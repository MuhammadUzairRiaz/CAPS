// CAPS Dynamics: velocity Verlet (Swope et al., J. Chem. Phys. 76, 637 (1982)); Bussi–Donadio–Parrinello velocity
// rescaling (J. Chem. Phys. 126, 014101 (2007)); Langevin dynamics with the BAOAB splitting (Leimkuhler and Matthews,
// Appl. Math. Res. Express 2013, 34); stochastic cell rescaling (Bernetti and Bussi, J. Chem. Phys. 153, 114107
// (2020)) and Berendsen coupling (J. Chem. Phys. 81, 3684 (1984)) for isotropic pressure control.
#include "caps/uff.hpp"
#include "caps/dynamics.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <memory>
#include <cstdio>
#include <random>

namespace caps {

namespace {

constexpr double kB = 0.0019872041;       // kcal/(mol·K)
constexpr double kAcc = 4.184e-4;         // (kcal/mol/Å)/(g/mol) → Å/fs²
constexpr double kAtm = 68568.415;        // kcal/(mol·Å³) → atm
constexpr double kNA = 6.02214076e23;

}  // namespace

Thermostat thermostat_from_string(const std::string& s) {
  if (s == "none" || s == "nve") return Thermostat::None;
  if (s == "bussi" || s == "csvr" || s == "v-rescale") return Thermostat::Bussi;
  if (s == "langevin" || s == "baoab") return Thermostat::Langevin;
  throw std::invalid_argument("unknown thermostat '" + s + "' (none, bussi, langevin)");
}

Barostat barostat_from_string(const std::string& s) {
  if (s == "none") return Barostat::None;
  if (s == "crescale" || s == "c-rescale") return Barostat::CRescale;
  if (s == "berendsen") return Barostat::Berendsen;
  throw std::invalid_argument("unknown barostat '" + s + "' (none, crescale, berendsen)");
}

const char* to_string(Thermostat t) {
  switch (t) {
    case Thermostat::None: return "none (NVE)";
    case Thermostat::Bussi: return "Bussi velocity rescaling";
    case Thermostat::Langevin: return "Langevin (BAOAB)";
  }
  return "?";
}

const char* to_string(Barostat b) {
  switch (b) {
    case Barostat::None: return "none";
    case Barostat::CRescale: return "stochastic cell rescaling";
    case Barostat::Berendsen: return "Berendsen";
  }
  return "?";
}

double kinetic_energy(const std::vector<double>& v, const std::vector<double>& m) {
  double k = 0;
  for (size_t i = 0; i < m.size(); ++i) k += m[i] * (v[3 * i] * v[3 * i] + v[3 * i + 1] * v[3 * i + 1] + v[3 * i + 2] * v[3 * i + 2]);
  return 0.5 * k / kAcc;
}

void run_dynamics(System& s, const DynamicsOptions& o, DynamicsReport* rep_out) {
  DynamicsReport rep;
  const size_t n = s.atoms.size();
  if (n < 2) throw FieldError("dynamics needs at least two atoms");
  // coarse-grained beads (every site element 0) move slowly enough for the model's own steps (SDK 10 fs, MARTINI 20–40 fs)
  const bool beads = std::all_of(s.atoms.begin(), s.atoms.end(), [](const Atom& a) { return a.element == 0; });
  const double dt_max = beads ? 50.0 : 5.0;
  if (o.dt <= 0 || o.dt > dt_max)
    throw std::invalid_argument(beads ? "time step must be in (0, 50] fs for coarse-grained beads" : "time step must be in (0, 5] fs");
  if (o.barostat != Barostat::None && !s.cell.valid()) throw FieldError("pressure control needs a periodic cell");
  if (o.thermostat == Thermostat::None && o.barostat == Barostat::CRescale)
    throw std::invalid_argument("stochastic cell rescaling needs a thermostat");
  if (o.anisotropic && o.barostat == Barostat::CRescale)
    throw std::invalid_argument("per-axis pressure coupling is available with the Berendsen barostat");
  if (o.deform_axis >= 0 && (o.deform_axis > 2 || !s.cell.valid())) throw std::invalid_argument("deformation needs a periodic cell and an axis 0, 1 or 2");
  if (o.deform_axis >= 0 && o.barostat != Barostat::None && (!o.anisotropic || o.couple_axis[o.deform_axis]))
    throw std::invalid_argument("the deformed axis cannot also follow the barostat: use per-axis coupling without that axis");
  if (o.respa < 1 || o.respa > 16) throw std::invalid_argument("r-RESPA inner steps must be 1 … 16");
  if (o.respa > 1 && o.thermostat == Thermostat::Langevin) throw std::invalid_argument("r-RESPA runs with the Bussi thermostat or none");

  if (o.field && o.field->atom_type.size() != n)
    throw FieldError("the assigned force field is for " + std::to_string(o.field->atom_type.size()) + " atoms, the structure has " +
                     std::to_string(n));
  const ForceField ff = o.field ? *o.field : default_forcefield(s);
  for (const auto& note : ff.notes) rep.notes.push_back(note);
  // r-RESPA (Tuckerman, Berne & Martyna 1992): the non-bonded forces every step, the bonded ones respa times per step
  const int respa = std::max(1, o.respa);
  EnergyOptions eo = o.energy;
  if (respa > 1) eo.parts = 2;
  Evaluator ev(ff, eo);
  std::unique_ptr<Evaluator> ev_fast;
  if (respa > 1) {
    EnergyOptions ef = o.energy;
    ef.parts = 1;
    ev_fast = std::make_unique<Evaluator>(ff, ef);
  }
  // virtual sites (massless, placed from their atoms) are not integrated: held, with a nominal mass for the arithmetic
  std::vector<double> m = ff.mass;
  std::vector<char> held(n, 0);
  for (const auto& vs : ff.vsites) held[vs.site] = 2, m[vs.site] = 1.0;
  double mtot = 0;
  for (size_t i = 0; i < n; ++i) mtot += held[i] == 2 ? 0.0 : m[i];
  size_t nheld = 0, nvs = ff.vsites.size();
  for (size_t i = 0; i < n && i < o.fixed.size(); ++i)
    if (o.fixed[i] && !held[i]) held[i] = 1, ++nheld;
  nheld += nvs;
  if (nheld + 1 >= n) throw std::invalid_argument("dynamics needs at least two atoms that move");
  // with held atoms momentum is not conserved: every free coordinate counts
  const double ndof = nheld ? 3.0 * double(n - nheld) : 3.0 * n - 3.0;
  if (nheld > nvs) rep.notes.push_back(std::to_string(nheld - nvs) + " atoms held in place");
  if (nvs) rep.notes.push_back(std::to_string(nvs) + " virtual sites placed from their atoms each step");

  Cell cell = s.cell;
  std::vector<double> x(3 * n), v(3 * n, 0.0), f, f_slow, f_fast;   // f_slow, f_fast: r-RESPA's two parts (f is their sum)
  for (size_t i = 0; i < n; ++i)
    for (int k = 0; k < 3; ++k) x[3 * i + k] = s.atoms[i].pos[k];

  std::mt19937_64 rng(o.seed);
  std::normal_distribution<double> gauss(0.0, 1.0);

  // Velocities: from the system, or Maxwell–Boltzmann at the target temperature without centre-of-mass motion.
  if (s.velocities.size() == n && !o.new_velocities) {
    for (size_t i = 0; i < n; ++i)
      for (int k = 0; k < 3; ++k) v[3 * i + k] = s.velocities[i][k];
    rep.notes.push_back("velocities taken from the structure");
  } else {
    for (size_t i = 0; i < n; ++i) {
      const double sd = std::sqrt(kB * o.temperature * kAcc / m[i]);
      for (int k = 0; k < 3; ++k) v[3 * i + k] = sd * gauss(rng);
    }
    double p[3] = {0, 0, 0}, mfree = 0;
    for (size_t i = 0; i < n; ++i)
      if (!held[i]) {
        mfree += m[i];
        for (int k = 0; k < 3; ++k) p[k] += m[i] * v[3 * i + k];
      }
    for (size_t i = 0; i < n; ++i)
      for (int k = 0; k < 3; ++k) v[3 * i + k] = held[i] ? 0.0 : v[3 * i + k] - p[k] / mfree;
    const double t0 = 2 * kinetic_energy(v, m) / (ndof * kB);
    if (t0 > 0)
      for (auto& q : v) q *= std::sqrt(o.temperature / t0);
    char b[96];
    std::snprintf(b, sizeof b, "velocities drawn at %.1f K (seed %llu)", o.temperature, static_cast<unsigned long long>(o.seed));
    rep.notes.push_back(b);
  }

  if (nheld)
    for (size_t i = 0; i < n; ++i)
      if (held[i]) v[3 * i] = v[3 * i + 1] = v[3 * i + 2] = 0;
  // steered pulling: the group, its mass, and where its centre starts along the pull direction
  std::vector<uint32_t> pulled;
  double mpull = 0, com0 = 0, pull_f = 0, pull_x = 0;
  int64_t cur_step = 0;
  const double plen = norm(o.pull_dir);
  const Vec3 pdir = plen > 0 ? o.pull_dir * (1 / plen) : Vec3{1, 0, 0};
  auto pull_com = [&] {
    double c = 0;
    for (uint32_t i : pulled) c += m[i] * (x[3 * i] * pdir[0] + x[3 * i + 1] * pdir[1] + x[3 * i + 2] * pdir[2]);
    return c / mpull;
  };
  for (size_t i = 0; i < n && i < o.pull_group.size(); ++i)
    if (o.pull_group[i] && !held[i]) pulled.push_back(uint32_t(i)), mpull += m[i];
  if (!pulled.empty()) com0 = pull_com();
  // forces, with none on held atoms and the pulling spring on the group
  auto hold = [&](std::vector<double>& F) {
    if (nheld)
      for (size_t i = 0; i < n; ++i)
        if (held[i]) F[3 * i] = F[3 * i + 1] = F[3 * i + 2] = 0;
  };
  EnergyTerms et_fast;
  auto compute_fast = [&] {   // r-RESPA inner step: bonded forces only
    et_fast = ev_fast->compute(x, cell, f_fast);
    hold(f_fast);
  };
  auto compute = [&] {
    EnergyTerms t = ev.compute(x, cell, respa > 1 ? f_slow : f);
    auto& F = respa > 1 ? f_slow : f;
    hold(F);
    if (!pulled.empty()) {
      pull_x = pull_com() - com0;
      const double anchor = o.pull_rate * double(cur_step) * o.dt * 1e-3;
      pull_f = o.pull_k * (anchor - pull_x);
      for (uint32_t i : pulled)
        for (int k = 0; k < 3; ++k) F[3 * i + k] += pdir[k] * pull_f * m[i] / mpull;
    }
    if (respa > 1) {   // both parts, for the energies, the virial and anyone reading f
      compute_fast();
      t.bond = et_fast.bond, t.angle = et_fast.angle, t.dihedral = et_fast.dihedral, t.improper = et_fast.improper;
      t.virial += et_fast.virial;
      for (int c = 0; c < 6; ++c) t.w[c] += et_fast.w[c];
      f.resize(f_slow.size());
      for (size_t k = 0; k < f.size(); ++k) f[k] = f_slow[k] + f_fast[k];
    }
    return t;
  };
  EnergyTerms et = compute();
  double bath = 0;   // energy taken out of the system by the thermostat and barostat, kcal/mol

  auto density = [&] { return cell.valid() ? mtot / kNA / (cell.volume() * 1e-24) : 0.0; };
  const bool ramp = o.temperature_end >= 0 && o.steps > 0;
  auto target_t = [&](int64_t step) {
    return ramp ? o.temperature + (o.temperature_end - o.temperature) * double(step) / double(o.steps) : o.temperature;
  };
  auto row = [&](int64_t step) {
    ThermoRow r;
    r.step = step + o.step_offset;
    r.time_ps = (step + o.step_offset) * o.dt * 1e-3;
    r.target_temperature = target_t(step);
    r.kinetic = kinetic_energy(v, m);
    r.temperature = 2 * r.kinetic / (ndof * kB);
    r.potential = et.total();
    r.total = r.potential + r.kinetic;
    r.volume = cell.valid() ? cell.volume() : 0;
    r.density = density();
    r.pressure = cell.valid() ? (2 * r.kinetic + et.virial) / (3 * r.volume) * kAtm : 0;
    if (cell.valid()) {
      double kt[6] = {0, 0, 0, 0, 0, 0};
      static const int ia[6] = {0, 1, 2, 0, 0, 1}, ib[6] = {0, 1, 2, 1, 2, 2};
      for (size_t i = 0; i < n; ++i)
        for (int c = 0; c < 6; ++c) kt[c] += m[i] * v[3 * i + ia[c]] * v[3 * i + ib[c]];
      for (int c = 0; c < 6; ++c) r.p[c] = (kt[c] / kAcc + et.w[c]) / r.volume * kAtm;
      r.lx = norm(cell.a);
      r.ly = norm(cell.b);
      r.lz = norm(cell.c);
    }
    r.conserved = r.total + bath + (o.barostat != Barostat::None ? o.pressure * r.volume / kAtm : 0.0);
    r.pull_force = pull_f;
    r.pull_disp = pull_x;
    return r;
  };
  auto emit = [&](int64_t step) {
    ThermoRow r = row(step);
    rep.thermo.push_back(r);
    if (o.progress && !o.progress(r)) throw DynamicsCancelled();
  };

  const double dt = o.dt;
  const double c1 = std::exp(-dt / std::max(o.tau_t, 1e-9)), c2 = std::sqrt(1 - c1 * c1);   // Langevin O step
  std::vector<double> sdv(n);
  for (size_t i = 0; i < n; ++i) sdv[i] = std::sqrt(kB * o.temperature * kAcc / m[i]);
  double kt_target = 0.5 * ndof * kB * o.temperature;
  double t_now = o.temperature;

  auto kick_with = [&](const std::vector<double>& F, double h) {
    for (size_t i = 0; i < n; ++i) {
      const double a = h * kAcc / m[i];
      v[3 * i] += a * F[3 * i];
      v[3 * i + 1] += a * F[3 * i + 1];
      v[3 * i + 2] += a * F[3 * i + 2];
    }
  };
  auto kick = [&](double h) { kick_with(f, h); };
  auto drift = [&](double h) {
    for (size_t k = 0; k < x.size(); ++k) x[k] += h * v[k];
  };

  // Bussi et al. 2007, appendix: the new kinetic energy is drawn from its canonical distribution in one step.
  auto bussi = [&](double h) {
    const double k = kinetic_energy(v, m);
    if (k <= 0) return;
    const double c = std::exp(-h / o.tau_t);
    const double r1 = gauss(rng);
    std::gamma_distribution<double> chi((ndof - 1) / 2, 2.0);   // Σ_{i≥2} R_i², chi-squared with ndof − 1
    const double sum = ndof > 1 ? chi(rng) : 0.0;
    const double knew = k + (1 - c) * (kt_target * (r1 * r1 + sum) / ndof - k) + 2 * r1 * std::sqrt(c * (1 - c) * kt_target * k / ndof);
    const double alpha = std::sqrt(std::max(knew, 0.0) / k);
    for (auto& q : v) q *= alpha;
    bath -= knew - k;
  };

  auto scale_cell = [&](double mu) {
    for (size_t i = 0; i < n; ++i)
      for (int k = 0; k < 3; ++k) x[3 * i + k] = cell.origin[k] + mu * (x[3 * i + k] - cell.origin[k]);
    cell.a = cell.a * mu;
    cell.b = cell.b * mu;
    cell.c = cell.c * mu;
  };

  // diagonal scaling of Cartesian axis k (atoms about the cell origin, and every cell vector's k component)
  auto scale_axis = [&](int k, double mu) {
    for (size_t i = 0; i < n; ++i) x[3 * i + k] = cell.origin[k] + mu * (x[3 * i + k] - cell.origin[k]);
    cell.a[k] *= mu;
    cell.b[k] *= mu;
    cell.c[k] *= mu;
  };
  auto deform = [&](int64_t step) {
    if (o.deform_axis < 0) return;
    const double t1 = step * dt * 1e-3, t0 = (step - 1) * dt * 1e-3;
    scale_axis(o.deform_axis, (1 + o.deform_rate * t1) / (1 + o.deform_rate * t0));
  };

  const auto t_start = std::chrono::steady_clock::now();
  emit(0);
  if (o.frame && o.frame_every > 0) o.frame(x, cell, o.step_offset);

  for (int64_t step = 1; step <= o.steps; ++step) {
    cur_step = step;
    if (ramp) {
      t_now = target_t(step);
      kt_target = 0.5 * ndof * kB * t_now;
      for (size_t i = 0; i < n; ++i) sdv[i] = std::sqrt(kB * t_now * kAcc / m[i]);
    }
    if (o.thermostat == Thermostat::Langevin) {
      kick(0.5 * dt);
      drift(0.5 * dt);
      const double k0 = kinetic_energy(v, m);
      for (size_t i = 0; i < n; ++i)
        for (int k = 0; k < 3; ++k) v[3 * i + k] = held[i] ? 0.0 : c1 * v[3 * i + k] + c2 * sdv[i] * gauss(rng);
      bath -= kinetic_energy(v, m) - k0;
      drift(0.5 * dt);
      deform(step);
      et = compute();
      kick(0.5 * dt);
    } else if (respa > 1) {
      // outer half kick with the non-bonded forces, `respa` velocity-Verlet steps with the bonded ones, outer half kick
      const double h = dt / respa;
      kick_with(f_slow, 0.5 * dt);
      for (int j = 0; j < respa; ++j) {
        kick_with(f_fast, 0.5 * h);
        drift(h);
        if (j == respa - 1) deform(step);
        if (j + 1 < respa) compute_fast();
        else et = compute();   // the last inner step: both parts at the new positions
        kick_with(f_fast, 0.5 * h);
      }
      kick_with(f_slow, 0.5 * dt);
      if (o.thermostat == Thermostat::Bussi) bussi(dt);
    } else {
      kick(0.5 * dt);
      drift(dt);
      deform(step);
      et = compute();
      kick(0.5 * dt);
      if (o.thermostat == Thermostat::Bussi) bussi(dt);
    }

    if (o.barostat != Barostat::None && o.anisotropic && step % std::max(1, o.barostat_every) == 0) {
      // Berendsen per axis: μ_k = exp(−β hp (P0 − P_kk) / (3τ)), the linear-strain form of the isotropic update
      const double hp = dt * std::max(1, o.barostat_every);
      const ThermoRow r = row(step);
      for (int k = 0; k < 3; ++k) {
        if (!o.couple_axis[k]) continue;
        const double deps = std::clamp(-o.compressibility / o.tau_p * (o.pressure - r.p[k]) * hp / 3, -0.0033, 0.0033);
        scale_axis(k, std::exp(deps));
      }
      et = compute();
    } else if (o.barostat != Barostat::None && step % std::max(1, o.barostat_every) == 0) {
      const double hp = dt * std::max(1, o.barostat_every);
      const double k = kinetic_energy(v, m);
      const double vol = cell.volume();
      const double p = (2 * k + et.virial) / (3 * vol) * kAtm;
      double deps = -o.compressibility / o.tau_p * (o.pressure - p) * hp;
      if (o.barostat == Barostat::CRescale)
        deps += std::sqrt(2 * kB * t_now / vol * kAtm * o.compressibility * hp / o.tau_p) * gauss(rng);
      deps = std::clamp(deps, -0.01, 0.01);   // at most 1% in volume per update
      const double mu = std::exp(deps / 3);
      scale_cell(mu);
      if (o.barostat == Barostat::CRescale) {
        for (auto& q : v) q /= mu;
        bath -= kinetic_energy(v, m) - k;
      }
      et = compute();
    }

    if (o.each_step) o.each_step(et, x, cell, step + o.step_offset);
    if (step % std::max(1, o.thermo_every) == 0 || step == o.steps) emit(step);
    if (o.frame && o.frame_every > 0 && (step % o.frame_every == 0 || step == o.steps)) o.frame(x, cell, step + o.step_offset);
    for (size_t k = 0; k < x.size(); k += 3 * 64)
      if (!std::isfinite(x[k])) throw FieldError("the simulation became unstable at step " + std::to_string(step) + "; lower the time step or relax first");
  }

  rep.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t_start).count();
  rep.steps = o.steps;
  rep.ns_per_day = rep.seconds > 0 ? o.steps * dt * 1e-6 / rep.seconds * 86400 : 0;
  rep.list_builds = ev.list_builds();
  char b[256];
  std::snprintf(b, sizeof b, "%lld steps of %.2f fs (%.1f ps) in %.1f s · %.2f ns/day · %d neighbour-list builds · %d threads",
                static_cast<long long>(o.steps), dt, o.steps * dt * 1e-3, rep.seconds, rep.ns_per_day, rep.list_builds, ev.threads());
  rep.notes.insert(rep.notes.begin(), b);
  std::snprintf(b, sizeof b, "thermostat %s (τ %.0f fs) · barostat %s", to_string(o.thermostat), o.tau_t, to_string(o.barostat));
  rep.notes.insert(rep.notes.begin() + 1, b);
  if (respa > 1) {
    std::snprintf(b, sizeof b, "r-RESPA: non-bonded forces every %.2f fs, bonded forces every %.3f fs (%d inner steps)", dt, dt / respa, respa);
    rep.notes.insert(rep.notes.begin() + 2, b);
  }
  if (o.barostat != Barostat::None) {
    std::snprintf(b, sizeof b, "barostat target %.1f atm, τ %.0f fs, compressibility %.2e atm⁻¹, every %d steps%s; LJ tail correction %s",
                  o.pressure, o.tau_p, o.compressibility, o.barostat_every, o.anisotropic ? ", each axis on its own" : "", o.energy.tail ? "on" : "off");
    rep.notes.push_back(b);
  }

  place_virtual_sites(ff, x, cell);
  for (size_t i = 0; i < n; ++i)
    for (int k = 0; k < 3; ++k) s.atoms[i].pos[k] = x[3 * i + k];
  s.velocities.assign(n, Vec3{0, 0, 0});
  for (size_t i = 0; i < n; ++i) s.velocities[i] = {v[3 * i], v[3 * i + 1], v[3 * i + 2]};
  for (size_t i = 0; i < n; ++i) s.atoms[i].charge = ff.charge[i];
  s.has_charges = true;
  s.cell = cell;
  s.unwrapped = true;
  if (rep_out) *rep_out = std::move(rep);
}

}  // namespace caps
