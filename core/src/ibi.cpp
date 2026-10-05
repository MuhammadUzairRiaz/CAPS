#include "caps/ibi.hpp"

#include <cstdio>
#include <algorithm>
#include <cmath>
#include <set>
#include <stdexcept>

#include "caps/dynamics.hpp"

namespace caps {

namespace {
constexpr double kB = 0.0019872043;   // kcal/mol/K
constexpr double kPiI = 3.14159265358979323846;
}  // namespace

std::vector<double> nonbonded_gr(const System& top, const std::vector<std::vector<Vec3>>& frames, const std::vector<Cell>& cells, double dr,
                                 size_t nbin) {
  const size_t n = top.atoms.size();
  const auto nb = top.neighbours();
  const auto mol = top.molecules();
  // within three bonds of each bead (the excluded intramolecular pairs)
  std::vector<std::set<uint32_t>> near(n);
  for (uint32_t x = 0; x < n; ++x) {
    std::vector<uint32_t> front{x};
    for (int step = 0; step < 3; ++step) {
      std::vector<uint32_t> next;
      for (uint32_t u : front)
        for (uint32_t w : nb[u])
          if (w != x && near[x].insert(w).second) next.push_back(w);
      front = next;
    }
  }
  std::vector<double> h(nbin, 0.0);
  double pairs = 0, vol = 0;
  const double rmax = dr * double(nbin);
  for (size_t f = 0; f < frames.size(); ++f) {
    const Cell& c = cells[f];
    vol += c.volume();
    const auto& x = frames[f];
    for (size_t a = 0; a < n; ++a)
      for (size_t b = a + 1; b < n; ++b) {
        if (mol[a] == mol[b] && near[a].count(uint32_t(b))) continue;
        pairs += 1;
        const double d = norm(c.minimum_image(x[b] - x[a]));
        const size_t k = size_t(d / dr);   // d just under rmax can round to nbin
        if (d < rmax && k < nbin) h[k] += 1;
      }
  }
  std::vector<double> g(nbin, 0.0);
  if (frames.empty() || pairs <= 0) return g;
  const double V = vol / double(frames.size()), per_frame = pairs / double(frames.size());
  for (size_t k = 0; k < nbin; ++k) {
    const double r0 = double(k) * dr, r1 = r0 + dr, shell = 4.0 / 3.0 * kPiI * (r1 * r1 * r1 - r0 * r0 * r0);
    g[k] = h[k] / double(frames.size()) / (per_frame * shell / V);
  }
  return g;
}

IbiResult run_ibi(const System& beads, const ForceField& ff0, const std::vector<double>& r, const std::vector<double>& g, const IbiOptions& o) {
  if (!beads.cell.valid()) throw std::invalid_argument("IBI needs a periodic cell of beads");
  if (r.size() < 10 || r.size() != g.size()) throw std::invalid_argument("IBI needs a target g(r) on an even grid");
  const double drt = r[1] - r[0];
  const double kT = kB * o.temperature;
  double rc = o.cutoff > 0 ? std::min(o.cutoff, r.back()) : std::min(r.back(), 15.0);
  const double half = 0.5 * std::min({norm(beads.cell.a), norm(beads.cell.b), norm(beads.cell.c)});
  rc = std::min(rc, half);
  const size_t nt = std::min(r.size(), size_t((rc - r[0]) / drt) + 1);   // target points inside the cut-off
  IbiResult res;
  res.r.assign(r.begin(), r.begin() + long(nt));
  res.target.assign(g.begin(), g.begin() + long(nt));
  // the start: −k_B T ln g where resolved; first resolved point: g ≥ 0.02
  size_t k0 = 0;
  while (k0 < nt && res.target[k0] < 0.02) ++k0;
  if (k0 + 5 >= nt) throw std::invalid_argument("the target g(r) is not resolved inside the cut-off");
  std::vector<double> U(nt, 0.0);
  for (size_t k = k0; k < nt; ++k) U[k] = -kT * std::log(std::max(res.target[k], 1e-6));
  auto close = [&](std::vector<double>& u) {
    // below the first resolved point: a straight wall at least as steep as 5 k_B T/Å
    const double slope = std::max(5 * kT, (u[k0] - u[k0 + 1]) / drt);
    for (size_t k = 0; k < k0; ++k) u[k] = u[k0] + slope * double(k0 - k) * drt;
    const double end = u[nt - 1];
    for (auto& v : u) v -= end;   // zero at the cut-off
  };
  close(U);
  // the table from U: linear in r between the target's points, the force by central differences
  auto table_of = [&](const std::vector<double>& u) {
    TabulatedPair tb;
    tb.r0 = 0.05;   // from close to zero: LAMMPS stops on a pair below a table's start
    tb.dr = o.dr;
    for (double x = tb.r0; x <= res.r.back() + 1e-9; x += tb.dr) {
      double e;
      if (x <= res.r.front()) e = u.front() + (u[0] - u[1]) / drt * (res.r.front() - x);
      else {
        const double s = (x - res.r.front()) / drt;
        const size_t k = std::min(nt - 2, size_t(s));
        const double w = s - double(k);
        e = u[k] * (1 - w) + u[k + 1] * w;
      }
      tb.e.push_back(e);
    }
    tb.f.assign(tb.e.size(), 0.0);
    for (size_t k = 0; k < tb.e.size(); ++k) {
      const size_t a = k == 0 ? 0 : k - 1, b = std::min(tb.e.size() - 1, k + 1);
      tb.f[k] = -(tb.e[b] - tb.e[a]) / (double(b - a) * tb.dr);
    }
    return tb;
  };
  // the model: the bonded terms as given, every bead pair on the table
  auto ff = std::make_shared<ForceField>(ff0);
  ff->pair_func.clear();
  ff->tables = {table_of(U)};
  const int ntype = int(ff->type_names.size());
  for (int a = 0; a < ntype; ++a)
    for (int b = a; b < ntype; ++b) ff->pair_func[{a, b}] = {kPairTable, 0, 0, 0};
  std::fill(ff->charge.begin(), ff->charge.end(), 0.0);
  if (ntype > 1) res.notes.push_back(std::to_string(ntype) + " bead types share one table (the pooled g(r))");
  System s = beads;
  s.velocities.clear();
  DynamicsOptions d;
  d.energy.cutoff = rc;
  d.energy.coulomb = false;
  d.energy.tail = false;
  d.dt = o.dt;
  d.temperature = o.temperature;
  d.thermostat = Thermostat::Bussi;
  d.tau_t = 200;
  d.seed = o.seed;
  d.frame_every = std::max(1, o.frame_every);
  d.thermo_every = std::max(1, o.frame_every);
  const int64_t nrun = std::max<int64_t>(10, std::llround(o.run_ps * 1000 / o.dt));
  const int64_t neq = std::llround(o.first_equilibrate_ps * 1000 / o.dt);
  for (int it = 0; it <= o.iterations; ++it) {
    if (o.progress && !o.progress("iteration " + std::to_string(it) + " of " + std::to_string(o.iterations), double(it) / (o.iterations + 1)))
      throw std::runtime_error("IBI cancelled");
    ff->tables = {table_of(U)};
    d.field = ff;
    d.new_velocities = it == 0;
    if (it == 0 && neq > 0) {   // the first run settles the start
      d.steps = neq;
      d.frame = nullptr;
      d.progress = nullptr;
      run_dynamics(s, d);
      d.new_velocities = false;
    }
    std::vector<std::vector<Vec3>> frames;
    std::vector<Cell> cells;
    double psum = 0;
    long np = 0;
    d.steps = nrun;
    d.frame = [&](const std::vector<double>& x, const Cell& c, int64_t step) {
      if (step == 0) return;
      std::vector<Vec3> p(x.size() / 3);
      for (size_t i = 0; i < p.size(); ++i) p[i] = {x[3 * i], x[3 * i + 1], x[3 * i + 2]};
      frames.push_back(std::move(p));
      cells.push_back(c);
    };
    d.progress = [&](const ThermoRow& row) { psum += row.pressure, ++np; return true; };
    run_dynamics(s, d);
    const auto gi = nonbonded_gr(s, frames, cells, drt, nt);
    IbiIteration rec;
    rec.iteration = it;
    rec.gr = gi;
    rec.pressure = np > 0 ? psum / double(np) : 0;
    double num = 0, den = 0;
    for (size_t k = 0; k < nt; ++k) num += (gi[k] - res.target[k]) * (gi[k] - res.target[k]), den += res.target[k] * res.target[k];
    rec.residual = den > 0 ? num / den : 0;
    res.history.push_back(rec);
    if (it == o.iterations) break;
    // the update
    std::vector<double> dU(nt, 0.0);
    for (size_t k = k0; k < nt; ++k)
      if (gi[k] > 1e-3 && res.target[k] > 1e-3) dU[k] = o.alpha * kT * std::log(gi[k] / res.target[k]);
    std::vector<double> sm = dU;
    for (size_t k = k0 + 1; k + 1 < nt; ++k) sm[k] = (dU[k - 1] + dU[k] + dU[k + 1]) / 3;
    for (size_t k = k0; k < nt; ++k) U[k] += sm[k];
    if (o.pressure_correction) {   // ΔV = A (1 − r/r_c), A = −0.1 k_B T sign(ΔP) min(1, 0.0003 |ΔP| / bar)
      const double dP = (rec.pressure - o.pressure) * 1.01325;   // atm → bar
      const double A = -0.1 * kT * (dP > 0 ? 1 : -1) * std::min(1.0, 0.0003 * std::fabs(dP));
      for (size_t k = k0; k < nt; ++k) U[k] += A * (1 - res.r[k] / res.r.back());
    }
    close(U);
  }
  ff->native_timestep = o.dt;
  ff->cutoff = rc;   // the table's end is the model's own cut-off: every run evaluates it whole, as LAMMPS does
  ff->name += " · non-bonded by IBI (tabulated)";   // the time step the model was refined and checked at (its wall is steeper than the WCA's)
  res.ff = ff;
  res.last = s;
  const auto tb = ff->tables.front();
  for (size_t k = 0; k < tb.e.size(); ++k) res.table_r.push_back(tb.r0 + double(k) * tb.dr), res.potential.push_back(tb.e[k]);
  char b[200];
  std::snprintf(b, sizeof b, "IBI: %d iterations of %.0f ps at %.0f K · residual %.4f → %.4f · cut-off %.1f Å", o.iterations, o.run_ps, o.temperature,
                res.history.front().residual, res.history.back().residual, rc);
  res.notes.insert(res.notes.begin(), b);
  return res;
}

}  // namespace caps
