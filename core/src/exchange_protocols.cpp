// CAPS Exchange, protocols (see exchange_protocols.hpp).
#include "caps/exchange_protocols.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace caps {

namespace {

constexpr double kAtmMPa = 0.101325;

void check_axis(int axis) {
  if (axis < 0 || axis > 2) throw std::invalid_argument("axis must be 0, 1 or 2 (x, y, z)");
}

Vec3& cell_vector(Cell& c, int axis) { return axis == 0 ? c.a : axis == 1 ? c.b : c.c; }

// Whether a bond crosses the plane at fraction `plane` of the cell along the axis (by the minimum image).
bool crosses(const System& s, const Bond& b, int axis, double plane) {
  const Vec3 pi = s.atoms[b.i].pos, d = s.cell.minimum_image(s.atoms[b.j].pos - pi);
  double fi = s.cell.to_fractional(pi)[axis];
  const double fj0 = s.cell.to_fractional(pi + d)[axis] - fi;
  fi -= std::floor(fi);
  const double lo = std::min(fi, fi + fj0), hi = std::max(fi, fi + fj0);
  return std::floor(hi - plane) - std::floor(lo - plane) > 0;
}

// Keeps the atoms marked, renumbering bonds; per-atom extras that no longer fit go.
void keep_atoms(System& s, const std::vector<char>& keep) {
  std::vector<int64_t> to(s.atoms.size(), -1);
  std::vector<Atom> atoms;
  std::vector<Vec3> vel;
  for (size_t i = 0; i < s.atoms.size(); ++i)
    if (keep[i]) {
      to[i] = int64_t(atoms.size());
      atoms.push_back(s.atoms[i]);
      if (s.velocities.size() == s.atoms.size()) vel.push_back(s.velocities[i]);
    }
  auto renumber = [&](std::vector<Bond>& bonds) {
    std::vector<Bond> out;
    for (const auto& b : bonds)
      if (to[b.i] >= 0 && to[b.j] >= 0) out.push_back({uint32_t(to[b.i]), uint32_t(to[b.j]), b.order});
    bonds = std::move(out);
  };
  renumber(s.bonds);
  renumber(s.typing_bonds);
  s.atoms = std::move(atoms);
  s.velocities = std::move(vel);
  s.topology = nullptr, s.forcefield = nullptr;
}

// One run with the protocol's settings and exchange moves; the force field the moves left becomes the protocol's.
DynamicsReport run(System& s, ProtocolOptions& o, DynamicsOptions d, ExchangeStats& total) {
  d.field = o.field;
  d.energy = o.energy;
  d.dt = o.dt;
  d.tau_t = o.tau_t;
  d.tau_p = o.tau_p;
  d.seed = o.seed++;
  auto st = std::make_shared<ExchangeStats>();
  if (o.exchange) {
    ExchangeOptions m = o.moves;
    m.seed = o.moves.seed++;
    d.topology_every = o.exchange_every;
    d.topology = exchange_callback(m, o.energy, st);
  }
  DynamicsReport r;
  run_dynamics(s, d, &r);
  if (r.field) o.field = r.field;
  total.attempts += st->attempts, total.accepted += st->accepted, total.no_candidate += st->no_candidate;
  total.unparameterised += st->unparameterised, total.energy += st->energy;
  total.local = st->local;
  if (!st->note.empty()) total.note = st->note;
  return r;
}

int64_t steps_for(double ps, double dt) { return std::max<int64_t>(1, int64_t(std::llround(ps * 1000.0 / dt))); }

}  // namespace

size_t cut_bonds_across(System& s, int axis, double plane) {
  check_axis(axis);
  if (!s.cell.valid()) throw std::invalid_argument("a cut needs a periodic cell");
  const size_t before = s.bonds.size();
  auto across = [&](const Bond& b) { return crosses(s, b, axis, plane); };
  s.bonds.erase(std::remove_if(s.bonds.begin(), s.bonds.end(), across), s.bonds.end());
  s.typing_bonds.erase(std::remove_if(s.typing_bonds.begin(), s.typing_bonds.end(), across), s.typing_bonds.end());
  const size_t cut = before - s.bonds.size();
  if (cut) s.topology = nullptr, s.forcefield = nullptr;   // the terms given with the file no longer match the bonds
  return cut;
}

void open_gap(System& s, int axis, double plane, double gap) {
  check_axis(axis);
  if (!s.cell.valid()) throw std::invalid_argument("a gap needs a periodic cell");
  Vec3& v = cell_vector(s.cell, axis);
  const double len = norm(v);
  const Vec3 dir = v * (1.0 / len);
  for (auto& a : s.atoms) {
    double f = s.cell.to_fractional(a.pos)[axis];
    f -= std::floor(f);
    if (f >= plane) a.pos = a.pos + dir * gap;
  }
  v = v * ((len + gap) / len);
}

size_t remove_sphere(System& s, const Vec3& centre, double radius) {
  std::vector<char> keep(s.atoms.size(), 1);
  size_t removed = 0;
  for (size_t i = 0; i < s.atoms.size(); ++i) {
    Vec3 d = s.atoms[i].pos - centre;
    if (s.cell.valid()) d = s.cell.minimum_image(d);
    if (norm(d) < radius) keep[i] = 0, ++removed;
  }
  if (removed) keep_atoms(s, keep);
  return removed;
}

size_t fragment(System& s, const std::vector<int>& axes) {
  size_t n = 0;
  for (int a : axes) n += cut_bonds_across(s, a, 0.5);
  return n;
}

StepStrainResult step_strain(System& s, ProtocolOptions& o, int axis, double strain, double hold_ps, int sample_every) {
  check_axis(axis);
  if (!s.cell.valid()) throw std::invalid_argument("a step strain needs a periodic cell");
  if (!(strain > -1)) throw std::invalid_argument("the strain must be above −1");
  // affine, keeping the volume: the axis stretched by 1 + ε, the other two by 1/√(1 + ε)
  std::vector<Vec3> frac;
  for (const auto& a : s.atoms) frac.push_back(s.cell.to_fractional(a.pos));
  const double lat = 1.0 / std::sqrt(1.0 + strain);
  for (int k = 0; k < 3; ++k) cell_vector(s.cell, k) = cell_vector(s.cell, k) * (k == axis ? 1.0 + strain : lat);
  for (size_t i = 0; i < s.atoms.size(); ++i) s.atoms[i].pos = s.cell.to_cartesian(frac[i]);
  StepStrainResult r;
  DynamicsOptions d;
  d.temperature = o.temperature;
  d.steps = steps_for(hold_ps, o.dt);
  d.thermo_every = std::max(1, sample_every);
  const int b1 = (axis + 1) % 3, b2 = (axis + 2) % 3;
  d.progress = [&](const ThermoRow& row) {
    r.t_ps.push_back(row.time_ps);
    r.stress.push_back((-row.p[axis] + 0.5 * (row.p[b1] + row.p[b2])) * kAtmMPa);
    return true;
  };
  run(s, o, d, r.stats);
  const auto kw = fit_kww(r.t_ps, r.stress, 0);
  r.kww_tau = kw[0], r.kww_beta = kw[1], r.kww_mean = kw[2];
  return r;
}

TvScanResult tv_scan(System& s, ProtocolOptions& o, double t_from, double t_to, double t_step, double ps_per_step) {
  if (!(t_step > 0) || !(ps_per_step > 0)) throw std::invalid_argument("the Tv scan needs a positive step and time per step");
  TvScanResult r;
  ExchangeStats total;
  const double dir = t_to >= t_from ? 1.0 : -1.0;
  for (double T = t_from; dir * (t_to - T) >= -1e-9; T += dir * t_step) {
    DynamicsOptions d;
    d.temperature = T;
    d.barostat = Barostat::CRescale;
    d.pressure = o.pressure;
    d.steps = steps_for(ps_per_step, o.dt);
    d.thermo_every = int(std::max<int64_t>(1, d.steps / 50));
    std::vector<double> dens;
    d.progress = [&](const ThermoRow& row) {
      if (row.step > d.steps / 2 && row.density > 0) dens.push_back(row.density);
      return true;
    };
    const int64_t before = total.accepted;
    run(s, o, d, total);
    double rho = 0;
    for (double x : dens) rho += x / double(dens.size());
    r.T.push_back(T);
    r.specific_volume.push_back(rho > 0 ? 1.0 / rho : 0.0);
    r.accepted_per_ps.push_back(double(total.accepted - before) / ps_per_step);
  }
  r.fit = fit_bilinear(r.T, r.specific_volume);
  std::vector<double> tau;
  for (double k : r.accepted_per_ps) tau.push_back(k > 0 ? 1.0 / k : 0.0);
  r.rate = fit_arrhenius(r.T, tau);
  return r;
}

HealResult heal_run(System& s, ProtocolOptions& o, const std::vector<int8_t>& side, int axis, double plane, double ps, int chunks) {
  if (chunks < 1) throw std::invalid_argument("a heal run needs one chunk or more");
  if (side.size() != s.atoms.size()) throw std::invalid_argument("one side per atom");
  HealResult r;
  double t = 0;
  for (int c = 0; c < chunks; ++c) {
    DynamicsOptions d;
    d.temperature = o.temperature;
    d.steps = steps_for(ps / chunks, o.dt);
    d.thermo_every = int(std::max<int64_t>(1, d.steps));
    if (o.npt) d.barostat = Barostat::CRescale, d.pressure = o.pressure;
    run(s, o, d, r.stats);
    t += ps / chunks;
    const auto ib = interface_bonds(s, side, axis, plane, 1.0);
    r.t_ps.push_back(t), r.bridges.push_back(double(ib.bridges)), r.crossing.push_back(double(ib.crossing_plane));
  }
  return r;
}

}  // namespace caps
