// CAPS torsion scan (see caps/torsion.hpp).
#include "caps/torsion.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "caps/analysis.hpp"
#include "caps/relax.hpp"

namespace caps {

namespace {
constexpr double kDeg = 3.14159265358979323846 / 180.0;
double wrap180(double a) { a = std::fmod(a + 180.0, 360.0); if (a < 0) a += 360.0; return a - 180.0; }
}  // namespace

std::vector<uint32_t> moving_side(const System& s, int b, int c) {
  const size_t n = s.atoms.size();
  if (b < 0 || c < 0 || size_t(b) >= n || size_t(c) >= n) throw std::invalid_argument("atom out of range");
  std::vector<std::vector<uint32_t>> nb(n);
  bool bonded = false;
  for (const auto& e : s.bonds) {
    nb[e.i].push_back(e.j), nb[e.j].push_back(e.i);
    bonded |= (int(e.i) == b && int(e.j) == c) || (int(e.i) == c && int(e.j) == b);
  }
  if (!bonded) throw std::invalid_argument("the middle two atoms of a torsion scan must be bonded");
  std::vector<char> seen(n, 0);
  std::vector<uint32_t> out{uint32_t(c)}, stack{uint32_t(c)};
  seen[size_t(c)] = 1;
  seen[size_t(b)] = 1;
  while (!stack.empty()) {
    const uint32_t a = stack.back();
    stack.pop_back();
    for (uint32_t q : nb[a]) {
      if (a == uint32_t(c) && q == uint32_t(b)) continue;
      if (q == uint32_t(b)) throw std::invalid_argument("the b–c bond is in a ring: its torsion cannot be scanned by rotation");
      if (seen[q]) continue;
      seen[q] = 1;
      out.push_back(q);
      stack.push_back(q);
    }
  }
  return out;
}

double dihedral_angle(const System& s, const std::array<int, 4>& at) {
  auto P = [&](int k) { return s.atoms[size_t(at[size_t(k)])].pos; };
  Vec3 b0 = P(0) - P(1), b1 = P(2) - P(1), b2 = P(3) - P(2);
  if (s.cell.valid()) b0 = s.cell.minimum_image(b0), b1 = s.cell.minimum_image(b1), b2 = s.cell.minimum_image(b2);
  const Vec3 u = b1 * (1.0 / norm(b1));
  const Vec3 v = b0 - u * dot(b0, u), w = b2 - u * dot(b2, u);
  return std::atan2(dot(cross(u, v), w), dot(v, w)) / kDeg;
}

void set_dihedral(System& s, const std::array<int, 4>& at, double phi, const std::vector<uint32_t>& moving) {
  const double delta = wrap180(phi - dihedral_angle(s, at)) * kDeg;
  const Vec3 B = s.atoms[size_t(at[1])].pos;
  Vec3 axis = s.atoms[size_t(at[2])].pos - B;
  if (s.cell.valid()) axis = s.cell.minimum_image(axis);
  axis = axis * (1.0 / norm(axis));
  const double c = std::cos(delta), sn = std::sin(delta);
  for (uint32_t i : moving) {
    Vec3 p = s.atoms[i].pos - B;
    if (s.cell.valid()) p = s.cell.minimum_image(p);
    // Rodrigues: rotation by delta about the axis through B
    const Vec3 r = p * c + cross(axis, p) * sn + axis * (dot(axis, p) * (1 - c));
    s.atoms[i].pos = B + r;
  }
}

std::string torsion_state(double phi) {
  const double a = std::fabs(wrap180(phi));
  if (a >= 150) return "trans";
  if (a >= 90) return phi > 0 ? "anticlinal+" : "anticlinal−";
  if (a >= 30) return phi > 0 ? "gauche+" : "gauche−";
  return "cis";
}

TorsionScanResult torsion_scan(const System& s0, const ForceField& ff, const TorsionScanOptions& o) {
  TorsionScanResult R;
  const auto& at = o.atoms;
  const size_t n = s0.atoms.size();
  for (int k : at) if (k < 0 || size_t(k) >= n) throw std::invalid_argument("choose four atoms for the torsion");
  if (at[0] == at[1] || at[1] == at[2] || at[2] == at[3] || at[0] == at[3]) throw std::invalid_argument("the four atoms must differ");
  if (!(o.step > 0)) throw std::invalid_argument("the step must be positive");
  if (ff.atom_type.size() != n && ff.mass.size() != n) throw std::invalid_argument("the force field does not match the structure");
  R.moving = moving_side(s0, at[1], at[2]);
  R.phi_start = dihedral_angle(s0, at);
  System s = s0;
  EnergyOptions eo = o.energy;
  if (!s.cell.valid()) eo.tail = false;
  Evaluator ev(ff, eo);
  std::vector<double> x(3 * n), f(3 * n);
  const int steps = std::max(1, int(std::floor((o.to - o.from) / o.step + 1e-9)));
  const bool full_turn = o.to - o.from >= 360 - 1e-9;
  const int count = full_turn ? steps : steps + 1;   // −180 and 180 are one geometry
  RelaxOptions ro;
  ro.field = std::shared_ptr<const ForceField>(&ff, [](const ForceField*) {});
  ro.pushoff = false;
  ro.ftol = o.ftol;
  ro.max_iterations = 2000;
  ro.fixed.assign(n, 0);
  for (int k : at) ro.fixed[size_t(k)] = 1;
  for (int k = 0; k < count; ++k) {
    const double phi = o.from + k * o.step;
    set_dihedral(s, at, phi, R.moving);
    for (size_t i = 0; i < n; ++i) for (int d = 0; d < 3; ++d) x[3 * i + d] = s.atoms[i].pos[size_t(d)];
    TorsionPoint p;
    if (o.relax) {
      const RelaxStage st = minimise(ev, x, s.cell, ro, "scan");
      (void)st;
      for (size_t i = 0; i < n; ++i) s.atoms[i].pos = {x[3 * i], x[3 * i + 1], x[3 * i + 2]};
    }
    p.terms = ev.compute(x, s.cell, f);
    p.energy = p.terms.total();
    p.phi = dihedral_angle(s, at);
    for (const auto& a : s.atoms) p.positions.push_back(a.pos);
    R.points.push_back(std::move(p));
    if (o.progress && !o.progress(k + 1, count)) { R.notes.push_back("stopped after " + std::to_string(k + 1) + " points"); break; }
  }
  if (R.points.empty()) return R;
  double emin = 1e300, emax = -1e300;
  for (const auto& p : R.points) emin = std::min(emin, p.energy), emax = std::max(emax, p.energy);
  R.minimum = emin;
  R.barrier = emax - emin;
  // minima of the sampled curve, periodic over a full turn; a parabola through each minimum and its neighbours
  const size_t m = R.points.size();
  for (size_t k = 0; k < m; ++k) {
    const bool has_prev = full_turn || k > 0, has_next = full_turn || k + 1 < m;
    const auto& c = R.points[k];
    const auto& prev = R.points[(k + m - 1) % m];
    const auto& next = R.points[(k + 1) % m];
    if (has_prev && prev.energy <= c.energy) continue;
    if (has_next && next.energy < c.energy) continue;
    TorsionConformer cf;
    cf.phi = c.phi, cf.energy = c.energy;
    if (has_prev && has_next) {
      const double e0 = prev.energy, e1 = c.energy, e2 = next.energy;
      const double den = e0 - 2 * e1 + e2;
      if (den > 1e-12) {
        const double t = 0.5 * (e0 - e2) / den;   // offset in steps
        cf.phi = wrap180(c.phi + t * o.step);
        cf.energy = e1 - 0.25 * (e0 - e2) * t;
      }
    }
    cf.energy -= emin;
    cf.state = torsion_state(cf.phi);
    R.conformers.push_back(cf);
  }
  std::sort(R.conformers.begin(), R.conformers.end(), [](const TorsionConformer& a, const TorsionConformer& b) { return a.energy < b.energy; });
  if (o.relax) R.notes.push_back("relaxed scan: the four torsion atoms held, the rest minimised at each step");
  return R;
}

std::array<int, 4> default_torsion(const System& s) {
  auto chains = backbones(s, 4);
  std::sort(chains.begin(), chains.end(), [](const auto& a, const auto& b) { return a.size() > b.size(); });
  for (const auto& ch : chains) {
    const size_t mid = ch.size() / 2;
    for (size_t off = 0; off < ch.size(); ++off) {
      for (int sign : {1, -1}) {
        const long k = long(mid) + sign * long(off) - 1;
        if (k < 0 || size_t(k) + 3 >= ch.size() + 0 || size_t(k + 3) >= ch.size()) continue;
        try {
          moving_side(s, int(ch[size_t(k) + 1]), int(ch[size_t(k) + 2]));
          return {int(ch[size_t(k)]), int(ch[size_t(k) + 1]), int(ch[size_t(k) + 2]), int(ch[size_t(k) + 3])};
        } catch (const std::exception&) {}
      }
    }
  }
  return {-1, -1, -1, -1};
}

}  // namespace caps
