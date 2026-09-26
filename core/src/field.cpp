// GAFF typing and parameters for hydrocarbons, and the energy / force evaluator used by Relax.
//
// Parameters: GAFF 1.81 (Wang, Wolf, Caldwell, Kollman and Case, J. Comput. Chem. 25, 1157 (2004)), gaff.dat values
// for the c3 / ca / hc / ha subset. AMBER forms and 1-4 scaling (LJ 1/2, Coulomb 1/1.2). Electrostatics use the damped
// shifted force sum (Fennell and Gezelter, J. Chem. Phys. 124, 234104 (2006)), a real-space method with no Ewald part.
#include "caps/field.hpp"
#include "caps/kspace.hpp"
#include "kspace_pool.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <map>
#include <set>

#include "caps/grow.hpp"
#include "parallel.hpp"

namespace caps {

namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr double kDeg = kPi / 180.0;
constexpr double kCoulomb = 332.06371;   // kcal·Å/(mol·e²)

struct LjParam { const char* t; double rstar, eps, mass; };
constexpr LjParam kLj[] = {
    {"c3", 1.9080, 0.1094, 12.01}, {"ca", 1.9080, 0.0860, 12.01}, {"hc", 1.4870, 0.0157, 1.008}, {"ha", 1.4590, 0.0150, 1.008}};

struct BondParam { const char *a, *b; double k, r0; };
constexpr BondParam kBonds[] = {
    {"c3", "c3", 303.1, 1.535}, {"c3", "ca", 323.5, 1.513}, {"ca", "ca", 478.4, 1.387},
    {"c3", "hc", 337.3, 1.092}, {"ca", "ha", 344.3, 1.087}};

struct AngleParam { const char *a, *b, *c; double k, t0; };
constexpr AngleParam kAngles[] = {
    {"c3", "c3", "c3", 63.21, 110.63}, {"c3", "c3", "hc", 46.37, 110.05}, {"hc", "c3", "hc", 39.43, 108.35},
    {"c3", "c3", "ca", 63.26, 112.09}, {"ca", "c3", "hc", 46.96, 110.15}, {"c3", "ca", "ca", 63.84, 120.63},
    {"ca", "ca", "ca", 67.18, 119.97}, {"ca", "ca", "ha", 48.46, 120.01}};

// Proper torsions: specific entries replace the generic X-b-c-X entry for the same central pair.
struct TorsionParam { const char *a, *b, *c, *d; double v; int n; double delta; };
constexpr TorsionParam kTorsions[] = {
    {"X", "c3", "c3", "X", 1.40 / 9, 3, 0},
    {"hc", "c3", "c3", "hc", 0.15, 3, 0},
    {"c3", "c3", "c3", "hc", 0.16, 3, 0},
    {"c3", "c3", "c3", "c3", 0.18, 3, 0}, {"c3", "c3", "c3", "c3", 0.25, 2, 180}, {"c3", "c3", "c3", "c3", 0.20, 1, 180},
    {"X", "c3", "ca", "X", 0.0, 2, 0},
    {"X", "ca", "ca", "X", 14.50 / 4, 2, 180}};

// Out-of-plane term for every aromatic carbon (X-X-ca-ha, applied to all ca centres), AMBER order: centre third.
constexpr double kImproperV = 1.1;

const LjParam& lj_of(const std::string& t) {
  for (const auto& p : kLj)
    if (t == p.t) return p;
  throw FieldError("no GAFF LJ parameters for type " + t);
}

}  // namespace

ForceField assign_gaff(const System& s) {
  ForceField ff;
  ff.name = "GAFF 1.81 (hydrocarbon subset: c3, ca, hc, ha)";
  const size_t n = s.atoms.size();
  const auto nb = s.neighbours();

  // Aromatic carbons: three-connected carbons in a six-membered ring of three-connected carbons.
  std::vector<char> tri(n, 0), aromatic(n, 0);
  for (size_t i = 0; i < n; ++i) tri[i] = s.atoms[i].element == 6 && nb[i].size() == 3;
  for (size_t i = 0; i < n; ++i) {
    if (!tri[i]) continue;
    // depth-limited search for a 6-cycle through i
    std::vector<uint32_t> path{static_cast<uint32_t>(i)};
    bool found = false;
    std::function<void()> dfs = [&]() {
      if (found) return;
      const uint32_t u = path.back();
      for (uint32_t v : nb[u]) {
        if (!tri[v]) continue;
        if (path.size() == 6) {
          if (v == i) found = true;
          continue;
        }
        if (std::find(path.begin(), path.end(), v) != path.end()) continue;
        path.push_back(v);
        dfs();
        path.pop_back();
        if (found) return;
      }
    };
    dfs();
    aromatic[i] = found;
  }

  ff.atom_type.resize(n);
  ff.why.resize(n);
  for (size_t i = 0; i < n; ++i) {
    const auto& a = s.atoms[i];
    if (a.element == 6) {
      if (nb[i].size() == 4) { ff.atom_type[i] = "c3"; ff.why[i] = "sp3 carbon (4 neighbours)"; }
      else if (aromatic[i]) { ff.atom_type[i] = "ca"; ff.why[i] = "aromatic carbon (6-ring of 3-connected carbons)"; }
      else throw FieldError("atom " + std::to_string(a.id) + ": carbon with " + std::to_string(nb[i].size()) +
                            " neighbours is not sp3 or aromatic; only c3 / ca are parameterised in this version");
    } else if (a.element == 1) {
      if (nb[i].size() != 1) throw FieldError("atom " + std::to_string(a.id) + ": hydrogen with " + std::to_string(nb[i].size()) + " bonds");
    } else {
      throw FieldError("atom " + std::to_string(a.id) + ": element " + std::to_string(a.element) +
                       " has no parameters yet (this version types C and H only)");
    }
  }
  for (size_t i = 0; i < n; ++i) {
    if (s.atoms[i].element != 1) continue;
    const uint32_t h = nb[i][0];
    if (ff.atom_type[h] == "c3") { ff.atom_type[i] = "hc"; ff.why[i] = "hydrogen on sp3 carbon"; }
    else if (ff.atom_type[h] == "ca") { ff.atom_type[i] = "ha"; ff.why[i] = "hydrogen on aromatic carbon"; }
    else throw FieldError("atom " + std::to_string(s.atoms[i].id) + ": hydrogen not bonded to carbon");
  }

  // Types (in table order, so numbering is stable), LJ, masses.
  std::map<std::string, int> tix;
  for (const auto& p : kLj)
    if (std::find(ff.atom_type.begin(), ff.atom_type.end(), p.t) != ff.atom_type.end()) {
      tix[p.t] = static_cast<int>(ff.type_names.size());
      ff.type_names.push_back(p.t);
      ff.lj.push_back({p.eps, 2.0 * p.rstar / std::pow(2.0, 1.0 / 6.0)});
    }
  ff.type_index.resize(n);
  for (size_t i = 0; i < n; ++i) {
    ff.type_index[i] = tix.at(ff.atom_type[i]);
    ff.mass.push_back(lj_of(ff.atom_type[i]).mass);
  }

  // Charges.
  if (s.has_charges) {
    for (const auto& a : s.atoms) ff.charge.push_back(a.charge);
    ff.notes.push_back("charges kept from the structure");
  } else {
    ff.charge = gasteiger_ch(s, aromatic);
    ff.notes.push_back("Gasteiger–Marsili charges computed");
  }

  const auto& T = ff.atom_type;
  auto match = [](const char* p, const std::string& t) { return std::string(p) == "X" || t == p; };

  for (const auto& b : s.bonds) {
    bool ok = false;
    for (const auto& p : kBonds)
      if ((T[b.i] == p.a && T[b.j] == p.b) || (T[b.i] == p.b && T[b.j] == p.a)) {
        ff.bonds.push_back({b.i, b.j, p.k, p.r0});
        ok = true;
        break;
      }
    if (!ok) throw FieldError("no GAFF bond parameters for " + T[b.i] + "-" + T[b.j]);
  }

  std::set<std::pair<uint32_t, uint32_t>> ex12, ex13;
  for (const auto& b : s.bonds) ex12.insert({std::min(b.i, b.j), std::max(b.i, b.j)});
  for (uint32_t j = 0; j < n; ++j)
    for (size_t x = 0; x < nb[j].size(); ++x)
      for (size_t y = x + 1; y < nb[j].size(); ++y) {
        const uint32_t i = nb[j][x], k = nb[j][y];
        bool ok = false;
        for (const auto& p : kAngles)
          if (T[j] == p.b && ((T[i] == p.a && T[k] == p.c) || (T[i] == p.c && T[k] == p.a))) {
            ff.angles.push_back({i, j, k, p.k, p.t0 * kDeg});
            ok = true;
            break;
          }
        if (!ok) throw FieldError("no GAFF angle parameters for " + T[i] + "-" + T[j] + "-" + T[k]);
        ex13.insert({std::min(i, k), std::max(i, k)});
      }

  std::set<std::pair<uint32_t, uint32_t>> p14;
  for (const auto& b : s.bonds) {
    const uint32_t j = b.i, k = b.j;
    for (uint32_t i : nb[j]) {
      if (i == k) continue;
      for (uint32_t l : nb[k]) {
        if (l == j || l == i) continue;
        // specific entries for this (a, b, c, d) first, the generic X-b-c-X entry otherwise
        auto fits = [&](const TorsionParam& p, bool generic) {
          if ((std::string(p.a) == "X") != generic) return false;
          return (match(p.a, T[i]) && T[j] == p.b && T[k] == p.c && match(p.d, T[l])) ||
                 (match(p.a, T[l]) && T[k] == p.b && T[j] == p.c && match(p.d, T[i]));
        };
        bool any = false;
        for (const auto& p : kTorsions)
          if (fits(p, false)) { any = true; if (p.v != 0) ff.dihedrals.push_back({i, j, k, l, p.v, p.n, p.delta * kDeg}); }
        if (!any) {
          for (const auto& p : kTorsions)
            if (fits(p, true)) { any = true; if (p.v != 0) ff.dihedrals.push_back({i, j, k, l, p.v, p.n, p.delta * kDeg}); }
        }
        if (!any) throw FieldError("no GAFF torsion parameters for " + T[i] + "-" + T[j] + "-" + T[k] + "-" + T[l]);
        const auto key = std::make_pair(std::min(i, l), std::max(i, l));
        if (!ex12.count(key) && !ex13.count(key)) p14.insert(key);
      }
    }
  }
  for (const auto& p : p14) ff.pairs14.push_back({p.first, p.second});

  for (uint32_t c = 0; c < n; ++c) {
    if (T[c] != "ca" || nb[c].size() != 3) continue;
    // the substituent (not in the ring) is the last atom; ring neighbours first
    std::vector<uint32_t> ring, other;
    for (uint32_t v : nb[c]) (T[v] == "ca" ? ring : other).push_back(v);
    if (ring.size() == 3) { other.push_back(ring.back()); ring.pop_back(); }
    if (ring.size() != 2 || other.size() != 1) continue;
    ff.impropers.push_back({ring[0], ring[1], c, other[0], kImproperV, 2, 180 * kDeg});
  }

  ff.excluded.assign(n, {});
  auto add_ex = [&](uint32_t a, uint32_t b) { ff.excluded[a].push_back(b); ff.excluded[b].push_back(a); };
  for (const auto& p : ex12) add_ex(p.first, p.second);
  for (const auto& p : ex13) if (!ex12.count(p)) add_ex(p.first, p.second);
  for (const auto& p : p14) add_ex(p.first, p.second);
  for (auto& e : ff.excluded) {
    std::sort(e.begin(), e.end());
    e.erase(std::unique(e.begin(), e.end()), e.end());
  }

  std::map<std::string, int> counts;
  for (const auto& t : T) counts[t]++;
  std::string c;
  for (const auto& [t, k] : counts) c += (c.empty() ? "" : ", ") + std::to_string(k) + " " + t;
  ff.notes.push_back("typed " + std::to_string(n) + " atoms: " + c);
  ff.notes.push_back(std::to_string(ff.bonds.size()) + " bonds, " + std::to_string(ff.angles.size()) + " angles, " +
                     std::to_string(ff.dihedrals.size()) + " torsion terms, " + std::to_string(ff.impropers.size()) + " impropers, " +
                     std::to_string(ff.pairs14.size()) + " 1-4 pairs");
  return ff;
}

// ---------------------------------------------------------------------------------------------------------------
// Evaluator

Evaluator::Evaluator(const ForceField& ff, const EnergyOptions& o)
    : ff_(ff), opt_(o), pool_(std::make_unique<ThreadPool>(o.threads > 0 ? o.threads : default_threads())) {
  std::set<std::pair<uint32_t, uint32_t>> p14;
  for (const auto& p : ff.pairs14) p14.insert({p[0], p[1]});
  for (uint32_t i = 0; i < ff.excluded.size(); ++i)
    for (uint32_t j : ff.excluded[i])
      if (j > i) excl_.push_back({i, j, p14.count({i, j}) ? ff.coul14 : 0.0});
  const size_t nt = ff.lj.size();
  eps_.resize(nt * nt);
  s6_.resize(nt * nt);
  for (size_t a = 0; a < nt; ++a)
    for (size_t b = 0; b < nt; ++b) {
      const PairType pt = mixed_pair(ff, int(a), int(b));
      eps_[a * nt + b] = pt.eps;
      s6_[a * nt + b] = std::pow(pt.sigma, 6);
    }
  form_.assign(nt * nt, 0);
  pa_.assign(nt * nt, 0.0);
  pb_.assign(nt * nt, 0.0);
  pc_.assign(nt * nt, 0.0);
  for (const auto& [ab, pf] : ff.pair_func) {
    for (const size_t tp : {size_t(ab.first) * nt + ab.second, size_t(ab.second) * nt + ab.first}) {
      form_[tp] = uint8_t(pf.form);
      pa_[tp] = pf.a;
      pb_[tp] = pf.b;
      pc_[tp] = pf.c;
      eps_[tp] = 1;   // non-zero: the pair is evaluated
    }
  }
  if (!ff.lj14_types.empty()) {
    if (ff.lj14_types.size() != nt) throw FieldError("lj14_types needs one entry per type");
    ForceField f14;
    f14.mixing = ff.mixing;
    f14.lj = ff.lj14_types;
    eps14_.resize(nt * nt);
    s614_.resize(nt * nt);
    for (size_t a = 0; a < nt; ++a)
      for (size_t b = 0; b < nt; ++b) {
        const PairType pt = mixed_pair(f14, int(a), int(b));
        eps14_[a * nt + b] = pt.eps;
        s614_[a * nt + b] = std::pow(pt.sigma, 6);
      }
  }
  type_count_.assign(ff.lj.size(), 0.0);
  for (int t : ff.type_index) type_count_[t] += 1;
  if (ff.pair_form != "lj12-6" && ff.pair_form != "lj9-6") throw FieldError("unknown pair form '" + ff.pair_form + "'");
  lj96_ = ff.pair_form == "lj9-6";
  cap_radii();
}

Evaluator::~Evaluator() = default;

int Evaluator::threads() const { return pool_->size(); }

void Evaluator::set_options(const EnergyOptions& o) {
  if (o.threads != opt_.threads) pool_ = std::make_unique<ThreadPool>(o.threads > 0 ? o.threads : default_threads());
  const bool relist = o.cutoff != opt_.cutoff || o.skin != opt_.skin;
  opt_ = o;
  if (relist) built_ = false;
  cap_radii();
}

// For each type pair, the radius inside which |F_LJ| exceeds the cap; the potential is continued linearly inside it.
void Evaluator::cap_radii() {
  const size_t nt2 = eps_.size();
  rcap2_.assign(nt2, 0.0);
  ecap_.assign(nt2, 0.0);
  fcap_.assign(nt2, 0.0);
  if (opt_.force_cap <= 0) return;
  for (size_t t = 0; t < nt2; ++t) {
    const double e = eps_[t], s6 = s6_[t];
    if (e <= 0 || form_[t] != 0) continue;   // the cap is for Lennard-Jones only
    const double sg = std::pow(s6, 1.0 / 6);
    auto force = [&](double r) {
      if (lj96_) { const double q = sg / r, q3 = q * q * q; return 18 * e / r * (q3 * q3 * q3 - q3 * q3); }
      const double q = s6 / std::pow(r, 6);
      return 24 * e / r * (2 * q * q - q);
    };
    double lo = 0.3 * sg, hi = lj96_ ? sg : std::pow(2 * s6, 1.0 / 6);   // F(hi) = 0
    if (force(lo) < opt_.force_cap) continue;                                   // cap never reached
    for (int it = 0; it < 80; ++it) {
      const double m = 0.5 * (lo + hi);
      (force(m) > opt_.force_cap ? lo : hi) = m;
    }
    const double r = hi, q = s6 / std::pow(r, 6);
    rcap2_[t] = r * r;
    ecap_[t] = lj96_ ? e * (2 * std::pow(sg / r, 9) - 3 * q) : 4 * e * (q * q - q);
    fcap_[t] = force(r);
  }
}

namespace {

// erfc(x) with exp(−x²) returned too: Abramowitz and Stegun 7.1.26 (|error| < 1.5e-7), the form LAMMPS coul/dsf uses.
inline double erfc_exp(double x, double& ex2) {
  constexpr double P = 0.3275911, A1 = 0.254829592, A2 = -0.284496736, A3 = 1.421413741, A4 = -1.453152027, A5 = 1.061405429;
  ex2 = std::exp(-x * x);
  const double t = 1.0 / (1.0 + P * x);
  return t * (A1 + t * (A2 + t * (A3 + t * (A4 + t * A5)))) * ex2;
}

// ---- class II terms ----------------------------------------------------------------------------------------------
// Each returns the energy, the force on each atom and Σ r·f. Internal coordinates and their gradients:
// a bond length, an angle at a vertex, the IUPAC dihedral φ (signed, as LAMMPS dihedral class2) and the Wilson
// out-of-plane angle χ of improper class2.

struct Class2Out {
  double e = 0;
  Vec3 f[4] = {};
  double virial = 0;
  double w[6] = {0, 0, 0, 0, 0, 0};   // virial tensor xx yy zz xy xz yz
  // adds the virial of force f acting at relative position d
  void vir(const Vec3& d, const Vec3& fv) {
    virial += dot(d, fv);
    w[0] += d[0] * fv[0]; w[1] += d[1] * fv[1]; w[2] += d[2] * fv[2];
    w[3] += 0.5 * (d[0] * fv[1] + d[1] * fv[0]); w[4] += 0.5 * (d[0] * fv[2] + d[2] * fv[0]); w[5] += 0.5 * (d[1] * fv[2] + d[2] * fv[1]);
  }
};

// θ between u and v (both from the vertex) and dθ/du, dθ/dv.
inline double angle_grad(const Vec3& u, const Vec3& v, Vec3& gu, Vec3& gv) {
  const double lu = norm(u), lv = norm(v);
  const double c = std::clamp(dot(u, v) / (lu * lv), -1.0, 1.0);
  const double sn = std::max(std::sqrt(1 - c * c), 1e-8);
  gu = (v * (1 / (lu * lv)) - u * (c / (lu * lu))) * (-1 / sn);
  gv = (u * (1 / (lu * lv)) - v * (c / (lv * lv))) * (-1 / sn);
  return std::acos(c);
}

// φ of i-j-k-l from b1 = x_j − x_i, b2 = x_k − x_j, b3 = x_l − x_k, and dφ/dx for the four atoms.
inline double dihedral_grad(const Vec3& b1, const Vec3& b2, const Vec3& b3, Vec3 g[4]) {
  const Vec3 m = cross(b1, b2), nn = cross(b2, b3);
  const double m2 = std::max(dot(m, m), 1e-24), n2 = std::max(dot(nn, nn), 1e-24), lb2 = norm(b2);
  g[0] = m * (-lb2 / m2);
  g[3] = nn * (lb2 / n2);
  const double p = dot(b1, b2) / (lb2 * lb2), q = dot(b3, b2) / (lb2 * lb2);
  g[1] = g[3] * q - g[0] * (1 + p);
  g[2] = g[0] * p - g[3] * (1 + q);
  return std::atan2(lb2 * dot(b1, nn), dot(m, nn));
}

// Wilson angle χ between the plane of a, b and the vector c (all from the centre): sin χ = (a × b)·c / (|a × b| |c|).
inline double wilson_grad(const Vec3& a, const Vec3& b, const Vec3& c, Vec3& ga, Vec3& gb, Vec3& gc) {
  const Vec3 p = cross(a, b);
  const double P = norm(p), C = norm(c), t = dot(p, c);
  const double s = std::clamp(t / (P * C), -1.0, 1.0);
  const double inv = 1 / std::max(std::sqrt(1 - s * s), 1e-8);
  ga = (cross(b, c) * (1 / (P * C)) - cross(b, p) * (s / (P * P))) * inv;
  gb = (cross(c, a) * (1 / (P * C)) - cross(p, a) * (s / (P * P))) * inv;
  gc = (p * (1 / (P * C)) - c * (s / (C * C))) * inv;
  return std::asin(s);
}

// d = x_j − x_i
inline Class2Out class2_bond(const Class2Bond& t, const Vec3& d) {
  Class2Out o;
  const double r = norm(d), dr = r - t.r0, dr2 = dr * dr;
  o.e = t.k2 * dr2 + t.k3 * dr2 * dr + t.k4 * dr2 * dr2;
  const double de = 2 * t.k2 * dr + 3 * t.k3 * dr2 + 4 * t.k4 * dr2 * dr;
  o.f[1] = r > 1e-12 ? d * (-de / r) : Vec3{0, 0, 0};   // a core and its shell on the same point: no force
  o.f[0] = o.f[1] * -1.0;
  o.vir(d, o.f[1]);
  return o;
}

// u = x_i − x_j, v = x_k − x_j
inline Class2Out class2_angle(const Class2Angle& t, const Vec3& u, const Vec3& v) {
  Class2Out o;
  Vec3 gu, gv;
  const double th = angle_grad(u, v, gu, gv), r1 = norm(u), r2 = norm(v);
  const double dt = th - t.theta0, dt2 = dt * dt, d1 = r1 - t.bb_r1, d2 = r2 - t.bb_r2, e1 = r1 - t.ba_r1, e2 = r2 - t.ba_r2;
  o.e = t.k2 * dt2 + t.k3 * dt2 * dt + t.k4 * dt2 * dt2 + t.bb_m * d1 * d2 + (t.ba_n1 * e1 + t.ba_n2 * e2) * dt;
  const double dEdt = 2 * t.k2 * dt + 3 * t.k3 * dt2 + 4 * t.k4 * dt2 * dt + t.ba_n1 * e1 + t.ba_n2 * e2;
  const double dEdr1 = t.bb_m * d2 + t.ba_n1 * dt, dEdr2 = t.bb_m * d1 + t.ba_n2 * dt;
  o.f[0] = (gu * dEdt + u * (dEdr1 / r1)) * -1.0;
  o.f[2] = (gv * dEdt + v * (dEdr2 / r2)) * -1.0;
  o.f[1] = (o.f[0] + o.f[2]) * -1.0;
  o.vir(u, o.f[0]); o.vir(v, o.f[2]);
  return o;
}

// b1 = x_j − x_i, b2 = x_k − x_j, b3 = x_l − x_k
inline Class2Out class2_dihedral(const Class2Dihedral& t, const Vec3& b1, const Vec3& b2, const Vec3& b3) {
  Class2Out o;
  Vec3 gp[4], ga1, gc1, ga2, gc2;
  const double phi = dihedral_grad(b1, b2, b3, gp);
  const double th1 = angle_grad(b1 * -1.0, b2, ga1, gc1);   // i-j-k: vectors j→i, j→k
  const double th2 = angle_grad(b2 * -1.0, b3, ga2, gc2);   // j-k-l: vectors k→j, k→l
  const double r1 = norm(b1), r2 = norm(b2), r3 = norm(b3);
  const double c1 = std::cos(phi), c2 = std::cos(2 * phi), c3 = std::cos(3 * phi);
  const double s1 = std::sin(phi), s2 = std::sin(2 * phi), s3 = std::sin(3 * phi);
  const double dr1 = r1 - t.ebt_r1, dr2 = r2 - t.mbt_r2, dr3 = r3 - t.ebt_r3;
  const double dt1 = th1 - t.at_theta1, dt2 = th2 - t.at_theta2;
  const double da1 = th1 - t.aat_theta1, da2 = th2 - t.aat_theta2;
  const double db1 = r1 - t.bb13_r1, db3 = r3 - t.bb13_r3;
  const double A = t.mbt[0] * c1 + t.mbt[1] * c2 + t.mbt[2] * c3, dA = -(t.mbt[0] * s1 + 2 * t.mbt[1] * s2 + 3 * t.mbt[2] * s3);
  const double B = t.ebt_b[0] * c1 + t.ebt_b[1] * c2 + t.ebt_b[2] * c3, dB = -(t.ebt_b[0] * s1 + 2 * t.ebt_b[1] * s2 + 3 * t.ebt_b[2] * s3);
  const double C = t.ebt_c[0] * c1 + t.ebt_c[1] * c2 + t.ebt_c[2] * c3, dC = -(t.ebt_c[0] * s1 + 2 * t.ebt_c[1] * s2 + 3 * t.ebt_c[2] * s3);
  const double D = t.at_d[0] * c1 + t.at_d[1] * c2 + t.at_d[2] * c3, dD = -(t.at_d[0] * s1 + 2 * t.at_d[1] * s2 + 3 * t.at_d[2] * s3);
  const double E = t.at_e[0] * c1 + t.at_e[1] * c2 + t.at_e[2] * c3, dE = -(t.at_e[0] * s1 + 2 * t.at_e[1] * s2 + 3 * t.at_e[2] * s3);
  o.e = t.k1 * (1 - std::cos(phi - t.phi1)) + t.k2 * (1 - std::cos(2 * phi - t.phi2)) + t.k3 * (1 - std::cos(3 * phi - t.phi3)) +
        dr2 * A + dr1 * B + dr3 * C + dt1 * D + dt2 * E + t.aat_m * da1 * da2 * c1 + t.bb13_n * db1 * db3;
  const double dEdphi = t.k1 * std::sin(phi - t.phi1) + 2 * t.k2 * std::sin(2 * phi - t.phi2) + 3 * t.k3 * std::sin(3 * phi - t.phi3) +
                        dr2 * dA + dr1 * dB + dr3 * dC + dt1 * dD + dt2 * dE - t.aat_m * da1 * da2 * s1;
  const double dEdr1 = B + t.bb13_n * db3, dEdr2 = A, dEdr3 = C + t.bb13_n * db1;
  const double dEdth1 = D + t.aat_m * da2 * c1, dEdth2 = E + t.aat_m * da1 * c1;
  Vec3 g[4];
  for (int a = 0; a < 4; ++a) g[a] = gp[a] * dEdphi;
  // bond lengths: r1 = |x_j − x_i|, r2 = |x_k − x_j|, r3 = |x_l − x_k|
  const Vec3 u1 = b1 * (dEdr1 / r1), u2 = b2 * (dEdr2 / r2), u3 = b3 * (dEdr3 / r3);
  g[0] = g[0] - u1; g[1] = g[1] + u1 - u2; g[2] = g[2] + u2 - u3; g[3] = g[3] + u3;
  // θ1 at j (i-j-k), θ2 at k (j-k-l)
  g[0] = g[0] + ga1 * dEdth1; g[2] = g[2] + gc1 * dEdth1; g[1] = g[1] - (ga1 + gc1) * dEdth1;
  g[1] = g[1] + ga2 * dEdth2; g[3] = g[3] + gc2 * dEdth2; g[2] = g[2] - (ga2 + gc2) * dEdth2;
  for (int a = 0; a < 4; ++a) o.f[a] = g[a] * -1.0;
  o.vir(b1, o.f[1]); o.vir(b1 + b2, o.f[2]); o.vir(b1 + b2 + b3, o.f[3]);
  return o;
}

// a = x_i − x_j, c = x_k − x_j, d = x_l − x_j (j the centre; LAMMPS A, B, C, D = i, j, k, l)
inline Class2Out class2_improper(const Class2Improper& t, const Vec3& a, const Vec3& c, const Vec3& d) {
  Class2Out o;
  Vec3 ga{0, 0, 0}, gc{0, 0, 0}, gd{0, 0, 0};
  if (t.kchi != 0) {
    Vec3 x1, x2, x3;
    double chi = 0;
    chi += wilson_grad(c, d, a, x1, x2, x3); gc = gc + x1; gd = gd + x2; ga = ga + x3;   // plane C-B-D, bond B-A
    chi += wilson_grad(d, a, c, x1, x2, x3); gd = gd + x1; ga = ga + x2; gc = gc + x3;   // plane D-B-A, bond B-C
    chi += wilson_grad(a, c, d, x1, x2, x3); ga = ga + x1; gc = gc + x2; gd = gd + x3;   // plane A-B-C, bond B-D
    chi /= 3;
    const double dc = chi - t.chi0;
    o.e = t.kchi * dc * dc;
    const double w = 2 * t.kchi * dc / 3;
    ga = ga * w; gc = gc * w; gd = gd * w;
  }
  Vec3 p1, p2, q1, q2, s1, s2;
  const double tABC = angle_grad(a, c, p1, p2), tABD = angle_grad(a, d, q1, q2), tCBD = angle_grad(c, d, s1, s2);
  const double d1 = tABC - t.theta1, d2 = tABD - t.theta2, d3 = tCBD - t.theta3;
  o.e += t.m1 * d1 * d3 + t.m2 * d1 * d2 + t.m3 * d2 * d3;
  const double e1 = t.m1 * d3 + t.m2 * d2, e2 = t.m2 * d1 + t.m3 * d3, e3 = t.m1 * d1 + t.m3 * d2;
  ga = ga + p1 * e1 + q1 * e2;
  gc = gc + p2 * e1 + s1 * e3;
  gd = gd + q2 * e2 + s2 * e3;
  o.f[0] = ga * -1.0; o.f[2] = gc * -1.0; o.f[3] = gd * -1.0;
  o.f[1] = (o.f[0] + o.f[2] + o.f[3]) * -1.0;
  o.vir(a, o.f[0]); o.vir(c, o.f[2]); o.vir(d, o.f[3]);
  return o;
}

// a, b, d: bonds from the centre. Returns forces in the order a, centre, b, d.
inline Class2Out inversion(const InversionTerm& t, const Vec3& a, const Vec3& b, const Vec3& d) {
  Class2Out o;
  Vec3 ga{0, 0, 0}, gb{0, 0, 0}, gd{0, 0, 0};
  const Vec3* v[3] = {&a, &b, &d};
  Vec3* g[3] = {&ga, &gb, &gd};
  for (int q = 0; q < 3; ++q) {
    // bond q against the plane of the other two
    const int p1 = (q + 1) % 3, p2 = (q + 2) % 3;
    Vec3 x1, x2, x3;
    const double chi = wilson_grad(*v[p1], *v[p2], *v[q], x1, x2, x3);
    const double w = std::fabs(chi), sg = chi < 0 ? -1.0 : 1.0, dw = w - t.w0;
    double c;
    if (t.form == 1) {   // 1 − cos ω
      o.e += t.kw / 3 * (1 - std::cos(w));
      c = t.kw / 3 * std::sin(w) * sg;
    } else if (t.form == 2) {   // UFF pyramidal centres: C0 + C1 cos ω + C2 cos 2ω, minimum at ω0 (LAMMPS improper fourier)
      const double C2 = 1, C1 = -4 * std::cos(t.w0), C0 = -(C1 * std::cos(t.w0) + C2 * std::cos(2 * t.w0));
      o.e += t.kw / 3 * (C0 + C1 * std::cos(w) + C2 * std::cos(2 * w));
      c = -t.kw / 3 * (C1 * std::sin(w) + 2 * C2 * std::sin(2 * w)) * sg;
    } else {
      o.e += t.kw / 3 * dw * dw;
      c = 2 * t.kw / 3 * dw * sg;
    }
    *g[p1] = *g[p1] + x1 * c; *g[p2] = *g[p2] + x2 * c; *g[q] = *g[q] + x3 * c;
  }
  o.f[0] = ga * -1.0; o.f[2] = gb * -1.0; o.f[3] = gd * -1.0;
  o.f[1] = (o.f[0] + o.f[2] + o.f[3]) * -1.0;
  o.vir(a, o.f[0]); o.vir(b, o.f[2]); o.vir(d, o.f[3]);
  return o;
}

inline bool orthogonal(const Cell& c) {
  return c.a[1] == 0 && c.a[2] == 0 && c.b[0] == 0 && c.b[2] == 0 && c.c[0] == 0 && c.c[1] == 0;
}

}  // namespace

void Evaluator::build(const std::vector<double>& x, const Cell& cell_in) {
  const size_t n = x.size() / 3;
  const double rcs = opt_.cutoff + opt_.skin;
  Cell cell = cell_in;
  if (!cell.valid()) {
    Vec3 lo{1e300, 1e300, 1e300}, hi{-1e300, -1e300, -1e300};
    for (size_t i = 0; i < n; ++i)
      for (int k = 0; k < 3; ++k) { lo[k] = std::min(lo[k], x[3 * i + k]); hi[k] = std::max(hi[k], x[3 * i + k]); }
    if (n == 0) lo = hi = {0, 0, 0};
    cell.origin = lo - Vec3{1, 1, 1};
    cell.a = {hi[0] - lo[0] + 2, 0, 0};
    cell.b = {0, hi[1] - lo[1] + 2, 0};
    cell.c = {0, 0, hi[2] - lo[2] + 2};
    cell.periodic = {false, false, false};
  }
  const double v = cell.volume();
  const double w[3] = {v / norm(cross(cell.b, cell.c)), v / norm(cross(cell.c, cell.a)), v / norm(cross(cell.a, cell.b))};
  int nb[3], reach[3];
  for (int k = 0; k < 3; ++k) {
    // bins of about half the list radius: a 5×5×5 stencil scans 15.6 (rc + skin)³ instead of 27 (rc + skin)³
    nb[k] = std::max(1, std::min(512, static_cast<int>(2 * w[k] / rcs)));
    reach[k] = cell.periodic[k] ? static_cast<int>(std::ceil(rcs / (w[k] / nb[k]))) : 1;
  }
  std::vector<Vec3> fr(n);
  std::vector<std::array<int, 3>> bin(n);
  std::vector<std::vector<uint32_t>> bins(static_cast<size_t>(nb[0]) * nb[1] * nb[2]);
  auto index = [&](int a, int b, int c) { return (static_cast<size_t>(a) * nb[1] + b) * nb[2] + c; };
  for (size_t i = 0; i < n; ++i) {
    Vec3 f = cell.to_fractional({x[3 * i], x[3 * i + 1], x[3 * i + 2]});
    for (int k = 0; k < 3; ++k) {
      if (cell.periodic[k]) f[k] -= std::floor(f[k]);
      f[k] = std::clamp(f[k], 0.0, 1.0 - 1e-12);
      bin[i][k] = std::min(nb[k] - 1, static_cast<int>(f[k] * nb[k]));
    }
    fr[i] = f;
    bins[index(bin[i][0], bin[i][1], bin[i][2])].push_back(static_cast<uint32_t>(i));
  }
  const double rcs2 = rcs * rcs;
  const int nth = pool_->size();
  std::vector<std::vector<uint32_t>> vi(nth), vj(nth);
  std::vector<std::vector<double>> vs(nth);
  pool_->run(n, [&](int t, size_t b, size_t en) {
    auto& PI = vi[t];
    auto& PJ = vj[t];
    auto& SH = vs[t];
    PI.clear(); PJ.clear(); SH.clear();
    for (size_t i = b; i < en; ++i) {
      const auto& ex = ff_.excluded[i];
      for (int dx = -reach[0]; dx <= reach[0]; ++dx)
        for (int dy = -reach[1]; dy <= reach[1]; ++dy)
          for (int dz = -reach[2]; dz <= reach[2]; ++dz) {
            int u[3] = {bin[i][0] + dx, bin[i][1] + dy, bin[i][2] + dz}, q[3] = {0, 0, 0};
            bool skip = false;
            for (int k = 0; k < 3; ++k) {
              if (u[k] < 0 || u[k] >= nb[k]) {
                if (!cell.periodic[k]) { skip = true; break; }
                q[k] = static_cast<int>(std::floor(static_cast<double>(u[k]) / nb[k]));
                u[k] -= q[k] * nb[k];
              }
            }
            if (skip) continue;
            const bool zero = q[0] == 0 && q[1] == 0 && q[2] == 0;
            for (uint32_t j : bins[index(u[0], u[1], u[2])]) {
              // Half list balanced over atoms: an unordered pair is kept from the side given by the parity of i + j.
              if (j != i && (j > i) != (((i + j) & 1) == 0)) continue;
              if (j == i) {
                // self images: keep one of each ±shift pair
                if (zero) continue;
                if (q[0] < 0 || (q[0] == 0 && (q[1] < 0 || (q[1] == 0 && q[2] < 0)))) continue;
              }
              const Vec3 df = fr[j] + Vec3{double(q[0]), double(q[1]), double(q[2])} - fr[i];
              const Vec3 d = cell.a * df[0] + cell.b * df[1] + cell.c * df[2];
              if (dot(d, d) >= rcs2) continue;
              // 1-2, 1-3 and 1-4 partners are excluded at their bonded (nearest) image only; further images interact
              if (j != i && dot(d, d) < 36.0 && std::binary_search(ex.begin(), ex.end(), j)) continue;
              PI.push_back(static_cast<uint32_t>(i));
              PJ.push_back(j);
              // d = x_j − x_i + shift with the actual (possibly unwrapped) positions
              for (int k = 0; k < 3; ++k) SH.push_back(d[k] - (x[3 * j + k] - x[3 * i + k]));
            }
          }
    }
  });
  pi_.clear();
  pj_.clear();
  shift_.clear();
  for (int t = 0; t < nth; ++t) {
    pi_.insert(pi_.end(), vi[t].begin(), vi[t].end());
    pj_.insert(pj_.end(), vj[t].begin(), vj[t].end());
    shift_.insert(shift_.end(), vs[t].begin(), vs[t].end());
  }
  // Each shift is a lattice vector; keep its integer coefficients so a change of cell (barostat, compression)
  // only rescales the shifts instead of forcing a rebuild.
  ishift_.assign(shift_.size(), 0);
  if (cell_in.valid())
    for (size_t p = 0; p < shift_.size(); p += 3) {
      const Vec3 fr0 = cell_in.to_fractional(cell_in.origin + Vec3{shift_[p], shift_[p + 1], shift_[p + 2]});
      for (int k = 0; k < 3; ++k) ishift_[p + k] = static_cast<int32_t>(std::lround(fr0[k]));
    }
  x0_ = x;
  cell0_ = cell_in;
  cells_ = cell_in;
  built_ = true;
  ++builds_;
}

EnergyTerms Evaluator::compute(const std::vector<double>& x, const Cell& cell, std::vector<double>& f) {
  const size_t n = x.size() / 3;
  f.assign(x.size(), 0.0);
  EnergyTerms e;

  // r-RESPA splits the forces: bonded terms alone need no pair list
  const bool nonb = opt_.parts & 2, bonded = opt_.parts & 1;
  if (nonb) {
    auto same = [](const Cell& p, const Cell& q) { return p.a == q.a && p.b == q.b && p.c == q.c && p.origin == q.origin; };
    bool rebuild = !built_ || x0_.size() != x.size();
    if (frozen_ && !rebuild) {
      if (!same(cell, cells_) && cell.valid()) {
        for (size_t p = 0; p < shift_.size(); p += 3) {
          const Vec3 v = cell.a * ishift_[p] + cell.b * ishift_[p + 1] + cell.c * ishift_[p + 2];
          shift_[p] = v[0]; shift_[p + 1] = v[1]; shift_[p + 2] = v[2];
        }
        cells_ = cell;
      }
    } else if (!rebuild) {
      if (same(cell, cell0_)) {
        const double lim = 0.25 * opt_.skin * opt_.skin;
        for (size_t i = 0; i < n && !rebuild; ++i) {
          const double dx = x[3 * i] - x0_[3 * i], dy = x[3 * i + 1] - x0_[3 * i + 1], dz = x[3 * i + 2] - x0_[3 * i + 2];
          rebuild = dx * dx + dy * dy + dz * dz > lim;
        }
      } else if (!cell.valid() || !cell0_.valid()) {
        rebuild = true;
      } else {
        // The cell was deformed since the build. Shrinking eats into the skin by (1 − s)(rc + skin), s the smallest
        // ratio of perpendicular widths; the rest of the skin covers moves relative to the affinely mapped positions.
        auto widths = [](const Cell& c) {
          const double v = c.volume();
          return Vec3{v / norm(cross(c.b, c.c)), v / norm(cross(c.c, c.a)), v / norm(cross(c.a, c.b))};
        };
        const Vec3 w0 = widths(cell0_), w1 = widths(cell);
        const double sr = std::min({w1[0] / w0[0], w1[1] / w0[1], w1[2] / w0[2], 1.0});
        const double allow = opt_.skin - (1 - sr) * (opt_.cutoff + opt_.skin);
        if (allow < 0.1 * opt_.skin) rebuild = true;
        const double lim = 0.25 * allow * allow;
        for (size_t i = 0; i < n && !rebuild; ++i) {
          const Vec3 m = cell.to_cartesian(cell0_.to_fractional({x0_[3 * i], x0_[3 * i + 1], x0_[3 * i + 2]}));
          const double dx = x[3 * i] - m[0], dy = x[3 * i + 1] - m[1], dz = x[3 * i + 2] - m[2];
          rebuild = dx * dx + dy * dy + dz * dz > lim;
        }
        if (!rebuild && !same(cell, cells_)) {
          for (size_t p = 0; p < shift_.size(); p += 3) {
            const Vec3 v = cell.a * ishift_[p] + cell.b * ishift_[p + 1] + cell.c * ishift_[p + 2];
            shift_[p] = v[0]; shift_[p + 1] = v[1]; shift_[p + 2] = v[2];
          }
          cells_ = cell;
        }
      }
    }
    if (rebuild) {
      if (frozen_) throw std::logic_error("frozen pair list: the configuration changed size");
      build(x, cell);
    }
  }

  const bool ortho = cell.valid() && orthogonal(cell);
  const double L[3] = {cell.a[0], cell.b[1], cell.c[2]};
  auto mi = [&](Vec3 d) {
    if (!cell.valid()) return d;
    if (ortho) {
      for (int k = 0; k < 3; ++k)
        if (cell.periodic[k]) d[k] -= L[k] * std::round(d[k] / L[k]);
      return d;
    }
    return cell.minimum_image(d);
  };
  auto pos = [&](uint32_t i) { return Vec3{x[3 * i], x[3 * i + 1], x[3 * i + 2]}; };

  // Non-bonded constants
  const size_t nt = ff_.lj.size();
  const double rc = opt_.cutoff, rc2 = rc * rc;
  const double a = opt_.dsf_alpha;
  double ex_rc;
  const double erfc_rc = erfc_exp(a * rc, ex_rc);
  const double dsf_e0 = erfc_rc / rc;
  const double dsf_f0 = erfc_rc / (rc * rc) + 2 * a / std::sqrt(kPi) * std::exp(-a * a * rc2) / rc;
  const bool coul = opt_.coulomb;
  // PME: the pair term is erfc(βr)/r with no shift; the reciprocal part, self energy and background come after
  const bool pme = coul && opt_.electrostatics == EnergyOptions::Electrostatics::PME && cell.valid() && cell.periodic[0] && cell.periodic[1] &&
                   cell.periodic[2];
  const double beta = pme ? ewald_beta(rc, opt_.ewald_rtol) : 0.0;
  const double b2pi = 2 * beta / std::sqrt(kPi);
  const double a2pi = 2 * a / std::sqrt(kPi);
  const bool capped = opt_.force_cap > 0;

  auto lj = [&](size_t tp, double r2, double scale, double& en) -> double {   // returns −(dE/dr)/r
    const double eps = eps_[tp];
    if (eps <= 0) return 0.0;
    if (form_[tp] == 1) {   // Buckingham
      const double r = std::sqrt(r2), ex = pa_[tp] * std::exp(-r / pb_[tp]), c6 = pc_[tp] / (r2 * r2 * r2);
      en += scale * (ex - c6);
      return scale * (ex / pb_[tp] - 6 * c6 / r) / r;
    }
    if (form_[tp] == 2) {   // Morse
      const double r = std::sqrt(r2), e1 = std::exp(-pb_[tp] * (r - pc_[tp]));
      en += scale * pa_[tp] * (e1 * e1 - 2 * e1);
      return scale * 2 * pa_[tp] * pb_[tp] * (e1 * e1 - e1) / r;
    }
    if (capped && r2 < rcap2_[tp]) {
      const double r = std::sqrt(r2), rcp = std::sqrt(rcap2_[tp]);
      en += scale * (ecap_[tp] + fcap_[tp] * (rcp - r));
      return scale * fcap_[tp] / r;
    }
    const double q = s6_[tp] / (r2 * r2 * r2);
    if (lj96_) {
      const double q9 = q * std::sqrt(q);   // (σ/r)⁹
      en += scale * eps * (2 * q9 - 3 * q);
      return scale * 18 * eps * (q9 - q) / r2;
    }
    en += scale * 4 * eps * (q * q - q);
    return scale * 24 * eps * (2 * q * q - q) / r2;
  };

  // Energy shift so that LJ is zero at the cut-off.
  // With the tail correction the LJ energy is truncated, not shifted (the tail term assumes the plain potential).
  auto lj_shift = [&](size_t tp) {
    if (opt_.tail) return 0.0;
    if (form_[tp] != 0) {
      double e0 = 0;
      lj(tp, rc2, 1.0, e0);
      return eps_[tp] > 0 ? e0 : 0.0;
    }
    const double q = s6_[tp] / (rc2 * rc2 * rc2);
    if (lj96_) return eps_[tp] * (2 * q * std::sqrt(q) - 3 * q);
    return 4 * eps_[tp] * (q * q - q);
  };

  const double* q = ff_.charge.data();
  const int nth = pool_->size();
  tf_.resize(nth);
  // per worker: bond, angle, dihedral, improper, vdW, Coulomb, virial, virial tensor (xx yy zz xy xz yz)
  std::vector<std::array<double, 13>> acc(nth);
  for (auto& v : acc) v.fill(0.0);

  const bool frozen = frozen_;
  // Phase 1: pair list. Every worker is called (possibly with an empty range) and clears its own force buffer.
  pool_->run(nonb ? pi_.size() : 0, [&](int t, size_t b, size_t en) {
    std::vector<double>& ft = tf_[t];
    ft.assign(x.size(), 0.0);
    double evdw = 0, ecoul = 0, vir = 0, wv[6] = {0, 0, 0, 0, 0, 0};
    for (size_t p = b; p < en; ++p) {
      const uint32_t i = pi_[p], j = pj_[p];
      const double dx = x[3 * j] - x[3 * i] + shift_[3 * p], dy = x[3 * j + 1] - x[3 * i + 1] + shift_[3 * p + 1],
                   dz = x[3 * j + 2] - x[3 * i + 2] + shift_[3 * p + 2];
      const double r2 = dx * dx + dy * dy + dz * dz;
      if (frozen ? !inside_[p] : r2 >= rc2) continue;
      const double w = i == j ? 0.5 : 1.0;   // a self image is shared between the atom and its copy
      const size_t tp = size_t(ff_.type_index[i]) * nt + ff_.type_index[j];
      double ev = 0;
      double fr = lj(tp, r2, 1.0, ev);
      ev -= lj_shift(tp);
      evdw += w * ev;
      if (coul && q[i] != 0 && q[j] != 0) {
        const double r = std::sqrt(r2), qq = kCoulomb * q[i] * q[j];
        double ex2;
        if (pme) {
          const double er = erfc_exp(beta * r, ex2);
          ecoul += w * qq * er / r;
          fr += qq * (er / r2 + b2pi * ex2 / r) / r;
        } else {
          const double er = erfc_exp(a * r, ex2);
          ecoul += w * qq * (er / r - dsf_e0 + dsf_f0 * (r - rc));
          fr += qq * (er / r2 + a2pi * ex2 / r - dsf_f0) / r;
        }
      }
      fr *= w;
      const double fx = fr * dx, fy = fr * dy, fz = fr * dz;
      ft[3 * j] += fx; ft[3 * j + 1] += fy; ft[3 * j + 2] += fz;
      ft[3 * i] -= fx; ft[3 * i + 1] -= fy; ft[3 * i + 2] -= fz;
      vir += fx * dx + fy * dy + fz * dz;
      wv[0] += fx * dx; wv[1] += fy * dy; wv[2] += fz * dz;
      wv[3] += fx * dy; wv[4] += fx * dz; wv[5] += fy * dz;   // pair forces are central: f ∥ d, so f_a d_b = f_b d_a
    }
    acc[t][4] += evdw;
    acc[t][5] += ecoul;
    acc[t][6] += vir;
    for (int c = 0; c < 6; ++c) acc[t][7 + c] += wv[c];
  });

  // Phase 2: bonded terms, 1-4 pairs and the electrostatics of bonded partners, as one index space split over workers.
  const size_t nb = ff_.bonds.size(), na = ff_.angles.size(), nd = ff_.dihedrals.size(), ni = ff_.impropers.size();
  const size_t n14 = ff_.pairs14.size(), nx = coul ? excl_.size() : 0, nh = ff_.impropers_harmonic.size();
  const size_t o1 = nb, o2 = o1 + na, o3 = o2 + nd, o4 = o3 + ni, o5 = o4 + n14, o6 = o5 + nx, o7 = o6 + nh;
  const size_t c1 = o7 + ff_.bonds2.size(), c2 = c1 + ff_.angles2.size(), c3 = c2 + ff_.dihedrals2.size();
  const size_t c4 = c3 + ff_.impropers2.size();
  const size_t c5 = c4 + ff_.inversions.size(), c6 = c5 + ff_.bonds_x.size(), c7 = c6 + ff_.angles_x.size();
  const size_t total = c7 + ff_.urey_bradley.size();
  pool_->run(total, [&](int t, size_t b, size_t en) {
    std::vector<double>& ft = tf_[t];
    auto& A = acc[t];
    auto add = [&](uint32_t i, const Vec3& v) { ft[3 * i] += v[0]; ft[3 * i + 1] += v[1]; ft[3 * i + 2] += v[2]; };
    // virial of force fv at relative position d, scalar and tensor
    auto V = [&](const Vec3& d, const Vec3& fv) {
      A[6] += dot(d, fv);
      A[7] += d[0] * fv[0]; A[8] += d[1] * fv[1]; A[9] += d[2] * fv[2];
      A[10] += 0.5 * (d[0] * fv[1] + d[1] * fv[0]); A[11] += 0.5 * (d[0] * fv[2] + d[2] * fv[0]); A[12] += 0.5 * (d[1] * fv[2] + d[2] * fv[1]);
    };
    // Torsions (proper and improper share the form).
    // φ follows IUPAC (atan2 of |b2| b1·n and m·n); Field.ForcesMatchFiniteDifferences checks the gradients.
    auto torsion = [&](const TorsionTerm& tt, double& en) {
      const Vec3 b1 = mi(pos(tt.j) - pos(tt.i)), b2 = mi(pos(tt.k) - pos(tt.j)), b3 = mi(pos(tt.l) - pos(tt.k));
      const Vec3 m = cross(b1, b2), nn = cross(b2, b3);
      const double m2 = dot(m, m), n2 = dot(nn, nn), lb2 = norm(b2);
      if (m2 < 1e-12 || n2 < 1e-12 || lb2 < 1e-9) return;
      const double phi = std::atan2(lb2 * dot(b1, nn), dot(m, nn));
      en += tt.v * (1 + std::cos(tt.n * phi - tt.delta));
      const double dEdphi = -tt.v * tt.n * std::sin(tt.n * phi - tt.delta);
      const Vec3 g1 = m * (-lb2 / m2), g4 = nn * (lb2 / n2);
      const double p = dot(b1, b2) / (lb2 * lb2), qd = dot(b3, b2) / (lb2 * lb2);
      const Vec3 g2 = g4 * qd - g1 * (1 + p);
      const Vec3 g3 = g1 * p - g4 * (1 + qd);
      const Vec3 f1 = g1 * -dEdphi, f2 = g2 * -dEdphi, f3 = g3 * -dEdphi, f4 = g4 * -dEdphi;
      add(tt.i, f1); add(tt.j, f2); add(tt.k, f3); add(tt.l, f4);
      // positions relative to atom i
      const Vec3 r2 = b1, r3 = b1 + b2, r4 = b1 + b2 + b3;
      V(r2, f2); V(r3, f3); V(r4, f4);
    };
    for (size_t k = b; k < en; ++k) {
      // 1-4 pairs and the electrostatics of bonded partners are non-bonded; everything else is bonded
      if ((k >= o4 && k < o6) ? !nonb : !bonded) continue;
      if (k < o1) {
        const auto& bd = ff_.bonds[k];
        const Vec3 d = mi(pos(bd.j) - pos(bd.i));
        const double r = norm(d), dr = r - bd.r0;
        A[0] += bd.k * dr * dr;
        const Vec3 fj = r > 1e-12 ? d * (-2 * bd.k * dr / r) : Vec3{0, 0, 0};
        add(bd.j, fj);
        add(bd.i, fj * -1.0);
        V(d, fj);
      } else if (k < o2) {
        const auto& an = ff_.angles[k - o1];
        const Vec3 u = mi(pos(an.i) - pos(an.j)), v = mi(pos(an.k) - pos(an.j));
        const double lu = norm(u), lv = norm(v);
        const double c = std::clamp(dot(u, v) / (lu * lv), -1.0, 1.0);
        const double th = std::acos(c), sn = std::max(std::sqrt(1 - c * c), 1e-8);
        const double dth = th - an.theta0;
        A[1] += an.kt * dth * dth;
        const double g = 2 * an.kt * dth / sn;   // −dE/dθ · (−1/sinθ) = dE/dθ / sinθ
        const Vec3 fi = (v * (1 / (lu * lv)) - u * (c / (lu * lu))) * g;
        const Vec3 fk = (u * (1 / (lu * lv)) - v * (c / (lv * lv))) * g;
        add(an.i, fi);
        add(an.k, fk);
        add(an.j, (fi + fk) * -1.0);
        V(u, fi); V(v, fk);
      } else if (k < o3) {
        torsion(ff_.dihedrals[k - o2], A[2]);
      } else if (k < o4) {
        torsion(ff_.impropers[k - o3], A[3]);
      } else if (k < o5) {
        // 1-4 pairs: LJ scaled (with the same cut-off shift), at the nearest image.
        const auto& pr = ff_.pairs14[k - o4];
        const uint32_t i = pr[0], j = pr[1];
        const Vec3 d = mi(pos(j) - pos(i));
        const size_t tp = size_t(ff_.type_index[i]) * nt + ff_.type_index[j];
        double ev = 0, fr = 0;
        if (!eps14_.empty()) {
          // separate 1-4 Lennard-Jones parameters (CHARMM, GROMOS)
          const double r2 = dot(d, d), qq = s614_[tp] / (r2 * r2 * r2), e14 = eps14_[tp];
          ev = ff_.lj14 * 4 * e14 * (qq * qq - qq);
          fr = ff_.lj14 * 24 * e14 * (2 * qq * qq - qq) / r2;
          if (!opt_.tail) {
            const double qc = s614_[tp] / (rc2 * rc2 * rc2);
            ev -= ff_.lj14 * 4 * e14 * (qc * qc - qc);
          }
        } else {
          fr = lj(tp, dot(d, d), ff_.lj14, ev);
          ev -= ff_.lj14 * lj_shift(tp);
        }
        A[4] += ev;
        const Vec3 fj = d * fr;
        add(j, fj);
        add(i, fj * -1.0);
        V(d, fj);
      } else if (k >= o7) {
        // Class II terms: energy as a function of internal coordinates, forces by the chain rule.
        Class2Out o;
        if (k < c1) {
          const auto& t = ff_.bonds2[k - o7];
          o = class2_bond(t, mi(pos(t.j) - pos(t.i)));
          add(t.i, o.f[0]); add(t.j, o.f[1]);
          A[0] += o.e;
        } else if (k < c2) {
          const auto& t = ff_.angles2[k - c1];
          o = class2_angle(t, mi(pos(t.i) - pos(t.j)), mi(pos(t.k) - pos(t.j)));
          add(t.i, o.f[0]); add(t.j, o.f[1]); add(t.k, o.f[2]);
          A[1] += o.e;
        } else if (k < c3) {
          const auto& t = ff_.dihedrals2[k - c2];
          o = class2_dihedral(t, mi(pos(t.j) - pos(t.i)), mi(pos(t.k) - pos(t.j)), mi(pos(t.l) - pos(t.k)));
          add(t.i, o.f[0]); add(t.j, o.f[1]); add(t.k, o.f[2]); add(t.l, o.f[3]);
          A[2] += o.e;
        } else if (k >= c7) {
          const auto& t = ff_.urey_bradley[k - c7];
          const Vec3 d = mi(pos(t.k) - pos(t.i));
          const double r = norm(d), dr = r - t.r0;
          o.e = t.kub * dr * dr;
          o.f[1] = d * (-2 * t.kub * dr / r);
          o.f[0] = o.f[1] * -1.0;
          o.vir(d, o.f[1]);
          add(t.i, o.f[0]); add(t.k, o.f[1]);
          A[1] += o.e;
        } else if (k >= c6) {
          const auto& t = ff_.angles_x[k - c6];
          const Vec3 u = mi(pos(t.i) - pos(t.j)), v = mi(pos(t.k) - pos(t.j));
          Vec3 gu, gv;
          const double th = angle_grad(u, v, gu, gv), c = std::cos(th), dc = c - std::cos(t.b);
          double dEdth;
          if (t.form == 2) {   // K (1 + cos θ): linear centres (LAMMPS angle cosine, DREIDING)
            o.e = t.a * (1 + c);
            dEdth = -t.a * std::sin(th);
          } else if (t.form == 3) {   // UFF K [C0 + C1 cos θ + C2 cos 2θ] with the C from θ0 (LAMMPS angle fourier)
            const double s0 = std::sin(t.b), c0 = std::cos(t.b);
            const double C2 = 1 / (4 * std::max(s0 * s0, 1e-8)), C1 = -4 * C2 * c0, C0 = C2 * (2 * c0 * c0 + 1);
            o.e = t.a * (C0 + C1 * c + C2 * std::cos(2 * th));
            dEdth = -t.a * (C1 * std::sin(th) + 2 * C2 * std::sin(2 * th));
          } else if (t.form > 10) {   // UFF K (1 − cos nθ) / n², n = form − 10; n = 1 is K (1 + cos θ) (LAMMPS cosine/periodic)
            const int nn = t.form - 10;
            if (nn == 1) {
              o.e = t.a * (1 + c);
              dEdth = -t.a * std::sin(th);
            } else {
              o.e = t.a * (1 - std::cos(nn * th)) / (nn * nn);
              dEdth = t.a * std::sin(nn * th) / nn;
            }
            if (c > 0.8660) {   // UFF's wall below 30°, where these periodic forms have a spurious minimum (as RDKit)
              const double w = std::exp(-20 * (th - t.b + 0.25));
              o.e += w;
              dEdth += -20 * w;
            }
          } else {
            o.e = t.a * dc * dc;
            dEdth = -2 * t.a * dc * std::sin(th);
          }
          o.f[0] = gu * -dEdth; o.f[2] = gv * -dEdth; o.f[1] = (o.f[0] + o.f[2]) * -1.0;
          o.vir(u, o.f[0]); o.vir(v, o.f[2]);
          add(t.i, o.f[0]); add(t.j, o.f[1]); add(t.k, o.f[2]);
          A[1] += o.e;
        } else if (k >= c5) {
          const auto& t = ff_.bonds_x[k - c5];
          const Vec3 d = mi(pos(t.j) - pos(t.i));
          const double r = norm(d);
          double dEdr = 0;
          if (t.form == 1) {
            const double e1 = std::exp(-t.b * (r - t.c));
            o.e = t.a * (1 - e1) * (1 - e1);
            dEdr = 2 * t.a * t.b * e1 * (1 - e1);
          } else {
            const double q = r * r - t.b * t.b;   // LAMMPS bond gromos: K (r² − r0²)², K = ¼ K_GROMOS
            o.e = t.a * q * q;
            dEdr = 4 * t.a * q * r;
          }
          o.f[1] = r > 1e-12 ? d * (-dEdr / r) : Vec3{0, 0, 0};
          o.f[0] = o.f[1] * -1.0;
          o.vir(d, o.f[1]);
          add(t.i, o.f[0]); add(t.j, o.f[1]);
          A[0] += o.e;
        } else if (k >= c4) {
          const auto& t = ff_.inversions[k - c4];
          o = inversion(t, mi(pos(t.a) - pos(t.c)), mi(pos(t.b) - pos(t.c)), mi(pos(t.d) - pos(t.c)));
          add(t.a, o.f[0]); add(t.c, o.f[1]); add(t.b, o.f[2]); add(t.d, o.f[3]);
          A[3] += o.e;
        } else {
          const auto& t = ff_.impropers2[k - c3];
          o = class2_improper(t, mi(pos(t.i) - pos(t.j)), mi(pos(t.k) - pos(t.j)), mi(pos(t.l) - pos(t.j)));
          add(t.i, o.f[0]); add(t.j, o.f[1]); add(t.k, o.f[2]); add(t.l, o.f[3]);
          A[3] += o.e;
        }
        A[6] += o.virial;
        for (int c = 0; c < 6; ++c) A[7 + c] += o.w[c];
      } else if (k >= o6) {
        // Harmonic improper k2 (χ − χ0)² on the i-j-k-l dihedral angle (LAMMPS improper_style harmonic).
        const auto& h = ff_.impropers_harmonic[k - o6];
        const Vec3 b1 = mi(pos(h.j) - pos(h.i)), b2 = mi(pos(h.k) - pos(h.j)), b3 = mi(pos(h.l) - pos(h.k));
        const Vec3 m = cross(b1, b2), nn = cross(b2, b3);
        const double m2 = dot(m, m), n2 = dot(nn, nn), lb2 = norm(b2);
        if (m2 < 1e-12 || n2 < 1e-12 || lb2 < 1e-9) continue;
        const double phi = std::atan2(lb2 * dot(b1, nn), dot(m, nn));
        double dchi = phi - h.chi0;
        dchi -= 2 * kPi * std::round(dchi / (2 * kPi));
        A[3] += h.k2 * dchi * dchi;
        const double dEdphi = 2 * h.k2 * dchi;
        const Vec3 g1 = m * (-lb2 / m2), g4 = nn * (lb2 / n2);
        const double p = dot(b1, b2) / (lb2 * lb2), qd = dot(b3, b2) / (lb2 * lb2);
        const Vec3 g2 = g4 * qd - g1 * (1 + p), g3 = g1 * p - g4 * (1 + qd);
        const Vec3 f1 = g1 * -dEdphi, f2 = g2 * -dEdphi, f3 = g3 * -dEdphi, f4 = g4 * -dEdphi;
        add(h.i, f1); add(h.j, f2); add(h.k, f3); add(h.l, f4);
        V(b1, f2); V(b1 + b2, f3); V(b1 + b2 + b3, f4);
      } else {
        // Electrostatics of bonded partners, as LAMMPS coul/dsf does it: the damped sum includes them and the bare
        // Coulomb term is removed in proportion (1 − factor): factor 0 for 1-2 and 1-3 pairs, 1/1.2 for 1-4 pairs.
        const auto& ex = excl_[k - o5];
        const uint32_t i = ex.i, j = ex.j;
        if (q[i] == 0 || q[j] == 0) continue;
        const Vec3 d = mi(pos(j) - pos(i));
        const double r2 = dot(d, d);
        if (r2 >= rc2) continue;
        const double r = std::sqrt(r2), qq = kCoulomb * q[i] * q[j];
        double ex2, fr;
        if (r < 0.5) {
          // a core and its shell (or any pair this close): erfc(κr)/r − 1/r cancels catastrophically as r → 0, so the
          // same energy is written −erf(κr)/r (exact, with its r → 0 limits) plus what the factor keeps of 1/r
          const double kap = pme ? beta : a, sp = 2 / std::sqrt(kPi);
          const double x = kap * r;
          const double g = r > 1e-8 ? std::erf(x) / r : sp * kap;                                   // erf(κr)/r
          const double gr = r > 1e-4 ? (sp * kap * std::exp(-x * x) * r - std::erf(x)) / (r2 * r)   // g'(r) / r
                                     : -2 * sp * kap * kap * kap / 3;
          const double keep = ex.factor;   // the part of the bare 1/r the exclusion keeps (0 for 1-2 pairs)
          const double shift = pme ? 0.0 : -dsf_e0 + dsf_f0 * (r - rc);
          const double inv = r > 1e-8 ? 1 / r : 0.0;
          A[5] += qq * (-g + keep * inv + shift);
          fr = qq * (gr + keep * inv * inv * inv - (pme ? 0.0 : dsf_f0) * inv);
          const Vec3 fj = d * fr;
          add(j, fj);
          add(i, fj * -1.0);
          V(d, fj);
          continue;
        }
        if (pme) {
          // Ewald: the reciprocal sum includes the pair; add its real-space part and remove (1 − factor) of 1/r
          const double er = erfc_exp(beta * r, ex2);
          A[5] += qq * er / r - (1 - ex.factor) * qq / r;
          fr = qq * ((er / r2 + b2pi * ex2 / r) - (1 - ex.factor) / r2) / r;
        } else {
          const double er = erfc_exp(a * r, ex2);
          A[5] += qq * (er / r - dsf_e0 + dsf_f0 * (r - rc)) - (1 - ex.factor) * qq / r;
          fr = qq * ((er / r2 + a2pi * ex2 / r - dsf_f0) - (1 - ex.factor) / r2) / r;
        }
        const Vec3 fj = d * fr;
        add(j, fj);
        add(i, fj * -1.0);
        V(d, fj);
      }
    }
  });

  // Reduce in worker order (deterministic for a given thread count); the reduction itself is split by atoms.
  pool_->run(x.size(), [&](int, size_t b, size_t en) {
    for (int t = 0; t < nth; ++t) {
      const std::vector<double>& ft = tf_[t];
      for (size_t k = b; k < en; ++k) f[k] += ft[k];
    }
  });
  for (int t = 0; t < nth; ++t) {
    e.bond += acc[t][0];
    e.angle += acc[t][1];
    e.dihedral += acc[t][2];
    e.improper += acc[t][3];
    e.vdw += acc[t][4];
    e.coulomb += acc[t][5];
    e.virial += acc[t][6];
    for (int c = 0; c < 6; ++c) e.w[c] += acc[t][7 + c];
  }
  if (!nonb) return e;   // r-RESPA inner step: bonded terms only
  if (ff_.sw.on) {   // Stillinger–Weber: pairs and triplets within aσ, from the pair list (its shifts carry the images)
    const auto& S = ff_.sw;
    const double rcs = S.a * S.sigma, rcs2 = rcs * rcs;
    struct Nb { uint32_t j; Vec3 d; double r; };
    std::vector<std::vector<Nb>> nbs(n);
    double esw = 0;
    auto push = [&](uint32_t i, const Vec3& fv) { f[3 * i] += fv[0]; f[3 * i + 1] += fv[1]; f[3 * i + 2] += fv[2]; };
    auto vir = [&](const Vec3& d, const Vec3& fv) {
      e.virial += dot(d, fv);
      e.w[0] += d[0] * fv[0]; e.w[1] += d[1] * fv[1]; e.w[2] += d[2] * fv[2];
      e.w[3] += 0.5 * (d[0] * fv[1] + d[1] * fv[0]); e.w[4] += 0.5 * (d[0] * fv[2] + d[2] * fv[0]); e.w[5] += 0.5 * (d[1] * fv[2] + d[2] * fv[1]);
    };
    for (size_t k = 0; k < pi_.size(); ++k) {
      const uint32_t i = pi_[k], j = pj_[k];
      if (!S.atom[i] || !S.atom[j]) continue;
      const Vec3 d = pos(j) - pos(i) + Vec3{shift_[3 * k], shift_[3 * k + 1], shift_[3 * k + 2]};
      const double r2 = dot(d, d);
      if (r2 >= rcs2) continue;
      const double r = std::sqrt(r2);
      nbs[i].push_back({j, d, r});
      nbs[j].push_back({i, d * -1.0, r});
      // two-body
      const double sr = S.sigma / r, srp = std::pow(sr, S.p), srq = std::pow(sr, S.q), ex = std::exp(S.sigma / (r - rcs));
      const double u = S.B * srp - srq;
      esw += S.A * S.eps * u * ex;
      const double dudr = (-S.p * S.B * srp + S.q * srq) / r;
      const double dEdr = S.A * S.eps * ex * (dudr - u * S.sigma / ((r - rcs) * (r - rcs)));
      const Vec3 fj = d * (-dEdr / r);
      push(j, fj);
      push(i, fj * -1.0);
      vir(d, fj);
    }
    // three-body: every pair of neighbours j, k of a vertex i
    for (uint32_t i = 0; i < n; ++i) {
      const auto& L = nbs[i];
      for (size_t a1 = 0; a1 < L.size(); ++a1)
        for (size_t a2 = a1 + 1; a2 < L.size(); ++a2) {
          const Nb& J = L[a1];
          const Nb& K = L[a2];
          const double c = dot(J.d, K.d) / (J.r * K.r), dc = c - S.cos0;
          const double ej = std::exp(S.gamma * S.sigma / (J.r - rcs)), ek = std::exp(S.gamma * S.sigma / (K.r - rcs));
          const double E3 = S.lambda * S.eps * dc * dc * ej * ek;
          esw += E3;
          const double dEdc = 2 * S.lambda * S.eps * dc * ej * ek;
          const double dEdrj = -E3 * S.gamma * S.sigma / ((J.r - rcs) * (J.r - rcs));
          const double dEdrk = -E3 * S.gamma * S.sigma / ((K.r - rcs) * (K.r - rcs));
          // ∂cos/∂d_j = d_k/(r_j r_k) − cos d_j / r_j²
          const Vec3 gj = K.d * (1 / (J.r * K.r)) - J.d * (c / (J.r * J.r));
          const Vec3 gk = J.d * (1 / (J.r * K.r)) - K.d * (c / (K.r * K.r));
          const Vec3 fj = (J.d * (dEdrj / J.r) + gj * dEdc) * -1.0;
          const Vec3 fk = (K.d * (dEdrk / K.r) + gk * dEdc) * -1.0;
          push(J.j, fj);
          push(K.j, fk);
          push(i, (fj + fk) * -1.0);
          vir(J.d, fj);
          vir(K.d, fk);
        }
    }
    e.vdw += esw;
  }
  if (pme) {
    // reciprocal part (Fortran), self energy −β/√π Σq², and the neutralising background −π Q²/(2Vβ²) for a net charge
    const PmeGrid grid = pme_grid(cell, beta, opt_.pme_spacing, opt_.pme_order);
    std::vector<double> fk(x.size(), 0.0);
    double vk[6];
    e.coulomb += pme_reciprocal(x, ff_.charge, cell, grid, fk, vk, *pool_);
    for (size_t k = 0; k < x.size(); ++k) f[k] += fk[k];
    e.virial += vk[0] + vk[1] + vk[2];
    for (int c = 0; c < 6; ++c) e.w[c] += vk[c];
    double q2 = 0, qt = 0;
    for (double c : ff_.charge) q2 += c * c, qt += c;
    e.coulomb -= kCoulomb * beta / std::sqrt(kPi) * q2;
    if (std::fabs(qt) > 1e-10) {
      const double eb = -kCoulomb * kPi * qt * qt / (2 * cell.volume() * beta * beta);
      e.coulomb += eb;
      e.virial += 3 * eb;   // E_bg ∝ 1/V ∝ λ⁻³, so −dE/dλ = 3 E_bg
      for (int c = 0; c < 3; ++c) e.w[c] += eb;
    }
  } else if (coul) {
    double q2 = 0;
    for (double c : ff_.charge) q2 += c * c;
    // Self energy: half the r → 0 limit of the damped shifted pair potential minus the bare 1/r, i.e. with the full
    // shift (potential and force-shift parts), as LAMMPS coul/dsf (e_shift = erfc(αrc)/rc + f0 rc). A constant: no force.
    e.coulomb -= (erfc_rc / (2 * rc) + dsf_f0 * rc / 2 + a / std::sqrt(kPi)) * q2 * kCoulomb;
  }

  // Tail corrections for a homogeneous fluid beyond the cut-off (Allen and Tildesley, 2nd ed., §2.8), as LAMMPS
  // pair_modify tail yes. Per ordered type pair: E = 8π N_a N_b ε [σ¹²/(9rc⁹) − σ⁶/(3rc³)] / (2V) summed twice, and
  // P = 16π N_a N_b ε [2σ¹²/(3rc⁹) − σ⁶/rc³] / (6V²), kept as a virial 3PV.
  if (opt_.tail && cell.valid()) {
    const double vol = cell.volume(), rc3 = rc2 * rc, rc9 = rc3 * rc3 * rc3;
    double et = 0, wt = 0;
    for (size_t a0 = 0; a0 < nt; ++a0)
      for (size_t b0 = 0; b0 < nt; ++b0) {
        const size_t tp = a0 * nt + b0;
        const double eps = eps_[tp], s6 = s6_[tp], s12 = s6 * s6, nn = type_count_[a0] * type_count_[b0];
        if (form_[tp] == 1) {
          // Buckingham, as LAMMPS pair buck: E = 2π N N [Aρ e^(−rc/ρ)(rc² + 2ρrc + 2ρ²) − C/(3rc³)] / V
          const double A = pa_[tp], rho = pb_[tp], Cc = pc_[tp], ex = std::exp(-rc / rho);
          et += nn * 2 * kPi * (A * rho * ex * (rc2 + 2 * rho * rc + 2 * rho * rho) - Cc / (3 * rc3));
          wt += nn * 2 * kPi * (A * ex * (rc3 + 3 * rho * rc2 + 6 * rho * rho * rc + 6 * rho * rho * rho) - 2 * Cc / rc3);
          continue;
        }
        if (form_[tp] == 2) continue;   // Morse: no tail term (LAMMPS pair morse)
        if (lj96_) {
          // lj/class2: E = 2π N_a N_b ε σ⁶ (σ³ − 3rc³) / (3rc⁶ V), P = 2π N_a N_b ε σ⁶ (σ³ − 2rc³) / (rc⁶ V²) per ordered pair
          const double s3 = std::sqrt(s6), rc6 = rc3 * rc3;
          et += nn * 2 * kPi * eps * s6 * (s3 - 3 * rc3) / (3 * rc6);
          wt += nn * 3 * 2 * kPi * eps * s6 * (s3 - 2 * rc3) / rc6;
          continue;
        }
        et += nn * 2 * kPi * 4 * eps * (s12 / (9 * rc9) - s6 / (3 * rc3));
        wt += nn * 2 * kPi * 4 * eps * (4 * s12 / (3 * rc9) - 2 * s6 / rc3);   // Σ r·f equivalent, so P_tail = wt / (3V)
      }
    e.vdw += et / vol;
    e.virial += wt / vol;
    for (int c = 0; c < 3; ++c) e.w[c] += wt / vol / 3;   // isotropic
  }
  return e;
}

void Evaluator::freeze_pairs(const std::vector<double>& x, const Cell& cell) {
  frozen_ = false;
  build(x, cell);
  const double rc2 = opt_.cutoff * opt_.cutoff;
  inside_.assign(pi_.size(), 0);
  for (size_t p = 0; p < pi_.size(); ++p) {
    const uint32_t i = pi_[p], j = pj_[p];
    double r2 = 0;
    for (int k = 0; k < 3; ++k) {
      const double d = x[3 * j + k] - x[3 * i + k] + shift_[3 * p + k];
      r2 += d * d;
    }
    inside_[p] = r2 < rc2;
  }
  frozen_ = true;
}

PairType mixed_pair(const ForceField& ff, int a, int b) {
  if (auto it = ff.pair_override.find({std::min(a, b), std::max(a, b)}); it != ff.pair_override.end()) return it->second;
  const double ea = ff.lj[a].eps, eb = ff.lj[b].eps, sa = ff.lj[a].sigma, sb = ff.lj[b].sigma;
  double e = std::sqrt(ea * eb), sg = 0.5 * (sa + sb);
  if (ff.mixing == "geometric") sg = std::sqrt(sa * sb);
  else if (ff.mixing == "sixthpower") {
    const double s6a = std::pow(sa, 6), s6b = std::pow(sb, 6);
    sg = std::pow(0.5 * (s6a + s6b), 1.0 / 6);
    e = s6a + s6b > 0 ? 2 * std::sqrt(ea * eb) * sa * sa * sa * sb * sb * sb / (s6a + s6b) : 0;
  }
  return {e, sg};
}

double pressure_atm(double virial, double volume) {
  constexpr double kAtm = 68568.415;   // kcal/(mol·Å³) → atm
  return volume > 0 ? virial / (3 * volume) * kAtm : 0.0;
}

}  // namespace caps
