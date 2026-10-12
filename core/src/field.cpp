// GAFF typing and parameters for hydrocarbons, and the energy / force evaluator used by Relax.
//
// Parameters: GAFF 1.81 (Wang, Wolf, Caldwell, Kollman and Case, J. Comput. Chem. 25, 1157 (2004)), gaff.dat values
// for the c3 / ca / hc / ha subset. AMBER forms and 1-4 scaling (LJ 1/2, Coulomb 1/1.2). Electrostatics use the damped
// shifted force sum (Fennell and Gezelter, J. Chem. Phys. 124, 234104 (2006)), a real-space method with no Ewald part.
#include "caps/field.hpp"
#include "caps/kspace.hpp"
#include "kspace_pool.hpp"
#include "nbnxm.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <functional>
#include <map>
#include <set>
#include <type_traits>

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
  if (ff.cutoff > 0) opt_.cutoff = ff.cutoff;
  if (ff.lj_shift) opt_.tail = false;   // Martini 3: shifted at the cut-off, no tail correction
  if (ff.lj_fsw) {   // CHARMM: switched to zero at the cut-off, no tail correction
    opt_.tail = false;
    if (ff.lj_inner <= 0 || ff.lj_inner >= opt_.cutoff) throw FieldError("the CHARMM force switch needs an inner radius below the cut-off");
  }
  // a relative permittivity εr (MARTINI 15, SDK 80) divides every Coulomb term: charges scaled by 1/√εr throughout
  if (ff.dielectric <= 0) throw FieldError("the relative permittivity must be positive");
  qeff_ = ff.charge;
  if (ff.dielectric != 1)
    for (auto& c : qeff_) c /= std::sqrt(ff.dielectric);
  std::map<std::pair<uint32_t, uint32_t>, double> p14;   // each 1-4 pair's Coulomb scale (its part's, when merged by group)
  for (size_t k = 0; k < ff.pairs14.size(); ++k) {
    const auto& p = ff.pairs14[k];
    const double c = k < ff.pairs14_coul.size() ? ff.pairs14_coul[k] : ff.coul14;
    p14[{std::min(p[0], p[1]), std::max(p[0], p[1])}] = c;
  }
  for (uint32_t i = 0; i < ff.excluded.size(); ++i)
    for (uint32_t j : ff.excluded[i])
      if (j > i) {
        const auto it = p14.find({i, j});
        excl_.push_back({i, j, it != p14.end() ? it->second : 0.0});
      }
  const size_t nt = ff.lj.size();
  eps_.resize(nt * nt);
  s6_.resize(nt * nt);
  for (size_t a = 0; a < nt; ++a)
    for (size_t b = 0; b < nt; ++b) {
      const PairType pt = mixed_pair(ff, int(a), int(b));
      eps_[a * nt + b] = pt.eps;
      s6_[a * nt + b] = std::pow(pt.sigma, 6);
    }
  skip_type_.assign(nt * nt, 0);
  for (const auto& [a, b] : ff.excluded_type_pairs) skip_type_[size_t(a) * nt + b] = skip_type_[size_t(b) * nt + a] = 1;
  form_.assign(nt * nt, 0);
  pa_.assign(nt * nt, 0.0);
  pb_.assign(nt * nt, 0.0);
  pc_.assign(nt * nt, 0.0);
  for (const auto& [ab, pf] : ff.pair_func) {
    if (pf.form == kPairTable && (pf.a < 0 || size_t(pf.a) >= ff.tables.size())) throw FieldError("a tabulated pair refers to a table that is not there");
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

void Evaluator::set_options(const EnergyOptions& o0) {
  EnergyOptions o = o0;
  if (ff_.cutoff > 0) o.cutoff = ff_.cutoff;   // the model's own cut-off
  if (ff_.lj_shift || ff_.lj_fsw) o.tail = false;
  if (o.threads != opt_.threads) pool_ = std::make_unique<ThreadPool>(o.threads > 0 ? o.threads : default_threads()), cl_built_ = false;
  const bool relist = o.cutoff != opt_.cutoff || o.skin != opt_.skin;
  opt_ = o;
  if (relist) built_ = false, cl_built_ = false;
  cap_radii();
}

// For each type pair, the radius inside which |F_LJ| exceeds the cap; the potential is continued linearly inside it.
void Evaluator::cap_radii() {
  const size_t nt2 = eps_.size();
  // lj/gromacs switch coefficients (LAMMPS pair lj/gromacs) and the SDK minima, for the cut-off in use
  gsw_.assign(5 * nt2, 0.0);
  rmin2_.assign(nt2, 0.0);
  emin_.assign(nt2, 0.0);
  const double rc = opt_.cutoff, r1 = ff_.lj_inner;
  for (size_t t = 0; t < nt2; ++t) {
    if (form_[t] == kPairGromacs) {
      if (r1 <= 0 || r1 >= rc) throw FieldError("lj/gromacs needs an inner radius below the cut-off");
      const double e = pa_[t], sg = pb_[t], s6 = std::pow(sg, 6), s12 = s6 * s6;
      const double l1 = 48 * e * s12, l2 = 24 * e * s6, l3 = 4 * e * s12, l4 = 4 * e * s6;
      const double r6i = 1 / std::pow(rc, 6), r8i = 1 / std::pow(rc, 8), tt = rc - r1, t2i = 1 / (tt * tt), t3i = t2i / tt, t3 = tt * tt * tt;
      const double a6 = (7 * r1 - 10 * rc) * r8i * t2i, b6 = (9 * rc - 7 * r1) * r8i * t3i;
      const double a12 = (13 * r1 - 16 * rc) * r6i * r8i * t2i, b12 = (15 * rc - 13 * r1) * r6i * r8i * t3i;
      const double c6 = r6i - t3 * (6 * a6 / 3 + 6 * b6 * tt / 4), c12 = r6i * r6i - t3 * (12 * a12 / 3 + 12 * b12 * tt / 4);
      double* g = &gsw_[5 * t];
      g[0] = l1 * a12 - l2 * a6;
      g[1] = l1 * b12 - l2 * b6;
      g[2] = -l3 * 12 * a12 / 3 + l4 * 6 * a6 / 3;
      g[3] = -l3 * 12 * b12 / 4 + l4 * 6 * b6 / 4;
      g[4] = -l3 * c12 + l4 * c6;
    } else if (form_[t] >= kPairSdk96 && form_[t] <= kPairSdk125) {
      static const double pw1[] = {9, 12, 12, 12}, pw2[] = {6, 4, 6, 5}, pre[] = {6.75, 2.59807621135332, 4.0, 3.20377984125109};
      const int k = form_[t] - kPairSdk96;
      const double rmin = pb_[t] * std::exp(std::log(pw1[k] / pw2[k]) / (pw1[k] - pw2[k])), q = pb_[t] / rmin;
      rmin2_[t] = rmin * rmin;
      emin_[t] = pre[k] * pa_[t] * (std::pow(q, pw1[k]) - std::pow(q, pw2[k]));
    }
  }
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

namespace {

// Whether a pair list built at (x0, cell0), its images now for cells, still serves (x, cell): every atom has moved less
// than half the skin (relative to the affinely mapped positions when the cell was deformed, shrinking eating into the
// skin). 0 it serves, 1 it serves once its images are rescaled to cell, 2 it must be rebuilt.
int list_state(const std::vector<double>& x, const Cell& cell, const std::vector<double>& x0, const Cell& cell0, const Cell& cells,
               double skin, double cutoff) {
  auto same = [](const Cell& p, const Cell& q) { return p.a == q.a && p.b == q.b && p.c == q.c && p.origin == q.origin; };
  const size_t n = x.size() / 3;
  if (x0.size() != x.size()) return 2;
  if (same(cell, cell0)) {
    const double lim = 0.25 * skin * skin;
    for (size_t i = 0; i < n; ++i) {
      const double dx = x[3 * i] - x0[3 * i], dy = x[3 * i + 1] - x0[3 * i + 1], dz = x[3 * i + 2] - x0[3 * i + 2];
      if (dx * dx + dy * dy + dz * dz > lim) return 2;
    }
    return 0;
  }
  if (!cell.valid() || !cell0.valid()) return 2;
  // Shrinking eats into the skin by (1 − s)(rc + skin), s the smallest ratio of perpendicular widths; the rest of the
  // skin covers moves relative to the affinely mapped positions.
  auto widths = [](const Cell& c) {
    const double v = c.volume();
    return Vec3{v / norm(cross(c.b, c.c)), v / norm(cross(c.c, c.a)), v / norm(cross(c.a, c.b))};
  };
  const Vec3 w0 = widths(cell0), w1 = widths(cell);
  const double sr = std::min({w1[0] / w0[0], w1[1] / w0[1], w1[2] / w0[2], 1.0});
  const double allow = skin - (1 - sr) * (cutoff + skin);
  if (allow < 0.1 * skin) return 2;
  const double lim = 0.25 * allow * allow;
  for (size_t i = 0; i < n; ++i) {
    const Vec3 m = cell.to_cartesian(cell0.to_fractional({x0[3 * i], x0[3 * i + 1], x0[3 * i + 2]}));
    const double dx = x[3 * i] - m[0], dy = x[3 * i + 1] - m[1], dz = x[3 * i + 2] - m[2];
    if (dx * dx + dy * dy + dz * dz > lim) return 2;
  }
  return same(cell, cells) ? 0 : 1;
}

}  // namespace

// The cluster kernel takes plain pairs only: 12-6 or 9-6 Lennard-Jones, PME, DSF or no charges, nothing that needs the
// atom pair list (Stillinger–Weber, a frozen list) or a branch per pair.
bool Evaluator::cluster_mode() const {
  if (!opt_.mixed || frozen_ || opt_.force_cap > 0 || ff_.lj_fsw || ff_.hbond.on() || ff_.sw.on) return false;
  if (opt_.alchemical && opt_.alchemical->size() == ff_.type_index.size()) return false;
  if (opt_.coulomb && (ff_.coul_rf || ff_.coul_gromacs)) return false;
  return std::all_of(form_.begin(), form_.end(), [](uint8_t v) { return v == 0; });
}

// Per type pair c12, c6 (12-6: 4εσ¹², 4εσ⁶; 9-6: 2εσ⁹, 3εσ⁶), and per-type factors when every pair's coefficients are
// their products (geometric mixing, as OPLS-AA): the kernel then needs no lookup per pair.
void Evaluator::cluster_coefficients() {
  const size_t nt = ff_.lj.size();
  cl_c12_.assign(nt * nt, 0.0f);
  cl_c6_.assign(nt * nt, 0.0f);
  std::vector<double> c12(nt * nt), c6(nt * nt);
  for (size_t tp = 0; tp < nt * nt; ++tp) {
    const double eps = eps_[tp] > 0 ? eps_[tp] : 0.0;
    c12[tp] = lj96_ ? 2 * eps * s6_[tp] * std::sqrt(s6_[tp]) : 4 * eps * s6_[tp] * s6_[tp];
    c6[tp] = lj96_ ? 3 * eps * s6_[tp] : 4 * eps * s6_[tp];
    cl_c12_[tp] = float(c12[tp]);
    cl_c6_[tp] = float(c6[tp]);
  }
  cl_f12_.assign(nt, 0.0);
  cl_f6_.assign(nt, 0.0);
  bool comb = true;
  for (size_t t = 0; t < nt; ++t) cl_f12_[t] = std::sqrt(c12[t * nt + t]), cl_f6_[t] = std::sqrt(c6[t * nt + t]);
  for (size_t a = 0; a < nt && comb; ++a)
    for (size_t b = 0; b < nt && comb; ++b) {
      const double p12 = cl_f12_[a] * cl_f12_[b], p6 = cl_f6_[a] * cl_f6_[b];
      comb = std::fabs(p12 - c12[a * nt + b]) <= 1e-7 * std::max(std::fabs(c12[a * nt + b]), 1e-300) &&
             std::fabs(p6 - c6[a * nt + b]) <= 1e-7 * std::max(std::fabs(c6[a * nt + b]), 1e-300);
    }
  if (!comb) cl_f12_.clear(), cl_f6_.clear();
  cl_coef_ = true;
}

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
    // bins are at least half the list radius wide: a pair up to rcs apart can be two bins away, periodic or not
    reach[k] = static_cast<int>(std::ceil(rcs / (w[k] / nb[k])));
  }
  // Atoms binned, then put in bin order (a counting sort, atoms of a bin in index order): slot s holds atom perm_[s],
  // and every bin is a contiguous range of slots. fr and bin are kept by slot.
  const size_t nbins = static_cast<size_t>(nb[0]) * nb[1] * nb[2];
  auto index = [&](int a, int b, int c) { return (static_cast<size_t>(a) * nb[1] + b) * nb[2] + c; };
  std::vector<Vec3> fra(n);
  std::vector<uint32_t> bina(n);
  std::vector<size_t> bstart(nbins + 1, 0);
  for (size_t i = 0; i < n; ++i) {
    Vec3 f = cell.to_fractional({x[3 * i], x[3 * i + 1], x[3 * i + 2]});
    int bi[3];
    for (int k = 0; k < 3; ++k) {
      if (cell.periodic[k]) f[k] -= std::floor(f[k]);
      f[k] = std::clamp(f[k], 0.0, 1.0 - 1e-12);
      bi[k] = std::min(nb[k] - 1, static_cast<int>(f[k] * nb[k]));
    }
    fra[i] = f;
    bina[i] = static_cast<uint32_t>(index(bi[0], bi[1], bi[2]));
    ++bstart[bina[i] + 1];
  }
  for (size_t b = 0; b < nbins; ++b) bstart[b + 1] += bstart[b];
  perm_.resize(n);
  iperm_.resize(n);
  {
    std::vector<size_t> fill(bstart.begin(), bstart.end() - 1);
    for (size_t i = 0; i < n; ++i) {
      const size_t s = fill[bina[i]]++;
      perm_[s] = static_cast<uint32_t>(i);
      iperm_[i] = static_cast<uint32_t>(s);
    }
  }
  std::vector<Vec3> fr(n);
  std::vector<std::array<int, 3>> bin(n);
  for (size_t s = 0; s < n; ++s) {
    const uint32_t i = perm_[s];
    fr[s] = fra[i];
    const size_t b = bina[i];
    bin[s] = {static_cast<int>(b / (size_t(nb[1]) * nb[2])), static_cast<int>(b / nb[2] % nb[1]), static_cast<int>(b % nb[2])};
  }
  // wrapped Cartesian positions by slot, and each image code's lattice vector: a candidate pair's separation is then
  // three subtractions (the shift taken once per stencil bin)
  std::vector<Vec3> cw(n);
  for (size_t s = 0; s < n; ++s) cw[s] = cell.a * fr[s][0] + cell.b * fr[s][1] + cell.c * fr[s][2];
  tslot_.resize(n);
  qslot_.resize(n);
  for (size_t s = 0; s < n; ++s) tslot_[s] = static_cast<uint32_t>(ff_.type_index[perm_[s]]), qslot_[s] = qeff_[perm_[s]];
  const double rcs2 = rcs * rcs;
  const size_t ntypes = ff_.lj.size();
  const bool skipping = !ff_.excluded_type_pairs.empty();
  // a neighbour bin u in [−reach, nb − 1 + reach] is the image q = ⌊u / nb⌋ of bin u − q nb: the codes of the images
  int qlo[3], span[3];
  for (int k = 0; k < 3; ++k) {
    qlo[k] = cell.periodic[k] ? -((reach[k] + nb[k] - 1) / nb[k]) : 0;
    span[k] = (cell.periodic[k] ? (nb[k] - 1 + reach[k]) / nb[k] : 0) - qlo[k] + 1;
  }
  if (size_t(span[0]) * span[1] * span[2] > 65536) throw FieldError("the cell is too small for the cut-off: more than 65536 periodic images in range");
  qint_.clear();
  for (int a = 0; a < span[0]; ++a)
    for (int b = 0; b < span[1]; ++b)
      for (int c = 0; c < span[2]; ++c) qint_.insert(qint_.end(), {qlo[0] + a, qlo[1] + b, qlo[2] + c});
  const int nth = pool_->size();
  // Every slot in a bin scans the same stencil: its bins in order of image code (then bin), worked out once per bin, so
  // a slot's pairs come out grouped by image, and within an image in slot order.
  const int sx = 2 * reach[0] + 1, sy = 2 * reach[1] + 1, sz = 2 * reach[2] + 1;
  struct Cand { uint16_t code; uint32_t bin; };
  std::vector<std::vector<uint32_t>> vj(nth);
  std::vector<std::vector<uint16_t>> vq(nth), vsc(nth);
  std::vector<std::vector<uint32_t>> vse(nth);
  std::vector<uint32_t> count(n), nseg(n);
  pool_->run(n, [&](int t, size_t b, size_t en) {
    auto& PJ = vj[t];
    auto& PQ = vq[t];
    auto& SC = vsc[t];
    auto& SE = vse[t];
    PJ.clear(); PQ.clear(); SC.clear(); SE.clear();
    std::vector<Cand> stencil;
    stencil.reserve(size_t(sx) * sy * sz);
    std::array<int, 3> sbin{-1, -1, -1};
    for (size_t i = b; i < en; ++i) {
      if (bin[i] != sbin) {
        sbin = bin[i];
        stencil.clear();
        for (int dx = -reach[0]; dx <= reach[0]; ++dx)
          for (int dy = -reach[1]; dy <= reach[1]; ++dy)
            for (int dz = -reach[2]; dz <= reach[2]; ++dz) {
              int u[3] = {sbin[0] + dx, sbin[1] + dy, sbin[2] + dz}, q[3] = {0, 0, 0};
              bool skip = false;
              for (int k = 0; k < 3; ++k) {
                if (u[k] < 0 || u[k] >= nb[k]) {
                  if (!cell.periodic[k]) { skip = true; break; }
                  q[k] = static_cast<int>(std::floor(static_cast<double>(u[k]) / nb[k]));
                  u[k] -= q[k] * nb[k];
                }
              }
              if (skip) continue;
              const uint16_t code = static_cast<uint16_t>(((q[0] - qlo[0]) * span[1] + (q[1] - qlo[1])) * span[2] + (q[2] - qlo[2]));
              stencil.push_back({code, static_cast<uint32_t>(index(u[0], u[1], u[2]))});
            }
        std::sort(stencil.begin(), stencil.end(), [](const Cand& p, const Cand& r) { return p.code != r.code ? p.code < r.code : p.bin < r.bin; });
      }
      const uint32_t ai = perm_[i];
      const auto& ex = ff_.excluded[ai];
      const size_t before = PJ.size();
      const size_t seg0 = SC.size();
      const uint32_t exlo = ex.empty() ? 1 : ex.front(), exhi = ex.empty() ? 0 : ex.back();
      for (const Cand& c : stencil) {
        const size_t code = c.code;
        const int q[3] = {qint_[3 * code], qint_[3 * code + 1], qint_[3 * code + 2]};
        const bool zero = q[0] == 0 && q[1] == 0 && q[2] == 0;
        const Vec3 xi = cw[i] - (cell.a * double(q[0]) + cell.b * double(q[1]) + cell.c * double(q[2]));
        const size_t pairs_before = PJ.size();
        for (size_t j = bstart[c.bin]; j < bstart[c.bin + 1]; ++j) {
          // Half list balanced over slots: an unordered pair is kept from the side given by the parity of i + j.
          if (j != i && (j > i) != (((i + j) & 1) == 0)) continue;
          if (j == i) {
            // self images: keep one of each ±shift pair
            if (zero) continue;
            if (q[0] < 0 || (q[0] == 0 && (q[1] < 0 || (q[1] == 0 && q[2] < 0)))) continue;
          }
          const Vec3 d = cw[j] - xi;
          if (dot(d, d) >= rcs2) continue;
          const uint32_t aj = perm_[j];
          // 1-2, 1-3 and 1-4 partners are excluded at their nearest image only (the bonded one, whatever the bond
          // length: coarse-grained bonds are several Å); further images interact, as LAMMPS's minimum-image check
          if (j != i && aj >= exlo && aj <= exhi && std::binary_search(ex.begin(), ex.end(), aj)) {
            Vec3 m = fr[j] - fr[i];
            for (int k = 0; k < 3; ++k)
              if (cell.periodic[k]) m[k] -= std::round(m[k]);
            const Vec3 dm = cell.a * m[0] + cell.b * m[1] + cell.c * m[2];
            if (dot(d, d) <= dot(dm, dm) + 1e-6) continue;
          }
          if (skipping && skip_type_[size_t(tslot_[i]) * ntypes + tslot_[j]]) continue;
          PJ.push_back(static_cast<uint32_t>(j));
          PQ.push_back(c.code);
        }
        if (PJ.size() == pairs_before) continue;
        if (SC.size() > seg0 && SC.back() == c.code) SE.back() = static_cast<uint32_t>(PJ.size() - before);
        else SC.push_back(c.code), SE.push_back(static_cast<uint32_t>(PJ.size() - before));
      }
      count[i] = static_cast<uint32_t>(PJ.size() - before);
      nseg[i] = static_cast<uint32_t>(SC.size() - seg0);
    }
  });
  first_.resize(n + 1);
  segoff_.resize(n + 1);
  first_[0] = 0;
  segoff_[0] = 0;
  for (size_t i = 0; i < n; ++i) first_[i + 1] = first_[i] + count[i], segoff_[i + 1] = segoff_[i] + nseg[i];
  pj_.resize(first_[n]);
  pq_.resize(first_[n]);
  segcode_.resize(segoff_[n]);
  segend_.resize(segoff_[n]);
  // every worker's slots are a contiguous range listed in order: its pairs go in one piece (the same split as above)
  pool_->run(n, [&](int t, size_t b, size_t) {
    std::copy(vj[t].begin(), vj[t].end(), pj_.begin() + std::ptrdiff_t(first_[b]));
    std::copy(vq[t].begin(), vq[t].end(), pq_.begin() + std::ptrdiff_t(first_[b]));
    std::copy(vsc[t].begin(), vsc[t].end(), segcode_.begin() + std::ptrdiff_t(segoff_[b]));
    std::copy(vse[t].begin(), vse[t].end(), segend_.begin() + std::ptrdiff_t(segoff_[b]));
  });
  // Each slot's offset to its binned position: a lattice vector (in the periodic directions) whose integer
  // coefficients are kept, so a change of cell (barostat, compression) only rescales the offsets and images instead
  // of forcing a rebuild.
  off_.resize(3 * n);
  ioff_.assign(3 * n, 0);
  for (size_t s = 0; s < n; ++s) {
    const uint32_t i = perm_[s];
    const Vec3 w = cell.origin + cell.a * fr[s][0] + cell.b * fr[s][1] + cell.c * fr[s][2];
    for (int k = 0; k < 3; ++k) off_[3 * s + k] = w[k] - x[3 * i + k];
  }
  if (cell_in.valid())
    for (size_t s = 0; s < n; ++s) {
      const Vec3 o = cell_in.to_fractional(cell_in.origin + Vec3{off_[3 * s], off_[3 * s + 1], off_[3 * s + 2]});
      for (int k = 0; k < 3; ++k) ioff_[3 * s + k] = static_cast<int32_t>(std::lround(o[k]));
    }
  qvec_.resize(qint_.size());
  for (size_t c = 0; c < qint_.size(); c += 3) {
    const Vec3 v = cell.a * qint_[c] + cell.b * qint_[c + 1] + cell.c * qint_[c + 2];
    for (int k = 0; k < 3; ++k) qvec_[c + k] = v[k];
  }
  x0_ = x;
  cell0_ = cell_in;
  cells_ = cell_in;
  built_ = true;
  ++builds_;
}

void Evaluator::rescale_images(const Cell& cell) {
  for (size_t i = 0; i < off_.size(); i += 3) {
    const Vec3 v = cell.a * ioff_[i] + cell.b * ioff_[i + 1] + cell.c * ioff_[i + 2];
    off_[i] = v[0]; off_[i + 1] = v[1]; off_[i + 2] = v[2];
  }
  for (size_t c = 0; c < qint_.size(); c += 3) {
    const Vec3 v = cell.a * qint_[c] + cell.b * qint_[c + 1] + cell.c * qint_[c + 2];
    qvec_[c] = v[0]; qvec_[c + 1] = v[1]; qvec_[c + 2] = v[2];
  }
  cells_ = cell;
}

const double* Evaluator::listed_positions(const std::vector<double>& x) {
  xs_.resize(x.size());
  for (size_t s = 0; s < perm_.size(); ++s) {
    const size_t i = 3 * size_t(perm_[s]);
    xs_[3 * s] = x[i] + off_[3 * s];
    xs_[3 * s + 1] = x[i + 1] + off_[3 * s + 1];
    xs_[3 * s + 2] = x[i + 2] + off_[3 * s + 2];
  }
  return xs_.data();
}

void Evaluator::ewald_table(double kappa, double s_max) {
  // nodes at most 0.025/κ² apart in s (and at least 1024 to the cut-off): the cubic's error is below 1e-12 of G
  const size_t m = std::max<size_t>(1024, size_t(std::ceil(s_max * kappa * kappa / 0.025)));
  const double h = s_max / double(m);
  if (kappa == gtab_kappa_ && h == gtab_h_ && !gtab_.empty()) return;
  const double c = 2 * kappa / std::sqrt(kPi), k2 = kappa * kappa;
  gtab_.assign(2 * (m + 2), 0.0);   // one node past s_max: a pair just inside the cut-off reads node m + 1
  for (size_t k = 0; k <= m + 1; ++k) {
    const double s = double(k) * h, x = k2 * s;
    double G, dG;
    if (x < 1) {
      // small κ²s, where the closed form cancels: G = c Σ (−x)ⁿ / (n! (2n + 1)), dG/ds = −c κ² Σ (−x)ⁿ / (n! (2n + 3))
      G = 0, dG = 0;
      double term = 1;   // (−x)ⁿ / n!
      for (int nn = 0; nn < 30; ++nn) {
        G += term / (2 * nn + 1);
        dG -= term / (2 * nn + 3);
        term *= -x / (nn + 1);
      }
      G *= c;
      dG *= c * k2;
    } else {
      const double r = std::sqrt(s), e = std::erf(kappa * r);
      G = e / r;
      dG = (c * std::exp(-x) - G) / (2 * s);
    }
    gtab_[2 * k] = G;
    gtab_[2 * k + 1] = h * dG;
  }
  gtab_kappa_ = kappa;
  gtab_h_ = h;
}

int apply_rigid_geometry(System& s, const ForceField& ff, size_t first, size_t last) {
  last = std::min(last, s.atoms.size());
  std::vector<std::vector<uint32_t>> nb(s.atoms.size());
  for (const auto& b : s.bonds)
    if (b.i >= first && b.i < last && b.j >= first && b.j < last) nb[b.i].push_back(b.j), nb[b.j].push_back(b.i);
  // the atoms reached from `from` without crossing the bond from–across; empty when the bond is in a ring
  auto side = [&](uint32_t from, uint32_t across) {
    std::vector<uint32_t> out{from}, todo{from};
    std::vector<char> seen(s.atoms.size(), 0);
    seen[from] = 1;
    while (!todo.empty()) {
      const uint32_t a = todo.back();
      todo.pop_back();
      for (uint32_t c : nb[a]) {
        if (a == from && c == across) continue;
        if (c == across) return std::vector<uint32_t>{};
        if (!seen[c]) seen[c] = 1, out.push_back(c), todo.push_back(c);
      }
    }
    return out;
  };
  auto in = [&](uint32_t a) { return a >= first && a < last; };
  int changed = 0;
  for (const auto& b : ff.bonds) {
    if (b.k != 0 || b.r0 <= 0 || !in(b.i) || !in(b.j)) continue;
    const Vec3 d = s.atoms[b.j].pos - s.atoms[b.i].pos;
    const double r = norm(d);
    if (r < 1e-9 || std::fabs(r - b.r0) < 1e-6) continue;
    const auto mv = side(b.j, b.i);
    if (mv.empty()) continue;
    const Vec3 shift = d * ((b.r0 - r) / r);
    for (uint32_t a : mv) s.atoms[a].pos = s.atoms[a].pos + shift;
    ++changed;
  }
  for (const auto& a : ff.angles) {
    if (a.kt != 0 || !in(a.i) || !in(a.j) || !in(a.k)) continue;
    const Vec3 u = s.atoms[a.i].pos - s.atoms[a.j].pos, w = s.atoms[a.k].pos - s.atoms[a.j].pos;
    const double nu = norm(u), nw = norm(w);
    if (nu < 1e-9 || nw < 1e-9) continue;
    const double th = std::acos(std::clamp(dot(u, w) / (nu * nw), -1.0, 1.0)), dth = a.theta0 - th;
    if (std::fabs(dth) < 1e-6) continue;
    Vec3 n = cross(u, w);
    if (norm(n) < 1e-9) {   // collinear now: any normal to the line
      n = cross(u, std::fabs(u[0]) < 0.9 * nu ? Vec3{1, 0, 0} : Vec3{0, 1, 0});
    }
    n = n * (1.0 / norm(n));
    const auto mv = side(a.k, a.j);
    if (mv.empty()) continue;
    // Rodrigues: turn by dth about n through the vertex (positive dth opens the angle)
    const double c = std::cos(dth), sn = std::sin(dth);
    const Vec3 o = s.atoms[a.j].pos;
    for (uint32_t m : mv) {
      const Vec3 v = s.atoms[m].pos - o;
      s.atoms[m].pos = o + v * c + cross(n, v) * sn + n * (dot(n, v) * (1 - c));
    }
    ++changed;
  }
  return changed;
}

void place_virtual_sites(const ForceField& ff, std::vector<double>& x, const Cell& cell) {
  for (const auto& v : ff.vsites) {
    const Vec3 ref{x[3 * v.from[0]], x[3 * v.from[0] + 1], x[3 * v.from[0] + 2]};
    Vec3 c{0, 0, 0};
    for (size_t k = 0; k < v.from.size(); ++k) {
      Vec3 d = Vec3{x[3 * v.from[k]], x[3 * v.from[k] + 1], x[3 * v.from[k] + 2]} - ref;
      if (cell.valid()) d = cell.minimum_image(d);
      c = c + (ref + d) * v.w[k];
    }
    if (v.c != 0 && v.from.size() == 3) {   // out of the plane: c (r_ij × r_ik)
      Vec3 u = Vec3{x[3 * v.from[1]], x[3 * v.from[1] + 1], x[3 * v.from[1] + 2]} - ref;
      Vec3 w = Vec3{x[3 * v.from[2]], x[3 * v.from[2] + 1], x[3 * v.from[2] + 2]} - ref;
      if (cell.valid()) u = cell.minimum_image(u), w = cell.minimum_image(w);
      c = c + cross(u, w) * v.c;
    }
    for (int k = 0; k < 3; ++k) x[3 * v.site + k] = c[k];
  }
}

EnergyTerms Evaluator::compute(const std::vector<double>& x, const Cell& cell, std::vector<double>& f) {
  if (ff_.vsites.empty()) return compute_placed(x, cell, f);
  xv_ = x;
  place_virtual_sites(ff_, xv_, cell);
  EnergyTerms e = compute_placed(xv_, cell, f);
  // the site's force goes back to its constructing atoms (the virial is unchanged by a linear combination);
  // in reverse order, so a site built on other sites hands its force to them before they hand theirs on
  for (auto it = ff_.vsites.rbegin(); it != ff_.vsites.rend(); ++it) {
    const auto& v = *it;
    const Vec3 F{f[3 * v.site], f[3 * v.site + 1], f[3 * v.site + 2]};
    for (size_t k = 0; k < v.from.size(); ++k)
      for (int c = 0; c < 3; ++c) f[3 * v.from[k] + c] += v.w[k] * F[c];
    if (v.c != 0 && v.from.size() == 3) {
      // s = c (u × w), u = x_j − x_i, w = x_k − x_i: F_u = c (w × F), F_w = c (F × u), and −(F_u + F_w) on i
      const auto at = [&](uint32_t a) { return Vec3{xv_[3 * a], xv_[3 * a + 1], xv_[3 * a + 2]}; };
      Vec3 u = at(v.from[1]) - at(v.from[0]), w = at(v.from[2]) - at(v.from[0]);
      if (cell.valid()) u = cell.minimum_image(u), w = cell.minimum_image(w);
      const Vec3 fu = cross(w, F) * v.c, fw = cross(F, u) * v.c;
      for (int c = 0; c < 3; ++c) f[3 * v.from[1] + c] += fu[c], f[3 * v.from[2] + c] += fw[c], f[3 * v.from[0] + c] -= fu[c] + fw[c];
      // the cross part grows as λ² when the atoms are scaled by λ: one more s ⊗ F in the virial than the placed site counted
      const Vec3 s = cross(u, w) * v.c;
      e.virial += dot(s, F);
      e.w[0] += s[0] * F[0], e.w[1] += s[1] * F[1], e.w[2] += s[2] * F[2];
      e.w[3] += 0.5 * (s[0] * F[1] + s[1] * F[0]), e.w[4] += 0.5 * (s[0] * F[2] + s[2] * F[0]), e.w[5] += 0.5 * (s[1] * F[2] + s[2] * F[1]);
    }
    for (int c = 0; c < 3; ++c) f[3 * v.site + c] = 0;
  }
  return e;
}

EnergyTerms Evaluator::compute_placed(const std::vector<double>& x, const Cell& cell, std::vector<double>& f) {
  const size_t n = x.size() / 3;
  f.assign(x.size(), 0.0);
  EnergyTerms e;

  // r-RESPA splits the forces: bonded terms alone need no pair list
  const bool nonb = opt_.parts & 2, bonded = opt_.parts & 1;
  const bool clmode = nonb && cluster_mode();
  if (clmode) {
    if (!cl_) cl_ = std::make_unique<nbnxm::List>();
    if (!cl_coef_) cluster_coefficients();
    const int st = cl_built_ ? list_state(x, cell, x0c_, cell0c_, cl_->cells, opt_.skin, opt_.cutoff) : 2;
    if (st == 2) {
      const std::vector<char> none;
      cl_->build(x, cell, opt_.cutoff + opt_.skin, ff_.excluded, ff_.type_index, ff_.lj.size(),
                 ff_.excluded_type_pairs.empty() ? none : skip_type_, qeff_, cl_f12_, cl_f6_, *pool_);
      x0c_ = x;
      cell0c_ = cell;
      cl_built_ = true;
      ++builds_;
    } else if (st == 1) {
      cl_->rescale(cell);
    }
    cl_->positions(x);
    // Dynamic pruning: the kernel takes the pairs within r_c + b of the outer list, pruned again once an atom has moved
    // b/2 since the last prune (b = 0.8 Å, or the skin when that is shorter)
    const double b = std::min(0.8, opt_.skin);
    if (st == 2 || !cl_pruned_ || list_state(x, cell, x0p_, cell0p_, cell, b, opt_.cutoff) == 2) {
      cl_->prune(opt_.cutoff + b, *pool_);
      x0p_ = x;
      cell0p_ = cell;
      cl_pruned_ = true;
    }
  } else if (nonb) {
    auto same = [](const Cell& p, const Cell& q) { return p.a == q.a && p.b == q.b && p.c == q.c && p.origin == q.origin; };
    bool rebuild = !built_ || x0_.size() != x.size();
    if (frozen_ && !rebuild) {
      if (!same(cell, cells_) && cell.valid()) rescale_images(cell);
    } else if (!rebuild) {
      const int st = list_state(x, cell, x0_, cell0_, cells_, opt_.skin, opt_.cutoff);
      rebuild = st == 2;
      if (st == 1) rescale_images(cell);
    }
    if (rebuild) {
      if (frozen_) throw std::logic_error("frozen pair list: the configuration changed size");
      build(x, cell);
    }
  }

  const bool ortho = cell.valid() && orthogonal(cell);
  const double L[3] = {cell.a[0], cell.b[1], cell.c[2]};
  const bool cvalid = cell.valid();   // a volume: once, not per call
  auto mi = [&](Vec3 d) {
    if (!cvalid) return d;
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
  // MARTINI: Coulomb / εr with GROMACS's force switch (coul/gromacs) replaces DSF / PME
  const bool gro = coul && ff_.coul_gromacs;
  const double ri = ff_.coul_inner, qscale = 1.0;   // εr is in the scaled charges
  double gc[5] = {0, 0, 0, 0, 0};
  if (gro) {
    const double r3i = 1 / (rc2 * rc), tt = rc - ri, t2i = 1 / (tt * tt), t3i = t2i / tt;
    const double a1 = (2 * ri - 5 * rc) * r3i * t2i, b1 = (4 * rc - 2 * ri) * r3i * t3i;
    gc[0] = a1; gc[1] = b1; gc[2] = -a1 / 3; gc[3] = -b1 / 4; gc[4] = 1 / rc - tt * tt * tt * (a1 / 3 + b1 * tt / 4);
  }
  // Martini 3: reaction field, 1/r + k_rf r² − c_rf, k_rf = (ε_rf − ε_r) / ((2 ε_rf + ε_r) r_c³) (ε_rf infinite: 1 / (2 r_c³)); the
  // relative permittivity ε_r is in the scaled charges
  const bool rf = coul && !gro && ff_.coul_rf;
  double krf = 0, crf = 0;
  if (rf) {
    const double er = ff_.dielectric, erf = ff_.eps_rf;
    krf = erf <= 0 ? 1 / (2 * rc2 * rc) : (erf - er) / ((2 * erf + er) * rc2 * rc);
    crf = 1 / rc + krf * rc2;
  }
  const bool pme = coul && !gro && !rf && opt_.electrostatics == EnergyOptions::Electrostatics::PME && cell.valid() && cell.periodic[0] && cell.periodic[1] &&
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
    if (form_[tp] == kPairCos2 || form_[tp] == kPairCos2Wca) {   // cosine-squared attraction, optionally with WCA
      const double r = std::sqrt(r2), ep = pa_[tp], sg = pb_[tp], rcp = pc_[tp];
      if (r >= rcp && !(form_[tp] == kPairCos2Wca && r < sg)) return 0.0;
      double ev = 0, dEdr = 0;
      if (r < sg) {
        if (rcp > sg) ev = -ep;   // the attraction's flat bottom (none when the pair is WCA only)
        if (form_[tp] == kPairCos2Wca) {
          const double s6 = std::pow(sg / r, 6);
          ev += ep * (s6 * s6 - 2 * s6 + 1);
          dEdr += -12 * ep * (s6 * s6 - s6) / r;
        }
      } else {
        const double w = kPi * (r - sg) / (2 * (rcp - sg)), c = std::cos(w);
        ev = -ep * c * c;
        dEdr = ep * 2 * c * std::sin(w) * kPi / (2 * (rcp - sg));
      }
      en += scale * ev;
      return -scale * dEdr / r;
    }
    if (form_[tp] == kPairTable) {   // tabulated: E and F = −dE/dr linearly interpolated
      const auto& tb = ff_.tables[size_t(pa_[tp])];
      const double r = std::sqrt(r2);
      const size_t m = tb.e.size();
      if (m < 2 || r >= tb.r0 + tb.dr * double(m - 1)) return 0.0;
      if (r < tb.r0) {   // below the table: a straight wall at the first force
        en += scale * (tb.e[0] + tb.f[0] * (tb.r0 - r));
        return scale * tb.f[0] / r;
      }
      const double u = (r - tb.r0) / tb.dr;
      const size_t k = std::min(m - 2, size_t(u));
      const double w = u - double(k);
      en += scale * (tb.e[k] * (1 - w) + tb.e[k + 1] * w);
      return scale * (tb.f[k] * (1 - w) + tb.f[k + 1] * w) / r;
    }
    if (form_[tp] == kPairGromacs) {   // LJ with GROMACS's force switch (lj/gromacs)
      const double e = pa_[tp], s6 = std::pow(pb_[tp], 6), r2i = 1 / r2, r6i = r2i * r2i * r2i;
      const double* g = &gsw_[5 * tp];
      double ev = r6i * (4 * e * s6 * s6 * r6i - 4 * e * s6) + g[4], fr = r6i * (48 * e * s6 * s6 * r6i - 24 * e * s6);
      if (r2 > ff_.lj_inner * ff_.lj_inner) {
        const double r = std::sqrt(r2), tl = r - ff_.lj_inner;
        ev += tl * tl * tl * (g[2] + g[3] * tl);
        fr += r * tl * tl * (g[0] + g[1] * tl);
      }
      en += scale * ev;
      return scale * fr * r2i;
    }
    if (form_[tp] >= kPairSdk96 && form_[tp] <= kPairSdk125) {   // SDK / SPICA C ε [(σ/r)^m − (σ/r)^n]
      static const double pw1[] = {9, 12, 12, 12}, pw2[] = {6, 4, 6, 5}, pre[] = {6.75, 2.59807621135332, 4.0, 3.20377984125109};
      const int k = form_[tp] - kPairSdk96;
      const double sr2 = pb_[tp] * pb_[tp] / r2, sm = std::pow(sr2, pw1[k] / 2), sn = std::pow(sr2, pw2[k] / 2), c = pre[k] * pa_[tp];
      en += scale * c * (sm - sn);
      return scale * c * (pw1[k] * sm - pw2[k] * sn) / r2;
    }
    if (capped && r2 < rcap2_[tp]) {
      const double r = std::sqrt(r2), rcp = std::sqrt(rcap2_[tp]);
      en += scale * (ecap_[tp] + fcap_[tp] * (rcp - r));
      return scale * fcap_[tp] / r;
    }
    if (ff_.lj_fsw && !lj96_) {
      // CHARMM force switch: below r_in the plain 12-6 plus a constant, between r_in and r_c
      //   E = A k12 (r⁻⁶ − r_c⁻⁶)² − B k6 (r⁻³ − r_c⁻³)²,  k12 = r_c⁶/(r_c⁶ − r_in⁶), k6 = r_c³/(r_c³ − r_in³)
      // (A = 4εσ¹², B = 4εσ⁶); the force is its exact derivative
      const double A = 4 * eps * s6_[tp] * s6_[tp], B = 4 * eps * s6_[tp], ri = ff_.lj_inner;
      const double ri3 = ri * ri * ri, rc3 = rc2 * rc, ri6 = ri3 * ri3, rc6 = rc3 * rc3, r6i = 1 / (r2 * r2 * r2);
      if (r2 <= ri * ri) {
        en += scale * (A * (r6i * r6i - 1 / (ri6 * rc6)) - B * (r6i - 1 / (ri3 * rc3)));
        return scale * (12 * A * r6i * r6i - 6 * B * r6i) / r2;
      }
      const double r3i = 1 / (r2 * std::sqrt(r2)), k12 = rc6 / (rc6 - ri6), k6 = rc3 / (rc3 - ri3);
      const double u = r6i - 1 / rc6, w = r3i - 1 / rc3;
      en += scale * (A * k12 * u * u - B * k6 * w * w);
      return scale * (12 * A * k12 * u * r6i - 6 * B * k6 * w * r3i) / r2;
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
    if (opt_.tail || ff_.lj_fsw || form_[tp] == kPairGromacs || form_[tp] == kPairCos2 || form_[tp] == kPairCos2Wca || form_[tp] == kPairTable) return 0.0;   // zero at their cut-offs by construction
    if (form_[tp] != 0) {
      double e0 = 0;
      lj(tp, rc2, 1.0, e0);
      return eps_[tp] > 0 ? e0 : 0.0;
    }
    const double q = s6_[tp] / (rc2 * rc2 * rc2);
    if (lj96_) return eps_[tp] * (2 * q * std::sqrt(q) - 3 * q);
    return 4 * eps_[tp] * (q * q - q);
  };

  // DREIDING hydrogen bond for a donor D, acceptor A and hydrogen H: d1 = D − H, d2 = A − H (minimum images).
  // E = S(r) ε [5 (σ/r)¹² − 6 (σ/r)¹⁰] cosⁿθ (r = |A − D|, θ the D–H···A angle), only for cos θ < cos of the cut-off; S the
  // switch from the inner to the outer radius (LAMMPS hbond/dreiding/lj). Forces are the exact derivatives.
  struct HbOut { double e; Vec3 fd, fa, fh; };
  const bool hbon = ff_.hbond.on() && nonb;
  const double hro2 = ff_.hbond.outer * ff_.hbond.outer;
  auto hb_one = [&](const std::array<double, 3>& pr, const Vec3& d1, const Vec3& d2, HbOut& o) -> bool {
    const auto& Hb = ff_.hbond;
    const Vec3 dA = d2 - d1;   // A − D
    const double r2 = dot(dA, dA);
    if (r2 >= hro2) return false;
    const double l1 = std::sqrt(dot(d1, d1)), l2 = std::sqrt(dot(d2, d2));
    if (l1 < 1e-12 || l2 < 1e-12) return false;
    const double c = std::clamp(dot(d1, d2) / (l1 * l2), -1.0, 1.0);
    if (!(c < Hb.cos_cut)) return false;
    const double eps = pr[0], sg = pr[1];
    const int nn = int(pr[2]);
    const double r = std::sqrt(r2), sr2 = sg * sg / r2, s10 = sr2 * sr2 * sr2 * sr2 * sr2, s12 = s10 * sr2;
    const double U = eps * (5 * s12 - 6 * s10), dU = eps * (-60 * s12 + 60 * s10) / r;
    double S = 1, dS = 0;
    const double ri2 = Hb.inner * Hb.inner;
    if (r2 > ri2) {
      const double den = (hro2 - ri2) * (hro2 - ri2) * (hro2 - ri2);
      S = (hro2 - r2) * (hro2 - r2) * (hro2 + 2 * r2 - 3 * ri2) / den;
      dS = -12 * r * (hro2 - r2) * (r2 - ri2) / den;
    }
    double cn = 1, cn1 = 0;
    for (int k = 0; k < nn; ++k) { cn1 = cn; cn *= c; }
    const double dEdr = (dS * U + S * dU) * cn, dEdc = S * U * nn * cn1;
    const Vec3 g1 = d2 * (1 / (l1 * l2)) - d1 * (c / (l1 * l1)), g2 = d1 * (1 / (l1 * l2)) - d2 * (c / (l2 * l2)), u = dA * (1 / r);
    o.e = S * U * cn;
    o.fd = u * dEdr - g1 * dEdc;      // x_D moves d1 (+) and r (along −u)
    o.fa = u * (-dEdr) - g2 * dEdc;   // x_A moves d2 (+) and r (along +u)
    o.fh = (g1 + g2) * dEdc;          // x_H moves d1 and d2 (−)
    return true;
  };

  const double* q = qeff_.data();
  const int nth = pool_->size();
  tf_.resize(nth);
  tfs_.resize(nth);
  // per worker: bond, angle, dihedral, improper, vdW, Coulomb, virial, virial tensor (xx yy zz xy xz yz)
  std::vector<std::array<double, 15>> acc(nth);   // + ∂U/∂λ of the alchemical pairs (LJ, Coulomb)
  // alchemical pairs: plain 12-6 for every type an alchemical atom has, pairwise electrostatics, no tail
  const std::vector<char>* alch = opt_.alchemical && opt_.alchemical->size() == n ? opt_.alchemical : nullptr;
  if (alch && nonb) {
    if (pme) throw FieldError("alchemical pairs need pairwise electrostatics (DSF or reaction field), not PME");
    if (opt_.tail) throw FieldError("alchemical pairs need the tail correction off");
    if (opt_.force_cap > 0 || ff_.lj_fsw || lj96_) throw FieldError("alchemical pairs take the plain 12-6 Lennard-Jones only");
    std::vector<char> at(nt, 0);
    for (size_t i = 0; i < n; ++i) if ((*alch)[i]) at[size_t(ff_.type_index[i])] = 1;
    for (size_t a = 0; a < nt; ++a)
      if (at[a])
        for (size_t b = 0; b < nt; ++b)
          if (form_[a * nt + b] != 0) throw FieldError("alchemical pairs take the plain 12-6 Lennard-Jones only (" + ff_.type_names[a] + " has another pair form)");
  }
  const double lam_lj = opt_.lambda_lj, lam_c = opt_.lambda_coul, sca = opt_.sc_alpha;
  // the soft-core pair: energy, −(dU/dr)/r and ∂U/∂λ at r² (σ⁶ = s6, ε)
  auto softcore = [&](double r2, double s6, double eps, double& u, double& fr_out, double& dudl) {
    const double r6 = r2 * r2 * r2, rs6 = sca * s6 * (1 - lam_lj) + r6, x = s6 / rs6;
    u = lam_lj * 4 * eps * (x * x - x);
    fr_out = 24 * eps * lam_lj * (2 * x - 1) * x * r6 / (rs6 * r2);   // −(dU/dr)/r with dr_sc⁶/dr = 6r⁵
    dudl = 4 * eps * (x * x - x) + lam_lj * 4 * eps * (2 * x - 1) * x * sca * s6 / rs6;
  };
  for (auto& v : acc) v.fill(0.0);

  const bool frozen = frozen_;
  // Phase 1: pair list. Every worker is called (possibly with an empty range) and clears its own force buffers; its
  // slots are a contiguous range holding about 1/nth of the pairs. Pair forces go to a buffer in slot order.
  const double* xs = nonb && !clmode ? listed_positions(x) : nullptr;
  const size_t npairs = nonb && !clmode ? pj_.size() : 0;
  // Plain 12-6 pairs (no other pair form, cap, switch, alchemy, hydrogen bonds or frozen list; PME, DSF or no
  // charges) take a copy of the loop with those branches compiled out, the image shift taken once per segment and
  // the two coefficients per type pair: 12-6 E = c12 r⁻¹² − cc6 r⁻⁶ (4ε σ¹², 4ε σ⁶); class II 9-6 E = c12 r⁻⁹ − cc6 r⁻⁶
  // (2ε σ⁹, 3ε σ⁶).
  const bool plain = nonb && !alch && !hbon && !frozen_ && !capped && !ff_.lj_fsw && !rf && !gro &&
                     std::all_of(form_.begin(), form_.end(), [](uint8_t v) { return v == 0; });
  std::vector<double> shift0, c12, cc6;
  if (plain) {
    shift0.resize(nt * nt);
    c12.resize(nt * nt);
    cc6.resize(nt * nt);
    for (size_t tp = 0; tp < shift0.size(); ++tp) {
      shift0[tp] = lj_shift(tp);
      const double eps = eps_[tp] > 0 ? eps_[tp] : 0.0;
      c12[tp] = lj96_ ? 2 * eps * s6_[tp] * std::sqrt(s6_[tp]) : 4 * eps * s6_[tp] * s6_[tp];
      cc6[tp] = lj96_ ? 3 * eps * s6_[tp] : 4 * eps * s6_[tp];
    }
  }
  // erfc(κr)/r (E) and −(dE/dr)/r (F) at r² = s, 1/r given: from the table inside the cut-off
  const double kappa = pme ? beta : a;
  if (nonb && coul && !rf && !gro) ewald_table(kappa, rc2);
  const double* gt = gtab_.data();
  const double ghi = gtab_h_ > 0 ? 1 / gtab_h_ : 0.0;
  auto ewald_ef = [gt, ghi](double s, double rinv, double& E, double& F) {
    const double u = s * ghi;
    const size_t k = size_t(u);
    const double t = u - double(k), t2 = t * t;
    const double* g = gt + 2 * k;   // G_k, h G'_k, G_k+1, h G'_k+1
    const double G = g[0] + t2 * (3 - 2 * t) * (g[2] - g[0]) + t * (1 - t) * ((1 - t) * g[1] - t * g[3]);
    const double dG = (6 * t * (1 - t) * (g[2] - g[0]) * ghi) + (1 - t) * (1 - 3 * t) * g[1] * ghi + t * (3 * t - 2) * g[3] * ghi;
    E = rinv - G;
    F = rinv * rinv * rinv + 2 * dG;
  };
  // kind: 0 any pair, 1 plain 12-6, 2 plain 9-6
  auto pair_phase = [&](int t, auto kind_c) {
    constexpr int kKind = decltype(kind_c)::value;
    constexpr bool kPlain = kKind != 0;
    std::vector<double>& ft = tf_[t];
    ft.assign(x.size(), 0.0);
    std::vector<double>& fs = tfs_[t];
    if (nonb) fs.assign(x.size(), 0.0);
    double evdw = 0, ecoul = 0, vir = 0, wv[6] = {0, 0, 0, 0, 0, 0}, dl_lj = 0, dl_c = 0;
    if (npairs == 0) return;
    const size_t ib = size_t(std::lower_bound(first_.begin(), first_.end() - 1, npairs * size_t(t) / size_t(nth)) - first_.begin());
    const size_t ie = size_t(std::lower_bound(first_.begin(), first_.end() - 1, npairs * size_t(t + 1) / size_t(nth)) - first_.begin());
    if constexpr (kPlain) {
      const uint32_t* ts = tslot_.data();
      const double* qs = qslot_.data();
      double* F = fs.data();
      const bool pmec = coul && pme, dsfc = coul && !pme;
      // Two passes per slot: the separations of its listed pairs, keeping those inside the cut-off without a branch,
      // then the interactions of those alone (a third or more of a list with a 2 Å skin is outside, at random).
      size_t most = 0;
      for (size_t si = ib; si < ie; ++si) most = std::max(most, first_[si + 1] - first_[si]);
      std::vector<uint32_t> bj(most + 1);
      std::vector<double> bd(4 * (most + 1));
      for (size_t si = ib; si < ie; ++si) {
        size_t m = 0;
        for (size_t k = segoff_[si]; k < segoff_[si + 1]; ++k) {
          const double* sv = &qvec_[3 * size_t(segcode_[k])];
          // d = xs_j − (xs_i − shift)
          const double xi = xs[3 * si] - sv[0], yi = xs[3 * si + 1] - sv[1], zi = xs[3 * si + 2] - sv[2];
          const size_t pb = k == segoff_[si] ? first_[si] : first_[si] + segend_[k - 1], pe = first_[si] + segend_[k];
          for (size_t p = pb; p < pe; ++p) {
            const uint32_t sj = pj_[p];
            const double dx = xs[3 * sj] - xi, dy = xs[3 * sj + 1] - yi, dz = xs[3 * sj + 2] - zi;
            const double r2 = dx * dx + dy * dy + dz * dz;
            bj[m] = sj;
            double* d = &bd[4 * m];
            d[0] = dx; d[1] = dy; d[2] = dz; d[3] = r2;
            m += r2 < rc2;
          }
        }
        double fxi = 0, fyi = 0, fzi = 0;
        const size_t ti = size_t(ts[si]) * nt;
        const double qi = kCoulomb * qs[si];
        for (size_t c = 0; c < m; ++c) {
          const uint32_t sj = bj[c];
          const double* d = &bd[4 * c];
          const double dx = d[0], dy = d[1], dz = d[2], r2 = d[3];
          const size_t tp = ti + ts[sj];
          const double r2i = 1 / r2, r6i = r2i * r2i * r2i, rinv = std::sqrt(r2i);
          double fr;
          if constexpr (kKind == 2) {
            const double r9i = r6i * r2i * rinv;
            evdw += c12[tp] * r9i - cc6[tp] * r6i - shift0[tp];
            fr = (9 * c12[tp] * r9i - 6 * cc6[tp] * r6i) * r2i;
          } else {
            evdw += r6i * (c12[tp] * r6i - cc6[tp]) - shift0[tp];
            fr = r6i * (12 * c12[tp] * r6i - 6 * cc6[tp]) * r2i;
          }
          if (pmec || dsfc) {
            const double qq = qi * qs[sj];
            double E, Fc;
            ewald_ef(r2, rinv, E, Fc);
            if (pmec) {
              ecoul += qq * E;
              fr += qq * Fc;
            } else {
              ecoul += qq * (E - dsf_e0 + dsf_f0 * (r2 * rinv - rc));
              fr += qq * (Fc - dsf_f0 * rinv);
            }
          }
          const double fx = fr * dx, fy = fr * dy, fz = fr * dz;
          F[3 * sj] += fx; F[3 * sj + 1] += fy; F[3 * sj + 2] += fz;
          fxi -= fx; fyi -= fy; fzi -= fz;
          wv[0] += fx * dx; wv[1] += fy * dy; wv[2] += fz * dz;
          wv[3] += fx * dy; wv[4] += fx * dz; wv[5] += fy * dz;   // pair forces are central: f ∥ d, so f_a d_b = f_b d_a
        }
        F[3 * si] += fxi; F[3 * si + 1] += fyi; F[3 * si + 2] += fzi;
      }
      vir = wv[0] + wv[1] + wv[2];
    } else {
    // the same two passes: the pairs that interact (inside the cut-off, or frozen in), then their terms
    size_t most = 0;
    for (size_t si = ib; si < ie; ++si) most = std::max(most, first_[si + 1] - first_[si]);
    std::vector<size_t> bp(most + 1);
    std::vector<double> bd(4 * (most + 1));
    for (size_t si = ib; si < ie; ++si) {
      const uint32_t i = perm_[si];
      const double xi = xs[3 * si], yi = xs[3 * si + 1], zi = xs[3 * si + 2];
      size_t m = 0;
      for (size_t p = first_[si]; p < first_[si + 1]; ++p) {
        const uint32_t sj = pj_[p];
        const double* sv = &qvec_[3 * size_t(pq_[p])];
        const double dx = xs[3 * sj] - xi + sv[0], dy = xs[3 * sj + 1] - yi + sv[1], dz = xs[3 * sj + 2] - zi + sv[2];
        const double r2 = dx * dx + dy * dy + dz * dz;
        bp[m] = p;
        double* d = &bd[4 * m];
        d[0] = dx; d[1] = dy; d[2] = dz; d[3] = r2;
        m += frozen ? inside_[p] != 0 : r2 < rc2;
      }
      double fxi = 0, fyi = 0, fzi = 0;
      const size_t ti = size_t(ff_.type_index[i]) * nt;
      for (size_t c = 0; c < m; ++c) {
        const size_t p = bp[c];
        const uint32_t sj = pj_[p], j = perm_[sj];
        const double* d = &bd[4 * c];
        const double dx = d[0], dy = d[1], dz = d[2], r2 = d[3];
        // a self image is listed once for each ±shift pair (the pair list keeps one), so it counts in full, as in LAMMPS
        const double w = 1.0;
        const size_t tp = ti + ff_.type_index[j];
        const bool ab = alch && (*alch)[i] != (*alch)[j];   // an alchemical pair: soft-core, scaled
        double ev = 0, fr = 0;
        if (ab) {
          if (eps_[tp] > 0) {
            double u, fsc, dl, uc, fcc, dlc;
            softcore(r2, s6_[tp], eps_[tp], u, fsc, dl);
            softcore(rc2, s6_[tp], eps_[tp], uc, fcc, dlc);   // shifted to zero at the cut-off, as the plain pairs
            ev = u - uc;
            fr = fsc;
            dl_lj += dl - dlc;
          }
        } else {
          fr = lj(tp, r2, 1.0, ev);
          ev -= lj_shift(tp);
        }
        evdw += w * ev;
        const double ecoul0 = ecoul, fr0 = fr;
        if (coul && q[i] != 0 && q[j] != 0) {
          const double qq = kCoulomb * q[i] * q[j];
          if (rf) {
            const double r = std::sqrt(r2);
            ecoul += w * qq * (1 / r + krf * r2 - crf);
            fr += qq * (1 / (r2 * r) - 2 * krf);
          } else if (gro) {
            const double r = std::sqrt(r2);
            const double qd = qq * qscale;
            double ec = qd * (1 / r - gc[4]), fc = qd / r;
            if (r > ri) {
              const double tc = r - ri;
              ec += qd * tc * tc * tc * (gc[2] + gc[3] * tc);
              fc += qd * r * tc * tc * (gc[0] + gc[1] * tc);
            }
            ecoul += w * ec;
            fr += fc / r2;
          } else if (r2 < rc2) {
            const double rinv = 1 / std::sqrt(r2);
            double E, Fc;
            ewald_ef(r2, rinv, E, Fc);
            if (pme) {
              ecoul += w * qq * E;
              fr += qq * Fc;
            } else {
              ecoul += w * qq * (E - dsf_e0 + dsf_f0 * (r2 * rinv - rc));
              fr += qq * (Fc - dsf_f0 * rinv);
            }
          } else {   // a frozen pair beyond the cut-off: the closed form
            const double r = std::sqrt(r2);
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
        }
        if (ab && ecoul != ecoul0) {   // the pair's electrostatics × λ_coul, its full value is ∂U/∂λ_coul
          const double de = ecoul - ecoul0;
          dl_c += de;
          ecoul = ecoul0 + lam_c * de;
          fr = fr0 + lam_c * (fr - fr0);
        }
        if (hbon && i != j && r2 < hro2) {   // hydrogen bonds, either atom the donor
          auto run = [&](uint32_t D, uint32_t Aa, const Vec3& dDA) {
            const auto it = ff_.hbond.param.find({ff_.type_index[D], ff_.type_index[Aa]});
            if (it == ff_.hbond.param.end()) return;
            for (uint32_t h : ff_.hbond.hyd[D]) {
              const Vec3 d1 = mi(pos(D) - pos(h)), d2 = dDA + d1;
              HbOut o;
              if (!hb_one(it->second, d1, d2, o)) continue;
              evdw += o.e;
              const size_t sD = iperm_[D], sA = iperm_[Aa], sH = iperm_[h];
              for (int c = 0; c < 3; ++c) { fs[3 * sD + c] += o.fd[c]; fs[3 * sA + c] += o.fa[c]; fs[3 * sH + c] += o.fh[c]; }
              vir += dot(d1, o.fd) + dot(d2, o.fa);
              wv[0] += d1[0] * o.fd[0] + d2[0] * o.fa[0]; wv[1] += d1[1] * o.fd[1] + d2[1] * o.fa[1]; wv[2] += d1[2] * o.fd[2] + d2[2] * o.fa[2];
              wv[3] += 0.5 * (d1[0] * o.fd[1] + d1[1] * o.fd[0] + d2[0] * o.fa[1] + d2[1] * o.fa[0]);
              wv[4] += 0.5 * (d1[0] * o.fd[2] + d1[2] * o.fd[0] + d2[0] * o.fa[2] + d2[2] * o.fa[0]);
              wv[5] += 0.5 * (d1[1] * o.fd[2] + d1[2] * o.fd[1] + d2[1] * o.fa[2] + d2[2] * o.fa[1]);
            }
          };
          if (!ff_.hbond.hyd[i].empty() && ff_.hbond.acceptor[j]) run(i, j, Vec3{dx, dy, dz});
          if (!ff_.hbond.hyd[j].empty() && ff_.hbond.acceptor[i]) run(j, i, Vec3{-dx, -dy, -dz});
        }
        fr *= w;
        const double fx = fr * dx, fy = fr * dy, fz = fr * dz;
        fs[3 * sj] += fx; fs[3 * sj + 1] += fy; fs[3 * sj + 2] += fz;
        fxi -= fx; fyi -= fy; fzi -= fz;
        vir += fx * dx + fy * dy + fz * dz;
        wv[0] += fx * dx; wv[1] += fy * dy; wv[2] += fz * dz;
        wv[3] += fx * dy; wv[4] += fx * dz; wv[5] += fy * dz;   // pair forces are central: f ∥ d, so f_a d_b = f_b d_a
      }
      fs[3 * si] += fxi; fs[3 * si + 1] += fyi; fs[3 * si + 2] += fzi;
    }
    }
    acc[t][4] += evdw;
    acc[t][5] += ecoul;
    acc[t][6] += vir;
    for (int c = 0; c < 6; ++c) acc[t][7 + c] += wv[c];
    acc[t][13] += dl_lj;
    acc[t][14] += dl_c;
  };

  // Mixed precision: the cluster kernel over each worker's i clusters (about 1/nth of the cluster pairs), float forces
  // by slot
  nbnxm::Params cp;
  if (clmode) {
    tfc_.resize(nth);
    cp.rc2 = float(rc2);
    cp.rc = float(rc);
    cp.lj = lj96_ ? 2 : 1;
    cp.comb = !cl_f12_.empty();
    cp.c12 = cl_c12_.data();
    cp.c6 = cl_c6_.data();
    cp.nt = nt;
    cp.shift = opt_.tail ? 0.0f : 1.0f;   // as lj_shift: truncated with the tail correction, shifted without
    cp.rcn = float(lj96_ ? std::pow(rc, -9.0) : std::pow(rc, -12.0));
    cp.rc6 = float(std::pow(rc, -6.0));
    cp.coul = coul && !rf && !gro ? (pme ? 1 : 2) : 0;
    cp.kappa = float(pme ? beta : a);
    cp.kc = float(2 * (pme ? beta : a) / std::sqrt(kPi));
    cp.e0 = float(dsf_e0);
    cp.f0 = float(dsf_f0);
    cp.kcoul = float(kCoulomb);
    cp.energies = energies_;
  }
  auto cluster_phase = [&](int t) {
    tf_[t].assign(x.size(), 0.0);
    const size_t ns = cl_->perm.size();
    std::vector<float>& fc = tfc_[t];
    fc.assign(3 * ns, 0.0f);
    nbnxm::Sums sm;
    nbnxm::kernel(*cl_, cp, t, fc.data(), fc.data() + ns, fc.data() + 2 * ns, sm);
    acc[t][4] += sm.evdw;
    acc[t][5] += sm.ecoul;
    acc[t][6] += sm.w[0] + sm.w[1] + sm.w[2];
    for (int c = 0; c < 6; ++c) acc[t][7 + c] += sm.w[c];
  };

  // Phase 2: bonded terms, 1-4 pairs and the electrostatics of bonded partners, as one index space split over workers.
  const size_t nb = ff_.bonds.size(), na = ff_.angles.size(), nd = ff_.dihedrals.size(), ni = ff_.impropers.size();
  const size_t n14 = ff_.pairs14.size(), nx = coul && !gro ? excl_.size() : 0, nh = ff_.impropers_harmonic.size();
  const size_t o1 = nb, o2 = o1 + na, o3 = o2 + nd, o4 = o3 + ni, o5 = o4 + n14, o6 = o5 + nx, o7 = o6 + nh;
  const size_t c1 = o7 + ff_.bonds2.size(), c2 = c1 + ff_.angles2.size(), c3 = c2 + ff_.dihedrals2.size();
  const size_t c4 = c3 + ff_.impropers2.size();
  const size_t c5 = c4 + ff_.inversions.size(), c6 = c5 + ff_.bonds_x.size(), c7 = c6 + ff_.angles_x.size();
  const size_t c8 = c7 + ff_.urey_bradley.size(), c9 = c8 + ff_.cbt.size();
  const size_t total = c9 + ff_.lj_pairs.size();
  auto bonded_phase = [&](int t, size_t b, size_t en) {
    std::vector<double>& ft = tf_[t];
    auto& A = acc[t];
    auto add = [&](uint32_t i, const Vec3& v) { ft[3 * i] += v[0]; ft[3 * i + 1] += v[1]; ft[3 * i + 2] += v[2]; };
    // virial of force fv at relative position d, scalar and tensor
    auto V = [&](const Vec3& d, const Vec3& fv) {
      A[6] += dot(d, fv);
      A[7] += d[0] * fv[0]; A[8] += d[1] * fv[1]; A[9] += d[2] * fv[2];
      A[10] += 0.5 * (d[0] * fv[1] + d[1] * fv[0]); A[11] += 0.5 * (d[0] * fv[2] + d[2] * fv[0]); A[12] += 0.5 * (d[1] * fv[2] + d[2] * fv[1]);
    };
    // Torsions (proper and improper share the form); φ follows IUPAC, Field.ForcesMatchFiniteDifferences checks the
    // gradients. Torsions of one quadruple i-j-k-l in a row (a Fourier series is one term each) share its geometry: φ from its cosine
    // and sine (no atan2), cos(nφ − δ) and sin(nφ − δ) by the recurrence for integral n, forces once for the summed dE/dφ.
    // Terms [t, lim) of T from t; t is left on the last one taken.
    auto torsion_group = [&](const std::vector<TorsionTerm>& T, size_t& t, size_t lim, double& en) {
      const TorsionTerm& t0 = T[t];
      auto same = [&](size_t u) { return T[u].i == t0.i && T[u].j == t0.j && T[u].k == t0.k && T[u].l == t0.l; };
      const Vec3 b1 = mi(pos(t0.j) - pos(t0.i)), b2 = mi(pos(t0.k) - pos(t0.j)), b3 = mi(pos(t0.l) - pos(t0.k));
      const Vec3 m = cross(b1, b2), nn = cross(b2, b3);
      const double m2 = dot(m, m), n2 = dot(nn, nn), lb2 = norm(b2);
      if (m2 < 1e-12 || n2 < 1e-12 || lb2 < 1e-9) {
        while (t + 1 < lim && same(t + 1)) ++t;
        return;
      }
      const double xc = dot(m, nn), yc = lb2 * dot(b1, nn), rc_ = std::sqrt(xc * xc + yc * yc);
      const double c1 = xc / rc_, s1 = yc / rc_;   // cos φ, sin φ (φ as atan2(y, x), IUPAC)
      double dEdphi = 0;
      for (;;) {
        const TorsionTerm& tt = T[t];
        double cn = 1, sn = 0;
        if (tt.n >= 0 && tt.n <= 12) {
          for (int r = 0; r < tt.n; ++r) {
            const double c = cn * c1 - sn * s1;
            sn = sn * c1 + cn * s1;
            cn = c;
          }
        } else {
          const double phi = std::atan2(yc, xc);
          cn = std::cos(tt.n * phi), sn = std::sin(tt.n * phi);
        }
        double cd = 1, sd = 0;   // δ = 0 or π nearly always
        if (tt.delta != 0) {
          if (std::fabs(tt.delta - kPi) < 1e-12) cd = -1;
          else cd = std::cos(tt.delta), sd = std::sin(tt.delta);
        }
        en += tt.v * (1 + cn * cd + sn * sd);
        dEdphi -= tt.v * tt.n * (sn * cd - cn * sd);
        if (t + 1 < lim && same(t + 1)) ++t;
        else break;
      }
      const Vec3 g1 = m * (-lb2 / m2), g4 = nn * (lb2 / n2);
      const double p = dot(b1, b2) / (lb2 * lb2), qd = dot(b3, b2) / (lb2 * lb2);
      const Vec3 g2 = g4 * qd - g1 * (1 + p);
      const Vec3 g3 = g1 * p - g4 * (1 + qd);
      const Vec3 f1 = g1 * -dEdphi, f2 = g2 * -dEdphi, f3 = g3 * -dEdphi, f4 = g4 * -dEdphi;
      add(t0.i, f1); add(t0.j, f2); add(t0.k, f3); add(t0.l, f4);
      const Vec3 r2 = b1, r3 = b1 + b2, r4 = b1 + b2 + b3;
      V(r2, f2); V(r3, f3); V(r4, f4);
    };
    for (size_t k = b; k < en; ++k) {
      if (k >= c8) {
        if (k < c9) {
          // combined bending–torsion: E = Σ a_n C^n (AB)^((3−n)/2) / (L1^(3/2) L2³ L3^(3/2)), A = |b1×b2|², B = |b2×b3|²,
          // C = (b1×b2)·(b2×b3), L the squared bond lengths
          if (!bonded) continue;
          const auto& t = ff_.cbt[k - c8];
          const Vec3 b1 = mi(pos(t.j) - pos(t.i)), b2 = mi(pos(t.k) - pos(t.j)), b3 = mi(pos(t.l) - pos(t.k));
          const double L1 = dot(b1, b1), L2 = dot(b2, b2), L3 = dot(b3, b3);
          const double d12 = dot(b1, b2), d23 = dot(b2, b3), d13 = dot(b1, b3);
          const double Aa = L1 * L2 - d12 * d12, Bb = L2 * L3 - d23 * d23, Cc = d12 * d23 - d13 * L2;
          const double AB = Aa * Bb;
          if (AB < 1e-24 || L1 < 1e-12 || L2 < 1e-12 || L3 < 1e-12) continue;
          const double D = std::pow(L1, 1.5) * L2 * L2 * L2 * std::pow(L3, 1.5);
          double e = 0, dC = 0, dAB = 0;   // dAB: ∂(Σ)/∂(AB)
          for (int n = 0; n <= 4; ++n) {
            if (t.a[n] == 0) continue;
            const double pw = std::pow(AB, 0.5 * (3 - n)), cn = n ? std::pow(Cc, n) : 1.0;
            e += t.a[n] * cn * pw;
            if (n) dC += t.a[n] * n * std::pow(Cc, n - 1) * pw;
            dAB += t.a[n] * cn * 0.5 * (3 - n) * pw / AB;
          }
          const double E = e / D;
          A[2] += E;
          const double dEdA = dAB * Bb / D, dEdB = dAB * Aa / D, dEdC = dC / D;
          const double dL1 = -1.5 * E / L1, dL2 = -3 * E / L2, dL3 = -1.5 * E / L3;
          // gradients with respect to the bond vectors
          const Vec3 gA1 = (b1 * L2 - b2 * d12) * 2.0, gA2 = (b2 * L1 - b1 * d12) * 2.0;
          const Vec3 gB2 = (b2 * L3 - b3 * d23) * 2.0, gB3 = (b3 * L2 - b2 * d23) * 2.0;
          const Vec3 gC1 = b2 * d23 - b3 * L2, gC3 = b2 * d12 - b1 * L2, gC2 = b1 * d23 + b3 * d12 - b2 * (2 * d13);
          const Vec3 g1 = gA1 * dEdA + gC1 * dEdC + b1 * (2 * dL1);
          const Vec3 g2 = gA2 * dEdA + gB2 * dEdB + gC2 * dEdC + b2 * (2 * dL2);
          const Vec3 g3 = gB3 * dEdB + gC3 * dEdC + b3 * (2 * dL3);
          const Vec3 f1 = g1, f2 = g2 - g1, f3 = g3 - g2, f4 = g3 * -1.0;
          add(t.i, f1); add(t.j, f2); add(t.k, f3); add(t.l, f4);
          V(b1, f2); V(b1 + b2, f3); V(b1 + b2 + b3, f4);
        } else {
          // an explicit LJ pair, no cut-off
          if (!nonb) continue;
          const auto& pr = ff_.lj_pairs[k - c9];
          const Vec3 d = mi(pos(pr.j) - pos(pr.i));
          const double r2 = dot(d, d);
          if (r2 < 1e-12) continue;
          const double sr2 = pr.sigma * pr.sigma / r2, sr6 = sr2 * sr2 * sr2;
          A[4] += 4 * pr.eps * (sr6 * sr6 - sr6);
          double fr = 24 * pr.eps * (2 * sr6 * sr6 - sr6) / r2;
          if (q[pr.i] != 0 && q[pr.j] != 0) {   // and the pair's plain Coulomb, scaled as 1-4 pairs are (GROMACS fudgeQQ)
            const double r = std::sqrt(r2), qq = kCoulomb * q[pr.i] * q[pr.j] * ff_.coul14;
            A[5] += qq / r;
            fr += qq / (r2 * r);
          }
          const Vec3 fj = d * fr;
          add(pr.j, fj);
          add(pr.i, fj * -1.0);
          V(d, fj);
        }
        continue;
      }
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
        size_t t = k - o2;
        torsion_group(ff_.dihedrals, t, std::min(en, o3) - o2, A[2]);
        k = o2 + t;
      } else if (k < o4) {
        size_t t = k - o3;
        torsion_group(ff_.impropers, t, std::min(en, o4) - o3, A[3]);
        k = o3 + t;
      } else if (k < o5) {
        // 1-4 pairs: LJ scaled (with the same cut-off shift), at the nearest image.
        const auto& pr = ff_.pairs14[k - o4];
        const uint32_t i = pr[0], j = pr[1];
        const Vec3 d = mi(pos(j) - pos(i));
        const size_t tp = size_t(ff_.type_index[i]) * nt + ff_.type_index[j];
        double ev = 0, fr = 0;
        const double lj14 = k - o4 < ff_.pairs14_lj.size() ? ff_.pairs14_lj[k - o4] : ff_.lj14;   // its part's scale (merged by group)
        if (!eps14_.empty()) {
          // separate 1-4 Lennard-Jones parameters (CHARMM, GROMOS)
          const double r2 = dot(d, d), qq = s614_[tp] / (r2 * r2 * r2), e14 = eps14_[tp];
          ev = lj14 * 4 * e14 * (qq * qq - qq);
          fr = lj14 * 24 * e14 * (2 * qq * qq - qq) / r2;
          if (ff_.lj_fsw) {   // the force switch's constant below r_in (LAMMPS dihedral charmmfsw adds the same)
            const double ri3 = std::pow(ff_.lj_inner, 3), rc3 = rc2 * std::sqrt(rc2), s6 = s614_[tp];
            ev -= lj14 * 4 * e14 * (s6 * s6 / (ri3 * ri3 * rc3 * rc3) - s6 / (ri3 * rc3));
          } else if (!opt_.tail) {
            const double qc = s614_[tp] / (rc2 * rc2 * rc2);
            ev -= lj14 * 4 * e14 * (qc * qc - qc);
          }
        } else {
          fr = lj(tp, dot(d, d), lj14, ev);
          ev -= lj14 * lj_shift(tp);
        }
        A[4] += ev;
        const Vec3 fj = d * fr;
        add(j, fj);
        add(i, fj * -1.0);
        V(d, fj);
        if (hbon && dot(d, d) < hro2) {   // a hydrogen bond between 1-4 partners (DREIDING counts 1-4 pairs in full)
          auto run = [&](uint32_t D, uint32_t Aa, const Vec3& dDA) {
            const auto it = ff_.hbond.param.find({ff_.type_index[D], ff_.type_index[Aa]});
            if (it == ff_.hbond.param.end()) return;
            for (uint32_t h : ff_.hbond.hyd[D]) {
              const Vec3 d1 = mi(pos(D) - pos(h)), d2 = dDA + d1;
              HbOut o;
              if (!hb_one(it->second, d1, d2, o)) continue;
              A[4] += lj14 * o.e;
              add(D, o.fd * lj14); add(Aa, o.fa * lj14); add(h, o.fh * lj14);
              V(d1, o.fd * lj14); V(d2, o.fa * lj14);
            }
          };
          if (!ff_.hbond.hyd[i].empty() && ff_.hbond.acceptor[j]) run(i, j, d);
          if (!ff_.hbond.hyd[j].empty() && ff_.hbond.acceptor[i]) run(j, i, d * -1.0);
        }
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
            if (c > 0.8660) {   // a wall below 30°, where these periodic forms have a spurious minimum (not in the paper)
              const double w = std::exp(-20 * (th - t.b + 0.25));
              o.e += w;
              dEdth += -20 * w;
            }
          } else if (t.form == 5) {   // restricted bending (Martini 3, GROMACS 10): K (cos θ − cos θ0)² / sin² θ
            const double sn = std::sin(th), s2 = std::max(sn * sn, 1e-12);
            o.e = t.a * dc * dc / s2;
            // d/dθ: −2K dc sin θ / sin² θ − 2K dc² cos θ / sin³ θ
            dEdth = -2 * t.a * dc * sn / s2 - 2 * t.a * dc * dc * c / (s2 * sn);
          } else if (t.form == 4) {   // SDK: K (θ − θ0)²
            o.e = t.a * (th - t.b) * (th - t.b);
            dEdth = 2 * t.a * (th - t.b);
          } else {
            o.e = t.a * dc * dc;
            dEdth = -2 * t.a * dc * std::sin(th);
          }
          o.f[0] = gu * -dEdth; o.f[2] = gv * -dEdth;
          if (t.form == 4) {   // the end atoms' SDK pair, repulsive part only: cut at its minimum, shifted to zero there
            const size_t tp = size_t(ff_.type_index[t.i]) * nt + ff_.type_index[t.k];
            const Vec3 d13 = v - u;   // i → k
            const double r2 = dot(d13, d13);
            if (r2 < rmin2_[tp]) {
              double e13 = 0;
              const double fr = lj(tp, r2, 1.0, e13);
              o.e += e13 - emin_[tp];
              o.f[2] = o.f[2] + d13 * fr;
              o.f[0] = o.f[0] - d13 * fr;
            }
          }
          o.f[1] = (o.f[0] + o.f[2]) * -1.0;
          o.vir(u, o.f[0]); o.vir(v, o.f[2]);
          add(t.i, o.f[0]); add(t.j, o.f[1]); add(t.k, o.f[2]);
          A[1] += o.e;
        } else if (k >= c5) {
          const auto& t = ff_.bonds_x[k - c5];
          const Vec3 d = mi(pos(t.j) - pos(t.i));
          const double r = norm(d);
          double dEdr = 0;
          if (t.form == 3) {   // FENE (LAMMPS bond fene): −½ K R0² ln(1 − (r/R0)²), plus WCA below 2^(1/6) σ when ε > 0
            const double q = r / t.b;
            if (q >= 1) throw FieldError("a FENE bond is stretched past R0 (" + std::to_string(r) + " Å)");
            o.e = -0.5 * t.a * t.b * t.b * std::log(1 - q * q);
            dEdr = t.a * r / (1 - q * q);
            if (t.c > 0 && r < std::pow(2.0, 1.0 / 6) * t.d) {
              const double s6 = std::pow(t.d / r, 6);
              o.e += 4 * t.c * (s6 * s6 - s6) + t.c;
              dEdr += -24 * t.c * (2 * s6 * s6 - s6) / r;
            }
          } else if (t.form == 1) {
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
        if (rf) {
          // reaction field (GROMACS, Verlet): an excluded pair within the cut-off keeps k_rf r² − c_rf, plus what its factor
          // keeps of 1/r (the 1-4 scale)
          const double keep = ex.factor;
          A[5] += qq * (keep / r + krf * r2 - crf);
          fr = qq * (keep / (r2 * r) - 2 * krf);
          const Vec3 fj = d * fr;
          add(j, fj);
          add(i, fj * -1.0);
          V(d, fj);
          continue;
        }
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
        if (r2 < rc2) {
          // the real-space part from the table (as the pairs), less (1 − factor) of 1/r
          double E, Fc;
          ewald_ef(r2, 1 / r, E, Fc);
          if (pme) {
            A[5] += qq * E - (1 - ex.factor) * qq / r;
            fr = qq * (Fc - (1 - ex.factor) / (r2 * r));
          } else {
            A[5] += qq * (E - dsf_e0 + dsf_f0 * (r - rc)) - (1 - ex.factor) * qq / r;
            fr = qq * (Fc - dsf_f0 / r - (1 - ex.factor) / (r2 * r));
          }
        } else if (pme) {
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
  };
  // one run for both phases: each worker lists its pairs, then its share of the bonded terms, into its own buffer
  pool_->run(size_t(nth), [&](int t, size_t, size_t) {
    if (clmode) cluster_phase(t);
    else if (plain && lj96_) pair_phase(t, std::integral_constant<int, 2>{});
    else if (plain) pair_phase(t, std::integral_constant<int, 1>{});
    else pair_phase(t, std::integral_constant<int, 0>{});
    bonded_phase(t, total * size_t(t) / size_t(nth), total * size_t(t + 1) / size_t(nth));
  });

  // Reduce in worker order (deterministic for a given thread count); the reduction itself is split by atoms. Pair
  // forces are in slot order.
  const bool slots = nonb && !clmode && iperm_.size() == n;
  pool_->run(n, [&](int, size_t b, size_t en) {
    for (int t = 0; t < nth; ++t) {
      const std::vector<double>& ft = tf_[t];
      for (size_t k = 3 * b; k < 3 * en; ++k) f[k] += ft[k];
      if (clmode) {
        const std::vector<float>& fc = tfc_[t];
        const size_t ns = cl_->perm.size();
        for (size_t i = b; i < en; ++i) {
          const size_t sl = cl_->iperm[i];
          f[3 * i] += double(fc[sl]); f[3 * i + 1] += double(fc[ns + sl]); f[3 * i + 2] += double(fc[2 * ns + sl]);
        }
        continue;
      }
      if (!slots) continue;
      const std::vector<double>& fs = tfs_[t];
      for (size_t i = b; i < en; ++i) {
        const size_t sl = 3 * size_t(iperm_[i]);
        f[3 * i] += fs[sl]; f[3 * i + 1] += fs[sl + 1]; f[3 * i + 2] += fs[sl + 2];
      }
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
    e.dudl_lj += acc[t][13];
    e.dudl_coul += acc[t][14];
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
    for (uint32_t si = 0; si < n; ++si) {
      const uint32_t i = perm_[si];
      for (size_t k = first_[si]; k < first_[si + 1]; ++k) {
        const uint32_t sj = pj_[k], j = perm_[sj];
        if (!S.atom[i] || !S.atom[j]) continue;
        const double* sv = &qvec_[3 * size_t(pq_[k])];
        const Vec3 d = Vec3{xs[3 * sj], xs[3 * sj + 1], xs[3 * sj + 2]} - Vec3{xs[3 * si], xs[3 * si + 1], xs[3 * si + 2]} + Vec3{sv[0], sv[1], sv[2]};
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
    e.coulomb += pme_reciprocal(x, qeff_, cell, grid, fk, vk, *pool_);
    for (size_t k = 0; k < x.size(); ++k) f[k] += fk[k];
    e.virial += vk[0] + vk[1] + vk[2];
    for (int c = 0; c < 6; ++c) e.w[c] += vk[c];
    double q2 = 0, qt = 0;
    for (double c : qeff_) q2 += c * c, qt += c;
    e.coulomb -= kCoulomb * beta / std::sqrt(kPi) * q2;
    if (std::fabs(qt) > 1e-10) {
      const double eb = -kCoulomb * kPi * qt * qt / (2 * cell.volume() * beta * beta);
      e.coulomb += eb;
      e.virial += 3 * eb;   // E_bg ∝ 1/V ∝ λ⁻³, so −dE/dλ = 3 E_bg
      for (int c = 0; c < 3; ++c) e.w[c] += eb;
    }
  } else if (rf) {
    // GROMACS's reaction-field self term: −½ c_rf Σ q²
    double q2 = 0;
    for (double c : qeff_) q2 += c * c;
    e.coulomb -= 0.5 * crf * q2 * kCoulomb;
  } else if (coul && !gro) {
    double q2 = 0;
    for (double c : qeff_) q2 += c * c;
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
        if (form_[tp] == kPairSdk126) {   // a 12-6 pair in a 9-6 system (another force field's group): its own 12-6 tail
          const double e6 = pa_[tp], sg6 = std::pow(pb_[tp], 6), sg12 = sg6 * sg6;
          et += nn * 2 * kPi * 4 * e6 * (sg12 / (9 * rc9) - sg6 / (3 * rc3));
          wt += nn * 2 * kPi * 4 * e6 * (4 * sg12 / (3 * rc9) - 2 * sg6 / rc3);
          continue;
        }
        if (form_[tp] >= kPairSdk96) continue;   // SDK (LAMMPS refuses a tail for lj/sdk) and lj/gromacs (zero at rc)
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
  const double* xs = listed_positions(x);
  inside_.assign(pj_.size(), 0);
  for (size_t i = 0; i + 1 < first_.size(); ++i)
    for (size_t p = first_[i]; p < first_[i + 1]; ++p) {
      const uint32_t j = pj_[p];
      const double* sv = &qvec_[3 * size_t(pq_[p])];
      double r2 = 0;
      for (int k = 0; k < 3; ++k) {
        const double d = xs[3 * j + k] - xs[3 * i + k] + sv[k];
        r2 += d * d;
      }
      inside_[p] = r2 < rc2;
    }
  frozen_ = true;
}

double torsion_energy(const ForceField& ff, const TorsionRefs& refs, const std::function<Vec3(uint32_t)>& pos, const Cell& cell) {
  const bool per = cell.valid();
  auto mi = [&](const Vec3& d) { return per ? cell.minimum_image(d) : d; };
  double e = 0;
  for (uint32_t k : refs.dihedrals) {   // v [1 + cos(nφ − δ)], φ as the Evaluator's (IUPAC)
    const auto& tt = ff.dihedrals[k];
    const Vec3 b1 = mi(pos(tt.j) - pos(tt.i)), b2 = mi(pos(tt.k) - pos(tt.j)), b3 = mi(pos(tt.l) - pos(tt.k));
    const Vec3 m = cross(b1, b2), nn = cross(b2, b3);
    const double lb2 = norm(b2);
    if (dot(m, m) < 1e-12 || dot(nn, nn) < 1e-12 || lb2 < 1e-9) continue;
    const double phi = std::atan2(lb2 * dot(b1, nn), dot(m, nn));
    e += tt.v * (1 + std::cos(tt.n * phi - tt.delta));
  }
  for (uint32_t k : refs.dihedrals2) {
    const auto& t = ff.dihedrals2[k];
    e += class2_dihedral(t, mi(pos(t.j) - pos(t.i)), mi(pos(t.k) - pos(t.j)), mi(pos(t.l) - pos(t.k))).e;
  }
  for (uint32_t k : refs.cbt) {   // combined bending–torsion, as in Evaluator::compute
    const auto& t = ff.cbt[k];
    const Vec3 b1 = mi(pos(t.j) - pos(t.i)), b2 = mi(pos(t.k) - pos(t.j)), b3 = mi(pos(t.l) - pos(t.k));
    const double L1 = dot(b1, b1), L2 = dot(b2, b2), L3 = dot(b3, b3);
    const double d12 = dot(b1, b2), d23 = dot(b2, b3), d13 = dot(b1, b3);
    const double Aa = L1 * L2 - d12 * d12, Bb = L2 * L3 - d23 * d23, Cc = d12 * d23 - d13 * L2, AB = Aa * Bb;
    if (AB < 1e-24 || L1 < 1e-12 || L2 < 1e-12 || L3 < 1e-12) continue;
    double s = 0;
    for (int n = 0; n <= 4; ++n)
      if (t.a[n] != 0) s += t.a[n] * (n ? std::pow(Cc, n) : 1.0) * std::pow(AB, 0.5 * (3 - n));
    e += s / (std::pow(L1, 1.5) * L2 * L2 * L2 * std::pow(L3, 1.5));
  }
  return e;
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
