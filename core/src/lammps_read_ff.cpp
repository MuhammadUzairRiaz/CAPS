// The force field of a LAMMPS data file CAPS wrote (write_lammps_data_ff), rebuilt from the file: pair and bonded
// coefficients in the styles the writer uses, charges and masses from Atoms and Masses, and from the "# caps-ff" lines
// under the title the force field's name and its 1-2 / 1-3 / 1-4 weights (special_bonds, which LAMMPS keeps in the input).
// Anything else (another program's file, a style not read back, what the writer says the file cannot hold) leaves the
// reason in System::forcefield_error, so runs ask for --ff instead of falling back to GAFF or UFF.
#include "lammps_read_ff.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>

#include "caps/field.hpp"

namespace caps {

namespace {

constexpr double kD2R = 0.017453292519943295;
constexpr double kPi = 3.141592653589793;

struct Unreadable : std::runtime_error {
  using std::runtime_error::runtime_error;
};

void need(bool ok, const std::string& why) {
  if (!ok) throw Unreadable(why);
}

bool is_number(const std::string& t) {
  if (t.empty()) return false;
  char* end = nullptr;
  std::strtod(t.c_str(), &end);
  return end && *end == '\0';
}

double number(const std::string& t, const std::string& where) {
  need(is_number(t), where + ": '" + t + "' is not a number");
  return std::strtod(t.c_str(), nullptr);
}

// One coefficient line: its style (the hybrid word, or the section's style) and numbers
struct Coef {
  std::string style;
  std::vector<double> v;
};

// "N Coeffs" rows by type id. Hybrid sections name each row's style; a class II cross section's "skip" rows are left out.
std::map<int, Coef> coefficients(const LammpsFileTerms& d, const std::string& section, const std::string& owner_style = "") {
  std::map<int, Coef> out;
  auto it = d.rows.find(section);
  if (it == d.rows.end()) return out;
  auto st = d.style.find(section);
  std::string style = st == d.style.end() ? owner_style : st->second;
  if (style.empty()) style = owner_style;
  for (const auto& r : it->second) {
    need(!r.empty() && is_number(r[0]), section + ": a row without its type number");
    const int t = int(number(r[0], section));
    Coef c;
    size_t k = 1;
    if (r.size() > 1 && !is_number(r[1])) c.style = r[k++];
    else c.style = style;
    if (c.style == "skip") continue;
    for (; k < r.size(); ++k) {
      if (!is_number(r[k])) {   // a word among the numbers: an improper fourier's "all" flag stays a number
        need(false, section + " type " + std::to_string(t) + ": '" + r[k] + "' (CAPS reads back the coefficients it writes, all numbers)");
      }
      c.v.push_back(number(r[k], section));
    }
    need(out.emplace(t, std::move(c)).second, section + ": type " + std::to_string(t) + " given twice");
  }
  return out;
}

const Coef& coef_of(const std::map<int, Coef>& m, int t, const std::string& what) {
  auto it = m.find(t);
  need(it != m.end(), what + " type " + std::to_string(t) + " has no coefficients in the file");
  return it->second;
}

void count(const Coef& c, size_t n, const std::string& what) {
  need(c.v.size() == n, what + " " + c.style + ": " + std::to_string(c.v.size()) + " coefficients, CAPS writes " + std::to_string(n));
}

bool close(double a, double b, double rel = 1e-9) { return std::fabs(a - b) <= rel * std::max({1.0, std::fabs(a), std::fabs(b)}); }

// The LJ / Buckingham / Morse base of a pair style ("lj/cut/coul/long" → "lj/cut"); "" for a Coulomb-only sub-style
std::string pair_base(const std::string& st) {
  if (st.rfind("coul/", 0) == 0 || st.rfind("tip4p/", 0) == 0) return "";
  const size_t p = st.find("/coul/");
  return p == std::string::npos ? st : st.substr(0, p);
}

std::string first_word(const std::string& s) {
  std::istringstream is(s);
  std::string w;
  is >> w;
  return w;
}

ForceField build(const System& s, const LammpsFileTerms& d, const std::unordered_map<int64_t, uint32_t>& index, const double lj14,
                 const double coul14, const bool keep13) {
  ForceField F;
  const size_t n = s.atoms.size();
  // types: 1..N from Masses, named by their comments
  int nt = 0;
  for (const auto& t : s.types) nt = std::max(nt, t.type);
  for (const auto& a : s.atoms) need(a.type >= 1 && a.type <= nt, "atom " + std::to_string(a.id) + " has type " + std::to_string(a.type) + ", which Masses does not list");
  need(s.has_charges, "no charges (atom_style full has them)");
  F.type_names.assign(size_t(nt), "");
  std::vector<double> tmass(size_t(nt), 0.0);
  for (const auto& t : s.types) {
    F.type_names[size_t(t.type - 1)] = first_word(t.label).empty() ? "type" + std::to_string(t.type) : first_word(t.label);
    tmass[size_t(t.type - 1)] = t.mass;
  }
  for (int t = 0; t < nt; ++t)
    if (F.type_names[size_t(t)].empty()) F.type_names[size_t(t)] = "type" + std::to_string(t + 1);
  for (const auto& a : s.atoms) {
    const size_t t = size_t(a.type - 1);
    F.type_index.push_back(int(t));
    F.atom_type.push_back(F.type_names[t]);
    F.why.push_back("type " + std::to_string(t + 1) + " of the data file");
    F.charge.push_back(a.charge);
    F.mass.push_back(tmass[t]);
  }
  F.lj14 = lj14;
  F.coul14 = coul14;
  F.keep13 = keep13;

  // pairs: every i-j pair written out (PairIJ Coeffs), so nothing is left to LAMMPS's mixing
  need(!d.rows.count("Pair Coeffs"), "Pair Coeffs (pairs mixed by LAMMPS; CAPS writes every pair in PairIJ Coeffs)");
  auto pij = d.rows.find("PairIJ Coeffs");
  need(pij != d.rows.end() && !pij->second.empty(), "no pair coefficients (PairIJ Coeffs)");
  const std::string pstyle = d.style.count("PairIJ Coeffs") ? first_word(d.style.at("PairIJ Coeffs")) : "";
  need(!pstyle.empty(), "PairIJ Coeffs names no pair style");
  const bool phybrid = pstyle.rfind("hybrid", 0) == 0;
  struct P { std::string base; std::vector<double> v; };
  std::map<std::pair<int, int>, P> pairs;
  bool class2 = false;
  for (const auto& r : pij->second) {
    need(r.size() >= 2 && is_number(r[0]) && is_number(r[1]), "PairIJ Coeffs: a row without its two types");
    int a = int(number(r[0], "PairIJ")) - 1, b = int(number(r[1], "PairIJ")) - 1;
    need(a >= 0 && b >= 0 && a < nt && b < nt, "PairIJ Coeffs: a type outside 1.." + std::to_string(nt));
    if (a > b) std::swap(a, b);
    size_t k = 2;
    std::string st = pstyle;
    if (phybrid) {
      need(r.size() > 2 && !is_number(r[2]), "PairIJ Coeffs (hybrid): a row without its sub-style");
      st = r[k++];
    }
    const std::string base = pair_base(st);
    if (base.empty()) continue;   // the Coulomb sub-style of a hybrid
    P p{base, {}};
    for (; k < r.size(); ++k) p.v.push_back(number(r[k], "PairIJ Coeffs"));
    const std::string where = "pair " + std::to_string(a + 1) + "-" + std::to_string(b + 1) + " (" + st + ")";
    if (base == "lj/cut" || base == "lj/class2") {
      need(p.v.size() == 2, where + ": ε σ only (a pair's own cut-off is not read back)");
      class2 = class2 || base == "lj/class2";
    } else if (base == "buck" || base == "morse") {
      need(p.v.size() == 3 || p.v.size() == 4, where + ": 3 coefficients and an optional cut-off");
    } else if (base == "born") {   // a shell model's Buckingham written as Born–Mayer–Huggins (σ = 0, D = 0)
      need(p.v.size() == 5 && p.v[2] == 0 && p.v[4] == 0, where + ": only Buckingham as Born–Mayer–Huggins (σ = 0, D = 0) is read back");
      p = {"buck", {p.v[0], p.v[1], p.v[3]}};
    } else {
      need(false, "pair style " + st + " (CAPS reads back lj/cut, lj/class2, buck, morse and their Coulomb forms)");
    }
    need(pairs.emplace(std::make_pair(a, b), std::move(p)).second, where + ": two van der Waals terms for one pair");
  }
  F.pair_form = class2 ? "lj9-6" : "lj12-6";
  const std::string lj_base = class2 ? "lj/class2" : "lj/cut";
  // per type: its own pair when that is the base Lennard-Jones (else zero); unlike pairs by the rule they follow, the rest explicit
  F.lj.assign(size_t(nt), PairType{0, 0});
  for (int t = 0; t < nt; ++t)
    if (auto it = pairs.find({t, t}); it != pairs.end() && it->second.base == lj_base) F.lj[size_t(t)] = {it->second.v[0], it->second.v[1]};
  std::map<std::pair<int, int>, PairType> lj_pairs;   // every base Lennard-Jones pair (absent pairs: zero)
  for (int a = 0; a < nt; ++a)
    for (int b = a; b < nt; ++b) {
      auto it = pairs.find({a, b});
      if (it == pairs.end()) {
        lj_pairs[{a, b}] = {0.0, 0.0};
      } else if (it->second.base == lj_base) {
        lj_pairs[{a, b}] = {it->second.v[0], it->second.v[1]};
      } else {
        const auto& v = it->second.v;
        PairFunc pf{0, 0, 0, 0};
        if (it->second.base == "lj/cut") pf = {kPairSdk126, v[0], v[1], 0};   // a 12-6 pair in a 9-6 system
        else pf = {it->second.base == "buck" ? 1 : 2, v[0], v[1], v[2], v.size() > 3 ? v[3] : 0.0};
        F.pair_func[{a, b}] = pf;
      }
    }
  auto matches = [&](const PairType& want, const PairType& got) {
    return close(want.eps, got.eps, 1e-10) && (want.eps == 0 || close(want.sigma, got.sigma, 1e-10));
  };
  size_t best_miss = size_t(-1);
  std::string best = "geometric";
  for (const char* rule : {"geometric", "arithmetic", "sixthpower"}) {
    F.mixing = rule;
    size_t miss = 0;
    for (const auto& [ab, pt] : lj_pairs)
      if (ab.first != ab.second && !matches(pt, mixed_pair(F, ab.first, ab.second))) ++miss;
    if (miss < best_miss) best_miss = miss, best = rule;
  }
  F.mixing = best;
  for (const auto& [ab, pt] : lj_pairs) {
    if (ab.first == ab.second) {
      if (!matches(pt, F.lj[size_t(ab.first)])) F.pair_override[ab] = pt;   // a type whose own pair is not Lennard-Jones
      continue;
    }
    if (!matches(pt, mixed_pair(F, ab.first, ab.second))) F.pair_override[ab] = pt;
  }

  auto atom = [&](int64_t id, const std::string& what) {
    auto it = index.find(id);
    need(it != index.end(), what + " names atom " + std::to_string(id) + ", which Atoms does not hold");
    return it->second;
  };
  // bonds (a bond with style zero, or none, is a bond with no term: it still excludes its pairs)
  const auto bc = coefficients(d, "Bond Coeffs");
  need(d.bond_type.size() == s.bonds.size(), "Bonds: rows that name atoms the file does not hold");
  for (size_t q = 0; q < s.bonds.size(); ++q) {
    const uint32_t i = s.bonds[q].i, j = s.bonds[q].j;
    if (bc.empty()) need(false, "Bonds without Bond Coeffs");
    const Coef& c = coef_of(bc, d.bond_type[q], "bond");
    if (c.style == "zero") continue;
    if (c.style == "harmonic") count(c, 2, "bond"), F.bonds.push_back({i, j, c.v[0], c.v[1]});
    else if (c.style == "class2") count(c, 4, "bond"), F.bonds2.push_back({i, j, c.v[0], c.v[1], c.v[2], c.v[3]});
    else if (c.style == "morse") count(c, 3, "bond"), F.bonds_x.push_back({i, j, 1, c.v[0], c.v[1], c.v[2]});
    else if (c.style == "gromos") count(c, 2, "bond"), F.bonds_x.push_back({i, j, 2, c.v[0], c.v[1], 0});
    else if (c.style == "fene") count(c, 4, "bond"), F.bonds_x.push_back({i, j, 3, c.v[0], c.v[1], c.v[2], c.v[3]});
    else need(false, "bond style " + c.style);
  }
  // angles (class II: BondBond and BondAngle beside the angle)
  const auto ac = coefficients(d, "Angle Coeffs"), abb = coefficients(d, "BondBond Coeffs", "class2"), aba = coefficients(d, "BondAngle Coeffs", "class2");
  for (const auto& t : d.angles) {
    need(t.ids.size() == 3, "Angles: three atoms each");
    const uint32_t i = atom(t.ids[0], "an angle"), j = atom(t.ids[1], "an angle"), k = atom(t.ids[2], "an angle");
    const Coef& c = coef_of(ac, t.type, "angle");
    const auto& v = c.v;
    if (c.style == "zero") continue;
    if (c.style == "harmonic") {
      count(c, 2, "angle");
      F.angles.push_back({i, j, k, v[0], v[1] * kD2R});
    } else if (c.style == "charmm") {   // Urey–Bradley joined to its angle
      count(c, 4, "angle");
      F.angles.push_back({i, j, k, v[0], v[1] * kD2R});
      if (v[2] != 0) F.urey_bradley.push_back({i, k, v[2], v[3]});
    } else if (c.style == "class2") {
      count(c, 4, "angle");
      const Coef& b1 = coef_of(abb, t.type, "BondBond"), &b2 = coef_of(aba, t.type, "BondAngle");
      count(b1, 3, "BondBond"), count(b2, 4, "BondAngle");
      F.angles2.push_back({i, j, k, v[0] * kD2R, v[1], v[2], v[3], b1.v[0], b1.v[1], b1.v[2], b2.v[0], b2.v[1], b2.v[2], b2.v[3]});
    } else if (c.style == "cosine/squared") {
      count(c, 2, "angle");
      F.angles_x.push_back({i, j, k, 1, v[0], v[1] * kD2R});
    } else if (c.style == "cosine/squared/restricted") {
      count(c, 2, "angle");
      F.angles_x.push_back({i, j, k, 5, v[0], v[1] * kD2R});
    } else if (c.style == "cosine") {
      count(c, 1, "angle");
      F.angles_x.push_back({i, j, k, 2, v[0], 0});
    } else if (c.style == "fourier") {   // UFF's general bend: θ0 from C1 = −4 C2 cos θ0, the C's as UFF sets them
      count(c, 4, "angle");
      need(v[3] > 0, "angle fourier with C2 ≤ 0 is not UFF's form");
      const double c0 = std::clamp(-v[2] / (4 * v[3]), -1.0, 1.0), s0 = 1 - c0 * c0;
      need(close(v[3], 1 / (4 * std::max(s0, 1e-8)), 1e-7) && close(v[1], v[3] * (2 * c0 * c0 + 1), 1e-7),
           "angle fourier type " + std::to_string(t.type) + " is not UFF's form (C0, C1, C2 from θ0)");
      F.angles_x.push_back({i, j, k, 3, v[0], std::acos(c0)});
    } else if (c.style == "cosine/periodic") {   // UFF's linear, trigonal, square and octahedral bends
      count(c, 3, "angle");
      const int m = int(std::lround(v[2])), B = int(std::lround(v[1]));
      need(m >= 1 && m <= 6 && B == (m == 1 ? 1 : (m % 2 ? -1 : 1)), "angle cosine/periodic type " + std::to_string(t.type) + " is not UFF's form");
      F.angles_x.push_back({i, j, k, 10 + m, 2 * v[0], 0});
    } else {
      need(false, "angle style " + c.style);
    }
  }
  // dihedrals (class II: five cross sections)
  const auto dc = coefficients(d, "Dihedral Coeffs");
  const auto mbt = coefficients(d, "MiddleBondTorsion Coeffs", "class2"), ebt = coefficients(d, "EndBondTorsion Coeffs", "class2"),
             at = coefficients(d, "AngleTorsion Coeffs", "class2"), aat = coefficients(d, "AngleAngleTorsion Coeffs", "class2"),
             bb13 = coefficients(d, "BondBond13 Coeffs", "class2");
  for (const auto& t : d.dihedrals) {
    need(t.ids.size() == 4, "Dihedrals: four atoms each");
    uint32_t q[4];
    for (int x = 0; x < 4; ++x) q[x] = atom(t.ids[size_t(x)], "a dihedral");
    const Coef& c = coef_of(dc, t.type, "dihedral");
    const auto& v = c.v;
    auto torsion = [&](double K, int m, double delta) { F.dihedrals.push_back({q[0], q[1], q[2], q[3], K, m, delta}); };
    if (c.style == "zero") continue;
    if (c.style == "fourier") {
      need(!v.empty() && v.size() == 1 + 3 * size_t(std::lround(v[0])), "dihedral fourier type " + std::to_string(t.type) + ": m then m (K n d)");
      for (size_t x = 0; x < size_t(std::lround(v[0])); ++x) torsion(v[1 + 3 * x], int(std::lround(v[2 + 3 * x])), v[3 + 3 * x] * kD2R);
    } else if (c.style == "opls") {   // ½K1(1+cos φ) + ½K2(1−cos 2φ) + ½K3(1+cos 3φ) + ½K4(1−cos 4φ)
      count(c, 4, "dihedral");
      for (int m = 1; m <= 4; ++m)
        if (v[size_t(m - 1)] != 0) torsion(v[size_t(m - 1)] / 2, m, m % 2 ? 0.0 : kPi);
    } else if (c.style == "harmonic") {   // K [1 + d cos(nφ)]
      count(c, 3, "dihedral");
      need(v[1] == 1 || v[1] == -1, "dihedral harmonic with d other than ±1");
      torsion(v[0], int(std::lround(v[2])), v[1] > 0 ? 0.0 : kPi);
    } else if (c.style == "charmm") {
      count(c, 4, "dihedral");
      need(v[3] == 0, "dihedral charmm with a 1-4 weight (CHARMM's 1-4 pairs through the torsions)");
      torsion(v[0], int(std::lround(v[1])), v[2] * kD2R);
    } else if (c.style == "class2") {
      count(c, 6, "dihedral");
      const Coef &m1 = coef_of(mbt, t.type, "MiddleBondTorsion"), &e1 = coef_of(ebt, t.type, "EndBondTorsion"), &a1 = coef_of(at, t.type, "AngleTorsion"),
                 &aa = coef_of(aat, t.type, "AngleAngleTorsion"), &b13 = coef_of(bb13, t.type, "BondBond13");
      count(m1, 4, "MiddleBondTorsion"), count(e1, 8, "EndBondTorsion"), count(a1, 8, "AngleTorsion"), count(aa, 3, "AngleAngleTorsion"), count(b13, 3, "BondBond13");
      Class2Dihedral x{};
      x.i = q[0], x.j = q[1], x.k = q[2], x.l = q[3];
      x.k1 = v[0], x.phi1 = v[1] * kD2R, x.k2 = v[2], x.phi2 = v[3] * kD2R, x.k3 = v[4], x.phi3 = v[5] * kD2R;
      for (int z = 0; z < 3; ++z) {
        x.mbt[z] = m1.v[size_t(z)];
        x.ebt_b[z] = e1.v[size_t(z)], x.ebt_c[z] = e1.v[size_t(3 + z)];
        x.at_d[z] = a1.v[size_t(z)], x.at_e[z] = a1.v[size_t(3 + z)];
      }
      x.mbt_r2 = m1.v[3];
      x.ebt_r1 = e1.v[6], x.ebt_r3 = e1.v[7];
      x.at_theta1 = a1.v[6] * kD2R, x.at_theta2 = a1.v[7] * kD2R;
      x.aat_m = aa.v[0], x.aat_theta1 = aa.v[1] * kD2R, x.aat_theta2 = aa.v[2] * kD2R;
      x.bb13_n = b13.v[0], x.bb13_r1 = b13.v[1], x.bb13_r3 = b13.v[2];
      F.dihedrals2.push_back(x);
    } else {
      need(false, "dihedral style " + c.style);
    }
  }
  // impropers (class II: AngleAngle beside them); DREIDING's umbrella as CAPS writes it, three permutations at K/3
  const auto ic = coefficients(d, "Improper Coeffs"), iaa = coefficients(d, "AngleAngle Coeffs", "class2");
  std::map<std::pair<uint32_t, std::array<uint32_t, 3>>, std::vector<std::pair<double, std::array<uint32_t, 3>>>> umbrella;
  std::vector<std::pair<uint32_t, std::array<uint32_t, 3>>> umbrella_order;
  for (const auto& t : d.impropers) {
    need(t.ids.size() == 4, "Impropers: four atoms each");
    uint32_t q[4];
    for (int x = 0; x < 4; ++x) q[x] = atom(t.ids[size_t(x)], "an improper");
    const Coef& c = coef_of(ic, t.type, "improper");
    const auto& v = c.v;
    if (c.style == "zero") continue;
    if (c.style == "cvff") {   // K [1 + d cos(nχ)]
      count(c, 3, "improper");
      need(v[1] == 1 || v[1] == -1, "improper cvff with d other than ±1");
      F.impropers.push_back({q[0], q[1], q[2], q[3], v[0], int(std::lround(v[2])), v[1] > 0 ? 0.0 : kPi});
    } else if (c.style == "harmonic") {
      count(c, 2, "improper");
      F.impropers_harmonic.push_back({q[0], q[1], q[2], q[3], v[0], v[1] * kD2R});
    } else if (c.style == "class2") {
      count(c, 2, "improper");
      const Coef& a = coef_of(iaa, t.type, "AngleAngle");
      count(a, 6, "AngleAngle");
      F.impropers2.push_back({q[0], q[1], q[2], q[3], v[0], v[1] * kD2R, a.v[0], a.v[1], a.v[2], a.v[3] * kD2R, a.v[4] * kD2R, a.v[5] * kD2R});
    } else if (c.style == "inversion/harmonic") {
      count(c, 2, "improper");
      F.inversions.push_back({q[0], q[1], q[2], q[3], v[0], v[1] * kD2R, 0});
    } else if (c.style == "umbrella") {
      count(c, 2, "improper");
      need(v[1] == 0, "improper umbrella with ω0 ≠ 0 (CAPS writes DREIDING's planar form only)");
      std::array<uint32_t, 3> o{q[1], q[2], q[3]};
      std::array<uint32_t, 3> key = o;
      std::sort(key.begin(), key.end());
      auto& g = umbrella[{q[0], key}];
      if (g.empty()) umbrella_order.push_back({q[0], key});
      g.push_back({v[0], o});
    } else if (c.style == "fourier") {   // centre first, all = 1: the three permutations summed
      count(c, 5, "improper");
      need(v[4] == 1, "improper fourier with all = 0 (CAPS writes all three permutations)");
      if (v[1] == 1 && v[2] == -1 && v[3] == 0) {   // K (1 − cos ω): DREIDING's planar form, UFF's own style
        F.inversions.push_back({q[0], q[1], q[2], q[3], 3 * v[0], 0.0, 1});
      } else {
        need(v[3] == 1, "improper fourier type " + std::to_string(t.type) + " is not UFF's form");
        const double w0 = std::acos(std::clamp(-v[2] / 4, -1.0, 1.0));
        need(close(v[1], -(v[2] * std::cos(w0) + std::cos(2 * w0)), 1e-7), "improper fourier type " + std::to_string(t.type) + " is not UFF's form");
        F.inversions.push_back({q[0], q[1], q[2], q[3], 3 * v[0], w0, 2});
      }
    } else {
      need(false, "improper style " + c.style);
    }
  }
  for (const auto& key : umbrella_order) {
    const auto& g = umbrella[key];
    need(g.size() == 3 && close(g[0].first, g[1].first) && close(g[0].first, g[2].first),
         "umbrella impropers about atom " + std::to_string(key.first + 1) + " are not the three permutations CAPS writes");
    // written centre, o[p+1], o[p+2], o[p] for p = 0: the first line's last atom is the inversion's a
    const auto& o = g[0].second;
    F.inversions.push_back({key.first, o[2], o[0], o[1], 3 * g[0].first, 0.0, 1});
  }

  // exclusions and 1-4 pairs as LAMMPS's special_bonds makes them from the bonds: 1-2, 1-3 (unless kept in full), and the
  // 1-4 pairs that are neither
  std::vector<std::vector<uint32_t>> nb(n);
  std::set<std::pair<uint32_t, uint32_t>> ex12, ex13, p14;
  for (const auto& b : s.bonds) {
    if (b.i == b.j) continue;
    nb[b.i].push_back(b.j), nb[b.j].push_back(b.i);
    ex12.insert({std::min(b.i, b.j), std::max(b.i, b.j)});
  }
  for (uint32_t j = 0; j < n; ++j)
    for (uint32_t i : nb[j])
      for (uint32_t k : nb[j])
        if (i < k && !ex12.count({i, k})) ex13.insert({i, k});
  for (const auto& [j, k] : ex12)
    for (uint32_t i : nb[j])
      for (uint32_t l : nb[k]) {
        if (i == k || l == j || i == l) continue;
        const auto key = std::make_pair(std::min(i, l), std::max(i, l));
        if (!ex12.count(key) && !ex13.count(key)) p14.insert(key);
      }
  if (keep13 && lj14 == 1 && coul14 == 1) p14.clear();   // 1-3 and 1-4 in full: ordinary pairs
  for (const auto& p : p14) F.pairs14.push_back({p.first, p.second});
  F.excluded.assign(n, {});
  auto add_ex = [&](uint32_t a, uint32_t b) { F.excluded[a].push_back(b), F.excluded[b].push_back(a); };
  for (const auto& p : ex12) add_ex(p.first, p.second);
  if (!keep13)
    for (const auto& p : ex13) add_ex(p.first, p.second);
  for (const auto& p : p14) add_ex(p.first, p.second);
  for (auto& e : F.excluded) {
    std::sort(e.begin(), e.end());
    e.erase(std::unique(e.begin(), e.end()), e.end());
  }
  return F;
}

}  // namespace

void rebuild_lammps_forcefield(System& s, const LammpsFileTerms& d, const std::unordered_map<int64_t, uint32_t>& index, const std::string& path) {
  static const char* kCoeffs[] = {"Pair Coeffs", "PairIJ Coeffs", "Bond Coeffs", "Angle Coeffs", "Dihedral Coeffs", "Improper Coeffs"};
  if (std::none_of(std::begin(kCoeffs), std::end(kCoeffs), [&](const char* k) { return d.rows.count(k) && !d.rows.at(k).empty(); })) return;
  const std::string file = path.substr(path.find_last_of("/\\") == std::string::npos ? 0 : path.find_last_of("/\\") + 1);
  s.forcefield_atoms = s.atoms.size();
  auto refuse = [&](const std::string& why) {
    s.forcefield_error = file + " holds force-field coefficients CAPS cannot use as they are: " + why + ". Give --ff FF.json for the run";
  };
  std::string name, left_out;
  std::vector<double> special;
  for (const auto& line : d.caps_ff) {
    std::istringstream is(line);
    std::string key;
    is >> key;
    std::string rest;
    std::getline(is, rest);
    rest.erase(0, rest.find_first_not_of(" \t"));
    while (!rest.empty() && std::isspace(static_cast<unsigned char>(rest.back()))) rest.pop_back();
    if (key == "name") name = rest;
    else if (key == "not-in-file") left_out = rest;
    else if (key == "special") {   // lj w12 w13 w14 coul w12 w13 w14
      std::istringstream ss(rest);
      std::string w;
      std::vector<std::string> t;
      while (ss >> w) t.push_back(w);
      if (t.size() == 8 && t[0] == "lj" && t[4] == "coul")
        for (int k : {1, 2, 3, 5, 6, 7}) special.push_back(is_number(t[size_t(k)]) ? std::strtod(t[size_t(k)].c_str(), nullptr) : NAN);
    }
  }
  if (special.size() != 6) {
    refuse("it does not say its 1-2, 1-3 and 1-4 weights (special_bonds; CAPS writes them under the title, so another program "
           "or an earlier CAPS wrote this file)");
    return;
  }
  if (!left_out.empty()) {
    refuse("its force field has " + left_out + ", which a data file does not hold");
    return;
  }
  if (special[0] != 0 || special[3] != 0 || special[1] != special[4] || (special[1] != 0 && special[1] != 1) || std::isnan(special[2]) || std::isnan(special[5])) {
    refuse("1-2 or 1-3 weights other than 0, or 0 / 1 for both");
    return;
  }
  try {
    ForceField F = build(s, d, index, special[2], special[5], special[1] == 1);
    F.name = name.empty() ? "the force field of " + file : name;
    F.notes.push_back("force field: " + F.name + ", rebuilt from " + file + " (its coefficients, charges and 1-4 weights)");
    s.forcefield = std::make_shared<const ForceField>(std::move(F));
    s.notes.push_back("force field " + s.forcefield->name + " read from the file's coefficients");
  } catch (const Unreadable& e) {
    refuse(e.what());
  }
}

}  // namespace caps
