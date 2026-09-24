// DL_FIELD force-field libraries (C. W. Yong, STFC Daresbury Laboratory): .par parameters, .sf atom types,
// .bci bond charge increments. Each family's columns are interpreted the way DL_FIELD writes them for LAMMPS
// (bench/ff/compare_dlfield.py checks CAPS against DL_FIELD's own LAMMPS output, term by term).
#include <tuple>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <map>
#include <set>
#include <sstream>

#include "caps/elements.hpp"
#include "caps/ffdef.hpp"

namespace caps {

namespace {

bool is_number(const std::string& w) {
  if (w.empty()) return false;
  char* end = nullptr;
  std::strtod(w.c_str(), &end);
  return end && *end == 0;
}

// Comments: '!' anywhere (CHARMM style, even glued to a number: "180.0!"), '#' at the start of a line or after
// white space ('#' also appears inside DL_F notation names such as "Li+#4_binary_oxide1").
std::string strip_comments(std::string line) {
  if (auto h = line.find('!'); h != std::string::npos) line = line.substr(0, h);
  for (size_t k = 0; k < line.size(); ++k)
    if (line[k] == '#' && (k == 0 || line[k - 1] == ' ' || line[k - 1] == '\t')) { line = line.substr(0, k); break; }
  return line;
}

std::string upper(std::string s) {
  for (auto& c : s) c = char(std::toupper(static_cast<unsigned char>(c)));
  return s;
}

// One data line of a section: atom keys, the numbers after them, and the remark.
struct Row {
  std::vector<std::string> keys;
  std::string kind;              // a word between the keys and the numbers (inorganic "buckingham", GROMOS none)
  std::vector<double> v;
  std::string remark;
  int line = 0;
};

struct ParFile {
  std::string potential, unit;
  std::map<std::string, std::vector<Row>> sec;                 // BOND, ANGLE, DIHEDRAL, IMPROPER, INVERSION, VDW, VDW_FIX, ...
  std::map<std::string, std::vector<std::string>> header;      // column names after the section keyword
  std::vector<std::pair<std::string, std::vector<std::string>>> equivalence;   // key > vdw_x bond_y ...
};

const std::map<std::string, int> kKeys = {{"BOND", 2}, {"ANGLE", 3}, {"DIHEDRAL", 4}, {"IMPROPER", 4}, {"INVERSION", 4},
                                          {"VDW", 1},  {"VDW_FIX", 2}, {"SHELL", 2}, {"UREY_BRADLEY", 3}};

ParFile read_par(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw FFError("cannot open " + path);
  ParFile pf;
  std::string line, cur;
  int lineno = 0;
  while (std::getline(in, line)) {
    ++lineno;
    line = strip_comments(line);
    std::istringstream ls(line);
    std::vector<std::string> t;
    for (std::string w; ls >> w;) t.push_back(w);
    if (t.empty()) continue;
    const std::string head = upper(t[0]);
    if (head == "UNIT" && t.size() > 1) { pf.unit = t[1]; continue; }
    if (head == "POTENTIAL" && t.size() > 1) { pf.potential = t[1]; continue; }
    if (head == "END") { cur.clear(); continue; }
    if (head == "TERMINATE") break;
    if (cur.empty()) {
      if (kKeys.count(head) || head == "EQUIVALENCE") {
        cur = head;
        pf.sec[cur];
        pf.header[cur] = std::vector<std::string>(t.begin() + 1, t.end());
      }
      continue;
    }
    if (cur == "EQUIVALENCE") {
      // key > vdw_a bond_b ...
      if (t.size() >= 3 && t[1] == ">") pf.equivalence.push_back({t[0], std::vector<std::string>(t.begin() + 2, t.end())});
      continue;
    }
    const int nk = kKeys.at(cur);
    if (int(t.size()) < nk + 1) continue;
    Row r;
    r.line = lineno;
    size_t k = 0;
    for (; k < size_t(nk); ++k) r.keys.push_back(t[k]);
    // VDW rows of inorganic / MISC files carry a second key before the potential kind
    if (cur == "VDW" && k < t.size() && !is_number(t[k]) && k + 1 < t.size() && !is_number(t[k + 1])) r.keys.push_back(t[k++]);
    if (k < t.size() && !is_number(t[k])) r.kind = t[k++];
    for (; k < t.size() && is_number(t[k]); ++k) r.v.push_back(std::stod(t[k]));
    for (; k < t.size(); ++k) r.remark += (r.remark.empty() ? "" : " ") + t[k];
    if (r.v.empty()) continue;
    pf.sec[cur].push_back(r);
  }
  return pf;
}

double energy_factor(const std::string& unit) {
  const std::string u = upper(unit);
  if (u == "KCAL/MOL" || u.empty()) return 1.0;
  if (u == "KJ/MOL") return 1.0 / 4.184;
  if (u == "EV") return 23.060548;
  if (u == "K") return 0.0019872036;
  throw FFError("unknown DL_FIELD energy unit '" + unit + "'");
}

// "X" is DL_FIELD's wildcard; every other key is a literal name (PCFF has types such as "h*")
std::string key(const std::string& k) { return k == "X" ? "*" : glob_escape(k); }

std::string rule_name(const std::vector<std::string>& keys) {
  std::string n;
  for (const auto& k : keys) n += (n.empty() ? "" : "-") + k;
  return n;
}

std::string ref_of(const Row& r) {
  return r.remark;
}

// Rules are matched last-first: wildcard rules go first (most wildcards first) so specific rules win, and within
// equal specificity the first entry in the file wins (it is placed last).
void order_rules(std::vector<FFRule>& rules) {
  auto wild = [](const FFRule& r) { return int(std::count(r.match.begin(), r.match.end(), std::string("*"))); };
  std::vector<std::pair<int, size_t>> idx;
  for (size_t i = 0; i < rules.size(); ++i) idx.push_back({wild(rules[i]), i});
  std::stable_sort(idx.begin(), idx.end(), [](const auto& a, const auto& b) {
    if (a.first != b.first) return a.first > b.first;
    return a.second > b.second;
  });
  std::vector<FFRule> out;
  for (const auto& [w, i] : idx) out.push_back(rules[i]);
  rules = std::move(out);
}

// Consecutive DIHEDRAL rows for the same keys form one torsion: DL_FIELD marks "more terms follow" with a
// negative multiplicity.
// DL_FIELD takes the rows after a "more terms follow" marker as further terms of the same torsion whatever keys
// they carry; rows whose keys differ from the first are library slips, listed in `slips`.
std::vector<std::string> g_slips;
std::vector<std::vector<const Row*>> torsion_groups(const std::vector<Row>& rows) {
  std::vector<std::vector<const Row*>> g;
  bool open = false;
  for (const auto& r : rows) {
    if (open && !g.empty()) {
      if (g.back().front()->keys != r.keys) g_slips.push_back(rule_name(g.back().front()->keys) + " continued by " + rule_name(r.keys) + " (line " + std::to_string(r.line) + ")");
      g.back().push_back(&r);
    } else {
      g.push_back({&r});
    }
    open = r.v.size() >= 2 && r.v[1] < 0;
  }
  return g;
}

// The same grouping when the multiplicity is in column `ncol` (AMBER: fsc Kchi gamma n ...).
std::vector<std::vector<const Row*>> torsion_groups_n(const std::vector<Row>& rows, size_t ncol) {
  std::vector<std::vector<const Row*>> g;
  bool open = false;
  for (const auto& r : rows) {
    if (open && !g.empty()) {
      if (g.back().front()->keys != r.keys) g_slips.push_back(rule_name(g.back().front()->keys) + " continued by " + rule_name(r.keys) + " (line " + std::to_string(r.line) + ")");
      g.back().push_back(&r);
    } else {
      g.push_back({&r});
    }
    open = r.v.size() > ncol && r.v[ncol] < 0;
  }
  return g;
}

}  // namespace

FFDef import_dlfield(const std::string& par, const std::string& sf, const std::string& bci) {
  g_slips.clear();
  const ParFile pf = read_par(par);
  const std::string fam = upper(pf.potential);
  const double E = energy_factor(pf.unit);
  FFDef ff;
  ff.name = pf.potential;
  ff.source = "DL_FIELD library: " + par;
  ff.references.push_back("C. W. Yong, DL_FIELD - a force field and model development tool for DL_POLY, STFC Daresbury Laboratory");
  ff.equivalence = "replace";

  // ---- family interpretation (functional forms and columns) ----
  enum class Fam { ClassII, CVFF, OPLS, AMBER, CHARMM, DREIDING, GROMOS, TRAPPE, TYPED } family;
  if (fam == "PCFF" || fam == "COMPASS") family = Fam::ClassII;
  else if (fam == "CVFF") family = Fam::CVFF;
  else if (fam.rfind("OPLS", 0) == 0) family = Fam::OPLS;
  else if (fam.rfind("AMBER", 0) == 0) family = Fam::AMBER;
  else if (fam.rfind("CHARMM", 0) == 0) family = Fam::CHARMM;
  else if (fam == "DREIDING") family = Fam::DREIDING;
  else if (fam == "G54A7" || fam.rfind("GROMOS", 0) == 0) family = Fam::GROMOS;
  else if (fam.rfind("TRAPPE", 0) == 0) family = Fam::TRAPPE;
  else if (fam.rfind("INORGANIC", 0) == 0 || fam == "MISC_FF") family = Fam::TYPED;
  else throw FFError("DL_FIELD family '" + pf.potential + "' is not interpreted yet (" + par + ")");

  auto rows = [&](const char* s) -> const std::vector<Row>& {
    static const std::vector<Row> none;
    auto it = pf.sec.find(s);
    return it == pf.sec.end() ? none : it->second;
  };
  auto need = [&](const Row& r, size_t n, const char* what) {
    if (r.v.size() < n) throw FFError(par + ":" + std::to_string(r.line) + ": " + what + " needs " + std::to_string(n) + " values");
  };
  // column of a named value in a section (the header lists the columns after the keys); `def` when not named
  auto col = [&](const char* s, std::initializer_list<const char*> names, int def) {
    auto it = pf.header.find(s);
    if (it == pf.header.end()) return def;
    for (size_t k = 0; k < it->second.size(); ++k)
      for (const char* n : names)
        if (it->second[k] == n) return int(k);
    return def;
  };
  auto rule = [&](const Row& r, const std::string& style, std::vector<double> p) {
    FFRule x;
    for (const auto& k : r.keys) x.match.push_back(key(k));
    x.name = rule_name(r.keys);
    x.style = style;
    x.params = std::move(p);
    x.comment = ref_of(r);
    return x;
  };

  if (family == Fam::ClassII) {
    // PCFF / COMPASS as DL_FIELD keeps them: the diagonal class II terms only (DL_POLY has no cross terms).
    ff.pair_style = "lj/class2/coul/long";
    ff.bond_style = "class2";
    ff.angle_style = "quartic";
    ff.dihedral_style = "fourier";
    ff.improper_style = "inversion";
    ff.mixing = "sixthpower";
    ff.special_lj[2] = ff.special_coul[2] = 1.0;
    ff.improper_order = "center2_sorted";
    ff.improper_matched_order = true;
    ff.notes.push_back("DL_FIELD's " + pf.potential + " has no class II cross terms (bond-bond, bond-angle, torsion couplings, "
                       "angle-angle): it is the diagonal part only, as DL_POLY can evaluate. For the full force field use "
                       "an Accelrys .frc conversion.");
    for (const auto& r : rows("BOND")) { need(r, 4, "BOND"); ff.bonds.push_back(rule(r, "class2", {r.v[0], r.v[1] * E, r.v[2] * E, r.v[3] * E})); }
    for (const auto& r : rows("ANGLE")) { need(r, 4, "ANGLE"); ff.angles.push_back(rule(r, "quartic", {r.v[0], r.v[1] * E, r.v[2] * E, r.v[3] * E})); }
    for (const auto& g : torsion_groups(rows("DIHEDRAL"))) {
      // V = K [1 − cos(nφ − φ0)] = K [1 + cos(nφ − φ0 − 180°)]
      std::vector<double> p = {0};
      for (const Row* r : g) {
        need(*r, 3, "DIHEDRAL");
        const int n = int(std::lround(std::fabs(r->v[1])));
        if (n == 0 && r->v[0] == 0) continue;   // n = 0 with K ≠ 0 is a constant term (TraPPE c0): kept
        p.push_back(r->v[0] * E); p.push_back(n); p.push_back(std::fmod(r->v[2] + 180.0, 360.0));
        p[0] += 1;
      }
      ff.dihedrals.push_back(rule(*g.front(), "fourier", p));
    }
    for (const auto& r : rows("INVERSION")) { need(r, 2, "INVERSION"); ff.impropers.push_back(rule(r, "inversion", {r.v[0] * E, r.v[1]})); }
    for (const auto& r : rows("VDW")) {   // Rmin eps: lj/class2 σ is the minimum
      need(r, 2, "VDW");
      FFRule x = rule(r, "", {r.v[1] * E, r.v[0]});
      ff.pairs.push_back(x);
    }
  } else if (family == Fam::OPLS) {
    // OPLS-AA family (Jorgensen): harmonic K (r − r0)², K (θ − θ0)², OPLS Fourier torsions V1..V3, improper as the
    // torsion ½V2 [1 − cos 2φ] with the centre third (AMBER order), geometric mixing, 1-4 scaled by ½.
    ff.pair_style = "lj/cut/coul/long";
    ff.bond_style = "harmonic";
    ff.angle_style = "harmonic";
    ff.dihedral_style = "opls";
    ff.improper_style = "cvff";
    ff.mixing = "geometric";
    ff.special_lj[2] = ff.special_coul[2] = 0.5;
    // DL_FIELD tries the outer atoms in permutation order from sorted, accepting a key read either way
    ff.improper_order = "center3_sorted";
    ff.improper_matched_order = true;
    ff.improper_reversible = true;
    const int bk = col("BOND", {"k", "K", "k0"}, 0), bb = col("BOND", {"b0"}, 1);
    const int ak = col("ANGLE", {"0.5K", "K", "0.5k", "k"}, 0), at = col("ANGLE", {"theta0"}, 1);
    const double kscale = E == 1.0 ? 1.0 : 0.5;   // kJ/mol libraries give the DL_POLY force constant k (E = ½ k Δ²)
    if (E != 1.0) ff.notes.push_back("force constants in " + pf.unit + " are taken as DL_POLY's k (E = ½kΔ²); unverified against DL_FIELD");
    for (const auto& r : rows("BOND")) { need(r, 2, "BOND"); ff.bonds.push_back(rule(r, "harmonic", {r.v[bk] * E * kscale, r.v[bb]})); }
    for (const auto& r : rows("ANGLE")) { need(r, 2, "ANGLE"); ff.angles.push_back(rule(r, "harmonic", {r.v[ak] * E * kscale, r.v[at]})); }
    for (const auto& r : rows("DIHEDRAL")) {
      need(r, 3, "DIHEDRAL");
      // V1 V2 V3 only: a fourth number on these rows is a flag of the library (the header lists V1 V2 V3)
      const bool v4 = col("DIHEDRAL", {"V4"}, -1) == 3;
      ff.dihedrals.push_back(rule(r, "opls", {r.v[0] * E, r.v[1] * E, r.v[2] * E, v4 && r.v.size() > 3 ? r.v[3] * E : 0.0}));
    }
    for (const auto& r : rows("IMPROPER")) {
      // ½V1 [1 + cos φ] + ½V2 [1 − cos 2φ] + ½V3 [1 + cos 3φ]: DL_FIELD keeps the V2 term (the others are zero in these files)
      need(r, 3, "IMPROPER");
      if (r.v[0] != 0 || r.v[2] != 0)
        ff.notes.push_back("improper " + rule_name(r.keys) + ": V1 / V3 terms are not converted (only V2)");
      ff.impropers.push_back(rule(r, "cvff", {0.5 * r.v[1] * E, -1, 2}));
    }
    const int vs = col("VDW", {"sigma/angstrom", "sigma"}, 0), ve = col("VDW", {"eps/kcal/mol", "eps/kJ/mol", "eps"}, 1);
    for (const auto& r : rows("VDW")) { need(r, 2, "VDW"); ff.pairs.push_back(rule(r, "", {r.v[ve] * E, r.v[vs]})); }
  } else if (family == Fam::AMBER) {
    // AMBER / GAFF (Cornell et al. 1995; Wang et al. 2004): K (r − r0)², K (θ − θ0)², torsions (Kchi / fsc) [1 + cos(nφ − γ)]
    // with 1-4 scaling SCEE / SCNB, impropers Kchi [1 + cos(nφ − γ)] with the centre third, σ = 2 (R/2) / 2^(1/6),
    // Lorentz–Berthelot mixing.
    ff.pair_style = "lj/cut/coul/long";
    ff.bond_style = "harmonic";
    ff.angle_style = "harmonic";
    ff.dihedral_style = "fourier";
    ff.improper_style = "cvff";
    ff.mixing = "arithmetic";
    ff.special_lj[2] = 0.5;
    ff.special_coul[2] = 0.8333333333;
    // centre third; outer atoms tried in bond-list order, key read either way (the order of the two wildcard atoms is a
    // convention: AMBER itself leaves it to the builder)
    ff.improper_order = "center3_sorted";
    ff.improper_matched_order = true;
    ff.improper_reversible = true;
    const int bk = col("BOND", {"k", "K"}, 0), bb = col("BOND", {"b0"}, 1);
    const int ak = col("ANGLE", {"0.5K", "K"}, 0), at = col("ANGLE", {"theta0"}, 1);
    for (const auto& r : rows("BOND")) { need(r, 2, "BOND"); ff.bonds.push_back(rule(r, "harmonic", {r.v[bk] * E, r.v[bb]})); }
    for (const auto& r : rows("ANGLE")) { need(r, 2, "ANGLE"); ff.angles.push_back(rule(r, "harmonic", {r.v[ak] * E, r.v[at]})); }
    std::set<std::pair<double, double>> scale14;
    for (const auto& g : torsion_groups_n(rows("DIHEDRAL"), 3)) {
      // fsc Kchi gamma n SCEE SCNB
      std::vector<double> p = {0};
      for (const Row* r : g) {
        need(*r, 4, "DIHEDRAL");
        const int n = int(std::lround(std::fabs(r->v[3])));
        if (r->v.size() >= 6) scale14.insert({r->v[4], r->v[5]});
        if (n == 0 || r->v[0] == 0) continue;
        p.push_back(r->v[1] / r->v[0] * E); p.push_back(n); p.push_back(r->v[2]);
        p[0] += 1;
      }
      ff.dihedrals.push_back(rule(*g.front(), "fourier", p));
    }
    if (scale14.size() > 1) {
      std::string l;
      for (const auto& [e, n] : scale14) l += " (" + std::to_string(e) + ", " + std::to_string(n) + ")";
      ff.notes.push_back("torsions carry different 1-4 scalings (SCEE, SCNB):" + l + "; CAPS applies 1/1.2 and 1/2 to all 1-4 pairs, as DL_FIELD does for LAMMPS");
    }
    for (const auto& r : rows("IMPROPER")) {   // Kchi gamma n
      need(r, 3, "IMPROPER");
      const double g = std::fmod(std::fabs(r.v[1]) + 360, 360);
      if (std::fabs(g) > 1e-6 && std::fabs(g - 180) > 1e-6) ff.notes.push_back("improper " + rule_name(r.keys) + ": phase " + std::to_string(r.v[1]) + " is not 0 or 180");
      ff.impropers.push_back(rule(r, "cvff", {r.v[0] * E, std::fabs(g - 180) < 1e-6 ? -1.0 : 1.0, std::fabs(r.v[2])}));
    }
    const int vr = col("VDW", {"R/2"}, 0), ve = col("VDW", {"Epsilon", "eps"}, 1);
    for (const auto& r : rows("VDW")) {
      need(r, 2, "VDW");
      ff.pairs.push_back(rule(r, "", {r.v[ve] * E, 2 * r.v[vr] / std::pow(2.0, 1.0 / 6)}));
    }
  } else if (family == Fam::CHARMM) {
    // CHARMM (Brooks et al.; MacKerell et al.): K (r − r0)², K (θ − θ0)² + K_UB (r13 − S0)², Σ K [1 + cos(nφ − δ)],
    // impropers K (ψ − ψ0)² with the centre first, σ = 2 (Rmin/2) / 2^(1/6) with separate 1-4 values, NBFIX pairs,
    // Lorentz–Berthelot mixing, 1-4 electrostatics and van der Waals unscaled (the 1-4 LJ uses its own parameters).
    ff.pair_style = "lj/charmm/coul/long";
    ff.bond_style = "harmonic";
    ff.angle_style = "charmm";
    ff.dihedral_style = "fourier";
    ff.improper_style = "harmonic";
    ff.mixing = "arithmetic";
    ff.special_lj[2] = ff.special_coul[2] = 1.0;
    ff.improper_order = "center1_sorted";
    ff.improper_matched_order = true;
    ff.improper_reversible = true;
    ff.improper_all_explicit = true;
    for (const auto& r : rows("BOND")) { need(r, 2, "BOND"); ff.bonds.push_back(rule(r, "harmonic", {r.v[0] * E, r.v[1]})); }
    for (const auto& r : rows("ANGLE")) {
      need(r, 2, "ANGLE");
      const double kub = r.v.size() > 2 ? r.v[2] * E : 0, s0 = r.v.size() > 3 ? r.v[3] : 0;
      ff.angles.push_back(kub != 0 ? rule(r, "charmm", {r.v[0] * E, r.v[1], kub, s0}) : rule(r, "harmonic", {r.v[0] * E, r.v[1]}));
    }
    for (const auto& g : torsion_groups(rows("DIHEDRAL"))) {   // Kchi n delta
      std::vector<double> p = {0};
      for (const Row* r : g) {
        need(*r, 3, "DIHEDRAL");
        const int n = int(std::lround(std::fabs(r->v[1])));
        if (n == 0 && r->v[0] == 0) continue;   // n = 0 with K ≠ 0 is a constant term (TraPPE c0): kept
        p.push_back(r->v[0] * E); p.push_back(n); p.push_back(r->v[2]);
        p[0] += 1;
      }
      ff.dihedrals.push_back(rule(*g.front(), "fourier", p));
    }
    for (const auto& r : rows("IMPROPER")) { need(r, 2, "IMPROPER"); ff.impropers.push_back(rule(r, "harmonic", {r.v[0] * E, r.v[1]})); }
    const double s6 = std::pow(2.0, 1.0 / 6);
    for (const auto& r : rows("VDW")) {   // Eps Rmin/2 [Eps_1-4 Rmin/2,1-4]
      need(r, 2, "VDW");
      FFRule x = rule(r, "", {std::fabs(r.v[0]) * E, 2 * r.v[1] / s6});
      if (r.v.size() >= 4 && (r.v[2] != 0 || r.v[3] != 0)) x.cross["lj14"] = {std::fabs(r.v[2]) * E, 2 * r.v[3] / s6};
      ff.pairs.push_back(x);
    }
    for (const auto& r : rows("VDW_FIX")) {   // NBFIX: a b lj Emin Rmin
      if (r.keys.size() < 2 || r.v.size() < 2) continue;
      ff.pairs.push_back(rule(r, "", {std::fabs(r.v[0]) * E, r.v[1] / s6}));
    }
    if (fam == "CHARMM19") ff.notes.push_back("CHARMM19 is a polar-hydrogen (united-atom) force field");
  } else if (family == Fam::DREIDING) {
    // DREIDING (Mayo, Olafson, Goddard 1990) as DL_FIELD keeps it: ½k (r − r0)², ½k (θ − θ0)² (Morse and cosine
    // variants kept as alternatives), ½V [1 − cos(n(φ − φ0))], planar inversion K (1 − cos ω), LJ from Rmin.
    ff.pair_style = "lj/cut/coul/long";
    ff.bond_style = "harmonic";
    ff.angle_style = "harmonic";
    ff.dihedral_style = "fourier";
    ff.improper_style = "planar";
    ff.mixing = "arithmetic";
    ff.special_lj[2] = ff.special_coul[2] = 1.0;
    ff.improper_order = "center1_sorted";
    ff.improper_matched_order = true;
    for (const auto& r : rows("BOND")) {   // k b0 De a
      need(r, 2, "BOND");
      FFRule x = rule(r, "harmonic", {0.5 * r.v[0] * E, r.v[1]});
      if (r.v.size() >= 4 && r.v[2] > 0) x.cross["morse"] = {r.v[2] * E, r.v[3], r.v[1]};
      ff.bonds.push_back(x);
    }
    for (const auto& r : rows("ANGLE")) {   // k C theta0
      need(r, 3, "ANGLE");
      // linear centres use K (1 + cos θ) with the full k (DREIDING's form for θ0 = 180°, as DL_FIELD writes it)
      FFRule x = std::fabs(r.v[2] - 180.0) < 1e-6 ? rule(r, "cosine", {r.v[0] * E}) : rule(r, "harmonic", {0.5 * r.v[0] * E, r.v[2]});
      if (r.v[1] > 0) x.cross["cosine/squared"] = {0.5 * r.v[1] * E, r.v[2]};
      ff.angles.push_back(x);
    }
    for (const auto& g : torsion_groups(rows("DIHEDRAL"))) {   // V n chi0: ½V [1 − cos(n(φ − φ0))] = ½V [1 + cos(nφ − nφ0 − 180)]
      std::vector<double> p = {0};
      for (const Row* r : g) {
        need(*r, 3, "DIHEDRAL");
        const int n = int(std::lround(std::fabs(r->v[1])));
        if (n == 0 && r->v[0] == 0) continue;   // n = 0 with K ≠ 0 is a constant term (TraPPE c0): kept
        p.push_back(0.5 * r->v[0] * E); p.push_back(n); p.push_back(std::fmod(std::fmod(n * r->v[2] + 180.0, 360.0) + 360.0, 360.0));
        p[0] += 1;
      }
      ff.dihedrals.push_back(rule(*g.front(), "fourier", p));
    }
    // DREIDING (Mayo et al. 1990): every planar centre (_2, _R) with three bonds has K (1 − cos ω), K = 40 kcal/mol;
    // specific entries of the library come after these generic ones and so take precedence
    for (const char* c : {"*_2", "*_R"}) {
      FFRule g;
      g.name = std::string(c) + " planar centre (DREIDING rule)";
      g.match = {c, "*", "*", "*"};
      g.style = "planar";
      g.params = {40.0, 0.0};   // kcal/mol
      g.comment = "DREIDING generic planar inversion";
      ff.impropers.push_back(g);
    }
    for (const auto& r : rows("INVERSION")) { need(r, 1, "INVERSION"); ff.impropers.push_back(rule(r, "planar", {r.v[0] * E, 0.0})); }
    const double s6 = std::pow(2.0, 1.0 / 6);
    for (const auto& r : rows("VDW")) { need(r, 2, "VDW"); ff.pairs.push_back(rule(r, "", {r.v[1] * E, r.v[0] / s6})); }
    ff.notes.push_back("DL_FIELD also places inversions where its molecule templates ask (e.g. amide N typed N_3); CAPS applies "
                       "DREIDING's rule: planar (_2, _R) centres only");
    ff.notes.push_back("DREIDING hydrogen-bond terms are not converted; Morse bonds and cosine angles are kept as alternatives in each rule");
  } else if (family == Fam::GROMOS) {
    // GROMOS 54A7 (Schmid et al. 2011), kJ/mol and nm: quartic bonds ¼kb (b² − b0²)² (func 2), cosine-harmonic angles
    // ½kθ (cos θ − cos θ0)² (func 2), k [1 + cos(nφ − δ)], harmonic impropers ½kξ (ξ − ξ0)², C6 / C12 with separate
    // 1-4 values (geometric combination), 1-4 electrostatics unscaled.
    ff.pair_style = "lj/cut/coul/long";
    ff.bond_style = "gromos";
    ff.angle_style = "cosine/squared";
    ff.dihedral_style = "fourier";
    ff.improper_style = "harmonic";
    ff.mixing = "geometric";
    ff.special_lj[2] = ff.special_coul[2] = 1.0;
    ff.improper_order = "center1_sorted";
    ff.improper_matched_order = true;
    ff.improper_reversible = true;
    const double nm = 10.0;   // nm → Å
    for (const auto& r : rows("BOND")) {   // func b0 k
      need(r, 3, "BOND");
      const int f = int(r.v[0]);
      if (f == 2) ff.bonds.push_back(rule(r, "gromos", {0.25 * r.v[2] * E / (nm * nm * nm * nm), r.v[1] * nm}));
      else ff.bonds.push_back(rule(r, "harmonic", {0.5 * r.v[2] * E / (nm * nm), r.v[1] * nm}));
    }
    for (const auto& r : rows("ANGLE")) {   // func theta0 k
      need(r, 3, "ANGLE");
      const int f = int(r.v[0]);
      if (f == 2) ff.angles.push_back(rule(r, "cosine/squared", {0.5 * r.v[2] * E, r.v[1]}));
      else ff.angles.push_back(rule(r, "harmonic", {0.5 * r.v[2] * E, r.v[1]}));
    }
    for (const auto& g : torsion_groups_n(rows("DIHEDRAL"), 3)) {   // func delta k n
      std::vector<double> p = {0};
      for (const Row* r : g) {
        need(*r, 4, "DIHEDRAL");
        const int n = int(std::lround(std::fabs(r->v[3])));
        if (n == 0 && r->v[2] == 0) continue;
        p.push_back(r->v[2] * E); p.push_back(n); p.push_back(r->v[1]);
        p[0] += 1;
      }
      ff.dihedrals.push_back(rule(*g.front(), "fourier", p));
    }
    for (const auto& r : rows("IMPROPER")) { need(r, 3, "IMPROPER"); ff.impropers.push_back(rule(r, "harmonic", {0.5 * r.v[2] * E, r.v[1]})); }
    for (const auto& r : rows("VDW")) {   // C6 C12 [C6_14 C12_14], kJ/mol nm^6 / nm^12
      need(r, 2, "VDW");
      auto lj = [&](double c6, double c12) -> std::vector<double> {
        if (c6 <= 0 || c12 <= 0) return {0.0, 0.0};
        return {c6 * c6 / (4 * c12) * E, std::pow(c12 / c6, 1.0 / 6) * nm};
      };
      FFRule x = rule(r, "", lj(r.v[0], r.v[1]));
      if (r.v.size() >= 4) x.cross["lj14"] = lj(r.v[2], r.v[3]);
      ff.pairs.push_back(x);
    }
    ff.notes.push_back("GROMOS exclusion lists beyond third neighbours (aromatic rings) are not converted");
  } else if (family == Fam::TRAPPE) {
    // TraPPE (Siepmann et al.), energies in K: bonds and angles in DL_POLY's ½k convention, torsions
    // Kchi [1 + cos(nφ − γ)] (EH: Kchi / fsc), impropers Kchi [1 + cos(nφ − γ)], σ ε, Lorentz–Berthelot, no 1-4 terms.
    ff.pair_style = "lj/cut/coul/long";
    ff.bond_style = "harmonic";
    ff.angle_style = "harmonic";
    ff.dihedral_style = "fourier";
    ff.improper_style = "cvff";
    ff.mixing = "arithmetic";
    ff.special_lj[2] = ff.special_coul[2] = 0.0;
    ff.improper_order = "center3_sorted";
    ff.improper_matched_order = true;
    ff.improper_reversible = true;
    const bool fsc = col("DIHEDRAL", {"fsc"}, -1) == 0;
    // DL_FIELD reads TraPPE-UA constants as DL_POLY's k (E = ½kΔ²) and TraPPE-EH ones as K (E = KΔ²)
    const double kf = fam == "TRAPPE_EH" ? 1.0 : 0.5;
    for (const auto& r : rows("BOND")) { need(r, 2, "BOND"); ff.bonds.push_back(rule(r, "harmonic", {kf * r.v[0] * E, r.v[1]})); }
    for (const auto& r : rows("ANGLE")) { need(r, 2, "ANGLE"); ff.angles.push_back(rule(r, "harmonic", {kf * r.v[0] * E, r.v[1]})); }
    const size_t o = fsc ? 1 : 0;
    for (const auto& g : torsion_groups_n(rows("DIHEDRAL"), o + 2)) {
      std::vector<double> p = {0};
      for (const Row* r : g) {
        need(*r, o + 3, "DIHEDRAL");
        const int n = int(std::lround(std::fabs(r->v[o + 2])));
        if (n == 0 && r->v[o] == 0) continue;   // n = 0 with K ≠ 0 is TraPPE's constant c0: kept
        p.push_back(r->v[o] / (fsc && r->v[0] != 0 ? r->v[0] : 1.0) * E); p.push_back(n); p.push_back(r->v[o + 1]);
        p[0] += 1;
      }
      ff.dihedrals.push_back(rule(*g.front(), "fourier", p));
    }
    for (const auto& r : rows("IMPROPER")) {   // Kchi gamma n
      need(r, 3, "IMPROPER");
      ff.impropers.push_back(rule(r, "cvff", {r.v[0] * E, std::fabs(std::fmod(std::fabs(r.v[1]), 360.0) - 180) < 1e-6 ? -1.0 : 1.0, std::fabs(r.v[2])}));
    }
    for (const auto& r : rows("VDW")) { need(r, 2, "VDW"); ff.pairs.push_back(rule(r, "", {r.v[1] * E, r.v[0]})); }
  } else if (family == Fam::TYPED) {
    // DL_FIELD's inorganic and miscellaneous libraries name the DL_POLY form on every row.
    ff.pair_style = "buck/coul/long";
    ff.bond_style = "harmonic";
    ff.angle_style = "harmonic";
    ff.dihedral_style = "fourier";
    ff.improper_style = "fourier";
    ff.mixing = "geometric";
    ff.special_lj[2] = ff.special_coul[2] = 1.0;
    ff.improper_order = "center1_sorted";
    ff.improper_matched_order = true;
    // DL_FIELD reads the first four characters of the form's name (so the library's "buckinghsm" is Buckingham)
    auto lower = [&](std::string k) {
      for (auto& c : k) c = char(std::tolower(static_cast<unsigned char>(c)));
      if (!k.empty() && k[0] == '-') k.erase(0, 1);
      static const std::map<std::string, std::string> names = {{"buck", "buckingham"}, {"mors", "morse"}, {"lj12", "lj126"},
                                                               {"harm", "harmonic"}, {"cosi", "cosine"}, {"cos3", "cos3"},
                                                               {"quar", "quartic"}, {"coul", "coulomb"}};
      auto it = names.find(k.substr(0, 4));
      if (it == names.end()) return k;
      if (k != it->second && k.size() > 4) ff.notes.push_back("form name '" + k + "' read as " + it->second + " (first four characters, as DL_FIELD)");
      return it->second;
    };
    std::set<std::string> skipped;
    for (const auto& r : rows("BOND")) {
      const std::string k = lower(r.kind);
      if (k == "harmonic" && r.v.size() >= 2) ff.bonds.push_back(rule(r, "harmonic", {0.5 * r.v[0] * E, r.v[1]}));
      else if (k == "morse" && r.v.size() >= 3) ff.bonds.push_back(rule(r, "morse", {r.v[0] * E, r.v[2], r.v[1]}));
      else if (k == "quartic" && r.v.size() >= 4) ff.bonds.push_back(rule(r, "class2", {r.v[0], r.v[1] * E, r.v[2] * E, r.v[3] * E}));
      else skipped.insert("bond " + r.kind);
    }
    for (const auto& r : rows("ANGLE")) {
      const std::string k = lower(r.kind);
      if (k == "harmonic" && r.v.size() >= 2) ff.angles.push_back(rule(r, "harmonic", {0.5 * r.v[0] * E, r.v[1]}));
      else if (k == "quartic" && r.v.size() >= 4) ff.angles.push_back(rule(r, "quartic", {r.v[0], r.v[1] * E, r.v[2] * E, r.v[3] * E}));
      else skipped.insert("angle " + r.kind);
    }
    auto torsion_rule = [&](const std::vector<const Row*>& g, bool improper) {
      std::vector<double> p = {0};
      for (const Row* r : g) {
        const std::string k = lower(r->kind);
        if (k == "cosine" && r->v.size() >= 3) {   // A delta m: A [1 + cos(mφ − δ)]
          const int n = int(std::lround(std::fabs(r->v[2])));
          if (n == 0 && r->v[0] == 0) continue;
          p.push_back(r->v[0] * E); p.push_back(n); p.push_back(r->v[1]); p[0] += 1;
        } else if (k == "cos3" && r->v.size() >= 3) {   // ½[A1(1 + cos φ) + A2(1 − cos 2φ) + A3(1 + cos 3φ)]
          const double t[3][3] = {{r->v[0], 1, 0}, {r->v[1], 2, 180}, {r->v[2], 3, 0}};
          for (const auto& q : t)
            if (q[0] != 0) { p.push_back(0.5 * q[0] * E); p.push_back(q[1]); p.push_back(q[2]); p[0] += 1; }
        } else {
          skipped.insert(std::string(improper ? "improper " : "dihedral ") + r->kind);
        }
      }
      return rule(*g.front(), "fourier", p);
    };
    std::vector<std::vector<const Row*>> dg;
    for (const auto& g : torsion_groups_n(rows("DIHEDRAL"), 2)) dg.push_back(g);
    for (const auto& g : dg) ff.dihedrals.push_back(torsion_rule(g, false));
    for (const auto& r : rows("IMPROPER")) ff.impropers.push_back(torsion_rule({&r}, true));
    for (const auto& r : rows("INVERSION")) {
      if (r.v.size() >= 2) {
        ff.impropers.push_back(rule(r, "inversion", {r.v[0] * E, r.v[1]}));
        ff.improper_order = "center2_sorted";
        ff.improper_style = "inversion";
      }
    }
    for (const auto& r : rows("SHELL")) {   // core shell k2 [k4]: ½k2 r² + (1/24) k4 r⁴
      if (r.keys.size() < 2 || r.v.empty()) continue;
      FFRule x = rule(r, "class2", {0.0, 0.5 * r.v[0] * E, 0.0, (r.v.size() > 1 ? r.v[1] : 0.0) * E / 24});
      x.comment = "core-shell spring · " + x.comment;
      ff.bonds.push_back(x);
    }
    for (const auto& r : rows("VDW")) {
      const std::string k = lower(r.kind);
      std::vector<std::string> keys = r.keys;
      const bool self = keys.size() == 1 || keys[1] == "X" || keys[1] == keys[0];
      FFRule x;
      if (k == "buckingham" && r.v.size() >= 3) x = rule(r, "buck", {r.v[0] * E, r.v[1], r.v[2] * E});
      else if (k == "lj" && r.v.size() >= 2) x = rule(r, "", {r.v[0] * E, r.v[1]});
      else if (k == "lj126" && r.v.size() >= 2) {
        const double A = r.v[0] * E, B = r.v[1] * E;
        x = rule(r, "", {A > 0 ? B * B / (4 * A) : 0.0, B > 0 ? std::pow(A / B, 1.0 / 6) : 0.0});
      } else if (k == "morse" && r.v.size() >= 3) x = rule(r, "morse", {r.v[0] * E, r.v[2], r.v[1]});
      else { skipped.insert("vdw " + r.kind); continue; }
      if (self) x.match = {x.match[0]};
      if (self && !x.style.empty()) x.match.push_back(x.match[0]);   // non-LJ forms need an explicit pair
      ff.pairs.push_back(x);
    }
    for (const auto& k : skipped) ff.notes.push_back("not converted: " + k + " (DL_POLY form not available in CAPS)");
  } else {
    // CVFF
    ff.pair_style = "lj/cut/coul/long";
    ff.bond_style = "harmonic";
    ff.angle_style = "harmonic";
    ff.dihedral_style = "fourier";
    ff.improper_style = "cvff";
    ff.mixing = "geometric";
    ff.special_lj[2] = ff.special_coul[2] = 1.0;
    ff.improper_order = "center2_sorted";
    ff.improper_matched_order = true;
    ff.wildcard_torsion_scaling = "msi2lmp";
    for (const auto& r : rows("BOND")) { need(r, 2, "BOND"); ff.bonds.push_back(rule(r, "harmonic", {r.v[1] * E, r.v[0]})); }
    for (const auto& r : rows("ANGLE")) { need(r, 2, "ANGLE"); ff.angles.push_back(rule(r, "harmonic", {r.v[1] * E, r.v[0]})); }
    for (const auto& g : torsion_groups(rows("DIHEDRAL"))) {
      // V = K [1 + cos(nφ − φ0)]
      std::vector<double> p = {0};
      for (const Row* r : g) {
        need(*r, 3, "DIHEDRAL");
        const int n = int(std::lround(std::fabs(r->v[1])));
        if (n == 0 && r->v[0] == 0) continue;   // n = 0 with K ≠ 0 is a constant term (TraPPE c0): kept
        p.push_back(r->v[0] * E); p.push_back(n); p.push_back(r->v[2]);
        p[0] += 1;
      }
      ff.dihedrals.push_back(rule(*g.front(), "fourier", p));
    }
    for (const auto& r : rows("IMPROPER")) {
      // V = K [1 + cos(nχ − χ0)], χ0 0 or 180 → LAMMPS cvff K d n
      need(r, 3, "IMPROPER");
      ff.impropers.push_back(rule(r, "cvff", {r.v[0] * E, std::cos(r.v[2] * 3.14159265358979323846 / 180) >= 0 ? 1.0 : -1.0, r.v[1]}));
    }
    for (const auto& r : rows("VDW")) {   // A/r¹² − B/r⁶
      need(r, 2, "VDW");
      const double A = r.v[0] * E, B = r.v[1] * E;
      const double eps = A > 0 ? B * B / (4 * A) : 0, sig = B > 0 ? std::pow(A / B, 1.0 / 6) : 0;
      ff.pairs.push_back(rule(r, "", {eps, sig}));
    }
  }
  for (auto* v : {&ff.bonds, &ff.angles, &ff.dihedrals, &ff.impropers}) order_rules(*v);
  ff.improper_max_neighbours = 3;   // DL_FIELD's out-of-plane terms are for planar (three-connected) centres

  // ---- atom types: .sf ATOM_TYPE (DL_FIELD type → key), and every key named in EQUIVALENCE or VDW ----
  std::map<std::string, size_t> tix;
  auto type_of = [&](const std::string& k) -> FFType& {
    auto [it, fresh] = tix.emplace(k, ff.types.size());
    if (fresh) { FFType t; t.name = k; ff.types.push_back(t); }
    return ff.types[it->second];
  };
  if (!sf.empty()) {
    std::ifstream in(sf);
    if (!in) throw FFError("cannot open " + sf);
    // every ATOM_TYPE entry per key first: DL_FIELD's libraries carry typos in the element or mass of single entries
    // (AMBER25 HC_alkyne → ha as C 12.0115, CVFF O_alcohol → oh at 12.0115, CGenFF N_nitro → NG2O1 at 12.007, ...)
    struct Entry { std::string dl; int z; double m; std::string remark; };
    std::map<std::string, std::vector<Entry>> entries;
    std::vector<std::string> order;
    std::string line;
    bool on = false;
    while (std::getline(in, line)) {
      line = strip_comments(line);
      std::istringstream ls(line);
      std::vector<std::string> t;
      for (std::string w; ls >> w;) t.push_back(w);
      if (t.empty()) continue;
      const std::string head = upper(t[0]);
      if (head == "ATOM_TYPE") { on = true; continue; }
      if (head == "END" || head == "MOLECULE" || head == "MOLECULE_TYPE") { on = false; continue; }
      if (!on || t.size() < 4 || !is_number(t[3])) continue;
      std::string remark;
      for (size_t k = 4; k < t.size(); ++k) remark += (remark.empty() ? "" : " ") + t[k];
      if (!entries.count(t[1])) order.push_back(t[1]);
      entries[t[1]].push_back({t[0], element_from_symbol(t[2]), std::stod(t[3]), remark});
    }
    // a mass equal to the standard mass of another element than the one given is a typo
    auto other_element_mass = [](double m, int z) {
      for (int e : {1, 6, 7, 8, 9, 14, 15, 16, 17})
        if (e != z && std::fabs(m - element(e).mass) < 0.005) return true;
      for (double dl : {1.00797, 12.0115, 14.0067, 15.9994})   // DL_FIELD's own values
        if (std::fabs(m - dl) < 1e-4 && std::fabs(m - element(z).mass) > 0.05) return true;
      return false;
    };
    for (const auto& key : order) {
      const auto& es = entries[key];
      FFType& ty = type_of(key);
      for (const auto& e : es)
        if (std::find(ty.aliases.begin(), ty.aliases.end(), e.dl) == ty.aliases.end()) ty.aliases.push_back(e.dl);
      // element: the majority (first on a tie)
      std::map<int, int> votes;
      for (const auto& e : es) votes[e.z]++;
      int z = es[0].z;
      for (const auto& [zz, n] : votes)
        if (n > votes[z]) z = zz;
      // mass: among entries of that element, the one nearest the element's standard mass
      const Entry* best = nullptr;
      for (const auto& e : es)
        if (e.z == z && (!best || std::fabs(e.m - element(z).mass) < std::fabs(best->m - element(z).mass))) best = &e;
      double m = best->m;
      ty.element = z;
      ty.description = best->remark;
      if (z > 0 && other_element_mass(m, z)) {
        ff.notes.push_back("type " + key + ": DL_FIELD gives mass " + std::to_string(m) + " for " + element(z).symbol + "; the element's mass " +
                           std::to_string(element(z).mass) + " is used");
        m = element(z).mass;
      }
      ty.mass = m;
      if (votes.size() > 1)
        ff.notes.push_back("type " + key + ": DL_FIELD entries disagree on the element; " + element(z).symbol + " (" + std::to_string(votes[z]) + " of " +
                           std::to_string(es.size()) + ") is used");
      bool differ = false;
      for (const auto& e : es) differ |= e.z == z && std::fabs(e.m - ty.mass) > 1e-3;
      if (differ) ff.notes.push_back("type " + key + ": DL_FIELD entries give different masses; " + std::to_string(ty.mass) + " is used");
    }
  }
  // "bond2_" etc. are second-choice equivalents, tried when the first lookup finds nothing
  const char* kinds[][2] = {{"vdw_", "vdw"}, {"bond2_", "bond2"}, {"angle2_", "angle2"}, {"dihedral2_", "dihedral2"},
                            {"torsion2_", "dihedral2"}, {"inv2_", "improper2"}, {"imp2_", "improper2"}, {"bond_", "bond"},
                            {"angle_", "angle"}, {"dihedral_", "dihedral"}, {"torsion_", "dihedral"},
                            {"inv_", "improper"}, {"imp_", "improper"}, {"inc_", "increment"}};
  // inv_ equivalences serve the impropers when the library's out-of-plane terms are INVERSION entries
  const bool inversions = family == Fam::ClassII || family == Fam::DREIDING || (!rows("INVERSION").empty() && rows("IMPROPER").empty());
  for (const auto& [k, targets] : pf.equivalence) {
    FFType& ty = type_of(k);
    for (const auto& w : targets)
      for (const auto& kd : kinds)
        if (w.rfind(kd[0], 0) == 0) {
          const std::string what = kd[1];
          if ((what == "improper" || what == "improper2") &&
              ((inversions && w.rfind("imp", 0) == 0) || (!inversions && w.rfind("inv", 0) == 0))) break;
          ty.equiv[what] = w.substr(std::string(kd[0]).size());
          break;
        }
  }
  // keys with parameters but no .sf entry still become types (element unknown until a structure gives it)
  for (const auto& r : rows("VDW"))
    if (r.keys.size() == 1) type_of(r.keys[0]);

  // ---- bond charge increments ----
  if (!bci.empty()) {
    std::ifstream in(bci);
    if (!in) throw FFError("cannot open " + bci);
    std::string line;
    bool on = false;
    while (std::getline(in, line)) {
      if (auto h = line.find('#'); h != std::string::npos) line = line.substr(0, h);
      std::istringstream ls(line);
      std::vector<std::string> t;
      for (std::string w; ls >> w;) t.push_back(w);
      if (t.empty()) continue;
      if (upper(t[0]) == "BOND_INCREMENTS") { on = true; continue; }
      if (upper(t[0]) == "END") { on = false; continue; }
      if (!on || t.size() < 6 || !is_number(t[4]) || !is_number(t[5])) continue;
      FFRule r;
      r.match = {glob_escape(t[2]), glob_escape(t[3])};
      r.name = t[2] + "-" + t[3];
      r.params = {std::stod(t[4]), std::stod(t[5])};
      r.comment = "ver " + t[0] + ", ref " + t[1];
      ff.bond_increments.push_back(r);
    }
    std::reverse(ff.bond_increments.begin(), ff.bond_increments.end());   // the first entry wins
  }
  // ---- automatic parameters (lib/supplementary/<ff>_*_auto): the Accelrys cff91_auto / cvff_auto tables ----
  {
    const std::string dir = par.substr(0, par.find_last_of('/') + 1) + "supplementary/";
    std::string prefix = pf.potential;
    for (auto& c : prefix) c = char(std::tolower(static_cast<unsigned char>(c)));
    auto table = [&](const std::string& name) {
      std::vector<std::vector<std::string>> out;
      std::ifstream in(dir + prefix + "_" + name);
      std::string line;
      while (std::getline(in, line)) {
        if (auto h = line.find('#'); h != std::string::npos) line = line.substr(0, h);
        if (auto h = line.find('!'); h != std::string::npos) line = line.substr(0, h);
        std::istringstream ls(line);
        std::vector<std::string> t;
        for (std::string w; ls >> w;) t.push_back(w);
        if (t.size() >= 3 && is_number(t[0])) out.push_back(t);
      }
      return out;
    };
    // Type NonB BondInct Bond AngleEnd AngleApex TorsionEnd TorsionCenter OOPEnd OOPCenter
    const auto eq = table("equivalence_auto");
    const char* cols[] = {"", "", "auto_bond", "auto_angle_end", "auto_angle_apex", "auto_torsion_end", "auto_torsion_center",
                          "auto_oop_end", "auto_oop_center"};
    for (const auto& t : eq) {
      if (t.size() < 12) continue;
      FFType& ty = type_of(t[2]);
      for (int c = 2; c < 9; ++c) ty.equiv[cols[c]] = t[3 + c];
    }
    // "*7" is a wildcard whose number sets precedence (the lowest number wins)
    auto auto_rule = [&](const std::vector<std::string>& t, int nk, const std::string& style, std::vector<double> p, int& prec) {
      FFRule x;
      prec = 0;
      for (int k = 0; k < nk; ++k) {
        const std::string& w = t[2 + k];
        if (w[0] == '*') {
          x.match.push_back("*");
          if (w.size() > 1) prec = std::max(prec, std::atoi(w.c_str() + 1));
        } else {
          x.match.push_back(glob_escape(w));
        }
      }
      x.name = "auto";
      for (int k = 0; k < nk; ++k) x.name += " " + t[2 + k];
      x.style = style;
      x.params = std::move(p);
      x.comment = "auto, ver " + t[0] + ", ref " + t[1];
      return x;
    };
    auto load = [&](const std::vector<std::vector<std::string>>& rows, int nk, auto make, std::vector<FFRule>& dst) {
      std::vector<std::pair<std::tuple<int, int, size_t>, FFRule>> v;
      for (size_t i = 0; i < rows.size(); ++i) {
        const auto& t = rows[i];
        if (int(t.size()) < 2 + nk + 2) continue;
        std::vector<double> num;
        for (size_t k = 2 + nk; k < t.size() && is_number(t[k]); ++k) num.push_back(std::stod(t[k]));
        int prec = 0;
        FFRule x = auto_rule(t, nk, "", {}, prec);
        if (!make(x, num)) continue;
        const int wild = int(std::count(x.match.begin(), x.match.end(), std::string("*")));
        v.push_back({{wild, prec, i}, x});
      }
      // last match wins: most wildcards first, then highest precedence number, then later entries
      std::stable_sort(v.begin(), v.end(), [](const auto& a, const auto& b) {
        const auto& [wa, pa, ia] = a.first;
        const auto& [wb, pb, ib] = b.first;
        if (wa != wb) return wa > wb;
        if (pa != pb) return pa > pb;
        return ia > ib;
      });
      for (auto& [k, x] : v) dst.push_back(x);
    };
    const bool c2 = family == Fam::ClassII;
    auto bonds = table("bond_auto");
    if (bonds.empty()) bonds = table("bond_quadratic_auto");
    load(bonds, 2, [&](FFRule& x, const std::vector<double>& n) {   // R0 K2: K2 (r − r0)²
      if (n.size() < 2) return false;
      x.style = c2 ? "class2" : "harmonic";
      x.params = c2 ? std::vector<double>{n[0], n[1] * E, 0, 0} : std::vector<double>{n[1] * E, n[0]};
      return true;
    }, ff.auto_bonds);
    load(table("angle_auto"), 3, [&](FFRule& x, const std::vector<double>& n) {   // Theta0 K2: K2 (θ − θ0)²
      if (n.size() < 2) return false;
      x.style = c2 ? "quartic" : "harmonic";
      x.params = c2 ? std::vector<double>{n[0], n[1] * E, 0, 0} : std::vector<double>{n[1] * E, n[0]};
      return true;
    }, ff.auto_angles);
    load(table("dihedral_auto"), 4, [&](FFRule& x, const std::vector<double>& n) {   // Kphi n Phi0: Kphi [1 + cos(nφ − φ0)]
      if (n.size() < 3) return false;
      x.style = "fourier";
      x.params = {1, n[0] * E, n[1], n[2]};
      if (n[1] == 0) x.params = {0};
      return true;
    }, ff.auto_dihedrals);
    if (!ff.auto_bonds.empty() || !ff.auto_angles.empty() || !ff.auto_dihedrals.empty())
      ff.notes.push_back("automatic parameters from " + dir + prefix + "_*_auto: " + std::to_string(ff.auto_bonds.size()) + " bonds, " +
                         std::to_string(ff.auto_angles.size()) + " angles, " + std::to_string(ff.auto_dihedrals.size()) + " torsions");
  }
  for (const auto& sl : g_slips) ff.notes.push_back("library slip: torsion " + sl + "; taken as one torsion, as DL_FIELD does");
  ff.notes.push_back("converted from DL_FIELD by CAPS (" + pf.potential + ", unit " + (pf.unit.empty() ? "kcal/mol" : pf.unit) +
                     "); equivalences replace the type name for each kind of parameter");
  return ff;
}

}  // namespace caps
