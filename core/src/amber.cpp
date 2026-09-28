// AMBER prmtop and inpcrd / rst7 readers (see amber.hpp). The prmtop is a list of %FLAG sections, each with a Fortran
// %FORMAT (20a4, 10I8, 5E16.8 …) read by its fixed field widths. Units as the file: charges × 18.2223 (e·√(kcal Å/mol)),
// bond k in kcal/(mol Å²) for k (r − r0)², angle k in kcal/(mol rad²) for k (θ − θ0)², torsions V [1 + cos(nφ − δ)],
// Lennard-Jones A/r¹² − B/r⁶ per type pair.
#include "caps/amber.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <map>
#include <set>
#include <sstream>

#include "caps/elements.hpp"
#include "caps/io.hpp"

namespace caps {

namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr double kChargeUnit = 18.2223;          // AMBER's internal charge unit per e (√332.0522173)
constexpr double kVelUnit = 20.455;              // AMBER velocity unit: Å per (1/20.455 ps)

std::string lower(std::string s) {
  for (auto& c : s) c = char(std::tolower(static_cast<unsigned char>(c)));
  return s;
}

std::string trim(const std::string& s) {
  const auto a = s.find_first_not_of(" \t\r\n");
  if (a == std::string::npos) return "";
  return s.substr(a, s.find_last_not_of(" \t\r\n") - a + 1);
}

struct Section {
  char kind = 'a';   // a text, I integer, E / F real
  int width = 80;
  std::vector<std::string> fields;
};

// %FORMAT(10I8) → I, 8; (5E16.8) → E, 16; (20a4) → a, 4; (1a80) → a, 80
void parse_format(const std::string& f, Section& s, const std::string& where) {
  const auto a = f.find('('), b = f.find(')');
  if (a == std::string::npos || b == std::string::npos) throw ReadError(where + ": bad %FORMAT " + f);
  const std::string in = f.substr(a + 1, b - a - 1);
  size_t k = 0;
  while (k < in.size() && std::isdigit(static_cast<unsigned char>(in[k]))) ++k;
  if (k >= in.size()) throw ReadError(where + ": bad %FORMAT " + f);
  s.kind = char(std::toupper(static_cast<unsigned char>(in[k])));
  if (s.kind == 'A') s.kind = 'a';
  size_t e = k + 1;
  while (e < in.size() && std::isdigit(static_cast<unsigned char>(in[e]))) ++e;
  s.width = std::stoi(in.substr(k + 1, e - k - 1));
  if (s.width <= 0) throw ReadError(where + ": bad %FORMAT " + f);
}

std::map<std::string, Section> read_sections(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw ReadError("cannot open " + path);
  std::map<std::string, Section> out;
  Section* cur = nullptr;
  std::string line, flag;
  size_t ln = 0;
  bool any = false;
  while (std::getline(in, line)) {
    ++ln;
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (line.rfind("%VERSION", 0) == 0 || line.rfind("%COMMENT", 0) == 0) continue;
    if (line.rfind("%FLAG", 0) == 0) {
      flag = trim(line.substr(5));
      cur = &out[flag];
      *cur = Section{};
      any = true;
      continue;
    }
    if (line.rfind("%FORMAT", 0) == 0) {
      if (!cur) throw ReadError(path + ":" + std::to_string(ln) + ": %FORMAT before any %FLAG");
      parse_format(line, *cur, path + ":" + std::to_string(ln));
      continue;
    }
    if (!cur) continue;
    // fixed-width fields; a text field keeps its spaces trimmed, a blank numeric tail ends the line
    for (size_t p = 0; p < line.size(); p += size_t(cur->width)) {
      const std::string f = line.substr(p, size_t(cur->width));
      if (cur->kind == 'a') cur->fields.push_back(trim(f));
      else {
        const std::string t = trim(f);
        if (t.empty()) break;
        cur->fields.push_back(t);
      }
    }
    if (cur->kind == 'a' && cur->width >= 80) {   // a title line: one field, whatever its length
      if (!cur->fields.empty()) cur->fields.back() = trim(line);
    }
  }
  if (!any) throw ReadError(path + ": not an AMBER topology (no %FLAG sections)");
  return out;
}

struct Reader {
  const std::map<std::string, Section>& S;
  std::string path;
  bool has(const std::string& f) const { return S.count(f) > 0; }
  const Section& sec(const std::string& f) const {
    auto it = S.find(f);
    if (it == S.end()) throw ReadError(path + ": no %FLAG " + f);
    return it->second;
  }
  std::vector<long> ints(const std::string& f, size_t want) const {
    const auto& s = sec(f);
    std::vector<long> v;
    v.reserve(s.fields.size());
    for (const auto& x : s.fields) {
      try { v.push_back(std::stol(x)); } catch (...) { throw ReadError(path + ": %FLAG " + f + ": '" + x + "' is not an integer"); }
    }
    if (want != size_t(-1) && v.size() < want)
      throw ReadError(path + ": %FLAG " + f + " has " + std::to_string(v.size()) + " values, " + std::to_string(want) + " expected");
    return v;
  }
  std::vector<double> reals(const std::string& f, size_t want) const {
    const auto& s = sec(f);
    std::vector<double> v;
    v.reserve(s.fields.size());
    for (auto x : s.fields) {
      std::replace(x.begin(), x.end(), 'D', 'E');   // Fortran double exponents
      try { v.push_back(std::stod(x)); } catch (...) { throw ReadError(path + ": %FLAG " + f + ": '" + x + "' is not a number"); }
    }
    if (want != size_t(-1) && v.size() < want)
      throw ReadError(path + ": %FLAG " + f + " has " + std::to_string(v.size()) + " values, " + std::to_string(want) + " expected");
    return v;
  }
  std::vector<std::string> texts(const std::string& f, size_t want) const {
    const auto& s = sec(f);
    if (want != size_t(-1) && s.fields.size() < want)
      throw ReadError(path + ": %FLAG " + f + " has " + std::to_string(s.fields.size()) + " entries, " + std::to_string(want) + " expected");
    return s.fields;
  }
};

Cell box_cell(double a, double b, double c, double al, double be, double ga) {
  const double ca = std::cos(al * kPi / 180), cb = std::cos(be * kPi / 180), cg = std::cos(ga * kPi / 180), sg = std::sin(ga * kPi / 180);
  Cell cell;
  cell.a = {a, 0, 0};
  cell.b = {b * cg, b * sg, 0};
  const double cx = c * cb, cy = c * (ca - cb * cg) / sg;
  cell.c = {cx, cy, std::sqrt(std::max(0.0, c * c - cx * cx - cy * cy))};
  return cell;
}

}  // namespace

bool is_amber_topology_path(const std::string& path) {
  std::string p = lower(path);
  if (p.size() > 3 && p.substr(p.size() - 3) == ".gz") p = p.substr(0, p.size() - 3);
  const std::string e = std::filesystem::path(p).extension().string();
  return e == ".prmtop" || e == ".parm7" || e == ".prmtp";
}

bool is_amber_coordinates_path(const std::string& path) {
  const std::string e = lower(std::filesystem::path(path).extension().string());
  return e == ".inpcrd" || e == ".rst7" || e == ".restrt" || e == ".crd" || e == ".rst" || e == ".ncrst" || e == ".mdcrd_rst";
}

AmberTopology read_amber_prmtop(const std::string& path) {
  const auto S = read_sections(path);
  const Reader R{S, path};
  const std::string name = std::filesystem::path(path).filename().string();
  auto refuse = [&](const std::string& what) { throw ReadError(path + ": " + what); };
  if (R.has("CTITLE") || R.has("FORCE_FIELD_TYPE") || R.has("CHARMM_UREY_BRADLEY"))
    refuse("a CHARMM topology converted by chamber (Urey–Bradley, harmonic impropers, CHARMM 1-4 parameters): open the CHARMM files (PSF + parameters) instead");
  for (const char* f : {"CHARMM_CMAP_COUNT", "CMAP_COUNT"})
    if (R.has(f)) refuse("CMAP correction maps (ff19SB, CHARMM) are not read: CAPS has no CMAP term, and leaving them out would change the energy");
  if (R.has("POLARIZABILITY") || R.has("DIPOLE_DAMP_FACTOR")) refuse("a polarisable topology (AMOEBA / induced dipoles) is not read");

  const auto P = R.ints("POINTERS", 30);
  const size_t natom = size_t(P[0]), ntypes = size_t(P[1]), nbonh = size_t(P[2]), mbona = size_t(P[3]), ntheth = size_t(P[4]),
               mtheta = size_t(P[5]), nphih = size_t(P[6]), mphia = size_t(P[7]), nnb = size_t(P[10]), nres = size_t(P[11]),
               numbnd = size_t(P[15]), numang = size_t(P[16]), nptra = size_t(P[17]), nphb = size_t(P[19]);
  const long ifpert = P[20], ifbox = P[27], ifcap = P.size() > 29 ? P[29] : 0, numextra = P.size() > 30 ? P[30] : 0;
  if (ifpert) refuse("a perturbation topology (IFPERT) is not read");
  if (ifcap) refuse("a solvent-cap topology (IFCAP) is not read");
  if (numextra > 0) refuse(std::to_string(numextra) + " extra points (TIP4P / TIP5P virtual sites): not read");
  if (natom == 0) refuse("no atoms");

  const auto names = R.texts("ATOM_NAME", natom);
  auto charge = R.reals("CHARGE", natom);
  const auto mass = R.reals("MASS", natom);
  const auto tindex = R.ints("ATOM_TYPE_INDEX", natom);
  const auto nexcl = R.ints("NUMBER_EXCLUDED_ATOMS", natom);
  const auto nbindex = R.ints("NONBONDED_PARM_INDEX", ntypes * ntypes);
  const auto reslab = R.texts("RESIDUE_LABEL", nres);
  const auto resptr = R.ints("RESIDUE_POINTER", nres);
  const auto bk = R.reals("BOND_FORCE_CONSTANT", numbnd), br = R.reals("BOND_EQUIL_VALUE", numbnd);
  const auto ak = R.reals("ANGLE_FORCE_CONSTANT", numang), at = R.reals("ANGLE_EQUIL_VALUE", numang);
  const auto dk = R.reals("DIHEDRAL_FORCE_CONSTANT", nptra), dn = R.reals("DIHEDRAL_PERIODICITY", nptra), dp = R.reals("DIHEDRAL_PHASE", nptra);
  std::vector<double> scee(nptra, 1.2), scnb(nptra, 2.0);
  const bool has_scale = R.has("SCEE_SCALE_FACTOR") && R.has("SCNB_SCALE_FACTOR");
  if (has_scale) scee = R.reals("SCEE_SCALE_FACTOR", nptra), scnb = R.reals("SCNB_SCALE_FACTOR", nptra);
  const size_t npair = ntypes * (ntypes + 1) / 2;
  const auto lja = R.reals("LENNARD_JONES_ACOEF", npair), ljb = R.reals("LENNARD_JONES_BCOEF", npair);
  if (nphb > 0 && R.has("HBOND_ACOEF")) {
    const auto ha = R.reals("HBOND_ACOEF", 0), hb = R.reals("HBOND_BCOEF", 0);
    for (size_t k = 0; k < std::min(ha.size(), hb.size()); ++k)
      if (ha[k] != 0 || hb[k] != 0) refuse("10-12 hydrogen-bond terms with non-zero coefficients are not read");
  }
  const auto amber_type = R.has("AMBER_ATOM_TYPE") ? R.texts("AMBER_ATOM_TYPE", natom) : names;
  const auto excl = R.ints("EXCLUDED_ATOMS_LIST", nnb);
  std::vector<long> atomic_number;
  if (R.has("ATOMIC_NUMBER")) atomic_number = R.ints("ATOMIC_NUMBER", natom);

  AmberTopology out;
  auto& notes = out.notes;
  System& s = out.system;
  s.source_format = "amber-prmtop";
  if (R.has("TITLE") && !R.sec("TITLE").fields.empty()) s.title = R.sec("TITLE").fields.front();
  else if (R.has("CTITLE") && !R.sec("CTITLE").fields.empty()) s.title = R.sec("CTITLE").fields.front();

  auto ff = std::make_shared<ForceField>();
  ForceField& F = *ff;
  F.name = "AMBER topology " + name;
  F.pair_form = "lj12-6";
  F.mixing = "arithmetic";
  F.native_pair = "lj/cut/coul/long";
  F.native_dihedral = "fourier";
  F.native_improper = "cvff";
  F.native_special = "amber";
  F.native_cutoff = 10;

  // types: an AMBER atom type with its Lennard-Jones index and mass (hydrogen mass repartitioning gives one type two masses)
  std::map<std::tuple<std::string, long, double>, int> tkey;
  std::map<std::string, int> name_uses;
  std::vector<long> type_lj;   // per CAPS type: the file's LJ type (0-based)
  F.atom_type.resize(natom);
  F.type_index.resize(natom);
  F.charge.resize(natom);
  F.mass.resize(natom);
  F.why.assign(natom, "type from the AMBER topology");
  for (size_t i = 0; i < natom; ++i) {
    const long lt = tindex[i] - 1;
    if (lt < 0 || size_t(lt) >= ntypes) refuse("atom " + std::to_string(i + 1) + " has Lennard-Jones type " + std::to_string(tindex[i]) + " of " + std::to_string(ntypes));
    const auto key = std::make_tuple(amber_type[i], lt, mass[i]);
    auto it = tkey.find(key);
    if (it == tkey.end()) {
      const int uses = ++name_uses[amber_type[i]];
      const std::string tn = uses == 1 ? amber_type[i] : amber_type[i] + "_" + std::to_string(uses);
      it = tkey.emplace(key, int(F.type_names.size())).first;
      F.type_names.push_back(tn);
      type_lj.push_back(lt);
    }
    F.type_index[i] = it->second;
    F.atom_type[i] = F.type_names[size_t(it->second)];
    F.charge[i] = charge[i] / kChargeUnit;
    F.mass[i] = mass[i];
  }
  for (const auto& [n, u] : name_uses)
    if (u > 1) notes.push_back("AMBER type " + n + " comes with " + std::to_string(u) + " masses or Lennard-Jones types: kept apart as " + n + ", " + n + "_2 …");

  // Lennard-Jones: A = 4εσ¹², B = 4εσ⁶ for each type with itself; a pair the arithmetic (Lorentz–Berthelot) rule does
  // not reproduce (NBFIX, off-diagonal edits) keeps its own coefficients
  auto ab = [&](long a, long b, double& A, double& B) {
    const long ix = nbindex[size_t(a) * ntypes + size_t(b)];
    if (ix < 0) { A = B = 0; return false; }   // a 10-12 pair (checked zero above)
    if (ix == 0 || size_t(ix) > npair) refuse("NONBONDED_PARM_INDEX " + std::to_string(ix) + " out of range");
    A = lja[size_t(ix - 1)], B = ljb[size_t(ix - 1)];
    return true;
  };
  auto eps_sigma = [](double A, double B) -> PairType {
    if (A <= 0 || B <= 0) return {0, 0};
    return {B * B / (4 * A), std::pow(A / B, 1.0 / 6)};
  };
  const size_t nt = F.type_names.size();
  F.lj.resize(nt);
  for (size_t t = 0; t < nt; ++t) {
    double A, B;
    ab(type_lj[t], type_lj[t], A, B);
    if ((A > 0) != (B > 0)) refuse("type " + F.type_names[t] + ": a purely repulsive or purely attractive Lennard-Jones pair (A " + std::to_string(A) + ", B " + std::to_string(B) + ") has no ε, σ");
    F.lj[t] = eps_sigma(A, B);
  }
  int overrides = 0;
  for (size_t a = 0; a < nt; ++a)
    for (size_t b = a + 1; b < nt; ++b) {
      double A, B;
      ab(type_lj[a], type_lj[b], A, B);
      const PairType &p = F.lj[a], &q = F.lj[b];
      const double e = std::sqrt(p.eps * q.eps), sg = 0.5 * (p.sigma + q.sigma);
      const double Am = 4 * e * std::pow(sg, 12), Bm = 4 * e * std::pow(sg, 6);
      auto close = [](double x, double y) { return std::fabs(x - y) <= 2e-6 * std::max({std::fabs(x), std::fabs(y), 1e-30}); };
      if (close(A, Am) && close(B, Bm)) continue;
      if ((A > 0) != (B > 0)) refuse("types " + F.type_names[a] + " and " + F.type_names[b] + ": a Lennard-Jones pair with only A or only B has no ε, σ");
      F.pair_override[{int(a), int(b)}] = A > 0 ? eps_sigma(A, B) : PairType{0, 1};
      ++overrides;
    }
  if (overrides) notes.push_back(std::to_string(overrides) + " type pairs with their own Lennard-Jones coefficients (not the arithmetic rule)");

  // atoms: residues, molecules
  s.atoms.resize(natom);
  std::vector<long> resstart(resptr.begin(), resptr.begin() + long(nres));
  for (size_t r = 0; r < nres; ++r) {
    const size_t a0 = size_t(resstart[r] - 1), a1 = r + 1 < nres ? size_t(resstart[r + 1] - 1) : natom;
    if (a0 > a1 || a1 > natom) refuse("RESIDUE_POINTER out of order");
    for (size_t i = a0; i < a1; ++i) s.atoms[i].resname = reslab[r], s.atoms[i].resid = long(r + 1);
  }
  std::vector<long> mol(natom, 0);
  if (ifbox > 0 && R.has("ATOMS_PER_MOLECULE")) {
    const auto sp = R.ints("SOLVENT_POINTERS", 3);
    const auto apm = R.ints("ATOMS_PER_MOLECULE", size_t(sp[1]));
    size_t i = 0;
    for (size_t m = 0; m < size_t(sp[1]); ++m)
      for (long k = 0; k < apm[m] && i < natom; ++k) mol[i++] = long(m + 1);
    if (i != natom) refuse("ATOMS_PER_MOLECULE covers " + std::to_string(i) + " of " + std::to_string(natom) + " atoms");
  }

  // bonded terms (coordinate indices: 3 × atom)
  const size_t bad = size_t(-1);
  auto atom = [&](long c) -> uint32_t {
    const long a = std::labs(c);
    if (a % 3 != 0 || size_t(a / 3) >= natom) return uint32_t(bad);
    return uint32_t(a / 3);
  };
  auto check = [&](uint32_t a, const char* what) {
    if (a == uint32_t(bad)) refuse(std::string("a ") + what + " refers to an atom out of range");
    return a;
  };
  std::set<std::pair<uint32_t, uint32_t>> bonded;
  auto bonds = [&](const char* flag, size_t n) {
    const auto v = R.ints(flag, 3 * n);
    for (size_t k = 0; k < n; ++k) {
      const uint32_t i = check(atom(v[3 * k]), "bond"), j = check(atom(v[3 * k + 1]), "bond");
      const long t = v[3 * k + 2] - 1;
      if (t < 0 || size_t(t) >= numbnd) refuse("bond type out of range");
      F.bonds.push_back({i, j, bk[size_t(t)], br[size_t(t)]});
      s.bonds.push_back({i, j, 1});
      bonded.insert({std::min(i, j), std::max(i, j)});
    }
  };
  bonds("BONDS_INC_HYDROGEN", nbonh);
  bonds("BONDS_WITHOUT_HYDROGEN", mbona);
  auto angles = [&](const char* flag, size_t n) {
    const auto v = R.ints(flag, 4 * n);
    for (size_t k = 0; k < n; ++k) {
      const long t = v[4 * k + 3] - 1;
      if (t < 0 || size_t(t) >= numang) refuse("angle type out of range");
      F.angles.push_back({check(atom(v[4 * k]), "angle"), check(atom(v[4 * k + 1]), "angle"), check(atom(v[4 * k + 2]), "angle"), ak[size_t(t)], at[size_t(t)]});
    }
  };
  angles("ANGLES_INC_HYDROGEN", ntheth);
  angles("ANGLES_WITHOUT_HYDROGEN", mtheta);
  // torsions: a negative third index — the 1-4 pair is not computed here (another term or a ring already counts it); a
  // negative fourth — an improper. AMBER impropers are periodic torsions too; one whose phase is not 0 or 180° is kept
  // with the proper torsions (the same energy, and an exact LAMMPS form)
  std::set<std::pair<uint32_t, uint32_t>> p14;
  std::set<std::pair<double, double>> scales;   // (SCEE, SCNB) of the torsions that carry a 1-4 pair
  int improper_as_proper = 0;
  auto torsions = [&](const char* flag, size_t n) {
    const auto v = R.ints(flag, 5 * n);
    for (size_t k = 0; k < n; ++k) {
      const uint32_t i = check(atom(v[5 * k]), "torsion"), j = check(atom(v[5 * k + 1]), "torsion"), c = check(atom(v[5 * k + 2]), "torsion"),
                     l = check(atom(v[5 * k + 3]), "torsion");
      const long t = v[5 * k + 4] - 1;
      if (t < 0 || size_t(t) >= nptra) refuse("torsion type out of range");
      const bool no14 = v[5 * k + 2] < 0, improper = v[5 * k + 3] < 0;
      const double per = dn[size_t(t)];
      if (std::fabs(per - std::round(per)) > 1e-6) refuse("a torsion with periodicity " + std::to_string(per) + " (not whole) is not read");
      const TorsionTerm term{i, j, c, l, dk[size_t(t)], int(std::lround(std::fabs(per))), dp[size_t(t)]};
      if (improper) {
        const double cs = std::cos(term.delta);
        if (std::fabs(std::fabs(cs) - 1) < 1e-9) F.impropers.push_back(term);
        else { F.dihedrals.push_back(term); ++improper_as_proper; }
      } else {
        F.dihedrals.push_back(term);
        if (!no14) {
          p14.insert({std::min(i, l), std::max(i, l)});
          scales.insert({scee[size_t(t)], scnb[size_t(t)]});
        }
      }
    }
  };
  torsions("DIHEDRALS_INC_HYDROGEN", nphih);
  torsions("DIHEDRALS_WITHOUT_HYDROGEN", mphia);
  if (improper_as_proper) notes.push_back(std::to_string(improper_as_proper) + " improper torsions with a phase other than 0 or 180° kept with the proper torsions (the same periodic form)");
  // zero-barrier torsion terms carry only their 1-4 pair: drop the term, keep the pair
  F.dihedrals.erase(std::remove_if(F.dihedrals.begin(), F.dihedrals.end(), [](const TorsionTerm& t) { return t.v == 0; }), F.dihedrals.end());
  if (scales.size() > 1) {
    std::string list;
    for (const auto& [e, v] : scales) list += (list.empty() ? "" : ", ") + std::to_string(e) + "/" + std::to_string(v);
    refuse("1-4 scaling differs between torsions (SCEE/SCNB " + list + ", e.g. GLYCAM with a protein force field): CAPS applies one 1-4 scaling to the whole structure");
  }
  const double sc_e = scales.empty() ? 1.2 : scales.begin()->first, sc_n = scales.empty() ? 2.0 : scales.begin()->second;
  if (sc_e <= 0 || sc_n <= 0) refuse("a 1-4 scale factor of zero");
  F.coul14 = 1 / sc_e;
  F.lj14 = 1 / sc_n;
  if (!has_scale) notes.push_back("no SCEE / SCNB sections (an older topology): the AMBER defaults 1.2 and 2.0");
  for (const auto& [a, b] : p14) F.pairs14.push_back({a, b});

  // exclusions: the file's list (1-2, 1-3 and 1-4 partners; a lone 0 marks an atom with none), made symmetric
  F.excluded.assign(natom, {});
  size_t pos = 0;
  for (size_t i = 0; i < natom; ++i) {
    const long ne = nexcl[i];
    if (ne < 0 || pos + size_t(ne) > excl.size()) refuse("EXCLUDED_ATOMS_LIST is shorter than NUMBER_EXCLUDED_ATOMS says");
    for (long k = 0; k < ne; ++k) {
      const long j = excl[pos + size_t(k)];
      if (j <= 0) continue;
      if (size_t(j) > natom) refuse("an excluded atom out of range");
      F.excluded[i].push_back(uint32_t(j - 1));
      F.excluded[size_t(j - 1)].push_back(uint32_t(i));
    }
    pos += size_t(ne);
  }
  auto add_ex = [&](uint32_t a, uint32_t b) { F.excluded[a].push_back(b), F.excluded[b].push_back(a); };
  int added = 0;
  auto excluded = [&](uint32_t a, uint32_t b) { return std::find(F.excluded[a].begin(), F.excluded[a].end(), b) != F.excluded[a].end(); };
  for (const auto& [a, b] : bonded) if (!excluded(a, b)) add_ex(a, b), ++added;
  for (const auto& pr : F.pairs14) if (!excluded(pr[0], pr[1])) add_ex(pr[0], pr[1]), ++added;
  for (auto& e : F.excluded) {
    std::sort(e.begin(), e.end());
    e.erase(std::unique(e.begin(), e.end()), e.end());
  }
  if (added) notes.push_back(std::to_string(added) + " bonded or 1-4 pairs missing from EXCLUDED_ATOMS_LIST added to the exclusions");

  // atoms of the structure: named by their types (as a GROMACS topology), elements from the file or the mass
  int from_mass = 0;
  for (size_t i = 0; i < natom; ++i) {
    Atom& a = s.atoms[i];
    a.id = int64_t(i + 1);
    a.mol = mol[i];
    a.type = F.type_index[i] + 1;
    a.name = F.atom_type[i];
    a.charge = F.charge[i];
    int z = atomic_number.empty() ? 0 : int(atomic_number[i]);
    if (z <= 0 || z > max_element()) {
      // the atom name's element when its mass fits it (hydrogens repartitioned to ~3 u stay H), else the nearest mass
      const int zn = element_from_name(names[i]);
      if (zn > 0 && (std::fabs(element(zn).mass - mass[i]) < 0.6 || (zn == 1 && mass[i] < 4.1))) z = zn;
      else z = element_from_mass(mass[i], 0.6);
      ++from_mass;
    }
    a.element = z;
  }
  if (from_mass && atomic_number.empty()) notes.push_back("no ATOMIC_NUMBER section: elements from the atom names and masses");
  s.types.clear();
  for (size_t t = 0; t < nt; ++t) {
    TypeInfo ti;
    ti.type = int(t + 1);
    ti.label = F.type_names[t];
    for (size_t i = 0; i < natom; ++i)
      if (size_t(F.type_index[i]) == t) { ti.mass = mass[i]; break; }
    s.types.push_back(ti);
  }
  s.bonds_from_file = true;
  s.has_charges = true;
  if (ifbox > 0) {
    if (R.has("BOX_DIMENSIONS")) {
      const auto b = R.reals("BOX_DIMENSIONS", 4);
      s.cell = box_cell(b[1], b[2], b[3], ifbox == 2 ? 109.4712206 : 90, b[0], ifbox == 2 ? 109.4712206 : 90);
    }
    if (ifbox == 2) notes.push_back("a truncated octahedron (IFBOX 2): the box as its triclinic cell");
  }
  if (mol.front() == 0) {   // no molecule list: bonded fragments
    int count = 0;
    const auto m = s.molecules(&count);
    for (size_t i = 0; i < natom; ++i) s.atoms[i].mol = m[i] + 1;
  }
  s.has_mol = true;

  // net charge as the file holds it (AMBER charges are × 18.2223, so the sum is rarely exactly whole)
  double q = 0;
  for (double x : F.charge) q += x;
  F.notes.push_back(std::string("from ") + name + ": " + std::to_string(natom) + " atoms, " + std::to_string(nt) + " types, " + std::to_string(F.bonds.size()) +
                    " bonds, " + std::to_string(F.angles.size()) + " angles, " + std::to_string(F.dihedrals.size()) + " torsion terms, " +
                    std::to_string(F.impropers.size()) + " impropers, " + std::to_string(F.pairs14.size()) + " 1-4 pairs (Coulomb × 1/" +
                    std::to_string(sc_e).substr(0, 4) + ", LJ × 1/" + std::to_string(sc_n).substr(0, 4) + ")");
  char buf[64];
  std::snprintf(buf, sizeof buf, "%.4f", q);
  F.notes.push_back(std::string("net charge ") + buf + " e");
  for (const auto& n : notes) F.notes.push_back(n);
  s.notes = notes;
  out.ff = ff;
  return out;
}

AmberCoordinates read_amber_coordinates(const std::string& path) {
  std::ifstream in(path, std::ios::binary);
  if (!in) throw ReadError("cannot open " + path);
  char magic[4] = {0, 0, 0, 0};
  in.read(magic, 3);
  if (std::string(magic, 3) == "CDF" || (magic[0] == '\x89' && magic[1] == 'H'))
    throw ReadError(path + ": a NetCDF restart — convert it to an ASCII restart (cpptraj: trajout name.rst7 restart) to open it");
  in.clear();
  in.seekg(0);
  AmberCoordinates c;
  std::string line;
  std::getline(in, line);
  if (!line.empty() && line.back() == '\r') line.pop_back();
  c.title = trim(line);
  if (!std::getline(in, line)) throw ReadError(path + ": no atom count line");
  size_t natom = 0;
  {
    std::istringstream is(line);
    long n = 0;
    if (!(is >> n) || n <= 0) throw ReadError(path + ": no atom count on line 2");
    natom = size_t(n);
  }
  // 6F12.7 records; a line may be shorter at the end of a block
  std::vector<double> v;
  v.reserve(natom * 3 + 6);
  size_t ln = 2;
  while (std::getline(in, line)) {
    ++ln;
    if (!line.empty() && line.back() == '\r') line.pop_back();
    for (size_t p = 0; p + 1 <= line.size(); p += 12) {
      const std::string f = trim(line.substr(p, 12));
      if (f.empty()) continue;
      try { v.push_back(std::stod(f)); } catch (...) { throw ReadError(path + ":" + std::to_string(ln) + ": '" + f + "' is not a number"); }
    }
  }
  const size_t n3 = 3 * natom;
  if (v.size() < n3) throw ReadError(path + ": " + std::to_string(v.size() / 3) + " coordinates for " + std::to_string(natom) + " atoms");
  const bool vel = v.size() >= 2 * n3, box = v.size() == n3 + 6 || v.size() == 2 * n3 + 6;
  if (!(v.size() == n3 || v.size() == n3 + 6 || v.size() == 2 * n3 || v.size() == 2 * n3 + 6))
    throw ReadError(path + ": " + std::to_string(v.size()) + " numbers after the header fit neither coordinates, velocities nor a box for " + std::to_string(natom) + " atoms");
  c.positions.resize(natom);
  for (size_t i = 0; i < natom; ++i) c.positions[i] = {v[3 * i], v[3 * i + 1], v[3 * i + 2]};
  if (vel) {
    c.velocities.resize(natom);
    for (size_t i = 0; i < natom; ++i) c.velocities[i] = Vec3{v[n3 + 3 * i], v[n3 + 3 * i + 1], v[n3 + 3 * i + 2]} * (kVelUnit / 1000.0);
  }
  if (box) {
    const double* b = &v[v.size() - 6];
    c.cell = box_cell(b[0], b[1], b[2], b[3], b[4], b[5]);
    c.has_box = true;
  }
  return c;
}

}  // namespace caps
