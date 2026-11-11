// Tripos mol2 reader and writer (the preferred input for force-field assignment: explicit bonds, bond orders,
// atom types and charges, so nothing is guessed from distances).
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <set>

#include "caps/elements.hpp"
#include "caps/io.hpp"
#include "caps/typing.hpp"
#include "io_util.hpp"

namespace caps {

namespace {

// Element from a SYBYL type ("C.ar", "N.am", "Cl") or, failing that, from the atom name ("C12", "HA3").
int element_of(const std::string& type, const std::string& name) {
  const auto dot = type.find('.');
  if (dot != std::string::npos) {
    if (int z = element_from_symbol(type.substr(0, dot))) return z;
  }
  if (type.size() <= 2 && !type.empty() && std::isupper(static_cast<unsigned char>(type[0])))
    if (int z = element_from_symbol(type)) return z;
  if (int z = element_from_name(name)) return z;
  return 0;
}

int order_code(const std::string& t) {
  if (t == "1") return 1;
  if (t == "2") return 2;
  if (t == "3") return 3;
  if (t == "ar") return 4;
  if (t == "am") return 5;
  return 0;   // "du", "un", "nc"
}

const char* order_text(int o) {
  switch (o) {
    case 2: return "2";
    case 3: return "3";
    case 4: return "ar";
    case 5: return "am";
    default: return "1";
  }
}

// The COMMENT line write_mol2 adds: the atom name column (2) holds Atom::name, a force-field type or an atom label,
// and the atom type column (6) the SYBYL type worked out for other programs. read_mol2 takes Atom::name back from
// column 2 when it sees it.
constexpr const char* kNameMarker = "CAPS: the atom name column holds each atom's CAPS type or label; the atom type column holds its SYBYL type";

bool known_element(int z) { return z > 0 && z <= max_element(); }

// Is t a Tripos SYBYL atom type for element z, in the exact spelling other readers expect?
bool is_sybyl_type(const std::string& t, int z) {
  if (!known_element(z)) return false;
  static const std::set<std::string> dotted = {"C.3",   "C.2",   "C.1",   "C.ar",  "C.cat", "N.3",   "N.2",   "N.1",   "N.ar",
                                               "N.am",  "N.pl3", "N.4",   "O.3",   "O.2",   "O.co2", "O.spc", "O.t3p", "S.3",
                                               "S.2",   "S.O",   "S.O2",  "P.3",   "H.spc", "H.t3p", "Cr.th", "Cr.oh", "Co.oh"};
  const std::string sym = element(z).symbol;
  const auto dot = t.find('.');
  if (dot != std::string::npos) return t.compare(0, dot, sym) == 0 && dotted.count(t) > 0;
  // C, N, O, S and P only come with a hybridisation suffix; the other elements are their bare symbol
  return t == sym && z != 6 && z != 7 && z != 8 && z != 15 && z != 16;
}

// SYBYL atom types from the structure: bond orders, aromaticity and formal charges as CAPS perceives them for typing
// (the file's own bond orders where it has them, otherwise from valences), a bare element symbol where SYBYL has
// nothing finer, "Du" for an unknown element.
std::vector<std::string> sybyl_types(const System& s) {
  const size_t n = s.atoms.size();
  std::vector<std::string> out(n);
  for (size_t i = 0; i < n; ++i) out[i] = known_element(s.atoms[i].element) ? element(s.atoms[i].element).symbol : "Du";
  bool bonds_ok = true;
  for (const auto& b : s.bonds) bonds_ok = bonds_ok && b.i < n && b.j < n && b.i != b.j;
  if (!bonds_ok) return out;
  Perception p;
  try {
    p = perceive(s);
  } catch (const std::exception&) {
    p = Perception{};
  }
  const bool have = p.nb.size() == n && p.order.size() == n && p.aromatic.size() == n && p.charge.size() == n;
  auto z_of = [&](uint32_t i) { return s.atoms[i].element; };
  // covalent neighbours (a coordinate bond to a metal carries order 0 and does not count)
  auto nbs = [&](uint32_t i) {
    std::vector<std::pair<uint32_t, int>> v;
    if (have)
      for (size_t k = 0; k < p.nb[i].size(); ++k)
        if (p.order[i][k] > 0) v.push_back({p.nb[i][k], p.order[i][k]});
    return v;
  };
  auto degree = [&](uint32_t i) { return nbs(i).size(); };
  auto max_order = [&](uint32_t i) {
    int m = 0;
    for (const auto& [o, ord] : nbs(i)) { (void)o; m = std::max(m, ord); }
    return m;
  };
  auto doubles = [&](uint32_t i) {
    int c = 0;
    for (const auto& [o, ord] : nbs(i)) { (void)o; c += ord == 2; }
    return c;
  };
  auto terminal_o = [&](uint32_t i) {
    int c = 0;
    for (const auto& [o, ord] : nbs(i)) { (void)ord; c += z_of(o) == 8 && degree(o) == 1; }
    return c;
  };
  // a carbonyl (or thiocarbonyl) carbon: C=O, C=S
  auto carbonyl_c = [&](uint32_t c) {
    if (z_of(c) != 6) return false;
    for (const auto& [o, ord] : nbs(c))
      if (ord == 2 && (z_of(o) == 8 || z_of(o) == 16)) return true;
    return false;
  };
  for (uint32_t i = 0; i < n; ++i) {
    const int z = z_of(i);
    if (!known_element(z) || !have) continue;
    const auto nb = nbs(i);
    const size_t d = nb.size();
    const int mo = max_order(i), dbl = doubles(i);
    if (z == 6) {
      if (p.aromatic[i]) out[i] = "C.ar";
      else if (mo == 3 || dbl >= 2) out[i] = "C.1";
      else if (dbl == 1) {
        // guanidinium / amidinium centre: three N, the positive charge on one of them
        int nn = 0, q = 0;
        for (const auto& [o, ord] : nb) { (void)ord; nn += z_of(o) == 7; q += z_of(o) == 7 ? p.charge[o] : 0; }
        out[i] = d == 3 && nn == 3 && q > 0 ? "C.cat" : "C.2";
      } else out[i] = d == 3 && p.charge[i] > 0 ? "C.cat" : "C.3";
    } else if (z == 7) {
      // a pyrrole-type ring N (three neighbours, no double bond in the Kekulé structure, neutral) is N.pl3, as SYBYL
      // types it; written as N.ar, readers cannot place its hydrogen or substituent and fail to kekulize the ring
      if (p.aromatic[i]) out[i] = d == 3 && dbl == 0 && p.charge[i] == 0 ? "N.pl3" : "N.ar";
      else if (mo == 3 || dbl >= 2) out[i] = "N.1";
      else if (d >= 4) out[i] = "N.4";
      else if (dbl == 1) out[i] = d == 3 ? "N.pl3" : "N.2";
      else if (p.charge[i] > 0) out[i] = "N.4";
      else {
        bool amide = false, conj = false;
        for (const auto& [o, ord] : nb) {
          (void)ord;
          if (carbonyl_c(o)) amide = true;
          else if (p.aromatic[o] || doubles(o) > 0) conj = true;
        }
        out[i] = d == 3 && amide ? "N.am" : (conj ? "N.pl3" : "N.3");
      }
    } else if (z == 8) {
      bool co2 = false;
      if (d == 1) {
        const uint32_t c = nb[0].first;
        // carboxylate (C with two terminal O and three neighbours) and phosphate O: SYBYL O.co2
        co2 = (z_of(c) == 6 && degree(c) == 3 && terminal_o(c) == 2) || (z_of(c) == 15 && terminal_o(c) >= 2);
        // a nitro, sulfonyl or sulfonate O is O.2 even where the Kekulé structure gives it a single bond and a charge
        if (!co2 && (z_of(c) == 7 || z_of(c) == 16) && terminal_o(c) >= 2) {
          out[i] = "O.2";
          continue;
        }
      }
      out[i] = co2 ? "O.co2" : (d == 1 && dbl == 1 ? "O.2" : "O.3");
    } else if (z == 16) {
      const int to = terminal_o(i);
      if (to == 1 && d == 3) out[i] = "S.O";
      else if (to >= 2) out[i] = "S.O2";
      else out[i] = dbl >= 1 && d <= 2 ? "S.2" : "S.3";
    } else if (z == 15) {
      out[i] = "P.3";
    }
  }
  return out;
}

// An atom name for column 2: Atom::name when it is one word, else nothing.
bool usable_name(const std::string& name) {
  if (name.empty()) return false;
  for (const char c : name)
    if (c == ' ' || c == '\t' || c == '\r' || c == '\n') return false;
  return true;
}

}  // namespace

System read_mol2(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw ReadError("cannot open " + path);
  System s;
  s.source_format = "mol2";
  std::string line, section;
  size_t base = 0;                 // first atom of the current MOLECULE record
  std::map<int64_t, uint32_t> index;   // atom id → index, per record
  int64_t mol_offset = 0, max_sub = 0;
  bool charges = false;
  int lineno = 0, rec_line = 0;
  // A file CAPS wrote keeps Atom::name in the atom name column (its COMMENT record says so); the names of the current
  // MOLECULE record wait here until the record ends.
  std::vector<std::string> labels;
  bool names_in_labels = false;
  auto end_record = [&]() {
    if (names_in_labels)
      for (size_t k = base; k < s.atoms.size(); ++k) s.atoms[k].name = labels[k - base];
    labels.clear();
    names_in_labels = false;
  };
  while (std::getline(in, line)) {
    ++lineno;
    const std::string t = trim(line);
    if (t.empty() || t[0] == '#') continue;
    if (t.rfind("@<TRIPOS>", 0) == 0) {
      section = t.substr(9);
      if (section == "MOLECULE") {
        end_record();
        base = s.atoms.size();
        index.clear();
        mol_offset += max_sub;
        max_sub = 0;
        rec_line = 0;
      }
      continue;
    }
    const auto w = split(t);
    if (section == "MOLECULE") {
      ++rec_line;
      if (rec_line == 1 && s.title.empty()) s.title = t;
      // name, counts, molecule type, charge type, then optional status bits and comment
      if (rec_line == 4) charges = lower(t) != "no_charges";
    } else if (section == "ATOM") {
      if (w.size() < 6) throw ReadError(path + ":" + std::to_string(lineno) + ": ATOM record needs id name x y z type");
      Atom a;
      a.id = std::stoll(w[0]);
      a.pos = {std::stod(w[2]), std::stod(w[3]), std::stod(w[4])};
      a.name = w[5];
      a.element = element_of(w[5], w[1]);
      const int64_t sub = w.size() > 6 ? std::stoll(w[6]) : 1;
      a.mol = mol_offset + sub;   // substructure (residue), kept per atom; molecules come from bonds
      a.resid = sub;
      max_sub = std::max(max_sub, sub);
      if (w.size() > 7) a.resname = w[7];
      if (w.size() > 8) a.charge = std::stod(w[8]);
      index[a.id] = uint32_t(s.atoms.size());
      a.id = int64_t(s.atoms.size()) + 1;
      s.atoms.push_back(a);
      labels.push_back(w[1]);
    } else if (section == "BOND") {
      if (w.size() < 3) continue;
      auto i = index.find(std::stoll(w[1])), j = index.find(std::stoll(w[2]));
      if (i == index.end() || j == index.end()) throw ReadError(path + ":" + std::to_string(lineno) + ": bond to an unknown atom");
      s.bonds.push_back({i->second, j->second, w.size() > 3 ? order_code(w[3]) : 0});
    } else if (section == "CRYSIN" && w.size() >= 6) {
      // a b c alpha beta gamma (Å, degrees): a along x, b in the xy plane
      const double a = std::stod(w[0]), b = std::stod(w[1]), c = std::stod(w[2]);
      const double al = std::stod(w[3]) * M_PI / 180, be = std::stod(w[4]) * M_PI / 180, ga = std::stod(w[5]) * M_PI / 180;
      s.cell.a = {a, 0, 0};
      s.cell.b = {b * std::cos(ga), b * std::sin(ga), 0};
      const double cx = c * std::cos(be), cy = c * (std::cos(al) - std::cos(be) * std::cos(ga)) / std::sin(ga);
      s.cell.c = {cx, cy, std::sqrt(std::max(0.0, c * c - cx * cx - cy * cy))};
    } else if (section == "COMMENT" && t == kNameMarker) {
      names_in_labels = true;
    }
  }
  end_record();
  if (s.atoms.empty()) throw ReadError(path + ": no @<TRIPOS>ATOM records");
  s.has_charges = charges;
  s.has_mol = false;   // substructures are residues, not molecules
  s.bonds_from_file = !s.bonds.empty();
  if (!s.cell.valid()) s.cell.periodic = {false, false, false};
  return s;
}

void write_mol2(const System& s, const std::string& path) {
  std::FILE* f = std::fopen(path.c_str(), "w");
  if (!f) throw ReadError("cannot write " + path);
  std::fprintf(f, "@<TRIPOS>MOLECULE\n%s\n%5zu %5zu %5d 0 0\n%s\n%s\n\n", s.title.empty() ? "CAPS" : s.title.c_str(), s.atoms.size(),
               s.bonds.size(), 1, s.atoms.size() > 1000 ? "BIOPOLYMER" : "SMALL", s.has_charges ? "USER_CHARGES" : "NO_CHARGES");
  std::fprintf(f, "@<TRIPOS>ATOM\n");
  // Column 6 is the SYBYL type other programs take the element from; Atom::name (a force-field type such as GAFF
  // "c3", or a label such as "C12" from a builder, a PDB or a CIF) goes in the atom name column, so a file read back
  // by CAPS keeps it (the COMMENT record below says so). An atom with no name gets its element symbol there, which
  // is what CAPS read back for it before.
  const std::vector<std::string> sybyl = sybyl_types(s);
  for (size_t i = 0; i < s.atoms.size(); ++i) {
    const Atom& a = s.atoms[i];
    const char* sym = known_element(a.element) ? element(a.element).symbol : "X";
    const std::string name = usable_name(a.name) ? a.name : std::string(sym);
    const std::string& type = is_sybyl_type(a.name, a.element) ? a.name : sybyl[i];
    std::fprintf(f, "%7zu %-8s %11.5f %11.5f %11.5f %-8s %5lld %-6s %10.6f\n", i + 1, name.c_str(), a.pos[0], a.pos[1], a.pos[2], type.c_str(),
                 static_cast<long long>(a.mol > 0 ? a.mol : 1), a.resname.empty() ? "MOL" : a.resname.c_str(), a.charge);
  }
  std::fprintf(f, "@<TRIPOS>BOND\n");
  for (size_t k = 0; k < s.bonds.size(); ++k)
    std::fprintf(f, "%6zu %5u %5u %s\n", k + 1, s.bonds[k].i + 1, s.bonds[k].j + 1, order_text(s.bonds[k].order));
  if (s.cell.valid()) {
    const double a = norm(s.cell.a), b = norm(s.cell.b), c = norm(s.cell.c);
    auto ang = [](const Vec3& u, const Vec3& v) { return std::acos(dot(u, v) / (norm(u) * norm(v))) * 180 / M_PI; };
    std::fprintf(f, "@<TRIPOS>CRYSIN\n%10.4f %10.4f %10.4f %8.3f %8.3f %8.3f 1 1\n", a, b, c, ang(s.cell.b, s.cell.c), ang(s.cell.a, s.cell.c),
                 ang(s.cell.a, s.cell.b));
  }
  std::fprintf(f, "@<TRIPOS>COMMENT\n%s\n", kNameMarker);
  std::fclose(f);
}

}  // namespace caps
