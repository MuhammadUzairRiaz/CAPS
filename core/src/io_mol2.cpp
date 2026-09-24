// Tripos mol2 reader and writer (the preferred input for force-field assignment: explicit bonds, bond orders,
// atom types and charges, so nothing is guessed from distances).
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>

#include "caps/elements.hpp"
#include "caps/io.hpp"
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
  while (std::getline(in, line)) {
    ++lineno;
    const std::string t = trim(line);
    if (t.empty() || t[0] == '#') continue;
    if (t.rfind("@<TRIPOS>", 0) == 0) {
      section = t.substr(9);
      if (section == "MOLECULE") {
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
      if (rec_line == 5) charges = lower(t) != "no_charges";
    } else if (section == "ATOM") {
      if (w.size() < 6) throw ReadError(path + ":" + std::to_string(lineno) + ": ATOM record needs id name x y z type");
      Atom a;
      a.id = std::stoll(w[0]);
      a.pos = {std::stod(w[2]), std::stod(w[3]), std::stod(w[4])};
      a.name = w[5];
      a.element = element_of(w[5], w[1]);
      const int64_t sub = w.size() > 6 ? std::stoll(w[6]) : 1;
      a.mol = mol_offset + sub;
      max_sub = std::max(max_sub, sub);
      if (w.size() > 7) a.resname = w[7];
      if (w.size() > 8) a.charge = std::stod(w[8]);
      index[a.id] = uint32_t(s.atoms.size());
      a.id = int64_t(s.atoms.size()) + 1;
      s.atoms.push_back(a);
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
    }
  }
  (void)base;
  if (s.atoms.empty()) throw ReadError(path + ": no @<TRIPOS>ATOM records");
  s.has_charges = charges;
  s.has_mol = true;
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
  std::map<int, int> count;
  for (size_t i = 0; i < s.atoms.size(); ++i) {
    const Atom& a = s.atoms[i];
    const char* sym = a.element > 0 ? element(a.element).symbol : "X";
    const std::string name = std::string(sym) + std::to_string(++count[a.element]);
    const std::string type = a.name.empty() ? sym : a.name;
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
  std::fclose(f);
}

}  // namespace caps
