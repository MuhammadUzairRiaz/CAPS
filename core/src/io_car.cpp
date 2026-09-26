// Materials Studio / Discover structures: the .car file (coordinates, cell, force-field type, element and charge per atom,
// one block per molecule) and, beside it, the .mdf file (the bonds, including those across the cell). INTERFACE, ClayFF,
// PCFF and COMPASS models are distributed this way, typed and charged.
#include <cmath>
#include <filesystem>
#include <fstream>
#include <map>
#include <set>
#include <cstdlib>
#include <algorithm>
#include <ctime>
#include <stdexcept>

#include "caps/elements.hpp"
#include "caps/io.hpp"
#include "io_util.hpp"

namespace caps {

namespace {

constexpr double kDeg = 3.14159265358979323846 / 180.0;

// A cell from a b c α β γ, as Materials Studio orients it: a along x, b in the xy plane.
Cell cell_from(double a, double b, double c, double al, double be, double ga) {
  Cell k;
  const double ca = std::cos(al * kDeg), cb = std::cos(be * kDeg), cg = std::cos(ga * kDeg), sg = std::sin(ga * kDeg);
  k.a = {a, 0, 0};
  k.b = {b * cg, b * sg, 0};
  const double cx = c * cb, cy = c * (ca - cb * cg) / sg;
  k.c = {cx, cy, std::sqrt(std::max(0.0, c * c - cx * cx - cy * cy))};
  return k;
}

}  // namespace

System read_car(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw ReadError("cannot open " + path);
  System s;
  s.source_format = "car";
  std::string line;
  bool pbc = false;
  int64_t mol = 1;
  int header = 0;
  // per molecule: "resname_resnum:name" → atom index, for the .mdf's connections
  std::vector<std::map<std::string, uint32_t>> keys(1);
  auto number = [](const std::string& x) {
    char* e = nullptr;
    std::strtod(x.c_str(), &e);
    return e && *e == 0 && !x.empty();
  };
  while (std::getline(in, line)) {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    const auto t = split(line);
    if (t.empty()) continue;
    if (line[0] == '!') { ++header; continue; }
    if (t[0].rfind("PBC=", 0) == 0) {
      pbc = t[0] == "PBC=ON" || t[0] == "PBC=2D";
      continue;
    }
    if (t[0] == "PBC" && t.size() >= 7) {
      s.cell = cell_from(std::stod(t[1]), std::stod(t[2]), std::stod(t[3]), std::stod(t[4]), std::stod(t[5]), std::stod(t[6]));
      continue;
    }
    if (t[0] == "end") {   // the end of a molecule (a second "end" closes the file)
      if (!keys.back().empty()) {
        ++mol;
        keys.emplace_back();
      }
      continue;
    }
    // an atom: name x y z residue number type element charge
    if (t.size() < 9 || !number(t[1]) || !number(t[2]) || !number(t[3]) || !number(t[8])) {
      if (s.title.empty()) s.title = trim(line);   // the title line after PBC=
      continue;
    }
    Atom a;
    a.id = int64_t(s.atoms.size()) + 1;
    a.mol = mol;
    a.pos = {std::stod(t[1]), std::stod(t[2]), std::stod(t[3])};
    a.resname = t[4];
    a.resid = std::atoll(t[5].c_str());
    a.name = t[6];   // the force-field type (ay2, c3, oy1 …): --names types by it
    a.element = element_from_symbol(t[7]);
    if (!a.element) a.element = element_from_name(t[0]);
    a.charge = std::stod(t[8]);
    keys.back()[t[4] + "_" + t[5] + ":" + t[0]] = uint32_t(s.atoms.size());
    s.atoms.push_back(a);
  }
  (void)header;
  if (s.atoms.empty()) throw ReadError(path + ": no atoms (read as a Materials Studio .car file)");
  if (!pbc) s.cell = Cell{};
  s.has_charges = true;
  s.has_mol = true;
  // the bonds: the .mdf beside the .car
  std::filesystem::path mdf = std::filesystem::path(path).replace_extension(".mdf");
  if (!std::filesystem::exists(mdf)) mdf = std::filesystem::path(path).replace_extension(".MDF");
  if (!std::filesystem::exists(mdf)) {
    s.notes.push_back("no .mdf beside the .car: bonds from distances");
    return s;
  }
  std::ifstream md(mdf);
  std::set<std::pair<uint32_t, uint32_t>> seen;
  size_t m = 0;
  int conn_col = 12, n_missing = 0;
  bool started = false, topology = true;   // atoms only in #topology (#atomset's @list lines also start "RES:NAME")
  while (std::getline(md, line)) {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    const auto t = split(line);
    if (t.empty() || line[0] == '!') continue;
    if (t[0][0] == '#' && t[0] != "#") topology = t[0] == "#topology";
    if (!topology) continue;
    if (t[0] == "@column" && t.size() >= 3 && t[2] == "connections") conn_col = std::atoi(t[1].c_str());
    if (t[0] == "@molecule") {
      if (started) ++m;
      started = true;
      continue;
    }
    if (t[0][0] == '@' || t[0][0] == '#' || m >= keys.size()) continue;
    const std::string key = t[0];
    const auto colon = key.find(':');
    if (colon == std::string::npos) continue;
    const std::string res = key.substr(0, colon);
    const auto self = keys[m].find(key);
    if (self == keys[m].end()) { ++n_missing; continue; }
    for (size_t c = size_t(conn_col); c < t.size(); ++c) {
      std::string ref = t[c];
      int order = 1;
      if (const auto sl = ref.find('/'); sl != std::string::npos) {
        const double o = std::atof(ref.substr(sl + 1).c_str());
        order = o == 1.5 ? 4 : int(std::lround(o));
        ref = ref.substr(0, sl);
      }
      if (const auto pc = ref.find('%'); pc != std::string::npos) ref = ref.substr(0, pc);   // a bond across the cell
      if (const auto hs = ref.find('#'); hs != std::string::npos) ref = ref.substr(0, hs);
      const std::string full = ref.find(':') == std::string::npos ? res + ":" + ref : ref;
      const auto other = keys[m].find(full);
      if (other == keys[m].end()) { ++n_missing; continue; }
      const uint32_t i = std::min(self->second, other->second), j = std::max(self->second, other->second);
      if (i == j || !seen.insert({i, j}).second) continue;
      s.bonds.push_back({i, j, order});
    }
  }
  s.bonds_from_file = true;   // the .mdf lists every bond (a metal or an ionic crystal: none)
  if (n_missing) s.notes.push_back(std::to_string(n_missing) + " .mdf entries name atoms the .car does not have");
  s.notes.push_back(std::to_string(s.bonds.size()) + " bonds from " + mdf.filename().string());
  return s;
}

// The same pair written: the cell in Materials Studio's orientation (a along x, b in the xy plane; positions keep their
// fractional coordinates), one block per molecule (bond-connected), each atom's name (its force-field type once one is
// assigned) as the potential type, its charge; the .mdf lists every atom's bonds, "%abc#1" where the partner is the
// image a b c cells over.
void write_car(const System& s, const std::string& path) {
  const std::string stem = std::filesystem::path(path).replace_extension("").string();
  std::ofstream car(stem + ".car"), mdf(stem + ".mdf");
  if (!car) throw std::runtime_error("cannot write " + stem + ".car");
  if (!mdf) throw std::runtime_error("cannot write " + stem + ".mdf");
  const bool pbc = s.cell.valid();
  Cell ms;
  double la = 0, lb = 0, lc = 0, al = 90, be = 90, ga = 90;
  if (pbc) {
    la = norm(s.cell.a), lb = norm(s.cell.b), lc = norm(s.cell.c);
    al = std::acos(std::clamp(dot(s.cell.b, s.cell.c) / (lb * lc), -1.0, 1.0)) / kDeg;
    be = std::acos(std::clamp(dot(s.cell.a, s.cell.c) / (la * lc), -1.0, 1.0)) / kDeg;
    ga = std::acos(std::clamp(dot(s.cell.a, s.cell.b) / (la * lb), -1.0, 1.0)) / kDeg;
    ms = cell_from(la, lb, lc, al, be, ga);
  }
  const size_t n = s.atoms.size();
  std::vector<Vec3> pos(n), frac(n);
  for (size_t i = 0; i < n; ++i) {
    if (pbc) {   // inside the cell, as Materials Studio writes them; bonds across it carry their images in the .mdf
      frac[i] = s.cell.to_fractional(s.atoms[i].pos);
      for (int k = 0; k < 3; ++k) frac[i][k] -= std::floor(frac[i][k]);
      pos[i] = ms.to_cartesian(frac[i]);
    } else {
      pos[i] = s.atoms[i].pos;
    }
  }
  int nmol = 0;
  const auto mol = s.molecules(&nmol);
  std::vector<std::vector<uint32_t>> members(size_t(std::max(nmol, 0)));
  for (size_t i = 0; i < n; ++i) members[size_t(mol[i])].push_back(uint32_t(i));
  const auto nb = s.neighbours();
  // atom names: element and a count within the residue (MS style), residue: the file's, else XXXX and the molecule
  std::vector<std::string> name(n), res(n);
  std::map<std::string, int> count;
  for (size_t m = 0; m < members.size(); ++m)
    for (uint32_t i : members[m]) {
      const auto& a = s.atoms[i];
      const std::string rn = a.resname.empty() ? "XXXX" : a.resname.substr(0, 4);
      const long rid = a.resid > 0 ? long(a.resid) : long(m + 1);
      res[i] = rn + "_" + std::to_string(rid);
      const std::string el = a.element ? element(a.element).symbol : "X";
      name[i] = el + std::to_string(++count[res[i] + ":" + el]);
    }
  std::time_t now = std::time(nullptr);
  char date[64];
  std::strftime(date, sizeof date, "%a %b %d %H:%M:%S %Y", std::localtime(&now));
  car << "!BIOSYM archive 3\n" << (pbc ? "PBC=ON" : "PBC=OFF") << "\n" << (s.title.empty() ? "CAPS structure" : s.title) << "\n!DATE " << date << "\n";
  char b[256];
  if (pbc) {
    std::snprintf(b, sizeof b, "PBC %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f (P1)\n", la, lb, lc, al, be, ga);
    car << b;
  }
  mdf << "!BIOSYM molecular_data 4\n\n!Date: " << date << "   written by CAPS\n\n#topology\n\n"
      << "@column 1 element\n@column 2 atom_type\n@column 3 charge_group\n@column 4 isotope\n@column 5 formal_charge\n"
         "@column 6 charge\n@column 7 switching_atom\n@column 8 oop_flag\n@column 9 chirality_flag\n@column 10 occupancy\n"
         "@column 11 xray_temp_factor\n@column 12 connections\n\n";
  for (size_t m = 0; m < members.size(); ++m) {
    mdf << "@molecule MOL" << m + 1 << "\n\n";
    for (uint32_t i : members[m]) {
      const auto& a = s.atoms[i];
      const std::string el = a.element ? element(a.element).symbol : "X";
      const std::string type = a.name.empty() ? "?" : a.name.substr(0, 7);
      const std::string rn = res[i].substr(0, res[i].find('_')), rid = res[i].substr(res[i].find('_') + 1);
      std::snprintf(b, sizeof b, "%-5s %14.9f %14.9f %14.9f %-4s %-6s %-7s %-2s %7.4f\n", name[i].c_str(), pos[i][0], pos[i][1], pos[i][2],
                    rn.c_str(), rid.c_str(), type.c_str(), el.c_str(), a.charge);
      car << b;
      std::string conn;
      for (uint32_t j : nb[i]) {
        conn += " " + (res[j] == res[i] ? name[j] : res[j] + ":" + name[j]);
        if (pbc) {   // the partner's image: where the bond reaches it
          const Vec3 df = frac[j] - frac[i];
          const int t[3] = {-int(std::lround(df[0])), -int(std::lround(df[1])), -int(std::lround(df[2]))};
          if (t[0] || t[1] || t[2]) conn += "%" + std::to_string(t[0]) + std::to_string(t[1]) + std::to_string(t[2]) + "#1";
        }
      }
      std::snprintf(b, sizeof b, "%-20s %-2s %-7s ?     0  0  %8.4f 0 0 8 1.0000  0.0000", (res[i] + ":" + name[i]).c_str(), el.c_str(), type.c_str(), a.charge);
      mdf << b << conn << "\n";
    }
    mdf << "\n";
    car << "end\n";
  }
  car << "end\n";
  mdf << "#end\n";
}

}  // namespace caps
