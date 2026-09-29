// GROMACS .gro (fixed columns, nm) and XYZ / extended XYZ readers.
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>

#include "caps/elements.hpp"
#include "caps/io.hpp"
#include "io_util.hpp"

namespace caps {

Trajectory read_gro(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw ReadError("cannot open " + path);
  Trajectory tr;
  std::string line;
  size_t lineno = 0;
  bool first = true;
  while (std::getline(in, line)) {
    ++lineno;
    const std::string title = trim(line);
    if (title.empty() && in.eof()) break;
    if (!std::getline(in, line)) break;
    ++lineno;
    const size_t n = std::stoull(trim(line));
    System s;
    s.title = title;
    s.source_format = "gro";
    s.atoms.resize(n);
    std::vector<Vec3> pos(n);
    int64_t prev_res = -1, mol = 0;
    std::string prev_resname;
    for (size_t i = 0; i < n; ++i) {
      if (!std::getline(in, line)) throw ReadError(path + ": file ends after " + std::to_string(i) + " of " + std::to_string(n) + " atoms");
      ++lineno;
      if (line.size() < 44) throw ReadError(path + ":" + std::to_string(lineno) + ": line shorter than the 44 fixed columns of .gro");
      Atom& a = s.atoms[i];
      const int64_t resnr = std::stoll(trim(line.substr(0, 5)));
      a.resname = trim(line.substr(5, 5));
      a.resid = resnr;
      a.name = trim(line.substr(10, 5));
      a.id = static_cast<int64_t>(i + 1);   // column 15-20 wraps at 99999; use order
      // Coordinates: fixed width 8 in the standard, but precision may vary: parse the rest by spacing.
      auto t = split(line.substr(20));
      if (t.size() < 3) throw ReadError(path + ":" + std::to_string(lineno) + ": missing coordinates");
      a.pos = {std::stod(t[0]) * 10.0, std::stod(t[1]) * 10.0, std::stod(t[2]) * 10.0};
      a.element = element_from_name(a.name);
      // Residue numbers reset or wrap; a residue change starts a new residue, not necessarily a molecule.
      if (resnr != prev_res || a.resname != prev_resname) { ++mol; prev_res = resnr; prev_resname = a.resname; }
      a.mol = mol;
      pos[i] = a.pos;
    }
    s.has_mol = false;   // residues are not molecules; molecules come from bonds
    if (!std::getline(in, line)) throw ReadError(path + ": missing box line");
    ++lineno;
    auto b = split(line);
    if (b.size() >= 3) {
      const double v[9] = {std::stod(b[0]), std::stod(b[1]), std::stod(b[2]), b.size() >= 9 ? std::stod(b[3]) : 0, b.size() >= 9 ? std::stod(b[4]) : 0,
                           b.size() >= 9 ? std::stod(b[5]) : 0, b.size() >= 9 ? std::stod(b[6]) : 0, b.size() >= 9 ? std::stod(b[7]) : 0,
                           b.size() >= 9 ? std::stod(b[8]) : 0};
      // v1(x) v2(y) v3(z) v1(y) v1(z) v2(x) v2(z) v3(x) v3(y), nm
      s.cell.a = Vec3{v[0], v[3], v[4]} * 10.0;
      s.cell.b = Vec3{v[5], v[1], v[6]} * 10.0;
      s.cell.c = Vec3{v[7], v[8], v[2]} * 10.0;
    }
    s.notes.push_back("molecules found from bonds; residue numbers kept per atom (a .top gives exact molecule boundaries)");
    int unknown = 0;
    for (const auto& a : s.atoms) unknown += a.element == 0;
    if (unknown) s.notes.push_back(std::to_string(unknown) + " atom names did not map to an element");
    if (first) { tr.topology = s; first = false; }
    tr.positions.push_back(std::move(pos));
    tr.cells.push_back(s.cell);
    tr.timesteps.push_back(static_cast<int64_t>(tr.positions.size() - 1));
  }
  if (tr.positions.empty()) throw ReadError(path + ": no frames found");
  return tr;
}

Trajectory read_xyz(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw ReadError("cannot open " + path);
  Trajectory tr;
  std::string line;
  bool first = true;
  while (std::getline(in, line)) {
    if (trim(line).empty()) continue;
    const size_t n = std::stoull(trim(line));
    std::string comment;
    std::getline(in, comment);
    System s;
    s.source_format = "xyz";
    s.title = trim(comment);
    // extended XYZ keys are not a title: Lattice="…", Properties=…, pbc="…" come out of it
    for (const char* key : {"Lattice=", "Properties=", "pbc="}) {
      auto k = s.title.find(key);
      if (k == std::string::npos) continue;
      size_t e = k + std::strlen(key);
      if (e < s.title.size() && s.title[e] == '"') e = s.title.find('"', e + 1);
      else e = s.title.find(' ', e);
      s.title.erase(k, e == std::string::npos ? std::string::npos : e + 1 - k);
      s.title = trim(s.title);
    }
    // extended XYZ: Lattice="ax ay az bx by bz cx cy cz"
    if (auto p = comment.find("Lattice=\""); p != std::string::npos) {
      auto q = comment.find('"', p + 9);
      auto v = split(comment.substr(p + 9, q - p - 9));
      if (v.size() == 9) {
        s.cell.a = {std::stod(v[0]), std::stod(v[1]), std::stod(v[2])};
        s.cell.b = {std::stod(v[3]), std::stod(v[4]), std::stod(v[5])};
        s.cell.c = {std::stod(v[6]), std::stod(v[7]), std::stod(v[8])};
      }
    }
    s.atoms.resize(n);
    std::vector<Vec3> pos(n);
    for (size_t i = 0; i < n; ++i) {
      if (!std::getline(in, line)) throw ReadError(path + ": file ends inside a frame");
      auto t = split(line);
      if (t.size() < 4) throw ReadError(path + ": atom line needs a symbol and x y z");
      Atom& a = s.atoms[i];
      a.id = static_cast<int64_t>(i + 1);
      a.name = t[0];
      a.element = element_from_symbol(t[0]);
      if (!a.element) a.element = element_from_name(t[0]);
      a.pos = {std::stod(t[1]), std::stod(t[2]), std::stod(t[3])};
      pos[i] = a.pos;
    }
    if (first) { tr.topology = s; first = false; }
    tr.positions.push_back(std::move(pos));
    tr.cells.push_back(s.cell);
    tr.timesteps.push_back(static_cast<int64_t>(tr.positions.size() - 1));
  }
  if (tr.positions.empty()) throw ReadError(path + ": no frames found");
  return tr;
}

void write_xyz_frame(std::ostream& out, const System& s, const std::string& extra) {
  out << s.atoms.size() << "\n";
  const Cell& c = s.cell;
  if (c.valid())
    out << "Lattice=\"" << c.a[0] << " " << c.a[1] << " " << c.a[2] << " " << c.b[0] << " " << c.b[1] << " " << c.b[2] << " " << c.c[0] << " "
        << c.c[1] << " " << c.c[2] << "\" Properties=species:S:1:pos:R:3" << (extra.empty() ? "" : " " + extra) << "\n";
  else
    out << (extra.empty() ? s.title : extra) << "\n";
  for (const auto& a : s.atoms) out << element(a.element).symbol << " " << a.pos[0] << " " << a.pos[1] << " " << a.pos[2] << "\n";
}

void write_xyz(const System& s, const std::string& path) {
  std::ofstream out(path);
  if (!out) throw ReadError("cannot write " + path);
  write_xyz_frame(out, s, "");
}

void write_gro(const System& s, const std::string& path) {
  std::ofstream out(path);
  if (!out) throw ReadError("cannot write " + path);
  write_gro_frame(out, s, s.title.empty() ? "CAPS structure" : s.title);
}

void write_gro_frame(std::ostream& out, const System& s, const std::string& title) {
  out << title << "\n" << s.atoms.size() << "\n";
  const auto mol = s.molecules();
  char b[96];
  for (size_t i = 0; i < s.atoms.size(); ++i) {
    const Atom& a = s.atoms[i];
    const long res = (a.resid > 0 ? long(a.resid) : s.has_mol && a.mol > 0 ? long(a.mol) : long(mol[i] + 1)) % 100000;
    std::string rn = a.resname.empty() ? "MOL" : a.resname.substr(0, 5);
    std::string an = a.name.empty() ? element(a.element).symbol : a.name.substr(0, 5);
    std::snprintf(b, sizeof b, "%5ld%-5s%5s%5ld%8.3f%8.3f%8.3f\n", res, rn.c_str(), an.c_str(), long((i + 1) % 100000), a.pos[0] / 10, a.pos[1] / 10, a.pos[2] / 10);
    out << b;
  }
  const Cell& c = s.cell;
  if (c.valid()) {
    std::snprintf(b, sizeof b, "%10.5f%10.5f%10.5f", c.a[0] / 10, c.b[1] / 10, c.c[2] / 10);
    out << b;
    if (std::fabs(c.a[1]) + std::fabs(c.a[2]) + std::fabs(c.b[0]) + std::fabs(c.b[2]) + std::fabs(c.c[0]) + std::fabs(c.c[1]) > 1e-9) {
      std::snprintf(b, sizeof b, "%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f", c.a[1] / 10, c.a[2] / 10, c.b[0] / 10, c.b[2] / 10, c.c[0] / 10, c.c[1] / 10);
      out << b;
    }
    out << "\n";
  } else {
    out << "   0.00000   0.00000   0.00000\n";
  }
}

}  // namespace caps
