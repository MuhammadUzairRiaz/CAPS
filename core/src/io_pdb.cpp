// PDB reader (ATOM/HETATM, CRYST1, CONECT, MODEL/ENDMDL) and writer. Fixed columns per the wwPDB v3.3 format.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <unordered_map>

#include "caps/elements.hpp"
#include "caps/io.hpp"
#include "io_util.hpp"

namespace caps {
namespace {

std::string col(const std::string& l, size_t a, size_t b) {   // 1-based inclusive columns
  if (l.size() < a) return "";
  return trim(l.substr(a - 1, std::min(b, l.size()) - a + 1));
}

Cell cryst1(double a, double b, double c, double al, double be, double ga) {
  const double ca = std::cos(al * M_PI / 180), cb = std::cos(be * M_PI / 180), cg = std::cos(ga * M_PI / 180), sg = std::sin(ga * M_PI / 180);
  Cell k;
  k.a = {a, 0, 0};
  k.b = {b * cg, b * sg, 0};
  const double cx = c * cb, cy = c * (ca - cb * cg) / sg;
  k.c = {cx, cy, std::sqrt(std::max(0.0, c * c - cx * cx - cy * cy))};
  return k;
}

}  // namespace

Trajectory read_pdb(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw ReadError("cannot open " + path);
  Trajectory tr;
  System cur;
  cur.source_format = "pdb";
  Cell cell;
  std::vector<Vec3> pos;
  std::unordered_map<long, uint32_t> serial;
  std::vector<std::pair<long, long>> conect;
  bool first_done = false;
  int64_t chain_mol = 0;
  std::string prev_chain = "\x01";
  auto finish = [&] {
    if (pos.empty()) return;
    if (!first_done) {
      cur.cell = cell;
      for (auto [a, b] : conect) {
        auto i = serial.find(a), j = serial.find(b);
        if (i != serial.end() && j != serial.end() && i->second < j->second) cur.bonds.push_back({i->second, j->second});
      }
      cur.bonds_from_file = !cur.bonds.empty();
      tr.topology = cur;
      first_done = true;
    }
    tr.positions.push_back(pos);
    tr.cells.push_back(cell);
    tr.timesteps.push_back(int64_t(tr.positions.size() - 1));
    pos.clear();
  };
  std::string l;
  size_t lineno = 0;
  while (std::getline(in, l)) {
    ++lineno;
    const std::string rec = l.substr(0, std::min<size_t>(6, l.size()));
    if (rec == "CRYST1") {
      cell = cryst1(std::stod(col(l, 7, 15)), std::stod(col(l, 16, 24)), std::stod(col(l, 25, 33)), std::stod(col(l, 34, 40)), std::stod(col(l, 41, 47)),
                    std::stod(col(l, 48, 54)));
      // PDB files for non-periodic molecules carry a unit cell of 1 Å; ignore it.
      if (norm(cell.a) <= 1.0 + 1e-9 && norm(cell.b) <= 1.0 + 1e-9) cell = Cell{};
    } else if (rec == "ATOM  " || rec == "HETATM") {
      if (l.size() < 54) throw ReadError(path + ":" + std::to_string(lineno) + ": ATOM record shorter than 54 columns");
      const Vec3 r{std::stod(col(l, 31, 38)), std::stod(col(l, 39, 46)), std::stod(col(l, 47, 54))};
      if (!first_done) {
        Atom a;
        a.id = static_cast<int64_t>(cur.atoms.size() + 1);
        a.name = col(l, 13, 16);
        a.resname = col(l, 18, 20);
        const std::string chain = col(l, 22, 22);
        if (chain != prev_chain) { ++chain_mol; prev_chain = chain; }
        a.mol = chain_mol;
        const std::string el = col(l, 77, 78);
        a.element = el.empty() ? element_from_name(a.name) : element_from_symbol(el);
        a.pos = r;
        const std::string ser = col(l, 7, 11);
        if (!ser.empty() && ser.find_first_not_of("0123456789") == std::string::npos) serial[std::stol(ser)] = uint32_t(cur.atoms.size());
        cur.atoms.push_back(a);
      }
      pos.push_back(r);
    } else if (rec == "CONECT" && !first_done) {
      const std::string a = col(l, 7, 11);
      if (a.empty()) continue;
      for (size_t c = 12; c + 4 <= l.size() && c <= 27; c += 5) {
        const std::string b = col(l, c, c + 4);
        if (!b.empty()) conect.emplace_back(std::stol(a), std::stol(b));
      }
    } else if (rec == "ENDMDL") {
      finish();
    }
  }
  finish();
  if (tr.positions.empty()) throw ReadError(path + ": no ATOM or HETATM records");
  if (!tr.topology.bonds_from_file) tr.topology.notes.push_back("no CONECT records; bonds will be perceived from distances");
  tr.topology.notes.push_back("molecules found from bonds; chain identifiers are kept per atom");
  tr.topology.has_mol = false;   // molecules follow bonds; chain ids stay per atom
  return tr;
}

void write_pdb(const System& s, const std::string& path) {
  std::FILE* f = std::fopen(path.c_str(), "w");
  if (!f) throw ReadError("cannot write " + path);
  if (s.cell.valid()) {
    const double a = norm(s.cell.a), b = norm(s.cell.b), c = norm(s.cell.c);
    const double al = std::acos(dot(s.cell.b, s.cell.c) / (b * c)) * 180 / M_PI, be = std::acos(dot(s.cell.a, s.cell.c) / (a * c)) * 180 / M_PI,
                 ga = std::acos(dot(s.cell.a, s.cell.b) / (a * b)) * 180 / M_PI;
    std::fprintf(f, "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1\n", a, b, c, al, be, ga);
  }
  const bool big = s.atoms.size() > 99999;
  for (size_t i = 0; i < s.atoms.size(); ++i) {
    const auto& at = s.atoms[i];
    const std::string nm = at.name.empty() ? element(at.element).symbol : at.name.substr(0, 4);
    const std::string res = at.resname.empty() ? "MOL" : at.resname.substr(0, 3);
    std::fprintf(f, "HETATM%5zu %-4s %3s A%4lld    %8.3f%8.3f%8.3f  1.00  0.00          %2s\n", big ? (i + 1) % 100000 : i + 1, nm.c_str(), res.c_str(),
                 static_cast<long long>(at.mol % 10000), at.pos[0], at.pos[1], at.pos[2], element(at.element).symbol);
  }
  if (!big && !s.bonds.empty()) {
    const auto nb = s.neighbours();
    for (size_t i = 0; i < nb.size(); ++i)
      for (size_t k = 0; k < nb[i].size(); k += 4) {
        std::fprintf(f, "CONECT%5zu", i + 1);
        for (size_t q = k; q < std::min(k + 4, nb[i].size()); ++q) std::fprintf(f, "%5u", nb[i][q] + 1);
        std::fprintf(f, "\n");
      }
  }
  std::fprintf(f, "END\n");
  std::fclose(f);
}

}  // namespace caps
