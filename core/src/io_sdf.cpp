// MDL molfile / SD file (V2000 and V3000 connection tables; every record one molecule) and VASP POSCAR / CONTCAR.
#include <cmath>
#include <filesystem>
#include <fstream>
#include <map>

#include "caps/elements.hpp"
#include "caps/io.hpp"
#include "io_util.hpp"

namespace caps {

namespace {

double num(const std::string& s, const std::string& what) {
  try {
    size_t k = 0;
    const double v = std::stod(s, &k);
    return v;
  } catch (...) {
    throw ReadError("could not read " + what + " from '" + s + "'");
  }
}

std::string field(const std::string& l, size_t from, size_t len) { return l.size() > from ? trim(l.substr(from, len)) : std::string(); }

// MDL bond type → CAPS order (1, 2, 3 single to triple; 4 aromatic)
int mdl_order(int t) { return t >= 1 && t <= 4 ? t : 0; }

}  // namespace

System read_sdf(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw ReadError("cannot open " + path);
  System s;
  s.source_format = "sdf";
  std::vector<std::string> lines;
  for (std::string l; std::getline(in, l);) {
    if (!l.empty() && l.back() == '\r') l.pop_back();
    lines.push_back(l);
  }
  size_t k = 0, records = 0;
  int64_t mol = 0;
  auto charge_of = [](int ccc) { return ccc >= 1 && ccc <= 7 && ccc != 4 ? double(4 - ccc) : 0.0; };   // 1 → +3 … 7 → −3
  while (k + 3 < lines.size()) {
    // header: name, program line, comment
    const std::string name = trim(lines[k]);
    if (records == 0) s.title = name;
    k += 3;
    const std::string counts = lines[k++];
    const size_t base = s.atoms.size();
    ++mol;
    if (counts.find("V3000") != std::string::npos) {
      bool in_atoms = false, in_bonds = false;
      std::map<int64_t, uint32_t> index;
      for (; k < lines.size(); ++k) {
        const std::string& l = lines[k];
        if (l.rfind("M  END", 0) == 0) { ++k; break; }
        if (l.rfind("M  V30 ", 0) != 0) continue;
        const std::string body = trim(l.substr(7));
        if (body.rfind("BEGIN ATOM", 0) == 0) { in_atoms = true; continue; }
        if (body.rfind("END ATOM", 0) == 0) { in_atoms = false; continue; }
        if (body.rfind("BEGIN BOND", 0) == 0) { in_bonds = true; continue; }
        if (body.rfind("END BOND", 0) == 0) { in_bonds = false; continue; }
        const auto t = split(body);
        if (in_atoms && t.size() >= 5) {
          Atom a;
          a.id = int64_t(s.atoms.size() + 1);
          a.mol = mol;
          a.name = t[1];
          a.element = element_from_symbol(t[1]);
          a.pos = {num(t[2], "x"), num(t[3], "y"), num(t[4], "z")};
          for (size_t q = 6; q < t.size(); ++q)
            if (t[q].rfind("CHG=", 0) == 0) a.charge = num(t[q].substr(4), "charge"), s.has_charges = true;
          index[int64_t(num(t[0], "atom index"))] = uint32_t(s.atoms.size());
          s.atoms.push_back(a);
        } else if (in_bonds && t.size() >= 4) {
          const auto i = index.find(int64_t(num(t[2], "bond atom"))), j = index.find(int64_t(num(t[3], "bond atom")));
          if (i == index.end() || j == index.end()) throw ReadError(path + ": a bond names an atom that is not there");
          s.bonds.push_back({i->second, j->second, mdl_order(int(num(t[1], "bond type")))});
        }
      }
    } else {
      const int na = int(num(field(counts, 0, 3), "the atom count")), nb = int(num(field(counts, 3, 3), "the bond count"));
      if (na < 0 || nb < 0 || k + size_t(na + nb) > lines.size()) throw ReadError(path + ": the counts line does not match the file");
      for (int i = 0; i < na; ++i) {
        const std::string& l = lines[k++];
        Atom a;
        a.id = int64_t(s.atoms.size() + 1);
        a.mol = mol;
        a.pos = {num(field(l, 0, 10), "x"), num(field(l, 10, 10), "y"), num(field(l, 20, 10), "z")};
        a.name = field(l, 31, 3);
        a.element = element_from_symbol(a.name);
        const std::string ccc = field(l, 36, 3);
        if (!ccc.empty() && ccc != "0") a.charge = charge_of(int(num(ccc, "charge code"))), s.has_charges = true;
        s.atoms.push_back(a);
      }
      for (int b = 0; b < nb; ++b) {
        const std::string& l = lines[k++];
        const int i = int(num(field(l, 0, 3), "bond atom")), j = int(num(field(l, 3, 3), "bond atom"));
        if (i < 1 || j < 1 || i > na || j > na) throw ReadError(path + ": a bond names an atom that is not there");
        s.bonds.push_back({uint32_t(base + size_t(i - 1)), uint32_t(base + size_t(j - 1)), mdl_order(int(num(field(l, 6, 3), "bond type")))});
      }
      // properties: M  CHG overrides the atom block's charges
      bool chg = false;
      for (; k < lines.size(); ++k) {
        const std::string& l = lines[k];
        if (l.rfind("M  END", 0) == 0) { ++k; break; }
        if (l.rfind("M  CHG", 0) == 0) {
          if (!chg)
            for (size_t q = base; q < s.atoms.size(); ++q) s.atoms[q].charge = 0;
          chg = true;
          const auto t = split(l.substr(6));
          for (size_t q = 1; q + 1 < t.size(); q += 2) {
            const int i = int(num(t[q], "charged atom"));
            if (i >= 1 && i <= na) s.atoms[base + size_t(i - 1)].charge = num(t[q + 1], "charge");
          }
          s.has_charges = true;
        }
      }
    }
    ++records;
    // data items up to the record end
    while (k < lines.size() && lines[k].rfind("$$$$", 0) != 0) ++k;
    if (k < lines.size()) ++k;
  }
  if (s.atoms.empty()) throw ReadError(path + ": no atoms found (read as an MDL molfile / SD file)");
  for (const auto& a : s.atoms)
    if (a.element == 0) throw ReadError(path + ": unknown element '" + a.name + "'");
  s.bonds_from_file = true;
  s.has_mol = true;
  if (records > 1) s.notes.push_back(std::to_string(records) + " molecules from the SD file, one molecule each");
  return s;
}

System read_poscar(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw ReadError("cannot open " + path);
  std::vector<std::string> lines;
  for (std::string l; std::getline(in, l);) {
    if (!l.empty() && l.back() == '\r') l.pop_back();
    lines.push_back(l);
  }
  if (lines.size() < 8) throw ReadError(path + ": too short for a POSCAR");
  System s;
  s.source_format = "poscar";
  s.title = trim(lines[0]);
  const double scale = num(split(lines[1]).at(0), "the scale factor");
  Vec3 v[3];
  for (int r = 0; r < 3; ++r) {
    const auto t = split(lines[size_t(2 + r)]);
    if (t.size() < 3) throw ReadError(path + ": lattice vector " + std::to_string(r + 1) + " is incomplete");
    v[r] = {num(t[0], "a lattice vector"), num(t[1], "a lattice vector"), num(t[2], "a lattice vector")};
  }
  s.cell.a = v[0];
  s.cell.b = v[1];
  s.cell.c = v[2];
  // a negative scale is the cell volume
  const double f = scale < 0 ? std::cbrt(-scale / s.cell.volume()) : scale;
  s.cell.a = s.cell.a * f, s.cell.b = s.cell.b * f, s.cell.c = s.cell.c * f;
  size_t k = 5;
  std::vector<std::string> species;
  auto t = split(lines[k]);
  if (!t.empty() && !std::isdigit(static_cast<unsigned char>(t[0][0]))) {   // VASP 5: the species line
    species = t;
    t = split(lines[++k]);
  } else {
    // VASP 4: the species may be named in the title line
    for (const auto& w : split(s.title))
      if (element_from_symbol(w) > 0) species.push_back(w);
  }
  std::vector<int> counts;
  for (const auto& w : t) counts.push_back(int(num(w, "an atom count")));
  if (species.size() != counts.size())
    throw ReadError(path + ": name the species (VASP 5's line of symbols after the lattice) — " + std::to_string(counts.size()) + " counts, " +
                    std::to_string(species.size()) + " names");
  std::string mode = lower(trim(lines.at(++k)));
  if (!mode.empty() && mode[0] == 's') mode = lower(trim(lines.at(++k)));   // Selective dynamics
  const bool cart = !mode.empty() && (mode[0] == 'c' || mode[0] == 'k');
  ++k;
  for (size_t sp = 0; sp < species.size(); ++sp) {
    const int e = element_from_symbol(species[sp].substr(0, species[sp].find_first_of("_/")));
    if (e == 0) throw ReadError(path + ": unknown element '" + species[sp] + "'");
    for (int i = 0; i < counts[sp]; ++i) {
      if (k >= lines.size()) throw ReadError(path + ": fewer positions than atoms");
      const auto p = split(lines[k++]);
      if (p.size() < 3) throw ReadError(path + ": a position is incomplete");
      const Vec3 r{num(p[0], "a position"), num(p[1], "a position"), num(p[2], "a position")};
      Atom a;
      a.id = int64_t(s.atoms.size() + 1);
      a.element = e;
      a.name = species[sp];
      a.pos = cart ? r * f : s.cell.to_cartesian(r);
      s.atoms.push_back(a);
    }
  }
  s.notes.push_back("VASP POSCAR: " + std::to_string(s.atoms.size()) + " atoms, " + (cart ? "Cartesian" : "direct (fractional)") + " positions");
  return s;
}

}  // namespace caps
