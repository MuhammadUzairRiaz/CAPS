// MDL molfile / SD file (V2000 and V3000 connection tables; every record one molecule) and VASP POSCAR / CONTCAR.
#include <algorithm>
#include <cctype>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <map>

#include "caps/elements.hpp"
#include "caps/spacegroup.hpp"
#include "caps/io.hpp"
#include "caps/typing.hpp"
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
int mdl_order(int t) { return t >= 1 && t <= 4 ? t : t == 9 ? kBondDative : 0; }   // V3000 9: coordination

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
  // molecules follow the bonds, which never cross records: a record that holds a whole cell (as write_sdf writes
  // one) splits into its chains, and a record never merges with the next
  s.has_mol = false;
  if (records > 1) s.notes.push_back(std::to_string(records) + " records from the SD file");
  return s;
}

void write_sdf(const System& s, const std::string& path) {
  std::ofstream out(path);
  if (!out) throw std::runtime_error("cannot write " + path);
  write_sdf_record(s, out, {});
}

void write_sdf_record(const System& s, std::ostream& out, const std::vector<std::pair<std::string, std::string>>& data) {
  const size_t na = s.atoms.size(), nb = s.bonds.size();
  // formal charges (MDL keeps formal charges, not partial ones): from the structure's perceived chemistry
  std::vector<int> formal(na, 0);
  try { formal = perceive(s).charge; } catch (...) {}
  auto order = [](int o) { return o >= 1 && o <= 4 ? o : 1; };   // V2000 has no coordination type: a single bond
  auto order3 = [](int o) { return o >= 1 && o <= 4 ? o : o == kBondDative ? 9 : 1; };   // V3000: 9 coordination
  const bool v3000 = na > 999 || nb > 999;
  // the header block: the name (at most 80 bytes, cut at a whole UTF-8 character), then the program line in the MDL
  // columns: initials (1-2), program (3-10), date and time (11-20, left blank so the file is reproducible) and the
  // dimensional code (21-22)
  std::string head = s.title.empty() ? std::string("CAPS") : s.title;
  for (char& ch : head) if (ch == '\n' || ch == '\r') ch = ' ';
  if (head.size() > 80) {
    size_t cut = 80;
    while (cut > 0 && (static_cast<unsigned char>(head[cut]) & 0xC0) == 0x80) --cut;
    head.resize(cut);
  }
  out << head << "\n  CAPS" << std::string(14, ' ') << "3D\n\n";   // 2 + 8 (CAPS padded) + 10 blank columns
  char b[160];
  if (!v3000) {
    std::snprintf(b, sizeof b, "%3zu%3zu  0  0  0  0  0  0  0  0999 V2000\n", na, nb);
    out << b;
    for (size_t i = 0; i < na; ++i) {
      const auto& a = s.atoms[i];
      std::snprintf(b, sizeof b, "%10.4f%10.4f%10.4f %-3s 0  0  0  0  0  0  0  0  0  0  0  0\n", a.pos[0], a.pos[1], a.pos[2], element(a.element).symbol);
      out << b;
    }
    for (const auto& bd : s.bonds) {
      std::snprintf(b, sizeof b, "%3u%3u%3d  0  0  0  0\n", bd.i + 1, bd.j + 1, order(bd.order));
      out << b;
    }
    std::vector<size_t> charged;
    for (size_t i = 0; i < na; ++i) if (formal[i] != 0) charged.push_back(i);
    for (size_t k = 0; k < charged.size(); k += 8) {
      const size_t m = std::min<size_t>(8, charged.size() - k);
      std::snprintf(b, sizeof b, "M  CHG%3zu", m);
      out << b;
      for (size_t q = k; q < k + m; ++q) { std::snprintf(b, sizeof b, " %3zu %3d", charged[q] + 1, formal[charged[q]]); out << b; }
      out << "\n";
    }
  } else {
    out << "  0  0  0     0  0            999 V3000\n";
    out << "M  V30 BEGIN CTAB\n";
    std::snprintf(b, sizeof b, "M  V30 COUNTS %zu %zu 0 0 0\n", na, nb);
    out << b << "M  V30 BEGIN ATOM\n";
    for (size_t i = 0; i < na; ++i) {
      const auto& a = s.atoms[i];
      std::snprintf(b, sizeof b, "M  V30 %zu %s %.4f %.4f %.4f 0", i + 1, element(a.element).symbol, a.pos[0], a.pos[1], a.pos[2]);
      out << b;
      if (formal[i] != 0) out << " CHG=" << formal[i];
      out << "\n";
    }
    out << "M  V30 END ATOM\nM  V30 BEGIN BOND\n";
    for (size_t k = 0; k < nb; ++k) {
      std::snprintf(b, sizeof b, "M  V30 %zu %d %u %u\n", k + 1, order3(s.bonds[k].order), s.bonds[k].i + 1, s.bonds[k].j + 1);
      out << b;
    }
    out << "M  V30 END BOND\nM  V30 END CTAB\n";
  }
  out << "M  END\n";
  // data items (Dalby et al. 1992): a header line >  <name>, the value, a blank line
  for (const auto& [name, value] : data) out << ">  <" << name << ">\n" << value << "\n\n";
  out << "$$$$\n";
}

// a CIF 1.1 data block code from a title: printable ASCII only (every other byte, whitespace and each byte of a
// multi-byte UTF-8 character, becomes '_'), runs of '_' collapsed, trimmed, at most 75 characters; "caps" when empty
std::string cif_block_name(const std::string& title) {
  std::string name;
  for (char ch : title) {
    const unsigned char u = static_cast<unsigned char>(ch);
    const char c = u >= 0x21 && u <= 0x7E ? ch : '_';
    if (c == '_' && !name.empty() && name.back() == '_') continue;
    name += c;
  }
  while (!name.empty() && name.front() == '_') name.erase(name.begin());
  if (name.size() > 75) name.resize(75);
  while (!name.empty() && name.back() == '_') name.pop_back();
  return name.empty() ? std::string("caps") : name;
}

// CIF: the cell and the atoms at their fractional coordinates, with the element as type symbol. The space group kept
// with the structure is written (its operations and the asymmetric unit) when every operation still maps the atoms
// onto themselves; otherwise the file is P 1 with every atom.
void write_cif(const System& s, const std::string& path) {
  if (!s.cell.valid()) throw std::runtime_error("a CIF needs a periodic cell: this structure has none");
  std::ofstream out(path);
  if (!out) throw std::runtime_error("cannot write " + path);
  const Cell& c = s.cell;
  const double la = norm(c.a), lb = norm(c.b), lc = norm(c.c);
  auto ang = [](const Vec3& u, const Vec3& v) { return std::acos(std::clamp(dot(u, v) / (norm(u) * norm(v)), -1.0, 1.0)) * 180 / M_PI; };
  const std::string name = cif_block_name(s.title);
  char b[200];
  // the formula in Hill order (C, H, then alphabetical; alphabetical without carbon)
  std::map<std::string, int> nel;
  for (const auto& a : s.atoms) ++nel[element(a.element).symbol];
  std::vector<std::string> order;
  const bool carbon = nel.count("C") > 0;
  if (carbon) { order.push_back("C"); if (nel.count("H")) order.push_back("H"); }
  for (const auto& [e, n] : nel) if (!(carbon && (e == "C" || e == "H"))) order.push_back(e);
  std::string formula;
  for (const auto& e : order) formula += (formula.empty() ? "" : " ") + e + (nel[e] > 1 ? std::to_string(nel[e]) : "");
  // both the current (_space_group_*) and the older (_symmetry_*) names, so every reader finds the space group
  // the space group kept with the structure, when every operation still maps the atoms onto themselves: the asymmetric
  // unit and the operations; P 1 otherwise
  const bool keep_sym = s.symmetry.ops.size() > 1 && symmetry_holds(s, s.symmetry, 0.05);
  const std::string hm = keep_sym && !s.symmetry.name.empty() ? s.symmetry.name : keep_sym ? std::string("?") : std::string("P 1");
  const int itn = keep_sym ? s.symmetry.number : 1;
  out << "# written by CAPS\ndata_" << name << "\n";
  // the full title (the block code keeps only its ASCII part), quoted so that the reader takes it back as the name
  // CIF 1.1 allows printable ASCII only: the common non-ASCII characters of a CAPS title in an ASCII form, any other
  // character as '?'
  std::string title;
  {
    static const std::pair<const char*, const char*> kAscii[] = {
        {"\u00B7", "."}, {"\u00D7", "x"}, {"\u00C5", "A"}, {"\u212B", "A"}, {"\u00B0", " deg"}, {"\u2013", "-"}, {"\u2014", "-"}, {"\u2212", "-"},
        {"\u03B1", "alpha"}, {"\u03B2", "beta"}, {"\u03B3", "gamma"}, {"\u03B4", "delta"}, {"\u03B5", "epsilon"}, {"\u03BB", "lambda"},
        {"\u03BC", "mu"}, {"\u03C0", "pi"}, {"\u03C3", "sigma"}, {"\u2192", "->"}, {"\u2018", "'"}, {"\u2019", "'"}, {"\u201C", "\""}, {"\u201D", "\""},
        {"\u00B2", "2"}, {"\u00B3", "3"}, {"\u2070", "0"}, {"\u00B9", "1"}, {"\u2080", "0"}, {"\u2081", "1"}, {"\u2082", "2"}, {"\u2083", "3"}, {"\u2084", "4"}};
    for (size_t i = 0; i < s.title.size();) {
      const unsigned char u = static_cast<unsigned char>(s.title[i]);
      if (u < 0x80) {
        title += (u == '\n' || u == '\r' || u == '\t') ? ' ' : (u >= 0x20 && u < 0x7F ? char(u) : '?');
        ++i;
        continue;
      }
      const size_t len = u >= 0xF0 ? 4 : u >= 0xE0 ? 3 : u >= 0xC0 ? 2 : 1;
      const std::string ch = s.title.substr(i, len);
      const char* rep = "?";
      for (const auto& [from, to] : kAscii)
        if (ch == from) { rep = to; break; }
      title += rep;
      i += len;
    }
  }
  if (!title.empty() && title.find_first_not_of(' ') != std::string::npos) {
    if (title.find('\'') == std::string::npos) out << "_chemical_name_common '" << title << "'\n";
    else if (title.find('"') == std::string::npos) out << "_chemical_name_common \"" << title << "\"\n";
    else out << "_chemical_name_common\n;" << title << "\n;\n";
  }
  out << "_chemical_formula_sum '" << formula << "'\n"
      << "_space_group_name_H-M_alt '" << hm << "'\n";
  if (itn > 0) out << "_space_group_IT_number " << itn << "\n";
  out << "_symmetry_space_group_name_H-M '" << hm << "'\n";
  if (itn > 0) out << "_symmetry_Int_Tables_number " << itn << "\n";
  std::snprintf(b, sizeof b, "_cell_length_a %.6f\n_cell_length_b %.6f\n_cell_length_c %.6f\n_cell_angle_alpha %.6f\n_cell_angle_beta %.6f\n_cell_angle_gamma %.6f\n_cell_volume %.4f\n",
                la, lb, lc, ang(c.b, c.c), ang(c.a, c.c), ang(c.a, c.b), c.volume());
  out << b << "loop_\n_space_group_symop_operation_xyz\n";
  if (keep_sym)
    for (const auto& op : s.symmetry.ops) out << "'" << symop_string(op) << "'\n";
  else
    out << "'x, y, z'\n";
  out << "loop_\n_atom_site_label\n_atom_site_type_symbol\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n";
  if (s.has_charges) out << "_atom_site_charge\n";
  std::map<int, int> count;
  std::vector<size_t> rows;
  if (keep_sym) rows = asymmetric_unit(s, s.symmetry, 0.05);
  else
    for (size_t i = 0; i < s.atoms.size(); ++i) rows.push_back(i);
  for (size_t r : rows) {
    const Atom& a = s.atoms[r];
    Vec3 f = c.to_fractional(a.pos);
    for (int k = 0; k < 3; ++k) f[size_t(k)] -= std::floor(f[size_t(k)]);
    const char* sym = element(a.element).symbol;
    std::snprintf(b, sizeof b, "%s%d %s %.6f %.6f %.6f", sym, ++count[a.element], sym, f[0], f[1], f[2]);
    out << b;
    if (s.has_charges) { std::snprintf(b, sizeof b, " %.5f", a.charge); out << b; }
    out << "\n";
  }
}

void write_poscar(const System& s, const std::string& path, const std::vector<char>& fixed, std::vector<uint32_t>* order) {
  if (!s.cell.valid()) throw std::runtime_error("a POSCAR needs a periodic cell: this structure has none");
  if (s.atoms.empty()) throw std::runtime_error("a POSCAR needs atoms");
  std::ofstream out(path);
  if (!out) throw std::runtime_error("cannot write " + path);
  // the species in the order each first appears, and the atoms grouped by them
  std::vector<int> species;
  for (const auto& a : s.atoms) if (std::find(species.begin(), species.end(), a.element) == species.end()) species.push_back(a.element);
  std::vector<uint32_t> ord;
  for (int e : species)
    for (uint32_t i = 0; i < s.atoms.size(); ++i) if (s.atoms[i].element == e) ord.push_back(i);
  std::string title = s.title.empty() ? std::string("written by CAPS") : s.title;
  for (char& ch : title) if (ch == '\n' || ch == '\r') ch = ' ';
  char b[200];
  out << title << "\n1.0\n";
  // rounded at 1e-12 (Å and fractions): the round-off of building the cell (3e-16, 0.4999999999999999) is not written
  auto clean = [](double x) { const double r = std::round(x * 1e12) / 1e12; return r == 0 ? 0.0 : r; };
  for (const Vec3& v : {s.cell.a, s.cell.b, s.cell.c}) {
    std::snprintf(b, sizeof b, "  %21.16f %21.16f %21.16f\n", clean(v[0]), clean(v[1]), clean(v[2]));
    out << b;
  }
  for (int e : species) out << "  " << element(e).symbol;
  out << "\n";
  for (int e : species) out << "  " << std::count_if(s.atoms.begin(), s.atoms.end(), [&](const Atom& a) { return a.element == e; });
  out << "\n";
  const bool sel = !fixed.empty() && std::any_of(fixed.begin(), fixed.end(), [](char f) { return f != 0; });
  if (sel) out << "Selective dynamics\n";
  out << "Direct\n";
  for (uint32_t i : ord) {
    Vec3 f = s.cell.to_fractional(s.atoms[i].pos);
    for (int k = 0; k < 3; ++k) {
      f[size_t(k)] = clean(f[size_t(k)] - std::floor(f[size_t(k)]));
      if (f[size_t(k)] >= 1.0) f[size_t(k)] = 0.0;   // −1e-17 wraps to 1.0
    }
    std::snprintf(b, sizeof b, "  %19.16f %19.16f %19.16f", f[0], f[1], f[2]);
    out << b;
    if (sel) {   // per coordinate: F where it is held
      const char f = i < fixed.size() ? fixed[i] : 0;
      for (int k = 0; k < 3; ++k) out << (f && holds_axis(f, k) ? "   F" : "   T");
    }
    out << "\n";
  }
  if (order) *order = ord;
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
