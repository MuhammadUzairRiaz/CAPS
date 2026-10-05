// CAPS literature many-body potentials for a group of atoms (see manybody.hpp).
#include "caps/manybody.hpp"

#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <map>
#include <set>
#include <sstream>

#include "caps/elements.hpp"
#include "caps/uff.hpp"

namespace caps {

namespace {

// words per entry of the element-triplet formats (LAMMPS's NPARAMS_PER_LINE: three elements, then the parameters)
const std::map<std::string, int> kTriplet = {{"tersoff", 17}, {"tersoff/mod", 20}, {"tersoff/mod/c", 21}, {"tersoff/zbl", 21}, {"sw", 14},
                                             {"vashishta", 17}, {"gw", 17}, {"gw/zbl", 21}};
const std::set<std::string> kSetfl = {"eam/alloy", "eam/fs"};
// the Brenner family: C and H, read by LAMMPS in metal units only (the export is then in metal units)
const std::set<std::string> kBrenner = {"airebo", "airebo/morse", "rebo"};
// not written: MEAM names its library entries freely (not by element), BOP / COMB / others need what CAPS does not map
const std::map<std::string, std::string> kRefused = {
    {"meam/c", "meam/c is LAMMPS's old name: use meam"},
    {"bop", "BOP tables are not checked or mapped by CAPS"},
    {"comb", "COMB equilibrates charges itself (fix qeq/comb), which CAPS does not set up"},
    {"comb3", "COMB3 equilibrates charges itself (fix qeq/comb), which CAPS does not set up"},
    {"lcbop", "LCBOP is not mapped by CAPS (use airebo or rebo for carbon)"},
    {"polymorphic", "polymorphic tables are not checked or mapped by CAPS"},
    {"edip", "EDIP is not mapped by CAPS"},
    {"extep", "ExTeP is not mapped by CAPS"}};

bool is_number(const std::string& w) {
  if (w.empty()) return false;
  char* end = nullptr;
  std::strtod(w.c_str(), &end);
  return end && *end == '\0';
}

std::string first_line(const std::string& path) {
  std::ifstream f(path);
  std::string l;
  std::getline(f, l);
  return l;
}

// the word after "key" on the line ("UNITS: metal" → metal)
std::string tag(const std::string& line, const std::string& key) {
  std::istringstream is(line);
  for (std::string w; is >> w;)
    if (w == key) {
      std::string v;
      is >> v;
      return v;
    }
  return "";
}

std::string after(const std::string& line, const std::string& key) {
  const size_t k = line.find(key);
  if (k == std::string::npos) return "";
  std::string v = line.substr(k + key.size());
  for (const char* next : {" DATE:", " UNITS:", " CONTRIBUTOR:", " COMMENT:"}) {
    const size_t e = v.find(next);
    if (e != std::string::npos) v = v.substr(0, e);
  }
  while (!v.empty() && (v.front() == ' ')) v.erase(v.begin());
  while (!v.empty() && (v.back() == ' ' || v.back() == '\r')) v.pop_back();
  return v;
}

std::string base_name(const std::string& p) { return std::filesystem::path(p).filename().string(); }

}  // namespace

std::vector<MeamEntry> meam_library(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw FieldError("cannot read the MEAM library " + path);
  // as LAMMPS reads it: comments out, quotes as separators, 19 words per entry
  std::vector<std::string> w;
  for (std::string l; std::getline(in, l);) {
    if (const size_t h = l.find('#'); h != std::string::npos) l.resize(h);
    for (char& c : l) if (c == '\'') c = ' ';
    std::istringstream is(l);
    for (std::string x; is >> x;) w.push_back(x);
  }
  if (w.empty() || w.size() % 19 != 0)
    throw FieldError(base_name(path) + " is not a MEAM library: its entries are 19 values each (name, lattice, z, atomic number, mass, α, b0–b3, a, Ec, A, t0–t3, ρ0, ibar)");
  std::vector<MeamEntry> r;
  std::set<std::string> seen;
  for (size_t k = 0; k < w.size(); k += 19) {
    if (!seen.insert(w[k]).second) continue;   // LAMMPS keeps the first of a repeated name
    MeamEntry e;
    e.name = w[k], e.lattice = w[k + 1];
    if (!is_number(w[k + 3]) || !is_number(w[k + 4])) throw FieldError(base_name(path) + ": entry " + w[k] + " has no atomic number and mass where they belong");
    e.z = std::stoi(w[k + 3]), e.mass = std::stod(w[k + 4]);
    r.push_back(e);
  }
  return r;
}

bool manybody_caps_converts(const std::string& style) { return kBrenner.count(style) > 0; }

namespace {

// The Brenner family's scalars in LAMMPS's read order (pair_airebo.cpp read_file): those with energy dimension — A_CC …
// A_HH (20–22), BIJc (23–31), the LJ ε (56–58), the torsion ε (62–64), AIREBO-M's Morse ε (65–67). The splines that
// follow (g, P, π^rc, T) are dimensionless and stay as they are.
bool brenner_energy(int k) { return (k >= 20 && k <= 31) || (k >= 56 && k <= 58) || (k >= 62 && k <= 67); }

std::string brenner_to_real(const std::string& text, const std::string& style) {
  const int scalars = style == "airebo/morse" ? 74 : 65;
  std::istringstream is(text);
  std::string out, l;
  int line = 0, k = 0;
  char b[64];
  while (std::getline(is, l)) {
    ++line;
    if (line == 1) {   // the units LAMMPS looks for on the first line
      const size_t u = l.find("UNITS: metal");
      if (u != std::string::npos) l.replace(u, 12, "UNITS: real");
      else l += " UNITS: real";
      out += l + " (converted by CAPS: A, B and the ε's × 23.060549; the splines are dimensionless)\n";
      continue;
    }
    std::string head = l.substr(0, l.find('#'));
    std::istringstream ls(head);
    std::string first;
    if (k < scalars && (ls >> first) && is_number(first)) {
      if (brenner_energy(k)) {
        std::snprintf(b, sizeof b, "%.17g", std::stod(first) * 23.060549);
        const size_t at = l.find(first);
        l.replace(at, first.size(), b);
      }
      ++k;
    }
    out += l + "\n";
  }
  if (k < scalars) throw FieldError("the " + style + " file has " + std::to_string(k) + " parameters before its splines, not " + std::to_string(scalars));
  return out;
}

}  // namespace

const std::vector<std::string>& manybody_styles() {
  static const std::vector<std::string> s = {"tersoff", "tersoff/mod", "tersoff/mod/c", "tersoff/zbl", "sw", "vashishta", "gw", "gw/zbl", "eam/alloy", "eam/fs",
                                             "airebo", "airebo/morse", "rebo", "meam"};
  return s;
}

ForceField manybody_part(const System& g, const ManyBodySpec& spec, std::vector<std::string>* notes) {
  const std::string& st = spec.style;
  if (auto it = kRefused.find(st); it != kRefused.end()) throw FieldError("pair style " + st + ": " + it->second);
  const bool meam = st == "meam";
  if (!kTriplet.count(st) && !kSetfl.count(st) && !kBrenner.count(st) && !meam) {
    std::string l;
    for (const auto& x : manybody_styles()) l += (l.empty() ? "" : ", ") + x;
    throw FieldError("pair style \"" + st + "\" is not one CAPS writes (" + l + ")");
  }
  if (spec.file.empty()) throw FieldError(st + ": no potential file. Give the file from LAMMPS's potentials folder, the NIST Interatomic Potentials Repository or the paper");
  std::ifstream in(spec.file);
  if (!in) throw FieldError(st + ": cannot read the potential file " + spec.file);
  if (g.atoms.empty()) throw FieldError(st + ": the group has no atoms");

  // the group's elements, in order of appearance: one type each
  std::vector<std::string> els;
  std::vector<int> zs;
  std::vector<int> tix(g.atoms.size());
  for (size_t i = 0; i < g.atoms.size(); ++i) {
    const int z = g.atoms[i].element;
    if (z <= 0) throw FieldError(st + ": atom " + std::to_string(i + 1) + " of the group has no element");
    const std::string sym = element(z).symbol;
    auto it = std::find(els.begin(), els.end(), sym);
    if (it == els.end()) els.push_back(sym), zs.push_back(z), it = els.end() - 1;
    tix[i] = int(it - els.begin());
  }

  // the file: its units, and an entry for what the group needs
  const std::string head = first_line(spec.file);
  const std::string file_units = tag(head, "UNITS:");
  std::string units = spec.units;
  if (!file_units.empty() && !units.empty() && file_units != units)
    throw FieldError(base_name(spec.file) + " says its units are " + file_units + ", not " + units + " as given");
  if (units.empty()) units = file_units;
  if (units.empty())
    throw FieldError(base_name(spec.file) + " does not say its units (no \"UNITS:\" on its first line), and LAMMPS would read it as it stands: "
                     "say whether its energies are in eV (metal, as most published files) or kcal/mol (real)");
  if (units != "metal" && units != "real") throw FieldError(base_name(spec.file) + ": units " + units + " (CAPS converts metal or real)");
  if (kBrenner.count(st) && units != "metal") throw FieldError(st + ": LAMMPS reads its file in metal units only, and " + base_name(spec.file) + " is in " + units);
  if (meam && units != "metal") throw FieldError("meam: LAMMPS reads MEAM files in metal units only, and " + base_name(spec.file) + " is in " + units);
  std::set<std::string> have;   // "A B C" triplets, or elements
  std::vector<std::string> meam_entry;   // per element of the group
  std::vector<std::string> meam_extract;   // the entries read, in order
  if (meam) {
    const auto lib = meam_library(spec.file);
    for (size_t t = 0; t < els.size(); ++t) {
      const MeamEntry* e = nullptr;
      auto it = std::find_if(spec.entries.begin(), spec.entries.end(), [&](const auto& kv) { return kv.first == els[t]; });
      if (it != spec.entries.end()) {
        for (const auto& x : lib) if (x.name == it->second) e = &x;
        if (!e) throw FieldError(base_name(spec.file) + " has no entry '" + it->second + "' (for " + els[t] + ")");
      } else {
        for (const auto& x : lib) if (x.z == zs[t]) { e = &x; break; }
        if (!e) throw FieldError(base_name(spec.file) + " has no entry for " + els[t] + " (atomic number " + std::to_string(zs[t]) + ")");
      }
      if (e->z != zs[t]) throw FieldError(base_name(spec.file) + ": entry '" + e->name + "' is atomic number " + std::to_string(e->z) + ", not " + els[t]);
      meam_entry.push_back(e->name);
    }
    // the entries read: those given, in their order (the parameter file's indices), then the group's others
    for (const auto& [el, name] : spec.entries) {
      bool known = false;
      for (const auto& x : lib) known = known || x.name == name;
      if (!known) throw FieldError(base_name(spec.file) + " has no entry '" + name + "' (for " + el + ")");
      if (std::find(meam_extract.begin(), meam_extract.end(), name) == meam_extract.end()) meam_extract.push_back(name);
    }
    for (const auto& n : meam_entry)
      if (std::find(meam_extract.begin(), meam_extract.end(), n) == meam_extract.end()) meam_extract.push_back(n);
    if (!spec.file2.empty()) {
      std::ifstream p2(spec.file2);
      if (!p2) throw FieldError("meam: cannot read the parameter file " + spec.file2);
      const std::string u2 = tag(first_line(spec.file2), "UNITS:");
      if (!u2.empty() && u2 != "metal") throw FieldError("meam: " + base_name(spec.file2) + " is in " + u2 + " units; MEAM is read in metal units only");
      // the highest element index it names: lattce(1,2), Ec(1,2), Cmin(1,1,2) …
      int top = 0;
      for (std::string l; std::getline(p2, l);) {
        if (const size_t h = l.find('#'); h != std::string::npos) l.resize(h);
        const size_t a = l.find('('), b = l.find(')');
        if (a == std::string::npos || b == std::string::npos || b < a) continue;
        std::string in = l.substr(a + 1, b - a - 1);
        for (char& c : in) if (c == ',') c = ' ';
        std::istringstream is(in);
        for (int k; is >> k;) top = std::max(top, k);
      }
      if (top > int(meam_extract.size()))
        throw FieldError(base_name(spec.file2) + " names element " + std::to_string(top) + ", but " + std::to_string(meam_extract.size()) +
                         " library entries are read: give every entry it refers to, in its order (SiC.meam: Si=Si C=C)");
    }
  } else if (kBrenner.count(st)) {   // CH.airebo, CH.airebo-m, CH.rebo: carbon and hydrogen
    std::string other;
    for (const auto& e : els) if (e != "C" && e != "H") other += (other.empty() ? "" : ", ") + e;
    if (!other.empty()) throw FieldError(st + " covers carbon and hydrogen only, not " + other);
    // the Brenner tables: a comment block, then the numbers; a Tersoff / setfl file here is a mistake
    std::string l;
    int numbers = 0;
    while (numbers < 20 && std::getline(in, l)) {
      if (l.empty() || l[0] == '#') continue;
      std::istringstream is(l);
      std::string w;
      if (!(is >> w)) continue;
      if (!is_number(w)) throw FieldError(base_name(spec.file) + " is not a " + st + " file (\"" + w + "\" where its table of numbers begins)");
      ++numbers;
    }
    if (numbers < 20) throw FieldError(base_name(spec.file) + " is not a " + st + " file (too short)");
  } else if (kTriplet.count(st)) {
    const int per = kTriplet.at(st);
    std::vector<std::string> words;
    for (std::string l; std::getline(in, l);) {
      if (const size_t h = l.find('#'); h != std::string::npos) l.resize(h);
      std::istringstream is(l);
      for (std::string w; is >> w;) words.push_back(w);
    }
    if (words.empty() || words.size() % size_t(per) != 0)
      throw FieldError(base_name(spec.file) + " is not a " + st + " file: its entries are " + std::to_string(per) + " words each (three elements, then " +
                       std::to_string(per - 3) + " numbers), and it has " + std::to_string(words.size()) + " words");
    for (size_t k = 0; k < words.size(); k += size_t(per)) {
      for (int j = 0; j < 3; ++j)
        if (is_number(words[k + size_t(j)])) throw FieldError(base_name(spec.file) + " is not a " + st + " file: entry " + std::to_string(k / size_t(per) + 1) + " does not start with three elements");
      for (int j = 3; j < per; ++j)
        if (!is_number(words[k + size_t(j)])) throw FieldError(base_name(spec.file) + " is not a " + st + " file: entry " + std::to_string(k / size_t(per) + 1) + " has \"" + words[k + size_t(j)] + "\" where a number belongs");
      have.insert(words[k] + " " + words[k + 1] + " " + words[k + 2]);
    }
    std::string missing;
    for (const auto& a : els)
      for (const auto& b : els)
        for (const auto& c : els)
          if (!have.count(a + " " + b + " " + c)) missing += (missing.empty() ? "" : ", ") + a + "-" + b + "-" + c;
    if (!missing.empty()) throw FieldError(base_name(spec.file) + " has no " + st + " entry for " + missing + ": LAMMPS needs one for every triplet of the group's elements");
  } else {   // setfl (eam/alloy, eam/fs): three comment lines, then "N El1 … ElN"
    std::string l;
    for (int k = 0; k < 4; ++k) std::getline(in, l);
    std::istringstream is(l);
    int n = 0;
    if (!(is >> n) || n <= 0) throw FieldError(base_name(spec.file) + " is not a setfl file (" + st + "): its fourth line should give the number of elements and their names");
    for (std::string w; is >> w;) have.insert(w);
    std::string missing;
    for (const auto& a : els) if (!have.count(a)) missing += (missing.empty() ? "" : ", ") + a;
    if (!missing.empty()) throw FieldError(base_name(spec.file) + " has no " + st + " functions for " + missing);
  }

  ForceField F;
  F.name = st + " (" + base_name(spec.file) + ")";
  const bool c96 = spec.pair_form == "lj9-6";
  if (!c96 && spec.pair_form != "lj12-6") throw FieldError("cross-pair form " + spec.pair_form + " (lj12-6 or lj9-6)");
  F.pair_form = spec.pair_form;
  F.mixing = c96 ? "sixthpower" : "geometric";   // the pairs among the group's types are zero anyway
  F.type_names = els;
  F.excluded.assign(g.atoms.size(), {});
  for (size_t t = 0; t < els.size(); ++t) {
    double x = 0, d = 0;
    if (!uff_vdw(zs[t], x, d)) throw FieldError("UFF has no van der Waals parameters for " + els[t] + " (the cross pairs): give them explicitly");
    // UFF's x is the minimum: 12-6 σ = x / 2^(1/6); the 9-6 form ε[2(σ/r)⁹ − 3(σ/r)⁶] has its minimum at σ itself
    F.lj.push_back({d, c96 ? x : x / std::pow(2.0, 1.0 / 6)});
  }
  // inside the group the potential does it all: no Lennard-Jones between its types
  for (size_t a = 0; a < els.size(); ++a)
    for (size_t b = a; b < els.size(); ++b) F.pair_override[{int(a), int(b)}] = {0.0, 0.5 * (F.lj[a].sigma + F.lj[b].sigma)};
  for (size_t i = 0; i < g.atoms.size(); ++i) {
    F.type_index.push_back(tix[i]);
    F.atom_type.push_back(els[size_t(tix[i])]);
    F.why.push_back(st + " (" + base_name(spec.file) + "): element " + els[size_t(tix[i])]);
    F.charge.push_back(0.0);
    F.mass.push_back(element(zs[size_t(tix[i])]).mass);
  }
  F.manybody.style = st;
  F.manybody.file = spec.file;
  F.manybody.units = units;
  F.manybody.tagged = !file_units.empty();
  F.manybody.element = els;
  F.manybody.citation = after(head, "CITATION:");
  F.manybody.metal_only = kBrenner.count(st) > 0 || meam;
  F.manybody.file2 = meam ? spec.file2 : "";
  F.manybody.entry = meam_entry;
  F.manybody.extract = meam_extract;
  F.manybody.args = !spec.args.empty() ? spec.args : (st == "airebo" || st == "airebo/morse") ? "3.0 1 1" : "";
  std::vector<std::string> said;
  said.push_back(std::to_string(g.atoms.size()) + " atoms by the " + st + " potential of " + base_name(spec.file) + (F.manybody.citation.empty() ? "" : " (" + F.manybody.citation + ")") +
                 ", energies in " + (F.manybody.metal_only ? "eV (LAMMPS reads it in metal units only: the LAMMPS files are written in metal units, eV, ps, bar)"
                                                        : units == "metal" ? "eV (LAMMPS converts them to kcal/mol)" : "kcal/mol") +
                 "; standard atomic masses, no charges, no bonded terms" + (F.manybody.args.empty() ? "" : "; pair_style " + st + " " + F.manybody.args));
  said.push_back(std::string("cross pairs with the rest: Lennard-Jones ") + (c96 ? "9-6 (the class II form of the other groups, with UFF's well depth D and minimum x) " : "12-6 ") +
                 "from UFF (Rappé et al. 1992) for " + [&] {
    std::string l;
    for (const auto& e : els) l += (l.empty() ? "" : ", ") + e;
    return l;
  }() + " unless given; CAPS does not evaluate " + st + ": in CAPS these atoms are held still, the LAMMPS inputs carry the potential");
  if (g.has_charges) {
    bool charged = false;
    for (const auto& a : g.atoms) charged = charged || a.charge != 0;
    if (charged) said.push_back("the structure's charges on these atoms were dropped: the " + st + " potential carries their interactions");
  }
  for (const auto& s : said) F.notes.push_back(s);
  if (notes) *notes = said;
  return F;
}

std::string manybody_file_name(const ManyBodyFile& mb, const std::string& units) {
  const std::filesystem::path p(mb.file);
  if (units == "real" && kBrenner.count(mb.style)) return p.stem().string() + "-real" + p.extension().string();   // CAPS's converted copy
  if (mb.tagged) return p.filename().string();
  return p.stem().string() + "-" + mb.units + p.extension().string();
}

std::string manybody_file2_name(const ManyBodyFile& mb) { return mb.file2.empty() ? "" : std::filesystem::path(mb.file2).filename().string(); }

std::string write_manybody_file(const ManyBodyFile& mb, const std::string& dir, const std::string& units) {
  const std::string name = manybody_file_name(mb, units);
  const std::filesystem::path out = std::filesystem::path(dir.empty() ? "." : dir) / name;
  std::error_code ec;
  if (!mb.file2.empty()) {   // MEAM's parameter file, as it is
    const std::filesystem::path o2 = std::filesystem::path(dir.empty() ? "." : dir) / manybody_file2_name(mb);
    if (!(std::filesystem::exists(o2) && std::filesystem::equivalent(o2, mb.file2, ec))) std::filesystem::copy_file(mb.file2, o2, std::filesystem::copy_options::overwrite_existing);
  }
  if (std::filesystem::exists(out) && std::filesystem::equivalent(out, mb.file, ec)) return name;   // the inputs are beside the file itself
  std::ifstream in(mb.file, std::ios::binary);
  if (!in) throw std::runtime_error("cannot read the potential file " + mb.file);
  std::stringstream ss;
  ss << in.rdbuf();
  std::string text = ss.str();
  if (units == "real" && kBrenner.count(mb.style)) text = brenner_to_real(text, mb.style);
  else if (!mb.tagged) {   // LAMMPS looks for the units on the first line only
    const size_t eol = text.find('\n');
    std::string first = text.substr(0, eol);
    const std::string add = "UNITS: " + mb.units;
    const bool comment_line = kSetfl.count(mb.style) || (!first.empty() && first.find_first_not_of(" \t") != std::string::npos && first[first.find_first_not_of(" \t")] == '#');
    if (comment_line) text = first + " " + add + (eol == std::string::npos ? "" : text.substr(eol));
    else text = "# " + add + " (the units given in CAPS; the file did not say)\n" + text;
  }
  std::ofstream f(out, std::ios::binary);
  if (!f) throw std::runtime_error("cannot write " + out.string());
  f << text;
  return name;
}

}  // namespace caps
