// CAPS literature many-body potentials for a group of atoms (see manybody.hpp).
#include "caps/manybody.hpp"

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
// read by LAMMPS in metal units only (no conversion to real)
const std::set<std::string> kMetalOnly = {"airebo", "airebo/morse", "rebo", "meam", "meam/c", "bop", "comb", "comb3", "lcbop", "polymorphic", "edip", "extep"};

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

const std::vector<std::string>& manybody_styles() {
  static const std::vector<std::string> s = {"tersoff", "tersoff/mod", "tersoff/mod/c", "tersoff/zbl", "sw", "vashishta", "gw", "gw/zbl", "eam/alloy", "eam/fs"};
  return s;
}

ForceField manybody_part(const System& g, const ManyBodySpec& spec, std::vector<std::string>* notes) {
  const std::string& st = spec.style;
  if (kMetalOnly.count(st))
    throw FieldError("pair style " + st + ": LAMMPS reads its file in metal units only (it does not convert it to real units), and CAPS writes its "
                     "inputs in real units (kcal/mol, Å, fs) for the force field of the rest. Use a potential LAMMPS converts: " +
                     [] {
                       std::string l;
                       for (const auto& x : manybody_styles()) l += (l.empty() ? "" : ", ") + x;
                       return l;
                     }());
  if (!kTriplet.count(st) && !kSetfl.count(st)) {
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
  std::set<std::string> have;   // "A B C" triplets, or elements
  if (kTriplet.count(st)) {
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
  std::vector<std::string> said;
  said.push_back(std::to_string(g.atoms.size()) + " atoms by the " + st + " potential of " + base_name(spec.file) + (F.manybody.citation.empty() ? "" : " (" + F.manybody.citation + ")") +
                 ", energies in " + (units == "metal" ? "eV (LAMMPS converts them to kcal/mol)" : "kcal/mol") + "; standard atomic masses, no charges, no bonded terms");
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

std::string manybody_file_name(const ManyBodyFile& mb) {
  const std::filesystem::path p(mb.file);
  if (mb.tagged) return p.filename().string();
  return p.stem().string() + "-" + mb.units + p.extension().string();
}

std::string write_manybody_file(const ManyBodyFile& mb, const std::string& dir) {
  const std::string name = manybody_file_name(mb);
  const std::filesystem::path out = std::filesystem::path(dir.empty() ? "." : dir) / name;
  std::error_code ec;
  if (std::filesystem::exists(out) && std::filesystem::equivalent(out, mb.file, ec)) return name;   // the inputs are beside the file itself
  std::ifstream in(mb.file, std::ios::binary);
  if (!in) throw std::runtime_error("cannot read the potential file " + mb.file);
  std::stringstream ss;
  ss << in.rdbuf();
  std::string text = ss.str();
  if (!mb.tagged) {   // LAMMPS looks for the units on the first line only
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
