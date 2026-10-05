// GROMACS topologies (.top / .itp) read as a CAPS structure with its explicit topology: every molecule's atoms (named
// by their types, with charges, residue names and numbers, masses), bonds, constraints (stiff bonds), angles,
// dihedrals, virtual sites and exclusions, repeated as [ molecules ] lists them. The non-bonded terms are not read:
// they come from the CAPS force field the structure is typed with (atom types by name). Coordinates come from a
// separate file (.gro, .pdb, a data file) with the same atoms in the same order.
#include <cctype>
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>

#include "caps/elements.hpp"
#include "caps/io.hpp"

namespace caps {

namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr double kDeg = kPi / 180;
constexpr double kKJ = 4.184;

struct MolType {
  std::string name;
  int nrexcl = 1;
  struct A { std::string type, resname, atomname; long resnr = 0; double charge = 0, mass = std::nan(""); };
  std::vector<A> atoms;
  ExplicitTopology topo;   // atom indices within the molecule
  std::vector<std::pair<uint32_t, uint32_t>> chem;   // bonds and constraints that count for nrexcl
};

std::vector<std::string> words(const std::string& l) {
  std::istringstream is(l);
  std::vector<std::string> w;
  for (std::string x; is >> x;) w.push_back(x);
  return w;
}

// the file with #include resolved and #ifdef / #ifndef / #else / #endif applied; includes that are not found
// (force-field files such as martini_v3.0.0.itp: their non-bonded terms are not read) are noted and skipped
void preprocess(const std::filesystem::path& file, std::set<std::string>& defines, std::vector<std::string>& out, std::vector<std::string>& notes, int depth) {
  if (depth > 20) throw ReadError(file.string() + ": #include nested too deeply");
  std::ifstream in(file);
  if (!in) throw ReadError("cannot open " + file.string());
  std::vector<bool> active;
  std::string line;
  while (std::getline(in, line)) {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    const std::string code = line.substr(0, line.find(';'));
    const auto w = words(code);
    const bool on = std::find(active.begin(), active.end(), false) == active.end();
    if (!w.empty() && w[0][0] == '#') {
      if (w[0] == "#ifdef" || w[0] == "#ifndef") {
        const bool def = w.size() > 1 && defines.count(w[1]);
        active.push_back(w[0] == "#ifdef" ? def : !def);
      } else if (w[0] == "#else") {
        if (active.empty()) throw ReadError(file.string() + ": #else without #ifdef");
        active.back() = !active.back();
      } else if (w[0] == "#endif") {
        if (active.empty()) throw ReadError(file.string() + ": #endif without #ifdef");
        active.pop_back();
      } else if (on && w[0] == "#define" && w.size() > 1) {
        defines.insert(w[1]);
      } else if (on && w[0] == "#include" && w.size() > 1) {
        std::string inc = w[1];
        inc.erase(std::remove(inc.begin(), inc.end(), '"'), inc.end());
        const std::filesystem::path p = file.parent_path() / inc;
        if (std::filesystem::exists(p)) preprocess(p, defines, out, notes, depth + 1);
        else notes.push_back("#include " + inc + " not found: skipped (non-bonded terms come from the CAPS force field)");
      }
      continue;
    }
    if (on) out.push_back(code);
  }
}

double num(const std::vector<std::string>& w, size_t k, const std::string& where) {
  if (k >= w.size()) throw ReadError(where + ": missing parameter");
  return std::stod(w[k]);
}

}  // namespace

System read_gromacs_topology(const std::string& path, std::vector<std::string>* notes_out) {
  std::vector<std::string> lines, notes;
  std::set<std::string> defines;
  preprocess(std::filesystem::path(path), defines, lines, notes, 0);
  std::map<std::string, MolType> types;
  std::vector<std::pair<std::string, long>> molecules;
  std::string sec, title;
  MolType* cur = nullptr;
  std::set<std::string> ignored;
  for (size_t ln = 0; ln < lines.size(); ++ln) {
    const auto w = words(lines[ln]);
    if (w.empty()) continue;
    if (w[0][0] == '[') {
      std::string s = lines[ln];
      s.erase(std::remove_if(s.begin(), s.end(), [](char c) { return c == '[' || c == ']' || std::isspace(static_cast<unsigned char>(c)); }), s.end());
      sec = s;
      continue;
    }
    const std::string where = path + " [ " + sec + " ]";
    auto idx = [&](size_t k) {
      const long i = std::stol(w.at(k)) - 1;
      if (!cur || i < 0 || size_t(i) >= cur->atoms.size()) throw ReadError(where + ": atom " + w.at(k) + " out of range");
      return uint32_t(i);
    };
    if (sec == "moleculetype") {
      MolType m;
      m.name = w[0];
      m.nrexcl = w.size() > 1 ? std::stoi(w[1]) : 1;
      types[m.name] = m;
      cur = &types[m.name];
    } else if (sec == "atoms") {
      if (!cur) throw ReadError(where + ": atoms outside a moleculetype");
      MolType::A a;
      if (w.size() < 5) throw ReadError(where + ": short atom line");
      a.type = w[1];
      a.resnr = std::stol(w[2]);
      a.resname = w[3];
      a.atomname = w[4];
      if (w.size() > 6) a.charge = std::stod(w[6]);
      if (w.size() > 7) a.mass = std::stod(w[7]);
      cur->atoms.push_back(a);
    } else if (sec == "bonds") {
      const int f = int(num(w, 2, where));
      if (f != 1 && f != 6) throw ReadError(where + ": bond function " + w[2] + " is not handled");
      cur->topo.bonds.push_back({idx(0), idx(1), num(w, 4, where) / (2 * kKJ * 100), num(w, 3, where) * 10, f == 6 ? "bond (type 6)" : "bond"});
      if (f == 1) cur->chem.push_back({idx(0), idx(1)});
    } else if (sec == "constraints") {
      cur->topo.bonds.push_back({idx(0), idx(1), 1e6 / (2 * kKJ * 100), num(w, 3, where) * 10, "constraint"});
      if (int(num(w, 2, where)) == 1) cur->chem.push_back({idx(0), idx(1)});
    } else if (sec == "angles") {
      const int f = int(num(w, 3, where));
      const int form = f == 1 ? 0 : f == 2 ? 1 : f == 10 ? 5 : -1;
      if (form < 0) throw ReadError(where + ": angle function " + w[3] + " is not handled");
      cur->topo.angles.push_back({idx(0), idx(1), idx(2), form, num(w, 5, where) / (2 * kKJ), num(w, 4, where) * kDeg, "angle"});
    } else if (sec == "dihedrals") {
      const int f = int(num(w, 4, where));
      if (f == 1 || f == 9 || f == 4)
        cur->topo.dihedrals.push_back({idx(0), idx(1), idx(2), idx(3), f, num(w, 6, where) / kKJ, num(w, 5, where) * kDeg, int(num(w, 7, where)), "dihedral"});
      else if (f == 2)
        cur->topo.dihedrals.push_back({idx(0), idx(1), idx(2), idx(3), 2, num(w, 6, where) / (2 * kKJ), num(w, 5, where) * kDeg, 0, "improper"});
      else if (f == 11) {   // combined bending–torsion: k a0 … a4
        ExplicitTopology::Dihedral d{idx(0), idx(1), idx(2), idx(3), 11, 0, 0, 0, "bending-torsion"};
        for (int n = 0; n < 5; ++n) d.c[size_t(n)] = num(w, 5, where) * num(w, size_t(6 + n), where) / kKJ;
        cur->topo.dihedrals.push_back(d);
      }
      else throw ReadError(where + ": dihedral function " + w[4] + " is not handled");
    } else if (sec == "exclusions") {
      for (size_t k = 1; k < w.size(); ++k) cur->topo.exclusions.push_back({idx(0), idx(k)});
    } else if (sec == "virtual_sitesn") {
      ExplicitTopology::VSite v;
      v.site = idx(0);
      const int f = int(num(w, 1, where));
      if (f == 3) {
        for (size_t k = 2; k + 1 < w.size(); k += 2) v.from.push_back(idx(k)), v.w.push_back(std::stod(w[k + 1]));
      } else {
        for (size_t k = 2; k < w.size(); ++k) v.from.push_back(idx(k));
        if (f == 1) v.w.assign(v.from.size(), 1.0);
        else if (f != 2) throw ReadError(where + ": virtual_sitesn function " + w[1] + " is not handled");
      }
      if (!v.w.empty()) {   // centre of geometry / weighted centre: weights as fractions
        double sum = 0;
        for (double x : v.w) sum += x;
        if (sum <= 0) throw ReadError(where + ": virtual site weights must be positive");
        for (double& x : v.w) x /= sum;
      }
      cur->topo.vsites.push_back(v);
    } else if (sec == "virtual_sites2" || sec == "virtual_sites3") {
      ExplicitTopology::VSite v;
      v.site = idx(0);
      const bool three = sec == "virtual_sites3";
      const int f = int(num(w, three ? 4 : 3, where));
      if (three && f == 4) {   // 3out: x_i + a r_ij + b r_ik + c (r_ij × r_ik), c in 1/nm
        const double a = num(w, 5, where), b = num(w, 6, where), c = num(w, 7, where);
        v.from = {idx(1), idx(2), idx(3)};
        v.w = {1 - a - b, a, b};
        v.c = c / 10;
        cur->topo.vsites.push_back(v);
        continue;
      }
      if (f != 1) throw ReadError(where + ": only linear (function 1) and out-of-plane (virtual_sites3 function 4) sites are handled");
      if (three) {
        const double a = num(w, 5, where), b = num(w, 6, where);
        v.from = {idx(1), idx(2), idx(3)};
        v.w = {1 - a - b, a, b};
      } else {
        const double a = num(w, 4, where);
        v.from = {idx(1), idx(2)};
        v.w = {1 - a, a};
      }
      cur->topo.vsites.push_back(v);
    } else if (sec == "system") {
      title += (title.empty() ? "" : " ") + lines[ln];
    } else if (sec == "molecules") {
      molecules.push_back({w[0], std::stol(w.at(1))});
    } else if (sec == "pairs" && w.size() >= 5 && int(num(w, 2, where)) == 1) {   // with its own σ, ε
      cur->topo.pairs.push_back({idx(0), idx(1), num(w, 4, where) / kKJ, num(w, 3, where) * 10});
    } else if (sec == "pairs" || sec == "cmap" || sec == "settles" || sec == "virtual_sites4") {
      if (!ignored.count(sec)) notes.push_back("[ " + sec + " ] not read");
      ignored.insert(sec);
    }
    // defaults, atomtypes, nonbond_params, pairtypes, position_restraints ...: the force field's or not needed
  }
  if (molecules.empty()) {   // an .itp alone: one copy of each molecule type
    for (const auto& [n, t] : types) molecules.push_back({n, 1});
    notes.push_back("no [ molecules ]: one copy of each molecule type");
  }
  System s;
  s.title = title;
  s.source_format = "gromacs-top";
  auto topo = std::make_shared<ExplicitTopology>();
  topo->source = "GROMACS topology " + std::filesystem::path(path).filename().string();
  std::map<std::string, int> tid;
  int64_t mol = 0;
  bool any_mass = false;
  for (const auto& [name, count] : molecules) {
    auto it = types.find(name);
    if (it == types.end()) throw ReadError(path + ": [ molecules ] names " + name + ", which no [ moleculetype ] defines");
    const MolType& m = it->second;
    // exclusions GROMACS makes from nrexcl: atoms within nrexcl chemical bonds (1-2 are excluded by the bonds themselves)
    std::vector<std::pair<uint32_t, uint32_t>> gen;
    if (m.nrexcl >= 2) {
      std::vector<std::vector<uint32_t>> nb(m.atoms.size());
      for (const auto& [a, b] : m.chem) nb[a].push_back(b), nb[b].push_back(a);
      for (uint32_t a = 0; a < m.atoms.size(); ++a) {
        std::map<uint32_t, int> d{{a, 0}};
        std::vector<uint32_t> front{a};
        for (int step = 1; step <= m.nrexcl; ++step) {
          std::vector<uint32_t> next;
          for (uint32_t u : front)
            for (uint32_t v : nb[u])
              if (!d.count(v)) d[v] = step, next.push_back(v);
          front = next;
        }
        for (const auto& [b, dist] : d)
          if (b > a && dist >= 2) gen.push_back({a, b});
      }
    }
    for (long c = 0; c < count; ++c) {
      ++mol;
      const uint32_t off = uint32_t(s.atoms.size());
      for (const auto& a : m.atoms) {
        Atom at;
        at.id = int64_t(s.atoms.size() + 1);
        at.mol = mol;
        at.name = a.type;
        at.resname = a.resname;
        at.resid = a.resnr;
        at.charge = a.charge;
        // an element when the mass is an element's and the type or name starts with its symbol; beads stay 0
        for (int z = 1; z <= 54 && !std::isnan(a.mass); ++z) {
          const auto& e = element(z);
          const std::string sym = e.symbol;
          if (std::fabs(e.mass - a.mass) < 0.02 && (a.atomname.rfind(sym, 0) == 0 || a.type.rfind(sym, 0) == 0)) { at.element = z; break; }
        }
        auto [ti, fresh] = tid.emplace(a.type, int(tid.size()) + 1);
        at.type = ti->second;
        if (fresh) s.types.push_back({ti->second, std::isnan(a.mass) ? 0.0 : a.mass, a.type});
        s.atoms.push_back(at);
        topo->masses.push_back(a.mass);
        any_mass = any_mass || !std::isnan(a.mass);
      }
      auto sh = [&](uint32_t i) { return i + off; };
      for (auto b : m.topo.bonds) {
        b.i = sh(b.i), b.j = sh(b.j);
        topo->bonds.push_back(b);
        s.bonds.push_back({b.i, b.j, 1});
      }
      for (auto a : m.topo.angles) a.i = sh(a.i), a.j = sh(a.j), a.k = sh(a.k), topo->angles.push_back(a);
      for (auto d : m.topo.dihedrals) d.i = sh(d.i), d.j = sh(d.j), d.k = sh(d.k), d.l = sh(d.l), topo->dihedrals.push_back(d);
      for (const auto& e : m.topo.exclusions) topo->exclusions.push_back({sh(e.first), sh(e.second)});
      for (auto pr : m.topo.pairs) pr.i = sh(pr.i), pr.j = sh(pr.j), topo->pairs.push_back(pr);
      for (const auto& e : gen) topo->exclusions.push_back({sh(e.first), sh(e.second)});
      for (auto v : m.topo.vsites) {
        v.site = sh(v.site);
        for (auto& f : v.from) f = sh(f);
        topo->vsites.push_back(v);
      }
    }
  }
  if (!any_mass) topo->masses.clear();
  topo->natoms = s.atoms.size();
  s.topology = topo;
  s.bonds_from_file = true;
  s.has_charges = true;
  s.has_mol = true;
  if (notes_out) *notes_out = notes;
  s.notes = notes;
  return s;
}

}  // namespace caps
