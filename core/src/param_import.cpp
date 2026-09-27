// Parameters from other programs' files, as CAPS force-field rules to merge over a library force field: AMBER frcmod
// (and the parameter sections of a .dat) and the GROMACS [ *types ] sections of a .itp / .top. Units become LAMMPS
// real (kcal/mol, Å, degrees) in the CAPS styles: bonds and angles harmonic K (x − x0)², dihedrals fourier, impropers
// cvff / fourier / harmonic, pairs ε σ.
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <map>
#include <sstream>

#include "caps/elements.hpp"
#include "caps/ffdef.hpp"

namespace caps {

namespace {

std::string trim_s(const std::string& s) {
  const auto b = s.find_first_not_of(" \t\r\n");
  if (b == std::string::npos) return "";
  return s.substr(b, s.find_last_not_of(" \t\r\n") - b + 1);
}

std::vector<std::string> words(const std::string& s) {
  std::istringstream in(s);
  std::vector<std::string> w;
  for (std::string t; in >> t;) w.push_back(t);
  return w;
}

double to_num(const std::string& s, const std::string& what, const std::string& line) {
  try {
    size_t k = 0;
    const double v = std::stod(s, &k);
    if (k != s.size()) throw 0;
    return v;
  } catch (...) {
    throw FFError("could not read " + what + " from '" + trim_s(line) + "'");
  }
}

std::string amber_type(std::string t) {
  t = trim_s(t);
  return t == "X" ? "*" : t;
}

// "c3-c3-ca" (types of up to two characters, may carry spaces) → the types, then the rest of the line
std::pair<std::vector<std::string>, std::string> amber_key(const std::string& line, int n) {
  std::vector<std::string> types;
  size_t pos = 0;
  for (int k = 0; k < n; ++k) {
    const size_t dash = k + 1 < n ? line.find('-', pos) : std::string::npos;
    if (k + 1 < n && dash == std::string::npos) throw FFError("expected " + std::to_string(n) + " types in '" + trim_s(line) + "'");
    if (k + 1 < n) {
      types.push_back(amber_type(line.substr(pos, dash - pos)));
      pos = dash + 1;
    } else {
      // the last type: up to two characters (then the numbers)
      const size_t end = std::min(line.size(), pos + 2);
      types.push_back(amber_type(line.substr(pos, end - pos)));
      pos = end;
    }
  }
  return {types, pos < line.size() ? line.substr(pos) : std::string()};
}

FFRule rule(const std::string& name, const std::vector<std::string>& match, const std::string& style, std::vector<double> params,
            const std::string& comment = "") {
  FFRule r;
  r.name = name;
  r.match = match;
  r.style = style;
  r.params = std::move(params);
  r.comment = comment;
  return r;
}

bool is_int(const std::string& s) { return !s.empty() && s.find_first_not_of("0123456789") == std::string::npos; }
bool is_number(const std::string& s) {
  char* end = nullptr;
  std::strtod(s.c_str(), &end);
  return end && *end == 0 && !s.empty();
}

std::string join(const std::vector<std::string>& t) {
  std::string r;
  for (const auto& x : t) r += (r.empty() ? "" : "-") + (x == "*" ? std::string("X") : x);
  return r;
}

}  // namespace

FFDef import_frcmod(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw FFError("cannot open " + path);
  FFDef d;
  d.name = "frcmod";
  d.source = path;
  d.dihedral_style = "fourier";
  d.improper_style = "cvff";
  std::string section;
  std::map<std::string, size_t> open_dihedral;   // a negative periodicity continues the same dihedral
  std::string line;
  std::getline(in, line);   // title
  while (std::getline(in, line)) {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    const std::string t = trim_s(line);
    const std::string head = t.substr(0, 4);
    if (head == "MASS" || head == "BOND" || head == "ANGL" || head == "DIHE" || head == "IMPR" || head == "NONB" || head == "HBON" || head == "CMAP" ||
        head == "IPOL") {
      section = head;
      continue;
    }
    if (t.empty() || t == "END") {
      if (section != "NONB") section.clear();
      continue;
    }
    if (section == "MASS") {
      const auto w = words(line);
      if (w.size() < 2) continue;
      FFType ty;
      ty.name = w[0];
      ty.mass = to_num(w[1], "a mass", line);
      ty.element = element_from_mass(ty.mass, 0.6);
      ty.description = "from " + path.substr(path.find_last_of("/\\") + 1);
      d.types.push_back(ty);
    } else if (section == "BOND") {
      auto [ty, rest] = amber_key(line, 2);
      const auto w = words(rest);
      if (w.size() < 2) throw FFError("a BOND line needs K and r0: '" + t + "'");
      d.bonds.push_back(rule(join(ty), ty, "harmonic", {to_num(w[0], "K", line), to_num(w[1], "r0", line)}));
    } else if (section == "ANGL") {
      auto [ty, rest] = amber_key(line, 3);
      const auto w = words(rest);
      if (w.size() < 2) throw FFError("an ANGL line needs K and θ0: '" + t + "'");
      d.angles.push_back(rule(join(ty), ty, "harmonic", {to_num(w[0], "K", line), to_num(w[1], "theta0", line)}, ""));
    } else if (section == "DIHE") {
      auto [ty, rest] = amber_key(line, 4);
      const auto w = words(rest);
      if (w.size() < 4) throw FFError("a DIHE line needs IDIVF, PK, phase and PN: '" + t + "'");
      const double idivf = to_num(w[0], "IDIVF", line), pk = to_num(w[1], "PK", line), phase = to_num(w[2], "phase", line),
                   pn = to_num(w[3], "PN", line);
      const double k = pk / (idivf > 0 ? idivf : 1.0);
      const std::string key = join(ty);
      auto it = open_dihedral.find(key);
      if (it != open_dihedral.end()) {   // the previous line of this dihedral had a negative periodicity
        auto& p = d.dihedrals[it->second].params;
        p[0] += 1;
        p.insert(p.end(), {k, std::fabs(pn), phase});
      } else {
        d.dihedrals.push_back(rule(key, ty, "fourier", {1, k, std::fabs(pn), phase}, ""));
      }
      if (pn < 0) open_dihedral[key] = it != open_dihedral.end() ? it->second : d.dihedrals.size() - 1;
      else open_dihedral.erase(key);
    } else if (section == "IMPR") {
      auto [ty, rest] = amber_key(line, 4);
      const auto w = words(rest);
      if (w.size() < 3) throw FFError("an IMPR line needs PK, phase and PN: '" + t + "'");
      const double pk = to_num(w[0], "PK", line), phase = to_num(w[1], "phase", line), pn = to_num(w[2], "PN", line);
      if (std::fabs(phase) > 1e-6 && std::fabs(std::fabs(phase) - 180) > 1e-6)
        d.impropers.push_back(rule(join(ty), ty, "fourier", {1, pk, std::fabs(pn), phase}, ""));
      else
        d.impropers.push_back(rule(join(ty), ty, "cvff", {pk, std::fabs(phase) > 90 ? -1.0 : 1.0, std::fabs(pn)}, ""));
    } else if (section == "NONB") {
      const auto w = words(line);
      if (w.size() < 3) continue;
      const double rstar = to_num(w[1], "R*", line), eps = to_num(w[2], "epsilon", line);
      d.pairs.push_back(rule(w[0], {w[0]}, "", {eps, 2 * rstar / std::pow(2.0, 1.0 / 6.0)}, "R*/2 " + w[1]));
    }
    // HBON (10-12 terms), CMAP and IPOL are not read
  }
  if (d.bonds.empty() && d.angles.empty() && d.dihedrals.empty() && d.impropers.empty() && d.pairs.empty() && d.types.empty())
    throw FFError(path + ": no MASS, BOND, ANGL, DIHE, IMPR or NONB parameters found (not an AMBER frcmod?)");
  return d;
}

FFDef import_gromacs_params(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw FFError("cannot open " + path);
  FFDef d;
  d.name = "gromacs";
  d.source = path;
  constexpr double kJ = 1 / 4.184;
  int comb = 2;
  std::string section;
  std::vector<std::string> skipped;
  for (std::string line; std::getline(in, line);) {
    const auto sc = line.find(';');
    std::string t = trim_s(sc == std::string::npos ? line : line.substr(0, sc));
    if (t.empty() || t[0] == '#') continue;
    if (t.front() == '[') {
      section = trim_s(t.substr(1, t.find(']') - 1));
      continue;
    }
    const auto w = words(t);
    if (section == "defaults" && w.size() >= 2) comb = int(to_num(w[1], "the combination rule", line));
    else if (section == "atomtypes" && w.size() >= 6) {
      // name [bond_type] [at.num] mass charge ptype V W
      // the particle type (A S V D) is always third from the end: σ and ε (or C6 C12) follow it
      const size_t p = w.size() - 3;
      if (p < 3 || !(w[p] == "A" || w[p] == "S" || w[p] == "V" || w[p] == "D")) throw FFError("could not read the atom type '" + t + "'");
      const double mass = to_num(w[p - 2], "a mass", line), v = to_num(w[p + 1], "V", line), wv = to_num(w[p + 2], "W", line);
      FFType ty;
      ty.name = w[0];
      ty.mass = mass;
      ty.element = p >= 4 && is_int(w[p - 3]) ? int(to_num(w[p - 3], "the atomic number", line)) : element_from_mass(mass, 0.6);
      ty.description = "from " + path.substr(path.find_last_of("/\\") + 1);
      d.types.push_back(ty);
      double sigma, eps;
      if (comb == 1) {   // C6, C12
        sigma = v > 0 && wv > 0 ? std::pow(wv / v, 1.0 / 6.0) : 0;
        eps = wv > 0 ? v * v / (4 * wv) : 0;
      } else {
        sigma = v;
        eps = wv;
      }
      d.pairs.push_back(rule(w[0], {w[0]}, "", {eps * kJ, sigma * 10}, ""));
    } else if (section == "bondtypes" && w.size() >= 5) {
      const int fn = int(to_num(w[2], "the function", line));
      if (fn != 1) { skipped.push_back("bond function " + std::to_string(fn) + " (" + w[0] + "-" + w[1] + ")"); continue; }
      d.bonds.push_back(rule(w[0] + "-" + w[1], {w[0], w[1]}, "harmonic", {0.5 * to_num(w[4], "kb", line) * kJ / 100, to_num(w[3], "b0", line) * 10}, ""));
    } else if (section == "angletypes" && w.size() >= 6) {
      const int fn = int(to_num(w[3], "the function", line));
      if (fn != 1 && fn != 5) { skipped.push_back("angle function " + std::to_string(fn) + " (" + w[0] + "-" + w[1] + "-" + w[2] + ")"); continue; }
      d.angles.push_back(rule(w[0] + "-" + w[1] + "-" + w[2], {w[0], w[1], w[2]}, "harmonic", {0.5 * to_num(w[5], "ktheta", line) * kJ, to_num(w[4], "theta0", line)}, ""));
      if (fn == 5 && w.size() >= 8 && to_num(w[7], "kub", line) != 0) skipped.push_back("Urey–Bradley part of " + w[0] + "-" + w[1] + "-" + w[2]);
    } else if (section == "dihedraltypes" && w.size() >= 5) {
      // four types, or two (the central pair) in older files
      const bool four = w.size() >= 6 && is_int(w[4]) && !is_number(w[3]);
      if (!four && !is_int(w[2])) throw FFError("could not read the dihedral type '" + t + "'");
      const size_t nt = four ? 4 : 2;
      const int fn = int(to_num(w[nt], "the function", line));
      std::vector<std::string> ty(w.begin(), w.begin() + long(nt));
      for (auto& x : ty) if (x == "X") x = "*";
      if (nt == 2) ty = {"*", ty[0], ty[1], "*"};
      const auto p = [&](size_t k) { return to_num(w.at(nt + 1 + k), "a dihedral parameter", line); };
      if (fn == 1 || fn == 9) {
        const std::string key = join(ty);
        // function 9 lists several terms for the same types: one fourier rule
        auto it = std::find_if(d.dihedrals.begin(), d.dihedrals.end(), [&](const FFRule& r) { return r.name == key; });
        if (fn == 9 && it != d.dihedrals.end()) {
          it->params[0] += 1;
          it->params.insert(it->params.end(), {p(1) * kJ, p(2), p(0)});
        } else {
          d.dihedrals.push_back(rule(key, ty, "fourier", {1, p(1) * kJ, p(2), p(0)}, ""));
        }
      } else if (fn == 4) {
        d.impropers.push_back(rule(join(ty), ty, "fourier", {1, p(1) * kJ, p(2), p(0)}, "GROMACS periodic improper"));
      } else if (fn == 2) {
        d.impropers.push_back(rule(join(ty), ty, "harmonic", {0.5 * p(1) * kJ, p(0)}, "GROMACS harmonic improper"));
      } else {
        skipped.push_back("dihedral function " + std::to_string(fn) + " (" + join(ty) + ")");
      }
    }
  }
  if (d.bonds.empty() && d.angles.empty() && d.dihedrals.empty() && d.impropers.empty() && d.pairs.empty())
    throw FFError(path + ": no [ atomtypes ], [ bondtypes ], [ angletypes ] or [ dihedraltypes ] parameters found");
  if (!skipped.empty()) {
    std::string s = "not imported (no CAPS form yet): ";
    for (size_t k = 0; k < skipped.size() && k < 6; ++k) s += (k ? ", " : "") + skipped[k];
    if (skipped.size() > 6) s += " and " + std::to_string(skipped.size() - 6) + " more";
    d.references.push_back(s);
  }
  return d;
}

}  // namespace caps
