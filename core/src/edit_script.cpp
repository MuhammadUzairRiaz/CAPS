// Written edits (caps edit, Document.edit): one operation per ';' or line, applied in order.
#include <cstdio>
#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "caps/edit.hpp"
#include "caps/elements.hpp"
#include "caps/lattice.hpp"
#include "caps/spacegroup.hpp"

namespace caps {

namespace {

std::vector<std::string> words(const std::string& t) {
  std::istringstream in(t);
  std::vector<std::string> w;
  for (std::string x; in >> x;) w.push_back(x);
  return w;
}

double number(const std::string& t, const char* what) {
  size_t used = 0;
  double v = 0;
  try { v = std::stod(t, &used); } catch (const std::exception&) { used = 0; }
  if (used != t.size() || !std::isfinite(v)) throw std::invalid_argument(std::string("not a number for ") + what + ": '" + t + "'");
  return v;
}

uint32_t atom(const System& s, const std::string& t) {
  const double v = number(t, "an atom");
  if (v < 1 || v != std::floor(v) || v > double(s.atoms.size()))
    throw std::invalid_argument("atom " + t + " does not exist (1 to " + std::to_string(s.atoms.size()) + ")");
  return uint32_t(v) - 1;
}

Vec3 vec(const std::string& t, const char* what) {
  std::vector<double> x;
  std::string cur;
  for (char c : t + ",") {
    if (c == ',') { x.push_back(number(cur, what)); cur.clear(); } else cur += c;
  }
  if (x.size() != 3) throw std::invalid_argument(std::string(what) + " needs x,y,z: '" + t + "'");
  return {x[0], x[1], x[2]};
}

std::vector<uint32_t> indices(const std::vector<char>& m) {
  std::vector<uint32_t> v;
  for (uint32_t i = 0; i < m.size(); ++i) if (m[i]) v.push_back(i);
  return v;
}

}  // namespace

std::vector<char> parse_selection(const System& s, const std::string& sel) {
  const size_t n = s.atoms.size();
  auto starts = [&](const char* p) { return sel.rfind(p, 0) == 0; };
  std::vector<char> m;
  if (sel == "all") m.assign(n, 1);
  else if (starts("element:")) {
    std::string e = sel.substr(8);
    std::replace(e.begin(), e.end(), ',', ' ');
    m = select_element(s, e);
  } else if (starts("smarts:")) m = select_smarts(s, sel.substr(7));
  else if (starts("type:")) m = select_type(s, sel.substr(5));
  else {
    m.assign(n, 0);
    std::string cur;
    for (char c : sel + ",") {
      if (c != ',') { cur += c; continue; }
      if (cur.empty()) continue;
      const auto dash = cur.find('-', 1);
      if (dash == std::string::npos) m[atom(s, cur)] = 1;
      else {
        const uint32_t a = atom(s, cur.substr(0, dash)), b = atom(s, cur.substr(dash + 1));
        if (b < a) throw std::invalid_argument("range " + cur + " runs backwards");
        for (uint32_t i = a; i <= b; ++i) m[i] = 1;
      }
      cur.clear();
    }
  }
  m.resize(n, 0);
  if (std::none_of(m.begin(), m.end(), [](char c) { return c != 0; })) throw std::invalid_argument("'" + sel + "' selects no atoms");
  return m;
}

std::vector<std::string> edit_script(System& s, const std::string& script) {
  std::vector<std::string> ops, report;
  std::string cur;
  for (char c : script + ";") {
    if (c == ';' || c == '\n') {
      // a comment: '#' opening the edit or after a space (SMARTS keeps its [#6])
      for (size_t h = cur.find('#'); h != std::string::npos; h = cur.find('#', h + 1))
        if (h == 0 || cur[h - 1] == ' ' || cur[h - 1] == '\t') { cur.erase(h); break; }
      const auto a = cur.find_first_not_of(" \t\r"), b = cur.find_last_not_of(" \t\r");
      if (a != std::string::npos) ops.push_back(cur.substr(a, b - a + 1));
      cur.clear();
    } else cur += c;
  }
  int k = 0;
  for (const auto& op : ops) {
    ++k;
    const auto w = words(op);
    const std::string& v = w[0];
    auto need = [&](size_t lo, size_t hi, const char* usage) {
      if (w.size() < lo || w.size() > hi) throw std::invalid_argument(std::string("usage: ") + usage);
    };
    try {
      std::string done;
      if (v == "element") {
        need(3, 3, "element SEL Sym");
        const int z = element_from_symbol(w[2]);
        if (z <= 0) throw std::invalid_argument("unknown element " + w[2]);
        const auto sel = indices(parse_selection(s, w[1]));
        for (uint32_t i : sel) set_element(s, i, z);
        done = std::to_string(sel.size()) + " atoms made " + element(z).symbol;
      } else if (v == "delete") {
        need(2, 2, "delete SEL");
        const auto m = parse_selection(s, w[1]);
        const auto n = std::count(m.begin(), m.end(), 1);
        delete_atoms(s, m);
        done = std::to_string(n) + " atoms deleted (the atoms after them are renumbered)";
      } else if (v == "bond") {
        need(3, 4, "bond I J [order]");
        const uint32_t i = atom(s, w[1]), j = atom(s, w[2]);
        if (i == j) throw std::invalid_argument("an atom cannot bond to itself");
        const int order = w.size() > 3 ? int(number(w[3], "the bond order")) : 1;
        if (order < 1 || order > 3) throw std::invalid_argument("bond order 1, 2 or 3");
        add_bond(s, i, j, order);
        done = "bond " + w[1] + "–" + w[2] + " (order " + std::to_string(order) + ")";
      } else if (v == "unbond") {
        need(3, 3, "unbond I J");
        if (!remove_bond(s, atom(s, w[1]), atom(s, w[2]))) throw std::invalid_argument("atoms " + w[1] + " and " + w[2] + " are not bonded");
        done = "bond " + w[1] + "–" + w[2] + " removed";
      } else if (v == "addh") {
        need(1, 2, "addh [SEL]");
        const int n = add_hydrogens(s, w.size() > 1 ? parse_selection(s, w[1]) : std::vector<char>{});
        done = std::to_string(n) + " hydrogens added";
      } else if (v == "attach") {
        need(3, 3, "attach I SMILES");
        const auto added = attach_fragment(s, atom(s, w[1]), w[2]);
        done = std::to_string(added.size()) + " atoms of " + w[2] + " attached to atom " + w[1];
      } else if (v == "length") {
        need(4, 4, "length I J Å");
        const double r = number(w[3], "the length");
        if (r <= 0) throw std::invalid_argument("a length must be positive");
        set_bond_length(s, atom(s, w[1]), atom(s, w[2]), r);
        done = "bond " + w[1] + "–" + w[2] + " = " + w[3] + " Å (the side of atom " + w[2] + " moved)";
      } else if (v == "angle") {
        need(5, 5, "angle I J K °");
        set_bond_angle(s, atom(s, w[1]), atom(s, w[2]), atom(s, w[3]), number(w[4], "the angle"));
        done = "angle " + w[1] + "–" + w[2] + "–" + w[3] + " = " + w[4] + "°";
      } else if (v == "torsion") {
        need(6, 6, "torsion I J K L °");
        set_torsion(s, atom(s, w[1]), atom(s, w[2]), atom(s, w[3]), atom(s, w[4]), number(w[5], "the dihedral"));
        done = "dihedral " + w[1] + "–" + w[2] + "–" + w[3] + "–" + w[4] + " = " + w[5] + "°";
      } else if (v == "invert") {
        need(2, 2, "invert I");
        invert_centre(s, atom(s, w[1]));
        done = "centre " + w[1] + " inverted";
      } else if (v == "config") {
        need(3, 3, "config I R|S");
        if (w[2] != "R" && w[2] != "S") throw std::invalid_argument("config takes R or S");
        if (!set_configuration(s, atom(s, w[1]), w[2])) throw std::invalid_argument("atom " + w[1] + " is not a stereocentre");
        done = "centre " + w[1] + " is " + w[2];
      } else if (v == "rotate") {
        need(4, 4, "rotate SEL x,y,z °");
        const Vec3 ax = vec(w[2], "the axis");
        if (norm(ax) == 0) throw std::invalid_argument("the axis is zero");
        const auto sel = indices(parse_selection(s, w[1]));
        rotate_atoms(s, sel, ax, number(w[3], "the angle"));
        done = std::to_string(sel.size()) + " atoms rotated " + w[3] + "° about their centre";
      } else if (v == "mirror") {
        need(3, 3, "mirror SEL nx,ny,nz");
        const Vec3 nrm = vec(w[2], "the plane normal");
        if (norm(nrm) == 0) throw std::invalid_argument("the normal is zero");
        const auto sel = indices(parse_selection(s, w[1]));
        mirror_atoms(s, sel, nrm);
        done = std::to_string(sel.size()) + " atoms mirrored (their stereocentres inverted)";
      } else if (v == "move") {
        need(3, 3, "move SEL dx,dy,dz");
        const Vec3 d = vec(w[2], "the shift");
        const auto sel = indices(parse_selection(s, w[1]));
        for (uint32_t i : sel) s.atoms[i].pos = s.atoms[i].pos + d;
        done = std::to_string(sel.size()) + " atoms moved by (" + w[2] + ") Å";
      } else if (v == "clean") {
        need(1, 2, "clean [SEL]");
        clean_up(s, w.size() > 1 ? parse_selection(s, w[1]) : std::vector<char>{});
        done = w.size() > 1 ? "the selection cleaned up (UFF, the rest held)" : "cleaned up (UFF)";
      } else if (v == "tacticity") {
        need(2, 2, "tacticity iso|syndio");
        if (w[1] != "iso" && w[1] != "syndio") throw std::invalid_argument("tacticity iso or syndio");
        const int n = set_tacticity(s, w[1] == "iso");
        done = std::to_string(n) + " centres inverted: " + (w[1] == "iso" ? "isotactic" : "syndiotactic");
      } else if (v == "supercell") {
        need(4, 4, "supercell na nb nc");
        const int na = int(number(w[1], "na")), nb = int(number(w[2], "nb")), nc = int(number(w[3], "nc"));
        if (na < 1 || nb < 1 || nc < 1) throw std::invalid_argument("repeats must be 1 or more");
        s = supercell(s, na, nb, nc);
        done = std::to_string(s.atoms.size()) + " atoms";
      } else if (v == "primitive") {
        need(1, 2, "primitive [tolerance Å]");
        int kk = 1;
        s = find_primitive_cell(s, w.size() > 1 ? number(w[1], "the tolerance") : 0.1, &kk);
        done = kk == 1 ? "already primitive" : std::to_string(s.atoms.size()) + " atoms (" + std::to_string(kk) + " lattice points folded)";
      } else if (v == "niggli") {
        need(1, 1, "niggli");
        NiggliResult r;
        s = niggli_cell(s, &r);
        char b[140];
        std::snprintf(b, sizeof b, "a %.4f b %.4f c %.4f Å, α %.2f β %.2f γ %.2f°", r.a, r.b, r.c, r.alpha, r.beta, r.gamma);
        done = b;
      } else if (v == "conventional") {
        need(1, 2, "conventional [tolerance Å]");
        ConventionalResult r;
        s = conventional_cell(s, w.size() > 1 ? number(w[1], "the tolerance") : 0.1, &r);
        done = r.hm + " (No. " + std::to_string(r.number) + "), " + std::to_string(s.atoms.size()) + " atoms";
      } else if (v == "redefine") {
        need(10, 10, "redefine m11 m12 m13 m21 m22 m23 m31 m32 m33 (column j: the new vector j in the old ones)");
        Mat3 m{};
        for (int q = 0; q < 9; ++q) m[size_t(q / 3)][size_t(q % 3)] = number(w[size_t(q + 1)], "a matrix entry");
        s = transform_cell(s, m);
        done = std::to_string(s.atoms.size()) + " atoms";
      } else if (v == "vacuum") {
        need(2, 2, "vacuum Å");
        SlabResult r;
        s = vacuum_slab(s, number(w[1], "the vacuum"), true, &r);
        char b[96];
        std::snprintf(b, sizeof b, "%.2f Å slab, %.2f Å vacuum along c", r.thickness, r.vacuum);
        done = b;
      } else if (v == "nanowire") {
        need(5, 8, "nanowire u v w radius [repeats] [cylinder|hexagonal|square] [vacuum]");
        WireOptions o;
        o.uvw = {int(number(w[1], "u")), int(number(w[2], "v")), int(number(w[3], "w"))};
        o.radius = number(w[4], "the radius");
        if (w.size() > 5) o.repeats = int(number(w[5], "the repeats"));
        if (w.size() > 6) o.shape = w[6];
        if (w.size() > 7) o.vacuum = number(w[7], "the vacuum");
        s = nanowire(s, o);
        done = std::to_string(s.atoms.size()) + " atoms";
      } else {
        throw std::invalid_argument("unknown edit '" + v + "' (element, delete, bond, unbond, addh, attach, length, angle, torsion, invert, config, rotate, mirror, move, clean, tacticity, supercell, primitive, niggli, conventional, redefine, vacuum, nanowire)");
      }
      report.push_back(std::to_string(k) + ". " + op + " → " + done);
    } catch (const std::exception& e) {
      throw std::invalid_argument("edit " + std::to_string(k) + " (" + op + "): " + e.what());
    }
  }
  return report;
}

}  // namespace caps
