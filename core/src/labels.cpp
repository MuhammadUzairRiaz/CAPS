#include "caps/labels.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <stdexcept>

#include "caps/analysis.hpp"
#include "caps/appearance.hpp"
#include "caps/elements.hpp"
#include "caps/field.hpp"
#include "caps/typing.hpp"

namespace caps {

const std::vector<LabelKind>& atom_label_kinds() {
  static const std::vector<LabelKind> k = {
      {"element", "Element symbol", "Identity"}, {"element_name", "Element name", "Identity"}, {"atomic_number", "Atomic number", "Identity"},
      {"index", "Index", "Identity"}, {"id", "File id", "Identity"}, {"name", "Atom name", "Identity"},
      {"type", "Force-field type", "Identity"}, {"type_number", "Type number", "Identity"},
      {"molecule", "Molecule", "Identity"}, {"molecule_formula", "Molecule composition (on its first atom)", "Identity"},
      {"chain_end", "Chain ends (the molecule's number at both ends of its backbone)", "Identity"},
      {"residue", "Residue", "Identity"}, {"colour", "Display colour", "Identity"},
      {"charge", "Partial charge", "Chemistry"}, {"formal_charge", "Formal charge", "Chemistry"},
      {"oxidation_state", "Oxidation state", "Chemistry"}, {"hybridisation", "Hybridisation", "Chemistry"},
      {"stereo", "R/S and E/Z", "Chemistry"},
      {"coordination", "Neighbours", "Chemistry"}, {"hydrogens", "Attached hydrogens", "Chemistry"},
      {"valence", "Bond-order sum", "Chemistry"}, {"ring", "Smallest ring", "Chemistry"}, {"aromatic", "Aromatic", "Chemistry"},
      {"electronegativity", "Pauling electronegativity", "Element"}, {"configuration", "Electron configuration", "Element"},
      {"mass", "Mass", "Element"}, {"mass_number", "Mass number (most common isotope)", "Element"},
      {"covalent_radius", "Covalent radius", "Element"}, {"vdw_radius", "van der Waals radius", "Element"},
      {"xyz", "Position x, y, z", "Geometry"}, {"x", "x", "Geometry"}, {"y", "y", "Geometry"}, {"z", "z", "Geometry"},
      {"fractional", "Fractional coordinates", "Geometry"}, {"image", "Periodic image", "Geometry"},
  };
  return k;
}

const std::vector<LabelKind>& bond_label_kinds() {
  static const std::vector<LabelKind> k = {
      {"order", "Bond order", "Chemistry"}, {"order_name", "Bond type (single, double …)", "Chemistry"},
      {"chemical", "Chemical type (C=O, C–H …)", "Chemistry"}, {"ring", "Ring or chain", "Chemistry"},
      {"rotatable", "Rotatable", "Chemistry"}, {"delta_en", "Electronegativity difference", "Chemistry"},
      {"length", "Length", "Geometry"}, {"midpoint", "Midpoint x, y, z", "Geometry"},
      {"index", "Index", "Identity"}, {"name", "Name (atom names)", "Identity"}, {"type", "Force-field type pair", "Identity"},
      {"ff_style", "Force-field style", "Force field"}, {"ff_r0", "Equilibrium length r₀", "Force field"},
      {"ff_k", "Force constant k", "Force field"}, {"energy", "Bond energy", "Force field"},
  };
  return k;
}

namespace {

std::string fmt(const char* f, double v) {
  char b[64];
  std::snprintf(b, sizeof b, f, v);
  return b;
}

std::string signed_int(int v) { return v > 0 ? "+" + std::to_string(v) : std::to_string(v); }

std::string sub_digits(const std::string& f) {   // C6H12O → C₆H₁₂O
  static const char* d[] = {"₀", "₁", "₂", "₃", "₄", "₅", "₆", "₇", "₈", "₉"};
  std::string o;
  for (char c : f) o += (c >= '0' && c <= '9') ? d[c - '0'] : std::string(1, c);
  return o;
}

std::string type_label(const System& s, uint32_t i, const std::vector<std::string>* types) {
  if (types && i < types->size() && !(*types)[i].empty()) return (*types)[i];
  for (const auto& t : s.types) if (t.type == s.atoms[i].type && !t.label.empty()) return t.label;
  return s.atoms[i].type ? std::to_string(s.atoms[i].type) : "";
}

bool metal(int z) {
  return (z >= 3 && z <= 4) || (z >= 11 && z <= 13) || (z >= 19 && z <= 31) || (z >= 37 && z <= 50) || (z >= 55 && z <= 84) || z >= 87;
}

// molecules: the file's ids, else bonded components
std::vector<int64_t> molecules(const System& s) {
  std::vector<int64_t> m(s.atoms.size());
  bool ids = false;
  for (size_t i = 0; i < s.atoms.size(); ++i) { m[i] = s.atoms[i].mol; ids = ids || m[i] != 0; }
  if (ids) return m;
  std::vector<uint32_t> p(s.atoms.size());
  for (uint32_t i = 0; i < p.size(); ++i) p[i] = i;
  auto find = [&](uint32_t x) { while (p[x] != x) x = p[x] = p[p[x]]; return x; };
  for (const auto& b : s.bonds) p[find(b.i)] = find(b.j);
  std::map<uint32_t, int64_t> id;
  for (uint32_t i = 0; i < p.size(); ++i) m[i] = id.emplace(find(i), int64_t(id.size()) + 1).first->second;
  return m;
}

}  // namespace

std::vector<std::string> atom_labels(const System& s, const std::string& kind, const std::vector<std::string>* types) {
  const size_t n = s.atoms.size();
  std::vector<std::string> out(n);
  auto each = [&](auto f) { for (uint32_t i = 0; i < n; ++i) out[i] = f(i, s.atoms[i]); return out; };
  if (kind == "element") return each([](uint32_t, const Atom& a) { return std::string(element(a.element).symbol); });
  if (kind == "element_name") return each([](uint32_t, const Atom& a) { return std::string(element_name(a.element)); });
  if (kind == "atomic_number") return each([](uint32_t, const Atom& a) { return a.element ? std::to_string(a.element) : std::string(); });
  if (kind == "index") return each([](uint32_t i, const Atom&) { return std::to_string(i + 1); });
  if (kind == "id") return each([](uint32_t i, const Atom& a) { return std::to_string(a.id ? a.id : int64_t(i) + 1); });
  if (kind == "name") return each([](uint32_t, const Atom& a) { return a.name; });
  if (kind == "type") return each([&](uint32_t i, const Atom&) { return type_label(s, i, types); });
  if (kind == "type_number") return each([](uint32_t, const Atom& a) { return a.type ? std::to_string(a.type) : std::string(); });
  if (kind == "residue") return each([](uint32_t, const Atom& a) { return a.resname.empty() && !a.resid ? std::string() : a.resname + (a.resid ? std::to_string(a.resid) : ""); });
  if (kind == "colour")
    return each([](uint32_t, const Atom& a) { char b[12]; std::snprintf(b, sizeof b, "#%06X", element(a.element).rgb & 0xFFFFFFu); return std::string(b); });
  if (kind == "charge") return each([](uint32_t, const Atom& a) { return fmt("%+.3f", a.charge); });
  if (kind == "electronegativity") return each([](uint32_t, const Atom& a) { const double e = pauling_electronegativity(a.element); return e > 0 ? fmt("%.2f", e) : std::string(); });
  if (kind == "configuration") return each([](uint32_t, const Atom& a) { return std::string(electron_configuration(a.element)); });
  if (kind == "mass")
    return each([&](uint32_t, const Atom& a) {
      for (const auto& t : s.types) if (t.type == a.type && a.type && t.mass > 0) return fmt("%.3f", t.mass);
      return fmt("%.3f", element(a.element).mass);
    });
  if (kind == "mass_number") return each([](uint32_t, const Atom& a) { const int m = most_common_mass_number(a.element); return m ? std::to_string(m) : std::string(); });
  if (kind == "covalent_radius") return each([](uint32_t, const Atom& a) { return fmt("%.2f Å", element(a.element).covalent); });
  if (kind == "vdw_radius") return each([](uint32_t, const Atom& a) { return fmt("%.2f Å", element(a.element).vdw); });
  if (kind == "x" || kind == "y" || kind == "z") {
    const int c = kind[0] - 'x';
    return each([c](uint32_t, const Atom& a) { return fmt("%.3f", a.pos[size_t(c)]); });
  }
  if (kind == "xyz") return each([](uint32_t, const Atom& a) { char b[96]; std::snprintf(b, sizeof b, "(%.2f, %.2f, %.2f)", a.pos[0], a.pos[1], a.pos[2]); return std::string(b); });
  if (kind == "fractional") {
    if (!s.cell.valid()) return out;
    return each([&](uint32_t, const Atom& a) { const Vec3 f = s.cell.to_fractional(a.pos); char b[96]; std::snprintf(b, sizeof b, "(%.3f, %.3f, %.3f)", f[0], f[1], f[2]); return std::string(b); });
  }
  if (kind == "image") return each([](uint32_t, const Atom& a) { char b[64]; std::snprintf(b, sizeof b, "%d %d %d", a.image[0], a.image[1], a.image[2]); return std::string(b); });
  if (kind == "chain_end") {
    // colour is never the only cue (design/boards/ColourVision): each chain's number at both ends of its backbone
    const auto m = molecules(s);
    for (const auto& b : backbones(s, 3)) {
      if (b.empty()) continue;
      out[b.front()] = std::to_string(m[b.front()]);
      out[b.back()] = std::to_string(m[b.back()]);
    }
    return out;
  }
  if (kind == "molecule" || kind == "molecule_formula") {
    const auto m = molecules(s);
    if (kind == "molecule") return each([&](uint32_t i, const Atom&) { return std::to_string(m[i]); });
    // Hill formula of each molecule, on its first atom
    std::map<int64_t, std::map<std::string, int>> count;
    std::map<int64_t, uint32_t> first;
    for (uint32_t i = 0; i < n; ++i) { ++count[m[i]][element(s.atoms[i].element).symbol]; first.emplace(m[i], i); }
    for (const auto& [mol, c] : count) {
      std::string f;
      auto put = [&](const std::string& e) { const auto it = c.find(e); if (it != c.end()) f += e + (it->second > 1 ? std::to_string(it->second) : ""); };
      const bool carbon = c.count("C") > 0;
      if (carbon) { put("C"); put("H"); }
      for (const auto& [e, k] : c) if (!(carbon && (e == "C" || e == "H"))) put(e);
      out[first[mol]] = sub_digits(f);
    }
    return out;
  }
  if (kind == "stereo") {
    const auto rs = stereo_labels(s);
    const auto ez = ez_labels(s);
    for (size_t i = 0; i < n; ++i) out[i] = rs[i].empty() ? ez[i] : rs[i];
    return out;
  }
  // the perceived chemistry
  const Perception p = perceive(s);
  auto order_sum = [&](uint32_t i) {
    int v = p.implicit_h[i];
    for (size_t k = 0; k < p.nb[i].size(); ++k) v += p.order[i][k];
    return v;
  };
  if (kind == "formal_charge") return each([&](uint32_t i, const Atom&) { return p.charge[i] ? signed_int(p.charge[i]) : std::string("0"); });
  if (kind == "coordination") return each([&](uint32_t i, const Atom&) { return std::to_string(p.nb[i].size() + size_t(p.implicit_h[i])); });
  if (kind == "hydrogens") return each([&](uint32_t i, const Atom&) { return std::to_string(p.hcount[i]); });
  if (kind == "valence") return each([&](uint32_t i, const Atom&) { return std::to_string(order_sum(i)); });
  if (kind == "ring") return each([&](uint32_t i, const Atom&) { return p.smallest_ring[i] ? std::to_string(p.smallest_ring[i]) + "-ring" : std::string(); });
  if (kind == "aromatic") return each([&](uint32_t i, const Atom&) { return p.aromatic[i] ? std::string("ar") : std::string(); });
  if (kind == "hybridisation")
    return each([&](uint32_t i, const Atom& a) -> std::string {
      const int z = a.element;
      if (z <= 2 || z == 10 || z == 18 || z == 36 || z == 54 || z == 86 || metal(z) || p.nb[i].size() + size_t(p.implicit_h[i]) == 0) return "";
      int doubles = 0, triples = 0;
      bool arom = false;
      for (size_t k = 0; k < p.nb[i].size(); ++k) {
        arom = arom || p.arom_bond[i][k];
        if (p.order[i][k] == 2) ++doubles;
        if (p.order[i][k] == 3) ++triples;
      }
      // by the neighbours as well (steric number): S in a sulfone or sulfate, with four, stays sp³
      const size_t nn = p.nb[i].size() + size_t(p.implicit_h[i]);
      if ((triples > 0 || doubles >= 2) && nn <= 2) return "sp";
      if ((doubles >= 1 || arom) && nn <= 3) return "sp²";
      if (z == 5 && nn == 3) return "sp²";
      return nn <= 4 ? "sp³" : "";
    });
  if (kind == "oxidation_state")
    return each([&](uint32_t i, const Atom& a) -> std::string {
      const double en = pauling_electronegativity(a.element);
      if (en <= 0) return "";
      int ox = p.charge[i];
      for (size_t k = 0; k < p.nb[i].size(); ++k) {
        const int zb = s.atoms[p.nb[i][k]].element;
        if (zb == a.element) continue;
        const double eb = pauling_electronegativity(zb);
        if (eb <= 0) return "";
        if (eb > en) ox += p.order[i][k];
        else if (eb < en) ox -= p.order[i][k];
      }
      if (p.implicit_h[i] > 0 && a.element != 1) ox += (pauling_electronegativity(1) > en ? 1 : -1) * p.implicit_h[i];
      return signed_int(ox);
    });
  throw std::invalid_argument("unknown atom label '" + kind + "'");
}

std::vector<BondLabel> bond_labels(const System& s, const std::string& kind, const std::vector<std::string>* types, const ForceField* ff) {
  std::vector<BondLabel> out;
  out.reserve(s.bonds.size());
  const bool cell = s.cell.valid();
  auto sep = [&](const Bond& b) { const Vec3 d = s.atoms[b.j].pos - s.atoms[b.i].pos; return cell ? s.cell.minimum_image(d) : d; };
  for (const auto& b : s.bonds) {
    BondLabel l;
    l.i = b.i, l.j = b.j;
    const Vec3 raw = s.atoms[b.j].pos - s.atoms[b.i].pos;
    l.crossing = cell && norm(raw - sep(b)) > 1e-3;
    out.push_back(l);
  }
  const bool needs_p = kind == "ring" || kind == "rotatable" || kind == "chemical" || kind == "order" || kind == "order_name";
  Perception p;
  if (needs_p) p = perceive(s);
  auto order_of = [&](const Bond& b) -> int {   // the file's order, else the perceived one (aromatic 4)
    if (b.order > 0) return b.order;
    for (size_t k = 0; k < p.nb[b.i].size(); ++k)
      if (p.nb[b.i][k] == b.j) return p.arom_bond[b.i][k] ? 4 : p.order[b.i][k];
    return 1;
  };
  auto in_ring = [&](const Bond& b) {
    for (size_t k = 0; k < p.nb[b.i].size(); ++k) if (p.nb[b.i][k] == b.j) return bool(p.ring_bond[b.i][k]);
    return false;
  };
  // the force field's term of each bond (harmonic, class II or another form), by its atom pair
  struct Term { std::string style; double r0 = 0, k2 = 0, k3 = 0, k4 = 0; };
  std::map<std::pair<uint32_t, uint32_t>, Term> term;
  auto key = [](uint32_t a, uint32_t c) { return std::make_pair(std::min(a, c), std::max(a, c)); };
  if (ff) {
    for (const auto& t : ff->bonds) term[key(t.i, t.j)] = {"harmonic", t.r0, t.k, 0, 0};
    for (const auto& t : ff->bonds2) term[key(t.i, t.j)] = {"class2", t.r0, t.k2, t.k3, t.k4};
    for (const auto& t : ff->bonds_x) term[key(t.i, t.j)] = {t.form == 3 ? "fene" : "form " + std::to_string(t.form)};
  }
  for (size_t k = 0; k < s.bonds.size(); ++k) {
    const Bond& b = s.bonds[k];
    const Atom &A = s.atoms[b.i], &B = s.atoms[b.j];
    const std::string ea = element(A.element).symbol, eb = element(B.element).symbol;
    const double r = norm(sep(b));
    std::string& t = out[k].text;
    if (kind == "index") t = std::to_string(k + 1);
    else if (kind == "length") t = fmt("%.3f Å", r);
    else if (kind == "midpoint") { const Vec3 m = A.pos + sep(b) * 0.5; char c[96]; std::snprintf(c, sizeof c, "(%.2f, %.2f, %.2f)", m[0], m[1], m[2]); t = c; }
    else if (kind == "name") t = (A.name.empty() ? ea + std::to_string(b.i + 1) : A.name) + "–" + (B.name.empty() ? eb + std::to_string(b.j + 1) : B.name);
    else if (kind == "type") t = type_label(s, b.i, types) + "–" + type_label(s, b.j, types);
    else if (kind == "delta_en") {
      const double x = pauling_electronegativity(A.element), y = pauling_electronegativity(B.element);
      if (x > 0 && y > 0) t = fmt("Δχ %.2f", std::fabs(x - y));
    } else if (kind == "order" || kind == "order_name" || kind == "chemical") {
      const int o = order_of(b);
      if (kind == "order") t = o == 4 ? "1.5" : o == 5 ? "am" : o == kBondDative ? "dative" : std::to_string(o);
      else if (kind == "order_name") t = o == 1 ? "single" : o == 2 ? "double" : o == 3 ? "triple" : o == 4 ? "aromatic" : o == 5 ? "amide" : o == kBondDative ? "dative" : "?";
      else {
        const char* sym = o == 2 ? "=" : o == 3 ? "≡" : o == 4 ? ":" : o == kBondDative ? "→" : "–";
        // as bonds are usually written: H last, else the less electronegative atom first (C=O, O–H, S=O, Si–O)
        auto first = [&](int za, int zb) {
          if (za == 1 || zb == 1) return zb == 1;
          const double xa = pauling_electronegativity(za), xb = pauling_electronegativity(zb);
          return xa != xb ? xa < xb : std::string(element(za).symbol) <= element(zb).symbol;
        };
        t = first(A.element, B.element) ? ea + sym + eb : eb + sym + ea;
      }
    } else if (kind == "ring") t = in_ring(b) ? "ring" : "chain";
    else if (kind == "rotatable") {
      auto heavy = [&](uint32_t a) { int h = 0; for (uint32_t x : p.nb[a]) h += s.atoms[x].element != 1; return h; };
      const bool rot = order_of(b) == 1 && !in_ring(b) && A.element != 1 && B.element != 1 && heavy(b.i) > 1 && heavy(b.j) > 1;
      t = rot ? "rot" : "";
    } else if (kind == "ff_style" || kind == "ff_r0" || kind == "ff_k" || kind == "energy") {
      const auto it = term.find(key(b.i, b.j));
      if (it == term.end()) continue;
      const Term& f = it->second;
      const double d = r - f.r0;
      if (kind == "ff_style") t = f.style;
      else if (f.style == "harmonic" || f.style == "class2") {
        if (kind == "ff_r0") t = fmt("%.3f Å", f.r0);
        else if (kind == "ff_k") t = fmt("%.1f", f.k2) + (f.style == "class2" ? " (K2)" : "");
        else t = fmt("%.3f kcal/mol", f.k2 * d * d + f.k3 * d * d * d + f.k4 * d * d * d * d);   // k (r − r0)², class II + K3, K4
      }
    } else throw std::invalid_argument("unknown bond label '" + kind + "'");
  }
  return out;
}

}  // namespace caps
