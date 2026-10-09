// 2D sheets, surface sites, terminations, VASP POSCAR conventions and the 2D validator (caps/slab2d.hpp).
#include "caps/slab2d.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <numeric>
#include <regex>
#include <set>
#include <sstream>
#include <stdexcept>

#include "caps/elements.hpp"
#include "npy_random.hpp"

namespace caps {

namespace {

std::string sym(int z) { return element(z).symbol; }
int el_of(const std::string& s) {
  const int z = element_from_symbol(s);
  if (z <= 0) throw std::invalid_argument("unknown element " + s);
  return z;
}
double py_mod1(double x) { return x - std::floor(x); }   // Python's x % 1.0 for floats

std::string fmt(const char* f, double v) { char b[64]; std::snprintf(b, sizeof b, f, v); return b; }

Json read_json_file(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("cannot read " + path);
  std::stringstream ss;
  ss << in.rdbuf();
  return Json::parse(ss.str());
}

// ASE Atoms.repeat((na, nb, 1)): images m0 (outer) × m1, all atoms per image, the cell's a and b multiplied
System repeat_ab(const System& s, int na, int nb) {
  System out = s;
  out.atoms.clear();
  out.bonds.clear();
  for (int m0 = 0; m0 < na; ++m0)
    for (int m1 = 0; m1 < nb; ++m1)
      for (const Atom& a0 : s.atoms) {
        Atom a = a0;
        const Vec3 off{m0 * s.cell.a[0] + m1 * s.cell.b[0], m0 * s.cell.a[1] + m1 * s.cell.b[1], m0 * s.cell.a[2] + m1 * s.cell.b[2]};
        a.pos = a.pos + off;
        a.id = int64_t(out.atoms.size() + 1);
        out.atoms.push_back(a);
      }
  out.cell.a = s.cell.a * double(na);
  out.cell.b = s.cell.b * double(nb);
  return out;
}

// in-plane minimum image of d (a, b periodic; z kept)
Vec3 mic_ab(const Cell& c, Vec3 d) {
  Vec3 best = d;
  double bl = 1e300;
  const Vec3 f = c.to_fractional(d + c.origin);
  const double r0 = std::round(f[0]), r1 = std::round(f[1]);
  for (int i = -1; i <= 1; ++i)
    for (int j = -1; j <= 1; ++j) {
      const Vec3 e = d - c.a * (r0 + i) - c.b * (r1 + j);
      const double l = e[0] * e[0] + e[1] * e[1] + e[2] * e[2];
      if (l < bl) bl = l, best = e;
    }
  return best;
}

struct Cluster { std::vector<uint32_t> idx; std::vector<double> z; };
// atoms split along z by gaps > gap (periodic): each cluster's atoms in z order, z unwrapped (caps_mxene periodic_z_clusters)
std::vector<Cluster> z_clusters(const System& s, double gap) {
  const size_t n = s.atoms.size();
  const double cz = s.cell.c[2];
  std::vector<double> z(n);
  for (size_t i = 0; i < n; ++i) z[i] = cz > 0 ? s.atoms[i].pos[2] - std::floor(s.atoms[i].pos[2] / cz) * cz : s.atoms[i].pos[2];
  std::vector<uint32_t> order(n);
  std::iota(order.begin(), order.end(), 0u);
  std::stable_sort(order.begin(), order.end(), [&](uint32_t a, uint32_t b) { return z[a] < z[b]; });
  if (n == 0) return {};
  std::vector<size_t> big;
  for (size_t k = 0; k < n; ++k) {
    const double next = k + 1 < n ? z[order[k + 1]] : z[order[0]] + cz;
    if (next - z[order[k]] > gap) big.push_back(k);
  }
  if (big.empty()) {
    Cluster c;
    for (uint32_t i : order) c.idx.push_back(i), c.z.push_back(z[i]);
    return {c};
  }
  const size_t start = (big.back() + 1) % n;
  std::rotate(order.begin(), order.begin() + long(start), order.end());
  std::vector<double> zu(n);
  for (size_t k = 0; k < n; ++k) zu[k] = z[order[k]];
  for (size_t k = 1; k < n; ++k) if (zu[k] < zu[0]) zu[k] += cz;
  std::vector<Cluster> out(1);
  out[0].idx.push_back(order[0]), out[0].z.push_back(zu[0]);
  for (size_t k = 1; k < n; ++k) {
    if (zu[k] - zu[k - 1] > gap) out.emplace_back();
    out.back().idx.push_back(order[k]);
    out.back().z.push_back(zu[k]);
  }
  return out;
}

std::map<std::string, int> composition(const System& s, const std::vector<uint32_t>* idx = nullptr) {
  std::map<std::string, int> c;
  if (idx) for (uint32_t i : *idx) ++c[sym(s.atoms[i].element)];
  else for (const auto& a : s.atoms) ++c[sym(a.element)];
  return c;
}

std::map<std::string, int> parse_formula(const std::string& f) {
  std::map<std::string, int> out;
  static const std::regex re("([A-Z][a-z]?)(\\d*)");
  for (auto it = std::sregex_iterator(f.begin(), f.end(), re); it != std::sregex_iterator(); ++it)
    out[(*it)[1].str()] += (*it)[2].str().empty() ? 1 : std::stoi((*it)[2].str());
  return out;
}

// is `c` (restricted to the elements of `core`) a positive multiple of core?
bool multiple_of(const std::map<std::string, int>& c, const std::map<std::string, int>& core) {
  if (core.empty()) return false;
  const auto& [e0, k0] = *core.begin();
  const auto it = c.find(e0);
  if (it == c.end() || it->second == 0 || it->second % k0) return false;
  const int m = it->second / k0;
  for (const auto& [e, k] : core) {
    const auto jt = c.find(e);
    if (jt == c.end() || jt->second != m * k) return false;
  }
  return true;
}

const std::set<std::string> kTermAtoms = {"O", "F", "H", "Cl", "Br", "I", "S", "Se", "Te"};

std::string json_dict(const std::map<std::string, int>& m) {
  std::string s = "{";
  bool first = true;
  for (const auto& [k, v] : m) { s += (first ? "'" : ", '") + k + "': " + std::to_string(v); first = false; }
  return s + "}";
}

// three-fold hollows of one outer layer (caps_mxene hollow_candidates): sorted by rounded fractional coordinates
std::vector<Vec3> hollows(const System& sup, const std::vector<uint32_t>& idxs, double* dnn_out, std::vector<Vec3>* bridges = nullptr) {
  System lay;
  lay.cell = sup.cell;
  lay.cell.c = {0, 0, 50.0};
  for (uint32_t i : idxs) { Atom a; a.element = 1; a.pos = {sup.atoms[i].pos[0], sup.atoms[i].pos[1], 0.0}; lay.atoms.push_back(a); }
  double cut = 4.5;
  auto nb0 = neighbours_pbc(lay, cut, {true, true, false});
  if (nb0.empty()) nb0 = neighbours_pbc(lay, cut = 2.0 * norm(sup.cell.a), {true, true, false});
  if (nb0.empty()) throw std::runtime_error("the outer layer has no in-plane neighbours");
  double dnn = 1e300;
  for (const auto& n : nb0) dnn = std::min(dnn, n.r);
  if (dnn_out) *dnn_out = dnn;
  const auto nbl = neighbours_pbc(lay, dnn * 1.15, {true, true, false});
  std::vector<std::vector<Vec3>> nb(lay.atoms.size());
  for (const auto& n : nbl) nb[n.i].push_back(n.d);
  auto key_of = [&](const Vec3& c) {
    const Vec3 f = lay.cell.to_fractional(c);
    std::array<double, 2> k{};
    for (int q = 0; q < 2; ++q) {
      double v = py_mod1(f[size_t(q)]);
      v = std::nearbyint(v * 1000.0) / 1000.0;
      v = py_mod1(py_mod1(v)) + 0.0;
      k[size_t(q)] = v;
    }
    return k;
  };
  std::map<std::array<double, 2>, Vec3> found, mids;
  for (size_t a = 0; a < lay.atoms.size(); ++a) {
    const auto& D = nb[a];
    for (size_t p = 0; p < D.size(); ++p) {
      if (bridges) { const Vec3 m = lay.atoms[a].pos + D[p] * 0.5; mids[key_of(m)] = m; }
      for (size_t q = p + 1; q < D.size(); ++q)
        if (std::fabs(norm(D[p] - D[q]) - dnn) < 0.15 * dnn) {
          const Vec3 c = lay.atoms[a].pos + (D[p] + D[q]) * (1.0 / 3.0);
          found[key_of(c)] = c;
        }
    }
  }
  std::vector<Vec3> out;
  for (const auto& [k, c] : found) out.push_back(c);
  if (bridges) { bridges->clear(); for (const auto& [k, m] : mids) bridges->push_back(m); }
  return out;
}

// distinct layer planes (z) of a set of atoms, outermost first for the face
std::vector<double> planes(const System& s, const std::string& face, double tol = 0.3) {
  std::vector<double> z;
  for (const auto& a : s.atoms) z.push_back(a.pos[2]);
  std::sort(z.begin(), z.end());
  std::vector<double> p;
  for (double v : z) if (p.empty() || v - p.back() > tol) p.push_back(v);
  if (face == "top") std::reverse(p.begin(), p.end());
  return p;
}

// the site class of an in-plane position on a face (caps_mxene classify_site, general): the atom directly beneath
// (within tol in-plane): in the outer layer → top, the second layer → hcp, deeper → fcc; none → other
std::pair<std::string, std::string> classify(const System& bare, const Vec3& xy, double zref, const std::string& face, double tol = 0.35) {
  int pick = -1;
  for (size_t k = 0; k < bare.atoms.size(); ++k) {
    const Vec3 d = mic_ab(bare.cell, Vec3{bare.atoms[k].pos[0] - xy[0], bare.atoms[k].pos[1] - xy[1], 0.0});
    if (std::hypot(d[0], d[1]) >= tol) continue;
    const double z = bare.atoms[k].pos[2];
    if (face == "top" ? z > zref + 0.1 : z < zref - 0.1) continue;
    if (pick < 0 || (face == "top" ? z > bare.atoms[size_t(pick)].pos[2] : z < bare.atoms[size_t(pick)].pos[2])) pick = int(k);
  }
  if (pick < 0) return {"other", ""};
  const double zp = bare.atoms[size_t(pick)].pos[2];
  const std::string el = sym(bare.atoms[size_t(pick)].element);
  if (std::fabs(zp - zref) < 0.2) return {"top", el};
  const auto pl = planes(bare, face);
  int layer = 0;
  for (size_t q = 0; q < pl.size(); ++q) if (std::fabs(pl[q] - zp) < 0.3) { layer = int(q) + 1; break; }
  return {layer == 2 ? "hcp" : "fcc", el};
}

std::string outer_element(const System& s, const std::string& face) {
  double zx = face == "top" ? -1e300 : 1e300;
  for (const auto& a : s.atoms) if (a.element != 1) zx = face == "top" ? std::max(zx, a.pos[2]) : std::min(zx, a.pos[2]);
  std::map<std::string, int> c;
  for (const auto& a : s.atoms) if (a.element != 1 && std::fabs(a.pos[2] - zx) < 0.2) ++c[sym(a.element)];
  std::string best;
  int bn = 0;
  for (const auto& [e, n] : c) if (n > bn) bn = n, best = e;
  return best;
}

std::vector<std::string> assign_species(size_t n, const std::vector<std::pair<std::string, double>>& spec, npy::Pcg64& rng) {
  double tot = 0;
  for (const auto& [k, v] : spec) tot += v;
  if (!(tot > 0)) throw std::invalid_argument("termination fractions must add up to more than 0");
  std::vector<double> raw;
  std::vector<int> cnt;
  int sum = 0;
  for (const auto& [k, v] : spec) {
    raw.push_back(double(n) * v / tot);
    cnt.push_back(int(std::floor(raw.back())));
    sum += cnt.back();
  }
  const int rem = int(n) - sum;
  std::vector<size_t> ord(spec.size());
  std::iota(ord.begin(), ord.end(), size_t(0));
  std::stable_sort(ord.begin(), ord.end(), [&](size_t a, size_t b) { return raw[a] - cnt[a] > raw[b] - cnt[b]; });
  for (int k = 0; k < rem && size_t(k) < ord.size(); ++k) ++cnt[ord[size_t(k)]];
  std::vector<std::string> labs;
  for (size_t k = 0; k < spec.size(); ++k) for (int c = 0; c < cnt[k]; ++c) labs.push_back(spec[k].first);
  rng.shuffle(labs);
  return labs;
}

}  // namespace

// ---------------------------------------------------------------- neighbours, wrapping
std::vector<Neighbour> neighbours_pbc(const System& s, double cutoff, std::array<bool, 3> pbc) {
  std::vector<Neighbour> out;
  const size_t n = s.atoms.size();
  int nr[3] = {0, 0, 0};
  if (s.cell.valid()) {
    const Vec3 v[3] = {s.cell.a, s.cell.b, s.cell.c};
    const double vol = std::fabs(dot(v[0], cross(v[1], v[2])));
    for (int k = 0; k < 3; ++k)
      if (pbc[size_t(k)]) {
        const double w = vol / norm(cross(v[(k + 1) % 3], v[(k + 2) % 3]));
        nr[k] = int(std::ceil(cutoff / w));
      }
  } else if (pbc[0] || pbc[1]) {
    // in-plane periodic with only a and b given
    const double area = norm(cross(s.cell.a, s.cell.b));
    if (area > 0) {
      nr[0] = pbc[0] ? int(std::ceil(cutoff / (area / norm(s.cell.b)))) : 0;
      nr[1] = pbc[1] ? int(std::ceil(cutoff / (area / norm(s.cell.a)))) : 0;
    }
  }
  const double c2 = cutoff * cutoff;
  for (size_t i = 0; i < n; ++i)
    for (size_t j = 0; j < n; ++j)
      for (int a = -nr[0]; a <= nr[0]; ++a)
        for (int b = -nr[1]; b <= nr[1]; ++b)
          for (int c = -nr[2]; c <= nr[2]; ++c) {
            if (i == j && a == 0 && b == 0 && c == 0) continue;
            const Vec3 d = s.atoms[j].pos + s.cell.a * double(a) + s.cell.b * double(b) + s.cell.c * double(c) - s.atoms[i].pos;
            const double r2 = dot(d, d);
            if (r2 < c2) out.push_back({uint32_t(i), uint32_t(j), d, std::sqrt(r2)});
          }
  return out;
}

void wrap_ase(System& s, std::array<bool, 3> pbc) {
  for (auto& a : s.atoms) {
    Vec3 f = s.cell.to_fractional(a.pos);
    for (int k = 0; k < 3; ++k)
      if (pbc[size_t(k)]) f[size_t(k)] = py_mod1(f[size_t(k)] + 1e-7) - 1e-7;
    a.pos = s.cell.to_cartesian(f);
  }
}

// ---------------------------------------------------------------- sheets
System sheet_from_layers(const SheetSpec& sp) {
  if (sp.layers.empty()) throw std::invalid_argument("a sheet needs at least one layer");
  if (!(sp.a > 0)) throw std::invalid_argument("a sheet needs its lattice constant a");
  System s;
  s.title = sp.name;
  if (sp.lattice == "rectangular") {
    s.cell.a = {sp.a, 0, 0};
    s.cell.b = {0, sp.b > 0 ? sp.b : sp.a, 0};
  } else {
    s.cell.a = {sp.a, 0.0, 0.0};
    s.cell.b = {-sp.a / 2, sp.a * std::sqrt(3.0) / 2, 0.0};
  }
  s.cell.c = {0, 0, sp.c};
  for (const auto& l : sp.layers) {
    double fx = l.fx, fy = l.fy;
    if (l.site == "A") fx = 0.0, fy = 0.0;
    else if (l.site == "B") fx = 2.0 / 3, fy = 1.0 / 3;
    else if (l.site == "C") fx = 1.0 / 3, fy = 2.0 / 3;
    else if (!l.site.empty()) throw std::invalid_argument("stacking site " + l.site + " (use A, B or C, or fx fy)");
    Atom a;
    a.element = el_of(l.element);
    a.name = l.element;
    a.id = int64_t(s.atoms.size() + 1);
    a.pos = {fx * s.cell.a[0] + fy * s.cell.b[0], fx * s.cell.a[1] + fy * s.cell.b[1], l.z};
    s.atoms.push_back(a);
  }
  s.notes.push_back("sheet " + sp.name + (sp.source.empty() ? "" : " (" + sp.source + ")"));
  return s;
}

std::vector<std::string> sheet_preset_names(const std::string& data_dir) {
  std::vector<std::string> out;
  const Json j = read_json_file(data_dir + "/sheets/sheets.json");
  for (const auto& e : j["sheets"].items()) out.push_back(e.text("name"));
  return out;
}

SheetSpec sheet_preset(const std::string& data_dir, const std::string& name, double a, std::vector<std::string>* notes) {
  auto low = [](std::string s) { for (char& c : s) c = char(std::tolower(static_cast<unsigned char>(c))); return s; };
  const Json j = read_json_file(data_dir + "/sheets/sheets.json");
  for (const auto& e : j["sheets"].items()) {
    if (low(e.text("name")) != low(name)) continue;
    SheetSpec s;
    s.name = e.text("name");
    s.lattice = e.text("lattice", "hexagonal");
    s.a = e.num("a", 0);
    s.b = e.num("b", 0);
    s.source = e.text("source");
    for (const auto& l : e["layers"].items()) s.layers.push_back({l.text("element"), l.text("site"), l.num("fx", 0), l.num("fy", 0), l.num("z", 0)});
    if (a > 0 && std::fabs(a - s.a) > 1e-12) {
      const double f = a / s.a;
      for (auto& l : s.layers) l.z *= f;
      s.b *= f;
      if (notes) notes->push_back("a = " + fmt("%.6f", a) + " Å instead of " + fmt("%.6f", s.a) + ": layer heights scaled by the same ratio (a start for the DFT cell relaxation)");
      s.a = a;
    }
    return s;
  }
  std::string list;
  for (const auto& n : sheet_preset_names(data_dir)) list += (list.empty() ? "" : ", ") + n;
  throw std::invalid_argument("no sheet preset " + name + " (presets: " + list + "; or give the layers, or isolate one from a structure)");
}

System isolate_layer(const System& src0, const IsolateOptions& o, IsolateReport* rep) {
  System src = src0;
  if (!o.remove.empty()) {
    std::set<int> rm;
    for (const auto& e : o.remove) rm.insert(el_of(e));
    std::vector<Atom> keep;
    for (const auto& a : src.atoms) if (!rm.count(a.element)) keep.push_back(a);
    if (rep) rep->notes.push_back(std::to_string(src.atoms.size() - keep.size()) + " atoms removed first");
    src.atoms = keep;
    src.bonds.clear();
  }
  if (!src.cell.valid()) throw std::invalid_argument("isolate: the structure needs a periodic cell");
  const auto core = parse_formula(o.formula);
  const auto cl = z_clusters(src, o.gap);
  std::vector<const Cluster*> good;
  for (const auto& c : cl) {
    const auto comp = composition(src, &c.idx);
    if (core.empty() ? true : multiple_of(comp, core) && comp.size() == core.size()) good.push_back(&c);
  }
  if (rep) rep->layers = int(cl.size()), rep->complete = int(good.size());
  if (good.empty()) throw std::runtime_error("no complete layer" + (o.formula.empty() ? std::string() : " of " + o.formula) + " found (" + std::to_string(cl.size()) + " layers along z)");
  if (o.index < 0 || size_t(o.index) >= good.size()) throw std::invalid_argument("layer " + std::to_string(o.index) + " of " + std::to_string(good.size()));
  const Cluster& c = *good[size_t(o.index)];
  const double zmin = *std::min_element(c.z.begin(), c.z.end());
  System out;
  out.cell = src.cell;
  out.cell.c = {0, 0, o.c};
  for (size_t k = 0; k < c.idx.size(); ++k) {
    Atom a = src.atoms[c.idx[k]];
    a.pos[2] = c.z[k] - zmin;
    a.id = int64_t(out.atoms.size() + 1);
    out.atoms.push_back(a);
  }
  if (rep) rep->formula = formula_ordered(out);
  out.notes.push_back("one layer (" + formula_ordered(out) + ") of " + std::to_string(cl.size()) + " along z");
  return out;
}

// ---------------------------------------------------------------- sites and terminations
std::vector<SurfaceSite> surface_sites(const System& slab, const std::string& element) {
  std::vector<SurfaceSite> out;
  for (const std::string face : {"top", "bottom"}) {
    const std::string el = element.empty() ? outer_element(slab, face) : element;
    if (el.empty()) continue;
    const int z = el_of(el);
    double zref = face == "top" ? -1e300 : 1e300;
    for (const auto& a : slab.atoms) if (a.element == z) zref = face == "top" ? std::max(zref, a.pos[2]) : std::min(zref, a.pos[2]);
    std::vector<uint32_t> idx;
    for (uint32_t i = 0; i < slab.atoms.size(); ++i) if (slab.atoms[i].element == z && std::fabs(slab.atoms[i].pos[2] - zref) < 0.2) idx.push_back(i);
    double dnn = 0;
    std::vector<Vec3> br;
    const auto hs = hollows(slab, idx, &dnn, &br);
    for (const Vec3& h : hs) {
      const auto [kind, beneath] = classify(slab, h, zref, face);
      out.push_back({{h[0], h[1], zref}, kind, face, beneath, dnn / std::sqrt(3.0)});
    }
    for (uint32_t i : idx) out.push_back({slab.atoms[i].pos, "top", face, el, 0.0});
    for (const Vec3& m : br) out.push_back({{m[0], m[1], zref}, "bridge", face, "", dnn / 2});
  }
  return out;
}

std::vector<TerminationKind> termination_library(const std::string& data_dir, const std::string& surf) {
  std::vector<TerminationKind> out;
  const int zs = element_from_symbol(surf);
  const Json j = read_json_file(data_dir + "/sheets/terminations.json");
  for (const auto& e : j["terminations"].items()) {
    TerminationKind k;
    k.name = e.text("name");
    k.atom = e.text("atom");
    k.tail = e.text("tail");
    const int za = el_of(k.atom);
    bool known = false;
    if (e.has("bond") && e["bond"].has(surf)) { k.bond = e["bond"][surf].number(); known = true; }
    else k.bond = (zs > 0 ? element(zs).covalent : 1.5) + element(za).covalent;
    if (e.has("range") && e["range"].has(surf)) k.lo = e["range"][surf][0].number(), k.hi = e["range"][surf][1].number();
    else k.lo = k.bond - 0.2, k.hi = k.bond + 0.3;
    if (!k.tail.empty()) k.tail_bond = e.has("tail_bond") ? e["tail_bond"].number() : element(za).covalent + element(el_of(k.tail)).covalent;
    k.source = known ? e.text("source") : "sum of covalent radii " + surf + "–" + k.atom + " (Cordero 2008): a starting guess";
    out.push_back(k);
  }
  return out;
}

std::vector<std::pair<std::string, double>> parse_fractions(const std::string& text) {
  std::vector<std::pair<std::string, double>> out;
  std::stringstream ss(text);
  for (std::string part; std::getline(ss, part, ',');) {
    part.erase(0, part.find_first_not_of(" \t"));
    part.erase(part.find_last_not_of(" \t") + 1);
    if (part.empty()) continue;
    const auto c = part.find(':');
    out.push_back({part.substr(0, c), c == std::string::npos ? 1.0 : std::stod(part.substr(c + 1))});
  }
  return out;
}

System terminate_slab(const System& slab, const TerminateOptions& o, TerminateReport* rep) {
  if (!slab.cell.valid()) throw std::invalid_argument("terminate: the slab needs a cell");
  if (o.na < 1 || o.nb < 1) throw std::invalid_argument("terminate: supercell must be ≥ 1 × 1");
  System sup = repeat_ab(slab, o.na, o.nb);
  const std::string surf = o.surface_element.empty() ? outer_element(sup, "top") : o.surface_element;
  const int zsurf = el_of(surf);
  const auto lib = termination_library(o.data_dir, surf);
  auto kind = [&](const std::string& n) {
    for (const auto& k : lib) if (k.name == n) return k;
    throw std::invalid_argument("unknown termination " + n + " (data/sheets/terminations.json)");
  };
  double zt_max = -1e300, zt_min = 1e300;
  for (const auto& a : sup.atoms) if (a.element == zsurf) zt_max = std::max(zt_max, a.pos[2]), zt_min = std::min(zt_min, a.pos[2]);
  npy::Pcg64 rng(o.seed);
  struct New { std::string el; Vec3 p; };
  std::vector<New> add;
  for (const std::string face : {"top", "bottom"}) {
    const FaceSpec& fs = face == "top" || !o.janus ? o.top : o.bottom;
    if (fs.fractions.empty()) continue;
    const double zref = face == "top" ? zt_max : zt_min;
    const double sgn = face == "top" ? 1.0 : -1.0;
    std::vector<uint32_t> idx;
    for (uint32_t i = 0; i < sup.atoms.size(); ++i) if (sup.atoms[i].element == zsurf && std::fabs(sup.atoms[i].pos[2] - zref) < 0.2) idx.push_back(i);
    double dnn = 0;
    std::vector<Vec3> br;
    const auto cands = hollows(sup, idx, &dnn, fs.site == "bridge" ? &br : nullptr);
    std::vector<Vec3> sites;
    double r = 0;
    if (fs.site == "top") { for (uint32_t i : idx) sites.push_back(sup.atoms[i].pos); }
    else if (fs.site == "bridge") { sites = br; r = dnn / 2; }
    else {
      for (const Vec3& c : cands) if (classify(sup, c, zref, face).first == fs.site) sites.push_back(c);
      r = dnn / std::sqrt(3.0);
    }
    if (sites.empty()) throw std::runtime_error("no '" + fs.site + "' sites on the " + face + " face");
    if (rep) rep->r = r;
    const auto labels = assign_species(sites.size(), fs.fractions, rng);
    for (size_t s = 0; s < sites.size(); ++s) {
      const TerminationKind k = kind(labels[s]);
      const double d = o.bond.count(k.name) ? o.bond.at(k.name) : k.bond;
      if (fs.site != "top" && d <= r) throw std::invalid_argument(surf + "–" + k.name + " bond " + fmt("%.3f", d) + " Å is not longer than the site's in-plane distance " + fmt("%.3f", r) + " Å");
      const double h = fs.site == "top" ? d : std::sqrt(d * d - r * r);
      const double z = zref + sgn * h;
      add.push_back({k.atom, {sites[s][0], sites[s][1], z}});
      if (!k.tail.empty()) {
        const double t = o.bond.count(k.name + ".tail") ? o.bond.at(k.name + ".tail") : k.tail_bond;
        add.push_back({k.tail, {sites[s][0], sites[s][1], z + sgn * t}});
      }
      if (rep) ++(face == "top" ? rep->top : rep->bottom)[k.name];
    }
  }
  System out = sup;
  for (const auto& n : add) {
    Atom a;
    a.element = el_of(n.el);
    a.name = n.el;
    a.pos = n.p;
    out.atoms.push_back(a);
  }
  out = group_by_species(out, o.order);
  double zmin = 1e300, zmax = -1e300;
  for (const auto& a : out.atoms) zmin = std::min(zmin, a.pos[2]), zmax = std::max(zmax, a.pos[2]);
  const double shift = o.vacuum / 2.0 - zmin;
  for (auto& a : out.atoms) a.pos[2] += shift;
  out.cell.c = {0, 0, (zmax - zmin) + o.vacuum};
  wrap_ase(out, {true, true, false});
  for (size_t i = 0; i < out.atoms.size(); ++i) out.atoms[i].id = int64_t(i + 1);
  out.title = formula_ordered(out, o.order);
  if (rep) rep->notes.push_back(std::to_string(o.na) + " × " + std::to_string(o.nb) + ", " + std::to_string(add.size()) + " termination atoms, vacuum " + fmt("%.1f", o.vacuum) + " Å");
  return out;
}

// ---------------------------------------------------------------- species order, formula, POSCAR
System group_by_species(const System& s, const std::vector<std::string>& order) {
  std::vector<int> species;
  for (const auto& e : order) {
    const int z = element_from_symbol(e);
    if (z > 0 && std::any_of(s.atoms.begin(), s.atoms.end(), [&](const Atom& a) { return a.element == z; })) species.push_back(z);
  }
  for (const auto& a : s.atoms) if (std::find(species.begin(), species.end(), a.element) == species.end()) species.push_back(a.element);
  System out = s;
  out.atoms.clear();
  out.bonds.clear();
  for (int z : species) for (const auto& a : s.atoms) if (a.element == z) out.atoms.push_back(a);
  return out;
}

std::string formula_ordered(const System& s, const std::vector<std::string>& order) {
  const auto c = composition(s);
  std::string f;
  auto put = [&](const std::string& e) { const auto it = c.find(e); if (it != c.end()) f += e + (it->second > 1 ? std::to_string(it->second) : ""); };
  for (const auto& e : order) put(e);
  for (const auto& [e, n] : c) if (std::find(order.begin(), order.end(), e) == order.end()) put(e);
  return f;
}

void write_vasp_poscar(const System& s, const std::string& path, const std::string& comment, const std::vector<bool>& fixed) {
  std::vector<std::pair<std::string, int>> runs;
  for (const auto& a : s.atoms) {
    const std::string e = sym(a.element);
    if (!runs.empty() && runs.back().first == e) ++runs.back().second;
    else runs.push_back({e, 1});
  }
  std::set<std::string> seen;
  for (const auto& r : runs) if (!seen.insert(r.first).second) throw std::runtime_error("POSCAR: species " + r.first + " is not contiguous (group the atoms by element first)");
  std::ofstream f(path);
  if (!f) throw std::runtime_error("cannot write " + path);
  char b[160];
  f << comment << "\n1.0\n";
  for (const Vec3& v : {s.cell.a, s.cell.b, s.cell.c}) { std::snprintf(b, sizeof b, "  %.10f  %.10f  %.10f\n", v[0], v[1], v[2]); f << b; }
  f << " ";
  for (const auto& r : runs) f << " " << (r == runs.front() ? "" : " ") << r.first;
  f << "\n ";
  for (const auto& r : runs) f << " " << (r == runs.front() ? "" : " ") << r.second;
  f << "\n";
  const bool sel = std::any_of(fixed.begin(), fixed.end(), [](bool x) { return x; });
  if (sel) f << "Selective dynamics\n";
  f << "Cartesian\n";
  for (size_t i = 0; i < s.atoms.size(); ++i) {
    const Vec3& p = s.atoms[i].pos;
    std::snprintf(b, sizeof b, "  %.10f  %.10f  %.10f", p[0], p[1], p[2]);
    f << b;
    if (sel) f << (i < fixed.size() && fixed[i] ? "   F   F   F" : "   T   T   T");
    f << "\n";
  }
}

System read_vasp_poscar(const std::string& path, std::map<std::string, std::string>* bad, std::string* comment) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("cannot open " + path);
  std::vector<std::string> L;
  for (std::string l; std::getline(in, l);) { if (!l.empty() && l.back() == '\r') l.pop_back(); L.push_back(l); }
  auto split = [](const std::string& s) { std::stringstream ss(s); std::vector<std::string> t; for (std::string w; ss >> w;) t.push_back(w); return t; };
  if (L.size() < 8) throw std::runtime_error(path + ": too short for a POSCAR");
  if (comment) *comment = L[0];
  System s;
  s.title = L[0];
  s.source_format = "poscar";
  const double scale = std::stod(split(L[1]).at(0));
  Vec3 v[3];
  for (int r = 0; r < 3; ++r) { const auto t = split(L[size_t(2 + r)]); v[r] = {std::stod(t.at(0)), std::stod(t.at(1)), std::stod(t.at(2))}; }
  const double f = scale < 0 ? std::cbrt(-scale / std::fabs(dot(v[0], cross(v[1], v[2])))) : scale;
  s.cell.a = v[0] * f, s.cell.b = v[1] * f, s.cell.c = v[2] * f;
  const auto tok = split(L[5]);
  if (tok.empty() || std::isdigit(static_cast<unsigned char>(tok[0][0])))
    throw std::runtime_error(path + ": a VASP 4 POSCAR without species names is not supported");
  std::vector<int> counts;
  for (const auto& w : split(L[6])) counts.push_back(std::stoi(w));
  if (counts.size() != tok.size()) throw std::runtime_error(path + ": " + std::to_string(tok.size()) + " species, " + std::to_string(counts.size()) + " counts");
  size_t k = 7;
  bool selective = false;
  if (!L.at(k).empty() && (std::tolower(static_cast<unsigned char>(L[k].find_first_not_of(' ') != std::string::npos ? L[k][L[k].find_first_not_of(' ')] : ' ')) == 's')) selective = true, ++k;
  const std::string mode = L.at(k);
  const char m0 = char(std::tolower(static_cast<unsigned char>(mode[mode.find_first_not_of(' ')])));
  const bool cart = m0 == 'c' || m0 == 'k';
  ++k;
  (void)selective;
  for (size_t sp = 0; sp < tok.size(); ++sp) {
    std::string lab = tok[sp], real = lab.substr(0, lab.find_first_of("_/"));
    if (element_from_symbol(real) <= 0 || real == "X") {
      if (lab == "0") real = "O";
      else throw std::runtime_error(path + ": unknown species label '" + lab + "'");
      if (bad) (*bad)[lab] = real;
    }
    for (int i = 0; i < counts[sp]; ++i) {
      const auto p = split(L.at(k++));
      const Vec3 r{std::stod(p.at(0)), std::stod(p.at(1)), std::stod(p.at(2))};
      Atom a;
      a.element = el_of(real);
      a.name = real;
      a.id = int64_t(s.atoms.size() + 1);
      a.pos = cart ? r * f : s.cell.to_cartesian(r);
      s.atoms.push_back(a);
    }
  }
  return s;
}

// ---------------------------------------------------------------- the 2D validator
std::string Validation2D::text(const std::string& name) const {
  std::ostringstream o;
  o << "=== " << (name.empty() ? info.text("name", "structure") : name) << " : " << status << " ===\n";
  auto get = [&](const char* k) { return info.has(k) ? info[k].dump(0) : std::string("-"); };
  o << "composition : " << get("composition") << "\n";
  if (info.has("cell")) o << "cell        : a=" << fmt("%.4f", info["cell"].num("a", 0)) << " b=" << fmt("%.4f", info["cell"].num("b", 0)) << " c=" << fmt("%.2f", info["cell"].num("c", 0))
                          << " A, gamma=" << fmt("%.1f", info["cell"].num("gamma", 0)) << " deg\n";
  o << "blocks      : " << get("blocks") << "   layers per element: " << get("layers_per_element") << "\n";
  o << "vacuum gap  : " << fmt("%.1f", info.num("vacuum_gap", 0)) << " A    min distance: " << fmt("%.2f", info.num("min_distance", 0)) << " A\n";
  o << "terminations: " << get("faces") << "   sites: " << get("sites") << "\n";
  if (info.has("layers")) {
    o << "layers (z in A):\n";
    for (const auto& l : info["layers"].items()) o << "   " << fmt("%8.3f", l["z"].number()) << "   " << l["composition"].dump(0) << "\n";
  }
  if (findings.empty()) o << "findings: none\n";
  else {
    o << "findings:\n";
    for (const auto& f : findings) o << "   [" << f.level << "] " << f.text << "\n";
  }
  return o.str();
}

Json Validation2D::json() const {
  Json j = Json::object();
  j["status"] = status;
  Json fs = Json::array();
  for (const auto& f : findings) {
    Json x = Json::object();
    x["level"] = f.level;
    x["text"] = f.text;
    Json at = Json::array();
    for (int a : f.atoms) at.push_back(a);
    x["atoms"] = at;
    fs.push_back(x);
  }
  j["findings"] = fs;
  j["info"] = info;
  return j;
}

Validation2D validate_2d(const System& s, const ValidateOptions& o) {
  Validation2D R;
  auto add = [&](const std::string& lvl, const std::string& t, std::vector<int> atoms = {}) { R.findings.push_back({lvl, t, std::move(atoms)}); };
  const size_t n = s.atoms.size();
  std::vector<std::string> el(n);
  for (size_t i = 0; i < n; ++i) el[i] = sym(s.atoms[i].element);
  const auto comp = composition(s);
  {
    Json c = Json::object();
    for (const auto& [e, k] : comp) c[e] = k;
    R.info["composition"] = c;
  }
  const double la = norm(s.cell.a), lb = norm(s.cell.b), lc = norm(s.cell.c);
  {
    Json c = Json::object();
    c["a"] = la, c["b"] = lb, c["c"] = lc;
    c["gamma"] = la > 0 && lb > 0 ? std::acos(std::clamp(dot(s.cell.a, s.cell.b) / (la * lb), -1.0, 1.0)) * 180 / M_PI : 0.0;
    R.info["cell"] = c;
  }
  if (std::fabs(s.cell.c[0]) > 1e-4 || std::fabs(s.cell.c[1]) > 1e-4) add("WARN", "c axis is not along z; vacuum/layer checks assume a slab normal to z");
  for (const auto& [lab, real] : o.bad_labels)
    add("ERROR", "species label '" + lab + "' is not an element symbol (digit zero instead of letter O?); VASP/POTCAR mapping will fail. Treated as '" + real + "' for this analysis");
  if (!o.comment.empty()) {
    // leading tokens of the comment that look like a formula (with a digit)
    std::map<std::string, int> ct;
    std::stringstream ss(o.comment);
    static const std::regex tokre("(?:[A-Z][a-z]?\\d*)+");
    for (std::string t; ss >> t;) {
      if (!std::regex_match(t, tokre) || t.find_first_of("0123456789") == std::string::npos) break;
      for (const auto& [e, k] : parse_formula(t)) ct[e] += k;
    }
    if (!ct.empty() && ct != comp) add("WARN", "comment line says " + json_dict(ct) + " but the atoms are " + json_dict(comp) + " (stale label; harmless to VASP but misleading)");
  }
  // the bare layer's formula
  std::map<std::string, int> core = parse_formula(o.core);
  if (core.empty() && !o.expect.empty()) {
    static const std::regex re("([A-Z][a-z]?)(\\d*)");
    int k = 0;
    for (auto it = std::sregex_iterator(o.expect.begin(), o.expect.end(), re); it != std::sregex_iterator() && k < 2; ++it, ++k)
      core[(*it)[1].str()] = (*it)[2].str().empty() ? 1 : std::stoi((*it)[2].str());
  }
  if (!o.expect.empty()) {
    const auto exp = parse_formula(o.expect);
    bool ok = !exp.empty() && comp.size() == exp.size();
    if (ok) {
      const auto& [e0, k0] = *exp.begin();
      const double mult = comp.count(e0) ? double(comp.at(e0)) / k0 : 0;
      for (const auto& [e, k] : exp) ok = ok && comp.count(e) && std::fabs(comp.at(e) - mult * k) < 1e-6;
      ok = ok && mult > 0;
    }
    if (!ok) add("ERROR", "composition " + json_dict(comp) + " does not match expected " + o.expect + " x N");
  }
  // blocks (z clusters), fragments, floaters
  const auto cl = z_clusters(s, 2.0);
  if (core.empty()) {
    // the largest block's elements other than terminations, reduced by their gcd
    const Cluster* big = nullptr;
    for (const auto& c : cl) if (!big || c.idx.size() > big->idx.size()) big = &c;
    if (big) {
      int g = 0;
      for (const auto& [e, k] : composition(s, &big->idx)) if (!kTermAtoms.count(e)) core[e] = k, g = std::gcd(g, k);
      if (g > 1) for (auto& [e, k] : core) k /= g;
    }
  }
  const Cluster* sheet = nullptr;
  int nsheets = 0;
  for (const auto& c : cl) {
    const auto cc = composition(s, &c.idx);
    bool has_core = false;
    for (const auto& [e, k] : core) has_core |= cc.count(e) > 0;
    std::vector<int> at(c.idx.begin(), c.idx.end());
    const double z0 = *std::min_element(c.z.begin(), c.z.end()), z1 = *std::max_element(c.z.begin(), c.z.end());
    if (!has_core)
      add("ERROR", json_dict(cc) + " at z = " + fmt("%.2f", z0) + "-" + fmt("%.2f", z1) + " A are not attached to any sheet (> 2 A vertical gap): floating in vacuum", at);
    else if (multiple_of(cc, core)) { ++nsheets; if (!sheet) sheet = &c; }
    else
      add("ERROR", "incomplete sheet fragment " + json_dict(cc) + " at z = " + fmt("%.2f", z0) + "-" + fmt("%.2f", z1) + " A (a sheet must be a multiple of " + o.core +
                       (o.core.empty() ? std::string("its formula") : "") + " in one contiguous block)", at);
  }
  R.info["blocks"] = nsheets;
  if (nsheets != 1) add("ERROR", "found " + std::to_string(nsheets) + " complete sheets (expected exactly 1)");
  // layers of every block
  {
    Json layers = Json::array();
    for (const auto& c : cl) {
      std::vector<size_t> ord(c.idx.size());
      std::iota(ord.begin(), ord.end(), size_t(0));
      std::stable_sort(ord.begin(), ord.end(), [&](size_t a, size_t b) { return c.z[a] < c.z[b]; });
      std::vector<size_t> cur{ord[0]};
      auto flush = [&] {
        double zm = 0;
        std::map<std::string, int> cc;
        for (size_t k : cur) zm += c.z[k], ++cc[el[c.idx[k]]];
        Json l = Json::object();
        l["z"] = std::round(zm / double(cur.size()) * 1000) / 1000;
        Json cj = Json::object();
        for (const auto& [e, k] : cc) cj[e] = k;
        l["composition"] = cj;
        layers.push_back(l);
      };
      for (size_t q = 1; q < ord.size(); ++q) {
        if (c.z[ord[q]] - c.z[cur.back()] > 0.35) { flush(); cur = {ord[q]}; }
        else cur.push_back(ord[q]);
      }
      flush();
    }
    R.info["layers"] = layers;
  }
  System sub;   // the sheet's atoms, z unwrapped
  std::vector<uint32_t> sub_of;
  if (sheet) {
    Json lpe = Json::object();
    bool bad = false;
    std::string have;
    for (const auto& [e, k] : core) {
      std::vector<double> zz;
      for (size_t q = 0; q < sheet->idx.size(); ++q) if (el[sheet->idx[q]] == e) zz.push_back(sheet->z[q]);
      std::sort(zz.begin(), zz.end());
      int d = zz.empty() ? 0 : 1;
      for (size_t q = 1; q < zz.size(); ++q) d += zz[q] - zz[q - 1] > 0.3;
      lpe[e] = d;
      bad |= d != k;
      have += (have.empty() ? "" : " and ") + std::to_string(d) + " " + e;
    }
    R.info["layers_per_element"] = lpe;
    std::string want;
    for (const auto& [e, k] : core) want += (want.empty() ? "" : " and ") + std::to_string(k) + " " + e;
    if (bad) add("ERROR", "sheet has " + have + " distinct layers (one sheet has " + want + ")");
    sub.cell = s.cell;
    for (size_t q = 0; q < sheet->idx.size(); ++q) {
      const uint32_t i = sheet->idx[q];
      if (!core.count(el[i])) continue;
      Atom a = s.atoms[i];
      a.pos[2] = sheet->z[q];
      sub.atoms.push_back(a);
      sub_of.push_back(i);
    }
  }
  // vacuum
  {
    std::vector<double> zs;
    for (const auto& a : s.atoms) zs.push_back(lc > 0 ? a.pos[2] - std::floor(a.pos[2] / lc) * lc : a.pos[2]);
    std::sort(zs.begin(), zs.end());
    double vac = 0;
    for (size_t q = 0; q < zs.size(); ++q) vac = std::max(vac, (q + 1 < zs.size() ? zs[q + 1] : zs[0] + lc) - zs[q]);
    R.info["vacuum_gap"] = vac;
    if (vac < 10) add("ERROR", "vacuum gap only " + fmt("%.1f", vac) + " A");
    else if (vac < 15) add("WARN", "vacuum gap " + fmt("%.1f", vac) + " A (< 15 A recommended)");
  }
  // distances
  const auto nb32 = neighbours_pbc(s, 3.2);
  std::vector<int> cnt(n, 0);
  double dmin = 1e300;
  std::vector<std::string> close;
  std::vector<int> close_atoms;
  for (const auto& p : nb32) {
    ++cnt[p.i];
    dmin = std::min(dmin, p.r);
    if (p.i < p.j) {
      const std::set<std::string> pr{el[p.i], el[p.j]};
      const double lim = pr == std::set<std::string>{"O", "H"} ? 0.9 : (pr.count("H") ? 1.0 : 1.4);
      if (p.r < lim) {
        if (close.size() < 3) close.push_back("(" + std::to_string(p.i + 1) + ", '" + el[p.i] + "', " + std::to_string(p.j + 1) + ", '" + el[p.j] + "', " + fmt("%.3f", p.r) + ")");
        close_atoms.push_back(int(p.i)), close_atoms.push_back(int(p.j));
      }
    }
  }
  for (size_t i = 0; i < n; ++i) if (!cnt[i]) add("ERROR", "atom #" + std::to_string(i + 1) + " (" + el[i] + ") has no neighbour within 3.2 A (isolated)", {int(i)});
  if (!close_atoms.empty()) {
    std::string e;
    for (const auto& c : close) e += (e.empty() ? "" : ", ") + c;
    add("ERROR", std::to_string(close_atoms.size() / 2) + " atom pairs too close, e.g. [" + e + "]", close_atoms);
  }
  R.info["min_distance"] = nb32.empty() ? 0.0 : dmin;
  // terminations: atoms bonded to the outer element of the sheet
  std::string surf;
  if (!sub.atoms.empty()) surf = outer_element(sub, "top");
  const int zsurf = surf.empty() ? 0 : element_from_symbol(surf);
  std::vector<TerminationKind> lib;
  if (zsurf && !o.data_dir.empty()) try { lib = termination_library(o.data_dir, surf); } catch (...) {}
  auto range_of = [&](const std::string& kind) -> std::pair<double, double> {
    for (const auto& k : lib) if (k.name == kind) return {k.lo, k.hi};
    return {0, 0};
  };
  double hi_max = 2.35;
  for (const auto& k : lib) hi_max = std::max(hi_max, k.hi);
  const auto nbt = neighbours_pbc(s, std::max(2.8, hi_max + 0.45));
  std::vector<std::vector<double>> near(n);
  for (const auto& p : nbt) if (s.atoms[p.j].element == zsurf) near[p.i].push_back(p.r);
  std::map<uint32_t, std::vector<std::pair<uint32_t, Vec3>>> hof;   // O/N → its H
  std::set<uint32_t> hbound;
  for (const auto& p : neighbours_pbc(s, 1.3))
    if ((el[p.i] == "O" || el[p.i] == "N") && el[p.j] == "H" && !core.count(el[p.i])) hof[p.i].push_back({p.j, p.d}), hbound.insert(p.j);
  for (size_t k = 0; k < n; ++k) {
    if (el[k] != "H" || hbound.count(uint32_t(k))) continue;
    double best = 1e300;
    for (size_t o2 = 0; o2 < n; ++o2)
      if (el[o2] == "O") best = std::min(best, norm(s.cell.minimum_image(s.atoms[o2].pos - s.atoms[k].pos)));
    add("ERROR", best < 1e299 ? "H atom #" + std::to_string(k + 1) + " is not bonded to any O (nearest O at " + fmt("%.2f", best) + " A; an O-H bond is about 0.97 A)"
                              : "H atom #" + std::to_string(k + 1) + " has no O", {int(k)});
  }
  double zc = 0;
  int nsurf = 0;
  for (size_t i = 0; i < n; ++i) if (s.atoms[i].element == zsurf) zc += s.atoms[i].pos[2], ++nsurf;
  zc = nsurf ? zc / nsurf : 0;
  std::map<std::string, std::map<std::string, int>> face_cnt{{"top", {}}, {"bottom", {}}}, site_cnt{{"top", {}}, {"bottom", {}}};
  for (size_t k = 0; k < n; ++k) {
    if (el[k] == "H" || core.count(el[k]) || !(kTermAtoms.count(el[k]) || el[k] == "N")) continue;
    const std::string kind = (el[k] == "O" || el[k] == "N") && hof.count(uint32_t(k)) ? el[k] + "H" : el[k];
    auto ds = near[k];
    std::sort(ds.begin(), ds.end());
    if (ds.empty()) { add("ERROR", kind + " atom #" + std::to_string(k + 1) + " has no " + surf + " within " + fmt("%.1f", std::max(2.8, hi_max + 0.45)) + " A: not bonded to the surface", {int(k)}); continue; }
    auto [lo, hi] = range_of(kind);
    if (hi > 0) {
      std::vector<double> bonded;
      for (double x : ds) if (x <= hi + 0.05) bonded.push_back(x);
      const double mx = bonded.empty() ? ds.back() : bonded.back();
      if (ds.front() < lo || mx > hi + 0.05) {
        std::string l;
        for (double x : ds) l += (l.empty() ? "" : ", ") + fmt("%.2f", x);
        add("WARN", kind + " atom #" + std::to_string(k + 1) + ": " + surf + "-X distances [" + l + "] outside " + fmt("%.2f", lo) + "-" + fmt("%.2f", hi) + " A", {int(k)});
      }
      if (bonded.size() != 1 && bonded.size() != 3)
        add("WARN", kind + " atom #" + std::to_string(k + 1) + " has " + std::to_string(bonded.size()) + " " + surf + " neighbours (hollow=3, top=1)", {int(k)});
    }
    const std::string face = s.atoms[k].pos[2] > zc ? "top" : "bottom";
    ++face_cnt[face][kind];
    if (hof.count(uint32_t(k))) {
      const auto& hs = hof.at(uint32_t(k));
      if (hs.size() > 1) add("ERROR", el[k] + " atom #" + std::to_string(k + 1) + " carries " + std::to_string(hs.size()) + " H atoms", {int(k)});
      for (const auto& [h, D] : hs) {
        const double L = norm(D);
        if (el[k] == "O" && !(L >= 0.93 && L <= 1.03))
          add("ERROR", "O-H bond #" + std::to_string(k + 1) + "-" + std::to_string(h + 1) + " is " + fmt("%.2f", L) + " A (expected 0.93-1.03 A)", {int(k), int(h)});
        const double up = (face == "top" ? 1.0 : -1.0) * D[2] / L;
        if (up < 0.9) add("WARN", el[k] + "-H #" + std::to_string(k + 1) + "-" + std::to_string(h + 1) + " tilts away from the surface normal (cos = " + fmt("%.2f", up) + ")", {int(k), int(h)});
      }
    }
  }
  if (!sub.atoms.empty()) {
    double zmid = 0, ztop = -1e300, zbot = 1e300;
    int ns = 0;
    for (const auto& a : sub.atoms) if (a.element == zsurf) zmid += a.pos[2], ++ns, ztop = std::max(ztop, a.pos[2]), zbot = std::min(zbot, a.pos[2]);
    zmid /= std::max(1, ns);
    for (size_t q = 0; q < sheet->idx.size(); ++q) {
      const uint32_t k = sheet->idx[q];
      if (core.count(el[k]) || el[k] == "H" || !(kTermAtoms.count(el[k]) || el[k] == "N")) continue;
      const std::string face = sheet->z[q] > zmid ? "top" : "bottom";
      ++site_cnt[face][classify(sub, s.atoms[k].pos, face == "top" ? ztop : zbot, face).first];
    }
  }
  auto to_json = [](const std::map<std::string, std::map<std::string, int>>& m) {
    Json j = Json::object();
    for (const auto& [f, c] : m) { Json x = Json::object(); for (const auto& [k, v] : c) x[k] = v; j[f] = x; }
    return j;
  };
  R.info["faces"] = to_json(face_cnt);
  R.info["sites"] = to_json(site_cnt);
  const bool any = !face_cnt["top"].empty() || !face_cnt["bottom"].empty();
  R.info["asymmetric"] = any && face_cnt["top"] != face_cnt["bottom"];
  if (any && face_cnt["top"] != face_cnt["bottom"]) add("WARN", "top and bottom faces differ: use the dipole correction (IDIPOL = 3)");
  if (!surf.empty()) R.info["surface_element"] = surf;
  for (const auto& f : R.findings) if (f.level == "ERROR") R.status = "FAIL";
  return R;
}

}  // namespace caps
