// A molecule on a 2D slab for DFT: anchors, orientations, symmetry-equivalent azimuths, placement, image distance and
// an independent audit (caps/adsorb_dft.hpp).
#include "caps/adsorb_dft.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <deque>
#include <fstream>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>

#include "caps/adsorption.hpp"
#include "caps/analysis.hpp"
#include "caps/crystal.hpp"
#include "caps/uff.hpp"
#include "caps/edit.hpp"
#include "caps/elements.hpp"

namespace caps {

namespace {

std::string fmt(const char* f, double v) { char b[64]; std::snprintf(b, sizeof b, f, v); return b; }

Json read_json_file(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("cannot read " + path);
  std::stringstream ss;
  ss << in.rdbuf();
  return Json::parse(ss.str());
}

// in-plane minimum image (a, b periodic), exact over the neighbouring images
Vec3 mic_ab(const Cell& c, Vec3 d) {
  const Vec3 f = c.to_fractional(d + c.origin);
  const double r0 = std::round(f[0]), r1 = std::round(f[1]);
  Vec3 best = d;
  double bl = 1e300;
  for (int i = -1; i <= 1; ++i)
    for (int j = -1; j <= 1; ++j) {
      const Vec3 e = d - c.a * (r0 + i) - c.b * (r1 + j);
      const double l = dot(e, e);
      if (l < bl) bl = l, best = e;
    }
  return best;
}

// closest contact molecule ↔ slab (any pair; non-H to non-H), in-plane periodic
std::pair<double, double> contacts(const System& slab, const std::vector<Vec3>& mol, const std::vector<int>& mol_el, const Cell& cell) {
  double any = 1e300, heavy = 1e300;
  for (size_t k = 0; k < mol.size(); ++k)
    for (const auto& a : slab.atoms) {
      const double d = norm(mic_ab(cell, a.pos - mol[k]));
      any = std::min(any, d);
      if (mol_el[k] != 1 && a.element != 1) heavy = std::min(heavy, d);
    }
  return {any, heavy};
}

// the unwrapped molecule to its own in-plane images (the eight neighbours)
double image_distance(const std::vector<Vec3>& m, const Vec3& a, const Vec3& b) {
  double best = 1e300;
  for (int i = -1; i <= 1; ++i)
    for (int j = -1; j <= 1; ++j) {
      if (i == 0 && j == 0) continue;
      const Vec3 v = a * double(i) + b * double(j);
      for (const Vec3& p : m)
        for (const Vec3& q : m) best = std::min(best, norm(p - (q + v)));
    }
  return best;
}

// rotation taking unit v onto unit t (Rodrigues; antiparallel: 180° about an axis ⟂ v)
std::array<Vec3, 3> rot_to(Vec3 v, Vec3 t) {
  v = v * (1.0 / norm(v));
  t = t * (1.0 / norm(t));
  const double c = dot(v, t);
  auto K = [](const Vec3& ax) { return std::array<Vec3, 3>{Vec3{0, -ax[2], ax[1]}, Vec3{ax[2], 0, -ax[0]}, Vec3{-ax[1], ax[0], 0}}; };
  auto mm = [](const std::array<Vec3, 3>& A, const std::array<Vec3, 3>& B) {
    std::array<Vec3, 3> C{};
    for (int r = 0; r < 3; ++r) for (int q = 0; q < 3; ++q) for (int k = 0; k < 3; ++k) C[r][q] += A[r][k] * B[k][q];
    return C;
  };
  std::array<Vec3, 3> I{Vec3{1, 0, 0}, Vec3{0, 1, 0}, Vec3{0, 0, 1}};
  if (c > 1 - 1e-9) return I;
  if (c < -1 + 1e-9) {
    Vec3 ax = cross(v, {1, 0, 0});
    ax = norm(ax) > 1e-6 ? ax * (1.0 / norm(ax)) : Vec3{0, 1, 0};
    const auto k = K(ax), kk = mm(k, k);
    for (int r = 0; r < 3; ++r) for (int q = 0; q < 3; ++q) I[r][q] += 2 * kk[r][q];
    return I;
  }
  const Vec3 ax = cross(v, t);
  const double s = norm(ax);
  const auto k = K(ax), kk = mm(k, k);
  for (int r = 0; r < 3; ++r) for (int q = 0; q < 3; ++q) I[r][q] += k[r][q] + kk[r][q] * ((1 - c) / (s * s));
  return I;
}
Vec3 rotate_by(const std::array<Vec3, 3>& R, const Vec3& p) { return {dot(R[0], p), dot(R[1], p), dot(R[2], p)}; }

std::vector<Vec3> orient(const System& mol, const std::string& mode, const std::vector<int>& anchor, double az) {
  std::vector<Vec3> p;
  Vec3 c{0, 0, 0};
  for (const auto& a : mol.atoms) c = c + a.pos;
  c = c * (1.0 / double(mol.atoms.size()));
  for (const auto& a : mol.atoms) p.push_back(a.pos - c);
  if (mode == "parallel" || mode == "upright") {
    // inertia frame: long axis (smallest moment) → x, flat normal (largest) → z; upright: long axis → z
    Vec3 com{0, 0, 0};
    double M = 0;
    for (size_t i = 0; i < p.size(); ++i) { const double m = element(mol.atoms[i].element).mass; com = com + p[i] * m; M += m; }
    com = com * (1.0 / M);
    double I[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    for (size_t i = 0; i < p.size(); ++i) {
      const double m = element(mol.atoms[i].element).mass;
      const Vec3 r = p[i] - com;
      const double r2 = dot(r, r);
      for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b) I[a][b] += m * ((a == b ? r2 : 0) - r[size_t(a)] * r[size_t(b)]);
    }
    double w[3], V[3][3];
    symmetric_eigen3(I, w, V);
    const Vec3 e0{V[0][0], V[1][0], V[2][0]}, e1{V[0][1], V[1][1], V[2][1]}, e2{V[0][2], V[1][2], V[2][2]};
    for (auto& q : p) q = mode == "parallel" ? Vec3{dot(q, e0), dot(q, e1), dot(q, e2)} : Vec3{dot(q, e1), dot(q, e2), dot(q, e0)};
    if (mode == "upright" && !anchor.empty()) {
      // the anchor end down
      double za = 0;
      for (int i : anchor) za += p[size_t(i)][2];
      if (za > 0) for (auto& q : p) q = {q[0], -q[1], -q[2]};
    }
  } else {
    Vec3 a{0, 0, 0};
    for (int i : anchor) a = a + p[size_t(i)];
    a = a * (1.0 / double(anchor.size()));
    if (norm(a) > 1e-9) {
      const auto R = rot_to(a, {0, 0, -1.0});
      for (auto& q : p) q = rotate_by(R, q);
    }
  }
  const double t = az * M_PI / 180;
  for (auto& q : p) q = {std::cos(t) * q[0] - std::sin(t) * q[1], std::sin(t) * q[0] + std::cos(t) * q[1], q[2]};
  return p;
}

std::vector<int> species_order_index(const System& s, const std::vector<std::string>& order) {
  std::vector<int> idx(s.atoms.size());
  std::iota(idx.begin(), idx.end(), 0);
  auto rank = [&](int z) {
    const std::string e = element(z).symbol;
    const auto it = std::find(order.begin(), order.end(), e);
    return it == order.end() ? 99 : int(it - order.begin());
  };
  std::stable_sort(idx.begin(), idx.end(), [&](int a, int b) { return rank(s.atoms[size_t(a)].element) < rank(s.atoms[size_t(b)].element); });
  return idx;
}


// slab + the molecule's positions m → the complex as written (vacuum centred, in-plane wrapped, species grouped with
// the slab's atoms first in each block) and its checks: contacts, the image distance of the unwrapped molecule (with the
// smallest passing supercell), the slab part re-validated
void finish_complex(AdsorbComplex& cx, const System& slab, const System& mol, const std::vector<Vec3>& m, const AdsorbSetOptions& o) {
  // slab + molecule, vacuum centred, in-plane wrapped, species grouped (slab atoms first in each block)
  System both = slab;
  std::vector<char> is_mol(slab.atoms.size(), 0);
  for (size_t k = 0; k < m.size(); ++k) {
    Atom a = mol.atoms[k];
    a.pos = m[k];
    a.mol = 0;
    both.atoms.push_back(a);
    is_mol.push_back(1);
  }
  double zmin = 1e300, zmax = -1e300;
  for (const auto& a : both.atoms) zmin = std::min(zmin, a.pos[2]), zmax = std::max(zmax, a.pos[2]);
  for (auto& a : both.atoms) a.pos[2] += o.vacuum / 2.0 - zmin;
  both.cell.c = {0, 0, (zmax - zmin) + o.vacuum};
  std::vector<Vec3> unwrapped;
  for (size_t i = slab.atoms.size(); i < both.atoms.size(); ++i) unwrapped.push_back(both.atoms[i].pos);
  wrap_ase(both, {true, true, false});
  const auto idx = species_order_index(both, o.order);
  System sorted = both;
  sorted.atoms.clear();
  for (int i : idx) { if (is_mol[size_t(i)]) cx.molecule.push_back(int(sorted.atoms.size())); sorted.atoms.push_back(both.atoms[size_t(i)]); }
  sorted.bonds.clear();
  for (size_t i = 0; i < sorted.atoms.size(); ++i) sorted.atoms[i].id = int64_t(i + 1);
  cx.system = sorted;
  // checks on the complex as written
  System slab_part = sorted;
  slab_part.atoms.clear();
  std::vector<Vec3> mpos;
  std::set<int> ms(cx.molecule.begin(), cx.molecule.end());
  for (size_t i = 0; i < sorted.atoms.size(); ++i) (ms.count(int(i)) ? mpos.push_back(sorted.atoms[i].pos) : slab_part.atoms.push_back(sorted.atoms[i]));
  std::vector<int> mel;
  for (int i : cx.molecule) mel.push_back(sorted.atoms[size_t(i)].element);
  std::tie(cx.any, cx.heavy) = contacts(slab_part, mpos, mel, sorted.cell);
  cx.image = image_distance(unwrapped, sorted.cell.a, sorted.cell.b);
  if (cx.any < o.dmin - 0.1) cx.issues.push_back({"ERROR", "closest molecule-slab contact " + fmt("%.2f", cx.any) + " A < " + fmt("%.2f", o.dmin) + " A", {}});
  if (cx.heavy < o.dheavy - 0.1) cx.issues.push_back({"ERROR", "closest heavy-atom molecule-slab contact " + fmt("%.2f", cx.heavy) + " A < " + fmt("%.2f", o.dheavy) + " A", {}});
  if (cx.image < o.image_warn) {
    // the smallest n × n multiple of the slab's cell that passes
    const double na = 1;
    for (int nn = 2; nn <= 12; ++nn) if (image_distance(unwrapped, sorted.cell.a * (nn / na), sorted.cell.b * (nn / na)) >= o.image_warn) { cx.suggested_supercell = nn; break; }
    const std::string tip = cx.suggested_supercell ? " (a cell " + std::to_string(cx.suggested_supercell) + " times larger in-plane passes)" : "";
    if (cx.image < o.image_error) cx.issues.push_back({"ERROR", "molecule overlaps its own periodic image (" + fmt("%.2f", cx.image) + " A): use a larger supercell" + tip, cx.molecule});
    else cx.issues.push_back({"WARN", "molecule only " + fmt("%.2f", cx.image) + " A from its periodic image (lateral interaction): consider a larger supercell" + tip, cx.molecule});
  }
  ValidateOptions vo;
  const auto rep = validate_2d(slab_part, vo);
  for (const auto& f : rep.findings) if (f.level == "ERROR") cx.issues.push_back({"ERROR", "slab: " + f.text, {}});
  for (const auto& f : cx.issues) if (f.level == "ERROR") cx.status = "FAIL";
}

}  // namespace

std::vector<AnchorDef> anchor_library(const std::string& data_dir) {
  const Json j = read_json_file(data_dir + "/sheets/anchors.json");
  std::vector<AnchorDef> out;
  for (const auto& e : j["anchors"].items()) out.push_back({e.text("name"), e.text("smarts")});
  return out;
}

std::map<std::string, std::vector<int>> find_anchors(const System& mol, const std::vector<AnchorDef>& defs) {
  std::map<std::string, std::vector<int>> out;
  const auto nb = mol.neighbours();
  for (const auto& d : defs) {
    const auto sel = select_smarts(mol, d.smarts);
    int first = -1;
    for (size_t i = 0; i < sel.size(); ++i) if (sel[i]) { first = int(i); break; }
    if (first < 0) continue;
    // the connected group of matched atoms holding the first match
    std::vector<int> grp{first};
    std::set<int> seen{first};
    for (size_t h = 0; h < grp.size(); ++h)
      for (uint32_t w : nb[size_t(grp[h])])
        if (sel[w] && seen.insert(int(w)).second) grp.push_back(int(w));
    std::sort(grp.begin(), grp.end());
    out[d.name] = grp;
  }
  return out;
}

std::vector<double> outer_layer_rotations(const System& slab, const Vec3& c, double tol) {
  double zt = -1e300;
  for (const auto& a : slab.atoms) zt = std::max(zt, a.pos[2]);
  std::vector<const Atom*> layer;
  for (const auto& a : slab.atoms) if (a.pos[2] > zt - 0.5) layer.push_back(&a);
  std::vector<double> out;
  for (int deg = 30; deg < 360; deg += 30) {
    const double t = deg * M_PI / 180;
    bool ok = true;
    for (const Atom* a : layer) {
      const Vec3 d = a->pos - c;
      const Vec3 r{c[0] + std::cos(t) * d[0] - std::sin(t) * d[1], c[1] + std::sin(t) * d[0] + std::cos(t) * d[1], a->pos[2]};
      bool hit = false;
      for (const Atom* b : layer)
        if (b->element == a->element && norm(mic_ab(slab.cell, b->pos - r)) < tol) { hit = true; break; }
      if (!hit) { ok = false; break; }
    }
    if (ok) out.push_back(deg);
  }
  return out;
}

System slab_from_relaxed(const System& contcar, int na, int nb, double vacuum, const std::vector<std::string>& order) {
  // one block along z: unwrap across the boundary, the largest gap is the vacuum
  System s = contcar;
  const double cz = s.cell.c[2];
  std::vector<double> z;
  for (const auto& a : s.atoms) z.push_back(a.pos[2] - std::floor(a.pos[2] / cz) * cz);
  std::vector<double> zs = z;
  std::sort(zs.begin(), zs.end());
  double gap = -1, cut = 0;
  for (size_t k = 0; k < zs.size(); ++k) {
    const double g = (k + 1 < zs.size() ? zs[k + 1] : zs[0] + cz) - zs[k];
    if (g > gap) gap = g, cut = zs[k];
  }
  for (size_t i = 0; i < s.atoms.size(); ++i) s.atoms[i].pos[2] = z[i] > cut + 1e-9 ? z[i] : z[i] + cz;
  // every atom within 2 Å of the next along z: one contiguous block
  std::vector<double> zz;
  for (const auto& a : s.atoms) zz.push_back(a.pos[2]);
  std::sort(zz.begin(), zz.end());
  for (size_t k = 1; k < zz.size(); ++k)
    if (zz[k] - zz[k - 1] > 2.0) throw std::runtime_error("the relaxed slab is not one contiguous block (atoms left the surface?)");
  TerminateOptions o;   // only the repeat / vacuum / order parts
  System rep = s;
  rep.atoms.clear();
  for (int m0 = 0; m0 < na; ++m0)
    for (int m1 = 0; m1 < nb; ++m1)
      for (const auto& a0 : s.atoms) { Atom a = a0; a.pos = a.pos + s.cell.a * double(m0) + s.cell.b * double(m1); rep.atoms.push_back(a); }
  rep.cell.a = s.cell.a * double(na);
  rep.cell.b = s.cell.b * double(nb);
  double zmin = 1e300, zmax = -1e300;
  for (const auto& a : rep.atoms) zmin = std::min(zmin, a.pos[2]), zmax = std::max(zmax, a.pos[2]);
  for (auto& a : rep.atoms) a.pos[2] += vacuum / 2.0 - zmin;
  rep.cell.c = {0, 0, (zmax - zmin) + vacuum};
  wrap_ase(rep, {true, true, false});
  rep.bonds.clear();
  (void)o;
  return group_by_species(rep, order);
}

AdsorbSet build_adsorption_set(const System& slab0, const System& mol, const std::map<std::string, std::vector<int>>& anchors, const AdsorbSetOptions& o) {
  if (mol.atoms.empty()) throw std::invalid_argument("adsorption set: no molecule");
  AdsorbSet out;
  out.slab = group_by_species(slab0, o.order);
  out.anchors = anchors;
  const System& slab = out.slab;
  // the free molecule in a box with padding on every side (Γ point)
  {
    System box = group_by_species(mol, o.order);
    Vec3 lo{1e300, 1e300, 1e300}, hi{-1e300, -1e300, -1e300};
    for (const auto& a : box.atoms) for (int k = 0; k < 3; ++k) lo[size_t(k)] = std::min(lo[size_t(k)], a.pos[size_t(k)]), hi[size_t(k)] = std::max(hi[size_t(k)], a.pos[size_t(k)]);
    for (auto& a : box.atoms) a.pos = a.pos - lo + Vec3{o.box_padding, o.box_padding, o.box_padding};
    box.cell.a = {hi[0] - lo[0] + 2 * o.box_padding, 0, 0};
    box.cell.b = {0, hi[1] - lo[1] + 2 * o.box_padding, 0};
    box.cell.c = {0, 0, hi[2] - lo[2] + 2 * o.box_padding};
    box.bonds.clear();
    out.molecule = box;
  }
  std::vector<int> mol_el;
  for (const auto& a : mol.atoms) mol_el.push_back(a.element);
  std::vector<std::string> modes = o.modes;
  if (modes.empty()) {
    for (const auto& [name, at] : anchors) modes.push_back("anchor:" + name);
    modes.push_back("parallel");
  }
  const Vec3 ctr = (slab.cell.a + slab.cell.b) * 0.5;
  double top = -1e300;
  for (const auto& a : slab.atoms) top = std::max(top, a.pos[2]);
  out.symmetry_rotations = outer_layer_rotations(slab, {ctr[0], ctr[1], 0});
  if (!out.symmetry_rotations.empty()) {
    std::string r;
    for (double d : out.symmetry_rotations) r += (r.empty() ? "" : ", ") + fmt("%.0f", d);
    out.notes.push_back("the outer layer is unchanged by rotations of " + r + "° about the placement point: azimuths differing by these are equivalent");
  }
  for (const auto& mode : modes) {
    std::vector<int> anchor;
    std::string tag = mode;
    if (mode.rfind("anchor:", 0) == 0) {
      const auto it = anchors.find(mode.substr(7));
      if (it == anchors.end()) throw std::invalid_argument("no anchor " + mode.substr(7) + " in this molecule");
      anchor = it->second;
      tag = mode.substr(7) + "_down";
    } else if (mode == "upright") {
      if (!anchors.empty()) anchor = anchors.begin()->second;
    } else if (mode != "parallel") throw std::invalid_argument("orientation " + mode + " (anchor:<name>, parallel, upright)");
    for (double az : o.azimuths) {
      AdsorbComplex cx;
      cx.mode = mode;
      cx.azimuth = az;
      cx.name = "complex_" + tag + "_az" + std::to_string(int(std::lround(az)));
      // symmetry-equivalent to an earlier azimuth of this mode?
      for (const auto& prev : out.complexes) {
        if (prev.mode != mode) continue;
        double d = std::fmod(std::fabs(az - prev.azimuth), 360.0);
        for (double r : out.symmetry_rotations)
          if (std::fabs(d - r) < 1e-6 || std::fabs(360 - d - r) < 1e-6) { cx.equivalent_to = prev.name; break; }
        if (!cx.equivalent_to.empty()) break;
      }
      if (o.skip_equivalent && !cx.equivalent_to.empty()) continue;
      // centre over the cell, lower until a contact limit is reached, then one step back
      std::vector<Vec3> base = orient(mol, mode, anchor, az);
      Vec3 mc{0, 0, 0};
      double zlo = 1e300;
      for (const auto& p : base) mc = mc + p, zlo = std::min(zlo, p[2]);
      mc = mc * (1.0 / double(base.size()));
      for (auto& p : base) p = {p[0] + ctr[0] - mc[0], p[1] + ctr[1] - mc[1], p[2]};
      double z = top + 14.0 - zlo;
      const double step = 0.05;
      std::vector<Vec3> m(base.size());
      for (;;) {
        for (size_t k = 0; k < base.size(); ++k) m[k] = base[k] + Vec3{0, 0, z};
        const auto [a, h] = contacts(slab, m, mol_el, slab.cell);
        if (a < o.dmin || h < o.dheavy) { z += step; break; }
        z -= step;
        if (z < -50) throw std::runtime_error("could not place the molecule");
      }
      for (size_t k = 0; k < base.size(); ++k) m[k] = base[k] + Vec3{0, 0, z};
      finish_complex(cx, slab, mol, m, o);
      out.complexes.push_back(std::move(cx));
    }
  }
  return out;
}

std::vector<AuditRow> audit_complexes(const System& slab, const System& molecule, const std::vector<std::pair<std::string, System>>& cxs, double dmin, double dheavy, double dimage) {
  std::map<int, int> sc;
  for (const auto& a : slab.atoms) ++sc[a.element];
  // the free molecule's bonds (covalent radii × 1.2), as element pairs in order
  auto bonds_of = [](const std::vector<Vec3>& p, const std::vector<int>& el) {
    std::multiset<std::pair<int, int>> b;
    for (size_t i = 0; i < p.size(); ++i)
      for (size_t j = i + 1; j < p.size(); ++j)
        if (norm(p[i] - p[j]) < 1.2 * (element(el[i]).covalent + element(el[j]).covalent)) b.insert({std::min(el[i], el[j]), std::max(el[i], el[j])});
    return b;
  };
  std::vector<Vec3> mp;
  std::vector<int> me;
  for (const auto& a : molecule.atoms) mp.push_back(a.pos), me.push_back(a.element);
  const auto ref_bonds = bonds_of(mp, me);
  std::vector<AuditRow> out;
  for (const auto& [name, at] : cxs) {
    AuditRow r;
    r.name = name;
    std::map<int, int> seen;
    std::vector<size_t> mol, sl;
    for (size_t i = 0; i < at.atoms.size(); ++i) (++seen[at.atoms[i].element] > sc[at.atoms[i].element] ? mol : sl).push_back(i);
    r.n_mol = int(mol.size());
    // whole along bonds < 1.9 Å (in-plane minimum image)
    std::vector<Vec3> m(mol.size());
    std::vector<char> done(mol.size(), 0);
    if (!mol.empty()) {
      m[0] = at.atoms[mol[0]].pos;
      done[0] = 1;
      std::deque<size_t> q{0};
      while (!q.empty()) {
        const size_t i = q.front();
        q.pop_front();
        for (size_t j = 0; j < mol.size(); ++j) {
          if (done[j]) continue;
          const Vec3 D = mic_ab(at.cell, at.atoms[mol[j]].pos - at.atoms[mol[i]].pos);
          if (norm(D) < 1.9) { m[j] = m[i] + D; done[j] = 1; q.push_back(j); }
        }
      }
    }
    const bool connected = std::all_of(done.begin(), done.end(), [](char c) { return c != 0; });
    double any = 1e300, heavy = 1e300;
    for (size_t k = 0; k < mol.size(); ++k)
      for (size_t j : sl) {
        const double d = norm(mic_ab(at.cell, at.atoms[j].pos - m[k]));
        any = std::min(any, d);
        if (at.atoms[mol[k]].element != 1 && at.atoms[j].element != 1) heavy = std::min(heavy, d);
      }
    r.any = any, r.heavy = heavy;
    r.image = image_distance(m, at.cell.a, at.cell.b);
    std::vector<int> el;
    for (size_t i : mol) el.push_back(at.atoms[i].element);
    r.intact = connected && int(mol.size()) == int(molecule.atoms.size()) && bonds_of(m, el) == ref_bonds;
    r.ok = r.intact && r.any >= dmin && r.heavy >= dheavy && r.image >= dimage;
    out.push_back(r);
  }
  return out;
}

std::vector<AdsorbComplex> prescreen_complexes(const System& slab0, const System& mol, const AdsorbSetOptions& o, const PrescreenOptions& p, std::vector<double>* energies) {
  const System slab = group_by_species(slab0, o.order);
  // slab (molecule 1, bonds from the crystal) + the molecule above its centre (molecule 2, its own bonds)
  System s = slab;
  s.bonds = crystal_bonds(slab);
  for (auto& a : s.atoms) a.mol = 1;
  double top = -1e300;
  for (const auto& a : slab.atoms) top = std::max(top, a.pos[2]);
  Vec3 c{0, 0, 0};
  double zlo = 1e300;
  for (const auto& a : mol.atoms) c = c + a.pos, zlo = std::min(zlo, a.pos[2]);
  c = c * (1.0 / double(mol.atoms.size()));
  const Vec3 ctr = (slab.cell.a + slab.cell.b) * 0.5;
  const uint32_t off = uint32_t(s.atoms.size());
  for (const auto& a0 : mol.atoms) {
    Atom a = a0;
    a.pos = {a0.pos[0] - c[0] + ctr[0], a0.pos[1] - c[1] + ctr[1], a0.pos[2] - zlo + top + 3.0};
    a.mol = 2;
    s.atoms.push_back(a);
  }
  for (const auto& b : mol.bonds) s.bonds.push_back({b.i + off, b.j + off, b.order});
  s.has_mol = true;
  const ForceField ff = assign_uff(s);
  AdsorptionOptions ao;
  ao.first_mobile_atom = int(off);
  ao.cycles = p.cycles;
  ao.steps = p.steps;
  ao.keep = std::max(1, p.keep) * 4;   // more than asked: near-duplicates are left out below
  ao.seed = p.seed;
  ao.coulomb = false;   // UFF without charges: the contact and dispersion picture, not electrostatics
  ao.z_lo = top + 1.0;
  ao.z_hi = top + p.window;
  AdsorptionReport rep;
  System run = s;
  locate_adsorption(run, ff, ao, &rep);
  std::vector<AdsorbComplex> out;
  std::vector<std::vector<Vec3>> kept;
  for (size_t k = 0; k < rep.configs.size() && int(out.size()) < p.keep; ++k) {
    // the molecule's atoms of this configuration, made whole along its bonds
    std::vector<Vec3> m(mol.atoms.size());
    for (size_t i = 0; i < m.size(); ++i) m[i] = rep.configs[k].positions[off + i];
    std::vector<char> done(m.size(), 0);
    done[0] = 1;
    std::deque<size_t> q{0};
    while (!q.empty()) {
      const size_t i = q.front();
      q.pop_front();
      for (size_t j = 0; j < m.size(); ++j)
        if (!done[j]) {
          const Vec3 d = mic_ab(slab.cell, m[j] - m[i]);
          if (norm(d) < 1.9) m[j] = m[i] + d, done[j] = 1, q.push_back(j);
        }
    }
    // a configuration within 0.5 Å RMS (minimum image) of one already kept is the same minimum
    bool same = false;
    for (const auto& prev : kept) {
      double ss = 0;
      for (size_t i = 0; i < m.size(); ++i) { const Vec3 d = mic_ab(slab.cell, m[i] - prev[i]); ss += dot(d, d); }
      same |= std::sqrt(ss / double(m.size())) < 0.5;
    }
    if (same) continue;
    kept.push_back(m);
    AdsorbComplex cx;
    cx.mode = "prescreen";
    cx.name = "complex_ff" + std::to_string(out.size() + 1);
    finish_complex(cx, slab, mol, m, o);
    out.push_back(std::move(cx));
    if (energies) energies->push_back(rep.configs[k].energy);
  }
  return out;
}

}  // namespace caps
