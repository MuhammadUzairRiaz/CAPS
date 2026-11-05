// CAPS layer stacks (see layers.hpp).
#include "caps/layers.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <stdexcept>

namespace caps {

namespace {

// a and b in the xy plane (any angle between them, right-handed), c along z
bool planar(const Cell& c) {
  const double tol = 1e-3;
  return c.valid() && std::fabs(c.a[2]) < tol && std::fabs(c.b[2]) < tol && std::fabs(c.c[0]) < tol && std::fabs(c.c[1]) < tol && c.c[2] > 0 &&
         c.a[0] * c.b[1] - c.a[1] * c.b[0] > tol;
}

// in-plane fractional coordinates of (x, y) relative to the origin in the cell (a, b)
void frac2(const Cell& c, double x, double y, double& u, double& v) {
  const double det = c.a[0] * c.b[1] - c.a[1] * c.b[0];
  u = (x * c.b[1] - y * c.b[0]) / det;
  v = (c.a[0] * y - c.a[1] * x) / det;
}

// the layer turned about z so a lies along +x, its molecules whole along z (centres back in the cell) and its atoms
// wrapped in the plane
System prepared(const System& in) {
  System s = in;
  const Vec3 o = s.cell.origin;
  const double th = std::atan2(s.cell.a[1], s.cell.a[0]);
  if (std::fabs(th) > 1e-12) {
    const double c = std::cos(-th), sn = std::sin(-th);
    auto rot = [&](const Vec3& v) { return Vec3{c * v[0] - sn * v[1], sn * v[0] + c * v[1], v[2]}; };
    for (auto& a : s.atoms) a.pos = o + rot(a.pos - o);
    s.cell.a = rot(s.cell.a), s.cell.b = rot(s.cell.b);
    s.cell.a[1] = 0;
  }
  const double Lz = s.cell.c[2];
  const size_t n = s.atoms.size();
  const auto nb = s.neighbours();
  std::vector<char> seen(n, 0);
  for (size_t r = 0; r < n; ++r) {   // z unwrapped along the bonds, one molecule at a time
    if (seen[r]) continue;
    std::vector<size_t> comp{r}, stack{r};
    seen[r] = 1;
    while (!stack.empty()) {
      const size_t u = stack.back();
      stack.pop_back();
      for (uint32_t v : nb[u]) {
        if (seen[v]) continue;
        double dz = s.atoms[v].pos[2] - s.atoms[u].pos[2];
        dz -= Lz * std::round(dz / Lz);
        s.atoms[v].pos[2] = s.atoms[u].pos[2] + dz;
        seen[v] = 1;
        stack.push_back(v), comp.push_back(v);
      }
    }
    double zc = 0;
    for (size_t i : comp) zc += s.atoms[i].pos[2];
    zc /= double(comp.size());
    const double shift = -Lz * std::floor((zc - o[2]) / Lz);
    for (size_t i : comp) s.atoms[i].pos[2] += shift;
  }
  for (auto& a : s.atoms) {
    double u, v;
    frac2(s.cell, a.pos[0] - o[0], a.pos[1] - o[1], u, v);
    u -= std::floor(u), v -= std::floor(v);
    a.pos[0] = o[0] + u * s.cell.a[0] + v * s.cell.b[0];
    a.pos[1] = o[1] + u * s.cell.a[1] + v * s.cell.b[1];
  }
  return s;
}

// repeats n of an edge L to come closest to the target T (the strain T / (n L) − 1 smallest)
int best_repeat(double L, double T, int max_repeat) {
  int best = 1;
  for (int k = 1; k <= max_repeat; ++k)
    if (std::fabs(T / (k * L) - 1) < std::fabs(T / (best * L) - 1)) best = k;
  return best;
}

}  // namespace

System stack_layers(const std::vector<StackLayerInput>& layers, const StackOptions& o, StackReport* report) {
  if (layers.size() < 2) throw std::invalid_argument("a stack needs two layers or more");
  for (const auto& l : layers) {
    if (!l.system || l.system->atoms.empty()) throw std::invalid_argument((l.name.empty() ? "a layer" : l.name) + " has no atoms");
    if (!planar(l.system->cell))
      throw std::invalid_argument((l.name.empty() ? "a layer" : l.name) + ": a and b must lie in the xy plane and c along z (a slab's cell)");
  }
  StackReport rep;
  const double max_rep = std::max(1, o.max_repeat);
  std::vector<System> prep;
  for (const auto& l : layers) prep.push_back(prepared(*l.system));
  // a flipped layer (180° about x): y → b_y − y, z → −z; its lattice is spanned by a and (−b_x, b_y)
  for (size_t li = 0; li < layers.size(); ++li)
    if (layers[li].flip) {
      System& s = prep[li];
      const Vec3 o0 = s.cell.origin;
      for (auto& a : s.atoms) a.pos[1] = 2 * o0[1] + s.cell.b[1] - a.pos[1], a.pos[2] = -a.pos[2];
      s.cell.b[0] = -s.cell.b[0];
      for (auto& a : s.atoms) {
        double u, v;
        frac2(s.cell, a.pos[0] - o0[0], a.pos[1] - o0[1], u, v);
        u -= std::floor(u), v -= std::floor(v);
        a.pos[0] = o0[0] + u * s.cell.a[0] + v * s.cell.b[0];
        a.pos[1] = o0[1] + u * s.cell.a[1] + v * s.cell.b[1];
      }
    }
  auto len = [](const Vec3& v) { return std::hypot(v[0], v[1]); };
  // the lateral cell: the first layer's, repeated (match both) where that lowers the worst strain of the others
  const Cell& c0 = prep[0].cell;
  double best_cost = 1e300;
  int n0[2] = {1, 1};
  for (int dir = 0; dir < 2; ++dir) {
    const double L0 = dir == 0 ? len(c0.a) : len(c0.b);
    best_cost = 1e300;
    const int top = o.match == "both" ? int(max_rep) : 1;
    for (int k = 1; k <= top; ++k) {
      const double T = k * L0;
      if (k > 1 && T > o.max_cell) break;
      double worst = 0;
      for (size_t li = 1; li < layers.size(); ++li) {
        const double L = dir == 0 ? len(prep[li].cell.a) : len(prep[li].cell.b);
        const int m = best_repeat(L, T, int(max_rep));
        worst = std::max(worst, std::fabs(T / (m * L) - 1));
      }
      if (worst < best_cost - 1e-4) best_cost = worst, n0[dir] = k;   // the smallest repeat that does as well
    }
  }
  Vec3 TA = c0.a * double(n0[0]), TB = c0.b * double(n0[1]);
  TA[2] = TB[2] = 0;
  if (o.match == "average") {   // every layer's repeats to that cell, then the mean of their own lateral vectors
    Vec3 sa = TA, sb = TB;
    for (size_t li = 1; li < layers.size(); ++li) {
      const Cell& c = prep[li].cell;
      sa = sa + c.a * double(best_repeat(len(c.a), len(TA), int(max_rep)));
      sb = sb + c.b * double(best_repeat(len(c.b), len(TB), int(max_rep)));
    }
    TA = sa * (1.0 / double(layers.size())), TB = sb * (1.0 / double(layers.size()));
    TA[2] = TB[2] = 0;
  }
  const double A = len(TA), B = len(TB);
  auto angle = [&](const Vec3& u, const Vec3& v) { return std::acos(std::clamp((u[0] * v[0] + u[1] * v[1]) / (len(u) * len(v)), -1.0, 1.0)) * 180 / 3.14159265358979323846; };
  System out;
  out.title = "layer stack";
  std::map<std::pair<std::string, long long>, int> type_of;   // (label, mass × 1e4) → merged type
  int64_t mol_off = 0, id = 0;
  double z = o.vacuum > 0 ? 0.0 : 0.5 * o.gap;
  for (size_t li = 0; li < layers.size(); ++li) {
    const System& s = prep[li];
    const Vec3 la = s.cell.a, lb = s.cell.b;
    const int na = li == 0 ? n0[0] : best_repeat(len(la), A, int(max_rep));
    const int nb = li == 0 ? n0[1] : best_repeat(len(lb), B, int(max_rep));
    // F maps the layer's repeated lateral vectors (na a, nb b) onto the stack's (TA, TB)
    const double p00 = na * la[0], p10 = na * la[1], p01 = nb * lb[0], p11 = nb * lb[1];
    const double dp = p00 * p11 - p01 * p10;
    const double F00 = (TA[0] * p11 - TB[0] * p10) / dp, F01 = (TB[0] * p00 - TA[0] * p01) / dp;
    const double F10 = (TA[1] * p11 - TB[1] * p10) / dp, F11 = (TB[1] * p00 - TA[1] * p01) / dp;
    const double J = F00 * F11 - F01 * F10;   // the layer's change of area
    double zlo = 1e300, zhi = -1e300;
    for (const auto& a : s.atoms) zlo = std::min(zlo, a.pos[2]), zhi = std::max(zhi, a.pos[2]);
    const double zs = o.keep_volume ? 1.0 / J : 1.0;   // the thickness scaled so the layer keeps its volume
    // the layer's types into the merged table
    std::map<int, int> tmap;
    for (const auto& t : s.types) {
      const auto key = std::make_pair(t.label, static_cast<long long>(std::llround(t.mass * 1e4)));
      auto it = type_of.find(key);
      if (it == type_of.end()) {
        const int nt = int(out.types.size()) + 1;
        it = type_of.emplace(key, nt).first;
        out.types.push_back({nt, t.mass, t.label});
      }
      tmap[t.type] = it->second;
    }
    // molecules: renumbered in the layer's order, each copy its own
    const auto molidx = s.molecules();
    int nmol = 0;
    for (int x : molidx) nmol = std::max(nmol, x + 1);
    const Vec3 org = s.cell.origin;
    for (int ia = 0; ia < na; ++ia)
      for (int ib = 0; ib < nb; ++ib) {
        const uint32_t base = uint32_t(out.atoms.size());
        for (size_t i = 0; i < s.atoms.size(); ++i) {
          Atom a = s.atoms[i];
          const double px = a.pos[0] - org[0] + ia * la[0] + ib * lb[0], py = a.pos[1] - org[1] + ia * la[1] + ib * lb[1];
          a.pos = {F00 * px + F01 * py + layers[li].shift_x, F10 * px + F11 * py + layers[li].shift_y, z + (a.pos[2] - zlo) * zs};
          a.image = {0, 0, 0};
          a.mol = mol_off + int64_t(molidx[i]) + 1;
          a.id = ++id;
          if (a.type != 0) { auto t = tmap.find(a.type); a.type = t == tmap.end() ? 0 : t->second; }
          out.atoms.push_back(a);
        }
        for (const auto& b : s.bonds) out.bonds.push_back({b.i + base, b.j + base, b.order});
        mol_off += nmol;
      }
    StackLayerReport lr;
    lr.name = (layers[li].name.empty() ? "layer " + std::to_string(li + 1) : layers[li].name) + (layers[li].flip ? " (flipped)" : "");
    lr.na = na, lr.nb = nb;
    lr.strain_a = A / (na * len(la)) - 1, lr.strain_b = B / (nb * len(lb)) - 1;
    lr.shear_deg = angle(TA, TB) - angle(la, lb);
    lr.thickness_scale = zs;
    lr.z_lo = z, lr.z_hi = z + (zhi - zlo) * zs;
    lr.atoms = s.atoms.size() * size_t(na * nb);
    rep.layers.push_back(lr);
    out.has_charges = out.has_charges || s.has_charges;
    z = lr.z_hi + o.gap;
  }
  const double C = o.vacuum > 0 ? z - o.gap + o.vacuum : z - 0.5 * o.gap;
  out.cell.origin = {0, 0, 0};
  out.cell.a = TA;
  out.cell.b = TB;
  out.cell.c = {0, 0, C};
  out.has_mol = true;
  out.bonds_from_file = true;
  // bonds across the plane follow the minimum image; the positions stay in the cell there
  for (auto& a : out.atoms) {
    double u, v;
    frac2(out.cell, a.pos[0], a.pos[1], u, v);
    u -= std::floor(u), v -= std::floor(v);
    a.pos[0] = u * TA[0] + v * TB[0], a.pos[1] = u * TA[1] + v * TB[1];
  }
  rep.a = A, rep.b = B, rep.c = C;
  char b[240];
  for (const auto& l : rep.layers) {
    std::snprintf(b, sizeof b, "%s: %d × %d · strain %+.2f %% × %+.2f %%%s%s · z %.1f–%.1f Å · %zu atoms", l.name.c_str(), l.na, l.nb, 100 * l.strain_a, 100 * l.strain_b,
                  std::fabs(l.shear_deg) > 0.005 ? (" · shear " + std::to_string(l.shear_deg).substr(0, 5) + "°").c_str() : "",
                  std::fabs(l.thickness_scale - 1) > 1e-6 ? (" · thickness × " + std::to_string(l.thickness_scale).substr(0, 6)).c_str() : "", l.z_lo, l.z_hi, l.atoms);
    rep.notes.push_back(b);
  }
  std::snprintf(b, sizeof b, "stack of %zu layers: cell %.2f × %.2f × %.2f Å · gap %.1f Å · %s", layers.size(), A, B, C, o.gap,
                o.vacuum > 0 ? ("vacuum " + std::to_string(int(std::lround(o.vacuum))) + " Å on top").c_str() : "periodic in z");
  rep.notes.insert(rep.notes.begin(), b);
  double worst = 0;
  for (const auto& l : rep.layers) worst = std::max({worst, std::fabs(l.strain_a), std::fabs(l.strain_b), std::fabs(l.shear_deg) * 3.14159265358979323846 / 180});
  if (worst > 0.05) rep.notes.push_back("a layer is strained by more than 5 %: relax (or pick other repeats) before dynamics; crystals should not be strained");
  rep.notes.push_back("the gaps are a start: relax, then equilibrate the stack at constant pressure along z");
  if (report) *report = rep;
  return out;
}

}  // namespace caps
