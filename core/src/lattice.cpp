#include "caps/lattice.hpp"

#include <cstdio>
#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>

#include "caps/analysis.hpp"
#include "caps/crystal.hpp"
#include "caps/spacegroup.hpp"

namespace caps {

namespace {

double det3(const Mat3& m) {
  return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
         m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
}
Mat3 mul(const Mat3& x, const Mat3& y) {
  Mat3 r{};
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      for (int k = 0; k < 3; ++k) r[size_t(i)][size_t(j)] += x[size_t(i)][size_t(k)] * y[size_t(k)][size_t(j)];
  return r;
}
Mat3 identity() { return {{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}}; }
// columns of m combine the old vectors
std::array<Vec3, 3> combine(const std::array<Vec3, 3>& v, const Mat3& m) {
  std::array<Vec3, 3> r{};
  for (int j = 0; j < 3; ++j) r[size_t(j)] = v[0] * m[0][size_t(j)] + v[1] * m[1][size_t(j)] + v[2] * m[2][size_t(j)];
  return r;
}
double angle(const Vec3& a, const Vec3& b) { return std::acos(std::clamp(dot(a, b) / (norm(a) * norm(b)), -1.0, 1.0)) * 180 / M_PI; }

}  // namespace

System transform_cell(const System& s, const Mat3& m, double tol) {
  if (!s.cell.valid()) throw std::invalid_argument("changing the lattice needs a periodic cell");
  const double d = det3(m);
  if (std::fabs(d) < 1e-6) throw std::invalid_argument("the matrix is singular (its determinant is 0)");
  const auto nv = combine({s.cell.a, s.cell.b, s.cell.c}, m);
  System out = s;
  out.cell.a = nv[0], out.cell.b = nv[1], out.cell.c = nv[2];
  out.atoms.clear();
  out.bonds.clear();
  out.velocities.clear();
  // the old-lattice translations that can reach the new cell: its corners in old fractional coordinates
  std::array<double, 3> lo{1e9, 1e9, 1e9}, hi{-1e9, -1e9, -1e9};
  for (int c = 0; c < 8; ++c) {
    const Vec3 corner = nv[0] * double(c & 1) + nv[1] * double((c >> 1) & 1) + nv[2] * double((c >> 2) & 1);
    const Vec3 f = s.cell.to_fractional(corner + s.cell.origin);
    for (int k = 0; k < 3; ++k) lo[size_t(k)] = std::min(lo[size_t(k)], f[size_t(k)]), hi[size_t(k)] = std::max(hi[size_t(k)], f[size_t(k)]);
  }
  const double expect = double(s.atoms.size()) * std::fabs(d);
  if (expect > 2e6) throw std::invalid_argument("the new cell would hold more than two million atoms");
  std::vector<Vec3> kept;
  const double eps = 1e-6;
  for (const auto& at : s.atoms) {
    Vec3 f0 = s.cell.to_fractional(at.pos);
    for (int k = 0; k < 3; ++k) f0[size_t(k)] -= std::floor(f0[size_t(k)]);
    for (int i = int(std::floor(lo[0])) - 1; i <= int(std::ceil(hi[0])); ++i)
      for (int j = int(std::floor(lo[1])) - 1; j <= int(std::ceil(hi[1])); ++j)
        for (int k = int(std::floor(lo[2])) - 1; k <= int(std::ceil(hi[2])); ++k) {
          const Vec3 r = s.cell.to_cartesian(f0 + Vec3{double(i), double(j), double(k)});
          Vec3 f = out.cell.to_fractional(r);
          bool in = true;
          for (int q = 0; q < 3; ++q) {
            if (f[size_t(q)] > 1 - eps) f[size_t(q)] -= 1;   // on the far face: the near face's image
            if (f[size_t(q)] < -eps || f[size_t(q)] >= 1 - eps) in = false;
          }
          if (!in) continue;
          const Vec3 p = out.cell.to_cartesian(f);
          bool dup = false;
          for (const auto& q : kept) if (norm(out.cell.minimum_image(q - p)) < tol) { dup = true; break; }
          if (dup) continue;
          kept.push_back(p);
          Atom b = at;
          b.pos = p;
          b.image = {0, 0, 0};
          b.id = int64_t(out.atoms.size() + 1);
          out.atoms.push_back(b);
        }
  }
  if (std::fabs(double(out.atoms.size()) - expect) > 0.5)
    throw std::invalid_argument("the matrix does not map the lattice onto itself: " + std::to_string(out.atoms.size()) + " atoms found, N·|det| = " +
                                std::to_string(expect));
  out.bonds = crystal_bonds(out);
  out.unwrapped = true;
  return out;
}

NiggliResult niggli_reduce(const Cell& cell) {
  NiggliResult r;
  r.m = identity();
  std::array<Vec3, 3> v{cell.a, cell.b, cell.c};
  const double eps = 1e-5 * std::pow(std::fabs(cell.volume()), 2.0 / 3.0);   // relative (Grosse-Kunstleve et al. 2004)
  auto step = [&](const Mat3& m) { v = combine(v, m); r.m = mul(r.m, m); };
  auto sgn = [](double x) { return x > 0 ? 1.0 : -1.0; };
  for (int it = 0; it < 10000; ++it) {
    r.iterations = it + 1;
    const double A = dot(v[0], v[0]), B = dot(v[1], v[1]), C = dot(v[2], v[2]);
    const double xi = 2 * dot(v[1], v[2]), eta = 2 * dot(v[0], v[2]), zeta = 2 * dot(v[0], v[1]);
    // A1
    if (A > B + eps || (std::fabs(A - B) < eps && std::fabs(xi) > std::fabs(eta) + eps)) { step({{{0, -1, 0}, {-1, 0, 0}, {0, 0, -1}}}); continue; }
    // A2
    if (B > C + eps || (std::fabs(B - C) < eps && std::fabs(eta) > std::fabs(zeta) + eps)) { step({{{-1, 0, 0}, {0, 0, -1}, {0, -1, 0}}}); continue; }
    // A3 / A4: all angles acute (+++) or all non-acute (−−−)
    int pos = 0, neg = 0;
    for (double x : {xi, eta, zeta}) pos += x > eps, neg += x < -eps;
    const bool plus = pos == 3 || (pos + neg == 3 && neg == 2 ? false : (pos == 1 && neg == 2 ? false : false)) || (pos == 1 && neg == 0 && false);
    const int zero = 3 - pos - neg;
    if (pos == 3 || (zero == 0 && pos == 1)) {   // the product ξηζ > 0: make all positive
      const double i = xi < -eps ? -1 : 1, j = eta < -eps ? -1 : 1, k = zeta < -eps ? -1 : 1;
      if (i < 0 || j < 0 || k < 0) { step({{{i, 0, 0}, {0, j, 0}, {0, 0, k}}}); continue; }
    } else {   // make all non-positive, keeping the determinant +1
      double f[3] = {1, 1, 1};
      int p = -1;
      const double g[3] = {xi, eta, zeta};
      for (int q = 0; q < 3; ++q) {
        if (g[q] > eps) f[q] = -1;
        else if (!(g[q] < -eps)) p = q;
      }
      if (f[0] * f[1] * f[2] < 0) {
        if (p < 0) throw std::runtime_error("Niggli reduction: no zero angle term to fix the sign (report this cell)");
        f[p] = -1;
      }
      if (f[0] < 0 || f[1] < 0 || f[2] < 0) { step({{{f[0], 0, 0}, {0, f[1], 0}, {0, 0, f[2]}}}); continue; }
    }
    (void)plus;
    // A5
    if (std::fabs(xi) > B + eps || (std::fabs(xi - B) < eps && 2 * eta < zeta - eps) || (std::fabs(xi + B) < eps && zeta < -eps)) {
      step({{{1, 0, 0}, {0, 1, -sgn(xi)}, {0, 0, 1}}});
      continue;
    }
    // A6
    if (std::fabs(eta) > A + eps || (std::fabs(eta - A) < eps && 2 * xi < zeta - eps) || (std::fabs(eta + A) < eps && zeta < -eps)) {
      step({{{1, 0, -sgn(eta)}, {0, 1, 0}, {0, 0, 1}}});
      continue;
    }
    // A7
    if (std::fabs(zeta) > A + eps || (std::fabs(zeta - A) < eps && 2 * xi < eta - eps) || (std::fabs(zeta + A) < eps && eta < -eps)) {
      step({{{1, -sgn(zeta), 0}, {0, 1, 0}, {0, 0, 1}}});
      continue;
    }
    // A8
    const double sum = xi + eta + zeta + A + B;
    if (sum < -eps || (std::fabs(sum) < eps && 2 * (A + eta) + zeta > eps)) {
      step({{{1, 0, 1}, {0, 1, 1}, {0, 0, 1}}});
      continue;
    }
    r.a = norm(v[0]), r.b = norm(v[1]), r.c = norm(v[2]);
    r.alpha = angle(v[1], v[2]), r.beta = angle(v[0], v[2]), r.gamma = angle(v[0], v[1]);
    return r;
  }
  throw std::runtime_error("Niggli reduction did not converge");
}

System find_primitive_cell(const System& s, double tol, int* kout) {
  if (!s.cell.valid()) throw std::invalid_argument("a primitive cell needs a periodic cell");
  const size_t n = s.atoms.size();
  if (n == 0) throw std::invalid_argument("no atoms");
  if (n > 4000) throw std::invalid_argument("more than 4000 atoms: find the primitive cell of a smaller cell");
  std::vector<Vec3> f(n);
  for (size_t i = 0; i < n; ++i) { f[i] = s.cell.to_fractional(s.atoms[i].pos); for (auto& x : f[i]) x -= std::floor(x); }
  auto dist = [&](Vec3 d) { for (auto& x : d) x -= std::round(x); return norm(s.cell.to_cartesian(d) - s.cell.origin); };
  // the rarest element's atoms give the fewest candidate translations
  std::map<int, int> count;
  for (const auto& a : s.atoms) ++count[a.element];
  int rare = s.atoms[0].element;
  for (const auto& [z, c] : count) if (c < count[rare]) rare = z;
  size_t i0 = 0;
  while (s.atoms[i0].element != rare) ++i0;
  std::vector<Vec3> trans{{0, 0, 0}};
  for (size_t j = 0; j < n; ++j) {
    if (j == i0 || s.atoms[j].element != rare) continue;
    Vec3 t = f[j] - f[i0];
    for (auto& x : t) x -= std::floor(x);
    bool ok = true;
    for (size_t a = 0; a < n && ok; ++a) {
      bool hit = false;
      for (size_t b = 0; b < n && !hit; ++b) hit = s.atoms[b].element == s.atoms[a].element && dist(f[a] + t - f[b]) < tol;
      ok = hit;
    }
    if (ok) trans.push_back(t);
  }
  const int k = int(trans.size());
  if (kout) *kout = k;
  if (k == 1) return s;
  if (n % size_t(k) != 0) throw std::runtime_error("the translations found do not divide the atoms evenly: loosen or tighten the tolerance");
  // lattice vectors: the translations and the cell vectors, shortest first
  struct L { Vec3 fr; double len; };
  std::vector<L> vs;
  for (const auto& t : trans)
    for (int i = -1; i <= 1; ++i)
      for (int j = -1; j <= 1; ++j)
        for (int q = -1; q <= 1; ++q) {
          const Vec3 fr = t + Vec3{double(i), double(j), double(q)};
          const double len = norm(s.cell.to_cartesian(fr) - s.cell.origin);
          if (len > 1e-6) vs.push_back({fr, len});
        }
  std::sort(vs.begin(), vs.end(), [](const L& x, const L& y) { return x.len < y.len; });
  const double target = 1.0 / k;   // the primitive volume as a fraction of the cell's
  auto det = [](const Vec3& x, const Vec3& y, const Vec3& z) { return dot(x, cross(y, z)); };
  const size_t lim = std::min<size_t>(vs.size(), 60);
  for (size_t a = 0; a < lim; ++a)
    for (size_t b = a + 1; b < lim; ++b)
      for (size_t c = b + 1; c < lim; ++c) {
        const double d = det(vs[a].fr, vs[b].fr, vs[c].fr);
        if (std::fabs(std::fabs(d) - target) > 1e-6) continue;
        Mat3 m{};
        const Vec3* v3[3] = {&vs[a].fr, &vs[b].fr, d > 0 ? &vs[c].fr : nullptr};
        Vec3 negc = vs[c].fr * -1.0;
        if (!v3[2]) v3[2] = &negc;   // right-handed
        for (int r = 0; r < 3; ++r)
          for (int col = 0; col < 3; ++col) m[size_t(r)][size_t(col)] = (*v3[col])[size_t(r)];
        System out = transform_cell(s, m, tol);
        out.notes.insert(out.notes.begin(), "primitive cell: " + std::to_string(out.atoms.size()) + " atoms (the cell held " + std::to_string(k) + " lattice points)");
        return out;
      }
  throw std::runtime_error("no primitive basis found among the shortest translations");
}

System niggli_cell(const System& s, NiggliResult* result) {
  const NiggliResult r = niggli_reduce(s.cell);
  if (result) *result = r;
  // the reduced matrix is integer and unimodular: the same atoms
  Mat3 m = r.m;
  for (auto& row : m) for (auto& x : row) x = std::round(x);
  System out = transform_cell(s, m);
  out.notes.insert(out.notes.begin(), "Niggli-reduced cell");
  return out;
}

System conventional_cell(const System& s, double tol, ConventionalResult* result) {
  if (!s.cell.valid()) throw std::invalid_argument("a conventional cell needs a periodic cell");
  // from the reduced cell, lattice vectors up to ±2 of each reduced vector
  NiggliResult nr;
  const System red = niggli_cell(find_primitive_cell(s, tol), &nr);   // a supercell or centred cell first comes back primitive
  const std::array<Vec3, 3> rv{red.cell.a, red.cell.b, red.cell.c};
  struct V { std::array<int, 3> n; Vec3 r; double len; };
  std::vector<V> vs;
  for (int i = -2; i <= 2; ++i)
    for (int j = -2; j <= 2; ++j)
      for (int k = -2; k <= 2; ++k) {
        if (i == 0 && j == 0 && k == 0) continue;
        const Vec3 r = rv[0] * double(i) + rv[1] * double(j) + rv[2] * double(k);
        vs.push_back({{i, j, k}, r, norm(r)});
      }
  std::sort(vs.begin(), vs.end(), [](const V& x, const V& y) { return x.len < y.len; });
  if (vs.size() > 80) vs.resize(80);
  auto near = [](double x, double t) { return std::fabs(x - t) < 1.5; };   // degrees
  struct Cand { Mat3 m; double vol, orth; };
  std::vector<Cand> cands;
  std::map<std::string, bool> seen;
  const double v0 = std::fabs(red.cell.volume());
  for (const auto& a : vs)
    for (const auto& b : vs) {
      if (&a == &b) continue;
      const double gab = angle(a.r, b.r);
      if (!(near(gab, 90) || near(gab, 120))) continue;
      for (const auto& c : vs) {
        if (&c == &a || &c == &b) continue;
        const double ga = angle(b.r, c.r), gb = angle(a.r, c.r);
        // orthogonal / hexagonal (c ⟂ a, b) or monoclinic (b ⟂ a, c; any β)
        const bool ortho = near(ga, 90) && near(gb, 90);
        const bool mono = near(gab, 90) && near(ga, 90) && gb > 90.5;
        if (!ortho && !mono) continue;
        Mat3 m{};
        for (int q = 0; q < 3; ++q) m[size_t(q)][0] = a.n[size_t(q)], m[size_t(q)][1] = b.n[size_t(q)], m[size_t(q)][2] = c.n[size_t(q)];
        const double d = det3(m);
        if (d < 0.5 || d > 4.5) continue;   // right-handed, up to 4× (F)
        const double vol = v0 * d;
        char key[160];
        std::snprintf(key, sizeof key, "%.2f %.2f %.2f %.1f %.1f %.1f %.0f", a.len, b.len, c.len, ga, gb, gab, d);
        if (seen.count(key)) continue;
        seen[key] = true;
        cands.push_back({m, vol, std::fabs(ga - 90) + std::fabs(gb - 90) + std::min(std::fabs(gab - 90), std::fabs(gab - 120))});
      }
    }
  std::sort(cands.begin(), cands.end(), [](const Cand& x, const Cand& y) { return x.vol != y.vol ? x.vol < y.vol : x.orth < y.orth; });
  if (cands.size() > 150) cands.resize(150);
  // the reduced cell itself is the triclinic answer
  ConventionalResult best;
  best.m = identity();
  System bestSys = red;
  try {
    const auto f = find_symmetry(red, tol);
    best.hm = f.hm, best.number = f.number, best.operations = f.operations;
  } catch (const std::exception&) {}
  double bestVol = v0;
  for (const auto& cd : cands) {
    System t;
    try { t = transform_cell(red, cd.m); } catch (const std::exception&) { continue; }
    if (t.atoms.size() > 2000) continue;
    SymmetryFound f;
    try { f = find_symmetry(t, tol); } catch (const std::exception&) { continue; }
    ++best.candidates;
    // symmetry per primitive volume: more operations per volume is more symmetric
    const double score = f.operations / cd.vol, cur = best.operations / bestVol;
    if (score > cur * (1 + 1e-9) || (std::fabs(score - cur) <= cur * 1e-9 && f.number > best.number)) {
      best.m = cd.m, best.hm = f.hm, best.number = f.number, best.operations = f.operations;
      bestSys = t;
      bestVol = cd.vol;
    }
  }
  if (result) *result = best;
  bestSys.notes.insert(bestSys.notes.begin(), "conventional cell: " + best.hm + " (No. " + std::to_string(best.number) + "), " + std::to_string(best.operations) +
                                                  " symmetry operations, from " + std::to_string(best.candidates) + " candidate cells");
  return bestSys;
}

System vacuum_slab(const System& s, double vacuum, bool centre, SlabResult* result) {
  if (!s.cell.valid()) throw std::invalid_argument("a vacuum slab needs a periodic cell");
  if (vacuum < 0) throw std::invalid_argument("the vacuum must be positive");
  System out = s;
  const Vec3 n = cross(s.cell.a, s.cell.b) * (1 / norm(cross(s.cell.a, s.cell.b)));   // the surface normal
  const double h = dot(s.cell.c, n);   // perpendicular height of the cell
  // fractional c of every atom, and the slab made contiguous across the largest gap
  std::vector<double> fc;
  for (const auto& a : s.atoms) { double f = s.cell.to_fractional(a.pos)[2]; fc.push_back(f - std::floor(f)); }
  std::vector<double> sorted = fc;
  std::sort(sorted.begin(), sorted.end());
  double gap = 0, cut = 0;
  for (size_t i = 0; i < sorted.size(); ++i) {
    const double next = i + 1 < sorted.size() ? sorted[i + 1] : sorted[0] + 1;
    if (next - sorted[i] > gap) gap = next - sorted[i], cut = sorted[i] + (next - sorted[i]) / 2;
  }
  double lo = 1e9, hi = -1e9;
  for (size_t i = 0; i < out.atoms.size(); ++i) {
    Vec3 f = s.cell.to_fractional(out.atoms[i].pos);
    f[2] = fc[i] < cut ? fc[i] + 1 : fc[i];   // atoms above the gap come first
    f[2] -= cut;
    out.atoms[i].pos = s.cell.to_cartesian(f);
    const double z = dot(out.atoms[i].pos - s.cell.origin, n);
    lo = std::min(lo, z), hi = std::max(hi, z);
  }
  const double thick = s.atoms.empty() ? 0 : hi - lo;
  const double newh = thick + vacuum;
  out.cell.c = s.cell.c * (newh / h);
  // centre the slab in the new height (or put its bottom at the cell's floor)
  const double shift = centre ? (newh - thick) / 2 - lo : -lo;
  const Vec3 along = s.cell.c * (1 / h);   // moving along c by one Å of height
  for (auto& a : out.atoms) a.pos = a.pos + along * shift;
  out.cell.periodic = {true, true, true};
  out.unwrapped = true;
  out.bonds = crystal_bonds(out);
  if (result) result->thickness = thick, result->vacuum = vacuum;
  char b[160];
  std::snprintf(b, sizeof b, "vacuum slab: %.2f Å of atoms and %.2f Å of vacuum along c", thick, vacuum);
  out.notes.insert(out.notes.begin(), b);
  return out;
}

System nanowire(const System& bulk, const WireOptions& o, WireResult* result) {
  if (!bulk.cell.valid()) throw std::invalid_argument("a nanowire is cut from a crystal (a periodic cell)");
  if (o.uvw == std::array<int, 3>{0, 0, 0}) throw std::invalid_argument("the wire direction [uvw] is zero");
  if (o.radius <= 0 || o.repeats < 1) throw std::invalid_argument("the radius and the repeats must be positive");
  const Vec3 t = bulk.cell.a * double(o.uvw[0]) + bulk.cell.b * double(o.uvw[1]) + bulk.cell.c * double(o.uvw[2]);
  const double period = norm(t);
  const Vec3 e = t * (1 / period);
  // two directions across the wire
  Vec3 p = std::fabs(e[0]) < 0.9 ? Vec3{1, 0, 0} : Vec3{0, 1, 0};
  p = p - e * dot(p, e);
  p = p * (1 / norm(p));
  const Vec3 q = cross(e, p);
  const double L = period * o.repeats;
  const double R = o.radius;
  auto inside = [&](double x, double y) {
    if (o.shape == "square") return std::fabs(x) <= R && std::fabs(y) <= R;
    if (o.shape == "hexagonal") {   // a hexagon of inscribed radius R, flat sides across y
      for (int k = 0; k < 6; ++k) {
        const double th = k * M_PI / 3;
        if (x * std::cos(th) + y * std::sin(th) > R + 1e-9) return false;
      }
      return true;
    }
    if (o.shape != "cylinder") throw std::invalid_argument("wire shape cylinder, hexagonal or square");
    return x * x + y * y <= R * R;
  };
  // lattice translations covering the wire's bounding box
  const double reach = std::sqrt(2.0) * R * 1.2 + L;
  std::array<int, 3> nmax{};
  const Vec3 rec[3] = {cross(bulk.cell.b, bulk.cell.c), cross(bulk.cell.c, bulk.cell.a), cross(bulk.cell.a, bulk.cell.b)};
  const double V = bulk.cell.volume();
  for (int k = 0; k < 3; ++k) nmax[size_t(k)] = int(std::ceil(reach * norm(rec[k]) / std::fabs(V))) + 1;
  if (double(nmax[0]) * nmax[1] * nmax[2] * 8 * double(bulk.atoms.size()) > 5e7) throw std::invalid_argument("the wire is too large");
  System out;
  out.title = bulk.title.empty() ? "nanowire" : bulk.title + " nanowire";
  out.types = bulk.types;
  const double side = 2 * R + 2 * o.vacuum;
  out.cell.a = p * side, out.cell.b = q * side, out.cell.c = e * L;
  out.cell.periodic = {true, true, true};
  const Vec3 centre = (p + q) * (side / 2);
  std::vector<Vec3> kept;
  for (const auto& at : bulk.atoms) {
    const Vec3 f0 = bulk.cell.to_fractional(at.pos);
    for (int i = -nmax[0]; i <= nmax[0]; ++i)
      for (int j = -nmax[1]; j <= nmax[1]; ++j)
        for (int k = -nmax[2]; k <= nmax[2]; ++k) {
          const Vec3 r = bulk.cell.to_cartesian(f0 + Vec3{double(i), double(j), double(k)}) - bulk.cell.origin;
          const double z = dot(r, e);
          if (z < -1e-6 || z >= L - 1e-6) continue;
          const double x = dot(r, p), y = dot(r, q);
          if (!inside(x, y)) continue;
          const Vec3 pos = centre + p * x + q * y + e * z;
          bool dup = false;
          for (const auto& w : kept) if (norm(out.cell.minimum_image(w - pos)) < 0.05) { dup = true; break; }
          if (dup) continue;
          kept.push_back(pos);
          Atom b = at;
          b.pos = pos;
          b.id = int64_t(out.atoms.size() + 1);
          b.image = {0, 0, 0};
          out.atoms.push_back(b);
        }
  }
  if (out.atoms.empty()) throw std::invalid_argument("no atoms inside the wire: raise the radius");
  out.bonds = crystal_bonds(out);
  out.unwrapped = true;
  if (result) result->period = period, result->atoms = out.atoms.size();
  char b[200];
  std::snprintf(b, sizeof b, "%s nanowire along [%d %d %d] (period %.3f Å × %d), radius %.1f Å, %.1f Å vacuum around", o.shape.c_str(), o.uvw[0], o.uvw[1], o.uvw[2],
                period, o.repeats, R, o.vacuum);
  out.notes.push_back(b);
  return out;
}


System cut_cluster(const System& s0, const Vec3& centre, double radius, bool any_atom, ClusterResult* result) {
  if (!s0.cell.valid()) throw std::invalid_argument("the structure has no periodic cell");
  System s = s0;
  if (!s.unwrapped) make_molecules_whole(s);
  int nm = 0;
  const auto mol = s.molecules(&nm);
  std::vector<std::vector<uint32_t>> members(size_t(std::max(nm, 0)));
  for (uint32_t i = 0; i < mol.size(); ++i) members[size_t(mol[i])].push_back(i);
  const double big = std::cbrt(s.cell.volume());
  std::vector<char> keep(s.atoms.size(), 0);
  ClusterResult R;
  for (const auto& m : members) {
    Vec3 c{0, 0, 0};
    for (uint32_t i : m) c = c + s.atoms[i].pos;
    c = c * (1.0 / double(m.size()));
    double ext = 0;   // the molecule's own size: one larger than the cell is a network
    for (uint32_t i : m) ext = std::max(ext, norm(s.atoms[i].pos - c));
    if (ext > big) {   // a network: each atom to its own nearest image
      for (uint32_t i : m) {
        const Vec3 d = s.cell.minimum_image(s.atoms[i].pos - centre);
        s.atoms[i].pos = centre + d;
        if (radius <= 0 || norm(d) <= radius) keep[i] = 1;
      }
      continue;
    }
    const Vec3 shift = s.cell.minimum_image(c - centre) - (c - centre);
    bool in = radius <= 0 || norm(c + shift - centre) <= radius;
    for (uint32_t i : m) {
      s.atoms[i].pos = s.atoms[i].pos + shift;
      if (any_atom && norm(s.atoms[i].pos - centre) <= radius) in = true;
    }
    if (in) {
      for (uint32_t i : m) keep[i] = 1;
      ++R.molecules;
    } else {
      ++R.dropped;
    }
  }
  // the kept atoms, the bonds among them
  std::vector<int64_t> newi(s.atoms.size(), -1);
  System out = s;
  out.atoms.clear(), out.bonds.clear(), out.velocities.clear();
  for (size_t i = 0; i < s.atoms.size(); ++i)
    if (keep[i]) newi[i] = int64_t(out.atoms.size()), out.atoms.push_back(s.atoms[i]);
  for (const auto& b : s.bonds)
    if (newi[b.i] >= 0 && newi[b.j] >= 0) out.bonds.push_back({uint32_t(newi[b.i]), uint32_t(newi[b.j]), b.order});
  out.cell = Cell{};
  out.unwrapped = true;
  for (const auto& a : out.atoms) R.radius = std::max(R.radius, norm(a.pos - centre));
  R.atoms = out.atoms.size();
  if (result) *result = R;
  return out;
}

}  // namespace caps
