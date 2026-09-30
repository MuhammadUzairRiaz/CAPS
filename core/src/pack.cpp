// CAPS Pack: rigid-body overlap minimisation with region constraints.
#include "caps/pack.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <deque>
#include <filesystem>
#include <fstream>
#include <numeric>
#include <random>
#include <sstream>

#include "caps/io.hpp"
#include "caps/relax.hpp"
#include "parallel.hpp"

namespace caps {

namespace {

using Mat3 = std::array<std::array<double, 3>, 3>;

Mat3 identity() { return {{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}}; }

Mat3 mul(const Mat3& a, const Mat3& b) {
  Mat3 c{};
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) c[i][j] = a[i][0] * b[0][j] + a[i][1] * b[1][j] + a[i][2] * b[2][j];
  return c;
}

Vec3 apply(const Mat3& m, const Vec3& v) {
  return {m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2], m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
          m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2]};
}

// Rotation matrix of the rotation vector w (Rodrigues).
Mat3 expmap(const Vec3& w) {
  const double th = norm(w);
  if (th < 1e-12) return identity();
  const Vec3 k = w * (1 / th);
  const double c = std::cos(th), s = std::sin(th), v = 1 - c;
  return {{{c + k[0] * k[0] * v, k[0] * k[1] * v - k[2] * s, k[0] * k[2] * v + k[1] * s},
           {k[1] * k[0] * v + k[2] * s, c + k[1] * k[1] * v, k[1] * k[2] * v - k[0] * s},
           {k[2] * k[0] * v - k[1] * s, k[2] * k[1] * v + k[0] * s, c + k[2] * k[2] * v}}};
}

// Left Jacobian of SO(3): exp([w + dw]) ≈ exp([J dw]) exp([w]).
Mat3 left_jacobian(const Vec3& w) {
  const double th = norm(w);
  Mat3 J = identity();
  if (th < 1e-8) return J;
  const double a = (1 - std::cos(th)) / (th * th), b = (th - std::sin(th)) / (th * th * th);
  const Mat3 W = {{{0, -w[2], w[1]}, {w[2], 0, -w[0]}, {-w[1], w[0], 0}}};
  const Mat3 W2 = mul(W, W);
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) J[i][j] += a * W[i][j] + b * W2[i][j];
  return J;
}

Mat3 random_rotation(std::mt19937_64& rng) {
  // uniform unit quaternion (Shoemake)
  std::uniform_real_distribution<double> u(0, 1);
  const double u1 = u(rng), u2 = 2 * M_PI * u(rng), u3 = 2 * M_PI * u(rng);
  const double a = std::sqrt(1 - u1), b = std::sqrt(u1);
  const double x = a * std::sin(u2), y = a * std::cos(u2), z = b * std::sin(u3), w = b * std::cos(u3);
  return {{{1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)},
           {2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)},
           {2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)}}};
}

Mat3 euler_xyz(const Vec3& ang) {
  const Mat3 rx = {{{1, 0, 0}, {0, std::cos(ang[0]), -std::sin(ang[0])}, {0, std::sin(ang[0]), std::cos(ang[0])}}};
  const Mat3 ry = {{{std::cos(ang[1]), 0, std::sin(ang[1])}, {0, 1, 0}, {-std::sin(ang[1]), 0, std::cos(ang[1])}}};
  const Mat3 rz = {{{std::cos(ang[2]), -std::sin(ang[2]), 0}, {std::sin(ang[2]), std::cos(ang[2]), 0}, {0, 0, 1}}};
  return mul(rz, mul(ry, rx));
}

// Penalty P = v² of a region with a safety margin m (inside regions shrink, outside regions grow), and its gradient.
double region_penalty(const Region& R, const Vec3& x, double m, Vec3* grad) {
  Vec3 g{0, 0, 0};
  double p = 0;
  switch (R.kind) {
    case Region::InsideBox:
      for (int k = 0; k < 3; ++k) {
        const double lo = R.a[k] + m, hi = R.b[k] - m;
        if (x[k] < lo) { const double v = lo - x[k]; p += v * v; g[k] -= 2 * v; }
        else if (x[k] > hi) { const double v = x[k] - hi; p += v * v; g[k] += 2 * v; }
      }
      break;
    case Region::OutsideBox: {
      double best = 1e300;
      int axis = -1, dir = 0;
      for (int k = 0; k < 3; ++k) {
        const double lo = R.a[k] - m, hi = R.b[k] + m;
        if (x[k] <= lo || x[k] >= hi) return 0.0;   // outside along this axis: satisfied
        if (x[k] - lo < best) { best = x[k] - lo; axis = k; dir = -1; }
        if (hi - x[k] < best) { best = hi - x[k]; axis = k; dir = +1; }
      }
      p = best * best;
      g[axis] = -2 * best * dir;   // moving along dir reduces the depth
      break;
    }
    case Region::InsideSphere: {
      const Vec3 d = x - R.a;
      const double r = norm(d), v = r - (R.r - m);
      if (v > 0) { p = v * v; if (r > 0) g = d * (2 * v / r); }
      break;
    }
    case Region::OutsideSphere: {
      const Vec3 d = x - R.a;
      const double r = norm(d), v = (R.r + m) - r;
      if (v > 0) { p = v * v; g = r > 0 ? d * (-2 * v / r) : Vec3{-2 * v, 0, 0}; }
      break;
    }
    case Region::InsideCylinder: {
      const Vec3 rel = x - R.a;
      const double t = dot(rel, R.b);
      const Vec3 rr = rel - R.b * t;
      const double dr = norm(rr);
      if (const double v = dr - (R.r - m); v > 0) { p += v * v; if (dr > 0) g = g + rr * (2 * v / dr); }
      if (const double v = m - t; v > 0) { p += v * v; g = g - R.b * (2 * v); }
      if (const double v = t - (R.length - m); v > 0) { p += v * v; g = g + R.b * (2 * v); }
      break;
    }
    case Region::OutsideCylinder: {
      const Vec3 rel = x - R.a;
      const double t = dot(rel, R.b);
      const Vec3 rr = rel - R.b * t;
      const double dr = norm(rr);
      const double d1 = (R.r + m) - dr, d2 = t + m, d3 = (R.length + m) - t;
      if (d1 <= 0 || d2 <= 0 || d3 <= 0) return 0.0;
      if (d1 <= d2 && d1 <= d3) { p = d1 * d1; g = dr > 0 ? rr * (-2 * d1 / dr) : Vec3{0, 0, 0}; }
      else if (d2 <= d3) { p = d2 * d2; g = R.b * (2 * d2); }
      else { p = d3 * d3; g = R.b * (-2 * d3); }
      break;
    }
    case Region::OverPlane: {
      const double v = (R.r + m) - dot(R.a, x);
      if (v > 0) { p = v * v; g = R.a * (-2 * v); }
      break;
    }
    case Region::BelowPlane: {
      const double v = dot(R.a, x) - (R.r - m);
      if (v > 0) { p = v * v; g = R.a * (2 * v); }
      break;
    }
    case Region::InsideEllipsoid:
    case Region::OutsideEllipsoid: {
      // packmol's f = Σ((x−a)/b)² − d, scaled to Å by the smallest semi-axis (√d b_min is the reach along it); the margin
      // shrinks (inside) or grows (outside) the ellipsoid by m along that axis
      const double bmin = std::max(1e-9, std::min({R.b[0], R.b[1], R.b[2]}));
      const double sd = std::sqrt(std::max(0.0, R.r));
      double q = 0;
      Vec3 dq{0, 0, 0};
      for (int k = 0; k < 3; ++k) {
        const double u = (x[k] - R.a[k]) / R.b[k];
        q += u * u;
        dq[k] = 2 * u / R.b[k];
      }
      const double rq = std::sqrt(q);
      const bool inside = R.kind == Region::InsideEllipsoid;
      const double lim = inside ? sd - m / bmin : sd + m / bmin;
      const double v = (inside ? rq - lim : lim - rq) * bmin;   // Å
      if (v > 0 && rq > 1e-12) {
        p = v * v;
        const Vec3 drq = dq * (0.5 / rq);   // ∇√q
        g = drq * (2 * v * bmin * (inside ? 1.0 : -1.0));
      } else if (v > 0) {
        p = v * v;   // at the centre of an outside ellipsoid: any direction
        g = Vec3{-2 * v, 0, 0};
      }
      break;
    }
  }
  if (grad) *grad = g;
  return p;
}

struct Instance {
  int item;
  bool fixed;
  Vec3 c;       // centre
  Mat3 R0;      // anchor rotation
  Vec3 w;       // rotation vector on top of the anchor
};

// Uniform cell grid over a fixed domain (periodic, or clamped at the edges).
struct Grid {
  Vec3 lo{0, 0, 0}, L{1, 1, 1};
  bool periodic = false;
  int n[3] = {1, 1, 1};
  std::vector<uint32_t> start, order, cell_of;
  void setup(const Vec3& lo_, const Vec3& L_, bool per, double cs) {
    lo = lo_;
    L = L_;
    periodic = per;
    for (int k = 0; k < 3; ++k) n[k] = std::max(1, std::min(1024, int(L[k] / cs)));
  }
  int index1(double x, int k) const {
    double f = (x - lo[k]) / L[k];
    if (periodic) f -= std::floor(f);
    const int c = int(f * n[k]);
    return std::clamp(c, 0, n[k] - 1);
  }
  size_t cell(const Vec3& x) const { return (size_t(index1(x[0], 0)) * n[1] + index1(x[1], 1)) * n[2] + index1(x[2], 2); }
  size_t ncells() const { return size_t(n[0]) * n[1] * n[2]; }
  void build(const std::vector<Vec3>& x) {
    const size_t nc = ncells();
    cell_of.resize(x.size());
    start.assign(nc + 1, 0);
    for (size_t i = 0; i < x.size(); ++i) { cell_of[i] = uint32_t(cell(x[i])); ++start[cell_of[i] + 1]; }
    for (size_t c = 0; c < nc; ++c) start[c + 1] += start[c];
    order.resize(x.size());
    std::vector<uint32_t> fill(start.begin(), start.end() - 1);
    for (size_t i = 0; i < x.size(); ++i) order[fill[cell_of[i]]++] = uint32_t(i);
  }
  Vec3 sep(const Vec3& a, const Vec3& b) const {   // b − a, minimum image when periodic
    Vec3 d = b - a;
    if (periodic)
      for (int k = 0; k < 3; ++k) d[k] -= L[k] * std::round(d[k] / L[k]);
    return d;
  }
  // neighbouring cells of (x, y, z), each once
  template <class F>
  void neighbours(int x, int y, int z, F&& f) const {
    int xs[3], ys[3], zs[3], cx = 0, cy = 0, cz = 0;
    auto gather = [&](int v, int k, int* out, int& c) {
      for (int d = -1; d <= 1; ++d) {
        int u = v + d;
        if (u < 0 || u >= n[k]) {
          if (!periodic) continue;
          u = (u + n[k]) % n[k];
        }
        bool dup = false;
        for (int q = 0; q < c; ++q) dup |= out[q] == u;
        if (!dup) out[c++] = u;
      }
    };
    gather(x, 0, xs, cx);
    gather(y, 1, ys, cy);
    gather(z, 2, zs, cz);
    for (int a = 0; a < cx; ++a)
      for (int b = 0; b < cy; ++b)
        for (int c = 0; c < cz; ++c) f((size_t(xs[a]) * n[1] + ys[b]) * n[2] + zs[c]);
  }
};

}  // namespace

double Region::violation(const Vec3& x) const { return std::sqrt(region_penalty(*this, x, 0.0, nullptr)); }

System pack(const std::vector<PackItem>& items, const PackOptions& o, PackReport* rep_out) {
  const auto t0 = std::chrono::steady_clock::now();
  PackReport rep;
  if (items.empty()) throw PackError("nothing to pack");
  if (o.tolerance <= 0) throw PackError("the tolerance must be positive");
  if (o.periodic) {
    if (!o.cell.valid()) throw PackError("periodic packing needs a cell");
    if (o.cell.b[0] != 0 || o.cell.c[0] != 0 || o.cell.c[1] != 0 || o.cell.a[1] != 0 || o.cell.a[2] != 0 || o.cell.b[2] != 0)
      throw PackError("periodic packing supports orthorhombic cells");
  }

  // Templates: coordinates relative to the geometric centre.
  struct Template { std::vector<Vec3> r0; Vec3 centre; };
  std::vector<Template> tpl(items.size());
  for (size_t t = 0; t < items.size(); ++t) {
    const auto& m = items[t].molecule;
    if (m.atoms.empty()) throw PackError("structure '" + items[t].name + "' has no atoms");
    if (items[t].count < 0) throw PackError("structure '" + items[t].name + "': negative number");
    Vec3 c{0, 0, 0};
    for (const auto& a : m.atoms) c = c + a.pos;
    c = c * (1.0 / m.atoms.size());
    tpl[t].centre = c;
    for (const auto& a : m.atoms) tpl[t].r0.push_back(a.pos - c);
  }

  const double tol = o.tolerance + 0.05;   // internal tolerance with a margin, so a zero penalty means d ≥ tolerance
  const double margin = 0.05;
  const double tol2 = tol * tol;
  const double wreg = 4 * tol2;            // region penalty on the same scale as the pair term

  // Domain for the grid.
  Vec3 lo, L;
  if (o.periodic) {
    lo = o.cell.origin;
    L = {o.cell.a[0], o.cell.b[1], o.cell.c[2]};
  } else {
    Vec3 a{1e300, 1e300, 1e300}, b{-1e300, -1e300, -1e300};
    bool any = false;
    auto grow = [&](const Vec3& p) { for (int k = 0; k < 3; ++k) { a[k] = std::min(a[k], p[k]); b[k] = std::max(b[k], p[k]); } any = true; };
    for (size_t t = 0; t < items.size(); ++t) {
      const auto& it = items[t];
      if (it.fixed) {
        for (const auto& p : it.molecule.atoms) grow(p.pos + (it.center ? it.position - tpl[t].centre : Vec3{0, 0, 0}));
        continue;
      }
      bool bounded = false;
      for (const auto& r : it.regions) {
        if (!r.atoms.empty()) continue;
        if (r.kind == Region::InsideEllipsoid) {
          const double sd = std::sqrt(std::max(0.0, r.r));
          grow(r.a - r.b * sd); grow(r.a + r.b * sd); bounded = true;
        }
        if (r.kind == Region::InsideBox) { grow(r.a); grow(r.b); bounded = true; }
        if (r.kind == Region::InsideSphere) { grow(r.a - Vec3{r.r, r.r, r.r}); grow(r.a + Vec3{r.r, r.r, r.r}); bounded = true; }
        if (r.kind == Region::InsideCylinder) {
          grow(r.a - Vec3{r.r, r.r, r.r}); grow(r.a + Vec3{r.r, r.r, r.r});
          const Vec3 e = r.a + r.b * r.length;
          grow(e - Vec3{r.r, r.r, r.r}); grow(e + Vec3{r.r, r.r, r.r});
          bounded = true;
        }
      }
      if (!bounded && it.count > 0)
        throw PackError("structure '" + it.name + "' needs an inside box, sphere or cylinder when the cell is not periodic");
    }
    if (!any) throw PackError("no region to pack into");
    lo = a - Vec3{tol, tol, tol};
    L = b - a + Vec3{2 * tol, 2 * tol, 2 * tol};
  }

  // Instances: fixed first (kept in item order in the output), each molecule a rigid body.
  std::mt19937_64 rng(o.seed);
  std::uniform_real_distribution<double> uni(0, 1);
  std::vector<Instance> inst;
  std::vector<int> first_atom;
  int natoms = 0;
  for (size_t t = 0; t < items.size(); ++t)
    for (int k = 0; k < items[t].count; ++k) {
      Instance in{int(t), items[t].fixed, tpl[t].centre, identity(), {0, 0, 0}};
      if (items[t].fixed) {
        in.R0 = euler_xyz(items[t].angles);
        if (items[t].center) in.c = items[t].position;
        else in.c = tpl[t].centre + items[t].position;
      }
      first_atom.push_back(natoms);
      natoms += int(tpl[t].r0.size());
      inst.push_back(in);
    }
  first_atom.push_back(natoms);
  rep.molecules = int(inst.size());
  rep.atoms = natoms;
  if (inst.empty()) throw PackError("no molecules requested");

  std::vector<int> mol_of(natoms);
  for (size_t m = 0; m < inst.size(); ++m)
    for (int i = first_atom[m]; i < first_atom[m + 1]; ++i) mol_of[i] = int(m);

  std::vector<int> movable;
  for (size_t m = 0; m < inst.size(); ++m)
    if (!inst[m].fixed) movable.push_back(int(m));

  // Positions from the rigid-body state.
  std::vector<Vec3> x(natoms), p(natoms);   // p: rotated body vector (x = c + p)
  auto place = [&](size_t m) {
    const Mat3 R = mul(expmap(inst[m].w), inst[m].R0);
    const auto& r0 = tpl[inst[m].item].r0;
    for (size_t k = 0; k < r0.size(); ++k) {
      p[first_atom[m] + k] = apply(R, r0[k]);
      x[first_atom[m] + k] = inst[m].c + p[first_atom[m] + k];
    }
  };

  auto emit = [&](const std::string& stage, int loop, double f, double dmin, int bad) {
    if (!o.progress) return;
    PackProgress pr{stage, loop, o.max_loops, f, dmin, bad};
    if (!o.progress(pr)) throw PackError("packing cancelled");
  };

  Grid grid;
  grid.setup(lo, L, o.periodic, tol);
  ThreadPool pool(o.threads > 0 ? o.threads : default_threads());

  // Sampling a centre that satisfies the regions (as a point) for item t.
  auto sample_centre = [&](int t) {
    Vec3 a = lo, b = lo + L;
    for (const auto& r : items[t].regions) {
      if (!r.atoms.empty()) continue;   // a region of some atoms does not bound the centre
      Vec3 ra, rb;
      bool box = true;
      if (r.kind == Region::InsideBox) { ra = r.a; rb = r.b; }
      else if (r.kind == Region::InsideSphere) { ra = r.a - Vec3{r.r, r.r, r.r}; rb = r.a + Vec3{r.r, r.r, r.r}; }
      else box = false;
      if (box)
        for (int k = 0; k < 3; ++k) { a[k] = std::max(a[k], ra[k]); b[k] = std::min(b[k], rb[k]); }
    }
    Vec3 c = a;
    for (int tries = 0; tries < 500; ++tries) {
      for (int k = 0; k < 3; ++k) c[k] = a[k] + uni(rng) * std::max(0.0, b[k] - a[k]);
      bool ok = true;
      for (const auto& r : items[t].regions) ok = ok && (!r.atoms.empty() || region_penalty(r, c, 0.0, nullptr) == 0.0);
      if (ok) break;
    }
    return c;
  };

  // Overlap score of molecule m's atoms against atoms already on the grid (with `present` marking placed molecules).
  std::vector<char> present(inst.size(), 0);
  std::vector<std::vector<uint32_t>> dyn(grid.ncells());   // insertion grid, filled as molecules are placed
  auto score = [&](size_t m) {
    double s = 0;
    for (int i = first_atom[m]; i < first_atom[m + 1]; ++i) {
      const Vec3& xi = x[i];
      const int cx = grid.index1(xi[0], 0), cy = grid.index1(xi[1], 1), cz = grid.index1(xi[2], 2);
      grid.neighbours(cx, cy, cz, [&](size_t c) {
        for (uint32_t j : dyn[c]) {
          if (mol_of[j] == int(m)) continue;
          const Vec3 d = grid.sep(xi, x[j]);
          const double r2 = dot(d, d);
          if (r2 < tol2) s += (tol2 - r2) * (tol2 - r2);
        }
      });
      for (const auto& r : items[inst[m].item].regions)
        if (r.applies(int(i) - first_atom[m])) s += wreg * region_penalty(r, xi, margin, nullptr);
    }
    return s;
  };
  auto add_dyn = [&](size_t m) {
    for (int i = first_atom[m]; i < first_atom[m + 1]; ++i) dyn[grid.cell(x[i])].push_back(uint32_t(i));
    present[m] = 1;
  };
  auto remove_dyn = [&](size_t m) {
    for (int i = first_atom[m]; i < first_atom[m + 1]; ++i) {
      auto& v = dyn[grid.cell(x[i])];
      v.erase(std::remove(v.begin(), v.end(), uint32_t(i)), v.end());
    }
    present[m] = 0;
  };
  auto insert = [&](size_t m) {
    double best = 1e300;
    Vec3 bc{0, 0, 0};
    Mat3 bR = identity();
    for (int k = 0; k < std::max(1, o.trials); ++k) {
      inst[m].c = sample_centre(inst[m].item);
      inst[m].R0 = random_rotation(rng);
      inst[m].w = {0, 0, 0};
      place(m);
      const double s = score(m);
      if (s < best) { best = s; bc = inst[m].c; bR = inst[m].R0; }
      if (s == 0) break;
    }
    inst[m].c = bc;
    inst[m].R0 = bR;
    place(m);
  };

  // Stage 1: fixed molecules, then random sequential insertion of the rest.
  for (size_t m = 0; m < inst.size(); ++m)
    if (inst[m].fixed) { place(m); add_dyn(m); }
  for (size_t k = 0; k < movable.size(); ++k) {
    insert(movable[k]);
    add_dyn(movable[k]);
    if (k % 2000 == 1999) emit("insertion", 0, 0, 0, int(movable.size() - k));
  }

  // Objective and gradient over the movable rigid bodies.
  const size_t nv = 6 * movable.size();
  std::vector<double> gx(size_t(natoms) * 3);
  std::vector<std::vector<double>> tg(pool.size());
  std::vector<double> molpen(inst.size());
  auto evaluate = [&](const std::vector<double>& v, std::vector<double>* grad, bool per_molecule) {
    ++rep.evaluations;
    for (size_t k = 0; k < movable.size(); ++k) {
      Instance& in = inst[movable[k]];
      in.c = {v[6 * k], v[6 * k + 1], v[6 * k + 2]};
      in.w = {v[6 * k + 3], v[6 * k + 4], v[6 * k + 5]};
    }
    pool.run(movable.size(), [&](int, size_t b, size_t e) { for (size_t k = b; k < e; ++k) place(movable[k]); });
    grid.build(x);
    const int nth = pool.size();
    std::vector<double> fs(nth, 0.0);
    if (per_molecule) std::fill(molpen.begin(), molpen.end(), 0.0);
    std::vector<std::vector<std::pair<int, double>>> mp(nth);
    pool.run(grid.ncells(), [&](int t, size_t cb, size_t ce) {
      auto& g = tg[t];
      g.assign(size_t(natoms) * 3, 0.0);
      double f = 0;
      for (size_t c = cb; c < ce; ++c) {
        const int cx = int(c / (size_t(grid.n[1]) * grid.n[2])), cy = int((c / grid.n[2]) % grid.n[1]), cz = int(c % grid.n[2]);
        for (uint32_t a = grid.start[c]; a < grid.start[c + 1]; ++a) {
          const uint32_t i = grid.order[a];
          const int mi = mol_of[i];
          grid.neighbours(cx, cy, cz, [&](size_t c2) {
            for (uint32_t b2 = grid.start[c2]; b2 < grid.start[c2 + 1]; ++b2) {
              const uint32_t j = grid.order[b2];
              if (j <= i || mol_of[j] == mi) continue;
              if (inst[mi].fixed && inst[mol_of[j]].fixed) continue;
              const Vec3 d = grid.sep(x[i], x[j]);
              const double r2 = dot(d, d);
              if (r2 >= tol2) continue;
              const double q = tol2 - r2;
              f += q * q;
              if (per_molecule) { mp[t].push_back({mi, 0.5 * q * q}); mp[t].push_back({mol_of[j], 0.5 * q * q}); }
              const double s = 4 * q;   // −∂(q²)/∂x_j = 4q d
              g[3 * j] -= s * d[0]; g[3 * j + 1] -= s * d[1]; g[3 * j + 2] -= s * d[2];
              g[3 * i] += s * d[0]; g[3 * i + 1] += s * d[1]; g[3 * i + 2] += s * d[2];
            }
          });
        }
      }
      fs[t] = f;
    });
    // regions and reduction per movable molecule
    std::vector<double> fr(nth, 0.0);
    pool.run(movable.size(), [&](int t, size_t b, size_t e) {
      for (size_t k = b; k < e; ++k) {
        const int m = movable[k];
        const auto& regs = items[inst[m].item].regions;
        Vec3 gc{0, 0, 0}, tq{0, 0, 0};
        double pm = 0;
        for (int i = first_atom[m]; i < first_atom[m + 1]; ++i) {
          Vec3 gi{0, 0, 0};
          for (int w = 0; w < nth; ++w) gi = gi + Vec3{tg[w][3 * i], tg[w][3 * i + 1], tg[w][3 * i + 2]};
          for (const auto& r : regs) {
            if (!r.applies(i - first_atom[m])) continue;
            Vec3 rg;
            const double pr = wreg * region_penalty(r, x[i], margin, &rg);
            if (pr > 0) { pm += pr; gi = gi + rg * wreg; }
          }
          gc = gc + gi;
          tq = tq + cross(p[i], gi);
        }
        fr[t] += pm;
        if (per_molecule) molpen[m] += pm;
        if (grad) {
          const Mat3 J = left_jacobian(inst[m].w);
          double* G = grad->data() + 6 * k;
          G[0] = gc[0]; G[1] = gc[1]; G[2] = gc[2];
          for (int a = 0; a < 3; ++a) G[3 + a] = J[0][a] * tq[0] + J[1][a] * tq[1] + J[2][a] * tq[2];   // Jᵀ τ
        }
      }
    });
    if (per_molecule)
      for (const auto& v2 : mp)
        for (const auto& [m, val] : v2) molpen[m] += val;
    double f = 0;
    for (int t = 0; t < nth; ++t) f += fs[t] + fr[t];
    return f;
  };

  auto pack_vars = [&]() {
    std::vector<double> v(nv);
    for (size_t k = 0; k < movable.size(); ++k) {
      const Instance& in = inst[movable[k]];
      for (int a = 0; a < 3; ++a) { v[6 * k + a] = in.c[a]; v[6 * k + 3 + a] = in.w[a]; }
    }
    return v;
  };
  // Fold the rotation vectors into the anchors (keeps them small; L-BFGS restarts afterwards).
  auto reanchor = [&]() {
    for (int m : movable) {
      inst[m].R0 = mul(expmap(inst[m].w), inst[m].R0);
      inst[m].w = {0, 0, 0};
    }
  };

  // L-BFGS (m = 10) with backtracking; steps limited to 1 Å in position and 0.3 rad in rotation per molecule.
  auto lbfgs = [&](int maxit, int loop) {
    std::vector<double> v = pack_vars(), g(nv), vn(nv), gn(nv), d(nv);
    double f = evaluate(v, &g, false);
    std::deque<std::vector<double>> S, Y;
    std::deque<double> rho;
    int it = 0;
    for (; it < maxit && f > 1e-14; ++it) {
      std::vector<double> q = g;
      std::vector<double> al(S.size());
      for (int j = int(S.size()) - 1; j >= 0; --j) {
        al[j] = rho[j] * std::inner_product(S[j].begin(), S[j].end(), q.begin(), 0.0);
        for (size_t k = 0; k < nv; ++k) q[k] -= al[j] * Y[j][k];
      }
      double gamma = 1e-3;
      if (!S.empty())
        gamma = std::inner_product(S.back().begin(), S.back().end(), Y.back().begin(), 0.0) /
                std::inner_product(Y.back().begin(), Y.back().end(), Y.back().begin(), 0.0);
      for (auto& z : q) z *= gamma;
      for (size_t j = 0; j < S.size(); ++j) {
        const double bb = rho[j] * std::inner_product(Y[j].begin(), Y[j].end(), q.begin(), 0.0);
        for (size_t k = 0; k < nv; ++k) q[k] += S[j][k] * (al[j] - bb);
      }
      for (size_t k = 0; k < nv; ++k) d[k] = -q[k];
      double slope = std::inner_product(g.begin(), g.end(), d.begin(), 0.0);
      if (slope >= 0) {
        for (size_t k = 0; k < nv; ++k) d[k] = -g[k] * 1e-3;
        slope = std::inner_product(g.begin(), g.end(), d.begin(), 0.0);
        S.clear(); Y.clear(); rho.clear();
      }
      double scale = 1;
      for (size_t k = 0; k < movable.size(); ++k) {
        const double dt = std::sqrt(d[6 * k] * d[6 * k] + d[6 * k + 1] * d[6 * k + 1] + d[6 * k + 2] * d[6 * k + 2]);
        const double dr = std::sqrt(d[6 * k + 3] * d[6 * k + 3] + d[6 * k + 4] * d[6 * k + 4] + d[6 * k + 5] * d[6 * k + 5]);
        if (dt * scale > 1.0) scale = 1.0 / dt;
        if (dr * scale > 0.3) scale = 0.3 / dr;
      }
      double a = scale, fn = 0;
      bool ok = false;
      for (int ls = 0; ls < 25; ++ls) {
        for (size_t k = 0; k < nv; ++k) vn[k] = v[k] + a * d[k];
        fn = evaluate(vn, &gn, false);
        if (fn <= f + 1e-4 * a * slope) { ok = true; break; }
        a *= 0.5;
      }
      if (!ok) { evaluate(v, &g, false); break; }
      std::vector<double> s(nv), y(nv);
      for (size_t k = 0; k < nv; ++k) { s[k] = vn[k] - v[k]; y[k] = gn[k] - g[k]; }
      const double sy = std::inner_product(s.begin(), s.end(), y.begin(), 0.0);
      if (sy > 1e-16) {
        S.push_back(std::move(s)); Y.push_back(std::move(y)); rho.push_back(1 / sy);
        if (S.size() > 10) { S.pop_front(); Y.pop_front(); rho.pop_front(); }
      }
      v.swap(vn);
      g.swap(gn);
      f = fn;
      bool big = false;
      for (int m : movable) big = big || norm(inst[m].w) > 0.5;
      if (big) {   // re-anchor rotations and restart the memory
        reanchor();
        v = pack_vars();
        f = evaluate(v, &g, false);
        S.clear(); Y.clear(); rho.clear();
      }
      if ((it & 15) == 0) emit("optimisation", loop, f, 0, 0);
    }
    rep.iterations += it;
    reanchor();
    return f;
  };

  // Stages 2–3: optimisation rounds; molecules still in violation are partly moved to new places between rounds.
  double f = 0;
  int loop = 0;
  for (; loop < o.max_loops && !movable.empty(); ++loop) {
    f = lbfgs(o.iterations, loop + 1);
    std::vector<double> v = pack_vars();
    f = evaluate(v, nullptr, true);
    std::vector<int> bad;
    for (int m : movable)
      if (molpen[m] > 1e-12) bad.push_back(m);
    emit("optimisation", loop + 1, f, 0, int(bad.size()));
    if (f <= 1e-14 || bad.empty()) break;
    if (loop + 1 == o.max_loops) break;
    // move the worst: a share of the molecules in violation, at least one
    std::sort(bad.begin(), bad.end(), [&](int a, int b) { return molpen[a] > molpen[b]; });
    const size_t nmove = std::max<size_t>(1, size_t(std::ceil(o.move_fraction * bad.size())));
    for (auto& c : dyn) c.clear();
    for (size_t m = 0; m < inst.size(); ++m) add_dyn(m);
    for (size_t k = 0; k < nmove && k < bad.size(); ++k) {
      remove_dyn(bad[k]);
      insert(bad[k]);
      add_dyn(bad[k]);
      ++rep.moved;
    }
  }
  rep.loops = loop + (movable.empty() ? 0 : 1);

  // Assemble the output.
  System out;
  out.title = "CAPS Pack";
  out.source_format = "caps-pack";
  out.cell = o.cell;
  out.unwrapped = true;
  out.has_mol = true;
  if (!o.cell.valid()) {
    out.cell.origin = lo;
    out.cell.a = {L[0], 0, 0};
    out.cell.b = {0, L[1], 0};
    out.cell.c = {0, 0, L[2]};
    out.cell.periodic = {false, false, false};
  }
  // types: each structure's types keep their numbers, offset past the previous structures'
  std::vector<int> type_offset(items.size(), 0);
  int toff = 0;
  for (size_t t = 0; t < items.size(); ++t) {
    type_offset[t] = toff;
    int tmax = 0;
    for (const auto& a : items[t].molecule.atoms) tmax = std::max(tmax, a.type);
    for (const auto& ti : items[t].molecule.types) {
      TypeInfo c = ti;
      c.type += toff;
      out.types.push_back(c);
    }
    toff += tmax;
  }
  // molecules: each placed structure keeps its own (a fixed host of ten chains stays ten molecules); the contact check
  // below is between placed structures
  std::vector<int> placed;
  int64_t next_mol = 1;
  for (size_t m = 0; m < inst.size(); ++m) {
    const auto& mol = items[inst[m].item].molecule;
    const uint32_t base = uint32_t(out.atoms.size());
    int nmol = 0;
    const auto inner = mol.molecules(&nmol);
    for (size_t k = 0; k < mol.atoms.size(); ++k) {
      Atom a = mol.atoms[k];
      a.pos = x[first_atom[m] + k];
      a.id = int64_t(out.atoms.size() + 1);
      a.mol = next_mol + (k < inner.size() ? inner[k] : 0);
      placed.push_back(int(m));
      if (a.type > 0) a.type += type_offset[inst[m].item];
      a.image = {0, 0, 0};
      out.atoms.push_back(a);
      if (mol.has_charges) out.has_charges = true;
    }
    for (const auto& b : mol.bonds) out.bonds.push_back({base + b.i, base + b.j, b.order});
    next_mol += std::max(1, nmol);
  }
  out.bonds_from_file = !out.bonds.empty();
  rep.molecules = int(next_mol - 1);   // molecules in the cell (a fixed host keeps its own), not placed structures

  // Verification with the requested tolerance.
  const auto [dmin, close] = intermolecular_contacts(out, o.tolerance, o.periodic, placed);
  rep.dmin = dmin;
  rep.close_pairs = close;
  double rv = 0;
  for (size_t m = 0; m < inst.size(); ++m) {
    if (inst[m].fixed) continue;
    for (const auto& r : items[inst[m].item].regions)
      for (int i = first_atom[m]; i < first_atom[m + 1]; ++i)
        if (r.applies(i - first_atom[m])) rv = std::max(rv, r.violation(x[i]));
  }
  rep.region_violation = rv;
  rep.success = close == 0 && rv <= 1e-6;
  rep.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  char b[256];
  std::snprintf(b, sizeof b, "%d molecules, %d atoms · %d rounds, %d iterations, %d molecules moved · %.2f s", rep.molecules, rep.atoms, rep.loops,
                rep.iterations, rep.moved, rep.seconds);
  rep.notes.push_back(b);
  std::snprintf(b, sizeof b, "smallest distance between molecules %.3f Å (tolerance %.2f Å) · %d pairs closer · largest region violation %.3g Å",
                rep.dmin, o.tolerance, rep.close_pairs, rep.region_violation);
  rep.notes.push_back(b);
  if (!o.cell.valid() || !o.periodic) rep.notes.push_back("not periodic: distances across the cell faces were not checked");
  emit("verification", rep.loops, f, rep.dmin, rep.close_pairs);
  if (rep_out) *rep_out = rep;
  if (!rep.success) {
    std::snprintf(b, sizeof b, "could not pack within the tolerance: %d pairs closer than %.2f Å (smallest %.3f Å), largest region violation %.3f Å "
                               "after %d rounds; lower the number of molecules or the tolerance, or enlarge the regions",
                  rep.close_pairs, o.tolerance, rep.dmin, rep.region_violation, rep.loops);
    throw PackError(b);
  }
  // stage 3: compression to the target density (the packed cell holds no overlaps to start from)
  if (o.compress_to > 0) {
    if (!o.periodic || !out.cell.valid()) throw PackError("compression needs a periodic cell (pbc)");
    for (const auto& m : inst)
      if (m.fixed) throw PackError("compression would move the fixed structures: pack them at the density you want instead");
    const double rho0 = out.density();
    if (o.compress_to > rho0 * 1.0001) {
      RelaxOptions ro;
      ro.target_density = o.compress_to;
      ro.compress_step = 0.06;
      ro.ftol = 2.0;
      ro.max_iterations = 2000;
      ro.energy.threads = o.threads;
      ro.progress = [&](const RelaxProgress& p) {
        if (!o.progress) return true;
        PackProgress pr{"compression", p.stage_index, p.stages, p.energy, 0, 0};
        return o.progress(pr);
      };
      RelaxReport rr;
      try {
        relax(out, ro, &rr);
      } catch (const RelaxCancelled&) {
        throw PackError("packing cancelled");
      }
      std::snprintf(b, sizeof b, "compressed from %.3f to %.3f g/cm³ (affine steps of 6 %%, push-off and minimisation with %s)", rho0, out.density(),
                    rr.field.c_str());
      if (rep_out) rep_out->notes.push_back(b);
      out.velocities.clear();
    } else if (rep_out) {
      std::snprintf(b, sizeof b, "already at %.3f g/cm³: no compression to %.3f", rho0, o.compress_to);
      rep_out->notes.push_back(b);
    }
  }
  return out;
}

std::pair<double, int> intermolecular_contacts(const System& s, double tolerance, bool periodic) { return intermolecular_contacts(s, tolerance, periodic, {}); }

std::pair<double, int> intermolecular_contacts(const System& s, double tolerance, bool periodic, const std::vector<int>& group) {
  const size_t n = s.atoms.size();
  const auto mol = group.size() == n ? group : s.molecules();
  const double rs = tolerance + 1.0;   // search radius: distances beyond it are reported as ≥ rs
  Grid g;
  Vec3 lo, L;
  const bool per = periodic && s.cell.valid();
  if (per) {
    lo = s.cell.origin;
    L = {s.cell.a[0], s.cell.b[1], s.cell.c[2]};
  } else {
    Vec3 a{1e300, 1e300, 1e300}, b{-1e300, -1e300, -1e300};
    for (const auto& at : s.atoms)
      for (int k = 0; k < 3; ++k) { a[k] = std::min(a[k], at.pos[k]); b[k] = std::max(b[k], at.pos[k]); }
    lo = a - Vec3{1, 1, 1};
    L = b - a + Vec3{2, 2, 2};
  }
  g.setup(lo, L, per, rs);
  std::vector<Vec3> x(n);
  for (size_t i = 0; i < n; ++i) x[i] = s.atoms[i].pos;
  g.build(x);
  double dmin2 = rs * rs;
  int close = 0;
  const double t2 = tolerance * tolerance;
  for (size_t c = 0; c < g.ncells(); ++c) {
    const int cx = int(c / (size_t(g.n[1]) * g.n[2])), cy = int((c / g.n[2]) % g.n[1]), cz = int(c % g.n[2]);
    for (uint32_t a = g.start[c]; a < g.start[c + 1]; ++a) {
      const uint32_t i = g.order[a];
      g.neighbours(cx, cy, cz, [&](size_t c2) {
        for (uint32_t b2 = g.start[c2]; b2 < g.start[c2 + 1]; ++b2) {
          const uint32_t j = g.order[b2];
          if (j <= i || mol[i] == mol[j]) continue;
          const Vec3 d = g.sep(x[i], x[j]);
          const double r2 = dot(d, d);
          dmin2 = std::min(dmin2, r2);
          if (r2 < t2) ++close;
        }
      });
    }
  }
  return {std::sqrt(dmin2), close};
}

std::vector<PackItem> read_packmol_input(const std::string& path, PackOptions& o, std::string* output) {
  std::ifstream f(path);
  if (!f) throw PackError("cannot open " + path);
  std::stringstream ss;
  ss << f.rdbuf();
  return parse_packmol_input(ss.str(), std::filesystem::path(path).parent_path().string(), o, output, path);
}

std::vector<PackItem> parse_packmol_input(const std::string& text, const std::string& base_dir, PackOptions& o, std::string* output,
                                          const std::string& name) {
  std::istringstream in(text);
  const std::filesystem::path base = base_dir;
  const std::string path = name;
  std::vector<PackItem> items;
  PackItem* cur = nullptr;
  std::vector<int> block_atoms;   // inside "atoms … end atoms"
  auto lower = [](std::string w) { std::transform(w.begin(), w.end(), w.begin(), [](unsigned char c) { return char(std::tolower(c)); }); return w; };
  std::string line;
  int lineno = 0;
  o.periodic = false;
  auto bad = [&](const std::string& why) { return PackError(path + ":" + std::to_string(lineno) + ": " + why); };
  while (std::getline(in, line)) {
    ++lineno;
    if (auto h = line.find('#'); h != std::string::npos) line = line.substr(0, h);
    std::istringstream ls(line);
    std::vector<std::string> t;
    for (std::string w; ls >> w;) t.push_back(w);
    if (t.empty()) continue;
    std::string k = t[0];
    std::transform(k.begin(), k.end(), k.begin(), [](unsigned char c) { return char(std::tolower(c)); });
    auto num = [&](size_t i) {
      if (i >= t.size()) throw bad("missing number after '" + t[0] + "'");
      try { return std::stod(t[i]); } catch (...) { throw bad("'" + t[i] + "' is not a number"); }
    };
    // a file name is the rest of the line (quotes around it allowed): folders with spaces ("/Applications/CAPS Studio.app/…")
    auto rest = [&]() {
      std::string r = line.substr(line.find(t[0]) + t[0].size());
      const auto a = r.find_first_not_of(" \t\r"), b = r.find_last_not_of(" \t\r");
      r = a == std::string::npos ? std::string() : r.substr(a, b - a + 1);
      if (r.size() >= 2 && (r.front() == '"' || r.front() == '\'') && r.back() == r.front()) r = r.substr(1, r.size() - 2);
      return r;
    };
    if (k == "structure") {
      if (t.size() < 2) throw bad("structure needs a file");
      const std::string file = rest();
      std::filesystem::path f = file;
      if (f.is_relative()) f = base / f;
      PackItem it;
      it.name = file;
      Trajectory tr = open_file(f.string());
      it.molecule = tr.frame(0);
      items.push_back(std::move(it));
      cur = &items.back();
    } else if (k == "end") {
      if (t.size() > 1 && lower(t[1]) == "atoms") { block_atoms.clear(); continue; }   // end of an atoms block
      cur = nullptr;
      block_atoms.clear();
    } else if (!cur) {
      if (k == "tolerance") o.tolerance = num(1);
      else if (k == "compress") o.compress_to = num(1);   // CAPS: pack loosely, then compress the cell to this density (g/cm³)
      else if (k == "seed") { const double s = num(1); o.seed = s < 0 ? uint64_t(std::random_device{}()) : uint64_t(s); }
      else if (k == "output") { if (output && t.size() > 1) *output = (base / rest()).string(); }
      else if (k == "pbc") {
        Vec3 a{0, 0, 0}, b;
        if (t.size() >= 7) { a = {num(1), num(2), num(3)}; b = {num(4), num(5), num(6)}; }
        else b = {num(1), num(2), num(3)};
        o.cell.origin = a;
        o.cell.a = {b[0] - a[0], 0, 0};
        o.cell.b = {0, b[1] - a[1], 0};
        o.cell.c = {0, 0, b[2] - a[2]};
        o.cell.periodic = {true, true, true};
        o.periodic = true;
      } else if (k == "nloop") o.max_loops = int(num(1));
      else if (k == "maxit") o.iterations = int(num(1));
      else if (k == "filetype" || k == "discale" || k == "precision" || k == "movefrac" || k == "movebadrandom" || k == "add_box_sides" ||
               k == "randominitialpoint" || k == "avoid_overlap" || k == "writeout" || k == "writebad" || k == "check" || k == "sidemax" ||
               k == "connect" || k == "add_amber_ter" || k == "chkgrad" || k == "iprint1" || k == "iprint2" || k == "fbins" || k == "nloop0")
        continue;   // accepted, no effect here
      else throw bad("'" + t[0] + "' is not supported by caps pack");
    } else {
      if (k == "number") cur->count = int(num(1));
      else if (k == "inside" || k == "outside" || k == "over" || k == "above" || k == "below") {
        if (t.size() < 2) throw bad("region kind missing");
        std::string r = t[1];
        std::transform(r.begin(), r.end(), r.begin(), [](unsigned char c) { return char(std::tolower(c)); });
        const bool inside = k == "inside";
        Region g;
        if ((k == "inside" || k == "outside") && r == "box") {
          g.kind = inside ? Region::InsideBox : Region::OutsideBox;
          g.a = {num(2), num(3), num(4)};
          g.b = {num(5), num(6), num(7)};
        } else if ((k == "inside" || k == "outside") && r == "cube") {
          g.kind = inside ? Region::InsideBox : Region::OutsideBox;
          g.a = {num(2), num(3), num(4)};
          const double d = num(5);
          g.b = g.a + Vec3{d, d, d};
        } else if ((k == "inside" || k == "outside") && r == "sphere") {
          g.kind = inside ? Region::InsideSphere : Region::OutsideSphere;
          g.a = {num(2), num(3), num(4)};
          g.r = num(5);
          if (t.size() > 6) throw bad("a sphere takes a centre and a radius (inside ellipsoid for semi-axes)");
        } else if ((k == "inside" || k == "outside") && r == "ellipsoid") {
          g.kind = inside ? Region::InsideEllipsoid : Region::OutsideEllipsoid;
          g.a = {num(2), num(3), num(4)};
          g.b = {num(5), num(6), num(7)};
          g.r = num(8);
          if (g.b[0] <= 0 || g.b[1] <= 0 || g.b[2] <= 0 || g.r <= 0) throw bad("ellipsoid semi-axes and d must be positive");
        } else if ((k == "inside" || k == "outside") && r == "cylinder") {
          g.kind = inside ? Region::InsideCylinder : Region::OutsideCylinder;
          g.a = {num(2), num(3), num(4)};
          Vec3 dir{num(5), num(6), num(7)};
          const double len = norm(dir);
          if (len <= 0) throw bad("cylinder direction is zero");
          g.b = dir * (1 / len);
          g.r = num(8);
          g.length = num(9);
        } else if ((k == "over" || k == "above" || k == "below") && r == "plane") {
          g.kind = k == "below" ? Region::BelowPlane : Region::OverPlane;
          Vec3 nrm{num(2), num(3), num(4)};
          const double len = norm(nrm);
          if (len <= 0) throw bad("plane normal is zero");
          g.a = nrm * (1 / len);
          g.r = num(5) / len;
        } else {
          throw bad("region '" + t[0] + " " + t[1] + "' is not supported by caps pack");
        }
        g.atoms = block_atoms;   // inside "atoms … end atoms": those atoms only
        cur->regions.push_back(g);
      } else if (k == "fixed") {
        cur->fixed = true;
        cur->position = {num(1), num(2), num(3)};
        cur->angles = {num(4), num(5), num(6)};
      } else if (k == "center" || k == "centerofmass") {
        cur->center = true;
      } else if (k == "atoms") {   // the regions until "end atoms" hold for these atoms (1-based in the file)
        block_atoms.clear();
        for (size_t q = 1; q < t.size(); ++q) {
          const int a = int(num(q)) - 1;
          if (a < 0 || size_t(a) >= cur->molecule.atoms.size()) throw bad("atoms: " + t[q] + " is not an atom of " + cur->name);
          block_atoms.push_back(a);
        }
        if (block_atoms.empty()) throw bad("atoms needs atom numbers");
        std::sort(block_atoms.begin(), block_atoms.end());
      } else if (k == "constrain_rotation") {
        throw bad("constrain_rotation is not supported by caps pack: orient the molecule with regions on some of its atoms (atoms … end atoms)");
      } else if (k == "resnumbers" || k == "chain" || k == "segid" || k == "changechains" || k == "radius" || k == "discale" ||
                 k == "maxmove" || k == "nloop" || k == "movebadrandom") {
        continue;
      } else {
        throw bad("'" + t[0] + "' inside a structure is not supported by caps pack");
      }
    }
  }
  if (items.empty()) throw PackError(path + ": no structures");
  for (auto& it : items)
    if (it.fixed && it.count != 1) it.count = 1;
  return items;
}

System insert_molecules(const System& host, const System& guest, int count, const PackOptions& o0, PackReport* report, const std::vector<Region>& regions) {
  if (!host.cell.valid()) throw PackError("inserting molecules needs a periodic cell");
  if (guest.atoms.empty() || count <= 0) return host;
  PackItem h;
  h.name = host.title.empty() ? "structure" : host.title;
  h.molecule = host;
  h.fixed = true;
  PackItem g;
  g.name = guest.title.empty() ? "guest" : guest.title;
  g.molecule = guest;
  g.count = count;
  g.regions = regions;
  PackOptions o = o0;
  o.cell = host.cell;
  o.periodic = true;
  System out = pack({h, g}, o, report);
  // the host's own molecules, names, charges and bond orders; the guests numbered after them
  int64_t top = 0;
  for (const auto& a : host.atoms) top = std::max(top, a.mol);
  const size_t nh = host.atoms.size(), ng = guest.atoms.size();
  for (size_t i = 0; i < out.atoms.size(); ++i) {
    if (i < nh) {
      const Atom& src = host.atoms[i];
      out.atoms[i].mol = src.mol;
      out.atoms[i].name = src.name;
      out.atoms[i].charge = src.charge;
      out.atoms[i].element = src.element;
    } else {
      out.atoms[i].mol = top + 1 + int64_t((i - nh) / ng);
    }
  }
  std::vector<Bond> bonds = host.bonds;
  for (int k = 0; k < count; ++k)
    for (const auto& b : guest.bonds) bonds.push_back({uint32_t(nh + size_t(k) * ng + b.i), uint32_t(nh + size_t(k) * ng + b.j), b.order});
  out.bonds = std::move(bonds);
  out.bonds_from_file = true;
  out.title = host.title;
  out.has_charges = host.has_charges;
  if (report) report->notes.insert(report->notes.begin(), std::to_string(count) + " × " + g.name + " inserted among " + std::to_string(nh) + " atoms held in place");
  return out;
}

}  // namespace caps
