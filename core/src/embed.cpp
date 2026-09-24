// Distance-bounds embedding and the molecule builder (see molecule.hpp).
#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <random>
#include <set>

#include "caps/elements.hpp"
#include "caps/ffdef.hpp"
#include "caps/molecule.hpp"
#include "caps/relax.hpp"
#include "caps/typing.hpp"
#include "caps/uff.hpp"

namespace caps {
namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr int D = 4;   // embedding dimensions (the fourth is squeezed out)

struct Pair {
  int i, j;
  double lb, ub;
  double w;
};

struct Chiral {
  int c, n1, n2, n3;
  int sign;       // +1 or −1 for V = (n1−c)·((n2−c)×(n3−c))
  double vmin;
};

struct Planar {
  int c, n0, n1, n2;   // for a double bond: n0 on c, c=n1, n2 on n1 (the four atoms coplanar)
};

// Limited-memory BFGS with a backtracking line search; returns the final value.
double lbfgs(std::vector<double>& x, const std::function<double(const std::vector<double>&, std::vector<double>&)>& fg, int maxit,
             double gtol) {
  const size_t n = x.size();
  const int m = 8;
  std::vector<std::vector<double>> S, Y;
  std::vector<double> rho;
  std::vector<double> g(n), gn(n), d(n), xn(n);
  double f = fg(x, g);
  for (int it = 0; it < maxit; ++it) {
    double gmax = 0;
    for (double v : g) gmax = std::max(gmax, std::fabs(v));
    if (gmax < gtol) break;
    // two-loop recursion
    d = g;
    std::vector<double> alpha(S.size());
    for (int k = int(S.size()) - 1; k >= 0; --k) {
      double a = 0;
      for (size_t t = 0; t < n; ++t) a += S[size_t(k)][t] * d[t];
      a *= rho[size_t(k)];
      alpha[size_t(k)] = a;
      for (size_t t = 0; t < n; ++t) d[t] -= a * Y[size_t(k)][t];
    }
    if (!S.empty()) {
      double sy = 0, yy = 0;
      for (size_t t = 0; t < n; ++t) sy += S.back()[t] * Y.back()[t], yy += Y.back()[t] * Y.back()[t];
      const double gamma = yy > 0 ? sy / yy : 1;
      for (double& v : d) v *= gamma;
    } else {
      double gn2 = 0;
      for (double v : g) gn2 += v * v;
      const double sc = 0.1 / std::sqrt(gn2 + 1e-30);
      for (double& v : d) v *= sc;
    }
    for (size_t k = 0; k < S.size(); ++k) {
      double b = 0;
      for (size_t t = 0; t < n; ++t) b += Y[k][t] * d[t];
      b *= rho[k];
      for (size_t t = 0; t < n; ++t) d[t] += S[k][t] * (alpha[k] - b);
    }
    double slope = 0;
    for (size_t t = 0; t < n; ++t) slope -= g[t] * d[t];
    if (slope >= 0) {   // not a descent direction: restart along −g
      S.clear(); Y.clear(); rho.clear();
      double gn2 = 0;
      for (double v : g) gn2 += v * v;
      const double sc = 0.1 / std::sqrt(gn2 + 1e-30);
      for (size_t t = 0; t < n; ++t) d[t] = g[t] * sc;
      slope = -gn2 * sc;
    }
    double step = 1, fn = 0;
    bool ok = false;
    for (int ls = 0; ls < 30; ++ls) {
      for (size_t t = 0; t < n; ++t) xn[t] = x[t] - step * d[t];
      fn = fg(xn, gn);
      if (fn <= f + 1e-4 * step * slope) { ok = true; break; }
      step *= 0.5;
    }
    if (!ok) break;
    std::vector<double> s(n), y(n);
    double sy = 0;
    for (size_t t = 0; t < n; ++t) s[t] = xn[t] - x[t], y[t] = gn[t] - g[t], sy += s[t] * y[t];
    if (sy > 1e-12) {
      if (int(S.size()) == m) { S.erase(S.begin()); Y.erase(Y.begin()); rho.erase(rho.begin()); }
      S.push_back(std::move(s));
      Y.push_back(std::move(y));
      rho.push_back(1 / sy);
    }
    const double df = f - fn;
    x.swap(xn);
    g.swap(gn);
    f = fn;
    if (df < 1e-12 * (1 + std::fabs(f))) break;
  }
  return f;
}

double radius(int z) {
  // single-bond covalent radii (Cordero et al. 2008; sp3 carbon)
  return z == 1 ? 0.31 : element(z).covalent > 0 ? element(z).covalent : 0.77;
}

struct Model {
  int n = 0;
  std::vector<std::vector<std::pair<int, int>>> adj;   // neighbour, bond order
  std::vector<int> hyb;                                 // 1 sp, 2 sp2, 3 sp3
  std::vector<Pair> pairs;
  std::vector<Chiral> chiral;
  std::vector<Planar> planar;
  std::vector<Planar> flat;                             // across double and aromatic bonds: dihedral 0 or 180
  std::vector<std::vector<int>> topo;                   // topological distance, capped at 4 (4 = ≥ 4)
};

int ring_size(const Model& M, int i, int j, int k) {
  // smallest ring through i-j-k: shortest path i → k avoiding j (breadth first, up to 8 atoms)
  std::vector<int> dist(size_t(M.n), -1);
  std::vector<int> q{i};
  dist[size_t(i)] = 0;
  for (size_t h = 0; h < q.size(); ++h) {
    const int u = q[h];
    if (dist[size_t(u)] >= 6) break;
    for (auto [w, o] : M.adj[size_t(u)]) {
      if (w == j || dist[size_t(w)] >= 0) continue;
      if (u == i && w == k) continue;   // not the direct i-k bond of a 3-ring... handled below
      dist[size_t(w)] = dist[size_t(u)] + 1;
      if (w == k) return dist[size_t(w)] + 2;
      q.push_back(w);
    }
  }
  for (auto [w, o] : M.adj[size_t(i)])
    if (w == k) return 3;
  return 0;
}

double bond_length(const MolGraph& g, const Model& M, int i, int j, int order) {
  double r = radius(g.atoms[size_t(i)].element) + radius(g.atoms[size_t(j)].element);
  if (order == 2) r -= 0.19;
  else if (order == 3) r -= 0.32;
  else if (order == 4) r -= 0.12;
  else if (M.hyb[size_t(i)] < 3 && M.hyb[size_t(j)] < 3 && g.atoms[size_t(i)].element != 1 && g.atoms[size_t(j)].element != 1)
    r -= 0.05;   // conjugated single bond
  return r;
}

double ideal_angle(const MolGraph& g, const Model& M, int i, int j, int k) {
  const int rs = ring_size(M, i, j, k);
  if (rs == 3) return 60;
  if (rs == 4) return 90;
  if (rs == 5) return 108;
  const int h = M.hyb[size_t(j)];
  const size_t nb = M.adj[size_t(j)].size();
  if (h == 1) return 180;
  // exocyclic angles at a small-ring atom
  int small = 0;
  for (auto [a, oa] : M.adj[size_t(j)])
    for (auto [b, ob] : M.adj[size_t(j)])
      if (a < b) {
        const int r = ring_size(M, a, j, b);
        if (r >= 3 && r <= 5 && !(a == i && b == k) && !(a == k && b == i)) small = small ? std::min(small, r) : r;
      }
  if (h == 2) {
    if (small && nb == 3) return (360.0 - (small == 3 ? 60 : small == 4 ? 90 : 108)) / 2;
    return 120;
  }
  if (nb > 4) return 90;
  if (small && rs == 0) return small == 3 ? 117 : small == 4 ? 113 : 111;
  (void)g;
  return 109.47;
}

// The angle i-j-k at a square-planar (Ni, Pd, Pt, Au with four neighbours), five- or six-coordinate centre j, from the
// order the neighbours are written in: square planar 1-3 and 2-4 trans; trigonal bipyramid 1 and 5 axial; octahedron
// 1 and 6 axial, 2-4 and 3-5 trans. 0 for other centres.
double polyhedral_angle(const MolGraph& g, const Model& M, int j, int i, int k) {
  const size_t nb = M.adj[size_t(j)].size();
  const int z = g.atoms[size_t(j)].element;
  const bool square = nb == 4 && (z == 28 || z == 46 || z == 78 || z == 79);
  if (!square && nb != 5 && nb != 6) return 0;
  // written order: the atom's order list (bracket hydrogens, -2, are its added hydrogens in order), else bond order
  std::vector<int> ord;
  for (int w : g.atoms[size_t(j)].order)
    if (w >= 0) ord.push_back(w);
  for (auto [w, o] : M.adj[size_t(j)])
    if (std::find(ord.begin(), ord.end(), w) == ord.end()) ord.push_back(w);
  const auto pi = std::find(ord.begin(), ord.end(), i) - ord.begin(), pk = std::find(ord.begin(), ord.end(), k) - ord.begin();
  const long a = std::min(pi, pk), b = std::max(pi, pk);
  if (square) return b - a == 2 ? 180 : 90;
  if (nb == 5) {
    const bool ax_a = a == 0 || a == 4, ax_b = b == 0 || b == 4;
    return ax_a && ax_b ? 180 : ax_a || ax_b ? 90 : 120;
  }
  if (a == 0 && b == 5) return 180;
  if (a == 0 || b == 5) return 90;
  return b - a == 2 ? 180 : 90;
}

Model make_model(const MolGraph& g) {
  Model M;
  M.n = int(g.atoms.size());
  const size_t n = size_t(M.n);
  M.adj.assign(n, {});
  for (const auto& b : g.bonds) {
    M.adj[size_t(b.a)].push_back({b.b, b.order});
    M.adj[size_t(b.b)].push_back({b.a, b.order});
  }
  // hybridisation
  M.hyb.assign(n, 3);
  for (size_t i = 0; i < n; ++i) {
    const int z = g.atoms[i].element;
    int dbl = 0, tri = 0, aro = 0;
    for (auto [w, o] : M.adj[i]) dbl += o == 2, tri += o == 3, aro += o == 4;
    if (z > 10 && !aro) continue;   // S, P, Si …: tetrahedral (sulfones, phosphates)
    if (tri || dbl >= 2) M.hyb[i] = 1;
    else if (dbl || aro) M.hyb[i] = 2;
  }
  for (size_t i = 0; i < n; ++i) {   // amide, aniline and enamine nitrogen: planar
    if (g.atoms[i].element != 7 || M.hyb[i] != 3 || M.adj[i].size() != 3 || g.atoms[i].charge > 0) continue;
    for (auto [w, o] : M.adj[i])
      if (M.hyb[size_t(w)] == 2) { M.hyb[i] = 2; break; }
  }
  // topological distances up to 3
  M.topo.assign(n, std::vector<int>(n, 4));
  for (size_t s = 0; s < n; ++s) {
    M.topo[s][s] = 0;
    std::vector<int> q{int(s)};
    for (size_t h = 0; h < q.size(); ++h) {
      const int u = q[h];
      if (M.topo[s][size_t(u)] >= 3) continue;
      for (auto [w, o] : M.adj[size_t(u)])
        if (M.topo[s][size_t(w)] == 4) { M.topo[s][size_t(w)] = M.topo[s][size_t(u)] + 1; q.push_back(w); }
    }
  }
  // bounds
  std::map<std::pair<int, int>, std::pair<double, double>> b12, b13, b14;
  std::map<std::pair<int, int>, double> len;
  auto key = [](int a, int b) { return a < b ? std::make_pair(a, b) : std::make_pair(b, a); };
  for (const auto& b : g.bonds) {
    const double r = bond_length(g, M, b.a, b.b, b.order);
    len[key(b.a, b.b)] = r;
    b12[key(b.a, b.b)] = {r * 0.99, r * 1.01};
  }
  auto L = [&](int a, int b) { return len.at(key(a, b)); };
  std::map<std::tuple<int, int, int>, double> ang;   // (i, j, k) i<k, radians
  for (int j = 0; j < M.n; ++j)
    for (auto [i, oi] : M.adj[size_t(j)])
      for (auto [k, ok] : M.adj[size_t(j)]) {
        if (i >= k) continue;
        const double th = ideal_angle(g, M, i, j, k) * kPi / 180;
        ang[{i, j, k}] = th;
        const double a = L(i, j), c = L(j, k);
        const double d = std::sqrt(std::max(0.0, a * a + c * c - 2 * a * c * std::cos(th)));
        auto kk = key(i, k);
        if (b12.count(kk)) continue;   // 3-ring: already a bond
        const double tol = M.hyb[size_t(j)] == 1 ? 0.01 : 0.03;
        auto it = b13.find(kk);
        // square-planar metals, trigonal bipyramids and octahedra: the angle from the neighbours' written order
        // (OpenSMILES @SP1, @TB1, @OH1 defaults)
        const double poly = polyhedral_angle(g, M, j, i, k);
        if (poly > 0) {
          const double dp = std::sqrt(std::max(0.0, a * a + c * c - 2 * a * c * std::cos(poly * kPi / 180)));
          if (it == b13.end()) b13[kk] = {dp - 0.05, dp + 0.05};
          continue;
        }
        if (it == b13.end()) b13[kk] = {d - tol, d + tol};
        else it->second = {std::min(it->second.first, d - tol), std::max(it->second.second, d + tol)};
      }
  auto A = [&](int i, int j, int k) { return i < k ? ang.at({i, j, k}) : ang.at({k, j, i}); };
  // E/Z: the side of each directional substituent, normalised to "seen from the double bond outwards"
  std::map<std::pair<int, int>, int> side;   // (double-bond atom, substituent) → +1 / −1
  for (const auto& b : g.bonds) {
    if (!b.dir) continue;
    side[{b.a, b.b}] = -b.dir;   // written a→b as '/': b is below a, seen from a
    side[{b.b, b.a}] = b.dir;
  }
  auto dist14 = [&](int i, int j, int k, int l, double phi) {
    const double a = L(i, j), b = L(j, k), c = L(k, l), t1 = A(i, j, k), t2 = A(j, k, l);
    const Vec3 pi{a * std::cos(t1), a * std::sin(t1), 0};
    const Vec3 pl{b - c * std::cos(t2), c * std::sin(t2) * std::cos(phi), c * std::sin(t2) * std::sin(phi)};
    return norm(pl - pi);
  };
  for (const auto& b : g.bonds) {
    const int j = b.a, k = b.b;
    for (auto [i, oi] : M.adj[size_t(j)]) {
      if (i == k) continue;
      for (auto [l, ol] : M.adj[size_t(k)]) {
        if (l == j || l == i) continue;
        auto kk = key(i, l);
        if (b12.count(kk) || b13.count(kk)) continue;
        double lo = dist14(i, j, k, l, 0), hi = dist14(i, j, k, l, kPi);
        if (b.order == 2 && ring_size(M, j, k, l) == 0 && ring_size(M, i, j, k) == 0) {
          // a specified double bond fixes cis or trans
          int si = 0, sl = 0;
          if (auto it = side.find({j, i}); it != side.end()) si = it->second;
          if (auto it = side.find({k, l}); it != side.end()) sl = -it->second;
          // a substituent without a mark is opposite to the marked one on the same atom
          if (!si) for (auto [x, ox] : M.adj[size_t(j)]) if (x != i && x != k) if (auto it = side.find({j, x}); it != side.end()) si = -it->second;
          if (!sl) for (auto [x, ox] : M.adj[size_t(k)]) if (x != l && x != j) if (auto it = side.find({k, x}); it != side.end()) sl = it->second;
          if (si && sl) {
            const double d = si == sl ? hi : lo;   // same side marks: trans
            lo = hi = d;
            lo -= 0.05, hi += 0.05;
          }
        }
        auto it = b14.find(kk);
        if (it == b14.end()) b14[kk] = {lo - 0.05, hi + 0.05};
        else it->second = {std::min(it->second.first, lo - 0.05), std::max(it->second.second, hi + 0.05)};
      }
    }
  }
  for (const auto& [k, v] : b12) M.pairs.push_back({k.first, k.second, v.first, v.second, 1.0});
  for (const auto& [k, v] : b13) M.pairs.push_back({k.first, k.second, v.first, v.second, 1.0});
  for (const auto& [k, v] : b14)
    if (!b12.count(k) && !b13.count(k)) M.pairs.push_back({k.first, k.second, v.first, v.second, 0.5});
  for (int i = 0; i < M.n; ++i)
    for (int j = i + 1; j < M.n; ++j) {
      if (M.topo[size_t(i)][size_t(j)] < 4 && (b12.count(key(i, j)) || b13.count(key(i, j)) || b14.count(key(i, j)))) continue;
      if (M.topo[size_t(i)][size_t(j)] < 3) continue;
      const double lb = 0.75 * (element(g.atoms[size_t(i)].element).vdw + element(g.atoms[size_t(j)].element).vdw);
      M.pairs.push_back({i, j, lb, 1e9, 0.3});
    }
  // chirality (four neighbours)
  for (int c = 0; c < M.n; ++c) {
    const MolAtom& a = g.atoms[size_t(c)];
    if (!a.chiral || a.order.size() != 4) continue;
    if (std::any_of(a.order.begin(), a.order.end(), [](int v) { return v < 0; })) continue;
    const double r1 = L(c, a.order[1]), r2 = L(c, a.order[2]), r3 = L(c, a.order[3]);
    M.chiral.push_back({c, a.order[1], a.order[2], a.order[3], a.chiral == 1 ? -1 : 1, 0.45 * r1 * r2 * r3});
  }
  for (int c = 0; c < M.n; ++c)
    if (M.hyb[size_t(c)] == 2 && M.adj[size_t(c)].size() == 3)
      M.planar.push_back({c, M.adj[size_t(c)][0].first, M.adj[size_t(c)][1].first, M.adj[size_t(c)][2].first});
  for (const auto& b : g.bonds) {
    if (b.order != 2 && b.order != 4) continue;
    if (M.hyb[size_t(b.a)] != 2 || M.hyb[size_t(b.b)] != 2) continue;   // cumulenes and S=O, P=O are not planar
    for (auto [i, oi] : M.adj[size_t(b.a)])
      for (auto [l, ol] : M.adj[size_t(b.b)])
        if (i != b.b && l != b.a && i != l) M.flat.push_back({b.a, i, b.b, l});
  }
  return M;
}

// The embedding error and its gradient. squeeze > 0 pulls the fourth coordinate to zero.
double error(const Model& M, const std::vector<double>& x, std::vector<double>& gr, double squeeze, int keep = 3) {
  std::fill(gr.begin(), gr.end(), 0.0);
  double e = 0;
  for (const auto& p : M.pairs) {
    const double* a = &x[size_t(p.i) * D];
    const double* b = &x[size_t(p.j) * D];
    double d[D], d2 = 0;
    for (int t = 0; t < D; ++t) d[t] = a[t] - b[t], d2 += d[t] * d[t];
    double coef = 0;
    if (d2 > p.ub * p.ub) {
      const double u2 = p.ub * p.ub, f = d2 / u2 - 1;
      e += p.w * f * f;
      coef = p.w * 2 * f * 2 / u2;
    } else if (d2 < p.lb * p.lb) {
      const double l2 = p.lb * p.lb, den = l2 + d2, f = 2 * l2 / den - 1;
      e += p.w * f * f;
      coef = p.w * 2 * f * (-2 * l2 / (den * den)) * 2;
    } else
      continue;
    for (int t = 0; t < D; ++t) {
      gr[size_t(p.i) * D + size_t(t)] += coef * d[t];
      gr[size_t(p.j) * D + size_t(t)] -= coef * d[t];
    }
  }
  auto P = [&](int i) { return Vec3{x[size_t(i) * D], x[size_t(i) * D + 1], x[size_t(i) * D + 2]}; };
  auto add = [&](int i, const Vec3& v, double s) { for (int t = 0; t < 3; ++t) gr[size_t(i) * D + size_t(t)] += s * v[t]; };
  for (const auto& c : M.chiral) {
    const Vec3 pc = P(c.c), a = P(c.n1) - pc, b = P(c.n2) - pc, d = P(c.n3) - pc;
    const double V = dot(a, cross(b, d)) * c.sign;
    if (V >= c.vmin) continue;
    const double f = c.vmin - V;
    e += 1.0 * f * f;
    const double s = -2.0 * f * c.sign;   // dE/dV
    const Vec3 ga = cross(b, d), gb = cross(d, a), gd = cross(a, b);
    add(c.n1, ga, s); add(c.n2, gb, s); add(c.n3, gd, s);
    add(c.c, ga + gb + gd, -s);
  }
  for (const auto& p : M.planar) {
    const Vec3 pc = P(p.c), a = P(p.n0) - pc, b = P(p.n1) - pc, d = P(p.n2) - pc;
    const double V = dot(a, cross(b, d));
    e += 0.5 * V * V;
    const double s = 1.0 * V;
    const Vec3 ga = cross(b, d), gb = cross(d, a), gd = cross(a, b);
    add(p.n0, ga, s); add(p.n1, gb, s); add(p.n2, gd, s);
    add(p.c, ga + gb + gd, -s);
  }
  for (const auto& p : M.flat) {
    const Vec3 pj = P(p.c), a = P(p.n0) - pj, b = P(p.n1) - pj, d = P(p.n2) - pj;
    const double V = dot(a, cross(b, d));
    e += 0.5 * V * V;
    const double s = 1.0 * V;
    const Vec3 ga = cross(b, d), gb = cross(d, a), gd = cross(a, b);
    add(p.n0, ga, s); add(p.n1, gb, s); add(p.n2, gd, s);
    add(p.c, ga + gb + gd, -s);
  }
  if (squeeze > 0)
    for (int i = 0; i < M.n; ++i)
      for (int t = keep; t < D; ++t) {
        const double w = x[size_t(i) * D + size_t(t)];
        e += squeeze * w * w;
        gr[size_t(i) * D + size_t(t)] += 2 * squeeze * w;
      }
  return e;
}

}  // namespace

std::vector<Vec3> depict(const MolGraph& g0, uint64_t seed) {
  // the written atoms only (implicit hydrogens are labels, not atoms)
  MolGraph g;
  const int nw = g0.heavy > 0 ? std::min<int>(g0.heavy, int(g0.atoms.size())) : int(g0.atoms.size());
  g.atoms.assign(g0.atoms.begin(), g0.atoms.begin() + nw);
  for (const auto& b : g0.bonds)
    if (b.a < nw && b.b < nw) g.bonds.push_back(b);
  Model M;
  M.n = nw;
  const size_t n = size_t(nw);
  if (!n) return {};
  M.adj.assign(n, {});
  for (const auto& b : g.bonds) M.adj[size_t(b.a)].push_back({b.b, b.order}), M.adj[size_t(b.b)].push_back({b.a, b.order});
  M.topo.assign(n, std::vector<int>(n, 4));
  for (size_t s0 = 0; s0 < n; ++s0) {
    M.topo[s0][s0] = 0;
    std::vector<int> q{int(s0)};
    for (size_t h = 0; h < q.size(); ++h) {
      const int u = q[h];
      if (M.topo[s0][size_t(u)] >= 3) continue;
      for (auto [w, o] : M.adj[size_t(u)])
        if (M.topo[s0][size_t(w)] == 4) { M.topo[s0][size_t(w)] = M.topo[s0][size_t(u)] + 1; q.push_back(w); }
    }
  }
  auto key = [](int a, int b) { return a < b ? std::make_pair(a, b) : std::make_pair(b, a); };
  std::map<std::pair<int, int>, std::pair<double, double>> bnd;
  std::map<std::tuple<int, int, int>, double> ang;
  for (const auto& b : g.bonds) {
    bnd[key(b.a, b.b)] = {0.99, 1.01};
    M.pairs.push_back({b.a, b.b, 0.99, 1.01, 2.0});
  }
  // angles: ring interiors, the rest shared out
  for (int j = 0; j < nw; ++j) {
    const auto& nb = M.adj[size_t(j)];
    const size_t k = nb.size();
    if (k < 2) continue;
    std::vector<std::tuple<int, int, int>> ps;   // (ring size or 0, a, b)
    for (size_t x = 0; x < k; ++x)
      for (size_t y = x + 1; y < k; ++y) {
        const int rs = ring_size(M, nb[x].first, j, nb[y].first);
        ps.push_back({rs > 0 && rs <= 8 ? rs : 0, nb[x].first, nb[y].first});
      }
    auto interior = [](int rs) { return 180.0 * (rs - 2) / rs; };
    auto set = [&](int a, int b, double deg) { ang[{std::min(a, b), j, std::max(a, b)}] = deg; };
    if (k == 2) {
      auto [rs, a, b] = ps[0];
      bool linear = false;
      for (auto [w, o] : nb) linear |= o == 3;
      if (nb[0].second == 2 && nb[1].second == 2) linear = true;
      set(a, b, rs ? interior(rs) : linear ? 180 : 120);
    } else if (k == 3) {
      std::sort(ps.begin(), ps.end(), [](const auto& u, const auto& v) {
        const int ru = std::get<0>(u) ? std::get<0>(u) : 99, rv = std::get<0>(v) ? std::get<0>(v) : 99;
        return ru < rv;
      });
      const int r0 = std::get<0>(ps[0]), r1 = std::get<0>(ps[1]);
      if (r0 && r1) {   // fused: two ring interiors, the third angle is what is left
        const double a0 = interior(r0), a1 = interior(r1);
        set(std::get<1>(ps[0]), std::get<2>(ps[0]), a0);
        set(std::get<1>(ps[1]), std::get<2>(ps[1]), a1);
        set(std::get<1>(ps[2]), std::get<2>(ps[2]), 360 - a0 - a1);
      } else if (r0) {
        const double a0 = interior(r0);
        set(std::get<1>(ps[0]), std::get<2>(ps[0]), a0);
        set(std::get<1>(ps[1]), std::get<2>(ps[1]), (360 - a0) / 2);
        set(std::get<1>(ps[2]), std::get<2>(ps[2]), (360 - a0) / 2);
      } else
        for (auto [rs, a, b] : ps) set(a, b, 120);
    } else {
      for (auto [rs, a, b] : ps) set(a, b, rs ? interior(rs) : -1);   // −1: 90 to 180, sorted out by the contacts
    }
  }
  for (const auto& [t, deg] : ang) {
    const auto [a, j, b] = t;
    auto kk = key(a, b);
    if (bnd.count(kk)) continue;
    if (deg < 0) { M.pairs.push_back({a, b, std::sqrt(2.0), 2.0, 0.5}); continue; }
    const double d = std::sqrt(2 - 2 * std::cos(deg * kPi / 180));
    M.pairs.push_back({a, b, d * 0.99, d * 1.01, 1.0});
    bnd[kk] = {d, d};
  }
  // 1-4: chains zig-zag (trans), double bonds as written, the rest between cis and trans
  std::map<std::pair<int, int>, int> side;
  for (const auto& b : g.bonds)
    if (b.dir) side[{b.a, b.b}] = -b.dir, side[{b.b, b.a}] = b.dir;
  auto A = [&](int a, int j, int b) { auto it = ang.find({std::min(a, b), j, std::max(a, b)}); return it == ang.end() || it->second < 0 ? 120.0 : it->second; };
  std::set<std::pair<int, int>> done14;
  for (const auto& b : g.bonds) {
    const int j = b.a, k = b.b;
    for (auto [i, oi] : M.adj[size_t(j)])
      for (auto [l, ol] : M.adj[size_t(k)]) {
        if (i == k || l == j || i == l) continue;
        auto kk = key(i, l);
        if (bnd.count(kk) || done14.count(kk)) continue;
        done14.insert(kk);
        const double t1 = A(i, j, k) * kPi / 180, t2 = A(j, k, l) * kPi / 180;
        const Vec3 pi{std::cos(t1), std::sin(t1), 0};
        auto at = [&](double phi) { return norm(Vec3{1 - std::cos(t2), std::sin(t2) * std::cos(phi), std::sin(t2) * std::sin(phi)} - pi); };
        const double cis = at(0), trans = at(kPi);
        const bool ring = ring_size(M, i, j, k) || ring_size(M, j, k, l);
        double lo = cis, hi = trans, w = 0.2;
        int si = 0, sl = 0;
        if (auto it = side.find({j, i}); it != side.end()) si = it->second;
        if (auto it = side.find({k, l}); it != side.end()) sl = -it->second;
        if (!si) for (auto [x, ox] : M.adj[size_t(j)]) if (x != i && x != k) if (auto it = side.find({j, x}); it != side.end()) si = -it->second;
        if (!sl) for (auto [x, ox] : M.adj[size_t(k)]) if (x != l && x != j) if (auto it = side.find({k, x}); it != side.end()) sl = it->second;
        if (b.order == 2 && si && sl) { lo = hi = si == sl ? trans : cis; w = 1.0; }
        else if (!ring && M.adj[size_t(j)].size() == 2 && M.adj[size_t(k)].size() == 2) { lo = hi = trans; w = 0.5; }
        else if (ring) continue;
        M.pairs.push_back({i, l, lo * 0.98, hi * 1.02, w});
      }
  }
  for (int i = 0; i < nw; ++i)
    for (int j = i + 1; j < nw; ++j)
      if (!bnd.count(key(i, j)) && !done14.count(key(i, j))) M.pairs.push_back({i, j, 1.0, 1e9, 0.5});
  // several starts; the lowest error wins
  std::mt19937_64 rng(seed * 0x9E3779B97F4A7C15ull + 777);
  const double side_len = 2.0 * std::max(2.0, std::sqrt(double(n)) * 1.2);
  std::uniform_real_distribution<double> U(-side_len / 2, side_len / 2);
  std::vector<double> best;
  double best_e = std::numeric_limits<double>::infinity();
  for (int attempt = 0; attempt < 8; ++attempt) {
    std::vector<double> x(n * D);
    for (double& v : x) v = U(rng);
    auto fg = [&M](double sq) { return [&M, sq](const std::vector<double>& xx, std::vector<double>& gg) { return error(M, xx, gg, sq, 2); }; };
    lbfgs(x, fg(0.0), 1500, 1e-6);
    lbfgs(x, fg(0.1), 1500, 1e-6);
    for (size_t i = 0; i < n; ++i) x[i * D + 2] = x[i * D + 3] = 0;
    const double e = lbfgs(x, fg(50.0), 2000, 1e-7);
    if (e < best_e) best_e = e, best = x;
  }
  std::vector<Vec3> p(n);
  Vec3 c{0, 0, 0};
  for (size_t i = 0; i < n; ++i) p[i] = {best[i * D], best[i * D + 1], 0}, c = c + p[i];
  c = c * (1.0 / double(n));
  // principal axis horizontal
  double sxx = 0, syy = 0, sxy = 0;
  for (auto& v : p) {
    v = v - c;
    sxx += v[0] * v[0], syy += v[1] * v[1], sxy += v[0] * v[1];
  }
  const double th = 0.5 * std::atan2(2 * sxy, sxx - syy);
  const double cs = std::cos(-th), sn = std::sin(-th);
  for (auto& v : p) v = {v[0] * cs - v[1] * sn, v[0] * sn + v[1] * cs, 0};
  return p;
}

std::vector<int> chirality_check(const MolGraph& g, const std::vector<Vec3>& pos) {
  std::vector<int> out;
  for (size_t c = 0; c < g.atoms.size(); ++c) {
    const MolAtom& a = g.atoms[c];
    if (!a.chiral || a.order.size() != 4 || std::any_of(a.order.begin(), a.order.end(), [](int v) { return v < 0; })) continue;
    const Vec3 p = pos[c];
    const double V = dot(pos[size_t(a.order[1])] - p, cross(pos[size_t(a.order[2])] - p, pos[size_t(a.order[3])] - p));
    out.push_back((V < 0) == (a.chiral == 1) ? 1 : -1);
  }
  return out;
}

std::vector<Vec3> embed(const MolGraph& g, const EmbedOptions& o) {
  if (g.atoms.empty()) throw std::runtime_error("nothing to embed");
  const Model M = make_model(g);
  const size_t n = size_t(M.n);
  std::mt19937_64 rng(o.seed * 0x9E3779B97F4A7C15ull + 12345);
  const double side = 2.0 * std::max(3.0, 2.5 * std::cbrt(double(n)));
  std::uniform_real_distribution<double> U(-side / 2, side / 2);
  std::vector<Vec3> best;
  double best_err = std::numeric_limits<double>::infinity();
  for (int attempt = 0; attempt < std::max(1, o.attempts); ++attempt) {
    std::vector<double> x(n * D);
    for (double& v : x) v = U(rng);
    auto fg = [&](double sq) {
      return [&M, sq](const std::vector<double>& xx, std::vector<double>& gg) { return error(M, xx, gg, sq); };
    };
    lbfgs(x, fg(0.0), 2000, 1e-5);          // four dimensions: room to untangle
    lbfgs(x, fg(0.2), 2000, 1e-5);          // squeeze into three
    for (size_t i = 0; i < n; ++i) x[i * D + 3] = 0;
    // last pass in three dimensions: the fourth coordinate is frozen at zero by a stiff squeeze
    const double e = lbfgs(x, fg(50.0), 3000, 1e-6);
    std::vector<Vec3> p(n);
    for (size_t i = 0; i < n; ++i) p[i] = {x[i * D], x[i * D + 1], x[i * D + 2]};
    // acceptance: bonds within 5 %, every specified centre the right way round
    bool ok = true;
    for (const auto& pr : M.pairs) {
      if (M.topo[size_t(pr.i)][size_t(pr.j)] != 1) continue;
      const double d = norm(p[size_t(pr.i)] - p[size_t(pr.j)]), r = 0.5 * (pr.lb + pr.ub);
      if (std::fabs(d - r) > 0.05 * r) { ok = false; break; }
    }
    for (int s : chirality_check(g, p)) ok = ok && s > 0;
    if (ok && e < best_err) {
      best_err = e;
      best = std::move(p);
      break;
    }
  }
  if (best.empty()) throw std::runtime_error("could not embed the molecule with its bond lengths and stereochemistry; check the SMILES");
  // centre at the origin
  Vec3 c{0, 0, 0};
  for (const auto& v : best) c = c + v;
  c = c * (1.0 / double(best.size()));
  for (auto& v : best) v = v - c;
  return best;
}

std::shared_ptr<const ForceField> molecule_forcefield(const MolGraph& g, const std::string& path, const std::string& charges,
                                                      std::vector<std::string>& notes, std::string* name, const std::vector<Vec3>* pos) {
  if (path.empty()) return nullptr;
  if (is_uff(path)) {
    try {
      auto ff = std::make_shared<ForceField>(assign_uff(molecule_system(g, pos ? *pos : std::vector<Vec3>(g.atoms.size()))));
      if (name) *name = "UFF";
      for (const auto& note : ff->notes)
        if (note.find("closest") != std::string::npos || note.find("more than six") != std::string::npos) notes.push_back(note);
      return ff;
    } catch (const std::exception& e) {
      notes.push_back(std::string("no UFF clean-up: ") + e.what());
      return nullptr;
    }
  }
  try {
    FFDef def = load_forcefield(path);
    if (name) *name = def.name;
    System s0 = molecule_system(g, std::vector<Vec3>(g.atoms.size()));
    const TypingResult t = assign_types(s0, def);
    int untyped = 0;
    for (const auto& x : t.types) untyped += x.empty();
    if (untyped) {
      notes.push_back(std::to_string(untyped) + " atoms have no " + def.name + " type; the geometry is the embedding's");
      return nullptr;
    }
    ParamReport pr;
    std::shared_ptr<const ForceField> ff;
    try {
      ff = std::make_shared<ForceField>(parameterize(s0, def, t.types, charges, &pr, true));
    } catch (const std::exception&) {
      if (charges == "types") throw;
      pr = ParamReport{};
      ff = std::make_shared<ForceField>(parameterize(s0, def, t.types, "types", &pr, true));   // Gasteiger covers C, H, N, O only
      notes.push_back("charges from the force field's types (Gasteiger covers C, H, N and O only)");
    }
    if (!pr.missing.empty())
      notes.push_back(std::to_string(pr.missing.size()) + " " + def.name + " terms are missing and left out of the clean-up (first: " + pr.missing.front() + ")");
    return ff;
  } catch (const std::exception& e) {
    notes.push_back(std::string("no force-field clean-up: ") + e.what());
    return nullptr;
  }
}

bool minimise_molecule(const MolGraph& g, const std::shared_ptr<const ForceField>& ff, double ftol, std::vector<Vec3>& pos, double* energy,
                       std::string* why) {
  if (!ff) return false;
  System s = molecule_system(g, pos);
  RelaxOptions ro;
  ro.field = ff;
  ro.pushoff = false;
  ro.ftol = ftol;
  ro.max_iterations = 20000;
  RelaxReport rep;
  try {
    relax(s, ro, &rep);
  } catch (const std::exception& e) {
    if (why) *why = e.what();
    return false;
  }
  std::vector<Vec3> p;
  for (const auto& a : s.atoms) p.push_back(a.pos);
  const auto chk = chirality_check(g, p);
  if (!std::all_of(chk.begin(), chk.end(), [](int v) { return v > 0; })) {
    if (why) *why = "minimisation inverted a centre";
    return false;
  }
  pos = std::move(p);
  if (energy) *energy = rep.final.total();
  return true;
}

BuildResult build_molecule(const std::string& smiles, const BuildOptions& o) {
  BuildResult R;
  R.graph = parse_smiles(smiles);
  add_hydrogens(R.graph);
  R.info = molecule_info(R.graph);
  if (!R.info.problems.empty()) throw std::runtime_error(R.info.problems.front());
  for (const auto& a : R.graph.atoms)
    if (a.chiral && a.order.size() != 4)
      R.notes.push_back("a stereocentre with three neighbours (lone pair) is built without its configuration");

  std::string ffname;
  std::shared_ptr<const ForceField> ff;
  const int nconf = std::max(1, o.conformers);
  for (int k = 0; k < nconf; ++k) {
    EmbedOptions eo;
    eo.seed = o.seed + uint64_t(k) * 7919;
    Conformer c;
    c.pos = embed(R.graph, eo);
    if (k == 0) {   // UFF places five-coordinate axial pairs from the first embedding (all embeddings share them)
      ff = molecule_forcefield(R.graph, o.forcefield, o.charges, R.notes, &ffname, &c.pos);
      if (!ff && !o.forcefield.empty() && !is_uff(o.forcefield)) {   // UFF covers every element the other cannot type
        std::vector<std::string> un;
        ff = molecule_forcefield(R.graph, "uff", o.charges, un, &ffname, &c.pos);
        if (ff) R.notes.push_back("cleaned up with UFF instead");
      }
      R.method = ff ? "CAPS distance-bounds embedding + " + ffname + " minimisation" : "CAPS distance-bounds embedding";
    }
    if (ff) {
      std::string why;
      if (minimise_molecule(R.graph, ff, o.ftol, c.pos, &c.energy, &why)) c.minimised = true;
      else R.notes.push_back("conformer " + std::to_string(k + 1) + (why == "minimisation inverted a centre" ? ": " + why + "; kept the embedding" : " not minimised: " + why));
    }
    R.conformers.push_back(std::move(c));
  }
  std::stable_sort(R.conformers.begin(), R.conformers.end(), [](const Conformer& a, const Conformer& b) {
    if (a.minimised != b.minimised) return a.minimised;
    return a.energy < b.energy;
  });
  // near-identical minima (same energy to 0.01 kcal/mol) are one conformer
  std::vector<Conformer> uniq;
  for (auto& c : R.conformers)
    if (uniq.empty() || !c.minimised || std::fabs(c.energy - uniq.back().energy) > 0.01) uniq.push_back(std::move(c));
  if (uniq.size() < R.conformers.size())
    R.notes.push_back(std::to_string(R.conformers.size() - uniq.size()) + " embeddings reached a minimum already found");
  R.conformers = std::move(uniq);
  R.system = molecule_system(R.graph, R.conformers.front().pos);
  return R;
}

}  // namespace caps
