// Exact Voronoi / radical tessellation, one convex cell at a time (see caps/voronoi.hpp).
#include "caps/voronoi.hpp"

#include <algorithm>
#include <cmath>
#include <map>

#include "caps/elements.hpp"
#include "parallel.hpp"

namespace caps {
namespace {

// A convex polyhedron around the origin (the atom): vertices, faces as vertex loops, each face's plane distance from
// the origin and the atom behind it (−1: the starting box).
struct Poly {
  std::vector<Vec3> v;
  std::vector<std::vector<int>> f;
  std::vector<int> nb;
  std::vector<double> h;

  explicit Poly(double L) {
    for (int z : {-1, 1})
      for (int k = 0; k < 4; ++k) {
        const int x = (k == 1 || k == 2) ? 1 : -1, y = k >= 2 ? 1 : -1;
        v.push_back({x * L, y * L, z * L});
      }
    f = {{0, 1, 2, 3}, {4, 5, 6, 7}, {0, 1, 5, 4}, {3, 2, 6, 7}, {0, 3, 7, 4}, {1, 2, 6, 5}};
    nb.assign(6, -1);
    h.assign(6, L);
  }

  double radius() const {
    double r = 0;
    for (const auto& loop : f)
      for (int k : loop) r = std::max(r, norm(v[size_t(k)]));
    return r;
  }

  // keeps the part with n·x ≤ d (n a unit vector); false when the plane does not cut
  bool clip(const Vec3& n, double d, int who) {
    const double eps = 1e-10 * std::max(1.0, std::fabs(d));
    std::vector<double> s(v.size(), 0.0);
    bool any = false;
    for (const auto& loop : f)
      for (int k : loop) {
        s[size_t(k)] = dot(v[size_t(k)], n) - d;
        if (s[size_t(k)] > eps) any = true;
      }
    if (!any) return false;
    std::map<std::pair<int, int>, int> cut;
    auto inter = [&](int a, int b) {
      const auto key = std::make_pair(std::min(a, b), std::max(a, b));
      auto it = cut.find(key);
      if (it != cut.end()) return it->second;
      const double t = s[size_t(a)] / (s[size_t(a)] - s[size_t(b)]);
      v.push_back(v[size_t(a)] + (v[size_t(b)] - v[size_t(a)]) * t);
      s.push_back(0.0);
      cut[key] = int(v.size()) - 1;
      return int(v.size()) - 1;
    };
    std::vector<std::vector<int>> nf;
    std::vector<int> nnb;
    std::vector<double> nh;
    std::vector<int> on;
    for (size_t q = 0; q < f.size(); ++q) {
      const auto& loop = f[q];
      std::vector<int> o;
      for (size_t k = 0; k < loop.size(); ++k) {
        const int a = loop[k], b = loop[(k + 1) % loop.size()];
        const double sa = s[size_t(a)], sb = s[size_t(b)];
        if (sa <= eps) o.push_back(a);
        if ((sa < -eps && sb > eps) || (sa > eps && sb < -eps)) o.push_back(inter(a, b));
      }
      if (o.size() < 3) continue;
      for (int k : o)
        if (std::fabs(s[size_t(k)]) <= eps) on.push_back(k);
      nf.push_back(std::move(o));
      nnb.push_back(nb[q]);
      nh.push_back(h[q]);
    }
    std::sort(on.begin(), on.end());
    on.erase(std::unique(on.begin(), on.end()), on.end());
    // the new face: the points on the plane, in order around their centre
    if (on.size() >= 3) {
      Vec3 c{0, 0, 0};
      for (int k : on) c = c + v[size_t(k)];
      c = c * (1.0 / double(on.size()));
      Vec3 e1{0, 0, 0};
      for (int k : on) {
        const Vec3 r = v[size_t(k)] - c;
        if (norm(r) > norm(e1)) e1 = r;
      }
      if (norm(e1) > 0) {
        e1 = e1 * (1 / norm(e1));
        const Vec3 e2 = cross(n, e1);
        std::sort(on.begin(), on.end(), [&](int a, int b) {
          const Vec3 ra = v[size_t(a)] - c, rb = v[size_t(b)] - c;
          return std::atan2(dot(ra, e2), dot(ra, e1)) < std::atan2(dot(rb, e2), dot(rb, e1));
        });
        nf.push_back(on);
        nnb.push_back(who);
        nh.push_back(d);
      }
    }
    f = std::move(nf);
    nb = std::move(nnb);
    h = std::move(nh);
    return true;
  }
};

double loop_area(const std::vector<Vec3>& v, const std::vector<int>& loop) {
  Vec3 a{0, 0, 0};
  for (size_t k = 1; k + 1 < loop.size(); ++k)
    a = a + cross(v[size_t(loop[k])] - v[size_t(loop[0])], v[size_t(loop[k + 1])] - v[size_t(loop[0])]);
  return 0.5 * norm(a);
}

}  // namespace

std::vector<VoronoiCell> voronoi_cells(const System& s, const VoronoiOptions& o) {
  const size_t n = s.atoms.size();
  std::vector<VoronoiCell> out(n);
  if (n == 0) return out;
  const Cell& c = s.cell;
  const bool per = c.valid();
  std::vector<double> w(n, 0.0);
  double wmin = 0;
  if (o.radical) {
    wmin = 1e300;
    for (size_t i = 0; i < n; ++i) w[i] = std::pow(element(s.atoms[i].element).vdw, 2), wmin = std::min(wmin, w[i]);
  }
  // bins in fractional coordinates (periodic) or over the bounding box
  Vec3 lo{0, 0, 0}, span{1, 1, 1};
  double width[3] = {1, 1, 1};
  if (per) {
    const double V = c.volume();
    width[0] = V / norm(cross(c.b, c.c)), width[1] = V / norm(cross(c.c, c.a)), width[2] = V / norm(cross(c.a, c.b));
  } else {
    Vec3 hi{-1e300, -1e300, -1e300};
    lo = {1e300, 1e300, 1e300};
    for (const auto& a : s.atoms)
      for (int d = 0; d < 3; ++d) lo[size_t(d)] = std::min(lo[size_t(d)], a.pos[size_t(d)]), hi[size_t(d)] = std::max(hi[size_t(d)], a.pos[size_t(d)]);
    for (int d = 0; d < 3; ++d) span[size_t(d)] = std::max(1e-6, hi[size_t(d)] - lo[size_t(d)]), width[d] = span[size_t(d)];
  }
  const double dens = per ? double(n) / c.volume() : double(n) / std::max(1.0, width[0] * width[1] * width[2]);
  const double target = std::clamp(1.2 * std::cbrt(1.0 / std::max(dens, 1e-6)), 1.0, 8.0);   // about one atom per bin
  int nb[3];
  for (int d = 0; d < 3; ++d) nb[d] = std::max(1, std::min(200, int(width[d] / target)));
  std::vector<Vec3> frac(n);
  for (size_t i = 0; i < n; ++i) {
    if (per) {
      Vec3 f = c.to_fractional(s.atoms[i].pos);
      for (int d = 0; d < 3; ++d) f[size_t(d)] -= std::floor(f[size_t(d)]);
      frac[i] = f;
    } else {
      for (int d = 0; d < 3; ++d) frac[i][size_t(d)] = std::clamp((s.atoms[i].pos[size_t(d)] - lo[size_t(d)]) / span[size_t(d)], 0.0, 1.0 - 1e-12);
    }
  }
  auto bin = [&](const Vec3& f, int d) { return std::clamp(int(f[size_t(d)] * nb[d]), 0, nb[d] - 1); };
  std::vector<std::vector<uint32_t>> grid(static_cast<size_t>(nb[0]) * size_t(nb[1]) * size_t(nb[2]));
  for (uint32_t i = 0; i < n; ++i) grid[(size_t(bin(frac[i], 0)) * size_t(nb[1]) + size_t(bin(frac[i], 1))) * size_t(nb[2]) + size_t(bin(frac[i], 2))].push_back(i);
  double binw = 1e300;
  for (int d = 0; d < 3; ++d) binw = std::min(binw, width[d] / nb[d]);
  const double L = per ? 2 * std::max({width[0], width[1], width[2]}) + 10 : 1e3;
  const int maxreach = per ? 1000 : std::max({nb[0], nb[1], nb[2]});

  // the cells are independent: split over the threads
  ThreadPool pool(o.threads > 0 ? o.threads : default_threads());
  pool.run(n, [&](int, size_t first, size_t last) {
  for (size_t i = first; i < last; ++i) {
    if (!o.only.empty() && !o.only[i]) continue;
    Poly P(L);
    const int b0[3] = {bin(frac[i], 0), bin(frac[i], 1), bin(frac[i], 2)};
    for (int reach = 0; reach <= maxreach; ++reach) {
      // the shell of bins at Chebyshev distance `reach`, with the periodic image each one stands for
      std::vector<std::pair<double, std::pair<Vec3, uint32_t>>> cand;
      for (int dx = -reach; dx <= reach; ++dx)
        for (int dy = -reach; dy <= reach; ++dy)
          for (int dz = -reach; dz <= reach; ++dz) {
            if (std::max({std::abs(dx), std::abs(dy), std::abs(dz)}) != reach) continue;
            int b[3] = {b0[0] + dx, b0[1] + dy, b0[2] + dz};
            int wrap[3] = {0, 0, 0};
            bool skip = false;
            for (int d = 0; d < 3; ++d) {
              if (per) {
                wrap[d] = int(std::floor(double(b[d]) / nb[d]));
                b[d] -= wrap[d] * nb[d];
              } else if (b[d] < 0 || b[d] >= nb[d]) skip = true;
            }
            if (skip) continue;
            for (uint32_t j : grid[(size_t(b[0]) * size_t(nb[1]) + size_t(b[1])) * size_t(nb[2]) + size_t(b[2])]) {
              Vec3 r;
              if (per) {
                const Vec3 df{frac[j][0] + wrap[0] - frac[i][0], frac[j][1] + wrap[1] - frac[i][1], frac[j][2] + wrap[2] - frac[i][2]};
                r = c.a * df[0] + c.b * df[1] + c.c * df[2];
              } else {
                r = s.atoms[j].pos - s.atoms[i].pos;
              }
              const double r2 = dot(r, r);
              if (r2 < 1e-12) continue;   // the atom itself (or one on top of it)
              cand.push_back({r2, {r, j}});
            }
          }
      std::sort(cand.begin(), cand.begin() + long(cand.size()), [](const auto& a, const auto& b) { return a.first < b.first; });
      for (const auto& [r2, rj] : cand) {
        const double r = std::sqrt(r2);
        const double d = (r2 + w[i] - w[rj.second]) / (2 * r);
        if (d > P.radius() + 1e-9) continue;   // cannot reach the cell
        P.clip(rj.first * (1 / r), d, int(rj.second));
      }
      // no atom beyond the shells searched can cut the cell any more
      const double R = reach * binw;   // every atom closer than this has been seen
      const double reach_plane = R > 0 ? (R * R - (w[i] - wmin)) / (2 * R) : 0;
      if (R > 0 && reach_plane >= P.radius()) break;
    }
    VoronoiCell& cell = out[i];
    cell.index.assign(13, 0);
    for (size_t q = 0; q < P.f.size(); ++q) {
      const double a = loop_area(P.v, P.f[q]);
      cell.volume += a * P.h[q] / 3;
      cell.area += a;
      if (P.nb[q] < 0) { cell.bounded = false; continue; }
      if (a < o.face_area_min || a <= 0) continue;
      int edges = 0;
      const auto& loop = P.f[q];
      for (size_t k = 0; k < loop.size(); ++k)
        if (norm(P.v[size_t(loop[(k + 1) % loop.size()])] - P.v[size_t(loop[k])]) > o.edge_min) ++edges;
      if (edges < 3) continue;
      ++cell.faces;
      cell.max_face_order = std::max(cell.max_face_order, edges);
      ++cell.index[size_t(std::min(edges, 12))];
      cell.neighbours.push_back(P.nb[q]);
    }
    if (!cell.bounded) cell.volume = 0;
  }
  });
  return out;
}

}  // namespace caps
