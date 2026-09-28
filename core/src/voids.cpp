#include "caps/voids.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <stdexcept>

#include "caps/elements.hpp"

namespace caps {

namespace {
double len(const Vec3& v) { return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]); }
Vec3 cross3(const Vec3& a, const Vec3& b) { return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]}; }
}  // namespace

VoidReport largest_voids(const System& s, const VoidOptions& o) {
  if (!s.cell.valid()) throw std::invalid_argument("voids need a periodic cell");
  VoidReport rep;
  const Cell& c = s.cell;
  const Vec3 L{len(c.a), len(c.b), len(c.c)};
  int n[3];
  for (int k = 0; k < 3; ++k) n[k] = std::max(4, int(std::ceil(L[k] / std::max(0.1, o.grid))));
  for (int k = 0; k < 3; ++k) rep.grid[k] = n[k];
  const size_t N = size_t(n[0]) * n[1] * n[2];
  auto index = [&](int i, int j, int k) { return (size_t(i) * n[1] + j) * n[2] + k; };
  auto point = [&](int i, int j, int k) { return c.to_cartesian({(i + 0.5) / n[0], (j + 0.5) / n[1], (k + 0.5) / n[2]}); };
  std::vector<float> field(N, float(o.reach));
  // how many grid steps along each axis cover a distance r: perpendicular widths of the cell
  const double V = c.volume();
  const Vec3 w{V / len(cross3(c.b, c.c)), V / len(cross3(c.c, c.a)), V / len(cross3(c.a, c.b))};
  for (size_t ai = 0; ai < s.atoms.size(); ++ai) {
    const auto& a = s.atoms[ai];
    const double R = ai < o.radii.size() ? o.radii[ai] : element(a.element).vdw > 0 ? element(a.element).vdw : 1.7;
    const double reach = R + o.reach;
    const Vec3 f = c.to_fractional(a.pos);
    int lo[3], hi[3];
    for (int k = 0; k < 3; ++k) {
      const double centre = (f[k] - std::floor(f[k])) * n[k] - 0.5;
      const int span = int(std::ceil(reach / w[k] * n[k])) + 1;
      lo[k] = int(std::floor(centre)) - span;
      hi[k] = int(std::ceil(centre)) + span;
      if (hi[k] - lo[k] + 1 > n[k]) { lo[k] = 0; hi[k] = n[k] - 1; }
    }
    for (int i = lo[0]; i <= hi[0]; ++i)
      for (int j = lo[1]; j <= hi[1]; ++j)
        for (int k = lo[2]; k <= hi[2]; ++k) {
          const int ii = ((i % n[0]) + n[0]) % n[0], jj = ((j % n[1]) + n[1]) % n[1], kk = ((k % n[2]) + n[2]) % n[2];
          const Vec3 d = c.minimum_image(point(ii, jj, kk) - a.pos);
          const double dist = len(d) - R;
          float& fv = field[index(ii, jj, kk)];
          if (dist < fv) fv = float(dist);
        }
  }
  size_t free0 = 0, freep = 0;
  for (float v : field) { free0 += v > 0; freep += v > o.probe; }
  rep.accessible_point = double(free0) / N;
  rep.accessible_probe = double(freep) / N;
  // local maxima of the field (26 neighbours, periodic), biggest first, kept when they do not overlap a bigger one
  std::vector<std::pair<float, size_t>> maxima;
  for (int i = 0; i < n[0]; ++i)
    for (int j = 0; j < n[1]; ++j)
      for (int k = 0; k < n[2]; ++k) {
        const float v = field[index(i, j, k)];
        if (v < o.min_radius) continue;
        bool top = true;
        for (int di = -1; di <= 1 && top; ++di)
          for (int dj = -1; dj <= 1 && top; ++dj)
            for (int dk = -1; dk <= 1 && top; ++dk) {
              if (!di && !dj && !dk) continue;
              if (field[index((i + di + n[0]) % n[0], (j + dj + n[1]) % n[1], (k + dk + n[2]) % n[2])] > v) top = false;
            }
        if (top) maxima.push_back({v, index(i, j, k)});
      }
  std::sort(maxima.begin(), maxima.end(), [](const auto& x, const auto& y) { return x.first > y.first; });
  for (const auto& [v, q] : maxima) {
    if (int(rep.spheres.size()) >= o.max_count) break;
    const int i = int(q / (size_t(n[1]) * n[2])), j = int((q / n[2]) % n[1]), k = int(q % n[2]);
    const Vec3 p = point(i, j, k);
    bool clear = true;
    for (const auto& sp : rep.spheres)
      if (len(c.minimum_image(p - sp.centre)) < sp.radius + v) { clear = false; break; }
    if (clear) rep.spheres.push_back({p, double(v)});
  }
  rep.largest = rep.spheres.empty() ? 0 : rep.spheres.front().radius;
  return rep;
}

Mesh void_mesh(const std::vector<VoidSphere>& v, int subdivisions) {
  // unit icosahedron, subdivided
  const double t = (1 + std::sqrt(5.0)) / 2;
  std::vector<Vec3> P = {{-1, t, 0}, {1, t, 0}, {-1, -t, 0}, {1, -t, 0}, {0, -1, t}, {0, 1, t}, {0, -1, -t}, {0, 1, -t},
                         {t, 0, -1}, {t, 0, 1}, {-t, 0, -1}, {-t, 0, 1}};
  std::vector<std::array<uint32_t, 3>> T = {{0, 11, 5}, {0, 5, 1}, {0, 1, 7}, {0, 7, 10}, {0, 10, 11}, {1, 5, 9}, {5, 11, 4}, {11, 10, 2}, {10, 7, 6}, {7, 1, 8},
                                            {3, 9, 4}, {3, 4, 2}, {3, 2, 6}, {3, 6, 8}, {3, 8, 9}, {4, 9, 5}, {2, 4, 11}, {6, 2, 10}, {8, 6, 7}, {9, 8, 1}};
  for (auto& p : P) { const double l = len(p); p = {p[0] / l, p[1] / l, p[2] / l}; }
  for (int s = 0; s < subdivisions; ++s) {
    std::map<std::pair<uint32_t, uint32_t>, uint32_t> mid;
    auto m = [&](uint32_t a, uint32_t b) {
      const auto key = std::minmax(a, b);
      if (auto it = mid.find(key); it != mid.end()) return it->second;
      Vec3 q = (P[a] + P[b]) * 0.5;
      const double l = len(q);
      P.push_back({q[0] / l, q[1] / l, q[2] / l});
      return mid[key] = uint32_t(P.size() - 1);
    };
    std::vector<std::array<uint32_t, 3>> T2;
    for (const auto& f : T) {
      const uint32_t a = m(f[0], f[1]), b = m(f[1], f[2]), cc = m(f[2], f[0]);
      T2.push_back({f[0], a, cc});
      T2.push_back({f[1], b, a});
      T2.push_back({f[2], cc, b});
      T2.push_back({a, b, cc});
    }
    T = std::move(T2);
  }
  Mesh out;
  for (const auto& sp : v) {
    const uint32_t base = uint32_t(out.vertices.size());
    for (const auto& p : P) {
      out.vertices.push_back(sp.centre + p * sp.radius);
      out.normals.push_back(p);
    }
    for (const auto& f : T) out.triangles.push_back({base + f[0], base + f[1], base + f[2]});
  }
  return out;
}

void write_voids_pdb(const System& s, const std::vector<VoidSphere>& v, const std::string& path) {
  std::ofstream f(path);
  if (!f) throw std::runtime_error("cannot write " + path);
  char b[128];
  if (s.cell.valid()) {
    const auto ang = [](const Vec3& x, const Vec3& y) { return std::acos(std::clamp(dot(x, y) / (len(x) * len(y)), -1.0, 1.0)) * 180 / 3.14159265358979; };
    std::snprintf(b, sizeof b, "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1           1\n", len(s.cell.a), len(s.cell.b), len(s.cell.c),
                  ang(s.cell.b, s.cell.c), ang(s.cell.a, s.cell.c), ang(s.cell.a, s.cell.b));
    f << b;
  }
  f << "REMARK   1 CAPS voids: largest empty spheres, radius (Å) in the B-factor column\n";
  for (size_t k = 0; k < v.size(); ++k) {
    std::snprintf(b, sizeof b, "HETATM%5zu  VO  VOI X%4zu    %8.3f%8.3f%8.3f  1.00%6.2f           X\n", k + 1, k + 1, v[k].centre[0], v[k].centre[1],
                  v[k].centre[2], v[k].radius);
    f << b;
  }
  f << "END\n";
}

}  // namespace caps
