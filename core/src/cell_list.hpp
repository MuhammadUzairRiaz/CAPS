// Internal: periodic cell list over fractional coordinates, shared by bonds and pair analyses.
#pragma once
#include <cstdint>
#include <algorithm>
#include <cmath>
#include <vector>

#include "caps/system.hpp"

namespace caps {

// Cell list over fractional coordinates. For an invalid cell a padded, non-periodic bounding box is used.
struct Grid {
  Cell cell;
  int n[3];
  std::vector<std::vector<uint32_t>> bins;
  std::vector<Vec3> frac;

  Grid(const System& s, double cutoff) {
    cell = s.cell;
    if (!cell.valid()) {
      Vec3 lo{1e300, 1e300, 1e300}, hi{-1e300, -1e300, -1e300};
      for (const auto& a : s.atoms)
        for (int k = 0; k < 3; ++k) { lo[k] = std::min(lo[k], a.pos[k]); hi[k] = std::max(hi[k], a.pos[k]); }
      if (s.atoms.empty()) lo = hi = {0, 0, 0};
      cell.origin = lo - Vec3{1, 1, 1};
      cell.a = {hi[0] - lo[0] + 2, 0, 0};
      cell.b = {0, hi[1] - lo[1] + 2, 0};
      cell.c = {0, 0, hi[2] - lo[2] + 2};
      cell.periodic = {false, false, false};
    }
    const double v = cell.volume();
    const double w[3] = {v / norm(cross(cell.b, cell.c)), v / norm(cross(cell.c, cell.a)), v / norm(cross(cell.a, cell.b))};
    for (int k = 0; k < 3; ++k) n[k] = std::max(1, std::min(512, static_cast<int>(w[k] / cutoff)));
    bins.assign(static_cast<size_t>(n[0]) * n[1] * n[2], {});
    frac.resize(s.atoms.size());
    for (size_t i = 0; i < s.atoms.size(); ++i) {
      Vec3 f = cell.to_fractional(s.atoms[i].pos);
      for (int k = 0; k < 3; ++k) {
        if (cell.periodic[k]) f[k] -= std::floor(f[k]);
        f[k] = std::clamp(f[k], 0.0, 1.0 - 1e-12);
      }
      frac[i] = f;
      bins[index(bin(f, 0), bin(f, 1), bin(f, 2))].push_back(static_cast<uint32_t>(i));
    }
  }
  int bin(const Vec3& f, int k) const { return std::min(n[k] - 1, static_cast<int>(f[k] * n[k])); }
  size_t index(int x, int y, int z) const { return (static_cast<size_t>(x) * n[1] + y) * n[2] + z; }

  // Minimum-image separation from fractional coordinates (exact for cells that are not strongly skewed).
  Vec3 sep(uint32_t i, uint32_t j) const {
    Vec3 d = frac[j] - frac[i];
    for (int k = 0; k < 3; ++k)
      if (cell.periodic[k]) d[k] -= std::round(d[k]);
    return cell.a * d[0] + cell.b * d[1] + cell.c * d[2];
  }

  template <class F>
  void for_neighbour_bins(int x, int y, int z, F&& f) const {
    int xs[3], ys[3], zs[3], cx = 0, cy = 0, cz = 0;
    auto gather = [&](int v, int k, int* out, int& c) {
      for (int d = -1; d <= 1; ++d) {
        int u = v + d;
        if (u < 0 || u >= n[k]) {
          if (!cell.periodic[k]) continue;
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
        for (int c = 0; c < cz; ++c) f(bins[index(xs[a], ys[b], zs[c])]);
  }
};


}  // namespace caps
