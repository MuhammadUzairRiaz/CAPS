// Lattice planes (hkl) cut to the cell, their spacing and the view down their normal (miller.hpp).
#include "caps/miller.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace caps {

namespace {

struct Recip {
  Vec3 g;        // h a* + k b* + l c*
  double len;    // |g| = 1 / d
};
Recip reciprocal(const Cell& cell, int h, int k, int l) {
  if (h == 0 && k == 0 && l == 0) throw std::invalid_argument("(0 0 0) is not a plane");
  const double V = dot(cell.a, cross(cell.b, cell.c));
  if (!cell.valid() || std::fabs(V) < 1e-9) throw std::invalid_argument("lattice planes need a periodic cell");
  const Vec3 as = cross(cell.b, cell.c) * (1 / V), bs = cross(cell.c, cell.a) * (1 / V), cs = cross(cell.a, cell.b) * (1 / V);
  const Vec3 g = as * double(h) + bs * double(k) + cs * double(l);
  return {g, norm(g)};
}

}  // namespace

std::string MillerPlanes::label() const {
  return "(" + std::to_string(hkl[0]) + " " + std::to_string(hkl[1]) + " " + std::to_string(hkl[2]) + ")";
}

double d_spacing(const Cell& cell, int h, int k, int l) { return 1.0 / reciprocal(cell, h, k, l).len; }

std::array<double, 2> view_along(const Vec3& n0) {
  const double L = norm(n0);
  if (L < 1e-12) return {0, 0};
  const Vec3 n = n0 * (1 / L);
  const double rho = std::sqrt(n[0] * n[0] + n[2] * n[2]);
  // the view turns by yaw about y, then by pitch about x (render.cpp View::rot); n must come out as +z (toward the eye)
  return {rho > 1e-12 ? std::atan2(-n[0], n[2]) : 0.0, std::atan2(n[1], rho)};
}

MillerPlanes miller_planes(const Cell& cell, int h, int k, int l, double shift, int max_planes) {
  const Recip r = reciprocal(cell, h, k, l);
  MillerPlanes p;
  p.hkl = {h, k, l};
  p.d = 1.0 / r.len;
  p.normal = r.g * (1.0 / r.len);
  const auto va = view_along(p.normal);
  p.yaw = va[0], p.pitch = va[1];
  const int hk[3] = {h, k, l};
  auto value = [&](const Vec3& f) { return h * f[0] + k * f[1] + l * f[2]; };
  double lo = 1e300, hi = -1e300;
  for (int c = 0; c < 8; ++c) {
    const double v = value({double(c & 1), double((c >> 1) & 1), double((c >> 2) & 1)});
    lo = std::min(lo, v), hi = std::max(hi, v);
  }
  const double eps = 1e-9;
  std::vector<double> levels;
  for (long m = long(std::floor(lo - shift)) - 1; m <= long(std::ceil(hi - shift)) + 1; ++m) {
    const double c = double(m) + shift;
    if (c > lo + eps && c < hi - eps) levels.push_back(c);
  }
  if (max_planes > 0 && int(levels.size()) > max_planes) {
    const double mid = (lo + hi) / 2;
    std::stable_sort(levels.begin(), levels.end(), [&](double x, double y) { return std::fabs(x - mid) < std::fabs(y - mid); });
    levels.resize(size_t(max_planes));
    std::sort(levels.begin(), levels.end());
  }
  auto cart = [&](const Vec3& f) { return cell.origin + cell.a * f[0] + cell.b * f[1] + cell.c * f[2]; };
  const Vec3 u0 = std::fabs(p.normal[0]) < 0.9 ? Vec3{1, 0, 0} : Vec3{0, 1, 0};
  const Vec3 eu = cross(p.normal, u0) * (1 / norm(cross(p.normal, u0))), ev = cross(p.normal, eu);
  for (double c : levels) {
    std::vector<Vec3> pts;   // fractional
    auto add = [&](const Vec3& f) {
      for (const auto& q : pts)
        if (std::fabs(q[0] - f[0]) < 1e-9 && std::fabs(q[1] - f[1]) < 1e-9 && std::fabs(q[2] - f[2]) < 1e-9) return;
      pts.push_back(f);
    };
    for (int ax = 0; ax < 3; ++ax)   // the 12 edges of the unit cube along each axis
      for (int e = 0; e < 4; ++e) {
        Vec3 s{0, 0, 0};
        s[(ax + 1) % 3] = double(e & 1), s[(ax + 2) % 3] = double((e >> 1) & 1);
        const double g0 = value(s);
        if (hk[ax] == 0) {
          if (std::fabs(g0 - c) < eps) { add(s); Vec3 t = s; t[ax] = 1; add(t); }
          continue;
        }
        const double t = (c - g0) / hk[ax];
        if (t < -eps || t > 1 + eps) continue;
        Vec3 f = s;
        f[ax] = std::clamp(t, 0.0, 1.0);
        add(f);
      }
    if (pts.size() < 3) continue;
    std::vector<Vec3> poly;
    Vec3 centre{0, 0, 0};
    for (const auto& f : pts) { poly.push_back(cart(f)); centre = centre + poly.back(); }
    centre = centre * (1.0 / double(poly.size()));
    std::sort(poly.begin(), poly.end(), [&](const Vec3& x, const Vec3& y) {
      return std::atan2(dot(x - centre, ev), dot(x - centre, eu)) < std::atan2(dot(y - centre, ev), dot(y - centre, eu));
    });
    p.polygons.push_back(std::move(poly));
    p.offsets.push_back(c * p.d);
  }
  return p;
}

Mesh miller_mesh(const MillerPlanes& p) {
  Mesh m;
  for (const auto& poly : p.polygons) {
    const uint32_t base = uint32_t(m.vertices.size());
    for (const auto& v : poly) m.vertices.push_back(v), m.normals.push_back(p.normal);
    for (size_t k = 1; k + 1 < poly.size(); ++k) m.triangles.push_back({base, base + uint32_t(k), base + uint32_t(k + 1)});
  }
  return m;
}

}  // namespace caps
