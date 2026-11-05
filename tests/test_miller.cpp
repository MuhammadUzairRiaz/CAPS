#include <gtest/gtest.h>

#include <cmath>

#include "caps/miller.hpp"

using namespace caps;

namespace {
Cell cell_of(double a, double b, double c, double al, double be, double ga) {   // degrees, a along x, b in the xy plane
  const double r = M_PI / 180, ca = std::cos(al * r), cb = std::cos(be * r), cg = std::cos(ga * r), sg = std::sin(ga * r);
  Cell x;
  x.a = {a, 0, 0};
  x.b = {b * cg, b * sg, 0};
  const double cx = c * cb, cy = c * (ca - cb * cg) / sg;
  x.c = {cx, cy, std::sqrt(c * c - cx * cx - cy * cy)};
  return x;
}
double area(const std::vector<Vec3>& p) {
  Vec3 s{0, 0, 0};
  for (size_t k = 0; k < p.size(); ++k) s = s + cross(p[k], p[(k + 1) % p.size()]);
  return norm(s) / 2;
}
}  // namespace

// d-spacings against the closed forms: cubic a/√(h²+k²+l²); hexagonal 1/d² = 4(h²+hk+k²)/(3a²) + l²/c²; triclinic
// 1/d² = hᵀ G⁻¹ h with G the metric tensor (International Tables B, ch. 1.1)
TEST(Miller, Spacings) {
  const Cell cub = cell_of(4.05, 4.05, 4.05, 90, 90, 90);
  EXPECT_NEAR(d_spacing(cub, 1, 1, 1), 4.05 / std::sqrt(3.0), 1e-12);
  EXPECT_NEAR(d_spacing(cub, 2, 0, 0), 4.05 / 2, 1e-12);
  EXPECT_NEAR(d_spacing(cub, 2, 2, 2), d_spacing(cub, 1, 1, 1) / 2, 1e-12);
  const Cell hex = cell_of(2.46, 2.46, 6.70, 90, 90, 120);
  for (auto [h, k, l] : {std::array{1, 0, 0}, std::array{1, 1, 0}, std::array{1, 0, 2}, std::array{2, -1, 3}}) {
    const double inv = 4.0 * (h * h + h * k + k * k) / (3 * 2.46 * 2.46) + l * l / (6.70 * 6.70);
    EXPECT_NEAR(d_spacing(hex, h, k, l), 1 / std::sqrt(inv), 1e-12) << h << k << l;
  }
  const Cell tri = cell_of(5.1, 6.3, 7.2, 81, 97, 106);
  const Vec3 A[3] = {tri.a, tri.b, tri.c};
  double G[3][3], Gi[3][3];
  for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) G[i][j] = dot(A[i], A[j]);
  const double det = G[0][0] * (G[1][1] * G[2][2] - G[1][2] * G[2][1]) - G[0][1] * (G[1][0] * G[2][2] - G[1][2] * G[2][0]) + G[0][2] * (G[1][0] * G[2][1] - G[1][1] * G[2][0]);
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) {
      const int i1 = (j + 1) % 3, i2 = (j + 2) % 3, j1 = (i + 1) % 3, j2 = (i + 2) % 3;
      Gi[i][j] = (G[i1][j1] * G[i2][j2] - G[i1][j2] * G[i2][j1]) / det;
    }
  const int hk[3] = {1, -2, 3};
  double inv = 0;
  for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) inv += hk[i] * Gi[i][j] * hk[j];
  EXPECT_NEAR(d_spacing(tri, 1, -2, 3), 1 / std::sqrt(inv), 1e-12);
  EXPECT_THROW(d_spacing(cub, 0, 0, 0), std::invalid_argument);
}

// The (111) planes of a cube of side a: through the origin, two cut the cell, triangles with corners at the cube's
// edge ends (side a√2); shifted by half a spacing, three: a triangle, the regular hexagon of side a/√2 through the centre
// (area 3√3 a²/4) and a triangle. Every corner lies on its plane, the offsets are d apart, the normal is (1,1,1)/√3
TEST(Miller, CubePlanes) {
  const double a = 4.0;
  const Cell cub = cell_of(a, a, a, 90, 90, 90);
  const auto p = miller_planes(cub, 1, 1, 1);
  ASSERT_EQ(p.polygons.size(), 2u);
  for (const auto& poly : p.polygons) {
    ASSERT_EQ(poly.size(), 3u);
    for (size_t k = 0; k < 3; ++k) EXPECT_NEAR(norm(poly[k] - poly[(k + 1) % 3]), a * std::sqrt(2.0), 1e-12);
  }
  EXPECT_NEAR(p.offsets[1] - p.offsets[0], p.d, 1e-12);
  for (int k = 0; k < 3; ++k) EXPECT_NEAR(p.normal[k], 1 / std::sqrt(3.0), 1e-12);
  const auto q = miller_planes(cub, 1, 1, 1, 0.5);
  ASSERT_EQ(q.polygons.size(), 3u);
  const auto& hexagon = q.polygons[1];
  ASSERT_EQ(hexagon.size(), 6u);
  for (size_t k = 0; k < 6; ++k) EXPECT_NEAR(norm(hexagon[k] - hexagon[(k + 1) % 6]), a / std::sqrt(2.0), 1e-12);
  EXPECT_NEAR(area(hexagon), 3 * std::sqrt(3.0) * a * a / 4, 1e-9);
  for (size_t n = 0; n < q.polygons.size(); ++n)
    for (const auto& v : q.polygons[n]) EXPECT_NEAR(dot(v, q.normal), q.offsets[n], 1e-12);
  // (100): planes on the faces only touch the cell; shifted by half, one square through the middle
  const auto f = miller_planes(cub, 1, 0, 0, 0.5);
  ASSERT_EQ(f.polygons.size(), 1u);
  EXPECT_NEAR(area(f.polygons[0]), a * a, 1e-12);
  EXPECT_EQ(miller_planes(cub, 4, 0, 0, 0.5, 2).polygons.size(), 2u);   // at most two, the middle ones
  EXPECT_EQ(miller_mesh(q).triangles.size(), 1u + 4u + 1u);
}

// The camera that looks down a normal turns it onto the eye (+z after yaw about y, then pitch about x), for any
// direction; straight up is the top view (pitch 90°)
TEST(Miller, ViewAlongTheNormal) {
  for (Vec3 n : {Vec3{1, 1, 1}, Vec3{-0.3, 0.2, -0.9}, Vec3{0, 0, 1}, Vec3{1, 0, 0}, Vec3{0.1, -1, 0.05}}) {
    n = n * (1 / norm(n));
    const auto v = view_along(n);
    const double cy = std::cos(v[0]), sy = std::sin(v[0]), cp = std::cos(v[1]), sp = std::sin(v[1]);
    const double x = n[0] * cy + n[2] * sy, z = -n[0] * sy + n[2] * cy, y = n[1];
    EXPECT_NEAR(x, 0, 1e-12);
    EXPECT_NEAR(y * cp - z * sp, 0, 1e-12);
    EXPECT_NEAR(y * sp + z * cp, 1, 1e-12);
  }
  EXPECT_NEAR(view_along({0, 1, 0})[1], M_PI / 2, 1e-12);
}
