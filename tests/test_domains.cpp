#include <gtest/gtest.h>

#include <cmath>
#include <functional>

#include "caps/domains.hpp"

using namespace caps;

namespace {
// a cubic 30 Å cell on a 0.5 Å grid, the field 1 where `inside` holds (minimum image to the given point) and 0 elsewhere
GridField field_of(const std::function<double(const Vec3&)>& f) {
  GridField g;
  g.cell.a = {30, 0, 0}, g.cell.b = {0, 30, 0}, g.cell.c = {0, 0, 30};
  g.n[0] = g.n[1] = g.n[2] = 60;
  for (int i = 0; i < 60; ++i)
    for (int j = 0; j < 60; ++j)
      for (int k = 0; k < 60; ++k) g.values.push_back(f({(i + 0.5) * 0.5, (j + 0.5) * 0.5, (k + 0.5) * 0.5}));
  return g;
}
double mind(const Vec3& a, const Vec3& b) {
  Vec3 d = a - b;
  for (auto& x : d) x -= 30 * std::round(x / 30);
  return norm(d);
}
}  // namespace

// Two spheres of radius 4 Å, 12 Å apart: two enclosed domains of volume 4/3 π 4³, equivalent radius and ellipsoid axes
// 4 Å (within the voxels' error), centroids at the centres, each other's nearest neighbour at 12 Å
TEST(Domains, TwoSpheres) {
  const Vec3 c1{9, 15, 15}, c2{21, 15, 15};
  const auto g = field_of([&](const Vec3& p) { return mind(p, c1) < 4 || mind(p, c2) < 4 ? 1.0 : 0.0; });
  const auto r = label_domains(g, 0.5, true);
  ASSERT_EQ(r.domains.size(), 2u);
  for (const auto& d : r.domains) {
    EXPECT_TRUE(d.enclosed);
    EXPECT_NEAR(d.volume, 4.0 / 3 * M_PI * 64, 0.03 * 268);
    EXPECT_NEAR(d.equivalent_radius, 4.0, 0.04);
    for (int k = 0; k < 3; ++k) EXPECT_NEAR(d.axes[size_t(k)], 4.0, 0.12);
    EXPECT_NEAR(d.nearest, 12.0, 1e-9);
    EXPECT_LT(std::min(mind(d.centroid, c1), mind(d.centroid, c2)), 1e-9);
  }
}

// A sphere cut by the cell's corner is one domain, enclosed, its centroid at the corner
TEST(Domains, SphereAcrossTheEdges) {
  const Vec3 c{0, 0, 0};
  const auto g = field_of([&](const Vec3& p) { return mind(p, c) < 5 ? 1.0 : 0.0; });
  const auto r = label_domains(g, 0.5, true);
  ASSERT_EQ(r.domains.size(), 1u);
  EXPECT_TRUE(r.domains[0].enclosed);
  EXPECT_LT(mind(r.domains[0].centroid, c), 1e-9);
  EXPECT_NEAR(r.domains[0].equivalent_radius, 5.0, 0.05);
}

// A cylinder along z percolates along c only; a slab normal to z along a and b; the space around the spheres (below
// the threshold) percolates along all three
TEST(Domains, Percolation) {
  const auto cyl = label_domains(field_of([](const Vec3& p) { return std::hypot(p[0] - 15, p[1] - 15) < 3 ? 1.0 : 0.0; }), 0.5, true);
  ASSERT_EQ(cyl.domains.size(), 1u);
  EXPECT_EQ(cyl.domains[0].percolates, (std::array<bool, 3>{false, false, true}));
  EXPECT_FALSE(cyl.domains[0].enclosed);
  const auto slab = label_domains(field_of([](const Vec3& p) { return std::fabs(p[2] - 15) < 2 ? 1.0 : 0.0; }), 0.5, true);
  ASSERT_EQ(slab.domains.size(), 1u);
  EXPECT_EQ(slab.domains[0].percolates, (std::array<bool, 3>{true, true, false}));
  const auto out = label_domains(field_of([](const Vec3& p) { return mind(p, {15, 15, 15}) < 4 ? 1.0 : 0.0; }), 0.5, false);
  ASSERT_EQ(out.domains.size(), 1u);
  EXPECT_EQ(out.domains[0].percolates, (std::array<bool, 3>{true, true, true}));
}

// A slice of a field linear in x and y (inside the cell, away from the wrap) is exact on a tilted plane
TEST(Domains, SliceOfALinearField) {
  const auto g = field_of([](const Vec3& p) { return 2 * p[0] - 0.5 * p[1] + 0.25 * p[2]; });
  const Vec3 u = Vec3{1, 1, 0} * (1 / std::sqrt(2.0)), v{0, 0, 1};
  const auto s = field_slice(g, {15, 15, 15}, u, v, 10, 8, 11, 9);
  for (int y = 0; y < 9; ++y)
    for (int x = 0; x < 11; ++x) {
      const Vec3 p = Vec3{15, 15, 15} + u * (10 * (x / 10.0 - 0.5)) + v * (8 * (y / 8.0 - 0.5));
      EXPECT_NEAR(s[size_t(y * 11 + x)], 2 * p[0] - 0.5 * p[1] + 0.25 * p[2], 1e-9);
    }
}

// The density field holds the group's mass: Σ ρ · voxel = m / (N_A · 10⁻²⁴)
TEST(Domains, FieldHoldsTheMass) {
  System s;
  s.cell.a = {20, 0, 0}, s.cell.b = {0, 18, 0}, s.cell.c = {0, 0, 22};
  for (int k = 0; k < 10; ++k) {
    Atom a;
    a.element = k % 2 ? 6 : 8, a.pos = {1.7 * k, 0.9 * k + 0.3, 19.5 - 0.4 * k};
    s.atoms.push_back(a);
  }
  std::vector<char> pick(10, 0);
  for (int k = 0; k < 10; k += 2) pick[size_t(k)] = 1;   // the oxygens
  const GridField g = density_field(s, pick, 0.7, 1.2);
  double sum = 0;
  for (double v : g.values) sum += v;
  const double vox = s.cell.volume() / double(g.values.size());
  EXPECT_NEAR(sum * vox / 1.66053906660, 5 * 15.999, 1e-9 * 80);
}
