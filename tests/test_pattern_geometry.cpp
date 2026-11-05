#include <gtest/gtest.h>

#include <cmath>

#include "caps/pattern_geometry.hpp"
#include "caps/properties.hpp"

using namespace caps;

namespace {
// CH3–S–S–CH3 as heavy atoms: S–S 2.05 Å, C–S 1.82 Å, both C–S–S angles 104°, the C–S–S–C torsion 90°; two frames,
// the second shifted across the cell edge
Trajectory disulfide() {
  const double a = 104 * M_PI / 180, phi = 90 * M_PI / 180;
  const Vec3 s1{0, 0, 0}, s2{2.05, 0, 0};
  const Vec3 c1 = s1 + Vec3{std::cos(a), std::sin(a), 0} * 1.82;
  const Vec3 c2 = s2 + Vec3{-std::cos(a), std::sin(a) * std::cos(phi), std::sin(a) * std::sin(phi)} * 1.82;
  Trajectory t;
  for (auto [z, p] : std::vector<std::pair<int, Vec3>>{{6, c1}, {16, s1}, {16, s2}, {6, c2}}) {
    Atom at;
    at.element = z, at.pos = p + Vec3{1, 1, 1};
    t.topology.atoms.push_back(at);
  }
  t.topology.bonds = {{0, 1, 1}, {1, 2, 1}, {2, 3, 1}};
  t.topology.cell.a = {20, 0, 0}, t.topology.cell.b = {0, 20, 0}, t.topology.cell.c = {0, 0, 20};
  std::vector<Vec3> f0, f1;
  for (const auto& at : t.topology.atoms) {
    f0.push_back(at.pos);
    Vec3 q = at.pos + Vec3{-1.5, 0, 0};
    q[0] -= 20 * std::floor(q[0] / 20);
    f1.push_back(q);
  }
  t.positions = {f0, f1};
  t.cells = {t.topology.cell, t.topology.cell};
  t.timesteps = {0, 1};
  return t;
}
}  // namespace

// Map numbers read from the text, in pattern order
TEST(PatternGeometry, MapNumbers) {
  std::string plain;
  EXPECT_EQ(smarts_map_numbers("C[S:1][S:2]Cl", &plain), (std::vector<int>{0, 1, 2, 0}));
  EXPECT_EQ(plain, "C[S][S]Cl");
  EXPECT_EQ(smarts_map_numbers("[O:2]c1[c:1]cccc1"), (std::vector<int>{2, 0, 1, 0, 0, 0, 0}));
}

// Distance, angle and torsion of the disulfide bridge, each match once, the same in a frame across the cell edge
TEST(PatternGeometry, DisulfideBridge) {
  const Trajectory t = disulfide();
  const auto d = pattern_geometry(t, {0, 1}, "[S:1][S:2]");
  ASSERT_EQ(d.matches.size(), 1u);
  EXPECT_EQ(d.kind, "distance");
  for (double v : d.values) EXPECT_NEAR(v, 2.05, 1e-12);
  const auto a = pattern_geometry(t, {0, 1}, "[C:1][S:2][S:3]");
  ASSERT_EQ(a.matches.size(), 2u);   // C–S–S from either end
  for (double v : a.values) EXPECT_NEAR(v, 104.0, 1e-9);
  const auto w = pattern_geometry(t, {0, 1}, "[C:1][S:2][S:3][C:4]");
  ASSERT_EQ(w.matches.size(), 1u);
  for (double v : w.values) EXPECT_NEAR(std::fabs(v), 90.0, 1e-9);
  EXPECT_NEAR(w.sd, 0.0, 1e-9);
  // through the analysis
  AnalyzeOptions o;
  o.pattern = "[S:1][S:2]";
  const auto p = analyze(t, {"pattern"}, o)[0];
  EXPECT_NEAR(p.value, 2.05, 1e-12);
  EXPECT_EQ(p.unit, "Å");
}

// A centroid point: a six-membered carbon ring of radius 1.4 Å and an O on the radial line of one carbon, 1.36 Å out:
// the ring centre to O is 2.76 Å, one match however the ring is walked
TEST(PatternGeometry, RingCentroid) {
  Trajectory t;
  for (int k = 0; k < 6; ++k) {
    Atom at;
    at.element = 6, at.pos = {10 + 1.4 * std::cos(k * M_PI / 3), 10 + 1.4 * std::sin(k * M_PI / 3), 10};
    t.topology.atoms.push_back(at);
  }
  Atom o;
  o.element = 8, o.pos = {10 + 2.76, 10, 10};
  t.topology.atoms.push_back(o);
  for (uint32_t k = 0; k < 6; ++k) t.topology.bonds.push_back({k, (k + 1) % 6, 1});
  t.topology.bonds.push_back({0, 6, 1});
  t.topology.cell.a = {20, 0, 0}, t.topology.cell.b = {0, 20, 0}, t.topology.cell.c = {0, 0, 20};
  t.positions = {{}};
  for (const auto& a : t.topology.atoms) t.positions[0].push_back(a.pos);
  t.cells = {t.topology.cell};
  t.timesteps = {0};
  const auto g = pattern_geometry(t, {0}, "[O:2]~[#6:1]1~[#6:1]~[#6:1]~[#6:1]~[#6:1]~[#6:1]~1");
  ASSERT_EQ(g.matches.size(), 1u);
  EXPECT_NEAR(g.values[0], 2.76, 1e-12);
}
