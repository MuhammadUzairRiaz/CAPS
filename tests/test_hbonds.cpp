#include <gtest/gtest.h>

#include <cmath>

#include "caps/properties.hpp"

using namespace caps;

namespace {
Trajectory frame_of(const System& s) {
  Trajectory t;
  t.topology = s;
  t.positions.push_back({});
  for (const auto& a : s.atoms) t.positions.back().push_back(a.pos);
  t.timesteps.push_back(0);
  return t;
}
void add(System& s, int z, Vec3 p, int64_t mol) {
  Atom a;
  a.element = z, a.pos = p, a.mol = mol;
  s.atoms.push_back(a);
}
// a water donating along x from the origin (O at 0, H at 0.96 Å) and an acceptor O placed from that H at 1.84 Å with
// the D–H···A angle `angle` (degrees)
System water_pair(double angle) {
  System s;
  add(s, 8, {0, 0, 0}, 1);
  add(s, 1, {0.96, 0, 0}, 1);
  add(s, 1, {-0.24, 0.93, 0}, 1);
  const double th = (180 - angle) * M_PI / 180;
  const Vec3 A = Vec3{0.96, 0, 0} + Vec3{std::cos(th), std::sin(th), 0} * 1.84;
  add(s, 8, A, 2);
  add(s, 1, A + Vec3{0.3, 0.9, 0.2}, 2);
  add(s, 1, A + Vec3{0.3, -0.9, -0.2}, 2);
  s.has_mol = true;
  s.bonds = {{0, 1}, {0, 2}, {3, 4}, {3, 5}};
  return s;
}
double hbonds(const System& s, const AnalyzeOptions& o) { return analyze(frame_of(s), {"hbonds"}, o)[0].value; }
}  // namespace

// The IUPAC form (Arunan et al. 2011): H···A below the sum of the van der Waals radii (1.84 Å < 1.20 + 1.52 Å) and the
// D–H···A angle at least 110°: 120° counts, 100° does not; an explicit 1.8 Å limit excludes the 1.84 Å contact
TEST(HBonds, IupacCriterion) {
  AnalyzeOptions o;
  o.hb_criterion = "ha";
  EXPECT_NEAR(hbonds(water_pair(180), o), 1.0, 1e-12);
  EXPECT_NEAR(hbonds(water_pair(120), o), 1.0, 1e-12);
  EXPECT_NEAR(hbonds(water_pair(100), o), 0.0, 1e-12);
  o.hb_ha = 1.8;
  EXPECT_NEAR(hbonds(water_pair(180), o), 0.0, 1e-12);
  // the D···A form keeps its own limits: at 160° at H the H–D···A angle is 13.2° (inside 30°), at 120° it is 40.3°
  AnalyzeOptions lc;
  EXPECT_NEAR(hbonds(water_pair(160), lc), 1.0, 1e-12);
  EXPECT_NEAR(hbonds(water_pair(120), lc), 0.0, 1e-12);
}

// Weak C–H···O: a methyl hydrogen pointing at a water oxygen 2.4 Å away counts only when carbon is a donor
TEST(HBonds, WeakCHDonors) {
  System s;
  add(s, 6, {0, 0, 0}, 1);
  add(s, 1, {1.09, 0, 0}, 1);
  add(s, 1, {-0.36, 1.03, 0}, 1);
  add(s, 1, {-0.36, -0.51, 0.89}, 1);
  add(s, 1, {-0.36, -0.51, -0.89}, 1);
  add(s, 8, {3.49, 0, 0}, 2);
  add(s, 1, {3.8, 0.9, 0}, 2);
  add(s, 1, {3.8, -0.9, 0}, 2);
  s.has_mol = true;
  s.bonds = {{0, 1}, {0, 2}, {0, 3}, {0, 4}, {5, 6}, {5, 7}};
  AnalyzeOptions o;
  o.hb_criterion = "ha";
  EXPECT_NEAR(hbonds(s, o), 0.0, 1e-12);
  o.hb_donors = "C,N,O";
  EXPECT_NEAR(hbonds(s, o), 1.0, 1e-12);
}

// Between two groups: water 1 donates to water 2 and water 3 donates to water 2; counted between molecules 1 and 2
// only, one bond; with the groups swapped, the same
TEST(HBonds, BetweenTwoGroups) {
  System s = water_pair(180);
  // a third water donating to molecule 2's oxygen from the other side
  const Vec3 A = s.atoms[3].pos;
  add(s, 8, A + Vec3{2.8, 0, 0}, 3);
  add(s, 1, A + Vec3{1.84, 0, 0}, 3);
  add(s, 1, A + Vec3{3.04, 0.93, 0}, 3);
  s.bonds.push_back({6, 7});
  s.bonds.push_back({6, 8});
  AnalyzeOptions o;
  EXPECT_NEAR(hbonds(s, o), 2.0, 1e-12);
  o.hb_group_a = "molecule 1", o.hb_group_b = "molecule 2";
  EXPECT_NEAR(hbonds(s, o), 1.0, 1e-12);
  std::swap(o.hb_group_a, o.hb_group_b);
  EXPECT_NEAR(hbonds(s, o), 1.0, 1e-12);
}
