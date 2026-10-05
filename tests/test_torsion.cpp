#include <gtest/gtest.h>

#include <cmath>

#include "caps/analysis.hpp"
#include "caps/molecule.hpp"
#include "caps/torsion.hpp"
#include "caps/uff.hpp"

using namespace caps;

namespace {

System butane(std::array<int, 4>& at) {
  BuildOptions o;
  o.forcefield = "uff";
  System s = build_molecule("CCCC", o).system;
  at = {0, 1, 2, 3};   // the carbons, as written
  return s;
}

}  // namespace

TEST(Torsion, RigidScanHitsEachAngle) {
  std::array<int, 4> at;
  System s = butane(at);
  const ForceField ff = assign_uff(s);
  TorsionScanOptions o;
  o.atoms = at;
  o.step = 30;
  const auto R = torsion_scan(s, ff, o);
  ASSERT_EQ(R.points.size(), 12u);   // −180 … 150: a full turn, 180 is −180
  for (size_t k = 0; k < R.points.size(); ++k) {
    double d = R.points[k].phi - (-180 + 30.0 * double(k));
    d = std::fmod(d + 540.0, 360.0) - 180.0;
    EXPECT_NEAR(d, 0, 1e-6) << k;
  }
  EXPECT_EQ(R.moving.size(), 7u);   // C3, its two H, C4 and its three H
}

TEST(Torsion, ButaneHasTransAndTwoGauche) {
  std::array<int, 4> at;
  System s = butane(at);
  const ForceField ff = assign_uff(s);
  TorsionScanOptions o;
  o.atoms = at;
  o.step = 10;
  o.relax = true;
  const auto R = torsion_scan(s, ff, o);
  ASSERT_GE(R.conformers.size(), 3u);
  EXPECT_EQ(R.conformers[0].state, "trans");
  EXPECT_NEAR(R.conformers[0].energy, 0, 0.05);   // the parabola may dip a little below the lowest sample
  // gauche+ and gauche− next, equal by symmetry, above trans but below the barrier
  int gp = -1, gm = -1;
  for (size_t k = 0; k < R.conformers.size(); ++k) {
    if (R.conformers[k].state == "gauche+") gp = int(k);
    if (R.conformers[k].state == "gauche−") gm = int(k);
  }
  ASSERT_GE(gp, 0);
  ASSERT_GE(gm, 0);
  EXPECT_NEAR(R.conformers[size_t(gp)].energy, R.conformers[size_t(gm)].energy, 0.05);
  EXPECT_NEAR(R.conformers[size_t(gp)].phi, -R.conformers[size_t(gm)].phi, 3.0);
  EXPECT_GT(R.conformers[size_t(gp)].energy, 0.1);
  EXPECT_GT(R.barrier, R.conformers[size_t(gp)].energy + 1.0);
  EXPECT_NEAR(std::fabs(R.conformers[size_t(gp)].phi), 65, 15);
}

TEST(Torsion, RingBondsAreRefused) {
  BuildOptions o;
  o.forcefield = "uff";
  const System ring = build_molecule("C1CCCCC1", o).system;
  EXPECT_THROW(moving_side(ring, 0, 1), std::invalid_argument);
  EXPECT_EQ(torsion_state(180), "trans");
  EXPECT_EQ(torsion_state(-60), "gauche−");
  EXPECT_EQ(torsion_state(120), "anticlinal+");
  EXPECT_EQ(torsion_state(5), "cis");
}

// Backbone torsions of decane set to all-trans (the zig-zag: every dihedral 180°) and to TG (alternately 180° and 60°).
TEST(Torsion, BackboneTorsionPatterns) {
  BuildOptions b;
  b.forcefield = "uff";
  System s = build_molecule("CCCCCCCCCC", b).system;
  const auto bb = backbones(s);
  ASSERT_EQ(bb.size(), 1u);
  ASSERT_EQ(bb[0].size(), 10u);
  std::vector<std::string> notes;
  EXPECT_EQ(set_backbone_torsions(s, {180}, {}, &notes), 7);
  for (size_t k = 0; k + 3 < bb[0].size(); ++k)
    EXPECT_NEAR(std::fabs(dihedral_angle(s, {int(bb[0][k]), int(bb[0][k + 1]), int(bb[0][k + 2]), int(bb[0][k + 3])})), 180.0, 1e-6);
  set_backbone_torsions(s, {180, 60});
  for (size_t k = 0; k + 3 < bb[0].size(); ++k) {
    const double phi = dihedral_angle(s, {int(bb[0][k]), int(bb[0][k + 1]), int(bb[0][k + 2]), int(bb[0][k + 3])});
    EXPECT_NEAR(k % 2 == 0 ? std::fabs(phi) : phi, k % 2 == 0 ? 180.0 : 60.0, 1e-6);
  }
}
