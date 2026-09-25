#include <gtest/gtest.h>

#include <cmath>

#include "caps/chipair.hpp"
#include "caps/elements.hpp"
#include "caps/uff.hpp"

using namespace caps;

TEST(ChiPair, ArgonPairsTouchAtTheBondiSumAndPackUnderTwelve) {
  ChiPairOptions o;
  o.a_smiles = "[Ar]", o.b_smiles = "[Ar]";
  o.samples = 2000, o.pack_trials = 200;
  o.temperatures = {300};
  const ChiPairResult r = chi_by_contacts(o);
  // one atom each: every contact is at 2 R_vdW, so every pair energy is the UFF Lennard-Jones energy there
  double x = 0, d = 0;
  ASSERT_TRUE(uff_vdw(18, x, d));
  const double rc = 2 * element(18).vdw, sigma = x / std::pow(2.0, 1.0 / 6), sr = std::pow(sigma / rc, 6);
  const double e = 4 * d * (sr * sr - sr);
  EXPECT_NEAR(r.aa.e_min, e, 1e-9);
  EXPECT_NEAR(r.aa.e_mean, e, 1e-9);
  EXPECT_NEAR(r.aa.e_t[0], e, 1e-9);
  // random sequential packing of equal spheres around one: below the kissing number 12
  EXPECT_GT(r.aa.z, 7.0);
  EXPECT_LT(r.aa.z, 12.0);
  // the two packings of one kind (different random streams) agree
  EXPECT_NEAR(r.aa.z, r.bb.z, 0.5);
}

TEST(ChiPair, SelfMixingControlIsZeroWithinItsErrorAndWaterIsNot) {
  ChiPairOptions o;
  o.a_smiles = "*CC(*)C", o.b_smiles = "*CC(*)C";
  o.samples = 200000, o.pack_trials = 1000;
  const ChiPairResult self = chi_by_contacts(o);
  EXPECT_LT(std::fabs(self.chi_report), 3 * self.chi_report_error + 0.05);
  o.b_smiles = "O";
  const ChiPairResult water = chi_by_contacts(o);
  EXPECT_GT(water.chi_report, 2.0);
  EXPECT_EQ(capped_unit("*CC(*)c1ccccc1"), "[H]CC([H])c1ccccc1");
  EXPECT_EQ(capped_unit("[*]CC([*:2])C[C@H]1CC1"), "[H]CC([H])C[C@H]1CC1");
}
