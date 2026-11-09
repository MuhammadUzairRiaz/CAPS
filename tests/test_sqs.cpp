#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <numeric>

#include "caps/crystal.hpp"
#include "caps/spacegroup.hpp"
#include "caps/sqs.hpp"

using namespace caps;

namespace {
System copper() { return read_cif(std::string(CAPS_SOURCE_DIR) + "/data/crystals/copper.cif"); }
}  // namespace

// The four sites of the conventional fcc cell: one substituent is one orbit of 4; two are one orbit of 6 (every pair is a
// nearest-neighbour pair, the L1_0 layering)
TEST(Sqs, FccConventionalOrbits) {
  EnumOptions o;
  o.from = 29, o.to = 79, o.count = 1;
  auto r = enumerate_occupations(copper(), o);
  ASSERT_EQ(r.configurations.size(), 1u);
  EXPECT_EQ(r.multiplicity[0], 4u);
  o.count = 2;
  r = enumerate_occupations(copper(), o);
  ASSERT_EQ(r.configurations.size(), 1u);
  EXPECT_EQ(r.multiplicity[0], 6u);
}

// The 2 × 2 × 2 fcc supercell (32 sites): a pair of substituents falls in the five classes of separation the cell can hold,
// (½ ½ 0), (1 0 0), (1 ½ ½), (1 1 0), (1 1 1) a, of 12, 3, 12, 3, 1 partners, i.e. orbits of 192, 48, 192, 48, 16 = C(32, 2)
TEST(Sqs, FccSupercellPairOrbits) {
  const System s = supercell(copper(), 2, 2, 2);
  EnumOptions o;
  o.from = 29, o.to = 79, o.count = 2;
  const auto r = enumerate_occupations(s, o);
  EXPECT_EQ(r.total, 496u);
  ASSERT_EQ(r.configurations.size(), 5u);
  std::vector<size_t> m = r.multiplicity;
  std::sort(m.begin(), m.end());
  EXPECT_EQ(m, (std::vector<size_t>{16, 48, 48, 192, 192}));
  EXPECT_EQ(std::accumulate(m.begin(), m.end(), size_t(0)), r.total);
  EXPECT_EQ(r.operations, 48 * 32);   // the point group times the 32 fcc translations the cell holds
  // three substituents: the orbits still add up to C(32, 3)
  o.count = 3;
  const auto r3 = enumerate_occupations(s, o);
  EXPECT_EQ(std::accumulate(r3.multiplicity.begin(), r3.multiplicity.end(), size_t(0)), 4960u);
}

// SQS of a 50 % alloy on the 32-site cell: the composition is exact and the pair correlations of the first shells are
// those of the random alloy, (1 − 2x)² = 0, recomputed here from the returned structure
TEST(Sqs, FccHalfAlloy) {
  const System s = supercell(copper(), 2, 2, 2);
  SqsOptions o;
  o.from = 29, o.to = 79, o.fraction = 0.5, o.shells = 2;
  const auto r = special_quasirandom(s, o);
  EXPECT_EQ(r.sites, 32);
  EXPECT_EQ(r.substituted, 16);
  int au = 0;
  for (const auto& a : r.structure.atoms) au += a.element == 79;
  EXPECT_EQ(au, 16);
  ASSERT_EQ(r.correlation.size(), 2u);
  EXPECT_NEAR(r.shell_distance[0], 3.6149 / std::sqrt(2.0), 1e-6);
  EXPECT_NEAR(r.shell_distance[1], 3.6149, 1e-6);
  // independent recount of Π_1
  double sum = 0;
  int pairs = 0;
  for (size_t i = 0; i < s.atoms.size(); ++i)
    for (size_t j = i + 1; j < s.atoms.size(); ++j) {
      if (std::fabs(norm(s.cell.minimum_image(s.atoms[j].pos - s.atoms[i].pos)) - r.shell_distance[0]) > 1e-3) continue;
      const int si = r.structure.atoms[i].element == 29 ? 1 : -1, sj = r.structure.atoms[j].element == 29 ? 1 : -1;
      sum += si * sj, ++pairs;
    }
  EXPECT_EQ(pairs, 32 * 12 / 2);
  EXPECT_NEAR(sum / pairs, r.correlation[0], 1e-12);
  std::printf("Pi = %.4f %.4f, objective %.4f\n", r.correlation[0], r.correlation[1], r.objective);
  EXPECT_NEAR(r.correlation[0], 0.0, 1e-12);
  EXPECT_NEAR(r.correlation[1], 0.0, 1e-12);
}
