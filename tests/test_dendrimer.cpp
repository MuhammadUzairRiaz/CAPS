#include <gtest/gtest.h>

#include <cmath>

#include "caps/dendrimer.hpp"

using namespace caps;

// PAMAM on ethylenediamine: generation 0 is C22H48N10O4 (four amidoamine units, eight NH hydrogens at the ends) and
// each generation doubles the units; every B site filled, so no linear units and both degrees of branching 1
TEST(Dendrimer, PamamFormulas) {
  DendrimerSpec sp;
  sp.preset = "pamam", sp.generation = 0, sp.relax = false;
  DendrimerReport r;
  build_dendrimer(sp, &r);
  EXPECT_EQ(r.formula, "C22H48N10O4");
  EXPECT_EQ(r.units, 4);
  sp.generation = 1;
  build_dendrimer(sp, &r);
  EXPECT_EQ(r.formula, "C62H128N26O12");   // + 8 units C5H10N2O
  EXPECT_EQ(r.units, 12);
  EXPECT_EQ(r.counts.per_shell, (std::vector<int>{4, 8}));
  EXPECT_EQ(r.counts.linear, 0);
  EXPECT_DOUBLE_EQ(r.counts.db_frey, 1.0);
  EXPECT_DOUBLE_EQ(r.counts.db_frechet, 1.0);
}

// PPI on 1,4-diaminobutane: DAB-Am-4 (generation 1) is C16H40N6, DAB-Am-8 C40H96N14; relaxed, no two atoms closer
// than 0.9 Å
TEST(Dendrimer, PpiFormulasAndNoOverlap) {
  DendrimerSpec sp;
  sp.preset = "ppi", sp.generation = 1, sp.relax = false;
  DendrimerReport r;
  build_dendrimer(sp, &r);
  EXPECT_EQ(r.formula, "C16H40N6");
  sp.generation = 2, sp.relax = true;
  const System s = build_dendrimer(sp, &r);
  EXPECT_EQ(r.formula, "C40H96N14");
  double closest = 1e9;
  for (size_t i = 0; i < s.atoms.size(); ++i)
    for (size_t j = i + 1; j < s.atoms.size(); ++j) closest = std::min(closest, norm(s.atoms[i].pos - s.atoms[j].pos));
  EXPECT_GT(closest, 0.9);
  EXPECT_EQ(s.bonds.size(), s.atoms.size() - 1);   // one molecule, no rings
}

// End groups cap every B site left: PAMAM generation 0 with acetyl ends is C38H64N10O12 (8 × C2H2O more)
TEST(Dendrimer, EndGroups) {
  DendrimerSpec sp;
  sp.preset = "pamam", sp.generation = 0, sp.end = "[*:1]C(C)=O", sp.relax = false;
  DendrimerReport r;
  build_dendrimer(sp, &r);
  EXPECT_EQ(r.ends, 8);
  EXPECT_EQ(r.formula, "C38H64N10O12");
}

// Flory's AB₂ statistics: a binary tree has one more terminal unit than dendritic ones (T = D + 1), the tree is drawn
// with exactly the size asked, and Frey's degree of branching tends to 1/2 (the AB₂ bound; units 1/4 T, 1/2 L, 1/4 D)
TEST(Dendrimer, FloryTreeStatistics) {
  for (uint64_t seed : {1u, 2u, 3u}) {
    const auto t = flory_tree(4000, 2, seed);
    ASSERT_EQ(t.parent.size(), 4000u);
    for (size_t k = 1; k < t.parent.size(); ++k) ASSERT_LT(t.parent[k], int(k));
    const auto c = branch_counts(t, 2);
    EXPECT_EQ(c.terminal, c.dendritic + 1);
    EXPECT_NEAR(c.db_frey, 0.5, 0.03) << seed;
    EXPECT_NEAR(double(c.linear) / 4000, 0.5, 0.03);
  }
}

// Slow monomer addition of AB₂: T = D + 1 again, and the degree of branching tends to 2/3 (the rate equations give
// T, L and D a third each)
TEST(Dendrimer, SlowAdditionTendsToTwoThirds) {
  for (uint64_t seed : {1u, 2u, 3u}) {
    const auto t = sma_tree(4000, 2, 0, seed);
    const auto c = branch_counts(t, 2);
    EXPECT_EQ(c.terminal, c.dendritic + 1);
    EXPECT_NEAR(c.db_frey, 2.0 / 3.0, 0.02) << seed;
    EXPECT_NEAR(double(c.terminal) / 4000, 1.0 / 3.0, 0.02);
  }
}

// A hyperbranched bis-MPA molecule of 20 units (no core, the focal acid kept) built in 3D by each growth: 20 bis-MPA
// less one water per ester bond
TEST(Dendrimer, HyperbranchedBuilds) {
  for (const char* g : {"sma", "flory"}) {
    DendrimerSpec sp;
    sp.unit = dendrimer_presets()[2].unit, sp.focal = "[*:2]O", sp.growth = g, sp.units = 20, sp.relax = false;
    DendrimerReport r;
    const System s = build_dendrimer(sp, &r);
    EXPECT_EQ(r.units, 20) << g;
    EXPECT_EQ(r.formula, "C100H162O61") << g;   // 20 bis-MPA (C5H10O4) − 19 H2O
    EXPECT_EQ(s.bonds.size(), s.atoms.size() - 1);
  }
}
