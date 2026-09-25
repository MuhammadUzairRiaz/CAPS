#include <gtest/gtest.h>

#include <cmath>

#include "caps/interactions.hpp"
#include "caps/solvate.hpp"

using namespace caps;

namespace {

void water(System& s, const Vec3& o, const Vec3& h1, const Vec3& h2, int64_t mol) {
  const uint32_t base = uint32_t(s.atoms.size());
  for (int k = 0; k < 3; ++k) {
    Atom a;
    a.element = k ? 1 : 8;
    a.mol = mol;
    a.pos = k == 0 ? o : k == 1 ? h1 : h2;
    a.charge = k ? 0.417 : -0.834;
    s.atoms.push_back(a);
  }
  s.bonds.push_back({base, base + 1, 1});
  s.bonds.push_back({base, base + 2, 1});
  s.has_mol = true;
}

}  // namespace

TEST(Interactions, WaterDimerHasOneHydrogenBond) {
  System s;
  // donor at the origin with an H pointing along +x; acceptor 2.9 Å away
  water(s, {0, 0, 0}, {0.9572, 0, 0}, {-0.24, 0.927, 0}, 1);
  water(s, {2.9, 0, 0}, {3.14, 0.927, 0}, {3.14, -0.46, 0.8}, 2);
  const auto R = find_interactions(s);
  EXPECT_EQ(R.n_hbonds, 1u);
  ASSERT_EQ(R.hbonds.size(), 1u);
  EXPECT_EQ(R.hbonds[0].donor, 0u);
  EXPECT_EQ(R.hbonds[0].acceptor, 3u);
  EXPECT_NEAR(R.hbonds[0].distance, 2.9, 1e-9);
  EXPECT_LT(R.hbonds[0].angle, 1e-6);
  EXPECT_EQ(R.n_clashes, 0u);
  EXPECT_NEAR(R.net_charge, 0, 1e-9);
}

TEST(Interactions, OverlappingMoleculesClash) {
  System s;
  water(s, {0, 0, 0}, {0.9572, 0, 0}, {-0.24, 0.927, 0}, 1);
  water(s, {1.2, 0.9, 0.3}, {2.1, 0.9, 0.3}, {0.96, 1.83, 0.3}, 2);   // O···O 1.53 Å
  const auto R = find_interactions(s);
  EXPECT_GE(R.n_clashes, 1u);
  ASSERT_FALSE(R.issues.empty());
  EXPECT_EQ(R.issues[0].level, "error");
  EXPECT_EQ(R.issues[0].fix, "push_apart");
}

TEST(Interactions, PackedWaterBox) {
  SolvateOptions o;
  o.edge = 20;
  o.ion_mode = 0;
  o.water_model = "TIP3P";
  System w = solvate(nullptr, o);
  const auto R = find_interactions(w);
  EXPECT_GT(R.n_hbonds, 10u);      // random orientations still give some
  // packed at the usual 2.0 Å: O···O pairs between 2.0 and 0.75 × 3.04 = 2.28 Å are clashes a relaxation removes
  for (const auto& c : R.clashes) EXPECT_GE(c.distance, 2.0 - 1e-6);
  EXPECT_EQ(R.molecules, int(w.atoms.size() / 3));
  // whole molecules straddle the periodic boundary; wrapping the atoms folds them in
  System wrapped = w;
  for (auto& a : wrapped.atoms) a.pos = wrapped.cell.wrap(a.pos);
  EXPECT_EQ(molecules_across_edge(wrapped), 0);
}
