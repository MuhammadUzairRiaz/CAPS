#include <gtest/gtest.h>

#include <cmath>

#include "caps/conformers.hpp"
#include "caps/molecule.hpp"
#include "caps/torsion.hpp"

using namespace caps;

namespace {
System from_smiles(const char* smi) {
  BuildOptions b;
  b.forcefield = "uff";
  return build_molecule(smi, b).system;
}
}  // namespace

// n-Butane: one rotor; the anti conformer lowest, gauche above it by a fraction of a kcal/mol (the force field's
// anti–gauche gap), and both found from random starts.
TEST(Conformers, ButaneAntiAndGauche) {
  const System s = from_smiles("CCCC");
  const auto rot = rotatable_bonds(s);
  ASSERT_EQ(rot.size(), 1u);
  ConformerSearchOptions o;
  o.trials = 20;
  const auto r = conformer_search(s, o);
  ASSERT_GE(r.conformers.size(), 2u);
  ASSERT_LE(r.conformers.size(), 3u);   // gauche± may superpose within the cut-off (four heavy atoms)
  System lo = s;
  for (size_t i = 0; i < lo.atoms.size(); ++i) lo.atoms[i].pos = r.conformers[0].pos[i];
  const std::array<int, 4> d{int(rot[0][0]), int(rot[0][1]), int(rot[0][2]), int(rot[0][3])};
  EXPECT_GT(std::fabs(dihedral_angle(lo, d)), 150.0);
  EXPECT_GT(r.conformers[1].relative, 0.2);
  EXPECT_LT(r.conformers[1].relative, 2.0);
  double z = 0;
  for (const auto& c : r.conformers) z += c.population;
  EXPECT_NEAR(z, 1.0, 1e-9);
  EXPECT_GT(r.conformers[0].population, r.conformers[1].population);
}

// Benzene: no rotors, one conformer; the anneal search runs and every snapshot is minimised.
TEST(Conformers, RigidAndAnneal) {
  const System s = from_smiles("c1ccccc1");
  ConformerSearchOptions o;
  const auto r = conformer_search(s, o);
  EXPECT_EQ(r.rotors, 0);
  EXPECT_EQ(r.conformers.size(), 1u);
  const System b = from_smiles("CCCC");
  o.method = "anneal";
  o.trials = 6;
  o.anneal_ps = 0.5;
  const auto a = conformer_search(b, o);
  EXPECT_EQ(a.minima.size(), 6u);
  EXPECT_GE(a.conformers.size(), 1u);
}
