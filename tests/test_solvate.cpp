#include <gtest/gtest.h>

#include <cmath>

#include "caps/pack.hpp"
#include "caps/peptide.hpp"
#include "caps/solvate.hpp"

using namespace caps;

TEST(Solvate, PlanCountsWaterAndSalt) {
  SolvateOptions o;
  o.edge = 30;
  o.concentration = 0.15;
  const SolvatePlan p = solvate_plan(nullptr, o);
  // 0.997 g/cm³ × 27 000 Å³ → 900 waters; 0.15 mol/L × 27 000 Å³ → 2 NaCl, each ion displacing a water
  EXPECT_EQ(p.cations, 2);
  EXPECT_EQ(p.anions, 2);
  EXPECT_EQ(p.solvent, 900 - 4);
  EXPECT_NEAR(p.density, 1.0, 0.02);
  o.salt = "CaCl2";
  o.concentration = 0.3;
  const SolvatePlan q = solvate_plan(nullptr, o);
  EXPECT_EQ(q.anions, 2 * q.cations);
  o.ion_mode = 3;
  o.cations = 3, o.anions = 6;
  EXPECT_EQ(solvate_plan(nullptr, o).cations, 3);
  o.solvent = "toluene";
  o.ion_mode = 0;
  const SolvatePlan t = solvate_plan(nullptr, o);
  // 0.862 g/cm³, 92.14 g/mol → 152 molecules in 27 000 Å³
  EXPECT_NEAR(t.solvent, 152, 1);
  EXPECT_NEAR(t.solvent_mass, 92.14, 0.05);
}

TEST(Solvate, NeutralisesAChargedSolute) {
  PeptideOptions po;
  po.sequence = "GKKKG";   // +3: NH3+ COO− cancel, three Lys+
  po.cleanup = false;
  const System pep = build_peptide(po);
  SolvateOptions o;
  o.shape = 2;
  o.padding = 6;
  o.ion_mode = 1;
  const SolvatePlan p = solvate_plan(&pep, o);
  EXPECT_NEAR(p.solute_charge, 3, 1e-9);
  EXPECT_EQ(p.cations, 0);
  EXPECT_EQ(p.anions, 3);
  EXPECT_GT(p.solute_volume, 200);
  EXPECT_LT(p.free_volume, p.box_volume);
  SolvateReport rep;
  const System s = solvate(&pep, o, &rep);
  EXPECT_TRUE(rep.pack.success);
  double q = 0;
  int cl = 0;
  for (const auto& a : s.atoms) q += a.charge, cl += a.element == 17;
  EXPECT_NEAR(q, 0, 1e-6);   // TIP4P/2005 water carries its M-site charge on O: each molecule neutral
  EXPECT_EQ(cl, 3);
  EXPECT_GE(rep.pack.dmin, o.tolerance - 1e-6);
  EXPECT_EQ(s.atoms.size(), pep.atoms.size() + size_t(3 * p.solvent) + 3);
}

TEST(Solvate, WaterModelsAndErrors) {
  SolvateOptions o;
  o.water_model = "SPC/E";
  const System w = solvent_molecule(o);
  ASSERT_EQ(w.atoms.size(), 3u);
  EXPECT_NEAR(norm(w.atoms[1].pos - w.atoms[0].pos), 1.0, 1e-9);
  EXPECT_NEAR(w.atoms[0].charge, -0.8476, 1e-9);
  o.water_model = "TIP5P";
  EXPECT_THROW(solvent_molecule(o), std::invalid_argument);
  o.water_model = "TIP3P";
  o.solvent = "mercury";
  EXPECT_THROW(solvate_plan(nullptr, o), std::invalid_argument);
  o.solvent = "water";
  o.shape = 2;
  EXPECT_THROW(solvate_plan(nullptr, o), std::invalid_argument);   // padding needs a solute
}
