#include <gtest/gtest.h>

#include <cmath>
#include <fstream>
#include <string>

#include "caps/amber.hpp"
#include "caps/field.hpp"
#include "caps/io.hpp"

using namespace caps;

namespace {

const std::string kAmber = std::string(CAPS_SOURCE_DIR) + "/tests/data/amber/";

EnergyTerms energy(const System& s, const ForceField& ff, double cutoff) {
  EnergyOptions eo;
  eo.cutoff = cutoff;
  eo.dsf_alpha = 0;   // undamped: DSF tends to the plain Coulomb sum as the cut-off grows
  Evaluator ev(ff, eo);
  std::vector<double> x, f;
  for (const auto& a : s.atoms) x.insert(x.end(), a.pos.begin(), a.pos.end());
  return ev.compute(x, s.cell, f);
}

// Coulomb without a cut-off: undamped DSF is E∞ + a/rc + b/rc², solved from three cut-offs
double coulomb_limit(const System& s, const ForceField& ff) {
  const double r1 = 500, r2 = 1000, r3 = 2000;
  const double e1 = energy(s, ff, r1).coulomb, e2 = energy(s, ff, r2).coulomb, e3 = energy(s, ff, r3).coulomb;
  // Lagrange extrapolation in u = 1/rc to u = 0
  const double u1 = 1 / r1, u2 = 1 / r2, u3 = 1 / r3;
  return e1 * (u2 * u3) / ((u1 - u2) * (u1 - u3)) + e2 * (u1 * u3) / ((u2 - u1) * (u2 - u3)) + e3 * (u1 * u2) / ((u3 - u1) * (u3 - u2));
}

}  // namespace

// AMBER prmtop + coordinates (intermol's test systems, MIT; tests/data/amber) against OpenMM 8.5.2 on the same files
// (Reference platform, NoCutoff, no constraints): every term as OpenMM's HarmonicBond, HarmonicAngle, PeriodicTorsion and
// NonbondedForce (LJ alone with the charges zeroed; Coulomb the rest), kcal/mol.
TEST(Amber, PhenolMatchesOpenMM) {
  const Trajectory t = open_file(kAmber + "phenol.prmtop");   // phenol.crd beside it
  const System s = t.frame(0);
  ASSERT_TRUE(s.forcefield);
  const ForceField& ff = *s.forcefield;
  EXPECT_EQ(s.atoms.size(), 13u);
  EXPECT_EQ(ff.type_names.size(), 4u);
  EXPECT_EQ(s.atoms[0].element, 6);
  EXPECT_NEAR(ff.lj14, 0.5, 1e-12);
  EXPECT_NEAR(ff.coul14, 1 / 1.2, 1e-12);
  const EnergyTerms e = energy(s, ff, 500);
  EXPECT_NEAR(e.bond, 0.178425, 2e-6);
  EXPECT_NEAR(e.angle, 0.018066, 2e-6);
  EXPECT_NEAR(e.dihedral + e.improper, 0.000272, 2e-6);
  EXPECT_NEAR(e.vdw, 3.519492, 2e-6);
  EXPECT_NEAR(coulomb_limit(s, ff), -15.577348, 1e-4);
}

TEST(Amber, ZincFingerMatchesOpenMM) {
  const Trajectory t = open_file(kAmber + "znf.rst", kAmber + "znf.prmtop");
  const System s = t.frame(0);
  ASSERT_TRUE(s.forcefield);
  const ForceField& ff = *s.forcefield;
  EXPECT_EQ(s.atoms.size(), 551u);
  const EnergyTerms e = energy(s, ff, 500);
  EXPECT_NEAR(e.bond, 26.394708, 2e-5);
  EXPECT_NEAR(e.angle, 122.824343, 2e-5);
  EXPECT_NEAR(e.dihedral + e.improper, 319.041935, 2e-5);
  EXPECT_NEAR(e.vdw, -48.406750, 1e-4);
  EXPECT_NEAR(coulomb_limit(s, ff), -1995.355369, 1e-3);
  EXPECT_GT(ff.impropers.size(), 0u);
  // the bonded terms as the file lists them: 1-4 pairs from the torsions, all excluded from the ordinary pairs
  for (const auto& p : ff.pairs14) {
    const auto& ex = ff.excluded[p[0]];
    EXPECT_TRUE(std::find(ex.begin(), ex.end(), p[1]) != ex.end());
  }
}

// What the reader refuses rather than reads in part, and the coordinates that do not fit
TEST(Amber, RefusesWhatItCannotRepresent) {
  const std::string dir = ::testing::TempDir();
  std::ifstream in(kAmber + "phenol.prmtop");
  const std::string text((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
  {
    std::ofstream o(dir + "cmap.prmtop");
    o << text << "%FLAG CHARMM_CMAP_COUNT\n%FORMAT(2I8)\n       0       0\n";
  }
  EXPECT_THROW(read_amber_prmtop(dir + "cmap.prmtop"), ReadError);
  {
    std::ofstream o(dir + "short.inpcrd");
    o << "phenol\n    13\n   1.0 2.0 3.0\n";
  }
  EXPECT_THROW(open_file(dir + "short.inpcrd", kAmber + "phenol.prmtop"), ReadError);
  EXPECT_THROW(open_file(dir + "cmap.prmtop"), ReadError);   // no coordinates beside it either
}
