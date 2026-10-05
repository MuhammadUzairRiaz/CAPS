#include <algorithm>
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

TEST(Amber, ThioredoxinMatchesOpenMM) {
  const Trajectory t = open_file(kAmber + "trx.inpcrd", kAmber + "trx.prmtop");
  const System s = t.frame(0);
  ASSERT_TRUE(s.forcefield);
  const ForceField& ff = *s.forcefield;
  EXPECT_EQ(s.atoms.size(), 1654u);
  const EnergyTerms e = energy(s, ff, 500);
  EXPECT_NEAR(e.bond, 631.899271, 5e-5);
  EXPECT_NEAR(e.angle, 898.254266, 5e-5);
  EXPECT_NEAR(e.dihedral + e.improper, 566.445419, 5e-5);
  EXPECT_NEAR(e.vdw, -419.538352, 1e-4);
  EXPECT_NEAR(coulomb_limit(s, ff), -1894.052213, 2e-3);
  EXPECT_GT(ff.impropers.size(), 0u);
  // no SCEE / SCNB in this older file: AMBER's defaults
  EXPECT_NEAR(ff.coul14, 1 / 1.2, 1e-12);
  // the bonded terms as the file lists them: 1-4 pairs from the torsions, all excluded from the ordinary pairs
  for (const auto& p : ff.pairs14) {
    const auto& ex = ff.excluded[p[0]];
    EXPECT_TRUE(std::find(ex.begin(), ex.end(), p[1]) != ex.end());
  }
}

// A force field CAPS cannot hold whole is not taken (the zinc finger's 12-6-4 ion terms, C/r⁴: OpenMM puts
// −89.855 kcal/mol in them); its structure still opens
TEST(Amber, TwelveSixFourIonsLeaveTheForceField) {
  const Trajectory t = open_file(kAmber + "znf.rst", kAmber + "znf.prmtop");
  EXPECT_EQ(t.topology.atoms.size(), 551u);
  EXPECT_FALSE(t.topology.forcefield);
  bool said = false;
  for (const auto& n : t.topology.notes) said = said || n.find("12-6-4") != std::string::npos;
  EXPECT_TRUE(said);
  EXPECT_TRUE(t.topology.has_charges);
  EXPECT_GT(t.topology.bonds.size(), 500u);
}

// What the reader leaves out (the force field) or refuses (the file), and the coordinates that do not fit
TEST(Amber, RefusesWhatItCannotRepresent) {
  const std::string dir = ::testing::TempDir();
  std::ifstream in(kAmber + "phenol.prmtop");
  const std::string text((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
  {
    std::ofstream o(dir + "cmap.prmtop");
    o << text << "%FLAG CHARMM_CMAP_COUNT\n%FORMAT(2I8)\n       0       0\n";
  }
  EXPECT_FALSE(read_amber_prmtop(dir + "cmap.prmtop").ff);   // CMAP: the structure without the file's force field
  {
    std::string pert = text;   // a perturbation topology (IFPERT, the 21st pointer) is refused outright
    const auto at = pert.find("%FLAG POINTERS");
    size_t line3 = at;   // the %FLAG and %FORMAT lines, then two lines of ten pointers
    for (int k = 0; k < 4; ++k) line3 = pert.find('\n', line3) + 1;
    pert.replace(line3, 8, "       1");
    std::ofstream o(dir + "pert.prmtop");
    o << pert;
  }
  EXPECT_THROW(read_amber_prmtop(dir + "pert.prmtop"), ReadError);
  {
    std::ofstream o(dir + "short.inpcrd");
    o << "phenol\n    13\n   1.0 2.0 3.0\n";
  }
  EXPECT_THROW(open_file(dir + "short.inpcrd", kAmber + "phenol.prmtop"), ReadError);
  EXPECT_THROW(open_file(dir + "cmap.prmtop"), ReadError);   // no coordinates beside it either
}

// AMBER trajectories and restarts written by other programs, read on the prmtop's atoms: the values MDAnalysis and ParmEd
// read back from the same files (tests/data/amber/NOTICE)
TEST(Amber, NetcdfTrajectoryAndRestart) {
  const Trajectory t = open_file(kAmber + "znf_box.nc", kAmber + "znf.prmtop");
  ASSERT_EQ(t.frames(), 2u);
  const Vec3 x1 = t.positions[1][5];
  EXPECT_NEAR(x1[0], 6.981757164001465, 1e-6);
  EXPECT_NEAR(x1[1], 8.309558868408203, 1e-6);
  EXPECT_NEAR(x1[2], -0.44021475315093994, 1e-6);
  // the box per frame: lengths 60 then 61, 62, 64 Å; angles 80, 85, 75°
  auto deg = [](const Vec3& u, const Vec3& v) { return std::acos(dot(u, v) / (norm(u) * norm(v))) * 180 / M_PI; };
  EXPECT_NEAR(norm(t.cells[0].a), 60, 1e-9);
  EXPECT_NEAR(norm(t.cells[1].a), 61, 1e-9);
  EXPECT_NEAR(norm(t.cells[1].b), 62, 1e-9);
  EXPECT_NEAR(norm(t.cells[1].c), 64, 1e-9);
  EXPECT_NEAR(deg(t.cells[1].b, t.cells[1].c), 80, 1e-9);
  EXPECT_NEAR(deg(t.cells[1].a, t.cells[1].c), 85, 1e-9);
  EXPECT_NEAR(deg(t.cells[1].a, t.cells[1].b), 75, 1e-9);
  // the first frame's velocities, 1.5 −0.25 2.0 Å/ps, as Å/fs
  ASSERT_EQ(t.topology.velocities.size(), t.topology.atoms.size());
  EXPECT_NEAR(t.topology.velocities[5][0], 1.5e-3, 1e-9);
  EXPECT_NEAR(t.topology.velocities[5][1], -0.25e-3, 1e-9);
  // a NetCDF restart (ParmEd): double coordinates, velocities and the box
  const AmberCoordinates c = read_amber_coordinates(kAmber + "znf.ncrst");
  ASSERT_EQ(c.positions.size(), 551u);
  EXPECT_NEAR(c.positions[7][0], 6.1131157875061035, 1e-9);
  EXPECT_NEAR(c.positions[7][2], -0.24049928784370422, 1e-9);
  ASSERT_EQ(c.velocities.size(), 551u);
  EXPECT_NEAR(c.velocities[7][2], 3.999999985080627e-3, 1e-12);
  EXPECT_TRUE(c.has_box);
  EXPECT_NEAR(norm(c.cell.b), 41, 1e-12);
  // opened as a structure: the restart with its topology beside it (another name: given together)
  const Trajectory r = open_file(kAmber + "znf.ncrst", kAmber + "znf.prmtop");
  EXPECT_EQ(r.frames(), 1u);
  EXPECT_NEAR(r.topology.cell.volume(), 40.0 * 41 * 42, 1e-6);
}

TEST(Amber, MdcrdFramesAndBoxLines) {
  const Trajectory t = open_file(kAmber + "znf.mdcrd", kAmber + "znf.prmtop");
  ASSERT_EQ(t.frames(), 4u);
  EXPECT_FALSE(t.cells[3].valid());
  // frame 4, first and last atoms as MDAnalysis reads them (8.3f: to 0.001 Å)
  EXPECT_NEAR(t.positions[3][0][0], 4.121, 1e-6);
  EXPECT_NEAR(t.positions[3][0][1], 8.316, 1e-6);
  EXPECT_NEAR(t.positions[3][0][2], -1.276, 1e-6);
  // a periodic topology: a box line after each frame — AMBER's 3F8.3, or wider columns separated by spaces
  const std::string dir = ::testing::TempDir();
  System top = open_file(kAmber + "phenol.prmtop").frame(0);
  top.cell.a = {30, 0, 0}, top.cell.b = {0, 30, 0}, top.cell.c = {0, 0, 30};
  for (int style = 0; style < 2; ++style) {
    const std::string f = dir + "box" + std::to_string(style) + ".mdcrd";
    {
      std::ofstream o(f);
      o << "phenol frames\n";
      for (int k = 0; k < 2; ++k) {
        char b[16];
        for (size_t i = 0; i < 3 * top.atoms.size(); ++i) {
          std::snprintf(b, sizeof b, "%8.3f", 1.0 + 0.5 * double(i) + k);
          o << b << ((i + 1) % 10 == 0 ? "\n" : "");
        }
        if ((3 * top.atoms.size()) % 10) o << "\n";
        if (style == 0) std::snprintf(b, sizeof b, "%8.3f", 31.0 + k), o << b << b << b << "\n";
        else o << "  31.000   32.000   33.000\n";
      }
    }
    const Trajectory m = read_amber_mdcrd(f, top);
    ASSERT_EQ(m.frames(), 2u);
    EXPECT_NEAR(m.positions[1][12][2], 1.0 + 0.5 * 38 + 1, 1e-9);
    EXPECT_NEAR(norm(m.cells[1].a), style == 0 ? 32.0 : 31.0, 1e-9);
    EXPECT_NEAR(norm(m.cells[1].c), style == 0 ? 32.0 : 33.0, 1e-9);
  }
  // a frame cut short
  {
    std::ofstream o(dir + "short.mdcrd");
    o << "cut\n   1.000   2.000\n";
  }
  EXPECT_THROW(read_amber_mdcrd(dir + "short.mdcrd", top), ReadError);
}

// The writer: a force field out as prmtop + inpcrd and read back gives the same energy, term by term (OpenMM on the
// written files agrees too: bench notes in the commit); a mixing rule other than AMBER's goes out as the pair table
TEST(Amber, WriteReadRoundTrip) {
  const std::string dir = ::testing::TempDir();
  for (const char* name : {"phenol", "trx"}) {
    const std::string crd = std::string(name) == "phenol" ? "phenol.crd" : "trx.inpcrd";
    const System s = open_file(kAmber + crd, kAmber + name + ".prmtop").frame(0);
    ASSERT_TRUE(s.forcefield);
    write_amber(s, *s.forcefield, dir + name + "_out");
    const System r = open_file(dir + name + "_out.inpcrd", dir + name + "_out.prmtop").frame(0);
    ASSERT_TRUE(r.forcefield);
    ASSERT_EQ(r.atoms.size(), s.atoms.size());
    const EnergyTerms a = energy(s, *s.forcefield, 50), b = energy(r, *r.forcefield, 50);
    EXPECT_NEAR(a.bond, b.bond, 1e-6 * std::max(1.0, std::fabs(a.bond))) << name;
    EXPECT_NEAR(a.angle, b.angle, 1e-6 * std::max(1.0, std::fabs(a.angle))) << name;
    EXPECT_NEAR(a.dihedral + a.improper, b.dihedral + b.improper, 1e-6 * std::max(1.0, std::fabs(a.dihedral))) << name;
    EXPECT_NEAR(a.vdw, b.vdw, 1e-6 * std::max(1.0, std::fabs(a.vdw))) << name;
    EXPECT_NEAR(a.coulomb, b.coulomb, 1e-5 * std::max(1.0, std::fabs(a.coulomb))) << name;   // charges through × 18.2223 and E16.8
    EXPECT_EQ(r.forcefield->pairs14.size(), s.forcefield->pairs14.size()) << name;
  }
  // geometric mixing (OPLS): each type pair's own coefficients in the table, read back as overrides
  const System s = open_file(kAmber + "phenol.prmtop").frame(0);
  ForceField g = *s.forcefield;
  g.mixing = "geometric";
  write_amber(s, g, dir + "geo");
  const System r = open_file(dir + "geo.inpcrd", dir + "geo.prmtop").frame(0);
  ASSERT_TRUE(r.forcefield);
  EXPECT_FALSE(r.forcefield->pair_override.empty());
  EXPECT_NEAR(energy(s, g, 50).vdw, energy(r, *r.forcefield, 50).vdw, 1e-6);
  // what AMBER cannot hold is refused with the reason
  ForceField c2 = *s.forcefield;
  c2.pair_form = "lj9-6";
  EXPECT_THROW(write_amber(s, c2, dir + "c2"), FieldError);
  ForceField ub = *s.forcefield;
  ub.urey_bradley.push_back({0, 2, 10.0, 2.4});
  EXPECT_THROW(write_amber(s, ub, dir + "ub"), FieldError);
}
