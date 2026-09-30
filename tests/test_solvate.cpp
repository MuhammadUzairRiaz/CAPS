#include <gtest/gtest.h>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <sstream>

#include "caps/pack.hpp"
#include "caps/peptide.hpp"
#include "caps/solvate.hpp"
#include "caps/io.hpp"
#include "caps/relax.hpp"
#include "caps/water.hpp"

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

// Water models: every one neutral; TIP4P/2005 puts M 0.1546 Å from O on the bisector with the model's geometry, a
// three-site model takes it away again; LAMMPS files leave M out with its charge on O and the tip4p style
TEST(Water, ModelsGeometryAndSites) {
  EXPECT_EQ(caps::water_models().size(), 11u);
  for (const auto& m : caps::water_models()) EXPECT_NEAR(2 * m.q_h + m.q_neg, 0.0, 1e-9) << m.name;
  caps::System s;
  {
    const caps::Trajectory w = caps::open_file(std::string(CAPS_SAMPLES) + "/water.pdb");
    for (int k = 0; k < 3; ++k) {
      caps::System one = w.frame(0);
      for (auto& a : one.atoms) a.pos = a.pos + caps::Vec3{4.0 * k, 0, 0}, a.mol = k + 1;
      const uint32_t base = uint32_t(s.atoms.size());
      for (auto& a : one.atoms) s.atoms.push_back(a);
      for (auto b : one.bonds) s.bonds.push_back({b.i + base, b.j + base, b.order});
    }
    s.has_mol = true;
    s.cell.a = {12, 0, 0}, s.cell.b = {0, 12, 0}, s.cell.c = {0, 0, 12};
  }
  const auto& t = caps::water_model("TIP4P/2005");
  ASSERT_EQ(caps::apply_water_model(s, t), 3u);
  ASSERT_EQ(s.atoms.size(), 12u);
  const auto ws = caps::find_waters(s);
  ASSERT_EQ(ws.size(), 3u);
  for (const auto& w : ws) {
    ASSERT_GE(w[3], 0);
    const auto O = s.atoms[size_t(w[0])].pos, H1 = s.atoms[size_t(w[1])].pos, H2 = s.atoms[size_t(w[2])].pos, M = s.atoms[size_t(w[3])].pos;
    EXPECT_NEAR(caps::norm(H1 - O), 0.9572, 1e-9);
    EXPECT_NEAR(std::acos(caps::dot(H1 - O, H2 - O) / (caps::norm(H1 - O) * caps::norm(H2 - O))) * 180 / M_PI, 104.52, 1e-7);
    EXPECT_NEAR(caps::norm(M - O), 0.1546, 1e-9);
    EXPECT_NEAR(s.atoms[size_t(w[3])].charge, -1.1128, 1e-12);
    EXPECT_EQ(s.atoms[size_t(w[0])].charge, 0.0);
  }
  std::vector<uint32_t> all(s.atoms.size());
  for (size_t i = 0; i < all.size(); ++i) all[i] = uint32_t(i);
  const caps::ForceField F = caps::water_forcefield(s, t, all);
  ASSERT_EQ(F.vsites.size(), 3u);
  EXPECT_NEAR(F.vsites[0].w[0] + F.vsites[0].w[1] + F.vsites[0].w[2], 1.0, 1e-12);
  // the LAMMPS files: 9 atoms, O with M's charge, the tip4p styles
  const auto dir = std::filesystem::temp_directory_path() / "caps_water_lmp";
  std::filesystem::create_directories(dir);
  caps::EnergyOptions e;
  e.electrostatics = caps::EnergyOptions::Electrostatics::PME;
  caps::write_lammps_data_ff(s, F, e, (dir / "w.data").string(), false);
  caps::write_lammps_input(s, F, e, "w.data", (dir / "w.in").string(), 0, true, {}, {}, nullptr);
  std::ifstream dat(dir / "w.data"), inp(dir / "w.in");
  std::stringstream ds, is;
  ds << dat.rdbuf(), is << inp.rdbuf();
  EXPECT_NE(ds.str().find("9 atoms"), std::string::npos);
  EXPECT_NE(ds.str().find("-1.112800"), std::string::npos) << "M's charge on O";
  EXPECT_NE(is.str().find("lj/cut/tip4p/long 1 2 1 1 0.154600"), std::string::npos) << is.str();
  EXPECT_NE(is.str().find("pppm/tip4p"), std::string::npos);
  std::filesystem::remove_all(dir);
  // back to three sites: the M sites leave
  ASSERT_EQ(caps::apply_water_model(s, caps::water_model("spce")), 3u);
  EXPECT_EQ(s.atoms.size(), 9u);
  EXPECT_NEAR(s.atoms[0].charge, -0.8476, 1e-12);
}
