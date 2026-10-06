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
  o.water_model = "TIP5P";   // the lone pairs' charge on O until the model is applied
  const System w5 = solvent_molecule(o);
  EXPECT_NEAR(w5.atoms[0].charge, -0.482, 1e-12);
  EXPECT_NEAR(w5.atoms[1].charge, 0.241, 1e-12);
  o.water_model = "nonsense";
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
  EXPECT_EQ(caps::water_models().size(), 12u);
  for (const auto& m : caps::water_models()) EXPECT_NEAR(2 * m.q_h + (m.sites == 5 ? 2 : 1) * m.q_neg, 0.0, 1e-9) << m.name;
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
  dat.close(); inp.close();   // Windows will not remove a file still open
  std::filesystem::remove_all(dir);
  // back to three sites: the M sites leave
  ASSERT_EQ(caps::apply_water_model(s, caps::water_model("spce")), 3u);
  EXPECT_EQ(s.atoms.size(), 9u);
  EXPECT_NEAR(s.atoms[0].charge, -0.8476, 1e-12);
}

// TIP5P: two lone pairs 0.70 Å from O at 109.47°, symmetric about the molecular plane and away from the hydrogens (GROMACS's
// 3out constants a = −0.344908, c = 6.4437 /nm); the forces handed back through the cross product match the energy's
// finite differences, and the virial matches dE/dλ for the atoms scaled by λ (the sites placed from them)
TEST(Water, Tip5pLonePairsForcesAndVirial) {
  caps::System s;
  {
    const caps::Trajectory w = caps::open_file(std::string(CAPS_SAMPLES) + "/water.pdb");
    const caps::Vec3 shift[2] = {{0, 0, 0}, {2.4, 1.3, 1.1}};
    for (int k = 0; k < 2; ++k) {
      caps::System one = w.frame(0);
      for (auto& a : one.atoms) {
        a.pos = a.pos + shift[k], a.mol = k + 1;
        if (k == 1) a.pos = caps::Vec3{a.pos[0], -a.pos[2] + 2.6, a.pos[1]};   // turned, so the pairs are not symmetric
      }
      const uint32_t base = uint32_t(s.atoms.size());
      for (auto& a : one.atoms) s.atoms.push_back(a);
      for (auto b : one.bonds) s.bonds.push_back({b.i + base, b.j + base, b.order});
    }
    s.has_mol = true;
  }
  const auto& m = caps::water_model("TIP5P");
  ASSERT_EQ(caps::apply_water_model(s, m), 2u);
  ASSERT_EQ(s.atoms.size(), 10u);
  const auto lp = caps::water_lp_coefficients(m);
  EXPECT_NEAR(lp[0], -0.344908, 2e-6);
  EXPECT_NEAR(lp[1] * 10, 6.4437, 1e-3);
  for (const auto& w : caps::find_waters(s)) {
    ASSERT_GE(w[4], 0);
    const auto O = s.atoms[size_t(w[0])].pos, H1 = s.atoms[size_t(w[1])].pos, H2 = s.atoms[size_t(w[2])].pos;
    const auto L1 = s.atoms[size_t(w[3])].pos, L2 = s.atoms[size_t(w[4])].pos;
    EXPECT_NEAR(caps::norm(L1 - O), 0.70, 1e-9);
    EXPECT_NEAR(caps::norm(L2 - O), 0.70, 1e-9);
    EXPECT_NEAR(std::acos(caps::dot(L1 - O, L2 - O) / 0.49) * 180 / M_PI, 109.47, 1e-6);
    EXPECT_LT(caps::dot(L1 + L2 - O * 2.0, H1 + H2 - O * 2.0), 0);   // the side away from the hydrogens
    EXPECT_NEAR(s.atoms[size_t(w[3])].charge, -0.241, 1e-12);
  }
  std::vector<uint32_t> all(s.atoms.size());
  for (size_t i = 0; i < all.size(); ++i) all[i] = uint32_t(i);
  const caps::ForceField F = caps::water_forcefield(s, m, all);
  ASSERT_EQ(F.vsites.size(), 4u);
  caps::EnergyOptions eo;
  eo.cutoff = 9.0;
  eo.electrostatics = caps::EnergyOptions::Electrostatics::DSF;
  caps::Evaluator ev(F, eo);
  std::vector<double> x(3 * s.atoms.size()), f;
  for (size_t i = 0; i < s.atoms.size(); ++i) for (int c = 0; c < 3; ++c) x[3 * i + c] = s.atoms[i].pos[c];
  const caps::Cell none;
  const auto e0 = ev.compute(x, none, f);
  const double h = 1e-5;
  double worst = 0;
  for (size_t i = 0; i < s.atoms.size(); ++i) {
    if (F.mass[i] == 0) continue;   // the sites carry no force of their own
    for (int c = 0; c < 3; ++c) {
      auto xp = x, xm = x;
      xp[3 * i + c] += h, xm[3 * i + c] -= h;
      std::vector<double> g;
      const double fd = -(ev.compute(xp, none, g).total() - ev.compute(xm, none, g).total()) / (2 * h);
      worst = std::max(worst, std::fabs(fd - f[3 * i + c]));
    }
  }
  // DSF's erfc is Abramowitz–Stegun 7.1.26 (|error| < 1.5e-7, as LAMMPS coul/dsf): energy and force agree to ~1e-5
  // (TIP4P/2005's linear sites in the same setup: 7.8e-5), independent of the step
  EXPECT_LT(worst, 5e-5) << "forces against finite differences";
  // virial: −dE/dλ at λ = 1 for every real atom scaled by λ (intramolecular terms included, as in Σ r·f)
  auto scaled = [&](double l) {
    auto y = x;
    for (auto& v : y) v *= l;
    std::vector<double> g;
    return ev.compute(y, none, g).total();
  };
  const double dEdl = (scaled(1 + h) - scaled(1 - h)) / (2 * h);
  EXPECT_NEAR(e0.virial, -dEdl, 1e-4 * std::max(1.0, std::fabs(dEdl)));   // without the λ² term: off by s·F
}
