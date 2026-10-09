// Coarse-grained analyses (caps/cg_analysis.hpp): the entanglement estimators against their formulas (Hoy, Foteinopoulou &
// Kröger 2009, Eqs. 4–7, 13), config.Z1 written and read back as paths, the tension analysis on a curve with known modulus and
// hardening modulus, ⟨P₂⟩ of aligned bonds, g₃ of a drifting melt, the AA ↔ CG time factor.
#include <gtest/gtest.h>

#include <cmath>
#include <filesystem>
#include <fstream>

#include "caps/cg_analysis.hpp"

using namespace caps;

TEST(CgAnalysis, EstimatorFormulas) {
  CgChainPaths p;
  p.lpp = {100, 120}, p.r2 = {2000, 2400}, p.z = {3, 5}, p.beads = {50, 50};
  const CgEntanglement e = entanglement_of(p);
  const double N = 50, R2 = 2200, L = 110, L2 = (100.0 * 100 + 120.0 * 120) / 2, Z = 4;
  EXPECT_NEAR(e.ne_s_coil, (N - 1) * R2 / (L * L), 1e-9);                 // (4)
  EXPECT_NEAR(e.ne_s_kink, N * (N - 1) / (Z * (N - 1) + N), 1e-9);        // (5)
  EXPECT_NEAR(e.ne_mod_s_kink, N / Z, 1e-9);                              // (6)
  EXPECT_NEAR(e.ne_mod_s_coil, (N - 1) / (L2 / R2 - 1), 1e-9);            // (7)
  EXPECT_NEAR(e.a_pp, R2 / L, 1e-9);
  // M-kink: ⟨Z⟩ = N/40 − 0.5 exactly → N_e = 40
  std::vector<CgEntanglement> sets;
  for (double n : {100.0, 200.0, 400.0}) { CgEntanglement s; s.N = n; s.z = n / 40 - 0.5; s.lpp = 10; sets.push_back(s); }
  EXPECT_NEAR(multi_estimators(sets, {}, 0).ne_m_kink, 40.0, 1e-9);
}

TEST(CgAnalysis, Z1ConfigAndPaths) {
  CgTopology t;
  t.chains = {{0, 1, 2}, {3, 4, 5, 6}};
  for (int i = 0; i < 7; ++i) t.kind.push_back("A"), t.mol.push_back(i < 3 ? 0 : 1);
  t.bonds = {{0, 1}, {1, 2}, {3, 4}, {4, 5}, {5, 6}};
  Cell c;
  c.a = {20, 0, 0}, c.b = {0, 20, 0}, c.c = {0, 0, 20};
  // the second chain crosses the box edge: written whole
  std::vector<Vec3> pos = {{1, 1, 1}, {4, 1, 1}, {7, 1, 1}, {18, 5, 5}, {1, 5, 5}, {4, 5, 5}, {7, 5, 5}};
  const auto dir = std::filesystem::temp_directory_path() / "caps_cg_analysis_test";
  std::filesystem::create_directories(dir);
  const std::string f = (dir / "config.Z1").string();
  write_z1_config(t, pos, c, f);
  std::ifstream in(f);
  std::string l1, l2, l3, l4;
  std::getline(in, l1), std::getline(in, l2), std::getline(in, l3);
  EXPECT_EQ(l1, "2");
  EXPECT_EQ(l3, "3 4");   // every chain's length on the third line
  for (int k = 0; k < 4; ++k) std::getline(in, l4);
  EXPECT_EQ(l4, "18.000000 5.000000 5.000000");
  std::getline(in, l4);
  EXPECT_EQ(l4, "21.000000 5.000000 5.000000");   // made whole across the edge
  // the same layout read as paths: chain 1 straight (Z = 1 interior node), chain 2 9 Å long
  const CgChainPaths p = read_z1_paths(f);
  ASSERT_EQ(p.lpp.size(), 2u);
  EXPECT_NEAR(p.lpp[0], 6, 1e-9);
  EXPECT_NEAR(p.lpp[1], 9, 1e-9);
  EXPECT_EQ(p.z[0], 1);
  EXPECT_EQ(p.z[1], 2);
  // a PPA dump read back: contour of the paths, ends from the start
  std::vector<Vec3> after = pos;
  after[1] = {4, 3, 1};
  const CgChainPaths q = paths_from_positions(t, pos, after, c);
  EXPECT_NEAR(q.lpp[0], 2 * std::sqrt(9.0 + 4.0), 1e-9);
  EXPECT_NEAR(q.r2[0], 36, 1e-9);
  std::filesystem::remove_all(dir);
}

TEST(CgAnalysis, TensionCurve) {
  // σ = E ε up to 5 % (E = 2000 MPa), a drop, then σ = 60 + G_R (λ² − 1/λ − h0) with G_R = 12 MPa after ε = 0.3
  CgStressStrain c;
  for (int k = 0; k <= 3000; ++k) {
    const double e = k * 0.001, lam = 1 + e;
    double s;
    if (e <= 0.05) s = 2000 * e;
    else if (e <= 0.15) s = 100 - 400 * (e - 0.05);
    else if (e <= 0.3) s = 60;
    else s = 60 + 12 * ((lam * lam - 1 / lam) - (1.3 * 1.3 - 1 / 1.3));
    c.strain.push_back(e), c.stress.push_back(s);
    for (auto* v : {&c.bond, &c.angle, &c.dihedral, &c.pair, &c.kinetic, &c.density}) v->push_back(0);
  }
  const CgTensionAnalysis a = analyse_tension(c, 0.02, 0.3);
  EXPECT_NEAR(a.modulus, 2000, 1);
  EXPECT_NEAR(a.yield_stress, 100, 3);
  EXPECT_NEAR(a.yield_strain, 0.05, 0.01);
  EXPECT_NEAR(a.softening, 40, 3);
  EXPECT_NEAR(a.hardening_modulus, 12, 0.05);
}

TEST(CgAnalysis, OrientationAndDynamics) {
  // bonds along z: ⟨P₂⟩ = 1; chains along z: large anisotropy
  CgTopology t;
  std::vector<Vec3> pos;
  for (int ch = 0; ch < 4; ++ch) {
    std::vector<int> c;
    for (int k = 0; k < 5; ++k) {
      c.push_back(int(pos.size()));
      t.kind.push_back("A"), t.mol.push_back(ch);
      if (k) t.bonds.push_back({int(pos.size()) - 1, int(pos.size())});
      pos.push_back({double(ch) * 5 + 1, 1.0 + 0.01 * k, double(k) * 3 + 1});
    }
    t.chains.push_back(c);
  }
  Cell cell;
  cell.a = {40, 0, 0}, cell.b = {0, 40, 0}, cell.c = {0, 0, 40};
  const CgOrientation o = orientation_of(t, pos, cell, 1.0);
  EXPECT_NEAR(o.p2, 1.0, 1e-3);
  EXPECT_GT(o.ree_anisotropy, 100);
  EXPECT_GT(o.void_fraction, 0.9);
  // the melt drifting rigidly along x at 0.01 Å/fs: g₁ = g₃ = (0.01 t)², g₂ = 0, P₁ = 1
  std::vector<std::vector<Vec3>> frames;
  std::vector<Cell> cells;
  std::vector<double> times;
  for (int f = 0; f < 30; ++f) {
    std::vector<Vec3> p = pos;
    for (auto& v : p) v[0] += 0.01 * 100 * f;
    frames.push_back(p), cells.push_back(cell), times.push_back(100.0 * f);
  }
  const CgDynamics D = cg_dynamics(t, frames, cells, times);
  ASSERT_GT(D.t.size(), 5u);
  for (size_t k = 0; k < D.t.size(); ++k) {
    EXPECT_NEAR(D.g3[k], std::pow(0.01 * D.t[k], 2), 1e-6);
    EXPECT_NEAR(D.g1[k], D.g3[k], 1e-6);
    EXPECT_NEAR(D.g2[k], 0, 1e-9);
    EXPECT_NEAR(D.p1[k], 1, 1e-9);
  }
  // time mapping: g₁ = t^0.6 against the same with times ×5
  std::vector<double> ta, ga, tc, gc;
  for (int k = 1; k <= 50; ++k) {
    const double tt = std::pow(10.0, k / 10.0);
    tc.push_back(tt), gc.push_back(std::pow(tt, 0.6));
    ta.push_back(5 * tt), ga.push_back(std::pow(tt, 0.6));
  }
  double spread = 1;
  EXPECT_NEAR(time_mapping(ta, ga, tc, gc, &spread), 5.0, 1e-6);
  EXPECT_LT(spread, 1e-9);
}
