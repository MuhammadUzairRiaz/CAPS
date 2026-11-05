#include <gtest/gtest.h>

#include <cmath>
#include <filesystem>
#include <fstream>

#include "caps/tension.hpp"

using namespace caps;

// A constant tensor: P_zz 200 atm, P_xx = P_yy 100 atm, L_z 60 Å, two interfaces: γ = 30 Å × 100 atm = 3000 atm·Å =
// 30.3975 mN/m exactly (1 atm·Å = 101325 Pa × 10⁻¹⁰ m); the error from five equal blocks is zero
TEST(Tension, ConstantTensor) {
  std::vector<TensionSample> s(50, TensionSample{200, 100, 100, 60});
  const auto r = surface_tension(s, 2, 0.2, 5);
  EXPECT_EQ(r.samples, 40);
  EXPECT_NEAR(r.gamma, 30.3975, 1e-12);
  EXPECT_NEAR(r.error, 0, 1e-12);
  EXPECT_NEAR(surface_tension(s, 1, 0.0, 5).gamma, 2 * 30.3975, 1e-12);   // one interface (a film on a wall)
}

// Alternating samples ±δ about the mean in blocks of two: the mean exact, the blocks' error zero; uneven blocks a
// positive error; the skipped share left out
TEST(Tension, BlocksAndSkip) {
  std::vector<TensionSample> s;
  for (int k = 0; k < 20; ++k) s.push_back({200.0 + (k % 2 ? 50 : -50), 100, 100, 60});
  const auto r = surface_tension(s, 2, 0.0, 5);
  EXPECT_NEAR(r.gamma, 30.3975, 1e-12);
  EXPECT_NEAR(r.error, 0, 1e-12);
  std::vector<TensionSample> t = s;
  for (int k = 0; k < 10; ++k) t[size_t(k)].pn += 40;   // the first half higher
  EXPECT_GT(surface_tension(t, 2, 0.0, 2).error, 0);
  EXPECT_NEAR(surface_tension(t, 2, 0.5, 5).gamma, 30.3975, 1e-12);   // the higher half skipped
}

// From thermo rows along each axis, and from a LAMMPS log with Pxx Pyy Pzz and Lz (or zhi − zlo)
TEST(Tension, FromRowsAndLog) {
  ThermoRow row;
  row.p[0] = 100, row.p[1] = 100, row.p[2] = 200, row.lx = 30, row.ly = 30, row.lz = 60;
  const auto z = tension_samples(std::vector<ThermoRow>(4, row), 2);
  EXPECT_NEAR(surface_tension(z, 2, 0).gamma, 30.3975, 1e-12);
  const auto x = tension_samples(std::vector<ThermoRow>(4, row), 0);   // normal x: P_xx − (P_yy + P_zz)/2 = −50 atm, L 30
  EXPECT_NEAR(surface_tension(x, 2, 0).gamma, 15 * -50 * 0.0101325, 1e-12);
  const ThermoLog log = parse_lammps_log("Step Temp Pxx Pyy Pzz Zlo Zhi\n0 300 100 100 200 -10 50\n100 300 100 100 200 -10 50\nLoop time of 1 on 1 procs\n");
  EXPECT_NEAR(surface_tension(tension_samples(log, 2), 2, 0).gamma, 30.3975, 1e-12);
  EXPECT_THROW(tension_samples(parse_lammps_log("Step Temp Press\n0 300 1\nLoop time of 1 on 1 procs\n"), 2), std::invalid_argument);
}

// CAPS's own thermo CSV (caps md --log) read back: the same γ as from the rows
TEST(Tension, FromCapsThermoCsv) {
  const auto p = (std::filesystem::temp_directory_path() / "caps_tension_thermo.csv").string();
  {
    std::ofstream f(p);
    f << "step,time_ps,stage,target_K,temperature_K,potential,kinetic,total,conserved,pressure_atm,pxx_atm,pyy_atm,pzz_atm,pxy_atm,pxz_atm,pyz_atm,lx_A,ly_A,lz_A,volume_A3,density_g_cm3\n";
    for (int k = 0; k < 6; ++k) f << k * 100 << ",0.1,md,300,300,-1,1,0,0,133,100,100,200,0,0,0,30,30,60,54000,0.8\n";
  }
  EXPECT_NEAR(surface_tension(tension_samples_from_file(p, 2), 2, 0).gamma, 30.3975, 1e-12);
}
