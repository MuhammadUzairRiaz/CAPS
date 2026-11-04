#include <gtest/gtest.h>

#include <cmath>
#include <memory>

#include <algorithm>
#include "caps/uff.hpp"
#include "caps/zmethod.hpp"

using namespace caps;

namespace {
// 256 Lennard-Jones argon atoms (ε/k_B = 119.8 K, σ = 3.405 Å) on an fcc lattice at reduced density rho
System lj_crystal(double rho) {
  const double sigma = 3.405, a = std::cbrt(4 / rho) * sigma;
  System s;
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j)
      for (int k = 0; k < 4; ++k)
        for (const auto& b : std::vector<Vec3>{{0, 0, 0}, {0.5, 0.5, 0}, {0.5, 0, 0.5}, {0, 0.5, 0.5}}) {
          Atom at;
          at.element = 18;
          at.mol = int64_t(s.atoms.size() + 1);
          at.pos = {(i + b[0]) * a, (j + b[1]) * a, (k + b[2]) * a};
          s.atoms.push_back(at);
        }
  s.cell.a = {4 * a, 0, 0}, s.cell.b = {0, 4 * a, 0}, s.cell.c = {0, 0, 4 * a};
  return s;
}

ZMethodOptions lj_options(const System& s) {
  ZMethodOptions o;
  auto ff = std::make_shared<ForceField>(default_forcefield(s));
  for (auto& t : ff->lj) t = {119.8 * 0.0019872043, 3.405};
  std::fill(ff->charge.begin(), ff->charge.end(), 0.0);
  o.field = ff;
  o.energy.cutoff = 2.5 * 3.405;
  o.energy.coulomb = false;
  o.energy.tail = false;
  o.dt = 5;
  return o;
}
}  // namespace

// The perfect lattice at rest given velocities at T0 shares the energy equally between kinetic and potential (virial
// theorem of a harmonic crystal): it settles at T0 / 2, and its Bragg order stays near 1 (Debye–Waller); the order
// parameter of the perfect cell is exactly 1 on fcc's {111} star
TEST(ZMethod, HarmonicLimitHalvesTheTemperature) {
  const System s = lj_crystal(1.0);
  ZMethodOptions o = lj_options(s);
  o.temperatures = {10};
  o.steps = 2000;
  const ZResult r = z_method(s, o);
  ASSERT_EQ(r.points.size(), 1u);
  EXPECT_NEAR(r.order0, 1.0, 1e-9);
  EXPECT_EQ(r.g.size(), 4u);   // the four {111} vectors, one of each ± pair
  const auto& p = r.points[0];
  std::printf("T0 10 K: T %.3f ± %.3f K, order %.4f, drift %.2e kcal/mol per atom\n", p.temperature, p.temperature_sd, p.order, p.drift);
  EXPECT_NEAR(p.temperature, 5.0, 0.15);
  EXPECT_GT(p.order, 0.95);
  EXPECT_FALSE(p.melted);
  EXPECT_LT(std::fabs(p.drift), 1e-4);
}

// The Z of the Lennard-Jones crystal at ρ* = 0.95: T rises with E on the crystal branch, the first melted run (found by
// bisecting the starting temperature) lies below the hottest crystal run in temperature, and its order has fallen to the
// liquid's; eleven runs of 4000 steps
TEST(ZMethod, LennardJonesCrystalTracesAZ) {
  const System s = lj_crystal(0.95);
  ZMethodOptions o = lj_options(s);
  o.temperatures = {100, 200, 250, 300, 350, 400, 450};
  o.bisect = 4;
  o.steps = 4000;
  const ZResult r = z_method(s, o);
  for (const auto& p : r.points)
    std::printf("T0 %5.0f  E %9.5f  T %7.2f  P %8.0f  order %.3f  msd %6.2f  %s\n", p.t0, p.energy, p.temperature, p.pressure, p.order, p.msd, p.melted ? "melted" : "crystal");
  for (const auto& n : r.notes) std::printf("%s\n", n.c_str());
  ASSERT_GT(r.t_melt, 0);
  ASSERT_GT(r.t_superheat, 0);
  EXPECT_LT(r.t_melt, r.t_superheat);
  EXPECT_GT(r.e_melt, r.e_superheat);
  // the crystal branch rises up to the superheating limit (a run past it may be caught while it melts)
  double last = 0;
  for (const auto& p : r.points) {
    if (p.energy > r.e_superheat) break;
    EXPECT_GT(p.temperature, last);
    last = p.temperature;
  }
  for (const auto& p : r.points)
    if (p.melted) EXPECT_LT(p.order, 0.2);
}
