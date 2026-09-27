#include <gtest/gtest.h>

#include <algorithm>
#include <cstdio>

#include "caps/dpd.hpp"

using namespace caps;

// A one-component DPD fluid: without repulsion the pressure is the ideal gas's ρkT exactly; at ρ = 3, a = 25 (the water
// model) p = 23.65 as Groot & Warren's simulations give, with kT held at 1; (p − ρkT)/(aρ²) approaches their high-density
// α = 0.101 (J. Chem. Phys. 107, 4423 (1997)) — 0.099 at ρ = 6
TEST(Dpd, EquationOfStateAndTemperature) {
  auto run = [](double rho, double a) {
    DpdOptions o;
    o.density = rho;
    o.species.push_back({"water", "A", int(rho * 216)});
    o.a["AA"] = a;
    o.steps = 4000;
    o.equilibration = 1000;
    o.frame_every = 0;
    return run_dpd(o);
  };
  const auto ideal = run(3, 0);
  EXPECT_NEAR(ideal.pressure / (3 * ideal.kT), 1.0, 0.005);
  const auto water = run(3, 25);
  EXPECT_NEAR(water.kT, 1.0, 0.01);
  EXPECT_NEAR(water.pressure, 23.65, 0.01 * 23.65);
  const auto dense = run(6, 25);
  EXPECT_NEAR((dense.pressure - 6 * dense.kT) / (25 * 36), 0.101, 0.003);
  DpdOptions d;
  EXPECT_NEAR(dpd_repulsion(d, 'A', 'A'), 25.0, 1e-12);   // 75 kT / ρ
}

// A symmetric diblock A5B5: with a_AB = 25 (χ = 0) the melt stays mixed; with χ = 4.3 (χN = 43, far past the
// order–disorder transition near χN ≈ 10) the blocks segregate into microdomains: the order parameter rises and S(q) of A
// peaks at a finite q
TEST(Dpd, DiblockMicrophaseSeparation) {
  DpdOptions o;
  o.species.push_back({"A5B5", "AAAAABBBBB", 154});
  o.steps = 15000;
  o.equilibration = 5000;
  o.frame_every = 5000;
  o.chi["AB"] = 0.0;
  const auto mixed = run_dpd(o);
  o.chi["AB"] = 4.3;
  const auto seg = run_dpd(o);
  std::printf("diblock: ψ mixed %.3f, segregated %.3f · q* %.2f (spacing %.2f r_c) · peak S %.1f\n", mixed.order, seg.order, seg.q_peak, seg.spacing,
              seg.sq.empty() ? 0.0 : *std::max_element(seg.sq.begin(), seg.sq.end()));
  EXPECT_NEAR(dpd_repulsion(o, 'A', 'B'), 25 + 4.3 / 0.286, 1e-9);
  EXPECT_LT(mixed.order, 0.35);
  EXPECT_GT(seg.order, 0.6);
  EXPECT_GT(*std::max_element(seg.sq.begin(), seg.sq.end()), 5 * *std::max_element(mixed.sq.begin(), mixed.sq.end()));
  EXPECT_EQ(seg.frames.frames(), 4u);
  EXPECT_EQ(seg.frames.topology.atoms.size(), 1540u);
}
