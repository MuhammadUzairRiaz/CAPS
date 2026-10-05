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

// Chain stiffness and domains: k_θ (1 + cos θ) straightens chains (the mean cos θ at the middle beads moves toward −1),
// and a strongly segregated A10/B10 blend ends as one A and one B domain spanning the cell, the largest holding nearly
// all of its type's cells.
TEST(Dpd, StiffnessAndDomains) {
  auto mean_cos = [](const DpdReport& r) {
    const auto& t = r.frames;
    const auto& x = t.positions.back();
    const auto nb = t.topology.neighbours();
    double s = 0;
    int k = 0;
    const double L = t.cells.back().a[0];
    auto wrap = [&](Vec3 d) { for (int q = 0; q < 3; ++q) d[q] -= L * std::round(d[q] / L); return d; };
    for (size_t j = 0; j < nb.size(); ++j)
      if (nb[j].size() == 2) {
        const Vec3 a = wrap(x[nb[j][0]] - x[j]), b = wrap(x[nb[j][1]] - x[j]);
        s += dot(a, b) / (norm(a) * norm(b)), ++k;
      }
    return s / k;
  };
  DpdOptions o;
  o.species.push_back({"A10", "AAAAAAAAAA", 65});
  o.steps = 3000, o.equilibration = 1000, o.frame_every = 3000;
  const auto flex = run_dpd(o);
  o.angle_k = 5;
  const auto stiff = run_dpd(o);
  EXPECT_LT(mean_cos(stiff), mean_cos(flex) - 0.3);
  DpdOptions b;
  b.species = {{"A10", "AAAAAAAAAA", 33}, {"B10", "BBBBBBBBBB", 33}};
  b.chi["AB"] = 6;
  b.steps = 8000, b.equilibration = 2000, b.frame_every = 0;
  const auto seg = run_dpd(b);
  EXPECT_GT(seg.order, 0.5);   // the cells on the interface keep ψ below 1
  EXPECT_LE(seg.domains_a, 3);
  EXPECT_GT(seg.largest_a, 0.85);
  EXPECT_GT(seg.largest_b, 0.85);
  std::printf("DPD: <cos> %.2f flexible, %.2f stiff · blend psi %.2f, A domains %d (largest %.2f), B %d (%.2f)\n", mean_cos(flex), mean_cos(stiff), seg.order,
              seg.domains_a, seg.largest_a, seg.domains_b, seg.largest_b);
}
