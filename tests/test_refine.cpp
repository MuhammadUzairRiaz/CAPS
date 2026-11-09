#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>

#include "caps/crystal.hpp"
#include "caps/powder.hpp"
#include "caps/refine.hpp"

using namespace caps;

namespace {
const std::string kData = std::string(CAPS_SOURCE_DIR) + "/data/crystals/";
RefinePhase phase(const std::string& file, const std::string& name) {
  RefinePhase p;
  p.name = name;
  p.cell = read_cif(kData + file);
  return p;
}
std::vector<double> grid() {
  std::vector<double> t;
  for (double x = 20; x <= 100.0001; x += 0.02) t.push_back(x);
  return t;
}
// a phase with its cell scaled to lattice constant a (cubic)
RefinePhase scaled(RefinePhase p, double a) {
  const double f = a / norm(p.cell.cell.a);
  p.cell.cell.a = p.cell.cell.a * f, p.cell.cell.b = p.cell.cell.b * f, p.cell.cell.c = p.cell.cell.c * f;
  for (auto& at : p.cell.atoms) at.pos = at.pos * f;
  return p;
}
}  // namespace

// The refinement's calculated pattern is the powder pattern of powder.hpp (same physics, normalised alike)
TEST(Refine, ForwardModelIsThePowderPattern) {
  const RefinePhase cu = phase("copper.cif", "Cu");
  RefineOptions o;
  o.W = 0.01, o.Y = 0.02;
  const auto t = grid();
  auto yc = calculated_pattern(t, {cu}, {1.0}, o);
  PowderOptions po;
  po.W = 0.01, po.Y = 0.02, po.two_theta_min = 20, po.two_theta_max = 100, po.step = 0.02;
  Trajectory tr;
  tr.topology = cu.cell;
  const PowderReport pr = powder_pattern(tr, {}, po);
  const double m1 = *std::max_element(yc.begin(), yc.end());
  double worst = 0;
  for (size_t i = 0; i < t.size() && i < pr.intensity.size(); ++i) worst = std::max(worst, std::fabs(100 * yc[i] / m1 - pr.intensity[i]));
  std::printf("largest difference %.2e (of 100)\n", worst);
  EXPECT_LT(worst, 0.05);
}

// Rietveld recovers the lattice, zero shift and widths that made a synthetic copper pattern (background added)
TEST(Refine, RietveldRecoversSyntheticCopper) {
  const RefinePhase cu = phase("copper.cif", "Cu");
  RefineOptions truth;
  truth.W = 0.012, truth.Y = 0.03, truth.zero = 0.03;
  const auto t = grid();
  auto y = calculated_pattern(t, {scaled(cu, 3.6175)}, {0.02}, truth);
  for (size_t i = 0; i < t.size(); ++i) y[i] += 40 + 0.3 * (t[i] - 20);
  MeasuredPattern d{t, y, {}};
  RefineOptions o;
  o.W = 0.02, o.Y = 0.02, o.zero = 0, o.refine_U = false;
  const RefineResult r = refine_pattern(d, {cu}, o);
  std::printf("a %.6f (true 3.6175 ± %.1e) · zero %.5f · W %.6f · Y %.6f · Rwp %.2e after %d iterations\n", r.phases[0].lattice[0], r.phases[0].lattice_esd[0], r.zero, r.W,
              r.Y, r.rwp, r.iterations);
  EXPECT_NEAR(r.phases[0].lattice[0], 3.6175, 1e-6);
  EXPECT_NEAR(r.phases[0].lattice[2], 3.6175, 1e-6);   // the cubic equality kept
  EXPECT_NEAR(r.zero, 0.03, 1e-6);
  EXPECT_NEAR(r.W, 0.012, 1e-6);
  EXPECT_NEAR(r.Y, 0.03, 1e-6);
  EXPECT_LT(r.rwp, 1e-4);
}

// Pawley recovers the lattice and the reflections' intensity ratio without the structure
TEST(Refine, PawleyRecoversIntensities) {
  const RefinePhase cu = phase("copper.cif", "Cu");
  RefineOptions truth;
  truth.W = 0.012, truth.Y = 0.03;
  const auto t = grid();
  auto y = calculated_pattern(t, {scaled(cu, 3.6175)}, {0.02}, truth);
  for (size_t i = 0; i < t.size(); ++i) y[i] += 40;
  MeasuredPattern d{t, y, {}};
  RefineOptions o;
  o.mode = "pawley", o.W = 0.02, o.Y = 0.02, o.refine_U = false, o.refine_zero = false;
  const RefineResult p = refine_pattern(d, {cu}, o);
  RefineOptions orr = o;
  orr.mode = "rietveld";
  const RefineResult r = refine_pattern(d, {cu}, orr);
  EXPECT_NEAR(p.phases[0].lattice[0], 3.6175, 1e-6);
  EXPECT_LT(p.rwp, 1e-5);
  // Pawley's intensity per group is the structure's Σ m|F|² × scale; Rietveld's listed one includes LP: their ratios agree
  ASSERT_GE(p.phases[0].intensity.size(), 2u);
  const double lam = 1.540593;   // Cu Kα1 (powder.hpp's table)
  auto lp = [&](double tt) {
    const double th = tt / 2 * M_PI / 180, s = std::sin(th), c = std::cos(th), c2 = std::cos(2 * th);
    return (1 + c2 * c2) / 2 / (s * s * c);
  };
  (void)lam;
  const double rp = p.phases[0].intensity[0] / p.phases[0].intensity[1];
  const double rr = (r.phases[0].intensity[0] / lp(r.phases[0].two_theta[0])) / (r.phases[0].intensity[1] / lp(r.phases[0].two_theta[1]));
  std::printf("(111)/(200) Pawley %.6f, from the structure %.6f\n", rp, rr);
  EXPECT_NEAR(rp, rr, 1e-5 * rr);
}

// Weight fractions of a copper + rock-salt mixture made with 30 % / 70 % come back from the scale factors (Hill & Howard)
TEST(Refine, QuantitativePhaseAnalysis) {
  const RefinePhase cu = phase("copper.cif", "Cu"), nacl = phase("rock-salt.cif", "NaCl");
  auto zmv = [](const RefinePhase& p) {
    double m = 0;
    for (const auto& a : p.cell.atoms) m += p.cell.mass_of(a);
    return m * p.cell.cell.volume();
  };
  const double s1 = 0.3 / zmv(cu) * 1e3, s2 = 0.7 / zmv(nacl) * 1e3;
  RefineOptions truth;
  truth.W = 0.012, truth.Y = 0.03;
  const auto t = grid();
  auto y = calculated_pattern(t, {cu, nacl}, {s1, s2}, truth);
  for (size_t i = 0; i < t.size(); ++i) y[i] += 20;
  MeasuredPattern d{t, y, {}};
  RefineOptions o;
  o.W = 0.015, o.Y = 0.025, o.refine_U = false, o.refine_zero = false;
  const RefineResult r = refine_pattern(d, {cu, nacl}, o);
  std::printf("weight fractions %.6f %.6f · Rwp %.2e\n", r.phases[0].weight_fraction, r.phases[1].weight_fraction, r.rwp);
  EXPECT_NEAR(r.phases[0].weight_fraction, 0.3, 1e-6);
  EXPECT_NEAR(r.phases[1].weight_fraction, 0.7, 1e-6);
}
