#include <gtest/gtest.h>

#include <cmath>
#include <random>
#include <stdexcept>
#include <vector>

#include "caps/exchange_analysis.hpp"

using namespace caps;

namespace {

constexpr double kR = 8.314462618e-3;   // kJ/(mol K)

Cell box(double x, double y, double z) {
  Cell c;
  c.a = {x, 0, 0};
  c.b = {0, y, 0};
  c.c = {0, 0, z};
  return c;
}

}  // namespace

TEST(Exchange, ArrheniusRecoversTheActivationEnergy) {
  std::vector<double> T, tau;
  for (double t = 300; t <= 400; t += 20) T.push_back(t), tau.push_back(1e-3 * std::exp(80.0 / (kR * t)));
  const auto f = fit_arrhenius(T, tau);
  ASSERT_TRUE(f.ok);
  EXPECT_NEAR(f.ea, 80.0, 1e-8);
  EXPECT_NEAR(f.ln_tau0, std::log(1e-3), 1e-8);
  EXPECT_NEAR(f.r2, 1.0, 1e-12);
  EXPECT_EQ(f.points, 6u);
  // a window leaves points out; reduced units give Ea in the temperature's unit
  EXPECT_EQ(fit_arrhenius(T, tau, false, 320, 380).points, 4u);
  const auto r = fit_arrhenius({0.5, 0.6, 0.8}, {std::exp(3.24 / 0.5), std::exp(3.24 / 0.6), std::exp(3.24 / 0.8)}, true);
  EXPECT_NEAR(r.ea, 3.24, 1e-10);
  // one temperature is not a fit
  EXPECT_FALSE(fit_arrhenius({300, 300}, {1, 2}).ok);
  const auto props = arrhenius_properties(f);
  ASSERT_EQ(props.size(), 1u);
  EXPECT_EQ(props[0].id, "arrhenius");
  EXPECT_NEAR(props[0].extra.at("Ea (kcal/mol)"), 80.0 / 4.184, 1e-8);
}

TEST(Exchange, CrossingTemperatureOfTwoArrheniusLines) {
  // ln τa = 5000/T − ln 1e6, ln τb = 2000/T − ln 1e2: equal where 3000/T = ln 1e4
  std::vector<double> T, a, b;
  for (double t = 250; t <= 450; t += 25) T.push_back(t), a.push_back(1e-6 * std::exp(5000 / t)), b.push_back(1e-2 * std::exp(2000 / t));
  EXPECT_NEAR(crossing_temperature(T, a, b), 3000 / std::log(1e4), 1e-9);
  EXPECT_TRUE(std::isnan(crossing_temperature({300, 400}, {1, 1}, {2, 3})));
}

TEST(Exchange, KwwFitOfAnExactStretchedExponential) {
  std::vector<double> t{0}, y{1};   // the fit normalises by the first point
  for (int k = 0; k < 400; ++k) t.push_back(0.01 * std::pow(1.03, k)), y.push_back(std::exp(-std::pow(t.back() / 5.0, 0.6)));
  const auto kw = fit_kww(t, y, 0, 0.05, 0.95);
  EXPECT_NEAR(kw[0], 5.0, 1e-6);
  EXPECT_NEAR(kw[1], 0.6, 1e-8);
  EXPECT_NEAR(kw[2], 5.0 / 0.6 * std::tgamma(1 / 0.6), 1e-5);
}

TEST(Exchange, DynamicModuliOfAMaxwellModel) {
  // G = exp(−t), G' = ω²/(1+ω²), G'' = ω/(1+ω²)
  std::vector<double> t{0}, g{1};
  for (int k = 0; k < 3000; ++k) t.push_back(1e-5 * std::pow(1.006, k)), g.push_back(std::exp(-t.back()));
  ASSERT_GT(t.back(), 50);
  const std::vector<double> w{0.05, 0.3, 1, 3, 20};
  std::vector<double> gp, gpp;
  dynamic_moduli(t, g, w, 0, gp, gpp);
  for (size_t k = 0; k < w.size(); ++k) {
    EXPECT_NEAR(gp[k], w[k] * w[k] / (1 + w[k] * w[k]), 2e-4) << "ω " << w[k];
    EXPECT_NEAR(gpp[k], w[k] / (1 + w[k] * w[k]), 2e-4) << "ω " << w[k];
  }
  // a solid: G = G_eq after the end adds G_eq to G' at low frequency
  std::vector<double> gp2, gpp2;
  dynamic_moduli({0, 1, 2}, {3, 3, 3}, {1e-3}, 3, gp2, gpp2);
  EXPECT_NEAR(gp2[0], 3, 1e-9);
  EXPECT_NEAR(gpp2[0], 0, 1e-9);
}

TEST(Exchange, MultiTauCorrelatorOfAnExponentialSignal) {
  // AR(1): ⟨x(0) x(k)⟩ = exp(−k/τ) for unit variance
  const double tau = 20, a = std::exp(-1 / tau);
  std::mt19937_64 rng(7);
  std::normal_distribution<double> n(0, 1);
  MultiTauCorrelator c(1, 16, 2);
  double x = 0;
  for (int k = 0; k < 400000; ++k) {
    x = a * x + std::sqrt(1 - a * a) * n(rng);
    c.add({x});
  }
  std::vector<double> lag, corr;
  c.result(lag, corr);
  ASSERT_GT(lag.size(), 40u);
  EXPECT_EQ(lag[0], 0);
  int checked = 0;
  for (size_t k = 0; k < lag.size(); ++k)
    if (lag[k] <= 80) {
      EXPECT_NEAR(corr[k], std::exp(-lag[k] / tau), 0.04) << "lag " << lag[k];
      ++checked;
    }
  EXPECT_GT(checked, 20);
  // lags grow geometrically past the first level
  EXPECT_EQ(lag[15], 15);
  EXPECT_EQ(lag[16], 16);
  EXPECT_EQ(lag[17], 18);
}

TEST(Exchange, StressRelaxationModulusUnitsAndDecay) {
  // independent AR(1) shear stresses of σ atm, no normal stress: Σ P°P° = 2·3σ² at t = 0, so G(0) = 0.6 V σ²/(k_B T)
  const double sigma = 100, tau = 20, a = std::exp(-1 / tau);
  std::mt19937_64 rng(11);
  std::normal_distribution<double> n(0, 1);
  std::vector<std::array<double, 6>> p(200000);
  std::array<double, 3> s{0, 0, 0};
  for (auto& q : p) {
    for (auto& v : s) v = a * v + std::sqrt(1 - a * a) * sigma * n(rng);
    q = {1.0, 1.0, 1.0, s[0], s[1], s[2]};
  }
  StressRelaxationOptions o;
  o.dt_fs = 10, o.volume = 1e5, o.temperature = 300, o.fit_from_ps = 0;
  const auto r = stress_relaxation(p, o);
  const double g0 = 0.6 * o.volume * 1e-30 * std::pow(sigma * 101325.0, 2) / (1.380649e-23 * o.temperature) * 1e-6;
  EXPECT_NEAR(r.g0 / g0, 1.0, 0.03);
  EXPECT_NEAR(r.kww_beta, 1.0, 0.1);
  EXPECT_NEAR(r.kww_tau, tau * o.dt_fs * 1e-3, 0.03);
  // the Green–Kubo integral of G0 e^(−t/τ) is G0 τ (MPa·ps → mPa·s)
  EXPECT_NEAR(r.viscosity / (g0 * tau * o.dt_fs * 1e-3 * 1e-3), 1.0, 0.08);
  ASSERT_FALSE(r.omega.empty());
  EXPECT_EQ(r.storage.size(), r.omega.size());
  EXPECT_THROW(stress_relaxation(p, StressRelaxationOptions{}), std::invalid_argument);   // no volume
  EXPECT_EQ(stress_relaxation_properties(r)[0].series.size(), 3u);
}

TEST(Exchange, TwoStepFitRecoversItsParameters) {
  std::vector<double> t, f;
  for (int k = 0; k < 300; ++k) {
    t.push_back(0.01 * std::pow(1.05, k));
    f.push_back(0.5 * (1 + std::exp(-t.back() / 1.0)) * std::exp(-std::pow(t.back() / 100.0, 0.6)));
  }
  const auto r = fit_two_step(t, f);
  ASSERT_TRUE(r.ok);
  EXPECT_NEAR(r.tau, 1.0, 0.02);
  EXPECT_NEAR(r.tau_d, 100.0, 2.0);
  EXPECT_NEAR(r.beta, 0.6, 0.01);
  EXPECT_LT(r.rms, 1e-4);
}

TEST(Exchange, BondMemoryCountsPartnersAndChanges) {
  // site 0 keeps its partner, site 1 alternates between two partners every frame, site 2 is always free
  std::vector<std::vector<int64_t>> pt;
  const int nf = 40;
  for (int f = 0; f < nf; ++f) pt.push_back({100, f % 2 ? 201 : 200, -1});
  const auto r = bond_memory(pt, 0.5);
  ASSERT_GE(r.c.size(), 3u);
  EXPECT_DOUBLE_EQ(r.c[0], 1.0);
  EXPECT_DOUBLE_EQ(r.c[1], 0.5);
  EXPECT_DOUBLE_EQ(r.c[2], 1.0);
  EXPECT_DOUBLE_EQ(r.t[1], 0.5);
  EXPECT_EQ(r.changes, size_t(nf - 1));
  EXPECT_NEAR(r.change_time, 2 * 0.5, 1e-12);   // two bonded sites per frame, one change per frame
  EXPECT_TRUE(std::isnan(r.tau_e));
  EXPECT_THROW(bond_memory({{1, 2}, {1}}, 1.0), std::invalid_argument);
  // independent exchange at a fixed rate per frame: c(t) = (1 − k)^t, 1/e at t = −1/ln(1 − k)
  std::mt19937_64 rng(3);
  std::uniform_real_distribution<double> u(0, 1);
  const double k = 0.02;
  std::vector<std::vector<int64_t>> ex(600, std::vector<int64_t>(400));
  std::vector<int64_t> cur(400);
  int64_t next = 0;
  for (auto& c : cur) c = next++;
  for (auto& fr : ex) {
    for (auto& c : cur)
      if (u(rng) < k) c = next++;
    fr = cur;
  }
  const auto m = bond_memory(ex, 1.0);
  EXPECT_NEAR(m.tau_e, -1 / std::log(1 - k), 2.0);
  EXPECT_NEAR(m.kww_beta, 1.0, 0.05);
}

TEST(Exchange, InterfaceBondsAndCrossingProfile) {
  System s;
  s.cell = box(10, 10, 20);
  auto atom = [&](double z) {
    Atom a;
    a.element = 6;
    a.pos = {5, 5, z};
    s.atoms.push_back(a);
  };
  atom(9.5), atom(10.5);   // across the cut at z = 10
  atom(19.5), atom(0.5);   // across the periodic boundary at z = 0
  atom(3.2), atom(3.8);    // inside the lower half
  s.bonds = {{0, 1}, {2, 3}, {4, 5}};
  const auto side = sides_from_positions(s, 2, 0.5);
  EXPECT_EQ(side[0], 0);
  EXPECT_EQ(side[1], 1);
  EXPECT_EQ(side[2], 1);
  EXPECT_EQ(side[3], 0);
  const auto r = interface_bonds(s, side, 2, 0.5, 1.0);
  EXPECT_EQ(r.bonds, 3u);
  EXPECT_EQ(r.bridges, 2u);
  EXPECT_EQ(r.crossing_plane, 1u);
  ASSERT_EQ(r.crossings.size(), 20u);
  EXPECT_DOUBLE_EQ(r.crossings[0], 1.0);    // z = 0
  EXPECT_DOUBLE_EQ(r.crossings[10], 1.0);   // z = 10
  double total = 0;
  for (double c : r.crossings) total += c;
  EXPECT_DOUBLE_EQ(total, 2.0);
  // an atom left out takes no part in bridges
  auto side2 = side;
  side2[1] = -1;
  EXPECT_EQ(interface_bonds(s, side2, 2, 0.5, 1.0).bridges, 1u);
}

TEST(Exchange, HealingEfficiencyDefinitions) {
  StressStrain intact, healed;
  for (int k = 0; k <= 100; ++k) intact.strain.push_back(0.001 * k), intact.stress.push_back(100 * 0.001 * k);
  for (int k = 0; k <= 50; ++k) healed.strain.push_back(0.001 * k), healed.stress.push_back(80 * 0.001 * k);
  const auto h = healing_efficiency(healed, intact, 0.03, 0.04);
  EXPECT_NEAR(h.toughness, 100 * (0.5 * 4 * 0.05) / (0.5 * 10 * 0.1), 1e-9);   // 20 %
  EXPECT_NEAR(h.peak, 40, 1e-9);
  EXPECT_NEAR(h.modulus, 80, 1e-9);
  EXPECT_NEAR(h.stress_at, 80, 1e-9);
  const auto h2 = healing_efficiency(healed, intact, 0.03, 0, 0.05);
  EXPECT_NEAR(h2.toughness, 80, 1e-9);   // both areas to strain 0.05
  EXPECT_TRUE(std::isnan(h2.stress_at));
  const auto props = healing_properties(h);
  EXPECT_EQ(props[0].id, "healing");
  EXPECT_NEAR(props[0].extra.at("peak-stress ratio (%)"), 40, 1e-9);
  EXPECT_THROW(healing_efficiency(StressStrain{{0}, {0}}, intact), std::invalid_argument);
}
