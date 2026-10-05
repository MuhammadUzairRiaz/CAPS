#include <gtest/gtest.h>

#include <cmath>
#include <cstdio>
#include <random>

#include "caps/pairmodel.hpp"
#include "caps/sorption.hpp"

using namespace caps;

namespace {
// 27 LJ host atoms on a 4 Å cubic lattice in a 12 Å cell, then one LJ sorbate atom (the template)
struct Host {
  System s;
  ForceField ff;
};
Host make_host(double eps_host, double eps_sorbate) {
  Host h;
  h.s.cell.a = {12, 0, 0}, h.s.cell.b = {0, 12, 0}, h.s.cell.c = {0, 0, 12};
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      for (int k = 0; k < 3; ++k) {
        Atom a;
        a.element = 6;
        a.pos = {4.0 * i, 4.0 * j, 4.0 * k};
        h.s.atoms.push_back(a);
        h.ff.type_index.push_back(0);
      }
  Atom g;
  g.element = 18;
  g.pos = {2, 2, 2};
  h.s.atoms.push_back(g);
  h.ff.type_index.push_back(1);
  h.ff.lj = {{eps_host, 3.0}, {eps_sorbate, 3.4}};
  h.ff.type_names = {"C", "Ar"};
  h.ff.charge.assign(h.s.atoms.size(), 0.0);
  h.ff.mass.assign(h.s.atoms.size(), 12.011);
  h.ff.mass.back() = 39.948;
  h.ff.name = "test LJ";
  return h;
}
}  // namespace

// Widom: ⟨e^(−βU)⟩ from random insertions equals the average over a fine grid of the cell (the exact integral)
TEST(Sorption, WidomEqualsTheGridIntegral) {
  Host h = make_host(0.1, 0.24);
  const PairModel pm(h.ff, 5.9, false);
  const double beta = 1 / (0.0019872043 * 300);
  double sum = 0;
  const int G = 60;
  for (int i = 0; i < G; ++i)
    for (int j = 0; j < G; ++j)
      for (int k = 0; k < G; ++k) {
        const Vec3 r{12.0 * (i + 0.5) / G, 12.0 * (j + 0.5) / G, 12.0 * (k + 0.5) / G};
        double e = 0;
        for (uint32_t b = 0; b < 27; ++b) {
          Vec3 d = h.s.atoms[b].pos - r;
          for (auto& x : d) x -= 12.0 * std::round(x / 12.0);
          e += pm.energy(27, b, dot(d, d));
        }
        sum += std::isfinite(e) ? std::exp(-beta * e) : 0.0;
      }
  const double exact = sum / (G * G * G);
  SorptionOptions o;
  o.template_first_atom = 27;
  o.insertions = 400000;
  o.cutoff = 5.9;
  o.coulomb = false;
  const auto r = sorption(h.s, h.ff, o);
  std::printf("widom %.5f ± %.5f, grid %.5f\n", r.widom_w, r.widom_error, exact);
  EXPECT_NEAR(r.widom_w, exact, 0.03 * exact) << "error " << r.widom_error;
  EXPECT_NEAR(r.mu_ex, -std::log(r.widom_w) / beta, 1e-9);
  EXPECT_NEAR(r.solubility, r.widom_w * 273.15 / 300, 1e-12);
}

// GCMC with no interactions at all is an ideal gas: ⟨N⟩ = βpV exactly; with the host, the low-pressure uptake is the
// Henry-law value p V ⟨W⟩ / kT from Widom
TEST(Sorption, GcmcIdealGasAndHenryLimit) {
  Host ideal = make_host(0.0, 0.0);
  SorptionOptions o;
  o.template_first_atom = 27;
  o.insertions = 0;
  o.coulomb = false;
  o.cutoff = 5.9;
  const double V = 1728e-30, kT = 1.380649e-23 * 300;
  const double p = 20 * kT / V / 1000;   // kPa for ⟨N⟩ = 20
  o.pressures_kpa = {p};
  o.steps = 200000;
  const auto r = sorption(ideal.s, ideal.ff, o);
  ASSERT_EQ(r.isotherm.size(), 1u);
  EXPECT_NEAR(r.isotherm[0].loading, 20.0, 0.6);
  Host h = make_host(0.1, 0.24);
  SorptionOptions w = o;
  w.insertions = 200000;
  w.pressures_kpa = {p / 10};   // about two molecules: the Henry regime
  w.steps = 300000;
  const auto rh = sorption(h.s, h.ff, w);
  const double henry = p / 10 * 1000 * V / kT * rh.widom_w;
  std::printf("ideal N %.3f (20) · henry N %.3f vs %.3f\n", r.isotherm[0].loading, rh.isotherm[0].loading, henry);
  EXPECT_NEAR(rh.isotherm[0].loading, henry, 0.08 * henry) << "Widom W " << rh.widom_w;
}

// The density map: for the ideal gas the sorbate fills the cell evenly — the map summed over the cell (density × voxel
// volume) is the mean loading, and every voxel holds about ⟨N⟩/V.
TEST(Sorption, DensityMapOfAnIdealGas) {
  Host ideal = make_host(0.0, 0.0);
  SorptionOptions o;
  o.template_first_atom = 27;
  o.insertions = 0;
  o.coulomb = false;
  o.cutoff = 5.9;
  o.map_grid = 6;
  const double V = 1728e-30, kT = 1.380649e-23 * 300;
  o.pressures_kpa = {20 * kT / V / 1000};
  o.steps = 200000;
  const auto r = sorption(ideal.s, ideal.ff, o);
  const auto& pt = r.isotherm[0];
  ASSERT_EQ(pt.grid, 6);
  ASSERT_EQ(pt.density.size(), 216u);
  const double vox = 1728.0 / 216;
  double tot = 0, mx = 0;
  for (float d : pt.density) tot += d * vox, mx = std::max(mx, double(d));
  EXPECT_NEAR(tot, pt.loading, 0.05 * pt.loading);
  EXPECT_LT(mx, 2.0 * pt.loading / 1728.0);   // no corner fills up: an even gas
}

// Mixtures: in an ideal gas each species takes y_i βpV exactly; in a host at low loading each follows its own Henry
// law, N_i = y_i p V ⟨W_i⟩ / kT, so the adsorption selectivity is the ratio of the Widom factors
TEST(Sorption, BinaryMixtureIdealAndHenry) {
  auto two = [](double eps_host, double eps_a, double eps_b) {
    Host h = make_host(eps_host, eps_a);
    Atom g;
    g.element = 18;
    g.pos = {6, 6, 6};
    h.s.atoms.push_back(g);
    h.ff.type_index.push_back(2);
    h.ff.lj.push_back({eps_b, 3.4});
    h.ff.type_names.push_back("Ar2");
    h.ff.charge.push_back(0.0);
    h.ff.mass.push_back(39.948);
    return h;
  };
  const double V = 1728e-30, kT = 1.380649e-23 * 300;
  const double p = 20 * kT / V / 1000;   // kPa for ⟨N⟩ = 20 in total
  SorptionOptions o;
  o.template_first_atom = 27;
  o.species_first_atom = {27, 28};
  o.mole_fractions = {0.3, 0.7};
  o.insertions = 0;
  o.coulomb = false;
  o.cutoff = 5.9;
  o.pressures_kpa = {p};
  o.steps = 300000;
  const auto ideal = sorption(two(0, 0, 0).s, two(0, 0, 0).ff, o);
  const auto& q = ideal.isotherm.at(0);
  ASSERT_EQ(q.species_loading.size(), 2u);
  std::printf("ideal mixture N %.2f + %.2f (6 + 14), total %.2f\n", q.species_loading[0], q.species_loading[1], q.loading);
  EXPECT_NEAR(q.species_loading[0], 6.0, 0.4);
  EXPECT_NEAR(q.species_loading[1], 14.0, 0.6);
  EXPECT_NEAR(q.loading, q.species_loading[0] + q.species_loading[1], 1e-9);
  EXPECT_NEAR(q.selectivity[1], 1.0, 0.1);

  Host h = two(0.1, 0.24, 0.12);
  SorptionOptions w = o;
  w.mole_fractions = {0.5, 0.5};
  w.insertions = 200000;
  w.pressures_kpa = {p};   // under one molecule in all: the Henry regime
  w.steps = 600000;
  const auto r = sorption(h.s, h.ff, w);
  ASSERT_EQ(r.species.size(), 2u);
  const auto& pt = r.isotherm.at(0);
  for (int i = 0; i < 2; ++i) {
    const double henry = 0.5 * p * 1000 * V / kT * r.species[size_t(i)].widom_w;
    std::printf("species %d: N %.3f ± %.3f vs Henry %.3f\n", i, pt.species_loading[size_t(i)], pt.species_error[size_t(i)], henry);
    EXPECT_NEAR(pt.species_loading[size_t(i)], henry, 0.08 * henry);   // ~5 % high: the sorbates attract each other
  }
  const double s_widom = r.species[1].widom_w / r.species[0].widom_w;
  std::printf("selectivity %.3f vs Widom ratio %.3f\n", pt.selectivity[1], s_widom);
  EXPECT_NEAR(pt.selectivity[1], s_widom, 0.08 * s_widom);
  EXPECT_DOUBLE_EQ(r.widom_w, r.species[0].widom_w);
}

// TraPPE CO2 (Potoff & Siepmann, AIChE J. 47, 1676 (2001)): rigid, linear, C=O 1.16 Å; C σ 2.80 Å ε/k 27.0 K q +0.70,
// O σ 3.05 Å ε/k 79.0 K q −0.35, Lorentz–Berthelot. Molecules m = 0, 1, … along x at the origin, atoms O C O.
namespace {
struct Co2 {
  System s;
  ForceField ff;
};
Co2 trappe_co2(int molecules, double box) {
  Co2 c;
  if (box > 0) c.s.cell.a = {box, 0, 0}, c.s.cell.b = {0, box, 0}, c.s.cell.c = {0, 0, box};
  const double kB = 0.0019872036;
  c.ff.lj = {{27.0 * kB, 2.80}, {79.0 * kB, 3.05}};
  c.ff.type_names = {"CO2C", "CO2O"};
  c.ff.mixing = "arithmetic";
  c.ff.name = "TraPPE CO2";
  for (int m = 0; m < molecules; ++m)
    for (int k = 0; k < 3; ++k) {
      Atom a;
      a.element = k == 1 ? 6 : 8;
      a.pos = {(k - 1) * 1.16, 0, 0};
      a.mol = m + 1;
      c.s.atoms.push_back(a);
      c.ff.type_index.push_back(k == 1 ? 0 : 1);
      c.ff.charge.push_back(k == 1 ? 0.70 : -0.35);
      c.ff.mass.push_back(k == 1 ? 12.011 : 15.9994);
    }
  return c;
}
// B2(T) = −2π N_A ∫ ⟨e^(−βu) − 1⟩_Ω r² dr over the centre–centre distance, orientations uniform (cm³/mol)
double b2_co2(const PairModel& pm, double T, double rmax, int nr, int norient, uint64_t seed) {
  std::mt19937_64 rng(seed);
  std::normal_distribution<double> g(0, 1);
  const double beta = 1 / (0.0019872043 * T);
  auto axis = [&] { Vec3 v{g(rng), g(rng), g(rng)}; return v * (1 / norm(v)); };
  double integral = 0;
  const double dr = rmax / nr;
  for (int i = 0; i < nr; ++i) {
    const double r = (i + 0.5) * dr;
    double f = 0;
    for (int k = 0; k < norient; ++k) {
      const Vec3 a = axis(), b = axis();
      const Vec3 pa[3] = {a * -1.16, {0, 0, 0}, a * 1.16}, pb[3] = {Vec3{r, 0, 0} - b * 1.16, {r, 0, 0}, Vec3{r, 0, 0} + b * 1.16};
      double u = 0;
      for (int p = 0; p < 3; ++p)
        for (int q = 0; q < 3; ++q) {
          const Vec3 d = pb[q] - pa[p];
          u += pm.energy(uint32_t(p), uint32_t(3 + q), dot(d, d));
        }
      f += std::isfinite(u) ? std::exp(-beta * u) - 1 : -1;
    }
    integral += f / norient * r * r * dr;
  }
  return -2 * 3.14159265358979323846 * 6.02214076e23 * integral * 1e-24;   // Å³ → cm³
}
}  // namespace

// The model's second virial coefficient at 300 K against experiment (−121 cm³/mol, Span & Wagner's equation of state),
// and GCMC of the pure gas at 10 bar against the virial expansion of the same (truncated) model: βf = ρ e^(2 B2 ρ)
TEST(Sorption, TrappeCo2VirialAndGcmc) {
  Co2 two = trappe_co2(2, 0);
  const PairModel full(two.ff, 40.0, true, 0.0);   // the whole potential (shifted at 40 Å: negligible for these multipoles)
  const double b2 = b2_co2(full, 300, 30.0, 300, 3000, 7);
  std::printf("TraPPE CO2 B2(300 K) = %.1f cm3/mol (experiment −121)\n", b2);
  EXPECT_NEAR(b2, -121.0, 0.15 * 121.0);
  // NIST's benchmark simulations of the same model (1000 rigid molecules, LJ 15 Å + tail, PPPM), 300 K: 0.1 mol/L at
  // 2.435 ± 0.0003 atm, 0.5 mol/L at 11.65 ± 0.0074 atm; (Z − 1)/ρ = B2 + B3 ρ
  const double RT = 0.0820573661 * 300;   // L atm / mol
  const double nist01 = (2.435 / (0.1 * RT) - 1) / 1e-4, nist05 = (11.65 / (0.5 * RT) - 1) / 5e-4;
  std::printf("NIST (Z − 1)/ρ: %.1f at 0.1 mol/L, %.1f at 0.5 mol/L\n", nist01, nist05);
  EXPECT_NEAR(b2, nist01, 4.0);
  {   // GCMC at the fugacity of NIST's 0.5 mol/L state: βμ_ex = 2 B2 ρ + 1.5 B3 ρ², B3 from the two NIST points
    const double b3 = (nist05 - nist01) / 4e-4, rho = 5e-4;
    const double f_kpa = 0.5 * RT * 101.325 * std::exp(2 * b2 * rho + 1.5 * b3 * rho * rho);
    Co2 g = trappe_co2(1, 60.0);
    SorptionOptions so;
    so.template_first_atom = 0, so.insertions = 0, so.cutoff = 12.0, so.temperature = 300, so.pressures_kpa = {f_kpa}, so.steps = 300000;
    const auto r5 = sorption(g.s, g.ff, so);
    const double got = r5.isotherm.at(0).loading / 216000.0 / 6.02214076e23 * 1e27;
    std::printf("GCMC at f = %.1f kPa: %.4f ± %.4f mol/L (NIST 0.5)\n", f_kpa, got, r5.isotherm[0].loading_error / 216000.0 / 6.02214076e23 * 1e27);
    EXPECT_NEAR(got, 0.5, 0.02);
  }
  // GCMC in an empty 60 Å box at f = 10 bar, 300 K, the sorption engine's own model (12 Å, DSF α 0.2)
  Co2 gas = trappe_co2(1, 60.0);
  const PairModel trunc(two.ff, 12.0, true, 0.2);
  const double b2t = b2_co2(trunc, 300, 12.0, 240, 3000, 9);
  SorptionOptions o;
  o.template_first_atom = 0;
  o.insertions = 0;
  o.cutoff = 12.0;
  o.pressures_kpa = {1000};
  o.steps = 400000;
  o.temperature = 300;
  const auto r = sorption(gas.s, gas.ff, o);
  const double V = 216000e-24;                                // cm³
  const double bf = 1000e3 / (1.380649e-23 * 300) * 1e-6;    // βf, molecules/cm³
  double rho = bf;
  for (int it = 0; it < 50; ++it) rho = bf * std::exp(-2 * b2t / 6.02214076e23 * rho);
  const double expect = rho * V, ideal = bf * V;
  std::printf("GCMC N %.2f ± %.2f · virial %.2f (B2 %.1f cm3/mol with the 12 Å model) · ideal gas %.2f\n", r.isotherm[0].loading,
              r.isotherm[0].loading_error, expect, b2t, ideal);
  EXPECT_NEAR(r.isotherm[0].loading, expect, std::max(3 * r.isotherm[0].loading_error, 0.02 * expect));
  EXPECT_GT(r.isotherm[0].loading, ideal);   // attraction: more than the ideal gas
}

// TraPPE N2 against NIST's benchmark simulations of the same model (1000 rigid molecules, LJ 15 Å + tail, PPPM): the gas
// at 110 K and 1.0 mol/L has p = 7.925 ± 0.006 atm, Z = 0.8780 ± 0.0007, so B2 + B3 ρ = (Z − 1)/ρ = −122.0 cm³/mol;
// B3 ρ is a few cm³/mol here. The model: σ 3.31 Å, ε/k 36.0 K on each N (−0.482 e), +0.964 e at the centre, N–N 1.10 Å.
TEST(Sorption, TrappeN2AgainstNistBenchmark) {
  ForceField ff;
  const double kB = 0.0019872036;
  ff.lj = {{36.0 * kB, 3.31}, {0, 0}};
  ff.type_names = {"N2N", "N2M"};
  ff.mixing = "arithmetic";
  const double off[3] = {-0.55, 0.0, 0.55};
  for (int m = 0; m < 2; ++m)
    for (int k = 0; k < 3; ++k) {
      ff.type_index.push_back(k == 1 ? 1 : 0);
      ff.charge.push_back(k == 1 ? 0.964 : -0.482);
      ff.mass.push_back(k == 1 ? 0 : 14.007);
    }
  const PairModel pm(ff, 40.0, true, 0.0);
  std::mt19937_64 rng(11);
  std::normal_distribution<double> g(0, 1);
  auto axis = [&] { Vec3 v{g(rng), g(rng), g(rng)}; return v * (1 / norm(v)); };
  const double T = 110, beta = 1 / (0.0019872043 * T), rmax = 30, nr = 300, dr = rmax / nr;
  double integral = 0;
  for (int i = 0; i < nr; ++i) {
    const double r = (i + 0.5) * dr;
    double f = 0;
    const int no = 4000;
    for (int k = 0; k < no; ++k) {
      const Vec3 a = axis(), b = axis();
      double u = 0;
      for (int p = 0; p < 3; ++p)
        for (int q = 0; q < 3; ++q) {
          const Vec3 d = Vec3{r, 0, 0} + b * off[q] - a * off[p];
          u += pm.energy(uint32_t(p), uint32_t(3 + q), dot(d, d));
        }
      f += std::isfinite(u) ? std::exp(-beta * u) - 1 : -1;
    }
    integral += f / no * r * r * dr;
  }
  const double b2 = -2 * 3.14159265358979323846 * 6.02214076e23 * integral * 1e-24;
  const double p_ideal = 1000 * 8.314462618 * T / 101325;   // atm at 1.0 mol/L
  const double nist = (7.925 / p_ideal - 1) / 1e-3;           // cm³/mol
  std::printf("TraPPE N2 B2(110 K) = %.1f cm3/mol · NIST (Z − 1)/ρ at 1 mol/L = %.1f (B2 + B3 ρ)\n", b2, nist);
  EXPECT_NEAR(b2, nist, 6.0);
  // GCMC of the gas at the fugacity of NIST's state, f = ρRT exp(βμ_ex), βμ_ex = 2 B2 ρ + 1.5 B3 ρ² with B3 ρ = nist − b2, in a 60 Å box:
  // the density comes back at 1.0 mol/L (the engine's 12 Å, DSF, no tail: within a few per cent)
  const double rho = 1e-3, b3rho = nist - b2;
  const double f_kpa = p_ideal * 101.325 * std::exp(2 * b2 * rho + 1.5 * b3rho * rho);
  System s;
  s.cell.a = {60, 0, 0}, s.cell.b = {0, 60, 0}, s.cell.c = {0, 0, 60};
  ForceField one = ff;
  one.type_index.resize(3), one.charge.resize(3), one.mass.resize(3);
  for (int k = 0; k < 3; ++k) {
    Atom a;
    a.element = k == 1 ? 0 : 7;
    a.pos = {off[k], 0, 0};
    s.atoms.push_back(a);
  }
  SorptionOptions o;
  o.template_first_atom = 0;
  o.insertions = 0;
  o.temperature = T;
  o.cutoff = 12.0;
  o.pressures_kpa = {f_kpa};
  o.steps = 300000;
  const auto r = sorption(s, one, o);
  const double got = r.isotherm.at(0).loading / 216000.0 / 6.02214076e23 * 1e27;   // mol/L
  std::printf("GCMC at f = %.1f kPa: %.4f ± %.4f mol/L (NIST 1.0)\n", f_kpa, got, r.isotherm[0].loading_error / 216000.0 / 6.02214076e23 * 1e27);
  EXPECT_NEAR(got, 1.0, 0.03);
}
