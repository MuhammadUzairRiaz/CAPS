#include <gtest/gtest.h>

#include <cmath>
#include <cstdio>

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
