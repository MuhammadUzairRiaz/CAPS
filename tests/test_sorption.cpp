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
