#include <gtest/gtest.h>

#include <cmath>
#include <cstdio>
#include <limits>

#include "caps/adsorption.hpp"

using namespace caps;

namespace {
// a fixed square lattice of LJ sites (2.5 Å apart, z = 10) in a 20 × 20 × 30 Å periodic cell, and `extra` adsorbate atoms
struct Slab {
  System s;
  ForceField ff;
};
Slab make_slab(int ads, bool diatomic) {
  Slab r;
  r.s.cell.a = {20, 0, 0}, r.s.cell.b = {0, 20, 0}, r.s.cell.c = {0, 0, 30};
  for (int i = 0; i < 8; ++i)
    for (int j = 0; j < 8; ++j) {
      Atom a;
      a.element = 6;
      a.pos = {2.5 * i, 2.5 * j, 10.0};
      r.s.atoms.push_back(a);
      r.ff.type_index.push_back(0);
    }
  for (int k = 0; k < ads; ++k) {
    Atom a;
    a.element = diatomic ? 7 : 18;
    a.pos = {3.0 + 4 * k, 5.0, 20.0};
    r.s.atoms.push_back(a);
    r.ff.type_index.push_back(1);
    if (diatomic) {
      Atom b = a;
      b.pos = {a.pos[0] + 1.1, 5.0, 20.0};
      r.s.atoms.push_back(b);
      r.ff.type_index.push_back(1);
      r.s.bonds.push_back({uint32_t(r.s.atoms.size() - 2), uint32_t(r.s.atoms.size() - 1), 3});
    }
  }
  r.ff.lj = {{0.07, 3.4}, {diatomic ? 0.072 : 0.238, diatomic ? 3.3 : 3.4}};
  r.ff.type_names = {"C", diatomic ? "N" : "Ar"};
  r.ff.charge.assign(r.s.atoms.size(), 0.0);
  r.ff.mass.assign(r.s.atoms.size(), 12.0);
  r.ff.name = "test LJ";
  return r;
}
}  // namespace

// One argon atom: the locator's lowest energy equals the minimum of a fine scan over the lattice's unit cell and heights
TEST(Adsorption, FindsTheScanMinimum) {
  Slab sl = make_slab(1, false);
  double emin = std::numeric_limits<double>::infinity();
  System t = sl.s;
  for (int i = 0; i <= 20; ++i)
    for (int j = 0; j <= i; ++j)   // the unit cell's symmetric half
      for (int k = 0; k <= 150; ++k) {
        t.atoms.back().pos = {1.25 * i / 10.0, 1.25 * j / 10.0, 12.0 + 3.0 * k / 150};
        emin = std::min(emin, adsorption_interaction(t, sl.ff, adsorbate_molecules(t, 64), 9.0, false));
      }
  AdsorptionOptions o;
  o.coulomb = false;
  o.cutoff = 9.0;
  o.steps = 4000;
  o.cycles = 2;
  o.z_lo = 10.5, o.z_hi = 25;
  o.first_mobile_atom = 64;
  AdsorptionReport r;
  locate_adsorption(sl.s, sl.ff, o, &r);
  EXPECT_LT(emin, -0.5);
  EXPECT_NEAR(r.adsorption_energy, emin, 0.002 * std::fabs(emin)) << "scan " << emin;
  EXPECT_NEAR(adsorption_interaction(sl.s, sl.ff, adsorbate_molecules(sl.s, 64), 9.0, false), r.adsorption_energy, 1e-9);
  ASSERT_EQ(r.components.size(), 1u);
  EXPECT_EQ(r.components[0].name, "Ar");
  EXPECT_NEAR(r.components[0].de_dn, r.adsorption_energy, 1e-9);
  EXPECT_FALSE(r.hist_counts.empty());
}

// Three rigid diatomics: each stays whole, the substrate stays put, and dE/dN sums as it must (each adsorbate–adsorbate
// pair counted for both molecules)
TEST(Adsorption, RigidMoleculesAndDeDn) {
  Slab sl = make_slab(3, true);
  const System before = sl.s;
  AdsorptionOptions o;
  o.coulomb = false;
  o.cutoff = 9.0;
  o.steps = 3000;
  o.cycles = 2;
  o.z_lo = 10.5, o.z_hi = 25;
  o.first_mobile_atom = 64;
  AdsorptionReport r;
  locate_adsorption(sl.s, sl.ff, o, &r);
  for (size_t i = 0; i < 64; ++i) EXPECT_EQ(sl.s.atoms[i].pos, before.atoms[i].pos);
  for (const auto& b : sl.s.bonds) EXPECT_NEAR(norm(sl.s.cell.minimum_image(sl.s.atoms[b.j].pos - sl.s.atoms[b.i].pos)), 1.1, 1e-9);
  ASSERT_EQ(r.components.size(), 1u);
  EXPECT_EQ(r.components[0].molecules, 3);
  EXPECT_NEAR(3 * r.components[0].de_dn, r.adsorbate_substrate + 2 * r.adsorbate_adsorbate, 1e-9);
  EXPECT_LT(r.adsorption_energy, 0);
  EXPECT_GE(r.configs.size(), 1u);
  EXPECT_NEAR(r.configs[0].energy, r.adsorption_energy, 1e-9);
}
