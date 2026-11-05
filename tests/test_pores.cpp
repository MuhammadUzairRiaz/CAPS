#include <gtest/gtest.h>

#include <cmath>

#include "caps/elements.hpp"
#include "caps/pores.hpp"

using namespace caps;

namespace {
// one carbon (radius set by hand) at (0.1, 0.1, 0.1) in a cubic cell of edge L: a simple cubic lattice of spheres
System lattice(double L) {
  System s;
  s.cell.a = {L, 0, 0}, s.cell.b = {0, L, 0}, s.cell.c = {0, 0, L};
  Atom a;
  a.element = 6, a.pos = {0.1, 0.1, 0.1};
  s.atoms = {a};
  return s;
}
}  // namespace

// A simple cubic lattice of spheres (L = 20 Å, R = 2 Å, probe 1 Å): the largest included sphere sits at the body centre,
// LCD = 2(√3 L/2 − R); the windows are the face centres between four spheres, PLD = 2(L/√2 − R) along every axis; the
// accessible surface is the whole probe sphere, 4π(R + r)², per volume and per mass by the definitions
TEST(Pores, SimpleCubicLattice) {
  const double L = 20, R = 2, r = 1;
  PoreNetworkOptions o;
  o.probe = r, o.grid = 0.2, o.radii = {R}, o.samples = 4000;
  const PoreNetworkReport p = pore_network(lattice(L), o);
  EXPECT_NEAR(p.lcd, 2 * (std::sqrt(3.0) * L / 2 - R), 1e-5);
  for (int a = 0; a < 3; ++a) EXPECT_NEAR(p.pld_axis[size_t(a)], 2 * (L / std::sqrt(2.0) - R), 1e-5) << a;
  EXPECT_NEAR(p.pld, 2 * (L / std::sqrt(2.0) - R), 1e-5);
  EXPECT_NEAR(p.dif, p.lcd, 1e-5);   // the free path runs through the body centre's cavity
  EXPECT_EQ(p.channel_dims, 3);
  ASSERT_EQ(p.clusters.size(), 1u);
  EXPECT_TRUE(p.clusters[0].channel);
  const double area = 4 * M_PI * (R + r) * (R + r);
  EXPECT_NEAR(p.asa, area, 1e-9 * area);
  EXPECT_EQ(p.nasa, 0.0);
  EXPECT_NEAR(p.asa_volume, area / (L * L * L) * 1e4, 1e-9);
  const double m = element(6).mass;
  EXPECT_NEAR(p.asa_mass, area * 1e-20 * 6.02214076e23 / m, 1e-6);
  EXPECT_NEAR(p.density, m / 6.02214076e23 / (L * L * L * 1e-24), 1e-12);
  // the probe-centre space: the cell less a sphere of R + r (to the grid's resolution)
  EXPECT_NEAR(p.av_fraction, 1 - 4.0 / 3 * M_PI * std::pow(R + r, 3) / (L * L * L), 2e-3);
}

// Spheres that close the windows (R = 15 Å > L/√2) leave a cavity at the body centre: one pocket, no channel, PLD 0, all
// surface non-accessible
TEST(Pores, ClosedPocket) {
  PoreNetworkOptions o;
  o.probe = 0.5, o.grid = 0.2, o.radii = {15.0}, o.samples = 500;
  const PoreNetworkReport p = pore_network(lattice(20), o);
  EXPECT_EQ(p.channel_dims, 0);
  EXPECT_EQ(p.pld, 0.0);
  ASSERT_EQ(p.clusters.size(), 1u);
  EXPECT_FALSE(p.clusters[0].channel);
  EXPECT_EQ(p.av_fraction, 0.0);
  EXPECT_GT(p.nav_fraction, 0.0);
  EXPECT_EQ(p.asa, 0.0);
  EXPECT_GT(p.nasa, 0.0);
  EXPECT_NEAR(p.lcd, 2 * (std::sqrt(3.0) * 10 - 15), 1e-5);
}

// Two overlapping probe spheres (R + r = 3 Å, centres 4 Å apart, alone in a large cell): the union's area, each sphere
// less a cap of height R + r − d/2
TEST(Pores, OverlappingSpheres) {
  System s;
  s.cell.a = {30, 0, 0}, s.cell.b = {0, 30, 0}, s.cell.c = {0, 0, 30};
  Atom a;
  a.element = 6, a.pos = {13, 15, 15};
  Atom b = a;
  b.pos = {17, 15, 15};
  s.atoms = {a, b};
  PoreNetworkOptions o;
  o.probe = 1, o.radii = {2, 2}, o.grid = 0.5, o.samples = 20000;
  const PoreNetworkReport p = pore_network(s, o);
  const double Rp = 3, h = Rp - 2;
  EXPECT_NEAR(p.asa, 2 * (4 * M_PI * Rp * Rp - 2 * M_PI * Rp * h), 2e-3 * 2 * 4 * M_PI * Rp * Rp);
}
