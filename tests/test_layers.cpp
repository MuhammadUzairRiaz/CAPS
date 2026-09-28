#include <gtest/gtest.h>

#include <cmath>

#include "caps/layers.hpp"
#include "caps/polymer.hpp"

using namespace caps;

namespace {
// a square-lattice slab of argon, 5 layers thick, in a 20 × 20 Å cell with vacuum
System slab() {
  System s;
  s.cell.a = {20, 0, 0}, s.cell.b = {0, 20, 0}, s.cell.c = {0, 0, 30};
  for (int i = 0; i < 5; ++i)
    for (int j = 0; j < 5; ++j)
      for (int k = 0; k < 5; ++k) {
        Atom a;
        a.element = 18;
        a.pos = {4.0 * i + 1, 4.0 * j + 1, 4.0 * k + 3};
        s.atoms.push_back(a);
      }
  return s;
}
}  // namespace

TEST(Layers, StackKeepsChainsWholeAndLayersApart) {
  const System base = slab();
  GrowOptions o;
  o.chains = 4;
  o.box = 21;
  o.seed = 3;
  ChainSpec c;
  c.units = {{"A", "*CC*"}};
  c.dp = 20;
  const System film = grow_chains(c, o);
  StackOptions so;
  so.gap = 2.5;
  so.vacuum = 10;
  StackReport rep;
  const System s = stack_layers({{"argon", &base}, {"PE", &film}}, so, &rep);
  ASSERT_EQ(rep.layers.size(), 2u);
  // 20 Å against 21 Å: one repeat each, the film squeezed by 4.8 %
  EXPECT_EQ(rep.layers[1].na, 1);
  EXPECT_NEAR(rep.layers[1].strain_a, 20.0 / 21.0 - 1, 1e-9);
  EXPECT_EQ(s.atoms.size(), base.atoms.size() + film.atoms.size());
  EXPECT_NEAR(rep.a, 20, 1e-9);
  // the film starts gap above the slab's top atoms; vacuum above the film
  double slab_top = -1e9, film_lo = 1e9, film_hi = -1e9;
  for (size_t i = 0; i < s.atoms.size(); ++i) {
    if (i < base.atoms.size()) slab_top = std::max(slab_top, s.atoms[i].pos[2]);
    else film_lo = std::min(film_lo, s.atoms[i].pos[2]), film_hi = std::max(film_hi, s.atoms[i].pos[2]);
  }
  EXPECT_NEAR(film_lo - slab_top, 2.5, 1e-9);
  EXPECT_NEAR(rep.c - film_hi, 10, 1e-9);
  // chains are whole along z: every bond short (x and y by minimum image)
  double longest = 0;
  for (const auto& b : s.bonds) {
    Vec3 d = s.atoms[b.j].pos - s.atoms[b.i].pos;
    d[0] -= rep.a * std::round(d[0] / rep.a);
    d[1] -= rep.b * std::round(d[1] / rep.b);
    longest = std::max(longest, norm(d));
  }
  EXPECT_LT(longest, 1.7);
  // molecules: 125 argon atoms, then the four chains
  int nm = 0;
  s.molecules(&nm);
  EXPECT_EQ(nm, 125 + 4);
  // a small film on a big slab is repeated: 10 Å cells tile 20 Å twice
  System small = film;
  for (auto& a : small.atoms) a.pos = a.pos * (10.0 / 21.0);
  small.cell.a = {10, 0, 0}, small.cell.b = {0, 10, 0}, small.cell.c = {0, 0, 10};
  const System t = stack_layers({{"argon", &base}, {"small", &small}}, so, &rep);
  EXPECT_EQ(rep.layers[1].na, 2);
  EXPECT_EQ(rep.layers[1].nb, 2);
  EXPECT_EQ(t.atoms.size(), base.atoms.size() + 4 * small.atoms.size());
  // a sheared cell is refused
  System bad = film;
  bad.cell.b = {3, 21, 0};
  EXPECT_THROW(stack_layers({{"argon", &base}, {"bad", &bad}}, so), std::invalid_argument);
}
