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
  // a cell whose c leans out of z is refused (in-plane obliqueness is fine: see ObliqueAndKeepVolume)
  System bad = film;
  bad.cell.c = {3, 0, 21};
  EXPECT_THROW(stack_layers({{"argon", &base}, {"bad", &bad}}, so), std::invalid_argument);
}

// Oblique cells: a hexagonal sheet (a = 2.46 Å, γ = 120°) under a 2.50 Å hexagonal sheet turned by 30° in its file: the
// stack is hexagonal on the first layer's cell (γ 120°, a along x), the second strained by 2.46/2.50 − 1 along both axes,
// no shear, every atom in the cell; with keep_volume the second layer's thickness grows by 1/(0.984)² so its area ×
// thickness is kept, while the first (unstrained) keeps its own
TEST(Layers, ObliqueAndKeepVolume) {
  auto sheet = [](double a, double turn_deg, int z_layers) {
    System s;
    const double t = turn_deg * M_PI / 180;
    auto rot = [&](Vec3 v) { return Vec3{std::cos(t) * v[0] - std::sin(t) * v[1], std::sin(t) * v[0] + std::cos(t) * v[1], v[2]}; };
    s.cell.a = rot({a, 0, 0}), s.cell.b = rot({-a / 2, a * std::sqrt(3.0) / 2, 0}), s.cell.c = {0, 0, 20};
    for (int k = 0; k < z_layers; ++k)
      for (Vec3 f : {Vec3{0, 0, 0}, Vec3{1.0 / 3, 2.0 / 3, 0}}) {
        Atom at;
        at.element = 6;
        at.pos = s.cell.a * f[0] + s.cell.b * f[1] + Vec3{0, 0, 2.0 + 3.35 * k};
        s.atoms.push_back(at);
      }
    return s;
  };
  const System s1 = sheet(2.46, 0, 2), s2 = sheet(2.50, 30, 3);
  StackOptions so;
  so.match = "first", so.max_repeat = 1, so.vacuum = 10;
  StackReport rep;
  const System t = stack_layers({{"A", &s1}, {"B", &s2}}, so, &rep);
  EXPECT_NEAR(t.cell.a[0], 2.46, 1e-12);
  EXPECT_NEAR(t.cell.a[1], 0, 1e-12);
  EXPECT_NEAR(std::acos(dot(t.cell.a, t.cell.b) / (norm(t.cell.a) * norm(t.cell.b))) * 180 / M_PI, 120, 1e-9);
  EXPECT_NEAR(rep.layers[1].strain_a, 2.46 / 2.50 - 1, 1e-12);
  EXPECT_NEAR(rep.layers[1].strain_b, 2.46 / 2.50 - 1, 1e-12);
  EXPECT_NEAR(rep.layers[1].shear_deg, 0, 1e-9);
  for (const auto& a : t.atoms) {
    const double det = t.cell.a[0] * t.cell.b[1] - t.cell.a[1] * t.cell.b[0];
    const double u = (a.pos[0] * t.cell.b[1] - a.pos[1] * t.cell.b[0]) / det, v = (t.cell.a[0] * a.pos[1] - t.cell.a[1] * a.pos[0]) / det;
    EXPECT_GE(u, -1e-12);
    EXPECT_LT(u, 1);
    EXPECT_GE(v, -1e-12);
    EXPECT_LT(v, 1);
  }
  const double thick = 2 * 3.35;
  EXPECT_NEAR(rep.layers[1].z_hi - rep.layers[1].z_lo, thick, 1e-12);
  so.keep_volume = true;
  const System k = stack_layers({{"A", &s1}, {"B", &s2}}, so, &rep);
  const double f = 2.46 / 2.50;
  EXPECT_NEAR(rep.layers[1].z_hi - rep.layers[1].z_lo, thick / (f * f), 1e-12);
  EXPECT_NEAR((rep.layers[1].z_hi - rep.layers[1].z_lo) * f * f, thick, 1e-12);   // area × thickness kept
  EXPECT_NEAR(rep.layers[0].z_hi - rep.layers[0].z_lo, 3.35, 1e-12);
  EXPECT_EQ(k.atoms.size(), 10u);
}

// Layer operations: the lattice mismatch shared (match "average": 20 and 21 Å cells meet at 20.5, ±2.5 % each), a layer
// flipped (its hydrogen face from the top to the bottom), and one shifted in the plane by a set amount.
TEST(Layers, AverageFlipAndShift) {
  System a = slab();
  System b = slab();
  b.cell.a = {21, 0, 0}, b.cell.b = {0, 21, 0};
  for (auto& at : b.atoms) at.pos[0] *= 21.0 / 20, at.pos[1] *= 21.0 / 20;
  StackOptions o;
  o.match = "average";
  o.max_repeat = 1;
  StackReport r;
  stack_layers({{"A", &a}, {"B", &b}}, o, &r);
  EXPECT_NEAR(r.a, 20.5, 1e-9);
  EXPECT_NEAR(r.layers[0].strain_a, 20.5 / 20 - 1, 1e-9);
  EXPECT_NEAR(r.layers[1].strain_a, 20.5 / 21 - 1, 1e-9);
  // a top face of hydrogens on the second layer; flipped, they end at its bottom
  System h = slab();
  for (int i = 0; i < 5; ++i) { Atom x; x.element = 1; x.pos = {4.0 * i + 1, 1, 4.0 * 4 + 3 + 1.0}; h.atoms.push_back(x); }
  StackOptions p;
  StackLayerInput top{"H on top", &h};
  top.flip = true;
  top.shift_x = 1.5;
  StackReport rf;
  const System st = stack_layers({{"base", &a}, top}, p, &rf);
  double hz = 1e300, ar_min = 1e300;
  for (size_t i = a.atoms.size(); i < st.atoms.size(); ++i) {
    if (st.atoms[i].element == 1) hz = std::min(hz, st.atoms[i].pos[2]);
    else ar_min = std::min(ar_min, st.atoms[i].pos[2]);
  }
  EXPECT_LT(hz, ar_min);   // the hydrogens below the flipped layer's argon
  EXPECT_NE(rf.layers[1].name.find("flipped"), std::string::npos);
  // the shift: the flipped layer's argon x coordinates are the base lattice's plus 1.5 Å (mod the cell)
  for (size_t i = a.atoms.size(); i < st.atoms.size(); ++i)
    if (st.atoms[i].element == 18) EXPECT_NEAR(std::fmod(st.atoms[i].pos[0] - 1 - 1.5 + 40, 4.0), 0.0, 1e-9);
}
