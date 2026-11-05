#include <gtest/gtest.h>

#include <cmath>

#include "caps/crystal.hpp"
#include "caps/nano.hpp"

using namespace caps;

namespace {
// a cubic 5 Å cell holding one C–C dimer (1.5 Å, along x) that the cell's edge cuts: one atom at x = 4.25, the other
// wrapped to x = 0.75 (bonded through the boundary)
System dimer_crystal() {
  System s;
  s.cell.a = {5, 0, 0}, s.cell.b = {0, 5, 0}, s.cell.c = {0, 0, 5};
  Atom a;
  a.element = 6, a.pos = {4.25, 2.5, 2.5};
  Atom b = a;
  b.pos = {0.75, 2.5, 2.5};
  s.atoms = {a, b};
  s.bonds = {{0, 1, 1}};
  return s;
}
// every carbon of the result bonded to exactly one other carbon 1.5 Å away
int partnered(const System& s) {
  int n = 0;
  for (size_t i = 0; i < s.atoms.size(); ++i) {
    int k = 0;
    for (size_t j = 0; j < s.atoms.size(); ++j)
      if (i != j && std::fabs(norm(s.atoms[j].pos - s.atoms[i].pos) - 1.5) < 1e-6) ++k;
    n += k == 1;
  }
  return n;
}
}  // namespace

// Whole molecules: every dimer of a sphere cut is complete, and the dimers kept are exactly those whose centres lie
// within the radius
TEST(MolecularCuts, SphereKeepsWholeDimers) {
  const System x = dimer_crystal();
  ParticleOptions o;
  o.shape = ParticleShape::Sphere;
  o.radius = 9.0;
  o.on_atom = false;   // centred on the cell centre (2.5, 2.5, 2.5)
  o.whole_molecules = true;
  NanoReport rep;
  const System s = nanoparticle(x, o, &rep);
  EXPECT_EQ(partnered(s), int(s.atoms.size()));
  int expect = 0;   // whole dimers centred at (5 + 5i, 2.5 + 5j, 2.5 + 5k), the cut's centre (2.5, 2.5, 2.5): within 9 Å
  for (int i = -3; i <= 3; ++i)
    for (int j = -3; j <= 3; ++j)
      for (int k = -3; k <= 3; ++k) expect += std::pow(2.5 + 5 * i, 2) + 25.0 * (j * j + k * k) <= 81.0;
  EXPECT_EQ(int(s.atoms.size()), 2 * expect);
}

// (100) as the base: the dimers, along x in the crystal, lie along z in the particle
TEST(MolecularCuts, BaseNormalToPlanes) {
  ParticleOptions o;
  o.shape = ParticleShape::Cube;
  o.radius = 10;
  o.on_atom = false;
  o.whole_molecules = true;
  o.hkl = {1, 0, 0};
  const System s = nanoparticle(dimer_crystal(), o);
  ASSERT_FALSE(s.atoms.empty());
  for (size_t i = 0; i < s.atoms.size(); ++i)
    for (size_t j = i + 1; j < s.atoms.size(); ++j) {
      const Vec3 d = s.atoms[j].pos - s.atoms[i].pos;
      if (std::fabs(norm(d) - 1.5) < 1e-6) EXPECT_NEAR(std::fabs(d[2]), 1.5, 1e-9);
    }
}

// A polymer crystal's chains are periodic networks: kept atom by atom, with a note
TEST(MolecularCuts, NetworksAreCutAtomByAtom) {
  const System pe = read_cif(std::string(CAPS_SOURCE_DIR) + "/data/crystals/polymers/pe-orthorhombic.cif");
  ParticleOptions o;
  o.shape = ParticleShape::Sphere;
  o.radius = 8;
  o.whole_molecules = true;
  NanoReport rep;
  const System s = nanoparticle(pe, o, &rep);
  EXPECT_FALSE(s.atoms.empty());
  bool noted = false;
  for (const auto& n : rep.notes) noted = noted || n.find("periodic network") != std::string::npos;
  EXPECT_TRUE(noted);
}
