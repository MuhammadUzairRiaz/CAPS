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

// A chain crystal cut with whole molecules gives the same atoms as a plain cut (its chains are one periodic network):
// no H is left behind by a carbon outside the cube
TEST(MolecularCuts, NetworkCutLeavesNoBareHydrogen) {
  const System pe = read_cif(std::string(CAPS_SOURCE_DIR) + "/data/crystals/polymers/pe-orthorhombic.cif");
  ParticleOptions o;
  o.shape = ParticleShape::Cube;
  o.radius = 10;
  const System plain = nanoparticle(pe, o);
  o.whole_molecules = true;
  const System whole = nanoparticle(pe, o);
  EXPECT_EQ(whole.atoms.size(), plain.atoms.size());
  for (size_t i = 0; i < whole.atoms.size(); ++i) {
    if (whole.atoms[i].element != 1) continue;
    double dmin = 1e300;
    for (size_t j = 0; j < whole.atoms.size(); ++j)
      if (j != i) dmin = std::min(dmin, norm(whole.atoms[j].pos - whole.atoms[i].pos));
    EXPECT_LT(dmin, 1.2);
  }
}

// A sphere cut from a chain crystal keeps no carbon cut off from every other carbon, so passivation makes no methane
TEST(MolecularCuts, NoDetachedCarbonsToPassivate) {
  const System pe = read_cif(std::string(CAPS_SOURCE_DIR) + "/data/crystals/polymers/pe-orthorhombic.cif");
  ParticleOptions o;
  o.shape = ParticleShape::Sphere;
  o.radius = 10;
  o.passivate = true;
  const System s = nanoparticle(pe, o);
  ASSERT_FALSE(s.atoms.empty());
  for (size_t i = 0; i < s.atoms.size(); ++i) {
    if (s.atoms[i].element != 6) continue;
    double dmin = 1e300;
    for (size_t j = 0; j < s.atoms.size(); ++j)
      if (j != i && s.atoms[j].element == 6) dmin = std::min(dmin, norm(s.atoms[j].pos - s.atoms[i].pos));
    EXPECT_LT(dmin, 1.7);
  }
}

// A passivated quartz channel keeps Si, O and H as three types with their own masses
TEST(MolecularCuts, PassivatedChannelKeepsSiliconType) {
  const System quartz = read_cif(std::string(CAPS_SOURCE_DIR) + "/data/crystals/alpha-quartz.cif");
  PoreOptions o;
  o.kind = PoreKind::Cylinder;
  o.crystal = &quartz;
  o.width = 14;
  o.wall = 6;
  o.length = 20;
  o.passivate = true;
  const System s = build_pore(o);
  ASSERT_EQ(s.types.size(), 3u);
  for (const auto& a : s.atoms) {
    const TypeInfo* t = nullptr;
    for (const auto& ti : s.types)
      if (ti.type == a.type) t = &ti;
    ASSERT_NE(t, nullptr);
    if (a.element == 14) EXPECT_NEAR(t->mass, 28.0855, 0.01);
    if (a.element == 1) EXPECT_NEAR(t->mass, 1.008, 0.01);
  }
}
