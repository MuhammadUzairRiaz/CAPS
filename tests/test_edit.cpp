#include <gtest/gtest.h>

#include <cmath>

#include "caps/appearance.hpp"
#include "caps/edit.hpp"
#include "caps/elements.hpp"
#include "caps/molecule.hpp"
#include "caps/polymer.hpp"

using namespace caps;

namespace {

System molecule(const std::string& smiles) {
  BuildOptions o;
  o.forcefield = "uff";
  return build_molecule(smiles, o).system;
}

double angle_deg(const Vec3& a, const Vec3& b, const Vec3& c) {
  const Vec3 u = a - b, v = c - b;
  return std::acos(dot(u, v) / (norm(u) * norm(v))) * 180 / M_PI;
}

System polystyrene(Tacticity t) {
  ChainSpec spec;
  spec.units.push_back({"styrene", "[*]CC([*])c1ccccc1"});
  spec.dp = 12;
  spec.tacticity = t;
  GrowOptions o;
  o.chains = 1;
  o.density = 0.3;
  o.seed = 3;
  return grow_chains(spec, o);
}

}  // namespace

TEST(Edit, HydrogensFillValencesTetrahedrally) {
  System s;
  for (int z : {6, 6, 8}) { Atom a; a.element = z; s.atoms.push_back(a); }
  s.atoms[1].pos = {1.52, 0, 0};
  s.atoms[2].pos = {2.0, 1.35, 0};
  s.bonds = {{0, 1, 1}, {1, 2, 1}};
  EXPECT_EQ(add_hydrogens(s), 6);   // CH3 CH2 OH
  ASSERT_EQ(s.atoms.size(), 9u);
  int per[3] = {0, 0, 0};
  for (const auto& b : s.bonds)
    for (uint32_t k : {b.i, b.j}) if (k < 3 && s.atoms[b.i + b.j - k].element == 1) ++per[k];
  EXPECT_EQ(per[0], 3);
  EXPECT_EQ(per[1], 2);
  EXPECT_EQ(per[2], 1);
  // every angle at the methyl carbon near tetrahedral
  std::vector<uint32_t> around;
  for (const auto& b : s.bonds) { if (b.i == 0) around.push_back(b.j); if (b.j == 0) around.push_back(b.i); }
  for (size_t a = 0; a < around.size(); ++a)
    for (size_t b = a + 1; b < around.size(); ++b)
      EXPECT_NEAR(angle_deg(s.atoms[around[a]].pos, s.atoms[0].pos, s.atoms[around[b]].pos), 109.47, 3.0);
  EXPECT_EQ(add_hydrogens(s), 0);   // nothing left to fill
  EXPECT_EQ(default_valence(7, 1), 4);
}

TEST(Edit, AtomsBondsElementsAndDeletion) {
  System s = molecule("C");
  const size_t n0 = s.atoms.size();   // CH4
  // replace one H by a carbon: ethane's skeleton, then fill it
  uint32_t h = 1;
  set_element(s, h, 6);
  EXPECT_EQ(s.atoms[h].element, 6);
  EXPECT_EQ(add_hydrogens(s), 3);
  EXPECT_EQ(s.atoms.size(), n0 + 3);
  const uint32_t o = add_atom(s, int(h), 8, 1);
  EXPECT_NEAR(norm(s.atoms[o].pos - s.atoms[h].pos), element(8).covalent + element(6).covalent, 1e-9);
  add_bond(s, 0, o, 1);
  EXPECT_TRUE(remove_bond(s, 0, o));
  std::vector<char> del(s.atoms.size(), 0);
  del[o] = 1;
  const size_t nb = s.bonds.size();
  delete_atoms(s, del);
  EXPECT_EQ(s.atoms.size(), n0 + 3);
  EXPECT_EQ(s.bonds.size(), nb - 1);
}

TEST(Edit, InvertingTheCentreOfAlanine) {
  System ala = molecule("N[C@@H](C)C(=O)O");   // L = S
  uint32_t ca = 1;
  auto lab = stereo_labels(ala);
  ASSERT_EQ(lab[ca], "S");
  invert_centre(ala, ca);
  lab = stereo_labels(ala);
  EXPECT_EQ(lab[ca], "R");
}

TEST(Edit, TacticityOfGrownPolystyrene) {
  EXPECT_EQ(tacticity(polystyrene(Tacticity::Isotactic)).label, "isotactic");
  EXPECT_EQ(tacticity(polystyrene(Tacticity::Syndiotactic)).label, "syndiotactic");
  System at = polystyrene(Tacticity::Atactic);
  const auto T = tacticity(at);
  ASSERT_EQ(T.chains.size(), 1u);
  EXPECT_GE(T.chains[0].centres.size(), 10u);
  EXPECT_EQ(T.m + T.r, int(T.chains[0].centres.size()) - 1);
  set_tacticity(at, true);
  EXPECT_EQ(tacticity(at).label, "isotactic");
  set_tacticity(at, false);
  EXPECT_EQ(tacticity(at).label, "syndiotactic");
}

TEST(Edit, Selections) {
  const System tol = molecule("Cc1ccccc1");
  const auto ring = select_smarts(tol, "c1ccccc1");
  EXPECT_EQ(std::count(ring.begin(), ring.end(), 1), 6);
  const auto c = select_element(tol, "C");
  EXPECT_EQ(std::count(c.begin(), c.end(), 1), 7);
  auto grown = select_grow(tol, ring, 1);
  EXPECT_EQ(std::count(grown.begin(), grown.end(), 1), 6 + 5 + 1);   // the ring, its 5 H and the methyl C
  const auto near = select_within(tol, ring, 1.2);
  EXPECT_EQ(std::count(near.begin(), near.end(), 1), 6 + 5);          // ring H within 1.2 Å, the methyl C (1.5 Å) not
  EXPECT_THROW(select_element(tol, "Xx"), EditError);
}
