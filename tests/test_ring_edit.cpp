#include <gtest/gtest.h>

#include <cmath>
#include <map>

#include "caps/edit.hpp"
#include "caps/elements.hpp"
#include "caps/molecule.hpp"
#include "caps/ring_edit.hpp"

using namespace caps;

namespace {
System molecule(const std::string& smiles) {
  MolGraph g = parse_smiles(smiles);
  add_hydrogens(g);
  return molecule_system(g, embed(g));
}
std::string formula(const System& s) {   // Hill order
  std::map<std::string, int> n;
  for (const auto& a : s.atoms) ++n[element(a.element).symbol];
  std::string f;
  for (const char* first : {"C", "H"})
    if (n.count(first)) { f += first + (n[first] > 1 ? std::to_string(n[first]) : ""); n.erase(first); }
  for (const auto& [el, c] : n) f += el + (c > 1 ? std::to_string(c) : "");
  return f;
}
double dist(const System& s, uint32_t i, uint32_t j) { return norm(s.atoms[i].pos - s.atoms[j].pos); }
double angle(const System& s, uint32_t i, uint32_t j, uint32_t k) {
  const Vec3 a = s.atoms[i].pos - s.atoms[j].pos, b = s.atoms[k].pos - s.atoms[j].pos;
  return std::acos(dot(a, b) / (norm(a) * norm(b))) * 180 / M_PI;
}
double dihedral(const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& d) {
  const Vec3 b1 = b - a, b2 = c - b, b3 = d - c;
  const Vec3 n1 = cross(b1, b2), n2 = cross(b2, b3);
  return std::atan2(dot(cross(n1, n2), b2 * (1 / norm(b2))), dot(n1, n2)) * 180 / M_PI;
}
int bonds_of(const System& s, uint32_t a) {
  int n = 0;
  for (const auto& b : s.bonds) n += b.i == a || b.j == a;
  return n;
}
// an ideal methane: C at the origin, the hydrogens at 1.09 Å along the cube diagonals
System methane() {
  System s;
  Atom c;
  c.element = 6, c.mol = 1;
  s.atoms.push_back(c);
  const double h = 1.09 / std::sqrt(3.0);
  for (Vec3 d : {Vec3{1, 1, 1}, Vec3{1, -1, -1}, Vec3{-1, 1, -1}, Vec3{-1, -1, 1}}) {
    Atom x;
    x.element = 1, x.mol = 1, x.pos = d * h;
    s.atoms.push_back(x);
    s.bonds.push_back({0, uint32_t(s.atoms.size() - 1), 1});
  }
  return s;
}
}  // namespace

// Ethane + a three-ring of one O on its C–C: ethylene oxide, C2H4O; before any clean-up the O is the apex of the
// equilateral triangle on the C–C bond (the regular polygon with side |CC|)
TEST(RingEdit, EpoxideFusedOnEthane) {
  System s = molecule("CC");
  const double cc = dist(s, 0, 1);
  RingSpec sp;
  sp.size = 3, sp.members = "O", sp.clean = false;
  const auto r = fuse_ring(s, 0, 1, sp);
  EXPECT_EQ(formula(s), "C2H4O");
  ASSERT_EQ(r.ring.size(), 3u);
  EXPECT_EQ(r.removed, 2);
  const uint32_t o = r.ring[2];
  EXPECT_EQ(s.atoms[o].element, 8);
  EXPECT_NEAR(dist(s, 0, o), cc, 1e-12);
  EXPECT_NEAR(dist(s, 1, o), cc, 1e-12);
  EXPECT_EQ(bonds_of(s, o), 2);
  EXPECT_EQ(bonds_of(s, 0), 4);
  EXPECT_EQ(bonds_of(s, 1), 4);
  // cleaned with UFF the ring stays a three-ring
  System t = molecule("CC");
  sp.clean = true;
  const auto r2 = fuse_ring(t, 0, 1, sp);
  EXPECT_LT(dist(t, 0, r2.ring[2]), 1.6);
  EXPECT_LT(angle(t, 0, r2.ring[2], 1), 70);
}

// A six-ring of carbons on ethane's C–C is cyclohexane, C6H12: cleaned, it is a chair (ring torsions alternating in
// sign, each 40–70°)
TEST(RingEdit, CyclohexaneFusedIsAChair) {
  System s = molecule("CC");
  RingSpec sp;
  sp.size = 6;
  const auto r = fuse_ring(s, 0, 1, sp);
  EXPECT_EQ(formula(s), "C6H12");
  ASSERT_EQ(r.ring.size(), 6u);
  double sign = 0;
  for (int k = 0; k < 6; ++k) {
    const double t = dihedral(s.atoms[r.ring[k]].pos, s.atoms[r.ring[(k + 1) % 6]].pos, s.atoms[r.ring[(k + 2) % 6]].pos, s.atoms[r.ring[(k + 3) % 6]].pos);
    EXPECT_GT(std::fabs(t), 40) << k;
    EXPECT_LT(std::fabs(t), 70) << k;
    if (k > 0) EXPECT_LT(t * sign, 0) << "ring torsions alternate in a chair";
    sign = t;
  }
}

// An aromatic six-ring "C C C N" fused onto a C–C bond of benzene: quinoline, C9H7N, the 11 ring bonds of the two
// rings aromatic, the nitrogen without hydrogen; before clean-up a regular hexagon with benzene's side
TEST(RingEdit, PyridineFusedOnBenzeneIsQuinoline) {
  System s = molecule("c1ccccc1");
  const double cc = dist(s, 0, 1);
  RingSpec sp;
  sp.size = 6, sp.members = "C C C N", sp.aromatic = true, sp.clean = false;
  const auto r = fuse_ring(s, 0, 1, sp);
  EXPECT_EQ(formula(s), "C9H7N");
  ASSERT_EQ(r.ring.size(), 6u);
  for (int k = 0; k < 6; ++k) EXPECT_NEAR(dist(s, r.ring[k], r.ring[(k + 1) % 6]), cc, 1e-9);
  for (int k = 0; k < 6; ++k) EXPECT_NEAR(angle(s, r.ring[k], r.ring[(k + 1) % 6], r.ring[(k + 2) % 6]), 120, 1e-7);
  int aromatic = 0;
  for (const auto& b : s.bonds) aromatic += b.order == 4;
  EXPECT_EQ(aromatic, 11);
  EXPECT_EQ(bonds_of(s, r.ring[5]), 2);
  EXPECT_THROW(fuse_ring(s, 0, 2, sp), EditError);   // not bonded
  sp.size = 5;
  EXPECT_THROW(fuse_ring(s, 0, 1, sp), EditError);   // aromatic rings of six only
}

// Spiro rings on a cyclohexane carbon: four carbons give spiro[3.5]nonane, C9H16, the new ring a square with side
// 2 r_cov(C) through the shared atom; a γ-lactam (N, C=O, C, C) gives 2-azaspiro[4.5]decan-3-one, C9H15NO, with one H
// on N and none on the carbonyl carbon
TEST(RingEdit, SpiroRings) {
  System s = molecule("C1CCCCC1");
  RingSpec sp;
  sp.size = 4, sp.clean = false;
  const auto r = spiro_ring(s, 0, sp);
  EXPECT_EQ(formula(s), "C9H16");
  ASSERT_EQ(r.ring.size(), 4u);
  const double L = 2 * element(6).covalent;
  for (int k = 0; k < 4; ++k) EXPECT_NEAR(dist(s, r.ring[k], r.ring[(k + 1) % 4]), L, 1e-9);
  for (int k = 0; k < 4; ++k) EXPECT_NEAR(angle(s, r.ring[k], r.ring[(k + 1) % 4], r.ring[(k + 2) % 4]), 90, 1e-7);
  EXPECT_EQ(bonds_of(s, r.ring[0]), 4);

  System t = molecule("C1CCCCC1");
  sp.size = 5, sp.members = "N C=O C C", sp.clean = true;
  const auto q = spiro_ring(t, 0, sp);
  EXPECT_EQ(formula(t), "C9H15NO");
  ASSERT_EQ(q.ring.size(), 5u);
  int nh = 0, ch = 0;
  for (const auto& b : t.bonds)
    for (auto [x, y] : {std::pair{b.i, b.j}, std::pair{b.j, b.i}}) {
      if (x == q.ring[1] && t.atoms[y].element == 1) ++nh;
      if (x == q.ring[2] && t.atoms[y].element == 1) ++ch;
    }
  EXPECT_EQ(nh, 1);
  EXPECT_EQ(ch, 0);
  EXPECT_THROW(spiro_ring(t, q.ring[2], sp), EditError);   // the carbonyl carbon has no room
}

// One more sulfur in dimethyl disulfide: dimethyl trisulfide, C2H6S3. Before clean-up both S–S bonds are 2 r_cov(S),
// the S–S–S angle tetrahedral, and each side kept its shape
TEST(RingEdit, InsertSulfurIntoDisulfide) {
  System s = molecule("CSSC");
  const System before = s;
  const auto r = insert_atom(s, 1, 2, 16);
  EXPECT_EQ(formula(s), "C2H6S3");
  ASSERT_EQ(r.ring.size(), 3u);
  const uint32_t x = r.ring[1];
  const double ss = 2 * element(16).covalent;
  EXPECT_NEAR(dist(s, 1, x), ss, 1e-12);
  EXPECT_NEAR(dist(s, x, 2), ss, 1e-12);
  EXPECT_NEAR(angle(s, 1, x, 2), std::acos(-1.0 / 3.0) * 180 / M_PI, 1e-9);
  EXPECT_NEAR(dist(s, 2, 3), dist(before, 2, 3), 1e-12);
  EXPECT_NEAR(dist(s, 0, 1), dist(before, 0, 1), 1e-12);
  // in a ring the shared atoms stay: cyclohexane + C is cycloheptane
  System c = molecule("C1CCCCC1");
  const Vec3 p0 = c.atoms[0].pos, p1 = c.atoms[1].pos;
  const auto q = insert_atom(c, 0, 1, 6);
  EXPECT_EQ(formula(c), "C7H14");
  EXPECT_EQ(c.atoms[0].pos, p0);
  EXPECT_EQ(c.atoms[1].pos, p1);
  EXPECT_NEAR(dist(c, 0, q.ring[1]), 2 * element(6).covalent, 1e-12);
}

// Two ideal methanes joined by a hydrogen each: ethane, C2H6, the C–C bond 2 r_cov(C), the piece rigid, turned to the
// staggered form (every H–C–C–H dihedral 60° or 180°); a torsion asked is set exactly
TEST(RingEdit, JoinTwoMethanesIsStaggeredEthane) {
  const System b = methane();
  EXPECT_EQ(default_leaving_atom(b), 1u);
  System a = methane();
  const auto r = join_pieces(a, b, 1, default_leaving_atom(b));
  EXPECT_EQ(formula(a), "C2H6");
  ASSERT_EQ(r.ring.size(), 2u);
  const uint32_t c1 = r.ring[0], c2 = r.ring[1];
  EXPECT_NEAR(dist(a, c1, c2), 2 * element(6).covalent, 1e-12);
  std::vector<uint32_t> h1, h2;
  for (const auto& bd : a.bonds)
    for (auto [x, y] : {std::pair{bd.i, bd.j}, std::pair{bd.j, bd.i}})
      if (a.atoms[y].element == 1) (x == c1 ? h1 : h2).push_back(y);
  ASSERT_EQ(h1.size(), 3u);
  ASSERT_EQ(h2.size(), 3u);
  for (uint32_t x : h2) EXPECT_NEAR(dist(a, c2, x), 1.09, 1e-12);
  for (uint32_t x : h1)
    for (uint32_t y : h2) {
      const double t = std::fabs(dihedral(a.atoms[x].pos, a.atoms[c1].pos, a.atoms[c2].pos, a.atoms[y].pos));
      EXPECT_TRUE(std::fabs(t - 60) < 1e-4 || std::fabs(t - 180) < 1e-4) << t;
    }
  System e = methane();
  JoinOptions o;
  o.torsion = 0;
  const auto q = join_pieces(e, b, 1, 1, o);
  // the reference atoms: the first neighbour of each carbon other than the leaving atom (methane's atom 3 → index 1
  // after atom 2 left; the piece's atom 3, its first hydrogen kept, is the second atom added)
  EXPECT_NEAR(dihedral(e.atoms[1].pos, e.atoms[q.ring[0]].pos, e.atoms[q.ring[1]].pos, e.atoms[q.added[1]].pos), 0, 1e-9);
  EXPECT_THROW(join_pieces(e, b, 0, 1), EditError);   // a carbon is no leaving atom
}
