#include <gtest/gtest.h>

#include <cmath>
#include <memory>

#include "caps/crystal.hpp"
#include "caps/molecule.hpp"
#include "caps/reaction_path.hpp"
#include "caps/relax.hpp"
#include "caps/uff.hpp"

using namespace caps;

namespace {
System molecule(const std::string& smiles) {
  MolGraph g = parse_smiles(smiles);
  add_hydrogens(g);
  return molecule_system(g, embed(g));
}
// a fixed shuffle of n indices: out[k] is the old index placed at k
std::vector<size_t> shuffle(size_t n) {
  std::vector<size_t> p(n);
  for (size_t k = 0; k < n; ++k) p[k] = (k * 7 + 3) % n;
  return p;
}
Vec3 rot(const Vec3& x, double t) { return {std::cos(t) * x[0] - std::sin(t) * x[1], std::sin(t) * x[0] + std::cos(t) * x[1], x[2]}; }
double dihedral(const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& d) {
  const Vec3 b1 = b - a, b2 = c - b, b3 = d - c;
  const Vec3 n1 = cross(b1, b2), n2 = cross(b2, b3);
  return std::atan2(dot(cross(n1, n2), b2 * (1 / norm(b2))), dot(n1, n2)) * 180 / M_PI;
}
}  // namespace

// A molecule listed in another order, turned and moved: every atom found again (the methyl hydrogens told apart by
// position), RMSD zero
TEST(ReactionPath, MatchShuffledRotatedMolecule) {
  const System a = molecule("CCO");
  const size_t n = a.atoms.size();
  const auto p = shuffle(n);
  System b = a;
  b.bonds.clear();
  for (size_t k = 0; k < n; ++k) {
    b.atoms[k] = a.atoms[p[k]];
    b.atoms[k].pos = rot(a.atoms[p[k]].pos, 1.1) + Vec3{3, -2, 5};
  }
  const AtomMatch m = match_atoms(a, b);
  for (size_t i = 0; i < n; ++i) EXPECT_EQ(p[size_t(m.map[i])], i);
  EXPECT_LT(m.rmsd, 1e-6);
  EXPECT_TRUE(m.by_graph);
}

// A crystal's atoms in another order, a few wrapped by a lattice vector and all moved slightly: found by minimum image
TEST(ReactionPath, MatchPeriodic) {
  const System a = read_cif(std::string(CAPS_SOURCE_DIR) + "/data/crystals/alpha-quartz.cif");
  const size_t n = a.atoms.size();
  const auto p = shuffle(n);
  System b = a;
  for (size_t k = 0; k < n; ++k) {
    b.atoms[k] = a.atoms[p[k]];
    b.atoms[k].pos = b.atoms[k].pos + Vec3{0.03 * std::sin(double(k)), 0.02, -0.04 * std::cos(double(k))} + (k % 3 == 0 ? a.cell.a : Vec3{0, 0, 0});
  }
  const AtomMatch m = match_atoms(a, b);
  for (size_t i = 0; i < n; ++i) EXPECT_EQ(p[size_t(m.map[i])], i);
  const System r = reorder_atoms(b, m.map, a);
  for (size_t i = 0; i < n; ++i) EXPECT_LT(norm(r.atoms[i].pos - a.atoms[i].pos), 0.06);
}

// A dimer turned by 90°: the straight line shortens the bond to 1.5 cos 45° at the middle; LST and IDPP keep it at 1.5 Å
// on every image (the pair distance is the same at both ends); the ends stay where they are
TEST(ReactionPath, InterpolationKeepsDistances) {
  System a;
  Atom c1, c2;
  c1.element = c2.element = 6;
  c1.pos = {-0.75, 0, 0}, c2.pos = {0.75, 0, 0};
  a.atoms = {c1, c2};
  System b = a;
  b.atoms[0].pos = {0, -0.75, 0}, b.atoms[1].pos = {0, 0.75, 0};
  PathOptions o;
  o.images = 5;
  o.method = "linear";
  const auto li = interpolate_path(a, b, o);
  EXPECT_NEAR(norm(li[3][1] - li[3][0]), 1.5 * std::cos(M_PI / 4), 1e-12);
  for (const char* method : {"lst", "idpp"}) {
    o.method = method, o.tolerance = 1e-6;
    PathReport rep;
    const auto path = interpolate_path(a, b, o, &rep);
    ASSERT_EQ(path.size(), 7u);
    for (const auto& img : path) EXPECT_NEAR(norm(img[1] - img[0]), 1.5, 1e-3) << method;
    EXPECT_EQ(path.front()[0], a.atoms[0].pos);
    EXPECT_EQ(path.back()[1], b.atoms[1].pos);
  }
}

// The climbing image finds the saddle of E = (x² − 1)² + 5(y − c(1 − x²))² + z² at (0, c, 0), E = 1, from a straight
// band along y = 0 that passes 5c² higher; both barriers are 1
TEST(ReactionPath, ClimbingImageFindsTheSaddle) {
  const double c = 0.5;
  const ImageEnergy E = [&](int, const std::vector<Vec3>& p, std::vector<Vec3>& f) {
    const double x = p[0][0], y = p[0][1], z = p[0][2], u = y - c * (1 - x * x);
    f.assign(1, Vec3{-(4 * x * (x * x - 1) + 10 * u * 2 * c * x), -(10 * u), -2 * z});
    return (x * x - 1) * (x * x - 1) + 5 * u * u + z * z;
  };
  std::vector<std::vector<Vec3>> path;
  for (int k = 0; k <= 8; ++k) path.push_back({Vec3{-1 + 0.25 * k, 0, 0}});
  NebOptions o;
  o.spring = 1, o.fmax = 1e-6, o.max_steps = 20000;
  const NebResult r = neb(path, E, o);
  ASSERT_TRUE(r.converged);
  ASSERT_GE(r.climbing, 0);
  const Vec3 s = r.images[size_t(r.climbing)][0];
  EXPECT_NEAR(s[0], 0, 1e-5);
  EXPECT_NEAR(s[1], c, 1e-5);
  EXPECT_NEAR(r.energy[size_t(r.climbing)], 1.0, 1e-9);
  EXPECT_NEAR(r.forward, 1.0, 1e-9);
  EXPECT_NEAR(r.reverse, 1.0, 1e-9);
}

// Ethane's methyl turned by 120° with UFF: the climbing image is the eclipsed form (H–C–C–H 0° by symmetry), above both
// ends, which have the same energy
TEST(ReactionPath, EthaneRotationWithTheForceField) {
  System a = molecule("CC");
  auto ff = std::make_shared<ForceField>(assign_uff(a));
  RelaxOptions ro;
  ro.field = ff, ro.ftol = 1e-4, ro.pushoff = false, ro.max_iterations = 20000;
  relax(a, ro);
  std::vector<size_t> hc1, hc2;   // the hydrogens on each carbon
  for (const auto& bd : a.bonds)
    for (int s = 0; s < 2; ++s) {
      const uint32_t x = s ? bd.j : bd.i, y = s ? bd.i : bd.j;
      if (a.atoms[x].element == 1) (y == 0 ? hc1 : hc2).push_back(x);
    }
  ASSERT_EQ(hc1.size(), 3u);
  ASSERT_EQ(hc2.size(), 3u);
  // B: the second methyl turned by 120° about the C–C axis
  System b = a;
  const Vec3 c0 = a.atoms[0].pos, axis = (a.atoms[1].pos - c0) * (1 / norm(a.atoms[1].pos - c0));
  const double t = 2 * M_PI / 3;
  for (size_t h : hc2) {
    const Vec3 v = a.atoms[h].pos - c0;
    b.atoms[h].pos = v * std::cos(t) + cross(axis, v) * std::sin(t) + axis * (dot(axis, v) * (1 - std::cos(t))) + c0;
  }
  PathOptions po;
  po.images = 7;
  const auto path = interpolate_path(a, b, po);
  EnergyOptions eo;
  eo.coulomb = false;
  NebOptions no;
  no.fmax = 0.01, no.max_steps = 4000;
  const NebResult r = neb_forcefield(a, path, *ff, eo, no);
  ASSERT_TRUE(r.converged);
  ASSERT_GE(r.climbing, 0);
  const auto& ts = r.images[size_t(r.climbing)];
  double best = 180;
  for (size_t h1 : hc1)
    for (size_t h2 : hc2) best = std::min(best, std::fabs(dihedral(ts[h1], ts[0], ts[1], ts[h2])));
  EXPECT_LT(best, 2.0);
  EXPECT_GT(r.forward, 0.5);
  EXPECT_NEAR(r.forward, r.reverse, 1e-3);
}
