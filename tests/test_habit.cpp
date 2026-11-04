#include <gtest/gtest.h>

#include <cmath>

#include "caps/habit.hpp"
#include "caps/uff.hpp"

using namespace caps;

namespace {
ForceField lj_ff(const System& s, double eps, double sig) {
  ForceField ff = default_forcefield(s);
  ff.type_names = {"Ar"};
  ff.lj.assign(1, {eps, sig});
  for (size_t i = 0; i < s.atoms.size(); ++i) ff.type_index[i] = 0, ff.atom_type[i] = "Ar";
  ff.charge.assign(s.atoms.size(), 0.0);
  ff.bonds.clear(), ff.angles.clear(), ff.dihedrals.clear(), ff.impropers.clear(), ff.inversions.clear();
  ff.pairs14.clear();
  for (auto& e : ff.excluded) e.clear();
  return ff;
}
System cubic(double a, const std::vector<Vec3>& frac) {
  System s;
  s.cell.a = {a, 0, 0}, s.cell.b = {0, a, 0}, s.cell.c = {0, 0, a};
  for (const auto& f : frac) {
    Atom at;
    at.element = 18;
    at.pos = s.cell.to_cartesian(f);
    s.atoms.push_back(at);
  }
  return s;
}
const HabitFace* form(const HabitResult& r, std::array<int, 3> h) {
  for (const auto& f : r.faces)
    for (const auto& m : f.members)
      if (m == h) return &f;
  return nullptr;
}
}  // namespace

// BFDH of fcc (Fm-3m, conventional cell): {111} at a/√3 first, then {100} and {110} at half their spacings (the F centring
// extinguishes (100) and (110)), with multiplicities 8, 6 and 12
TEST(Habit, BfdhFccWithExtinctions) {
  const double a = 5.26;
  const System s = cubic(a, {{0, 0, 0}, {0, 0.5, 0.5}, {0.5, 0, 0.5}, {0.5, 0.5, 0}});
  HabitOptions o;
  o.max_index = 3;
  const HabitResult r = crystal_habit(s, o);
  std::printf("%s, %d operations; first forms:", r.space_group.c_str(), r.operations);
  for (size_t i = 0; i < 4 && i < r.faces.size(); ++i) std::printf(" {%d%d%d} %.4f", r.faces[i].hkl[0], r.faces[i].hkl[1], r.faces[i].hkl[2], r.faces[i].d_eff);
  std::printf("\n");
  EXPECT_EQ(r.operations, 192);
  ASSERT_GE(r.faces.size(), 3u);
  EXPECT_EQ(r.faces[0].members.size(), 8u);
  EXPECT_NEAR(r.faces[0].d_eff, a / std::sqrt(3.0), 1e-9);
  EXPECT_EQ(r.faces[1].members.size(), 6u);
  EXPECT_EQ(r.faces[1].order, 2);
  EXPECT_NEAR(r.faces[1].d_eff, a / 2, 1e-9);
  EXPECT_EQ(r.faces[2].members.size(), 12u);
  EXPECT_NEAR(r.faces[2].d_eff, a / (2 * std::sqrt(2.0)), 1e-9);
  // the BFDH shape: {111} and {100} only, the {111} the larger share
  EXPECT_GT(r.faces[0].area_bfdh, r.faces[1].area_bfdh);
  EXPECT_NEAR(r.faces[0].area_bfdh + r.faces[1].area_bfdh, 1.0, 1e-9);
}

// The Wulff construction of a cube truncated by {111}: areas and volume as geometry gives them
TEST(Habit, WulffTruncatedCube) {
  std::vector<Vec3> n;
  std::vector<double> d;
  for (int ax = 0; ax < 3; ++ax)
    for (int sg : {-1, 1}) {
      Vec3 v{0, 0, 0};
      v[size_t(ax)] = sg;
      n.push_back(v), d.push_back(1.0);
    }
  const WulffShape cube = wulff_shape(n, d);
  EXPECT_NEAR(cube.area, 24, 1e-9);
  EXPECT_NEAR(cube.volume, 8, 1e-9);
  EXPECT_NEAR(cube.aspect, 1, 1e-9);
  const double t = 1.5, c = 3 - std::sqrt(3.0) * t;   // corners cut c along each edge
  for (int x : {-1, 1})
    for (int y : {-1, 1})
      for (int z : {-1, 1}) n.push_back(Vec3{double(x), double(y), double(z)} * (1 / std::sqrt(3.0))), d.push_back(t);
  const WulffShape w = wulff_shape(n, d);
  double sq = 0, tri = 0;
  for (size_t k = 0; k < w.facets.size(); ++k) (w.facet_plane[k] < 6 ? sq : tri) += w.facet_area[k];
  std::printf("squares %.9f (exact %.9f), triangles %.9f (%.9f), volume %.9f (%.9f), %zu vertices\n", sq, 6 * (4 - 2 * c * c), tri,
              8 * std::sqrt(3.0) / 2 * c * c, w.volume, 8 - 8 * c * c * c / 6, w.vertices.size());
  EXPECT_NEAR(sq, 6 * (4 - 2 * c * c), 1e-9);
  EXPECT_NEAR(tri, 8 * std::sqrt(3.0) / 2 * c * c, 1e-9);
  EXPECT_NEAR(w.volume, 8 - 8 * c * c * c / 6, 1e-9);
  EXPECT_EQ(w.vertices.size(), 24u);
  EXPECT_EQ(w.facets.size(), 14u);
  // two parallel planes bound nothing
  EXPECT_FALSE(wulff_shape({{1, 0, 0}, {-1, 0, 0}}, {1, 1}).ok());
}

// Simple cubic with nearest-neighbour bonds u only (the cut-off between a and a√2): E_latt = 3u per atom; the attachment
// energies of {100}, {110}, {111} are u, 2u, 3u and the surface energies −u/2a², −u/√2a², −√3u/2a² (bonds broken per area)
TEST(Habit, SimpleCubicBrokenBonds) {
  const double a = 3.9, eps = 0.24, sig = 3.4;
  const System s = cubic(a, {{0, 0, 0}});
  HabitOptions o;
  o.max_index = 1;
  o.attachment = o.surface = true;
  o.energy.cutoff = 4.7, o.energy.coulomb = false, o.energy.tail = false;
  o.field_for = [&](const System& x) { return lj_ff(x, eps, sig); };
  const HabitResult r = crystal_habit(s, o);
  // CAPS's Lennard-Jones without tail corrections is shifted to zero at the cut-off: the bond is φ(a) − φ(rc)
  auto phi = [&](double r) { const double x = std::pow(sig / r, 6); return 4 * eps * (x * x - x); };
  const double u = phi(a) - phi(o.energy.cutoff);
  const double J = 4184.0 / 6.02214076e23 / 1e-20;
  const HabitFace *f100 = form(r, {1, 0, 0}), *f110 = form(r, {1, 1, 0}), *f111 = form(r, {1, 1, 1});
  ASSERT_TRUE(f100 && f110 && f111);
  std::printf("u %.6f · E_latt %.6f · E_att %.6f %.6f %.6f · γ %.6f %.6f %.6f J/m² (exact %.6f %.6f %.6f)\n", u, r.e_latt, f100->e_att, f110->e_att, f111->e_att,
              f100->gamma, f110->gamma, f111->gamma, -u / (2 * a * a) * J, -u / (std::sqrt(2.0) * a * a) * J, -std::sqrt(3.0) * u / (2 * a * a) * J);
  for (const auto& n : r.notes) std::printf("note: %s\n", n.c_str());
  EXPECT_NEAR(r.e_latt, 3 * u, 1e-9);
  EXPECT_NEAR(f100->e_att, u, 1e-9);
  EXPECT_NEAR(f110->e_att, 2 * u, 1e-9);
  EXPECT_NEAR(f111->e_att, 3 * u, 1e-9);
  EXPECT_NEAR(f100->gamma, -u / (2 * a * a) * J, 1e-9);
  EXPECT_NEAR(f110->gamma, -u / (std::sqrt(2.0) * a * a) * J, 1e-9);
  EXPECT_NEAR(f111->gamma, -std::sqrt(3.0) * u / (2 * a * a) * J, 1e-9);
  // growth ∝ |E_att|: the {100} cube alone (the {110}, {111} planes lie outside it)
  EXPECT_NEAR(f100->area_growth, 1.0, 1e-9);
  EXPECT_NEAR(r.growth.volume, 8 * std::pow(std::fabs(u), 3), 1e-6 * std::pow(std::fabs(u), 3));
}

// A particle carved in the cube shape: the lattice points within the cube, (2·2 + 1)³ = 125
TEST(Habit, ParticleCube) {
  const double a = 3.9;
  const System s = cubic(a, {{0, 0, 0}});
  std::vector<Vec3> n;
  std::vector<double> d;
  for (int ax = 0; ax < 3; ++ax)
    for (int sg : {-1, 1}) {
      Vec3 v{0, 0, 0};
      v[size_t(ax)] = sg;
      n.push_back(v), d.push_back(1.0);
    }
  const System p = carve_particle(s, n, d, 2.5 * a * std::sqrt(3.0));   // half side 2.5 a
  EXPECT_EQ(p.atoms.size(), 125u);
  EXPECT_FALSE(p.cell.valid());
}
