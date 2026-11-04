#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>

#include "caps/crystal.hpp"
#include "caps/molecule.hpp"
#include "caps/powder_index.hpp"
#include "caps/spacegroup.hpp"

using namespace caps;

namespace {
MeasuredPattern synthetic(const System& s, double from, double to, double W) {
  RefineOptions o;
  o.wavelength = 1.540593, o.W = W, o.Y = 0.0;
  MeasuredPattern d;
  for (double t = from; t <= to + 1e-9; t += 0.01) d.two_theta.push_back(t);
  RefinePhase p;
  p.cell = s;
  d.y = calculated_pattern(d.two_theta, {p}, {1e-3}, o);
  for (auto& v : d.y) v += 10;
  return d;
}
}  // namespace

// Copper's lines index on a primitive cubic cell of its lattice constant, ranked first by M20
TEST(PowderIndex, CopperIsCubic) {
  const System cu = read_cif(std::string(CAPS_SOURCE_DIR) + "/data/crystals/copper.cif");
  const auto pk = find_peaks(synthetic(cu, 30, 140, 0.002));
  std::vector<double> tt;
  for (const auto& p : pk) tt.push_back(p.two_theta);
  ASSERT_GE(tt.size(), 6u);
  IndexOptions o;
  o.max_length = 12;
  o.systems = {"cubic", "tetragonal", "orthorhombic"};
  const IndexResult r = index_pattern(tt, o);
  ASSERT_FALSE(r.solutions.empty());
  const auto& s = r.solutions[0];
  std::printf("%zu lines; best %s a %.5f (M20 %.1f, F_N %.1f, %d indexed), %.2f s\n", tt.size(), s.system.c_str(), s.cell[0], s.M20, s.FN, s.indexed, r.seconds);
  EXPECT_EQ(s.system, "cubic");
  EXPECT_NEAR(s.cell[0], 3.6149, 1e-3);
}

// Orthorhombic polyethylene's lines give its cell (2.534, 4.93, 7.40 Å) back
TEST(PowderIndex, PolyethyleneIsOrthorhombic) {
  const System pe = read_cif(std::string(CAPS_SOURCE_DIR) + "/data/crystals/polymers/pe-orthorhombic.cif");
  const auto pk = find_peaks(synthetic(pe, 10, 75, 0.002));
  std::vector<double> tt;
  for (const auto& p : pk) tt.push_back(p.two_theta);
  IndexOptions o;
  o.max_length = 9;
  o.systems = {"orthorhombic"};
  o.lines = 14;
  const IndexResult r = index_pattern(tt, o);
  ASSERT_FALSE(r.solutions.empty());
  auto c = r.solutions[0].cell;
  std::array<double, 3> l{c[0], c[1], c[2]};
  std::sort(l.begin(), l.end());
  std::printf("%zu lines; best a b c %.4f %.4f %.4f (M20 %.1f), %.2f s\n", tt.size(), l[0], l[1], l[2], r.solutions[0].M20, r.seconds);
  EXPECT_NEAR(l[0], 2.534, 3e-3);
  EXPECT_NEAR(l[1], 4.93, 3e-3);
  EXPECT_NEAR(l[2], 7.40, 3e-3);
}

// A dichloroethane crystal in P212121 made from a known placement: annealing finds a structure whose pattern fits it
TEST(PowderIndex, SolvesASyntheticMolecularCrystal) {
  System mol = build_molecule("ClCCCl").system;
  const std::array<double, 6> cp{9.0, 10.0, 11.0, 90, 90, 90};
  const Cell cell = cell_parameters(cp[0], cp[1], cp[2], cp[3], cp[4], cp[5]);
  // the truth: the molecule's centroid at (0.2, 0.3, 0.1), as built, and its symmetry copies
  Vec3 c0{0, 0, 0};
  for (const auto& a : mol.atoms) c0 = c0 + a.pos;
  c0 = c0 * (1.0 / double(mol.atoms.size()));
  const auto ops = hall_operations(find_space_group("P212121")->hall);
  System truth;
  truth.cell = cell;
  for (const auto& op : ops)
    for (const auto& a : mol.atoms) {
      const Vec3 f = cell.to_fractional(a.pos - c0 + cell.to_cartesian({0.2, 0.3, 0.1}));
      Vec3 g;
      for (int r = 0; r < 3; ++r) g[size_t(r)] = op.R[r][0] * f[0] + op.R[r][1] * f[1] + op.R[r][2] * f[2] + op.t[r];
      Atom b = a;
      b.pos = cell.to_cartesian(g);
      truth.atoms.push_back(b);
    }
  const MeasuredPattern d = synthetic(truth, 5, 50, 0.01);
  SolveOptions o;
  o.space_group = "P212121";
  o.refine.wavelength = 1.540593, o.refine.W = 0.01, o.refine.Y = 0, o.refine.background = 2;
  o.steps = 6000, o.restarts = 2, o.contact = 2.5;
  {
    RefinePhase tp;
    tp.cell = truth, tp.lattice = false;
    RefineOptions ro = o.refine;
    ro.refine_U = ro.refine_V = ro.refine_W = ro.refine_X = ro.refine_Y = ro.refine_zero = false;
    std::printf("the truth's Rwp %.2e\n", refine_pattern(d, {tp}, ro).rwp);
  }
  const SolveResult r = solve_structure(d, mol, cp, o);
  std::printf("Rwp %.4f after %d moves (%d accepted, %d torsions)\n", r.rwp, r.steps, r.accepted, r.torsions);
  EXPECT_LT(r.rwp, 0.05);
  EXPECT_EQ(r.crystal.atoms.size(), truth.atoms.size());
}
