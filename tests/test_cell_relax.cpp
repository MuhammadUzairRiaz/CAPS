#include <gtest/gtest.h>

#include <cmath>

#include "caps/cell_relax.hpp"
#include "caps/uff.hpp"

using namespace caps;

namespace {
constexpr double kEps = 0.2385, kSig = 3.405, kRc = 8.5;
ForceField argon_ff(const System& s) {
  ForceField ff = default_forcefield(s);
  ff.type_names = {"Ar"};
  ff.lj.assign(1, {kEps, kSig});
  for (size_t i = 0; i < s.atoms.size(); ++i) ff.type_index[i] = 0, ff.atom_type[i] = "Ar";
  ff.bonds.clear(), ff.angles.clear(), ff.dihedrals.clear(), ff.impropers.clear(), ff.inversions.clear(), ff.bonds_x.clear(), ff.angles_x.clear();
  ff.pairs14.clear();
  for (auto& e : ff.excluded) e.clear();
  ff.charge.assign(s.atoms.size(), 0.0);
  return ff;
}
EnergyOptions lj(bool tail = false) {
  EnergyOptions e;
  e.cutoff = kRc, e.coulomb = false, e.tail = tail;
  return e;
}
// conventional fcc argon in the cell of these parameters (four atoms at the face centres, fractional)
System argon(const std::array<double, 6>& p) {
  System s;
  s.cell = standard_cell(p);
  for (Vec3 f : {Vec3{0, 0, 0}, Vec3{0.5, 0.5, 0}, Vec3{0.5, 0, 0.5}, Vec3{0, 0.5, 0.5}}) {
    Atom a;
    a.element = 18;
    a.pos = s.cell.a * f[0] + s.cell.b * f[1] + s.cell.c * f[2];
    s.atoms.push_back(a);
  }
  return s;
}
double energy(const std::array<double, 6>& p, bool tail = false) {
  System s = argon(p);
  CellRelaxOptions o;
  o.max_steps = 0;
  return relax_cell(s, argon_ff, lj(tail), o).energy_before;
}
}  // namespace

// The cell gradient against central differences of the energy at fixed fractional coordinates, in a sheared,
// stretched cell: every ∂E/∂p, lengths per Å and angles per radian
TEST(CellRelax, GradientMatchesDifferences) {
  const std::array<double, 6> p{5.6, 5.3, 5.5, 88.0, 97.0, 92.0};
  for (bool tail : {false, true}) {   // with the long-range tail too (its energy depends on the volume)
    System s = argon(p);
    CellRelaxOptions o;
    o.max_steps = 0;
    const auto r = relax_cell(s, argon_ff, lj(tail), o);
    for (int k = 0; k < 6; ++k) {
      const double d = k < 3 ? 1e-4 : 1e-3;   // Å, degrees
      auto pp = p, pm = p;
      pp[size_t(k)] += d, pm[size_t(k)] -= d;
      const double fd = (energy(pp, tail) - energy(pm, tail)) / (2 * d) / (k < 3 ? 1.0 : M_PI / 180);
      EXPECT_NEAR(r.gradient[size_t(k)], fd, 1e-5 * std::max(1.0, std::fabs(fd))) << k << " tail " << tail;
    }
  }
}

// Fcc argon relaxed from a = 5.6 Å with every parameter free stays cubic at the minimum of E(a) found by a golden-section
// search over the cubic lattice alone, with no pressure left
TEST(CellRelax, CubicArgonReachesTheMinimum) {
  double lo = 5.0, hi = 5.8;
  const double gr = (std::sqrt(5.0) - 1) / 2;
  for (int it = 0; it < 80; ++it) {
    const double x1 = hi - gr * (hi - lo), x2 = lo + gr * (hi - lo);
    if (energy({x1, x1, x1, 90, 90, 90}) < energy({x2, x2, x2, 90, 90, 90})) hi = x2; else lo = x1;
  }
  const double amin = (lo + hi) / 2;
  System s = argon({5.6, 5.6, 5.6, 90, 90, 90});
  CellRelaxOptions o;
  o.ftol = 1e-5;
  const auto r = relax_cell(s, argon_ff, lj(), o);
  ASSERT_TRUE(r.converged);
  for (int k = 0; k < 3; ++k) EXPECT_NEAR(r.after[size_t(k)], amin, 1e-4);
  for (int k = 3; k < 6; ++k) EXPECT_NEAR(r.after[size_t(k)], 90, 1e-6);
  EXPECT_NEAR(r.pressure_after, 0, 1.0);   // atm
  EXPECT_LT(r.energy_after, r.energy_before);
}

// Held parameters stay exactly as they were (c and β of a monoclinic start), the free ones carry no stress
TEST(CellRelax, HeldParametersStay) {
  const std::array<double, 6> p{5.6, 5.3, 5.5, 90.0, 97.0, 90.0};
  System s = argon(p);
  CellRelaxOptions o;
  o.hold = {false, false, true, false, true, false};
  o.ftol = 1e-5;
  const auto r = relax_cell(s, argon_ff, lj(), o);
  ASSERT_TRUE(r.converged);
  EXPECT_NEAR(r.after[2], 5.5, 1e-12);
  EXPECT_NEAR(r.after[4], 97.0, 1e-9);
  for (int k : {0, 1, 3, 5}) EXPECT_LT(std::fabs(r.gradient[size_t(k)]), 1e-3) << k;
  EXPECT_GT(std::fabs(r.gradient[4]), 1e-2);   // β is held against a force
  // the cell keeps its orientation: a stays along x
  EXPECT_NEAR(s.cell.a[1], 0, 1e-12);
  EXPECT_NEAR(s.cell.a[2], 0, 1e-12);
}
