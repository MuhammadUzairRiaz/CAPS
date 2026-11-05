#include <gtest/gtest.h>

#include <cmath>

#include "caps/field.hpp"
#include "caps/uff.hpp"

using namespace caps;

namespace {
// a dipolar slab: +1 charges on a 4 × 4 square net (3 Å) at z = 10, −1 charges on the net shifted by half a cell at
// z = 12.5 (extra: one more +1 at z = 11 for a net charge); a 12 × 12 Å cell, Lz long
System slab(double Lz, double extra_z = -1, double shift = 0) {
  System s;
  s.cell.a = {12, 0, 0}, s.cell.b = {0, 12, 0}, s.cell.c = {0, 0, Lz};
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j) {
      Atom p, m;
      p.element = 11, p.charge = 1, p.pos = {3.0 * i, 3.0 * j, 10 + shift};
      m.element = 17, m.charge = -1, m.pos = {3.0 * i + 1.5, 3.0 * j + 1.5, 12.5 + shift};
      s.atoms.push_back(p), s.atoms.push_back(m);
    }
  if (extra_z > 0) {
    Atom q;
    q.element = 11, q.charge = 1, q.pos = {1.5, 0.0, extra_z + shift};
    s.atoms.push_back(q);
  }
  s.has_charges = true;
  return s;
}
ForceField ions(const System& s) {
  ForceField ff = default_forcefield(s);
  for (auto& lj : ff.lj) lj = {0, 1};
  ff.bonds.clear(), ff.angles.clear(), ff.dihedrals.clear(), ff.impropers.clear(), ff.inversions.clear();
  ff.pairs14.clear();
  for (auto& e : ff.excluded) e.clear();
  ff.charge.clear();
  for (const auto& a : s.atoms) ff.charge.push_back(a.charge);
  return ff;
}
double energy(const System& s, double slab_factor, std::vector<double>* f = nullptr) {
  const ForceField ff = ions(s);
  EnergyOptions o;
  o.cutoff = 6, o.tail = false, o.electrostatics = EnergyOptions::Electrostatics::PME, o.ewald_rtol = 1e-8, o.pme_spacing = 0.3, o.pme_order = 6;
  o.slab = slab_factor;
  Evaluator ev(ff, o);
  std::vector<double> x, ff_;
  for (const auto& a : s.atoms) x.insert(x.end(), a.pos.begin(), a.pos.end());
  const double e = ev.compute(x, s.cell, ff_).coulomb;
  if (f) *f = ff_;
  return e;
}
}  // namespace

// The dipolar slab: with the correction the energy hardly depends on the vacuum (Lz 40 or 80 Å, factor 3), without it
// the images' dipole–dipole energy (∝ 1/V) moves it by far more
TEST(SlabEwald, VacuumIndependence) {
  const double c40 = energy(slab(40), 3), c80 = energy(slab(80), 3);
  const double u40 = energy(slab(40), 0), u80 = energy(slab(80), 0);
  EXPECT_LT(std::fabs(c40 - c80), 0.02 * std::fabs(u40 - u80)) << c40 << " " << c80 << " " << u40 << " " << u80;
  EXPECT_LT(std::fabs(c40 - c80) / std::fabs(c80), 1e-3);
  // the lengthening factor matters little once the vacuum is wide (3 against 5)
  EXPECT_LT(std::fabs(energy(slab(40), 5) - c40) / std::fabs(c40), 1e-3);
}

// A slab with a net charge moved 3 Å along z inside its vacuum: the corrected energy (with the net-charge terms) is
// unchanged; the forces are −∂E/∂z of an atom by central differences
TEST(SlabEwald, ChargedSlabTranslationAndForces) {
  const double e0 = energy(slab(50, 11), 3), e1 = energy(slab(50, 11, 3.0), 3);
  EXPECT_NEAR(e0, e1, 1e-6 * std::fabs(e0));
  System s = slab(50, 11);
  std::vector<double> f;
  energy(s, 3, &f);
  for (size_t i : {0u, 5u, 32u}) {
    const double h = 1e-4;
    System p = s, m = s;
    p.atoms[i].pos[2] += h, m.atoms[i].pos[2] -= h;
    const double fd = -(energy(p, 3) - energy(m, 3)) / (2 * h);
    EXPECT_NEAR(f[3 * i + 2], fd, 1e-4 * std::max(1.0, std::fabs(fd))) << i;
  }
}

// Against LAMMPS (11 Feb 2026 development build, kspace_style ewald 1e-10, pair_style coul/long 6.0, units real) on the charged slab in a
// 12 × 12 × 40 Å cell: boundary p p f with kspace_modify slab 3.0 gives 2452.6201 kcal/mol, the plain periodic sum
// 1914.1118; CAPS's PME (rtol 1e-8, grid 0.3 Å, order 6) within 2e-3 of both
TEST(SlabEwald, AgreesWithLammps) {
  const System s = slab(40, 11);
  EXPECT_NEAR(energy(s, 3), 2452.6201, 2e-3);
  EXPECT_NEAR(energy(s, 0), 1914.1118, 2e-3);
}
