#include <gtest/gtest.h>

#include <cmath>

#include "caps/dpd.hpp"

using namespace caps;

namespace {
// the start alone: no steps, positions in r_c
DpdReport start_of(DpdOptions o) {
  o.steps = 0, o.equilibration = 0, o.frame_every = 0, o.rc_angstrom = 1.0;
  return run_dpd(o);
}
double rmin(const Vec3& x, const Vec3& c, double L) {
  Vec3 d = x - c;
  for (auto& q : d) q -= L * std::round(q / L);
  return norm(d);
}
}  // namespace

// A core–shell droplet in a solvent: core beads (A, monomers) within 3 r_c of the centre, shell dimers (BB) between 3
// and 5 r_c, the solvent (C) outside 5 r_c — every bead in its own region
TEST(DpdRegions, CoreShellDroplet) {
  DpdOptions o;
  o.box = 16;
  o.regions = {{"core", "droplet", {0.5, 0.5, 0.5}, 3.0, 0, 3, 2, 2}, {"shell", "shell", {0.5, 0.5, 0.5}, 5.0, 3.0, 3, 2, 1}};
  o.species = {{"core", "A", 300, "core", ""}, {"shell", "BB", 600, "shell", ""}, {"solvent", "C", 9000, "", ""}};
  const auto r = start_of(o);
  const auto& t = r.frames;
  const Vec3 c{8, 8, 8};
  for (size_t i = 0; i < t.topology.atoms.size(); ++i) {
    const double d = rmin(t.positions[0][i], c, 16);
    const int z = t.topology.atoms[i].element;   // A → C, B → O, C → N
    if (z == 6) EXPECT_LE(d, 3.0) << i;
    else if (z == 8) { EXPECT_GE(d, 3.0) << i; EXPECT_LE(d, 5.0) << i; }
    else EXPECT_GT(d, 5.0) << i;
  }
}

// Priority: a droplet (priority 1) inside a slab (priority 0) keeps the slab's beads out of it; an empty void
// (a region no species names) stays empty
TEST(DpdRegions, PriorityAndVoid) {
  DpdOptions o;
  o.box = 12;
  o.regions = {{"slab", "slab", {0.5, 0.5, 0.5}, 0, 0, 6.0, 0, 0}, {"drop", "droplet", {0.5, 0.5, 0.5}, 2.0, 0, 3, 2, 1},
               {"void", "droplet", {0.5, 0.5, 0.05}, 1.5, 0, 3, 2, 2}};
  o.species = {{"slab", "A", 1000, "slab", ""}, {"drop", "B", 60, "drop", ""}, {"rest", "C", 2500, "", ""}};
  const auto r = start_of(o);
  const auto& t = r.frames;
  for (size_t i = 0; i < t.topology.atoms.size(); ++i) {
    const Vec3& x = t.positions[0][i];
    const int z = t.topology.atoms[i].element;
    const double d = rmin(x, {6, 6, 6}, 12);
    if (z == 6) { EXPECT_GT(d, 2.0); EXPECT_LE(std::fabs(x[0] - 6), 3.0); }
    if (z == 8) EXPECT_LE(d, 2.0);
    EXPECT_GT(rmin(x, {6, 6, 0.6}, 12), 1.5) << "a bead in the void";
  }
}

// Amphiphiles on a droplet: each molecule's head (first bead) on the surface at 4 r_c, its tail 0.7 r_c per bead inward
// along the radius
TEST(DpdRegions, AmphiphilesAlongTheNormal) {
  DpdOptions o;
  o.box = 14;
  o.regions = {{"drop", "droplet", {0.5, 0.5, 0.5}, 4.0, 0, 3, 2, 0}};
  o.species = {{"surfactant", "HTTT", 40, "drop", "in"}, {"rest", "C", 6000, "", ""}};
  // head H and tail T are not among A–D: they map to P in the frames; the molecule order gives each bead's place
  const auto r = start_of(o);
  const auto& P = r.frames.positions[0];
  const Vec3 c{7, 7, 7};
  for (int m = 0; m < 40; ++m) {
    const Vec3& head = P[size_t(4 * m)];
    EXPECT_NEAR(rmin(head, c, 14), 4.0, 1e-9);
    for (int k = 1; k < 4; ++k) EXPECT_NEAR(rmin(P[size_t(4 * m + k)], c, 14), 4.0 - 0.7 * k, 1e-9);
  }
}
