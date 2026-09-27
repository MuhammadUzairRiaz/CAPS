#include <gtest/gtest.h>

#include <cmath>
#include <filesystem>
#include <fstream>

#include "caps/elements.hpp"
#include "caps/io.hpp"
#include "caps/pack.hpp"
#include "caps/molecule.hpp"

using namespace caps;

namespace {

System water() {
  System w;
  auto add = [&](int el, Vec3 p) { Atom a; a.element = el; a.pos = p; a.id = int64_t(w.atoms.size() + 1); a.name = el == 8 ? "OW" : "HW"; w.atoms.push_back(a); };
  add(8, {0, 0, 0});
  add(1, {0.957, 0, 0});
  add(1, {-0.240, 0.927, 0});
  w.bonds = {{0, 1}, {0, 2}};
  return w;
}

Cell cube(double l) {
  Cell c;
  c.a = {l, 0, 0};
  c.b = {0, l, 0};
  c.c = {0, 0, l};
  return c;
}

Region box(Vec3 lo, Vec3 hi, Region::Kind k = Region::InsideBox) {
  Region r;
  r.kind = k;
  r.a = lo;
  r.b = hi;
  return r;
}

}  // namespace

TEST(Pack, RegionViolations) {
  const Region b = box({0, 0, 0}, {10, 10, 10});
  EXPECT_EQ(b.violation({5, 5, 5}), 0);
  EXPECT_NEAR(b.violation({12, 5, 5}), 2, 1e-12);
  const Region ob = box({0, 0, 0}, {10, 10, 10}, Region::OutsideBox);
  EXPECT_NEAR(ob.violation({5, 5, 1}), 1, 1e-12);   // one Å from the nearest face
  EXPECT_EQ(ob.violation({-1, 5, 5}), 0);
  Region s;
  s.kind = Region::InsideSphere;
  s.a = {0, 0, 0};
  s.r = 5;
  EXPECT_NEAR(s.violation({8, 0, 0}), 3, 1e-12);
  s.kind = Region::OutsideSphere;
  EXPECT_NEAR(s.violation({3, 0, 0}), 2, 1e-12);
  Region c;
  c.kind = Region::InsideCylinder;
  c.a = {0, 0, 0};
  c.b = {0, 0, 1};
  c.r = 2;
  c.length = 10;
  EXPECT_EQ(c.violation({1, 1, 5}), 0);
  EXPECT_NEAR(c.violation({4, 0, 5}), 2, 1e-12);
  EXPECT_NEAR(c.violation({0, 0, 12}), 2, 1e-12);
  Region p;
  p.kind = Region::OverPlane;
  p.a = {0, 0, 1};
  p.r = 3;
  EXPECT_NEAR(p.violation({0, 0, 1}), 2, 1e-12);
  EXPECT_EQ(p.violation({0, 0, 4}), 0);
  p.kind = Region::BelowPlane;
  EXPECT_NEAR(p.violation({0, 0, 4}), 1, 1e-12);
}

TEST(Pack, PeriodicWaterBoxMeetsTolerance) {
  PackItem w{"water", water(), 400, {box({0, 0, 0}, {23, 23, 23})}};
  PackOptions o;
  o.cell = cube(23);
  o.seed = 3;
  PackReport r;
  const System s = pack({w}, o, &r);
  EXPECT_TRUE(r.success);
  ASSERT_EQ(s.atoms.size(), 1200u);
  EXPECT_EQ(s.bonds.size(), 800u);
  const auto [dmin, close] = intermolecular_contacts(s, 2.0, true);
  EXPECT_GE(dmin, 2.0);
  EXPECT_EQ(close, 0);
  // rigid: every O–H stays 0.957 Å
  for (size_t m = 0; m < 400; ++m) EXPECT_NEAR(norm(s.atoms[3 * m + 1].pos - s.atoms[3 * m].pos), 0.957, 1e-9);
  int nm = 0;
  s.molecules(&nm);
  EXPECT_EQ(nm, 400);
}

TEST(Pack, RegionsAndFixedMolecules) {
  // A fixed "wall" molecule, water inside a sphere, water outside the sphere but inside a box (not periodic).
  System wall;
  for (int i = 0; i < 5; ++i)
    for (int j = 0; j < 5; ++j) {
      Atom a;
      a.element = 6;
      a.pos = {10 + 3.0 * i, 10 + 3.0 * j, 2};
      a.id = int64_t(wall.atoms.size() + 1);
      wall.atoms.push_back(a);
    }
  PackItem fixed{"wall", wall, 1};
  fixed.fixed = true;
  Region sph;
  sph.kind = Region::InsideSphere;
  sph.a = {15, 15, 15};
  sph.r = 7;
  Region out = sph;
  out.kind = Region::OutsideSphere;
  out.r = 9;
  PackItem inner{"inner", water(), 40, {sph}};
  PackItem outer{"outer", water(), 150, {box({0, 0, 0}, {30, 30, 30}), out}};
  PackOptions o;
  o.periodic = false;
  o.seed = 5;
  PackReport r;
  const System s = pack({fixed, inner, outer}, o, &r);
  EXPECT_TRUE(r.success);
  EXPECT_EQ(r.region_violation, 0);
  for (int i = 0; i < 25; ++i) EXPECT_EQ(s.atoms[i].pos, wall.atoms[i].pos);   // fixed molecule untouched
  for (int i = 25; i < 25 + 120; ++i) EXPECT_LE(norm(s.atoms[i].pos - sph.a), 7.0);
  for (size_t i = 25 + 120; i < s.atoms.size(); ++i) {
    EXPECT_GE(norm(s.atoms[i].pos - sph.a), 9.0);
    for (int k = 0; k < 3; ++k) { EXPECT_GE(s.atoms[i].pos[k], 0.0); EXPECT_LE(s.atoms[i].pos[k], 30.0); }
  }
  EXPECT_GE(intermolecular_contacts(s, 2.0, false).first, 2.0);
}

TEST(Pack, ImpossibleRequestFailsLoudly) {
  PackItem w{"water", water(), 600, {box({0, 0, 0}, {15, 15, 15})}};
  PackOptions o;
  o.cell = cube(15);
  o.max_loops = 4;
  PackReport r;
  EXPECT_THROW(pack({w}, o, &r), PackError);
  EXPECT_FALSE(r.success);
  EXPECT_GT(r.close_pairs, 0);
}

TEST(Pack, SameSeedSamePacking) {
  PackItem w{"water", water(), 150, {box({0, 0, 0}, {18, 18, 18})}};
  PackOptions o;
  o.cell = cube(18);
  o.seed = 9;
  const System a = pack({w}, o), b = pack({w}, o);
  for (size_t i = 0; i < a.atoms.size(); ++i) ASSERT_EQ(a.atoms[i].pos, b.atoms[i].pos);
}

TEST(Pack, ReadsPackmolInput) {
  const auto dir = std::filesystem::temp_directory_path() / "caps_pack_inp";
  std::filesystem::create_directories(dir);
  write_xyz(water(), (dir / "w.xyz").string());
  {
    std::ofstream f(dir / "in.inp");
    f << "tolerance 2.5\nseed 4\nfiletype xyz\noutput out.xyz\npbc 0 0 0 20 20 30\n"
         "structure w.xyz\n  number 12\n  inside cube 0. 0. 0. 20.\n  below plane 0 0 2 30  # z < 15\nend structure\n"
         "structure w.xyz\n  number 1\n  fixed 5 5 25 0 0 0\n  center\nend structure\n";
  }
  PackOptions o;
  std::string out;
  const auto items = read_packmol_input((dir / "in.inp").string(), o, &out);
  ASSERT_EQ(items.size(), 2u);
  EXPECT_EQ(o.tolerance, 2.5);
  EXPECT_EQ(o.seed, 4u);
  EXPECT_TRUE(o.periodic);
  EXPECT_NEAR(o.cell.c[2], 30, 1e-12);
  EXPECT_EQ(out, (dir / "out.xyz").string());
  EXPECT_EQ(items[0].count, 12);
  ASSERT_EQ(items[0].regions.size(), 2u);
  EXPECT_EQ(items[0].regions[0].b, (Vec3{20, 20, 20}));
  EXPECT_NEAR(items[0].regions[1].r, 15, 1e-12);   // normalised plane
  EXPECT_TRUE(items[1].fixed && items[1].center);
  const System s = pack(items, o);
  for (size_t i = 0; i < 36; ++i) EXPECT_LE(s.atoms[i].pos[2], 15.0);
  {
    std::ofstream f(dir / "bad.inp");
    f << "structure w.xyz\n  number 3\n  inside ellipsoid 0 0 0 1 1 1 2\nend structure\n";
  }
  EXPECT_THROW(read_packmol_input((dir / "bad.inp").string(), o, &out), PackError);
}

// Stage 3: hexane packed loosely, then compressed to 0.6 g/cm³ with push-off: the density reached, molecules whole
// and no atoms of different molecules closer than 1.5 Å
TEST(Pack, CompressionToTargetDensity) {
  BuildOptions bo;
  bo.forcefield = "uff";
  const System hex = build_molecule("CCCCCC", bo).system;
  const auto dir = std::filesystem::temp_directory_path();
  write_xyz(hex, (dir / "caps_hexane.xyz").string());
  PackOptions o;
  const auto items = parse_packmol_input("tolerance 2.0\nseed 3\npbc 30 30 30\ncompress 0.6\nstructure caps_hexane.xyz\n  number 40\n  inside box 0 0 0 30 30 30\nend structure\n",
                                         dir.string(), o, nullptr);
  EXPECT_DOUBLE_EQ(o.compress_to, 0.6);
  PackReport r;
  const System s = pack(items, o, &r);
  EXPECT_NEAR(s.density(), 0.6, 1e-6);
  EXPECT_TRUE(std::any_of(r.notes.begin(), r.notes.end(), [](const std::string& n) { return n.rfind("compressed from", 0) == 0; }));
  double dmin = 1e9;
  for (size_t i = 0; i < s.atoms.size(); ++i)
    for (size_t j = i + 1; j < s.atoms.size(); ++j)
      if (s.atoms[i].mol != s.atoms[j].mol) dmin = std::min(dmin, norm(s.cell.minimum_image(s.atoms[i].pos - s.atoms[j].pos)));
  EXPECT_GT(dmin, 1.5);
}
