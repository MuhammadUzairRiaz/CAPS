#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <set>

#include "caps/crystal.hpp"
#include "caps/spacegroup.hpp"

using namespace caps;

namespace {

const std::string kCrystals = std::string(CAPS_SOURCE_DIR) + "/data/crystals/";

size_t order(const std::string& q) {
  const auto* sg = find_space_group(q);
  if (!sg) return 0;
  return hall_operations(sg->hall).size();
}

std::set<std::string> op_set(const std::vector<SymOp>& ops) {
  std::set<std::string> out;
  for (auto op : ops) {
    for (double& t : op.t) { t -= std::floor(t + 1e-9); if (std::fabs(t - 1) < 1e-9) t = 0; }
    out.insert(symop_string(op));
  }
  return out;
}

CrystalSite site(const std::string& label, int z, double x, double y, double w) {
  CrystalSite s;
  s.label = label, s.element = z, s.frac = {x, y, w};
  return s;
}

}  // namespace

TEST(SpaceGroup, TableHasAllSettings) {
  EXPECT_EQ(space_group_settings().size(), 530u);
  std::set<int> numbers;
  for (const auto& s : space_group_settings()) numbers.insert(s.number);
  EXPECT_EQ(numbers.size(), 230u);
}

TEST(SpaceGroup, LookupByKeyNumberAndSymbol) {
  ASSERT_NE(find_space_group("Pnma"), nullptr);
  EXPECT_EQ(find_space_group("Pnma")->number, 62);
  EXPECT_EQ(find_space_group("P 21/c")->number, 14);
  EXPECT_EQ(find_space_group("fd-3m")->number, 227);
  EXPECT_EQ(find_space_group("227")->key, "227:2");
  EXPECT_EQ(find_space_group("166")->key, "166:h");
  EXPECT_EQ(find_space_group("nonsense"), nullptr);
  EXPECT_EQ(crystal_system(14), "monoclinic");
  EXPECT_EQ(crystal_system(194), "hexagonal");
  EXPECT_EQ(crystal_system(166), "trigonal");
}

TEST(SpaceGroup, GroupOrdersFromHallSymbols) {
  EXPECT_EQ(order("P 1"), 1u);
  EXPECT_EQ(order("P -1"), 2u);
  EXPECT_EQ(order("P 21/c"), 4u);
  EXPECT_EQ(order("C 2/c"), 8u);
  EXPECT_EQ(order("Pnma"), 8u);
  EXPECT_EQ(order("P 42/m n m"), 16u);
  EXPECT_EQ(order("P 32 2 1"), 6u);
  EXPECT_EQ(order("R -3 m:h"), 36u);
  EXPECT_EQ(order("R -3 m:r"), 12u);
  EXPECT_EQ(order("P 63 m c"), 12u);
  EXPECT_EQ(order("P 63/m m c"), 24u);
  EXPECT_EQ(order("F m -3 m"), 192u);
  EXPECT_EQ(order("227:2"), 192u);
  EXPECT_EQ(order("227:1"), 192u);
  EXPECT_EQ(order("I a -3 d"), 96u);
  EXPECT_EQ(order("I m -3 m"), 96u);
  // every setting closes into a group whose order divides by its point group's
  for (const auto& s : space_group_settings()) {
    const auto ops = hall_operations(s.hall);
    ASSERT_FALSE(ops.empty()) << s.key;
    EXPECT_TRUE(ops.size() <= 192) << s.key;
  }
}

TEST(SpaceGroup, QuartzOperationsMatchTheCif) {
  const auto ops = hall_operations(find_space_group("P 32 2 1")->hall);
  std::vector<SymOp> cif;
  for (const char* s : {"x,y,z", "-y,x-y,z+2/3", "-x+y,-x,z+1/3", "y,x,-z", "x-y,-y,-z+1/3", "-x,-x+y,-z+2/3"}) cif.push_back(parse_symop(s));
  EXPECT_EQ(op_set(ops), op_set(cif));
}

TEST(SpaceGroup, SymopRoundTrip) {
  const SymOp op = parse_symop("-y+1/2, x-y, z+3/4");
  EXPECT_EQ(op_set({parse_symop(symop_string(op))}), op_set({op}));
}

TEST(SpaceGroup, BuildsTheCrystalsOfTheLibrary) {
  CrystalReport rep;
  CrystalSpec quartz;
  quartz.space_group = "P 32 2 1";
  quartz.a = quartz.b = 4.9134, quartz.c = 5.4052, quartz.gamma = 120;
  quartz.sites = {site("Si1", 14, 0.4697, 0, 2.0 / 3), site("O1", 8, 0.4135, 0.2669, 0.785767)};
  System q = build_crystal(quartz, &rep);
  EXPECT_EQ(q.atoms.size(), 9u);
  EXPECT_EQ(rep.number, 154);
  EXPECT_EQ(rep.system, "trigonal");
  EXPECT_EQ(rep.multiplicity, (std::vector<int>{3, 6}));
  EXPECT_NEAR(rep.volume, 113.0, 0.1);
  EXPECT_EQ(read_cif(kCrystals + "alpha-quartz.cif").atoms.size(), q.atoms.size());

  CrystalSpec rutile;
  rutile.space_group = "P 42/m n m";
  rutile.a = rutile.b = 4.5937, rutile.c = 2.9587;
  rutile.sites = {site("Ti1", 22, 0, 0, 0), site("O1", 8, 0.3048, 0.3048, 0)};
  EXPECT_EQ(build_crystal(rutile, &rep).atoms.size(), 6u);
  EXPECT_EQ(rep.multiplicity, (std::vector<int>{2, 4}));

  CrystalSpec diamond;
  diamond.space_group = "F d -3 m";
  diamond.a = diamond.b = diamond.c = 3.567;
  diamond.sites = {site("C1", 6, 0.125, 0.125, 0.125)};
  System d = build_crystal(diamond, &rep);
  EXPECT_EQ(d.atoms.size(), 8u);
  EXPECT_EQ(rep.multiplicity, (std::vector<int>{8}));
  EXPECT_EQ(d.bonds.size(), 16u);   // 8 atoms × 4 neighbours / 2 in the periodic cell

  CrystalSpec salt;
  salt.space_group = "F m -3 m";
  salt.a = salt.b = salt.c = 5.64;
  salt.sites = {site("Na1", 11, 0, 0, 0), site("Cl1", 17, 0.5, 0.5, 0.5)};
  salt.supercell = {2, 2, 2};
  EXPECT_EQ(build_crystal(salt, &rep).atoms.size(), 64u);
  EXPECT_EQ(rep.atoms_per_cell, 8u);

  CrystalSpec zno;
  zno.space_group = "P 63 m c";
  zno.a = zno.b = 3.25, zno.c = 5.207, zno.gamma = 120;
  zno.sites = {site("Zn1", 30, 1.0 / 3, 2.0 / 3, 0), site("O1", 8, 1.0 / 3, 2.0 / 3, 0.382)};
  EXPECT_EQ(build_crystal(zno, &rep).atoms.size(), 4u);
}

TEST(SpaceGroup, FindsTheSymmetryOfLibraryCrystals) {
  const auto f = [](const std::string& name) { return find_symmetry(read_cif(kCrystals + name), 0.05); };
  EXPECT_EQ(f("rock-salt.cif").number, 225);
  EXPECT_EQ(f("copper.cif").number, 225);
  EXPECT_EQ(f("diamond.cif").number, 227);
  EXPECT_EQ(f("alpha-iron.cif").number, 229);
  EXPECT_EQ(f("beta-brass.cif").number, 221);
  EXPECT_EQ(f("rutile.cif").number, 136);
  EXPECT_EQ(f("alpha-quartz.cif").number, 154);
  const auto salt = f("rock-salt.cif");
  EXPECT_EQ(salt.sites.size(), 2u);
  EXPECT_EQ(salt.operations, 192);
}

TEST(SpaceGroup, BuiltCrystalRoundTripsThroughFindSymmetry) {
  CrystalSpec rutile;
  rutile.space_group = "P 42/m n m";
  rutile.a = rutile.b = 4.5937, rutile.c = 2.9587;
  rutile.sites = {site("Ti1", 22, 0, 0, 0), site("O1", 8, 0.3048, 0.3048, 0)};
  const auto found = find_symmetry(build_crystal(rutile));
  EXPECT_EQ(found.number, 136);
  ASSERT_EQ(found.sites.size(), 2u);
  CrystalSpec again = rutile;
  again.space_group = found.key;
  again.sites = found.sites;
  EXPECT_EQ(build_crystal(again).atoms.size(), 6u);
}

TEST(SpaceGroup, PrimitiveCellsAndSupercells) {
  const System salt = read_cif(kCrystals + "rock-salt.cif");
  const System p = primitive_cell(salt, 'F');
  EXPECT_EQ(p.atoms.size(), salt.atoms.size() / 4);
  EXPECT_NEAR(std::fabs(p.cell.volume()), std::fabs(salt.cell.volume()) / 4, 1e-6);
  const System fe = read_cif(kCrystals + "alpha-iron.cif");
  EXPECT_EQ(primitive_cell(fe, 'I').atoms.size(), 1u);
  EXPECT_EQ(supercell(fe, 2, 3, 1).atoms.size(), fe.atoms.size() * 6);
  EXPECT_EQ(primitive_cell(fe, 'P').atoms.size(), fe.atoms.size());
}

TEST(SpaceGroup, SymmetrizeSnapsSitesOntoSpecialPositions) {
  CrystalSpec zno;
  zno.space_group = "P 63 m c";
  zno.a = zno.b = 3.25, zno.c = 5.207, zno.gamma = 120;
  zno.sites = {site("Zn1", 30, 0.33, 0.67, 0), site("O1", 8, 0.3335, 0.6668, 0.382)};
  CrystalReport rep;
  EXPECT_GT(build_crystal(zno, &rep).atoms.size(), 4u);   // near-misses of the special position give extra atoms
  int moved = 0;
  const CrystalSpec fixed = symmetrize_sites(zno, 0.3, &moved);
  EXPECT_EQ(moved, 2);
  EXPECT_NEAR(fixed.sites[0].frac[0], 1.0 / 3, 1e-9);
  EXPECT_NEAR(fixed.sites[0].frac[1], 2.0 / 3, 1e-9);
  EXPECT_NEAR(fixed.sites[1].frac[2], 0.382, 1e-9);   // z is free on 2b
  EXPECT_EQ(build_crystal(fixed).atoms.size(), 4u);
  // a general position stays where it is
  CrystalSpec q;
  q.space_group = "P 32 2 1";
  q.a = q.b = 4.9134, q.c = 5.4052, q.gamma = 120;
  q.sites = {site("O1", 8, 0.4135, 0.2669, 0.785767)};
  symmetrize_sites(q, 0.3, &moved);
  EXPECT_EQ(moved, 0);
}
