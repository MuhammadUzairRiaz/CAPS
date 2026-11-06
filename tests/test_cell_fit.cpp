#include <gtest/gtest.h>

#include <filesystem>

#include "caps/cell_fit.hpp"
#include "caps/io.hpp"
#include "caps/spacegroup.hpp"

using namespace caps;

namespace {
System pe() {
  const auto p = std::filesystem::path(CAPS_SOURCE_DIR) / "data/crystals/polymers/pe-orthorhombic.cif";
  return open_file(p.string()).frame(0);
}
System two_atoms(double box) {
  System s;
  s.cell.a = {box, 0, 0}, s.cell.b = {0, box, 0}, s.cell.c = {0, 0, box};
  for (double x : {0.0, 1.5}) {
    Atom a;
    a.element = 6, a.pos = {x, 0, 0};
    s.atoms.push_back(a);
  }
  s.bonds.push_back({0, 1});
  return s;
}
}  // namespace

// The cells at which polyethylene's OPLS energy per cell reaches its converged value are the ones that fit (measured
// against 4×4×10): c at 4 repeats so 1-4 partners along the chain stop tying, b at 2 so hydrogens read their own chain
TEST(CellFit, PolyethyleneCells) {
  const System c = pe();
  for (const std::array<int, 3> n : {std::array<int, 3>{1, 1, 1}, {1, 1, 2}, {1, 1, 3}, {2, 2, 3}, {3, 3, 3}, {1, 1, 4}, {1, 1, 8}, {2, 1, 8}})
    EXPECT_FALSE(cell_fit(supercell(c, n[0], n[1], n[2])).fits) << n[0] << n[1] << n[2];
  for (const std::array<int, 3> n : {std::array<int, 3>{1, 2, 4}, {2, 2, 4}, {3, 3, 4}, {2, 2, 8}, {3, 3, 8}})
    EXPECT_TRUE(cell_fit(supercell(c, n[0], n[1], n[2])).fits) << n[0] << n[1] << n[2];
}

TEST(CellFit, CopiesToFit) {
  const auto n = copies_to_fit(pe());
  EXPECT_TRUE(cell_fit(supercell(pe(), n[0], n[1], n[2])).fits);
  EXPECT_EQ(n, (std::array<int, 3>{1, 2, 4}));   // whose energy per cell is already the 4×4×10 value
  EXPECT_NE(cell_fit_message(pe()).find("1×2×4"), std::string::npos) << cell_fit_message(pe());
}

// A molecule in a wide box, no bonds at all, and no periodic directions all pass; two images of one neighbour do not
TEST(CellFit, SimpleCases) {
  EXPECT_TRUE(cell_fit(two_atoms(20)).fits);
  System noperiod = two_atoms(2.8);
  noperiod.cell.periodic = {false, false, false};
  EXPECT_TRUE(cell_fit(noperiod).fits);
  System unbonded = two_atoms(5);
  unbonded.bonds.clear();
  EXPECT_TRUE(cell_fit(unbonded).fits);
  const CellFit f = cell_fit(two_atoms(2.8));   // 1.5 Å one way, 1.3 Å the other
  EXPECT_FALSE(f.fits);
  EXPECT_GE(f.double_images, 1u);
  EXPECT_TRUE(f.axis[0]);
}
