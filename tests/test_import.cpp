#include <gtest/gtest.h>

#include <cmath>
#include <filesystem>
#include <fstream>

#include "caps/import.hpp"

using namespace caps;

namespace {

// Benzene and ethane without bonds (plain XYZ); with a lattice, the ethane straddles the cell's x face.
std::string write_xyz(const std::string& name, bool lattice) {
  const auto path = (std::filesystem::temp_directory_path() / name).string();
  std::ofstream out(path);
  out << 20 << "\n";
  out << (lattice ? "Lattice=\"20 0 0 0 20 0 0 0 20\" Properties=species:S:1:pos:R:3" : "benzene and ethane") << "\n";
  for (int k = 0; k < 6; ++k) {
    const double t = k * M_PI / 3;
    out << "C " << 10 + 1.39 * std::cos(t) << " " << 10 + 1.39 * std::sin(t) << " 5\n";
  }
  for (int k = 0; k < 6; ++k) {
    const double t = k * M_PI / 3;
    out << "H " << 10 + 2.47 * std::cos(t) << " " << 10 + 2.47 * std::sin(t) << " 5\n";
  }
  // ethane along x: C at 19.3 and (19.3 + 1.54) − 20 = 0.84 when wrapped
  const double c1 = 19.3, c2 = lattice ? 0.84 : 20.84;
  out << "C " << c1 << " 10 15\nC " << c2 << " 10 15\n";
  for (int s : {-1, 1}) out << "H " << c1 - 0.36 << " " << 10 + s * 0.89 << " " << 15 + 0.51 << "\n";
  out << "H " << c1 - 0.36 << " 10 " << 15 - 1.03 << "\n";
  for (int s : {-1, 1}) out << "H " << c2 + 0.36 << " " << 10 + s * 0.89 << " " << 15 - 0.51 << "\n";
  out << "H " << c2 + 0.36 << " 10 " << 15 + 1.03 << "\n";
  return path;
}

}  // namespace

TEST(Import, PerceivesBondsOrdersAndMolecules) {
  const auto path = write_xyz("caps_import_plain.xyz", false);
  const Trajectory t = import_file(path, "", {});
  EXPECT_EQ(t.topology.bonds.size(), 19u);   // 6 + 6 in benzene, 7 in ethane
  int aromatic = 0, single = 0;
  for (const auto& b : t.topology.bonds) { aromatic += b.order == 4; single += b.order == 1; }
  EXPECT_EQ(aromatic, 6);
  EXPECT_EQ(single, 13);
  int n = 0;
  t.topology.molecules(&n);
  EXPECT_EQ(n, 2);
  EXPECT_TRUE(t.topology.has_mol);
  EXPECT_EQ(t.topology.atoms[0].mol, 1);
  EXPECT_EQ(t.topology.atoms[12].mol, 2);
}

TEST(Import, NoBondsAndTightTolerance) {
  const auto path = write_xyz("caps_import_plain2.xyz", false);
  ImportOptions none;
  none.bonds = ImportOptions::None;
  EXPECT_TRUE(import_file(path, "", none).topology.bonds.empty());
  ImportOptions tight;
  tight.tolerance = 0.0;   // with bare covalent radii only benzene's 1.39 Å ring bonds are short enough (C–H 1.08 > 1.07,
                           // ethane's C–C 1.54 > 1.52)
  const auto t = import_file(path, "", tight);
  EXPECT_EQ(t.topology.bonds.size(), 6u);
  ImportOptions from_file;
  from_file.bonds = ImportOptions::FromFile;
  EXPECT_TRUE(import_file(path, "", from_file).topology.bonds.empty());
}

TEST(Import, UnwrapsAcrossTheCellAndCanDropIt) {
  const auto path = write_xyz("caps_import_cell.xyz", true);
  const Trajectory t = import_file(path, "", {});
  ASSERT_TRUE(t.topology.cell.valid());
  EXPECT_EQ(t.topology.bonds.size(), 19u);
  const Vec3 d = t.topology.atoms[13].pos - t.topology.atoms[12].pos;
  EXPECT_NEAR(std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]), 1.54, 0.01);   // whole, not 18.46 Å apart
  EXPECT_TRUE(t.topology.unwrapped);

  ImportOptions nocell;
  nocell.use_cell = false;
  const Trajectory u = import_file(path, "", nocell);
  EXPECT_FALSE(u.topology.cell.valid());
  EXPECT_LT(u.topology.bonds.size(), 19u);   // the C–C bond across the face is gone
}

TEST(Import, PreviewHasAFragmentAndCounts) {
  const auto path = write_xyz("caps_import_prev.xyz", true);
  const ImportPreview p = import_preview(path, {});
  EXPECT_EQ(p.atoms, 20u);
  EXPECT_EQ(p.bonds_in_file, 0u);
  EXPECT_EQ(p.bonds, 19u);
  EXPECT_EQ(p.molecules, 2u);
  EXPECT_EQ(p.aromatic, 6);
  EXPECT_EQ(p.cell, "20 × 20 × 20 Å, orthogonal");
  EXPECT_EQ(p.fragment_heavy, 6);            // benzene's ring (the first carbon's component)
  EXPECT_EQ(p.fragment.atoms.size(), 12u);
  EXPECT_EQ(p.fragment.bonds.size(), 12u);
  EXPECT_FALSE(p.file.head.empty());
}
