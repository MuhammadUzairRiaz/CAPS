#include <gtest/gtest.h>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <tuple>
#include <vector>

#include "caps/analysis.hpp"
#include "caps/import.hpp"
#include "caps/io.hpp"

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

// Materials Studio .car/.mdf: coordinates, the cell from a b c α β γ (a along x, b in the xy plane), the force-field type
// as each atom's name, its charge; bonds from the .mdf, one across the cell (%100) and one with an order (/1.0), each once.
TEST(Import, MaterialsStudioCarMdf) {
  const Trajectory t = open_file(std::string(CAPS_SOURCE_DIR) + "/tests/data/car/water_pbc.car");
  const System& s = t.topology;
  ASSERT_EQ(s.atoms.size(), 6u);
  EXPECT_EQ(s.source_format, "car");
  EXPECT_EQ(s.title, "Two SPC waters, one across the cell");
  EXPECT_EQ(s.atoms[0].name, "o*");
  EXPECT_EQ(s.atoms[1].name, "h*");
  EXPECT_EQ(s.atoms[0].element, 8);
  EXPECT_NEAR(s.atoms[3].charge, -0.82, 1e-12);
  EXPECT_TRUE(s.has_charges);
  EXPECT_NEAR(s.cell.a[0], 10.0, 1e-9);
  EXPECT_NEAR(s.cell.b[0], 12.0 * std::cos(120.0 * M_PI / 180), 1e-9);
  EXPECT_NEAR(s.cell.b[1], 12.0 * std::sin(120.0 * M_PI / 180), 1e-9);
  EXPECT_NEAR(s.cell.c[2], 14.0, 1e-9);
  // O-H twice per water, listed from both ends in the .mdf; #atomset's "@list subset" lines (WAT_1:H1 H2) are no atoms
  ASSERT_EQ(s.bonds.size(), 4u);
  EXPECT_TRUE(s.bonds_from_file);
  std::set<std::pair<uint32_t, uint32_t>> b;
  for (const auto& x : s.bonds) b.insert({x.i, x.j});
  EXPECT_TRUE(b.count({0, 1}) && b.count({0, 2}) && b.count({3, 4}) && b.count({3, 5}));
}

// .car/.mdf written back: the same atoms (types, charges), bonds and cell lengths and angles; atoms inside the cell
TEST(Import, MaterialsStudioCarRoundTrip) {
  const System a = open_file(std::string(CAPS_SOURCE_DIR) + "/tests/data/car/water_pbc.car").topology;
  const std::string out = (std::filesystem::temp_directory_path() / "caps_roundtrip.car").string();
  write_car(a, out);
  EXPECT_TRUE(std::filesystem::exists((std::filesystem::temp_directory_path() / "caps_roundtrip.mdf")));
  const System b = open_file(out).topology;
  ASSERT_EQ(b.atoms.size(), a.atoms.size());
  for (size_t i = 0; i < a.atoms.size(); ++i) {
    EXPECT_EQ(b.atoms[i].name, a.atoms[i].name);
    EXPECT_NEAR(b.atoms[i].charge, a.atoms[i].charge, 1e-4);
    const Vec3 f = b.cell.to_fractional(b.atoms[i].pos);
    for (int k = 0; k < 3; ++k) EXPECT_TRUE(f[k] >= -1e-9 && f[k] < 1 + 1e-9);
  }
  EXPECT_EQ(b.bonds.size(), a.bonds.size());
  EXPECT_NEAR(norm(b.cell.b), norm(a.cell.b), 1e-6);
  EXPECT_NEAR(dot(b.cell.a, b.cell.b), dot(a.cell.a, a.cell.b), 1e-4);
}

// Bond perception: a long contact between two atoms bonded to a common third closes a ring across it, it is no bond
// (Si···Si over an edge-sharing Si2O2 ring of amorphous silica, 2.6 Å); cyclopropane keeps its three C–C bonds.
TEST(Import, PerceptionLeavesRingContacts) {
  System s;
  const double y = std::sqrt(1.65 * 1.65 - 1.3 * 1.3);
  for (auto [z, x, yy] : std::vector<std::tuple<int, double, double>>{{14, -1.3, 0}, {14, 1.3, 0}, {8, 0, y}, {8, 0, -y}}) {
    Atom a;
    a.element = z;
    a.pos = {x, yy, 0};
    s.atoms.push_back(a);
  }
  const auto b = perceive_bonds(s);
  EXPECT_EQ(b.size(), 4u);
  for (const auto& x : b) EXPECT_FALSE(s.atoms[x.i].element == 14 && s.atoms[x.j].element == 14);
  System c;
  for (int k = 0; k < 3; ++k) {
    Atom a;
    a.element = 6;
    a.pos = {1.51 / std::sqrt(3.0) * std::cos(2 * M_PI * k / 3), 1.51 / std::sqrt(3.0) * std::sin(2 * M_PI * k / 3), 0};
    c.atoms.push_back(a);
  }
  EXPECT_EQ(perceive_bonds(c).size(), 3u);
}

namespace {
const std::string kTraj = std::string(CAPS_SOURCE_DIR) + "/tests/data/traj/";

// frames of a multi-frame .gro as written by gmx (nm → Å)
std::vector<std::vector<caps::Vec3>> gro_frames(const std::string& path) {
  std::ifstream in(path);
  std::vector<std::string> L;
  for (std::string l; std::getline(in, l);) L.push_back(l);
  std::vector<std::vector<caps::Vec3>> out;
  for (size_t k = 0; k + 1 < L.size();) {
    const size_t n = std::stoul(L[k + 1]);
    std::vector<caps::Vec3> f;
    for (size_t i = 0; i < n; ++i) {
      const auto& a = L[k + 2 + i];
      f.push_back({std::stod(a.substr(20, 8)) * 10, std::stod(a.substr(28, 8)) * 10, std::stod(a.substr(36, 8)) * 10});
    }
    out.push_back(f);
    k += n + 3;
  }
  return out;
}
}  // namespace

// GROMACS xtc (compressed, with water's run-length coding) and trr from gmx trjconv: every coordinate back to the .gro's
TEST(Import, XtcAndTrrMatchGromacs) {
  const auto ref = gro_frames(kTraj + "water3.gro");
  ASSERT_EQ(ref.size(), 3u);
  for (const char* ext : {"xtc", "trr"}) {
    const caps::Trajectory t = caps::open_file(kTraj + "water3." + ext, kTraj + "water3.gro");
    ASSERT_EQ(t.frames(), 3u) << ext;
    ASSERT_EQ(t.topology.atoms.size(), ref[0].size());
    EXPECT_NEAR(caps::norm(t.cells[0].a), 16.0, 1e-4);
    double worst = 0;
    for (size_t f = 0; f < 3; ++f)
      for (size_t i = 0; i < ref[f].size(); ++i) worst = std::max(worst, caps::norm(t.positions[f][i] - ref[f][i]));
    EXPECT_LT(worst, 1e-4) << ext;
    EXPECT_EQ(t.topology.atoms[0].element, 8);
  }
  EXPECT_THROW(caps::open_file(kTraj + "water3.xtc"), caps::ReadError);   // coordinates only: needs its structure
}

// DCD from LAMMPS (dump dcd, a tilted cell): coordinates and cell as LAMMPS's own text dump of the same steps
TEST(Import, DcdMatchesLammps) {
  const caps::Trajectory t = caps::open_file(kTraj + "water.dcd", kTraj + "water.data");
  const caps::Trajectory d = caps::open_file(kTraj + "water.lammpstrj", kTraj + "water.data");
  ASSERT_EQ(t.frames(), d.frames());
  ASSERT_EQ(t.frames(), 3u);
  EXPECT_NEAR(t.cells[0].volume(), d.cells[0].volume(), 1e-3);
  double worst = 0;
  for (size_t f = 0; f < t.frames(); ++f)
    for (size_t i = 0; i < t.topology.atoms.size(); ++i)
      worst = std::max(worst, caps::norm(t.cells[f].minimum_image(t.positions[f][i] - d.positions[f][i])));
  EXPECT_LT(worst, 1e-4);
}

// MDL molfile V2000 (charges from M  CHG) and V3000, and a VASP 5 POSCAR in direct coordinates
TEST(Import, SdfAndPoscar) {
  const auto dir = std::filesystem::temp_directory_path();
  {
    std::ofstream o(dir / "caps_acetate.sdf");
    o << "acetate\n  test\n\n  4  3  0  0  0  0  0  0  0  0999 V2000\n"
         "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
         "    1.5000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
         "    2.1000    1.1000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
         "    2.1000   -1.1000    0.0000 O   0  5  0  0  0  0  0  0  0  0  0  0\n"
         "  1  2  1  0\n  2  3  2  0\n  2  4  1  0\nM  CHG  1   4  -1\nM  END\n$$$$\n"
         "methane\n  test\n\n  0  0  0     0  0            999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 1 0 0 0 0\nM  V30 BEGIN ATOM\n"
         "M  V30 1 C 5.0 5.0 5.0 0\nM  V30 END ATOM\nM  V30 END CTAB\nM  END\n$$$$\n";
  }
  const caps::System s = caps::open_file((dir / "caps_acetate.sdf").string()).topology;
  ASSERT_EQ(s.atoms.size(), 5u);
  EXPECT_EQ(s.bonds.size(), 3u);
  EXPECT_EQ(s.bonds[1].order, 2);
  EXPECT_DOUBLE_EQ(s.atoms[3].charge, -1.0);
  EXPECT_EQ(s.atoms[4].mol, 2);
  EXPECT_DOUBLE_EQ(s.atoms[4].pos[2], 5.0);
  {
    std::ofstream o(dir / "POSCAR");
    o << "rutile\n1.0\n4.594 0 0\n0 4.594 0\n0 0 2.959\nTi O\n2 4\nDirect\n0 0 0\n0.5 0.5 0.5\n0.3048 0.3048 0\n0.6952 0.6952 0\n"
         "0.8048 0.1952 0.5\n0.1952 0.8048 0.5\n";
  }
  const caps::Trajectory p = caps::open_file((dir / "POSCAR").string());
  ASSERT_EQ(p.topology.atoms.size(), 6u);
  EXPECT_EQ(p.topology.atoms[0].element, 22);
  EXPECT_EQ(p.topology.atoms[5].element, 8);
  EXPECT_NEAR(p.topology.atoms[1].pos[2], 1.4795, 1e-9);
  EXPECT_NEAR(p.topology.cell.volume(), 4.594 * 4.594 * 2.959, 1e-9);
}
