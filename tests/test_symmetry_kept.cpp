#include <gtest/gtest.h>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <sstream>

#include "caps/crystal.hpp"
#include "caps/io.hpp"
#include "caps/spacegroup.hpp"

using namespace caps;

namespace {
std::string tmp(const std::string& name) { return (std::filesystem::temp_directory_path() / name).string(); }
std::string slurp(const std::string& p) {
  std::ifstream f(p);
  std::stringstream b;
  b << f.rdbuf();
  return b.str();
}
// every atom of a has an atom of the same element at the same place (modulo the cell) in b, and the counts agree
bool same_atoms(const System& a, const System& b, double tol = 1e-4) {
  if (a.atoms.size() != b.atoms.size()) return false;
  for (const auto& x : a.atoms) {
    bool hit = false;
    for (const auto& y : b.atoms) {
      if (x.element != y.element) continue;
      Vec3 d = a.cell.to_fractional(x.pos) - b.cell.to_fractional(y.pos);
      for (auto& v : d) v -= std::round(v);
      hit = hit || norm(a.cell.to_cartesian(d) - a.cell.to_cartesian({0, 0, 0})) < tol;
    }
    if (!hit) return false;
  }
  return true;
}
size_t count_lines(const std::string& text, const std::string& prefix) {
  std::istringstream in(text);
  size_t n = 0;
  for (std::string l; std::getline(in, l);) n += l.rfind(prefix, 0) == 0;
  return n;
}
System rock_salt() {
  CrystalSpec spec;
  spec.space_group = "Fm-3m";
  spec.a = spec.b = spec.c = 5.64;
  spec.sites = {{"Na", 11, {0, 0, 0}, 1}, {"Cl", 17, {0.5, 0.5, 0.5}, 1}};
  return build_crystal(spec);
}
System quartz() { return read_cif(std::string(CAPS_SOURCE_DIR) + "/data/crystals/alpha-quartz.cif"); }
}  // namespace

// Rock salt built in F m -3 m keeps its 192 operations; the CIF lists them with the two sites, and reading it back gives
// the same eight atoms and the same group
TEST(SymmetryKept, BuiltCrystalWritesItsAsymmetricUnit) {
  const System s = rock_salt();
  ASSERT_EQ(s.symmetry.ops.size(), 192u);
  EXPECT_EQ(s.symmetry.number, 225);
  EXPECT_TRUE(symmetry_holds(s, s.symmetry));
  EXPECT_EQ(asymmetric_unit(s, s.symmetry).size(), 2u);
  const std::string p = tmp("caps_symkept_nacl.cif");
  write_cif(s, p);
  const std::string text = slurp(p);
  EXPECT_NE(text.find("_space_group_IT_number 225"), std::string::npos);
  EXPECT_EQ(count_lines(text, "Na1 ") + count_lines(text, "Cl1 "), 2u);
  EXPECT_EQ(count_lines(text, "Na2 "), 0u);
  const System r = read_cif(p);
  EXPECT_TRUE(same_atoms(s, r));
  EXPECT_EQ(r.symmetry.ops.size(), 192u);
  EXPECT_EQ(r.symmetry.number, 225);
}

// Quartz read from a CIF keeps P 32 2 1 (6 operations) and writes back as its two sites; with the group stripped it is
// found again, in the structure's own coordinates; an edit that breaks it writes P 1
TEST(SymmetryKept, ReadFoundAndWritten) {
  const System q = quartz();
  ASSERT_EQ(q.symmetry.ops.size(), 6u);
  EXPECT_TRUE(symmetry_holds(q, q.symmetry));
  EXPECT_EQ(asymmetric_unit(q, q.symmetry).size(), 2u);
  const std::string p = tmp("caps_symkept_quartz.cif");
  write_cif(q, p);
  EXPECT_TRUE(same_atoms(q, read_cif(p)));
  System bare = q;
  bare.symmetry = {};
  const Symmetry g = symmetry_of(bare, 0.05);
  EXPECT_EQ(g.ops.size(), 6u);
  EXPECT_TRUE(symmetry_holds(bare, g, 0.05));
  System broken = q;
  broken.atoms[0].pos = broken.atoms[0].pos + Vec3{0.3, 0, 0};
  EXPECT_FALSE(symmetry_holds(broken, broken.symmetry));
  const std::string p1 = tmp("caps_symkept_broken.cif");
  write_cif(broken, p1);
  EXPECT_NE(slurp(p1).find("_space_group_IT_number 1\n"), std::string::npos);
}

// The cell follows the crystal system: cubic from a = 5, b = 6, c = 7, β = 95° becomes a = b = c = √((25 + 36 + 49)/3),
// all angles 90°; hexagonal quartz from γ = 118° becomes γ = 120° with a = b; the atoms keep their fractional coordinates
TEST(SymmetryKept, CellByCrystalSystem) {
  System s = rock_salt();
  set_cell_symmetric(s, 5, 6, 7, 90, 95, 90);
  const double a = std::sqrt((25.0 + 36 + 49) / 3);
  EXPECT_NEAR(norm(s.cell.a), a, 1e-9);
  EXPECT_NEAR(norm(s.cell.b), a, 1e-9);
  EXPECT_NEAR(norm(s.cell.c), a, 1e-9);
  EXPECT_NEAR(dot(s.cell.a, s.cell.c), 0, 1e-9);
  EXPECT_NEAR(dot(s.cell.b, s.cell.c), 0, 1e-9);
  EXPECT_NEAR(dot(s.cell.a, s.cell.b), 0, 1e-9);
  EXPECT_TRUE(symmetry_holds(s, s.symmetry, 1e-6));

  System q = quartz();
  const Vec3 f0 = q.cell.to_fractional(q.atoms[3].pos);
  set_cell_symmetric(q, 4.9, 5.0, 5.4, 91, 89, 118);
  EXPECT_NEAR(norm(q.cell.a), norm(q.cell.b), 1e-9);
  EXPECT_NEAR(std::acos(dot(q.cell.a, q.cell.b) / (norm(q.cell.a) * norm(q.cell.b))) * 180 / M_PI, 120, 1e-9);
  EXPECT_NEAR(dot(q.cell.a, q.cell.c), 0, 1e-9);
  EXPECT_NEAR(dot(q.cell.b, q.cell.c), 0, 1e-9);
  EXPECT_NEAR(norm(q.cell.to_fractional(q.atoms[3].pos) - f0), 0, 1e-12);
  EXPECT_TRUE(symmetry_holds(q, q.symmetry, 1e-5));
}

// (the quartz file gives z = 0.666667 for 2/3: its images agree to 4·10⁻⁶ Å)
// Moves mirrored: an oxygen of quartz (general position, 6 images) moved by a vector takes its five images with it and
// lands where asked; a silicon (site symmetry 2) moves along its twofold axis only, its three images with it; the group
// holds after both
TEST(SymmetryKept, MovesMirroredToImages) {
  System q = quartz();
  size_t o = 0, si = 0;
  for (size_t i = 0; i < q.atoms.size(); ++i) (q.atoms[i].element == 8 ? o : si) = i;
  const Vec3 to = q.atoms[o].pos + Vec3{0.05, 0.02, -0.03};
  EXPECT_EQ(move_atom_symmetric(q, o, to), 6u);
  EXPECT_NEAR(norm(q.atoms[o].pos - to), 0, 1e-9);
  EXPECT_TRUE(symmetry_holds(q, q.symmetry, 1e-5));
  const Vec3 before = q.atoms[si].pos;
  EXPECT_EQ(move_atom_symmetric(q, si, before + Vec3{0.04, 0.03, 0.02}), 3u);
  EXPECT_TRUE(symmetry_holds(q, q.symmetry, 1e-5));
  EXPECT_GT(norm(q.atoms[si].pos - before), 1e-4);   // the allowed part of the move
}
