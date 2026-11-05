#include <gtest/gtest.h>

#include <cmath>
#include <map>
#include <string>

#include "caps/crystal.hpp"

using namespace caps;

namespace {
std::string cif(const std::string& group, const std::string& rows, bool disorder_columns) {
  std::string t = "data_test\n_cell_length_a 10\n_cell_length_b 10\n_cell_length_c 10\n_cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90\n"
                  "_symmetry_space_group_name_H-M '" + group + "'\nloop_\n_symmetry_equiv_pos_as_xyz\n'x,y,z'\n";
  if (group == "P -1") t += "'-x,-y,-z'\n";
  t += "loop_\n_atom_site_label\n_atom_site_type_symbol\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n_atom_site_occupancy\n";
  if (disorder_columns) t += "_atom_site_disorder_assembly\n_atom_site_disorder_group\n";
  return t + rows;
}
int count(const System& s, const std::string& label) {
  int n = 0;
  for (const auto& a : s.atoms) n += a.name == label;
  return n;
}
}  // namespace

// Disorder groups: the major group by default, a named group on request, and a random supercell in which every copy
// of the assembly keeps one whole group, the groups in the ratio of their occupancies
TEST(CifDisorder, GroupsMajorNamedAndRandom) {
  const std::string rows = "O1 O 0 0 0 1 . .\n"
                           "C1A C 0.5 0.5 0.5 0.7 A 1\nC2A C 0.6 0.5 0.5 0.7 A 1\n"
                           "C1B C 0.5 0.6 0.5 0.3 A 2\nC2B C 0.5 0.7 0.5 0.3 A 2\n";
  const std::string text = cif("P 1", rows, true);
  CifDisorderReport rep;
  const System major = parse_cif(text, "t", CifDisorder{}, &rep);
  EXPECT_EQ(major.atoms.size(), 3u);
  EXPECT_EQ(count(major, "C1A"), 1);
  EXPECT_EQ(count(major, "C1B"), 0);
  EXPECT_EQ(rep.assemblies, 1);
  EXPECT_EQ(rep.groups, 2);
  CifDisorder g;
  g.mode = "group", g.group = "2";
  const System minor = parse_cif(text, "t", g);
  EXPECT_EQ(count(minor, "C1B"), 1);
  EXPECT_EQ(count(minor, "C2B"), 1);
  EXPECT_EQ(count(minor, "C1A"), 0);
  CifDisorder r;
  r.mode = "random", r.supercell = {10, 10, 10}, r.seed = 4;
  const System rnd = parse_cif(text, "t", r, &rep);
  EXPECT_EQ(count(rnd, "O1"), 1000);
  // per cell: both atoms of one group, never a mixture
  std::map<int, std::pair<int, int>> per;   // cell → (group 1 atoms, group 2 atoms)
  for (const auto& a : rnd.atoms) {
    if (a.element != 6) continue;
    const int c = int(a.pos[0] / 10) * 100 + int(a.pos[1] / 10) * 10 + int(a.pos[2] / 10);
    (a.name.back() == 'A' ? per[c].first : per[c].second) += 1;
  }
  int g1 = 0;
  for (const auto& [c, n] : per) {
    EXPECT_TRUE((n.first == 2 && n.second == 0) || (n.first == 0 && n.second == 2)) << "cell " << c;
    g1 += n.first == 2;
  }
  EXPECT_EQ(per.size(), 1000u);
  EXPECT_NEAR(g1 / 1000.0, 0.7, 3 * std::sqrt(0.21 / 1000));
  std::printf("random fill: group 1 in %d of 1000 cells (0.7 expected)\n", g1);
}

// With symmetry: a disordered pair in P-1 whose atoms fall on both sides of the cell boundary stays whole in every
// copy of the random supercell (each atom has its partner 1.2 Å away)
TEST(CifDisorder, CopiesStayWholeAcrossTheCellEdge) {
  const std::string rows = "C1A C 0.10 0.30 0.30 0.6 A 1\nC2A C -0.02 0.30 0.30 0.6 A 1\n"
                           "C1B C 0.10 0.40 0.30 0.4 A 2\nC2B C -0.02 0.40 0.30 0.4 A 2\n";
  CifDisorder r;
  r.mode = "random", r.supercell = {6, 6, 6}, r.seed = 9;
  const System s = parse_cif(cif("P -1", rows, true), "t", r);
  ASSERT_FALSE(s.atoms.empty());
  for (size_t i = 0; i < s.atoms.size(); ++i) {
    const std::string partner = s.atoms[i].name == "C1A" ? "C2A" : s.atoms[i].name == "C2A" ? "C1A" : s.atoms[i].name == "C1B" ? "C2B" : "C1B";
    double best = 1e9;
    for (size_t j = 0; j < s.atoms.size(); ++j)
      if (s.atoms[j].name == partner) best = std::min(best, norm(s.cell.minimum_image(s.atoms[j].pos - s.atoms[i].pos)));
    EXPECT_NEAR(best, 1.2, 1e-6) << s.atoms[i].name << " " << i;
  }
  EXPECT_EQ(s.atoms.size(), size_t(2 * 2 * 216));   // two molecule copies (the inversion) per cell, two atoms each
}

// Places without groups: a mixed Mg/Fe site keeps Mg by default and is filled in proportion at random, one atom per
// site; a lone 0.25 site is left out by default and filled a quarter of the time; two half-occupied split positions
// 0.3 Å apart give one atom, not two overlapping ones
TEST(CifDisorder, MixedSitesVacanciesAndSplitPositions) {
  const std::string mixed = cif("P 1", "M1 Mg 0 0 0 0.6\nM2 Fe 0 0 0 0.4\nO1 O 0.5 0.5 0.5 1\nNa1 Na 0.25 0.25 0.25 0.25\n"
                                           "C1 C 0.3 0.8 0.3 0.5\nC2 C 0.33 0.8 0.3 0.5\n", false);
  CifDisorderReport rep;
  const System major = parse_cif(mixed, "t", CifDisorder{}, &rep);
  EXPECT_EQ(count(major, "M1"), 1);
  EXPECT_EQ(count(major, "M2"), 0);
  EXPECT_EQ(count(major, "Na1"), 0);
  EXPECT_EQ(count(major, "C1") + count(major, "C2"), 1);
  EXPECT_EQ(major.atoms.size(), 3u);
  CifDisorder r;
  r.mode = "random", r.supercell = {10, 10, 10}, r.seed = 2;
  const System s = parse_cif(mixed, "t", r, &rep);
  const int mg = count(s, "M1"), fe = count(s, "M2"), na = count(s, "Na1"), c = count(s, "C1") + count(s, "C2");
  EXPECT_EQ(mg + fe, 1000);   // every metal site filled once
  EXPECT_NEAR(fe / 1000.0, 0.4, 3 * std::sqrt(0.24 / 1000));
  EXPECT_NEAR(na / 1000.0, 0.25, 3 * std::sqrt(0.1875 / 1000));
  EXPECT_EQ(c, 1000);         // the split pair: one carbon per cell
  EXPECT_NEAR(rep.expected[26], 400.0, 1e-6);
  EXPECT_EQ(rep.obtained[26], fe);
  std::printf("random fill: Fe %d / 1000 (400), Na %d / 1000 (250)\n", fe, na);
}
