#include <gtest/gtest.h>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <sstream>

#include "caps/io.hpp"
#include "caps/table_tools.hpp"

using namespace caps;

// 1 … 100 by Freedman–Diaconis: quartiles 25.75 and 75.25 (linear interpolation), h = 2 · 49.5 / 100^(1/3) = 21.33,
// so ⌈99 / h⌉ = 5 bins of 19.8 with 20 values each; Sturges gives ⌈log₂ 100⌉ + 1 = 8 bins
TEST(TableTools, BinRules) {
  std::vector<double> v;
  for (int i = 1; i <= 100; ++i) v.push_back(i);
  const auto fd = histogram(v, "fd");
  ASSERT_EQ(fd.counts.size(), 5u);
  EXPECT_NEAR(fd.width, 19.8, 1e-12);
  for (int c : fd.counts) EXPECT_EQ(c, 20);
  EXPECT_EQ(histogram(v, "sturges").counts.size(), 8u);
  // Scott on 1 … 10: s = 3.0277, h = 3.49 s / 10^(1/3) = 4.90, two bins of 4.5 holding five values each
  std::vector<double> w;
  for (int i = 1; i <= 10; ++i) w.push_back(i);
  const auto sc = histogram(w, "scott");
  ASSERT_EQ(sc.counts.size(), 2u);
  EXPECT_EQ(sc.counts[0], 5);
  EXPECT_EQ(sc.counts[1], 5);
  // four bins over [0, 1]: the last closed on the right; a NaN left out and counted
  const auto b = histogram({0, 0.25, 0.5, 0.75, 1, NAN}, "bins", 4);
  EXPECT_EQ(b.counts, (std::vector<int>{1, 1, 1, 2}));
  EXPECT_EQ(b.left_out, 1);
  EXPECT_EQ(b.n, 5);
  // all equal: no width, one bin
  const auto e = histogram({3, 3, 3}, "fd");
  ASSERT_EQ(e.counts.size(), 1u);
  EXPECT_EQ(e.counts[0], 3);
}

// Ascending by the first key, ties descending by the second, the row with no value last
TEST(TableTools, MultiColumnSort) {
  const auto o = sort_rows({{2, 1, 2, NAN, 1}, {5, 3, 4, 1, 9}}, {false, true});
  EXPECT_EQ(o, (std::vector<size_t>{4, 1, 0, 2, 3}));
  const auto d = sort_rows({{2, 1, 2, NAN, 1}}, {true});
  EXPECT_EQ(d, (std::vector<size_t>{0, 2, 1, 4, 3}));   // stable among equals, not finite still last
}

// An SD file of the rows with a structure: one record, its data items as >  <name> / value / blank line, read back
TEST(TableTools, SdfWithDataItems) {
  const auto dir = std::filesystem::temp_directory_path() / "caps_table_tools";
  std::filesystem::create_directories(dir);
  const auto xyz = (dir / "methane.xyz").string(), sdf = (dir / "rows.sdf").string();
  {
    std::ofstream f(xyz);
    f << "5\nmethane\nC 0 0 0\nH 0.629 0.629 0.629\nH 0.629 -0.629 -0.629\nH -0.629 0.629 -0.629\nH -0.629 -0.629 0.629\n";
  }
  std::vector<std::string> notes;
  const int n = write_sdf_table({{"row 1", xyz, {{"Tg (K)", "350.5"}, {"name", "methane"}}}, {"row 2", "", {{"Tg (K)", "1"}}}}, sdf, &notes);
  EXPECT_EQ(n, 1);
  ASSERT_EQ(notes.size(), 1u);
  std::ifstream in(sdf);
  std::stringstream ss;
  ss << in.rdbuf();
  const std::string t = ss.str();
  EXPECT_EQ(t.rfind("row 1\n", 0), 0u);
  EXPECT_NE(t.find("M  END\n>  <Tg (K)>\n350.5\n\n>  <name>\nmethane\n\n$$$$\n"), std::string::npos);
  const Trajectory back = open_file(sdf);
  EXPECT_EQ(back.topology.atoms.size(), 5u);
  EXPECT_THROW(write_sdf_table({{"x", xyz, {{"a<b", "1"}}}}, sdf), std::invalid_argument);
}
