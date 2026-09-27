#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <string>

#include "caps/dlpoly.hpp"
#include "caps/polymer.hpp"
#include "caps/uff.hpp"

using namespace caps;

// A small PE cell: one molecular type per distinct chain, every section counted, each 1-4 pair scaled exactly once,
// the atoms in CONFIG inside the cell about its centre
TEST(Dlpoly, FieldConfigControl) {
  ChainSpec c;
  c.units.push_back({"*CC*", "*CC*"});
  c.dp = 6;
  GrowOptions g;
  g.chains = 3;
  g.density = 0.3;
  g.seed = 5;
  const System s = grow_chains(c, g);
  const ForceField ff = default_forcefield(s);
  const auto dir = std::filesystem::temp_directory_path() / "caps_dlpoly_test";
  std::filesystem::remove_all(dir);
  const auto notes = write_dlpoly(s, ff, dir.string());
  ASSERT_FALSE(notes.empty());
  std::ifstream f(dir / "FIELD");
  std::string l;
  std::map<std::string, int> head;
  std::set<std::pair<int, int>> scaled;
  int mols = 0, nummols = 0, dih_scaled = 0;
  std::string section;
  while (std::getline(f, l)) {
    std::istringstream w(l);
    std::string k;
    w >> k;
    if (k == "molecule") ++mols;
    if (k == "nummols") { int x; w >> x; nummols += x; }
    if (k == "bonds" || k == "angles" || k == "dihedrals" || k == "atoms" || k == "vdw") { int x; w >> x; head[k] += x; section = k; continue; }
    if (section == "dihedrals" && (k == "cos" || k == "cos3")) {
      int a, b, cc, d;
      w >> a >> b >> cc >> d;
      double p[5] = {0, 0, 0, 0, 0};
      for (double& x : p) w >> x;
      const double e14 = p[3], v14 = p[4];   // cos3: A1 A2 A3 then the scales; cos: A δ m then the scales
      if (e14 != 0 || v14 != 0) {
        ++dih_scaled;
        EXPECT_TRUE(scaled.insert({std::min(a, d), std::max(a, d)}).second) << l;   // once per 1-4 pair
        EXPECT_NEAR(e14, ff.coul14, 1e-7);
        EXPECT_NEAR(v14, ff.lj14, 1e-7);
      }
    }
  }
  EXPECT_EQ(nummols, 3);
  EXPECT_GE(mols, 1);
  EXPECT_GT(head["bonds"], 0);
  EXPECT_GT(head["vdw"], 0);
  EXPECT_GT(dih_scaled, 0);
  // CONFIG: the header and one position per atom, inside the cell about its centre
  std::ifstream cf(dir / "CONFIG");
  std::getline(cf, l);
  int lev = -1, imcon = -1;
  size_t natms = 0;
  cf >> lev >> imcon >> natms;
  EXPECT_EQ(lev, 0);
  EXPECT_EQ(imcon, 1);
  EXPECT_EQ(natms, s.atoms.size());
  EXPECT_TRUE(std::filesystem::exists(dir / "CONTROL"));
  std::filesystem::remove_all(dir);
}
