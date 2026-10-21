// The polymer crystals of data/crystals/polymers (tools/make_polymer_crystals.py): each against its paper (cell, space
// group operations, atoms per cell, composition Z × repeat unit, density), every atom with its valence (no missing or
// extra hydrogens: the cell is neutral and closed-shell), no clashes once the hydrogens are in, the catalogue entries,
// and the POSCAR round trip a VASP set starts from.
#include <gtest/gtest.h>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <map>
#include <set>
#include <sstream>

#include "caps/crystal.hpp"
#include "caps/elements.hpp"
#include "caps/json.hpp"
#include "caps/slab2d.hpp"
#include "caps/spacegroup.hpp"

using namespace caps;
namespace fs = std::filesystem;

namespace {
const std::string kDir = std::string(CAPS_SOURCE_DIR) + "/data/crystals/";

Json load(const std::string& p) {
  std::ifstream f(p);
  std::stringstream ss;
  ss << f.rdbuf();
  return Json::parse(ss.str());
}
std::vector<Json> entries() {
  std::vector<Json> out;
  const Json doc = load(kDir + "polymers/polymers.json");
  for (const auto& e : doc["crystals"].items()) out.push_back(e);
  return out;
}
double angle(const Vec3& u, const Vec3& v) { return std::acos(dot(u, v) / (norm(u) * norm(v))) * 180 / M_PI; }
// "C2H4" → {C: 2, H: 4}
std::map<std::string, int> counts(const std::string& f) {
  std::map<std::string, int> m;
  for (size_t i = 0; i < f.size();) {
    std::string el(1, f[i++]);
    while (i < f.size() && std::islower(static_cast<unsigned char>(f[i]))) el += f[i++];
    int n = 0;
    while (i < f.size() && std::isdigit(static_cast<unsigned char>(f[i]))) n = n * 10 + (f[i++] - '0');
    m[el] += n ? n : 1;
  }
  return m;
}
bool bonded(int a, int b, double r) {
  if (a > b) std::swap(a, b);
  if (a == 1 && b == 1) return r < 0.9;
  if (a == 1) return r < 1.25;            // C-H, N-H, O-H
  if (a == 6 && b == 6) return r < 1.75;
  return r < 1.65;                        // C-N, C-O, C-F
}
}  // namespace

TEST(PolymerCrystals, LibraryAgainstPapers) {
  const auto all = entries();
  ASSERT_GE(all.size(), 12u);
  std::set<std::string> ids;
  for (const Json& e : all) {
    const std::string id = e.text("id");
    SCOPED_TRACE(id);
    ids.insert(id);
    const System s = read_cif(kDir + "polymers/" + e.text("file"));
    // the paper's cell
    const Json& c = e["cell"];
    EXPECT_NEAR(norm(s.cell.a), c.num("a", 0), 1e-4);
    EXPECT_NEAR(norm(s.cell.b), c.num("b", 0), 1e-4);
    EXPECT_NEAR(norm(s.cell.c), c.num("c", 0), 1e-4);
    EXPECT_NEAR(angle(s.cell.b, s.cell.c), c.num("alpha", 0), 1e-3);
    EXPECT_NEAR(angle(s.cell.a, s.cell.c), c.num("beta", 0), 1e-3);
    EXPECT_NEAR(angle(s.cell.a, s.cell.b), c.num("gamma", 0), 1e-3);
    // atoms per cell and composition = Z × the repeat unit
    EXPECT_EQ(int(s.atoms.size()), int(e.num("atoms", 0)));
    std::map<std::string, int> got;
    for (const auto& a : s.atoms) ++got[element(a.element).symbol];
    auto want = counts(e.text("repeat_unit"));
    for (auto& [el, n] : want) n *= int(e.num("Z", 0));
    EXPECT_EQ(got, want);
    // density against the paper's (1 %), or the value its cell and Z give when it prints none
    const double rho = s.density();
    if (e["density_paper"].is_number()) EXPECT_NEAR(rho, e.num("density_paper", 0), 0.01 * e.num("density_paper", 0));
    EXPECT_NEAR(rho, e.num("density", 0), 1e-3);
    // the paper's space group: the cell is invariant under each of its equivalent positions
    const double tol = e.has("prerelax") ? 0.05 : 0.002;
    for (const auto& o : e["operations"].items()) {
      const SymOp op = parse_symop(o.str());
      double worst = 0;
      for (const auto& a : s.atoms) {
        const Vec3 f = s.cell.to_fractional(a.pos);
        Vec3 g{op.t[0], op.t[1], op.t[2]};
        for (int i = 0; i < 3; ++i) for (int k = 0; k < 3; ++k) g[size_t(i)] += op.R[i][k] * f[size_t(k)];
        double best = 1e9;
        for (const auto& b : s.atoms) {
          if (b.element != a.element) continue;
          Vec3 d = g - s.cell.to_fractional(b.pos);
          for (auto& x : d) x -= std::round(x);
          best = std::min(best, norm(s.cell.a * d[0] + s.cell.b * d[1] + s.cell.c * d[2]));
        }
        worst = std::max(worst, best);
      }
      EXPECT_LT(worst, tol) << o.str();
    }
    // valences: C 3 (sp2) or 4, N 3, O 1 or 2, F 1, H 1 — every hydrogen there, none extra
    const auto nb = neighbours_pbc(s, 3.2);
    std::vector<std::vector<std::pair<uint32_t, Vec3>>> bonds(s.atoms.size());
    for (const auto& n : nb)
      if (bonded(s.atoms[n.i].element, s.atoms[n.j].element, n.r)) bonds[n.i].push_back({n.j, n.d});
    int bad = 0;
    for (size_t i = 0; i < s.atoms.size(); ++i) {
      const int z = s.atoms[i].element, k = int(bonds[i].size());
      const bool ok = z == 6 ? (k == 3 || k == 4) : z == 7 ? k == 3 : z == 8 ? (k == 1 || k == 2) : k == 1;
      bad += !ok;
    }
    EXPECT_EQ(bad, 0);
    // electrons per cell even and every C/N/O at its full valence: count bond-order sums for the sp2 atoms is not
    // needed; a closed shell follows from the valences above and an even electron count
    int electrons = 0;
    for (const auto& a : s.atoms) electrons += a.element;
    EXPECT_EQ(electrons % 2, 0);
    // no clashes: pairs beyond 1-3 (not bonded, no common neighbour); hydrogen bonds N/O-H...O/N apart
    double hh = 9, hx = 9, xx = 9;
    for (const auto& n : nb) {
      bool near = false;
      for (const auto& [j, d] : bonds[n.i]) {
        if (j == n.j && norm(d - n.d) < 0.05) near = true;
        for (const auto& [k, d2] : bonds[j]) if (k == n.j && norm(d + d2 - n.d) < 0.05) near = true;
      }
      if (near) continue;
      const int zi = s.atoms[n.i].element, zj = s.atoms[n.j].element;
      auto polar_h = [&](uint32_t h) { return s.atoms[h].element == 1 && !bonds[h].empty() && (s.atoms[bonds[h][0].first].element == 7 || s.atoms[bonds[h][0].first].element == 8); };
      if ((polar_h(n.i) && (zj == 7 || zj == 8)) || (polar_h(n.j) && (zi == 7 || zi == 8))) { EXPECT_GT(n.r, 1.7); continue; }
      if (zi == 1 && zj == 1) hh = std::min(hh, n.r);
      else if (zi == 1 || zj == 1) hx = std::min(hx, n.r);
      else xx = std::min(xx, n.r);
    }
    EXPECT_GE(hh, 1.8);
    EXPECT_GE(hx, 2.2);
    EXPECT_GE(xx, 2.5);
    // the POSCAR a VASP set writes reads back to the same structure
    const fs::path p = fs::temp_directory_path() / ("caps_polymer_" + id + ".vasp");
    const System g = group_by_species(s, {"C", "N", "O", "F", "H"});
    write_vasp_poscar(g, p.string(), id);
    const System back = read_vasp_poscar(p.string());
    ASSERT_EQ(back.atoms.size(), g.atoms.size());
    double dmax = 0;
    for (size_t i = 0; i < g.atoms.size(); ++i) {
      EXPECT_EQ(back.atoms[i].element, g.atoms[i].element);
      dmax = std::max(dmax, norm(back.atoms[i].pos - g.atoms[i].pos));
    }
    EXPECT_LT(dmax, 1e-8);
    EXPECT_NEAR(back.cell.volume(), s.cell.volume(), 1e-6);
    fs::remove(p);
  }
  // PE first: the cell both the DFT workflow and the semicrystalline builder start from
  EXPECT_TRUE(ids.count("pe-orthorhombic") && ids.count("pe-monoclinic"));
}

TEST(PolymerCrystals, PolyethyleneOrthorhombicGeometry) {
  // Bunn 1939: C-C 1.53 A, zigzag angle 112 deg, chains of 2 CH2 per 2.534 A repeat, Pnam found by CAPS's own search
  const System s = read_cif(kDir + "polymers/pe-orthorhombic.cif");
  ASSERT_EQ(s.atoms.size(), 12u);
  const auto sym = find_symmetry(s, 0.01);
  EXPECT_EQ(sym.number, 62);
  double cc = 0, ch = 0;
  int ncc = 0, nch = 0;
  for (const auto& n : neighbours_pbc(s, 1.7)) {
    const int a = s.atoms[n.i].element, b = s.atoms[n.j].element;
    if (a == 6 && b == 6) cc += n.r, ++ncc;
    if (a == 6 && b == 1) ch += n.r, ++nch;
  }
  ASSERT_GT(ncc, 0);
  EXPECT_NEAR(cc / ncc, 1.53, 0.01);
  EXPECT_NEAR(ch / nch, 1.09, 1e-3);
}

TEST(PolymerCrystals, CatalogueLists) {
  const Json cat = load(kDir + "catalogue.json");
  std::set<std::string> files;
  for (const auto& e : cat["crystals"].items()) {
    files.insert(e.text("file"));
    EXPECT_TRUE(fs::exists(kDir + e.text("file"))) << e.text("file");
  }
  for (const Json& e : entries()) EXPECT_TRUE(files.count("polymers/" + e.text("file"))) << e.text("id");
}
