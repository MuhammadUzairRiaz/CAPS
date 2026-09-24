#include <gtest/gtest.h>

#include <algorithm>
#include <tuple>
#include <array>
#include <cmath>
#include <filesystem>
#include <random>
#include <set>

#include "caps/ffdef.hpp"
#include "caps/field.hpp"
#include "caps/grow.hpp"
#include "caps/io.hpp"

using namespace caps;

namespace {

System small_cell(int chains, int dp, double density) {
  GrowOptions g;
  g.chains = chains;
  g.dp = dp;
  g.density = density;
  g.seed = 5;
  return grow(g);
}

std::vector<double> flat(const System& s) {
  std::vector<double> x;
  for (const auto& a : s.atoms) x.insert(x.end(), a.pos.begin(), a.pos.end());
  return x;
}

// The GAFF topology of a small polystyrene cell with every term replaced by a class II term with random,
// all non-zero coefficients (reference values near the actual geometry), and 9-6 LJ.
ForceField random_class2(const System& s, uint64_t seed) {
  ForceField ff = assign_gaff(s);
  std::mt19937_64 rng(seed);
  std::uniform_real_distribution<double> u(-1, 1);
  auto r = [&](double lo, double hi) { return lo + (hi - lo) * 0.5 * (u(rng) + 1); };
  constexpr double D = 3.14159265358979323846 / 180;
  for (const auto& b : ff.bonds) ff.bonds2.push_back({b.i, b.j, b.r0 + r(-0.05, 0.05), r(200, 400), r(-600, -200), r(200, 600)});
  for (const auto& a : ff.angles)
    ff.angles2.push_back({a.i, a.j, a.k, a.theta0 + r(-5, 5) * D, r(30, 60), r(-20, -5), r(-10, 10), r(-20, 20), r(1.0, 1.6), r(1.0, 1.6),
                          r(-20, 20), r(-20, 20), r(1.0, 1.6), r(1.0, 1.6)});
  std::set<std::tuple<uint32_t, uint32_t, uint32_t, uint32_t>> seen;
  for (const auto& t : ff.dihedrals) {
    if (!seen.insert({t.i, t.j, t.k, t.l}).second) continue;
    Class2Dihedral d{};
    d.i = t.i; d.j = t.j; d.k = t.k; d.l = t.l;
    d.k1 = r(-1, 1); d.phi1 = r(-180, 180) * D; d.k2 = r(-1, 1); d.phi2 = r(-180, 180) * D; d.k3 = r(-1, 1); d.phi3 = r(-180, 180) * D;
    for (int q = 0; q < 3; ++q) {
      d.mbt[q] = r(-10, 10); d.ebt_b[q] = r(-5, 5); d.ebt_c[q] = r(-5, 5); d.at_d[q] = r(-2, 2); d.at_e[q] = r(-2, 2);
    }
    d.mbt_r2 = r(1.3, 1.6); d.ebt_r1 = r(1.0, 1.6); d.ebt_r3 = r(1.0, 1.6);
    d.at_theta1 = r(100, 125) * D; d.at_theta2 = r(100, 125) * D;
    d.aat_m = r(-15, 15); d.aat_theta1 = r(100, 125) * D; d.aat_theta2 = r(100, 125) * D;
    d.bb13_n = r(-10, 10); d.bb13_r1 = r(1.0, 1.6); d.bb13_r3 = r(1.0, 1.6);
    ff.dihedrals2.push_back(d);
  }
  const auto nb = s.neighbours();
  for (uint32_t c = 0; c < s.atoms.size(); ++c)
    if (nb[c].size() >= 3)
      ff.impropers2.push_back({nb[c][0], c, nb[c][1], nb[c][2], r(5, 20), r(-5, 5) * D, r(-5, 5), r(-5, 5), r(-5, 5), r(100, 120) * D,
                               r(100, 120) * D, r(100, 120) * D});
  for (uint32_t c = 0; c < s.atoms.size(); ++c)
    if (nb[c].size() >= 3) ff.inversions.push_back({c, nb[c][0], nb[c][1], nb[c][2], r(5, 40), r(0, 20) * D});
  ff.bonds.clear();
  ff.angles.clear();
  ff.dihedrals.clear();
  ff.impropers.clear();
  ff.pair_form = "lj9-6";
  ff.mixing = "sixthpower";
  ff.lj14 = 1.0;
  ff.coul14 = 1.0;
  return ff;
}


// −dE/dε_ab of an affine strain of positions and cell (symmetric ε; Voigt order xx yy zz xy xz yz).
std::array<double, 6> strain_derivative(Evaluator& ev, const std::vector<double>& x, const Cell& cell, double h = 1e-6) {
  static const int ia[6] = {0, 1, 2, 0, 0, 1}, ib[6] = {0, 1, 2, 1, 2, 2};
  std::array<double, 6> out{};
  for (int v = 0; v < 6; ++v) {
    auto energy = [&](double e) {
      double m[3][3] = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
      if (ia[v] == ib[v]) m[ia[v]][ia[v]] += e;
      else { m[ia[v]][ib[v]] += e / 2; m[ib[v]][ia[v]] += e / 2; }
      auto map = [&](const Vec3& r) { return Vec3{m[0][0] * r[0] + m[0][1] * r[1] + m[0][2] * r[2], m[1][0] * r[0] + m[1][1] * r[1] + m[1][2] * r[2],
                                                  m[2][0] * r[0] + m[2][1] * r[1] + m[2][2] * r[2]}; };
      Cell c = cell;
      c.a = map(cell.a); c.b = map(cell.b); c.c = map(cell.c);
      std::vector<double> y = x, g;
      for (size_t i = 0; i < y.size(); i += 3) {
        const Vec3 r = map(Vec3{x[i] - cell.origin[0], x[i + 1] - cell.origin[1], x[i + 2] - cell.origin[2]});
        for (int k = 0; k < 3; ++k) y[i + k] = cell.origin[k] + r[k];
      }
      return ev.compute(y, c, g).total();
    };
    out[v] = -(energy(h) - energy(-h)) / (2 * h);
  }
  return out;
}
}  // namespace

TEST(Class2, ForcesMatchFiniteDifferences) {
  const System s = small_cell(2, 3, 0.25);
  const ForceField ff = random_class2(s, 11);
  ASSERT_FALSE(ff.impropers2.empty());
  std::vector<double> x = flat(s);
  std::mt19937_64 rng(3);
  std::normal_distribution<double> nd(0, 0.05);
  for (auto& v : x) v += nd(rng);
  for (double cap : {0.0, 20.0}) {
    EnergyOptions eo;
    eo.force_cap = cap;
    Evaluator ev(ff, eo);
    std::vector<double> f, g;
    const EnergyTerms e0 = ev.compute(x, s.cell, f);
    EXPECT_NE(e0.dihedral, 0.0);
    EXPECT_NE(e0.improper, 0.0);
    const double h = 1e-5;
    double worst = 0;
    for (size_t k = 0; k < x.size(); k += 3) {
      std::vector<double> xp = x, xm = x;
      xp[k] += h;
      xm[k] -= h;
      const double fd = -(ev.compute(xp, s.cell, g).total() - ev.compute(xm, s.cell, g).total()) / (2 * h);
      worst = std::max(worst, std::fabs(fd - f[k]) / std::max(1.0, std::fabs(f[k])));
    }
    EXPECT_LT(worst, 2e-4) << "force cap " << cap;
  }
}

TEST(Class2, VirialMatchesVolumeDerivative) {
  const System s = small_cell(2, 4, 0.35);
  const ForceField ff = random_class2(s, 12);
  EnergyOptions eo;
  eo.tail = false;
  Evaluator ev(ff, eo);
  std::vector<double> x = flat(s), f;
  const double w = ev.compute(x, s.cell, f).virial;
  auto energy_at = [&](double sc) {
    std::vector<double> y = x;
    Cell c = s.cell;
    for (size_t i = 0; i < y.size(); ++i) y[i] = c.origin[i % 3] + sc * (y[i] - c.origin[i % 3]);
    c.a = c.a * sc; c.b = c.b * sc; c.c = c.c * sc;
    std::vector<double> g;
    return ev.compute(y, c, g).total();
  };
  const double h = 1e-6, dEds = (energy_at(1 + h) - energy_at(1 - h)) / (2 * h);
  EXPECT_NEAR(w, -dEds, 1e-4 * std::max(1.0, std::fabs(w)));
  const EnergyTerms et = ev.compute(x, s.cell, f);
  const auto d = strain_derivative(ev, x, s.cell);
  for (int v = 0; v < 6; ++v) EXPECT_NEAR(et.w[v], d[v], 2e-4 * std::max(1.0, std::fabs(w))) << "class II virial component " << v;
}

TEST(Class2, SixthPowerMixing) {
  ForceField ff;
  ff.mixing = "sixthpower";
  ff.lj = {{0.1, 3.0}, {0.4, 4.0}};
  const PairType p = mixed_pair(ff, 0, 1);
  const double s6a = std::pow(3.0, 6), s6b = std::pow(4.0, 6);
  EXPECT_NEAR(p.sigma, std::pow(0.5 * (s6a + s6b), 1.0 / 6), 1e-12);
  EXPECT_NEAR(p.eps, 2 * std::sqrt(0.04) * 27 * 64 / (s6a + s6b), 1e-12);
}

// moltemplate's COMPASS file converts with its cross terms, bond increments and centre-second impropers.
TEST(FFDef, ImportsCompassFromMoltemplate) {
  const std::string lt = std::string(std::getenv("HOME") ? std::getenv("HOME") : "") + "/moltemplate/moltemplate/force_fields/compass_published.lt";
  if (!std::filesystem::exists(lt)) GTEST_SKIP() << "moltemplate not installed";
  const FFDef ff = import_moltemplate(lt);
  EXPECT_EQ(ff.bond_style, "class2");
  EXPECT_EQ(ff.improper_order, "center2_sorted");
  EXPECT_EQ(ff.mixing, "sixthpower");
  EXPECT_DOUBLE_EQ(ff.special_lj[2], 1.0);
  EXPECT_FALSE(ff.bond_increments.empty());
  bool angle_cross = false, dihedral_cross = false, improper_cross = false;
  for (const auto& r : ff.angles) angle_cross = angle_cross || (r.cross.count("bb") && r.cross.count("ba") && r.params.size() == 4);
  for (const auto& r : ff.dihedrals)
    dihedral_cross = dihedral_cross || (r.cross.size() == 5 && r.cross.at("ebt").size() == 8 && r.params.size() == 6);
  for (const auto& r : ff.impropers) improper_cross = improper_cross || (r.cross.count("aa") && r.cross.at("aa").size() == 6);
  EXPECT_TRUE(angle_cross);
  EXPECT_TRUE(dihedral_cross);
  EXPECT_TRUE(improper_cross);
  // JSON round trip keeps the cross terms
  const std::string path = (std::filesystem::temp_directory_path() / "caps_compass_test.json").string();
  save_forcefield(ff, path);
  const FFDef back = load_forcefield(path);
  ASSERT_EQ(back.dihedrals.size(), ff.dihedrals.size());
  EXPECT_EQ(back.dihedrals[3].cross, ff.dihedrals[3].cross);
  EXPECT_EQ(back.bond_increments.size(), ff.bond_increments.size());
  std::filesystem::remove(path);
}

TEST(FFDef, GlobEscapeKeepsLiteralStars) {
  EXPECT_TRUE(glob_match("h*", "hc"));
  EXPECT_FALSE(glob_match(glob_escape("h*"), "hc"));
  EXPECT_TRUE(glob_match(glob_escape("h*"), "h*"));
  EXPECT_TRUE(glob_match("*~b" + glob_escape("c=") + "~*", "c=1~bc=~a"));
}

TEST(IO, Mol2RoundTripKeepsTypesChargesAndBondOrders) {
  System s;
  s.title = "formaldehyde";
  Atom c, o, h1, h2;
  c.element = 6; c.name = "c"; c.charge = 0.5; c.pos = {0, 0, 0};
  o.element = 8; o.name = "o"; o.charge = -0.5; o.pos = {1.21, 0, 0};
  h1.element = 1; h1.name = "h4"; h1.pos = {-0.55, 0.94, 0};
  h2.element = 1; h2.name = "h4"; h2.pos = {-0.55, -0.94, 0};
  s.atoms = {c, o, h1, h2};
  s.bonds = {{0, 1, 2}, {0, 2, 1}, {0, 3, 1}};
  s.has_charges = true;
  const std::string path = (std::filesystem::temp_directory_path() / "caps_test.mol2").string();
  write_mol2(s, path);
  const System r = read_mol2(path);
  ASSERT_EQ(r.atoms.size(), 4u);
  EXPECT_EQ(r.atoms[1].name, "o");
  EXPECT_EQ(r.atoms[1].element, 8);
  EXPECT_NEAR(r.atoms[0].charge, 0.5, 1e-9);
  ASSERT_EQ(r.bonds.size(), 3u);
  EXPECT_EQ(r.bonds[0].order, 2);
  EXPECT_EQ(detect_format(path), "mol2");
  std::filesystem::remove(path);
}

// DL_FIELD's PCFF converts with equivalences, bond increments and the cff91_auto tables.
TEST(FFDef, ImportsDlfieldPcff) {
  const std::string par = "~/project/dl_f_4.13/lib/PCFF.par";
  if (!std::filesystem::exists(par)) GTEST_SKIP() << "DL_FIELD not installed";
  const FFDef ff = import_dlfield(par, "~/project/dl_f_4.13/lib/PCFF.sf", "~/project/dl_f_4.13/lib/PCFF.bci");
  EXPECT_EQ(ff.mixing, "sixthpower");
  EXPECT_FALSE(ff.auto_dihedrals.empty());
  EXPECT_FALSE(ff.bond_increments.empty());
  const FFType* hc = ff.type("hc");
  ASSERT_NE(hc, nullptr);
  EXPECT_EQ(hc->equiv.at("vdw"), "h");
  EXPECT_EQ(hc->element, 1);
  EXPECT_NEAR(hc->mass, 1.00797, 1e-6);   // not the 12.0115 typo of HC_benzyl
}

// Morse / GROMOS bonds, cosine-squared and linear-cosine angles, Urey–Bradley, planar inversions, Buckingham pairs
// and separate 1-4 Lennard-Jones: forces are the gradient of the energy.
TEST(FieldForms, ForcesMatchFiniteDifferences) {
  const System s = small_cell(2, 3, 0.25);
  ForceField ff = assign_gaff(s);
  std::mt19937_64 rng(21);
  std::uniform_real_distribution<double> u(0, 1);
  auto r = [&](double lo, double hi) { return lo + (hi - lo) * u(rng); };
  constexpr double D = 3.14159265358979323846 / 180;
  for (size_t k = 0; k < ff.bonds.size(); ++k) {
    const auto& b = ff.bonds[k];
    ff.bonds_x.push_back(k % 2 ? BondX{b.i, b.j, 1, r(50, 100), r(1.5, 2.5), b.r0} : BondX{b.i, b.j, 2, r(20, 60), b.r0, 0});
  }
  for (size_t k = 0; k < ff.angles.size(); ++k) {
    const auto& a = ff.angles[k];
    ff.angles_x.push_back({a.i, a.j, a.k, k % 3 == 0 ? 2 : 1, r(20, 60), a.theta0 + r(-5, 5) * D});
    ff.urey_bradley.push_back({a.i, a.k, r(5, 30), r(2.0, 2.6)});
  }
  const auto nb = s.neighbours();
  for (uint32_t c = 0; c < s.atoms.size(); ++c)
    if (nb[c].size() == 3) ff.inversions.push_back({c, nb[c][0], nb[c][1], nb[c][2], r(10, 50), 0.0, 1});
  const int nt = int(ff.lj.size());
  ff.pair_func[{0, nt - 1}] = {1, 5000.0, 0.3, 20.0};
  for (const auto& p : ff.lj) ff.lj14_types.push_back({p.eps * 0.5, p.sigma * 0.9});
  ff.lj14 = 1.0;
  std::vector<double> x = flat(s);
  std::normal_distribution<double> nd(0, 0.05);
  for (auto& v : x) v += nd(rng);
  // shifted pair energies (tail off): with tail corrections LJ is truncated at the cut-off, a step a finite
  // difference would see when a pair sits within h of it
  EnergyOptions shifted;
  shifted.tail = false;
  Evaluator ev(ff, shifted);
  std::vector<double> f, g;
  ev.compute(x, s.cell, f);
  const double h = 1e-5;
  double worst = 0;
  for (size_t k = 0; k < x.size(); k += 3) {
    std::vector<double> xp = x, xm = x;
    xp[k] += h;
    xm[k] -= h;
    const double fd = -(ev.compute(xp, s.cell, g).total() - ev.compute(xm, s.cell, g).total()) / (2 * h);
    worst = std::max(worst, std::fabs(fd - f[k]) / std::max(1.0, std::fabs(f[k])));
  }
  EXPECT_LT(worst, 2e-4);
  // virial = −dE/ds for an affine scaling
  EnergyOptions nt0;
  nt0.tail = false;
  Evaluator e2(ff, nt0);
  const double w = e2.compute(x, s.cell, f).virial;
  auto energy_at = [&](double sc) {
    std::vector<double> y = x;
    Cell c = s.cell;
    for (size_t i = 0; i < y.size(); ++i) y[i] = c.origin[i % 3] + sc * (y[i] - c.origin[i % 3]);
    c.a = c.a * sc; c.b = c.b * sc; c.c = c.c * sc;
    return e2.compute(y, c, g).total();
  };
  EXPECT_NEAR(w, -(energy_at(1 + 1e-6) - energy_at(1 - 1e-6)) / 2e-6, 1e-4 * std::max(1.0, std::fabs(w)));
  const EnergyTerms et = e2.compute(x, s.cell, f);
  const auto d = strain_derivative(e2, x, s.cell);
  for (int v = 0; v < 6; ++v) EXPECT_NEAR(et.w[v], d[v], 2e-4 * std::max(1.0, std::fabs(w))) << "virial component " << v;
}

// CHARMM libraries keep their separate 1-4 van der Waals and Urey–Bradley terms through save / load.
TEST(FFDef, DlfieldCharmmKeepsOneFourAndUreyBradley) {
  const std::string par = "~/project/dl_f_4.13/lib/CHARMM36_cgenff.par";
  if (!std::filesystem::exists(par)) GTEST_SKIP() << "DL_FIELD not installed";
  const FFDef ff = import_dlfield(par, "~/project/dl_f_4.13/lib/CHARMM36_cgenff.sf");
  const std::string path = (std::filesystem::temp_directory_path() / "caps_cgenff_test.json").string();
  save_forcefield(ff, path);
  const FFDef back = load_forcefield(path);
  std::filesystem::remove(path);
  bool lj14 = false, ub = false;
  for (const auto& r : back.pairs) lj14 = lj14 || (r.match.size() == 1 && r.match[0] == "OG2D1" && r.cross.count("lj14"));
  for (const auto& r : back.angles) ub = ub || (r.style == "charmm" && r.params.size() == 4 && r.params[2] > 0);
  EXPECT_TRUE(lj14);
  EXPECT_TRUE(ub);
  EXPECT_DOUBLE_EQ(back.special_lj[2], 1.0);
}
