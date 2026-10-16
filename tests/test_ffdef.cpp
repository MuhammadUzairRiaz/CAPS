#include <gtest/gtest.h>

#include <algorithm>
#include <cstdlib>
#include <tuple>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <map>
#include <random>
#include <set>

#include "caps/analysis.hpp"
#include "caps/crystal.hpp"
#include "caps/ffdef.hpp"
#include "caps/field.hpp"
#include "caps/grow.hpp"
#include "caps/io.hpp"
#include "caps/martini3_molecules.hpp"
#include "caps/molecule.hpp"
#include "caps/resolution.hpp"
#include "caps/spacegroup.hpp"
#include "caps/typing.hpp"
#include "caps/uff.hpp"

using namespace caps;

namespace {

// The first of these library ids whose file is in data/forcefields, else "" (a model being rebuilt from its own papers)
std::string library_file(std::initializer_list<const char*> ids) {
  for (const char* id : ids) {
    const std::string p = std::string(CAPS_SOURCE_DIR) + "/data/forcefields/" + id + ".json";
    if (std::filesystem::exists(p)) return p;
  }
  return "";
}

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

// Morse / GROMOS bonds, cosine-squared and linear-cosine angles, Urey–Bradley, planar inversions, Buckingham and Rydberg
// pairs and separate 1-4 Lennard-Jones: forces are the gradient of the energy.
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
    if (k % 4 == 3)   // the screened harmonic bend (Tilocca's SM1): its radial part adds forces along the arms
      ff.angles_x.push_back({a.i, a.j, a.k, 6, r(200, 400), a.theta0 + r(-15, 15) * D, r(0.8, 1.2), r(0.8, 1.2)});
    else
      ff.angles_x.push_back({a.i, a.j, a.k, k % 3 == 0 ? 2 : 1, r(20, 60), a.theta0 + r(-5, 5) * D});
    ff.urey_bradley.push_back({a.i, a.k, r(5, 30), r(2.0, 2.6)});
  }
  const auto nb = s.neighbours();
  for (uint32_t c = 0; c < s.atoms.size(); ++c)
    if (nb[c].size() == 3) ff.inversions.push_back({c, nb[c][0], nb[c][1], nb[c][2], r(10, 50), 0.0, 1});
  const int nt = int(ff.lj.size());
  ff.pair_func[{0, nt - 1}] = {1, 5000.0, 0.3, 20.0};
  ASSERT_GT(nt, 2);
  ff.pair_func[{1, nt - 1}] = {kPairRydberg, 3000.0, 0.35, -150.0};   // (A + C r) e^(−r/ρ), C < 0 as PbTiO3's Ti-O
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

// mW water: the Stillinger–Weber two- and three-body forces and the virial are the energy's derivatives (sites on a
// jittered lattice at liquid density, pairs and triplets across the periodic faces)
TEST(FieldForms, StillingerWeberForcesMatchFiniteDifferences) {
  // mW water from Molinero & Moore 2009 (bench/ff/make_water_models.py)
  const std::string mw = library_file({"mw-molinero2009"});
  if (mw.empty()) GTEST_SKIP() << "no mW water library yet (mw-molinero2009)";
  FFDef def = load_forcefield(mw);
  System s;
  const int n = 4;
  const double a = 3.1;   // 64 sites in 12.4 Å: 1.00 g/cm³
  s.cell.a = {n * a, 0, 0};
  s.cell.b = {0, n * a, 0};
  s.cell.c = {0, 0, n * a};
  std::mt19937_64 rng(7);
  std::normal_distribution<double> nd(0, 0.35);
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      for (int k = 0; k < n; ++k) {
        Atom at;
        at.element = 8;
        at.name = "MW";
        at.pos = {i * a + nd(rng), j * a + nd(rng), k * a + nd(rng)};
        s.atoms.push_back(at);
      }
  ParamReport rep;
  const ForceField ff = parameterize(s, def, std::vector<std::string>(s.atoms.size(), "MW"), "types", &rep, false);
  ASSERT_TRUE(ff.sw.on);
  EnergyOptions o;
  o.coulomb = false;
  o.tail = false;
  Evaluator ev(ff, o);
  std::vector<double> x = flat(s), f, g;
  const EnergyTerms et = ev.compute(x, s.cell, f);
  EXPECT_LT(et.vdw, -100.0);   // bound: about −6 kcal/mol per site
  const double h = 1e-5;
  double worst = 0;
  for (size_t k = 0; k < x.size(); ++k) {
    std::vector<double> xp = x, xm = x;
    xp[k] += h;
    xm[k] -= h;
    const double fd = -(ev.compute(xp, s.cell, g).total() - ev.compute(xm, s.cell, g).total()) / (2 * h);
    worst = std::max(worst, std::fabs(fd - f[k]) / std::max(1.0, std::fabs(f[k])));
  }
  EXPECT_LT(worst, 1e-5);
  const auto d = strain_derivative(ev, x, s.cell);
  for (int v = 0; v < 6; ++v) EXPECT_NEAR(et.w[v], d[v], 1e-5 * std::max(1.0, std::fabs(et.virial))) << "virial component " << v;
}

// Coarse-grained forms: MARTINI (lj/gromacs, coul/gromacs with εr 15, cosine/squared angles, 1-3 pairs kept) and SDK
// (lj/sdk 9-6 / 12-4, angle sdk with its 1-3 repulsion): forces and virial are the energy's derivatives, for bead
// molecules and water beads in a periodic cell
namespace {
System cg_box(const FFDef& def, const std::vector<std::pair<std::string, int>>& parts, double edge, uint64_t seed, double clear = 4.0) {
  System box;
  box.cell.a = {edge, 0, 0};
  box.cell.b = {0, edge, 0};
  box.cell.c = {0, 0, edge};
  std::mt19937_64 rng(seed);
  std::uniform_real_distribution<double> u(0, edge);
  for (const auto& [text, count] : parts)
    for (int c = 0; c < count; ++c) {
      const System m = build_bead_molecule(text, def, seed + c);
      Vec3 t{0, 0, 0};
      for (int tr = 0; tr < 5000; ++tr) {   // a place `clear` Å from every bead already in the box
        t = {u(rng), u(rng), u(rng)};
        bool ok = true;
        for (const auto& a : m.atoms)
          for (const auto& b : box.atoms)
            ok = ok && norm(box.cell.minimum_image(a.pos + t - b.pos)) > clear;
        if (ok) break;
      }
      const uint32_t off = uint32_t(box.atoms.size());
      const int64_t mol = box.atoms.empty() ? 1 : box.atoms.back().mol + 1;
      for (auto a : m.atoms) {
        a.pos = a.pos + t;
        a.mol = mol;   // molecules by id: a virtual site is bonded to nothing
        a.id = int64_t(box.atoms.size()) + 1;
        box.atoms.push_back(a);
      }
      for (auto b : m.bonds) box.bonds.push_back({b.i + off, b.j + off, 1});
    }
  box.bonds_from_file = true;
  box.has_mol = true;
  box.has_charges = true;
  return box;
}

void check_cg_forces(const FFDef& def, System s, const std::string& charges) {
  std::string ch = charges;
  prepare_for_forcefield(s, def, ch);
  const TypingResult tr = assign_types(s, def);
  ASSERT_EQ(tr.untyped, 0);
  ParamReport rep;
  const ForceField ff = parameterize(s, def, tr.types, ch, &rep, false);
  EnergyOptions o;
  o.tail = false;
  Evaluator ev(ff, o);
  std::vector<double> x, f, g;
  for (const auto& a : s.atoms) x.insert(x.end(), {a.pos[0], a.pos[1], a.pos[2]});
  const EnergyTerms et = ev.compute(x, s.cell, f);
  EXPECT_TRUE(std::isfinite(et.total()));
  double worst = 0;
  for (size_t k = 0; k < x.size(); ++k) {
    std::vector<double> xp = x, xm = x;
    xp[k] += 1e-5;
    xm[k] -= 1e-5;
    const double fd = -(ev.compute(xp, s.cell, g).total() - ev.compute(xm, s.cell, g).total()) / 2e-5;
    const double err = std::fabs(fd - f[k]) / std::max(1.0, std::fabs(f[k]));
    if (err > 1e-4) std::printf("  %s atom %zu (%s) comp %zu: fd %.8g f %.8g\n", def.name.c_str(), k / 3, s.atoms[k / 3].name.c_str(), k % 3, fd, f[k]);
    worst = std::max(worst, err);
  }
  EXPECT_LT(worst, 1e-4) << def.name;
  const auto d = strain_derivative(ev, x, s.cell);
  for (int v = 0; v < 6; ++v) EXPECT_NEAR(et.w[v], d[v], 1e-4 * std::max(1.0, std::fabs(et.virial))) << def.name << " virial " << v;
}
}  // namespace

TEST(FieldForms, CoarseGrainedForcesMatchFiniteDifferences) {
  const std::string dir = std::string(CAPS_SOURCE_DIR) + "/data/forcefields/";
  // Cooke-Kremer-Deserno from the paper (cooke-deserno2005). TODO: MARTINI 2, its amino acids and SDK once rebuilt from
  // their papers (not in hand yet; expected ids martini2-marrink2007, martini2-aminoacids-monticelli2008, sdk-shinoda2007)
  if (const std::string f = library_file({"martini2-marrink2007"}); !f.empty()) {
    const FFDef martini = load_forcefield(f);
    check_cg_forces(martini, cg_box(martini, {{"DPPC", 3}, {"NA+", 2}, {"CL-", 2}, {"W", 30}}, 36.0, 3), "keep");
  }
  if (const std::string f = library_file({"sdk-shinoda2007"}); !f.empty()) {
    const FFDef sdk = load_forcefield(f);
    check_cg_forces(sdk, cg_box(sdk, {{"[NC][PH][GL]([EST1][CM][CM][CT2])[EST2][CM][CM][CT2]", 3}, {"[W]", 30}}, 32.0, 5), "types");
  }
  if (const std::string f = library_file({"martini2-aminoacids-monticelli2008"}); !f.empty()) {
    const FFDef aa = load_forcefield(f);   // harmonic impropers, GROMACS type-2 order
    check_cg_forces(aa, cg_box(aa, {{"TRP", 2}, {"HIS", 2}, {"PHE", 2}, {"TYR", 2}, {"LYS", 2}, {"ASP", 2}, {"[P4]", 30}}, 36.0, 9), "keep");
  }
  if (const std::string f = library_file({"cooke-deserno2005"}); !f.empty()) {
    const FFDef cd = load_forcefield(f);
    check_cg_forces(cd, cg_box(cd, {{"lipid", 40}}, 12.0, 7, 1.2), "types");
  }
  // Martini 3: reaction field (εr 15, εrf ∞), LJ shifted at 11 Å, every pair from the table; its molecules with their
  // own topology: virtual sites (nucleobases, the weighted centre), exclusions, restricted bending (TXE), impropers
  const FFDef m3 = load_forcefield(dir + "martini3.json");
  check_cg_forces(m3, cg_box(m3, {{"THYM", 2}, {"GUAN", 2}, {"TXE", 2}, {"POPS", 2}, {"DIM", 2}, {"NA", 4}, {"CL", 4}, {"W", 30}}, 30.0, 11), "auto");
}

// Martini 3's pair table and molecule templates: pairs as martini_v3.0.0.itp gives them (no mixing), a molecule's own
// topology when built, and recognised again in the structure written to a data file and read back
TEST(CoarseGrained, Martini3TemplatesAndPairTable) {
  const FFDef m3 = load_forcefield(std::string(CAPS_SOURCE_DIR) + "/data/forcefields/martini3.json");
  ASSERT_TRUE(m3.molecule_templates);
  EXPECT_GT(bead_template_list(m3).size(), 200u);
  const System w = build_bead_molecule("THYM", m3, 1);
  ASSERT_TRUE(w.topology);
  EXPECT_EQ(w.atoms.size(), 5u);
  EXPECT_EQ(w.topology->vsites.size(), 1u);
  System box = cg_box(m3, {{"W", 20}, {"NA", 2}, {"CL", 2}, {"THYM", 2}, {"POPC", 2}}, 30.0, 2);
  std::map<std::string, int> tix;   // one numeric type per bead name, labelled, as CAPS Pack writes a box
  for (auto& a : box.atoms) {
    auto [it, fresh] = tix.emplace(a.name, int(tix.size()) + 1);
    a.type = it->second;
    if (fresh)
      for (const auto& t : m3.types)
        if (t.name == a.name) box.types.push_back({it->second, t.mass, a.name});
  }
  const std::string path = (std::filesystem::temp_directory_path() / "caps_m3_box.data").string();
  write_lammps_data(box, path);
  System s = open_file(path).frame(0);
  std::string ch = "auto";
  const std::string note = prepare_for_forcefield(s, m3, ch);
  ASSERT_TRUE(s.topology) << note;
  EXPECT_EQ(ch, "keep");
  EXPECT_NE(note.find("28 molecules recognised"), std::string::npos) << note;
  const TypingResult tr = assign_types(s, m3);
  EXPECT_EQ(tr.untyped, 0);
  ParamReport rep;
  const ForceField ff = parameterize(s, m3, tr.types, ch, &rep, false);
  EXPECT_TRUE(ff.coul_rf);
  EXPECT_TRUE(ff.lj_shift);
  EXPECT_EQ(ff.vsites.size(), 2u);
  double q = 0;
  for (double c : ff.charge) q += c;
  EXPECT_NEAR(q, 0.0, 1e-12);
  auto ti = [&](const std::string& n) { return int(std::find(ff.type_names.begin(), ff.type_names.end(), n) - ff.type_names.begin()); };
  // martini_v3.0.0.itp: W–W σ 0.47 nm, ε 4.65 kJ/mol; TQ5 (Na+) – W σ 0.385 nm, ε 11.46 kJ/mol (not a mixing rule's)
  const int iw = ti("W"), iq = ti("TQ5");
  ASSERT_LT(iw, int(ff.type_names.size()));
  ASSERT_LT(iq, int(ff.type_names.size()));
  EXPECT_NEAR(ff.lj[iw].sigma, 4.7, 1e-12);
  EXPECT_NEAR(ff.lj[iw].eps, 4.65 / 4.184, 1e-12);
  const PairType& p = ff.pair_override.at({std::min(iw, iq), std::max(iw, iq)});
  EXPECT_NEAR(p.sigma, 3.85, 1e-12);
  EXPECT_NEAR(p.eps, 11.46 / 4.184, 1e-12);
  std::filesystem::remove(path);
}

TEST(CoarseGrained, BeadSmilesParse) {
  const BeadMolecule m = parse_bead_smiles("[Q0+1][Qa-1][Na]([Na][C1])[C1].[SC4]1[SC4][SC4]1[Qa1-1]");
  ASSERT_EQ(m.type.size(), 10u);
  EXPECT_EQ(m.type[0], "Q0");
  EXPECT_EQ(m.charge[0], 1.0);
  EXPECT_EQ(m.type[1], "Qa");
  EXPECT_EQ(m.charge[1], -1.0);
  EXPECT_EQ(m.type[9], "Qa1");
  EXPECT_EQ(m.charge[9], -1.0);
  EXPECT_EQ(m.bonds.size(), 9u);   // 5 in the lipid fragment, 3 in the ring, 1 to Qa1 ('.' separates)
  EXPECT_THROW(parse_bead_smiles("[C1]1[C1]"), std::invalid_argument);
  const System s = build_beads("[C1][C1][C1][C1]");
  ASSERT_EQ(s.atoms.size(), 4u);
  for (const auto& b : s.bonds) EXPECT_NEAR(norm(s.atoms[b.j].pos - s.atoms[b.i].pos), 4.7, 0.05);
}

// SDK: all-atom alkanes tile into CT / CT2 ends and CM between for every chain length; DMPC maps onto Shinoda's
// thirteen beads with the esters told apart by their GL bond lengths
TEST(CoarseGrained, SdkMapsAllAtomStructures) {
  // TODO: SDK rebuilt from Shinoda, DeVane & Klein 2007 (expected id sdk-shinoda2007)
  const std::string f = library_file({"sdk-shinoda2007"});
  if (f.empty()) GTEST_SKIP() << "no SDK library yet (sdk-shinoda2007)";
  const FFDef sdk = load_forcefield(f);
  auto map = [&](const std::string& smiles) {
    System s = build_molecule(smiles, {}).system;
    std::string ch = "auto";
    prepare_for_forcefield(s, sdk, ch);
    std::map<std::string, int> n;
    for (const auto& a : s.atoms) n[a.name]++;
    return n;
  };
  EXPECT_EQ(map(std::string(12, 'C')), (std::map<std::string, int>{{"CM", 2}, {"CT", 2}}));
  EXPECT_EQ(map(std::string(13, 'C')), (std::map<std::string, int>{{"CM", 3}, {"CT2", 2}}));
  EXPECT_EQ(map(std::string(14, 'C')), (std::map<std::string, int>{{"CM", 3}, {"CT", 1}, {"CT2", 1}}));
  EXPECT_EQ(map("CCCCCCCCCCCCCC(=O)OCC(COP(=O)([O-])OCC[N+](C)(C)C)OC(=O)CCCCCCCCCCCCC"),
            (std::map<std::string, int>{{"CM", 6}, {"CT", 2}, {"EST1", 1}, {"EST2", 1}, {"GL", 1}, {"NC", 1}, {"PH", 1}}));
  // the LAMMPS SDK examples' own topologies: C12E8 (OA, 8 EO, 3 CM, CT2) and SDS (SO4, 3 CM, CT)
  EXPECT_EQ(map("CCCCCCCCCCCCOCCOCCOCCOCCOCCOCCOCCOCCO"), (std::map<std::string, int>{{"CM", 3}, {"CT2", 1}, {"EO", 8}, {"OA", 1}}));
  EXPECT_EQ(map("CCCCCCCCCCCCOS(=O)(=O)[O-]"), (std::map<std::string, int>{{"CM", 3}, {"CT", 1}, {"SO4", 1}}));
}

// A structure in vacuum sees every pair within the cut-off, however wide it is (the neighbour bins are half the list
// radius wide, so partners can sit two bins apart): the same energy as in a periodic box too large for images
TEST(FieldForms, VacuumPairsMatchALargePeriodicBox) {
  System s;
  for (int k = 0; k < 30; ++k) {
    Atom a;
    a.element = 18;
    a.name = "Ar";
    a.pos = {4.5 * k, 0.3 * (k % 3), 0};
    s.atoms.push_back(a);
  }
  const ForceField ff = assign_uff(s);
  EnergyOptions o;
  o.tail = false;
  o.coulomb = false;
  std::vector<double> x = flat(s), f, g;
  Evaluator vac(ff, o);
  const double ev = vac.compute(x, Cell{}, f).total();
  Cell big;
  big.origin = {-100, -100, -100};
  big.a = {400, 0, 0};
  big.b = {0, 400, 0};
  big.c = {0, 0, 400};
  Evaluator per(ff, o);
  const double ep = per.compute(x, big, g).total();
  EXPECT_NEAR(ev, ep, 1e-9 * std::max(1.0, std::fabs(ep)));
  for (size_t k = 0; k < f.size(); ++k) EXPECT_NEAR(f[k], g[k], 1e-9) << k;
}

// Martini 3 small molecules from all-atom structures, by graph: CAPS's own toluene and anthracene (its atom names, not
// CHARMM's) become TOLU and ANTR; anthracene's sites built on a site (R1 and R7 on R4, itself between R2 and R6) are
// placed in order and hand their forces back in reverse: forces and virial are the energy's derivatives
TEST(CoarseGrained, Martini3SmallMoleculesByGraph) {
  const std::string dir = std::string(CAPS_SOURCE_DIR) + "/data/";
  const FFDef m3 = load_forcefield(dir + "forcefields/martini3.json");
  System box;
  box.cell.a = {40, 0, 0};
  box.cell.b = {0, 40, 0};
  box.cell.c = {0, 0, 40};
  int k = 0;
  for (const char* smi : {"Cc1ccccc1", "c1ccc2cc3ccccc3cc2c1", "c1ccc2[nH]ccc2c1"}) {
    System s = build_molecule(smi).system;
    int64_t mol = ++k;
    const uint32_t off = uint32_t(box.atoms.size());
    for (auto a : s.atoms) {
      a.pos = a.pos + Vec3{10.0 * k, 12.0 * k, 20};
      a.mol = mol;
      box.atoms.push_back(a);
    }
    for (auto b : s.bonds) box.bonds.push_back({b.i + off, b.j + off, b.order});
  }
  box.bonds_from_file = true;
  std::vector<std::string> un, notes;
  const System cg = martini3_small_molecules(box, dir + "martini/martini3-small-molecules.json", 1e6, &un, &notes);
  EXPECT_TRUE(un.empty());
  ASSERT_FALSE(notes.empty());
  EXPECT_NE(notes[0].find("1 ANTR, 1 INDO, 1 TOLU"), std::string::npos) << notes[0];
  EXPECT_EQ(cg.atoms.size(), 3u + 7u + 5u);
  ASSERT_TRUE(cg.topology);
  EXPECT_EQ(cg.topology->vsites.size(), 3u + 1u);   // ANTR's three, INDO's one
  // the same through the force field's typing (all-atom in, beads out)
  System s = box;
  std::string ch = "auto";
  const std::string note = prepare_for_forcefield(s, m3, ch);
  EXPECT_EQ(s.atoms.size(), cg.atoms.size()) << note;
  check_cg_forces(m3, cg, "keep");
}

// A GROMACS topology read back: the .itp CAPS writes for Martini 3 small molecules (toluene, anthracene with its virtual
// sites) gives the same explicit topology again (every bond and constraint, angle, dihedral and exclusion)
TEST(GromacsTopology, ReadsMartini3ItpBack) {
  const std::string dir = std::string(CAPS_SOURCE_DIR) + "/data/";
  System box;
  box.cell.a = {40, 0, 0};
  box.cell.b = {0, 40, 0};
  box.cell.c = {0, 0, 40};
  int k = 0;
  for (const char* smi : {"Cc1ccccc1", "c1ccc2cc3ccccc3cc2c1"}) {
    System s = build_molecule(smi).system;
    int64_t mol = ++k;
    const uint32_t off = uint32_t(box.atoms.size());
    for (auto a : s.atoms) {
      a.pos = a.pos + Vec3{10.0 * k, 12.0 * k, 20};
      a.mol = mol;
      box.atoms.push_back(a);
    }
    for (auto b : s.bonds) box.bonds.push_back({b.i + off, b.j + off, b.order});
  }
  box.bonds_from_file = true;
  const System cg = martini3_small_molecules(box, dir + "martini/martini3-small-molecules.json");
  ASSERT_TRUE(cg.topology);
  const std::string path = (std::filesystem::temp_directory_path() / "caps_m3_readback.itp").string();
  {
    std::ofstream f(path);
    f << martini3_itp(cg);
  }
  const System T = read_gromacs_topology(path);
  std::filesystem::remove(path);
  ASSERT_TRUE(T.topology);
  EXPECT_EQ(T.atoms.size(), cg.atoms.size());
  auto bonds = [](const ExplicitTopology& t) {
    std::multiset<std::tuple<uint32_t, uint32_t, long, long>> s;
    for (const auto& b : t.bonds) s.insert({std::min(b.i, b.j), std::max(b.i, b.j), std::lround(b.r0 * 100), std::lround(b.k * 100)});
    return s;
  };
  EXPECT_EQ(bonds(*T.topology), bonds(*cg.topology));
  EXPECT_EQ(T.topology->angles.size(), cg.topology->angles.size());
  EXPECT_EQ(T.topology->dihedrals.size(), cg.topology->dihedrals.size());
  EXPECT_EQ(T.topology->vsites.size(), cg.topology->vsites.size());
}

// GROMACS's combined bending–torsion (dihedral function 11, Martini 3 polymer backbones): k sin³θ1 sin³θ2 Σ a_n cos^n φ,
// energies as GROMACS 2026 gives them for six random conformations (k 1.30929464, a 6.51130136 −4.6419125 10.24767684
// −7.64602475 12.3432375 kJ/mol), then forces and virial against finite differences, an explicit LJ pair with charges too
TEST(FieldForms, CombinedBendingTorsionAndExplicitPairs) {
  ForceField ff;
  ff.type_names = {"X"};
  ff.type_index = {0, 0, 0, 0};
  ff.lj = {{0, 0}};
  ff.mass = {72, 72, 72, 72};
  ff.charge = {0, 0, 0, 0};
  ff.excluded.assign(4, {});
  const double k = 1.30929464, a[5] = {6.51130136, -4.6419125, 10.24767684, -7.64602475, 12.3432375};
  CbtTorsion t{0, 1, 2, 3, {}};
  for (int n = 0; n < 5; ++n) t.a[n] = k * a[n] / 4.184;
  ff.cbt.push_back(t);
  EnergyOptions o;
  o.tail = false;
  o.coulomb = false;
  Cell none;
  const double frames[6][12] = {{2.0, 2.0, 2.0, 1.734, 2.045, 1.868, 1.801, 2.126, 1.587, 1.573, 2.284, 1.474},
                                {2.0, 2.0, 2.0, 1.858, 2.264, 1.984, 2.135, 2.245, 2.099, 1.936, 2.322, 2.309},
                                {2.0, 2.0, 2.0, 2.023, 2.244, 2.173, 1.769, 2.394, 2.226, 1.674, 2.17, 2.401},
                                {2.0, 2.0, 2.0, 1.981, 2.15, 2.259, 2.114, 2.411, 2.194, 2.284, 2.38, 2.44},
                                {2.0, 2.0, 2.0, 2.172, 1.818, 1.835, 2.017, 2.072, 1.8, 2.178, 1.819, 1.809},
                                {2.0, 2.0, 2.0, 1.834, 1.783, 2.124, 1.89, 2.052, 2.245, 2.063, 2.196, 2.443}};
  const double gromacs_kj[6] = {33.835560, 7.251500, 5.617820, 21.503141, 0.011181, 1.107190};
  {
    Evaluator ev(ff, o);
    for (int f = 0; f < 6; ++f) {
      std::vector<double> x(frames[f], frames[f] + 12), g;
      for (double& v : x) v *= 10;   // nm → Å
      EXPECT_NEAR(ev.compute(x, none, g).total() * 4.184, gromacs_kj[f], 2e-5 * std::max(1.0, gromacs_kj[f])) << f;
    }
  }
  // forces and virial: the CBT term and a charged explicit pair between the ends, in a periodic cell
  ff.lj_pairs.push_back({0, 3, 0.8, 4.6});
  ff.charge = {0.5, 0, 0, -1};
  ff.coul14 = 1.0;
  o.coulomb = true;
  Cell c;
  c.a = {40, 0, 0};
  c.b = {0, 40, 0};
  c.c = {0, 0, 40};
  Evaluator ev(ff, o);
  std::vector<double> x = {20, 20, 20, 17.3, 20.4, 18.7, 18.0, 21.3, 15.9, 15.7, 22.8, 14.7}, f, g;
  const EnergyTerms et = ev.compute(x, c, f);
  for (size_t kk = 0; kk < x.size(); ++kk) {
    std::vector<double> xp = x, xm = x;
    xp[kk] += 1e-5;
    xm[kk] -= 1e-5;
    const double fd = -(ev.compute(xp, c, g).total() - ev.compute(xm, c, g).total()) / 2e-5;
    EXPECT_NEAR(f[kk], fd, 1e-4 * std::max(1.0, std::fabs(fd))) << kk;
  }
  const auto d = strain_derivative(ev, x, c);
  for (int v = 0; v < 6; ++v) EXPECT_NEAR(et.w[v], d[v], 1e-4 * std::max(1.0, std::fabs(et.virial))) << "virial " << v;
}

// CVFF, PCFF and COMPASS from their .frc files (bench/ff/convert_frc.py; msi2lmp's assignment, checked against msi2lmp
// by bench/ff/check_msi2lmp.py): a polystyrene melt is typed and fully parameterised, the class II files with their
// cross terms and Wilson / angle-angle out-of-plane terms, CVFF with its improper torsions; each keeps its own LAMMPS
// styles for the native export.
TEST(FFDef, FrcConversionsParameteriseFully) {
  const System melt = read_lammps_data(std::string(CAPS_SOURCE_DIR) + "/samples/ps_melt.data");
  for (const char* id : {"pcff-frc", "cvff-frc", "compass-frc"}) {
    System s = melt;
    const FFDef def = load_forcefield(std::string(CAPS_SOURCE_DIR) + "/data/forcefields/" + id + ".json");
    EXPECT_EQ(def.equivalence, "fallback") << id;
    EXPECT_EQ(def.oop_scheme, "msi2lmp") << id;
    std::string ch = "types";
    prepare_for_forcefield(s, def, ch);
    const TypingResult tr = assign_types(s, def);
    ASSERT_EQ(tr.untyped, 0) << id;
    ParamReport rep;
    const ForceField ff = parameterize(s, def, tr.types, ch, &rep, false);
    EXPECT_TRUE(rep.missing.empty()) << id;
    if (std::string(id) == "cvff-frc") {
      EXPECT_EQ(ff.native_dihedral, "harmonic");
      EXPECT_EQ(ff.native_improper, "cvff");
      EXPECT_FALSE(ff.impropers.empty());   // the aromatic carbons' out-of-plane torsions
      EXPECT_TRUE(ff.dihedrals2.empty());
    } else {
      EXPECT_EQ(ff.native_dihedral, "class2") << id;
      ASSERT_FALSE(ff.dihedrals2.empty()) << id;
      ASSERT_FALSE(ff.impropers2.empty()) << id;
      bool mbt = false, bb = false, aa = false;
      for (const auto& d : ff.dihedrals2) mbt = mbt || d.mbt[0] != 0 || d.mbt[1] != 0 || d.mbt[2] != 0;
      for (const auto& a : ff.angles2) bb = bb || a.bb_m != 0;
      // the angle-angle reference angles in LAMMPS improper class2's order: θ1 = i-j-k, θ2 = i-j-l, θ3 = k-j-l, each the θ0
      // of that angle's own term (msi2lmp swaps θ2 and θ3; with them swapped liquid decane's improper energy came out
      // 19.59 kcal/mol in LAMMPS against 6.91)
      std::map<std::array<uint32_t, 3>, double> t0;
      for (const auto& a : ff.angles2) t0[{a.j, std::min(a.i, a.k), std::max(a.i, a.k)}] = a.theta0;
      auto at = [&](uint32_t a, uint32_t v, uint32_t b) {
        auto it = t0.find({v, std::min(a, b), std::max(a, b)});
        return it == t0.end() ? -1.0 : it->second;
      };
      size_t distinct = 0;
      for (const auto& q : ff.impropers2) {
        aa = aa || q.m1 != 0 || q.m2 != 0 || q.m3 != 0;
        EXPECT_NEAR(q.theta1, at(q.i, q.j, q.k), 1e-12) << id;
        EXPECT_NEAR(q.theta2, at(q.i, q.j, q.l), 1e-12) << id;
        EXPECT_NEAR(q.theta3, at(q.k, q.j, q.l), 1e-12) << id;
        distinct += std::fabs(q.theta2 - q.theta3) > 1e-6;
      }
      EXPECT_GT(distinct, 0u) << id;   // the check has terms where the order matters
      EXPECT_TRUE(mbt && bb && aa) << id;
    }
  }
}

// DREIDING's torsion rules (Mayo, Olafson, Goddard 1990), which give every torsion of the paper's library: polystyrene's backbone CH–ring
// carbon torsions are case (b), V = 1 kcal/mol, n = 6, φ0 = 0, over the 2 × 3 torsions about the bond; the export keeps
// DREIDING's own LAMMPS styles and special_bonds keyword.
TEST(FFDef, DreidingTorsionRulesFillWhatTheFileLacks) {
  System s = read_lammps_data(std::string(CAPS_SOURCE_DIR) + "/samples/ps_melt.data");
  const FFDef def = load_forcefield(std::string(CAPS_SOURCE_DIR) + "/data/forcefields/dreiding1990.json");
  EXPECT_EQ(def.torsion_rules, "dreiding1990");
  std::string ch = "gasteiger";
  prepare_for_forcefield(s, def, ch);
  const TypingResult tr = assign_types(s, def);
  ASSERT_EQ(tr.untyped, 0);
  ParamReport rep;
  const ForceField ff = parameterize(s, def, tr.types, ch, &rep, false);
  EXPECT_TRUE(rep.missing.empty());
  EXPECT_GT(rep.used["dihedral DREIDING rule (b)"], 0);
  bool found = false;
  for (const auto& d : ff.dihedrals)
    if (ff.atom_type[d.j] == "C_3" && ff.atom_type[d.k] == "C_R" && d.n == 6) {
      EXPECT_NEAR(d.v, 0.5 / 6, 1e-12);
      EXPECT_NEAR(d.delta, M_PI, 1e-12);   // ½V[1 − cos 6φ] = ½V[1 + cos(6φ − 180°)]
      found = true;
    }
  EXPECT_TRUE(found);
  EXPECT_EQ(ff.native_dihedral, "charmm");
  EXPECT_EQ(ff.native_special, "dreiding");
}

// CHARMM's force switch (lj/charmmfsw): two Lennard-Jones atoms, ε 0.1, σ 3.5, switched from 10 to 12 Å. Below 10 Å the
// plain 12-6 plus LAMMPS's constant, continuous at 10 Å, zero at 12 Å, and the force the exact derivative throughout.
TEST(FieldForms, CharmmForceSwitch) {
  ForceField ff;
  ff.type_names = {"X"};
  ff.type_index = {0, 0};
  ff.lj = {{0.1, 3.5}};
  ff.mass = {12, 12};
  ff.charge = {0, 0};
  ff.excluded.assign(2, {});
  ff.lj_fsw = true;
  ff.lj_inner = 10;
  ff.cutoff = 12;
  EnergyOptions o;
  o.coulomb = false;
  Cell none;
  Evaluator ev(ff, o);
  auto energy = [&](double r, double* force = nullptr) {
    std::vector<double> x = {0, 0, 0, r, 0, 0}, g;
    const double e = ev.compute(x, none, g).total();
    if (force) *force = g[3];
    return e;
  };
  const double A = 4 * 0.1 * std::pow(3.5, 12), B = 4 * 0.1 * std::pow(3.5, 6);
  EXPECT_NEAR(energy(5), A * (std::pow(5, -12) - std::pow(120, -6)) - B * (std::pow(5, -6) - std::pow(120, -3)), 1e-12);   // LAMMPS's form
  EXPECT_NEAR(energy(10 - 1e-7), energy(10 + 1e-7), 1e-9);
  EXPECT_NEAR(energy(12 - 1e-6), 0, 1e-12);
  for (double r : {4.0, 9.5, 10.5, 11.0, 11.9}) {
    double f = 0;
    energy(r, &f);
    const double h = 1e-5, dedr = (energy(r + h) - energy(r - h)) / (2 * h);
    EXPECT_NEAR(f, -dedr, 1e-7 * std::max(1.0, std::fabs(dedr))) << r;
  }
}

// DREIDING's hydrogen bond (hbond/dreiding/lj): donor O, its hydrogen, an acceptor O; ε 4, σ 2.75, cos⁴θ, switched from 6
// to 6.5 Å, only past 90°. Forces and virial are the exact derivatives, inside and in the switching region; no term
// below 90°.
TEST(FieldForms, DreidingHydrogenBond) {
  ForceField ff;
  ff.type_names = {"O_3_hd", "H_HB", "O_2_ha"};
  ff.type_index = {0, 1, 2};
  ff.lj = {{0, 0}, {0, 0}, {0, 0}};
  ff.mass = {16, 1, 16};
  ff.charge = {0, 0, 0};
  ff.bonds.push_back({0, 1, 0, 0.97});   // no bond energy: only the hydrogen's link to its donor
  ff.excluded.assign(3, {});
  ff.excluded[0] = {1};
  ff.excluded[1] = {0};
  ff.hbond.hyd = {{1}, {}, {}};
  ff.hbond.acceptor = {0, 0, 1};
  ff.hbond.param[{0, 2}] = {4.0, 2.75, 4};
  ff.hbond.htype[{0, 2}] = 1;
  ff.hbond.cos_cut = 0;
  EnergyOptions o;
  o.coulomb = false;
  o.tail = false;
  Cell none;
  Evaluator ev(ff, o);
  auto energy = [&](const std::vector<double>& x, std::vector<double>* f = nullptr) {
    std::vector<double> g;
    const double e = ev.compute(x, none, g).total();
    if (f) *f = g;
    return e;
  };
  // D at the origin, H along +x bent a little, A along +x at r
  for (double r : {2.8, 3.5, 6.2, 6.4}) {
    const std::vector<double> x = {0, 0, 0, 0.97, 0.12, 0.05, r, -0.3, 0.2};
    std::vector<double> f;
    const double e = energy(x, &f);
    EXPECT_LT(e, r < 6 ? -0.0 : 1.0) << r;
    for (size_t k = 0; k < 9; ++k) {
      auto xp = x, xm = x;
      const double h = 1e-6;
      xp[k] += h;
      xm[k] -= h;
      EXPECT_NEAR(f[k], -(energy(xp) - energy(xm)) / (2 * h), 1e-6 * std::max(1.0, std::fabs(f[k]))) << r << " " << k;
    }
  }
  // the hydrogen on the far side (D–H···A below 90°): no term
  EXPECT_NEAR(energy({0, 0, 0, -0.97, 0, 0, 2.9, 0, 0}), 0, 1e-12);
}

// IFF's silica and metals by rule (bench/ff/convert_frc.py IFF_EXTRA): orthosilicic acid's Si is sc4, its silanol O oc24
// and H hoy, with the charges IFF states (+1.1, -0.675, +0.4) from bond increments; quartz is sc4 / oc23 (+1.1 / -0.55);
// copper's perceived Cu-Cu bonds are dropped (IFF's metals are Lennard-Jones atoms). Both IFF variants.
TEST(FFDef, IffSilicaAndMetalsByRule) {
  for (const char* id : {"iff-pcff", "iff-cvff"}) {
    const FFDef def = load_forcefield(std::string(CAPS_SOURCE_DIR) + "/data/forcefields/" + id + ".json");
    System s = build_molecule("O[Si](O)(O)O").system;
    const TypingResult tr = assign_types(s, def);
    ASSERT_EQ(tr.untyped, 0) << id;
    ParamReport rep;
    const ForceField ff = parameterize(s, def, tr.types, "types", &rep, false);
    EXPECT_TRUE(rep.missing.empty()) << id;
    for (size_t i = 0; i < s.atoms.size(); ++i) {
      const int z = s.atoms[i].element;
      EXPECT_EQ(tr.types[i], z == 14 ? "sc4" : z == 8 ? "oc24" : "hoy") << id;
      EXPECT_NEAR(ff.charge[i], z == 14 ? 1.1 : z == 8 ? -0.675 : 0.4, 1e-12) << id;
    }
    System q = read_cif(std::string(CAPS_SOURCE_DIR) + "/data/crystals/alpha-quartz.cif");
    q.bonds = perceive_bonds(q);
    ASSERT_EQ(q.bonds.size(), 2 * q.atoms.size() / 3 * 2) << id;   // four Si-O per Si
    const TypingResult tq = assign_types(q, def);
    const ForceField fq = parameterize(q, def, tq.types, "types", nullptr, false);
    for (size_t i = 0; i < q.atoms.size(); ++i) {
      EXPECT_EQ(tq.types[i], q.atoms[i].element == 14 ? "sc4" : "oc23") << id;
      EXPECT_NEAR(fq.charge[i], q.atoms[i].element == 14 ? 1.1 : -0.55, 1e-12) << id;
    }
    System cu = read_cif(std::string(CAPS_SOURCE_DIR) + "/data/crystals/copper.cif");
    cu.bonds = perceive_bonds(cu);
    std::string ch = "types";
    const std::string note = prepare_for_forcefield(cu, def, ch);
    EXPECT_TRUE(cu.bonds.empty()) << id;
    EXPECT_NE(note.find("does not bond Cu"), std::string::npos) << note;
    EXPECT_EQ(assign_types(cu, def).types[0], "Cu") << id;
  }
}

// DREIDING (Mayo, Olafson & Goddard 1990): acetanilide (an aramid's amide on a ring) is typed and every term is found
TEST(FFDef, BondOrderVariantsKeepTorsions) {
  System s = build_molecule("CC(=O)Nc1ccccc1").system;
  const FFDef def = load_forcefield(std::string(CAPS_SOURCE_DIR) + "/data/forcefields/dreiding1990.json");
  std::string ch = "gasteiger";
  if (needs_prepare(def)) prepare_for_forcefield(s, def, ch);
  const TypingResult tr = assign_types(s, def);
  ASSERT_EQ(tr.untyped, 0);
  ParamReport rep;
  (void)parameterize(s, def, tr.types, ch, &rep, false);
  EXPECT_TRUE(rep.missing.empty()) << (rep.missing.empty() ? "" : rep.missing.front());
}

// AMBER frcmod: harmonic terms as written, multi-term dihedrals (negative periodicity), IDIVF, impropers, R*/2 → σ
TEST(ParamImport, Frcmod) {
  const auto path = (std::filesystem::temp_directory_path() / "caps_test.frcmod").string();
  {
    std::ofstream o(path);
    o << "remark\nMASS\nzz 12.01         0.878               my carbon\n\nBOND\nzz-c3  303.1   1.5350\n\nANGL\nc3-zz-hc   46.37   110.05\n\n"
         "DIHE\nX -zz-c3-X    9    1.400       0.000           3.000\nhc-zz-c3-hc   1    0.150       0.000          -3.000\nhc-zz-c3-hc   1    0.250     180.000           1.000\n\n"
         "IMPR\nX -X -zz-o    1.1      180.0         2.0\n\nNONB\n  zz          1.9080  0.1094\n\n";
  }
  const caps::FFDef d = caps::import_frcmod(path);
  ASSERT_EQ(d.types.size(), 1u);
  EXPECT_EQ(d.types[0].element, 6);
  ASSERT_EQ(d.bonds.size(), 1u);
  EXPECT_EQ(d.bonds[0].match, (std::vector<std::string>{"zz", "c3"}));
  EXPECT_DOUBLE_EQ(d.bonds[0].params[0], 303.1);
  EXPECT_DOUBLE_EQ(d.angles[0].params[1], 110.05);
  ASSERT_EQ(d.dihedrals.size(), 2u);
  EXPECT_EQ(d.dihedrals[0].match[0], "*");
  EXPECT_NEAR(d.dihedrals[0].params[1], 1.4 / 9, 1e-12);                        // PK / IDIVF
  EXPECT_EQ(d.dihedrals[1].params, (std::vector<double>{2, 0.15, 3, 0, 0.25, 1, 180}));   // two terms in one rule
  EXPECT_EQ(d.impropers[0].params, (std::vector<double>{1.1, -1, 2}));
  EXPECT_NEAR(d.pairs[0].params[1], 3.39966950842, 1e-9);                        // 2 R*/2 / 2^(1/6)
  EXPECT_DOUBLE_EQ(d.pairs[0].params[0], 0.1094);
}

// GROMACS [ *types ]: GROMACS's own AMBER ff14SB gives back the AMBER numbers (kJ/mol, nm → kcal/mol, Å; ½k → K)
TEST(ParamImport, GromacsAmberTypes) {
  const std::string dir = "/opt/homebrew/share/gromacs/top/amber14sb.ff/";
  if (!std::filesystem::exists(dir + "ffbonded.itp")) GTEST_SKIP() << "GROMACS's amber14sb.ff not installed";
  const caps::FFDef b = caps::import_gromacs_params(dir + "ffbonded.itp");
  auto find = [](const std::vector<caps::FFRule>& v, std::vector<std::string> m) {
    for (const auto& r : v) if (r.match == m) return &r;
    return static_cast<const caps::FFRule*>(nullptr);
  };
  const auto* ctct = find(b.bonds, {"CT", "CT"});
  ASSERT_NE(ctct, nullptr);
  EXPECT_NEAR(ctct->params[0], 310.0, 0.05);
  EXPECT_NEAR(ctct->params[1], 1.526, 1e-9);
  const auto* hch = find(b.angles, {"HC", "CT", "HC"});
  ASSERT_NE(hch, nullptr);
  EXPECT_NEAR(hch->params[0], 35.0, 0.01);
  EXPECT_NEAR(hch->params[1], 109.5, 1e-9);
  const auto* xctctx = find(b.dihedrals, {"*", "CT", "CT", "*"});
  ASSERT_NE(xctctx, nullptr);
  EXPECT_NEAR(xctctx->params[1], 1.4 / 9, 1e-4);
  EXPECT_EQ(xctctx->params[2], 3);
  const caps::FFDef nb = caps::import_gromacs_params(dir + "ffnonbonded.itp");
  const auto* ct = find(nb.pairs, {"CT"});
  ASSERT_NE(ct, nullptr);
  EXPECT_NEAR(ct->params[1], 3.39967, 1e-4);
  EXPECT_NEAR(ct->params[0], 0.1094, 1e-4);
}

// CHARMM19 / GROMOS style: one torsion per central bond (the end of highest atomic number; every H where only H is
// there, one where torsion_single_h names the centre: an amide NH2); every 1-4 pair kept
TEST(FFDef, OneTorsionPerBond) {
  const std::string path = (std::filesystem::temp_directory_path() / "caps_one_per_bond.json").string();
  std::ofstream(path) << R"({"format": "caps-forcefield", "format_version": 1, "name": "one per bond", "units": "real",
 "styles": {"pair": "lj/cut/coul/long", "bond": "harmonic", "angle": "harmonic", "dihedral": "fourier", "improper": "harmonic"},
 "mixing": "arithmetic", "special_lj": [0, 0, 1], "special_coul": [0, 0, 1], "cutoff": 10, "torsion_selection": "one_per_bond",
 "torsion_single_h": ["[#7X3;$(*-[#6]=[#8X1])]"],
 "atom_types": [{"name": "CT", "element": "C", "mass": 12.011}, {"name": "C", "element": "C", "mass": 12.011},
                {"name": "O", "element": "O", "mass": 15.999}, {"name": "N", "element": "N", "mass": 14.007},
                {"name": "H", "element": "H", "mass": 1.008}, {"name": "HC", "element": "H", "mass": 1.008}],
 "pairs": [{"name": "*", "match": ["*"], "params": [0.1, 3.0]}],
 "bonds": [{"name": "X-X", "match": ["*", "*"], "params": [300, 1.2]}],
 "angles": [{"name": "X-X-X", "match": ["*", "*", "*"], "params": [50, 115]}],
 "dihedrals": [{"name": "X-C-N-X", "match": ["*", "C", "N", "*"], "params": [1, 8.2, 2, 180]},
               {"name": "X-CT-C-X", "match": ["*", "CT", "C", "*"], "params": [1, 1.0, 3, 0]}],
 "impropers": []})";
  FFDef def = load_forcefield(path);
  // acetamide CH3-C(=O)-NH2
  System s;
  const std::vector<std::pair<int, Vec3>> at = {{6, {0, 0, 0}}, {6, {1.52, 0, 0}}, {8, {2.1, 1.05, 0}}, {7, {2.2, -1.15, 0}},
                                                {1, {3.2, -1.15, 0}}, {1, {1.7, -2.05, 0}}, {1, {-0.4, 1.0, 0}},
                                                {1, {-0.4, -0.5, 0.87}}, {1, {-0.4, -0.5, -0.87}}};
  for (const auto& [z, p] : at) {
    Atom a;
    a.element = z;
    a.pos = p;
    s.atoms.push_back(a);
  }
  s.bonds = {{0, 1, 1}, {1, 2, 2}, {1, 3, 1}, {3, 4, 1}, {3, 5, 1}, {0, 6, 1}, {0, 7, 1}, {0, 8, 1}};
  const std::vector<std::string> types = {"CT", "C", "O", "N", "H", "H", "HC", "HC", "HC"};
  ParamReport rep;
  ForceField ff = parameterize(s, def, types, "keep", &rep, false);
  // CT-C: the three methyl H with O; C-N: O with the first N-H only
  ASSERT_EQ(ff.dihedrals.size(), 4u);
  int amide = 0;
  for (const auto& d : ff.dihedrals) {
    const bool cn = (d.j == 1 && d.k == 3) || (d.j == 3 && d.k == 1);
    if (cn) {
      ++amide;
      EXPECT_TRUE((d.i == 2 && d.l == 4) || (d.i == 4 && d.l == 2));
    } else {
      EXPECT_TRUE(d.i == 2 || d.l == 2);
    }
  }
  EXPECT_EQ(amide, 1);
  def.torsions_one_per_bond = false;   // every torsion: 3 x 2 about CT-C, 2 x 2 about C-N
  ff = parameterize(s, def, types, "keep", &rep, false);
  EXPECT_EQ(ff.dihedrals.size(), 10u);
  std::filesystem::remove(path);
}

// ClayFF from its paper (Cygan, Liang & Kalinichev, J. Phys. Chem. B 108, 1255 (2004); bench/ff/make_clayff_cygan2004.py):
// every charge, Do and Ro of Table 1 and every term of Table 2 (p 1259), as an independent transcription here; Do, Ro as
// Lennard-Jones epsilon = Do, sigma = Ro / 2^(1/6), Lorentz-Berthelot (Eqs. 4-5). Brucite (data/crystals) is typed mgh /
// oh / ho with only its O-H bonds kept and three Mg-O-H bends per hydroxyl through contacts; SPC water o* / h* with its
// H-O-H bend; pyrophyllite's formula unit Al2Si4O10(OH)2 and brucite are neutral with the paper's charges.
TEST(FFDef, ClayffCygan2004AsPublished) {
  const FFDef def = load_forcefield(std::string(CAPS_SOURCE_DIR) + "/data/forcefields/clayff-cygan2004.json");
  EXPECT_EQ(def.mixing, "arithmetic");
  struct Row { const char* sym; double q, D, R; };
  const Row table1[] = {{"h*", 0.41, 0, 0}, {"ho", 0.425, 0, 0}, {"o*", -0.82, 0.1554, 3.5532}, {"oh", -0.95, 0.1554, 3.5532},
                        {"ob", -1.05, 0.1554, 3.5532}, {"obos", -1.1808, 0.1554, 3.5532}, {"obts", -1.1688, 0.1554, 3.5532},
                        {"obss", -1.2996, 0.1554, 3.5532}, {"ohs", -1.0808, 0.1554, 3.5532}, {"st", 2.1, 1.8405e-6, 3.7064},
                        {"ao", 1.575, 1.3298e-6, 4.7943}, {"at", 1.575, 1.8405e-6, 3.7064}, {"mgo", 1.36, 9.0298e-7, 5.9090},
                        {"mgh", 1.05, 9.0298e-7, 5.9090}, {"cao", 1.36, 5.0298e-6, 6.2484}, {"cah", 1.05, 5.0298e-6, 6.2428},
                        {"feo", 1.575, 9.0298e-6, 5.5070}, {"lio", 0.525, 9.0298e-6, 4.7257}, {"Na", 1, 0.1301, 2.6378},
                        {"K", 1, 0.1, 3.7423}, {"Cs", 1, 0.1, 4.3002}, {"Ca", 2, 0.1, 3.2237}, {"Ba", 2, 0.047, 4.2840},
                        {"Cl", -1, 0.1001, 4.9388}};
  ASSERT_EQ(def.types.size(), std::size(table1));
  for (const Row& r : table1) {
    const FFType* t = def.type(r.sym);
    ASSERT_NE(t, nullptr) << r.sym;
    EXPECT_DOUBLE_EQ(t->charge, r.q) << r.sym;
    const auto p = std::find_if(def.pairs.begin(), def.pairs.end(), [&](const FFRule& x) { return x.name == r.sym; });
    ASSERT_NE(p, def.pairs.end()) << r.sym;
    EXPECT_DOUBLE_EQ(p->params[0], r.D) << r.sym;
    EXPECT_NEAR(p->params[1] * std::pow(2.0, 1.0 / 6), r.R, 1e-12) << r.sym;
  }
  ASSERT_EQ(def.bonds.size(), 3u);
  for (const auto& b : def.bonds) {
    EXPECT_DOUBLE_EQ(b.params[0], 554.1349) << b.name;
    EXPECT_DOUBLE_EQ(b.params[1], 1.0) << b.name;
  }
  ASSERT_EQ(def.angles.size(), 19u);   // water's H-O-H and Metal-O-H for nine metals and two hydroxyl oxygens
  for (const auto& a : def.angles) {
    EXPECT_DOUBLE_EQ(a.params[0], a.name == "h*-o*-h*" ? 45.7696 : 30.0) << a.name;
    EXPECT_DOUBLE_EQ(a.params[1], 109.47) << a.name;
  }
  auto q = [&](const char* sym) { return def.type(sym)->charge; };
  EXPECT_NEAR(2 * q("ao") + 4 * q("st") + 10 * q("ob") + 2 * (q("oh") + q("ho")), 0, 1e-12);   // pyrophyllite

  System s = supercell(read_cif(std::string(CAPS_SOURCE_DIR) + "/data/crystals/brucite.cif"), 4, 4, 3);
  s.bonds = perceive_bonds(s);
  std::string ch = "types";
  prepare_for_forcefield(s, def, ch);
  const TypingResult tr = assign_types(s, def);
  ASSERT_EQ(tr.untyped, 0);
  for (size_t i = 0; i < s.atoms.size(); ++i)
    EXPECT_EQ(tr.types[i], s.atoms[i].element == 12 ? "mgh" : s.atoms[i].element == 8 ? "oh" : "ho") << i;
  ParamReport rep;
  const ForceField ff = parameterize(s, def, tr.types, ch, &rep, false);
  EXPECT_TRUE(rep.missing.empty());
  EXPECT_EQ(ff.bonds.size(), 96u);          // the O-H bonds only
  EXPECT_EQ(ff.angles.size(), 3 * 96u);     // three Mg around each hydroxyl O
  for (const auto& a : ff.angles) {
    EXPECT_EQ(tr.types[a.i], "mgh");
    EXPECT_NEAR(a.kt, 30.0, 1e-12);
  }
  double total = 0;
  for (double c : ff.charge) total += c;
  EXPECT_NEAR(total, 0, 1e-9);

  System w = build_molecule("O").system;
  const TypingResult tw = assign_types(w, def);
  ASSERT_EQ(tw.untyped, 0);
  const ForceField fw = parameterize(w, def, tw.types, "types", nullptr, false);
  ASSERT_EQ(fw.angles.size(), 1u);
  EXPECT_NEAR(fw.angles[0].kt, 45.7696, 1e-12);
  EXPECT_EQ(tw.types[size_t(fw.angles[0].j)], "o*");
}

// Sanders, Leslie & Catlow's O-Si-O bends by distance: Si and O unbonded, the O shells within 1.9 A of each Si bent
TEST(FFDef, ThreeBodyBendsByDistanceOnShells) {
  const FFDef def = load_forcefield(std::string(CAPS_SOURCE_DIR) + "/data/forcefields/sanders1984-sio2.json");
  System q = read_cif(std::string(CAPS_SOURCE_DIR) + "/data/crystals/alpha-quartz.cif");
  q.bonds = perceive_bonds(q);
  std::string ch = "types";
  prepare_for_forcefield(q, def, ch);
  ASSERT_EQ(q.atoms.size(), 15u);   // 3 Si, 6 O cores, 6 O shells
  EXPECT_EQ(q.bonds.size(), 6u);    // the core-shell springs only
  const TypingResult tr = assign_types(q, def);
  ASSERT_EQ(tr.untyped, 0);
  ParamReport rep;
  const ForceField ff = parameterize(q, def, tr.types, "types", &rep, false);
  ASSERT_EQ(ff.angles.size(), 18u);   // six per tetrahedron
  for (const auto& a : ff.angles) {
    EXPECT_EQ(tr.types[a.j], "Si_s84");
    EXPECT_EQ(tr.types[a.i], "O_s84s");
    EXPECT_EQ(tr.types[a.k], "O_s84s");
    EXPECT_NEAR(a.kt, 2.097 * 23.060548, 1e-9);
  }
}

// Tilocca, de Leeuw & Cormack's SM1: the screened O-Si-O bend on the shells within 2.5 A, 1/2 k (θ − θ0)² exp(−(r1 + r2)/ρ);
// alpha-quartz's 18 bends sum to 0.147118 kcal/mol (an independent evaluation of the paper's eq. 4 from the same coordinates)
TEST(FFDef, ScreenedBendsOfTheShellModel) {
  const FFDef def = load_forcefield(std::string(CAPS_SOURCE_DIR) + "/data/forcefields/tilocca2006-sm.json");
  System q = read_cif(std::string(CAPS_SOURCE_DIR) + "/data/crystals/alpha-quartz.cif");
  q.bonds = perceive_bonds(q);
  std::string ch = "types";
  prepare_for_forcefield(q, def, ch);
  const TypingResult tr = assign_types(q, def);
  ASSERT_EQ(tr.untyped, 0);
  const ForceField ff = parameterize(q, def, tr.types, "types", nullptr, false);
  ASSERT_EQ(ff.angles_x.size(), 18u);
  EXPECT_TRUE(ff.angles.empty());
  for (const auto& a : ff.angles_x) EXPECT_EQ(a.form, 6);
  EnergyOptions o;
  o.tail = false;
  std::vector<double> f;
  const EnergyTerms e = Evaluator(ff, o).compute(flat(q), q.cell, f);
  EXPECT_NEAR(e.angle, 0.147118, 2e-6);
}

// Xiao, Edwards & Grater's carbonate: no C-O bond term and no exclusions; three O-C-O bends and one O1-C-O2-O3 plane term
// (K (1 - cos 2psi), the centre second, the O sorted by index) per carbonate, found by distance
TEST(FFDef, CarbonateTermsByDistanceWithoutBonds) {
  const FFDef def = load_forcefield(std::string(CAPS_SOURCE_DIR) + "/data/forcefields/xiao2011-caco3.json");
  System q = read_cif(std::string(CAPS_SOURCE_DIR) + "/data/crystals/calcite.cif");
  q.bonds = perceive_bonds(q);
  std::string ch = "types";
  prepare_for_forcefield(q, def, ch);
  EXPECT_TRUE(q.bonds.empty());
  const TypingResult tr = assign_types(q, def);
  ASSERT_EQ(tr.untyped, 0);
  const ForceField ff = parameterize(q, def, tr.types, "types", nullptr, false);
  ASSERT_EQ(ff.angles.size(), 18u);
  ASSERT_EQ(ff.impropers.size(), 6u);
  for (const auto& t : ff.impropers) {
    EXPECT_EQ(tr.types[t.j], "Cm_x");
    EXPECT_TRUE(t.i < t.k && t.k < t.l);
    EXPECT_EQ(t.n, 2);
    EXPECT_NEAR(t.v, 28.9 / 4.184, 1e-12);
    EXPECT_NEAR(t.delta, 3.14159265358979323846, 1e-12);   // v [1 + cos(2φ − π)] = v (1 − cos 2φ)
  }
}

// Schroder et al.'s hydroxyls: the O-H bond carries no term and no exclusion (the Morse is a pair within 1.4 A), yet the
// hydroxyl O and H are typed from it (the dropped bond kept for typing)
TEST(FFDef, DroppedBondsStillTypeAndPairCutOff) {
  const FFDef def = load_forcefield(std::string(CAPS_SOURCE_DIR) + "/data/forcefields/schroder1992-zeolite.json");
  System s = build_molecule("O[Si](O)(O)O").system;
  std::string ch = "types";
  prepare_for_forcefield(s, def, ch);
  EXPECT_TRUE(s.bonds.empty());
  EXPECT_EQ(s.typing_bonds.size(), 8u);
  const TypingResult tr = assign_types(s, def);
  ASSERT_EQ(tr.untyped, 0);
  for (size_t i = 0; i < s.atoms.size(); ++i)
    EXPECT_EQ(tr.types[i], s.atoms[i].element == 14 ? "Si_sz" : s.atoms[i].element == 8 ? "Oh_sz" : "H_sz");
  const ForceField ff = parameterize(s, def, tr.types, "types", nullptr, false);
  bool morse = false;
  for (const auto& [k, p] : ff.pair_func)
    if (p.form == 2) {
      morse = true;
      EXPECT_NEAR(p.rc, 1.4, 1e-12);
    }
  EXPECT_TRUE(morse);
  EXPECT_EQ(ff.angles.size(), 6u);
}

// Zeolite sorbates from their papers (bench/ff/make_zeolite_sorbates.py): O'Malley & Catlow's all-atom alkane neutral by
// its C types (CH3 -0.3, CH2 -0.2, CH -0.1); Plant et al.'s methanol and Peguin et al.'s HFA134a typed with their tables'
// charges; a Buckingham term with its own cut-off (10 A, under the 13 A global one) adds nothing beyond it, not -shift
TEST(FFDef, ZeoliteSorbatesFromTheirPapers) {
  const std::string dir = std::string(CAPS_SOURCE_DIR) + "/data/forcefields/";
  auto typed = [&](const std::string& ff, const std::string& smiles, std::map<std::string, int>& count) {
    const FFDef def = load_forcefield(dir + ff + ".json");
    System s = build_molecule(smiles).system;
    std::string ch = "types";
    prepare_for_forcefield(s, def, ch);
    const TypingResult tr = assign_types(s, def);
    EXPECT_EQ(tr.untyped, 0) << ff;
    for (const auto& t : tr.types) ++count[t];
    const ForceField f = parameterize(s, def, tr.types, "types", nullptr, false);
    double q = 0;
    for (double x : f.charge) q += x;
    EXPECT_NEAR(q, 0.0, 1e-9) << ff;
    return f;
  };
  std::map<std::string, int> c;
  typed("omalley2013-alkanes", "CC(C)CCC(C)(C)C", c);
  EXPECT_EQ(c["C3_om"], 5);
  EXPECT_EQ(c["C2_om"], 2);
  EXPECT_EQ(c["C1_om"], 1);
  EXPECT_EQ(c["C0_om"], 1);
  EXPECT_EQ(c["H_om"], 20);
  c.clear();
  const ForceField m = typed("plant2007-methanol", "CO", c);
  EXPECT_EQ(c["C_pm"], 1);
  EXPECT_EQ(c["Om_pm"], 1);
  EXPECT_EQ(c["Hm_pm"], 1);
  EXPECT_EQ(c["H_pm"], 3);
  ASSERT_EQ(m.dihedrals.size(), 3u);
  c.clear();
  typed("peguin2009-hfa134a", "FCC(F)(F)F", c);
  EXPECT_EQ(c["C3_hfa"], 1);
  EXPECT_EQ(c["C2_hfa"], 1);
  EXPECT_EQ(c["F_hfa"], 4);
  EXPECT_EQ(c["H_hfa"], 2);

  const FFDef def = load_forcefield(dir + "omalley2013-alkanes.json");
  System s;
  for (int k = 0; k < 2; ++k) {
    Atom a;
    a.element = 8;
    a.name = "O";
    a.pos = {11.0 * k, 0, 0};
    s.atoms.push_back(a);
  }
  const ForceField ff = parameterize(s, def, {"O_om", "O_om"}, "types", nullptr, false);
  EnergyOptions o;
  o.tail = false;
  o.coulomb = false;
  o.cutoff = 13;
  std::vector<double> f;
  EXPECT_NEAR(Evaluator(ff, o).compute(flat(s), Cell{}, f).vdw, 0.0, 1e-12);
}

// TraPPE-EH from the TraPPE database (bench/ff/make_trappe_eh.py): each model's own charges by whole-molecule rules, so
// phenol and pyridine are typed as themselves and neutral; an alkane's H is the C-H bond-centre site (C-M 0.55 A); a
// structure that is no TraPPE-EH model stays untyped; 1-4 pairs carry nothing
TEST(FFDef, TrappeEhModelsTypedAsThemselves) {
  const FFDef def = load_forcefield(std::string(CAPS_SOURCE_DIR) + "/data/forcefields/trappe-eh-umn.json");
  for (const std::string smi : {"Oc1ccccc1", "c1ccncc1", "CCCCCCCC", "CC(N)=O", "Clc1cccs1"}) {
    System s = build_molecule(smi).system;
    std::string ch = "types";
    prepare_for_forcefield(s, def, ch);
    const TypingResult tr = assign_types(s, def);
    ASSERT_EQ(tr.untyped, 0) << smi;
    const ForceField ff = parameterize(s, def, tr.types, "types", nullptr, false);
    double q = 0;
    for (double x : ff.charge) q += x;
    EXPECT_NEAR(q, 0.0, 1e-9) << smi;
    if (smi == "CCCCCCCC")
      for (const auto& b : ff.bonds)
        if (s.atoms[b.i].element == 1 || s.atoms[b.j].element == 1) EXPECT_NEAR(b.r0, 0.55, 1e-12);
    if (smi == "Clc1cccs1")   // the database's slip (the ring S typed as furan's O) corrected
      for (size_t i = 0; i < s.atoms.size(); ++i)
        if (s.atoms[i].element == 16) EXPECT_EQ(tr.types[i].substr(0, 3), "S7_");
  }
  System t = build_molecule("Cc1ccccc1").system;   // toluene: no TraPPE-EH model
  std::string ch = "types";
  prepare_for_forcefield(t, def, ch);
  EXPECT_GT(assign_types(t, def).untyped, 0);
}

// TraPPE-UA's transferable hydrocarbon terms (bench/ff/make_trappe.py, from Martin & Siepmann 1998 and 1999 and Wick,
// Martin & Siepmann 2000): polymer oligomers whose torsions no database molecule has are parameterised whole, with the
// papers' numbers; cis and trans C=C get their own harmonic term; a carbon in a three-membered ring stays untyped
TEST(FFDef, TrappeUaTransferableHydrocarbonTerms) {
  const FFDef def = load_forcefield(std::string(CAPS_SOURCE_DIR) + "/data/forcefields/trappe-ua-umn.json");
  constexpr double KB = 0.0019872043;   // kcal/mol/K
  const double kPi = std::acos(-1.0);
  struct Want { const char* t[4]; int n; double v, delta; };
  struct Case { const char* smi; std::vector<Want> tors; double k_ct, chi0; int n_ct; };
  const std::vector<Case> cases = {
      // polypropylene: CHx-(CH2)-(CH)-CHy, 1999 Table 2 (c3 441.27 K)
      {"CC(C)CC(C)CC(C)CC(C)C", {{{"CH3_2", "CH_4", "CH2_3", "CH_4"}, 3, 441.27 * KB, 0}}, 0, 0, 0},
      // polystyrene: CHx(aro)-C(aro)-(CH)-CHx, 2000 Table 4 (e0 167 K, e1 5 pi / 3), and the backbone CH-CH2-CH-C(aro)
      {"CC(c1ccccc1)CC(c1ccccc1)C", {{{"CH_88", "C_89", "CH_4", "CH2_3"}, 2, 167 * KB, 240},
                                      {{"C_89", "CH_4", "CH2_3", "CH_4"}, 1, 428.73 * KB, 0}}, 0, 0, 0},
      // polyisobutylene: CHx-(CH2)-(C)-CHy, 1999 Table 2 (c3 461.29 K)
      {"CC(C)(C)CC(C)(C)CC(C)(C)C", {{{"CH3_2", "C_5", "CH2_3", "C_5"}, 3, 461.29 * KB, 0}}, 0, 0, 0},
      // cis- and trans-1,4-polybutadiene: CHx-(CH2)-(CH)=CHy, 2000 Table 4 (c1 86.36 K); cis / trans about CH=CH
      {"CC/C=C\\CC/C=C\\CC", {{{"CH2_3", "CH2_3", "CH_30", "CH_30"}, 1, 86.36 * KB, 0}}, 24800.0 / 2 * KB, 0.0, 2},
      {"CC/C=C/CC/C=C/CC", {{{"CH2_3", "CH2_3", "CH_30", "CH_30"}, 1, 86.36 * KB, 0}}, 26800.0 / 2 * KB, kPi, 2},
  };
  for (const auto& c : cases) {
    System s = build_molecule(c.smi).system;
    std::string ch = "types";
    prepare_for_forcefield(s, def, ch);
    const TypingResult tr = assign_types(s, def);
    ASSERT_EQ(tr.untyped, 0) << c.smi;
    ParamReport rep;
    const ForceField ff = parameterize(s, def, tr.types, ch, &rep, false);
    EXPECT_TRUE(rep.missing.empty()) << c.smi;
    for (const auto& b : ff.bonds) {
      const bool dbl = tr.types[b.i].rfind("CH_30", 0) == 0 && tr.types[b.j].rfind("CH_30", 0) == 0;
      const bool aro = tr.types[b.i].find("_8") != std::string::npos && tr.types[b.j].find("_8") != std::string::npos;
      EXPECT_NEAR(b.r0, dbl ? 1.33 : aro ? 1.40 : 1.54, 1e-12) << c.smi << " " << tr.types[b.i] << "-" << tr.types[b.j];
    }
    for (const auto& w : c.tors) {
      int found = 0;
      for (const auto& d : ff.dihedrals) {
        const std::array<std::string, 4> t{tr.types[d.i], tr.types[d.j], tr.types[d.k], tr.types[d.l]};
        const bool fw = t[0] == w.t[0] && t[1] == w.t[1] && t[2] == w.t[2] && t[3] == w.t[3];
        const bool bw = t[3] == w.t[0] && t[2] == w.t[1] && t[1] == w.t[2] && t[0] == w.t[3];
        if ((fw || bw) && d.n == w.n) {
          EXPECT_NEAR(d.v, w.v, 1e-8) << c.smi << " " << w.t[0] << "-" << w.t[1] << "-" << w.t[2] << "-" << w.t[3];
          EXPECT_NEAR(std::fmod(d.delta * 180 / kPi + 360, 360.0), w.delta, 1e-6) << c.smi;
          ++found;
        }
      }
      EXPECT_GT(found, 0) << c.smi << " " << w.t[0] << "-" << w.t[1] << "-" << w.t[2] << "-" << w.t[3];
    }
    ASSERT_EQ(int(ff.impropers_harmonic.size()), c.n_ct) << c.smi;
    for (const auto& h : ff.impropers_harmonic) {
      EXPECT_NEAR(h.k2, c.k_ct, 1e-8) << c.smi;
      EXPECT_NEAR(h.chi0, c.chi0, 1e-12) << c.smi;
    }
  }
  System r = build_molecule("CC1CC1CCC").system;   // a cyclopropane ring: TraPPE-UA has no bends for it
  std::string ch = "types";
  prepare_for_forcefield(r, def, ch);
  EXPECT_EQ(assign_types(r, def).untyped, 3);
}

// TraPPE-UA's alcohol, ether and ring terms (bench/ff/make_trappe.py, from Chen, Potoff & Siepmann 2001, Stubbs, Potoff &
// Siepmann 2004 and Keasler et al. 2012): poly(vinyl alcohol), poly(vinyl methyl ether) and poly(vinyl cyclohexane)
// oligomers parameterised whole and neutral with the tables' numbers; the ring torsions as a series in psi = phi - 180;
// a bridged ring (norbornane) stays untyped
TEST(FFDef, TrappeUaOxygenAndRingTerms) {
  const FFDef def = load_forcefield(std::string(CAPS_SOURCE_DIR) + "/data/forcefields/trappe-ua-umn.json");
  constexpr double KB = 0.0019872043;
  const double kPi = std::acos(-1.0);
  struct Want { const char* t[4]; int n; double v, delta; };
  struct Case { const char* smi; std::vector<Want> tors; const char* bend[3]; double th0, kt; };
  const std::vector<Case> cases = {
      // PVA: CHx-(CH)-CHy 112 deg at the hydroxyl CH (2001 Table 1); CHx-(CH)-(O)-H c1 197.33 K
      {"CC(O)CC(O)CC(O)C", {{{"CH2_3", "CH_10", "O_6", "H_7"}, 1, 197.33 * KB, 0}}, {"CH2_3", "CH_10", "CH2_3"}, 112.0, 62500.0 / 2 * KB},
      // PVME: CHx-CHy-O-CHz c1 725.35 K; CHx-(CHy)-O 112 deg, 50300 K (2004 Table 1)
      {"CC(OC)CC(OC)CC(OC)C", {{{"CH2_3", "CH_18", "O_20", "CH3_16"}, 1, 725.35 * KB, 0}}, {"CH2_3", "CH_18", "O_20"}, 112.0, 50300.0 / 2 * KB},
      // PVCH: the six-membered ring's C-C-C-C, c1 6840 K in psi (2012 Table 3); its C-C-C 114 deg (Table 2)
      {"CC(C1CCCCC1)CC(C1CCCCC1)C", {{{"CH2_140", "CH2_140", "CH2_140", "CH2_140"}, 1, 6840 * KB, 180},
                                      {{"CH2_140", "CH2_140", "CH2_140", "CH2_140"}, 3, 63 * KB, 180}},
       {"CH2_140", "CH_4", "CH2_140"}, 114.0, 62500.0 / 2 * KB},
  };
  for (const auto& c : cases) {
    System s = build_molecule(c.smi).system;
    std::string ch = "types";
    prepare_for_forcefield(s, def, ch);
    const TypingResult tr = assign_types(s, def);
    ASSERT_EQ(tr.untyped, 0) << c.smi;
    ParamReport rep;
    const ForceField ff = parameterize(s, def, tr.types, ch, &rep, false);
    EXPECT_TRUE(rep.missing.empty()) << c.smi;
    double q = 0;
    for (double x : ff.charge) q += x;
    EXPECT_NEAR(q, 0.0, 1e-9) << c.smi;
    for (const auto& w : c.tors) {
      int found = 0;
      for (const auto& d : ff.dihedrals) {
        const std::array<std::string, 4> t{tr.types[d.i], tr.types[d.j], tr.types[d.k], tr.types[d.l]};
        const bool fw = t[0] == w.t[0] && t[1] == w.t[1] && t[2] == w.t[2] && t[3] == w.t[3];
        const bool bw = t[3] == w.t[0] && t[2] == w.t[1] && t[1] == w.t[2] && t[0] == w.t[3];
        if ((fw || bw) && d.n == w.n) {
          EXPECT_NEAR(d.v, w.v, 1e-8) << c.smi << " " << w.t[1] << "-" << w.t[2];
          EXPECT_NEAR(std::fmod(d.delta * 180 / kPi + 720, 360.0), w.delta, 1e-6) << c.smi;
          ++found;
        }
      }
      EXPECT_GT(found, 0) << c.smi << " " << w.t[0] << "-" << w.t[1] << "-" << w.t[2] << "-" << w.t[3];
    }
    int nb = 0;
    for (const auto& a : ff.angles) {
      const std::string t0 = tr.types[a.i], t1 = tr.types[a.j], t2 = tr.types[a.k];
      if (t1 == c.bend[1] && ((t0 == c.bend[0] && t2 == c.bend[2]) || (t0 == c.bend[2] && t2 == c.bend[0]))) {
        EXPECT_NEAR(a.theta0 * 180 / kPi, c.th0, 1e-9) << c.smi;
        EXPECT_NEAR(a.kt, c.kt, 1e-9) << c.smi;
        ++nb;
      }
    }
    EXPECT_GT(nb, 0) << c.smi;
  }
  System r = build_molecule("C1CC2CCC1C2").system;   // norbornane: no TraPPE-UA paper fits a bridged ring
  std::string ch = "types";
  prepare_for_forcefield(r, def, ch);
  EXPECT_EQ(assign_types(r, def).untyped, 3);
}

// Inorganic sets from their papers (bench/ff/make_catlow1976_mgo.py, make_catlow1985_rutile.py, make_khan1998_zro2.py,
// make_sayle2000.py, make_sepliarsky2005.py): each crystal typed whole, every pair and spring found, neutral, and the
// tables' numbers in place (eV to kcal/mol); the quartic core-shell springs as class II K4 = c4/24
TEST(FFDef, InorganicSetsFromTheirPapers) {
  constexpr double EV = 23.060548;
  const std::string dir = std::string(CAPS_SOURCE_DIR) + "/data/forcefields/";
  struct Case { const char* ff; const char* group; double a, c; std::vector<CrystalSite> sites; };
  const std::vector<Case> cases = {
      {"catlow1976-mgo", "F m -3 m", 4.212, 4.212, {{"Mg1", 12, {0, 0, 0}}, {"O1", 8, {0.5, 0.5, 0.5}}}},
      {"khan1998-zro2", "F m -3 m", 5.07, 5.07, {{"Zr1", 40, {0, 0, 0}}, {"O1", 8, {0.25, 0.25, 0.25}}}},
      {"sayle2000-alkaline-earth-oxides", "F m -3 m", 5.52, 5.52, {{"Ba1", 56, {0, 0, 0}}, {"O1", 8, {0.5, 0.5, 0.5}}}},
      {"catlow1985-rutile", "P 42/m n m", 4.594, 2.959, {{"Ti1", 22, {0, 0, 0}}, {"O1", 8, {0.3053, 0.3053, 0}}}},
      {"sepliarsky2005-perovskites", "P m -3 m", 4.0, 4.0, {{"K1", 19, {0, 0, 0}}, {"Nb1", 41, {0.5, 0.5, 0.5}}, {"O1", 8, {0.5, 0.5, 0}}}},
      {"sepliarsky2005-perovskites", "P m -3 m", 3.905, 3.905, {{"Sr1", 38, {0, 0, 0}}, {"Ti1", 22, {0.5, 0.5, 0.5}}, {"O1", 8, {0.5, 0.5, 0}}}},
  };
  for (const auto& c : cases) {
    const FFDef def = load_forcefield(dir + c.ff + ".json");
    CrystalSpec spec;
    spec.space_group = c.group;
    spec.a = spec.b = c.a;
    spec.c = c.c;
    spec.sites = c.sites;
    spec.supercell = {2, 2, 2};
    System s = build_crystal(spec);
    s.bonds = perceive_bonds(s);
    std::string ch = "types";
    prepare_for_forcefield(s, def, ch);
    const TypingResult tr = assign_types(s, def);
    ASSERT_EQ(tr.untyped, 0) << c.ff << " " << c.sites[0].label;
    ParamReport rep;
    const ForceField ff = parameterize(s, def, tr.types, ch, &rep, false);
    EXPECT_TRUE(rep.missing.empty()) << c.ff << " " << c.sites[0].label;
    double q = 0;
    for (double x : ff.charge) q += x;
    EXPECT_NEAR(q, 0.0, 1e-9) << c.ff << " " << c.sites[0].label;
  }
  auto pair = [&](const FFDef& d, const std::string& name) {
    const auto p = std::find_if(d.pairs.begin(), d.pairs.end(), [&](const FFRule& x) { return x.name == name; });
    EXPECT_NE(p, d.pairs.end()) << name;
    return p == d.pairs.end() ? std::vector<double>(4, 0.0) : p->params;
  };
  auto bond = [&](const FFDef& d, const std::string& name) {
    const auto p = std::find_if(d.bonds.begin(), d.bonds.end(), [&](const FFRule& x) { return x.name == name; });
    EXPECT_NE(p, d.bonds.end()) << name;
    return p == d.bonds.end() ? std::vector<double>(4, 0.0) : p->params;
  };
  const FFDef mgo = load_forcefield(dir + "catlow1976-mgo.json");   // Table 1, model 1
  EXPECT_NEAR(pair(mgo, "Mg_c76m1-O_c76m1s")[0], 1152.0 * EV, 1e-9);
  EXPECT_NEAR(pair(mgo, "O_c76m1s-O_c76m1s")[2], 28.96 * EV, 1e-9);
  EXPECT_NEAR(bond(mgo, "O_c76m1c-O_c76m1s")[1], 38.03 / 2 * EV, 1e-9);
  EXPECT_DOUBLE_EQ(mgo.type("O_c76m1s")->charge, -2.62);
  const FFDef rut = load_forcefield(dir + "catlow1985-rutile.json");   // Table IVa: K 37.3 / 86.4, quartic 1.1e5
  EXPECT_NEAR(bond(rut, "Ti_cfrc-Ti_cfrs")[1], 37.3 / 2 * EV, 1e-9);
  EXPECT_NEAR(bond(rut, "O_cfrc-O_cfrs")[3], 1.1e5 / 24 * EV, 1e-6);
  EXPECT_NEAR(rut.type("Ti_cfrc")->charge + rut.type("Ti_cfrs")->charge, 4.0, 1e-12);
  const FFDef zr = load_forcefield(dir + "khan1998-zro2.json");
  EXPECT_NEAR(pair(zr, "Zr_k98s-O_k98s")[0], 985.869 * EV, 1e-9);
  EXPECT_DOUBLE_EQ(zr.type("Zr_k98s")->charge, 1.35);
  const FFDef sw = load_forcefield(dir + "sayle2000-alkaline-earth-oxides.json");
  EXPECT_NEAR(pair(sw, "Mg_sw-O_sw")[1], 0.2945, 1e-12);
  EXPECT_NEAR(pair(sw, "O_sw-O_sw")[0], 22764.0 * EV, 1e-9);
  const FFDef sp = load_forcefield(dir + "sepliarsky2005-perovskites.json");
  EXPECT_NEAR(pair(sp, "K_knts-O_knts")[0], 124872.44 * EV, 1e-6);
  EXPECT_NEAR(bond(sp, "O_bstc-O_bsts")[3], 4000.0 / 24 * EV, 1e-9);
  EXPECT_EQ(sp.type("Pb_ptc"), nullptr);   // PbTiO3 withheld: Table 1's Pb-Ti row makes the crystal unbound (make_sepliarsky2005.py)
  EXPECT_DOUBLE_EQ(sp.type("Nb_kntc")->charge, -2.9878);
}

// Bandura & Kubicki's hydroxylated rutile (110) (tests/data/rutile110_hydroxylated.pdb, bench/ff/make_bandura2003.py
// --test-slab): Table 9's surface sites from the neighbours, Ti-O bonds dropped, the O-H Morse bonds kept, Ti-O-H bends
// by distance, Table 10's surface pairs
TEST(FFDef, BanduraKubickiSurfaceSites) {
  const FFDef def = load_forcefield(std::string(CAPS_SOURCE_DIR) + "/data/forcefields/bandura2003-tio2-water.json");
  Trajectory t = read_pdb(std::string(CAPS_SOURCE_DIR) + "/tests/data/rutile110_hydroxylated.pdb");
  System s = t.topology;
  s.bonds = perceive_bonds(s);
  std::string ch = "types";
  prepare_for_forcefield(s, def, ch);
  const TypingResult tr = assign_types(s, def);
  ASSERT_EQ(tr.untyped, 0);
  std::map<std::string, int> n;
  for (const auto& x : tr.types) n[x]++;
  // two surfaces of 2 x 3 cells: 6 TiV and 6 TiVI each; on top one TiV holds an OH, one a water; one bridging O protonated
  EXPECT_EQ(n["TiV_bk"], 11);
  EXPECT_EQ(n["TiVIoh_bk"], 1);
  EXPECT_EQ(n["TiVI_bk"], 12);
  EXPECT_EQ(n["Ti_bk"], 12);
  EXPECT_EQ(n["Os_bk"], 24);   // the in-plane O of the two surface Ti layers
  EXPECT_EQ(n["Ob_bk"], 11);
  EXPECT_EQ(n["Obh_bk"], 1);
  EXPECT_EQ(n["Ooh_bk"], 1);
  EXPECT_EQ(n["Ow_bk"], 1);
  EXPECT_EQ(n["O_bk"], 36);
  ParamReport rep;
  const ForceField ff = parameterize(s, def, tr.types, ch, &rep, false);
  EXPECT_TRUE(rep.missing.empty());
  EXPECT_EQ(ff.bonds_x.size(), 4u);   // Morse: the two hydroxyls' O-H and water's two
  ASSERT_EQ(ff.angles.size(), 4u);    // water's H-O-H; Ti-Obh-H twice, Ti-Ooh-H once
  for (const auto& a : ff.angles) {
    const bool water = tr.types[a.j] == "Ow_bk";
    EXPECT_NEAR(a.kt, water ? 103.045 / 2 : 14.136 / 2, 1e-12);
    if (!water) EXPECT_EQ(s.atoms[a.i].element == 22 ? s.atoms[a.k].element : s.atoms[a.i].element, 1);
  }
  const FFRule* os = nullptr;
  for (const auto& r : def.pairs)
    if (r.name == "TiV_bk-Os_bk") os = &r;
  ASSERT_NE(os, nullptr);
  EXPECT_DOUBLE_EQ(os->params[0], 449223.3);
}

// CL&P's deep-eutectic-solvent set from Padua's des.ff (bench/ff/convert_fftool.py, bench/typing/make_clandp_des_rules.py):
// each DES component typed whole, at its integer charge with des.ff's type charges, every bonded term found; the
// components des.ff cannot make neutral stay untyped; des.ff's kJ/mol numbers in place (k/2 constants, V2/2 impropers);
// a four-connected P on carbons is perceived as phosphonium (+1)
TEST(FFDef, ClandpDesComponentsTypedWhole) {
  const FFDef def = load_forcefield(std::string(CAPS_SOURCE_DIR) + "/data/forcefields/clandp-des.json");
  auto typed = [&](const std::string& smi, std::map<std::string, int>* count = nullptr) {
    System s = build_molecule(smi).system;
    std::string ch = "types";
    prepare_for_forcefield(s, def, ch);
    const TypingResult tr = assign_types(s, def);
    if (count)
      for (const auto& t : tr.types) ++(*count)[t];
    if (tr.untyped) return std::make_pair(tr.untyped, 0.0);
    ParamReport rep;
    const ForceField f = parameterize(s, def, tr.types, "types", &rep, false);
    EXPECT_TRUE(rep.missing.empty()) << smi << ": " << (rep.missing.empty() ? "" : rep.missing[0]);
    double q = 0;
    for (double x : f.charge) q += x;
    int formal = 0;
    for (int c : perceive(s).charge) formal += c;
    EXPECT_NEAR(q, formal, 1e-9) << smi;
    return std::make_pair(0, q);
  };
  const std::vector<std::pair<std::string, int>> whole = {
      {"C[N+](C)(C)CCO", 1}, {"CCCC[N+](CCCC)(CCCC)CCCC", 1}, {"CC[NH+](CC)CCO", 1}, {"CCCC[P+](CCCC)(CCCC)CCCC", 1},
      {"C[N+](C)(C)CC(=O)[O-]", 0}, {"[Cl-]", -1}, {"[Br-]", -1}, {"CCC(=O)[O-]", -1}, {"NC(N)=O", 0}, {"CCO", 0},
      {"OCCO", 0}, {"OCC(O)CO", 0}, {"CC(=O)O", 0}, {"OC(=O)CCC(=O)O", 0}, {"CC(=O)CCC(=O)O", 0}, {"CC(C)=O", 0}};
  for (const auto& [smi, q] : whole) {
    const auto r = typed(smi);
    EXPECT_EQ(r.first, 0) << smi;
    EXPECT_NEAR(r.second, q, 1e-9) << smi;
  }
  for (const std::string smi : {"CO", "OC(=O)CC(=O)O", "OC(=O)C(=O)O", "CC(O)C(=O)O", "CC[NH3+]"})
    EXPECT_GT(typed(smi).first, 0) << smi;
  std::map<std::string, int> c;   // cholinium as Padua's cholinium.xyz labels it
  typed("C[N+](C)(C)CCO", &c);
  EXPECT_EQ(c, (std::map<std::string, int>{{"N4", 1}, {"C1", 3}, {"C1A", 1}, {"COL", 1}, {"H1", 11}, {"HC", 2}, {"OH", 1}, {"HO", 1}}));
  c.clear();
  typed("OCC(O)CO", &c);
  EXPECT_EQ(c, (std::map<std::string, int>{{"CTO", 2}, {"CSO", 1}, {"H1O", 5}, {"OHT", 3}, {"HOT", 3}}));
  c.clear();
  typed("C[N+](C)(C)CC(=O)[O-]", &c);
  EXPECT_EQ(c["N4B"], 1);
  EXPECT_EQ(c["C2C"], 1);
  auto rule = [&](const std::vector<FFRule>& v, const std::string& name) {
    const auto p = std::find_if(v.begin(), v.end(), [&](const FFRule& x) { return x.name == name; });
    EXPECT_NE(p, v.end()) << name;
    return p == v.end() ? std::vector<double>(4, 0.0) : p->params;
  };
  EXPECT_NEAR(rule(def.dihedrals, "OH-CT-CT-NT")[0], -44.0515 / 4.184, 1e-12);
  EXPECT_NEAR(rule(def.bonds, "C-N")[0], 4100.3 / 2 / 4.184, 1e-12);
  EXPECT_NEAR(rule(def.impropers, "N-N-C-O")[0], 87.864 / 2 / 4.184, 1e-12);
  EXPECT_DOUBLE_EQ(def.type("CU")->charge, 0.142);
}

// OPLS-UA: the 2024 file's united-atom types, with the hydrocarbon terms it lacks from Jorgensen, Madura & Swenson 1984
// (bench/ff/make_oplsua2024.py, bench/ff/open/oplsua2024_supplement.json). Polypropylene, 1,4-polybutadiene and
// cis-1,4-polyisoprene oligomers are complete and neutral; the 2024 file's own terms are kept where it has them
TEST(FFDef, OplsUa1984HydrocarbonTerms) {
  const FFDef def = load_forcefield(std::string(CAPS_SOURCE_DIR) + "/data/forcefields/oplsua2024.json");
  auto apply = [&](const std::string& smi) {
    System s = build_molecule(smi).system;
    std::string ch = "types";
    prepare_for_forcefield(s, def, ch);
    const TypingResult tr = assign_types(s, def);
    EXPECT_EQ(tr.untyped, 0) << smi;
    ParamReport rep;
    const ForceField f = parameterize(s, def, tr.types, ch, &rep, false);
    EXPECT_TRUE(rep.missing.empty()) << smi << ": " << (rep.missing.empty() ? "" : rep.missing[0]);
    double q = 0;
    for (double x : f.charge) q += x;
    EXPECT_NEAR(q, 0.0, 1e-12) << smi;
    return rep.used;
  };
  const auto pp = apply("CC(C)CC(C)CC(C)CC(C)C");
  EXPECT_GT(pp.count("dihedral C3-CH-C2-CH"), 0u);   // Table II isopentane, per torsion
  EXPECT_GT(pp.count("dihedral C2-CH-C2-CH"), 0u);
  const auto pb = apply("C/C=C\\CC/C=C/CC/C=C\\C");
  EXPECT_GT(pb.count("dihedral C2-C2-C8-C8"), 0u);   // Table II 1-butene
  EXPECT_GT(pb.count("dihedral C2-C8-C8-C2"), 0u);   // the 2024 file's alkene V2
  EXPECT_GT(pb.count("dihedral C8-C2-C2-C8"), 0u);   // Table II other n-alkanes
  EXPECT_GT(pb.count("angle C2-C2-C8"), 0u);         // Table I C-C4-C
  const auto pi = apply("C/C=C(C)\\CC/C=C(C)\\CC/C=C(C)\\C");
  EXPECT_GT(pi.count("angle C2-C8-C7"), 0u);         // Table I C-C3-C
  EXPECT_GT(pi.count("dihedral C2-C2-C7-C3"), 0u);   // zero: the C= carbon's single-bonded substituent
  apply("CCCCCCCC");
  auto rule = [&](const std::vector<FFRule>& v, const std::string& name) {
    const auto p = std::find_if(v.begin(), v.end(), [&](const FFRule& x) { return x.name == name; });
    EXPECT_NE(p, v.end()) << name;
    return p == v.end() ? std::vector<double>(4, 0.0) : p->params;
  };
  EXPECT_EQ(rule(def.dihedrals, "C3-CH-C2-CH"), (std::vector<double>{1.526, -0.533, 1.7265, 0.0}));
  EXPECT_EQ(rule(def.dihedrals, "C2-C2-C8-C8"), (std::vector<double>{0.343, -0.436, -1.121, 0.0}));
  EXPECT_EQ(rule(def.dihedrals, "C8-C2-C2-C8"), (std::vector<double>{1.411, -0.271, 3.145, 0.0}));
  EXPECT_EQ(rule(def.angles, "C2-C2-C8"), (std::vector<double>{63.0, 112.0}));
  EXPECT_EQ(rule(def.dihedrals, "C2-C2-C2-C2"), (std::vector<double>{-3.4, 1.25, 3.1, 0.0}));   // the 2024 file's, kept
  EXPECT_EQ(std::count_if(def.dihedrals.begin(), def.dihedrals.end(), [](const FFRule& r) { return r.name == "C2-C2-C2-C2"; }), 1);
}

// OPLS-UA alcohols and ethers: the terms the 2024 file lacks from Jorgensen 1986 (Tables II, III) and Briggs, Matsui &
// Jorgensen 1990 (Table II, p. 960). Poly(vinyl alcohol), polyoxytetramethylene and the papers' own alcohols and diethyl
// ether are complete and neutral; poly(ethylene oxide)'s O-C-C-O torsion is in neither paper and stays missing
TEST(FFDef, OplsUaAlcoholEtherTerms) {
  const FFDef def = load_forcefield(std::string(CAPS_SOURCE_DIR) + "/data/forcefields/oplsua2024.json");
  auto apply = [&](const std::string& smi, ParamReport& rep) {
    System s = build_molecule(smi).system;
    std::string ch = "types";
    prepare_for_forcefield(s, def, ch);
    const TypingResult tr = assign_types(s, def);
    EXPECT_EQ(tr.untyped, 0) << smi;
    const ForceField f = parameterize(s, def, tr.types, ch, &rep, false);
    double q = 0;
    for (double x : f.charge) q += x;
    EXPECT_NEAR(q, 0.0, 1e-9) << smi;
  };
  for (const std::string smi : {"CC(O)CC(O)CC(O)CC(O)C", "CCCCOCCCCOCCCCOCCCC", "CCO.CCCO.CC(C)O.CC(C)(C)O.CCOCC"}) {
    ParamReport rep;
    apply(smi, rep);
    EXPECT_TRUE(rep.missing.empty()) << smi << ": " << (rep.missing.empty() ? "" : rep.missing[0]);
  }
  // PEO: none of the papers has O-C-C-O; CAPS's fit to a 1,2-dimethoxyethane scan (RMS 0.04 kcal/mol) completes it
  {
    ParamReport peo;
    EXPECT_NO_THROW(apply("CCOCCOCCOCC", peo));
    EXPECT_TRUE(peo.missing.empty());
  }
  auto rule = [&](const std::vector<FFRule>& v, const std::string& name) {
    const auto p = std::find_if(v.begin(), v.end(), [&](const FFRule& x) { return x.name == name; });
    EXPECT_NE(p, v.end()) << name;
    return p == v.end() ? std::vector<double>(4, 0.0) : p->params;
  };
  EXPECT_EQ(rule(def.dihedrals, "C2-C2-OS-C2"), (std::vector<double>{2.8827, -0.6507, 2.2183, 0.0}));   // 1990, EME / DEE
  EXPECT_EQ(rule(def.dihedrals, "CH-CH-OH-HO"), (std::vector<double>{0.784, -0.125, 0.3455, 0.0}));     // 1986 2-propanol
  EXPECT_EQ(rule(def.dihedrals, "CH-C2-CH-OH"), (std::vector<double>{0.702, -0.212, 3.06, 0.0}));       // 1986 1-propanol
  EXPECT_EQ(rule(def.angles, "C2-OS-CH"), (std::vector<double>{100.0, 112.0}));                          // 1990 Table II
  EXPECT_EQ(rule(def.dihedrals, "C2-C2-OH-HO"), (std::vector<double>{0.3, 0.0, 1.3, 0.0}));              // the 2024 file's, kept
  EXPECT_EQ(std::count_if(def.dihedrals.begin(), def.dihedrals.end(), [](const FFRule& r) { return r.name == "C2-C2-OH-HO"; }), 1);
  EXPECT_EQ(rule(def.dihedrals, "OS-C2-C2-OS"), (std::vector<double>{2.0063, -1.8941, 4.0419, -0.1512}));  // CAPS fit (QM), DME
}

// OPLS-UA sulfur (Jorgensen, J. Phys. Chem. 90, 6379 (1986)): thiols, sulfides and disulfides typed, neutral and complete,
// the C-S and S-S rotations as Table III prints them
TEST(FFDef, OplsUaSulfurTerms) {
  const FFDef def = load_forcefield(std::string(CAPS_SOURCE_DIR) + "/data/forcefields/oplsua2024.json");
  for (const std::string smi : {"CCS", "CCSC", "CCSCC", "CSSC", "CCSSCC", "CCCSCCC"}) {
    System s = build_molecule(smi).system;
    std::string ch = "types";
    prepare_for_forcefield(s, def, ch);
    const TypingResult tr = assign_types(s, def);
    EXPECT_EQ(tr.untyped, 0) << smi;
    ParamReport rep;
    const ForceField f = parameterize(s, def, tr.types, ch, &rep, false);
    double q = 0;
    for (double x : f.charge) q += x;
    EXPECT_NEAR(q, 0.0, 1e-9) << smi;
    EXPECT_TRUE(rep.missing.empty()) << smi << ": " << (rep.missing.empty() ? "" : rep.missing[0]);
  }
  auto rule = [&](const std::string& name) {
    const auto p = std::find_if(def.dihedrals.begin(), def.dihedrals.end(), [&](const FFRule& x) { return x.name == name; });
    EXPECT_NE(p, def.dihedrals.end()) << name;
    return p == def.dihedrals.end() ? std::vector<double>(4, 0.0) : p->params;
  };
  EXPECT_EQ(rule("C3-C2-SH-HS"), (std::vector<double>{-0.2067, -0.0733, 1.2167, 0.0}));   // ethanethiol
  EXPECT_EQ(rule("C2-S-C2-C3"), (std::vector<double>{1.46, -1.073, 2.31, 0.0}));          // diethyl sulfide
  EXPECT_EQ(rule("C3-S-S-C3"), (std::vector<double>{1.582, -8.67, 2.008, 0.0}));          // dimethyl disulfide
}

// OPLS-DES from Doherty & Acevedo's paper and SI only (bench/ff/make_opls_des_doherty2018.py,
// bench/typing/make_oplsdes_doherty2018_rules.py): the paper's nine molecules typed as Schemes S1 and S2 type them, choline
// and chloride at +-0.8 e (ion_charge_scale), the HBDs neutral, every bonded term found; molecules outside the set
// untyped; numbers as the SI's tables print them
TEST(FFDef, OplsDesDoherty2018FromThePaper) {
  const FFDef def = load_forcefield(std::string(CAPS_SOURCE_DIR) + "/data/forcefields/oplsdes-doherty2018.json");
  EXPECT_DOUBLE_EQ(def.ion_charge_scale, 0.8);
  auto typed = [&](const std::string& smi, std::map<std::string, int>* count = nullptr) {
    System s = build_molecule(smi).system;
    std::string ch = "types";
    prepare_for_forcefield(s, def, ch);
    const TypingResult tr = assign_types(s, def);
    if (count)
      for (const auto& t : tr.types) ++(*count)[t];
    if (tr.untyped) return std::make_pair(tr.untyped, 0.0);
    ParamReport rep;
    const ForceField f = parameterize(s, def, tr.types, "types", &rep, false);
    EXPECT_TRUE(rep.missing.empty()) << smi << ": " << (rep.missing.empty() ? "" : rep.missing[0]);
    double q = 0;
    for (double x : f.charge) q += x;
    int formal = 0;
    for (int c : perceive(s).charge) formal += c;
    EXPECT_NEAR(q, 0.8 * formal, 1e-9) << smi;
    return std::make_pair(0, q);
  };
  const std::vector<std::pair<std::string, double>> whole = {
      {"C[N+](C)(C)CCO", 0.8}, {"[Cl-]", -0.8}, {"OCCO", 0}, {"OCC(O)CO", 0}, {"CC(=O)CCC(=O)O", 0}, {"OC(=O)CC(=O)O", 0},
      {"OC(=O)C(=O)O", 0}, {"Oc1ccccc1", 0}, {"NC(N)=O", 0}};
  for (const auto& [smi, q] : whole) {
    const auto r = typed(smi);
    EXPECT_EQ(r.first, 0) << smi;
    EXPECT_NEAR(r.second, q, 1e-9) << smi;
  }
  for (const std::string smi : {"CC(O)CO", "CC(=O)O", "OC(=O)CCC(=O)O", "C[N+](C)(C)C", "C[N+](C)(C)CC(=O)[O-]", "[Br-]", "CNC(N)=O"})
    EXPECT_GT(typed(smi).first, 0) << smi;
  std::map<std::string, int> c;   // Scheme S1
  typed("C[N+](C)(C)CCO", &c);
  EXPECT_EQ(c, (std::map<std::string, int>{{"NA", 1}, {"CA", 3}, {"HA", 9}, {"CS", 1}, {"HS", 2}, {"CW", 1}, {"HW", 2}, {"OY", 1}, {"HY", 1}}));
  c.clear();   // Scheme S2
  typed("OCC(O)CO", &c);
  EXPECT_EQ(c, (std::map<std::string, int>{{"CM_gly", 1}, {"OM_gly", 1}, {"HM_gly", 1}, {"HZ_gly", 1}, {"CB_gly", 2}, {"HC_gly", 4},
                                           {"OH_gly", 2}, {"HO_gly", 2}}));
  c.clear();
  typed("CC(=O)CCC(=O)O", &c);
  EXPECT_EQ(c["CZ_lev"], 1);
  EXPECT_EQ(c["CF_lev"], 1);
  EXPECT_EQ(c["OC_lev"], 1);
  EXPECT_EQ(c["OB_lev"], 1);
  auto rule = [&](const std::vector<FFRule>& v, const std::string& name) {
    const auto p = std::find_if(v.begin(), v.end(), [&](const FFRule& x) { return x.name == name; });
    EXPECT_NE(p, v.end()) << name;
    return p == v.end() ? std::vector<double>(4, 0.0) : p->params;
  };
  EXPECT_DOUBLE_EQ(def.type("NA")->charge, 0.791);                         // Table S1
  EXPECT_DOUBLE_EQ(def.type("OM_gly")->charge, -0.5447);                   // Table S4
  EXPECT_EQ(rule(def.dihedrals, "OY-CW-CS-NA"), (std::vector<double>{-6.0, -5.0, 3.2, 0.0}));   // Table S3
  EXPECT_EQ(rule(def.bonds, "CS-HS"), (std::vector<double>{340.0, 1.0805}));                  // Table S2
  EXPECT_EQ(rule(def.angles, "CD_ox-CD_ox-OH_ox"), (std::vector<double>{70.96, 118.03}));     // Table S5
  EXPECT_EQ(rule(def.angles, "OH_lev-CB_lev-CZ_lev"), (std::vector<double>{70.0, 108.0}));    // Table S5's HO-CB-CZ
  EXPECT_EQ(rule(def.impropers, "*-CD_mal-*-*"), (std::vector<double>{10.5, -1.0, 2.0}));     // Table S6, V2 21
  EXPECT_EQ(rule(def.impropers, "*-NA-*-*"), (std::vector<double>{1.0, -1.0, 2.0}));          // Table S3, V2 2
}
