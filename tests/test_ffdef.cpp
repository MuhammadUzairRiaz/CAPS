#include <gtest/gtest.h>

#include <algorithm>
#include <tuple>
#include <array>
#include <cmath>
#include <filesystem>
#include <map>
#include <random>
#include <set>

#include "caps/ffdef.hpp"
#include "caps/field.hpp"
#include "caps/grow.hpp"
#include "caps/io.hpp"
#include "caps/martini_protein.hpp"
#include "caps/molecule.hpp"
#include "caps/resolution.hpp"
#include "caps/typing.hpp"
#include "caps/uff.hpp"

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

// mW water: the Stillinger–Weber two- and three-body forces and the virial are the energy's derivatives (sites on a
// jittered lattice at liquid density, pairs and triplets across the periodic faces)
TEST(FieldForms, StillingerWeberForcesMatchFiniteDifferences) {
  FFDef def = load_forcefield(std::string(CAPS_SOURCE_DIR) + "/data/forcefields/mw-moltemplate.json");
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
  const FFDef martini = load_forcefield(dir + "martini-moltemplate.json");
  check_cg_forces(martini, cg_box(martini, {{"DPPC", 3}, {"NA+", 2}, {"CL-", 2}, {"W", 30}}, 36.0, 3), "keep");
  const FFDef sdk = load_forcefield(dir + "sdk-moltemplate.json");
  check_cg_forces(sdk, cg_box(sdk, {{"[NC][PH][GL]([EST1][CM][CM][CT2])[EST2][CM][CM][CT2]", 3}, {"[W]", 30}}, 32.0, 5), "types");
  const FFDef aa = load_forcefield(dir + "martini-aminoacids.json");   // harmonic impropers, GROMACS type-2 order
  check_cg_forces(aa, cg_box(aa, {{"TRP", 2}, {"HIS", 2}, {"PHE", 2}, {"TYR", 2}, {"LYS", 2}, {"ASP", 2}, {"[P4]", 30}}, 36.0, 9), "keep");
  const FFDef cd = load_forcefield(dir + "cooke-deserno-moltemplate.json");
  check_cg_forces(cd, cg_box(cd, {{"lipid", 40}}, 12.0, 7, 1.2), "types");
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
  const FFDef sdk = load_forcefield(std::string(CAPS_SOURCE_DIR) + "/data/forcefields/sdk-moltemplate.json");
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
