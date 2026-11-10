#include <gtest/gtest.h>

#include <cmath>

#include "caps/crystal.hpp"
#include "caps/dielectric.hpp"
#include "caps/elements.hpp"
#include "caps/phonons.hpp"
#include "caps/qha.hpp"
#include "caps/uff.hpp"

using namespace caps;

namespace {
constexpr double kPi = 3.14159265358979323846;
// a force field of ions only (no bonded terms, no exclusions), types by element
ForceField ions(const System& s, const std::vector<std::pair<int, double>>& charge) {
  ForceField ff = default_forcefield(s);
  ff.type_names.clear();
  for (const auto& [z, q] : charge) ff.type_names.push_back(element(z).symbol);
  ff.lj.assign(charge.size(), {0, 0});
  for (size_t i = 0; i < s.atoms.size(); ++i)
    for (size_t t = 0; t < charge.size(); ++t)
      if (s.atoms[i].element == charge[t].first) ff.type_index[i] = int(t), ff.atom_type[i] = ff.type_names[t], ff.charge[i] = charge[t].second;
  ff.bonds.clear(), ff.angles.clear(), ff.dihedrals.clear(), ff.impropers.clear(), ff.inversions.clear(), ff.bonds_x.clear(), ff.angles_x.clear();
  ff.pairs14.clear();
  for (auto& e : ff.excluded) e.clear();
  return ff;
}
ForceField salt_ff(const System& s) {
  ForceField ff = ions(s, {{11, 1.0}, {17, -1.0}});
  ff.pair_func[{0, 1}] = {1, 26990.0, 0.3216, 0.0};
  ff.pair_func[{1, 1}] = {1, 50000.0, 0.33, 1500.0};
  return ff;
}
System salt(double a) {
  System s;
  s.cell.a = {0, a / 2, a / 2}, s.cell.b = {a / 2, 0, a / 2}, s.cell.c = {a / 2, a / 2, 0};
  for (auto [z, f] : std::vector<std::pair<int, Vec3>>{{11, {0, 0, 0}}, {17, {0.5, 0.5, 0.5}}}) {
    Atom at;
    at.element = z;
    at.pos = s.cell.to_cartesian(f);
    s.atoms.push_back(at);
  }
  return s;
}
EnergyOptions ewald() {
  EnergyOptions e;
  e.cutoff = 8.0, e.tail = false, e.electrostatics = EnergyOptions::Electrostatics::PME, e.ewald_rtol = 1e-7, e.pme_spacing = 0.5;
  return e;
}
}  // namespace

// Rigid-ion rock salt: ε⁰ − 1 = q²/(ε₀ V μ ω_TO²) with ω_TO the transverse optical frequency the phonon module finds at Γ;
// ε(ω → 0) = ε⁰; ε^∞ = 1; no piezoelectricity (centrosymmetric)
TEST(Dielectric, RockSaltIsOneOscillator) {
  DielectricOptions o;
  o.energy = ewald();
  o.field_for = salt_ff;
  o.w_max = 400, o.points = 400, o.damping = 1;
  const DielectricResult r = crystal_dielectric(salt(5.65), o);
  ASSERT_EQ(r.modes.size(), 3u);
  PhononOptions po;
  po.field_for = salt_ff, po.energy = ewald(), po.relax = false, po.step = 0.002;
  po.path = {{0, 0, 0}, {0.001, 0, 0}}, po.labels = {"G", "x"}, po.path_points = 2, po.mesh = {1, 1, 1}, po.temperatures = {300};
  System cell = salt(5.65);
  cell.cell = r.cell_vectors;
  for (size_t i = 0; i < 2; ++i) cell.atoms[i].pos = r.positions[i];
  const PhononResult ph = phonons(cell, po);
  const double wto = ph.bands[0].back();   // cm⁻¹
  const double mu = 22.98977 * 35.453 / (22.98977 + 35.453) * 1e-3 / 6.02214076e23;
  const double w = 2 * kPi * 2.99792458e10 * wto;
  const double expect = 1 + std::pow(1.602176634e-19, 2) / (8.8541878128e-12 * r.volume * 1e-30 * mu * w * w);
  std::printf("ε⁰ %.6f (one oscillator at the phonon module's ω_TO %.3f cm⁻¹: %.6f) · modes %.3f %.3f %.3f cm⁻¹ · ε(0) %.6f\n", r.eps0_principal[0], wto, expect,
              r.modes[0], r.modes[1], r.modes[2], r.eps_re[0][0]);
  for (int k = 0; k < 3; ++k) EXPECT_NEAR(r.eps0_principal[size_t(k)], expect, 2e-4 * expect);
  for (double m : r.modes) EXPECT_NEAR(m, wto, 2e-3 * wto);
  EXPECT_NEAR(r.eps_re[0][0], r.eps0_principal[0], 1e-7 * r.eps0_principal[0]);
  EXPECT_NEAR(r.epsinf_principal[0], 1.0, 1e-12);
  double emax = 0;
  for (const auto& row : r.e)
    for (double v : row) emax = std::max(emax, std::fabs(v));
  EXPECT_LT(emax, 1e-6);
}

// A shell on the anion: ε^∞ > 1 and ε(ω) far above the modes reaches it; ε(0) = ε⁰
TEST(Dielectric, ShellModelHighFrequencyLimit) {
  auto shelled = [](double a) {
    System s = salt(a);
    Atom sh = s.atoms[1];
    s.atoms.push_back(sh);   // the Cl shell at its core
    return s;
  };
  auto ff_for = [](const System& s) {
    ForceField ff = salt_ff(s);
    // the shells: every third atom of each copy (Na, Cl core, Cl shell); type 2 "Cls", charge −2.5, core +1.5
    ff.type_names.push_back("Cls");
    ff.lj.push_back({0, 0});
    for (size_t i = 2; i < s.atoms.size(); i += 3) {
      ff.type_index[i] = 2, ff.atom_type[i] = "Cls", ff.charge[i] = -2.5, ff.charge[i - 1] = 1.5;
      ff.bonds.push_back({uint32_t(i - 1), uint32_t(i), 1500.0, 0.0});
      ff.excluded[i - 1].push_back(uint32_t(i)), ff.excluded[i].push_back(uint32_t(i - 1));
    }
    // the short-range pairs act on the shell
    ff.pair_func.clear();
    ff.pair_func[{0, 2}] = {1, 26990.0, 0.3216, 0.0};
    ff.pair_func[{2, 2}] = {1, 50000.0, 0.33, 1500.0};
    return ff;
  };
  DielectricOptions o;
  o.energy = ewald();
  o.field_for = ff_for;
  o.shell = {0, 0, 1};
  o.w_max = 20000, o.points = 200, o.damping = 1, o.piezo = false;
  const DielectricResult r = crystal_dielectric(shelled(5.65), o);
  std::printf("ε^∞ %.6f · ε⁰ %.6f · ε(%.0f cm⁻¹) %.6f · ε(0) %.6f · n %.4f · max force %.2e\n", r.epsinf_principal[0], r.eps0_principal[0], r.w.back(), r.eps_re[0].back(),
              r.eps_re[0][0], r.refractive[0], r.max_force);
  EXPECT_GT(r.epsinf_principal[0], 1.05);
  EXPECT_GT(r.eps0_principal[0], r.epsinf_principal[0]);
  EXPECT_NEAR(r.eps_re[0].back(), r.epsinf_principal[0], 2e-3 * r.epsinf_principal[0]);
  EXPECT_NEAR(r.eps_re[0][0], r.eps0_principal[0], 1e-6 * r.eps0_principal[0]);
}

// Zinc blende (rigid ±2 ions): the proper piezoelectric e14 from the internal-strain forces equals the polarisation of the
// ions relaxed under a small shear (beyond the affine motion), and the other components vanish by symmetry
TEST(Dielectric, ZincBlendePiezoelectricByDirectStrain) {
  auto zns_ff = [](const System& s) {
    ForceField ff = ions(s, {{30, 2.0}, {16, -2.0}});
    ff.pair_func[{0, 1}] = {1, 25000.0, 0.35, 0.0};
    ff.pair_func[{1, 1}] = {1, 50000.0, 0.33, 1500.0};
    return ff;
  };
  DielectricOptions o;
  o.energy = ewald();
  o.field_for = zns_ff;
  o.points = 10;
  const System zns = read_cif(std::string(CAPS_SOURCE_DIR) + "/data/crystals/sphalerite.cif");
  const DielectricResult r = crystal_dielectric(zns, o);
  std::printf("modes %.1f … %.1f cm⁻¹, max force %.1e; e14 %.6f e25 %.6f e36 %.6f C/m², e11 %.1e\n", r.modes.front(), r.modes.back(), r.max_force, r.e[0][3], r.e[1][4], r.e[2][5],
              r.e[0][0]);
  ASSERT_GT(r.modes.front(), 0) << "zinc blende is not a minimum of this model";
  EXPECT_GT(std::fabs(r.e[0][3]), 1e-3);
  EXPECT_NEAR(r.e[1][4], r.e[0][3], 1e-6 + 1e-4 * std::fabs(r.e[0][3]));
  EXPECT_NEAR(r.e[2][5], r.e[0][3], 1e-6 + 1e-4 * std::fabs(r.e[0][3]));
  EXPECT_LT(std::fabs(r.e[0][0]), 1e-6);
  // direct: shear γ_yz, the ions relaxed in the strained cell
  const double g = 2e-3;
  System st = zns;
  st.cell = r.cell_vectors;
  auto F = [&](const Vec3& v) { return Vec3{v[0], v[1] + g / 2 * v[2], v[2] + g / 2 * v[1]}; };
  st.cell.a = F(st.cell.a), st.cell.b = F(st.cell.b), st.cell.c = F(st.cell.c);
  for (size_t i = 0; i < st.atoms.size(); ++i) st.atoms[i].pos = F(r.positions[i]);
  DielectricOptions od = o;
  od.relax_cell = false, od.piezo = false;
  const DielectricResult rs = crystal_dielectric(st, od);
  double P = 0;
  const ForceField ff = zns_ff(st);
  for (size_t i = 0; i < st.atoms.size(); ++i) P += ff.charge[i] * st.cell.minimum_image(rs.positions[i] - st.atoms[i].pos)[0];
  const double e14 = P * 1.602176634e-19 * 1e-10 / (st.cell.volume() * 1e-30) / g;
  std::printf("e14 direct %.6f, linear response %.6f C/m²\n", e14, r.e[0][3]);
  EXPECT_NEAR(e14, r.e[0][3], 1e-3 * std::fabs(r.e[0][3]) + 1e-5);
}

// Phonons of a shell model (adiabatic): the shells condensed out at every k, three branches per ion, and the Γ optical
// modes those the dielectric module finds by its own condensation. A shell vibrating as a particle of its own gave core-
// shell modes of thousands of cm⁻¹ (MgO: 4063) and their zero-point energy, which bent the quasi-harmonic results.
TEST(Dielectric, ShellModelPhononsCondenseTheShells) {
  System s = salt(5.65);
  s.atoms.push_back(s.atoms[1]);   // the Cl shell at its core
  auto ff_for = [](const System& x) {
    ForceField ff = salt_ff(x);
    ff.type_names.push_back("Cls");
    ff.lj.push_back({0, 0});
    for (size_t i = 2; i < x.atoms.size(); i += 3) {
      ff.type_index[i] = 2, ff.atom_type[i] = "Cls", ff.charge[i] = -2.5, ff.charge[i - 1] = 1.5;
      ff.bonds.push_back({uint32_t(i - 1), uint32_t(i), 1500.0, 0.0});
      ff.excluded[i - 1].push_back(uint32_t(i)), ff.excluded[i].push_back(uint32_t(i - 1));
    }
    ff.pair_func.clear();
    ff.pair_func[{0, 2}] = {1, 26990.0, 0.3216, 0.0};
    ff.pair_func[{2, 2}] = {1, 50000.0, 0.33, 1500.0};
    return ff;
  };
  DielectricOptions d;
  d.energy = ewald();
  d.field_for = ff_for;
  d.shell = {0, 0, 1};
  d.points = 50, d.piezo = false;
  d.relax_cell = false;   // the phonons keep the cell as given
  const DielectricResult dr = crystal_dielectric(s, d);
  PhononOptions o;
  o.energy = ewald();
  o.field_for = ff_for;
  o.mesh = {4, 4, 4};
  o.path = {{0, 0, 0}, {0.5, 0, 0.5}};
  o.labels = {"G", "X"};
  o.path_points = 10;
  const PhononResult pr = phonons(s, o);
  ASSERT_FALSE(pr.bands.empty());
  ASSERT_EQ(pr.bands[0].size(), 6u);   // two ions, not three particles
  for (int k = 0; k < 3; ++k) EXPECT_NEAR(pr.bands[0][size_t(k)], 0, 1.0);   // the translations
  ASSERT_EQ(dr.modes.size(), 3u);
  for (int k = 0; k < 3; ++k) EXPECT_NEAR(pr.bands[0][size_t(3 + k)], dr.modes[size_t(k)], 0.01 * dr.modes[size_t(k)]) << k;
  for (const auto& b : pr.bands)
    for (double w : b) EXPECT_LT(w, 400);   // no core-shell mode
}

// Quasi-harmonic run of a charged shell-model crystal from the default (damped shifted force) Coulomb sum: switched to
// Ewald, since the cut-off sum puts kinks in G(V) as shells of ions cross it, and the core-shell modes condensed out. The
// fit is then smooth, the crystal expands as it warms, faster when hotter, and softens (in MgO the bulk modulus had risen
// with temperature, from both faults)
TEST(Dielectric, ShellModelQuasiHarmonicIsSmooth) {
  System s = salt(5.65);
  s.atoms.push_back(s.atoms[1]);
  auto ff_for = [](const System& x) {
    ForceField ff = salt_ff(x);
    ff.type_names.push_back("Cls");
    ff.lj.push_back({0, 0});
    for (size_t i = 2; i < x.atoms.size(); i += 3) {
      ff.type_index[i] = 2, ff.atom_type[i] = "Cls", ff.charge[i] = -2.5, ff.charge[i - 1] = 1.5;
      ff.bonds.push_back({uint32_t(i - 1), uint32_t(i), 1500.0, 0.0});
      ff.excluded[i - 1].push_back(uint32_t(i)), ff.excluded[i].push_back(uint32_t(i - 1));
    }
    ff.pair_func.clear();
    ff.pair_func[{0, 2}] = {1, 26990.0, 0.3216, 0.0};
    ff.pair_func[{2, 2}] = {1, 50000.0, 0.33, 1500.0};
    return ff;
  };
  QhaOptions o;
  o.phonon.field_for = ff_for;
  o.phonon.energy.cutoff = 8.0;   // the default electrostatics: damped shifted force
  o.phonon.mesh = {4, 4, 4};
  o.points = 7, o.range = 0.03;
  o.temperatures = {0, 200, 400, 600};
  const QhaResult r = quasi_harmonic(s, o);
  bool ewald = false;
  for (const auto& n : r.notes) ewald = ewald || n.find("Ewald") != std::string::npos;
  EXPECT_TRUE(ewald);
  EXPECT_LT(r.fit_rms, 0.01);
  ASSERT_EQ(r.rows.size(), 4u);
  for (size_t k = 1; k < r.rows.size(); ++k) {
    EXPECT_GT(r.rows[k].volume, r.rows[k - 1].volume) << k;
    EXPECT_LT(r.rows[k].B, r.rows[k - 1].B) << k;
    if (k > 1) EXPECT_GT(r.rows[k].alpha_V, r.rows[k - 1].alpha_V) << k;
  }
}

