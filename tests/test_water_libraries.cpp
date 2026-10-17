#include <gtest/gtest.h>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <tuple>
#include <vector>

#include "caps/ffdef.hpp"
#include "caps/field.hpp"
#include "caps/molecule.hpp"
#include "caps/relax.hpp"
#include "caps/typing.hpp"
#include "caps/water.hpp"

using namespace caps;

namespace {
constexpr double kKB = 0.0019872043;   // kcal/mol/K
const double kPi = std::acos(-1.0);

std::string lib(const std::string& id) { return std::string(CAPS_SOURCE_DIR) + "/data/forcefields/" + id + ".json"; }

// n waters of the given geometry, 3.1 Å apart along x, in a 3.1 n x 12 x 12 Å cell (or no cell)
System waters(int n, double r_oh, double theta_deg, bool cell) {
  System s;
  const double h = theta_deg * kPi / 360;
  for (int k = 0; k < n; ++k) {
    const Vec3 o{3.1 * k + 1.0, 6.0 + 0.3 * k, 6.0};
    const uint32_t b = uint32_t(s.atoms.size());
    Atom a;
    a.element = 8, a.pos = o, a.mol = k + 1, a.name = "O";
    s.atoms.push_back(a);
    a.element = 1, a.name = "H";
    a.pos = o + Vec3{r_oh * std::sin(h), 0, r_oh * std::cos(h)};
    s.atoms.push_back(a);
    a.pos = o + Vec3{-r_oh * std::sin(h), 0, r_oh * std::cos(h)};
    s.atoms.push_back(a);
    s.bonds.push_back({b, b + 1, 1});
    s.bonds.push_back({b, b + 2, 1});
  }
  for (size_t i = 0; i < s.atoms.size(); ++i) s.atoms[i].id = int64_t(i) + 1;
  s.has_mol = true;
  if (cell) s.cell.a = {3.1 * n, 0, 0}, s.cell.b = {0, 12, 0}, s.cell.c = {0, 0, 12};
  return s;
}

struct Typed { System s; ForceField ff; std::vector<std::string> types; };

Typed assign(const std::string& id, System s) {
  const FFDef def = load_forcefield(lib(id));
  std::string ch = "types";
  prepare_for_forcefield(s, def, ch);
  const TypingResult tr = assign_types(s, def);
  EXPECT_EQ(tr.untyped, 0) << id;
  ParamReport rep;
  ForceField ff = parameterize(s, def, tr.types, ch, &rep, false);
  EXPECT_TRUE(rep.missing.empty()) << id;
  return {s, ff, tr.types};
}

// the published numbers (bench/ff/make_water_models.py names where each was read)
struct Model { const char* id; const char* caps_model; double q_h, r_oh, theta, sigma, eps, d_om, d_ol, theta_l, eps_h, sigma_h; };
// SPC and SPC/E: -(A/r)^6 + (B/r)^12, A 0.37122 (kJ/mol)^(1/6) nm, B 0.3428 (kJ/mol)^(1/12) nm (Berendsen 1987, eq 18);
// TIP3P and TIP4P: A/r^12 - C/r^6 in kcal/mol and A (Jorgensen 1983, Table I)
const double spc_c6 = std::pow(0.37122, 6), spc_c12 = std::pow(0.3428, 12);
const double spc_sigma = 10 * std::pow(spc_c12 / spc_c6, 1.0 / 6), spc_eps = spc_c6 * spc_c6 / (4 * spc_c12) / 4.184;
const double rmin = std::pow(2.0, 1.0 / 6);
const std::vector<Model> kModels = {
    {"spc-berendsen1981", "spc", 0.41, 1.0, 109.47, spc_sigma, spc_eps, 0, 0, 0, 0, 0},
    {"spce-berendsen1987", "spce", 0.4238, 1.0, 109.47, spc_sigma, spc_eps, 0, 0, 0, 0, 0},
    {"tip3p-jorgensen1983", "tip3p", 0.417, 0.9572, 104.52, std::pow(582.0e3 / 595.0, 1.0 / 6), 595.0 * 595.0 / (4 * 582.0e3), 0, 0, 0, 0, 0},
    {"tip3p-charmm-mackerell1998", "tip3p-charmm", 0.417, 0.9572, 104.52, 2 * 1.7682 / rmin, 0.1521, 0, 0, 0, 0.046, 2 * 0.2245 / rmin},
    {"tip3p-pme-price2004", "tip3p-ew", 0.415, 0.9572, 104.52, 3.188, 0.102, 0, 0, 0, 0, 0},
    {"tip4p-jorgensen1983", "tip4p", 0.520, 0.9572, 104.52, std::pow(600.0e3 / 610.0, 1.0 / 6), 610.0 * 610.0 / (4 * 600.0e3), 0.15, 0, 0, 0, 0},
    {"tip4pew-horn2004", "tip4p-ew", 0.52422, 0.9572, 104.52, 3.16435, 0.162750, 0.125, 0, 0, 0, 0},
    {"tip4p2005-abascal2005", "tip4p2005", 0.5564, 0.9572, 104.52, 3.1589, 93.2 * kKB, 0.1546, 0, 0, 0, 0},
    {"tip4pice-abascal2005", "tip4p-ice", 0.5897, 0.9572, 104.52, 3.1668, 106.1 * kKB, 0.1577, 0, 0, 0, 0},
    {"tip5p-mahoney2000", "tip5p", 0.241, 0.9572, 104.52, 3.12, 0.16, 0, 0.70, 109.47, 0, 0},
    {"spcfw-wu2006", "spcfw", 0.41, 1.012, 113.24, 3.165492, 0.1554253, 0, 0, 0, 0, 0},
    {"opc-izadi2014", "opc", 0.6791, 0.8724, 103.6, 3.16655, 0.89036 / 4.184, 0.1594, 0, 0, 0, 0},
    {"opc3-izadi2016", "opc3", 0.447585, 0.97888, 109.47, 3.17427, 0.68369 / 4.184, 0, 0, 0, 0, 0},
};
}  // namespace

// Every water library built from its source: typed waters carry the model's charges (each molecule neutral), its
// geometry in the O-H bond and H-O-H angle, its Lennard-Jones, and its massless sites where the model puts them (M d_OM
// from O on the bisector, TIP5P's lone pairs d_OL from O at L-O-L, away from the hydrogens); CAPS's built-in water
// models (water_models(), the Solvation page and water groups) carry the same numbers
TEST(WaterLibraries, ModelsAsPublished) {
  for (const auto& m : kModels) {
    const int extra = m.d_om > 0 ? 1 : m.d_ol > 0 ? 2 : 0;
    Typed t = assign(m.id, waters(2, m.r_oh, m.theta, false));
    ASSERT_EQ(t.s.atoms.size(), size_t(2 * (3 + extra))) << m.id;
    ASSERT_EQ(t.ff.vsites.size(), size_t(2 * extra)) << m.id;
    for (int w = 0; w < 2; ++w) {
      double q = 0;
      int o = -1;
      std::vector<int> h, x;
      for (size_t i = 0; i < t.s.atoms.size(); ++i) {
        if (t.s.atoms[i].mol != w + 1) continue;
        q += t.ff.charge[i];
        if (t.s.atoms[i].element == 8) o = int(i);
        if (t.s.atoms[i].element == 1) h.push_back(int(i));
        if (t.s.atoms[i].element == 0) x.push_back(int(i));
      }
      ASSERT_GE(o, 0);
      ASSERT_EQ(h.size(), 2u);
      ASSERT_EQ(int(x.size()), extra) << m.id;
      EXPECT_NEAR(q, 0.0, 1e-12) << m.id;
      EXPECT_NEAR(t.ff.charge[size_t(h[0])], m.q_h, 1e-12) << m.id;
      EXPECT_NEAR(t.ff.charge[size_t(o)], extra ? 0.0 : -2 * m.q_h, 1e-12) << m.id;
      const auto& lo = t.ff.lj[size_t(t.ff.type_index[size_t(o)])];
      const auto& lh = t.ff.lj[size_t(t.ff.type_index[size_t(h[0])])];
      EXPECT_NEAR(lo.sigma, m.sigma, 1e-12) << m.id;
      EXPECT_NEAR(lo.eps, m.eps, 1e-12) << m.id;
      EXPECT_NEAR(lh.eps, m.eps_h, 1e-12) << m.id;
      if (m.eps_h > 0) EXPECT_NEAR(lh.sigma, m.sigma_h, 1e-12) << m.id;
      const Vec3 O = t.s.atoms[size_t(o)].pos, H1 = t.s.atoms[size_t(h[0])].pos, H2 = t.s.atoms[size_t(h[1])].pos;
      if (m.d_om > 0) {
        const Vec3 M = t.s.atoms[size_t(x[0])].pos;
        EXPECT_NEAR(norm(M - O), m.d_om, 1e-9) << m.id;
        EXPECT_NEAR(norm(cross(M - O, (H1 - O) + (H2 - O))), 0.0, 1e-9) << m.id << ": M on the bisector";
        EXPECT_NEAR(t.ff.charge[size_t(x[0])], -2 * m.q_h, 1e-12) << m.id;
      }
      if (m.d_ol > 0) {
        const Vec3 L1 = t.s.atoms[size_t(x[0])].pos, L2 = t.s.atoms[size_t(x[1])].pos;
        EXPECT_NEAR(norm(L1 - O), m.d_ol, 1e-9) << m.id;
        EXPECT_NEAR(norm(L2 - O), m.d_ol, 1e-9) << m.id;
        EXPECT_NEAR(std::acos(dot(L1 - O, L2 - O) / (m.d_ol * m.d_ol)) * 180 / kPi, m.theta_l, 1e-7) << m.id;
        EXPECT_LT(dot((L1 - O) + (L2 - O), (H1 - O) + (H2 - O)), 0) << m.id << ": away from the hydrogens";
        EXPECT_NEAR(t.ff.charge[size_t(x[0])], -m.q_h, 1e-12) << m.id;
      }
    }
    for (const auto& b : t.ff.bonds) EXPECT_NEAR(b.r0, m.r_oh, 1e-12) << m.id;
    for (const auto& a : t.ff.angles) EXPECT_NEAR(a.theta0 * 180 / kPi, m.theta, 1e-9) << m.id;
    // the built-in model of the same name: the same numbers
    const WaterModel& w = water_model(m.caps_model);
    EXPECT_EQ(w.sites, 3 + extra) << m.id;
    EXPECT_NEAR(w.q_h, m.q_h, 1e-12) << m.id;
    EXPECT_NEAR(w.r_oh, m.r_oh, 1e-12) << m.id;
    EXPECT_NEAR(w.theta, m.theta, 1e-12) << m.id;
    EXPECT_NEAR(w.sigma_o, m.sigma, 1e-12) << m.id;
    EXPECT_NEAR(w.eps_o, m.eps, 1e-12) << m.id;
    EXPECT_NEAR(w.eps_h, m.eps_h, 1e-12) << m.id;
    EXPECT_NEAR(w.d_om, m.d_om, 1e-12) << m.id;
    EXPECT_NEAR(w.d_ol, m.d_ol, 1e-12) << m.id;
  }
}

// The four- and five-site libraries export: LAMMPS's tip4p styles with the model's d_OM (M leaves the data file, its
// charge on O), GROMACS's virtual sites (TIP5P's lone pairs as 3out sites); the energy of the typed waters equals that
// of CAPS's built-in model on the same structure
TEST(WaterLibraries, SitesExportAndEnergy) {
  const auto dir = std::filesystem::temp_directory_path() / "caps_water_libraries";
  std::filesystem::create_directories(dir);
  EnergyOptions e;
  e.electrostatics = EnergyOptions::Electrostatics::PME;
  e.cutoff = 6.0;
  {
    Typed t = assign("tip4p2005-abascal2005", waters(4, 0.9572, 104.52, true));
    write_lammps_data_ff(t.s, t.ff, e, (dir / "w.data").string(), false);
    write_lammps_input(t.s, t.ff, e, "w.data", (dir / "w.in").string(), 0, true, {}, {}, nullptr);
    std::ifstream dat(dir / "w.data"), inp(dir / "w.in");
    std::stringstream ds, is;
    ds << dat.rdbuf(), is << inp.rdbuf();
    EXPECT_NE(ds.str().find("12 atoms"), std::string::npos);
    EXPECT_NE(ds.str().find("-1.112800"), std::string::npos) << "M's charge on O";
    EXPECT_NE(is.str().find("0.154600"), std::string::npos) << is.str();
    EXPECT_NE(is.str().find("pppm/tip4p"), std::string::npos);
    // the same waters through the built-in model: the same energy
    System s2 = waters(4, 0.9572, 104.52, true);
    const WaterModel& wm = water_model("tip4p2005");
    apply_water_model(s2, wm);
    std::vector<uint32_t> all(s2.atoms.size());
    for (size_t i = 0; i < all.size(); ++i) all[i] = uint32_t(i);
    const ForceField F2 = water_forcefield(s2, wm, all);
    auto energy = [&](const System& s, const ForceField& ff) {
      Evaluator ev(ff, e);
      std::vector<double> x(3 * s.atoms.size()), f;
      for (size_t i = 0; i < s.atoms.size(); ++i) for (int c = 0; c < 3; ++c) x[3 * i + c] = s.atoms[i].pos[c];
      return ev.compute(x, s.cell, f).total();
    };
    EXPECT_NEAR(energy(t.s, t.ff), energy(s2, F2), 1e-6);
  }
  {
    Typed t = assign("tip5p-mahoney2000", waters(4, 0.9572, 104.52, true));
    const std::string stem = (dir / "t5").string();
    write_gromacs(t.s, t.ff, e, stem);
    std::string all;
    for (const char* ext : {".top", ".itp"}) {
      std::ifstream f(stem + ext);
      std::stringstream ss;
      ss << f.rdbuf();
      all += ss.str();
    }
    EXPECT_NE(all.find("virtual_sites3"), std::string::npos) << all;
    EXPECT_EQ(t.ff.vsites.size(), 8u);
    for (const auto& v : t.ff.vsites) EXPECT_NE(v.c, 0.0);
  }
  std::filesystem::remove_all(dir);
}

// mW (Molinero & Moore 2009): one Stillinger-Weber site per water at its oxygen, the paper's lambda, epsilon and sigma
TEST(WaterLibraries, MonatomicWater) {
  Typed t = assign("mw-molinero2009", waters(3, 0.9572, 104.52, true));
  ASSERT_EQ(t.s.atoms.size(), 3u);
  for (const auto& ty : t.types) EXPECT_EQ(ty, "MW");
  EXPECT_NEAR(t.ff.mass[0], 15.9994 + 2 * 1.008, 1e-9);
  const FFDef def = load_forcefield(lib("mw-molinero2009"));
  ASSERT_EQ(def.pairs.size(), 1u);
  const auto& p = def.pairs[0].params;
  ASSERT_GE(p.size(), 10u);
  EXPECT_DOUBLE_EQ(p[0], 6.189);
  EXPECT_DOUBLE_EQ(p[1], 2.3925);
  EXPECT_DOUBLE_EQ(p[3], 23.15);
}

// Noble gases (Rutkai, Thol, Span & Vrabec 2017, Table 1): Ne, Ar, Kr, Xe as Lennard-Jones atoms with the paper's
// epsilon/k_B and sigma, unlike pairs by Lorentz-Berthelot
TEST(GasLibraries, NobleGasesAsPublished) {
  System s;
  const char* el[4] = {"Ne", "Ar", "Kr", "Xe"};
  const int z[4] = {10, 18, 36, 54};
  const double ek[4] = {33.921, 116.79, 162.58, 226.14}, sg[4] = {2.801, 3.3952, 3.6274, 3.949};
  for (int k = 0; k < 4; ++k) {
    Atom a;
    a.element = z[k], a.pos = {5.0 * k, 0, 0}, a.mol = k + 1, a.name = el[k], a.id = k + 1;
    s.atoms.push_back(a);
  }
  s.has_mol = true;
  Typed t = assign("noble-gases-rutkai2017", s);
  for (int k = 0; k < 4; ++k) {
    EXPECT_EQ(t.types[size_t(k)], el[k]);
    const auto& lj = t.ff.lj[size_t(t.ff.type_index[size_t(k)])];
    EXPECT_NEAR(lj.eps, ek[k] * kKB, 1e-12) << el[k];
    EXPECT_NEAR(lj.sigma, sg[k], 1e-12) << el[k];
    EXPECT_EQ(t.ff.charge[size_t(k)], 0.0);
  }
}

// The Cooke-Kremer-Deserno lipid (Phys. Rev. E 72, 011506 (2005)): the paper's FENE bonds (k 30, r_inf 1.5), the
// straightening spring (k_bend 10 for (1/2) k (r - 4)^2) and its pair terms (WCA b 0.95 / 1, cos^2 attraction of width 1.6)
TEST(CgLibraries, CookeDesernoAsPublished) {
  const FFDef def = load_forcefield(lib("cooke-deserno2005"));
  const double r6 = std::pow(2.0, 1.0 / 6);
  ASSERT_EQ(def.pairs.size(), 3u);
  for (const auto& p : def.pairs) {
    const bool tt = p.name == "T-T";
    EXPECT_NEAR(p.params[1], (tt ? 1.0 : 0.95) * r6, 1e-12) << p.name;
    EXPECT_NEAR(p.params[2], (tt ? r6 + 1.6 : 0.95 * r6), 1e-12) << p.name;
  }
  int fene = 0;
  for (const auto& b : def.bonds) {
    if (b.style == "fene") {
      EXPECT_DOUBLE_EQ(b.params[0], 30.0);
      EXPECT_DOUBLE_EQ(b.params[1], 1.5);
      ++fene;
    } else {
      EXPECT_DOUBLE_EQ(b.params[0], 5.0);
      EXPECT_DOUBLE_EQ(b.params[1], 4.0);
    }
  }
  EXPECT_EQ(fene, 2);
}

// CO2 (EPM2, Harris & Yung 1995, Table 1): the charges, the paper's own C-O cross term, the 1.149 A bond and the bend; and
// H2 (Buch 1994; Darkrim & Levesque 1998): the two hydrogens with a massless centre midway carrying the Lennard-Jones
// (and Darkrim-Levesque's -2q against q on each proton, q = 2 Q / d^2 from the paper's quadrupole 0.615e-26 esu cm^2)
TEST(GasLibraries, Co2AndH2AsPublished) {
  {
    System s;
    const Vec3 c{5, 5, 5};
    for (int k = -1; k <= 1; ++k) {
      Atom a;
      a.element = k == 0 ? 6 : 8, a.pos = c + Vec3{1.149 * k, 0, 0}, a.mol = 1, a.id = int64_t(s.atoms.size()) + 1;
      s.atoms.push_back(a);
    }
    s.bonds = {{1, 0, 2}, {1, 2, 2}};
    s.has_mol = true;
    Typed t = assign("epm2-harris1995", s);
    EXPECT_NEAR(t.ff.charge[1], 0.6512, 1e-12);
    EXPECT_NEAR(t.ff.charge[0] + t.ff.charge[1] + t.ff.charge[2], 0.0, 1e-12);
    const FFDef def = load_forcefield(lib("epm2-harris1995"));
    bool cross = false;
    for (const auto& p : def.pairs)
      if (p.match.size() == 2) {
        cross = true;
        EXPECT_NEAR(p.params[0], 47.588 * kKB, 1e-12);
        EXPECT_NEAR(p.params[1], 2.892, 1e-12);
      }
    EXPECT_TRUE(cross);
    for (const auto& b : t.ff.bonds) EXPECT_NEAR(b.r0, 1.149, 1e-12);
    ASSERT_EQ(t.ff.angles.size(), 1u);
    EXPECT_NEAR(t.ff.angles[0].theta0 * 180 / kPi, 180.0, 1e-9);
  }
  const double e_esu = 4.803204712570263e-10, q_dl = 2 * 0.615e-26 / std::pow(0.741e-8, 2) / e_esu;
  for (const auto& [id, ek, sg, q] : std::vector<std::tuple<std::string, double, double, double>>{
           {"h2-buch1994", 34.2, 2.96, 0.0}, {"h2-darkrim1998", 36.7, 2.958, q_dl}}) {
    System s;
    for (int k = 0; k < 2; ++k) {
      Atom a;
      a.element = 1, a.pos = {5 + 0.741 * k, 5, 5}, a.mol = 1, a.id = k + 1;
      s.atoms.push_back(a);
    }
    s.bonds = {{0, 1, 1}};
    s.has_mol = true;
    Typed t = assign(id, s);
    ASSERT_EQ(t.s.atoms.size(), 3u) << id;
    ASSERT_EQ(t.ff.vsites.size(), 1u) << id;
    const uint32_t m = t.ff.vsites[0].site;
    EXPECT_NEAR(norm(t.s.atoms[m].pos - Vec3{5.3705, 5, 5}), 0.0, 1e-9) << id;
    const auto& lj = t.ff.lj[size_t(t.ff.type_index[m])];
    EXPECT_NEAR(lj.eps, ek * kKB, 1e-12) << id;
    EXPECT_NEAR(lj.sigma, sg, 1e-12) << id;
    EXPECT_NEAR(t.ff.charge[m], -2 * q, 1e-12) << id;
    EXPECT_NEAR(t.ff.charge[0], q, 1e-12) << id;
  }
  EXPECT_NEAR(q_dl, 0.4664, 1e-4);
}

// L-OPLS (Siu, Pluhackova & Boeckmann 2012, Table 2) on OPLS-AA 2024: a chain's CH3, CH2 and alkene CH with their
// L-OPLS charges and Lennard-Jones (each group neutral), the CH2 next to a carboxyl left to OPLS-AA 2024, and the
// refitted CT-CT-CT-CT torsion (Fourier 2.697394, -0.896807, 0.74567 kJ/mol)
TEST(OrganicLibraries, LoplsAsPublished) {
  std::map<std::string, int> seen;
  for (const char* smi : {"CCCCCCCCCCCCCCCC", "CCCCCCCC/C=C\\CCCCCCCC(=O)O"}) {
    Typed t = assign("lopls-siu2012", build_molecule(smi).system);
    double q = 0;
    for (size_t i = 0; i < t.s.atoms.size(); ++i) {
      q += t.ff.charge[i];
      seen[t.types[i]]++;
      const std::string& ty = t.types[i];
      const auto& lj = t.ff.lj[size_t(t.ff.type_index[i])];
      if (ty == "lopls_135") EXPECT_NEAR(t.ff.charge[i], -0.222, 1e-12);
      if (ty == "lopls_136") EXPECT_NEAR(t.ff.charge[i], -0.148, 1e-12);
      if (ty == "lopls_142") EXPECT_NEAR(t.ff.charge[i], -0.160, 1e-12);
      if (ty == "lopls_144") EXPECT_NEAR(t.ff.charge[i], 0.160, 1e-12);
      if (ty == "lopls_140b") EXPECT_NEAR(lj.eps, 0.110 / 4.184, 1e-12);
      if (ty == "lopls_136") EXPECT_NEAR(lj.sigma, 3.5, 1e-12);
    }
    EXPECT_NEAR(q, 0.0, 1e-9) << smi;
    if (std::string(smi).size() == 16) {
      int ct = 0;
      for (const auto& d : t.ff.dihedrals) {
        if (t.types[d.i][0] == 'H' || t.types[d.l][0] == 'H' || t.types[d.j] == "lopls_135" || t.types[d.k] == "lopls_135") continue;
        if (t.s.atoms[d.i].element != 6 || t.s.atoms[d.l].element != 6) continue;
        // opls K1..K3 as cosine terms: n = 1, 2, 3 with V = K (kcal/mol)
        if (d.n == 1) {
          EXPECT_NEAR(d.v, 2.697394 / 4.184 / 2, 1e-9);
          ++ct;
        }
      }
      EXPECT_GT(ct, 0);
    }
  }
  EXPECT_GT(seen["lopls_135"], 0);
  EXPECT_GT(seen["lopls_136"], 0);
  EXPECT_GT(seen["lopls_140"], 0);
  EXPECT_GT(seen["lopls_140b"], 0);
  EXPECT_GT(seen["lopls_142"], 0);
  EXPECT_GT(seen["lopls_144"], 0);
  EXPECT_GT(seen["opls_136"], 0);   // the CH2 next to the carboxyl stays OPLS-AA 2024's
}
