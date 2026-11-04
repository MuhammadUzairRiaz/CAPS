#include <gtest/gtest.h>

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <map>
#include <memory>
#include <sstream>

#include "caps/dynamics.hpp"
#include "caps/relax.hpp"
#include "caps/restraints.hpp"
#include "caps/uff.hpp"

using namespace caps;

namespace {
// the four atoms gmx 2026.3 evaluated (nm × 10), in a 30 Å cell
System four_atoms() {
  System s;
  const Vec3 p[4] = {{10.0, 10.0, 10.0}, {11.8, 10.0, 10.0}, {11.8, 11.5, 10.0}, {12.06, 11.5, 8.52}};
  for (int i = 0; i < 4; ++i) {
    Atom a;
    a.element = 18;
    a.mol = i + 1;
    a.pos = p[i];
    s.atoms.push_back(a);
  }
  s.has_mol = true;
  s.cell.a = {30, 0, 0}, s.cell.b = {0, 30, 0}, s.cell.c = {0, 0, 30};
  return s;
}
std::vector<double> coords(const System& s) {
  std::vector<double> x;
  for (const auto& a : s.atoms) x.insert(x.end(), a.pos.begin(), a.pos.end());
  return x;
}
double energy(const std::vector<Restraint>& rs, const std::vector<double>& x, const Cell& c) {
  std::vector<double> f(x.size(), 0.0);
  return add_restraint_forces(rs, x, c, f).energy;
}
constexpr double kJ = 4.184;
}  // namespace

// The same restraints GROMACS 2026.3 evaluated (bonds function 10 low 0.10 up1 0.15 nm, k 1000; an angle 100°, k 300;
// a dihedral restraint φ0 60°, Δφ 10°, k 200; a position restraint k_x 500; a flat-bottomed sphere r 0.05 nm, k 700):
// CAPS's energies equal gmx's 0.449999, 4.56926, 685.059, 0.0999998 and 0.315001 kJ/mol to its single precision
TEST(Restraints, EnergiesEqualGromacs) {
  const System s = four_atoms();
  const auto x = coords(s);
  Restraint d;
  d.kind = Restraint::Kind::Distance, d.atoms = {0, 1, 0, 0}, d.lower = 1.0, d.upper = 1.5, d.k = 1000 / (kJ * 100);
  EXPECT_NEAR(energy({d}, x, s.cell) * kJ, 0.449999, 2e-6);
  Restraint a;
  a.kind = Restraint::Kind::Angle, a.atoms = {0, 1, 2, 0}, a.lower = a.upper = 100, a.k = 300 / kJ;
  EXPECT_NEAR(energy({a}, x, s.cell) * kJ, 4.56926, 2e-5);
  Restraint t;
  t.kind = Restraint::Kind::Torsion, t.atoms = {0, 1, 2, 3}, t.lower = 50, t.upper = 70, t.k = 200 / kJ;
  EXPECT_NEAR(energy({t}, x, s.cell) * kJ, 685.059, 2e-3);
  EXPECT_NEAR(restraint_value(t, x, s.cell), -99.96, 0.01);   // IUPAC sign
  Restraint p;
  p.kind = Restraint::Kind::Position, p.atoms = {3, 0, 0, 0}, p.anchor = {12.26, 11.8, 8.52}, p.k = 500 / (kJ * 100), p.axes = 1;
  EXPECT_NEAR(energy({p}, x, s.cell) * kJ, 0.0999998, 2e-7);
  Restraint fb;
  fb.kind = Restraint::Kind::Position, fb.atoms = {2, 0, 0, 0}, fb.anchor = {11.8, 12.3, 10.0}, fb.upper = 0.5, fb.k = 700 / (kJ * 100);
  EXPECT_NEAR(energy({fb}, x, s.cell) * kJ, 0.315001, 2e-6);
}

// Forces are minus the gradient of the energy (central differences) for every kind, harmonic and flat-bottomed, across
// the cell boundary; the virial of the internal restraints is Σ r·f
TEST(Restraints, ForcesAreGradients) {
  System s = four_atoms();
  for (auto& a : s.atoms) a.pos[0] += 19.0;   // atoms 2–4 cross x = 30
  auto x = coords(s);
  for (size_t i = 0; i < x.size(); i += 3) x[i] = std::fmod(x[i], 30.0);
  std::vector<Restraint> rs(7);
  rs[0].kind = Restraint::Kind::Distance, rs[0].atoms = {0, 1, 0, 0}, rs[0].lower = rs[0].upper = 2.2, rs[0].k = 30;
  rs[1].kind = Restraint::Kind::Distance, rs[1].atoms = {0, 2, 0, 0}, rs[1].lower = 1.0, rs[1].upper = 1.5, rs[1].k = 20;
  rs[2].kind = Restraint::Kind::Angle, rs[2].atoms = {0, 1, 2, 0}, rs[2].lower = 100, rs[2].upper = 110, rs[2].k = 60;
  rs[3].kind = Restraint::Kind::Torsion, rs[3].atoms = {0, 1, 2, 3}, rs[3].lower = rs[3].upper = 150, rs[3].k = 40;
  rs[4].kind = Restraint::Kind::Torsion, rs[4].atoms = {3, 2, 1, 0}, rs[4].lower = 170, rs[4].upper = 200, rs[4].k = 40;   // across ±180°
  rs[5].kind = Restraint::Kind::Position, rs[5].atoms = {3, 0, 0, 0}, rs[5].anchor = {1.0, 11.0, 9.0}, rs[5].k = 15, rs[5].axes = 5;
  rs[6].kind = Restraint::Kind::Position, rs[6].atoms = {2, 0, 0, 0}, rs[6].anchor = {29.5, 12.0, 10.5}, rs[6].k = 15, rs[6].upper = 0.3;
  for (const auto& r : rs) {
    std::vector<double> f(x.size(), 0.0);
    const auto t = add_restraint_forces({r}, x, s.cell, f);
    EXPECT_GT(t.energy, 1e-3) << to_string(r.kind);
    for (size_t i = 0; i < x.size(); ++i) {
      auto xp = x, xm = x;
      xp[i] += 1e-5, xm[i] -= 1e-5;
      const double num = -(energy({r}, xp, s.cell) - energy({r}, xm, s.cell)) / 2e-5;
      EXPECT_NEAR(f[i], num, 1e-5 * std::max(1.0, std::fabs(num))) << to_string(r.kind) << " coordinate " << i;
    }
  }
  // virial: with the atoms made whole (no image crossing), Σ r·f
  const System w = four_atoms();
  const auto xw = coords(w);
  std::vector<double> f(xw.size(), 0.0);
  const auto t = add_restraint_forces({rs[0], rs[2], rs[3]}, xw, w.cell, f);
  double vir = 0;
  for (size_t i = 0; i < xw.size(); ++i) vir += xw[i] * f[i];
  EXPECT_NEAR(t.virial, vir, 1e-9 * std::max(1.0, std::fabs(vir)));
}

// The text form: atoms from 1, bounds LO:HI, force constants, ramps, position selections with axes and radius; and back
TEST(Restraints, TextForm) {
  const System s = four_atoms();
  const auto rs = parse_restraints("distance 1 2 2.5 20 ramp 0; angle 1 2 3 100:110\ntorsion 1 2 3 4 180; position 3-4 5 axes xz radius 0.2", s);
  ASSERT_EQ(rs.size(), 5u);
  EXPECT_EQ(rs[0].kind, Restraint::Kind::Distance);
  EXPECT_EQ(rs[0].atoms[1], 1u);
  EXPECT_DOUBLE_EQ(rs[0].k, 20);
  EXPECT_DOUBLE_EQ(rs[0].k_end, 0);
  EXPECT_DOUBLE_EQ(rs[1].lower, 100);
  EXPECT_DOUBLE_EQ(rs[1].upper, 110);
  EXPECT_DOUBLE_EQ(rs[1].k, 50);
  EXPECT_EQ(rs[3].axes, 5);
  EXPECT_DOUBLE_EQ(rs[4].upper, 0.2);
  EXPECT_EQ(rs[4].anchor, s.atoms[3].pos);
  const auto back = parse_restraints(restraints_text(rs), s);
  ASSERT_EQ(back.size(), rs.size());
  for (size_t k = 0; k < rs.size(); ++k) {
    EXPECT_DOUBLE_EQ(back[k].upper, rs[k].upper);
    EXPECT_EQ(back[k].atoms, rs[k].atoms);
  }
  EXPECT_THROW(parse_restraints("distance 1 9 2.0", s), std::invalid_argument);
  EXPECT_THROW(parse_restraints("bend 1 2 3 90", s), std::invalid_argument);
}

// Dynamics: an NVE run holds its total energy with the restraint energy in the potential, and a minimisation ends at the
// restrained distance
TEST(Restraints, DynamicsConservesAndMinimisationHolds) {
  System s;
  for (int i = 0; i < 2; ++i) {
    Atom a;
    a.element = 18, a.mol = i + 1, a.pos = {10.0 + 7.0 * i, 10, 10};
    s.atoms.push_back(a);
  }
  s.has_mol = true;
  s.cell.a = {40, 0, 0}, s.cell.b = {0, 40, 0}, s.cell.c = {0, 0, 40};
  RelaxOptions ro;
  ro.pushoff = false, ro.ftol = 1e-4, ro.max_iterations = 2000;
  ro.restraint_terms = parse_restraints("distance 1 2 5.0 100", s);
  System m = s;
  relax(m, ro);
  EXPECT_NEAR(norm(m.atoms[1].pos - m.atoms[0].pos), 5.0, 0.01);
  DynamicsOptions o;
  o.thermostat = Thermostat::None;
  o.temperature = 300, o.steps = 2000, o.dt = 1, o.frame_every = 0, o.thermo_every = 100;
  o.restraints = parse_restraints("distance 1 2 4.0 20; position 1 2", s);
  o.energy.coulomb = false;
  DynamicsReport rep;
  run_dynamics(s, o, &rep);
  ASSERT_GT(rep.thermo.size(), 3u);
  double emax = 0, rmax = 0;
  for (const auto& r : rep.thermo) emax = std::max(emax, std::fabs(r.total - rep.thermo[0].total)), rmax = std::max(rmax, r.restraint);
  EXPECT_GT(rmax, 1.0);   // 7 Å held at 4: the restraint stores energy
  EXPECT_LT(emax, 1e-3 * rmax);
  std::printf("NVE with restraints: largest |ΔE| %.2e kcal/mol, restraint energy up to %.2f\n", emax, rmax);
}

namespace {
bool have(const char* tool) { return std::system((std::string("command -v ") + tool + " >/dev/null 2>&1").c_str()) == 0; }
std::string slurp(const std::filesystem::path& p) {
  std::ifstream f(p);
  return std::string((std::istreambuf_iterator<char>(f)), {});
}
// "Energies (kJ/mol)" block of a GROMACS log: names in 15-character columns over their values
std::map<std::string, double> gmx_energies(const std::string& log) {
  std::map<std::string, double> e;
  const auto at = log.find("   Energies (kJ/mol)");
  if (at == std::string::npos) return e;
  std::istringstream in(log.substr(at));
  std::string line, head;
  std::getline(in, line);
  while (std::getline(in, head) && !head.empty() && std::getline(in, line)) {
    for (size_t c = 0; c + 15 <= head.size() + 1; c += 15) {
      std::string name = head.substr(c, 15);
      name.erase(0, name.find_first_not_of(' '));
      if (name.empty() || c >= line.size()) continue;
      e[name] = std::stod(line.substr(c, 15));
    }
  }
  return e;
}
}  // namespace

// The GROMACS files CAPS writes, run by gmx (when installed): its restraint energies equal CAPS's; the restrained atoms'
// molecules become their own molecule types
TEST(Restraints, GromacsExportRunsTheSameRestraints) {
  if (!have("gmx")) GTEST_SKIP() << "gmx not installed";
  System s = four_atoms();
  const ForceField ff = default_forcefield(s);
  std::vector<Restraint> rs(5);
  rs[0].kind = Restraint::Kind::Distance, rs[0].atoms = {0, 1, 0, 0}, rs[0].lower = 1.0, rs[0].upper = 1.5, rs[0].k = 2.5;
  rs[1].kind = Restraint::Kind::Angle, rs[1].atoms = {0, 1, 2, 0}, rs[1].lower = rs[1].upper = 100, rs[1].k = 70;
  rs[2].kind = Restraint::Kind::Torsion, rs[2].atoms = {0, 1, 2, 3}, rs[2].lower = 50, rs[2].upper = 70, rs[2].k = 40;
  rs[3].kind = Restraint::Kind::Position, rs[3].atoms = {3, 0, 0, 0}, rs[3].anchor = {12.26, 11.8, 8.52}, rs[3].k = 1.2, rs[3].axes = 3;
  rs[4].kind = Restraint::Kind::Position, rs[4].atoms = {2, 0, 0, 0}, rs[4].anchor = {11.8, 12.3, 10.0}, rs[4].upper = 0.5, rs[4].k = 1.7;
  const auto dir = std::filesystem::temp_directory_path() / "caps_restraints_gmx";
  std::filesystem::remove_all(dir);
  std::filesystem::create_directories(dir);
  const std::string stem = (dir / "r").string();
  EnergyOptions e;
  e.cutoff = 10, e.coulomb = false;
  const auto notes = write_gromacs(s, ff, e, stem, rs);
  const std::string itp = slurp(stem + ".itp");
  EXPECT_NE(itp.find("[ position_restraints ]"), std::string::npos);
  const std::string cmd = "cd '" + dir.string() + "' && gmx grompp -f r.mdp -c r.gro -r r_ref.gro -p r.top -o r.tpr -maxwarn 10 > grompp.txt 2>&1 && gmx mdrun -s r.tpr -deffnm r -nt 1 > mdrun.txt 2>&1";
  ASSERT_EQ(std::system(cmd.c_str()), 0) << slurp(dir / "grompp.txt");
  const auto ge = gmx_energies(slurp(dir / "r.log"));
  const auto x = coords(s);
  auto caps_kj = [&](std::initializer_list<int> which) {
    std::vector<Restraint> sub;
    for (int k : which) sub.push_back(rs[size_t(k)]);
    return energy(sub, x, s.cell) * kJ;
  };
  ASSERT_TRUE(ge.count("Restraint Pot.") && ge.count("Dih. Rest.") && ge.count("Position Rest.") && ge.count("Flat-b. P-R.")) << slurp(dir / "r.log");
  EXPECT_NEAR(ge.at("Restraint Pot."), caps_kj({0}), 1e-4 * std::max(1.0, caps_kj({0})));
  EXPECT_NEAR(ge.at("Angle"), caps_kj({1}), 1e-4 * std::max(1.0, caps_kj({1})));
  EXPECT_NEAR(ge.at("Dih. Rest."), caps_kj({2}), 1e-4 * std::max(1.0, caps_kj({2})));
  EXPECT_NEAR(ge.at("Position Rest."), caps_kj({3}), 1e-4 * std::max(1.0, caps_kj({3})));
  EXPECT_NEAR(ge.at("Flat-b. P-R."), caps_kj({4}), 1e-4 * std::max(1.0, caps_kj({4})));
  std::printf("gmx: distance %.6f angle %.6f dihedral %.6f position %.6f flat %.6f kJ/mol\n", ge.at("Restraint Pot."), ge.at("Angle"), ge.at("Dih. Rest."),
              ge.at("Position Rest."), ge.at("Flat-b. P-R."));
}

// The LAMMPS input CAPS writes, run by LAMMPS (when installed): fix restrain gives K (r − r0)² = ½ k Δr² for a bond and
// an angle, and K [1 − cos(φ − φ0)] for the torsion with φ0 written 180° over for LAMMPS's reckoning
TEST(Restraints, LammpsExportRunsTheSameRestraints) {
  if (!have("lmp")) GTEST_SKIP() << "lmp not installed";
  System s = four_atoms();
  const ForceField ff = default_forcefield(s);
  std::vector<Restraint> rs(3);
  rs[0].kind = Restraint::Kind::Distance, rs[0].atoms = {0, 1, 0, 0}, rs[0].lower = rs[0].upper = 1.5, rs[0].k = 20;
  rs[1].kind = Restraint::Kind::Angle, rs[1].atoms = {0, 1, 2, 0}, rs[1].lower = rs[1].upper = 100, rs[1].k = 6;
  rs[2].kind = Restraint::Kind::Torsion, rs[2].atoms = {0, 1, 2, 3}, rs[2].lower = rs[2].upper = -60, rs[2].k = 4;
  const auto dir = std::filesystem::temp_directory_path() / "caps_restraints_lmp";
  std::filesystem::remove_all(dir);
  std::filesystem::create_directories(dir);
  EnergyOptions e;
  e.cutoff = 10, e.coulomb = false;
  LammpsStyle st;
  st.restraints = rs;
  write_lammps_data_ff(s, ff, e, (dir / "r.data").string(), false, st);
  LammpsRun run;
  run.kind = LammpsRun::Kind::None;
  write_lammps_input(s, ff, e, "r.data", (dir / "r.in").string(), 0, true, run, st);
  {
    std::ofstream in(dir / "r.in", std::ios::app);
    in << "thermo_style custom step f_restraints\nthermo_modify format float %.10f\nrun 0\n";
  }
  const std::string cmd = "cd '" + dir.string() + "' && lmp -in r.in -log none > out.txt 2>&1";
  ASSERT_EQ(std::system(cmd.c_str()), 0) << slurp(dir / "out.txt");
  const std::string out = slurp(dir / "out.txt");
  const auto at = out.find("f_restraints");
  ASSERT_NE(at, std::string::npos) << out;
  std::istringstream ls(out.substr(out.find('\n', at) + 1));
  long step = -1;
  double lmp_e = 0;
  ls >> step >> lmp_e;
  const auto x = coords(s);
  const double dphi = (restraint_value(rs[2], x, s.cell) - rs[2].lower) * M_PI / 180;
  const double expect = energy({rs[0], rs[1]}, x, s.cell) + rs[2].k * (1 - std::cos(dphi));
  EXPECT_NEAR(lmp_e, expect, 1e-6 * std::max(1.0, expect)) << out;
  std::printf("lmp fix restrain %.8f, CAPS %.8f kcal/mol\n", lmp_e, expect);
}
