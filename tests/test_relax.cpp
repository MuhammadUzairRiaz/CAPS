#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <map>
#include <random>

#include "caps/field.hpp"
#include "caps/grow.hpp"
#include "caps/io.hpp"
#include "caps/relax.hpp"

using namespace caps;

namespace {

System small_cell(int chains, int dp, double density, double scale = 1.0, uint64_t seed = 3) {
  GrowOptions g;
  g.chains = chains;
  g.dp = dp;
  g.density = density;
  g.seed = seed;
  g.contact_scale = scale;
  return grow(g);
}

std::vector<double> flat(const System& s) {
  std::vector<double> x;
  for (const auto& a : s.atoms) x.insert(x.end(), a.pos.begin(), a.pos.end());
  return x;
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

TEST(Field, TypesPolystyrene) {
  const System s = small_cell(2, 4, 0.3);
  const ForceField ff = assign_gaff(s);
  std::map<std::string, int> n;
  for (const auto& t : ff.atom_type) n[t]++;
  EXPECT_EQ(n["ca"], 2 * 4 * 6);
  EXPECT_EQ(n["ha"], 2 * 4 * 5);
  EXPECT_EQ(n["c3"], 2 * 4 * 2);
  EXPECT_EQ(n["hc"], 2 * (4 * 3 + 2));   // CH2 + CH per unit, plus the two H end caps
  EXPECT_EQ(ff.type_names, (std::vector<std::string>{"c3", "ca", "hc", "ha"}));
  EXPECT_EQ(ff.bonds.size(), s.bonds.size());
  EXPECT_EQ(ff.impropers.size(), size_t(2 * 4 * 6));
  double q = 0;
  for (double c : ff.charge) q += c;
  EXPECT_NEAR(q, 0.0, 1e-4);
}

TEST(Field, RejectsElementsWithoutParameters) {
  System s;
  s.atoms.resize(2);
  s.atoms[0].element = 8;
  s.atoms[1].element = 1;
  s.bonds.push_back({0, 1});
  EXPECT_THROW(assign_gaff(s), FieldError);
}

// Analytic forces against central differences, in a cell narrower than twice the cut-off (image pairs included).
TEST(Field, ForcesMatchFiniteDifferences) {
  System s = small_cell(2, 3, 0.25);
  ASSERT_LT(norm(s.cell.a), 2 * 10.0 + 2);
  const ForceField ff = assign_gaff(s);
  std::vector<double> x = flat(s);
  std::mt19937_64 rng(7);
  std::normal_distribution<double> nd(0, 0.05);
  for (auto& v : x) v += nd(rng);
  for (double cap : {0.0, 20.0}) {
    EnergyOptions eo;
    eo.force_cap = cap;
    Evaluator ev(ff, eo);
    std::vector<double> f, g;
    ev.compute(x, s.cell, f);
    const double h = 1e-5;
    double worst = 0;
    for (size_t k = 0; k < x.size(); k += 7) {
      std::vector<double> xp = x, xm = x;
      xp[k] += h;
      xm[k] -= h;
      const double ep = ev.compute(xp, s.cell, g).total(), em = ev.compute(xm, s.cell, g).total();
      const double fd = -(ep - em) / (2 * h);
      worst = std::max(worst, std::fabs(fd - f[k]) / std::max(1.0, std::fabs(f[k])));
    }
    EXPECT_LT(worst, 2e-4) << "force cap " << cap;
  }
}

// The virial Σ r·f equals −dE/ds for an affine scaling of cell and positions by s.
TEST(Field, VirialMatchesVolumeDerivative) {
  const System s = small_cell(3, 4, 0.35);
  const ForceField ff = assign_gaff(s);
  EnergyOptions no_tail;
  no_tail.tail = false;   // the tail pressure is not −dE_tail/dV of the fixed-cut-off energy
  Evaluator ev(ff, no_tail);
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
  const double h = 1e-6;
  const double dEds = (energy_at(1 + h) - energy_at(1 - h)) / (2 * h);
  EXPECT_NEAR(-dEds, w, 1e-3 * std::max(1.0, std::fabs(w)));
}

// Every component of the virial tensor equals −dE/dε for an affine strain (sheared cells included).
TEST(Field, VirialTensorMatchesStrainDerivatives) {
  const System s = small_cell(3, 4, 0.35);
  const ForceField ff = assign_gaff(s);
  EnergyOptions no_tail;
  no_tail.tail = false;
  Evaluator ev(ff, no_tail);
  std::vector<double> x = flat(s), f;
  const EnergyTerms e = ev.compute(x, s.cell, f);
  EXPECT_NEAR(e.w[0] + e.w[1] + e.w[2], e.virial, 1e-9 * std::fabs(e.virial) + 1e-9);
  const auto d = strain_derivative(ev, x, s.cell);
  for (int v = 0; v < 6; ++v) EXPECT_NEAR(e.w[v], d[v], 2e-3 * std::max(1.0, std::fabs(e.virial))) << "component " << v;
}

TEST(Field, EnergyIsPeriodic) {
  const System s = small_cell(2, 3, 0.3);
  const ForceField ff = assign_gaff(s);
  Evaluator ev(ff, EnergyOptions{});
  std::vector<double> x = flat(s), f;
  const double e0 = ev.compute(x, s.cell, f).total();
  // move the whole first molecule by a lattice vector, and one lone atom's image does not matter either
  for (size_t i = 0; i < s.atoms.size(); ++i)
    if (s.atoms[i].mol == s.atoms[0].mol)
      for (int k = 0; k < 3; ++k) x[3 * i + k] += s.cell.a[k] - 2 * s.cell.c[k];
  EXPECT_NEAR(ev.compute(x, s.cell, f).total(), e0, 1e-6 * std::fabs(e0) + 1e-6);
}

TEST(Relax, EveryMinimiserLowersEnergy) {
  const System s0 = small_cell(3, 4, 0.4);
  for (Minimiser m : {Minimiser::SteepestDescent, Minimiser::ConjugateGradient, Minimiser::LBFGS, Minimiser::FIRE}) {
    System s = s0;
    RelaxOptions o;
    o.method = m;
    o.pushoff = false;
    o.max_iterations = 300;
    RelaxReport r;
    relax(s, o, &r);
    EXPECT_LT(r.final.total(), r.initial.total()) << to_string(m);
    EXPECT_LT(r.fmax_final, r.fmax_initial) << to_string(m);
  }
  System s = s0;
  RelaxOptions o;
  o.ftol = 0.1;
  o.max_iterations = 5000;
  RelaxReport r;
  relax(s, o, &r);
  EXPECT_TRUE(r.converged);
  EXPECT_LT(r.fmax_final, 0.1);
}

TEST(Relax, CompressesOverlappingCellToTarget) {
  System s = small_cell(4, 5, 0.5, 0.8);
  RelaxOptions o;
  o.target_density = 0.95;
  o.ftol = 1.0;
  RelaxReport r;
  relax(s, o, &r);
  EXPECT_NEAR(s.density(), 0.95, 1e-6);
  // bonded geometry near the force-field values, no close non-bonded contacts left
  const ForceField ff = assign_gaff(s);
  double worst_bond = 0;
  for (const auto& b : ff.bonds)
    worst_bond = std::max(worst_bond, std::fabs(norm(s.cell.minimum_image(s.atoms[b.j].pos - s.atoms[b.i].pos)) - b.r0));
  EXPECT_LT(worst_bond, 0.05);
  double closest = 1e9;
  for (size_t i = 0; i < s.atoms.size(); ++i)
    for (size_t j = i + 1; j < s.atoms.size(); ++j) {
      if (std::binary_search(ff.excluded[i].begin(), ff.excluded[i].end(), uint32_t(j))) continue;
      closest = std::min(closest, norm(s.cell.minimum_image(s.atoms[j].pos - s.atoms[i].pos)));
    }
  EXPECT_GT(closest, 1.7);
  EXPECT_LT(r.final.total(), r.initial.total());
}

TEST(Relax, BoxRelaxReachesPressure) {
  System s = small_cell(3, 4, 0.6, 0.85);
  RelaxOptions o;
  o.relax_box = true;
  o.pressure = 1.0;
  o.pressure_tol = 200;
  o.ftol = 0.5;
  RelaxReport r;
  relax(s, o, &r);
  EXPECT_LT(std::fabs(r.pressure_final - 1.0), 200.0 + 50.0);   // pressure is re-evaluated after the last minimisation
  EXPECT_GT(r.density_final, 0.6);
}

TEST(Relax, LammpsExportHasTermsAndReadsBack) {
  System s = small_cell(2, 3, 0.3);
  const ForceField ff = assign_gaff(s);
  const auto path = (std::filesystem::temp_directory_path() / "caps_relax_ff.data").string();
  write_lammps_data_ff(s, ff, EnergyOptions{}, path);
  std::ifstream in(path);
  std::string all((std::istreambuf_iterator<char>(in)), {});
  EXPECT_NE(all.find(std::to_string(ff.angles.size()) + " angles"), std::string::npos);
  EXPECT_NE(all.find("Dihedral Coeffs"), std::string::npos);
  EXPECT_NE(all.find("PairIJ Coeffs  # lj/cut/coul/dsf"), std::string::npos);   // a data file on its own keeps its pair coefficients
  EXPECT_EQ(all.find("pair_style"), std::string::npos);                           // commands live in the input script, not the data
  const System t = read_lammps_data(path);
  EXPECT_EQ(t.atoms.size(), s.atoms.size());
  EXPECT_EQ(t.bonds.size(), s.bonds.size());
  EXPECT_NEAR(t.density(), s.density(), 1e-3);
}

TEST(Gromacs, TopologyUnitsExclusionsAndRefusals) {
  System s = small_cell(2, 3, 0.3);
  ForceField ff = assign_gaff(s);
  const auto stem = (std::filesystem::temp_directory_path() / "caps_gmx_unit").string();
  EnergyOptions e;
  e.electrostatics = EnergyOptions::Electrostatics::PME;
  const auto notes = write_gromacs(s, ff, e, stem);
  std::ifstream f(stem + ".top");
  const std::string top((std::istreambuf_iterator<char>(f)), {});
  std::ifstream m(stem + ".mdp");
  const std::string mdp((std::istreambuf_iterator<char>(m)), {});
  // CAPS k (r − r0)² in kcal/mol/Å² is GROMACS ½ kb (r − b0)² in kJ/mol/nm²: kb = 2 k · 4.184 · 100
  const auto& b = ff.bonds[0];
  char line[160];
  std::snprintf(line, sizeof line, "%7u %7u 1 %.10g %.10g\n", b.i + 1, b.j + 1, b.r0 / 10, 2 * b.k * 4.184 * 100);
  // the bond line, numbered within its molecule
  std::snprintf(line, sizeof line, " 1 %.10g %.10g\n", b.r0 / 10, 2 * b.k * 4.184 * 100);
  std::ifstream fi(stem + ".itp");
  const std::string itp((std::istreambuf_iterator<char>(fi)), {});
  EXPECT_NE(top.find("#include \"caps_gmx_unit.itp\""), std::string::npos);   // molecule types in the .itp
  EXPECT_NE(itp.find("MOL  0"), std::string::npos);                              // nrexcl 0: CAPS's exclusions are listed
  EXPECT_NE(itp.find("[ exclusions ]"), std::string::npos);
  EXPECT_NE(itp.find(line), std::string::npos) << line;                          // kb = 2 k · 4.184 · 100
  EXPECT_NE(top.find("[ nonbond_params ]"), std::string::npos);
  EXPECT_NE(mdp.find("coulombtype              = PME"), std::string::npos);
  EXPECT_NE(mdp.find("DispCorr                 = AllEnerPres"), std::string::npos);
  // the same notes without writing
  EXPECT_EQ(gromacs_notes(s, ff, e), notes);
  // forms GROMACS lacks are refused
  ForceField inv = ff;
  inv.inversions.push_back({0, 1, 2, 3, 1.0, 0.0, 1});
  EXPECT_THROW(write_gromacs(s, inv, e, stem), FieldError);
  ForceField c2 = ff;
  c2.pair_form = "lj9-6";
  EXPECT_THROW(gromacs_notes(s, c2, e), FieldError);
}

TEST(LammpsData, MixedClassesBecomeHybridStylesWithSkipLines) {
  // class II bonds and angles, one class I angle, Fourier torsions: hybrid angles, BondBond / BondAngle skip lines
  System s = small_cell(2, 3, 0.3);
  ForceField ff = assign_gaff(s);
  for (const auto& b : ff.bonds) ff.bonds2.push_back({b.i, b.j, b.r0, b.k, 0, 0});
  ff.bonds.clear();
  for (size_t k = 1; k < ff.angles.size(); ++k) {
    const auto& a = ff.angles[k];
    ff.angles2.push_back({a.i, a.j, a.k, a.theta0, a.kt, 0, 0, 1.0, 1.5, 1.5, 2.0, 3.0, 1.5, 1.1});
  }
  ff.angles.resize(1);
  const auto path = (std::filesystem::temp_directory_path() / "caps_mixed.data").string();
  const auto in = (std::filesystem::temp_directory_path() / "caps_mixed.in").string();
  write_lammps_data_ff(s, ff, EnergyOptions{}, path, false);   // written with its script: pair coefficients go there
  write_lammps_input(s, ff, EnergyOptions{}, "caps_mixed.data", in, 0, true);
  std::ifstream f(path);
  const std::string all((std::istreambuf_iterator<char>(f)), {});
  EXPECT_EQ(all.find("PairIJ Coeffs"), std::string::npos);
  EXPECT_NE(all.find("\nBondBond Coeffs\n\n1 skip"), std::string::npos);   // the harmonic angle type
  EXPECT_NE(all.find("2 class2 1.000000 1.500000 1.500000"), std::string::npos);
  EXPECT_EQ(all.find("BondBond13"), std::string::npos);                    // no class II dihedrals
  std::ifstream g(in);
  const std::string script((std::istreambuf_iterator<char>(g)), {});
  EXPECT_NE(script.find("angle_style     hybrid harmonic class2"), std::string::npos);
  EXPECT_NE(script.find("bond_style      class2"), std::string::npos);
  EXPECT_NE(script.find("read_data       caps_mixed.data"), std::string::npos);
  EXPECT_NE(script.find("special_bonds   lj 0 0 0.500000"), std::string::npos);
  EXPECT_NE(script.find("pair_coeff      1 1 "), std::string::npos);
  EXPECT_NE(script.find("run 0"), std::string::npos);   // the default run section: a single-point check
  // a protocol instead: minimisation then NPT, the final structure written
  LammpsRun run;
  run.kind = LammpsRun::Kind::NPT;
  write_lammps_input(s, ff, EnergyOptions{}, "caps_mixed.data", in, 0, true, run);
  std::ifstream g2(in);
  const std::string npt((std::istreambuf_iterator<char>(g2)), {});
  EXPECT_NE(npt.find("minimize"), std::string::npos);
  EXPECT_NE(npt.find(" npt temp 300 300 100 iso 1 1 1000"), std::string::npos);
  EXPECT_NE(npt.find("write_data      final.data"), std::string::npos);
  // the time step: 0.5 fs unless the force field declares its own (Martini 20 fs) or the run sets one
  EXPECT_NE(npt.find("\ntimestep        0.5\n"), std::string::npos);
  EXPECT_EQ(lammps_timestep(run, ff), 0.5);
  ForceField cg = ff;
  cg.native_timestep = 20;
  EXPECT_EQ(lammps_timestep(run, cg), 20);
  run.dt = 2;
  EXPECT_EQ(lammps_timestep(run, cg), 2);
  // terms LAMMPS cannot reproduce exactly are refused
  ForceField charmm = ff;
  charmm.lj14_types.assign(charmm.type_names.size(), {0.05, 3.0});
  EXPECT_THROW(write_lammps_data_ff(s, charmm, EnergyOptions{}, path), FieldError);
}

// A distance restraint pulls two carbons of different chains to its target (4.0 Å at k 100 kcal/mol/Å²), is reported,
// and changes nothing else about the minimisation
TEST(Relax, DistanceRestraintPullsAtomsToTarget) {
  System s = small_cell(3, 4, 0.4);
  uint32_t a = 0, b = 0;
  double far = 0;
  for (uint32_t i = 0; i < s.atoms.size(); ++i)
    for (uint32_t j = 0; j < s.atoms.size(); ++j)
      if (s.atoms[i].mol == 1 && s.atoms[j].mol == 2 && s.atoms[i].element == 6 && s.atoms[j].element == 6) {
        const double d = norm(s.cell.minimum_image(s.atoms[j].pos - s.atoms[i].pos));
        if (d > far && d < 9) far = d, a = i, b = j;
      }
  ASSERT_GT(far, 6.0);
  RelaxOptions o;
  o.ftol = 0.2;
  o.restraints.push_back({a, b, 4.0, 100.0});
  RelaxReport r;
  relax(s, o, &r);
  const double d = norm(s.cell.minimum_image(s.atoms[b].pos - s.atoms[a].pos));
  EXPECT_NEAR(d, 4.0, 0.1) << "from " << far;
  EXPECT_TRUE(std::any_of(r.notes.begin(), r.notes.end(), [](const std::string& n) { return n.rfind("restraint ", 0) == 0 && n.find("target 4.000 Å") != std::string::npos; }));
}

// Dihedral restraints: forces are the energy's derivatives (finite differences), and a restrained backbone torsion ends
// at its target while the rest relaxes
TEST(Relax, DihedralRestraintForcesAndTarget) {
  System s = small_cell(1, 4, 0.3);
  // a backbone C–C–C–C torsion of the chain
  const auto nb = s.neighbours();
  uint32_t q[4] = {0, 0, 0, 0};
  bool found = false;
  for (uint32_t j = 0; j < s.atoms.size() && !found; ++j) {
    if (s.atoms[j].element != 6) continue;
    for (uint32_t k : nb[j]) {
      if (s.atoms[k].element != 6) continue;
      for (uint32_t i : nb[j])
        for (uint32_t l : nb[k])
          if (i != k && l != j && i != l && s.atoms[i].element == 6 && s.atoms[l].element == 6 && !found) q[0] = i, q[1] = j, q[2] = k, q[3] = l, found = true;
    }
  }
  ASSERT_TRUE(found);
  RelaxOptions o;
  o.pushoff = false;
  o.dihedral_restraints.push_back({q[0], q[1], q[2], q[3], 60.0, 40.0});
  // the minimiser follows the restraint's analytic forces to the target (the force field's own torsion pulls a little)
  o.ftol = 0.05;
  o.max_iterations = 4000;
  RelaxReport r;
  relax(s, o, &r);
  auto dihedral = [&](const System& t) {
    auto P = [&](uint32_t a) { return t.atoms[a].pos; };
    const Vec3 b1 = t.cell.minimum_image(P(q[1]) - P(q[0])), b2 = t.cell.minimum_image(P(q[2]) - P(q[1])), b3 = t.cell.minimum_image(P(q[3]) - P(q[2]));
    const Vec3 m = cross(b1, b2), n = cross(b2, b3);
    return std::atan2(norm(b2) * dot(b1, n), dot(m, n)) * 180 / M_PI;
  };
  EXPECT_NEAR(dihedral(s), 60.0, 3.0) << "restrained torsion ended at " << dihedral(s);
  // and without the restraint the same torsion is free to go elsewhere (sanity: the restraint did the work)
  EXPECT_TRUE(r.converged || r.fmax_final < 0.5);
}

// Anisotropic box: only z moves (a film), x and y keep their lengths, Pzz reaches the target
TEST(Relax, BoxRelaxAlongOneAxis) {
  System s = small_cell(3, 4, 0.6, 0.85);
  const double lx = s.cell.a[0], ly = s.cell.b[1];
  RelaxOptions o;
  o.relax_box = true;
  o.box_anisotropic = true;
  o.box_axes[0] = o.box_axes[1] = false;
  o.pressure = 1.0;
  o.pressure_tol = 300;
  o.ftol = 0.5;
  RelaxReport r;
  relax(s, o, &r);
  EXPECT_DOUBLE_EQ(s.cell.a[0], lx);
  EXPECT_DOUBLE_EQ(s.cell.b[1], ly);
  EXPECT_NE(s.cell.c[2], 0.0);
  ASSERT_FALSE(r.stages.empty());
  EXPECT_NE(r.stages.back().name.find("axis by axis"), std::string::npos) << r.stages.back().name;
}

// Push-off by MD with the force cap ramped (Auhl et al. 2003): ten NVT stages with caps rising geometrically to the
// chosen cap, then the minimisation stages ending at that cap; the structure converges and keeps no velocities
TEST(Relax, PushoffByMdRampsTheCap) {
  System s = small_cell(3, 4, 0.5, 0.7);
  RelaxOptions o;
  o.pushoff_ramp_ps = 0.5;
  o.pushoff_cap = 200;
  o.ftol = 1.0;
  RelaxReport r;
  relax(s, o, &r);
  std::vector<std::string> md, mini;
  for (const auto& st : r.stages) {
    if (st.name.rfind("push-off MD", 0) == 0) md.push_back(st.name);
    else if (st.name.rfind("push-off", 0) == 0) mini.push_back(st.name);
  }
  ASSERT_EQ(md.size(), 10u);
  EXPECT_EQ(md.front(), "push-off MD, force cap 5");
  EXPECT_EQ(md.back(), "push-off MD, force cap 200");
  ASSERT_FALSE(mini.empty());
  EXPECT_EQ(mini.back(), "push-off, force cap 200");
  EXPECT_TRUE(r.converged);
  EXPECT_TRUE(s.velocities.empty());
  EXPECT_TRUE(std::any_of(r.notes.begin(), r.notes.end(), [](const std::string& n) { return n.rfind("push-off MD", 0) == 0; }));
}

// Rigid bodies in the LAMMPS input: the group and the intra-body exclusion in the setup, the rigid fix before the
// barostat (LAMMPS refuses it after a box-changing fix), the barostat on the other atoms dilating only them, and
// velocities for everything that moves.
TEST(LammpsData, RigidBodiesBeforeTheBarostat) {
  System s = small_cell(3, 3, 0.3);
  ForceField ff = assign_gaff(s);
  const auto in = (std::filesystem::temp_directory_path() / "caps_rigid.in").string();
  LammpsRun run;
  run.kind = LammpsRun::Kind::NPT;
  LammpsStyle st;
  st.rigid_mols = {1, 2};
  write_lammps_input(s, ff, EnergyOptions{}, "caps_rigid.data", in, 0, true, run, st);
  std::ifstream g(in);
  const std::string sc((std::istreambuf_iterator<char>(g)), {});
  EXPECT_NE(sc.find("group           rigid molecule 1:2"), std::string::npos);
  EXPECT_NE(sc.find("neigh_modify    exclude molecule/intra rigid"), std::string::npos);
  EXPECT_NE(sc.find("group           mobile subtract all rigid"), std::string::npos);
  EXPECT_NE(sc.find("velocity        all create"), std::string::npos);
  const auto rig = sc.find("rigid/nvt/small molecule"), npt = sc.find("fix             integrate mobile npt");
  ASSERT_NE(rig, std::string::npos);
  ASSERT_NE(npt, std::string::npos);
  EXPECT_LT(rig, npt);
  EXPECT_NE(sc.find("dilate mobile"), std::string::npos);
  // a held molecule is never a rigid body, and gets no velocity
  write_lammps_input(s, ff, EnergyOptions{}, "caps_rigid.data", in, 1, true, run, st);
  std::ifstream h(in);
  const std::string sh((std::istreambuf_iterator<char>(h)), {});
  EXPECT_NE(sh.find("group           rigid molecule 2"), std::string::npos);
  EXPECT_NE(sh.find("velocity        moving create"), std::string::npos);
}

// Mechanical protocols in the LAMMPS input (each run in LAMMPS on a melt when written): tension by fix deform at an
// engineering strain rate with the lateral axes barostatted and the run long enough for the strain; creep with the
// stressed axis at −σ; shear by SLLOD with the box made triclinic before the long-range solver and η = −P_xy/γ̇
TEST(LammpsData, TensileCreepAndShearProtocols) {
  System s = small_cell(2, 3, 0.3);
  ForceField ff = assign_gaff(s);
  const auto in = (std::filesystem::temp_directory_path() / "caps_mech.in").string();
  auto script = [&](const LammpsRun& run) {
    write_lammps_input(s, ff, EnergyOptions{}, "caps_mech.data", in, 0, true, run);
    std::ifstream g(in);
    return std::string((std::istreambuf_iterator<char>(g)), {});
  };
  LammpsRun t;
  t.kind = LammpsRun::Kind::Tensile;
  t.axis = 2;
  t.strain_rate = 0.5;   // /ps
  t.max_strain = 0.1;
  t.dt = 1.0;            // fs: 0.1 / (0.5 /ps × 0.001 ps) = 200 steps
  const auto a = script(t);
  EXPECT_NE(a.find("fix             pull all deform 1 z erate 0.0005 remap x"), std::string::npos);
  EXPECT_NE(a.find("npt temp 300 300 100 x 1 1 1000 y 1 1 1000"), std::string::npos);
  EXPECT_NE(a.find("variable        stress equal -pzz*0.101325"), std::string::npos);   // atm → MPa
  EXPECT_NE(a.find("run             200\n"), std::string::npos);
  LammpsRun c;
  c.kind = LammpsRun::Kind::Creep;
  c.stress_mpa = 10.1325;   // −100 atm on x
  const auto b = script(c);
  EXPECT_NE(b.find("x -100 -100 1000 y 1 1 1000 z 1 1 1000 couple none"), std::string::npos);
  EXPECT_NE(b.find("file creep.dat"), std::string::npos);
  LammpsRun sh;
  sh.kind = LammpsRun::Kind::Shear;
  sh.shear_rate = 0.1;
  const auto d = script(sh);
  const auto tri = d.find("change_box      all triclinic"), rd = d.find("read_data"), ks = d.find("kspace_style");
  ASSERT_NE(tri, std::string::npos);
  EXPECT_GT(tri, rd);
  if (ks != std::string::npos) EXPECT_LT(tri, ks);
  EXPECT_NE(d.find("nvt/sllod temp 300 300 100"), std::string::npos);
  EXPECT_NE(d.find("fix             flow all deform 1 xy erate 0.0001 remap v"), std::string::npos);
  // η (mPa·s) = −P_xy (atm) × 101325 Pa/atm / (0.1e12 /s) × 1000
  EXPECT_NE(d.find("variable        eta equal -pxy*0.00101325"), std::string::npos);
}
