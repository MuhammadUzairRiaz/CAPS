#include <gtest/gtest.h>

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
  EXPECT_NE(all.find("pair_style lj/cut/coul/dsf"), std::string::npos);
  const System t = read_lammps_data(path);
  EXPECT_EQ(t.atoms.size(), s.atoms.size());
  EXPECT_EQ(t.bonds.size(), s.bonds.size());
  EXPECT_NEAR(t.density(), s.density(), 1e-3);
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
  write_lammps_data_ff(s, ff, EnergyOptions{}, path);
  write_lammps_input(s, ff, EnergyOptions{}, path, in);
  std::ifstream f(path);
  const std::string all((std::istreambuf_iterator<char>(f)), {});
  EXPECT_NE(all.find("angle_style hybrid harmonic class2"), std::string::npos);
  EXPECT_NE(all.find("bond_style class2"), std::string::npos);
  EXPECT_NE(all.find("\nBondBond Coeffs\n\n1 skip"), std::string::npos);   // the harmonic angle type
  EXPECT_NE(all.find("2 class2 1 1.5 1.5"), std::string::npos);
  EXPECT_EQ(all.find("BondBond13"), std::string::npos);                    // no class II dihedrals
  std::ifstream g(in);
  const std::string script((std::istreambuf_iterator<char>(g)), {});
  EXPECT_NE(script.find("read_data " + path), std::string::npos);
  EXPECT_NE(script.find("special_bonds lj 0 0 0.5"), std::string::npos);
  // terms LAMMPS cannot reproduce exactly are refused
  ForceField charmm = ff;
  charmm.lj14_types.assign(charmm.type_names.size(), {0.05, 3.0});
  EXPECT_THROW(write_lammps_data_ff(s, charmm, EnergyOptions{}, path), FieldError);
}
