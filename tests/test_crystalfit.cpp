#include <gtest/gtest.h>

#include <cmath>

#include "caps/crystalfit.hpp"
#include "caps/ffdef.hpp"
#include "caps/uff.hpp"

using namespace caps;

namespace {
// primitive fcc cell of lattice constant a with atoms at the given fractional positions
System fcc(double a, const std::vector<std::pair<int, Vec3>>& atoms) {
  System s;
  s.cell.a = {0, a / 2, a / 2}, s.cell.b = {a / 2, 0, a / 2}, s.cell.c = {a / 2, a / 2, 0};
  for (const auto& [el, f] : atoms) {
    Atom at;
    at.element = el;
    at.pos = s.cell.to_cartesian(f);
    s.atoms.push_back(at);
  }
  return s;
}
// names every atom's type by the given name for its element
void name_types(ForceField& ff, const System& s, const std::vector<std::pair<int, std::string>>& names) {
  ff.type_names.clear();
  for (const auto& [el, n] : names) ff.type_names.push_back(n);
  ff.lj.assign(names.size(), {0, 0});
  for (size_t i = 0; i < s.atoms.size(); ++i)
    for (size_t t = 0; t < names.size(); ++t)
      if (s.atoms[i].element == names[t].first) ff.type_index[i] = int(t), ff.atom_type[i] = names[t].second;
  // ions and atoms only: no bonded terms (UFF would bond Na to Cl by distance), no exclusions
  ff.bonds.clear(), ff.angles.clear(), ff.dihedrals.clear(), ff.impropers.clear(), ff.inversions.clear(), ff.bonds_x.clear(), ff.angles_x.clear();
  ff.pairs14.clear();
  for (auto& e : ff.excluded) e.clear();
}

constexpr double kEps = 0.2385, kSig = 3.405, kRc = 8.5;
ForceField argon_ff(const System& s, double eps, double sig) {
  ForceField ff = default_forcefield(s);
  name_types(ff, s, {{18, "Ar"}});
  ff.lj[0] = {eps, sig};
  ff.charge.assign(s.atoms.size(), 0.0);
  return ff;
}
EnergyOptions lj_energy() {
  EnergyOptions e;
  e.cutoff = kRc, e.coulomb = false, e.tail = false;
  return e;
}

// rigid-ion rock salt: Buckingham Na–Cl and Cl–Cl, charges ±1, Ewald sums
constexpr double kA = 26990.0, kRho = 0.3216, kAcl = 50000.0, kRhoCl = 0.33, kCcl = 1500.0;
ForceField salt_ff(const System& s, double A) {
  ForceField ff = default_forcefield(s);
  name_types(ff, s, {{11, "Na"}, {17, "Cl"}});
  ff.charge.assign(s.atoms.size(), 0.0);
  for (size_t i = 0; i < s.atoms.size(); ++i) ff.charge[i] = s.atoms[i].element == 11 ? 1.0 : -1.0;
  ff.pair_func[{0, 1}] = {1, A, kRho, 0.0};
  ff.pair_func[{1, 1}] = {1, kAcl, kRhoCl, kCcl};
  return ff;
}
EnergyOptions salt_energy() {
  EnergyOptions e;
  e.cutoff = 8.0, e.tail = false;
  e.electrostatics = EnergyOptions::Electrostatics::PME;
  e.ewald_rtol = 1e-7;
  return e;
}
System salt(double a) { return fcc(a, {{11, {0, 0, 0}}, {17, {0.5, 0.5, 0.5}}}); }

// strains positions and cell by I + ε (ε symmetric; e = xx yy zz yz xz xy)
System strained(const System& s, const std::array<double, 6>& e) {
  auto F = [&](const Vec3& v) {
    return Vec3{v[0] + e[0] * v[0] + e[5] * v[1] + e[4] * v[2], v[1] + e[5] * v[0] + e[1] * v[1] + e[3] * v[2], v[2] + e[4] * v[0] + e[3] * v[1] + e[2] * v[2]};
  };
  System t = s;
  t.cell.a = F(s.cell.a), t.cell.b = F(s.cell.b), t.cell.c = F(s.cell.c);
  for (auto& at : t.atoms) at.pos = F(at.pos);
  return t;
}
}  // namespace

// The stress (−W/V from the virial tensor, Ewald included) is the strain derivative of the energy per volume: a sheared
// rock salt with one ion off its site, stresses against central differences of the energy over ± strain
TEST(CrystalFit, StressIsTheStrainDerivative) {
  System base = strained(salt(5.5), {0, 0, 0, 0, 0, 0.02});
  base.atoms[1].pos = base.atoms[1].pos + Vec3{0.05, -0.03, 0.02};
  const double d = 1e-5;
  FitOptions o;
  o.energy = salt_energy();
  auto add = [&](const System& s) {
    FitStructure fs;
    fs.cell = s;
    fs.field_for = [](const System& x) { return salt_ff(x, kA); };
    o.structures.push_back(fs);
  };
  add(base);
  for (int k : {0, 5}) {   // xx and xy
    std::array<double, 6> e{};
    e[size_t(k)] = k < 3 ? d : d / 2;
    add(strained(base, e));
    e[size_t(k)] = -e[size_t(k)];
    add(strained(base, e));
  }
  const FitResult r = fit_to_crystals(o);
  const double V = base.cell.volume();
  const double kGPa = 4184.0 / 6.02214076e23 / 1e-30 / 1e9;
  const double sxx = (r.structures[1].energy_initial - r.structures[2].energy_initial) / (2 * d) / V * kGPa;
  const double sxy = (r.structures[3].energy_initial - r.structures[4].energy_initial) / (2 * d) / V * kGPa;
  std::printf("σxx %.6f (virial) %.6f (energy); σxy %.6f %.6f GPa; largest force %.3f kcal/mol/Å\n", r.structures[0].stress_initial[0], sxx,
              r.structures[0].stress_initial[5], sxy, r.structures[0].force_initial);
  EXPECT_NEAR(r.structures[0].stress_initial[0], sxx, 1e-4 * std::fabs(sxx) + 1e-5);
  EXPECT_NEAR(r.structures[0].stress_initial[5], sxy, 1e-4 * std::fabs(sxy) + 1e-5);
  EXPECT_GT(std::fabs(sxy), 0.05);
  EXPECT_GT(r.structures[0].force_initial, 1.0);
}

// Relaxed Lennard-Jones fcc: the volume and C11, C12, C44 equal the lattice sums at the analytic minimum (central
// pair forces at zero stress: C12 = C44, the Cauchy relation)
TEST(CrystalFit, LennardJonesFccIsExact) {
  // analytic: E(a) per atom, its minimum, then C_ijkl = (1/2Ω) Σ (φ'' − φ'/r) R_i R_j R_k R_l / r²
  auto lattice = [](double a, auto&& visit) {
    const System s = fcc(a, {});
    for (int i = -8; i <= 8; ++i)
      for (int j = -8; j <= 8; ++j)
        for (int l = -8; l <= 8; ++l) {
          if (!i && !j && !l) continue;
          const Vec3 R = s.cell.a * double(i) + s.cell.b * double(j) + s.cell.c * double(l);
          if (norm(R) < kRc) visit(R);
        }
  };
  auto dEda = [&](double a) {
    double g = 0;
    lattice(a, [&](const Vec3& R) {
      const double r = norm(R), sr6 = std::pow(kSig / r, 6);
      g += 0.5 * 4 * kEps * (-12 * sr6 * sr6 + 6 * sr6) / r * (r / a);
    });
    return g;
  };
  double lo = 5.0, hi = 5.6;
  for (int k = 0; k < 200; ++k) {
    const double mid = 0.5 * (lo + hi);
    (dEda(mid) > 0 ? hi : lo) = mid;
  }
  const double a0 = 0.5 * (lo + hi), omega = a0 * a0 * a0 / 4;
  double c11 = 0, c12 = 0, c44 = 0;
  lattice(a0, [&](const Vec3& R) {
    const double r = norm(R), sr6 = std::pow(kSig / r, 6);
    const double d1 = 4 * kEps * (-12 * sr6 * sr6 + 6 * sr6) / r, d2 = 4 * kEps * (156 * sr6 * sr6 - 42 * sr6) / (r * r);
    const double w = (d2 - d1 / r) / (r * r) / (2 * omega);
    c11 += w * R[0] * R[0] * R[0] * R[0], c12 += w * R[0] * R[0] * R[1] * R[1], c44 += w * R[1] * R[2] * R[1] * R[2];
  });
  const double kGPa = 4184.0 / 6.02214076e23 / 1e-30 / 1e9;
  c11 *= kGPa, c12 *= kGPa, c44 *= kGPa;

  FitOptions o;
  o.energy = lj_energy();
  o.relaxed = true;
  FitStructure fs;
  fs.cell = fcc(5.4, {{18, {0, 0, 0}}});
  fs.field_for = [](const System& x) { return argon_ff(x, kEps, kSig); };
  o.structures.push_back(fs);
  for (const char* k : {"volume", "C11", "C12", "C44", "bulk"}) o.targets.push_back({0, k, 0, 1});
  const FitResult r = fit_to_crystals(o);
  const auto& ob = r.observations;
  std::printf("a %.6f (analytic %.6f) · C11 %.5f (%.5f) C12 %.5f (%.5f) C44 %.5f (%.5f) GPa\n", std::cbrt(4 * ob[0].final), a0, ob[1].final, c11,
              ob[2].final, c12, ob[3].final, c44);
  EXPECT_NEAR(ob[0].final, omega, 1e-6 * omega);
  EXPECT_NEAR(ob[1].final, c11, 1e-4 * c11);
  EXPECT_NEAR(ob[2].final, c12, 1e-4 * c11);
  EXPECT_NEAR(ob[3].final, c44, 1e-4 * c11);
  EXPECT_NEAR(ob[4].final, (c11 + 2 * c12) / 3, 1e-4 * c11);
  EXPECT_NEAR(r.structures[0].final[3], 60.0, 1e-6);   // the primitive cell keeps its angles
}

// Relaxed fitting recovers the parameters that made the data: ε and σ of argon from its cell and C11
TEST(CrystalFit, RelaxedFitRecoversLennardJones) {
  FitOptions truth;
  truth.energy = lj_energy();
  truth.relaxed = true;
  FitStructure fs;
  fs.cell = fcc(5.3, {{18, {0, 0, 0}}});
  fs.field_for = [](const System& x) { return argon_ff(x, kEps, kSig); };
  truth.structures.push_back(fs);
  truth.targets.push_back({0, "C11", 0, 1});
  const FitResult t = fit_to_crystals(truth);
  const double c11 = t.observations[0].final;

  FitOptions o;
  o.energy = lj_energy();
  o.relaxed = true;
  fs.cell.cell = t.structures[0].cell;
  fs.cell.atoms[0].pos = t.structures[0].positions[0];
  fs.field_for = [](const System& x) { return argon_ff(x, 1.15 * kEps, 0.96 * kSig); };
  o.structures.push_back(fs);
  o.targets.push_back({0, "C11", c11, 1});
  o.parameters.push_back({"lj", {"Ar"}, "epsilon"});
  o.parameters.push_back({"lj", {"Ar"}, "sigma"});
  o.parameters.back().lower = 3.0;
  const FitResult r = fit_to_crystals(o);
  std::printf("ε %.7f σ %.7f after %d iterations (%s), S %.3g → %.3g\n", r.final[0], r.final[1], r.iterations, r.stopped_by.c_str(), r.ss_initial, r.ss_final);
  EXPECT_NEAR(r.final[0], kEps, 1e-5 * kEps);
  EXPECT_NEAR(r.final[1], kSig, 1e-5 * kSig);
  ASSERT_EQ(r.terms.size(), 1u);
  EXPECT_EQ(r.terms[0].kind, "lj");
  EXPECT_NEAR(r.terms[0].values[1], r.final[1], 1e-12);
}

// Conventional fitting recovers the Buckingham A that holds the observed rock salt at zero stress (Ewald sums)
TEST(CrystalFit, ConventionalFitRecoversBuckingham) {
  FitOptions truth;
  truth.energy = salt_energy();
  truth.relaxed = true;
  FitStructure fs;
  fs.cell = salt(5.6);
  fs.field_for = [](const System& x) { return salt_ff(x, kA); };
  truth.structures.push_back(fs);
  const FitResult t = fit_to_crystals(truth);
  System obs = fs.cell;
  obs.cell = t.structures[0].cell;
  for (size_t i = 0; i < obs.atoms.size(); ++i) obs.atoms[i].pos = t.structures[0].positions[i];

  FitOptions o;
  o.energy = salt_energy();
  fs.cell = obs;
  fs.field_for = [](const System& x) { return salt_ff(x, 0.8 * kA); };
  o.structures.push_back(fs);
  o.parameters.push_back({"pair", {"Na", "Cl"}, "A"});
  const FitResult r = fit_to_crystals(o);
  std::printf("a %.5f Å · A %.4f (true %.4f) · stress %.3g → %.3g GPa\n", std::sqrt(2.0) * t.structures[0].final[0], r.final[0], kA,
              r.structures[0].stress_initial[0], r.structures[0].stress_final[0]);
  EXPECT_NEAR(r.final[0], kA, 1e-5 * kA);
  EXPECT_LT(std::fabs(r.structures[0].stress_final[0]), 1e-5);
  ASSERT_EQ(r.terms.size(), 1u);
  EXPECT_EQ(r.terms[0].style, "buck");
  EXPECT_NEAR(r.terms[0].values[1], kRho, 1e-12);
}

// A balanced charge: the balancing type takes what keeps the cell neutral; the terms write back into a definition
TEST(CrystalFit, BalancedChargeAndWriteBack) {
  FitOptions o;
  o.energy = salt_energy();
  o.iterations = 0;
  FitStructure fs;
  fs.cell = salt(5.6);
  fs.field_for = [](const System& x) { return salt_ff(x, kA); };
  o.structures.push_back(fs);
  FitParameter q{"charge", {"Na"}, "q"};
  q.lower = q.upper = 0.9;
  q.balance = "Cl";
  o.parameters.push_back(q);
  o.parameters.push_back({"pair", {"Cl", "Cl"}, "C"});
  o.parameters.push_back({"pair", {"Na", "Na"}, "sigma"});
  const FitResult r = fit_to_crystals(o);
  ASSERT_EQ(r.terms.size(), 4u);
  EXPECT_NEAR(r.terms[0].values[0], 0.9, 1e-12);
  EXPECT_EQ(r.terms[1].types[0], "Cl");
  EXPECT_NEAR(r.terms[1].values[0], -0.9, 1e-12);

  FFDef def;
  for (const char* n : {"Na", "Cl"}) {
    FFType t;
    t.name = n;
    def.types.push_back(t);
  }
  FFRule cl;
  cl.name = "Cl-Cl", cl.match = {"Cl", "Cl"}, cl.style = "buck", cl.params = {1, 2, 3};
  def.pairs.push_back(cl);
  const auto notes = apply_fitted_terms(def, r.terms);
  EXPECT_NEAR(def.types[0].charge, 0.9, 1e-12);
  EXPECT_NEAR(def.types[1].charge, -0.9, 1e-12);
  ASSERT_EQ(def.pairs.size(), 2u);   // Cl-Cl updated, Na-Na appended (a Lennard-Jones pair)
  EXPECT_NEAR(def.pairs[0].params[2], kCcl, 1e-12);
  EXPECT_EQ(def.pairs[1].match, (std::vector<std::string>{"Na", "Na"}));
  EXPECT_TRUE(def.pairs[1].style.empty());
}
