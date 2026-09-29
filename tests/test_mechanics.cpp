#include <gtest/gtest.h>

#include <cmath>
#include <memory>
#include <random>
#include <vector>

#include "caps/mechanics.hpp"
#include "caps/grow.hpp"
#include "caps/relax.hpp"

using namespace caps;

namespace {

constexpr double kGPa = 4184.0 / 6.02214076e23 / 1e-30 / 1e9;
constexpr int VA[6] = {0, 1, 2, 1, 0, 0}, VB[6] = {0, 1, 2, 2, 2, 1};

// Argon fcc crystal, n × n × n conventional cells of edge a; one LJ type, no bonds.
System argon(int n, double a) {
  System s;
  const double basis[4][3] = {{0, 0, 0}, {0.5, 0.5, 0}, {0.5, 0, 0.5}, {0, 0.5, 0.5}};
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      for (int k = 0; k < n; ++k)
        for (const auto& b : basis) {
          Atom at;
          at.element = 18;
          at.pos = {(i + b[0]) * a, (j + b[1]) * a, (k + b[2]) * a};
          s.atoms.push_back(at);
        }
  s.cell.a = {n * a, 0, 0};
  s.cell.b = {0, n * a, 0};
  s.cell.c = {0, 0, n * a};
  return s;
}

std::shared_ptr<ForceField> argon_ff(size_t n) {
  auto ff = std::make_shared<ForceField>();
  ff->name = "argon";
  ff->type_names = {"Ar"};
  ff->lj = {{0.238, 3.405}};
  ff->atom_type.assign(n, "Ar");
  ff->why.assign(n, "");
  ff->type_index.assign(n, 0);
  ff->charge.assign(n, 0.0);
  ff->mass.assign(n, 39.948);
  ff->excluded.assign(n, {});
  return ff;
}

EnergyOptions lj_only() {
  EnergyOptions e;
  e.cutoff = 8.5;
  e.coulomb = false;
  e.tail = false;
  return e;
}

// Born matrix of a pair potential by the lattice sum (1/V) Σ_pairs (φ'' − φ'/r) d_a d_b d_c d_d / r² (GPa).
Mat6 lattice_born(const System& s, double eps, double sig, double rc) {
  Mat6 B{};
  const double L = s.cell.a[0], V = s.cell.volume();
  for (size_t i = 0; i < s.atoms.size(); ++i)
    for (size_t j = i + 1; j < s.atoms.size(); ++j) {
      Vec3 d = s.atoms[j].pos - s.atoms[i].pos;
      for (int k = 0; k < 3; ++k) d[k] -= L * std::round(d[k] / L);
      const double r = norm(d);
      if (r >= rc) continue;
      const double sr6 = std::pow(sig / r, 6);
      const double d1 = 4 * eps * (-12 * sr6 * sr6 + 6 * sr6) / r;         // φ'
      const double d2 = 4 * eps * (156 * sr6 * sr6 - 42 * sr6) / (r * r);  // φ''
      const double c = (d2 - d1 / r) / (r * r);
      for (int I = 0; I < 6; ++I)
        for (int J = 0; J < 6; ++J) B[I][J] += c * d[VA[I]] * d[VB[I]] * d[VA[J]] * d[VB[J]] / V * kGPa;
    }
  return B;
}

}  // namespace

TEST(Mechanics, BornMatrixMatchesLatticeSum) {
  const System s = argon(4, 5.40);   // slightly expanded: under tension, so stress terms matter
  auto ff = argon_ff(s.atoms.size());
  Evaluator ev(*ff, lj_only());
  std::vector<double> x;
  for (const auto& a : s.atoms) x.insert(x.end(), a.pos.begin(), a.pos.end());
  const Mat6 B = born_matrix(ev, x, s.cell, 1e-5);
  const Mat6 R = lattice_born(s, 0.238, 3.405, 8.5);
  for (int I = 0; I < 6; ++I)
    for (int J = 0; J < 6; ++J) EXPECT_NEAR(B[I][J], R[I][J], 1e-4 * std::fabs(R[0][0])) << I << J;
  // cubic symmetry: C11 = C22 = C33, C12 = C13 = C23, C44 = C55 = C66; for central forces C12 = C44 (Cauchy relation)
  EXPECT_NEAR(B[0][1], B[3][3], 1e-4 * B[0][0]);
}

// A Bravais lattice has no internal relaxation: the static constants equal the Born term plus the stress terms of
// the Cauchy stress, C_IJ = C^B_IJ + (E_J σ0 + σ0 E_J)_I − σ0_I tr E_J.
TEST(Mechanics, StaticConstantsOfBravaisLattice) {
  const System s = argon(4, 5.40);
  auto ff = argon_ff(s.atoms.size());
  StaticElasticOptions o;
  o.field = ff;
  o.energy = lj_only();
  o.strain = 1e-5;
  const ElasticResult r = static_elastic({s}, o);
  const Mat6 B = lattice_born(s, 0.238, 3.405, 8.5);
  const double p0 = r.prestress[0];   // hydrostatic for a cubic lattice
  EXPECT_NEAR(r.prestress[1], p0, 1e-6);
  EXPECT_GT(p0, 0.01);   // tension, GPa
  for (int I = 0; I < 6; ++I)
    for (int J = 0; J < 6; ++J) {
      // σ0 = p0 I: (E σ0 + σ0 E) = 2 p0 E, E_J in Voigt: 1 on the diagonal (J ≤ 3), ½ for shear; tr E_J = [J ≤ 3]
      const double e = I == J ? (J < 3 ? 1.0 : 0.5) : 0.0;
      const double expect = B[I][J] + 2 * p0 * e - (I < 3 ? p0 : 0.0) * (J < 3 ? 1.0 : 0.0);
      EXPECT_NEAR(r.C[I][J], expect, 2e-3 * B[0][0]) << I << J;
    }
  EXPECT_GT(r.K_hill, 0);
}

// At low temperature the fluctuation constants of a crystal approach its static constants.
TEST(Mechanics, FluctuationConstantsOfColdCrystal) {
  // zero-pressure lattice for this cut-off: a ≈ 5.30 Å
  const System s0 = argon(3, 5.30);
  auto ff = argon_ff(s0.atoms.size());
  EnergyOptions e = lj_only();
  e.cutoff = 7.5;
  Trajectory t;
  t.topology = s0;
  System s = s0;
  DynamicsOptions d;
  d.field = ff;
  d.energy = e;
  d.dt = 5;
  d.steps = 8000;
  d.temperature = 5;
  d.thermostat = Thermostat::Langevin;
  d.tau_t = 1000;
  d.seed = 3;
  d.frame_every = 20;
  d.frame = [&](const std::vector<double>& x, const Cell& c, int64_t step) {
    if (step < 2000) return;
    std::vector<Vec3> p(x.size() / 3);
    for (size_t i = 0; i < p.size(); ++i) p[i] = {x[3 * i], x[3 * i + 1], x[3 * i + 2]};
    t.positions.push_back(p);
    t.cells.push_back(c);
    t.timesteps.push_back(step);
  };
  run_dynamics(s, d);
  std::vector<size_t> fr(t.frames());
  for (size_t k = 0; k < fr.size(); ++k) fr[k] = k;
  FluctuationOptions fo;
  fo.ff = ff.get();
  fo.energy = e;
  fo.temperature = 5;
  const ElasticResult fl = fluctuation_elastic(t, fr, fo);
  StaticElasticOptions so;
  so.field = ff;
  so.energy = e;
  const ElasticResult st = static_elastic({s0}, so);
  EXPECT_NEAR(fl.C[0][0], st.C[0][0], 0.06 * st.C[0][0]);
  EXPECT_NEAR(fl.C[0][1], st.C[0][1], 0.10 * st.C[0][1]);
  EXPECT_NEAR(fl.C[3][3], st.C[3][3], 0.10 * st.C[3][3]);
  EXPECT_LT(fl.fluct[0][0], 0.2 * fl.born[0][0]);   // a cold crystal: fluctuations are a small correction
}

TEST(Mechanics, IsotropicAveragesOfIsotropicSolid) {
  ElasticResult r;
  const double lam = 2.0, mu = 1.0;
  for (int I = 0; I < 3; ++I)
    for (int J = 0; J < 3; ++J) r.C[I][J] = lam + (I == J ? 2 * mu : 0);
  for (int I = 3; I < 6; ++I) r.C[I][I] = mu;
  isotropic_averages(r);
  EXPECT_NEAR(r.K_voigt, lam + 2 * mu / 3, 1e-12);
  EXPECT_NEAR(r.K_reuss, r.K_voigt, 1e-12);
  EXPECT_NEAR(r.G_hill, mu, 1e-12);
  EXPECT_NEAR(r.E_hill, mu * (3 * lam + 2 * mu) / (lam + mu), 1e-12);
  EXPECT_NEAR(r.nu_hill, lam / (2 * (lam + mu)), 1e-12);
}

TEST(Mechanics, BilinearFitFindsHinge) {
  std::mt19937 rng(5);
  std::normal_distribution<double> noise(0, 2e-4);
  std::vector<double> T, v;
  for (double t = 200; t <= 500; t += 20) {
    T.push_back(t);
    v.push_back(0.95 + (t < 373 ? 2.0e-4 : 5.5e-4) * (t - 373) + noise(rng));
  }
  const BilinearFit f = fit_bilinear(T, v);
  ASSERT_TRUE(f.ok);
  EXPECT_NEAR(f.tg, 373, 8);
  EXPECT_GT(f.tg_err, 0);
  EXPECT_LT(f.tg_err, 15);
  EXPECT_NEAR(f.slope_low, 2.0e-4, 0.3e-4);
  EXPECT_NEAR(f.slope_high, 5.5e-4, 0.3e-4);
  EXPECT_TRUE(f.note.empty()) << f.note;
}

// Uniaxial tension of a cold crystal along [100] with free lateral faces: E100 = (C11 − C12)(C11 + 2C12)/(C11 + C12),
// ν = C12/(C11 + C12).
TEST(Mechanics, TensileModulusOfColdCrystal) {
  auto ff = argon_ff(108);
  EnergyOptions e = lj_only();
  e.cutoff = 7.9;   // between the fcc shells at 7.44 and 8.32 Å, so no shell crosses the cut-off while straining
  // zero-pressure lattice constant for this cut-off (bisection on the 0 K virial pressure)
  double lo = 5.1, hi = 5.5;
  for (int it = 0; it < 50; ++it) {
    const double a = 0.5 * (lo + hi);
    const System t = argon(3, a);
    Evaluator ev(*ff, e);
    std::vector<double> x, f;
    for (const auto& at : t.atoms) x.insert(x.end(), at.pos.begin(), at.pos.end());
    (ev.compute(x, t.cell, f).virial > 0 ? lo : hi) = a;   // repulsive (compressed): expand
  }
  const System s0 = argon(3, 0.5 * (lo + hi));
  StaticElasticOptions so;
  so.field = ff;
  so.energy = e;
  const ElasticResult st = static_elastic({s0}, so);
  const double c11 = st.C[0][0], c12 = st.C[0][1];
  EXPECT_NEAR(st.prestress[0], 0.0, 1e-4);
  System s = s0;
  TensileOptions o;
  o.field = ff;
  o.energy = e;
  o.axis = 0;
  o.rate = 5e-4;
  o.max_strain = 0.02;
  o.temperature = 2;
  o.dt = 5;
  o.tau_p = 500;
  o.compressibility = 1.0 / (st.K_hill / 1.01325e-4);   // atm⁻¹, so the lateral response is fast and stable
  o.thermostat = Thermostat::Langevin;
  o.tau_t = 500;
  o.fit_strain = 0.02;
  o.sample_every = 20;
  const TensileResult r = run_tensile(s, o);
  const double E = (c11 - c12) * (c11 + 2 * c12) / (c11 + c12), nu = c12 / (c11 + c12);
  EXPECT_NEAR(r.modulus, E, 0.12 * E);
  EXPECT_NEAR(r.poisson, nu, 0.1);
  EXPECT_GT(r.curve.size(), 50u);
  // the three directions averaged: a cubic crystal gives the same E along each, so the mean is E and the spread small
  System s3 = s0;
  o.axis = 3;
  const TensileResult r3 = run_tensile(s3, o);
  EXPECT_NEAR(r3.modulus, E, 0.12 * E);
  EXPECT_LT(r3.modulus_err, 0.1 * E);
  EXPECT_NEAR(r3.poisson, nu, 0.1);
  EXPECT_EQ(r3.curve.size(), r.curve.size());
  EXPECT_NEAR(r3.curve.back().strain, r.curve.back().strain, 1e-12);
  EXPECT_NE(r3.method.find("x, y and z"), std::string::npos);
}

// Green–Kubo viscosity: an Ornstein–Uhlenbeck stress (variance σ², relaxation τ) in the three off-diagonal components
// integrates to η = V/(10 kT) · 6 σ² τ; the full protocol runs on a small cell and gives a positive viscosity
TEST(Mechanics, GreenKuboViscosity) {
  std::mt19937_64 rng(7);
  std::normal_distribution<double> g(0, 1);
  const double sigma = 300.0, tau = 50.0, dt = 4.0, V = 30000.0, T = 300.0;   // atm, fs, fs, Å³, K
  const double a = std::exp(-dt / tau), b = sigma * std::sqrt(1 - a * a);
  std::vector<std::array<double, 6>> p(400000);
  double x[3] = {0, 0, 0};
  for (auto& row : p) {
    for (double& xi : x) xi = a * xi + b * g(rng);
    row = {0, 0, 0, x[0], x[1], x[2]};
  }
  const ViscosityResult r = green_kubo_viscosity(p, dt, V, T, 1.0, 5);
  const double expect = 101325.0 * 101325.0 * V * 1e-30 * 6 * sigma * sigma * tau * 1e-15 / (1.380649e-23 * T) / 10 * 1000;
  EXPECT_NEAR(r.eta, expect, 0.06 * expect);
  EXPECT_TRUE(r.plateau);
  EXPECT_GT(r.error, 0);
  EXPECT_LT(r.error, 0.1 * expect);
  EXPECT_NEAR(r.acf.front(), 1.0, 1e-12);
  // the protocol, briefly
  GrowOptions go;
  go.chains = 3;
  go.dp = 4;
  go.density = 0.5;
  go.seed = 2;
  System s = grow(go);
  RelaxOptions ro;
  ro.ftol = 2;
  relax(s, ro);
  ViscosityOptions vo;
  vo.ps = 4;
  vo.equilibrate_ps = 1;
  vo.corr_ps = 0.5;
  vo.sample_every = 2;
  const ViscosityResult m = viscosity_green_kubo(s, vo);
  EXPECT_TRUE(std::isfinite(m.eta));
  EXPECT_EQ(m.t_ps.size(), m.running.size());
}
