#include <gtest/gtest.h>

#include <cmath>

#include "caps/cell_relax.hpp"
#include "caps/qha.hpp"
#include "caps/uff.hpp"

using namespace caps;

namespace {
constexpr double kB = 1.987204e-3;   // kcal/mol/K
constexpr double kGPa = 6.947695;    // GPa per kcal/mol/Å³
// A model with one Einstein frequency for all n modes and a Grüneisen γ: θ(V) = θ₀ (V/V₀)^(−γ), E(V) = ½ B₀ (V − V₀)²/V₀
struct Model {
  double V0 = 100, B0 = 5.0, theta0 = 200, gamma = 2, n = 12;   // Å³, kcal/mol/Å³, K, —, modes
  double theta(double V) const { return theta0 * std::pow(V / V0, -gamma); }
  double Fvib(double V, double T) const {
    const double e = kB * theta(V);
    return n * (0.5 * e + (T > 0 ? kB * T * std::log(-std::expm1(-theta(V) / T)) : 0.0));
  }
  double Cv(double V, double T) const {
    if (T <= 0) return 0;
    const double x = theta(V) / T, ex = std::exp(-x);
    return n * 1000 * kB * x * x * ex / ((1 - ex) * (1 - ex));
  }
  double S(double V, double T) const {
    if (T <= 0) return 0;
    const double x = theta(V) / T;
    return n * 1000 * kB * (x / std::expm1(x) - std::log(-std::expm1(-x)));
  }
  double E(double V) const { return 0.5 * B0 * (V - V0) * (V - V0) / V0; }
  double G(double V, double T) const { return E(V) + Fvib(V, T); }
  // the exact minimum by bisection on central differences of G
  double Vmin(double T) const {
    double a = 0.9 * V0, b = 1.1 * V0;
    for (int k = 0; k < 200; ++k) {
      const double m = (a + b) / 2, h = 1e-3;
      if (G(m + h, T) > G(m - h, T)) b = m; else a = m;
    }
    return (a + b) / 2;
  }
};
std::vector<QhaSample> model_samples(const Model& M, const std::vector<double>& Ts, double range, int pts) {
  std::vector<QhaSample> out;
  for (int k = 0; k < pts; ++k) {
    QhaSample s;
    const double e = -range + 2 * range * k / (pts - 1);
    s.strain = {e, e, e};
    s.volume = M.V0 * std::pow(1 + e, 3);
    s.energy = M.E(s.volume);
    for (double T : Ts) {
      const double d = qha_step(T);
      const std::array<double, 3> t3 = T > 0 ? std::array<double, 3>{T - d, T, T + d} : std::array<double, 3>{0, d, 2 * d};
      s.F.push_back({M.Fvib(s.volume, t3[0]), M.Fvib(s.volume, t3[1]), M.Fvib(s.volume, t3[2])});
      s.Cv.push_back({M.Cv(s.volume, t3[0]), M.Cv(s.volume, t3[1]), M.Cv(s.volume, t3[2])});
      s.S.push_back({M.S(s.volume, t3[0]), M.S(s.volume, t3[1]), M.S(s.volume, t3[2])});
    }
    out.push_back(s);
  }
  return out;
}

constexpr double kEps = 0.2385, kSig = 3.405, kRc = 13;   // a shell crossing a short cut-off puts a kink in E(V)
ForceField argon_ff(const System& s) {
  ForceField ff = default_forcefield(s);
  ff.type_names = {"Ar"};
  ff.lj.assign(1, {kEps, kSig});
  for (size_t i = 0; i < s.atoms.size(); ++i) ff.type_index[i] = 0, ff.atom_type[i] = "Ar";
  ff.bonds.clear(), ff.angles.clear(), ff.dihedrals.clear(), ff.impropers.clear(), ff.inversions.clear(), ff.bonds_x.clear(), ff.angles_x.clear();
  ff.pairs14.clear();
  for (auto& e : ff.excluded) e.clear();
  ff.charge.assign(s.atoms.size(), 0.0);
  return ff;
}
System argon(double a) {
  System s;
  s.cell.a = {a, 0, 0}, s.cell.b = {0, a, 0}, s.cell.c = {0, 0, a};
  for (Vec3 f : {Vec3{0, 0, 0}, Vec3{0.5, 0.5, 0}, Vec3{0.5, 0, 0.5}, Vec3{0, 0.5, 0.5}}) {
    Atom x;
    x.element = 18, x.pos = f * a;
    s.atoms.push_back(x);
  }
  return s;
}
QhaOptions argon_options(const std::string& mode) {
  QhaOptions o;
  o.phonon.field_for = argon_ff;
  o.phonon.energy.cutoff = kRc, o.phonon.energy.coulomb = false, o.phonon.energy.tail = false;
  o.phonon.mesh = {6, 6, 6};
  o.mode = mode, o.range = 0.04;
  o.temperatures = {0, 20, 40, 60};
  return o;
}
}  // namespace

// The fit against the exact model: V(T) from an independent bisection on G; the Grüneisen relation α B_T V = γ C_v that the
// model obeys exactly (γ constant, B_T from G''); C_v and S at the exact minimum; α = 0 at T = 0; C_p − C_v = T V B α²
TEST(Qha, ExactGruneisenModel) {
  const Model M;
  const std::vector<double> Ts{0, 50, 150, 300, 500};
  const auto r = qha_analyze(model_samples(M, Ts, 0.03, 9), Ts, "volume", {true, true, true}, {4.6416, 4.6416, 4.6416}, M.V0, 0);
  ASSERT_EQ(r.rows.size(), Ts.size());
  for (const auto& row : r.rows) {
    const double Vx = M.Vmin(row.T);
    EXPECT_NEAR(row.volume / Vx, 1, 2e-5) << row.T;
    EXPECT_NEAR(row.Cv, M.Cv(Vx, row.T), 2e-3 * std::max(1.0, M.Cv(Vx, row.T))) << row.T;
    EXPECT_NEAR(row.S, M.S(Vx, row.T), 2e-3 * std::max(1.0, M.S(Vx, row.T))) << row.T;
    const double h = 0.01, G2 = (M.G(Vx + h, row.T) - 2 * M.G(Vx, row.T) + M.G(Vx - h, row.T)) / (h * h);
    const double Bx = Vx * G2 * kGPa;
    EXPECT_NEAR(row.B, Bx, 2e-3 * Bx) << row.T;
    if (row.T > 0) {
      const double lhs = row.alpha_V * (row.B / kGPa) * row.volume * 1000;   // cal/mol/K
      EXPECT_NEAR(lhs, M.gamma * row.Cv, 5e-3 * M.gamma * row.Cv) << row.T;
      EXPECT_NEAR(row.Cp - row.Cv, row.T * row.volume * (row.B / kGPa) * row.alpha_V * row.alpha_V * 1000, 1e-6);
      EXPECT_GT(row.Cp, row.Cv);
      EXPECT_NEAR(row.alpha[0], row.alpha_V / 3, 2e-3 * row.alpha_V);
    } else {
      EXPECT_EQ(row.alpha_V, 0);
      EXPECT_GT(row.volume, M.V0);   // zero-point expansion
    }
    EXPECT_FALSE(row.extrapolated);
  }
}

// A pressure shifts the minimum to the root of dF/dV = −P, within the fit's accuracy
TEST(Qha, PressureShiftsTheMinimum) {
  const Model M;
  const std::vector<double> Ts{300};
  const double Patm = 5000, P = Patm * 1.458397e-5;
  const auto r = qha_analyze(model_samples(M, Ts, 0.03, 9), Ts, "volume", {true, true, true}, {1, 1, 1}, M.V0, Patm);
  double a = 0.9 * M.V0, b = 1.1 * M.V0;
  for (int k = 0; k < 200; ++k) {
    const double m = (a + b) / 2, h = 1e-3;
    if (M.G(m + h, 300) + P * (m + h) > M.G(m - h, 300) + P * (m - h)) b = m; else a = m;
  }
  EXPECT_NEAR(r.rows[0].volume / ((a + b) / 2), 1, 2e-5);
  EXPECT_LT(r.rows[0].volume, M.Vmin(300));
}

// Fcc argon (Lennard-Jones, phonons on a 6 × 6 × 6 mesh): zero-point expansion above the static cell, expansion on heating,
// C_v rising towards 3R per atom (classical limit), C_p above C_v, and the two modes agreeing: axes mode finds the same
// volume and α_a = α_b = α_c = α_V/3 for this cubic crystal
TEST(Qha, ArgonVolumeAndAxes) {
  const System a0 = argon(5.3);
  const auto v = quasi_harmonic(a0, argon_options("volume"));
  ASSERT_EQ(v.rows.size(), 4u);
  EXPECT_GT(v.rows[0].volume, v.volume0);
  for (size_t k = 1; k < v.rows.size(); ++k) {
    EXPECT_GT(v.rows[k].volume, v.rows[k - 1].volume);
    EXPECT_GT(v.rows[k].alpha_V, 0);
    EXPECT_GT(v.rows[k].Cv, v.rows[k - 1].Cv);
    EXPECT_GT(v.rows[k].Cp, v.rows[k].Cv);
    EXPECT_LT(v.rows[k].B, v.rows[k - 1].B);
  }
  const double dulong = 4 * 3 * 1.987204;   // cal/mol/K per cell of four atoms
  EXPECT_GT(v.rows[3].Cv, 0.85 * dulong);
  EXPECT_LT(v.rows[3].Cv, dulong);
  const auto x = quasi_harmonic(a0, argon_options("axes"));
  for (size_t k = 0; k < x.rows.size(); ++k) {
    EXPECT_NEAR(x.rows[k].volume / v.rows[k].volume, 1, 2e-3) << k;
    if (k > 0) {
      for (int d = 0; d < 3; ++d) EXPECT_NEAR(x.rows[k].alpha[size_t(d)], x.rows[k].alpha_V / 3, 0.02 * x.rows[k].alpha_V) << k;
      EXPECT_NEAR(x.rows[k].alpha_V / v.rows[k].alpha_V, 1, 0.1) << k;
      EXPECT_NEAR(x.rows[k].B / v.rows[k].B, 1, 0.1) << k;
    }
  }
}

// A locked axis stays at its length in axes mode
TEST(Qha, LockedAxis) {
  QhaOptions o = argon_options("axes");
  o.axes = {true, true, false};
  o.temperatures = {40};
  const auto r = quasi_harmonic(argon(5.3), o);
  EXPECT_DOUBLE_EQ(r.rows[0].length[2], r.lengths0[2]);
  EXPECT_EQ(r.rows[0].alpha[2], 0);
  EXPECT_GT(r.rows[0].alpha[0], 0);
  EXPECT_EQ(r.samples.size(), 1u + 2 * 4 + 4);
}
