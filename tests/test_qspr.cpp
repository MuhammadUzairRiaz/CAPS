#include <gtest/gtest.h>

#include <cmath>
#include <map>

#include "caps/molecule.hpp"
#include "caps/qspr.hpp"

using namespace caps;

namespace {
std::map<std::string, double> desc(const std::string& smiles) {
  std::map<std::string, double> m;
  for (const auto& d : molecular_descriptors(build_molecule(smiles).system)) m[d.name] = d.value;
  return m;
}
}  // namespace

// Graph indices against their published values: Wiener's hexanes (35, 32, 31, 28, 29), Zagreb, Randić and the Kier
// shape indices of n-hexane, Balaban's J of n-hexane (2.3391), Kier–Hall valence χ of ethanol and chloromethane
TEST(Qspr, GraphIndicesArePublishedValues) {
  const auto hx = desc("CCCCCC");
  EXPECT_EQ(hx.at("W"), 35);
  EXPECT_EQ(desc("CC(C)CCC").at("W"), 32);
  EXPECT_EQ(desc("CCC(C)CC").at("W"), 31);
  EXPECT_EQ(desc("CC(C)(C)CC").at("W"), 28);
  EXPECT_EQ(desc("CC(C)C(C)C").at("W"), 29);
  EXPECT_EQ(hx.at("M1"), 18);
  EXPECT_EQ(hx.at("M2"), 16);
  EXPECT_NEAR(hx.at("1chi"), 2 / std::sqrt(2.0) + 1.5, 1e-12);
  EXPECT_NEAR(hx.at("kappa1"), 6, 1e-12);
  EXPECT_NEAR(hx.at("kappa2"), 5, 1e-12);
  EXPECT_NEAR(hx.at("kappa3"), 3.0 * 16 / 9, 1e-12);
  const double J = 5 * (2 / std::sqrt(165.0) + 2 / std::sqrt(99.0) + 1 / 9.0);
  EXPECT_NEAR(hx.at("J"), J, 1e-12);
  EXPECT_NEAR(J, 2.3391, 1e-4);
  const auto et = desc("CCO");
  EXPECT_NEAR(et.at("1chiv"), 1 / std::sqrt(2.0) + 1 / std::sqrt(10.0), 1e-12);
  EXPECT_NEAR(et.at("0chiv"), 1 + 1 / std::sqrt(2.0) + 1 / std::sqrt(5.0), 1e-12);
  EXPECT_EQ(et.at("hbond_donors"), 1);
  EXPECT_EQ(et.at("hbond_acceptors"), 1);
  EXPECT_NEAR(desc("CCl").at("0chiv"), 1 + 1 / std::sqrt(7.0 / 9.0), 1e-12);   // Cl: δv = (7 − 0)/(17 − 7 − 1)
  std::printf("n-hexane: W %.0f, 1χ %.4f, κ %.3f %.3f %.3f, J %.4f; ethanol 1χv %.4f\n", hx.at("W"), hx.at("1chi"), hx.at("kappa1"), hx.at("kappa2"),
              hx.at("kappa3"), hx.at("J"), et.at("1chiv"));
}

// A single atom: its van der Waals volume is the sphere's and its shadow fills π/4 of the enclosing square
TEST(Qspr, SphereVolumeAndShadow) {
  System s;
  Atom a;
  a.element = 18;
  s.atoms.push_back(a);
  std::map<std::string, double> m;
  for (const auto& d : molecular_descriptors(s)) m[d.name] = d.value;
  const double r = 1.88;
  EXPECT_NEAR(m.at("vdw_volume"), 4.0 / 3.0 * M_PI * r * r * r, 3e-3 * m.at("vdw_volume"));
  EXPECT_NEAR(m.at("vdw_area"), 4 * M_PI * r * r, 1e-2 * m.at("vdw_area"));
  EXPECT_NEAR(m.at("shadow_XY_ratio"), M_PI / 4, 3e-3);
}

namespace {
// deterministic pseudo-random numbers in (−1, 1)
double rnd(int i, int k) { return std::fmod(std::sin(12.9898 * i + 78.233 * k) * 43758.5453, 1.0); }
}  // namespace

// MLR reproduces an exact plane; PLS with every component is OLS; leave-one-out PRESS equals Σ (e/(1 − h))²; Cook's
// distances equal the change of the fitted values when each row is left out
TEST(Qspr, LinearModelsAndDiagnostics) {
  const int n = 12;
  std::vector<double> y(n), x1(n), x2(n), x3(n);
  for (int i = 0; i < n; ++i) x1[i] = rnd(i, 1), x2[i] = rnd(i, 2), x3[i] = rnd(i, 3), y[i] = 1 + 2 * x1[i] - 3 * x2[i];
  ModelOptions o;
  ModelResult r = fit_model(y, {x1, x2}, {"a", "b"}, o);
  EXPECT_NEAR(r.terms[0].coef, 1, 1e-10);
  EXPECT_NEAR(r.terms[1].coef, 2, 1e-10);
  EXPECT_NEAR(r.terms[2].coef, -3, 1e-10);
  EXPECT_NEAR(r.r2, 1, 1e-12);
  for (int i = 0; i < n; ++i) y[i] += 0.3 * rnd(i, 4) + 0.5 * x3[i];
  o.folds = n;
  r = fit_model(y, {x1, x2, x3}, {"a", "b", "c"}, o);
  double press = 0;
  for (int i = 0; i < n; ++i) press += std::pow(r.residuals[size_t(i)] / (1 - r.leverage[size_t(i)]), 2);
  double sst = 0, my = 0;
  for (double v : y) my += v / n;
  for (double v : y) sst += (v - my) * (v - my);
  EXPECT_NEAR(r.q2, 1 - press / sst, 1e-10);
  // Cook's distance by deletion
  double s2 = 0;
  for (double e : r.residuals) s2 += e * e;
  s2 /= double(n - 4);
  for (int i : {0, 5, 11}) {
    std::vector<double> yy;
    std::vector<std::vector<double>> xx(3);
    for (int j = 0; j < n; ++j)
      if (j != i) yy.push_back(y[j]), xx[0].push_back(x1[j]), xx[1].push_back(x2[j]), xx[2].push_back(x3[j]);
    ModelOptions oo;
    const ModelResult d = fit_model(yy, xx, {"a", "b", "c"}, oo);
    double sum = 0;
    for (int j = 0; j < n; ++j) {
      const double f = d.terms[0].coef + d.terms[1].coef * x1[j] + d.terms[2].coef * x2[j] + d.terms[3].coef * x3[j];
      sum += std::pow(r.fitted[size_t(j)] - f, 2);
    }
    EXPECT_NEAR(r.cooks[size_t(i)], sum / (4 * s2), 1e-10);
  }
  ModelOptions op;
  op.method = "pls";
  op.components = 3;
  const ModelResult p = fit_model(y, {x1, x2, x3}, {"a", "b", "c"}, op);
  for (size_t k = 0; k < 4; ++k) EXPECT_NEAR(p.terms[k].coef, r.terms[k].coef, 1e-9);
  std::printf("MLR R² %.4f q²(LOO) %.4f; PLS(3) coefficients equal OLS\n", r.r2, r.q2);
}

// Stepwise AIC keeps the one column y depends on; y-randomisation leaves the shuffled fits far below it
TEST(Qspr, StepwiseAndRandomisation) {
  const int n = 30;
  std::vector<double> y(n), x1(n), x2(n), x3(n);
  for (int i = 0; i < n; ++i) x1[i] = rnd(i, 1), x2[i] = rnd(i, 2), x3[i] = rnd(i, 3), y[i] = 4 * x1[i] + 0.05 * rnd(i, 5);
  ModelOptions o;
  o.stepwise = true, o.randomise = 20, o.folds = 5;
  const ModelResult r = fit_model(y, {x1, x2, x3}, {"a", "b", "c"}, o);
  ASSERT_GE(r.terms.size(), 2u);
  EXPECT_EQ(r.terms[1].name, "a");
  double mx = 0;
  for (double v : r.random_r2) mx = std::max(mx, v);
  EXPECT_EQ(r.random_r2.size(), 20u);
  EXPECT_LT(mx, 0.5);
  EXPECT_GT(r.q2, 0.99);
  std::printf("stepwise kept %zu term(s); R² %.4f, q² %.4f; shuffled R² at most %.3f\n", r.terms.size() - 1, r.r2, r.q2, mx);
}

// GFA finds the hinge y = 3 ⟨x − 0.5⟩₊ (knot at a data value); a small network fits a smooth curve
TEST(Qspr, GfaHingeAndNetwork) {
  const int n = 21;
  std::vector<double> y(n), x1(n), x2(n);
  for (int i = 0; i < n; ++i) x1[i] = i / 20.0, x2[i] = rnd(i, 7), y[i] = 3 * std::max(0.0, x1[i] - 0.5);
  ModelOptions o;
  o.method = "gfa";
  o.gfa_generations = 4000;
  const ModelResult g = fit_model(y, {x1, x2}, {"x", "z"}, o);
  std::printf("GFA: %s (LOF %.3g)\n", g.equation.c_str(), g.lof);
  EXPECT_LT(g.rmse, 1e-8);
  bool hinge = false;
  for (const auto& t : g.terms) hinge = hinge || t.name == "<x-0.5>+";
  EXPECT_TRUE(hinge);
  std::vector<double> ys(30), xs(30);
  for (int i = 0; i < 30; ++i) xs[i] = i / 29.0 * 3, ys[i] = std::sin(2 * xs[i]);
  ModelOptions on;
  on.method = "nn", on.nn_hidden = 6, on.nn_decay = 1e-6;
  const ModelResult nnr = fit_model(ys, {xs}, {"x"}, on);
  std::printf("network: R² %.6f\n", nnr.r2);
  EXPECT_GT(nnr.r2, 0.999);
}

// Ward's method on 0, 1, 10, 11: merges at 0.5, 0.5 and 100 (the within sum of squares added), two clusters; k-means
// agrees; Fisher's direction is S_W⁻¹(μ1 − μ2) and separates two clouds
TEST(Qspr, ClusteringAndDiscriminant) {
  const std::vector<std::vector<double>> pts = {{0}, {1}, {10}, {11}};
  const RowClusters w = cluster_rows(pts, "ward", 2, false);
  ASSERT_EQ(w.merges.size(), 3u);
  EXPECT_NEAR(w.merges[0][2], 0.5, 1e-12);
  EXPECT_NEAR(w.merges[1][2], 0.5, 1e-12);
  EXPECT_NEAR(w.merges[2][2], 100, 1e-12);
  EXPECT_EQ(w.label[0], w.label[1]);
  EXPECT_EQ(w.label[2], w.label[3]);
  EXPECT_NE(w.label[0], w.label[2]);
  EXPECT_NEAR(w.within_ss, 1.0, 1e-12);
  const RowClusters km = cluster_rows(pts, "kmeans", 2, false);
  EXPECT_EQ(km.label[0], km.label[1]);
  EXPECT_NE(km.label[0], km.label[2]);
  // two classes with a shared within-class scatter
  std::vector<std::vector<double>> X;
  std::vector<std::string> cls;
  for (int i = 0; i < 20; ++i) {
    const double a = rnd(i, 1), b = rnd(i, 2);
    X.push_back({a + 0.5 * b, b}), cls.push_back("A");
    X.push_back({a + 0.5 * b + 3, b + 1}), cls.push_back("B");
  }
  const DiscriminantResult d = discriminant(X, cls);
  EXPECT_EQ(d.accuracy, 1.0);
  EXPECT_EQ(d.accuracy_loo, 1.0);
  // S_W⁻¹ (μ_A − μ_B)
  double mA[2] = {0, 0}, mB[2] = {0, 0}, S[2][2] = {{0, 0}, {0, 0}};
  for (size_t i = 0; i < X.size(); ++i)
    for (int k = 0; k < 2; ++k) (cls[i] == "A" ? mA : mB)[k] += X[i][size_t(k)] / 20;
  for (size_t i = 0; i < X.size(); ++i) {
    const double* m = cls[i] == "A" ? mA : mB;
    for (int a = 0; a < 2; ++a)
      for (int b = 0; b < 2; ++b) S[a][b] += (X[i][size_t(a)] - m[a]) * (X[i][size_t(b)] - m[b]);
  }
  const double det = S[0][0] * S[1][1] - S[0][1] * S[1][0];
  const double dx = mA[0] - mB[0], dy = mA[1] - mB[1];
  const double v0 = (S[1][1] * dx - S[0][1] * dy) / det, v1 = (-S[1][0] * dx + S[0][0] * dy) / det;
  const auto& u = d.directions[0];
  const double cosang = std::fabs(u[0] * v0 + u[1] * v1) / (std::hypot(u[0], u[1]) * std::hypot(v0, v1));
  EXPECT_NEAR(cosang, 1.0, 1e-9);
}
