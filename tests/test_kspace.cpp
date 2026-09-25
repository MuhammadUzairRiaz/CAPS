#include <gtest/gtest.h>

#include <cmath>
#include <complex>
#include <random>
#include <vector>

#include "caps/kspace.hpp"

using namespace caps;

namespace {
constexpr double kPi = 3.14159265358979323846;

Cell cubic(double L) {
  Cell c;
  c.a = {L, 0, 0};
  c.b = {0, L, 0};
  c.c = {0, 0, L};
  return c;
}

// NaCl rock salt, n × n × n conventional cells of edge a0 (8 ions each)
void rock_salt(int n, double a0, std::vector<double>& x, std::vector<double>& q) {
  const double base[8][4] = {{0, 0, 0, 1}, {0.5, 0.5, 0, 1}, {0.5, 0, 0.5, 1}, {0, 0.5, 0.5, 1},
                             {0.5, 0, 0, -1}, {0, 0.5, 0, -1}, {0, 0, 0.5, -1}, {0.5, 0.5, 0.5, -1}};
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      for (int k = 0; k < n; ++k)
        for (const auto& b : base) {
          x.insert(x.end(), {(i + b[0]) * a0, (j + b[1]) * a0, (k + b[2]) * a0});
          q.push_back(b[3]);
        }
}

// full Ewald energy: real space by direct sum over images, self term, reciprocal from `recip`
template <class R>
double ewald_total(const std::vector<double>& x, const std::vector<double>& q, double L, double beta, double rc, R recip) {
  const size_t n = q.size();
  double er = 0;
  for (size_t i = 0; i < n; ++i)
    for (size_t j = i + 1; j < n; ++j) {
      double d[3];
      for (int k = 0; k < 3; ++k) {
        d[k] = x[3 * j + k] - x[3 * i + k];
        d[k] -= L * std::round(d[k] / L);
      }
      const double r = std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
      if (r < rc) er += kCoulombConstant * q[i] * q[j] * std::erfc(beta * r) / r;
    }
  double q2 = 0;
  for (double c : q) q2 += c * c;
  return er - kCoulombConstant * beta / std::sqrt(kPi) * q2 + recip();
}
}  // namespace

TEST(Kspace, FftMatchesTheDirectTransform) {
  for (int k1 : {6, 7, 10, 36, 45, 50}) {
    const int k2 = 5, k3 = 4;
    std::mt19937 rng(3);
    std::uniform_real_distribution<double> U(-1, 1);
    std::vector<std::complex<double>> a(size_t(k1 * k2 * k3));
    for (auto& v : a) v = {U(rng), U(rng)};
    auto b = a;
    fft3d(b, k1, k2, k3, +1);
    double worst = 0;
    for (int m1 = 0; m1 < k1; ++m1)
      for (int m2 = 0; m2 < k2; ++m2)
        for (int m3 = 0; m3 < k3; ++m3) {
          std::complex<double> s = 0;
          for (int i = 0; i < k1; ++i)
            for (int j = 0; j < k2; ++j)
              for (int k = 0; k < k3; ++k)
                s += a[size_t(i + k1 * (j + k2 * k))] * std::polar(1.0, 2 * kPi * (double(m1 * i) / k1 + double(m2 * j) / k2 + double(m3 * k) / k3));
          worst = std::max(worst, std::abs(s - b[size_t(m1 + k1 * (m2 + k2 * m3))]));
        }
    EXPECT_LT(worst, 1e-10) << k1;
  }
  EXPECT_EQ(fft_good_size(31), 32);
  EXPECT_EQ(fft_good_size(49), 49);   // 7²: radix 7 is FFT-friendly too
  EXPECT_EQ(fft_good_size(13), 14);
  EXPECT_EQ(fft_good_size(43), 45);
}

TEST(Kspace, MadelungConstantOfRockSalt) {
  // E per ion pair = −M kC / r0 with the nearest-neighbour distance r0 = a0/2 and M = 1.747565
  std::vector<double> x, q;
  const double a0 = 5.64, L = 2 * a0;
  rock_salt(2, a0, x, q);
  const Cell cell = cubic(L);
  const double rc = 5.5, beta = ewald_beta(rc, 1e-8);
  const int kmax[3] = {12, 12, 12};
  std::vector<double> f(x.size(), 0.0);
  double vir[6];
  const double e = ewald_total(x, q, L, beta, rc, [&] { return ewald_reciprocal(x, q, cell, beta, kmax, f, vir); });
  const double M = -e / (q.size() / 2) * (a0 / 2) / kCoulombConstant;
  EXPECT_NEAR(M, 1.747565, 2e-5);
  // PME on a fine grid agrees
  PmeGrid g = pme_grid(cell, beta, 0.3, 6);
  std::vector<double> f2(x.size(), 0.0);
  const double e2 = ewald_total(x, q, L, beta, rc, [&] { return pme_reciprocal(x, q, cell, g, f2, vir); });
  EXPECT_NEAR(e2, e, 1e-4 * std::fabs(e));
}

TEST(Kspace, PmeMatchesEwaldForForcesEnergyAndVirial) {
  // random neutral charges in a triclinic cell
  std::mt19937 rng(11);
  std::uniform_real_distribution<double> U(0, 1);
  Cell cell;
  cell.a = {20, 0, 0};
  cell.b = {3, 19, 0};
  cell.c = {-2, 4, 21};
  const int n = 200;
  std::vector<double> x, q;
  for (int i = 0; i < n; ++i) {
    const double s[3] = {U(rng), U(rng), U(rng)};
    for (int k = 0; k < 3; ++k) x.push_back(s[0] * cell.a[k] + s[1] * cell.b[k] + s[2] * cell.c[k]);
    q.push_back(i % 2 ? 0.4 : -0.4);
  }
  const double beta = ewald_beta(9.0, 1e-5);
  const int kmax[3] = {14, 14, 14};
  std::vector<double> fe(x.size(), 0.0), fp(x.size(), 0.0);
  double ve[6], vp[6];
  const double ee = ewald_reciprocal(x, q, cell, beta, kmax, fe, ve);
  const PmeGrid g = pme_grid(cell, beta, 0.8, 6);
  const double ep = pme_reciprocal(x, q, cell, g, fp, vp);
  EXPECT_NEAR(ep, ee, 1e-4 * std::fabs(ee) + 1e-4);
  double rms = 0, err = 0;
  for (size_t k = 0; k < x.size(); ++k) rms += fe[k] * fe[k], err += (fp[k] - fe[k]) * (fp[k] - fe[k]);
  EXPECT_LT(std::sqrt(err / rms), 1e-3);
  for (int c = 0; c < 6; ++c) EXPECT_NEAR(vp[c], ve[c], 1e-3 * std::fabs(ee) + 1e-3) << c;
  // forces are minus the gradient of the reciprocal energy (finite differences on a few coordinates)
  for (int k : {0, 7, 301}) {
    auto xp = x, xm = x;
    const double h = 1e-5;
    xp[size_t(k)] += h;
    xm[size_t(k)] -= h;
    std::vector<double> t(x.size(), 0.0);
    double vv[6];
    const double dE = (ewald_reciprocal(xp, q, cell, beta, kmax, t, vv) - ewald_reciprocal(xm, q, cell, beta, kmax, t, vv)) / (2 * h);
    EXPECT_NEAR(fe[size_t(k)], -dE, 1e-5 * (1 + std::fabs(dE))) << k;
  }
  // trace of the virial = −dE/dλ under x → λx, cell → λ cell
  auto scaled = [&](double lam) {
    auto y = x;
    for (double& v : y) v *= lam;
    Cell c2 = cell;
    c2.a = c2.a * lam, c2.b = c2.b * lam, c2.c = c2.c * lam;
    std::vector<double> t(x.size(), 0.0);
    double vv[6];
    return ewald_reciprocal(y, q, c2, beta, kmax, t, vv);
  };
  const double eps = 1e-5;
  const double dEdl = (scaled(1 + eps) - scaled(1 - eps)) / (2 * eps);
  EXPECT_NEAR(ve[0] + ve[1] + ve[2], -dEdl, 1e-5 * std::fabs(ee) + 1e-6);
}
