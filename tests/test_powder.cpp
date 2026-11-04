#include <gtest/gtest.h>

#include <cmath>
#include <cstdio>
#include <random>

#include "caps/powder.hpp"
#include "caps/properties.hpp"

using namespace caps;

namespace {
constexpr double kPi = 3.14159265358979323846;
// a cubic cell of edge a with atoms (element, fractional position), repeated n times along each axis
Trajectory cubic(double a, const std::vector<std::pair<int, Vec3>>& basis, int n = 1) {
  Trajectory t;
  t.topology.cell.a = {a * n, 0, 0}, t.topology.cell.b = {0, a * n, 0}, t.topology.cell.c = {0, 0, a * n};
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      for (int k = 0; k < n; ++k)
        for (const auto& [z, f] : basis) {
          Atom at;
          at.element = z;
          at.pos = {a * (i + f[0]), a * (j + f[1]), a * (k + f[2])};
          t.topology.atoms.push_back(at);
        }
  return t;
}
const std::vector<Vec3> kFcc = {{0, 0, 0}, {0.5, 0.5, 0}, {0.5, 0, 0.5}, {0, 0.5, 0.5}};
double lp(double th) { const double c2 = std::cos(2 * th); return (1 + c2 * c2) / 2 / (std::sin(th) * std::sin(th) * std::cos(th)); }
const PowderPeak* find(const PowderReport& r, int h, int k, int l) {
  for (const auto& p : r.peaks) {
    int a[3] = {std::abs(p.h), std::abs(p.k), std::abs(p.l)}, b[3] = {std::abs(h), std::abs(k), std::abs(l)};
    std::sort(a, a + 3), std::sort(b, b + 3);
    if (a[0] == b[0] && a[1] == b[1] && a[2] == b[2]) return &p;
  }
  return nullptr;
}
// the reflection of exactly ±(h k l)
const PowderPeak* exact(const PowderReport& r, int h, int k, int l) {
  for (const auto& p : r.peaks)
    if ((p.h == h && p.k == k && p.l == l) || (p.h == -h && p.k == -k && p.l == -l)) return &p;
  return nullptr;
}
}  // namespace

// fcc copper: reflections where Bragg's law puts them, the multiplicities of the cubic forms, the mixed-parity ones absent,
// and every intensity m · |4 f(q)|² · LP(θ) exactly
TEST(Powder, FccCopperIsExact) {
  const double a = 3.615, lam = 1.5406;
  std::vector<std::pair<int, Vec3>> basis;
  for (const auto& f : kFcc) basis.push_back({29, f});
  const auto t = cubic(a, basis);
  PowderOptions o;
  o.wavelength = lam;
  o.two_theta_max = 120;   // 400 at 117°
  o.profile = "gaussian";
  const auto r = powder_pattern(t, {}, o);
  const int want[6][4] = {{1, 1, 1, 8}, {2, 0, 0, 6}, {2, 2, 0, 12}, {3, 1, 1, 24}, {2, 2, 2, 8}, {4, 0, 0, 6}};
  ASSERT_EQ(r.reflections, 6);
  auto expect_i = [&](int h, int k, int l, int m) {
    const double q = 2 * kPi * std::sqrt(double(h * h + k * k + l * l)) / a, th = std::asin(lam * q / (4 * kPi));
    const double f = 4 * xray_f(29, q);
    return m * f * f * lp(th);
  };
  const double i111 = expect_i(1, 1, 1, 8);
  for (const auto& w : want) {
    const auto* p = find(r, w[0], w[1], w[2]);
    ASSERT_NE(p, nullptr) << w[0] << w[1] << w[2];
    const double d = a / std::sqrt(double(w[0] * w[0] + w[1] * w[1] + w[2] * w[2]));
    EXPECT_NEAR(p->d, d, 1e-9);
    EXPECT_NEAR(p->two_theta, 2 * std::asin(lam / (2 * d)) * 180 / kPi, 1e-9);
    EXPECT_EQ(p->multiplicity, w[3]);
    EXPECT_NEAR(p->intensity / 100, expect_i(w[0], w[1], w[2], w[3]) / i111 * (find(r, 1, 1, 1)->intensity / 100), 1e-9) << w[0] << w[1] << w[2];
  }
  EXPECT_EQ(find(r, 1, 0, 0), nullptr);
  EXPECT_EQ(find(r, 1, 1, 0), nullptr);
  EXPECT_EQ(find(r, 2, 1, 0), nullptr);
}

// rock salt: the odd reflections carry f_Na − f_Cl, the even ones f_Na + f_Cl
TEST(Powder, RockSaltAlternates) {
  const double a = 5.640, lam = 1.5406;
  std::vector<std::pair<int, Vec3>> basis;
  for (const auto& f : kFcc) basis.push_back({11, f}), basis.push_back({17, Vec3{f[0] + 0.5, f[1], f[2]}});
  const auto t = cubic(a, basis);
  PowderOptions o;
  o.wavelength = lam;
  o.two_theta_max = 60;
  const auto r = powder_pattern(t, {}, o);
  const auto *p111 = find(r, 1, 1, 1), *p200 = find(r, 2, 0, 0);
  ASSERT_TRUE(p111 && p200);
  auto q = [&](int s) { return 2 * kPi * std::sqrt(double(s)) / a; };
  const double fo = 4 * (xray_f(11, q(3)) - xray_f(17, q(3))), fe = 4 * (xray_f(11, q(4)) + xray_f(17, q(4)));
  const double ratio = 8 * fo * fo * lp(std::asin(lam * q(3) / (4 * kPi))) / (6 * fe * fe * lp(std::asin(lam * q(4) / (4 * kPi))));
  std::printf("NaCl I(111)/I(200) %.5f (analytic %.5f)\n", p111->intensity / p200->intensity, ratio);
  EXPECT_NEAR(p111->intensity / p200->intensity, ratio, 1e-9);
  EXPECT_NEAR(p111->f, std::fabs(fo), 1e-9);
}

// the Kα doublet: Kα2 at its own Bragg angle with half the intensity; a Lorentzian of Scherrer width Kλ / (L cos θ)
TEST(Powder, DoubletAndScherrerWidth) {
  std::vector<std::pair<int, Vec3>> basis;
  for (const auto& f : kFcc) basis.push_back({29, f});
  const auto t = cubic(3.615, basis);
  PowderOptions o;
  o.anode = "Cu";
  o.profile = "lorentzian";
  o.W = 0;
  o.size = 150;
  o.step = 0.002;
  o.two_theta_min = 40, o.two_theta_max = 46;
  const auto r = powder_pattern(t, {}, o);
  ASSERT_EQ(r.wavelengths.size(), 2u);
  const auto* p = find(r, 1, 1, 1);
  ASSERT_NE(p, nullptr);
  const double th = p->two_theta / 2 * kPi / 180;
  const double scherrer = 0.9 * r.wavelengths[0] / (150 * std::cos(th)) * 180 / kPi;
  EXPECT_NEAR(p->fwhm, scherrer, 1e-9);
  // the pattern: the peak's half-maximum width, and the Kα2 peak's height
  const double th2 = std::asin(r.wavelengths[1] / r.wavelengths[0] * std::sin(th)) * 2 * 180 / kPi;
  size_t i1 = 0, i2 = 0;
  for (size_t i = 0; i < r.two_theta.size(); ++i) {
    if (std::fabs(r.two_theta[i] - p->two_theta) < std::fabs(r.two_theta[i1] - p->two_theta)) i1 = i;
    if (std::fabs(r.two_theta[i] - th2) < std::fabs(r.two_theta[i2] - th2)) i2 = i;
  }
  // heights: Kα1 peak = I/(πH/2) + Kα2 tail there; solve the two-Lorentzian heights from the two points
  const double H = p->fwhm, d = th2 - p->two_theta;
  const double tail = 1 / (1 + 4 * d * d / (H * H));   // one Lorentzian's height fraction at the other's centre
  const double y1 = r.intensity[i1], y2 = r.intensity[i2];
  // y1 = A1 + A2 tail', y2 = A2 + A1 tail (equal widths, near enough at Δθ ≈ 0.1°)
  const double A1 = (y1 - y2 * tail) / (1 - tail * tail), A2 = (y2 - y1 * tail) / (1 - tail * tail);
  std::printf("Cu 111: FWHM %.4f° (Scherrer %.4f°) · Kα2/Kα1 heights %.4f\n", p->fwhm, scherrer, A2 / A1);
  EXPECT_NEAR(A2 / A1, 0.5, 0.02);
}

// March–Dollase: along [001] with r, the (001) family (α = 0) and the (100) family (α = 90°) change by r⁻³ and r^(3/2)
TEST(Powder, MarchDollaseFactor) {
  Trajectory t;
  t.topology.cell.a = {4, 0, 0}, t.topology.cell.b = {0, 4, 0}, t.topology.cell.c = {0, 0, 6};
  Atom at;
  at.element = 6;
  at.pos = {0, 0, 0};
  t.topology.atoms.push_back(at);
  PowderOptions o;
  o.wavelength = 1.5406;
  o.two_theta_max = 30;
  const auto r0 = powder_pattern(t, {}, o);
  o.preferred = {0, 0, 1};
  o.march_r = 0.5;
  const auto r1 = powder_pattern(t, {}, o);
  auto ratio = [&](const PowderReport& r) {
    const auto *c = exact(r, 0, 0, 1), *a = exact(r, 1, 0, 0) ? exact(r, 1, 0, 0) : exact(r, 0, 1, 0);
    EXPECT_TRUE(c && a);
    return c && a ? c->intensity / a->intensity : 0.0;
  };
  EXPECT_NEAR(ratio(r1) / ratio(r0), std::pow(0.5, -4.5), 1e-9);
}

// the frames' thermal motion: random displacements of σ per axis give the Bragg intensities of the still crystal
// with B = 8π²σ² (Debye–Waller), up to the diffuse part (1/N of it at the Bragg positions)
TEST(Powder, ThermalFramesMatchDebyeWaller) {
  std::vector<std::pair<int, Vec3>> basis;
  for (const auto& f : kFcc) basis.push_back({29, f});
  auto t = cubic(3.615, basis, 4);
  const double sigma = 0.09;
  std::mt19937_64 rng(3);
  std::normal_distribution<double> N(0.0, sigma);
  for (int f = 0; f < 40; ++f) {
    std::vector<Vec3> p;
    for (const auto& a : t.topology.atoms) p.push_back(a.pos + Vec3{N(rng), N(rng), N(rng)});
    t.positions.push_back(p), t.cells.push_back(t.topology.cell), t.timesteps.push_back(f);
  }
  std::vector<size_t> fr(40);
  for (size_t k = 0; k < 40; ++k) fr[k] = k;
  PowderOptions o;
  o.wavelength = 1.5406;
  o.two_theta_max = 120;
  o.max_peaks = -1;
  const auto hot = powder_pattern(t, fr, o);
  Trajectory still = t;
  still.positions.clear(), still.cells.clear(), still.timesteps.clear();
  o.b_iso = 8 * kPi * kPi * sigma * sigma;
  const auto cold = powder_pattern(still, {}, o);
  for (auto hkl : std::vector<std::array<int, 3>>{{8, 0, 0}, {8, 8, 0}, {12, 4, 4}}) {   // 200, 220, 311 of the unit cell in the 4×4×4 cell
    const auto *a = find(hot, hkl[0], hkl[1], hkl[2]), *b = find(cold, hkl[0], hkl[1], hkl[2]);
    ASSERT_TRUE(a && b);
    const double ra = a->intensity / find(hot, 4, 4, 4)->intensity, rb = b->intensity / find(cold, 4, 4, 4)->intensity;
    std::printf("Cu (%d%d%d)/(444) frames %.4f · B %.4f\n", hkl[0] / 4, hkl[1] / 4, hkl[2] / 4, ra, rb);
    EXPECT_NEAR(ra, rb, 0.04 * rb);
  }
}
