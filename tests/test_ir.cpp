#include <gtest/gtest.h>

#include <cmath>

#include "caps/properties.hpp"

using namespace caps;

namespace {
constexpr double kPi = 3.14159265358979323846;
constexpr double kNu0 = 1500, kA = 0.05, kQ = 1.0, kL = 20.0, kT = 300.0;   // cm⁻¹, Å, e, Å, K
constexpr double kDtFs = 1.0;
constexpr int kFrames = 20000;

// one ±q pair whose separation oscillates at ν̃0: M = q x(t), x = x0 + A cos(ω0 t), with or without exact velocities
Trajectory oscillating_dipole(bool velocities) {
  Trajectory t;
  System s;
  s.cell.a = {kL, 0, 0}, s.cell.b = {0, kL, 0}, s.cell.c = {0, 0, kL};
  for (double q : {-kQ, kQ}) {
    Atom a;
    a.element = q < 0 ? 17 : 11;
    a.charge = q;
    s.atoms.push_back(a);
  }
  s.has_charges = true;
  const double w = 2 * kPi * 2.99792458e10 * kNu0 * 1e-15;   // rad/fs
  for (int k = 0; k < kFrames; ++k) {
    const double tt = k * kDtFs, x = 2.5 + kA * std::cos(w * tt), xd = -kA * w * std::sin(w * tt);
    t.positions.push_back({Vec3{10 - x / 2, 10, 10}, Vec3{10 + x / 2, 10, 10}});
    t.cells.push_back(s.cell);
    t.timesteps.push_back(k);
    if (velocities) t.velocities.push_back({Vec3{-xd / 2, 0, 0}, Vec3{xd / 2, 0, 0}});
  }
  s.atoms[0].pos = t.positions[0][0], s.atoms[1].pos = t.positions[0][1];
  t.topology = s;
  return t;
}

double area(const Series& s) {
  double a = 0;
  for (size_t k = 1; k < s.x.size(); ++k) a += 0.5 * (s.y[k] + s.y[k - 1]) * (s.x[k] - s.x[k - 1]);
  return a;
}
}  // namespace

// The band of an oscillating dipole sits at its frequency, and its area is the classical sum rule
// ∫ α n dω = π β ⟨Ṁ²⟩ / (6 c V ε₀), ⟨Ṁ²⟩ = (q A ω0)² / 2, whether Ṁ comes from velocities or from differences of M
TEST(Infrared, OscillatingDipoleSumRule) {
  AnalyzeOptions o;
  o.temperature = kT;
  o.frame_ps = kDtFs * 1e-3;
  o.ir_max_wavenumber = 3000;
  constexpr double kE = 1.602176634e-19, kEps0 = 8.8541878128e-12, kKB = 1.380649e-23, kC = 2.99792458e8;
  const double w0 = 2 * kPi * kC * 100 * kNu0;                    // rad/s
  const double mdot2 = std::pow(kQ * kE * kA * 1e-10 * w0, 2) / 2;   // (C·m/s)²
  const double V = kL * kL * kL * 1e-30;
  const double expect = kPi / (kKB * kT) * mdot2 / (6 * kC * V * kEps0) / 100 / (2 * kPi * kC * 100);   // cm⁻² (α n in cm⁻¹ over ν̃ in cm⁻¹)
  for (bool vel : {true, false}) {
    const auto r = analyze(oscillating_dipole(vel), {"ir"}, o);
    ASSERT_EQ(r.size(), 1u);
    const Property& p = r[0];
    ASSERT_FALSE(p.series.empty()) << (p.notes.empty() ? "" : p.notes[0]);
    const double a = area(p.series[0]);
    std::printf("%s: strongest band %.2f cm⁻¹, area %.6g (sum rule %.6g) cm⁻², resolution %.2f cm⁻¹\n", vel ? "velocities" : "differences", p.value, a, expect,
                p.extra.at("resolution (cm⁻¹)"));
    EXPECT_NEAR(p.value, kNu0, p.extra.at("resolution (cm⁻¹)") / 2);
    EXPECT_NEAR(a, expect, 2e-3 * expect);
  }
}

// The standard and Schofield spectra are the harmonic one times 2 tanh(x/2)/x and 2 sinh(x/2)/x, x = βħω
TEST(Infrared, QuantumCorrectionFactors) {
  AnalyzeOptions o;
  o.temperature = kT;
  o.frame_ps = kDtFs * 1e-3;
  o.ir_max_wavenumber = 3000;
  const auto r = analyze(oscillating_dipole(true), {"ir"}, o);
  const Property& p = r[0];
  ASSERT_GE(p.series.size(), 3u);
  size_t k = 0;
  while (p.series[0].x[k] < kNu0) ++k;
  const double x = 1.054571817e-34 * 2 * kPi * 2.99792458e10 * p.series[0].x[k] / (1.380649e-23 * kT);
  EXPECT_NEAR(p.series[1].y[k] / p.series[0].y[k], 2 * std::tanh(x / 2) / x, 1e-12);
  EXPECT_NEAR(p.series[2].y[k] / p.series[0].y[k], 2 * std::sinh(x / 2) / x, 1e-12);
  std::printf("at %.1f cm⁻¹ (βħω = %.3f): standard/harmonic %.4f, Schofield/harmonic %.4f\n", p.series[0].x[k], x, p.series[1].y[k] / p.series[0].y[k],
              p.series[2].y[k] / p.series[0].y[k]);
}
