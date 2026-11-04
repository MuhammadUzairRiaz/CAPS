#include <gtest/gtest.h>

#include <cmath>
#include <cstdio>

#include "caps/ddft.hpp"

using namespace caps;

namespace {
constexpr double kPi = 3.14159265358979323846;
// the discrete Gaussian chain's correlation (1/N) Σ_{s∈I, s'∈J} K^|s−s'|, K = exp(−q²a²/6)
double S(const std::string& seq, char I, char J, double q, double a) {
  const double K = std::exp(-q * q * a * a / 6);
  double s = 0;
  for (size_t i = 0; i < seq.size(); ++i)
    for (size_t j = 0; j < seq.size(); ++j)
      if (seq[i] == I && seq[j] == J) s += std::pow(K, std::fabs(double(i) - double(j)));
  return s / double(seq.size());
}
}  // namespace

// The chain functional: a weak cosine potential on a homopolymer melt gives δθ = −ε S_D(q) cos(qx), S_D the discrete
// chain's Debye sum (1/N) Σ K^|s−s'|
TEST(Ddft, ChainResponseIsTheDebyeSum) {
  DdftOptions o;
  o.molecules = {{"A", "A8", 1.0}};
  o.grid = {32, 1, 1};
  const double eps = 1e-4, q = 2 * kPi * 3 / 32;
  std::vector<double> U(32);
  for (int i = 0; i < 32; ++i) U[size_t(i)] = eps * std::cos(q * i);
  const auto th = ddft_chain_densities(o, U);
  const double want = -eps * S("AAAAAAAA", 'A', 'A', q, o.bond);
  double got = 0;
  for (int i = 0; i < 32; ++i) got += (th[size_t(i)] - 1.0) * std::cos(q * i);
  got *= 2.0 / 32;
  std::printf("Debye response %.6e (exact %.6e)\n", got, want);
  EXPECT_NEAR(got, want, 2e-3 * std::fabs(want));
}

// A symmetric blend A8 + B8 (nearly incompressible): below the RPA spinodal at the box's longest wave every ripple dies,
// well above it (3 χ_s: the one unstable wave grows as e^(Γτ), Γ ≈ q² θ̄ (χK − 2/S)) it grows
TEST(Ddft, BlendSpinodal) {
  const double a = 1.1543, q = 2 * kPi / 16;
  const double Sq = S("AAAAAAAA", 'A', 'A', q, a);   // per species, φ = ½: φ S_D(q) each
  const double K = std::exp(-q * q * a * a / 6);
  const double chis = 1 / (0.5 * Sq) / K / 2 * 2;    // 2χ K = 1/(φ_A S) + 1/(φ_B S)
  auto run = [&](double chi) {
    DdftOptions o;
    o.molecules = {{"A", "A8", 0.5}, {"B", "B8", 0.5}};
    o.chi = {{"AB", chi}};
    o.grid = {16, 16, 1};
    o.kappa = 50, o.noise = 0, o.dt = 0.05, o.steps = 2000, o.frame_every = 2000;
    const auto r = ddft_run(o);
    return std::make_pair(r.order.front()[0], r.order.back()[0]);
  };
  const auto lo = run(0.7 * chis), hi = run(3.0 * chis);
  std::printf("χ_s %.4f · at 0.7 χ_s: %.3e → %.3e · at 3 χ_s: %.3e → %.3e\n", chis, lo.first, lo.second, hi.first, hi.second);
  EXPECT_LT(lo.second, lo.first);
  EXPECT_GT(hi.second, 10 * hi.first);
}

// A symmetric diblock A4B4 above its spinodal starts to form lamellae at the RPA wavelength 2π/q* (the early, linear
// stage; they coarsen later), q* the maximum of the incompressible structure factor (S_AA S_BB − S_AB²)/(S_AA + S_BB + 2 S_AB)
TEST(Ddft, DiblockLamellaeAtTheRpaSpacing) {
  const double a = 1.1543;
  double best = 0, qs = 0;
  for (double q = 0.05; q < 3; q += 0.001) {
    const double aa = S("AAAABBBB", 'A', 'A', q, a), bb = S("AAAABBBB", 'B', 'B', q, a), ab = S("AAAABBBB", 'A', 'B', q, a);
    const double s = (aa * bb - ab * ab) / (aa + bb + 2 * ab);
    if (s > best) best = s, qs = q;
  }
  DdftOptions o;
  o.molecules = {{"AB", "A4B4", 1.0}};
  o.chi = {{"AB", 1.6 / (2 * best)}};   // 1.6 × χ_s
  o.grid = {32, 32, 1};
  o.kappa = 20, o.noise = 200, o.dt = 0.01, o.steps = 300, o.frame_every = 300;
  const auto r = ddft_run(o);
  const double want = 2 * kPi / qs, got = r.dominant_spacing.back();
  std::printf("A4B4: χ_s %.3f, RPA spacing %.2f · lamellae %.2f · order %.4f → %.4f · %d inversion iterations at most · %.1f s\n", 1 / (2 * best), want, got,
              r.order.front()[0], r.order.back()[0], r.worst_inverse, r.seconds);
  EXPECT_NEAR(got, want, 0.12 * want);
  EXPECT_GT(r.order.back()[0], 1e-3);
}
