// Polymer statistics: chain-length draws, the terminal copolymer model, stereo sequences, Flory–Huggins, Ewald sizes.
#include <gtest/gtest.h>

#include <cmath>
#include <numeric>

#include "caps/kspace.hpp"
#include "caps/polymer.hpp"
#include "caps/polystats.hpp"

using namespace caps;

TEST(Polystats, SchulzZimmDrawsMatchTheTargetForLargeSamples) {
  const auto L = draw_chain_lengths("schulz-zimm", 40, 1.1, 20000, 7);
  double s1 = 0, s2 = 0;
  for (int n : L) s1 += n, s2 += double(n) * n;
  const double nn = s1 / double(L.size()), pdi = s2 / s1 / nn;
  EXPECT_NEAR(nn, 40, 0.3);
  EXPECT_NEAR(pdi, 1.1, 0.01);
  EXPECT_EQ(L, draw_chain_lengths("schulz-zimm", 40, 1.1, 20000, 7));   // deterministic for a seed
  const auto F = draw_chain_lengths("flory", 50, 2, 20000, 3);
  s1 = s2 = 0;
  for (int n : F) s1 += n, s2 += double(n) * n;
  EXPECT_NEAR(s2 / s1 / (s1 / double(F.size())), 2 - 1.0 / 50, 0.05);
  // the Gamma density integrates to 1
  double area = 0;
  for (double x = 0.05; x < 200; x += 0.1) area += chain_length_pdf("schulz-zimm", 40, 1.1, x) * 0.1;
  EXPECT_NEAR(area, 1, 1e-3);
}

TEST(Polystats, TerminalModelMatchesMayoLewis) {
  const auto m = copolymer_terminal(0.52, 0.46, 0.5);   // styrene / MMA
  EXPECT_NEAR(m.F1, 0.510, 5e-4);
  EXPECT_NEAR(m.azeotrope, 0.529, 5e-4);
  EXPECT_NEAR(m.run_a, 1.52, 5e-3);
  EXPECT_NEAR(m.run_b, 1.46, 5e-3);
  EXPECT_LT(copolymer_terminal(2, 0.5, 0.5).azeotrope, 0);   // one ratio above 1, one below: no azeotrope
  // long chains drawn with the terminal model reach the model's composition
  ChainSpec spec;
  spec.units.resize(2);
  spec.sequence = Sequence::Terminal;
  spec.dp = 40000, spec.r1 = 0.52, spec.r2 = 0.46, spec.weights = {0.5, 0.5};
  const auto seq = chain_sequence(spec, 12);
  EXPECT_NEAR(double(std::count(seq.begin(), seq.end(), 0)) / double(seq.size()), m.F1, 0.01);
}

TEST(Polystats, StereoModelsAndFits) {
  const auto b = stereo_bernoulli(0.5);
  EXPECT_NEAR(b.mm, 0.25, 1e-12);
  EXPECT_NEAR(b.mr, 0.5, 1e-12);
  EXPECT_NEAR(std::accumulate(b.pentads.begin(), b.pentads.end(), 0.0), 1, 1e-12);
  EXPECT_NEAR(b.pentads[0], 0.0625, 1e-12);   // mmmm = Pm⁴
  const auto iso = stereo_bernoulli(0.8);
  double rms = 1;
  EXPECT_NEAR(fit_bernoulli(iso.pentads, &rms).pm, 0.8, 1e-6);
  EXPECT_LT(rms, 1e-8);
  const auto mk = stereo_markov(0.2, 0.6);
  const auto fit = fit_markov(mk.pentads, &rms);
  EXPECT_NEAR(fit.p_mr, 0.2, 1e-5);
  EXPECT_NEAR(fit.p_rm, 0.6, 1e-5);
  const auto c = count_stereo("mmrmr");
  EXPECT_EQ(c.m, 3);
  EXPECT_EQ(c.mm, 1);
  EXPECT_EQ(c.mr, 3);
  EXPECT_EQ(c.pentad_total, 2);
}

TEST(Polystats, FloryHugginsCriticalPointSpinodalAndBinodal) {
  const auto c = blend_critical(100, 200);
  EXPECT_NEAR(c.chi_c, 0.01457, 1e-5);
  EXPECT_NEAR(c.phi_c, 0.5858, 1e-4);
  double lo, hi;
  ASSERT_TRUE(blend_spinodal(100, 200, 0.03, lo, hi));
  EXPECT_NEAR(lo, 0.1857, 1e-4);
  EXPECT_NEAR(hi, 0.8977, 1e-4);
  ASSERT_TRUE(blend_binodal(100, 200, 0.03, lo, hi));
  EXPECT_NEAR(lo, 0.0385, 1e-4);
  EXPECT_NEAR(hi, 0.9932, 1e-4);
  EXPECT_FALSE(blend_binodal(100, 200, 0.01, lo, hi));   // below χc: one phase
  // symmetric blend: binodal symmetric about ½
  ASSERT_TRUE(blend_binodal(50, 50, 0.05, lo, hi));
  EXPECT_NEAR(lo + hi, 1, 1e-9);
}

TEST(Polystats, EwaldAndSolventEstimates) {
  EXPECT_NEAR(ewald_beta(12, 1e-5), 0.2603, 1e-4);
  EXPECT_EQ(pme_mesh_size(45.3, 1.2), 40);
  EXPECT_EQ(pme_mesh_size(16.0, 1.2), 14);   // 2·7
  EXPECT_NEAR(hildebrand_chi(106.3, 18.2, 18.6, 298.15), 0.347, 1e-3);
}
