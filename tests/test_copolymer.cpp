#include <gtest/gtest.h>

#include <array>
#include <cmath>

#include "caps/polymer.hpp"
#include "caps/polystats.hpp"

using namespace caps;

// Two units: the stationary composition of the terminal chain is the Mayo–Lewis composition, and the mean run lengths
// are 1 / (1 − P_ii)
TEST(Copolymer, TerminalModelIsMayoLewis) {
  for (double r1 : {0.1, 0.5, 2.0, 8.0})
    for (double r2 : {0.05, 0.4, 1.0, 3.0})
      for (double f1 : {0.1, 0.5, 0.85}) {
        const auto m = terminal_model({f1, 1 - f1}, {1, r1, r2, 1});
        EXPECT_NEAR(m.composition[0], mayo_lewis(r1, r2, f1), 1e-12);
        const auto c = copolymer_terminal(r1, r2, f1);
        EXPECT_NEAR(m.run_length[0], c.run_a, 1e-9 * c.run_a);
        EXPECT_NEAR(m.run_length[1], c.run_b, 1e-9 * c.run_b);
      }
}

// Three units: the composition equals the Alfrey–Goldfinger equation (J. Chem. Phys. 12, 322 (1944))
TEST(Copolymer, TernaryIsAlfreyGoldfinger) {
  const double r12 = 0.4, r13 = 1.8, r21 = 0.7, r23 = 0.25, r31 = 3.0, r32 = 0.6;
  for (const auto& M : std::vector<std::array<double, 3>>{{0.2, 0.3, 0.5}, {0.6, 0.1, 0.3}, {1, 1, 1}}) {
    const double m1 = M[0], m2 = M[1], m3 = M[2];
    const double d1 = m1 * (m1 / (r31 * r21) + m2 / (r21 * r32) + m3 / (r31 * r23)) * (m1 + m2 / r12 + m3 / r13);
    const double d2 = m2 * (m1 / (r12 * r31) + m2 / (r12 * r32) + m3 / (r32 * r13)) * (m2 + m1 / r21 + m3 / r23);
    const double d3 = m3 * (m1 / (r13 * r21) + m2 / (r23 * r12) + m3 / (r13 * r23)) * (m3 + m1 / r31 + m2 / r32);
    const double s = d1 + d2 + d3;
    const auto m = terminal_model({m1, m2, m3}, {1, r12, r13, r21, 1, r23, r31, r32, 1});
    EXPECT_NEAR(m.composition[0], d1 / s, 1e-12);
    EXPECT_NEAR(m.composition[1], d2 / s, 1e-12);
    EXPECT_NEAR(m.composition[2], d3 / s, 1e-12);
  }
}

// Penultimate model: the composition of Merz, Alfrey & Goldfinger (J. Polym. Sci. 1, 75 (1946)), x = [M1]/[M2]:
// d[M1]/d[M2] = (1 + r1' x (r1 x + 1)/(r1' x + 1)) / (1 + r2' (r2 + x)/(x (r2' + x))); equal ratios give the terminal model
TEST(Copolymer, PenultimateComposition) {
  const double r1 = 0.3, r1p = 1.5, r2 = 0.8, r2p = 0.1;
  for (double f1 : {0.2, 0.5, 0.7}) {
    const double x = f1 / (1 - f1);
    const double ratio = (1 + r1p * x * (r1 * x + 1) / (r1p * x + 1)) / (1 + r2p * (r2 + x) / (x * (r2p + x)));
    const auto m = penultimate_model(f1, r1, r1p, r2, r2p);
    EXPECT_NEAR(m.composition[0], ratio / (1 + ratio), 1e-12);
    const auto t = penultimate_model(f1, r1, r1, r2, r2);
    EXPECT_NEAR(t.composition[0], mayo_lewis(r1, r2, f1), 1e-12);
    const auto tt = terminal_model({f1, 1 - f1}, {1, r1, r2, 1});
    for (size_t k = 0; k < 8; ++k) EXPECT_NEAR(t.triads[k], tt.triads[k], 1e-12);
  }
}

// Drift in a batch: the integrated feed equals the closed form of Meyer & Lowry (J. Polym. Sci. A 3, 2843 (1965)),
// X = 1 − (f1/f1⁰)^α (f2/f2⁰)^β ((f1⁰ − δ)/(f1 − δ))^γ, α = r2/(1 − r2), β = r1/(1 − r1),
// γ = (1 − r1 r2)/((1 − r1)(1 − r2)), δ = (1 − r2)/(2 − r1 − r2); the cumulative composition closes the mass balance
TEST(Copolymer, DriftIsMeyerLowry) {
  const double r1 = 0.5, r2 = 0.2, f10 = 0.6;
  const double a = r2 / (1 - r2), b = r1 / (1 - r1), g = (1 - r1 * r2) / ((1 - r1) * (1 - r2)), d = (1 - r2) / (2 - r1 - r2);
  const auto pts = composition_drift({f10, 1 - f10}, {1, r1, r2, 1}, 0.9, 46);
  ASSERT_EQ(pts.size(), 46u);
  for (const auto& p : pts) {
    const double f1 = p.feed[0];
    const double X = 1 - std::pow(f1 / f10, a) * std::pow((1 - f1) / (1 - f10), b) * std::pow((f10 - d) / (f1 - d), g);
    EXPECT_NEAR(X, p.conversion, 1e-8) << "f1 " << f1;
    if (p.conversion > 0) EXPECT_NEAR(p.conversion * p.cumulative[0] + (1 - p.conversion) * f1, f10, 1e-12);
  }
  // below the azeotrope δ = 0.615 the copolymer is richer in A than the feed, which loses A
  EXPECT_LT(pts.back().feed[0], f10);
}

// Built chains carry the model: composition, dyads, triads and run lengths of 300 chains of 300 units (three units)
// within their sampling error; the given matrix and the penultimate model likewise; "along" drift makes the head richer
// in the reactive unit than the tail
TEST(Copolymer, BuiltChainsFollowTheModel) {
  ChainSpec spec;
  spec.units.resize(3);
  spec.sequence = Sequence::Terminal;
  spec.dp = 300;
  spec.weights = {0.3, 0.3, 0.4};
  spec.reactivity = {1, 0.4, 1.8, 0.7, 1, 0.25, 3.0, 0.6, 1};
  std::vector<std::vector<int>> chains;
  for (int c = 0; c < 300; ++c) chains.push_back(chain_sequence(spec, 1 + 101 * uint64_t(c)));
  const auto st = sequence_statistics(chains, 3);
  const auto m = terminal_model(spec.weights, spec.reactivity);
  for (int i = 0; i < 3; ++i) {
    EXPECT_NEAR(st.composition[size_t(i)], m.composition[size_t(i)], 0.008);
    EXPECT_NEAR(st.run_length[size_t(i)], m.run_length[size_t(i)], 0.04 * m.run_length[size_t(i)]);
  }
  for (size_t k = 0; k < 9; ++k) EXPECT_NEAR(st.dyads[k], m.dyads[k], 0.006);
  for (size_t k = 0; k < 27; ++k) EXPECT_NEAR(st.triads[k], m.triads[k], 0.005);
  ChainSpec mk = spec;
  mk.sequence = Sequence::Markov;
  mk.transition = {0.8, 0.1, 0.1, 0.2, 0.5, 0.3, 0.05, 0.05, 0.9};
  chains.clear();
  for (int c = 0; c < 300; ++c) chains.push_back(chain_sequence(mk, 7 + 101 * uint64_t(c)));
  const auto sm = sequence_statistics(chains, 3);
  const auto tm = transition_model(mk.transition, 3);
  for (size_t k = 0; k < 9; ++k) EXPECT_NEAR(sm.dyads[k], tm.dyads[k], 0.01);
  ChainSpec pe;
  pe.units.resize(2);
  pe.sequence = Sequence::Penultimate;
  pe.dp = 300;
  pe.weights = {0.5, 0.5};
  pe.penultimate = {0.3, 1.5, 0.8, 0.1};
  chains.clear();
  for (int c = 0; c < 300; ++c) chains.push_back(chain_sequence(pe, 3 + 101 * uint64_t(c)));
  const auto sp = sequence_statistics(chains, 2);
  const auto pm = penultimate_model(0.5, 0.3, 1.5, 0.8, 0.1);
  for (size_t k = 0; k < 8; ++k) EXPECT_NEAR(sp.triads[k], pm.triads[k], 0.006);
  // drift along the chain: r1 = 5, r2 = 0.2 from an equal feed to 90 % conversion; A is used up early
  ChainSpec dr;
  dr.units.resize(2);
  dr.sequence = Sequence::Terminal;
  dr.dp = 200;
  dr.weights = {0.5, 0.5};
  dr.r1 = 5, dr.r2 = 0.2;
  dr.drift = "along", dr.conversion = 0.9;
  double head = 0, tail = 0;
  for (int c = 0; c < 200; ++c) {
    const auto s = chain_sequence(dr, 11 + 101 * uint64_t(c));
    for (int i = 0; i < 20; ++i) head += s[size_t(i)] == 0, tail += s[size_t(180 + i)] == 0;
  }
  head /= 4000, tail /= 4000;
  const auto table = composition_drift({0.5, 0.5}, {1, 5, 0.2, 1}, 0.9, 201);
  EXPECT_NEAR(head, table.front().instantaneous[0], 0.04);
  EXPECT_LT(tail, 0.5 * head);
  std::printf("drift along: A at the head %.3f (model %.3f), at the tail %.3f\n", head, table.front().instantaneous[0], tail);
}
