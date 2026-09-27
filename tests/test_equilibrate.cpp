#include <gtest/gtest.h>

#include <cmath>

#include "caps/analysis.hpp"
#include "caps/equilibrate.hpp"
#include "caps/grow.hpp"
#include "caps/relax.hpp"

using namespace caps;

namespace {

const System& small_relaxed() {
  static const System s = [] {
    GrowOptions g;
    g.chains = 3;
    g.dp = 4;
    g.density = 0.4;
    g.seed = 5;
    System c = grow(g);
    RelaxOptions r;
    r.target_density = 0.9;
    r.ftol = 1.0;
    relax(c, r);
    return c;
  }();
  return s;
}

}  // namespace

TEST(Protocol, Larsen21MatchesThePublishedSchedule) {
  ProtocolParams p;
  p.t_max = 600;
  const auto s = larsen21(p);
  ASSERT_EQ(s.size(), 21u);
  EXPECT_DOUBLE_EQ(protocol_ps(s), 1560.0);
  const double fractions[] = {0.02, 0.6, 1.0, 0.5, 0.1, 0.01};
  int k = 0;
  for (size_t i = 2; i < 20; i += 3) {
    EXPECT_EQ(s[i].ensemble, Ensemble::NPT);
    EXPECT_NEAR(s[i].pressure, fractions[k++] * 49346.2, 1e-6);
    EXPECT_EQ(s[i].t_start, 300);
    EXPECT_EQ(s[i - 2].t_start, 600);   // heat before each pressure step
  }
  EXPECT_EQ(s[20].ensemble, Ensemble::NPT);
  EXPECT_EQ(s[20].ps, 800);
  EXPECT_EQ(s[20].pressure, 1);
  p.time_scale = 0.1;
  EXPECT_NEAR(protocol_ps(larsen21(p)), 156.0, 1e-9);
}

TEST(Protocol, TextRoundTripAndErrors) {
  ProtocolParams p;
  for (const auto& name : {"larsen21", "annealing", "pushoff"}) {
    const auto a = protocol_by_name(name, p);
    const auto b = parse_protocol(protocol_text(a));
    ASSERT_EQ(a.size(), b.size()) << name;
    for (size_t i = 0; i < a.size(); ++i) {
      EXPECT_EQ(a[i].ensemble, b[i].ensemble);
      EXPECT_NEAR(a[i].ps, b[i].ps, 1e-9);
      EXPECT_NEAR(a[i].t_start, b[i].t_start, 1e-9);
      EXPECT_NEAR(a[i].t_end, b[i].t_end, 1e-9);
      EXPECT_NEAR(a[i].pressure, b[i].pressure, 0.006);   // the text keeps 0.01 atm
      EXPECT_NEAR(a[i].force_cap, b[i].force_cap, 1e-9);
      EXPECT_EQ(a[i].label, b[i].label);
    }
  }
  const auto s = parse_protocol("NPT 2 ns T 450 P 1.01325 bar  # warm\nnvt 500 fs T 300 to 350\n\n# comment only\n");
  ASSERT_EQ(s.size(), 2u);
  EXPECT_EQ(s[0].ps, 2000);
  EXPECT_NEAR(s[0].pressure, 1.0, 1e-9);
  EXPECT_EQ(s[0].label, "warm");
  EXPECT_EQ(s[1].t_end, 350);
  EXPECT_NEAR(s[1].ps, 0.5, 1e-12);
  EXPECT_THROW(parse_protocol("npt 5 ps T 300"), std::invalid_argument);         // no pressure
  EXPECT_THROW(parse_protocol("nvt T 300"), std::invalid_argument);              // no duration
  EXPECT_THROW(parse_protocol("md 5 ps"), std::invalid_argument);                // unknown ensemble
  EXPECT_THROW(parse_protocol("nvt 5 ps T 300 Q 2"), std::invalid_argument);     // unknown keyword
  EXPECT_THROW(parse_protocol(""), std::invalid_argument);
}

TEST(Protocol, AnnealingCycles) {
  ProtocolParams p;
  p.cycles = 2;
  const auto s = annealing(p);
  ASSERT_EQ(s.size(), 9u);
  EXPECT_EQ(s[1].t_start, 300);
  EXPECT_EQ(s[1].t_end, 600);
  EXPECT_EQ(s[3].t_end, 300);
  for (const auto& st : s) EXPECT_EQ(st.ensemble, Ensemble::NPT);
}

TEST(Equilibrate, RunsStagesContinuously) {
  System s = small_relaxed();
  EquilibrateOptions o;
  o.stages = parse_protocol("nvt 1 ps T 400\nnpt 1 ps T 300 P 1 atm\nnvt 1 ps T 300 to 500 # ramp\nnvt 0.5 ps T 300 cap 20");
  o.thermo_ps = 0.1;
  o.frame_ps = 0.5;
  int frames = 0;
  o.frame = [&](const std::vector<double>&, const Cell&, int64_t) { ++frames; };
  EquilibrateReport r;
  equilibrate(s, o, &r);
  ASSERT_EQ(r.stages.size(), 4u);
  EXPECT_NEAR(r.ps, 3.5, 1e-9);
  for (size_t k = 1; k < r.thermo.size(); ++k) ASSERT_GT(r.thermo[k].step, r.thermo[k - 1].step);
  EXPECT_EQ(r.thermo.back().step, 3500);
  // the ramp ends at 500 K
  bool seen = false;
  for (const auto& row : r.thermo)
    if (row.step == 3000) { EXPECT_NEAR(row.target_temperature, 500, 1e-9); seen = true; }
  EXPECT_TRUE(seen);
  EXPECT_EQ(frames, 8);   // start, then every 0.5 ps and stage ends without repeats
  EXPECT_EQ(s.velocities.size(), s.atoms.size());
  EXPECT_EQ(r.stages[2].label, "ramp");
}

TEST(Equilibrate, PushoffTamesOverlaps) {
  GrowOptions g;
  g.chains = 3;
  g.dp = 4;
  g.density = 0.6;
  g.contact_scale = 0.75;
  g.seed = 2;
  System s = grow(g);
  RelaxOptions r;
  r.pushoff = false;
  r.max_iterations = 50;   // barely minimised: close contacts remain
  relax(s, r);
  ProtocolParams p;
  p.time_scale = 0.1;
  EquilibrateOptions o;
  o.stages = pushoff(p);
  o.md.dt = 0.5;
  EquilibrateReport rep;
  equilibrate(s, o, &rep);
  EXPECT_TRUE(std::isfinite(rep.thermo.back().potential));
  EXPECT_LT(std::fabs(rep.stages.back().temperature - 300), 60);
}

TEST(Equilibrate, ConvergenceChecks) {
  {
    System s = small_relaxed();
    EquilibrateOptions o;
    o.stages = parse_protocol("npt 2 ps T 300 P 1 atm");
    o.until_converged = true;
    o.block_ps = 1;
    o.min_blocks = 3;
    o.max_blocks = 6;
    o.tol_density = 0.2;
    o.tol_energy = 1.0;
    o.tol_rg = 0.5;
    o.tol_internal = 0.5;
    EquilibrateReport r;
    equilibrate(s, o, &r);
    EXPECT_TRUE(r.converged);
    EXPECT_EQ(r.blocks, 3);
    ASSERT_EQ(r.checks.size(), 4u);
    EXPECT_EQ(r.checks[3].quantity, "internal distances");
    for (const auto& c : r.checks) EXPECT_TRUE(c.ok) << c.quantity;
  }
  {
    // the same, but the internal distances must match a target far from the chains' own: never converged
    System s = small_relaxed();
    EquilibrateOptions o;
    o.stages = parse_protocol("npt 2 ps T 300 P 1 atm");
    o.until_converged = true;
    o.block_ps = 1;
    o.min_blocks = 3;
    o.max_blocks = 4;
    o.tol_density = 0.2;
    o.tol_energy = 1.0;
    o.tol_rg = 0.5;
    o.tol_internal = 0.05;
    o.internal_target.assign(50, 25.0);   // C_n = 25 everywhere: no real chain
    EquilibrateReport r;
    equilibrate(s, o, &r);
    EXPECT_FALSE(r.converged);
    ASSERT_EQ(r.checks.size(), 4u);
    EXPECT_FALSE(r.checks[3].ok);
    EXPECT_GT(r.checks[3].change, 0.5);
  }
  {
    System s = small_relaxed();
    EquilibrateOptions o;
    o.stages = parse_protocol("nvt 1 ps T 300");
    o.until_converged = true;
    o.block_ps = 0.5;
    o.max_blocks = 3;
    o.tol_density = 1e-12;
    o.tol_energy = 1e-12;
    o.tol_rg = 1e-12;
    EquilibrateReport r;
    equilibrate(s, o, &r);
    EXPECT_FALSE(r.converged);
    EXPECT_EQ(r.blocks, 3);
    EXPECT_EQ(r.checks[0].blocks.size(), 3u);
  }
}

TEST(Chains, BackbonesAndInternalDistances) {
  GrowOptions g;
  g.chains = 2;
  g.dp = 10;
  g.density = 0.2;
  g.seed = 4;
  const System s = grow(g);
  const auto bb = backbones(s);
  ASSERT_EQ(bb.size(), 2u);
  EXPECT_EQ(bb[0].size(), 20u);   // two backbone carbons per styrene unit
  for (uint32_t i : bb[0]) EXPECT_EQ(s.atoms[i].element, 6);
  const auto d = internal_distances(s);
  EXPECT_EQ(d.chains, 2);
  ASSERT_FALSE(d.ratio.empty());
  EXPECT_NEAR(d.ratio[0], 1.0, 1e-12);
  EXPECT_NEAR(std::sqrt(d.b2), 1.54, 0.02);
  EXPECT_GT(d.ratio.back(), 1.5);
}
