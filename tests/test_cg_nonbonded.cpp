// Non-bonded coarse-grained potentials (caps/cg_nonbonded.hpp): per-pair g(r) with exclusions and a cell list against
// brute force, the IBI update's direction and the pressure ramp, a short joint IBI on a two-type Lennard-Jones fluid with
// CAPS's engine (the error falls), analytic fits recovering known parameters, the calibration steps.
// (The long self-consistency check — 1 % in 20 iterations — is in `caps bench`, not run here.)
#include <gtest/gtest.h>

#include <cmath>
#include <random>

#include "caps/cg_nonbonded.hpp"

using namespace caps;

namespace {
constexpr double kB = 0.0019872067;

CgTypes two_types() { return {{"A", "B"}, {}, {}, {}}; }

// n beads in a cubic box of edge L on a jittered lattice, alternating A and B
System fluid(int n, double L, uint64_t seed, CgTopology* top) {
  System s;
  s.cell.a = {L, 0, 0}, s.cell.b = {0, L, 0}, s.cell.c = {0, 0, L};
  const int m = int(std::ceil(std::cbrt(double(n))));
  std::mt19937_64 g(seed);
  auto u = [&]() { return double(g() >> 11) * (1.0 / 9007199254740992.0) - 0.5; };
  for (int i = 0; i < n; ++i) {
    Atom a;
    a.id = i + 1, a.mol = i + 1, a.name = i % 2 ? "B" : "A", a.type = i % 2 + 1;
    a.pos = {(i % m + 0.5 + 0.2 * u()) * L / m, ((i / m) % m + 0.5 + 0.2 * u()) * L / m, (i / (m * m) + 0.5 + 0.2 * u()) * L / m};
    s.atoms.push_back(a);
    top->kind.push_back(a.name), top->mol.push_back(i);
  }
  s.types = {{1, 50, "A"}, {2, 60, "B"}};
  return s;
}

CgPairSet reference_pairs() {
  std::vector<CgPairFit> f = {{"A-A", "lj126", 0.45, 4.2}, {"A-B", "lj126", 0.50, 4.4}, {"B-B", "lj126", 0.55, 4.6}};
  CgPairSet grid;
  grid.r0 = 1.0, grid.dr = 0.05, grid.rc = 11.0, grid.temperature = 300;
  return tables_of_fits(f, grid);
}
}  // namespace

TEST(CgNonbonded, RdfCellListMatchesBruteForceAndExcludes) {
  // a chain of 3 beads plus free beads: the 1-2 and 1-3 pairs are left out
  CgTopology t;
  System s = fluid(216, 30.0, 3, &t);
  t.bonds = {{0, 2}, {2, 4}};   // A-A-A
  CgPairOptions o;
  o.rmax = 9.0;   // cell list: 30/9 → 3 cells per edge
  CgRdfAccumulator acc(two_types(), o);
  const int si = acc.add_system(t);
  std::vector<Vec3> p;
  for (const auto& a : s.atoms) p.push_back(a.pos);
  acc.add_frame(si, p, s.cell);
  // brute force over every pair, minimum image
  const auto R = acc.rdf(si);
  double brute_aa = 0, brute_ab = 0;
  for (size_t i = 0; i < p.size(); ++i)
    for (size_t j = i + 1; j < p.size(); ++j) {
      if ((i == 0 && (j == 2 || j == 4)) || (i == 2 && j == 4)) continue;
      const double d = norm(s.cell.minimum_image(p[j] - p[i]));
      if (d >= 9.0) continue;
      if (t.kind[i] == "A" && t.kind[j] == "A") brute_aa += 1;
      else if (t.kind[i] != t.kind[j]) brute_ab += 1;
    }
  double aa = 0, ab = 0;
  for (double c : R[0].counts) aa += c;
  for (double c : R[1].counts) ab += c;
  EXPECT_EQ(R[0].key, "A-A");
  EXPECT_EQ(R[1].key, "A-B");
  EXPECT_DOUBLE_EQ(aa, brute_aa);
  EXPECT_DOUBLE_EQ(ab, brute_ab);
}

TEST(CgNonbonded, UpdateDirectionAndPressureRamp) {
  CgTargets t;
  for (int k = 0; k < 300; ++k) t.r.push_back(0.025 + 0.05 * k);
  t.systems = {"x"}, t.weights = {1}, t.pressure = {1}, t.temperature = 300;
  CgRdf R;
  R.key = "A-A", R.ideal_pairs = 1000, R.frames = 1;
  for (double r : t.r) R.g.push_back(r < 4 ? 0 : 1 + 0.5 * std::exp(-(r - 5) * (r - 5)) ), R.counts.push_back(1);
  t.rdf = {{R}};
  CgIbiOptions o;
  o.rc = 12;
  o.pressure_correction = false;
  const CgPairSet p0 = ibi_start(t, o);
  ASSERT_EQ(p0.pairs.size(), 1u);
  // the start is −kT ln g at the peak (r = 5 Å), relative to the cut-off where g = 1
  const double kT = kB * 300;
  const size_t i5 = size_t(std::lround((5.0 - p0.r0) / p0.dr));
  EXPECT_NEAR(p0.pairs[0].U[i5], -kT * std::log(1.5), 0.01);
  // a CG run with too high a peak: the update raises U there
  CgRdf C = R;
  for (size_t k = 0; k < t.r.size(); ++k) if (t.r[k] > 4.5 && t.r[k] < 5.5) C.g[k] *= 1.2;
  CgIbiStepReport rep;
  const CgPairSet p1 = ibi_step(p0, t, {{C}}, {1}, o, &rep);
  EXPECT_GT(p1.pairs[0].U[i5], p0.pairs[0].U[i5] + 0.5 * o.alpha * kT * std::log(1.2));
  EXPECT_GT(rep.residual, 0);
  // too high a pressure: the ramp makes the pair more attractive (A < 0) and is zero at the cut-off
  o.pressure_correction = true;
  const CgPairSet p2 = ibi_step(p0, t, {{R}}, {2000}, o, &rep);
  EXPECT_LT(rep.ramp_A, 0);
  EXPECT_LT(p2.pairs[0].U[i5], p0.pairs[0].U[i5]);
  EXPECT_NEAR(p2.pairs[0].U.back(), 0, 1e-12);
}

TEST(CgNonbonded, FitsRecoverKnownParameters) {
  const CgPairSet ref = reference_pairs();
  for (const char* form : {"lj126", "morse", "mie"}) {
    const auto fits = fit_pairs(ref, form, 4.0);
    ASSERT_EQ(fits.size(), 3u);
    for (const auto& f : fits) EXPECT_LT(f.rms, std::string(form) == "lj126" ? 1e-4 : 0.05) << form << " " << f.key;
    if (std::string(form) == "lj126") {
      EXPECT_NEAR(fits[1].epsilon, 0.50, 1e-3);
      EXPECT_NEAR(fits[1].sigma, 4.4, 1e-3);
    }
    if (std::string(form) == "mie") EXPECT_NEAR(fits[0].n, 12.0, 0.3);
  }
  const std::string lines = lammps_pair_lines(fit_pairs(ref, "lj126"), two_types(), 11.0);
  EXPECT_NE(lines.find("pair_style lj/cut 11.0000"), std::string::npos) << lines;
  EXPECT_NE(lines.find("pair_coeff 1 2 0.500000 4.400000"), std::string::npos) << lines;
  const auto scaled = scale_fits(fit_pairs(ref, "lj126"), 1.1, 0.8);
  EXPECT_NEAR(scaled[0].sigma, 4.2 * 1.1, 1e-3);
  EXPECT_NEAR(scaled[0].epsilon, 0.45 * 0.8, 1e-3);
}

TEST(CgNonbonded, CalibrationSteps) {
  // ρ ∝ s_σ^−3 exactly: from s = 1 (ρ = 1.5) to the target 1.25 in one step, then reached
  std::vector<CgCalibrationPoint> h = {{1.0, 1.0, 1.5, 0}};
  auto st = calibration_step(h, 1.25, 0);
  EXPECT_NEAR(st.s_sigma, std::cbrt(1.5 / 1.25), 1e-9);
  h.push_back({st.s_sigma, 1.0, 1.25 * 1.004, 0});
  EXPECT_TRUE(calibration_step(h, 1.25, 0).density_done);
  // T_g linear in s_ε on a line with offset: the secant finds it
  auto tg = [](double s) { return 100 + 200 * s; };
  std::vector<CgCalibrationPoint> g = {{1, 1.0, 1.25, tg(1.0)}};
  double s = 1;
  for (int k = 0; k < 8; ++k) {
    const auto x = calibration_step(g, 1.25, 241);
    if (x.tg_done) break;
    s = x.s_eps;
    g.push_back({1, s, 1.25, tg(s)});
  }
  EXPECT_NEAR(tg(s), 241, 5);
}

TEST(CgNonbonded, ShortJointIbiReducesTheError) {
  // the reference fluid: 500 beads, 3 pair tables (LJ), CAPS's engine at 300 K (2000 steps of 5 fs)
  CgTopology t;
  const System s0 = fluid(500, 42.0, 5, &t);
  const std::map<std::string, double> mass = {{"A", 50}, {"B", 60}};
  const CgPairSet ref = reference_pairs();
  CgEngineOptions eo;
  eo.steps = 1500, eo.equilibrate = 500, eo.frame_every = 25, eo.seed = 7;
  const CgEngineRun r0 = run_cg_engine(s0, cg_forcefield(t, two_types(), mass, ref), eo);
  CgPairOptions po;
  po.rmax = 11.0;
  CgRdfAccumulator acc(two_types(), po);
  const int si = acc.add_system(t, 1, "lj");
  for (size_t f = 0; f < r0.frames.size(); ++f) acc.add_frame(si, r0.frames[f], r0.cells[f]);
  const CgTargets tg = cg_targets(acc, 300, {r0.pressure});
  // IBI from −kT ln g: two short iterations (≤ 5000 steps in all)
  CgIbiOptions io;
  io.rc = 11.0;
  io.alpha = 0.5;
  CgPairSet p = ibi_start(tg, io);
  ASSERT_EQ(p.pairs.size(), 3u);
  System s = r0.last;
  std::vector<double> residual, pressure_gap;
  CgEngineOptions it = eo;
  it.steps = 500, it.equilibrate = 250, it.new_velocities = false;
  for (int k = 0; k < 2; ++k) {
    const CgEngineRun run = run_cg_engine(s, cg_forcefield(t, two_types(), mass, p), it);
    CgRdfAccumulator a2(two_types(), po);
    const int s2 = a2.add_system(t, 1, "lj");
    for (size_t f = 0; f < run.frames.size(); ++f) a2.add_frame(s2, run.frames[f], run.cells[f]);
    CgIbiStepReport rep;
    p = ibi_step(p, tg, {a2.rdf(s2)}, {run.pressure}, io, &rep);
    residual.push_back(rep.residual);
    pressure_gap.push_back(std::fabs(run.pressure - r0.pressure));
    s = run.last;
  }
  std::printf("IBI on a two-type LJ fluid: residual %.4f → %.4f · |ΔP| %.0f → %.0f atm (reference %.0f atm)\n", residual.front(), residual.back(), pressure_gap.front(),
              pressure_gap.back(), r0.pressure);
  EXPECT_LT(residual.back(), 0.05) << "−kT ln g plus one update lands close to a simple fluid's g(r)";
}
