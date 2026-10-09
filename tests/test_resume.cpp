// Runs stopped from outside and continued (caps/live.hpp request_stop: what SIGUSR1 / SIGTERM do in a cluster job near
// its time limit): molecular dynamics (a checkpoint of the step, the state kept), an equilibration protocol (stage and
// step), a cooling scan (the temperatures done), a tensile pull (the curve and L0) and a recipe (the stages done). Each
// resumed run reaches the same end as one that was never stopped would: the full step count, every temperature, the
// final strain, every stage.
#include <gtest/gtest.h>

#include <cmath>
#include <filesystem>
#include <fstream>

#include "caps/dynamics.hpp"
#include "caps/equilibrate.hpp"
#include "caps/grow.hpp"
#include "caps/live.hpp"
#include "caps/mechanics.hpp"
#include "caps/recipe.hpp"
#include "caps/relax.hpp"

using namespace caps;

namespace {
const System& cell() {
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
struct StopReset { ~StopReset() { request_stop(false); } };
}  // namespace

TEST(Resume, DynamicsStopsWithACheckpointAndContinues) {
  StopReset reset;
  System s = cell();
  DynamicsOptions d;
  d.steps = 400;
  d.thermo_every = 10;
  d.progress = [](const ThermoRow& r) { if (r.step == 150) request_stop(true); return true; };
  int checkpoints = 0;
  int64_t at = -1;
  d.checkpoint = [&](const std::vector<double>&, const std::vector<double>& v, const Cell&, int64_t step) { ++checkpoints, at = step; EXPECT_EQ(v.size(), 3 * s.atoms.size()); };
  try {
    run_dynamics(s, d);
    FAIL() << "the run did not stop";
  } catch (const DynamicsInterrupted& e) {
    EXPECT_EQ(e.step, 150);
    EXPECT_EQ(checkpoints, 1);
    EXPECT_EQ(at, 150);
  }
  EXPECT_EQ(s.velocities.size(), s.atoms.size());   // the state of the stop is kept
  request_stop(false);
  DynamicsOptions r = d;
  r.progress = nullptr;
  r.step_offset = 150;
  r.steps = 250;
  r.new_velocities = false;
  DynamicsReport rep;
  run_dynamics(s, r, &rep);
  EXPECT_EQ(rep.thermo.back().step, 400);
}

TEST(Resume, EquilibrationGoesOnFromItsStageAndStep) {
  StopReset reset;
  EquilibrateOptions e;
  auto make_stage = [](const char* label, double ps, double t0, double t1 = -1) { Stage x; x.label = label, x.ensemble = Ensemble::NVT, x.ps = ps, x.t_start = t0, x.t_end = t1; return x; };
  e.stages = {make_stage("hot", 0.2, 500), make_stage("cool", 0.3, 500, 300), make_stage("hold", 0.1, 300)};
  e.thermo_ps = 0.01;
  int stage = -1;
  int64_t step = -1;
  e.checkpoint = [&](const std::vector<double>&, const std::vector<double>&, const Cell&, int k, int64_t at) { stage = k, step = at; };
  e.progress = [](int k, int, const std::string&, const ThermoRow& r) { if (k == 2 && r.step >= 300) request_stop(true); return true; };
  System s = cell();
  EXPECT_THROW(equilibrate(s, e), DynamicsInterrupted);
  ASSERT_EQ(stage, 1);   // the second stage, 0-based
  EXPECT_GT(step, 0);
  EXPECT_LT(step, 300);
  request_stop(false);
  EquilibrateOptions r = e;
  r.progress = nullptr;
  r.start_stage = stage;
  r.start_step = step;
  EquilibrateReport rep;
  equilibrate(s, r, &rep);
  EXPECT_EQ(rep.thermo.back().step, 600);   // 200 + 300 + 100 steps of 1 fs: the protocol's end
  EXPECT_EQ(rep.stages.size(), 2u);         // the stage continued and the last
  // the ramp continues from where it was: the first row of the continued stage is near its target then
  const double frac = double(step) / 300;
  EXPECT_NEAR(rep.thermo.front().target_temperature, 500 + (300 - 500) * frac, 1.0);
}

TEST(Resume, CoolingScanSkipsTheTemperaturesDone) {
  StopReset reset;
  CoolingOptions c;
  c.t_start = 400, c.t_end = 300, c.t_step = 50, c.ps_per_step = 0.2, c.equilibrate_ps = 0.1;
  std::vector<CoolingPoint> done;
  System state;
  c.after_point = [&](const System& now, const std::vector<CoolingPoint>& pts) { done = pts, state = now; if (pts.size() == 1) request_stop(true); };
  System s = cell();
  EXPECT_THROW(run_cooling(s, c), DynamicsInterrupted);
  ASSERT_EQ(done.size(), 1u);
  EXPECT_DOUBLE_EQ(done[0].temperature, 400);
  request_stop(false);
  CoolingOptions r = c;
  r.after_point = nullptr;
  r.done = done;
  const CoolingResult res = run_cooling(state, r);
  ASSERT_EQ(res.points.size(), 3u);
  EXPECT_DOUBLE_EQ(res.points[0].density, done[0].density);   // kept, not run again
  EXPECT_DOUBLE_EQ(res.points[2].temperature, 300);
}

TEST(Resume, TensilePullContinuesFromTheStrainReached) {
  StopReset reset;
  TensileOptions t;
  t.rate = 0.05, t.max_strain = 0.06, t.equilibrate_ps = 0.05, t.sample_every = 20;
  std::vector<TensilePoint> curve;
  double L0[3] = {0, 0, 0};
  System state;
  System s = cell();
  t.checkpoint = [&](const std::vector<double>& x, const std::vector<double>& v, const Cell& c, const std::vector<TensilePoint>& cu, const double l0[3]) {
    curve = cu;
    for (int k = 0; k < 3; ++k) L0[k] = l0[k];
    state = s;
    state.cell = c;
    state.velocities.resize(state.atoms.size());
    for (size_t i = 0; i < state.atoms.size(); ++i) state.atoms[i].pos = {x[3 * i], x[3 * i + 1], x[3 * i + 2]}, state.velocities[i] = {v[3 * i], v[3 * i + 1], v[3 * i + 2]};
  };
  t.progress = [](const TensilePoint& p) { if (p.strain > 0.025) request_stop(true); return true; };
  EXPECT_THROW(run_tensile(s, t), DynamicsInterrupted);
  ASSERT_FALSE(curve.empty());
  const double e0 = curve.back().strain;
  EXPECT_GT(e0, 0.02);
  request_stop(false);
  TensileResume rs;
  rs.curve = curve;
  for (int k = 0; k < 3; ++k) rs.L0[k] = L0[k];
  TensileOptions r = t;
  r.progress = nullptr;
  r.checkpoint = nullptr;
  r.resume = &rs;
  const TensileResult res = run_tensile(state, r);
  EXPECT_NEAR(res.curve.back().strain, 0.06, 0.002);   // the pull went on to its end at the same rate
  EXPECT_NEAR(state.cell.a[0] / L0[0] - 1, 0.06, 0.002);   // measured from the original L0
  for (size_t i = 1; i < res.curve.size(); ++i) EXPECT_GE(res.curve[i].strain, res.curve[i - 1].strain - 1e-9);
}

TEST(Resume, RecipeGoesOnFromTheStagesDone) {
  StopReset reset;
  const auto dir = std::filesystem::temp_directory_path() / "caps_resume_recipe";
  std::filesystem::remove_all(dir);
  const Json r = Json::parse(R"({"recipe":1,"name":"tiny","build":{"polymer":{"smiles":"*CC*","dp":4,"chains":2}},"grow":{"density":0.5,"seed":1},
                                "relax":{"fmax":5},"md":{"ps":0.3,"temperature":300,"ensemble":"nvt"}})");
  RecipeOptions o;
  o.out_dir = dir.string();
  o.checkpoint_dir = dir.string();
  o.sha256 = "test-recipe";
  o.progress = [](const RecipeEvent& e) { if (e.name == "md" && e.status == "running" && e.fraction > 0.3) request_stop(true); };
  EXPECT_THROW(run_recipe(r, o), RecipeError);
  ASSERT_TRUE(std::filesystem::exists(dir / "recipe.state.json"));
  request_stop(false);
  std::vector<std::string> events;
  RecipeOptions again = o;
  again.resume = true;
  again.progress = [&](const RecipeEvent& e) { if (e.status == "done") events.push_back(e.name + ":" + e.detail); };
  const RecipeResult res = run_recipe(r, again);
  EXPECT_FALSE(res.system.atoms.empty());
  bool from_checkpoint = false, md_done = false;
  for (const auto& e : events) from_checkpoint |= e == "relax:from the checkpoint", md_done |= e.rfind("md:", 0) == 0 && e.find("checkpoint") == std::string::npos;
  EXPECT_TRUE(from_checkpoint);
  EXPECT_TRUE(md_done);
  // another recipe cannot take this checkpoint
  RecipeOptions other = again;
  other.sha256 = "another";
  EXPECT_THROW(run_recipe(r, other), RecipeError);
}
