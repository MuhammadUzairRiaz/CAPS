// Live output of long runs (caps/live.hpp): the thermo CSV's fixed columns, progress lines read while a run writes them
// (a line cut off mid-write is left for the next read), frames appended one by one and read back as a LAMMPS dump, and
// the thread count a cluster job gives (CAPS_THREADS, then SLURM_CPUS_PER_TASK).
#include <gtest/gtest.h>

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <sstream>

#include "caps/io.hpp"
#include "caps/live.hpp"

using namespace caps;

namespace {
std::string slurp(const std::string& p) {
  std::ifstream f(p);
  std::stringstream ss;
  ss << f.rdbuf();
  return ss.str();
}
void set_env(const char* k, const char* v) {
#ifdef _WIN32
  _putenv_s(k, v ? v : "");
#else
  if (v) ::setenv(k, v, 1);
  else ::unsetenv(k);
#endif
}
std::filesystem::path tmpdir(const char* name) {
  auto d = std::filesystem::temp_directory_path() / name;
  std::filesystem::remove_all(d);
  std::filesystem::create_directories(d);
  return d;
}
}  // namespace

TEST(Live, ThermoColumnsFixedAndRowsMatch) {
  const std::string h = thermo_header();
  const auto cols = std::count(h.begin(), h.end(), ',') + 1;
  ThermoRow r;
  r.step = 500, r.time_ps = 0.5, r.temperature = 301.25, r.pressure = -12.5, r.density = 1.0123, r.lx = 30;
  const std::string row = thermo_csv_row(r, "nvt, hot");
  EXPECT_EQ(std::count(row.begin(), row.end(), ',') + 1, cols);   // a comma in the stage name does not add a column
  EXPECT_EQ(row.rfind("500,0.5000,nvt; hot,", 0), 0u);
  EXPECT_NE(h.find("density_g_cm3"), std::string::npos);
}

TEST(Live, ProgressLinesReadWhileWritten) {
  const auto d = tmpdir("caps_live_progress");
  const std::string f = (d / "progress.jsonl").string();
  LiveOutput live("md");
  live.open_progress(f);
  ThermoRow r;
  r.step = 100, r.temperature = 300;
  live.progress("md", 0, 0, 0.1, &r);
  r.step = 200;
  Json extra = Json::object();
  extra["stress_MPa"] = 12.5;
  live.progress("md", 0, 0, 0.2, &r, extra);
  std::string text = slurp(f);
  text += "{\"command\":\"md\",\"step\":3";   // a line still being written
  size_t used = 0;
  const auto lines = read_progress_lines(text, &used);
  ASSERT_EQ(lines.size(), 3u);   // the start line and two reports
  EXPECT_EQ(lines[0]["state"].str(), "running");
  EXPECT_EQ(lines[2]["step"].number(), 200);
  EXPECT_EQ(lines[2]["stress_MPa"].number(), 12.5);
  EXPECT_EQ(used, slurp(f).size());   // the partial line is left for the next read
  live.done(true);
  const auto all = read_progress_lines(slurp(f));
  EXPECT_EQ(all.back()["state"].str(), "finished");
}

TEST(Live, FramesAppendedReadBackAsDump) {
  const auto d = tmpdir("caps_live_frames");
  System s;
  s.cell.a = {20, 0, 0}, s.cell.b = {0, 20, 0}, s.cell.c = {0, 0, 20};
  for (int i = 0; i < 3; ++i) {
    Atom a;
    a.id = i + 1, a.mol = 1, a.type = 1, a.element = 6, a.pos = {1.0 + i, 2, 3};
    s.atoms.push_back(a);
  }
  const std::string f = (d / "traj.lammpstrj").string();
  LiveOutput live;
  live.open_frames(f, s);
  for (int k = 0; k < 4; ++k) {
    std::vector<double> x;
    for (int i = 0; i < 3; ++i) x.insert(x.end(), {1.0 + i + 0.1 * k, 2.0, 3.0});
    live.frame(x, s.cell, 100 * k);
  }
  EXPECT_EQ(live.frames_written(), 4);
  const Trajectory t = read_lammps_dump(f, &s);
  ASSERT_EQ(t.frames(), 4u);
  EXPECT_EQ(t.timesteps.back(), 300);
  EXPECT_NEAR(t.positions[3][2][0], 3.3, 1e-4);
}

TEST(Live, ThreadsFromTheJob) {
  set_env("CAPS_THREADS", nullptr);
  set_env("SLURM_CPUS_PER_TASK", nullptr);
  EXPECT_EQ(env_threads(), 0);
  set_env("SLURM_CPUS_PER_TASK", "52");
  EXPECT_EQ(env_threads(), 52);
  set_env("CAPS_THREADS", "26");
  EXPECT_EQ(env_threads(), 26);   // CAPS_THREADS first
  set_env("CAPS_THREADS", "x");
  EXPECT_EQ(env_threads(), 52);   // not a number: the job's cores
  set_env("CAPS_THREADS", nullptr);
  set_env("SLURM_CPUS_PER_TASK", nullptr);
  set_env("SLURM_JOB_ID", nullptr);
  set_env("CAPS_JOB", nullptr);
  EXPECT_EQ(default_progress_file(), "");
  set_env("CAPS_JOB", "md-1");
  EXPECT_EQ(default_progress_file(), "progress.jsonl");
  set_env("CAPS_JOB", nullptr);
}
