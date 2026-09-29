#include <gtest/gtest.h>

#include <cmath>

#include <filesystem>

#include "caps/io.hpp"
#include "caps/peptide.hpp"
#include "caps/trajectory.hpp"

using namespace caps;

TEST(Trajectory, ReadsBothLogStyles) {
  const std::string one =
      "LAMMPS (2 Aug 2023)\n"
      "Step Temp Press Density\n"
      "       0   300.0   1.0   1.05\n"
      "     100   301.5   -3.0  1.04\n"
      "Loop time of 1.0 on 1 procs for 100 steps\n"
      "run 100\n"
      "   Step          Temp          Press    \n"
      "     100   302.0   2.0\n"
      "     200   299.0   5.0\n"
      "Loop time of 1.0 on 1 procs for 100 steps\n";
  const ThermoLog a = parse_lammps_log(one);
  ASSERT_EQ(a.columns.size(), 4u);
  EXPECT_EQ(a.columns[1], "Temp");
  EXPECT_EQ(a.rows.size(), 4u);
  EXPECT_EQ(a.run_starts, (std::vector<size_t>{0, 2}));
  EXPECT_DOUBLE_EQ(a.rows[3][1], 299.0);
  EXPECT_TRUE(std::isnan(a.rows[3][3]));   // the second run printed no density
  const std::string multi =
      "------------ Step              0 ----- CPU =            0 (sec) -------------\n"
      "TotEng   =    -25356.2057 KinEng   =     21444.8313 Temp     =       299.0397 \n"
      "Volume   =    307995.0335\n"
      "------------ Step             50 ----- CPU =     1.894152 (sec) -------------\n"
      "TotEng   =    -25330.0307 KinEng   =     21501.0009 Temp     =       299.8229 \n"
      "Volume   =    308031.6762\n"
      "Loop time of 3.8862 on 4 procs for 100 steps with 32000 atoms\n";
  const ThermoLog b = parse_lammps_log(multi);
  EXPECT_EQ(b.rows.size(), 2u);
  EXPECT_EQ(b.columns.size(), 5u);   // Step TotEng KinEng Temp Volume
  EXPECT_DOUBLE_EQ(b.rows[1][0], 50);
  EXPECT_NEAR(b.rows[1][3], 299.8229, 1e-9);
  EXPECT_THROW(parse_lammps_log("nothing here\n"), std::runtime_error);
}

TEST(Trajectory, SeriesFollowTheChainAndTheLog) {
  PeptideOptions po;
  po.sequence = "GAGAGAGA";
  po.structure = "PPPPPPPP";
  po.cleanup = false;
  System pep = build_peptide(po);
  pep.cell.a = {60, 0, 0}, pep.cell.b = {0, 60, 0}, pep.cell.c = {0, 0, 60};
  Trajectory t;
  t.topology = pep;
  for (int k = 0; k < 3; ++k) {
    std::vector<Vec3> p;
    for (const auto& a : pep.atoms) p.push_back(a.pos * (1.0 + 0.1 * k));   // the chain stretches
    t.positions.push_back(p);
    t.cells.push_back(pep.cell);
    t.timesteps.push_back(1000 * k);
  }
  ThermoLog log = parse_lammps_log("Step Temp\n 0 300\n 1000 310\n 2000 320\nLoop time of 1\n");
  TrajectorySeriesOptions o;
  o.log = &log;
  o.dt_fs = 2.0;
  const DataTable T = trajectory_series(t, o);
  ASSERT_EQ(T.rows.size(), 3u);
  auto col = [&](const std::string& n) { return size_t(std::find(T.columns.begin(), T.columns.end(), n) - T.columns.begin()); };
  EXPECT_NEAR(T.rows[2][col("Time (ps)")], 4.0, 1e-12);
  EXPECT_NEAR(T.rows[1][col("Ree (Å)")] / T.rows[0][col("Ree (Å)")], 1.1, 1e-9);
  EXPECT_NEAR(T.rows[2][col("Rg (Å)")] / T.rows[0][col("Rg (Å)")], 1.2, 1e-9);
  EXPECT_DOUBLE_EQ(T.rows[1][col("Temp")], 310);
  const ChainEnds e = chain_ends(pep, 0);
  EXPECT_GE(e.first, 0);
  EXPECT_NE(e.first, e.last);
  // smoothing over three frames gives the middle frame
  const auto s = smoothed_positions(t, 1, 3);
  EXPECT_NEAR(norm(s[5] - t.positions[1][5]), 0.0, 1e-9);
}

// Reading a frame selection: 20 frames written as a dump, a DCD and an XYZ; frames 2 … 15 every 4 read back as
// 2, 6, 10, 14 with their own positions and timesteps (the dump passes over the others unread and stops after 15)
TEST(Trajectory, ReadsAFrameSelection) {
  const std::string data = std::string(CAPS_SAMPLES) + "/ps_melt.data";
  const Trajectory src = open_file(data);
  Trajectory t;
  t.topology = src.topology;
  for (int k = 0; k < 20; ++k) {
    auto p = src.positions.front();
    p[0][0] += 0.05 * k;
    t.positions.push_back(p);
    t.cells.push_back(src.topology.cell);
    t.timesteps.push_back(100 * k);
  }
  const auto dir = std::filesystem::temp_directory_path() / "caps_frame_selection";
  std::filesystem::create_directories(dir);
  OpenProgress sel;
  sel.frames = {2, 15, 4};
  for (const char* ext : {".lammpstrj", ".dcd", ".xyz"}) {
    const std::string path = (dir / (std::string("t") + ext)).string();
    write_trajectory(t, path, 1.0);
    const Trajectory r = open_file(path, std::string(ext) == ".xyz" ? "" : data, sel);
    ASSERT_EQ(r.frames(), 4u) << ext;
    for (size_t j = 0; j < 4; ++j) {
      const size_t k = 2 + 4 * j;
      EXPECT_NEAR(r.positions[j][0][0], t.positions[k][0][0], 2e-3) << ext << " frame " << k;
      EXPECT_NEAR(r.positions[j][5][1], t.positions[k][5][1], 2e-3) << ext;
    }
    EXPECT_NEAR(r.topology.atoms[0].pos[0], t.positions[2][0][0], 2e-3) << ext << ": the topology follows the first kept frame";
    if (std::string(ext) == ".lammpstrj") EXPECT_EQ(r.timesteps, (std::vector<int64_t>{200, 600, 1000, 1400}));
    if (std::string(ext) == ".xyz") EXPECT_EQ(r.timesteps, (std::vector<int64_t>{2, 6, 10, 14})) << "XYZ frames are numbered as in the file";
    if (std::string(ext) == ".lammpstrj") EXPECT_EQ(r.frames_read, 16u) << "reading stops after the last frame asked for";
    OpenProgress past;
    past.frames = {30, SIZE_MAX, 1};
    EXPECT_THROW(open_file(path, std::string(ext) == ".xyz" ? "" : data, past), std::exception) << ext;
  }
  // the same frames from those held in memory
  Trajectory m = t;
  m.select({2, 15, 4});
  ASSERT_EQ(m.frames(), 4u);
  EXPECT_EQ(m.timesteps, (std::vector<int64_t>{200, 600, 1000, 1400}));
  std::filesystem::remove_all(dir);
}
