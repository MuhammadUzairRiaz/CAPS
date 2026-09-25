// CAPS trajectory player data (design/boards/Trajectory): per-frame series for plots linked to the frame, a LAMMPS
// log's thermo columns matched to the frames by timestep, and the end-to-end vector of one chain.
#pragma once
#include <functional>
#include <string>
#include <vector>

#include "caps/pipeline.hpp"
#include "caps/system.hpp"

namespace caps {

// Thermo output of a LAMMPS log: the columns (Step first) and one row per printed step, all runs in order. Both the
// one-line style ("Step Temp Press …" tables up to "Loop time") and the multi style ("---- Step N ----" blocks of
// "Name = value" pairs) are read. run_starts: the row where each run begins.
struct ThermoLog {
  std::vector<std::string> columns;
  std::vector<std::vector<double>> rows;
  std::vector<size_t> run_starts;
};
ThermoLog parse_lammps_log(const std::string& text);
ThermoLog read_lammps_log(const std::string& path);

struct TrajectorySeriesOptions {
  int molecule = 0;            // the chain for Rg and Ree (molecule id; 0: the largest)
  double dt_fs = 1.0;          // fs per timestep, for the time axis
  const ThermoLog* log = nullptr;
  int stride = 1;
  std::function<bool(int, int)> progress;   // (done, total) → false stops
};

// One row per frame: Frame, Timestep, Time (ps), Density (g/cm³) and Volume (Å³) with a cell, Rg and Ree (Å) of the
// chain, then the log's columns at the frame's timestep (the nearest printed step; NaN when none is within half the
// print interval).
DataTable trajectory_series(const Trajectory& t, const TrajectorySeriesOptions& o);

// The molecule id used for the chain (the largest when molecule is 0) and its backbone's first and last atoms.
struct ChainEnds { int molecule = 0; int first = -1, last = -1; };
ChainEnds chain_ends(const System& s, int molecule);

// Positions averaged over `window` frames centred on frame k (each atom unwrapped to frame k's image); window ≤ 1
// returns frame k's positions.
std::vector<Vec3> smoothed_positions(const Trajectory& t, size_t k, int window);

}  // namespace caps
