// CAPS live output: the files a long run writes while it goes, so a run on a cluster can be followed from outside
// (caps job status / tail, the Studio's Jobs page) and a run cut off by a time limit leaves what it had done.
//
//   thermo CSV       one fixed column set for every MD-based command (thermo_header()), a row per report, flushed
//   frames           a LAMMPS dump (id mol type xu yu zu) appended frame by frame, flushed
//   progress.jsonl   one JSON object per report: {"command","stage","stage_index","stages","step","time_ps","fraction",
//                    "T","P","rho","pe","ke","etotal","volume","lx","ly","lz", …extra}, flushed; the last complete line
//                    is the state of the run (a reader drops a line cut off mid-write)
//
// Threads: run_threads() is what a command uses when --threads is not given: CAPS_THREADS, then SLURM_CPUS_PER_TASK
// (the cores a SLURM job was given), else the automatic choice (caps/parallel.hpp).
#pragma once
#include <cstdint>
#include <fstream>
#include <memory>
#include <string>
#include <vector>

#include "caps/dynamics.hpp"
#include "caps/json.hpp"
#include "caps/system.hpp"

namespace caps {

// threads from the environment: CAPS_THREADS, then SLURM_CPUS_PER_TASK; 0 when neither is set (or not a number ≥ 1)
int env_threads();

// "progress.jsonl" when the command runs as a job (SLURM_JOB_ID or CAPS_JOB set), else "" (no file unless asked)
std::string default_progress_file();

// the thermo CSV's columns, the same for md, equilibrate, tg, tensile, pull and recipes
std::string thermo_header();
std::string thermo_csv_row(const ThermoRow& r, const std::string& stage = "");

class LiveOutput {
 public:
  explicit LiveOutput(std::string command = "") : command_(std::move(command)) {}
  void open_thermo(const std::string& path);
  void open_frames(const std::string& path, const System& topology);
  void open_progress(const std::string& path);
  bool thermo_open() const { return thermo_ != nullptr; }
  bool frames_open() const { return frames_ != nullptr; }
  bool progress_open() const { return progress_ != nullptr; }

  void thermo(const ThermoRow& r, const std::string& stage = "");
  void frame(const std::vector<double>& x, const Cell& c, int64_t step);
  // one progress line; stage_index from 1 (0: a single-stage command); fraction of the whole command (0 … 1, < 0: unknown);
  // the thermo values when r is given; extra keys merged in (stress for tensile, temperature step for tg …)
  void progress(const std::string& stage, int stage_index, int stages, double fraction, const ThermoRow* r = nullptr,
                const Json& extra = Json());
  // the last line: {"state": "finished"|"failed", …}
  void done(bool ok, const std::string& note = "");
  int64_t frames_written() const { return nframes_; }

 private:
  std::string command_;
  std::unique_ptr<std::ofstream> thermo_, frames_, progress_;
  System topology_;
  int64_t nframes_ = 0;
};

// One frame in LAMMPS dump format (the frame writer of write_lammps_dump)
void write_lammps_dump_frame(std::ostream& out, const System& topology, const std::vector<Vec3>& pos, const Cell& c, int64_t step);

// The complete lines of a chunk of a progress file, parsed in order; a last line cut off mid-write is left for the next
// read: consumed is how many bytes were used (the next read starts there).
std::vector<Json> read_progress_lines(const std::string& text, size_t* consumed = nullptr);

}  // namespace caps
