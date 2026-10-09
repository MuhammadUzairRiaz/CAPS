// CAPS Bench: the validation suite that runs with CAPS alone (design board Bench). Each table measures one claim on
// the shipped samples and says pass, fail or info; tables that need other programs or long equilibrated runs are
// listed as not run with the reason, so the suite never reports what it did not measure.
//
//   T1  forces and virial against finite differences of the energy
//   T2  Pack on water boxes: wall time, closest contact, success
//   T4  molecular-dynamics throughput (ns/day)
//   T5  NVE energy conservation (drift, kT/ns/atom)
//   T6  properties of equilibrated cells (the user's) against experiment: density
//   T7  chain statistics of the same cells against the literature: C∞
//   T8  thread scaling of the force evaluation
//   T9  reproducibility: the same run twice is bit-identical
//   T11 rendering time
//   T12 molecule builder: stereochemistry and round trips
//   T13 coarse-grained mapping of grown copolyesters; T14 bonded Boltzmann inversion of an ideal chain; T15 IBI
//   self-consistency on a Lennard-Jones fluid; T16, T17 Kremer–Grest entanglement and hardening (not run: cluster runs);
//   T18 API regressions
#pragma once
#include <functional>
#include <string>
#include <vector>

namespace caps {

struct BenchRow {
  std::vector<std::string> cells;
  std::string status;           // pass | fail | info
  std::string file;             // the structure the row was measured on (Open run), when there is one
};

struct BenchTable {
  std::string id, title, scope;
  std::vector<std::string> columns;
  std::vector<BenchRow> rows;
  std::string status = "not run";   // pass | fail | info | not run
  std::string note;             // how it was measured, or why it was not run
  double seconds = 0;
};

struct BenchOptions {
  std::string samples;          // the samples directory (ps_melt.data, water.pdb)
  std::string forcefields;      // data/forcefields (T12 cleans molecules up with GAFF2 when given)
  int repeats = 3;
  bool quick = false;           // shorter runs (tests, a first look)
  // T6 / T7: equilibrated cells to measure, named by their reference id (ps-atactic.data, with ps-atactic.lammpstrj or
  // .dcd beside it for frames), and the reference ranges (data/reference/polymers.json)
  std::string cells;
  std::string reference;
  // (table id, what, fraction of the table done) → false cancels
  std::function<bool(const std::string&, const std::string&, double)> progress;
};

std::vector<std::string> bench_ids();                  // T1 … T12
BenchTable bench_describe(const std::string& id);      // title, scope and columns, not run
BenchTable run_bench(const std::string& id, const BenchOptions& o);

std::string bench_markdown(const std::vector<BenchTable>& tables);
std::string bench_csv(const BenchTable& t);
std::string bench_latex(const std::vector<BenchTable>& tables);

}  // namespace caps
