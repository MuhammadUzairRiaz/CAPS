// CAPS Equilibrate: published equilibration protocols run as chained Dynamics stages, with convergence checks.
#pragma once
#include <functional>
#include <string>
#include <vector>

#include "caps/dynamics.hpp"
#include "caps/system.hpp"

namespace caps {

enum class Ensemble { NVE, NVT, NPT };

struct Stage {
  std::string label;            // shown in reports, e.g. "3 · compress 0.02 Pmax"
  Ensemble ensemble = Ensemble::NVT;
  double ps = 10;               // duration
  double t_start = 300, t_end = -1;   // K; t_end ≥ 0 ramps the thermostat target
  double pressure = 1;          // atm, NPT only
  double force_cap = 0;         // kcal/mol/Å; > 0 caps LJ forces (push-off)
};

struct ProtocolParams {
  double t_final = 300;         // K
  double t_max = 600;           // K; choose about Tg + 200–300 K for the polymer
  double p_final = 1;           // atm
  double p_max = 49346.2;       // atm (5 × 10⁴ bar, the Larsen et al. value)
  double time_scale = 1.0;      // multiplies every stage duration (< 1 gives a shortened protocol, labelled as such)
  // Annealing
  int cycles = 3;
  double t_low = 300, t_high = 600;
  double ramp_ps = 50, hold_ps = 50;
  // Push-off
  std::vector<double> caps{5, 10, 20, 50, 100, 200, 500};
  double cap_ps = 5;
};

// The 21-step compression and decompression protocol (Larsen, Lin, Hart and Colina, Macromolecules 44, 6944 (2011)):
// alternating NVT at t_max and t_final with NPT compressions to 0.02, 0.6, 1, 0.5, 0.1 and 0.01 Pmax, then a long
// NPT stage at t_final and p_final. 1.56 ns at time_scale 1.
std::vector<Stage> larsen21(const ProtocolParams& p);
// Simulated annealing: NPT ramps between t_low and t_high with holds, `cycles` times, ending at t_low.
std::vector<Stage> annealing(const ProtocolParams& p);
// Molecular-dynamics push-off: NVT at t_final with LJ forces capped at increasing values, then uncapped NVT and NPT
// (after Auhl, Everaers, Grest, Kremer and Plimpton, J. Chem. Phys. 119, 12718 (2003)).
std::vector<Stage> pushoff(const ProtocolParams& p);
std::vector<Stage> protocol_by_name(const std::string& name, const ProtocolParams& p);   // larsen21 | annealing | pushoff

// Plain-text protocols, one stage per line, '#' starts a label or comment:
//   nvt 50 ps T 600                 npt 50 ps T 300 P 987 atm (or bar)
//   nvt 20 ps T 300 to 600          nvt 5 ps T 300 cap 20        nve 10 ps
std::vector<Stage> parse_protocol(const std::string& text);
std::string protocol_text(const std::vector<Stage>& stages);
double protocol_ps(const std::vector<Stage>& stages);

struct ConvergenceCheck {
  std::string quantity;         // "density", "potential energy per atom", "mean Rg"
  std::vector<double> blocks;   // block means
  double change = 0;            // last relative (or per-atom) change between blocks
  double tolerance = 0;
  bool ok = false;
};

struct EquilibrateOptions {
  std::vector<Stage> stages;
  DynamicsOptions md;           // dt, thermostat / barostat kinds and times, seed, energy options; per-stage fields are set
  double frame_ps = 10;         // record a frame this often (and at the end of each stage)
  double thermo_ps = 0.5;
  // After the stages: production NPT at the last stage's conditions in blocks, until the checks pass.
  bool until_converged = false;
  double block_ps = 20;
  int min_blocks = 3, max_blocks = 20;
  double tol_density = 0.005;   // relative change of block means
  double tol_energy = 0.005;    // kcal/mol per atom
  double tol_rg = 0.02;         // relative
  std::function<bool(int stage, int stages, const std::string& label, const ThermoRow& row)> progress;
  std::function<void(const std::vector<double>& x, const Cell& c, int64_t step)> frame;
};

struct StageSummary {
  std::string label;
  double ps = 0;
  double temperature = 0, pressure = 0, density = 0, potential = 0;   // means over the second half of the stage
  double density_end = 0;
};

struct EquilibrateReport {
  std::vector<StageSummary> stages;
  std::vector<ThermoRow> thermo;
  std::vector<ConvergenceCheck> checks;
  bool converged = false;
  int blocks = 0;
  double ps = 0, seconds = 0;
  std::vector<std::string> notes;
};

void equilibrate(System& s, const EquilibrateOptions& o, EquilibrateReport* report = nullptr);

}  // namespace caps
