// CAPS iterative Boltzmann inversion (Reith, Pütz & Müller-Plathe, J. Comput. Chem. 24, 1624 (2003)): a tabulated
// non-bonded bead–bead potential refined until the coarse-grained model's g(r) matches the target (the mapped all-atom
// g(r) of cg_map, or any other).
//
//   start      V₀(r) = −k_B T ln g_target(r) where g_target > 0, continued below as a straight wall at its slope there,
//              shifted to zero at the cut-off
//   iterate    a CG MD run (NVT, Bussi) from the last configuration; its g(r) over the non-bonded pairs (other chains, or
//              more than three bonds apart along one, as cg_map's); V ← V + α k_B T ln(g_i / g_target) where both are
//              resolved, smoothed over three points, zero at the cut-off; optionally the linear pressure correction
//              ΔV = A (1 − r/r_c), A = −0.1 k_B T sign(ΔP) min(1, 0.0003 |ΔP| / bar)
//   residual   ∫ (g_i − g_target)² dr / ∫ g_target² dr per iteration
// One table for every bead pair (the pooled g(r) of cg_map); the bonded terms stay as given.
#pragma once
#include <cstdint>
#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "caps/field.hpp"
#include "caps/system.hpp"

namespace caps {

struct IbiOptions {
  double temperature = 300;      // K
  int iterations = 8;
  double run_ps = 20, first_equilibrate_ps = 5, dt = 5;   // fs per step
  double cutoff = 0;             // Å, the table's end (0: the target's range, at most 15 Å)
  double dr = 0.05;              // Å, table spacing
  double alpha = 0.3;            // damping of each update
  bool pressure_correction = false;
  double pressure = 1.0;         // atm, the correction's target
  int frame_every = 50;          // steps between g(r) samples
  uint64_t seed = 1;
  std::function<bool(const std::string& stage, double fraction)> progress;
};

struct IbiIteration {
  int iteration = 0;
  double residual = 0;
  double pressure = 0;           // atm, mean over the run
  std::vector<double> gr;        // on the target's r
};

struct IbiResult {
  std::shared_ptr<ForceField> ff;            // the model with the refined table
  std::vector<double> r, target;             // the target g(r)
  std::vector<double> table_r, potential;    // the final V(r), kcal/mol
  std::vector<IbiIteration> history;
  System last;                               // the last configuration
  std::vector<std::string> notes;
};

// The non-bonded pooled g(r) of beads (other molecules, or more than three bonds apart in one) over frames.
std::vector<double> nonbonded_gr(const System& top, const std::vector<std::vector<Vec3>>& frames, const std::vector<Cell>& cells, double dr,
                                 size_t nbin);

// beads: the CG structure (bonds, molecules, a periodic cell); ff: its bead model (bonded terms kept, non-bonded
// replaced by the table); r, g: the target g(r) on an even grid (bin centres).
IbiResult run_ibi(const System& beads, const ForceField& ff, const std::vector<double>& r, const std::vector<double>& g, const IbiOptions& o);

}  // namespace caps
