// CAPS coarse-grained melts (design/boards/CoarseGrained): Kremer–Grest bead-spring chains (J. Chem. Phys. 92, 5057
// (1990)) — M chains of N beads at a reduced density ρσ³ in a cubic box L = (M N / ρ)^{1/3} σ — built as random walks
// with a bond of 0.97 σ that never folds straight back (bond angles above 60°; a bending stiffness k_θ makes straighter
// walks more likely), and written for LAMMPS in reduced units with the standard preparation: a soft-core push-off whose
// strength is ramped up (Auhl et al. 2003), then FENE bonds (K = 30 ε/σ², R₀ = 1.5 σ) and the WCA repulsion (LJ cut at
// 2^{1/6} σ and shifted), Langevin thermostat at T = 1.
//
// In the structure, σ = 1 Å numerically (positions and box in σ); beads carry element 6 for drawing and mass 1.
#pragma once
#include <cstdint>
#include <string>

#include "caps/system.hpp"

namespace caps {

struct KgOptions {
  int chains = 50;            // M
  int beads = 100;            // N
  double density = 0.85;      // ρσ³
  double k_theta = 0.0;       // bending stiffness, ε (0: flexible)
  double bond = 0.97;         // σ, the FENE–WCA bond length at T = 1
  uint64_t seed = 1;
};

struct KgReport {
  double box = 0;             // σ
  double closest = 0;         // closest non-bonded pair, σ (overlaps are expected before the push-off)
  double mean_r2 = 0;         // ⟨R_ee²⟩ / (N − 1) of the built chains, σ² per bond
};

System kremer_grest(const KgOptions& o, KgReport* rep = nullptr);

// LAMMPS data (atom_style angle; bond type 1; angle type 1 when k_θ > 0) and an input deck: push-off with pair soft, then
// the Kremer–Grest force field; `stem` gets .data and .in appended.
void write_kg_lammps(const System& s, const KgOptions& o, const std::string& stem, double pushoff_steps = 20000, double run_steps = 100000);

}  // namespace caps
