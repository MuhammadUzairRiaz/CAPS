// CAPS CBMC: configurational-bias Monte Carlo regrowth of chain ends on a built, typed cell.
//
// Siepmann & Frenkel, Mol. Phys. 75, 59 (1992); Frenkel & Smit, Understanding Molecular Simulation (2002), ch. 13.
// A move picks a chain, one of its ends and 1 … max_torsions of the rotatable backbone bonds nearest that end, and
// regrows the end torsion by torsion, the innermost first. At each bond `trials` torsion angles, uniform on the circle,
// are each weighted by exp(−u/kT), u the energy of the atoms that torsion places: van der Waals (shifted to zero at the
// cut-off, 1-4 pairs scaled as the force field says) and damped-shifted-force electrostatics with every atom already
// placed, plus the force field's torsion terms about that bond (class II cross terms and bending–torsion included).
// One is drawn by its weight. The old end is retraced the same way, its own torsion being one of the trials, and the
// new end is accepted with probability min(1, W_new / W_old), W the Rosenbluth weight Π_bonds Σ_trials exp(−u/kT):
// detailed balance holds for uniform torsion trials at fixed bond lengths and angles. Bond lengths, angles, impropers,
// ring and double-bond geometry and every configuration (R/S, E/Z) are untouched; carbonyl C–N and C–O bonds (amides,
// esters), double, aromatic and ring bonds are held.
#pragma once

#include <cstdint>
#include <functional>
#include <string>
#include <vector>

#include "caps/field.hpp"
#include "caps/system.hpp"

namespace caps {

struct CbmcOptions {
  int moves = 1000;             // attempted regrowths
  int trials = 8;               // k, torsion trials per bond
  int max_torsions = 4;         // rotatable bonds regrown per move, at most (1 … this, uniformly)
  int max_atoms = 120;          // an end moved by a bond may hold at most this many atoms
  double temperature = 300;     // K
  double cutoff = 9.0;          // Å, pair energies (at most half the cell's narrowest width)
  bool coulomb = true;          // damped shifted force (Fennell & Gezelter 2006), as Dynamics' default
  double dsf_alpha = 0.2;       // Å⁻¹
  uint64_t seed = 1;
  // about every 1 % of the moves: (moves done, accepted); false stops
  std::function<bool(int, int)> progress;
  // with each progress call: the positions now (whole molecules)
  std::function<void(const std::vector<Vec3>&)> snapshot;
};

struct CbmcReport {
  int attempted = 0, accepted = 0;
  int chains = 0;               // chain ends that can be regrown
  double acceptance = 0;
  double energy_change = 0;     // kcal/mol: the trial-energy model's change over the accepted moves
  double r2_before = 0, r2_after = 0;   // mean squared end-to-end distance of the backbones, Å²
  double cutoff = 0;            // Å, as used
  double seconds = 0;
  std::vector<std::string> notes;
};

// Regrows chain ends of s (positions are made whole per molecule and changed in place; velocities are dropped).
// ff must be the force field assigned to s (per-atom types and charges). Throws std::invalid_argument when the force
// field has pair forms other than Lennard-Jones, hydrogen bonds, virtual sites or explicit pairs, or when no chain end
// can be regrown.
void cbmc_regrow(System& s, const ForceField& ff, const CbmcOptions& o, CbmcReport* report = nullptr);

}  // namespace caps
