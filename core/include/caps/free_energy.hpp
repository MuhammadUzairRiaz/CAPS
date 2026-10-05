// CAPS solvation free energy by thermodynamic integration (C5): the solute's interactions with everything else switched
// off in windows — its electrostatics first (λ_coul 1 → 0, the Lennard-Jones on), then its Lennard-Jones (λ_lj 1 → 0,
// soft-core: Beutler et al. 1994) — sampling ⟨∂U/∂λ⟩ in each window by MD from the previous window's end, and
//   ΔG_solv = ∫₀¹ ⟨∂U/∂λ_coul⟩ dλ_coul + ∫₀¹ ⟨∂U/∂λ_lj⟩ dλ_lj   (trapezoid rule over the windows)
// the free energy of taking the solute from the gas into the solvent (decoupling: its own internal interactions stay
// on in both states, so no gas-phase leg is needed). Pairwise electrostatics (DSF) and no tail correction, as the
// soft-core pairs need.
#pragma once
#include <cstdint>
#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "caps/field.hpp"
#include "caps/system.hpp"

namespace caps {

struct SolvationOptions {
  std::shared_ptr<const ForceField> field;   // null: CAPS's default force field
  EnergyOptions energy;                      // DSF electrostatics; tail is turned off
  std::vector<char> solute;                  // per atom: 1 the solute
  std::vector<double> coul_windows{1.0, 0.75, 0.5, 0.25, 0.0};
  std::vector<double> lj_windows{1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.15, 0.1, 0.05, 0.0};
  double temperature = 300.0, dt = 1.0, tau_t = 100.0;
  double ps = 20.0, equilibrate_ps = 5.0;   // per window: sampled, and run first
  int blocks = 5;
  uint64_t seed = 1;
  std::function<bool(const std::string& stage, double fraction)> progress;
};

struct TiWindow {
  std::string leg;            // "coulomb" or "lj"
  double lambda = 0;
  double dudl = 0, error = 0; // ⟨∂U/∂λ⟩ and its block error, kcal/mol
};

struct SolvationResult {
  double dg = 0, dg_err = 0;              // kcal/mol: ΔG_solv (gas → solvent)
  double dg_coul = 0, dg_lj = 0;          // the two legs
  std::vector<TiWindow> windows;
  std::vector<std::string> notes;
};

SolvationResult solvation_free_energy(const System& s, const SolvationOptions& o);

}  // namespace caps
