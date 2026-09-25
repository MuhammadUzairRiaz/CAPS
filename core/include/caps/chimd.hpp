// CAPS χ from molecular dynamics (design/boards/SolventScreen, BlendPhase): the Flory–Huggins interaction parameter from
// the energy of mixing. Three cells are built and run with the same force field (the built-in typing: GAFF for C and H,
// UFF otherwise) — component A alone, component B alone and a mixture — each relaxed, then NPT at the temperature; the
// cohesive energy density (Σ E of each molecule alone − E of the cell, per volume, as Analyze's CED) is taken over ten
// frames of the production run. With volume fractions φ from the pure components' volumes,
//
//     Δe_mix = φ_A CED_A + φ_B CED_B − CED_mix,      χ = V_ref Δe_mix / (R T φ_A φ_B)
//
// (Flory–Huggins; V_ref the solvent's molecular volume, or for two polymers the geometric mean of the repeat units').
// This is the enthalpic χ only: no entropic part, no contact statistics; short runs make it noisy (the standard error
// comes from block averages of the three energies). A test of trends, not a replacement for measured χ.
//
// Measured resolution (natural rubber mixed with itself, where χ must be 0; 6 chains of 10 units per cell): χ = −5.0
// after 4 ps of NPT per cell, −2.9 after 50 ps, −1.2 ± 0.5 after 300 ps (31 CPU-minutes). Real χ of interest are
// 0.01–1, so the Studio does not offer this; use it from Python or C with long runs and larger cells, and run the
// self-mixing control beside it.
#pragma once
#include <functional>
#include <string>
#include <vector>

#include "caps/polymer.hpp"

namespace caps {

struct ChiMdOptions {
  ChainSpec polymer;                 // component A
  int chains = 6;                    // chains of A in its pure cell (the mixture takes half)
  // component B: a solvent from SMILES, or a second polymer when b_polymer is set
  std::string solvent_smiles;
  int solvent_molecules = 0;         // in B's pure cell; 0: the same mass as A's chains
  bool b_polymer = false;
  ChainSpec polymer_b;
  int chains_b = 6;
  double temperature = 300;          // K
  double pressure = 1.0;             // atm
  double eq_ps = 10, prod_ps = 10;   // NPT equilibration and production per cell
  double dt = 1.0;                   // fs
  double cutoff = 10.0;              // Å
  uint64_t seed = 1;
  int threads = 0;
  std::function<bool(const std::string& stage, double fraction)> progress;   // return false to cancel
};

struct ChiMdCell {
  std::string name;
  int atoms = 0, molecules = 0;
  double energy = 0, volume = 0;     // ⟨E_pot⟩ kcal/mol, ⟨V⟩ Å³
  double e_density = 0, e_error = 0; // cohesive energy density, kcal/mol/Å³, and its standard error
  double density = 0;                // g/cm³
};

struct ChiMdResult {
  ChiMdCell a, b, mix;
  double phi_a = 0;                  // volume fraction of A in the mixture
  double v_ref = 0;                  // Å³
  double de_mix = 0;                 // kcal/mol/Å³
  double chi = 0, chi_error = 0;
  std::vector<std::string> notes;
};

struct ChiMdCancelled : std::runtime_error {
  ChiMdCancelled() : std::runtime_error("χ run cancelled") {}
};

ChiMdResult chi_by_md(const ChiMdOptions& o);

}  // namespace caps
