// CAPS sorption: how much of a gas or small molecule a fixed host (a polymer or elastomer cell, a filler, a porous
// crystal) takes up, by test-particle insertion (Widom, J. Chem. Phys. 39, 2808 (1963)) and grand-canonical Monte Carlo
// (Frenkel & Smit, Understanding Molecular Simulation, ch. 5).
//
// The host stays fixed; the sorbate is a molecule (the template: one copy typed together with the host, its atoms the
// last of the structure), rigid unless `flexible` (below). Energies are the force field's non-bonded ones (pairmodel.hpp).
//   Widom: ⟨W⟩ = ⟨e^(−βΔU)⟩ over random places and orientations of one sorbate in the empty host; μ_ex = −kT ln⟨W⟩; at low
//   pressure the uptake is p V ⟨W⟩ / kT molecules, so the Henry constant K_H = ⟨W⟩ V / (R T m_host) and the solubility
//   coefficient S = ⟨W⟩ T₀ / T  in cm³(STP) per cm³ of host per atm (T₀ = 273.15 K).
//   GCMC at each pressure (ideal-gas reservoir, fugacity = pressure): insertions and deletions accepted with
//   min(1, βfV/(N+1) e^(−βΔU)) and min(1, N/(βfV) e^(−βΔU)), translations and rotations with min(1, e^(−βΔU)); the
//   loading is the average N over the production steps, its error from ten blocks, and the isosteric heat from the
//   fluctuations Q_st = kT − (⟨UN⟩ − ⟨U⟩⟨N⟩) / (⟨N²⟩ − ⟨N⟩²).
//   Mixtures (C6): several sorbate templates one after the other, gas mole fractions y_i, each species' fugacity y_i p
//   (ideal reservoir); an insertion or deletion picks a species at random and is accepted with that species' f_i and N_i;
//   per species loadings, Widom terms and isosteric heats q_i = kT − Σ_j cov(U, N_j) [cov(N, N)⁻¹]_ji (the
//   multicomponent fluctuation formula); the adsorption selectivity S_i/0 = (x_i / x_0) / (y_i / y_0).
//
// Flexible sorbates (flexible = true): a sorbate with rotatable bonds is cut there into rigid pieces (bond lengths,
// angles, rings and double bonds keep the template's geometry) and grown piece by piece by configurational bias
// (Siepmann & Frenkel, Mol. Phys. 75, 59 (1992); Smit & Siepmann, J. Phys. Chem. 98, 8442 (1994)): the first piece at
// k random places and orientations, each later one at k torsion angles uniform on the circle, each trial weighted by
// e^(−βu) (u: the new atoms with the host, the other sorbates and the atoms already placed, 1-4 pairs scaled, plus the
// torsion terms about that bond), one drawn by its weight. W = Π_steps (1/k) Σ_trials e^(−βu) is the Rosenbluth weight.
// The same growth with no host gives the isolated molecule's ⟨W_IG⟩ (the ideal gas reference), so
//   Widom: μ_ex = −kT ln(⟨W⟩ / ⟨W_IG⟩), and the Henry constant and solubility from ⟨W⟩ / ⟨W_IG⟩ as above;
//   GCMC: insertion accepted with min(1, βfV/(N+1) W / ⟨W_IG⟩), deletion (the molecule retraced, its own conformation
//   one of the trials) with min(1, N/(βfV) ⟨W_IG⟩ / W); partial regrowth beyond a random rotatable bond with
//   min(1, W_new / W_old); twists of one rotatable bond and rigid translations and rotations by Metropolis;
//   isosteric heat q = kT + ⟨u_intra⟩_gas − cov(U, N) / var(N), U holding the intramolecular energy too and
//   ⟨u_intra⟩_gas = ⟨u W⟩ / ⟨W⟩ over the isolated growths.
// Fixed loading (fixed_loading, molecules per species): canonical Monte Carlo at that loading (the molecules grown in
// by configurational bias first, then translations, rotations, regrowths and twists): the mean energy and its parts, the
// adsorption energy per molecule (⟨U⟩ − Σ N_i ⟨u_intra,i⟩_gas) / N, the distribution of single molecules' binding
// energies (with the host and the other sorbates), the lowest-energy configuration met and a density map.
#pragma once

#include <cstdint>
#include <functional>
#include <string>
#include <vector>

#include "caps/field.hpp"
#include "caps/system.hpp"

namespace caps {

struct SorptionOptions {
  int template_first_atom = -1;      // the sorbate template's atoms: this one to the last (required)
  double temperature = 300;          // K
  int insertions = 100000;           // Widom (0: none)
  std::vector<double> pressures_kpa; // GCMC isotherm points (empty: Widom only)
  int steps = 200000;                // GCMC steps per pressure (the first quarter equilibrates)
  double cutoff = 12.0;              // Å
  bool coulomb = true;
  double dsf_alpha = 0.2;
  uint64_t seed = 1;
  // a density map of the sorbate in the host (C7): its molecules' centres sampled over the production steps of each
  // pressure on a grid of map_grid points along each cell edge (0: no map)
  int map_grid = 0;
  // mixtures: the first atom of each species' template (ascending; each runs to the next or to the end; empty: one
  // species from template_first_atom) and their gas mole fractions (normalised; empty: equal)
  std::vector<int> species_first_atom;
  std::vector<double> mole_fractions;
  // flexible sorbates: rotatable bonds grown by configurational bias with `trials` trials per step; the isolated
  // molecule's Rosenbluth weight from ig_samples growths; twists up to twist_max_deg
  bool flexible = false;
  int trials = 8;
  int ig_samples = 20000;
  double twist_max_deg = 45;
  // fixed loading: molecules per species (empty: none); grown in with at most loading_attempts tries, then `steps`
  // canonical Monte Carlo steps (the first quarter equilibrates)
  std::vector<int> fixed_loading;
  int loading_attempts = 100000;
  int energy_bins = 60;
  std::function<bool(const std::string& stage, double fraction)> progress;   // false stops
};

struct IsothermPoint {
  double pressure_kpa = 0;
  double loading = 0, loading_error = 0;   // molecules per cell
  double mol_per_kg = 0;                   // per kg of host
  double cm3stp_per_cm3 = 0;               // per cm³ of cell
  double heat = 0;                         // isosteric heat, kcal/mol (positive: exothermic)
  double acceptance_insert = 0, acceptance_delete = 0;
  double acceptance_regrow = 0, acceptance_twist = 0;   // flexible sorbates
  std::vector<std::vector<Vec3>> molecules;   // the sorbate molecules at the end of the run
  // the density map (map_grid > 0): sorbate centres per Å³, averaged over the samples, index (i·g + j)·g + k along a, b, c
  std::vector<float> density;
  int grid = 0;
  // per species (one entry for a pure gas): molecules per cell, errors, mol/kg, isosteric heats, selectivity over species 0
  std::vector<double> species_loading, species_error, species_mol_per_kg, species_heat, selectivity;
};

struct SpeciesWidom {
  double fraction = 1;                     // gas mole fraction
  double widom_w = 0, widom_error = 0, mu_ex = 0, henry_mol_kg_kpa = 0, solubility = 0;   // widom_w: ⟨W⟩ / ⟨W_IG⟩ when flexible
  int rotatable = 0;                       // rotatable bonds grown (0: rigid)
  double w_ig = 1, w_ig_error = 0;         // the isolated molecule's ⟨W_IG⟩
  double u_intra_gas = 0;                  // ⟨u_intra⟩ of the isolated molecule, kcal/mol
};

// The canonical run at a fixed loading
struct FixedLoadingResult {
  std::vector<int> molecules;              // per species, as loaded
  int requested = 0, loaded = 0;
  double energy = 0, energy_error = 0;     // ⟨U⟩ per cell, kcal/mol
  double host_sorbate = 0, sorbate_sorbate = 0, intramolecular = 0;   // its parts (intramolecular: torsions and pairs across rotatable bonds)
  double adsorption_energy = 0;            // per molecule, kcal/mol (negative: bound)
  std::vector<double> hist_edges, hist_counts;   // single molecules' binding energies (kcal/mol), counts
  double lowest_energy = 0;                // the lowest U met in production
  std::vector<std::vector<Vec3>> lowest;   // its molecules
  std::vector<std::vector<Vec3>> molecules_end;
  std::vector<float> density;
  int grid = 0;
  double acceptance_translate = 0, acceptance_rotate = 0, acceptance_regrow = 0, acceptance_twist = 0;
};

struct SorptionReport {
  double widom_w = 0, widom_error = 0;     // ⟨e^(−βΔU)⟩ and its standard error (blocks)
  double mu_ex = 0;                        // kcal/mol
  double henry_mol_kg_kpa = 0;             // mol per kg of host per kPa
  double solubility = 0;                   // cm³(STP) / (cm³ atm)
  double host_mass = 0;                    // g/mol (of the cell's host)
  double volume = 0;                       // Å³
  std::vector<IsothermPoint> isotherm;
  std::vector<SpeciesWidom> species;       // per species (the top-level Widom values are species 0's)
  bool fixed_ran = false;
  FixedLoadingResult fixed;
  double seconds = 0;
  std::vector<std::string> notes;
};

// s: host followed by one sorbate template (typed together in ff); s is not changed.
SorptionReport sorption(const System& s, const ForceField& ff, const SorptionOptions& o);

}  // namespace caps
