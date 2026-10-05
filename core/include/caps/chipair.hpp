// CAPS χ from pair contacts (design/boards/SolventScreen, BlendPhase): the Flory–Huggins interaction parameter as a
// function of temperature from the energies of molecule pairs in contact and coordination numbers, after Fan, Olafson,
// Blanco & Hsu, Macromolecules 25, 3667 (1992) and Blanco, J. Comput. Chem. 12, 237 (1991).
//
// Each component is one rigid molecule: a solvent, or a polymer's repeat unit with its two ends capped by hydrogen
// (*CC(*)c1ccccc1 → ethylbenzene). For each pair of kinds (AA, AB, BB) CAPS draws random relative orientations and
// directions, slides the second molecule along the direction until the two van der Waals surfaces touch (Bondi radii
// × contact_scale; the contact distance is solved exactly from every atom pair), and records the non-bonded energy of the
// pair (Lennard-Jones and Coulomb, ε = 1, the same force field for both). The energies are Boltzmann-averaged at each
// temperature:
//
//     ⟨E_ij⟩_T = Σ E e^(−E/RT) / Σ e^(−E/RT)
//
// The coordination number Z_ij is how many j molecules fit around one i molecule at contact without overlapping each
// other: j molecules are added at random contacts, each kept when it overlaps none already placed, until `pack_misses`
// placements in a row fail; the mean over `pack_trials` shells. Then
//
//     ΔE_mix(T) = ½ [Z_AB ⟨E_AB⟩ + Z_BA ⟨E_BA⟩ − Z_AA ⟨E_AA⟩ − Z_BB ⟨E_BB⟩],      χ(T) = ΔE_mix / RT
//
// and χ(T) is fitted to A + B/T. The lattice site is one molecule of each kind (no volume correction between kinds of
// different size), and the molecules keep one conformer: a screen for trends and rankings, as the paper says, not a
// replacement for measured χ. Standard errors come from ten blocks of the samples.
//
// Force field: one for both molecules — a library force field with typing rules (the Studio uses GAFF2), else GAFF's
// built-in C and H subset when it types both, else UFF; charges Gasteiger–Marsili (QEq with UFF where Gasteiger has none).
#pragma once
#include <cstdint>
#include <stdexcept>
#include <functional>
#include <string>
#include <vector>

#include "caps/system.hpp"

namespace caps {

struct ChiPairOptions {
  std::string a_smiles, b_smiles;   // molecules or repeat units (* ends capped with H)
  // A caps-forcefield JSON with typing rules (GAFF2, OPLS-AA …) for both molecules, Gasteiger charges; empty: GAFF's
  // built-in C and H subset when it types both, else UFF with Gasteiger / QEq charges. Only the Lennard-Jones types and
  // charges are used, so missing bonded parameters do not matter.
  std::string forcefield;
  int samples = 1000000;            // contact configurations per pair of kinds
  int pack_trials = 5000;           // shells packed per coordination number
  int pack_misses = 1000;           // failed placements in a row that end a shell
  double contact_scale = 1.0;       // contact at scale × (R_i + R_j), Bondi radii
  std::vector<double> temperatures{250, 275, 300, 325, 350, 375, 400};   // K
  double report_temperature = 298.15;
  uint64_t seed = 1;
  int threads = 0;
  std::function<bool(const std::string& stage, double fraction)> progress;   // return false to cancel
};

struct ChiPairKind {
  std::string name;                 // "A–A", "A–B", "B–A", "B–B"
  double z = 0, z_error = 0;        // coordination number
  double e_min = 0, e_mean = 0;     // kcal/mol: the lowest sampled and the plain mean
  std::vector<double> e_t, e_t_error;   // ⟨E⟩_T per temperature, kcal/mol
  std::vector<double> hist_e, hist_p;   // energy histogram (bin centres, probability)
};

struct ChiPairResult {
  std::string a_name, b_name;       // the molecules as built (formula)
  int a_atoms = 0, b_atoms = 0;
  std::string forcefield;
  ChiPairKind aa, ab, ba, bb;
  std::vector<double> temperatures, chi, chi_error;
  double fit_a = 0, fit_b = 0;      // χ(T) = A + B/T
  double chi_report = 0, chi_report_error = 0;   // at report_temperature
  std::vector<std::string> notes;
};

struct ChiPairCancelled : std::runtime_error {
  ChiPairCancelled() : std::runtime_error("χ run cancelled") {}
};

// The SMILES a repeat unit becomes: every * replaced by [H].
std::string capped_unit(const std::string& unit_smiles);

ChiPairResult chi_by_contacts(const ChiPairOptions& o);

}  // namespace caps
