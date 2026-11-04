// CAPS force-field fitting to crystals: least-squares refinement of force-field parameters against observed crystal
// structures and properties (Gale, Philos. Mag. B 73, 3 (1996)).
//
//  Parameters   any pair parameter (Buckingham A, ρ, C; Morse D, α, r0; Lennard-Jones ε, σ of a pair or of a type), bond
//               parameter (harmonic, class II, Morse: core–shell springs included) or type charge, by type names; a
//               charge can carry a balancing type that keeps every structure's total charge.
//  Structures   each crystal is evaluated in its periodic copies (a supercell at least twice the cut-off wide, the copies
//               moving together), so a cell of any size sees every neighbour within the cut-off. Shells (positions not
//               observed) are relaxed in every evaluation.
//  Structure    conventional fitting: the observed structure should be a stationary point; its residual forces
//  targets      (kcal/mol/Å, the shells relaxed) and stresses (GPa) are fitted to zero. Relaxed fitting: the structure is
//               relaxed (atoms and the six strains) at every evaluation and its cell lengths, angles and coordinates are
//               fitted to the observed ones (Gale 1996: the gradients alone weigh the curvature wrongly; relaxed fitting
//               fits the displacements themselves).
//  Properties   lattice energy per formula unit, cell lengths, angles and volume, elastic constants C_IJ (relaxed ions,
//               central differences of the stress over ± strain) and the Hill bulk and shear moduli.
//  Least squares  S = Σ w (calc − obs)², minimised by Levenberg–Marquardt (Marquardt, J. Soc. Ind. Appl. Math. 11, 431
//               (1963)) with a forward-difference Jacobian; standard errors from s² (JᵀJ)⁻¹ when there are more
//               observations than parameters (meaningful when each w is 1/σ²).
// Stress σ = −W/V (W the virial tensor Σ r_a f_b, tensile positive); Voigt order xx yy zz yz xz xy, engineering shear.
#pragma once

#include <array>
#include <functional>
#include <limits>
#include <string>
#include <vector>

#include "caps/field.hpp"
#include "caps/json.hpp"
#include "caps/system.hpp"

namespace caps {

struct FFDef;

struct FitParameter {
  std::string kind;                  // "pair", "lj" (a type's own ε, σ), "bond", "charge"
  std::vector<std::string> types;    // pair, bond: two type names; lj, charge: one
  // pair: A rho C (Buckingham, Rydberg), D alpha r0 (Morse), epsilon sigma (Lennard-Jones);
  // lj: epsilon sigma; bond: k r0 (harmonic), r0 k2 k3 k4 (class II), D alpha r0 (Morse); charge: q
  std::string name;
  double lower = -std::numeric_limits<double>::infinity();
  double upper = std::numeric_limits<double>::infinity();
  std::string balance;               // charge: this type's charge follows, so every structure keeps its total charge
  std::string label() const;         // "pair Mg O A"
};

struct FitTarget {
  size_t structure = 0;
  // a b c (Å), alpha beta gamma (degrees), volume (Å³ per cell), energy (kcal/mol per formula unit), C11 … C66 (GPa,
  // any of the 21), bulk shear (Hill, GPa)
  std::string kind;
  double value = 0, weight = 1;
};

struct FitStructure {
  std::string name;
  System cell;                                           // the observed structure, a periodic cell
  std::function<ForceField(const System&)> field_for;    // the force field of any supercell of it; required
  std::vector<char> shell;                               // per atom of the cell: a shell (empty: none)
  bool fit_structure = true;                             // the observed structure is a target
  double formula_units = 1;                              // energy targets are per formula unit
  std::array<int, 3> supercell{0, 0, 0};                 // 0: each width ≥ 2 × cut-off
};

struct FitOptions {
  std::vector<FitStructure> structures;
  std::vector<FitParameter> parameters;
  std::vector<FitTarget> targets;
  bool relaxed = false;              // relaxed fitting (else conventional: zero gradients at the observed structure)
  EnergyOptions energy;
  // weights of the structure terms: per Å² (lengths, coordinates), per degree² (angles), per (kcal/mol/Å)² (forces),
  // per GPa² (stresses)
  double w_length = 100, w_angle = 1, w_coord = 100, w_force = 1, w_stress = 1;
  int iterations = 50;
  double tolerance = 1e-6;           // stop when an accepted step lowers S by less than this fraction
  double fd_step = 1e-4;             // relative step of the Jacobian's differences
  double relax_ftol = 1e-6;          // kcal/mol/Å (and the strain gradient per cell length)
  int relax_iterations = 5000;
  double elastic_strain = 1e-4;
  size_t max_atoms = 15000;          // the automatic supercell stays below this
  std::function<bool(const std::string&, double)> progress;   // false stops
};

struct FitObservation {
  size_t structure = 0;
  std::string what, unit;
  double target = 0, initial = 0, final = 0, weight = 0;
};

struct FitStructureResult {
  std::string name;
  std::array<int, 3> supercell{1, 1, 1};
  size_t supercell_atoms = 0;
  std::array<double, 6> observed{}, initial{}, final{};   // a b c (Å) alpha beta gamma (degrees)
  double volume_observed = 0, volume_final = 0;
  double energy_initial = 0, energy_final = 0;            // kcal/mol per formula unit
  double shift_initial = 0, shift_final = 0;              // Å, RMS of the atoms' (not shells') displacement from observed
  double force_initial = 0, force_final = 0;              // kcal/mol/Å, largest at the observed structure (conventional)
  std::array<double, 6> stress_initial{}, stress_final{}; // GPa at the observed structure (conventional)
  bool elastic = false;
  std::array<std::array<double, 6>, 6> C{};               // GPa, final
  double bulk = 0, shear = 0;                             // Hill, GPa, final
  Cell cell;                                              // final cell (relaxed fitting) or the observed one
  std::vector<Vec3> positions;                            // final positions of the cell's atoms
  int relax_steps = 0;                                    // of the final relaxation
};

// A fitted term in full (the fitted parameter with the others of its term), for writing back.
struct FittedTerm {
  std::string kind;                  // pair, lj, bond, charge
  std::vector<std::string> types;
  std::string style;                 // buck, morse, rydberg, lj/cut; harmonic, class2, morse; "" for lj and charge
  std::vector<double> values;        // in the force-field file's order
};

struct FitResult {
  std::vector<std::string> names;
  std::vector<double> initial, final, error;   // error: standard error (NaN when not estimated)
  std::vector<FitObservation> observations;
  std::vector<FitStructureResult> structures;
  std::vector<FittedTerm> terms;
  double ss_initial = 0, ss_final = 0;
  std::vector<double> history;                 // S after each iteration
  int residuals = 0, iterations = 0, evaluations = 0;
  bool converged = false;
  std::string stopped_by;
  std::vector<std::string> notes;
  double seconds = 0;
};

// Fit the parameters (none: evaluate the structures and targets once).
FitResult fit_to_crystals(const FitOptions& o);

// Writes the fitted terms into a force-field definition: an explicit rule of the same types is updated, else a rule
// for exactly those types is appended (the last matching rule wins); fitted charges go on the types. Returns notes.
std::vector<std::string> apply_fitted_terms(FFDef& def, const std::vector<FittedTerm>& terms);

// The fit's settings from JSON (everything but the structures): "parameters": [{kind, types, name, lower, upper, balance}],
// "targets": [{structure, kind, value, weight, unit}] (an energy's unit: kcal/mol, kJ/mol or eV, per formula unit),
// "relaxed", "iterations", "tolerance", "fd_step", "relax_ftol", "elastic_strain",
// "weights": {length, angle, coord, force, stress}.
void fit_options_from_json(const Json& j, FitOptions& o);
// The report as JSON (parameters, observations, structures with their cells and elastic constants, terms, notes).
Json fit_result_json(const FitResult& r);
// The fittable terms of a force field (a structure's): pairs between the types present (their form and parameter names),
// the types' own Lennard-Jones, bonds and charges, as {kind, types, style, names, values}.
Json fittable_terms(const ForceField& ff);

// Cell lengths and angles (a b c in Å, alpha beta gamma in degrees).
std::array<double, 6> cell_parameters(const Cell& c);

}  // namespace caps
