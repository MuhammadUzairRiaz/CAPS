// CAPS Relax: energy minimisation, soft push-off, compression to a target density and box relaxation.
#pragma once
#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "caps/field.hpp"
#include "caps/system.hpp"

namespace caps {

enum class Minimiser { SteepestDescent, ConjugateGradient, LBFGS, FIRE };

Minimiser minimiser_from_string(const std::string& s);   // "sd" | "cg" | "lbfgs" | "fire"
const char* to_string(Minimiser m);

struct RelaxProgress {
  std::string stage;       // "push-off 20", "compress 0.85", "minimise", "box"
  int stage_index = 0, stages = 0;
  int iteration = 0;
  double energy = 0, fmax = 0, density = 0, pressure = 0;
};

struct RelaxOptions {
  Minimiser method = Minimiser::LBFGS;
  double ftol = 0.5;              // stop when the largest atomic force is below this, kcal/mol/Å
  double etol = 1e-8;             // or when the relative energy change per step is below this
  int max_iterations = 5000;      // per minimisation
  EnergyOptions energy;
  // The force field to use (from CAPS Field); null: CAPS's built-in GAFF typing of C and H.
  std::shared_ptr<const ForceField> field;
  // Per atom: 1 = held in place during minimisation (a substrate under a film). Compression and box relaxation still
  // scale every atom.
  std::vector<char> fixed;
  // Distance restraints: k (r − r0)² between atoms i and j (minimum image), added to the energy that is minimised
  // (pull two groups to a contact distance, hold a hydrogen bond, close a gap). Reported at the end.
  struct Restraint { uint32_t i = 0, j = 0; double r0 = 0, k = 10; };   // Å, kcal/mol/Å²
  std::vector<Restraint> restraints;

  // Soft push-off before minimising (Auhl et al., J. Chem. Phys. 119, 12718 (2003)): LJ forces capped, the cap raised
  // stage by stage. Needed after Grow with a contact scale below 1, or for any structure with overlaps.
  bool pushoff = true;
  std::vector<double> pushoff_caps{5, 20, 100, 500};   // kcal/mol/Å

  // Compression: scale the cell and all positions affinely, stage by stage, minimising after each stage.
  double target_density = 0.0;    // g/cm³, 0 = keep the cell
  double compress_step = 0.06;    // relative density increase per stage

  // Isotropic box relaxation towards a pressure (0 K mechanical equilibrium, not a thermal density).
  bool relax_box = false;
  double pressure = 1.0;          // atm
  double pressure_tol = 100.0;    // atm
  int box_cycles = 40;

  // Called often; return false to cancel (throws RelaxCancelled).
  std::function<bool(const RelaxProgress&)> progress;
  // Called after each stage with the positions (3N) and cell, for a frame-by-frame record.
  std::function<void(const std::vector<double>&, const Cell&, const std::string&)> snapshot;
};

struct RelaxStage {
  std::string name;
  int iterations = 0;
  double energy = 0, fmax = 0, density = 0, pressure = 0;
  std::string stopped_by;   // "force", "energy", "iterations", "line search"
};

struct RelaxReport {
  std::string field;
  EnergyTerms initial, final;
  double fmax_initial = 0, fmax_final = 0;
  double density_initial = 0, density_final = 0, pressure_final = 0;
  int iterations = 0, evaluations = 0, list_builds = 0;
  std::vector<RelaxStage> stages;
  std::vector<std::string> notes;
  bool converged = false;
};

struct RelaxCancelled : std::runtime_error {
  RelaxCancelled() : std::runtime_error("relax cancelled") {}
};

// Relax `s` in place (positions, cell, and charges when they were computed). Throws FieldError when the structure
// cannot be typed, RelaxCancelled when the progress callback returns false.
void relax(System& s, const RelaxOptions& o, RelaxReport* report = nullptr);

// One minimisation with a given evaluator; positions x (3N) are updated. Exposed for tests and benches.
RelaxStage minimise(Evaluator& ev, std::vector<double>& x, const Cell& cell, const RelaxOptions& o, const std::string& name,
                    int* evaluations = nullptr);

// Largest atomic force magnitude.
double max_force(const std::vector<double>& f);

// LAMMPS data file with the force field (units real, atom_style full): every term CAPS evaluates, in the LAMMPS style
// with the same energy; kinds that mix forms (class II with class I, Morse with harmonic bonds, ...) become hybrid
// styles with "skip" lines in the class II sections. The header lists the LAMMPS commands that go with it. Throws
// FieldError for terms LAMMPS cannot reproduce exactly (separate 1-4 LJ parameters, non-planar cvff impropers).
void write_lammps_data_ff(const System& s, const ForceField& ff, const EnergyOptions& e, const std::string& path);
// A LAMMPS input script for that data file: styles, special_bonds, read_data, run 0 with every energy term.
// held_mol > 0: that molecule is held in place (group, zero velocity, fix setforce), as CAPS holds an interface's surface.
void write_lammps_input(const System& s, const ForceField& ff, const EnergyOptions& e, const std::string& data_path, const std::string& path,
                        int64_t held_mol = 0);

// GROMACS files with the force field: STEM.top (every term in the GROMACS function with the same energy; every
// Lennard-Jones type pair and 1-4 pair written out, CAPS's exclusions listed), STEM.gro (nm, 8 decimals; molecules made
// whole) and STEM.mdp (a single-point run with the matching cut-offs, tail and electrostatics). Returns notes where
// GROMACS cannot compute exactly the same (DSF becomes PME; no cell). Throws FieldError for forms GROMACS lacks (class
// II, inversions, 9-6 or Buckingham pairs).
std::vector<std::string> write_gromacs(const System& s, const ForceField& ff, const EnergyOptions& e, const std::string& stem);
// The notes alone (nothing written); throws as write_gromacs does.
std::vector<std::string> gromacs_notes(const System& s, const ForceField& ff, const EnergyOptions& e);
// The .mdp non-bonded settings alone (cut-offs, modifiers, dispersion correction, electrostatics).
std::string gromacs_mdp(const System& s, const ForceField& ff, const EnergyOptions& e);

}  // namespace caps
