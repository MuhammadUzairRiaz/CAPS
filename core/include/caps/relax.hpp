// CAPS Relax: energy minimisation, soft push-off, compression to a target density and box relaxation.
#pragma once
#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "caps/constraints.hpp"
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
  // Dihedral restraints: k (φ − φ0)² on the i-j-k-l dihedral (IUPAC sign, φ0 in degrees, k in kcal/mol/rad²): hold a
  // backbone torsion trans or gauche, set a side group's orientation while the rest relaxes.
  struct DihedralRestraint { uint32_t i = 0, j = 0, k = 0, l = 0; double phi0 = 180, kphi = 50; };
  std::vector<DihedralRestraint> dihedral_restraints;

  // Soft push-off before minimising (Auhl et al., J. Chem. Phys. 119, 12718 (2003)): LJ forces capped, the cap raised
  // stage by stage. Needed after Grow with a contact scale below 1, or for any structure with overlaps.
  bool pushoff = true;
  std::vector<double> pushoff_caps{5, 20, 100, 500};   // kcal/mol/Å
  // The same push-off by molecular dynamics first (Auhl et al.): NVT at pushoff_temperature for pushoff_ramp_ps, the
  // LJ force cap raised geometrically in pushoff_ramp_segments steps from the first cap to pushoff_cap (λ ramp), so
  // overlapping chains move apart by motion rather than by the steepest path down. 0 ps: minimisation stages only.
  // pushoff_cap > 0 also ends the minimisation stages there.
  double pushoff_ramp_ps = 0;
  double pushoff_cap = 0;         // kcal/mol/Å, 0 = the last of pushoff_caps
  double pushoff_temperature = 300;
  int pushoff_ramp_segments = 10;

  // Compression: scale the cell and all positions affinely, stage by stage, minimising after each stage.
  double target_density = 0.0;    // g/cm³, 0 = keep the cell
  double compress_step = 0.06;    // relative density increase per stage

  // Isotropic box relaxation towards a pressure (0 K mechanical equilibrium, not a thermal density).
  bool relax_box = false;
  double pressure = 1.0;          // atm
  double pressure_tol = 100.0;    // atm
  int box_cycles = 40;
  // Anisotropic: each Cartesian axis with box_axes[k] set follows its own diagonal pressure P_kk (a film or slab: z
  // alone, or x and y at a fixed thickness); the others keep their length. Orthorhombic cells.
  bool box_anisotropic = false;
  bool box_axes[3] = {true, true, true};

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
// A minimisation that ends with an angle stuck at 180° where its force field wants a bent one (an sp3 centre pushed flat
// through, as an opened ring or a moved hydrogen can leave it) sits on a saddle: the angle's force is zero by symmetry there
// and the minimiser stops, but any motion throws the atoms apart. relax() finds such angles afterwards, moves each centre
// 0.15 Å off the line and minimises again (up to three times), and says so in the notes.
void relax(System& s, const RelaxOptions& o, RelaxReport* report = nullptr);
// The vertices of angle terms with θ0 below 150° held above min_deg; each displaced 0.15 Å perpendicular to its two
// partners, in a direction drawn from seed; atoms held in place (fixed) never move. Returns how many were moved.
int kick_linear_angles(System& s, const ForceField& ff, uint64_t seed, double min_deg = 172.0, const std::vector<char>* fixed = nullptr);

// One minimisation with a given evaluator; positions x (3N) are updated. Exposed for tests and benches.
RelaxStage minimise(Evaluator& ev, std::vector<double>& x, const Cell& cell, const RelaxOptions& o, const std::string& name,
                    int* evaluations = nullptr);

// Largest atomic force magnitude.
double max_force(const std::vector<double>& f);

// LAMMPS data file with the force field (units real, atom_style full): every term CAPS evaluates, in the LAMMPS style
// with the same energy; kinds that mix forms (class II with class I, Morse with harmonic bonds, ...) become hybrid
// styles with "skip" lines in the class II sections. The header lists the LAMMPS commands that go with it. Throws
// FieldError for terms LAMMPS cannot reproduce exactly (separate 1-4 LJ parameters, non-planar cvff impropers).
// pair_coeffs false: the pair coefficients go in the input script instead (write_lammps_input with pair_coeffs), and
// the data file holds the structure, masses and bonded coefficients only.
// The description a CAPS file's title line is written with: an old CAPS header line read back loses its program,
// style and force-field parts ("structure" when nothing is left).
std::string export_title(std::string title, const std::string& ffname);

// How the LAMMPS files are written. exact (native false): the styles that give exactly CAPS's energy (damped shifted
// force or Ewald as CAPS computes, torsions as Fourier sums): for checking CAPS against LAMMPS. native: the force
// field's own styles, as its file declares them, for production runs: its dihedral style (OPLS K1–K4, CHARMM), long-range
// Coulomb by PPPM, its own cut-off. hybrid: every style written as "hybrid <sub-style>", each coefficient line naming its
// sub-style (moltemplate's layout); styles that need more than one sub-style are hybrid whatever this says.
struct LammpsStyle {
  bool native = false;
  bool hybrid = false;
  std::string coulomb = "auto";   // native: auto (PPPM when the force field is long-range and the cell periodic), pppm, ewald, dsf, cut
  double cutoff = 0;              // native: Å; 0 the force field's
  double kspace_accuracy = 1e-4;  // native: relative accuracy of PPPM / Ewald
  int tail = -1;                  // native: 1 pair_modify tail yes, 0 no, −1 as the energy options say
  // > 0: the data file holds only the first write_atoms atoms (and the terms among them) while its type and coefficient
  // tables cover the whole system — a fix bond/react cell whose post-reaction types come from a reacted copy appended
  size_t write_atoms = 0;
  // LAMMPS units: real (kcal/mol, fs, atm), metal (eV, ps, bar: every energy parameter divided by 23.060549, LAMMPS's own
  // factor), or auto: metal when the force field has a many-body potential LAMMPS reads in metal units only (AIREBO,
  // REBO), else real
  std::string units = "auto";
  // Named atom sets written as LAMMPS groups (a composite's filler and matrix; the Field page's groups): by atom type
  // when their types are their own (so per-group pair styles, fixes and computes can name them), else by molecule or
  // atom id, said in a comment. The input says where the first group's types are numbered.
  struct Group { std::string name; std::vector<uint32_t> atoms; };
  std::vector<Group> groups;
  // Rigid bodies (filler particles, rigid molecules): these molecule ids move as rigid bodies in the LAMMPS run (group
  // rigid, pairs inside a body excluded, fix rigid/nvt/small molecule; the others integrate on their own, and an NPT
  // barostat dilates only them). CAPS's own runs have no rigid bodies.
  std::vector<int64_t> rigid_mols;
  // Four-site water (TIP4P family): set by the writers when they leave the M sites out; LAMMPS's tip4p styles place
  // them d_OM (qdist, Å) along the bisector from the O of type o_type, with its H type, O–H bond and H–O–H angle types
  double tip4p_qdist = 0;
};

// Whether a LAMMPS export of this force field is written in metal units (see LammpsStyle::units); throws when real units
// are asked for a potential LAMMPS reads in metal units only.
bool lammps_metal_units(const ForceField& ff, const LammpsStyle& style = {});
// The force field with every energy parameter in eV (kcal/mol ÷ 23.060549): pair ε and the other pair forms' energies,
// bond, angle, torsion, improper and class II cross-term constants, DREIDING hydrogen bonds, Stillinger–Weber ε.
ForceField forcefield_in_metal_units(const ForceField& ff);

// The groups of LammpsStyle::groups as LAMMPS group commands, with the comment saying how the types are numbered.
std::string lammps_group_lines(const System& s, const ForceField& ff, const std::vector<LammpsStyle::Group>& groups);
void write_lammps_data_ff(const System& s, const ForceField& ff, const EnergyOptions& e, const std::string& path, bool pair_coeffs = true,
                          const LammpsStyle& style = {});
// The numbering a LAMMPS data file of s would use (1-based): each atom's type, and every bond, angle, dihedral and
// improper with its type — for molecule templates that must agree with the data file.
struct LammpsTerms {
  std::vector<int> atom_type;
  std::vector<std::pair<int, std::vector<uint32_t>>> bonds, angles, dihedrals, impropers;
  std::vector<std::string> type_names;
};
LammpsTerms lammps_terms(const System& s, const ForceField& ff, const EnergyOptions& e, const LammpsStyle& style = {});
// The same, or the structure alone (atoms, types, charges, bonds; no coefficients) when LAMMPS has no form for the force
// field (Martini 3's reaction field and virtual sites): returns why, or "" when the coefficients were written.
std::string write_lammps_data_or_structure(const System& s, const ForceField& ff, const EnergyOptions& e, const std::string& path);

// How a LAMMPS input ends. Check: a single point with every energy term printed (the parity benches). Otherwise a
// protocol: an optional conjugate-gradient minimisation, then NVT or NPT (Nosé–Hoover, LAMMPS's standard; CAPS's own
// runs use Bussi and stochastic cell rescaling), with thermo output, a dump and the final structure written.
struct LammpsRun {
  // Tensile: the box stretched along axis at a constant engineering strain rate (fix deform erate), the other two axes
  // at pressure (Nosé–Hoover), stress_strain.dat with the engineering strain and the tensile stress −P_axis (MPa).
  // Creep: a constant true tensile stress stress_mpa along axis (its pressure −σ, the others at pressure; NPT, axes
  // uncoupled), creep.dat with time (ps) and strain. Shear: planar Couette flow by SLLOD (fix nvt/sllod, the xy tilt at
  // shear_rate, remap v), viscosity.dat with η = −⟨P_xy⟩/γ̇ (mPa·s) averaged in blocks.
  // Protocol: a CAPS equilibration protocol (protocol: its text, as protocol_text writes it) as LAMMPS stages, each
  // its own fix and run, every length scaled by the input's ${scale} (default 1), then production_ps of NPT at the last
  // stage's temperature and pressure with the density averaged (density.dat).
  enum class Kind { Check, None, Minimize, NVT, NPT, Tensile, Creep, Shear, Protocol } kind = Kind::Check;
  std::string protocol;            // Protocol: the stages' text
  double production_ps = 0;        // Protocol: NPT after the stages (0: none)
  int axis = 0;                    // 0 x, 1 y, 2 z (tensile, creep)
  double strain_rate = 1e-3;       // 1/ps, engineering (tensile)
  double max_strain = 0;           // > 0: the tensile run's steps from it (strain / (rate·dt))
  double stress_mpa = 50;          // MPa, tension positive (creep)
  double shear_rate = 0.01;        // 1/ps (shear)
  bool minimize_first = true;
  double temperature = 300, pressure = 1.0;   // K, atm
  double dt = 0;                               // fs; 0: the force field's own (lammps_timestep)
  int64_t steps = 100000;
  double tdamp = 100, pdamp = 1000;            // fs
  int thermo_every = 1000, dump_every = 5000;
  uint64_t seed = 4928459;
  // bonds to hydrogen (and rigid water) or every bond held at its length: fix shake in LAMMPS, constraints in GROMACS
  ConstraintMode constraints = ConstraintMode::None;
};
// The time step engine inputs are written with: the run's when set, else the force field's own (Martini 20 fs), else 2 fs
// with constraints, else 0.5 fs.
double lammps_timestep(const LammpsRun& run, const ForceField& ff);

// A LAMMPS input script for that data file: units, styles, special_bonds, read_data (as data_path is given), the pair
// coefficients when pair_coeffs (every i-j pair written out: nothing left to LAMMPS's mixing), neighbour settings, then
// the run section. held_mol > 0: that molecule is held in place (group, zero velocity, fix setforce), as CAPS holds an
// interface's surface.
void write_lammps_input(const System& s, const ForceField& ff, const EnergyOptions& e, const std::string& data_path, const std::string& path,
                        int64_t held_mol = 0, bool pair_coeffs = false, const LammpsRun& run = {}, const LammpsStyle& style = {},
                        std::vector<std::string>* notes = nullptr);

// The LAMMPS fix shake line for these constraints on group (empty when there is nothing to hold), with the bond and angle
// type numbers of the data file written with the same style.
std::string lammps_shake_fix(const System& s, const ForceField& ff, const EnergyOptions& e, ConstraintMode mode, const std::string& group = "all",
                             const LammpsStyle& style = {});

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
