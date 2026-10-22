#pragma once
// Lattice thermal conductivity of a crystal with VASP and phono3py, the way the group VASP guide runs it (ch. 7–8): CAPS
// writes the inputs and the batch scripts of every step and reads the results; phono3py itself runs on the cluster
// (CAPS does not need it installed).
//
//  plan        The supercells (fc3 and fc2) from the cell heights by a stated rule, the number of displaced supercells
//              (the full finite-difference set counted with phono3py's own rules: site symmetry, least displacements
//              along the axes, +/- where the site symmetry does not supply it, second displacements per symmetry-reduced
//              pair; or N random displacements fitted with symfc), the k-meshes, the q-point mesh, and the cost in
//              core-hours BEFORE anything is submitted.
//  setup       POSCAR-unitcell, INCAR (single point: IBRION −1, ISIF 0, EDIFF 1E-8, PREC Accurate, LREAL .FALSE.), the
//              supercell KPOINTS, POTCAR, and run_step1.sh … run_step7.sh (displacements; fc3 and fc2 job arrays with the
//              vasprun.xml check and restart; check and quota cleanup; FORCES_FC3 / FORCES_FC2; fc3.hdf5 / fc2.hdf5;
//              kappa), run_dryrun.sh (a few displacements through the whole chain before any array, which also sets the
//              arrays' time limits from the measured time), cleanup.sh, analyze_kappa.py.
//  check       The displacement folders of a set: done, failed (vasprun.xml incomplete), not run; the array line to rerun.
//  analyze     kappa-m*.hdf5 (through analyze_kappa.py, h5py) or its extract: the κ tensor per temperature, along the
//              chain against the transverse mean, cumulative κ against frequency and mean free path, lifetimes and group
//              velocities, as CSV tables.
#include <array>
#include <string>
#include <vector>

#include "caps/json.hpp"
#include "caps/system.hpp"

namespace caps {

// The lattice vectors along which a covalent molecule is endless (a polymer crystal's chain axis): 0 = a, 1 = b, 2 = c;
// empty for a molecular crystal of finite molecules, three for a covalent network.
std::vector<int> chain_axes(const System& s);

// n_i = max(1, ceil(length × (1 − slack) / h_i)) with h_i the cell's height along axis i (the distance between the
// faces a_j × a_k): the smallest supercell whose images are at least `length` apart, 5 % short accepted.
std::array<int, 3> supercell_for_length(const System& s, double length, double slack = 0.05);

struct SymOpFrac { std::array<std::array<int, 3>, 3> R{}; Vec3 t{0, 0, 0}; };
// The space-group operations of a structure in its own fractional coordinates (rotations with entries −1…1 that keep
// the metric; translations from the atom mapping), tolerance in Å.
std::vector<SymOpFrac> crystal_operations(const System& s, double tol = 1e-3);

struct DisplacementCount {
  int first = 0;          // single displacements (phono3py's first supercells; also the fc2 set without --dim-fc2)
  int pairs = 0;          // second displacements over all first ones
  int fc3() const { return first + pairs; }   // supercells phono3py -d writes for fc3
  int fc2 = 0;            // supercells for fc2 (--dim-fc2)
  int independent_atoms = 0;
  int operations = 0;     // of the unit cell
};
// phono3py's full finite-difference set for these supercells; cutoff_pair > 0 (Å): pairs farther apart are left out, as
// phono3py --cutoff-pair leaves them out.
DisplacementCount count_displacements(const System& unit, const std::array<int, 3>& dim_fc3, const std::array<int, 3>& dim_fc2, double tol = 1e-3,
                                      double cutoff_pair = 0);

struct Phono3pyOptions {
  std::string label;
  std::string data_dir, pp_dir, profile = "tud-barnard", account;
  std::array<int, 3> dim_fc3{0, 0, 0}, dim_fc2{0, 0, 0};   // 0 = the length rule
  double fc3_length = 10.0, fc2_length = 15.0;             // Å
  std::string mode = "auto";      // full | random | both | auto (random when the full fc3 set exceeds random_above)
  int random_above = 3000;
  int rd = 0;                     // random fc3 supercells (0 = full set / 10, at least 50)
  bool rd_auto = false;           // phono3py --rd auto chooses the number (the estimate keeps full / 10)
  double cutoff_pair = 0;         // Å, full set: phono3py --cutoff-pair (0 = none)
  int random_seed = 1;
  double kdens = 50.0;            // Å, the supercells' k-meshes
  std::array<int, 3> mesh{0, 0, 0};                        // 0 = 74 Å / |a_i| across, 405 Å / |a_i| along the chain
  int chain_axis = -2;            // −2 = found from the bonds, −1 = none, 0–2
  std::vector<double> temperatures{100, 200, 300, 400};
  std::string method = "br";      // br | lbte
  double encut = 0;               // eV; 0 = relax4's (chain.json beside the structure), else 2 × ENMAX
  int ivdw = 12;
  int dry_run = 2;                // fc3 supercells of the dry run (plus one fc2)
  int array_limit = 0;            // %N on the arrays (0 = none)
  // 2D sheets: --dim / --dim-fc2 / k / q-mesh 1 along c, the dipole correction when polar, κ rescaled by h / thickness
  std::string dimensionality;     // "2d" | "3d" | "" (3d unless the structure has a vacuum gap: then an error)
  std::string dipole = "auto";    // 2d: auto | on | off
  double thickness = 0;           // Å, 2d: the sheet's effective thickness (e.g. the bulk parent's interlayer distance)
};

struct Phono3pyPlan {
  std::string label;
  std::string dimensionality = "3d";
  double cell_height = 0, thickness = 0;   // 2d: Å
  int atoms = 0, electrons = 0;
  std::array<int, 3> dim_fc3{1, 1, 1}, dim_fc2{1, 1, 1};
  int atoms_fc3 = 0, atoms_fc2 = 0;
  DisplacementCount full;
  std::string mode;               // full | random
  int fc3_supercells = 0;         // what this set computes (full count or --rd N)
  int rd = 0;
  std::array<int, 3> kpoints_fc3{1, 1, 1}, kpoints_fc2{1, 1, 1}, mesh{1, 1, 1};
  int chain_axis = -1;
  double encut = 0;
  double core_h_fc3 = 0, core_h_fc2 = 0;   // per supercell
  double core_h_total = 0;
  int ranks_fc3 = 0, ranks_fc2 = 0;
  std::vector<std::string> notes;
  Json json() const;
  std::string text() const;
};
// Valence electrons of the standard PAW datasets (H 1, C 4, N 5, O 6, F 7; data/dft/enmax.json "zval").
int valence_electrons(const System& s, const std::string& data_dir);
// The cost model: core-hours of one single point = 4.33 × (N_e / 288)² — calibrated on the guide's PE (a 144-atom fc3
// supercell, 288 valence electrons, about 10 min on 26 cores; it then gives its 432-atom fc2 supercells about 1.5 h, as
// the guide reports 1–2 h). An order of magnitude, not a promise: the dry run measures it.
double single_point_core_hours(int electrons);
Phono3pyPlan plan_phono3py(const System& unit, const Phono3pyOptions& o, const std::string& mode);

struct Phono3pyResult { std::vector<std::string> dirs; std::vector<Phono3pyPlan> plans; std::vector<std::string> notes; };
// One set per mode (both: full/ and random/ below dir).
Phono3pyResult write_phono3py_set(const System& unit, const std::string& dir, const Phono3pyOptions& o);

// The references for a set (keys of the built-in bibliography): VASP, PAW, PBE and D3(BJ), phono3py, symfc unless every
// displacement is computed (mode full), the direct LBTE solution for method lbte.
std::vector<std::string> phono3py_cites(const std::string& mode, const std::string& method);

// The displacement folders of a set (after run_step1.sh): done / failed / missing, and the array spec to rerun.
Json phono3py_check(const std::string& dir);

// Analysis of an extract written by analyze_kappa.py (kappa_summary.json + modes_T<T>.csv): per temperature the tensor
// (xx yy zz yz xz xy, W/m/K), along the chain against the transverse mean, cumulative κ against frequency and against
// mean free path (CSV), lifetimes τ = 1 / (4π Γ) and group velocities; writes the tables into out_dir. A 2D sheet:
// phono3py's κ is per cell volume (vacuum included), so the sheet conductance G = κ h (h the cell height along c) is
// reported, and with a thickness t also κ h / t, the sheet's own conductivity.
Json analyze_kappa_extract(const std::string& extract_dir, const std::string& out_dir, int chain_axis, const std::string& dimensionality = "3d",
                           double cell_height = 0, double thickness = 0);

}  // namespace caps
