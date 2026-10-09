#pragma once
// Reading VASP output (the DFT surface & adsorption workbench): OUTCAR, OSZICAR and vasp.out (OSZICAR is written in
// large buffered blocks; the more complete of the two is used), XDATCAR, DOSCAR (LORBIT 11), LOCPOT/CHGCAR grids,
// Bader's ACF.dat; and the checks built on them: a stage's status, a running job's progress with an ETA, and the
// physical health of an adsorption set.
#include <array>
#include <string>
#include <vector>

#include "caps/json.hpp"
#include "caps/system.hpp"

namespace caps {

struct Outcar {
  bool exists = false;
  bool finished = false;          // "General timing and accounting"
  bool relaxed = false;           // "reached required accuracy"
  bool ediff_reached = false;     // "EDIFF is reached" (a single point converged)
  bool scf_failed = false;        // "was not achieved"
  int nions = 0, nelm = 0, nsw = -1, ibrion = -100;
  std::vector<double> e0;         // energy(sigma->0) of every ionic step (eV)
  std::vector<std::vector<Vec3>> forces;   // TOTAL-FORCE blocks (eV/Å)
  std::vector<double> fmax;       // largest force per block
  std::vector<std::array<double, 6>> stress;   // "in kB" XX YY ZZ XY YZ ZX
  std::vector<double> magnetization, efermi, dipole;   // per step (dipole: z component, e·Å)
  std::vector<std::string> warnings;
  double seconds = 0;             // Elapsed time
  std::vector<double> freq_real, freq_imag;   // cm⁻¹ (IBRION 5–8)
  double encut = 0;
};
Outcar read_outcar(const std::string& path);

struct IonicStep { double f = 0, e0 = 0; int scf = 0; double temperature = 0, etot = 0; };
// the F= lines (and MD's T= E= F=) with the SCF iterations before each
std::vector<IonicStep> read_oszicar(const std::string& path);
// OSZICAR or vasp.out of a folder, whichever has more steps
std::vector<IonicStep> read_ionic_steps(const std::string& dir, std::string* used = nullptr);

std::vector<System> read_xdatcar(const std::string& path);

struct Doscar {
  double efermi = 0;
  std::vector<double> energy;     // E − E_F (eV)
  std::vector<double> total;      // states/eV (spin channels summed)
  std::vector<std::vector<double>> site;   // per atom, every orbital and spin summed
};
Doscar read_doscar(const std::string& path);

struct VaspGrid {
  System system;
  std::array<int, 3> n{0, 0, 0};
  std::vector<double> values;     // as in the file, x fastest
  double at(int i, int j, int k) const { return values[size_t(i) + size_t(n[0]) * (size_t(j) + size_t(n[1]) * size_t(k))]; }
};
// LOCPOT (values in eV) or CHGCAR (values ρ·V; divide by the volume for e/Å³): the first grid of the file
VaspGrid read_vasp_grid(const std::string& path);
void write_chgcar(const VaspGrid& g, const std::string& path);

std::vector<double> read_acf(const std::string& path);   // Bader electrons per atom (column 5)

// species line and counts of a POSCAR/CONTCAR
std::vector<std::pair<std::string, int>> poscar_species(const std::string& path);

// ---------------------------------------------------------------- checks
// One stage folder: finished, converged (relaxation or single point), SCF failure, ionic steps, final E0, max force,
// magnetisation, warnings, the in-plane lattice change of a cell stage (02_cell2 must move a by < 0.002 Å), and the
// CONTCAR re-validated (slabs). Folders of a case (01_cell … 04_static) are expanded.
Json check_runs(const std::vector<std::string>& dirs, bool complex, const std::string& data_dir);
// A case folder's job: stages done and running (from the newest slurm-*.out), the live files in the workspace it names,
// steps, last energies, max force against the target, time per step, and an ETA from a fit of log(fmax).
Json run_progress(const std::string& case_dir);
// Every complex of an adsorption set against E_ref = E(slab) + E(molecule): ENERGY (|ΔE| > 3 eV), JUMP (> 2 eV
// between steps), DIPOLE (> 1 e·Å), SCF (NELM hit), GEOM (atoms < 0.8 Å, molecule bonds broken, desorbed > 6 Å, below
// the termination plane), HANG (> 20 min without output while running) — each with what to do.
Json health_check(const std::string& set_root);

// ETA of a relaxation from its max forces: steps left to reach target by a straight-line fit of log(fmax) over the
// last 8 steps (−1 when the force is not falling)
int steps_to_target(const std::vector<double>& fmax, double target = 0.01);

}  // namespace caps
