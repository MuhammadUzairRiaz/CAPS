#pragma once
// VASP calculation sets (the DFT surface & adsorption workbench): INCAR stages with a reason for every tag, Γ-centred
// k-meshes of one density, POTCARs assembled from the user's licensed PAW directory (never shipped) and checked by
// their TITEL lines, parallelisation and job scripts from editable cluster profiles (stages chained, resume, live backup,
// hang watchdog, scratch workspaces, long-term storage with links), the job tools (submit list, update, reset, cleanup,
// store), convergence and lattice-scan sets with their collectors, and derived sets (charge, charge-density difference,
// partial-Hessian frequencies, AIMD).
#include <map>
#include <string>
#include <vector>

#include "caps/adsorb_dft.hpp"
#include "caps/json.hpp"
#include "caps/system.hpp"

namespace caps {

// ---------------------------------------------------------------- settings shared by every run that is compared
struct VaspSettings {
  double encut = 0;               // eV; 0 = 1.3 × the largest ENMAX of the POTCARs (520 without POTCARs)
  std::string prec = "Accurate", algo = "Normal";
  int nelm = 120;
  int ismear = 0;
  double sigma = 0.05;
  std::string gga = "PE";
  int ivdw = 12;                  // 0 none, 11 D3 zero damping, 12 D3-BJ, 13 D4 …
  bool lasph = true;
  int lmaxmix = 0;                // 0 = 4 with d elements (Sc–Zn, Y–Cd, Hf–Hg), else 2
  int ispin = 1;
  std::string lreal;              // "" = .FALSE. up to 40 atoms, else Auto (one value for every run of a set)
  std::string ediff_relax = "1E-6", ediff_static = "1E-7";
  std::string static_out = "all"; // 04_static: all (CHGCAR, AECCAR, LOCPOT), pot (LOCPOT), none (energy + DOS)
  int dipole = -1;                // −1 = when the faces differ (validator); 0 off; 1 on (IDIPOL = 3 energy-only, no LDIPOL)
  double kdens = 45.0;            // Å: N_i = ceil(kdens / |a_i|), 1 along the vacuum
  int nsw_cell = 200, nsw_relax = 0;   // 0 = 300 up to 100 atoms, else 600
};

struct ClusterRule { int max_atoms = 0; int ranks = 0, kpar = 1, ncore = 1; std::string time; };
struct ClusterProfile {
  std::string name, about, scheduler = "slurm", account, licenses, modules, vasp = "vasp_std";
  int cores_per_node = 0, mem_per_cpu = 0;
  std::string ws_allocate, ws_find, ws_release;
  std::vector<ClusterRule> parallel;   // first rule whose max_atoms (0 = any) is not exceeded
  ClusterRule gamma;
  std::vector<std::pair<int, int>> hang;   // (max_atoms or 0, seconds)
};
std::vector<ClusterProfile> cluster_profiles(const std::string& data_dir);   // the shipped ones, then ~/CAPS/dft-clusters.json
std::string user_clusters_file();
ClusterProfile cluster_profile(const std::string& data_dir, const std::string& name);
ClusterRule parallel_rule(const ClusterProfile& p, int atoms, bool gamma_only);
int hang_seconds(const ClusterProfile& p, int atoms);

// ---------------------------------------------------------------- the files
// "cell" (01_cell, 02_cell2), "relax" (03_relax), "static" (04_static); dipole_z < 0: no dipole correction
struct IncarLine { std::string key, value, why; };
std::vector<IncarLine> incar_lines(const std::string& kind, int atoms, bool gamma_only, double encut, double dipole_z, const VaspSettings& s,
                                   const ClusterRule& par, const std::vector<int>& elements);
std::string incar_text(const std::string& system, const std::vector<IncarLine>& lines);
// The reason for a tag (Studio tooltips, --help)
std::string incar_why(const std::string& key);

std::array<int, 3> kmesh(const System& s, double kdens, bool gamma_only = false);
std::string kpoints_text(const System& s, double kdens, bool gamma_only = false);

struct PotcarInfo {
  std::vector<std::string> species;    // POSCAR order
  std::vector<std::string> datasets;   // folder names used
  std::vector<std::string> titel;      // TITEL lines read
  std::vector<double> enmax, zval;
  double encut = 0;                    // 1.3 × max ENMAX
};
std::map<std::string, std::string> potcar_map(const std::string& data_dir);
// Concatenates pp_dir/<dataset>/POTCAR in species order into out (when out is not ""), reading ENMAX, ZVAL and the
// TITEL lines; throws when a file is missing or the TITEL order differs from the species.
PotcarInfo assemble_potcar(const std::vector<std::string>& species, const std::string& pp_dir, const std::map<std::string, std::string>& map, const std::string& out = "");
// Reads an existing POTCAR: TITEL element order (throws if it differs from `species`), ENMAX, ZVAL.
PotcarInfo verify_potcar(const std::string& potcar, const std::vector<std::string>& species);
std::string make_potcar_script(const std::map<std::string, std::string>& map);

std::string job_script(const ClusterProfile& p, const std::string& name, const std::vector<std::string>& stages, int atoms, bool gamma_only);
std::string aimd_job_script(const ClusterProfile& p, const std::string& name, int segments, int nsw);

// ---------------------------------------------------------------- sets
const std::vector<std::string>& stages_slab();    // 01_cell 02_cell2 03_relax 04_static
const std::vector<std::string>& stages_fixed();   // 03_relax 04_static

struct VaspSetOptions {
  std::string label;
  std::vector<std::string> stages;      // empty: stages_fixed()
  VaspSettings settings;
  std::string profile = "slurm-workspace";
  std::string pp_dir;                   // the licensed PAW directory ("" = make_potcar.sh only)
  bool gamma_only = false;
  bool asymmetric = false;              // the validator found different faces (dipole −1 then switches it on)
  std::string report_text;              // written as report.txt (the validator's)
  Json report_json;                     // and report.json
  std::string data_dir;
};
struct VaspSetResult {
  std::string dir;
  std::vector<std::string> files;
  double encut = 0;
  std::array<int, 3> kpoints{1, 1, 1};
  ClusterRule parallel;
  std::vector<std::string> notes;
};
// One case folder: POSCAR (species grouped, "<formula> CAPS <label>"), INCAR.cell / .relax / .static for its stages,
// KPOINTS, job.slurm, make_potcar.sh, POTCAR (with pp_dir), report.txt/json.
VaspSetResult write_vasp_set(const System& s, const std::string& dir, const VaspSetOptions& o);
// slab/, molecule/ (Γ only), complex_* — one settings set for all (LREAL Auto, EDIFF 1E-5 for the relaxations, the
// dipole correction everywhere) so the energies can be subtracted; failing complexes are not written unless forced.
std::vector<VaspSetResult> write_adsorption_set(const AdsorbSet& set, const std::string& root, VaspSetOptions o, bool force = false);

// convergence: conv_encut_<E>/ and conv_k_<N>/ single points from a case folder; the collector picks the cheapest
// setting for which it and every more expensive one pass (ENCUT: stress ≤ 1 kB, force ≤ 0.01 eV/Å; k: E0 ≤ 1 meV/atom,
// stress ≤ 1 kB, against the most expensive).
std::vector<std::string> write_convergence_set(const std::string& case_dir, const std::vector<int>& encuts = {400, 450, 520, 600, 700, 800},
                                               const std::vector<int>& kmeshes = {6, 9, 12, 15, 18});
struct ConvRow { std::string value; double e_atom = 0, stress = 0, fmax = 0; double de_mev = 0, ds = 0, df = 0; bool pass = false; };
struct ConvResult { std::vector<ConvRow> encut, k; std::string encut_choice, k_choice; };
ConvResult collect_convergence(const std::string& case_dir);
// ...and the same choice from given rows (value, E0/atom, in-plane stress, max force), most expensive last
std::string convergence_choice(std::vector<ConvRow>& rows, bool is_encut);

// lattice-constant scan: ascan_<a>/ fixed-cell relaxations; fit: parabola minimum, and Birch–Murnaghan (3rd order,
// linear in V^(−2/3)) on the cell volume
std::vector<std::string> write_lattice_scan(const std::string& case_dir, double lo, double hi, double step);
struct ScanFit { std::vector<std::pair<double, double>> points; double a_parabola = 0, a_bm = 0, v0 = 0, b0_gpa = 0; bool inside = false; };
ScanFit fit_lattice_scan(const std::string& case_dir);
ScanFit fit_scan_points(const std::vector<std::pair<double, double>>& a_e, const System& cell_at_a0, double a0);

// derived sets from a finished complex
std::string write_charge_set(const std::string& complex_dir);                       // <cx>/charge: all outputs
std::vector<std::string> write_cdd_set(const std::string& complex_dir, const System& slab_poscar, const std::string& pp_dir, const std::map<std::string, std::string>& map);
// partial Hessian: Selective dynamics frees `free_atoms` (empty: the nitrile N, its C, the α-C and its H)
std::string write_freq_set(const std::string& case_dir, std::vector<int> free_atoms = {});
std::string write_aimd_set(const System& start, const std::string& out, const VaspSetOptions& o, double temperature = 300, int nsw = 1000, int segments = 6);

// ---------------------------------------------------------------- job tools (run on the cluster's files)
struct ToolAction { std::string path, what; long long bytes = 0; };
// stage folders below `root` that a job may still run (not DONE), skipping `queued` working directories
std::vector<std::string> submittable(const std::vector<std::string>& case_dirs, const std::vector<std::string>& queued);
int update_job_scripts(const std::vector<std::string>& case_dirs, const ClusterProfile& p);   // keeps name, ranks, time, stages, hang
std::string reset_stage(const std::string& case_dir, const std::string& stage);              // renames the real folder <stage>_bad[n]
std::vector<ToolAction> cleanup_plan(const std::string& root, bool keep_best);               // only finished stages
std::vector<ToolAction> store_plan(const std::string& root, const std::string& store);       // stage folders → storage + links
void apply_plan(const std::vector<ToolAction>& plan, const std::string& kind, const std::string& root = "", const std::string& store = "");

}  // namespace caps
