// CAPS non-bonded coarse-grained potentials (the coarse-graining workflow's "Non-bonded" stage).
//
//   targets      g_ij(r) per bead-type pair and per system, from mapped all-atom frames, over the non-bonded bead pairs (the
//                model's own exclusions: 1-2 and 1-3 by default, LAMMPS special_bonds lj 0 0 1), normalised by the ideal
//                pair density N_i N_j / V
//   IBI          one table per pair, U₀ = −k_B T ln g_target (a straight wall below the first resolved point, zero at the
//                cut-off), refined from CG runs: U ← U + α k_B T ln(g_CG / g_target) (Reith, Pütz & Müller-Plathe,
//                J. Comput. Chem. 24, 1624 (2003)). Several systems are fitted jointly: a pair that occurs in several
//                (B–B in PBS, PBSA and PBAT) gets the average of their corrections weighted by each system's statistics
//                at that distance and its weight, so one table set serves them all. The pressure correction of the same
//                paper, ΔU = A (1 − r/r_c) with A = −0.1 k_B T sign(ΔP) min(1, 0.0003 |ΔP| / bar), uses the weighted mean
//                pressure error of the systems. Convergence per pair: ∫(g_CG − g_target)² dr / ∫ g_target² dr and the
//                largest |g_CG − g_target|.
//   analytic     a table fitted by LJ 12-6, LJ 9-6 (class2), Morse or Mie(n, m) over its well and the repulsion up to
//                a few k_B T (Levenberg–Marquardt), for models that must transfer beyond the state point
//   calibration  a fitted set scaled to a target density (all σ by one factor) and a target T_g (all ε by one factor), the
//                ratios between pairs kept, by secant steps on the CG runs' results
//
// The CG runs are LAMMPS's (pair_style table and the bonded tables of caps/cg_bonded.hpp: the decks and a loop script are
// written) or CAPS's own engine for small systems (bonded terms as harmonic springs from the tables' wells; dihedrals
// left out, said so).
#pragma once
#include <map>
#include <string>
#include <vector>

#include "caps/cg_bonded.hpp"
#include "caps/field.hpp"
#include "caps/json.hpp"

namespace caps {

// all unordered bead-type pairs of a type list, "A-A", "A-B", … (canonical: the smaller name first)
std::vector<std::string> cg_pair_keys(const CgTypes& t);

struct CgPairOptions {
  double rmax = 15;      // Å, the g(r) range (at most half the shortest cell edge is used)
  double dr = 0.05;      // Å
  int exclude = 3;       // bead pairs fewer than this many bonds apart are left out: 2 → 1-2, 3 → 1-2 and 1-3, 4 → also 1-4
};

struct CgRdf {
  std::string key;
  std::vector<double> g, counts;    // per bin: g(r) and the raw pair count (all frames)
  double ideal_pairs = 0;           // N_i N_j (or N_i (N_i − 1)/2) per frame
  long frames = 0;
};

class CgRdfAccumulator {
 public:
  CgRdfAccumulator(const CgTypes& types, const CgPairOptions& o);
  int add_system(const CgTopology& t, double weight = 1, const std::string& name = "");
  void add_frame(int s, const std::vector<Vec3>& beads, const Cell& cell);   // a periodic cell is needed
  std::vector<CgRdf> rdf(int s) const;                                      // per pair key, in cg_pair_keys order
  std::vector<double> r() const;                                            // bin centres
  int systems() const { return int(sys_.size()); }
  const std::string& name(int s) const { return sys_.at(size_t(s)).name; }
  double weight(int s) const { return sys_.at(size_t(s)).weight; }
  const std::vector<std::string>& keys() const { return keys_; }
  const CgPairOptions& options() const { return o_; }

 private:
  struct Sys { CgTopology t; double weight = 1; std::string name; std::vector<int> tix; std::vector<std::vector<int>> excl;
               std::vector<std::vector<double>> h; double volume = 0; long frames = 0; std::vector<double> ntype; };
  CgTypes types_;
  CgPairOptions o_;
  std::vector<std::string> keys_;
  std::vector<Sys> sys_;
  double rmax_used_ = 0;
};

// Targets per system (what IBI matches), saved and read back as targets.json.
struct CgTargets {
  std::vector<double> r;
  std::vector<std::string> systems;
  std::vector<double> weights;
  std::vector<std::vector<CgRdf>> rdf;   // [system][pair]
  std::vector<double> pressure;          // atm, the AA pressure of each system (the IBI target pressure)
  double temperature = 300;
};
Json cg_targets_json(const CgTargets& t);
CgTargets cg_targets_from_json(const Json& j);
CgTargets cg_targets(const CgRdfAccumulator& acc, double temperature, const std::vector<double>& pressure = {});

struct CgPairTable {
  std::string key;
  std::vector<double> U, F;              // on r0, r0 + dr, … rc (kcal/mol, kcal/mol/Å)
};
struct CgPairSet {
  double r0 = 0.5, dr = 0.05, rc = 15;   // Å, the table grid
  double temperature = 300;
  int iteration = 0;
  std::vector<CgPairTable> pairs;
  std::string source;                    // "IBI", "LJ 12-6 fit" …
  size_t points() const { return pairs.empty() ? 0 : pairs.front().U.size(); }
};
Json cg_pairs_json(const CgPairSet& p);
CgPairSet cg_pairs_from_json(const Json& j);

struct CgIbiOptions {
  double alpha = 0.2;                    // damping of each update
  double rc = 15;                        // Å (at most the targets' range)
  double table_dr = 0.05;                // Å
  double smooth = 1;                     // Gaussian smoothing of each update, σ in bins
  bool pressure_correction = true;
  double ramp = 0.1;                     // A = −ramp k_B T sign(ΔP) min(1, 0.0003 |ΔP|/bar)
};
// The start: −k_B T ln g of the targets pooled over the systems (weighted by their statistics).
CgPairSet ibi_start(const CgTargets& t, const CgIbiOptions& o);

struct CgIbiPairReport { std::string key; double residual = 0, max_dev = 0; };
struct CgIbiStepReport {
  std::vector<CgIbiPairReport> pairs;
  std::vector<double> pressure;          // atm, per system (the CG runs')
  double ramp_A = 0;                     // kcal/mol, the pressure correction's amplitude
  double residual = 0;                   // the largest pair residual
};
// One joint update from the CG runs' g(r) (current[s][pair], the targets' grid) and pressures (atm).
CgPairSet ibi_step(const CgPairSet& cur, const CgTargets& t, const std::vector<std::vector<CgRdf>>& current, const std::vector<double>& pressure,
                   const CgIbiOptions& o, CgIbiStepReport* rep = nullptr);
// The residual of CG g(r)s against the targets without updating (a converged check).
CgIbiStepReport ibi_compare(const CgTargets& t, const std::vector<std::vector<CgRdf>>& current, const std::vector<double>& pressure);

// Files: pairs.table (LAMMPS pair_style table, one section per pair) and pair.in (pair_style and pair_coeff numbered by the
// type list, ${PAIR} for the folder), pairs.json.
std::vector<std::string> write_pairs(const CgPairSet& p, const CgTypes& types, const std::string& dir);

// ---------------------------------------------------------------- analytic fits
struct CgPairFit {
  std::string key, form;                 // lj126 | lj96 | morse | mie
  double epsilon = 0, sigma = 0;         // kcal/mol, Å (Morse: D0, r0; and alpha in a)
  double a = 0, n = 12, m = 6;           // Morse α (1/Å); Mie exponents
  double rms = 0;                        // kcal/mol over the fitted range
  double lo = 0, hi = 0;                 // the fitted range
};
// fits each pair of a table set; the range: from where U falls below `repulsion` k_B T to the cut-off
std::vector<CgPairFit> fit_pairs(const CgPairSet& p, const std::string& form, double repulsion = 4.0);
double pair_energy(const CgPairFit& f, double r);
// the analytic set as tables on p's grid (for runs or a further IBI) and its LAMMPS lines
CgPairSet tables_of_fits(const std::vector<CgPairFit>& fits, const CgPairSet& grid);
std::string lammps_pair_lines(const std::vector<CgPairFit>& fits, const CgTypes& types, double rc);

// ---------------------------------------------------------------- thermomechanical calibration
// σ and ε scale factors from the runs so far (secant on density(σ-scale) and T_g(ε-scale), log-linear; the first step
// from the power laws ρ ∝ s_σ⁻³ and T_g ∝ s_ε).
struct CgCalibrationPoint { double s_sigma = 1, s_eps = 1, density = 0, tg = 0; };
struct CgCalibrationStep { double s_sigma = 1, s_eps = 1; bool density_done = false, tg_done = false; std::string note; };
CgCalibrationStep calibration_step(const std::vector<CgCalibrationPoint>& history, double target_density, double target_tg,
                                   double density_tol = 0.01, double tg_tol = 5.0);
std::vector<CgPairFit> scale_fits(const std::vector<CgPairFit>& fits, double s_sigma, double s_eps);

// ---------------------------------------------------------------- CAPS's engine (small systems, local checks)
// The bead force field: each pair on its table (pairs), bonds and angles as harmonic springs at the bonded tables' wells
// (k from a parabola through the well; dihedrals left out, said in the notes) — or none when the topology has no bonds;
// beads fewer than `exclude` bonds apart excluded (special_bonds lj 0 0 1 for 3). Masses: each kind's from `masses`.
ForceField cg_forcefield(const CgTopology& t, const CgTypes& types, const std::map<std::string, double>& masses, const CgPairSet& pairs,
                         const CgBondedResult* bonded = nullptr, int exclude = 3);
struct CgEngineOptions {
  double temperature = 300, dt = 5;     // K, fs
  int64_t steps = 2000, equilibrate = 500, frame_every = 50;
  bool npt = false;
  double pressure = 1;                   // atm (NPT)
  uint64_t seed = 1;
  bool new_velocities = true;
};
struct CgEngineRun {
  std::vector<std::vector<Vec3>> frames;
  std::vector<Cell> cells;
  double pressure = 0, density = 0;      // means over the sampled part
  System last;
};
CgEngineRun run_cg_engine(const System& beads, const ForceField& ff, const CgEngineOptions& o);

}  // namespace caps
