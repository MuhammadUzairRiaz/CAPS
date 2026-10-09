// CAPS bonded coarse-grained potentials by tabulated Boltzmann inversion (the coarse-graining workflow's "Bonded" stage):
// bead–bead distributions gathered over many frames of many systems mapped alike (caps/cg_rules.hpp), each system's
// histograms normalised and weighted before they are pooled, then inverted per type
//
//   bonds       U(r) = −k_B T ln[P(r) / r²]
//   angles      U(θ) = −k_B T ln[P(θ) / sin θ]
//   dihedrals   U(φ) = −k_B T ln P(φ)                (IUPAC: trans = 180°)
//
// (Tschöp, Kremer, Batoulis, Bürger & Hahn, Acta Polym. 49, 61 (1998); Reith, Pütz & Müller-Plathe, J. Comput. Chem. 24,
// 1624 (2003)). Only where the smoothed distribution exceeds a fraction of its maximum (5 % by default) is it inverted;
// outside, the potential continues with its edge slope plus a stiff quadratic wall (value and slope continuous), so sparse
// tails — distorted contacts of an unrelaxed start — cannot leave soft spots that let beads collapse. Dihedrals are
// periodic: unsampled ranges take a smooth barrier through the gap, joined at both edges.
//
// Each table also carries how well it is sampled: the samples split in two halves — odd and even molecules, or the frames'
// halves for a single molecule — each inverted alone, and the largest difference over the inverted range (k_B T). Iterative refinement (bonded IBI): with the same distributions from a CG run,
// U ← U + α k_B T ln(P_CG / P_target) where both are sampled, which corrects for the non-bonded 1–3 and 1–4 terms that
// shift the bonded distributions.
//
// Outputs: LAMMPS bond_style table, angle_style table (degrees) and dihedral_style table/cut (with its aat switch, so the
// dihedral turns off as a neighbouring angle straightens), GROMACS tabulated bonded files (table_b<n>.xvg, table_a<n>.xvg,
// table_d<n>.xvg), the coefficient lines numbered by the shared type list, and the distributions (JSON) for plots.
#pragma once
#include <array>
#include <cstdint>
#include <map>
#include <string>
#include <vector>

#include "caps/cg_rules.hpp"
#include "caps/json.hpp"
#include "caps/system.hpp"

namespace caps {

// The bead topology a mapping file describes (no all-atom structure needed): kinds, bonds, angles, dihedrals, molecules.
struct CgTopology {
  std::vector<std::string> kind;
  std::vector<int> mol;
  std::vector<std::pair<int, int>> bonds;
  std::vector<std::array<int, 3>> angles;
  std::vector<std::array<int, 4>> dihedrals;
  std::vector<std::vector<int>> chains;
  size_t beads() const { return kind.size(); }
};
CgTopology cg_topology_from_map(const Json& map_json);
CgTopology cg_topology_of(const CgMapping& m);

struct CgBondedOptions {
  double temperature = 300;          // K
  double bond_bin = 0.02, bond_max = 20;      // Å
  double angle_bin = 1.0;            // degrees
  double dihedral_bin = 5.0;         // degrees
  double threshold = 0.05;           // invert where the smoothed P exceeds this share of its maximum
  double smooth = 1;                 // Gaussian smoothing, σ in bins at least (0: none); widened to Silverman's bandwidth
                                     // 0.9 sd n^(−1/5) when the samples are few
  double bond_wall = 20;             // kcal/mol/Å², the walls' curvature (at least the inverted well's)
  double angle_wall = 0.005;         // kcal/mol/deg²
  double dihedral_cap = 0;           // kcal/mol: the gap barrier's height above the sampled maximum (0: 2 k_B T)
  // tables
  double bond_table_lo = 0;          // Å (0: the sampled minimum − 2 Å, at least 0.5)
  double bond_table_hi = 0;          // Å (0: the sampled maximum + 3 Å)
  double bond_table_dr = 0.01;       // Å
  double angle_table_dt = 0.5;       // degrees (0 … 180)
  double dihedral_table_dp = 1.0;    // degrees (−180 … 180)
  double aat_k = 1.0, aat_theta1 = 170, aat_theta2 = 180;   // dihedral_style table/cut switch (degrees)
};

// One pooled histogram per type (weights already applied), with the halves for the sampling check.
struct CgBondedHistogram {
  int kind = 0;                      // 0 bond, 1 angle, 2 dihedral
  std::string key;
  double lo = 0, bin = 0;            // first edge, bin width (Å or degrees)
  std::vector<double> h, half[2];    // pooled; first and second halves of each system's frames
  long count = 0;
};

class CgBondedAccumulator {
 public:
  CgBondedAccumulator(const CgTypes& types, const CgBondedOptions& o);
  // a system: its topology and weight (relative; each system's histograms are normalised before pooling)
  int add_system(const CgTopology& t, double weight = 1, const std::string& name = "");
  // a frame of system s (bead positions in bead order; the cell for minimum images, or none for unwrapped beads)
  void add_frame(int s, const std::vector<Vec3>& beads, const Cell& cell);
  // frames per system (the halves are split by frame count)
  void set_expected_frames(int s, size_t frames);
  std::vector<CgBondedHistogram> pooled() const;
  size_t frames() const;
  const CgTypes& types() const { return types_; }
  const CgBondedOptions& options() const { return o_; }
  std::vector<std::string> notes() const;

 private:
  struct Sys { CgTopology t; double weight = 1; std::string name; size_t frames = 0, expected = 0;
               std::vector<int> bond_t, angle_t, dih_t; std::vector<std::vector<double>> h[3][2]; };
  CgTypes types_;
  CgBondedOptions o_;
  std::vector<Sys> sys_;
  std::vector<std::string> notes_;
  size_t nbins(int kind) const;
};

struct CgBondedTable {
  int kind = 0;
  std::string key;
  std::vector<double> x, U, F;       // the table (Å or degrees; kcal/mol; −dU/dx per Å or per degree)
  std::vector<double> xh, P, Uinv;   // the histogram's centres, P (normalised, Jacobian removed), U where inverted (NaN elsewhere)
  double lo = 0, hi = 0;             // the inverted range
  double x0 = 0;                     // the potential's minimum
  double k_harmonic = 0;             // curvature of a parabola through the well (U < 2 k_B T): E = k (x − x0)², per Å² or deg²
  double half_diff = 0;              // largest difference between the halves' potentials over the inverted range (k_B T)
  long count = 0;
  bool sampled = false;
};
struct CgBondedResult {
  std::vector<CgBondedTable> bonds, angles, dihedrals;
  double temperature = 0;
  size_t frames = 0;
  std::vector<std::string> notes;
};

CgBondedResult invert_bonded(const CgBondedAccumulator& acc);

// One bonded IBI step: `current` the tables in use, `cg` the same distributions sampled from a CG run with them;
// U ← U + α k_B T ln(P_cg / P_target) where both are above the threshold, walls rebuilt. Each table's residual
// (∫|P_cg − P_target| over the range) goes to `residual` by key.
CgBondedResult refine_bonded(const CgBondedResult& current, const CgBondedAccumulator& cg, double alpha = 0.5,
                             std::map<std::string, double>* residual = nullptr);

// Files in dir: bonds.table, angles.table, dihedrals.table (LAMMPS), bonded.in (styles and coefficients numbered by
// `types`, ${BONDED} for the folder), gromacs/table_b<n>.xvg, table_a<n>.xvg, table_d<n>.xvg (n = type number − 1) and
// bonded.json (every distribution and potential). Returns the files written.
std::vector<std::string> write_bonded(const CgBondedResult& r, const CgTypes& types, const CgBondedOptions& o, const std::string& dir);
Json bonded_json(const CgBondedResult& r);
// Tables read back from bonded.json (for a refinement step).
CgBondedResult bonded_from_json(const Json& j);

// The internal coordinates of a bonded term (Å; degrees; IUPAC dihedral, trans = 180°).
double cg_angle_deg(const Vec3& a, const Vec3& b, const Vec3& c);
double cg_dihedral_deg(const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& d);

}  // namespace caps
