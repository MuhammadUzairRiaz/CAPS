// CAPS entanglement, mechanics and dynamics analyses of coarse-grained melts (the coarse-graining workflow's "Analyse" stage).
//
// Entanglements
//   primitive paths    CAPS's own PPA (caps/entangle.hpp) on the beads, or a LAMMPS PPA deck for large systems (chain ends
//                      held, intra-chain pairs off, bonds with zero rest length, a quench) whose dump is read back here;
//                      Z1+ (the Z1+ package of Kröger, Dietz, Hoy & Luap, Mendeley Data (2022), the method of Kröger,
//                      Comput. Phys. Commun. 168, 209 (2005)): config.Z1 written, the program run when it is installed, its
//                      shortest paths read back (Z = interior nodes per chain)
//   estimators         Hoy, Foteinopoulou & Kröger, Phys. Rev. E 80, 031803 (2009), N beads per chain:
//                        classical S-coil  N_e = (N − 1) ⟨R²⟩ / ⟨L_pp⟩²                   (4)
//                        classical S-kink  N_e = N (N − 1) / (⟨Z⟩ (N − 1) + N)            (5)
//                        modified S-kink   N_e = N / ⟨Z⟩                                  (6)
//                        modified S-coil   N_e = (N − 1) (⟨L_pp²⟩ / ⟨R²⟩ − 1)⁻¹           (7)
//                        M-kink            1 / N_e = d⟨Z⟩ / dN over several chain lengths  (13)
//                        M-coil            C(x)/x at x = N_e equals d/dN (⟨L_pp⟩² / ((N − 1) l₀²)), C(x) the
//                                          characteristic ratio of a chain of x beads      (15)
//                      M_e = N_e × the mean bead mass; Z per chain = N / N_e
// Mechanics (LAMMPS decks, CAPS's analysis of what they print)
//   tension            uniaxial stress (lateral axes at P, NPT) or constant volume, z stretched at a constant engineering
//                      rate; σ = −(P_zz − (P_xx + P_yy)/2) and its bond / angle / dihedral / pair / kinetic parts (compute
//                      stress/atom, reduced); frames at even strain steps
//   analysis           modulus (linear fit up to `fit_strain`), yield (first maximum of the smoothed stress), strain
//                      softening (yield − the minimum after it), strain hardening modulus G_R from σ = σ₀ + G_R (λ² − 1/λ)
//                      after the softening minimum (Hoy & Robbins, J. Polym. Sci. B 44, 3487 (2006)); per frame ⟨P₂⟩ of the
//                      bonds along the pull, the end-to-end anisotropy ⟨R_z²⟩ / (½ (⟨R_x²⟩ + ⟨R_y²⟩)), the empty share of a
//                      grid (points farther than a probe radius from every bead) and Z through Z1+
// Dynamics
//   g₁ (the inner half of each chain), g₂ (relative to the chain's centre of mass), g₃ (the centres of mass), bond P₁(t) and
//   ⟨R_ee(t)·R_ee(0)⟩/⟨R_ee²⟩ over all time origins at logarithmic lags; D from g₃'s last decade, τ_R where the end-to-end
//   correlation falls to 1/e, τ_e where g₁'s local exponent first falls below 3/8 (between Rouse ½ and reptation ¼)
//   time mapping       the factor s with t_AA = s t_CG from the overlap of two g₁ curves (the geometric mean of the time
//                      ratios at equal g₁)
#pragma once
#include <string>
#include <vector>

#include "caps/cg_bonded.hpp"
#include "caps/system.hpp"

namespace caps {

// ---------------------------------------------------------------- entanglements
struct CgChainPaths {
  std::vector<double> lpp, r2;    // per chain: primitive path contour length (Å), squared end-to-end distance (Å²)
  std::vector<int> z;             // per chain: kinks (Z1+), empty when not known
  std::vector<int> beads;         // per chain
  int steps = 0;                  // CAPS's PPA: minimisation steps, whether it converged, the bead diameter used
  bool converged = true;
  double sigma = 0;
};
// Z1 format: the number of chains, the box edges, every chain's length on the third line, then x y z per bead (chains
// made whole along their bonds)
void write_z1_config(const CgTopology& t, const std::vector<Vec3>& pos, const Cell& cell, const std::string& path);
// Z1+'s shortest paths (Z1+SP.dat): node counts on the third line or before each chain, then the nodes; Z = nodes − 2
CgChainPaths read_z1_paths(const std::string& path);
// primitive paths from a PPA (CAPS's or LAMMPS's: the beads after the quench, ends where they were)
CgChainPaths paths_from_positions(const CgTopology& t, const std::vector<Vec3>& start, const std::vector<Vec3>& after, const Cell& cell);
// CAPS's PPA on one frame of beads
CgChainPaths caps_ppa(const CgTopology& t, const std::vector<Vec3>& pos, const Cell& cell, double sigma = 0, int max_steps = 200000);

struct CgEntanglement {
  double N = 0;                     // beads per chain (mean)
  double r2 = 0, lpp = 0, lpp2 = 0, z = -1;   // ⟨R²⟩, ⟨L_pp⟩, ⟨L_pp²⟩, ⟨Z⟩ (−1: not known)
  double ne_s_coil = 0, ne_s_kink = 0, ne_mod_s_kink = 0, ne_mod_s_coil = 0;
  double a_pp = 0;                  // tube step ⟨R²⟩ / ⟨L_pp⟩ (Å)
  int chains = 0;
};
CgEntanglement entanglement_of(const CgChainPaths& p);
struct CgMultiEstimate { double ne_m_kink = 0, ne_m_coil = 0; std::string note; };
// M-kink and M-coil from several chain lengths; `internal` the ⟨R²(n)⟩/n of the longest chains (n = 1, 2, …) and l0² the
// mean squared bond length, for C(x)
CgMultiEstimate multi_estimators(const std::vector<CgEntanglement>& sets, const std::vector<double>& internal, double l0sq);

// LAMMPS PPA deck: ends held (a group file written beside it), intra-chain pairs excluded, harmonic bonds with zero rest
// length, WCA beads of diameter sigma, a quench; the paths dumped for paths_from_positions
std::string ppa_deck(const CgTopology& t, double sigma, const std::string& ends_file);
std::string ppa_ends(const CgTopology& t);   // "group ends id …" lines

// ---------------------------------------------------------------- mechanics
struct CgTensionDeck {
  std::string mode = "stress";      // stress | volume
  double rate = 1e-7;               // engineering strain rate, 1/fs
  double max_strain = 3.0;
  double T = 300, P = 1, dt = 10;
  int frames = 60;                  // dumps at even strain steps
  int print_every = 0;              // steps between stress lines (0: 500 lines over the run)
  int exclude = 3;
};
std::string tension_deck(const CgTensionDeck& d);

struct CgStressStrain {
  std::vector<double> strain, stress, bond, angle, dihedral, pair, kinetic, density;   // MPa; g/cm³
};
CgStressStrain read_stress_strain(const std::string& path);

struct CgTensionAnalysis {
  double modulus = 0;               // MPa, linear fit up to fit_strain
  double yield_stress = 0, yield_strain = 0;
  double softening = 0;             // MPa, yield − the minimum after it
  double min_strain = 0;            // where that minimum is
  double hardening_modulus = 0, hardening_rms = 0, hardening_from = 0;   // G_R (MPa) and its fit from strain hardening_from
  std::string note;
};
CgTensionAnalysis analyse_tension(const CgStressStrain& c, double fit_strain = 0.02, double hardening_from = -1);

struct CgOrientation { double strain = 0, p2 = 0, ree_anisotropy = 0, void_fraction = -1; };
// one frame: ⟨P₂⟩ of the bonds along z, R_ee anisotropy, the empty share of a 1 Å grid at the probe radius (≤ 0: skip)
CgOrientation orientation_of(const CgTopology& t, const std::vector<Vec3>& pos, const Cell& cell, double probe = 0);

// ---------------------------------------------------------------- dynamics
struct CgDynamics {
  std::vector<double> t;            // lag, fs
  std::vector<double> g1, g2, g3;   // Å²
  std::vector<double> p1, ree;      // bond and end-to-end autocorrelations
  double D = 0;                     // Å²/fs, from g₃
  double tau_R = 0, tau_e = 0;      // fs (0: not reached)
  std::string note;
};
// frames: bead positions (unwrapped, or whole with minimum images per frame), times in fs
CgDynamics cg_dynamics(const CgTopology& t, const std::vector<std::vector<Vec3>>& frames, const std::vector<Cell>& cells, const std::vector<double>& times);
// s with t_AA = s t_CG (geometric mean over the overlap of the two g₁ curves); spread: the ratios' standard deviation in ln s
double time_mapping(const std::vector<double>& t_aa, const std::vector<double>& g_aa, const std::vector<double>& t_cg, const std::vector<double>& g_cg,
                    double* spread = nullptr, int* points = nullptr);

}  // namespace caps
