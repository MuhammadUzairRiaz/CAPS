// CAPS polymer statistics (design/boards/Copolymer, Tacticity, BlendPhase, SolventScreen, Electrostatics): the small
// analytic models the builders and settings pages show beside what was built.
//
//  copolymer   the terminal (Mayo–Lewis) model: instantaneous composition, azeotrope, mean run lengths
//  stereo      Bernoulli and first-order Markov dyad models: triads, the ten pentads, fits to measured pentads
//  blend       Flory–Huggins free energy of mixing for a binary blend: critical point, spinodal, binodal
//  solvent     the Hildebrand estimate χ ≈ V (δs − δp)² / RT + 0.34
//  ewald       the splitting parameter from a tolerance and FFT-friendly mesh sizes
#pragma once
#include <array>
#include <string>
#include <vector>

namespace caps {

// ---- copolymer (Mayo & Lewis, J. Am. Chem. Soc. 66, 1594 (1944)); low conversion, no composition drift
struct CopolymerModel {
  double r1 = 1, r2 = 1, f1 = 0.5;
  double F1 = 0.5;                 // A in the copolymer formed from feed f1
  double paa = 0.5, pbb = 0.5;     // P(A follows A), P(B follows B)
  double run_a = 2, run_b = 2;     // mean run lengths 1 / (1 − p)
  double azeotrope = -1;           // feed with F1 = f1, or −1 when there is none in (0, 1)
};
double mayo_lewis(double r1, double r2, double f1);
CopolymerModel copolymer_terminal(double r1, double r2, double f1);

// ---- stereo sequence statistics (Bovey; Randall). Dyads are 'm' or 'r'.
constexpr int kPentadCount = 10;
// mmmm mmmr rmmr mmrr mmrm rmrm rmrr rrrr rrrm mrrm (reversed sequences counted with their pair)
const std::array<const char*, kPentadCount>& pentad_names();
struct StereoModel {
  std::string kind = "bernoulli";  // bernoulli | markov
  double pm = 0.5;                 // meso fraction
  double p_mr = 0.5, p_rm = 0.5;   // Markov: P(r after m), P(m after r); Bernoulli: 1 − pm and pm
  double mm = 0.25, mr = 0.5, rr = 0.25;
  std::array<double, kPentadCount> pentads{};
};
StereoModel stereo_bernoulli(double pm);
StereoModel stereo_markov(double p_mr, double p_rm);
struct StereoCounts {
  int m = 0, r = 0, mm = 0, mr = 0, rr = 0, pentad_total = 0;
  std::array<int, kPentadCount> pentads{};
};
StereoCounts count_stereo(const std::string& dyads);
// Least-squares fits to measured pentad fractions (normalised first); rms is the fit's residual.
StereoModel fit_bernoulli(const std::array<double, kPentadCount>& measured, double* rms = nullptr);
StereoModel fit_markov(const std::array<double, kPentadCount>& measured, double* rms = nullptr);
// A dyad string drawn from the model (dp units → dp − 1 dyads).
std::string draw_dyads(const StereoModel& m, int dp, uint64_t seed);

// ---- Flory–Huggins blend (Flory 1942; Huggins 1942): f/kT per site = φ/NA ln φ + (1−φ)/NB ln(1−φ) + χ φ(1−φ), φ of A
struct BlendCritical { double chi_c = 0, phi_c = 0.5; };
BlendCritical blend_critical(double na, double nb);
// Spinodal compositions at χ (∂²f/∂φ² = 0); false below χc.
bool blend_spinodal(double na, double nb, double chi, double& lo, double& hi);
// Coexisting compositions at χ (equal exchange chemical potential and equal grand potential); false below χc or when
// Newton fails. guess_lo / guess_hi: a previous solution to continue from (0 for none).
bool blend_binodal(double na, double nb, double chi, double& lo, double& hi, double guess_lo = 0, double guess_hi = 0);

// ---- Hildebrand χ with the entropic 0.34 (V in cm³/mol, δ in MPa½, T in K)
double hildebrand_chi(double v_solvent, double delta_solvent, double delta_polymer, double t);

// ---- Ewald
// Mesh points along an edge: the smallest size ≥ edge / spacing (and ≥ order + 1) whose factors are 2, 3, 5 and 7.
int pme_mesh_size(double edge, double spacing, int order = 4);

}  // namespace caps
