// CAPS whole-pattern refinement: a measured powder pattern fitted by the crystal phases' calculated patterns.
//
//  Pawley      each reflection group's intensity a free (non-negative) parameter, with the lattice, zero shift and peak
//              widths (Pawley, J. Appl. Cryst. 14, 357 (1981)): the cell and profile without a structure.
//  Rietveld    the intensities from each phase's structure, y_c(2θ) = b(2θ) + Σ_p S_p Σ_hkl m |F|² LP P_MD Ω(2θ − 2θ_hkl)
//              (Rietveld, J. Appl. Cryst. 2, 65 (1969)), F with Cromer–Mann form factors (or neutron lengths) and a
//              per-phase isotropic B; Ω the Thompson–Cox–Hastings pseudo-Voigt of powder.hpp, LP and March–Dollase as
//              there.
//  Least squares  S = Σ w_i (y_i − y_c,i)², w_i = 1/σ_i² (σ_i² = y_i when not given). The linear parameters (scales,
//              Pawley intensities, the Chebyshev background) are solved exactly at every trial of the others (zero shift,
//              lattice with the cell's equalities kept, U V W X Y, B, March r) by Levenberg–Marquardt.
//  Results     Rp, Rwp, Rexp, χ² = (Rwp/Rexp)²; the lattice; Scherrer size from Y and Lorentzian strain from X; weight
//              fractions W_p = S_p (ZMV)_p / Σ_i S_i (ZMV)_i (Hill & Howard, J. Appl. Cryst. 20, 467 (1987)), ZM the
//              cell's mass and V its volume.
// Atom positions, occupancies and rigid-body or torsion parameters are not refined.
#pragma once

#include <array>
#include <functional>
#include <string>
#include <vector>

#include "caps/json.hpp"
#include "caps/system.hpp"

namespace caps {

struct MeasuredPattern {
  std::vector<double> two_theta, y, sigma;   // degrees; counts; sigma empty: √y
};
// Two or three columns (2θ, intensity[, σ]); lines starting with # or ! and non-numeric lines are skipped.
MeasuredPattern read_pattern(const std::string& path);

struct RefinePhase {
  std::string name;
  System cell;                         // the structure (its cell and atoms)
  // which nonlinear parameters to refine for this phase
  bool lattice = true, b_iso = false, march = false;
  double b = 0;                        // Å², start
  std::array<int, 3> preferred{0, 0, 0};
  double march_r = 1;
};

struct RefineOptions {
  std::string mode = "rietveld";       // rietveld | pawley (one phase)
  std::string radiation = "xray";
  std::string anode = "Cu";
  double wavelength = 0;
  bool kalpha2 = true;
  double ka2_ratio = 0.5;
  double monochromator_2theta = 0;
  bool polarised = false;
  int background = 6;                  // Chebyshev terms in 2θ
  // profile starts and which to refine
  double U = 0, V = 0, W = 0.01, X = 0, Y = 0.02, zero = 0;
  bool refine_U = true, refine_V = false, refine_W = true, refine_X = false, refine_Y = true, refine_zero = true;
  double two_theta_min = 0, two_theta_max = 0;   // the range fitted (0: the data's)
  int iterations = 50;
  std::function<bool(const std::string&, double)> progress;
};

struct RefinedPhase {
  std::string name;
  std::array<double, 6> lattice{};     // a b c (Å) alpha beta gamma (degrees)
  std::array<double, 6> lattice_esd{};
  double scale = 0, b = 0, march_r = 1;
  double weight_fraction = 0;          // Hill & Howard
  double volume = 0, mass = 0;         // cell, Å³ and g/mol
  std::vector<std::array<int, 3>> hkl; // reflection groups (one member), with
  std::vector<double> two_theta, intensity;   // their positions (Kα1) and integrated intensities (Pawley: refined)
};

struct RefineResult {
  std::vector<double> two_theta, y_obs, y_calc, background;   // the fitted range
  std::vector<RefinedPhase> phases;
  double U = 0, V = 0, W = 0, X = 0, Y = 0, zero = 0;
  double size = 0, strain = 0;         // Å from Y (Scherrer K 0.9), Lorentzian microstrain from X
  double rp = 0, rwp = 0, rexp = 0, chi2 = 0;
  int iterations = 0, parameters = 0;
  std::vector<std::pair<std::string, double>> refined;   // name, value of every refined nonlinear parameter
  std::vector<std::string> notes;
};

RefineResult refine_pattern(const MeasuredPattern& data, const std::vector<RefinePhase>& phases, const RefineOptions& o);

// The calculated pattern alone (no fitting): every phase at the given scale (phases[k].b, preferred, march_r used), for
// making test data and for checks.
std::vector<double> calculated_pattern(const std::vector<double>& two_theta, const std::vector<RefinePhase>& phases, const std::vector<double>& scales,
                                       const RefineOptions& o);

Json refine_json(const RefineResult& r);

}  // namespace caps
