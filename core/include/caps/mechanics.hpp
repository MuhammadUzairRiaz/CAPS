// CAPS mechanics and thermal transitions: elastic constants (static strain with minimisation, and stress fluctuations
// over an NVT trajectory), uniaxial deformation (stress–strain) and the glass transition from a stepwise cooling run.
// Units: stresses and moduli in GPa (stress–strain curve in MPa), strains dimensionless, temperatures in K.
// Voigt order 1 xx, 2 yy, 3 zz, 4 yz, 5 xz, 6 xy; shear strains are engineering strains (ε4 = 2ε_yz).
#pragma once
#include <array>
#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "caps/dynamics.hpp"
#include "caps/field.hpp"
#include "caps/properties.hpp"
#include "caps/system.hpp"

namespace caps {

using Mat6 = std::array<std::array<double, 6>, 6>;

struct ElasticResult {
  Mat6 C{};                      // GPa
  Mat6 err{};                    // standard error (over configurations or blocks); NaN when not estimated
  Mat6 born{}, fluct{}, kinetic{};   // fluctuation method: C = ⟨C_Born⟩ − fluctuation + kinetic (GPa)
  std::array<double, 6> prestress{};   // stress of the unstrained state (σ = −P), GPa
  // isotropic averages (GPa, ν dimensionless): Voigt (uniform strain), Reuss (uniform stress), Hill (their mean)
  double K_voigt = 0, G_voigt = 0, K_reuss = 0, G_reuss = 0, K_hill = 0, G_hill = 0, E_hill = 0, nu_hill = 0, lambda_hill = 0;
  double asymmetry = 0;          // max |C_IJ − C_JI| / max |C|, before symmetrising (static method)
  int configurations = 0;
  std::string method;
  std::vector<std::string> notes;
};

// Voigt, Reuss and Hill averages of C (sets the isotropic fields).
void isotropic_averages(ElasticResult& r);

// Static elastic constants (Theodorou and Suter, Macromolecules 19, 139 (1986)): each configuration is minimised at
// fixed cell; each of the six strains is applied as ±strain (pure, symmetric, affine), the atoms are re-minimised and
// C_IJ = Δσ_I / (2 strain). Averaged over the configurations.
struct StaticElasticOptions {
  std::shared_ptr<const ForceField> field;   // null: built-in GAFF of C and H
  EnergyOptions energy;
  double strain = 1e-4;
  double ftol = 1e-4;            // kcal/mol/Å, minimiser force tolerance
  int max_iterations = 50000;
  std::function<bool(const std::string& what, double fraction)> progress;   // return false to cancel
};
ElasticResult static_elastic(const std::vector<System>& configs, const StaticElasticOptions& o);

// Elastic constants from stress fluctuations at constant volume and temperature (Lutsko, J. Appl. Phys. 65, 2991
// (1989); Clavier et al., Mol. Sim. 43, 1413 (2017)):
//   C_IJ = ⟨C^B_IJ⟩ − (V / kT) cov(σ_I, σ_J) + N kT / V (δ_IJ + [I, J ≤ 3])
// with σ the virial stress. The Born term C^B = (1/V) ∂²U/∂η_I∂η_J (Lagrangian strain) is the central difference of
// the second Piola–Kirchhoff virial stress S = −F⁻¹ W F⁻ᵀ under F = √(I + 2η), with the pair set frozen (no impulsive
// cut-off term, as analytic Born terms). Frames must come from an NVT run at `temperature`.
struct FluctuationOptions {
  const ForceField* ff = nullptr;   // required: the force field of the trajectory's atoms
  EnergyOptions energy;
  double temperature = 300.0;
  double strain = 1e-5;          // finite-difference step of the Born term
  int blocks = 5;
  std::function<bool(const std::string& what, double fraction)> progress;
};
ElasticResult fluctuation_elastic(const Trajectory& t, const std::vector<size_t>& frames, const FluctuationOptions& o);

// The same on the fly: NVT dynamics of `s` at `temperature` for `ps`, the virial stress sampled at every step (the
// covariance needs very many samples: for bonded polymers the fluctuation term cancels ~99 % of the Born term) and the
// Born term every `born_every` steps; block errors over time. `s` ends as the last configuration.
struct FluctuationRunOptions {
  std::shared_ptr<const ForceField> field;
  EnergyOptions energy;
  double temperature = 300.0, ps = 100.0, equilibrate_ps = 10.0, dt = 1.0, tau_t = 1000.0;
  Thermostat thermostat = Thermostat::Langevin;
  int born_every = 100;
  double strain = 1e-5;
  int blocks = 5;
  uint64_t seed = 1;
  bool new_velocities = false;
  std::function<bool(const std::string& what, double fraction)> progress;
};
ElasticResult fluctuation_run(System& s, const FluctuationRunOptions& o);

// Born matrix of one configuration (GPa), as fluctuation_elastic computes it (tests and benches).
Mat6 born_matrix(Evaluator& ev, const std::vector<double>& x, const Cell& cell, double strain);
// Virial stress σ = −W / V of one configuration, Voigt, kcal/mol/Å³.
std::array<double, 6> virial_stress(const EnergyTerms& e, double volume);

// Uniaxial deformation at constant engineering strain rate (MD, as LAMMPS fix deform erate), the lateral axes held at
// the target pressure by per-axis Berendsen coupling or kept fixed. Stress σ = −P_aa (true stress).
struct TensilePoint {
  double strain = 0, stress = 0;       // MPa
  double lateral1 = 0, lateral2 = 0;   // strains of the two other axes
  double temperature = 0, time_ps = 0;
};
struct TensileOptions {
  std::shared_ptr<const ForceField> field;
  EnergyOptions energy;
  int axis = 0;                  // 0 x, 1 y, 2 z
  double rate = 1e-3;            // engineering strain rate, 1/ps (1e-3/ps = 10⁹ s⁻¹)
  double max_strain = 0.2;
  double temperature = 300.0, dt = 1.0, tau_t = 100.0;
  Thermostat thermostat = Thermostat::Bussi;
  bool lateral_pressure = true;  // false: lateral axes fixed (uniaxial strain)
  double pressure = 1.0, tau_p = 1000.0, compressibility = 4.5e-5;
  double equilibrate_ps = 20;    // NPT first (each axis at `pressure`, Berendsen per axis), so the pull starts stress-free
  int sample_every = 50;         // steps between stress samples
  double fit_strain = 0.02;      // modulus and Poisson ratio from 0 to this strain
  uint64_t seed = 1;
  bool new_velocities = false;
  std::function<bool(const TensilePoint&)> progress;   // return false to cancel
  std::function<void(const std::vector<double>&, const Cell&, int64_t)> frame;
  int frame_every = 0;
};
struct TensileResult {
  std::vector<TensilePoint> curve;
  std::vector<double> smooth;    // stress smoothed over ±0.5 % strain, MPa
  double modulus = 0, modulus_err = 0;   // GPa
  double poisson = 0, poisson_err = 0;
  double yield_stress = 0, yield_strain = 0;   // 0.2 % offset on the smoothed curve (0 when not reached)
  double peak_stress = 0, peak_strain = 0;
  std::string method;
  std::vector<std::string> notes;
};
TensileResult run_tensile(System& s, const TensileOptions& o);
// The analysis part (modulus, Poisson ratio, yield, peak) of a curve.
void analyse_tensile(TensileResult& r, double fit_strain, bool lateral);

// Glass transition: stepwise cooling at constant pressure; at each temperature the density is averaged over the second
// part of the hold, and Tg is the hinge of a continuous two-line fit of specific volume against T.
struct BilinearFit {
  bool ok = false;
  double tg = 0, tg_err = 0;     // K (error: bootstrap of residuals)
  double slope_low = 0, slope_high = 0, value_at_tg = 0;   // dv/dT below and above, v(Tg)
  double alpha_low = 0, alpha_high = 0;   // volumetric expansion (1/v) dv/dT at Tg, 1/K
  double rss = 0;
  std::string note;
};
// y(T) = y0 + b_low min(T − Tg, 0) + b_high max(T − Tg, 0), least squares with the hinge free.
BilinearFit fit_bilinear(const std::vector<double>& T, const std::vector<double>& y);

struct CoolingPoint {
  double temperature = 0, density = 0, density_err = 0, specific_volume = 0, potential = 0;
};
struct CoolingOptions {
  std::shared_ptr<const ForceField> field;
  EnergyOptions energy;
  double t_start = 500, t_end = 200, t_step = 20;   // K
  double ps_per_step = 100;      // hold at each temperature
  double equilibrate_ps = -1;    // unsampled run at t_start first (the structure may come from another temperature); −1: one hold
  double average_from = 0.5;     // fraction of each hold discarded before averaging
  double dt = 1.0, tau_t = 100.0, tau_p = 1000.0, pressure = 1.0, compressibility = 4.5e-5;
  Barostat barostat = Barostat::CRescale;
  uint64_t seed = 1;
  bool new_velocities = false;
  std::function<bool(const ThermoRow&, int step_index, int steps)> progress;   // return false to cancel
};
struct CoolingResult {
  std::vector<CoolingPoint> points;
  BilinearFit fit;
  std::vector<ThermoRow> thermo;
  std::string method;
  std::vector<std::string> notes;
};
CoolingResult run_cooling(System& s, const CoolingOptions& o);

// Results as Analyze properties (cards with method, errors, notes and curves).
//   elastic: cij (value: mean of C11, C22, C33; every C_IJ in extra), youngs, bulk, shear, poisson (Hill averages);
//            ids get the suffix given ("" for static strain, "_fluct" for fluctuations)
//   tensile: tensile_modulus (with Poisson ratio, yield and peak in extra; stress–strain curves), yield
//   cooling: tg (with expansion coefficients; specific volume and density against T, and the two-line fit)
std::vector<Property> elastic_properties(const ElasticResult& r, const std::string& suffix);
std::vector<Property> tensile_properties(const TensileResult& r);
std::vector<Property> cooling_properties(const CoolingResult& r);

}  // namespace caps
