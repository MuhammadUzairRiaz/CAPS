// CAPS thermal conductivity by molecular dynamics: Green–Kubo (heat-flux autocorrelation, equilibrium NVE) and
// reverse non-equilibrium MD (Müller-Plathe, J. Chem. Phys. 106, 6082 (1997): kinetic energy moved from a cold slab to
// a hot slab by exchanging velocities, κ from the imposed flux and the temperature gradient that answers it).
//
// The heat flux (Irving–Kirkwood, as LAMMPS compute heat/flux with compute centroid/stress/atom):
//   Q = Σ_i e_i v_i + Σ_terms Σ_{k in term} (r_k − r_c) (F_k · v_k)
// e_i the atom's kinetic plus potential energy (pairs split half and half, a bonded term equally among its atoms),
// F_k the force the term puts on atom k and r_c the centroid of the term's atoms (positions taken whole about its first
// atom). For pairs this is ½ r_ij F_ij·(v_i + v_j); for angles, torsions and impropers the centroid form is the one that
// is right for many-body terms (Surblys et al., Phys. Rev. E 99, 051301 (2019); Boone, Babaei & Wilmer, J. Chem. Theory
// Comput. 15, 5579 (2019)): splitting their virial equally among the atoms, as compute stress/atom does, gives a wrong
// flux. The DSF self energy of each atom is part of its e_i (as LAMMPS's pe/atom); the Lennard-Jones tail correction,
// a constant of the whole cell, is not. On a single frame the flux agrees with LAMMPS compute heat/flux on compute
// centroid/stress/atom (argon and GAFF polyethylene with DSF, checked to 1e-6).
//
// Supported terms: 12-6 Lennard-Jones pairs (with the force field's mixing, overrides, 1-4 scaling or 1-4 parameters,
// shift when the tail correction is off), damped shifted force Coulomb (Fennell & Gezelter 2006, as CAPS's pair and
// exclusion terms), harmonic bonds and angles, Urey–Bradley, Fourier torsions and impropers, harmonic impropers. Other
// forms (class II cross terms, PME, Buckingham, many-body potentials, virtual sites) throw FieldError: write the LAMMPS
// inputs instead (thermal_lammps_input), where compute centroid/stress/atom does the same with LAMMPS's own terms.
//
// κ (W m⁻¹ K⁻¹): Green–Kubo κ_αα = 1/(V k_B T²) ∫₀^τ ⟨Q_α(0) Q_α(t)⟩ dt, the running integral's mean over the last
// third of τ (errors from blocks); NEMD κ = J / |dT/dz| with J = E_exchanged / (2 A t) (the flux leaves the hot slab in
// both directions of the periodic cell) and the gradient from straight-line fits to both halves of the temperature
// profile, the slabs next to the hot and cold slabs left out.
#pragma once
#include <array>
#include <cstdint>
#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "caps/dynamics.hpp"
#include "caps/field.hpp"
#include "caps/properties.hpp"
#include "caps/system.hpp"

namespace caps {

// (kcal/mol)/(Å fs K) in W m⁻¹ K⁻¹
constexpr double kKappaSI = 4184.0 / 6.02214076e23 / (1e-10 * 1e-15);

class HeatFlux {
 public:
  // Throws FieldError for terms it does not support (above).
  HeatFlux(const ForceField& ff, const EnergyOptions& e);
  // Q (kcal/mol·Å/fs) for positions x (3N, Å) and velocities v (3N, Å/fs); masses from the force field.
  Vec3 flux(const std::vector<double>& x, const std::vector<double>& v, const Cell& cell);
  // Per-atom potential energies (kcal/mol), centroid virials W_k,αβ = Σ (r_k − r_c)_α F_k,β (row-major 3×3, kcal/mol)
  // and the total forces (3N): Σ e_i is the potential energy without constant terms, Σ W_k the virial tensor.
  void per_atom(const std::vector<double>& x, const Cell& cell, std::vector<double>& e, std::vector<std::array<double, 9>>& w,
                std::vector<double>* forces = nullptr);
  const std::vector<double>& masses() const { return mass_; }
  // the constant energy left out of the per-atom sum (the Lennard-Jones tail correction; 0 without it)
  double constant_energy() const { return econst_; }

 private:
  const ForceField& ff_;
  EnergyOptions e_;
  std::vector<double> mass_, q_, self_;
  size_t nt_ = 0;
  std::vector<double> eps_, s6_, eps14_, s614_;
  double rc_ = 10, econst_ = 0;
  bool coul_ = false, shift_ = false;
};

struct ThermalOptions {
  std::shared_ptr<const ForceField> field;
  EnergyOptions energy;
  std::string method = "gk";        // gk | nemd
  // 3: a bulk material (κ along x, y, z; the cell's volume). 2: a sheet or film periodic in x and y with vacuum along z
  // (never guessed from the structure: the caller says which): only the in-plane flux, κ over the volume A·thickness
  // (thickness 0: the atoms' z extent plus their van der Waals radii, said in the notes), the sheet conductance κ·t too
  int dimensionality = 3;
  double thickness = 0;             // Å, 2D only
  double temperature = 300.0;       // K
  double equilibrate_ps = 10.0;     // NVT (Bussi) before the sampled run
  double ps = 100.0;                // the sampled NVE run (gk) or the NEMD run
  double dt = 1.0;                  // fs
  double tau_t = 100.0;             // fs, equilibration thermostat
  // Green–Kubo
  int sample_every = 5;             // steps between heat-flux samples
  double corr_ps = 5.0;             // correlation time integrated
  int blocks = 5;
  // NEMD (Müller-Plathe)
  int axis = 2;                     // 0 x, 1 y, 2 z: the flux direction (gk: reported for all three)
  int slabs = 20;                   // even; slab 0 cold, slab slabs/2 hot
  int swap_every = 100;             // steps between exchanges
  double transient = 0.3;           // fraction of the NEMD run before the profile and flux are averaged
  uint64_t seed = 1;
  bool new_velocities = false;
  std::function<bool(const std::string& what, double fraction)> progress;
};

struct ThermalResult {
  std::string method;
  std::array<double, 3> kappa{0, 0, 0}, error{0, 0, 0};   // W/(m K); NEMD: kappa[axis]
  double kappa_mean = 0, kappa_mean_error = 0;           // GK: the mean of the three
  int axis = -1;
  int dimensionality = 3;
  double thickness = 0;              // Å, the effective thickness κ is normalised by (2D)
  double sheet_conductance = 0;      // 2D: κ·t, W/K (the sheet's thermal conductance per unit width and length)
  // Green–Kubo: running integrals and normalised autocorrelations per axis
  std::vector<double> t_ps;
  std::array<std::vector<double>, 3> running, acf;
  bool plateau = true;
  // NEMD: slab centres (Å along the axis from the cold slab), temperatures, flux and gradient
  std::vector<double> slab_x, slab_t;
  double flux = 0, gradient = 0;     // W/m², K/m
  double exchanged = 0;              // kcal/mol moved from cold to hot over the averaged part
  int exchanges = 0;
  double temperature = 0;            // mean temperature of the sampled part, K
  std::vector<std::string> notes;
};

// The Green–Kubo integral of heat-flux samples Q (kcal/mol·Å/fs), dt_fs apart, in a cell of volume V (Å³) at T (K).
ThermalResult green_kubo_kappa(const std::vector<Vec3>& q, double dt_fs, double volume, double temperature, double corr_ps, int blocks);
// A Müller-Plathe temperature profile and exchanged energy to κ: slab temperatures (slab 0 cold, n/2 hot), slab
// width (Å), the exchanged energy (kcal/mol) over time (fs), cross-section area (Å²).
ThermalResult nemd_kappa(const std::vector<double>& slab_t, double slab_width, double exchanged, double time_fs, double area);

// The effective thickness of a sheet: the atoms' extent along z plus the van der Waals radii (Bondi) of the lowest and
// highest atoms (Å).
double sheet_thickness(const System& s);

// Runs the method on s (NVT equilibration, then NVE): s ends as the last configuration.
ThermalResult thermal_conductivity(System& s, const ThermalOptions& o);
std::vector<Property> thermal_properties(const ThermalResult& r);
std::string thermal_summary(const ThermalResult& r);

// LAMMPS input sections (appended to a CAPS input set up with read_data and the force field): Green–Kubo with compute
// heat/flux on compute centroid/stress/atom (the many-body-correct form) and fix ave/correlate, or Müller-Plathe NEMD with
// fix thermal/conductivity and the slab temperature profile (chunk/atom). κ is written by the input itself.
struct ThermalLammps {
  std::string method = "gk";
  double temperature = 300, dt = 1.0;
  double equilibrate_ps = 100, ps = 2000;
  int sample_every = 5;
  double corr_ps = 10;
  int axis = 2, slabs = 20, swap_every = 100;
  int dimensionality = 3;
  double thickness = 0;   // Å (2D: the GK volume is lx·ly·thickness; NEMD's cross-section ly·thickness or lx·thickness)
  bool kspace = false;   // the force field runs with PPPM: centroid/stress/atom then needs a pair-only flux (a note in the input)
};
std::string thermal_lammps_run(const ThermalLammps& q);

}  // namespace caps
