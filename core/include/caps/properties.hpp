// CAPS Analyze: properties from a trajectory (or one structure) — structure, chains, dynamics, cohesive energy and free
// volume. Every property carries its method, its estimate of uncertainty (block averages over frames where it applies)
// and the curves behind it.
//
// Properties (ids):
//   density      mass / volume over the frames                                  g/cm³
//   rdf          g(r) between element sets, averaged over frames; first peak      Å
//   sq           total structure factor S(q) (Faber–Ziman): direct reciprocal-lattice sum at low q, g(r) transform above; first peak   Å⁻¹
//   xray         X-ray I(q) (Cromer–Mann form factors of every element to Cf, International Tables Vol. C)
//   electron     electron-diffraction I(q) (Peng et al. 1996 elastic scattering factors, International Tables Vol. C)
//   neutron      neutron S(q) (coherent scattering lengths)
//   rg           radius of gyration of molecules, √⟨Rg²⟩                          Å
//   ree          backbone end-to-end distance, √⟨R²⟩                            Å
//   cn           characteristic ratio C_n = ⟨R²(n)⟩ / (n ⟨b²⟩) and C∞ (extrapolated in 1/n)
//   persistence  persistence length: ⟨R·b₁⟩/|b| (Flory) and the exponential decay of bond correlations   Å
//   msd          mean-square displacement of atoms and molecule centres (all time origins, drift removed)
//   diffusion    D from the linear part of the molecule-centre MSD (Einstein), 10⁻⁵ cm²/s
//   relaxation   end-to-end vector and backbone bond (P2) autocorrelations, KWW fits; τ   ps
//   ced          cohesive energy density (E_isolated − E_bulk)/V with the force field, and δ = √CED   J/cm³, MPa^½
//   ffv          free volume by probe insertion: accessible fraction and Bondi FFV
//   entanglements  primitive-path analysis (Everaers 2004): N_e (modified S-coil), M_e, tube step, plateau modulus   bonds
//   psd          pore size distribution (largest sphere containing each free point, Gelb & Gubbins)   Å
//   cij_fluct    elastic constants from stress fluctuations of an NVT trajectory (needs ff and temperature); adds
//                youngs_fluct, bulk_fluct, shear_fluct, poisson_fluct (see mechanics.hpp)   GPa
#pragma once
#include <algorithm>
#include <functional>
#include <limits>
#include <map>
#include <string>
#include <vector>

#include "caps/field.hpp"
#include "caps/system.hpp"

namespace caps {

struct Series {
  std::string label, x_label, y_label;
  std::vector<double> x, y;
};

struct Property {
  std::string id, name, unit, method;
  double value = std::numeric_limits<double>::quiet_NaN();
  double error = std::numeric_limits<double>::quiet_NaN();   // standard error (block averages), NaN when not estimated
  std::map<std::string, double> extra;                        // further numbers (named with their unit)
  std::vector<Series> series;
  std::vector<std::string> notes;
};

struct AnalyzeOptions {
  // frames: first, last (inclusive; −1 = the last), stride
  long first = 0, last = -1, stride = 1;
  double frame_ps = 0;          // time between stored frames; 0: from the timesteps × timestep_fs
  double timestep_fs = 1.0;     // MD time step of the timesteps in the file
  int blocks = 5;               // block averages for the standard errors
  // structure
  int elem_a = 0, elem_b = 0;   // g(r) element pair (0 = any)
  bool inter_only = false;
  double rdf_rmax = 0;          // 0: half the shortest cell width (at most 15 Å)
  double rdf_dr = 0.02;
  double qmax = 25.0, dq = 0.02;   // Å⁻¹
  double q_direct = 4.0;           // S(q) by the direct reciprocal-lattice sum up to here (Å⁻¹); 0: from g(r) only
  // neutron contrast: which hydrogens scatter as deuterium (b = 6.671 fm) — 0 none, 1 every H, 2 H on aliphatic
  // carbons (a deuterated backbone), 3 H on aromatic carbons (deuterated rings), 4 H on O and N (exchanged in D₂O)
  int deuterate = 0;
  // dynamics
  double fit_from = 0.2, fit_to = 0.5;   // Einstein fit window, as fractions of the run
  // cohesive energy: the force field of the structure (atoms in trajectory order), and its energy settings
  const ForceField* ff = nullptr;
  EnergyOptions energy;
  // elastic constants from stress fluctuations (cij_fluct): the temperature of the NVT run the frames come from
  double temperature = 0;       // K
  // response functions from fluctuations (fluct): the pressure of the NPT run the frames come from
  double pressure = 1.0;        // atm
  // interfaces: bin width of the density profile along z (Å); a molecule left out of the chain analyses (the substrate)
  double zbin = 0.5;
  int64_t exclude_mol = 0;
  // interfaces and orientation: the direction profiles, adhesion and Herman's f run along — the normal of the cell face
  // spanned by the other two axes (0 = a: x in an orthogonal cell, 1 = b: y, 2 = c: z)
  int axis = 2;
  // the surface or filler of the interface properties (zprofile, adhesion, interaction): these molecule ids; empty: molecule 1
  std::vector<int64_t> surface_mols;
  bool is_surface(int64_t mol) const {
    return surface_mols.empty() ? mol == 1 : std::find(surface_mols.begin(), surface_mols.end(), mol) != surface_mols.end();
  }
  int ppa_frames = 3;           // entanglements: primitive paths of this many frames, spread over the chosen ones
  // free volume
  double probe = 0.0;           // probe radius, Å (0: points outside every van der Waals sphere)
  std::string radii = "bondi";  // atom radii: bondi (Bondi 1964), uff (x/2, Rappé 1992), forcefield (half the assigned LJ minimum)
  double grid = 0.4;            // Å
  int threads = 0;
  std::function<bool(const std::string& what, double fraction)> progress;   // return false to cancel
};

// The frames an analysis uses.
std::vector<size_t> analysis_frames(const Trajectory& t, const AnalyzeOptions& o);
// Time of each stored frame, ps.
std::vector<double> frame_times(const Trajectory& t, const AnalyzeOptions& o);

// Runs the requested properties (ids above). Unknown ids throw std::invalid_argument; a property that cannot be
// computed for this input (no chains, no force field, one frame for dynamics) comes back with a note and no value.
std::vector<Property> analyze(const Trajectory& t, const std::vector<std::string>& ids, const AnalyzeOptions& o);

// Scattering data: neutron coherent scattering length (fm) and X-ray form factor f(q) (electrons; q in Å⁻¹).
double neutron_b(int z);                // z = kDeuterium: ²H
constexpr int kDeuterium = 1001;
// The hydrogens a deuteration choice (AnalyzeOptions::deuterate) marks, by atom.
std::vector<char> deuterated_hydrogens(const System& s, int pattern);
double xray_f(int z, double q);
double electron_f(int z, double q);   // Å; Peng et al. 1996
// Atom radii for free volume and voids: "bondi" (Bondi 1964), "uff" (x/2, Rappé 1992), "forcefield" (half the Lennard-Jones
// minimum of ff, which must be assigned to s). Throws std::invalid_argument otherwise.
std::vector<double> free_volume_radii(const System& s, const std::string& kind, const ForceField* ff = nullptr);
// The scattering weight of element z at q for kind xray | electron | neutron.
double scatter_w(const std::string& kind, int z, double q);

// Results as JSON: [{id, name, value, error, unit, method, extra{...}, notes[...], series[{label, x_label, y_label, x[], y[]}]}]
// (NaN as null).
std::string properties_json(const std::vector<Property>& props);

// Part of a force field: the atoms given (in that order) and the terms among them only.
ForceField subset_forcefield(const ForceField& ff, const std::vector<uint32_t>& atoms);

}  // namespace caps
