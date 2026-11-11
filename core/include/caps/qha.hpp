#pragma once
// Quasi-harmonic free energy and thermal expansion of a crystal: the harmonic phonons (phonons.hpp) at a set of strained
// cells, each with its atoms relaxed, give G(strain, T) = E_static + F_vib + P V; a polynomial fitted to G at each
// temperature has its minimum at the cell of that temperature (the numerical route, as in Togo & Tanaka, Scr. Mater. 108,
// 1 (2015); free-energy minimisation of crystals: Gale, J. Phys. Chem. B 102, 5423 (1998), there with analytical
// derivatives).
//   volume  the cell scaled as a whole: G as a polynomial (degree up to 4) in V/V₀ − 1; for a cubic or isotropic crystal
//   axes    a, b and c scaled on their own, angles kept: G as a polynomial in the free strains with every term to third
//           order, x_i⁴ and x_i²x_j², from each axis at ±δ, ±2δ (δ = range / 2), pairs of axes at (±δ, ±δ) and all three at
//           ±δ, ±2δ: 29 cells for three axes; for anisotropic crystals (a polymer's negative expansion along its chains)
// At the minimum: α_V and α_a,b,c by central differences over ±ΔT, B_T = V / (vᵀ H⁻¹ v) with H the Hessian of G in the
// strains and v = ∂V/∂strain, C_v and S from fits of their own, C_p = C_v + T V B_T α_V².
// Charges are summed by Ewald (PME) whatever the options say (the damped shifted force cut-off puts kinks in G(V)), and
// the Lennard-Jones term is shifted at the cut-off with no tail term whatever the options say: the energies of a few per
// cent strain are small, and a truncated potential jumps each time a pair crosses the cut-off (a tenth of the signal in
// polyethylene). Limits: no anharmonicity beyond the volume dependence of the modes, so results degrade towards melting.
// The volume mode strains a, b and c alike, so for an anisotropic crystal (a polymer chain crystal such as PE) its B_T is a
// Voigt-like average over stiff and soft axes; the axes mode is the right choice there.
// A cell with imaginary mesh modes is left out of every fit (its F_vib lacks those modes and jumps against its neighbours),
// and the stable cells set the range: a minimum beyond them is marked extrapolated. When too few stable cells remain, all
// cells are used and every row is marked unstable. Every fit is least squares with more cells than terms (volume: degree
// 2 from 4 cells, 3 from 6, 4 from 7; axes: the highest-order terms go first when the cells cannot carry them).
// --points (QhaOptions::points) spaces the volume grid only; the axes grid is the fixed design above.
#include <array>
#include <functional>
#include <string>
#include <vector>

#include "caps/phonons.hpp"
#include "caps/system.hpp"

namespace caps {

struct QhaOptions {
  PhononOptions phonon;                       // field_for, energy, step, supercell, mesh, relax (its temperatures are set here)
  std::string mode = "volume";                // volume | axes
  std::array<bool, 3> axes{true, true, true}; // axes: which of a, b, c may change
  double range = 0.02;                        // linear strain at the ends of the grid
  int points = 7;                             // volume: cells on the grid (at least 4)
  std::vector<double> temperatures{0, 100, 200, 298.15, 400};   // K
  double pressure = 0;                        // atm
  bool relax_cell = true;                     // the static cell relaxed first (the grid centred on it)
  std::function<bool(const std::string&, double)> progress;   // false stops
};

// One strained cell. F, Cv, S per requested temperature and j = 0, 1, 2 (T − ΔT, T, T + ΔT, see qha_step).
struct QhaSample {
  std::array<double, 3> strain{0, 0, 0};      // linear strains of a, b, c
  double volume = 0, energy = 0;              // Å³, kcal/mol per cell (static)
  std::vector<std::array<double, 3>> F, Cv, S;   // kcal/mol, cal/mol/K, cal/mol/K per cell (vibrational)
  int imaginary = 0;
  double min_wavenumber = 0;
};
struct QhaRow {
  double T = 0;
  std::array<double, 3> strain{0, 0, 0};
  std::array<double, 3> length{0, 0, 0};      // a, b, c, Å
  double volume = 0, G = 0;                   // Å³, kcal/mol per cell
  double alpha_V = 0;                         // 1/K
  std::array<double, 3> alpha{0, 0, 0};       // along a, b, c, 1/K
  double B = 0;                               // isothermal bulk modulus, GPa
  double Cv = 0, Cp = 0, S = 0;               // cal/mol/K per cell
  bool extrapolated = false;                  // the minimum (or one at T ± ΔT) lies outside the stable strains sampled
  bool unstable = false;                      // G is not a minimum there, or the fit had to use cells with imaginary modes
};
struct QhaResult {
  std::string mode;
  std::vector<QhaRow> rows;
  std::vector<QhaSample> samples;
  std::array<double, 3> lengths0{0, 0, 0};    // the centre cell's a, b, c
  double volume0 = 0;
  double fit_rms = 0;                         // kcal/mol, the largest over the temperatures
  std::string field;
  std::vector<std::string> notes;
};

// The ΔT used about a requested temperature (1 K, less near zero); at T = 0 the three are 0, ΔT, 2ΔT and α = 0 (third law).
double qha_step(double T);
// The fit alone. variable[k]: strain k free (mode "axes"); mode "volume" treats the cell as one variable.
QhaResult qha_analyze(const std::vector<QhaSample>& samples, const std::vector<double>& temperatures, const std::string& mode,
                      const std::array<bool, 3>& variable, const std::array<double, 3>& lengths0, double volume0, double pressure_atm);
QhaResult quasi_harmonic(const System& cell, const QhaOptions& o);

}  // namespace caps
