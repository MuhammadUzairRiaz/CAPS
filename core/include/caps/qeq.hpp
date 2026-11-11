// CAPS QEq: charge equilibration (Rappé & Goddard, J. Phys. Chem. 95, 3358 (1991)) for any element, with the
// electronegativities χ and idempotentials J of the UFF / QEq table.
//
// Minimises E(q) = Σ χ_i q_i + ½ Σ J_i q_i² + Σ_{i<j} q_i q_j K(r_ij) at fixed total charge, with the shielded Coulomb
// kernel of Ohno and Klopman, K = 14.40 / √(r² + a_ij²) eV, a_ij = 7.20 (1/J_i + 1/J_j) Å, which tends to the
// point-charge 1/r at long range and stays finite at bond lengths; a 7th-order taper brings it to zero at the cut-off
// (summed over every periodic image within the cut-off, so a cell narrower than twice the cut-off keeps all of it). The linear system is solved by Jacobi-preconditioned conjugate gradients, twice
// (H s = −χ, H t = 1; q = s + μ t with μ fixing the total charge), as LAMMPS fix qeq does. Like the original QEq it
// gives ionic crystals large charges (quartz O about −1.6 e); molecules come out close to chemical intuition.
//
// For periodic fillers (silica, ZnO, clays, CaCO3) and the original model:
//  kernel "slater"  the Coulomb integral between the atoms' normalised ns Slater densities |φ|², φ ∝ r^(n−1) e^(−ζr),
//                   n the element's period, ζ = λ (2n + 1) / (2R) with R the QEq radius and λ = 0.4913 (Rappé & Goddard);
//                   tends to 1/r at long range, finite at r = 0. Computed exactly (incomplete gamma functions, then a
//                   radial quadrature) and tabulated per element pair. "point" is the bare 1/r (for checks).
//  ewald            (the default for a cell periodic along a, b and c, when it is affordable: by default a cell whose
//                   reciprocal sum would cost more than 2·10⁷ atom × k-vector terms per step takes the tapered kernel over
//                   the images within the cut-off instead, with a note; asked for explicitly, Ewald is always used)
//                   the long-range 1/r summed over all images by Ewald (as Wilmer, Kim & Snurr's EQeq): real space
//                   erfc(αr)/r, reciprocal space (4π/V) Σ_k e^(−k²/4α²)/k² cos k·r, self −2α/√π; the kernel's short-range
//                   difference from 1/r is added for pairs within the cut-off. The matrix is never stored: the reciprocal
//                   part is applied as structure factors.
//  hydrogen         J_H(q) = J_H⁰ (1 + q / ζ_H), ζ_H the hydrogen Slater exponent in bohr⁻¹ (Rappé & Goddard's charge-
//                   dependent hydrogen), solved self-consistently with mixing.
//  bounds           per element q_min ≤ q ≤ q_max: atoms beyond a bound are held there and the rest re-equalised, until
//                   the Karush–Kuhn–Tucker conditions hold (as LAMMPS qeq/ctip bounds charges).
#pragma once
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "caps/json.hpp"
#include "caps/system.hpp"

namespace caps {

struct QEqOptions {
  double total_charge = 0.0;   // e
  double cutoff = 10.0;        // Å (pairs over every periodic image within it, however narrow the cell)
  double tolerance = 1e-8;     // CG residual (relative)
  int max_iterations = 1000;
  std::string kernel = "ohno";   // ohno | slater | point
  bool ewald = true;             // cells periodic along a, b and c; false: the tapered kernel over the images within the cut-off
  bool ewald_auto = true;        // true: the default, Ewald skipped for a cell too large to sum it quickly; false: as asked
  double ewald_accuracy = 1e-7;  // erfc(α r_c) and the reciprocal cut-off's Gaussian
  double lambda = 0.4913;        // Slater exponent scale
  bool hydrogen = false;         // charge-dependent hydrogen hardness
  std::map<int, std::pair<double, double>> bounds;   // element → (q_min, q_max)
  int scf_iterations = 200;
  double scf_tolerance = 1e-6;   // e, largest charge change
};

struct QEqReport {
  int iterations = 0;
  double residual = 0;
  double cutoff = 0;
  double mu = 0;                 // the equalised electronegativity (eV)
  int scf = 0;                   // self-consistency passes (hydrogen, bounds)
  int bounded = 0;               // atoms held at a bound
  double alpha = 0;              // Ewald α (1/Å), 0 without Ewald
  int kvectors = 0;
  bool converged = true;         // false: the self-consistent passes, the bounds or the CG solve did not converge (see notes)
  std::vector<std::string> notes;
};

// Charges (e) for every atom of s. Throws std::invalid_argument for an element without QEq parameters.
std::vector<double> qeq_charges(const System& s, const QEqOptions& o = {}, QEqReport* report = nullptr);

// QEq parameters of an element from the UFF table: χ (eV), J (eV) and the QEq radius (Å); false when unknown.
bool qeq_parameters(int z, double& chi, double& J, double& radius);

// Options from JSON (ewald absent: the default, on for a cell periodic along a, b and c where affordable; given: as asked): {kernel, ewald, accuracy, lambda, hydrogen, cutoff, total, bounds: {"O": [qmin, qmax], …}}.
QEqOptions qeq_options_from_json(const Json& j);
// Bounds from text: "O:-1.2:0, Si:0:2.4".
std::map<int, std::pair<double, double>> qeq_bounds_from_text(const std::string& t);

// The Coulomb integral (hartree, for unit charges) between normalised ns Slater densities with principal quantum numbers
// n1, n2 and exponents ζ1, ζ2 (bohr⁻¹), centres R bohr apart.
double slater_coulomb_integral(int n1, double zeta1, int n2, double zeta2, double R);
// The period of an element (the Slater principal quantum number).
int principal_quantum_number(int z);
// The Coulomb energy ½ qᵀ K q (eV) of fixed charges under the options' kernel, cut-off and Ewald (no χ, no J): a check of
// the summation (a rock-salt crystal's Madelung energy with kernel "point").
double qeq_coulomb_energy(const System& s, const std::vector<double>& q, const QEqOptions& o = {});

}  // namespace caps
