// CAPS QEq: charge equilibration (Rappé & Goddard, J. Phys. Chem. 95, 3358 (1991)) for any element, with the
// electronegativities χ and idempotentials J of the UFF / QEq table.
//
// Minimises E(q) = Σ χ_i q_i + ½ Σ J_i q_i² + Σ_{i<j} q_i q_j K(r_ij) at fixed total charge, with the shielded Coulomb
// kernel of Ohno and Klopman, K = 14.40 / √(r² + a_ij²) eV, a_ij = 7.20 (1/J_i + 1/J_j) Å, which tends to the
// point-charge 1/r at long range and stays finite at bond lengths; a 7th-order taper brings it to zero at the cut-off
// (minimum image in periodic cells). The linear system is solved by Jacobi-preconditioned conjugate gradients, twice
// (H s = −χ, H t = 1; q = s + μ t with μ fixing the total charge), as LAMMPS fix qeq does. Like the original QEq it
// gives ionic crystals large charges (quartz O about −1.6 e); molecules come out close to chemical intuition.
#pragma once
#include <string>
#include <vector>

#include "caps/system.hpp"

namespace caps {

struct QEqOptions {
  double total_charge = 0.0;   // e
  double cutoff = 10.0;        // Å (reduced to half the narrowest cell width when the cell is smaller)
  double tolerance = 1e-8;     // CG residual (relative)
  int max_iterations = 1000;
};

struct QEqReport {
  int iterations = 0;
  double residual = 0;
  double cutoff = 0;
  std::vector<std::string> notes;
};

// Charges (e) for every atom of s. Throws std::invalid_argument for an element without QEq parameters.
std::vector<double> qeq_charges(const System& s, const QEqOptions& o = {}, QEqReport* report = nullptr);

// QEq parameters of an element from the UFF table: χ (eV), J (eV) and the QEq radius (Å); false when unknown.
bool qeq_parameters(int z, double& chi, double& J, double& radius);

}  // namespace caps
