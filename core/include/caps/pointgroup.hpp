#pragma once
// Molecular point groups: the symmetry operations a molecule has within a tolerance, its Schoenflies point group, and
// the molecule made exactly symmetric.
//
// Candidate axes are the principal axes of inertia, the directions of atoms, of midpoints of like atoms and of normals
// to pairs of like atoms (the strategy of Beruski & Vidal, J. Comput. Chem. 35, 290 (2014)); an operation holds when
// every atom lands within the tolerance on an atom of its element. C_n (n ≤ 8), mirror planes, the inversion and S_2n
// are tested; the group follows the usual flowchart: linear (C∞v, D∞h), cubic (T, Td, Th, O, Oh) and icosahedral (I, Ih)
// from two or more C₃ axes, then Dnh, Dnd, Dn from n C₂ axes perpendicular to the principal one, then Cnh, Cnv, S2n,
// Cn, and Cs, Ci, C1. Symmetrisation closes the operations into the group and replaces each atom by the mean of its
// images, g⁻¹ x_π(i) over the group, repeated with the axes found again until the deviation is below 10⁻¹⁰ Å.
#include <array>
#include <string>
#include <vector>

#include "caps/system.hpp"

namespace caps {

struct SymmetryElement {
  std::string kind;            // "C", "S", "sigma", "i"
  int n = 1;                   // order (C_n, S_n)
  Vec3 axis{0, 0, 1};          // the axis, or the plane's normal
};
struct PointGroupReport {
  std::string group;           // Schoenflies: "C2v", "Td", "D6h", "C*v" for C∞v, "D*h" for D∞h
  int order = 1;               // operations in the group (0 for the linear groups)
  std::vector<SymmetryElement> elements;   // one per distinct axis or plane (the highest order on each axis)
  double max_deviation = 0;    // Å: the largest distance of an atom's image from its partner, over the group
  Vec3 centre{0, 0, 0};        // the centre of mass
  std::vector<std::string> notes;
};
// tol in Å. The structure is a molecule (the cell, if any, is ignored).
PointGroupReport point_group(const System& s, double tol = 0.1);
// Makes s exactly symmetric under the point group found within tol (positions only); returns the group after it.
PointGroupReport symmetrize_molecule(System& s, double tol = 0.1);

}  // namespace caps
