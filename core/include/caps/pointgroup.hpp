#pragma once
// Molecular point groups: the symmetry operations a molecule has within a tolerance, its Schoenflies point group, and
// the molecule made exactly symmetric.
//
// Candidate axes are the principal axes of inertia, the directions of atoms, of midpoints of like atoms and of normals
// to pairs of like atoms (the strategy of Beruski & Vidal, J. Comput. Chem. 35, 290 (2014)); an operation holds when
// every atom lands within the tolerance on an atom of its element (a C_n only when its permutation has order n).
// C_n (n ≤ 8), mirror planes, the inversion and S_2n are tested. Linear molecules give C∞v or D∞h; otherwise the
// operations are closed on the atoms' permutations, every operation of the closed group must hold within the tolerance
// (else the largest subgroup found whose operations all hold is kept, with a note), and the Schoenflies name is read
// from that group's order, rotation orders, reflections and inversion, so name, order and miss describe one group.
// Symmetrisation replaces each atom by the mean of its images, g⁻¹ x_π(i) over the group, repeated with the axes found
// again until the deviation is below 10⁻¹⁰ Å.
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
  double max_deviation = 0;    // Å: the largest distance of an atom's image from its partner, over the group (≤ tol)
  Vec3 centre{0, 0, 0};        // the centre of mass
  std::vector<std::string> notes;
};
// tol in Å. The structure is a molecule (the cell, if any, is ignored).
PointGroupReport point_group(const System& s, double tol = 0.1);
// Makes s exactly symmetric under the point group found within tol (positions only); returns the group after it.
PointGroupReport symmetrize_molecule(System& s, double tol = 0.1);

}  // namespace caps
