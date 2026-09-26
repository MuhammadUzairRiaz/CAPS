// CAPS model resolution (design/boards/ModelResolution): the same structure at three resolutions, mass conserved.
//
//  united_atom    every hydrogen bonded to carbon folded into its carbon (CH3, CH2, CH sites; masses and charges summed;
//                 polar hydrogens kept, as TraPPE-UA does)
//  coarse_grain   each molecule's backbone cut into beads of `per_bead` backbone atoms; side groups and hydrogens join
//                 the bead of the backbone atom they hang on; a bead sits at its atoms' centre of mass and carries their
//                 mass (molecules without a backbone become one bead each)
//
// Both return new structures (the original is untouched). Going back up: united-atom → all-atom is add_hydrogens then
// a relax; coarse-grained → all-atom is backmap onto the all-atom structure the beads came from.
#pragma once
#include <string>
#include <vector>

#include "caps/system.hpp"

namespace caps {

struct ResolutionReport {
  int sites = 0;             // atoms, united sites or beads
  int hydrogens = 0;         // explicit hydrogens left
  double mass = 0;           // g/mol, the sum of the sites' masses (IUPAC standard atomic weights underneath)
  std::vector<int> site_of;  // for each original atom, its site in the result
  std::vector<std::string> notes;
};

// All-atom mass and counts of a structure (by element, standard atomic weights).
ResolutionReport all_atom_summary(const System& s);
// hosts: the elements whose hydrogens fold in (carbon by default; oxygen for a one-site water such as mW, site "OH2")
System united_atom(const System& s, ResolutionReport* rep = nullptr, const std::vector<int>& hosts = {6});
System coarse_grain(const System& s, int per_bead = 5, ResolutionReport* rep = nullptr);

// Backmapping (multiscale equilibration: coarse-grain here, run the beads elsewhere, come back): the beads of
// coarse_grain(all_atom, per_bead) have moved to `beads` (same count and order, e.g. read from a LAMMPS dump). Each
// bead's atoms are carried along rigidly: translated with the bead and turned by the rotation that best maps the bead
// and its bonded beads from their old places to their new ones (Horn's quaternion fit; a chain end turns with its one
// neighbour, a lone bead only moves). Bonds between beads come out stretched or squeezed: relax afterwards (push-off).
// The cell is the beads'. Throws when the bead counts differ.
struct BackmapReport {
  int beads = 0, atoms = 0;
  double rms_turn = 0;       // degrees, the beads' rotations
  double worst_bond = 0;     // Å, the longest bond between two beads' atoms before a relax
  std::vector<std::string> notes;
};
System backmap(const System& all_atom, const System& beads, int per_bead = 5, BackmapReport* rep = nullptr);

}  // namespace caps
