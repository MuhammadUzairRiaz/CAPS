// CAPS model resolution (design/boards/ModelResolution): the same structure at three resolutions, mass conserved.
//
//  united_atom    every hydrogen bonded to carbon folded into its carbon (CH3, CH2, CH sites; masses and charges summed;
//                 polar hydrogens kept, as TraPPE-UA does)
//  coarse_grain   each molecule's backbone cut into beads of `per_bead` backbone atoms; side groups and hydrogens join
//                 the bead of the backbone atom they hang on; a bead sits at its atoms' centre of mass and carries their
//                 mass (molecules without a backbone become one bead each)
//
// Both return new structures (the original is untouched). Going back up: united-atom → all-atom is add_hydrogens then
// a relax; coarse-grained → all-atom needs a backmap (not built here).
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
System united_atom(const System& s, ResolutionReport* rep = nullptr);
System coarse_grain(const System& s, int per_bead = 5, ResolutionReport* rep = nullptr);

}  // namespace caps
