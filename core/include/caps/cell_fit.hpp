#pragma once
// Is a periodic cell large enough for its bonds? Force fields read a pair's separation at the nearest image and exclude or
// scale 1-2, 1-3 and 1-4 partners there, and a bond list holds one entry per pair of atoms. Both assumptions fail in a
// small crystal cell, which is common: a polyethylene cell repeats its chain every 2.5 Å, so an atom is bonded to two images
// of one neighbour and partners three bonds apart sit closer through the next cell. Energies are then wrong by far more
// than any effect being studied, and no warning follows from the numbers alone. The test walks every 1-3 and 1-4 path and
// compares it with the image the engines would use; it asks nothing about widths, so a wide cell holding a stretched
// network fails and a narrow cell of small molecules (a noble gas, a dense molecular crystal) passes. Coarse-grained beads
// (element 0) are left out: their models exclude what they choose, often bonded pairs alone.
#include <array>
#include <cstddef>
#include <string>

#include "caps/system.hpp"
#include "caps/typing.hpp"

namespace caps {

struct CellFit {
  bool fits = true;
  size_t double_images = 0;   // bonded pairs (or an atom and itself) within bonding distance at more than one image
  size_t wrong_images = 0;    // 1-3 and 1-4 partners reached by bonds at an image other than the nearest (or a tie)
  std::array<bool, 3> axis{false, false, false};   // the cell directions at fault
};
// field: only the partners it excludes or scales count (a bead-spring model excludes bonded pairs alone); none: all of them
CellFit cell_fit(const System& s, const ForceField* field = nullptr);

// The smallest repetition that fits, found by growing the axes at fault; throws std::invalid_argument past max_atoms.
std::array<int, 3> copies_to_fit(const System& s, size_t max_atoms = 15000);

// What is wrong and what to do, for error messages and reports ("" when the cell fits).
std::string cell_fit_message(const System& s, const ForceField* field = nullptr);

// assign_types for a cell too small for its bonds: the cell's bonds hold one entry per pair of atoms, so a carbon of a
// polyethylene cell bonded to two images of one neighbour would be typed with three neighbours (an alkene). The copies
// that fit (copies_to_fit) are typed instead and each atom takes the type of its own copy (the first: a supercell keeps
// atom i of that copy at i), with a note naming the copies. A cell that fits is typed as it is. `copies`, when given,
// receives the repetition used ({1, 1, 1} for a fit). Throws std::invalid_argument when no repetition within 15000 atoms
// fits, with cell_fit_message's reason.
TypingResult assign_types_fitted(const System& s, const FFDef& ff, std::array<int, 3>* copies = nullptr);

}  // namespace caps
