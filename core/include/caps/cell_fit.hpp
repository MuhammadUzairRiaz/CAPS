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

}  // namespace caps
