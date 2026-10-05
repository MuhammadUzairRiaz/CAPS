#pragma once
// Clipboard pieces and stamps (design/boards/Stamp): atoms copied out of a structure with their bonds, types, charges,
// names and molecules, as a small JSON document the Studio keeps in its tray; a stamp places a piece into a structure
// at a point, turned about an axis, and pushed clear of the atoms already there.
#include <string>
#include <vector>

#include "caps/system.hpp"

namespace caps {

// {"format":"caps-piece","name","formula","atoms":[{"el","x","y","z","q","type","label","mass","name","res","mol"}],
//  "bonds":[[i,j,order]]}: positions about the piece's centre; molecules renumbered from 1 in order of appearance.
std::string piece_json(const System& s, const std::vector<size_t>& atoms, const std::string& name);

struct StampResult {
  std::vector<size_t> added;   // the new atoms' indices
  double shift = 0.0;          // Å the piece was pushed to clear the structure
  double closest = 0.0;        // Å, its closest contact after
};
// Adds the piece centred at `at`, turned by `degrees` about `axis` (through its centre), then moved the shortest
// distance (up to 12 Å) that keeps every atom at least `clear` Å from the structure (minimum image); with no such place
// (a large piece in a dense cell) the roomiest one found, `closest` below `clear` (Minimise clears the contacts). Types are
// matched by label (else new types), molecule ids continue after the structure's. Throws for a piece that is not one.
StampResult stamp_piece(System& s, const std::string& piece, const Vec3& at, const Vec3& axis, double degrees, double clear);

}  // namespace caps
