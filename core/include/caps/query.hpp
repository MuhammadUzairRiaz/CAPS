// Selection queries (design/boards/SmartSelect): one grammar for the Studio's query bar, the command palette, macros and
// Python (caps_select mode "query").
//
//   smarts "c1ccccc1"          atoms matched by a SMARTS pattern
//   element C N O              by element          type c3 ca     by force-field type label
//   chain 1-4 | molecule 2,5   molecules (1-based, in file order)      index 1-20   atoms (1-based)
//   ring 5                     the atoms of the 5th ring (SSSR, in the order perception finds them)
//   stereo R | S | *           tetrahedral stereocentres: CIP R, S, or all four-different centres
//   within 5.0 of <query>      atoms within the distance (Å, minimum image) of the query's atoms
//   sel | selection | all | none
//   and · or · not · ( )       combine (not binds tightest, then and, then or)
//
// CIP labels rank the four substituents by atomic number sphere by sphere (duplicate atoms for double and triple bonds,
// rules 1a/1b); a centre whose substituents tie within 16 spheres (a backbone CH between two long identical chain
// ends) is a stereocentre (*) without an R/S label.
#pragma once
#include <string>
#include <vector>

#include "caps/system.hpp"

namespace caps {

struct QueryResult {
  std::vector<char> atoms;
  int rings = 0;   // SSSR rings whose atoms are all selected
};
QueryResult select_query(const System& s, const std::string& query, const std::vector<char>& current = {});

// CIP labels of the tetrahedral centres: 'R', 'S', '*' (four different neighbours but a tie), 0 elsewhere.
std::vector<char> cip_labels(const System& s);

}  // namespace caps
