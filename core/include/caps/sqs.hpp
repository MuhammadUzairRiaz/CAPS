// CAPS substitutional disorder: special quasirandom structures and symmetry-distinct occupations of a sublattice.
//
//  Sublattice   the sites of element `from` in a periodic cell (a supercell chosen by the user); a fraction x of them take
//               element `to`. Occupation σ_i = +1 (from) or −1 (to).
//  SQS          the arrangement whose pair correlation functions Π_m = ⟨σ_i σ_j⟩ over the first neighbour shells m (pairs
//               of the sublattice at equal distance) best match the random alloy's (1 − 2x)² (Zunger, Wei, Ferreira &
//               Bernard, Phys. Rev. Lett. 65, 353 (1990)): Σ_m w_m |Π_m − (1 − 2x)²|, w_m = 1/m, minimised by simulated
//               annealing of swaps (Kirkpatrick et al., Science 220, 671 (1983)).
//  Enumeration  every occupation with k substituted sites, one per orbit of the cell's symmetry operations (those that
//               map the sublattice onto itself, found from the undecorated cell): the canonical member is the
//               lexicographically smallest image; the orbit sizes are reported (they sum to C(n, k)).
#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "caps/json.hpp"
#include "caps/system.hpp"

namespace caps {

struct SqsOptions {
  int from = 0, to = 0;              // elements (atomic numbers); to = 0: vacancies are not allowed here
  double fraction = 0.5;             // of the sublattice that takes `to`
  int shells = 3;                    // pair shells matched
  int steps = 20000;
  int restarts = 4;
  uint64_t seed = 1;
};
struct SqsResult {
  System structure;                  // the cell with the substitution made
  std::vector<double> shell_distance, correlation;   // Å; Π_m of the result
  double target = 0;                 // (1 − 2x)²
  double objective = 0;
  int sites = 0, substituted = 0;
  std::vector<std::string> notes;
};
SqsResult special_quasirandom(const System& cell, const SqsOptions& o);

struct EnumOptions {
  int from = 0, to = 0;
  int count = 1;                     // substituted sites
  size_t max_configurations = 200000;   // C(n, k) above this is refused
  double tolerance = 0.1;            // Å, symmetry
};
struct EnumResult {
  std::vector<std::vector<int>> configurations;   // the substituted sites (indices into the sublattice), canonical, one per orbit
  std::vector<size_t> multiplicity;               // orbit sizes
  std::vector<uint32_t> sites;                    // the sublattice: atom indices of the cell
  int operations = 1;
  size_t total = 0;                               // C(n, k)
  std::vector<std::string> notes;
};
EnumResult enumerate_occupations(const System& cell, const EnumOptions& o);

// the cell with configuration `which` of an enumeration made (its sites take element `to`)
System occupied(const System& cell, const EnumResult& r, size_t which, int to);
// written by the extension: .cif, .vasp or POSCAR, otherwise PDB
void write_disordered(const System& s, const std::string& path);
Json sqs_json(const SqsResult& r);
Json enumeration_json(const EnumResult& r);

}  // namespace caps
