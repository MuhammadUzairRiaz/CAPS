#pragma once
// Bond rules (design/boards/BondRules): bonds from what the structure actually holds. Per pair of elements, a histogram
// of the distances between their atoms (0.8 – 3.2 Å by default); the cut-off goes in the gap between bonded and
// non-bonded pairs, where it can be seen, and an ionic pair (Pauling electronegativity difference ≥ 1.7: Zn–O, Na–Cl)
// starts as never bonded.
#include <string>
#include <vector>

#include "caps/system.hpp"

namespace caps {

struct PairHistogram {
  int za = 0, zb = 0;             // za ≤ zb
  std::vector<int> counts;        // atom pairs per bin
  double lo = 0.8, bin = 0.04;    // Å
  size_t bonded = 0;              // the structure's bonds between these elements now
  double covalent = 0;            // the sum of the covalent radii, Å
  double suggested = 0;           // the cut-off in the first gap after the bonded peak (else the radii + 0.45 Å)
  bool ionic = false;
};
std::vector<PairHistogram> pair_histograms(const System& s, double lo = 0.8, double hi = 3.2, double bin = 0.04);

struct BondRule {
  int za = 0, zb = 0;
  double max = 0;                 // Å: bonded within it
  bool never = false;
};
// Bonds by the rules: a pair of elements with a rule bonds within its max (or never); the others within their covalent
// radii + 0.45 Å, as perception does. Pairs closer than min_distance are left alone (overlaps, not bonds).
std::vector<Bond> bonds_by_rules(const System& s, const std::vector<BondRule>& rules, double min_distance = 0.4);

}  // namespace caps
