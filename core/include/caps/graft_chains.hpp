// CAPS grafting chains onto an existing structure: chains of their own units grown from chosen atoms of a backbone, a
// cell or a molecule (NR-g-PMMA, PP-g-MA, compatibilisers, long-chain branches added afterwards). (Silane coupling
// agents on silanols are edit.hpp's graft_silanes.)
//
//  Sites      atoms given (0-based; a hydrogen stands for the atom it sits on), each optionally standing for every atom
//             equivalent to it: the same element, the same number of hydrogens, the same neighbour elements and the
//             same residue name. Each candidate is kept with `probability`, at least `min_spacing` apart (minimum
//             image), up to `max_grafts` in all.
//  Grafts     each site gives up `per_site` of its hydrogens (the one pointing into most room first); a chain of the
//             graft spec grows from there along the hydrogen's direction, the host held fixed (grow_chains with the host
//             as substrate and the sites as anchors). A host without a cell is put in a box with room around it.
//  Molecules  renumbered by connectivity afterwards: a graft belongs to the molecule it hangs on.
#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "caps/grow.hpp"
#include "caps/polymer.hpp"
#include "caps/system.hpp"

namespace caps {

struct ChainGraftOptions {
  std::vector<uint32_t> sites;         // atoms of the host (0-based)
  bool equivalent = false;             // each site stands for every atom equivalent to it
  double probability = 1.0;            // per candidate site
  int per_site = 1;                    // grafts on one site atom (hydrogens it gives up)
  int max_grafts = 0;                  // 0: no limit
  double min_spacing = 0;              // Å between sites
  double padding = 0;                  // a host without a cell: room around it (0: from the graft length)
  GrowOptions grow;                    // seed, trials, method …; cell, chains and anchors are set here
};

struct ChainGraftReport {
  int candidates = 0;                  // sites before the draw
  int grafted = 0;
  std::vector<uint32_t> sites;         // host atoms carrying grafts (once per graft), host numbering
  GrowReport grow;
  std::vector<std::string> notes;
};

// The atoms equivalent to `atom` (itself included), by the rule above.
std::vector<uint32_t> equivalent_atoms(const System& s, uint32_t atom);

System graft_chains(const System& host, const ChainSpec& graft, const ChainGraftOptions& o, ChainGraftReport* report = nullptr);

}  // namespace caps
