// Typing by example: force-field types given by hand (or seeded automatically) on a small example — a repeat unit's
// head, body and tail in a short chain, a copolymer's units and their junctions — carried over to any structure made
// of the same chemistry by matching each atom's chemical environment.
//
// An atom's environment signature at radius r is a hash of its element, heavy degree, hydrogen count and smallest ring
// (from the bonds alone: files differ in whether they give bond orders, and with explicit hydrogens the counts imply
// the unsaturation), refined r times with its neighbours' signatures (Weisfeiler–Lehman / Morgan
// refinement, as extended-connectivity fingerprints: Rogers & Hahn, J. Chem. Inf. Model. 50, 742 (2010)). Two atoms
// with the same signature at radius r have the same bonded surroundings to r bonds.
//
// learn_types picks the smallest radius (from min_radius) at which no two example atoms with one signature carry
// different types; atoms of a structure take the type of their signature at that radius, and where the structure has
// an environment the example lacks (a branch point, a crosslink, a chain end not in the example), a shorter radius is
// tried when it is still unambiguous in the example; otherwise the atom is left untyped and reported.
#pragma once
#include <cstdint>
#include <map>
#include <string>
#include <vector>

#include "caps/system.hpp"

namespace caps {

// Signatures at radius 0 … radius (result[r][atom]).
std::vector<std::vector<uint64_t>> environment_signatures(const System& s, int radius);

struct ExampleTypes {
  int radius = 3;                                              // the radius the table is exact at
  std::vector<std::map<uint64_t, std::string>> by_radius;       // radius → signature → type (ambiguous ones left out)
  std::vector<std::string> conflicts;                          // environments the example types two ways even at the largest radius
  size_t environments = 0;                                     // distinct environments at `radius`
};

// Learns from example atoms with types ("" = not given; those atoms teach nothing). Throws when nothing is typed.
ExampleTypes learn_types(const System& example, const std::vector<std::string>& types, int min_radius = 2, int max_radius = 6);

struct ExampleMatch {
  std::vector<std::string> types;   // per atom, "" when no environment of the example fits
  size_t exact = 0, shorter = 0, unmatched = 0;
};
ExampleMatch apply_types(const System& s, const ExampleTypes& t);

// Atoms of the example with the same environment as `atom` at `radius` (for "assign to all equivalent atoms").
std::vector<uint32_t> equivalent_atoms(const System& s, uint32_t atom, int radius);

}  // namespace caps
