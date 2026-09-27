// CAPS functionalisation: chemical groups grafted covalently onto fillers — the sidewalls, ends and edges of nanotubes and
// sheets (graphene, h-BN), fullerenes and any sp² framework, or in place of the hydrogens of an edge or surface.
//
// A group is a SMILES with one attachment point (*), or a preset name (hydroxyl, carboxyl, amine …). The sites follow
// a pattern:
//   random  a share (or a count) of the free sidewall sites drawn at random
//   all     every free sidewall site, taken in turn, at least min_spacing apart (an ordered pattern)
//   band    as random, only within a window along the axis (a tube's z, a sheet's x), from … to as fractions
//   helix   as all, only on a helix round a tube: angle 2π z / pitch + phase, within 20°
//   ends    the hydrogens of a tube's open ends (the rim atoms), replaced
//   edges   the hydrogens of a flake's edges, replaced
//   atoms   the atoms given (their hydrogen replaced when they have one, else on their free side)
// A sidewall group sits on the site's outward normal (the normal of its three bonds, turned away from the filler's
// centre: outside a tube or particle, on top of a sheet; side "inner" turns it in, "both" draws each at random). Sites
// are at least min_spacing apart, so no two neighbours both carry a group. The grafted sp² atom keeps its position:
// relax afterwards to let it pyramidalise (sp³).
#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "caps/system.hpp"

namespace caps {

struct FunctionalizeOptions {
  std::string group = "hydroxyl";   // preset name or SMILES with one *
  std::string pattern = "random";   // random | all | band | helix | ends | edges | atoms
  std::string elements;             // site elements, e.g. "C" or "B,N" ("" any sp² framework atom)
  double fraction = 0.05;           // of the eligible sites (random, band)
  int count = 0;                    // exactly this many (overrides fraction)
  double min_spacing = 3.0;         // Å between grafted sites
  double from = 0.0, to = 1.0;      // band: window along the axis, as fractions of the filler's length
  double pitch = 20.0, phase = 0.0; // helix: Å per turn, degrees at z = 0
  std::string side = "outer";       // outer | inner | both
  std::vector<uint32_t> atoms;      // pattern "atoms"
  uint64_t seed = 1;
};

struct FunctionalizeReport {
  size_t eligible = 0, grafted = 0, added_atoms = 0;
  double degree = 0;                // grafted sites per 100 framework atoms
  std::vector<uint32_t> sites;
  std::vector<std::string> notes;
};

// The preset groups: hydroxyl, carboxyl, amine, methyl, fluoro, phenyl, nitrophenyl, amide, ester, hydroxymethyl,
// vinyl, thiol, aminopropyl, octadecylamide, peg3; any other string is taken as SMILES.
std::string functional_group_smiles(const std::string& name);
const std::vector<std::string>& functional_group_names();

FunctionalizeReport functionalize(System& s, const FunctionalizeOptions& o);

}  // namespace caps
