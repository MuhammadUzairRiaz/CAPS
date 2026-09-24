// CAPS polymer builder: chains of any repeat unit, grown into a periodic cell.
//
// A repeat unit is a SMILES with two attachment points, head first: "*CC(*)c1ccccc1" (styrene), "*CCO*" (ethylene
// oxide), "*O[Si](C)(C)*" (dimethylsiloxane). Each unit is embedded once (molecule builder, cleaned with a force field
// when one is given) with carbons in place of the attachment points, and turned into internal coordinates: bond,
// angle and torsion of every atom along a spanning tree from the head. A chain then grows unit by unit (NeRF
// placement), as Grow does for polystyrene: each step tries torsions for the link bond, the bond before it and the
// rotatable single bonds inside the unit (sp3–sp3: trans and gauche; next to an sp2 atom: 30° steps; ester and amide
// links: trans) and keeps the trial whose closest non-bonded contact is roomiest, against every atom already in the
// cell and across the periodic faces. Tacticity mirrors whole units: a meso dyad repeats the previous unit's
// configuration, a racemo dyad mirrors it. Chains end with hydrogen caps.
#pragma once
#include <cstdint>
#include <string>
#include <vector>

#include "caps/grow.hpp"
#include "caps/molecule.hpp"
#include "caps/system.hpp"

namespace caps {

struct RepeatUnit {
  std::string name;
  std::string smiles;   // two attachment points (* or [*:1] head, [*:2] tail)
};

// Homopolymer: unit A only. Alternating: A B C … A B C. Block: blocks of each unit in turn (lengths in `blocks`).
// Random: each position drawn with `weights` (Bernoulli / multinomial). Gradient: the share of A falls linearly
// from 1 at the head to 0 at the tail (two units; with more, the last unit rises). Pattern: `pattern` repeated
// ("AAB", letters A, B, C … for the units).
enum class Sequence { Homopolymer, Alternating, Block, Random, Gradient, Pattern };
Sequence sequence_from_string(const std::string& s);   // homopolymer | alternating | block | random | gradient | pattern
const char* to_string(Sequence s);

struct ChainSpec {
  std::vector<RepeatUnit> units;   // A, B, …
  Sequence sequence = Sequence::Homopolymer;
  int dp = 20;                     // units per chain
  std::vector<int> blocks;         // Block: lengths of the blocks, cycling through the units (A_n B_m …)
  std::vector<double> weights;     // Random: relative share of each unit (equal when empty)
  std::string pattern;             // Pattern: e.g. "AAB"
  Tacticity tacticity = Tacticity::Atactic;
  double pm = 0.5;                 // Atactic: probability of a meso dyad
  std::string forcefield;          // caps-forcefield JSON with typing rules for the unit templates (optional)
};

struct UnitInfo {
  std::string formula;             // of the repeat unit (without the attachment points)
  double mass = 0;                 // g/mol
  int atoms = 0;                   // with hydrogens
  int head = -1, tail = -1;        // atom indices in the SMILES
  int stereocentres = 0;           // atoms that tacticity acts on (backbone atoms with two different side groups)
  std::string head_element, tail_element;
};

// Parses and checks a repeat unit (two attachment points, each on one atom). Throws with the reason.
UnitInfo repeat_unit_info(const std::string& smiles);

// Unit index for each position of a chain (the sequence kind; Random draws from the seed).
std::vector<int> chain_sequence(const ChainSpec& spec, uint64_t seed);

// The chain as a molecule graph with hydrogen end caps (for its formula, mass and SMILES).
MolGraph chain_graph(const ChainSpec& spec, const std::vector<int>& sequence);
double chain_mass(const ChainSpec& spec, const std::vector<int>& sequence);

// Grows o.chains chains of spec (spec.dp units each) into a periodic cubic cell of edge o.box, or at o.density.
// Uses o.seed, o.trials, o.accept, o.contact_scale, o.max_restarts, o.progress; o.polymer and o.dp are ignored.
// Positions are unwrapped; atoms carry element names, molecule ids per chain and no charges (assign a force field in
// Field). Throws GrowError when a chain cannot be placed.
System grow_chains(const ChainSpec& spec, const GrowOptions& o, GrowReport* report = nullptr);

}  // namespace caps
