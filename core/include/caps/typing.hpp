// CAPS atom typing: chemical perception (bond orders, formal charges, rings, aromaticity) and SMARTS typing rules.
//
// A force field's types carry SMARTS rules (FFType.smarts). An atom gets the type whose rule matches it with the atom
// as the pattern's first atom; when several match, a type listed in another matching type's `overrides` drops out,
// then the highest `priority` wins, then the earliest rule in the file. Rules may name other types with %type
// (as in foyer): typing repeats until no type changes.
//
// Supported SMARTS (Daylight): organic-subset and bracket atoms; primitives * a A #n element (aromatic lower case),
// D X H h v R r x (with counts), charges + - +n -n, recursive $(...), %type; logic ! & , ; ; bonds - = # : ~ @ and
// their negation; branches; ring closures (digits and %nn). Chirality and isotopes are ignored. CAPS extension:
// {AR1}..{AR5}, antechamber's ring classes (Perception::ar_class), for GAFF.
//
// Conjugated type pairs (FFDef::typing_pairs, e.g. GAFF cc/cd): after typing, atoms carrying either type of any pair
// are coloured along their bonds (double bond: the other type of the pair, single: the same), starting from the first
// type at the lowest-numbered atom of each conjugated system, as antechamber does. With pair_mode "double_same"
// (CGenFF's CG2DC1 / CG2DC2) a double bond keeps the type and a conjugated single bond switches it.
#pragma once
#include <string>
#include <vector>

#include "caps/ffdef.hpp"
#include "caps/system.hpp"

namespace caps {

// Chemistry of a structure as the typing rules see it.
struct Perception {
  std::vector<std::vector<uint32_t>> nb;     // neighbours
  std::vector<std::vector<int>> order;       // bond order to each neighbour (1 2 3; aromatic bonds keep a Kekulé order)
  std::vector<std::vector<bool>> arom_bond;  // aromatic bond to each neighbour
  std::vector<bool> aromatic;
  std::vector<int> charge;                   // formal charge
  std::vector<int> hcount;                   // attached hydrogens (explicit, plus implicit)
  std::vector<int> implicit_h;               // hydrogens a united-atom site carries (CH2, CH3, CR1 ...), 0 otherwise
  bool united_atom = false;                  // no hydrogen on any carbon, and carbon sites named for their hydrogens
  std::vector<int> ring_count;               // number of SSSR rings containing the atom (SMARTS R)
  std::vector<int> smallest_ring;            // size of the smallest ring (SMARTS r), 0 if acyclic
  std::vector<int> ring_bonds;               // ring bonds at the atom (SMARTS x)
  std::vector<std::vector<bool>> ring_bond;  // bond to each neighbour is in a ring
  std::vector<std::vector<uint32_t>> rings;  // SSSR
  std::vector<int> ar_class;                 // antechamber ring class 1..5 (SMARTS extension {AR1}..{AR5}), 0 acyclic
  std::vector<std::string> notes;
  int bond_order(uint32_t a, uint32_t b) const;
};

// Bond orders come from the file when it gives them (mol2), otherwise from valences. Hydrogens are explicit, except in
// a united-atom structure (no hydrogen bonded to any carbon): there a carbon or sulfur site named for the hydrogens it
// carries (CH, CH0-CH4, CH1E-CH3E, CR1, CR1E, SH1E; as united_atom() and united-atom force fields name them) has them
// as implicit hydrogens, counted by SMARTS H and X and in the valences.
Perception perceive(const System& s);
// Hydrogens a site named like a united atom carries (-1: not such a name).
int united_atom_hydrogens(const std::string& name, int element);

// A compiled SMARTS pattern.
class Smarts {
 public:
  explicit Smarts(const std::string& text);
  ~Smarts();
  Smarts(Smarts&&) noexcept;
  Smarts& operator=(Smarts&&) noexcept;
  // Does the pattern match with its first atom on `atom`? `types` (may be empty) serves %type references.
  bool matches(const System& s, const Perception& p, uint32_t atom, const std::vector<std::string>& types) const;
  bool uses_types() const;
  const std::string& text() const;

 private:
  struct Impl;
  Impl* d_;
};

struct TypingResult {
  std::vector<std::string> types;             // one per atom ("" when no rule matched)
  std::vector<std::string> why;               // the rule that assigned it, or the reason there is none
  std::vector<int> rule;                      // index into FFDef::typing of that rule, -1 when none matched
  std::vector<std::vector<std::string>> candidates;   // every matching type, before precedence
  std::vector<std::string> notes;
  int untyped = 0, ambiguous = 0;
};

// Types every atom with the force field's rules. Throws FFError when the force field has no typing rules.
TypingResult assign_types(const System& s, const FFDef& ff);

// Bond-order variants (FFDef::type_variants), run by assign_types after the rules.
void refine_bond_order_variants(const System& s, const Perception& p, const FFDef& ff, std::vector<std::string>& types, std::vector<std::string>& why);

}  // namespace caps
