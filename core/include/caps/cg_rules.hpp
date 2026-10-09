// CAPS chemistry-aware coarse-grained mapping (the coarse-graining workflow's "Mapping" stage): beads from rules about the
// chemistry, not from counting atoms, so several polymers that share chemistry share beads.
//
//   cut       bonds to cut, as SMARTS whose first and last atoms are bonded (that bond is cut: "[CX3](=O)-[OX2;!H1]" cuts
//             every ester C(=O)–O, the carbonyl O only giving context): the connected fragments left are the beads, end
//             groups staying in their fragment
//   beads     fragment SMARTS → bead name: an exact cover of each molecule's heavy atoms (caps/resolution.hpp
//             map_to_beads), hydrogens with their heavy atoms
//   explicit  an atom → bead index list
//
// Beads are named by an ordered list of (name, SMARTS): the first whose pattern embeds in the fragment's heavy atoms names
// it (so "S" = "[CX3](=O)[CH2][CH2][CX3]=O" names a succinate bead with or without an acid end). Every fragment also
// gets a class: its formula, aromaticity and a canonical key (a Weisfeiler–Lehman hash of the heavy-atom graph labelled
// by element, aromaticity and hydrogen count). A fragment no rule names is reported with its class and stops the mapping
// (strict, the default), or is named by its class (F1, F2 … in the order met) when strict is off; it is never merged
// into a neighbour.
//
// A bead sits at the centre of mass of its atoms (default), their centre of geometry, or on one atom (the first atom of
// the first embedding of a SMARTS in the fragment). Chains are ordered from end to end (the end with the lower first
// atom index first). Types (bead, bond, angle, dihedral) are numbered by a shared list so several systems get the same
// numbering and can share table files; without one, the list is the mapping's own, sorted.
#pragma once
#include <cstdint>
#include <functional>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "caps/json.hpp"
#include "caps/system.hpp"

namespace caps {

struct CgNameRule { std::string name, smarts; };

struct CgRules {
  std::vector<std::string> cut;                // bond SMARTS (cut form)
  std::vector<CgNameRule> beads;               // fragment SMARTS → name (exact-cover form)
  std::vector<int> explicit_bead;              // atom → bead (explicit form; −1 not allowed)
  std::vector<CgNameRule> names;               // naming rules, first match wins (cut and explicit forms)
  std::string position = "com";                // com | cog | atom:<SMARTS>
  bool strict = true;                          // an unnamed fragment stops the mapping
};

// The shared type list: "B", "A-B", "A-B-A", "A-B-S-B" (bond, angle and dihedral keys in canonical order: the smaller of
// the forward and the reversed sequence).
struct CgTypes {
  std::vector<std::string> beads, bonds, angles, dihedrals;
  bool empty() const { return beads.empty(); }
};
std::string cg_key(const std::vector<std::string>& seq);   // canonical key of a bonded sequence
CgTypes cg_types_from_json(const Json& j);
Json cg_types_json(const CgTypes& t);
// The union of several lists, each part sorted (a list for several systems mapped alike).
CgTypes cg_types_union(const std::vector<CgTypes>& lists);

struct CgFragmentClass {
  std::string key;          // canonical key (hex)
  std::string formula;      // C4H8O2
  bool aromatic = false;
  std::string name;         // bead name given to it
  int count = 0;
};

struct CgMapping {
  std::vector<std::vector<uint32_t>> bead_atoms;   // atoms of each bead (all-atom indices)
  std::vector<std::string> bead_kind;              // its name
  std::vector<std::string> bead_class;             // its class key
  std::vector<int> bead_mol;                       // its molecule (0-based, as System::molecules numbers them)
  std::vector<double> bead_mass;                   // g/mol
  std::vector<int> bead_anchor;                    // position "atom:…": the atom it sits on, else −1
  std::vector<std::vector<int>> chains;            // beads of each molecule in chain order
  std::vector<std::pair<int, int>> bonds;          // bead bonds (from the atom bonds between beads)
  std::vector<std::array<int, 3>> angles;
  std::vector<std::array<int, 4>> dihedrals;
  std::vector<CgFragmentClass> classes;
  std::string position = "com";
  int atoms = 0;                                   // all-atom count it was made for
  double mass_aa = 0, mass_cg = 0;
  bool branched = false;
  std::vector<std::string> notes;

  size_t beads() const { return bead_atoms.size(); }
  std::vector<int> bead_of() const;               // atom → bead
};

// Throws std::invalid_argument with the unnamed fragments (their class, formula and count) under strict naming.
CgMapping cg_mapping(const System& aa, const CgRules& rules);

// map.json: {"format": "caps-cg-map", "version": 1, "position", "atoms", "beads": [{"kind", "class", "mol", "mass",
// "atoms": [1-based ids …], "anchor"}], "chains": [[bead indices …]], "bonds", "classes"}; reading it back gives the same
// mapping (angles and dihedrals rebuilt from the bonds).
Json cg_mapping_json(const CgMapping& m, const System& aa);
CgMapping cg_mapping_from_json(const Json& j, const System& aa);
CgRules cg_rules_from_json(const Json& j);

// Bead positions for one frame of the all-atom structure: each bead made whole about its first atom (minimum image),
// then its centre of mass / geometry / anchor atom.
std::vector<Vec3> cg_positions(const CgMapping& m, const System& aa, const std::vector<Vec3>& pos, const Cell& cell);

// The beads as a structure (element 0, named by kind, typed by `types` — the mapping's own sorted list when empty;
// throws when a kind or bonded sequence is not in a given list), bonded, angles and dihedrals in System::topology when
// the structure keeps an explicit topology, molecule ids from the all-atom molecules.
System cg_structure(const CgMapping& m, const System& aa, const std::vector<Vec3>& beads, const Cell& cell, const CgTypes& types = {});
// The mapping's own type list (sorted).
CgTypes cg_types_of(const CgMapping& m);

// What a mapping holds: per kind the count, per molecule the sequence (composition of kinds), mass conservation.
struct CgMappingSummary {
  std::map<std::string, int> kinds;
  std::map<std::string, double> fraction;          // share of each kind among the beads
  double mass_error = 0;                           // |Σ bead mass − Σ atom mass|
  int chains = 0, min_beads = 0, max_beads = 0;
  std::vector<std::string> sequences;              // the first few chains as kind strings (B S B A …)
};
CgMappingSummary cg_mapping_summary(const CgMapping& m, int sequences = 3);

}  // namespace caps

namespace caps {

// A LAMMPS data file of the beads (atom_style full, charges 0): Masses, Atoms, Bonds, Angles, Dihedrals, every type
// numbered by `types` (the mapping's own sorted list when empty) and named in a comment, so tables can be shared.
void write_cg_lammps_data(const CgMapping& m, const System& cg, const CgTypes& types, const std::string& path);

// Streaming trajectory mapping (one frame in memory): a LAMMPS text dump of the all-atom structure — "id" and either
// xu yu zu, x y z (with ix iy iz when present), or xs ys zs / xsu ysu zsu; any column order, atoms matched by id — mapped
// bead by bead and written as a LAMMPS dump "id mol type xu yu zu" of the beads. Without unwrapped coordinates or image
// flags the beads are made whole along each chain's bead bonds (minimum image from its first bead).
struct CgTrajectoryOptions {
  size_t first = 0, last = SIZE_MAX, stride = 1;   // frames kept (0-based, inclusive)
  CgTypes types;                                  // bead type numbers (empty: the mapping's own)
  std::function<bool(double fraction, size_t frames)> progress;   // false stops (the frames written are kept)
};
struct CgTrajectoryReport {
  size_t frames_read = 0, frames_written = 0;
  double seconds = 0;
  bool unwrapped_input = false, images = false, stopped = false;
  std::vector<std::string> notes;
};
CgTrajectoryReport map_lammps_dump(const CgMapping& m, const System& aa, const std::string& dump, const std::string& out,
                                   const CgTrajectoryOptions& o = {});

}  // namespace caps
