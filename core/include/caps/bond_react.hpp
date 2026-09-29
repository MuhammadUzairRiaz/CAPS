// CAPS reactions as a LAMMPS fix bond/react set (REACTER: Gissinger, Jensen & Wise, Polymer 128, 211 (2017);
// Macromolecules 53, 9953 (2020)): for every template, pre- and post-reaction molecule templates and a map file
// (InitiatorIDs, EdgeIDs, DeleteIDs, Equivalences) cut from real reaction sites of the structure, typed with the
// structure's force field before and after the reaction — the post types, charges and new bonded terms are the ones the
// force field gives the reacted structure — a data file whose type and coefficient tables hold every type the reactions
// create, and an input script with the fix. LAMMPS matches a template by its atom types and bonds, so each distinct
// chemical environment of a reactive site (a chain end nearby, a different neighbour type) needs its own template: the
// most frequent ones are written, and the report says how many of the candidate sites they cover.
#pragma once
#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "caps/field.hpp"
#include "caps/react.hpp"
#include "caps/relax.hpp"
#include "caps/system.hpp"

namespace caps {

struct BondReactOptions {
  int radius = 3;               // bonds from the reacting atoms kept in a template (edge atoms beyond)
  int max_variants = 6;         // templates per reaction (the most frequent environments)
  double survey_capture = 0;    // Å: candidate sites for the survey (0: the larger of each template's capture and 8 Å)
  bool keep_byproducts = false; // byproduct atoms stay as molecules (else DeleteIDs)
  bool between_chains = false;  // molecule inter: initiators in different molecules (the data file's molecule ids)
  std::vector<double> weights;  // relative rates per template: prob = template probability × weight / largest weight
  int nevery = 100;             // steps between reaction checks
  // Rmax of every reaction in LAMMPS, Å (0: the template's capture). LAMMPS stabilises only the reacting atoms, so a bond
  // formed from far apart (with a proton jumping to its new partner) can blow the run up; REACTER's own examples use about
  // 3 Å and let diffusion bring groups together, where CAPS's cycles relax the whole cell and can reach further
  double rmax = 3.5;
  double temperature = 300;     // K
  int64_t steps = 100000;
  double timestep = 0;          // fs (0: the force field's own)
  double stabilization_xmax = 0.03;   // Å per step: nve/limit on the reacting atoms (stabilization yes)
  uint64_t seed = 12345;
  EnergyOptions energy;
  LammpsStyle style;
};

struct BondReactReport {
  std::vector<std::string> files;
  std::vector<std::string> notes;
  struct Variant { std::string reaction, name; int sites = 0, pre_atoms = 0, edge = 0, deleted = 0; };
  std::vector<Variant> variants;
  int candidates = 0, covered = 0;
};

// field: the force field of a structure (the user's assignment re-run on it); required. The files are written to
// dir/stem.data, dir/stem.in, dir/stem_<reaction>_<k>_pre.mol, …_post.mol, …_map.txt.
BondReactReport write_bond_react(const System& s, const std::vector<ReactionTemplate>& templates,
                                 const std::function<std::shared_ptr<const ForceField>(const System&)>& field, const std::string& dir,
                                 const std::string& stem, const BondReactOptions& o = {});

// The other way: a fix bond/react set (pre- and post-reaction molecule templates and the map file, CAPS's or anyone's)
// as a CAPS reaction template. Elements come from the masses: the molecule file's Masses section, else masses_from (a
// LAMMPS data file whose Masses section numbers the same types). The pattern is the pre-reaction template: every atom by
// element, hydrogen count and — away from the edge — its number of bonds; bonds formed and broken from comparing the two
// templates through the Equivalences; DeleteIDs deleted; a hydrogen that leaves one partner for another becomes a move.
// LAMMPS matches by atom type, CAPS by element and connectivity: the notes say what the pattern keeps.
ReactionTemplate read_bond_react(const std::string& pre, const std::string& post, const std::string& map, const std::string& masses_from = "",
                                 const std::string& name = "", double capture = 0, std::vector<std::string>* notes = nullptr);
// A template as the text parse_templates reads.
std::string template_text(const ReactionTemplate& t);

}  // namespace caps
