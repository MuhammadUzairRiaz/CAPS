// CAPS reactions as a LAMMPS fix bond/react set (REACTER: Gissinger, Jensen & Wise, Polymer 128, 211 (2017);
// Macromolecules 53, 9953 (2020)): for every template, pre- and post-reaction molecule templates and a map file
// (InitiatorIDs, EdgeIDs, DeleteIDs, Equivalences) cut from real reaction sites of the structure, typed with the
// structure's force field before and after the reaction — the post types, charges and new bonded terms are the ones the
// force field gives the reacted structure — a data file whose type and coefficient tables hold every type the reactions
// create, and an input script with the fix. LAMMPS matches a template by its atom types and bonds, so each distinct
// chemical environment of a reactive site (a chain end nearby, a different neighbour type) needs its own template: the
// most frequent ones are written, and the report says how many of the candidate sites they cover.
#pragma once
#include <cstdint>
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
  // steps each reaction's atoms stay under nve/limit after it happens (react … stabilize_steps; REACTER's default is 60,
  // too few for a hydrogen that moves to a new partner a few Å away)
  int stabilize_steps = 200;
  // a hydrogen that changes partner (an acid H onto an epoxide O): the reaction waits until it is within this distance of
  // its new partner (a distance constraint in the map file; 0: none), so the new O–H bond is not formed stretched
  double h_transfer_max = 3.5;
  uint64_t seed = 12345;
  EnergyOptions energy;
  LammpsStyle style;
  // Survey after a virtual cure: CAPS React runs on scratch copies of the structure to these conversions (the same
  // templates, between_chains and weights; relaxed after each cycle when survey_relax), and the candidate sites of every
  // template are surveyed in the structure and in each copy — so a step whose groups appear only as the cure goes on (a
  // half-ester's acid, a maleate's second acid) gets its templates. The data file stays the unreacted structure.
  std::vector<double> survey_after;
  bool survey_relax = true;
  // Atom types split by component (type_group: one per atom of the structure, −1 none; names in the order their types are
  // numbered): in the data file, every template and the input's groups (split_types_by_group)
  std::vector<int> type_group;
  std::vector<std::string> type_group_names;
  // Crosslink-density targets: the input runs in chunks of check_every steps, sums the reactions of link_reactions (globs
  // over the reaction names written: "enr_acid_ester_*"; empty: the second steps where a template has two, else every
  // reaction), writes crosslink_progress.dat (step, time, each reaction, links, crosslink density = links / limiting,
  // ν) and stem_XL<target>.data (and a restart) when each target is reached; stops at the last target or max_steps.
  // targets as fractions (0.2) or percent (20).
  std::vector<double> targets;
  double limiting = 0;
  std::vector<std::string> link_reactions;
  int64_t check_every = 1000, max_steps = 2000000;
  // no reaction for stall_chunks chunks: every Rmax grows by rmax_step, up to rmax_limit (0: never)
  int stall_chunks = 0;
  double rmax_step = 0.5, rmax_limit = 0;
  // molecule ids during the cure: "reset" (LAMMPS's default: the bonded pieces after each reaction), "keep"
  // (reset_mol_ids no: the data file's molecules stay), "molmap" (reset_mol_ids molmap, LAMMPS 2 Apr 2025 or later: the
  // chains keep their ids and a small molecule takes the id of the chain it first bonds to — so molecule inter means
  // different original chains, and a crosslinker cannot close back on the chain it hangs from)
  std::string mol_ids = "reset";
};

struct BondReactReport {
  std::vector<std::string> files;
  std::vector<std::string> notes;
  struct Variant { std::string reaction, name, step; int sites = 0, pre_atoms = 0, edge = 0, deleted = 0; };
  std::vector<Variant> variants;
  int candidates = 0, covered = 0;
  // coverage per reaction step ("first": on a small molecule not yet bonded; "second": on one bonded at its other end; ""
  // when a reaction has one step), over the structure and every surveyed copy
  struct Step { std::string reaction, step; int candidates = 0, covered = 0, variants = 0; };
  std::vector<Step> steps;
  struct Frame { double conversion = 0; int reactions = 0, atoms = 0; };   // the virtual-cure copies surveyed
  std::vector<Frame> frames;
  std::vector<std::string> link_reactions;   // the reactions counted as links in the input (targets)
  int h_transfers = 0;                       // hydrogen-transfer distance constraints written
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
