// CAPS React: bond-forming reactions from atom-mapped templates, run as cycles of react → retype → relax
// (the Polymatic cycle: Abbott, Hart and Colina, Theor. Chem. Acc. 132, 1334 (2013)), with distance capture and
// probability as in REACTER (Gissinger, Jensen and Wise, Polymer 128, 211 (2017)), and network analysis.
#pragma once
#include <cstdint>
#include <functional>
#include <stdexcept>
#include <string>
#include <vector>

#include "caps/field.hpp"
#include "caps/system.hpp"

namespace caps {

struct ReactError : std::runtime_error {
  using std::runtime_error::runtime_error;
};

// One atom of the pre-reaction pattern.
struct TemplateAtom {
  int map = 0;                  // map number (1-based, unique in the template)
  int element = 0;              // atomic number
  int h_min = -1, h_max = -1;   // hydrogen neighbours allowed (−1: any)
  int degree = -1;              // exact number of neighbours (−1: any)
  bool ring3 = false;           // in a three-membered ring (epoxide)
  bool not_aromatic = false;    // not in an aromatic six-ring
  std::vector<int> bonded;      // map numbers of pattern atoms it is bonded to (defined earlier)
};

// Atom-mapped reaction template. Text form (one statement per line, '#' comments):
//   reaction epoxy_amine
//   atom 1 C ring3                   atom 4 N H=2           atom 5 H bonded 4
//   atom 2 O ring3 bonded 1          atom 3 C ring3 bonded 1 2
//   initiators 4 1                   # the pair whose distance is tested
//   capture 4.5                      # Å
//   probability 1.0
//   min_path 4                       # initiators ≥ 4 bonds apart, or in different molecules (0: different molecules only)
//   form 4 1 · break 1 2 · delete 3 · move 5 2   (move: the atom leaves its partners and bonds to the second atom)
//   sites 1 2 3                      # the group counted for conversion: distinct matches of these map atoms (default: first initiator)
struct ReactionTemplate {
  std::string name;
  std::vector<TemplateAtom> atoms;
  int init_a = 0, init_b = 0;
  double capture = 4.5, probability = 1.0;
  int min_path = 4;
  std::vector<std::pair<int, int>> form, brk, move;
  std::vector<int> remove;
  std::vector<int> sites;       // map atoms that make up one counted reactive group
  std::string text;             // the source, for reports
};

std::vector<ReactionTemplate> parse_templates(const std::string& text);   // one or more "reaction" blocks
std::vector<std::string> builtin_template_names();                        // "cc_crosslink", "epoxy_amine_primary", …
std::string builtin_template(const std::string& name);                    // text of a built-in template

// A matched reaction site: system atom index for each template atom (same order as ReactionTemplate::atoms).
struct Match {
  int reaction = 0;
  std::vector<uint32_t> atoms;
  double distance = 0;
};

// All sites where a template's pattern matches with its initiators closer than the capture distance
// (minimum image). Sorted by distance.
std::vector<Match> find_matches(const System& s, const ReactionTemplate& t, int reaction_index = 0);
// Number of distinct reactive groups (matches of the template's site atoms) in the structure, for conversion.
int count_sites(const System& s, const ReactionTemplate& t);

// Apply non-overlapping matches (each atom reacts at most once per call). Deleted atoms are removed and the
// indices compacted; molecules are recomputed from bonds. Returns the number applied.
int apply_matches(System& s, const std::vector<ReactionTemplate>& t, const std::vector<Match>& m);

struct ClusterStats {
  int clusters = 0;
  double largest_fraction = 0;  // mass of the largest connected cluster / total
  double reduced_mw = 0;        // weight-average mass excluding the largest cluster, g/mol
  double mw = 0;                // weight-average mass, g/mol
};
ClusterStats cluster_stats(const System& s);

// Flory–Stockmayer critical conversion of the A groups for an A_fA + B_fB step-growth network with ratio
// r = (B groups)/(A groups) ≤ 1 taken on the A side: α_c = 1 / sqrt(r (fA − 1)(fB − 1)).
double flory_stockmayer(double r, double fa, double fb);

struct CycleRow {
  int cycle = 0, reactions = 0, total = 0;
  double conversion = 0;        // reactions so far / initial counted sites
  ClusterStats clusters;
  double energy = 0;            // after relaxation, kcal/mol
  int atoms = 0;
};

struct ReactOptions {
  std::vector<ReactionTemplate> templates;
  uint64_t seed = 1;
  int max_cycles = 100;
  int max_per_cycle = 10;       // reactions per cycle (Polymatic forms one bond per cycle; more is faster)
  double target_conversion = 1.0;
  int stall_cycles = 3;         // stop after this many cycles without a reaction (after dynamics, if any)
  bool relax = true;            // minimise after each cycle (needs every atom to be typed)
  int relax_iterations = 500;
  double relax_ftol = 2.0;
  double md_ps = 0;             // NVT after each cycle, ps (0: none)
  double temperature = 300;
  // REACTER-style (after Gissinger, Jensen & Wise, Polymer 128, 211 (2017)): one continuous NVT run, reactions checked
  // every md_ps; each reacted site (its atoms, those within two bonds and all atoms within 5 Å) is stabilised by a local minimisation
  // with capped forces while the rest of the cell is held and keeps its velocities — in place of LAMMPS's nve/limit
  // stabilisation. No global minimisation between checks. Stops at the cycle limit or the target conversion.
  bool during_md = false;
  EnergyOptions energy;
  std::function<bool(const CycleRow&)> progress;   // return false to cancel
  std::function<void(const System&, int cycle)> frame;
};

struct ReactReport {
  std::vector<CycleRow> cycles;
  int initial_sites = 0, reactions = 0;
  double gel_conversion = -1;   // conversion where the reduced weight-average mass peaked (−1: not seen)
  std::vector<std::string> notes;
  double seconds = 0;
};

// Runs cycles of find → react → retype → relax (→ dynamics) until the target conversion, the cycle limit, or no
// more reactions. `s` is changed in place.
void react(System& s, const ReactOptions& o, ReactReport* report = nullptr);

// The template editor (design/boards/ReactionTemplate): the pattern before and after the reaction as 2D drawings with
// map numbers, what changes, and checks. JSON {name, pre: {atoms: [{map, symbol, x, y, reacting}], bonds: [[m1, m2]]},
// post: {…}, changes: [{kind: formed|broken|deleted|moved, text}], checks: [{ok, text}], initiators: [a, b], capture,
// probability, min_path}.
std::string template_view(const ReactionTemplate& t);

}  // namespace caps
