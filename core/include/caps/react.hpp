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
  bool aromatic = false;        // in an aromatic six-ring (a phenol's ring carbon, not a carboxyl carbon)
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
//   charges keep                     # after the reaction: charges kept, a deleted atom's to its partner (default: forcefield)
//   byproduct 3 4                    # atoms that leave as a small molecule (H2, H2O, HCl): their bonds to the rest break,
//                                    # bonds among them (form 3 4) stay; the run removes them or keeps them as molecules
//   weight 2                         # relative rate against the other templates of a run (selection by weights)
struct ReactionTemplate {
  std::string name;
  std::vector<TemplateAtom> atoms;
  int init_a = 0, init_b = 0;
  double capture = 4.5, probability = 1.0;
  int min_path = 4;
  std::vector<std::pair<int, int>> form, brk, move;
  std::vector<int> remove;
  std::vector<int> sites;       // map atoms that make up one counted reactive group
  std::vector<int> byproduct;   // map atoms that leave as a small molecule (removed or kept, as the run chooses)
  double weight = 1.0;          // relative rate when a run selects among templates by weight
  bool keep_charges = false;    // "charges keep": the atoms keep their charges, a deleted atom's joins its partner (net charge
                                // conserved); default "charges forcefield": recomputed for the new chemistry
  std::string text;             // the source, for reports
};

std::vector<ReactionTemplate> parse_templates(const std::string& text);   // one or more "reaction" blocks
// A reaction scheme as atom-mapped SMILES — reactants and products, one entry per molecule — turned into a template. Map
// numbers tag the atoms whose bonds change; a mapped atom absent from the products leaves (with the hydrogens it takes:
// HCl, HBr); small separate products (water, an alcohol) are byproducts; hydrogens that change partner are moves. The
// pattern is the changing atoms by element, bonds and hydrogens and their neighbours by element and bonds — the model
// compound's other atoms stand for the chain. Two molecules meet per step: a scheme of three at once is refused (write it as
// steps). Throws ReactError with the reason.
ReactionTemplate template_from_scheme(const std::vector<std::string>& reactants, const std::vector<std::string>& products,
                                      const std::string& name, std::vector<std::string>* notes = nullptr);
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
// allow (optional): whether initiators a and b may react (a run's chain rule); capture > 0 replaces the template's.
std::vector<Match> find_matches(const System& s, const ReactionTemplate& t, int reaction_index = 0,
                                const std::function<bool(uint32_t, uint32_t)>& allow = {}, double capture = 0);
// Number of distinct reactive groups (matches of the template's site atoms) in the structure, for conversion.
int count_sites(const System& s, const ReactionTemplate& t);
// The distinct reactive groups themselves (each the sorted atoms matching the template's site atoms).
std::vector<std::vector<uint32_t>> site_groups(const System& s, const ReactionTemplate& t);

// The reactive sites of each chain for a set of templates (the chains as a run takes them: molecules of 30 atoms or more
// and a fifth of the largest): chain id (molecule, from 1) → sites, its repeat units (residues) and its mass (g/mol).
struct ChainSites {
  int64_t chain = 0;
  int sites = 0, units = 0, atoms = 0;
  double mass = 0;
};
// chains (optional, one per atom): the chains a reaction run started from (its chains_after), so chains joined by links are
// still counted one by one; empty: the molecules of the structure as it is
std::vector<ChainSites> chain_sites(const System& s, const std::vector<ReactionTemplate>& templates, const std::vector<int64_t>& chains = {});

// Apply non-overlapping matches (each atom reacts at most once per call). Deleted atoms are removed and the
// indices compacted; molecules are recomputed from bonds. Returns the number applied.
// keep_byproducts: a template's byproduct atoms stay as their own molecule (else removed); byproducts counts those kept or
// removed; tag (optional): a per-atom label carried through the compaction (the run's original chain of each atom).
// carry (optional): another per-atom label carried through the compaction unchanged (an atom's index in an earlier structure).
int apply_matches(System& s, const std::vector<ReactionTemplate>& t, const std::vector<Match>& m, bool keep_byproducts = false,
                  int* byproducts = nullptr, std::vector<int64_t>* tag = nullptr, std::vector<int64_t>* carry = nullptr);

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
  double max_force = 0;         // largest force after the cycle's relaxation, kcal/mol/Å (0: not relaxed)
  int atoms = 0;
  int crosslinks = 0;           // links between different chains so far
  double capture = 0;           // the capture distance this cycle used (Å; auto capture raises it)
  int target = 0;               // the crosslink target as links (0: a conversion target)
  double density = 0, degree = 0;   // ν so far (mol/m³), DC so far (%)
};

// What the run aims for: the conversion of the counted sites, or a number of links between chains given as a count, per
// chain, a crosslink density or a molecular weight between crosslinks.
// DegreePercent: the degree of crosslinking DC = 2 N_links / N_monomers × 100 % (Vasilev, Lorenz & Breitkopf, Polymers 13,
// 315 (2021); Alamfard et al., Polymers 15, 2058 (2023)), monomers counted as the chains' residues (the repeat units Grow numbers).
enum class ReactTarget { Conversion = 0, Crosslinks = 1, PerChain = 2, Density = 3, Mc = 4, DegreePercent = 5 };

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
  // a cycle that fails (a relaxation that cannot converge, dynamics that blow up) leaves the structure as the last
  // completed cycle did and stops there: the report says which cycle and why (failed_cycle, failure). False: throw.
  bool keep_on_failure = true;
  // the force field for the structure as the reactions leave it (the user's assignment re-run on the product: types,
  // charges and parameters for the new bonds); null: the built-in default (GAFF for C and H, UFF otherwise)
  std::function<std::shared_ptr<const ForceField>(const System&)> retype;
  std::string field_name;       // for the report
  // bonds only between different chains: each atom keeps the chain it started in; a small molecule (a curative, a
  // crosslinker) belongs to the chains it has bonded to, so ENR(A)–MAH cannot close back onto chain A
  bool between_chains = false;
  // each atom's chain from an earlier run (the report's chains_after), so a cure in several runs keeps the chains it started
  // from; empty (or another atom count): the molecules of the structure as it is now
  std::vector<int64_t> chains;
  // a label per atom carried through the run unchanged (e.g. each atom's index in the structure the run started from);
  // the report's carry_after holds it for the atoms at the end. Empty: nothing carried
  std::vector<int64_t> carry;
  bool keep_byproducts = false; // byproduct atoms kept as molecules (else removed)
  // at most this many of each chain's reactive sites react (a template's site atoms on that chain; 0: no limit) — e.g. 2 of
  // the 5 epoxides of a 10-unit ENR-50 chain; a site reacting in the template's step on a small molecule is not counted
  int sites_per_chain = 0;
  // several templates: 0 closest pairs first whatever the template (by distance); 1 by relative weights (weights[k], else
  // the template's weight) — each pick chooses a template in proportion to its weight, then its closest free pair
  int selection = 0;
  std::vector<double> weights;
  // auto capture: when a cycle finds no pair, every template's capture grows by capture_step, up to capture_max
  bool auto_capture = false;
  double capture_step = 0.5, capture_max = 8.0;
  ReactTarget target = ReactTarget::Conversion;
  double target_value = 0;      // links, links per chain, mol/m³, g/mol (target_conversion for Conversion)
  EnergyOptions energy;
  std::function<bool(const CycleRow&)> progress;   // return false to cancel
  std::function<void(const System&, int cycle)> frame;
  // after every cycle, for a live view: the structure, its row, each atom's chain of the start (negative: a byproduct) and
  // the atoms of the links formed between chains so far
  std::function<void(const System&, const CycleRow&, const std::vector<int64_t>& chain, const std::vector<char>& linked)> live;
};

// One link between two chains: where it is (each chain and repeat unit — the residue numbers Grow gives — 0 when the chain
// carries none) and the molecule that bridges them (a crosslinker's start id, 0 for a direct bond).
struct LinkRecord {
  int cycle = 0;
  std::string reaction;
  int64_t chain_a = 0, chain_b = 0, unit_a = 0, unit_b = 0, via = 0;
  std::string via_name;         // the bridge's element formula at the start (C4H2O3 for maleic anhydride)
};

struct ReactReport {
  std::vector<CycleRow> cycles;
  int initial_sites = 0, reactions = 0;
  double gel_conversion = -1;   // conversion where the reduced weight-average mass peaked (−1: not seen)
  std::vector<std::string> notes;
  double seconds = 0;
  int failed_cycle = 0;         // > 0: this cycle failed; the structure is the one after cycle failed_cycle − 1
  std::string failure;
  // the network: chains (molecules of the start with ≥ 30 atoms and ≥ 20 % of the largest), links between different chains
  // and within one; crosslink density ν = links / (V N_A); strands 2 × links (tetrafunctional junctions); Mc = chain mass / strands
  int chains = 0, crosslinks = 0, intrachain = 0, byproducts = 0;
  int target_crosslinks = 0;    // the target as a number of links (0: a conversion target)
  double volume = 0, chain_mass = 0;   // Å³, g/mol
  double density = 0, per_chain = 0, mc = 0;   // mol/m³, 2 × links / chains, g/mol
  int monomers = 0;             // repeat units of the chains (residues), 0 when the chains carry no residue numbers
  double degree = 0;            // DC = 2 × links / monomers × 100 %
  std::string field;            // the force field that relaxed the network
  std::vector<int64_t> chains_after;   // each atom's chain at the end (for the next run's ReactOptions::chains)
  std::vector<int64_t> carry_after;    // ReactOptions::carry for the atoms at the end
  std::vector<LinkRecord> links;       // every link between chains, in the order they formed
};

// Runs cycles of find → react → retype → relax (→ dynamics) until the target conversion, the cycle limit, or no
// more reactions. `s` is changed in place.
void react(System& s, const ReactOptions& o, ReactReport* report = nullptr);

// The template editor (design/boards/ReactionTemplate): the pattern before and after the reaction as 2D drawings with
// map numbers, what changes, and checks. JSON {name, pre: {atoms: [{map, symbol, x, y, reacting}], bonds: [[m1, m2]]},
// post: {…}, changes: [{kind: formed|broken|deleted|moved, text}], checks: [{ok, text}], initiators: [a, b], capture,
// probability, min_path}.
std::string template_view(const ReactionTemplate& t);

// The template as an atom-mapped reaction SMARTS (reactants>>products): the query atoms with their constraints
// ([#6;X4;!H0;A:1] — element, connections, hydrogens, three-ring, aliphatic), the products by element and map number;
// deleted atoms are absent from the products (removed, as RDKit reads it). The pattern's bonds are connectivity only, written
// '~' (any bond), except a formed bond between atoms that keep their bond count (a substitution): single, '-'.
std::string reaction_smarts(const ReactionTemplate& t);

}  // namespace caps
