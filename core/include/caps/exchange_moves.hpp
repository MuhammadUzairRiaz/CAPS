// CAPS Exchange, moves: Monte Carlo moves that break and make bonds during a run (stage E2 of
// docs/studies/vitrimers/caps_module_proposal.md), on the TopologyEditor of stage E1. Three kinds cover the exchange
// chemistries:
//
//   Swap   associative exchange of partners, as LAMMPS fix bond/swap (Sides et al. 2004) and the hybrid MD/MC vitrimer
//          models (Branham-Ferrari et al. 2026, Mamat et al. 2026): bonds a1–a2 and b1–b2, all four atoms in group A
//          and every distance a1–a2, b1–b2, a1–b2, b1–a2 within the capture distance, become a1–b2 and b1–a2. With
//          labels, a1 and b1 must carry the same one and so must a2 and b2 (as LAMMPS checks molecule IDs: a bead's place
//          along its chain keeps chain lengths; one label for every sticker lets any two stickers exchange).
//   Shift  associative exchange with valences (addition–elimination; multivalent linkers, polysulfide sulfur, the
//          pivot–residue model of Rao, Xia & Ni 2024): a pivot of group A bonded to a residue r0 of group B takes a
//          residue r1 of B within the capture distance that has a valence left, and lets r0 go.
//   Template  any exchange a React template describes (the schemes of data/reactions, data/exchange): its sites found by
//          chemistry within the template's capture distance, its bonds broken and made as it says, accepted with the
//          ratio of the matches before and after (the way back is the same template on the products). A hydrogen the
//          template moves goes onto its new partner at a covalent bond length, on the side it came from (as React places
//          it), and needs the relaxation below. Types follow roles: the template matched the other way on the products
//          names the atom that plays each role after the change, and it takes the type and charge of the atom that played
//          the role before (the oxygen an ester leaves becomes an alcohol's; atoms outside the pattern are paired through
//          their neighbours' roles). With ExchangeOptions::parameterize, types come from the run's own assignment of the
//          product's neighbourhood instead (roles need not be equivalent: an ester moving from a tertiary carbon to a
//          secondary one), and charges still follow roles. min_path 0, React's "different molecules", reads here as at least 4 bonds apart, so
//          a network exchanges within itself. Templates that delete atoms run in LAMMPS (fix bond/react) instead.
//   Bind   dissociative exchange (Diels–Alder, hindered ureas, reversible stickers): a site of A and a site of B (or two of
//          A when B is empty) within the capture distance bond, each while it has a valence left, or a bond between them
//          breaks; the binding constant K = exp(ε/k_B T) weighs bonded against free.
//
// Every move is Metropolis–Hastings with its proposal counts, so detailed balance holds: acceptance
// exp(−E_a/k_B T) · min[1, (q_reverse/q_forward) · K^±1 · exp(−ΔU/k_B T)], the barrier E_a (the exchange chemistry's
// activation energy) slowing every move alike without biasing the equilibrium. ΔU comes from the atoms within three bonds
// of the change (in the graph before and after it), whose energy difference equals the whole system's: every bonded
// term and every exclusion a change touches lies among them, and the rest (pairs with other atoms, k-space sums, tail
// corrections) is the same before and after. Force fields with many-body terms, hydrogen bonds or virtual sites use the
// whole system's energy instead.
//
// Relaxed moves (relax_steps > 0; all-atom templates): a new bond made at the capture distance is stretched by 1 to 3 Å,
// so its energy turns down nearly every move. With relax_steps, the atoms of the change, the atoms bonded to them and
// their hydrogens relax by steepest descent after it, the rest held in place: every bonded term within three bonds of
// them, and every atom within the cut-off of them, enters that energy. ΔU is then the relaxed energy against the energy
// before, and an accepted move keeps the relaxed positions. This is the relaxed bond-exchange criterion of the
// atomistic vitrimer studies (Yang et al. 2016, Wick et al. 2025), not exact detailed balance: the positions change.
#pragma once
#include <cstdint>
#include <functional>
#include <memory>
#include <random>
#include <string>
#include <vector>

#include "caps/dynamics.hpp"
#include "caps/field.hpp"
#include "caps/react.hpp"
#include "caps/system.hpp"
#include "caps/topology_edit.hpp"

namespace caps {

enum class ExchangeKind { Swap, Shift, Bind, Template };
ExchangeKind exchange_kind_from_string(const std::string& s);   // "swap" | "shift" | "bind" | "template"
const char* to_string(ExchangeKind k);

struct ExchangeOptions {
  ExchangeKind kind = ExchangeKind::Swap;
  std::vector<uint32_t> group_a;     // swap: the exchangeable atoms; shift: the pivots; bind: the sites (and their partners when group_b is empty)
  std::vector<uint32_t> group_b;     // shift: the residues; bind: the partner sites
  std::vector<int64_t> labels;       // swap: per atom; a1 and b1 must share theirs, and a2 and b2 (empty: any)
  int valence_a = 1, valence_b = 1;  // shift: bonds a residue (b) can hold; bind: dynamic bonds a site of a / b can hold
  double capture = 1.5;              // Å
  double barrier = 0;                // kcal/mol, E_a
  double binding = 0;                // kcal/mol, ε (bind): K = exp(ε/k_B T)
  int attempts = 1;                  // moves tried at each call
  uint64_t seed = 1;
  bool local = true;                 // ΔU from the atoms near the change (false: the whole system's energy)
  int relax_steps = 0;               // > 0: relaxed moves (see above), this many steepest-descent steps; needed by templates that move a hydrogen
  int relax_depth = 3;               // bonds from the change within which atoms relax (with the hydrogens on them)
  // The run's force-field assignment. Each change's neighbourhood, six bonds deep with hydrogen caps on the cut bonds
  // (none across a double or aromatic bond), is typed with it: the atoms within three bonds of the change take its types,
  // and the terms of combinations of types the structure never had are learned (TopologyEditor::learn; polysulfides
  // exchanging into a disulfide make the first C–S–S–C torsion). An assignment with an overlay or with CAPS ML
  // completion (caps ff complete: a torsion the force field lacks fitted to the ML potential) fills what the force field
  // itself has no term for. Null, or an assignment that fails: the change is turned down as one the force field cannot
  // describe (the reason in ExchangeStats::note).
  std::function<ForceField(const System&)> parameterize;
  std::vector<ReactionTemplate> templates;   // template: the exchanges to try (one picked at random each move)
};

struct ExchangeStats {
  int64_t attempts = 0, accepted = 0;
  int64_t no_candidate = 0;          // tries with no partner within the capture distance
  int64_t unparameterised = 0;       // tries the force field could not describe (rejected)
  int64_t relaxed = 0;               // tries whose ΔU came after a local relaxation
  int64_t learned = 0;               // combinations of types learned from the force-field assignment
  int64_t new_types = 0;             // types the assignment gave that the run's force field had for no atom before
  double energy = 0;                 // ΔU of the accepted moves, summed, kcal/mol
  bool local = true;                 // ΔU was local (false: the whole system's, because the force field needs it)
  std::string note;
};

class ExchangeEngine {
 public:
  // s and ff are the run's structure and force field (edited in place; they must outlive the engine). eo: the run's
  // energy settings (cut-off, electrostatics), for the local energies.
  ExchangeEngine(System& s, ForceField& ff, const ExchangeOptions& o, const EnergyOptions& eo);
  ~ExchangeEngine();
  // o.attempts moves at positions x (3N, Å) in `cell` at `temperature` (K). potential(): the whole system's potential
  // energy with the force field as it stands (needed when the energies are not local; may be null otherwise). Returns
  // true when a move was accepted.
  bool run(const std::vector<double>& x, const Cell& cell, double temperature, const std::function<double()>& potential);
  // The same, and an accepted relaxed move writes its relaxed positions into x (relaxed moves need this form).
  bool run(std::vector<double>& x, const Cell& cell, double temperature, const std::function<double()>& potential);
  // ΔU of a change (bonds broken and made) at x, without keeping it: local when the engine is, else from potential().
  double delta(const std::vector<TopologyEditor::Pair>& broken, const std::vector<TopologyEditor::Pair>& made, const std::vector<double>& x,
               const Cell& cell, const std::function<double()>& potential);
  const ExchangeStats& stats() const { return st_; }
  TopologyEditor& editor() { return ed_; }

 private:
  bool try_swap(const std::vector<double>& x, const Cell& cell, double beta, const std::function<double()>& potential);
  bool try_shift(const std::vector<double>& x, const Cell& cell, double beta, const std::function<double()>& potential);
  bool try_template(const std::vector<double>& x, const Cell& cell, double beta, const std::function<double()>& potential);
  bool try_bind(const std::vector<double>& x, const Cell& cell, double beta, const std::function<double()>& potential);
  bool decide(const std::vector<TopologyEditor::Pair>& broken, const std::vector<TopologyEditor::Pair>& made, double log_ratio, double beta,
              const std::vector<double>& x, const Cell& cell, const std::function<double()>& potential);
  // applies the change and returns ΔU (throws TopologyError, with nothing changed, when the force field cannot describe it)
  double apply(const std::vector<TopologyEditor::Pair>& broken, const std::vector<TopologyEditor::Pair>& made, const std::vector<double>& x,
               const Cell& cell, const std::function<double()>& potential);
  double distance(const std::vector<double>& x, const Cell& cell, uint32_t i, uint32_t j) const;
  // ed_.change, learning the terms of new combinations of types first when it has to (o_.parameterize)
  void edit(const std::vector<TopologyEditor::Pair>& broken, const std::vector<TopologyEditor::Pair>& made, const std::vector<TopologyEditor::Retype>& retype);
  double cluster_energy(const std::vector<uint32_t>& atoms, const std::vector<double>& x, const Cell& cell) const;
  ForceField cluster_field(const std::vector<uint32_t>& atoms) const;   // the terms among these atoms, renumbered
  double apply_relaxed(const std::vector<TopologyEditor::Pair>& broken, const std::vector<TopologyEditor::Pair>& made, const std::vector<double>& x,
                       const Cell& cell);

  System& s_;
  ForceField& ff_;
  ExchangeOptions o_;
  EnergyOptions eo_;
  TopologyEditor ed_;
  std::vector<char> in_a_, in_b_;
  std::mt19937_64 rng_;
  ExchangeStats st_;
  std::unique_ptr<ForceField> shell_;   // the force field without atoms and terms: the base of each cluster's
  std::vector<double>* xw_ = nullptr;   // the run's positions, when they may change (run with a mutable x)
  std::vector<TopologyEditor::Pair> moves_;   // the trial's moved hydrogens and their new partners
  std::vector<TopologyEditor::Retype> retype_, unretype_;   // the trial's types following roles, and the way back
  std::vector<std::pair<uint32_t, Vec3>> pending_;   // a relaxed trial's positions, written when it is accepted
};

// A DynamicsOptions::topology callback running an ExchangeEngine on the run's own structure and force field (made at the
// first call). stats, when given, follows the engine's.
std::function<bool(TopologyContext&)> exchange_callback(const ExchangeOptions& o, const EnergyOptions& eo,
                                                        std::shared_ptr<ExchangeStats> stats = nullptr);

// Labels for swaps along linear chains, as LAMMPS fix bond/swap reads molecule IDs: each atom's place along its chain
// counted from one end (1, 2, 3 … N: swaps keep every chain's length, one partner per bead), or from the nearer end when
// `from_both_ends` (1, 2, 3, 3, 2, 1: swaps may happen at either end). Atoms not on a linear chain get 0.
std::vector<int64_t> chain_position_labels(const System& s, bool from_both_ends = false);

}  // namespace caps
