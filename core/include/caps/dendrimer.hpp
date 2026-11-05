#pragma once
// Dendrimers and hyperbranched molecules from AB_n branching units, each a SMILES with its sites marked:
//   [*:1]  the unit's A site (the one bond toward the core or the focal unit),
//   [*:2]  its n B sites (where the next units hang).
// A core is a SMILES with f B sites [*:2]; an end group a SMILES with one A site [*:1] that caps every B site left
// open. Without an end group an open site stays a hydrogen (NH₂ ends of PAMAM, OH ends of bis-MPA). Units are joined
// in 3D through their sites (ring_edit.hpp join_pieces: the two site atoms leave, the new bond at the covalent sum,
// each unit turned about it for the most room), then UFF relaxes the whole molecule.
//
//   perfect  every B site of every shell carries a unit, `shells` shells around the core;
//   sma      slow monomer addition: units one at a time, each onto a B site drawn uniformly from all open ones (on the
//            core too); for AB₂ the degree of branching tends to 2/3 (rate equations: the terminal, linear and dendritic
//            fractions all tend to 1/3; Hanselmann, Hölter & Frey, Macromolecules 31, 3790 (1998) for the method);
//   flory    Flory's statistics of AB_n polycondensation (J. Am. Chem. Soc. 74, 2718 (1952): equal reactivity, no
//            cycles): seen from the focal unit a molecule is a Galton–Watson tree in which each B site reacts
//            independently; a molecule of exactly N units is that tree conditioned on its size, drawn exactly by the
//            sequence-and-rotation method (Devroye, SIAM J. Comput. 41, 1 (2012); the cycle lemma of Dvoretzky &
//            Motzkin, Duke Math. J. 14, 305 (1947)); for AB₂ the degree of branching tends to 1/2 (Hölter, Burgath &
//            Frey, Acta Polym. 48, 30 (1997): it cannot exceed 0.5 in random AB₂ polycondensation). No core.
//
// Degrees of branching (AB₂): Frey's DB = 2D / (2D + L) (Hölter et al. 1997) and Fréchet's DB = (D + T) / (D + T + L)
// (Hawker, Lee & Fréchet, J. Am. Chem. Soc. 113, 4583 (1991)), D, L, T the units with 2, 1, 0 of their B sites
// carrying a unit (end groups do not count). For n > 2 the report gives the units by the number of B sites used.
//
// Presets (core, unit, focal group; open sites stay H): pamam (ethylenediamine core, generation g = g + 1 shells; Tomalia et al.,
// Polym. J. 17, 117 (1985)), ppi (1,4-diaminobutane core, generation g = g shells; de Brabander-van den Berg & Meijer,
// Angew. Chem. Int. Ed. 32, 1308 (1993)), bismpa (trimethylolpropane core, 2,2-bis(hydroxymethyl)propionic acid
// ester units, generation g = g shells; Ihre, Hult & Söderlind, J. Am. Chem. Soc. 118, 6388 (1996)).
#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

#include "caps/system.hpp"

namespace caps {

struct DendrimerPreset {
  std::string name, title, core, unit, focal;
  int shell_offset = 0;   // shells = generation + offset
};
const std::vector<DendrimerPreset>& dendrimer_presets();

struct DendrimerSpec {
  std::string preset;           // fills core and unit (a core or unit given here wins)
  std::string core, unit, end;
  std::string focal;            // without a core: a group on the focal unit's A site, a SMILES with one [*:2] ("[*:2]O"
                                // makes bis-MPA's acid); empty: a hydrogen
  int generation = -1;          // with a preset: the family's generation number (≥ 0); else −1 and `shells`
  int shells = 2;               // perfect: unit shells around the core
  std::string growth = "perfect";   // perfect | sma | flory
  int units = 30;               // sma, flory: units in the molecule (the core not counted)
  uint64_t seed = 1;
  bool relax = true;
  int relax_steps = 2000;
};

// The branching tree: unit k hangs on B site slot[k] of unit parent[k] (parent −1: on core site slot[k], or the focal
// unit when there is no core); units are in the order they are built (every parent before its children).
struct DendrimerTree {
  std::vector<int> parent, slot;
};
DendrimerTree perfect_tree(int core_sites, int n, int shells);
DendrimerTree sma_tree(int units, int n, int core_sites, uint64_t seed);
DendrimerTree flory_tree(int units, int n, uint64_t seed);

struct BranchCounts {
  std::vector<int> by_used;     // units with k of their n B sites carrying a unit, k = 0 … n
  int dendritic = 0, linear = 0, terminal = 0;   // AB₂
  double db_frey = NAN, db_frechet = NAN;          // AB₂
  std::vector<int> per_shell;   // units per shell (shell 1 on the core or the focal unit)
};
BranchCounts branch_counts(const DendrimerTree& t, int n);

struct DendrimerReport {
  int units = 0, shells = 0, n = 0, core_sites = 0, ends = 0;
  BranchCounts counts;
  double db_limit = NAN;        // the large-N degree of branching (Frey) of the growth, AB₂
  std::string formula;
  double mass = 0;              // g/mol
  std::vector<std::string> notes;
};
System build_dendrimer(const DendrimerSpec& spec, DendrimerReport* report = nullptr);

}  // namespace caps
