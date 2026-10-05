// CAPS conformer search: the low-energy minima of a molecule (or a few) under its force field, in vacuum.
//
//  Torsions   every rotatable bond (acyclic, single, a heavy atom beyond each end) set at random to one of its three
//             staggered positions (60°, 180°, 300°) with a little noise, then minimised (L-BFGS); the first trial is the
//             structure as it is.
//  Anneal     NVT dynamics at a high temperature with snapshots quenched to their minima (a quench / anneal search, the
//             distribution of minima it visits as a by-product).
//  Clusters   the minima in energy order; one joins the first lower conformer within `rmsd` Å over the heavy atoms after
//             optimal superposition (Horn 1987), else starts a conformer of its own. Conformers above the lowest by more
//             than `window` kcal/mol are dropped. Populations are Boltzmann weights of the minimised energies (no
//             vibrational entropy) at `temperature`.
#pragma once
#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "caps/field.hpp"
#include "caps/system.hpp"

namespace caps {

struct ConformerSearchOptions {
  std::shared_ptr<const ForceField> field;   // null: CAPS's default force field
  EnergyOptions energy;                      // the cell is ignored: molecules in vacuum
  std::string method = "torsions";           // torsions | anneal
  int trials = 50;
  double ftol = 0.05;                        // kcal/mol/Å
  double window = 10.0;                      // kcal/mol above the lowest kept
  double rmsd = 0.5;                         // Å, heavy atoms
  double temperature = 298.15;               // K, populations
  double anneal_temperature = 1000.0;        // K
  double anneal_ps = 1.0;                    // ps between quenched snapshots
  uint64_t seed = 1;
  std::function<bool(double)> progress;      // fraction → false cancels
};

struct ConformerHit {
  std::vector<Vec3> pos;
  double energy = 0;         // kcal/mol
  double relative = 0;       // above the lowest
  double population = 0;     // Boltzmann, sums to 1
  int found = 0;             // minima that fell into it
};

struct ConformerSearchResult {
  std::vector<ConformerHit> conformers;   // lowest first
  std::vector<double> minima;             // every trial's minimised energy, kcal/mol, in the order found
  int rotors = 0;
  std::vector<std::array<uint32_t, 4>> rotor_atoms;
  std::string field, method;
  std::vector<std::string> notes;
};

// The rotatable bonds as dihedrals a-b-c-d (b–c the bond; a and d heavy neighbours).
std::vector<std::array<uint32_t, 4>> rotatable_bonds(const System& s);

ConformerSearchResult conformer_search(const System& s, const ConformerSearchOptions& o);

}  // namespace caps
