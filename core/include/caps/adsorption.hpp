// CAPS adsorption locator: the low-energy places of adsorbate molecules on a substrate (a surface slab, a crystal, a
// polymer or filler surface), by Monte Carlo simulated annealing (Kirkpatrick, Gelatt & Vecchi, Science 220, 671 (1983)),
// the method of Materials Studio's Adsorption Locator.
//
// The substrate stays fixed; each adsorbate molecule moves as a rigid body (translation, rotation about its centre, a
// jump to a random place and orientation) under the Metropolis rule at a temperature that falls geometrically from
// t_high to t_low in each cycle; each cycle ends with a zero-temperature polish. The energy is the force field's
// non-bonded interaction between molecules (van der Waals shifted at the cut-off, damped-shifted-force electrostatics:
// pairmodel.hpp); the adsorbates' internal energy is constant (rigid) and left out, as is the substrate's.
//
// Reported for the lowest configuration: the adsorption energy (adsorbate–substrate plus adsorbate–adsorbate
// interaction, rigid: no deformation term), each component's dE/dN (the interaction of one of its molecules with all
// else, averaged over its molecules), the lowest configurations of the cycles, and the distribution of the energies
// sampled at the cold end of the cycles.
#pragma once

#include <cstdint>
#include <functional>
#include <string>
#include <vector>

#include "caps/field.hpp"
#include "caps/system.hpp"

namespace caps {

struct AdsorptionOptions {
  std::vector<int> mobile;          // molecules (System::molecules() components, from 0) that move; empty: first_mobile_atom
  int first_mobile_atom = -1;       // the molecules of this atom and every later one move (the adsorbates put after the
                                    // substrate); < 0: every molecule but the largest
  int cycles = 3;
  int steps = 20000;                // Monte Carlo steps per cycle
  double t_high = 1e4, t_low = 100; // K
  double cutoff = 12.0;             // Å
  bool coulomb = true;
  double dsf_alpha = 0.2;           // Å⁻¹
  double z_lo = 0, z_hi = 0;        // the adsorbates' atoms stay between these heights (Å); both 0: anywhere in the cell
  bool randomise = true;            // place the adsorbates at random before the first cycle
  int keep = 10;                    // lowest configurations kept
  uint64_t seed = 1;
  std::function<bool(int cycle, int step, double best)> progress;   // about every 1 %; false stops
};

struct AdsorptionConfig {
  double energy = 0;                // kcal/mol
  int cycle = 0;
  std::vector<Vec3> positions;      // every atom
};

struct AdsorptionComponent {
  std::string name;                 // the molecule's formula
  int molecules = 0;
  double de_dn = 0;                 // kcal/mol
};

struct AdsorptionReport {
  double adsorption_energy = 0;     // kcal/mol: the lowest configuration's total interaction (rigid)
  double adsorbate_substrate = 0, adsorbate_adsorbate = 0;
  std::vector<AdsorptionComponent> components;
  std::vector<AdsorptionConfig> configs;           // lowest first
  std::vector<double> hist_edges, hist_counts;     // energies sampled in the cold fifth of each cycle
  double acceptance = 0;
  long steps = 0;
  double seconds = 0;
  std::vector<std::string> notes;
};

// Moves the adsorbates of s to the lowest configuration found. ff: the force field assigned to s (substrate and
// adsorbates typed together). Throws std::invalid_argument without a mobile molecule or with pair forms other than LJ.
void locate_adsorption(System& s, const ForceField& ff, const AdsorptionOptions& o, AdsorptionReport* report = nullptr);

// The molecules (System::molecules() ids) holding atom first_atom or a later one; first_atom < 0: all but the largest.
std::vector<int> adsorbate_molecules(const System& s, int first_atom);

// The interaction energy (kcal/mol) of the mobile molecules with everything else in s, as the locator counts it.
double adsorption_interaction(const System& s, const ForceField& ff, const std::vector<int>& mobile, double cutoff, bool coulomb, double dsf_alpha = 0.2);

}  // namespace caps
