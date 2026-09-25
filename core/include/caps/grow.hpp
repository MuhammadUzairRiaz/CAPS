// CAPS Grow: all-atom polymer chains built from internal coordinates and grown inside a periodic cell.
#pragma once
#include <array>
#include <cstdint>
#include <functional>
#include <stdexcept>
#include <string>
#include <vector>

#include "caps/system.hpp"

namespace caps {

enum class Tacticity { Atactic, Isotactic, Syndiotactic };

Tacticity tacticity_from_string(const std::string& s);   // "atactic" | "isotactic" | "syndiotactic"
const char* to_string(Tacticity t);

struct GrowOptions {
  std::string polymer = "polystyrene";   // the only monomer in this slice
  int chains = 10;
  int dp = 8;                            // degree of polymerisation (monomer units per chain)
  Tacticity tacticity = Tacticity::Atactic;
  uint64_t seed = 1;
  double box = 0.0;                      // cubic cell edge, Å; 0 = derive from density
  double density = 0.0;                  // target g/cm^3 when box == 0
  bool curve = true;                     // allow gauche backbone torsions (false: near-trans only)
  int trials = 120;                      // torsion / ring trials per growth step (the roomiest is kept)
  double comfortable = 0.0;              // stop trying once a trial clears every limit by this much, Å
  bool escalate = false;                 // back up further (and try harder) where a step keeps failing
  int max_backtracks = 0;                // per chain start; 0 = 400, or 40 × dp when escalating
  int max_restarts = 400;                // new start points per chain before giving up
  double accept = -0.05;                 // worst allowed (distance − limit), Å
  // Called once per growth round with (chains finished, chains, restarts so far); return false to cancel.
  std::function<bool(int, int, int)> progress;
  double contact_scale = 1.0;            // scales the contact limits (C–C 3.0, C–H 2.45, H–H 2.0 Å); < 1 needs Relax afterwards
  // Films and interfaces (the general polymer builder, grow_chains): an orthorhombic cell, a height range for the
  // chains, and fixed atoms they must avoid.
  std::array<double, 3> cell{0, 0, 0};   // edges x, y, z (Å), used when all three are > 0 (box and density are then ignored)
  double z_lo = 0, z_hi = 0;             // when z_hi > z_lo, every chain atom stays between these heights
  const System* substrate = nullptr;     // fixed atoms (a slab) the chains avoid; they come first in the result
  // a spherical region (blend droplets): when sphere_radius > 0 every chain atom stays inside the sphere about
  // sphere_centre, or outside it with sphere_outside (minimum image in the periodic cell)
  double sphere_radius = 0;
  std::array<double, 3> sphere_centre{0, 0, 0};
  bool sphere_outside = false;
  // grow_chains: when a chain cannot be placed, try again at contact scales 0.85, 0.75, 0.7, 0.6 of the full limits
  // (quaternary backbones such as polyisobutylene and methacrylates, dense films); Relax with push-off afterwards
  bool auto_scale = false;
};

struct GrowReport {
  int chains_placed = 0;
  int restarts = 0;
  int backtracks = 0;
  double worst_margin = 0.0;             // smallest (distance − limit) over accepted non-bonded pairs, Å
  double box = 0.0;
  double density = 0.0;
  std::vector<std::string> notes;
};

struct GrowError : std::runtime_error {
  using std::runtime_error::runtime_error;
};

// Mass of one chain (g/mol) for a monomer and degree of polymerisation, H end caps included.
double chain_mass(const std::string& polymer, int dp);

// Cell edge (Å) that gives the target density for the requested chains.
double box_for_density(const GrowOptions& o);

// Grow all chains into a periodic cubic cell. Positions are unwrapped (each chain continuous).
// Atoms carry GAFF-style names (c3, ca, hc, ha), types 1–4, molecule ids and Gasteiger–Marsili charges.
// Throws GrowError when a chain cannot be placed within the limits.
System grow(const GrowOptions& o, GrowReport* report = nullptr);

// Gasteiger–Marsili partial charges (H, C, N, O, F, Cl, Br, I, sp³ S; hybridisation from the perceived bonds). One
// charge per atom; throws std::invalid_argument for other elements.
std::vector<double> gasteiger_ch(const System& s, const std::vector<char>& aromatic, int iterations = 6);

}  // namespace caps
