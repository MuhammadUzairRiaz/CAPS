// CAPS Grow: all-atom polymer chains built from internal coordinates and grown inside a periodic cell.
#pragma once
#include <array>
#include <cstdint>
#include <cmath>
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
  int threads = 0;                       // grow_chains: threads scoring the trials of a step (0: up to 8; 1: serial)
  int lookahead = 1;                     // grow_chains: bonds ahead that need room (1: the next head; 2–4 also points further along)
  // Called once per growth round with (chains finished, chains, restarts so far); return false to cancel.
  std::function<bool(int, int, int)> progress;
  double contact_scale = 1.0;            // scales the contact limits (C–C 3.0, C–H 2.45, H–H 2.0 Å); < 1 needs Relax afterwards
  // Films and interfaces (the general polymer builder, grow_chains): an orthorhombic cell, a height range for the
  // chains, and fixed atoms they must avoid.
  std::array<double, 3> cell{0, 0, 0};   // edges x, y, z (Å), used when all three are > 0 (box and density are then ignored)
  double z_lo = 0, z_hi = 0;             // when z_hi > z_lo, every chain atom stays between these heights
  const System* substrate = nullptr;     // fixed atoms (a slab) the chains avoid; they come first in the result
  // grafted chains (brushes): chain k < anchors.size() starts bonded to substrate atom anchors[k].atom, its first unit's
  // head leaving along anchors[k].dir (the site's outward direction); that head takes no end cap. The atoms within three
  // bonds of the anchor count as bonded neighbours, not contacts. The substrate must have the site's own valence free
  // (a silanol's H already taken off).
  struct Anchor { uint32_t atom = 0; std::array<double, 3> dir{0, 0, 1}; };
  std::vector<Anchor> anchors;
  // a spherical region (blend droplets): when sphere_radius > 0 every chain atom stays inside the sphere about
  // sphere_centre, or outside it with sphere_outside (minimum image in the periodic cell)
  double sphere_radius = 0;
  std::array<double, 3> sphere_centre{0, 0, 0};
  bool sphere_outside = false;
  // Region shapes (grow_chains, when no cell is given; the edges come from the density):
  //   slab_thickness > 0: a film of that thickness between slab_vacuum / 2 of vacuum above and below (square in x, y)
  //   cylinder_radius > 0: a cylinder along z through the cell centre — chains inside it (a pore), or outside it with
  //   cylinder_outside (around a fibre); a cubic cell unless cylinder_length > 0
  double slab_thickness = 0, slab_vacuum = 0;
  double cylinder_radius = 0, cylinder_length = 0;
  bool cylinder_outside = false;
  // Growth method (grow_chains): 0 best of the trials by contact margin (roomy but stretched chains); 1 Rosenbluth: a
  // trial drawn with probability ∝ exp(−E/kT), E soft-sphere overlap + torsion energy (after Theodorou & Suter 1985);
  // 2 Rosenbluth with UFF Lennard-Jones between atoms (by element, cut at 6 Å) + the same torsion energy
  // (configurational-bias growth; there is no Monte Carlo acceptance step). The torsion energy about sp3–sp3 bonds is
  // Jorgensen's butane potential (J. Am. Chem. Soc. 106, 6638 (1984): trans 0, gauche 0.86, barrier 4.6 kcal/mol).
  // Trials under the contact limits still fail. The chains' Rosenbluth weights are reported.
  int method = 0;
  double method_temperature = 450;       // K, of the Boltzmann factors
  // Orientation (grow_chains; design/boards/Grow "Orientation"): an aligning field on each unit's backbone chord (the
  // previous unit's tail to this one's), energy −s P₂(cos θ) in kT with θ its angle to orient_axis — drawn chains, fibres,
  // an oriented amorphous start. Rosenbluth methods weight each trial by exp(s P₂) besides exp(−E/kT); the best-of-k
  // method adds 0.25 s P₂ Å to a trial's contact margin (among trials within the limits). s = 0 (or no axis): isotropic.
  std::array<double, 3> orient_axis{0, 0, 0};
  double orient_strength = 0;
  // grow_chains: when a chain cannot be placed, try again at contact scales 0.85, 0.75, 0.7, 0.6 of the full limits
  // (quaternary backbones such as polyisobutylene and methacrylates, dense films); Relax with push-off afterwards
  bool auto_scale = false;
  // A live view while growing (grow_chains; design/boards/GrowAllAtom): about every snapshot_seconds, the chains so far
  // (atoms and bonds, molecule = chain, no end caps) and where the growth stands.
  struct Live {
    int chains_done = 0, chains = 0;
    long units = 0, units_total = 0;
    int restarts = 0;
    double worst_margin = 0;   // smallest (distance − limit) over accepted pairs so far, Å
    double density = 0;        // of the atoms placed so far, g/cm³
    double acceptance = 1;     // growth steps placed at the first attempt / all attempts (a backtrack counts as a failed one)
    double ln_w = std::nan(""); // Rosenbluth methods: mean ln W of the chains finished so far
  };
  std::function<void(const System&, const Live&)> snapshot;
  double snapshot_seconds = 0.25;
};

struct GrowReport {
  int chains_placed = 0;
  int restarts = 0;
  int backtracks = 0;
  double worst_margin = 0.0;             // smallest (distance − limit) over accepted non-bonded pairs, Å
  double box = 0.0;
  double density = 0.0;
  double ln_rosenbluth = 0.0;            // Rosenbluth methods: mean over chains of ln W (W = Π_steps Σ w / k)
  double orientation = 0.0;              // ⟨P₂⟩ of the units' backbone chords against orient_axis (z when none is given)
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
