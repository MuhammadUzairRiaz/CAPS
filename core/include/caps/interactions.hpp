// CAPS interactions and checks (design/boards/Interactions): hydrogen bonds, close contacts and clashes of a frame, and
// the problems a structure has before a simulation, each with a fix.
//
//  Hydrogen bonds  donor D (N, O or F) bonded to H, acceptor A (N, O or F) in another molecule or four or more bonds
//                  away: R(D···A) ≤ 3.5 Å and angle H–D···A ≤ 30° (Luzar & Chandler, Nature 379, 55 (1996)).
//  Contacts        non-bonded pairs (not 1-2 or 1-3) closer than the sum of van der Waals radii − 0.4 Å.
//  Clashes         such pairs closer than 0.75 × the sum of the radii.
//  Checks          clashes; molecules cut by the cell edge (wrap); atoms over their valence or lacking hydrogens; the
//                  net charge.
#pragma once
#include <string>
#include <vector>

#include "caps/system.hpp"

namespace caps {

struct InteractionOptions {
  double hb_distance = 3.5;      // Å, D···A
  double hb_angle = 30.0;        // degrees, H–D···A
  double contact_margin = 0.4;   // Å below the vdW sum
  double clash_factor = 0.75;    // × the vdW sum
  size_t max_pairs = 20000;      // listed pairs of each kind (the counts are complete)
};

struct HBond { uint32_t donor, hydrogen, acceptor; double distance, angle; };
struct Contact { uint32_t i, j; double distance, limit; };

struct Issue {
  std::string level;    // error, warn, ok
  std::string title, detail;
  std::string fix;      // "push_apart", "wrap", "add_h", "" none
  std::vector<uint32_t> atoms;
};

struct InteractionReport {
  std::vector<HBond> hbonds;
  std::vector<Contact> contacts, clashes;
  size_t n_hbonds = 0, n_contacts = 0, n_clashes = 0;
  int molecules = 0;
  double net_charge = 0;
  std::vector<Issue> issues;
};

InteractionReport find_interactions(const System& s, const InteractionOptions& o = {});

// Molecules the cell edge cuts (atoms outside the cell): what Wrap would fold.
int molecules_across_edge(const System& s);

}  // namespace caps
