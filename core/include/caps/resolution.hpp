// CAPS model resolution (design/boards/ModelResolution): the same structure at three resolutions, mass conserved.
//
//  united_atom    every hydrogen bonded to carbon folded into its carbon (CH3, CH2, CH sites; masses and charges summed;
//                 polar hydrogens kept, as TraPPE-UA does)
//  coarse_grain   each molecule's backbone cut into beads of `per_bead` backbone atoms; side groups and hydrogens join
//                 the bead of the backbone atom they hang on; a bead sits at its atoms' centre of mass and carries their
//                 mass (molecules without a backbone become one bead each)
//
// Both return new structures (the original is untouched). Going back up: united-atom → all-atom is add_hydrogens then
// a relax; coarse-grained → all-atom is backmap onto the all-atom structure the beads came from.
#pragma once
#include <functional>
#include <map>
#include <string>
#include <vector>

#include "caps/system.hpp"

namespace caps {

struct ResolutionReport {
  int sites = 0;             // atoms, united sites or beads
  int hydrogens = 0;         // explicit hydrogens left
  double mass = 0;           // g/mol, the sum of the sites' masses (IUPAC standard atomic weights underneath)
  std::vector<int> site_of;  // for each original atom, its site in the result
  std::vector<std::string> notes;
};

// All-atom mass and counts of a structure (by element, standard atomic weights).
ResolutionReport all_atom_summary(const System& s);
// hosts: the elements whose hydrogens fold in (carbon by default; oxygen for a one-site water such as mW, site "OH2")
System united_atom(const System& s, ResolutionReport* rep = nullptr, const std::vector<int>& hosts = {6});
System coarse_grain(const System& s, int per_bead = 5, ResolutionReport* rep = nullptr);

// Backmapping (multiscale equilibration: coarse-grain here, run the beads elsewhere, come back): the beads of
// coarse_grain(all_atom, per_bead) have moved to `beads` (same count and order, e.g. read from a LAMMPS dump). Each
// bead's atoms are carried along rigidly: translated with the bead and turned by the rotation that best maps the bead
// and its bonded beads from their old places to their new ones (Horn's quaternion fit; a chain end turns with its one
// neighbour, a lone bead only moves). Bonds between beads come out stretched or squeezed: relax afterwards (push-off).
// The cell is the beads'. Throws when the bead counts differ.
struct BackmapReport {
  int beads = 0, atoms = 0;
  double rms_turn = 0;       // degrees, the beads' rotations
  double worst_bond = 0;     // Å, the longest bond between two beads' atoms before a relax
  std::vector<std::string> notes;
};
System backmap(const System& all_atom, const System& beads, int per_bead = 5, BackmapReport* rep = nullptr);

// Coarse-grained molecules written as bead SMILES, the template notation of the MARTINI and SDK sources: [TYPE] or
// [TYPE±q] beads (Qa-1, Qd+1), branches, ring closures 1–9, %nn and %nnn; '.' separates molecules. Throws
// std::invalid_argument on a malformed string.
struct BeadMolecule {
  std::vector<std::string> type;
  std::vector<double> charge;
  std::vector<std::pair<int, int>> bonds;
};
BeadMolecule parse_bead_smiles(const std::string& text);

// One site per bead, named by its type (element 0: a bead, not an element), its charge from the suffix, bonded as
// written; each type's mass from `mass` (else default_mass). Coordinates: breadth-first placement at the bond lengths
// (`bond_length`, else default_bond), then a spring model (bonds at their lengths, beads two bonds apart at least
// 1.4 bond lengths apart, others 0.9 of the longest bond) — a start for a relax with the force field.
struct BeadBuildOptions {
  std::function<double(const std::string&, const std::string&)> bond_length;   // Å, 0 = unknown
  std::function<double(const std::string&)> mass;                             // g/mol, 0 = unknown
  double default_bond = 4.7, default_mass = 72.0;
  uint64_t seed = 1;
};
System build_beads(const std::string& text, const BeadBuildOptions& o = {});

// All-atom → beads by fragment rules (a coarse-grained force field's typing file, "beads" and "bead_groups"). Every
// heavy atom of a molecule is covered by exactly one fragment (an exact cover: at chain ends the fragments are chosen so
// the chain tiles, as SDK's CT / CT2 ends with CM between); hydrogens join their heavy atom; each bead sits at its
// fragment's centre of mass, named by its type (element 0), with the fragment's mass and summed charge. A rule with
// several types (SDK's two ester beads) is resolved per molecule by the bond lengths the force field gives to the bead
// they share (`bond_length`, Å). Groups gather whole small molecules (SDK's W: three waters, nearest first). Beads are
// bonded where their fragments are. Molecules that cannot be covered keep their atoms (reported, left for typing to
// flag).
struct BeadRule {
  std::vector<std::string> types;   // one, or alternatives told apart by bond length
  std::string smarts, description;  // heavy atoms of the fragment
};
struct BeadGroup {
  std::string type, molecule, description;   // molecule: SMARTS for the one heavy atom of each small molecule
  int count = 1;
};
struct BeadMapReport {
  int beads = 0, molecules_mapped = 0;
  std::map<std::string, int> by_type;
  std::vector<uint32_t> uncovered;   // atoms of the input left as they were
  std::vector<std::string> notes;
};
System map_to_beads(const System& s, const std::vector<BeadRule>& rules, const std::vector<BeadGroup>& groups,
                    const std::function<double(const std::string&, const std::string&)>& bond_length = {}, BeadMapReport* rep = nullptr);

}  // namespace caps
