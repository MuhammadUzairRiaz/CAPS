// CAPS structure editing (design/boards/Main builder tools, ElementPicker, SelectionStereo, AddHydrogens): the small
// operations the Studio's builder tools make, each keeping atoms, bonds and types consistent.
//
//  add_atom        a new atom bonded to another at the covalent-radius bond length, in the ideal direction for the
//                  parent's hybridisation (sp³ tetrahedral, sp² trigonal, sp linear) that is farthest from its bonds.
//  add_hydrogens   the hydrogens each atom lacks for its valence (formal charge counted), placed the same way.
//  invert_centre   the configuration of a tetrahedral centre changed by swapping its two smallest acyclic branches
//                  (a 180° rotation about the bisector of their bonds).
//  tacticity       vinyl-polymer stereocentres along each backbone (a CH with one side group), their relative
//                  configuration, m/r dyads and mm/mr/rr triads; set_tacticity inverts centres to make a chain iso- or
//                  syndiotactic.
//  selections      by SMARTS (the pattern's first atom), element, type, charge range, distance from a set, or grown
//                  along bonds.
#pragma once
#include <stdexcept>
#include <string>
#include <vector>

#include "caps/system.hpp"

namespace caps {

struct EditError : std::runtime_error {
  using std::runtime_error::runtime_error;
};

// The usual valence of an element with a formal charge (N⁺ 4, O⁻ 1 …); 0 for elements without one.
int default_valence(int z, int charge = 0);

// The atom's type follows its element (a type labelled with the symbol is added when missing).
void set_element(System& s, uint32_t atom, int z);
// Returns the new atom's index. bonded_to < 0: an isolated atom at `at`.
// geometry: 0 from the parent's neighbours (sp³ unless it has a double or triple bond), 3 sp³, 2 sp², 1 sp.
uint32_t add_atom(System& s, int bonded_to, int z, int order = 1, int geometry = 0, Vec3 at = {0, 0, 0});
void add_bond(System& s, uint32_t i, uint32_t j, int order = 1);
bool remove_bond(System& s, uint32_t i, uint32_t j);
// Removes the flagged atoms and their bonds; the rest keep their order.
void delete_atoms(System& s, const std::vector<char>& remove);
// Adds the hydrogens the flagged atoms (empty: all) lack; returns how many.
int add_hydrogens(System& s, const std::vector<char>& atoms = {});
// What add_hydrogens would do, by kind of atom ("aromatic C with 2 C neighbours": atoms, H to add).
struct HydrogenPlanRow {
  std::string label;
  int atoms = 0, hydrogens = 0;
};
std::vector<HydrogenPlanRow> hydrogen_plan(const System& s, const std::vector<char>& atoms = {});
void invert_centre(System& s, uint32_t centre);

struct TacticityChain {
  std::vector<uint32_t> centres;   // stereocentres in chain order
  std::vector<int> sign;           // relative configuration (+1 / −1) of each
  std::string dyads;               // "m" or "r" between consecutive centres
};
struct TacticityReport {
  std::vector<TacticityChain> chains;
  int m = 0, r = 0, mm = 0, mr = 0, rr = 0;
  std::string label;               // isotactic, syndiotactic, atactic, or "" without centres
};
TacticityReport tacticity(const System& s);
// Makes every chain isotactic (iso = true) or syndiotactic; returns the centres inverted.
int set_tacticity(System& s, bool iso);

// Attaches a fragment written as SMILES with * attachment points ("*C(=O)O*"): its `which`-th attachment point goes on
// the target (replacing one of the target's hydrogens when it has one, else in its free direction); the fragment is
// rolled about the new bond to keep clear of the structure; other attachment points become hydrogens. Returns the
// fragment's atoms in the structure.
std::vector<uint32_t> attach_fragment(System& s, uint32_t target, const std::string& smiles, int which = 0, bool replace_h = true);
// Attachment points of a fragment's SMILES: for each *, the index (among the written atoms) of the atom it hangs on.
std::vector<int> fragment_attach_atoms(const std::string& smiles);

// Fuses a benzene ring onto the bond i–j: a hydrogen of each atom on the same side goes, four carbons (with a hydrogen
// each) complete a regular hexagon on that side, in the plane of the bond and the two hydrogens, and the ring's six
// bonds (i–j included) become aromatic. Returns the new atoms' indices (after the two hydrogens are removed).
std::vector<uint32_t> fuse_benzene(System& s, uint32_t i, uint32_t j);

// Protonation of amino-acid residues at a pH (the Add hydrogens tool, before add_hydrogens): formal charges on the
// titratable atoms, found by residue name and chemistry (not atom names): Asp and Glu carboxylates −1 above pKa 3.9 and
// 4.3, His +1 below 6.0, Cys (not in a disulfide) −1 above 8.3, Tyr −1 above 10.1, Lys +1 below 10.5, Arg +1 below
// 12.5 (the peptide builder's values), the N-terminus +1 below 8.0 and the C-terminus −1 above 3.1 — model pKa
// values, no shifts from the environment. Needs bond orders (orders_from_geometry for heavy-atom files). Returns the
// charged sites; `notes` gets what was set.
int protonate_residues(System& s, double ph, std::vector<std::string>* notes = nullptr);
// protonate_residues, then add_hydrogens, then histidine's ring N–H that bond orders cannot tell (aromatic imidazole):
// both ring nitrogens of His+, the one farther from the backbone (HIE) of neutral His. Returns the hydrogens added.
int add_hydrogens_at_ph(System& s, double ph, const std::vector<char>& atoms = {}, std::vector<std::string>* notes = nullptr);

// Minimises the flagged atoms (the rest held) with UFF; push-off first for overlaps.
void clean_up(System& s, const std::vector<char>& atoms = {}, double ftol = 0.5);

std::vector<char> select_smarts(const System& s, const std::string& pattern);
std::vector<char> select_element(const System& s, const std::string& symbols);   // "C", "C N O"
std::vector<char> select_type(const System& s, const std::string& label);         // type label or number
std::vector<char> select_charge(const System& s, double lo, double hi);
std::vector<char> select_within(const System& s, const std::vector<char>& from, double distance);
std::vector<char> select_grow(const System& s, const std::vector<char>& from, int steps = 1);

}  // namespace caps
