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

// Exact internal coordinates by moving one side of the structure (the side of the last atom, which must not be in a
// ring with the bond that moves): the bond i–j to r Å, the angle i–j–k to θ degrees, the dihedral i–j–k–l to φ degrees.
// Throws when the chosen bond is in a ring (the two sides are one).
void set_bond_length(System& s, uint32_t i, uint32_t j, double r);
void set_bond_angle(System& s, uint32_t i, uint32_t j, uint32_t k, double theta_deg);
void set_torsion(System& s, uint32_t i, uint32_t j, uint32_t k, uint32_t l, double phi_deg);
// Rigid rotation of the atoms by `degrees` about `axis` through their centre; reflection of the atoms through the
// plane with normal `normal` through their centre (a mirror image: every stereocentre among them inverts).
void rotate_atoms(System& s, const std::vector<uint32_t>& atoms, const Vec3& axis, double degrees);
// The bonded neighbours of `centre` placed at the ideal directions of a coordination geometry — linear, trigonal
// (planar), tetrahedral, square_planar, trigonal_bipyramidal, square_pyramidal, octahedral — each keeping its bond length
// and carrying its own substituents (a ligand in a chelate ring moves alone). The ideal set is turned to move the
// ligands least; fewer ligands than sites fill the sites nearest them. Returns the largest angle a ligand moved (degrees).
double set_coordination(System& s, uint32_t centre, const std::string& geometry);
// The ideal unit directions of a coordination geometry (empty for an unknown name).
std::vector<Vec3> coordination_directions(const std::string& geometry);
void mirror_atoms(System& s, const std::vector<uint32_t>& atoms, const Vec3& normal);
// The centre made R or S (CIP, from the 3D geometry) by inverting it when it is the other; false when it is not a
// stereocentre.
bool set_configuration(System& s, uint32_t centre, const std::string& rs);

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
// direction: where the attaching atom goes from the target (a site's outward normal) when no hydrogen is replaced.
std::vector<uint32_t> attach_fragment(System& s, uint32_t target, const std::string& smiles, int which = 0, bool replace_h = true, const Vec3* direction = nullptr);

// Silane coupling agents grafted onto surface silanols (silica fillers in rubber: TESPT / Si69 couples silica to the
// sulfur-cured matrix): each chosen Si–O–H loses its H and the silane's silicon bonds to that oxygen (the condensation's
// ethanol goes), the rest of the silane pointing away from the surface. Sites are drawn at random (seed) among the
// silanols, at least min_spacing apart; count > 0 grafts that many, else fraction of the silanols.
struct GraftOptions {
  std::string smiles = "*[Si](OCC)(OCC)CCCSSSSCCC[Si](OCC)(OCC)OCC";   // the silane with * where it bonds (TESPT)
  std::string name = "TESPT";
  int count = 0;
  double fraction = 0.25;
  double min_spacing = 5.0;   // Å between grafted oxygens
  uint64_t seed = 1;
};
struct GraftReport { size_t silanols = 0, grafted = 0, added_atoms = 0; std::vector<std::string> notes; };
// the preset silanes by name: TESPT (Si69), TESPD (Si75), MPTES, APTES, VTES, OCTEO
std::string silane_smiles(const std::string& name);
GraftReport graft_silanes(System& s, const GraftOptions& o);
// Thiolate capping of a metal particle or surface (gold, silver, copper, platinum, palladium): alkanethiolates bound
// through S in three-fold hollow sites of the surface metal atoms (2.45 Å from each of the three; on top of an atom where
// the surface has no hollow), sites at least min_spacing apart (0: √3 × the metal's nearest-neighbour distance, the
// √3 × √3 R30° packing of thiolates on Au(111)), the tail pointing out of the surface. No metal–S bonds are written:
// the S carries the thiolate; relax with a force field that has metal–S terms (INTERFACE) before dynamics.
struct ThiolateOptions {
  std::string smiles = "*SCCCCCC";   // the ligand with * where it binds (hexanethiolate)
  std::string name = "hexanethiolate";
  double fraction = 1.0;             // of the sites that fit at the spacing
  double min_spacing = 0.0;          // Å between S atoms, 0 = √3 · d(M–M)
  uint64_t seed = 1;
};
struct ThiolateReport { size_t surface_atoms = 0, ligands = 0, hollow = 0, on_top = 0, added_atoms = 0; std::vector<std::string> notes; };
// the preset ligands by name: C6 (hexanethiolate), C12 (dodecanethiolate), C18, MPA (3-mercaptopropionate, acid form),
// MUA (11-mercaptoundecanoic acid), MHA (6-mercaptohexanol)
std::string thiolate_smiles(const std::string& name);
ThiolateReport cap_thiolates(System& s, const ThiolateOptions& o);

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

// Hydrogens on phosphorus become hydroxyls (a nucleic-acid strand's 3′ end, capped with H by the grower, becomes a
// 3′-phosphate P–OH): the H turns into O at 1.61 Å and gets its own H. Returns the groups changed.
int hydroxylate_phosphorus(System& s);

// Minimises the flagged atoms (the rest held) with UFF; push-off first for overlaps.
void clean_up(System& s, const std::vector<char>& atoms = {}, double ftol = 0.5);

std::vector<char> select_smarts(const System& s, const std::string& pattern);
std::vector<char> select_element(const System& s, const std::string& symbols);   // "C", "C N O"
std::vector<char> select_type(const System& s, const std::string& label);         // type label or number
std::vector<char> select_charge(const System& s, double lo, double hi);
std::vector<char> select_within(const System& s, const std::vector<char>& from, double distance);
std::vector<char> select_grow(const System& s, const std::vector<char>& from, int steps = 1);

// A written list of edits applied in order (the CLI's `caps edit`), separated by ';' or new lines ('#' after a space
// starts a comment),
// atoms numbered from 1 as the structure stands at that edit (a delete renumbers the atoms after it):
//   element SEL Sym · delete SEL · bond I J [order] · unbond I J · addh [SEL] · attach I SMILES · length I J Å ·
//   angle I J K ° · torsion I J K L ° · invert I · config I R|S · rotate SEL x,y,z ° · mirror SEL nx,ny,nz ·
//   move SEL dx,dy,dz · clean [SEL] · tacticity iso|syndio
//   lattice (lattice.hpp): supercell na nb nc · primitive [tol] · niggli · conventional [tol] · redefine m11 … m33 ·
//   vacuum Å · nanowire u v w radius [repeats] [shape] [vacuum]
// SEL: "3,5-9", "all", "element:C,N", "smarts:PATTERN", "type:LABEL". Returns one line per edit saying what it did;
// throws on the first edit that cannot be done, naming it.
std::vector<std::string> edit_script(System& s, const std::string& script);
// The atoms a SEL names (see edit_script).
std::vector<char> parse_selection(const System& s, const std::string& sel);

}  // namespace caps
