// CAPS molecule builder: SMILES to a 3D structure.
//
//  parse_smiles   OpenSMILES: organic subset and bracket atoms (isotope, chirality @/@@, hydrogen count, charge,
//                 class), bonds - = # $ : / \, branches, ring closures (digits and %nn), dot-disconnected parts.
//                 Implicit hydrogens by the organic-subset valences; they are added as explicit atoms.
//  embed          Distance-bounds embedding (Crippen & Havel; Blaney & Dixon): 1-2 and 1-3 distances from bond
//                 lengths and hybridisation angles, 1-4 ranges from cis and trans, contact lower bounds, with
//                 tetrahedral chirality, sp2 planarity and double-bond E/Z terms; minimised from random
//                 coordinates in four dimensions, then squeezed into three (as in ETKDG, Riniker & Landrum,
//                 J. Chem. Inf. Model. 55, 2562 (2015), without its experimental torsion preferences).
//  build_molecule Several embeddings, each minimised with a force field from CAPS Field when one is given (GAFF2 by
//                 default in the Studio, or UFF for any element); conformers ranked by energy.
#pragma once
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

#include <memory>

#include "caps/field.hpp"
#include "caps/system.hpp"

namespace caps {

struct SmilesError : std::runtime_error {
  size_t position;
  SmilesError(const std::string& what, size_t pos) : std::runtime_error(what), position(pos) {}
};

struct MolAtom {
  int element = 6;
  bool aromatic = false;
  bool bracket = false;
  int charge = 0;
  int isotope = 0;
  int hcount = -1;          // hydrogens written in a bracket; -1 = implicit (organic subset)
  int chiral = 0;           // 0 none, 1 @ (anticlockwise), 2 @@ (clockwise)
  int map = 0;              // atom class [C:1]
  std::vector<int> order;   // neighbours in the order written (for chirality); -2 marks a bracket's own hydrogen
};

struct MolBond {
  int a = 0, b = 0;
  int order = 1;            // 1, 2, 3, 4 aromatic (mol2 code), quadruple bonds are kept as 3
  int dir = 0;              // +1 '/', -1 '\' as written from a to b; 0 none
};

struct MolGraph {
  std::string smiles;
  std::vector<MolAtom> atoms;
  std::vector<MolBond> bonds;
  int heavy = 0;            // atoms written in the SMILES (hydrogens added later come after them)
  int parts = 1;            // dot-disconnected components
};

MolGraph parse_smiles(const std::string& smiles);
// Adds the implicit and bracket hydrogens as atoms bonded to their parent; the chirality neighbour lists follow.
void add_hydrogens(MolGraph& g);
// Valence problems (e.g. "atom 3 (C) has 5 bonds"); empty when the structure is chemically sensible.
std::vector<std::string> valence_problems(const MolGraph& g);

struct MolInfo {
  std::string formula;      // Hill order
  double mass = 0;          // g/mol
  int atoms = 0, heavy = 0, bonds = 0, rings = 0;
  int stereocentres = 0;    // specified tetrahedral centres
  int stereo_bonds = 0;     // specified E/Z double bonds
  int charge = 0;
  std::vector<std::string> problems;
};
MolInfo molecule_info(const MolGraph& g);   // g with hydrogens added

struct EmbedOptions {
  uint64_t seed = 1;
  int attempts = 30;        // random starts before giving up
};
// 3D coordinates for g (hydrogens added) that respect its bond lengths, angles, chirality and E/Z. Throws
// std::runtime_error when no embedding meets the constraints.
std::vector<Vec3> embed(const MolGraph& g, const EmbedOptions& o = {});
// Signed volumes of the specified tetrahedral centres: +1 when the geometry matches the SMILES, -1 when inverted.
std::vector<int> chirality_check(const MolGraph& g, const std::vector<Vec3>& pos);

// 2D drawing coordinates (z = 0, bond length 1, long axis horizontal) for the atoms written in the SMILES (the first
// g.heavy; added hydrogens are not drawn): ring polygons, fused rings sharing edges, zig-zag chains, E/Z as written.
std::vector<Vec3> depict(const MolGraph& g, uint64_t seed = 1);
// SMILES for the written atoms of g (as parse_smiles leaves it): depth-first, ring closures numbered from 1, bracket
// atoms where needed, tetrahedral chirality and E/Z carried over. parse_smiles(write_smiles(g)) is the same molecule.
std::string write_smiles(const MolGraph& g);

struct BuildOptions {
  int conformers = 1;
  uint64_t seed = 1;
  std::string forcefield;   // caps-forcefield JSON with typing rules, or "uff"; empty: embedding only
  std::string charges = "gasteiger";
  double ftol = 0.05;       // kcal/mol/Å
  // after each minimisation, a rotor search: every rotatable bond (acyclic single, a heavy neighbour on both ends) tried
  // at its three staggered positions with the force field's own torsions and contacts, greedily over two passes, then
  // minimised again; kept when lower. A force-field alternative to ETKDG's CSD torsion table, which CAPS does not carry.
  bool rotor_search = false;
  bool implicit_hydrogens = true;   // false: the SMILES as written, heavy atoms only (united-atom models)
};

struct Conformer {
  std::vector<Vec3> pos;
  double energy = 0;        // kcal/mol (force field), or the embedding error without one
  bool minimised = false;
};

struct BuildResult {
  MolGraph graph;
  MolInfo info;
  System system;            // the lowest conformer: atoms, bonds with orders, one molecule, no cell
  std::vector<Conformer> conformers;   // lowest first
  std::string method;
  std::vector<std::string> notes;
};

BuildResult build_molecule(const std::string& smiles, const BuildOptions& o = {});
// A force field for a molecule (hydrogens added) from a caps-forcefield JSON with typing rules, or UFF when path is
// "uff" (then pos, an embedding, sets the axial pairs of five-coordinate centres); null (with the reason in notes)
// when it cannot type or parameterise every atom.
std::shared_ptr<const ForceField> molecule_forcefield(const MolGraph& g, const std::string& path, const std::string& charges,
                                                      std::vector<std::string>& notes, std::string* name = nullptr,
                                                      const std::vector<Vec3>* pos = nullptr);
// Minimises pos with ff (no cell); false (and why) when it fails or inverts a specified centre.
bool minimise_molecule(const MolGraph& g, const std::shared_ptr<const ForceField>& ff, double ftol, std::vector<Vec3>& pos,
                       double* energy = nullptr, std::string* why = nullptr);
// The graph as a System (element, bonds with orders, molecule 1, names like C1 H9) at the given positions.
System molecule_system(const MolGraph& g, const std::vector<Vec3>& pos);

}  // namespace caps
