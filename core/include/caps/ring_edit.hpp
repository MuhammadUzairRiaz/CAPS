#pragma once
// Rings, insertions and joins in 3D (the sketcher's ring tool beyond benzene):
//
//  fuse_ring    a ring of 3 to 8 atoms sharing the bond i–j: a hydrogen of each (on the same side) goes, the new atoms
//               complete a regular polygon with side |ij| in the plane of the bond and that side (epoxide on C–C, a
//               lactam or pyridine ring onto a chain or ring);
//  spiro_ring   a ring of 3 to 8 atoms sharing one atom: two of its hydrogens go (or, with one or none, its two open
//               tetrahedral directions), the polygon lies in the plane of those two directions, symmetric about their
//               bisector, its side the mean covalent sum of the ring bonds (spiro siloxanes, spiro lactams);
//  insert_atom  an atom put into the bond i–j (one more S in a polysulfide bridge, O into Si–Si, CH2 into a chain):
//               i–X and X–j at covalent sums, the angle i–X–j tetrahedral; an open chain gives way (the smaller side
//               moves out along the bond), in a ring X sits off the bond on the side away from the ring;
//  join_pieces  a second structure bonded to the first by two leaving atoms (a terminal atom of each, usually a
//               hydrogen): the leaving atoms go, the piece is turned so its bond points back along the first's, the
//               new bond at the covalent sum, and the piece turned about the bond to the torsion asked, or to the one
//               that keeps it farthest from the first structure.
//
// Ring members are written one token per new atom in ring order (from j round to i for a fused ring, from the shared
// atom round for a spiro ring), separated by spaces or commas: an element symbol, optionally "=O" or "=S" for a doubly
// bonded atom outside the ring ("C=O" for the carbonyl of a lactam). Hydrogens fill the new atoms' valences. Covalent
// radii are Cordero et al. (Dalton Trans. 2008, 2832), the double-bond shortening is CAPS's add_atom rule. The start
// geometry is a regular polygon (puckered for 6- to 8-rings so a minimiser can leave the plane); clean = true relaxes
// the new atoms and the atoms they touch with UFF.
#include <cmath>
#include <string>
#include <vector>

#include "caps/system.hpp"

namespace caps {

struct RingSpec {
  int size = 6;                  // ring atoms, 3–8, the shared ones included
  std::string members;           // the new atoms (empty: all carbon)
  bool aromatic = false;         // fused six-rings only: all ring bonds aromatic (i–j too)
  bool hydrogens = true;
  bool clean = true;
};
struct RingEditResult {
  std::vector<uint32_t> ring;    // the ring's atoms in order (shared atoms included), indices after the edit
  std::vector<uint32_t> added;   // every new atom, indices after the edit
  int removed = 0;               // leaving atoms deleted
  std::vector<std::string> notes;
};
RingEditResult fuse_ring(System& s, uint32_t i, uint32_t j, const RingSpec& spec);
RingEditResult spiro_ring(System& s, uint32_t centre, const RingSpec& spec);
RingEditResult insert_atom(System& s, uint32_t i, uint32_t j, int z, bool hydrogens = true, bool clean = false);

struct JoinOptions {
  int order = 1;                 // the new bond
  double torsion = NAN;          // degrees: the dihedral (first-structure neighbour)–pa–pb–(piece neighbour); NaN: the clearest
  bool clean = false;
};
// a's leaving atom la and b's leaving atom lb must each have exactly one bond. b's atoms are appended (types matched by
// label, molecule of a's bonding atom); a's cell is kept.
RingEditResult join_pieces(System& a, const System& b, uint32_t la, uint32_t lb, const JoinOptions& o = {});
// The leaving atom CAPS picks in a piece when none is named: the first hydrogen bonded to the piece's first heavy atom
// (then any terminal atom); throws when the piece has none.
uint32_t default_leaving_atom(const System& b);

}  // namespace caps
