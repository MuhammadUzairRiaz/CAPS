// Labels for the 3D view: one text per atom or per bond, of a chosen property (design/boards/Appearance, the Materials
// Studio label dialog's objects and properties). Every value comes from the structure, its perceived chemistry
// (typing.hpp perceive: formal charges, bond orders, rings, aromaticity) or tabulated element data (elements.hpp).
//
//  Derived quantities, and how:
//   hybridisation      from bond orders: a triple bond or two double bonds sp, one double or an aromatic bond sp²,
//                      otherwise sp³ (boron with three single bonds sp²); none for H, noble gases and metals
//   oxidation state    the IUPAC bond rule: each bond to a more electronegative partner (Pauling) +order, to a less
//                      electronegative one −order, to the same element 0, plus the formal charge (Kekulé orders for
//                      aromatic bonds; implicit united-atom hydrogens count as bonds to H)
//   unpaired electrons the valence left open: the element's usual valence at its formal charge minus the bond-order
//                      sum (a radical); main-group elements only
#pragma once
#include <string>
#include <vector>

#include "caps/system.hpp"

namespace caps {

struct ForceField;

struct LabelKind {
  std::string id, title, group;
};
const std::vector<LabelKind>& atom_label_kinds();
const std::vector<LabelKind>& bond_label_kinds();

// One text per atom ("" where the property does not apply). `types`: the force-field type of each atom from an
// assignment (else the structure's own type labels).
std::vector<std::string> atom_labels(const System& s, const std::string& kind, const std::vector<std::string>* types = nullptr);

struct BondLabel {
  uint32_t i = 0, j = 0;
  std::string text;
  bool crossing = false;   // the two atoms as stored lie in different periodic images (the bond crosses the cell)
};
// One text per bond of s.bonds. `ff`: the assignment's terms, for the force-field kinds (style, r0, k, energy).
std::vector<BondLabel> bond_labels(const System& s, const std::string& kind, const std::vector<std::string>* types = nullptr, const ForceField* ff = nullptr);

}  // namespace caps
