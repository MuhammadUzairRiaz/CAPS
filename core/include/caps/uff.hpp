// CAPS UFF: the Universal Force Field (Rappé, Casewit, Colwell, Goddard & Skiff, J. Am. Chem. Soc. 114, 10024 (1992))
// for every element from H to Lr, built from the atoms' elements and bonds alone, for geometry clean-up and
// minimisation of any molecule (polymers with Si, P, S, halogens, metals ...).
//
//  Typing     element + hybridisation (1 linear, 2 trigonal, R resonant, 3 tetrahedral, 4 square planar,
//             5 trigonal bipyramidal, 6 octahedral) + oxidation state, from the perceived bond orders, lone pairs and
//             conjugation; metals take the geometry of their coordination number. When the exact label is not in the
//             table, the closest one of the same element is used (and reported).
//  Bonds      k (r − r0)², r0 = ri + rj + rBO − rEN, k = 332.06 Zi Zj / r0³ (UFF's ½ k is folded in).
//  Angles     K [C0 + C1 cos θ + C2 cos 2θ] at general centres; K (1 − cos nθ) / n² at linear (n = 1, as 1 + cos θ),
//             trigonal (n = 3) and square-planar / octahedral (n = 4) centres; three- and four-ring sp2 special cases.
//  Torsions   ½ V [1 − cos(n φ0) cos(n φ)] shared over the torsions about a bond, with UFF's sp3–sp3, sp2–sp2,
//             sp2–sp3 and group-16 rules.
//  Inversions K [C0 + C1 cos ω + C2 cos 2ω] at sp2 C, N, O (K = 6, 50 for C bonded to sp2 O) and pyramidal P, As, Sb, Bi.
//  van der Waals  D [(x/r)¹² − 2 (x/r)⁶], geometric means of x and D; 1-2 and 1-3 pairs excluded, 1-4 at full strength.
//  No charges (UFF's clean-up use, as in RDKit and Open Babel).
//
// The functional forms and special cases follow RDKit's implementation of UFF, whose parameter table CAPS carries.
#pragma once
#include <string>
#include <vector>

#include "caps/ffdef.hpp"
#include "caps/field.hpp"
#include "caps/system.hpp"

namespace caps {

struct UffOptions {
  bool keep_charges = false;          // use the system's charges (with Coulomb) instead of UFF's neutral atoms
  std::vector<std::string> labels;    // per atom: a UFF label set by hand ("" or missing: the typer's)
};

// UFF for every atom of `s`. Throws FieldError only for atoms with no element or an element past Lr.
ForceField assign_uff(const System& s, const UffOptions& o = {});

// The UFF label for each atom of `s` (as assign_uff types them) and the reason for each.
std::vector<std::string> uff_types(const System& s, std::vector<std::string>* why = nullptr);

// The force field used when none is given: CAPS's built-in GAFF for saturated and aromatic hydrocarbons, UFF for
// anything else (double bonds, heteroatoms, metals).
ForceField default_forcefield(const System& s);

// True for the force-field names that mean UFF: "uff" (any case), or a path whose file name is uff / uff.json.
bool is_uff(const std::string& name_or_path);

// UFF as a force-field description (every label as a type, with its element, mass and parameters in the
// description; no typing rules: assign_uff types), for pages that list a force field's types.
FFDef uff_definition();

// Number of UFF atom labels in the table (127) and whether a label is one of them.
int uff_label_count();
bool uff_has_label(const std::string& label);

}  // namespace caps
