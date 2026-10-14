// Martini 3 molecules: small molecules mapped from all-atom structures (data/martini/martini3-small-molecules.json:
// vermouth-martinize's martini3001 blocks and mappings, bench/ff/convert_vermouth_martini3_small.py), molecules given as
// GROMACS terms (data/martini/martini3-molecules.json), and the topology as a GROMACS .itp.
#pragma once
#include <cstdint>
#include <functional>
#include <string>
#include <vector>

#include "caps/json.hpp"
#include "caps/system.hpp"

namespace caps {

// The topology as a GROMACS .itp (stiff bonds from constraints written back as constraints)
std::string martini3_itp(const System& beads, double constraint_kj = 1e6);

// Martini 3 small molecules: every molecule of an all-atom structure whose heavy atoms and bonds (and hydrogen counts,
// when it has hydrogens) are one of the model's residues, whatever its atom names, becomes that molecule's beads
// (mass-weighted centres of the mapped atoms, or geometric centres with `geometric`; virtual sites placed) with its
// explicit topology. Molecules that match none are listed in `unmatched` (formulas) and left out.
System martini3_small_molecules(const System& aa, const std::string& data_path, double constraint_kj = 1e6,
                                std::vector<std::string>* unmatched = nullptr, std::vector<std::string>* notes = nullptr, bool geometric = false);

// An all-atom structure in Martini 3: its molecules as small molecules (beads at the geometric centres of their mapped
// atoms, hydrogens included: the mapping Martini 3's small-molecule parameters reproduce, bench/ff/check_martini3_small.py);
// water and ions are left out (noted). Molecules that match no small molecule are refused (std::invalid_argument) with
// refuse_unmatched, as typing must not drop part of a structure; otherwise left out and listed in the notes.
System martini3_all_atom(const System& aa, const std::string& small_path, double constraint_kj = 1e6,
                         std::vector<std::string>* notes = nullptr, bool geometric = true, bool refuse_unmatched = true);

// A molecule given as GROMACS terms (data/martini/martini3-molecules.json: atoms with type, name, charge and optional mass;
// bonds, constraints, angles, dihedrals, exclusions and virtual sites by 0-based index, nm, kJ/mol, degrees, the .itp
// functions): beads named by their types, placed by the bead builder at the bond lengths, with the explicit topology.
// Constraints become stiff bonds of constraint_kj (kJ/mol/nm²); type_mass gives each type's mass for mass-weighted sites.
System gromacs_molecule(const Json& mol, const std::function<double(const std::string&)>& type_mass, double constraint_kj, uint64_t seed = 1);

}  // namespace caps
