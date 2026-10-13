// The force field of a LAMMPS data file CAPS wrote, rebuilt from its coefficient sections and charges (io_lammps.cpp).
#pragma once
#include <cstdint>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

#include "caps/system.hpp"

namespace caps {

// What read_lammps_data keeps of a data file for its force field: every coefficient section's rows (tokens, comments
// stripped) and style comment, each bonded term's type and atom ids, the bond types of System::bonds, and the
// "# caps-ff" lines CAPS writes under the title (the force field's name, its 1-2 / 1-3 / 1-4 weights, or what the file
// cannot hold).
struct LammpsFileTerms {
  std::map<std::string, std::vector<std::vector<std::string>>> rows;   // "PairIJ Coeffs", "Bond Coeffs", "BondBond Coeffs" ...
  std::map<std::string, std::string> style;                            // section → its style comment ("harmonic", "hybrid")
  std::vector<int> bond_type;                                          // per System::bonds entry
  struct Term { int type; std::vector<int64_t> ids; };
  std::vector<Term> angles, dihedrals, impropers;
  std::vector<std::string> caps_ff;                                    // the text after "# caps-ff"
};

// Sets s.forcefield from the file's coefficients when CAPS wrote them all (the "# caps-ff" record and styles it reads
// back), else s.forcefield_error with the reason. A file with no coefficient sections is a structure: neither is set.
void rebuild_lammps_forcefield(System& s, const LammpsFileTerms& d, const std::unordered_map<int64_t, uint32_t>& index, const std::string& path);

}  // namespace caps
