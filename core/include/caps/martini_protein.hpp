// Martini 2.2 proteins (data/martini/martini22-protein.json, converted from vermouth-martinize by
// bench/ff/convert_vermouth_martini22.py): an all-atom protein mapped onto Martini beads with the model's explicit
// topology, as martinize writes it; and DSSP secondary structure (Kabsch & Sander, Biopolymers 22, 2577 (1983)).
#pragma once
#include <string>
#include <vector>

#include "caps/system.hpp"

namespace caps {

// A residue of an all-atom protein: its atoms, name, chain position.
struct ProteinResidue {
  std::string name;              // as in the file (HIS, ASP ...)
  int64_t resid = 0;
  int chain = 0;                 // index of the chain (bonded run of residues)
  std::vector<uint32_t> atoms;
  int n = -1, ca = -1, c = -1, o = -1;   // backbone atoms, −1 when missing
};

// Residues in file order, chains split where the peptide bond C(i)–N(i+1) is missing.
std::vector<ProteinResidue> protein_residues(const System& s);

// DSSP secondary structure, one letter per residue of protein_residues: H G I E B T S, 'C' for none (DSSP's blank).
// Amide hydrogens are placed as DSSP places them (N + the unit vector from the previous O to C), the given ones ignored.
std::string dssp(const System& s);

// vermouth's conversion of DSSP letters to Martini's (1 2 3 at helix ends and short helices, H, E, T, S, C).
std::string dssp_to_martini(const std::string& dssp_letters, const std::string& data_path);

struct MartiniProteinReport {
  std::string dssp, cg_ss;       // one letter per residue
  int residues = 0, beads = 0, chains = 0, disulfides = 0;
  std::vector<std::string> residue_names;   // the Martini residue each became (HSE, ASP0 ...)
  std::vector<std::string> notes;
};

// The protein's beads (named by their Martini bead types, element 0, at the mass-weighted centre of the atoms the
// mapping lists, residue names and ids kept) with the explicit Martini 2.2 topology in System::topology. `ss`: one letter
// per residue in DSSP's alphabet, or empty for DSSP on the structure. Throws std::invalid_argument for residues the
// model has no mapping for or missing backbone atoms. Other molecules (water, ions ...) are not converted: they must be
// absent (the caller splits them off).
System martini22_protein(const System& aa, const std::string& ss, const std::string& data_path, MartiniProteinReport* rep = nullptr);

// The topology as a GROMACS .itp (nm, kJ/mol, the Martini functions; stiff bonds written back as constraints when their
// force constant is the constraint one), for comparison with martinize2.
std::string martini_itp(const System& beads, const std::string& data_path);

}  // namespace caps
