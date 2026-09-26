// Martini 2.2 proteins (data/martini/martini22-protein.json, converted from vermouth-martinize by
// bench/ff/convert_vermouth_martini22.py): an all-atom protein mapped onto Martini beads with the model's explicit
// topology, as martinize writes it; and DSSP secondary structure (Kabsch & Sander, Biopolymers 22, 2577 (1983)).
#pragma once
#include <functional>
#include <string>
#include <vector>

#include "caps/json.hpp"
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

// DSSP secondary structure, one letter per residue of protein_residues that has backbone atoms (waters, ions and
// ligands are skipped): H G I E B T S, 'C' for none (DSSP's blank).
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

// Martini 3 proteins as martinize2 builds them (data/martini/martini3-protein.json: vermouth-martinize's martini3001 force
// field and mappings, converted by bench/ff/convert_vermouth_martini3.py). The options are martinize2's.
struct Martini3Options {
  std::string ss;                  // DSSP letters per residue, one letter for all (-ss C), "" for DSSP on the structure, "-" for none
  bool scfix = true;               // side-chain fix dihedrals and angles (martinize2's default; -noscfix turns it off)
  bool neutral_termini = false;    // -nt
  bool extdih = false;             // dihedrals rather than elastic bonds for extended regions (-extdih)
  bool disulfides = true;          // bonded SG pairs become disulfide bridges (-cys auto)
  std::vector<std::pair<int64_t, int64_t>> idr;   // disordered regions by residue number (-id-regions)
  bool elastic = false;            // the rubber-band elastic network (-elastic) between backbone beads
  double ef = 700, el = 0, eu = 0.9, ea = 0, ep = 1, es = 0, em = 0;   // kJ/mol/nm², nm; decay exp(−ea (d − es)^ep)
  int ermd = -1;                   // minimum residue separation (−1: vermouth's default, 2)
  std::string eunit = "molecule";  // molecule, chain or all
  double constraint_kj = 1e6;      // constraints as stiff bonds (kJ/mol/nm², Martini 3's own stiff_fc)
  // small molecules' beads at the geometric centre of their mapped atoms, hydrogens included: the mapping Martini 3's
  // small-molecule parameters reproduce (bench/ff/check_martini3_small.py); false: mass-weighted, as vermouth centres
  bool small_geometric = true;
  // molecules that are neither amino acids, Martini 3 small molecules, water nor ions: refused (typing must not drop
  // part of a structure); false: left out with a note (caps martini)
  bool refuse_unmatched = true;
};
System martini3_protein(const System& aa, const Martini3Options& o, const std::string& data_path, MartiniProteinReport* rep = nullptr);
// The topology as a GROMACS .itp (stiff bonds from constraints written back as constraints), for comparison with martinize2
std::string martini3_itp(const System& beads, double constraint_kj = 1e6);
// Martini 3 small molecules (data/martini/martini3-small-molecules.json: vermouth's martini3001 blocks and mappings,
// the CHARMM residues they are written for; bench/ff/convert_vermouth_martini3_small.py): every molecule of an
// all-atom structure whose heavy atoms and bonds (and hydrogen counts, when it has hydrogens) are one of those residues,
// whatever its atom names, becomes that molecule's beads (mass-weighted centres of the mapped atoms, virtual sites
// placed) with its explicit topology. Molecules that match none are listed in `unmatched` (formulas) and left out.
System martini3_small_molecules(const System& aa, const std::string& data_path, double constraint_kj = 1e6,
                                std::vector<std::string>* unmatched = nullptr, std::vector<std::string>* notes = nullptr, bool geometric = false);
// An all-atom structure in Martini 3: its amino acids as martini3_protein, its other molecules as small molecules;
// water and ions are left out (noted).
System martini3_all_atom(const System& aa, const Martini3Options& o, const std::string& protein_path, const std::string& small_path,
                         MartiniProteinReport* rep = nullptr);

// Whether a model file is a vermouth model (Martini 3) rather than the Martini 2.2 tables
bool is_martini3_model(const std::string& path);

// A molecule given as GROMACS terms (data/martini/martini3-molecules.json: atoms with type, name, charge and optional mass;
// bonds, constraints, angles, dihedrals, exclusions and virtual sites by 0-based index, nm, kJ/mol, degrees, the .itp
// functions): beads named by their types, placed by the bead builder at the bond lengths, with the explicit topology.
// Constraints become stiff bonds of constraint_kj (kJ/mol/nm²); type_mass gives each type's mass for mass-weighted sites.
System gromacs_molecule(const Json& mol, const std::function<double(const std::string&)>& type_mass, double constraint_kj, uint64_t seed = 1);

// The topology as a GROMACS .itp (nm, kJ/mol, the Martini functions; stiff bonds written back as constraints when their
// force constant is the constraint one), for comparison with martinize2.
std::string martini_itp(const System& beads, const std::string& data_path);

}  // namespace caps
