// CAPS peptide builder (design/boards/BioBuilder): a peptide from its one-letter sequence and a secondary structure per
// residue, all-atom with hydrogens.
//
//  Backbone      N, CA, C, O placed by NeRF from φ, ψ, ω with Engh & Huber geometry (N–CA 1.458, CA–C 1.525, C–N
//                1.329, C=O 1.231 Å; N–CA–C 111.2°, CA–C–N 116.2°, C–N–CA 121.7°). α-helix φ/ψ −57/−47 (Pauling,
//                Corey & Branson 1951), β-strand −139/135, PPII −75/145; coil drawn from the αR, β and PPII basins.
//                Proline keeps φ −65°.
//  Side chains   Each residue's free amino acid (L, written N[C@@H](R)C(=O)O with its charge at the pH) is embedded
//                from SMILES and superposed on the residue's N, CA, C; its side chain and hydrogens are taken over.
//  Termini       NH3+, NH2 or acetyl (ACE); COO−, COOH or N-methyl amide (NME). Other hydrogens placed by geometry.
//  Clean-up      UFF minimisation of the whole peptide (optional).
#pragma once
#include <array>
#include <cstdint>
#include <string>
#include <vector>

#include "caps/system.hpp"

namespace caps {

struct PeptideOptions {
  std::string sequence;            // one-letter codes; spaces, digits and line breaks are ignored
  std::string structure;           // per residue H (α-helix), E (β-strand), P (PPII), C (coil); missing residues are C
  std::array<double, 3> helix{-57, -47, 180}, strand{-139, 135, 180}, ppii{-75, 145, 180};   // φ, ψ, ω (degrees)
  std::string n_term = "NH3+";     // NH3+ | NH2 | ACE
  std::string c_term = "COO-";     // COO- | COOH | NME
  double ph = 7.0;                 // side-chain charges: D E below pKa 3.9/4.3 neutral, H+ below 6.0, K+ below 10.5, R+ below
                                   // 12.5, C− above 8.3, Y− above 10.1
  bool neutral = false;            // every side chain uncharged (overrides pH)
  bool cleanup = true;             // UFF minimisation after placement
  double ftol = 1.0;               // kcal/mol/Å for the clean-up
  uint64_t seed = 1;               // coil torsions
};

struct PeptideReport {
  std::string smiles, formula, structure;   // structure: the per-residue letters used
  int residues = 0, charge = 0;
  size_t atoms = 0;
  double mass = 0;                           // g/mol
  std::vector<std::string> notes;
};

// The peptide (one molecule, no cell), atoms named as in PDB files (N CA C O CB … H HA …) with three-letter residue names.
System build_peptide(const PeptideOptions& o, PeptideReport* report = nullptr);

// The sequence of the first record of a FASTA text (or of plain one-letter text).
std::string parse_fasta(const std::string& text);
// One-letter code to its three-letter name ("ALA"); empty when unknown.
std::string residue_name(char code);
// φ, ψ, ω of a structure letter (H E P; C and others give the coil's centre −70, 140, 180).
std::array<double, 3> structure_torsions(const PeptideOptions& o, char ss);

}  // namespace caps
