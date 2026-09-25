// One molecule described (design/boards/MoleculeInspector) and solvent-accessible surface areas
// (design/boards/SurfaceArea).
#pragma once
#include <string>
#include <vector>

#include "caps/system.hpp"

namespace caps {

struct MoleculeInfo {
  int molecule = 0;                    // 0-based, in the order of System::molecules
  int atoms = 0, bonds = 0;
  std::string formula;                 // Hill order
  std::string smiles;                  // written from the perceived bonds (not canonical)
  double mass = 0;                     // g/mol, standard atomic weights
  double monoisotopic = 0;             // u, the most abundant isotopes; NaN when an element is not tabulated
  double dbe = 0;                      // double-bond equivalents (rings + π bonds)
  double inertia[3] = {0, 0, 0};       // principal moments, amu·Å², ascending
  double inertia_defect = 0;           // I_C − I_A − I_B (0 for a planar molecule)
  double rg = 0;                       // mass-weighted radius of gyration, Å
  double net_charge = 0;               // e
  bool has_charges = false;
  double dipole = 0;                   // debye, from the charges about the centre of mass (NaN without charges)
  int rotatable = 0;                   // acyclic single bonds between two non-terminal heavy atoms
  int rings = 0;
};
// The molecule containing atom `atom` (positions made whole across the cell).
MoleculeInfo molecule_info(const System& s, uint32_t atom);

struct SasaResult {
  std::vector<double> area;            // Å² per atom
  double total = 0;
  double probe = 1.4;
  int points = 200;
};
// Shrake & Rupley (1973): points on a golden spiral at r_vdW (Bondi) + probe, each counted when no other atom's sphere
// covers it; periodic cells use the minimum image.
SasaResult sasa(const System& s, double probe = 1.4, int points = 200);

}  // namespace caps
