// CAPS import (design/boards/ImportDialog): opening a structure file with the choices the file leaves open — where bonds
// come from (perceived from distances with a tolerance, read from the file, or none), bond orders and aromaticity
// written into the bonds, molecules by connectivity, molecules made whole across periodic boundaries, the cell kept or
// dropped — and a preview of those choices on a small fragment before the whole file is read.
#pragma once
#include <string>
#include <vector>

#include "caps/io.hpp"
#include "caps/system.hpp"

namespace caps {

struct ImportOptions {
  enum Bonds { Perceive = 0, FromFile = 1, None = 2 };
  int bonds = Perceive;
  double tolerance = 0.45;   // Å added to the sum of covalent radii (Cordero 2008)
  bool bond_orders = true;   // orders (1 2 3, 4 aromatic) from valences written into the bonds
  bool split = true;         // molecule ids from bond-connected components
  bool unwrap = true;        // molecules made whole across the cell's faces
  bool use_cell = true;      // false: the file's cell is dropped (a non-periodic structure)
};

// Applies the options to a trajectory read by open_file (bonds_in_file: how many bonds the file itself gave).
void apply_import(Trajectory& t, const ImportOptions& o, size_t bonds_in_file);

Trajectory import_file(const std::string& path, const std::string& topology_path, const ImportOptions& o);

struct ImportPreview {
  FileInspection file;
  size_t atoms = 0, bonds_in_file = 0, bonds = 0, molecules = 0;
  int single = 0, dbl = 0, triple = 0, aromatic = 0;
  std::string cell;                 // "36.2 × 36.2 × 36.2 Å, orthogonal", or "" without one
  System fragment;                  // the first `heavy` connected heavy atoms (carbons first) with their hydrogens
  int fragment_heavy = 0;
  std::vector<std::string> notes;
};
// Frame 0 only, so the preview is quick for long trajectories.
ImportPreview import_preview(const std::string& path, const ImportOptions& o, int heavy = 10);

// Bond orders of a structure without hydrogens (heavy atoms only) from its geometry: hybridisation from bond angles
// (terminal atoms: the bond length), flat sp² six-rings aromatic (4), flat five-rings conjugated, then double and triple
// bonds paired shortest first; a planar N with three bonds keeps its lone pair.
void orders_from_geometry(System& s);

}  // namespace caps
