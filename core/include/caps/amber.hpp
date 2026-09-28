// AMBER topologies (prmtop / parm7) and coordinates (inpcrd / rst7 / crd): the structure with the force field the
// topology file carries, term by term — no typing, no library: the charges, masses, Lennard-Jones A/B tables, bonds,
// angles, torsions (1-4 scaling from SCEE / SCNB) and exclusions exactly as the file gives them (AMBER file formats,
// ambermd.org/FileFormats.php). CHARMM topologies converted by chamber, CMAP, extra points and polarisable
// topologies are refused rather than read in part.
#pragma once
#include <memory>
#include <string>
#include <vector>

#include "caps/field.hpp"
#include "caps/system.hpp"

namespace caps {

struct AmberTopology {
  System system;                             // atoms (names = AMBER atom types), residues, molecules, bonds, box
  std::shared_ptr<const ForceField> ff;      // the file's force field for these atoms
  std::vector<std::string> notes;
};

AmberTopology read_amber_prmtop(const std::string& path);

// ASCII restart / coordinate file: positions (Å), velocities when present (converted to Å/fs), the box when present
// (lengths and angles). Throws for a NetCDF restart.
struct AmberCoordinates {
  std::vector<Vec3> positions, velocities;
  Cell cell;
  bool has_box = false;
  std::string title;
};
AmberCoordinates read_amber_coordinates(const std::string& path);

bool is_amber_topology_path(const std::string& path);      // .prmtop, .parm7, .prmtop.gz …
bool is_amber_coordinates_path(const std::string& path);   // .inpcrd, .rst7, .restrt, .crd, .rst

}  // namespace caps
