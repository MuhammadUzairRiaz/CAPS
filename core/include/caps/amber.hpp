// AMBER topologies (prmtop / parm7) and coordinates (inpcrd / rst7 / crd): the structure with the force field the
// topology file carries, term by term — no typing, no library: the charges, masses, Lennard-Jones A/B tables, bonds,
// angles, torsions (1-4 scaling from SCEE / SCNB) and exclusions exactly as the file gives them (AMBER file formats,
// ambermd.org/FileFormats.php). A force field CAPS cannot hold whole — chamber (CHARMM) topologies, CMAP, 12-6-4 ions,
// extra points, polarisable terms, mixed 1-4 scaling — is not taken at all (ff null, a note says why): the structure
// still opens with its atoms, residues, bonds, charges and masses. Perturbation and cap topologies are refused.
#pragma once
#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "caps/field.hpp"
#include "caps/system.hpp"

namespace caps {

struct AmberTopology {
  System system;                             // atoms (names = AMBER atom types), residues, molecules, bonds, box
  std::shared_ptr<const ForceField> ff;      // the file's force field for these atoms, or null (see notes)
  std::vector<std::string> notes;
};

AmberTopology read_amber_prmtop(const std::string& path);

// Restart / coordinate file, ASCII (inpcrd / rst7) or NetCDF (ncrst, the AMBERRESTART convention): positions (Å),
// velocities when present (converted to Å/fs), the box when present (lengths and angles).
struct AmberCoordinates {
  std::vector<Vec3> positions, velocities;
  Cell cell;
  bool has_box = false;
  std::string title;
};
AmberCoordinates read_amber_coordinates(const std::string& path);

// AMBER trajectories on a topology with the same atoms (a prmtop, or any structure): NetCDF (the AMBER convention,
// ambermd.org/netcdf/nctraj.xhtml; NetCDF-3 classic or 64-bit offset — a NetCDF-4 / HDF5 file is refused with how to
// convert it) with coordinates, velocities and the box per frame; ASCII mdcrd (10F8.3, a box line after each frame
// when the topology has a box). Frames are numbered 0, 1, 2 …; the times, where the file has them, go in the notes.
Trajectory read_amber_netcdf(const std::string& path, const System& topology, size_t max_frames = 0,
                             const std::function<bool(double, const Trajectory&)>& progress = {}, const FrameSelection& frames = {});
Trajectory read_amber_mdcrd(const std::string& path, const System& topology, size_t max_frames = 0,
                            const std::function<bool(double, const Trajectory&)>& progress = {}, const FrameSelection& frames = {});

// The structure with its force field as an AMBER topology and restart (STEM.prmtop, STEM.inpcrd) for AMBER, OpenMM and
// ParmEd: charges, masses, a Lennard-Jones A/B table for every type pair (whatever the mixing rule: each pair as the force
// field mixes it), harmonic bonds and angles, periodic torsions and impropers, the 1-4 pairs carried by the torsions
// (SCEE = 1/coul14, SCNB = 1/lj14; a 1-4 pair no torsion term reaches gets a zero-barrier term, as LEaP writes them),
// the exclusions, residues, molecules and the box. A force field AMBER's functional forms cannot hold (class II, Urey–
// Bradley, harmonic impropers, separate 1-4 Lennard-Jones, non-LJ pairs, virtual sites, many-body terms, a relative
// permittivity) is refused with the reason. No generalised-Born radii are written. Returns notes.
std::vector<std::string> write_amber(const System& s, const ForceField& ff, const std::string& stem);

bool is_amber_topology_path(const std::string& path);      // .prmtop, .parm7, .prmtop.gz …
bool is_amber_coordinates_path(const std::string& path);   // .inpcrd, .rst7, .restrt, .crd, .rst

}  // namespace caps
