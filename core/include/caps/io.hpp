#pragma once
#include <functional>
#include <stdexcept>
#include <string>

#include "caps/system.hpp"

namespace caps {

struct ReadError : std::runtime_error {
  using std::runtime_error::runtime_error;
};

// Formats: "lammps-data", "lammps-dump", "gro", "xyz", "pdb", "mol2", "cif", "car", "sdf", "poscar", and the binary
// trajectories "xtc", "trr", "dcd" (by extension). Detected from content, then extension.
std::string detect_format(const std::string& path);

System read_lammps_data(const std::string& path);
Trajectory read_lammps_dump(const std::string& path, const System* topology = nullptr);
// The same, stopping after max_frames frames (0: all); progress gets the share of the file read and the frames so far
// after each frame, and stops reading when it returns false (the frames read so far are kept, with a note).
Trajectory read_lammps_dump(const std::string& path, const System* topology, size_t max_frames,
                            const std::function<bool(double, const Trajectory&)>& progress);
Trajectory read_gro(const std::string& path);
// A GROMACS topology (.top or .itp; #include and #ifdef resolved, force-field includes that are not found skipped):
// the atoms of every molecule [ molecules ] lists, named by their types, with charges, residues, masses, and the
// explicit topology (bonds, constraints as stiff bonds, angles, dihedrals, virtual sites, exclusions) in
// System::topology. No coordinates: open_file(COORDS, TOP) joins them.
System read_gromacs_topology(const std::string& path, std::vector<std::string>* notes = nullptr);
Trajectory read_xyz(const std::string& path);
// Binary trajectories (coordinates only; the atoms come from topology, which must have as many): GROMACS .xtc
// (compressed; nm → Å) and .trr (single or double precision; velocities of the first frame kept), CHARMM / NAMD / LAMMPS
// .dcd (either byte order; the cell of each frame when the file has it). max_frames and progress as read_lammps_dump.
Trajectory read_xtc(const std::string& path, const System& topology, size_t max_frames = 0,
                    const std::function<bool(double, const Trajectory&)>& progress = {});
Trajectory read_trr(const std::string& path, const System& topology, size_t max_frames = 0,
                    const std::function<bool(double, const Trajectory&)>& progress = {});
Trajectory read_dcd(const std::string& path, const System& topology, size_t max_frames = 0,
                    const std::function<bool(double, const Trajectory&)>& progress = {});
// MDL molfile / SD file (V2000 and V3000): atoms, bonds with orders, formal charges (M  CHG); each record a molecule.
System read_sdf(const std::string& path);
// VASP POSCAR / CONTCAR (VASP 5 species line, or VASP 4 with the species in the title): cell (scale or volume),
// direct or Cartesian positions, selective dynamics skipped.
System read_poscar(const std::string& path);
Trajectory read_pdb(const std::string& path);
// Tripos mol2: atoms with their type (Atom::name: SYBYL "C.ar" or a force-field type such as GAFF "ca"), atom
// name kept in Atom::resname's companion field label, charges, substructures as molecules, bonds with orders,
// and the CRYSIN cell. Several MOLECULE records are read as one system.
System read_mol2(const std::string& path);
// Materials Studio .car (coordinates, cell, force-field types, charges; molecules) with the bonds of the .mdf beside it
// (bonds across the cell included). Each atom's name is its force-field type.
System read_car(const std::string& path);
// Writes STEM.car and STEM.mdf (path's extension replaced): the cell as Materials Studio orients it, atom names as the
// potential types, charges, one block per molecule, the bonds (with their cell images) in the .mdf.
void write_car(const System& s, const std::string& path);

// Open any supported file. A topology (LAMMPS data) may be given for dumps.
Trajectory open_file(const std::string& path, const std::string& topology_path = "");

// Staged opening (the Studio's progressive open, design/boards/VisLoading): report is called as each stage finishes —
// 0 format detected, 1 frame 0 read, 2 topology joined — and during stage 3 with the share of the file read. Returning
// false stops reading; the frames read so far are kept. max_frames > 0 stops after that many frames, so frame 0 can be
// shown at once and the rest read in the background.
struct OpenProgress {
  std::function<bool(int stage, double fraction, const std::string& detail)> report;
  size_t max_frames = 0;
};
Trajectory open_file(const std::string& path, const std::string& topology_path, const OpenProgress& progress);

// What a file holds before it is opened (design/boards/OpenLammps, OpenGromacs): format from content, first lines,
// dump columns and what they map to, types with masses and elements (from the file or the topology), frames (counted
// by scanning, without reading coordinates) and where bonds will come from.
struct FileColumn { std::string name, maps_to, kind; bool used = true; };
struct FileType { int type = 0; std::string label; double mass = 0; std::string element; };
struct FileInspection {
  std::string format, format_name;
  std::vector<std::string> head;       // first lines
  std::vector<FileColumn> columns;     // dump columns
  std::vector<FileType> types;
  size_t atoms = 0, frames = 0, bytes = 0;
  std::string bonds_from, units;
  std::vector<std::string> notes;
};
FileInspection inspect_file(const std::string& path, const std::string& topology_path = "", int head_lines = 40);

void write_lammps_data(const System& s, const std::string& path);
// All frames as a LAMMPS text dump with unwrapped coordinates (id mol type xu yu zu).
void write_lammps_dump(const Trajectory& t, const std::string& path);
void write_xyz(const System& s, const std::string& path);
// GROMACS .gro: one residue per molecule (its residue name, else MOL), positions and box in nm.
void write_gro(const System& s, const std::string& path);
void write_pdb(const System& s, const std::string& path);
// Tripos mol2 with Atom::name as the atom type (force-field types or SYBYL), charges and bond orders.
void write_mol2(const System& s, const std::string& path);
// MDL SD file (V2000; V3000 above 999 atoms or bonds), the formal charges perceived from the structure
void write_sdf(const System& s, const std::string& path);
// CIF in P 1: the cell and every atom's fractional coordinates (throws without a periodic cell)
void write_cif(const System& s, const std::string& path);
// DCD as LAMMPS writes it (every frame with its cell); dt_fs is the time step recorded in the header
void write_dcd(const Trajectory& t, const std::string& path, double dt_fs = 1.0);
// A moltemplate system (.lt text) from the LAMMPS data and input files CAPS writes for a structure with its force field:
// In Init (units, styles, special bonds), Data Masses, In Settings (pair and every coefficient, class II cross terms as
// bb / ba / mbt / ebt / at / aat / bb13 / aa), Data Boundary, and the atoms (unwrapped), bonds, angles, dihedrals and
// impropers, the types named after the force field's. moltemplate.sh gives back the same energies.
std::string lammps_to_moltemplate(const std::string& data_path, const std::string& input_path, const std::string& title = "");

}  // namespace caps
