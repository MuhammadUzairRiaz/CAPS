#pragma once
#include <functional>
#include <stdexcept>
#include <string>

#include "caps/system.hpp"

namespace caps {

struct ReadError : std::runtime_error {
  using std::runtime_error::runtime_error;
};

// Formats: "lammps-data", "lammps-dump", "gro", "xyz", "pdb", "mol2". Detected from content, then extension.
std::string detect_format(const std::string& path);

System read_lammps_data(const std::string& path);
Trajectory read_lammps_dump(const std::string& path, const System* topology = nullptr);
// The same, stopping after max_frames frames (0: all); progress gets the share of the file read and the frames so far
// after each frame, and stops reading when it returns false (the frames read so far are kept, with a note).
Trajectory read_lammps_dump(const std::string& path, const System* topology, size_t max_frames,
                            const std::function<bool(double, const Trajectory&)>& progress);
Trajectory read_gro(const std::string& path);
Trajectory read_xyz(const std::string& path);
Trajectory read_pdb(const std::string& path);
// Tripos mol2: atoms with their type (Atom::name: SYBYL "C.ar" or a force-field type such as GAFF "ca"), atom
// name kept in Atom::resname's companion field label, charges, substructures as molecules, bonds with orders,
// and the CRYSIN cell. Several MOLECULE records are read as one system.
System read_mol2(const std::string& path);

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

void write_lammps_data(const System& s, const std::string& path);
// All frames as a LAMMPS text dump with unwrapped coordinates (id mol type xu yu zu).
void write_lammps_dump(const Trajectory& t, const std::string& path);
void write_xyz(const System& s, const std::string& path);
// GROMACS .gro: one residue per molecule (its residue name, else MOL), positions and box in nm.
void write_gro(const System& s, const std::string& path);
void write_pdb(const System& s, const std::string& path);
// Tripos mol2 with Atom::name as the atom type (force-field types or SYBYL), charges and bond orders.
void write_mol2(const System& s, const std::string& path);

}  // namespace caps
