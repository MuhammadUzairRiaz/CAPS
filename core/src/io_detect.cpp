#include <filesystem>
#include <fstream>

#include "caps/analysis.hpp"
#include "caps/io.hpp"
#include "io_util.hpp"

namespace caps {

std::string detect_format(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw ReadError("cannot open " + path);
  std::string l1, l2, l3;
  std::getline(in, l1);
  std::getline(in, l2);
  std::getline(in, l3);
  if (l1.rfind("ITEM:", 0) == 0) return "lammps-dump";
  for (const auto* l : {&l1, &l2, &l3})
    if (l->find("@<TRIPOS>") != std::string::npos) return "mol2";
  for (const auto* l : {&l1, &l2, &l3})
    for (const char* rec : {"CRYST1", "ATOM  ", "HETATM", "HEADER", "MODEL ", "REMARK", "COMPND", "TITLE "})
      if (l->rfind(rec, 0) == 0) return "pdb";
  const std::string ext = lower(std::filesystem::path(path).extension().string());
  // XYZ: first line is only an integer, and the third line starts with a symbol.
  auto t1 = split(l1);
  if (t1.size() == 1 && t1[0].find_first_not_of("0123456789") == std::string::npos) {
    auto t3 = split(l3);
    if (!t3.empty() && std::isalpha(static_cast<unsigned char>(t3[0][0])) && t3.size() >= 4) return "xyz";
  }
  // GRO: second line is only an integer, third line has fixed columns.
  auto t2 = split(l2);
  if (t2.size() == 1 && t2[0].find_first_not_of("0123456789") == std::string::npos && l3.size() >= 44) return "gro";
  if (ext == ".mol2") return "mol2";
  if (ext == ".gro") return "gro";
  if (ext == ".pdb" || ext == ".ent") return "pdb";
  if (ext == ".xyz" || ext == ".extxyz") return "xyz";
  if (ext == ".lammpstrj" || ext == ".dump") return "lammps-dump";
  if (ext == ".data" || ext == ".lmp" || ext == ".lammps") return "lammps-data";
  // A data file has a header with "atoms" on an early line.
  std::string line;
  for (int k = 0; k < 20 && std::getline(in, line); ++k) {
    auto t = split(line);
    if (t.size() >= 2 && t[1] == "atoms") return "lammps-data";
  }
  throw ReadError(path + ": format not recognised (supported: LAMMPS data and dump, GROMACS .gro, PDB, mol2, XYZ)");
}

Trajectory open_file(const std::string& path, const std::string& topology_path) {
  const std::string fmt = detect_format(path);
  Trajectory tr;
  if (fmt == "lammps-data") {
    tr.topology = read_lammps_data(path);
    std::vector<Vec3> p;
    for (const auto& a : tr.topology.atoms) p.push_back(a.pos);
    tr.positions.push_back(std::move(p));
    tr.cells.push_back(tr.topology.cell);
    tr.timesteps.push_back(0);
  } else if (fmt == "lammps-dump") {
    if (!topology_path.empty()) {
      System top = read_lammps_data(topology_path);
      tr = read_lammps_dump(path, &top);
    } else {
      tr = read_lammps_dump(path);
    }
  } else if (fmt == "mol2") {
    tr.topology = read_mol2(path);
    std::vector<Vec3> p;
    for (const auto& a : tr.topology.atoms) p.push_back(a.pos);
    tr.positions.push_back(std::move(p));
    tr.cells.push_back(tr.topology.cell);
    tr.timesteps.push_back(0);
  } else if (fmt == "gro") {
    tr = read_gro(path);
  } else if (fmt == "pdb") {
    tr = read_pdb(path);
  } else {
    tr = read_xyz(path);
  }
  if (tr.topology.bonds.empty()) {
    System f0 = tr.frame(0);
    tr.topology.bonds = perceive_bonds(f0);
    tr.topology.notes.push_back(std::to_string(tr.topology.bonds.size()) + " bonds perceived from distances (none in file)");
  }
  return tr;
}

}  // namespace caps
