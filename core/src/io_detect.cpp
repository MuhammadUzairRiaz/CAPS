#include <cstdio>
#include <filesystem>
#include <fstream>

#include "caps/analysis.hpp"
#include "caps/crystal.hpp"
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
  {
    const std::string e = lower(std::filesystem::path(path).extension().string());
    if (e == ".cif") return "cif";
  }
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
  throw ReadError(path + ": format not recognised (supported: LAMMPS data and dump, GROMACS .gro, PDB, mol2, XYZ, CIF)");
}

Trajectory open_file(const std::string& path, const std::string& topology_path) { return open_file(path, topology_path, OpenProgress{}); }

namespace {
const char* format_name(const std::string& f) {
  if (f == "lammps-dump") return "LAMMPS dump";
  if (f == "lammps-data") return "LAMMPS data";
  if (f == "gro") return "GROMACS .gro";
  if (f == "pdb") return "PDB";
  if (f == "mol2") return "Tripos mol2";
  if (f == "cif") return "CIF";
  return "XYZ";
}
}  // namespace

Trajectory open_file(const std::string& path, const std::string& topology_path, const OpenProgress& progress) {
  const auto report = [&](int stage, double f, const std::string& detail) { return !progress.report || progress.report(stage, f, detail); };
  const std::string fmt = detect_format(path);
  report(0, 1, format_name(fmt));
  bool told = false;
  const auto tell = [&](const Trajectory& t, bool joined) {
    told = true;
    const System& s = t.topology;
    char b[160];
    if (s.cell.valid())
      std::snprintf(b, sizeof b, "%zu atoms · box %.1f × %.1f × %.1f Å", s.atoms.size(), norm(s.cell.a), norm(s.cell.b), norm(s.cell.c));
    else
      std::snprintf(b, sizeof b, "%zu atoms · no cell", s.atoms.size());
    report(1, 1, b);
    if (joined) std::snprintf(b, sizeof b, "%s · %zu bonds · %zu atom types", std::filesystem::path(topology_path).filename().string().c_str(), s.bonds.size(), s.types.size());
    else if (!s.bonds.empty()) std::snprintf(b, sizeof b, "%zu bonds in the file", s.bonds.size());
    else std::snprintf(b, sizeof b, "no topology file · bonds from distances");
    report(2, 1, b);
  };
  Trajectory tr;
  if (fmt == "lammps-data") {
    tr.topology = read_lammps_data(path);
    std::vector<Vec3> p;
    for (const auto& a : tr.topology.atoms) p.push_back(a.pos);
    tr.positions.push_back(std::move(p));
    tr.cells.push_back(tr.topology.cell);
    tr.timesteps.push_back(0);
  } else if (fmt == "lammps-dump") {
    System top;
    if (!topology_path.empty()) top = read_lammps_data(topology_path);
    tr = read_lammps_dump(path, topology_path.empty() ? nullptr : &top, progress.max_frames, [&](double f, const Trajectory& t) {
      if (!told) tell(t, !topology_path.empty());
      return report(3, f, std::to_string(t.frames()) + " frames");
    });
  } else if (fmt == "mol2") {
    tr.topology = read_mol2(path);
    std::vector<Vec3> p;
    for (const auto& a : tr.topology.atoms) p.push_back(a.pos);
    tr.positions.push_back(std::move(p));
    tr.cells.push_back(tr.topology.cell);
    tr.timesteps.push_back(0);
  } else if (fmt == "cif") {
    tr.topology = read_cif(path);
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
  if (!told) tell(tr, fmt == "lammps-dump" && !topology_path.empty());
  report(3, 1, std::to_string(tr.frames()) + " frames");
  if (tr.topology.bonds.empty()) {
    System f0 = tr.frame(0);
    tr.topology.bonds = perceive_bonds(f0);
    tr.topology.notes.push_back(std::to_string(tr.topology.bonds.size()) + " bonds perceived from distances (none in file)");
  }
  return tr;
}

}  // namespace caps
