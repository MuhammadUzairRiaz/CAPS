#include <cctype>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <map>

#include "caps/analysis.hpp"
#include "caps/crystal.hpp"
#include "caps/elements.hpp"
#include "caps/io.hpp"
#include "io_util.hpp"

namespace caps {

std::string detect_format(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw ReadError("cannot open " + path);
  {   // binary trajectories and names that say what they are
    const std::string e = lower(std::filesystem::path(path).extension().string());
    const std::string stem = lower(std::filesystem::path(path).filename().string());
    if (e == ".xtc") return "xtc";
    if (e == ".trr") return "trr";
    if (e == ".dcd") return "dcd";
    if (e == ".sdf" || e == ".sd" || e == ".mol" || e == ".mdl") return "sdf";
    if (e == ".poscar" || e == ".vasp" || stem.rfind("poscar", 0) == 0 || stem.rfind("contcar", 0) == 0) return "poscar";
  }
  std::string l1, l2, l3;
  std::getline(in, l1);
  std::getline(in, l2);
  std::getline(in, l3);
  if (l1.rfind("ITEM:", 0) == 0) return "lammps-dump";
  {
    const std::string e = lower(std::filesystem::path(path).extension().string());
    if (e == ".cif") return "cif";
    if (e == ".car") return "car";
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
  // an MDL molfile: the counts line (fourth) ends with V2000 or V3000
  if (l1.size() || l2.size()) {
    std::string l4;
    std::getline(in, l4);
    if (l4.find("V2000") != std::string::npos || l4.find("V3000") != std::string::npos) return "sdf";
  }
  throw ReadError(path + ": format not recognised (supported: LAMMPS data, dump and DCD, GROMACS .gro/.top, .xtc and .trr, PDB, mol2, SDF/MOL, XYZ, CIF, VASP POSCAR, Materials Studio .car/.mdf)");
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
  if (f == "car") return "Materials Studio .car";
  if (f == "xtc") return "GROMACS .xtc";
  if (f == "trr") return "GROMACS .trr";
  if (f == "dcd") return "DCD";
  if (f == "sdf") return "MDL molfile / SD";
  if (f == "poscar") return "VASP POSCAR";
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
  } else if (fmt == "car") {
    tr.topology = read_car(path);
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
  } else if (fmt == "xtc" || fmt == "trr" || fmt == "dcd") {
    // coordinates only: the atoms from the structure or topology given with them
    if (topology_path.empty())
      throw ReadError(path + ": a " + format_name(fmt) + " trajectory holds coordinates only — open it with its structure (.gro, .pdb, .data, .top …)");
    const std::string tl0 = lower(std::filesystem::path(topology_path).extension().string());
    System top;
    if (tl0 == ".top" || tl0 == ".itp") top = read_gromacs_topology(topology_path);
    else {
      const Trajectory st = open_file(topology_path);
      top = st.frame(0);
      top.bonds = st.topology.bonds;
    }
    const auto prog = [&](double f, const Trajectory& t) {
      if (!told) tell(t, true);
      return report(3, f, std::to_string(t.frames()) + " frames");
    };
    tr = fmt == "xtc" ? read_xtc(path, top, progress.max_frames, prog) : fmt == "trr" ? read_trr(path, top, progress.max_frames, prog)
                      : read_dcd(path, top, progress.max_frames, prog);
    for (size_t i = 0; i < tr.topology.atoms.size(); ++i) tr.topology.atoms[i].pos = tr.positions.front()[i];
    if (!top.bonds.empty()) tr.topology.bonds_from_file = true;
  } else if (fmt == "sdf" || fmt == "poscar") {
    tr.topology = fmt == "sdf" ? read_sdf(path) : read_poscar(path);
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
  if (tr.topology.atoms.empty()) throw ReadError(path + ": no atoms found (read as " + format_name(fmt) + ")");
  // a GROMACS topology: its atoms (types, charges, residues), bonds and explicit terms, the coordinates from the file
  const std::string tl = topology_path.size() > 4 ? topology_path.substr(topology_path.size() - 4) : std::string();
  const bool gmx_top = tl == ".top" || tl == ".itp";
  if (gmx_top) {
    System T = read_gromacs_topology(topology_path);
    if (T.atoms.size() != tr.topology.atoms.size())
      throw ReadError(topology_path + " has " + std::to_string(T.atoms.size()) + " atoms, " + path + " " + std::to_string(tr.topology.atoms.size()));
    for (size_t i = 0; i < T.atoms.size(); ++i) {
      T.atoms[i].pos = tr.topology.atoms[i].pos;
      T.atoms[i].image = tr.topology.atoms[i].image;
    }
    T.cell = tr.topology.cell;
    T.title = tr.topology.title.empty() ? T.title : tr.topology.title;
    T.velocities = tr.topology.velocities;
    T.unwrapped = tr.topology.unwrapped;
    T.notes.insert(T.notes.begin(), tr.topology.notes.begin(), tr.topology.notes.end());
    T.notes.push_back("topology from " + std::filesystem::path(topology_path).filename().string() + ": " + std::to_string(T.bonds.size()) + " bonds, " +
                      std::to_string(T.topology->angles.size()) + " angles, " + std::to_string(T.topology->dihedrals.size()) + " dihedrals, " +
                      std::to_string(T.topology->vsites.size()) + " virtual sites");
    tr.topology = std::move(T);
  }
  if (!told) tell(tr, (fmt == "lammps-dump" || gmx_top) && !topology_path.empty());
  report(3, 1, std::to_string(tr.frames()) + " frames");
  if (tr.topology.bonds.empty() && !tr.topology.bonds_from_file) {   // a file that declares its bonds (0 too) keeps them
    System f0 = tr.frame(0);
    tr.topology.bonds = perceive_bonds(f0);
    tr.topology.notes.push_back(std::to_string(tr.topology.bonds.size()) + " bonds perceived from distances (none in file)");
  }
  return tr;
}

FileInspection inspect_file(const std::string& path, const std::string& topology_path, int head_lines) {
  FileInspection r;
  r.format = detect_format(path);
  r.format_name = format_name(r.format);
  r.bytes = size_t(std::filesystem::file_size(path));
  std::ifstream in(path, std::ios::binary);
  std::string line;
  const bool binary = r.format == "xtc" || r.format == "trr" || r.format == "dcd";
  if (binary) {
    r.notes.push_back(std::string(format_name(r.format)) + " holds coordinates only: the atoms, bonds and types come from the topology or structure given with it" +
                      (topology_path.empty() ? " (none chosen yet)" : ""));
    head_lines = 0;
  }
  // first lines, and the dump's header for its columns
  for (int k = 0; k < head_lines && std::getline(in, line); ++k) {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    r.head.push_back(line.size() > 200 ? line.substr(0, 200) + " …" : line);
  }
  if (r.format == "lammps-dump") {
    in.clear();
    in.seekg(0);
    size_t natoms = 0;
    bool got_cols = false;
    while (std::getline(in, line)) {
      if (line.rfind("ITEM: NUMBER OF ATOMS", 0) == 0 && std::getline(in, line)) natoms = std::stoull(trim(line));
      else if (line.rfind("ITEM: ATOMS", 0) == 0 && !got_cols) {
        got_cols = true;
        for (const auto& c : split(line.substr(11))) {
          FileColumn col{c, "", "float", true};
          if (c == "id") col = {c, "Particle Identifier", "int", true};
          else if (c == "mol") col = {c, "Molecule Identifier", "int", true};
          else if (c == "type") col = {c, "Particle Type", "int", true};
          else if (c == "element") col = {c, "Element", "text", true};
          else if (c == "q") col = {c, "Charge", "float", true};
          else if (c == "x" || c == "y" || c == "z") col.maps_to = "Position." + std::string(1, char(std::toupper(c[0]))) + " (wrapped)";
          else if (c == "xu" || c == "yu" || c == "zu") col.maps_to = "Position." + std::string(1, char(std::toupper(c[0]))) + " (unwrapped)";
          else if (c == "xs" || c == "ys" || c == "zs") col.maps_to = "Position." + std::string(1, char(std::toupper(c[0]))) + " (scaled)";
          else if (c == "ix" || c == "iy" || c == "iz") col = {c, "Periodic image", "int", true};
          else { col.maps_to = "not read"; col.used = false; }
          r.columns.push_back(col);
        }
        // skip this frame's atoms quickly
        for (size_t i = 0; i < natoms && std::getline(in, line); ++i) {}
        r.frames = 1;
      } else if (line.rfind("ITEM: TIMESTEP", 0) == 0 && got_cols) {
        ++r.frames;
      }
    }
    r.atoms = natoms;
    size_t unused = 0;
    for (const auto& c : r.columns) unused += !c.used;
    if (unused) r.notes.push_back(std::to_string(unused) + " column(s) are not read (CAPS keeps positions, ids, molecules, types and charges)");
  }
  // types: from the data file (the file itself or the topology)
  const std::string data = r.format == "lammps-data" ? path : topology_path;
  if (!data.empty()) {
    try {
      const System t = read_lammps_data(data);
      for (const auto& ti : t.types) {
        FileType ft{ti.type, ti.label, ti.mass, ""};
        for (const auto& a : t.atoms) if (a.type == ti.type && a.element) { ft.element = element(a.element).symbol; break; }
        r.types.push_back(ft);
      }
      if (r.atoms == 0) r.atoms = t.atoms.size();
      if (r.frames == 0) r.frames = 1;
      r.bonds_from = t.bonds.empty() ? "perceived from distances (the data file has none)" :
                     std::to_string(t.bonds.size()) + " bonds from " + std::filesystem::path(data).filename().string();
      r.units = "LAMMPS real (Å, fs, kcal/mol, e)";
    } catch (const std::exception& e) {
      r.notes.push_back(std::string("topology: ") + e.what());
    }
  }
  if (r.types.empty() || r.format != "lammps-dump") {
    // other formats: open to count atoms and types (small files); frames as read
    try {
      const Trajectory t = open_file(path, topology_path);
      r.atoms = t.topology.atoms.size();
      r.frames = t.frames();
      if (r.types.empty()) {
        std::map<std::string, FileType> by;
        for (const auto& a : t.topology.atoms) {
          const std::string key = a.name.empty() ? element(a.element).symbol : a.name;
          auto& ft = by[key];
          ft.label = key;
          ft.type = a.type;
          ft.mass = element(a.element).mass;
          ft.element = element(a.element).symbol;
        }
        for (auto& [k, v] : by) r.types.push_back(v);
      }
      if (r.bonds_from.empty()) r.bonds_from = t.topology.bonds_from_file ? std::to_string(t.topology.bonds.size()) + " bonds from the file" : "perceived from distances";
      if (r.units.empty()) r.units = r.format == "gro" ? "nm in the file, converted to Å" : "Å";
    } catch (const std::exception& e) {
      r.notes.push_back(e.what());
    }
  }
  if (r.bonds_from.empty()) r.bonds_from = "perceived from distances (no topology file)";
  return r;
}

}  // namespace caps
