// CAPS semicrystalline cells: lamellar stacks of a polymer crystal with amorphous layers between the lamellae, the
// cells the papers on semicrystalline polyethylene simulate (Balijepalli & Rutledge, J. Chem. Phys. 109, 6523 (1998);
// in 't Veld & Rutledge, Macromolecules 36, 7358 (2003); Lu et al., J. Appl. Phys. 123, 015107 (2018); Kumar et al.,
// Macromolecules 50, 1206 (2017)).
//
// The stack normal is z. A lamella of nc unit cells along the chain axis c (= z) holds straight crystalline stems; each
// stem runs through the whole lamella and leaves it at both faces. The amorphous layer between the lamella's top face
// and the bottom face of its periodic image above holds the chain segments that join the stems' ends:
//   loop:   leaves a face and comes back into the same face at another stem;
//   bridge: leaves one face and enters the other (it crosses the amorphous layer);
//   tail:   leaves a face and ends free in the amorphous layer.
// Every chain is stem, segment, stem, …, with a tail at both ends; no chain closes on itself. Counts and fractions are
// per amorphous segment (a loop and a bridge use two stem ends, a tail one).
#pragma once
#include <array>
#include <cstdint>
#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "caps/field.hpp"
#include "caps/system.hpp"

namespace caps {

// A polymer crystal's unit cell for the lamellar builder: the atoms of one cell (c is the chain axis, stems along it),
// with the source the numbers come from.
struct PolymerCrystalCell {
  std::string id, name, source;
  System cell;   // atoms of one unit cell (fractional positions in [0, 1)), orthorhombic or monoclinic about c
};

// A unit cell from CAPS's polymer crystal library (data_dir/crystals/polymers/polymers.json, each entry a CIF from its
// paper: "pe-orthorhombic" (alias "pe"), Bunn, Trans. Faraday Soc. 35, 482 (1939), Pnam, a 7.40, b 4.93, c 2.534 Å;
// "pe-monoclinic", "ipp-alpha", "pvdf-beta" …), or a CIF / structure file with a cell when id_or_path is a file.
PolymerCrystalCell polymer_crystal_cell(const std::string& id_or_path, const std::string& data_dir);

// Stems of a crystal: each chain whose bonds (perceived by distance, through the cell's faces) repeat with a translation
// of one cell edge (c for the lamellar builder; a or b reported too). For each stem, the atoms of one cell's repeat and
// its backbone atoms (the shortest path of atoms with two heavy neighbours or more from one cell's chain to the next) in
// order along the axis; helical when the backbone spreads across the axis in two directions (not a planar zig-zag).
struct CrystalStem {
  std::vector<uint32_t> atoms;      // atoms of the stem in one unit cell
  std::vector<uint32_t> backbone;   // backbone atoms in one cell, ordered along the axis
  Vec3 axis{0, 0, 1};               // unit vector of the stem's translation
  int lattice_axis = -1;            // the cell edge the chain repeats along: 0 a, 1 b, 2 c (−1: not a chain)
  std::string formula;              // of one cell's repeat ("C2H4")
  bool helical = false;             // the backbone turns about the axis (more than a planar zig-zag)
};
std::vector<CrystalStem> find_stems(const System& cell, double bond_tolerance = 0.45);

struct SemicrystallineOptions {
  std::string polymer = "pe";      // a library id (pe = pe-orthorhombic) or a structure file of the unit cell
  std::string data_dir;            // CAPS's data folder (holds crystals/polymers)
  int nx = 4, ny = 6;              // lateral repeats of the unit cell (a, b)
  int nc = 0;                      // cells along c in a lamella; 0: from crystal_nm
  double crystal_nm = 5.0;         // lamella thickness (used when nc == 0)
  double amorphous_nm = 4.0;       // amorphous layer thickness
  // target fractions of the amorphous segments (normalised); a loop or bridge takes two stem ends, a tail one
  double loops = 0.5, bridges = 0.2, tails = 0.3;
  // where a loop re-enters its face: "adjacent" (mostly a neighbouring stem: 70 % of loops at the nearest free stems,
  // the rest anywhere on the face) or "random" (anywhere on the face: the switchboard model)
  std::string reentry = "adjacent";
  double density = 0.855;          // g/cm³ of the amorphous layer (melt-like PE at 300 K, as Lu et al. 2018)
  uint64_t seed = 1;
  bool united_atom = false;        // carbons only (a united-atom force field folds the H later)
  double contact = 3.0;            // Å, closest approach of non-bonded backbone carbons while growing
  int trials = 24;                 // torsion trials per growth step (dense near a closure)
  std::function<bool(const std::string& stage, double fraction)> progress;   // return false to cancel
};

// Topology and structure of a lamellar cell (from the builder, or read from any frame by analyse_semicrystalline).
struct SemicrystallineReport {
  int nx = 0, ny = 0, nc = 0;
  int stems = 0, chains = 0;
  int loops = 0, bridges = 0, tails = 0;
  int loops_top = 0, loops_bottom = 0, tails_top = 0, tails_bottom = 0;
  std::vector<int> bridge_units, loop_units, tail_units;   // backbone atoms (CH2 units) of each segment
  std::vector<int> chain_units;                              // backbone atoms per chain
  double lx = 0, ly = 0, lz = 0;                             // cell edges, Å
  double crystal_thickness = 0, amorphous_thickness = 0;     // Å
  double amorphous_density = 0, density = 0;                 // g/cm³ (amorphous: the layer's atoms over its volume)
  double crystallinity_p2 = 0;     // fraction of backbone atoms whose local P2 (chord over ±3 neighbours) exceeds p2_cut
  double crystallinity_mass = 0;   // mass fraction of the atoms in the crystal group
  double p2_cut = 0.75;
  double tilt_deg = 0, tilt_sd = 0;   // angle of the stems (end to end) to the stack normal
  std::vector<double> z, rho, p2;     // profiles along the stack normal: bin centres (Å), density (g/cm³), mean P2
  size_t atoms = 0;
  double worst_contact = 0;        // closest non-bonded backbone pair as grown, Å
  int eased_steps = 0, strained_closures = 0;
  std::vector<std::string> notes;
  // analyse_semicrystalline only: each atom's class, the backbone carbons by their segment and every other atom as the
  // carbon it is bonded to (hydrogens): kNone (a chain without a stem, or not a chain), kStem, kLoop, kBridge, kTail
  enum AtomClass : int8_t { kNone = 0, kStem = 1, kLoop = 2, kBridge = 3, kTail = 4 };
  std::vector<int8_t> atom_class;
  std::vector<int32_t> atom_segment;   // the segment (stem, loop, bridge or tail) each atom belongs to, 0-based; −1 none
};

// Builds the lamellar cell: crystal from the unit cell (stems along z), amorphous segments grown from the stems' ends
// at the amorphous density, avoiding the crystal and each other (carbons; hydrogens added after), each closing as a
// loop, a bridge or a tail in the requested fractions. Atoms carry the tags "crystal" and "amorphous", molecule ids per
// chain; positions wrapped in the cell. Throws std::invalid_argument for options it cannot meet (atom budget, a
// density too low for the topology, a cell whose stems are not polymethylene) and std::runtime_error when cancelled.
System build_semicrystalline(const SemicrystallineOptions& o, SemicrystallineReport* report = nullptr);

// Push-off and minimisation of a built cell: the crystal held while the amorphous layer is pushed apart (LJ force cap
// raised stage by stage), then everything minimised. field null: the built-in GAFF of C and H.
void relax_semicrystalline(System& s, std::shared_ptr<const ForceField> field, const EnergyOptions& e, double fmax = 1.0,
                           std::function<bool(const std::string&, double)> progress = {});

// The topology and structure of a lamellar cell from one frame (stack normal z, one lamella per cell): backbone
// carbons with at most two carbon neighbours, a backbone atom crystalline when the P2 of its chord (r_{i+1} − r_{i−1})
// against the stems' mean direction, averaged over ±3 neighbours, exceeds p2_cut; stems are runs of at least
// min_stem crystalline atoms that cross the lamella's middle plane; the segments between stems are loops (both ends at
// one face), bridges (one end at each) or tails. bins: density / P2 profile bins along z.
SemicrystallineReport analyse_semicrystalline(const System& s, double p2_cut = 0.75, int min_stem = 8, int bins = 60);

// The report as JSON (counts, fractions, distributions, profiles).
std::string semicrystalline_json(const SemicrystallineReport& r);
// A short human summary (several lines).
std::string semicrystalline_summary(const SemicrystallineReport& r);

// Cluster inputs for the equilibration of a built cell (LAMMPS, data file already written with the force field):
// minimisation, NVT at 300 K with the crystal tethered (fix spring/self), then NPT (x, y, z coupled separately) with the
// tether released, density and box written. Returns the input text; data_name is the data file as the input reads it.
struct SemicrystallineEquilibration {
  double temperature = 300, pressure = 1;   // K, atm
  double restrained_ps = 200, free_ps = 2000;
  double spring = 10;                        // kcal/mol/Å², the crystal's tether in the first stage
  double dt = 1;                             // fs
};
std::string semicrystalline_lammps_run(const SemicrystallineEquilibration& q);

}  // namespace caps
