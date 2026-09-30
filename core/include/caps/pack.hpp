// CAPS Pack: packing rigid molecules into regions without overlaps (a native replacement for packmol).
//
// Objective, after Martínez, Andrade, Birgin and Martínez, J. Comput. Chem. 30, 2157 (2009):
//   f = Σ_{i<j, different molecules} max(0, d² − |x_i − x_j|²)²  +  w Σ_i g_region(x_i)²
// over rigid bodies (centre + rotation), minimised with L-BFGS, with poorly placed molecules moved between rounds.
#pragma once
#include <cstdint>
#include <functional>
#include <stdexcept>
#include <algorithm>
#include <string>
#include <vector>

#include "caps/system.hpp"

namespace caps {

struct Region {
  enum Kind { InsideBox, OutsideBox, InsideSphere, OutsideSphere, InsideCylinder, OutsideCylinder, OverPlane, BelowPlane,
              InsideEllipsoid, OutsideEllipsoid };   // ellipsoid: a centre, b semi-axes, r = d (packmol: Σ((x−a)/b)² ≤ d)
  Kind kind = InsideBox;
  Vec3 a{0, 0, 0};   // box: low corner · sphere: centre · cylinder: base point · plane: normal (a·x = r)
  Vec3 b{0, 0, 0};   // box: high corner · cylinder: axis direction (unit)
  double r = 0;      // sphere / cylinder radius · plane offset
  double length = 0; // cylinder length along the axis
  std::vector<int> atoms;   // the molecule's atoms it holds for (0-based, sorted; packmol "atoms … end atoms"); empty: all
  bool applies(int atom) const { return atoms.empty() || std::binary_search(atoms.begin(), atoms.end(), atom); }
  // Violation distance (Å) of a point; 0 when it satisfies the region.
  double violation(const Vec3& x) const;
};

struct PackItem {
  std::string name;
  System molecule;              // one molecule (its atoms, bonds, types, charges)
  int count = 1;
  std::vector<Region> regions;  // all must hold for every atom
  bool fixed = false;           // placed as given (optionally moved by `center` / `angles`), never optimised
  bool center = false;          // fixed: put the molecule's centre at `position`
  Vec3 position{0, 0, 0};
  Vec3 angles{0, 0, 0};         // fixed: rotations about x, y, z (radians, packmol convention)
};

struct PackProgress {
  std::string stage;            // "insertion", "optimisation", "verification", "compression"
  int loop = 0, loops = 0;
  double penalty = 0;
  double dmin = 0;              // smallest intermolecular distance so far, Å
  int bad = 0;                  // molecules still in violation
};

struct PackOptions {
  double tolerance = 2.0;       // Å, minimum distance between atoms of different molecules
  Cell cell;                    // the output cell; with periodic = true pairs use the minimum image
  bool periodic = true;
  uint64_t seed = 1;
  int max_loops = 60;           // optimisation rounds (each followed by moving the worst molecules)
  int iterations = 400;         // L-BFGS iterations per round
  int trials = 16;              // placements tried per molecule during insertion
  double move_fraction = 0.05;  // share of molecules in violation moved between rounds
  int threads = 0;              // 0 = automatic
  // Stage 3 (a dense cell): pack loosely in the given cell, then compress it to this density (g/cm³) by affine steps with
  // capped-force push-off and minimisation (Relax's compression, the built-in force field); 0 = off. Periodic, no fixed
  // structures. packmol text: "compress 0.9".
  double compress_to = 0;
  std::function<bool(const PackProgress&)> progress;   // return false to cancel
};

struct PackReport {
  bool success = false;
  double dmin = 0;              // smallest distance between atoms of different molecules, Å
  int close_pairs = 0;          // pairs closer than the tolerance
  double region_violation = 0;  // largest distance of an atom outside its regions, Å
  int molecules = 0, atoms = 0;
  int loops = 0, iterations = 0, evaluations = 0, moved = 0;
  double seconds = 0;
  std::vector<std::string> notes;
};

struct PackError : std::runtime_error {
  using std::runtime_error::runtime_error;
};

// Packs the items. Throws PackError (with the report filled when given) when the tolerance or a region cannot be met:
// a cell below tolerance is never returned.
System pack(const std::vector<PackItem>& items, const PackOptions& o, PackReport* report = nullptr);

// packmol input files (the common subset): tolerance, seed, output, filetype, pbc, structure … end structure with
// number, inside / outside box, cube, sphere, cylinder, over / below plane, fixed, center. Relative structure paths
// are resolved against `base_dir`. Returns the items and fills options (cell, tolerance, seed) and the output path.
// CAPS Pack input, the readable form (the Studio writes it; packmol files are read too):
//   cell      40 40 40                     Å, periodic (or: cell 0 0 0 to 40 40 40)
//   distance  2.0                          Å: no two atoms of different molecules closer
//   seed      1                            (or: seed new)
//   compress  0.95                         g/cm³ after packing (optional) · save FILE · loops N · iterations N
//   molecule  water.pdb                    a structure file (spaces allowed), then its lines, then `end`
//     count   100
//     in      box from 0 0 0 to 40 40 40   · cube from x y z size d · sphere at x y z radius r
//             ellipsoid at x y z axes a b c scale d · cylinder from x y z along u v w radius r length l
//     not in  (the same regions)           outside
//     above   plane normal a b c at d      a·x ≥ d (below: a·x ≤ d)
//     fixed   at x y z rotated a b c degrees     placed as given, never moved
//     centred                              the fixed position is its centre
//     atoms   1 2 3 in sphere at …         a region for some of its atoms only
//   end
bool is_caps_pack_input(const std::string& text);
// The CAPS form as packmol's (throws PackError naming the line), and packmol's as the CAPS form.
std::string caps_pack_to_packmol(const std::string& text);
std::string packmol_to_caps_pack(const std::string& text);
std::vector<PackItem> read_packmol_input(const std::string& path, PackOptions& o, std::string* output);
// The same from text; `name` labels error messages.
std::vector<PackItem> parse_packmol_input(const std::string& text, const std::string& base_dir, PackOptions& o, std::string* output,
                                          const std::string& name = "input");

// Inserts `count` copies of `guest` (one molecule) into the free space of a periodic structure, which stays where it is
// (curatives, sulfur donors, solvent into a polymer cell). The host keeps its molecules, bonds and bond orders; the
// guests are numbered after them. Uses o.tolerance, o.seed and the optimisation settings; the cell is the host's.
System insert_molecules(const System& host, const System& guest, int count, const PackOptions& o, PackReport* report = nullptr,
                        const std::vector<Region>& regions = {});   // every guest atom must satisfy these (a pore)

// Smallest distance between atoms of different molecules (minimum image when the cell is periodic) and the number of
// such pairs closer than `tolerance`. Uses a cell list.
std::pair<double, int> intermolecular_contacts(const System& s, double tolerance, bool periodic);
// The same between groups (one int per atom; pairs within a group are not counted), e.g. the structures Pack placed.
std::pair<double, int> intermolecular_contacts(const System& s, double tolerance, bool periodic, const std::vector<int>& group);

}  // namespace caps
