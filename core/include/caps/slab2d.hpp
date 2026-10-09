#pragma once
// 2D sheets and terminated surfaces for DFT (the DFT surface & adsorption workbench), general for any layered material:
// a sheet built from a list of layers (element, stacking position or fractional x y, height) — MXenes, other carbides
// and nitrides, hexagonal 2D materials — or ONE layer cut out of a bulk, stacked or MAX-type structure; sites on the
// outer layer of any slab (hollows classified fcc / hcp by the atom beneath, top, bridge); terminations or adatoms on
// those sites from a library (bond length → height, mixed and Janus faces, a seed that gives the same pattern on every
// platform); the POSCAR conventions VASP needs; and a validator that knows what a terminated 2D slab must look like.
// Presets (data/sheets/sheets.json) and the termination library (data/sheets/terminations.json) are data, not code.
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "caps/json.hpp"
#include "caps/system.hpp"

namespace caps {

// ---------------------------------------------------------------- sheets from layers
struct SheetLayer {
  std::string element;
  std::string site;               // A (0,0), B (2/3,1/3), C (1/3,2/3) of a hexagonal cell; "" = use fx, fy
  double fx = 0, fy = 0;          // fractional in-plane position (when site is "")
  double z = 0;                   // Å above the first layer
};
struct SheetSpec {
  std::string name;
  std::string lattice = "hexagonal";   // hexagonal (a1 = a x, a2 = a(−½, √3/2)) | rectangular (a x, b y)
  double a = 0, b = 0;            // Å (b: rectangular only)
  std::vector<SheetLayer> layers; // one atom per layer and cell
  double c = 30.0;                // cell height (Å)
  std::string source;             // where the geometry comes from
};
System sheet_from_layers(const SheetSpec& s);
// A preset by name from data/sheets/sheets.json (case-insensitive), with `a` (> 0) overriding its lattice constant
// (heights scaled by the same ratio, said in notes). Throws for an unknown name; preset_names lists them.
SheetSpec sheet_preset(const std::string& data_dir, const std::string& name, double a = 0, std::vector<std::string>* notes = nullptr);
std::vector<std::string> sheet_preset_names(const std::string& data_dir);

struct IsolateOptions {
  int index = 0;                  // which complete layer (from the bottom)
  double gap = 2.0;               // Å: atoms farther apart than this along z belong to different layers
  std::string formula;            // the layer's composition up to a multiple, e.g. Ti3C2 ("" = any)
  std::vector<std::string> remove;   // elements taken out first (a MAX phase's A element: Al, Si …)
  double c = 30.0;                // cell height of the result
};
struct IsolateReport { int layers = 0, complete = 0; std::string formula; std::vector<std::string> notes; };
// ONE complete layer out of a bulk, stacked or MAX-type structure: atoms clustered along z (periodic: a layer crossing
// the cell boundary is joined), the `remove` elements taken out first, a layer complete when its composition is a
// multiple of `formula`. Throws when there is none.
System isolate_layer(const System& src, const IsolateOptions& o, IsolateReport* rep = nullptr);

// ---------------------------------------------------------------- sites and terminations
struct SurfaceSite {
  Vec3 pos;                       // in the outer layer's plane (Å)
  std::string kind;               // fcc (above the 3rd layer) | hcp (above the 2nd layer) | top | bridge | hollow
  std::string face;               // top | bottom
  std::string beneath;            // element directly beneath (hollows, top)
  double r = 0;                   // in-plane distance to the outer-layer neighbours (Å)
};
// The sites of both faces: three-fold hollows (centroids of nearest-neighbour triangles of the outer layer), top
// (the outer atoms) and bridge (midpoints of nearest-neighbour pairs). The outer layer is the heavy atoms within 0.2 Å
// of the extreme z (of `element` only when given). Sorted by rounded fractional coordinates (−0 made 0): the same
// order on every platform.
std::vector<SurfaceSite> surface_sites(const System& slab, const std::string& element = "");

struct TerminationKind {
  std::string name;               // O, OH, F, Cl, NH …
  std::string atom;               // bonded to the surface
  std::string tail;               // H of OH / NH ("" when none), along the outward normal
  double bond = 0;                // surface atom – X (Å)
  double tail_bond = 0;           // X – tail (Å)
  double lo = 0, hi = 0;          // the validator's range for the surface–X distance (Å)
  std::string source;
};
// The library for a surface element (data/sheets/terminations.json: values given per surface element where they are
// known, else the sum of covalent radii, Cordero 2008, as a starting guess — said in `source`).
std::vector<TerminationKind> termination_library(const std::string& data_dir, const std::string& surface_element);

struct FaceSpec {
  std::vector<std::pair<std::string, double>> fractions;   // {O 0.5, OH 0.25, F 0.25}; empty: a bare face
  std::string site = "fcc";
};
struct TerminateOptions {
  FaceSpec top, bottom;
  bool janus = false;             // false: the bottom takes the top's species, site and counts (no net dipole)
  int na = 1, nb = 1;             // in-plane supercell
  double vacuum = 20.0;           // Å, centred along z
  uint64_t seed = 1;              // arrangement of a mixed face: numpy default_rng(seed).shuffle, top face first
  std::string surface_element;    // the outer layer's element ("" = the most common heavy element of the outer layer)
  std::map<std::string, double> bond;       // overrides: "O" → 2.05, "OH" → 2.2, "OH.tail" → 0.97
  std::vector<std::string> order = {"Ti", "C", "N", "O", "F", "H"};   // species order of the output (POTCAR order)
  std::string data_dir;
};
struct TerminateReport {
  std::map<std::string, int> top, bottom;
  double r = 0;                   // in-plane distance of the sites used (Å)
  std::vector<std::string> notes;
};
// The slab repeated na × nb, each face's sites filled (largest-remainder counts, shuffled with the seed), heights
// h = √(d² − r²) from the bond (error if d ≤ r; top sites: h = d), tails along the outward normal, vacuum centred,
// positions wrapped in-plane, atoms grouped by element in `order`.
System terminate_slab(const System& slab, const TerminateOptions& o, TerminateReport* rep = nullptr);
std::vector<std::pair<std::string, double>> parse_fractions(const std::string& text);   // "O:0.5,OH:0.25,F:0.25" or "O"

// ---------------------------------------------------------------- species order and POSCAR (VASP 5)
System group_by_species(const System& s, const std::vector<std::string>& order = {"Ti", "C", "N", "O", "F", "H"});
std::string formula_ordered(const System& s, const std::vector<std::string>& order = {"Ti", "C", "N", "O", "F", "H"});
// Cartesian, %.10f, species contiguous (throws if not), the comment as given; `fixed` true = F F F (Selective dynamics).
void write_vasp_poscar(const System& s, const std::string& path, const std::string& comment, const std::vector<bool>& fixed = {});
// Tolerant reader: a label that is not an element (the digit 0 for O) is mapped and reported in `bad`.
System read_vasp_poscar(const std::string& path, std::map<std::string, std::string>* bad = nullptr, std::string* comment = nullptr);
// Periodic pairs within `cutoff` (images included; i == j only for a non-zero image): j, the vector i → j (Å), |d|.
struct Neighbour { uint32_t i, j; Vec3 d; double r; };
std::vector<Neighbour> neighbours_pbc(const System& s, double cutoff, std::array<bool, 3> pbc = {true, true, true});
// Positions wrapped into the cell along the periodic axes the way ASE does it (fractions in [−1e-7, 1 − 1e-7)).
void wrap_ase(System& s, std::array<bool, 3> pbc = {true, true, false});

// ---------------------------------------------------------------- the 2D validator
struct Finding {
  std::string level;              // ERROR | WARN
  std::string text;
  std::vector<int> atoms;         // 0-based, for the 3D view
};
struct Validation2D {
  std::string status = "PASS";
  std::vector<Finding> findings;
  Json info = Json::object();
  std::string text(const std::string& name = "") const;
  Json json() const;
};
struct ValidateOptions {
  std::string expect;             // expected formula per cell, e.g. Ti3C2O2 (× N)
  std::string core;               // the bare layer's formula, e.g. Ti3C2 ("" = from the expected formula's
                                  // first two elements, else the largest block's composition)
  std::string comment;            // the POSCAR comment, compared with the atoms
  std::map<std::string, std::string> bad_labels;
  std::string data_dir;
};
// PASS/FAIL with findings: non-element labels, composition × N, complete blocks (one expected) and their distinct
// layers, incomplete fragments, atoms floating > 2 Å from the block, vacuum (< 10 Å error, < 15 warn), close contacts
// (1.4 Å; 0.9 O–H; 1.0 with H), isolated atoms (none within 3.2 Å), surface–termination distances in range and
// coordination (hollow 3, top 1), each H on exactly one O within 0.93–1.03 Å pointing outward, site class per face,
// top/bottom asymmetry (→ dipole correction), a comment-line formula that disagrees with the atoms.
Validation2D validate_2d(const System& s, const ValidateOptions& o = {});

}  // namespace caps
