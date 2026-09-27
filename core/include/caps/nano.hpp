// CAPS nanostructures: graphene sheets, carbon nanotubes and nanoparticles cut from crystals — the fillers of rubber
// composites — and a filler embedded in a grown polymer matrix.
//
//  graphene_sheet  Rectangular graphene (zigzag along x), one or more AB-stacked layers 3.35 Å apart; periodic in the
//                  plane, or a flake with its edge carbons capped by hydrogen.
//  nanotube        (n, m) tube rolled from graphene (or several concentric armchair / zigzag walls) (Saito, Dresselhaus & Dresselhaus, 1998): diameter
//                  a √(n² + nm + m²) / π, a = √3 × C–C; periodic along z, or a finite tube with hydrogen-capped ends.
//  nanoparticle    A sphere, cube, octahedron, cuboctahedron, truncated octahedron or icosahedron (the size is the
//                  circumscribed radius; an icosahedron cut from the bulk keeps its lattice: it is not multiply twinned), or a fibre
//                  (a cylinder periodic along the crystal's c axis) cut from a bulk
//                  crystal, centred on an atom or on the cell centre; isolated atoms dropped; optionally passivated
//                  (silanols on silica) as the surface builder does, the new groups pointing away from the centre.
//  embed_filler    The filler held fixed at the centre of a periodic cell (its periodic axes kept) while polymer chains
//                  grow around it to the matrix density: filled rubber, CNT– or graphene–rubber composites.
#pragma once
#include <array>
#include <string>
#include <vector>

#include "caps/polymer.hpp"
#include "caps/system.hpp"

namespace caps {

struct NanoReport {
  double diameter = 0;        // tube diameter, or particle diameter (Å)
  double chiral_angle = 0;    // degrees
  double translation = 0;     // tube period |T| (Å)
  int atoms_per_period = 0;
  int capped = 0;             // hydrogens on edges / ends
  int added_h = 0, added_oh = 0;
  std::vector<std::string> notes;
};

struct SheetOptions {
  double lx = 20.0, ly = 20.0;   // Å, rounded to whole rectangular cells (2.46 × 4.26 Å at C–C 1.42)
  double cc = 1.42;               // C–C, Å
  int layers = 1;                 // AB stacked, 3.35 Å apart
  bool periodic = true;           // periodic in the plane; false: a flake with hydrogen on the edge carbons
  double vacuum = 15.0;           // Å above and below (and around a flake)
};
System graphene_sheet(const SheetOptions& o, NanoReport* rep = nullptr);

struct NanotubeOptions {
  int n = 10, m = 10;
  double length = 25.0;           // Å, rounded to whole periods
  double cc = 1.42;
  bool periodic = true;           // periodic along z; false: finite, ends capped with hydrogen
  double vacuum = 10.0;           // Å around the tube
  // multi-walled: concentric armchair or zigzag walls about wall_spacing apart ((5,5)@(10,10)@(15,15) at 3.4 Å)
  int walls = 1;
  double wall_spacing = 3.4;      // Å, the target; whole (n, m) steps give 3.39 Å (armchair) or 3.52 Å (zigzag)
};
System nanotube(const NanotubeOptions& o, NanoReport* rep = nullptr);
// Diameter (Å), chiral angle (degrees) and period |T| (Å) of an (n, m) tube.
std::array<double, 3> nanotube_geometry(int n, int m, double cc = 1.42);

// Fibre: a cylinder along the crystal's c axis (which must be normal to a and b), periodic along it — a glass, carbon or
// steel fibre to embed in a rubber matrix.
enum class ParticleShape { Sphere, Cube, Octahedron, Cuboctahedron, Fibre, TruncatedOctahedron, Icosahedron };
const char* to_string(ParticleShape s);
ParticleShape particle_shape_from_string(const std::string& s);

struct ParticleOptions {
  ParticleShape shape = ParticleShape::Sphere;
  double radius = 12.0;           // circumscribed radius (fibre: its radius), Å
  double length = 20.0;           // fibre: length along the axis, Å (whole cells; the fibre is periodic along it)
  bool on_atom = true;            // centre on the atom nearest the cell centre; false: on the cell centre
  bool passivate = false;
  double vacuum = 10.0;
};
System nanoparticle(const System& bulk, const ParticleOptions& o, NanoReport* rep = nullptr);

struct FillerMatrixOptions {
  double density = 0.9;           // matrix density, g/cm³ (the chains' mass over the cell volume left by the filler)
  int chains = 10;
  std::array<bool, 3> keep_axis{false, false, false};   // cell edges taken from the filler (a tube's z, a sheet's x and y)
  GrowOptions grow;               // seed, contact scale, trials, progress …; auto_scale is on by default
};
struct FillerReport {
  double filler_mass_fraction = 0, filler_volume_fraction = 0, density = 0;
  std::array<double, 3> cell{0, 0, 0};
  std::vector<std::string> notes;
};
// The filler (molecule 1, held while growing) centred in a new periodic cell with `chains` chains of spec around it.
System embed_filler(const System& filler, const ChainSpec& spec, const FillerMatrixOptions& o, FillerReport* rep = nullptr);

// Volume (Å³) inside the van der Waals spheres of the atoms (a 0.4 Å grid over their bounding box).
double occupied_volume(const System& s);

// Pores (design/boards/SlitPore): a fluid confined in
//   Slit       two graphene walls (1–3 AB-stacked sheets each) a width H apart (carbon centre to carbon centre),
//              periodic in the plane; periodic along z (the walls meet through the boundary) or with vacuum above
//   Cylinder   a cylindrical channel of diameter D along c, carved from a crystal block periodic in all directions,
//              the carved surface optionally passivated (silanols pointing into the pore)
//   Framework  a crystal's supercell as it is (zeolite, MOF): the fluid goes into its free space
// The walls are molecule 1 (hold them in Relax and Dynamics); the fluid is packed by CAPS Pack inside the pore only.
enum class PoreKind { Slit, Cylinder, Framework };
struct PoreOptions {
  PoreKind kind = PoreKind::Slit;
  double width = 10.0;            // slit H (C–C centres) or cylinder diameter, Å
  int layers = 1;                 // slit: sheets per wall
  double lx = 26.0, ly = 22.0;    // slit: in-plane size (whole graphene cells)
  bool vacuum = false;            // slit: vacuum above the upper wall (non-periodic in z)
  double vacuum_gap = 20.0;
  const System* crystal = nullptr;   // cylinder / framework: the bulk crystal (a CIF)
  double wall = 6.0;              // cylinder: material around the channel, Å
  double length = 20.0;           // cylinder: along c (whole cells)
  int repeat[3] = {2, 2, 2};      // framework: supercell
  bool passivate = false;         // cylinder: silanols / H on the carved surface
  const System* fluid = nullptr;  // one molecule; null: an empty pore
  int count = 0;
  double tolerance = 2.0;         // Å between fluid and wall atoms (and between fluid molecules)
  uint64_t seed = 1;
};
struct PoreReport {
  int wall_atoms = 0, fluid_molecules = 0;
  double width = 0;               // as built (rounded to the lattice)
  double pore_volume = 0;         // Å³ of the region the fluid may occupy
  double fluid_density = 0;       // g/cm³ over that region
  double dmin = 0;                // closest fluid–wall or fluid–fluid contact, Å
  std::vector<std::string> notes;
};
System build_pore(const PoreOptions& o, PoreReport* rep = nullptr);

}  // namespace caps
