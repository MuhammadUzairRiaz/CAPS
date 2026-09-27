// CAPS appearance (design/boards/Appearance): what the Studio draws besides balls and sticks.
//
//  stereo_labels     R/S of tetrahedral centres by the Cahn–Ingold–Prelog rules (rule 1a: atomic number, explored
//                    sphere by sphere through a hierarchical digraph with duplicate atoms for double and triple bonds),
//                    read off the 3D geometry.
//  surface_mesh      A molecular surface by marching tetrahedra on a grid: solvent-accessible (van der Waals radii plus
//                    the probe; Lee & Richards 1971), van der Waals, or solvent-excluded (the accessible surface pulled
//                    back by the probe radius; Connolly 1983, on the grid's distance field).
//  surface_potential The Coulomb potential of the atoms' partial charges at the mesh vertices (kcal/mol/e, ε = 1),
//                    for colouring the surface.
//  polyhedra         Coordination polyhedra: the convex hull of the bonded neighbours of 4- to 8-coordinated centres
//                    (SiO4 tetrahedra in silica, TiO6 octahedra in rutile).
#pragma once
#include <array>
#include <cstdint>
#include <string>
#include <vector>

#include "caps/system.hpp"

namespace caps {

// "R", "S" or "" for every atom; `only` (when not empty) limits the centres considered.
std::vector<std::string> stereo_labels(const System& s, const std::vector<char>& only = {});
// "E" or "Z" on both atoms of every stereogenic acyclic double bond (CIP ranks of the two substituents at each end, the
// configuration read off the 3D geometry); bonds, when given, receives the double bonds labelled.
std::vector<std::string> ez_labels(const System& s, std::vector<std::pair<uint32_t, uint32_t>>* bonds = nullptr);

struct Mesh {
  std::vector<Vec3> vertices, normals;
  std::vector<std::array<uint32_t, 3>> triangles;
  std::vector<unsigned> colours;      // per vertex 0xRRGGBB (empty: one colour)
  std::vector<int32_t> nearest;       // per vertex: the nearest atom (for colouring by atom)
  double area() const;                // Å²
};

enum class SurfaceKind { Accessible, VanDerWaals, Excluded };

struct SurfaceOptions {
  SurfaceKind kind = SurfaceKind::Accessible;
  double probe = 1.4;                 // Å (water)
  double spacing = 0.6;               // Å grid
  std::vector<char> atoms;            // which atoms the surface wraps (empty: all)
};

Mesh surface_mesh(const System& s, const SurfaceOptions& o);
std::vector<double> surface_potential(const System& s, const Mesh& m, const std::vector<char>& atoms = {});

// The hull triangles (as a mesh, coloured by the centre's element) of the centres flagged in `centres` (empty: all).
Mesh polyhedra(const System& s, const std::vector<char>& centres = {});

// A smooth tube path through each backbone of the flagged atoms (points every ~0.5 Å), for the ribbon style.
std::vector<std::vector<Vec3>> ribbon_paths(const System& s, const std::vector<char>& atoms);

}  // namespace caps
