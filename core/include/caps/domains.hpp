// CAPS domains of a density field: phase domains of blends and thermoplastic elastomers, filler-rich regions.
//
//  Field      the mass density of a group of atoms (a selection query; all atoms by default) on a grid spanning the
//             cell, each atom a Gaussian of width sigma normalised on the grid, so the field holds the group's mass.
//  Domains    voxels above (or below) a threshold, labelled into connected regions through the faces of the voxels and
//             across the periodic boundaries; while a region is flooded each voxel keeps the image it was reached in, and
//             a region that reaches one of its own voxels in another image percolates along those axes. Per region: the
//             volume, the equivalent radius (3V/4π)^⅓, the centroid and the semi-axes √(5λ_k) of the ellipsoid with the
//             same gyration tensor (exact for a uniform ellipsoid), and the distance to the nearest other enclosed
//             region's centroid.
//  Field      the histogram of the voxel values, and a slice on any plane by trilinear interpolation (periodic).
#pragma once

#include <array>
#include <string>
#include <vector>

#include "caps/pipeline.hpp"
#include "caps/system.hpp"

namespace caps {

struct DomainOptions {
  std::string group;            // selection query (query.hpp); "" all atoms
  double grid = 0.8;            // Å
  double sigma = 1.5;           // Å
  double threshold = 0;         // g/cm³ (or the field's unit); 0: the field's mean
  bool above = true;            // domains where the field is above the threshold (false: below)
  int min_voxels = 1;           // smaller regions are left out
};

struct Domain {
  size_t voxels = 0;
  double volume = 0;            // Å³
  double equivalent_radius = 0; // Å
  Vec3 centroid{0, 0, 0};       // Å (in the cell)
  Vec3 axes{0, 0, 0};           // fitted-ellipsoid semi-axes, largest first, Å
  std::array<bool, 3> percolates{false, false, false};   // along a, b, c
  bool enclosed = true;         // percolates along no axis
  double nearest = 0;           // Å to the nearest other enclosed region's centroid (0: none)
};

struct DomainResult {
  GridField field;
  std::vector<int> label;       // per voxel: domain index, −1 none
  std::vector<Domain> domains;  // largest first
  double threshold = 0, fraction = 0;   // the level used; the share of the cell in domains
  std::vector<double> hist_x, hist_y;   // voxel values: bin centres, fraction of voxels
  std::vector<std::string> notes;
};

// The group's Gaussian mass density (g/cm³) on a grid of about `grid` Å.
GridField density_field(const System& s, const std::vector<char>& atoms, double grid, double sigma);
// Domains of a field at a level.
DomainResult label_domains(const GridField& g, double threshold, bool above, int min_voxels = 1);
// Both: the group's field and its domains.
DomainResult density_domains(const System& s, const DomainOptions& o);
// The field on a plane through `point` with unit vectors u, v in it: nx × ny samples over width × height, centred on
// the point; row-major (y outer). Periodic trilinear interpolation.
std::vector<double> field_slice(const GridField& g, const Vec3& point, const Vec3& u, const Vec3& v, double width, double height, int nx, int ny);

}  // namespace caps
