// CAPS voids (design/boards/FreeVolume): the empty space of a periodic cell as a distance field — at every grid point
// the distance to the nearest atom's van der Waals surface (Bondi radii) — its local maxima taken as the centres of
// the largest empty spheres (non-overlapping, biggest first), and the share of the cell a probe's centre can reach.
#pragma once
#include <string>
#include <vector>

#include "caps/appearance.hpp"
#include "caps/system.hpp"

namespace caps {

struct VoidSphere {
  Vec3 centre{0, 0, 0};
  double radius = 0;   // Å, to the nearest van der Waals surface
};

struct VoidOptions {
  double grid = 0.5;         // Å between field points
  double probe = 1.4;        // Å: the accessible share is for a probe of this radius
  int max_count = 40;        // spheres kept
  double min_radius = 1.0;   // Å
  double reach = 10.0;       // Å: distances beyond this from every surface count as this
};

struct VoidReport {
  std::vector<VoidSphere> spheres;   // largest first
  double accessible_point = 0;       // share of the cell outside every van der Waals sphere (a point probe)
  double accessible_probe = 0;       // share a probe centre of radius `probe` can reach
  int grid[3] = {0, 0, 0};
  double largest = 0;                // Å
};

VoidReport largest_voids(const System& s, const VoidOptions& o = {});
// Icospheres for drawing (subdivision 2: 320 triangles each).
Mesh void_mesh(const std::vector<VoidSphere>& v, int subdivisions = 2);
// HETATM records named VOI at the centres, the radius in the B-factor column (Å), the cell in CRYST1.
void write_voids_pdb(const System& s, const std::vector<VoidSphere>& v, const std::string& path);

}  // namespace caps
