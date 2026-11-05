#pragma once
// Lattice planes (hkl) of a cell, for the viewer: the planes h u + k v + l w = m + s (u, v, w fractional, m integer,
// s a shift in units of the spacing) that cross the cell, each cut to the cell as a polygon; their spacing
// d_hkl = 1 / |h a* + k b* + l c*| (a*, b*, c* the reciprocal vectors without 2π: a* = b × c / V …), their unit normal
// along h a* + k b* + l c*, and the camera that looks straight down that normal. The indices are kept as given: (222)
// is spaced half as far as (111), as in diffraction.
#include <array>
#include <string>
#include <vector>

#include "caps/appearance.hpp"
#include "caps/system.hpp"

namespace caps {

struct MillerPlanes {
  std::array<int, 3> hkl{0, 0, 1};
  double d = 0;                                   // Å
  Vec3 normal{0, 0, 1};                           // unit, Cartesian
  std::vector<double> offsets;                    // Å from the cell origin along the normal, one per plane
  std::vector<std::vector<Vec3>> polygons;        // each plane cut to the cell, corners in order (Cartesian, Å)
  double yaw = 0, pitch = 0;                      // radians: the view down the normal (caps::Camera)
  std::string label() const;                      // "(1 1 -1)" with bars written as minus signs
};
// shift: where the planes sit, in units of d (0: through the origin). max_planes: at most this many, those nearest the
// middle of the cell (0: every plane that crosses it). Throws for (000) or a cell that is not periodic.
MillerPlanes miller_planes(const Cell& cell, int h, int k, int l, double shift = 0.0, int max_planes = 0);
// Triangles of the polygons (both faces' normals along the plane normal) for the view.
Mesh miller_mesh(const MillerPlanes& p);
// The camera angles (Camera::yaw, pitch) that look straight down n (n toward the viewer).
std::array<double, 2> view_along(const Vec3& n);
// d_hkl alone.
double d_spacing(const Cell& cell, int h, int k, int l);

}  // namespace caps
