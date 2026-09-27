// Exact Voronoi tessellation, cell by cell (as voro++ builds it, Rycroft, Chaos 19, 041111 (2009)): each atom's cell
// starts as a box and is clipped by the bisecting plane of every neighbour — nearest first, until no farther neighbour
// can reach the cell — or, radical (power) tessellation, by the plane where d² − r_i² = d² − r_j² (Gellatly & Finney
// 1982). Periodic cells use every image within the search radius.
#pragma once
#include <vector>

#include "caps/system.hpp"

namespace caps {

struct VoronoiOptions {
  bool radical = false;          // power tessellation weighted by the van der Waals radii (Bondi)
  double face_area_min = 0.0;    // Å²: smaller faces are not counted as faces (their volume stays in the cell)
  double edge_min = 1e-6;        // Å: shorter edges are not counted in a face's edge count
  std::vector<char> only;        // the atoms whose cells to build (empty: every atom); all atoms still bound them
  int threads = 0;               // 0: automatic
};

struct VoronoiCell {
  double volume = 0;             // Å³
  double area = 0;               // Å², the cell's surface
  int faces = 0;                 // faces counted (area ≥ face_area_min): the Voronoi coordination
  int max_face_order = 0;        // the most edges of any face
  std::vector<int> index;        // index[k]: faces with k edges (k = 0 … 12), the Voronoi index ⟨n3, n4, n5, n6, …⟩
  std::vector<int> neighbours;   // the atom behind each counted face
  bool bounded = true;           // false: no cell (an isolated atom without a periodic cell)
};

std::vector<VoronoiCell> voronoi_cells(const System& s, const VoronoiOptions& o = {});

}  // namespace caps
