// CAPS pore network (gas barrier of butyl and halobutyl rubbers, porous fillers): the geometry a spherical probe sees in a
// periodic cell, after Willems, Rycroft, Kazi, Meza & Haranczyk, Microporous Mesoporous Mater. 149, 134 (2012).
//
//  Field      the distance from every grid point to the nearest van der Waals surface (caps/voids.hpp).
//  Clusters   the space a probe's centre can reach (distance > probe radius), labelled into connected clusters across the
//             periodic boundaries (caps/domains.hpp); a cluster that reaches its own image percolates along those axes and
//             is a channel, one that does not is a pocket (inaccessible from outside).
//  Diameters  LCD, the largest included sphere (twice the field's maximum); PLD, the pore-limiting diameter along each
//             axis, twice the largest radius r for which the points farther than r from every surface still percolate
//             (found by bisection); Dif, the largest included sphere along that free path.
//  Surface    the accessible surface area: points spread evenly (a Fibonacci lattice) over every atom's sphere of radius
//             R + probe, kept when no other such sphere covers them, and counted as accessible when they open into a
//             channel and as non-accessible when they open into a pocket; per cell, per volume (m²/cm³), per mass (m²/g).
#pragma once
#include <array>
#include <string>
#include <vector>

#include "caps/json.hpp"
#include "caps/pipeline.hpp"
#include "caps/system.hpp"

namespace caps {

struct PoreNetworkOptions {
  double probe = 1.4;          // Å, probe radius (1.4: a water-sized sphere)
  double grid = 0.25;          // Å between field points
  int samples = 2000;          // surface points per atom
  std::vector<double> radii;   // per atom, Å; empty: Bondi's van der Waals radii
};

struct PoreCluster {
  size_t voxels = 0;
  double volume = 0;                          // Å³ of probe-centre space
  std::array<bool, 3> percolates{false, false, false};
  bool channel = false;                       // percolates along at least one axis
  double largest = 0;                         // Å, diameter of the largest sphere included in it
};

struct PoreNetworkReport {
  double lcd = 0, pld = 0, dif = 0;           // Å; pld the largest over the three axes, dif along that free path
  std::array<double, 3> pld_axis{0, 0, 0};    // Å, 0 where nothing percolates along that axis
  int channel_dims = 0;                       // axes some channel percolates along
  double density = 0;                         // g/cm³
  double av_fraction = 0, nav_fraction = 0;   // probe-centre space in channels and in pockets, share of the cell
  double av = 0, nav = 0;                     // the same, cm³/g
  double geometric_fraction = 0;              // share of the cell outside every van der Waals sphere (a point probe)
  double asa = 0, nasa = 0;                   // Å² per cell: accessible and non-accessible surface
  double asa_volume = 0, nasa_volume = 0;     // m²/cm³
  double asa_mass = 0, nasa_mass = 0;         // m²/g
  std::vector<PoreCluster> clusters;          // largest first
  GridField field;                            // the distance field
  std::vector<std::string> notes;
};

PoreNetworkReport pore_network(const System& s, const PoreNetworkOptions& o = {});
// The report as JSON (without the field): {lcd, pld, dif, pld_axis, channel_dims, density, av_fraction, nav_fraction,
// av, nav, geometric_fraction, asa, nasa, asa_volume, nasa_volume, asa_mass, nasa_mass, clusters [{volume, voxels,
// percolates, channel, largest}], grid, notes}.
Json pore_json(const PoreNetworkReport& r);

}  // namespace caps
