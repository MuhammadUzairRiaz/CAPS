// CAPS interfaces: a polymer film grown onto a crystal slab (see polymer.hpp, build_interface).
#include <algorithm>
#include <cmath>
#include <cstdio>

#include "caps/polymer.hpp"

namespace caps {

System build_interface(const System& slab, const ChainSpec& spec, const InterfaceOptions& o, GrowReport* report) {
  if (!slab.cell.valid() || slab.atoms.empty()) throw GrowError("the slab has no atoms or no cell");
  const Cell& c = slab.cell;
  if (std::fabs(c.a[1]) > 1e-6 || std::fabs(c.a[2]) > 1e-6 || std::fabs(c.b[0]) > 1e-6 || std::fabs(c.b[2]) > 1e-6 || std::fabs(c.c[0]) > 1e-6 ||
      std::fabs(c.c[1]) > 1e-6)
    throw GrowError("the slab's surface cell is not rectangular; cleave it with the orthogonal option");
  if (o.film <= 2) throw GrowError("the film must be thicker than 2 Å");
  const double Lx = c.a[0], Ly = c.b[1];
  // the slab at the bottom of the cell, wrapped laterally into it
  System sub = slab;
  double zmin = 1e300, zmax = -1e300;
  for (const auto& a : sub.atoms) zmin = std::min(zmin, a.pos[2]), zmax = std::max(zmax, a.pos[2]);
  for (auto& a : sub.atoms) {
    a.pos[0] -= Lx * std::floor((a.pos[0] - c.origin[0]) / Lx);
    a.pos[1] -= Ly * std::floor((a.pos[1] - c.origin[1]) / Ly);
    a.pos[0] -= c.origin[0], a.pos[1] -= c.origin[1];
    a.pos[2] += 0.5 - zmin;
  }
  const double top = zmax - zmin + 0.5;
  const double z_lo = top + std::max(0.0, o.gap), z_hi = z_lo + o.film;
  const double Lz = o.vacuum > 0 ? z_hi + o.vacuum : z_hi + std::max(0.0, o.gap) + 0.5;
  GrowOptions g = o.grow;
  g.cell = {Lx, Ly, Lz};
  g.z_lo = z_lo, g.z_hi = z_hi;
  g.substrate = &sub;
  if (o.chains > 0) g.chains = o.chains;
  else {
    // chains for the film density, from the mean chain mass over a few sequences
    double m = 0;
    for (int k = 0; k < 8; ++k) m += chain_mass(spec, chain_sequence(spec, g.seed + uint64_t(k) * 101));
    m /= 8;
    g.chains = std::max(1, int(std::lround(o.density * Lx * Ly * o.film * 0.602214076 / m)));
  }
  // dense films between hard walls: the contact scale steps down when a chain cannot be placed
  g.auto_scale = o.auto_scale;
  GrowReport rep;
  System s = grow_chains(spec, g, &rep);
  s.title = "CAPS interface: " + (slab.title.empty() ? std::string("slab") : slab.title) + " + " + std::to_string(g.chains) + " chains";
  char b[200];
  std::snprintf(b, sizeof b, "interface: slab %.2f Å thick · film %.1f–%.1f Å (%.1f Å) · %s above · cell %.2f × %.2f × %.2f Å", top - 0.5, z_lo, z_hi, o.film,
                o.vacuum > 0 ? (std::to_string(int(std::lround(o.vacuum))) + " Å vacuum").c_str() : "the slab's periodic image", Lx, Ly, Lz);
  rep.notes.insert(rep.notes.begin(), b);
  s.notes = rep.notes;
  if (report) *report = rep;
  return s;
}

}  // namespace caps
