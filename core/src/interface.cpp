// CAPS interfaces: a polymer film grown onto a crystal slab (see polymer.hpp, build_interface).
#include <algorithm>
#include <array>
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

System grow_blend(const std::vector<BlendComponent>& comps, const BlendOptions& o, BlendReport* report) {
  if (comps.empty()) throw GrowError("a blend needs at least one component");
  if (o.density <= 0) throw GrowError("give the growth density");
  if (o.morphology == BlendMorphology::Slabs && comps.size() != 2) throw GrowError("the two-slab start takes exactly two components");
  BlendReport R;
  // mean chain mass per component
  std::vector<double> mass(comps.size());
  for (size_t k = 0; k < comps.size(); ++k) {
    double m = 0;
    for (int q = 0; q < 8; ++q) m += chain_mass(comps[k].spec, chain_sequence(comps[k].spec, o.grow.seed + uint64_t(q) * 101));
    mass[k] = m / 8;
  }
  // chain counts from the weight fractions, the first component setting the scale
  double wsum = 0;
  for (const auto& c : comps) wsum += std::max(0.0, c.weight);
  if (wsum <= 0) throw GrowError("the weight fractions add up to zero");
  const int n0 = comps[0].chains > 0 ? comps[0].chains : std::max(1, o.chains);
  std::vector<int> n(comps.size());
  n[0] = n0;
  for (size_t k = 1; k < comps.size(); ++k)
    n[k] = comps[k].chains > 0 ? comps[k].chains
                               : std::max(1, int(std::lround(double(n0) * mass[0] / mass[k] * comps[k].weight / std::max(1e-12, comps[0].weight))));
  double mtot = 0;
  for (size_t k = 0; k < comps.size(); ++k) mtot += n[k] * mass[k];
  const double vol = mtot / (o.density * 0.602214076);
  std::array<double, 3> L;
  if (o.morphology == BlendMorphology::Slabs) {
    const double a = std::cbrt(vol / 2);
    L = {a, a, 2 * a};
  } else {
    const double a = std::cbrt(vol);
    L = {a, a, a};
  }
  System s;
  bool first = true;
  int64_t mol_before = 0;
  for (size_t k = 0; k < comps.size(); ++k) {
    GrowOptions g = o.grow;
    g.cell = L;
    g.chains = n[k];
    g.auto_scale = true;
    g.substrate = first ? nullptr : &s;
    if (o.morphology == BlendMorphology::Slabs) {   // each component in its own half along z
      g.z_lo = k == 0 ? 0.5 : L[2] / 2 + 0.5;
      g.z_hi = k == 0 ? L[2] / 2 - 0.5 : L[2] - 0.5;
    }
    GrowReport gr;
    s = grow_chains(comps[k].spec, g, &gr);
    int64_t top = 0;
    for (const auto& a : s.atoms) top = std::max(top, a.mol);
    R.molecules.push_back({mol_before + 1, top});
    mol_before = top;
    for (const auto& note : gr.notes)
      if (note.find("contact scale") != std::string::npos) R.notes.push_back(comps[k].spec.units.empty() ? note : comps[k].spec.units[0].name + ": " + note);
    first = false;
  }
  R.chains = n;
  for (size_t k = 0; k < comps.size(); ++k) R.weight_fraction.push_back(n[k] * mass[k] / mtot);
  std::string comp;
  for (size_t k = 0; k < comps.size(); ++k) {
    char b[160];
    std::snprintf(b, sizeof b, "%s%s: %d chains, %.1f wt %% (molecules %lld–%lld)", k ? " · " : "",
                  comps[k].spec.units.empty() ? "component" : comps[k].spec.units[0].name.c_str(), n[k], 100 * R.weight_fraction[k],
                  static_cast<long long>(R.molecules[k].first), static_cast<long long>(R.molecules[k].second));
    comp += b;
  }
  char b[200];
  std::snprintf(b, sizeof b, "blend of %zu components, %s start · cell %.2f × %.2f × %.2f Å · %.3f g/cm³", comps.size(),
                o.morphology == BlendMorphology::Slabs ? "two-slab" : "mixed", L[0], L[1], L[2], o.density);
  R.notes.insert(R.notes.begin(), comp);
  R.notes.insert(R.notes.begin(), b);
  s.title = "CAPS blend: " + comp;
  s.notes = R.notes;
  if (report) *report = R;
  return s;
}

}  // namespace caps
