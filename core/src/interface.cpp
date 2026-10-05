// CAPS interfaces: a polymer film grown onto a crystal slab (see polymer.hpp, build_interface).
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <random>

#include "caps/rng.hpp"
#include "caps/polymer.hpp"
#include "caps/elements.hpp"

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

System build_brush(const System& slab, const ChainSpec& spec, const BrushOptions& o, BrushReport* report) {
  if (!slab.cell.valid() || slab.atoms.empty()) throw GrowError("the slab has no atoms or no cell");
  const Cell& c = slab.cell;
  if (std::fabs(c.a[1]) > 1e-6 || std::fabs(c.a[2]) > 1e-6 || std::fabs(c.b[0]) > 1e-6 || std::fabs(c.b[2]) > 1e-6 || std::fabs(c.c[0]) > 1e-6 ||
      std::fabs(c.c[1]) > 1e-6)
    throw GrowError("the slab's surface cell is not rectangular; cleave it with the orthogonal option");
  if (o.film <= 5) throw GrowError("give the brush more than 5 Å to grow into");
  BrushReport R;
  const double Lx = c.a[0], Ly = c.b[1];
  R.area = Lx * Ly / 100.0;
  const int zs = element_from_symbol(o.site);
  if (zs <= 0) throw GrowError("unknown site element '" + o.site + "'");
  double ztop = -1e300;
  for (const auto& a : slab.atoms) ztop = std::max(ztop, a.pos[2]);
  // sites: the element near the top with a hydrogen pointing up
  const auto nb = slab.neighbours();
  struct Site { uint32_t atom, h; Vec3 dir; };
  std::vector<Site> sites;
  for (uint32_t i = 0; i < slab.atoms.size(); ++i) {
    if (slab.atoms[i].element != zs || slab.atoms[i].pos[2] < ztop - 2.5) continue;
    for (uint32_t h : nb[i])
      if (slab.atoms[h].element == 1) {
        const Vec3 d = slab.cell.minimum_image(slab.atoms[h].pos - slab.atoms[i].pos);
        if (d[2] > 0.3 * norm(d)) { sites.push_back({i, h, d * (1.0 / norm(d))}); break; }
      }
  }
  R.sites = int(sites.size());
  if (sites.empty()) throw GrowError("no graft sites: no " + o.site + " near the top face carries a hydrogen pointing up (hydroxylate the surface first)");
  const int want = o.chains > 0 ? o.chains : std::max(1, int(std::lround(o.density * R.area)));
  // drawn at random, at least min_spacing apart in the surface plane
  std::mt19937_64 rng(o.grow.seed);
  caps::shuffle(sites.begin(), sites.end(), rng);
  std::vector<Site> pick;
  for (const auto& s : sites) {
    if (int(pick.size()) >= want) break;
    bool ok = true;
    for (const auto& p : pick) {
      double dx = slab.atoms[s.atom].pos[0] - slab.atoms[p.atom].pos[0], dy = slab.atoms[s.atom].pos[1] - slab.atoms[p.atom].pos[1];
      dx -= Lx * std::round(dx / Lx), dy -= Ly * std::round(dy / Ly);
      if (dx * dx + dy * dy < o.min_spacing * o.min_spacing) { ok = false; break; }
    }
    if (ok) pick.push_back(s);
  }
  if (int(pick.size()) < want)
    R.notes.push_back("only " + std::to_string(pick.size()) + " of " + std::to_string(want) + " chains: the sites run out at " + std::to_string(o.min_spacing).substr(0, 4) +
                      " Å spacing (" + std::to_string(sites.size()) + " sites)");
  // the slab without the sites' hydrogens, at the bottom of the cell
  std::vector<char> drop(slab.atoms.size(), 0);
  for (const auto& p : pick) drop[p.h] = 1;
  std::vector<int64_t> newi(slab.atoms.size(), -1);
  System sub = slab;
  sub.atoms.clear(), sub.bonds.clear(), sub.velocities.clear();
  double zmin = 1e300;
  for (const auto& a : slab.atoms) zmin = std::min(zmin, a.pos[2]);
  for (size_t i = 0; i < slab.atoms.size(); ++i)
    if (!drop[i]) {
      newi[i] = int64_t(sub.atoms.size());
      Atom a = slab.atoms[i];
      a.pos[0] -= c.origin[0] + Lx * std::floor((a.pos[0] - c.origin[0]) / Lx);
      a.pos[1] -= c.origin[1] + Ly * std::floor((a.pos[1] - c.origin[1]) / Ly);
      a.pos[2] += 0.5 - zmin;
      sub.atoms.push_back(a);
    }
  for (const auto& b : slab.bonds)
    if (newi[b.i] >= 0 && newi[b.j] >= 0) sub.bonds.push_back({uint32_t(newi[b.i]), uint32_t(newi[b.j]), b.order});
  const double top = ztop - zmin + 0.5;
  GrowOptions g = o.grow;
  g.cell = {Lx, Ly, top + o.film + o.vacuum};
  g.z_lo = top - 4.0, g.z_hi = top + o.film;   // the heads sit at their sites, below the top hydrogens; the substrate contacts keep the chains out of the slab
  g.substrate = &sub;
  g.chains = int(pick.size());
  g.anchors.clear();
  for (const auto& p : pick) g.anchors.push_back({uint32_t(newi[p.atom]), {p.dir[0], p.dir[1], p.dir[2]}});
  g.auto_scale = true;
  GrowReport rep;
  System s = grow_chains(spec, g, &rep);
  R.grafted = int(pick.size());
  R.sigma = R.grafted / R.area;
  s.title = "CAPS brush: " + std::to_string(R.grafted) + " chains on " + (slab.title.empty() ? std::string("a surface") : slab.title);
  char b[220];
  std::snprintf(b, sizeof b, "brush: %d chains grafted to %s sites (%d found) · σ = %.3f chains/nm² on %.1f nm² · film %.0f Å, %.0f Å vacuum", R.grafted,
                o.site.c_str(), R.sites, R.sigma, R.area, o.film, o.vacuum);
  R.notes.insert(R.notes.begin(), b);
  for (const auto& n : rep.notes) R.notes.push_back(n);
  s.notes = R.notes;
  if (report) *report = R;
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
  // droplet: the minor component first, inside a sphere holding its share of the volume
  std::vector<size_t> order(comps.size());
  for (size_t k = 0; k < order.size(); ++k) order[k] = k;
  size_t minor = 0;
  if (o.morphology == BlendMorphology::Droplet) {
    for (size_t k = 1; k < comps.size(); ++k)
      if (n[k] * mass[k] < n[minor] * mass[minor]) minor = k;
    std::stable_partition(order.begin(), order.end(), [&](size_t k) { return k == minor; });
  }
  const double rdrop = std::cbrt(3 * n[minor] * mass[minor] / (o.density * 0.602214076) / (4 * 3.14159265358979));
  System s;
  bool first = true;
  int64_t mol_before = 0;
  std::vector<std::pair<int64_t, int64_t>> mols(comps.size());
  for (size_t k : order) {
    GrowOptions g = o.grow;
    g.cell = L;
    g.chains = n[k];
    g.auto_scale = true;
    g.substrate = first ? nullptr : &s;
    if (o.morphology == BlendMorphology::Slabs) {   // each component in its own half along z
      g.z_lo = k == 0 ? 0.5 : L[2] / 2 + 0.5;
      g.z_hi = k == 0 ? L[2] / 2 - 0.5 : L[2] - 0.5;
    }
    if (o.morphology == BlendMorphology::Droplet) {
      g.sphere_radius = rdrop;
      g.sphere_centre = {L[0] / 2, L[1] / 2, L[2] / 2};
      g.sphere_outside = k != minor;
    }
    GrowReport gr;
    s = grow_chains(comps[k].spec, g, &gr);
    int64_t top = 0;
    for (const auto& a : s.atoms) top = std::max(top, a.mol);
    mols[k] = {mol_before + 1, top};
    mol_before = top;
    for (const auto& note : gr.notes)
      if (note.find("contact scale") != std::string::npos) R.notes.push_back(comps[k].spec.units.empty() ? note : comps[k].spec.units[0].name + ": " + note);
    first = false;
  }
  R.molecules = mols;
  if (o.morphology == BlendMorphology::Droplet) {
    char d[160];
    std::snprintf(d, sizeof d, "droplet of %s: radius %.1f Å at the cell centre", comps[minor].spec.units.empty() ? "the minor component" : comps[minor].spec.units[0].name.c_str(), rdrop);
    R.notes.push_back(d);
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
                o.morphology == BlendMorphology::Slabs ? "two-slab" : o.morphology == BlendMorphology::Droplet ? "droplet" : "mixed", L[0], L[1], L[2], o.density);
  R.notes.insert(R.notes.begin(), comp);
  R.notes.insert(R.notes.begin(), b);
  s.title = "CAPS blend: " + comp;
  s.notes = R.notes;
  if (report) *report = R;
  return s;
}

}  // namespace caps
