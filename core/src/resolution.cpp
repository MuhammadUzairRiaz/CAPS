// CAPS model resolution: see caps/resolution.hpp.
#include "caps/resolution.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>

#include "caps/analysis.hpp"
#include "caps/elements.hpp"

namespace caps {

namespace {
std::vector<std::vector<uint32_t>> adjacency(const System& s) {
  std::vector<std::vector<uint32_t>> nb(s.atoms.size());
  for (const auto& b : s.bonds) nb[b.i].push_back(b.j), nb[b.j].push_back(b.i);
  return nb;
}
// A type per label with its mass (the site's mass: a united CH2 weighs C + 2 H).
int type_of(System& out, std::map<std::string, int>& types, const std::string& label, double mass) {
  auto it = types.find(label);
  if (it != types.end()) return it->second;
  TypeInfo t;
  t.type = int(out.types.size()) + 1;
  t.label = label;
  t.mass = mass;
  out.types.push_back(t);
  return types[label] = t.type;
}
}  // namespace

ResolutionReport all_atom_summary(const System& s) {
  ResolutionReport r;
  r.sites = int(s.atoms.size());
  for (const auto& a : s.atoms) {
    r.mass += element(a.element).mass;
    r.hydrogens += a.element == 1 ? 1 : 0;
  }
  return r;
}

System united_atom(const System& s, ResolutionReport* rep) {
  const size_t n = s.atoms.size();
  const auto nb = adjacency(s);
  std::vector<int> host(n, -1);   // the carbon a hydrogen folds into
  for (size_t i = 0; i < n; ++i)
    if (s.atoms[i].element == 1 && nb[i].size() == 1 && s.atoms[nb[i][0]].element == 6) host[i] = int(nb[i][0]);
  System out;
  out.title = s.title.empty() ? "united-atom" : s.title + " (united-atom)";
  out.cell = s.cell;
  out.has_mol = s.has_mol;
  out.unwrapped = s.unwrapped;
  out.bonds_from_file = true;
  std::map<std::string, int> types;
  ResolutionReport r;
  r.site_of.assign(n, -1);
  std::vector<int> hcount(n, 0);
  std::vector<double> hcharge(n, 0.0);
  for (size_t i = 0; i < n; ++i)
    if (host[i] >= 0) ++hcount[size_t(host[i])], hcharge[size_t(host[i])] += s.atoms[i].charge;
  for (size_t i = 0; i < n; ++i) {
    if (host[i] >= 0) continue;
    Atom a = s.atoms[i];
    a.id = int64_t(out.atoms.size() + 1);
    const int h = hcount[i];
    const double mass = element(a.element).mass + h * element(1).mass;
    std::string label = element(a.element).symbol;
    if (a.element == 6 && h > 0) label = h == 1 ? "CH" : "CH" + std::to_string(h);
    a.name = label;
    a.charge += hcharge[i];
    a.type = type_of(out, types, label, mass);
    r.site_of[i] = int(out.atoms.size());
    out.atoms.push_back(a);
    r.mass += mass;
    r.hydrogens += a.element == 1 ? 1 : 0;
  }
  for (size_t i = 0; i < n; ++i)
    if (host[i] >= 0) r.site_of[i] = r.site_of[size_t(host[i])];
  for (const auto& b : s.bonds) {
    if (host[b.i] >= 0 || host[b.j] >= 0) continue;
    out.bonds.push_back({uint32_t(r.site_of[b.i]), uint32_t(r.site_of[b.j]), b.order});
  }
  r.sites = int(out.atoms.size());
  int folded = 0;
  for (size_t i = 0; i < n; ++i) folded += host[i] >= 0 ? 1 : 0;
  r.notes.push_back(std::to_string(folded) + " hydrogens on carbon folded into their carbons; " + std::to_string(r.hydrogens) + " polar hydrogens kept");
  if (rep) *rep = std::move(r);
  return out;
}

System coarse_grain(const System& s, int per_bead, ResolutionReport* rep) {
  if (per_bead < 1) throw std::invalid_argument("coarse_grain: at least one backbone atom per bead");
  const size_t n = s.atoms.size();
  const auto nb = adjacency(s);
  int nm = 0;
  const auto mol = s.molecules(&nm);
  // bead of each backbone atom, then every other atom takes the bead of the nearest backbone atom along bonds
  std::vector<int> bead(n, -1);
  int beads = 0;
  std::vector<int> mol_of_bead;
  std::vector<char> has(size_t(std::max(nm, 1)), 0);
  for (const auto& path : backbones(s, 2)) {
    for (size_t k = 0; k < path.size(); ++k) {
      if (k % size_t(per_bead) == 0) {
        // a short tail joins the bead before it (fewer than half a bead)
        if (k > 0 && path.size() - k < size_t(per_bead + 1) / 2) { bead[path[k]] = beads - 1; continue; }
        mol_of_bead.push_back(mol[path[k]]);
        ++beads;
      }
      bead[path[k]] = beads - 1;
    }
    has[size_t(mol[path[0]])] = 1;
  }
  std::vector<uint32_t> queue;
  for (uint32_t i = 0; i < n; ++i) if (bead[i] >= 0) queue.push_back(i);
  for (size_t h = 0; h < queue.size(); ++h)
    for (uint32_t w : nb[queue[h]])
      if (bead[w] < 0) bead[w] = bead[queue[h]], queue.push_back(w);
  std::vector<int> mol_bead(size_t(std::max(nm, 1)), -1);   // molecules without a backbone: one bead
  for (uint32_t i = 0; i < n; ++i) {
    if (bead[i] >= 0) continue;
    int& b = mol_bead[size_t(mol[i])];
    if (b < 0) b = beads++, mol_of_bead.push_back(mol[i]);
    bead[i] = b;
  }
  // centre of mass per bead, unwrapped about the bead's first atom (minimum image in a periodic cell)
  std::vector<double> mass(size_t(beads), 0.0), charge(size_t(beads), 0.0);
  std::vector<Vec3> com(size_t(beads), Vec3{0, 0, 0});
  std::vector<int> first(size_t(beads), -1);
  const bool periodic = s.cell.valid();
  for (uint32_t i = 0; i < n; ++i) {
    const auto b = size_t(bead[i]);
    if (first[b] < 0) first[b] = int(i);
    Vec3 p = s.atoms[i].pos;
    if (periodic) p = s.atoms[size_t(first[b])].pos + s.cell.minimum_image(p - s.atoms[size_t(first[b])].pos);
    const double m = element(s.atoms[i].element).mass;
    com[b] = com[b] + p * m;
    mass[b] += m;
    charge[b] += s.atoms[i].charge;
  }
  System out;
  out.title = s.title.empty() ? "coarse-grained" : s.title + " (coarse-grained)";
  out.cell = s.cell;
  out.has_mol = true;
  out.bonds_from_file = true;
  std::map<std::string, int> types;
  ResolutionReport r;
  for (int b = 0; b < beads; ++b) {
    Atom a;
    a.id = b + 1;
    a.element = 6;   // drawn as a carbon-sized sphere; the type carries the bead's mass
    a.mol = int64_t(mol_of_bead[size_t(b)] + 1);
    a.pos = com[size_t(b)] * (1.0 / std::max(1e-12, mass[size_t(b)]));
    a.charge = charge[size_t(b)];
    char label[32];
    std::snprintf(label, sizeof label, "B%.0f", std::round(mass[size_t(b)]));
    a.name = label;
    a.type = type_of(out, types, label, mass[size_t(b)]);
    out.atoms.push_back(a);
    r.mass += mass[size_t(b)];
  }
  std::map<std::pair<int, int>, bool> seen;
  for (const auto& bd : s.bonds) {
    int x = bead[bd.i], y = bead[bd.j];
    if (x == y) continue;
    if (x > y) std::swap(x, y);
    if (seen.emplace(std::pair{x, y}, true).second) out.bonds.push_back({uint32_t(x), uint32_t(y), 1});
  }
  r.sites = beads;
  r.site_of = bead;
  r.notes.push_back(std::to_string(per_bead) + " backbone atoms per bead, centre of mass; " + std::to_string(beads) + " beads");
  if (rep) *rep = std::move(r);
  return out;
}

}  // namespace caps
