// CAPS layer stacks (see layers.hpp).
#include "caps/layers.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <stdexcept>

namespace caps {

namespace {

bool rectangular(const Cell& c) {
  const double tol = 1e-3;
  return c.valid() && std::fabs(c.a[1]) < tol && std::fabs(c.a[2]) < tol && std::fabs(c.b[0]) < tol && std::fabs(c.b[2]) < tol &&
         std::fabs(c.c[0]) < tol && std::fabs(c.c[1]) < tol && c.a[0] > 0 && c.b[1] > 0 && c.c[2] > 0;
}

// the layer with its molecules whole along z (centres back in the cell) and atoms wrapped in x and y
System prepared(const System& in) {
  System s = in;
  const Vec3 o = s.cell.origin;
  const double L[3] = {s.cell.a[0], s.cell.b[1], s.cell.c[2]};
  const size_t n = s.atoms.size();
  const auto nb = s.neighbours();
  std::vector<char> seen(n, 0);
  for (size_t r = 0; r < n; ++r) {   // z unwrapped along the bonds, one molecule at a time
    if (seen[r]) continue;
    std::vector<size_t> comp{r}, stack{r};
    seen[r] = 1;
    while (!stack.empty()) {
      const size_t u = stack.back();
      stack.pop_back();
      for (uint32_t v : nb[u]) {
        if (seen[v]) continue;
        double dz = s.atoms[v].pos[2] - s.atoms[u].pos[2];
        dz -= L[2] * std::round(dz / L[2]);
        s.atoms[v].pos[2] = s.atoms[u].pos[2] + dz;
        seen[v] = 1;
        stack.push_back(v), comp.push_back(v);
      }
    }
    double zc = 0;
    for (size_t i : comp) zc += s.atoms[i].pos[2];
    zc /= double(comp.size());
    const double shift = -L[2] * std::floor((zc - o[2]) / L[2]);
    for (size_t i : comp) s.atoms[i].pos[2] += shift;
  }
  for (auto& a : s.atoms)
    for (int k = 0; k < 2; ++k) a.pos[size_t(k)] -= L[k] * std::floor((a.pos[size_t(k)] - o[size_t(k)]) / L[k]);
  return s;
}

// repeats n of an edge L to come closest to the target T (the strain T / (n L) − 1 smallest)
int best_repeat(double L, double T, int max_repeat) {
  int best = 1;
  for (int k = 1; k <= max_repeat; ++k)
    if (std::fabs(T / (k * L) - 1) < std::fabs(T / (best * L) - 1)) best = k;
  return best;
}

}  // namespace

System stack_layers(const std::vector<StackLayerInput>& layers, const StackOptions& o, StackReport* report) {
  if (layers.size() < 2) throw std::invalid_argument("a stack needs two layers or more");
  for (const auto& l : layers) {
    if (!l.system || l.system->atoms.empty()) throw std::invalid_argument((l.name.empty() ? "a layer" : l.name) + " has no atoms");
    if (!rectangular(l.system->cell))
      throw std::invalid_argument((l.name.empty() ? "a layer" : l.name) + ": the cell must be rectangular (a along x, b along y, c along z) — make the slab orthogonal");
  }
  StackReport rep;
  const double max_rep = std::max(1, o.max_repeat);
  // the lateral cell: the first layer's, repeated (match both) where that lowers the worst strain of the others
  const Cell& c0 = layers[0].system->cell;
  double best_cost = 1e300;
  int n0[2] = {1, 1};
  for (int dir = 0; dir < 2; ++dir) {
    const double L0 = dir == 0 ? c0.a[0] : c0.b[1];
    best_cost = 1e300;
    const int top = o.match == "both" ? int(max_rep) : 1;
    for (int k = 1; k <= top; ++k) {
      const double T = k * L0;
      if (k > 1 && T > o.max_cell) break;
      double worst = 0;
      for (size_t li = 1; li < layers.size(); ++li) {
        const double L = dir == 0 ? layers[li].system->cell.a[0] : layers[li].system->cell.b[1];
        const int m = best_repeat(L, T, int(max_rep));
        worst = std::max(worst, std::fabs(T / (m * L) - 1));
      }
      if (worst < best_cost - 1e-4) best_cost = worst, n0[dir] = k;   // the smallest repeat that does as well
    }
  }
  const double A = n0[0] * c0.a[0], B = n0[1] * c0.b[1];
  System out;
  out.title = "layer stack";
  std::map<std::pair<std::string, long long>, int> type_of;   // (label, mass × 1e4) → merged type
  int64_t mol_off = 0, id = 0;
  double z = o.vacuum > 0 ? 0.0 : 0.5 * o.gap;
  for (size_t li = 0; li < layers.size(); ++li) {
    const System s = prepared(*layers[li].system);
    const double La = s.cell.a[0], Lb = s.cell.b[1];
    const int na = li == 0 ? n0[0] : best_repeat(La, A, int(max_rep));
    const int nb = li == 0 ? n0[1] : best_repeat(Lb, B, int(max_rep));
    const double fa = A / (na * La), fb = B / (nb * Lb);
    double zlo = 1e300, zhi = -1e300;
    for (const auto& a : s.atoms) zlo = std::min(zlo, a.pos[2]), zhi = std::max(zhi, a.pos[2]);
    const double dz = z - zlo;
    // the layer's types into the merged table
    std::map<int, int> tmap;
    for (const auto& t : s.types) {
      const auto key = std::make_pair(t.label, static_cast<long long>(std::llround(t.mass * 1e4)));
      auto it = type_of.find(key);
      if (it == type_of.end()) {
        const int nt = int(out.types.size()) + 1;
        it = type_of.emplace(key, nt).first;
        out.types.push_back({nt, t.mass, t.label});
      }
      tmap[t.type] = it->second;
    }
    // molecules: renumbered in the layer's order, each copy its own
    std::map<int64_t, int64_t> mol_seen;
    const auto molidx = s.molecules();
    int nmol = 0;
    for (int x : molidx) nmol = std::max(nmol, x + 1);
    const Vec3 org = s.cell.origin;
    for (int ia = 0; ia < na; ++ia)
      for (int ib = 0; ib < nb; ++ib) {
        const uint32_t base = uint32_t(out.atoms.size());
        for (size_t i = 0; i < s.atoms.size(); ++i) {
          Atom a = s.atoms[i];
          const double x = (a.pos[0] - org[0] + ia * La) * fa, y = (a.pos[1] - org[1] + ib * Lb) * fb;
          a.pos = {x, y, a.pos[2] + dz};
          a.image = {0, 0, 0};
          a.mol = mol_off + int64_t(molidx[i]) + 1;
          a.id = ++id;
          if (a.type != 0) { auto t = tmap.find(a.type); a.type = t == tmap.end() ? 0 : t->second; }
          out.atoms.push_back(a);
        }
        for (const auto& b : s.bonds) out.bonds.push_back({b.i + base, b.j + base, b.order});
        mol_off += nmol;
      }
    StackLayerReport lr;
    lr.name = layers[li].name.empty() ? "layer " + std::to_string(li + 1) : layers[li].name;
    lr.na = na, lr.nb = nb;
    lr.strain_a = fa - 1, lr.strain_b = fb - 1;
    lr.z_lo = z, lr.z_hi = z + (zhi - zlo);
    lr.atoms = s.atoms.size() * size_t(na * nb);
    rep.layers.push_back(lr);
    out.has_charges = out.has_charges || s.has_charges;
    z = lr.z_hi + o.gap;
  }
  const double C = o.vacuum > 0 ? z - o.gap + o.vacuum : z - 0.5 * o.gap;
  out.cell.origin = {0, 0, 0};
  out.cell.a = {A, 0, 0};
  out.cell.b = {0, B, 0};
  out.cell.c = {0, 0, C};
  out.has_mol = true;
  out.bonds_from_file = true;
  // bonds across x and y follow the minimum image; the positions stay in the cell there
  for (auto& a : out.atoms)
    for (int k = 0; k < 2; ++k) {
      const double L = k == 0 ? A : B;
      a.pos[size_t(k)] -= L * std::floor(a.pos[size_t(k)] / L);
    }
  rep.a = A, rep.b = B, rep.c = C;
  char b[240];
  for (const auto& l : rep.layers) {
    std::snprintf(b, sizeof b, "%s: %d × %d · strain %+.2f %% × %+.2f %% · z %.1f–%.1f Å · %zu atoms", l.name.c_str(), l.na, l.nb, 100 * l.strain_a, 100 * l.strain_b, l.z_lo, l.z_hi,
                  l.atoms);
    rep.notes.push_back(b);
  }
  std::snprintf(b, sizeof b, "stack of %zu layers: cell %.2f × %.2f × %.2f Å · gap %.1f Å · %s", layers.size(), A, B, C, o.gap,
                o.vacuum > 0 ? ("vacuum " + std::to_string(int(std::lround(o.vacuum))) + " Å on top").c_str() : "periodic in z");
  rep.notes.insert(rep.notes.begin(), b);
  double worst = 0;
  for (const auto& l : rep.layers) worst = std::max({worst, std::fabs(l.strain_a), std::fabs(l.strain_b)});
  if (worst > 0.05) rep.notes.push_back("a layer is strained by more than 5 %: relax (or pick other repeats) before dynamics; crystals should not be strained");
  rep.notes.push_back("the gaps are a start: relax, then equilibrate the stack at constant pressure along z");
  if (report) *report = rep;
  return out;
}

}  // namespace caps
