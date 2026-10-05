#include <algorithm>
#include <cmath>

#include "caps/analysis.hpp"
#include "caps/elements.hpp"
#include "cell_list.hpp"

namespace caps {

std::vector<Bond> perceive_bonds(const System& s, const BondOptions& opt) {
  double rmax = 0;
  for (const auto& a : s.atoms) rmax = std::max(rmax, element(a.element).covalent);
  const double cutoff = 2 * rmax + opt.tolerance;
  Grid g(s, cutoff);
  std::vector<Bond> out;
  for (uint32_t i = 0; i < s.atoms.size(); ++i) {
    const Vec3& f = g.frac[i];
    const double ri = element(s.atoms[i].element).covalent;
    const bool hi = s.atoms[i].element == 1;
    g.for_neighbour_bins(g.bin(f, 0), g.bin(f, 1), g.bin(f, 2), [&](const std::vector<uint32_t>& bin) {
      for (uint32_t j : bin) {
        if (j <= i) continue;
        if (opt.skip_hh && hi && s.atoms[j].element == 1) continue;
        const double lim = ri + element(s.atoms[j].element).covalent + opt.tolerance;
        const Vec3 d = g.sep(i, j);
        const double r2 = dot(d, d);
        if (r2 < lim * lim && r2 > opt.min_distance * opt.min_distance) out.push_back({i, j});
      }
    });
  }
  std::sort(out.begin(), out.end(), [](const Bond& x, const Bond& y) { return x.i != y.i ? x.i < y.i : x.j < y.j; });
  // A long contact (beyond the covalent radii + 0.25 Å) between two atoms bonded to a common third closes a ring across
  // it, it is not a bond: Si···Si over an edge-sharing Si₂O₂ ring of amorphous silica (2.6 Å, a Si–Si bond is 2.35 Å).
  // Cyclopropane's C–C (radii + 0.0 Å) is kept.
  std::vector<std::vector<uint32_t>> nb(s.atoms.size());
  for (const auto& b : out) {
    nb[b.i].push_back(b.j);
    nb[b.j].push_back(b.i);
  }
  std::vector<Bond> kept;
  kept.reserve(out.size());
  for (const auto& b : out) {
    const double lim = element(s.atoms[b.i].element).covalent + element(s.atoms[b.j].element).covalent + 0.25;
    const Vec3 d = g.sep(b.i, b.j);
    bool ring = false;
    if (dot(d, d) > lim * lim)
      for (uint32_t k : nb[b.i])
        if (k != b.j && std::find(nb[b.j].begin(), nb[b.j].end(), k) != nb[b.j].end()) ring = true;
    if (!ring) kept.push_back(b);
  }
  return kept;
}

std::vector<Vec3> whole_positions(const System& s) {
  std::vector<Vec3> p(s.atoms.size());
  for (size_t i = 0; i < s.atoms.size(); ++i) p[i] = s.atoms[i].pos;
  if (!s.cell.valid() || s.bonds.empty()) return p;
  const auto nb = s.neighbours();
  std::vector<char> done(s.atoms.size(), 0);
  std::vector<uint32_t> stack;
  for (uint32_t root = 0; root < s.atoms.size(); ++root) {
    if (done[root]) continue;
    done[root] = 1;
    stack.push_back(root);
    while (!stack.empty()) {
      const uint32_t i = stack.back();
      stack.pop_back();
      for (uint32_t j : nb[i]) {
        if (done[j]) continue;
        p[j] = p[i] + s.cell.minimum_image(p[j] - p[i]);
        done[j] = 1;
        stack.push_back(j);
      }
    }
  }
  return p;
}

void make_molecules_whole(System& s) {
  if (!s.cell.valid() || s.bonds.empty()) return;
  const auto p = whole_positions(s);
  for (size_t i = 0; i < s.atoms.size(); ++i) s.atoms[i].pos = p[i];
  s.unwrapped = true;
}

}  // namespace caps
