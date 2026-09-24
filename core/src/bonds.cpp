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
  return out;
}

void make_molecules_whole(System& s) {
  if (!s.cell.valid() || s.bonds.empty()) return;
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
        s.atoms[j].pos = s.atoms[i].pos + s.cell.minimum_image(s.atoms[j].pos - s.atoms[i].pos);
        done[j] = 1;
        stack.push_back(j);
      }
    }
  }
  s.unwrapped = true;
}

}  // namespace caps
