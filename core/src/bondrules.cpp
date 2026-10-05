#include "caps/bondrules.hpp"

#include <algorithm>
#include <cmath>
#include <map>

#include "caps/elements.hpp"
#include "cell_list.hpp"

namespace caps {

namespace {
std::pair<int, int> key(int a, int b) { return a <= b ? std::pair{a, b} : std::pair{b, a}; }
bool is_ionic(int a, int b) {
  const double ea = pauling_electronegativity(a), eb = pauling_electronegativity(b);
  return ea > 0 && eb > 0 && std::abs(ea - eb) >= 1.7;
}
}  // namespace

std::vector<PairHistogram> pair_histograms(const System& s, double lo, double hi, double bin) {
  const int nb = std::max(1, int(std::ceil((hi - lo) / bin)));
  std::map<std::pair<int, int>, PairHistogram> by;
  std::map<int, int> present;
  for (const auto& a : s.atoms) ++present[a.element];
  for (const auto& [za, _] : present)
    for (const auto& [zb, __] : present)
      if (za <= zb) {
        auto& h = by[{za, zb}];
        h.za = za, h.zb = zb, h.lo = lo, h.bin = bin;
        h.counts.assign(size_t(nb), 0);
        h.covalent = element(za).covalent + element(zb).covalent;
        h.ionic = is_ionic(za, zb);
      }
  Grid g(s, hi);
  for (uint32_t i = 0; i < s.atoms.size(); ++i) {
    const Vec3& f = g.frac[i];
    g.for_neighbour_bins(g.bin(f, 0), g.bin(f, 1), g.bin(f, 2), [&](const std::vector<uint32_t>& bn) {
      for (uint32_t j : bn) {
        if (j <= i) continue;
        const Vec3 d = g.sep(i, j);
        const double r = std::sqrt(dot(d, d));
        if (r < lo || r >= hi) continue;
        ++by[key(s.atoms[i].element, s.atoms[j].element)].counts[size_t((r - lo) / bin)];
      }
    });
  }
  for (const auto& b : s.bonds)
    if (b.i < s.atoms.size() && b.j < s.atoms.size()) ++by[key(s.atoms[b.i].element, s.atoms[b.j].element)].bonded;
  std::vector<PairHistogram> out;
  for (auto& [k, h] : by) {
    // the gap after the bonded distances: the first empty stretch at least 0.16 Å wide (a single empty bin between an
    // aromatic 1.39 Å and a single 1.54 Å C–C is not one) below 2.6 Å; the cut-off 0.1 Å into it
    h.suggested = h.covalent + 0.45;
    int first = -1;
    for (int t = 0; t < nb; ++t) if (h.counts[size_t(t)] > 0) { first = t; break; }
    const int wide = std::max(1, int(std::ceil(0.16 / bin - 1e-9)));
    if (first >= 0 && lo + first * bin < h.covalent + 0.5)
      for (int t = first + 1; t + wide <= nb && lo + t * bin < 2.6; ++t) {
        bool empty = true;
        for (int u = t; u < t + wide; ++u) empty = empty && h.counts[size_t(u)] == 0;
        if (empty) { h.suggested = std::min(lo + t * bin + 0.1, h.covalent + 0.45); break; }
      }
    out.push_back(h);
  }
  // the pairs the structure has most of first
  std::sort(out.begin(), out.end(), [](const PairHistogram& a, const PairHistogram& b) {
    long ca = 0, cb = 0;
    for (int c : a.counts) ca += c;
    for (int c : b.counts) cb += c;
    return a.bonded != b.bonded ? a.bonded > b.bonded : ca > cb;
  });
  return out;
}

std::vector<Bond> bonds_by_rules(const System& s, const std::vector<BondRule>& rules, double min_distance) {
  std::map<std::pair<int, int>, BondRule> rule;
  double reach = 0;
  for (const auto& r : rules) { rule[key(r.za, r.zb)] = r; if (!r.never) reach = std::max(reach, r.max); }
  double rmax = 0;
  for (const auto& a : s.atoms) rmax = std::max(rmax, element(a.element).covalent);
  reach = std::max(reach, 2 * rmax + 0.45);
  Grid g(s, reach);
  std::vector<Bond> out;
  for (uint32_t i = 0; i < s.atoms.size(); ++i) {
    const Vec3& f = g.frac[i];
    g.for_neighbour_bins(g.bin(f, 0), g.bin(f, 1), g.bin(f, 2), [&](const std::vector<uint32_t>& bn) {
      for (uint32_t j : bn) {
        if (j <= i) continue;
        const int za = s.atoms[i].element, zb = s.atoms[j].element;
        double lim = element(za).covalent + element(zb).covalent + 0.45;
        if (auto it = rule.find(key(za, zb)); it != rule.end()) {
          if (it->second.never) continue;
          lim = it->second.max;
        }
        const Vec3 d = g.sep(i, j);
        const double r2 = dot(d, d);
        if (r2 <= lim * lim && r2 > min_distance * min_distance) out.push_back({i, j});
      }
    });
  }
  std::sort(out.begin(), out.end(), [](const Bond& x, const Bond& y) { return x.i != y.i ? x.i < y.i : x.j < y.j; });
  return out;
}

}  // namespace caps
