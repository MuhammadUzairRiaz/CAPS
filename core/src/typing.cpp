// Atom typing: chemical perception and SMARTS rules (see caps/typing.hpp).
#include "caps/typing.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <functional>
#include <deque>
#include <map>
#include <memory>
#include <set>

#include "caps/elements.hpp"

namespace caps {

// ---------------------------------------------------------------------------------------------------------------------
// Perception

int Perception::bond_order(uint32_t a, uint32_t b) const {
  for (size_t k = 0; k < nb[a].size(); ++k)
    if (nb[a][k] == b) return order[a][k];
  return 0;
}

namespace {

// Valence an atom needs (sum of bond orders) for a neutral closed-shell structure, from its element and degree.
// Hypervalent S / P / Se take the valence their number of neighbours implies (sulfoxide 4, sulfone 6, phosphate 5).
int target_valence(int z, int degree) {
  switch (z) {
    case 1: case 9: case 17: case 35: case 53: return 1;
    case 5: case 13: return 3;
    case 6: case 14: case 32: return 4;
    case 7: case 15: case 33:
      if (z == 7) return degree >= 4 ? 4 : 3;
      return degree <= 3 ? 3 : 5;
    case 8: return 2;
    case 16: case 34: return degree <= 2 ? 2 : degree == 3 ? 4 : 6;
    default: return -1;   // metals, noble gases: no covalent valence rule
  }
}

// Formal charge of an isolated ion from its group.
int ion_charge(int z) {
  switch (z) {
    case 3: case 11: case 19: case 37: case 55: return 1;
    case 4: case 12: case 20: case 38: case 56: case 30: return 2;
    case 13: return 3;
    case 9: case 17: case 35: case 53: return -1;
    case 8: case 16: return -2;
    default: return 0;
  }
}

// Smallest set of smallest rings: the shortest cycle through each bond, then a cycle basis over GF(2) of those,
// smallest first.
std::vector<std::vector<uint32_t>> sssr(const std::vector<std::vector<uint32_t>>& nb) {
  const size_t n = nb.size();
  std::map<std::pair<uint32_t, uint32_t>, size_t> edge_id;
  for (uint32_t i = 0; i < n; ++i)
    for (uint32_t j : nb[i])
      if (i < j) edge_id.emplace(std::make_pair(i, j), edge_id.size());
  // components for the cycle rank E − V + C
  std::vector<int> comp(n, -1);
  int ncomp = 0;
  for (uint32_t s = 0; s < n; ++s) {
    if (comp[s] >= 0) continue;
    std::vector<uint32_t> st = {s};
    comp[s] = ncomp;
    while (!st.empty()) {
      uint32_t a = st.back();
      st.pop_back();
      for (uint32_t b : nb[a])
        if (comp[b] < 0) { comp[b] = ncomp; st.push_back(b); }
    }
    ++ncomp;
  }
  const long rank = long(edge_id.size()) - long(n) + ncomp;
  if (rank <= 0) return {};
  // shortest cycle through each edge (i, j): BFS from j to i without the edge
  std::vector<std::vector<uint32_t>> cand;
  std::set<std::vector<uint32_t>> seen;
  for (const auto& [e, id] : edge_id) {
    (void)id;
    const uint32_t i = e.first, j = e.second;
    std::vector<int64_t> prev(n, -2);
    std::vector<uint32_t> q = {j};
    prev[j] = -1;
    bool found = false;
    for (size_t h = 0; h < q.size() && !found; ++h) {
      const uint32_t a = q[h];
      for (uint32_t b : nb[a]) {
        if (a == j && b == i) continue;
        if (prev[b] != -2) continue;
        prev[b] = a;
        if (b == i) { found = true; break; }
        q.push_back(b);
      }
    }
    if (!found) continue;
    std::vector<uint32_t> cyc;
    for (int64_t a = i; a != -1; a = prev[a]) cyc.push_back(uint32_t(a));
    if (cyc.size() > 24) continue;
    auto key = cyc;
    std::sort(key.begin(), key.end());
    if (seen.insert(key).second) cand.push_back(cyc);
  }
  // also the shortest cycles through each atom pair of fused systems are covered by the edge cycles for the rings
  // chemists care about; pick an independent basis, smallest first
  std::stable_sort(cand.begin(), cand.end(), [](const auto& a, const auto& b) { return a.size() < b.size(); });
  std::vector<std::vector<char>> basis;   // reduced edge vectors
  std::vector<size_t> pivot;
  std::vector<std::vector<uint32_t>> out;
  const size_t ne = edge_id.size();
  for (const auto& c : cand) {
    std::vector<char> v(ne, 0);
    for (size_t k = 0; k < c.size(); ++k) {
      uint32_t a = c[k], b = c[(k + 1) % c.size()];
      v[edge_id.at({std::min(a, b), std::max(a, b)})] ^= 1;
    }
    for (size_t r = 0; r < basis.size(); ++r)
      if (v[pivot[r]])
        for (size_t x = 0; x < ne; ++x) v[x] ^= basis[r][x];
    auto it = std::find(v.begin(), v.end(), 1);
    if (it == v.end()) continue;
    pivot.push_back(size_t(it - v.begin()));
    basis.push_back(v);
    out.push_back(c);
    if (long(out.size()) == rank) break;
  }
  return out;
}

}  // namespace

int united_atom_hydrogens(const std::string& name, int element) {
  if (element == 16) return name == "SH1E" || name == "SH" ? 1 : -1;
  if (element == 8) return name == "OH2" ? 2 : name == "OH" ? 1 : -1;   // a one-site water (mW) or hydroxyl site
  if (element != 6 || name.size() < 2 || name[0] != 'C') return -1;
  std::string x = name;
  if (x.size() > 2 && (x.back() == 'E' || x.back() == 'p')) x.pop_back();   // CHARMM19 CH2E, GROMOS CH3p
  if (x == "CH") return 1;
  if (x == "CR1") return 1;                                                    // aromatic CH (GROMOS, CHARMM19)
  if (x.size() == 3 && x[1] == 'H' && x[2] >= '0' && x[2] <= '4') return x[2] - '0';
  return -1;
}

Perception perceive(const System& s) {
  const size_t n = s.atoms.size();
  Perception p;
  p.nb.assign(n, {});
  p.order.assign(n, {});
  p.arom_bond.assign(n, {});
  p.aromatic.assign(n, false);
  p.charge.assign(n, 0);
  p.hcount.assign(n, 0);
  // bonds: known orders are kept, the others (unknown, aromatic, amide) are solved from valences
  bool any_order = false;
  for (const auto& b : s.bonds) any_order = any_order || b.order > 0;
  struct E { uint32_t i, j; int order; bool fixed; bool arom_in; };
  std::vector<E> edges;
  for (const auto& b : s.bonds) {
    if (b.i == b.j || b.i >= n || b.j >= n) continue;
    if (b.order == kBondDative) { edges.push_back({b.i, b.j, 0, true, false}); continue; }   // a coordinate bond: no valence
    const bool fixed = any_order && (b.order == 1 || b.order == 2 || b.order == 3 || b.order == 5);
    edges.push_back({b.i, b.j, fixed ? (b.order == 5 ? 1 : b.order) : 1, fixed, b.order == 4});
  }
  std::vector<std::vector<size_t>> inc(n);
  for (size_t e = 0; e < edges.size(); ++e) { inc[edges[e].i].push_back(e); inc[edges[e].j].push_back(e); }
  std::vector<int> z(n), deg(n);
  for (size_t i = 0; i < n; ++i) { z[i] = s.atoms[i].element; deg[i] = int(inc[i].size()); }
  // united-atom sites: only when no carbon has an explicit hydrogen (a heavy-atom PDB's "CH2" is a ring carbon's name)
  p.implicit_h.assign(n, 0);
  {
    bool ch = false, named = false;
    for (const auto& e : edges) ch = ch || (z[e.i] == 6 && z[e.j] == 1) || (z[e.j] == 6 && z[e.i] == 1);
    std::vector<char> has_h(n, 0);   // an atom with an explicit hydrogen never gets implicit ones (an all-atom water "OH2")
    for (const auto& e : edges) {
      if (z[e.j] == 1) has_h[e.i] = 1;
      if (z[e.i] == 1) has_h[e.j] = 1;
    }
    if (!ch)
      for (size_t i = 0; i < n; ++i) {
        const int h = has_h[i] ? -1 : united_atom_hydrogens(s.atoms[i].name, z[i]);
        if (h >= 0) p.implicit_h[i] = h, named = true;
      }
    p.united_atom = named;
  }
  // charged centres the valences alone cannot tell: nitro / N-oxide N (N+ with O−), four-connected N+
  std::vector<int> target(n, -1);
  for (size_t i = 0; i < n; ++i) {
    target[i] = target_valence(z[i], deg[i] + p.implicit_h[i]);
    if (target[i] >= 0) target[i] = std::max(0, target[i] - p.implicit_h[i]);   // implicit hydrogens: single bonds already there
    const bool lone = deg[i] == 0 && p.implicit_h[i] == 0;   // a united-atom methane is not an ion
    if (lone && target[i] >= 0 && z[i] != 1) { target[i] = -1; p.charge[i] = ion_charge(z[i]); }
    if (lone && target[i] < 0) p.charge[i] = ion_charge(z[i]);
    if (z[i] == 7 && deg[i] == 4) p.charge[i] = 1;
    if (z[i] == 7 && deg[i] == 3) {
      int term_o = 0;
      for (size_t e : inc[i]) {
        const uint32_t o = edges[e].i == i ? edges[e].j : edges[e].i;
        if (z[o] == 8 && deg[o] == 1) ++term_o;
      }
      if (term_o >= 2) { target[i] = 4; p.charge[i] = 1; }   // nitro, nitrate: N+(=O)O−
    }
    // azide R-N=N+=N−: the middle N of an N-N-N chain ending in a one-connected N
    if (z[i] == 7 && deg[i] == 2) {
      int term_n = 0, n_nb = 0;
      for (size_t e : inc[i]) {
        const uint32_t o = edges[e].i == i ? edges[e].j : edges[e].i;
        if (z[o] == 7) ++n_nb;
        if (z[o] == 7 && deg[o] == 1) ++term_n;
      }
      if (term_n == 1 && n_nb == 2) { target[i] = 4; p.charge[i] = 1; }
    }
    // oxonium (hydronium, protonated ethers): three-connected O+
    if (z[i] == 8 && deg[i] == 3) { target[i] = 3; p.charge[i] = 1; }
  }
  // free valence left after one bond order per bond (and the fixed extra orders)
  std::vector<int> fv(n, 0);
  for (size_t i = 0; i < n; ++i) {
    if (target[i] < 0) continue;
    int used = 0;
    for (size_t e : inc[i]) used += edges[e].order;
    fv[i] = std::max(0, target[i] - used);
  }
  auto flexible = [&](size_t e) { return !edges[e].fixed; };
  auto other = [&](size_t e, size_t a) { return edges[e].i == a ? edges[e].j : edges[e].i; };
  // raise bond orders until no atom has free valence: forced moves first, then a search over the rest (Kekulé
  // structures of rings), keeping the assignment that satisfies the most atoms
  std::function<void()> propagate = [&] {
    for (bool changed = true; changed;) {
      changed = false;
      for (size_t a = 0; a < n; ++a) {
        if (fv[a] <= 0) continue;
        std::vector<size_t> opts;
        for (size_t e : inc[a])
          if (flexible(e) && fv[other(e, a)] > 0 && edges[e].order < 3) opts.push_back(e);
        if (opts.size() == 1) {
          const size_t e = opts[0];
          const size_t b = other(e, a);
          const int add = std::min({fv[a], fv[b], 3 - edges[e].order});
          edges[e].order += add;
          fv[a] -= add;
          fv[b] -= add;
          changed = true;
        } else if (opts.empty()) {
          continue;
        }
      }
    }
  };
  // forced moves, then a depth-first search over the remaining choices (Kekulé structures of rings), keeping the
  // assignment that leaves the least free valence
  auto solve = [&] {
    for (auto& e : edges)
      if (!e.fixed) e.order = 1;
    for (size_t i = 0; i < n; ++i) {
      fv[i] = 0;
      if (target[i] < 0) continue;
      int used = 0;
      for (size_t e : inc[i]) used += edges[e].order;
      fv[i] = std::max(0, target[i] - used);
    }
    propagate();
    int best_left = 1 << 30;
    std::vector<int> best_order, best_fv;
    long budget = 200000;
    int dropped = 0;   // valence left unsatisfied on purpose along the current branch
    std::function<void()> search = [&] {
      if (--budget < 0) return;
      size_t a = n;
      for (size_t i = 0; i < n && a == n; ++i)
        if (fv[i] > 0)
          for (size_t e : inc[i])
            if (flexible(e) && fv[other(e, i)] > 0 && edges[e].order < 3) { a = i; break; }
      if (a == n) {
        int left = dropped;
        for (size_t i = 0; i < n; ++i) left += fv[i];
        if (left < best_left) {
          best_left = left;
          best_order.clear();
          for (const auto& x : edges) best_order.push_back(x.order);
          best_fv = fv;
        }
        return;
      }
      for (size_t e : inc[a]) {
        const size_t b = other(e, a);
        if (!flexible(e) || fv[b] <= 0 || edges[e].order >= 3) continue;
        std::vector<int> so;
        for (const auto& x : edges) so.push_back(x.order);
        const auto sf = fv;
        edges[e].order += 1;
        fv[a] -= 1;
        fv[b] -= 1;
        propagate();
        search();
        for (size_t k = 0; k < edges.size(); ++k) edges[k].order = so[k];
        fv = sf;
        if (best_left == 0) return;
      }
      // leave this atom's valence unsatisfied (a charge or radical) as the last option
      const int keep = fv[a];
      fv[a] = 0;
      dropped += keep;
      search();
      dropped -= keep;
      fv[a] = keep;
    };
    search();
    if (!best_order.empty()) {
      for (size_t k = 0; k < edges.size(); ++k) edges[k].order = best_order[k];
      fv = best_fv;
    }
    if (budget < 0) p.notes.push_back("bond-order search stopped early; some bond orders may be wrong");
  };
  solve();
  // Bond orders from the file are trusted only when they are consistent: a file that writes some double bonds but
  // leaves others as single (CHARMM topologies write aromatic rings with BOND) leaves atoms short of valence; its
  // single bonds are then solved again (its double and triple bonds kept).
  auto short_valence = [&] {
    int l = 0;
    for (size_t i = 0; i < n; ++i) {
      if (target[i] < 0) continue;
      int used = 0;
      for (size_t e : inc[i]) used += edges[e].order;
      l += std::max(0, target[i] - used);
    }
    return l;
  };
  if (any_order && short_valence() > 0) {
    const int before = short_valence();
    std::vector<E> kept = edges;
    for (auto& e : edges)
      if (e.fixed && e.order == 1) e.fixed = false;
    solve();
    if (short_valence() < before) p.notes.push_back("bond orders in the file leave atoms short of valence; single bonds re-solved");
    else { edges = kept; solve(); }
  }
  // Cations the valences alone cannot place (iminium, amidinium, guanidinium, imidazolium, protonated Schiff bases):
  // when valence is left unsatisfied, try each three-connected neutral N next to an sp2-capable carbon as N+ (four
  // bonds) and keep the one that leaves the least; repeat for polycations.
  auto leftover = [&] {
    int l = 0;
    for (size_t i = 0; i < n; ++i) {
      if (target[i] < 0) continue;
      int used = 0;
      for (size_t e : inc[i]) used += edges[e].order;
      l += std::max(0, target[i] - used);
    }
    return l;
  };
  for (int round = 0; round < 4 && leftover() > 0; ++round) {
    const int before = leftover();
    int best = before;
    size_t best_n = n;
    for (size_t c = 0; c < n; ++c) {
      if (z[c] != 7 || deg[c] != 3 || p.charge[c] != 0 || target[c] != 3) continue;
      bool next_to_sp2 = false;
      for (size_t e : inc[c]) {
        const size_t o = other(e, c);
        next_to_sp2 = next_to_sp2 || (z[o] == 6 && deg[o] == 3 && flexible(e));
      }
      if (!next_to_sp2) continue;
      target[c] = 4;
      solve();
      const int l = leftover();
      if (l < best) { best = l; best_n = c; }
      target[c] = 3;
    }
    if (best_n == n) { solve(); break; }
    target[best_n] = 4;
    p.charge[best_n] = 1;
    solve();
  }
  // what is still free becomes a formal charge on O / S / N / C (carboxylates, thiolates, amide anions, carbanions)
  for (size_t i = 0; i < n; ++i) {
    int left = 0;
    if (target[i] >= 0) {
      int used = 0;
      for (size_t e : inc[i]) used += edges[e].order;
      left = target[i] - used;
    }
    if (left > 0 && (z[i] == 8 || z[i] == 16 || z[i] == 7 || z[i] == 6 || z[i] == 34)) p.charge[i] -= left;
  }
  // formal charges given in a mol2 file win over the valence estimate only where the file is explicit (no field for
  // formal charges in mol2): nothing to do
  for (size_t i = 0; i < n; ++i)
    for (size_t e : inc[i]) {
      const uint32_t o = uint32_t(other(e, i));
      p.nb[i].push_back(o);
      p.order[i].push_back(edges[e].order);
      if (z[o] == 1) ++p.hcount[i];
    }
  for (size_t i = 0; i < n; ++i) p.hcount[i] += p.implicit_h[i];
  // rings
  p.rings = sssr(p.nb);
  p.ring_count.assign(n, 0);
  p.smallest_ring.assign(n, 0);
  p.ring_bonds.assign(n, 0);
  p.ring_bond.assign(n, {});
  for (size_t i = 0; i < n; ++i) { p.ring_bond[i].assign(p.nb[i].size(), false); p.arom_bond[i].assign(p.nb[i].size(), false); }
  std::set<std::pair<uint32_t, uint32_t>> rbonds;
  for (const auto& r : p.rings) {
    for (size_t k = 0; k < r.size(); ++k) {
      const uint32_t a = r[k], b = r[(k + 1) % r.size()];
      rbonds.insert({std::min(a, b), std::max(a, b)});
      ++p.ring_count[a];
      if (p.smallest_ring[a] == 0 || int(r.size()) < p.smallest_ring[a]) p.smallest_ring[a] = int(r.size());
    }
  }
  for (size_t i = 0; i < n; ++i)
    for (size_t k = 0; k < p.nb[i].size(); ++k) {
      const uint32_t o = p.nb[i][k];
      if (rbonds.count({std::min<uint32_t>(uint32_t(i), o), std::max<uint32_t>(uint32_t(i), o)})) {
        p.ring_bond[i][k] = true;
        ++p.ring_bonds[i];
      }
    }
  // aromaticity (Hückel 4n+2 over each SSSR ring, then over pairs of fused rings): a ring atom gives 1 π electron
  // for a double bond in a ring, 0 for an exocyclic double bond to an electronegative atom, 2 for a lone pair
  // (N, P with three single bonds; O, S, Se with two; carbanions), and breaks the ring otherwise
  auto in_set = [](const std::vector<uint32_t>& r, uint32_t a) { return std::find(r.begin(), r.end(), a) != r.end(); };
  auto electrons = [&](const std::vector<uint32_t>& r, uint32_t a) -> int {
    int dbl_ring = 0, dbl_exo = 0, dbl_exo_c = 0;
    for (size_t k = 0; k < p.nb[a].size(); ++k) {
      if (p.order[a][k] != 2) continue;
      if (p.ring_bond[a][k]) ++dbl_ring;
      else if (z[p.nb[a][k]] == 6) ++dbl_exo_c;
      else ++dbl_exo;
    }
    (void)r;
    if (dbl_ring == 1) return 1;
    if (dbl_ring > 1) return -1;
    for (size_t k = 0; k < p.nb[a].size(); ++k)
      if (p.order[a][k] == 3) return -1;
    if (dbl_exo == 1 && (z[a] == 6 || z[a] == 7 || z[a] == 15 || z[a] == 16)) return 0;
    if (dbl_exo_c) return -1;
    const int d = int(p.nb[a].size());
    if ((z[a] == 7 || z[a] == 15) && d == 3 && p.charge[a] == 0) return 2;
    if ((z[a] == 7 || z[a] == 15) && d == 2 && p.charge[a] == -1) return 2;
    if ((z[a] == 8 || z[a] == 16 || z[a] == 34) && d == 2 && p.charge[a] == 0) return 2;
    if (z[a] == 6 && p.charge[a] == -1 && d == 3) return 2;
    if (z[a] == 6 && p.charge[a] == 1 && d == 3) return 0;
    if (z[a] == 5 && d == 3) return 0;
    return -1;
  };
  auto huckel = [&](const std::vector<uint32_t>& r) {
    int pi = 0;
    for (uint32_t a : r) {
      const int e = electrons(r, a);
      if (e < 0) return false;
      pi += e;
    }
    return pi >= 2 && (pi - 2) % 4 == 0;
  };
  std::vector<std::vector<uint32_t>> arom_rings;
  std::vector<bool> ring_arom(p.rings.size(), false);
  for (size_t r = 0; r < p.rings.size(); ++r)
    if (p.rings[r].size() >= 5 && p.rings[r].size() <= 7 && huckel(p.rings[r])) { ring_arom[r] = true; arom_rings.push_back(p.rings[r]); }
  // fused pairs sharing one bond: the envelope ring (e.g. azulene-like or Kekulé forms that leave one ring alone)
  for (size_t r1 = 0; r1 < p.rings.size(); ++r1)
    for (size_t r2 = r1 + 1; r2 < p.rings.size(); ++r2) {
      if (ring_arom[r1] && ring_arom[r2]) continue;
      std::vector<uint32_t> shared;
      for (uint32_t a : p.rings[r1])
        if (in_set(p.rings[r2], a)) shared.push_back(a);
      if (shared.size() != 2) continue;
      std::vector<uint32_t> env;
      for (uint32_t a : p.rings[r1]) env.push_back(a);
      for (uint32_t a : p.rings[r2])
        if (!in_set(env, a)) env.push_back(a);
      if (huckel(env)) {
        ring_arom[r1] = ring_arom[r2] = true;
        arom_rings.push_back(p.rings[r1]);
        arom_rings.push_back(p.rings[r2]);
      }
    }
  for (const auto& r : arom_rings)
    for (size_t k = 0; k < r.size(); ++k) {
      const uint32_t a = r[k], b = r[(k + 1) % r.size()];
      p.aromatic[a] = p.aromatic[b] = true;
      for (size_t q = 0; q < p.nb[a].size(); ++q)
        if (p.nb[a][q] == b) p.arom_bond[a][q] = true;
      for (size_t q = 0; q < p.nb[b].size(); ++q)
        if (p.nb[b][q] == a) p.arom_bond[b][q] = true;
    }
  // aromatic bonds read from the file (mol2 "ar") where perception found no aromatic ring: keep the file's word
  for (const auto& e : edges)
    if (e.arom_in && !(p.aromatic[e.i] && p.aromatic[e.j])) {
      p.notes.push_back("bond " + std::to_string(e.i + 1) + "-" + std::to_string(e.j + 1) +
                        " is aromatic in the file but not by Hückel's rule; typed as aromatic");
      p.aromatic[e.i] = p.aromatic[e.j] = true;
      for (size_t q = 0; q < p.nb[e.i].size(); ++q)
        if (p.nb[e.i][q] == e.j) p.arom_bond[e.i][q] = true;
      for (size_t q = 0; q < p.nb[e.j].size(); ++q)
        if (p.nb[e.j][q] == e.i) p.arom_bond[e.j][q] = true;
    }
  // antechamber's ring classes (GAFF typing): AR1 pure aromatic (six-membered, C / N, aromatic by Hückel, no
  // exocyclic double bond); AR2 planar ring with two consecutive single bonds and at least two double bonds (pyrrole,
  // furan, cyclopentadiene); AR3 planar ring with an exocyclic double bond (pyridone, quinone); AR5 ring of sp3 carbons
  // only; AR4 any other ring. An atom takes the lowest class of its rings.
  p.ar_class.assign(n, 0);
  for (const auto& r : p.rings) {
    const size_t m = r.size();
    auto ring_order = [&](size_t k) { return p.bond_order(r[k], r[(k + 1) % m]); };
    bool planar = true, all_cn = true, all_sp3c = true, exo_double = false, all_arom = true;
    for (uint32_t a : r) {
      const int d = int(p.nb[a].size());
      const int e = z[a];
      bool pl = (e == 6 && d == 3) || (e == 7 && (d == 2 || d == 3)) || ((e == 8 || e == 16 || e == 34) && d == 2) ||
                (e == 15 && (d == 2 || d == 3)) || (e == 5 && d == 3);
      planar = planar && pl;
      all_cn = all_cn && (e == 6 || e == 7);
      all_sp3c = all_sp3c && e == 6 && d == 4;
      all_arom = all_arom && p.aromatic[a];
      for (size_t q = 0; q < p.nb[a].size(); ++q)
        if (p.order[a][q] == 2 && std::find(r.begin(), r.end(), p.nb[a][q]) == r.end()) exo_double = true;
    }
    for (size_t k = 0; k < m && all_arom; ++k)
      for (size_t q = 0; q < p.nb[r[k]].size(); ++q)
        if (p.nb[r[k]][q] == r[(k + 1) % m] && !p.arom_bond[r[k]][q]) all_arom = false;
    int cls = 4;
    if (all_sp3c) cls = 5;
    else if (planar && all_arom && all_cn && m == 6 && !exo_double) cls = 1;
    else if (planar && exo_double) cls = 3;
    else if (planar) {
      int doubles = 0;
      bool two_singles = false;
      for (size_t k = 0; k < m; ++k) {
        if (ring_order(k) == 2) ++doubles;
        if (ring_order(k) == 1 && ring_order((k + 1) % m) == 1) two_singles = true;
      }
      if (doubles >= 2 && two_singles) cls = 2;
      else if (all_arom) cls = 2;   // aromatic by Hückel but not a pure six-membered C / N ring
    }
    for (uint32_t a : r)
      if (p.ar_class[a] == 0 || cls < p.ar_class[a]) p.ar_class[a] = cls;
  }
  return p;
}

// ---------------------------------------------------------------------------------------------------------------------
// SMARTS

class SmartsNode;

namespace {

enum class Op { And, Or, AndLow, Not, Prim };
enum class Prim { True, Aromatic, Aliphatic, Element, Degree, Connect, HCount, ImplicitH, Valence, RingCount, RingSize,
                  RingBonds, Charge, Recursive, Type, Chiral, ArClass };

struct AtomExpr {
  Op op = Op::Prim;
  Prim prim = Prim::True;
  int value = 0;         // element / count (−1: "any non-zero" for R r x)
  int arom = -1;         // for Element: 1 aromatic, 0 aliphatic, −1 either
  std::string name;      // %type
  std::shared_ptr<SmartsNode> rec;
  std::vector<std::unique_ptr<AtomExpr>> kids;
};

enum class BPrim { Default, Single, Double, Triple, Aromatic, Any, Ring, Up };
struct BondExpr {
  Op op = Op::Prim;
  BPrim prim = BPrim::Default;
  std::vector<std::unique_ptr<BondExpr>> kids;
};

}  // namespace

class SmartsNode {
 public:
  struct PBond { int a, b; std::unique_ptr<BondExpr> expr; };
  std::vector<std::unique_ptr<AtomExpr>> atoms;
  std::vector<PBond> bonds;
  std::vector<std::vector<std::pair<int, int>>> adj;   // (other atom, bond index)
  bool uses_types = false;
};

namespace {

struct Parser {
  const std::string& t;
  size_t i = 0;
  size_t bracket_start = std::string::npos;   // position after the current '[': H there is a hydrogen atom
  explicit Parser(const std::string& s) : t(s) {}
  [[noreturn]] void fail(const std::string& why) const {
    throw FFError("SMARTS '" + t + "': " + why + " at position " + std::to_string(i + 1));
  }
  char peek(size_t k = 0) const { return i + k < t.size() ? t[i + k] : '\0'; }
  int number(int def) {
    if (!std::isdigit(static_cast<unsigned char>(peek()))) return def;
    int v = 0;
    while (std::isdigit(static_cast<unsigned char>(peek()))) v = v * 10 + (t[i++] - '0');
    return v;
  }
  static std::unique_ptr<AtomExpr> prim(Prim p, int v = 0) {
    auto e = std::make_unique<AtomExpr>();
    e->prim = p;
    e->value = v;
    return e;
  }
  static std::unique_ptr<AtomExpr> elem(int z, int arom) {
    auto e = prim(Prim::Element, z);
    e->arom = arom;
    return e;
  }
  // element symbol at the cursor (two letters first when they name an element), aromatic when lower case
  bool element_symbol(std::unique_ptr<AtomExpr>& out, bool bracket) {
    const char c = peek();
    if (bracket && std::isupper(static_cast<unsigned char>(c))) {
      if (std::islower(static_cast<unsigned char>(peek(1)))) {
        const std::string two = t.substr(i, 2);
        if (int z = element_from_symbol(two); z > 0) { i += 2; out = elem(z, 0); return true; }
      }
      if (int z = element_from_symbol(std::string(1, c)); z > 0) { ++i; out = elem(z, 0); return true; }
      return false;
    }
    if (bracket && std::islower(static_cast<unsigned char>(c))) {
      for (const char* a : {"se", "as"})
        if (t.compare(i, 2, a) == 0) { i += 2; out = elem(element_from_symbol(a), 1); return true; }
      for (const char* a : {"c", "n", "o", "s", "p", "b"})
        if (c == a[0]) { ++i; out = elem(element_from_symbol(a), 1); return true; }
      return false;
    }
    // organic subset
    for (const char* a : {"Cl", "Br"})
      if (t.compare(i, 2, a) == 0) { i += 2; out = elem(element_from_symbol(a), 0); return true; }
    for (const char* a : {"B", "C", "N", "O", "P", "S", "F", "I"})
      if (c == a[0]) { ++i; out = elem(element_from_symbol(a), 0); return true; }
    for (const char* a : {"c", "n", "o", "s", "p", "b"})
      if (c == a[0]) { ++i; out = elem(element_from_symbol(a), 1); return true; }
    return false;
  }
  std::unique_ptr<AtomExpr> primitive(bool first) {
    const char c = peek();
    if (c == '$') {
      ++i;
      if (peek() != '(') fail("'$' needs '('");
      int depth = 0;
      const size_t start = i + 1;
      for (; i < t.size(); ++i) {
        if (t[i] == '(') ++depth;
        else if (t[i] == ')' && --depth == 0) break;
      }
      if (depth) fail("unbalanced $(");
      auto e = prim(Prim::Recursive);
      const std::string sub = t.substr(start, i - start);
      ++i;
      e->rec = std::make_shared<SmartsNode>();
      compile(sub, *e->rec);
      return e;
    }
    if (c == '%') {
      ++i;
      const size_t start = i;
      while (i < t.size() && std::string(";,&]!)").find(t[i]) == std::string::npos) ++i;
      if (i == start) fail("'%' needs a type name");
      auto e = prim(Prim::Type);
      e->name = t.substr(start, i - start);
      return e;
    }
    if (c == '*') { ++i; return prim(Prim::True); }
    if (c == '{') {   // CAPS extension: {AR1} .. {AR5}, antechamber's ring classes
      const size_t close = t.find('}', i);
      if (close == std::string::npos || t.compare(i + 1, 2, "AR") != 0 || close != i + 4 || !std::isdigit(static_cast<unsigned char>(t[i + 3])))
        fail("expected {AR1} .. {AR5}");
      const int v = t[i + 3] - '0';
      i = close + 1;
      return prim(Prim::ArClass, v);
    }
    if (c == '#') { ++i; const int z = number(-1); if (z < 0) fail("'#' needs a number"); return elem(z, -1); }
    if (c == '+' || c == '-') {
      const int sign = c == '+' ? 1 : -1;
      int v = 0;
      while (peek() == c) { ++i; ++v; }
      if (v == 1 && std::isdigit(static_cast<unsigned char>(peek()))) v = number(1);
      return prim(Prim::Charge, sign * v);
    }
    if (c == '@') { while (peek() == '@') ++i; return prim(Prim::Chiral); }
    // two-letter element symbols first ([Hg], [Xe], [Ru]); then H is a hydrogen atom only as the first primitive
    // ([H], [H+]) and a hydrogen count after it ([CH3]); D X R are always counts
    if (std::isupper(static_cast<unsigned char>(c)) && std::islower(static_cast<unsigned char>(peek(1))))
      if (int z = element_from_symbol(t.substr(i, 2)); z > 0) { i += 2; return elem(z, 0); }
    if (c == 'H' && first && i == bracket_start && !std::isdigit(static_cast<unsigned char>(peek(1)))) {
      ++i;
      return elem(1, 0);
    }
    std::unique_ptr<AtomExpr> e;
    if (c != 'H' && c != 'D' && c != 'X' && c != 'R' && element_symbol(e, true)) return e;
    switch (c) {
      case 'a': ++i; return prim(Prim::Aromatic);
      case 'A': ++i; return prim(Prim::Aliphatic);
      case 'D': ++i; return prim(Prim::Degree, number(1));
      case 'X': ++i; return prim(Prim::Connect, number(1));
      case 'H': ++i; return prim(Prim::HCount, number(1));
      case 'h': ++i; return prim(Prim::ImplicitH, number(-1));
      case 'v': ++i; return prim(Prim::Valence, number(1));
      case 'R': ++i; return prim(Prim::RingCount, number(-1));
      case 'r': ++i; return prim(Prim::RingSize, number(-1));
      case 'x': ++i; return prim(Prim::RingBonds, number(-1));
      default: break;
    }
    fail(std::string("unknown atom primitive '") + c + "'");
  }
  std::unique_ptr<AtomExpr> a_not(bool first) {
    if (peek() == '!') {
      ++i;
      auto e = std::make_unique<AtomExpr>();
      e->op = Op::Not;
      e->kids.push_back(a_not(false));
      return e;
    }
    return primitive(first);
  }
  std::unique_ptr<AtomExpr> a_and(bool first) {
    auto e = a_not(first);
    while (true) {
      if (peek() == '&') ++i;
      else if (peek() == ']' || peek() == ',' || peek() == ';' || peek() == '\0') break;
      auto x = std::make_unique<AtomExpr>();
      x->op = Op::And;
      x->kids.push_back(std::move(e));
      x->kids.push_back(a_not(false));
      e = std::move(x);
    }
    return e;
  }
  std::unique_ptr<AtomExpr> a_or(bool first) {
    auto e = a_and(first);
    while (peek() == ',') {
      ++i;
      auto x = std::make_unique<AtomExpr>();
      x->op = Op::Or;
      x->kids.push_back(std::move(e));
      x->kids.push_back(a_and(true));
      e = std::move(x);
    }
    return e;
  }
  std::unique_ptr<AtomExpr> a_low() {
    auto e = a_or(true);
    while (peek() == ';') {
      ++i;
      auto x = std::make_unique<AtomExpr>();
      x->op = Op::AndLow;
      x->kids.push_back(std::move(e));
      x->kids.push_back(a_or(true));
      e = std::move(x);
    }
    return e;
  }
  bool bond_char(char c) const { return c == '-' || c == '=' || c == '#' || c == ':' || c == '~' || c == '@' || c == '!' ||
                                        c == '/' || c == '\\' || c == ',' || c == ';' || c == '&'; }
  std::unique_ptr<BondExpr> bond() {
    // precedence as for atoms: ! > & (implicit) > , > ;
    std::function<std::unique_ptr<BondExpr>()> b_not, b_and, b_or, b_low;
    b_not = [&]() -> std::unique_ptr<BondExpr> {
      if (peek() == '!') {
        ++i;
        auto e = std::make_unique<BondExpr>();
        e->op = Op::Not;
        e->kids.push_back(b_not());
        return e;
      }
      auto e = std::make_unique<BondExpr>();
      switch (peek()) {
        case '-': e->prim = BPrim::Single; break;
        case '=': e->prim = BPrim::Double; break;
        case '#': e->prim = BPrim::Triple; break;
        case ':': e->prim = BPrim::Aromatic; break;
        case '~': e->prim = BPrim::Any; break;
        case '@': e->prim = BPrim::Ring; break;
        case '/': case '\\': e->prim = BPrim::Single; break;
        default: fail("bad bond");
      }
      ++i;
      return e;
    };
    b_and = [&]() {
      auto e = b_not();
      while (true) {
        if (peek() == '&') ++i;
        else if (!(peek() == '-' || peek() == '=' || peek() == '#' || peek() == ':' || peek() == '~' || peek() == '@' || peek() == '!')) break;
        auto x = std::make_unique<BondExpr>();
        x->op = Op::And;
        x->kids.push_back(std::move(e));
        x->kids.push_back(b_not());
        e = std::move(x);
      }
      return e;
    };
    b_or = [&]() {
      auto e = b_and();
      while (peek() == ',') {
        ++i;
        auto x = std::make_unique<BondExpr>();
        x->op = Op::Or;
        x->kids.push_back(std::move(e));
        x->kids.push_back(b_and());
        e = std::move(x);
      }
      return e;
    };
    b_low = [&]() {
      auto e = b_or();
      while (peek() == ';') {
        ++i;
        auto x = std::make_unique<BondExpr>();
        x->op = Op::AndLow;
        x->kids.push_back(std::move(e));
        x->kids.push_back(b_or());
        e = std::move(x);
      }
      return e;
    };
    if (!bond_char(peek())) return std::make_unique<BondExpr>();   // default: single or aromatic
    return b_low();
  }
  static void add_bond(SmartsNode& g, int a, int b, std::unique_ptr<BondExpr> e) {
    g.bonds.push_back({a, b, std::move(e)});
    const int k = int(g.bonds.size()) - 1;
    g.adj[a].push_back({b, k});
    g.adj[b].push_back({a, k});
  }
  static bool uses_types(const AtomExpr& e) {
    if (e.prim == Prim::Type && e.op == Op::Prim) return true;
    if (e.rec && e.rec->uses_types) return true;
    for (const auto& k : e.kids)
      if (uses_types(*k)) return true;
    return false;
  }
  static void compile(const std::string& text, SmartsNode& g) {
    Parser P(text);
    std::vector<int> stack;
    int prev = -1;
    std::unique_ptr<BondExpr> pending;
    std::map<int, std::pair<int, std::unique_ptr<BondExpr>>> ring;
    while (P.i < text.size()) {
      const char c = P.peek();
      if (c == '(') {
        if (prev < 0) P.fail("branch before an atom");
        stack.push_back(prev);
        ++P.i;
        continue;
      }
      if (c == ')') {
        if (stack.empty()) P.fail("unbalanced ')'");
        prev = stack.back();
        stack.pop_back();
        ++P.i;
        continue;
      }
      if (c == '.') P.fail("disconnected patterns are not supported");
      if (P.bond_char(c) && c != ',' && c != ';' && c != '&') {
        if (prev < 0) P.fail("bond before an atom");
        pending = P.bond();
        continue;
      }
      if (std::isdigit(static_cast<unsigned char>(c)) || c == '%') {
        if (prev < 0) P.fail("ring closure before an atom");
        int d;
        if (c == '%') {
          ++P.i;
          d = P.number(-1);
          if (d < 0) P.fail("'%' ring closure needs digits");
        } else {
          d = c - '0';
          ++P.i;
        }
        auto it = ring.find(d);
        if (it == ring.end()) {
          ring[d] = {prev, std::move(pending)};
        } else {
          auto e = pending ? std::move(pending) : std::move(it->second.second);
          if (!e) e = std::make_unique<BondExpr>();
          add_bond(g, it->second.first, prev, std::move(e));
          ring.erase(it);
        }
        pending.reset();
        continue;
      }
      std::unique_ptr<AtomExpr> a;
      if (c == '[') {
        ++P.i;
        P.bracket_start = P.i;
        a = P.a_low();
        if (P.peek() != ']') P.fail("expected ']'");
        ++P.i;
      } else if (c == '*') {
        ++P.i;
        a = prim(Prim::True);
      } else if (c == 'a' || c == 'A') {
        ++P.i;
        a = prim(c == 'a' ? Prim::Aromatic : Prim::Aliphatic);
      } else if (!P.element_symbol(a, false)) {
        P.fail(std::string("unexpected '") + c + "'");
      }
      g.atoms.push_back(std::move(a));
      g.adj.emplace_back();
      const int cur = int(g.atoms.size()) - 1;
      if (prev >= 0) add_bond(g, prev, cur, pending ? std::move(pending) : std::make_unique<BondExpr>());
      pending.reset();
      prev = cur;
    }
    if (!ring.empty()) P.fail("unclosed ring bond");
    if (!stack.empty()) P.fail("unbalanced '('");
    if (g.atoms.empty()) P.fail("empty pattern");
    for (const auto& a : g.atoms) g.uses_types = g.uses_types || uses_types(*a);
  }
};

struct Ctx {
  const System& s;
  const Perception& p;
  const std::vector<std::string>& types;
};

bool match_node(const SmartsNode& g, const Ctx& c, uint32_t anchor);

bool eval_atom(const AtomExpr& e, const Ctx& c, uint32_t a) {
  switch (e.op) {
    case Op::Not: return !eval_atom(*e.kids[0], c, a);
    case Op::And: case Op::AndLow: return eval_atom(*e.kids[0], c, a) && eval_atom(*e.kids[1], c, a);
    case Op::Or: return eval_atom(*e.kids[0], c, a) || eval_atom(*e.kids[1], c, a);
    case Op::Prim: break;
  }
  const auto& p = c.p;
  switch (e.prim) {
    case Prim::True: case Prim::Chiral: return true;
    case Prim::Aromatic: return p.aromatic[a];
    case Prim::Aliphatic: return !p.aromatic[a];
    case Prim::Element:
      if (c.s.atoms[a].element != e.value) return false;
      return e.arom < 0 || (e.arom == 1) == bool(p.aromatic[a]);
    case Prim::Degree: return int(p.nb[a].size()) == e.value;
    case Prim::Connect: return int(p.nb[a].size()) + (a < p.implicit_h.size() ? p.implicit_h[a] : 0) == e.value;   // implicit: united atoms
    case Prim::HCount: return p.hcount[a] == e.value;
    case Prim::ImplicitH: return (a < p.implicit_h.size() ? p.implicit_h[a] : 0) == e.value;
    case Prim::Valence: {
      int v = a < p.implicit_h.size() ? p.implicit_h[a] : 0;
      for (int o : p.order[a]) v += o;
      return v == e.value;
    }
    case Prim::RingCount: return e.value < 0 ? p.ring_count[a] > 0 : p.ring_count[a] == e.value;
    case Prim::RingSize: {
      if (e.value < 0) return p.smallest_ring[a] > 0;
      if (e.value == 0) return p.smallest_ring[a] == 0;
      // r<n>: in a ring of that size (Daylight: smallest; here any SSSR ring, so fused atoms match both sizes)
      for (const auto& r : p.rings)
        if (int(r.size()) == e.value && std::find(r.begin(), r.end(), a) != r.end()) return true;
      return false;
    }
    case Prim::RingBonds: return e.value < 0 ? p.ring_bonds[a] > 0 : p.ring_bonds[a] == e.value;
    case Prim::Charge: return p.charge[a] == e.value;
    case Prim::Recursive: return match_node(*e.rec, c, a);
    case Prim::Type: return a < c.types.size() && c.types[a] == e.name;
    case Prim::ArClass: return p.ar_class[a] == e.value;
  }
  return false;
}

bool eval_bond(const BondExpr& e, const Ctx& c, uint32_t a, size_t k) {
  switch (e.op) {
    case Op::Not: return !eval_bond(*e.kids[0], c, a, k);
    case Op::And: case Op::AndLow: return eval_bond(*e.kids[0], c, a, k) && eval_bond(*e.kids[1], c, a, k);
    case Op::Or: return eval_bond(*e.kids[0], c, a, k) || eval_bond(*e.kids[1], c, a, k);
    case Op::Prim: break;
  }
  const bool ar = c.p.arom_bond[a][k];
  const int o = c.p.order[a][k];
  switch (e.prim) {
    case BPrim::Default: return ar || o == 1;
    case BPrim::Single: case BPrim::Up: return !ar && o == 1;
    case BPrim::Double: return !ar && o == 2;
    case BPrim::Triple: return o == 3;
    case BPrim::Aromatic: return ar;
    case BPrim::Any: return true;
    case BPrim::Ring: return c.p.ring_bond[a][k];
  }
  return false;
}

bool match_node(const SmartsNode& g, const Ctx& c, uint32_t anchor) {
  const size_t m = g.atoms.size();
  if (!eval_atom(*g.atoms[0], c, anchor)) return false;
  std::vector<int64_t> map(m, -1);
  std::vector<char> used(c.s.atoms.size(), 0);
  map[0] = anchor;
  used[anchor] = 1;
  // pattern atoms in parse order: each atom k > 0 is bonded to an earlier one (its first bond)
  std::function<bool(size_t)> go = [&](size_t k) -> bool {
    if (k == m) return true;
    int parent = -1, pb = -1;
    for (auto [o, b] : g.adj[k])
      if (size_t(o) < k) { parent = o; pb = b; break; }
    if (parent < 0) return false;
    const uint32_t pa = uint32_t(map[parent]);
    for (size_t q = 0; q < c.p.nb[pa].size(); ++q) {
      const uint32_t t = c.p.nb[pa][q];
      if (used[t]) continue;
      if (!eval_bond(*g.bonds[pb].expr, c, pa, q)) continue;
      if (!eval_atom(*g.atoms[k], c, t)) continue;
      // bonds to other already-mapped atoms (ring closures)
      bool ok = true;
      for (auto [o, b] : g.adj[k]) {
        if (b == pb || size_t(o) >= k) continue;
        const uint32_t ta = uint32_t(map[o]);
        size_t kk = SIZE_MAX;
        for (size_t z = 0; z < c.p.nb[t].size(); ++z)
          if (c.p.nb[t][z] == ta) kk = z;
        if (kk == SIZE_MAX || !eval_bond(*g.bonds[b].expr, c, t, kk)) { ok = false; break; }
      }
      if (!ok) continue;
      map[k] = t;
      used[t] = 1;
      if (go(k + 1)) return true;
      used[t] = 0;
      map[k] = -1;
    }
    return false;
  };
  return go(1);
}

// Every complete embedding with its first atom on `anchor`, avoiding atoms marked in `taken` (may be null).
void match_all(const SmartsNode& g, const Ctx& c, uint32_t anchor, const std::vector<char>* taken, std::vector<std::vector<uint32_t>>& out) {
  const size_t m = g.atoms.size();
  if ((taken && (*taken)[anchor]) || !eval_atom(*g.atoms[0], c, anchor)) return;
  std::vector<int64_t> map(m, -1);
  std::vector<char> used(c.s.atoms.size(), 0);
  if (taken) used = *taken;
  map[0] = anchor;
  used[anchor] = 1;
  std::function<void(size_t)> go = [&](size_t k) {
    if (k == m) {
      out.emplace_back(map.begin(), map.end());
      return;
    }
    int parent = -1, pb = -1;
    for (auto [o, b] : g.adj[k])
      if (size_t(o) < k) { parent = o; pb = b; break; }
    if (parent < 0) return;
    const uint32_t pa = uint32_t(map[parent]);
    for (size_t q = 0; q < c.p.nb[pa].size(); ++q) {
      const uint32_t t = c.p.nb[pa][q];
      if (used[t] || !eval_bond(*g.bonds[pb].expr, c, pa, q) || !eval_atom(*g.atoms[k], c, t)) continue;
      bool ok = true;
      for (auto [o, b] : g.adj[k]) {
        if (b == pb || size_t(o) >= k) continue;
        const uint32_t ta = uint32_t(map[o]);
        size_t kk = SIZE_MAX;
        for (size_t z = 0; z < c.p.nb[t].size(); ++z)
          if (c.p.nb[t][z] == ta) kk = z;
        if (kk == SIZE_MAX || !eval_bond(*g.bonds[b].expr, c, t, kk)) { ok = false; break; }
      }
      if (!ok) continue;
      map[k] = t;
      used[t] = 1;
      go(k + 1);
      used[t] = 0;
      map[k] = -1;
    }
  };
  go(1);
}

}  // namespace

struct Smarts::Impl {
  std::string text;
  SmartsNode g;
};

Smarts::Smarts(const std::string& text) : d_(new Impl) {
  d_->text = text;
  try {
    Parser::compile(text, d_->g);
  } catch (...) {
    delete d_;
    throw;
  }
}
Smarts::~Smarts() { delete d_; }
Smarts::Smarts(Smarts&& o) noexcept : d_(o.d_) { o.d_ = nullptr; }
Smarts& Smarts::operator=(Smarts&& o) noexcept {
  std::swap(d_, o.d_);
  return *this;
}
bool Smarts::matches(const System& s, const Perception& p, uint32_t atom, const std::vector<std::string>& types) const {
  Ctx c{s, p, types};
  return match_node(d_->g, c, atom);
}
std::vector<std::vector<uint32_t>> Smarts::embeddings(const System& s, const Perception& p, uint32_t atom, const std::vector<char>* taken) const {
  static const std::vector<std::string> none;
  Ctx c{s, p, none};
  std::vector<std::vector<uint32_t>> out;
  match_all(d_->g, c, atom, taken, out);
  return out;
}
size_t Smarts::size() const { return d_->g.atoms.size(); }
bool Smarts::uses_types() const { return d_->g.uses_types; }
const std::string& Smarts::text() const { return d_->text; }

// ---------------------------------------------------------------------------------------------------------------------
// Typing

TypingResult assign_types(const System& s, const FFDef& ff, const std::vector<std::string>* context) {
  struct Rule { const TypingRule* r; const FFType* t; size_t order; Smarts sm; };
  std::map<std::string, const FFType*> byname;
  for (const auto& t : ff.types) byname[t.name] = &t;
  std::vector<Rule> rules;
  std::deque<FFType> absent;   // types the rules name but this file lacks (typing_unknown_untyped): names only
  for (size_t k = 0; k < ff.typing.size(); ++k) {
    const auto& tr = ff.typing[k];
    auto it = byname.find(tr.type);
    if (it == byname.end() && ff.typing_unknown_untyped) {
      absent.push_back(FFType{});
      absent.back().name = tr.type;
      it = byname.emplace(tr.type, &absent.back()).first;
    }
    if (it == byname.end()) throw FFError("typing rule for " + tr.type + ": no such type in " + ff.name);
    try {
      rules.push_back({&tr, it->second, k, Smarts(tr.smarts)});
    } catch (const FFError& e) {
      throw FFError("typing rule for " + tr.type + ": " + e.what());
    }
  }
  if (rules.empty()) throw FFError(ff.name + " has no typing rules");
  // the first atom of each rule decides the element cheaply
  const Perception p = perceive(s);
  const size_t n = s.atoms.size();
  TypingResult r;
  r.notes = p.notes;
  r.types.assign(n, "");
  r.why.assign(n, "");
  r.candidates.assign(n, {});
  bool any_type_refs = false;
  for (const auto& x : rules) any_type_refs = any_type_refs || x.sm.uses_types();
  // conditions on the whole structure: which elements it holds
  std::set<int> present;
  for (const auto& at : s.atoms) present.insert(at.element);
  auto applies = [&](const Rule& x, uint32_t a) {
    for (int z : x.r->needs_elements)
      if (!present.count(z)) return false;
    for (int z : x.r->no_elements)
      if (present.count(z)) return false;
    return x.r->atom_name.empty() || s.atoms[a].name == x.r->atom_name;
  };
  for (int pass = 0; pass < (any_type_refs && !context ? 10 : 1); ++pass) {
    std::vector<std::string> next(n);
    std::vector<std::string> why(n);
    std::vector<int> rule_of(n, -1);
    std::vector<std::vector<std::string>> cands(n);
    int untyped = 0, ambiguous = 0;
    for (uint32_t a = 0; a < n; ++a) {
      std::vector<const Rule*> hit;
      for (const auto& x : rules)   // the SMARTS decides the element (library element columns have typos)
        if (applies(x, a) && x.sm.matches(s, p, a, context ? *context : r.types)) hit.push_back(&x);
      for (const auto* h : hit)
        if (std::find(cands[a].begin(), cands[a].end(), h->t->name) == cands[a].end()) cands[a].push_back(h->t->name);
      if (hit.empty()) { ++untyped; why[a] = "no rule matches"; continue; }
      // drop types that another matching type overrides
      std::set<std::string> over;
      for (const auto* h : hit)
        for (const auto& o : h->r->overrides) over.insert(o);
      std::vector<const Rule*> keep;
      for (const auto* h : hit)
        if (!over.count(h->t->name)) keep.push_back(h);
      if (keep.empty()) keep = hit;   // a cycle of overrides: fall back to priority
      std::stable_sort(keep.begin(), keep.end(), [](const Rule* x, const Rule* y) {
        if (x->r->priority != y->r->priority) return x->r->priority > y->r->priority;
        return x->order < y->order;
      });
      next[a] = keep[0]->t->name;
      rule_of[a] = int(keep[0]->order);
      why[a] = keep[0]->t->name + " ← " + keep[0]->sm.text();
      if (!ff.typing_ordered && keep.size() > 1 && keep[1]->r->priority == keep[0]->r->priority && keep[1]->t != keep[0]->t) {
        ++ambiguous;
        why[a] += "  (also matched";
        for (size_t q = 1; q < keep.size() && keep[q]->r->priority == keep[0]->r->priority; ++q)
          if (keep[q]->t != keep[0]->t) why[a] += " " + keep[q]->t->name;
        why[a] += "; first rule taken)";
      }
    }
    const bool same = next == r.types;
    r.types = std::move(next);
    r.why = std::move(why);
    r.rule = std::move(rule_of);
    r.candidates = std::move(cands);
    r.untyped = untyped;
    r.ambiguous = ambiguous;
    if (same) break;
  }
  // conjugated pairs: colour each conjugated system along its bonds (GAFF: double bond → the other member, single →
  // the same; CGenFF pair_mode double_same: the reverse)
  if (!ff.typing_pairs.empty()) {
    std::map<std::string, std::pair<size_t, int>> member;   // type → (pair, 0 first / 1 second)
    for (size_t k = 0; k < ff.typing_pairs.size(); ++k) {
      member.emplace(ff.typing_pairs[k].first, std::make_pair(k, 0));
      member.emplace(ff.typing_pairs[k].second, std::make_pair(k, 1));
    }
    std::vector<int> colour(n, -1);
    int conflicts = 0;
    for (uint32_t s0 = 0; s0 < n; ++s0) {
      if (colour[s0] >= 0 || !member.count(r.types[s0])) continue;
      colour[s0] = 0;
      std::vector<uint32_t> queue = {s0};
      for (size_t h = 0; h < queue.size(); ++h) {
        const uint32_t a = queue[h];
        for (size_t q = 0; q < p.nb[a].size(); ++q) {
          const uint32_t b = p.nb[a][q];
          const int o = p.order[a][q];
          if (!member.count(r.types[b]) || (o != 1 && o != 2)) continue;
          const bool differ = (o == 2) != ff.typing_pairs_double_same;
          const int want = differ ? 1 - colour[a] : colour[a];
          if (colour[b] < 0) {
            colour[b] = want;
            queue.push_back(b);
          } else if (colour[b] != want) {
            ++conflicts;
          }
        }
      }
    }
    for (uint32_t a = 0; a < n; ++a) {
      if (colour[a] < 0) continue;
      const auto [k, which] = member.at(r.types[a]);
      const std::string t = colour[a] == 0 ? ff.typing_pairs[k].first : ff.typing_pairs[k].second;
      if (t != r.types[a]) {
        r.why[a] += "  (" + t + " for the conjugation pattern)";
        r.types[a] = t;
      }
    }
    if (conflicts) r.notes.push_back(std::to_string(conflicts / 2) + " conjugated bonds could not follow the cc/cd pattern (odd ring)");
  }
  refine_bond_order_variants(s, p, ff, r.types, r.why);
  // DREIDING hydrogen bonds: an N, O or F carrying a hydrogen takes its donor variant (T_hd) and that hydrogen the
  // hydrogen-bond type; one without takes its acceptor variant (T_ha), where the force field has them
  if (!ff.hbonds.terms.empty()) {
    std::set<std::string> names;
    for (const auto& t : ff.types) names.insert(t.name);
    const std::string hb = ff.hbonds.terms.front().hydrogen;
    int nd = 0, na = 0;
    for (size_t a = 0; a < s.atoms.size(); ++a) {
      const int z = s.atoms[a].element;
      if ((z != 7 && z != 8 && z != 9) || r.types[a].empty()) continue;
      std::vector<uint32_t> hs;
      for (uint32_t b : p.nb[a])
        if (s.atoms[b].element == 1 && r.types[b] == "H") hs.push_back(b);
      if (!hs.empty() && names.count(r.types[a] + "_hd") && names.count(hb)) {
        r.types[a] += "_hd";
        r.why[a] += "  (hydrogen-bond donor)";
        for (uint32_t b : hs) {
          r.types[b] = hb;
          r.why[b] += "  (on a hydrogen-bond donor)";
        }
        ++nd;
      } else if (names.count(r.types[a] + "_ha")) {
        r.types[a] += "_ha";
        r.why[a] += "  (hydrogen-bond acceptor)";
        ++na;
      }
    }
    if (nd + na) r.notes.push_back(std::to_string(nd) + " hydrogen-bond donors, " + std::to_string(na) + " acceptors");
  }
  return r;
}

}  // namespace caps
