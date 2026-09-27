#include "caps/example_typing.hpp"

#include <algorithm>
#include <stdexcept>

namespace caps {

namespace {
uint64_t mix(uint64_t h, uint64_t v) {   // splitmix-style combine
  h ^= v + 0x9E3779B97F4A7C15ull + (h << 6) + (h >> 2);
  h ^= h >> 31;
  h *= 0xBF58476D1CE4E5B9ull;
  return h ^ (h >> 29);
}
int order_of(int o) { return o == 5 ? 1 : o <= 0 ? 1 : o; }
}  // namespace

std::vector<std::vector<uint64_t>> environment_signatures(const System& s, int radius) {
  const size_t n = s.atoms.size();
  std::vector<std::vector<std::pair<uint32_t, int>>> nb(n);
  for (const auto& b : s.bonds) {
    nb[b.i].push_back({b.j, order_of(b.order)});
    nb[b.j].push_back({b.i, order_of(b.order)});
  }
  // the smallest ring through each atom (0: none; up to 8), from the bonds alone
  std::vector<int> ring(n, 0);
  for (size_t a = 0; a < n; ++a) {
    // shortest cycle through a: BFS from a, a cycle closes where two branches from different first neighbours meet
    std::vector<int> dist(n, -1), root(n, -1);
    std::vector<uint32_t> q{uint32_t(a)};
    dist[a] = 0;
    int best = 0;
    for (size_t h = 0; h < q.size() && !best; ++h) {
      const uint32_t u = q[h];
      if (dist[u] >= 4) break;
      for (auto [w, o] : nb[u]) {
        if (dist[w] < 0) {
          dist[w] = dist[u] + 1, root[w] = u == a ? int(w) : root[u];
          q.push_back(w);
        } else if (w != a && root[w] != root[u] && u != a) {
          const int len = dist[u] + dist[w] + 1;
          if (len <= 8 && (!best || len < best)) best = len;
        }
      }
    }
    ring[a] = best;
  }
  std::vector<std::vector<uint64_t>> out;
  std::vector<uint64_t> cur(n);
  for (size_t i = 0; i < n; ++i) {
    int heavy = 0, h = 0;
    for (auto [j, o] : nb[i]) (s.atoms[j].element == 1 ? h : heavy)++;
    uint64_t v = mix(uint64_t(s.atoms[i].element), uint64_t(heavy));
    v = mix(v, uint64_t(h));
    v = mix(v, uint64_t(ring[i]));
    cur[i] = v;
  }
  out.push_back(cur);
  for (int r = 1; r <= radius; ++r) {
    std::vector<uint64_t> next(n);
    for (size_t i = 0; i < n; ++i) {
      std::vector<uint64_t> parts;
      for (auto [j, o] : nb[i]) parts.push_back(cur[j]);   // bonds only: files differ in whether they give orders
      std::sort(parts.begin(), parts.end());
      uint64_t v = mix(cur[i], 0x51ED27ull + uint64_t(r));
      for (uint64_t p : parts) v = mix(v, p);
      next[i] = v;
    }
    cur = next;
    out.push_back(cur);
  }
  return out;
}

ExampleTypes learn_types(const System& example, const std::vector<std::string>& types, int min_radius, int max_radius) {
  if (types.size() != example.atoms.size()) throw std::invalid_argument("one type per example atom");
  if (std::all_of(types.begin(), types.end(), [](const std::string& t) { return t.empty(); })) throw std::invalid_argument("no atom of the example has a type");
  const auto sig = environment_signatures(example, max_radius);
  ExampleTypes t;
  t.by_radius.resize(size_t(max_radius) + 1);
  std::vector<bool> ambiguous_at(size_t(max_radius) + 1, false);
  for (int r = 0; r <= max_radius; ++r) {
    std::map<uint64_t, std::string> m;
    std::map<uint64_t, bool> bad;
    for (size_t i = 0; i < types.size(); ++i) {
      if (types[i].empty()) continue;
      auto [it, fresh] = m.emplace(sig[size_t(r)][i], types[i]);
      if (!fresh && it->second != types[i]) bad[sig[size_t(r)][i]] = true;
    }
    for (const auto& [k, v] : bad) m.erase(k), ambiguous_at[size_t(r)] = true;
    t.by_radius[size_t(r)] = std::move(m);
  }
  t.radius = max_radius;
  for (int r = std::max(0, min_radius); r <= max_radius; ++r)
    if (!ambiguous_at[size_t(r)]) { t.radius = r; break; }
  if (ambiguous_at[size_t(t.radius)]) {   // two types for one environment even at the largest radius
    std::map<uint64_t, std::vector<size_t>> by;
    for (size_t i = 0; i < types.size(); ++i) if (!types[i].empty()) by[sig[size_t(max_radius)][i]].push_back(i);
    for (const auto& [k, idx] : by) {
      std::string first = types[idx.front()];
      for (size_t i : idx)
        if (types[i] != first) {
          t.conflicts.push_back("atoms " + std::to_string(idx.front() + 1) + " and " + std::to_string(i + 1) + " have the same surroundings but types " + first +
                                " and " + types[i]);
          break;
        }
    }
  }
  t.environments = t.by_radius[size_t(t.radius)].size();
  return t;
}

ExampleMatch apply_types(const System& s, const ExampleTypes& t) {
  const int rmax = int(t.by_radius.size()) - 1;
  const auto sig = environment_signatures(s, rmax);
  ExampleMatch m;
  m.types.assign(s.atoms.size(), "");
  for (size_t i = 0; i < s.atoms.size(); ++i) {
    // the exact radius first, then larger ones (more specific, still consistent), then shorter ones
    bool done = false;
    for (int r = t.radius; r <= rmax && !done; ++r)
      if (auto it = t.by_radius[size_t(r)].find(sig[size_t(r)][i]); it != t.by_radius[size_t(r)].end()) m.types[i] = it->second, done = true, ++m.exact;
    for (int r = t.radius - 1; r >= 1 && !done; --r)
      if (auto it = t.by_radius[size_t(r)].find(sig[size_t(r)][i]); it != t.by_radius[size_t(r)].end()) m.types[i] = it->second, done = true, ++m.shorter;
    if (!done) ++m.unmatched;
  }
  return m;
}

std::vector<uint32_t> equivalent_atoms(const System& s, uint32_t atom, int radius) {
  const auto sig = environment_signatures(s, radius);
  std::vector<uint32_t> out;
  if (atom >= s.atoms.size()) return out;
  for (size_t i = 0; i < s.atoms.size(); ++i)
    if (sig[size_t(radius)][i] == sig[size_t(radius)][atom]) out.push_back(uint32_t(i));
  return out;
}

}  // namespace caps
