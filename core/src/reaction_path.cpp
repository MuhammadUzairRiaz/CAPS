// CAPS reaction paths: atom matching, path interpolation, climbing-image NEB (see caps/reaction_path.hpp).
#include "caps/reaction_path.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <map>
#include <stdexcept>
#include <tuple>

#include "caps/analysis.hpp"
#include "caps/superpose.hpp"

namespace caps {

namespace {

constexpr double kBohr = 0.529177210903;   // Å

// Hungarian method (Kuhn 1955) for a square cost matrix, by potentials (O(n³)): row r goes to column out[r].
std::vector<int> hungarian(const std::vector<std::vector<double>>& c) {
  const int n = int(c.size());
  const double inf = std::numeric_limits<double>::infinity();
  std::vector<double> u(size_t(n + 1), 0), v(size_t(n + 1), 0), minv(size_t(n + 1));
  std::vector<int> p(size_t(n + 1), 0), way(size_t(n + 1), 0);
  std::vector<char> used(size_t(n + 1));
  for (int i = 1; i <= n; ++i) {
    p[0] = i;
    int j0 = 0;
    std::fill(minv.begin(), minv.end(), inf);
    std::fill(used.begin(), used.end(), 0);
    do {
      used[size_t(j0)] = 1;
      const int i0 = p[size_t(j0)];
      double delta = inf;
      int j1 = 0;
      for (int j = 1; j <= n; ++j) {
        if (used[size_t(j)]) continue;
        const double cur = c[size_t(i0 - 1)][size_t(j - 1)] - u[size_t(i0)] - v[size_t(j)];
        if (cur < minv[size_t(j)]) minv[size_t(j)] = cur, way[size_t(j)] = j0;
        if (minv[size_t(j)] < delta) delta = minv[size_t(j)], j1 = j;
      }
      for (int j = 0; j <= n; ++j) {
        if (used[size_t(j)]) u[size_t(p[size_t(j)])] += delta, v[size_t(j)] -= delta;
        else minv[size_t(j)] -= delta;
      }
      j0 = j1;
    } while (p[size_t(j0)] != 0);
    do {
      const int j1 = way[size_t(j0)];
      p[size_t(j0)] = p[size_t(j1)];
      j0 = j1;
    } while (j0);
  }
  std::vector<int> out(size_t(n), -1);
  for (int j = 1; j <= n; ++j)
    if (p[size_t(j)] > 0) out[size_t(p[size_t(j)] - 1)] = j - 1;
  return out;
}

using Adjacency = std::vector<std::vector<uint32_t>>;

// Each atom's bonded neighbours: the structure's bonds, or bonds perceived from the distances when it has none.
Adjacency adjacency(const System& s) {
  Adjacency nb(s.atoms.size());
  const auto bonds = s.bonds.empty() ? perceive_bonds(s) : s.bonds;
  for (const auto& bd : bonds) nb[bd.i].push_back(bd.j), nb[bd.j].push_back(bd.i);
  return nb;
}

bool bonded(const Adjacency& nb, uint32_t x, uint32_t y) { return std::find(nb[x].begin(), nb[x].end(), y) != nb[x].end(); }

// Bonded-environment labels of both structures from one dictionary (element, then the sorted neighbours' labels, three
// rounds), so equal environments get equal numbers in A and in B.
void environment_labels(const System& a, const System& b, const Adjacency& na, const Adjacency& nb, std::vector<int>& la, std::vector<int>& lb) {
  std::map<std::string, int> dict;
  la.resize(a.atoms.size()), lb.resize(b.atoms.size());
  for (size_t i = 0; i < a.atoms.size(); ++i) la[i] = a.atoms[i].element;
  for (size_t i = 0; i < b.atoms.size(); ++i) lb[i] = b.atoms[i].element;
  for (int round = 0; round < 3; ++round) {
    auto step = [&](const std::vector<int>& l, const Adjacency& nbr) {
      std::vector<int> out(l.size());
      for (size_t i = 0; i < l.size(); ++i) {
        std::vector<int> ns;
        for (uint32_t j : nbr[i]) ns.push_back(l[j]);
        std::sort(ns.begin(), ns.end());
        std::string key = std::to_string(l[i]) + ":";
        for (int x : ns) key += std::to_string(x) + ",";
        out[i] = dict.emplace(key, int(dict.size()) + 1000).first->second;
      }
      return out;
    };
    la = step(la, na), lb = step(lb, nb);
  }
}

// Does the map send every atom to one of its own environment and every bond of A onto a bond of B? (With equal bond
// counts, that makes it a graph isomorphism.)
bool keeps_bonds(const Adjacency& na, const Adjacency& nb, const std::vector<int>& la, const std::vector<int>& lb, const std::vector<int>& map) {
  for (size_t i = 0; i < map.size(); ++i) {
    if (map[i] < 0 || la[i] != lb[size_t(map[i])]) return false;
    for (uint32_t j : na[i])
      if (!bonded(nb, uint32_t(map[i]), uint32_t(map[j]))) return false;
  }
  return true;
}

// The bond-preserving map nearest the preferred one: A's atoms are visited breadth-first through each molecule (from
// its atom of the rarest environment), each given a free atom of B with its environment that is bonded to the
// partners of its already placed neighbours, the preferred partner first, then the nearest. Chronological backtracking
// within `budget` tries; empty when no such map is found.
std::vector<int> graph_map(const Adjacency& na, const Adjacency& nb, const std::vector<int>& la, const std::vector<int>& lb,
                           const std::vector<int>& pref, const std::function<double(size_t, size_t)>& d2, long budget, bool& over_budget) {
  const size_t n = na.size();
  long tries = 0;
  std::map<int, std::vector<uint32_t>> class_b;
  for (size_t y = 0; y < n; ++y) class_b[lb[y]].push_back(uint32_t(y));
  auto class_size = [&](uint32_t i) {
    const auto it = class_b.find(la[i]);
    return it == class_b.end() ? size_t(0) : it->second.size();
  };
  std::vector<uint32_t> order;
  order.reserve(n);
  std::vector<int> parent(n, -1);   // the earlier visited neighbour whose partner's neighbours are the candidates
  std::vector<char> seen(n, 0);
  for (size_t s = 0; s < n; ++s) {
    if (seen[s]) continue;
    std::vector<uint32_t> comp{uint32_t(s)};
    seen[s] = 1;
    for (size_t q = 0; q < comp.size(); ++q)
      for (uint32_t j : na[comp[q]])
        if (!seen[j]) seen[j] = 1, comp.push_back(j);
    uint32_t root = comp[0];
    for (uint32_t x : comp)
      if (class_size(x) < class_size(root)) root = x;
    if (class_size(root) == 0) return {};
    const size_t start = order.size();
    order.push_back(root), seen[root] = 2;
    for (size_t q = start; q < order.size(); ++q)
      for (uint32_t j : na[order[q]])
        if (seen[j] != 2) seen[j] = 2, parent[j] = int(order[q]), order.push_back(j);
  }
  std::vector<int> map(n, -1), inv(n, -1);
  std::vector<std::vector<uint32_t>> cand(n);
  std::vector<size_t> next(n, 0);
  std::vector<char> full(n, 0);   // a molecule's first atom tries its preferred partner alone before the whole class
  auto sort_cand = [&](size_t k) {
    const uint32_t i = order[k];
    std::sort(cand[k].begin(), cand[k].end(), [&](uint32_t x, uint32_t y) {
      const bool px = int(x) == pref[i], py = int(y) == pref[i];
      if (px != py) return px;
      return d2(i, x) < d2(i, y);
    });
  };
  auto whole_class = [&](size_t k, int skip) {
    const uint32_t i = order[k];
    cand[k].clear(), next[k] = 0, full[k] = 1;
    for (uint32_t y : class_b.at(la[i]))
      if (inv[y] < 0 && int(y) != skip) cand[k].push_back(y);
    sort_cand(k);
  };
  auto fill = [&](size_t k) {
    const uint32_t i = order[k];
    cand[k].clear(), next[k] = 0;
    if (parent[i] >= 0) {
      for (uint32_t y : nb[size_t(map[size_t(parent[i])])])
        if (inv[y] < 0 && lb[y] == la[i]) cand[k].push_back(y);
      full[k] = 1;
      sort_cand(k);
    } else if (pref[i] >= 0 && inv[size_t(pref[i])] < 0 && lb[size_t(pref[i])] == la[i]) {
      cand[k].push_back(uint32_t(pref[i]));
      full[k] = 0;
    } else {
      whole_class(k, -1);
    }
  };
  auto fits = [&](uint32_t i, uint32_t y) {
    if (inv[y] >= 0 || lb[y] != la[i]) return false;
    int placed = 0, used = 0;
    for (uint32_t j : na[i])
      if (map[j] >= 0) {
        ++placed;
        if (!bonded(nb, y, uint32_t(map[j]))) return false;
      }
    for (uint32_t z : nb[y]) used += inv[z] >= 0;
    return placed == used;   // no bond of B among placed atoms that A lacks
  };
  size_t k = 0;
  if (n) fill(0);
  while (k < n) {
    const uint32_t i = order[k];
    if (map[i] >= 0) inv[size_t(map[i])] = -1, map[i] = -1;   // back here: the last choice failed further on
    bool placed = false;
    for (;;) {
      if (next[k] == cand[k].size()) {
        if (full[k]) break;
        whole_class(k, pref[i]);   // the preferred partner was tried; now the rest of the class
        continue;
      }
      const uint32_t y = cand[k][next[k]++];
      if (++tries > budget) {
        over_budget = true;
        return {};
      }
      if (fits(i, y)) {
        map[i] = int(y), inv[y] = int(i), placed = true;
        break;
      }
    }
    if (placed) {
      if (++k < n) fill(k);
    } else {
      if (k == 0) return {};
      --k;
    }
  }
  return map;
}

Vec3 rotate(const std::array<std::array<double, 3>, 3>& R, const Vec3& x) {
  return {R[0][0] * x[0] + R[0][1] * x[1] + R[0][2] * x[2], R[1][0] * x[0] + R[1][1] * x[1] + R[1][2] * x[2], R[2][0] * x[0] + R[2][1] * x[1] + R[2][2] * x[2]};
}

}  // namespace

AtomMatch match_atoms(const System& a, const System& b) {
  const size_t n = a.atoms.size();
  if (n != b.atoms.size()) throw std::invalid_argument("matching atoms: A has " + std::to_string(n) + " atoms, B " + std::to_string(b.atoms.size()));
  if (n == 0) throw std::invalid_argument("matching atoms: no atoms");
  AtomMatch R;
  // the classes: environments when both structures have the same ones, elements otherwise
  std::vector<int> la, lb;
  const Adjacency na = adjacency(a), nb = adjacency(b);
  environment_labels(a, b, na, nb, la, lb);
  auto counts = [](const std::vector<int>& l) { std::map<int, int> c; for (int x : l) ++c[x]; return c; };
  R.by_graph = counts(la) == counts(lb);
  if (!R.by_graph) {
    for (size_t i = 0; i < n; ++i) la[i] = a.atoms[i].element, lb[i] = b.atoms[i].element;
    if (counts(la) != counts(lb)) throw std::invalid_argument("matching atoms: A and B do not have the same elements");
    R.notes.push_back("the bonded environments differ (bonds made or broken): atoms matched by element and position");
  }
  std::map<int, std::pair<std::vector<size_t>, std::vector<size_t>>> cls;
  for (size_t i = 0; i < n; ++i) cls[la[i]].first.push_back(i), cls[lb[i]].second.push_back(i);
  const bool periodic = a.cell.valid() && b.cell.valid();
  bool greedy_noted = false;
  auto assign = [&](const std::function<double(size_t, size_t)>& d2) {
    std::vector<int> map(n, -1);
    for (const auto& [lab, ij] : cls) {
      const auto& [ia, ib] = ij;
      const size_t m = ia.size();
      if (m <= 1500) {
        std::vector<std::vector<double>> c(m, std::vector<double>(m));
        for (size_t x = 0; x < m; ++x)
          for (size_t y = 0; y < m; ++y) c[x][y] = d2(ia[x], ib[y]);
        const auto h = hungarian(c);
        for (size_t x = 0; x < m; ++x) map[ia[x]] = int(ib[size_t(h[x])]);
      } else {   // large classes: nearest free partner, closest pairs first
        std::vector<std::tuple<double, size_t, size_t>> pairs;
        for (size_t x = 0; x < m; ++x)
          for (size_t y = 0; y < m; ++y) {
            const double v = d2(ia[x], ib[y]);
            if (v < 9.0) pairs.emplace_back(v, x, y);
          }
        std::sort(pairs.begin(), pairs.end());
        std::vector<char> ua(m, 0), ub(m, 0);
        for (const auto& [v, x, y] : pairs)
          if (!ua[x] && !ub[y]) ua[x] = ub[y] = 1, map[ia[x]] = int(ib[y]);
        for (size_t x = 0, y = 0; x < m; ++x)
          if (!ua[x]) {
            while (ub[y]) ++y;
            ub[y] = 1, map[ia[x]] = int(ib[y]);
          }
        if (!greedy_noted) R.notes.push_back("classes above 1500 atoms matched greedily (nearest first), not by the Hungarian method");
        greedy_noted = true;
      }
    }
    return map;
  };
  // When the environments agree, the map must also keep every bond: a position-only assignment within one class can
  // give a hydrogen of one carbon the partner of a hydrogen on the next (butane turned end to end), and B in that
  // order has bonds two angstroms long. graph_map turns the assignment into the nearest map that keeps the bonds.
  const long budget = 200 * long(n) + 100000;
  bool over_budget = false, graph_failed = false;
  auto keep_bonds = [&](std::vector<int> m, const std::function<double(size_t, size_t)>& d2) {
    // an exhausted search proves no such map exists, and a spent budget is not spent again on every iteration
    if (!R.by_graph || graph_failed) return m;
    auto g = graph_map(na, nb, la, lb, m, d2, budget, over_budget);
    if (g.empty()) graph_failed = true;
    return g.empty() ? m : g;
  };
  bool best_valid = false;
  if (periodic) {
    const auto d2 = [&](size_t i, size_t j) {
      const Vec3 d = a.cell.minimum_image(b.atoms[j].pos - a.atoms[i].pos);
      return dot(d, d);
    };
    R.map = keep_bonds(assign(d2), d2);
    best_valid = R.by_graph && keeps_bonds(na, nb, la, lb, R.map);
    double s = 0;
    for (size_t i = 0; i < n; ++i) {
      const Vec3 d = a.cell.minimum_image(b.atoms[size_t(R.map[i])].pos - a.atoms[i].pos);
      s += dot(d, d);
    }
    R.rmsd = std::sqrt(s / double(n));
  } else {
    // four proper orientations of B's principal axes onto A's as starts; assignment and fit repeated until stable
    Vec3 ca{0, 0, 0}, cb{0, 0, 0};
    for (size_t i = 0; i < n; ++i) ca = ca + a.atoms[i].pos, cb = cb + b.atoms[i].pos;
    ca = ca * (1.0 / double(n)), cb = cb * (1.0 / double(n));
    auto axes = [&](const System& s, const Vec3& c, double V[3][3]) {
      double G[3][3] = {};
      for (const auto& at : s.atoms) {
        const Vec3 d = at.pos - c;
        for (int p = 0; p < 3; ++p)
          for (int q = 0; q < 3; ++q) G[p][q] += d[size_t(p)] * d[size_t(q)];
      }
      double w[3];
      symmetric_eigen3(G, w, V);
    };
    double VA[3][3], VB[3][3];
    axes(a, ca, VA), axes(b, cb, VB);
    double best = std::numeric_limits<double>::infinity();
    if (R.by_graph && n >= 3) {   // B in its file order, when that keeps the bonds, wins a tie
      std::vector<int> id(n);
      for (size_t i = 0; i < n; ++i) id[i] = int(i);
      if (keeps_bonds(na, nb, la, lb, id)) {
        std::vector<Vec3> ref(n), mov(n);
        for (size_t i = 0; i < n; ++i) ref[i] = a.atoms[i].pos, mov[i] = b.atoms[i].pos;
        best = superpose(ref, mov).rmsd, R.map = id, best_valid = true;
      }
    }
    const int signs[4][3] = {{1, 1, 1}, {1, -1, -1}, {-1, 1, -1}, {-1, -1, 1}};
    for (const auto& sg : signs) {
      std::array<std::array<double, 3>, 3> Rm{};
      for (int p = 0; p < 3; ++p)
        for (int q = 0; q < 3; ++q)
          for (int k = 0; k < 3; ++k) Rm[size_t(p)][size_t(q)] += VA[p][k] * sg[k] * VB[q][k];
      const double det = Rm[0][0] * (Rm[1][1] * Rm[2][2] - Rm[1][2] * Rm[2][1]) - Rm[0][1] * (Rm[1][0] * Rm[2][2] - Rm[1][2] * Rm[2][0]) +
                         Rm[0][2] * (Rm[1][0] * Rm[2][1] - Rm[1][1] * Rm[2][0]);
      if (det < 0)   // the eigenvector sets differ in handedness: flip the first axis
        for (int p = 0; p < 3; ++p)
          for (int q = 0; q < 3; ++q) Rm[size_t(p)][size_t(q)] -= 2 * VA[p][0] * sg[0] * VB[q][0];
      std::vector<Vec3> pb(n);
      for (size_t j = 0; j < n; ++j) pb[j] = rotate(Rm, b.atoms[j].pos - cb) + ca;
      std::vector<int> map;
      double rmsd = std::numeric_limits<double>::infinity();
      const auto d2 = [&](size_t i, size_t j) {
        const Vec3 d = pb[j] - a.atoms[i].pos;
        return dot(d, d);
      };
      for (int it = 0; it < 30; ++it) {
        auto next = keep_bonds(assign(d2), d2);
        std::vector<Vec3> ref(n), mov(n);
        for (size_t i = 0; i < n; ++i) ref[i] = a.atoms[i].pos, mov[i] = b.atoms[size_t(next[i])].pos;
        if (n >= 3) {
          const Superposition sp = superpose(ref, mov);
          for (size_t j = 0; j < n; ++j) pb[j] = sp.apply(b.atoms[j].pos);
          rmsd = sp.rmsd;
        } else {
          double s = 0;
          for (size_t i = 0; i < n; ++i) { const Vec3 d = pb[size_t(next[i])] - a.atoms[i].pos; s += dot(d, d); }
          rmsd = std::sqrt(s / double(n));
        }
        const bool same = next == map;
        map = std::move(next);
        if (same) break;
      }
      const bool valid = R.by_graph && keeps_bonds(na, nb, la, lb, map);
      if ((valid && !best_valid) || (valid == best_valid && rmsd < best - 1e-9)) best = rmsd, R.map = map, best_valid = valid;
    }
    R.rmsd = best;
  }
  // A's bonds that B, in the matched order, stretches past 1.5 times: expected where a reaction breaks bonds, a wrong
  // match otherwise (an energy on A's bonds then sees them as huge strains).
  {
    int stretched = 0;
    std::string first;
    for (size_t i = 0; i < n; ++i)
      for (uint32_t j : na[i]) {
        if (j <= i) continue;
        Vec3 da = a.atoms[j].pos - a.atoms[i].pos, db = b.atoms[size_t(R.map[j])].pos - b.atoms[size_t(R.map[i])].pos;
        if (periodic) da = a.cell.minimum_image(da), db = b.cell.minimum_image(db);
        const double ra = norm(da), rb = norm(db);
        if (rb > 1.5 * ra) {
          if (++stretched <= 3) {
            char e[96];
            std::snprintf(e, sizeof e, "%s%zu-%u %.3f → %.3f Å", first.empty() ? "" : ", ", i + 1, j + 1, ra, rb);
            first += e;
          }
        }
      }
    if (R.by_graph && !best_valid) {
      R.notes.push_back(over_budget
                            ? "no atom map that keeps every bond was found within " + std::to_string(budget) + " tries: the closest map by position is used"
                            : "no atom map keeps every bond (the bond graphs differ): the closest map by position is used");
    }
    if (stretched > 0)
      R.notes.push_back(std::to_string(stretched) + " of A's bonds are over 1.5 times longer in B as matched (" + first +
                        (stretched > 3 ? ", ..." : "") + (R.by_graph ? "): check the match, or keep B's file order" : "): bonds the reaction breaks, or atoms matched to the wrong partner"));
  }
  for (size_t i = 0; i < n; ++i) R.reordered += R.map[i] != int(i);
  char buf[160];
  std::snprintf(buf, sizeof buf, "%d of %zu atoms in another order; RMSD %.4f Å after %s", R.reordered, n, R.rmsd,
                periodic ? "minimum-image matching" : "superposition");
  R.notes.insert(R.notes.begin(), buf);
  return R;
}

System reorder_atoms(const System& b, const std::vector<int>& map, const System& a) {
  if (map.size() != b.atoms.size() || a.atoms.size() != b.atoms.size()) throw std::invalid_argument("reorder: one partner per atom");
  System out = b;
  out.bonds.clear();
  std::vector<int> inv(map.size(), -1);
  for (size_t i = 0; i < map.size(); ++i) {
    out.atoms[i] = b.atoms[size_t(map[i])];
    inv[size_t(map[i])] = int(i);
    if (a.cell.valid()) out.atoms[i].pos = a.atoms[i].pos + a.cell.minimum_image(out.atoms[i].pos - a.atoms[i].pos);
  }
  for (const auto& bd : b.bonds) out.bonds.push_back({uint32_t(inv[bd.i]), uint32_t(inv[bd.j]), bd.order});
  return out;
}

// ---------------------------------------------------------------- the band

NebResult neb(std::vector<std::vector<Vec3>> path, const ImageEnergy& energy, const NebOptions& o) {
  const size_t m = path.size();
  if (m < 3) throw std::invalid_argument("NEB: at least one image between the ends");
  const size_t n = path[0].size();
  for (const auto& p : path)
    if (p.size() != n) throw std::invalid_argument("NEB: every image needs the same atoms");
  NebResult R;
  std::vector<double> E(m);
  std::vector<std::vector<Vec3>> F(m, std::vector<Vec3>(n)), G(m, std::vector<Vec3>(n, Vec3{0, 0, 0}));   // true and NEB forces
  E[0] = energy(0, path[0], F[0]);
  E[m - 1] = energy(int(m - 1), path[m - 1], F[m - 1]);
  auto sub = [&](const std::vector<Vec3>& x, const std::vector<Vec3>& y) {
    std::vector<Vec3> d(n);
    for (size_t k = 0; k < n; ++k) d[k] = x[k] - y[k];
    return d;
  };
  auto dotv = [&](const std::vector<Vec3>& x, const std::vector<Vec3>& y) {
    double s = 0;
    for (size_t k = 0; k < n; ++k) s += dot(x[k], y[k]);
    return s;
  };
  bool climbing = false;
  auto forces = [&]() {
    for (size_t i = 1; i + 1 < m; ++i) E[i] = energy(int(i), path[i], F[i]);
    size_t top = 1;
    for (size_t i = 1; i + 1 < m; ++i)
      if (E[i] > E[top]) top = i;
    // Climbing (Henkelman, Uberuaga, Jonsson, J. Chem. Phys. 113, 9901 (2000)) inverts the force along the tangent so the
    // image walks uphill to the saddle. That is only bounded when the image is a maximum along the band: E above both
    // neighbours (the ends count). If the unrelaxed ends lie above every interior image, the highest interior image is
    // not a maximum, the inverted force has nothing to stop it, and the image runs off. Such an image keeps the spring force.
    const bool top_is_max = E[top] > E[top - 1] && E[top] > E[top + 1];
    R.climbing = climbing && top_is_max ? int(top) : -1;
    double fmax = 0;
    for (size_t i = 1; i + 1 < m; ++i) {
      const auto tp = sub(path[i + 1], path[i]), tm = sub(path[i], path[i - 1]);
      std::vector<Vec3> tau(n);
      const double vp = E[i + 1], v0 = E[i], vm = E[i - 1];
      if (vp > v0 && v0 > vm) tau = tp;
      else if (vp < v0 && v0 < vm) tau = tm;
      else {   // at an extremum: the two neighbours' directions weighted by the energy differences
        const double dmax = std::max(std::fabs(vp - v0), std::fabs(vm - v0)), dmin = std::min(std::fabs(vp - v0), std::fabs(vm - v0));
        const double wp = vp > vm ? dmax : dmin, wm = vp > vm ? dmin : dmax;
        for (size_t k = 0; k < n; ++k) tau[k] = tp[k] * wp + tm[k] * wm;
      }
      const double tn = std::sqrt(dotv(tau, tau));
      if (tn > 0)
        for (auto& t : tau) t = t * (1 / tn);
      const double fpar = dotv(F[i], tau);
      if (R.climbing >= 0 && i == top) {
        for (size_t k = 0; k < n; ++k) G[i][k] = F[i][k] - tau[k] * (2 * fpar);
      } else {
        const double spring = o.spring * (std::sqrt(dotv(tp, tp)) - std::sqrt(dotv(tm, tm)));
        for (size_t k = 0; k < n; ++k) G[i][k] = F[i][k] - tau[k] * fpar + tau[k] * spring;
      }
      for (size_t k = 0; k < n; ++k) fmax = std::max(fmax, norm(G[i][k]));
    }
    return fmax;
  };
  // FIRE over every intermediate image (the parameters of caps/relax)
  const int nmin = 5;
  const double finc = 1.1, fdec = 0.5, astart = 0.1, falpha = 0.99, dtmax = 0.5, maxstep = 0.1;
  double dt = 0.05, alpha = astart;
  int npos = 0;
  std::vector<std::vector<Vec3>> v(m, std::vector<Vec3>(n, Vec3{0, 0, 0}));
  auto band_length = [&]() {
    double L = 0;
    for (size_t i = 1; i < m; ++i) {
      const auto dd = sub(path[i], path[i - 1]);
      L += std::sqrt(dotv(dd, dd));
    }
    return L;
  };
  const double length0 = band_length();
  bool diverged = false;
  double fmax = forces();
  int step = 0;
  for (; step < o.max_steps; ++step) {
    // A band whose length has grown tenfold, or whose energies are no longer finite, has left the path between the ends:
    // stop and say so rather than report a meaningless barrier.
    bool finite = std::isfinite(fmax);
    for (double e : E) finite = finite && std::isfinite(e);
    if (!finite || (length0 > 0 && band_length() > 10 * length0)) {
      diverged = true;
      break;
    }
    if (o.climb && !climbing && (fmax < std::max(o.climb_after, o.fmax) || step > o.max_steps / 2)) {
      climbing = true;
      fmax = forces();
      for (auto& x : v) std::fill(x.begin(), x.end(), Vec3{0, 0, 0});
      npos = 0, dt = 0.05, alpha = astart;
    }
    if (fmax < o.fmax && (climbing || !o.climb)) {
      R.converged = true;
      break;
    }
    double p = 0, vv = 0, ff = 0;
    for (size_t i = 1; i + 1 < m; ++i) p += dotv(G[i], v[i]), vv += dotv(v[i], v[i]), ff += dotv(G[i], G[i]);
    if (p > 0) {
      const double scale = alpha * std::sqrt(vv) / std::sqrt(std::max(ff, 1e-300));
      for (size_t i = 1; i + 1 < m; ++i)
        for (size_t k = 0; k < n; ++k) v[i][k] = v[i][k] * (1 - alpha) + G[i][k] * scale;
      if (++npos > nmin) dt = std::min(dt * finc, dtmax), alpha *= falpha;
    } else {
      for (auto& x : v) std::fill(x.begin(), x.end(), Vec3{0, 0, 0});
      dt *= fdec, alpha = astart, npos = 0;
    }
    double md = 0;
    for (size_t i = 1; i + 1 < m; ++i)
      for (size_t k = 0; k < n; ++k) {
        v[i][k] = v[i][k] + G[i][k] * dt;
        md = std::max(md, norm(v[i][k]) * dt);
      }
    const double sc = md > maxstep ? maxstep / md : 1.0;
    for (size_t i = 1; i + 1 < m; ++i)
      for (size_t k = 0; k < n; ++k) path[i][k] = path[i][k] + v[i][k] * (dt * sc);
    fmax = forces();
  }
  R.images = path;
  R.energy = E;
  R.fmax = fmax;
  R.steps = step;
  R.arc.assign(m, 0.0);
  for (size_t i = 1; i < m; ++i) {
    const auto dd = sub(path[i], path[i - 1]);
    R.arc[i] = R.arc[i - 1] + std::sqrt(dotv(dd, dd));
  }
  const double top = *std::max_element(E.begin() + 1, E.end() - 1);   // the highest intermediate image
  R.forward = top - E.front(), R.reverse = top - E.back();
  if (diverged)
    R.notes.push_back("the band diverged (grew past ten times its starting length or lost finite energies) and was stopped: "
                      "the barrier below is not meaningful; relax the ends or use more images");
  if (top < std::max(E.front(), E.back()))
    R.notes.push_back(o.climb && R.climbing < 0
                          ? "no image lies above both ends, so no image climbed: no barrier along this path (are the ends minima? relax them first)"
                          : "no image lies above both ends: no barrier along this path (are the ends minima?)");
  char b[200];
  std::snprintf(b, sizeof b, "%zu images, %d steps, largest NEB force %.4g%s; barrier %.6g from A, %.6g from B", m, step, fmax,
                R.converged ? "" : diverged ? " (not converged, band diverged)" : " (not converged)", R.forward, R.reverse);
  R.notes.push_back(b);
  return R;
}

// ---------------------------------------------------------------- interpolation

std::vector<std::vector<Vec3>> interpolate_path(const System& a, const System& b, const PathOptions& o, PathReport* report) {
  const size_t n = a.atoms.size();
  if (n != b.atoms.size()) throw std::invalid_argument("path: A and B have different atom counts");
  for (size_t i = 0; i < n; ++i)
    if (a.atoms[i].element != b.atoms[i].element) throw std::invalid_argument("path: atom " + std::to_string(i + 1) + " is another element in B (match the atoms first)");
  if (o.images < 1) throw std::invalid_argument("path: at least one image");
  if (o.method != "linear" && o.method != "lst" && o.method != "idpp") throw std::invalid_argument("path: linear, lst or idpp, not " + o.method);
  PathReport R;
  const bool periodic = a.cell.valid();
  std::vector<Vec3> A(n), B(n);
  for (size_t i = 0; i < n; ++i) {
    A[i] = a.atoms[i].pos;
    B[i] = periodic ? A[i] + a.cell.minimum_image(b.atoms[i].pos - A[i]) : b.atoms[i].pos;
  }
  if (!periodic && n >= 3) {   // a molecule: B turned and moved onto A, so the path carries no rigid motion
    const Superposition sp = superpose(A, B);
    for (auto& x : B) x = sp.apply(x);
    char line[120];
    std::snprintf(line, sizeof line, "B superposed on A (RMSD %.4f Å)", sp.rmsd);
    R.notes.push_back(line);
  }
  const int m = o.images + 2;
  std::vector<std::vector<Vec3>> path(static_cast<size_t>(m), std::vector<Vec3>(n));
  for (int k = 0; k < m; ++k) {
    const double f = double(k) / double(m - 1);
    for (size_t i = 0; i < n; ++i) path[size_t(k)][i] = A[i] * (1 - f) + B[i] * f;
  }
  auto sep = [&](const Vec3& x, const Vec3& y) { return periodic ? a.cell.minimum_image(y - x) : y - x; };
  if (o.method != "linear") {
    if (n > 3000) throw std::invalid_argument("path: lst and idpp use every atom pair; " + std::to_string(n) + " atoms is above 3000");
    std::vector<double> dA(n * n), dB(n * n);
    for (size_t i = 0; i < n; ++i)
      for (size_t j = i + 1; j < n; ++j) dA[i * n + j] = norm(sep(A[i], A[j])), dB[i * n + j] = norm(sep(B[i], B[j]));
    const std::vector<std::vector<Vec3>> li = path;
    // S_κ = Σ_{i<j} (d^κ − d)² / d⁴ with distances in bohr (+ γ Σ |r − r^LI|² for lst); f = −∂S/∂r per Å
    auto objective = [&](int k, const std::vector<Vec3>& x, std::vector<Vec3>& f, bool lst) {
      const double t = double(k) / double(m - 1);
      f.assign(n, Vec3{0, 0, 0});
      double S = 0;
      for (size_t i = 0; i < n; ++i)
        for (size_t j = i + 1; j < n; ++j) {
          const Vec3 r = sep(x[i], x[j]);
          const double rl = norm(r);
          if (rl < 1e-9) continue;
          const double d = rl / kBohr, target = (dA[i * n + j] + t * (dB[i * n + j] - dA[i * n + j])) / kBohr;
          const double d4 = d * d * d * d, e = target - d;
          S += e * e / d4;
          const double g = (-2 * e / d4 - 4 * e * e / (d4 * d)) / kBohr;   // ∂S/∂|r| in Å⁻¹ units
          const Vec3 u = r * (1 / rl);
          f[j] = f[j] - u * g, f[i] = f[i] + u * g;
        }
      if (lst) {
        const double gamma = 1e-6;   // bohr⁻⁴
        for (size_t i = 0; i < n; ++i) {
          const Vec3 dd = (x[i] - li[size_t(k)][i]) * (1 / kBohr);
          S += gamma * dot(dd, dd);
          f[i] = f[i] - dd * (2 * gamma / kBohr);
        }
      }
      return S;
    };
    if (o.method == "idpp") {
      NebOptions no;
      no.spring = 1.0;   // equal springs on the objective surface, as the paper
      no.climb = false;
      no.fmax = o.tolerance;
      no.max_steps = o.max_steps;
      const auto r = neb(path, [&](int k, const std::vector<Vec3>& x, std::vector<Vec3>& f) { return objective(k, x, f, false); }, no);
      path = r.images;
      R.steps = r.steps;
      if (!r.converged) R.notes.push_back("the IDPP band stopped at the step limit (largest force " + std::to_string(r.fmax) + ")");
    } else {   // each image on its own: FIRE on S_LST from its linear image
      for (int k = 1; k + 1 < m; ++k) {
        std::vector<Vec3> x = li[size_t(k)], f, v(n, Vec3{0, 0, 0});
        objective(k, x, f, true);
        double dt = 0.05, alpha = 0.1;
        int npos = 0, it = 0;
        for (; it < o.max_steps; ++it) {
          double fm = 0, p = 0, vv = 0, ff = 0;
          for (size_t i = 0; i < n; ++i) fm = std::max(fm, norm(f[i])), p += dot(f[i], v[i]), vv += dot(v[i], v[i]), ff += dot(f[i], f[i]);
          if (fm < o.tolerance) break;
          if (p > 0) {
            const double scale = alpha * std::sqrt(vv / std::max(ff, 1e-300));
            for (size_t i = 0; i < n; ++i) v[i] = v[i] * (1 - alpha) + f[i] * scale;
            if (++npos > 5) dt = std::min(dt * 1.1, 0.5), alpha *= 0.99;
          } else {
            std::fill(v.begin(), v.end(), Vec3{0, 0, 0});
            dt *= 0.5, alpha = 0.1, npos = 0;
          }
          double md = 0;
          for (size_t i = 0; i < n; ++i) v[i] = v[i] + f[i] * dt, md = std::max(md, norm(v[i]) * dt);
          const double sc = md > 0.1 ? 0.1 / md : 1.0;
          for (size_t i = 0; i < n; ++i) x[i] = x[i] + v[i] * (dt * sc);
          objective(k, x, f, true);
        }
        R.steps = std::max(R.steps, it);
        if (it >= o.max_steps) R.notes.push_back("LST image " + std::to_string(k) + " stopped at the step limit");
        path[size_t(k)] = x;
      }
    }
  }
  if (n <= 3000)
    for (const auto& img : path) {
      double dmin = std::numeric_limits<double>::infinity();
      for (size_t i = 0; i < n; ++i)
        for (size_t j = i + 1; j < n; ++j) dmin = std::min(dmin, norm(sep(img[i], img[j])));
      R.min_distance.push_back(n > 1 ? dmin : 0.0);
    }
  char line[160];
  std::snprintf(line, sizeof line, "%d images by %s interpolation%s", o.images, o.method.c_str(), periodic ? " (minimum image)" : "");
  R.notes.insert(R.notes.begin(), line);
  if (report) *report = R;
  return path;
}

NebResult neb_forcefield(const System& a, std::vector<std::vector<Vec3>> path, const ForceField& ff, const EnergyOptions& eo, const NebOptions& o) {
  Evaluator ev(ff, eo);
  std::vector<double> x, f;
  const ImageEnergy en = [&](int, const std::vector<Vec3>& p, std::vector<Vec3>& out) {
    x.resize(3 * p.size());
    for (size_t i = 0; i < p.size(); ++i) x[3 * i] = p[i][0], x[3 * i + 1] = p[i][1], x[3 * i + 2] = p[i][2];
    const EnergyTerms t = ev.compute(x, a.cell, f);
    out.resize(p.size());
    for (size_t i = 0; i < p.size(); ++i) out[i] = {f[3 * i], f[3 * i + 1], f[3 * i + 2]};
    return t.total();
  };
  NebResult r = neb(std::move(path), en, o);
  r.notes.push_back("energies in kcal/mol from the force field");
  return r;
}

}  // namespace caps
