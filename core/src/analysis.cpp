#include "caps/analysis.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "caps/elements.hpp"
#include "cell_list.hpp"

namespace caps {

void symmetric_eigen3(const double Ain[3][3], double w[3], double V[3][3]) {
  double A[3][3];
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) { A[i][j] = Ain[i][j]; V[i][j] = i == j; }
  for (int sweep = 0; sweep < 60; ++sweep) {
    int p = 0, q = 1;
    double big = std::fabs(A[0][1]);
    if (std::fabs(A[0][2]) > big) { p = 0; q = 2; big = std::fabs(A[0][2]); }
    if (std::fabs(A[1][2]) > big) { p = 1; q = 2; big = std::fabs(A[1][2]); }
    if (big < 1e-14) break;
    const double th = 0.5 * std::atan2(2 * A[p][q], A[q][q] - A[p][p]);
    const double c = std::cos(th), s = std::sin(th);
    for (int k = 0; k < 3; ++k) {
      const double akp = A[k][p], akq = A[k][q];
      A[k][p] = c * akp - s * akq;
      A[k][q] = s * akp + c * akq;
    }
    for (int k = 0; k < 3; ++k) {
      const double apk = A[p][k], aqk = A[q][k];
      A[p][k] = c * apk - s * aqk;
      A[q][k] = s * apk + c * aqk;
    }
    for (int k = 0; k < 3; ++k) {
      const double vkp = V[k][p], vkq = V[k][q];
      V[k][p] = c * vkp - s * vkq;
      V[k][q] = s * vkp + c * vkq;
    }
  }
  int order[3] = {0, 1, 2};
  std::sort(order, order + 3, [&](int x, int y) { return A[x][x] < A[y][y]; });
  double Vs[3][3];
  for (int k = 0; k < 3; ++k) {
    w[k] = A[order[k]][order[k]];
    for (int r = 0; r < 3; ++r) Vs[r][k] = V[r][order[k]];
  }
  for (int r = 0; r < 3; ++r)
    for (int k = 0; k < 3; ++k) V[r][k] = Vs[r][k];
}

double measure(const System& s, const std::vector<uint32_t>& idx) {
  if (idx.size() < 2 || idx.size() > 4) throw std::invalid_argument("measure needs 2, 3 or 4 atoms");
  for (auto i : idx)
    if (i >= s.atoms.size()) throw std::out_of_range("atom index out of range");
  // Chain of minimum-image steps so the atoms form one connected geometry.
  std::vector<Vec3> p{s.atoms[idx[0]].pos};
  for (size_t k = 1; k < idx.size(); ++k) p.push_back(p.back() + s.cell.minimum_image(s.atoms[idx[k]].pos - s.atoms[idx[k - 1]].pos));
  if (p.size() == 2) return norm(p[1] - p[0]);
  if (p.size() == 3) {
    const Vec3 a = p[0] - p[1], b = p[2] - p[1];
    return std::acos(std::clamp(dot(a, b) / (norm(a) * norm(b)), -1.0, 1.0)) * 180.0 / M_PI;
  }
  const Vec3 b1 = p[1] - p[0], b2 = p[2] - p[1], b3 = p[3] - p[2];
  const Vec3 n1 = cross(b1, b2), n2 = cross(b2, b3);
  // IUPAC: positive when the near bond turns clockwise onto the far one, looking along b2
  return std::atan2(dot(cross(n1, n2), b2 * (1.0 / norm(b2))), dot(n1, n2)) * 180.0 / M_PI;
}

std::vector<MoleculeShape> molecule_shapes(const System& s) {
  int nm = 0;
  const auto mol = s.molecules(&nm);
  std::vector<MoleculeShape> out(nm);
  std::vector<double> m(nm, 0);
  std::vector<Vec3> com(nm, {0, 0, 0});
  for (size_t i = 0; i < s.atoms.size(); ++i) {
    const double w = s.mass_of(s.atoms[i]);
    m[mol[i]] += w;
    com[mol[i]] = com[mol[i]] + s.atoms[i].pos * w;
    out[mol[i]].atoms++;
  }
  std::vector<std::array<double, 6>> G(nm, {0, 0, 0, 0, 0, 0});
  for (int k = 0; k < nm; ++k) com[k] = com[k] * (1.0 / m[k]);
  for (size_t i = 0; i < s.atoms.size(); ++i) {
    const int k = mol[i];
    const double w = s.mass_of(s.atoms[i]);
    const Vec3 d = s.atoms[i].pos - com[k];
    G[k][0] += w * d[0] * d[0]; G[k][1] += w * d[1] * d[1]; G[k][2] += w * d[2] * d[2];
    G[k][3] += w * d[0] * d[1]; G[k][4] += w * d[0] * d[2]; G[k][5] += w * d[1] * d[2];
  }
  for (int k = 0; k < nm; ++k) {
    const double A[3][3] = {{G[k][0] / m[k], G[k][3] / m[k], G[k][4] / m[k]},
                            {G[k][3] / m[k], G[k][1] / m[k], G[k][5] / m[k]},
                            {G[k][4] / m[k], G[k][5] / m[k], G[k][2] / m[k]}};
    double w[3], V[3][3];
    symmetric_eigen3(A, w, V);
    auto& o = out[k];
    o.molecule = k;
    o.mass = m[k];
    o.com = com[k];
    for (int j = 0; j < 3; ++j) o.lambda[j] = w[j], o.axis[j] = {V[0][j], V[1][j], V[2][j]};
    const double t = w[0] + w[1] + w[2];
    o.rg = std::sqrt(std::max(0.0, t));
    o.kappa2 = t > 0 ? 1 - 3 * (w[0] * w[1] + w[1] * w[2] + w[2] * w[0]) / (t * t) : 0;
  }
  return out;
}

std::vector<std::pair<double, double>> rdf(const System& s, int ea, int eb, double rmax, double dr, bool inter_only) {
  // !(x > 0) also rejects NaN; a bin count past a million is a mistyped --dr, not a histogram
  if (!(dr > 0) || !std::isfinite(dr)) throw std::invalid_argument("rdf: the bin width dr must be a number greater than 0 (got " + std::to_string(dr) + ")");
  if (!(rmax > 0) || !std::isfinite(rmax)) throw std::invalid_argument("rdf: r_max must be a number greater than 0 (got " + std::to_string(rmax) + ")");
  if (rmax / dr > 1e6) throw std::invalid_argument("rdf: r_max / dr gives more than a million bins; use a larger dr");
  const int nb = static_cast<int>(rmax / dr);
  std::vector<double> h(nb, 0);
  std::vector<uint32_t> A, B;
  for (uint32_t i = 0; i < s.atoms.size(); ++i) {
    if (!ea || s.atoms[i].element == ea) A.push_back(i);
    if (!eb || s.atoms[i].element == eb) B.push_back(i);
  }
  // Large systems: an evenly strided subset of the centres (every neighbour still counts); the average over 20 000
  // centres is already smooth, and the cost stays bounded for million-atom cells.
  if (A.size() > kRdfMaxCentres) {
    const double step = double(A.size()) / double(kRdfMaxCentres);
    std::vector<uint32_t> sub;
    sub.reserve(kRdfMaxCentres);
    for (size_t k = 0; k < kRdfMaxCentres; ++k) sub.push_back(A[size_t(k * step)]);
    A = std::move(sub);
  }
  const auto mol = inter_only ? s.molecules() : std::vector<int>();
  std::vector<char> inB(s.atoms.size(), 0);
  for (uint32_t j : B) inB[j] = 1;
  // Cell list when rmax fits inside half the cell (exact minimum image); all pairs otherwise.
  const double v = s.cell.valid() ? s.cell.volume() : 0;
  const double wmin = v > 0 ? std::min({v / norm(cross(s.cell.b, s.cell.c)), v / norm(cross(s.cell.c, s.cell.a)), v / norm(cross(s.cell.a, s.cell.b))}) : 0;
  if (v > 0 && rmax < wmin / 2) {
    Grid g(s, rmax);
    for (uint32_t i : A) {
      const Vec3& f = g.frac[i];
      g.for_neighbour_bins(g.bin(f, 0), g.bin(f, 1), g.bin(f, 2), [&](const std::vector<uint32_t>& bin) {
        for (uint32_t j : bin) {
          if (j == i || !inB[j]) continue;
          if (inter_only && mol[i] == mol[j]) continue;
          const double r = norm(g.sep(i, j));
          if (r < rmax) { const int k = static_cast<int>(r / dr); if (k < nb) h[k] += 1; }   // rmax / dr rounds down
        }
      });
    }
  } else if (v > 0) {
    // rmax reaches past half the cell: every periodic image within rmax counts, not only the nearest one (a small
    // crystal cell has most of its shells in the neighbouring cells, and an atom's own images are neighbours too)
    const Cell& c = s.cell;
    const Vec3 ra = cross(c.b, c.c) * (1 / v), rb = cross(c.c, c.a) * (1 / v), rc = cross(c.a, c.b) * (1 / v);
    const Vec3 rec[3] = {ra, rb, rc}, lat[3] = {c.a, c.b, c.c};
    double reach[3], images = 1;
    for (int k = 0; k < 3; ++k) {
      reach[k] = c.periodic[size_t(k)] ? rmax * norm(rec[k]) : 0;   // rmax over the spacing of the lattice planes
      images *= 2 * std::ceil(reach[k]) + 1;
    }
    if (double(A.size()) * double(B.size()) * images > 4e9)   // a few seconds per frame
      throw std::invalid_argument("rdf: r_max = " + std::to_string(rmax) + " Å reaches about " + std::to_string(long(images)) +
                                  " periodic images of this cell for every pair, too many to count; use a smaller r_max");
    for (uint32_t i : A)
      for (uint32_t j : B) {
        if (inter_only && mol[i] == mol[j]) continue;
        Vec3 d = s.atoms[j].pos - s.atoms[i].pos;
        int lo[3], hi[3];
        for (int k = 0; k < 3; ++k) {
          double f = dot(d, rec[k]);
          if (c.periodic[size_t(k)]) {
            const double sh = std::round(f);
            f -= sh;
            d = d - lat[k] * sh;
          }
          // images t with |f + t| ≤ reach along this axis (a point within rmax is within rmax of every plane through i)
          lo[k] = int(std::ceil(-reach[k] - f));
          hi[k] = int(std::floor(reach[k] - f));
          if (!c.periodic[size_t(k)]) lo[k] = hi[k] = 0;
        }
        for (int t0 = lo[0]; t0 <= hi[0]; ++t0)
          for (int t1 = lo[1]; t1 <= hi[1]; ++t1)
            for (int t2 = lo[2]; t2 <= hi[2]; ++t2) {
              if (i == j && t0 == 0 && t1 == 0 && t2 == 0) continue;
              const double r = norm(d + lat[0] * double(t0) + lat[1] * double(t1) + lat[2] * double(t2));
              if (r < rmax) { const int k = static_cast<int>(r / dr); if (k < nb) h[k] += 1; }   // rmax / dr rounds down
            }
      }
  } else {
    // no cell: plain distances (the density below is zero, so g(r) needs a periodic cell to mean anything)
    for (uint32_t i : A)
      for (uint32_t j : B) {
        if (i == j) continue;
        if (inter_only && mol[i] == mol[j]) continue;
        const double r = norm(s.atoms[j].pos - s.atoms[i].pos);
        if (r < rmax) { const int k = static_cast<int>(r / dr); if (k < nb) h[k] += 1; }   // rmax / dr rounds down
      }
  }
  const double rho = s.cell.valid() ? B.size() / s.cell.volume() : 0;
  std::vector<std::pair<double, double>> g(nb);
  for (int k = 0; k < nb; ++k) {
    const double r0 = k * dr, r1 = r0 + dr;
    const double shell = 4.0 / 3.0 * M_PI * (r1 * r1 * r1 - r0 * r0 * r0);
    g[k] = {r0 + dr / 2, (A.empty() || rho == 0) ? 0 : h[k] / (A.size() * rho * shell)};
  }
  return g;
}

namespace {

// Bonds that lie on a cycle, by Tarjan's bridge search (iterative). in_ring[i] is set for atoms with a cycle bond.
std::vector<char> ring_atoms(const std::vector<std::vector<uint32_t>>& nb, const std::vector<char>& use) {
  const size_t n = nb.size();
  std::vector<int> disc(n, -1), low(n, 0);
  std::vector<char> ring(n, 0);
  int timer = 0;
  struct Frame { uint32_t v; int parent; size_t k; };
  for (uint32_t r = 0; r < n; ++r) {
    if (!use[r] || disc[r] >= 0) continue;
    std::vector<Frame> st{{r, -1, 0}};
    disc[r] = low[r] = timer++;
    while (!st.empty()) {
      Frame& f = st.back();
      if (f.k < nb[f.v].size()) {
        const uint32_t w = nb[f.v][f.k++];
        if (!use[w] || int(w) == f.parent) continue;
        if (disc[w] < 0) {
          disc[w] = low[w] = timer++;
          st.push_back({w, int(f.v), 0});
        } else {
          low[f.v] = std::min(low[f.v], disc[w]);
        }
      } else {
        const uint32_t v = f.v;
        const int p = f.parent;
        st.pop_back();
        if (p >= 0) {
          low[p] = std::min(low[p], low[v]);
          if (low[v] <= disc[p]) { ring[v] = 1; ring[p] = 1; }   // edge p–v is not a bridge: it is on a cycle
        }
      }
    }
  }
  return ring;
}

// Longest shortest path (double breadth-first search) in the subgraph of atoms with use[i] from a start atom.
std::vector<uint32_t> diameter_path(const std::vector<std::vector<uint32_t>>& nb, const std::vector<char>& use, uint32_t start) {
  auto bfs = [&](uint32_t s0, std::vector<int>& parent) {
    std::vector<uint32_t> q{s0};
    parent.assign(nb.size(), -2);
    parent[s0] = -1;
    uint32_t last = s0;
    for (size_t h = 0; h < q.size(); ++h) {
      last = q[h];
      for (uint32_t w : nb[last])
        if (use[w] && parent[w] == -2) { parent[w] = int(last); q.push_back(w); }
    }
    return last;
  };
  std::vector<int> parent;
  const uint32_t a = bfs(start, parent);
  const uint32_t b = bfs(a, parent);
  std::vector<uint32_t> path;
  for (int v = int(b); v >= 0; v = parent[v]) path.push_back(uint32_t(v));
  return path;
}

// Sum of the minimum-image bond vectors around a closed loop of bonded atoms (last atom bonded to the first). Zero for
// an ordinary ring; a lattice vector when the loop is a chain that closes on itself through the periodic cell.
Vec3 loop_closure(const System& s, const std::vector<uint32_t>& loop) {
  Vec3 sum{0, 0, 0};
  for (size_t k = 0; k < loop.size(); ++k) {
    const uint32_t a = loop[k], b = loop[(k + 1) % loop.size()];
    sum = sum + s.cell.minimum_image(s.atoms[b].pos - s.atoms[a].pos);
  }
  return sum;
}

// A loop whose closure is at least this long (Å) runs through the cell; a real ring closes to within rounding.
constexpr double kPeriodicClosure = 1.0;

}  // namespace

std::vector<std::vector<uint32_t>> backbones(const System& s, int min_atoms) {
  const size_t n = s.atoms.size();
  const auto nb = s.neighbours();
  std::vector<char> heavy(n);
  for (size_t i = 0; i < n; ++i) heavy[i] = s.atoms[i].element != 1;
  const auto ring = ring_atoms(nb, heavy);
  std::vector<char> chain(n);
  for (size_t i = 0; i < n; ++i) chain[i] = heavy[i] && !ring[i];
  const auto mol = s.molecules();
  std::vector<std::vector<uint32_t>> out;
  // Chains bonded through the cell (a crystal's infinite chains): in a component of heavy atoms whose cycle atoms form
  // one simple cycle (every cycle atom with two cycle neighbours; side groups such as iPP's methyls hang off it) and
  // whose bond vectors add up to a lattice vector, the cycle is walked in order and cut once, and the whole component
  // is used up. A ring or a cyclic molecule closes on itself and is left to the ordinary search below.
  std::vector<char> in_cycle(n, 0);
  if (s.cell.valid()) {
    std::vector<char> seen(n, 0);
    auto on_cycle = [&](uint32_t w) { return heavy[w] && ring[w]; };
    auto cycle_deg = [&](uint32_t a) { int k = 0; for (uint32_t w : nb[a]) k += on_cycle(w) ? 1 : 0; return k; };
    for (uint32_t i = 0; i < n; ++i) {
      if (!heavy[i] || seen[i]) continue;
      std::vector<uint32_t> comp{i};
      seen[i] = 1;
      for (size_t h = 0; h < comp.size(); ++h)
        for (uint32_t w : nb[comp[h]])
          if (heavy[w] && !seen[w]) { seen[w] = 1; comp.push_back(w); }
      std::vector<uint32_t> cyc;
      for (uint32_t a : comp) if (ring[a]) cyc.push_back(a);
      if (cyc.size() < 3 || !std::all_of(cyc.begin(), cyc.end(), [&](uint32_t a) { return cycle_deg(a) == 2; })) continue;
      std::vector<uint32_t> path{cyc[0]};
      uint32_t prev = cyc[0], cur = cyc[0];
      for (uint32_t w : nb[cur]) if (on_cycle(w)) { cur = w; break; }
      while (cur != cyc[0] && path.size() <= cyc.size()) {
        path.push_back(cur);
        uint32_t next = cur;
        for (uint32_t w : nb[cur]) if (on_cycle(w) && w != prev) { next = w; break; }
        prev = cur;
        cur = next;
      }
      if (path.size() != cyc.size()) continue;
      const Vec3 closure = loop_closure(s, path);
      if (dot(closure, closure) < kPeriodicClosure * kPeriodicClosure) continue;
      for (uint32_t a : comp) in_cycle[a] = 1;
      if (int(path.size()) >= min_atoms) out.push_back(std::move(path));
    }
  }
  // The other heavy atoms, grouped by molecule and, within a molecule, by bonded piece: a crystal or data file that puts
  // every chain in one molecule ID still gives one backbone per chain. Paths stay inside their group.
  std::vector<char> grouped(n, 0), in_group(n, 0), group_chain(n, 0), seen(n, 0);
  for (uint32_t i = 0; i < n; ++i) {
    if (!heavy[i] || in_cycle[i] || grouped[i]) continue;
    std::vector<uint32_t> m{i};
    grouped[i] = 1;
    for (size_t h = 0; h < m.size(); ++h)
      for (uint32_t w : nb[m[h]])
        if (heavy[w] && !in_cycle[w] && !grouped[w] && mol[w] == mol[i]) { grouped[w] = 1; m.push_back(w); }
    if (int(m.size()) < min_atoms) continue;
    for (uint32_t a : m) { in_group[a] = 1; group_chain[a] = chain[a]; }
    // longest non-ring path over the group's non-ring components
    std::vector<uint32_t> best;
    for (uint32_t a : m) {
      if (!group_chain[a] || seen[a]) continue;
      auto p = diameter_path(nb, group_chain, a);
      // mark the component
      std::vector<uint32_t> q{a};
      seen[a] = 1;
      for (size_t h = 0; h < q.size(); ++h)
        for (uint32_t w : nb[q[h]])
          if (group_chain[w] && !seen[w]) { seen[w] = 1; q.push_back(w); }
      if (p.size() > best.size()) best = std::move(p);
    }
    const auto full = diameter_path(nb, in_group, m[0]);
    if (best.size() * 2 < full.size()) best = full;   // rings in the main chain
    for (uint32_t a : m) { in_group[a] = 0; group_chain[a] = 0; }
    if (int(best.size()) >= min_atoms) out.push_back(std::move(best));
  }
  return out;
}

InternalDistances internal_distances(const System& s) {
  InternalDistances r;
  const auto bb = backbones(s);
  r.chains = int(bb.size());
  if (bb.empty()) return r;
  size_t lmax = 0;
  for (const auto& p : bb) lmax = std::max(lmax, p.size());
  std::vector<double> sum(lmax, 0.0);
  std::vector<long> cnt(lmax, 0);
  double b2 = 0;
  long nb = 0;
  std::vector<Vec3> x;
  std::vector<std::vector<uint32_t>> nbr;
  for (const auto& p : bb) {
    // Consecutive backbone atoms are bonded: chain minimum-image bond vectors so a chain that crosses the cell
    // (or runs through it as an infinite crystal chain) is measured unwrapped, with no bond spanning the box.
    x.assign(1, s.atoms[p[0]].pos);
    for (size_t k = 1; k < p.size(); ++k) x.push_back(x.back() + s.cell.minimum_image(s.atoms[p[k]].pos - s.atoms[p[k - 1]].pos));
    for (size_t i = 0; i < p.size(); ++i)
      for (size_t j = i + 1; j < p.size(); ++j) {
        const Vec3 d = x[j] - x[i];
        sum[j - i] += dot(d, d);
        ++cnt[j - i];
      }
    r.r2_end += dot(x.back() - x.front(), x.back() - x.front());
    // a chain that closes on itself through the cell (last atom bonded to the first) is infinite: its end-to-end
    // distance is one period of the crystal, not a coil size
    if (p.size() >= 3) {
      if (nbr.empty()) nbr = s.neighbours();
      if (std::find(nbr[p.front()].begin(), nbr[p.front()].end(), p.back()) != nbr[p.front()].end()) {
        const Vec3 c = loop_closure(s, p);
        if (dot(c, c) >= kPeriodicClosure * kPeriodicClosure) ++r.periodic;
      }
    }
  }
  b2 = cnt[1] ? sum[1] / cnt[1] : 0;
  nb = cnt[1];
  (void)nb;
  r.b2 = b2;
  r.r2_end /= bb.size();
  for (size_t k = 1; k < lmax; ++k)
    if (cnt[k] > 0 && b2 > 0) {
      r.n.push_back(int(k));
      r.ratio.push_back(sum[k] / cnt[k] / (k * b2));
    }
  return r;
}

}  // namespace caps
