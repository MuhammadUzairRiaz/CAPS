// CAPS sorption: Widom insertion, grand-canonical and fixed-loading Monte Carlo in a fixed host, rigid or flexible
// sorbates (see sorption.hpp).
#include "caps/rng.hpp"
#include "caps/sorption.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <limits>
#include <map>
#include <random>
#include <stdexcept>

#include "caps/pairmodel.hpp"
#include "caps/typing.hpp"

namespace caps {

namespace {

constexpr double kB = 0.0019872043;          // kcal/(mol·K)
constexpr double kBoltzmannJ = 1.380649e-23; // J/K
constexpr double kAvogadro = 6.02214076e23;
constexpr double kPi = 3.14159265358979323846;
constexpr double kInf = std::numeric_limits<double>::infinity();

struct Quat { double w, x, y, z; };
Vec3 rotate(const Quat& q, const Vec3& v) {   // q v q*
  const Vec3 u{q.x, q.y, q.z};
  const Vec3 t = cross(u, v) * 2.0;
  return v + t * q.w + cross(u, t);
}
// p turned by angle t about the axis through c along unit u (Rodrigues)
Vec3 rotate_about(const Vec3& p, const Vec3& c, const Vec3& u, double t) {
  const Vec3 v = p - c;
  const double ct = std::cos(t), st = std::sin(t);
  return c + v * ct + cross(u, v) * st + u * (dot(u, v) * (1 - ct));
}
// ln Σ exp(l) over finite entries (−∞ when none)
double log_sum_exp(const std::vector<double>& l) {
  double m = -kInf;
  for (double v : l) m = std::max(m, v);
  if (!std::isfinite(m)) return -kInf;
  double s = 0;
  for (double v : l) s += std::exp(v - m);
  return m + std::log(s);
}
Vec3 centroid(const std::vector<Vec3>& P) {
  Vec3 c{0, 0, 0};
  for (const auto& p : P) c = c + p;
  return c * (1.0 / double(std::max<size_t>(1, P.size())));
}

// The fixed host atoms binned (bins at least rc/2 wide, two either side), and the pair model.
// A cut-off wider than half the cell's narrowest width reaches more than one periodic image of an atom, and the
// minimum image then drops pairs inside the cut-off (a shorter cut-off instead would bias the Henry constant and the
// loadings). So past that width every image within the cut-off counts: the host's images are ghost copies binned in
// space (a sorbate atom folded into the cell meets them at plain distances), and sorbate pairs are lattice sums over
// the translations T with |d + T| < rc. The energy per cell of the periodic system is then
// E = sum_i u_host(i) + 1/2 sum_{i,j} sum_T' u(r_ij + T) (Allen & Tildesley, Computer Simulation of Liquids, sec. 1.5),
// so a molecule also meets its own images: 1/2 sum_{k,q} sum_{T != 0} u(r_kq + T), see own().
struct Host {
  const System& s;
  const PairModel& pm;
  bool ortho = false;
  double L[3] = {0, 0, 0};
  int nb[3] = {1, 1, 1};
  std::vector<std::vector<uint32_t>> bins;
  std::vector<std::vector<int>> stencil;
  bool images = false;
  double w[3] = {0, 0, 0}, wmin = 0;   // the cell's widths across each pair of faces (Å), and the narrowest
  Vec3 rec[3];                         // fractional height of a separation d across faces k: dot(rec[k], d)
  std::vector<Vec3> gpos;              // the host's images within the cut-off of the cell
  std::vector<uint32_t> gatom;         // ... and the atom each copies
  Vec3 glo{0, 0, 0};
  double gw[3] = {1, 1, 1};
  int gn[3] = {1, 1, 1};
  std::vector<std::vector<uint32_t>> gbins;
  Host(const System& sys, const PairModel& p, size_t n_host) : s(sys), pm(p) {
    const Cell& c = sys.cell;
    ortho = std::fabs(c.a[1]) + std::fabs(c.a[2]) + std::fabs(c.b[0]) + std::fabs(c.b[2]) + std::fabs(c.c[0]) + std::fabs(c.c[1]) < 1e-9;
    L[0] = c.a[0], L[1] = c.b[1], L[2] = c.c[2];
    const double V = c.volume(), det = dot(c.a, cross(c.b, c.c));
    w[0] = V / norm(cross(c.b, c.c)), w[1] = V / norm(cross(c.c, c.a)), w[2] = V / norm(cross(c.a, c.b));
    rec[0] = cross(c.b, c.c) * (1.0 / det), rec[1] = cross(c.c, c.a) * (1.0 / det), rec[2] = cross(c.a, c.b) * (1.0 / det);
    wmin = std::min({w[0], w[1], w[2]});
    images = p.cutoff() > 0.5 * wmin;
    if (images) build_images(n_host);
    for (int d = 0; d < 3; ++d) nb[d] = std::max(1, int(std::floor(2 * w[d] / p.cutoff())));
    bins.resize(size_t(nb[0] * nb[1] * nb[2]));
    for (uint32_t i = 0; i < n_host; ++i) bins[size_t(bin_of(sys.atoms[i].pos))].push_back(i);
    stencil.resize(bins.size());
    for (int a = 0; a < nb[0]; ++a)
      for (int b = 0; b < nb[1]; ++b)
        for (int cc = 0; cc < nb[2]; ++cc) {
          auto offs = [&](int d, int v) {
            std::vector<int> r;
            if (nb[d] < 5) for (int k = 0; k < nb[d]; ++k) r.push_back(k);
            else for (int dv = -2; dv <= 2; ++dv) r.push_back(((v + dv) % nb[d] + nb[d]) % nb[d]);
            return r;
          };
          for (int i : offs(0, a))
            for (int j : offs(1, b))
              for (int k : offs(2, cc)) stencil[size_t((a * nb[1] + b) * nb[2] + cc)].push_back((i * nb[1] + j) * nb[2] + k);
        }
  }
  void build_images(size_t n_host) {
    const double rc = pm.cutoff();
    double per_pair = 1;   // a bounded budget: a cut-off many times the cell's width needs millions of images
    for (int d = 0; d < 3; ++d) per_pair *= 1 + 2 * rc / w[d];
    const double est = per_pair * double(n_host);
    if (est > 4e6 || per_pair > 2e5) {
      char b[320];
      std::snprintf(b, sizeof b, "sorption: a %.1f Å cut-off in a cell %.2f Å across needs about %.0f periodic images; use a shorter cut-off or a larger cell",
                    rc, wmin, std::max(est, per_pair));
      throw std::invalid_argument(b);
    }
    int n[3];
    double reach[3];
    for (int d = 0; d < 3; ++d) n[d] = int(std::ceil(rc / w[d])), reach[d] = rc / w[d];
    for (uint32_t a = 0; a < n_host; ++a) {
      Vec3 f = s.cell.to_fractional(s.atoms[a].pos);
      for (int d = 0; d < 3; ++d) f[size_t(d)] -= std::floor(f[size_t(d)]);
      for (int i = -n[0]; i <= n[0]; ++i)
        for (int j = -n[1]; j <= n[1]; ++j)
          for (int k = -n[2]; k <= n[2]; ++k) {
            const Vec3 g{f[0] + i, f[1] + j, f[2] + k};
            bool keep = true;   // within the cut-off of the cell (a point in it is at most rc / w away across each pair of faces)
            for (int d = 0; d < 3 && keep; ++d) keep = g[size_t(d)] > -reach[d] && g[size_t(d)] < 1 + reach[d];
            if (keep) gpos.push_back(s.cell.to_cartesian(g)), gatom.push_back(a);
          }
    }
    // bins at least rc wide over the images' box; a point outside the box searches the bins at its edge
    Vec3 hi{-1e30, -1e30, -1e30};
    glo = {1e30, 1e30, 1e30};
    for (const auto& q : gpos)
      for (size_t d = 0; d < 3; ++d) glo[d] = std::min(glo[d], q[d]), hi[d] = std::max(hi[d], q[d]);
    if (gpos.empty()) glo = {0, 0, 0};
    for (int d = 0; d < 3; ++d) {
      const double ext = gpos.empty() ? 0 : hi[size_t(d)] - glo[size_t(d)];
      gn[d] = std::clamp(int(std::floor(ext / rc)), 1, 64);
      gw[d] = std::max(rc, ext / gn[d]);
    }
    gbins.assign(size_t(gn[0] * gn[1] * gn[2]), {});
    for (uint32_t g = 0; g < gpos.size(); ++g) {
      int id[3];
      gbin(gpos[g], id);
      gbins[size_t((id[0] * gn[1] + id[1]) * gn[2] + id[2])].push_back(g);
    }
  }
  void gbin(const Vec3& p, int id[3]) const {
    for (int d = 0; d < 3; ++d) id[d] = std::clamp(int(std::floor((p[size_t(d)] - glo[size_t(d)]) / gw[d])), 0, gn[d] - 1);
  }
  // pm.energy of atoms ta and tb over the lattice translations T with |d + T| inside the cut-off; skip0: T = 0 left out
  double lattice_sum(uint32_t ta, uint32_t tb, const Vec3& d, bool skip0) const {
    const double rc = pm.cutoff();
    const Cell& c = s.cell;
    int lo[3], hi[3];
    for (int t = 0; t < 3; ++t) {
      const double f = dot(rec[t], d);   // (f + i) w_t must lie within ±rc
      lo[t] = int(std::ceil(-rc / w[t] - f)), hi[t] = int(std::floor(rc / w[t] - f));
    }
    double e = 0;
    for (int i = lo[0]; i <= hi[0]; ++i)
      for (int j = lo[1]; j <= hi[1]; ++j)
        for (int k = lo[2]; k <= hi[2]; ++k) {
          if (skip0 && i == 0 && j == 0 && k == 0) continue;
          const Vec3 r = d + c.a * double(i) + c.b * double(j) + c.c * double(k);
          e += pm.energy(ta, tb, dot(r, r));
        }
    return e;
  }
  // two sorbate atoms (force-field indices ta, tb) a separation d apart, every image within the cut-off
  double pair(uint32_t ta, uint32_t tb, const Vec3& d) const {
    if (images) return lattice_sum(ta, tb, d, false);
    const Vec3 m = mi(d);
    return pm.energy(ta, tb, dot(m, m));
  }
  // a molecule (atoms tpl at P) with its own periodic images, each pair of images once; 0 while the cut-off reaches
  // one image only
  double own(const std::vector<uint32_t>& tpl, const std::vector<Vec3>& P) const {
    if (!images) return 0;
    double e = 0;
    for (size_t k = 0; k < tpl.size(); ++k) {
      e += 0.5 * lattice_sum(tpl[k], tpl[k], Vec3{0, 0, 0}, true);
      for (size_t q = k + 1; q < tpl.size(); ++q) e += lattice_sum(tpl[k], tpl[q], P[q] - P[k], true);
    }
    return e;
  }
  Vec3 mi(Vec3 d) const {
    if (!ortho) return s.cell.minimum_image(d);
    for (int k = 0; k < 3; ++k) d[size_t(k)] -= L[k] * std::round(d[size_t(k)] / L[k]);
    return d;
  }
  int bin_of(const Vec3& p) const {
    const Vec3 f = s.cell.to_fractional(p);
    int id[3];
    for (int d = 0; d < 3; ++d) id[d] = std::min(nb[d] - 1, int((f[size_t(d)] - std::floor(f[size_t(d)])) * nb[d]));
    return (id[0] * nb[1] + id[1]) * nb[2] + id[2];
  }
  // one atom (force-field index t) at p with the host
  double atom(uint32_t t, const Vec3& p) const {
    double e = 0;
    if (images) {
      const Vec3 q = s.cell.wrap(p);
      int id[3];
      gbin(q, id);
      for (int a = std::max(0, id[0] - 1); a <= std::min(gn[0] - 1, id[0] + 1); ++a)
        for (int b = std::max(0, id[1] - 1); b <= std::min(gn[1] - 1, id[1] + 1); ++b)
          for (int c = std::max(0, id[2] - 1); c <= std::min(gn[2] - 1, id[2] + 1); ++c)
            for (uint32_t g : gbins[size_t((a * gn[1] + b) * gn[2] + c)]) {
              const Vec3 d = gpos[g] - q;
              e += pm.energy(t, gatom[g], dot(d, d));
            }
      return e;
    }
    for (int bn : stencil[size_t(bin_of(p))])
      for (uint32_t b : bins[size_t(bn)]) {
        const Vec3 d = mi(s.atoms[b].pos - p);
        e += pm.energy(t, b, dot(d, d));
      }
    return e;
  }
  // the template's atoms (indices tpl) at positions P with the host
  double energy(const std::vector<uint32_t>& tpl, const std::vector<Vec3>& P) const {
    double e = 0;
    for (size_t k = 0; k < tpl.size(); ++k) {
      e += atom(tpl[k], P[k]);
      if (!std::isfinite(e)) return e;
    }
    return e;
  }
  double between(const std::vector<uint32_t>& tp, const std::vector<Vec3>& P, const std::vector<uint32_t>& tq, const std::vector<Vec3>& Q) const {
    double e = 0;
    for (size_t k = 0; k < tp.size(); ++k)
      for (size_t q = 0; q < tq.size(); ++q) e += pair(tp[k], tq[q], Q[q] - P[k]);
    return e;
  }
};

// A sorbate species: its template atoms and, when flexible, its growth plan. Local atom indices 0 … m−1 follow the
// template's order; the root piece is placed first, then one piece per rotatable bond (parents before children).
struct Species {
  uint32_t first = 0;                 // the template's first force-field atom index
  std::vector<uint32_t> tpl;          // force-field atom indices
  std::vector<Vec3> shape;            // the template about its centroid
  std::vector<int> root;              // local atoms of the root piece
  struct Step { int near = 0, far = 0; std::vector<int> atoms, subtree; TorsionRefs refs; };
  std::vector<Step> steps;            // one per rotatable bond, parents first
  std::vector<std::vector<int>> below; // per step: the steps in its subtree (itself first), in growth order
  struct IntraPair { int a, b; bool p14; };
  std::vector<IntraPair> pairs;       // pairs across rotatable bonds that interact (not 1-2 or 1-3)
  std::vector<char> rel;              // m × m: 0 interacting, 1 excluded, 2 1-4
  bool flexible() const { return !steps.empty(); }
};

Species make_species(const System& s, const ForceField& ff, size_t a, size_t b, bool flexible) {
  Species S;
  S.first = uint32_t(a);
  const size_t m = b - a;
  Vec3 c0{0, 0, 0};
  for (size_t i = a; i < b; ++i) S.tpl.push_back(uint32_t(i)), c0 = c0 + s.atoms[i].pos;
  c0 = c0 * (1.0 / double(m));
  for (size_t i = a; i < b; ++i) S.shape.push_back(s.atoms[i].pos - c0);
  S.root.resize(m);
  for (size_t i = 0; i < m; ++i) S.root[i] = int(i);
  if (!flexible || m < 4) return S;
  // the template alone, its bonds renumbered
  System sub;
  for (size_t i = a; i < b; ++i) sub.atoms.push_back(s.atoms[i]);
  for (const auto& bd : s.bonds)
    if (bd.i >= a && bd.i < b && bd.j >= a && bd.j < b) {
      Bond nb = bd;
      nb.i = uint32_t(bd.i - a), nb.j = uint32_t(bd.j - a);
      sub.bonds.push_back(nb);
    }
  const auto nbr = sub.neighbours();
  const Perception P = perceive(sub);
  auto el = [&](uint32_t x) { return sub.atoms[x].element; };
  auto carbonyl = [&](uint32_t c) {
    if (el(c) != 6) return false;
    for (uint32_t w : P.nb[c])
      if ((el(w) == 8 || el(w) == 16) && P.bond_order(c, w) == 2) return true;
    return false;
  };
  auto rotatable = [&](uint32_t x, uint32_t y) {   // as CBMC regrowth: single, acyclic, not amide or ester, both ends inner
    if (P.bond_order(x, y) != 1) return false;
    for (size_t k = 0; k < P.nb[x].size(); ++k)
      if (P.nb[x][k] == y && (P.ring_bond[x][k] || P.arom_bond[x][k])) return false;
    if (nbr[x].size() < 2 || nbr[y].size() < 2) return false;
    auto hetero = [&](uint32_t z) { return el(z) == 7 || el(z) == 8 || el(z) == 16; };
    if ((carbonyl(x) && hetero(y)) || (carbonyl(y) && hetero(x))) return false;
    return true;
  };
  // rigid pieces: bonded atoms not across a rotatable bond
  std::vector<int> piece(m, -1);
  int np = 0;
  for (size_t i = 0; i < m; ++i) {
    if (piece[i] >= 0) continue;
    std::vector<uint32_t> st{uint32_t(i)};
    piece[i] = np;
    while (!st.empty()) {
      const uint32_t u = st.back();
      st.pop_back();
      for (uint32_t w : nbr[u])
        if (piece[w] < 0 && !rotatable(u, w)) piece[w] = np, st.push_back(w);
    }
    ++np;
  }
  if (np < 2) return S;
  std::vector<std::vector<int>> members(static_cast<size_t>(np));
  for (size_t i = 0; i < m; ++i) members[size_t(piece[i])].push_back(int(i));
  // the root: the largest piece (the first of equals)
  int rp = 0;
  for (int p = 1; p < np; ++p) if (members[size_t(p)].size() > members[size_t(rp)].size()) rp = p;
  S.root = members[size_t(rp)];
  // breadth first over the pieces through the rotatable bonds
  std::vector<char> done(size_t(np), 0);
  done[size_t(rp)] = 1;
  std::vector<int> queue{rp};
  for (size_t qi = 0; qi < queue.size(); ++qi)
    for (int u : members[size_t(queue[qi])])
      for (uint32_t w : nbr[size_t(u)])
        if (!done[size_t(piece[w])]) {
          done[size_t(piece[w])] = 1;
          queue.push_back(piece[w]);
          Species::Step st;
          st.near = u, st.far = int(w), st.atoms = members[size_t(piece[w])];
          S.steps.push_back(std::move(st));
        }
  // each step's subtree: the atoms on far's side
  for (auto& st : S.steps) {
    std::vector<char> seen(m, 0);
    seen[size_t(st.near)] = 1, seen[size_t(st.far)] = 1;
    std::vector<int> stack{st.far};
    st.subtree.push_back(st.far);
    while (!stack.empty()) {
      const int u = stack.back();
      stack.pop_back();
      for (uint32_t w : nbr[size_t(u)])
        if (!seen[w]) seen[w] = 1, st.subtree.push_back(int(w)), stack.push_back(int(w));
    }
  }
  S.below.resize(S.steps.size());
  for (size_t k = 0; k < S.steps.size(); ++k) {
    std::vector<char> in(m, 0);
    for (int x : S.steps[k].subtree) in[size_t(x)] = 1;
    for (size_t j = k; j < S.steps.size(); ++j)
      if (in[size_t(S.steps[j].far)]) S.below[k].push_back(int(j));
  }
  // torsion terms about each rotatable bond
  std::map<std::pair<uint32_t, uint32_t>, size_t> step_of;
  for (size_t k = 0; k < S.steps.size(); ++k) {
    const uint32_t x = uint32_t(a) + uint32_t(S.steps[k].near), y = uint32_t(a) + uint32_t(S.steps[k].far);
    step_of[{std::min(x, y), std::max(x, y)}] = k;
  }
  auto key = [](uint32_t x, uint32_t y) { return std::make_pair(std::min(x, y), std::max(x, y)); };
  for (size_t k = 0; k < ff.dihedrals.size(); ++k)
    if (auto it = step_of.find(key(ff.dihedrals[k].j, ff.dihedrals[k].k)); it != step_of.end()) S.steps[it->second].refs.dihedrals.push_back(uint32_t(k));
  for (size_t k = 0; k < ff.dihedrals2.size(); ++k)
    if (auto it = step_of.find(key(ff.dihedrals2[k].j, ff.dihedrals2[k].k)); it != step_of.end()) S.steps[it->second].refs.dihedrals2.push_back(uint32_t(k));
  for (size_t k = 0; k < ff.cbt.size(); ++k)
    if (auto it = step_of.find(key(ff.cbt[k].j, ff.cbt[k].k)); it != step_of.end()) S.steps[it->second].refs.cbt.push_back(uint32_t(k));
  // exclusions: the force field's (1-2, 1-3, 1-4 partners, 1-4 pairs scaled), else from the bonds
  S.rel.assign(m * m, 0);
  if (ff.excluded.size() == s.atoms.size()) {
    for (size_t i = 0; i < m; ++i)
      for (uint32_t g : ff.excluded[a + i])
        if (g >= a && g < b) S.rel[i * m + (g - a)] = 1;
    for (const auto& p : ff.pairs14)
      if (p[0] >= a && p[0] < b && p[1] >= a && p[1] < b) S.rel[(p[0] - a) * m + (p[1] - a)] = S.rel[(p[1] - a) * m + (p[0] - a)] = 2;
  } else {
    for (size_t i = 0; i < m; ++i) {
      std::vector<int> d(m, -1);
      d[i] = 0;
      std::vector<uint32_t> fr{uint32_t(i)};
      for (int depth = 1; depth <= 3; ++depth) {
        std::vector<uint32_t> nx;
        for (uint32_t u : fr)
          for (uint32_t w : nbr[u])
            if (d[w] < 0) d[w] = depth, nx.push_back(w);
        fr = std::move(nx);
      }
      for (size_t j = 0; j < m; ++j)
        if (d[j] > 0) S.rel[i * m + j] = d[j] == 3 ? 2 : 1;
    }
  }
  for (size_t i = 0; i < m; ++i)
    for (size_t j = i + 1; j < m; ++j)
      if (piece[i] != piece[j] && S.rel[i * m + j] != 1) S.pairs.push_back({int(i), int(j), S.rel[i * m + j] == 2});
  return S;
}

struct E3 {
  double host = 0, mut = 0, intra = 0;   // with the host, with the other sorbates, within the molecule
  double sum() const { return host + mut + intra; }
  E3& operator+=(const E3& o) { host += o.host, mut += o.mut, intra += o.intra; return *this; }
  E3& operator-=(const E3& o) { host -= o.host, mut -= o.mut, intra -= o.intra; return *this; }
};

}  // namespace

SorptionReport sorption(const System& s, const ForceField& ff, const SorptionOptions& o) {
  const auto t0 = std::chrono::steady_clock::now();
  SorptionReport rep;
  const size_t n = s.atoms.size();
  if (!s.cell.valid()) throw std::invalid_argument("sorption: the host needs a periodic cell");
  if (o.template_first_atom < 0 || size_t(o.template_first_atom) >= n) throw std::invalid_argument("sorption: no sorbate template (its atoms come last)");
  if (ff.type_index.size() != n || ff.charge.size() != n) throw std::invalid_argument("sorption: the force field was assigned to another structure");
  if (!ff.pair_func.empty() || ff.hbond.on() || ff.sw.on) throw std::invalid_argument("sorption: the energies need Lennard-Jones pairs");
  if (o.temperature <= 0) throw std::invalid_argument("sorption: T > 0");
  if (o.flexible && o.trials < 1) throw std::invalid_argument("sorption: at least one trial per growth step");
  const size_t nh = size_t(o.template_first_atom);
  if (!(o.cutoff > 0)) throw std::invalid_argument("sorption: the cut-off must be positive");
  // the cut-off holds as asked: past half the cell's narrowest width every periodic image within it counts (see Host)
  const PairModel pm(ff, o.cutoff, o.coulomb, o.dsf_alpha);
  const Host H(s, pm, nh);
  if (H.images) {
    char b[320];
    std::snprintf(b, sizeof b, "the %.1f Å cut-off is wider than half the cell (%.2f Å across): every periodic image within it counts (%zu host images)",
                  o.cutoff, H.wmin, H.gpos.size());
    rep.notes.push_back(b);
    if (!o.pressures_kpa.empty() || !o.fixed_loading.empty()) {
      int sc[3];
      for (int d = 0; d < 3; ++d) sc[d] = std::max(1, int(std::ceil(2 * o.cutoff / H.w[d] - 1e-9)));
      std::snprintf(b, sizeof b, "each sorbed molecule also meets its own periodic images, so its neighbours repeat with the cell; for loadings free of that "
                    "make a %dx%dx%d supercell first (the Henry constant, one molecule alone, is unaffected)", sc[0], sc[1], sc[2]);
      rep.notes.push_back(b);
    }
  }
  // a molecule's own images: part of the periodic system's energy once sorbed (GCMC, fixed loading), but not of the
  // Henry constant, which is one molecule alone in the infinite host
  bool count_own = false;
  // the species: each template's atoms, its shape about its centre and, when flexible, its growth plan
  std::vector<int> starts = o.species_first_atom.empty() ? std::vector<int>{o.template_first_atom} : o.species_first_atom;
  for (size_t k = 0; k < starts.size(); ++k)
    if (starts[k] < int(nh) || size_t(starts[k]) >= n || (k > 0 && starts[k] <= starts[k - 1]))
      throw std::invalid_argument("sorption: the species' templates must follow the host in order");
  const size_t ns = starts.size();
  std::vector<Species> sp_(ns);
  std::vector<std::vector<uint32_t>> tpls(ns);
  std::vector<std::vector<Vec3>> shapes(ns);
  bool any_flex = false;
  for (size_t k = 0; k < ns; ++k) {
    const size_t a = size_t(starts[k]), b = k + 1 < ns ? size_t(starts[k + 1]) : n;
    sp_[k] = make_species(s, ff, a, b, o.flexible);
    tpls[k] = sp_[k].tpl, shapes[k] = sp_[k].shape;
    any_flex = any_flex || sp_[k].flexible();
  }
  std::vector<double> y(ns, 1.0 / double(ns));
  if (!o.mole_fractions.empty()) {
    if (o.mole_fractions.size() != ns) throw std::invalid_argument("sorption: one mole fraction per species");
    double t = 0;
    for (double v : o.mole_fractions) { if (v < 0) throw std::invalid_argument("sorption: mole fractions ≥ 0"); t += v; }
    if (t <= 0) throw std::invalid_argument("sorption: the mole fractions add up to zero");
    for (size_t k = 0; k < ns; ++k) y[k] = o.mole_fractions[k] / t;
  }
  if (!o.fixed_loading.empty() && o.fixed_loading.size() != ns) throw std::invalid_argument("sorption: one fixed loading per species");
  const auto& tpl = tpls[0];
  const double beta = 1 / (kB * o.temperature);
  const double V = s.cell.volume();
  rep.volume = V;
  for (size_t i = 0; i < nh; ++i) rep.host_mass += i < ff.mass.size() ? ff.mass[i] : s.mass_of(s.atoms[i]);
  std::mt19937_64 rng(o.seed);
  caps::UniformReal<double> U(0.0, 1.0);
  auto random_quat = [&] {   // uniform on SO(3) (Shoemake 1992)
    const double u1 = U(rng), u2 = 2 * kPi * U(rng), u3 = 2 * kPi * U(rng);
    return Quat{std::sqrt(u1) * std::cos(u3), std::sqrt(1 - u1) * std::sin(u2), std::sqrt(1 - u1) * std::cos(u2), std::sqrt(u1) * std::sin(u3)};
  };
  auto place = [&](size_t sp, const Vec3& c, const Quat& q) {
    std::vector<Vec3> P(shapes[sp].size());
    for (size_t k = 0; k < P.size(); ++k) P[k] = c + rotate(q, shapes[sp][k]);
    return P;
  };
  auto random_point = [&] { return s.cell.to_cartesian(Vec3{U(rng), U(rng), U(rng)}); };
  bool stopped = false;
  const Cell no_cell;

  // ---- the sorbed molecules (GCMC and fixed loading)
  struct Mol { size_t sp; std::vector<Vec3> P; };
  std::vector<Mol> mols;
  std::vector<std::vector<size_t>> of(ns);   // the molecules of each species (indices into mols)
  auto remove = [&](size_t m) {   // mols[m] out: the last takes its place
    const size_t last = mols.size() - 1;
    auto& av = of[mols[m].sp];
    av.erase(std::find(av.begin(), av.end(), m));
    if (m != last) {
      auto& bv = of[mols[last].sp];
      *std::find(bv.begin(), bv.end(), last) = m;
      mols[m] = std::move(mols[last]);
    }
    mols.pop_back();
  };
  // one atom of species sp (local index a) at p with the other sorbed molecules (skip: its own)
  auto mut_atom = [&](size_t sp, int a, const Vec3& p, size_t skip) {
    double e = 0;
    for (size_t j = 0; j < mols.size(); ++j) {
      if (j == skip) continue;
      const auto& tq = tpls[mols[j].sp];
      for (size_t q = 0; q < tq.size(); ++q) e += H.pair(tpls[sp][size_t(a)], tq[q], mols[j].P[q] - p);
    }
    return e;
  };
  // the energy within a flexible molecule that its conformation changes: pairs across rotatable bonds, torsions about them
  auto intra_full = [&](size_t sp, const std::vector<Vec3>& X) {
    const Species& S = sp_[sp];
    double e = 0;
    for (const auto& p : S.pairs) {
      const Vec3 d = X[size_t(p.b)] - X[size_t(p.a)];
      e += pm.energy(S.tpl[size_t(p.a)], S.tpl[size_t(p.b)], dot(d, d), p.p14);
    }
    for (const auto& st : S.steps)
      if (!st.refs.dihedrals.empty() || !st.refs.dihedrals2.empty() || !st.refs.cbt.empty())
        e += torsion_energy(ff, st.refs, [&](uint32_t g) { return X[g - S.first]; }, no_cell);
    return e;
  };
  // a whole molecule with the host and the others (and within itself when flexible)
  auto mol_energy = [&](size_t sp, const std::vector<Vec3>& X, size_t skip, bool with_intra) {
    E3 e;
    e.host = H.energy(tpls[sp], X);
    if (!std::isfinite(e.host)) return e;
    for (size_t j = 0; j < mols.size(); ++j) if (j != skip) e.mut += H.between(tpls[sp], X, tpls[mols[j].sp], mols[j].P);
    if (count_own) e.mut += H.own(tpls[sp], X);
    if (with_intra && sp_[sp].flexible()) e.intra = intra_full(sp, X);
    return e;
  };

  // ---- configurational-bias growth. New: X is grown (from the template's shape when root) and changed; retrace: X is
  // the molecule as it is (its own placement the first trial of each step) and stays. Returns ln W, W = Π (1/k) Σ e^(−βu),
  // and the energy of the placement kept (new) or of X (retrace). external = false: no host, no other molecules.
  const int kt = std::max(1, o.trials);
  struct Grown { double lnw = -kInf; E3 e; };
  std::vector<char> placed;
  std::vector<Vec3> Z;
  // step energy of the local atoms idx at positions Z (the molecule's others as placed)
  auto step_u = [&](const Species& S, size_t sp, const std::vector<int>& idx, const TorsionRefs* refs, size_t skip, bool external) {
    E3 e;
    const size_t m = S.tpl.size();
    for (int a : idx) {
      if (external) {
        e.host += H.atom(S.tpl[size_t(a)], Z[size_t(a)]);
        if (!std::isfinite(e.host)) return e;
        e.mut += mut_atom(sp, a, Z[size_t(a)], skip);
        if (count_own && H.images) {   // the step's share of the molecule's own images (as Host::own, each pair once)
          e.mut += 0.5 * H.lattice_sum(S.tpl[size_t(a)], S.tpl[size_t(a)], Vec3{0, 0, 0}, true);
          for (size_t b2 = 0; b2 < m; ++b2)
            if (placed[b2] || (std::find(idx.begin(), idx.end(), int(b2)) != idx.end() && int(b2) > a))
              e.mut += H.lattice_sum(S.tpl[size_t(a)], S.tpl[b2], Z[b2] - Z[size_t(a)], true);
        }
      }
      for (size_t b2 = 0; b2 < m; ++b2) {
        if (!placed[b2] || S.rel[size_t(a) * m + b2] == 1) continue;
        const Vec3 d = Z[b2] - Z[size_t(a)];
        e.intra += pm.energy(S.tpl[size_t(a)], S.tpl[b2], dot(d, d), S.rel[size_t(a) * m + b2] == 2);
      }
    }
    if (refs && (!refs->dihedrals.empty() || !refs->dihedrals2.empty() || !refs->cbt.empty()))
      e.intra += torsion_energy(ff, *refs, [&](uint32_t g) { return Z[g - S.first]; }, no_cell);
    return e;
  };
  auto grow = [&](size_t sp, std::vector<Vec3>& X, bool root, int from_step, bool retrace, size_t skip, bool external) {
    const Species& S = sp_[sp];
    const size_t m = S.tpl.size();
    Grown g;
    g.lnw = 0;
    placed.assign(m, root ? 0 : 1);
    std::vector<int> order;   // the torsion steps grown
    if (root) for (size_t k = 0; k < S.steps.size(); ++k) order.push_back(int(k));
    else {
      order = S.below[size_t(from_step)];
      for (int a : S.steps[size_t(from_step)].subtree) placed[size_t(a)] = 0;
    }
    std::vector<double> lw(static_cast<size_t>(kt));
    std::vector<E3> ue(static_cast<size_t>(kt));
    const double lnk = std::log(double(kt));
    if (root) {
      if (!retrace) X = S.shape;
      const Vec3 cen = centroid(X);
      std::vector<Vec3> cs(static_cast<size_t>(kt));
      std::vector<Quat> qs(static_cast<size_t>(kt));
      for (int t = 0; t < kt; ++t) {
        Z = X;
        if (!(retrace && t == 0)) {
          cs[size_t(t)] = random_point(), qs[size_t(t)] = random_quat();
          for (int a : S.root) Z[size_t(a)] = cs[size_t(t)] + rotate(qs[size_t(t)], X[size_t(a)] - cen);
        }
        ue[size_t(t)] = external ? step_u(S, sp, S.root, nullptr, skip, true) : E3{};
        const double u = ue[size_t(t)].sum();
        lw[size_t(t)] = std::isfinite(u) ? -beta * u : -kInf;
      }
      const double lse = log_sum_exp(lw);
      if (!std::isfinite(lse)) { g.lnw = -kInf; return g; }
      g.lnw += lse - lnk;
      int pick = 0;
      if (!retrace) {
        const double r = U(rng);
        double acc = 0;
        pick = kt - 1;
        for (int t = 0; t < kt; ++t) { acc += std::exp(lw[size_t(t)] - lse); if (r < acc) { pick = t; break; } }
        while (!std::isfinite(lw[size_t(pick)])) pick = (pick + kt - 1) % kt;
        for (auto& x : X) x = cs[size_t(pick)] + rotate(qs[size_t(pick)], x - cen);
      }
      g.e += ue[size_t(pick)];
      for (int a : S.root) placed[size_t(a)] = 1;
    }
    for (int k : order) {
      const auto& st = S.steps[size_t(k)];
      const Vec3 axis = X[size_t(st.far)] - X[size_t(st.near)];
      const Vec3 u = axis * (1 / norm(axis));
      const Vec3 c = X[size_t(st.far)];
      std::vector<double> ang(size_t(kt), 0.0);
      for (int t = 0; t < kt; ++t) {
        ang[size_t(t)] = retrace && t == 0 ? 0.0 : 2 * kPi * U(rng);
        // the whole side turns: the torsions about this bond reach the next pieces' first atoms, which this angle
        // already fixes (they sit on their own bonds' axes); pair energies count only this piece's atoms
        Z = X;
        for (int a : st.subtree) Z[size_t(a)] = rotate_about(X[size_t(a)], c, u, ang[size_t(t)]);
        ue[size_t(t)] = step_u(S, sp, st.atoms, &st.refs, skip, external);
        const double uu = ue[size_t(t)].sum();
        lw[size_t(t)] = std::isfinite(uu) ? -beta * uu : -kInf;
      }
      const double lse = log_sum_exp(lw);
      if (!std::isfinite(lse)) { g.lnw = -kInf; return g; }
      g.lnw += lse - lnk;
      int pick = 0;
      if (!retrace) {
        const double r = U(rng);
        double acc = 0;
        pick = kt - 1;
        for (int t = 0; t < kt; ++t) { acc += std::exp(lw[size_t(t)] - lse); if (r < acc) { pick = t; break; } }
        while (!std::isfinite(lw[size_t(pick)])) pick = (pick + kt - 1) % kt;
        for (int a : st.subtree) X[size_t(a)] = rotate_about(X[size_t(a)], c, u, ang[size_t(pick)]);   // the side further out rides along
      }
      g.e += ue[size_t(pick)];
      for (int a : st.atoms) placed[size_t(a)] = 1;
    }
    return g;
  };

  // ---- the isolated molecule: ⟨W_IG⟩ and ⟨u_intra⟩ from growths with no host
  std::vector<double> w_ig(ns, 1.0), w_ig_err(ns, 0.0), u_gas(ns, 0.0);
  for (size_t sp = 0; sp < ns && !stopped; ++sp) {
    if (!sp_[sp].flexible()) continue;
    const int nblk = 10, per = std::max(1, o.ig_samples / nblk);
    std::vector<double> blk(nblk, 0.0);
    double sw = 0, swu = 0;
    std::vector<Vec3> X;
    for (int b = 0; b < nblk; ++b) {
      double sum = 0;
      for (int i = 0; i < per; ++i) {
        const Grown g = grow(sp, X, true, 0, false, size_t(-1), false);
        const double w = std::isfinite(g.lnw) ? std::exp(g.lnw) : 0.0;
        sum += w, sw += w, swu += w * g.e.intra;
      }
      blk[size_t(b)] = sum / per;
    }
    double mm = 0, v = 0;
    for (double x : blk) mm += x;
    mm /= nblk;
    for (double x : blk) v += (x - mm) * (x - mm);
    w_ig[sp] = mm, w_ig_err[sp] = std::sqrt(v / (nblk - 1) / nblk), u_gas[sp] = sw > 0 ? swu / sw : 0.0;
    if (!(mm > 0)) throw std::invalid_argument("sorption: the isolated sorbate could not be grown (every conformation overlaps itself)");
    if (o.progress && !o.progress("isolated molecule, species " + std::to_string(sp + 1), 1.0)) stopped = true;
  }

  // ---- Widom, per species
  const double Vm3 = V * 1e-30, mkg = rep.host_mass / kAvogadro * 1e-3;
  for (size_t sp = 0; sp < ns && o.insertions > 0; ++sp) {
    const int nblk = 10;
    const int per = std::max(1, o.insertions / nblk);
    std::vector<double> blk(nblk, 0.0);
    const bool flex = sp_[sp].flexible();
    std::vector<Vec3> X;
    for (int b = 0; b < nblk && !stopped; ++b) {
      double sum = 0;
      for (int i = 0; i < per; ++i) {
        if (flex) {   // a growth by configurational bias: its Rosenbluth weight
          const Grown g = grow(sp, X, true, 0, false, size_t(-1), true);
          sum += std::isfinite(g.lnw) ? std::exp(g.lnw) : 0.0;
        } else {
          const auto P = place(sp, random_point(), random_quat());
          const double e = H.energy(tpls[sp], P);
          sum += std::isfinite(e) ? std::exp(-beta * e) : 0.0;
        }
      }
      blk[size_t(b)] = sum / per;
      if (o.progress && !o.progress(ns > 1 ? "Widom insertion, species " + std::to_string(sp + 1) : std::string("Widom insertion"), (double(sp) + double(b + 1) / nblk) / double(ns)))
        stopped = true;
    }
    double mm = 0, v = 0;
    for (double x : blk) mm += x;
    mm /= nblk;
    for (double x : blk) v += (x - mm) * (x - mm);
    double err = std::sqrt(v / (nblk - 1) / nblk);
    if (flex) {   // ⟨W⟩ / ⟨W_IG⟩, the errors combined
      const double r = mm / w_ig[sp];
      err = r * std::sqrt((mm > 0 ? (err / mm) * (err / mm) : 0.0) + (w_ig_err[sp] / w_ig[sp]) * (w_ig_err[sp] / w_ig[sp]));
      mm = r;
    }
    SpeciesWidom sw;
    sw.fraction = y[sp];
    sw.widom_w = mm;
    sw.widom_error = err;
    sw.mu_ex = mm > 0 ? -std::log(mm) / beta : kInf;
    sw.henry_mol_kg_kpa = mkg > 0 ? mm * Vm3 / (kBoltzmannJ * kAvogadro * o.temperature) / mkg * 1000.0 : 0.0;
    sw.solubility = mm * 273.15 / o.temperature;
    rep.species.push_back(sw);
  }
  if (!rep.species.empty()) {
    const auto& f = rep.species.front();
    rep.widom_w = f.widom_w, rep.widom_error = f.widom_error, rep.mu_ex = f.mu_ex, rep.henry_mol_kg_kpa = f.henry_mol_kg_kpa, rep.solubility = f.solubility;
  } else
    for (size_t sp = 0; sp < ns; ++sp) rep.species.push_back(SpeciesWidom{y[sp]});
  for (size_t sp = 0; sp < ns; ++sp)
    rep.species[sp].rotatable = int(sp_[sp].steps.size()), rep.species[sp].w_ig = w_ig[sp], rep.species[sp].w_ig_error = w_ig_err[sp], rep.species[sp].u_intra_gas = u_gas[sp];

  // ---- the moves that keep N: translation, rotation, regrowth beyond a rotatable bond, twist of one rotatable bond.
  // Each changes cur by the accepted energy change; kind 0 … 3; returns whether accepted (tries counted by the caller).
  const double twist = std::clamp(o.twist_max_deg, 1.0, 180.0) * kPi / 180;
  auto keep_n_move = [&](int kind, size_t m, E3& cur) {
    const size_t sp = mols[m].sp;
    const Species& S = sp_[sp];
    if (kind <= 1) {   // rigid-body moves: the conformation (and so the intramolecular energy) is kept
      std::vector<Vec3> Q = mols[m].P;
      if (kind == 0) {
        const Vec3 d{(2 * U(rng) - 1) * 1.0, (2 * U(rng) - 1) * 1.0, (2 * U(rng) - 1) * 1.0};
        for (auto& q : Q) q = q + d;
      } else {
        const Vec3 c = centroid(Q);
        const Vec3 ax = [&] { const double z = 2 * U(rng) - 1, t = 2 * kPi * U(rng), rr = std::sqrt(std::max(0.0, 1 - z * z)); return Vec3{rr * std::cos(t), rr * std::sin(t), z}; }();
        const double a = (2 * U(rng) - 1) * 30 * kPi / 180;
        const Quat q{std::cos(a / 2), ax[0] * std::sin(a / 2), ax[1] * std::sin(a / 2), ax[2] * std::sin(a / 2)};
        for (auto& x : Q) x = c + rotate(q, x - c);
      }
      const E3 eo = mol_energy(sp, mols[m].P, m, false), en = mol_energy(sp, Q, m, false);
      const double de = en.sum() - eo.sum();
      if (std::isfinite(en.host) && (de <= 0 || U(rng) < std::exp(-beta * de))) {
        mols[m].P = std::move(Q);
        cur += en, cur -= eo;
        return true;
      }
      return false;
    }
    if (!S.flexible()) return false;
    const int k = int(U(rng) * double(S.steps.size())) % int(S.steps.size());
    if (kind == 2) {   // regrowth beyond bond k: retrace the old side, grow a new one, min(1, W_new / W_old)
      std::vector<Vec3> Xo = mols[m].P, Xn = mols[m].P;
      const Grown go = grow(sp, Xo, false, k, true, m, true);
      const Grown gn = grow(sp, Xn, false, k, false, m, true);
      if (!std::isfinite(gn.lnw)) return false;
      if (std::isfinite(go.lnw) && !(std::log(std::max(U(rng), 1e-300)) < gn.lnw - go.lnw)) return false;
      mols[m].P = std::move(Xn);
      cur += gn.e, cur -= go.e;
      return true;
    }
    // twist: bond k turned by up to ±twist, Metropolis
    const auto& st = S.steps[size_t(k)];
    std::vector<Vec3> Q = mols[m].P;
    const Vec3 axis = Q[size_t(st.far)] - Q[size_t(st.near)];
    const Vec3 u = axis * (1 / norm(axis));
    const Vec3 c = Q[size_t(st.far)];
    const double t = (2 * U(rng) - 1) * twist;
    for (int a : st.subtree) Q[size_t(a)] = rotate_about(Q[size_t(a)], c, u, t);
    const E3 eo = mol_energy(sp, mols[m].P, m, true), en = mol_energy(sp, Q, m, true);
    const double de = en.sum() - eo.sum();
    if (std::isfinite(en.host) && std::isfinite(en.intra) && (de <= 0 || U(rng) < std::exp(-beta * de))) {
      mols[m].P = std::move(Q);
      cur += en, cur -= eo;
      return true;
    }
    return false;
  };
  auto map_sample = [&](std::vector<double>& dens, int g) {   // the molecules' centres on the map (every species)
    for (const auto& mol : mols) {
      Vec3 fr = s.cell.to_fractional(centroid(mol.P));
      int ix[3];
      for (int k = 0; k < 3; ++k) ix[k] = std::clamp(int(std::floor((fr[k] - std::floor(fr[k])) * g)), 0, g - 1);
      dens[(size_t(ix[0]) * g + size_t(ix[1])) * g + size_t(ix[2])] += 1;
    }
  };

  // ---- GCMC at each pressure (from here on molecules are sorbed: their own images count)
  count_own = true;
  for (size_t pi = 0; pi < o.pressures_kpa.size() && !stopped; ++pi) {
    const double p = o.pressures_kpa[pi];
    std::vector<double> bfV(ns);   // β f_i V per species, ideal-gas fugacity y_i p
    for (size_t k = 0; k < ns; ++k) bfV[k] = y[k] * p * 1000.0 * V * 1e-30 / (kBoltzmannJ * o.temperature);
    mols.clear();
    for (auto& v : of) v.clear();
    E3 cur;
    long ins_try = 0, ins_ok = 0, del_try = 0, del_ok = 0, reg_try = 0, reg_ok = 0, tw_try = 0, tw_ok = 0;
    const int equil = o.steps / 4;
    // sums for the averages and the fluctuations: N_i, N_i N_j, U, U N_i
    std::vector<double> sN(ns, 0.0), sUN(ns, 0.0), sNN(ns * ns, 0.0);
    double sU = 0;
    long samples = 0;
    const int nblk = 10;
    std::vector<double> blkC(nblk, 0.0), blkN(size_t(nblk) * ns, 0.0);
    const int report_every = std::max(1, o.steps / 50);
    const int g = std::clamp(o.map_grid, 0, 96);
    std::vector<double> dens(size_t(g) * g * g, 0.0);
    long map_samples = 0;
    const int map_every = std::max(1, (o.steps - equil) / 2000);
    // move shares: rigid sorbates as before (¼ each); with flexible ones, regrowth and twists too
    const double c_ins = any_flex ? 0.2 : 0.25, c_del = any_flex ? 0.4 : 0.5, c_tr = any_flex ? 0.55 : 0.75, c_rot = any_flex ? 0.7 : 1.0, c_reg = 0.9;
    std::vector<Vec3> X;
    for (int step = 0; step < o.steps && !stopped; ++step) {
      const double r = U(rng);
      const size_t N = mols.size();
      if (r < c_ins) {   // insert one molecule of a species picked at random
        ++ins_try;
        const size_t sp = std::min(ns - 1, size_t(U(rng) * double(ns)));
        if (sp_[sp].flexible()) {
          const Grown gr = grow(sp, X, true, 0, false, size_t(-1), true);
          if (std::isfinite(gr.lnw) && U(rng) < bfV[sp] / double(of[sp].size() + 1) * std::exp(gr.lnw) / w_ig[sp]) {
            of[sp].push_back(mols.size());
            mols.push_back({sp, X});
            cur += gr.e;
            ++ins_ok;
          }
        } else {
          auto P = place(sp, random_point(), random_quat());
          const E3 e = mol_energy(sp, P, size_t(-1), false);
          const double dU = e.sum();
          if (std::isfinite(dU) && U(rng) < bfV[sp] / double(of[sp].size() + 1) * std::exp(-beta * dU)) {
            of[sp].push_back(mols.size());
            mols.push_back({sp, std::move(P)});
            cur += e;
            ++ins_ok;
          }
        }
      } else if (r < c_del) {   // delete one of a species picked at random
        const size_t sp = std::min(ns - 1, size_t(U(rng) * double(ns)));
        const size_t Ns = of[sp].size();
        if (Ns > 0) {
          ++del_try;
          const size_t m = of[sp][size_t(U(rng) * double(Ns)) % Ns];
          if (sp_[sp].flexible()) {
            std::vector<Vec3> Xo = mols[m].P;
            const Grown gr = grow(sp, Xo, true, 0, true, m, true);
            const double acc = bfV[sp] <= 0 ? kInf : double(Ns) / bfV[sp] * w_ig[sp] * std::exp(-gr.lnw);
            if (U(rng) < acc) { cur -= gr.e; remove(m); ++del_ok; }
          } else {
            const E3 e = mol_energy(sp, mols[m].P, m, false);
            if (bfV[sp] <= 0 || U(rng) < double(Ns) / bfV[sp] * std::exp(beta * e.sum())) { cur -= e; remove(m); ++del_ok; }
          }
        }
      } else if (N > 0) {   // keep N: translate, rotate, regrow, twist
        const size_t m = size_t(U(rng) * double(N)) % N;
        int kind = r < c_tr ? 0 : r < c_rot ? 1 : r < c_reg ? 2 : 3;
        if (kind >= 2 && !sp_[mols[m].sp].flexible()) kind = kind == 2 ? 0 : 1;   // a rigid molecule: a rigid-body move instead
        const bool ok = keep_n_move(kind, m, cur);
        if (kind == 2) ++reg_try, reg_ok += ok;
        if (kind == 3) ++tw_try, tw_ok += ok;
      }
      if (step >= equil) {
        const double Ucur = cur.sum();
        sU += Ucur, ++samples;
        const int b = std::min(nblk - 1, int(double(step - equil) / double(o.steps - equil) * nblk));
        blkC[size_t(b)] += 1;
        for (size_t i = 0; i < ns; ++i) {
          const double Ni = double(of[i].size());
          sN[i] += Ni, sUN[i] += Ucur * Ni, blkN[size_t(b) * ns + i] += Ni;
          for (size_t j = 0; j < ns; ++j) sNN[i * ns + j] += Ni * double(of[j].size());
        }
        if (g > 0 && (step - equil) % map_every == 0) { map_sample(dens, g); ++map_samples; }
      }
      if (o.progress && (step + 1) % report_every == 0) {
        char st[80];
        std::snprintf(st, sizeof st, "GCMC at %g kPa · %zu molecules", p, mols.size());
        if (!o.progress(st, (pi + double(step + 1) / o.steps) / o.pressures_kpa.size())) stopped = true;
      }
    }
    IsothermPoint pt;
    pt.pressure_kpa = p;
    if (samples > 0) {
      const double inv = 1.0 / double(samples), mU = sU * inv;
      std::vector<double> mN(ns), cUN(ns), cNN(ns * ns);
      for (size_t i = 0; i < ns; ++i) mN[i] = sN[i] * inv;
      for (size_t i = 0; i < ns; ++i) {
        cUN[i] = sUN[i] * inv - mU * mN[i];
        for (size_t j = 0; j < ns; ++j) cNN[i * ns + j] = sNN[i * ns + j] * inv - mN[i] * mN[j];
      }
      pt.species_heat.assign(ns, 0.0);
      {   // q_i = kT + ⟨u_intra,i⟩_gas − Σ_j cov(U, N_j) [cov(N, N)⁻¹]_ji: solve cov(N, N) z = cov(U, N) (symmetric), Gauss with
          // pivoting; species never adsorbed left out
        std::vector<size_t> act;
        for (size_t i = 0; i < ns; ++i) if (cNN[i * ns + i] > 1e-9) act.push_back(i);
        const size_t m = act.size();
        std::vector<double> A(m * m), z(m);
        for (size_t a = 0; a < m; ++a) {
          z[a] = cUN[act[a]];
          for (size_t b = 0; b < m; ++b) A[a * m + b] = cNN[act[a] * ns + act[b]];
        }
        bool ok = true;
        for (size_t c = 0; c < m && ok; ++c) {
          size_t piv = c;
          for (size_t r2 = c + 1; r2 < m; ++r2) if (std::fabs(A[r2 * m + c]) > std::fabs(A[piv * m + c])) piv = r2;
          if (std::fabs(A[piv * m + c]) < 1e-12) { ok = false; break; }
          for (size_t k = 0; k < m; ++k) std::swap(A[c * m + k], A[piv * m + k]);
          std::swap(z[c], z[piv]);
          for (size_t r2 = 0; r2 < m; ++r2)
            if (r2 != c) {
              const double f = A[r2 * m + c] / A[c * m + c];
              for (size_t k = 0; k < m; ++k) A[r2 * m + k] -= f * A[c * m + k];
              z[r2] -= f * z[c];
            }
        }
        if (ok)
          for (size_t a = 0; a < m; ++a) pt.species_heat[act[a]] = 1 / beta + u_gas[act[a]] - z[a] / A[a * m + a];
      }
      double tot = 0;
      for (double v : mN) tot += v;
      pt.loading = tot;
      pt.species_loading = mN;
      pt.species_error.assign(ns, 0.0);
      // the error of the total and of each species: ten blocks
      auto blk_err = [&](auto value) {
        double bm = 0, bv = 0;
        int nbk = 0;
        for (int b = 0; b < nblk; ++b) if (blkC[size_t(b)] > 0) bm += value(b), ++nbk;
        bm /= std::max(1, nbk);
        for (int b = 0; b < nblk; ++b) if (blkC[size_t(b)] > 0) { const double x = value(b) - bm; bv += x * x; }
        return nbk > 1 ? std::sqrt(bv / (nbk - 1) / nbk) : 0.0;
      };
      pt.loading_error = blk_err([&](int b) { double t = 0; for (size_t i = 0; i < ns; ++i) t += blkN[size_t(b) * ns + i]; return t / blkC[size_t(b)]; });
      for (size_t i = 0; i < ns; ++i) pt.species_error[i] = blk_err([&](int b) { return blkN[size_t(b) * ns + i] / blkC[size_t(b)]; });
      if (ns == 1) {   // one species: the single-component formula
        pt.heat = pt.species_heat[0];
      } else {         // the adsorbed phase's mean: Σ x_i q_i
        pt.heat = 0;
        if (tot > 0) for (size_t i = 0; i < ns; ++i) pt.heat += mN[i] / tot * pt.species_heat[i];
      }
      pt.mol_per_kg = mkg > 0 ? tot / kAvogadro / mkg : 0;
      pt.species_mol_per_kg.resize(ns);
      for (size_t i = 0; i < ns; ++i) pt.species_mol_per_kg[i] = mkg > 0 ? mN[i] / kAvogadro / mkg : 0;
      pt.cm3stp_per_cm3 = tot / (V * 1e-24) * 22413.969 / kAvogadro;   // molecules per cm³ → cm³(STP) per cm³
      pt.selectivity.assign(ns, 0.0);
      for (size_t i = 0; i < ns; ++i)   // (x_i / x_0) / (y_i / y_0)
        if (mN[0] > 0 && y[i] > 0) pt.selectivity[i] = (mN[i] / mN[0]) / (y[i] / y[0]);
    }
    pt.acceptance_insert = ins_try ? double(ins_ok) / ins_try : 0;
    pt.acceptance_delete = del_try ? double(del_ok) / del_try : 0;
    pt.acceptance_regrow = reg_try ? double(reg_ok) / reg_try : 0;
    pt.acceptance_twist = tw_try ? double(tw_ok) / tw_try : 0;
    for (const auto& mol : mols) pt.molecules.push_back(mol.P);
    if (g > 0 && map_samples > 0) {   // per Å³: counts over the samples and the voxel's volume
      const double vox = V / (double(g) * g * g);
      pt.grid = g;
      pt.density.resize(dens.size());
      for (size_t k = 0; k < dens.size(); ++k) pt.density[k] = float(dens[k] / (double(map_samples) * vox));
    }
    rep.isotherm.push_back(std::move(pt));
  }

  // ---- fixed loading: canonical Monte Carlo at a set number of molecules
  if (!o.fixed_loading.empty() && !stopped) {
    auto& F = rep.fixed;
    rep.fixed_ran = true;
    mols.clear();
    for (auto& v : of) v.clear();
    E3 cur;
    F.molecules.assign(ns, 0);
    for (size_t k = 0; k < ns; ++k) F.requested += std::max(0, o.fixed_loading[k]);
    // grown in one by one (configurational bias, k trials even for rigid molecules), each species in turn
    std::vector<Vec3> X;
    int attempts = 0;
    for (size_t k = 0; k < ns && !stopped; ++k)
      while (F.molecules[k] < o.fixed_loading[k] && attempts < o.loading_attempts) {
        ++attempts;
        const Grown gr = grow(k, X, true, 0, false, size_t(-1), true);
        if (!std::isfinite(gr.lnw)) continue;
        of[k].push_back(mols.size());
        mols.push_back({k, X});
        cur += gr.e;
        ++F.molecules[k], ++F.loaded;
        if (o.progress && F.loaded % 10 == 0 && !o.progress("loading " + std::to_string(F.loaded) + " of " + std::to_string(F.requested), 0.0)) stopped = true;
      }
    if (F.loaded < F.requested) {
      char b[200];
      std::snprintf(b, sizeof b, "fixed loading: %d of %d molecules fitted in %d attempts", F.loaded, F.requested, attempts);
      rep.notes.push_back(b);
    }
    const int equil = o.steps / 4, nblk = 10;
    const int sample_every = std::max(1, (o.steps - equil) / 2000), report_every = std::max(1, o.steps / 50);
    long tries[4] = {0, 0, 0, 0}, oks[4] = {0, 0, 0, 0};
    double sE = 0, sH = 0, sM = 0, sI = 0;
    long samples = 0;
    std::vector<double> blk(nblk, 0.0), blkC(nblk, 0.0), binding;
    const int g = std::clamp(o.map_grid, 0, 96);
    std::vector<double> dens(size_t(g) * g * g, 0.0);
    long map_samples = 0;
    F.lowest_energy = kInf;
    for (int step = 0; step < o.steps && !stopped && !mols.empty(); ++step) {
      const double r = U(rng);
      const size_t m = size_t(U(rng) * double(mols.size())) % mols.size();
      int kind = any_flex ? (r < 0.3 ? 0 : r < 0.6 ? 1 : r < 0.85 ? 2 : 3) : (r < 0.5 ? 0 : 1);
      if (kind >= 2 && !sp_[mols[m].sp].flexible()) kind = kind == 2 ? 0 : 1;
      ++tries[kind];
      oks[kind] += keep_n_move(kind, m, cur);
      if (step >= equil) {
        const double e = cur.sum();
        sE += e, sH += cur.host, sM += cur.mut, sI += cur.intra, ++samples;
        const int b = std::min(nblk - 1, int(double(step - equil) / double(o.steps - equil) * nblk));
        blk[size_t(b)] += e, blkC[size_t(b)] += 1;
        if (e < F.lowest_energy) {
          F.lowest_energy = e;
          F.lowest.clear();
          for (const auto& mol : mols) F.lowest.push_back(mol.P);
        }
        if ((step - equil) % sample_every == 0) {
          for (size_t j = 0; j < mols.size(); ++j) {   // each molecule's binding energy with the host and the others
            const E3 eb = mol_energy(mols[j].sp, mols[j].P, j, false);
            binding.push_back(eb.host + eb.mut);
          }
          if (g > 0) { map_sample(dens, g); ++map_samples; }
        }
      }
      if (o.progress && (step + 1) % report_every == 0) {
        char st[80];
        std::snprintf(st, sizeof st, "fixed loading · %d molecules", F.loaded);
        if (!o.progress(st, double(step + 1) / o.steps)) stopped = true;
      }
    }
    if (samples > 0) {
      const double inv = 1.0 / double(samples);
      F.energy = sE * inv, F.host_sorbate = sH * inv, F.sorbate_sorbate = sM * inv, F.intramolecular = sI * inv;
      double bm = 0, bv = 0;
      int nbk = 0;
      for (int b = 0; b < nblk; ++b) if (blkC[size_t(b)] > 0) bm += blk[size_t(b)] / blkC[size_t(b)], ++nbk;
      bm /= std::max(1, nbk);
      for (int b = 0; b < nblk; ++b) if (blkC[size_t(b)] > 0) { const double x = blk[size_t(b)] / blkC[size_t(b)] - bm; bv += x * x; }
      F.energy_error = nbk > 1 ? std::sqrt(bv / (nbk - 1) / nbk) : 0.0;
      double ref = 0;   // the molecules' own intramolecular energy in the gas
      for (size_t k = 0; k < ns; ++k) ref += F.molecules[k] * u_gas[k];
      if (F.loaded > 0) F.adsorption_energy = (F.energy - ref) / F.loaded;
    }
    if (!binding.empty()) {
      double lo = *std::min_element(binding.begin(), binding.end()), hi = *std::max_element(binding.begin(), binding.end());
      if (hi - lo < 1e-9) hi = lo + 1;
      const int nb = std::clamp(o.energy_bins, 5, 400);
      F.hist_edges.resize(size_t(nb) + 1), F.hist_counts.assign(size_t(nb), 0.0);
      for (int k = 0; k <= nb; ++k) F.hist_edges[size_t(k)] = lo + (hi - lo) * k / nb;
      for (double e : binding) F.hist_counts[size_t(std::min(nb - 1, int((e - lo) / (hi - lo) * nb)))] += 1;
    }
    if (!std::isfinite(F.lowest_energy)) F.lowest_energy = 0;
    for (const auto& mol : mols) F.molecules_end.push_back(mol.P);
    F.acceptance_translate = tries[0] ? double(oks[0]) / tries[0] : 0;
    F.acceptance_rotate = tries[1] ? double(oks[1]) / tries[1] : 0;
    F.acceptance_regrow = tries[2] ? double(oks[2]) / tries[2] : 0;
    F.acceptance_twist = tries[3] ? double(oks[3]) / tries[3] : 0;
    if (g > 0 && map_samples > 0) {
      const double vox = V / (double(g) * g * g);
      F.grid = g;
      F.density.resize(dens.size());
      for (size_t k = 0; k < dens.size(); ++k) F.density[k] = float(dens[k] / (double(map_samples) * vox));
    }
  }

  if (stopped) rep.notes.push_back("stopped before the end");
  {
    char b[300];
    int rot = 0;
    for (const auto& S : sp_) rot += int(S.steps.size());
    const char* kind = rot > 0 ? "flexible" : "rigid";
    if (ns == 1)
      std::snprintf(b, sizeof b, "host %zu atoms (%.0f g/mol per cell), cell %.0f Å³ · sorbate %zu atoms, %s · %.0f K", nh, rep.host_mass, V, tpl.size(), kind, o.temperature);
    else
      std::snprintf(b, sizeof b, "host %zu atoms (%.0f g/mol per cell), cell %.0f Å³ · %zu sorbates, %s, ideal gas mixture · %.0f K", nh, rep.host_mass, V, ns, kind, o.temperature);
    rep.notes.insert(rep.notes.begin(), b);
    if (rot > 0) {
      std::snprintf(b, sizeof b, "flexible sorbates: %d rotatable bond%s grown by configurational bias, %d trials per step (Siepmann & Frenkel 1992); bond lengths, angles and rings as the template",
                    rot, rot == 1 ? "" : "s", kt);
      rep.notes.push_back(b);
    } else if (o.flexible)
      rep.notes.push_back("no rotatable bond in the sorbate: treated as rigid");
    rep.notes.push_back(std::string("energies: the force field's van der Waals") + (pm.coulomb() ? " and damped shifted force electrostatics" : "") +
                        "; the host is held fixed; ideal-gas reservoir (fugacity = " + std::string(ns > 1 ? "y_i × pressure)" : "pressure)"));
  }
  rep.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  return rep;
}

}  // namespace caps
