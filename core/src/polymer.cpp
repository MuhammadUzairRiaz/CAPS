// CAPS polymer builder: chains of any repeat unit grown into a periodic cell (see polymer.hpp).
#include "caps/polymer.hpp"
#include "caps/uff.hpp"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <array>
#include <cmath>
#include <functional>
#include <map>
#include <optional>
#include <numeric>
#include <random>
#include <set>
#include <sstream>

#include "caps/elements.hpp"

namespace caps {
namespace {

constexpr double kPi = 3.14159265358979323846;
Vec3 unitv(const Vec3& v) { return v * (1.0 / norm(v)); }

// NeRF: position of d given a-b-c, |cd| = bond, angle bcd = ang (rad), torsion abcd = tor (rad).
Vec3 place(const Vec3& a, const Vec3& b, const Vec3& c, double bond, double ang, double tor) {
  const Vec3 bc = unitv(c - b);
  const Vec3 n = unitv(cross(b - a, bc));
  const Vec3 m = cross(n, bc);
  const Vec3 d2{-bond * std::cos(ang), bond * std::sin(ang) * std::cos(tor), bond * std::sin(ang) * std::sin(tor)};
  return c + bc * d2[0] + m * d2[1] + n * d2[2];
}

double angle_at(const Vec3& a, const Vec3& b, const Vec3& c) {
  const Vec3 u = unitv(a - b), v = unitv(c - b);
  return std::acos(std::clamp(dot(u, v), -1.0, 1.0));
}

double dihedral(const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& d) {
  const Vec3 b1 = b - a, b2 = c - b, b3 = d - c;
  const Vec3 n1 = cross(b1, b2), n2 = cross(b2, b3);
  const Vec3 m = cross(n1, unitv(b2));
  return -std::atan2(dot(m, n2), dot(n1, n2));   // the sign convention of place(): place(..., dihedral(...)) rebuilds d
}

// torsion kinds: 0 sp3–sp3 (trans, gauche±), 1 next to an sp2 atom (30° steps), 2 conjugated (0 or 180°),
// 3 ester / amide link (trans)
int torsion_kind(bool sp2a, bool donor_a, bool sp2b, bool donor_b) {
  if (sp2a && sp2b) return 2;
  if ((sp2a && donor_b) || (sp2b && donor_a)) return 3;
  if (sp2a || sp2b) return 1;
  return 0;
}

struct Template {
  std::string name;
  int n = 0;
  std::vector<int> z;
  std::vector<bool> sp2, donor;
  std::vector<std::array<int, 3>> ref;   // parent, grandparent, great-grandparent; < 0: -1 previous tail, -2 its parent, -3 …
  std::vector<double> bond, ang, tor;
  std::vector<int> group;                // sampled torsion group per atom, -1 fixed (template torsion)
  std::vector<double> offset;            // torsion − the group primary's
  std::vector<int> group_kind;           // per group; the link group's kind is set per step
  int link_group = -1;                   // the group about the link bond (head – previous tail)
  int tail = 0;                          // unit index (the head is 0)
  double tail_angle = 1.95;              // angle(next head, tail, tail's parent)
  double head_bond = 1.53, tail_bond = 1.53;   // head / tail to a carbon in the template
  std::vector<std::array<int, 3>> bonds; // a, b, order
  std::vector<int> parent;               // tree parent (unit index, -1 for the head)
  int head_extra_h = 0, tail_extra_h = 0;   // bracket end atoms that need a cap hydrogen in the chain graph
  // where the next unit bonds: fixed by the tail's own substituents when it has a child in the tree (the direction
  // in a frame from the tail to its parent r1 and child r2, r1 may be -1 = the unit's link partner); otherwise sampled
  bool tail_fixed = false;
  int r1 = -1, r2 = -1;
  double tf[3] = {0, 0, 0};
};

// orthonormal frame at t from its neighbours a (e1 toward a) and b
void frame(const Vec3& t, const Vec3& a, const Vec3& b, Vec3& e1, Vec3& e2, Vec3& e3) {
  e1 = unitv(a - t);
  const Vec3 w = b - t;
  e2 = unitv(w - e1 * dot(w, e1));
  e3 = cross(e1, e2);
}

struct Parsed {
  MolGraph g;
  int d1 = -1, d2 = -1, head = -1, tail = -1;
};

Parsed parse_unit(const std::string& smiles) {
  Parsed p;
  p.g = parse_smiles(smiles);
  std::vector<int> dummies;
  for (size_t i = 0; i < p.g.atoms.size(); ++i)
    if (p.g.atoms[i].element == 0) dummies.push_back(int(i));
  if (dummies.size() != 2) throw std::runtime_error("a repeat unit needs exactly two attachment points (*), head first; this one has " + std::to_string(dummies.size()));
  // [*:1] / [*:2] fix the order; otherwise the first * written is the head
  if (p.g.atoms[size_t(dummies[0])].map == 2 && p.g.atoms[size_t(dummies[1])].map == 1) std::swap(dummies[0], dummies[1]);
  p.d1 = dummies[0];
  p.d2 = dummies[1];
  auto nb = [&](int d) {
    std::vector<int> v;
    for (const auto& b : p.g.bonds) {
      if (b.a == d) v.push_back(b.b);
      if (b.b == d) v.push_back(b.a);
    }
    return v;
  };
  const auto n1 = nb(p.d1), n2 = nb(p.d2);
  if (n1.size() != 1 || n2.size() != 1) throw std::runtime_error("each attachment point must be bonded to one atom");
  if (p.g.atoms[size_t(n1[0])].element == 0 || p.g.atoms[size_t(n2[0])].element == 0) throw std::runtime_error("the attachment points are bonded to each other");
  for (const auto& b : p.g.bonds)
    if ((b.a == p.d1 || b.b == p.d1 || b.a == p.d2 || b.b == p.d2) && b.order != 1)
      throw std::runtime_error("the links to the attachment points must be single bonds");
  p.head = n1[0];
  p.tail = n2[0];
  if (p.g.parts != 1) throw std::runtime_error("a repeat unit must be one connected piece");
  return p;
}

// The unit with carbons in place of the attachment points, hydrogens added.
MolGraph capped(const Parsed& p) {
  MolGraph g = p.g;
  for (int d : {p.d1, p.d2}) {
    auto& a = g.atoms[size_t(d)];
    a.element = 6;
    a.bracket = false;
    a.hcount = -1;
    a.aromatic = false;
    a.chiral = 0;
  }
  add_hydrogens(g);
  return g;
}

Template make_template(const std::string& name, const std::string& smiles, const std::string& ff, uint64_t seed, std::vector<std::string>& notes) {
  const Parsed P = parse_unit(smiles);
  MolGraph g = capped(P);
  const auto problems = valence_problems(g);
  if (!problems.empty()) throw std::runtime_error(name + ": " + problems.front());
  EmbedOptions eo;
  eo.seed = seed;
  std::vector<Vec3> pos = embed(g, eo);
  if (!ff.empty()) {
    std::vector<std::string> fn;
    auto f = molecule_forcefield(g, ff, "gasteiger", fn, nullptr, &pos);
    if (!f && !is_uff(ff)) {   // units the chosen force field cannot type (silicones, phosphazenes …): UFF
      std::vector<std::string> un;
      f = molecule_forcefield(g, "uff", "gasteiger", un, nullptr, &pos);
      if (f) fn.clear();
    }
    std::string why;
    if (!f || !minimise_molecule(g, f, 0.05, pos, nullptr, &why))
      notes.push_back(name + ": unit geometry from the embedding" + (fn.empty() ? (why.empty() ? "" : " (" + why + ")") : " (" + fn.front() + ")"));
  }
  const size_t N = g.atoms.size();
  std::vector<std::vector<std::pair<int, int>>> adj(N);
  for (const auto& b : g.bonds) adj[size_t(b.a)].push_back({b.b, b.order}), adj[size_t(b.b)].push_back({b.a, b.order});
  // leave out the stand-in carbons and their hydrogens
  std::vector<bool> keep(N, true);
  for (int d : {P.d1, P.d2}) {
    keep[size_t(d)] = false;
    for (auto [w, o] : adj[size_t(d)])
      if (g.atoms[size_t(w)].element == 1) keep[size_t(w)] = false;
  }
  // spanning tree from the head: heavy atoms before hydrogens
  std::vector<int> order, parent_old(N, -2), newidx(N, -1);
  std::vector<int> q{P.head};
  parent_old[size_t(P.head)] = -1;
  for (size_t h = 0; h < q.size(); ++h) {
    const int u = q[h];
    order.push_back(u);
    std::vector<int> nb;
    for (auto [w, o] : adj[size_t(u)])
      if (keep[size_t(w)] && parent_old[size_t(w)] == -2) nb.push_back(w);
    std::stable_sort(nb.begin(), nb.end(), [&](int a, int b) { return (g.atoms[size_t(a)].element == 1) < (g.atoms[size_t(b)].element == 1); });
    for (int w : nb) parent_old[size_t(w)] = u, q.push_back(w);
  }
  for (size_t k = 0; k < order.size(); ++k) newidx[size_t(order[k])] = int(k);
  Template T;
  T.name = name;
  T.n = int(order.size());
  // stand-ins for the previous unit: the head's carbon and two of its hydrogens
  std::vector<int> dh;
  for (auto [w, o] : adj[size_t(P.d1)])
    if (g.atoms[size_t(w)].element == 1) dh.push_back(w);
  const Vec3 vm1 = pos[size_t(P.d1)], vm2 = pos[size_t(dh.at(0))], vm3 = pos[size_t(dh.at(1))];
  auto X = [&](int local) -> Vec3 {
    if (local >= 0) return pos[size_t(order[size_t(local)])];
    return local == -1 ? vm1 : local == -2 ? vm2 : vm3;
  };
  T.parent.resize(size_t(T.n));
  for (int k = 0; k < T.n; ++k) T.parent[size_t(k)] = parent_old[size_t(order[size_t(k)])] < 0 ? -1 : newidx[size_t(parent_old[size_t(order[size_t(k)])])];
  auto up = [&](int local) { return local >= 0 ? T.parent[size_t(local)] : local - 1; };
  for (int k = 0; k < T.n; ++k) {
    const int old = order[size_t(k)];
    T.z.push_back(g.atoms[size_t(old)].element);
    bool pi = false;
    for (auto [w, o] : adj[size_t(old)]) pi |= o == 2 || o == 4 || o == 3;
    T.sp2.push_back(pi);
    T.donor.push_back(!pi && (g.atoms[size_t(old)].element == 7 || g.atoms[size_t(old)].element == 8 || g.atoms[size_t(old)].element == 16));
    const int p = up(k), gp = up(p), ggp = up(gp);
    T.ref.push_back({p, gp, ggp});
    T.bond.push_back(norm(X(k) - X(p)));
    T.ang.push_back(angle_at(X(k), X(p), X(gp)));
    T.tor.push_back(dihedral(X(k), X(p), X(gp), X(ggp)));
  }
  T.tail = newidx[size_t(P.tail)];
  T.head_bond = norm(pos[size_t(P.head)] - pos[size_t(P.d1)]);
  T.tail_bond = norm(pos[size_t(P.tail)] - pos[size_t(P.d2)]);
  T.tail_angle = angle_at(pos[size_t(P.d2)], pos[size_t(P.tail)], X(T.parent[size_t(T.tail)]));
  for (int k = 0; k < T.n; ++k)
    if (T.parent[size_t(k)] == T.tail) { T.r2 = k; break; }
  if (T.r2 >= 0) {
    T.tail_fixed = true;
    T.r1 = T.parent[size_t(T.tail)];
    Vec3 e1, e2, e3;
    frame(X(T.tail), X(T.r1), X(T.r2), e1, e2, e3);
    const Vec3 d = unitv(pos[size_t(P.d2)] - X(T.tail));
    T.tf[0] = dot(d, e1), T.tf[1] = dot(d, e2), T.tf[2] = dot(d, e3);
  }
  // rings: a bond is in a ring when its ends stay connected without it
  auto in_ring = [&](int a, int b) {
    std::vector<int> seen(N, 0), st{a};
    seen[size_t(a)] = 1;
    while (!st.empty()) {
      const int u = st.back();
      st.pop_back();
      for (auto [w, o] : adj[size_t(u)]) {
        if ((u == a && w == b) || seen[size_t(w)]) continue;
        if (w == b) return true;
        seen[size_t(w)] = 1;
        st.push_back(w);
      }
    }
    return false;
  };
  auto order_of = [&](int a, int b) {
    for (auto [w, o] : adj[size_t(a)])
      if (w == b) return o;
    return 0;
  };
  // torsion groups: atoms turning about the same bond share one sampled value
  std::map<std::pair<int, int>, int> groups;
  T.group.assign(size_t(T.n), -1);
  T.offset.assign(size_t(T.n), 0.0);
  std::vector<int> primary;
  for (int k = 1; k < T.n; ++k) {
    const int p = T.ref[size_t(k)][0], gp = T.ref[size_t(k)][1];
    int kind;
    if (gp < 0) kind = -2;   // the link bond: decided at each step from the previous unit's tail
    else {
      const int op = order[size_t(p)], og = order[size_t(gp)];
      if (order_of(op, og) != 1 || in_ring(op, og)) continue;
      const bool sp = T.sp2[size_t(p)], sg = T.sp2[size_t(gp)];
      if (sp && sg) continue;                                    // conjugated: the template's (planar) torsion
      if ((sp && T.donor[size_t(gp)]) || (sg && T.donor[size_t(p)])) continue;   // ester, amide: as built
      kind = torsion_kind(sp, T.donor[size_t(p)], sg, T.donor[size_t(gp)]);
    }
    auto [it, fresh] = groups.emplace(std::make_pair(p, gp), int(T.group_kind.size()));
    if (fresh) {
      T.group_kind.push_back(kind);
      primary.push_back(k);
      if (kind == -2) T.link_group = it->second;
    }
    T.group[size_t(k)] = it->second;
    T.offset[size_t(k)] = T.tor[size_t(k)] - T.tor[size_t(primary[size_t(it->second)])];
  }
  for (const auto& b : g.bonds)
    if (keep[size_t(b.a)] && keep[size_t(b.b)]) T.bonds.push_back({newidx[size_t(b.a)], newidx[size_t(b.b)], b.order});
  T.head_extra_h = P.g.atoms[size_t(P.head)].bracket ? 1 : 0;
  T.tail_extra_h = P.g.atoms[size_t(P.tail)].bracket ? 1 : 0;
  return T;
}

// ---- the cell: every atom placed so far, all chains, with a periodic cell list
struct Cell3 {
  double L[3] = {0, 0, 0}, cs[3] = {4, 4, 4};
  int nc[3] = {1, 1, 1};
  std::vector<std::vector<int>> bins;
  std::vector<Vec3> x;        // by id
  std::vector<int> z, chain, local;   // chain −1: a fixed substrate atom
  std::vector<char> alive;
  void init(const std::array<double, 3>& edge, double reach) {
    for (int k = 0; k < 3; ++k) {
      L[k] = edge[size_t(k)];
      nc[k] = std::max(1, int(std::floor(L[k] / std::max(reach, 1.0))));
      cs[k] = L[k] / nc[k];
    }
    bins.assign(size_t(nc[0]) * size_t(nc[1]) * size_t(nc[2]), {});
  }
  int wi(int k, double v) const { return std::clamp(int(std::floor((v - L[k] * std::floor(v / L[k])) / cs[k])), 0, nc[k] - 1); }
  int bin_of(const Vec3& p) const { return (wi(0, p[0]) * nc[1] + wi(1, p[1])) * nc[2] + wi(2, p[2]); }
  int add(const Vec3& p, int zz, int ch, int loc) {
    const int id = int(x.size());
    x.push_back(p), z.push_back(zz), chain.push_back(ch), local.push_back(loc), alive.push_back(1);
    bins[size_t(bin_of(p))].push_back(id);
    return id;
  }
  void kill(int id) {
    alive[size_t(id)] = 0;
    auto& b = bins[size_t(bin_of(x[size_t(id)]))];
    b.erase(std::remove(b.begin(), b.end(), id), b.end());
  }
  Vec3 mi(Vec3 d) const {
    for (int k = 0; k < 3; ++k) d[k] -= L[k] * std::round(d[k] / L[k]);
    return d;
  }
  template <class F>
  void near(const Vec3& p, F&& f) const {
    int lo[3], hi[3];
    for (int k = 0; k < 3; ++k) {
      if (nc[k] < 3) lo[k] = 0, hi[k] = nc[k] - 1;   // every bin along a short axis
      else lo[k] = wi(k, p[k]) - 1, hi[k] = wi(k, p[k]) + 1;
    }
    for (int i = lo[0]; i <= hi[0]; ++i)
      for (int j = lo[1]; j <= hi[1]; ++j)
        for (int k = lo[2]; k <= hi[2]; ++k) {
          const int xi = (i + nc[0]) % nc[0], yj = (j + nc[1]) % nc[1], zk = (k + nc[2]) % nc[2];
          for (int id : bins[size_t((xi * nc[1] + yj) * nc[2] + zk)]) f(id);
        }
  }
};

struct ChainState {
  std::vector<int> seq;
  std::vector<char> mirror;
  std::vector<Vec3> pos;      // local (unwrapped); 0..2 are ghosts standing in for the start
  std::vector<int> z, tparent, gid, unit_of;
  std::vector<char> sp2, donor, backbone;
  std::vector<std::vector<int>> adj;
  std::vector<int> unit_start;   // local index of each placed unit's head
  int backtracks = 0, fails = 0, starts = 0;
  bool done = false;
  // branched molecules (star, comb, branched): a chain hung on an atom of another, grown once that chain (and the arm
  // before it on the same molecule) is complete
  int mol = 0;                         // molecule, 0-based
  int parent = -1, anchor_unit = 0, after = -1;
  bool core = false;                   // a star arm: on the parent's head atom, the head cap's valence first
  int anchor = -1;                     // the parent's atom the head bonds to (local), resolved when the arm starts
  Vec3 anchor_dir{0, 0, 0};
  bool resolved = false, started = false, cap_taken = false;   // cap_taken: an arm took this chain's head cap valence
  int drop_h = -1;                     // the parent's hydrogen this arm replaced (−1: the head cap's valence)
  std::set<int> tried;                 // hydrogens given back after the arm found no room there
  int reserve = -1;                    // cell id of a stand-in carbon keeping the head's place until the arm starts
  std::set<int> dropped;               // local hydrogens replaced by arms
  std::map<int, int> pdist;            // cell id → bond distance from the anchor, for the molecule's atoms near it
};

}  // namespace

Sequence sequence_from_string(const std::string& s) {
  if (s == "alternating") return Sequence::Alternating;
  if (s == "block") return Sequence::Block;
  if (s == "random") return Sequence::Random;
  if (s == "gradient") return Sequence::Gradient;
  if (s == "terminal") return Sequence::Terminal;
  if (s == "pattern") return Sequence::Pattern;
  return Sequence::Homopolymer;
}

Architecture architecture_from_string(const std::string& s) {
  if (s.empty() || s == "linear") return Architecture::Linear;
  if (s == "star") return Architecture::Star;
  if (s == "comb") return Architecture::Comb;
  if (s == "branched") return Architecture::Branched;
  throw std::invalid_argument("unknown architecture '" + s + "' (linear, star, comb, branched)");
}

const char* to_string(Architecture a) {
  switch (a) {
    case Architecture::Star: return "star";
    case Architecture::Comb: return "comb";
    case Architecture::Branched: return "branched";
    default: return "linear";
  }
}

const char* to_string(Sequence s) {
  switch (s) {
    case Sequence::Alternating: return "alternating";
    case Sequence::Block: return "block";
    case Sequence::Random: return "random";
    case Sequence::Gradient: return "gradient";
    case Sequence::Terminal: return "terminal";
    case Sequence::Pattern: return "pattern";
    default: return "homopolymer";
  }
}

UnitInfo repeat_unit_info(const std::string& smiles) {
  const Parsed P = parse_unit(smiles);
  MolGraph g = capped(P);
  const auto problems = valence_problems(g);
  if (!problems.empty()) throw std::runtime_error(problems.front());
  UnitInfo u;
  u.head = P.head;
  u.tail = P.tail;
  u.head_element = element(P.g.atoms[size_t(P.head)].element).symbol;
  u.tail_element = element(P.g.atoms[size_t(P.tail)].element).symbol;
  // the unit alone: drop the stand-in carbons and their hydrogens
  std::map<int, int> count;
  std::vector<std::vector<int>> adj(g.atoms.size());
  for (const auto& b : g.bonds) adj[size_t(b.a)].push_back(b.b), adj[size_t(b.b)].push_back(b.a);
  std::vector<bool> keep(g.atoms.size(), true);
  for (int d : {P.d1, P.d2}) {
    keep[size_t(d)] = false;
    for (int w : adj[size_t(d)])
      if (g.atoms[size_t(w)].element == 1) keep[size_t(w)] = false;
  }
  for (size_t i = 0; i < g.atoms.size(); ++i)
    if (keep[i]) ++count[g.atoms[i].element], u.mass += element(g.atoms[i].element).mass, ++u.atoms;
  // stereocentres tacticity acts on: sp3 atoms with four different neighbour kinds, two of them backbone
  for (size_t i = 0; i < P.g.atoms.size(); ++i) {
    if (!keep[i] || g.atoms[i].element != 6) continue;
    std::set<std::string> kinds;
    int heavy = 0;
    for (int w : adj[i]) {
      std::string k = element(g.atoms[size_t(w)].element).symbol;
      if (!keep[size_t(w)] && g.atoms[size_t(w)].element == 6) k = "link";
      int deg = 0;
      for (int x : adj[size_t(w)]) deg += g.atoms[size_t(x)].element != 1;
      kinds.insert(k + std::to_string(deg));
      heavy += g.atoms[size_t(w)].element != 1;
    }
    if (adj[i].size() == 4 && kinds.size() >= 3 && heavy >= 3 && (int(i) == P.head || int(i) == P.tail || true)) {
      // a CH with three different heavy neighbours (backbone in, backbone out, side group)
      int h = 0;
      for (int w : adj[i]) h += g.atoms[size_t(w)].element == 1;
      if (h == 1 && kinds.size() == 4) ++u.stereocentres;
    }
  }
  std::ostringstream f;
  const bool carbon = count.count(6) > 0;
  auto put = [&](const std::string& e, int k) { f << e; if (k > 1) f << k; };
  if (carbon) {
    put("C", count[6]);
    if (count.count(1)) put("H", count[1]);
  }
  std::vector<std::pair<std::string, int>> rest;
  for (const auto& [z, k] : count)
    if (!carbon || (z != 6 && z != 1)) rest.push_back({element(z).symbol, k});
  std::sort(rest.begin(), rest.end());
  for (const auto& [e, k] : rest) put(e, k);
  u.formula = f.str();
  return u;
}

std::vector<int> chain_sequence(const ChainSpec& spec, uint64_t seed) {
  const int n = std::max(1, spec.dp), nu = std::max<int>(1, int(spec.units.size()));
  std::vector<int> s(size_t(n), 0);
  std::mt19937_64 rng(seed * 0x9E3779B97F4A7C15ull + 99);
  switch (spec.sequence) {
    case Sequence::Homopolymer: break;
    case Sequence::Alternating:
      for (int i = 0; i < n; ++i) s[size_t(i)] = i % nu;
      break;
    case Sequence::Block: {
      std::vector<int> b = spec.blocks;
      if (b.empty()) b.assign(size_t(nu), std::max(1, n / nu));
      int i = 0, k = 0;
      while (i < n) {
        for (int j = 0; j < std::max(1, b[size_t(k % int(b.size()))]) && i < n; ++j) s[size_t(i++)] = k % nu;
        ++k;
      }
      break;
    }
    case Sequence::Random: {
      std::vector<double> w = spec.weights;
      w.resize(size_t(nu), w.empty() ? 1.0 : 0.0);
      if (std::accumulate(w.begin(), w.end(), 0.0) <= 0) w.assign(size_t(nu), 1.0);
      std::discrete_distribution<int> D(w.begin(), w.end());
      for (int i = 0; i < n; ++i) s[size_t(i)] = D(rng);
      break;
    }
    case Sequence::Gradient: {
      std::uniform_real_distribution<double> U(0, 1);
      for (int i = 0; i < n; ++i) {
        const double f = n > 1 ? double(i) / (n - 1) : 0.0;   // 0 at the head, 1 at the tail
        s[size_t(i)] = U(rng) < f ? nu - 1 : 0;
      }
      break;
    }
    case Sequence::Terminal: {
      // f1: A in the feed; P(A→A) = r1 f1 / (r1 f1 + f2), P(B→B) = r2 f2 / (r2 f2 + f1); the first unit from the
      // instantaneous copolymer composition F1 (Mayo & Lewis 1944)
      double f1 = spec.weights.size() >= 2 ? spec.weights[0] / std::max(1e-12, spec.weights[0] + spec.weights[1]) : spec.weights.size() == 1 ? spec.weights[0] : 0.5;
      f1 = std::clamp(f1, 0.0, 1.0);
      const double f2 = 1 - f1, r1 = std::max(0.0, spec.r1), r2 = std::max(0.0, spec.r2);
      const double paa = r1 * f1 + f2 > 0 ? r1 * f1 / (r1 * f1 + f2) : 0, pbb = r2 * f2 + f1 > 0 ? r2 * f2 / (r2 * f2 + f1) : 0;
      const double den = r1 * f1 * f1 + 2 * f1 * f2 + r2 * f2 * f2, F1 = den > 0 ? (r1 * f1 * f1 + f1 * f2) / den : f1;
      auto u01 = [&] { return double(rng() >> 11) * (1.0 / 9007199254740992.0); };
      for (int i = 0; i < n; ++i) {
        if (i == 0) s[0] = u01() < F1 ? 0 : 1 % nu;
        else s[size_t(i)] = s[size_t(i - 1)] == 0 ? (u01() < paa ? 0 : 1 % nu) : (u01() < pbb ? 1 % nu : 0);
      }
      break;
    }
    case Sequence::Pattern: {
      std::vector<int> p;
      for (char c : spec.pattern)
        if (c >= 'A' && c <= 'Z' && c - 'A' < nu) p.push_back(c - 'A');
        else if (c >= 'a' && c <= 'z' && c - 'a' < nu) p.push_back(c - 'a');
      if (p.empty()) p.push_back(0);
      for (int i = 0; i < n; ++i) s[size_t(i)] = p[size_t(i) % p.size()];
      break;
    }
  }
  return s;
}

MolGraph chain_graph(const ChainSpec& spec, const std::vector<int>& seq) {
  std::vector<Parsed> units;
  for (const auto& u : spec.units) units.push_back(parse_unit(u.smiles));
  MolGraph g;
  int prev_tail = -1;
  for (size_t k = 0; k < seq.size(); ++k) {
    const Parsed& P = units.at(size_t(seq[k]));
    std::vector<int> map(P.g.atoms.size(), -1);
    for (size_t i = 0; i < P.g.atoms.size(); ++i) {
      if (int(i) == P.d1 || int(i) == P.d2) continue;
      MolAtom a = P.g.atoms[i];
      a.order.clear();
      a.chiral = 0;
      map[i] = int(g.atoms.size());
      g.atoms.push_back(a);
    }
    for (const auto& b : P.g.bonds)
      if (map[size_t(b.a)] >= 0 && map[size_t(b.b)] >= 0) g.bonds.push_back({map[size_t(b.a)], map[size_t(b.b)], b.order, 0});
    const int head = map[size_t(P.head)], tail = map[size_t(P.tail)];
    if (prev_tail >= 0) g.bonds.push_back({prev_tail, head, 1, 0});
    else if (g.atoms[size_t(head)].bracket) g.atoms[size_t(head)].hcount += 1;   // the head cap
    prev_tail = tail;
    if (k + 1 == seq.size() && g.atoms[size_t(tail)].bracket) g.atoms[size_t(tail)].hcount += 1;   // the tail cap
  }
  g.heavy = int(g.atoms.size());
  g.smiles = write_smiles(g);
  add_hydrogens(g);
  return g;
}

double chain_mass(const ChainSpec& spec, const std::vector<int>& seq) {
  double m = 0;
  for (const auto& a : chain_graph(spec, seq).atoms) m += element(a.element).mass;
  return m;
}

namespace {
System grow_chains_once(const ChainSpec& spec, const GrowOptions& o, GrowReport* report) {
  if (spec.units.empty()) throw GrowError("no repeat unit");
  GrowReport rep;
  std::vector<Template> T;
  for (size_t k = 0; k < spec.units.size(); ++k) {
    try {
      T.push_back(make_template(spec.units[k].name.empty() ? spec.units[k].smiles : spec.units[k].name, spec.units[k].smiles, spec.forcefield,
                                o.seed + k, rep.notes));
    } catch (const std::exception& e) {
      throw GrowError(e.what());
    }
  }
  const int nchains = std::max(1, o.chains);
  std::mt19937_64 rng(o.seed * 0x9E3779B97F4A7C15ull + 7);
  std::uniform_real_distribution<double> U(0, 1);
  std::normal_distribution<double> Nd(0, 1);

  // sequences, tacticity, box
  std::vector<ChainState> C(static_cast<size_t>(nchains));
  double mass = 0;
  auto make_sequence = [&](ChainState& ch, int c, int dp) {
    if (dp > 0) {   // polydisperse or an arm: this chain's own length
      ChainSpec one = spec;
      one.dp = std::max(1, dp);
      ch.seq = chain_sequence(one, o.seed + uint64_t(c) * 101);
    } else {
      ch.seq = chain_sequence(spec, o.seed + uint64_t(c) * 101);
    }
    ch.mirror.resize(ch.seq.size());
    for (size_t k = 0; k < ch.seq.size(); ++k) {
      if (k == 0) ch.mirror[k] = U(rng) < 0.5;
      else if (spec.tacticity == Tacticity::Isotactic) ch.mirror[k] = ch.mirror[k - 1];
      else if (spec.tacticity == Tacticity::Syndiotactic) ch.mirror[k] = !ch.mirror[k - 1];
      else if (spec.p_mr >= 0 && spec.p_rm >= 0) {   // Markov: the next dyad depends on the previous one
        const double a = std::min(1.0, spec.p_mr), b = std::min(1.0, spec.p_rm);
        const bool prev_m = k >= 2 ? ch.mirror[k - 1] == ch.mirror[k - 2] : U(rng) < (a + b > 0 ? b / (a + b) : 0.5);
        const bool m = prev_m ? U(rng) >= a : U(rng) < b;
        ch.mirror[k] = m ? ch.mirror[k - 1] : !ch.mirror[k - 1];
      } else ch.mirror[k] = U(rng) < std::clamp(spec.pm, 0.0, 1.0) ? ch.mirror[k - 1] : !ch.mirror[k - 1];
    }
    mass += chain_mass(spec, ch.seq);
  };
  for (int c = 0; c < nchains; ++c) {
    C[size_t(c)].mol = c;
    make_sequence(C[size_t(c)], c, size_t(c) < spec.chain_dp.size() ? spec.chain_dp[size_t(c)] : 0);
  }
  // arms of branched molecules: each replaces a hydrogen (or the head cap) of its parent and loses its own head cap
  const Architecture arch = spec.architecture;
  int n_arms = 0;
  if (arch != Architecture::Linear) {
    std::mt19937_64 brng(o.seed * 0xD1B54A32D192ED03ull + 11);
    auto add_arm = [&](int parent, int unit, bool core, int dp, int after) {
      ChainState a;
      a.mol = C[size_t(parent)].mol;
      a.parent = parent, a.anchor_unit = unit, a.core = core, a.after = after;
      const int idx = int(C.size());
      make_sequence(a, idx, dp);
      mass -= 2 * element(1).mass;
      C.push_back(std::move(a));
      ++n_arms;
      return idx;
    };
    if (arch == Architecture::Star && (spec.arms < 3 || spec.arms > 4)) throw GrowError("a star has 3 or 4 arms on its core carbon");
    for (int m = 0; m < nchains; ++m) {
      const int len = int(C[size_t(m)].seq.size());
      int prev = -1;
      if (arch == Architecture::Star)
        for (int j = 1; j < spec.arms; ++j) prev = add_arm(m, 0, true, int(C[size_t(m)].seq.size()), prev);
      else if (arch == Architecture::Comb)
        for (int u = std::max(1, spec.spacing) - 1; u < len; u += std::max(1, spec.spacing)) prev = add_arm(m, u, false, std::max(1, spec.arm_dp), prev);
      else
        for (int u = 1; u + 1 < len; ++u)
          if (std::uniform_real_distribution<double>(0, 1)(brng) < spec.branch_probability) prev = add_arm(m, u, false, std::max(1, spec.arm_dp), prev);
    }
  }
  const bool ortho = o.cell[0] > 0 && o.cell[1] > 0 && o.cell[2] > 0;
  const bool film = o.z_hi > o.z_lo;
  std::array<double, 3> Lv = o.cell;
  if (!ortho) {
    double L = o.box;
    if (L <= 0) {
      if (o.density <= 0) throw GrowError("give a box edge or a density");
      L = std::cbrt(mass / (o.density * 0.602214076));
    }
    Lv = {L, L, L};
  }
  if (film && (o.z_lo < 0 || o.z_hi > Lv[2])) throw GrowError("the film heights lie outside the cell");
  const double L = Lv[0];
  const double vol = Lv[0] * Lv[1] * (film ? o.z_hi - o.z_lo : Lv[2]);
  rep.box = L;
  rep.density = mass / (0.602214076 * vol);
  const double scale = o.contact_scale > 0 ? o.contact_scale : 1.0;
  // contact limits as the polystyrene grower's (C–C 3.0, C–H 2.45, H–H 2.0 Å), from Bondi radii for other elements
  // (a negative element is an arm's reserved head: a carbon with room for its hydrogens, 0.5 Å more)
  auto limit = [&](int a, int b) {
    const double extra = a < 0 || b < 0 ? 0.5 : 0.0;
    a = std::abs(a), b = std::abs(b);
    return scale * (0.88 * (element(a).vdw + element(b).vdw) - 0.08 * ((a == 1) + (b == 1)) + extra);
  };
  Cell3 cell;
  cell.init(Lv, scale * 0.86 * 2 * 2.3);
  if (o.substrate)
    for (size_t i = 0; i < o.substrate->atoms.size(); ++i) cell.add(o.substrate->atoms[i].pos, o.substrate->atoms[i].element, -1, int(i));
  // a film: heights outside [z_lo, z_hi] count as contacts
  const bool sphere = o.sphere_radius > 0;
  const Vec3 sc{o.sphere_centre[0], o.sphere_centre[1], o.sphere_centre[2]};
  auto region = [&](const Vec3& p) {
    double m = film ? std::min(p[2] - o.z_lo, o.z_hi - p[2]) : 1e9;
    if (sphere) {
      const double d = norm(cell.mi(p - sc));
      m = std::min(m, o.sphere_outside ? d - o.sphere_radius : o.sphere_radius - d);
    }
    return m;
  };
  const bool gauche = o.curve;
  const int trials = std::max(4, o.trials);
  rep.worst_margin = 1e9;

  // arms' anchors, all of a parent's at once when it is complete: the atom each arm's head bonds to and the valence it
  // takes (a hydrogen, which is dropped, or the head cap's); a stand-in carbon keeps each head's place, so the arms
  // grown first leave room for those after them
  // room for an arm's head carbon at P + d·1.53 and for the chain beyond it at P + d·2.9 (distance − contact limit, the
  // parent's atoms within two bonds of P left out)
  auto anchor_room = [&](int parent, int P, const Vec3& d) {
    const auto& par = C[size_t(parent)];
    std::set<int> near_p{P};
    for (int u : par.adj[size_t(P)]) {
      near_p.insert(u);
      for (int v : par.adj[size_t(u)]) near_p.insert(v);
    }
    double room = 1e9;
    for (double r : {1.53, 2.9}) {
      const Vec3 x = par.pos[size_t(P)] + d * r;
      cell.near(x, [&](int id) {
        if (cell.chain[size_t(id)] == parent && near_p.count(cell.local[size_t(id)])) return;
        room = std::min(room, norm(cell.mi(x - cell.x[size_t(id)])) - limit(6, cell.z[size_t(id)]));
      });
    }
    return room;
  };
  // one arm's anchor: a star arm on the parent's head atom (the head cap's valence first, then its hydrogens); a side
  // chain on the hydrogen, of the head or tail atom of its unit or the units beside it, with the most room
  auto resolve_one = [&](int e) {
    auto& ch = C[size_t(e)];
    auto& par = C[size_t(ch.parent)];
    int P = -1, h = -1;
    if (ch.core && par.parent < 0 && !par.cap_taken) {
      P = par.unit_start[0];
      par.cap_taken = true;
      ch.anchor_dir = unitv(par.pos[2] - par.pos[size_t(P)]);
    } else {
      std::set<int> taken;   // units already carrying a side chain of this parent
      for (size_t f = 0; f < C.size(); ++f)
        if (int(f) != e && C[f].parent == ch.parent && C[f].resolved && !C[f].core) taken.insert(C[f].anchor_unit);
      const int nu = int(par.unit_start.size());
      std::vector<int> units{ch.anchor_unit};
      if (!ch.core)
        for (int du : {-1, 1})
          if (ch.anchor_unit + du >= 0 && ch.anchor_unit + du < nu) units.push_back(ch.anchor_unit + du);
      double best_room = -1e9;
      int best_unit = ch.anchor_unit;
      for (int u : units) {
        if (!ch.core && taken.count(u)) continue;
        const int base = par.unit_start[size_t(u)];
        const Template& tp = T[size_t(par.seq[size_t(u)])];
        for (int cand : ch.core ? std::vector<int>{base} : std::vector<int>{base, base + tp.tail})
          for (int w : par.adj[size_t(cand)]) {
            if (w < 3 || par.z[size_t(w)] != 1 || par.dropped.count(w) || ch.tried.count(w)) continue;
            // the unit asked for first: a neighbour must have clearly more room to take its place
            const double room = anchor_room(ch.parent, cand, unitv(par.pos[size_t(w)] - par.pos[size_t(cand)])) - (u == ch.anchor_unit ? 0.0 : 0.3);
            if (room > best_room) best_room = room, P = cand, h = w, best_unit = u;
          }
      }
      if (P < 0)
        throw GrowError(ch.core ? "a star of " + std::to_string(spec.arms) + " arms needs " + std::to_string(spec.arms - 2) + " hydrogens on the head atom (" +
                                      element(par.z[size_t(par.unit_start[0])]).symbol + ") of the unit; choose a unit whose head is CH₂ or CH₃, or 3 arms"
                                : "no room for a side chain near unit " + std::to_string(ch.anchor_unit + 1) + " (no free hydrogen on its head or tail atom)");
      ch.anchor_unit = best_unit;
      ch.drop_h = h;
      par.dropped.insert(h);
      if (par.gid[size_t(h)] >= 0) cell.kill(par.gid[size_t(h)]), par.gid[size_t(h)] = -1;
      ch.anchor_dir = unitv(par.pos[size_t(h)] - par.pos[size_t(P)]);
    }
    ch.anchor = P;
    ch.reserve = cell.add(par.pos[size_t(P)] + ch.anchor_dir * 1.53, -6, -3, e);
    ch.resolved = true;
  };
  // all of a parent's arms at once when it is complete, so stand-in carbons keep each head's place and the arms grown
  // first leave room for those after them
  auto resolve_anchors = [&](int parent) {
    for (size_t e = 0; e < C.size(); ++e)
      if (C[e].parent == parent && !C[e].resolved) resolve_one(int(e));
  };
  // the parent lost the atoms its arms hang on (it backed up or restarted): the arms' places are chosen again later
  auto unresolve_arms = [&](int parent) {
    auto& par = C[size_t(parent)];
    bool any = false;
    for (auto& e : C)
      if (e.parent == parent && e.resolved) {
        if (e.reserve >= 0) cell.kill(e.reserve), e.reserve = -1;
        e.resolved = false, e.drop_h = -1, e.tried.clear();
        any = true;
      }
    if (any) par.dropped.clear(), par.cap_taken = false;
  };
  // an arm that finds no room gives its hydrogen back and takes the next best site (at most three times)
  auto reanchor = [&](int e) {
    auto& ch = C[size_t(e)];
    auto& par = C[size_t(ch.parent)];
    if (ch.core || ch.drop_h < 0 || ch.tried.size() >= 3) return false;
    ch.tried.insert(ch.drop_h);
    par.dropped.erase(ch.drop_h);
    par.gid[size_t(ch.drop_h)] = cell.add(par.pos[size_t(ch.drop_h)], 1, ch.parent, ch.drop_h);
    for (size_t i = 3; i < ch.gid.size(); ++i)
      if (ch.gid[i] >= 0) cell.kill(ch.gid[i]), ch.gid[i] = -1;
    ch.resolved = false;
    ch.started = false;
    resolve_one(e);
    ++rep.restarts;
    ch.starts = 0;
    return true;
  };
  // an arm starting: its stand-in goes, and the bond distances from its anchor to the molecule's atoms within four
  // bonds (other arms' stand-ins count as bonded to their anchors) exempt or soften those contacts
  auto begin_arm = [&](int ci) {
    auto& ch = C[size_t(ci)];
    if (!ch.resolved) resolve_anchors(ch.parent);
    if (ch.reserve >= 0) cell.kill(ch.reserve), ch.reserve = -1;
    std::map<std::pair<int, int>, int> seen{{{ch.parent, ch.anchor}, 0}};
    std::vector<std::pair<int, int>> q{{ch.parent, ch.anchor}};
    for (size_t qi = 0; qi < q.size(); ++qi) {
      const auto [c, i] = q[qi];
      const int d = seen[{c, i}];
      if (d == 4) continue;
      std::vector<std::pair<int, int>> nb;
      for (int w : C[size_t(c)].adj[size_t(i)])
        if (w >= 3 && !C[size_t(c)].dropped.count(w)) nb.push_back({c, w});
      if (C[size_t(c)].parent >= 0 && i == 3) nb.push_back({C[size_t(c)].parent, C[size_t(c)].anchor});
      for (size_t e = 0; e < C.size(); ++e)
        if (C[e].parent == c && C[e].anchor == i && int(e) != ci && C[e].done) nb.push_back({int(e), 3});
      for (const auto& w : nb)
        if (seen.emplace(w, d + 1).second) q.push_back(w);
    }
    ch.pdist.clear();
    for (const auto& [node, d] : seen)
      if (const int g = C[size_t(node.first)].gid[size_t(node.second)]; g >= 0) ch.pdist[g] = d;
    for (const auto& e : C)
      if (e.reserve >= 0 && e.mol == ch.mol)
        if (auto it = seen.find({e.parent, e.anchor}); it != seen.end() && it->second < 4) ch.pdist[e.reserve] = it->second + 1;
    ch.started = true;
  };

  auto start_chain = [&](int ci) {
    auto& ch = C[size_t(ci)];
    if (n_arms) unresolve_arms(ci);
    for (size_t i = 3; i < ch.gid.size(); ++i)
      if (ch.gid[i] >= 0) cell.kill(ch.gid[i]);
    ch.pos.clear(), ch.z.clear(), ch.tparent.clear(), ch.gid.clear(), ch.unit_of.clear(), ch.sp2.clear(), ch.donor.clear(), ch.backbone.clear(), ch.adj.clear(), ch.unit_start.clear();
    if (ch.parent >= 0) {   // an arm: the ghosts are its anchor atom and the two atoms before it in the parent
      if (!ch.started) begin_arm(ci);
      const auto& par = C[size_t(ch.parent)];
      // references for the first unit's torsions: a heavy atom bonded to the anchor (never the start ghost, which sits
      // where a head-cap arm's head goes) and one bonded to that
      int a1 = -1, a0 = -1;
      for (int w : par.adj[size_t(ch.anchor)])
        if (w >= 3 && par.z[size_t(w)] != 1 && (a1 < 0 || par.backbone[size_t(w)])) a1 = w;
      if (a1 >= 0)
        for (int w : par.adj[size_t(a1)])
          if (w >= 3 && w != ch.anchor && !par.dropped.count(w) && (a0 < 0 || par.z[size_t(w)] != 1)) a0 = w;
      const Vec3 P = par.pos[size_t(ch.anchor)], g1 = a1 >= 0 ? par.pos[size_t(a1)] : P - ch.anchor_dir * 1.53;
      Vec3 g0 = a0 >= 0 ? par.pos[size_t(a0)] : g1;
      if (norm(cross(g0 - g1, P - g1)) < 1e-6) {   // a straight or missing reference: any perpendicular
        Vec3 w{Nd(rng), Nd(rng), Nd(rng)};
        const Vec3 u = unitv(P - g1);
        g0 = g1 + unitv(w - u * dot(w, u)) * 1.53;
      }
      for (const Vec3& p : {g0, g1, P}) {
        ch.pos.push_back(p), ch.z.push_back(6), ch.gid.push_back(-1), ch.unit_of.push_back(-1), ch.sp2.push_back(0), ch.donor.push_back(0), ch.backbone.push_back(1), ch.adj.push_back({});
      }
      ch.tparent = {-1, 0, 1};
      ++ch.starts;
      ch.fails = 0;
      ch.backtracks = 0;
      return;
    }
    // three ghosts: a start point and a random frame (the head bonds to ghost 2)
    Vec3 s{U(rng) * Lv[0], U(rng) * Lv[1], film ? o.z_lo + 1 + U(rng) * std::max(0.0, o.z_hi - o.z_lo - 2) : U(rng) * Lv[2]};
    for (int tries = 0; sphere && region(s) < 1.0 && tries < 1000; ++tries)   // a start inside the allowed region
      s = {U(rng) * Lv[0], U(rng) * Lv[1], film ? o.z_lo + 1 + U(rng) * std::max(0.0, o.z_hi - o.z_lo - 2) : U(rng) * Lv[2]};
    Vec3 u{Nd(rng), Nd(rng), Nd(rng)};
    u = unitv(u);
    Vec3 w{Nd(rng), Nd(rng), Nd(rng)};
    w = unitv(w - u * dot(w, u));
    const Vec3 g2 = s, g1 = s + u * 1.53, g0 = g1 + unitv(u * 0.4 + w) * 1.53;
    for (const Vec3& p : {g0, g1, g2}) {
      ch.pos.push_back(p), ch.z.push_back(6), ch.gid.push_back(-1), ch.unit_of.push_back(-1), ch.sp2.push_back(0), ch.donor.push_back(0), ch.backbone.push_back(1), ch.adj.push_back({});
    }
    ch.tparent = {-1, 0, 1};
    ++ch.starts;
    ch.fails = 0;
    ch.backtracks = 0;
  };

  auto remove_units = [&](int ci, int count) {
    auto& ch = C[size_t(ci)];
    if (n_arms && count >= int(ch.unit_start.size())) unresolve_arms(ci);   // a star's core unit goes too
    for (int r = 0; r < count && !ch.unit_start.empty(); ++r) {
      const int from = ch.unit_start.back();
      ch.unit_start.pop_back();
      for (int i = int(ch.pos.size()) - 1; i >= from; --i) {
        if (ch.gid[size_t(i)] >= 0) cell.kill(ch.gid[size_t(i)]);
        for (int w : ch.adj[size_t(i)]) {
          auto& a = ch.adj[size_t(w)];
          a.erase(std::remove(a.begin(), a.end(), i), a.end());
        }
      }
      ch.pos.resize(size_t(from)), ch.z.resize(size_t(from)), ch.tparent.resize(size_t(from)), ch.gid.resize(size_t(from));
      ch.unit_of.resize(size_t(from)), ch.sp2.resize(size_t(from)), ch.donor.resize(size_t(from)), ch.backbone.resize(size_t(from)), ch.adj.resize(size_t(from));
    }
  };

  auto draw = [&](int kind, double t0) {
    switch (kind) {
      case 0: {
        // mostly staggered (trans, gauche±), sometimes anywhere: crowded units (methacrylates) need the room between
        if (gauche && U(rng) < 0.3) return (U(rng) * 2 - 1) * kPi;
        const double c = gauche ? std::array<double, 3>{kPi, kPi / 3, -kPi / 3}[size_t(U(rng) * 3) % 3] : kPi;
        return c + Nd(rng) * 15 * kPi / 180;
      }
      case 1: return t0 + std::floor(U(rng) * 12) * kPi / 6 + Nd(rng) * 5 * kPi / 180;
      case 2: return (U(rng) < 0.5 ? 0.0 : kPi) + Nd(rng) * 8 * kPi / 180;
      default: return kPi + Nd(rng) * 8 * kPi / 180;
    }
  };

  // direction of the free valence at unit u's tail (none when the tail has no child: then it is sampled)
  auto free_valence = [&](int ci, int u) -> std::optional<Vec3> {
    const auto& ch = C[size_t(ci)];
    const Template& t = T[size_t(ch.seq[size_t(u)])];
    if (!t.tail_fixed) return std::nullopt;
    const int base = ch.unit_start[size_t(u)];
    auto abs = [&](int local) { return local >= 0 ? base + local : ch.tparent[size_t(base)]; };
    Vec3 e1, e2, e3;
    frame(ch.pos[size_t(base + t.tail)], ch.pos[size_t(abs(t.r1))], ch.pos[size_t(abs(t.r2))], e1, e2, e3);
    const double c3 = ch.mirror[size_t(u)] ? -t.tf[2] : t.tf[2];
    return unitv(e1 * t.tf[0] + e2 * t.tf[1] + e3 * c3);
  };

  // one growth step of chain ci: place its next unit (best of the trials) or back up
  auto advance = [&](int ci) {
    auto& ch = C[size_t(ci)];
    const int k = int(ch.unit_start.size());
    const Template& t = T[size_t(ch.seq[size_t(k)])];
    const bool mir = ch.mirror[size_t(k)];
    const int base = int(ch.pos.size());
    const int prev_tail = k == 0 ? 2 : ch.unit_start.back() + T[size_t(ch.seq[size_t(k - 1)])].tail;
    const double link = k == 0 ? t.head_bond : t.head_bond + T[size_t(ch.seq[size_t(k - 1)])].tail_bond - 1.53;
    const double link_angle = k == 0 ? 1.95 : T[size_t(ch.seq[size_t(k - 1)])].tail_angle;
    auto absref = [&](int local) {
      if (local >= 0) return base + local;
      int a = prev_tail;
      for (int s = -1; s > local; --s) a = ch.tparent[size_t(a)];
      return a;
    };
    const int pt = prev_tail, ptp = ch.tparent[size_t(pt)];
    std::optional<Vec3> next_dir;
    if (k > 0) next_dir = free_valence(ci, k - 1);
    else if (ch.parent >= 0) next_dir = ch.anchor_dir;   // an arm's head takes the parent's freed valence
    const int root_kind = torsion_kind(ch.sp2[size_t(pt)], ch.donor[size_t(pt)], ch.sp2[size_t(ptp)], ch.donor[size_t(ptp)]);
    const int link_kind = torsion_kind(t.sp2[0], t.donor[0], ch.sp2[size_t(pt)], ch.donor[size_t(pt)]);
    // atoms within three bonds of each new atom: the unit's own bonds plus the link, walked from the new atom
    std::vector<std::vector<int>> nadj(size_t(t.n));
    for (const auto& b : t.bonds) nadj[size_t(b[0])].push_back(b[1]), nadj[size_t(b[1])].push_back(b[0]);
    // bond distance (up to four) from each new atom to chain-local atoms: within three bonds exempt, four (1-5 pairs,
    // such as a side group and the next backbone carbon, always close on a crowded backbone) at 85 % of the limit
    std::vector<std::map<int, int>> excl(size_t(t.n));
    std::map<int, int> near_tail;                   // bond distance + 1 from the tail: the look-ahead atom's
    for (int a = 0; a < t.n; ++a) {
      const int depth = 4;
      std::vector<std::pair<int, int>> st{{base + a, 0}};
      std::set<int> seen{base + a};
      std::map<int, int> dist{{base + a, 0}};
      while (!st.empty()) {
        // breadth first, so each atom's distance is its shortest
        std::sort(st.begin(), st.end(), [](const auto& x, const auto& y) { return x.second > y.second; });
        auto [u, d] = st.back();
        st.pop_back();
        if (d == depth) continue;
        std::vector<int> next;
        if (u >= base) {
          for (int w : nadj[size_t(u - base)]) next.push_back(base + w);
          if (u == base) next.push_back(pt);
        } else {
          next = ch.adj[size_t(u)];
          if (u == pt) next.push_back(base);
        }
        for (int w : next)
          if (w >= 3 && seen.insert(w).second) st.push_back({w, d + 1}), dist[w] = d + 1;
      }
      if (a == t.tail)
        for (const auto& [w, d] : dist)
          if (d <= 3) near_tail[w] = d + 1;
      excl[size_t(a)] = std::move(dist);
    }
    // an arm: bond distance from each new atom to the arm's head (then + 1 to the anchor), for the parent's atoms nearby
    std::vector<int> dhead(size_t(t.n), 99);
    int dhead_next = 99;
    if (ch.parent >= 0) {
      for (int a = 0; a < t.n; ++a)
        if (auto it = excl[size_t(a)].find(3); it != excl[size_t(a)].end()) dhead[size_t(a)] = it->second;
      if (auto it = near_tail.find(3); it != near_tail.end()) dhead_next = it->second;
    }
    std::vector<Vec3> best, trial(size_t(t.n));
    double best_m = -1e9;
    std::vector<double> gv(t.group_kind.size()), best_gv;
    double best_root = 0;
    // places the unit's atoms for a root torsion and group torsions, and returns the worst contact margin (Å); stops
    // counting once the margin is below `floor`
    auto score = [&](double root_t, const std::vector<double>& gv, double floor) {
      const double best_m = floor;
      auto P = [&](int local) -> Vec3 {
        const int a = absref(local);
        return a >= base ? trial[size_t(a - base)] : ch.pos[size_t(a)];
      };
      for (int a = 0; a < t.n; ++a) {
        const auto& r = t.ref[size_t(a)];
        if (a == 0 && next_dir) {   // the previous tail's free valence
          trial[0] = ch.pos[size_t(pt)] + *next_dir * link;
          continue;
        }
        double tor, bond = t.bond[size_t(a)], ang = t.ang[size_t(a)];
        if (a == 0) tor = root_t, bond = link, ang = link_angle;
        else if (t.group[size_t(a)] >= 0) tor = gv[size_t(t.group[size_t(a)])] + t.offset[size_t(a)];
        else tor = t.tor[size_t(a)];
        if (mir) tor = -tor;
        trial[size_t(a)] = place(P(r[2]), P(r[1]), P(r[0]), bond, ang, tor);
      }
      // worst margin against everything placed, and within the unit beyond three bonds
      double worst = 1e9;
      for (int a = 0; a < t.n && worst > best_m; ++a) {
        const Vec3& x = trial[size_t(a)];
        worst = std::min(worst, region(x));
        cell.near(x, [&](int id) {
          double f = 1.0;
          if (cell.chain[size_t(id)] == ci) {
            const auto it = excl[size_t(a)].find(cell.local[size_t(id)]);
            if (it != excl[size_t(a)].end()) {
              if (it->second <= 3) return;
              f = 0.85;
            }
          } else if (dhead[size_t(a)] <= 3) {   // across an arm's junction
            const auto it = ch.pdist.find(id);
            if (it != ch.pdist.end()) {
              const int tot = dhead[size_t(a)] + 1 + it->second;   // branch points are crowded, as quaternary carbons are
              if (tot <= 3) return;
              f = tot == 4 ? 0.85 : 0.9;
            }
          }
          const double d = norm(cell.mi(x - cell.x[size_t(id)]));
          const double m = d - f * limit(t.z[size_t(a)], cell.z[size_t(id)]);
          if (m < worst) worst = m;
        });
        for (int b = a + 1; b < t.n; ++b)
          if (auto it = excl[size_t(a)].find(base + b); it == excl[size_t(a)].end() || it->second > 3) {
            const double f = it == excl[size_t(a)].end() ? 1.0 : 0.85;
            const double m = norm(cell.mi(x - trial[size_t(b)])) - f * limit(t.z[size_t(a)], t.z[size_t(b)]);
            if (m < worst) worst = m;
          }
      }
      // look ahead: the next unit's head must have room at this unit's free valence (else the chain folds into itself)
      if (worst > best_m && t.tail_fixed && k + 1 < int(ch.seq.size())) {
        auto tp = [&](int local) { return local >= 0 ? trial[size_t(local)] : ch.pos[size_t(pt)]; };
        Vec3 e1, e2, e3;
        frame(trial[size_t(t.tail)], tp(t.r1), tp(t.r2), e1, e2, e3);
        const Vec3 dir = unitv(e1 * t.tf[0] + e2 * t.tf[1] + e3 * (mir ? -t.tf[2] : t.tf[2]));
        const Vec3 nx = trial[size_t(t.tail)] + dir * 1.53;
        worst = std::min(worst, region(nx));
        const auto& ex = near_tail;
        cell.near(nx, [&](int id) {
          double f = 1.0;
          if (cell.chain[size_t(id)] == ci) {
            const auto it = ex.find(cell.local[size_t(id)]);
            if (it != ex.end()) {
              if (it->second <= 3) return;
              f = 0.85;
            }
          } else if (dhead_next <= 3) {
            const auto it = ch.pdist.find(id);
            if (it != ch.pdist.end()) {
              const int tot = dhead_next + 1 + it->second;
              if (tot <= 3) return;
              f = tot == 4 ? 0.85 : 0.9;
            }
          }
          const double m = norm(cell.mi(nx - cell.x[size_t(id)])) - f * limit(6, cell.z[size_t(id)]);
          if (m < worst) worst = m;
        });
        for (int b = 0; b < t.n; ++b)
          if (auto it = ex.find(base + b); it == ex.end() || it->second > 3) {
            const double m = norm(cell.mi(nx - trial[size_t(b)])) - (it == ex.end() ? 1.0 : 0.85) * limit(6, t.z[size_t(b)]);
            if (m < worst) worst = m;
          }
      }
      return worst;
    };
    for (int tr = 0; tr < trials; ++tr) {
      const double root_t = draw(root_kind, 0);
      for (size_t gi = 0; gi < gv.size(); ++gi) gv[gi] = draw(int(gi) == t.link_group ? link_kind : t.group_kind[gi], 0);
      // the template's own torsion is the reference for groups next to an sp2 atom
      for (int a = 1; a < t.n; ++a)
        if (t.group[size_t(a)] >= 0 && t.group_kind[size_t(t.group[size_t(a)])] == 1 && t.offset[size_t(a)] == 0.0) gv[size_t(t.group[size_t(a)])] += t.tor[size_t(a)];
      const double worst = score(root_t, gv, best_m);
      if (worst > best_m) {
        best_m = worst;
        best = trial;
        best_gv = gv;
        best_root = root_t;
        if (o.comfortable > 0 && worst >= o.comfortable) break;
      }
    }
    // just short of the limits (long flexible units, crowded junctions): nudge single torsions of the best trial
    if (best_m < 0.05 && best_m > -0.8 && !best.empty()) {
      for (int it = 0; it < 160 && best_m < 0.05; ++it) {
        std::vector<double> g = best_gv;
        double r = best_root;
        const int moves = U(rng) < 0.3 ? 2 : 1;
        for (int mv = 0; mv < moves; ++mv) {
          const size_t pick = size_t(U(rng) * double(g.size() + 1));
          const double step = Nd(rng) * 15 * kPi / 180;
          if (pick >= g.size()) r += step;
          else g[pick] += step;
        }
        const double worst = score(r, g, best_m);
        if (worst > best_m) {
          best_m = worst;
          best = trial;
          best_gv = std::move(g);
          best_root = r;
        }
      }
    }
    if (best_m < o.accept) {
      ++ch.fails;
      ++rep.backtracks;
      const int limit_bt = o.max_backtracks > 0 ? o.max_backtracks : 40 * std::max(1, spec.dp);
      if (++ch.backtracks > limit_bt || k == 0) {
        if (ch.parent >= 0 && ch.starts > std::min(o.max_restarts, 60) && reanchor(ci)) {
          start_chain(ci);
          return;
        }
        if (ch.starts > o.max_restarts)
          throw GrowError((ch.parent >= 0 ? "an arm of molecule " + std::to_string(ch.mol + 1) : "chain " + std::to_string(ci + 1)) + " could not be placed after " +
                          std::to_string(o.max_restarts) + " restarts; lower the density or the contact scale");
        ++rep.restarts;
        start_chain(ci);
        return;
      }
      remove_units(ci, std::min(k, 1 + ch.fails / 4));
      return;
    }
    ch.fails = std::max(0, ch.fails - 1);
    rep.worst_margin = std::min(rep.worst_margin, best_m);
    ch.unit_start.push_back(base);
    for (int a = 0; a < t.n; ++a) {
      const int idx = base + a;
      ch.pos.push_back(best[size_t(a)]);
      ch.z.push_back(t.z[size_t(a)]);
      ch.tparent.push_back(a == 0 ? pt : base + t.parent[size_t(a)]);
      ch.unit_of.push_back(k);
      ch.sp2.push_back(t.sp2[size_t(a)]);
      ch.donor.push_back(t.donor[size_t(a)]);
      ch.backbone.push_back(a == 0 || a == t.tail);
      ch.adj.push_back({});
      ch.gid.push_back(cell.add(best[size_t(a)], t.z[size_t(a)], ci, idx));
    }
    for (const auto& b : t.bonds) ch.adj[size_t(base + b[0])].push_back(base + b[1]), ch.adj[size_t(base + b[1])].push_back(base + b[0]);
    if (k > 0) ch.adj[size_t(base)].push_back(pt), ch.adj[size_t(pt)].push_back(base);
    if (k == 0 && arch == Architecture::Star && ch.parent < 0) resolve_anchors(ci);   // the core's places, before the chain folds back
    if (int(ch.unit_start.size()) == int(ch.seq.size())) {
      ch.done = true;
      if (n_arms) resolve_anchors(ci);   // its arms' places are kept from now on
    }
  };

  for (int c = 0; c < nchains; ++c) start_chain(c);
  const int nstates = int(C.size());
  int finished = 0;
  long units_total = 0;
  for (const auto& ch : C) units_total += long(ch.seq.size());
  auto last_snap = std::chrono::steady_clock::now() - std::chrono::hours(1);
  auto snapshot = [&] {   // the chains so far, for a live view
    System p;
    p.cell.a = {Lv[0], 0, 0}, p.cell.b = {0, Lv[1], 0}, p.cell.c = {0, 0, Lv[2]};
    p.has_mol = true;
    GrowOptions::Live live;
    live.chains = nchains, live.units_total = units_total, live.restarts = rep.restarts, live.worst_margin = rep.worst_margin;
    double placed = 0;
    std::vector<std::vector<uint32_t>> at(C.size());
    for (int c = 0; c < nstates; ++c) {
      const auto& ch = C[size_t(c)];
      live.chains_done += ch.done && ch.parent < 0 ? 1 : 0;
      live.units += long(ch.unit_start.size());
      auto& m = at[size_t(c)];
      m.assign(ch.pos.size(), UINT32_MAX);
      for (size_t i = 3; i < ch.pos.size(); ++i) {
        if (ch.dropped.count(int(i))) continue;
        Atom a;
        a.element = ch.z[i], a.pos = ch.pos[i], a.mol = ch.mol + 1, a.id = int64_t(p.atoms.size() + 1);
        m[i] = uint32_t(p.atoms.size());
        p.atoms.push_back(a);
        placed += element(ch.z[i]).mass;
      }
      for (size_t i = 3; i < ch.adj.size(); ++i)
        for (int j : ch.adj[i])
          if (j > int(i) && j >= 3 && m[i] != UINT32_MAX && m[size_t(j)] != UINT32_MAX) p.bonds.push_back({m[i], m[size_t(j)], 1});
      if (ch.parent >= 0 && ch.pos.size() > 3 && at[size_t(ch.parent)].size() > size_t(ch.anchor)) p.bonds.push_back({at[size_t(ch.parent)][size_t(ch.anchor)], m[3], 1});
    }
    p.bonds_from_file = true;
    live.density = placed / (6.02214076e23 * vol * 1e-24);
    o.snapshot(p, live);
  };
  while (finished < nstates) {
    for (int c = 0; c < nstates; ++c) {
      auto& ch = C[size_t(c)];
      if (ch.done) continue;
      if (ch.parent >= 0) {   // an arm waits for its parent and the arm before it
        if (!C[size_t(ch.parent)].done || (ch.after >= 0 && !C[size_t(ch.after)].done)) continue;
        if (!ch.started) start_chain(c);
      }
      advance(c);
    }
    finished = int(std::count_if(C.begin(), C.end(), [](const ChainState& s) { return s.done; }));
    if (o.progress && !o.progress(finished, nstates, rep.restarts)) throw GrowError("cancelled");
    if (o.snapshot && (finished == nstates || std::chrono::steady_clock::now() - last_snap > std::chrono::duration<double>(o.snapshot_seconds))) {
      snapshot();
      last_snap = std::chrono::steady_clock::now();
    }
  }

  // assemble: per chain its atoms (without the ghosts), the bonds, and a hydrogen cap on each end
  System s;
  s.title = "CAPS Grow: " + std::to_string(nchains) + " chains × " + std::to_string(spec.dp) + " units";
  s.cell.origin = {0, 0, 0};
  s.cell.a = {Lv[0], 0, 0};
  s.cell.b = {0, Lv[1], 0};
  s.cell.c = {0, 0, Lv[2]};
  s.unwrapped = true;
  s.has_mol = true;
  std::map<int, int> type_of;
  auto type_for = [&](int z) {
    auto it = type_of.find(z);
    if (it != type_of.end()) return it->second;
    const int t = int(type_of.size()) + 1;
    type_of[z] = t;
    TypeInfo ti;
    ti.type = t;
    ti.mass = element(z).mass;
    ti.label = element(z).symbol;
    s.types.push_back(ti);
    return t;
  };
  auto add = [&](int z, const Vec3& p, int mol) {
    Atom a;
    a.id = int64_t(s.atoms.size() + 1);
    a.mol = mol;
    a.element = z;
    a.type = type_for(z);
    a.name = element(z).symbol;
    a.pos = p;
    s.atoms.push_back(a);
    return uint32_t(s.atoms.size() - 1);
  };
  // the substrate first, keeping its molecule ids (a slab or filler is molecule 1; an earlier blend component keeps its
  // chains); the new chains are numbered after them
  // residue names: the unit's name, three letters upper case (PDB), else U + its letter
  std::vector<std::string> unit_code;
  for (size_t u = 0; u < spec.units.size(); ++u) {
    std::string code;
    for (char ch : spec.units[u].name)
      if (std::isalnum(static_cast<unsigned char>(ch)) && code.size() < 3) code += char(std::toupper(static_cast<unsigned char>(ch)));
    unit_code.push_back(code.empty() ? std::string("U") + char('A' + int(u % 26)) : code);
  }
  int mol0 = 0;
  if (o.substrate) {
    int64_t top = 0;
    for (const auto& a : o.substrate->atoms) top = std::max(top, a.mol);
    for (const auto& a : o.substrate->atoms) add(a.element, a.pos, top > 0 ? int(a.mol) : 1);
    for (const auto& b : o.substrate->bonds) s.bonds.push_back(b);
    mol0 = int(std::max<int64_t>(top, 1));
  }
  // the atoms molecule by molecule: each main chain, then its arms (residues numbered on along the molecule)
  std::vector<std::vector<uint32_t>> maps(C.size());
  std::vector<int> resid_off(size_t(nchains), 0);
  std::vector<int> order;
  for (int m = 0; m < nchains; ++m) {
    order.push_back(m);
    for (int c = nchains; c < nstates; ++c)
      if (C[size_t(c)].mol == m) order.push_back(c);
  }
  for (int c : order) {
    auto& ch = C[size_t(c)];
    const int molid = ch.mol + 1 + mol0, roff = resid_off[size_t(ch.mol)];
    auto& map = maps[size_t(c)];
    map.assign(ch.pos.size(), UINT32_MAX);
    for (size_t i = 3; i < ch.pos.size(); ++i) {
      if (ch.dropped.count(int(i))) continue;   // a hydrogen an arm replaced
      map[i] = add(ch.z[i], ch.pos[i], molid);
      if (const int k = ch.unit_of[i]; k >= 0) {   // residues: one per repeat unit, numbered along the chain
        s.atoms[map[i]].resid = roff + k + 1;
        s.atoms[map[i]].resname = unit_code[size_t(ch.seq[size_t(k)])];
      }
    }
    for (size_t k = 0; k < ch.unit_start.size(); ++k) {
      const Template& t = T[size_t(ch.seq[k])];
      const int base = ch.unit_start[k];
      for (const auto& b : t.bonds)
        if (map[size_t(base + b[0])] != UINT32_MAX && map[size_t(base + b[1])] != UINT32_MAX) s.bonds.push_back({map[size_t(base + b[0])], map[size_t(base + b[1])], b[2]});
      if (k > 0) {
        const int pt = ch.unit_start[k - 1] + T[size_t(ch.seq[k - 1])].tail;
        s.bonds.push_back({map[size_t(pt)], map[size_t(base)], 1});
      }
    }
    // caps: the head's toward the start ghost (an arm's head bonds to its anchor instead), the tail's where a next
    // unit would go
    const int head = ch.unit_start.front();
    if (ch.parent >= 0) {
      s.bonds.push_back({maps[size_t(ch.parent)][size_t(ch.anchor)], map[size_t(head)], 1});
    } else if (!ch.cap_taken) {
      const Vec3 hp = ch.pos[size_t(head)] + unitv(ch.pos[2] - ch.pos[size_t(head)]) * 1.09;
      s.bonds.push_back({map[size_t(head)], add(1, hp, molid), 1});
      s.atoms.back().resid = roff + 1, s.atoms.back().resname = unit_code[size_t(ch.seq.front())];
    }
    const Template& tl = T[size_t(ch.seq.back())];
    const int tail = ch.unit_start.back() + tl.tail;
    const int tp = ch.tparent[size_t(tail)], tgp = ch.tparent[size_t(tp)];
    const auto fv = free_valence(c, int(ch.unit_start.size()) - 1);
    const Vec3 tpos = fv ? ch.pos[size_t(tail)] + *fv * 1.09 : place(ch.pos[size_t(tgp)], ch.pos[size_t(tp)], ch.pos[size_t(tail)], 1.09, tl.tail_angle, kPi);
    s.bonds.push_back({map[size_t(tail)], add(1, tpos, molid), 1});
    s.atoms.back().resid = roff + int64_t(ch.seq.size()), s.atoms.back().resname = unit_code[size_t(ch.seq.back())];
    resid_off[size_t(ch.mol)] += int(ch.seq.size());
  }
  s.bonds_from_file = true;
  rep.chains_placed = nchains;
  char cb[96];
  if (ortho) std::snprintf(cb, sizeof cb, "cell %.2f × %.2f × %.2f Å", Lv[0], Lv[1], Lv[2]);
  else std::snprintf(cb, sizeof cb, "box %.3f Å", L);
  std::string units_text = std::to_string(nchains) + " chains × " + std::to_string(spec.dp) + " units";
  if (arch == Architecture::Star) units_text = std::to_string(nchains) + " stars of " + std::to_string(spec.arms) + " arms × " + std::to_string(spec.dp) + " units";
  else if (arch != Architecture::Linear)
    units_text = std::to_string(nchains) + (arch == Architecture::Comb ? " combs" : " branched chains") + ": backbones of " + std::to_string(spec.dp) + " units, " +
                 std::to_string(n_arms) + " side chains of " + std::to_string(std::max(1, spec.arm_dp)) + " units";
  if (!spec.chain_dp.empty()) {
    int lo = 1 << 30, hi = 0;
    long total = 0;
    for (const auto& ch : C) lo = std::min(lo, int(ch.seq.size())), hi = std::max(hi, int(ch.seq.size())), total += long(ch.seq.size());
    units_text = std::to_string(nchains) + " chains, " + std::to_string(lo) + "–" + std::to_string(hi) + " units (Σ " + std::to_string(total) + ")";
  }
  rep.notes.insert(rep.notes.begin(), units_text + " · " + std::to_string(s.atoms.size()) + " atoms · " + cb +
                                          " · " + std::to_string(rep.density).substr(0, 5) + " g/cm³" + (film ? " in the film" : ""));
  if (o.substrate) rep.notes.push_back(std::to_string(o.substrate->atoms.size()) + " substrate atoms kept fixed while growing (molecule 1)");
  if (report) *report = rep;
  return s;
}

}  // namespace

System grow_chains(const ChainSpec& spec, const GrowOptions& o, GrowReport* report) {
  if (!o.auto_scale) return grow_chains_once(spec, o, report);
  std::vector<double> scales = {o.contact_scale > 0 ? o.contact_scale : 1.0};
  for (double x : {0.85, 0.75, 0.7, 0.6})
    if (x < scales.front() - 1e-9) scales.push_back(x);
  for (size_t k = 0; k < scales.size(); ++k) {
    GrowOptions g = o;
    g.auto_scale = false;
    g.contact_scale = scales[k];
    if (k + 1 < scales.size()) g.max_restarts = std::min(o.max_restarts, 10);   // give up early on the stricter limits
    try {
      System s = grow_chains_once(spec, g, report);
      if (k > 0 && report) {
        char n[200];
        std::snprintf(n, sizeof n, "grown at contact scale %.2f (%.2f was too crowded); relax with push-off before dynamics", scales[k], scales[0]);
        report->notes.push_back(n);
      }
      return s;
    } catch (const GrowError& e) {
      if (k + 1 == scales.size() || std::string(e.what()) == "cancelled") throw;
    }
  }
  throw GrowError("no contact scale worked");
}

// ---------------------------------------------------------------- chain-length distributions

namespace {
struct Draw {
  std::mt19937_64 rng;
  double u01() { return (double(rng() >> 11) + 0.5) * (1.0 / 9007199254740992.0); }
  double normal() { return std::sqrt(-2 * std::log(u01())) * std::cos(2 * M_PI * u01()); }
  // Marsaglia & Tsang (2000); shape < 1 by the boost x = G(shape + 1) · U^(1/shape)
  double gamma(double k) {
    if (k < 1) return gamma(k + 1) * std::pow(u01(), 1 / k);
    const double d = k - 1.0 / 3, c = 1 / std::sqrt(9 * d);
    for (;;) {
      double x, v;
      do { x = normal(); v = 1 + c * x; } while (v <= 0);
      v = v * v * v;
      const double u = u01();
      if (u < 1 - 0.0331 * x * x * x * x || std::log(u) < 0.5 * x * x + d * (1 - v + std::log(v))) return d * v;
    }
  }
  int poisson(double lam) {   // Knuth for small means, a normal approximation above 60
    if (lam > 60) return std::max(0, int(std::lround(lam + std::sqrt(lam) * normal())));
    const double L = std::exp(-lam);
    int k = 0;
    double p = 1;
    do { ++k; p *= u01(); } while (p > L);
    return k - 1;
  }
};
}  // namespace

std::vector<int> draw_chain_lengths(const std::string& dist, double nn, double pdi, int count, uint64_t seed) {
  if (count < 1) return {};
  nn = std::max(2.0, nn);
  Draw d{std::mt19937_64(seed * 0x9E3779B97F4A7C15ull + 31)};
  std::vector<int> out;
  for (int i = 0; i < count; ++i) {
    double x = nn;
    if (dist == "schulz-zimm") {
      const double k = pdi > 1.0001 ? 1 / (pdi - 1) : 1e4;
      x = d.gamma(k) * nn / k;
    } else if (dist == "flory") {   // geometric with mean nn: P(N) = p (1 − p)^(N − 1)
      const double p = 1 / nn;
      x = 1 + std::floor(std::log(d.u01()) / std::log(1 - p));
    } else if (dist == "poisson") {
      x = 1 + d.poisson(nn - 1);
    } else if (dist != "monodisperse") {
      throw std::invalid_argument("chain lengths: monodisperse, schulz-zimm, flory or poisson");
    }
    out.push_back(std::max(2, int(std::lround(x))));
  }
  return out;
}

double chain_length_pdf(const std::string& dist, double nn, double pdi, double n) {
  if (n <= 0) return 0;
  if (dist == "schulz-zimm") {
    const double k = pdi > 1.0001 ? 1 / (pdi - 1) : 1e4, th = nn / k;
    return std::exp((k - 1) * std::log(n) - n / th - std::lgamma(k) - k * std::log(th));
  }
  if (dist == "flory") return (1 / nn) * std::pow(1 - 1 / nn, n - 1);
  if (dist == "poisson") { const double lam = nn - 1, m = std::round(n) - 1; return m < 0 ? 0 : std::exp(m * std::log(lam) - lam - std::lgamma(m + 1)); }
  return std::fabs(n - nn) < 0.5 ? 1.0 : 0.0;
}

}  // namespace caps
