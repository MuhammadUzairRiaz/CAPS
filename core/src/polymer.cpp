// CAPS polymer builder: chains of any repeat unit grown into a periodic cell (see polymer.hpp).
#include "caps/polymer.hpp"

#include <algorithm>
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
    const auto f = molecule_forcefield(g, ff, "gasteiger", fn);
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
  double L = 0, cs = 4.0;
  int nc = 1;
  std::vector<std::vector<int>> bins;
  std::vector<Vec3> x;        // by id
  std::vector<int> z, chain, local;
  std::vector<char> alive;
  void init(double edge, double reach) {
    L = edge;
    nc = std::max(1, int(std::floor(L / std::max(reach, 1.0))));
    cs = L / nc;
    bins.assign(size_t(nc) * nc * nc, {});
  }
  int bin_of(const Vec3& p) const {
    auto w = [&](double v) {
      int i = int(std::floor((v - L * std::floor(v / L)) / cs));
      return std::clamp(i, 0, nc - 1);
    };
    return (w(p[0]) * nc + w(p[1])) * nc + w(p[2]);
  }
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
    for (int k = 0; k < 3; ++k) d[k] -= L * std::round(d[k] / L);
    return d;
  }
  template <class F>
  void near(const Vec3& p, F&& f) const {
    if (nc < 3) {
      for (size_t id = 0; id < x.size(); ++id)
        if (alive[id]) f(int(id));
      return;
    }
    auto w = [&](double v) { return std::clamp(int(std::floor((v - L * std::floor(v / L)) / cs)), 0, nc - 1); };
    const int bx = w(p[0]), by = w(p[1]), bz = w(p[2]);
    for (int i = -1; i <= 1; ++i)
      for (int j = -1; j <= 1; ++j)
        for (int k = -1; k <= 1; ++k) {
          const int xi = (bx + i + nc) % nc, yj = (by + j + nc) % nc, zk = (bz + k + nc) % nc;
          for (int id : bins[size_t((xi * nc + yj) * nc + zk)]) f(id);
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
};

}  // namespace

Sequence sequence_from_string(const std::string& s) {
  if (s == "alternating") return Sequence::Alternating;
  if (s == "block") return Sequence::Block;
  if (s == "random") return Sequence::Random;
  if (s == "gradient") return Sequence::Gradient;
  if (s == "pattern") return Sequence::Pattern;
  return Sequence::Homopolymer;
}

const char* to_string(Sequence s) {
  switch (s) {
    case Sequence::Alternating: return "alternating";
    case Sequence::Block: return "block";
    case Sequence::Random: return "random";
    case Sequence::Gradient: return "gradient";
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

System grow_chains(const ChainSpec& spec, const GrowOptions& o, GrowReport* report) {
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
  for (int c = 0; c < nchains; ++c) {
    auto& ch = C[size_t(c)];
    ch.seq = chain_sequence(spec, o.seed + uint64_t(c) * 101);
    ch.mirror.resize(ch.seq.size());
    for (size_t k = 0; k < ch.seq.size(); ++k) {
      if (k == 0) ch.mirror[k] = U(rng) < 0.5;
      else if (spec.tacticity == Tacticity::Isotactic) ch.mirror[k] = ch.mirror[k - 1];
      else if (spec.tacticity == Tacticity::Syndiotactic) ch.mirror[k] = !ch.mirror[k - 1];
      else ch.mirror[k] = U(rng) < std::clamp(spec.pm, 0.0, 1.0) ? ch.mirror[k - 1] : !ch.mirror[k - 1];
    }
    mass += chain_mass(spec, ch.seq);
  }
  double L = o.box;
  if (L <= 0) {
    if (o.density <= 0) throw GrowError("give a box edge or a density");
    L = std::cbrt(mass / (o.density * 0.602214076));
  }
  rep.box = L;
  rep.density = mass / (0.602214076 * L * L * L);
  const double scale = o.contact_scale > 0 ? o.contact_scale : 1.0;
  // contact limits as the polystyrene grower's (C–C 3.0, C–H 2.45, H–H 2.0 Å), from Bondi radii for other elements
  auto limit = [&](int a, int b) { return scale * (0.88 * (element(a).vdw + element(b).vdw) - 0.08 * ((a == 1) + (b == 1))); };
  Cell3 cell;
  cell.init(L, scale * 0.86 * 2 * 2.3);
  const bool gauche = o.curve;
  const int trials = std::max(4, o.trials);
  rep.worst_margin = 1e9;

  auto start_chain = [&](int ci) {
    auto& ch = C[size_t(ci)];
    for (size_t i = 3; i < ch.gid.size(); ++i)
      if (ch.gid[i] >= 0) cell.kill(ch.gid[i]);
    ch.pos.clear(), ch.z.clear(), ch.tparent.clear(), ch.gid.clear(), ch.unit_of.clear(), ch.sp2.clear(), ch.donor.clear(), ch.backbone.clear(), ch.adj.clear(), ch.unit_start.clear();
    // three ghosts: a start point and a random frame (the head bonds to ghost 2)
    const Vec3 s{U(rng) * L, U(rng) * L, U(rng) * L};
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
    std::vector<Vec3> best, trial(size_t(t.n));
    double best_m = -1e9;
    int dbg_a = -1, dbg_b = -1;
    std::vector<double> gv(t.group_kind.size());
    for (int tr = 0; tr < trials; ++tr) {
      const double root_t = draw(root_kind, 0);
      for (size_t gi = 0; gi < gv.size(); ++gi) gv[gi] = draw(int(gi) == t.link_group ? link_kind : t.group_kind[gi], 0);
      // the template's own torsion is the reference for groups next to an sp2 atom
      for (int a = 1; a < t.n; ++a)
        if (t.group[size_t(a)] >= 0 && t.group_kind[size_t(t.group[size_t(a)])] == 1 && t.offset[size_t(a)] == 0.0) gv[size_t(t.group[size_t(a)])] += t.tor[size_t(a)];
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
      int wa = -1, wb = -1;
      for (int a = 0; a < t.n && worst > best_m; ++a) {
        const Vec3& x = trial[size_t(a)];
        cell.near(x, [&](int id) {
          double f = 1.0;
          if (cell.chain[size_t(id)] == ci) {
            const auto it = excl[size_t(a)].find(cell.local[size_t(id)]);
            if (it != excl[size_t(a)].end()) {
              if (it->second <= 3) return;
              f = 0.85;
            }
          }
          const double d = norm(cell.mi(x - cell.x[size_t(id)]));
          const double m = d - f * limit(t.z[size_t(a)], cell.z[size_t(id)]);
          if (m < worst) worst = m, wa = a, wb = cell.chain[size_t(id)] == ci ? cell.local[size_t(id)] : -1000 - cell.chain[size_t(id)];
        });
        for (int b = a + 1; b < t.n; ++b)
          if (auto it = excl[size_t(a)].find(base + b); it == excl[size_t(a)].end() || it->second > 3) {
            const double f = it == excl[size_t(a)].end() ? 1.0 : 0.85;
            const double m = norm(cell.mi(x - trial[size_t(b)])) - f * limit(t.z[size_t(a)], t.z[size_t(b)]);
            if (m < worst) worst = m, wa = a, wb = base + b;
          }
      }
      // look ahead: the next unit's head must have room at this unit's free valence (else the chain folds into itself)
      if (worst > best_m && t.tail_fixed && k + 1 < int(ch.seq.size())) {
        auto tp = [&](int local) { return local >= 0 ? trial[size_t(local)] : ch.pos[size_t(pt)]; };
        Vec3 e1, e2, e3;
        frame(trial[size_t(t.tail)], tp(t.r1), tp(t.r2), e1, e2, e3);
        const Vec3 dir = unitv(e1 * t.tf[0] + e2 * t.tf[1] + e3 * (mir ? -t.tf[2] : t.tf[2]));
        const Vec3 nx = trial[size_t(t.tail)] + dir * 1.53;
        const auto& ex = near_tail;
        cell.near(nx, [&](int id) {
          double f = 1.0;
          if (cell.chain[size_t(id)] == ci) {
            const auto it = ex.find(cell.local[size_t(id)]);
            if (it != ex.end()) {
              if (it->second <= 3) return;
              f = 0.85;
            }
          }
          const double m = norm(cell.mi(nx - cell.x[size_t(id)])) - f * limit(6, cell.z[size_t(id)]);
          if (m < worst) worst = m, wa = -2, wb = cell.local[size_t(id)];
        });
        for (int b = 0; b < t.n; ++b)
          if (auto it = ex.find(base + b); it == ex.end() || it->second > 3) {
            const double m = norm(cell.mi(nx - trial[size_t(b)])) - (it == ex.end() ? 1.0 : 0.85) * limit(6, t.z[size_t(b)]);
            if (m < worst) worst = m, wa = -2, wb = base + b;
          }
      }
      if (worst > best_m) {
        best_m = worst;
        dbg_a = wa, dbg_b = wb;
        best = trial;
        if (o.comfortable > 0 && worst >= o.comfortable) break;
      }
    }
    if (best_m < o.accept) {
      ++ch.fails;
      ++rep.backtracks;
      const int limit_bt = o.max_backtracks > 0 ? o.max_backtracks : 40 * std::max(1, spec.dp);
      if (++ch.backtracks > limit_bt || k == 0) {
        if (ch.starts > o.max_restarts) throw GrowError("chain " + std::to_string(ci + 1) + " could not be placed after " + std::to_string(o.max_restarts) +
                                                        " restarts; lower the density or the contact scale");
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
    if (int(ch.unit_start.size()) == int(ch.seq.size())) ch.done = true;
  };

  for (int c = 0; c < nchains; ++c) start_chain(c);
  int finished = 0;
  while (finished < nchains) {
    for (int c = 0; c < nchains; ++c)
      if (!C[size_t(c)].done) advance(c);
    finished = int(std::count_if(C.begin(), C.end(), [](const ChainState& s) { return s.done; }));
    if (o.progress && !o.progress(finished, nchains, rep.restarts)) throw GrowError("cancelled");
  }

  // assemble: per chain its atoms (without the ghosts), the bonds, and a hydrogen cap on each end
  System s;
  s.title = "CAPS Grow: " + std::to_string(nchains) + " chains × " + std::to_string(spec.dp) + " units";
  s.cell.origin = {0, 0, 0};
  s.cell.a = {L, 0, 0};
  s.cell.b = {0, L, 0};
  s.cell.c = {0, 0, L};
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
  for (int c = 0; c < nchains; ++c) {
    auto& ch = C[size_t(c)];
    std::vector<uint32_t> map(ch.pos.size());
    for (size_t i = 3; i < ch.pos.size(); ++i) map[i] = add(ch.z[i], ch.pos[i], c + 1);
    for (size_t k = 0; k < ch.unit_start.size(); ++k) {
      const Template& t = T[size_t(ch.seq[k])];
      const int base = ch.unit_start[k];
      for (const auto& b : t.bonds) s.bonds.push_back({map[size_t(base + b[0])], map[size_t(base + b[1])], b[2]});
      if (k > 0) {
        const int pt = ch.unit_start[k - 1] + T[size_t(ch.seq[k - 1])].tail;
        s.bonds.push_back({map[size_t(pt)], map[size_t(base)], 1});
      }
    }
    // caps: the head's toward the start ghost, the tail's where a next unit would go
    const int head = ch.unit_start.front();
    const Vec3 hp = ch.pos[size_t(head)] + unitv(ch.pos[2] - ch.pos[size_t(head)]) * 1.09;
    s.bonds.push_back({map[size_t(head)], add(1, hp, c + 1), 1});
    const Template& tl = T[size_t(ch.seq.back())];
    const int tail = ch.unit_start.back() + tl.tail;
    const int tp = ch.tparent[size_t(tail)], tgp = ch.tparent[size_t(tp)];
    const auto fv = free_valence(c, int(ch.unit_start.size()) - 1);
    const Vec3 tpos = fv ? ch.pos[size_t(tail)] + *fv * 1.09 : place(ch.pos[size_t(tgp)], ch.pos[size_t(tp)], ch.pos[size_t(tail)], 1.09, tl.tail_angle, kPi);
    s.bonds.push_back({map[size_t(tail)], add(1, tpos, c + 1), 1});
  }
  s.bonds_from_file = true;
  rep.chains_placed = nchains;
  rep.notes.insert(rep.notes.begin(), std::to_string(nchains) + " chains × " + std::to_string(spec.dp) + " units · " + std::to_string(s.atoms.size()) + " atoms · box " +
                                          std::to_string(L).substr(0, 6) + " Å · " + std::to_string(rep.density).substr(0, 5) + " g/cm³");
  if (report) *report = rep;
  return s;
}

}  // namespace caps
