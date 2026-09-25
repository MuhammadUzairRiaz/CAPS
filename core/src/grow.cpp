// CAPS Grow — polystyrene grown atom group by atom group inside a periodic cubic cell.
//
// Geometry (internal coordinates, NeRF placement): backbone C–C 1.53 Å, C–C–C 114°, C(H)–C(ipso) 1.51 Å,
// planar phenyl ring C–C 1.39 Å, aromatic C–H 1.08 Å, aliphatic C–H 1.09 Å, tetrahedral substituents.
// Each step places the next backbone carbon together with the substituents of the current one and keeps the
// trial (backbone torsion 180° or ±60°, ±15° jitter; ring rotation) whose worst non-bonded margin is largest.
// Margins are distance minus a contact limit (C–C 3.0, C–H 2.45, H–H 2.0 Å) against every atom already in the
// cell, other chains included, using the minimum image; pairs within three bonds are exempt. Poor steps back
// track a few units; a chain that cannot continue restarts elsewhere.
#include "caps/grow.hpp"

#include <cstdio>
#include <algorithm>
#include <cmath>
#include <random>

#include "caps/elements.hpp"
#include "caps/typing.hpp"

namespace caps {
namespace {

constexpr double kPi = 3.14159265358979323846;
double rad(double d) { return d * kPi / 180.0; }
Vec3 unit(const Vec3& v) { return v * (1.0 / norm(v)); }

// NeRF: position of d given a-b-c, |cd| = bond, angle bcd = ang (deg), torsion abcd = tor (deg).
Vec3 place(const Vec3& a, const Vec3& b, const Vec3& c, double bond, double ang, double tor) {
  const Vec3 bc = unit(c - b);
  const Vec3 n = unit(cross(b - a, bc));
  const Vec3 m = cross(n, bc);
  const double A = rad(ang), T = rad(tor);
  const double d0 = -bond * std::cos(A), d1 = bond * std::sin(A) * std::cos(T), d2 = bond * std::sin(A) * std::sin(T);
  return c + bc * d0 + m * d1 + n * d2;
}

// Unit vector of a tetrahedral substituent on an atom with neighbour directions u1, u2; s = ±1 picks the side.
Vec3 substituent(const Vec3& u1, const Vec3& u2, int s) {
  const Vec3 bis = unit(u1 + u2);
  const Vec3 nn = unit(cross(u1, u2));
  return unit(bis * -0.575 + nn * (0.818 * s));
}

double base_limit(char a, char b) {
  if (a == 'C' && b == 'C') return 3.0;
  if (a == 'H' && b == 'H') return 2.0;
  return 2.45;
}

struct Cand {
  char e;
  Vec3 p;
  std::vector<int> placed_nb;   // bonded atoms already in the cell
  std::vector<int> local_nb;    // bonded partners among this step's candidates
};

// Everything placed so far, all chains, with a periodic cell list.
struct Cell3 {
  double L;
  int n;
  double cs;
  std::vector<Vec3> p;
  std::vector<char> e;
  std::vector<std::vector<int>> adj;
  std::vector<std::vector<int>> bins;

  double scale = 1.0;
  double limit(char a, char b) const { return base_limit(a, b) * scale; }

  explicit Cell3(double L_, double scale_) : L(L_), scale(scale_) {
    n = std::max(1, int(L / 3.0));
    cs = L / n;
    bins.assign(size_t(n) * n * n, {});
  }
  int k1(double x) const { return ((int(std::floor(x / cs)) % n) + n) % n; }
  size_t key(const Vec3& q) const { return (size_t(k1(q[0])) * n + k1(q[1])) * n + k1(q[2]); }
  double dist(const Vec3& a, const Vec3& b) const {
    double s = 0;
    for (int k = 0; k < 3; ++k) {
      double d = b[k] - a[k];
      d -= L * std::round(d / L);
      s += d * d;
    }
    return std::sqrt(s);
  }
  int add(char el, const Vec3& q, const std::vector<int>& bonded) {
    const int idx = int(p.size());
    p.push_back(q);
    e.push_back(el);
    adj.emplace_back(bonded);
    for (int b : bonded) adj[b].push_back(idx);
    bins[key(q)].push_back(idx);
    return idx;
  }
  void remove(int idx) {
    auto& b = bins[key(p[size_t(idx)])];
    b.erase(std::remove(b.begin(), b.end(), idx), b.end());
    for (int j : adj[size_t(idx)]) {
      auto& a = adj[size_t(j)];
      a.erase(std::remove(a.begin(), a.end(), idx), a.end());
    }
    adj[size_t(idx)].clear();
  }
  void truncate(size_t size) {
    for (size_t i = size; i < p.size(); ++i) {
      auto& b = bins[key(p[i])];
      b.erase(std::remove_if(b.begin(), b.end(), [&](int t) { return size_t(t) >= size; }), b.end());
    }
    p.resize(size);
    e.resize(size);
    adj.resize(size);
    for (auto& a : adj) a.erase(std::remove_if(a.begin(), a.end(), [&](int t) { return size_t(t) >= size; }), a.end());
  }
  template <class F>
  void around(const Vec3& q, F&& f) const {
    const int x = k1(q[0]), y = k1(q[1]), z = k1(q[2]);
    int xs[3], ys[3], zs[3], cx = 0, cy = 0, cz = 0;
    auto gather = [&](int v, int* out, int& c) {
      for (int d = -1; d <= 1; ++d) {
        const int u = ((v + d) % n + n) % n;
        bool dup = false;
        for (int t = 0; t < c; ++t) dup |= out[t] == u;
        if (!dup) out[c++] = u;
      }
    };
    gather(x, xs, cx);
    gather(y, ys, cy);
    gather(z, zs, cz);
    for (int a = 0; a < cx; ++a)
      for (int b = 0; b < cy; ++b)
        for (int c = 0; c < cz; ++c)
          for (int t : bins[(size_t(xs[a]) * n + ys[b]) * n + zs[c]]) f(t);
  }

  // Worst (distance − limit) over the candidate set, exempting pairs within three bonds.
  double score(const std::vector<Cand>& cand) const {
    double worst = 9.0;
    std::vector<int> excl, frontier, next;
    for (const auto& c : cand) {
      excl = c.placed_nb;
      for (int l : c.local_nb) excl.insert(excl.end(), cand[size_t(l)].placed_nb.begin(), cand[size_t(l)].placed_nb.end());
      frontier = excl;
      for (int hop = 0; hop < 2; ++hop) {
        next.clear();
        for (int x : frontier)
          for (int y : adj[size_t(x)])
            if (std::find(excl.begin(), excl.end(), y) == excl.end()) { excl.push_back(y); next.push_back(y); }
        frontier.swap(next);
      }
      around(c.p, [&](int t) {
        const double v = dist(c.p, p[size_t(t)]) - limit(c.e, e[size_t(t)]);
        if (v < worst && std::find(excl.begin(), excl.end(), t) == excl.end()) worst = v;
      });
    }
    return worst;
  }
};

struct Rng {
  std::mt19937_64 g;
  explicit Rng(uint64_t s) : g(s) {}
  double uniform(double a, double b) { return std::uniform_real_distribution<double>(a, b)(g); }
  bool coin() { return std::uniform_int_distribution<int>(0, 1)(g) == 1; }
  int below(int n) { return std::uniform_int_distribution<int>(0, n - 1)(g); }
  Vec3 direction() {
    std::normal_distribution<double> N(0, 1);
    return unit(Vec3{N(g), N(g), N(g)});
  }
};

struct Sub {
  char e;
  Vec3 p;
  bool aromatic;
};

// One polystyrene chain: backbone, substituents per backbone atom, virtual neighbours for the end caps.
struct ChainGeom {
  std::vector<Vec3> B;
  std::vector<std::vector<Sub>> S;
  Vec3 first_virt, last_virt;
};

void ring(const Vec3& bk, const Vec3& d, double phi, std::vector<Sub>& out) {
  const Vec3 cen = bk + d * (1.51 + 1.39);
  const Vec3 ref = unit(cross(d, std::fabs(d[2]) < 0.9 ? Vec3{0, 0, 1} : Vec3{1, 0, 0}));
  const Vec3 w = unit(ref * std::cos(phi) + cross(d, ref) * std::sin(phi));
  for (int m = 0; m < 6; ++m) {
    const double a = kPi + m * kPi / 3;
    out.push_back({'C', cen + d * (1.39 * std::cos(a)) + w * (1.39 * std::sin(a)), true});
  }
  for (int m = 1; m < 6; ++m) {
    const double a = kPi + m * kPi / 3;
    out.push_back({'H', cen + d * (2.47 * std::cos(a)) + w * (2.47 * std::sin(a)), false});
  }
}

std::vector<Sub> substituents(int j, const std::vector<int>& sides, const Vec3& prev, const Vec3& bj, const Vec3& nxt, double phi) {
  const Vec3 u1 = unit(prev - bj), u2 = unit(nxt - bj);
  std::vector<Sub> out;
  if (j % 2 == 1) {
    ring(bj, substituent(u1, u2, sides[size_t(j)]), phi, out);
    out.push_back({'H', bj + substituent(u1, u2, -sides[size_t(j)]) * 1.09, false});
  } else {
    for (int s : {1, -1}) out.push_back({'H', bj + substituent(u1, u2, s) * 1.09, false});
  }
  return out;
}

// One chain growing in the shared cell. All chains advance one backbone step per round, so late chains are
// not left with only the gaps between finished ones. A chain backs up only its own atoms.
struct ChainGrower {
  const GrowOptions& o;
  Cell3& cell;
  Rng& rnd;
  int nb;
  std::vector<int> sides;
  std::vector<Vec3> B;
  std::vector<int> bid;
  std::vector<std::vector<Sub>> S;
  std::vector<std::vector<int>> added;   // atom indices added at each step (index nb holds the start atoms)
  std::vector<size_t> snapB;
  Vec3 first_virt{}, last_virt{};
  int jj = 0;
  int backtracks = 0;
  std::vector<int> fails;   // failed attempts per backbone step since the chain started
  bool started = false, done = false;

  ChainGrower(const GrowOptions& o_, Cell3& c, Rng& r) : o(o_), cell(c), rnd(r), nb(2 * o_.dp) {}

  void clear() {
    for (auto& v : added)
      for (int t : v) cell.remove(t);
    added.assign(static_cast<size_t>(nb) + 1, {});
    snapB.assign(static_cast<size_t>(nb) + 1, 0);
    S.assign(static_cast<size_t>(nb), {});
    B.clear();
    bid.clear();
    fails.assign(static_cast<size_t>(nb), 0);
    jj = 0;
    backtracks = 0;
    started = done = false;
  }

  bool start() {
    clear();
    sides.assign(static_cast<size_t>(nb), 0);
    for (int k = 1; k < nb; k += 2)
      sides[size_t(k)] = o.tacticity == Tacticity::Isotactic ? 1 : o.tacticity == Tacticity::Syndiotactic ? ((k / 2) % 2 == 0 ? 1 : -1) : (rnd.coin() ? 1 : -1);
    // The emptiest of several random points, then the best of several orientations for three backbone atoms.
    Vec3 c0{};
    double room = -1;
    for (int k = 0; k < 48; ++k) {
      const Vec3 q{rnd.uniform(0, cell.L), rnd.uniform(0, cell.L), rnd.uniform(0, cell.L)};
      double dmin = 9.0;
      cell.around(q, [&](int t) { dmin = std::min(dmin, cell.dist(q, cell.p[size_t(t)])); });
      if (dmin > room) { room = dmin; c0 = q; }
      if (room >= 9.0) break;
    }
    double best = -1e300;
    for (int k = 0; k < 24; ++k) {
      const Vec3 ex = rnd.direction();
      const Vec3 ey = unit(cross(ex, rnd.direction()));
      std::vector<Vec3> T = {c0, c0 + ex * 1.53};
      T.push_back(T[1] + ex * (1.53 * std::cos(rad(66))) + ey * (1.53 * std::sin(rad(66))));
      std::vector<Cand> c = {{'C', T[0], {}, {1}}, {'C', T[1], {}, {0, 2}}, {'C', T[2], {}, {1}}};
      const double sc = cell.score(c);
      if (sc > best) { best = sc; B = T; }
      if (sc >= 1.0) break;
    }
    if (best < 0) return false;
    first_virt = place(B[2], B[1], B[0], 1.53, 114.0, 180.0);
    auto& init = added[size_t(nb)];
    bid.push_back(cell.add('C', B[0], {}));
    bid.push_back(cell.add('C', B[1], {bid[0]}));
    bid.push_back(cell.add('C', B[2], {bid[1]}));
    init = bid;
    started = true;
    return true;
  }

  // One growth step (or a back-up). Returns false when this start has failed.
  bool step(double& worst_accepted, int& total_backtracks) {
    snapB[size_t(jj)] = B.size();
    const Vec3 prev = jj > 0 ? B[size_t(jj) - 1] : first_virt;
    const bool need_next = int(B.size()) <= jj + 1;   // next backbone atom (a look-ahead past the end at jj = nb − 1)
    double best_sc = -1e300;
    Vec3 best_next{};
    std::vector<Sub> best_sub;
    std::vector<Cand> best_cand;
    // Steps that failed before get more trials.
    const int trials = o.escalate ? o.trials * (1 + std::min(3, fails[size_t(jj)])) : o.trials;
    for (int trial = 0; trial < trials; ++trial) {
      Vec3 nxt;
      if (need_next) {   // jj >= 2 here
        const double t = ((!o.curve || rnd.uniform(0, 1) < 0.5) ? 180.0 : (rnd.coin() ? 60.0 : -60.0)) + rnd.uniform(-15, 15);
        nxt = place(B[size_t(jj) - 2], B[size_t(jj) - 1], B[size_t(jj)], 1.53, 114.0, t);
      } else {
        nxt = B[size_t(jj) + 1];
      }
      auto sb = substituents(jj, sides, prev, B[size_t(jj)], nxt, rnd.uniform(0, kPi));
      std::vector<Cand> cand;
      if (jj % 2 == 1) {
        for (int m = 0; m < 6; ++m) cand.push_back({'C', sb[size_t(m)].p, m == 0 ? std::vector<int>{bid[size_t(jj)]} : std::vector<int>{}, {(m + 5) % 6, (m + 1) % 6}});
        for (int m = 1; m < 6; ++m) cand.push_back({'H', sb[size_t(5 + m)].p, {}, {m}});
        cand.push_back({'H', sb[11].p, {bid[size_t(jj)]}, {}});
      } else {
        for (const auto& x : sb) cand.push_back({'H', x.p, {bid[size_t(jj)]}, {}});
      }
      if (need_next) cand.push_back({'C', nxt, {bid[size_t(jj)]}, {}});
      const double sc = cell.score(cand);
      if (sc > best_sc) { best_sc = sc; best_next = nxt; best_sub = sb; best_cand = cand; }
      if (sc >= o.comfortable) break;   // plenty of room: stop looking
    }
    if (best_sc < o.accept) {
      const int budget = o.max_backtracks > 0 ? o.max_backtracks : (o.escalate ? std::max(400, 40 * o.dp) : 400);
      if (jj >= 4 && backtracks < budget) {
        ++backtracks;
        ++total_backtracks;
        // Back up 2–5 steps, further each time the same step keeps failing.
        const int f = ++fails[size_t(jj)];
        const int to = jj - std::min(jj - 2, 2 + rnd.below(4) + (o.escalate ? 2 * (f - 1) : 0));
        for (int k = to; k < nb; ++k) {
          for (int t : added[size_t(k)]) cell.remove(t);
          added[size_t(k)].clear();
          S[size_t(k)].clear();
        }
        B.resize(snapB[size_t(to)]);
        bid.resize(snapB[size_t(to)]);
        jj = to;
        return true;
      }
      return false;
    }
    worst_accepted = std::min(worst_accepted, best_sc);
    std::vector<int> ids;
    for (const auto& c : best_cand) {
      std::vector<int> nbrs = c.placed_nb;
      for (int l : c.local_nb)
        if (size_t(l) < ids.size()) nbrs.push_back(ids[size_t(l)]);
      ids.push_back(cell.add(c.e, c.p, nbrs));
    }
    added[size_t(jj)] = ids;
    if (need_next) { B.push_back(best_next); bid.push_back(ids.back()); }
    S[size_t(jj)] = best_sub;
    if (++jj == nb) finish();
    return true;
  }

  void finish() {
    last_virt = B.size() > size_t(nb) ? B[size_t(nb)] : place(B[size_t(nb) - 3], B[size_t(nb) - 2], B[size_t(nb) - 1], 1.53, 114.0, 180.0);
    if (B.size() > size_t(nb)) {   // the look-ahead atom past the chain end is not part of the chain
      auto& last = added[size_t(nb) - 1];
      cell.remove(last.back());
      last.pop_back();
      B.resize(size_t(nb));
    }
    done = true;
  }

  ChainGeom geometry() const { return {B, S, first_virt, last_virt}; }
};

}  // namespace

Tacticity tacticity_from_string(const std::string& s) {
  if (s == "isotactic" || s == "iso") return Tacticity::Isotactic;
  if (s == "syndiotactic" || s == "syndio") return Tacticity::Syndiotactic;
  if (s == "atactic" || s == "a") return Tacticity::Atactic;
  throw std::invalid_argument("tacticity must be atactic, isotactic or syndiotactic");
}

const char* to_string(Tacticity t) {
  switch (t) {
    case Tacticity::Isotactic: return "isotactic";
    case Tacticity::Syndiotactic: return "syndiotactic";
    default: return "atactic";
  }
}

double chain_mass(const std::string& polymer, int dp) {
  if (polymer != "polystyrene") throw std::invalid_argument("unknown polymer '" + polymer + "' (this build has polystyrene)");
  // C8H8 per unit plus two H end caps.
  return dp * (8 * element(6).mass + 8 * element(1).mass) + 2 * element(1).mass;
}

double box_for_density(const GrowOptions& o) {
  if (o.density <= 0) throw std::invalid_argument("give a box edge or a positive target density");
  constexpr double kNA = 6.02214076e23;
  const double grams = o.chains * chain_mass(o.polymer, o.dp) / kNA;
  return std::cbrt(grams / o.density) * 1e8;   // cm -> Å
}

std::vector<double> gasteiger_ch(const System& s, const std::vector<char>& aromatic, int iterations) {
  // Gasteiger & Marsili, Tetrahedron 36, 3219 (1980): χ = a + b q + c q² per atom type (element and hybridisation),
  // charge moved along each bond by (χ_j − χ_i) / χ⁺ of the less electronegative atom, damped by ½ per iteration.
  // Parameters of the paper's Table 1 (H, C, N, O, F, Cl, Br, I) and the S sp³ extension; other elements are refused.
  struct P { double a, b, c; };
  static const P H{7.17, 6.24, -0.56}, C3{7.98, 9.18, 1.88}, C2{8.79, 9.32, 1.51}, C1{10.39, 9.45, 0.73},
      N3{11.54, 10.82, 1.36}, N2{12.87, 11.15, 0.85}, N1{15.68, 11.70, -0.27}, O3{14.18, 12.92, 1.39}, O2{17.07, 13.79, 0.47},
      F{14.66, 13.85, 2.31}, Cl{11.00, 9.69, 1.35}, Br{10.08, 8.47, 1.16}, I{9.90, 7.96, 0.96}, S3{10.14, 9.13, 1.38};
  const size_t n = s.atoms.size();
  const Perception pc = perceive(s);
  std::vector<const P*> t(n);
  for (size_t i = 0; i < n; ++i) {
    // hybridisation from the perceived bonds: a triple bond or two double bonds → sp, a double or aromatic bond → sp²
    int doubles = 0, triples = 0;
    bool arom = (i < aromatic.size() && aromatic[i]) || pc.aromatic[i];
    for (size_t k = 0; k < pc.nb[i].size(); ++k) {
      if (pc.arom_bond[i][k]) arom = true;
      else if (pc.order[i][k] == 2) ++doubles;
      else if (pc.order[i][k] == 3) ++triples;
    }
    const int hyb = triples > 0 || doubles > 1 ? 1 : doubles > 0 || arom ? 2 : 3;
    switch (s.atoms[i].element) {
      case 1: t[i] = &H; break;
      case 6: t[i] = hyb == 1 ? &C1 : hyb == 2 ? &C2 : &C3; break;
      case 7: t[i] = hyb == 1 ? &N1 : hyb == 2 ? &N2 : &N3; break;
      case 8: t[i] = hyb == 3 ? &O3 : &O2; break;
      case 9: t[i] = &F; break;
      case 17: t[i] = &Cl; break;
      case 35: t[i] = &Br; break;
      case 53: t[i] = &I; break;
      case 16:
        if (hyb != 3) throw std::invalid_argument("Gasteiger–Marsili parameters cover sp³ sulfur only; use QEq charges for S=O or thiophene");
        t[i] = &S3;
        break;
      default:
        throw std::invalid_argument(std::string("no Gasteiger–Marsili parameters for ") + element(s.atoms[i].element).symbol +
                                    " (H, C, N, O, F, Cl, Br, I and sp³ S); use QEq charges");
    }
  }
  const auto& nb = pc.nb;
  std::vector<double> q(n, 0.0), chi(n), dq(n);
  auto plus = [&](const P* p) { return p == &H ? 20.02 : p->a + p->b + p->c; };   // χ⁺ of the cation; H takes 20.02
  for (int k = 0; k < iterations; ++k) {
    for (size_t i = 0; i < n; ++i) chi[i] = t[i]->a + t[i]->b * q[i] + t[i]->c * q[i] * q[i];
    std::fill(dq.begin(), dq.end(), 0.0);
    const double damp = std::pow(0.5, k + 1);
    for (size_t i = 0; i < n; ++i)
      for (uint32_t j : nb[i]) dq[i] += (chi[j] - chi[i]) / plus(chi[j] > chi[i] ? t[i] : t[j]) * damp;
    for (size_t i = 0; i < n; ++i) q[i] += dq[i];
  }
  return q;
}

System grow(const GrowOptions& o, GrowReport* report) {
  if (o.polymer != "polystyrene") throw std::invalid_argument("unknown polymer '" + o.polymer + "' (this build has polystyrene)");
  if (o.chains < 1 || o.dp < 2) throw std::invalid_argument("need at least one chain of at least two units");
  const double L = o.box > 0 ? o.box : box_for_density(o);
  if (L < 12.0) throw std::invalid_argument("box edge below 12 Å is too small for a periodic cell with 3 Å contact limits");

  if (o.contact_scale < 0.5 || o.contact_scale > 1.2) throw std::invalid_argument("contact scale must be between 0.5 and 1.2");
  Cell3 cell(L, o.contact_scale);
  Rng rnd(o.seed);
  GrowReport rep;
  rep.box = L;
  double worst = 9.0;
  std::vector<ChainGrower> growers;
  growers.reserve(size_t(o.chains));
  for (int k = 0; k < o.chains; ++k) growers.emplace_back(o, cell, rnd);
  std::vector<int> restarts(size_t(o.chains), 0);
  // Rounds: every unfinished chain starts or advances one step.
  for (bool pending = true; pending;) {
    pending = false;
    if (o.progress) {
      int finished = 0;
      for (const auto& g : growers) finished += g.done;
      if (!o.progress(finished, o.chains, rep.restarts)) throw GrowError("cancelled");
    }
    for (int k = 0; k < o.chains; ++k) {
      auto& g = growers[size_t(k)];
      if (g.done) continue;
      pending = true;
      const bool ok = g.started ? g.step(worst, rep.backtracks) : g.start();
      if (!ok) {
        g.clear();
        ++rep.restarts;
        if (++restarts[size_t(k)] > o.max_restarts)
          throw GrowError("could not grow chain " + std::to_string(k + 1) + " of " + std::to_string(o.chains) + " in a " + std::to_string(L).substr(0, 5) +
                          " Å box after " + std::to_string(o.max_restarts) + " restarts; lower the density, the number of chains or the contact scale");
      }
    }
  }
  std::vector<ChainGeom> chains;
  for (const auto& g : growers) chains.push_back(g.geometry());
  rep.chains_placed = o.chains;

  // Assemble: per chain, backbone first, then substituents per backbone atom, then the two end-cap H.
  System s;
  s.title = "polystyrene " + std::to_string(o.chains) + " x DP " + std::to_string(o.dp) + " (" + to_string(o.tacticity) + ")";
  s.source_format = "caps-grow";
  s.cell.a = {L, 0, 0};
  s.cell.b = {0, L, 0};
  s.cell.c = {0, 0, L};
  s.has_mol = s.has_charges = s.bonds_from_file = s.unwrapped = true;
  s.types = {{1, element(6).mass, "c3"}, {2, element(6).mass, "ca"}, {3, element(1).mass, "hc"}, {4, element(1).mass, "ha"}};
  std::vector<char> aromatic;
  const int nb = 2 * o.dp;
  for (size_t ci = 0; ci < chains.size(); ++ci) {
    const auto& g = chains[ci];
    const uint32_t base = uint32_t(s.atoms.size());
    auto add = [&](int el, const Vec3& p, int type, const char* name, bool ar) {
      Atom a;
      a.id = int64_t(s.atoms.size() + 1);
      a.mol = int64_t(ci + 1);
      a.element = el;
      a.type = type;
      a.name = name;
      a.resname = "STY";
      a.pos = p;
      s.atoms.push_back(a);
      aromatic.push_back(ar);
    };
    for (int k = 0; k < nb; ++k) add(6, g.B[size_t(k)], 1, "c3", false);
    uint32_t x = base + uint32_t(nb);
    for (int k = 0; k < nb; ++k) {
      if (k % 2 == 1) {
        for (int m = 0; m < 6; ++m) add(6, g.S[size_t(k)][size_t(m)].p, 2, "ca", true);
        for (int m = 1; m < 6; ++m) add(1, g.S[size_t(k)][size_t(5 + m)].p, 4, "ha", false);
        add(1, g.S[size_t(k)][11].p, 3, "hc", false);
        for (int m = 0; m < 6; ++m) s.bonds.push_back({x + uint32_t(m), x + uint32_t((m + 1) % 6)});
        s.bonds.push_back({base + uint32_t(k), x});
        for (int m = 1; m < 6; ++m) s.bonds.push_back({x + uint32_t(m), x + 6 + uint32_t(m) - 1});
        s.bonds.push_back({base + uint32_t(k), x + 11});
        x += 12;
      } else {
        add(1, g.S[size_t(k)][0].p, 3, "hc", false);
        add(1, g.S[size_t(k)][1].p, 3, "hc", false);
        s.bonds.push_back({base + uint32_t(k), x});
        s.bonds.push_back({base + uint32_t(k), x + 1});
        x += 2;
      }
    }
    add(1, g.B[0] + unit(g.first_virt - g.B[0]) * 1.09, 3, "hc", false);
    add(1, g.B[size_t(nb) - 1] + unit(g.last_virt - g.B[size_t(nb) - 1]) * 1.09, 3, "hc", false);
    s.bonds.push_back({base, x});
    s.bonds.push_back({base + uint32_t(nb) - 1, x + 1});
    for (int k = 0; k + 1 < nb; ++k) s.bonds.push_back({base + uint32_t(k), base + uint32_t(k) + 1});
  }
  const auto q = gasteiger_ch(s, aromatic);
  for (size_t i = 0; i < s.atoms.size(); ++i) s.atoms[i].charge = q[i];

  rep.worst_margin = worst;
  rep.density = s.density();
  char buf[160];
  std::snprintf(buf, sizeof buf, "grown %d chains in a %.2f Å cell: density %.4f g/cm³, %d restarts, %d backtracks, worst contact margin %+.3f Å",
                o.chains, L, rep.density, rep.restarts, rep.backtracks, worst);
  rep.notes.push_back(buf);
  if (o.contact_scale < 1.0) {
    std::snprintf(buf, sizeof buf, "contact limits scaled to %.0f %% (C–C %.2f, C–H %.2f, H–H %.2f Å): minimise before dynamics", o.contact_scale * 100,
                  3.0 * o.contact_scale, 2.45 * o.contact_scale, 2.0 * o.contact_scale);
    rep.notes.push_back(buf);
  }
  s.notes = rep.notes;
  if (report) *report = rep;
  return s;
}

}  // namespace caps
