// CAPS powder indexing and direct-space structure solution (see powder_index.hpp).
#include "caps/powder_index.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <functional>
#include <random>
#include <set>
#include <stdexcept>

#include "caps/crystalfit.hpp"
#include "caps/spacegroup.hpp"

namespace caps {

namespace {
constexpr double kPi = 3.14159265358979323846;
constexpr double kDeg = 180.0 / kPi;
using I3 = std::array<int, 3>;
}  // namespace

// ---------------------------------------------------------------- peaks
std::vector<Peak> find_peaks(const MeasuredPattern& p, const PeakOptions& o) {
  const size_t n = p.y.size();
  if (n < 10) throw std::invalid_argument("peaks: too few points");
  const int h = std::max(0, o.smooth);
  std::vector<double> s(n, 0.0), bg(n, 0.0);
  for (size_t i = 0; i < n; ++i) {
    double sum = 0;
    int c = 0;
    for (int k = -h; k <= h; ++k) {
      const long j = long(i) + k;
      if (j >= 0 && j < long(n)) sum += p.y[size_t(j)], ++c;
    }
    s[i] = sum / c;
  }
  // background: the lowest smoothed value within ±1° 2θ, then smoothed again
  const double step = (p.two_theta.back() - p.two_theta.front()) / double(n - 1);
  const int w = std::max(3, int(1.0 / step));
  std::vector<double> mn(n);
  for (size_t i = 0; i < n; ++i) {
    double m = 1e300;
    for (long j = std::max(0L, long(i) - w); j <= std::min(long(n) - 1, long(i) + w); ++j) m = std::min(m, s[size_t(j)]);
    mn[i] = m;
  }
  for (size_t i = 0; i < n; ++i) {
    double sum = 0;
    int c = 0;
    for (long j = std::max(0L, long(i) - w); j <= std::min(long(n) - 1, long(i) + w); ++j) sum += mn[size_t(j)], ++c;
    bg[i] = sum / c;
  }
  double top = 0;
  for (size_t i = 0; i < n; ++i) top = std::max(top, s[i] - bg[i]);
  std::vector<Peak> pk;
  for (size_t i = 1; i + 1 < n; ++i) {
    const double v = s[i] - bg[i];
    if (v < o.min_height * top || s[i] < s[i - 1] || s[i] <= s[i + 1]) continue;
    // a parabola through the three points
    const double a = s[i - 1], b = s[i], c = s[i + 1], den = a - 2 * b + c;
    const double off = den != 0 ? 0.5 * (a - c) / den : 0;
    pk.push_back({p.two_theta[i] + std::clamp(off, -1.0, 1.0) * step, v});
  }
  if (o.drop_kalpha2) {
    std::vector<Peak> keep;
    for (const auto& q : pk) {
      bool partner = false;
      for (const auto& r : pk) {
        if (r.two_theta >= q.two_theta || r.height <= q.height) continue;
        const double sn = std::sin(r.two_theta / 2 / kDeg) * o.lambda2 / o.lambda1;
        if (sn >= 1) continue;
        const double t2 = 2 * std::asin(sn) * kDeg;
        if (std::fabs(t2 - q.two_theta) < std::max(0.05, 3 * step)) partner = true;
      }
      if (!partner) keep.push_back(q);
    }
    pk.swap(keep);
  }
  return pk;
}

// ---------------------------------------------------------------- indexing
namespace {
struct Box { double lo[4], hi[4]; };   // a, b, c (Å), β (degrees)
enum class Sys { cubic, tetragonal, hexagonal, orthorhombic, monoclinic };

Sys system_of(const std::string& s) {
  if (s == "cubic") return Sys::cubic;
  if (s == "tetragonal") return Sys::tetragonal;
  if (s == "hexagonal") return Sys::hexagonal;
  if (s == "orthorhombic") return Sys::orthorhombic;
  if (s == "monoclinic") return Sys::monoclinic;
  throw std::invalid_argument("index: crystal system '" + s + "' not known (cubic, tetragonal, hexagonal, orthorhombic or monoclinic)");
}

// A box is handed to least squares once no free dimension moves the top line's Q by more than this many tolerances
// across it: the box centre is then within a few tolerances of the cell (less on the lower lines, which lie further
// apart), so the nearest-index assignment converges; finer leaves only multiply the least-squares runs.
constexpr double kLeafTolerances = 8.0;
// a cell indexing every line with at least this M20 ends the search of the systems with more free parameters
constexpr double kStrongM20 = 100.0;

struct Indexer {
  const IndexOptions& o;
  std::string name;
  Sys sys = Sys::cubic;
  int m = 1;                            // reciprocal metric parameters (linear in Q)
  std::vector<int> free_dims;           // the box dimensions searched: a; a, c; a, b, c; a, b, c, β
  std::vector<double> Q, tol;           // observed Q (1/Å²) and their tolerances
  double qmax = 0;                      // the top line's Q plus its tolerance
  std::vector<I3> hkl;                  // the indices to try
  std::vector<std::array<double, 4>> co;   // per index: Q = Σ co[k] x[k] with the system's metric parameters x
  int H = 0;
  std::vector<int> lookup;              // (h + H, k, l) → index into hkl, or −1
  std::vector<IndexSolution> found;
  long tests = 0;
  double vol_limit = 0;                 // once a cell is found: boxes beyond about twice its volume hold only its supercells
  bool stopped = false;                 // the box budget ran out
  double frac = 0;                      // progress through the start boxes
  std::vector<std::pair<double, double>> rr;   // scratch: the Q ranges of the kept indices

  Indexer(const IndexOptions& opt, const std::string& s) : o(opt), name(s), sys(system_of(s)) {
    m = sys == Sys::cubic ? 1 : sys == Sys::orthorhombic ? 3 : sys == Sys::monoclinic ? 4 : 2;
    free_dims = sys == Sys::cubic ? std::vector<int>{0} : sys == Sys::orthorhombic ? std::vector<int>{0, 1, 2}
                : sys == Sys::monoclinic ? std::vector<int>{0, 1, 2, 3} : std::vector<int>{0, 2};
  }

  void add_indices() {
    H = int(o.max_length * std::sqrt(qmax)) + 1;
    lookup.assign(size_t(2 * H + 1) * size_t(H + 1) * size_t(H + 1), -1);
    for (int h = (sys == Sys::monoclinic ? -H : 0); h <= H; ++h)
      for (int k = 0; k <= H; ++k)
        for (int l = 0; l <= H; ++l) {
          if (!h && !k && !l) continue;
          if (sys == Sys::cubic && !(h >= k && k >= l)) continue;
          if ((sys == Sys::tetragonal || sys == Sys::hexagonal) && h < k) continue;
          // the monoclinic h < 0 half with l = 0 repeats h > 0 (Q is even in h when l = 0)
          if (sys == Sys::monoclinic && h < 0 && l == 0) continue;
          const double h2 = double(h) * h, k2 = double(k) * k, l2 = double(l) * l;
          std::array<double, 4> c{0, 0, 0, 0};
          switch (sys) {
            case Sys::cubic: c = {h2 + k2 + l2, 0, 0, 0}; break;
            case Sys::tetragonal: c = {h2 + k2, l2, 0, 0}; break;
            case Sys::hexagonal: c = {4.0 / 3.0 * (h2 + double(h) * k + k2), l2, 0, 0}; break;
            case Sys::orthorhombic: c = {h2, k2, l2, 0}; break;
            case Sys::monoclinic: c = {h2, k2, l2, double(h) * l}; break;
          }
          lookup[(size_t(h + H) * size_t(H + 1) + size_t(k)) * size_t(H + 1) + size_t(l)] = int(hkl.size());
          hkl.push_back({h, k, l});
          co.push_back(c);
        }
  }
  // the indices that can fall below qmax somewhere in a start box: Q ≥ (h/a)² + (k/b)² + (l/c)² bounded by each index
  // alone (for the monoclinic h·l term too: Q over β is never below max(|h|/a, |l|/c)²)
  std::vector<int> start_candidates(const Box& b) const {
    const double s = std::sqrt(qmax);
    const double A = b.hi[0], B = sys == Sys::orthorhombic || sys == Sys::monoclinic ? b.hi[1] : b.hi[0],
                 C = sys == Sys::cubic ? b.hi[0] : b.hi[2];
    // the hexagonal metric 4/3 (h² + hk + k²) ≥ h² and ≥ k² (h, k ≥ 0)
    const int Ha = std::min(H, int(A * s) + 1), Kb = std::min(H, int(B * s) + 1), Lc = std::min(H, int(C * s) + 1);
    std::vector<int> out;
    for (int h = (sys == Sys::monoclinic ? -Ha : 0); h <= Ha; ++h)
      for (int k = 0; k <= Kb; ++k)
        for (int l = 0; l <= Lc; ++l) {
          const int id = lookup[(size_t(h + H) * size_t(H + 1) + size_t(k)) * size_t(H + 1) + size_t(l)];
          if (id >= 0) out.push_back(id);
        }
    return out;
  }

  // the box's constants: reciprocal lengths at the box's ends, |cos β| and 1/sin² β at its ends
  struct BoxK { double ilo[3], ihi[3], g0, g1, w0, w1; };
  BoxK consts(const Box& b) const {
    BoxK k{};
    for (int d = 0; d < 3; ++d) k.ilo[d] = 1 / b.hi[d], k.ihi[d] = 1 / b.lo[d];
    if (sys == Sys::monoclinic) {
      // β in [90°, 180°): |cos β| grows with β
      k.g0 = -std::cos(b.lo[3] / kDeg), k.g1 = -std::cos(b.hi[3] / kDeg);
      k.w0 = 1 / (1 - k.g0 * k.g0), k.w1 = 1 / (1 - k.g1 * k.g1);
    }
    return k;
  }
  // the exact Q range of index id over the box
  std::pair<double, double> qrange(int id, const BoxK& k) const {
    const auto& c = co[size_t(id)];
    switch (sys) {
      case Sys::cubic: return {c[0] * k.ilo[0] * k.ilo[0], c[0] * k.ihi[0] * k.ihi[0]};
      case Sys::tetragonal:
      case Sys::hexagonal:
        return {c[0] * k.ilo[0] * k.ilo[0] + c[1] * k.ilo[2] * k.ilo[2], c[0] * k.ihi[0] * k.ihi[0] + c[1] * k.ihi[2] * k.ihi[2]};
      case Sys::orthorhombic:
        return {c[0] * k.ilo[0] * k.ilo[0] + c[1] * k.ilo[1] * k.ilo[1] + c[2] * k.ilo[2] * k.ilo[2],
                c[0] * k.ihi[0] * k.ihi[0] + c[1] * k.ihi[1] * k.ihi[1] + c[2] * k.ihi[2] * k.ihi[2]};
      case Sys::monoclinic: break;
    }
    // monoclinic (b unique, β ≥ 90°): Q = k²/b² + (u² + v² + 2 s u v γ)/(1 − γ²), u = |h|/a, v = |l|/c, γ = |cos β|,
    // s the sign of h·l. With s ≥ 0 Q grows with u, v and γ. With s < 0 Q is convex in (u, v) and, in γ, falls to
    // max(u, v)² (at γ = min(u, v)/max(u, v)) and rises again: its maximum is at a corner, its minimum on a γ face's
    // edges or, when that ratio can lie inside the box, no lower than max(u0, v0)².
    const double kb_lo = c[1] * k.ilo[1] * k.ilo[1], kb_hi = c[1] * k.ihi[1] * k.ihi[1];
    const double ah = std::sqrt(c[0]), al = std::sqrt(c[2]);
    const double u0 = ah * k.ilo[0], u1 = ah * k.ihi[0], v0 = al * k.ilo[2], v1 = al * k.ihi[2];
    if (c[3] >= 0) return {kb_lo + (u0 * u0 + v0 * v0 + 2 * u0 * v0 * k.g0) * k.w0, kb_hi + (u1 * u1 + v1 * v1 + 2 * u1 * v1 * k.g1) * k.w1};
    auto g = [](double u, double v, double ga) { return u * u + v * v - 2 * u * v * ga; };
    double hi = 0, lo = 1e300;
    const double gs[2] = {k.g0, k.g1}, ws[2] = {k.w0, k.w1};
    for (int e = 0; e < 2; ++e) {
      const double ga = gs[e], w = ws[e];
      hi = std::max({hi, g(u0, v0, ga) * w, g(u0, v1, ga) * w, g(u1, v0, ga) * w, g(u1, v1, ga) * w});
      lo = std::min({lo, g(u0, std::clamp(u0 * ga, v0, v1), ga) * w, g(u1, std::clamp(u1 * ga, v0, v1), ga) * w,
                     g(std::clamp(v0 * ga, u0, u1), v0, ga) * w, g(std::clamp(v1 * ga, u0, u1), v1, ga) * w});
    }
    const double t0 = u0 / v1, t1 = u1 / v0;   // the range of u/v
    double rmin, rmax;
    if (t0 <= 1 && t1 >= 1) rmax = 1, rmin = std::min(t0, 1 / t1);
    else if (t1 < 1) rmin = t0, rmax = t1;
    else rmin = 1 / t1, rmax = 1 / t0;
    if (rmax >= k.g0 && rmin <= k.g1) lo = std::min(lo, std::max(u0, v0) * std::max(u0, v0));
    return {kb_lo + lo, kb_hi + hi};
  }
  // every line within its tolerance of some index's Q range over the box (up to o.unindexed lines not); keep: the
  // candidates that can still fall below qmax inside the box
  bool passes(const Box& b, const std::vector<int>& cand, std::vector<int>& keep) {
    ++tests;
    const BoxK k = consts(b);
    rr.clear(), keep.clear();
    for (int id : cand) {
      const auto r = qrange(id, k);
      if (r.first <= qmax) rr.push_back(r), keep.push_back(id);
    }
    int miss = 0;
    for (size_t i = 0; i < Q.size(); ++i) {
      bool ok = false;
      for (const auto& r : rr)
        if (Q[i] >= r.first - tol[i] && Q[i] <= r.second + tol[i]) { ok = true; break; }
      if (!ok && ++miss > o.unindexed) return false;
    }
    return true;
  }
  double vol_min(const Box& b) const {
    switch (sys) {
      case Sys::cubic: return b.lo[0] * b.lo[0] * b.lo[0];
      case Sys::tetragonal: return b.lo[0] * b.lo[0] * b.lo[2];
      case Sys::hexagonal: return std::sqrt(3.0) / 2 * b.lo[0] * b.lo[0] * b.lo[2];
      case Sys::orthorhombic: return b.lo[0] * b.lo[1] * b.lo[2];
      case Sys::monoclinic: break;
    }
    return b.lo[0] * b.lo[1] * b.lo[2] * std::min(std::sin(b.lo[3] / kDeg), std::sin(b.hi[3] / kDeg));
  }
  // how far the top line's Q can move across the box along each free dimension (an estimate for the split and the leaf)
  void spreads(const Box& b, double s[4]) const {
    const double qt = Q.back();
    for (int d = 0; d < 3; ++d) {
      const double r = b.hi[d] / b.lo[d];
      s[d] = qt * (r * r - 1);
    }
    s[3] = 0;
    if (sys == Sys::monoclinic) {
      // |dQ/dβ| / Q is at most about 1/(1 − γ) + 2γ/(1 − γ²) (an h·l reflection and an h00 one), γ = |cos β|
      const double g1 = -std::cos(b.hi[3] / kDeg);
      s[3] = qt * (b.hi[3] - b.lo[3]) / kDeg * (1 / (1 - g1) + 2 * g1 / (1 - g1 * g1));
    }
  }
  // least squares on Q with the nearest indices; the cell, M20 and F_N
  void refine_leaf(const Box& b, const std::vector<int>& cand) {
    if (cand.empty()) return;
    double p[4];
    for (int k = 0; k < 4; ++k) p[k] = 0.5 * (b.lo[k] + b.hi[k]);
    double x[4] = {0, 0, 0, 0};   // the reciprocal metric parameters
    const double sb = std::sin(p[3] / kDeg), cb = std::cos(p[3] / kDeg);
    switch (sys) {
      case Sys::cubic: x[0] = 1 / (p[0] * p[0]); break;
      case Sys::tetragonal:
      case Sys::hexagonal: x[0] = 1 / (p[0] * p[0]), x[1] = 1 / (p[2] * p[2]); break;
      case Sys::orthorhombic: x[0] = 1 / (p[0] * p[0]), x[1] = 1 / (p[1] * p[1]), x[2] = 1 / (p[2] * p[2]); break;
      case Sys::monoclinic:
        x[0] = 1 / (p[0] * p[0] * sb * sb), x[1] = 1 / (p[1] * p[1]), x[2] = 1 / (p[2] * p[2] * sb * sb), x[3] = -2 * cb / (p[0] * p[2] * sb * sb);
        break;
    }
    const size_t mm = size_t(m);
    auto qcalc = [&](int id) {
      const auto& c = co[size_t(id)];
      return c[0] * x[0] + c[1] * x[1] + c[2] * x[2] + c[3] * x[3];
    };
    std::vector<int> as(Q.size(), -1);
    std::vector<char> ok(Q.size());
    for (int it = 0; it < 4; ++it) {
      for (size_t i = 0; i < Q.size(); ++i) {
        double best = 1e300;
        for (int id : cand) {
          const double d = std::fabs(qcalc(id) - Q[i]);
          if (d < best) best = d, as[i] = id;
        }
        // the box's centre can sit a few tolerances off: the first passes take every line's nearest reflection
        ok[i] = it < 2 || best <= tol[i];
      }
      double A[16] = {0}, rhs[4] = {0};
      for (size_t i = 0; i < Q.size(); ++i) {
        if (!ok[i] || as[i] < 0) continue;
        const auto& c = co[size_t(as[i])];
        for (size_t a = 0; a < mm; ++a) {
          rhs[a] += c[a] * Q[i];
          for (size_t bb = 0; bb < mm; ++bb) A[a * mm + bb] += c[a] * c[bb];
        }
      }
      // solve (Gauss–Jordan)
      for (size_t c = 0; c < mm; ++c) {
        size_t piv = c;
        for (size_t i = c + 1; i < mm; ++i) if (std::fabs(A[i * mm + c]) > std::fabs(A[piv * mm + c])) piv = i;
        if (std::fabs(A[piv * mm + c]) < 1e-300) return;
        for (size_t j = 0; j < mm; ++j) std::swap(A[c * mm + j], A[piv * mm + j]);
        std::swap(rhs[c], rhs[piv]);
        const double d = A[c * mm + c];
        for (size_t j = 0; j < mm; ++j) A[c * mm + j] /= d;
        rhs[c] /= d;
        for (size_t i = 0; i < mm; ++i)
          if (i != c) {
            const double f = A[i * mm + c];
            for (size_t j = 0; j < mm; ++j) A[i * mm + j] -= f * A[c * mm + j];
            rhs[i] -= f * rhs[c];
          }
      }
      for (size_t k = 0; k < mm; ++k) x[k] = rhs[k];
    }
    for (size_t k = 0; k < mm; ++k) if (!std::isfinite(x[k])) return;
    // only the monoclinic h·l term x[3] may be negative
    for (size_t k = 0; k < std::min<size_t>(mm, 3); ++k) if (x[k] <= 0) return;
    IndexSolution s;
    s.system = name;
    switch (sys) {
      case Sys::cubic: s.cell = {1 / std::sqrt(x[0]), 1 / std::sqrt(x[0]), 1 / std::sqrt(x[0]), 90, 90, 90}; break;
      case Sys::tetragonal: s.cell = {1 / std::sqrt(x[0]), 1 / std::sqrt(x[0]), 1 / std::sqrt(x[1]), 90, 90, 90}; break;
      case Sys::hexagonal: s.cell = {1 / std::sqrt(x[0]), 1 / std::sqrt(x[0]), 1 / std::sqrt(x[1]), 90, 90, 120}; break;
      case Sys::orthorhombic: s.cell = {1 / std::sqrt(x[0]), 1 / std::sqrt(x[1]), 1 / std::sqrt(x[2]), 90, 90, 90}; break;
      case Sys::monoclinic: {
        const double cosb = -x[3] / (2 * std::sqrt(x[0] * x[2]));
        if (std::fabs(cosb) >= 1) return;
        const double beta = std::acos(cosb) * kDeg, sn = std::sin(beta / kDeg);
        s.cell = {1 / (std::sqrt(x[0]) * sn), 1 / std::sqrt(x[1]), 1 / (std::sqrt(x[2]) * sn), 90, beta, 90};
        break;
      }
    }
    const auto& c = s.cell;
    const double ca = std::cos(c[3] / kDeg), cbb = std::cos(c[4] / kDeg), cg = std::cos(c[5] / kDeg);
    s.volume = c[0] * c[1] * c[2] * std::sqrt(std::max(0.0, 1 - ca * ca - cbb * cbb - cg * cg + 2 * ca * cbb * cg));
    if (s.volume > o.max_volume || c[0] > o.max_length * 1.05 || c[1] > o.max_length * 1.05 || c[2] > o.max_length * 1.05) return;
    // figures of merit
    double eps = 0, d2t = 0;
    int nind = 0;
    for (size_t i = 0; i < Q.size(); ++i) {
      double best = 1e300;
      int bh = -1;
      for (int id : cand) {
        const double d = std::fabs(qcalc(id) - Q[i]);
        if (d < best) best = d, bh = id;
      }
      const bool in = bh >= 0 && best <= tol[i];
      s.hkl.push_back(in ? hkl[size_t(bh)] : I3{0, 0, 0});
      const double sn2 = std::sqrt(std::max(0.0, Q[i])) * o.wavelength / 2;
      s.two_theta_obs.push_back(2 * std::asin(std::min(1.0, sn2)) * kDeg);
      const double qc = bh >= 0 ? qcalc(bh) : 0, snc = std::sqrt(std::max(0.0, qc)) * o.wavelength / 2;
      s.two_theta_calc.push_back(2 * std::asin(std::min(1.0, snc)) * kDeg);
      // each line's error floored at a tenth of its tolerance (what a peak search can place): exact data would otherwise
      // give any cell that indexes the lines an unbounded merit
      if (in) ++nind, eps += std::max(best, 0.1 * tol[i]), d2t += std::max(std::fabs(s.two_theta_obs.back() - s.two_theta_calc.back()), 0.1 * o.tolerance);
    }
    if (nind < int(Q.size()) - o.unindexed || nind == 0) return;
    s.indexed = nind;
    // N: the distinct calculated lines up to the last observed one, over every index (the refined cell may lie outside the box)
    std::set<long> distinct;
    for (size_t id = 0; id < hkl.size(); ++id) {
      const double q = qcalc(int(id));
      if (q > 0 && q <= Q.back() * (1 + 1e-9)) distinct.insert(std::lround(q * 1e7));
    }
    const double N = double(std::max<size_t>(1, distinct.size()));
    s.M20 = Q.back() / (2 * (eps / nind) * N);
    s.FN = (double(Q.size()) / N) / std::max(1e-9, d2t / nind);
    // a cell already found (within 0.2 %) keeps the better merit
    for (auto& f : found)
      if (f.system == s.system && std::fabs(f.volume - s.volume) < 2e-3 * s.volume) {
        bool same = true;
        for (int k = 0; k < 6; ++k) same = same && std::fabs(f.cell[size_t(k)] - s.cell[size_t(k)]) < 2e-3 * std::max(1.0, f.cell[size_t(k)]);
        if (same) {
          if (s.M20 > f.M20) f = s;
          return;
        }
      }
    found.push_back(s);
  }
  void poll() {
    if (o.progress && !o.progress(name, frac)) throw std::runtime_error("index stopped");
    if (o.max_boxes > 0 && tests >= o.max_boxes) stopped = true;
  }
  void search(const Box& b, int depth, const std::vector<int>& cand) {
    if (stopped) return;
    if ((tests & 0xFFFF) == 0xFFFF) {
      poll();
      if (stopped) return;
    }
    if (vol_min(b) > o.max_volume) return;
    if (vol_limit > 0 && vol_min(b) > vol_limit) return;
    // monoclinic: only the reduced setting of the ac net, a ≤ c and c |cos β| ≤ a/2 (every lattice has one, with
    // β ≤ 120°); the other settings are the same lattices again
    if (sys == Sys::monoclinic && (b.lo[0] > b.hi[2] || -std::cos(b.lo[3] / kDeg) * b.lo[2] > 0.5 * b.hi[0])) return;
    std::vector<int> keep;
    if (!passes(b, cand, keep)) return;
    double s[4];
    spreads(b, s);
    int dsplit = -1;
    double smax = 0;
    for (int d : free_dims)
      if (s[d] > smax) smax = s[d], dsplit = d;
    if (dsplit < 0 || smax <= kLeafTolerances * tol.back() || depth > 80) {
      const size_t before = found.size();
      refine_leaf(b, keep);
      if (found.size() > before && found.back().M20 >= 5) {
        const double v = found.back().volume * 2.05;
        vol_limit = vol_limit > 0 ? std::min(vol_limit, v) : v;
      }
      return;
    }
    // halve the dimension that moves Q the most
    Box l = b, u = b;
    const double mid = 0.5 * (b.lo[dsplit] + b.hi[dsplit]);
    l.hi[dsplit] = mid, u.lo[dsplit] = mid;
    search(l, depth + 1, keep);
    search(u, depth + 1, keep);
  }
};
}  // namespace

IndexResult index_pattern(const std::vector<double>& two_theta_in, const IndexOptions& o) {
  const auto t0 = std::chrono::steady_clock::now();
  IndexResult res;
  std::vector<double> tt = two_theta_in;
  std::sort(tt.begin(), tt.end());
  if (tt.size() > size_t(o.lines)) tt.resize(size_t(o.lines));
  if (tt.size() < 4) throw std::invalid_argument("index: four lines or more");
  if (!(o.max_length > 2)) throw std::invalid_argument("index: the longest cell edge must be above 2 Å");
  for (const auto& sys : o.systems) (void)system_of(sys);
  res.lines = tt;
  // DICVOL's convention: a cell of higher symmetry that indexes every line bounds the volume searched in the lower
  // systems that follow (beyond about twice its volume they only hold its supercells)
  // and a cell that indexes every line with a high merit ends the search of systems with more free parameters (they
  // would add only its supercells and pseudo-symmetric cells, and with few lines a family of cells that all fit)
  double vol_cap = 0;
  std::string cap_note, strong_note;
  int strong_m = 0;
  for (const auto& sys : o.systems) {
    if (o.progress && !o.progress(sys, 0)) throw std::runtime_error("index stopped");
    Indexer ix(o, sys);
    if (strong_m > 0 && ix.m > strong_m) {
      res.notes.push_back(sys + ": not searched, " + strong_note + " (--systems " + sys + " alone searches it)");
      continue;
    }
    for (double t : tt) {
      const double th = t / 2 / kDeg, s = std::sin(th);
      ix.Q.push_back(4 * s * s / (o.wavelength * o.wavelength));
      ix.tol.push_back(4 / (o.wavelength * o.wavelength) * s * std::cos(th) * (o.tolerance / kDeg));
    }
    ix.qmax = ix.Q.back() + ix.tol.back();
    ix.add_indices();
    if (vol_cap > 0) ix.vol_limit = vol_cap;
    const double lmin = 2.0, step = 0.5;
    const double bmin = std::max(90.0, o.min_beta), bmax = std::min(179.0, o.max_beta);
    std::vector<Box> start;
    for (double a = lmin; a < o.max_length; a += step) {
      if (ix.sys == Sys::cubic) { start.push_back({{a, a, a, 90}, {a + step, a + step, a + step, 90}}); continue; }
      for (double c = lmin; c < o.max_length; c += step) {
        if (ix.sys == Sys::tetragonal || ix.sys == Sys::hexagonal) { start.push_back({{a, a, c, 90}, {a + step, a + step, c + step, 90}}); continue; }
        for (double b = (ix.sys == Sys::orthorhombic ? a : lmin); b < (ix.sys == Sys::orthorhombic ? c + step : o.max_length); b += step) {
          if (ix.sys == Sys::orthorhombic) {
            if (b > c) continue;
            start.push_back({{a, b, c, 90}, {a + step, b + step, c + step, 90}});
          } else {
            if (a > c) continue;   // a ≤ c for the monoclinic setting (the other order is the same lattice)
            for (double be = bmin; be < bmax; be += 5) start.push_back({{a, b, c, be}, {a + step, b + step, c + step, std::min(be + 5, bmax)}});
          }
        }
      }
    }
    // the smallest cells first (the indexing sought is the smallest cell that indexes the lines)
    std::sort(start.begin(), start.end(), [&](const Box& x, const Box& y) { return ix.vol_min(x) < ix.vol_min(y); });
    double reached = 0;
    for (size_t k = 0; k < start.size() && !ix.stopped; ++k) {
      if (ix.vol_limit > 0 && ix.vol_min(start[k]) > ix.vol_limit) break;
      if (ix.vol_min(start[k]) > o.max_volume) break;
      ix.frac = double(k) / double(start.size());
      if (k % 200 == 0) ix.poll();
      if (ix.stopped) break;
      reached = ix.vol_min(start[k]);
      ix.search(start[k], 0, ix.start_candidates(start[k]));
    }
    for (auto& f : ix.found) res.solutions.push_back(f);
    char b[320];
    std::snprintf(b, sizeof b, "%s: %zu cells (%ld boxes tested)", sys.c_str(), ix.found.size(), ix.tests);
    res.notes.push_back(b);
    if (vol_cap > 0) {
      std::snprintf(b, sizeof b, "%s: cells above %.1f Å³ not searched (%s)", sys.c_str(), vol_cap, cap_note.c_str());
      res.notes.push_back(b);
    }
    if (ix.stopped) {
      std::snprintf(b, sizeof b, "%s: search stopped after %ld boxes (the budget), cells up to about %.0f Å³ searched; give fewer systems, "
                    "a smaller --max-length or more lines", sys.c_str(), ix.tests, reached);
      res.notes.push_back(b);
    }
    for (const auto& f : ix.found)
      if (f.indexed == int(tt.size()) && f.M20 >= kStrongM20 && strong_m == 0) {
        strong_m = ix.m;
        std::snprintf(b, sizeof b, "the %s cell %.4f %.4f %.4f indexes every line with M20 %.1f", f.system.c_str(), f.cell[0], f.cell[1], f.cell[2], f.M20);
        strong_note = b;
      }
    for (const auto& f : ix.found)
      if (f.indexed == int(tt.size()) && f.M20 >= 10 && (vol_cap <= 0 || 2.05 * f.volume < vol_cap)) {
        vol_cap = 2.05 * f.volume;
        std::snprintf(b, sizeof b, "about twice the volume of the %s cell %.4f %.4f %.4f that indexes every line, M20 %.1f", f.system.c_str(),
                      f.cell[0], f.cell[1], f.cell[2], f.M20);
        cap_note = b;
      }
  }
  // ranked by lines indexed, then by M20 with the higher symmetry preferred when the merits are close (a factor of two per
  // step from monoclinic to cubic): a cell of lower symmetry with a free parameter more can always fit as well
  auto rank = [](const std::string& s) { return s == "cubic" ? 4 : s == "hexagonal" ? 3 : s == "tetragonal" ? 2 : s == "orthorhombic" ? 1 : 0; };
  std::sort(res.solutions.begin(), res.solutions.end(), [&](const IndexSolution& a, const IndexSolution& b) {
    if (a.indexed != b.indexed) return a.indexed > b.indexed;
    return a.M20 * std::pow(2.0, rank(a.system)) > b.M20 * std::pow(2.0, rank(b.system));
  });
  if (res.solutions.size() > size_t(o.solutions)) res.solutions.resize(size_t(o.solutions));
  res.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  return res;
}

// ---------------------------------------------------------------- direct-space solution
namespace {
struct Torsion { size_t i, j; std::vector<size_t> moving; };

std::vector<Torsion> torsions_of(const System& m) {
  const size_t n = m.atoms.size();
  std::vector<std::vector<size_t>> nb(n);
  for (const auto& b : m.bonds) nb[b.i].push_back(b.j), nb[b.j].push_back(b.i);
  std::vector<Torsion> out;
  for (const auto& b : m.bonds) {
    const size_t i = b.i, j = b.j;
    if (m.atoms[i].element == 1 || m.atoms[j].element == 1) continue;
    int hi = 0, hj = 0;
    for (size_t k : nb[i]) hi += k != j && m.atoms[k].element != 1;
    for (size_t k : nb[j]) hj += k != i && m.atoms[k].element != 1;
    if (!hi || !hj) continue;
    // the side of j without the bond; a ring bond reaches i again
    std::vector<char> seen(n, 0);
    std::vector<size_t> st{j};
    seen[j] = 1;
    bool ring = false;
    while (!st.empty()) {
      const size_t v = st.back();
      st.pop_back();
      for (size_t w : nb[v]) {
        if (v == j && w == i) continue;
        if (w == i) { ring = true; break; }
        if (!seen[w]) seen[w] = 1, st.push_back(w);
      }
      if (ring) break;
    }
    if (ring) continue;
    Torsion t{i, j, {}};
    for (size_t k = 0; k < n; ++k) if (seen[k]) t.moving.push_back(k);
    if (t.moving.size() * 2 > n) {   // move the smaller side
      std::vector<size_t> other;
      for (size_t k = 0; k < n; ++k) if (!seen[k]) other.push_back(k);
      t = Torsion{j, i, other};
    }
    out.push_back(t);
  }
  return out;
}

Vec3 rotate(const Vec3& v, const Vec3& axis, double ang) {
  const double c = std::cos(ang), s = std::sin(ang);
  return v * c + cross(axis, v) * s + axis * (dot(axis, v) * (1 - c));
}
}  // namespace

SolveResult solve_structure(const MeasuredPattern& data, const System& mol, const std::array<double, 6>& cp, const SolveOptions& o) {
  SolveResult res;
  if (mol.atoms.empty()) throw std::invalid_argument("solve: no molecule");
  if (o.steps < 1) throw std::invalid_argument("solve: steps must be 1 or more");
  const SpaceGroupSetting* g = find_space_group(o.space_group);
  if (!g) throw std::invalid_argument("solve: space group '" + o.space_group + "' not known");
  const auto ops = hall_operations(g->hall);
  const Cell cell = cell_parameters(cp[0], cp[1], cp[2], cp[3], cp[4], cp[5]);
  const size_t n = mol.atoms.size();
  // the molecule about its centroid
  Vec3 c0{0, 0, 0};
  for (const auto& a : mol.atoms) c0 = c0 + a.pos;
  c0 = c0 * (1.0 / double(n));
  std::vector<Vec3> base(n);
  for (size_t i = 0; i < n; ++i) base[i] = mol.atoms[i].pos - c0;
  const std::vector<Torsion> tors = o.torsions ? torsions_of(mol) : std::vector<Torsion>{};
  res.torsions = int(tors.size());
  std::vector<size_t> heavy;
  for (size_t i = 0; i < n; ++i) if (mol.atoms[i].element != 1) heavy.push_back(i);
  // the contact distance as a fraction of each lattice plane spacing: how many cells away an image can still touch
  const double vol = cell.volume();
  const double reach[3] = {o.contact * norm(cross(cell.b, cell.c)) / vol, o.contact * norm(cross(cell.c, cell.a)) / vol,
                           o.contact * norm(cross(cell.a, cell.b)) / vol};
  struct St { Vec3 f; Vec3 axis; double ang; std::vector<double> tor; };
  auto geometry = [&](const St& s) {
    std::vector<Vec3> x = base;
    for (size_t k = 0; k < tors.size(); ++k) {
      const auto& t = tors[k];
      Vec3 ax = x[t.j] - x[t.i];
      ax = ax * (1 / norm(ax));
      for (size_t m : t.moving) x[m] = x[t.j] + rotate(x[m] - x[t.j], ax, s.tor[k]);
    }
    Vec3 cc{0, 0, 0};
    for (const auto& v : x) cc = cc + v;
    cc = cc * (1.0 / double(n));
    const Vec3 ctr = cell.to_cartesian(s.f);
    for (auto& v : x) v = ctr + rotate(v - cc, s.axis, s.ang);
    return x;
  };
  auto crystal = [&](const std::vector<Vec3>& x, bool* clash) {
    System cr;
    cr.cell = cell;
    std::vector<std::vector<Vec3>> copies, fracs;
    for (const auto& op : ops) {
      std::vector<Vec3> f(n);
      Vec3 fc{0, 0, 0};
      for (size_t i = 0; i < n; ++i) {
        const Vec3 fr = cell.to_fractional(x[i]);
        for (int r = 0; r < 3; ++r) f[i][size_t(r)] = op.R[r][0] * fr[0] + op.R[r][1] * fr[1] + op.R[r][2] * fr[2] + op.t[r];
        fc = fc + f[i];
      }
      fc = fc * (1.0 / double(n));
      const Vec3 sh{std::floor(fc[0]), std::floor(fc[1]), std::floor(fc[2])};
      std::vector<Vec3> cart(n);
      for (size_t i = 0; i < n; ++i) f[i] = f[i] - sh, cart[i] = cell.to_cartesian(f[i]);
      copies.push_back(cart);
      fracs.push_back(f);
    }
    if (clash) {
      // every pair of heavy atoms, over every lattice image within the contact distance: those of two copies, and a
      // copy's atoms against its own images in the neighbouring cells (a short axis can bring a molecule onto itself)
      *clash = false;
      const double c2 = o.contact * o.contact;
      for (size_t a = 0; a < fracs.size() && !*clash; ++a)
        for (size_t b = a; b < fracs.size() && !*clash; ++b)
          for (size_t ii = 0; ii < heavy.size() && !*clash; ++ii)
            for (size_t jj = (a == b ? ii : 0); jj < heavy.size() && !*clash; ++jj) {
              Vec3 d = fracs[b][heavy[jj]] - fracs[a][heavy[ii]];
              for (int k = 0; k < 3; ++k) d[size_t(k)] -= std::round(d[size_t(k)]);
              int lo[3], hi[3];
              for (int k = 0; k < 3; ++k) lo[k] = int(std::ceil(-reach[k] - d[size_t(k)])), hi[k] = int(std::floor(reach[k] - d[size_t(k)]));
              for (int l0 = lo[0]; l0 <= hi[0] && !*clash; ++l0)
                for (int l1 = lo[1]; l1 <= hi[1] && !*clash; ++l1)
                  for (int l2 = lo[2]; l2 <= hi[2]; ++l2) {
                    if (a == b && !l0 && !l1 && !l2) continue;   // the molecule's own bonds and contacts
                    const Vec3 v = cell.a * (d[0] + l0) + cell.b * (d[1] + l1) + cell.c * (d[2] + l2);
                    if (dot(v, v) < c2) { *clash = true; break; }
                  }
            }
    }
    for (const auto& cp2 : copies)
      for (size_t i = 0; i < n; ++i) {
        Atom at = mol.atoms[i];
        at.pos = cp2[i];
        cr.atoms.push_back(at);
      }
    return cr;
  };
  // the peak profile: the trials are scored with a fixed profile, so it is fitted to the data first (Pawley, the given
  // cell held); a profile far from the data's moves the minimum of Rwp away from the true structure
  RefineOptions ro = o.refine;
  ro.progress = nullptr;
  if (o.fit_profile) {
    try {
      RefineOptions po = ro;
      po.mode = "pawley";
      RefinePhase ph;
      // any placement gives the cell and the space group's absences
      ph.cell = crystal(geometry(St{{0.13, 0.27, 0.41}, {0, 0, 1}, 0.0, std::vector<double>(tors.size(), 0.0)}), nullptr), ph.lattice = false;
      if (o.progress && !o.progress("profile", 0)) throw std::runtime_error("solve stopped");
      const RefineResult pr = refine_pattern(data, {ph}, po);
      if (std::isfinite(pr.U) && std::isfinite(pr.V) && std::isfinite(pr.W) && std::isfinite(pr.X) && std::isfinite(pr.Y) && std::isfinite(pr.zero)) {
        ro.U = pr.U, ro.V = pr.V, ro.W = pr.W, ro.X = pr.X, ro.Y = pr.Y, ro.zero = pr.zero;
        char b[240];
        std::snprintf(b, sizeof b, "peak profile from a Pawley fit in the cell (Rwp %.4f): U %.6f V %.6f W %.6f X %.6f Y %.6f, zero %.4f°", pr.rwp, ro.U,
                      ro.V, ro.W, ro.X, ro.Y, ro.zero);
        res.notes.push_back(b);
      } else {
        res.notes.push_back("the Pawley fit of the peak profile did not converge: the given profile used");
      }
    } catch (const std::runtime_error& e) {
      if (std::string(e.what()) == "solve stopped") throw;
      res.notes.push_back(std::string("the peak profile could not be fitted (") + e.what() + "): the given profile used");
    } catch (const std::exception& e) {
      res.notes.push_back(std::string("the peak profile could not be fitted (") + e.what() + "): the given profile used");
    }
  }
  res.profile = {ro.U, ro.V, ro.W, ro.X, ro.Y, ro.zero};
  ro.refine_U = ro.refine_V = ro.refine_W = ro.refine_X = ro.refine_Y = ro.refine_zero = false;
  auto rwp_of = [&](const System& cr) {
    RefinePhase ph;
    ph.cell = cr, ph.lattice = false;
    return refine_pattern(data, {ph}, ro).rwp;
  };
  std::mt19937_64 rng(o.seed);
  std::uniform_real_distribution<double> U(0.0, 1.0);
  std::normal_distribution<double> N(0.0, 1.0);
  auto rand_axis = [&] {
    Vec3 a{N(rng), N(rng), N(rng)};
    return a * (1 / std::max(1e-12, norm(a)));
  };
  auto rand_state = [&] {
    St s{{U(rng), U(rng), U(rng)}, rand_axis(), 2 * kPi * U(rng), std::vector<double>(tors.size(), 0.0)};
    for (auto& t : s.tor) t = 2 * kPi * U(rng);
    return s;
  };
  // the best clash-free state (only such a state is ever a result)
  bool found = false;
  double best = 0;
  St best_s;
  auto consider = [&](double v, const St& s) {
    if (!found || v < best || (std::isnan(best) && !std::isnan(v))) best = v, best_s = s, found = true;
  };
  // --steps counts every move: shared among the restarts
  const int restarts = std::max(1, std::min(o.restarts, o.steps));
  const int total = o.steps;
  int done = 0;
  for (int rs = 0; rs < restarts; ++rs) {
    const int steps = total / restarts + (rs < total % restarts ? 1 : 0);
    // a random start free of contacts (a few hundred tries; a crowded cell may have none, and the moves then look for one)
    St s = rand_state();
    bool clash = false;
    System cr0 = crystal(geometry(s), &clash);
    for (int tries = 0; clash && tries < 300; ++tries) s = rand_state(), cr0 = crystal(geometry(s), &clash);
    double cur = 1e3;   // a clashing start: any free state below is taken
    if (!clash) cur = rwp_of(cr0), consider(cur, s);
    for (int it = 0; it < steps; ++it, ++done) {
      if (o.progress && it % 50 == 0 && !o.progress("annealing", double(done) / total)) throw std::runtime_error("solve stopped");
      // the last fifth: downhill only, with small moves (a polish of the best region)
      const bool polish = it >= steps * 4 / 5;
      const double T = polish ? 0 : o.t_start * std::pow(o.t_end / o.t_start, double(it) / std::max(1, steps * 4 / 5 - 1));
      // moves shrink as it cools
      const double sc = polish ? 0.03 : std::max(0.05, std::sqrt(T / o.t_start));
      St t = s;
      const int kind = int(U(rng) * (tors.empty() ? 2 : 3));
      if (kind == 0) {
        for (int k = 0; k < 3; ++k) t.f[size_t(k)] += N(rng) * 0.5 * sc / cp[size_t(k)];
      } else if (kind == 1) {
        // compose a small rotation with the current one (axis-angle by quaternions)
        const Vec3 ax = rand_axis();
        const double da = N(rng) * 0.4 * sc;
        const double qw1 = std::cos(t.ang / 2), qw2 = std::cos(da / 2);
        const Vec3 qv1 = t.axis * std::sin(t.ang / 2), qv2 = ax * std::sin(da / 2);
        const double w = qw2 * qw1 - dot(qv2, qv1);
        const Vec3 v = qv1 * qw2 + qv2 * qw1 + cross(qv2, qv1);
        const double nv = norm(v);
        t.ang = 2 * std::atan2(nv, w);
        t.axis = nv > 1e-12 ? v * (1 / nv) : Vec3{0, 0, 1};
      } else {
        t.tor[size_t(U(rng) * double(tors.size())) % tors.size()] += N(rng) * 0.8 * sc;
      }
      bool cl = false;
      const System cr = crystal(geometry(t), &cl);
      if (cl) continue;
      const double v = rwp_of(cr);
      if (v < cur || (T > 0 && U(rng) < std::exp(-(v * v - cur * cur) / T))) {
        s = t, cur = v, ++res.accepted;
        consider(v, t);
      }
      if (it % 100 == 0) res.rwp_history.push_back(cur);
    }
  }
  if (!found) {
    char b[400];
    std::snprintf(b, sizeof b,
                  "solve: no placement of the molecule keeps the heavy atoms of its %zu symmetry copies in %s (and of their lattice images) "
                  "%.2f Å apart in this cell (V = %.1f Å³, %.1f Å³ per molecule) within %d moves; the cell is likely too small for %zu molecules, "
                  "or give more --steps or another --seed",
                  ops.size(), g->hm.c_str(), o.contact, vol, vol / double(ops.size()), done, ops.size());
    throw std::runtime_error(b);
  }
  res.steps = done;
  res.rwp = best;
  res.crystal = crystal(geometry(best_s), nullptr);
  res.notes.push_back(std::to_string(ops.size()) + " symmetry copies in " + g->hm + ", " + std::to_string(tors.size()) + " torsions, " + std::to_string(done) +
                      " moves (" + std::to_string(restarts) + (restarts == 1 ? " run)" : " runs)"));
  return res;
}

}  // namespace caps
