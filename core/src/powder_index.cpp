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

struct Indexer {
  const IndexOptions& o;
  std::string sys;
  std::vector<double> Q, tol;           // observed Q (1/Å²) and their tolerances
  std::vector<I3> hkl;                  // the indices to try
  std::vector<IndexSolution> found;
  int dims = 0;                         // free parameters
  long tests = 0;

  // the Q range of reflection h over the box
  std::pair<double, double> qrange(const I3& h, const Box& b) const {
    const double h2 = double(h[0]) * h[0], k2 = double(h[1]) * h[1], l2 = double(h[2]) * h[2];
    if (sys == "cubic") {
      const double N = h2 + k2 + l2;
      return {N / (b.hi[0] * b.hi[0]), N / (b.lo[0] * b.lo[0])};
    }
    if (sys == "tetragonal") return {(h2 + k2) / (b.hi[0] * b.hi[0]) + l2 / (b.hi[2] * b.hi[2]), (h2 + k2) / (b.lo[0] * b.lo[0]) + l2 / (b.lo[2] * b.lo[2])};
    if (sys == "hexagonal") {
      const double M = 4.0 / 3.0 * (h2 + double(h[0]) * h[1] + k2);
      return {M / (b.hi[0] * b.hi[0]) + l2 / (b.hi[2] * b.hi[2]), M / (b.lo[0] * b.lo[0]) + l2 / (b.lo[2] * b.lo[2])};
    }
    if (sys == "orthorhombic")
      return {h2 / (b.hi[0] * b.hi[0]) + k2 / (b.hi[1] * b.hi[1]) + l2 / (b.hi[2] * b.hi[2]), h2 / (b.lo[0] * b.lo[0]) + k2 / (b.lo[1] * b.lo[1]) + l2 / (b.lo[2] * b.lo[2])};
    // monoclinic (b unique): A = 1/(a² sin²β), B = 1/b², C = 1/(c² sin²β), D = −2 cos β / (a c sin²β) ≥ 0 for β ≥ 90°
    const double s1 = std::sin(b.lo[3] / kDeg), s2 = std::sin(b.hi[3] / kDeg);
    const double smax2 = std::max(s1, s2) * std::max(s1, s2), smin2 = std::min(s1, s2) * std::min(s1, s2);
    const double c1 = std::fabs(std::cos(b.lo[3] / kDeg)), c2 = std::fabs(std::cos(b.hi[3] / kDeg));
    const double Alo = 1 / (b.hi[0] * b.hi[0] * smax2), Ahi = 1 / (b.lo[0] * b.lo[0] * smin2);
    const double Clo = 1 / (b.hi[2] * b.hi[2] * smax2), Chi = 1 / (b.lo[2] * b.lo[2] * smin2);
    const double Dlo = 2 * std::min(c1, c2) / (b.hi[0] * b.hi[2] * smax2), Dhi = 2 * std::max(c1, c2) / (b.lo[0] * b.lo[2] * smin2);
    const double hl = double(h[0]) * h[2];
    double lo = h2 * Alo + k2 / (b.hi[1] * b.hi[1]) + l2 * Clo, hi = h2 * Ahi + k2 / (b.lo[1] * b.lo[1]) + l2 * Chi;
    if (hl >= 0) lo += hl * Dlo, hi += hl * Dhi;
    else lo += hl * Dhi, hi += hl * Dlo;
    return {lo, hi};
  }
  bool passes(const Box& b) {
    ++tests;
    std::vector<std::pair<double, double>> rr;
    const double qmax = Q.back() + tol.back();
    for (const auto& h : hkl) {
      const auto r = qrange(h, b);
      if (r.first <= qmax) rr.push_back(r);
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
    if (sys == "cubic") return b.lo[0] * b.lo[0] * b.lo[0];
    if (sys == "tetragonal") return b.lo[0] * b.lo[0] * b.lo[2];
    if (sys == "hexagonal") return std::sqrt(3.0) / 2 * b.lo[0] * b.lo[0] * b.lo[2];
    if (sys == "orthorhombic") return b.lo[0] * b.lo[1] * b.lo[2];
    return b.lo[0] * b.lo[1] * b.lo[2] * std::min(std::sin(b.lo[3] / kDeg), std::sin(b.hi[3] / kDeg));
  }
  // least squares on Q with the nearest indices; the cell, M20 and F_N
  void refine_leaf(const Box& b) {
    double p[4];
    for (int k = 0; k < 4; ++k) p[k] = 0.5 * (b.lo[k] + b.hi[k]);
    // the reciprocal metric parameters (linear in Q) for the system
    auto coeffs = [&](const I3& h) -> std::vector<double> {
      const double h2 = double(h[0]) * h[0], k2 = double(h[1]) * h[1], l2 = double(h[2]) * h[2];
      if (sys == "cubic") return {h2 + k2 + l2};
      if (sys == "tetragonal") return {h2 + k2, l2};
      if (sys == "hexagonal") return {4.0 / 3.0 * (h2 + double(h[0]) * h[1] + k2), l2};
      if (sys == "orthorhombic") return {h2, k2, l2};
      return {h2, k2, l2, double(h[0]) * h[2]};
    };
    std::vector<double> x;   // parameters
    const double sb = std::sin(p[3] / kDeg), cb = std::cos(p[3] / kDeg);
    if (sys == "cubic") x = {1 / (p[0] * p[0])};
    else if (sys == "tetragonal" || sys == "hexagonal") x = {1 / (p[0] * p[0]), 1 / (p[2] * p[2])};
    else if (sys == "orthorhombic") x = {1 / (p[0] * p[0]), 1 / (p[1] * p[1]), 1 / (p[2] * p[2])};
    else x = {1 / (p[0] * p[0] * sb * sb), 1 / (p[1] * p[1]), 1 / (p[2] * p[2] * sb * sb), -2 * cb / (p[0] * p[2] * sb * sb)};
    const size_t m = x.size();
    auto qcalc = [&](const I3& h) {
      const auto c = coeffs(h);
      double q = 0;
      for (size_t k = 0; k < m; ++k) q += c[k] * x[k];
      return q;
    };
    std::vector<I3> as(Q.size());
    std::vector<char> ok(Q.size());
    for (int it = 0; it < 4; ++it) {
      for (size_t i = 0; i < Q.size(); ++i) {
        double best = 1e300;
        for (const auto& h : hkl) {
          const double d = std::fabs(qcalc(h) - Q[i]);
          if (d < best) best = d, as[i] = h;
        }
        // the box's centre can sit a few tolerances off: the first passes take every line's nearest reflection
        ok[i] = it < 2 || best <= tol[i];
      }
      std::vector<double> A(m * m, 0.0), rhs(m, 0.0);
      for (size_t i = 0; i < Q.size(); ++i) {
        if (!ok[i]) continue;
        const auto c = coeffs(as[i]);
        for (size_t a = 0; a < m; ++a) {
          rhs[a] += c[a] * Q[i];
          for (size_t bb = 0; bb < m; ++bb) A[a * m + bb] += c[a] * c[bb];
        }
      }
      // solve (Gauss–Jordan)
      bool sing = false;
      for (size_t c = 0; c < m && !sing; ++c) {
        size_t piv = c;
        for (size_t i = c + 1; i < m; ++i) if (std::fabs(A[i * m + c]) > std::fabs(A[piv * m + c])) piv = i;
        if (std::fabs(A[piv * m + c]) < 1e-300) { sing = true; break; }
        for (size_t j = 0; j < m; ++j) std::swap(A[c * m + j], A[piv * m + j]);
        std::swap(rhs[c], rhs[piv]);
        const double d = A[c * m + c];
        for (size_t j = 0; j < m; ++j) A[c * m + j] /= d;
        rhs[c] /= d;
        for (size_t i = 0; i < m; ++i)
          if (i != c) {
            const double f = A[i * m + c];
            for (size_t j = 0; j < m; ++j) A[i * m + j] -= f * A[c * m + j];
            rhs[i] -= f * rhs[c];
          }
      }
      if (sing) return;
      x = rhs;
    }
    for (double v : x) if (!std::isfinite(v)) return;
    IndexSolution s;
    s.system = sys;
    if (sys == "cubic") s.cell = {1 / std::sqrt(x[0]), 1 / std::sqrt(x[0]), 1 / std::sqrt(x[0]), 90, 90, 90};
    else if (sys == "tetragonal") s.cell = {1 / std::sqrt(x[0]), 1 / std::sqrt(x[0]), 1 / std::sqrt(x[1]), 90, 90, 90};
    else if (sys == "hexagonal") s.cell = {1 / std::sqrt(x[0]), 1 / std::sqrt(x[0]), 1 / std::sqrt(x[1]), 90, 90, 120};
    else if (sys == "orthorhombic") s.cell = {1 / std::sqrt(x[0]), 1 / std::sqrt(x[1]), 1 / std::sqrt(x[2]), 90, 90, 90};
    else {
      if (x[0] <= 0 || x[1] <= 0 || x[2] <= 0) return;   // only the h*l term x[3] may be negative
      const double cosb = -x[3] / (2 * std::sqrt(x[0] * x[2]));
      if (std::fabs(cosb) >= 1) return;
      const double beta = std::acos(cosb) * kDeg, sn = std::sin(beta / kDeg);
      s.cell = {1 / (std::sqrt(x[0]) * sn), 1 / std::sqrt(x[1]), 1 / (std::sqrt(x[2]) * sn), 90, beta, 90};
    }
    for (double v : x) if (v <= 0 && sys != "monoclinic") return;
    const auto& c = s.cell;
    const double ca = std::cos(c[3] / kDeg), cbb = std::cos(c[4] / kDeg), cg = std::cos(c[5] / kDeg);
    s.volume = c[0] * c[1] * c[2] * std::sqrt(std::max(0.0, 1 - ca * ca - cbb * cbb - cg * cg + 2 * ca * cbb * cg));
    if (s.volume > o.max_volume || c[0] > o.max_length * 1.05 || c[1] > o.max_length * 1.05 || c[2] > o.max_length * 1.05) return;
    // figures of merit
    double eps = 0, d2t = 0;
    int nind = 0;
    for (size_t i = 0; i < Q.size(); ++i) {
      double best = 1e300;
      I3 bh{0, 0, 0};
      for (const auto& h : hkl) {
        const double d = std::fabs(qcalc(h) - Q[i]);
        if (d < best) best = d, bh = h;
      }
      const bool in = best <= tol[i];
      s.hkl.push_back(in ? bh : I3{0, 0, 0});
      const double sn2 = std::sqrt(std::max(0.0, Q[i])) * o.wavelength / 2;
      s.two_theta_obs.push_back(2 * std::asin(std::min(1.0, sn2)) * kDeg);
      const double qc = qcalc(bh), snc = std::sqrt(std::max(0.0, qc)) * o.wavelength / 2;
      s.two_theta_calc.push_back(2 * std::asin(std::min(1.0, snc)) * kDeg);
      // each line's error floored at a tenth of its tolerance (what a peak search can place): exact data would otherwise
      // give any cell that indexes the lines an unbounded merit
      if (in) ++nind, eps += std::max(best, 0.1 * tol[i]), d2t += std::max(std::fabs(s.two_theta_obs.back() - s.two_theta_calc.back()), 0.1 * o.tolerance);
    }
    if (nind < int(Q.size()) - o.unindexed || nind == 0) return;
    s.indexed = nind;
    std::set<long> distinct;
    for (const auto& h : hkl) {
      const double q = qcalc(h);
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
  double vol_limit = 0;   // once a cell is found: boxes beyond about twice its volume hold only its supercells
  void search(const Box& b, int depth) {
    if (vol_min(b) > o.max_volume) return;
    if (vol_limit > 0 && vol_min(b) > vol_limit) return;
    if (!passes(b)) return;
    double wmax = 0;
    for (int k = 0; k < 3; ++k) wmax = std::max(wmax, b.hi[k] - b.lo[k]);
    const bool beta_wide = sys == "monoclinic" && b.hi[3] - b.lo[3] > 0.02;
    if ((wmax < 0.02 && !beta_wide) || depth > 30) {
      const size_t before = found.size();
      refine_leaf(b);
      if (found.size() > before && found.back().M20 >= 5) {
        const double v = found.back().volume * 2.05;
        vol_limit = vol_limit > 0 ? std::min(vol_limit, v) : v;
      }
      return;
    }
    // halve the widest free dimension(s)
    std::vector<int> dimsplit;
    if (sys == "cubic") dimsplit = {0};
    else if (sys == "tetragonal" || sys == "hexagonal") dimsplit = {0, 2};
    else if (sys == "orthorhombic") dimsplit = {0, 1, 2};
    else dimsplit = {0, 1, 2, 3};
    std::vector<Box> boxes{b};
    for (int d : dimsplit) {
      const double w = b.hi[d] - b.lo[d];
      if ((d < 3 && w < 0.02) || (d == 3 && w < 0.02)) continue;
      std::vector<Box> nb;
      for (const auto& x : boxes) {
        Box l = x, u = x;
        const double mid = 0.5 * (x.lo[d] + x.hi[d]);
        l.hi[d] = mid, u.lo[d] = mid;
        nb.push_back(l), nb.push_back(u);
      }
      boxes.swap(nb);
    }
    for (const auto& x : boxes) search(x, depth + 1);
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
  res.lines = tt;
  for (const auto& sys : o.systems) {
    if (o.progress && !o.progress(sys, 0)) throw std::runtime_error("index stopped");
    Indexer ix{o, sys, {}, {}, {}, {}, 0, 0};
    for (double t : tt) {
      const double th = t / 2 / kDeg, s = std::sin(th);
      ix.Q.push_back(4 * s * s / (o.wavelength * o.wavelength));
      ix.tol.push_back(4 / (o.wavelength * o.wavelength) * s * std::cos(th) * (o.tolerance / kDeg));
    }
    const double qmax = ix.Q.back() + ix.tol.back();
    const int H = int(o.max_length * std::sqrt(qmax)) + 1;
    for (int h = (sys == "monoclinic" ? -H : 0); h <= H; ++h)
      for (int k = 0; k <= H; ++k)
        for (int l = 0; l <= H; ++l) {
          if (!h && !k && !l) continue;
          if (sys == "cubic" && !(h >= k && k >= l)) continue;
          if ((sys == "tetragonal") && h < k) continue;
          if (sys == "hexagonal" && h < k) continue;
          ix.hkl.push_back({h, k, l});
        }
    const double lmin = 2.0, step = 0.5;
    std::vector<Box> start;
    for (double a = lmin; a < o.max_length; a += step) {
      if (sys == "cubic") { start.push_back({{a, a, a, 90}, {a + step, a + step, a + step, 90}}); continue; }
      for (double c = lmin; c < o.max_length; c += step) {
        if (sys == "tetragonal" || sys == "hexagonal") { start.push_back({{a, a, c, 90}, {a + step, a + step, c + step, 90}}); continue; }
        for (double b = (sys == "orthorhombic" ? a : lmin); b < (sys == "orthorhombic" ? c + step : o.max_length); b += step) {
          if (sys == "orthorhombic") {
            if (b > c) continue;
            start.push_back({{a, b, c, 90}, {a + step, b + step, c + step, 90}});
          } else {
            if (a > c) continue;   // a ≤ c for the monoclinic setting (the other order is the same lattice)
            for (double be = o.min_beta; be < o.max_beta; be += 5) start.push_back({{a, b, c, be}, {a + step, b + step, c + step, std::min(be + 5, o.max_beta)}});
          }
        }
      }
    }
    // the smallest cells first (the indexing sought is the smallest cell that indexes the lines)
    std::sort(start.begin(), start.end(), [&](const Box& x, const Box& y) { return ix.vol_min(x) < ix.vol_min(y); });
    for (size_t k = 0; k < start.size(); ++k) {
      if (ix.vol_limit > 0 && ix.vol_min(start[k]) > ix.vol_limit) break;
      if (o.progress && k % 200 == 0 && !o.progress(sys, double(k) / double(start.size()))) throw std::runtime_error("index stopped");
      ix.search(start[k], 0);
    }
    for (auto& f : ix.found) res.solutions.push_back(f);
    char b[160];
    std::snprintf(b, sizeof b, "%s: %zu cells (%ld boxes tested)", sys.c_str(), ix.found.size(), ix.tests);
    res.notes.push_back(b);
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
    std::vector<std::vector<Vec3>> copies;
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
      for (size_t i = 0; i < n; ++i) cart[i] = cell.to_cartesian(f[i] - sh);
      copies.push_back(cart);
    }
    if (clash) {
      *clash = false;
      const double c2 = o.contact * o.contact;
      for (size_t a = 0; a < copies.size() && !*clash; ++a)
        for (size_t b = a; b < copies.size() && !*clash; ++b)
          for (size_t i = 0; i < n && !*clash; ++i) {
            if (mol.atoms[i].element == 1) continue;
            for (size_t j = 0; j < n; ++j) {
              if (mol.atoms[j].element == 1 || (a == b && i == j)) continue;
              if (a == b) continue;   // the molecule's own atoms (images of itself are checked through other copies' cells)
              const Vec3 d = cell.minimum_image(copies[b][j] - copies[a][i]);
              if (dot(d, d) < c2) { *clash = true; break; }
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
  RefineOptions ro = o.refine;
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
  double best = 1e300;
  St best_s;
  const int total = std::max(1, o.restarts) * o.steps;
  int done = 0;
  for (int rs = 0; rs < std::max(1, o.restarts); ++rs) {
    St s{{U(rng), U(rng), U(rng)}, rand_axis(), 2 * kPi * U(rng), std::vector<double>(tors.size(), 0.0)};
    for (auto& t : s.tor) t = 2 * kPi * U(rng);
    bool clash = false;
    double cur = rwp_of(crystal(geometry(s), &clash));
    if (clash) cur = 1e3;
    // the starting state is a candidate too, so best_s is always a full state even when no move is accepted
    // (rs == 0 also covers a NaN Rwp, e.g. an all-zero pattern, which no comparison would accept)
    if (rs == 0 || cur < best) best = cur, best_s = s;
    for (int it = 0; it < o.steps; ++it, ++done) {
      if (o.progress && it % 50 == 0 && !o.progress("annealing", double(done) / total)) throw std::runtime_error("solve stopped");
      // the last fifth: downhill only, with small moves (a polish of the best region)
      const bool polish = it >= o.steps * 4 / 5;
      const double T = polish ? 0 : o.t_start * std::pow(o.t_end / o.t_start, double(it) / std::max(1, o.steps * 4 / 5 - 1));
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
        if (v < best) best = v, best_s = t;
      }
      if (it % 100 == 0) res.rwp_history.push_back(cur);
    }
  }
  res.steps = done;
  res.rwp = best;
  res.crystal = crystal(geometry(best_s), nullptr);
  res.notes.push_back(std::to_string(ops.size()) + " symmetry copies in " + g->hm + ", " + std::to_string(tors.size()) + " torsions, " + std::to_string(done) + " moves");
  return res;
}

}  // namespace caps
