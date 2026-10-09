// CAPS bonded coarse-grained potentials by tabulated Boltzmann inversion (see caps/cg_bonded.hpp).
#include "caps/cg_bonded.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <limits>
#include <numeric>
#include <set>
#include <stdexcept>

namespace caps {

namespace {
constexpr double kB = 0.0019872067;   // kcal/(mol K)
constexpr double kPi = 3.14159265358979323846;
constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

std::vector<double> gauss_smooth(const std::vector<double>& y, double s, bool periodic) {
  if (s <= 0 || y.empty()) return y;
  const int w = int(std::ceil(3 * s));
  std::vector<double> k(size_t(2 * w + 1));
  for (int i = -w; i <= w; ++i) k[size_t(i + w)] = std::exp(-0.5 * (i / s) * (i / s));
  const int n = int(y.size());
  std::vector<double> out(y.size(), 0.0);
  for (int i = 0; i < n; ++i) {
    double a = 0, wsum = 0;
    for (int j = -w; j <= w; ++j) {
      int q = i + j;
      if (periodic) q = ((q % n) + n) % n;
      else if (q < 0 || q >= n) continue;   // renormalised at the edges (no padding bias)
      a += k[size_t(j + w)] * y[size_t(q)];
      wsum += k[size_t(j + w)];
    }
    out[size_t(i)] = a / wsum;
  }
  return out;
}

// natural cubic spline through (x, y) — x increasing — evaluated with its derivative
struct Spline {
  std::vector<double> x, y, m;   // m: second derivatives
  bool periodic = false;
  double period = 0;
  void build(std::vector<double> xs, std::vector<double> ys, bool per = false, double p = 0) {
    x = std::move(xs), y = std::move(ys), periodic = per, period = p;
    const size_t n = x.size();
    m.assign(n, 0.0);
    if (n < 3) return;
    if (!periodic) {   // tridiagonal, natural ends
      std::vector<double> a(n, 0), b(n, 1), c(n, 0), d(n, 0);
      for (size_t i = 1; i + 1 < n; ++i) {
        const double h0 = x[i] - x[i - 1], h1 = x[i + 1] - x[i];
        a[i] = h0 / 6, b[i] = (h0 + h1) / 3, c[i] = h1 / 6;
        d[i] = (y[i + 1] - y[i]) / h1 - (y[i] - y[i - 1]) / h0;
      }
      for (size_t i = 1; i < n; ++i) { const double w = a[i] / b[i - 1]; b[i] -= w * c[i - 1]; d[i] -= w * d[i - 1]; }
      m[n - 1] = d[n - 1] / b[n - 1];
      for (size_t i = n - 1; i-- > 0;) m[i] = (d[i] - c[i] * m[i + 1]) / b[i];
      m[0] = m[n - 1] = 0;
    } else {   // even spacing on a circle: solve the cyclic system by iteration (diagonally dominant)
      const double h = x[1] - x[0];
      std::vector<double> d(n);
      for (size_t i = 0; i < n; ++i) d[i] = 6 * (y[(i + 1) % n] - 2 * y[i] + y[(i + n - 1) % n]) / (h * h);
      for (int it = 0; it < 500; ++it) {
        double change = 0;
        for (size_t i = 0; i < n; ++i) {
          const double v = (d[i] - m[(i + 1) % n] - m[(i + n - 1) % n]) / 4;
          change = std::max(change, std::fabs(v - m[i]));
          m[i] = v;
        }
        if (change < 1e-12) break;
      }
    }
  }
  // value and derivative at t
  std::pair<double, double> at(double t) const {
    const size_t n = x.size();
    if (n == 1) return {y[0], 0};
    size_t i;
    double x0, x1, y0, y1, m0, m1;
    if (periodic) {
      const double h = x[1] - x[0];
      double u = std::fmod(t - x[0], period);
      if (u < 0) u += period;
      i = std::min(n - 1, size_t(u / h));
      const size_t j = (i + 1) % n;
      x0 = x[0] + double(i) * h, x1 = x0 + h, y0 = y[i], y1 = y[j], m0 = m[i], m1 = m[j];
      t = x[0] + u;
    } else {
      i = size_t(std::upper_bound(x.begin(), x.end(), t) - x.begin());
      i = std::clamp<size_t>(i, 1, n - 1) - 1;
      x0 = x[i], x1 = x[i + 1], y0 = y[i], y1 = y[i + 1], m0 = m[i], m1 = m[i + 1];
    }
    const double h = x1 - x0, A = (x1 - t) / h, B = (t - x0) / h;
    const double v = A * y0 + B * y1 + ((A * A * A - A) * m0 + (B * B * B - B) * m1) * h * h / 6;
    const double dv = (y1 - y0) / h + (-(3 * A * A - 1) * m0 + (3 * B * B - 1) * m1) * h / 6;
    return {v, dv};
  }
};

// least-squares parabola through (x, y): returns the curvature c2 of y = c0 + c1 x + c2 x²
double parabola_c2(const std::vector<double>& x, const std::vector<double>& y) {
  const size_t n = x.size();
  if (n < 3) return 0;
  double mx = 0;
  for (double v : x) mx += v;
  mx /= double(n);
  double S[5] = {0, 0, 0, 0, 0}, T[3] = {0, 0, 0};
  for (size_t i = 0; i < n; ++i) {
    const double u = x[i] - mx;
    double p = 1;
    for (int k = 0; k < 5; ++k) { S[k] += p; if (k < 3) T[k] += p * y[i]; p *= u; }
  }
  // normal equations [S0 S1 S2; S1 S2 S3; S2 S3 S4] c = T
  const double M[3][3] = {{S[0], S[1], S[2]}, {S[1], S[2], S[3]}, {S[2], S[3], S[4]}};
  auto det3 = [](const double A[3][3]) {
    return A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) - A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0]) + A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);
  };
  const double D = det3(M);
  if (std::fabs(D) < 1e-300) return 0;
  double M2[3][3];
  for (int r = 0; r < 3; ++r) for (int c = 0; c < 3; ++c) M2[r][c] = c == 2 ? T[r] : M[r][c];
  return det3(M2) / D;
}

const char* kind_name(int k) { return k == 0 ? "bond" : k == 1 ? "angle" : "dihedral"; }

}  // namespace

// ------------------------------------------------------------------------------------------------ geometry
double cg_angle_deg(const Vec3& a, const Vec3& b, const Vec3& c) {
  const Vec3 u = a - b, v = c - b;
  const double cs = dot(u, v) / (norm(u) * norm(v));
  return std::acos(std::clamp(cs, -1.0, 1.0)) * 180.0 / kPi;
}
double cg_dihedral_deg(const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& d) {
  const Vec3 b1 = b - a, b2 = c - b, b3 = d - c;
  const Vec3 n1 = cross(b1, b2), n2 = cross(b2, b3);
  // IUPAC sign (Blondel & Karplus 1996): positive when, seen along b → c, a turns clockwise onto d; trans = ±180
  return std::atan2(norm(b2) * dot(b1, n2), dot(n1, n2)) * 180.0 / kPi;
}

// ------------------------------------------------------------------------------------------------ topology
CgTopology cg_topology_from_map(const Json& j) {
  if (j.text("format") != "caps-cg-map") throw std::invalid_argument("not a CAPS mapping (format caps-cg-map)");
  CgTopology t;
  for (const auto& b : j["beads"].items()) t.kind.push_back(b["kind"].str()), t.mol.push_back(int(b.num("mol", 0)));
  for (const auto& b : j["bonds"].items()) t.bonds.push_back({int(b[0].number()), int(b[1].number())});
  if (j.has("chains"))
    for (const auto& c : j["chains"].items()) {
      std::vector<int> v;
      for (const auto& x : c.items()) v.push_back(int(x.number()));
      t.chains.push_back(v);
    }
  std::vector<std::vector<int>> nbr(t.kind.size());
  for (const auto& [x, y] : t.bonds) nbr[size_t(x)].push_back(y), nbr[size_t(y)].push_back(x);
  for (auto& v : nbr) std::sort(v.begin(), v.end());
  for (size_t j2 = 0; j2 < nbr.size(); ++j2)
    for (size_t a = 0; a < nbr[j2].size(); ++a)
      for (size_t c = a + 1; c < nbr[j2].size(); ++c) t.angles.push_back({nbr[j2][a], int(j2), nbr[j2][c]});
  for (const auto& [b, c] : t.bonds)
    for (int a : nbr[size_t(b)])
      if (a != c)
        for (int d : nbr[size_t(c)])
          if (d != b && d != a) t.dihedrals.push_back({a, b, c, d});
  return t;
}

CgTopology cg_topology_of(const CgMapping& m) {
  CgTopology t;
  t.kind = m.bead_kind;
  t.mol = m.bead_mol;
  t.bonds = m.bonds;
  t.angles = m.angles;
  t.dihedrals = m.dihedrals;
  t.chains = m.chains;
  return t;
}

// ------------------------------------------------------------------------------------------------ accumulation
CgBondedAccumulator::CgBondedAccumulator(const CgTypes& types, const CgBondedOptions& o) : types_(types), o_(o) {
  if (o.temperature <= 0) throw std::invalid_argument("give the temperature of the inversion");
  if (o.bond_bin <= 0 || o.angle_bin <= 0 || o.dihedral_bin <= 0) throw std::invalid_argument("histogram bins must be positive");
}

size_t CgBondedAccumulator::nbins(int kind) const {
  if (kind == 0) return size_t(std::ceil(o_.bond_max / o_.bond_bin));
  if (kind == 1) return size_t(std::lround(180.0 / o_.angle_bin));
  return size_t(std::lround(360.0 / o_.dihedral_bin));
}

int CgBondedAccumulator::add_system(const CgTopology& t, double weight, const std::string& name) {
  if (weight <= 0) throw std::invalid_argument("a system's weight must be positive");
  Sys s;
  s.t = t;
  s.weight = weight;
  s.name = name;
  auto idx = [&](const std::vector<std::string>& list, const std::string& key, const char* what) {
    const auto it = std::find(list.begin(), list.end(), key);
    if (it == list.end()) throw std::invalid_argument(std::string("the shared type list has no ") + what + " type " + key + (name.empty() ? "" : " (" + name + ")"));
    return int(it - list.begin());
  };
  for (const auto& [x, y] : t.bonds) s.bond_t.push_back(idx(types_.bonds, cg_key({t.kind[size_t(x)], t.kind[size_t(y)]}), "bond"));
  for (const auto& q : t.angles) s.angle_t.push_back(idx(types_.angles, cg_key({t.kind[size_t(q[0])], t.kind[size_t(q[1])], t.kind[size_t(q[2])]}), "angle"));
  for (const auto& q : t.dihedrals)
    s.dih_t.push_back(idx(types_.dihedrals, cg_key({t.kind[size_t(q[0])], t.kind[size_t(q[1])], t.kind[size_t(q[2])], t.kind[size_t(q[3])]}), "dihedral"));
  const size_t nt[3] = {types_.bonds.size(), types_.angles.size(), types_.dihedrals.size()};
  for (int k = 0; k < 3; ++k)
    for (int h = 0; h < 2; ++h) s.h[k][h].assign(nt[k], std::vector<double>(nbins(k), 0.0));
  sys_.push_back(std::move(s));
  return int(sys_.size()) - 1;
}

void CgBondedAccumulator::set_expected_frames(int s, size_t frames) { sys_.at(size_t(s)).expected = frames; }

void CgBondedAccumulator::add_frame(int si, const std::vector<Vec3>& p, const Cell& cell) {
  Sys& s = sys_.at(size_t(si));
  if (p.size() != s.t.beads()) throw std::invalid_argument("a frame of " + s.name + " has " + std::to_string(p.size()) + " beads, its mapping " + std::to_string(s.t.beads()));
  // contiguous halves when the frame count is known, else alternate frames
  const int half = s.expected > 1 ? (s.frames < s.expected / 2 ? 0 : 1) : int(s.frames % 2);
  const bool pbc = cell.valid();
  auto vec = [&](int a, int b) { const Vec3 d = p[size_t(b)] - p[size_t(a)]; return pbc ? cell.minimum_image(d) : d; };
  const double nb0 = double(nbins(0)), nb1 = double(nbins(1)), nb2 = double(nbins(2));
  for (size_t k = 0; k < s.t.bonds.size(); ++k) {
    const double r = norm(vec(s.t.bonds[k].first, s.t.bonds[k].second));
    const double b = std::floor(r / o_.bond_bin);
    if (b >= 0 && b < nb0) s.h[0][half][size_t(s.bond_t[k])][size_t(b)] += 1;
  }
  for (size_t k = 0; k < s.t.angles.size(); ++k) {
    const auto& q = s.t.angles[k];
    const Vec3 u = vec(q[1], q[0]), v = vec(q[1], q[2]);
    const double th = std::acos(std::clamp(dot(u, v) / (norm(u) * norm(v)), -1.0, 1.0)) * 180.0 / kPi;
    const double b = std::min(nb1 - 1, std::floor(th / o_.angle_bin));
    s.h[1][half][size_t(s.angle_t[k])][size_t(b)] += 1;
  }
  for (size_t k = 0; k < s.t.dihedrals.size(); ++k) {
    const auto& q = s.t.dihedrals[k];
    const Vec3 b1 = vec(q[0], q[1]), b2 = vec(q[1], q[2]), b3 = vec(q[2], q[3]);
    const Vec3 n1 = cross(b1, b2), n2 = cross(b2, b3);
    const double ph = std::atan2(norm(b2) * dot(b1, n2), dot(n1, n2)) * 180.0 / kPi;   // IUPAC, as cg_dihedral_deg
    const double b = std::min(nb2 - 1, std::floor((ph + 180.0) / o_.dihedral_bin));
    s.h[2][half][size_t(s.dih_t[k])][size_t(std::max(0.0, b))] += 1;
  }
  ++s.frames;
}

size_t CgBondedAccumulator::frames() const {
  size_t n = 0;
  for (const auto& s : sys_) n += s.frames;
  return n;
}

std::vector<std::string> CgBondedAccumulator::notes() const { return notes_; }

std::vector<CgBondedHistogram> CgBondedAccumulator::pooled() const {
  std::vector<CgBondedHistogram> out;
  const std::vector<std::string>* lists[3] = {&types_.bonds, &types_.angles, &types_.dihedrals};
  const double lo[3] = {0, 0, -180}, bin[3] = {o_.bond_bin, o_.angle_bin, o_.dihedral_bin};
  for (int k = 0; k < 3; ++k)
    for (size_t t = 0; t < lists[k]->size(); ++t) {
      CgBondedHistogram H;
      H.kind = k;
      H.key = (*lists[k])[t];
      H.lo = lo[k];
      H.bin = bin[k];
      const size_t nb = nbins(k);
      H.h.assign(nb, 0.0);
      H.half[0].assign(nb, 0.0);
      H.half[1].assign(nb, 0.0);
      double wsum = 0, w0 = 0, w1 = 0;
      for (const auto& s : sys_) {
        double n0 = 0, n1 = 0;
        for (size_t b = 0; b < nb; ++b) n0 += s.h[k][0][t][b], n1 += s.h[k][1][t][b];
        if (n0 + n1 <= 0) continue;
        H.count += long(n0 + n1);
        wsum += s.weight;
        for (size_t b = 0; b < nb; ++b) H.h[b] += s.weight * (s.h[k][0][t][b] + s.h[k][1][t][b]) / (n0 + n1);
        if (n0 > 0) { w0 += s.weight; for (size_t b = 0; b < nb; ++b) H.half[0][b] += s.weight * s.h[k][0][t][b] / n0; }
        if (n1 > 0) { w1 += s.weight; for (size_t b = 0; b < nb; ++b) H.half[1][b] += s.weight * s.h[k][1][t][b] / n1; }
      }
      if (wsum > 0) for (auto& v : H.h) v /= wsum;
      if (w0 > 0) for (auto& v : H.half[0]) v /= w0;
      if (w1 > 0) for (auto& v : H.half[1]) v /= w1;
      out.push_back(std::move(H));
    }
  return out;
}

// ------------------------------------------------------------------------------------------------ inversion
namespace {

struct Inverted {
  std::vector<double> xc, P, U;   // bin centres, P (Jacobian removed, max 1), U (NaN outside the inverted range)
  int ilo = -1, ihi = -1;          // inverted bins
};

Inverted invert_hist(const std::vector<double>& h, int kind, double lo, double bin, const CgBondedOptions& o, double kT) {
  Inverted r;
  const size_t n = h.size();
  r.xc.resize(n);
  for (size_t i = 0; i < n; ++i) r.xc[i] = lo + (double(i) + 0.5) * bin;
  const std::vector<double> s = gauss_smooth(h, o.smooth, kind == 2);
  r.P.assign(n, 0.0);
  for (size_t i = 0; i < n; ++i) {
    const double jac = kind == 0 ? r.xc[i] * r.xc[i] : kind == 1 ? std::sin(r.xc[i] * kPi / 180.0) : 1.0;
    r.P[i] = jac > 1e-12 ? s[i] / jac : 0.0;
  }
  const double pmax = *std::max_element(r.P.begin(), r.P.end());
  r.U.assign(n, kNaN);
  if (pmax <= 0) return r;
  for (auto& v : r.P) v /= pmax;
  for (size_t i = 0; i < n; ++i)
    if (r.P[i] > o.threshold) { if (r.ilo < 0) r.ilo = int(i); r.ihi = int(i); }
  for (int i = r.ilo; i <= r.ihi; ++i)
    if (r.P[size_t(i)] > 0) r.U[size_t(i)] = -kT * std::log(r.P[size_t(i)]);
  // empty bins inside the range: linear in U between their sampled neighbours
  for (int i = r.ilo; i <= r.ihi; ++i)
    if (std::isnan(r.U[size_t(i)])) {
      int a = i - 1, b = i + 1;
      while (b <= r.ihi && std::isnan(r.U[size_t(b)])) ++b;
      const double t = double(i - a) / double(b - a);
      r.U[size_t(i)] = r.U[size_t(a)] + t * (r.U[size_t(b)] - r.U[size_t(a)]);
    }
  return r;
}

// the table of a non-periodic coordinate from its inverted bins: a spline inside, walls outside
void table_open(CgBondedTable& T, const Inverted& v, int kind, const CgBondedOptions& o, double kT, double tlo, double thi, double dx) {
  std::vector<double> xs, us;
  for (int i = v.ilo; i <= v.ihi; ++i) xs.push_back(v.xc[size_t(i)]), us.push_back(v.U[size_t(i)]);
  T.lo = xs.front(), T.hi = xs.back();
  Spline sp;
  sp.build(xs, us);
  // the well: its minimum and a parabola through U < 2 kT
  const size_t imin = size_t(std::min_element(us.begin(), us.end()) - us.begin());
  T.x0 = xs[imin];
  std::vector<double> wx, wy;
  for (size_t i = 0; i < xs.size(); ++i) if (us[i] - us[imin] < 2 * kT) wx.push_back(xs[i]), wy.push_back(us[i]);
  T.k_harmonic = std::max(0.0, parabola_c2(wx, wy));
  const double wall = std::max(T.k_harmonic, kind == 0 ? o.bond_wall : o.angle_wall);
  auto [ulo, slo] = sp.at(T.lo);
  auto [uhi, shi] = sp.at(T.hi);
  slo = std::min(slo, 0.0);   // the walls only rise outward
  shi = std::max(shi, 0.0);
  const int n = int(std::floor((thi - tlo) / dx + 1e-9)) + 1;
  T.x.resize(size_t(n)), T.U.resize(size_t(n)), T.F.resize(size_t(n));
  for (int i = 0; i < n; ++i) {
    const double x = tlo + i * dx;
    double u, du;
    if (x < T.lo) u = ulo + slo * (x - T.lo) + wall * (x - T.lo) * (x - T.lo), du = slo + 2 * wall * (x - T.lo);
    else if (x > T.hi) u = uhi + shi * (x - T.hi) + wall * (x - T.hi) * (x - T.hi), du = shi + 2 * wall * (x - T.hi);
    else std::tie(u, du) = sp.at(x);
    T.x[size_t(i)] = x, T.U[size_t(i)] = u, T.F[size_t(i)] = -du;
  }
  const double umin = *std::min_element(T.U.begin(), T.U.end());
  for (auto& u : T.U) u -= umin;
}

// a periodic coordinate (dihedral): sampled bins kept, gaps bridged by a smooth barrier, a periodic spline through all
void table_periodic(CgBondedTable& T, const Inverted& v, const CgBondedOptions& o, double kT) {
  const size_t n = v.xc.size();
  std::vector<double> u(n, kNaN);
  for (size_t i = 0; i < n; ++i) if (v.P[i] > o.threshold) u[i] = -kT * std::log(v.P[i]);
  double umax = 0;
  bool any = false;
  for (double x : u) if (!std::isnan(x)) umax = std::max(umax, x), any = true;
  if (!any) return;
  const double cap = o.dihedral_cap > 0 ? o.dihedral_cap : 2 * kT;
  // each gap (circularly): between its sampled edges a + (b − a) t + cap sin(π t)
  for (size_t i = 0; i < n; ++i) {
    if (!std::isnan(u[i]) || std::isnan(u[(i + n - 1) % n])) continue;
    size_t j = i;
    int len = 0;
    while (std::isnan(u[j % n])) ++j, ++len;
    const double a = u[(i + n - 1) % n], b = u[j % n];
    for (int k = 0; k < len; ++k) {
      const double t = double(k + 1) / double(len + 1);
      u[(i + size_t(k)) % n] = a + (b - a) * t + (std::max(a, b) - std::min(a, b) + cap) * std::sin(kPi * t);
    }
  }
  T.lo = -180, T.hi = 180;
  Spline sp;
  sp.build(v.xc, u, true, 360.0);
  const size_t imin = size_t(std::min_element(u.begin(), u.end()) - u.begin());
  T.x0 = v.xc[imin];
  const int m = int(std::lround(360.0 / o.dihedral_table_dp));
  T.x.resize(size_t(m)), T.U.resize(size_t(m)), T.F.resize(size_t(m));
  for (int i = 0; i < m; ++i) {
    const double x = -180.0 + i * o.dihedral_table_dp;   // −180 … < 180 (LAMMPS: less than a full turn)
    const auto [val, d] = sp.at(x);
    T.x[size_t(i)] = x, T.U[size_t(i)] = val, T.F[size_t(i)] = -d;
  }
  const double umin = *std::min_element(T.U.begin(), T.U.end());
  for (auto& x : T.U) x -= umin;
}

double half_difference(const CgBondedHistogram& H, const CgBondedOptions& o, double kT, int ilo, int ihi) {
  const Inverted a = invert_hist(H.half[0], H.kind, H.lo, H.bin, o, kT), b = invert_hist(H.half[1], H.kind, H.lo, H.bin, o, kT);
  if (a.ilo < 0 || b.ilo < 0) return kNaN;
  // compare where both halves are inverted, each shifted to its mean over that common range
  std::vector<std::pair<double, double>> pairs;
  for (int i = std::max({ilo, a.ilo, b.ilo}); i <= std::min({ihi, a.ihi, b.ihi}); ++i) pairs.push_back({a.U[size_t(i)], b.U[size_t(i)]});
  if (pairs.empty()) return kNaN;
  double ma = 0, mb = 0;
  for (const auto& [x, y] : pairs) ma += x, mb += y;
  ma /= double(pairs.size()), mb /= double(pairs.size());
  double worst = 0;
  for (const auto& [x, y] : pairs) worst = std::max(worst, std::fabs((x - ma) - (y - mb)));
  return worst / kT;
}

CgBondedTable make_table(const CgBondedHistogram& H, const CgBondedOptions& o, double kT, std::vector<std::string>& notes) {
  CgBondedTable T;
  T.kind = H.kind;
  T.key = H.key;
  T.count = H.count;
  const Inverted v = invert_hist(H.h, H.kind, H.lo, H.bin, o, kT);
  T.xh = v.xc;
  T.P = v.P;
  T.Uinv = v.U;
  if (v.ilo < 0 || H.count == 0) {
    notes.push_back(std::string("no samples for ") + kind_name(H.kind) + " " + H.key + ": no table");
    return T;
  }
  T.sampled = true;
  if (H.kind == 0) {
    const double tlo = o.bond_table_lo > 0 ? o.bond_table_lo : std::max(0.5, v.xc[size_t(v.ilo)] - 2.0);
    const double thi = o.bond_table_hi > 0 ? o.bond_table_hi : v.xc[size_t(v.ihi)] + 3.0;
    table_open(T, v, 0, o, kT, tlo, thi, o.bond_table_dr);
  } else if (H.kind == 1) {
    table_open(T, v, 1, o, kT, 0.0, 180.0, o.angle_table_dt);
  } else {
    table_periodic(T, v, o, kT);
  }
  T.half_diff = half_difference(H, o, kT, v.ilo, v.ihi);
  if (H.count < 1000) notes.push_back(std::string(kind_name(H.kind)) + " " + H.key + ": only " + std::to_string(H.count) + " samples");
  return T;
}

}  // namespace

CgBondedResult invert_bonded(const CgBondedAccumulator& acc) {
  CgBondedResult r;
  const auto& o = acc.options();
  r.temperature = o.temperature;
  r.frames = acc.frames();
  const double kT = kB * o.temperature;
  for (const auto& H : acc.pooled()) {
    CgBondedTable T = make_table(H, o, kT, r.notes);
    (H.kind == 0 ? r.bonds : H.kind == 1 ? r.angles : r.dihedrals).push_back(std::move(T));
  }
  for (const auto& n : acc.notes()) r.notes.push_back(n);
  return r;
}

CgBondedResult refine_bonded(const CgBondedResult& cur, const CgBondedAccumulator& cg, double alpha, std::map<std::string, double>* residual) {
  CgBondedResult r = cur;
  const auto& o = cg.options();
  const double kT = kB * o.temperature;
  const auto pooled = cg.pooled();
  auto find = [&](int kind, const std::string& key) -> const CgBondedHistogram* {
    for (const auto& h : pooled) if (h.kind == kind && h.key == key) return &h;
    return nullptr;
  };
  for (auto* list : {&r.bonds, &r.angles, &r.dihedrals})
    for (auto& T : *list) {
      if (!T.sampled) continue;
      const CgBondedHistogram* H = find(T.kind, T.key);
      if (!H || H->count == 0) { r.notes.push_back(std::string("the CG run has no ") + kind_name(T.kind) + " " + T.key + ": kept"); continue; }
      const Inverted c = invert_hist(H->h, T.kind, H->lo, H->bin, o, kT);
      if (c.P.size() != T.P.size()) throw std::invalid_argument("the CG run's histograms use other bins than the tables");
      // residual and the correction on the bins, then the table rebuilt from the corrected inverted potential
      double res = 0, norm_t = 0;
      Inverted v;
      v.xc = T.xh, v.P = T.P, v.U = T.Uinv;
      for (size_t i = 0; i < T.P.size(); ++i) {
        res += std::fabs(c.P[i] - T.P[i]);
        norm_t += T.P[i];
        if (!std::isnan(v.U[i]) && c.P[i] > o.threshold && T.P[i] > o.threshold) v.U[i] += alpha * kT * std::log(c.P[i] / T.P[i]);
      }
      if (residual) (*residual)[std::string(kind_name(T.kind)) + " " + T.key] = norm_t > 0 ? res / norm_t : 0;
      v.ilo = -1, v.ihi = -1;
      for (size_t i = 0; i < v.U.size(); ++i) if (!std::isnan(v.U[i])) { if (v.ilo < 0) v.ilo = int(i); v.ihi = int(i); }
      if (v.ilo < 0) continue;
      CgBondedTable N = T;
      N.Uinv = v.U;
      if (T.kind == 0) table_open(N, v, 0, o, kT, T.x.front(), T.x.back(), T.x.size() > 1 ? T.x[1] - T.x[0] : o.bond_table_dr);
      else if (T.kind == 1) table_open(N, v, 1, o, kT, 0.0, 180.0, o.angle_table_dt);
      else {
        // periodic: the corrected U back to a pseudo-distribution so the gap logic applies again
        Inverted q = v;
        for (size_t i = 0; i < q.U.size(); ++i) q.P[i] = std::isnan(q.U[i]) ? 0.0 : std::exp(-q.U[i] / kT);
        const double pm = *std::max_element(q.P.begin(), q.P.end());
        for (auto& x : q.P) x /= pm;
        table_periodic(N, q, o, kT);
      }
      T = std::move(N);
    }
  return r;
}

// ------------------------------------------------------------------------------------------------ files
Json bonded_json(const CgBondedResult& r) {
  Json j = Json::object();
  j["format"] = "caps-cg-bonded";
  j["version"] = 1;
  j["temperature"] = r.temperature;
  j["frames"] = double(r.frames);
  auto arr = [](const std::vector<double>& v) {
    Json a = Json::array();
    for (double x : v) a.push_back(std::isnan(x) ? Json() : Json(x));
    return a;
  };
  auto list = [&](const std::vector<CgBondedTable>& ts) {
    Json a = Json::array();
    for (const auto& T : ts) {
      Json e = Json::object();
      e["key"] = T.key;
      e["kind"] = kind_name(T.kind);
      e["sampled"] = T.sampled;
      e["count"] = double(T.count);
      e["x0"] = T.x0;
      e["k_harmonic"] = T.k_harmonic;
      e["lo"] = T.lo;
      e["hi"] = T.hi;
      e["half_diff_kT"] = std::isnan(T.half_diff) ? Json() : Json(T.half_diff);
      e["x"] = arr(T.x);
      e["U"] = arr(T.U);
      e["F"] = arr(T.F);
      e["xh"] = arr(T.xh);
      e["P"] = arr(T.P);
      e["Uinv"] = arr(T.Uinv);
      a.push_back(e);
    }
    return a;
  };
  j["bonds"] = list(r.bonds);
  j["angles"] = list(r.angles);
  j["dihedrals"] = list(r.dihedrals);
  Json notes = Json::array();
  for (const auto& n : r.notes) notes.push_back(n);
  j["notes"] = notes;
  return j;
}

CgBondedResult bonded_from_json(const Json& j) {
  if (j.text("format") != "caps-cg-bonded") throw std::invalid_argument("not a CAPS bonded set (format caps-cg-bonded)");
  CgBondedResult r;
  r.temperature = j.num("temperature", 300);
  r.frames = size_t(j.num("frames", 0));
  auto vec = [](const Json& a) {
    std::vector<double> v;
    for (const auto& x : a.items()) v.push_back(x.is_null() ? kNaN : x.number());
    return v;
  };
  auto list = [&](const char* k, int kind, std::vector<CgBondedTable>& out) {
    for (const auto& e : j[k].items()) {
      CgBondedTable T;
      T.kind = kind;
      T.key = e["key"].str();
      T.sampled = e["sampled"].boolean();
      T.count = long(e.num("count", 0));
      T.x0 = e.num("x0", 0);
      T.k_harmonic = e.num("k_harmonic", 0);
      T.lo = e.num("lo", 0);
      T.hi = e.num("hi", 0);
      T.half_diff = e["half_diff_kT"].is_null() ? kNaN : e["half_diff_kT"].number();
      T.x = vec(e["x"]), T.U = vec(e["U"]), T.F = vec(e["F"]), T.xh = vec(e["xh"]), T.P = vec(e["P"]), T.Uinv = vec(e["Uinv"]);
      out.push_back(std::move(T));
    }
  };
  list("bonds", 0, r.bonds);
  list("angles", 1, r.angles);
  list("dihedrals", 2, r.dihedrals);
  return r;
}

std::vector<std::string> write_bonded(const CgBondedResult& r, const CgTypes& types, const CgBondedOptions& o, const std::string& dir) {
  namespace fs = std::filesystem;
  fs::create_directories(dir);
  fs::create_directories(fs::path(dir) / "gromacs");
  std::vector<std::string> files;
  char line[256];
  auto open = [&](const std::string& name) {
    const std::string p = (fs::path(dir) / name).string();
    files.push_back(p);
    std::ofstream f(p);
    if (!f) throw std::runtime_error("cannot write " + p);
    return f;
  };
  const std::string head = "# CAPS coarse-grained bonded potentials, Boltzmann-inverted at " + std::to_string(int(std::lround(r.temperature))) + " K from " +
                           std::to_string(r.frames) + " frame(s); kcal/mol\n";
  {
    auto f = open("bonds.table");
    f << head << "# r (Å), E, −dE/dr (kcal/mol/Å)\n";
    for (const auto& T : r.bonds) {
      if (!T.sampled) continue;
      f << "\n" << T.key << "\nN " << T.x.size() << "\n\n";
      for (size_t i = 0; i < T.x.size(); ++i) { std::snprintf(line, sizeof line, "%zu %.5f %.6f %.6f\n", i + 1, T.x[i], T.U[i], T.F[i]); f << line; }
    }
  }
  {
    auto f = open("angles.table");
    f << head << "# θ (degrees), E, −dE/dθ (kcal/mol per degree, as angle_style table reads it)\n";
    for (const auto& T : r.angles) {
      if (!T.sampled) continue;
      f << "\n" << T.key << "\nN " << T.x.size() << "\n\n";
      for (size_t i = 0; i < T.x.size(); ++i) { std::snprintf(line, sizeof line, "%zu %.3f %.6f %.6f\n", i + 1, T.x[i], T.U[i], T.F[i]); f << line; }
    }
  }
  {
    auto f = open("dihedrals.table");
    f << head << "# φ (degrees, IUPAC: trans = 180), E; forces from LAMMPS's cyclic spline (NOF)\n";
    for (const auto& T : r.dihedrals) {
      if (!T.sampled) continue;
      f << "\n" << T.key << "\nN " << T.x.size() << " DEGREES NOF\n\n";
      for (size_t i = 0; i < T.x.size(); ++i) { std::snprintf(line, sizeof line, "%zu %.3f %.6f\n", i + 1, T.x[i], T.U[i]); f << line; }
    }
  }
  {
    auto f = open("bonded.in");
    f << "# CAPS coarse-grained bonded potentials (" << int(std::lround(r.temperature)) << " K): include after read_data; type numbers follow the shared type list\n";
    f << "# variable BONDED is the folder of the tables (e.g. -var BONDED bonded)\n";
    size_t nb = 0, na = 0, nd = 0;
    for (const auto& T : r.bonds) nb += T.sampled;
    for (const auto& T : r.angles) na += T.sampled;
    for (const auto& T : r.dihedrals) nd += T.sampled;
    size_t bmax = 0, amax = 0;
    for (const auto& T : r.bonds) bmax = std::max(bmax, T.x.size());
    for (const auto& T : r.angles) amax = std::max(amax, T.x.size());
    if (nb) f << "bond_style     table linear " << bmax << "\n";
    if (na) f << "angle_style    table linear " << amax << "\n";
    if (nd) f << "dihedral_style table/cut spline " << std::max<size_t>(360, size_t(std::lround(360.0 / o.dihedral_table_dp))) << "\n";
    auto num = [](const std::vector<std::string>& list, const std::string& k) { return size_t(std::find(list.begin(), list.end(), k) - list.begin()) + 1; };
    for (const auto& T : r.bonds) if (T.sampled) f << "bond_coeff     " << num(types.bonds, T.key) << " ${BONDED}/bonds.table " << T.key << "\n";
    for (const auto& T : r.angles) if (T.sampled) f << "angle_coeff    " << num(types.angles, T.key) << " ${BONDED}/angles.table " << T.key << "\n";
    for (const auto& T : r.dihedrals)
      if (T.sampled) {
        std::snprintf(line, sizeof line, "dihedral_coeff %zu aat %g %g %g ${BONDED}/dihedrals.table %s\n", num(types.dihedrals, T.key), o.aat_k, o.aat_theta1, o.aat_theta2, T.key.c_str());
        f << line;
      }
    // terms with no samples: said, not silently left out
    for (const auto* list : {&r.bonds, &r.angles, &r.dihedrals})
      for (const auto& T : *list) if (!T.sampled) f << "# no samples for " << kind_name(T.kind) << " " << T.key << ": give it a potential before a run\n";
  }
  // GROMACS: x in nm (bonds) or degrees (angles, dihedrals), kJ/mol; derivatives per nm or per degree
  constexpr double kJ = 4.184;
  auto xvg = [&](const std::string& name, const CgBondedTable& T, double xs, double fs) {
    auto f = open("gromacs/" + name);
    f << "# CAPS tabulated " << kind_name(T.kind) << " " << T.key << " (x, V kJ/mol, -dV/dx)\n";
    for (size_t i = 0; i < T.x.size(); ++i) { std::snprintf(line, sizeof line, "%.6f %.6e %.6e\n", T.x[i] * xs, T.U[i] * kJ, T.F[i] * fs); f << line; }
  };
  for (const auto& T : r.bonds) if (T.sampled) xvg("table_b" + std::to_string(size_t(std::find(types.bonds.begin(), types.bonds.end(), T.key) - types.bonds.begin())) + ".xvg", T, 0.1, kJ * 10);
  for (const auto& T : r.angles) if (T.sampled) xvg("table_a" + std::to_string(size_t(std::find(types.angles.begin(), types.angles.end(), T.key) - types.angles.begin())) + ".xvg", T, 1, kJ);
  for (const auto& T : r.dihedrals)
    if (T.sampled) {
      // GROMACS dihedral tables run −180 … 180 inclusive
      CgBondedTable G = T;
      G.x.push_back(180.0), G.U.push_back(T.U.front()), G.F.push_back(T.F.front());
      xvg("table_d" + std::to_string(size_t(std::find(types.dihedrals.begin(), types.dihedrals.end(), T.key) - types.dihedrals.begin())) + ".xvg", G, 1, kJ);
    }
  {
    auto f = open("gromacs/README.txt");
    f << "Tabulated bonded interactions for GROMACS (bond type 8, angle type 8, dihedral type 8 with table number n = the file's\n"
         "number; mdrun -tableb table_b0.xvg …). The GROMACS manual describes this format but notes that recent versions have\n"
         "tabulated bonded interactions disabled: check your version before relying on them. LAMMPS reads the tables in the\n"
         "folder above.\n";
  }
  {
    auto f = open("bonded.json");
    f << bonded_json(r).dump(0) << "\n";
  }
  return files;
}

}  // namespace caps
