// CAPS whole-pattern refinement (see refine.hpp).
#include "caps/refine.hpp"

#include "caps/measured.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <fstream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>

#include "caps/crystal.hpp"
#include "caps/crystalfit.hpp"
#include "caps/powder.hpp"
#include "caps/properties.hpp"
#include "caps/spacegroup.hpp"

namespace caps {

namespace {
constexpr double kPi = 3.14159265358979323846;
constexpr double kDeg = 180.0 / kPi;
using cd = std::complex<double>;

std::string fmt(double v) {
  char b[48];
  std::snprintf(b, sizeof b, "%.6g", v);
  return b;
}

// the lattice parameters as free variables: equal lengths share one, 90° and 120° angles stay, equal free angles share one
struct LatticeMap {
  std::array<int, 6> var{};    // per parameter: its variable, −1 fixed
  std::array<double, 6> fixed{};
  std::vector<double> start;
  std::array<double, 6> params(const std::vector<double>& v) const {
    std::array<double, 6> p{};
    for (int k = 0; k < 6; ++k) p[size_t(k)] = var[size_t(k)] < 0 ? fixed[size_t(k)] : v[size_t(var[size_t(k)])];
    return p;
  }
};
LatticeMap lattice_map(const std::array<double, 6>& p0) {
  LatticeMap m;
  for (int k = 0; k < 6; ++k) m.var[size_t(k)] = -2, m.fixed[size_t(k)] = p0[size_t(k)];
  for (int k = 0; k < 3; ++k) {
    if (m.var[size_t(k)] != -2) continue;
    m.var[size_t(k)] = int(m.start.size());
    m.start.push_back(p0[size_t(k)]);
    for (int q = k + 1; q < 3; ++q)
      if (std::fabs(p0[size_t(q)] - p0[size_t(k)]) < 1e-4 * p0[size_t(k)]) m.var[size_t(q)] = m.var[size_t(k)];
  }
  for (int k = 3; k < 6; ++k) {
    if (m.var[size_t(k)] != -2) continue;
    const double a = p0[size_t(k)];
    if (std::fabs(a - 90) < 1e-3 || std::fabs(a - 120) < 1e-3 || std::fabs(a - 60) < 1e-3) { m.var[size_t(k)] = -1; continue; }
    m.var[size_t(k)] = int(m.start.size());
    m.start.push_back(a);
    for (int q = k + 1; q < 6; ++q)
      if (std::fabs(p0[size_t(q)] - a) < 1e-3) m.var[size_t(q)] = m.var[size_t(k)];
  }
  return m;
}

// one reflection group: symmetry- (or metric-) equivalent vectors of one |G| at the starting lattice
struct Group {
  std::array<int, 3> hkl{};
  std::vector<std::array<int, 3>> members;   // half space; each counts twice (±)
  std::vector<std::vector<cd>> S;            // per member, per element: Σ_j e^{2πi h·f_j}
};

struct PhaseModel {
  const RefinePhase* src = nullptr;
  std::vector<int> elems;
  std::vector<Group> groups;
  LatticeMap lat;
  double mass = 0;
};

std::array<Vec3, 3> reciprocal(const Cell& c) {
  const double V = c.volume();
  return {cross(c.b, c.c) * (2 * kPi / V), cross(c.c, c.a) * (2 * kPi / V), cross(c.a, c.b) * (2 * kPi / V)};
}
Vec3 gv(const std::array<Vec3, 3>& B, const std::array<int, 3>& h) { return B[0] * double(h[0]) + B[1] * double(h[1]) + B[2] * double(h[2]); }

PhaseModel build_phase(const RefinePhase& ph, double qmax, bool neutron) {
  PhaseModel m;
  m.src = &ph;
  const System& s = ph.cell;
  if (!s.cell.valid() || s.atoms.empty()) throw std::invalid_argument("refine: phase " + ph.name + " needs a periodic cell with atoms");
  std::vector<int> eo(s.atoms.size());
  std::vector<Vec3> f(s.atoms.size());
  for (size_t j = 0; j < s.atoms.size(); ++j) {
    const int z = s.atoms[j].element;
    auto it = std::find(m.elems.begin(), m.elems.end(), z);
    if (it == m.elems.end()) m.elems.push_back(z), it = m.elems.end() - 1;
    eo[j] = int(it - m.elems.begin());
    f[j] = s.cell.to_fractional(s.atoms[j].pos);
    m.mass += s.mass_of(s.atoms[j]);
  }
  for (int z : m.elems) (void)scatter_w(neutron ? "neutron" : "xray", z, 1.0);
  m.lat = lattice_map(cell_parameters(s.cell));
  const auto B = reciprocal(s.cell);
  const int hm = int(qmax * norm(s.cell.a) / (2 * kPi)) + 1, km = int(qmax * norm(s.cell.b) / (2 * kPi)) + 1, lm = int(qmax * norm(s.cell.c) / (2 * kPi)) + 1;
  struct V { double q; std::array<int, 3> h; std::vector<cd> S; double f2; };
  std::vector<V> vs;
  double f2max = 0;
  for (int h = 0; h <= hm; ++h)
    for (int k = (h == 0 ? 0 : -km); k <= km; ++k)
      for (int l = (h == 0 && k == 0 ? 1 : -lm); l <= lm; ++l) {
        const std::array<int, 3> hk{h, k, l};
        const double q = norm(gv(B, hk));
        if (q > qmax) continue;
        V v{q, hk, std::vector<cd>(m.elems.size(), 0.0), 0};
        for (size_t j = 0; j < s.atoms.size(); ++j) v.S[size_t(eo[j])] += std::polar(1.0, 2 * kPi * (h * f[j][0] + k * f[j][1] + l * f[j][2]));
        cd F = 0;
        for (size_t e = 0; e < m.elems.size(); ++e) F += scatter_w(neutron ? "neutron" : "xray", m.elems[e], q) * v.S[e];
        v.f2 = std::norm(F);
        f2max = std::max(f2max, v.f2);
        vs.push_back(std::move(v));
      }
  std::sort(vs.begin(), vs.end(), [](const V& a, const V& b) { return a.q < b.q; });
  for (auto& v : vs) {
    if (v.f2 <= 1e-10 * f2max) continue;   // extinguished (an absence of the structure)
    if (m.groups.empty() || v.q - norm(gv(B, m.groups.back().members[0])) > 1e-6 * v.q + 1e-9) {
      Group g;
      g.hkl = v.h;
      m.groups.push_back(g);
    }
    Group& g = m.groups.back();
    g.members.push_back(v.h), g.S.push_back(v.S);
    if (v.h > g.hkl) g.hkl = v.h;
  }
  return m;
}

// everything that moves between trials
struct Params {
  double U, V, W, X, Y, zero;
  std::vector<std::vector<double>> lat;   // per phase: its lattice variables
  std::vector<double> b, r;               // per phase
};

struct Engine {
  const RefineOptions& o;
  std::vector<double> tt;                 // the fitted points
  std::vector<double> lam, wt;            // lines
  bool neutron = false;
  std::vector<PhaseModel> ph;
  double tmin = 0, tmax = 0;

  std::pair<double, double> shape(double th, const Params& p) const {
    const double tn = std::tan(th), cs = std::cos(th);
    const double hg2 = std::max(0.0, p.U * tn * tn + p.V * tn + p.W);
    const double hg = std::sqrt(hg2), hl = std::max(0.0, p.X * tn + p.Y / cs);
    double H = std::pow(std::pow(hg, 5) + 2.69269 * std::pow(hg, 4) * hl + 2.42843 * std::pow(hg, 3) * hl * hl + 4.47163 * hg * hg * std::pow(hl, 3) +
                            0.07842 * hg * std::pow(hl, 4) + std::pow(hl, 5), 0.2);
    const double r = H > 0 ? hl / H : 0;
    const double eta = std::clamp(1.36603 * r - 0.47719 * r * r + 0.11116 * r * r * r, 0.0, 1.0);
    return {std::max(H, 1e-6), eta};
  }
  double lp(double th) const {
    const double s = std::sin(th), c = std::cos(th), c2t = std::cos(2 * th);
    const double lor = 1 / (s * s * c);
    if (neutron || o.polarised) return lor;
    const double cm = std::cos(o.monochromator_2theta / kDeg), cm2 = cm * cm;
    return lor * (1 + cm2 * c2t * c2t) / (1 + cm2);
  }
  // adds intensity I of a peak at 2θ c (degrees) with (H, η) into col
  void peak(std::vector<double>& col, double c, double H, double eta, double I) const {
    const double c4l2 = 4 * std::log(2.0), gn = 2 * std::sqrt(std::log(2.0) / kPi), half = 40 * H;
    auto lo = std::lower_bound(tt.begin(), tt.end(), c - half), hi = std::upper_bound(tt.begin(), tt.end(), c + half);
    for (auto it = lo; it != hi; ++it) {
      const double x = *it - c, u = x * x / (H * H);
      col[size_t(it - tt.begin())] += I * (eta * 2 / (kPi * H) / (1 + 4 * u) + (1 - eta) * gn / H * std::exp(-c4l2 * u));
    }
  }
  // per phase: the pattern of each group at unit intensity (Pawley) or of the whole phase at unit scale (Rietveld)
  void columns(size_t pi, const Params& p, bool per_group, std::vector<std::vector<double>>& cols, std::vector<double>* pos = nullptr,
               std::vector<double>* inten = nullptr) const {
    const PhaseModel& m = ph[pi];
    const auto lp6 = m.lat.params(p.lat[pi]);
    const Cell cell = cell_parameters(lp6[0], lp6[1], lp6[2], lp6[3], lp6[4], lp6[5]);
    const auto B = reciprocal(cell);
    const bool po = (m.src->preferred[0] || m.src->preferred[1] || m.src->preferred[2]);
    Vec3 pref{0, 0, 1};
    if (po) pref = gv(B, m.src->preferred), pref = pref * (1 / norm(pref));
    if (!per_group) cols.push_back(std::vector<double>(tt.size(), 0.0));
    for (const auto& g : m.groups) {
      const double q = norm(gv(B, g.members[0]));
      double I = 1;
      if (!per_group) {
        const double s2 = (q / (4 * kPi)) * (q / (4 * kPi)), dw = std::exp(-p.b[pi] * s2);
        std::vector<double> w(m.elems.size());
        for (size_t e = 0; e < m.elems.size(); ++e) w[e] = scatter_w(neutron ? "neutron" : "xray", m.elems[e], q) * dw;
        I = 0;
        for (size_t k = 0; k < g.members.size(); ++k) {
          cd F = 0;
          for (size_t e = 0; e < m.elems.size(); ++e) F += w[e] * g.S[k][e];
          double pmd = 1;
          if (po) {
            const Vec3 G = gv(B, g.members[k]);
            const double c = dot(G, pref) / std::max(1e-12, norm(G)), c2 = c * c, r = p.r[pi];
            pmd = std::pow(r * r * c2 + (1 - c2) / r, -1.5);
          }
          I += 2 * std::norm(F) * pmd;
        }
      }
      if (per_group) cols.push_back(std::vector<double>(tt.size(), 0.0));
      std::vector<double>& col = cols.back();
      double first_tt = 0, first_I = 0;
      for (size_t li = 0; li < lam.size(); ++li) {
        const double sn = lam[li] * q / (4 * kPi);
        if (sn >= 1) continue;
        const double th = std::asin(sn), c = 2 * th * kDeg + p.zero;
        const auto [H, eta] = shape(th, p);
        const double In = I * wt[li] * lp(th);
        if (li == 0) first_tt = 2 * th * kDeg, first_I = In;
        if (c < tmin - 40 * H || c > tmax + 40 * H) continue;
        peak(col, c, H, eta, In);
      }
      if (pos) pos->push_back(first_tt);
      if (inten) inten->push_back(first_I);
    }
  }
  std::vector<std::vector<double>> background() const {
    std::vector<std::vector<double>> c(size_t(std::max(0, o.background)), std::vector<double>(tt.size()));
    for (size_t i = 0; i < tt.size(); ++i) {
      const double x = tmax > tmin ? 2 * (tt[i] - tmin) / (tmax - tmin) - 1 : 0;
      double t0 = 1, t1 = x;
      for (size_t k = 0; k < c.size(); ++k) {
        c[k][i] = k == 0 ? 1 : k == 1 ? x : 2 * x * t1 - t0;
        if (k >= 2) t0 = t1, t1 = c[k][i];
      }
    }
    return c;
  }
};

// weighted least squares on columns; nonneg[k]: column k must have a non-negative coefficient (an active-set loop)
std::vector<double> linear_fit(const std::vector<std::vector<double>>& cols, const std::vector<double>& y, const std::vector<double>& w,
                               const std::vector<char>& nonneg) {
  const size_t m = cols.size();
  std::vector<char> on(m, 1);
  std::vector<double> c(m, 0.0);
  for (int it = 0; it < 100; ++it) {
    std::vector<size_t> idx;
    for (size_t k = 0; k < m; ++k) if (on[k]) idx.push_back(k);
    const size_t n = idx.size();
    std::vector<double> A(n * n, 0.0), b(n, 0.0);
    for (size_t a = 0; a < n; ++a) {
      const auto& ca = cols[idx[a]];
      for (size_t i = 0; i < y.size(); ++i) b[a] += w[i] * ca[i] * y[i];
      for (size_t bb = a; bb < n; ++bb) {
        const auto& cb = cols[idx[bb]];
        double s = 0;
        for (size_t i = 0; i < y.size(); ++i) s += w[i] * ca[i] * cb[i];
        A[a * n + bb] = A[bb * n + a] = s;
      }
    }
    // Cholesky with a tiny ridge for empty columns
    double tr = 0;
    for (size_t a = 0; a < n; ++a) tr = std::max(tr, A[a * n + a]);
    for (size_t a = 0; a < n; ++a) A[a * n + a] += 1e-12 * tr + 1e-300;
    std::vector<double> L(n * n, 0.0);
    for (size_t i = 0; i < n; ++i)
      for (size_t j = 0; j <= i; ++j) {
        double s = A[i * n + j];
        for (size_t k = 0; k < j; ++k) s -= L[i * n + k] * L[j * n + k];
        L[i * n + j] = i == j ? std::sqrt(std::max(s, 1e-300)) : s / L[j * n + j];
      }
    std::vector<double> z(n), x(n);
    for (size_t i = 0; i < n; ++i) {
      double s = b[i];
      for (size_t k = 0; k < i; ++k) s -= L[i * n + k] * z[k];
      z[i] = s / L[i * n + i];
    }
    for (size_t i = n; i-- > 0;) {
      double s = z[i];
      for (size_t k = i + 1; k < n; ++k) s -= L[k * n + i] * x[k];
      x[i] = s / L[i * n + i];
    }
    std::fill(c.begin(), c.end(), 0.0);
    bool neg = false;
    for (size_t a = 0; a < n; ++a) {
      c[idx[a]] = x[a];
      if (nonneg[idx[a]] && x[a] < 0) on[idx[a]] = 0, neg = true;
    }
    if (!neg) break;
  }
  return c;
}
}  // namespace

MeasuredPattern read_pattern(const std::string& path) {
  // columns, XRDML, GSAS, UXD, Rigaku or JCAMP-DX (caps/measured.hpp)
  const MeasuredCurve c = read_measured(path);
  MeasuredPattern p;
  p.two_theta = c.x, p.y = c.y, p.sigma = c.sigma;
  if (p.two_theta.size() < 10) throw std::invalid_argument(path + ": fewer than ten points (2θ intensity [σ])");
  for (size_t i = 1; i < p.two_theta.size(); ++i)
    if (p.two_theta[i] <= p.two_theta[i - 1]) throw std::invalid_argument(path + ": 2θ must increase");
  return p;
}

namespace {
Engine make_engine(const std::vector<double>& two_theta, const std::vector<RefinePhase>& phases, const RefineOptions& o) {
  Engine e{o, two_theta, {}, {}, o.radiation == "neutron", {}, two_theta.front(), two_theta.back()};
  if (o.radiation != "xray" && o.radiation != "neutron") throw std::invalid_argument("refine: radiation xray or neutron");
  if (o.wavelength > 0) e.lam = {o.wavelength}, e.wt = {1};
  else {
    const auto ka = anode_kalpha(o.anode);
    e.lam = {ka[0]}, e.wt = {1};
    // a neutron source is monochromatic: no K-alpha2 line (as in powder.cpp)
    if (o.kalpha2 && !e.neutron) e.lam.push_back(ka[1]), e.wt.push_back(o.ka2_ratio);
  }
  const double lam_min = *std::min_element(e.lam.begin(), e.lam.end());
  const double th_max = std::min(89.999, 0.5 * (e.tmax + 2)) / kDeg;
  const double qmax = 4 * kPi * std::sin(th_max) / lam_min * 1.03;
  for (const auto& p : phases) e.ph.push_back(build_phase(p, qmax, e.neutron));
  return e;
}
}  // namespace

std::vector<double> calculated_pattern(const std::vector<double>& two_theta, const std::vector<RefinePhase>& phases, const std::vector<double>& scales,
                                       const RefineOptions& o) {
  if (two_theta.size() < 2) throw std::invalid_argument("refine: no points");
  const Engine e = make_engine(two_theta, phases, o);
  Params p{o.U, o.V, o.W, o.X, o.Y, o.zero, {}, {}, {}};
  for (size_t k = 0; k < phases.size(); ++k) p.lat.push_back(e.ph[k].lat.start), p.b.push_back(phases[k].b), p.r.push_back(phases[k].march_r);
  std::vector<double> y(two_theta.size(), 0.0);
  for (size_t k = 0; k < phases.size(); ++k) {
    std::vector<std::vector<double>> c;
    e.columns(k, p, false, c);
    for (size_t i = 0; i < y.size(); ++i) y[i] += (k < scales.size() ? scales[k] : 1.0) * c[0][i];
  }
  return y;
}

RefineResult refine_pattern(const MeasuredPattern& data0, const std::vector<RefinePhase>& phases, const RefineOptions& o) {
  RefineResult res;
  if (phases.empty()) throw std::invalid_argument("refine: no phases");
  const bool pawley = o.mode == "pawley";
  if (!pawley && o.mode != "rietveld") throw std::invalid_argument("refine: mode rietveld or pawley");
  if (pawley && phases.size() != 1) throw std::invalid_argument("refine: Pawley fits one phase");
  // the range
  MeasuredPattern d;
  for (size_t i = 0; i < data0.two_theta.size(); ++i) {
    const double t = data0.two_theta[i];
    if ((o.two_theta_min > 0 && t < o.two_theta_min) || (o.two_theta_max > 0 && t > o.two_theta_max)) continue;
    d.two_theta.push_back(t), d.y.push_back(data0.y[i]);
    if (!data0.sigma.empty()) d.sigma.push_back(data0.sigma[i]);
  }
  if (d.two_theta.size() < 10) throw std::invalid_argument("refine: fewer than ten points in the range");
  std::vector<double> w(d.y.size());
  for (size_t i = 0; i < w.size(); ++i) {
    const double s2 = !d.sigma.empty() ? d.sigma[i] * d.sigma[i] : std::max(1.0, d.y[i]);
    w[i] = 1 / std::max(s2, 1e-300);
  }
  Engine e = make_engine(d.two_theta, phases, o);
  // the nonlinear parameters
  struct PRef { std::string name; double* (*get)(Params&, size_t, size_t); size_t a, b; double lo, hi; };
  Params p{o.U, o.V, o.W, o.X, o.Y, o.zero, {}, {}, {}};
  for (size_t k = 0; k < phases.size(); ++k) p.lat.push_back(e.ph[k].lat.start), p.b.push_back(phases[k].b), p.r.push_back(phases[k].march_r);
  std::vector<PRef> refs;
  const double inf = std::numeric_limits<double>::infinity();
  if (o.refine_zero) refs.push_back({"zero", [](Params& q, size_t, size_t) { return &q.zero; }, 0, 0, -2.0, 2.0});
  // U alone is kept ≥ 0: a negative U with V fixed takes H_G² below zero at high angles, where the width is clipped and the
  // fit loses its gradient (Caglioti's U < 0 needs V refined with it)
  if (o.refine_U) refs.push_back({"U", [](Params& q, size_t, size_t) { return &q.U; }, 0, 0, o.refine_V ? -inf : 0.0, inf});
  if (o.refine_V) refs.push_back({"V", [](Params& q, size_t, size_t) { return &q.V; }, 0, 0, -inf, inf});
  if (o.refine_W) refs.push_back({"W", [](Params& q, size_t, size_t) { return &q.W; }, 0, 0, 0.0, inf});
  if (o.refine_X) refs.push_back({"X", [](Params& q, size_t, size_t) { return &q.X; }, 0, 0, 0.0, inf});
  if (o.refine_Y) refs.push_back({"Y", [](Params& q, size_t, size_t) { return &q.Y; }, 0, 0, 0.0, inf});
  for (size_t k = 0; k < phases.size(); ++k) {
    const std::string nm = phases[k].name.empty() ? "phase " + std::to_string(k + 1) : phases[k].name;
    if (phases[k].lattice)
      for (size_t v = 0; v < p.lat[k].size(); ++v) {
        int which = 0;
        for (int q = 0; q < 6; ++q) if (e.ph[k].lat.var[size_t(q)] == int(v)) { which = q; break; }
        static const char* pn[6] = {"a", "b", "c", "alpha", "beta", "gamma"};
        refs.push_back({nm + " " + pn[which], [](Params& q, size_t a, size_t b) { return &q.lat[a][b]; }, k, v, which < 3 ? 0.5 : 30.0, which < 3 ? 1e4 : 150.0});
      }
    if (phases[k].b_iso && !pawley) refs.push_back({nm + " B", [](Params& q, size_t a, size_t) { return &q.b[a]; }, k, 0, -2.0, 50.0});
    if (phases[k].march && !pawley) refs.push_back({nm + " March r", [](Params& q, size_t a, size_t) { return &q.r[a]; }, k, 0, 0.05, 20.0});
  }
  const size_t np = refs.size();
  const auto bg = e.background();
  // the linear solve at a trial: residuals, the coefficients
  auto evaluate = [&](Params& q, std::vector<double>& coef, std::vector<double>* ycalc, std::vector<std::vector<double>>* phase_cols) {
    std::vector<std::vector<double>> cols = bg;
    std::vector<char> nonneg(bg.size(), 0);
    std::vector<size_t> start;
    for (size_t k = 0; k < e.ph.size(); ++k) {
      start.push_back(cols.size());
      e.columns(k, q, pawley, cols);
    }
    nonneg.resize(cols.size(), 1);
    coef = linear_fit(cols, d.y, w, nonneg);
    std::vector<double> r(d.y.size()), yc(d.y.size(), 0.0);
    for (size_t c = 0; c < cols.size(); ++c)
      if (coef[c] != 0)
        for (size_t i = 0; i < yc.size(); ++i) yc[i] += coef[c] * cols[c][i];
    for (size_t i = 0; i < yc.size(); ++i) r[i] = std::sqrt(w[i]) * (d.y[i] - yc[i]);
    if (ycalc) *ycalc = yc;
    if (phase_cols) *phase_cols = cols;
    return r;
  };
  std::vector<double> coef;
  std::vector<double> res0 = evaluate(p, coef, nullptr, nullptr);
  auto ss_of = [](const std::vector<double>& v) { double s = 0; for (double x : v) s += x * x; return s; };
  double ss = ss_of(res0);
  // Levenberg–Marquardt on the nonlinear parameters (scaled by their starts)
  std::vector<double> scale(np);
  for (size_t k = 0; k < np; ++k) {
    const double v = *refs[k].get(p, refs[k].a, refs[k].b);
    scale[k] = std::max(std::fabs(v), refs[k].name == "zero" ? 0.01 : refs[k].name == "U" || refs[k].name == "V" ? 0.001 : 1e-3);
  }
  double lambda = 1e-3;
  int it = 0;
  std::vector<double> J;
  bool done = np == 0;
  for (; it < o.iterations && !done; ++it) {
    if (o.progress && !o.progress("iteration " + std::to_string(it + 1), double(it) / o.iterations)) throw std::runtime_error("refine stopped");
    const size_t m = res0.size();
    J.assign(m * np, 0.0);
    for (size_t k = 0; k < np; ++k) {
      Params q = p;
      double* v = refs[k].get(q, refs[k].a, refs[k].b);
      double h = 1e-6 * std::max(std::fabs(*v), scale[k]);
      if (*v + h > refs[k].hi) h = -h;
      *v += h;
      std::vector<double> cf;
      const auto rk = evaluate(q, cf, nullptr, nullptr);
      for (size_t i = 0; i < m; ++i) J[i * np + k] = (rk[i] - res0[i]) / h * scale[k];
    }
    std::vector<double> A(np * np, 0.0), g(np, 0.0);
    for (size_t i = 0; i < m; ++i)
      for (size_t a = 0; a < np; ++a) {
        g[a] -= J[i * np + a] * res0[i];
        for (size_t b = 0; b < np; ++b) A[a * np + b] += J[i * np + a] * J[i * np + b];
      }
    bool accepted = false;
    for (int tries = 0; tries < 12 && !accepted; ++tries) {
      std::vector<double> M = A, dlt = g;
      for (size_t a = 0; a < np; ++a) M[a * np + a] += lambda * std::max(A[a * np + a], 1e-12);
      // solve M dlt = g (Gauss–Jordan)
      bool sing = false;
      for (size_t c = 0; c < np && !sing; ++c) {
        size_t piv = c;
        for (size_t i = c + 1; i < np; ++i) if (std::fabs(M[i * np + c]) > std::fabs(M[piv * np + c])) piv = i;
        if (std::fabs(M[piv * np + c]) < 1e-300) { sing = true; break; }
        for (size_t j = 0; j < np; ++j) std::swap(M[c * np + j], M[piv * np + j]);
        std::swap(dlt[c], dlt[piv]);
        const double dd = M[c * np + c];
        for (size_t j = 0; j < np; ++j) M[c * np + j] /= dd;
        dlt[c] /= dd;
        for (size_t i = 0; i < np; ++i)
          if (i != c) {
            const double f = M[i * np + c];
            for (size_t j = 0; j < np; ++j) M[i * np + j] -= f * M[c * np + j];
            dlt[i] -= f * dlt[c];
          }
      }
      if (sing) { lambda *= 10; continue; }
      Params q = p;
      double dmax = 0;
      for (size_t k = 0; k < np; ++k) {
        double* v = refs[k].get(q, refs[k].a, refs[k].b);
        const double nv = std::clamp(*v + dlt[k] * scale[k], refs[k].lo, refs[k].hi);
        dmax = std::max(dmax, std::fabs(nv - *v) / scale[k]);
        *v = nv;
      }
      if (dmax < 1e-12) { done = true; break; }
      std::vector<double> cf;
      const auto rn = evaluate(q, cf, nullptr, nullptr);
      const double sn = ss_of(rn);
      if (sn < ss) {
        const double gain = (ss - sn) / ss;
        p = q, res0 = rn, ss = sn, coef = cf, accepted = true;
        lambda = std::max(lambda / 10, 1e-12);
        if (gain < 1e-9) done = true;
      } else {
        lambda *= 10;
      }
    }
    if (!accepted) done = true;
  }
  res.iterations = it;
  // esd of the nonlinear parameters: s² (JᵀJ)⁻¹ at the end (J of the last iteration)
  std::vector<double> esd(np, std::numeric_limits<double>::quiet_NaN());
  const size_t nlin = coef.size();
  const double dof = double(res0.size()) - double(np) - double(nlin);
  if (np > 0 && dof > 0 && !J.empty()) {
    std::vector<double> A(np * np, 0.0);
    for (size_t i = 0; i < res0.size(); ++i)
      for (size_t a = 0; a < np; ++a)
        for (size_t b = 0; b < np; ++b) A[a * np + b] += J[i * np + a] * J[i * np + b];
    // invert by Gauss–Jordan
    std::vector<double> Inv(np * np, 0.0);
    for (size_t a = 0; a < np; ++a) Inv[a * np + a] = 1;
    bool ok = true;
    for (size_t c = 0; c < np && ok; ++c) {
      size_t piv = c;
      for (size_t i = c + 1; i < np; ++i) if (std::fabs(A[i * np + c]) > std::fabs(A[piv * np + c])) piv = i;
      if (std::fabs(A[piv * np + c]) < 1e-300) { ok = false; break; }
      for (size_t j = 0; j < np; ++j) std::swap(A[c * np + j], A[piv * np + j]), std::swap(Inv[c * np + j], Inv[piv * np + j]);
      const double dd = A[c * np + c];
      for (size_t j = 0; j < np; ++j) A[c * np + j] /= dd, Inv[c * np + j] /= dd;
      for (size_t i = 0; i < np; ++i)
        if (i != c) {
          const double f = A[i * np + c];
          for (size_t j = 0; j < np; ++j) A[i * np + j] -= f * A[c * np + j], Inv[i * np + j] -= f * Inv[c * np + j];
        }
    }
    if (ok)
      for (size_t k = 0; k < np; ++k) esd[k] = std::sqrt(std::max(0.0, ss / dof * Inv[k * np + k])) * scale[k];
  }
  // results
  std::vector<double> yc;
  std::vector<std::vector<double>> cols;
  evaluate(p, coef, &yc, &cols);
  res.two_theta = d.two_theta, res.y_obs = d.y, res.y_calc = yc;
  res.background.assign(d.y.size(), 0.0);
  for (size_t c = 0; c < bg.size(); ++c)
    for (size_t i = 0; i < d.y.size(); ++i) res.background[i] += coef[c] * bg[c][i];
  res.U = p.U, res.V = p.V, res.W = p.W, res.X = p.X, res.Y = p.Y, res.zero = p.zero;
  const double lam1 = e.lam[0];
  res.size = p.Y > 1e-5 ? 0.9 * lam1 / (p.Y / kDeg) : 0;   // 0: no size broadening found
  res.strain = p.X / kDeg / 4;
  double sumw = 0, sumabs = 0, sumy = 0;
  for (size_t i = 0; i < d.y.size(); ++i) sumw += w[i] * d.y[i] * d.y[i], sumabs += std::fabs(d.y[i] - yc[i]), sumy += std::fabs(d.y[i]);
  res.parameters = int(np + nlin);
  res.rp = sumy > 0 ? sumabs / sumy : 0;
  res.rwp = std::sqrt(ss / sumw);
  res.rexp = std::sqrt(std::max(0.0, dof) / sumw);
  res.chi2 = res.rexp > 0 ? (res.rwp / res.rexp) * (res.rwp / res.rexp) : 0;
  for (size_t k = 0; k < np; ++k) res.refined.push_back({refs[k].name, *refs[k].get(p, refs[k].a, refs[k].b)});
  // phases
  double zmv = 0;
  size_t col = bg.size();
  for (size_t k = 0; k < phases.size(); ++k) {
    RefinedPhase rp;
    rp.name = phases[k].name.empty() ? "phase " + std::to_string(k + 1) : phases[k].name;
    const auto l6 = e.ph[k].lat.params(p.lat[k]);
    rp.lattice = l6;
    for (int q = 0; q < 6; ++q) {
      rp.lattice_esd[size_t(q)] = 0;
      if (e.ph[k].lat.var[size_t(q)] < 0) continue;
      for (size_t r = 0; r < np; ++r)
        if (refs[r].a == k && refs[r].get(p, k, refs[r].b) == &p.lat[k][size_t(e.ph[k].lat.var[size_t(q)])]) rp.lattice_esd[size_t(q)] = esd[r];
    }
    rp.b = p.b[k], rp.march_r = p.r[k];
    const Cell cell = cell_parameters(l6[0], l6[1], l6[2], l6[3], l6[4], l6[5]);
    rp.volume = cell.volume(), rp.mass = e.ph[k].mass;
    std::vector<std::vector<double>> unit;
    std::vector<double> pos, inten;
    e.columns(k, p, false, unit, &pos, &inten);
    for (const auto& g : e.ph[k].groups) rp.hkl.push_back(g.hkl);
    rp.two_theta = pos;
    if (pawley) {
      rp.intensity.assign(e.ph[k].groups.size(), 0.0);
      for (size_t g = 0; g < e.ph[k].groups.size(); ++g) rp.intensity[g] = coef[col + g];
      col += e.ph[k].groups.size();
      rp.scale = 1;
    } else {
      rp.scale = coef[col++];
      rp.intensity = inten;
      for (auto& v : rp.intensity) v *= rp.scale;
      zmv += rp.scale * rp.mass * rp.volume;
    }
    res.phases.push_back(rp);
  }
  if (!pawley && zmv > 0)
    for (auto& rp : res.phases) rp.weight_fraction = rp.scale * rp.mass * rp.volume / zmv;
  res.notes.push_back(std::string(pawley ? "Pawley" : "Rietveld") + " fit of " + std::to_string(d.y.size()) + " points (" + fmt(d.two_theta.front()) + "–" +
                      fmt(d.two_theta.back()) + "° 2θ): " + std::to_string(np) + " nonlinear and " + std::to_string(nlin) + " linear parameters, " +
                      std::to_string(it) + " iterations");
  if (d.sigma.empty()) res.notes.push_back("weights 1/y (counting statistics): give a third column σ when the data are not counts");
  res.notes.push_back("atom positions and occupancies are not refined");
  return res;
}

Json refine_json(const RefineResult& r) {
  auto arr = [](const std::vector<double>& v) { Json a = Json::array(); for (double x : v) a.push_back(std::isfinite(x) ? Json(x) : Json(nullptr)); return a; };
  Json j = Json::object(), ph = Json::array(), rf = Json::array(), nt = Json::array();
  for (const auto& p : r.phases) {
    Json x = Json::object(), hk = Json::array();
    for (const auto& h : p.hkl) {
      Json t = Json::array();
      for (int q : h) t.push_back(Json(q));
      hk.push_back(t);
    }
    x["name"] = p.name, x["lattice"] = arr(std::vector<double>(p.lattice.begin(), p.lattice.end()));
    x["lattice_esd"] = arr(std::vector<double>(p.lattice_esd.begin(), p.lattice_esd.end()));
    x["scale"] = p.scale, x["b"] = p.b, x["march_r"] = p.march_r, x["weight_fraction"] = p.weight_fraction, x["volume"] = p.volume, x["mass"] = p.mass;
    x["hkl"] = hk, x["two_theta"] = arr(p.two_theta), x["intensity"] = arr(p.intensity);
    ph.push_back(x);
  }
  for (const auto& [n, v] : r.refined) {
    Json x = Json::object();
    x["name"] = n, x["value"] = v;
    rf.push_back(x);
  }
  for (const auto& n : r.notes) nt.push_back(n);
  j["two_theta"] = arr(r.two_theta), j["y_obs"] = arr(r.y_obs), j["y_calc"] = arr(r.y_calc), j["background"] = arr(r.background);
  j["phases"] = ph, j["refined"] = rf, j["U"] = r.U, j["V"] = r.V, j["W"] = r.W, j["X"] = r.X, j["Y"] = r.Y, j["zero"] = r.zero;
  j["size"] = r.size, j["strain"] = r.strain, j["rp"] = r.rp, j["rwp"] = r.rwp, j["rexp"] = r.rexp, j["chi2"] = r.chi2;
  j["iterations"] = double(r.iterations), j["parameters"] = double(r.parameters), j["notes"] = nt;
  return j;
}

}  // namespace caps
