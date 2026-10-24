// CAPS Exchange, analysis (see exchange_analysis.hpp).
#include "caps/exchange_analysis.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <numeric>
#include <stdexcept>

namespace caps {

namespace {

constexpr double NaN = std::numeric_limits<double>::quiet_NaN();
constexpr double kGasKJ = 8.314462618e-3;   // kJ/(mol K)
constexpr double kBoltzmann = 1.380649e-23; // J/K
constexpr double kAtmPa = 101325.0;
constexpr double kInvE = 0.36787944117144233;

std::string fmt(const char* f, double v) {
  char b[64];
  std::snprintf(b, sizeof b, f, v);
  return b;
}

// Least squares y = a + b x with standard errors of a and b and R².
struct Line { double a = NaN, b = NaN, sa = NaN, sb = NaN, r2 = NaN; };
Line line_fit(const std::vector<double>& x, const std::vector<double>& y) {
  Line l;
  const size_t n = x.size();
  if (n < 2) return l;
  const double mx = std::accumulate(x.begin(), x.end(), 0.0) / n, my = std::accumulate(y.begin(), y.end(), 0.0) / n;
  double sxx = 0, sxy = 0, syy = 0;
  for (size_t i = 0; i < n; ++i) { sxx += (x[i] - mx) * (x[i] - mx); sxy += (x[i] - mx) * (y[i] - my); syy += (y[i] - my) * (y[i] - my); }
  if (sxx <= 0) return l;
  l.b = sxy / sxx;
  l.a = my - l.b * mx;
  double sse = 0;
  for (size_t i = 0; i < n; ++i) { const double r = y[i] - l.a - l.b * x[i]; sse += r * r; }
  l.r2 = syy > 0 ? 1 - sse / syy : 1.0;
  if (n > 2) {
    const double s2 = sse / double(n - 2);
    l.sb = std::sqrt(s2 / sxx);
    l.sa = std::sqrt(s2 * (1.0 / n + mx * mx / sxx));
  }
  return l;
}

// Nelder–Mead on a function of n variables; returns the best point found.
template <class F>
std::vector<double> nelder_mead(F f, std::vector<double> x0, const std::vector<double>& step, int iters = 4000, double tol = 1e-14) {
  const size_t n = x0.size();
  std::vector<std::vector<double>> s(n + 1, x0);
  for (size_t i = 0; i < n; ++i) s[i + 1][i] += step[i];
  std::vector<double> v(n + 1);
  for (size_t i = 0; i <= n; ++i) v[i] = f(s[i]);
  for (int it = 0; it < iters; ++it) {
    std::vector<size_t> o(n + 1);
    std::iota(o.begin(), o.end(), 0);
    std::sort(o.begin(), o.end(), [&](size_t a, size_t b) { return v[a] < v[b]; });
    {
      const auto s2 = s;
      const auto v2 = v;
      for (size_t i = 0; i <= n; ++i) s[i] = s2[o[i]], v[i] = v2[o[i]];
    }
    if (std::fabs(v[n] - v[0]) <= tol * (std::fabs(v[0]) + tol)) break;
    std::vector<double> c(n, 0.0);
    for (size_t i = 0; i < n; ++i)
      for (size_t k = 0; k < n; ++k) c[k] += s[i][k] / double(n);
    auto along = [&](double t) {
      std::vector<double> p(n);
      for (size_t k = 0; k < n; ++k) p[k] = c[k] + t * (s[n][k] - c[k]);
      return p;
    };
    const auto xr = along(-1.0);
    const double fr = f(xr);
    if (fr < v[0]) {
      const auto xe = along(-2.0);
      const double fe = f(xe);
      if (fe < fr) s[n] = xe, v[n] = fe;
      else s[n] = xr, v[n] = fr;
    } else if (fr < v[n - 1]) {
      s[n] = xr, v[n] = fr;
    } else {
      const auto xc = fr < v[n] ? along(-0.5) : along(0.5);
      const double fc = f(xc);
      if (fc < std::min(fr, v[n])) {
        s[n] = xc, v[n] = fc;
      } else {
        for (size_t i = 1; i <= n; ++i) {
          for (size_t k = 0; k < n; ++k) s[i][k] = s[0][k] + 0.5 * (s[i][k] - s[0][k]);
          v[i] = f(s[i]);
        }
      }
    }
  }
  return s[std::min_element(v.begin(), v.end()) - v.begin()];
}

// The thickness of the cell along an axis (the distance between the faces the other two edges span), Å.
double thickness(const Cell& c, int axis) {
  const Vec3 n = axis == 0 ? cross(c.b, c.c) : axis == 1 ? cross(c.c, c.a) : cross(c.a, c.b);
  const double a = norm(n);
  return a > 0 ? c.volume() / a : 0;
}

}  // namespace

// ---------------------------------------------------------------- Arrhenius

ArrheniusFit fit_arrhenius(const std::vector<double>& T, const std::vector<double>& tau, bool reduced, double t_min, double t_max) {
  if (T.size() != tau.size()) throw std::invalid_argument("fit_arrhenius: as many relaxation times as temperatures are needed");
  ArrheniusFit f;
  f.reduced = reduced;
  std::vector<double> x, y;
  for (size_t i = 0; i < T.size(); ++i) {
    if (!(T[i] > 0) || !(tau[i] > 0) || !std::isfinite(tau[i])) continue;
    if (t_min > 0 && T[i] < t_min) continue;
    if (t_max > 0 && T[i] > t_max) continue;
    f.T.push_back(T[i]), f.tau.push_back(tau[i]);
    x.push_back(1.0 / T[i]), y.push_back(std::log(tau[i]));
  }
  f.points = x.size();
  if (f.points < 2 || *std::max_element(f.T.begin(), f.T.end()) == *std::min_element(f.T.begin(), f.T.end())) {
    f.note = "an Arrhenius fit needs relaxation times at two temperatures or more";
    return f;
  }
  const Line l = line_fit(x, y);
  const double r = reduced ? 1.0 : kGasKJ;
  f.ea = l.b * r, f.ea_err = l.sb * r;
  f.ln_tau0 = l.a, f.ln_tau0_err = l.sa;
  f.r2 = l.r2;
  f.ok = std::isfinite(f.ea);
  if (f.points == 2) f.note = "two temperatures: the line goes through both, so there is no error estimate";
  return f;
}

double crossing_temperature(const std::vector<double>& T, const std::vector<double>& tau_a, const std::vector<double>& tau_b) {
  if (T.size() != tau_a.size() || T.size() != tau_b.size()) throw std::invalid_argument("crossing_temperature: one τ_a and one τ_b per temperature");
  std::vector<std::pair<double, double>> d;   // (1/T, ln τ_a − ln τ_b)
  for (size_t i = 0; i < T.size(); ++i)
    if (T[i] > 0 && tau_a[i] > 0 && tau_b[i] > 0) d.push_back({1.0 / T[i], std::log(tau_a[i]) - std::log(tau_b[i])});
  std::sort(d.begin(), d.end(), [](const auto& p, const auto& q) { return p.first > q.first; });   // ascending T
  for (size_t k = 0; k + 1 < d.size(); ++k) {
    if (d[k].second == 0) return 1.0 / d[k].first;
    if (d[k].second * d[k + 1].second < 0) {
      const double x = d[k].first + d[k].second * (d[k + 1].first - d[k].first) / (d[k].second - d[k + 1].second);
      return 1.0 / x;
    }
  }
  if (!d.empty() && d.back().second == 0) return 1.0 / d.back().first;
  return NaN;
}

std::vector<Property> arrhenius_properties(const ArrheniusFit& f, const std::string& what) {
  Property p{"arrhenius", "Activation energy of " + what, f.reduced ? "ε (k_B = 1)" : "kJ/mol", "", NaN, NaN, {}, {}, {}};
  p.method = "Arrhenius fit ln " + what + " = ln " + what + "₀ + Ea/(" + std::string(f.reduced ? "k_B" : "R") +
             " T), least squares in 1/T; error: standard error of the slope";
  if (!f.note.empty()) p.notes.push_back(f.note);
  if (!f.ok) return {p};
  p.value = f.ea, p.error = f.ea_err;
  p.extra["ln " + what + "0"] = f.ln_tau0;
  if (std::isfinite(f.ln_tau0_err)) p.extra["ln " + what + "0 standard error"] = f.ln_tau0_err;
  p.extra["R²"] = f.r2;
  p.extra["points"] = double(f.points);
  if (!f.reduced) p.extra["Ea (kcal/mol)"] = f.ea / 4.184;
  const double xs = f.reduced ? 1.0 : 1000.0;
  Series d{"ln " + what, f.reduced ? "1/T" : "1000/T (1/K)", "ln " + what, {}, {}}, line{"fit", d.x_label, d.y_label, {}, {}};
  for (size_t i = 0; i < f.T.size(); ++i) d.x.push_back(xs / f.T[i]), d.y.push_back(std::log(f.tau[i]));
  const auto [lo, hi] = std::minmax_element(f.T.begin(), f.T.end());
  for (double t : {*lo, *hi}) line.x.push_back(xs / t), line.y.push_back(f.ln_tau0 + f.ea / ((f.reduced ? 1.0 : kGasKJ) * t));
  p.series.push_back(std::move(d));
  p.series.push_back(std::move(line));
  return {p};
}

// ---------------------------------------------------------------- multi-tau correlator

MultiTauCorrelator::MultiTauCorrelator(size_t channels, int p, int m, int levels)
    : nc_(channels), p_(std::max(2, p)), m_(std::max(2, m)), nlev_(levels > 0 ? levels : 40) {
  if (p_ % m_) throw std::invalid_argument("multi-tau correlator: p must be a multiple of m");
  buf_.assign(nlev_, std::vector<std::vector<double>>(p_, std::vector<double>(nc_, 0.0)));
  filled_.assign(nlev_, 0);
  sum_.assign(nlev_, std::vector<double>(p_, 0.0));
  count_.assign(nlev_, std::vector<double>(p_, 0.0));
  acc_.assign(nlev_, std::vector<double>(nc_, 0.0));
  nacc_.assign(nlev_, 0);
}

void MultiTauCorrelator::add(const std::vector<double>& x) {
  if (x.size() != nc_) throw std::invalid_argument("multi-tau correlator: wrong number of channels");
  add_at(0, x);
}

void MultiTauCorrelator::add_at(int l, const std::vector<double>& x) {
  if (l >= nlev_) return;
  auto& b = buf_[l];
  std::rotate(b.begin(), b.end() - 1, b.end());
  b[0] = x;
  filled_[l] = std::min(filled_[l] + 1, p_);
  for (int j = l == 0 ? 0 : p_ / m_; j < filled_[l]; ++j) {
    double s = 0;
    for (size_t c = 0; c < nc_; ++c) s += b[0][c] * b[j][c];
    sum_[l][j] += s;
    count_[l][j] += 1;
  }
  for (size_t c = 0; c < nc_; ++c) acc_[l][c] += x[c];
  if (++nacc_[l] == m_) {
    std::vector<double> avg(nc_);
    for (size_t c = 0; c < nc_; ++c) avg[c] = acc_[l][c] / m_, acc_[l][c] = 0;
    nacc_[l] = 0;
    add_at(l + 1, avg);
  }
}

void MultiTauCorrelator::result(std::vector<double>& lag, std::vector<double>& corr) const {
  lag.clear(), corr.clear();
  double scale = 1;
  for (int l = 0; l < nlev_; ++l, scale *= m_)
    for (int j = l == 0 ? 0 : p_ / m_; j < p_; ++j)
      if (count_[l][j] > 0) lag.push_back(j * scale), corr.push_back(sum_[l][j] / count_[l][j]);
}

// ---------------------------------------------------------------- stress relaxation

void dynamic_moduli(const std::vector<double>& t, const std::vector<double>& g, const std::vector<double>& omega, double g_inf,
                    std::vector<double>& storage, std::vector<double>& loss) {
  // integration by parts of iω∫G e^(−iωt) dt with G piecewise linear (slopes s_k) from t = 0 and g_inf after the end:
  // G* = G(0) + (g_inf − G(t_N)) e^(−iωt_N) + Σ_k s_k (e^(−iωt_k) − e^(−iωt_{k+1}))/(iω)
  storage.assign(omega.size(), NaN), loss.assign(omega.size(), NaN);
  if (t.size() < 2 || g.size() != t.size()) return;
  using C = std::complex<double>;
  const C I(0, 1);
  for (size_t w = 0; w < omega.size(); ++w) {
    const double om = omega[w];
    if (!(om > 0)) continue;
    C sum = g[0] + (g_inf - g.back()) * std::exp(-I * om * t.back());
    for (size_t k = 0; k + 1 < t.size(); ++k) {
      const double dt = t[k + 1] - t[k];
      if (dt <= 0) continue;
      const double s = (g[k + 1] - g[k]) / dt;
      sum += s * (std::exp(-I * om * t[k]) - std::exp(-I * om * t[k + 1])) / (I * om);
    }
    storage[w] = sum.real(), loss[w] = sum.imag();
  }
}

std::array<double, 3> fit_kww(const std::vector<double>& t, const std::vector<double>& y, double t_from, double lo, double hi) {
  size_t i0 = 0;
  while (i0 < t.size() && t[i0] < t_from) ++i0;
  if (i0 >= t.size() || !(y[i0] > 0)) return {NaN, NaN, NaN};
  std::vector<double> x, z;
  for (size_t k = i0; k < t.size(); ++k) {
    const double c = y[k] / y[i0];
    if (c > lo && c < hi && t[k] > 0) x.push_back(std::log(t[k])), z.push_back(std::log(-std::log(c)));
  }
  if (x.size() < 2) return {NaN, NaN, NaN};
  const Line l = line_fit(x, z);
  if (!(l.b > 0)) return {NaN, NaN, NaN};
  const double beta = l.b, tau = std::exp(-l.a / l.b);
  return {tau, beta, tau / beta * std::tgamma(1.0 / beta)};
}

StressRelaxation stress_relaxation(const std::vector<std::array<double, 6>>& p, const StressRelaxationOptions& o) {
  if (!(o.volume > 0)) throw std::invalid_argument("stress relaxation: the cell volume (Å³) is needed");
  if (!(o.temperature > 0)) throw std::invalid_argument("stress relaxation: the temperature (K) is needed");
  if (!(o.dt_fs > 0)) throw std::invalid_argument("stress relaxation: the time between samples (fs) is needed");
  StressRelaxation r;
  r.samples = p.size();
  if (p.size() < 4) { r.notes.push_back("fewer than four pressure-tensor samples"); return r; }
  std::array<double, 6> mean{};
  if (o.subtract_mean)
    for (const auto& q : p)
      for (int k = 0; k < 6; ++k) mean[k] += q[k] / double(p.size());
  // the traceless symmetric tensor as six channels whose products sum to Σ_αβ P°_αβ P°_αβ
  MultiTauCorrelator corr(6, o.p, o.m);
  const double r2 = std::sqrt(2.0);
  std::vector<double> ch(6);
  for (const auto& q0 : p) {
    std::array<double, 6> q;
    for (int k = 0; k < 6; ++k) q[k] = q0[k] - mean[k];
    const double tr = (q[0] + q[1] + q[2]) / 3.0;
    ch = {q[0] - tr, q[1] - tr, q[2] - tr, r2 * q[3], r2 * q[4], r2 * q[5]};
    corr.add(ch);
  }
  std::vector<double> lag, c;
  corr.result(lag, c);
  const double factor = o.volume * 1e-30 / (10.0 * kBoltzmann * o.temperature) * kAtmPa * kAtmPa * 1e-6;   // atm² → MPa
  const double max_lag = std::max(4.0, o.max_lag_fraction * double(p.size()));
  for (size_t k = 0; k < lag.size(); ++k)
    if (lag[k] <= max_lag) r.t_ps.push_back(lag[k] * o.dt_fs * 1e-3), r.g.push_back(factor * c[k]);
  r.g0 = r.g.empty() ? 0 : r.g[0];
  // the Green–Kubo integral, to integrate_ps or to the first lag where G reaches zero (beyond it the integral adds noise)
  for (size_t k = 0; k + 1 < r.t_ps.size(); ++k) {
    if (o.integrate_ps > 0 ? r.t_ps[k] >= o.integrate_ps : r.g[k] <= 0) break;
    r.viscosity += 0.5 * (r.g[k] + r.g[k + 1]) * (r.t_ps[k + 1] - r.t_ps[k]);
    r.integrated_to = r.t_ps[k + 1];
  }
  r.viscosity *= 1e-3;   // MPa·ps → mPa·s
  const auto kw = fit_kww(r.t_ps, r.g, o.fit_from_ps);
  r.kww_tau = kw[0], r.kww_beta = kw[1], r.kww_mean = kw[2];
  if (r.t_ps.size() > 2) {
    const double w0 = 1.0 / r.t_ps.back(), w1 = 1.0 / r.t_ps[1];
    const int n = std::max(2, o.omega_points);
    for (int k = 0; k < n; ++k) r.omega.push_back(w0 * std::pow(w1 / w0, double(k) / (n - 1)));
    dynamic_moduli(r.t_ps, r.g, r.omega, o.g_inf, r.storage, r.loss);
  }
  size_t i0 = 0;
  while (i0 < r.t_ps.size() && r.t_ps[i0] < o.fit_from_ps) ++i0;
  if (i0 < r.g.size() && r.g[i0] > 0 && r.g.back() > 0.05 * r.g[i0])
    r.notes.push_back("G(t) has not decayed to 5 % of its value at " + fmt("%g", o.fit_from_ps) +
                      " ps by the longest lag: the run is shorter than the stress relaxation, and the viscosity and fits are lower bounds");
  if (std::isnan(r.kww_tau)) r.notes.push_back("no KWW fit: too few points between 95 % and 5 % of G at the start of the fit");
  if (o.g_inf != 0) r.notes.push_back("G'(ω), G''(ω) with G = " + fmt("%g", o.g_inf) + " MPa after the last lag");
  return r;
}

std::vector<Property> stress_relaxation_properties(const StressRelaxation& r) {
  Property p{"stress_relaxation", "Stress relaxation modulus", "MPa", "", NaN, NaN, {}, {}, {}};
  p.method = "G(t) = V/(10 k_B T) Σ_αβ ⟨P°_αβ(0) P°_αβ(t)⟩ over the traceless symmetric pressure tensor (fluctuations), multi-tau "
             "correlator (Ramírez et al. 2010); G'(ω), G''(ω) by the exact transform of the piecewise-linear G(t); value: G(0)";
  p.notes = r.notes;
  if (r.t_ps.empty()) return {p};
  p.value = r.g0;
  p.extra["samples"] = double(r.samples);
  p.extra["Green–Kubo viscosity (mPa·s)"] = r.viscosity;
  p.extra["viscosity integrated to (ps)"] = r.integrated_to;
  if (std::isfinite(r.kww_tau)) {
    p.extra["KWW τ (ps)"] = r.kww_tau;
    p.extra["KWW β"] = r.kww_beta;
    p.extra["KWW mean time τΓ(1/β)/β (ps)"] = r.kww_mean;
  }
  p.series = {Series{"G(t)", "t (ps)", "G (MPa)", r.t_ps, r.g}, Series{"G'(ω)", "ω (rad/ps)", "G' (MPa)", r.omega, r.storage},
              Series{"G''(ω)", "ω (rad/ps)", "G'' (MPa)", r.omega, r.loss}};
  return {p};
}

// ---------------------------------------------------------------- bond memory

TwoStepFit fit_two_step(const std::vector<double>& t, const std::vector<double>& f) {
  TwoStepFit r;
  std::vector<double> x, y;
  for (size_t k = 0; k < t.size() && k < f.size(); ++k)
    if (t[k] > 0 && std::isfinite(f[k])) x.push_back(t[k]), y.push_back(f[k]);
  if (x.size() < 4) return r;
  auto model = [](double tt, double tau, double td, double b) { return 0.5 * (1 + std::exp(-tt / tau)) * std::exp(-std::pow(tt / td, b)); };
  auto sse = [&](const std::vector<double>& v) {
    const double b = v[2];
    if (b < 0.1 || b > 2.0) return 1e30;
    const double tau = std::exp(v[0]), td = std::exp(v[1]);
    double s = 0;
    for (size_t k = 0; k < x.size(); ++k) { const double d = model(x[k], tau, td, b) - y[k]; s += d * d; }
    return s;
  };
  // starting points: τ_d where f falls below 1/(2e) (else the last time), τ several fractions of it
  double td0 = x.back();
  for (size_t k = 0; k < x.size(); ++k)
    if (y[k] < 0.5 * kInvE) { td0 = x[k]; break; }
  std::vector<double> best;
  double bv = 1e300;
  for (double frac : {1e-3, 1e-2, 1e-1, 0.5})
    for (double b0 : {0.5, 0.9}) {
      const auto v = nelder_mead(sse, {std::log(frac * td0), std::log(td0), b0}, {1.0, 1.0, 0.2});
      const double e = sse(v);
      if (e < bv) bv = e, best = v;
    }
  if (best.empty() || bv >= 1e29) return r;
  r.ok = true;
  r.tau = std::exp(best[0]), r.tau_d = std::exp(best[1]), r.beta = best[2];
  r.rms = std::sqrt(bv / double(x.size()));
  return r;
}

BondMemory bond_memory(const std::vector<std::vector<int64_t>>& partners, double frame_dt) {
  BondMemory r;
  r.frames = partners.size();
  if (partners.empty()) { r.notes.push_back("no frames"); return r; }
  r.sites = partners[0].size();
  for (const auto& f : partners)
    if (f.size() != r.sites) throw std::invalid_argument("bond memory: every frame must list the same sites");
  if (r.frames < 2 || r.sites == 0) { r.notes.push_back("bond memory needs two frames or more and one site or more"); return r; }
  // lags: every lag to 32, then about 5 % apart (the cost stays near 100 lags × frames × sites)
  std::vector<size_t> lags;
  for (size_t l = 0; l < r.frames; l = l < 32 ? l + 1 : std::max(l + 1, size_t(double(l) * 1.05))) lags.push_back(l);
  for (size_t l : lags) {
    double same = 0, bonded = 0;
    for (size_t o = 0; o + l < r.frames; ++o)
      for (size_t i = 0; i < r.sites; ++i)
        if (partners[o][i] >= 0) {
          bonded += 1;
          if (partners[o + l][i] == partners[o][i]) same += 1;
        }
    if (bonded == 0) continue;
    r.t.push_back(double(l) * frame_dt), r.c.push_back(same / bonded);
  }
  double bonded_pairs = 0;
  for (size_t f = 0; f + 1 < r.frames; ++f)
    for (size_t i = 0; i < r.sites; ++i) {
      if (partners[f][i] >= 0) bonded_pairs += 1;
      if (partners[f][i] != partners[f + 1][i]) ++r.changes;
    }
  if (r.changes > 0) r.change_time = bonded_pairs * frame_dt / double(r.changes);
  for (size_t k = 1; k < r.c.size(); ++k)
    if (r.c[k] < kInvE) {
      const double a = r.c[k - 1], b = r.c[k];
      r.tau_e = r.t[k - 1] + (a - kInvE) / (a - b) * (r.t[k] - r.t[k - 1]);
      break;
    }
  const auto kw = fit_kww(r.t, r.c, 0);
  r.kww_tau = kw[0], r.kww_beta = kw[1], r.kww_mean = kw[2];
  if (!r.c.empty() && r.c.back() < 0.9) r.two_step = fit_two_step(r.t, r.c);
  if (r.c.empty()) r.notes.push_back("no site is bonded in any frame");
  else if (r.c.back() > kInvE)
    r.notes.push_back("the partners have not decorrelated to 1/e by the longest lag: the run is shorter than the bond lifetime");
  r.notes.push_back("a site that leaves its partner and comes back to it counts as still bonded to it");
  return r;
}

std::vector<Property> bond_memory_properties(const BondMemory& r, const std::string& time_unit) {
  Property p{"bond_memory", "Bond lifetime (partner autocorrelation)", time_unit, "", NaN, NaN, {}, {}, {}};
  p.method = "c(t): exchangeable sites bonded at a time origin that are bonded to the same partner a time t later, over all "
             "origins; value: c = 1/e; KWW fit on 0.05 < c < 0.95 and the two-step form ½[1 + e^(−t/τ)] exp[−(t/τ_d)^β] "
             "(Branham-Ferrari et al. 2026)";
  p.notes = r.notes;
  p.value = r.tau_e;
  p.extra["sites"] = double(r.sites);
  p.extra["frames"] = double(r.frames);
  p.extra["partner changes between frames"] = double(r.changes);
  if (std::isfinite(r.change_time)) p.extra["mean time between partner changes seen (" + time_unit + ")"] = r.change_time;
  if (std::isfinite(r.kww_tau)) {
    p.extra["KWW τ (" + time_unit + ")"] = r.kww_tau;
    p.extra["KWW β"] = r.kww_beta;
    p.extra["KWW mean time (" + time_unit + ")"] = r.kww_mean;
  }
  if (r.two_step.ok) {
    p.extra["two-step τ (" + time_unit + ")"] = r.two_step.tau;
    p.extra["two-step τ_d (" + time_unit + ")"] = r.two_step.tau_d;
    p.extra["two-step β"] = r.two_step.beta;
    p.extra["two-step rms residual"] = r.two_step.rms;
  }
  p.series.push_back(Series{"c(t)", "t (" + time_unit + ")", "same partner", r.t, r.c});
  return {p};
}

// ---------------------------------------------------------------- healing

std::vector<int8_t> sides_from_positions(const System& s, int axis, double plane) {
  if (axis < 0 || axis > 2) throw std::invalid_argument("sides: axis must be 0, 1 or 2");
  std::vector<int8_t> side(s.atoms.size(), 0);
  if (s.cell.valid()) {
    for (size_t i = 0; i < s.atoms.size(); ++i) {
      double f = s.cell.to_fractional(s.atoms[i].pos)[axis];
      f -= std::floor(f);
      side[i] = f < plane ? 0 : 1;
    }
    return side;
  }
  double lo = 1e300, hi = -1e300;
  for (const auto& a : s.atoms) lo = std::min(lo, a.pos[axis]), hi = std::max(hi, a.pos[axis]);
  for (size_t i = 0; i < s.atoms.size(); ++i) side[i] = s.atoms[i].pos[axis] < lo + plane * (hi - lo) ? 0 : 1;
  return side;
}

InterfaceBonds interface_bonds(const System& s, const std::vector<int8_t>& side, int axis, double plane, double bin) {
  if (axis < 0 || axis > 2) throw std::invalid_argument("interface bonds: axis must be 0, 1 or 2");
  if (side.size() != s.atoms.size()) throw std::invalid_argument("interface bonds: one side per atom");
  if (!(bin > 0)) throw std::invalid_argument("interface bonds: the bin width must be positive");
  InterfaceBonds r;
  r.bonds = s.bonds.size();
  const bool periodic = s.cell.valid();
  // positions along the axis as fractions of the cell (periodic) or of the span of the atoms (no cell)
  double length = 0, lo = 0;
  if (periodic) {
    length = thickness(s.cell, axis);
  } else {
    double hi = -1e300;
    lo = 1e300;
    for (const auto& a : s.atoms) lo = std::min(lo, a.pos[axis]), hi = std::max(hi, a.pos[axis]);
    length = std::max(hi - lo, 1e-9);
  }
  const int nb = std::max(1, int(std::lround(length / bin)));
  std::vector<double> planes(nb);
  for (int k = 0; k < nb; ++k) planes[k] = double(k) / nb, r.z.push_back(planes[k] * length);
  r.crossings.assign(nb, 0.0);
  // the number of planes q + n (n an integer) with a < q + n ≤ b
  auto crossed = [&](double a, double b, double q) { return periodic ? std::floor(b - q) - std::floor(a - q) : (a < q && q <= b ? 1.0 : 0.0); };
  for (const auto& b : s.bonds) {
    const size_t i = b.i, j = b.j;
    if (i >= s.atoms.size() || j >= s.atoms.size()) continue;
    if (side[i] >= 0 && side[j] >= 0 && side[i] != side[j]) ++r.bridges;
    double fi, fj;
    if (periodic) {
      const Vec3 d = s.cell.minimum_image(s.atoms[j].pos - s.atoms[i].pos);
      const double f0 = s.cell.to_fractional(s.atoms[i].pos)[axis];
      fi = f0 - std::floor(f0);
      fj = fi + (s.cell.to_fractional(s.atoms[i].pos + d)[axis] - f0);
    } else {
      fi = (s.atoms[i].pos[axis] - lo) / length, fj = (s.atoms[j].pos[axis] - lo) / length;
    }
    const double a = std::min(fi, fj), c = std::max(fi, fj);
    for (int k = 0; k < nb; ++k) r.crossings[k] += crossed(a, c, planes[k]);
    if (crossed(a, c, plane) > 0) ++r.crossing_plane;
  }
  return r;
}

std::vector<Property> interface_properties(const InterfaceBonds& r, int axis) {
  const std::string ax = axis == 0 ? "x" : axis == 1 ? "y" : "z";
  Property p{"interface_bonds", "Bonds across the interface", "bonds", "", NaN, NaN, {}, {}, {}};
  p.method = "bridges: bonds joining an atom that started on one side of the cut to one that started on the other; B(" + ax +
             "): bonds crossing each plane normal to " + ax + " (minimum image)";
  p.value = double(r.bridges);
  p.extra["bonds in the structure"] = double(r.bonds);
  p.extra["bonds crossing the cut plane"] = double(r.crossing_plane);
  p.series.push_back(Series{"B(" + ax + ")", ax + " (Å)", "bonds crossing", r.z, r.crossings});
  return {p};
}

namespace {

void check_curve(const StressStrain& c, const char* what) {
  if (c.strain.size() != c.stress.size() || c.strain.size() < 2)
    throw std::invalid_argument(std::string("healing efficiency: the ") + what + " curve needs two points or more, one stress per strain");
}

double interp(const StressStrain& c, double e) {
  if (e <= c.strain.front()) return c.stress.front();
  for (size_t k = 1; k < c.strain.size(); ++k)
    if (e <= c.strain[k]) {
      const double d = c.strain[k] - c.strain[k - 1];
      return d > 0 ? c.stress[k - 1] + (e - c.strain[k - 1]) / d * (c.stress[k] - c.stress[k - 1]) : c.stress[k];
    }
  return NaN;   // beyond the curve's end
}

double area(const StressStrain& c, double up_to) {
  double a = 0;
  for (size_t k = 0; k + 1 < c.strain.size(); ++k) {
    double e0 = c.strain[k], e1 = c.strain[k + 1], s0 = c.stress[k], s1 = c.stress[k + 1];
    if (up_to > 0) {
      if (e0 >= up_to) break;
      if (e1 > up_to) s1 = s0 + (up_to - e0) / (e1 - e0) * (s1 - s0), e1 = up_to;
    }
    a += 0.5 * (s0 + s1) * (e1 - e0);
  }
  return a;
}

double initial_slope(const StressStrain& c, double upto) {
  std::vector<double> x, y;
  for (size_t k = 0; k < c.strain.size(); ++k)
    if (c.strain[k] >= 0 && c.strain[k] <= upto) x.push_back(c.strain[k]), y.push_back(c.stress[k]);
  return line_fit(x, y).b;
}

double stress_near(const StressStrain& c, double e, double w) {
  if (w <= 0) return interp(c, e);
  double s = 0;
  int n = 0;
  for (size_t k = 0; k < c.strain.size(); ++k)
    if (std::fabs(c.strain[k] - e) <= w) s += c.stress[k], ++n;
  return n ? s / n : interp(c, e);
}

}  // namespace

HealingEfficiency healing_efficiency(const StressStrain& healed, const StressStrain& intact, double modulus_strain, double at_strain,
                                     double up_to, double at_window) {
  check_curve(healed, "healed");
  check_curve(intact, "intact");
  HealingEfficiency h;
  h.up_to = up_to, h.modulus_strain = modulus_strain, h.at_strain = at_strain;
  const double ai = area(intact, up_to);
  if (ai > 0) h.toughness = 100.0 * area(healed, up_to) / ai;
  const double pi = *std::max_element(intact.stress.begin(), intact.stress.end());
  if (pi > 0) h.peak = 100.0 * *std::max_element(healed.stress.begin(), healed.stress.end()) / pi;
  if (modulus_strain > 0) {
    const double ei = initial_slope(intact, modulus_strain), eh = initial_slope(healed, modulus_strain);
    if (ei > 0 && std::isfinite(eh)) h.modulus = 100.0 * eh / ei;
    else h.notes.push_back("no initial slope: fewer than two points with strain up to " + fmt("%g", modulus_strain));
  }
  if (at_strain > 0) {
    const double si = stress_near(intact, at_strain, at_window), sh = stress_near(healed, at_strain, at_window);
    if (std::isfinite(si) && si != 0 && std::isfinite(sh)) h.stress_at = 100.0 * sh / si;
    else h.notes.push_back("a curve does not reach strain " + fmt("%g", at_strain));
  }
  if (up_to > 0 && (healed.strain.back() < up_to || intact.strain.back() < up_to))
    h.notes.push_back("a curve ends before strain " + fmt("%g", up_to) + ": its area runs to its own end");
  if (up_to <= 0) h.notes.push_back("areas run to each curve's own last point (failure, or the end of the test)");
  return h;
}

std::vector<Property> healing_properties(const HealingEfficiency& h) {
  Property p{"healing", "Healing efficiency (toughness ratio)", "%", "", NaN, NaN, {}, {}, {}};
  p.method = "healed against intact stress–strain curves: toughness ratio ∫σ dε (Chen et al. 2023), peak-stress ratio, initial "
             "modulus ratio (straight-line fit to strain " + fmt("%g", h.modulus_strain) + ", Yuan et al. 2026) and the stress at a "
             "set strain (Wick et al. 2025); the definitions differ, so each is named";
  p.value = h.toughness;
  p.notes = h.notes;
  if (std::isfinite(h.peak)) p.extra["peak-stress ratio (%)"] = h.peak;
  if (std::isfinite(h.modulus)) p.extra["initial modulus ratio (%)"] = h.modulus;
  if (std::isfinite(h.stress_at)) p.extra["stress ratio at strain " + fmt("%g", h.at_strain) + " (%)"] = h.stress_at;
  if (h.up_to > 0) p.extra["areas to strain"] = h.up_to;
  return {p};
}

}  // namespace caps
