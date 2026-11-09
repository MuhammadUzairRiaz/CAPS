// CAPS structure–property models (see qspr.hpp).
#include "caps/qspr.hpp"
#include "caps/rng.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <functional>
#include <map>
#include <numeric>
#include <random>
#include <set>
#include <stdexcept>

#include "caps/charges.hpp"
#include "caps/elements.hpp"
#include "caps/molinfo.hpp"
#include "caps/normal_modes.hpp"
#include "caps/typing.hpp"

namespace caps {

namespace {
const double kNaN = std::numeric_limits<double>::quiet_NaN();

// valence electrons of a main-group element (s and p blocks); −1 otherwise
int valence_electrons(int z) {
  if (z <= 2) return z;
  if (z <= 10) return z - 2;
  if (z <= 18) return z - 10;
  if (z == 19 || z == 20) return z - 18;
  if (z >= 31 && z <= 36) return z - 28;
  if (z == 37 || z == 38) return z - 36;
  if (z >= 49 && z <= 54) return z - 46;
  if (z == 55 || z == 56) return z - 54;
  if (z >= 81 && z <= 86) return z - 78;
  return -1;
}

std::string fmt(double v) {
  char b[64];
  std::snprintf(b, sizeof b, "%.6g", v);
  return b;
}

// ---- dense linear algebra (small)
using Mat = std::vector<std::vector<double>>;
// inverse by Gauss–Jordan with partial pivoting; false when singular
bool invert(Mat a, Mat& inv) {
  const size_t n = a.size();
  inv.assign(n, std::vector<double>(n, 0.0));
  for (size_t i = 0; i < n; ++i) inv[i][i] = 1;
  double scale = 0;
  for (const auto& r : a)
    for (double v : r) scale = std::max(scale, std::fabs(v));
  for (size_t c = 0; c < n; ++c) {
    size_t p = c;
    for (size_t i = c + 1; i < n; ++i) if (std::fabs(a[i][c]) > std::fabs(a[p][c])) p = i;
    if (std::fabs(a[p][c]) <= 1e-12 * std::max(1.0, scale)) return false;
    std::swap(a[c], a[p]), std::swap(inv[c], inv[p]);
    const double d = a[c][c];
    for (size_t j = 0; j < n; ++j) a[c][j] /= d, inv[c][j] /= d;
    for (size_t i = 0; i < n; ++i)
      if (i != c && a[i][c] != 0) {
        const double m = a[i][c];
        for (size_t j = 0; j < n; ++j) a[i][j] -= m * a[c][j], inv[i][j] -= m * inv[c][j];
      }
  }
  return true;
}

struct Ols {
  std::vector<double> beta, fitted, resid;
  Mat inv;   // (XᵀX)⁻¹
  double sse = 0;
  bool ok = false;
};
// least squares with the design matrix's columns given (no intercept added)
Ols ols(const std::vector<double>& y, const std::vector<std::vector<double>>& cols) {
  Ols r;
  const size_t n = y.size(), p = cols.size();
  Mat A(p, std::vector<double>(p, 0.0));
  std::vector<double> b(p, 0.0);
  for (size_t i = 0; i < p; ++i) {
    for (size_t j = i; j < p; ++j) {
      double s = 0;
      for (size_t k = 0; k < n; ++k) s += cols[i][k] * cols[j][k];
      A[i][j] = A[j][i] = s;
    }
    for (size_t k = 0; k < n; ++k) b[i] += cols[i][k] * y[k];
  }
  if (!invert(A, r.inv)) return r;
  r.beta.assign(p, 0.0);
  for (size_t i = 0; i < p; ++i)
    for (size_t j = 0; j < p; ++j) r.beta[i] += r.inv[i][j] * b[j];
  r.fitted.assign(n, 0.0), r.resid.assign(n, 0.0);
  for (size_t k = 0; k < n; ++k) {
    for (size_t i = 0; i < p; ++i) r.fitted[k] += r.beta[i] * cols[i][k];
    r.resid[k] = y[k] - r.fitted[k];
    r.sse += r.resid[k] * r.resid[k];
  }
  r.ok = true;
  return r;
}

double mean(const std::vector<double>& v) { return v.empty() ? 0 : std::accumulate(v.begin(), v.end(), 0.0) / double(v.size()); }
double ss_about_mean(const std::vector<double>& v) {
  const double m = mean(v);
  double s = 0;
  for (double x : v) s += (x - m) * (x - m);
  return s;
}

// ---- terms of linear models: a function of the columns and its name
struct Term {
  int kind = 0;        // 0 constant, 1 linear x_a, 2 square x_a², 3 product x_a x_b, 4 hinge+ <x_a − t>₊, 5 hinge− <t − x_a>₊
  int a = 0, b = 0;
  double t = 0;
  double eval(const std::vector<double>& row) const {
    switch (kind) {
      case 0: return 1;
      case 1: return row[size_t(a)];
      case 2: return row[size_t(a)] * row[size_t(a)];
      case 3: return row[size_t(a)] * row[size_t(b)];
      case 4: return std::max(0.0, row[size_t(a)] - t);
      default: return std::max(0.0, t - row[size_t(a)]);
    }
  }
  std::string name(const std::vector<std::string>& nm) const {
    switch (kind) {
      case 0: return "1";
      case 1: return nm[size_t(a)];
      case 2: return nm[size_t(a)] + "^2";
      case 3: return nm[size_t(a)] + "*" + nm[size_t(b)];
      case 4: return "<" + nm[size_t(a)] + "-" + fmt(t) + ">+";
      default: return "<" + fmt(t) + "-" + nm[size_t(a)] + ">+";
    }
  }
  bool operator<(const Term& o) const { return std::tie(kind, a, b, t) < std::tie(o.kind, o.a, o.b, o.t); }
  bool operator==(const Term& o) const { return kind == o.kind && a == o.a && b == o.b && t == o.t; }
};

// rows (row-major) from columns
std::vector<std::vector<double>> to_rows(const std::vector<std::vector<double>>& x, size_t n) {
  std::vector<std::vector<double>> r(n, std::vector<double>(x.size()));
  for (size_t j = 0; j < x.size(); ++j)
    for (size_t i = 0; i < n; ++i) r[i][j] = x[j][i];
  return r;
}
std::vector<std::vector<double>> design(const std::vector<Term>& terms, const std::vector<std::vector<double>>& rows) {
  std::vector<std::vector<double>> c(terms.size(), std::vector<double>(rows.size()));
  for (size_t t = 0; t < terms.size(); ++t)
    for (size_t i = 0; i < rows.size(); ++i) c[t][i] = terms[t].eval(rows[i]);
  return c;
}

// A fitted model that can predict rows; refit(rows, y) fits the same form on other data (cross-validation).
struct Fitted {
  std::function<double(const std::vector<double>&)> predict;
};
using Refit = std::function<Fitted(const std::vector<std::vector<double>>&, const std::vector<double>&)>;

Fitted linear_refit(const std::vector<Term>& terms, const std::vector<std::vector<double>>& rows, const std::vector<double>& y) {
  const Ols f = ols(y, design(terms, rows));
  std::vector<double> beta = f.ok ? f.beta : std::vector<double>(terms.size(), kNaN);
  return {[terms, beta](const std::vector<double>& r) {
    double s = 0;
    for (size_t t = 0; t < terms.size(); ++t) s += beta[t] * terms[t].eval(r);
    return s;
  }};
}

// ---- PLS1 by NIPALS on autoscaled x and y: the regression coefficients in the original units (intercept first)
std::vector<double> pls_coefficients(const std::vector<std::vector<double>>& rows, const std::vector<double>& y, int A) {
  const size_t n = rows.size(), p = rows.empty() ? 0 : rows[0].size();
  std::vector<double> mx(p, 0.0), sx(p, 0.0);
  for (size_t j = 0; j < p; ++j) {
    for (size_t i = 0; i < n; ++i) mx[j] += rows[i][j];
    mx[j] /= double(n);
    for (size_t i = 0; i < n; ++i) sx[j] += (rows[i][j] - mx[j]) * (rows[i][j] - mx[j]);
    sx[j] = std::sqrt(sx[j] / double(n > 1 ? n - 1 : 1));
    if (sx[j] == 0) sx[j] = 1;
  }
  const double my = mean(y);
  double sy = std::sqrt(ss_about_mean(y) / double(n > 1 ? n - 1 : 1));
  if (sy == 0) sy = 1;
  Mat X(n, std::vector<double>(p));
  std::vector<double> yy(n);
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < p; ++j) X[i][j] = (rows[i][j] - mx[j]) / sx[j];
    yy[i] = (y[i] - my) / sy;
  }
  Mat W, P;
  std::vector<double> C;
  for (int a = 0; a < A; ++a) {
    std::vector<double> w(p, 0.0);
    for (size_t j = 0; j < p; ++j)
      for (size_t i = 0; i < n; ++i) w[j] += X[i][j] * yy[i];
    double nw = 0;
    for (double v : w) nw += v * v;
    nw = std::sqrt(nw);
    if (nw < 1e-14) break;
    for (double& v : w) v /= nw;
    std::vector<double> t(n, 0.0);
    for (size_t i = 0; i < n; ++i)
      for (size_t j = 0; j < p; ++j) t[i] += X[i][j] * w[j];
    double tt = 0, yt = 0;
    for (size_t i = 0; i < n; ++i) tt += t[i] * t[i], yt += yy[i] * t[i];
    std::vector<double> pl(p, 0.0);
    for (size_t j = 0; j < p; ++j) {
      for (size_t i = 0; i < n; ++i) pl[j] += X[i][j] * t[i];
      pl[j] /= tt;
    }
    const double c = yt / tt;
    for (size_t i = 0; i < n; ++i) {
      for (size_t j = 0; j < p; ++j) X[i][j] -= t[i] * pl[j];
      yy[i] -= c * t[i];
    }
    W.push_back(w), P.push_back(pl), C.push_back(c);
  }
  const size_t a = W.size();
  // B = W (PᵀW)⁻¹ c
  Mat PW(a, std::vector<double>(a, 0.0)), inv;
  for (size_t i = 0; i < a; ++i)
    for (size_t j = 0; j < a; ++j)
      for (size_t k = 0; k < p; ++k) PW[i][j] += P[i][k] * W[j][k];
  std::vector<double> B(p, 0.0);
  if (a > 0 && invert(PW, inv)) {
    std::vector<double> g(a, 0.0);
    for (size_t i = 0; i < a; ++i)
      for (size_t j = 0; j < a; ++j) g[i] += inv[i][j] * C[j];
    for (size_t k = 0; k < p; ++k)
      for (size_t i = 0; i < a; ++i) B[k] += W[i][k] * g[i];
  }
  std::vector<double> out(p + 1, 0.0);
  out[0] = my;
  for (size_t k = 0; k < p; ++k) {
    out[k + 1] = B[k] * sy / sx[k];
    out[0] -= out[k + 1] * mx[k];
  }
  return out;
}
Fitted pls_refit(const std::vector<std::vector<double>>& rows, const std::vector<double>& y, int A) {
  const std::vector<double> b = pls_coefficients(rows, y, A);
  return {[b](const std::vector<double>& r) {
    double s = b[0];
    for (size_t k = 0; k < r.size(); ++k) s += b[k + 1] * r[k];
    return s;
  }};
}

// ---- neural network: one tanh hidden layer, inputs and output autoscaled; loss ½Σ(ŷ − y)² + ½λ Σ w² over the weights
struct Net {
  size_t p = 0, H = 0;
  std::vector<double> mx, sx;
  double my = 0, sy = 1;
  std::vector<double> w;   // H×p input weights, H biases, H output weights, 1 output bias
  double raw(const std::vector<double>& z) const {
    double out = w[H * p + 2 * H];
    for (size_t h = 0; h < H; ++h) {
      double a = w[H * p + h];
      for (size_t j = 0; j < p; ++j) a += w[h * p + j] * z[j];
      out += w[H * p + H + h] * std::tanh(a);
    }
    return out;
  }
  double predict(const std::vector<double>& r) const {
    std::vector<double> z(p);
    for (size_t j = 0; j < p; ++j) z[j] = (r[j] - mx[j]) / sx[j];
    return my + sy * raw(z);
  }
};

double net_loss(const Net& net, const std::vector<std::vector<double>>& Z, const std::vector<double>& t, double lambda, std::vector<double>& g) {
  const size_t H = net.H, p = net.p;
  g.assign(net.w.size(), 0.0);
  double L = 0;
  std::vector<double> act(H);
  for (size_t i = 0; i < Z.size(); ++i) {
    double out = net.w[H * p + 2 * H];
    for (size_t h = 0; h < H; ++h) {
      double a = net.w[H * p + h];
      for (size_t j = 0; j < p; ++j) a += net.w[h * p + j] * Z[i][j];
      act[h] = std::tanh(a);
      out += net.w[H * p + H + h] * act[h];
    }
    const double e = out - t[i];
    L += 0.5 * e * e;
    g[H * p + 2 * H] += e;
    for (size_t h = 0; h < H; ++h) {
      g[H * p + H + h] += e * act[h];
      const double da = e * net.w[H * p + H + h] * (1 - act[h] * act[h]);
      g[H * p + h] += da;
      for (size_t j = 0; j < p; ++j) g[h * p + j] += da * Z[i][j];
    }
  }
  for (size_t k = 0; k < H * p; ++k) L += 0.5 * lambda * net.w[k] * net.w[k], g[k] += lambda * net.w[k];
  for (size_t h = 0; h < H; ++h) {
    const double v = net.w[H * p + H + h];
    L += 0.5 * lambda * v * v, g[H * p + H + h] += lambda * v;
  }
  return L;
}

// L-BFGS with backtracking on f(x) (gradient g)
template <class F>
double lbfgs(std::vector<double>& x, F f, int iters) {
  const size_t n = x.size();
  std::vector<double> g, gn, xn;
  double fx = f(x, g);
  std::vector<std::vector<double>> S, Y;
  std::vector<double> rho;
  for (int it = 0; it < iters; ++it) {
    std::vector<double> q = g, al(S.size());
    for (size_t k = S.size(); k-- > 0;) {
      double d = 0;
      for (size_t v = 0; v < n; ++v) d += S[k][v] * q[v];
      al[k] = rho[k] * d;
      for (size_t v = 0; v < n; ++v) q[v] -= al[k] * Y[k][v];
    }
    double gam = 1e-2;
    if (!S.empty()) {
      double sy = 0, yy = 0;
      for (size_t v = 0; v < n; ++v) sy += S.back()[v] * Y.back()[v], yy += Y.back()[v] * Y.back()[v];
      if (yy > 0) gam = sy / yy;
    }
    for (double& v : q) v *= gam;
    for (size_t k = 0; k < S.size(); ++k) {
      double d = 0;
      for (size_t v = 0; v < n; ++v) d += Y[k][v] * q[v];
      const double be = rho[k] * d;
      for (size_t v = 0; v < n; ++v) q[v] += S[k][v] * (al[k] - be);
    }
    double slope = 0, gmax = 0;
    for (size_t v = 0; v < n; ++v) slope -= q[v] * g[v], gmax = std::max(gmax, std::fabs(g[v]));
    if (gmax < 1e-10) break;
    if (slope >= 0) {
      S.clear(), Y.clear(), rho.clear();
      q = g;
      slope = 0;
      for (size_t v = 0; v < n; ++v) q[v] *= 1e-2, slope -= q[v] * g[v];
    }
    double step = 1, fn = fx;
    bool ok = false;
    for (int ls = 0; ls < 40; ++ls) {
      xn = x;
      for (size_t v = 0; v < n; ++v) xn[v] -= step * q[v];
      fn = f(xn, gn);
      if (fn <= fx + 1e-4 * step * slope) { ok = true; break; }
      step *= 0.5;
    }
    if (!ok) break;
    std::vector<double> s(n), yv(n);
    double sy = 0;
    for (size_t v = 0; v < n; ++v) s[v] = xn[v] - x[v], yv[v] = gn[v] - g[v], sy += s[v] * yv[v];
    if (sy > 1e-16) {
      S.push_back(s), Y.push_back(yv), rho.push_back(1 / sy);
      if (S.size() > 10) S.erase(S.begin()), Y.erase(Y.begin()), rho.erase(rho.begin());
    }
    const bool small = std::fabs(fx - fn) < 1e-12 * std::max(1.0, std::fabs(fx));
    x.swap(xn), g.swap(gn), fx = fn;
    if (small) break;
  }
  return fx;
}

Net train_net(const std::vector<std::vector<double>>& rows, const std::vector<double>& y, const ModelOptions& o, uint64_t seed) {
  Net net;
  const size_t n = rows.size();
  net.p = rows[0].size(), net.H = size_t(std::max(1, o.nn_hidden));
  net.mx.assign(net.p, 0), net.sx.assign(net.p, 0);
  for (size_t j = 0; j < net.p; ++j) {
    for (size_t i = 0; i < n; ++i) net.mx[j] += rows[i][j];
    net.mx[j] /= double(n);
    for (size_t i = 0; i < n; ++i) net.sx[j] += (rows[i][j] - net.mx[j]) * (rows[i][j] - net.mx[j]);
    net.sx[j] = std::sqrt(net.sx[j] / double(std::max<size_t>(1, n - 1)));
    if (net.sx[j] == 0) net.sx[j] = 1;
  }
  net.my = mean(y);
  net.sy = std::sqrt(ss_about_mean(y) / double(std::max<size_t>(1, n - 1)));
  if (net.sy == 0) net.sy = 1;
  std::vector<std::vector<double>> Z(n, std::vector<double>(net.p));
  std::vector<double> t(n);
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < net.p; ++j) Z[i][j] = (rows[i][j] - net.mx[j]) / net.sx[j];
    t[i] = (y[i] - net.my) / net.sy;
  }
  std::mt19937_64 rng(seed);
  caps::Normal<double> N(0.0, 1.0);
  double best = std::numeric_limits<double>::infinity();
  std::vector<double> bw;
  const size_t nw = net.H * net.p + 2 * net.H + 1;
  for (int r = 0; r < std::max(1, o.nn_restarts); ++r) {
    std::vector<double> w(nw);
    for (size_t k = 0; k < nw; ++k) w[k] = N(rng) / std::sqrt(double(net.p));
    Net trial = net;
    auto f = [&](const std::vector<double>& x, std::vector<double>& g) {
      trial.w = x;
      return net_loss(trial, Z, t, o.nn_decay, g);
    };
    const double L = lbfgs(w, f, 3000);
    if (L < best) best = L, bw = w;
  }
  net.w = bw;
  return net;
}

// ---- GFA
struct GfaModel {
  std::vector<Term> basis;   // without the constant
  double lof = std::numeric_limits<double>::infinity();
  std::vector<double> beta;
};
void gfa_score(GfaModel& m, const std::vector<std::vector<double>>& rows, const std::vector<double>& y, double d) {
  std::vector<Term> all{Term{}};
  all.insert(all.end(), m.basis.begin(), m.basis.end());
  const Ols f = ols(y, design(all, rows));
  const double M = double(y.size()), c = double(m.basis.size()), p = double(m.basis.size());
  const double den = 1 - (c + d * p) / M;
  m.lof = (!f.ok || den <= 0) ? std::numeric_limits<double>::infinity() : (f.sse / M) / (den * den);
  m.beta = f.ok ? f.beta : std::vector<double>();
}

// a better model: a lower LOF, or as good (to rounding: an exact fit) with fewer terms
bool gfa_better(const GfaModel& a, const GfaModel& b) {
  if (!std::isfinite(a.lof) || !std::isfinite(b.lof)) return a.lof < b.lof;   // a singular model never ties with a real one
  const double tol = 1e-9 * std::max(std::fabs(a.lof), std::fabs(b.lof)) + 1e-24;
  if (std::fabs(a.lof - b.lof) <= tol) return a.basis.size() < b.basis.size();
  return a.lof < b.lof;
}

GfaModel run_gfa(const std::vector<std::vector<double>>& rows, const std::vector<double>& y, const ModelOptions& o, uint64_t seed,
                 std::vector<GfaModel>* top) {
  std::mt19937_64 rng(seed);
  const size_t n = rows.size(), p = rows[0].size();
  auto rnd = [&](size_t k) { return size_t(caps::UniformInt<size_t>(0, k - 1)(rng)); };
  auto random_term = [&] {
    Term t;
    const int kinds[4] = {1, 2, 4, 5};
    t.kind = kinds[rnd(4)];
    t.a = int(rnd(p));
    if (t.kind == 4 || t.kind == 5) t.t = rows[rnd(n)][size_t(t.a)];
    return t;
  };
  auto key = [](const GfaModel& m) {
    std::vector<Term> b = m.basis;
    std::sort(b.begin(), b.end());
    return b;
  };
  std::vector<GfaModel> pop;
  std::set<std::vector<Term>> seen;
  const size_t maxt = size_t(std::max(1, o.gfa_max_terms));
  for (int tries = 0; pop.size() < size_t(o.gfa_population) && tries < o.gfa_population * 20; ++tries) {
    GfaModel m;
    const size_t k = 1 + rnd(std::min<size_t>(3, maxt));
    while (m.basis.size() < k) {
      const Term t = random_term();
      if (std::find(m.basis.begin(), m.basis.end(), t) == m.basis.end()) m.basis.push_back(t);
    }
    if (!seen.insert(key(m)).second) continue;
    gfa_score(m, rows, y, o.gfa_smoothing);
    pop.push_back(m);
  }
  auto pick = [&] {   // binary tournament
    const size_t a = rnd(pop.size()), b = rnd(pop.size());
    return gfa_better(pop[a], pop[b]) ? a : b;
  };
  for (int g = 0; g < o.gfa_generations; ++g) {
    const GfaModel& A = pop[pick()];
    const GfaModel& B = pop[pick()];
    GfaModel c;
    const size_t ca = rnd(A.basis.size() + 1), cb = rnd(B.basis.size() + 1);
    for (size_t i = 0; i < ca; ++i) c.basis.push_back(A.basis[i]);
    for (size_t i = cb; i < B.basis.size(); ++i)
      if (std::find(c.basis.begin(), c.basis.end(), B.basis[i]) == c.basis.end()) c.basis.push_back(B.basis[i]);
    if (rnd(2) == 0) {   // mutation: a new term
      const Term t = random_term();
      if (std::find(c.basis.begin(), c.basis.end(), t) == c.basis.end()) c.basis.push_back(t);
    }
    if (c.basis.empty() || c.basis.size() > maxt) continue;
    if (!seen.insert(key(c)).second) continue;
    gfa_score(c, rows, y, o.gfa_smoothing);
    const auto worst = std::max_element(pop.begin(), pop.end(), [](const GfaModel& a, const GfaModel& b) { return gfa_better(a, b); });
    if (gfa_better(c, *worst)) *worst = c;
  }
  std::sort(pop.begin(), pop.end(), gfa_better);
  if (top) top->assign(pop.begin(), pop.begin() + std::min<size_t>(5, pop.size()));
  return pop.front();
}

// ---- the whole fit for one y (no validation): the model's description, fitted values and the refit of its form
struct Core {
  ModelResult r;
  Refit refit;
};

Core fit_core(const std::vector<double>& y, const std::vector<std::vector<double>>& rows, const std::vector<std::string>& names, const ModelOptions& o,
              uint64_t seed) {
  Core c;
  ModelResult& r = c.r;
  const size_t n = y.size(), p = names.size();
  r.method = o.method;
  auto finish_linear = [&](const std::vector<Term>& terms) {
    const Ols f = ols(y, design(terms, rows));
    if (!f.ok) throw std::invalid_argument("the columns are linearly dependent (or there are too few rows)");
    const size_t k = terms.size();
    const double s2 = n > k ? f.sse / double(n - k) : kNaN;
    for (size_t t = 0; t < k; ++t) r.terms.push_back({terms[t].name(names), f.beta[t], std::sqrt(s2 * f.inv[t][t])});
    r.fitted = f.fitted, r.residuals = f.resid;
    // leverage, standardised residuals, Cook's distance
    const auto D = design(terms, rows);
    r.leverage.assign(n, 0.0), r.std_residual.assign(n, kNaN), r.cooks.assign(n, kNaN);
    for (size_t i = 0; i < n; ++i) {
      double h = 0;
      for (size_t a = 0; a < k; ++a)
        for (size_t b = 0; b < k; ++b) h += D[a][i] * f.inv[a][b] * D[b][i];
      r.leverage[i] = h;
      if (std::isfinite(s2) && s2 > 0 && h < 1) {
        r.std_residual[i] = f.resid[i] / std::sqrt(s2 * (1 - h));
        r.cooks[i] = f.resid[i] * f.resid[i] / (double(k) * s2) * h / ((1 - h) * (1 - h));
      }
    }
    std::string eq = "y = ";
    for (size_t t = 0; t < k; ++t) eq += (t ? (f.beta[t] < 0 ? " − " : " + ") : "") + fmt(t ? std::fabs(f.beta[t]) : f.beta[t]) + (terms[t].kind ? " " + terms[t].name(names) : "");
    r.equation = eq;
    c.refit = [terms](const std::vector<std::vector<double>>& rr, const std::vector<double>& yy) { return linear_refit(terms, rr, yy); };
  };
  if (o.method == "mlr") {
    std::vector<Term> terms{Term{}};
    if (!o.stepwise) {
      for (size_t j = 0; j < p; ++j) terms.push_back({1, int(j), 0, 0});
    } else {
      std::vector<Term> cand;
      for (size_t j = 0; j < p; ++j) cand.push_back({1, int(j), 0, 0});
      if (o.squares)
        for (size_t j = 0; j < p; ++j) cand.push_back({2, int(j), 0, 0});
      if (o.products)
        for (size_t a = 0; a < p; ++a)
          for (size_t b = a + 1; b < p; ++b) cand.push_back({3, int(a), int(b), 0});
      auto aic = [&](const std::vector<Term>& t) {
        if (t.size() >= n) return std::numeric_limits<double>::infinity();
        const Ols f = ols(y, design(t, rows));
        if (!f.ok) return std::numeric_limits<double>::infinity();
        return double(n) * std::log(std::max(f.sse, 1e-300) / double(n)) + 2.0 * double(t.size() + 1);
      };
      double cur = aic(terms);
      for (int it = 0; it < 200; ++it) {
        double best = cur;
        std::vector<Term> bt;
        for (const auto& t : cand) {
          if (std::find(terms.begin(), terms.end(), t) != terms.end()) continue;
          auto tr = terms;
          tr.push_back(t);
          const double a = aic(tr);
          if (a < best - 1e-9) best = a, bt = tr;
        }
        for (size_t q = 1; q < terms.size(); ++q) {
          auto tr = terms;
          tr.erase(tr.begin() + long(q));
          const double a = aic(tr);
          if (a < best - 1e-9) best = a, bt = tr;
        }
        if (bt.empty()) break;
        terms = bt, cur = best;
      }
      r.notes.push_back("stepwise by AIC from " + std::to_string(cand.size()) + " candidate terms: " + std::to_string(terms.size() - 1) + " kept");
    }
    finish_linear(terms);
  } else if (o.method == "gfa") {
    std::vector<GfaModel> top;
    const GfaModel m = run_gfa(rows, y, o, seed, &top);
    std::vector<Term> terms{Term{}};
    terms.insert(terms.end(), m.basis.begin(), m.basis.end());
    finish_linear(terms);
    r.lof = m.lof;
    for (const auto& t : top) {
      std::string eq;
      for (size_t i = 0; i < t.basis.size(); ++i) eq += (i ? " + " : "") + t.basis[i].name(names);
      r.gfa_models.push_back({eq, t.lof});
    }
  } else if (o.method == "pls") {
    const int amax = int(std::min<size_t>({p, n > 2 ? n - 2 : 1, 10}));
    int A = o.components;
    if (A <= 0) {   // by cross-validation
      const int k = o.folds > 0 ? o.folds : int(n);
      double bestq = -std::numeric_limits<double>::infinity();
      for (int a = 1; a <= amax; ++a) {
        double press = 0;
        for (int f = 0; f < std::min<int>(k, int(n)); ++f) {
          std::vector<std::vector<double>> tr;
          std::vector<double> ty;
          for (size_t i = 0; i < n; ++i)
            if (int(i % size_t(std::min<int>(k, int(n)))) != f) tr.push_back(rows[i]), ty.push_back(y[i]);
          const Fitted m = pls_refit(tr, ty, a);
          for (size_t i = 0; i < n; ++i)
            if (int(i % size_t(std::min<int>(k, int(n)))) == f) press += std::pow(y[i] - m.predict(rows[i]), 2);
        }
        const double q2 = 1 - press / ss_about_mean(y);
        r.q2_by_components.push_back(q2);
        if (q2 > bestq + 1e-9) bestq = q2, A = a;
      }
    }
    A = std::max(1, std::min(A, int(p)));
    r.components = A;
    const std::vector<double> b = pls_coefficients(rows, y, A);
    r.terms.push_back({"1", b[0]});
    for (size_t j = 0; j < p; ++j) r.terms.push_back({names[j], b[j + 1]});
    r.fitted.resize(n), r.residuals.resize(n);
    for (size_t i = 0; i < n; ++i) {
      double s = b[0];
      for (size_t j = 0; j < p; ++j) s += b[j + 1] * rows[i][j];
      r.fitted[i] = s, r.residuals[i] = y[i] - s;
    }
    std::string eq = "y = " + fmt(b[0]);
    for (size_t j = 0; j < p; ++j) eq += (b[j + 1] < 0 ? " − " : " + ") + fmt(std::fabs(b[j + 1])) + " " + names[j];
    r.equation = eq + " (" + std::to_string(A) + " components)";
    c.refit = [A](const std::vector<std::vector<double>>& rr, const std::vector<double>& yy) { return pls_refit(rr, yy, A); };
  } else if (o.method == "nn") {
    const Net net = train_net(rows, y, o, seed);
    r.fitted.resize(n), r.residuals.resize(n);
    for (size_t i = 0; i < n; ++i) r.fitted[i] = net.predict(rows[i]), r.residuals[i] = y[i] - r.fitted[i];
    r.equation = "one hidden layer of " + std::to_string(net.H) + " tanh units, " + std::to_string(net.w.size()) + " weights";
    if (net.w.size() >= n) r.notes.push_back("as many weights as rows or more: the network can fit noise; judge it by its cross-validated q²");
    ModelOptions oo = o;
    c.refit = [oo, seed](const std::vector<std::vector<double>>& rr, const std::vector<double>& yy) {
      const Net m = train_net(rr, yy, oo, seed);
      return Fitted{[m](const std::vector<double>& q) { return m.predict(q); }};
    };
  } else {
    throw std::invalid_argument("model '" + o.method + "': mlr, pls, gfa or nn");
  }
  // statistics of the fit
  const double sst = ss_about_mean(y);
  double sse = 0;
  for (double e : r.residuals) sse += e * e;
  r.r2 = sst > 0 ? 1 - sse / sst : kNaN;
  r.rmse = std::sqrt(sse / double(n));
  const double k = o.method == "pls" ? double(r.components) : double(r.terms.size() > 0 ? r.terms.size() - 1 : 0);
  if (o.method != "nn" && double(n) - k - 1 > 0) {
    r.r2_adj = 1 - (1 - r.r2) * (double(n) - 1) / (double(n) - k - 1);
    if (k > 0 && sse > 0) r.f = ((sst - sse) / k) / (sse / (double(n) - k - 1));
  } else {
    r.r2_adj = kNaN;
  }
  return c;
}

double cross_validate(const Core& c, const std::vector<double>& y, const std::vector<std::vector<double>>& rows, int folds, std::vector<double>* pred) {
  const size_t n = y.size();
  const size_t k = size_t(std::min<int>(folds, int(n)));
  std::vector<double> pr(n, kNaN);
  for (size_t f = 0; f < k; ++f) {
    std::vector<std::vector<double>> tr;
    std::vector<double> ty;
    for (size_t i = 0; i < n; ++i)
      if (i % k != f) tr.push_back(rows[i]), ty.push_back(y[i]);
    const Fitted m = c.refit(tr, ty);
    for (size_t i = 0; i < n; ++i)
      if (i % k == f) pr[i] = m.predict(rows[i]);
  }
  double press = 0;
  for (size_t i = 0; i < n; ++i) press += (y[i] - pr[i]) * (y[i] - pr[i]);
  if (pred) *pred = pr;
  return press;
}
}  // namespace

ModelResult fit_model(const std::vector<double>& y, const std::vector<std::vector<double>>& x, const std::vector<std::string>& names, const ModelOptions& o) {
  const size_t n = y.size();
  if (x.size() != names.size()) throw std::invalid_argument("model: a name per column");
  if (x.empty()) throw std::invalid_argument("model: no columns");
  for (const auto& c : x)
    if (c.size() != n) throw std::invalid_argument("model: every column as long as y");
  if (n < 3) throw std::invalid_argument("model: three rows or more");
  for (double v : y) if (!std::isfinite(v)) throw std::invalid_argument("model: y has a missing value");
  for (const auto& c : x)
    for (double v : c) if (!std::isfinite(v)) throw std::invalid_argument("model: a column has a missing value");
  const auto rows = to_rows(x, n);
  Core c = fit_core(y, rows, names, o, o.seed);
  ModelResult r = c.r;
  if (o.folds > 1) {
    const double press = cross_validate(c, y, rows, o.folds, &r.predicted_cv);
    r.q2 = 1 - press / ss_about_mean(y);
    r.rmse_cv = std::sqrt(press / double(n));
    r.notes.push_back((o.folds >= int(n) ? std::string("leave-one-out") : std::to_string(o.folds) + "-fold") +
                      " cross-validation: the model's form refitted in each fold" + ((o.method == "mlr" && o.stepwise) || o.method == "gfa" ? " (its terms as selected on all rows)" : ""));
  }
  if (o.randomise > 0) {
    std::mt19937_64 rng(o.seed + 7);
    for (int k = 0; k < o.randomise; ++k) {
      std::vector<double> ys = y;
      caps::shuffle(ys.begin(), ys.end(), rng);
      try {
        Core cs = fit_core(ys, rows, names, o, o.seed + uint64_t(k) + 11);
        r.random_r2.push_back(cs.r.r2);
        if (o.folds > 1) r.random_q2.push_back(1 - cross_validate(cs, ys, rows, o.folds, nullptr) / ss_about_mean(ys));
      } catch (const std::invalid_argument&) {
      }
    }
    if (!r.random_r2.empty()) {
      const double mx = *std::max_element(r.random_r2.begin(), r.random_r2.end());
      char b[200];
      std::snprintf(b, sizeof b, "y-randomisation: %zu shuffled fits, R² mean %.3f, largest %.3f (the model's %.3f)", r.random_r2.size(), mean(r.random_r2), mx, r.r2);
      r.notes.push_back(b);
    }
  }
  return r;
}

// ---------------------------------------------------------------- clustering
RowClusters cluster_rows(const std::vector<std::vector<double>>& rows0, const std::string& method, int k, bool scale, uint64_t seed, int restarts) {
  RowClusters r;
  const size_t n = rows0.size();
  if (n == 0) throw std::invalid_argument("cluster: no rows");
  const size_t p = rows0[0].size();
  if (k < 1 || size_t(k) > n) throw std::invalid_argument("cluster: k between 1 and the row count");
  std::vector<std::vector<double>> X = rows0;
  if (scale)
    for (size_t j = 0; j < p; ++j) {
      double m = 0, s = 0;
      for (size_t i = 0; i < n; ++i) m += X[i][j];
      m /= double(n);
      for (size_t i = 0; i < n; ++i) s += (X[i][j] - m) * (X[i][j] - m);
      s = std::sqrt(s / double(std::max<size_t>(1, n - 1)));
      for (size_t i = 0; i < n; ++i) X[i][j] = s > 0 ? (X[i][j] - m) / s : 0;
    }
  auto d2 = [&](const std::vector<double>& a, const std::vector<double>& b) {
    double s = 0;
    for (size_t j = 0; j < p; ++j) s += (a[j] - b[j]) * (a[j] - b[j]);
    return s;
  };
  if (method == "kmeans") {
    std::mt19937_64 rng(seed);
    double best = std::numeric_limits<double>::infinity();
    for (int rs = 0; rs < std::max(1, restarts); ++rs) {
      std::vector<size_t> idx(n);
      std::iota(idx.begin(), idx.end(), 0);
      caps::shuffle(idx.begin(), idx.end(), rng);
      std::vector<std::vector<double>> C;
      for (int c = 0; c < k; ++c) C.push_back(X[idx[size_t(c)]]);
      std::vector<int> lab(n, -1);
      for (int it = 0; it < 1000; ++it) {
        bool changed = false;
        for (size_t i = 0; i < n; ++i) {
          int bc = 0;
          for (int c = 1; c < k; ++c) if (d2(X[i], C[size_t(c)]) < d2(X[i], C[size_t(bc)])) bc = c;
          if (lab[i] != bc) lab[i] = bc, changed = true;
        }
        if (!changed) break;
        for (int c = 0; c < k; ++c) {
          std::vector<double> m(p, 0.0);
          int cnt = 0;
          for (size_t i = 0; i < n; ++i)
            if (lab[i] == c) { ++cnt; for (size_t j = 0; j < p; ++j) m[j] += X[i][j]; }
          if (cnt) {
            for (double& v : m) v /= cnt;
            C[size_t(c)] = m;
          }
        }
      }
      double w = 0;
      for (size_t i = 0; i < n; ++i) w += d2(X[i], C[size_t(lab[i])]);
      if (w < best - 1e-12) best = w, r.label = lab;
    }
  } else {
    const int kind = method == "single" ? 0 : method == "complete" ? 1 : method == "average" ? 2 : method == "ward" ? 3 : -1;
    if (kind < 0) throw std::invalid_argument("cluster method '" + method + "': single, complete, average, ward or kmeans");
    // D: distances (Ward: the increase of the within sum of squares on merging, ½|a − b|² for two rows)
    std::vector<std::vector<double>> D(n, std::vector<double>(n, 0.0));
    for (size_t i = 0; i < n; ++i)
      for (size_t j = i + 1; j < n; ++j) D[i][j] = D[j][i] = kind == 3 ? 0.5 * d2(X[i], X[j]) : std::sqrt(d2(X[i], X[j]));
    std::vector<int> id(n), size(n, 1);
    std::iota(id.begin(), id.end(), 0);
    std::vector<bool> alive(n, true);
    std::vector<int> parent(2 * n, -1);
    int next = int(n);
    for (size_t step = 0; step + 1 < n; ++step) {
      size_t a = 0, b = 0;
      double best = std::numeric_limits<double>::infinity();
      for (size_t i = 0; i < n; ++i)
        if (alive[i])
          for (size_t j = i + 1; j < n; ++j)
            if (alive[j] && D[i][j] < best) best = D[i][j], a = i, b = j;
      r.merges.push_back({double(std::min(id[a], id[b])), double(std::max(id[a], id[b])), best});
      parent[size_t(id[a])] = parent[size_t(id[b])] = next;
      // Lance–Williams: the merged cluster takes slot a
      const double na = size[a], nb = size[b];
      for (size_t m = 0; m < n; ++m) {
        if (!alive[m] || m == a || m == b) continue;
        const double nm = size[m];
        double v;
        switch (kind) {
          case 0: v = std::min(D[a][m], D[b][m]); break;
          case 1: v = std::max(D[a][m], D[b][m]); break;
          case 2: v = (na * D[a][m] + nb * D[b][m]) / (na + nb); break;
          default: v = ((na + nm) * D[a][m] + (nb + nm) * D[b][m] - nm * D[a][b]) / (na + nb + nm);
        }
        D[a][m] = D[m][a] = v;
      }
      alive[b] = false;
      size[a] += size[b];
      id[a] = next++;
    }
    // the partition into k: the first n − k merges
    std::vector<int> uf(n);
    std::iota(uf.begin(), uf.end(), 0);
    std::function<int(int)> find = [&](int x) { return uf[size_t(x)] == x ? x : uf[size_t(x)] = find(uf[size_t(x)]); };
    std::vector<int> rep(2 * n, -1);   // cluster id → a row in it
    for (size_t i = 0; i < n; ++i) rep[i] = int(i);
    for (size_t s = 0; s + size_t(k) < n; ++s) {
      const int A = rep[size_t(r.merges[s][0])], B = rep[size_t(r.merges[s][1])];
      uf[size_t(find(A))] = find(B);
      rep[n + s] = B;
    }
    std::map<int, int> num;
    r.label.assign(n, 0);
    for (size_t i = 0; i < n; ++i) {
      const int root = find(int(i));
      if (!num.count(root)) num[root] = int(num.size());
      r.label[i] = num[root];
    }
  }
  // centroids (original units) and the within sum of squares (the clustered columns)
  r.centroids.assign(size_t(k), std::vector<double>(p, 0.0));
  std::vector<std::vector<double>> cs(size_t(k), std::vector<double>(p, 0.0));
  std::vector<int> cnt(size_t(k), 0);
  for (size_t i = 0; i < n; ++i) {
    ++cnt[size_t(r.label[i])];
    for (size_t j = 0; j < p; ++j) r.centroids[size_t(r.label[i])][j] += rows0[i][j], cs[size_t(r.label[i])][j] += X[i][j];
  }
  for (int c = 0; c < k; ++c)
    for (size_t j = 0; j < p; ++j)
      if (cnt[size_t(c)]) r.centroids[size_t(c)][j] /= cnt[size_t(c)], cs[size_t(c)][j] /= cnt[size_t(c)];
  for (size_t i = 0; i < n; ++i) r.within_ss += d2(X[i], cs[size_t(r.label[i])]);
  return r;
}

// ---------------------------------------------------------------- discriminant analysis
namespace {
struct Lda {
  std::vector<std::vector<double>> dirs;
  std::vector<double> eig;
  std::vector<std::vector<double>> means;   // per class, projected
};
Lda lda_fit(const std::vector<std::vector<double>>& X, const std::vector<int>& cls, int C) {
  const size_t n = X.size(), p = X[0].size();
  std::vector<std::vector<double>> mu(size_t(C), std::vector<double>(p, 0.0));
  std::vector<double> all(p, 0.0);
  std::vector<int> cnt(size_t(C), 0);
  for (size_t i = 0; i < n; ++i) {
    ++cnt[size_t(cls[i])];
    for (size_t j = 0; j < p; ++j) mu[size_t(cls[i])][j] += X[i][j], all[j] += X[i][j];
  }
  for (int c = 0; c < C; ++c)
    for (size_t j = 0; j < p; ++j) if (cnt[size_t(c)]) mu[size_t(c)][j] /= cnt[size_t(c)];
  for (double& v : all) v /= double(n);
  Mat Sw(p, std::vector<double>(p, 0.0)), Sb(p, std::vector<double>(p, 0.0));
  for (size_t i = 0; i < n; ++i)
    for (size_t a = 0; a < p; ++a)
      for (size_t b = 0; b < p; ++b) Sw[a][b] += (X[i][a] - mu[size_t(cls[i])][a]) * (X[i][b] - mu[size_t(cls[i])][b]);
  for (int c = 0; c < C; ++c)
    for (size_t a = 0; a < p; ++a)
      for (size_t b = 0; b < p; ++b) Sb[a][b] += cnt[size_t(c)] * (mu[size_t(c)][a] - all[a]) * (mu[size_t(c)][b] - all[b]);
  // Cholesky of S_W (a tiny ridge when singular)
  double tr = 0;
  for (size_t a = 0; a < p; ++a) tr += Sw[a][a];
  for (size_t a = 0; a < p; ++a) Sw[a][a] += 1e-10 * std::max(tr, 1e-300);
  Mat L(p, std::vector<double>(p, 0.0));
  for (size_t i = 0; i < p; ++i)
    for (size_t j = 0; j <= i; ++j) {
      double s = Sw[i][j];
      for (size_t k = 0; k < j; ++k) s -= L[i][k] * L[j][k];
      if (i == j) L[i][i] = std::sqrt(std::max(s, 1e-300));
      else L[i][j] = s / L[j][j];
    }
  Mat Li;
  invert(L, Li);
  // M = L⁻¹ S_B L⁻ᵀ
  std::vector<double> M(p * p, 0.0);
  for (size_t a = 0; a < p; ++a)
    for (size_t b = 0; b < p; ++b) {
      double s = 0;
      for (size_t c = 0; c < p; ++c)
        for (size_t d = 0; d < p; ++d) s += Li[a][c] * Sb[c][d] * Li[b][d];
      M[a * p + b] = s;
    }
  std::vector<double> w;
  symmetric_eigen(M, p, w);   // ascending; eigenvectors in the columns of M
  Lda r;
  const size_t keep = std::min<size_t>(p, size_t(std::max(1, C - 1)));
  for (size_t q = 0; q < keep; ++q) {
    const size_t col = p - 1 - q;
    std::vector<double> u(p), v(p, 0.0);
    for (size_t a = 0; a < p; ++a) u[a] = M[a * p + col];
    for (size_t a = 0; a < p; ++a)
      for (size_t b = 0; b < p; ++b) v[a] += Li[b][a] * u[b];   // L⁻ᵀ u
    r.dirs.push_back(v), r.eig.push_back(w[col]);
  }
  for (int c = 0; c < C; ++c) {
    std::vector<double> m;
    for (const auto& d : r.dirs) {
      double s = 0;
      for (size_t j = 0; j < p; ++j) s += d[j] * mu[size_t(c)][j];
      m.push_back(s);
    }
    r.means.push_back(m);
  }
  return r;
}
int lda_class(const Lda& m, const std::vector<double>& x) {
  std::vector<double> z;
  for (const auto& d : m.dirs) {
    double s = 0;
    for (size_t j = 0; j < x.size(); ++j) s += d[j] * x[j];
    z.push_back(s);
  }
  int best = 0;
  double bd = std::numeric_limits<double>::infinity();
  for (size_t c = 0; c < m.means.size(); ++c) {
    double dd = 0;
    for (size_t q = 0; q < z.size(); ++q) dd += (z[q] - m.means[c][q]) * (z[q] - m.means[c][q]);
    if (dd < bd) bd = dd, best = int(c);
  }
  return best;
}
}  // namespace

DiscriminantResult discriminant(const std::vector<std::vector<double>>& rows, const std::vector<std::string>& cls) {
  DiscriminantResult r;
  const size_t n = rows.size();
  if (n != cls.size() || n < 3) throw std::invalid_argument("discriminant: a class per row, three rows or more");
  std::map<std::string, int> idx;
  for (const auto& c : cls)
    if (!idx.count(c)) idx[c] = int(idx.size()), r.classes.push_back(c);
  if (r.classes.size() < 2) throw std::invalid_argument("discriminant: two classes or more");
  std::vector<int> y(n);
  for (size_t i = 0; i < n; ++i) y[i] = idx[cls[i]];
  const int C = int(r.classes.size());
  const Lda m = lda_fit(rows, y, C);
  r.directions = m.dirs, r.eigenvalues = m.eig;
  int ok = 0, okl = 0;
  for (size_t i = 0; i < n; ++i) {
    r.predicted.push_back(lda_class(m, rows[i]));
    ok += r.predicted.back() == y[i];
    std::vector<std::vector<double>> tr;
    std::vector<int> ty;
    for (size_t j = 0; j < n; ++j)
      if (j != i) tr.push_back(rows[j]), ty.push_back(y[j]);
    r.predicted_loo.push_back(lda_class(lda_fit(tr, ty, C), rows[i]));
    okl += r.predicted_loo.back() == y[i];
  }
  r.accuracy = double(ok) / double(n), r.accuracy_loo = double(okl) / double(n);
  return r;
}

// ---------------------------------------------------------------- descriptors
std::vector<Descriptor> molecular_descriptors(const System& s, std::string* charges_used) {
  std::vector<Descriptor> out;
  auto add = [&](const std::string& n, double v, const std::string& u, const std::string& w) { out.push_back({n, u, w, v}); };
  const size_t N = s.atoms.size();
  if (N == 0) throw std::invalid_argument("descriptors: no atoms");
  const Perception pc = perceive(s);
  // heavy-atom graph
  std::vector<int> heavy_index(N, -1);
  std::vector<size_t> heavy;
  for (size_t i = 0; i < N; ++i)
    if (s.atoms[i].element > 1) heavy_index[i] = int(heavy.size()), heavy.push_back(i);
  const size_t A = heavy.size();
  std::vector<std::vector<int>> g(A);
  std::vector<int> hcount(A, 0);
  for (size_t a = 0; a < A; ++a) {
    const size_t i = heavy[a];
    for (uint32_t j : pc.nb[i]) {
      if (heavy_index[j] >= 0) g[a].push_back(heavy_index[j]);
    }
    hcount[a] = pc.hcount[i];
  }
  std::vector<double> delta(A), dv(A);
  bool dv_ok = true;
  for (size_t a = 0; a < A; ++a) {
    delta[a] = double(g[a].size());
    const int z = s.atoms[heavy[a]].element, zv = valence_electrons(z);
    if (zv < 0) { dv_ok = false, dv[a] = kNaN; continue; }
    dv[a] = double(zv - hcount[a]) / double(z - zv - 1);
  }
  // distances (BFS)
  std::vector<std::vector<int>> d(A, std::vector<int>(A, -1));
  bool connected = true;
  for (size_t a = 0; a < A; ++a) {
    std::vector<int> q{int(a)};
    d[a][a] = 0;
    for (size_t h = 0; h < q.size(); ++h)
      for (int b : g[size_t(q[h])])
        if (d[a][size_t(b)] < 0) d[a][size_t(b)] = d[a][size_t(q[h])] + 1, q.push_back(b);
    for (size_t b = 0; b < A; ++b) connected = connected && d[a][b] >= 0;
  }
  double W = 0;
  for (size_t a = 0; a < A; ++a)
    for (size_t b = a + 1; b < A; ++b) W += d[a][b];
  size_t q_bonds = 0;
  double M1 = 0, M2 = 0;
  for (size_t a = 0; a < A; ++a) {
    M1 += delta[a] * delta[a];
    for (int b : g[a])
      if (size_t(b) > a) ++q_bonds, M2 += delta[a] * delta[size_t(b)];
  }
  // paths of length 1..3 (each once) and the connectivity indices
  double P[4] = {double(A), 0, 0, 0}, chi[4] = {0, 0, 0, 0}, chiv[4] = {0, 0, 0, 0};
  for (size_t a = 0; a < A; ++a)
    if (delta[a] > 0) chi[0] += 1 / std::sqrt(delta[a]), chiv[0] += dv[a] > 0 ? 1 / std::sqrt(dv[a]) : kNaN;
  std::function<void(std::vector<int>&)> walk = [&](std::vector<int>& path) {
    const size_t m = path.size() - 1;
    if (m >= 1 && path.front() < path.back()) {
      P[m] += 1;
      double c = 1, cv = 1;
      for (int v : path) c /= std::sqrt(delta[size_t(v)]), cv /= std::sqrt(dv[size_t(v)]);
      chi[m] += c, chiv[m] += cv;
    }
    if (m == 3) return;
    for (int nb : g[size_t(path.back())])
      if (std::find(path.begin(), path.end(), nb) == path.end()) {
        path.push_back(nb);
        walk(path);
        path.pop_back();
      }
  };
  for (size_t a = 0; a < A; ++a) {
    std::vector<int> path{int(a)};
    walk(path);
  }
  const double Ad = double(A);
  const double k1 = P[1] > 0 ? Ad * (Ad - 1) * (Ad - 1) / (P[1] * P[1]) : kNaN;
  const double k2 = P[2] > 0 ? (Ad - 1) * (Ad - 2) * (Ad - 2) / (P[2] * P[2]) : kNaN;
  const double k3 = P[3] > 0 ? (A % 2 ? (Ad - 1) * (Ad - 3) * (Ad - 3) : (Ad - 3) * (Ad - 2) * (Ad - 2)) / (P[3] * P[3]) : kNaN;
  // Balaban J
  int comps = 0;
  {
    std::vector<int> seen(A, 0);
    for (size_t a = 0; a < A; ++a) {
      if (seen[a]) continue;
      ++comps;
      std::vector<int> st{int(a)};
      seen[a] = 1;
      while (!st.empty()) {
        const int v = st.back();
        st.pop_back();
        for (int b : g[size_t(v)]) if (!seen[size_t(b)]) seen[size_t(b)] = 1, st.push_back(b);
      }
    }
  }
  const double mu = double(q_bonds) - Ad + comps;
  double J = kNaN;
  if (connected && q_bonds > 0) {
    std::vector<double> sd(A, 0.0);
    for (size_t a = 0; a < A; ++a)
      for (size_t b = 0; b < A; ++b) sd[a] += d[a][b];
    double sum = 0;
    for (size_t a = 0; a < A; ++a)
      for (int b : g[a]) if (size_t(b) > a) sum += 1 / std::sqrt(sd[a] * sd[size_t(b)]);
    J = double(q_bonds) / (mu + 1) * sum;
  }
  add("W", connected ? W : kNaN, "", "Wiener index: half the sum of the topological distances between heavy atoms");
  add("M1", M1, "", "first Zagreb index Σ δ²");
  add("M2", M2, "", "second Zagreb index Σ over bonds δ_i δ_j");
  for (int m = 0; m <= 3; ++m) add(std::to_string(m) + "chi", chi[m], "", "connectivity index of order " + std::to_string(m) + " (paths)");
  for (int m = 0; m <= 3; ++m) add(std::to_string(m) + "chiv", dv_ok ? chiv[m] : kNaN, "", "valence connectivity index of order " + std::to_string(m) + " (paths)");
  add("kappa1", k1, "", "Kier shape index of order 1");
  add("kappa2", k2, "", "Kier shape index of order 2");
  add("kappa3", k3, "", "Kier shape index of order 3");
  add("J", J, "", "Balaban index");
  for (int m = 1; m <= 3; ++m) add("P" + std::to_string(m), P[m], "", "paths of length " + std::to_string(m));
  add("heavy_atoms", Ad, "", "atoms other than hydrogen");
  add("rings", double(pc.rings.size()), "", "smallest set of smallest rings");
  int arom = 0, hbd = 0, hba = 0;
  for (size_t i = 0; i < N; ++i) {
    arom += pc.aromatic[i] && s.atoms[i].element > 1;
    const int z = s.atoms[i].element;
    if (z == 7 || z == 8) hba += 1, hbd += pc.hcount[i];
  }
  add("aromatic_atoms", double(arom), "", "aromatic atoms");
  add("hbond_donors", double(hbd), "", "hydrogens on N and O (Lipinski)");
  add("hbond_acceptors", double(hba), "", "N and O atoms (Lipinski)");
  const MoleculeInfo mi = molecule_info(s, 0);
  add("rotatable_bonds", double(mi.rotatable), "", "acyclic single bonds between non-terminal heavy atoms");
  add("mass", mi.mass, "g/mol", "molecular mass");
  add("Rg", mi.rg, "Å", "mass-weighted radius of gyration");
  for (int k = 0; k < 3; ++k) add("I" + std::to_string(k + 1), mi.inertia[k], "amu·Å²", "principal moment of inertia (ascending)");
  // ---- 3D: van der Waals volume (grid), area, shadows
  std::vector<double> rad(N);
  Vec3 lo{1e300, 1e300, 1e300}, hi{-1e300, -1e300, -1e300};
  for (size_t i = 0; i < N; ++i) {
    rad[i] = element(s.atoms[i].element).vdw;
    for (int k = 0; k < 3; ++k) lo[size_t(k)] = std::min(lo[size_t(k)], s.atoms[i].pos[size_t(k)] - rad[i]), hi[size_t(k)] = std::max(hi[size_t(k)], s.atoms[i].pos[size_t(k)] + rad[i]);
  }
  {
    const double h = 0.1;
    int nx = int(std::ceil((hi[0] - lo[0]) / h)), ny = int(std::ceil((hi[1] - lo[1]) / h)), nz = int(std::ceil((hi[2] - lo[2]) / h));
    double vol = kNaN;
    if (double(nx) * ny * nz < 6e7) {
      size_t inside = 0;
      for (int ix = 0; ix < nx; ++ix) {
        const double x = lo[0] + (ix + 0.5) * h;
        std::vector<size_t> near;
        for (size_t i = 0; i < N; ++i) if (std::fabs(s.atoms[i].pos[0] - x) < rad[i]) near.push_back(i);
        for (int iy = 0; iy < ny; ++iy) {
          const double y = lo[1] + (iy + 0.5) * h;
          std::vector<size_t> near2;
          for (size_t i : near) if (std::fabs(s.atoms[i].pos[1] - y) < rad[i]) near2.push_back(i);
          if (near2.empty()) continue;
          for (int iz = 0; iz < nz; ++iz) {
            const double z = lo[2] + (iz + 0.5) * h;
            for (size_t i : near2) {
              const double dx = s.atoms[i].pos[0] - x, dy = s.atoms[i].pos[1] - y, dz = s.atoms[i].pos[2] - z;
              if (dx * dx + dy * dy + dz * dz < rad[i] * rad[i]) { ++inside; break; }
            }
          }
        }
      }
      vol = double(inside) * h * h * h;
    }
    add("vdw_volume", vol, "Å³", "van der Waals volume (Bondi radii, 0.1 Å grid)");
  }
  System sn = s;
  sn.cell = Cell{};
  const SasaResult vdw = sasa(sn, 0.0, 960);
  add("vdw_area", vdw.total, "Å²", "van der Waals surface area");
  {
    // principal axes (mass-weighted), longest first; shadows on XY, YZ, XZ
    Vec3 c{0, 0, 0};
    double m = 0;
    for (const auto& a : s.atoms) { const double w = s.mass_of(a); c = c + a.pos * w, m += w; }
    c = c * (1 / m);
    std::vector<double> I(9, 0.0), ev;
    for (const auto& a : s.atoms) {
      const double w = s.mass_of(a);
      const Vec3 r = a.pos - c;
      const double r2 = dot(r, r);
      for (int p = 0; p < 3; ++p)
        for (int q = 0; q < 3; ++q) I[size_t(3 * p + q)] += w * ((p == q ? r2 : 0) - r[size_t(p)] * r[size_t(q)]);
    }
    symmetric_eigen(I, 3, ev);   // ascending: column 0 the smallest moment = the longest axis
    std::vector<Vec3> P2(N);
    for (size_t i = 0; i < N; ++i) {
      const Vec3 r = s.atoms[i].pos - c;
      for (int k = 0; k < 3; ++k) P2[i][size_t(k)] = r[0] * I[size_t(0 * 3 + k)] + r[1] * I[size_t(1 * 3 + k)] + r[2] * I[size_t(2 * 3 + k)];
    }
    auto shadow = [&](int u, int v, double& area, double& rect) {
      double ulo = 1e300, uhi = -1e300, vlo = 1e300, vhi = -1e300;
      for (size_t i = 0; i < N; ++i) {
        ulo = std::min(ulo, P2[i][size_t(u)] - rad[i]), uhi = std::max(uhi, P2[i][size_t(u)] + rad[i]);
        vlo = std::min(vlo, P2[i][size_t(v)] - rad[i]), vhi = std::max(vhi, P2[i][size_t(v)] + rad[i]);
      }
      const double h = 0.05;
      const int nu = int(std::ceil((uhi - ulo) / h)), nv = int(std::ceil((vhi - vlo) / h));
      size_t cnt = 0;
      for (int a = 0; a < nu; ++a) {
        const double x = ulo + (a + 0.5) * h;
        for (int b = 0; b < nv; ++b) {
          const double y = vlo + (b + 0.5) * h;
          for (size_t i = 0; i < N; ++i) {
            const double dx = P2[i][size_t(u)] - x, dy = P2[i][size_t(v)] - y;
            if (dx * dx + dy * dy < rad[i] * rad[i]) { ++cnt; break; }
          }
        }
      }
      area = double(cnt) * h * h;
      rect = (uhi - ulo) * (vhi - vlo);
    };
    const char* nm[3] = {"XY", "YZ", "XZ"};
    const int uv[3][2] = {{0, 1}, {1, 2}, {0, 2}};
    for (int k = 0; k < 3; ++k) {
      double ar = 0, re = 1;
      shadow(uv[k][0], uv[k][1], ar, re);
      add(std::string("shadow_") + nm[k], ar, "Å²", std::string("area of the molecule's shadow on the principal ") + nm[k] + " plane (X the longest axis)");
      add(std::string("shadow_") + nm[k] + "_ratio", ar / re, "", std::string("the ") + nm[k] + " shadow over its enclosing rectangle");
    }
  }
  // ---- charged partial surface areas
  std::vector<double> q(N, 0.0);
  std::string used;
  if (s.has_charges) {
    for (size_t i = 0; i < N; ++i) q[i] = s.atoms[i].charge;
    used = "the structure's own";
  } else {
    try {
      q = compute_charges(s, "gasteiger").q;
      used = "Gasteiger–Marsili";
    } catch (const std::exception&) {
      q.clear();
      used = "none";
    }
  }
  if (charges_used) *charges_used = used;
  const SasaResult sa = sasa(sn, 1.4, 500);
  double ppsa1 = 0, pnsa1 = 0, ppsa3 = 0, pnsa3 = 0, qp = 0, qn = 0;
  if (!q.empty())
    for (size_t i = 0; i < N; ++i) {
      if (q[i] > 0) ppsa1 += sa.area[i], ppsa3 += q[i] * sa.area[i], qp += q[i];
      if (q[i] < 0) pnsa1 += sa.area[i], pnsa3 += q[i] * sa.area[i], qn += q[i];
    }
  const double T = sa.total;
  const double ok = q.empty() ? kNaN : 1.0;
  const double PP[3] = {ppsa1 * ok, ppsa1 * qp * ok, ppsa3 * ok}, PN[3] = {pnsa1 * ok, pnsa1 * qn * ok, pnsa3 * ok};
  add("SASA", T, "Å²", "solvent-accessible area (1.4 Å probe)");
  for (int k = 0; k < 3; ++k) {
    const std::string K = std::to_string(k + 1);
    add("PPSA" + K, PP[k], k ? "e·Å²" : "Å²", "partial positive surface area " + K);
    add("PNSA" + K, PN[k], k ? "e·Å²" : "Å²", "partial negative surface area " + K);
    add("DPSA" + K, PP[k] - PN[k], k ? "e·Å²" : "Å²", "difference PPSA − PNSA " + K);
    add("FPSA" + K, PP[k] / T, k ? "e" : "", "PPSA" + K + " / SASA");
    add("FNSA" + K, PN[k] / T, k ? "e" : "", "PNSA" + K + " / SASA");
    add("WPSA" + K, PP[k] * T / 1000, k ? "e·Å⁴" : "Å⁴", "PPSA" + K + " · SASA / 1000");
    add("WNSA" + K, PN[k] * T / 1000, k ? "e·Å⁴" : "Å⁴", "PNSA" + K + " · SASA / 1000");
  }
  return out;
}

Json descriptors_json(const std::vector<Descriptor>& d) {
  Json a = Json::array();
  for (const auto& x : d) {
    Json o = Json::object();
    o["name"] = x.name, o["unit"] = x.unit, o["what"] = x.what, o["value"] = std::isfinite(x.value) ? Json(x.value) : Json(nullptr);
    a.push_back(o);
  }
  return a;
}

Json model_json(const ModelResult& r) {
  auto val = [](double v) { return std::isfinite(v) ? Json(v) : Json(nullptr); };
  auto arr = [&](const std::vector<double>& v) { Json a = Json::array(); for (double x : v) a.push_back(val(x)); return a; };
  Json j = Json::object(), ts = Json::array(), gm = Json::array(), nt = Json::array();
  for (const auto& t : r.terms) {
    Json o = Json::object();
    o["name"] = t.name, o["coef"] = t.coef, o["se"] = val(t.se);
    ts.push_back(o);
  }
  for (const auto& [e, l] : r.gfa_models) {
    Json o = Json::object();
    o["terms"] = e, o["lof"] = val(l);
    gm.push_back(o);
  }
  for (const auto& n : r.notes) nt.push_back(n);
  j["method"] = r.method, j["equation"] = r.equation, j["terms"] = ts, j["fitted"] = arr(r.fitted), j["residuals"] = arr(r.residuals);
  j["r2"] = val(r.r2), j["r2_adj"] = val(r.r2_adj), j["rmse"] = val(r.rmse), j["f"] = val(r.f), j["components"] = double(r.components);
  j["q2"] = val(r.q2), j["rmse_cv"] = val(r.rmse_cv), j["predicted_cv"] = arr(r.predicted_cv), j["q2_by_components"] = arr(r.q2_by_components);
  j["lof"] = val(r.lof), j["gfa_models"] = gm, j["random_r2"] = arr(r.random_r2), j["random_q2"] = arr(r.random_q2);
  j["leverage"] = arr(r.leverage), j["std_residual"] = arr(r.std_residual), j["cooks"] = arr(r.cooks), j["notes"] = nt;
  return j;
}

}  // namespace caps
