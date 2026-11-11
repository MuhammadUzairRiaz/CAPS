// CAPS QEq charge equilibration (see caps/qeq.hpp).
#include "caps/qeq.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <stdexcept>

#include "caps/elements.hpp"
#include "cell_list.hpp"

namespace caps {

namespace {

constexpr double kCoul = 14.399645;            // e²/(4πε0) in eV Å
constexpr double kBohr = 0.529177210903;       // Å
constexpr double kHartree = 27.211386245988;   // eV

double dotv(const std::vector<double>& a, const std::vector<double>& b) {
  double s = 0;
  for (size_t i = 0; i < a.size(); ++i) s += a[i] * b[i];
  return s;
}

// The regularised lower incomplete gamma function P(a, x) for integer a ≥ 1.
double lower_p(int a, double x) {
  if (x <= 0) return 0;
  if (x < a + 12) {   // e^(−x) Σ_{j≥a} x^j / j!
    double term = std::exp(-x);
    for (int j = 1; j <= a; ++j) term *= x / j;
    double sum = 0;
    for (int j = a; j < a + 400; ++j) {
      sum += term;
      term *= x / (j + 1);
      if (term < 1e-18 * sum && j > x) break;
    }
    return sum;
  }
  double term = std::exp(-x), sum = 0;   // 1 − e^(−x) Σ_{j<a} x^j / j!
  for (int j = 0; j < a; ++j) {
    sum += term;
    term *= x / (j + 1);
  }
  return 1 - sum;
}

// The potential Φ(s) of a normalised ns Slater density ρ ∝ r^(2n−2) e^(−2ζr) (atomic units), and G(u) = ∫₀ᵘ s Φ(s) ds.
// With x = 2ζs: Φ(s) = P(2n+1, x)/s + (ζ/n) Q(2n, x); G from ∫ s^m e^(−2ζs) ds = m! P(m+1, x)/(2ζ)^(m+1).
double slater_phi(int n, double z, double s) {
  if (s < 1e-12) return z / n;
  const double x = 2 * z * s;
  return lower_p(2 * n + 1, x) / s + (z / n) * (1 - lower_p(2 * n, x));
}
double slater_g(int n, double z, double u) {
  const double x = 2 * z * u;
  double a = 0, b = 0;
  for (int k = 0; k <= 2 * n; ++k) a += lower_p(k + 1, x);
  for (int k = 0; k <= 2 * n - 1; ++k) b += (k + 1) * lower_p(k + 2, x);
  return u - a / (2 * z) + b / (4 * n * z);
}

// Gauss–Legendre nodes and weights on [−1, 1].
void gauss_legendre(int m, std::vector<double>& x, std::vector<double>& w) {
  x.assign(size_t(m), 0), w.assign(size_t(m), 0);
  for (int i = 0; i < m; ++i) {
    double t = std::cos(M_PI * (i + 0.75) / (m + 0.5)), dp = 1;
    for (int it = 0; it < 100; ++it) {
      double p0 = 1, p1 = t;
      for (int k = 2; k <= m; ++k) {
        const double p2 = ((2 * k - 1) * t * p1 - (k - 1) * p0) / k;
        p0 = p1, p1 = p2;
      }
      dp = m * (t * p1 - p0) / (t * t - 1);
      const double dt = p1 / dp;
      t -= dt;
      if (std::fabs(dt) < 1e-15) break;
    }
    x[size_t(i)] = t, w[size_t(i)] = 2 / ((1 - t * t) * dp * dp);
  }
}

// A Slater kernel per element pair: K(r) in eV on a grid of h Å, the exact kCoul/r once the densities no longer overlap.
struct SlaterTable {
  double h = 0.01;
  std::vector<double> v;   // K at r = k h
  double far = 0;          // from here on K = kCoul / r
  double at(double r) const {
    if (r >= far) return kCoul / r;
    const double u = r / h;
    const int k = std::clamp(int(u), 1, int(v.size()) - 3);
    const double t = u - k;   // Catmull–Rom through k−1 … k+2
    const double a = v[size_t(k - 1)], b = v[size_t(k)], c = v[size_t(k + 1)], d = v[size_t(k + 2)];
    return b + 0.5 * t * (c - a + t * (2 * a - 5 * b + 4 * c - d + t * (3 * (b - c) + d - a)));
  }
};

SlaterTable slater_table(int n1, double z1, int n2, double z2, double rmax) {
  SlaterTable T;
  T.far = 1e300;
  // beyond both densities' extent (4πr²ρ below e^(−55) of its peak) the integral is 1/r
  const double reach = kBohr * ((2 * n1 + 60) / (2 * z1) + (2 * n2 + 60) / (2 * z2));
  rmax = std::min(rmax, reach);
  for (int k = 0;; ++k) {
    const double r = k * T.h;
    const double J = kHartree * slater_coulomb_integral(n1, z1, n2, z2, r / kBohr);
    T.v.push_back(J);
    if (k > 4 && std::fabs(J - kCoul / r) < 1e-11 * kCoul / r) {   // the overlap is gone: three more points, then 1/r
      for (int e = 1; e <= 3; ++e) T.v.push_back(kCoul / ((k + e) * T.h));
      T.far = r;
      break;
    }
    if (r > rmax + 3 * T.h) {
      T.far = r;
      break;
    }
  }
  return T;
}

// The linear operator H (eV/e²) of the equalisation; the reciprocal part is applied as structure factors, never stored.
struct Operator {
  std::vector<double> diag;
  std::vector<double> pre;   // the full diagonal (with the reciprocal self term), for the Jacobi preconditioner
  std::vector<std::vector<std::pair<uint32_t, double>>> off;   // each pair once (j > i)
  bool ewald = false;
  int H[3] = {0, 0, 0};
  std::vector<std::complex<double>> ph[3];   // per atom e^(2πi h f_a), h = −H…H
  struct K { int h, k, l; double w; };
  std::vector<K> ks;                         // half space; w already doubled
  void apply(const std::vector<double>& x, std::vector<double>& y) const {
    const size_t n = diag.size();
    for (size_t i = 0; i < n; ++i) y[i] = diag[i] * x[i];
    for (size_t i = 0; i < n; ++i)
      for (const auto& [j, v] : off[i]) y[i] += v * x[j], y[j] += v * x[i];
    if (!ewald) return;
    const size_t w0 = size_t(2 * H[0] + 1), w1 = size_t(2 * H[1] + 1), w2 = size_t(2 * H[2] + 1);
    std::vector<std::complex<double>> z(n);
    for (const auto& k : ks) {
      std::complex<double> S = 0;
      for (size_t i = 0; i < n; ++i) {
        z[i] = ph[0][i * w0 + size_t(k.h + H[0])] * ph[1][i * w1 + size_t(k.k + H[1])] * ph[2][i * w2 + size_t(k.l + H[2])];
        S += z[i] * x[i];
      }
      const double wc = k.w * kCoul;
      for (size_t i = 0; i < n; ++i) y[i] += wc * (z[i].real() * S.real() + z[i].imag() * S.imag());
    }
  }
};

// Jacobi-preconditioned conjugate gradients for H x = b over the free atoms; returns the iterations.
int cg(const Operator& A, const std::vector<char>& free, const std::vector<double>& b, std::vector<double>& x, double tol, int maxit, double& res) {
  const size_t n = b.size();
  std::vector<double> r(n, 0), z(n, 0), p(n, 0), q(n, 0);
  auto mask = [&](std::vector<double>& v) { for (size_t i = 0; i < n; ++i) if (!free[i]) v[i] = 0; };
  mask(x);
  A.apply(x, q);
  for (size_t i = 0; i < n; ++i) r[i] = free[i] ? b[i] - q[i] : 0;
  std::vector<double> bm = b;
  mask(bm);
  const double bn = std::sqrt(std::max(dotv(bm, bm), 1e-300));
  for (size_t i = 0; i < n; ++i) z[i] = r[i] / A.pre[i];
  p = z;
  double rz = dotv(r, z);
  int it = 0;
  for (; it < maxit; ++it) {
    res = std::sqrt(dotv(r, r)) / bn;
    if (res < tol) break;
    A.apply(p, q);
    mask(q);
    const double alpha = rz / dotv(p, q);
    for (size_t i = 0; i < n; ++i) x[i] += alpha * p[i], r[i] -= alpha * q[i];
    for (size_t i = 0; i < n; ++i) z[i] = r[i] / A.pre[i];
    const double rz2 = dotv(r, z);
    const double beta = rz2 / rz;
    rz = rz2;
    for (size_t i = 0; i < n; ++i) p[i] = z[i] + beta * p[i];
  }
  res = std::sqrt(dotv(r, r)) / bn;
  return it;
}

struct Built {
  Operator A;
  std::vector<double> chi, J, zeta;   // per atom; zeta the Slater exponent (bohr⁻¹)
  double rc = 0;
};

// χ, J and the operator for s under o (J on the diagonal unless coulomb_only).
Built build(const System& s, const QEqOptions& o, QEqReport& R, bool coulomb_only) {
  const size_t n = s.atoms.size();
  Built B;
  B.chi.resize(n), B.J.resize(n), B.zeta.resize(n);
  if (o.kernel != "ohno" && o.kernel != "slater" && o.kernel != "point") throw std::invalid_argument("QEq kernel ohno, slater or point, not " + o.kernel);
  const bool slater = o.kernel == "slater";
  if (o.lambda <= 0) throw std::invalid_argument("QEq λ must be positive");
  for (size_t i = 0; i < n; ++i) {
    double rad = 0;
    if (!qeq_parameters(s.atoms[i].element, B.chi[i], B.J[i], rad))
      throw std::invalid_argument(std::string("no QEq parameters for ") + element(s.atoms[i].element).symbol);
    B.zeta[i] = o.lambda * (2 * principal_quantum_number(s.atoms[i].element) + 1) / (2 * rad / kBohr);
  }
  const bool periodic = s.cell.valid() && s.cell.periodic[0] && s.cell.periodic[1] && s.cell.periodic[2];
  const bool ewald = o.ewald && periodic;
  if (o.ewald && !periodic) R.notes.push_back("Ewald needs a cell periodic along a, b and c: not used");
  if (ewald && !(o.ewald_accuracy > 0 && o.ewald_accuracy < 1)) throw std::invalid_argument("Ewald accuracy between 0 and 1");
  double rc = o.cutoff > 0 ? o.cutoff : 10.0;
  double widths[3] = {0, 0, 0};
  if (s.cell.valid()) {
    const double v = s.cell.volume();
    widths[0] = v / norm(cross(s.cell.b, s.cell.c)), widths[1] = v / norm(cross(s.cell.c, s.cell.a)), widths[2] = v / norm(cross(s.cell.a, s.cell.b));
    const double w = std::min({widths[0], widths[1], widths[2]});
    if (rc > w / 2 && !ewald) {   // with Ewald the real space runs over explicit images instead
      rc = w / 2;
      char b[120];
      std::snprintf(b, sizeof b, "cut-off reduced to %.2f Å, half the narrowest cell width", rc);
      R.notes.push_back(b);
    }
  }
  B.rc = rc;
  std::map<std::pair<int, int>, SlaterTable> tables;   // per element pair
  if (slater) {
    std::map<int, double> zeta_of;
    for (size_t i = 0; i < n; ++i) zeta_of[s.atoms[i].element] = B.zeta[i];
    for (auto a = zeta_of.begin(); a != zeta_of.end(); ++a)
      for (auto b = a; b != zeta_of.end(); ++b)
        tables[{a->first, b->first}] = slater_table(principal_quantum_number(a->first), a->second, principal_quantum_number(b->first), b->second, rc);
  }
  Operator& A = B.A;
  A.diag.assign(n, 0.0);
  if (!coulomb_only) A.diag = B.J;
  A.off.assign(n, {});
  const double alpha = ewald ? std::sqrt(-std::log(o.ewald_accuracy)) / rc : 0;
  const double rc2 = rc * rc;
  auto taper = [&](double r) {
    const double x = r / rc;
    return 1 + x * x * x * x * (-35 + x * (84 + x * (-70 + x * 20)));
  };
  auto kernel = [&](size_t i, size_t j, double r) {   // eV
    if (o.kernel == "ohno") {
      const double a = 0.5 * kCoul * (1 / B.J[i] + 1 / B.J[j]);
      return kCoul / std::sqrt(r * r + a * a);
    }
    if (slater) return tables.at(std::minmax(s.atoms[i].element, s.atoms[j].element)).at(r);
    return kCoul / r;
  };
  if (ewald) {
    // real space over every periodic image within rc (the cell may be narrower than 2 rc): images of the atoms binned in
    // fractional coordinates over [−f, 1 + f] per axis, f = rc / width, bins at least f wide so that neighbouring bins
    // hold every image within rc; each (i, j > i, image) once, an atom's own images onto the diagonal
    double fcut[3];
    int N[3], nb[3];
    double bw[3];
    for (int a = 0; a < 3; ++a) {
      fcut[a] = rc / widths[a];
      N[a] = int(std::ceil(fcut[a]));
      nb[a] = std::max(1, int(std::floor((1 + 2 * fcut[a]) / fcut[a])));
      bw[a] = (1 + 2 * fcut[a]) / nb[a];
    }
    struct Img { uint32_t j; Vec3 r; };
    std::vector<std::vector<Img>> bins(size_t(nb[0]) * size_t(nb[1]) * size_t(nb[2]));
    auto bidx = [&](int x, int y, int z) { return (size_t(x) * size_t(nb[1]) + size_t(y)) * size_t(nb[2]) + size_t(z); };
    std::vector<Vec3> fr(n);
    for (size_t i = 0; i < n; ++i) {
      fr[i] = s.cell.to_fractional(s.atoms[i].pos);
      for (int a = 0; a < 3; ++a) fr[i][size_t(a)] -= std::floor(fr[i][size_t(a)]);
    }
    auto binof = [&](const Vec3& f, int a) { return std::clamp(int((f[size_t(a)] + fcut[a]) / bw[a]), 0, nb[a] - 1); };
    for (uint32_t j = 0; j < n; ++j)
      for (int u = -N[0]; u <= N[0]; ++u)
        for (int v = -N[1]; v <= N[1]; ++v)
          for (int w = -N[2]; w <= N[2]; ++w) {
            const Vec3 f = fr[j] + Vec3{double(u), double(v), double(w)};
            bool in = true;
            for (int a = 0; a < 3; ++a) in = in && f[size_t(a)] >= -fcut[a] && f[size_t(a)] <= 1 + fcut[a];
            if (in) bins[bidx(binof(f, 0), binof(f, 1), binof(f, 2))].push_back({j, s.cell.to_cartesian(f) - s.cell.origin});
          }
    std::vector<double> acc(n, 0.0);
    std::vector<uint32_t> touched;
    for (uint32_t i = 0; i < n; ++i) {
      const Vec3 ri = s.cell.to_cartesian(fr[i]) - s.cell.origin;
      const int bx = binof(fr[i], 0), by = binof(fr[i], 1), bz = binof(fr[i], 2);
      touched.clear();
      for (int x = std::max(0, bx - 1); x <= std::min(nb[0] - 1, bx + 1); ++x)
        for (int y = std::max(0, by - 1); y <= std::min(nb[1] - 1, by + 1); ++y)
          for (int z = std::max(0, bz - 1); z <= std::min(nb[2] - 1, bz + 1); ++z)
            for (const auto& im : bins[bidx(x, y, z)]) {
              if (im.j < i) continue;
              const Vec3 d = im.r - ri;
              const double r2 = dot(d, d);
              if (r2 >= rc2 || r2 < 1e-12) continue;   // r ≈ 0: the atom itself
              const double r = std::sqrt(r2);
              const double v = kCoul * std::erfc(alpha * r) / r + (kernel(i, im.j, r) - kCoul / r) * (o.kernel == "ohno" ? taper(r) : 1.0);
              if (acc[im.j] == 0.0) touched.push_back(im.j);
              acc[im.j] += v;
              if (acc[im.j] == 0.0) acc[im.j] = 1e-300;   // keep it marked
            }
      for (uint32_t j : touched) {
        if (j == i) A.diag[i] += acc[j];
        else A.off[i].push_back({j, acc[j]});
        acc[j] = 0.0;
      }
    }
  }
  const Grid g(s, rc);
  for (uint32_t i = 0; i < n && !ewald; ++i) {
    const Vec3& f = g.frac[i];
    g.for_neighbour_bins(g.bin(f, 0), g.bin(f, 1), g.bin(f, 2), [&](const std::vector<uint32_t>& bin) {
      for (uint32_t j : bin) {
        if (j <= i) continue;
        const Vec3 d = g.sep(i, j);
        const double r2 = dot(d, d);
        if (r2 >= rc2) continue;
        const double r = std::sqrt(r2);
        A.off[i].push_back({j, taper(r) * kernel(i, j, r)});
      }
    });
  }
  if (ewald) {
    A.ewald = true;
    const double V = s.cell.volume();
    const double kmax = 2 * alpha * std::sqrt(-std::log(o.ewald_accuracy));
    // k = h b1 + k b2 + l b3 gives k·a = 2πh, so |k| ≤ kmax needs |h| ≤ kmax |a| / 2π (the lattice length, not the plane spacing)
    const double lens[3] = {norm(s.cell.a), norm(s.cell.b), norm(s.cell.c)};
    for (int k = 0; k < 3; ++k) A.H[k] = int(std::ceil(kmax * lens[k] / (2 * M_PI)));
    const Vec3 b1 = cross(s.cell.b, s.cell.c) * (2 * M_PI / V), b2 = cross(s.cell.c, s.cell.a) * (2 * M_PI / V), b3 = cross(s.cell.a, s.cell.b) * (2 * M_PI / V);
    double self = 0;
    for (int h = 0; h <= A.H[0]; ++h)
      for (int k = -A.H[1]; k <= A.H[1]; ++k)
        for (int l = -A.H[2]; l <= A.H[2]; ++l) {
          if (h == 0 && (k < 0 || (k == 0 && l <= 0))) continue;
          const Vec3 kv = b1 * double(h) + b2 * double(k) + b3 * double(l);
          const double k2 = dot(kv, kv);
          if (k2 > kmax * kmax) continue;
          const double w = 2 * (4 * M_PI / V) * std::exp(-k2 / (4 * alpha * alpha)) / k2;   // ×2: k and −k
          A.ks.push_back({h, k, l, w});
          self += w;
        }
    for (int a = 0; a < 3; ++a) {
      const size_t wd = size_t(2 * A.H[a] + 1);
      A.ph[a].resize(n * wd);
      for (size_t i = 0; i < n; ++i) {
        const Vec3 fr = s.cell.to_fractional(s.atoms[i].pos);   // unwrapped is fine: integer h
        for (int h = -A.H[a]; h <= A.H[a]; ++h) {
          const double t = 2 * M_PI * h * fr[size_t(a)];
          A.ph[a][i * wd + size_t(h + A.H[a])] = {std::cos(t), std::sin(t)};
        }
      }
    }
    for (size_t i = 0; i < n; ++i) A.diag[i] -= kCoul * 2 * alpha / std::sqrt(M_PI);   // apply() holds the reciprocal i = j term
    A.pre = A.diag;
    for (size_t i = 0; i < n; ++i) A.pre[i] += kCoul * self;
    R.alpha = alpha;
    R.kvectors = int(A.ks.size());
  }
  if (A.pre.empty()) A.pre = A.diag;
  return B;
}

}  // namespace

QEqOptions qeq_options_from_json(const Json& j) {
  QEqOptions o;
  o.kernel = j.text("kernel", o.kernel);
  o.ewald = j.has("ewald") && j["ewald"].kind() == Json::Bool && j["ewald"].boolean();
  o.ewald_accuracy = j.num("accuracy", o.ewald_accuracy);
  o.lambda = j.num("lambda", o.lambda);
  o.hydrogen = j.has("hydrogen") && j["hydrogen"].kind() == Json::Bool && j["hydrogen"].boolean();
  o.cutoff = j.num("cutoff", o.cutoff);
  o.total_charge = j.num("total", o.total_charge);
  if (j.has("bounds") && j["bounds"].is_object())
    for (const auto& [sym, v] : j["bounds"].members()) {
      const int z = element_from_symbol(sym);
      if (z <= 0) throw std::invalid_argument("QEq bounds: unknown element " + sym);
      const auto lim = v.items();
      if (lim.size() != 2) throw std::invalid_argument("QEq bounds of " + sym + ": [minimum, maximum]");
      o.bounds[z] = {lim[0].number(), lim[1].number()};
    }
  return o;
}

std::map<int, std::pair<double, double>> qeq_bounds_from_text(const std::string& t) {
  std::map<int, std::pair<double, double>> b;
  size_t at = 0;
  while (at < t.size()) {
    size_t e = t.find(',', at);
    if (e == std::string::npos) e = t.size();
    std::string item = t.substr(at, e - at);
    at = e + 1;
    item.erase(0, item.find_first_not_of(" \t"));
    item.erase(item.find_last_not_of(" \t") + 1);
    if (item.empty()) continue;
    const size_t c1 = item.find(':'), c2 = c1 == std::string::npos ? c1 : item.find(':', c1 + 1);
    if (c2 == std::string::npos) throw std::invalid_argument("QEq bounds 'O:-1.2:0, Si:0:2.4'");
    const int z = element_from_symbol(item.substr(0, c1));
    if (z <= 0) throw std::invalid_argument("QEq bounds: unknown element " + item.substr(0, c1));
    b[z] = {std::stod(item.substr(c1 + 1, c2 - c1 - 1)), std::stod(item.substr(c2 + 1))};
  }
  return b;
}

int principal_quantum_number(int z) {
  return z <= 2 ? 1 : z <= 10 ? 2 : z <= 18 ? 3 : z <= 36 ? 4 : z <= 54 ? 5 : z <= 86 ? 6 : 7;
}

double slater_coulomb_integral(int n1, double z1, int n2, double z2, double R) {
  if (n1 < 1 || n2 < 1 || !(z1 > 0) || !(z2 > 0)) throw std::invalid_argument("Slater integral: n ≥ 1 and ζ > 0");
  static const auto rule = [] {
    std::pair<std::vector<double>, std::vector<double>> g;
    gauss_legendre(16, g.first, g.second);
    return g;
  }();
  const auto& gx = rule.first;
  const auto& gw = rule.second;
  // J(R) = ∫ 4πr² ρ1(r) ⟨Φ2⟩(r) dr, ⟨Φ2⟩ the average of Φ2(|r − R|) over the sphere of radius r: [G(r+R) − G(|r−R|)] / (2rR)
  const double lnorm = (2 * n1 + 1) * std::log(2 * z1) - std::lgamma(2.0 * n1 + 1);
  auto w1 = [&](double r) { return std::exp(lnorm + 2 * n1 * std::log(r) - 2 * z1 * r); };
  auto avg = [&](double r) {
    if (R < 1e-6) return slater_phi(n2, z2, r);
    return (slater_g(n2, z2, r + R) - slater_g(n2, z2, std::fabs(r - R))) / (2 * r * R);
  };
  auto panel = [&](double a, double b) {
    double s = 0;
    const double c = 0.5 * (a + b), h = 0.5 * (b - a);
    for (size_t k = 0; k < gx.size(); ++k) {
      const double r = c + h * gx[k];
      s += gw[k] * w1(r) * avg(r);
    }
    return s * h;
  };
  const double tail = (2 * n1 + 60) / (2 * z1);   // 4πr²ρ1 below e^(−55) of its peak beyond
  const double edge = std::min(R, tail);          // the kink of ⟨Φ2⟩ at r = R sits on a panel boundary
  double J = 0;
  const int p1 = 12, p2 = 48;
  if (edge > 0)
    for (int k = 0; k < p1; ++k) J += panel(edge * k / p1, edge * (k + 1) / p1);
  for (int k = 0; k < p2; ++k) J += panel(edge + tail * k / p2, edge + tail * (k + 1) / p2);
  return J;
}

double qeq_coulomb_energy(const System& s, const std::vector<double>& q, const QEqOptions& o) {
  if (q.size() != s.atoms.size()) throw std::invalid_argument("one charge per atom");
  QEqReport R;
  const Built B = build(s, o, R, true);
  std::vector<double> y(q.size());
  B.A.apply(q, y);
  return 0.5 * dotv(q, y);
}

std::vector<double> qeq_charges(const System& s, const QEqOptions& o, QEqReport* report) {
  const size_t n = s.atoms.size();
  QEqReport R;
  if (n == 0) return {};
  Built B = build(s, o, R, false);
  Operator& A = B.A;
  R.cutoff = B.rc;
  const std::vector<double> D0 = A.diag;   // the diagonal before the hydrogen correction
  std::vector<int> hyd;
  if (o.hydrogen)
    for (size_t i = 0; i < n; ++i)
      if (s.atoms[i].element == 1) hyd.push_back(int(i));
  std::vector<double> lo(n, -1e300), hi(n, 1e300);
  for (size_t i = 0; i < n; ++i) {
    const auto it = o.bounds.find(s.atoms[i].element);
    if (it == o.bounds.end()) continue;
    if (it->second.first > it->second.second)
      throw std::invalid_argument(std::string("QEq bounds of ") + element(s.atoms[i].element).symbol + ": the minimum is above the maximum");
    lo[i] = it->second.first, hi[i] = it->second.second;
  }
  std::vector<char> pinned(n, 0);
  std::vector<double> q(n, 0.0), qn(n), y(n), xs(n, 0.0), xt(n, 0.0);
  const bool iterate = !hyd.empty() || !o.bounds.empty();
  int iters = 0, pass = 0;
  double res = 0, mu = 0, change = 0;
  bool set_changed = false;
  // the free atoms equalised at the total with the held ones fixed: qn (held atoms keep q), mu
  auto equalise = [&]() {
    std::vector<char> free(n);
    std::vector<double> qp(n, 0.0);
    double pinned_sum = 0;
    for (size_t i = 0; i < n; ++i) {
      free[i] = !pinned[i];
      if (pinned[i]) qp[i] = q[i], pinned_sum += q[i];
    }
    A.apply(qp, y);
    std::vector<double> bs(n), bt(n, 1.0);
    for (size_t i = 0; i < n; ++i) bs[i] = -B.chi[i] - y[i];
    double r1 = 0, r2 = 0;
    iters += cg(A, free, bs, xs, o.tolerance, o.max_iterations, r1);
    iters += cg(A, free, bt, xt, o.tolerance, o.max_iterations, r2);
    res = std::max(r1, r2);
    double ss = 0, st = 0;
    for (size_t i = 0; i < n; ++i)
      if (free[i]) ss += xs[i], st += xt[i];
    if (st == 0) throw std::invalid_argument("QEq: every atom is held at a bound, so the total charge cannot be met");
    mu = (o.total_charge - pinned_sum - ss) / st;
    for (size_t i = 0; i < n; ++i) qn[i] = free[i] ? xs[i] + mu * xt[i] : q[i];
  };
  for (; pass < (iterate ? std::max(1, o.scf_iterations) : 1); ++pass) {
    for (int i : hyd)   // J_H(q) = J_H⁰ (1 + q/ζ_H), kept positive
    {
      const double dJ = B.J[size_t(i)] * std::max(-0.9, q[size_t(i)] / B.zeta[size_t(i)]);
      A.pre[size_t(i)] += D0[size_t(i)] + dJ - A.diag[size_t(i)];
      A.diag[size_t(i)] = D0[size_t(i)] + dJ;
    }
    equalise();
    const double mix = (!hyd.empty() && pass > 0 && !set_changed) ? 0.5 : 1.0;   // hydrogen: damped steps after the first
    change = 0;
    for (size_t i = 0; i < n; ++i) {
      const double v = q[i] + mix * (qn[i] - q[i]);
      change = std::max(change, std::fabs(v - q[i]));
      q[i] = v;
    }
    if (!iterate) break;
    // bounds: hold atoms beyond a bound there; release held atoms the Karush–Kuhn–Tucker conditions no longer keep
    set_changed = false;
    if (!o.bounds.empty()) {
      std::vector<double> gq(n);
      A.apply(q, gq);
      for (size_t i = 0; i < n; ++i) {
        const double gi = B.chi[i] + gq[i];   // ∂E/∂q_i; free atoms share μ
        if (!pinned[i] && (q[i] > hi[i] + 1e-10 || q[i] < lo[i] - 1e-10)) {
          q[i] = std::clamp(q[i], lo[i], hi[i]);
          pinned[i] = 1, set_changed = true;
        } else if (pinned[i] && ((q[i] >= hi[i] && gi - mu > 1e-8) || (q[i] <= lo[i] && gi - mu < -1e-8))) {
          pinned[i] = 0, set_changed = true;
        }
      }
    }
    if (!set_changed && pass > 0 && change < o.scf_tolerance) break;
    if (!set_changed && hyd.empty()) break;   // bounds alone: the active set is settled
  }
  // the pass limit reached on a pass that held or released atoms: the clamped charges no longer add up to the total,
  // so equalise the free atoms once more with the final held set
  const bool unsettled = set_changed;
  if (unsettled) {
    equalise();
    q = qn;
  }
  R.iterations = iters;
  R.residual = res;
  R.mu = mu;   // H q = −χ + μ·1, so χ_i + (H q)_i = ∂E/∂q_i = μ: the equalised electronegativity is μ
  R.scf = std::min(pass + 1, std::max(1, o.scf_iterations));
  for (size_t i = 0; i < n; ++i) R.bounded += pinned[i];
  char b[320];
  std::snprintf(b, sizeof b, "QEq (Rappé & Goddard 1991), %s kernel%s, cut-off %.1f Å · %d CG iterations · residual %.1e · total %.4f e",
                o.kernel == "ohno" ? "Ohno–Klopman" : o.kernel == "slater" ? "ns Slater" : "point-charge", R.alpha > 0 ? ", Ewald sum" : "", B.rc,
                R.iterations, R.residual, o.total_charge);
  R.notes.insert(R.notes.begin(), b);
  if (R.alpha > 0) {
    std::snprintf(b, sizeof b, "Ewald α %.4f Å⁻¹, %d k-vectors (half space), accuracy %.0e", R.alpha, R.kvectors, o.ewald_accuracy);
    R.notes.push_back(b);
  }
  if (o.kernel == "slater") {
    std::snprintf(b, sizeof b, "Slater exponents ζ = λ(2n+1)/(2R), λ = %.4f", o.lambda);
    R.notes.push_back(b);
  }
  if (!hyd.empty()) R.notes.push_back(std::to_string(hyd.size()) + " hydrogen(s) with the charge-dependent hardness J_H(q) = J_H⁰(1 + q/ζ_H)");
  if (iterate) {
    if (!hyd.empty()) std::snprintf(b, sizeof b, "%d self-consistent pass(es), last change %.1e e", R.scf, change);
    else std::snprintf(b, sizeof b, "%d pass(es) to settle the atoms held at bounds", R.scf);
    R.notes.push_back(b);
    if (!hyd.empty() && change >= o.scf_tolerance) R.notes.push_back("the self-consistent passes did not converge");
  }
  if (R.bounded) R.notes.push_back(std::to_string(R.bounded) + " atom(s) held at a charge bound");
  if (unsettled) R.notes.push_back("the atoms held at bounds did not settle within the pass limit; free atoms re-equalised to the total once more");
  if (R.residual > 1e-5) R.notes.push_back("conjugate gradients did not converge fully");
  if (report) *report = R;
  return q;
}

}  // namespace caps
