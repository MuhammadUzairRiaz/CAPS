// CAPS QEq charge equilibration (see caps/qeq.hpp).
#include "caps/qeq.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>

#include "caps/elements.hpp"
#include "cell_list.hpp"

namespace caps {

namespace {

constexpr double kCoul = 14.399645;   // e²/(4πε0) in eV Å

struct Sparse {
  std::vector<double> diag;
  std::vector<std::vector<std::pair<uint32_t, double>>> off;   // each pair stored once (j > i)
  void mul(const std::vector<double>& x, std::vector<double>& y) const {
    const size_t n = diag.size();
    for (size_t i = 0; i < n; ++i) y[i] = diag[i] * x[i];
    for (size_t i = 0; i < n; ++i)
      for (const auto& [j, v] : off[i]) y[i] += v * x[j], y[j] += v * x[i];
  }
};

double dotv(const std::vector<double>& a, const std::vector<double>& b) {
  double s = 0;
  for (size_t i = 0; i < a.size(); ++i) s += a[i] * b[i];
  return s;
}

// Jacobi-preconditioned conjugate gradients for H x = b; returns iterations, leaves the relative residual in res.
int cg(const Sparse& H, const std::vector<double>& b, std::vector<double>& x, double tol, int maxit, double& res) {
  const size_t n = b.size();
  std::vector<double> r(n), z(n), p(n), q(n);
  H.mul(x, q);
  for (size_t i = 0; i < n; ++i) r[i] = b[i] - q[i];
  const double bn = std::sqrt(std::max(dotv(b, b), 1e-300));
  for (size_t i = 0; i < n; ++i) z[i] = r[i] / H.diag[i];
  p = z;
  double rz = dotv(r, z);
  int it = 0;
  for (; it < maxit; ++it) {
    res = std::sqrt(dotv(r, r)) / bn;
    if (res < tol) break;
    H.mul(p, q);
    const double alpha = rz / dotv(p, q);
    for (size_t i = 0; i < n; ++i) x[i] += alpha * p[i], r[i] -= alpha * q[i];
    for (size_t i = 0; i < n; ++i) z[i] = r[i] / H.diag[i];
    const double rz2 = dotv(r, z);
    const double beta = rz2 / rz;
    rz = rz2;
    for (size_t i = 0; i < n; ++i) p[i] = z[i] + beta * p[i];
  }
  res = std::sqrt(dotv(r, r)) / bn;
  return it;
}

}  // namespace

std::vector<double> qeq_charges(const System& s, const QEqOptions& o, QEqReport* report) {
  const size_t n = s.atoms.size();
  QEqReport R;
  if (n == 0) return {};
  std::vector<double> chi(n), J(n);
  for (size_t i = 0; i < n; ++i) {
    double r;
    if (!qeq_parameters(s.atoms[i].element, chi[i], J[i], r))
      throw std::invalid_argument(std::string("no QEq parameters for ") + element(s.atoms[i].element).symbol);
  }
  double rc = o.cutoff > 0 ? o.cutoff : 10.0;
  if (s.cell.valid()) {
    const double v = s.cell.volume();
    const double w = std::min({v / norm(cross(s.cell.b, s.cell.c)), v / norm(cross(s.cell.c, s.cell.a)), v / norm(cross(s.cell.a, s.cell.b))});
    if (rc > w / 2) {
      rc = w / 2;
      char b[120];
      std::snprintf(b, sizeof b, "cut-off reduced to %.2f Å, half the narrowest cell width", rc);
      R.notes.push_back(b);
    }
  }
  R.cutoff = rc;
  // the matrix: J on the diagonal, the tapered shielded Coulomb kernel off it
  Sparse H;
  H.diag = J;
  H.off.assign(n, {});
  const Grid g(s, rc);
  const double rc2 = rc * rc;
  auto taper = [&](double r) {
    const double x = r / rc;
    return 1 + x * x * x * x * (-35 + x * (84 + x * (-70 + x * 20)));
  };
  for (uint32_t i = 0; i < n; ++i) {
    const Vec3& f = g.frac[i];
    g.for_neighbour_bins(g.bin(f, 0), g.bin(f, 1), g.bin(f, 2), [&](const std::vector<uint32_t>& bin) {
      for (uint32_t j : bin) {
        if (j <= i) continue;
        const Vec3 d = g.sep(i, j);
        const double r2 = dot(d, d);
        if (r2 >= rc2) continue;
        const double a = 0.5 * kCoul * (1 / J[i] + 1 / J[j]);
        const double r = std::sqrt(r2);
        H.off[i].push_back({j, taper(r) * kCoul / std::sqrt(r2 + a * a)});
      }
    });
  }
  std::vector<double> bs(n), bt(n, 1.0), xs(n, 0.0), xt(n, 0.0);
  for (size_t i = 0; i < n; ++i) bs[i] = -chi[i];
  double r1 = 0, r2 = 0;
  const int i1 = cg(H, bs, xs, o.tolerance, o.max_iterations, r1);
  const int i2 = cg(H, bt, xt, o.tolerance, o.max_iterations, r2);
  double ss = 0, st = 0;
  for (size_t i = 0; i < n; ++i) ss += xs[i], st += xt[i];
  const double mu = (o.total_charge - ss) / st;
  std::vector<double> q(n);
  for (size_t i = 0; i < n; ++i) q[i] = xs[i] + mu * xt[i];
  R.iterations = i1 + i2;
  R.residual = std::max(r1, r2);
  char b[200];
  std::snprintf(b, sizeof b, "QEq (Rappé & Goddard 1991), Ohno–Klopman shielding, cut-off %.1f Å · %d CG iterations · residual %.1e · total %.4f e", rc, R.iterations,
                R.residual, o.total_charge);
  R.notes.insert(R.notes.begin(), b);
  if (R.residual > 1e-5) R.notes.push_back("conjugate gradients did not converge fully");
  if (report) *report = R;
  return q;
}

}  // namespace caps
