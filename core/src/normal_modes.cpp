#include "caps/normal_modes.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

#include "caps/uff.hpp"

namespace caps {

// ---------------------------------------------------------------- symmetric eigenproblem

void symmetric_eigen(std::vector<double>& a, size_t nn, std::vector<double>& w) {
  const int n = int(nn);
  auto V = [&](int i, int j) -> double& { return a[size_t(i) * nn + size_t(j)]; };
  std::vector<double> d(nn), e(nn);
  w.assign(nn, 0);
  if (n == 0) return;
  // Householder reduction to tridiagonal form (the transformations accumulated in V)
  for (int j = 0; j < n; ++j) d[j] = V(n - 1, j);
  for (int i = n - 1; i > 0; --i) {
    double scale = 0, h = 0;
    for (int k = 0; k < i; ++k) scale += std::fabs(d[k]);
    if (scale == 0) {
      e[i] = d[i - 1];
      for (int j = 0; j < i; ++j) d[j] = V(i - 1, j), V(i, j) = 0, V(j, i) = 0;
    } else {
      for (int k = 0; k < i; ++k) d[k] /= scale, h += d[k] * d[k];
      double f = d[i - 1], g = std::sqrt(h);
      if (f > 0) g = -g;
      e[i] = scale * g;
      h -= f * g;
      d[i - 1] = f - g;
      for (int j = 0; j < i; ++j) e[j] = 0;
      for (int j = 0; j < i; ++j) {
        f = d[j];
        V(j, i) = f;
        g = e[j] + V(j, j) * f;
        for (int k = j + 1; k <= i - 1; ++k) g += V(k, j) * d[k], e[k] += V(k, j) * f;
        e[j] = g;
      }
      f = 0;
      for (int j = 0; j < i; ++j) e[j] /= h, f += e[j] * d[j];
      const double hh = f / (h + h);
      for (int j = 0; j < i; ++j) e[j] -= hh * d[j];
      for (int j = 0; j < i; ++j) {
        f = d[j], g = e[j];
        for (int k = j; k <= i - 1; ++k) V(k, j) -= f * e[k] + g * d[k];
        d[j] = V(i - 1, j);
        V(i, j) = 0;
      }
    }
    d[i] = h;
  }
  for (int i = 0; i < n - 1; ++i) {
    V(n - 1, i) = V(i, i);
    V(i, i) = 1;
    const double h = d[i + 1];
    if (h != 0) {
      for (int k = 0; k <= i; ++k) d[k] = V(k, i + 1) / h;
      for (int j = 0; j <= i; ++j) {
        double g = 0;
        for (int k = 0; k <= i; ++k) g += V(k, i + 1) * V(k, j);
        for (int k = 0; k <= i; ++k) V(k, j) -= g * d[k];
      }
    }
    for (int k = 0; k <= i; ++k) V(k, i + 1) = 0;
  }
  for (int j = 0; j < n; ++j) d[j] = V(n - 1, j), V(n - 1, j) = 0;
  V(n - 1, n - 1) = 1;
  e[0] = 0;
  // implicit QL on the tridiagonal matrix
  for (int i = 1; i < n; ++i) e[i - 1] = e[i];
  e[n - 1] = 0;
  double f = 0, tst1 = 0;
  const double eps = std::ldexp(1.0, -52);
  for (int l = 0; l < n; ++l) {
    tst1 = std::max(tst1, std::fabs(d[l]) + std::fabs(e[l]));
    int m = l;
    while (m < n && std::fabs(e[m]) > eps * tst1) ++m;
    if (m == n) m = n - 1;
    if (m > l) {
      int iter = 0;
      do {
        if (++iter > 200) throw std::runtime_error("eigenvalues did not converge");
        double g = d[l];
        double p = (d[l + 1] - g) / (2 * e[l]);
        double r = std::hypot(p, 1.0);
        if (p < 0) r = -r;
        d[l] = e[l] / (p + r);
        d[l + 1] = e[l] * (p + r);
        const double dl1 = d[l + 1];
        double h = g - d[l];
        for (int i = l + 2; i < n; ++i) d[i] -= h;
        f += h;
        p = d[m];
        double c = 1, c2 = 1, c3 = 1, s = 0, s2 = 0;
        const double el1 = e[l + 1];
        for (int i = m - 1; i >= l; --i) {
          c3 = c2, c2 = c, s2 = s;
          g = c * e[i];
          h = c * p;
          r = std::hypot(p, e[i]);
          e[i + 1] = s * r;
          s = e[i] / r;
          c = p / r;
          p = c * d[i] - s * g;
          d[i + 1] = h + s * (c * g + s * d[i]);
          for (int k = 0; k < n; ++k) {
            h = V(k, i + 1);
            V(k, i + 1) = s * V(k, i) + c * h;
            V(k, i) = c * V(k, i) - s * h;
          }
        }
        p = -s * s2 * c3 * el1 * e[l] / dl1;
        e[l] = s * p;
        d[l] = c * p;
      } while (std::fabs(e[l]) > eps * tst1);
    }
    d[l] += f;
    e[l] = 0;
  }
  // ascending, the vectors with them
  std::vector<int> ord(nn);
  std::iota(ord.begin(), ord.end(), 0);
  std::sort(ord.begin(), ord.end(), [&](int x, int y) { return d[x] < d[y]; });
  std::vector<double> b(a.size());
  for (int j = 0; j < n; ++j) {
    w[j] = d[ord[j]];
    for (int i = 0; i < n; ++i) b[size_t(i) * nn + j] = V(i, ord[j]);
  }
  a.swap(b);
}

// ---------------------------------------------------------------- normal modes

namespace {
constexpr double kWavenumber = 108.5913;      // cm⁻¹ per √(kcal/mol/Å²/(g/mol)): √(4.184e-4 fs⁻²) / (2π c)
constexpr double kCm1ToKcal = 2.859144e-3;    // kcal/mol per cm⁻¹ (h c N_A)
constexpr double kC2 = 1.4387769;             // cm·K, h c / k_B
constexpr double kR = 1.98720426;             // cal/mol/K
}  // namespace

NormalModesResult normal_modes(const System& s, const NormalModesOptions& o) {
  NormalModesResult r;
  const size_t n = s.atoms.size();
  if (n < 2) throw std::invalid_argument("normal modes need at least two atoms");
  std::vector<char> mv = o.moving.size() == n ? o.moving : std::vector<char>(n, 1);
  std::vector<uint32_t> idx;
  for (size_t i = 0; i < n; ++i)
    if (mv[i]) idx.push_back(uint32_t(i));
  if (idx.empty()) throw std::invalid_argument("no atoms to vibrate");
  if (idx.size() > o.max_atoms)
    throw std::invalid_argument(std::to_string(idx.size()) + " atoms vibrate: the dense Hessian takes at most " + std::to_string(o.max_atoms) +
                                " (select a molecule or a region, the rest held)");
  const ForceField ff = o.field ? *o.field : default_forcefield(s);
  r.field = ff.name;
  const size_t m = idx.size(), dim = 3 * m;
  Evaluator ev(ff, o.energy);
  std::vector<double> x(3 * n), f;
  for (size_t i = 0; i < n; ++i)
    for (int k = 0; k < 3; ++k) x[3 * i + k] = s.atoms[i].pos[k];
  ev.compute(x, s.cell, f);
  for (uint32_t i : idx) r.max_force = std::max(r.max_force, std::sqrt(f[3 * i] * f[3 * i] + f[3 * i + 1] * f[3 * i + 1] + f[3 * i + 2] * f[3 * i + 2]));
  ev.freeze_pairs(x, s.cell);
  // Hessian of the moving coordinates: H_ab = −∂F_a/∂x_b by central differences
  std::vector<double> H(dim * dim);
  std::vector<double> fp, fm;
  for (size_t b = 0; b < dim; ++b) {
    if (o.progress && !o.progress(double(b) / dim)) throw std::runtime_error("normal modes cancelled");
    const size_t xb = 3 * idx[b / 3] + b % 3;
    const double x0 = x[xb];
    x[xb] = x0 + o.step;
    ev.compute(x, s.cell, fp);
    x[xb] = x0 - o.step;
    ev.compute(x, s.cell, fm);
    x[xb] = x0;
    for (size_t a = 0; a < dim; ++a) {
      const size_t xa = 3 * idx[a / 3] + a % 3;
      H[a * dim + b] = -(fp[xa] - fm[xa]) / (2 * o.step);
    }
  }
  ev.unfreeze();
  std::vector<double> sm(m);
  for (size_t k = 0; k < m; ++k) sm[k] = std::sqrt(ff.mass[idx[k]] > 0 ? ff.mass[idx[k]] : 1.0);
  for (size_t a = 0; a < dim; ++a)
    for (size_t b = a; b < dim; ++b) {
      const double v = 0.5 * (H[a * dim + b] + H[b * dim + a]) / (sm[a / 3] * sm[b / 3]);
      H[a * dim + b] = H[b * dim + a] = v;
    }
  // rigid motions, in mass-weighted coordinates: translations always; rotations without a periodic cell; none when
  // some atoms are held (the held ones pin the rest)
  std::vector<std::vector<double>> tr;
  if (m == n) {
    for (int ax = 0; ax < 3; ++ax) {
      std::vector<double> v(dim, 0);
      for (size_t k = 0; k < m; ++k) v[3 * k + ax] = sm[k];
      tr.push_back(v);
    }
    if (!s.cell.valid()) {
      double M = 0;
      Vec3 c{0, 0, 0};
      for (size_t k = 0; k < m; ++k) {
        const double mk = sm[k] * sm[k];
        M += mk;
        for (int q = 0; q < 3; ++q) c[q] += mk * x[3 * idx[k] + q];
      }
      for (int q = 0; q < 3; ++q) c[q] /= M;
      for (int ax = 0; ax < 3; ++ax) {
        std::vector<double> v(dim, 0);
        for (size_t k = 0; k < m; ++k) {
          const double rx = x[3 * idx[k]] - c[0], ry = x[3 * idx[k] + 1] - c[1], rz = x[3 * idx[k] + 2] - c[2];
          const double cr[3][3] = {{0, -rz, ry}, {rz, 0, -rx}, {-ry, rx, 0}};   // e_ax × r, column ax
          for (int q = 0; q < 3; ++q) v[3 * k + q] = sm[k] * cr[q][ax];
        }
        tr.push_back(v);
      }
    }
  }
  // orthonormalise (a linear molecule loses one rotation)
  std::vector<std::vector<double>> basis;
  for (auto v : tr) {
    for (const auto& u : basis) {
      const double p = std::inner_product(v.begin(), v.end(), u.begin(), 0.0);
      for (size_t q = 0; q < dim; ++q) v[q] -= p * u[q];
    }
    const double nv = std::sqrt(std::inner_product(v.begin(), v.end(), v.begin(), 0.0));
    if (nv > 1e-6 * std::sqrt(double(dim))) {
      for (double& q : v) q /= nv;
      basis.push_back(v);
    }
  }
  r.projected = int(basis.size());
  if (!basis.empty()) {   // H ← P H P with P = 1 − Σ v vᵀ
    std::vector<double> T(dim * dim);
    auto project_rows = [&](std::vector<double>& A) {   // A ← A P
      for (size_t a = 0; a < dim; ++a) {
        double* row = &A[a * dim];
        for (const auto& u : basis) {
          const double p = std::inner_product(row, row + dim, u.begin(), 0.0);
          for (size_t q = 0; q < dim; ++q) row[q] -= p * u[q];
        }
      }
    };
    project_rows(H);
    for (size_t a = 0; a < dim; ++a)
      for (size_t b = 0; b < dim; ++b) T[b * dim + a] = H[a * dim + b];
    project_rows(T);   // (H P)ᵀ P = P H P (H symmetric)
    for (size_t a = 0; a < dim; ++a)
      for (size_t b = a; b < dim; ++b) H[a * dim + b] = H[b * dim + a] = 0.5 * (T[a * dim + b] + T[b * dim + a]);
  }
  std::vector<double> lam;
  symmetric_eigen(H, dim, lam);
  // drop the projected rigid motions: the eigenvectors lying in their space
  std::vector<char> rigid(dim, 0);
  if (!basis.empty()) {
    std::vector<std::pair<double, size_t>> ov;
    for (size_t j = 0; j < dim; ++j) {
      double o2 = 0;
      for (const auto& u : basis) {
        double p = 0;
        for (size_t q = 0; q < dim; ++q) p += H[q * dim + j] * u[q];
        o2 += p * p;
      }
      ov.push_back({o2, j});
    }
    std::sort(ov.rbegin(), ov.rend());
    for (size_t k = 0; k < basis.size(); ++k) rigid[ov[k].second] = 1;
  }
  bool charged = false;
  for (uint32_t i : idx) charged = charged || std::fabs(ff.charge.size() == n ? ff.charge[i] : 0.0) > 1e-9;
  double imax = 0;
  for (size_t j = 0; j < dim; ++j) {
    if (rigid[j]) continue;
    const double l = lam[j];
    const double nu = (l >= 0 ? 1 : -1) * kWavenumber * std::sqrt(std::fabs(l));
    std::vector<double> mode(3 * n, 0);
    double inv_mu = 0, norm = 0;
    Vec3 dmu{0, 0, 0};
    for (size_t a = 0; a < dim; ++a) {
      const double q = H[a * dim + j], c = q / sm[a / 3];
      mode[3 * idx[a / 3] + a % 3] = c;
      inv_mu += c * c;
      norm += c * c;
      if (charged) dmu[a % 3] += ff.charge[idx[a / 3]] * c;
    }
    norm = std::sqrt(norm);
    for (double& c : mode) c /= norm;
    r.wavenumber.push_back(nu);
    r.mode.push_back(std::move(mode));
    r.reduced_mass.push_back(1 / inv_mu);
    const double I = dmu[0] * dmu[0] + dmu[1] * dmu[1] + dmu[2] * dmu[2];
    r.ir.push_back(I);
    imax = std::max(imax, I);
    if (nu < -10) ++r.imaginary;
  }
  for (double& I : r.ir) I = imax > 0 ? I / imax : 0;
  // quantum harmonic thermodynamics of the real modes above 10 cm⁻¹
  const double T = o.temperature;
  int soft = 0;
  for (double nu : r.wavenumber) {
    if (nu < 10) { soft += nu > -10; continue; }
    const double xq = kC2 * nu / T, ex = std::exp(xq);
    r.zpe += 0.5 * kCm1ToKcal * nu;
    r.e_vib += kCm1ToKcal * nu * (0.5 + 1 / (ex - 1));
    r.s_vib += kR * (xq / (ex - 1) - std::log(1 - 1 / ex));
    r.cv_vib += kR * xq * xq * ex / ((ex - 1) * (ex - 1));
  }
  if (r.max_force > 0.1)
    r.notes.push_back("the largest force is " + std::to_string(r.max_force).substr(0, 6) +
                      " kcal/mol/Å: not at a minimum, so the low modes are unreliable — minimise tightly first (ftol ≤ 0.01)");
  if (r.imaginary > 0) r.notes.push_back(std::to_string(r.imaginary) + " imaginary mode(s) (shown negative): a saddle point or an unconverged minimum");
  if (soft > 0) r.notes.push_back(std::to_string(soft) + " mode(s) below 10 cm⁻¹ left out of the thermodynamics");
  if (m < n) r.notes.push_back(std::to_string(n - m) + " atoms held: their coordinates are fixed, no rigid motions projected");
  return r;
}

std::vector<std::vector<Vec3>> mode_frames(const System& s, const NormalModesResult& r, size_t k, double amplitude, int frames) {
  if (k >= r.mode.size()) throw std::out_of_range("no mode " + std::to_string(k + 1));
  const auto& v = r.mode[k];
  double big = 0;
  for (size_t i = 0; i < s.atoms.size(); ++i) big = std::max(big, std::sqrt(v[3 * i] * v[3 * i] + v[3 * i + 1] * v[3 * i + 1] + v[3 * i + 2] * v[3 * i + 2]));
  const double scale = big > 0 ? amplitude / big : 0;
  std::vector<std::vector<Vec3>> out;
  for (int fr = 0; fr < std::max(2, frames); ++fr) {
    const double a = scale * std::sin(2 * M_PI * fr / std::max(2, frames));
    std::vector<Vec3> p(s.atoms.size());
    for (size_t i = 0; i < s.atoms.size(); ++i)
      for (int q = 0; q < 3; ++q) p[i][q] = s.atoms[i].pos[q] + a * v[3 * i + q];
    out.push_back(std::move(p));
  }
  return out;
}

}  // namespace caps
