// C++ side of the k-space kernels (see kspace.hpp and fortran/caps_kspace.f90).
#include "caps/kspace.hpp"

#include <array>

#include "kspace_pool.hpp"

#include <cmath>
#include <stdexcept>

extern "C" {
void caps_f_pme(int n, const double* x, const double* q, const double* cell, const int* grid, int p, double beta, double kcoul, double* energy,
                double* f, double* vir);
void caps_f_ewald(int n, const double* x, const double* q, const double* cell, const int* kmax, double beta, double kcoul, double* energy,
                  double* f, double* vir);
void caps_f_fft3d(double* re, double* im, int k1, int k2, int k3, int sign);
int caps_f_good_size(int n);
void caps_f_bspline_moduli(int k, int p, double* bsp);
void caps_f_pme_spread(int n, const double* x, const double* q, const double* cell, const int* grid, int p, std::complex<double>* g, double* th,
                       double* dth, int* idx);
void caps_f_fft_lines(std::complex<double>* g, int k1, int k2, int k3, int axis, int sign, int first, int last);
void caps_f_pme_influence(std::complex<double>* g, const int* grid, const double* cell, const double* b1, const double* b2, const double* b3,
                          double beta, double kcoul, int first, int last, double* energy, double* vir);
void caps_f_pme_gather(int n, const double* q, const double* cell, const int* grid, int p, const std::complex<double>* g, const double* th,
                       const double* dth, const int* idx, int first, int last, double* f);
}

namespace caps {

double ewald_beta(double cutoff, double rtol) {
  // erfc(β rc) decreases with β: bisection on [0, 20/rc]
  double lo = 0, hi = 20.0 / cutoff;
  for (int it = 0; it < 100; ++it) {
    const double mid = 0.5 * (lo + hi);
    (std::erfc(mid * cutoff) > rtol ? lo : hi) = mid;
  }
  return 0.5 * (lo + hi);
}

PmeGrid pme_grid(const Cell& cell, double beta, double spacing, int order) {
  PmeGrid g;
  g.beta = beta;
  g.order = order;
  const Vec3 e[3] = {cell.a, cell.b, cell.c};
  for (int d = 0; d < 3; ++d) g.k[d] = caps_f_good_size(std::max(order + 1, int(std::ceil(norm(e[d]) / spacing))));
  return g;
}

namespace {
void cell_array(const Cell& c, double out[9]) {
  const Vec3 e[3] = {c.a, c.b, c.c};
  for (int r = 0; r < 3; ++r)
    for (int k = 0; k < 3; ++k) out[3 * r + k] = e[r][k];
}
}  // namespace

double pme_reciprocal(const std::vector<double>& x, const std::vector<double>& q, const Cell& cell, const PmeGrid& g, std::vector<double>& f,
                      double vir[6]) {
  const int n = int(q.size());
  if (x.size() != 3 * q.size() || f.size() != x.size()) throw std::invalid_argument("pme: sizes of positions, charges and forces differ");
  double c[9], e = 0;
  cell_array(cell, c);
  caps_f_pme(n, x.data(), q.data(), c, g.k, g.order, g.beta, kCoulombConstant, &e, f.data(), vir);
  return e;
}

double ewald_reciprocal(const std::vector<double>& x, const std::vector<double>& q, const Cell& cell, double beta, const int kmax[3],
                        std::vector<double>& f, double vir[6]) {
  double c[9], e = 0;
  cell_array(cell, c);
  caps_f_ewald(int(q.size()), x.data(), q.data(), c, kmax, beta, kCoulombConstant, &e, f.data(), vir);
  return e;
}

double pme_reciprocal(const std::vector<double>& x, const std::vector<double>& q, const Cell& cell, const PmeGrid& g, std::vector<double>& f,
                      double vir[6], ThreadPool& pool) {
  const int n = int(q.size()), p = g.order, k1 = g.k[0], k2 = g.k[1], k3 = g.k[2];
  if (x.size() != 3 * q.size() || f.size() != x.size()) throw std::invalid_argument("pme: sizes of positions, charges and forces differ");
  double c[9];
  cell_array(cell, c);
  std::vector<std::complex<double>> grid(size_t(k1) * size_t(k2) * size_t(k3));
  std::vector<double> th(size_t(p) * 3 * size_t(n)), dth(th.size());
  std::vector<int> idx(th.size());
  caps_f_pme_spread(n, x.data(), q.data(), c, g.k, p, grid.data(), th.data(), dth.data(), idx.data());
  // lines of each axis split over the workers (disjoint ranges: no two write the same grid point)
  auto transform = [&](int sign) {
    const int lines[3] = {k2 * k3, k1 * k3, k1 * k2};
    for (int axis = 1; axis <= 3; ++axis)
      pool.run(size_t(lines[axis - 1]), [&](int, size_t b, size_t e) {
        if (b < e) caps_f_fft_lines(grid.data(), k1, k2, k3, axis, sign, int(b), int(e));
      });
  };
  transform(+1);
  std::vector<double> b1(static_cast<size_t>(k1)), b2(static_cast<size_t>(k2)), b3(static_cast<size_t>(k3));
  caps_f_bspline_moduli(k1, p, b1.data());
  caps_f_bspline_moduli(k2, p, b2.data());
  caps_f_bspline_moduli(k3, p, b3.data());
  std::vector<std::array<double, 7>> part(static_cast<size_t>(pool.size()));
  for (auto& a : part) a.fill(0.0);
  pool.run(size_t(k3), [&](int t, size_t b, size_t e) {
    if (b < e) caps_f_pme_influence(grid.data(), g.k, c, b1.data(), b2.data(), b3.data(), g.beta, kCoulombConstant, int(b), int(e), &part[size_t(t)][0],
                                     &part[size_t(t)][1]);
  });
  transform(-1);
  pool.run(size_t(n), [&](int, size_t b, size_t e) {
    if (b < e) caps_f_pme_gather(n, q.data(), c, g.k, p, grid.data(), th.data(), dth.data(), idx.data(), int(b), int(e), f.data());
  });
  double energy = 0;
  for (int k = 0; k < 6; ++k) vir[k] = 0;
  for (const auto& a : part) {   // worker order: deterministic for a given thread count
    energy += a[0];
    for (int k = 0; k < 6; ++k) vir[k] += a[1 + k];
  }
  return energy;
}

void fft3d(std::vector<std::complex<double>>& data, int k1, int k2, int k3, int sign) {
  std::vector<double> re(data.size()), im(data.size());
  for (size_t i = 0; i < data.size(); ++i) re[i] = data[i].real(), im[i] = data[i].imag();
  caps_f_fft3d(re.data(), im.data(), k1, k2, k3, sign);
  for (size_t i = 0; i < data.size(); ++i) data[i] = {re[i], im[i]};
}

int fft_good_size(int n) { return caps_f_good_size(n); }

}  // namespace caps
