// Stand-alone check of the cluster pair list and kernels against a direct double-precision sum, for building them on
// another instruction set than the test suite's (x86-64 under emulation, or -DCAPS_SIMD4_SCALAR):
//   c++ -std=c++20 -O2 -Icore/include -Icore/src tests/nbnxm_check.cpp core/src/nbnxm.cpp core/src/system.cpp core/src/elements.cpp \
//       core/src/nbnxm_avx2.cpp          (that one with -mavx2 -mfma for the eight-lane kernel on x86-64; CAPS_SIMD=avx2|sse2 picks)
// Random charged Lennard-Jones particles in a triclinic cell, a few excluded pairs; 12-6 and 9-6, Ewald, DSF (the fits)
// and Ewald beyond the fits (erfc), shifted and not, coefficients as factors and from the table; prints the errors.
#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

#include "nbnxm.hpp"
#include "parallel.hpp"

using namespace caps;

namespace {

struct Case {
  const char* name;
  int lj;        // 1 12-6, 2 9-6
  int coul;      // 0 none, 1 Ewald, 2 DSF
  double kappa;
  bool shift;
};

}  // namespace

int main() {
  const double rc = 8.0, skin = 2.0, kc = 332.06371;
  Cell cell;
  cell.a = {24, 0, 0};
  cell.b = {3, 22, 0};
  cell.c = {-2, 4, 23};
  const size_t n = 1100, nt = 3;
  std::mt19937 rng(7);
  std::uniform_real_distribution<double> U(0, 1);
  std::vector<double> x;
  std::vector<int> type(n);
  std::vector<double> q(n);
  // on a jittered lattice, so no two particles start closer than about 2 Å
  const int m = 11;
  for (size_t i = 0; i < n; ++i) {
    const double f[3] = {(double(i % m) + 0.3 * U(rng)) / m, (double(i / m % m) + 0.3 * U(rng)) / m, (double(i / (m * m)) + 0.3 * U(rng)) / 10.0};
    for (int k = 0; k < 3; ++k) x.push_back(f[0] * cell.a[k] + f[1] * cell.b[k] + f[2] * cell.c[k]);
    type[i] = int(i % nt);
    q[i] = (i % 2 ? 0.35 : -0.35) * (1 + 0.1 * type[i]);
  }
  // excluded: each particle with the next (a chain), sorted
  std::vector<std::vector<uint32_t>> ex(n);
  for (size_t i = 0; i + 1 < n; i += 2) ex[i].push_back(uint32_t(i + 1)), ex[i + 1].push_back(uint32_t(i));
  const double eps[nt] = {0.07, 0.03, 0.11}, sig[nt] = {3.5, 2.5, 3.2};

  int fails = 0;
  std::printf("kernel: %s\n", nbnxm::instruction_set());
  const Case cases[] = {{"12-6 Ewald", 1, 1, 0.32, false},     {"9-6 DSF shifted", 2, 2, 0.2, true},
                        {"12-6 Ewald erfc", 1, 1, 0.6, false}, {"9-6 no charges", 2, 0, 0.0, false}};
  for (const Case& cs : cases) {
    // geometric coefficients: 12-6 4εσ¹², 4εσ⁶; 9-6 2εσ⁹, 3εσ⁶
    std::vector<double> f12(nt), f6(nt);
    std::vector<float> c12(nt * nt), c6(nt * nt);
    for (size_t a = 0; a < nt; ++a) {
      f12[a] = std::sqrt(cs.lj == 2 ? 2 * eps[a] * std::pow(sig[a], 9) : 4 * eps[a] * std::pow(sig[a], 12));
      f6[a] = std::sqrt((cs.lj == 2 ? 3 : 4) * eps[a] * std::pow(sig[a], 6));
    }
    for (size_t a = 0; a < nt; ++a)
      for (size_t b = 0; b < nt; ++b) c12[a * nt + b] = float(f12[a] * f12[b]), c6[a * nt + b] = float(f6[a] * f6[b]);
    const double pw = cs.lj == 2 ? 9 : 12, sp = 2 / std::sqrt(M_PI);
    const double dsf_e0 = std::erfc(cs.kappa * rc) / rc;
    const double dsf_f0 = std::erfc(cs.kappa * rc) / (rc * rc) + sp * cs.kappa * std::exp(-cs.kappa * cs.kappa * rc * rc) / rc;

    // direct sum: every pair and image within rc, excluded pairs at their nearest image left out
    double evdw = 0, ecoul = 0, w[6] = {0, 0, 0, 0, 0, 0};
    std::vector<double> fref(3 * n, 0.0);
    for (size_t i = 0; i < n; ++i)
      for (size_t j = i; j < n; ++j)
        for (int ia = -1; ia <= 1; ++ia)
          for (int ib = -1; ib <= 1; ++ib)
            for (int ic = -1; ic <= 1; ++ic) {
              if (j == i && (ia > 0 || (ia == 0 && (ib > 0 || (ib == 0 && ic >= 0))))) continue;   // self: one of ± images
              Vec3 d;
              for (int k = 0; k < 3; ++k) d[k] = x[3 * j + k] - x[3 * i + k] + ia * cell.a[k] + ib * cell.b[k] + ic * cell.c[k];
              const double r2 = dot(d, d);
              if (r2 >= rc * rc) continue;
              if (j != i && !ex[i].empty() && ex[i][0] == j) {
                const Vec3 dm = cell.minimum_image(Vec3{x[3 * j] - x[3 * i], x[3 * j + 1] - x[3 * i + 1], x[3 * j + 2] - x[3 * i + 2]});
                if (r2 <= dot(dm, dm) + 1e-6) continue;
              }
              const double r = std::sqrt(r2), a12 = f12[size_t(type[i])] * f12[size_t(type[j])], a6 = f6[size_t(type[i])] * f6[size_t(type[j])];
              evdw += a12 * std::pow(r, -pw) - a6 * std::pow(r, -6.0);
              if (cs.shift) evdw -= a12 * std::pow(rc, -pw) - a6 * std::pow(rc, -6.0);
              double fr = (pw * a12 * std::pow(r, -pw) - 6 * a6 * std::pow(r, -6.0)) / r2;
              if (cs.coul) {
                const double qq = kc * q[i] * q[j], er = std::erfc(cs.kappa * r), E = er / r;
                const double F = (er / r + sp * cs.kappa * std::exp(-cs.kappa * cs.kappa * r2)) / r2;
                if (cs.coul == 1) ecoul += qq * E, fr += qq * F;
                else ecoul += qq * (E - dsf_e0 + dsf_f0 * (r - rc)), fr += qq * (F - dsf_f0 / r);
              }
              for (int k = 0; k < 3; ++k) fref[3 * j + k] += fr * d[k], fref[3 * i + k] -= fr * d[k];
              w[0] += fr * d[0] * d[0]; w[1] += fr * d[1] * d[1]; w[2] += fr * d[2] * d[2];
              w[3] += fr * d[0] * d[1]; w[4] += fr * d[0] * d[2]; w[5] += fr * d[1] * d[2];
            }

    for (int threads : {1, 3}) {
      ThreadPool pool(threads);
      nbnxm::List L;
      L.build(x, cell, rc + skin, ex, type, nt, {}, q, f12, f6, pool);
      L.positions(x);
      L.prune(rc + 0.8, pool);
      for (const bool comb : {true, false}) {
        nbnxm::Params p;
        p.rc2 = float(rc * rc);
        p.rc = float(rc);
        p.lj = cs.lj;
        p.comb = comb;
        p.c12 = c12.data();
        p.c6 = c6.data();
        p.nt = nt;
        p.shift = cs.shift ? 1.0f : 0.0f;
        p.rcn = float(std::pow(rc, -pw));
        p.rc6 = float(std::pow(rc, -6.0));
        p.coul = cs.coul;
        p.kappa = float(cs.kappa);
        p.kc = float(sp * cs.kappa);
        p.e0 = float(dsf_e0);
        p.f0 = float(dsf_f0);
        p.kcoul = float(kc);
        if (!comb) L.P12.assign(L.P12.size(), 0.0f), L.P6.assign(L.P6.size(), 0.0f);
        const size_t ns = L.perm.size();
        std::vector<double> f(3 * n, 0.0);
        nbnxm::Sums sum;
        for (int t = 0; t < threads; ++t) {
          std::vector<float> fb(3 * ns, 0.0f);
          nbnxm::kernel(L, p, t, fb.data(), fb.data() + ns, fb.data() + 2 * ns, sum);
          for (size_t i = 0; i < n; ++i)
            for (int k = 0; k < 3; ++k) f[3 * i + k] += fb[size_t(k) * ns + L.iperm[i]];
        }
        double rms = 0, err = 0;
        for (size_t k = 0; k < f.size(); ++k) rms += fref[k] * fref[k], err += (f[k] - fref[k]) * (f[k] - fref[k]);
        const double ev = std::fabs(sum.evdw - evdw) / std::fabs(evdw), ef = std::sqrt(err / rms);
        const double ec = cs.coul ? std::fabs(sum.ecoul - ecoul) / std::fabs(ecoul) : std::fabs(sum.ecoul);
        double ew = 0;
        for (int c = 0; c < 6; ++c) ew = std::max(ew, std::fabs(sum.w[c] - w[c]) / (std::fabs(w[0] + w[1] + w[2]) / 3));
        // erfc by A&S (beyond the fits) is good to 1.5e-7 absolute only
        const double tol_e = cs.kappa * rc > 4 ? 2e-5 : 2e-6;
        const bool ok = ev < 2e-6 && ec < tol_e && ef < 2e-6 && ew < 2e-5;
        fails += !ok;
        std::printf("%-17s %d threads, %-7s: vdW %.1e  Coulomb %.1e  forces %.1e  virial %.1e  %s\n", cs.name, threads, comb ? "factors" : "table", ev,
                    ec, ef, ew, ok ? "ok" : "FAIL");
      }
    }
  }
  return fails;
}
