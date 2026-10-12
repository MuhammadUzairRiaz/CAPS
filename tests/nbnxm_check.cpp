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
  // every other excluded pair a 1-4 pair keeping half its Coulomb (for the kernel's corrections)
  const double f14 = 0.5;
  std::vector<std::vector<uint32_t>> p14(n);
  for (size_t i = 0; i + 1 < n; i += 4) p14[i].push_back(uint32_t(i + 1)), p14[i + 1].push_back(uint32_t(i));

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

    // direct sum: every pair and image within rc, excluded pairs at their nearest image left out; their Ewald / DSF
    // corrections (Coulomb less the part of the bare 1/r the pair drops) kept apart
    double evdw = 0, ecoul = 0, w[6] = {0, 0, 0, 0, 0, 0}, ecx = 0, wx[6] = {0, 0, 0, 0, 0, 0};
    std::vector<double> fref(3 * n, 0.0), fx(3 * n, 0.0);
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
                if (r2 <= dot(dm, dm) + 1e-6) {
                  if (cs.coul) {
                    const double r = std::sqrt(r2), qq = kc * q[i] * q[j], er = std::erfc(cs.kappa * r), keep = p14[i].empty() ? 0.0 : f14;
                    double E = er / r - (1 - keep) / r, F = (er / r + sp * cs.kappa * std::exp(-cs.kappa * cs.kappa * r2)) / r2 - (1 - keep) / (r2 * r);
                    if (cs.coul == 2) E += -dsf_e0 + dsf_f0 * (r - rc), F -= dsf_f0 / r;
                    ecx += qq * E;
                    const double fr = qq * F;
                    for (int k = 0; k < 3; ++k) fx[3 * j + k] += fr * d[k], fx[3 * i + k] -= fr * d[k];
                    wx[0] += fr * d[0] * d[0]; wx[1] += fr * d[1] * d[1]; wx[2] += fr * d[2] * d[2];
                    wx[3] += fr * d[0] * d[1]; wx[4] += fr * d[0] * d[2]; wx[5] += fr * d[1] * d[2];
                  }
                  continue;
                }
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
      // the coefficients three ways: per-type factors, a per (type, j cluster) table, gathered from the type-pair table
      // and the excluded pairs' corrections left to the caller, or taken by the kernel
      for (const int mode : {0, 1, 2})
      for (const bool exk : {false, true}) {
        const bool comb = mode == 0;
        nbnxm::List L;
        const std::vector<double> none;
        L.build(x, cell, rc + skin, ex, type, nt, {}, q, comb ? f12 : none, comb ? f6 : none, pool, mode == 1 ? c12.data() : nullptr,
                mode == 1 ? c6.data() : nullptr, exk ? &p14 : nullptr);
        L.positions(x, pool);
        L.prune(rc + 0.8, pool);
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
        p.f14 = float(f14);
        const size_t ns = L.perm.size();
        std::vector<double> f(3 * n, 0.0);
        nbnxm::Sums sum;
        std::vector<double> P(3 * ns, 0.0);   // pair forces by slot, for the virial
        for (size_t t = 0; t < L.tb.size(); ++t) {
          std::vector<float> fb(3 * ns, 0.0f);
          nbnxm::kernel(L, p, int(t), fb.data(), fb.data() + ns, fb.data() + 2 * ns, sum);
          for (size_t sl = 0; sl < ns; ++sl)
            for (int k = 0; k < 3; ++k) P[3 * sl + size_t(k)] += fb[size_t(k) * ns + sl];
          for (size_t i = 0; i < n; ++i)
            for (int k = 0; k < 3; ++k) f[3 * i + k] += fb[size_t(k) * ns + L.iperm[i]];
        }
        double W[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
        nbnxm::pair_virial(L, P.data(), 0, ns, W);
        nbnxm::shift_virial(L, sum.fsh, W);
        const double wv[6] = {W[0], W[4], W[8], 0.5 * (W[1] + W[3]), 0.5 * (W[2] + W[6]), 0.5 * (W[5] + W[7])};
        double rms = 0, err = 0;
        const double wexp[6] = {w[0] + (exk ? wx[0] : 0), w[1] + (exk ? wx[1] : 0), w[2] + (exk ? wx[2] : 0),
                                w[3] + (exk ? wx[3] : 0), w[4] + (exk ? wx[4] : 0), w[5] + (exk ? wx[5] : 0)};
        for (size_t k = 0; k < f.size(); ++k) {
          const double fe = fref[k] + (exk ? fx[k] : 0.0);
          rms += fe * fe, err += (f[k] - fe) * (f[k] - fe);
        }
        const double ecexp = ecoul + (exk ? ecx : 0.0);
        const double ev = std::fabs(sum.evdw - evdw) / std::fabs(evdw), ef = std::sqrt(err / rms);
        // relative to the Coulomb energy without the corrections (they are negative, so the corrected total is smaller)
        const double ec = cs.coul ? std::fabs(sum.ecoul - ecexp) / std::fabs(ecoul) : std::fabs(sum.ecoul);
        double ew = 0;
        for (int c = 0; c < 6; ++c) ew = std::max(ew, std::fabs(wv[c] - wexp[c]) / (std::fabs(wexp[0] + wexp[1] + wexp[2]) / 3));
        // erfc by A&S (beyond the fits) is good to 1.5e-7 absolute only
        const double tol_e = cs.kappa * rc > 4 ? 2e-5 : 2e-6;
        const bool ok = ev < 2e-6 && ec < tol_e && ef < 2e-6 && ew < 2e-5;
        fails += !ok;
        std::printf("%-17s %d threads, %-7s %s: vdW %.1e  Coulomb %.1e  forces %.1e  virial %.1e  %s\n", cs.name, threads,
                    mode == 0 ? "factors" : mode == 1 ? "table" : "gather", exk ? "excl" : "    ", ev, ec, ef, ew, ok ? "ok" : "FAIL");
      }
    }
  }
  return fails;
}
