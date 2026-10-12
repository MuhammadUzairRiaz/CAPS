// The cluster kernel with AVX2 and FMA (x86-64): eight lanes hold two i atoms against the same four j atoms (the 2×N
// layout of GROMACS), so the list is the one the four-lane kernel uses. Built with -mavx2 -mfma and called only on a
// processor that has both (kernel_avx2 says whether it ran).
#include "nbnxm.hpp"
#include "nbnxm_fit.hpp"

#if defined(__AVX2__) && defined(__FMA__)
#include <immintrin.h>
#include <type_traits>
#include <vector>

namespace caps::nbnxm {

namespace {

inline __m256 s8(float a) { return _mm256_set1_ps(a); }
inline __m256 dup4(const float* p) { return _mm256_broadcast_ps(reinterpret_cast<const __m128*>(p)); }   // p[0..3] twice
inline __m256 pair(float a, float b) { return _mm256_set_m128(_mm_set1_ps(b), _mm_set1_ps(a)); }        // a ×4 | b ×4
inline __m256 mask8(unsigned bits) {
  const __m256i lanes = _mm256_setr_epi32(1, 2, 4, 8, 16, 32, 64, 128);
  return _mm256_castsi256_ps(_mm256_cmpeq_epi32(_mm256_and_si256(_mm256_set1_epi32(int(bits)), lanes), lanes));
}
inline __m128 fold(__m256 a) { return _mm_add_ps(_mm256_castps256_ps128(a), _mm256_extractf128_ps(a, 1)); }   // the two halves summed
inline float hsum4(__m128 a) {
  __m128 s = _mm_add_ps(a, _mm_movehl_ps(a, a));
  return _mm_cvtss_f32(_mm_add_ss(s, _mm_shuffle_ps(s, s, 1)));
}
inline double hsum8(__m256 a) { return double(hsum4(fold(a))); }
template <size_t N>
inline __m256 horner(const float (&c)[N], __m256 t) {
  __m256 p = s8(c[N - 1]);
  for (size_t k = N - 1; k-- > 0;) p = _mm256_fmadd_ps(p, t, s8(c[k]));
  return p;
}
inline __m256 rsqrt8(__m256 a) { return _mm256_div_ps(s8(1.0f), _mm256_sqrt_ps(a)); }
// e^x, x ≤ 0 (as simd4::exp_neg)
inline __m256 exp_neg8(__m256 x) {
  const __m256 lo = s8(-87.0f);
  const __m256 ok = _mm256_cmp_ps(lo, x, _CMP_LT_OQ);
  x = _mm256_blendv_ps(lo, x, ok);
  const __m256 n = _mm256_round_ps(_mm256_mul_ps(x, s8(1.44269504f)), _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);
  __m256 r = _mm256_fnmadd_ps(n, s8(0.693359375f), x);
  r = _mm256_fnmadd_ps(n, s8(-2.12194440e-4f), r);
  __m256 p = s8(1.0f / 720);
  p = _mm256_fmadd_ps(p, r, s8(1.0f / 120));
  p = _mm256_fmadd_ps(p, r, s8(1.0f / 24));
  p = _mm256_fmadd_ps(p, r, s8(1.0f / 6));
  p = _mm256_fmadd_ps(p, r, s8(0.5f));
  p = _mm256_fmadd_ps(p, r, s8(1.0f));
  p = _mm256_fmadd_ps(p, r, s8(1.0f));
  const __m256 sc = _mm256_castsi256_ps(_mm256_slli_epi32(_mm256_add_epi32(_mm256_cvtps_epi32(n), _mm256_set1_epi32(127)), 23));
  return _mm256_and_ps(ok, _mm256_mul_ps(p, sc));
}

template <int kLj, bool kComb, int kCoul, bool kE>
void kernel8(const List& L, const Params& p, int tw, float* fxo, float* fyo, float* fzo, Sums& out) {
  const __m256 rc2 = s8(p.rc2), one = s8(1.0f), rcn = s8(p.rcn), rc6 = s8(p.rc6);
  const bool shifted = p.shift != 0;
  const __m256 kap2 = s8(-p.kappa * p.kappa), epk = s8(kErfcP * p.kappa), kc = s8(p.kc), e0 = s8(p.e0), f0 = s8(p.f0), rcv = s8(p.rc);
  const __m256 kap = s8(p.kappa), kap3 = s8(p.kappa * p.kappa * p.kappa), ku = s8(p.kappa * p.kappa / kFitUMax);
  const __m256 twelve = s8(kLj == 2 ? 9.0f : 12.0f), six = s8(6.0f);
  const float* X = L.X.data();
  const float* Y = L.Y.data();
  const float* Z = L.Z.data();
  const float* Q = L.Q.data();
  const float* P12 = L.P12.data();
  const float* P6 = L.P6.data();
  const uint32_t* T = L.T.data();
  const size_t ncodes = L.qint.size() / 3;
  std::vector<__m128> fsh(3 * ncodes, _mm_setzero_ps());
  size_t zero_code = ncodes;
  for (size_t c = 0; c < ncodes; ++c)
    if (L.qint[3 * c] == 0 && L.qint[3 * c + 1] == 0 && L.qint[3 * c + 2] == 0) zero_code = c;
  const size_t cb = L.tb[size_t(tw)], ce = L.te[size_t(tw)];
  const Entry* ents = L.pe[size_t(tw)].data();
  const uint32_t* pfi = L.pf[size_t(tw)].data();
  const bool ct = !L.CT12.empty();
  const size_t nsl = L.perm.size();
  for (size_t ci = cb; ci < ce; ++ci) {
    const size_t si = 4 * ci;
    const double* oi = &L.org[3 * ci];
    // h = 0: i atoms 0 | 1, h = 1: 2 | 3
    __m256 xi[2], yi[2], zi[2], qi[2], p12i[2], p6i[2], fxa[2], fya[2], fza[2];
    for (int h = 0; h < 2; ++h) {
      xi[h] = pair(X[si + 2 * h], X[si + 2 * h + 1]);
      yi[h] = pair(Y[si + 2 * h], Y[si + 2 * h + 1]);
      zi[h] = pair(Z[si + 2 * h], Z[si + 2 * h + 1]);
      qi[h] = pair(Q[si + 2 * h] * p.kcoul, Q[si + 2 * h + 1] * p.kcoul);
      p12i[h] = pair(P12[si + 2 * h], P12[si + 2 * h + 1]);
      p6i[h] = pair(P6[si + 2 * h], P6[si + 2 * h + 1]);
      fxa[h] = fya[h] = fza[h] = _mm256_setzero_ps();
    }
    __m256 ev = _mm256_setzero_ps(), ec = _mm256_setzero_ps();
    for (size_t e = pfi[ci - cb]; e < pfi[ci - cb + 1]; ++e) {
      const Entry& en = ents[e];
      const size_t sj = 4 * size_t(en.cj);
      const double* sq = &L.qvec[3 * size_t(en.code)];
      const double* oj = &L.org[3 * size_t(en.cj)];
      const __m256 sx = s8(float(sq[0] + oj[0] - oi[0])), sy = s8(float(sq[1] + oj[1] - oi[1])), sz = s8(float(sq[2] + oj[2] - oi[2]));
      const __m256 xj = dup4(X + sj), yj = dup4(Y + sj), zj = dup4(Z + sj);
      __m256 qj = _mm256_setzero_ps(), p12j = qj, p6j = qj;
      if constexpr (kCoul != 0) qj = dup4(Q + sj);
      if constexpr (kComb) p12j = dup4(P12 + sj), p6j = dup4(P6 + sj);
      __m256 fjx = _mm256_setzero_ps(), fjy = fjx, fjz = fjx;
      for (int h = 0; h < 2; ++h) {
        const unsigned bits = (en.mask >> (8 * h)) & 0xffu;
        if (!bits) continue;
        const __m256 dx = _mm256_sub_ps(xj, _mm256_sub_ps(xi[h], sx));
        const __m256 dy = _mm256_sub_ps(yj, _mm256_sub_ps(yi[h], sy));
        const __m256 dz = _mm256_sub_ps(zj, _mm256_sub_ps(zi[h], sz));
        __m256 r2 = _mm256_fmadd_ps(dx, dx, _mm256_fmadd_ps(dy, dy, _mm256_mul_ps(dz, dz)));
        const __m256 m = _mm256_and_ps(mask8(bits), _mm256_cmp_ps(r2, rc2, _CMP_LT_OQ));
        r2 = _mm256_blendv_ps(one, r2, m);
        const __m256 rinv = rsqrt8(r2), r2i = _mm256_mul_ps(rinv, rinv), r6i = _mm256_mul_ps(_mm256_mul_ps(r2i, r2i), r2i);
        __m256 c12, c6;
        if constexpr (kComb) {
          c12 = _mm256_mul_ps(p12i[h], p12j);
          c6 = _mm256_mul_ps(p6i[h], p6j);
        } else if (ct) {
          const float* a12 = L.CT12.data() + size_t(T[si + 2 * h]) * nsl + sj;
          const float* b12 = L.CT12.data() + size_t(T[si + 2 * h + 1]) * nsl + sj;
          const float* a6 = L.CT6.data() + size_t(T[si + 2 * h]) * nsl + sj;
          const float* b6 = L.CT6.data() + size_t(T[si + 2 * h + 1]) * nsl + sj;
          c12 = _mm256_set_m128(_mm_loadu_ps(b12), _mm_loadu_ps(a12));
          c6 = _mm256_set_m128(_mm_loadu_ps(b6), _mm_loadu_ps(a6));
        } else {
          const int ra = int(T[si + 2 * h] * p.nt), rb = int(T[si + 2 * h + 1] * p.nt);
          const __m256i idx = _mm256_setr_epi32(ra + int(T[sj]), ra + int(T[sj + 1]), ra + int(T[sj + 2]), ra + int(T[sj + 3]), rb + int(T[sj]),
                                                rb + int(T[sj + 1]), rb + int(T[sj + 2]), rb + int(T[sj + 3]));
          c12 = _mm256_i32gather_ps(p.c12, idx, 4);
          c6 = _mm256_i32gather_ps(p.c6, idx, 4);
        }
        __m256 rn = _mm256_mul_ps(r6i, r6i);
        if constexpr (kLj == 2) rn = _mm256_mul_ps(_mm256_mul_ps(r6i, r2i), rinv);
        const __m256 Aa = _mm256_mul_ps(c12, rn), Bb = _mm256_mul_ps(c6, r6i);
        __m256 fr = _mm256_mul_ps(_mm256_fnmadd_ps(six, Bb, _mm256_mul_ps(twelve, Aa)), r2i);
        if constexpr (kE) {
          __m256 vlj = _mm256_sub_ps(Aa, Bb);
          if (shifted) vlj = _mm256_sub_ps(vlj, _mm256_fnmadd_ps(c6, rc6, _mm256_mul_ps(c12, rcn)));
          ev = _mm256_add_ps(ev, _mm256_and_ps(m, vlj));
        }
        if constexpr (kCoul != 0) {
          const __m256 qq = _mm256_mul_ps(qi[h], qj);
          if constexpr (kCoul == 1 || kCoul == 2) {
            const __m256 t = _mm256_mul_ps(ku, r2);
            const __m256 Fc = _mm256_fmadd_ps(kap3, _mm256_div_ps(horner(kGfP, t), horner(kGfQ, t)), _mm256_mul_ps(rinv, r2i));
            __m256 E = _mm256_setzero_ps();
            if constexpr (kE) E = _mm256_fnmadd_ps(kap, _mm256_div_ps(horner(kGvP, t), horner(kGvQ, t)), rinv);
            if constexpr (kCoul == 1) {
              fr = _mm256_fmadd_ps(qq, Fc, fr);
              if constexpr (kE) ec = _mm256_add_ps(ec, _mm256_and_ps(m, _mm256_mul_ps(qq, E)));
            } else {
              fr = _mm256_fmadd_ps(qq, _mm256_fnmadd_ps(f0, rinv, Fc), fr);
              if constexpr (kE)
                ec = _mm256_add_ps(ec, _mm256_and_ps(m, _mm256_mul_ps(qq, _mm256_fmadd_ps(f0, _mm256_sub_ps(_mm256_mul_ps(r2, rinv), rcv), _mm256_sub_ps(E, e0)))));
            }
          } else {
            const __m256 r = _mm256_mul_ps(r2, rinv);
            const __m256 ex = exp_neg8(_mm256_mul_ps(kap2, r2));
            const __m256 t = _mm256_div_ps(one, _mm256_fmadd_ps(epk, r, one));
            __m256 poly = s8(kErfcA[4]);
            for (int k = 3; k >= 0; --k) poly = _mm256_fmadd_ps(poly, t, s8(kErfcA[k]));
            const __m256 E = _mm256_mul_ps(_mm256_mul_ps(_mm256_mul_ps(poly, t), ex), rinv);
            const __m256 Fc = _mm256_mul_ps(_mm256_fmadd_ps(kc, ex, E), r2i);
            if constexpr (kCoul == 3) {
              fr = _mm256_fmadd_ps(qq, Fc, fr);
              if constexpr (kE) ec = _mm256_add_ps(ec, _mm256_and_ps(m, _mm256_mul_ps(qq, E)));
            } else {
              fr = _mm256_fmadd_ps(qq, _mm256_fnmadd_ps(f0, rinv, Fc), fr);
              if constexpr (kE) ec = _mm256_add_ps(ec, _mm256_and_ps(m, _mm256_mul_ps(qq, _mm256_fmadd_ps(f0, _mm256_sub_ps(r, rcv), _mm256_sub_ps(E, e0)))));
            }
          }
        }
        fr = _mm256_and_ps(m, fr);
        const __m256 fx = _mm256_mul_ps(fr, dx), fy = _mm256_mul_ps(fr, dy), fz = _mm256_mul_ps(fr, dz);
        fjx = _mm256_add_ps(fjx, fx); fjy = _mm256_add_ps(fjy, fy); fjz = _mm256_add_ps(fjz, fz);
        fxa[h] = _mm256_sub_ps(fxa[h], fx); fya[h] = _mm256_sub_ps(fya[h], fy); fza[h] = _mm256_sub_ps(fza[h], fz);
      }
      const __m128 gx = fold(fjx), gy = fold(fjy), gz = fold(fjz);
      _mm_storeu_ps(fxo + sj, _mm_add_ps(_mm_loadu_ps(fxo + sj), gx));
      _mm_storeu_ps(fyo + sj, _mm_add_ps(_mm_loadu_ps(fyo + sj), gy));
      _mm_storeu_ps(fzo + sj, _mm_add_ps(_mm_loadu_ps(fzo + sj), gz));
      if (en.code != zero_code) {
        __m128* fs = &fsh[3 * size_t(en.code)];
        fs[0] = _mm_add_ps(fs[0], gx); fs[1] = _mm_add_ps(fs[1], gy); fs[2] = _mm_add_ps(fs[2], gz);
      }
    }
    for (int h = 0; h < 2; ++h) {
      // lanes 0–3 belong to i atom 2h, lanes 4–7 to 2h + 1
      fxo[si + 2 * h] += hsum4(_mm256_castps256_ps128(fxa[h])); fxo[si + 2 * h + 1] += hsum4(_mm256_extractf128_ps(fxa[h], 1));
      fyo[si + 2 * h] += hsum4(_mm256_castps256_ps128(fya[h])); fyo[si + 2 * h + 1] += hsum4(_mm256_extractf128_ps(fya[h], 1));
      fzo[si + 2 * h] += hsum4(_mm256_castps256_ps128(fza[h])); fzo[si + 2 * h + 1] += hsum4(_mm256_extractf128_ps(fza[h], 1));
    }
    out.evdw += hsum8(ev);
    out.ecoul += hsum8(ec);
  }
  if (out.fsh.size() < 3 * ncodes) out.fsh.resize(3 * ncodes, 0.0);
  for (size_t c = 0; c < ncodes; ++c) {
    if (c == zero_code) continue;
    out.fsh[3 * c] += hsum4(fsh[3 * c]);
    out.fsh[3 * c + 1] += hsum4(fsh[3 * c + 1]);
    out.fsh[3 * c + 2] += hsum4(fsh[3 * c + 2]);
  }
}

template <int kLj, bool kComb, bool kE>
void by_coul8(const List& L, const Params& p, int tw, float* fx, float* fy, float* fz, Sums& out) {
  const bool fit = p.kappa * p.kappa * p.rc2 <= kFitUMax;
  if (p.coul == 1) fit ? kernel8<kLj, kComb, 1, kE>(L, p, tw, fx, fy, fz, out) : kernel8<kLj, kComb, 3, kE>(L, p, tw, fx, fy, fz, out);
  else if (p.coul == 2) fit ? kernel8<kLj, kComb, 2, kE>(L, p, tw, fx, fy, fz, out) : kernel8<kLj, kComb, 4, kE>(L, p, tw, fx, fy, fz, out);
  else kernel8<kLj, kComb, 0, kE>(L, p, tw, fx, fy, fz, out);
}

template <int kLj, bool kComb>
void by_energy8(const List& L, const Params& p, int tw, float* fx, float* fy, float* fz, Sums& out) {
  if (p.energies) by_coul8<kLj, kComb, true>(L, p, tw, fx, fy, fz, out);
  else by_coul8<kLj, kComb, false>(L, p, tw, fx, fy, fz, out);
}

}  // namespace

bool kernel_avx2(const List& L, const Params& p, int t, float* fx, float* fy, float* fz, Sums& out) {
  if (p.lj == 2) {
    if (p.comb) by_energy8<2, true>(L, p, t, fx, fy, fz, out);
    else by_energy8<2, false>(L, p, t, fx, fy, fz, out);
  } else {
    if (p.comb) by_energy8<1, true>(L, p, t, fx, fy, fz, out);
    else by_energy8<1, false>(L, p, t, fx, fy, fz, out);
  }
  return true;
}

bool avx2_built() { return true; }

}  // namespace caps::nbnxm

#else

namespace caps::nbnxm {
bool kernel_avx2(const List&, const Params&, int, float*, float*, float*, Sums&) { return false; }
bool avx2_built() { return false; }
}  // namespace caps::nbnxm

#endif
