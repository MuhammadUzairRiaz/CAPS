// Four single-precision lanes for the cluster pair kernel: NEON on 64-bit Arm, SSE2 on x86-64 (every x86-64 processor
// has it, so installers stay portable), plain arrays elsewhere. Only what the kernel needs.
#pragma once

#include <cmath>
#include <cstdint>
#include <cstring>

#if defined(CAPS_SIMD4_SCALAR)   // the plain version everywhere (to check the others against)
#elif defined(__aarch64__) || defined(_M_ARM64)
#include <arm_neon.h>
#define CAPS_SIMD4_NEON 1
#elif defined(__SSE2__) || defined(_M_X64)
#include <emmintrin.h>
#define CAPS_SIMD4_SSE 1
#endif

namespace caps::simd4 {

#if defined(CAPS_SIMD4_NEON)

struct F4 { float32x4_t v; };
struct M4 { uint32x4_t v; };   // lane masks: all ones or all zeros

inline F4 load(const float* p) { return {vld1q_f32(p)}; }
inline void store(float* p, F4 a) { vst1q_f32(p, a.v); }
inline F4 set1(float a) { return {vdupq_n_f32(a)}; }
inline F4 zero() { return {vdupq_n_f32(0.0f)}; }
inline F4 operator+(F4 a, F4 b) { return {vaddq_f32(a.v, b.v)}; }
inline F4 operator-(F4 a, F4 b) { return {vsubq_f32(a.v, b.v)}; }
inline F4 operator*(F4 a, F4 b) { return {vmulq_f32(a.v, b.v)}; }
inline F4 operator/(F4 a, F4 b) { return {vdivq_f32(a.v, b.v)}; }
inline F4 fma(F4 a, F4 b, F4 c) { return {vfmaq_f32(c.v, a.v, b.v)}; }    // a b + c
inline F4 fnma(F4 a, F4 b, F4 c) { return {vfmsq_f32(c.v, a.v, b.v)}; }   // c − a b
inline M4 operator<(F4 a, F4 b) { return {vcltq_f32(a.v, b.v)}; }
inline M4 operator&(M4 a, M4 b) { return {vandq_u32(a.v, b.v)}; }
inline F4 select(M4 m, F4 a, F4 b) { return {vbslq_f32(m.v, a.v, b.v)}; }   // m ? a : b
inline F4 masked(M4 m, F4 a) { return {vreinterpretq_f32_u32(vandq_u32(m.v, vreinterpretq_u32_f32(a.v)))}; }
inline M4 mask_bits(unsigned bits) {   // lane l set when bit l is (bits < 16): one load from a table
  alignas(16) static const uint32_t table[16][4] = {
      {0, 0, 0, 0}, {~0u, 0, 0, 0}, {0, ~0u, 0, 0}, {~0u, ~0u, 0, 0}, {0, 0, ~0u, 0}, {~0u, 0, ~0u, 0}, {0, ~0u, ~0u, 0}, {~0u, ~0u, ~0u, 0},
      {0, 0, 0, ~0u}, {~0u, 0, 0, ~0u}, {0, ~0u, 0, ~0u}, {~0u, ~0u, 0, ~0u}, {0, 0, ~0u, ~0u}, {~0u, 0, ~0u, ~0u}, {0, ~0u, ~0u, ~0u}, {~0u, ~0u, ~0u, ~0u}};
  return {vld1q_u32(table[bits & 15u])};
}
inline bool any(M4 m) { return vmaxvq_u32(m.v) != 0; }
inline unsigned bits(M4 m) {   // bit l set for lane l
  static const uint32_t lanes[4] = {1, 2, 4, 8};
  return vaddvq_u32(vandq_u32(m.v, vld1q_u32(lanes)));
}
inline float hsum(F4 a) { return vaddvq_f32(a.v); }
// 1/√a, a > 0: the estimate and two Newton steps (full single precision)
inline F4 rsqrt(F4 a) { return {vdivq_f32(vdupq_n_f32(1.0f), vsqrtq_f32(a.v))}; }
// 2ⁿ for integral n held as floats, |n| ≤ 126
inline F4 pow2i(F4 n) { return {vreinterpretq_f32_s32(vshlq_n_s32(vaddq_s32(vcvtq_s32_f32(n.v), vdupq_n_s32(127)), 23))}; }
inline F4 round_near(F4 a) { return {vrndnq_f32(a.v)}; }
template <int L> inline F4 lane(F4 a) { return {vdupq_laneq_f32(a.v, L)}; }   // lane L in every lane
inline F4 gather(const float* base, const uint32_t* idx) {
  float32x4_t r = vdupq_n_f32(0.0f);
  r = vld1q_lane_f32(base + idx[0], r, 0);
  r = vld1q_lane_f32(base + idx[1], r, 1);
  r = vld1q_lane_f32(base + idx[2], r, 2);
  r = vld1q_lane_f32(base + idx[3], r, 3);
  return {r};
}

#elif defined(CAPS_SIMD4_SSE)

struct F4 { __m128 v; };
struct M4 { __m128 v; };

inline F4 load(const float* p) { return {_mm_loadu_ps(p)}; }
inline void store(float* p, F4 a) { _mm_storeu_ps(p, a.v); }
inline F4 set1(float a) { return {_mm_set1_ps(a)}; }
inline F4 zero() { return {_mm_setzero_ps()}; }
inline F4 operator+(F4 a, F4 b) { return {_mm_add_ps(a.v, b.v)}; }
inline F4 operator-(F4 a, F4 b) { return {_mm_sub_ps(a.v, b.v)}; }
inline F4 operator*(F4 a, F4 b) { return {_mm_mul_ps(a.v, b.v)}; }
inline F4 operator/(F4 a, F4 b) { return {_mm_div_ps(a.v, b.v)}; }
inline F4 fma(F4 a, F4 b, F4 c) { return {_mm_add_ps(_mm_mul_ps(a.v, b.v), c.v)}; }
inline F4 fnma(F4 a, F4 b, F4 c) { return {_mm_sub_ps(c.v, _mm_mul_ps(a.v, b.v))}; }
inline M4 operator<(F4 a, F4 b) { return {_mm_cmplt_ps(a.v, b.v)}; }
inline M4 operator&(M4 a, M4 b) { return {_mm_and_ps(a.v, b.v)}; }
inline F4 select(M4 m, F4 a, F4 b) { return {_mm_or_ps(_mm_and_ps(m.v, a.v), _mm_andnot_ps(m.v, b.v))}; }
inline F4 masked(M4 m, F4 a) { return {_mm_and_ps(m.v, a.v)}; }
inline M4 mask_bits(unsigned bits) {
  const __m128i lanes = _mm_set_epi32(8, 4, 2, 1);
  const __m128i b = _mm_and_si128(_mm_set1_epi32(int(bits)), lanes);
  return {_mm_castsi128_ps(_mm_cmpeq_epi32(b, lanes))};
}
inline bool any(M4 m) { return _mm_movemask_ps(m.v) != 0; }
inline unsigned bits(M4 m) { return unsigned(_mm_movemask_ps(m.v)); }
inline float hsum(F4 a) {
  __m128 s = _mm_add_ps(a.v, _mm_movehl_ps(a.v, a.v));
  s = _mm_add_ss(s, _mm_shuffle_ps(s, s, 1));
  return _mm_cvtss_f32(s);
}
inline F4 rsqrt(F4 a) { return {_mm_div_ps(_mm_set1_ps(1.0f), _mm_sqrt_ps(a.v))}; }
inline F4 pow2i(F4 n) { return {_mm_castsi128_ps(_mm_slli_epi32(_mm_add_epi32(_mm_cvtps_epi32(n.v), _mm_set1_epi32(127)), 23))}; }
inline F4 round_near(F4 a) { return {_mm_cvtepi32_ps(_mm_cvtps_epi32(a.v))}; }   // round to nearest (the default mode)
template <int L> inline F4 lane(F4 a) { return {_mm_shuffle_ps(a.v, a.v, _MM_SHUFFLE(L, L, L, L))}; }
inline F4 gather(const float* base, const uint32_t* idx) { return {_mm_set_ps(base[idx[3]], base[idx[2]], base[idx[1]], base[idx[0]])}; }

#else

struct F4 { float v[4]; };
struct M4 { uint32_t v[4]; };

#define CAPS_F4_LANES(expr) F4 r; for (int l = 0; l < 4; ++l) r.v[l] = (expr); return r
inline F4 load(const float* p) { CAPS_F4_LANES(p[l]); }
inline void store(float* p, F4 a) { for (int l = 0; l < 4; ++l) p[l] = a.v[l]; }
inline F4 set1(float a) { CAPS_F4_LANES(a); }
inline F4 zero() { CAPS_F4_LANES(0.0f); }
inline F4 operator+(F4 a, F4 b) { CAPS_F4_LANES(a.v[l] + b.v[l]); }
inline F4 operator-(F4 a, F4 b) { CAPS_F4_LANES(a.v[l] - b.v[l]); }
inline F4 operator*(F4 a, F4 b) { CAPS_F4_LANES(a.v[l] * b.v[l]); }
inline F4 operator/(F4 a, F4 b) { CAPS_F4_LANES(a.v[l] / b.v[l]); }
inline F4 fma(F4 a, F4 b, F4 c) { CAPS_F4_LANES(a.v[l] * b.v[l] + c.v[l]); }
inline F4 fnma(F4 a, F4 b, F4 c) { CAPS_F4_LANES(c.v[l] - a.v[l] * b.v[l]); }
inline M4 operator<(F4 a, F4 b) { M4 r; for (int l = 0; l < 4; ++l) r.v[l] = a.v[l] < b.v[l] ? 0xffffffffu : 0u; return r; }
inline M4 operator&(M4 a, M4 b) { M4 r; for (int l = 0; l < 4; ++l) r.v[l] = a.v[l] & b.v[l]; return r; }
inline F4 select(M4 m, F4 a, F4 b) { CAPS_F4_LANES(m.v[l] ? a.v[l] : b.v[l]); }
inline F4 masked(M4 m, F4 a) { CAPS_F4_LANES(m.v[l] ? a.v[l] : 0.0f); }
inline M4 mask_bits(unsigned bits) { M4 r; for (int l = 0; l < 4; ++l) r.v[l] = (bits >> l) & 1u ? 0xffffffffu : 0u; return r; }
inline bool any(M4 m) { return (m.v[0] | m.v[1] | m.v[2] | m.v[3]) != 0; }
inline unsigned bits(M4 m) { return (m.v[0] & 1u) | (m.v[1] & 2u) | (m.v[2] & 4u) | (m.v[3] & 8u); }
inline float hsum(F4 a) { return (a.v[0] + a.v[1]) + (a.v[2] + a.v[3]); }
inline F4 rsqrt(F4 a) { CAPS_F4_LANES(1.0f / std::sqrt(a.v[l])); }
inline F4 pow2i(F4 n) { CAPS_F4_LANES(std::ldexp(1.0f, int(n.v[l]))); }
inline F4 round_near(F4 a) { CAPS_F4_LANES(std::nearbyint(a.v[l])); }
template <int L> inline F4 lane(F4 a) { CAPS_F4_LANES(a.v[L]); }
inline F4 gather(const float* base, const uint32_t* idx) { CAPS_F4_LANES(base[idx[l]]); }
#undef CAPS_F4_LANES

#endif

// e^x for x ≤ 0 (x below −87 gives 0): x = n ln2 + r, |r| ≤ ln2/2, e^r by its degree-6 Taylor polynomial (relative error
// below 2e-7), times 2ⁿ
inline F4 exp_neg(F4 x) {
  const F4 lo = set1(-87.0f);
  const M4 ok = lo < x;
  x = select(ok, x, lo);
  const F4 n = round_near(x * set1(1.44269504f));
  // Cody–Waite: ln 2 in two parts, so x − n ln2 keeps its digits
  F4 r = fnma(n, set1(0.693359375f), x);
  r = fnma(n, set1(-2.12194440e-4f), r);
  F4 p = set1(1.0f / 720);
  p = fma(p, r, set1(1.0f / 120));
  p = fma(p, r, set1(1.0f / 24));
  p = fma(p, r, set1(1.0f / 6));
  p = fma(p, r, set1(0.5f));
  p = fma(p, r, set1(1.0f));
  p = fma(p, r, set1(1.0f));
  return masked(ok, p * pow2i(n));
}

}  // namespace caps::simd4
