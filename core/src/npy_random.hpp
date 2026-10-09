// NumPy's default random generator, bit for bit: np.random.default_rng(seed) (SeedSequence → PCG64, XSL-RR 128/64)
// and Generator.shuffle (Fisher–Yates with masked rejection, random_interval). CAPS uses it where a structure must be
// identical to one made with the numpy-based workflow for the same seed (MXene mixed terminations), on every platform.
#pragma once
#include <cstdint>
#include <utility>
#include <vector>

namespace caps::npy {

class Pcg64 {
 public:
  explicit Pcg64(uint64_t seed) {
    // SeedSequence(entropy = seed): the entropy as 32-bit words, mixed into a pool of four
    std::vector<uint32_t> ent;
    do { ent.push_back(uint32_t(seed & 0xffffffffu)); seed >>= 32; } while (seed);
    uint32_t pool[4];
    uint32_t hc = 0x43b0d7e5u;
    auto hashmix = [&](uint32_t v) {
      v ^= hc;
      hc *= 0x931e8875u;
      v *= hc;
      v ^= v >> 16;
      return v;
    };
    auto mix = [](uint32_t x, uint32_t y) {
      uint32_t r = 0xca01f9ddu * x - 0x4973f715u * y;
      r ^= r >> 16;
      return r;
    };
    for (size_t i = 0; i < 4; ++i) pool[i] = hashmix(i < ent.size() ? ent[i] : 0u);
    for (size_t s = 0; s < 4; ++s)
      for (size_t d = 0; d < 4; ++d)
        if (s != d) pool[d] = mix(pool[d], hashmix(pool[s]));
    for (size_t s = 4; s < ent.size(); ++s)
      for (size_t d = 0; d < 4; ++d) pool[d] = mix(pool[d], hashmix(ent[s]));
    // generate_state(4, uint64): eight 32-bit words from the pool, paired little-endian
    uint32_t w[8];
    uint32_t hb = 0x8b51f9ddu;
    for (int i = 0; i < 8; ++i) {
      uint32_t v = pool[i % 4];
      v ^= hb;
      hb *= 0x58f38dedu;
      v *= hb;
      v ^= v >> 16;
      w[i] = v;
    }
    uint64_t val[4];
    for (int k = 0; k < 4; ++k) val[k] = uint64_t(w[2 * k]) | (uint64_t(w[2 * k + 1]) << 32);
    // pcg64_set_seed: state seed (val0:val1), increment (val2:val3), pcg_setseq_128_srandom_r
    const U128 initstate{val[0], val[1]}, initseq{val[2], val[3]};
    state_ = {0, 0};
    inc_ = shl1_or1(initseq);
    step();
    state_ = add(state_, initstate);
    step();
  }

  uint64_t next64() {
    step();
    const uint64_t x = state_.hi ^ state_.lo;
    const unsigned rot = unsigned(state_.hi >> 58);
    return (x >> rot) | (x << ((64 - rot) & 63));
  }
  uint32_t next32() {
    if (has32_) { has32_ = false; return buf32_; }
    const uint64_t n = next64();
    has32_ = true;
    buf32_ = uint32_t(n >> 32);
    return uint32_t(n & 0xffffffffu);
  }
  // random_interval(max): uniform in [0, max] by masked rejection
  uint64_t interval(uint64_t max) {
    if (max == 0) return 0;
    uint64_t mask = max;
    mask |= mask >> 1; mask |= mask >> 2; mask |= mask >> 4; mask |= mask >> 8; mask |= mask >> 16; mask |= mask >> 32;
    uint64_t v;
    if (max <= 0xffffffffull) { while ((v = (next32() & mask)) > max) {} }
    else { while ((v = (next64() & mask)) > max) {} }
    return v;
  }
  template <class T> void shuffle(std::vector<T>& x) {
    for (size_t i = x.size(); i-- > 1;) std::swap(x[i], x[size_t(interval(i))]);
  }

 private:
  struct U128 { uint64_t hi, lo; };
  static U128 add(U128 a, U128 b) { const uint64_t lo = a.lo + b.lo; return {a.hi + b.hi + (lo < a.lo ? 1 : 0), lo}; }
  static U128 mul(U128 a, U128 b) {
    // low 128 bits of a × b
    const uint64_t a0 = a.lo & 0xffffffffu, a1 = a.lo >> 32, b0 = b.lo & 0xffffffffu, b1 = b.lo >> 32;
    const uint64_t p00 = a0 * b0, p01 = a0 * b1, p10 = a1 * b0, p11 = a1 * b1;
    const uint64_t mid = (p00 >> 32) + (p01 & 0xffffffffu) + (p10 & 0xffffffffu);
    const uint64_t lo = (p00 & 0xffffffffu) | (mid << 32);
    const uint64_t hi = p11 + (p01 >> 32) + (p10 >> 32) + (mid >> 32) + a.hi * b.lo + a.lo * b.hi;
    return {hi, lo};
  }
  static U128 shl1_or1(U128 a) { return {(a.hi << 1) | (a.lo >> 63), (a.lo << 1) | 1u}; }
  void step() { state_ = add(mul(state_, kMult), inc_); }
  static constexpr U128 kMult{0x2360ED051FC65DA4ull, 0x4385DF649FCCF645ull};
  U128 state_{0, 0}, inc_{0, 0};
  bool has32_ = false;
  uint32_t buf32_ = 0;
};

}  // namespace caps::npy
