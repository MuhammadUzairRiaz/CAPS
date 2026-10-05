// CAPS random distributions, the same numbers on every platform. std::mt19937_64 is specified exactly by the standard,
// but std::uniform_real_distribution, normal_distribution, uniform_int_distribution, shuffle, discrete_distribution and
// gamma_distribution are not: libc++ (macOS), libstdc++ (Linux, MinGW) and MSVC turn one seed into different structures.
// These keep the standard interfaces with algorithms fixed here, so a seed in a provenance record reproduces the run on
// any machine.
//   uniform real: the engine's top 53 bits (64-bit engines) scaled to [0, 1), then to [a, b)
//   uniform int:  rejection on the largest multiple of the range below 2⁶⁴ (no modulo bias)
//   normal:       Marsaglia's polar method, the second value kept for the next call
//   shuffle:      Fisher–Yates from the back, swapping with a uniform index in [0, i]
//   discrete:     the cumulative weights, a uniform real searched
//   gamma:        Marsaglia & Tsang (2000) for shape ≥ 1, boosted by U^(1/shape) below 1
#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iterator>
#include <limits>
#include <numeric>
#include <utility>
#include <vector>

namespace caps {

// a double in [0, 1) from the engine's top 53 bits (two 32-bit draws for a 32-bit engine)
template <class G>
inline double uniform01(G& g) {
  if constexpr (sizeof(typename G::result_type) >= 8) {
    return double(uint64_t(g() - G::min()) >> 11) * 0x1.0p-53;
  } else {
    const uint64_t a = uint64_t(g() - G::min()) >> 5, b = uint64_t(g() - G::min()) >> 6;   // 27 + 26 bits
    return double(a * 67108864u + b) * 0x1.0p-53;
  }
}

template <class G>
inline uint64_t random_u64(G& g) {
  if constexpr (sizeof(typename G::result_type) >= 8) return uint64_t(g() - G::min());
  else return (uint64_t(g() - G::min()) << 32) | uint64_t(g() - G::min());
}

template <class T = double>
class UniformReal {
 public:
  explicit UniformReal(T a = 0, T b = 1) : a_(a), b_(b) {}
  template <class G> T operator()(G& g) const { return a_ + (b_ - a_) * T(uniform01(g)); }
  T a() const { return a_; }
  T b() const { return b_; }
 private:
  T a_, b_;
};

template <class T = int>
class UniformInt {
 public:
  explicit UniformInt(T a = 0, T b = std::numeric_limits<T>::max()) : a_(a), b_(b) {}
  template <class G> T operator()(G& g) const {
    const uint64_t range = uint64_t(b_) - uint64_t(a_);   // b − a, two's complement for signed types
    if (range == std::numeric_limits<uint64_t>::max()) return T(random_u64(g));
    const uint64_t n = range + 1, limit = std::numeric_limits<uint64_t>::max() - std::numeric_limits<uint64_t>::max() % n;
    uint64_t x;
    do x = random_u64(g); while (x >= limit);
    return T(uint64_t(a_) + x % n);
  }
 private:
  T a_, b_;
};

template <class T = double>
class Normal {
 public:
  explicit Normal(T mean = 0, T sd = 1) : m_(mean), s_(sd) {}
  template <class G> T operator()(G& g) {
    if (have_) { have_ = false; return m_ + s_ * spare_; }
    double u, v, q;
    do {
      u = 2 * uniform01(g) - 1, v = 2 * uniform01(g) - 1, q = u * u + v * v;
    } while (q >= 1 || q == 0);
    const double f = std::sqrt(-2 * std::log(q) / q);
    spare_ = T(v * f), have_ = true;
    return m_ + s_ * T(u * f);
  }
  T mean() const { return m_; }
  T stddev() const { return s_; }
  void reset() { have_ = false; }
 private:
  T m_, s_, spare_ = 0;
  bool have_ = false;
};

template <class It, class G>
inline void shuffle(It first, It last, G& g) {
  const auto n = std::distance(first, last);
  for (auto i = n - 1; i > 0; --i) {
    const auto j = UniformInt<int64_t>(0, int64_t(i))(g);
    using std::swap;
    swap(*(first + i), *(first + j));
  }
}

template <class T = int>
class Discrete {
 public:
  template <class It> Discrete(It first, It last) {
    double t = 0;
    for (auto it = first; it != last; ++it) t += std::max(0.0, double(*it)), cum_.push_back(t);
  }
  template <class G> T operator()(G& g) const {
    if (cum_.empty() || cum_.back() <= 0) return T(0);
    const double x = uniform01(g) * cum_.back();
    const auto it = std::upper_bound(cum_.begin(), cum_.end(), x);
    return T(std::min<std::ptrdiff_t>(std::distance(cum_.begin(), it), std::ptrdiff_t(cum_.size()) - 1));
  }
 private:
  std::vector<double> cum_;
};

// shape k, scale θ: mean kθ
template <class T = double>
class Gamma {
 public:
  explicit Gamma(T shape = 1, T scale = 1) : k_(shape), th_(scale) {}
  template <class G> T operator()(G& g) {
    if (k_ <= 0) return T(0);
    if (k_ < 1) {   // Γ(k) = Γ(k + 1) U^(1/k)
      Gamma up(k_ + 1, th_);
      double u;
      do u = uniform01(g); while (u == 0);
      return up(g) * T(std::pow(u, 1.0 / double(k_)));
    }
    const double d = double(k_) - 1.0 / 3.0, c = 1.0 / std::sqrt(9 * d);
    for (;;) {
      double x, v;
      do {
        x = n_(g);
        v = 1 + c * x;
      } while (v <= 0);
      v = v * v * v;
      double u;
      do u = uniform01(g); while (u == 0);
      if (u < 1 - 0.0331 * x * x * x * x || std::log(u) < 0.5 * x * x + d * (1 - v + std::log(v))) return T(d * v) * th_;
    }
  }
 private:
  T k_, th_;
  Normal<double> n_;
};

}  // namespace caps
