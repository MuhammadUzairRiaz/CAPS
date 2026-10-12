// Internal: a small persistent thread pool for static-partition loops. Results are deterministic for a given
// thread count because each worker owns a fixed slice and reductions run in worker order.
#pragma once
#include <algorithm>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <functional>
#include <mutex>
#include <thread>
#include <vector>

namespace caps {

class ThreadPool {
 public:
  explicit ThreadPool(int n) : n_(std::max(1, n)) {
    for (int t = 1; t < n_; ++t) workers_.emplace_back([this, t] { loop(t); });
  }
  ~ThreadPool() {
    {
      std::lock_guard<std::mutex> l(m_);
      stop_ = true;
      gen_.fetch_add(1, std::memory_order_release);
    }
    cv_.notify_all();
    for (auto& w : workers_) w.join();
  }
  int size() const { return n_; }

  // Runs f(worker, begin, end) over [0, count) split into size() contiguous slices; the caller runs slice 0.
  // Every worker is called once per run, so per-worker buffers can be reset inside f.
  // An MD step makes a dozen runs a few hundred µs apart: workers spin a short while for the next one before they
  // sleep, and the caller spins for the last worker, so a run costs no system call while the work keeps coming.
  void run(size_t count, const std::function<void(int, size_t, size_t)>& f) {
    if (n_ == 1) { f(0, 0, count); return; }   // otherwise every worker is called, possibly with an empty range
    job_ = &f;
    count_ = count;
    pending_.store(n_ - 1, std::memory_order_relaxed);
    bool wake;
    {
      std::lock_guard<std::mutex> l(m_);
      gen_.fetch_add(1, std::memory_order_release);
      wake = sleepers_ > 0;
    }
    if (wake) cv_.notify_all();
    slice(0);
    spin([this] { return pending_.load(std::memory_order_acquire) == 0; });
    if (pending_.load(std::memory_order_acquire) != 0) {
      std::unique_lock<std::mutex> l(m_);
      done_.wait(l, [this] { return pending_.load(std::memory_order_acquire) == 0; });
    }
    job_ = nullptr;
  }

  // Runs f(slice, begin, end) over [0, count) cut into `slices` fixed slices, each taken by whichever worker is free (a
  // fast core takes more): what a slice writes, and any sum kept per slice, is the same whichever worker ran it.
  void run_slices(size_t count, size_t slices, const std::function<void(size_t, size_t, size_t)>& f) {
    slices = std::max<size_t>(1, slices);
    if (n_ == 1) {
      for (size_t s = 0; s < slices; ++s) f(s, count * s / slices, count * (s + 1) / slices);
      return;
    }
    std::atomic<size_t> next{0};
    run(size_t(n_), [&](int, size_t, size_t) {
      for (size_t s; (s = next.fetch_add(1, std::memory_order_relaxed)) < slices;) f(s, count * s / slices, count * (s + 1) / slices);
    });
  }
  size_t slices() const { return n_ == 1 ? 1 : 4 * size_t(n_); }   // the usual number: 4 per worker

 private:
  // spins until done() or about 100 µs have passed
  template <class Done>
  static void spin(Done done) {
    const auto t0 = std::chrono::steady_clock::now();
    for (int k = 1; !done(); ++k) {
#if defined(__x86_64__) || defined(_M_X64) || defined(__i386__)
      __builtin_ia32_pause();
#elif defined(__aarch64__) || defined(__arm__)
      __asm__ __volatile__("yield");
#endif
      if (k % 64 == 0 && std::chrono::steady_clock::now() - t0 > std::chrono::microseconds(100)) return;
    }
  }
  void slice(int t) {
    const size_t b = count_ * t / n_, e = count_ * (t + 1) / n_;
    (*job_)(t, b, e);
  }
  void loop(int t) {
    size_t seen = 0;
    for (;;) {
      size_t g = seen;
      spin([&] { return (g = gen_.load(std::memory_order_acquire)) != seen; });
      if (g == seen) {
        std::unique_lock<std::mutex> l(m_);
        ++sleepers_;
        cv_.wait(l, [&] { return gen_.load(std::memory_order_acquire) != seen; });
        --sleepers_;
        g = gen_.load(std::memory_order_acquire);
      }
      seen = g;
      if (stop_) return;
      slice(t);
      if (pending_.fetch_sub(1, std::memory_order_acq_rel) == 1) {
        { std::lock_guard<std::mutex> l(m_); }   // the caller is spinning, or waits under the lock: no lost wake-up
        done_.notify_one();
      }
    }
  }
  int n_;
  std::vector<std::thread> workers_;
  std::mutex m_;
  std::condition_variable cv_, done_;
  const std::function<void(int, size_t, size_t)>* job_ = nullptr;
  size_t count_ = 0;
  std::atomic<size_t> gen_{0};
  std::atomic<int> pending_{0};
  int sleepers_ = 0;
  bool stop_ = false;
};

// Default worker count: set_max_threads, else a cluster job's CAPS_THREADS or SLURM_CPUS_PER_TASK, else performance cores
// where the OS says, otherwise hardware threads, at most 16.
int max_threads();   // caps/config.hpp: 0 = automatic

int env_threads();   // caps/live.hpp: CAPS_THREADS, then SLURM_CPUS_PER_TASK (0: neither)

inline int default_threads() {
  if (const int n = max_threads(); n > 0) return n;
  if (const int n = env_threads(); n > 0) return n;   // a cluster job's cores (not capped: a node has many)
  const unsigned h = std::thread::hardware_concurrency();
  return static_cast<int>(std::clamp(h == 0 ? 1u : h, 1u, 16u));
}

}  // namespace caps
