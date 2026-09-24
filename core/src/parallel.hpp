// Internal: a small persistent thread pool for static-partition loops. Results are deterministic for a given
// thread count because each worker owns a fixed slice and reductions run in worker order.
#pragma once
#include <algorithm>
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
      ++gen_;
    }
    cv_.notify_all();
    for (auto& w : workers_) w.join();
  }
  int size() const { return n_; }

  // Runs f(worker, begin, end) over [0, count) split into size() contiguous slices; the caller runs slice 0.
  // Every worker is called once per run, so per-worker buffers can be reset inside f.
  void run(size_t count, const std::function<void(int, size_t, size_t)>& f) {
    if (n_ == 1) { f(0, 0, count); return; }   // otherwise every worker is called, possibly with an empty range
    {
      std::lock_guard<std::mutex> l(m_);
      job_ = &f;
      count_ = count;
      pending_ = n_ - 1;
      ++gen_;
    }
    cv_.notify_all();
    slice(0);
    std::unique_lock<std::mutex> l(m_);
    done_.wait(l, [this] { return pending_ == 0; });
    job_ = nullptr;
  }

 private:
  void slice(int t) {
    const size_t b = count_ * t / n_, e = count_ * (t + 1) / n_;
    (*job_)(t, b, e);
  }
  void loop(int t) {
    size_t seen = 0;
    for (;;) {
      {
        std::unique_lock<std::mutex> l(m_);
        cv_.wait(l, [&] { return gen_ != seen; });
        seen = gen_;
        if (stop_) return;
      }
      slice(t);
      {
        std::lock_guard<std::mutex> l(m_);
        if (--pending_ == 0) done_.notify_one();
      }
    }
  }
  int n_;
  std::vector<std::thread> workers_;
  std::mutex m_;
  std::condition_variable cv_, done_;
  const std::function<void(int, size_t, size_t)>* job_ = nullptr;
  size_t count_ = 0, gen_ = 0;
  int pending_ = 0;
  bool stop_ = false;
};

// Default worker count: performance cores where the OS says, otherwise hardware threads, at most 16.
int max_threads();   // caps/config.hpp: 0 = automatic

inline int default_threads() {
  if (const int n = max_threads(); n > 0) return n;
  const unsigned h = std::thread::hardware_concurrency();
  return static_cast<int>(std::clamp(h == 0 ? 1u : h, 1u, 16u));
}

}  // namespace caps
