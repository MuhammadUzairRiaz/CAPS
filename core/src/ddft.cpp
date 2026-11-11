// CAPS density-field dynamics: Gaussian-chain dynamic mean-field density functional on a grid (see ddft.hpp).
#include "caps/rng.hpp"
#include "caps/ddft.hpp"

#include <algorithm>
#include <cctype>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdio>
#include <random>
#include <stdexcept>

#include "caps/kspace.hpp"

namespace caps {

namespace {
constexpr double kPi = 3.14159265358979323846;
using cd = std::complex<double>;

// "A8B8" or "AAAABBBB" → the letters, one per bead
std::string expand(const std::string& s) {
  std::string out;
  for (size_t i = 0; i < s.size();) {
    const char c = s[i];
    if (!std::isalpha(static_cast<unsigned char>(c))) throw std::invalid_argument("density fields: a sequence is letters with counts, e.g. A8B8");
    size_t j = i + 1;
    int n = 0;
    while (j < s.size() && std::isdigit(static_cast<unsigned char>(s[j]))) n = n * 10 + (s[j] - '0'), ++j;
    out.append(size_t(j > i + 1 ? n : 1), char(std::toupper(static_cast<unsigned char>(c))));
    i = j;
  }
  if (out.empty()) throw std::invalid_argument("density fields: an empty sequence");
  return out;
}

struct Model {
  int nx = 1, ny = 1, nz = 1;
  size_t M = 1;
  std::vector<std::string> types;
  std::vector<std::vector<int>> seq;     // per molecule, the type of each bead
  std::vector<double> frac;              // per molecule, share of all beads
  std::vector<double> mean;              // per type, θ̄
  std::vector<double> kernel;            // exp(−q²a²/6) per grid point
  std::vector<double> q2;                // q² per grid point
  std::vector<double> rinv;              // per grid point, R(q)⁻¹ (nt × nt): R_IJ = Σ_m (f_m/N_m) Σ_{s∈I,s'∈J} K^|s−s'| (0 at q = 0)
  std::vector<cd> buf;

  explicit Model(const DdftOptions& o) {
    nx = std::max(1, o.grid[0]), ny = std::max(1, o.grid[1]), nz = std::max(1, o.grid[2]);
    M = size_t(nx) * size_t(ny) * size_t(nz);
    if (o.molecules.empty()) throw std::invalid_argument("density fields: no molecules");
    double ft = 0;
    for (const auto& m : o.molecules) {
      if (m.fraction < 0) throw std::invalid_argument("density fields: fractions ≥ 0");
      ft += m.fraction;
    }
    if (ft <= 0) throw std::invalid_argument("density fields: the fractions add up to zero");
    for (const auto& m : o.molecules) {
      const std::string e = expand(m.sequence);
      std::vector<int> s;
      for (char c : e) {
        auto it = std::find(types.begin(), types.end(), std::string(1, c));
        if (it == types.end()) types.push_back(std::string(1, c)), it = types.end() - 1;
        s.push_back(int(it - types.begin()));
      }
      seq.push_back(std::move(s));
      frac.push_back(m.fraction / ft);
    }
    const size_t nt = types.size();
    mean.assign(nt, 0.0);
    for (size_t m = 0; m < seq.size(); ++m)
      for (int t : seq[m]) mean[size_t(t)] += frac[m] / double(seq[m].size());
    kernel.resize(M), q2.resize(M);
    auto wave = [](int i, int n) { const int m = i <= n / 2 ? i : i - n; return 2 * kPi * double(m) / double(n); };
    for (int k = 0; k < nz; ++k)
      for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i) {
          const double qx = nx > 1 ? wave(i, nx) : 0, qy = ny > 1 ? wave(j, ny) : 0, qz = nz > 1 ? wave(k, nz) : 0;
          const size_t c = size_t(i) + size_t(nx) * (size_t(j) + size_t(ny) * size_t(k));
          q2[c] = qx * qx + qy * qy + qz * qz;
          kernel[c] = std::exp(-q2[c] * o.bond * o.bond / 6);
        }
    // the inversion's Newton matrix: the chains' RPA response δθ_I = −Σ_J R_IJ(q) δU_J, inverted per wave vector (types
    // of one block copolymer are strongly coupled: the diagonal alone oscillates)
    rinv.assign(M * nt * nt, 0.0);
    std::vector<double> R(nt * nt), A(nt * nt);
    for (size_t c = 1; c < M; ++c) {
      std::fill(R.begin(), R.end(), 0.0);
      for (size_t m = 0; m < seq.size(); ++m) {
        const auto& s = seq[m];
        const double w = frac[m] / double(s.size());
        for (size_t a = 0; a < s.size(); ++a)
          for (size_t b = 0; b < s.size(); ++b) R[size_t(s[a]) * nt + size_t(s[b])] += w * std::pow(kernel[c], std::fabs(double(a) - double(b)));
      }
      // Gauss–Jordan on [R | I]
      A = R;
      std::vector<double> inv(nt * nt, 0.0);
      for (size_t i = 0; i < nt; ++i) inv[i * nt + i] = 1;
      bool ok = true;
      for (size_t col = 0; col < nt && ok; ++col) {
        size_t piv = col;
        for (size_t row = col + 1; row < nt; ++row) if (std::fabs(A[row * nt + col]) > std::fabs(A[piv * nt + col])) piv = row;
        if (std::fabs(A[piv * nt + col]) < 1e-14) { ok = false; break; }
        for (size_t k = 0; k < nt; ++k) std::swap(A[col * nt + k], A[piv * nt + k]), std::swap(inv[col * nt + k], inv[piv * nt + k]);
        const double d = A[col * nt + col];
        for (size_t k = 0; k < nt; ++k) A[col * nt + k] /= d, inv[col * nt + k] /= d;
        for (size_t row = 0; row < nt; ++row)
          if (row != col) {
            const double f = A[row * nt + col];
            for (size_t k = 0; k < nt; ++k) A[row * nt + k] -= f * A[col * nt + k], inv[row * nt + k] -= f * inv[col * nt + k];
          }
      }
      if (ok) std::copy(inv.begin(), inv.end(), rinv.begin() + long(c * nt * nt));
    }
    buf.resize(M);
  }
  // K ⋆ f (periodic), in place
  void convolve(std::vector<double>& f) {
    for (size_t c = 0; c < M; ++c) buf[c] = f[c];
    fft3d(buf, nx, ny, nz, -1);
    for (size_t c = 0; c < M; ++c) buf[c] *= kernel[c];
    fft3d(buf, nx, ny, nz, 1);
    const double inv = 1.0 / double(M);
    for (size_t c = 0; c < M; ++c) f[c] = buf[c].real() * inv;
  }
  // θ of every type for the potentials U (type-major)
  std::vector<double> densities(const std::vector<double>& U) {
    const size_t nt = types.size();
    std::vector<double> theta(nt * M, 0.0);
    std::vector<std::vector<double>> W(nt, std::vector<double>(M));
    for (size_t t = 0; t < nt; ++t)
      for (size_t c = 0; c < M; ++c) W[t][c] = std::exp(-std::clamp(U[t * M + c], -60.0, 60.0));
    for (size_t m = 0; m < seq.size(); ++m) {
      const auto& s = seq[m];
      const size_t N = s.size();
      std::vector<std::vector<double>> qf(N, std::vector<double>(M)), qb(N, std::vector<double>(M));
      // each step rescaled by its mean, the scale kept as a logarithm (strong fields would otherwise underflow e^(−N U))
      std::vector<double> lf(N, 0.0), lb(N, 0.0);
      auto rescale = [&](std::vector<double>& q, double& l) {
        double m = 0;
        for (double v : q) m += v;
        m /= double(M);
        if (m > 0 && std::isfinite(m)) { for (auto& v : q) v /= m; l += std::log(m); }
      };
      qf[0] = W[size_t(s[0])];
      rescale(qf[0], lf[0]);
      for (size_t k = 1; k < N; ++k) {
        qf[k] = qf[k - 1];
        lf[k] = lf[k - 1];
        convolve(qf[k]);
        for (size_t c = 0; c < M; ++c) qf[k][c] *= W[size_t(s[k])][c];
        rescale(qf[k], lf[k]);
      }
      qb[N - 1] = W[size_t(s[N - 1])];
      rescale(qb[N - 1], lb[N - 1]);
      for (size_t k = N - 1; k-- > 0;) {
        qb[k] = qb[k + 1];
        lb[k] = lb[k + 1];
        convolve(qb[k]);
        for (size_t c = 0; c < M; ++c) qb[k][c] *= W[size_t(s[k])][c];
        rescale(qb[k], lb[k]);
      }
      // Q = ⟨q_N⟩ = e^(lf_N): bead s gets (f/N) q_s q†_s e^(lf_s + lb_s − lf_N) / w_s
      for (size_t k = 0; k < N; ++k) {
        const size_t t = size_t(s[k]);
        const double norm = frac[m] / double(N) * std::exp(lf[k] + lb[k] - lf[N - 1]);
        for (size_t c = 0; c < M; ++c) theta[t * M + c] += norm * qf[k][c] * qb[k][c] / W[t][c];
      }
    }
    return theta;
  }
  // U such that densities(U) = target, from U (warm start); returns the iterations taken
  int invert(std::vector<double>& U, const std::vector<double>& target, int iterations, double tol) {
    const size_t nt = types.size();
    for (int it = 0; it < iterations; ++it) {
      const auto th = densities(U);
      double err = 0;
      for (size_t i = 0; i < th.size(); ++i) err = std::max(err, std::fabs(th[i] - target[i]));
      if (err < tol) return it;
      // δÛ(q) = λ R(q)⁻¹ ê(q), e = θ̄ ln(θ[U] / θ*): too dense where U is too low, so raise it there. The residual is
      // logarithmic (U ← U + λ ln(θ[U] / θ*) where the chain is local, R → θ̄ at large q): a cell far below its mean answers
      // a field change by about −θ δU, not −θ̄ δU, and a linear residual would step it θ/θ̄ too short each iteration
      std::vector<std::vector<cd>> e(nt);
      for (size_t t = 0; t < nt; ++t) {
        for (size_t c = 0; c < M; ++c) buf[c] = mean[t] * std::log(std::max(th[t * M + c], 1e-300) / std::max(target[t * M + c], 1e-300));
        fft3d(buf, nx, ny, nz, -1);
        e[t] = buf;
      }
      for (size_t t = 0; t < nt; ++t) {
        for (size_t c = 0; c < M; ++c) {
          cd v(0, 0);
          for (size_t u = 0; u < nt; ++u) v += rinv[c * nt * nt + t * nt + u] * e[u][c];
          buf[c] = v;
        }
        fft3d(buf, nx, ny, nz, 1);
        const double inv = 0.8 / double(M);
        for (size_t c = 0; c < M; ++c) U[t * M + c] += std::clamp(buf[c].real() * inv, -2.0, 2.0);   // a Newton step far from the linear regime is capped
      }
    }
    return iterations;
  }
};
}  // namespace

std::vector<double> ddft_chain_densities(const DdftOptions& o, const std::vector<double>& U) {
  Model m(o);
  if (U.size() != m.types.size() * m.M) throw std::invalid_argument("density fields: one potential per type and grid point");
  return m.densities(U);
}

DdftResult ddft_run(const DdftOptions& o) {
  const auto t0 = std::chrono::steady_clock::now();
  Model md(o);
  DdftResult r;
  r.types = md.types;
  r.grid = {md.nx, md.ny, md.nz};
  const size_t nt = md.types.size(), M = md.M;
  if (o.dt <= 0 || o.steps < 0) throw std::invalid_argument("density fields: a positive time step");
  // χ per type pair
  std::vector<double> chi(nt * nt, 0.0);
  for (const auto& [k, v] : o.chi) {
    if (k.size() != 2) throw std::invalid_argument("density fields: χ keys are two bead letters, e.g. AB");
    const auto a = std::find(md.types.begin(), md.types.end(), std::string(1, char(std::toupper(static_cast<unsigned char>(k[0]))))),
               b = std::find(md.types.begin(), md.types.end(), std::string(1, char(std::toupper(static_cast<unsigned char>(k[1])))));
    if (a == md.types.end() || b == md.types.end()) continue;   // a type no molecule has
    const size_t i = size_t(a - md.types.begin()), j = size_t(b - md.types.begin());
    if (i != j) chi[i * nt + j] = chi[j * nt + i] = v;
  }
  // the sharpest interface the χ allow, a / √(6χ) (strong segregation): narrower than the grid, the fields cannot follow it
  double chimax = 0;
  for (double v : chi) chimax = std::max(chimax, v);
  const double width = chimax > 0 ? o.bond / std::sqrt(6 * chimax) : 1e9;
  if (width < 0.35) {
    char b[240];
    std::snprintf(b, sizeof b, "the interface width a/√(6χ) = %.2f grid spacings is narrower than the grid: strong segregation may diverge (lower χ, or more beads per chain with a smaller χ each)", width);
    r.notes.push_back(b);
  }
  // the start: uniform, with a small random ripple when there is no noise to break the symmetry
  std::mt19937_64 rng(o.seed);
  caps::Normal<double> G(0.0, 1.0);
  std::vector<double> theta(nt * M), U(nt * M, 0.0);
  for (size_t t = 0; t < nt; ++t)
    for (size_t c = 0; c < M; ++c) theta[t * M + c] = md.mean[t];
  if (o.noise <= 0) {
    for (size_t t = 0; t < nt; ++t) {
      std::vector<double> d(M);
      double s = 0;
      for (auto& x : d) s += x = 1e-3 * md.mean[t] * G(rng);
      for (size_t c = 0; c < M; ++c) theta[t * M + c] += d[c] - s / double(M);
    }
    r.worst_inverse = md.invert(U, theta, o.inverse_iterations, o.inverse_tol);
  }
  const int dims[3] = {md.nx, md.ny, md.nz};
  const int active = (md.nx > 1) + (md.ny > 1) + (md.nz > 1);   // faces per cell: 2 × active
  auto idx = [&](int i, int j, int k) { return size_t(i) + size_t(md.nx) * (size_t(j) + size_t(md.ny) * size_t(k)); };
  auto record = [&](double tau) {
    r.times.push_back(tau);
    std::vector<float> f(theta.size());
    for (size_t i = 0; i < theta.size(); ++i) f[i] = float(theta[i]);
    r.frames.push_back(std::move(f));
    std::vector<double> ord(nt, 0.0);
    for (size_t t = 0; t < nt; ++t) {
      for (size_t c = 0; c < M; ++c) { const double d = theta[t * M + c] - md.mean[t]; ord[t] += d * d; }
      ord[t] /= double(M);
    }
    r.order.push_back(ord);
    double inter = 0, comp = 0;
    std::vector<std::vector<double>> sm(nt);
    for (size_t t = 0; t < nt; ++t) { sm[t].assign(theta.begin() + long(t * M), theta.begin() + long((t + 1) * M)); md.convolve(sm[t]); }
    for (size_t a = 0; a < nt; ++a)
      for (size_t b = a + 1; b < nt; ++b)
        if (chi[a * nt + b] != 0)
          for (size_t c = 0; c < M; ++c) inter += chi[a * nt + b] * theta[a * M + c] * sm[b][c];
    for (size_t c = 0; c < M; ++c) {
      double s = -1;
      for (size_t t = 0; t < nt; ++t) s += theta[t * M + c];
      comp += 0.5 * o.kappa * s * s;
    }
    r.interaction.push_back(inter / double(M)), r.compressibility.push_back(comp / double(M));
    // the strongest Fourier mode of type 0
    for (size_t c = 0; c < M; ++c) md.buf[c] = theta[c] - md.mean[0];
    fft3d(md.buf, md.nx, md.ny, md.nz, -1);
    size_t best = 0;
    double bv = 0;
    for (size_t c = 1; c < M; ++c) if (std::norm(md.buf[c]) > bv) bv = std::norm(md.buf[c]), best = c;
    r.dominant_spacing.push_back(best > 0 && md.q2[best] > 0 ? 2 * kPi / std::sqrt(md.q2[best]) : 0.0);
  };
  record(0);
  // dynamics
  std::vector<double> mu(nt * M), dtheta(nt * M);
  int clamped = 0;
  const int every = std::max(1, o.frame_every);
  for (int step = 1; step <= o.steps; ++step) {
    if (o.progress && step % 10 == 0 && !o.progress(double(step) / o.steps)) { r.notes.push_back("stopped at step " + std::to_string(step)); break; }
    // μ_I = −U_I + Σ_J χ_IJ (K ⋆ θ_J) + κ (Σθ − 1)
    std::vector<std::vector<double>> sm(nt);
    for (size_t t = 0; t < nt; ++t) { sm[t].assign(theta.begin() + long(t * M), theta.begin() + long((t + 1) * M)); md.convolve(sm[t]); }
    for (size_t c = 0; c < M; ++c) {
      double s = -1;
      for (size_t t = 0; t < nt; ++t) s += theta[t * M + c];
      for (size_t t = 0; t < nt; ++t) {
        double v = -U[t * M + c] + o.kappa * s;
        for (size_t u = 0; u < nt; ++u) v += chi[t * nt + u] * sm[u][c];
        mu[t * M + c] = v;
      }
    }
    // fluxes on the faces: J = −θ_face ∇μ + ζ
    std::fill(dtheta.begin(), dtheta.end(), 0.0);
    for (size_t t = 0; t < nt; ++t)
      for (int d = 0; d < 3; ++d) {
        if (dims[d] == 1) continue;
        for (int k = 0; k < md.nz; ++k)
          for (int j = 0; j < md.ny; ++j)
            for (int i = 0; i < md.nx; ++i) {
              const size_t c = idx(i, j, k);
              const size_t n = d == 0 ? idx((i + 1) % md.nx, j, k) : d == 1 ? idx(i, (j + 1) % md.ny, k) : idx(i, j, (k + 1) % md.nz);
              const double tf = std::max(0.0, 0.5 * (theta[t * M + c] + theta[t * M + n]));
              double J = -tf * (mu[t * M + n] - mu[t * M + c]);
              // the Gaussian noise has no bound: a tail draw would empty a cell. Through any one face it may move at
              // most the share 1/(2 × active) of the poorer neighbour, so the noise alone never takes a cell below zero
              if (o.noise > 0) {
                const double cap = std::min(theta[t * M + c], theta[t * M + n]) / (2 * active * o.dt);
                J += std::clamp(std::sqrt(2 * tf / (o.noise * o.dt)) * G(rng), -cap, cap);
              }
              dtheta[t * M + c] -= o.dt * J;
              dtheta[t * M + n] += o.dt * J;
            }
      }
    // the compressibility term implicitly (a stabilising correction in Fourier space): with s = Σθ, the explicit step
    // θ* gets − Δτ θ̄_I λ(q) κ (s_new − s_old), λ the discrete Laplacian's eigenvalue, s_new solved from its sum; the
    // stiff total-density mode (D ≈ κ) then has no time-step limit
    if (o.kappa > 0) {
      std::vector<cd> sn(M, cd(0, 0)), ss(M, cd(0, 0));
      std::vector<std::vector<cd>> ts(nt, std::vector<cd>(M));
      for (size_t t = 0; t < nt; ++t) {
        for (size_t c = 0; c < M; ++c) md.buf[c] = theta[t * M + c] + dtheta[t * M + c];
        fft3d(md.buf, md.nx, md.ny, md.nz, -1);
        ts[t] = md.buf;
        for (size_t c = 0; c < M; ++c) ss[c] += md.buf[c];
        for (size_t c = 0; c < M; ++c) md.buf[c] = theta[t * M + c];
        fft3d(md.buf, md.nx, md.ny, md.nz, -1);
        for (size_t c = 0; c < M; ++c) sn[c] += md.buf[c];
      }
      double mt = 0;
      for (double v : md.mean) mt += v;
      for (int k = 0; k < md.nz; ++k)
        for (int j = 0; j < md.ny; ++j)
          for (int i = 0; i < md.nx; ++i) {
            const size_t c = idx(i, j, k);
            const double lam = (md.nx > 1 ? 2 - 2 * std::cos(2 * kPi * i / md.nx) : 0.0) + (md.ny > 1 ? 2 - 2 * std::cos(2 * kPi * j / md.ny) : 0.0) +
                               (md.nz > 1 ? 2 - 2 * std::cos(2 * kPi * k / md.nz) : 0.0);
            const double g = o.dt * lam * o.kappa;
            const cd snew = (ss[c] + g * mt * sn[c]) / (1 + g * mt);
            for (size_t t = 0; t < nt; ++t) ts[t][c] -= g * md.mean[t] * (snew - sn[c]);
          }
      for (size_t t = 0; t < nt; ++t) {
        fft3d(ts[t], md.nx, md.ny, md.nz, 1);
        for (size_t c = 0; c < M; ++c) dtheta[t * M + c] = ts[t][c].real() / double(M) - theta[t * M + c];
      }
    }
    for (size_t i = 0; i < theta.size(); ++i) {
      theta[i] += dtheta[i];
      if (theta[i] < 1e-8) theta[i] = 1e-8, ++clamped;
    }
    r.worst_inverse = std::max(r.worst_inverse, md.invert(U, theta, o.inverse_iterations, o.inverse_tol));
    for (size_t i = 0; i < theta.size(); ++i)
      if (!std::isfinite(theta[i]) || !std::isfinite(U[i])) {
        char b[300];
        if (width < 0.5)
          std::snprintf(b, sizeof b, "density fields: the fields diverged at τ %.3g: the interface width a/√(6χ) = %.2f grid spacings is narrower than the grid (lower χ, or more beads per chain with a smaller χ each)",
                        step * o.dt, width);
        else
          std::snprintf(b, sizeof b, "density fields: the fields diverged at τ %.3g; take a smaller time step (now %.3g)", step * o.dt, o.dt);
        throw std::invalid_argument(b);
      }
    if (step % every == 0 || step == o.steps) record(step * o.dt);
  }
  {
    char b[300];
    std::string mols;
    for (size_t m = 0; m < md.seq.size(); ++m) {
      char f[32];
      std::snprintf(f, sizeof f, "%.3g", md.frac[m]);
      mols += (m ? ", " : "") + o.molecules[m].sequence + " (" + std::to_string(md.seq[m].size()) + " beads, " + f + ")";
    }
    std::snprintf(b, sizeof b, "grid %d × %d × %d, bond %.4g · κ_H %.3g · Ω %.3g · Δτ %.3g × %d steps · inversion: at most %d iterations a step", md.nx, md.ny, md.nz, o.bond,
                  o.kappa, o.noise, o.dt, o.steps, r.worst_inverse);
    r.notes.push_back("molecules: " + mols);
    r.notes.push_back(b);
    if (r.worst_inverse >= o.inverse_iterations) r.notes.push_back("the potentials did not always match the densities within the tolerance: take a smaller time step");
    if (clamped > 0) r.notes.push_back(std::to_string(clamped) + " cell values went below zero and were held at 10⁻⁸: take a smaller time step");
  }
  r.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  return r;
}

}  // namespace caps
