// Quasi-harmonic free energy and thermal expansion (qha.hpp).
#include "caps/qha.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <stdexcept>

#include "caps/cell_relax.hpp"

namespace caps {

namespace {

constexpr double kAtmToKcalA3 = 101325.0 * 1e-30 * 6.02214076e23 / 4184.0;   // 1 atm in kcal/mol/Å³
constexpr double kKcalA3ToGPa = 4184.0 / 6.02214076e23 / 1e-30 / 1e9;        // 1 kcal/mol/Å³ in GPa

using Exp = std::array<int, 3>;

// a polynomial in up to three variables, x scaled by `scale` so the normal equations stay well conditioned
struct Fit {
  std::vector<Exp> basis;
  std::vector<double> c;
  double scale = 1;
  double value(const std::array<double, 3>& x) const {
    double v = 0;
    for (size_t b = 0; b < basis.size(); ++b) {
      double t = c[b];
      for (int k = 0; k < 3; ++k) t *= std::pow(x[size_t(k)] / scale, basis[b][size_t(k)]);
      v += t;
    }
    return v;
  }
  // gradient and Hessian in x
  void derivatives(const std::array<double, 3>& x, std::array<double, 3>& g, std::array<std::array<double, 3>, 3>& H) const {
    g = {0, 0, 0};
    H = {};
    std::array<double, 3> u{x[0] / scale, x[1] / scale, x[2] / scale};
    auto pw = [](double a, int e) { return e <= 0 ? 1.0 : std::pow(a, e); };
    for (size_t b = 0; b < basis.size(); ++b) {
      const Exp& e = basis[b];
      for (int i = 0; i < 3; ++i) {
        if (e[size_t(i)] > 0) {
          double t = c[b] * e[size_t(i)] * pw(u[size_t(i)], e[size_t(i)] - 1);
          for (int k = 0; k < 3; ++k) if (k != i) t *= pw(u[size_t(k)], e[size_t(k)]);
          g[size_t(i)] += t / scale;
        }
        for (int j = i; j < 3; ++j) {
          double t = c[b];
          if (i == j) {
            if (e[size_t(i)] < 2) continue;
            t *= e[size_t(i)] * (e[size_t(i)] - 1) * pw(u[size_t(i)], e[size_t(i)] - 2);
            for (int k = 0; k < 3; ++k) if (k != i) t *= pw(u[size_t(k)], e[size_t(k)]);
          } else {
            if (e[size_t(i)] < 1 || e[size_t(j)] < 1) continue;
            t *= e[size_t(i)] * pw(u[size_t(i)], e[size_t(i)] - 1) * e[size_t(j)] * pw(u[size_t(j)], e[size_t(j)] - 1);
            const int k = 3 - i - j;
            t *= pw(u[size_t(k)], e[size_t(k)]);
          }
          H[size_t(i)][size_t(j)] += t / (scale * scale);
          if (i != j) H[size_t(j)][size_t(i)] += t / (scale * scale);
        }
      }
    }
  }
};

// solves A z = b by Gaussian elimination with pivoting; false when singular
bool solve(std::vector<std::vector<double>> A, std::vector<double> b, std::vector<double>& z) {
  const size_t n = b.size();
  for (size_t i = 0; i < n; ++i) {
    size_t p = i;
    for (size_t r = i + 1; r < n; ++r) if (std::fabs(A[r][i]) > std::fabs(A[p][i])) p = r;
    if (std::fabs(A[p][i]) < 1e-14) return false;
    std::swap(A[i], A[p]), std::swap(b[i], b[p]);
    for (size_t r = i + 1; r < n; ++r) {
      const double f = A[r][i] / A[i][i];
      for (size_t c = i; c < n; ++c) A[r][c] -= f * A[i][c];
      b[r] -= f * b[i];
    }
  }
  z.assign(n, 0);
  for (size_t i = n; i-- > 0;) {
    double t = b[i];
    for (size_t c = i + 1; c < n; ++c) t -= A[i][c] * z[c];
    z[i] = t / A[i][i];
  }
  return true;
}

Fit fit_poly(const std::vector<std::array<double, 3>>& x, const std::vector<double>& y, const std::vector<Exp>& basis, double scale, double* rms = nullptr) {
  Fit f;
  f.basis = basis, f.scale = scale;
  const size_t m = basis.size();
  std::vector<std::vector<double>> A(m, std::vector<double>(m, 0.0));
  std::vector<double> b(m, 0.0);
  std::vector<double> row(m);
  // about the mean value: the fit is of small differences on a large static energy
  double mean = 0;
  for (double v : y) mean += v;
  mean /= double(y.size());
  for (size_t s = 0; s < x.size(); ++s) {
    for (size_t k = 0; k < m; ++k) {
      double t = 1;
      for (int d = 0; d < 3; ++d) t *= std::pow(x[s][size_t(d)] / scale, basis[k][size_t(d)]);
      row[k] = t;
    }
    for (size_t i = 0; i < m; ++i) {
      b[i] += row[i] * (y[s] - mean);
      for (size_t j = 0; j < m; ++j) A[i][j] += row[i] * row[j];
    }
  }
  if (!solve(A, b, f.c)) throw std::runtime_error("quasi-harmonic: the strained cells do not determine the fit");
  for (size_t k = 0; k < m; ++k)
    if (basis[k] == Exp{0, 0, 0}) f.c[k] += mean;
  if (rms) {
    double v = 0;
    for (size_t s = 0; s < x.size(); ++s) { const double d = f.value(x[s]) - y[s]; v += d * d; }
    *rms = std::sqrt(v / double(x.size()));
  }
  return f;
}

// true when the cells at x fix every coefficient of the basis (the normal equations are not singular)
bool determines(const std::vector<std::array<double, 3>>& x, const std::vector<Exp>& basis, double scale) {
  const size_t m = basis.size();
  std::vector<std::vector<double>> A(m, std::vector<double>(m, 0.0));
  std::vector<double> row(m), z;
  for (const auto& p : x) {
    for (size_t k = 0; k < m; ++k) {
      double t = 1;
      for (int d = 0; d < 3; ++d) t *= std::pow(p[size_t(d)] / scale, basis[k][size_t(d)]);
      row[k] = t;
    }
    for (size_t i = 0; i < m; ++i)
      for (size_t j = 0; j < m; ++j) A[i][j] += row[i] * row[j];
  }
  return solve(A, std::vector<double>(m, 0.0), z);
}

}  // namespace

double qha_step(double T) { return T >= 2 ? 1.0 : T > 0 ? T / 2 : 1.0; }

QhaResult qha_analyze(const std::vector<QhaSample>& samples, const std::vector<double>& temperatures, const std::string& mode,
                      const std::array<bool, 3>& variable, const std::array<double, 3>& lengths0, double volume0, double pressure_atm) {
  QhaResult r;
  r.mode = mode;
  r.samples = samples;
  r.lengths0 = lengths0, r.volume0 = volume0;
  if (samples.empty()) throw std::invalid_argument("quasi-harmonic: no strained cells");
  const bool vol = mode == "volume";
  if (!vol && mode != "axes") throw std::invalid_argument("quasi-harmonic: mode volume or axes");
  // the fit's variables: x₀ = V/V₀ − 1 (volume), else the free linear strains
  std::vector<int> free;
  if (vol) free = {0};
  else for (int k = 0; k < 3; ++k) if (variable[size_t(k)]) free.push_back(k);
  if (free.empty()) throw std::invalid_argument("quasi-harmonic: no axis is free to change");
  // The cells the fit uses: those with no imaginary mesh mode. An unstable mode is left out of F_vib, where its soft
  // neighbours at the next strain add a large negative k T ln(h nu / k T): F_vib at that cell jumps up by an amount that
  // grows with T, steepens G there and holds the minimum inside the grid (PE at +2 % strain: B_T rising from 37 to
  // 100 GPa on heating, alpha falling). Only when too few stable cells remain for any fit are all cells used, and then
  // every row is marked unstable.
  auto coords = [&](const std::vector<size_t>& idx) {
    std::vector<std::array<double, 3>> out;
    for (size_t s : idx) {
      std::array<double, 3> x{0, 0, 0};
      if (vol) x[0] = samples[s].volume / volume0 - 1;
      else for (int k : free) x[size_t(k)] = samples[s].strain[size_t(k)];
      out.push_back(x);
    }
    return out;
  };
  // A least-squares basis for the cells at x, with at least one more cell than terms: a polynomial through every point
  // carries any noise in F_vib straight into the curvature (and so into B_T). Empty when no fit is possible.
  auto basis_for = [&](const std::vector<std::array<double, 3>>& x, double scale) {
    std::vector<Exp> b{{0, 0, 0}};
    const int n = int(x.size());
    if (vol) {
      // G as a polynomial in V/V0 - 1, degree up to 4: 4 cells a quadratic, 6 a cubic, 7 and more a quartic
      if (n < 4) return std::vector<Exp>{};
      const int deg = std::max(2, std::min(4, n - 3));
      for (int e = 1; e <= deg; ++e) b.push_back({e, 0, 0});
      return determines(x, b, scale) ? b : std::vector<Exp>{};
    }
    // every term up to third order, with x_i^4 and x_i^2 x_j^2 (crystals are strongly anharmonic over a few per cent);
    // when the cells cannot carry them, the x_i^2 x_j^2, then the x_i^4, then the mixed and then the pure cubic terms go,
    // each set whole so no axis is favoured; the quadratic surface always stays
    const size_t nf = free.size();
    std::vector<std::pair<int, Exp>> terms;   // (the order they go in, 4 = kept), exponents
    auto term = [&](int a, int bb, int c) {   // exponents of the free variables, in order
      Exp e{0, 0, 0};
      e[size_t(free[0])] = a;
      if (nf > 1) e[size_t(free[1])] = bb;
      if (nf > 2) e[size_t(free[2])] = c;
      const int tot = a + bb + c, nz = (a > 0) + (bb > 0) + (c > 0);
      const int cls = tot == 4 ? (nz > 1 ? 0 : 1) : tot == 3 ? (nz > 1 ? 2 : 3) : 4;
      terms.push_back({cls, e});
    };
    for (int a = 0; a <= 4; ++a)
      for (int bb = 0; bb <= (nf > 1 ? 4 : 0); ++bb)
        for (int c = 0; c <= (nf > 2 ? 4 : 0); ++c) {
          const int tot = a + bb + c;
          if (tot == 0) continue;
          const bool quartic = tot == 4 && ((a == 4 || bb == 4 || c == 4) || (a == 2 && bb == 2) || (a == 2 && c == 2) || (bb == 2 && c == 2));
          if (tot <= 3 || quartic) term(a, bb, c);
        }
    for (int drop = 0; drop <= 4; ++drop) {   // drop: the classes below it are left out
      b.assign(1, Exp{0, 0, 0});
      for (const auto& t : terms) if (t.first >= drop) b.push_back(t.second);
      if (int(b.size()) < n && determines(x, b, scale)) return b;
    }
    return std::vector<Exp>{};
  };
  auto span_of = [](const std::vector<std::array<double, 3>>& x) {
    double v = 0;
    for (const auto& p : x) for (double c : p) v = std::max(v, std::fabs(c));
    return v;
  };
  std::vector<size_t> use, every;
  for (size_t s = 0; s < samples.size(); ++s) {
    every.push_back(s);
    if (samples[s].imaginary == 0) use.push_back(s);
  }
  std::vector<std::array<double, 3>> X = coords(use);
  double span = span_of(X);
  std::vector<Exp> basis = span > 0 ? basis_for(X, span) : std::vector<Exp>{};
  const size_t left_out = samples.size() - use.size();
  const bool fallback = basis.empty() && left_out > 0;   // too few stable cells: all of them, every row marked
  if (fallback) {
    use = every, X = coords(use), span = span_of(X);
    if (span > 0) basis = basis_for(X, span);
  }
  if (span <= 0) throw std::invalid_argument("quasi-harmonic: the strained cells are all the same");
  if (basis.empty())
    throw std::invalid_argument(vol ? "quasi-harmonic: at least four volumes" : "quasi-harmonic: too few strained cells for the fit");
  double lo[3] = {1e300, 1e300, 1e300}, hi[3] = {-1e300, -1e300, -1e300};
  for (const auto& x : X)
    for (int k = 0; k < 3; ++k) lo[k] = std::min(lo[k], x[size_t(k)]), hi[k] = std::max(hi[k], x[size_t(k)]);
  const double P = pressure_atm * kAtmToKcalA3;
  auto volume_at = [&](const std::array<double, 3>& x) {
    if (vol) return volume0 * (1 + x[0]);
    return volume0 * (1 + x[0]) * (1 + x[1]) * (1 + x[2]);
  };
  struct Min { std::array<double, 3> x{0, 0, 0}; double V = 0, G = 0; std::array<std::array<double, 3>, 3> H{}; };
  // the minimum of G at one of the three temperatures about a requested one
  auto minimum = [&](size_t it, int j) {
    std::vector<double> y;
    size_t best = 0;
    for (size_t u = 0; u < use.size(); ++u) {
      const QhaSample& s = samples[use[u]];
      y.push_back(s.energy + s.F[it][size_t(j)] + P * s.volume);
      if (y[u] < y[best]) best = u;
    }
    double rms = 0;
    const Fit f = fit_poly(X, y, basis, span, &rms);
    r.fit_rms = std::max(r.fit_rms, rms);
    Min m;
    m.x = X[best];
    std::array<double, 3> g;
    for (int it2 = 0; it2 < 60; ++it2) {   // Newton on the free variables, steps limited to the grid's span
      f.derivatives(m.x, g, m.H);
      const size_t n = free.size();
      std::vector<std::vector<double>> A(n, std::vector<double>(n));
      std::vector<double> b(n), z;
      for (size_t a = 0; a < n; ++a) {
        b[a] = -g[size_t(free[a])];
        for (size_t c = 0; c < n; ++c) A[a][c] = m.H[size_t(free[a])][size_t(free[c])];
      }
      if (!solve(A, b, z)) break;
      double step = 0;
      for (double v : z) step = std::max(step, std::fabs(v));
      const double lim = step > span ? span / step : 1.0;
      for (size_t a = 0; a < n; ++a) m.x[size_t(free[a])] += lim * z[a];
      if (step < 1e-12) break;
    }
    f.derivatives(m.x, g, m.H);
    m.V = volume_at(m.x);
    m.G = f.value(m.x);
    return m;
  };
  auto lengths_of = [&](const std::array<double, 3>& x) {
    std::array<double, 3> L = lengths0;
    if (vol) { const double q = std::cbrt(1 + x[0]); for (auto& l : L) l *= q; }
    else for (int k = 0; k < 3; ++k) L[size_t(k)] *= 1 + x[size_t(k)];
    return L;
  };
  bool any_out = false, unstable = false;
  const double margin = 0.05 * span;
  auto outside = [&](const std::array<double, 3>& x) {
    for (int k : free) if (x[size_t(k)] < lo[k] - margin || x[size_t(k)] > hi[k] + margin) return true;
    return false;
  };
  for (size_t it = 0; it < temperatures.size(); ++it) {
    const double T = temperatures[it], dT = qha_step(T);
    QhaRow row;
    row.T = T;
    const Min m = minimum(it, T > 0 ? 1 : 0);
    row.volume = m.V, row.G = m.G, row.length = lengths_of(m.x);
    if (vol) { const double q = std::cbrt(1 + m.x[0]) - 1; row.strain = {q, q, q}; }
    else row.strain = m.x;
    row.extrapolated = outside(m.x);
    // the curvature must be a minimum along every free variable
    for (int k : free) if (m.H[size_t(k)][size_t(k)] <= 0) row.unstable = true;
    if (fallback) row.unstable = true;
    if (T > 0) {
      const Min a = minimum(it, 0), b = minimum(it, 2);
      if (outside(a.x) || outside(b.x)) row.extrapolated = true;   // alpha leans on the minima at T +/- dT too
      row.alpha_V = (b.V - a.V) / (2 * dT * m.V);
      const auto La = lengths_of(a.x), Lb = lengths_of(b.x);
      for (int k = 0; k < 3; ++k) row.alpha[size_t(k)] = (Lb[size_t(k)] - La[size_t(k)]) / (2 * dT * row.length[size_t(k)]);
    }
    any_out = any_out || row.extrapolated;
    // B_T = V / (vᵀ H⁻¹ v), v = ∂V/∂x
    {
      const size_t n = free.size();
      std::vector<std::vector<double>> A(n, std::vector<double>(n));
      std::vector<double> v(n), z;
      for (size_t i = 0; i < n; ++i) {
        v[i] = vol ? volume0 : m.V / (1 + m.x[size_t(free[i])]);
        for (size_t c = 0; c < n; ++c) A[i][c] = m.H[size_t(free[i])][size_t(free[c])];
      }
      if (solve(A, v, z)) {
        double q = 0;
        for (size_t i = 0; i < n; ++i) q += v[i] * z[i];
        if (q > 0) row.B = m.V / q * kKcalA3ToGPa;
        else row.unstable = true;   // G curves down along some combination of the strains
      } else {
        row.unstable = true;
      }
    }
    unstable = unstable || row.unstable;
    // C_v and S at the minimum, from fits of their own over the strained cells
    std::vector<double> yc, ys;
    for (size_t u : use) yc.push_back(samples[u].Cv[it][T > 0 ? 1 : 0]), ys.push_back(samples[u].S[it][T > 0 ? 1 : 0]);
    row.Cv = fit_poly(X, yc, basis, span).value(m.x);
    row.S = fit_poly(X, ys, basis, span).value(m.x);
    row.Cp = row.Cv + T * m.V * (row.B / kKcalA3ToGPa) * row.alpha_V * row.alpha_V * 1000;   // cal/mol/K
    r.rows.push_back(row);
  }
  int imag = 0;
  double wmin = 0;
  for (const auto& q : samples) imag = std::max(imag, q.imaginary), wmin = std::min(wmin, q.imaginary > 0 ? q.min_wavenumber : 0.0);
  char b[400];
  if (left_out > 0 && !fallback) {
    std::snprintf(b, sizeof b, "%zu of %zu strained cells have imaginary mesh modes (up to %d, down to %.1f cm⁻¹) and are left out of the fit: "
                  "their F_vib lacks those modes", left_out, samples.size(), imag, wmin);
    r.notes.push_back(b);
  }
  if (fallback) {
    std::snprintf(b, sizeof b, "%zu of %zu strained cells have imaginary mesh modes (up to %d, down to %.1f cm⁻¹) and too few stable cells "
                  "remain for a fit: all cells are used and every row is marked unstable", left_out, samples.size(), imag, wmin);
    r.notes.push_back(b);
  }
  if (any_out) {
    if (left_out > 0 && !fallback)
      r.notes.push_back("a minimum lies beyond the cells without imaginary modes: the lattice is not stable there along this strain path, so "
                        "the rows marked extrapolated are not a stable state" + std::string(vol ? " (for an anisotropic crystal try the axes mode)" : ""));
    else
      r.notes.push_back("a minimum lies outside the strains sampled: the cell there is extrapolated; widen the range");
  }
  if (unstable) r.notes.push_back("the fitted free energy is not a minimum along every axis at some temperature, or rests on unstable cells: rows marked unstable");
  return r;
}

QhaResult quasi_harmonic(const System& cell0, const QhaOptions& o) {
  if (!cell0.cell.valid()) throw std::invalid_argument("quasi-harmonic: a periodic cell is needed");
  if (!o.phonon.field_for) throw std::invalid_argument("quasi-harmonic: a force field is needed");
  if (o.range <= 0 || o.range > 0.2) throw std::invalid_argument("quasi-harmonic: a strain range between 0 and 0.2");
  if (o.temperatures.empty()) throw std::invalid_argument("quasi-harmonic: give the temperatures");
  for (double T : o.temperatures) if (T < 0) throw std::invalid_argument("quasi-harmonic: temperatures from 0 K");
  const bool vol = o.mode == "volume";
  if (!vol && o.mode != "axes") throw std::invalid_argument("quasi-harmonic: mode volume or axes");
  if (vol && o.points < 4) throw std::invalid_argument("quasi-harmonic: at least four volumes");
  auto stop = [&](const std::string& what, double f) {
    if (o.progress && !o.progress(what, f)) throw std::runtime_error("quasi-harmonic cancelled");
  };
  System base = cell0;
  std::vector<std::string> notes;
  // the free energy surface must be smooth in the strains: plain truncation with a tail term jumps each time a pair crosses
  // the cut-off, which swamps the small energies of a few per cent strain; the shifted form (no tail) is continuous
  QhaOptions oq = o;
  if (oq.phonon.energy.tail) notes.push_back("Lennard-Jones shifted at the cut-off, no tail term: the surface must be smooth in the strains (a truncated potential jumps as pairs cross the cut-off)");
  oq.phonon.energy.tail = false;
  // the same for charges: the damped shifted force sum is a cut-off method, and its curvature jumps each time a shell of
  // ions crosses the cut-off as the cell is strained (MgO: the bulk modulus rose with temperature); Ewald (PME) is smooth
  if (oq.phonon.energy.electrostatics == EnergyOptions::Electrostatics::DSF) {
    const ForceField f0 = oq.phonon.field_for(base);
    if (std::any_of(f0.charge.begin(), f0.charge.end(), [](double q) { return q != 0; })) {
      oq.phonon.energy.electrostatics = EnergyOptions::Electrostatics::PME;
      oq.phonon.energy.ewald_rtol = std::min(oq.phonon.energy.ewald_rtol, 1e-6);
      oq.phonon.energy.pme_spacing = std::min(oq.phonon.energy.pme_spacing, 0.6);
      notes.push_back("Coulomb by Ewald (PME) for the free-energy surface: the damped shifted force sum's cut-off puts kinks in G(V)");
    }
  }
  if (o.relax_cell) {
    stop("relaxing the cell", 0);
    CellRelaxOptions co;
    co.pressure = o.pressure, co.ftol = 0.02;
    if (!vol) for (int k = 0; k < 3; ++k) co.hold[size_t(k)] = !o.axes[size_t(k)];
    const CellRelaxReport cr = relax_cell(base, oq.phonon.field_for, oq.phonon.energy, co);
    char b[200];
    std::snprintf(b, sizeof b, "the grid is centred on the relaxed static cell: V %.2f → %.2f Å³", cr.volume_before, cr.volume_after);
    notes.push_back(b);
    if (!cr.converged) notes.push_back("the static cell did not converge fully: the grid is centred near it");
  }
  // the strains
  std::vector<std::array<double, 3>> grid;
  if (vol) {
    for (int k = 0; k < o.points; ++k) {
      const double e = -o.range + 2 * o.range * k / (o.points - 1);
      grid.push_back({e, e, e});
    }
  } else {
    const double d = o.range / 2;
    grid.push_back({0, 0, 0});
    std::vector<int> fr;
    for (int k = 0; k < 3; ++k) if (o.axes[size_t(k)]) fr.push_back(k);
    if (fr.empty()) throw std::invalid_argument("quasi-harmonic: no axis is free to change");
    for (int k : fr)
      for (double e : {-2 * d, -d, d, 2 * d}) {
        std::array<double, 3> s{0, 0, 0};
        s[size_t(k)] = e;
        grid.push_back(s);
      }
    for (size_t i = 0; i < fr.size(); ++i)
      for (size_t j = i + 1; j < fr.size(); ++j)
        for (double e1 : {-d, d})
          for (double e2 : {-d, d}) {
            std::array<double, 3> s{0, 0, 0};
            s[size_t(fr[i])] = e1, s[size_t(fr[j])] = e2;
            grid.push_back(s);
          }
    if (fr.size() == 3)
      for (double e : {-2 * d, -d, d, 2 * d}) grid.push_back({e, e, e});   // the hydrostatic path
    // the cell count is fixed by that design; the points option only spaces the volume grid
    if (o.points != QhaOptions{}.points)
      notes.push_back("the points option applies to the volume mode only: the axes mode uses its own grid of " + std::to_string(grid.size()) +
                      " cells (the centre, 4 along each free axis, 4 for each pair of free axes, 4 along the diagonal when all three are free)");
  }
  // three temperatures about each requested one (0 K: 0, ΔT, 2ΔT)
  std::vector<double> temps;
  for (double T : o.temperatures) {
    const double dT = qha_step(T);
    if (T > 0) temps.insert(temps.end(), {T - dT, T, T + dT});
    else temps.insert(temps.end(), {0.0, dT, 2 * dT});
  }
  const std::array<double, 3> L0{norm(base.cell.a), norm(base.cell.b), norm(base.cell.c)};
  const double V0 = base.cell.volume();
  std::string field;
  std::vector<QhaSample> samples;
  for (size_t g = 0; g < grid.size(); ++g) {
    const double f0 = double(g) / double(grid.size());
    stop("cell " + std::to_string(g + 1) + " of " + std::to_string(grid.size()), f0);
    // scale the lattice with the atoms affinely
    System s = base;
    const std::array<double, 3>& e = grid[g];
    Cell c = base.cell;
    c.a = base.cell.a * (1 + e[0]), c.b = base.cell.b * (1 + e[1]), c.c = base.cell.c * (1 + e[2]);
    for (auto& a : s.atoms) {
      const Vec3 f = base.cell.to_fractional(a.pos);
      a.pos = c.origin + c.a * f[0] + c.b * f[1] + c.c * f[2];
    }
    s.cell = c;
    PhononOptions po = oq.phonon;
    po.temperatures = temps;
    po.path_points = 2;                 // only the mesh is used
    po.progress = [&](const std::string& what, double f) {
      if (!o.progress) return true;
      return o.progress("cell " + std::to_string(g + 1) + "/" + std::to_string(grid.size()) + ": " + what, f0 + f / double(grid.size()));
    };
    const PhononResult pr = phonons(s, po);
    if (pr.T.size() != temps.size()) throw std::runtime_error("quasi-harmonic: the phonon thermodynamics are incomplete");
    field = pr.field;
    QhaSample q;
    q.strain = e;
    q.volume = c.volume();
    q.energy = pr.static_energy;
    q.imaginary = pr.imaginary, q.min_wavenumber = pr.min_wavenumber;
    for (size_t t = 0; t < o.temperatures.size(); ++t) {
      q.F.push_back({pr.F[3 * t], pr.F[3 * t + 1], pr.F[3 * t + 2]});
      q.Cv.push_back({pr.Cv[3 * t], pr.Cv[3 * t + 1], pr.Cv[3 * t + 2]});
      q.S.push_back({pr.S[3 * t], pr.S[3 * t + 1], pr.S[3 * t + 2]});
    }
    samples.push_back(std::move(q));
  }
  QhaResult r = qha_analyze(samples, o.temperatures, o.mode, o.axes, L0, V0, o.pressure);
  r.field = field;
  r.notes.insert(r.notes.begin(), notes.begin(), notes.end());
  stop("done", 1);
  return r;
}

}  // namespace caps
