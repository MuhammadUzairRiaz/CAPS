// CAPS mechanics and thermal transitions (see mechanics.hpp for the methods and references).
#include "caps/uff.hpp"
#include "caps/mechanics.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <limits>
#include <map>
#include <numeric>
#include <random>
#include <stdexcept>
#include <tuple>

#include "caps/relax.hpp"

namespace caps {

namespace {

constexpr double kGPa = 4184.0 / 6.02214076e23 / 1e-30 / 1e9;   // kcal/(mol·Å³) → GPa
constexpr double kB = 0.0019872041;                              // kcal/(mol·K)
constexpr double kAtmMPa = 0.101325;                             // atm → MPa
constexpr double NaN = std::numeric_limits<double>::quiet_NaN();

using M3 = std::array<std::array<double, 3>, 3>;

// Voigt index → Cartesian pair, and → index of the virial tensor (xx yy zz xy xz yz)
constexpr int VA[6] = {0, 1, 2, 1, 0, 0}, VB[6] = {0, 1, 2, 2, 2, 1}, VW[6] = {0, 1, 2, 5, 4, 3};

M3 identity() { return {{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}}; }

M3 mul(const M3& a, const M3& b) {
  M3 c{};
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      for (int k = 0; k < 3; ++k) c[i][j] += a[i][k] * b[k][j];
  return c;
}

M3 inverse(const M3& m) {
  const double det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
                     m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
  M3 r;
  r[0][0] = (m[1][1] * m[2][2] - m[1][2] * m[2][1]) / det;
  r[0][1] = (m[0][2] * m[2][1] - m[0][1] * m[2][2]) / det;
  r[0][2] = (m[0][1] * m[1][2] - m[0][2] * m[1][1]) / det;
  r[1][0] = (m[1][2] * m[2][0] - m[1][0] * m[2][2]) / det;
  r[1][1] = (m[0][0] * m[2][2] - m[0][2] * m[2][0]) / det;
  r[1][2] = (m[0][2] * m[1][0] - m[0][0] * m[1][2]) / det;
  r[2][0] = (m[1][0] * m[2][1] - m[1][1] * m[2][0]) / det;
  r[2][1] = (m[0][1] * m[2][0] - m[0][0] * m[2][1]) / det;
  r[2][2] = (m[0][0] * m[1][1] - m[0][1] * m[1][0]) / det;
  return r;
}

// Square root of a symmetric positive-definite matrix (Denman–Beavers iteration).
M3 sqrtm(const M3& a) {
  M3 y = a, z = identity();
  for (int it = 0; it < 40; ++it) {
    const M3 yi = inverse(y), zi = inverse(z);
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) {
        y[i][j] = 0.5 * (y[i][j] + zi[i][j]);
        z[i][j] = 0.5 * (z[i][j] + yi[i][j]);
      }
  }
  return y;
}

// Symmetric strain tensor for Voigt component v of magnitude e (shear: engineering, so the tensor holds e/2).
M3 voigt_strain(int v, double e) {
  M3 m{};
  if (v < 3) m[v][v] = e;
  else { m[VA[v]][VB[v]] = e / 2; m[VB[v]][VA[v]] = e / 2; }
  return m;
}

Vec3 mapv(const M3& f, const Vec3& r) {
  return {f[0][0] * r[0] + f[0][1] * r[1] + f[0][2] * r[2], f[1][0] * r[0] + f[1][1] * r[1] + f[1][2] * r[2],
          f[2][0] * r[0] + f[2][1] * r[1] + f[2][2] * r[2]};
}

// Maps positions and cell by the deformation gradient F about the cell origin.
void deform(const M3& f, const std::vector<double>& x, const Cell& c, std::vector<double>& y, Cell& d) {
  y.resize(x.size());
  for (size_t i = 0; i < x.size(); i += 3) {
    const Vec3 r = mapv(f, Vec3{x[i] - c.origin[0], x[i + 1] - c.origin[1], x[i + 2] - c.origin[2]});
    for (int k = 0; k < 3; ++k) y[i + k] = c.origin[k] + r[k];
  }
  d = c;
  d.a = mapv(f, c.a);
  d.b = mapv(f, c.b);
  d.c = mapv(f, c.c);
}

M3 virial_matrix(const EnergyTerms& e) {
  return {{{e.w[0], e.w[3], e.w[4]}, {e.w[3], e.w[1], e.w[5]}, {e.w[4], e.w[5], e.w[2]}}};
}

std::vector<double> positions(const System& s) {
  std::vector<double> x;
  x.reserve(3 * s.atoms.size());
  for (const auto& a : s.atoms) x.insert(x.end(), a.pos.begin(), a.pos.end());
  return x;
}

std::string fmt(double v, int digits) {
  char b[32];
  std::snprintf(b, sizeof b, "%.*g", digits, v);
  return b;
}

// Least-squares line y = a + b x with the standard error of b.
void linfit(const std::vector<double>& x, const std::vector<double>& y, double& a, double& b, double& berr) {
  const size_t n = x.size();
  a = b = 0;
  berr = NaN;
  if (n < 2) return;
  double mx = 0, my = 0;
  for (size_t i = 0; i < n; ++i) { mx += x[i]; my += y[i]; }
  mx /= n;
  my /= n;
  double sxx = 0, sxy = 0;
  for (size_t i = 0; i < n; ++i) { sxx += (x[i] - mx) * (x[i] - mx); sxy += (x[i] - mx) * (y[i] - my); }
  if (sxx <= 0) return;
  b = sxy / sxx;
  a = my - b * mx;
  if (n > 2) {
    double rss = 0;
    for (size_t i = 0; i < n; ++i) { const double r = y[i] - a - b * x[i]; rss += r * r; }
    berr = std::sqrt(rss / (n - 2) / sxx);
  }
}

bool cancelled(const std::function<bool(const std::string&, double)>& p, const std::string& what, double f) {
  return p && !p(what, f);
}

}  // namespace

// ---------------------------------------------------------------- isotropic averages

void isotropic_averages(ElasticResult& r) {
  const auto& C = r.C;
  r.K_voigt = ((C[0][0] + C[1][1] + C[2][2]) + 2 * (C[0][1] + C[0][2] + C[1][2])) / 9;
  r.G_voigt = ((C[0][0] + C[1][1] + C[2][2]) - (C[0][1] + C[0][2] + C[1][2]) + 3 * (C[3][3] + C[4][4] + C[5][5])) / 15;
  // compliance S = C⁻¹ (Gauss–Jordan, 6 × 6)
  double a[6][12];
  for (int i = 0; i < 6; ++i)
    for (int j = 0; j < 12; ++j) a[i][j] = j < 6 ? C[i][j] : (j - 6 == i ? 1 : 0);
  bool singular = false;
  for (int c = 0; c < 6 && !singular; ++c) {
    int p = c;
    for (int i = c + 1; i < 6; ++i) if (std::fabs(a[i][c]) > std::fabs(a[p][c])) p = i;
    if (std::fabs(a[p][c]) < 1e-12) { singular = true; break; }
    for (int j = 0; j < 12; ++j) std::swap(a[c][j], a[p][j]);
    const double d = a[c][c];
    for (int j = 0; j < 12; ++j) a[c][j] /= d;
    for (int i = 0; i < 6; ++i)
      if (i != c) {
        const double m = a[i][c];
        for (int j = 0; j < 12; ++j) a[i][j] -= m * a[c][j];
      }
  }
  if (singular) {
    r.K_reuss = r.G_reuss = NaN;
  } else {
    auto S = [&](int i, int j) { return a[i][j + 6]; };
    r.K_reuss = 1 / ((S(0, 0) + S(1, 1) + S(2, 2)) + 2 * (S(0, 1) + S(0, 2) + S(1, 2)));
    r.G_reuss = 15 / (4 * (S(0, 0) + S(1, 1) + S(2, 2)) - 4 * (S(0, 1) + S(0, 2) + S(1, 2)) + 3 * (S(3, 3) + S(4, 4) + S(5, 5)));
  }
  r.K_hill = std::isnan(r.K_reuss) ? r.K_voigt : 0.5 * (r.K_voigt + r.K_reuss);
  r.G_hill = std::isnan(r.G_reuss) ? r.G_voigt : 0.5 * (r.G_voigt + r.G_reuss);
  const double K = r.K_hill, G = r.G_hill;
  r.E_hill = 9 * K * G / (3 * K + G);
  r.nu_hill = (3 * K - 2 * G) / (2 * (3 * K + G));
  r.lambda_hill = K - 2 * G / 3;
}

// ---------------------------------------------------------------- stress helpers

std::array<double, 6> virial_stress(const EnergyTerms& e, double volume) {
  std::array<double, 6> s{};
  for (int v = 0; v < 6; ++v) s[v] = -e.w[VW[v]] / volume;
  return s;
}

Mat6 born_matrix(Evaluator& ev, const std::vector<double>& x, const Cell& cell, double h) {
  Mat6 B{};
  const double v0 = cell.volume();
  ev.freeze_pairs(x, cell);
  std::vector<double> y, f;
  Cell d;
  try {
    for (int J = 0; J < 6; ++J) {
      std::array<double, 6> sp{}, sm{};
      for (int sg : {1, -1}) {
        const M3 eta = voigt_strain(J, sg * h);
        M3 m = identity();
        for (int i = 0; i < 3; ++i)
          for (int j = 0; j < 3; ++j) m[i][j] += 2 * eta[i][j];
        const M3 F = sqrtm(m), Fi = inverse(F);
        deform(F, x, cell, y, d);
        const M3 S = mul(mul(Fi, virial_matrix(ev.compute(y, d, f))), Fi);   // −S (second Piola–Kirchhoff virial)
        auto& out = sg > 0 ? sp : sm;
        for (int I = 0; I < 6; ++I) out[I] = -S[VA[I]][VB[I]];
      }
      for (int I = 0; I < 6; ++I) B[I][J] = (sp[I] - sm[I]) / (2 * h) / v0 * kGPa;
    }
  } catch (...) {
    ev.unfreeze();
    throw;
  }
  ev.unfreeze();
  for (int I = 0; I < 6; ++I)
    for (int J = I + 1; J < 6; ++J) B[I][J] = B[J][I] = 0.5 * (B[I][J] + B[J][I]);
  return B;
}

// ---------------------------------------------------------------- static strain

ElasticResult static_elastic(const std::vector<System>& configs, const StaticElasticOptions& o) {
  ElasticResult res;
  if (configs.empty()) throw std::invalid_argument("no configurations");
  std::vector<Mat6> all;
  std::array<double, 6> pre{};
  RelaxOptions ro;
  ro.method = Minimiser::LBFGS;
  ro.ftol = o.ftol;
  ro.etol = 0;
  ro.max_iterations = o.max_iterations;
  int unconverged = 0;
  std::map<std::string, int> stops;
  double asym = 0;
  for (size_t q = 0; q < configs.size(); ++q) {
    const System& s = configs[q];
    if (!s.cell.valid()) throw std::invalid_argument("elastic constants need a periodic cell");
    const ForceField ff = o.field ? *o.field : default_forcefield(s);
    if (ff.atom_type.size() != s.atoms.size()) throw FieldError("the force field does not match the structure");
    Evaluator ev(ff, o.energy);
    std::vector<double> x = positions(s), f;
    const double steps = 13.0 * configs.size();
    const double base = 13.0 * q;
    if (cancelled(o.progress, "minimising configuration " + std::to_string(q + 1), base / steps)) throw std::runtime_error("elastic constants cancelled");
    // minimise, then freeze the pair set there: the energy is then smooth in the strain (no pair crosses the truncated
    // cut-off), so the minimiser can reach tight tolerances and ± strains see the same interactions
    minimise(ev, x, s.cell, ro, "minimise");
    ev.freeze_pairs(x, s.cell);
    auto st = minimise(ev, x, s.cell, ro, "minimise, pairs frozen");
    if (st.stopped_by != "force") { ++unconverged; stops[st.stopped_by]++; }
    const auto s0 = virial_stress(ev.compute(x, s.cell, f), s.cell.volume());
    for (int v = 0; v < 6; ++v) pre[v] += s0[v] * kGPa / configs.size();
    Mat6 C{};
    int k = 0;
    for (int J = 0; J < 6; ++J) {
      std::array<double, 6> sp{}, sm{};
      for (int sg : {1, -1}) {
        if (cancelled(o.progress, "strain " + std::to_string(J + 1) + (sg > 0 ? "+" : "−") + " of configuration " + std::to_string(q + 1),
                      (base + 1 + k++) / steps))
          throw std::runtime_error("elastic constants cancelled");
        M3 F = identity();
        const M3 e = voigt_strain(J, sg * o.strain);
        for (int i = 0; i < 3; ++i)
          for (int j = 0; j < 3; ++j) F[i][j] += e[i][j];
        std::vector<double> y;
        Cell d;
        deform(F, x, s.cell, y, d);
        auto sj = minimise(ev, y, d, ro, "strained");
        if (sj.stopped_by != "force") { ++unconverged; stops[sj.stopped_by]++; }
        const auto sv = virial_stress(ev.compute(y, d, f), d.volume());
        (sg > 0 ? sp : sm) = sv;
      }
      for (int I = 0; I < 6; ++I) C[I][J] = (sp[I] - sm[I]) / (2 * o.strain) * kGPa;
    }
    ev.unfreeze();
    double big = 0, dif = 0;
    for (int I = 0; I < 6; ++I)
      for (int J = 0; J < 6; ++J) { big = std::max(big, std::fabs(C[I][J])); dif = std::max(dif, std::fabs(C[I][J] - C[J][I])); }
    asym = std::max(asym, big > 0 ? dif / big : 0.0);
    for (int I = 0; I < 6; ++I)
      for (int J = I + 1; J < 6; ++J) C[I][J] = C[J][I] = 0.5 * (C[I][J] + C[J][I]);
    all.push_back(C);
  }
  const size_t n = all.size();
  for (int I = 0; I < 6; ++I)
    for (int J = 0; J < 6; ++J) {
      double m = 0;
      for (const auto& c : all) m += c[I][J];
      m /= n;
      double v = 0;
      for (const auto& c : all) v += (c[I][J] - m) * (c[I][J] - m);
      res.C[I][J] = m;
      res.err[I][J] = n > 1 ? std::sqrt(v / (n - 1) / n) : NaN;
    }
  res.prestress = pre;
  res.asymmetry = asym;
  res.configurations = int(n);
  isotropic_averages(res);
  res.method = "static strain (Theodorou & Suter 1986): each configuration minimised at fixed cell (L-BFGS, |F| < " + fmt(o.ftol, 2) +
               " kcal/mol/Å, pair set frozen at the minimum), ±" + fmt(o.strain, 2) + " pure strain in each Voigt component, atoms re-minimised, C_IJ = Δσ_I/Δε_J; " +
               std::to_string(n) + " configuration" + (n == 1 ? "" : "s");
  if (unconverged) {
    std::string why;
    for (const auto& [k, v] : stops) why += (why.empty() ? "" : ", ") + std::to_string(v) + " by " + k;
    res.notes.push_back(std::to_string(unconverged) + " minimisations stopped before the force tolerance (" + why + "); the constants carry that noise");
  }
  if (asym > 0.05)
    res.notes.push_back("C is not symmetric to " + fmt(100 * asym, 2) +
                        " %: at large strains (above ~1e-3) amorphous cells rearrange plastically during the re-minimisation, at very small ones the "
                        "force tolerance limits the precision; try 1e-4 to 1e-3");
  double pmax = 0;
  for (double p : pre) pmax = std::max(pmax, std::fabs(p));
  if (pmax > 0.05 * std::max(1e-9, res.K_hill))
    res.notes.push_back("the minimised cell is under stress (largest component " + fmt(pmax, 3) + " GPa); relax the box first for zero-pressure constants");
  if (n == 1) res.notes.push_back("one configuration: amorphous cells vary; average several (frames of an equilibrated trajectory) for an error estimate");
  return res;
}

// ---------------------------------------------------------------- stress fluctuations

ElasticResult fluctuation_elastic(const Trajectory& t, const std::vector<size_t>& frames, const FluctuationOptions& o) {
  ElasticResult res;
  if (!o.ff) throw std::invalid_argument("the fluctuation method needs the force field of the trajectory");
  if (frames.size() < 2) throw std::invalid_argument("the fluctuation method needs frames of an NVT run");
  if (!t.topology.cell.valid()) throw std::invalid_argument("elastic constants need a periodic cell");
  if (o.temperature <= 0) throw std::invalid_argument("give the temperature of the NVT run");
  Evaluator ev(*o.ff, o.energy);
  const size_t n = frames.size(), natoms = t.topology.atoms.size();
  std::vector<Mat6> born(n);
  std::vector<std::array<double, 6>> sig(n);
  std::vector<double> vol(n);
  std::vector<double> f;
  for (size_t k = 0; k < n; ++k) {
    if (cancelled(o.progress, "Born term and stress", double(k) / n)) throw std::runtime_error("elastic constants cancelled");
    const System s = t.frame(frames[k]);
    const std::vector<double> x = positions(s);
    vol[k] = s.cell.volume();
    sig[k] = virial_stress(ev.compute(x, s.cell, f), vol[k]);
    born[k] = born_matrix(ev, x, s.cell, o.strain);
  }
  const double vmean = std::accumulate(vol.begin(), vol.end(), 0.0) / n;
  const double vmin = *std::min_element(vol.begin(), vol.end()), vmax = *std::max_element(vol.begin(), vol.end());
  if ((vmax - vmin) / vmean > 1e-6) res.notes.push_back("the cell changes between frames: the fluctuation formula holds at constant volume (NVT)");
  const double kT = kB * o.temperature;
  // C over a set of frames
  auto estimate = [&](size_t b, size_t e, Mat6& Cb, Mat6& Fl) {
    const size_t m = e - b;
    std::array<double, 6> mean{};
    for (size_t k = b; k < e; ++k)
      for (int I = 0; I < 6; ++I) mean[I] += sig[k][I] / m;
    for (int I = 0; I < 6; ++I)
      for (int J = 0; J < 6; ++J) {
        double bm = 0, cov = 0;
        for (size_t k = b; k < e; ++k) { bm += born[k][I][J] / m; cov += (sig[k][I] - mean[I]) * (sig[k][J] - mean[J]) / m; }
        Cb[I][J] = bm;
        Fl[I][J] = vmean / kT * cov * kGPa;
      }
  };
  Mat6 Bm, Fm, K{};
  estimate(0, n, Bm, Fm);
  const double kin = natoms * kT / vmean * kGPa;
  for (int I = 0; I < 6; ++I)
    for (int J = 0; J < 6; ++J) K[I][J] = kin * ((I == J ? 1 : 0) + (I < 3 && J < 3 ? 1 : 0));
  for (int I = 0; I < 6; ++I)
    for (int J = 0; J < 6; ++J) res.C[I][J] = Bm[I][J] - Fm[I][J] + K[I][J];
  res.born = Bm;
  res.fluct = Fm;
  res.kinetic = K;
  // block errors
  const size_t nb = std::min<size_t>(std::max(2, o.blocks), n / 2);
  for (auto& row : res.err) row.fill(NaN);
  if (nb >= 2) {
    std::vector<Mat6> cs;
    for (size_t b = 0; b < nb; ++b) {
      Mat6 cb, fb;
      estimate(b * n / nb, (b + 1) * n / nb, cb, fb);
      for (int I = 0; I < 6; ++I)
        for (int J = 0; J < 6; ++J) cb[I][J] = cb[I][J] - fb[I][J] + K[I][J];
      cs.push_back(cb);
    }
    for (int I = 0; I < 6; ++I)
      for (int J = 0; J < 6; ++J) {
        double m = 0, v = 0;
        for (const auto& c : cs) m += c[I][J] / nb;
        for (const auto& c : cs) v += (c[I][J] - m) * (c[I][J] - m);
        res.err[I][J] = std::sqrt(v / (nb - 1) / nb);
      }
  }
  for (int v = 0; v < 6; ++v) {
    double m = 0;
    for (const auto& s : sig) m += s[v] / n;
    res.prestress[v] = m * kGPa - (v < 3 ? natoms * kT / vmean * kGPa : 0.0);   // total stress: virial minus kinetic pressure
  }
  res.configurations = int(n);
  isotropic_averages(res);
  res.method = "stress fluctuations at constant volume (Lutsko 1989): ⟨Born⟩ − (V/kT) cov(σ) + NkT/V over " + std::to_string(n) +
               " frames at " + fmt(o.temperature, 4) + " K; Born term by central differences (±" + fmt(o.strain, 2) +
               " Lagrangian strain) of the virial with the pair set frozen";
  // convergence hint: fluctuation and Born terms nearly cancel in bonded systems
  double ratio = 0;
  for (int I = 0; I < 3; ++I) ratio = std::max(ratio, Bm[I][I] != 0 ? Fm[I][I] / Bm[I][I] : 0.0);
  if (ratio > 0.8)
    res.notes.push_back("the fluctuation term cancels " + fmt(100 * ratio, 3) + " % of the Born term, so the variance must be known to a fraction of " +
                        fmt(100 * (1 - std::min(ratio, 0.999)), 2) + " %: " + std::to_string(n) +
                        " saved frames are far too few for a bonded system (the run protocol samples every step; caps elastic DATA --method fluct-run)");
  return res;
}

ViscosityResult green_kubo_viscosity(const std::vector<std::array<double, 6>>& p, double dt_fs, double volume, double temperature,
                                     double corr_ps, int blocks) {
  ViscosityResult r;
  const size_t n = p.size();
  if (n < 10) throw std::invalid_argument("too few pressure samples for Green–Kubo");
  const size_t L = std::min(n / 2, size_t(std::max(2.0, corr_ps * 1000 / dt_fs)));
  // the traceless symmetric tensor: off-diagonals, and the diagonal minus its trace/3
  std::vector<std::array<double, 6>> q(n);
  for (size_t k = 0; k < n; ++k) {
    const double tr = (p[k][0] + p[k][1] + p[k][2]) / 3;
    q[k] = {p[k][0] - tr, p[k][1] - tr, p[k][2] - tr, p[k][3], p[k][4], p[k][5]};
  }
  // Σ_αβ P°αβ(0) P°αβ(t) over all nine elements: the diagonal once, each off-diagonal twice
  auto acf_of = [&](size_t from, size_t to) {
    std::vector<double> c(L, 0.0);
    for (size_t lag = 0; lag < L; ++lag) {
      double sum = 0;
      size_t m = 0;
      for (size_t k = from; k + lag < to; ++k, ++m) {
        const auto& a = q[k];
        const auto& b = q[k + lag];
        sum += a[0] * b[0] + a[1] * b[1] + a[2] * b[2] + 2 * (a[3] * b[3] + a[4] * b[4] + a[5] * b[5]);
      }
      c[lag] = m ? sum / double(m) : 0;
    }
    return c;
  };
  // atm² · Å³ · fs / (J) → Pa·s: (101325 Pa)² · 1e-30 m³ · 1e-15 s / (kB T); ×1000 for mPa·s
  const double kT = 1.380649e-23 * temperature;
  const double to_mpas = 101325.0 * 101325.0 * volume * 1e-30 * dt_fs * 1e-15 / kT / 10.0 * 1000.0;
  auto integrate = [&](const std::vector<double>& c) {
    std::vector<double> run(L, 0.0);
    for (size_t k = 1; k < L; ++k) run[k] = run[k - 1] + 0.5 * (c[k - 1] + c[k]) * to_mpas;   // trapezoid
    return run;
  };
  auto plateau = [&](const std::vector<double>& run) {
    double m = 0;
    size_t c = 0;
    for (size_t k = 2 * L / 3; k < L; ++k) m += run[k], ++c;
    return c ? m / double(c) : run.back();
  };
  const auto c = acf_of(0, n);
  r.running = integrate(c);
  r.eta = plateau(r.running);
  for (size_t k = 0; k < L; ++k) r.t_ps.push_back(double(k) * dt_fs / 1000), r.acf.push_back(c[0] != 0 ? c[k] / c[0] : 0);
  // blocks of the run, each integrated alone: the spread of their plateaus
  const int nb = std::max(2, blocks);
  std::vector<double> e;
  for (int b = 0; b < nb; ++b) {
    const size_t from = n * size_t(b) / size_t(nb), to = n * size_t(b + 1) / size_t(nb);
    if (to - from < 2 * L) continue;
    e.push_back(plateau(integrate(acf_of(from, to))));
  }
  if (e.size() >= 2) {
    double m = 0, v = 0;
    for (double x : e) m += x;
    m /= double(e.size());
    for (double x : e) v += (x - m) * (x - m);
    r.error = std::sqrt(v / double(e.size() - 1) / double(e.size()));
  } else {
    r.error = std::numeric_limits<double>::quiet_NaN();
    r.notes.push_back("the run is shorter than " + std::to_string(2 * nb) + " correlation windows: no block error");
  }
  // levelled off: the last third varies by less than a fifth of its mean
  double lo = 1e300, hi = -1e300;
  for (size_t k = 2 * L / 3; k < L; ++k) lo = std::min(lo, r.running[k]), hi = std::max(hi, r.running[k]);
  r.plateau = std::fabs(r.eta) > 0 && (hi - lo) < 0.2 * std::fabs(r.eta);
  if (!r.plateau)
    r.notes.push_back("the running integral has not levelled off within " + fmt(corr_ps, 3) +
                      " ps: the stress has not relaxed (a polymer melt needs far longer runs and correlation windows); the value is a lower bound at best");
  return r;
}

ViscosityResult viscosity_green_kubo(System& s, const ViscosityOptions& o) {
  if (!s.cell.valid()) throw std::invalid_argument("viscosity needs a periodic cell");
  const ForceField ff = o.field ? *o.field : default_forcefield(s);
  DynamicsOptions d;
  d.field = std::make_shared<const ForceField>(ff);
  d.energy = o.energy;
  d.dt = o.dt;
  d.temperature = o.temperature;
  d.thermostat = Thermostat::NoseHoover;   // deterministic and weakly coupled: the dynamics the correlations need
  d.tau_t = o.tau_t;
  d.seed = o.seed;
  d.frame_every = 0;
  d.new_velocities = o.new_velocities || s.velocities.size() != s.atoms.size();
  const int64_t neq = std::llround(o.equilibrate_ps * 1000 / o.dt), nrun = std::llround(o.ps * 1000 / o.dt);
  if (neq > 0) {
    d.steps = neq;
    d.thermo_every = 1000;
    d.progress = [&](const ThermoRow& r) { return !cancelled(o.progress, "equilibrating at " + fmt(o.temperature, 4) + " K", double(r.step) / double(neq + nrun)); };
    run_dynamics(s, d);
    d.new_velocities = false;
  }
  std::vector<std::array<double, 6>> p;
  d.steps = nrun;
  d.step_offset = neq;
  d.thermo_every = std::max(1, o.sample_every);
  d.progress = [&](const ThermoRow& r) {
    if (r.step > neq) p.push_back({r.p[0], r.p[1], r.p[2], r.p[3], r.p[4], r.p[5]});
    return !cancelled(o.progress, "sampling the pressure tensor", double(r.step) / double(neq + nrun));
  };
  DynamicsReport rep;
  run_dynamics(s, d, &rep);
  ViscosityResult r = green_kubo_viscosity(p, o.dt * std::max(1, o.sample_every), s.cell.volume(), o.temperature, o.corr_ps, o.blocks);
  r.notes.insert(r.notes.begin(), fmt(o.ps, 4) + " ps NVT (Nosé–Hoover, τ " + fmt(o.tau_t, 4) + " fs) after " + fmt(o.equilibrate_ps, 4) +
                                  " ps; the pressure tensor every " + fmt(o.dt * std::max(1, o.sample_every), 3) + " fs (" + std::to_string(p.size()) + " samples)");
  return r;
}

std::vector<Property> viscosity_properties(const ViscosityResult& r) {
  Property q;
  q.id = "viscosity";
  q.name = "Shear viscosity (Green–Kubo)";
  q.unit = "mPa·s";
  q.value = r.eta;
  q.error = r.error;
  q.method = "Green–Kubo integral of the traceless pressure tensor's autocorrelation (Daivis & Evans 1994), the plateau of the running integral";
  q.notes = r.notes;
  Series run{"running integral", "t (ps)", "η (mPa·s)", r.t_ps, r.running};
  Series acf{"stress autocorrelation", "t (ps)", "⟨P°P°⟩(t)/⟨P°P°⟩(0)", r.t_ps, r.acf};
  q.series = {run, acf};
  return {q};
}

ElasticResult fluctuation_run(System& s, const FluctuationRunOptions& o) {
  ElasticResult res;
  if (!s.cell.valid()) throw std::invalid_argument("elastic constants need a periodic cell");
  const ForceField ff = o.field ? *o.field : default_forcefield(s);
  auto shared = std::make_shared<const ForceField>(ff);
  DynamicsOptions d;
  d.field = shared;
  d.energy = o.energy;
  d.dt = o.dt;
  d.temperature = o.temperature;
  d.thermostat = o.thermostat;
  d.tau_t = o.tau_t;
  d.seed = o.seed;
  d.frame_every = 0;
  d.thermo_every = 1000;
  d.new_velocities = o.new_velocities || s.velocities.size() != s.atoms.size();
  // equilibration at the temperature, not sampled
  const int64_t neq = static_cast<int64_t>(std::llround(o.equilibrate_ps * 1000 / o.dt));
  const int64_t nrun = static_cast<int64_t>(std::llround(o.ps * 1000 / o.dt));
  if (neq > 0) {
    d.steps = neq;
    d.progress = [&](const ThermoRow& r) { return !cancelled(o.progress, "equilibrating at " + fmt(o.temperature, 4) + " K", double(r.step) / (neq + nrun)); };
    run_dynamics(s, d);
    d.new_velocities = false;
  }
  // sampling: per block, Σσ, Σσσ and the Born sum
  const int nb = std::max(2, o.blocks);
  struct Block { double n = 0, nb = 0; std::array<double, 6> s{}; Mat6 ss{}, born{}; };
  std::vector<Block> blocks(nb);
  Evaluator ev(ff, o.energy);
  const double V = s.cell.volume();
  int64_t k = 0;
  d.steps = nrun;
  d.step_offset = neq;
  d.progress = [&](const ThermoRow& r) { return !cancelled(o.progress, "sampling stress fluctuations", double(r.step) / (neq + nrun)); };
  d.each_step = [&](const EnergyTerms& e, const std::vector<double>& x, const Cell& c, int64_t) {
    Block& b = blocks[std::min<int64_t>(nb - 1, k * nb / std::max<int64_t>(1, nrun))];
    const auto sg = virial_stress(e, c.volume());
    b.n += 1;
    for (int I = 0; I < 6; ++I) {
      b.s[I] += sg[I];
      for (int J = 0; J < 6; ++J) b.ss[I][J] += sg[I] * sg[J];
    }
    if (k % std::max(1, o.born_every) == 0) {
      const Mat6 B = born_matrix(ev, x, c, o.strain);
      b.nb += 1;
      for (int I = 0; I < 6; ++I)
        for (int J = 0; J < 6; ++J) b.born[I][J] += B[I][J];
    }
    ++k;
  };
  run_dynamics(s, d);
  const double kT = kB * o.temperature;
  const double kin = s.atoms.size() * kT / V * kGPa;
  Mat6 K{};
  for (int I = 0; I < 6; ++I)
    for (int J = 0; J < 6; ++J) K[I][J] = kin * ((I == J ? 1 : 0) + (I < 3 && J < 3 ? 1 : 0));
  auto estimate = [&](const std::vector<const Block*>& bs, Mat6& Bm, Mat6& Fm) {
    double n = 0, nbn = 0;
    std::array<double, 6> m{};
    Mat6 ss{}, bo{};
    for (const Block* b : bs) {
      n += b->n;
      nbn += b->nb;
      for (int I = 0; I < 6; ++I) {
        m[I] += b->s[I];
        for (int J = 0; J < 6; ++J) { ss[I][J] += b->ss[I][J]; bo[I][J] += b->born[I][J]; }
      }
    }
    for (int I = 0; I < 6; ++I) m[I] /= n;
    for (int I = 0; I < 6; ++I)
      for (int J = 0; J < 6; ++J) {
        Bm[I][J] = bo[I][J] / std::max(1.0, nbn);
        Fm[I][J] = V / kT * (ss[I][J] / n - m[I] * m[J]) * kGPa;
      }
  };
  std::vector<const Block*> all;
  for (const auto& b : blocks) all.push_back(&b);
  Mat6 Bm, Fm;
  estimate(all, Bm, Fm);
  for (int I = 0; I < 6; ++I)
    for (int J = 0; J < 6; ++J) res.C[I][J] = Bm[I][J] - Fm[I][J] + K[I][J];
  res.born = Bm;
  res.fluct = Fm;
  res.kinetic = K;
  std::vector<Mat6> cs;
  for (const auto& b : blocks) {
    Mat6 bb, ff2;
    estimate({&b}, bb, ff2);
    Mat6 c;
    for (int I = 0; I < 6; ++I)
      for (int J = 0; J < 6; ++J) c[I][J] = bb[I][J] - ff2[I][J] + K[I][J];
    cs.push_back(c);
  }
  for (int I = 0; I < 6; ++I)
    for (int J = 0; J < 6; ++J) {
      double m = 0, v = 0;
      for (const auto& c : cs) m += c[I][J] / nb;
      for (const auto& c : cs) v += (c[I][J] - m) * (c[I][J] - m);
      res.err[I][J] = std::sqrt(v / (nb - 1) / nb);
    }
  double n = 0;
  for (const auto& b : blocks) n += b.n;
  for (int v = 0; v < 6; ++v) {
    double m = 0;
    for (const auto& b : blocks) m += b.s[v];
    res.prestress[v] = m / n * kGPa - (v < 3 ? kin : 0.0);
  }
  res.configurations = static_cast<int>(n);
  isotropic_averages(res);
  res.method = "stress fluctuations on the fly (Lutsko 1989; Clavier et al. 2017): NVT at " + fmt(o.temperature, 4) + " K (" + to_string(o.thermostat) +
               "), " + fmt(o.equilibrate_ps, 4) + " ps equilibration, " + fmt(o.ps, 5) + " ps sampled with the virial stress at every step and the Born term every " +
               std::to_string(o.born_every) + " steps; errors from " + std::to_string(nb) + " time blocks";
  double ratio = 0;
  for (int I = 0; I < 3; ++I) ratio = std::max(ratio, Bm[I][I] != 0 ? Fm[I][I] / Bm[I][I] : 0.0);
  if (ratio > 0.8)
    res.notes.push_back("the fluctuation term cancels " + fmt(100 * ratio, 3) + " % of the Born term; the constants converge only as the run grows (check the errors)");
  // run length for a 0.1 GPa error on the diagonal (errors fall as 1/√t)
  const double e11 = (res.err[0][0] + res.err[1][1] + res.err[2][2]) / 3;
  if (std::isfinite(e11) && e11 > 0.1)
    res.notes.push_back("block error of C11–C33 " + fmt(e11, 2) + " GPa after " + fmt(o.ps, 4) + " ps; about " + fmt(o.ps * (e11 / 0.1) * (e11 / 0.1) / 1000, 2) +
                        " ns would bring it to 0.1 GPa. For bonded polymers the static method is far cheaper");
  return res;
}

// ---------------------------------------------------------------- uniaxial deformation

void analyse_tensile(TensileResult& r, double fit_strain, bool lateral) {
  const auto& c = r.curve;
  r.smooth.assign(c.size(), 0.0);
  for (size_t i = 0; i < c.size(); ++i) {
    double s = 0;
    int k = 0;
    for (size_t j = 0; j < c.size(); ++j)
      if (std::fabs(c[j].strain - c[i].strain) <= 0.005) { s += c[j].stress; ++k; }
    r.smooth[i] = k ? s / k : c[i].stress;
  }
  std::vector<double> e, s, l;
  for (const auto& p : c)
    if (p.strain <= fit_strain) { e.push_back(p.strain); s.push_back(p.stress); l.push_back(0.5 * (p.lateral1 + p.lateral2)); }
  double a, b, be;
  linfit(e, s, a, b, be);
  r.modulus = b / 1000;
  r.modulus_err = std::isnan(be) ? NaN : be / 1000;
  if (lateral) {
    linfit(e, l, a, b, be);
    r.poisson = -b;
    r.poisson_err = be;
  } else {
    r.poisson = r.poisson_err = NaN;
  }
  r.peak_stress = r.peak_strain = 0;
  for (size_t i = 0; i < c.size(); ++i)
    if (r.smooth[i] > r.peak_stress) { r.peak_stress = r.smooth[i]; r.peak_strain = c[i].strain; }
  r.yield_stress = r.yield_strain = 0;
  const double E = r.modulus * 1000;
  if (E > 0)
    for (size_t i = 0; i < c.size(); ++i)
      if (c[i].strain > std::max(0.004, fit_strain) && r.smooth[i] <= E * (c[i].strain - 0.002)) {
        r.yield_stress = r.smooth[i];
        r.yield_strain = c[i].strain;
        break;
      }
}

namespace {
// mean and sample standard deviation of the finite values (NaN when none)
std::pair<double, double> mean_sd(const std::vector<double>& v) {
  double sum = 0, n = 0;
  for (double x : v) if (std::isfinite(x)) sum += x, ++n;
  if (n == 0) return {std::nan(""), std::nan("")};
  const double m = sum / n;
  double ss = 0;
  for (double x : v) if (std::isfinite(x)) ss += (x - m) * (x - m);
  return {m, n > 1 ? std::sqrt(ss / (n - 1)) : std::nan("")};
}
}  // namespace

// The same test along x, y and z from the same start (an amorphous cell is isotropic only on average): the curves
// averaged point by point (the same rate gives the same strains), E, ν, yield and peak as mean ± s.d. over the axes
TensileResult run_tensile_averaged(System& s, const TensileOptions& o) {
  std::array<TensileResult, 3> r;
  System last;
  for (int k = 0; k < 3; ++k) {
    System c = s;
    TensileOptions ok = o;
    ok.axis = k;
    ok.seed = o.seed + uint64_t(k);
    r[size_t(k)] = run_tensile(c, ok);
    if (k == 2) last = std::move(c);
  }
  TensileResult res;
  size_t n = std::min({r[0].curve.size(), r[1].curve.size(), r[2].curve.size()});
  double worst = 0;
  for (size_t i = 0; i < n; ++i) {
    TensilePoint p;
    for (const auto& x : r) {
      p.strain += x.curve[i].strain / 3, p.stress += x.curve[i].stress / 3;
      p.lateral1 += x.curve[i].lateral1 / 3, p.lateral2 += x.curve[i].lateral2 / 3;
      p.temperature += x.curve[i].temperature / 3, p.time_ps += x.curve[i].time_ps / 3;
      worst = std::max(worst, std::fabs(x.curve[i].strain - r[0].curve[i].strain));
    }
    res.curve.push_back(p);
    double sm = 0;
    for (const auto& x : r) sm += (i < x.smooth.size() ? x.smooth[i] : x.curve[i].stress) / 3;
    res.smooth.push_back(sm);
  }
  if (worst > 1e-9) res.notes.push_back("the three curves' strains differ by up to " + std::to_string(worst) + ": averaged sample by sample");
  auto field = [&](auto get) { std::vector<double> v; for (const auto& x : r) v.push_back(get(x)); return mean_sd(v); };
  std::tie(res.modulus, res.modulus_err) = field([](const TensileResult& x) { return x.modulus; });
  std::tie(res.poisson, res.poisson_err) = field([](const TensileResult& x) { return x.poisson; });
  res.yield_stress = field([](const TensileResult& x) { return x.yield_strain > 0 ? x.yield_stress : std::nan(""); }).first;
  res.yield_strain = field([](const TensileResult& x) { return x.yield_strain > 0 ? x.yield_strain : std::nan(""); }).first;
  if (!std::isfinite(res.yield_strain)) res.yield_stress = 0, res.yield_strain = 0;
  res.peak_stress = field([](const TensileResult& x) { return x.peak_stress; }).first;
  res.peak_strain = field([](const TensileResult& x) { return x.peak_strain; }).first;
  res.method = r[0].method + "; along x, y and z from the same start, averaged (± s.d. over the three directions)";
  char b[256];
  std::snprintf(b, sizeof b, "E by axis: x %.3f, y %.3f, z %.3f GPa; ν: x %.3f, y %.3f, z %.3f", r[0].modulus, r[1].modulus, r[2].modulus,
                r[0].poisson, r[1].poisson, r[2].poisson);
  res.notes.push_back(b);
  for (int k = 0; k < 3; ++k)
    for (const auto& note : r[size_t(k)].notes) res.notes.push_back(std::string(1, "xyz"[k]) + ": " + note);
  res.notes.push_back("the structure kept is the one pulled along z");
  s = std::move(last);
  return res;
}

TensileResult run_tensile(System& s, const TensileOptions& o) {
  if (o.axis == 3) return run_tensile_averaged(s, o);
  TensileResult res;
  if (!s.cell.valid()) throw std::invalid_argument("a tensile test needs a periodic cell");
  if (o.axis < 0 || o.axis > 2) throw std::invalid_argument("axis must be 0 (x), 1 (y), 2 (z) or 3 (averaged over the three)");
  if (o.rate <= 0 || o.max_strain <= 0) throw std::invalid_argument("strain rate and final strain must be positive");
  DynamicsOptions d;
  d.field = o.field;
  d.energy = o.energy;
  d.dt = o.dt;
  d.steps = static_cast<int64_t>(std::ceil(o.max_strain / o.rate / (o.dt * 1e-3)));
  d.temperature = o.temperature;
  d.thermostat = o.thermostat;
  d.tau_t = o.tau_t;
  d.barostat = o.lateral_pressure ? Barostat::Berendsen : Barostat::None;
  d.anisotropic = true;
  for (int k = 0; k < 3; ++k) d.couple_axis[k] = k != o.axis;
  d.pressure = o.pressure;
  d.tau_p = o.tau_p;
  d.compressibility = o.compressibility;
  d.barostat_every = 10;
  d.deform_axis = o.axis;
  d.deform_rate = o.rate;
  d.thermo_every = std::max(1, o.sample_every);
  d.frame_every = o.frame_every;
  d.frame = o.frame;
  d.seed = o.seed;
  d.new_velocities = o.new_velocities;
  if (o.equilibrate_ps > 0) {
    DynamicsOptions e = d;
    e.deform_axis = -1;
    e.barostat = Barostat::Berendsen;
    for (int k = 0; k < 3; ++k) e.couple_axis[k] = true;
    e.steps = static_cast<int64_t>(std::llround(o.equilibrate_ps * 1000 / o.dt));
    e.frame_every = 0;
    e.frame = nullptr;
    e.progress = nullptr;
    run_dynamics(s, e);
    d.new_velocities = false;
  }
  double L0[3] = {0, 0, 0};
  bool first = true;
  d.progress = [&](const ThermoRow& r) {
    const double L[3] = {r.lx, r.ly, r.lz};
    if (first) { for (int k = 0; k < 3; ++k) L0[k] = L[k]; first = false; }
    TensilePoint p;
    p.strain = L[o.axis] / L0[o.axis] - 1;
    p.stress = -r.p[o.axis] * kAtmMPa;
    const int a1 = (o.axis + 1) % 3, a2 = (o.axis + 2) % 3;
    p.lateral1 = L[a1] / L0[a1] - 1;
    p.lateral2 = L[a2] / L0[a2] - 1;
    p.temperature = r.temperature;
    p.time_ps = r.time_ps;
    res.curve.push_back(p);
    return !o.progress || o.progress(p);
  };
  DynamicsReport rep;
  run_dynamics(s, d, &rep);
  analyse_tensile(res, o.fit_strain, o.lateral_pressure);
  const char* ax = o.axis == 0 ? "x" : o.axis == 1 ? "y" : "z";
  res.method = std::string("uniaxial deformation along ") + ax + " at " + fmt(o.rate * 1e12, 3) + " s⁻¹ (engineering) to " +
               fmt(100 * o.max_strain, 3) + " %, " + fmt(o.temperature, 4) + " K (" + to_string(o.thermostat) + "), lateral axes " +
               (o.lateral_pressure ? "at " + fmt(o.pressure, 3) + " atm (Berendsen per axis)" : std::string("fixed")) +
               "; " + fmt(o.equilibrate_ps, 4) + " ps NPT first (each axis at the target pressure); stress −P_" + ax + ax + "; modulus and Poisson ratio fitted from 0 to " + fmt(100 * o.fit_strain, 2) + " % strain";
  // the run log (steps, speed, typing) is not a warning; only the first line (speed) goes into the method
  if (!rep.notes.empty()) res.method += "; " + rep.notes.front();
  if (o.rate > 5e-3) res.notes.push_back("strain rates above ~10⁹ s⁻¹ raise the modulus and yield stress of glassy polymers; experiments are 10⁻³–10⁰ s⁻¹");
  if (!res.curve.empty() && res.yield_strain == 0) res.notes.push_back("no yield point (0.2 % offset) within the deformation");
  return res;
}

// ---------------------------------------------------------------- glass transition

BilinearFit fit_bilinear(const std::vector<double>& T0, const std::vector<double>& y0) {
  BilinearFit fit;
  const size_t n = T0.size();
  if (n != y0.size() || n < 5) { fit.note = "a two-line fit needs at least five temperatures"; return fit; }
  std::vector<size_t> ord(n);
  std::iota(ord.begin(), ord.end(), 0);
  std::sort(ord.begin(), ord.end(), [&](size_t a, size_t b) { return T0[a] < T0[b]; });
  std::vector<double> T(n), Y(n);
  for (size_t i = 0; i < n; ++i) { T[i] = T0[ord[i]]; Y[i] = y0[ord[i]]; }

  // least squares for a fixed hinge: y = c0 + c1 min(T − h, 0) + c2 max(T − h, 0)
  auto solve = [&](double h, const std::vector<double>& y, double c[3]) {
    double A[3][4] = {};
    for (size_t i = 0; i < n; ++i) {
      const double u[3] = {1, std::min(T[i] - h, 0.0), std::max(T[i] - h, 0.0)};
      for (int a = 0; a < 3; ++a) {
        for (int b = 0; b < 3; ++b) A[a][b] += u[a] * u[b];
        A[a][3] += u[a] * y[i];
      }
    }
    for (int k = 0; k < 3; ++k) {
      int p = k;
      for (int i = k + 1; i < 3; ++i) if (std::fabs(A[i][k]) > std::fabs(A[p][k])) p = i;
      for (int j = 0; j < 4; ++j) std::swap(A[k][j], A[p][j]);
      if (std::fabs(A[k][k]) < 1e-300) return std::numeric_limits<double>::max();
      for (int i = 0; i < 3; ++i)
        if (i != k) {
          const double m = A[i][k] / A[k][k];
          for (int j = 0; j < 4; ++j) A[i][j] -= m * A[k][j];
        }
    }
    for (int k = 0; k < 3; ++k) c[k] = A[k][3] / A[k][k];
    double rss = 0;
    for (size_t i = 0; i < n; ++i) {
      const double r = y[i] - (c[0] + c[1] * std::min(T[i] - h, 0.0) + c[2] * std::max(T[i] - h, 0.0));
      rss += r * r;
    }
    return rss;
  };
  // hinge between the second and the second-to-last temperature: a scan, then golden-section refinement
  auto best_hinge = [&](const std::vector<double>& y, int scan) {
    const double lo = T[1], hi = T[n - 2];
    double c[3], bh = lo, br = std::numeric_limits<double>::max();
    for (int k = 0; k <= scan; ++k) {
      const double h = lo + (hi - lo) * k / scan, r = solve(h, y, c);
      if (r < br) { br = r; bh = h; }
    }
    double a = std::max(lo, bh - (hi - lo) / scan), b = std::min(hi, bh + (hi - lo) / scan);
    const double g = 0.5 * (std::sqrt(5.0) - 1);
    double x1 = b - g * (b - a), x2 = a + g * (b - a), f1 = solve(x1, y, c), f2 = solve(x2, y, c);
    for (int it = 0; it < 60; ++it) {
      if (f1 < f2) { b = x2; x2 = x1; f2 = f1; x1 = b - g * (b - a); f1 = solve(x1, y, c); }
      else { a = x1; x1 = x2; f1 = f2; x2 = a + g * (b - a); f2 = solve(x2, y, c); }
    }
    const double h = 0.5 * (a + b);
    return solve(h, y, c) <= br ? h : bh;
  };
  const double h = best_hinge(Y, 400);
  double c[3];
  fit.rss = solve(h, Y, c);
  fit.tg = h;
  fit.value_at_tg = c[0];
  fit.slope_low = c[1];
  fit.slope_high = c[2];
  fit.alpha_low = c[1] / c[0];
  fit.alpha_high = c[2] / c[0];
  fit.ok = true;
  // bootstrap of residuals
  std::vector<double> fitted(n), res(n);
  for (size_t i = 0; i < n; ++i) {
    fitted[i] = c[0] + c[1] * std::min(T[i] - h, 0.0) + c[2] * std::max(T[i] - h, 0.0);
    res[i] = Y[i] - fitted[i];
  }
  std::mt19937_64 rng(12345);
  std::uniform_int_distribution<size_t> pick(0, n - 1);
  std::vector<double> hs;
  for (int b = 0; b < 200; ++b) {
    std::vector<double> yb(n);
    for (size_t i = 0; i < n; ++i) yb[i] = fitted[i] + res[pick(rng)];
    hs.push_back(best_hinge(yb, 100));
  }
  double m = 0, v = 0;
  for (double x : hs) m += x / hs.size();
  for (double x : hs) v += (x - m) * (x - m) / (hs.size() - 1);
  fit.tg_err = std::sqrt(v);
  if (std::fabs(c[2] - c[1]) < 0.1 * std::max(std::fabs(c[1]), std::fabs(c[2])))
    fit.note = "the two slopes differ by less than 10 %: no clear transition in this range";
  else if (h - T[0] < 1.5 * (T[1] - T[0]) || T[n - 1] - h < 1.5 * (T[n - 1] - T[n - 2]))
    fit.note = "the hinge sits at the edge of the range: extend the temperatures";
  return fit;
}

BilinearFit fit_two_ranges(const std::vector<double>& T, const std::vector<double>& y, double glassy_max, double rubbery_min) {
  BilinearFit fit;
  const size_t n = T.size();
  if (n != y.size() || n < 4) { fit.note = "a two-range fit needs at least four temperatures"; return fit; }
  const auto [tmin, tmax] = std::minmax_element(T.begin(), T.end());
  if (glassy_max <= 0) glassy_max = *tmin + (*tmax - *tmin) / 3.0;
  if (rubbery_min <= 0) rubbery_min = *tmax - (*tmax - *tmin) / 3.0;
  if (glassy_max >= rubbery_min) { fit.note = "the glassy range must end below where the rubbery range starts"; return fit; }
  std::vector<size_t> lo, hi;
  for (size_t i = 0; i < n; ++i) {
    if (T[i] <= glassy_max + 1e-9) lo.push_back(i);
    if (T[i] >= rubbery_min - 1e-9) hi.push_back(i);
  }
  if (lo.size() < 2 || hi.size() < 2) { fit.note = "each range needs at least two temperatures"; return fit; }
  auto line = [&](const std::vector<size_t>& idx, const std::vector<double>& Y, double& a, double& b) {
    double st = 0, sy = 0, stt = 0, sty = 0;
    for (size_t i : idx) st += T[i], sy += Y[i], stt += T[i] * T[i], sty += T[i] * Y[i];
    const double m = double(idx.size()), den = m * stt - st * st;
    b = den != 0 ? (m * sty - st * sy) / den : 0;
    a = (sy - b * st) / m;
  };
  auto cross_at = [&](const std::vector<double>& Y, double& tg, double& bl, double& bh, double& al, double& ah) {
    line(lo, Y, al, bl);
    line(hi, Y, ah, bh);
    tg = bl != bh ? (ah - al) / (bl - bh) : 0;
  };
  double al, bl, ah, bh, tg;
  cross_at(y, tg, bl, bh, al, ah);
  if (bl == bh) { fit.note = "the two lines are parallel: no crossing"; return fit; }
  fit.ok = true;
  fit.tg = tg;
  fit.slope_low = bl;
  fit.slope_high = bh;
  fit.value_at_tg = al + bl * tg;
  if (fit.value_at_tg != 0) fit.alpha_low = bl / fit.value_at_tg, fit.alpha_high = bh / fit.value_at_tg;
  // residuals of the points in the ranges, resampled onto the fitted lines
  std::vector<double> res(n, 0), fitted(n, 0);
  for (size_t i : lo) fitted[i] = al + bl * T[i], res[i] = y[i] - fitted[i];
  for (size_t i : hi) fitted[i] = ah + bh * T[i], res[i] = y[i] - fitted[i];
  for (size_t i : lo) fit.rss += res[i] * res[i];
  for (size_t i : hi) fit.rss += res[i] * res[i];
  std::vector<size_t> used(lo);
  used.insert(used.end(), hi.begin(), hi.end());
  std::mt19937_64 rng(12345);
  std::uniform_int_distribution<size_t> pick(0, used.size() - 1);
  std::vector<double> tgs;
  for (int b = 0; b < 400; ++b) {
    std::vector<double> yb = y;
    for (size_t i : used) yb[i] = fitted[i] + res[used[pick(rng)]];
    double t2, b1, b2, a1, a2;
    cross_at(yb, t2, b1, b2, a1, a2);
    if (b1 != b2 && std::isfinite(t2)) tgs.push_back(t2);
  }
  if (tgs.size() > 10) {
    double m = std::accumulate(tgs.begin(), tgs.end(), 0.0) / double(tgs.size()), v = 0;
    for (double x : tgs) v += (x - m) * (x - m);
    fit.tg_err = std::sqrt(v / double(tgs.size() - 1));
  }
  if (lo.size() < 3 || hi.size() < 3) {   // two points fix a line exactly: nothing for the bootstrap to resample
    fit.tg_err = std::numeric_limits<double>::quiet_NaN();
    fit.note = "a range with two temperatures fixes its line exactly: no error estimate — use at least three per range";
  }
  if (tg < glassy_max || tg > rubbery_min) {
    // lines that cross outside the gap between the ranges (nearly parallel, or noise) give no transition
    fit.ok = false;
    fit.note = "the lines cross at " + fmt(tg, 4) + " K, outside the gap between the ranges (" + fmt(glassy_max, 4) + "–" + fmt(rubbery_min, 4) +
               " K): no Tg — longer holds, a wider scan or other ranges";
  }
  return fit;
}

CoolingResult run_cooling(System& s, const CoolingOptions& o) {
  CoolingResult res;
  if (!s.cell.valid()) throw std::invalid_argument("a cooling run needs a periodic cell");
  if (o.t_step <= 0 || o.ps_per_step <= 0) throw std::invalid_argument("temperature step and hold time must be positive");
  std::vector<double> temps;
  const double dir = o.t_end < o.t_start ? -1 : 1;
  for (double T = o.t_start; dir * (T - o.t_end) <= 1e-9; T += dir * o.t_step) temps.push_back(T);
  int64_t offset = 0;
  const int64_t steps = static_cast<int64_t>(std::llround(o.ps_per_step * 1000 / o.dt));
  bool first_velocities = o.new_velocities;
  // equilibrate at the starting temperature before the first sampled hold
  const double eq_ps = o.equilibrate_ps < 0 ? o.ps_per_step : o.equilibrate_ps;
  if (eq_ps > 0) {
    DynamicsOptions d;
    d.field = o.field;
    d.energy = o.energy;
    d.dt = o.dt;
    d.steps = static_cast<int64_t>(std::llround(eq_ps * 1000 / o.dt));
    d.temperature = o.t_start;
    d.thermostat = Thermostat::Bussi;
    d.tau_t = o.tau_t;
    d.barostat = o.barostat;
    d.pressure = o.pressure;
    d.tau_p = o.tau_p;
    d.compressibility = o.compressibility;
    d.seed = o.seed + 1000;
    d.new_velocities = first_velocities;
    d.thermo_every = std::max<int64_t>(1, d.steps / 50);
    d.frame_every = 0;
    if (o.progress) d.progress = [&](const ThermoRow& r) { return o.progress(r, -1, int(temps.size())); };
    run_dynamics(s, d);
    offset += d.steps;
    first_velocities = false;
  }
  for (size_t k = 0; k < temps.size(); ++k) {
    DynamicsOptions d;
    d.field = o.field;
    d.energy = o.energy;
    d.dt = o.dt;
    d.steps = steps;
    d.temperature = temps[k];
    d.thermostat = Thermostat::Bussi;
    d.tau_t = o.tau_t;
    d.barostat = o.barostat;
    d.pressure = o.pressure;
    d.tau_p = o.tau_p;
    d.compressibility = o.compressibility;
    d.seed = o.seed + k;
    d.new_velocities = k == 0 && first_velocities;
    d.thermo_every = std::max<int64_t>(1, steps / 200);
    d.frame_every = 0;
    d.step_offset = offset;
    const int index = int(k), total = int(temps.size());
    if (o.progress) d.progress = [&, index, total](const ThermoRow& r) { return o.progress(r, index, total); };
    DynamicsReport rep;
    run_dynamics(s, d, &rep);
    offset += steps;
    // average over the second part of the hold
    std::vector<double> dens, pot;
    for (size_t i = 0; i < rep.thermo.size(); ++i)
      if (i >= size_t(o.average_from * rep.thermo.size())) { dens.push_back(rep.thermo[i].density); pot.push_back(rep.thermo[i].potential); }
    for (const auto& r : rep.thermo) res.thermo.push_back(r);
    CoolingPoint p;
    p.temperature = temps[k];
    p.density = std::accumulate(dens.begin(), dens.end(), 0.0) / std::max<size_t>(1, dens.size());
    p.potential = std::accumulate(pot.begin(), pot.end(), 0.0) / std::max<size_t>(1, pot.size());
    // standard error from 5 blocks
    const size_t nb = std::min<size_t>(5, dens.size());
    if (nb >= 2) {
      std::vector<double> bm;
      for (size_t b = 0; b < nb; ++b) {
        double sum = 0;
        size_t c = 0;
        for (size_t i = b * dens.size() / nb; i < (b + 1) * dens.size() / nb; ++i) { sum += dens[i]; ++c; }
        bm.push_back(sum / std::max<size_t>(1, c));
      }
      double v = 0;
      for (double x : bm) v += (x - p.density) * (x - p.density);
      p.density_err = std::sqrt(v / (nb - 1) / nb);
    }
    p.specific_volume = p.density > 0 ? 1 / p.density : 0;
    res.points.push_back(p);
  }
  std::vector<double> T, v;
  const double atoms = double(std::max<size_t>(1, s.atoms.size()));
  for (const auto& p : res.points) { T.push_back(p.temperature); v.push_back(o.property == 1 ? p.potential / atoms : p.specific_volume); }
  res.property = o.property;
  res.fitted = v;
  res.fit = o.fit == 1 ? fit_two_ranges(T, v, o.glassy_max, o.rubbery_min) : fit_bilinear(T, v);
  const std::string what = o.property == 1 ? "potential energy per atom" : "specific volume";
  const std::string how = o.fit == 1 ? "the crossing of straight lines through the glassy and the rubbery range" : "the hinge of a continuous two-line fit";
  const double rate = o.t_step / o.ps_per_step;   // K/ps
  res.method = "stepwise cooling " + fmt(o.t_start, 4) + " → " + fmt(o.t_end, 4) + " K in " + fmt(o.t_step, 3) + " K steps of " +
               fmt(o.ps_per_step, 4) + " ps after " + fmt(eq_ps, 4) + " ps at " + fmt(o.t_start, 4) + " K (NPT, " + fmt(o.pressure, 3) + " atm, " + to_string(o.barostat) + "); density averaged over the last " +
               fmt(100 * (1 - o.average_from), 3) + " % of each hold; Tg = " + how + " of " + what + ", error by bootstrap";
  res.notes.push_back("effective cooling rate " + fmt(rate * 1e12, 3) + " K/s: simulated Tg sits above the experimental value, roughly 3 K per decade of rate (Williams–Landel–Ferry)");
  if (!res.fit.note.empty()) res.notes.push_back(res.fit.note);
  return res;
}

PullResult run_pull(System& s, const PullOptions& o) {
  if (!s.cell.valid()) throw FieldError("a pull test needs a periodic cell");
  if (o.rate <= 0 || o.distance <= 0) throw std::invalid_argument("give a positive pulling rate and distance");
  const size_t n = s.atoms.size();
  std::vector<char> held(n, 0), group(n, 0);
  size_t nh = 0, ng = 0;
  for (size_t i = 0; i < n; ++i) {
    if (s.atoms[i].mol == o.surface_mol) held[i] = 1, ++nh;
    else group[i] = 1, ++ng;
  }
  if (!nh || !ng) throw std::invalid_argument("needs a surface (molecule " + std::to_string(o.surface_mol) + ") and a film of other molecules");
  PullResult R;
  R.area = norm(cross(s.cell.a, s.cell.b));
  // a film between the surface and the surface's periodic image slides on (or leaves) two interfaces
  int faces = 1;
  const bool fibre = !o.normal && o.axis == 2;
  if (fibre) {   // pull-out along a fibre: the interface is its side, 2π R L (R: the outermost held atoms from the axis)
    double cx = 0, cy = 0;
    for (size_t i = 0; i < n; ++i)
      if (held[i]) cx += s.atoms[i].pos[0], cy += s.atoms[i].pos[1];
    cx /= double(nh), cy /= double(nh);
    double rmax = 0;
    for (size_t i = 0; i < n; ++i)
      if (held[i]) rmax = std::max(rmax, std::hypot(s.atoms[i].pos[0] - cx, s.atoms[i].pos[1] - cy));
    R.area = 2 * 3.14159265358979 * rmax * norm(s.cell.c);
  } else {
    const double Lz = std::fabs(dot(s.cell.c, cross(s.cell.a, s.cell.b))) / R.area;
    auto h = [&](const Vec3& r) { double f = s.cell.to_fractional(r)[2]; return (f - std::floor(f)) * Lz; };
    double stop = -1e300, sbot = 1e300, ftop = -1e300;
    for (size_t i = 0; i < n; ++i) {
      const double z = h(s.atoms[i].pos);
      if (held[i]) stop = std::max(stop, z), sbot = std::min(sbot, z);
      else ftop = std::max(ftop, z);
    }
    if (ftop > stop && Lz + sbot - ftop < 8.0) faces = 2;
  }
  DynamicsOptions d;
  d.field = o.field ? o.field : std::make_shared<ForceField>(default_forcefield(s));
  d.energy = o.energy;
  d.dt = o.dt;
  d.temperature = o.temperature;
  d.tau_t = o.tau_t;
  d.thermostat = Thermostat::Bussi;
  d.fixed = held;
  d.seed = o.seed;
  d.frame_every = 0;
  if (o.relax_first) {   // a freshly grown film has close contacts: push-off and minimise, the surface held
    RelaxOptions ro;
    ro.field = d.field;
    ro.energy = o.energy;
    ro.fixed = held;
    ro.pushoff = true;
    ro.ftol = 2.0;
    ro.max_iterations = 3000;
    relax(s, ro);
  }
  if (o.equilibrate_ps > 0) {   // settle at the temperature first, nothing pulled
    d.steps = std::max<int64_t>(1, std::llround(o.equilibrate_ps * 1000 / o.dt));
    d.new_velocities = true;
    d.thermo_every = int(d.steps);
    run_dynamics(s, d);
    d.new_velocities = false;
  }
  d.pull_group = group;
  d.pull_dir = o.normal ? Vec3{0, 0, 1} : o.axis == 1 ? Vec3{0, 1, 0} : o.axis == 2 ? Vec3{0, 0, 1} : Vec3{1, 0, 0};
  d.pull_k = o.spring;
  d.pull_rate = o.rate;
  d.steps = std::max<int64_t>(1, std::llround(o.distance / o.rate * 1000 / o.dt));
  d.thermo_every = std::max(1, o.sample_every);
  d.progress = [&](const ThermoRow& r) {
    PullPoint p{r.time_ps, r.pull_disp, r.pull_force, r.temperature};
    R.curve.push_back(p);
    return !o.progress || o.progress(p);
  };
  DynamicsReport dr;
  run_dynamics(s, d, &dr);
  // smoothed force against displacement, the peak, and the work
  const size_t m = R.curve.size();
  R.smooth.resize(m);
  for (size_t i = 0; i < m; ++i) {
    double sum = 0;
    int c = 0;
    for (size_t j = 0; j < m; ++j)
      if (std::fabs(R.curve[j].displacement - R.curve[i].displacement) <= 0.5) sum += R.curve[j].force, ++c;
    R.smooth[i] = c ? sum / c : R.curve[i].force;
  }
  for (size_t i = 0; i < m; ++i)
    if (R.smooth[i] > R.peak_force) R.peak_force = R.smooth[i], R.peak_displacement = R.curve[i].displacement;
  double w = 0;
  for (size_t i = 1; i < m; ++i) w += 0.5 * (R.curve[i].force + R.curve[i - 1].force) * (R.curve[i].displacement - R.curve[i - 1].displacement);
  R.strength = R.peak_force / (faces * R.area) * 6947.7;   // kcal/mol/Å per Å² → MPa
  R.work = w / (faces * R.area) * 694.77;                    // kcal/mol per Å² → mJ/m²
  R.interfaces = faces;
  char b[300];
  std::snprintf(b, sizeof b, "%s pull of %zu atoms from a held surface of %zu atoms · spring %.1f kcal/mol/Å² at %.2f Å/ps over %.1f Å · %.0f K · area %.1f Å²",
                o.normal ? "normal (+z)" : fibre ? "fibre pull-out (z)" : o.axis == 1 ? "shear (y)" : "shear (x)", ng, nh, o.spring, o.rate, o.distance, o.temperature, R.area);
  R.method = b;
  R.notes.push_back("steered MD: the force depends on the pulling rate (far faster than experiment); compare systems at the same rate");
  if (faces == 2) R.notes.push_back("the film touches the surface and its periodic image: strength and work are per interface (two share the force)");
  if (o.normal) R.notes.push_back("normal pulls need vacuum above the film (a film between the surface and its periodic image is pushed into the image)");
  for (const auto& note : dr.notes) R.notes.push_back(note);
  return R;
}

}  // namespace caps
