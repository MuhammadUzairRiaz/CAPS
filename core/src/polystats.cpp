// CAPS polymer statistics: see caps/polystats.hpp.
#include "caps/polystats.hpp"

#include <algorithm>
#include <cmath>
#include <random>

#include "caps/kspace.hpp"

namespace caps {

// ---------------------------------------------------------------- copolymer

double mayo_lewis(double r1, double r2, double f1) {
  const double f2 = 1 - f1, den = r1 * f1 * f1 + 2 * f1 * f2 + r2 * f2 * f2;
  return den > 0 ? (r1 * f1 * f1 + f1 * f2) / den : f1;
}

CopolymerModel copolymer_terminal(double r1, double r2, double f1) {
  CopolymerModel m;
  m.r1 = r1 = std::max(0.0, r1), m.r2 = r2 = std::max(0.0, r2), m.f1 = f1 = std::clamp(f1, 0.0, 1.0);
  const double f2 = 1 - f1;
  m.F1 = mayo_lewis(r1, r2, f1);
  m.paa = r1 * f1 + f2 > 0 ? r1 * f1 / (r1 * f1 + f2) : 0;
  m.pbb = r2 * f2 + f1 > 0 ? r2 * f2 / (r2 * f2 + f1) : 0;
  m.run_a = m.paa < 1 ? 1 / (1 - m.paa) : INFINITY;
  m.run_b = m.pbb < 1 ? 1 / (1 - m.pbb) : INFINITY;
  // F1 = f1 away from 0 and 1: f1 = (1 − r2) / (2 − r1 − r2), both ratios on the same side of 1
  if (std::fabs(2 - r1 - r2) > 1e-12) {
    const double az = (1 - r2) / (2 - r1 - r2);
    if (az > 1e-9 && az < 1 - 1e-9) m.azeotrope = az;
  }
  return m;
}

// ---------------------------------------------------------------- stereo sequences

const std::array<const char*, kPentadCount>& pentad_names() {
  static const std::array<const char*, kPentadCount> n = {"mmmm", "mmmr", "rmmr", "mmrr", "mmrm", "rmrm", "rmrr", "rrrr", "rrrm", "mrrm"};
  return n;
}

namespace {
int pentad_index(const char* s) {   // s: four dyads
  const auto& n = pentad_names();
  const std::string a(s, 4), b(a.rbegin(), a.rend());
  for (int k = 0; k < kPentadCount; ++k)
    if (a == n[size_t(k)] || b == n[size_t(k)]) return k;
  return -1;
}

// Probabilities of a model's sequences: stationary start, then the transition matrix
void fill_model(StereoModel& m) {
  const double a = std::clamp(m.p_mr, 0.0, 1.0), b = std::clamp(m.p_rm, 0.0, 1.0);
  m.pm = a + b > 0 ? b / (a + b) : 0.5;
  auto P = [&](const char* s, int n) {
    double p = s[0] == 'm' ? m.pm : 1 - m.pm;
    for (int i = 1; i < n; ++i) {
      const bool from_m = s[i - 1] == 'm', to_m = s[i] == 'm';
      p *= from_m ? (to_m ? 1 - a : a) : (to_m ? b : 1 - b);
    }
    return p;
  };
  m.mm = P("mm", 2), m.rr = P("rr", 2), m.mr = P("mr", 2) + P("rm", 2);
  m.pentads.fill(0);
  char s[5] = {0};
  for (int bits = 0; bits < 16; ++bits) {
    for (int i = 0; i < 4; ++i) s[i] = (bits >> (3 - i)) & 1 ? 'r' : 'm';
    m.pentads[size_t(pentad_index(s))] += P(s, 4);
  }
}
}  // namespace

StereoModel stereo_bernoulli(double pm) {
  StereoModel m;
  m.kind = "bernoulli";
  pm = std::clamp(pm, 0.0, 1.0);
  m.p_mr = 1 - pm, m.p_rm = pm;
  fill_model(m);
  m.pm = pm;
  return m;
}

StereoModel stereo_markov(double p_mr, double p_rm) {
  StereoModel m;
  m.kind = "markov";
  m.p_mr = std::clamp(p_mr, 0.0, 1.0), m.p_rm = std::clamp(p_rm, 0.0, 1.0);
  fill_model(m);
  return m;
}

StereoCounts count_stereo(const std::string& d) {
  StereoCounts c;
  for (char x : d) (x == 'm' ? c.m : c.r)++;
  for (size_t i = 1; i < d.size(); ++i) {
    const std::string t = d.substr(i - 1, 2);
    if (t == "mm") ++c.mm;
    else if (t == "rr") ++c.rr;
    else ++c.mr;
  }
  for (size_t i = 3; i < d.size(); ++i) {
    const int k = pentad_index(d.c_str() + i - 3);
    if (k >= 0) ++c.pentads[size_t(k)], ++c.pentad_total;
  }
  return c;
}

namespace {
std::array<double, kPentadCount> normalised(const std::array<double, kPentadCount>& x) {
  double s = 0;
  for (double v : x) s += std::max(0.0, v);
  std::array<double, kPentadCount> out{};
  for (int k = 0; k < kPentadCount; ++k) out[size_t(k)] = s > 0 ? std::max(0.0, x[size_t(k)]) / s : 0;
  return out;
}
double residual(const StereoModel& m, const std::array<double, kPentadCount>& y) {
  double s = 0;
  for (int k = 0; k < kPentadCount; ++k) s += (m.pentads[size_t(k)] - y[size_t(k)]) * (m.pentads[size_t(k)] - y[size_t(k)]);
  return std::sqrt(s / kPentadCount);
}
}  // namespace

StereoModel fit_bernoulli(const std::array<double, kPentadCount>& measured, double* rms) {
  const auto y = normalised(measured);
  double best = 1e300, bp = 0.5;
  for (int i = 0; i <= 1000; ++i) {   // grid, then golden-section refinement around the best point
    const double p = i / 1000.0, r = residual(stereo_bernoulli(p), y);
    if (r < best) best = r, bp = p;
  }
  double lo = std::max(0.0, bp - 0.001), hi = std::min(1.0, bp + 0.001);
  const double g = (std::sqrt(5.0) - 1) / 2;
  for (int it = 0; it < 60; ++it) {
    const double a = hi - g * (hi - lo), b = lo + g * (hi - lo);
    if (residual(stereo_bernoulli(a), y) < residual(stereo_bernoulli(b), y)) hi = b;
    else lo = a;
  }
  auto m = stereo_bernoulli((lo + hi) / 2);
  if (rms) *rms = residual(m, y);
  return m;
}

StereoModel fit_markov(const std::array<double, kPentadCount>& measured, double* rms) {
  const auto y = normalised(measured);
  double best = 1e300, ba = 0.5, bb = 0.5;
  for (int i = 0; i <= 200; ++i)
    for (int j = 0; j <= 200; ++j) {
      const double r = residual(stereo_markov(i / 200.0, j / 200.0), y);
      if (r < best) best = r, ba = i / 200.0, bb = j / 200.0;
    }
  double step = 0.005;   // pattern search refinement
  while (step > 1e-7) {
    bool moved = false;
    for (auto [da, db] : {std::pair{step, 0.0}, {-step, 0.0}, {0.0, step}, {0.0, -step}}) {
      const double a = std::clamp(ba + da, 0.0, 1.0), b = std::clamp(bb + db, 0.0, 1.0), r = residual(stereo_markov(a, b), y);
      if (r < best) best = r, ba = a, bb = b, moved = true;
    }
    if (!moved) step /= 2;
  }
  auto m = stereo_markov(ba, bb);
  if (rms) *rms = best;
  return m;
}

std::string draw_dyads(const StereoModel& m, int dp, uint64_t seed) {
  std::mt19937_64 rng(seed * 0x9E3779B97F4A7C15ull + 5);
  auto u = [&] { return double(rng() >> 11) * (1.0 / 9007199254740992.0); };
  std::string d;
  for (int i = 0; i + 1 < dp; ++i) {
    if (i == 0) d += u() < m.pm ? 'm' : 'r';
    else if (d.back() == 'm') d += u() < m.p_mr ? 'r' : 'm';
    else d += u() < m.p_rm ? 'm' : 'r';
  }
  return d;
}

// ---------------------------------------------------------------- Flory–Huggins blend

BlendCritical blend_critical(double na, double nb) {
  const double a = 1 / std::sqrt(na), b = 1 / std::sqrt(nb);
  return {0.5 * (a + b) * (a + b), std::sqrt(nb) / (std::sqrt(na) + std::sqrt(nb))};
}

bool blend_spinodal(double na, double nb, double chi, double& lo, double& hi) {
  // 1/(NA φ) + 1/(NB (1 − φ)) = 2χ  →  2χ NA NB φ² + (NA − NB − 2χ NA NB) φ + NB = 0
  const double A = 2 * chi * na * nb, B = na - nb - 2 * chi * na * nb, C = nb, disc = B * B - 4 * A * C;
  if (A <= 0 || disc < 0) return false;
  const double s = std::sqrt(disc);
  lo = (-B - s) / (2 * A), hi = (-B + s) / (2 * A);
  return lo > 0 && hi < 1 && lo < hi;
}

bool blend_binodal(double na, double nb, double chi, double& lo, double& hi, double glo, double ghi) {
  const auto crit = blend_critical(na, nb);
  if (chi <= crit.chi_c * (1 + 1e-9)) return false;
  // Everything in the logit x = ln(φ/(1 − φ)): φ, 1 − φ, ln φ and ln(1 − φ) straight from x, so compositions within
  // 1e-300 of 0 or 1 (strong segregation, χN ≫ 1) keep their precision.
  struct Q { double p, q, lp, lq; };
  auto at = [](double x) {
    const double sp = x > 0 ? x + std::log1p(std::exp(-x)) : std::log1p(std::exp(x));   // softplus(x) = ln(1 + eˣ)
    return Q{1 / (1 + std::exp(-x)), 1 / (1 + std::exp(x)), x - sp, -sp};
  };
  // f'(φ) and the grand potential g = f − φ f'
  auto fp = [&](const Q& s) { return (s.lp + 1) / na - (s.lq + 1) / nb + chi * (s.q - s.p); };
  auto g = [&](const Q& s) { return s.p / na * s.lp + s.q / nb * s.lq + chi * s.p * s.q - s.p * fp(s); };
  double x1, x2;   // logits of the two compositions
  if (glo > 0 && ghi < 1 && glo < ghi) x1 = std::log(glo / (1 - glo)), x2 = std::log(ghi / (1 - ghi));
  else {
    double s1, s2;
    if (!blend_spinodal(na, nb, chi, s1, s2)) return false;
    const double p1 = std::max(1e-12, crit.phi_c - std::sqrt(3.0) * (crit.phi_c - s1));   // Landau: binodal ≈ √3 × spinodal width
    const double p2 = std::min(1 - 1e-12, crit.phi_c + std::sqrt(3.0) * (s2 - crit.phi_c));
    x1 = std::log(p1 / (1 - p1)), x2 = std::log(p2 / (1 - p2));
  }
  for (int it = 0; it < 3000; ++it) {   // steps of at most 2 in x: χN of a few thousand is reached
    const Q a = at(x1), b = at(x2);
    const double F1 = fp(a) - fp(b), F2 = g(a) - g(b);
    if (std::fabs(F1) < 1e-13 && std::fabs(F2) < 1e-13) break;
    // d/dx = dφ/dx · d/dφ with dφ/dx = φ(1 − φ); dg/dφ = −φ f''; f''·φ(1 − φ) formed without the huge factor
    auto fd = [&](const Q& s) { return s.q / na + s.p / nb - 2 * chi * s.p * s.q; };   // f'' φ (1 − φ)
    const double J11 = fd(a), J12 = -fd(b), J21 = -a.p * fd(a), J22 = b.p * fd(b);
    const double det = J11 * J22 - J12 * J21;
    if (std::fabs(det) < 1e-300) return false;
    double dx1 = (-F1 * J22 + F2 * J12) / det, dx2 = (-F2 * J11 + F1 * J21) / det;
    const double big = std::max(std::fabs(dx1), std::fabs(dx2));
    if (big > 2) dx1 *= 2 / big, dx2 *= 2 / big;
    x1 += dx1, x2 += dx2;
  }
  if (x1 > x2) std::swap(x1, x2);
  const Q a = at(x1), b = at(x2);
  lo = a.p, hi = b.p;
  return std::fabs(fp(a) - fp(b)) < 1e-9 && std::fabs(g(a) - g(b)) < 1e-9 && x2 - x1 > 1e-5;
}

// ---------------------------------------------------------------- solvents, Ewald

double hildebrand_chi(double v, double ds, double dp, double t) {
  return v * (ds - dp) * (ds - dp) / (8.314462618 * t) + 0.34;
}

int pme_mesh_size(double edge, double spacing, int order) {
  return fft_good_size(std::max(order + 1, int(std::ceil(edge / std::max(1e-6, spacing)))));
}

std::vector<double> ris_cn(const RisModel& m, double T, int nmax) {
  constexpr double kPiR = 3.14159265358979323846, kRk = 0.0019872043;
  const double l = m.bond, th = kPiR - m.angle_deg * kPiR / 180;   // θ: the supplement of the bond angle
  const double sg = std::exp(-m.e_sigma / (kRk * T)), om = std::exp(-m.e_omega / (kRk * T));
  const double U[3][3] = {{1, sg, sg}, {1, sg, sg * om}, {1, sg * om, sg}};   // rows: bond i−1 in t g+ g−, columns: bond i
  const double phis[3] = {0, m.gauche_deg * kPiR / 180, -m.gauche_deg * kPiR / 180};   // trans at 0 (Flory's convention)
  // T(φ): the frame of bond i+1 in the frame of bond i
  auto Tm = [&](double ph) {
    const double c = std::cos(th), s = std::sin(th), cp = std::cos(ph), sp = std::sin(ph);
    return std::array<std::array<double, 3>, 3>{{{c, s, 0}, {s * cp, -c * cp, sp}, {s * sp, -c * sp, -cp}}};
  };
  // generator matrix for one state: rows / columns (1, 3 vector, 1) = 5
  using M5 = std::array<std::array<double, 5>, 5>;
  auto G = [&](double ph) {
    const auto t = Tm(ph);
    M5 g{};
    g[0][0] = 1;
    for (int k = 0; k < 3; ++k) g[0][1 + k] = 2 * l * t[0][k];   // 2 lᵀ T with l = (l, 0, 0)
    g[0][4] = l * l;
    for (int a = 0; a < 3; ++a)
      for (int b = 0; b < 3; ++b) g[1 + a][1 + b] = t[a][b];
    g[1][4] = l;
    g[4][4] = 1;
    return g;
  };
  std::array<M5, 3> Gs{G(phis[0]), G(phis[1]), G(phis[2])};
  // the first bond's generator row [1, 2 l₁ᵀ T₁, l²] (its torsion is undefined: taken trans), carried in state t; J* = (1, 0, 0)
  std::vector<double> v(15, 0.0);
  {
    const auto t1 = Tm(0);
    v[0] = 1, v[4] = l * l;
    for (int k = 0; k < 3; ++k) v[size_t(1 + k)] = 2 * l * t1[0][size_t(k)];
  }
  std::array<double, 3> zrow{1, 0, 0};
  std::vector<double> out{1.0};   // C_1 = 1
  double z = 1;
  for (int n = 2; n <= nmax; ++n) {
    // close with the last bond's column [l², l_n, 1] (l_n = (l, 0, 0)), summed over the states (J)
    double r2 = 0;
    for (int st = 0; st < 3; ++st) r2 += v[size_t(5 * st)] * l * l + v[size_t(5 * st + 1)] * l + v[size_t(5 * st + 4)];
    z = zrow[0] + zrow[1] + zrow[2];
    out.push_back(r2 / z / (double(n) * l * l));
    // one more rotatable bond: v ← v · (U ⊗ E5) · diag(G_t, G_g+, G_g−); zrow ← zrow · U
    std::vector<double> w(15, 0.0);
    std::array<double, 3> zn{0, 0, 0};
    for (int a = 0; a < 3; ++a)
      for (int b = 0; b < 3; ++b) {
        zn[size_t(b)] += zrow[size_t(a)] * U[a][b];
        for (int p = 0; p < 5; ++p) {
          const double x = v[size_t(5 * a + p)] * U[a][b];
          if (x == 0) continue;
          for (int q = 0; q < 5; ++q) w[size_t(5 * b + q)] += x * Gs[size_t(b)][size_t(p)][size_t(q)];
        }
      }
    // the partition sums grow geometrically: rescale both by the same factor (the ratio is what counts)
    const double sc = zn[0] + zn[1] + zn[2];
    for (auto& x : w) x /= sc;
    for (auto& x : zn) x /= sc;
    v = w;
    zrow = zn;
  }
  (void)z;
  return out;
}

}  // namespace caps
