// CAPS Analyze: properties from trajectories (see caps/properties.hpp).
#include "caps/properties.hpp"

#include <algorithm>
#include <array>
#include <complex>
#include <cstdio>
#include <cmath>
#include <memory>
#include <tuple>
#include <numeric>
#include <set>
#include <stdexcept>

#include "caps/analysis.hpp"
#include "caps/elements.hpp"
#include "caps/json.hpp"
#include "caps/mechanics.hpp"

namespace caps {

namespace {

constexpr double kNA = 6.02214076e23;
constexpr double kPi = 3.14159265358979323846;
const double NaN = std::numeric_limits<double>::quiet_NaN();

// Mean and standard error of the mean from `blocks` consecutive blocks (Flyvbjerg–Petersen style, one level).
std::pair<double, double> block_mean(const std::vector<double>& v, int blocks) {
  if (v.empty()) return {NaN, NaN};
  const double mean = std::accumulate(v.begin(), v.end(), 0.0) / double(v.size());
  const int nb = std::min<int>(blocks, int(v.size()));
  if (nb < 2) return {mean, NaN};
  std::vector<double> bm(nb, 0.0);
  const size_t per = v.size() / size_t(nb);
  for (int b = 0; b < nb; ++b) {
    double s = 0;
    for (size_t k = size_t(b) * per; k < size_t(b + 1) * per; ++k) s += v[k];
    bm[size_t(b)] = s / double(per);
  }
  double var = 0;
  const double bmean = std::accumulate(bm.begin(), bm.end(), 0.0) / nb;
  for (double x : bm) var += (x - bmean) * (x - bmean);
  double err = std::sqrt(var / (nb - 1) / nb);
  if (err < 1e-12 * std::fabs(mean)) err = 0;   // identical blocks: rounding noise only
  return {mean, err};
}

// Least-squares line y = a + b x; returns {a, b}.
std::pair<double, double> linfit(const std::vector<double>& x, const std::vector<double>& y) {
  const size_t n = x.size();
  if (n < 2) return {NaN, NaN};
  double sx = 0, sy = 0, sxx = 0, sxy = 0;
  for (size_t i = 0; i < n; ++i) { sx += x[i]; sy += y[i]; sxx += x[i] * x[i]; sxy += x[i] * y[i]; }
  const double d = n * sxx - sx * sx;
  if (std::fabs(d) < 1e-300) return {NaN, NaN};
  const double b = (n * sxy - sx * sy) / d;
  return {(sy - b * sx) / n, b};
}

double shortest_width(const Cell& c) {
  if (!c.valid()) return 0;
  const double v = c.volume();
  return std::min({v / norm(cross(c.b, c.c)), v / norm(cross(c.c, c.a)), v / norm(cross(c.a, c.b))});
}

bool cancelled(const AnalyzeOptions& o, const std::string& what, double f) { return o.progress && !o.progress(what, f); }

struct Cancel : std::runtime_error {
  Cancel() : std::runtime_error("analysis cancelled") {}
};

// Unwrapped positions of every selected frame: from the file when it is unwrapped, otherwise by following each atom
// with minimum-image steps between consecutive frames (frames must be close enough that no atom moves half a cell).
std::vector<std::vector<Vec3>> unwrapped_positions(const Trajectory& t, const std::vector<size_t>& fr, std::vector<std::string>& notes) {
  std::vector<std::vector<Vec3>> out;
  out.reserve(fr.size());
  const bool file_unwrapped = t.topology.unwrapped;
  for (size_t k = 0; k < fr.size(); ++k) {
    const System s = t.frame(fr[k]);
    std::vector<Vec3> p(s.atoms.size());
    for (size_t i = 0; i < p.size(); ++i) p[i] = s.atoms[i].pos;
    if (!file_unwrapped && k > 0 && s.cell.valid())
      for (size_t i = 0; i < p.size(); ++i) p[i] = out.back()[i] + s.cell.minimum_image(p[i] - out.back()[i]);
    out.push_back(std::move(p));
  }
  if (!file_unwrapped) notes.push_back("positions unwrapped by minimum-image steps between frames (the file stores wrapped positions)");
  return out;
}

// ---------------------------------------------------------------- structure

Property density_prop(const Trajectory& t, const std::vector<size_t>& fr, const std::vector<double>& times, const AnalyzeOptions& o) {
  Property p{"density", "Density", "g/cm³", "mass / cell volume over the frames; standard error from block averages", NaN, NaN, {}, {}, {}};
  const double mass = t.topology.total_mass();
  Series s{"density", "time (ps)", "density (g/cm³)", {}, {}};
  std::vector<double> v;
  for (size_t k = 0; k < fr.size(); ++k) {
    const Cell& c = fr[k] < t.cells.size() ? t.cells[fr[k]] : t.topology.cell;
    if (!c.valid()) continue;
    const double rho = mass / kNA / (c.volume() * 1e-24);
    v.push_back(rho);
    s.x.push_back(times[k]);
    s.y.push_back(rho);
  }
  if (v.empty()) { p.notes.push_back("no periodic cell: density is undefined"); return p; }
  std::tie(p.value, p.error) = block_mean(v, o.blocks);
  p.series.push_back(std::move(s));
  return p;
}

// Frame-averaged partial g(r) between two element sets (0 = any).
std::vector<std::pair<double, double>> avg_rdf(const Trajectory& t, const std::vector<size_t>& fr, int ea, int eb, double rmax, double dr, bool inter,
                                               const AnalyzeOptions& o, const std::string& what) {
  std::vector<std::pair<double, double>> acc;
  for (size_t k = 0; k < fr.size(); ++k) {
    if (cancelled(o, what, double(k) / fr.size())) throw Cancel();
    const System s = t.frame(fr[k]);
    const auto g = rdf(s, ea, eb, rmax, dr, inter);
    if (acc.empty()) acc = g;
    else
      for (size_t b = 0; b < g.size(); ++b) acc[b].second += g[b].second;
  }
  for (auto& x : acc) x.second /= double(fr.size());
  return acc;
}

double rdf_rmax(const Trajectory& t, const AnalyzeOptions& o) {
  const double w = shortest_width(t.topology.cell.valid() ? t.topology.cell : (t.cells.empty() ? Cell{} : t.cells[0]));
  double r = o.rdf_rmax > 0 ? o.rdf_rmax : std::min(15.0, 0.5 * w);
  if (w > 0) r = std::min(r, 0.5 * w - 1e-6);
  return r;
}

Property rdf_prop(const Trajectory& t, const std::vector<size_t>& fr, const AnalyzeOptions& o) {
  Property p{"rdf", "Radial distribution", "Å", "", NaN, NaN, {}, {}, {}};
  if (!t.topology.cell.valid()) { p.notes.push_back("g(r) needs a periodic cell"); return p; }
  const double rmax = rdf_rmax(t, o);
  const auto g = avg_rdf(t, fr, o.elem_a, o.elem_b, rmax, o.rdf_dr, o.inter_only, o, "g(r)");
  auto name = [](int z) { return z ? std::string(element(z).symbol) : std::string("all"); };
  p.method = "g(r) " + name(o.elem_a) + "–" + name(o.elem_b) + (o.inter_only ? ", between molecules" : "") + ", averaged over " +
             std::to_string(fr.size()) + " frames, bins of " + std::to_string(o.rdf_dr).substr(0, 5) + " Å";
  Series s{"g(r)", "r (Å)", "g(r)", {}, {}}, n{"n(r)", "r (Å)", "coordination number", {}, {}};
  // coordination number n(r) = 4π ρ_B ∫ g r² dr
  size_t nb = 0;
  for (const auto& a : t.topology.atoms) nb += !o.elem_b || a.element == o.elem_b;
  const double rho = nb / t.topology.cell.volume();
  double cn = 0;
  for (const auto& [r, v] : g) {
    s.x.push_back(r);
    s.y.push_back(v);
    cn += 4 * kPi * rho * v * r * r * o.rdf_dr;
    n.x.push_back(r + o.rdf_dr / 2);
    n.y.push_back(cn);
  }
  // first peak: the highest bin after g first exceeds 1 (after the exclusion zone)
  size_t start = 0;
  while (start < g.size() && g[start].second < 1.0) ++start;
  size_t best = start;
  for (size_t b = start; b < g.size() && g[b].first < g[start].first + 1.5; ++b)
    if (g[b].second > g[best].second) best = b;
  if (start < g.size()) {
    p.value = g[best].first;
    p.extra["first peak height"] = g[best].second;
    // first minimum after the peak and the coordination number there
    size_t mn = best;
    for (size_t b = best; b < g.size() && g[b].first < g[best].first + 3.0; ++b)
      if (g[b].second < g[mn].second) mn = b;
    p.extra["first minimum (Å)"] = g[mn].first;
    p.extra["coordination number at the first minimum"] = n.y[mn];
  }
  p.series.push_back(std::move(s));
  p.series.push_back(std::move(n));
  return p;
}

// Direct S(q) on the cell's reciprocal lattice: S(k) = 1 + (|Σ w_j e^{ik·r_j}|²/N − ⟨w²⟩)/⟨w⟩² (Faber–Ziman weights), averaged
// over |k| shells of width dq, for |k| ≤ qcut. Exact for the periodic model (no g(r) truncation); frames are spread so that
// frames × k-vectors × atoms stays within a budget. Returns shell centres and values (empty shells skipped) and the frames used.
std::vector<size_t> spread_frames(const std::vector<size_t>& fr, size_t n);

std::string fmt(double v, int digits) {
  char b[32];
  std::snprintf(b, sizeof b, "%.*g", digits, v);
  return b;
}

struct DirectSq { std::vector<double> q, s; size_t frames = 0; double qcut = 0; };
DirectSq direct_sq(const Trajectory& t, const std::vector<size_t>& fr0, const AnalyzeOptions& o, const std::string& kind, double qcut) {
  DirectSq out;
  const size_t n = t.topology.atoms.size();
  const double V = t.topology.cell.volume();
  // half-space k-vectors up to qcut: about qcut³ V / (12 π²); keep frames × k × atoms ≲ 3·10⁹
  const double budget = 3e9;
  double nk = qcut * qcut * qcut * V / (12 * kPi * kPi);
  if (nk * n > budget) { qcut = std::cbrt(budget / n * 12 * kPi * kPi / V); nk = qcut * qcut * qcut * V / (12 * kPi * kPi); }
  out.qcut = qcut;
  const size_t nfr = std::max<size_t>(1, std::min<size_t>(20, size_t(budget / std::max(1.0, nk * n))));
  const auto fr = spread_frames(fr0, nfr);
  out.frames = fr.size();
  const size_t nb = size_t(qcut / o.dq) + 1;
  std::vector<double> sum(nb, 0.0), cnt(nb, 0.0), ksum(nb, 0.0);
  std::vector<int> el(n);
  for (size_t j = 0; j < n; ++j) el[j] = t.topology.atoms[j].element;
  for (size_t q = 0; q < fr.size(); ++q) {
    if (cancelled(o, kind == "sq" ? "S(q)" : "scattering", double(q) / fr.size())) throw Cancel();
    const System f = t.frame(fr[q]);
    const Cell& c = f.cell;
    const double vol = c.volume();
    const Vec3 b1 = cross(c.b, c.c) * (2 * kPi / vol), b2 = cross(c.c, c.a) * (2 * kPi / vol), b3 = cross(c.a, c.b) * (2 * kPi / vol);
    const int hm = int(qcut * norm(c.a) / (2 * kPi)) + 1, km = int(qcut * norm(c.b) / (2 * kPi)) + 1, lm = int(qcut * norm(c.c) / (2 * kPi)) + 1;
    // e^{i h (b1·r)} etc. for every atom by recurrence
    using cd = std::complex<double>;
    auto powers = [&](const Vec3& bv, int m) {
      std::vector<cd> e(size_t(2 * m + 1) * n);
      for (size_t j = 0; j < n; ++j) {
        const cd e1 = std::polar(1.0, dot(bv, f.atoms[j].pos));
        const cd e1c = std::conj(e1);
        e[size_t(m) * n + j] = 1.0;
        for (int h = 1; h <= m; ++h) {
          e[size_t(m + h) * n + j] = e[size_t(m + h - 1) * n + j] * e1;
          e[size_t(m - h) * n + j] = e[size_t(m - h + 1) * n + j] * e1c;
        }
      }
      return e;
    };
    const auto E1 = powers(b1, hm), E2 = powers(b2, km), E3 = powers(b3, lm);
    std::map<int, double> wz;
    std::vector<cd> ab(n);
    for (int h = 0; h <= hm; ++h)
      for (int k = (h == 0 ? 0 : -km); k <= km; ++k) {
        const Vec3 hk = b1 * double(h) + b2 * double(k);
        if (norm(hk) > qcut + norm(b3) * lm) continue;
        for (size_t j = 0; j < n; ++j) ab[j] = E1[size_t(hm + h) * n + j] * E2[size_t(km + k) * n + j];
        for (int l = (h == 0 && k == 0 ? 1 : -lm); l <= lm; ++l) {
          const double kn = norm(hk + b3 * double(l));
          if (kn > qcut || kn <= 0) continue;
          // weights at this |k|
          double wsum = 0, w2 = 0;
          if (kind == "sq") { wsum = 1; w2 = 1; }
          else wz.clear();   // form factors depend on |k|
          cd A = 0;
          if (kind == "sq")
            for (size_t j = 0; j < n; ++j) A += ab[j] * E3[size_t(lm + l) * n + j];
          else {
            double sw = 0, sw2 = 0;
            for (size_t j = 0; j < n; ++j) {
              auto it = wz.find(el[j]);
              if (it == wz.end()) it = wz.emplace(el[j], kind == "xray" ? xray_f(el[j], kn) : neutron_b(el[j])).first;
              A += it->second * ab[j] * E3[size_t(lm + l) * n + j];
              sw += it->second;
              sw2 += it->second * it->second;
            }
            wsum = sw / n;
            w2 = sw2 / n;
          }
          const double S = 1 + (std::norm(A) / n - w2) / (wsum * wsum);
          const size_t bi = size_t(kn / o.dq);
          if (bi < nb) { sum[bi] += S; cnt[bi] += 1; ksum[bi] += kn; }
        }
      }
  }
  // adjacent shells merge until each holds at least kMinShell k-vectors per frame (at low q a shell has only a few); the point
  // sits at the mean |k| of its vectors
  constexpr double kMinShell = 24;
  double S = 0, C = 0, K = 0;
  for (size_t bi = 0; bi < nb; ++bi) {
    S += sum[bi]; C += cnt[bi]; K += ksum[bi];
    if (C >= kMinShell * fr.size() || (bi + 1 == nb && C > 0)) {
      out.q.push_back(K / C);
      out.s.push_back(S / C);
      S = C = K = 0;
    }
  }
  return out;
}

// Faber–Ziman total structure factor from element partials, weights w(z, q). Lorch window against truncation ripples.
Property scattering_prop(const Trajectory& t, const std::vector<size_t>& fr, const AnalyzeOptions& o, const std::string& kind) {
  const std::string title = kind == "sq" ? "Structure factor S(q)" : kind == "xray" ? "X-ray scattering" : "Neutron scattering";
  Property p{kind, title, "Å⁻¹", "", NaN, NaN, {}, {}, {}};
  if (!t.topology.cell.valid()) { p.notes.push_back("the structure factor needs a periodic cell"); return p; }
  const double rmax = rdf_rmax(t, o), dr = o.rdf_dr;
  std::map<int, double> count;
  for (const auto& a : t.topology.atoms) count[a.element] += 1;
  const double natoms = double(t.topology.atoms.size());
  const double rho = natoms / t.topology.cell.volume();
  std::vector<int> el;
  for (const auto& [z, c] : count) el.push_back(z);
  // partial g_ab(r) for a ≤ b (number-weighted S(q) uses the all-atom g(r) directly)
  std::map<std::pair<int, int>, std::vector<std::pair<double, double>>> g;
  if (kind == "sq") g[{0, 0}] = avg_rdf(t, fr, 0, 0, rmax, dr, false, o, "S(q)");
  else
    for (size_t a = 0; a < el.size(); ++a)
      for (size_t b = a; b < el.size(); ++b)
        g[{el[a], el[b]}] = avg_rdf(t, fr, el[a], el[b], rmax, dr, false, o, title);
  auto partial_sq = [&](const std::vector<std::pair<double, double>>& gab, double q) {
    // S_ab(q) − 1 = 4πρ ∫ r² (g − 1) sin(qr)/(qr) L(r) dr
    double s = 0;
    for (const auto& [r, v] : gab) {
      const double x = kPi * r / rmax, lorch = x > 0 ? std::sin(x) / x : 1.0;
      s += r * r * (v - 1) * (q > 0 ? std::sin(q * r) / (q * r) : 1.0) * lorch * dr;
    }
    return 4 * kPi * rho * s;
  };
  Series s{kind == "sq" ? "S(q)" : kind == "xray" ? "I(q) X-ray (Faber–Ziman)" : "S(q) neutron (Faber–Ziman)", "q (Å⁻¹)",
           kind == "sq" ? "S(q)" : "S(q), normalised", {}, {}};
  const double qmin = 2 * kPi / rmax;
  // low q: exact reciprocal-lattice sum on the cell; above it the g(r) transform (short-range order, truncation harmless)
  const DirectSq D = o.q_direct > 0 ? direct_sq(t, fr, o, kind, std::min(o.q_direct, o.qmax)) : DirectSq{};
  for (size_t k = 0; k < D.q.size(); ++k) { s.x.push_back(D.q[k]); s.y.push_back(D.s[k]); }
  const double qstart = D.q.empty() ? std::max(o.dq, 0.5 * qmin) : D.qcut + o.dq;
  for (double q = qstart; q <= o.qmax + 1e-9; q += o.dq) {
    double F = 0;
    if (kind == "sq") F = partial_sq(g[{0, 0}], q);
    else {
      double wsum = 0;
      for (int z : el) wsum += count[z] / natoms * (kind == "xray" ? xray_f(z, q) : neutron_b(z));
      for (size_t a = 0; a < el.size(); ++a)
        for (size_t b = a; b < el.size(); ++b) {
          const double ca = count[el[a]] / natoms, cb = count[el[b]] / natoms;
          const double wa = kind == "xray" ? xray_f(el[a], q) : neutron_b(el[a]), wb = kind == "xray" ? xray_f(el[b], q) : neutron_b(el[b]);
          // Faber–Ziman: F = Σ_ab c_a c_b w_a w_b [S_ab − 1]; g_ab is normalised by ρ_b = c_b ρ, so S_ab − 1 uses the total ρ.
          // Unlike pairs appear twice in the double sum.
          F += (a == b ? 1.0 : 2.0) * ca * cb * wa * wb * partial_sq(g[{el[a], el[b]}], q);
        }
      F /= wsum * wsum;
    }
    s.x.push_back(q);
    s.y.push_back(1 + F);
  }
  // first peak above q = 0.5 Å⁻¹ (below it the finite box and the window dominate)
  size_t best = 0;
  bool found = false;
  for (size_t k = 0; k < s.x.size(); ++k)
    if (s.x[k] >= 0.5 && (!found || s.y[k] > s.y[best])) {
      if (s.x[k] > 3.0 && found) break;
      best = k;
      found = true;
    }
  if (found) {
    p.value = s.x[best];
    p.extra["first peak height"] = s.y[best];
    p.extra["d-spacing 2π/q (Å)"] = 2 * kPi / s.x[best];
  }
  const std::string w = kind == "sq" ? "number-weighted" : kind == "xray" ? "X-ray (Cromer–Mann form factors)" : "neutron (coherent scattering lengths)";
  if (!D.q.empty())
    p.method = w + " Faber–Ziman S(q): direct sum over the cell's reciprocal lattice up to " + fmt(D.qcut, 3) + " Å⁻¹ (" +
               std::to_string(D.frames) + " frames, |k| shells of " + fmt(o.dq, 3) + " Å⁻¹ merged to ≥ 24 vectors), above it from g(r) up to " + fmt(rmax, 4) +
               " Å with a Lorch window; the smallest q of this cell is " + fmt(2 * kPi / std::cbrt(t.topology.cell.volume()), 3) + " Å⁻¹";
  else
    p.method = w + " Faber–Ziman S(q) from g(r) up to " + fmt(rmax, 4) + " Å with a Lorch window; q below " + fmt(qmin, 3) +
               " Å⁻¹ is not resolved by this cell";
  if (kind == "neutron") {
    double b1 = 0, b2 = 0;
    for (const auto& [z, c] : count) { b1 += c / natoms * neutron_b(z); b2 += c / natoms * neutron_b(z) * neutron_b(z); }
    if (b1 * b1 < 0.25 * b2)
      p.notes.push_back("the mean scattering length is small (" + fmt(b1, 3) + " fm: hydrogen scatters with b < 0), so the normalised S(q) is "
                        "magnified; measurements on such materials usually use deuterated samples");
  }
  p.series.push_back(std::move(s));
  return p;
}

// ---------------------------------------------------------------- chains

struct ChainFrames {
  std::vector<std::vector<uint32_t>> bb;   // backbones (heavy atoms, in order)
  std::vector<std::vector<Vec3>> pos;      // whole-molecule positions per frame
};

ChainFrames chain_frames(const Trajectory& t, const std::vector<size_t>& fr) {
  ChainFrames c;
  System s0 = t.frame(fr.empty() ? 0 : fr[0]);
  c.bb = backbones(s0);
  for (size_t k : fr) {
    System s = t.frame(k);
    if (!s.unwrapped) make_molecules_whole(s);
    std::vector<Vec3> p(s.atoms.size());
    for (size_t i = 0; i < p.size(); ++i) p[i] = s.atoms[i].pos;
    c.pos.push_back(std::move(p));
  }
  return c;
}

Property rg_prop(const Trajectory& t, const std::vector<size_t>& fr, const AnalyzeOptions& o) {
  Property p{"rg", "Radius of gyration", "Å", "√⟨Rg²⟩ over molecules with at least 4 heavy atoms and frames (mass-weighted)", NaN, NaN, {}, {}, {}};
  std::vector<double> per_frame, all;
  for (size_t k : fr) {
    System s = t.frame(k);
    if (!s.unwrapped) make_molecules_whole(s);
    const auto shapes = molecule_shapes(s);
    const auto bb = backbones(s);
    std::set<int> chains;
    const auto mol = s.molecules();
    for (const auto& b : bb) chains.insert(mol[b[0]]);
    double sum = 0;
    int n = 0;
    for (const auto& m : shapes)
      if (chains.count(m.molecule)) { sum += m.rg * m.rg; all.push_back(m.rg); ++n; }
    if (n) per_frame.push_back(sum / n);
  }
  if (per_frame.empty()) { p.notes.push_back("no chains (molecules with at least 4 heavy atoms)"); return p; }
  const auto [m2, e2] = block_mean(per_frame, o.blocks);
  p.value = std::sqrt(m2);
  p.error = std::isnan(e2) ? NaN : e2 / (2 * p.value);
  // distribution
  Series h{"P(Rg)", "Rg (Å)", "probability density", {}, {}};
  const double mx = *std::max_element(all.begin(), all.end()) * 1.05 + 1e-9, bw = mx / 40;
  std::vector<double> hist(40, 0.0);
  for (double x : all) hist[std::min<size_t>(39, size_t(x / bw))] += 1;
  for (size_t b = 0; b < 40; ++b) { h.x.push_back((b + 0.5) * bw); h.y.push_back(hist[b] / (all.size() * bw)); }
  p.series.push_back(std::move(h));
  p.extra["chains"] = double(all.size()) / double(fr.size());
  return p;
}

Property ree_prop(const Trajectory& t, const std::vector<size_t>& fr, const AnalyzeOptions& o, const ChainFrames& c) {
  Property p{"ree", "End-to-end distance", "Å", "√⟨R²⟩ of the backbones (first to last heavy atom)", NaN, NaN, {}, {}, {}};
  if (c.bb.empty()) { p.notes.push_back("no chains"); return p; }
  std::vector<double> per_frame, all;
  for (const auto& pos : c.pos) {
    double s = 0;
    for (const auto& b : c.bb) {
      const Vec3 d = pos[b.back()] - pos[b.front()];
      s += dot(d, d);
      all.push_back(norm(d));
    }
    per_frame.push_back(s / c.bb.size());
  }
  const auto [m2, e2] = block_mean(per_frame, o.blocks);
  p.value = std::sqrt(m2);
  p.error = std::isnan(e2) ? NaN : e2 / (2 * p.value);
  Series h{"P(R)", "R (Å)", "probability density", {}, {}};
  const double mx = *std::max_element(all.begin(), all.end()) * 1.05 + 1e-9, bw = mx / 40;
  std::vector<double> hist(40, 0.0);
  for (double x : all) hist[std::min<size_t>(39, size_t(x / bw))] += 1;
  for (size_t b = 0; b < 40; ++b) { h.x.push_back((b + 0.5) * bw); h.y.push_back(hist[b] / (all.size() * bw)); }
  p.series.push_back(std::move(h));
  (void)t;
  (void)fr;
  return p;
}

Property cn_prop(const AnalyzeOptions& o, const ChainFrames& c) {
  Property p{"cn", "Characteristic ratio C∞", "", "", NaN, NaN, {}, {}, {}};
  if (c.bb.empty()) { p.notes.push_back("no chains"); return p; }
  size_t lmax = 0;
  for (const auto& b : c.bb) lmax = std::max(lmax, b.size());
  std::vector<double> sum(lmax, 0.0);
  std::vector<double> cnt(lmax, 0.0);
  std::vector<double> ratio_last;   // C_N (whole backbone) per frame
  for (const auto& pos : c.pos) {
    double bsum = 0, bn = 0, frame_r2 = 0, frame_nb2 = 0;
    for (const auto& b : c.bb) {
      for (size_t i = 0; i + 1 < b.size(); ++i) {
        const Vec3 d = pos[b[i + 1]] - pos[b[i]];
        bsum += dot(d, d);
        bn += 1;
      }
      for (size_t i = 0; i < b.size(); ++i)
        for (size_t j = i + 1; j < b.size(); ++j) {
          const Vec3 d = pos[b[j]] - pos[b[i]];
          sum[j - i] += dot(d, d);
          cnt[j - i] += 1;
        }
      const Vec3 e = pos[b.back()] - pos[b.front()];
      frame_r2 += dot(e, e);
      frame_nb2 += double(b.size() - 1);
    }
    const double b2 = bsum / bn;
    ratio_last.push_back(frame_r2 / (frame_nb2 * b2));
  }
  const double b2 = sum[1] / cnt[1];
  Series s{"C_n", "n (backbone bonds)", "C_n = ⟨R²(n)⟩ / (n ⟨b²⟩)", {}, {}};
  for (size_t n = 1; n < lmax; ++n)
    if (cnt[n] > 0) { s.x.push_back(double(n)); s.y.push_back(sum[n] / cnt[n] / (double(n) * b2)); }
  // C∞: C_n = C∞ (1 − a/n) fitted in 1/n over the upper half of n (enough pairs per n: at least a quarter of the chains)
  std::vector<double> xi, yi;
  const double nmax = s.x.empty() ? 0 : s.x.back();
  for (size_t k = 0; k < s.x.size(); ++k)
    if (s.x[k] >= nmax * 0.3 && s.x[k] <= nmax * 0.8) { xi.push_back(1.0 / s.x[k]); yi.push_back(s.y[k]); }
  const auto [cinf, slope] = linfit(xi, yi);
  p.value = cinf;
  const auto [cN, eN] = block_mean(ratio_last, o.blocks);
  p.extra["C_N (whole backbone)"] = cN;
  if (!std::isnan(eN)) p.extra["C_N standard error"] = eN;
  p.extra["⟨b²⟩ (Å²)"] = b2;
  p.extra["backbone bonds per chain"] = nmax;
  p.method = "C_n from all internal distances of the backbones (Theodorou & Suter 1985); C∞ from C_n = C∞(1 − a/n) fitted for n = " +
             std::to_string(int(nmax * 0.3)) + "…" + std::to_string(int(nmax * 0.8));
  if (nmax < 50) p.notes.push_back("short chains (" + std::to_string(int(nmax)) + " backbone bonds): C∞ is an extrapolation");
  (void)slope;
  p.series.push_back(std::move(s));
  return p;
}

Property persistence_prop(const AnalyzeOptions& o, const ChainFrames& c) {
  Property p{"persistence", "Persistence length", "Å", "", NaN, NaN, {}, {}, {}};
  if (c.bb.empty()) { p.notes.push_back("no chains"); return p; }
  size_t lmax = 0;
  for (const auto& b : c.bb) lmax = std::max(lmax, b.size());
  std::vector<double> corr(lmax, 0.0), cnt(lmax, 0.0);
  std::vector<double> flory;   // ⟨R · b1⟩ / |b| per frame, both chain ends
  double bl = 0, bn = 0;
  for (const auto& pos : c.pos) {
    double fs = 0, fn = 0;
    for (const auto& b : c.bb) {
      std::vector<Vec3> u;
      for (size_t i = 0; i + 1 < b.size(); ++i) {
        const Vec3 d = pos[b[i + 1]] - pos[b[i]];
        bl += norm(d);
        bn += 1;
        u.push_back(d * (1.0 / norm(d)));
      }
      for (size_t i = 0; i < u.size(); ++i)
        for (size_t j = i; j < u.size(); ++j) { corr[j - i] += dot(u[i], u[j]); cnt[j - i] += 1; }
      const Vec3 R = pos[b.back()] - pos[b.front()];
      fs += dot(R, u.front()) + dot(R, u.back());
      fn += 2;
    }
    flory.push_back(fs / fn);
  }
  const double lb = bl / bn;
  Series s{"⟨cos θ(k)⟩", "k (bonds apart)", "⟨u_i · u_{i+k}⟩", {}, {}};
  for (size_t k = 0; k < lmax; ++k)
    if (cnt[k] > 0) { s.x.push_back(double(k)); s.y.push_back(corr[k] / cnt[k]); }
  const auto [lf, ef] = block_mean(flory, o.blocks);
  p.value = lf;
  p.error = ef;
  // exponential decay length from ln⟨cos⟩ while it is above 0.1 (vinyl backbones alternate: fit pairs k even)
  std::vector<double> xi, yi;
  for (size_t k = 0; k < s.x.size(); ++k)
    if (s.y[k] > 0.1 && s.x[k] > 0) { xi.push_back(s.x[k]); yi.push_back(std::log(s.y[k])); }
    else if (s.x[k] > 0 && s.y[k] <= 0.1) break;
  const auto fit = linfit(xi, yi);
  if (!std::isnan(fit.second) && fit.second < 0) p.extra["decay length of bond correlations (Å)"] = -lb / fit.second;
  p.extra["mean backbone bond (Å)"] = lb;
  p.method = "Flory: projection of the end-to-end vector on the first (and last) backbone bond, ⟨R·u₁⟩; also the exponential decay "
             "length of backbone bond correlations";
  p.series.push_back(std::move(s));
  return p;
}

// ---------------------------------------------------------------- dynamics

// Multiple-origin average of f(origin, lag) over all origins; lags 0 … maxlag.
template <class F>
std::vector<double> origins_average(size_t nf, size_t maxlag, F&& f) {
  std::vector<double> out(maxlag + 1, 0.0);
  for (size_t lag = 0; lag <= maxlag; ++lag) {
    double s = 0;
    for (size_t o = 0; o + lag < nf; ++o) s += f(o, o + lag);
    out[lag] = s / double(nf - lag);
  }
  return out;
}

struct Dyn {
  std::vector<std::vector<Vec3>> x;    // unwrapped atoms, per frame (drift removed)
  std::vector<std::vector<Vec3>> com;  // molecule centres, per frame
  std::vector<double> t;               // ps, from the first frame
};

Dyn dynamics_data(const Trajectory& t, const std::vector<size_t>& fr, const std::vector<double>& times, std::vector<std::string>& notes) {
  Dyn d;
  d.x = unwrapped_positions(t, fr, notes);
  const System& top = t.topology;
  int nm = 0;
  const auto mol = top.molecules(&nm);
  std::vector<double> m(top.atoms.size());
  double mtot = 0;
  for (size_t i = 0; i < m.size(); ++i) { m[i] = top.mass_of(top.atoms[i]); mtot += m[i]; }
  // remove the drift of the whole system's centre of mass
  Vec3 c0{0, 0, 0};
  for (size_t k = 0; k < d.x.size(); ++k) {
    Vec3 c{0, 0, 0};
    for (size_t i = 0; i < m.size(); ++i) c = c + d.x[k][i] * m[i];
    c = c * (1.0 / mtot);
    if (k == 0) c0 = c;
    const Vec3 shift = c - c0;
    for (auto& p : d.x[k]) p = p - shift;
  }
  // molecule centres from unwrapped atoms: whole molecules need the first frame whole
  System s0 = t.frame(fr[0]);
  if (!s0.unwrapped) make_molecules_whole(s0);
  std::vector<Vec3> offset(m.size());
  for (size_t i = 0; i < m.size(); ++i) offset[i] = s0.atoms[i].pos - d.x[0][i];
  std::vector<double> mm(nm, 0.0);
  for (size_t i = 0; i < m.size(); ++i) mm[mol[i]] += m[i];
  for (const auto& xs : d.x) {
    std::vector<Vec3> c(nm, Vec3{0, 0, 0});
    for (size_t i = 0; i < m.size(); ++i) c[mol[i]] = c[mol[i]] + (xs[i] + offset[i]) * m[i];
    for (int k = 0; k < nm; ++k) c[k] = c[k] * (1.0 / mm[k]);
    d.com.push_back(std::move(c));
  }
  for (double x : times) d.t.push_back(x - times[0]);
  return d;
}

Property msd_prop(const Dyn& d, const AnalyzeOptions& o, Property* diffusion) {
  Property p{"msd", "Mean-square displacement", "Å²", "", NaN, NaN, {}, {}, {}};
  const size_t nf = d.x.size();
  if (nf < 3) { p.notes.push_back("the MSD needs a trajectory (at least 3 frames)"); return p; }
  const size_t maxlag = nf - 1;
  const size_t na = d.x[0].size(), nm = d.com[0].size();
  auto msd_of = [&](const std::vector<std::vector<Vec3>>& X, size_t n) {
    return origins_average(nf, maxlag, [&](size_t a, size_t b) {
      double s = 0;
      for (size_t i = 0; i < n; ++i) { const Vec3 v = X[b][i] - X[a][i]; s += dot(v, v); }
      return s / double(n);
    });
  };
  if (cancelled(o, "MSD", 0.0)) throw Cancel();
  const auto ma = msd_of(d.x, na);
  if (cancelled(o, "MSD", 0.5)) throw Cancel();
  const auto mc = msd_of(d.com, nm);
  Series sa{"atoms", "t (ps)", "MSD (Å²)", {}, {}}, sc{"molecule centres", "t (ps)", "MSD (Å²)", {}, {}};
  for (size_t l = 0; l <= maxlag; ++l) {
    sa.x.push_back(d.t[l]); sa.y.push_back(ma[l]);
    sc.x.push_back(d.t[l]); sc.y.push_back(mc[l]);
  }
  p.value = mc.back();
  p.extra["atom MSD at the longest lag (Å²)"] = ma.back();
  p.method = "all time origins, " + std::to_string(nf) + " frames over " + std::to_string(d.t.back()).substr(0, 7) +
             " ps; the system's centre-of-mass drift removed; value: molecule-centre MSD at the longest lag";
  // Einstein fit on the molecule centres over [fit_from, fit_to] of the run (averages still good there)
  if (diffusion) {
    Property& D = *diffusion;
    D = Property{"diffusion", "Diffusion coefficient", "10⁻⁵ cm²/s", "", NaN, NaN, {}, {}, {}};
    std::vector<double> xi, yi, lx, ly;
    for (size_t l = 1; l <= maxlag; ++l) {
      const double f = d.t[l] / d.t.back();
      if (f >= o.fit_from && f <= o.fit_to) { xi.push_back(d.t[l]); yi.push_back(mc[l]); if (mc[l] > 0) { lx.push_back(std::log(d.t[l])); ly.push_back(std::log(mc[l])); } }
    }
    const auto [a, b] = linfit(xi, yi);
    const auto [la, beta] = linfit(lx, ly);
    (void)a; (void)la;
    if (std::isnan(b)) {
      D.method = "Einstein: slope/6 of the molecule-centre MSD";
      D.notes.push_back("too few frames in the fit window (" + std::to_string(int(o.fit_from * 100)) + "–" + std::to_string(int(o.fit_to * 100)) +
                        " % of the run); store more frames");
    } else {
      D.value = b / 6.0 * 1e-4 * 1e5;   // Å²/ps → cm²/s (1e-4), in units of 1e-5 cm²/s
      D.extra["log-log slope β"] = beta;
      // error: slopes of the two halves of the window
      std::vector<double> x1(xi.begin(), xi.begin() + xi.size() / 2), y1(yi.begin(), yi.begin() + yi.size() / 2);
      std::vector<double> x2(xi.begin() + xi.size() / 2, xi.end()), y2(yi.begin() + yi.size() / 2, yi.end());
      const double b1 = linfit(x1, y1).second, b2 = linfit(x2, y2).second;
      if (!std::isnan(b1) && !std::isnan(b2)) D.error = std::fabs(b1 - b2) / 6.0 * 10.0 / 2;
      D.method = "Einstein: slope/6 of the molecule-centre MSD, fitted from " + std::to_string(int(o.fit_from * 100)) + " to " +
                 std::to_string(int(o.fit_to * 100)) + " % of the run";
      if (beta < 0.9)
        D.notes.push_back("the MSD is not yet linear in time (log-log slope " + std::to_string(beta).substr(0, 4) +
                          " < 0.9): sub-diffusive, so D is an upper bound; run longer");
    }
  }
  p.series.push_back(std::move(sa));
  p.series.push_back(std::move(sc));
  return p;
}

// KWW fit C(t) = exp(−(t/τ)^β) on 0.05 < C < 0.95; returns {τ_KWW, β, integral correlation time}.
std::tuple<double, double, double> kww(const std::vector<double>& t, const std::vector<double>& c) {
  std::vector<double> x, y;
  for (size_t k = 1; k < t.size(); ++k)
    if (c[k] > 0.05 && c[k] < 0.95 && t[k] > 0) { x.push_back(std::log(t[k])); y.push_back(std::log(-std::log(c[k]))); }
  const auto [a, b] = linfit(x, y);
  if (std::isnan(b) || b <= 0) return {NaN, NaN, NaN};
  const double beta = b, tau = std::exp(-a / b);
  return {tau, beta, tau / beta * std::tgamma(1.0 / beta)};
}

Property relaxation_prop(const Trajectory& t, const std::vector<size_t>& fr, const Dyn& d, const AnalyzeOptions& o) {
  Property p{"relaxation", "Relaxation times", "ps", "", NaN, NaN, {}, {}, {}};
  const size_t nf = d.x.size();
  if (nf < 3) { p.notes.push_back("relaxation needs a trajectory"); return p; }
  const System s0 = t.frame(fr[0]);
  const auto bb = backbones(s0);
  if (bb.empty()) { p.notes.push_back("no chains"); return p; }
  // end-to-end vectors (unwrapped atoms keep chains whole) and backbone bond vectors
  auto ree = [&](size_t k, const std::vector<uint32_t>& b) { return d.x[k][b.back()] - d.x[k][b.front()]; };
  const auto cee = origins_average(nf, nf - 1, [&](size_t a, size_t b) {
    double s = 0, n = 0;
    for (const auto& c : bb) { s += dot(ree(a, c), ree(b, c)); n += dot(ree(a, c), ree(a, c)); }
    return s / n;
  });
  if (cancelled(o, "relaxation", 0.5)) throw Cancel();
  const auto cp2 = origins_average(nf, nf - 1, [&](size_t a, size_t b) {
    double s = 0;
    long n = 0;
    for (const auto& c : bb)
      for (size_t i = 0; i + 1 < c.size(); ++i) {
        const Vec3 u = d.x[a][c[i + 1]] - d.x[a][c[i]], v = d.x[b][c[i + 1]] - d.x[b][c[i]];
        const double cs = dot(u, v) / (norm(u) * norm(v));
        s += 1.5 * cs * cs - 0.5;
        ++n;
      }
    return s / double(n);
  });
  Series se{"end-to-end ⟨R(t)·R(0)⟩/⟨R²⟩", "t (ps)", "autocorrelation", {}, {}}, sb{"backbone bonds P₂(t)", "t (ps)", "autocorrelation", {}, {}};
  for (size_t l = 0; l < nf; ++l) { se.x.push_back(d.t[l]); se.y.push_back(cee[l]); sb.x.push_back(d.t[l]); sb.y.push_back(cp2[l]); }
  const auto [tb, bb_beta, tb_int] = kww(sb.x, sb.y);
  const auto [te, be, te_int] = kww(se.x, se.y);
  p.value = tb_int;
  p.extra["segmental τ_KWW (ps)"] = tb;
  p.extra["segmental β"] = bb_beta;
  if (!std::isnan(te_int)) {
    p.extra["end-to-end τ (ps), integral of the KWW fit"] = te_int;
    p.extra["end-to-end β"] = be;
  }
  p.method = "backbone bond P₂ autocorrelation (segmental) and end-to-end vector autocorrelation, all time origins; "
             "KWW fits C = exp(−(t/τ)^β) where 0.05 < C < 0.95; value: segmental correlation time ∫C dt";
  if (se.y.back() > 0.5)
    p.notes.push_back("the end-to-end vector has not decorrelated (C = " + std::to_string(se.y.back()).substr(0, 4) +
                      " at the longest lag): its relaxation time is far longer than the run; the fit extrapolates");
  p.series.push_back(std::move(sb));
  p.series.push_back(std::move(se));
  return p;
}

// ---------------------------------------------------------------- cohesive energy

Property ced_prop(const Trajectory& t, const std::vector<size_t>& fr, const AnalyzeOptions& o, Property* delta) {
  Property p{"ced", "Cohesive energy density", "J/cm³", "", NaN, NaN, {}, {}, {}};
  *delta = Property{"delta", "Hildebrand solubility parameter", "MPa^½", "", NaN, NaN, {}, {}, {}};
  if (!o.ff) {
    p.notes.push_back("needs a force field (assign one in Field, or --ff)");
    delta->notes = p.notes;
    return p;
  }
  if (!t.topology.cell.valid()) { p.notes.push_back("needs a periodic cell"); delta->notes = p.notes; return p; }
  const ForceField& ff = *o.ff;
  int nm = 0;
  const auto mol = t.topology.molecules(&nm);
  std::vector<std::vector<uint32_t>> members(nm);
  for (uint32_t i = 0; i < mol.size(); ++i) members[mol[i]].push_back(i);
  Evaluator bulk(ff, o.energy);
  EnergyOptions iso_e = o.energy;
  iso_e.tail = false;
  std::vector<std::unique_ptr<Evaluator>> iso;
  std::vector<ForceField> iso_ff;
  iso_ff.reserve(nm);
  for (int k = 0; k < nm; ++k) iso_ff.push_back(subset_forcefield(ff, members[k]));
  for (int k = 0; k < nm; ++k) iso.push_back(std::make_unique<Evaluator>(iso_ff[k], iso_e));
  std::vector<double> ced, coh;
  Series s{"CED", "time (ps)", "CED (J/cm³)", {}, {}};
  const auto times = frame_times(t, o);
  for (size_t q = 0; q < fr.size(); ++q) {
    if (cancelled(o, "cohesive energy", double(q) / fr.size())) throw Cancel();
    System f = t.frame(fr[q]);
    if (!f.unwrapped) make_molecules_whole(f);
    std::vector<double> x, g;
    for (const auto& a : f.atoms) x.insert(x.end(), a.pos.begin(), a.pos.end());
    const double eb = bulk.compute(x, f.cell, g).total();
    double ei = 0;
    for (int k = 0; k < nm; ++k) {
      std::vector<double> xm;
      for (uint32_t i : members[k]) xm.insert(xm.end(), f.atoms[i].pos.begin(), f.atoms[i].pos.end());
      ei += iso[k]->compute(xm, Cell{}, g).total();
    }
    const double ecoh = ei - eb;   // kcal/mol per cell, positive when the bulk is bound
    const double v = f.cell.volume();
    const double c = ecoh * 4184.0 / kNA / (v * 1e-24);
    ced.push_back(c);
    coh.push_back(ecoh);
    s.x.push_back(times[fr[q]] - times[fr[0]]);
    s.y.push_back(c);
  }
  std::tie(p.value, p.error) = block_mean(ced, o.blocks);
  p.extra["cohesive energy per molecule (kcal/mol)"] = std::accumulate(coh.begin(), coh.end(), 0.0) / coh.size() / nm;
  p.method = "(Σ E_isolated molecule − E_bulk) / V with " + ff.name + " (same cut-off and electrostatics; molecules in vacuum), averaged over " +
             std::to_string(fr.size()) + " frames";
  p.series.push_back(std::move(s));
  delta->value = p.value > 0 ? std::sqrt(p.value) : NaN;
  delta->error = p.value > 0 && !std::isnan(p.error) ? p.error / (2 * delta->value) : NaN;
  delta->method = "δ = √CED";
  return p;
}

// ---------------------------------------------------------------- free volume

struct FreeGrid {
  int n[3] = {0, 0, 0};
  std::vector<float> dist;   // distance from each grid point to the nearest van der Waals surface (Å); < 0 inside
  double cell_volume = 0;
};

// Distance to the nearest atom surface at each grid point, for fractional coordinates on an n0 × n1 × n2 grid.
FreeGrid free_grid(const System& s, double spacing, double search) {
  FreeGrid G;
  const Cell& c = s.cell;
  G.cell_volume = c.volume();
  const double la = norm(c.a), lb = norm(c.b), lc = norm(c.c);
  G.n[0] = std::max(4, int(std::ceil(la / spacing)));
  G.n[1] = std::max(4, int(std::ceil(lb / spacing)));
  G.n[2] = std::max(4, int(std::ceil(lc / spacing)));
  const size_t np = size_t(G.n[0]) * G.n[1] * G.n[2];
  G.dist.assign(np, float(search));
  // bins of atoms in fractional space, bin width ≥ search
  int nb[3];
  const double w[3] = {G.cell_volume / norm(cross(c.b, c.c)), G.cell_volume / norm(cross(c.c, c.a)), G.cell_volume / norm(cross(c.a, c.b))};
  for (int k = 0; k < 3; ++k) nb[k] = std::max(1, int(w[k] / search));
  std::vector<std::vector<uint32_t>> bins(size_t(nb[0]) * nb[1] * nb[2]);
  std::vector<Vec3> frac(s.atoms.size());
  std::vector<double> rad(s.atoms.size());
  for (uint32_t i = 0; i < s.atoms.size(); ++i) {
    Vec3 f = c.to_fractional(s.atoms[i].pos);
    for (int k = 0; k < 3; ++k) f[k] -= std::floor(f[k]);
    frac[i] = f;
    rad[i] = element(s.atoms[i].element).vdw;
    int b[3];
    for (int k = 0; k < 3; ++k) b[k] = std::min(nb[k] - 1, int(f[k] * nb[k]));
    bins[(size_t(b[0]) * nb[1] + b[1]) * nb[2] + b[2]].push_back(i);
  }
  for (int i0 = 0; i0 < G.n[0]; ++i0)
    for (int i1 = 0; i1 < G.n[1]; ++i1)
      for (int i2 = 0; i2 < G.n[2]; ++i2) {
        const Vec3 f{(i0 + 0.5) / G.n[0], (i1 + 0.5) / G.n[1], (i2 + 0.5) / G.n[2]};
        const Vec3 r = c.to_cartesian(f);
        int b[3];
        for (int k = 0; k < 3; ++k) b[k] = std::min(nb[k] - 1, int(f[k] * nb[k]));
        double best = search;
        for (int d0 = -1; d0 <= 1; ++d0)
          for (int d1 = -1; d1 <= 1; ++d1)
            for (int d2 = -1; d2 <= 1; ++d2) {
              const int q0 = ((b[0] + d0) % nb[0] + nb[0]) % nb[0], q1 = ((b[1] + d1) % nb[1] + nb[1]) % nb[1], q2 = ((b[2] + d2) % nb[2] + nb[2]) % nb[2];
              if ((nb[0] < 3 && d0 > 0 && nb[0] + d0 > 2) || (nb[1] < 3 && d1 > 0 && nb[1] + d1 > 2) || (nb[2] < 3 && d2 > 0 && nb[2] + d2 > 2)) continue;
              for (uint32_t a : bins[(size_t(q0) * nb[1] + q1) * nb[2] + q2]) {
                const double dd = norm(c.minimum_image(s.atoms[a].pos - r)) - rad[a];
                if (dd < best) best = dd;
              }
            }
        G.dist[(size_t(i0) * G.n[1] + i1) * G.n[2] + i2] = float(best);
      }
  (void)frac;
  return G;
}

std::vector<size_t> spread_frames(const std::vector<size_t>& fr, size_t n) {
  if (fr.size() <= n) return fr;
  std::vector<size_t> out;
  for (size_t k = 0; k < n; ++k) out.push_back(fr[k * (fr.size() - 1) / (n - 1)]);
  return out;
}

Property ffv_prop(const Trajectory& t, const std::vector<size_t>& fr0, const AnalyzeOptions& o) {
  Property p{"ffv", "Free volume", "", "", NaN, NaN, {}, {}, {}};
  if (!t.topology.cell.valid()) { p.notes.push_back("needs a periodic cell"); return p; }
  const auto fr = spread_frames(fr0, 5);
  std::vector<double> acc, occ;
  Series s{"accessible fraction", "probe radius (Å)", "fraction of the cell", {}, {}};
  const std::vector<double> probes = {0.0, 0.5, 1.0, 1.28, 1.4, 1.82};   // point, He (1.28 Å, Bondi), H2O-sized, N2 / CO2 kinetic radii
  std::vector<double> sums(probes.size(), 0.0);
  for (size_t q = 0; q < fr.size(); ++q) {
    if (cancelled(o, "free volume", double(q) / fr.size())) throw Cancel();
    const System f = t.frame(fr[q]);
    const FreeGrid G = free_grid(f, o.grid, 4.0);
    const double np = double(G.dist.size());
    size_t free_p = 0, free0 = 0;
    for (float d : G.dist) { free_p += d > o.probe; free0 += d > 0; }
    acc.push_back(free_p / np);
    occ.push_back(1.0 - free0 / np);
    for (size_t k = 0; k < probes.size(); ++k) {
      size_t n = 0;
      for (float d : G.dist) n += d > probes[k];
      sums[k] += n / np;
    }
  }
  for (size_t k = 0; k < probes.size(); ++k) { s.x.push_back(probes[k]); s.y.push_back(sums[k] / fr.size()); }
  std::tie(p.value, p.error) = block_mean(acc, o.blocks);
  const double vw = std::accumulate(occ.begin(), occ.end(), 0.0) / occ.size();
  p.extra["van der Waals occupied fraction"] = vw;
  p.extra["Bondi FFV = 1 − 1.3 V_w/V"] = 1 - 1.3 * vw;
  p.method = "probe insertion on a " + std::to_string(o.grid).substr(0, 4) + " Å grid, Bondi van der Waals radii, probe radius " +
             std::to_string(o.probe).substr(0, 4) + " Å, " + std::to_string(fr.size()) + " frames; value: fraction of the cell a probe centre can reach";
  p.series.push_back(std::move(s));
  return p;
}

Property psd_prop(const Trajectory& t, const std::vector<size_t>& fr0, const AnalyzeOptions& o) {
  Property p{"psd", "Pore size distribution", "Å", "", NaN, NaN, {}, {}, {}};
  if (!t.topology.cell.valid()) { p.notes.push_back("needs a periodic cell"); return p; }
  const auto fr = spread_frames(fr0, 3);
  std::vector<double> hist;
  const double bw = 0.2;
  double mean_d = 0, npts = 0, largest = 0;
  for (size_t q = 0; q < fr.size(); ++q) {
    if (cancelled(o, "pore sizes", double(q) / fr.size())) throw Cancel();
    const System f = t.frame(fr[q]);
    const FreeGrid G = free_grid(f, o.grid, 6.0);
    const int n0 = G.n[0], n1 = G.n[1], n2 = G.n[2];
    // the largest sphere that contains each free point without overlapping an atom: centres in descending radius
    std::vector<size_t> order;
    for (size_t k = 0; k < G.dist.size(); ++k) if (G.dist[k] > 0) order.push_back(k);
    std::sort(order.begin(), order.end(), [&](size_t a, size_t b) { return G.dist[a] > G.dist[b]; });
    std::vector<float> best(G.dist.size(), 0.f);
    const Cell& c = f.cell;
    const double h[3] = {norm(c.a) / n0, norm(c.b) / n1, norm(c.c) / n2};
    for (size_t k : order) {
      const float r = G.dist[k];
      if (best[k] >= r) continue;   // a larger sphere already covers the centre: its own sphere lies mostly inside it
      const int i0 = int(k / (size_t(n1) * n2)), i1 = int(k / n2 % n1), i2 = int(k % n2);
      const int e0 = int(r / h[0]) + 1, e1 = int(r / h[1]) + 1, e2 = int(r / h[2]) + 1;
      const Vec3 rc = c.to_cartesian({(i0 + 0.5) / n0, (i1 + 0.5) / n1, (i2 + 0.5) / n2});
      for (int d0 = -e0; d0 <= e0; ++d0)
        for (int d1 = -e1; d1 <= e1; ++d1)
          for (int d2 = -e2; d2 <= e2; ++d2) {
            const int j0 = ((i0 + d0) % n0 + n0) % n0, j1 = ((i1 + d1) % n1 + n1) % n1, j2 = ((i2 + d2) % n2 + n2) % n2;
            const size_t j = (size_t(j0) * n1 + j1) * n2 + j2;
            if (best[j] >= r || G.dist[j] <= 0) continue;
            const Vec3 rj = c.to_cartesian({(i0 + d0 + 0.5) / n0, (i1 + d1 + 0.5) / n1, (i2 + d2 + 0.5) / n2});
            if (norm(rj - rc) <= r) best[j] = r;
          }
    }
    // the largest cavity off the grid: pattern search on the exact distance field from the best grid points
    for (size_t m = 0; m < std::min<size_t>(order.size(), 30); ++m) {
      const size_t k = order[m];
      const int i0 = int(k / (size_t(n1) * n2)), i1 = int(k / n2 % n1), i2 = int(k % n2);
      Vec3 x = c.to_cartesian({(i0 + 0.5) / n0, (i1 + 0.5) / n1, (i2 + 0.5) / n2});
      auto surface = [&](const Vec3& r) {
        double d = std::numeric_limits<double>::max();
        for (const auto& a : f.atoms) d = std::min(d, norm(c.minimum_image(a.pos - r)) - element(a.element).vdw);
        return d;
      };
      double fx = surface(x);
      for (double step = o.grid; step > 1e-3;) {
        bool moved = false;
        for (int ax = 0; ax < 3 && !moved; ++ax)
          for (double sgn : {-1.0, 1.0}) {
            Vec3 y = x;
            y[ax] += sgn * step;
            const double fy = surface(y);
            if (fy > fx) { x = y; fx = fy; moved = true; break; }
          }
        if (!moved) step *= 0.5;
      }
      largest = std::max(largest, 2 * fx);
    }
    for (size_t k = 0; k < best.size(); ++k)
      if (G.dist[k] > 0) {
        const double dia = 2.0 * best[k];
        const size_t b = size_t(dia / bw);
        if (hist.size() <= b) hist.resize(b + 1, 0.0);
        hist[b] += 1;
        mean_d += dia;
        npts += 1;
      }
  }
  Series s{"PSD", "pore diameter (Å)", "fraction of free volume per Å", {}, {}};
  for (size_t b = 0; b < hist.size(); ++b) { s.x.push_back((b + 0.5) * bw); s.y.push_back(npts > 0 ? hist[b] / (npts * bw) : 0); }
  if (npts > 0) {
    p.value = mean_d / npts;
    size_t mode = 0;
    for (size_t b = 1; b < hist.size(); ++b) if (hist[b] > hist[mode]) mode = b;
    p.extra["most frequent diameter (Å)"] = (mode + 0.5) * bw;
    p.extra["largest diameter (Å)"] = largest;   // refined off the grid
  }
  p.method = "each free point takes the diameter of the largest atom-free sphere containing it (Gelb & Gubbins 1999), Bondi radii, " +
             std::to_string(o.grid).substr(0, 4) + " Å grid, " + std::to_string(fr.size()) + " frames; value: volume-weighted mean diameter";
  p.series.push_back(std::move(s));
  return p;
}

}  // namespace

// ---------------------------------------------------------------- scattering data

double neutron_b(int z) {
  // coherent scattering lengths, fm (NIST Center for Neutron Research, natural isotopic abundance)
  static const std::map<int, double> b = {{1, -3.739}, {2, 3.26}, {3, -1.90}, {5, 5.30}, {6, 6.646}, {7, 9.36}, {8, 5.803}, {9, 5.654},
                                          {11, 3.63}, {12, 5.375}, {13, 3.449}, {14, 4.1491}, {15, 5.13}, {16, 2.847}, {17, 9.577}, {19, 3.67},
                                          {20, 4.70}, {22, -3.438}, {26, 9.45}, {29, 7.718}, {30, 5.68}, {35, 6.795}, {47, 5.922}, {53, 5.28},
                                          {78, 9.60}, {79, 7.63}};
  auto it = b.find(z);
  if (it == b.end()) throw std::invalid_argument(std::string("no neutron scattering length for ") + element(z).symbol);
  return it->second;
}

double xray_f(int z, double q) {
  // Cromer–Mann coefficients (International Tables for Crystallography, Vol. C, Table 6.1.1.4): a1 b1 a2 b2 a3 b3 a4 b4 c
  static const std::map<int, std::array<double, 9>> cm = {
      {1, {0.489918, 20.6593, 0.262003, 7.74039, 0.196767, 49.5519, 0.049879, 2.20159, 0.001305}},
      {6, {2.31000, 20.8439, 1.02000, 10.2075, 1.58860, 0.568700, 0.865000, 51.6512, 0.215600}},
      {7, {12.2126, 0.005700, 3.13220, 9.89330, 2.01250, 28.9975, 1.16630, 0.582600, -11.529}},
      {8, {3.04850, 13.2771, 2.28680, 5.70110, 1.54630, 0.323900, 0.867000, 32.9089, 0.250800}},
      {9, {3.53920, 10.2825, 2.64120, 4.29440, 1.51700, 0.261500, 1.02430, 26.1476, 0.277600}},
      {11, {4.76260, 3.28500, 3.17360, 8.84220, 1.26740, 0.313600, 1.11280, 129.424, 0.676000}},
      {14, {6.29150, 2.43860, 3.03530, 32.3337, 1.98910, 0.678500, 1.54100, 81.6937, 1.14070}},
      {15, {6.43450, 1.90670, 4.17910, 27.1570, 1.78000, 0.526000, 1.49080, 68.1645, 1.11490}},
      {16, {6.90530, 1.46790, 5.20340, 22.2151, 1.43790, 0.253600, 1.58630, 56.1720, 0.866900}},
      {17, {11.4604, 0.010400, 7.19640, 1.16620, 6.25560, 18.5194, 1.64550, 47.7784, -9.5574}},
      {35, {17.1789, 2.17230, 5.23580, 16.5796, 5.63770, 0.260900, 3.98510, 41.4328, 2.95570}},
      {53, {20.1472, 4.34700, 18.9949, 0.381400, 7.51380, 27.7660, 2.27350, 66.8776, 4.07120}}};
  auto it = cm.find(z);
  if (it == cm.end()) throw std::invalid_argument(std::string("no X-ray form factor for ") + element(z).symbol);
  const auto& a = it->second;
  const double s2 = (q / (4 * kPi)) * (q / (4 * kPi));   // (sin θ / λ)²
  return a[0] * std::exp(-a[1] * s2) + a[2] * std::exp(-a[3] * s2) + a[4] * std::exp(-a[5] * s2) + a[6] * std::exp(-a[7] * s2) + a[8];
}

// ---------------------------------------------------------------- force-field subset

ForceField subset_forcefield(const ForceField& ff, const std::vector<uint32_t>& atoms) {
  std::vector<int64_t> map(ff.charge.size(), -1);
  for (size_t k = 0; k < atoms.size(); ++k) map[atoms[k]] = int64_t(k);
  auto in = [&](std::initializer_list<uint32_t> v) {
    for (uint32_t x : v)
      if (map[x] < 0) return false;
    return true;
  };
  auto m = [&](uint32_t x) { return uint32_t(map[x]); };
  ForceField s;
  s.name = ff.name;
  s.type_names = ff.type_names;
  s.lj = ff.lj;
  s.lj14_types = ff.lj14_types;
  s.pair_func = ff.pair_func;
  s.pair_override = ff.pair_override;
  s.pair_form = ff.pair_form;
  s.mixing = ff.mixing;
  s.lj14 = ff.lj14;
  s.coul14 = ff.coul14;
  for (uint32_t a : atoms) {
    s.atom_type.push_back(ff.atom_type[a]);
    s.why.push_back(a < ff.why.size() ? ff.why[a] : "");
    s.type_index.push_back(ff.type_index[a]);
    s.charge.push_back(ff.charge[a]);
    s.mass.push_back(ff.mass[a]);
  }
  for (const auto& t : ff.bonds) if (in({t.i, t.j})) s.bonds.push_back({m(t.i), m(t.j), t.k, t.r0});
  for (const auto& t : ff.angles) if (in({t.i, t.j, t.k})) s.angles.push_back({m(t.i), m(t.j), m(t.k), t.kt, t.theta0});
  for (const auto& t : ff.dihedrals) if (in({t.i, t.j, t.k, t.l})) s.dihedrals.push_back({m(t.i), m(t.j), m(t.k), m(t.l), t.v, t.n, t.delta});
  for (const auto& t : ff.impropers) if (in({t.i, t.j, t.k, t.l})) s.impropers.push_back({m(t.i), m(t.j), m(t.k), m(t.l), t.v, t.n, t.delta});
  for (const auto& t : ff.impropers_harmonic) if (in({t.i, t.j, t.k, t.l})) s.impropers_harmonic.push_back({m(t.i), m(t.j), m(t.k), m(t.l), t.k2, t.chi0});
  for (const auto& t : ff.inversions) if (in({t.c, t.a, t.b, t.d})) s.inversions.push_back({m(t.c), m(t.a), m(t.b), m(t.d), t.kw, t.w0, t.form});
  for (const auto& t : ff.bonds_x) if (in({t.i, t.j})) s.bonds_x.push_back({m(t.i), m(t.j), t.form, t.a, t.b, t.c});
  for (const auto& t : ff.angles_x) if (in({t.i, t.j, t.k})) s.angles_x.push_back({m(t.i), m(t.j), m(t.k), t.form, t.a, t.b});
  for (const auto& t : ff.urey_bradley) if (in({t.i, t.k})) s.urey_bradley.push_back({m(t.i), m(t.k), t.kub, t.r0});
  for (auto t : ff.bonds2) if (in({t.i, t.j})) { t.i = m(t.i); t.j = m(t.j); s.bonds2.push_back(t); }
  for (auto t : ff.angles2) if (in({t.i, t.j, t.k})) { t.i = m(t.i); t.j = m(t.j); t.k = m(t.k); s.angles2.push_back(t); }
  for (auto t : ff.dihedrals2) if (in({t.i, t.j, t.k, t.l})) { t.i = m(t.i); t.j = m(t.j); t.k = m(t.k); t.l = m(t.l); s.dihedrals2.push_back(t); }
  for (auto t : ff.impropers2) if (in({t.i, t.j, t.k, t.l})) { t.i = m(t.i); t.j = m(t.j); t.k = m(t.k); t.l = m(t.l); s.impropers2.push_back(t); }
  for (const auto& p : ff.pairs14) if (in({p[0], p[1]})) s.pairs14.push_back({m(p[0]), m(p[1])});
  s.excluded.resize(atoms.size());
  for (size_t k = 0; k < atoms.size(); ++k)
    if (atoms[k] < ff.excluded.size())
      for (uint32_t j : ff.excluded[atoms[k]])
        if (map[j] >= 0) s.excluded[k].push_back(uint32_t(map[j]));
  for (auto& e : s.excluded) std::sort(e.begin(), e.end());
  return s;
}

// ---------------------------------------------------------------- JSON

std::string properties_json(const std::vector<Property>& props) {
  auto num = [](double v) { return std::isfinite(v) ? Json(v) : Json(); };
  Json arr = Json::array();
  for (const auto& p : props) {
    Json o = Json::object();
    o["id"] = p.id;
    o["name"] = p.name;
    o["value"] = num(p.value);
    o["error"] = num(p.error);
    o["unit"] = p.unit;
    o["method"] = p.method;
    Json ex = Json::object();
    for (const auto& [k, v] : p.extra) ex[k] = num(v);
    o["extra"] = ex;
    Json notes = Json::array();
    for (const auto& n : p.notes) notes.push_back(n);
    o["notes"] = notes;
    Json ser = Json::array();
    for (const auto& s : p.series) {
      Json j = Json::object();
      j["label"] = s.label;
      j["x_label"] = s.x_label;
      j["y_label"] = s.y_label;
      Json x = Json::array(), y = Json::array();
      for (double v : s.x) x.push_back(num(v));
      for (double v : s.y) y.push_back(num(v));
      j["x"] = x;
      j["y"] = y;
      ser.push_back(j);
    }
    o["series"] = ser;
    arr.push_back(o);
  }
  return arr.dump(0);
}

// ---------------------------------------------------------------- driver

std::vector<size_t> analysis_frames(const Trajectory& t, const AnalyzeOptions& o) {
  const long nf = long(t.frames());
  const long last = o.last < 0 || o.last >= nf ? nf - 1 : o.last;
  std::vector<size_t> fr;
  for (long k = std::max(0L, o.first); k <= last; k += std::max(1L, o.stride)) fr.push_back(size_t(k));
  if (fr.empty() && nf > 0) fr.push_back(size_t(nf - 1));
  return fr;
}

std::vector<double> frame_times(const Trajectory& t, const AnalyzeOptions& o) {
  std::vector<double> tt(t.frames());
  for (size_t k = 0; k < tt.size(); ++k)
    tt[k] = o.frame_ps > 0 ? k * o.frame_ps : (k < t.timesteps.size() ? double(t.timesteps[k]) * o.timestep_fs * 1e-3 : double(k));
  return tt;
}

std::vector<Property> analyze(const Trajectory& t, const std::vector<std::string>& ids, const AnalyzeOptions& o) {
  static const std::set<std::string> known = {"density", "rdf", "sq", "xray", "neutron", "rg", "ree", "cn", "persistence", "msd", "diffusion",
                                              "relaxation", "ced", "delta", "ffv", "psd", "cij_fluct"};
  for (const auto& id : ids)
    if (!known.count(id)) throw std::invalid_argument("unknown property '" + id + "'");
  const auto fr = analysis_frames(t, o);
  const auto all_times = frame_times(t, o);
  std::vector<double> times;
  for (size_t k : fr) times.push_back(all_times[k]);
  auto want = [&](const char* id) { return std::find(ids.begin(), ids.end(), id) != ids.end(); };
  std::vector<Property> out;
  auto run = [&](auto&& f) {
    try { f(); } catch (const Cancel&) { throw; }
  };
  run([&] {
    if (want("density")) out.push_back(density_prop(t, fr, times, o));
    if (want("rdf")) out.push_back(rdf_prop(t, fr, o));
    for (const char* k : {"sq", "xray", "neutron"})
      if (want(k)) out.push_back(scattering_prop(t, fr, o, k));
    if (want("rg")) out.push_back(rg_prop(t, fr, o));
    if (want("ree") || want("cn") || want("persistence")) {
      const ChainFrames c = chain_frames(t, fr);
      if (want("ree")) out.push_back(ree_prop(t, fr, o, c));
      if (want("cn")) out.push_back(cn_prop(o, c));
      if (want("persistence")) out.push_back(persistence_prop(o, c));
    }
    if (want("msd") || want("diffusion") || want("relaxation")) {
      std::vector<std::string> notes;
      if (fr.size() < 3) {
        for (const char* k : {"msd", "diffusion", "relaxation"})
          if (want(k)) {
            Property p;
            p.id = k;
            p.name = k;
            p.notes.push_back("needs a trajectory of at least 3 frames");
            out.push_back(p);
          }
      } else {
        const Dyn d = dynamics_data(t, fr, times, notes);
        if (want("msd") || want("diffusion")) {
          Property D;
          Property M = msd_prop(d, o, &D);
          for (const auto& n : notes) M.notes.push_back(n);
          if (want("msd")) out.push_back(M);
          if (want("diffusion")) out.push_back(D);
        }
        if (want("relaxation")) out.push_back(relaxation_prop(t, fr, d, o));
      }
    }
    if (want("ced") || want("delta")) {
      Property delta;
      Property c = ced_prop(t, fr, o, &delta);
      if (want("ced")) out.push_back(c);
      if (want("delta")) out.push_back(delta);
    }
    if (want("cij_fluct")) {
      if (!o.ff || o.temperature <= 0 || fr.size() < 2 || !t.topology.cell.valid()) {
        Property p;
        p.id = "cij_fluct";
        p.name = "Elastic constants (fluctuations)";
        p.unit = "GPa";
        p.notes.push_back(!o.ff ? "needs the force field of the trajectory"
                          : o.temperature <= 0 ? "give the temperature of the NVT run the frames come from"
                          : fr.size() < 2 ? "needs many frames of an NVT run" : "needs a periodic cell");
        out.push_back(p);
      } else {
        FluctuationOptions fo;
        fo.ff = o.ff;
        fo.energy = o.energy;
        fo.temperature = o.temperature;
        fo.blocks = o.blocks;
        if (o.progress) fo.progress = [&](const std::string& w, double f) { if (!o.progress(w, f)) throw Cancel(); return true; };
        for (auto& p : elastic_properties(fluctuation_elastic(t, fr, fo), "_fluct")) out.push_back(std::move(p));
      }
    }
    if (want("ffv")) out.push_back(ffv_prop(t, fr, o));
    if (want("psd")) out.push_back(psd_prop(t, fr, o));
  });
  return out;
}

}  // namespace caps
