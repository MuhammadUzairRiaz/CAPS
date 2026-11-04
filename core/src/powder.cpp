// CAPS powder pattern from the cell's reciprocal lattice (see powder.hpp).
#include "caps/powder.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdio>
#include <map>
#include <stdexcept>
#include <tuple>

#include "caps/properties.hpp"
#include "parallel.hpp"

namespace caps {

namespace {
constexpr double kPi = 3.14159265358979323846;
constexpr double kDeg = 180.0 / kPi;
std::string fmt(double v, int prec = 4) {
  char b[48];
  std::snprintf(b, sizeof b, "%.*g", prec, v);
  return b;
}
}  // namespace

std::array<double, 2> anode_kalpha(const std::string& anode) {
  // Kα1, Kα2 (Å) as in pymatgen's diffraction module (its table of X-ray wavelengths); Cu as the instrument defaults
  // (Kα1 1.54056, Kα2 1.54439)
  static const std::map<std::string, std::array<double, 2>> lines = {
      {"Cu", {1.54056, 1.54439}}, {"Mo", {0.70930, 0.71359}}, {"Co", {1.78896, 1.79285}},
      {"Fe", {1.93604, 1.93998}}, {"Cr", {2.28970, 2.29361}}, {"Ag", {0.559421, 0.563813}}};
  const auto it = lines.find(anode);
  if (it == lines.end()) throw std::invalid_argument("powder: unknown anode '" + anode + "' (Cu, Mo, Co, Fe, Cr, Ag, or give a wavelength)");
  return it->second;
}

PowderReport powder_pattern(const Trajectory& t, const std::vector<size_t>& frames_in, const PowderOptions& o) {
  const auto t0 = std::chrono::steady_clock::now();
  PowderReport rep;
  const bool neutron = o.radiation == "neutron";
  if (!neutron && o.radiation != "xray") throw std::invalid_argument("powder: radiation is xray or neutron");
  if (!(o.two_theta_max > o.two_theta_min) || o.two_theta_min < 0 || o.two_theta_max >= 180 || o.step <= 0)
    throw std::invalid_argument("powder: a 2θ range inside 0–180° and a positive step");
  // the lines
  if (o.wavelength > 0) rep.wavelengths = {o.wavelength}, rep.weights = {1.0};
  else {
    const auto ka = anode_kalpha(o.anode);
    rep.wavelengths = {ka[0]}, rep.weights = {1.0};
    if (o.kalpha2 && !neutron && o.ka2_ratio > 0) rep.wavelengths.push_back(ka[1]), rep.weights.push_back(o.ka2_ratio);
  }
  const double lam_min = *std::min_element(rep.wavelengths.begin(), rep.wavelengths.end());
  std::vector<size_t> frames = frames_in;
  if (frames.empty()) frames.push_back(0);
  const size_t nf_have = std::max<size_t>(1, t.frames());
  for (size_t f : frames)
    if (f >= nf_have) throw std::invalid_argument("powder: no frame " + std::to_string(f + 1));
  auto frame_of = [&](size_t k) { return t.frames() == 0 ? t.topology : t.frame(k); };
  const System s0 = frame_of(frames[0]);
  if (!s0.cell.valid()) throw std::invalid_argument("powder: the structure needs a periodic cell");
  const size_t n = s0.atoms.size();
  if (n == 0) throw std::invalid_argument("powder: no atoms");
  // elements and their weights
  std::vector<int> elem_of(n), elems;
  for (size_t j = 0; j < n; ++j) {
    const int z = s0.atoms[j].element;
    auto it = std::find(elems.begin(), elems.end(), z);
    if (it == elems.end()) elems.push_back(z), it = elems.end() - 1;
    elem_of[j] = int(it - elems.begin());
  }
  for (int z : elems) (void)scatter_w(neutron ? "neutron" : "xray", z, 1.0);   // an element with no data throws here
  // the reflections: half the reciprocal lattice up to q_max (a little beyond, for cells that change between frames)
  const double th_max = std::min(89.999, 0.5 * (o.two_theta_max + std::fabs(o.zero_shift))) / kDeg;
  const double qmax = 4 * kPi * std::sin(th_max) / lam_min;
  auto recip = [](const Cell& c) {
    const double vol = c.volume();
    return std::array<Vec3, 3>{cross(c.b, c.c) * (2 * kPi / vol), cross(c.c, c.a) * (2 * kPi / vol), cross(c.a, c.b) * (2 * kPi / vol)};
  };
  const auto B0 = recip(s0.cell);
  const int hm = int(qmax * norm(s0.cell.a) / (2 * kPi)) + 1, km = int(qmax * norm(s0.cell.b) / (2 * kPi)) + 1, lm = int(qmax * norm(s0.cell.c) / (2 * kPi)) + 1;
  struct HK { int h, k; std::vector<int> ls; };
  std::vector<HK> hks;
  size_t nvec = 0;
  for (int h = 0; h <= hm; ++h)
    for (int k = (h == 0 ? 0 : -km); k <= km; ++k) {
      HK g{h, k, {}};
      for (int l = (h == 0 && k == 0 ? 1 : -lm); l <= lm; ++l)
        if (norm(B0[0] * double(h) + B0[1] * double(k) + B0[2] * double(l)) <= qmax * 1.02) g.ls.push_back(l);
      if (!g.ls.empty()) nvec += g.ls.size(), hks.push_back(std::move(g));
    }
  if (double(nvec) * double(n) * double(frames.size()) > 4e10)
    throw std::invalid_argument("powder: " + std::to_string(nvec) + " reflections × " + std::to_string(n) + " atoms × " + std::to_string(frames.size()) +
                                " frames is too much: lower 2θ max or use fewer frames");
  // per reflection (h k l in hks order): Σ |F|², Σ |G|, Σ G (for the orientation)
  std::vector<size_t> offset(hks.size() + 1, 0);
  for (size_t g = 0; g < hks.size(); ++g) offset[g + 1] = offset[g] + hks[g].ls.size();
  std::vector<double> f2(nvec, 0.0), gsum(nvec, 0.0);
  std::vector<Vec3> gvec(nvec, Vec3{0, 0, 0});
  const int nth = o.threads > 0 ? o.threads : default_threads();
  ThreadPool pool(nth);
  using cd = std::complex<double>;
  const double biso = std::max(0.0, o.b_iso);
  for (size_t fi = 0; fi < frames.size(); ++fi) {
    if (o.progress && !o.progress(double(fi) / double(frames.size()))) throw std::invalid_argument("powder: stopped");
    const System s = fi == 0 ? s0 : frame_of(frames[fi]);
    if (s.atoms.size() != n) throw std::invalid_argument("powder: the frames hold different atoms");
    const auto B = recip(s.cell);
    auto powers = [&](const Vec3& bv, int m) {   // e^{i h (b·r)} for every atom, h = −m … m, by recurrence
      std::vector<cd> e(size_t(2 * m + 1) * n);
      for (size_t j = 0; j < n; ++j) {
        const cd e1 = std::polar(1.0, dot(bv, s.atoms[j].pos)), e1c = std::conj(e1);
        e[size_t(m) * n + j] = 1.0;
        for (int h = 1; h <= m; ++h) {
          e[size_t(m + h) * n + j] = e[size_t(m + h - 1) * n + j] * e1;
          e[size_t(m - h) * n + j] = e[size_t(m - h + 1) * n + j] * e1c;
        }
      }
      return e;
    };
    const auto E1 = powers(B[0], hm), E2 = powers(B[1], km), E3 = powers(B[2], lm);
    pool.run(hks.size(), [&](int, size_t a, size_t b) {
      std::vector<cd> ab(n);
      std::vector<double> w(elems.size());
      for (size_t g = a; g < b; ++g) {
        const auto& hk = hks[g];
        for (size_t j = 0; j < n; ++j) ab[j] = E1[size_t(hm + hk.h) * n + j] * E2[size_t(km + hk.k) * n + j];
        for (size_t li = 0; li < hk.ls.size(); ++li) {
          const int l = hk.ls[li];
          const Vec3 G = B[0] * double(hk.h) + B[1] * double(hk.k) + B[2] * double(l);
          const double q = norm(G), s2 = (q / (4 * kPi)) * (q / (4 * kPi));
          const double dw = std::exp(-biso * s2);
          for (size_t e = 0; e < elems.size(); ++e) w[e] = scatter_w(neutron ? "neutron" : "xray", elems[e], q) * dw;
          cd A = 0;
          for (size_t j = 0; j < n; ++j) A += w[size_t(elem_of[j])] * ab[j] * E3[size_t(lm + l) * n + j];
          const size_t idx = offset[g] + li;
          f2[idx] += std::norm(A), gsum[idx] += q, gvec[idx] = gvec[idx] + G;
        }
      }
    });
  }
  const double nfr = double(frames.size());
  // reflections of equal |G| merged: multiplicity (both signs), ⟨|F|²⟩ and the orientation factor per vector
  struct Ref { double q, f2sum, posum; int count; int h, k, l; std::vector<Vec3> dirs; };
  std::vector<size_t> order(nvec);
  for (size_t i = 0; i < nvec; ++i) order[i] = i, gsum[i] /= nfr, f2[i] /= nfr;
  std::sort(order.begin(), order.end(), [&](size_t x, size_t y) { return gsum[x] < gsum[y]; });
  std::vector<int> hh(nvec), kk(nvec), ll(nvec);
  for (size_t g = 0; g < hks.size(); ++g)
    for (size_t li = 0; li < hks[g].ls.size(); ++li) hh[offset[g] + li] = hks[g].h, kk[offset[g] + li] = hks[g].k, ll[offset[g] + li] = hks[g].ls[li];
  // the preferred direction: the reciprocal vector of its hkl (the normal of those planes)
  const bool po = (o.preferred[0] || o.preferred[1] || o.preferred[2]) && std::fabs(o.march_r - 1) > 1e-12;
  Vec3 pref{0, 0, 1};
  if (po) {
    pref = B0[0] * double(o.preferred[0]) + B0[1] * double(o.preferred[1]) + B0[2] * double(o.preferred[2]);
    pref = pref * (1 / norm(pref));
    if (o.march_r <= 0) throw std::invalid_argument("powder: the March–Dollase r must be positive");
  }
  auto march = [&](const Vec3& G) {
    if (!po) return 1.0;
    const double c = dot(G, pref) / std::max(1e-12, norm(G)), c2 = c * c, r = o.march_r;
    return std::pow(r * r * c2 + (1 - c2) / r, -1.5);
  };
  std::vector<Ref> refs;
  double f2max = 0;
  for (double v : f2) f2max = std::max(f2max, v);
  for (size_t ii = 0; ii < nvec; ++ii) {
    const size_t i = order[ii];
    if (gsum[i] > qmax) continue;
    const double tol = 1e-6 * gsum[i] + 1e-9;
    if (refs.empty() || gsum[i] - refs.back().q > tol) refs.push_back({gsum[i], 0, 0, 0, hh[i], kk[i], ll[i], {}});
    Ref& R = refs.back();
    const Vec3 G = gvec[i] * (1.0 / nfr);
    R.f2sum += 2 * f2[i];
    R.posum += 2 * f2[i] * march(G);
    R.count += 2;
    R.dirs.push_back(G);
    // the label: the largest of ±(h k l)
    auto key = [](int a, int b, int c) { return std::make_tuple(a, b, c); };
    if (key(hh[i], kk[i], ll[i]) > key(R.h, R.k, R.l)) R.h = hh[i], R.k = kk[i], R.l = ll[i];
    if (key(-hh[i], -kk[i], -ll[i]) > key(R.h, R.k, R.l)) R.h = -hh[i], R.k = -kk[i], R.l = -ll[i];
  }
  // absences: reflections whose structure factor vanishes
  refs.erase(std::remove_if(refs.begin(), refs.end(), [&](const Ref& R) { return R.f2sum / R.count <= 1e-10 * std::max(1e-300, f2max); }), refs.end());
  rep.vectors = int(2 * nvec), rep.frames = int(frames.size()), rep.reflections = int(refs.size());

  // widths at θ for line λ: (H, η) of the peak shape in degrees 2θ, and the crystallite size along G
  const bool aniso = o.size_abc[0] > 0 && o.size_abc[1] > 0 && o.size_abc[2] > 0;
  const Vec3 ua = s0.cell.a * (1 / norm(s0.cell.a)), ub = s0.cell.b * (1 / norm(s0.cell.b)), uc = s0.cell.c * (1 / norm(s0.cell.c));
  auto size_along = [&](const Ref& R) {
    if (aniso) {
      double tot = 0;   // an ellipsoid: the mean of 1/L² over the equivalent directions
      for (const auto& d : R.dirs) {
        const Vec3 nhat = d * (1 / std::max(1e-12, norm(d)));
        const double a1 = dot(nhat, ua) / o.size_abc[0], a2 = dot(nhat, ub) / o.size_abc[1], a3 = dot(nhat, uc) / o.size_abc[2];
        tot += a1 * a1 + a2 * a2 + a3 * a3;
      }
      return 1 / std::sqrt(tot / double(R.dirs.size()));
    }
    return o.size;
  };
  bool clamped = false;
  const double hmin = 2 * o.step;
  auto shape = [&](double th, double lam, double L) {
    const double tn = std::tan(th), cs = std::cos(th);
    const double strain = 4 * o.strain * tn * kDeg;
    const double hg2 = std::max(0.0, o.U * tn * tn + o.V * tn + o.W + strain * strain);
    const double hg = std::sqrt(hg2);
    const double hl = std::max(0.0, o.X * tn + o.Y / cs + (L > 0 ? o.scherrer_k * lam / (L * cs) * kDeg : 0.0));
    double H = 0, eta = 0;
    if (o.profile == "gaussian") H = std::sqrt(hg2 + hl * hl), eta = 0;
    else if (o.profile == "lorentzian") H = hg + hl, eta = 1;
    else {
      H = std::pow(std::pow(hg, 5) + 2.69269 * std::pow(hg, 4) * hl + 2.42843 * std::pow(hg, 3) * hl * hl + 4.47163 * hg * hg * std::pow(hl, 3) +
                       0.07842 * hg * std::pow(hl, 4) + std::pow(hl, 5), 0.2);
      if (o.profile == "pseudo-voigt") eta = std::clamp(o.eta, 0.0, 1.0);
      else if (o.profile == "tch") {
        const double r = H > 0 ? hl / H : 0;
        eta = std::clamp(1.36603 * r - 0.47719 * r * r + 0.11116 * r * r * r, 0.0, 1.0);
      } else throw std::invalid_argument("powder: profile tch, pseudo-voigt, gaussian or lorentzian");
    }
    if (H < hmin) H = hmin, clamped = true;
    return std::make_pair(H, eta);
  };
  auto lp = [&](double th) {
    const double s = std::sin(th), c = std::cos(th), c2t = std::cos(2 * th);
    const double lor = 1 / (s * s * c);
    if (neutron) return lor;
    if (o.polarised) return lor;
    const double cm = std::cos(o.monochromator_2theta / kDeg), cm2 = cm * cm;
    return lor * (1 + cm2 * c2t * c2t) / (1 + cm2);
  };

  // the pattern and the list
  const size_t np = size_t(std::floor((o.two_theta_max - o.two_theta_min) / o.step + 1e-9)) + 1;
  rep.two_theta.resize(np), rep.intensity.assign(np, 0.0), rep.q.resize(np);
  for (size_t i = 0; i < np; ++i) {
    rep.two_theta[i] = o.two_theta_min + o.step * double(i);
    rep.q[i] = 4 * kPi * std::sin(0.5 * rep.two_theta[i] / kDeg) / rep.wavelengths[0];
  }
  const double c4l2 = 4 * std::log(2.0), gn = 2 * std::sqrt(std::log(2.0) / kPi);
  std::vector<PowderPeak> list;
  for (const auto& R : refs) {
    const double L = size_along(R);
    for (size_t li = 0; li < rep.wavelengths.size(); ++li) {
      const double lam = rep.wavelengths[li], sn = lam * R.q / (4 * kPi);
      if (sn >= 1) continue;
      const double th = std::asin(sn), tt = 2 * th * kDeg;
      const double I = rep.weights[li] * R.posum * lp(th);
      const auto [H, eta] = shape(th, lam, L);
      if (li == 0 && tt >= o.two_theta_min && tt <= o.two_theta_max) {
        PowderPeak p;
        p.h = R.h, p.k = R.k, p.l = R.l, p.multiplicity = R.count, p.q = R.q, p.d = 2 * kPi / R.q, p.two_theta = tt;
        p.f = std::sqrt(R.f2sum / R.count), p.intensity = I, p.fwhm = H;
        list.push_back(p);
      }
      const double c = tt + o.zero_shift, half = 30 * H;
      if (c + half < o.two_theta_min || c - half > o.two_theta_max) continue;
      const size_t i0 = size_t(std::max(0.0, std::ceil((c - half - o.two_theta_min) / o.step)));
      const size_t i1 = std::min(np - 1, size_t(std::max(0.0, std::floor((c + half - o.two_theta_min) / o.step))));
      for (size_t i = i0; i <= i1 && i < np; ++i) {
        const double x = rep.two_theta[i] - c, u = x * x / (H * H);
        const double g = gn / H * std::exp(-c4l2 * u), lz = 2 / (kPi * H) / (1 + 4 * u);
        rep.intensity[i] += I * (eta * lz + (1 - eta) * g);
      }
    }
  }
  double imax = 0;
  for (const auto& p : list) imax = std::max(imax, p.intensity);
  for (auto& p : list) p.intensity = imax > 0 ? 100 * p.intensity / imax : 0;
  double ymax = *std::max_element(rep.intensity.begin(), rep.intensity.end());
  if (ymax > 0)
    for (auto& v : rep.intensity) v = 100 * v / ymax;
  std::stable_sort(list.begin(), list.end(), [](const PowderPeak& a, const PowderPeak& b) {
    return std::fabs(a.intensity - b.intensity) > 1e-9 ? a.intensity > b.intensity : a.two_theta < b.two_theta;
  });
  if (o.max_peaks >= 0 && list.size() > size_t(o.max_peaks)) list.resize(size_t(o.max_peaks));
  rep.peaks = std::move(list);

  // notes
  std::string lines = neutron ? "neutrons" : "X-ray";
  for (size_t li = 0; li < rep.wavelengths.size(); ++li) lines += (li ? " + " : " λ ") + fmt(rep.wavelengths[li], 6) + " Å" + (li ? " (× " + fmt(rep.weights[li], 3) + ")" : "");
  rep.notes.push_back(lines + " · " + std::to_string(rep.reflections) + " reflections (" + std::to_string(rep.vectors) + " vectors) from " + std::to_string(rep.frames) +
                      " frame" + (rep.frames == 1 ? "" : "s") + " · " + fmt(o.two_theta_min) + "–" + fmt(o.two_theta_max) + "° 2θ");
  rep.notes.push_back(std::string(neutron ? "coherent scattering lengths" : "Cromer–Mann form factors") + (biso > 0 ? ", B " + fmt(biso) + " Å²" : "") +
                      "; Lorentz" + (neutron || o.polarised ? "" : "–polarisation" + std::string(o.monochromator_2theta > 0 ? " with a monochromator at 2θ " + fmt(o.monochromator_2theta) + "°" : "")) +
                      "; " + o.profile + " peaks" + (o.size > 0 || aniso ? ", Scherrer size" : "") + (o.strain > 0 ? ", microstrain " + fmt(o.strain) : "") +
                      (po ? ", March–Dollase r " + fmt(o.march_r) : ""));
  if (clamped) rep.notes.push_back("some peaks were narrower than two steps: drawn at " + fmt(hmin) + "° so the pattern samples them");
  if (rep.frames > 1) rep.notes.push_back("|F|² averaged over the frames: their thermal motion damps the reflections and adds diffuse scattering");
  rep.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  return rep;
}

}  // namespace caps
