// CAPS thermal conductivity by MD (see thermal.hpp).
#include "caps/thermal.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <numeric>
#include <sstream>
#include <stdexcept>

#include "caps/elements.hpp"
#include "caps/uff.hpp"
#include "cell_list.hpp"

namespace caps {

namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr double kCoulombK = 332.06371;            // kcal Å / (mol e²), as LAMMPS real units
constexpr double kBoltz = 0.0019872043;            // kcal/mol/K
constexpr double kMv2 = 1.0 / 4.184e-4;            // g/mol · Å²/fs² in kcal/mol

inline double erfc_e(double x, double& ex2) {   // Abramowitz and Stegun 7.1.26, as CAPS's evaluator and LAMMPS coul/dsf
  const double t = 1.0 / (1.0 + 0.3275911 * x);
  ex2 = std::exp(-x * x);
  return t * (0.254829592 + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429)))) * ex2;
}

std::string fmtd(double v, int p = 3) {
  char b[64];
  std::snprintf(b, sizeof b, "%.*g", p, v);
  return b;
}

}  // namespace

HeatFlux::HeatFlux(const ForceField& ff, const EnergyOptions& e) : ff_(ff), e_(e) {
  auto no = [](const std::string& what) {
    throw FieldError("heat flux: " + what + " are not supported by CAPS's own heat flux; write the LAMMPS inputs (caps thermal --lammps-input), where compute centroid/stress/atom handles them");
  };
  if (!ff.bonds_x.empty()) no("bonds other than harmonic");
  if (!ff.angles_x.empty()) no("angles other than harmonic");
  if (!ff.bonds2.empty() || !ff.angles2.empty() || !ff.dihedrals2.empty() || !ff.impropers2.empty()) no("class II terms");
  if (!ff.inversions.empty()) no("inversion impropers");
  if (!ff.cbt.empty()) no("combined bending-torsion terms");
  if (!ff.lj_pairs.empty()) no("explicit LJ pairs");
  if (!ff.pair_func.empty()) no("pair forms other than 12-6 Lennard-Jones");
  if (ff.pair_form != "lj12-6") no("9-6 Lennard-Jones pairs");
  if (!ff.vsites.empty()) no("virtual sites");
  if (ff.sw.on || ff.manybody.on() || ff.hbond.on()) no("many-body and hydrogen-bond potentials");
  if (ff.coul_rf || ff.coul_gromacs || ff.lj_fsw || ff.keep13 || !ff.excluded_type_pairs.empty()) no("reaction-field, switched or MARTINI non-bonded settings");
  const size_t n = ff.type_index.size();
  mass_ = ff.mass;
  if (mass_.size() != n) throw FieldError("heat flux: the force field carries no masses");
  if (ff.cutoff > 0) e_.cutoff = ff.cutoff;
  rc_ = e_.cutoff;
  shift_ = !e_.tail;
  q_ = ff.charge;
  if (ff.dielectric != 1) for (auto& c : q_) c /= std::sqrt(ff.dielectric);
  coul_ = e_.coulomb && std::any_of(q_.begin(), q_.end(), [](double c) { return c != 0; });
  if (coul_ && e_.electrostatics == EnergyOptions::Electrostatics::PME)
    throw FieldError("heat flux: PME has no per-atom heat flux in CAPS; run with damped shifted force electrostatics (--dsf, the default)");
  nt_ = ff.lj.size();
  eps_.resize(nt_ * nt_), s6_.resize(nt_ * nt_);
  for (size_t a = 0; a < nt_; ++a)
    for (size_t b = 0; b < nt_; ++b) {
      const PairType pt = mixed_pair(ff, int(a), int(b));
      eps_[a * nt_ + b] = pt.eps, s6_[a * nt_ + b] = std::pow(pt.sigma, 6);
    }
  if (!ff.lj14_types.empty()) {
    ForceField f14;
    f14.mixing = ff.mixing;
    f14.lj = ff.lj14_types;
    eps14_.resize(nt_ * nt_), s614_.resize(nt_ * nt_);
    for (size_t a = 0; a < nt_; ++a)
      for (size_t b = 0; b < nt_; ++b) {
        const PairType pt = mixed_pair(f14, int(a), int(b));
        eps14_[a * nt_ + b] = pt.eps, s614_[a * nt_ + b] = std::pow(pt.sigma, 6);
      }
  }
  if (coul_) {   // DSF self energy of each atom (no force), in e_i as LAMMPS's pe/atom has it
    const double a = e_.dsf_alpha;
    double ex;
    const double erc = erfc_e(a * rc_, ex);
    const double f0 = erc / (rc_ * rc_) + 2 * a / std::sqrt(kPi) * std::exp(-a * a * rc_ * rc_) / rc_;
    self_.resize(q_.size());
    for (size_t i = 0; i < q_.size(); ++i) self_[i] = -(erc / (2 * rc_) + f0 * rc_ / 2 + a / std::sqrt(kPi)) * q_[i] * q_[i] * kCoulombK;
  }
}

void HeatFlux::per_atom(const std::vector<double>& x, const Cell& cell, std::vector<double>& e, std::vector<std::array<double, 9>>& w,
                        std::vector<double>* forces) {
  const size_t n = ff_.type_index.size();
  if (x.size() != 3 * n) throw FieldError("heat flux: positions do not match the force field");
  e.assign(n, 0.0);
  w.assign(n, std::array<double, 9>{});
  if (forces) forces->assign(3 * n, 0.0);
  if (!self_.empty()) e = self_;
  auto P = [&](uint32_t i) { return Vec3{x[3 * i], x[3 * i + 1], x[3 * i + 2]}; };
  auto mi = [&](const Vec3& d) { return cell.valid() ? cell.minimum_image(d) : d; };
  // a term over atoms ids with relative positions rel (about the first atom) and forces f: centroid virial, energy split
  auto term = [&](const uint32_t* ids, const Vec3* rel, const Vec3* f, int m, double en) {
    Vec3 c{0, 0, 0};
    for (int k = 0; k < m; ++k) c = c + rel[k];
    c = c * (1.0 / m);
    for (int k = 0; k < m; ++k) {
      const Vec3 d = rel[k] - c;
      auto& W = w[ids[k]];
      for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b) W[3 * a + b] += d[a] * f[k][b];
      e[ids[k]] += en / m;
      if (forces) for (int a = 0; a < 3; ++a) (*forces)[3 * ids[k] + a] += f[k][a];
    }
  };
  auto pair = [&](uint32_t i, uint32_t j, const Vec3& d /* r_j − r_i */, double en, double fr /* −(dE/dr)/r */) {
    const uint32_t ids[2] = {i, j};
    const Vec3 rel[2] = {{0, 0, 0}, d};
    const Vec3 fj = d * fr;
    const Vec3 f[2] = {fj * -1.0, fj};
    term(ids, rel, f, 2, en);
  };
  const double rc = rc_, rc2 = rc * rc;
  const double a = e_.dsf_alpha;
  double exrc;
  const double erc = erfc_e(a * rc, exrc);
  const double e0 = erc / rc, f0 = erc / rc2 + 2 * a / std::sqrt(kPi) * std::exp(-a * a * rc2) / rc;
  const double a2pi = 2 * a / std::sqrt(kPi);
  auto lj = [&](double eps, double s6, double r2, double scale, double& en) {
    const double q = s6 / (r2 * r2 * r2);
    en += scale * 4 * eps * (q * q - q);
    if (shift_) {
      const double qc = s6 / (rc2 * rc2 * rc2);
      en -= scale * 4 * eps * (qc * qc - qc);
    }
    return scale * 24 * eps * (2 * q * q - q) / r2;
  };
  // --- non-bonded pairs (not excluded) inside the cut-off
  if (cell.valid()) {
    const double v = cell.volume();
    const double wd[3] = {v / norm(cross(cell.b, cell.c)), v / norm(cross(cell.c, cell.a)), v / norm(cross(cell.a, cell.b))};
    if (std::min({wd[0], wd[1], wd[2]}) < 2 * rc) throw FieldError("heat flux: the cell must be at least twice the cut-off across (" + fmtd(2 * rc) + " Å)");
  }
  System tmp;
  tmp.cell = cell;
  tmp.atoms.resize(n);
  for (size_t i = 0; i < n; ++i) tmp.atoms[i].pos = P(uint32_t(i));
  const Grid g(tmp, rc);
  for (int bx = 0; bx < g.n[0]; ++bx)
    for (int by = 0; by < g.n[1]; ++by)
      for (int bz = 0; bz < g.n[2]; ++bz)
        for (uint32_t i : g.bins[g.index(bx, by, bz)])
          g.for_neighbour_bins(bx, by, bz, [&](const std::vector<uint32_t>& bin) {
            for (uint32_t j : bin) {
              if (j <= i) continue;
              const Vec3 d = g.sep(i, j);
              const double r2 = dot(d, d);
              if (r2 >= rc2) continue;
              const auto& ex = ff_.excluded.size() == n ? ff_.excluded[i] : std::vector<uint32_t>{};
              if (std::binary_search(ex.begin(), ex.end(), j)) continue;
              const size_t tp = size_t(ff_.type_index[i]) * nt_ + ff_.type_index[j];
              double en = 0, fr = 0;
              if (eps_[tp] > 0) fr = lj(eps_[tp], s6_[tp], r2, 1.0, en);
              if (coul_ && q_[i] != 0 && q_[j] != 0) {
                const double r = std::sqrt(r2), qq = kCoulombK * q_[i] * q_[j];
                double ex2;
                const double er = erfc_e(a * r, ex2);
                en += qq * (er / r - e0 + f0 * (r - rc));
                fr += qq * (er / r2 + a2pi * ex2 / r - f0) / r;
              }
              if (en != 0 || fr != 0) pair(i, j, d, en, fr);
            }
          });
  // --- 1-4 Lennard-Jones
  for (size_t k = 0; k < ff_.pairs14.size(); ++k) {
    const uint32_t i = ff_.pairs14[k][0], j = ff_.pairs14[k][1];
    const Vec3 d = mi(P(j) - P(i));
    const double r2 = dot(d, d);
    const size_t tp = size_t(ff_.type_index[i]) * nt_ + ff_.type_index[j];
    const double s = k < ff_.pairs14_lj.size() ? ff_.pairs14_lj[k] : ff_.lj14;
    double en = 0;
    const double fr = eps14_.empty() ? lj(eps_[tp], s6_[tp], r2, s, en) : lj(eps14_[tp], s614_[tp], r2, s, en);
    pair(i, j, d, en, fr);
  }
  // --- electrostatics of excluded partners (DSF: the damped sum includes them, (1 − factor) of 1/r removed)
  if (coul_) {
    std::map<std::pair<uint32_t, uint32_t>, double> p14;
    for (size_t k = 0; k < ff_.pairs14.size(); ++k) {
      const auto& p = ff_.pairs14[k];
      p14[{std::min(p[0], p[1]), std::max(p[0], p[1])}] = k < ff_.pairs14_coul.size() ? ff_.pairs14_coul[k] : ff_.coul14;
    }
    for (uint32_t i = 0; i < ff_.excluded.size(); ++i)
      for (uint32_t j : ff_.excluded[i]) {
        if (j <= i || q_[i] == 0 || q_[j] == 0) continue;
        const Vec3 d = mi(P(j) - P(i));
        const double r2 = dot(d, d);
        if (r2 >= rc2) continue;
        const auto it = p14.find({i, j});
        const double factor = it == p14.end() ? 0.0 : it->second;
        const double r = std::sqrt(r2), qq = kCoulombK * q_[i] * q_[j];
        double ex2;
        const double er = erfc_e(a * r, ex2);
        const double en = qq * (er / r - e0 + f0 * (r - rc)) - (1 - factor) * qq / r;
        const double fr = qq * ((er / r2 + a2pi * ex2 / r - f0) - (1 - factor) / r2) / r;
        pair(i, j, d, en, fr);
      }
  }
  // --- bonded terms
  for (const auto& bd : ff_.bonds) {
    const Vec3 d = mi(P(bd.j) - P(bd.i));
    const double r = norm(d), dr = r - bd.r0;
    pair(bd.i, bd.j, d, bd.k * dr * dr, r > 1e-12 ? -2 * bd.k * dr / r : 0.0);
  }
  for (const auto& ub : ff_.urey_bradley) {
    const Vec3 d = mi(P(ub.k) - P(ub.i));
    const double r = norm(d), dr = r - ub.r0;
    pair(ub.i, ub.k, d, ub.kub * dr * dr, r > 1e-12 ? -2 * ub.kub * dr / r : 0.0);
  }
  for (const auto& an : ff_.angles) {
    const Vec3 u = mi(P(an.i) - P(an.j)), v = mi(P(an.k) - P(an.j));
    const double lu = norm(u), lv = norm(v);
    const double c = std::clamp(dot(u, v) / (lu * lv), -1.0, 1.0);
    const double th = std::acos(c), sn = std::max(std::sqrt(1 - c * c), 1e-8);
    const double dth = th - an.theta0;
    const double gg = 2 * an.kt * dth / sn;
    const Vec3 fi = (v * (1 / (lu * lv)) - u * (c / (lu * lu))) * gg;
    const Vec3 fk = (u * (1 / (lu * lv)) - v * (c / (lv * lv))) * gg;
    const uint32_t ids[3] = {an.j, an.i, an.k};
    const Vec3 rel[3] = {{0, 0, 0}, u, v};
    const Vec3 f[3] = {(fi + fk) * -1.0, fi, fk};
    term(ids, rel, f, 3, an.kt * dth * dth);
  }
  // torsions: a run of terms on one quadruple shares its geometry
  auto torsions = [&](const std::vector<TorsionTerm>& T) {
    for (size_t t = 0; t < T.size();) {
      const TorsionTerm& t0 = T[t];
      size_t end = t + 1;
      while (end < T.size() && T[end].i == t0.i && T[end].j == t0.j && T[end].k == t0.k && T[end].l == t0.l) ++end;
      const Vec3 b1 = mi(P(t0.j) - P(t0.i)), b2 = mi(P(t0.k) - P(t0.j)), b3 = mi(P(t0.l) - P(t0.k));
      const Vec3 m = cross(b1, b2), nn = cross(b2, b3);
      const double m2 = dot(m, m), n2 = dot(nn, nn), lb2 = norm(b2);
      if (m2 < 1e-12 || n2 < 1e-12 || lb2 < 1e-9) { t = end; continue; }
      const double phi = std::atan2(lb2 * dot(b1, nn), dot(m, nn));
      double en = 0, dEdphi = 0;
      for (size_t u = t; u < end; ++u) {
        en += T[u].v * (1 + std::cos(T[u].n * phi - T[u].delta));
        dEdphi -= T[u].v * T[u].n * std::sin(T[u].n * phi - T[u].delta);
      }
      const Vec3 g1 = m * (-lb2 / m2), g4 = nn * (lb2 / n2);
      const double p = dot(b1, b2) / (lb2 * lb2), qd = dot(b3, b2) / (lb2 * lb2);
      const Vec3 g2 = g4 * qd - g1 * (1 + p), g3 = g1 * p - g4 * (1 + qd);
      const uint32_t ids[4] = {t0.i, t0.j, t0.k, t0.l};
      const Vec3 rel[4] = {{0, 0, 0}, b1, b1 + b2, b1 + b2 + b3};
      const Vec3 f[4] = {g1 * -dEdphi, g2 * -dEdphi, g3 * -dEdphi, g4 * -dEdphi};
      term(ids, rel, f, 4, en);
      t = end;
    }
  };
  torsions(ff_.dihedrals);
  torsions(ff_.impropers);
  for (const auto& h : ff_.impropers_harmonic) {
    const Vec3 b1 = mi(P(h.j) - P(h.i)), b2 = mi(P(h.k) - P(h.j)), b3 = mi(P(h.l) - P(h.k));
    const Vec3 m = cross(b1, b2), nn = cross(b2, b3);
    const double m2 = dot(m, m), n2 = dot(nn, nn), lb2 = norm(b2);
    if (m2 < 1e-12 || n2 < 1e-12 || lb2 < 1e-9) continue;
    const double chi = std::atan2(lb2 * dot(b1, nn), dot(m, nn));
    double dchi = chi - h.chi0;
    dchi -= 2 * kPi * std::round(dchi / (2 * kPi));
    const double dEdphi = 2 * h.k2 * dchi;
    const Vec3 g1 = m * (-lb2 / m2), g4 = nn * (lb2 / n2);
    const double p = dot(b1, b2) / (lb2 * lb2), qd = dot(b3, b2) / (lb2 * lb2);
    const Vec3 g2 = g4 * qd - g1 * (1 + p), g3 = g1 * p - g4 * (1 + qd);
    const uint32_t ids[4] = {h.i, h.j, h.k, h.l};
    const Vec3 rel[4] = {{0, 0, 0}, b1, b1 + b2, b1 + b2 + b3};
    const Vec3 f[4] = {g1 * -dEdphi, g2 * -dEdphi, g3 * -dEdphi, g4 * -dEdphi};
    term(ids, rel, f, 4, h.k2 * dchi * dchi);
  }
}

Vec3 HeatFlux::flux(const std::vector<double>& x, const std::vector<double>& v0, const Cell& cell) {
  std::vector<double> e;
  std::vector<std::array<double, 9>> w;
  per_atom(x, cell, e, w);
  // velocities about the centre of mass: a drift of the whole cell carries energy but no heat
  Vec3 vc{0, 0, 0};
  double mt = 0;
  for (size_t i = 0; i < e.size(); ++i) {
    for (int k = 0; k < 3; ++k) vc[k] += mass_[i] * v0[3 * i + size_t(k)];
    mt += mass_[i];
  }
  vc = vc * (1.0 / mt);
  std::vector<double> v(v0);
  for (size_t i = 0; i < e.size(); ++i)
    for (int k = 0; k < 3; ++k) v[3 * i + size_t(k)] -= vc[k];
  Vec3 q{0, 0, 0};
  for (size_t i = 0; i < e.size(); ++i) {
    const double vx = v[3 * i], vy = v[3 * i + 1], vz = v[3 * i + 2];
    const double ei = e[i] + 0.5 * mass_[i] * (vx * vx + vy * vy + vz * vz) * kMv2;
    const auto& W = w[i];
    q[0] += ei * vx + W[0] * vx + W[1] * vy + W[2] * vz;
    q[1] += ei * vy + W[3] * vx + W[4] * vy + W[5] * vz;
    q[2] += ei * vz + W[6] * vx + W[7] * vy + W[8] * vz;
  }
  return q;
}

// ---------------------------------------------------------------------------------------------------------------

ThermalResult green_kubo_kappa(const std::vector<Vec3>& q0, double dt_fs, double volume, double temperature, double corr_ps, int blocks) {
  ThermalResult r;
  r.method = "gk";
  const size_t n = q0.size();
  if (n < 10) throw std::invalid_argument("Green–Kubo needs more heat-flux samples (a longer run)");
  size_t L = std::min(n / 2, size_t(std::max(2.0, std::round(corr_ps * 1000 / dt_fs))));
  // no mean taken off: the flux's mean is zero at rest (velocities about the centre of mass), and taking off the sample mean
  // of a run only a few correlation times long pulls the integral down (the autocovariance of a mean-free series sums to
  // about zero over the run)
  const std::vector<Vec3>& q = q0;
  const double pre = dt_fs / (volume * kBoltz * temperature * temperature) * kKappaSI;
  auto acf_of = [&](size_t a, size_t b, int k, std::vector<double>& c) {   // samples [a, b), axis k, lags 0..L
    c.assign(L + 1, 0.0);
    for (size_t lag = 0; lag <= L; ++lag) {
      double s = 0;
      size_t cnt = 0;
      for (size_t t = a; t + lag < b; ++t) s += q[t][k] * q[t + lag][k], ++cnt;
      c[lag] = cnt ? s / double(cnt) : 0;
    }
  };
  auto integral = [&](const std::vector<double>& c, std::vector<double>* run) {
    double s = 0, tail = 0;
    size_t nt = 0;
    if (run) run->assign(c.size(), 0.0);
    for (size_t lag = 1; lag < c.size(); ++lag) {
      s += 0.5 * (c[lag - 1] + c[lag]) * pre;
      if (run) (*run)[lag] = s;
      if (lag >= 2 * (c.size() - 1) / 3) tail += s, ++nt;
    }
    return nt ? tail / double(nt) : s;
  };
  r.t_ps.resize(L + 1);
  for (size_t lag = 0; lag <= L; ++lag) r.t_ps[lag] = lag * dt_fs * 1e-3;
  const int nb = std::clamp(blocks, 2, int(n / (2 * L + 2) + 1));
  for (int k = 0; k < 3; ++k) {
    std::vector<double> c;
    acf_of(0, n, k, c);
    r.kappa[size_t(k)] = integral(c, &r.running[size_t(k)]);
    r.acf[size_t(k)].resize(c.size());
    for (size_t lag = 0; lag < c.size(); ++lag) r.acf[size_t(k)][lag] = c[0] != 0 ? c[lag] / c[0] : 0;
    // block errors
    std::vector<double> kb;
    for (int bkk = 0; bkk < nb; ++bkk) {
      const size_t a = n * size_t(bkk) / size_t(nb), b = n * size_t(bkk + 1) / size_t(nb);
      if (b - a < L + 2) continue;
      std::vector<double> cb;
      acf_of(a, b, k, cb);
      kb.push_back(integral(cb, nullptr));
    }
    if (kb.size() >= 2) {
      const double m = std::accumulate(kb.begin(), kb.end(), 0.0) / double(kb.size());
      double s2 = 0;
      for (double x : kb) s2 += (x - m) * (x - m);
      r.error[size_t(k)] = std::sqrt(s2 / double(kb.size() - 1) / double(kb.size()));
    }
    // plateau: the running integral's last third moves less than its block error (or 20 %)
    const auto& run = r.running[size_t(k)];
    const size_t a3 = 2 * L / 3;
    double lo = 1e300, hi = -1e300;
    for (size_t lag = a3; lag <= L; ++lag) lo = std::min(lo, run[lag]), hi = std::max(hi, run[lag]);
    if (hi - lo > std::max(0.2 * std::fabs(r.kappa[size_t(k)]), 2 * r.error[size_t(k)])) r.plateau = false;
  }
  r.kappa_mean = (r.kappa[0] + r.kappa[1] + r.kappa[2]) / 3;
  r.kappa_mean_error = std::sqrt(r.error[0] * r.error[0] + r.error[1] * r.error[1] + r.error[2] * r.error[2]) / 3;
  r.notes.push_back(std::to_string(n) + " heat-flux samples " + fmtd(dt_fs) + " fs apart; the integral to " + fmtd(L * dt_fs * 1e-3) +
                    " ps, κ the running integral's mean over its last third; errors from " + std::to_string(nb) + " blocks");
  if (!r.plateau) r.notes.push_back("the running integral has not levelled off over its last third: run longer (the noise of a Green–Kubo integral falls as 1/√(run length))");
  return r;
}

ThermalResult nemd_kappa(const std::vector<double>& T, double width, double exchanged, double time_fs, double area) {
  ThermalResult r;
  r.method = "nemd";
  const int n = int(T.size());
  if (n < 6 || n % 2) throw std::invalid_argument("NEMD needs an even number of slabs (at least 6)");
  r.slab_t = T;
  for (int k = 0; k < n; ++k) r.slab_x.push_back((k + 0.5) * width);
  // two straight lines: slabs 2 … n/2 − 2 (cold → hot) and n/2 + 2 … n − 2 (hot → cold), the slopes' magnitudes averaged
  auto slope = [&](int a, int b) {
    double sx = 0, sy = 0, sxx = 0, sxy = 0;
    int m = 0;
    for (int k = a; k <= b; ++k) {
      const double x = k * width, y = T[size_t(k)];
      sx += x, sy += y, sxx += x * x, sxy += x * y, ++m;
    }
    const double d = m * sxx - sx * sx;
    return d != 0 ? (m * sxy - sx * sy) / d : 0.0;
  };
  const int h = n / 2;
  const double g1 = slope(2, h - 2), g2 = slope(h + 2, n - 2);
  const double g = 0.5 * (std::fabs(g1) + std::fabs(g2));   // K/Å
  const double J = exchanged / (2 * area * time_fs);        // kcal/mol/(Å² fs)
  r.exchanged = exchanged;
  r.gradient = g * 1e10;
  r.flux = J * 4184.0 / 6.02214076e23 / (1e-20 * 1e-15);
  const double kappa = g > 0 ? J / g * kKappaSI : 0;
  r.kappa = {kappa, kappa, kappa};
  r.kappa_mean = kappa;
  // an error from the two halves' slopes
  if (g1 != 0 && g2 != 0) {
    const double k1 = J / std::fabs(g1) * kKappaSI, k2 = J / std::fabs(g2) * kKappaSI;
    r.kappa_mean_error = 0.5 * std::fabs(k1 - k2);
    r.error = {r.kappa_mean_error, r.kappa_mean_error, r.kappa_mean_error};
  }
  if (std::fabs(g1) > 0 && std::fabs(g2) > 0 && std::fabs(std::fabs(g1) - std::fabs(g2)) > 0.25 * g)
    r.notes.push_back("the two halves' gradients differ by more than 25 %: the profile is not yet steady or the run is noisy");
  return r;
}

double sheet_thickness(const System& s) {
  if (s.atoms.empty()) return 0;
  // z unwrapped about the first atom
  const double Lz = s.cell.valid() ? s.cell.c[2] : 0;
  size_t lo = 0, hi = 0;
  std::vector<double> z(s.atoms.size());
  for (size_t i = 0; i < s.atoms.size(); ++i) {
    double d = s.atoms[i].pos[2] - s.atoms[0].pos[2];
    if (Lz > 0) d -= Lz * std::round(d / Lz);
    z[i] = d;
    if (z[i] < z[lo]) lo = i;
    if (z[i] > z[hi]) hi = i;
  }
  return z[hi] - z[lo] + element(s.atoms[lo].element).vdw + element(s.atoms[hi].element).vdw;
}

ThermalResult thermal_conductivity(System& s, const ThermalOptions& o) {
  if (!s.cell.valid()) throw std::invalid_argument("thermal conductivity needs a periodic cell");
  if (o.method != "gk" && o.method != "nemd") throw std::invalid_argument("method: gk or nemd");
  if (o.dimensionality != 2 && o.dimensionality != 3) throw std::invalid_argument("dimensionality: 2 (a sheet, vacuum along z) or 3 (bulk)");
  const bool two = o.dimensionality == 2;
  double thick = 0;
  if (two) {
    if (std::fabs(s.cell.c[0]) > 1e-6 || std::fabs(s.cell.c[1]) > 1e-6 || std::fabs(s.cell.a[2]) > 1e-6 || std::fabs(s.cell.b[2]) > 1e-6)
      throw std::invalid_argument("2D: the sheet lies in the xy plane with the cell's c along z");
    thick = o.thickness > 0 ? o.thickness : sheet_thickness(s);
    const double extent = sheet_thickness(s);
    if (s.cell.c[2] - extent < 2 * o.energy.cutoff)
      throw std::invalid_argument("2D: the vacuum along z (" + fmtd(s.cell.c[2] - extent) + " Å) must exceed twice the cut-off, so the sheet does not see its image");
    if (o.method == "nemd" && o.axis == 2) throw std::invalid_argument("2D: the NEMD flux runs in the plane (axis x or y)");
  }
  const auto field = o.field ? o.field : std::make_shared<const ForceField>(default_forcefield(s));
  HeatFlux hf(*field, o.energy);   // checks the terms (both methods: NEMD needs none, but the same limits keep them alike)
  const auto& m = hf.masses();
  DynamicsOptions d;
  d.field = field;
  d.energy = o.energy;
  d.dt = o.dt;
  d.temperature = o.temperature;
  d.tau_t = o.tau_t;
  d.seed = o.seed;
  d.frame_every = 0;
  d.new_velocities = o.new_velocities || s.velocities.size() != s.atoms.size();
  const int64_t neq = std::llround(o.equilibrate_ps * 1000 / o.dt), nrun = std::llround(o.ps * 1000 / o.dt);
  auto stop = [&](const std::string& w, double f) { return o.progress && !o.progress(w, f); };
  if (neq > 0) {
    d.thermostat = Thermostat::Bussi;
    d.steps = neq;
    d.thermo_every = 100;
    d.progress = [&](const ThermoRow& r) { return !stop("equilibrating at " + fmtd(o.temperature, 4) + " K (NVT)", double(r.step) / double(neq + nrun)); };
    run_dynamics(s, d);
    d.new_velocities = false;
  }
  // the sampled run: NVE
  d.thermostat = Thermostat::None;
  d.steps = nrun;
  d.step_offset = neq;
  d.thermo_every = 100;
  double tsum = 0;
  int tn = 0;
  d.progress = [&](const ThermoRow& r) {
    if (r.step > neq) tsum += r.temperature, ++tn;
    return !stop(o.method == "gk" ? "sampling the heat flux (NVE)" : "exchanging kinetic energy (NVE)", double(r.step) / double(neq + nrun));
  };
  ThermalResult R;
  if (o.method == "gk") {
    std::vector<Vec3> q;
    const int every = std::max(1, o.sample_every);
    d.after_step = [&](const std::vector<double>& x, std::vector<double>& v, const Cell& cell, int64_t step) {
      if ((step - neq) % every == 0) q.push_back(hf.flux(x, v, cell));
    };
    run_dynamics(s, d);
    const double vol = two ? norm(cross(s.cell.a, s.cell.b)) * thick : s.cell.volume();
    R = green_kubo_kappa(q, o.dt * every, vol, o.temperature, o.corr_ps, o.blocks);
    if (two) {   // in-plane components only
      R.kappa[2] = R.error[2] = 0;
      R.running[2].clear(), R.acf[2].clear();
      R.kappa_mean = 0.5 * (R.kappa[0] + R.kappa[1]);
      R.kappa_mean_error = 0.5 * std::sqrt(R.error[0] * R.error[0] + R.error[1] * R.error[1]);
    }
    R.notes.insert(R.notes.begin(), fmtd(o.equilibrate_ps, 4) + " ps NVT (Bussi) then " + fmtd(o.ps, 4) + " ps NVE at " + fmtd(o.dt) + " fs steps");
  } else {
    const int ax = std::clamp(o.axis, 0, 2);
    const int ns = std::max(6, o.slabs + (o.slabs % 2));
    const Vec3 edge = ax == 0 ? s.cell.a : ax == 1 ? s.cell.b : s.cell.c;
    const double L = norm(edge);
    if (std::fabs(edge[ax] - L) > 1e-6 * L) throw std::invalid_argument("NEMD: the flux axis must be a cell edge along x, y or z");
    const double area = two ? s.cell.volume() / L / s.cell.c[2] * thick : s.cell.volume() / L;
    const size_t n = s.atoms.size();
    std::vector<double> tsl(size_t(ns), 0.0);
    std::vector<double> cnt(size_t(ns), 0.0);
    const int64_t skip = std::llround(o.transient * double(nrun));
    double moved = 0, moved_all = 0;
    int swaps = 0;
    const int every = std::max(1, o.sample_every);
    auto slab = [&](const std::vector<double>& x, size_t i) {
      double f = (x[3 * i + size_t(ax)] - s.cell.origin[size_t(ax)]) / L;
      f -= std::floor(f);
      return std::min(ns - 1, int(f * ns));
    };
    d.after_step = [&](const std::vector<double>& x, std::vector<double>& v, const Cell&, int64_t step) {
      const int64_t k = step - neq;
      if (k % std::max(1, o.swap_every) == 0) {
        // the hottest atom of the cold slab gives its velocity to the coldest atom of that mass in the hot slab
        size_t hot = SIZE_MAX, cold = SIZE_MAX;
        double kh = -1, kc = 1e300;
        auto ke = [&](size_t i) { return 0.5 * m[i] * (v[3 * i] * v[3 * i] + v[3 * i + 1] * v[3 * i + 1] + v[3 * i + 2] * v[3 * i + 2]) * kMv2; };
        for (size_t i = 0; i < n; ++i)
          if (slab(x, i) == 0 && ke(i) > kh) kh = ke(i), hot = i;
        if (hot != SIZE_MAX)
          for (size_t i = 0; i < n; ++i)
            if (slab(x, i) == ns / 2 && std::fabs(m[i] - m[hot]) < 1e-6 && ke(i) < kc) kc = ke(i), cold = i;
        if (hot != SIZE_MAX && cold != SIZE_MAX && kh > kc) {
          for (int c = 0; c < 3; ++c) std::swap(v[3 * hot + size_t(c)], v[3 * cold + size_t(c)]);
          moved_all += kh - kc;
          if (k > skip) moved += kh - kc, ++swaps;
        }
      }
      if (k > skip && k % every == 0) {
        std::vector<double> ek(size_t(ns), 0.0), nk(size_t(ns), 0.0);
        for (size_t i = 0; i < n; ++i) {
          const int b = slab(x, i);
          ek[size_t(b)] += m[i] * (v[3 * i] * v[3 * i] + v[3 * i + 1] * v[3 * i + 1] + v[3 * i + 2] * v[3 * i + 2]) * kMv2;
          nk[size_t(b)] += 1;
        }
        for (int b = 0; b < ns; ++b)
          if (nk[size_t(b)] > 0) tsl[size_t(b)] += ek[size_t(b)] / (3 * nk[size_t(b)] * kBoltz), cnt[size_t(b)] += 1;
      }
    };
    run_dynamics(s, d);
    for (int b = 0; b < ns; ++b) tsl[size_t(b)] = cnt[size_t(b)] > 0 ? tsl[size_t(b)] / cnt[size_t(b)] : 0;
    const double t_avg = double(nrun - skip) * o.dt;
    R = nemd_kappa(tsl, L / ns, moved, t_avg, area);
    R.axis = ax;
    R.exchanges = swaps;
    R.notes.insert(R.notes.begin(), fmtd(o.equilibrate_ps, 4) + " ps NVT (Bussi) then " + fmtd(o.ps, 4) + " ps NVE with an exchange every " +
                                        std::to_string(o.swap_every) + " steps along " + std::string(1, char('x' + ax)) + " (" + std::to_string(ns) +
                                        " slabs); profile and flux over the last " + fmtd(100 * (1 - o.transient), 3) + " %");
  }
  R.temperature = tn ? tsum / tn : 0;
  R.dimensionality = o.dimensionality;
  if (two) {
    R.thickness = thick;
    R.sheet_conductance = R.kappa_mean * thick * 1e-10;
    R.notes.push_back("2D: in-plane heat flux only; κ normalised by an effective thickness of " + fmtd(thick, 4) + " Å" +
                      (o.thickness > 0 ? " (given)" : " (the atoms' z extent plus their van der Waals radii; give the bulk parent's interlayer spacing with --thickness)") +
                      "; the sheet conductance κ·t = " + fmtd(R.sheet_conductance, 4) + " W/K does not depend on that choice");
  }
  return R;
}

std::vector<Property> thermal_properties(const ThermalResult& r) {
  Property q;
  q.id = "thermal_conductivity";
  q.unit = "W/(m·K)";
  q.notes = r.notes;
  if (r.method == "gk") {
    q.name = "Thermal conductivity (Green–Kubo)";
    q.value = r.kappa_mean;
    q.error = r.kappa_mean_error;
    q.method = "Green–Kubo integral of the heat-flux autocorrelation (centroid form for bonded terms), mean of x, y, z";
    for (int k = 0; k < 3; ++k) {
      q.series.push_back(Series{std::string("running integral ") + char('x' + k), "t (ps)", "κ (W/m/K)", r.t_ps, r.running[size_t(k)]});
      q.series.push_back(Series{std::string("autocorrelation ") + char('x' + k), "t (ps)", "⟨QQ⟩(t)/⟨QQ⟩(0)", r.t_ps, r.acf[size_t(k)]});
    }
  } else {
    q.name = "Thermal conductivity (NEMD, Müller-Plathe)";
    q.value = r.kappa_mean;
    q.error = r.kappa_mean_error;
    q.method = "reverse NEMD: kinetic energy exchanged between a cold and a hot slab, κ = J / |dT/dz|";
    q.series.push_back(Series{"temperature profile", "position (Å)", "T (K)", r.slab_x, r.slab_t});
  }
  return {q};
}

std::string thermal_summary(const ThermalResult& r) {
  char b[512];
  if (r.method == "gk" && r.dimensionality == 2) {
    std::snprintf(b, sizeof b, "κ (Green–Kubo, 2D, t = %.3g Å) x %.4g ± %.2g · y %.4g ± %.2g · in-plane mean %.4g ± %.2g W/(m K) · sheet conductance %.4g W/K at %.1f K%s\n",
                  r.thickness, r.kappa[0], r.error[0], r.kappa[1], r.error[1], r.kappa_mean, r.kappa_mean_error, r.sheet_conductance, r.temperature,
                  r.plateau ? "" : " (no plateau yet)");
  } else if (r.method == "gk") {
    std::snprintf(b, sizeof b, "κ (Green–Kubo) x %.4g ± %.2g · y %.4g ± %.2g · z %.4g ± %.2g · mean %.4g ± %.2g W/(m K) at %.1f K%s\n", r.kappa[0], r.error[0],
                  r.kappa[1], r.error[1], r.kappa[2], r.error[2], r.kappa_mean, r.kappa_mean_error, r.temperature, r.plateau ? "" : " (no plateau yet)");
  } else {
    if (r.dimensionality == 2) {
      std::snprintf(b, sizeof b, "κ (NEMD, 2D, t = %.3g Å, along %c) %.4g ± %.2g W/(m K) · sheet conductance %.4g W/K · %d exchanges · %.1f K\n", r.thickness,
                    char('x' + std::max(0, r.axis)), r.kappa_mean, r.kappa_mean_error, r.sheet_conductance, r.exchanges, r.temperature);
      return b;
    }
    std::snprintf(b, sizeof b, "κ (NEMD, along %c) %.4g ± %.2g W/(m K) · flux %.4g W/m² · gradient %.4g K/m · %d exchanges · %.1f K\n", char('x' + std::max(0, r.axis)),
                  r.kappa_mean, r.kappa_mean_error, r.flux, r.gradient, r.exchanges, r.temperature);
  }
  return b;
}

std::string thermal_lammps_run(const ThermalLammps& q) {
  std::ostringstream o;
  const long neq = std::lround(q.equilibrate_ps * 1000 / q.dt), nrun = std::lround(q.ps * 1000 / q.dt);
  const char ax = char('x' + std::clamp(q.axis, 0, 2));
  o << "\n# ---- CAPS thermal conductivity (" << (q.method == "gk" ? "Green-Kubo" : "Muller-Plathe NEMD") << "), units real\n"
    << "variable        T equal " << q.temperature << "\n"
    << "variable        dt equal " << q.dt << "\n"
    << "timestep        ${dt}\n"
    << "velocity        all create ${T} 4928459 dist gaussian mom yes rot yes\n"
    << "thermo          1000\nthermo_style    custom step temp press pe ke etotal\n"
    << "fix             equil all nvt temp ${T} ${T} 100.0\n"
    << "run             " << neq << "\n"
    << "unfix           equil\n"
    << "reset_timestep  0\n"
    << "fix             nve all nve\n";
  if (q.dimensionality == 2)
    o << "# 2D sheet (vacuum along z): no barostat here; an NPT stage before this would couple x and y only (x P P D y P P D)\n"
      << (q.kspace ? "kspace_modify   slab 3.0   # with boundary p p f set before read_data\n" : "");
  if (q.kspace)
    o << "# note: compute centroid/stress/atom does not take kspace terms (LAMMPS manual, Restrictions): with PPPM the flux below\n"
         "# misses the long-range part; for kappa write the inputs with coul/dsf (CAPS's default electrostatics)\n";
  if (q.method == "gk") {
    const long s = std::max(1, q.sample_every), p = std::max(2L, std::lround(q.corr_ps * 1000 / q.dt / s));
    o << "# heat flux: kinetic + potential energy per atom, centroid virial (correct for angle, dihedral and improper terms)\n"
      << "compute         ke all ke/atom\ncompute         pe all pe/atom\n"
      << "compute         stress all centroid/stress/atom NULL virial\n"
      << "compute         flux all heat/flux ke pe stress\n"
      << "variable        Jx equal c_flux[1]/vol\nvariable        Jy equal c_flux[2]/vol\nvariable        Jz equal c_flux[3]/vol\n"
      << "variable        s equal " << s << "\nvariable        p equal " << p << "\nvariable        d equal $p*$s\n"
      << "fix             JJ all ave/correlate $s $p $d c_flux[1] c_flux[2] c_flux[3] type auto file J0Jt.dat ave running\n"
      << "# kappa = V/(kB T^2) * integral, with the flux per volume: (kcal/mol)/(A fs K) -> W/(m K)\n"
      << "variable        kB equal 0.0019872043\n"
      << "variable        conv equal 4184.0/6.02214076e23/1.0e-25\n"
      << (q.dimensionality == 2 ? "# 2D sheet: the volume is the area times the effective thickness " + fmtd(q.thickness, 4) + " A (in-plane kappa only)\n"
                                  "variable        V equal lx*ly*" + fmtd(q.thickness, 6) + "\n"
                                : std::string("variable        V equal vol\n"))
      << "variable        scale equal ${conv}*${s}*${dt}/(${kB}*${T}*${T}*${V})\n"
      << "variable        kxx equal trap(f_JJ[3])*${scale}\nvariable        kyy equal trap(f_JJ[4])*${scale}\nvariable        kzz equal trap(f_JJ[5])*${scale}\n"
      << "# f_JJ exists only every $d steps: thermo then too, and the run a whole number of windows\n"
      << "thermo_style    custom step temp pe etotal v_kxx v_kyy v_kzz\nthermo          ${d}\n"
      << "run             " << (nrun + p * s - 1) / (p * s) * (p * s) << "\n"
      << "variable        k equal (v_kxx+v_kyy+v_kzz)/3.0\n"
      << "print           \"kappa (W/m/K): x ${kxx} y ${kyy} z ${kzz} mean ${k}\" file kappa.txt\n";
  } else {
    o << "# reverse NEMD (Muller-Plathe): kinetic energy swapped every " << q.swap_every << " steps between slab 1 (cold) and slab " << (q.slabs / 2 + 1)
      << " (hot)\n"
      << "compute         ke all ke/atom\nvariable        temp atom c_ke/(1.5*0.0019872043)\n"
      << "compute         layers all chunk/atom bin/1d " << ax << " lower " << 1.0 / q.slabs << " units reduced\n"
      << "fix             profile all ave/chunk 10 100 1000 layers v_temp file profile.mp\n"
      << "fix             mp all thermal/conductivity " << q.swap_every << " " << ax << " " << q.slabs << "\n"
      << "thermo_style    custom step temp pe etotal f_mp\n"
      << "run             " << nrun << "\n"
      << "# kappa = (f_mp / (2 A t)) / |dT/d" << ax << "|: f_mp the energy moved (kcal/mol), A the cross-section, t the time (fs),\n"
      << (q.dimensionality == 2 ? "# 2D sheet: A = (the in-plane edge across the flux) x " + fmtd(q.thickness, 4) + " A, the effective thickness\n" : std::string())
      << "# the gradient from profile.mp (straight lines through both halves, leaving out the slabs next to the hot and cold ones)\n"
      << "print           \"exchanged energy (kcal/mol) $(f_mp) over " << nrun * q.dt << " fs\" file exchanged.txt\n";
  }
  o << "write_data      thermal_final.data\n";
  return o.str();
}

}  // namespace caps
