// CAPS sorption: Widom insertion and grand-canonical Monte Carlo in a fixed host (see sorption.hpp).
#include "caps/sorption.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <limits>
#include <random>
#include <stdexcept>

#include "caps/pairmodel.hpp"

namespace caps {

namespace {

constexpr double kB = 0.0019872043;          // kcal/(mol·K)
constexpr double kBoltzmannJ = 1.380649e-23; // J/K
constexpr double kAvogadro = 6.02214076e23;
constexpr double kPi = 3.14159265358979323846;

struct Quat { double w, x, y, z; };
Vec3 rotate(const Quat& q, const Vec3& v) {   // q v q*
  const Vec3 u{q.x, q.y, q.z};
  const Vec3 t = cross(u, v) * 2.0;
  return v + t * q.w + cross(u, t);
}

// The fixed host atoms binned (bins at least rc/2 wide, two either side), and the pair model
struct Host {
  const System& s;
  const PairModel& pm;
  bool ortho = false;
  double L[3] = {0, 0, 0};
  int nb[3] = {1, 1, 1};
  std::vector<std::vector<uint32_t>> bins;
  std::vector<std::vector<int>> stencil;
  Host(const System& sys, const PairModel& p, size_t n_host) : s(sys), pm(p) {
    const Cell& c = sys.cell;
    ortho = std::fabs(c.a[1]) + std::fabs(c.a[2]) + std::fabs(c.b[0]) + std::fabs(c.b[2]) + std::fabs(c.c[0]) + std::fabs(c.c[1]) < 1e-9;
    L[0] = c.a[0], L[1] = c.b[1], L[2] = c.c[2];
    const double V = c.volume();
    const double w[3] = {V / norm(cross(c.b, c.c)), V / norm(cross(c.c, c.a)), V / norm(cross(c.a, c.b))};
    for (int d = 0; d < 3; ++d) nb[d] = std::max(1, int(std::floor(2 * w[d] / p.cutoff())));
    bins.resize(size_t(nb[0] * nb[1] * nb[2]));
    for (uint32_t i = 0; i < n_host; ++i) bins[size_t(bin_of(sys.atoms[i].pos))].push_back(i);
    stencil.resize(bins.size());
    for (int a = 0; a < nb[0]; ++a)
      for (int b = 0; b < nb[1]; ++b)
        for (int cc = 0; cc < nb[2]; ++cc) {
          auto offs = [&](int d, int v) {
            std::vector<int> r;
            if (nb[d] < 5) for (int k = 0; k < nb[d]; ++k) r.push_back(k);
            else for (int dv = -2; dv <= 2; ++dv) r.push_back(((v + dv) % nb[d] + nb[d]) % nb[d]);
            return r;
          };
          for (int i : offs(0, a))
            for (int j : offs(1, b))
              for (int k : offs(2, cc)) stencil[size_t((a * nb[1] + b) * nb[2] + cc)].push_back((i * nb[1] + j) * nb[2] + k);
        }
  }
  Vec3 mi(Vec3 d) const {
    if (!ortho) return s.cell.minimum_image(d);
    for (int k = 0; k < 3; ++k) d[size_t(k)] -= L[k] * std::round(d[size_t(k)] / L[k]);
    return d;
  }
  int bin_of(const Vec3& p) const {
    const Vec3 f = s.cell.to_fractional(p);
    int id[3];
    for (int d = 0; d < 3; ++d) id[d] = std::min(nb[d] - 1, int((f[size_t(d)] - std::floor(f[size_t(d)])) * nb[d]));
    return (id[0] * nb[1] + id[1]) * nb[2] + id[2];
  }
  // the template's atoms (indices tpl) at positions P with the host
  double energy(const std::vector<uint32_t>& tpl, const std::vector<Vec3>& P) const {
    double e = 0;
    for (size_t k = 0; k < tpl.size(); ++k) {
      for (int bn : stencil[size_t(bin_of(P[k]))])
        for (uint32_t b : bins[size_t(bn)]) {
          const Vec3 d = mi(s.atoms[b].pos - P[k]);
          e += pm.energy(tpl[k], b, dot(d, d));
        }
      if (!std::isfinite(e)) return e;
    }
    return e;
  }
  double between(const std::vector<uint32_t>& tpl, const std::vector<Vec3>& P, const std::vector<Vec3>& Q) const {
    double e = 0;
    for (size_t k = 0; k < tpl.size(); ++k)
      for (size_t q = 0; q < tpl.size(); ++q) {
        const Vec3 d = mi(Q[q] - P[k]);
        e += pm.energy(tpl[k], tpl[q], dot(d, d));
      }
    return e;
  }
};

}  // namespace

SorptionReport sorption(const System& s, const ForceField& ff, const SorptionOptions& o) {
  const auto t0 = std::chrono::steady_clock::now();
  SorptionReport rep;
  const size_t n = s.atoms.size();
  if (!s.cell.valid()) throw std::invalid_argument("sorption: the host needs a periodic cell");
  if (o.template_first_atom < 0 || size_t(o.template_first_atom) >= n) throw std::invalid_argument("sorption: no sorbate template (its atoms come last)");
  if (ff.type_index.size() != n || ff.charge.size() != n) throw std::invalid_argument("sorption: the force field was assigned to another structure");
  if (!ff.pair_func.empty() || ff.hbond.on() || ff.sw.on) throw std::invalid_argument("sorption: the energies need Lennard-Jones pairs");
  if (o.temperature <= 0) throw std::invalid_argument("sorption: T > 0");
  const size_t nh = size_t(o.template_first_atom);
  double rc = o.cutoff;
  {
    const double V = s.cell.volume();
    const double wmin = std::min({V / norm(cross(s.cell.b, s.cell.c)), V / norm(cross(s.cell.c, s.cell.a)), V / norm(cross(s.cell.a, s.cell.b))});
    if (rc > 0.5 * wmin - 1e-6) {
      rc = 0.5 * wmin - 1e-6;
      char b[160];
      std::snprintf(b, sizeof b, "cut-off reduced to %.2f Å (half the cell's narrowest width)", rc);
      rep.notes.push_back(b);
    }
  }
  const PairModel pm(ff, rc, o.coulomb, o.dsf_alpha);
  const Host H(s, pm, nh);
  std::vector<uint32_t> tpl;
  Vec3 c0{0, 0, 0};
  for (size_t i = nh; i < n; ++i) tpl.push_back(uint32_t(i)), c0 = c0 + s.atoms[i].pos;
  c0 = c0 * (1.0 / double(tpl.size()));
  std::vector<Vec3> shape;   // the template about its centre
  for (uint32_t i : tpl) shape.push_back(s.atoms[i].pos - c0);
  const double beta = 1 / (kB * o.temperature);
  const double V = s.cell.volume();
  rep.volume = V;
  for (size_t i = 0; i < nh; ++i) rep.host_mass += i < ff.mass.size() ? ff.mass[i] : s.mass_of(s.atoms[i]);
  std::mt19937_64 rng(o.seed);
  std::uniform_real_distribution<double> U(0.0, 1.0);
  auto random_quat = [&] {   // uniform on SO(3) (Shoemake 1992)
    const double u1 = U(rng), u2 = 2 * kPi * U(rng), u3 = 2 * kPi * U(rng);
    return Quat{std::sqrt(u1) * std::cos(u3), std::sqrt(1 - u1) * std::sin(u2), std::sqrt(1 - u1) * std::cos(u2), std::sqrt(u1) * std::sin(u3)};
  };
  auto place = [&](const Vec3& c, const Quat& q) {
    std::vector<Vec3> P(shape.size());
    for (size_t k = 0; k < shape.size(); ++k) P[k] = c + rotate(q, shape[k]);
    return P;
  };
  auto random_point = [&] { return s.cell.to_cartesian(Vec3{U(rng), U(rng), U(rng)}); };
  bool stopped = false;

  // ---- Widom
  if (o.insertions > 0) {
    const int nblk = 10;
    const int per = std::max(1, o.insertions / nblk);
    std::vector<double> blk(nblk, 0.0);
    for (int b = 0; b < nblk && !stopped; ++b) {
      double sum = 0;
      for (int i = 0; i < per; ++i) {
        const auto P = place(random_point(), random_quat());
        const double e = H.energy(tpl, P);
        sum += std::isfinite(e) ? std::exp(-beta * e) : 0.0;
      }
      blk[size_t(b)] = sum / per;
      if (o.progress && !o.progress("Widom insertion", double(b + 1) / nblk)) stopped = true;
    }
    double m = 0, v = 0;
    for (double x : blk) m += x;
    m /= nblk;
    for (double x : blk) v += (x - m) * (x - m);
    rep.widom_w = m;
    rep.widom_error = std::sqrt(v / (nblk - 1) / nblk);
    rep.mu_ex = m > 0 ? -std::log(m) / beta : std::numeric_limits<double>::infinity();
    const double Vm3 = V * 1e-30, mkg = rep.host_mass / kAvogadro * 1e-3;
    rep.henry_mol_kg_kpa = mkg > 0 ? m * Vm3 / (kBoltzmannJ * kAvogadro * o.temperature) / mkg * 1000.0 : 0.0;
    rep.solubility = m * 273.15 / o.temperature;
  }

  // ---- GCMC at each pressure
  for (size_t pi = 0; pi < o.pressures_kpa.size() && !stopped; ++pi) {
    const double p = o.pressures_kpa[pi];
    const double bfV = p * 1000.0 * V * 1e-30 / (kBoltzmannJ * o.temperature);   // β f V, ideal-gas fugacity
    std::vector<std::vector<Vec3>> mols;
    auto e_of = [&](const std::vector<Vec3>& P, size_t skip) {
      double e = H.energy(tpl, P);
      if (!std::isfinite(e)) return e;
      for (size_t j = 0; j < mols.size(); ++j) if (j != skip) e += H.between(tpl, P, mols[j]);
      return e;
    };
    double Ucur = 0;
    long ins_try = 0, ins_ok = 0, del_try = 0, del_ok = 0;
    const int equil = o.steps / 4;
    double sN = 0, sN2 = 0, sU = 0, sUN = 0;
    long samples = 0;
    const int nblk = 10;
    std::vector<double> blkN(nblk, 0.0), blkC(nblk, 0.0);
    const int report_every = std::max(1, o.steps / 50);
    for (int step = 0; step < o.steps && !stopped; ++step) {
      const double r = U(rng);
      const size_t N = mols.size();
      if (r < 0.25) {   // insert
        ++ins_try;
        auto P = place(random_point(), random_quat());
        const double dU = e_of(P, size_t(-1));
        if (std::isfinite(dU) && U(rng) < bfV / double(N + 1) * std::exp(-beta * dU)) { mols.push_back(std::move(P)); Ucur += dU; ++ins_ok; }
      } else if (r < 0.5) {   // delete
        if (N > 0) {
          ++del_try;
          const size_t m = size_t(U(rng) * double(N)) % N;
          const double dU = -e_of(mols[m], m);
          if (U(rng) < double(N) / bfV * std::exp(-beta * dU)) { mols[m] = std::move(mols.back()); mols.pop_back(); Ucur += dU; ++del_ok; }
        }
      } else if (N > 0) {   // translate or rotate one molecule
        const size_t m = size_t(U(rng) * double(N)) % N;
        std::vector<Vec3> Q = mols[m];
        if (r < 0.75) {
          const Vec3 d{(2 * U(rng) - 1) * 1.0, (2 * U(rng) - 1) * 1.0, (2 * U(rng) - 1) * 1.0};
          for (auto& q : Q) q = q + d;
        } else {
          Vec3 c{0, 0, 0};
          for (const auto& q : Q) c = c + q;
          c = c * (1.0 / double(Q.size()));
          const Vec3 ax = [&] { const double z = 2 * U(rng) - 1, t = 2 * kPi * U(rng), rr = std::sqrt(std::max(0.0, 1 - z * z)); return Vec3{rr * std::cos(t), rr * std::sin(t), z}; }();
          const double a = (2 * U(rng) - 1) * 30 * kPi / 180;
          const Quat q{std::cos(a / 2), ax[0] * std::sin(a / 2), ax[1] * std::sin(a / 2), ax[2] * std::sin(a / 2)};
          for (auto& x : Q) x = c + rotate(q, x - c);
        }
        const double e_old = e_of(mols[m], m), e_new = e_of(Q, m);
        if (std::isfinite(e_new) && (e_new <= e_old || U(rng) < std::exp(-beta * (e_new - e_old)))) { mols[m] = std::move(Q); Ucur += e_new - e_old; }
      }
      if (step >= equil) {
        const double Nn = double(mols.size());
        sN += Nn, sN2 += Nn * Nn, sU += Ucur, sUN += Ucur * Nn, ++samples;
        const int b = std::min(nblk - 1, int(double(step - equil) / double(o.steps - equil) * nblk));
        blkN[size_t(b)] += Nn, blkC[size_t(b)] += 1;
      }
      if (o.progress && (step + 1) % report_every == 0) {
        char st[80];
        std::snprintf(st, sizeof st, "GCMC at %g kPa · %zu molecules", p, mols.size());
        if (!o.progress(st, (pi + double(step + 1) / o.steps) / o.pressures_kpa.size())) stopped = true;
      }
    }
    IsothermPoint pt;
    pt.pressure_kpa = p;
    if (samples > 0) {
      const double mN = sN / samples, mN2 = sN2 / samples, mU = sU / samples, mUN = sUN / samples;
      pt.loading = mN;
      double bm = 0, bv = 0;
      int nbk = 0;
      for (int b = 0; b < nblk; ++b) if (blkC[size_t(b)] > 0) bm += blkN[size_t(b)] / blkC[size_t(b)], ++nbk;
      bm /= std::max(1, nbk);
      for (int b = 0; b < nblk; ++b) if (blkC[size_t(b)] > 0) { const double x = blkN[size_t(b)] / blkC[size_t(b)] - bm; bv += x * x; }
      pt.loading_error = nbk > 1 ? std::sqrt(bv / (nbk - 1) / nbk) : 0;
      const double varN = mN2 - mN * mN;
      pt.heat = varN > 1e-9 ? 1 / beta - (mUN - mU * mN) / varN : 0;
      const double mkg = rep.host_mass / kAvogadro * 1e-3;
      pt.mol_per_kg = mkg > 0 ? mN / kAvogadro / mkg : 0;
      pt.cm3stp_per_cm3 = mN / (V * 1e-24) * 22413.969 / kAvogadro;   // molecules per cm³ → cm³(STP) per cm³
    }
    pt.acceptance_insert = ins_try ? double(ins_ok) / ins_try : 0;
    pt.acceptance_delete = del_try ? double(del_ok) / del_try : 0;
    pt.molecules = mols;
    rep.isotherm.push_back(std::move(pt));
  }
  if (stopped) rep.notes.push_back("stopped before the end");
  {
    char b[240];
    std::snprintf(b, sizeof b, "host %zu atoms (%.0f g/mol per cell), cell %.0f Å³ · sorbate %zu atoms, rigid · %.0f K", nh, rep.host_mass, V, tpl.size(), o.temperature);
    rep.notes.insert(rep.notes.begin(), b);
    rep.notes.push_back(std::string("energies: the force field's van der Waals") + (pm.coulomb() ? " and damped shifted force electrostatics" : "") +
                        "; the host is held fixed; ideal-gas reservoir (fugacity = pressure)");
  }
  rep.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  return rep;
}

}  // namespace caps
