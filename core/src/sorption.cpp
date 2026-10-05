// CAPS sorption: Widom insertion and grand-canonical Monte Carlo in a fixed host (see sorption.hpp).
#include "caps/rng.hpp"
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
  double between(const std::vector<uint32_t>& tp, const std::vector<Vec3>& P, const std::vector<uint32_t>& tq, const std::vector<Vec3>& Q) const {
    double e = 0;
    for (size_t k = 0; k < tp.size(); ++k)
      for (size_t q = 0; q < tq.size(); ++q) {
        const Vec3 d = mi(Q[q] - P[k]);
        e += pm.energy(tp[k], tq[q], dot(d, d));
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
  // the species: each template's atoms and its shape about its centre
  std::vector<int> starts = o.species_first_atom.empty() ? std::vector<int>{o.template_first_atom} : o.species_first_atom;
  for (size_t k = 0; k < starts.size(); ++k)
    if (starts[k] < int(nh) || size_t(starts[k]) >= n || (k > 0 && starts[k] <= starts[k - 1]))
      throw std::invalid_argument("sorption: the species' templates must follow the host in order");
  const size_t ns = starts.size();
  std::vector<std::vector<uint32_t>> tpls(ns);
  std::vector<std::vector<Vec3>> shapes(ns);
  for (size_t k = 0; k < ns; ++k) {
    const size_t a = size_t(starts[k]), b = k + 1 < ns ? size_t(starts[k + 1]) : n;
    Vec3 c0{0, 0, 0};
    for (size_t i = a; i < b; ++i) tpls[k].push_back(uint32_t(i)), c0 = c0 + s.atoms[i].pos;
    c0 = c0 * (1.0 / double(b - a));
    for (size_t i = a; i < b; ++i) shapes[k].push_back(s.atoms[i].pos - c0);
  }
  std::vector<double> y(ns, 1.0 / double(ns));
  if (!o.mole_fractions.empty()) {
    if (o.mole_fractions.size() != ns) throw std::invalid_argument("sorption: one mole fraction per species");
    double t = 0;
    for (double v : o.mole_fractions) { if (v < 0) throw std::invalid_argument("sorption: mole fractions ≥ 0"); t += v; }
    if (t <= 0) throw std::invalid_argument("sorption: the mole fractions add up to zero");
    for (size_t k = 0; k < ns; ++k) y[k] = o.mole_fractions[k] / t;
  }
  const auto& tpl = tpls[0];
  const double beta = 1 / (kB * o.temperature);
  const double V = s.cell.volume();
  rep.volume = V;
  for (size_t i = 0; i < nh; ++i) rep.host_mass += i < ff.mass.size() ? ff.mass[i] : s.mass_of(s.atoms[i]);
  std::mt19937_64 rng(o.seed);
  caps::UniformReal<double> U(0.0, 1.0);
  auto random_quat = [&] {   // uniform on SO(3) (Shoemake 1992)
    const double u1 = U(rng), u2 = 2 * kPi * U(rng), u3 = 2 * kPi * U(rng);
    return Quat{std::sqrt(u1) * std::cos(u3), std::sqrt(1 - u1) * std::sin(u2), std::sqrt(1 - u1) * std::cos(u2), std::sqrt(u1) * std::sin(u3)};
  };
  auto place = [&](size_t sp, const Vec3& c, const Quat& q) {
    std::vector<Vec3> P(shapes[sp].size());
    for (size_t k = 0; k < P.size(); ++k) P[k] = c + rotate(q, shapes[sp][k]);
    return P;
  };
  auto random_point = [&] { return s.cell.to_cartesian(Vec3{U(rng), U(rng), U(rng)}); };
  bool stopped = false;

  // ---- Widom, per species
  const double Vm3 = V * 1e-30, mkg = rep.host_mass / kAvogadro * 1e-3;
  for (size_t sp = 0; sp < ns && o.insertions > 0; ++sp) {
    const int nblk = 10;
    const int per = std::max(1, o.insertions / nblk);
    std::vector<double> blk(nblk, 0.0);
    for (int b = 0; b < nblk && !stopped; ++b) {
      double sum = 0;
      for (int i = 0; i < per; ++i) {
        const auto P = place(sp, random_point(), random_quat());
        const double e = H.energy(tpls[sp], P);
        sum += std::isfinite(e) ? std::exp(-beta * e) : 0.0;
      }
      blk[size_t(b)] = sum / per;
      if (o.progress && !o.progress(ns > 1 ? "Widom insertion, species " + std::to_string(sp + 1) : std::string("Widom insertion"), (double(sp) + double(b + 1) / nblk) / double(ns)))
        stopped = true;
    }
    double m = 0, v = 0;
    for (double x : blk) m += x;
    m /= nblk;
    for (double x : blk) v += (x - m) * (x - m);
    SpeciesWidom sw;
    sw.fraction = y[sp];
    sw.widom_w = m;
    sw.widom_error = std::sqrt(v / (nblk - 1) / nblk);
    sw.mu_ex = m > 0 ? -std::log(m) / beta : std::numeric_limits<double>::infinity();
    sw.henry_mol_kg_kpa = mkg > 0 ? m * Vm3 / (kBoltzmannJ * kAvogadro * o.temperature) / mkg * 1000.0 : 0.0;
    sw.solubility = m * 273.15 / o.temperature;
    rep.species.push_back(sw);
  }
  if (!rep.species.empty()) {
    const auto& f = rep.species.front();
    rep.widom_w = f.widom_w, rep.widom_error = f.widom_error, rep.mu_ex = f.mu_ex, rep.henry_mol_kg_kpa = f.henry_mol_kg_kpa, rep.solubility = f.solubility;
  } else
    for (size_t sp = 0; sp < ns; ++sp) rep.species.push_back(SpeciesWidom{y[sp]});

  // ---- GCMC at each pressure
  struct Mol { size_t sp; std::vector<Vec3> P; };
  for (size_t pi = 0; pi < o.pressures_kpa.size() && !stopped; ++pi) {
    const double p = o.pressures_kpa[pi];
    std::vector<double> bfV(ns);   // β f_i V per species, ideal-gas fugacity y_i p
    for (size_t k = 0; k < ns; ++k) bfV[k] = y[k] * p * 1000.0 * V * 1e-30 / (kBoltzmannJ * o.temperature);
    std::vector<Mol> mols;
    std::vector<std::vector<size_t>> of(ns);   // the molecules of each species (indices into mols)
    auto e_of = [&](size_t sp, const std::vector<Vec3>& P, size_t skip) {
      double e = H.energy(tpls[sp], P);
      if (!std::isfinite(e)) return e;
      for (size_t j = 0; j < mols.size(); ++j) if (j != skip) e += H.between(tpls[sp], P, tpls[mols[j].sp], mols[j].P);
      return e;
    };
    auto remove = [&](size_t m) {   // mols[m] out: the last takes its place
      const size_t last = mols.size() - 1;
      auto& a = of[mols[m].sp];
      a.erase(std::find(a.begin(), a.end(), m));
      if (m != last) {
        auto& b = of[mols[last].sp];
        *std::find(b.begin(), b.end(), last) = m;
        mols[m] = std::move(mols[last]);
      }
      mols.pop_back();
    };
    double Ucur = 0;
    long ins_try = 0, ins_ok = 0, del_try = 0, del_ok = 0;
    const int equil = o.steps / 4;
    // sums for the averages and the fluctuations: N_i, N_i N_j, U, U N_i
    std::vector<double> sN(ns, 0.0), sUN(ns, 0.0), sNN(ns * ns, 0.0);
    double sU = 0;
    long samples = 0;
    const int nblk = 10;
    std::vector<double> blkC(nblk, 0.0), blkN(size_t(nblk) * ns, 0.0);
    const int report_every = std::max(1, o.steps / 50);
    const int g = std::clamp(o.map_grid, 0, 96);
    std::vector<double> dens(size_t(g) * g * g, 0.0);
    long map_samples = 0;
    const int map_every = std::max(1, (o.steps - equil) / 2000);
    for (int step = 0; step < o.steps && !stopped; ++step) {
      const double r = U(rng);
      const size_t N = mols.size();
      if (r < 0.25) {   // insert one molecule of a species picked at random
        ++ins_try;
        const size_t sp = std::min(ns - 1, size_t(U(rng) * double(ns)));
        auto P = place(sp, random_point(), random_quat());
        const double dU = e_of(sp, P, size_t(-1));
        if (std::isfinite(dU) && U(rng) < bfV[sp] / double(of[sp].size() + 1) * std::exp(-beta * dU)) {
          of[sp].push_back(mols.size());
          mols.push_back({sp, std::move(P)});
          Ucur += dU;
          ++ins_ok;
        }
      } else if (r < 0.5) {   // delete one of a species picked at random
        const size_t sp = std::min(ns - 1, size_t(U(rng) * double(ns)));
        const size_t Ns = of[sp].size();
        if (Ns > 0) {
          ++del_try;
          const size_t m = of[sp][size_t(U(rng) * double(Ns)) % Ns];
          const double dU = -e_of(sp, mols[m].P, m);
          if (bfV[sp] <= 0 || U(rng) < double(Ns) / bfV[sp] * std::exp(-beta * dU)) { remove(m); Ucur += dU; ++del_ok; }
        }
      } else if (N > 0) {   // translate or rotate one molecule
        const size_t m = size_t(U(rng) * double(N)) % N;
        const size_t sp = mols[m].sp;
        std::vector<Vec3> Q = mols[m].P;
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
        const double e_old = e_of(sp, mols[m].P, m), e_new = e_of(sp, Q, m);
        if (std::isfinite(e_new) && (e_new <= e_old || U(rng) < std::exp(-beta * (e_new - e_old)))) { mols[m].P = std::move(Q); Ucur += e_new - e_old; }
      }
      if (step >= equil) {
        sU += Ucur, ++samples;
        const int b = std::min(nblk - 1, int(double(step - equil) / double(o.steps - equil) * nblk));
        blkC[size_t(b)] += 1;
        for (size_t i = 0; i < ns; ++i) {
          const double Ni = double(of[i].size());
          sN[i] += Ni, sUN[i] += Ucur * Ni, blkN[size_t(b) * ns + i] += Ni;
          for (size_t j = 0; j < ns; ++j) sNN[i * ns + j] += Ni * double(of[j].size());
        }
        if (g > 0 && (step - equil) % map_every == 0) {   // the molecules' centres on the map (every species)
          for (const auto& mol : mols) {
            Vec3 cc{0, 0, 0};
            for (const auto& q : mol.P) cc = cc + q;
            Vec3 fr = s.cell.to_fractional(cc * (1.0 / double(mol.P.size())));
            int ix[3];
            for (int k = 0; k < 3; ++k) ix[k] = std::clamp(int(std::floor((fr[k] - std::floor(fr[k])) * g)), 0, g - 1);
            dens[(size_t(ix[0]) * g + size_t(ix[1])) * g + size_t(ix[2])] += 1;
          }
          ++map_samples;
        }
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
      const double inv = 1.0 / double(samples), mU = sU * inv;
      std::vector<double> mN(ns), cUN(ns), cNN(ns * ns);
      for (size_t i = 0; i < ns; ++i) mN[i] = sN[i] * inv;
      for (size_t i = 0; i < ns; ++i) {
        cUN[i] = sUN[i] * inv - mU * mN[i];
        for (size_t j = 0; j < ns; ++j) cNN[i * ns + j] = sNN[i * ns + j] * inv - mN[i] * mN[j];
      }
      pt.species_heat.assign(ns, 0.0);
      {   // q_i = kT − Σ_j cov(U, N_j) [cov(N, N)⁻¹]_ji: solve cov(N, N) z = cov(U, N) (symmetric), Gauss with pivoting; species never adsorbed left out
        std::vector<size_t> act;
        for (size_t i = 0; i < ns; ++i) if (cNN[i * ns + i] > 1e-9) act.push_back(i);
        const size_t m = act.size();
        std::vector<double> A(m * m), z(m);
        for (size_t a = 0; a < m; ++a) {
          z[a] = cUN[act[a]];
          for (size_t b = 0; b < m; ++b) A[a * m + b] = cNN[act[a] * ns + act[b]];
        }
        bool ok = true;
        for (size_t c = 0; c < m && ok; ++c) {
          size_t piv = c;
          for (size_t r2 = c + 1; r2 < m; ++r2) if (std::fabs(A[r2 * m + c]) > std::fabs(A[piv * m + c])) piv = r2;
          if (std::fabs(A[piv * m + c]) < 1e-12) { ok = false; break; }
          for (size_t k = 0; k < m; ++k) std::swap(A[c * m + k], A[piv * m + k]);
          std::swap(z[c], z[piv]);
          for (size_t r2 = 0; r2 < m; ++r2)
            if (r2 != c) {
              const double f = A[r2 * m + c] / A[c * m + c];
              for (size_t k = 0; k < m; ++k) A[r2 * m + k] -= f * A[c * m + k];
              z[r2] -= f * z[c];
            }
        }
        if (ok)
          for (size_t a = 0; a < m; ++a) pt.species_heat[act[a]] = 1 / beta - z[a] / A[a * m + a];
      }
      double tot = 0;
      for (double v : mN) tot += v;
      pt.loading = tot;
      pt.species_loading = mN;
      pt.species_error.assign(ns, 0.0);
      // the error of the total and of each species: ten blocks
      auto blk_err = [&](auto value) {
        double bm = 0, bv = 0;
        int nbk = 0;
        for (int b = 0; b < nblk; ++b) if (blkC[size_t(b)] > 0) bm += value(b), ++nbk;
        bm /= std::max(1, nbk);
        for (int b = 0; b < nblk; ++b) if (blkC[size_t(b)] > 0) { const double x = value(b) - bm; bv += x * x; }
        return nbk > 1 ? std::sqrt(bv / (nbk - 1) / nbk) : 0.0;
      };
      pt.loading_error = blk_err([&](int b) { double t = 0; for (size_t i = 0; i < ns; ++i) t += blkN[size_t(b) * ns + i]; return t / blkC[size_t(b)]; });
      for (size_t i = 0; i < ns; ++i) pt.species_error[i] = blk_err([&](int b) { return blkN[size_t(b) * ns + i] / blkC[size_t(b)]; });
      if (ns == 1) {   // one species: the single-component formula, as before
        pt.heat = pt.species_heat[0];
      } else {         // the adsorbed phase's mean: Σ x_i q_i
        pt.heat = 0;
        if (tot > 0) for (size_t i = 0; i < ns; ++i) pt.heat += mN[i] / tot * pt.species_heat[i];
      }
      pt.mol_per_kg = mkg > 0 ? tot / kAvogadro / mkg : 0;
      pt.species_mol_per_kg.resize(ns);
      for (size_t i = 0; i < ns; ++i) pt.species_mol_per_kg[i] = mkg > 0 ? mN[i] / kAvogadro / mkg : 0;
      pt.cm3stp_per_cm3 = tot / (V * 1e-24) * 22413.969 / kAvogadro;   // molecules per cm³ → cm³(STP) per cm³
      pt.selectivity.assign(ns, 0.0);
      for (size_t i = 0; i < ns; ++i)   // (x_i / x_0) / (y_i / y_0)
        if (mN[0] > 0 && y[i] > 0) pt.selectivity[i] = (mN[i] / mN[0]) / (y[i] / y[0]);
    }
    pt.acceptance_insert = ins_try ? double(ins_ok) / ins_try : 0;
    pt.acceptance_delete = del_try ? double(del_ok) / del_try : 0;
    for (const auto& mol : mols) pt.molecules.push_back(mol.P);
    if (g > 0 && map_samples > 0) {   // per Å³: counts over the samples and the voxel's volume
      const double vox = V / (double(g) * g * g);
      pt.grid = g;
      pt.density.resize(dens.size());
      for (size_t k = 0; k < dens.size(); ++k) pt.density[k] = float(dens[k] / (double(map_samples) * vox));
    }
    rep.isotherm.push_back(std::move(pt));
  }
  if (stopped) rep.notes.push_back("stopped before the end");
  {
    char b[240];
    if (ns == 1)
      std::snprintf(b, sizeof b, "host %zu atoms (%.0f g/mol per cell), cell %.0f Å³ · sorbate %zu atoms, rigid · %.0f K", nh, rep.host_mass, V, tpl.size(), o.temperature);
    else
      std::snprintf(b, sizeof b, "host %zu atoms (%.0f g/mol per cell), cell %.0f Å³ · %zu sorbates, rigid, ideal gas mixture · %.0f K", nh, rep.host_mass, V, ns, o.temperature);
    rep.notes.insert(rep.notes.begin(), b);
    rep.notes.push_back(std::string("energies: the force field's van der Waals") + (pm.coulomb() ? " and damped shifted force electrostatics" : "") +
                        "; the host is held fixed; ideal-gas reservoir (fugacity = " + std::string(ns > 1 ? "y_i × pressure)" : "pressure)"));
  }
  rep.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  return rep;
}

}  // namespace caps
