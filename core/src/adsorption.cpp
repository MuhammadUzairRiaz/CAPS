// CAPS adsorption locator (see adsorption.hpp).
#include "caps/adsorption.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <map>
#include <random>
#include <set>
#include <stdexcept>

#include "caps/elements.hpp"
#include "caps/pairmodel.hpp"

namespace caps {

namespace {

constexpr double kB = 0.0019872043;   // kcal/(mol·K)
constexpr double kPi = 3.14159265358979323846;

Vec3 rotate(const Vec3& v, const Vec3& u, double t) {   // v about unit axis u by t (Rodrigues)
  const double c = std::cos(t), s = std::sin(t);
  return v * c + cross(u, v) * s + u * (dot(u, v) * (1 - c));
}

// The substrate (fixed atoms) and the mobile molecules, with the energies the locator counts.
struct Model {
  const System& s;
  const ForceField& ff;
  PairModel pm;
  const Cell& cell;
  bool per = false, ortho = false;
  double L[3] = {0, 0, 0};
  bool pdim[3] = {false, false, false};
  std::vector<std::vector<uint32_t>> mol;   // atoms of each mobile molecule
  std::vector<uint32_t> fixed;
  // bins of the fixed atoms, at least rc/2 wide, two either side (fewer than five: every bin once)
  int nb[3] = {1, 1, 1};
  std::vector<std::vector<uint32_t>> bins;
  std::vector<std::vector<int>> stencil;

  Model(const System& sys, const ForceField& f, const std::vector<int>& mobile, double rc, bool coul, double alpha)
      : s(sys), ff(f), pm(f, rc, coul, alpha), cell(sys.cell) {
    per = cell.valid();
    ortho = per && std::fabs(cell.a[1]) + std::fabs(cell.a[2]) + std::fabs(cell.b[0]) + std::fabs(cell.b[2]) + std::fabs(cell.c[0]) + std::fabs(cell.c[1]) < 1e-9;
    L[0] = cell.a[0], L[1] = cell.b[1], L[2] = cell.c[2];
    for (int d = 0; d < 3; ++d) pdim[d] = cell.periodic[size_t(d)];
    int nm = 0;
    const auto comp = sys.molecules(&nm);
    std::map<int, size_t> slot;
    for (size_t k = 0; k < mobile.size(); ++k) slot[mobile[k]] = k;
    mol.resize(mobile.size());
    for (uint32_t i = 0; i < sys.atoms.size(); ++i) {
      const auto it = slot.find(comp[i]);
      if (it == slot.end()) fixed.push_back(i);
      else mol[it->second].push_back(i);
    }
    if (per) {
      const double V = cell.volume();
      const double w[3] = {V / norm(cross(cell.b, cell.c)), V / norm(cross(cell.c, cell.a)), V / norm(cross(cell.a, cell.b))};
      for (int d = 0; d < 3; ++d)
        if (pdim[d]) nb[d] = std::max(1, int(std::floor(2 * w[d] / rc)));
    }
    bins.resize(size_t(nb[0] * nb[1] * nb[2]));
    for (uint32_t i : fixed) bins[size_t(bin_of(sys.atoms[i].pos))].push_back(i);
    stencil.resize(bins.size());
    for (int a = 0; a < nb[0]; ++a)
      for (int b = 0; b < nb[1]; ++b)
        for (int c = 0; c < nb[2]; ++c) {
          auto offs = [&](int d, int v) {
            std::vector<int> r;
            if (nb[d] < 5) for (int k = 0; k < nb[d]; ++k) r.push_back(k);
            else for (int dv = -2; dv <= 2; ++dv) r.push_back(((v + dv) % nb[d] + nb[d]) % nb[d]);
            return r;
          };
          for (int i : offs(0, a))
            for (int j : offs(1, b))
              for (int k : offs(2, c)) stencil[size_t((a * nb[1] + b) * nb[2] + c)].push_back((i * nb[1] + j) * nb[2] + k);
        }
  }
  Vec3 mi(Vec3 d) const {
    if (!per) return d;
    if (!ortho) return cell.minimum_image(d);
    for (int k = 0; k < 3; ++k)
      if (pdim[k]) d[size_t(k)] -= L[k] * std::round(d[size_t(k)] / L[k]);
    return d;
  }
  int bin_of(const Vec3& p) const {
    if (!per) return 0;
    const Vec3 f = cell.to_fractional(p);
    int id[3];
    for (int d = 0; d < 3; ++d) id[d] = std::min(nb[d] - 1, int((f[size_t(d)] - std::floor(f[size_t(d)])) * nb[d]));
    return (id[0] * nb[1] + id[1]) * nb[2] + id[2];
  }
  // molecule m at positions P (one per atom of m) with the substrate
  double with_fixed(size_t m, const std::vector<Vec3>& P) const {
    double e = 0;
    for (size_t k = 0; k < mol[m].size(); ++k) {
      const uint32_t a = mol[m][k];
      for (int bn : stencil[size_t(bin_of(P[k]))])
        for (uint32_t b : bins[size_t(bn)]) {
          const Vec3 d = mi(s.atoms[b].pos - P[k]);
          e += pm.energy(a, b, dot(d, d));
        }
      if (!std::isfinite(e)) return e;
    }
    return e;
  }
  double between(size_t m, const std::vector<Vec3>& P, size_t n, const std::vector<Vec3>& Q) const {
    double e = 0;
    for (size_t k = 0; k < mol[m].size(); ++k)
      for (size_t q = 0; q < mol[n].size(); ++q) {
        const Vec3 d = mi(Q[q] - P[k]);
        e += pm.energy(mol[m][k], mol[n][q], dot(d, d));
      }
    return e;
  }
};

std::string formula(const System& s, const std::vector<uint32_t>& atoms) {   // Hill order
  std::map<std::string, int> n;
  for (uint32_t a : atoms) ++n[element(s.atoms[a].element).symbol];
  std::string f;
  auto put = [&](const std::string& el) {
    const auto it = n.find(el);
    if (it == n.end()) return;
    f += el + (it->second > 1 ? std::to_string(it->second) : "");
    n.erase(it);
  };
  if (n.count("C")) { put("C"); put("H"); }
  for (const auto& [el, c] : std::map<std::string, int>(n)) put(el);
  return f;
}

std::vector<int> default_mobile(const System& s) {   // every molecule but the largest
  int nm = 0;
  const auto comp = s.molecules(&nm);
  std::vector<int> size(size_t(nm), 0);
  for (int c : comp) ++size[size_t(c)];
  const int big = int(std::max_element(size.begin(), size.end()) - size.begin());
  std::vector<int> out;
  for (int m = 0; m < nm; ++m) if (m != big) out.push_back(m);
  return out;
}

}  // namespace

std::vector<int> adsorbate_molecules(const System& s, int first_atom) {
  if (first_atom < 0) return default_mobile(s);
  int nm = 0;
  const auto comp = s.molecules(&nm);
  std::set<int> m;
  for (size_t i = size_t(first_atom); i < s.atoms.size(); ++i) m.insert(comp[i]);
  return {m.begin(), m.end()};
}

double adsorption_interaction(const System& s, const ForceField& ff, const std::vector<int>& mobile_in, double cutoff, bool coulomb, double dsf_alpha) {
  const auto mobile = mobile_in.empty() ? default_mobile(s) : mobile_in;
  Model M(s, ff, mobile, cutoff, coulomb, dsf_alpha);
  std::vector<std::vector<Vec3>> P(M.mol.size());
  for (size_t m = 0; m < M.mol.size(); ++m) for (uint32_t a : M.mol[m]) P[m].push_back(s.atoms[a].pos);
  double e = 0;
  for (size_t m = 0; m < M.mol.size(); ++m) {
    e += M.with_fixed(m, P[m]);
    for (size_t n = m + 1; n < M.mol.size(); ++n) e += M.between(m, P[m], n, P[n]);
  }
  return e;
}

void locate_adsorption(System& s, const ForceField& ff, const AdsorptionOptions& o, AdsorptionReport* report) {
  const auto t0 = std::chrono::steady_clock::now();
  AdsorptionReport rep;
  const size_t n = s.atoms.size();
  if (ff.type_index.size() != n || ff.charge.size() != n) throw std::invalid_argument("adsorption: the force field was assigned to another structure");
  if (!ff.pair_func.empty() || ff.hbond.on() || ff.sw.on || !ff.vsites.empty())
    throw std::invalid_argument("adsorption: the energies need Lennard-Jones pairs (no other pair forms, hydrogen bonds or virtual sites)");
  if (o.steps < 10 || o.cycles < 1 || o.t_high <= 0 || o.t_low <= 0 || o.t_low > o.t_high) throw std::invalid_argument("adsorption: at least 10 steps and one cycle, 0 < t_low ≤ t_high");
  const auto mobile = o.mobile.empty() ? adsorbate_molecules(s, o.first_mobile_atom) : o.mobile;
  if (mobile.empty()) throw std::invalid_argument("adsorption: no adsorbate molecule (the structure is one molecule)");
  double rc = o.cutoff;
  if (s.cell.valid()) {
    const double V = s.cell.volume();
    double wmin = 1e30;
    const double w[3] = {V / norm(cross(s.cell.b, s.cell.c)), V / norm(cross(s.cell.c, s.cell.a)), V / norm(cross(s.cell.a, s.cell.b))};
    for (int d = 0; d < 3; ++d) if (s.cell.periodic[size_t(d)]) wmin = std::min(wmin, w[d]);
    if (rc > 0.5 * wmin - 1e-6) {
      rc = 0.5 * wmin - 1e-6;
      char b[160];
      std::snprintf(b, sizeof b, "cut-off reduced to %.2f Å (half the cell's narrowest width)", rc);
      rep.notes.push_back(b);
    }
  }
  Model M(s, ff, mobile, rc, o.coulomb, o.dsf_alpha);
  const size_t nm = M.mol.size();
  for (size_t m = 0; m < nm; ++m)
    if (M.mol[m].empty()) throw std::invalid_argument("adsorption: molecule " + std::to_string(mobile[m]) + " does not exist");
  if (M.fixed.empty()) rep.notes.push_back("no substrate: every molecule moves (a cluster search)");

  std::mt19937_64 rng(o.seed);
  std::uniform_real_distribution<double> U(0.0, 1.0);
  std::vector<std::vector<Vec3>> P(nm);
  for (size_t m = 0; m < nm; ++m) for (uint32_t a : M.mol[m]) P[m].push_back(s.atoms[a].pos);
  const bool zband = o.z_hi > o.z_lo;
  auto inside = [&](const std::vector<Vec3>& Q) {
    if (!zband) return true;
    for (const auto& q : Q) if (q[2] < o.z_lo || q[2] > o.z_hi) return false;
    return true;
  };
  auto centre = [](const std::vector<Vec3>& Q) { Vec3 c{0, 0, 0}; for (const auto& q : Q) c = c + q; return c * (1.0 / double(Q.size())); };
  auto fold = [&](std::vector<Vec3>& Q) {   // the centre back into the cell (whole molecule shifted)
    if (!M.per) return;
    const Vec3 c = centre(Q), d = s.cell.wrap(c) - c;
    for (auto& q : Q) q = q + d;
  };
  auto random_axis = [&] {
    const double z = 2 * U(rng) - 1, t = 2 * kPi * U(rng), r = std::sqrt(std::max(0.0, 1 - z * z));
    return Vec3{r * std::cos(t), r * std::sin(t), z};
  };
  auto random_place = [&](const std::vector<Vec3>& Q) {   // a uniform place and a uniform orientation
    const Vec3 c = centre(Q);
    const Vec3 u = random_axis();
    const double ang = 2 * kPi * U(rng);
    Vec3 target;
    if (M.per) {
      Vec3 f{U(rng), U(rng), U(rng)};
      target = s.cell.to_cartesian(f);
      if (zband) target[2] = o.z_lo + (o.z_hi - o.z_lo) * U(rng);
    } else {
      Vec3 lo{1e30, 1e30, 1e30}, hi{-1e30, -1e30, -1e30};
      for (uint32_t b : M.fixed) for (int k = 0; k < 3; ++k) lo[size_t(k)] = std::min(lo[size_t(k)], s.atoms[b].pos[size_t(k)]), hi[size_t(k)] = std::max(hi[size_t(k)], s.atoms[b].pos[size_t(k)]);
      if (M.fixed.empty()) lo = c - Vec3{10, 10, 10}, hi = c + Vec3{10, 10, 10};
      for (int k = 0; k < 3; ++k) target[size_t(k)] = lo[size_t(k)] - 5 + (hi[size_t(k)] - lo[size_t(k)] + 10) * U(rng);
      if (zband) target[2] = o.z_lo + (o.z_hi - o.z_lo) * U(rng);
    }
    std::vector<Vec3> R(Q.size());
    for (size_t k = 0; k < Q.size(); ++k) R[k] = target + rotate(Q[k] - c, u, ang);
    return R;
  };
  auto mol_energy = [&](size_t m, const std::vector<Vec3>& Q) {
    double e = M.with_fixed(m, Q);
    if (!std::isfinite(e)) return e;
    for (size_t k = 0; k < nm; ++k) if (k != m) e += M.between(m, Q, k, P[k]);
    return e;
  };
  auto total = [&] {
    double e = 0;
    for (size_t m = 0; m < nm; ++m) {
      e += M.with_fixed(m, P[m]);
      for (size_t k = m + 1; k < nm; ++k) e += M.between(m, P[m], k, P[k]);
    }
    return e;
  };
  // start: each molecule inserted in turn at the best of 50 random places
  if (o.randomise)
    for (size_t m = 0; m < nm; ++m) {
      std::vector<Vec3> best = P[m];
      double be = std::numeric_limits<double>::infinity();
      for (int t = 0; t < 50; ++t) {
        auto Q = random_place(P[m]);
        if (!inside(Q)) continue;
        fold(Q);
        const double e = mol_energy(m, Q);
        if (e < be) be = e, best = Q;
      }
      P[m] = best;
    }

  std::vector<AdsorptionConfig> finals;
  std::vector<std::vector<Vec3>> bestP = P;
  double best_e = total();
  std::vector<double> samples;
  long acc = 0, att = 0;
  const int every = std::max(1, o.steps / 100);
  const double ratio = std::log(o.t_high / o.t_low);
  bool stopped = false;
  for (int cyc = 0; cyc < o.cycles && !stopped; ++cyc) {
    if (cyc > 0) P = bestP;
    double cur = total();
    double cyc_best = cur;
    std::vector<std::vector<Vec3>> cyc_bestP = P;
    auto step = [&](double T, double dmax, double amax, double pjump, bool cold) {
      const size_t m = size_t(U(rng) * double(nm)) % nm;
      const double r = U(rng);
      std::vector<Vec3> Q = P[m];
      if (r < pjump) Q = random_place(P[m]);
      else if (r < pjump + 0.5 * (1 - pjump)) {   // translate
        const Vec3 d{(2 * U(rng) - 1) * dmax, (2 * U(rng) - 1) * dmax, (2 * U(rng) - 1) * dmax};
        for (auto& q : Q) q = q + d;
      } else {   // rotate about the centre
        const Vec3 c = centre(Q), u = random_axis();
        const double a = (2 * U(rng) - 1) * amax;
        for (auto& q : Q) q = c + rotate(q - c, u, a);
      }
      ++att;
      if (!inside(Q)) return;
      fold(Q);
      const double e_old = mol_energy(m, P[m]), e_new = mol_energy(m, Q);
      if (!std::isfinite(e_new)) return;
      const double dE = e_new - (std::isfinite(e_old) ? e_old : 1e300);
      if (dE <= 0 || (!cold && U(rng) < std::exp(-dE / (kB * T)))) {
        P[m] = std::move(Q);
        cur += std::isfinite(e_old) ? dE : 0;
        if (!std::isfinite(e_old)) cur = total();
        ++acc;
        if (cur < best_e) best_e = cur, bestP = P;
      }
    };
    for (int k = 0; k < o.steps; ++k) {
      const double f = 1.0 - double(k) / std::max(1, o.steps - 1);    // 1 hot … 0 cold
      const double T = o.t_low * std::exp(ratio * f);
      step(T, 0.2 + 2.8 * f, (5 + 55 * f) * kPi / 180, 0.3 * f, false);   // jumps fade out as it cools
      if (cur < cyc_best) cyc_best = cur, cyc_bestP = P;
      if (k >= o.steps * 4 / 5 && k % 10 == 0) samples.push_back(cur);
      if ((k + 1) % every == 0 && o.progress && !o.progress(cyc, k + 1, best_e)) { stopped = true; break; }
    }
    P = cyc_bestP, cur = cyc_best;   // the cycle's lowest configuration, polished
    {   // polish downhill with steps shrinking from 0.3 to 0.005 Å (and 10° to 0.2°)
      const int np = std::max(400, o.steps / 5);
      for (int k = 0; k < np; ++k) {
        const double f = std::pow(0.005 / 0.3, double(k) / (np - 1));
        step(0, 0.3 * f, 10 * f * kPi / 180, 0, true);
      }
    }
    cur = total();
    if (cur < best_e) best_e = cur, bestP = P;
    AdsorptionConfig c;
    c.energy = cur, c.cycle = cyc + 1;
    c.positions.resize(n);
    for (size_t i = 0; i < n; ++i) c.positions[i] = s.atoms[i].pos;
    for (size_t m = 0; m < nm; ++m) for (size_t k = 0; k < M.mol[m].size(); ++k) c.positions[M.mol[m][k]] = P[m][k];
    finals.push_back(std::move(c));
  }
  if (stopped) rep.notes.push_back("stopped before the last cycle ended");
  // the lowest configuration into the structure
  P = bestP;
  for (size_t m = 0; m < nm; ++m) for (size_t k = 0; k < M.mol[m].size(); ++k) s.atoms[M.mol[m][k]].pos = P[m][k];
  s.velocities.clear();
  {
    AdsorptionConfig b;
    b.energy = best_e, b.cycle = 0;
    b.positions.resize(n);
    for (size_t i = 0; i < n; ++i) b.positions[i] = s.atoms[i].pos;
    finals.push_back(std::move(b));
  }
  std::sort(finals.begin(), finals.end(), [](const auto& a, const auto& b) { return a.energy < b.energy; });
  for (const auto& c : finals) {
    if (!rep.configs.empty() && std::fabs(c.energy - rep.configs.back().energy) < 1e-6) continue;   // the same minimum twice
    rep.configs.push_back(c);
    if (int(rep.configs.size()) >= std::max(1, o.keep)) break;
  }
  // energies of the lowest configuration, by kind and per component
  double as = 0, aa = 0;
  std::vector<double> per(nm, 0.0);
  for (size_t m = 0; m < nm; ++m) {
    const double f = M.with_fixed(m, P[m]);
    as += f, per[m] += f;
    for (size_t k = m + 1; k < nm; ++k) {
      const double e = M.between(m, P[m], k, P[k]);
      aa += e, per[m] += e, per[k] += e;
    }
  }
  rep.adsorption_energy = as + aa, rep.adsorbate_substrate = as, rep.adsorbate_adsorbate = aa;
  std::map<std::string, std::pair<int, double>> comp;
  std::vector<std::string> order;
  for (size_t m = 0; m < nm; ++m) {
    const std::string f = formula(s, M.mol[m]);
    if (!comp.count(f)) order.push_back(f);
    auto& c = comp[f];
    ++c.first, c.second += per[m];
  }
  for (const auto& f : order) rep.components.push_back({f, comp[f].first, comp[f].second / comp[f].first});
  if (!samples.empty()) {
    const double lo = *std::min_element(samples.begin(), samples.end()), hi = *std::max_element(samples.begin(), samples.end());
    const int nb = 30;
    const double w = hi > lo ? (hi - lo) / nb : 1.0;
    for (int b = 0; b <= nb; ++b) rep.hist_edges.push_back(lo + b * w);
    rep.hist_counts.assign(size_t(nb), 0.0);
    for (double x : samples) rep.hist_counts[size_t(std::clamp(int((x - lo) / w), 0, nb - 1))] += 1;
  }
  rep.steps = att;
  rep.acceptance = att ? double(acc) / double(att) : 0;
  {
    char b[240];
    std::snprintf(b, sizeof b, "%zu adsorbate molecule%s on %zu fixed atoms · %d cycle%s of %d steps, %.0f → %.0f K · lowest %.3f kcal/mol",
                  nm, nm == 1 ? "" : "s", M.fixed.size(), o.cycles, o.cycles == 1 ? "" : "s", o.steps, o.t_high, o.t_low, best_e);
    rep.notes.insert(rep.notes.begin(), b);
    rep.notes.push_back(std::string("energies: the force field's van der Waals") + (M.pm.coulomb() ? " and damped shifted force electrostatics" : "") +
                        " between molecules; rigid adsorbates (no deformation energy)");
  }
  rep.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  if (report) *report = std::move(rep);
}

}  // namespace caps
