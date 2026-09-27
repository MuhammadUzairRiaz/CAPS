// CAPS CBMC: configurational-bias regrowth of chain ends (see cbmc.hpp).
#include "caps/cbmc.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <limits>
#include <map>
#include <random>
#include <stdexcept>
#include <unordered_set>

#include "caps/analysis.hpp"
#include "caps/typing.hpp"

namespace caps {

namespace {

constexpr double kCoulomb = 332.06371;   // kcal·Å/(mol·e²), as the Evaluator
constexpr double kB = 0.0019872043;      // kcal/(mol·K)
constexpr double kPi = 3.14159265358979323846;

// p rotated by angle t about the axis through c along unit u (Rodrigues)
Vec3 rotate_about(const Vec3& p, const Vec3& c, const Vec3& u, double t) {
  const Vec3 v = p - c;
  const double ct = std::cos(t), st = std::sin(t);
  return c + v * ct + cross(u, v) * st + u * (dot(u, v) * (1 - ct));
}

// ln Σ exp(l) over finite entries (−∞ when none)
double log_sum_exp(const std::vector<double>& l) {
  double m = -std::numeric_limits<double>::infinity();
  for (double v : l) m = std::max(m, v);
  if (!std::isfinite(m)) return -std::numeric_limits<double>::infinity();
  double s = 0;
  for (double v : l) s += std::exp(v - m);
  return m + std::log(s);
}

struct EndBond { uint32_t near, far; std::vector<uint32_t> down; };   // rotating about near→far moves `down` (far excluded)
struct End { std::vector<EndBond> bonds; };                           // outermost bond first

}  // namespace

void cbmc_regrow(System& s, const ForceField& ff, const CbmcOptions& o, CbmcReport* report) {
  const auto t0 = std::chrono::steady_clock::now();
  CbmcReport rep;
  const size_t n = s.atoms.size();
  if (ff.type_index.size() != n || ff.charge.size() != n) throw std::invalid_argument("CBMC: the force field was assigned to another structure");
  if (!ff.pair_func.empty() || ff.hbond.on() || !ff.vsites.empty() || !ff.lj_pairs.empty() || ff.sw.on || !ff.excluded_type_pairs.empty())
    throw std::invalid_argument("CBMC: the trial energies need Lennard-Jones pairs only (no other pair forms, hydrogen bonds, virtual sites or explicit pairs)");
  if (o.trials < 2 || o.max_torsions < 1 || o.moves < 0 || o.temperature <= 0) throw std::invalid_argument("CBMC: at least 2 trials, 1 torsion, T > 0");
  make_molecules_whole(s);
  const Cell& cell = s.cell;
  const bool per = cell.valid();
  auto mi = [&](const Vec3& d) { return per ? cell.minimum_image(d) : d; };

  // the cut-off: at most half the narrowest periodic width (one image per pair)
  double rc = o.cutoff;
  double width[3] = {1e30, 1e30, 1e30};
  if (per) {
    const double V = cell.volume();
    width[0] = V / norm(cross(cell.b, cell.c));
    width[1] = V / norm(cross(cell.c, cell.a));
    width[2] = V / norm(cross(cell.a, cell.b));
    double wmin = 1e30;
    for (int d = 0; d < 3; ++d)
      if (cell.periodic[size_t(d)]) wmin = std::min(wmin, width[d]);
    if (rc > 0.5 * wmin - 1e-6) {
      rc = 0.5 * wmin - 1e-6;
      char b[160];
      std::snprintf(b, sizeof b, "cut-off reduced to %.2f Å (half the cell's narrowest width)", rc);
      rep.notes.push_back(b);
    }
  }
  if (rc < 3) throw std::invalid_argument("CBMC: the cell is too small for a 3 Å cut-off");
  rep.cutoff = rc;
  const double rc2 = rc * rc;

  // ---- pair energies: LJ (12-6 or class II 9-6) shifted to zero at rc, 1-4 scaled; DSF Coulomb
  const size_t nt = ff.lj.size();
  const bool lj96 = ff.pair_form == "lj9-6";
  std::vector<double> eps(nt * nt), sig(nt * nt), eps14, sig14;
  auto lj_raw = [&](double r2, double e, double sg) {
    if (lj96) {
      const double t = sg / std::sqrt(r2), t3 = t * t * t;
      return e * (2 * t3 * t3 * t3 - 3 * t3 * t3);
    }
    const double q = sg * sg / r2, q3 = q * q * q;
    return 4 * e * (q3 * q3 - q3);
  };
  std::vector<double> shift(nt * nt), shift14;
  for (size_t a = 0; a < nt; ++a)
    for (size_t b = 0; b < nt; ++b) {
      const PairType p = mixed_pair(ff, int(a), int(b));
      eps[a * nt + b] = p.eps, sig[a * nt + b] = p.sigma;
      shift[a * nt + b] = lj_raw(rc2, p.eps, p.sigma);
    }
  if (!ff.lj14_types.empty()) {   // separate 1-4 parameters (CHARMM, GROMOS), mixed by the force field's rule
    eps14.resize(nt * nt), sig14.resize(nt * nt), shift14.resize(nt * nt);
    for (size_t a = 0; a < nt; ++a)
      for (size_t b = 0; b < nt; ++b) {
        const auto& A = ff.lj14_types[a];
        const auto& B = ff.lj14_types[b];
        const double sg = ff.mixing == "geometric" ? std::sqrt(A.sigma * B.sigma) : 0.5 * (A.sigma + B.sigma);
        eps14[a * nt + b] = std::sqrt(A.eps * B.eps), sig14[a * nt + b] = sg;
        shift14[a * nt + b] = lj_raw(rc2, eps14[a * nt + b], sg);
      }
  }
  const double alpha = o.dsf_alpha, erfc_c = std::erfc(alpha * rc) / rc;
  const double dsf_f = erfc_c / rc + 2 * alpha / std::sqrt(kPi) * std::exp(-alpha * alpha * rc2) / rc;
  const bool coul = o.coulomb && std::any_of(ff.charge.begin(), ff.charge.end(), [](double q) { return q != 0; });
  const double qscale = kCoulomb / (ff.dielectric > 0 ? ff.dielectric : 1.0);

  // exclusions: the force field's (1-2, 1-3 and 1-4 partners) with its 1-4 pairs, else from the bonds
  const auto nb = s.neighbours();
  std::vector<std::vector<uint32_t>> excl(n);
  std::unordered_set<uint64_t> is14;
  auto key = [&](uint32_t a, uint32_t b) { return a < b ? uint64_t(a) * n + b : uint64_t(b) * n + a; };
  if (ff.excluded.size() == n) {
    excl = ff.excluded;
    for (const auto& p : ff.pairs14) is14.insert(key(p[0], p[1]));
  } else {
    for (uint32_t a = 0; a < n; ++a) {
      std::map<uint32_t, int> d{{a, 0}};
      std::vector<uint32_t> fr{a};
      for (int depth = 1; depth <= 3; ++depth) {
        std::vector<uint32_t> nx;
        for (uint32_t u : fr)
          for (uint32_t w : nb[u])
            if (!d.count(w)) d[w] = depth, nx.push_back(w);
        fr = std::move(nx);
      }
      for (const auto& [w, dd] : d)
        if (dd > 0) {
          excl[a].push_back(w);
          if (dd == 3) is14.insert(key(a, w));
        }
      std::sort(excl[a].begin(), excl[a].end());
    }
  }
  auto pair_energy = [&](uint32_t a, uint32_t b, const Vec3& xa, const Vec3& xb) {
    const Vec3 d = mi(xb - xa);
    const double r2 = dot(d, d);
    if (r2 >= rc2) return 0.0;
    double fl = 1, fq = 1;
    bool p14 = false;
    if (std::binary_search(excl[a].begin(), excl[a].end(), b)) {
      if (!is14.count(key(a, b))) return 0.0;
      fl = ff.lj14, fq = ff.coul14, p14 = true;
    }
    if (r2 < 1e-6) return std::numeric_limits<double>::infinity();
    const size_t tp = size_t(ff.type_index[a]) * nt + size_t(ff.type_index[b]);
    double e = p14 && !eps14.empty() ? fl * (lj_raw(r2, eps14[tp], sig14[tp]) - shift14[tp]) : fl * (lj_raw(r2, eps[tp], sig[tp]) - shift[tp]);
    if (coul && fq != 0) {
      const double r = std::sqrt(r2), qq = ff.charge[a] * ff.charge[b];
      if (qq != 0) e += fq * qscale * qq * (std::erfc(alpha * r) / r - erfc_c + dsf_f * (r - rc));
    }
    return e;
  };

  // ---- torsion terms by central bond
  std::map<std::pair<uint32_t, uint32_t>, TorsionRefs> tors;
  auto bk = [](uint32_t a, uint32_t b) { return std::make_pair(std::min(a, b), std::max(a, b)); };
  for (size_t k = 0; k < ff.dihedrals.size(); ++k) tors[bk(ff.dihedrals[k].j, ff.dihedrals[k].k)].dihedrals.push_back(uint32_t(k));
  for (size_t k = 0; k < ff.dihedrals2.size(); ++k) tors[bk(ff.dihedrals2[k].j, ff.dihedrals2[k].k)].dihedrals2.push_back(uint32_t(k));
  for (size_t k = 0; k < ff.cbt.size(); ++k) tors[bk(ff.cbt[k].j, ff.cbt[k].k)].cbt.push_back(uint32_t(k));
  const TorsionRefs none;
  auto tors_of = [&](uint32_t a, uint32_t b) -> const TorsionRefs& {
    const auto it = tors.find(bk(a, b));
    return it == tors.end() ? none : it->second;
  };

  // ---- rotatable backbone bonds at each chain end
  const Perception P = perceive(s);
  auto el = [&](uint32_t a) { return s.atoms[a].element; };
  auto carbonyl = [&](uint32_t c) {
    if (el(c) != 6) return false;
    for (uint32_t w : P.nb[c])
      if ((el(w) == 8 || el(w) == 16) && P.bond_order(c, w) == 2) return true;
    return false;
  };
  auto rotatable = [&](uint32_t a, uint32_t b) {
    if (P.bond_order(a, b) != 1) return false;
    for (size_t k = 0; k < P.nb[a].size(); ++k)
      if (P.nb[a][k] == b && (P.ring_bond[a][k] || P.arom_bond[a][k])) return false;
    if (nb[a].size() < 2 || nb[b].size() < 2) return false;
    auto hetero = [&](uint32_t x) { return el(x) == 7 || el(x) == 8 || el(x) == 16; };
    if ((carbonyl(a) && hetero(b)) || (carbonyl(b) && hetero(a))) return false;   // amides, esters, thioesters
    return true;
  };
  // the atoms on far's side of the bond near–far (far itself left out); empty when the bond is in a ring or the side is too big
  auto downstream = [&](uint32_t near, uint32_t far) {
    std::vector<uint32_t> out;
    std::vector<uint32_t> st{far};
    std::unordered_set<uint32_t> seen{far, near};
    while (!st.empty()) {
      const uint32_t u = st.back();
      st.pop_back();
      for (uint32_t w : nb[u]) {
        if (u == far && w == near) continue;
        if (w == near) return std::vector<uint32_t>{};   // a ring through the bond
        if (seen.insert(w).second) {
          out.push_back(w);
          st.push_back(w);
          if (int(out.size()) > o.max_atoms) return std::vector<uint32_t>{};
        }
      }
    }
    return out;
  };
  const auto bbs = backbones(s, 4);
  std::vector<End> ends;
  auto heavy_deg = [&](uint32_t a) { int k = 0; for (uint32_t w : nb[a]) k += el(w) != 1; return k; };
  for (const auto& path : bbs) {
    for (int side = 0; side < 2; ++side) {
      std::vector<uint32_t> q(path.begin(), path.end());
      if (side == 0) std::reverse(q.begin(), q.end());   // the end last
      if (heavy_deg(q.back()) != 1) continue;             // only a molecule's own end
      End e;
      for (int i = int(q.size()) - 2; i >= 0 && int(e.bonds.size()) < o.max_torsions; --i) {
        if (!rotatable(q[size_t(i)], q[size_t(i) + 1])) continue;
        auto d = downstream(q[size_t(i)], q[size_t(i) + 1]);
        if (d.empty()) break;   // a ring, or the side is too big: the bonds further in are not an end
        e.bonds.push_back({q[size_t(i)], q[size_t(i) + 1], std::move(d)});
      }
      if (!e.bonds.empty()) ends.push_back(std::move(e));
    }
  }
  rep.chains = int(ends.size());
  if (ends.empty()) throw std::invalid_argument("CBMC: no chain end with a rotatable backbone bond (chains of at least four heavy atoms, with a free end)");

  std::vector<Vec3> x(n);
  for (size_t i = 0; i < n; ++i) x[i] = s.atoms[i].pos;
  auto r2_mean = [&] {
    double t = 0;
    for (const auto& p : bbs) { const Vec3 d = x[p.back()] - x[p.front()]; t += dot(d, d); }
    return bbs.empty() ? 0.0 : t / double(bbs.size());
  };
  rep.r2_before = r2_mean();

  // ---- cell list over fractional coordinates (bins at least rc wide; fewer than three: every bin once)
  int nbin[3] = {1, 1, 1};
  if (per)
    for (int d = 0; d < 3; ++d)
      if (cell.periodic[size_t(d)]) nbin[d] = std::max(1, int(std::floor(width[d] / rc)));
  auto bin_of = [&](const Vec3& p) {
    if (!per) return 0;
    const Vec3 f = cell.to_fractional(p);
    int id[3];
    for (int d = 0; d < 3; ++d) {
      const double g = f[size_t(d)] - std::floor(f[size_t(d)]);
      id[d] = std::min(nbin[d] - 1, int(g * nbin[d]));
    }
    return (id[0] * nbin[1] + id[1]) * nbin[2] + id[2];
  };
  std::vector<std::vector<uint32_t>> bins(size_t(nbin[0] * nbin[1] * nbin[2]));
  std::vector<int> abin(n);
  for (uint32_t i = 0; i < n; ++i) bins[size_t(abin[i] = bin_of(x[i]))].push_back(i);
  std::vector<std::vector<int>> nbins(bins.size());   // each bin's neighbour bins, each once
  for (int a = 0; a < nbin[0]; ++a)
    for (int b = 0; b < nbin[1]; ++b)
      for (int c = 0; c < nbin[2]; ++c) {
        std::vector<int> list;
        auto offs = [&](int d, int v) {
          std::vector<int> r;
          if (nbin[d] < 3) for (int k = 0; k < nbin[d]; ++k) r.push_back(k);
          else for (int dv = -1; dv <= 1; ++dv) r.push_back(((v + dv) % nbin[d] + nbin[d]) % nbin[d]);
          return r;
        };
        for (int i : offs(0, a))
          for (int j : offs(1, b))
            for (int k : offs(2, c)) list.push_back((i * nbin[1] + j) * nbin[2] + k);
        nbins[size_t((a * nbin[1] + b) * nbin[2] + c)] = std::move(list);
      }

  // ---- moves
  std::mt19937_64 rng(o.seed);
  std::uniform_real_distribution<double> U(0.0, 1.0);
  const double beta = 1 / (kB * o.temperature);
  const int k = o.trials;
  std::vector<int> mark(n, -1);          // growth step of each moved atom this move; −1 fixed
  std::vector<uint32_t> moved;           // atoms moved (the innermost bond's side), for restoring and marks
  std::vector<Vec3> saved;
  const int every = std::max(1, o.moves / 100);
  for (int mv = 0; mv < o.moves; ++mv) {
    const End& e = ends[size_t(U(rng) * double(ends.size())) % ends.size()];
    const int m = 1 + int(U(rng) * double(e.bonds.size())) % int(e.bonds.size());
    // growth order: innermost of the m bonds first
    std::vector<const EndBond*> G;
    for (int g = m - 1; g >= 0; --g) G.push_back(&e.bonds[size_t(g)]);
    moved = G[0]->down;
    // S_g = down(G[g]) \ down(G[g+1]): step g places them
    for (int g = 0; g < m; ++g) {
      std::unordered_set<uint32_t> next;
      if (g + 1 < m) next.insert(G[size_t(g) + 1]->down.begin(), G[size_t(g) + 1]->down.end());
      for (uint32_t a : G[size_t(g)]->down)
        if (!next.count(a)) mark[a] = g;
    }
    std::vector<std::vector<uint32_t>> S(static_cast<size_t>(m));
    for (uint32_t a : moved) S[size_t(mark[a])].push_back(a);
    saved.resize(moved.size());
    for (size_t i = 0; i < moved.size(); ++i) saved[i] = x[moved[i]];

    // the energy of step g's atoms at positions y, with everything placed before them
    std::vector<Vec3> y;
    auto step_energy = [&](int g) {
      const auto& Sg = S[size_t(g)];
      double E = 0;
      for (size_t ia = 0; ia < Sg.size(); ++ia) {
        const uint32_t a = Sg[ia];
        const Vec3& xa = y[ia];
        for (int bn : nbins[size_t(bin_of(xa))])
          for (uint32_t b : bins[size_t(bn)]) {
            if (mark[b] >= 0) continue;   // moved atoms below
            E += pair_energy(a, b, xa, x[b]);
          }
        for (uint32_t b : moved)
          if (mark[b] >= 0 && mark[b] < g) E += pair_energy(a, b, xa, x[b]);
        if (!std::isfinite(E)) return E;
      }
      const EndBond& B = *G[size_t(g)];
      const auto& refs = tors_of(B.near, B.far);
      if (!refs.dihedrals.empty() || !refs.dihedrals2.empty() || !refs.cbt.empty()) {
        E += torsion_energy(ff, refs, [&](uint32_t a) -> Vec3 {
          if (mark[a] == g)
            for (size_t ia = 0; ia < Sg.size(); ++ia)
              if (Sg[ia] == a) return y[ia];
          return x[a];
        }, cell);
      }
      return E;
    };
    auto trial_positions = [&](int g, double t) {
      const EndBond& B = *G[size_t(g)];
      const Vec3 axis = x[B.far] - x[B.near];
      const Vec3 u = axis * (1 / norm(axis));
      const auto& Sg = S[size_t(g)];
      y.resize(Sg.size());
      for (size_t ia = 0; ia < Sg.size(); ++ia) y[ia] = t == 0 ? x[Sg[ia]] : rotate_about(x[Sg[ia]], x[B.far], u, t);
    };
    // old: retrace with the present torsion as one of the trials
    double lnw_old = 0, e_old = 0;
    std::vector<double> lw(static_cast<size_t>(k)), ue(static_cast<size_t>(k));
    for (int g = 0; g < m; ++g) {
      for (int t = 0; t < k; ++t) {
        trial_positions(g, t == 0 ? 0.0 : 2 * kPi * U(rng));
        ue[size_t(t)] = step_energy(g);
        lw[size_t(t)] = std::isfinite(ue[size_t(t)]) ? -beta * ue[size_t(t)] : -std::numeric_limits<double>::infinity();
      }
      lnw_old += log_sum_exp(lw);
      e_old += ue[0];
    }
    // new: grow, drawing each torsion by its weight
    double lnw_new = 0, e_new = 0;
    bool dead = false;
    for (int g = 0; g < m && !dead; ++g) {
      std::vector<double> angle(static_cast<size_t>(k));
      for (int t = 0; t < k; ++t) {
        angle[size_t(t)] = 2 * kPi * U(rng);
        trial_positions(g, angle[size_t(t)]);
        ue[size_t(t)] = step_energy(g);
        lw[size_t(t)] = std::isfinite(ue[size_t(t)]) ? -beta * ue[size_t(t)] : -std::numeric_limits<double>::infinity();
      }
      const double lse = log_sum_exp(lw);
      if (!std::isfinite(lse)) { dead = true; break; }
      lnw_new += lse;
      double r = U(rng), acc = 0;
      int pick = k - 1;
      for (int t = 0; t < k; ++t) {
        acc += std::exp(lw[size_t(t)] - lse);
        if (r < acc) { pick = t; break; }
      }
      while (!std::isfinite(lw[size_t(pick)])) pick = (pick + k - 1) % k;
      e_new += ue[size_t(pick)];
      // turn this bond's whole side (the atoms placed later ride along)
      const EndBond& B = *G[size_t(g)];
      const Vec3 axis = x[B.far] - x[B.near];
      const Vec3 u = axis * (1 / norm(axis));
      for (uint32_t a : B.down) x[a] = rotate_about(x[a], x[B.far], u, angle[size_t(pick)]);
    }
    ++rep.attempted;
    const bool accept = !dead && std::isfinite(lnw_old) ? std::log(std::max(U(rng), 1e-300)) < lnw_new - lnw_old : !dead;
    if (accept) {
      ++rep.accepted;
      rep.energy_change += e_new - e_old;
      for (uint32_t a : moved) {   // re-bin the moved atoms
        const int nbn = bin_of(x[a]);
        if (nbn != abin[a]) {
          auto& v = bins[size_t(abin[a])];
          v.erase(std::find(v.begin(), v.end(), a));
          bins[size_t(nbn)].push_back(a);
          abin[a] = nbn;
        }
      }
    } else {
      for (size_t i = 0; i < moved.size(); ++i) x[moved[i]] = saved[i];
    }
    for (uint32_t a : moved) mark[a] = -1;
    if ((mv + 1) % every == 0 && o.snapshot) o.snapshot(x);
    if (o.progress && (mv + 1) % every == 0 && !o.progress(mv + 1, rep.accepted)) {
      rep.notes.push_back("stopped after " + std::to_string(mv + 1) + " moves");
      break;
    }
  }
  for (size_t i = 0; i < n; ++i) s.atoms[i].pos = x[i];
  s.velocities.clear();
  rep.r2_after = r2_mean();
  rep.acceptance = rep.attempted ? double(rep.accepted) / rep.attempted : 0;
  {
    char b[240];
    std::snprintf(b, sizeof b, "%d of %d regrowths accepted (%.1f %%) · %d chain ends · %d trials per torsion, up to %d torsions · %.0f K · ⟨R²⟩ %.1f → %.1f Å²",
                  rep.accepted, rep.attempted, 100 * rep.acceptance, rep.chains, k, o.max_torsions, o.temperature, rep.r2_before, rep.r2_after);
    rep.notes.insert(rep.notes.begin(), b);
    rep.notes.push_back(std::string("trial energies: the force field's van der Waals") + (coul ? ", damped shifted force electrostatics" : "") + " and torsions about each regrown bond");
  }
  rep.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  if (report) *report = std::move(rep);
}

}  // namespace caps
