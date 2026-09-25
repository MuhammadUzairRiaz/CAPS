// Primitive-path analysis (Everaers et al., Science 2004): see entangle.hpp.
#include "caps/entangle.hpp"

#include <algorithm>
#include <cmath>

namespace caps {
namespace {

// Minimum image by rounding fractional components (the cells CAPS builds are orthorhombic or mildly skewed).
struct Box {
  Cell cell;
  bool valid = false;
  Vec3 ra, rb, rc;   // rows of the inverse cell matrix
  Vec3 to_frac(const Vec3& d) const { return {dot(d, ra), dot(d, rb), dot(d, rc)}; }
  Vec3 image(const Vec3& d) const {
    if (!valid) return d;
    Vec3 f = to_frac(d);
    for (int k = 0; k < 3; ++k)
      if (cell.periodic[k]) f[k] -= std::round(f[k]);
    return cell.a * f[0] + cell.b * f[1] + cell.c * f[2];
  }
};

struct NeighbourList {
  std::vector<std::pair<uint32_t, uint32_t>> pairs;
  std::vector<Vec3> at_build;
};

// Pairs of beads on different chains closer than rc (a cell list over the wrapped positions).
void build_pairs(const Box& box, const std::vector<Vec3>& x, const std::vector<int>& chain, double rc, NeighbourList& nl) {
  nl.pairs.clear();
  nl.at_build = x;
  const size_t n = x.size();
  Cell c = box.cell;
  Vec3 ra = box.ra, rb = box.rb, rcv = box.rc;
  std::array<bool, 3> per = c.periodic;
  if (!box.valid) {
    Vec3 lo{1e300, 1e300, 1e300}, hi{-1e300, -1e300, -1e300};
    for (const auto& p : x)
      for (int k = 0; k < 3; ++k) lo[k] = std::min(lo[k], p[k]), hi[k] = std::max(hi[k], p[k]);
    c.origin = lo;
    c.a = {hi[0] - lo[0] + 1e-6, 0, 0};
    c.b = {0, hi[1] - lo[1] + 1e-6, 0};
    c.c = {0, 0, hi[2] - lo[2] + 1e-6};
    ra = {1 / c.a[0], 0, 0};
    rb = {0, 1 / c.b[1], 0};
    rcv = {0, 0, 1 / c.c[2]};
    per = {false, false, false};
  }
  const double v = std::abs(dot(c.a, cross(c.b, c.c)));
  const double w[3] = {v / norm(cross(c.b, c.c)), v / norm(cross(c.c, c.a)), v / norm(cross(c.a, c.b))};
  int nb[3];
  for (int k = 0; k < 3; ++k) {
    nb[k] = std::max(1, std::min(256, int(w[k] / rc)));
    if (nb[k] < 3) nb[k] = 1;   // fewer than three bins: one bin, every pair checked
  }
  std::vector<std::vector<uint32_t>> bins(size_t(nb[0]) * nb[1] * nb[2]);
  std::vector<std::array<int, 3>> where(n);
  for (size_t i = 0; i < n; ++i) {
    const Vec3 d = x[i] - c.origin;
    double f[3] = {dot(d, ra), dot(d, rb), dot(d, rcv)};
    for (int k = 0; k < 3; ++k) {
      if (per[k]) f[k] -= std::floor(f[k]);
      where[i][k] = std::clamp(int(f[k] * nb[k]), 0, nb[k] - 1);
    }
    bins[(size_t(where[i][0]) * nb[1] + where[i][1]) * nb[2] + where[i][2]].push_back(uint32_t(i));
  }
  const double rc2 = rc * rc;
  for (size_t i = 0; i < n; ++i) {
    for (int dx = -1; dx <= 1; ++dx) {
      if (nb[0] == 1 && dx) continue;
      int bx = where[i][0] + dx;
      if (bx < 0 || bx >= nb[0]) { if (!per[0]) continue; bx = (bx + nb[0]) % nb[0]; }
      for (int dy = -1; dy <= 1; ++dy) {
        if (nb[1] == 1 && dy) continue;
        int by = where[i][1] + dy;
        if (by < 0 || by >= nb[1]) { if (!per[1]) continue; by = (by + nb[1]) % nb[1]; }
        for (int dz = -1; dz <= 1; ++dz) {
          if (nb[2] == 1 && dz) continue;
          int bz = where[i][2] + dz;
          if (bz < 0 || bz >= nb[2]) { if (!per[2]) continue; bz = (bz + nb[2]) % nb[2]; }
          for (uint32_t j : bins[(size_t(bx) * nb[1] + by) * nb[2] + bz]) {
            if (j <= i || chain[j] == chain[i]) continue;
            const Vec3 d = box.image(x[j] - x[i]);
            if (dot(d, d) < rc2) nl.pairs.emplace_back(uint32_t(i), j);
          }
        }
      }
    }
  }
}

}  // namespace

PrimitivePaths primitive_paths(const System& frame, const std::vector<std::vector<uint32_t>>& backbones, const PrimitivePathOptions& o) {
  PrimitivePaths out;
  std::vector<Vec3> x;
  std::vector<int> chain;
  std::vector<char> fixed;
  std::vector<std::pair<uint32_t, uint32_t>> bonds;
  std::vector<size_t> start;
  double bl = 0, bn = 0;
  for (size_t c = 0; c < backbones.size(); ++c) {
    const auto& b = backbones[c];
    if (b.size() < 2) continue;
    start.push_back(x.size());
    for (size_t k = 0; k < b.size(); ++k) {
      if (k) {
        bonds.emplace_back(uint32_t(x.size() - 1), uint32_t(x.size()));
        bl += norm(frame.atoms[b[k]].pos - frame.atoms[b[k - 1]].pos);
        bn += 1;
      }
      x.push_back(frame.atoms[b[k]].pos);
      chain.push_back(int(start.size() - 1));
      fixed.push_back(k == 0 || k + 1 == b.size());
    }
  }
  start.push_back(x.size());
  const size_t nchains = start.size() - 1;
  if (!nchains) return out;
  const double sigma = o.sigma > 0 ? o.sigma : bl / bn / 0.97;
  out.sigma = sigma;

  Box box;
  box.cell = frame.cell;
  box.valid = frame.cell.valid();
  if (box.valid) {
    const Cell& c = frame.cell;
    const double v = dot(c.a, cross(c.b, c.c));
    box.ra = cross(c.b, c.c) * (1 / v);
    box.rb = cross(c.c, c.a) * (1 / v);
    box.rc = cross(c.a, c.b) * (1 / v);
  }

  // Everaers et al.'s potentials in reduced units of σ: FENE bonds (K = 30 ε/σ², R₀ = 1.5 σ, no rest length, so they pull
  // the path tight) and the WCA repulsion between beads of different chains (LJ cut at 2^{1/6} σ and shifted). σ is
  // set so the mean backbone bond is 0.97 σ, the Kremer–Grest bond.
  const double K = 30.0 / (sigma * sigma), R0 = 1.5 * sigma, eps = o.repulsion;
  const double rc_wca = std::pow(2.0, 1.0 / 6.0) * sigma;
  const size_t n = x.size();
  std::vector<Vec3> f(n), vel(n, Vec3{0, 0, 0});
  NeighbourList nl;
  const double skin = 0.3 * sigma, rlist = rc_wca + skin;
  build_pairs(box, x, chain, rlist, nl);
  const double k_bond = K;
  auto forces = [&] {
    std::fill(f.begin(), f.end(), Vec3{0, 0, 0});
    for (const auto& [i, j] : bonds) {
      const Vec3 d = x[j] - x[i];   // unwrapped along the chain
      const double r2 = std::min(dot(d, d), 0.98 * R0 * R0);
      const double g = K / (1 - r2 / (R0 * R0));   // −dU/dr / r
      f[i] = f[i] + d * g;
      f[j] = f[j] - d * g;
    }
    const double s2 = sigma * sigma;
    for (const auto& [i, j] : nl.pairs) {
      const Vec3 d = box.image(x[j] - x[i]);
      const double r2 = dot(d, d);
      if (r2 >= rc_wca * rc_wca) continue;
      const double ir2 = s2 / std::max(r2, 0.25 * s2), ir6 = ir2 * ir2 * ir2;
      const double g = 24 * eps * ir6 * (2 * ir6 - 1) / std::max(r2, 0.25 * s2);   // −dU/dr / r
      f[i] = f[i] - d * g;
      f[j] = f[j] + d * g;
    }
    double fmax = 0;
    for (size_t i = 0; i < n; ++i) {
      if (fixed[i]) { f[i] = {0, 0, 0}; continue; }
      fmax = std::max(fmax, norm(f[i]));
    }
    return fmax;
  };
  auto mean_lpp = [&] {
    double s = 0;
    for (const auto& [i, j] : bonds) s += norm(x[j] - x[i]);
    return s / double(nchains);
  };

  // FIRE (Bitzek et al., PRL 97, 170201, 2006), run per chain: each chain keeps its own time step and mixing, so one
  // chain that overshoots does not stop the others. Steps are capped at 0.05 σ so no bead jumps through another chain.
  // The run stops when the forces vanish or the mean L_pp changes by less than 10⁻⁴ of itself over 1000 steps.
  struct Fire { double dt, alpha = 0.1; int since_neg = 0; };
  const double dt0 = o.dt * sigma, dt_max = 10 * dt0, max_move = 0.05 * sigma, tol = o.force_tol * k_bond * sigma;
  std::vector<Fire> fire(o.per_chain ? nchains : 1, Fire{dt0});
  double fmax = forces();
  out.trace_step.push_back(0);
  out.trace_lpp.push_back(mean_lpp());
  double l_check = out.trace_lpp.back();
  int step = 0;
  for (; step < o.max_steps; ++step) {
    if (fmax < tol) { out.converged = true; break; }
    double moved = 0;
    for (size_t g = 0; g < fire.size(); ++g) {
      const size_t lo = o.per_chain ? start[g] : 0, hi = o.per_chain ? start[g + 1] : n;
      Fire& F = fire[g];
      double P = 0, vv = 0, ff = 0;
      for (size_t i = lo; i < hi; ++i) { P += dot(f[i], vel[i]); vv += dot(vel[i], vel[i]); ff += dot(f[i], f[i]); }
      const double vn = std::sqrt(vv), fn = std::sqrt(ff);
      if (P > 0) {
        for (size_t i = lo; i < hi; ++i) vel[i] = vel[i] * (1 - F.alpha) + f[i] * (F.alpha * vn / std::max(fn, 1e-300));
        if (++F.since_neg > 5) { F.dt = std::min(F.dt * 1.1, dt_max); F.alpha *= 0.99; }
      } else {
        for (size_t i = lo; i < hi; ++i) vel[i] = {0, 0, 0};
        F.dt = std::max(F.dt * 0.5, 0.1 * dt0);
        F.alpha = 0.1;
        F.since_neg = 0;
      }
      for (size_t i = lo; i < hi; ++i) {
        if (fixed[i]) continue;
        vel[i] = vel[i] + f[i] * F.dt;
        Vec3 dx = vel[i] * F.dt;
        const double m = norm(dx);
        if (m > max_move) dx = dx * (max_move / m);
        x[i] = x[i] + dx;
        moved = std::max(moved, norm(x[i] - nl.at_build[i]));
      }
    }
    if (moved > 0.5 * skin) build_pairs(box, x, chain, rlist, nl);
    fmax = forces();
    if ((step + 1) % o.record_every == 0) {
      out.trace_step.push_back(step + 1);
      out.trace_lpp.push_back(mean_lpp());
      if ((step + 1) % 1000 == 0) {
        if (std::abs(out.trace_lpp.back() - l_check) < 1e-4 * out.trace_lpp.back()) { out.converged = true; ++step; break; }
        l_check = out.trace_lpp.back();
      }
      if (o.progress && !o.progress(double(step + 1) / o.max_steps)) { out.stopped = true; ++step; break; }
    }
  }
  out.steps = step;
  if (out.trace_step.back() != step) { out.trace_step.push_back(step); out.trace_lpp.push_back(mean_lpp()); }

  for (size_t c = 0; c < nchains; ++c) {
    std::vector<Vec3> p(x.begin() + long(start[c]), x.begin() + long(start[c + 1]));
    double L = 0;
    for (size_t k = 1; k < p.size(); ++k) L += norm(p[k] - p[k - 1]);
    const Vec3 R = p.back() - p.front();
    out.lpp.push_back(L);
    out.r2.push_back(dot(R, R));
    out.bonds.push_back(double(p.size() - 1));
    out.paths.push_back(std::move(p));
  }
  return out;
}

EntanglementEstimate entanglement_estimate(const PrimitivePaths& p) {
  EntanglementEstimate e;
  e.chains = int(p.lpp.size());
  if (!e.chains) return e;
  for (size_t c = 0; c < p.lpp.size(); ++c) {
    e.nb += p.bonds[c];
    e.r2 += p.r2[c];
    e.lpp += p.lpp[c];
    e.lpp2 += p.lpp[c] * p.lpp[c];
  }
  e.nb /= e.chains;
  e.r2 /= e.chains;
  e.lpp /= e.chains;
  e.lpp2 /= e.chains;
  if (e.lpp > 0) e.a_pp = e.r2 / e.lpp;
  // straight paths (L_pp² ≈ R²) mean no entanglements: N_e is left at 0 then
  if (e.lpp * e.lpp > 1.0001 * e.r2) e.ne_coil = e.nb * e.r2 / (e.lpp * e.lpp);
  if (e.lpp2 > 1.0001 * e.r2) {
    e.ne_mscoil = e.nb / (e.lpp2 / e.r2 - 1);
    e.z = e.nb / e.ne_mscoil;
  }
  return e;
}

}  // namespace caps
