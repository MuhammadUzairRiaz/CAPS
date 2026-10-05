// CAPS DPD (see dpd.hpp): Groot–Warren dissipative particle dynamics in a cubic periodic box, reduced units.
#include "caps/dpd.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <numeric>
#include <array>
#include <complex>
#include <cstdio>
#include <random>
#include <set>
#include <stdexcept>

namespace caps {

namespace {

constexpr double kPi = 3.14159265358979323846;

std::string pair_key(char i, char j) { return i <= j ? std::string{i, j} : std::string{j, i}; }

int element_of_type(char t) {
  switch (t) { case 'A': return 6; case 'B': return 8; case 'C': return 7; case 'D': return 16; default: return 15; }
}

}  // namespace

double dpd_repulsion(const DpdOptions& o, char i, char j) {
  auto given = [&](char x, char y, double& v) {
    for (const auto& k : {std::string{x, y}, std::string{y, x}}) {
      const auto it = o.a.find(k);
      if (it != o.a.end()) { v = it->second; return true; }
    }
    return false;
  };
  double v = 0;
  if (given(i, j, v)) return v;
  auto like = [&](char x) { double w; return given(x, x, w) ? w : 75.0 / o.density; };
  const double base = 0.5 * (like(i) + like(j));
  if (i == j) return like(i);
  for (const auto& k : {std::string{i, j}, std::string{j, i}}) {
    const auto it = o.chi.find(k);
    if (it != o.chi.end()) return base + it->second / 0.286;   // Groot & Warren at ρ = 3
  }
  return base;
}

DpdReport run_dpd(const DpdOptions& o) {
  const auto t0 = std::chrono::steady_clock::now();
  DpdReport rep;
  if (o.density <= 0 || o.dt <= 0 || o.gamma <= 0) throw std::invalid_argument("DPD: density, dt and γ must be positive");
  // beads, their types and bonds, molecule by molecule
  std::vector<char> type;
  std::vector<int> mol;
  std::vector<std::pair<int, int>> bonds;
  std::vector<std::array<int, 3>> angles;
  std::set<char> used;
  for (const auto& sp : o.species) {
    if (sp.sequence.empty() || sp.count <= 0) continue;
    for (char c : sp.sequence) if (c < 'A' || c > 'Z') throw std::invalid_argument("DPD: bead types are the letters A–Z ('" + sp.sequence + "')");
    for (int k = 0; k < sp.count; ++k) {
      const int first = int(type.size());
      for (size_t b = 0; b < sp.sequence.size(); ++b) {
        type.push_back(sp.sequence[b]);
        mol.push_back(rep.molecules);
        used.insert(sp.sequence[b]);
        if (b > 0) bonds.push_back({first + int(b) - 1, first + int(b)});
        if (b > 1 && o.angle_k != 0) angles.push_back({first + int(b) - 2, first + int(b) - 1, first + int(b)});
      }
      ++rep.molecules;
    }
  }
  const int n = int(type.size());
  if (n < 2) throw std::invalid_argument("DPD: no beads (give species with sequences and counts)");
  rep.beads = n;
  rep.types.assign(used.begin(), used.end());
  if (!o.chi.empty() && std::fabs(o.density - 3.0) > 1e-9) rep.notes.push_back("χ → a_ij uses Groot & Warren's mapping for ρ = 3; at another density give a_ij");
  const double L = std::cbrt(n / o.density);
  rep.box = L;
  if (L < 3) throw std::invalid_argument("DPD: the box is under 3 r_c: add beads");
  // repulsion table
  std::map<std::pair<char, char>, double> A;
  for (char x : used) for (char y : used) A[{x, y}] = dpd_repulsion(o, x, y);
  std::vector<double> aij(26 * 26, 0.0);
  for (const auto& [k, v] : A) aij[size_t(k.first - 'A') * 26 + size_t(k.second - 'A')] = v;
  // start: chains as random walks of 0.7 r_c steps from random places; Maxwell velocities at kT = 1, no net momentum
  std::mt19937_64 rng(o.seed);
  std::uniform_real_distribution<double> U(0.0, 1.0);
  std::normal_distribution<double> G(0.0, 1.0);
  const size_t nn = static_cast<size_t>(n);
  std::vector<Vec3> r(nn), v(nn), f(nn), fold(nn);
  for (int i = 0; i < n; ++i) {
    if (i > 0 && mol[size_t(i)] == mol[size_t(i - 1)]) {
      const double z = 2 * U(rng) - 1, ph = 2 * kPi * U(rng), s = std::sqrt(1 - z * z);
      r[size_t(i)] = r[size_t(i - 1)] + Vec3{s * std::cos(ph), s * std::sin(ph), z} * 0.7;
    } else {
      r[size_t(i)] = {L * U(rng), L * U(rng), L * U(rng)};
    }
    v[size_t(i)] = {G(rng), G(rng), G(rng)};
  }
  Vec3 p{0, 0, 0};
  for (const auto& x : v) p = p + x;
  for (auto& x : v) x = x - p * (1.0 / n);
  auto wrapd = [&](Vec3 d) { for (auto& x : d) x -= L * std::round(x / L); return d; };
  // cell list with cells ≥ r_c
  const int nc = std::max(3, int(std::floor(L)));
  const double cw = L / nc;
  std::vector<int> head(static_cast<size_t>(nc * nc * nc)), next(nn);
  const double sigma = std::sqrt(2 * o.gamma * 1.0), sdt = 1 / std::sqrt(o.dt);
  const double sq3 = std::sqrt(3.0);
  double virial = 0;
  auto forces = [&](const std::vector<Vec3>& vel) {
    std::fill(f.begin(), f.end(), Vec3{0, 0, 0});
    std::fill(head.begin(), head.end(), -1);
    for (int i = 0; i < n; ++i) {
      int c[3];
      for (int k = 0; k < 3; ++k) {
        double x = r[size_t(i)][size_t(k)] - L * std::floor(r[size_t(i)][size_t(k)] / L);
        c[k] = std::min(nc - 1, int(x / cw));
      }
      const int id = (c[0] * nc + c[1]) * nc + c[2];
      next[size_t(i)] = head[size_t(id)], head[size_t(id)] = i;
    }
    virial = 0;
    for (int cx = 0; cx < nc; ++cx)
      for (int cy = 0; cy < nc; ++cy)
        for (int cz = 0; cz < nc; ++cz) {
          const int id = (cx * nc + cy) * nc + cz;
          for (int dx = -1; dx <= 1; ++dx)
            for (int dy = -1; dy <= 1; ++dy)
              for (int dz = -1; dz <= 1; ++dz) {
                const int jd = (((cx + dx + nc) % nc) * nc + (cy + dy + nc) % nc) * nc + (cz + dz + nc) % nc;
                if (jd < id) continue;   // each cell pair once
                for (int i = head[size_t(id)]; i >= 0; i = next[size_t(i)])
                  for (int j = (jd == id ? next[size_t(i)] : head[size_t(jd)]); j >= 0; j = next[size_t(j)]) {
                    const Vec3 d = wrapd(r[size_t(i)] - r[size_t(j)]);
                    const double r2 = dot(d, d);
                    if (r2 >= 1.0 || r2 < 1e-12) continue;
                    const double rr = std::sqrt(r2), w = 1 - rr;
                    const Vec3 e = d * (1 / rr);
                    const double a = aij[size_t(type[size_t(i)] - 'A') * 26 + size_t(type[size_t(j)] - 'A')];
                    const double theta = sq3 * (2 * U(rng) - 1);   // unit variance
                    const double fc = a * w, fd = -o.gamma * w * w * dot(e, vel[size_t(i)] - vel[size_t(j)]), fr = sigma * w * theta * sdt;
                    const Vec3 F = e * (fc + fd + fr);
                    f[size_t(i)] = f[size_t(i)] + F;
                    f[size_t(j)] = f[size_t(j)] - F;
                    virial += fc * rr;
                  }
              }
        }
    for (const auto& [i, j] : bonds) {
      const Vec3 d = wrapd(r[size_t(i)] - r[size_t(j)]);
      const double rr = norm(d);
      if (rr < 1e-12) continue;
      const double fb = -o.bond_k * (rr - o.bond_r0);
      const Vec3 F = d * (fb / rr);
      f[size_t(i)] = f[size_t(i)] + F;
      f[size_t(j)] = f[size_t(j)] - F;
      virial += fb * rr;
    }
    // stiffness: E = k (1 + cos θ), θ at the middle bead (LAMMPS angle_style cosine); F_i = −k ∂cos θ/∂r_i (the middle
    // bead takes the rest)
    for (const auto& [i, j, k] : angles) {
      const Vec3 a = wrapd(r[size_t(i)] - r[size_t(j)]), b = wrapd(r[size_t(k)] - r[size_t(j)]);
      const double la = norm(a), lb = norm(b);
      if (la < 1e-12 || lb < 1e-12) continue;
      const double c = dot(a, b) / (la * lb);
      const Vec3 Fi = (b * (1 / (la * lb)) - a * (c / (la * la))) * -o.angle_k;
      const Vec3 Fk = (a * (1 / (la * lb)) - b * (c / (lb * lb))) * -o.angle_k;
      f[size_t(i)] = f[size_t(i)] + Fi;
      f[size_t(k)] = f[size_t(k)] + Fk;
      f[size_t(j)] = f[size_t(j)] - Fi - Fk;
      virial += dot(Fi, a) + dot(Fk, b);
    }
  };
  // segregation of the first two types over cells of about 2 r_c
  const char ta = rep.types.empty() ? 'A' : rep.types[0], tb = rep.types.size() > 1 ? rep.types[1] : ta;
  auto order_now = [&] {
    if (ta == tb) return 0.0;
    const int g = std::max(1, int(std::floor(L / 2)));
    std::vector<int> na(size_t(g * g * g), 0), nb(size_t(g * g * g), 0);
    for (int i = 0; i < n; ++i) {
      int c[3];
      for (int k = 0; k < 3; ++k) c[k] = std::min(g - 1, int((r[size_t(i)][size_t(k)] - L * std::floor(r[size_t(i)][size_t(k)] / L)) / L * g));
      const int id = (c[0] * g + c[1]) * g + c[2];
      if (type[size_t(i)] == ta) ++na[size_t(id)];
      else if (type[size_t(i)] == tb) ++nb[size_t(id)];
    }
    double num = 0, den = 0;
    for (size_t k = 0; k < na.size(); ++k) num += std::abs(na[k] - nb[k]), den += na[k] + nb[k];
    return den > 0 ? num / den : 0.0;
  };
  // integrate (Groot & Warren's modified velocity Verlet, λ = 0.65)
  const double lambda = 0.65, dt = o.dt, V = L * L * L;
  forces(v);
  std::vector<double> kts, ps;
  std::vector<Vec3> vt(nn);
  const long every = std::max(1L, o.steps / 100);
  auto store_frame = [&](long step) {
    std::vector<Vec3> x(nn);
    for (int i = 0; i < n; ++i) x[size_t(i)] = r[size_t(i)] * o.rc_angstrom;
    rep.frames.positions.push_back(std::move(x));
    Cell c;
    c.a = {L * o.rc_angstrom, 0, 0}, c.b = {0, L * o.rc_angstrom, 0}, c.c = {0, 0, L * o.rc_angstrom};
    rep.frames.cells.push_back(c);
    rep.frames.timesteps.push_back(step);
  };
  store_frame(0);
  for (long step = 1; step <= o.steps; ++step) {
    for (int i = 0; i < n; ++i) {
      r[size_t(i)] = r[size_t(i)] + v[size_t(i)] * dt + f[size_t(i)] * (0.5 * dt * dt);
      vt[size_t(i)] = v[size_t(i)] + f[size_t(i)] * (lambda * dt);
      fold[size_t(i)] = f[size_t(i)];
    }
    forces(vt);
    double ke = 0;
    for (int i = 0; i < n; ++i) {
      v[size_t(i)] = v[size_t(i)] + (fold[size_t(i)] + f[size_t(i)]) * (0.5 * dt);
      ke += dot(v[size_t(i)], v[size_t(i)]);
    }
    const double kT = ke / (3.0 * (n - 1));
    if (step > o.equilibration && step % 10 == 0) {
      kts.push_back(kT);
      ps.push_back(n / V * kT + virial / (3 * V));
    }
    if (step % every == 0) {
      rep.order_series.push_back({double(step), order_now()});
      if (o.progress && !o.progress(step, kT)) { rep.notes.push_back("stopped at step " + std::to_string(step)); break; }
    }
    if (o.frame_every > 0 && step % o.frame_every == 0) store_frame(step);
  }
  auto block = [](const std::vector<double>& x, double& mean, double& err) {
    mean = err = 0;
    if (x.empty()) return;
    for (double y : x) mean += y;
    mean /= double(x.size());
    const int nb = 10;
    if (x.size() < size_t(2 * nb)) return;
    const size_t per = x.size() / nb;
    std::vector<double> bm;
    for (int b = 0; b < nb; ++b) { double s = 0; for (size_t k = 0; k < per; ++k) s += x[size_t(b) * per + k]; bm.push_back(s / double(per)); }
    double v2 = 0;
    for (double y : bm) v2 += (y - mean) * (y - mean);
    err = std::sqrt(v2 / (nb - 1) / nb);
  };
  block(kts, rep.kT, rep.kT_error);
  block(ps, rep.pressure, rep.pressure_error);
  rep.order = order_now();
  // domains of the first two types: cells about r_c wide, connected through faces (periodic)
  if (ta != tb) {
    const int g = std::max(2, int(std::floor(L)));
    const size_t nc = size_t(g) * g * g;
    std::vector<int> na(nc, 0), nb(nc, 0);
    for (int i = 0; i < n; ++i) {
      int c[3];
      for (int k = 0; k < 3; ++k) c[k] = std::min(g - 1, int((r[size_t(i)][size_t(k)] - L * std::floor(r[size_t(i)][size_t(k)] / L)) / L * g));
      const size_t id = size_t((c[0] * g + c[1]) * g + c[2]);
      if (type[size_t(i)] == ta) ++na[id];
      else if (type[size_t(i)] == tb) ++nb[id];
    }
    auto clusters = [&](bool first, int& count, double& mean_beads, double& largest) {
      std::vector<int> lab(nc, -1);
      std::vector<int> sizes, beads;
      size_t mine = 0;
      for (size_t s0 = 0; s0 < nc; ++s0) {
        const bool in = first ? na[s0] > nb[s0] : nb[s0] > na[s0];
        if (!in) continue;
        ++mine;
        if (lab[s0] >= 0) continue;
        const int L0 = int(sizes.size());
        sizes.push_back(0), beads.push_back(0);
        std::vector<size_t> st{s0};
        lab[s0] = L0;
        while (!st.empty()) {
          const size_t c = st.back();
          st.pop_back();
          ++sizes[size_t(L0)];
          beads[size_t(L0)] += first ? na[c] : nb[c];
          const int x = int(c / (size_t(g) * g)), y = int((c / size_t(g)) % size_t(g)), z = int(c % size_t(g));
          const int nbx[6][3] = {{1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}};
          for (const auto& d : nbx) {
            const size_t q = size_t((((x + d[0] + g) % g) * g + (y + d[1] + g) % g) * g + (z + d[2] + g) % g);
            const bool qin = first ? na[q] > nb[q] : nb[q] > na[q];
            if (qin && lab[q] < 0) lab[q] = L0, st.push_back(q);
          }
        }
      }
      count = int(sizes.size());
      mean_beads = count > 0 ? double(std::accumulate(beads.begin(), beads.end(), 0)) / count : 0;
      largest = mine > 0 && count > 0 ? double(*std::max_element(sizes.begin(), sizes.end())) / double(mine) : 0;
    };
    clusters(true, rep.domains_a, rep.domain_a_size, rep.largest_a);
    clusters(false, rep.domains_b, rep.domain_b_size, rep.largest_b);
  }
  // S(q) of the first type: q = 2π n / L, radially averaged over shells of unit width in |n|
  {
    std::vector<Vec3> ra;
    for (int i = 0; i < n; ++i) if (type[size_t(i)] == ta) ra.push_back(r[size_t(i)]);
    const int nmax = std::min(8, std::max(2, int(L / 2)));
    std::vector<double> sum(size_t(nmax + 1), 0.0), cnt(size_t(nmax + 1), 0.0);
    for (int hx = 0; hx <= nmax; ++hx)
      for (int hy = -nmax; hy <= nmax; ++hy)
        for (int hz = -nmax; hz <= nmax; ++hz) {
          if (hx == 0 && (hy < 0 || (hy == 0 && hz <= 0))) continue;   // half the vectors (S(−q) = S(q)), no q = 0
          const double hn = std::sqrt(double(hx * hx + hy * hy + hz * hz));
          const int shell = int(std::lround(hn));
          if (shell < 1 || shell > nmax) continue;
          const Vec3 q{2 * kPi * hx / L, 2 * kPi * hy / L, 2 * kPi * hz / L};
          std::complex<double> rho{0, 0};
          for (const auto& x : ra) rho += std::polar(1.0, dot(q, x));
          sum[size_t(shell)] += std::norm(rho) / double(ra.size());
          cnt[size_t(shell)] += 1;
        }
    double best = -1;
    for (int k = 1; k <= nmax; ++k) {
      if (cnt[size_t(k)] == 0) continue;
      rep.q.push_back(2 * kPi * k / L);
      rep.sq.push_back(sum[size_t(k)] / cnt[size_t(k)]);
      if (rep.sq.back() > best) best = rep.sq.back(), rep.q_peak = rep.q.back();
    }
    rep.spacing = rep.q_peak > 0 ? 2 * kPi / rep.q_peak : 0;
  }
  // the frames' topology: beads as atoms, bonds, molecules
  {
    System& t = rep.frames.topology;
    t.title = "DPD";
    for (int i = 0; i < n; ++i) {
      Atom a;
      a.id = i + 1;
      a.mol = mol[size_t(i)] + 1;
      a.element = element_of_type(type[size_t(i)]);
      a.name = std::string(1, type[size_t(i)]);
      a.type = type[size_t(i)] - 'A' + 1;
      a.pos = rep.frames.positions.front()[size_t(i)];
      t.atoms.push_back(a);
    }
    for (const auto& [i, j] : bonds) t.bonds.push_back({uint32_t(i), uint32_t(j), 1});
    for (char c : rep.types) { TypeInfo ti; ti.type = c - 'A' + 1; ti.mass = 1; ti.label = std::string(1, c); t.types.push_back(ti); }
    t.cell = rep.frames.cells.front();
    t.has_mol = true;
    t.unwrapped = true;
    t.source_format = "caps-dpd";
  }
  char b[240];
  std::snprintf(b, sizeof b, "%d beads in %d molecules, box %.2f r_c (ρ = %.2f) · kT %.3f ± %.3f · p %.2f ± %.2f · ψ %.2f · S(q) peak at q* = %.2f (spacing %.2f r_c)", n,
                rep.molecules, L, o.density, rep.kT, rep.kT_error, rep.pressure, rep.pressure_error, rep.order, rep.q_peak, rep.spacing);
  rep.notes.insert(rep.notes.begin(), b);
  std::string at;
  for (char x : used) for (char y : used) if (x <= y) { char q[48]; std::snprintf(q, sizeof q, "%sa_%c%c = %.2f", at.empty() ? "" : ", ", x, y, A[{x, y}]); at += q; }
  rep.notes.push_back("repulsions " + at + " (kT/r_c) · γ " + std::to_string(o.gamma).substr(0, 4) + " · dt " + std::to_string(o.dt).substr(0, 5));
  rep.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  return rep;
}

}  // namespace caps
