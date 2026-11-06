// CAPS force-field phonons: supercell force constants, dynamical matrix, dispersion, DOS, thermodynamics (see phonons.hpp).
#include "caps/phonons.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <stdexcept>

#include "caps/elements.hpp"
#include "caps/normal_modes.hpp"
#include "caps/spacegroup.hpp"
#include "parallel.hpp"

namespace caps {

namespace {
constexpr double kPi = 3.14159265358979323846;
constexpr double kWavenumber = 108.5913;      // cm⁻¹ per √(kcal/mol/Å²/(g/mol)), as normal_modes
constexpr double kCm1ToKcal = 2.859144e-3;    // kcal/mol per cm⁻¹
constexpr double kB = 0.0019872043;           // kcal/(mol·K)
constexpr double kC = 2.99792458e10;          // cm/s
using cd = std::complex<double>;

struct PairFC {
  size_t a = 0, b = 0;         // atoms of the cell (a: displaced; b: the supercell atom's site)
  double phi[3][3] = {};       // Φ_aα,Jβ
  std::vector<Vec3> d;         // r_J + L − r_a, the nearest image(s)
};
}  // namespace

PhononResult phonons(const System& s, const PhononOptions& o) {
  PhononResult r;
  if (!s.cell.valid()) throw std::invalid_argument("phonons: the crystal needs a periodic cell");
  if (!o.field_for) throw std::invalid_argument("phonons: no force field");
  const size_t N = s.atoms.size(), n = 3 * N;
  if (N == 0) throw std::invalid_argument("phonons: no atoms");
  auto stop = [&](const std::string& w, double f) {
    if (o.progress && !o.progress(w, f)) throw std::runtime_error("phonons stopped");
  };
  // ---- the supercell: each width at least twice the cut-off, unless given
  const double V = s.cell.volume();
  const double width[3] = {V / norm(cross(s.cell.b, s.cell.c)), V / norm(cross(s.cell.c, s.cell.a)), V / norm(cross(s.cell.a, s.cell.b))};
  std::array<int, 3> sc = o.supercell;
  const bool automatic = sc[0] <= 0 || sc[1] <= 0 || sc[2] <= 0;
  for (int d = 0; d < 3; ++d)
    if (sc[size_t(d)] <= 0) sc[size_t(d)] = std::max(1, int(std::ceil(2 * o.energy.cutoff / width[d] - 1e-9)));
  bool cut = false;
  if (automatic)
    while (size_t(sc[0]) * size_t(sc[1]) * size_t(sc[2]) * N > o.max_atoms && std::max({sc[0], sc[1], sc[2]}) > 1) {
      --sc[size_t(std::max_element(sc.begin(), sc.end()) - sc.begin())];
      cut = true;
    }
  if (cut) r.notes.push_back("the supercell was cut to stay under " + std::to_string(o.max_atoms) + " atoms: pairs beyond half its width are folded in");
  r.supercell = sc;
  const System sup = supercell(s, sc[0], sc[1], sc[2]);
  const size_t Ns = sup.atoms.size();
  r.supercell_atoms = Ns;
  const ForceField ff = o.field_for(sup);
  if (ff.type_index.size() != Ns) throw std::invalid_argument("phonons: the force field does not describe the supercell");
  r.field = ff.name;
  std::vector<double> mass(N);
  for (size_t a = 0; a < N; ++a) mass[a] = a < ff.mass.size() && ff.mass[a] > 0 ? ff.mass[a] : s.mass_of(s.atoms[a]);
  // ---- force constants by central differences (the cell's own atoms are the supercell's first N)
  Evaluator ev(ff, o.energy);
  std::vector<double> x(3 * Ns), f0, fp, fm;
  for (size_t i = 0; i < Ns; ++i)
    for (int k = 0; k < 3; ++k) x[3 * i + size_t(k)] = sup.atoms[i].pos[size_t(k)];
  // the supercell positions from the cell's own (copy g of atom a: sup index g·N + a, shifted by its lattice offset)
  const size_t ng = Ns / N;
  std::vector<Vec3> offset(ng);
  for (size_t g = 0; g < ng; ++g) offset[g] = sup.atoms[g * N].pos - s.atoms[0].pos;
  std::vector<Vec3> X(N);
  for (size_t a = 0; a < N; ++a) X[a] = s.atoms[a].pos;
  auto place = [&] {
    for (size_t g = 0; g < ng; ++g)
      for (size_t a = 0; a < N; ++a)
        for (int k = 0; k < 3; ++k) x[3 * (g * N + a) + size_t(k)] = X[a][size_t(k)] + offset[g][size_t(k)];
  };
  // each atom's force in the periodic crystal: the mean over its copies (equal by symmetry, summed for precision)
  auto cell_forces = [&](std::vector<Vec3>& F) {
    place();
    ev.compute(x, sup.cell, f0);
    F.assign(N, Vec3{0, 0, 0});
    for (size_t g = 0; g < ng; ++g)
      for (size_t a = 0; a < N; ++a)
        for (int k = 0; k < 3; ++k) F[a][size_t(k)] += f0[3 * (g * N + a) + size_t(k)] / double(ng);
    double fmax = 0;
    for (const auto& v : F) fmax = std::max(fmax, norm(v));
    return fmax;
  };
  std::vector<Vec3> F;
  double fmax = cell_forces(F);
  if (o.relax && fmax > o.relax_ftol) {   // FIRE on the cell's coordinates (masses as the atoms'; Å, fs)
    constexpr double kAcc = 4.184e-4;     // Å/fs² per (kcal/mol/Å)/(g/mol)
    std::vector<Vec3> v(N, Vec3{0, 0, 0});
    double dt = 0.5, alpha = 0.1;
    int positive = 0;
    const std::vector<Vec3> X0 = X;
    for (int it = 0; it < o.relax_iterations && fmax > o.relax_ftol; ++it) {
      if (it % 50 == 0) stop("relaxing the atoms", 0);
      double P = 0, vn = 0, fn = 0;
      for (size_t a = 0; a < N; ++a) P += dot(F[a], v[a]), vn += dot(v[a], v[a]), fn += dot(F[a], F[a]);
      if (P > 0) {
        const double sv = std::sqrt(vn), sf = std::sqrt(fn) > 0 ? std::sqrt(fn) : 1;
        for (size_t a = 0; a < N; ++a) v[a] = v[a] * (1 - alpha) + F[a] * (alpha * sv / sf);
        if (++positive > 5) dt = std::min(dt * 1.1, 4.0), alpha *= 0.99;
      } else {
        for (auto& w : v) w = Vec3{0, 0, 0};
        dt *= 0.5, alpha = 0.1, positive = 0;
      }
      for (size_t a = 0; a < N; ++a) {
        v[a] = v[a] + F[a] * (dt * kAcc / mass[a]);
        Vec3 dx = v[a] * dt;
        const double l = norm(dx);
        if (l > 0.1) dx = dx * (0.1 / l);
        X[a] = X[a] + dx;
      }
      fmax = cell_forces(F);
      r.relax_steps = it + 1;
    }
    for (size_t a = 0; a < N; ++a) r.relax_shift = std::max(r.relax_shift, norm(X[a] - X0[a]));
    char b[200];
    std::snprintf(b, sizeof b, "atoms relaxed in the supercell's periodic copies: %d steps, largest move %.3f Å, forces below %.2g kcal/mol/Å", r.relax_steps, r.relax_shift, fmax);
    r.notes.push_back(b);
  }
  place();
  r.max_force = fmax;
  r.positions = X;
  r.static_energy = ev.compute(x, sup.cell, f0).total() / double(ng);
  ev.freeze_pairs(x, sup.cell);
  std::vector<std::vector<double>> Phi(n, std::vector<double>(3 * Ns));   // row aα, column Jβ
  for (size_t a = 0; a < N; ++a)
    for (int al = 0; al < 3; ++al) {
      stop("force constants", double(3 * a + size_t(al)) / double(n));
      const size_t xa = 3 * a + size_t(al);
      const double x0 = x[xa];
      x[xa] = x0 + o.step;
      ev.compute(x, sup.cell, fp);
      x[xa] = x0 - o.step;
      ev.compute(x, sup.cell, fm);
      x[xa] = x0;
      for (size_t c = 0; c < 3 * Ns; ++c) Phi[xa][c] = -(fp[c] - fm[c]) / (2 * o.step);
    }
  ev.unfreeze();
  // acoustic sum rule on the self term
  for (size_t a = 0; a < N; ++a)
    for (int al = 0; al < 3; ++al)
      for (int be = 0; be < 3; ++be) {
        double sum = 0;
        for (size_t J = 0; J < Ns; ++J)
          if (J != a) sum += Phi[3 * a + size_t(al)][3 * J + size_t(be)];
        Phi[3 * a + size_t(al)][3 * a + size_t(be)] = -sum;
      }
  // ---- pairs with their nearest supercell images
  const Cell& C = sup.cell;
  double phimax = 0;
  for (const auto& row : Phi) for (double v : row) phimax = std::max(phimax, std::fabs(v));
  std::vector<PairFC> pairs;
  for (size_t a = 0; a < N; ++a)
    for (size_t J = 0; J < Ns; ++J) {
      PairFC p;
      p.a = a, p.b = J % N;
      double big = 0;
      for (int al = 0; al < 3; ++al)
        for (int be = 0; be < 3; ++be) big = std::max(big, std::fabs(p.phi[al][be] = Phi[3 * a + size_t(al)][3 * J + size_t(be)]));
      if (big <= 1e-9 * phimax && J != a) continue;
      const Vec3 pj{x[3 * J], x[3 * J + 1], x[3 * J + 2]}, pa{x[3 * a], x[3 * a + 1], x[3 * a + 2]};
      Vec3 fr = C.to_fractional(pj - pa + C.origin);
      for (int k = 0; k < 3; ++k) fr[size_t(k)] -= std::round(fr[size_t(k)]);
      double best = 1e300;
      std::vector<Vec3> cand;
      for (int i = -1; i <= 1; ++i)
        for (int j = -1; j <= 1; ++j)
          for (int k = -1; k <= 1; ++k) {
            const Vec3 d = C.a * (fr[0] + i) + C.b * (fr[1] + j) + C.c * (fr[2] + k);
            const double l = norm(d);
            if (l < best - 1e-4) best = l, cand.clear();
            if (l <= best + 1e-4) cand.push_back(d);
          }
      p.d = std::move(cand);
      pairs.push_back(std::move(p));
    }
  // ---- the dynamical matrix at k (Cartesian, Å⁻¹): wavenumbers ascending, and (optionally) each column's weight per atom
  const double vol = s.cell.volume();
  const Vec3 B0 = cross(s.cell.b, s.cell.c) * (2 * kPi / vol), B1 = cross(s.cell.c, s.cell.a) * (2 * kPi / vol), B2 = cross(s.cell.a, s.cell.b) * (2 * kPi / vol);
  auto kcart = [&](const std::array<double, 3>& f) { return B0 * f[0] + B1 * f[1] + B2 * f[2]; };
  auto solve = [&](const Vec3& k, std::vector<double>& w, std::vector<std::vector<double>>* atom_weight) {
    std::vector<cd> D(n * n, cd(0, 0));
    for (const auto& p : pairs) {
      cd ph(0, 0);
      for (const auto& d : p.d) ph += std::polar(1.0, dot(k, d));
      ph /= double(p.d.size());
      const double im = 1 / std::sqrt(mass[p.a] * mass[p.b]);
      for (int al = 0; al < 3; ++al)
        for (int be = 0; be < 3; ++be) D[(3 * p.a + size_t(al)) * n + 3 * p.b + size_t(be)] += p.phi[al][be] * im * ph;
    }
    // Hermitian part, then the real symmetric form of twice the size
    const size_t m2 = 2 * n;
    std::vector<double> A(m2 * m2);
    for (size_t i = 0; i < n; ++i)
      for (size_t j = 0; j < n; ++j) {
        const cd h = 0.5 * (D[i * n + j] + std::conj(D[j * n + i]));
        A[i * m2 + j] = h.real(), A[(i + n) * m2 + j + n] = h.real();
        A[i * m2 + j + n] = -h.imag(), A[(i + n) * m2 + j] = h.imag();
      }
    std::vector<double> ev2;
    symmetric_eigen(A, m2, ev2);
    w.resize(n);
    for (size_t i = 0; i < n; ++i) {
      const double lam = 0.5 * (ev2[2 * i] + ev2[2 * i + 1]);
      w[i] = (lam < 0 ? -1.0 : 1.0) * kWavenumber * std::sqrt(std::fabs(lam));
    }
    if (atom_weight) {   // per column (each mode appears twice, weight ½): its wavenumber last, then |e|² per atom
      atom_weight->assign(m2, std::vector<double>(N + 1, 0.0));
      for (size_t c = 0; c < m2; ++c) {
        auto& aw = (*atom_weight)[c];
        aw[N] = (ev2[c] < 0 ? -1.0 : 1.0) * kWavenumber * std::sqrt(std::fabs(ev2[c]));
        for (size_t row = 0; row < n; ++row) {
          const double u = A[row * m2 + c], v = A[(row + n) * m2 + c];
          aw[row / 3] += 0.5 * (u * u + v * v);
        }
      }
    }
  };
  // ---- the path
  std::vector<std::array<double, 3>> pts = o.path;
  std::vector<std::string> lab = o.labels;
  if (pts.size() < 2) {
    pts = {{0, 0, 0}, {0.5, 0, 0}, {0.5, 0.5, 0}, {0, 0.5, 0}, {0, 0, 0}, {0, 0, 0.5}, {0.5, 0.5, 0.5}};
    lab = {"Γ", "X", "S", "Y", "Γ", "Z", "R"};
  }
  lab.resize(pts.size());
  double total = 0;
  std::vector<double> seg(pts.size() - 1);
  for (size_t i = 0; i + 1 < pts.size(); ++i) total += seg[i] = norm(kcart(pts[i + 1]) - kcart(pts[i]));
  if (total <= 0) throw std::invalid_argument("phonons: the path has no length");
  std::vector<std::pair<double, std::array<double, 3>>> kp;   // (x, fractional k)
  double xs = 0;
  for (size_t i = 0; i + 1 < pts.size(); ++i) {
    r.tick_x.push_back(xs), r.tick_labels.push_back(lab[i]);
    const int m = std::max(2, int(std::round(o.path_points * seg[i] / total)));
    for (int t = 0; t <= m; ++t) {
      if (t == 0 && i > 0) continue;
      const double u = double(t) / m;
      kp.push_back({xs + u * seg[i], {pts[i][0] + u * (pts[i + 1][0] - pts[i][0]), pts[i][1] + u * (pts[i + 1][1] - pts[i][1]), pts[i][2] + u * (pts[i + 1][2] - pts[i][2])}});
    }
    xs += seg[i];
  }
  r.tick_x.push_back(xs), r.tick_labels.push_back(lab.back());
  ThreadPool pool(default_threads());
  r.path_x.resize(kp.size()), r.bands.resize(kp.size());
  stop("dispersion", 0);
  pool.run(kp.size(), [&](int, size_t a, size_t b) {
    for (size_t i = a; i < b; ++i) { r.path_x[i] = kp[i].first; solve(kcart(kp[i].second), r.bands[i], nullptr); }
  });
  // ---- the mesh: DOS, per element, thermodynamics
  std::array<int, 3> mesh = o.mesh;
  const double len[3] = {norm(s.cell.a), norm(s.cell.b), norm(s.cell.c)};
  for (int d = 0; d < 3; ++d)
    if (mesh[size_t(d)] <= 0) mesh[size_t(d)] = std::clamp(int(std::ceil(40.0 / len[d])), 1, 30);
  r.mesh = mesh;
  const size_t nk = size_t(mesh[0]) * size_t(mesh[1]) * size_t(mesh[2]);
  std::vector<std::vector<double>> wk(nk);
  std::vector<std::vector<std::vector<double>>> aw(nk);
  stop("density of states", 0);
  pool.run(nk, [&](int, size_t a, size_t b) {
    for (size_t q = a; q < b; ++q) {
      const int i = int(q / size_t(mesh[1] * mesh[2])), j = int(q / size_t(mesh[2]) % size_t(mesh[1])), k = int(q % size_t(mesh[2]));
      solve(kcart({double(i) / mesh[0], double(j) / mesh[1], double(k) / mesh[2]}), wk[q], &aw[q]);
    }
  });
  double wmin = 0, wmax = 0;
  for (size_t q = 0; q < nk; ++q)
    for (size_t i = 0; i < n; ++i) {
      if (q == 0 && i < 3) continue;   // Γ's acoustic modes
      wmin = std::min(wmin, wk[q][i]), wmax = std::max(wmax, wk[q][i]);
      r.imaginary += wk[q][i] < -10;
    }
  for (const auto& b : r.bands) wmin = std::min(wmin, b.size() > 3 ? b[3] : 0.0);
  r.min_wavenumber = wmin;
  // the elements
  std::vector<int> el_of(N), els;
  for (size_t a = 0; a < N; ++a) {
    const int z = s.atoms[a].element;
    auto it = std::find(els.begin(), els.end(), z);
    if (it == els.end()) els.push_back(z), it = els.end() - 1;
    el_of[a] = int(it - els.begin());
  }
  for (int z : els) r.pdos_labels.push_back(element(z).symbol);
  const double sig = std::max(0.1, o.dos_sigma), st = std::max(0.05, o.dos_step);
  const double lo = std::floor((std::min(0.0, wmin) - 4 * sig) / st) * st, hi = wmax + 4 * sig;
  const size_t nb = size_t((hi - lo) / st) + 1;
  r.dos_w.resize(nb), r.dos.assign(nb, 0.0), r.pdos.assign(els.size(), std::vector<double>(nb, 0.0));
  for (size_t i = 0; i < nb; ++i) r.dos_w[i] = lo + st * double(i);
  const double g0 = 1 / (sig * std::sqrt(2 * kPi)) / double(nk);
  for (size_t q = 0; q < nk; ++q)
    for (const auto& col : aw[q]) {
      const double w0 = col[N];
      const size_t i0 = size_t(std::max(0.0, (w0 - 5 * sig - lo) / st)), i1 = std::min(nb - 1, size_t(std::max(0.0, (w0 + 5 * sig - lo) / st)));
      for (size_t i = i0; i <= i1; ++i) {
        const double u = (r.dos_w[i] - w0) / sig, g = g0 * std::exp(-0.5 * u * u);
        double tot = 0;
        for (size_t a = 0; a < N; ++a) tot += col[a], r.pdos[size_t(el_of[a])][i] += col[a] * g;
        r.dos[i] += tot * g;
      }
    }
  // thermodynamics per cell
  for (double T : o.temperatures) {
    if (T < 0) continue;
    double zpe = 0, U = 0, F = 0, Cv = 0;
    for (size_t q = 0; q < nk; ++q)
      for (size_t i = 0; i < n; ++i) {
        const double w0 = wk[q][i];
        if (w0 < 0.5) continue;   // Γ's acoustic modes and unstable ones
        const double e = w0 * kCm1ToKcal;
        zpe += 0.5 * e;
        if (T == 0) { U += 0.5 * e, F += 0.5 * e; continue; }   // the ground state: no entropy, no heat capacity
        const double xx = e / (kB * T);
        U += e * (0.5 + 1 / std::expm1(xx));
        F += 0.5 * e + kB * T * std::log(-std::expm1(-xx));
        const double ex = std::exp(-xx);
        Cv += kB * xx * xx * ex / ((1 - ex) * (1 - ex));
      }
    r.T.push_back(T), r.zpe.push_back(zpe / double(nk)), r.U.push_back(U / double(nk)), r.F.push_back(F / double(nk));
    r.S.push_back(T > 0 ? (U - F) / double(nk) / T * 1000 : 0.0), r.Cv.push_back(Cv / double(nk) * 1000);
  }
  // sound speeds: the acoustic slopes at 1 % of the way to the zone face along a*, b*, c*
  for (int d = 0; d < 3; ++d) {
    std::array<double, 3> f{0, 0, 0};
    f[size_t(d)] = 0.01;
    const Vec3 k = kcart(f);
    std::vector<double> w;
    solve(k, w, nullptr);
    for (int b = 0; b < 3; ++b) r.sound[size_t(d)][size_t(b)] = 2 * kPi * kC * std::max(0.0, w[size_t(b)]) / (norm(k) * 1e10);
  }
  {
    char b[300];
    std::snprintf(b, sizeof b, "%s · supercell %d × %d × %d (%zu atoms), ±%.3g Å · %zu force-constant pairs · mesh %d × %d × %d", ff.name.c_str(), sc[0], sc[1], sc[2], Ns, o.step,
                  pairs.size(), mesh[0], mesh[1], mesh[2]);
    r.notes.insert(r.notes.begin(), b);
    if (r.max_force > 0.5) {
      std::snprintf(b, sizeof b, "the largest force at the reference is %.2f kcal/mol/Å: minimise the crystal (with its cell) first, or the low modes are off", r.max_force);
      r.notes.push_back(b);
    }
    if (r.imaginary > 0) r.notes.push_back(std::to_string(r.imaginary) + " mesh modes below −10 cm⁻¹: not at a minimum, or the supercell is too small");
  }
  return r;
}

}  // namespace caps
