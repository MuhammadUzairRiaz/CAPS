// CAPS dielectric, piezoelectric and optical response of a crystal (see dielectric.hpp).
#include "caps/dielectric.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <memory>
#include <stdexcept>

#include "caps/crystalfit.hpp"
#include "caps/normal_modes.hpp"
#include "caps/spacegroup.hpp"

namespace caps {

namespace {
constexpr double kPi = 3.14159265358979323846;
constexpr double kE = 1.602176634e-19, kEps0 = 8.8541878128e-12, kNA = 6.02214076e23, kC = 2.99792458e8;
constexpr double kKcalA2 = 4184.0 / 6.02214076e23 / 1e-20;   // J/m² per kcal/mol/Å²
using cd = std::complex<double>;

// the cell's particles in their periodic copies (copy g of particle a is supercell atom g·N + a)
struct Copies {
  size_t N = 0;
  std::array<int, 3> sc{1, 1, 1};
  std::vector<std::array<int, 3>> off;
  ForceField ff;
  std::unique_ptr<Evaluator> ev;
  // the mean force on each particle of the cell over its copies (kcal/mol/Å), 3N
  std::vector<double> forces(const std::vector<Vec3>& X, const Cell& c) {
    const size_t ng = off.size();
    std::vector<double> x(3 * N * ng), f;
    for (size_t g = 0; g < ng; ++g) {
      const Vec3 o = c.a * double(off[g][0]) + c.b * double(off[g][1]) + c.c * double(off[g][2]);
      for (size_t i = 0; i < N; ++i)
        for (int k = 0; k < 3; ++k) x[3 * (g * N + i) + size_t(k)] = X[i][size_t(k)] + o[size_t(k)];
    }
    Cell sup;
    sup.a = c.a * double(sc[0]), sup.b = c.b * double(sc[1]), sup.c = c.c * double(sc[2]);
    ev->compute(x, sup, f);
    std::vector<double> F(3 * N, 0.0);
    for (size_t g = 0; g < ng; ++g)
      for (size_t i = 0; i < 3 * N; ++i) F[i] += f[3 * g * N + i] / double(ng);
    return F;
  }
};

// pseudo-inverse of a symmetric matrix (eigenvalues below tol × the largest dropped)
std::vector<double> pinv(std::vector<double> A, size_t n, double tol = 1e-9) {
  std::vector<double> w;
  symmetric_eigen(A, n, w);
  double mx = 0;
  for (double v : w) mx = std::max(mx, std::fabs(v));
  std::vector<double> P(n * n, 0.0);
  for (size_t k = 0; k < n; ++k) {
    if (std::fabs(w[k]) <= tol * mx) continue;
    for (size_t i = 0; i < n; ++i)
      for (size_t j = 0; j < n; ++j) P[i * n + j] += A[i * n + k] * A[j * n + k] / w[k];
  }
  return P;
}
}  // namespace

DielectricResult crystal_dielectric(const System& s0, const DielectricOptions& o) {
  DielectricResult r;
  if (!s0.cell.valid()) throw std::invalid_argument("dielectric: the crystal needs a periodic cell");
  if (!o.field_for) throw std::invalid_argument("dielectric: no force field");
  const size_t N = s0.atoms.size();
  auto stop = [&](const std::string& w, double f) {
    if (o.progress && !o.progress(w, f)) throw std::runtime_error("dielectric stopped");
  };
  std::vector<char> shell(N, 0);
  for (size_t i = 0; i < N && i < o.shell.size(); ++i) shell[i] = o.shell[i];
  // ---- relaxation (cell at zero stress) and the relaxed-ion elastic constants, by the crystal fitter with no parameters
  stop("relaxing", 0);
  System s = s0;
  {
    FitOptions fo;
    fo.energy = o.energy;
    fo.relaxed = o.relax && o.relax_cell;
    FitStructure fs;
    fs.cell = s0, fs.field_for = o.field_for, fs.shell = shell;
    fo.structures.push_back(fs);
    for (int I = 1; I <= 6; ++I)
      for (int J = I; J <= 6; ++J) fo.targets.push_back({0, "C" + std::to_string(I) + std::to_string(J), 0, 0});
    const FitResult fr = fit_to_crystals(fo);
    r.C = fr.structures[0].C;
    if (fo.relaxed) {
      s.cell = fr.structures[0].cell;
      for (size_t i = 0; i < N; ++i) s.atoms[i].pos = fr.structures[0].positions[i];
    }
    for (const auto& nt : fr.notes)
      if (nt.rfind("Ewald grid", 0) != 0) r.notes.push_back(nt);   // the fitter's own grid note does not concern this
  }
  // ---- the copies
  Copies cp;
  cp.N = N;
  const double V0 = s.cell.volume();
  const double width[3] = {V0 / norm(cross(s.cell.b, s.cell.c)), V0 / norm(cross(s.cell.c, s.cell.a)), V0 / norm(cross(s.cell.a, s.cell.b))};
  for (int d = 0; d < 3; ++d) cp.sc[size_t(d)] = std::max(1, int(std::ceil(2.02 * o.energy.cutoff / width[d] - 1e-9)));
  if (size_t(cp.sc[0]) * size_t(cp.sc[1]) * size_t(cp.sc[2]) * N > 15000) throw std::invalid_argument("dielectric: the copies to twice the cut-off exceed 15,000 atoms");
  const System sup = supercell(s, cp.sc[0], cp.sc[1], cp.sc[2]);
  cp.ff = o.field_for(sup);
  if (cp.ff.type_index.size() != sup.atoms.size()) throw std::invalid_argument("dielectric: the force field does not describe the copies");
  EnergyOptions eo = o.energy;
  if (eo.electrostatics == EnergyOptions::Electrostatics::PME) eo.pme_spacing = std::min(eo.pme_spacing, 0.5);
  cp.ev = std::make_unique<Evaluator>(cp.ff, eo);
  {
    const size_t ng = sup.atoms.size() / N;
    cp.off.resize(ng);
    for (size_t g = 0; g < ng; ++g) {
      const Vec3 f = s.cell.to_fractional(sup.atoms[g * N].pos - s.atoms[0].pos);
      for (int k = 0; k < 3; ++k) cp.off[g][size_t(k)] = int(std::lround(f[size_t(k)]));
    }
  }
  std::vector<double> q(N), mass(N);
  double qsum = 0;
  for (size_t i = 0; i < N; ++i) q[i] = cp.ff.charge[i], mass[i] = cp.ff.mass[i] > 0 ? cp.ff.mass[i] : s.mass_of(s.atoms[i]), qsum += q[i];
  if (std::fabs(qsum) > 1e-6) r.notes.push_back("the cell carries a net charge of " + std::to_string(qsum) + " e");
  std::vector<Vec3> X(N);
  for (size_t i = 0; i < N; ++i) X[i] = s.atoms[i].pos;
  const size_t n = 3 * N;
  auto hessian = [&](const std::vector<Vec3>& Xh) {
    std::vector<double> K(n * n, 0.0);
    for (size_t a = 0; a < N; ++a)
      for (int k = 0; k < 3; ++k) {
        stop("force constants", double(3 * a + size_t(k)) / double(n));
        auto Xp = Xh, Xm = Xh;
        Xp[a][size_t(k)] += o.step, Xm[a][size_t(k)] -= o.step;
        const auto fp = cp.forces(Xp, s.cell), fm = cp.forces(Xm, s.cell);
        for (size_t i = 0; i < n; ++i) K[i * n + 3 * a + size_t(k)] = -(fp[i] - fm[i]) / (2 * o.step);
      }
    for (size_t i = 0; i < n; ++i)
      for (size_t j = i + 1; j < n; ++j) K[i * n + j] = K[j * n + i] = 0.5 * (K[i * n + j] + K[j * n + i]);
    return K;
  };
  std::vector<double> K = hessian(X);
  // ---- Newton steps on the particles at this cell (the start's force constants), then the force constants again
  if (o.relax) {
    const std::vector<double> K0p = pinv(K, n);
    for (int it = 0; it < 8; ++it) {
      const auto F = cp.forces(X, s.cell);
      double fmax = 0;
      for (double v : F) fmax = std::max(fmax, std::fabs(v));
      if (fmax < 1e-7) break;
      for (size_t i = 0; i < n; ++i) {
        double dx = 0;
        for (size_t j = 0; j < n; ++j) dx += K0p[i * n + j] * F[j];
        X[i / 3][i % 3] += dx;
      }
    }
    K = hessian(X);
  }
  {
    const auto F = cp.forces(X, s.cell);
    for (size_t i = 0; i < N; ++i) r.max_force = std::max(r.max_force, std::sqrt(F[3 * i] * F[3 * i] + F[3 * i + 1] * F[3 * i + 1] + F[3 * i + 2] * F[3 * i + 2]));
    if (r.max_force > 1e-2) r.notes.push_back("the particles are not at a minimum (largest force " + std::to_string(r.max_force) + " kcal/mol/Å): the response is not meaningful");
  }
  r.positions = X, r.cell_vectors = s.cell, r.cell = cell_parameters(s.cell), r.volume = V0;
  const double V = V0 * 1e-30;
  // ---- ε⁰ over every particle
  const std::vector<double> Kp = pinv(K, n);
  auto response = [&](const std::vector<size_t>& idx, const std::vector<double>& Kinv) {   // 1 + Σ q q K⁻¹ / ε₀V over idx
    const size_t m = idx.size();
    Mat3x3 out{};
    for (int a = 0; a < 3; ++a)
      for (int b = 0; b < 3; ++b) {
        double sum = 0;
        for (size_t i = 0; i < m; ++i)
          for (size_t j = 0; j < m; ++j) sum += q[idx[i]] * q[idx[j]] * Kinv[(3 * i + size_t(a)) * 3 * m + 3 * j + size_t(b)];
        out[size_t(a)][size_t(b)] = (a == b ? 1.0 : 0.0) + sum * kE * kE / kKcalA2 / (kEps0 * V);
      }
    return out;
  };
  std::vector<size_t> all(N);
  for (size_t i = 0; i < N; ++i) all[i] = i;
  r.eps0 = response(all, Kp);
  // ---- ε^∞ over the shells (cores held)
  std::vector<size_t> sh, co;
  for (size_t i = 0; i < N; ++i) (shell[i] ? sh : co).push_back(i);
  r.shells = !sh.empty();
  const size_t ns = sh.size(), nc = co.size();
  std::vector<double> Kssi;
  if (ns) {
    std::vector<double> Kss(9 * ns * ns, 0.0);
    for (size_t i = 0; i < ns; ++i)
      for (size_t j = 0; j < ns; ++j)
        for (int a = 0; a < 3; ++a)
          for (int b = 0; b < 3; ++b) Kss[(3 * i + size_t(a)) * 3 * ns + 3 * j + size_t(b)] = K[(3 * sh[i] + size_t(a)) * n + 3 * sh[j] + size_t(b)];
    Kssi = pinv(Kss, 3 * ns);
    r.epsinf = response(sh, Kssi);
  } else {
    for (int a = 0; a < 3; ++a) r.epsinf[size_t(a)][size_t(a)] = 1;
    r.notes.push_back("rigid ions: ε^∞ = 1 and no refractive indices (a shell model gives the electrons' part)");
  }
  auto principal = [](const Mat3x3& M, std::array<double, 3>& val, std::array<Vec3, 3>* axes) {
    std::vector<double> A(9), w;
    for (int a = 0; a < 3; ++a)
      for (int b = 0; b < 3; ++b) A[size_t(3 * a + b)] = M[size_t(a)][size_t(b)];
    symmetric_eigen(A, 3, w);
    for (int k = 0; k < 3; ++k) {
      val[size_t(k)] = w[size_t(k)];
      if (axes) (*axes)[size_t(k)] = Vec3{A[size_t(0 * 3 + k)], A[size_t(1 * 3 + k)], A[size_t(2 * 3 + k)]};
    }
  };
  principal(r.eps0, r.eps0_principal, &r.axes);
  principal(r.epsinf, r.epsinf_principal, nullptr);
  for (int k = 0; k < 3; ++k) r.refractive[size_t(k)] = r.shells ? std::sqrt(std::max(0.0, r.epsinf_principal[size_t(k)])) : 0.0;
  // ---- the cores' Γ modes with the shells condensed out, oscillator strengths, ε(ω) and reflectivity
  {
    const size_t m3 = 3 * nc;
    std::vector<double> Keff(m3 * m3), Z(m3 * 3, 0.0);
    for (size_t i = 0; i < m3; ++i)
      for (size_t j = 0; j < m3; ++j) Keff[i * m3 + j] = K[(3 * co[i / 3] + i % 3) * n + 3 * co[j / 3] + j % 3];
    for (size_t i = 0; i < nc; ++i)
      for (int a = 0; a < 3; ++a) Z[(3 * i + size_t(a)) * 3 + size_t(a)] = q[co[i]];
    std::vector<double> mc(nc);
    for (size_t i = 0; i < nc; ++i) mc[i] = mass[co[i]];
    if (ns) {
      for (size_t k = 0; k < ns; ++k) {   // each shell's mass to its nearest core
        size_t best = 0;
        double bd = 1e300;
        for (size_t i = 0; i < nc; ++i) {
          const double dd = norm(s.cell.minimum_image(X[sh[k]] - X[co[i]]));
          if (dd < bd) bd = dd, best = i;
        }
        mc[best] += mass[sh[k]];
      }
      const size_t s3 = 3 * ns;
      std::vector<double> T(m3 * s3, 0.0);   // K_cs K_ss⁻¹
      for (size_t i = 0; i < m3; ++i)
        for (size_t j = 0; j < s3; ++j) {
          double sum = 0;
          for (size_t k = 0; k < s3; ++k) sum += K[(3 * co[i / 3] + i % 3) * n + 3 * sh[k / 3] + k % 3] * Kssi[k * s3 + j];
          T[i * s3 + j] = sum;
        }
      for (size_t i = 0; i < m3; ++i) {
        for (size_t j = 0; j < m3; ++j) {
          double sum = 0;
          for (size_t k = 0; k < s3; ++k) sum += T[i * s3 + k] * K[(3 * sh[k / 3] + k % 3) * n + 3 * co[j / 3] + j % 3];
          Keff[i * m3 + j] -= sum;
        }
        for (int b = 0; b < 3; ++b) {
          double sum = 0;
          for (size_t k = 0; k < ns; ++k) sum += T[i * s3 + 3 * k + size_t(b)] * q[sh[k]];
          Z[i * 3 + size_t(b)] -= sum;
        }
      }
    }
    std::vector<double> D(Keff);
    for (size_t i = 0; i < m3; ++i)
      for (size_t j = 0; j < m3; ++j) D[i * m3 + j] /= std::sqrt(mc[i / 3] * mc[j / 3]);
    std::vector<double> lam;
    symmetric_eigen(D, m3, lam);
    std::vector<size_t> order(m3);
    for (size_t k = 0; k < m3; ++k) order[k] = k;
    std::sort(order.begin(), order.end(), [&](size_t a, size_t b) { return std::fabs(lam[a]) < std::fabs(lam[b]); });
    struct Osc { double w2; Vec3 p; };   // ω² (rad²/s²), mode charge (C/√kg)
    std::vector<Osc> osc;
    for (size_t kk = 3; kk < m3; ++kk) {
      const size_t k = order[kk];
      const double w2 = lam[k] * 4184.0 / (1e-3 * 1e-20);
      r.modes.push_back((lam[k] < 0 ? -1 : 1) * std::sqrt(std::fabs(w2)) / (2 * kPi * kC * 100));
      Vec3 p{0, 0, 0};
      for (size_t i = 0; i < m3; ++i)
        for (int b = 0; b < 3; ++b) p[size_t(b)] += Z[i * 3 + size_t(b)] * kE * D[i * m3 + k] / std::sqrt(mc[i / 3] * 1e-3 / kNA);
      if (w2 > 0) osc.push_back({w2, p});
    }
    std::sort(r.modes.begin(), r.modes.end());
    const double wmax = o.w_max > 0 ? o.w_max : 1.2 * (r.modes.empty() ? 1000.0 : std::max(1.0, r.modes.back()));
    const double gam = 2 * kPi * kC * 100 * o.damping;
    for (int k = 0; k <= o.points; ++k) {
      const double wn = wmax * k / o.points, w = 2 * kPi * kC * 100 * wn;
      cd E[3][3];
      for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b) E[a][b] = r.epsinf[size_t(a)][size_t(b)];
      for (const auto& m : osc) {
        const cd den = cd(m.w2 - w * w, -gam * w);
        for (int a = 0; a < 3; ++a)
          for (int b = 0; b < 3; ++b) E[a][b] += m.p[size_t(a)] * m.p[size_t(b)] / (kEps0 * V) / den;
      }
      r.w.push_back(wn);
      for (int ax = 0; ax < 3; ++ax) {
        const Vec3& u = r.axes[size_t(ax)];
        cd e = 0;
        for (int a = 0; a < 3; ++a)
          for (int b = 0; b < 3; ++b) e += u[size_t(a)] * E[a][b] * u[size_t(b)];
        const cd rt = std::sqrt(e), rr = (rt - 1.0) / (rt + 1.0);
        r.eps_re[ax].push_back(e.real()), r.eps_im[ax].push_back(e.imag()), r.reflectivity[ax].push_back(std::norm(rr));
      }
    }
  }
  // ---- piezoelectric constants from the internal-strain forces
  if (o.piezo) {
    stop("internal strain", 0.9);
    r.piezo = true;
    for (int J = 0; J < 6; ++J) {
      std::vector<double> Fs[2];
      for (int sg = 0; sg < 2; ++sg) {
        const double dd = (sg ? -o.strain : o.strain) * (J < 3 ? 1.0 : 0.5);
        double e6[6] = {0, 0, 0, 0, 0, 0};
        e6[J] = dd;
        auto Fm = [&](const Vec3& v) {
          return Vec3{v[0] + e6[0] * v[0] + e6[5] * v[1] + e6[4] * v[2], v[1] + e6[5] * v[0] + e6[1] * v[1] + e6[3] * v[2], v[2] + e6[4] * v[0] + e6[3] * v[1] + e6[2] * v[2]};
        };
        Cell c2 = s.cell;
        c2.a = Fm(s.cell.a), c2.b = Fm(s.cell.b), c2.c = Fm(s.cell.c);
        std::vector<Vec3> X2(N);
        for (size_t i = 0; i < N; ++i) X2[i] = Fm(X[i]);
        Fs[sg] = cp.forces(X2, c2);
      }
      // Θ = ∂F/∂ε, δu = K⁺Θ (Å per unit strain), e_αJ = (1/V) Σ q δu
      std::vector<double> th(n), du(n, 0.0);
      for (size_t i = 0; i < n; ++i) th[i] = (Fs[0][i] - Fs[1][i]) / (2 * o.strain);
      for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < n; ++j) du[i] += Kp[i * n + j] * th[j];
      for (int a = 0; a < 3; ++a) {
        double P = 0;
        for (size_t i = 0; i < N; ++i) P += q[i] * du[3 * i + size_t(a)];
        r.e[size_t(a)][size_t(J)] = P * kE * 1e-10 / V;   // C/m²
      }
    }
    // d = e S (S = C⁻¹ in 1/GPa) → pC/N
    double A[6][12];
    for (int i = 0; i < 6; ++i)
      for (int j = 0; j < 12; ++j) A[i][j] = j < 6 ? r.C[size_t(i)][size_t(j)] : (j - 6 == i ? 1 : 0);
    bool ok = true;
    for (int c = 0; c < 6 && ok; ++c) {
      int p = c;
      for (int i = c + 1; i < 6; ++i) if (std::fabs(A[i][c]) > std::fabs(A[p][c])) p = i;
      if (std::fabs(A[p][c]) < 1e-12) { ok = false; break; }
      for (int j = 0; j < 12; ++j) std::swap(A[c][j], A[p][j]);
      const double dd = A[c][c];
      for (int j = 0; j < 12; ++j) A[c][j] /= dd;
      for (int i = 0; i < 6; ++i)
        if (i != c) {
          const double m = A[i][c];
          for (int j = 0; j < 12; ++j) A[i][j] -= m * A[c][j];
        }
    }
    if (ok)
      for (int a = 0; a < 3; ++a)
        for (int j = 0; j < 6; ++j) {
          double sum = 0;
          for (int k = 0; k < 6; ++k) sum += r.e[size_t(a)][size_t(k)] * A[k][j + 6];
          r.d[size_t(a)][size_t(j)] = sum * 1e3;
        }
    else
      r.notes.push_back("the elastic constants are singular: no d coefficients");
  }
  stop("done", 1);
  return r;
}

}  // namespace caps
