// Cell relaxation with chosen lattice parameters held (cell_relax.hpp).
#include "caps/cell_relax.hpp"
#include "caps/cell_fit.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "caps/spacegroup.hpp"

namespace caps {

namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr double kDeg = kPi / 180.0;
constexpr double kAtmToKcalA3 = 101325.0 * 1e-30 * 6.02214076e23 / 4184.0;   // 1 atm in kcal/mol/Å³
constexpr double kKcalA3ToGPa = 4184.0 / 6.02214076e23 / 1e-30 / 1e9;        // 1 kcal/mol/Å³ in GPa

using M3 = std::array<std::array<double, 3>, 3>;
M3 cols(const Cell& c) { return {{{c.a[0], c.b[0], c.c[0]}, {c.a[1], c.b[1], c.c[1]}, {c.a[2], c.b[2], c.c[2]}}}; }
M3 mul(const M3& A, const M3& B) {
  M3 C{};
  for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) for (int k = 0; k < 3; ++k) C[i][j] += A[i][k] * B[k][j];
  return C;
}
Vec3 mul(const M3& A, const Vec3& v) { return {A[0][0] * v[0] + A[0][1] * v[1] + A[0][2] * v[2], A[1][0] * v[0] + A[1][1] * v[1] + A[1][2] * v[2], A[2][0] * v[0] + A[2][1] * v[1] + A[2][2] * v[2]}; }
Vec3 mul_t(const M3& A, const Vec3& v) { return {A[0][0] * v[0] + A[1][0] * v[1] + A[2][0] * v[2], A[0][1] * v[0] + A[1][1] * v[1] + A[2][1] * v[2], A[0][2] * v[0] + A[1][2] * v[1] + A[2][2] * v[2]}; }
double det(const M3& A) {
  return A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) - A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0]) + A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);
}
M3 inverse(const M3& A) {
  const double d = det(A);
  if (std::fabs(d) < 1e-14) throw std::runtime_error("a singular cell");
  M3 I{};
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) {
      const int i1 = (j + 1) % 3, i2 = (j + 2) % 3, j1 = (i + 1) % 3, j2 = (i + 2) % 3;
      I[i][j] = (A[i1][j1] * A[i2][j2] - A[i1][j2] * A[i2][j1]) / d;
    }
  return I;
}
Cell cell_of(const M3& h, const Vec3& origin) {
  Cell c;
  c.origin = origin;
  c.a = {h[0][0], h[1][0], h[2][0]}, c.b = {h[0][1], h[1][1], h[2][1]}, c.c = {h[0][2], h[1][2], h[2][2]};
  return c;
}

}  // namespace

std::array<double, 6> lattice_parameters(const Cell& c) {
  auto ang = [](const Vec3& x, const Vec3& y) { return std::acos(std::clamp(dot(x, y) / (norm(x) * norm(y)), -1.0, 1.0)) / kDeg; };
  return {norm(c.a), norm(c.b), norm(c.c), ang(c.b, c.c), ang(c.a, c.c), ang(c.a, c.b)};
}

Cell standard_cell(const std::array<double, 6>& p) {
  const double ca = std::cos(p[3] * kDeg), cb = std::cos(p[4] * kDeg), cg = std::cos(p[5] * kDeg), sg = std::sin(p[5] * kDeg);
  if (p[0] <= 0 || p[1] <= 0 || p[2] <= 0 || sg <= 1e-9) throw std::invalid_argument("not a cell: lengths must be positive, γ between 0° and 180°");
  const double cy = (ca - cb * cg) / sg, z2 = 1 - cb * cb - cy * cy;
  if (z2 <= 1e-12) throw std::invalid_argument("not a cell: the angles do not close");
  Cell c;
  c.a = {p[0], 0, 0};
  c.b = {p[1] * cg, p[1] * sg, 0};
  c.c = {p[2] * cb, p[2] * cy, p[2] * std::sqrt(z2)};
  return c;
}

CellRelaxReport relax_cell(System& s, const std::function<ForceField(const System&)>& field_for, const EnergyOptions& eo,
                           const CellRelaxOptions& o) {
  if (!s.cell.valid()) throw std::invalid_argument("cell relaxation needs a periodic cell");
  if (s.atoms.empty()) throw std::invalid_argument("no atoms");
  if (!field_for) throw std::invalid_argument("cell relaxation needs a force field");
  const size_t N = s.atoms.size();
  CellRelaxReport rep;
  const std::array<double, 6> p0 = lattice_parameters(s.cell);
  rep.before = p0;
  const M3 h0 = cols(s.cell);
  const M3 R0 = mul(h0, inverse(cols(standard_cell(p0))));   // the orientation kept
  const Vec3 origin = s.cell.origin;
  const M3 h0i = inverse(h0);
  std::vector<Vec3> frac(N);
  for (size_t i = 0; i < N; ++i) frac[i] = mul(h0i, s.atoms[i].pos - origin);
  // copies so the cell seen is at least twice the cut-off wide (3 % margin for shrinking)
  const double V0 = std::fabs(det(h0));
  const double width[3] = {V0 / norm(cross(s.cell.b, s.cell.c)), V0 / norm(cross(s.cell.c, s.cell.a)), V0 / norm(cross(s.cell.a, s.cell.b))};
  std::array<int, 3> sc{1, 1, 1};
  for (int d = 0; d < 3; ++d) sc[size_t(d)] = std::max(1, int(std::ceil(2.06 * eo.cutoff / width[d] - 1e-9)));
  const auto fit = copies_to_fit(s);   // a short cut-off still needs the bonds read at their own images
  for (size_t d = 0; d < 3; ++d) sc[d] = std::max(sc[d], fit[d]);
  if (size_t(sc[0]) * size_t(sc[1]) * size_t(sc[2]) * N > 15000)
    throw std::invalid_argument("the copies this cut-off needs exceed 15000 atoms: use a larger cell or a shorter cut-off");
  rep.copies = sc;
  const size_t ng = size_t(sc[0]) * size_t(sc[1]) * size_t(sc[2]);
  const System sup = ng > 1 ? supercell(s, sc[0], sc[1], sc[2]) : s;
  const ForceField ff = field_for(sup);
  if (ff.type_index.size() != sup.atoms.size()) throw std::invalid_argument("the force field does not describe the copies of the cell");
  std::vector<std::array<int, 3>> off(ng);
  for (size_t g = 0; g < ng; ++g) {
    const Vec3 f = mul(h0i, sup.atoms[g * N].pos - s.atoms[0].pos);
    for (int k = 0; k < 3; ++k) off[g][size_t(k)] = int(std::lround(f[size_t(k)]));
  }
  if (ng > 1) rep.notes.push_back("evaluated as " + std::to_string(sc[0]) + "×" + std::to_string(sc[1]) + "×" + std::to_string(sc[2]) + " copies (the cut-off is wider than half the cell)");
  // The search runs on a continuous enthalpy. With the tail correction on, the Evaluator truncates the pair energy at
  // the cut-off, so the energy steps by V(rc) whenever a pair crosses it, while the forces and the virial (which give
  // the search direction) never see the step. A line search on that energy stalls as soon as a descent step pushes a
  // pair out across rc. So the search evaluates with the pair energy shifted to zero at the cut-off (tail off, which
  // leaves forces and virial unchanged) and adds the tail term itself: E_t = e_t / V with e_t fixed by the composition
  // and the cut-off. The caller's options give the reported energies, enthalpy and stress.
  Evaluator ev_report(ff, eo);
  EnergyOptions eo_search = eo;
  eo_search.tail = false;
  Evaluator ev(ff, eo_search);
  double tail_et = 0;   // E_t × V of the copies, from the first evaluation with the caller's options
  const double P = o.pressure * kAtmToKcalA3;
  const double L = std::cbrt(V0);

  // the variables: movable atoms' fractional coordinates × L, then the free lattice parameters (angles in rad × L)
  std::vector<size_t> mv;
  for (size_t i = 0; i < N; ++i)
    if (o.move_atoms && !(i < o.fixed.size() && o.fixed[i])) mv.push_back(i);
  std::vector<int> free;
  for (int k = 0; k < 6; ++k) if (!o.hold[size_t(k)]) free.push_back(k);
  const size_t na = 3 * mv.size(), nv = na + free.size();
  auto params_of = [&](const std::vector<double>& z) {
    std::array<double, 6> p = p0;
    for (size_t q = 0; q < free.size(); ++q) {
      const int k = free[q];
      p[size_t(k)] = k < 3 ? z[na + q] : z[na + q] / L / kDeg;
    }
    return p;
  };
  auto h_of = [&](const std::array<double, 6>& p) { return mul(R0, cols(standard_cell(p))); };
  struct Eval {
    double E = 0, H = 0, V = 0;
    std::array<double, 6> W{};
    std::vector<Vec3> g;   // ∂E/∂x per atom of the cell
  };
  auto evaluate = [&](const std::vector<double>& z, M3& h, std::vector<Vec3>& fr, bool report) {
    const auto p = params_of(z);
    h = h_of(p);
    fr = frac;
    for (size_t m = 0; m < mv.size(); ++m)
      for (int k = 0; k < 3; ++k) fr[mv[m]][size_t(k)] = z[3 * m + size_t(k)] / L;
    std::vector<double> x(3 * N * ng), f;
    for (size_t g = 0; g < ng; ++g) {
      const Vec3 shift = mul(h, Vec3{double(off[g][0]), double(off[g][1]), double(off[g][2])});
      for (size_t i = 0; i < N; ++i) {
        const Vec3 r = origin + mul(h, fr[i]) + shift;
        for (int k = 0; k < 3; ++k) x[3 * (g * N + i) + size_t(k)] = r[size_t(k)];
      }
    }
    M3 hs = h;
    for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) hs[i][j] *= double(sc[size_t(j)]);
    const Cell sup_cell = cell_of(hs, origin);
    const EnergyTerms t = (report ? ev_report : ev).compute(x, sup_cell, f);
    Eval e;
    for (int k = 0; k < 6; ++k) e.W[size_t(k)] = t.w[k] / double(ng);
    // the tail term: the virial has the fluid's tail pressure; the gradient of the energy needs E_t's own volume
    // derivative, a virial of E_t on each diagonal component
    if (report) {
      tail_et = t.tail_energy * sup_cell.volume();
      e.E = t.total() / double(ng);
      for (int k = 0; k < 3; ++k) e.W[size_t(k)] += (t.tail_energy - t.tail_virial / 3) / double(ng);
    } else {
      const double et = tail_et / sup_cell.volume();
      e.E = (t.total() + et) / double(ng);
      for (int k = 0; k < 3; ++k) e.W[size_t(k)] += et / double(ng);
    }
    e.V = std::fabs(det(h));
    e.H = e.E + P * e.V;
    e.g.assign(N, Vec3{0, 0, 0});
    for (size_t g = 0; g < ng; ++g)
      for (size_t i = 0; i < N; ++i)
        for (int k = 0; k < 3; ++k) e.g[i][size_t(k)] -= f[3 * (g * N + i) + size_t(k)] / double(ng);
    if (!std::isfinite(e.E)) throw std::runtime_error("the energy is not finite");
    return e;
  };
  // ∂H/∂p for every lattice parameter: (P V I − W) h⁻ᵀ contracted with ∂h/∂p (∂h/∂p by central differences of the
  // geometry alone)
  auto cell_gradient = [&](const Eval& e, const M3& h, const std::array<double, 6>& p) {
    const M3 Wm = {{{e.W[0], e.W[3], e.W[4]}, {e.W[3], e.W[1], e.W[5]}, {e.W[4], e.W[5], e.W[2]}}};
    const M3 hi = inverse(h);
    // e.V = |det h| and ∂|det h|/∂h = |det h| h⁻ᵀ for either handedness, so no sign factor here
    M3 D{};
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) {
        double v = 0;
        for (int k = 0; k < 3; ++k) v += ((i == k ? P * e.V : 0.0) - Wm[i][k]) * hi[j][k];
        D[i][j] = v;
      }
    std::array<double, 6> gp{};
    for (int k = 0; k < 6; ++k) {
      const double step = k < 3 ? 1e-6 * p[size_t(k)] : 1e-6;   // Å, or degrees
      std::array<double, 6> pp = p, pm = p;
      pp[size_t(k)] += step, pm[size_t(k)] -= step;
      const M3 hp = h_of(pp), hm = h_of(pm);
      double v = 0;
      for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) v += D[i][j] * (hp[i][j] - hm[i][j]) / (2 * step);
      gp[size_t(k)] = k < 3 ? v : v / kDeg;   // per Å, or per radian
    }
    return gp;
  };
  M3 h;
  std::vector<Vec3> fr;
  auto value = [&](const std::vector<double>& z, std::vector<double>& gz, double& gmax, Eval* out, bool report = false) {
    const Eval e = evaluate(z, h, fr, report);
    gz.assign(nv, 0.0);
    gmax = 0;
    for (size_t m = 0; m < mv.size(); ++m) {
      const Vec3 gs = mul_t(h, e.g[mv[m]]);   // ∂E/∂s
      for (int k = 0; k < 3; ++k) gz[3 * m + size_t(k)] = gs[size_t(k)] / L;
      gmax = std::max(gmax, norm(e.g[mv[m]]));
    }
    if (!free.empty()) {
      const auto gp = cell_gradient(e, h, params_of(z));
      for (size_t q = 0; q < free.size(); ++q) {
        const int k = free[q];
        gz[na + q] = k < 3 ? gp[size_t(k)] : gp[size_t(k)] / L;
        gmax = std::max(gmax, std::fabs(gz[na + q]));
      }
    }
    if (out) *out = e;
    return e.H;
  };
  std::vector<double> z(nv, 0.0), g;
  for (size_t m = 0; m < mv.size(); ++m) for (int k = 0; k < 3; ++k) z[3 * m + size_t(k)] = frac[mv[m]][size_t(k)] * L;
  for (size_t q = 0; q < free.size(); ++q) z[na + q] = free[q] < 3 ? p0[size_t(free[q])] : p0[size_t(free[q])] * kDeg * L;
  double gmax = 0;
  Eval e0;
  value(z, g, gmax, &e0, true);   // the caller's energy, and the tail term's e_t
  rep.energy_before = e0.E, rep.volume_before = e0.V;
  double H = value(z, g, gmax, nullptr);
  std::vector<std::vector<double>> S, Y;
  std::vector<double> rho;
  int it = 0;
  bool stalled = false;
  bool done = gmax < o.ftol;
  for (; it < o.max_steps && !done; ++it) {
    if (o.progress && !o.progress(it, H, gmax)) throw std::runtime_error("cell relaxation cancelled");
    std::vector<double> q = g, al(S.size());
    for (size_t k = S.size(); k-- > 0;) {
      double d = 0;
      for (size_t v = 0; v < nv; ++v) d += S[k][v] * q[v];
      al[k] = rho[k] * d;
      for (size_t v = 0; v < nv; ++v) q[v] -= al[k] * Y[k][v];
    }
    double gamma = 1.0 / 100.0;   // Å²/(kcal/mol) until a curvature is known
    if (!S.empty()) {
      double sy = 0, yy = 0;
      for (size_t v = 0; v < nv; ++v) sy += S.back()[v] * Y.back()[v], yy += Y.back()[v] * Y.back()[v];
      if (yy > 0) gamma = sy / yy;
    }
    for (auto& v : q) v *= gamma;
    for (size_t k = 0; k < S.size(); ++k) {
      double d = 0;
      for (size_t v = 0; v < nv; ++v) d += Y[k][v] * q[v];
      const double be = rho[k] * d;
      for (size_t v = 0; v < nv; ++v) q[v] += S[k][v] * (al[k] - be);
    }
    std::vector<double> p(nv);
    double slope = 0, pmax = 0;
    for (size_t v = 0; v < nv; ++v) p[v] = -q[v], slope += p[v] * g[v], pmax = std::max(pmax, std::fabs(p[v]));
    if (slope >= 0) {
      S.clear(), Y.clear(), rho.clear();
      slope = 0, pmax = 0;
      for (size_t v = 0; v < nv; ++v) p[v] = -g[v] / 100.0, slope += p[v] * g[v], pmax = std::max(pmax, std::fabs(p[v]));
    }
    const bool fresh = S.empty();   // a steepest-descent direction, no curvature history behind it
    double step = pmax > 0.2 ? 0.2 / pmax : 1.0;
    bool ok = false;
    std::vector<double> zn(nv), gn;
    double Hn = H, gnmax = gmax;
    for (int ls = 0; ls < 40; ++ls) {
      for (size_t v = 0; v < nv; ++v) zn[v] = z[v] + step * p[v];
      try {
        Hn = value(zn, gn, gnmax, nullptr);
        if (Hn <= H + 1e-4 * step * slope) { ok = true; break; }
      } catch (const std::invalid_argument&) {   // the step made no cell: shorter
      }
      step *= 0.5;
    }
    // a failed L-BFGS direction gets one retry along steepest descent; a failed steepest descent ends the run
    if (!ok) { if (fresh) { stalled = true; break; } S.clear(), Y.clear(), rho.clear(); continue; }
    std::vector<double> sv(nv), yv(nv);
    double sy = 0;
    for (size_t v = 0; v < nv; ++v) sv[v] = zn[v] - z[v], yv[v] = gn[v] - g[v], sy += sv[v] * yv[v];
    if (sy > 1e-12) {
      S.push_back(sv), Y.push_back(yv), rho.push_back(1.0 / sy);
      if (S.size() > 12) S.erase(S.begin()), Y.erase(Y.begin()), rho.erase(rho.begin());
    }
    z = zn, g = gn, H = Hn, gmax = gnmax;
    done = gmax < o.ftol;
  }
  Eval ef;
  value(z, g, gmax, &ef, true);
  rep.steps = it, rep.converged = done;
  for (size_t i = 0; i < N; ++i) s.atoms[i].pos = origin + mul(h, fr[i]);
  s.cell = cell_of(h, origin);
  rep.after = lattice_parameters(s.cell);
  const auto gp = cell_gradient(ef, h, rep.after);
  rep.gradient = gp;
  rep.energy_after = ef.E, rep.enthalpy_after = ef.H, rep.volume_after = ef.V;
  rep.stress = {-ef.W[0] / ef.V * kKcalA3ToGPa, -ef.W[1] / ef.V * kKcalA3ToGPa, -ef.W[2] / ef.V * kKcalA3ToGPa,
                -ef.W[5] / ef.V * kKcalA3ToGPa, -ef.W[4] / ef.V * kKcalA3ToGPa, -ef.W[3] / ef.V * kKcalA3ToGPa};
  rep.pressure_after = (ef.W[0] + ef.W[1] + ef.W[2]) / (3 * ef.V) / kAtmToKcalA3;
  char nb[200];
  if (!done && stalled)
    std::snprintf(nb, sizeof nb, "stopped after %d steps: no lower enthalpy along the search at this precision (largest gradient %.4f, asked %.4f)", it, gmax, o.ftol);
  else if (!done)
    std::snprintf(nb, sizeof nb, "not converged in %d steps (largest gradient %.4f, asked %.4f)", it, gmax, o.ftol);
  if (!done) rep.notes.push_back(nb);
  return rep;
}

}  // namespace caps
