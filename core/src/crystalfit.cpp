// CAPS force-field fitting to crystals (see crystalfit.hpp).
#include "caps/crystalfit.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>

#include "caps/ffdef.hpp"
#include "caps/mechanics.hpp"
#include "caps/spacegroup.hpp"

namespace caps {

namespace {
constexpr double kPi = 3.14159265358979323846;
constexpr double kGPa = 4184.0 / 6.02214076e23 / 1e-30 / 1e9;   // GPa per kcal/mol/Å³
const double kNaN = std::numeric_limits<double>::quiet_NaN();

using M3 = std::array<std::array<double, 3>, 3>;

M3 identity() { return {{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}}; }
Vec3 mul(const M3& m, const Vec3& v) {
  return {m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2], m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
          m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2]};
}
Vec3 mul_t(const M3& m, const Vec3& v) {   // mᵀ v
  return {m[0][0] * v[0] + m[1][0] * v[1] + m[2][0] * v[2], m[0][1] * v[0] + m[1][1] * v[1] + m[2][1] * v[2],
          m[0][2] * v[0] + m[1][2] * v[1] + m[2][2] * v[2]};
}
M3 inverse(const M3& m) {
  const double d = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
                   m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
  M3 r;
  r[0][0] = (m[1][1] * m[2][2] - m[1][2] * m[2][1]) / d;
  r[0][1] = (m[0][2] * m[2][1] - m[0][1] * m[2][2]) / d;
  r[0][2] = (m[0][1] * m[1][2] - m[0][2] * m[1][1]) / d;
  r[1][0] = (m[1][2] * m[2][0] - m[1][0] * m[2][2]) / d;
  r[1][1] = (m[0][0] * m[2][2] - m[0][2] * m[2][0]) / d;
  r[1][2] = (m[0][2] * m[1][0] - m[0][0] * m[1][2]) / d;
  r[2][0] = (m[1][0] * m[2][1] - m[1][1] * m[2][0]) / d;
  r[2][1] = (m[0][1] * m[2][0] - m[0][0] * m[2][1]) / d;
  r[2][2] = (m[0][0] * m[1][1] - m[0][1] * m[1][0]) / d;
  return r;
}
// symmetric deformation I + ε from tensor components e = xx yy zz yz xz xy (not engineering shear)
M3 deformation(const std::array<double, 6>& e) {
  M3 F = identity();
  F[0][0] += e[0], F[1][1] += e[1], F[2][2] += e[2];
  F[1][2] += e[3], F[2][1] += e[3];
  F[0][2] += e[4], F[2][0] += e[4];
  F[0][1] += e[5], F[1][0] += e[5];
  return F;
}

// the state of one crystal: cell vectors and the cell's atom positions (Cartesian)
struct State {
  Vec3 a{}, b{}, c{};
  std::vector<Vec3> X;
};
Cell cell_of(const State& s) {
  Cell c;
  c.a = s.a, c.b = s.b, c.c = s.c;
  return c;
}
State deformed(const State& s, const M3& F) {
  State t;
  t.a = mul(F, s.a), t.b = mul(F, s.b), t.c = mul(F, s.c);
  t.X.resize(s.X.size());
  for (size_t i = 0; i < s.X.size(); ++i) t.X[i] = mul(F, s.X[i]);
  return t;
}

struct Evaluation {
  double E = 0;                  // per cell
  std::vector<Vec3> g;           // ∂E/∂X per cell atom
  std::array<double, 6> W{};     // per cell, xx yy zz xy xz yz
  double V = 0;
};

// One crystal in its periodic copies.
struct Crystal {
  size_t N = 0;
  std::array<int, 3> sc{1, 1, 1};
  std::vector<std::array<int, 3>> off;   // per copy g: the lattice offset (copy g of atom a is supercell atom g·N + a)
  ForceField base;                       // the supercell's force field as given
  std::vector<char> movable_all, shells_only;
  std::vector<char> shell;
  double Q0 = 0;                         // total charge of the supercell as given

  Evaluation evaluate(Evaluator& ev, const State& s) const {
    const size_t ng = off.size();
    std::vector<double> x(3 * N * ng), f;
    for (size_t g = 0; g < ng; ++g) {
      const Vec3 o = s.a * double(off[g][0]) + s.b * double(off[g][1]) + s.c * double(off[g][2]);
      for (size_t i = 0; i < N; ++i)
        for (int k = 0; k < 3; ++k) x[3 * (g * N + i) + size_t(k)] = s.X[i][size_t(k)] + o[size_t(k)];
    }
    Cell sup;
    sup.a = s.a * double(sc[0]), sup.b = s.b * double(sc[1]), sup.c = s.c * double(sc[2]);
    const EnergyTerms t = ev.compute(x, sup, f);
    Evaluation r;
    r.E = t.total() / double(ng);
    r.g.assign(N, Vec3{0, 0, 0});
    for (size_t g = 0; g < ng; ++g)
      for (size_t i = 0; i < N; ++i)
        for (int k = 0; k < 3; ++k) r.g[i][size_t(k)] -= f[3 * (g * N + i) + size_t(k)] / double(ng);
    for (int k = 0; k < 6; ++k) r.W[size_t(k)] = t.w[k] / double(ng);
    r.V = cell_of(s).volume();
    if (!std::isfinite(r.E)) throw std::runtime_error("the energy is not finite");
    return r;
  }
};

// stress (GPa, Voigt xx yy zz yz xz xy, tensile positive)
std::array<double, 6> stress(const Evaluation& e) {
  return {-e.W[0] / e.V * kGPa, -e.W[1] / e.V * kGPa, -e.W[2] / e.V * kGPa, -e.W[5] / e.V * kGPa, -e.W[4] / e.V * kGPa, -e.W[3] / e.V * kGPa};
}

// L-BFGS over the movable atoms' reference positions and (with cell) the six symmetric strains, scaled by a cell length.
// Returns the steps taken; s is updated.
int relax_state(const Crystal& c, Evaluator& ev, State& s, bool cell, const std::vector<char>& movable, double ftol, int max_it,
                bool* converged) {
  const size_t N = c.N;
  std::vector<size_t> mv;
  for (size_t i = 0; i < N; ++i) if (movable[i]) mv.push_back(i);
  const size_t nv = 3 * mv.size() + (cell ? 6 : 0);
  if (converged) *converged = true;
  if (nv == 0) return 0;
  const State s0 = s;
  const double L = std::cbrt(cell_of(s0).volume());
  // z → state: X = F u for the movable atoms (u their reference positions), affine for the others
  auto state_of = [&](const std::vector<double>& z, M3& F) {
    std::array<double, 6> e{};
    if (cell)
      for (int k = 0; k < 6; ++k) e[size_t(k)] = z[3 * mv.size() + size_t(k)] / L;
    F = deformation(e);
    State t = deformed(s0, F);
    for (size_t m = 0; m < mv.size(); ++m) t.X[mv[m]] = mul(F, Vec3{z[3 * m], z[3 * m + 1], z[3 * m + 2]});
    return t;
  };
  auto value = [&](const std::vector<double>& z, std::vector<double>& gz, double& gmax) {
    M3 F;
    const State t = state_of(z, F);
    const Evaluation e = c.evaluate(ev, t);
    gz.assign(nv, 0.0);
    gmax = 0;
    for (size_t m = 0; m < mv.size(); ++m) {
      const Vec3 gu = mul_t(F, e.g[mv[m]]);
      for (int k = 0; k < 3; ++k) gz[3 * m + size_t(k)] = gu[size_t(k)];
      gmax = std::max(gmax, norm(e.g[mv[m]]));
    }
    if (cell) {
      // ∂E/∂F = −W F⁻ᵀ, W symmetric; the symmetric strain's off-diagonal component takes both entries
      const M3 W = {{{e.W[0], e.W[3], e.W[4]}, {e.W[3], e.W[1], e.W[5]}, {e.W[4], e.W[5], e.W[2]}}};
      const M3 Fi = inverse(F);
      M3 D{};
      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
          for (int k = 0; k < 3; ++k) D[size_t(i)][size_t(j)] -= W[size_t(i)][size_t(k)] * Fi[size_t(j)][size_t(k)];
      const double ge[6] = {D[0][0], D[1][1], D[2][2], D[1][2] + D[2][1], D[0][2] + D[2][0], D[0][1] + D[1][0]};
      for (int k = 0; k < 6; ++k) {
        gz[3 * mv.size() + size_t(k)] = ge[k] / L;
        gmax = std::max(gmax, std::fabs(ge[k]) / L);
      }
    }
    return e.E;
  };
  std::vector<double> z(nv, 0.0), g, zn, gn;
  for (size_t m = 0; m < mv.size(); ++m)
    for (int k = 0; k < 3; ++k) z[3 * m + size_t(k)] = s0.X[mv[m]][size_t(k)];
  double gmax = 0;
  double E = value(z, g, gmax);
  const int mem = 12;
  std::vector<std::vector<double>> S, Y;
  std::vector<double> rho;
  int it = 0;
  bool done = gmax < ftol;
  for (; it < max_it && !done; ++it) {
    // two-loop recursion
    std::vector<double> q = g, al(S.size());
    for (size_t k = S.size(); k-- > 0;) {
      double d = 0;
      for (size_t v = 0; v < nv; ++v) d += S[k][v] * q[v];
      al[k] = rho[k] * d;
      for (size_t v = 0; v < nv; ++v) q[v] -= al[k] * Y[k][v];
    }
    double gamma = 1.0 / 100.0;   // Å²/(kcal/mol): a stiffness of 100 kcal/mol/Å² until a curvature is known
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
    if (slope >= 0) {   // not a descent direction: steepest descent
      S.clear(), Y.clear(), rho.clear();
      slope = 0, pmax = 0;
      for (size_t v = 0; v < nv; ++v) p[v] = -g[v] / 100.0, slope += p[v] * g[v], pmax = std::max(pmax, std::fabs(p[v]));
    }
    double step = pmax > 0.2 ? 0.2 / pmax : 1.0;
    bool ok = false;
    double En = E, gnmax = gmax;
    for (int ls = 0; ls < 30; ++ls) {
      zn = z;
      for (size_t v = 0; v < nv; ++v) zn[v] += step * p[v];
      try {
        En = value(zn, gn, gnmax);
      } catch (const std::exception&) {
        step *= 0.5;
        continue;
      }
      // Armijo, or (at the energy's rounding level) a smaller gradient
      const double tiny = 1e-12 * std::max(1.0, std::fabs(E));
      if (En <= E + 1e-4 * step * slope || (std::fabs(En - E) < tiny && gnmax < gmax)) {
        ok = true;
        break;
      }
      step *= 0.5;
    }
    if (!ok) {
      if (S.empty()) break;   // already steepest: nothing more to gain at this precision
      S.clear(), Y.clear(), rho.clear();
      continue;
    }
    std::vector<double> sv(nv), yv(nv);
    double sy = 0;
    for (size_t v = 0; v < nv; ++v) sv[v] = zn[v] - z[v], yv[v] = gn[v] - g[v], sy += sv[v] * yv[v];
    if (sy > 1e-16) {
      S.push_back(sv), Y.push_back(yv), rho.push_back(1.0 / sy);
      if (S.size() > size_t(mem)) S.erase(S.begin()), Y.erase(Y.begin()), rho.erase(rho.begin());
    }
    z.swap(zn), g.swap(gn), E = En, gmax = gnmax;
    done = gmax < ftol;
  }
  M3 F;
  s = state_of(z, F);
  if (converged) *converged = done;
  if (cell && std::fabs(cell_of(s).volume() / cell_of(s0).volume() - 1) > 0.5)
    throw std::runtime_error("the cell changed volume by more than half: no stable minimum near the structure");
  return it;
}

// ---- parameters on a force field
int type_index(const ForceField& F, const std::string& n) {
  for (size_t t = 0; t < F.type_names.size(); ++t)
    if (F.type_names[t] == n) return int(t);
  return -1;
}
bool same_pair(const ForceField& F, uint32_t i, uint32_t j, const std::vector<std::string>& t) {
  const std::string& a = F.atom_type[i];
  const std::string& b = F.atom_type[j];
  return (a == t[0] && b == t[1]) || (a == t[1] && b == t[0]);
}
int slot(const std::vector<std::string>& names, const std::string& n) {
  for (size_t k = 0; k < names.size(); ++k)
    if (names[k] == n) return int(k);
  return -1;
}
std::string join(const std::vector<std::string>& v) {
  std::string s;
  for (const auto& x : v) s += (s.empty() ? "" : " ") + x;
  return s;
}

// The term a parameter belongs to in F: its style, its parameters' names and current values. With set_k ≥ 0 it first
// writes v into parameter set_k. False when F has no such term.
struct Access {
  std::string style;
  std::vector<std::string> names;
  std::vector<double> values;
};
bool access(ForceField& F, const FitParameter& p, Access& out, int set_k = -1, double v = 0) {
  if (p.kind == "pair") {
    if (p.types.size() != 2) throw std::invalid_argument(p.label() + ": a pair parameter names two types");
    const int a = type_index(F, p.types[0]), b = type_index(F, p.types[1]);
    if (a < 0 || b < 0) return false;
    const std::pair<int, int> key{std::min(a, b), std::max(a, b)};
    if (auto it = F.pair_func.find(key); it != F.pair_func.end()) {
      PairFunc& f = it->second;
      double* field[3] = {&f.a, &f.b, &f.c};
      if (f.form == 1) out.style = "buck", out.names = {"A", "rho", "C"};
      else if (f.form == 2) out.style = "morse", out.names = {"D", "alpha", "r0"};
      else if (f.form == kPairRydberg) out.style = "rydberg", out.names = {"A", "rho", "C"};
      else throw std::invalid_argument(p.label() + ": this pair's form cannot be fitted");
      if (set_k >= 0 && set_k < 3) *field[set_k] = v;
      out.values = {f.a, f.b, f.c};
      return true;
    }
    PairType pt = mixed_pair(F, a, b);
    if (auto it = F.pair_override.find(key); it != F.pair_override.end()) pt = it->second;
    out.style = "lj/cut", out.names = {"epsilon", "sigma"};
    if (set_k == 0 || set_k == 1) {
      (set_k == 0 ? pt.eps : pt.sigma) = v;
      F.pair_override[key] = pt;
    }
    out.values = {pt.eps, pt.sigma};
    return true;
  }
  if (p.kind == "lj") {
    if (p.types.size() != 1) throw std::invalid_argument(p.label() + ": names one type");
    const int a = type_index(F, p.types[0]);
    if (a < 0) return false;
    PairType& pt = F.lj[size_t(a)];
    out.style = "", out.names = {"epsilon", "sigma"};
    if (set_k == 0 || set_k == 1) (set_k == 0 ? pt.eps : pt.sigma) = v;
    out.values = {pt.eps, pt.sigma};
    return true;
  }
  if (p.kind == "bond") {
    if (p.types.size() != 2) throw std::invalid_argument(p.label() + ": a bond parameter names two types");
    bool found = false;
    for (auto& t : F.bonds)
      if (same_pair(F, t.i, t.j, p.types)) {
        found = true, out.style = "harmonic", out.names = {"k", "r0"};
        if (set_k == 0) t.k = v;
        if (set_k == 1) t.r0 = v;
        out.values = {t.k, t.r0};
      }
    for (auto& t : F.bonds2)
      if (same_pair(F, t.i, t.j, p.types)) {
        found = true, out.style = "class2", out.names = {"r0", "k2", "k3", "k4"};
        double* field[4] = {&t.r0, &t.k2, &t.k3, &t.k4};
        if (set_k >= 0 && set_k < 4) *field[set_k] = v;
        out.values = {t.r0, t.k2, t.k3, t.k4};
      }
    for (auto& t : F.bonds_x)
      if (t.form == 1 && same_pair(F, t.i, t.j, p.types)) {
        found = true, out.style = "morse", out.names = {"D", "alpha", "r0"};
        double* field[3] = {&t.a, &t.b, &t.c};
        if (set_k >= 0 && set_k < 3) *field[set_k] = v;
        out.values = {t.a, t.b, t.c};
      }
    return found;
  }
  if (p.kind == "charge") {
    if (p.types.size() != 1) throw std::invalid_argument(p.label() + ": names one type");
    bool found = false;
    for (size_t i = 0; i < F.atom_type.size(); ++i)
      if (F.atom_type[i] == p.types[0]) {
        found = true;
        if (set_k == 0) F.charge[i] = v;
        out.values = {F.charge[i]};
      }
    out.style = "", out.names = {"q"};
    return found;
  }
  throw std::invalid_argument("parameter kind '" + p.kind + "': pair, lj, bond or charge");
}

// every parameter onto a copy of the crystal's force field
ForceField with_parameters(const Crystal& c, const std::vector<FitParameter>& ps, const std::vector<double>& v) {
  ForceField F = c.base;
  for (size_t k = 0; k < ps.size(); ++k) {
    Access a;
    if (!access(F, ps[k], a)) continue;
    access(F, ps[k], a, slot(a.names, ps[k].name), v[k]);
    if (ps[k].kind == "charge" && !ps[k].balance.empty()) {
      size_t nb = 0;
      double others = 0;
      for (size_t i = 0; i < F.atom_type.size(); ++i) {
        if (F.atom_type[i] == ps[k].balance) ++nb;
        else others += F.charge[i];
      }
      if (nb == 0) throw std::invalid_argument(ps[k].label() + ": no atom of the balancing type " + ps[k].balance);
      const double qb = (c.Q0 - others) / double(nb);
      for (size_t i = 0; i < F.atom_type.size(); ++i)
        if (F.atom_type[i] == ps[k].balance) F.charge[i] = qb;
    }
  }
  return F;
}

bool is_elastic_kind(const std::string& k) {
  return k == "bulk" || k == "shear" || (k.size() == 3 && k[0] == 'C' && k[1] >= '1' && k[1] <= '6' && k[2] >= '1' && k[2] <= '6');
}

// everything computed for one structure at one parameter set
struct Props {
  State st;                         // final state (relaxed, or the observed one with its shells relaxed)
  std::array<double, 6> cellp{};
  double V = 0, E = 0;              // per cell
  std::vector<Vec3> dX;             // Å, per atom (relaxed fitting): displacement from observed, drift removed
  std::vector<Vec3> force;          // kcal/mol/Å at the observed structure (conventional)
  std::array<double, 6> sigma{};    // GPa (conventional)
  bool elastic = false;
  std::array<std::array<double, 6>, 6> C{};
  double K = 0, G = 0;
  int steps = 0;
  bool converged = true;
};

Props properties(const Crystal& c, const State& observed, const State& start, const ForceField& F, const FitOptions& o, bool want_elastic) {
  Evaluator ev(F, o.energy);
  Props p;
  p.st = o.relaxed ? start : observed;
  if (o.relaxed) p.steps = relax_state(c, ev, p.st, true, c.movable_all, o.relax_ftol, o.relax_iterations, &p.converged);
  else p.steps = relax_state(c, ev, p.st, false, c.shells_only, o.relax_ftol, o.relax_iterations, &p.converged);
  const Evaluation e = c.evaluate(ev, p.st);
  p.E = e.E, p.V = e.V;
  p.cellp = cell_parameters(cell_of(p.st));
  if (o.relaxed) {
    // fractional displacement from observed, wrapped, its mean over the atoms (not shells) removed, in Å
    const M3 h = {{{p.st.a[0], p.st.b[0], p.st.c[0]}, {p.st.a[1], p.st.b[1], p.st.c[1]}, {p.st.a[2], p.st.b[2], p.st.c[2]}}};
    const M3 h0 = {{{observed.a[0], observed.b[0], observed.c[0]}, {observed.a[1], observed.b[1], observed.c[1]}, {observed.a[2], observed.b[2], observed.c[2]}}};
    const M3 hi = inverse(h), h0i = inverse(h0);
    std::vector<Vec3> df(c.N);
    Vec3 mean{0, 0, 0};
    size_t n = 0;
    for (size_t i = 0; i < c.N; ++i) {
      const Vec3 f1 = mul(hi, p.st.X[i]), f0 = mul(h0i, observed.X[i]);
      for (int k = 0; k < 3; ++k) {
        double d = f1[size_t(k)] - f0[size_t(k)];
        d -= std::round(d);
        df[i][size_t(k)] = d;
      }
      if (!c.shell[i]) mean = mean + df[i], ++n;
    }
    if (n) mean = mean * (1.0 / double(n));
    p.dX.assign(c.N, Vec3{0, 0, 0});
    for (size_t i = 0; i < c.N; ++i) {
      Vec3 d = df[i] - mean;
      for (auto& x : d) x -= std::round(x);
      p.dX[i] = mul(h, d);
    }
  } else {
    p.force.resize(c.N);
    for (size_t i = 0; i < c.N; ++i) p.force[i] = e.g[i] * -1.0;
    p.sigma = stress(e);
  }
  if (want_elastic) {
    // the atoms relaxed in the cell (already, in relaxed fitting), then ± each strain with the atoms relaxed again
    State ref = p.st;
    if (!o.relaxed) relax_state(c, ev, ref, false, c.movable_all, o.relax_ftol, o.relax_iterations, nullptr);
    const double d = o.elastic_strain;
    for (int J = 0; J < 6; ++J) {
      std::array<double, 6> sp[2];
      for (int sgn = 0; sgn < 2; ++sgn) {
        std::array<double, 6> e6{};
        e6[size_t(J)] = (sgn ? -d : d) * (J < 3 ? 1.0 : 0.5);
        State t = deformed(ref, deformation(e6));
        relax_state(c, ev, t, false, c.movable_all, o.relax_ftol, o.relax_iterations, nullptr);
        sp[sgn] = stress(c.evaluate(ev, t));
      }
      for (int I = 0; I < 6; ++I) p.C[size_t(I)][size_t(J)] = (sp[0][size_t(I)] - sp[1][size_t(I)]) / (2 * d);
    }
    for (int I = 0; I < 6; ++I)
      for (int J = I + 1; J < 6; ++J) p.C[size_t(I)][size_t(J)] = p.C[size_t(J)][size_t(I)] = 0.5 * (p.C[size_t(I)][size_t(J)] + p.C[size_t(J)][size_t(I)]);
    ElasticResult er;
    for (int I = 0; I < 6; ++I)
      for (int J = 0; J < 6; ++J) er.C[size_t(I)][size_t(J)] = p.C[size_t(I)][size_t(J)];
    isotropic_averages(er);
    p.K = er.K_hill, p.G = er.G_hill;
    p.elastic = true;
  }
  return p;
}

double target_value(const Props& p, const FitStructure& fs, const std::string& k) {
  static const std::map<std::string, int> cp = {{"a", 0}, {"b", 1}, {"c", 2}, {"alpha", 3}, {"beta", 4}, {"gamma", 5}};
  if (auto it = cp.find(k); it != cp.end()) return p.cellp[size_t(it->second)];
  if (k == "volume") return p.V;
  if (k == "energy") return p.E / fs.formula_units;
  if (k == "bulk") return p.K;
  if (k == "shear") return p.G;
  if (is_elastic_kind(k)) return p.C[size_t(k[1] - '1')][size_t(k[2] - '1')];
  throw std::invalid_argument("target '" + k + "': a b c alpha beta gamma volume energy C11…C66 bulk shear");
}
std::string unit_of(const std::string& k) {
  if (k == "a" || k == "b" || k == "c") return "Å";
  if (k == "alpha" || k == "beta" || k == "gamma") return "°";
  if (k == "volume") return "Å³";
  if (k == "energy") return "kcal/mol";
  return "GPa";
}

// small dense solve (Gauss–Jordan with pivoting); false when singular
bool solve(std::vector<double> A, std::vector<double>& b, size_t n) {
  for (size_t c = 0; c < n; ++c) {
    size_t p = c;
    for (size_t i = c + 1; i < n; ++i) if (std::fabs(A[i * n + c]) > std::fabs(A[p * n + c])) p = i;
    if (std::fabs(A[p * n + c]) < 1e-300) return false;
    if (p != c) {
      for (size_t j = 0; j < n; ++j) std::swap(A[c * n + j], A[p * n + j]);
      std::swap(b[c], b[p]);
    }
    const double d = A[c * n + c];
    for (size_t j = 0; j < n; ++j) A[c * n + j] /= d;
    b[c] /= d;
    for (size_t i = 0; i < n; ++i)
      if (i != c) {
        const double m = A[i * n + c];
        if (m == 0) continue;
        for (size_t j = 0; j < n; ++j) A[i * n + j] -= m * A[c * n + j];
        b[i] -= m * b[c];
      }
  }
  return true;
}
}  // namespace

std::string FitParameter::label() const { return kind + " " + join(types) + " " + name; }

std::array<double, 6> cell_parameters(const Cell& c) {
  auto ang = [](const Vec3& u, const Vec3& v) { return std::acos(std::clamp(dot(u, v) / (norm(u) * norm(v)), -1.0, 1.0)) * 180 / kPi; };
  return {norm(c.a), norm(c.b), norm(c.c), ang(c.b, c.c), ang(c.a, c.c), ang(c.a, c.b)};
}

namespace {
constexpr double kFitPmeSpacing = 0.5;
FitResult fit_impl(const FitOptions& o, FitResult r);
}  // namespace

FitResult fit_to_crystals(const FitOptions& o) {
  FitResult r;
  if (o.structures.empty()) throw std::invalid_argument("fit: no structures");
  // a fit compares nearby cells: the Ewald grid's discretisation error must not move the minimum (at 1 Å spacing it
  // shifts a rock-salt lattice constant by 5e-4 relative, at 0.5 Å by under 1e-6)
  FitOptions oo = o;
  if (oo.energy.electrostatics == EnergyOptions::Electrostatics::PME && oo.energy.coulomb && oo.energy.pme_spacing > kFitPmeSpacing) {
    char b[160];
    std::snprintf(b, sizeof b, "Ewald grid spacing %.2f Å → %.2f Å for the fit", oo.energy.pme_spacing, kFitPmeSpacing);
    r.notes.push_back(b);
    oo.energy.pme_spacing = kFitPmeSpacing;
  }
  return fit_impl(oo, std::move(r));
}

namespace {
FitResult fit_impl(const FitOptions& o, FitResult r) {
  const auto t0 = std::chrono::steady_clock::now();
  auto stop = [&](const std::string& w, double f) {
    if (o.progress && !o.progress(w, f)) throw std::runtime_error("fit stopped");
  };
  // ---- the crystals
  const size_t ns = o.structures.size();
  std::vector<Crystal> cr(ns);
  std::vector<State> observed(ns);
  std::vector<bool> want_el(ns, false);
  for (const auto& t : o.targets) {
    if (t.structure >= ns) throw std::invalid_argument("fit: a target names structure " + std::to_string(t.structure + 1) + " of " + std::to_string(ns));
    if (t.weight < 0) throw std::invalid_argument("fit: a negative weight");
    if (is_elastic_kind(t.kind)) want_el[t.structure] = true;
  }
  for (size_t si = 0; si < ns; ++si) {
    const FitStructure& fs = o.structures[si];
    const System& s = fs.cell;
    const std::string nm = fs.name.empty() ? "structure " + std::to_string(si + 1) : fs.name;
    if (!s.cell.valid()) throw std::invalid_argument("fit: " + nm + " needs a periodic cell");
    if (!fs.field_for) throw std::invalid_argument("fit: " + nm + " has no force field");
    if (s.atoms.empty()) throw std::invalid_argument("fit: " + nm + " has no atoms");
    if (fs.formula_units <= 0) throw std::invalid_argument("fit: formula units must be positive");
    Crystal& c = cr[si];
    c.N = s.atoms.size();
    const double V = s.cell.volume();
    const double width[3] = {V / norm(cross(s.cell.b, s.cell.c)), V / norm(cross(s.cell.c, s.cell.a)), V / norm(cross(s.cell.a, s.cell.b))};
    // relaxed fitting may shrink the cell a little: a 3 % margin over twice the cut-off
    std::array<int, 3> sc = fs.supercell;
    const bool automatic = sc[0] <= 0 || sc[1] <= 0 || sc[2] <= 0;
    for (int d = 0; d < 3; ++d)
      if (sc[size_t(d)] <= 0) sc[size_t(d)] = std::max(1, int(std::ceil(2.06 * o.energy.cutoff / width[d] - 1e-9)));
    if (automatic && size_t(sc[0]) * size_t(sc[1]) * size_t(sc[2]) * c.N > o.max_atoms)
      throw std::invalid_argument("fit: " + nm + " needs a supercell above " + std::to_string(o.max_atoms) + " atoms for this cut-off");
    c.sc = sc;
    const System sup = supercell(s, sc[0], sc[1], sc[2]);
    c.base = fs.field_for(sup);
    if (c.base.type_index.size() != sup.atoms.size()) throw std::invalid_argument("fit: the force field does not describe the supercell of " + nm);
    const size_t ng = sup.atoms.size() / c.N;
    const M3 h = {{{s.cell.a[0], s.cell.b[0], s.cell.c[0]}, {s.cell.a[1], s.cell.b[1], s.cell.c[1]}, {s.cell.a[2], s.cell.b[2], s.cell.c[2]}}};
    const M3 hi = inverse(h);
    c.off.resize(ng);
    for (size_t g = 0; g < ng; ++g) {
      const Vec3 f = mul(hi, sup.atoms[g * c.N].pos - s.atoms[0].pos);
      for (int k = 0; k < 3; ++k) c.off[g][size_t(k)] = int(std::lround(f[size_t(k)]));
    }
    c.shell.assign(c.N, 0);
    for (size_t i = 0; i < c.N && i < fs.shell.size(); ++i) c.shell[i] = fs.shell[i];
    c.movable_all.assign(c.N, 1);
    c.shells_only = c.shell;
    for (double q : c.base.charge) c.Q0 += q;
    observed[si].a = s.cell.a, observed[si].b = s.cell.b, observed[si].c = s.cell.c;
    for (const auto& a : s.atoms) observed[si].X.push_back(a.pos);
    FitStructureResult sr;
    sr.name = nm;
    sr.supercell = sc;
    sr.supercell_atoms = sup.atoms.size();
    sr.observed = cell_parameters(s.cell);
    sr.volume_observed = V;
    r.structures.push_back(sr);
  }
  // ---- the parameters: initial values from the first structure that has each
  const size_t np = o.parameters.size();
  std::vector<double> p(np), scale(np);
  for (size_t k = 0; k < np; ++k) {
    const FitParameter& fp = o.parameters[k];
    bool found = false;
    for (size_t si = 0; si < ns; ++si) {
      Access a;
      ForceField F = cr[si].base;
      if (!access(F, fp, a)) continue;
      const int s = slot(a.names, fp.name);
      if (s < 0) throw std::invalid_argument(fp.label() + ": this " + (a.style.empty() ? fp.kind : a.style) + " term has " + join(a.names));
      if (!found) {
        p[k] = a.values[size_t(s)];
        found = true;
      } else if (std::fabs(a.values[size_t(s)] - p[k]) > 1e-9 * std::max(1.0, std::fabs(p[k]))) {
        r.notes.push_back(fp.label() + " differs between structures at the start; all take the first one's value");
      }
    }
    if (!found) throw std::invalid_argument(fp.label() + ": no structure has this term");
    if (fp.lower > fp.upper) throw std::invalid_argument(fp.label() + ": lower bound above upper");
    p[k] = std::clamp(p[k], fp.lower, fp.upper);
    scale[k] = std::fabs(p[k]) > 0 ? std::fabs(p[k]) : 1.0;
    r.names.push_back(fp.label());
  }
  r.initial = p;
  // ---- residuals
  std::vector<State> start(observed);   // relaxed fitting starts each relaxation from the last accepted one
  int evals = 0;
  auto compute = [&](const std::vector<double>& pv, std::vector<Props>& props) {
    props.clear();
    for (size_t si = 0; si < ns; ++si) {
      const ForceField F = with_parameters(cr[si], o.parameters, pv);
      props.push_back(properties(cr[si], observed[si], start[si], F, o, want_el[si]));
    }
    ++evals;
  };
  auto residuals = [&](const std::vector<Props>& props) {
    std::vector<double> res;
    for (size_t si = 0; si < ns; ++si) {
      const Props& q = props[si];
      if (!o.structures[si].fit_structure) continue;
      if (o.relaxed) {
        const auto& ob = r.structures[si].observed;
        for (int k = 0; k < 3; ++k) res.push_back(std::sqrt(o.w_length) * (q.cellp[size_t(k)] - ob[size_t(k)]));
        for (int k = 3; k < 6; ++k) res.push_back(std::sqrt(o.w_angle) * (q.cellp[size_t(k)] - ob[size_t(k)]));
        for (size_t i = 0; i < cr[si].N; ++i)
          if (!cr[si].shell[i])
            for (int k = 0; k < 3; ++k) res.push_back(std::sqrt(o.w_coord) * q.dX[i][size_t(k)]);
      } else {
        for (size_t i = 0; i < cr[si].N; ++i)
          if (!cr[si].shell[i])
            for (int k = 0; k < 3; ++k) res.push_back(std::sqrt(o.w_force) * q.force[i][size_t(k)]);
        for (int k = 0; k < 6; ++k) res.push_back(std::sqrt(o.w_stress) * q.sigma[size_t(k)]);
      }
    }
    for (const auto& t : o.targets)
      res.push_back(std::sqrt(t.weight) * (target_value(props[t.structure], o.structures[t.structure], t.kind) - t.value));
    return res;
  };
  auto ss_of = [](const std::vector<double>& v) {
    double s = 0;
    for (double x : v) s += x * x;
    return s;
  };
  auto record = [&](const std::vector<Props>& props, bool initial) {
    for (size_t si = 0; si < ns; ++si) {
      const Props& q = props[si];
      FitStructureResult& sr = r.structures[si];
      (initial ? sr.initial : sr.final) = q.cellp;
      (initial ? sr.energy_initial : sr.energy_final) = q.E / o.structures[si].formula_units;
      double shift = 0, fmax = 0;
      size_t n = 0;
      for (size_t i = 0; i < cr[si].N; ++i) {
        if (cr[si].shell[i]) continue;
        ++n;
        if (o.relaxed) shift += dot(q.dX[i], q.dX[i]);
        else fmax = std::max(fmax, norm(q.force[i]));
      }
      (initial ? sr.shift_initial : sr.shift_final) = n ? std::sqrt(shift / double(n)) : 0;
      (initial ? sr.force_initial : sr.force_final) = fmax;
      (initial ? sr.stress_initial : sr.stress_final) = q.sigma;
      if (!initial) {
        sr.volume_final = q.V;
        sr.elastic = q.elastic, sr.C = q.C, sr.bulk = q.K, sr.shear = q.G;
        sr.cell = cell_of(q.st);
        sr.positions = q.st.X;
        sr.relax_steps = q.steps;
        if (!q.converged)
          r.notes.push_back(sr.name + ": the last relaxation stopped before every gradient was below " + std::to_string(o.relax_ftol));
      }
    }
  };
  stop("evaluating", 0);
  std::vector<Props> props;
  compute(p, props);
  std::vector<double> res = residuals(props);
  r.residuals = int(res.size());
  double ss = ss_of(res);
  r.ss_initial = ss;
  record(props, true);
  std::vector<double> obs_initial;
  for (const auto& t : o.targets) obs_initial.push_back(target_value(props[t.structure], o.structures[t.structure], t.kind));
  if (o.relaxed)
    for (size_t si = 0; si < ns; ++si) start[si] = props[si].st;
  // ---- Levenberg–Marquardt in the scaled parameters z = p / scale
  const size_t m = res.size();
  std::vector<double> J;   // m × np
  bool j_current = false;
  auto jacobian = [&](const std::vector<double>& pv, const std::vector<double>& r0) {
    J.assign(m * np, 0.0);
    for (size_t k = 0; k < np; ++k) {
      stop("derivatives", double(k) / double(np));
      std::vector<double> q = pv;
      double h = o.fd_step * std::max(std::fabs(pv[k]), scale[k]);
      if (q[k] + h > o.parameters[k].upper) h = -h;
      q[k] += h;
      std::vector<Props> pr;
      compute(q, pr);
      const std::vector<double> rk = residuals(pr);
      if (rk.size() != m) throw std::runtime_error("fit: the residual count changed");
      for (size_t i = 0; i < m; ++i) J[i * np + k] = (rk[i] - r0[i]) / h * scale[k];
    }
    j_current = true;
  };
  double lambda = 1e-3;
  int it = 0;
  r.stopped_by = np == 0 ? "no parameters" : "iterations";
  if (np > 0 && ss < 1e-24) r.stopped_by = "exact", r.converged = true;
  for (; np > 0 && it < o.iterations && !r.converged; ++it) {
    stop("iteration " + std::to_string(it + 1), double(it) / o.iterations);
    jacobian(p, res);
    std::vector<double> A(np * np, 0.0), gvec(np, 0.0);
    for (size_t i = 0; i < m; ++i)
      for (size_t a = 0; a < np; ++a) {
        gvec[a] -= J[i * np + a] * res[i];
        for (size_t b = 0; b < np; ++b) A[a * np + b] += J[i * np + a] * J[i * np + b];
      }
    bool accepted = false;
    for (int tries = 0; tries < 12 && !accepted && !r.converged; ++tries) {
      std::vector<double> M = A, d = gvec;
      for (size_t a = 0; a < np; ++a) M[a * np + a] += lambda * std::max(A[a * np + a], 1e-12);
      if (!solve(M, d, np)) { lambda *= 10; continue; }
      std::vector<double> q = p;
      double dmax = 0;
      for (size_t a = 0; a < np; ++a) {
        q[a] = std::clamp(p[a] + d[a] * scale[a], o.parameters[a].lower, o.parameters[a].upper);
        dmax = std::max(dmax, std::fabs(q[a] - p[a]) / scale[a]);
      }
      if (dmax < 1e-12) { r.stopped_by = "step", r.converged = true; break; }
      std::vector<Props> pr;
      double ssn = std::numeric_limits<double>::infinity();
      std::vector<double> rn;
      try {
        compute(q, pr);
        rn = residuals(pr);
        ssn = ss_of(rn);
      } catch (const std::runtime_error& e) {   // no stable structure at this trial: a smaller step
        if (std::string(e.what()) == "fit stopped") throw;
      }
      if (ssn < ss) {
        const double gain = (ss - ssn) / ss;
        p = q, res = rn, props = pr, ss = ssn;
        j_current = false;
        accepted = true;
        lambda = std::max(lambda / 10, 1e-12);
        if (o.relaxed)
          for (size_t si = 0; si < ns; ++si) start[si] = props[si].st;
        if (ss < 1e-24) r.stopped_by = "exact", r.converged = true;
        else if (gain < o.tolerance) r.stopped_by = "tolerance", r.converged = true;
      } else {
        lambda *= 10;
      }
    }
    r.history.push_back(ss);
    if (!accepted && !r.converged) {
      r.stopped_by = "no step lowers S";
      r.converged = true;
      ++it;
      break;
    }
  }
  r.iterations = it;
  // ---- standard errors: s² (JᵀJ)⁻¹
  r.error.assign(np, kNaN);
  if (np > 0 && m > np) {
    if (!j_current) jacobian(p, res);
    std::vector<double> A(np * np, 0.0);
    for (size_t i = 0; i < m; ++i)
      for (size_t a = 0; a < np; ++a)
        for (size_t b = 0; b < np; ++b) A[a * np + b] += J[i * np + a] * J[i * np + b];
    const double s2 = ss / double(m - np);
    for (size_t k = 0; k < np; ++k) {
      std::vector<double> e(np, 0.0);
      e[k] = 1;
      if (solve(A, e, np) && e[k] >= 0) r.error[k] = std::sqrt(s2 * e[k]) * scale[k];
    }
  } else if (np > 0) {
    r.notes.push_back("no standard errors: as many parameters as observations, or more");
  }
  r.final = p;
  r.ss_final = ss;
  record(props, false);
  for (size_t ti = 0; ti < o.targets.size(); ++ti) {
    const auto& t = o.targets[ti];
    FitObservation ob;
    ob.structure = t.structure, ob.what = t.kind, ob.unit = unit_of(t.kind), ob.target = t.value, ob.weight = t.weight;
    ob.initial = obs_initial[ti];
    ob.final = target_value(props[t.structure], o.structures[t.structure], t.kind);
    r.observations.push_back(ob);
  }
  // ---- the fitted terms in full (a balanced charge brings its balancing type)
  std::set<std::string> seen;
  for (size_t k = 0; k < np; ++k) {
    const FitParameter& fp = o.parameters[k];
    std::vector<FitParameter> also{fp};
    if (fp.kind == "charge" && !fp.balance.empty()) {
      FitParameter b = fp;
      b.types = {fp.balance}, b.balance.clear();
      also.push_back(b);
    }
    for (const auto& q : also) {
      const std::string key = q.kind + " " + join(q.types);
      if (seen.count(key)) continue;
      for (size_t si = 0; si < ns; ++si) {
        ForceField F = with_parameters(cr[si], o.parameters, p);
        Access a;
        if (!access(F, q, a)) continue;
        seen.insert(key);
        r.terms.push_back({q.kind, q.types, a.style, a.values});
        break;
      }
    }
  }
  r.evaluations = evals;
  r.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  stop("done", 1);
  return r;
}
}  // namespace

void fit_options_from_json(const Json& j, FitOptions& o) {
  auto strs = [](const Json& a) {
    std::vector<std::string> v;
    if (a.is_array())
      for (const auto& x : a.items()) v.push_back(x.str());
    return v;
  };
  if (j.has("parameters") && j["parameters"].is_array())
    for (const auto& p : j["parameters"].items()) {
      FitParameter fp;
      fp.kind = p.text("kind"), fp.name = p.text("name"), fp.balance = p.text("balance");
      if (p.has("types")) fp.types = strs(p["types"]);
      if (p.has("lower")) fp.lower = p["lower"].number();
      if (p.has("upper")) fp.upper = p["upper"].number();
      o.parameters.push_back(fp);
    }
  if (j.has("targets") && j["targets"].is_array())
    for (const auto& t : j["targets"].items()) {
      FitTarget ft;
      ft.structure = size_t(t.num("structure", 0));
      ft.kind = t.text("kind");
      ft.value = t.num("value", 0);
      ft.weight = t.num("weight", 1);
      const std::string u = t.text("unit", "");
      if (ft.kind == "energy" && !u.empty()) {
        if (u == "eV") ft.value *= 23.060547830619;   // kcal/mol per eV (CODATA 2018)
        else if (u == "kJ/mol") ft.value /= 4.184;
        else if (u != "kcal/mol") throw std::invalid_argument("energy unit '" + u + "': kcal/mol, kJ/mol or eV");
      }
      o.targets.push_back(ft);
    }
  if (j.has("relaxed")) o.relaxed = j["relaxed"].boolean();
  o.iterations = int(j.num("iterations", o.iterations));
  o.tolerance = j.num("tolerance", o.tolerance);
  o.fd_step = j.num("fd_step", o.fd_step);
  o.relax_ftol = j.num("relax_ftol", o.relax_ftol);
  o.elastic_strain = j.num("elastic_strain", o.elastic_strain);
  if (j.has("weights") && j["weights"].is_object()) {
    const Json& w = j["weights"];
    o.w_length = w.num("length", o.w_length), o.w_angle = w.num("angle", o.w_angle), o.w_coord = w.num("coord", o.w_coord);
    o.w_force = w.num("force", o.w_force), o.w_stress = w.num("stress", o.w_stress);
  }
}

Json fit_result_json(const FitResult& r) {
  auto arr = [](const auto& v) {
    Json a = Json::array();
    for (double x : v) a.push_back(Json(std::isfinite(x) ? x : 0.0));
    return a;
  };
  Json j = Json::object();
  Json ps = Json::array();
  for (size_t k = 0; k < r.names.size(); ++k) {
    Json p = Json::object();
    p["name"] = r.names[k], p["initial"] = r.initial[k], p["final"] = r.final[k];
    p["error"] = std::isfinite(r.error[k]) ? Json(r.error[k]) : Json(nullptr);
    ps.push_back(p);
  }
  j["parameters"] = ps;
  Json ob = Json::array();
  for (const auto& o : r.observations) {
    Json x = Json::object();
    x["structure"] = double(o.structure), x["what"] = o.what, x["unit"] = o.unit, x["target"] = o.target, x["initial"] = o.initial, x["final"] = o.final, x["weight"] = o.weight;
    ob.push_back(x);
  }
  j["observations"] = ob;
  Json st = Json::array();
  for (const auto& s : r.structures) {
    Json x = Json::object();
    x["name"] = s.name, x["supercell"] = arr(std::vector<double>{double(s.supercell[0]), double(s.supercell[1]), double(s.supercell[2])});
    x["supercell_atoms"] = double(s.supercell_atoms);
    x["observed"] = arr(s.observed), x["initial"] = arr(s.initial), x["final"] = arr(s.final);
    x["volume_observed"] = s.volume_observed, x["volume_final"] = s.volume_final;
    x["energy_initial"] = s.energy_initial, x["energy_final"] = s.energy_final;
    x["shift_initial"] = s.shift_initial, x["shift_final"] = s.shift_final, x["force_initial"] = s.force_initial, x["force_final"] = s.force_final;
    x["stress_initial"] = arr(s.stress_initial), x["stress_final"] = arr(s.stress_final);
    if (s.elastic) {
      Json C = Json::array();
      for (const auto& row : s.C) C.push_back(arr(row));
      x["C"] = C, x["bulk"] = s.bulk, x["shear"] = s.shear;
    }
    st.push_back(x);
  }
  j["structures"] = st;
  Json ts = Json::array();
  for (const auto& t : r.terms) {
    Json x = Json::object(), ty = Json::array();
    for (const auto& n : t.types) ty.push_back(n);
    x["kind"] = t.kind, x["types"] = ty, x["style"] = t.style, x["values"] = arr(t.values);
    ts.push_back(x);
  }
  j["terms"] = ts;
  j["ss_initial"] = r.ss_initial, j["ss_final"] = r.ss_final, j["history"] = arr(r.history);
  j["residuals"] = double(r.residuals), j["iterations"] = double(r.iterations), j["evaluations"] = double(r.evaluations);
  j["converged"] = r.converged, j["stopped_by"] = r.stopped_by, j["seconds"] = r.seconds;
  Json nt = Json::array();
  for (const auto& n : r.notes) nt.push_back(n);
  j["notes"] = nt;
  return j;
}

Json fittable_terms(const ForceField& ff0) {
  ForceField F = ff0;
  Json out = Json::array();
  std::set<int> present(F.type_index.begin(), F.type_index.end());
  auto add = [&](const FitParameter& p) {
    Access a;
    if (!access(F, p, a)) return;
    Json x = Json::object(), ty = Json::array(), nm = Json::array(), va = Json::array();
    for (const auto& t : p.types) ty.push_back(t);
    for (const auto& n : a.names) nm.push_back(n);
    for (double v : a.values) va.push_back(v);
    x["kind"] = p.kind, x["types"] = ty, x["style"] = a.style, x["names"] = nm, x["values"] = va;
    out.push_back(x);
  };
  for (int a : present)
    for (int b : present)
      if (a <= b) {
        const std::pair<int, int> key{a, b};
        const bool explicit_pair = F.pair_func.count(key) || F.pair_override.count(key);
        if (!explicit_pair) continue;
        if (auto it = F.pair_func.find(key); it != F.pair_func.end() && it->second.form != 1 && it->second.form != 2 && it->second.form != kPairRydberg) continue;
        add({"pair", {F.type_names[size_t(a)], F.type_names[size_t(b)]}, ""});
      }
  for (int a : present)
    if (F.lj[size_t(a)].eps != 0) add({"lj", {F.type_names[size_t(a)]}, ""});
  std::set<std::pair<std::string, std::string>> bonds;
  auto bond = [&](uint32_t i, uint32_t j) {
    std::string x = F.atom_type[i], y = F.atom_type[j];
    if (y < x) std::swap(x, y);
    bonds.insert({x, y});
  };
  for (const auto& t : F.bonds) bond(t.i, t.j);
  for (const auto& t : F.bonds2) bond(t.i, t.j);
  for (const auto& t : F.bonds_x) if (t.form == 1) bond(t.i, t.j);
  for (const auto& [x, y] : bonds) add({"bond", {x, y}, ""});
  std::set<std::string> charged;
  for (size_t i = 0; i < F.atom_type.size(); ++i)
    if (F.charge[i] != 0) charged.insert(F.atom_type[i]);
  for (const auto& t : charged) add({"charge", {t}, ""});
  return out;
}

std::vector<std::string> apply_fitted_terms(FFDef& def, const std::vector<FittedTerm>& terms) {
  std::vector<std::string> notes;
  auto equiv = [&](const std::string& t, const char* kind) {
    for (const auto& ty : def.types)
      if (ty.name == t) {
        auto it = ty.equiv.find(kind);
        return it == ty.equiv.end() ? t : it->second;
      }
    return t;
  };
  auto exact = [](const FFRule& r, const std::vector<std::string>& n) {
    if (r.match.size() != n.size()) return false;
    if (n.size() == 1) return r.match[0] == n[0];
    return (r.match[0] == n[0] && r.match[1] == n[1]) || (r.match[0] == n[1] && r.match[1] == n[0]);
  };
  for (const auto& t : terms) {
    if (t.kind == "charge") {
      bool done = false;
      for (auto& ty : def.types)
        if (ty.name == t.types[0]) ty.charge = t.values[0], done = true;
      notes.push_back(done ? "charge of " + t.types[0] + " set" : "no type " + t.types[0] + " to carry its charge");
      continue;
    }
    const bool bond = t.kind == "bond";
    std::vector<std::string> n;
    for (const auto& x : t.types) n.push_back(equiv(x, bond ? "bond" : "vdw"));
    auto& rules = bond ? def.bonds : def.pairs;
    FFRule* hit = nullptr;
    for (auto& rr : rules)
      if (exact(rr, n)) hit = &rr;   // the last exact rule is the one assigned
    const std::string style = t.kind == "lj" || t.style == "lj/cut" ? "" : t.style;
    if (hit) {
      hit->params = t.values;
      if (!style.empty()) hit->style = style;
      hit->comment = hit->comment.empty() ? "refitted to crystal data by CAPS"
                                          : "refitted to crystal data by CAPS (before: " + hit->comment + ")";
      notes.push_back(t.kind + " " + join(n) + ": rule " + hit->name + " updated");
    } else {
      FFRule rr;
      rr.name = join(n) + " (fitted)";
      rr.match = n;
      rr.style = style;
      rr.params = t.values;
      rr.comment = "fitted to crystal data by CAPS";
      rules.push_back(rr);
      notes.push_back(t.kind + " " + join(n) + ": rule appended");
    }
  }
  return notes;
}

}  // namespace caps
