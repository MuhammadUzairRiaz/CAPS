// CAPS's machine-learned potential through libtorch (see caps/mlpot.hpp). The compiled MACE model takes one dict of
// tensors: positions, the one-hot elements (node_attrs, in the model's own element order), the directed edges within its
// cutoff (both directions, no self-edges), their periodic shifts (zero here: one molecule), batch / ptr for a single
// graph, and the training labels and weights, which inference does not use but the compiled graph expects. The input
// was checked against MACE's own Python calculator: the same energy and forces to 1e-12 eV.
#include "caps/mlpot.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#ifdef CAPS_WITH_TORCH
#include <ATen/Parallel.h>
#include <torch/script.h>

#include <map>
#include <mutex>
#endif

namespace caps::ml {

#ifndef CAPS_WITH_TORCH

bool available() { return false; }
void set_threads(int) {}
struct Potential::Impl {};
Potential::Potential(const std::string&) {
  throw std::runtime_error("this CAPS was built without libtorch: configure with -DCAPS_WITH_TORCH=ON "
                           "-DTorch_DIR=<torch>/share/cmake/Torch to use the machine-learned potential");
}
Potential::~Potential() = default;
Result Potential::evaluate(const std::vector<int>&, const std::vector<std::array<double, 3>>&) const { return {}; }
double Potential::cutoff() const { return 0; }
const std::vector<int>& Potential::elements() const { static const std::vector<int> none; return none; }
bool Potential::covers(const std::vector<int>&) const { return false; }

#else

bool available() { return true; }

void set_threads(int n) {
  if (n > 0) {
    at::set_num_threads(n);
    try { at::set_num_interop_threads(n); } catch (const std::exception&) {}   // only settable before its first use
  }
}

struct Potential::Impl {
  mutable torch::jit::Module module;
  mutable std::mutex lock;    // one evaluation at a time: a TorchScript module is not reentrant
  double r_max = 0;
  std::vector<int> z;         // the model's elements, in its node_attrs order
  std::map<int, int> index;   // atomic number -> column of node_attrs
};

Potential::Potential(const std::string& path) : p_(std::make_unique<Impl>()) {
  try {
    p_->module = torch::jit::load(path, torch::kCPU);
  } catch (const std::exception& e) {
    throw std::runtime_error("cannot load the machine-learned potential " + path + ": " + e.what());
  }
  p_->module.eval();
  if (!p_->module.hasattr("r_max") || !p_->module.hasattr("atomic_numbers"))
    throw std::runtime_error(path + " is not a compiled MACE model (no r_max / atomic_numbers)");
  p_->r_max = p_->module.attr("r_max").toTensor().item<double>();
  auto an = p_->module.attr("atomic_numbers").toTensor().to(torch::kLong);
  for (int64_t k = 0; k < an.size(0); ++k) {
    p_->z.push_back(static_cast<int>(an[k].item<int64_t>()));
    p_->index[p_->z.back()] = static_cast<int>(k);
  }
}

Potential::~Potential() = default;
double Potential::cutoff() const { return p_->r_max; }
const std::vector<int>& Potential::elements() const { return p_->z; }

bool Potential::covers(const std::vector<int>& z) const {
  for (int e : z)
    if (!p_->index.count(e)) return false;
  return true;
}

Result Potential::evaluate(const std::vector<int>& z, const std::vector<std::array<double, 3>>& pos) const {
  const int64_t n = static_cast<int64_t>(z.size());
  if (n == 0 || pos.size() != z.size()) throw std::runtime_error("ml::Potential: no atoms, or positions and elements differ in number");
  if (!covers(z)) throw std::runtime_error("the machine-learned potential does not cover an element of this structure");
  const auto f64 = torch::TensorOptions().dtype(torch::kFloat64);
  const auto i64 = torch::TensorOptions().dtype(torch::kLong);
  const int64_t ne = static_cast<int64_t>(p_->z.size());
  auto positions = torch::empty({n, 3}, f64);
  auto node_attrs = torch::zeros({n, ne}, f64);
  auto pa = positions.accessor<double, 2>();
  auto na = node_attrs.accessor<double, 2>();
  for (int64_t i = 0; i < n; ++i) {
    for (int k = 0; k < 3; ++k) pa[i][k] = pos[i][k];
    na[i][p_->index.at(z[i])] = 1.0;
  }
  // directed edges within the cutoff (all pairs: fragments are small)
  std::vector<int64_t> snd, rcv;
  const double rc2 = p_->r_max * p_->r_max;
  for (int64_t i = 0; i < n; ++i)
    for (int64_t j = 0; j < n; ++j) {
      if (i == j) continue;
      double d2 = 0;
      for (int k = 0; k < 3; ++k) d2 += (pos[i][k] - pos[j][k]) * (pos[i][k] - pos[j][k]);
      if (d2 < rc2) {
        snd.push_back(i);
        rcv.push_back(j);
      }
    }
  const int64_t e = static_cast<int64_t>(snd.size());
  auto edge_index = torch::empty({2, e}, i64);
  auto ea = edge_index.accessor<int64_t, 2>();
  for (int64_t k = 0; k < e; ++k) {
    ea[0][k] = snd[k];
    ea[1][k] = rcv[k];
  }
  c10::Dict<std::string, torch::Tensor> d;
  d.insert("positions", positions);
  d.insert("node_attrs", node_attrs);
  d.insert("edge_index", edge_index);
  d.insert("shifts", torch::zeros({e, 3}, f64));
  d.insert("unit_shifts", torch::zeros({e, 3}, f64));
  d.insert("cell", torch::zeros({3, 3}, f64));
  d.insert("rcell", torch::zeros({3, 3}, f64));
  d.insert("volume", torch::zeros({1}, f64));
  d.insert("pbc", torch::zeros({1, 3}, torch::TensorOptions().dtype(torch::kBool)));
  d.insert("batch", torch::zeros({n}, i64));
  d.insert("ptr", torch::tensor({int64_t(0), n}, i64));
  d.insert("head", torch::zeros({1}, i64));
  // labels and weights: unused in inference, expected by the compiled graph
  d.insert("energy", torch::zeros({1}, f64));
  d.insert("forces", torch::zeros({n, 3}, f64));
  d.insert("stress", torch::zeros({1, 3, 3}, f64));
  d.insert("virials", torch::zeros({1, 3, 3}, f64));
  d.insert("dipole", torch::zeros({1, 3}, f64));
  d.insert("polarizability", torch::zeros({1, 3, 3}, f64));
  d.insert("charges", torch::zeros({n}, f64));
  d.insert("density_coefficients", torch::zeros({n, 1}, f64));
  d.insert("total_charge", torch::zeros({1}, f64));
  d.insert("total_spin", torch::ones({1}, f64));
  d.insert("elec_temp", torch::zeros({1}, f64));
  d.insert("fermi_level", torch::zeros({1}, f64));
  d.insert("external_field", torch::zeros({1, 3}, f64));
  for (const char* w : {"weight", "energy_weight", "forces_weight", "stress_weight", "virials_weight", "charges_weight"})
    d.insert(w, torch::ones({1}, f64));
  d.insert("dipole_weight", torch::ones({1, 3}, f64));
  d.insert("polarizability_weight", torch::ones({1, 3, 3}, f64));

  std::lock_guard<std::mutex> g(p_->lock);
  std::vector<torch::jit::IValue> args{d, false, true};   // training = false, compute_force = true
  auto out = p_->module.forward(args).toGenericDict();
  Result r;
  r.energy = out.at("energy").toTensor().item<double>();
  auto f = out.at("forces").toTensor().to(torch::kFloat64).contiguous();
  auto fa = f.accessor<double, 2>();
  r.forces.resize(n);
  for (int64_t i = 0; i < n; ++i) r.forces[i] = {fa[i][0], fa[i][1], fa[i][2]};
  return r;
}

#endif

// ------------------------------------------------------------------ constrained relaxation (both builds: it only calls
// evaluate, which a build without libtorch never reaches because no Potential can be made)
namespace {
using V = std::array<double, 3>;
V sub(const V& a, const V& b) { return {a[0] - b[0], a[1] - b[1], a[2] - b[2]}; }
double dot(const V& a, const V& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
V cross(const V& a, const V& b) { return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]}; }
V scale(const V& a, double s) { return {a[0] * s, a[1] * s, a[2] * s}; }
double norm(const V& a) { return std::sqrt(dot(a, a)); }
constexpr double kPi = 3.14159265358979323846;

// value and gradient (d value / d x of each atom of the Hold) of a bond, angle or dihedral
double coordinate(const Hold& h, const std::vector<V>& x, std::array<V, 4>* g) {
  const auto& q = h.idx;
  if (h.n == 2) {
    V d = sub(x[q[0]], x[q[1]]);
    double r = norm(d);
    if (g) { (*g)[0] = scale(d, 1 / r); (*g)[1] = scale(d, -1 / r); }
    return r;
  }
  if (h.n == 3) {
    V u = sub(x[q[0]], x[q[1]]), v = sub(x[q[2]], x[q[1]]);
    double lu = norm(u), lv = norm(v), c = dot(u, v) / (lu * lv);
    c = std::max(-1.0, std::min(1.0, c));
    double th = std::acos(c), s = std::max(1e-8, std::sin(th));
    if (g) {
      V gi = scale(sub(scale(v, 1 / (lu * lv)), scale(u, c / (lu * lu))), -1 / s);
      V gk = scale(sub(scale(u, 1 / (lu * lv)), scale(v, c / (lv * lv))), -1 / s);
      (*g)[0] = gi; (*g)[2] = gk; (*g)[1] = scale({gi[0] + gk[0], gi[1] + gk[1], gi[2] + gk[2]}, -1);
    }
    return th;
  }
  // dihedral i-j-k-l (Blondel & Karplus, J. Comput. Chem. 17, 1132 (1996))
  V b1 = sub(x[q[1]], x[q[0]]), b2 = sub(x[q[2]], x[q[1]]), b3 = sub(x[q[3]], x[q[2]]);
  V m = cross(b1, b2), n = cross(b2, b3);
  double lb2 = norm(b2), mm = std::max(1e-12, dot(m, m)), nn = std::max(1e-12, dot(n, n));
  double phi = std::atan2(lb2 * dot(b1, n), dot(m, n));
  if (g) {
    V gi = scale(m, -lb2 / mm), gl = scale(n, lb2 / nn);
    double p = dot(b1, b2) / (lb2 * lb2), r = dot(b3, b2) / (lb2 * lb2);
    // the middle atoms: d(phi)/dx_j = -(1 + p) g_i + r g_l, d(phi)/dx_k = p g_i - (1 + r) g_l (the four sum to zero);
    // checked against finite differences to 3e-10
    V gj = sub(scale(gl, r), scale(gi, 1 + p)), gk = sub(scale(gi, p), scale(gl, 1 + r));
    (*g)[0] = gi; (*g)[1] = gj; (*g)[2] = gk; (*g)[3] = gl;
  }
  return phi;
}

double wrap(double d) {
  while (d > kPi) d -= 2 * kPi;
  while (d <= -kPi) d += 2 * kPi;
  return d;
}
}  // namespace

double hold_value(const Hold& h, const std::vector<std::array<double, 3>>& pos) { return coordinate(h, pos, nullptr); }

RelaxResult relax(const Potential& pot, const std::vector<int>& z, std::vector<std::array<double, 3>> x,
                  const std::vector<Hold>& holds, double fmax, int max_steps, double /*k*/) {
  // Held coordinates are kept exactly, as ASE's FixInternals does: the force is projected onto the space where no held
  // coordinate changes, and after each step the geometry is brought back onto the constraints by Newton iterations.
  // The minimiser is L-BFGS on the projected force with ASE's settings (memory 100, initial curvature 70 eV/A^2, no
  // atom moves more than 0.2 A in a step).
  const size_t n = x.size(), m = holds.size(), dim = 3 * n;
  auto flat = [&](const std::vector<V>& a) { std::vector<double> f(dim); for (size_t i = 0; i < n; ++i) for (int d = 0; d < 3; ++d) f[3 * i + d] = a[i][d]; return f; };
  // the gradient of each held coordinate as a 3N vector
  auto grads = [&](const std::vector<V>& at) {
    std::vector<std::vector<double>> G(m, std::vector<double>(dim, 0.0));
    for (size_t c = 0; c < m; ++c) {
      std::array<V, 4> g{};
      coordinate(holds[c], at, &g);
      for (int a = 0; a < holds[c].n; ++a)
        for (int d = 0; d < 3; ++d) G[c][3 * holds[c].idx[a] + d] += g[a][d];
    }
    return G;
  };
  // solve the m x m system (Gram matrix) by Gauss elimination with partial pivoting; m is small
  auto solve = [](std::vector<std::vector<double>> A, std::vector<double> b) {
    const size_t k = b.size();
    for (size_t c = 0; c < k; ++c) {
      size_t p = c;
      for (size_t r = c + 1; r < k; ++r) if (std::fabs(A[r][c]) > std::fabs(A[p][c])) p = r;
      std::swap(A[c], A[p]); std::swap(b[c], b[p]);
      if (std::fabs(A[c][c]) < 1e-14) continue;
      for (size_t r = c + 1; r < k; ++r) {
        const double f = A[r][c] / A[c][c];
        for (size_t q = c; q < k; ++q) A[r][q] -= f * A[c][q];
        b[r] -= f * b[c];
      }
    }
    std::vector<double> x(k, 0.0);
    for (size_t c = k; c-- > 0;) {
      double s = b[c];
      for (size_t q = c + 1; q < k; ++q) s -= A[c][q] * x[q];
      x[c] = std::fabs(A[c][c]) < 1e-14 ? 0.0 : s / A[c][c];
    }
    return x;
  };
  auto project = [&](std::vector<double>& v, const std::vector<std::vector<double>>& G) {   // v minus its part along G
    if (!m) return;
    std::vector<std::vector<double>> A(m, std::vector<double>(m));
    std::vector<double> b(m);
    for (size_t a = 0; a < m; ++a) {
      for (size_t c = 0; c < m; ++c) { double s = 0; for (size_t q = 0; q < dim; ++q) s += G[a][q] * G[c][q]; A[a][c] = s; }
      double s = 0; for (size_t q = 0; q < dim; ++q) s += G[a][q] * v[q]; b[a] = s;
    }
    const auto l = solve(A, b);
    for (size_t c = 0; c < m; ++c) for (size_t q = 0; q < dim; ++q) v[q] -= l[c] * G[c][q];
  };
  auto enforce = [&](std::vector<V>& at) {   // Newton steps back onto every held value
    for (int it = 0; it < 50 && m; ++it) {
      std::vector<double> dev(m);
      double worst = 0;
      for (size_t c = 0; c < m; ++c) {
        const double v = coordinate(holds[c], at, nullptr);
        dev[c] = holds[c].n == 4 ? wrap(v - holds[c].target) : v - holds[c].target;
        worst = std::max(worst, std::fabs(dev[c]));
      }
      if (worst < 1e-10) break;
      const auto G = grads(at);
      std::vector<std::vector<double>> A(m, std::vector<double>(m));
      for (size_t a = 0; a < m; ++a)
        for (size_t c = 0; c < m; ++c) { double s = 0; for (size_t q = 0; q < dim; ++q) s += G[a][q] * G[c][q]; A[a][c] = s; }
      const auto mu = solve(A, dev);
      for (size_t c = 0; c < m; ++c)
        for (size_t i = 0; i < n; ++i) for (int d = 0; d < 3; ++d) at[i][d] -= mu[c] * G[c][3 * i + d];
    }
  };
  enforce(x);
  auto force = [&](const std::vector<V>& at, double* e) {   // the projected force, flat
    Result r = pot.evaluate(z, at);
    *e = r.energy;
    std::vector<double> f = flat(r.forces);
    project(f, grads(at));
    return f;
  };
  // L-BFGS (ASE's: two-loop recursion, H0 = 1 / 70, steps scaled so no atom moves more than 0.2 A)
  const size_t mem = 100;
  const double h0 = 1.0 / 70.0, maxstep = 0.2;
  std::vector<std::vector<double>> S, Y;
  std::vector<double> rho;
  double e = 0;
  std::vector<double> f = force(x, &e), x0, f0;
  RelaxResult out;
  for (int step = 0; step <= max_steps; ++step) {
    double worst = 0;
    for (size_t i = 0; i < n; ++i) worst = std::max(worst, std::sqrt(f[3 * i] * f[3 * i] + f[3 * i + 1] * f[3 * i + 1] + f[3 * i + 2] * f[3 * i + 2]));
    out.steps = step;
    if (worst < fmax) { out.converged = true; break; }
    if (step == max_steps) break;
    const std::vector<double> xf = flat(x);
    if (!x0.empty()) {
      std::vector<double> s(dim), y(dim);
      double sy = 0;
      for (size_t q = 0; q < dim; ++q) { s[q] = xf[q] - x0[q]; y[q] = f0[q] - f[q]; sy += s[q] * y[q]; }
      if (sy > 1e-12) {
        S.push_back(s); Y.push_back(y); rho.push_back(1.0 / sy);
        if (S.size() > mem) { S.erase(S.begin()); Y.erase(Y.begin()); rho.erase(rho.begin()); }
      }
    }
    // direction = -H g with g = -f
    std::vector<double> q(dim);
    for (size_t t = 0; t < dim; ++t) q[t] = -f[t];
    std::vector<double> al(S.size());
    for (size_t k = S.size(); k-- > 0;) {
      double a = 0; for (size_t t = 0; t < dim; ++t) a += S[k][t] * q[t];
      al[k] = rho[k] * a;
      for (size_t t = 0; t < dim; ++t) q[t] -= al[k] * Y[k][t];
    }
    for (size_t t = 0; t < dim; ++t) q[t] *= h0;
    for (size_t k = 0; k < S.size(); ++k) {
      double b = 0; for (size_t t = 0; t < dim; ++t) b += Y[k][t] * q[t];
      b *= rho[k];
      for (size_t t = 0; t < dim; ++t) q[t] += S[k][t] * (al[k] - b);
    }
    std::vector<double> dr(dim);
    for (size_t t = 0; t < dim; ++t) dr[t] = -q[t];
    project(dr, grads(x));
    double big = 0;
    for (size_t i = 0; i < n; ++i) big = std::max(big, std::sqrt(dr[3 * i] * dr[3 * i] + dr[3 * i + 1] * dr[3 * i + 1] + dr[3 * i + 2] * dr[3 * i + 2]));
    const double sc = big > maxstep ? maxstep / big : 1.0;
    x0 = xf;
    f0 = f;
    for (size_t i = 0; i < n; ++i) for (int d = 0; d < 3; ++d) x[i][d] += sc * dr[3 * i + d];
    enforce(x);
    f = force(x, &e);
  }
  out.energy = e;
  out.pos = x;
  for (const auto& h : holds) out.held.push_back(coordinate(h, x, nullptr));
  return out;
}

}  // namespace caps::ml
