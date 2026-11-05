// Molecular point groups and symmetrisation (pointgroup.hpp).
#include "caps/pointgroup.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "caps/elements.hpp"

namespace caps {

namespace {

constexpr double kPi = 3.14159265358979323846;
using M3 = std::array<std::array<double, 3>, 3>;

M3 identity() { return {{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}}; }
M3 mul(const M3& A, const M3& B) {
  M3 C{};
  for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) for (int k = 0; k < 3; ++k) C[i][j] += A[i][k] * B[k][j];
  return C;
}
Vec3 apply(const M3& A, const Vec3& v) { return {A[0][0] * v[0] + A[0][1] * v[1] + A[0][2] * v[2], A[1][0] * v[0] + A[1][1] * v[1] + A[1][2] * v[2], A[2][0] * v[0] + A[2][1] * v[1] + A[2][2] * v[2]}; }
Vec3 apply_t(const M3& A, const Vec3& v) { return {A[0][0] * v[0] + A[1][0] * v[1] + A[2][0] * v[2], A[0][1] * v[0] + A[1][1] * v[1] + A[2][1] * v[2], A[0][2] * v[0] + A[1][2] * v[1] + A[2][2] * v[2]}; }
double det(const M3& A) {
  return A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) - A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0]) + A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);
}
M3 rotation(const Vec3& u, double t) {
  const double c = std::cos(t), s = std::sin(t), C = 1 - c;
  return {{{c + u[0] * u[0] * C, u[0] * u[1] * C - u[2] * s, u[0] * u[2] * C + u[1] * s},
           {u[1] * u[0] * C + u[2] * s, c + u[1] * u[1] * C, u[1] * u[2] * C - u[0] * s},
           {u[2] * u[0] * C - u[1] * s, u[2] * u[1] * C + u[0] * s, c + u[2] * u[2] * C}}};
}
M3 reflection(const Vec3& n) {
  M3 R = identity();
  for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) R[i][j] -= 2 * n[i] * n[j];
  return R;
}
M3 improper(const Vec3& u, int n) { return mul(reflection(u), rotation(u, 2 * kPi / n)); }
M3 inversion() { return {{{-1, 0, 0}, {0, -1, 0}, {0, 0, -1}}}; }

// symmetric 3×3 eigen decomposition by cyclic Jacobi (eigenvectors as columns of V)
void jacobi(M3 A, std::array<double, 3>& w, M3& V) {
  V = identity();
  for (int sweep = 0; sweep < 60; ++sweep) {
    double off = std::fabs(A[0][1]) + std::fabs(A[0][2]) + std::fabs(A[1][2]);
    if (off < 1e-15) break;
    for (int p = 0; p < 2; ++p)
      for (int q = p + 1; q < 3; ++q) {
        if (std::fabs(A[p][q]) < 1e-300) continue;
        const double th = (A[q][q] - A[p][p]) / (2 * A[p][q]);
        const double t = (th >= 0 ? 1 : -1) / (std::fabs(th) + std::sqrt(th * th + 1));
        const double c = 1 / std::sqrt(t * t + 1), s = t * c;
        for (int k = 0; k < 3; ++k) {
          const double akp = A[k][p], akq = A[k][q];
          A[k][p] = c * akp - s * akq, A[k][q] = s * akp + c * akq;
        }
        for (int k = 0; k < 3; ++k) {
          const double apk = A[p][k], aqk = A[q][k];
          A[p][k] = c * apk - s * aqk, A[q][k] = s * apk + c * aqk;
        }
        for (int k = 0; k < 3; ++k) {
          const double vkp = V[k][p], vkq = V[k][q];
          V[k][p] = c * vkp - s * vkq, V[k][q] = s * vkp + c * vkq;
        }
      }
  }
  for (int k = 0; k < 3; ++k) w[size_t(k)] = A[k][k];
}

struct Mol {
  std::vector<Vec3> x;      // about the centre of mass
  std::vector<int> z;
  Vec3 centre{0, 0, 0};
};
Mol molecule_of(const System& s) {
  Mol m;
  double M = 0;
  for (const auto& a : s.atoms) {
    const double w = a.element > 0 ? element(a.element).mass : 1.0;
    m.centre = m.centre + a.pos * w, M += w;
  }
  m.centre = m.centre * (1 / M);
  for (const auto& a : s.atoms) m.x.push_back(a.pos - m.centre), m.z.push_back(a.element);
  return m;
}

// the permutation the operation makes (atom i lands on perm[i]) and the largest miss; empty when it does not hold
std::vector<int> matches(const Mol& m, const M3& R, double tol, double* dev) {
  const size_t n = m.x.size();
  std::vector<int> perm(n, -1);
  std::vector<char> used(n, 0);
  double worst = 0;
  for (size_t i = 0; i < n; ++i) {
    const Vec3 y = apply(R, m.x[i]);
    int best = -1;
    double bd = tol;
    for (size_t j = 0; j < n; ++j) {
      if (m.z[j] != m.z[i] || used[j]) continue;
      const double d = norm(y - m.x[j]);
      if (d < bd) bd = d, best = int(j);
    }
    if (best < 0) return {};
    perm[i] = best, used[size_t(best)] = 1, worst = std::max(worst, bd);
  }
  if (dev) *dev = worst;
  return perm;
}
bool holds(const Mol& m, const M3& R, double tol) { return !matches(m, R, tol, nullptr).empty(); }

Vec3 canonical(Vec3 v) {
  const double L = norm(v);
  v = v * (1 / L);
  for (int k = 0; k < 3; ++k)
    if (std::fabs(v[size_t(k)]) > 1e-6) { if (v[size_t(k)] < 0) v = v * -1.0; break; }
  return v;
}
void add_direction(std::vector<Vec3>& dirs, const Vec3& v) {
  if (norm(v) < 1e-3) return;
  const Vec3 u = canonical(v);
  for (const auto& d : dirs) if (std::fabs(dot(d, u)) > 1 - 1e-5) return;
  dirs.push_back(u);
}

struct Found {
  std::vector<std::pair<Vec3, int>> rot;   // C_n axes, the highest n on each
  std::vector<Vec3> mirrors;
  std::vector<std::pair<Vec3, int>> imp;   // S_n
  bool inv = false;
  bool linear = false;
  Vec3 linear_axis{0, 0, 1};
};

Found find_elements(const Mol& m, double tol) {
  Found f;
  const size_t n = m.x.size();
  // inertia
  M3 I{};
  for (size_t i = 0; i < n; ++i) {
    const double w = m.z[i] > 0 ? element(m.z[i]).mass : 1.0;
    const Vec3& r = m.x[i];
    const double r2 = dot(r, r);
    for (int a = 0; a < 3; ++a)
      for (int b = 0; b < 3; ++b) I[a][b] += w * ((a == b ? r2 : 0) - r[size_t(a)] * r[size_t(b)]);
  }
  std::array<double, 3> w;
  M3 V;
  jacobi(I, w, V);
  f.inv = holds(m, inversion(), tol);
  // linear: every atom within tol of the axis of least inertia
  {
    int k = 0;
    for (int q = 1; q < 3; ++q) if (w[size_t(q)] < w[size_t(k)]) k = q;
    const Vec3 ax{V[0][k], V[1][k], V[2][k]};
    bool lin = n >= 2;
    for (const auto& r : m.x) if (norm(r - ax * dot(r, ax)) > tol) { lin = false; break; }
    if (lin || n == 1) { f.linear = n >= 2; f.linear_axis = ax; return f; }
  }
  std::vector<Vec3> dirs;
  for (int k = 0; k < 3; ++k) add_direction(dirs, {V[0][k], V[1][k], V[2][k]});
  for (const auto& r : m.x) add_direction(dirs, r);
  const bool pairs = n <= 300;
  if (pairs)
    for (size_t i = 0; i < n; ++i)
      for (size_t j = i + 1; j < n; ++j) {
        if (m.z[i] != m.z[j] || std::fabs(norm(m.x[i]) - norm(m.x[j])) > 2 * tol) continue;
        add_direction(dirs, m.x[i] + m.x[j]);
        add_direction(dirs, cross(m.x[i], m.x[j]));
        add_direction(dirs, m.x[i] - m.x[j]);   // a mirror's normal between two like atoms
      }
  if (n <= 60)   // normals of planes through three like atoms at one distance (C3 axes through faces)
    for (size_t i = 0; i < n; ++i)
      for (size_t j = i + 1; j < n; ++j) {
        if (m.z[i] != m.z[j] || std::fabs(norm(m.x[i]) - norm(m.x[j])) > 2 * tol) continue;
        for (size_t k = j + 1; k < n; ++k) {
          if (m.z[k] != m.z[i] || std::fabs(norm(m.x[i]) - norm(m.x[k])) > 2 * tol) continue;
          add_direction(dirs, cross(m.x[j] - m.x[i], m.x[k] - m.x[i]));
        }
      }
  // an element found again on a nearly parallel candidate (noise) is the same element: within about 5°
  auto known = [](const std::vector<Vec3>& v, const Vec3& d) {
    for (const auto& e : v) if (std::fabs(dot(e, d)) > 0.996) return true;
    return false;
  };
  std::vector<Vec3> rot_axes;
  for (const auto& d : dirs) {
    if (!known(rot_axes, d))
      for (int k = 8; k >= 2; --k)
        if (holds(m, rotation(d, 2 * kPi / k), tol)) { f.rot.push_back({d, k}), rot_axes.push_back(d); break; }
    if (!known(f.mirrors, d) && holds(m, reflection(d), tol)) f.mirrors.push_back(d);
  }
  for (const auto& [d, k] : f.rot) {   // S_2n, else S_n (n ≥ 3), on each rotation axis: the highest that holds
    if (holds(m, improper(d, 2 * k), tol)) f.imp.push_back({d, 2 * k});
    else if (k >= 3 && holds(m, improper(d, k), tol)) f.imp.push_back({d, k});
  }
  return f;
}

std::string classify(const Found& f, Vec3* principal) {
  if (f.linear) { if (principal) *principal = f.linear_axis; return f.inv ? "D*h" : "C*v"; }
  int c3 = 0, nmax = 1;
  bool c4 = false, c5 = false;
  for (const auto& [d, k] : f.rot) {
    if (k % 3 == 0) ++c3;
    if (k == 4 || k == 8) c4 = true;
    if (k == 5) c5 = true;
    nmax = std::max(nmax, k);
  }
  if (c3 >= 2) {
    if (c5) return f.inv ? "Ih" : "I";
    if (c4) return f.inv ? "Oh" : "O";
    if (f.inv) return "Th";
    return f.mirrors.empty() ? "T" : "Td";
  }
  if (f.rot.empty()) {
    if (!f.mirrors.empty()) return "Cs";
    if (f.inv) return "Ci";
    return "C1";
  }
  Vec3 p{0, 0, 1};
  for (const auto& [d, k] : f.rot) if (k == nmax) { p = d; break; }
  if (principal) *principal = p;
  int perp_c2 = 0;
  for (const auto& [d, k] : f.rot) if (k % 2 == 0 && std::fabs(dot(d, p)) < 0.02) ++perp_c2;
  bool sh = false, sv = false;
  for (const auto& d : f.mirrors) {
    if (std::fabs(dot(d, p)) > 0.98) sh = true;
    if (std::fabs(dot(d, p)) < 0.02) sv = true;
  }
  const std::string N = std::to_string(nmax);
  if (perp_c2 >= nmax) {
    if (sh) return "D" + N + "h";
    if (sv) return "D" + N + "d";
    return "D" + N;
  }
  if (sh) return "C" + N + "h";
  if (sv) return "C" + N + "v";
  for (const auto& [d, k] : f.imp) if (k == 2 * nmax && std::fabs(dot(d, p)) > 0.98) return "S" + std::to_string(2 * nmax);
  return "C" + N;
}

// the group as permutations of the atoms: the found elements' powers, closed under composition (exact, whatever the
// noise in the axes), each with the matrix of its first finding (or of a product)
struct PermGroup {
  std::vector<std::vector<int>> perm;
  std::vector<int> sign;   // det: a planar molecule's own mirror is the identity permutation, told apart by it
  std::vector<M3> approx;
};
PermGroup group_of(const Mol& m, const Found& f, double tol) {
  PermGroup G;
  const size_t n = m.x.size();
  std::vector<int> id(n);
  for (size_t i = 0; i < n; ++i) id[i] = int(i);
  G.perm.push_back(id), G.sign.push_back(1), G.approx.push_back(identity());
  auto add = [&](const std::vector<int>& p, const M3& R) {
    if (p.empty()) return false;
    const int sg = det(R) > 0 ? 1 : -1;
    for (size_t q = 0; q < G.perm.size(); ++q) if (G.sign[q] == sg && G.perm[q] == p) return false;
    G.perm.push_back(p), G.sign.push_back(sg), G.approx.push_back(R);
    return true;
  };
  auto add_op = [&](const M3& R) { add(matches(m, R, tol, nullptr), R); };
  for (const auto& [d, k] : f.rot) for (int j = 1; j < k; ++j) add_op(rotation(d, 2 * kPi * j / k));
  for (const auto& d : f.mirrors) add_op(reflection(d));
  for (const auto& [d, k] : f.imp) add_op(improper(d, k));
  if (f.inv) add_op(inversion());
  for (bool grew = true; grew && G.perm.size() < 240;) {
    grew = false;
    const size_t c = G.perm.size();
    for (size_t a = 0; a < c && G.perm.size() < 240; ++a)
      for (size_t b = 0; b < c && G.perm.size() < 240; ++b) {
        std::vector<int> p(n);
        for (size_t i = 0; i < n; ++i) p[i] = G.perm[a][size_t(G.perm[b][i])];   // b, then a
        grew = add(p, mul(G.approx[a], G.approx[b])) || grew;
      }
  }
  return G;
}

// the orthogonal matrix that best carries each x_i onto x_perm(i) (polar factor of Σ x_π xᵀ, steadied by the
// approximate matrix where the atoms span no direction, as the normal of a planar molecule)
M3 fitted(const Mol& m, const std::vector<int>& perm, const M3& approx) {
  M3 H{};
  double scale = 0;
  for (size_t i = 0; i < m.x.size(); ++i) {
    const Vec3& a = m.x[size_t(perm[i])];
    const Vec3& b = m.x[i];
    for (int r = 0; r < 3; ++r) for (int c = 0; c < 3; ++c) H[r][c] += a[size_t(r)] * b[size_t(c)];
    scale += dot(b, b);
  }
  for (int r = 0; r < 3; ++r) for (int c = 0; c < 3; ++c) H[r][c] += 1e-6 * (scale + 1e-12) * approx[r][c];
  // R = H (HᵀH)^(−1/2)
  M3 HtH{};
  for (int r = 0; r < 3; ++r) for (int c = 0; c < 3; ++c) for (int k = 0; k < 3; ++k) HtH[r][c] += H[k][r] * H[k][c];
  std::array<double, 3> w;
  M3 V;
  jacobi(HtH, w, V);
  M3 S{};   // (HᵀH)^(−1/2) = V diag(1/√w) Vᵀ
  for (int r = 0; r < 3; ++r)
    for (int c = 0; c < 3; ++c)
      for (int k = 0; k < 3; ++k) S[r][c] += V[r][k] * V[c][k] / std::sqrt(std::max(w[size_t(k)], 1e-300));
  return mul(H, S);
}

PointGroupReport report_of(const Mol& m, const Found& f, double tol) {
  PointGroupReport r;
  r.centre = m.centre;
  r.group = classify(f, nullptr);
  if (f.linear) {
    r.order = 0;
    r.elements.push_back({"C", 0, f.linear_axis});
    if (f.inv) r.elements.push_back({"i", 2, {0, 0, 0}});
    double dev = 0;
    for (const auto& x : m.x) dev = std::max(dev, norm(x - f.linear_axis * dot(x, f.linear_axis)));
    r.max_deviation = dev;
    return r;
  }
  for (const auto& [d, k] : f.rot) r.elements.push_back({"C", k, d});
  for (const auto& [d, k] : f.imp) r.elements.push_back({"S", k, d});
  for (const auto& d : f.mirrors) r.elements.push_back({"sigma", 1, d});
  if (f.inv) r.elements.push_back({"i", 2, {0, 0, 0}});
  const PermGroup G = group_of(m, f, tol);
  r.order = int(G.perm.size());
  double worst = 0;
  for (size_t g = 0; g < G.perm.size(); ++g) {
    const M3 R = fitted(m, G.perm[g], G.approx[g]);
    for (size_t i = 0; i < m.x.size(); ++i) worst = std::max(worst, norm(apply(R, m.x[i]) - m.x[size_t(G.perm[g][i])]));
  }
  r.max_deviation = worst;
  return r;
}

}  // namespace

PointGroupReport point_group(const System& s, double tol) {
  if (s.atoms.empty()) throw std::invalid_argument("no atoms");
  if (tol <= 0) throw std::invalid_argument("the tolerance must be positive");
  const Mol m = molecule_of(s);
  return report_of(m, find_elements(m, tol), tol);
}

PointGroupReport symmetrize_molecule(System& s, double tol) {
  if (s.atoms.empty()) throw std::invalid_argument("no atoms");
  Mol m = molecule_of(s);
  const std::string group0 = report_of(m, find_elements(m, tol), tol).group;
  for (int it = 0; it < 30; ++it) {
    const Found f = find_elements(m, tol);
    if (f.linear) {   // onto the axis; D∞h: each atom's partner through the centre mirrored exactly
      std::vector<Vec3> y(m.x.size());
      for (size_t i = 0; i < m.x.size(); ++i) y[i] = f.linear_axis * dot(m.x[i], f.linear_axis);
      if (f.inv) {
        double dv = 0;
        const auto perm = matches(m, inversion(), tol, &dv);
        if (!perm.empty())
          for (size_t i = 0; i < y.size(); ++i) y[i] = (y[i] - y[size_t(perm[i])]) * 0.5;
      }
      m.x = y;
      continue;
    }
    const PermGroup G = group_of(m, f, tol);
    std::vector<Vec3> y(m.x.size(), Vec3{0, 0, 0});
    for (size_t g = 0; g < G.perm.size(); ++g) {
      const M3 R = fitted(m, G.perm[g], G.approx[g]);
      for (size_t i = 0; i < m.x.size(); ++i) y[i] = y[i] + apply_t(R, m.x[size_t(G.perm[g][i])]);
    }
    double moved = 0;
    for (size_t i = 0; i < y.size(); ++i) {
      y[i] = y[i] * (1.0 / double(G.perm.size()));
      moved = std::max(moved, norm(y[i] - m.x[i]));
    }
    m.x = y;
    if (moved < 1e-12) break;
  }
  for (size_t i = 0; i < s.atoms.size(); ++i) s.atoms[i].pos = m.centre + m.x[i];
  PointGroupReport r = point_group(s, tol);
  if (r.group != group0) r.notes.push_back("the group changed from " + group0 + " while symmetrising");
  return r;
}

}  // namespace caps
