// Molecular point groups and symmetrisation (pointgroup.hpp).
#include "caps/pointgroup.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
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

// The cycle lengths of a permutation: its order is their least common multiple (capped, as a guard)
long perm_order(const std::vector<int>& p) {
  std::vector<char> seen(p.size(), 0);
  long L = 1;
  for (size_t i = 0; i < p.size(); ++i) {
    if (seen[i]) continue;
    long len = 0;
    for (size_t j = i; !seen[j]; j = size_t(p[j])) seen[j] = 1, ++len;
    long a = L, b = len;
    while (b) { const long t = a % b; a = b, b = t; }
    L = L / a * len;
    if (L > 100000) return 100000;
  }
  return L;
}
// a rotation C_k really of order k: every atom is on the axis (a 1-cycle) or in a cycle of exactly k atoms, and some are.
// A ring of six then never passes as C7 or C8 at a loose tolerance, nor a near-linear molecule as any C_k.
bool cycles_of(const std::vector<int>& p, int k) {
  if (p.empty()) return false;
  std::vector<char> seen(p.size(), 0);
  bool some = false;
  for (size_t i = 0; i < p.size(); ++i) {
    if (seen[i]) continue;
    int len = 0;
    for (size_t j = i; !seen[j]; j = size_t(p[j])) seen[j] = 1, ++len;
    if (len != 1 && len != k) return false;
    some = some || len == k;
  }
  return some;
}
// the order of an operation given as its permutation and determinant: a planar molecule's own mirror permutes nothing
int element_order(const std::vector<int>& p, int sign) {
  const long L = perm_order(p);
  return int(sign < 0 && L % 2 == 1 ? 2 * L : L);
}

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
  {
    const auto p = matches(m, inversion(), tol, nullptr);
    f.inv = !p.empty() && perm_order(p) <= 2;
  }
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
        if (cycles_of(matches(m, rotation(d, 2 * kPi / k), tol, nullptr), k)) { f.rot.push_back({d, k}), rot_axes.push_back(d); break; }
    if (!known(f.mirrors, d)) {
      const auto p = matches(m, reflection(d), tol, nullptr);
      if (!p.empty() && perm_order(p) <= 2) f.mirrors.push_back(d);
    }
  }
  // S_n holds when its permutation has the order S_n has: n for even n, 2n for odd n
  auto s_holds = [&](const Vec3& d, int n) {
    const auto p = matches(m, improper(d, n), tol, nullptr);
    return !p.empty() && element_order(p, -1) == (n % 2 == 0 ? n : 2 * n);
  };
  for (const auto& [d, k] : f.rot) {   // S_2n, else S_n (n ≥ 3), on each rotation axis: the highest that holds
    if (s_holds(d, 2 * k)) f.imp.push_back({d, 2 * k});
    else if (k >= 3 && s_holds(d, k)) f.imp.push_back({d, k});
  }
  return f;
}

// the orthogonal matrix that best carries each x_i onto x_perm(i) (polar factor of Σ x_π xᵀ, steadied by the
// approximate matrix where the atoms span no direction, as the normal of a planar molecule)
M3 fit_once(const Mol& m, const std::vector<int>& perm, const M3& approx) {
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

// The small term in fit_once fixes a planar or linear molecule's rotation but pulls the result towards `approx`; refitted
// from its own answer that pull vanishes (it was about 1e-7 Å from the first estimate).
M3 fitted(const Mol& m, const std::vector<int>& perm, const M3& approx) {
  M3 R = fit_once(m, perm, approx);
  for (int pass = 0; pass < 2; ++pass) R = fit_once(m, perm, R);
  return R;
}

// The group as permutations of the atoms with their determinants (a planar molecule's own mirror is the identity
// permutation, told apart by its sign), closed under composition. Permutations compose exactly whatever the noise, so
// the closure can reach operations that no single test saw: each operation is therefore checked here against the
// tolerance. It keeps the product of its generators' ideal matrices, the orthogonal matrix fitted to its permutation,
// and its miss: the largest distance of an atom's image from its partner under the better of the two.
struct PermGroup {
  std::vector<std::vector<int>> perm;
  std::vector<int> sign;
  std::vector<M3> approx;
  std::vector<M3> fit;
  std::vector<double> miss;
  double worst = 0;
  bool contains(const std::vector<int>& p, int sg) const {
    for (size_t q = 0; q < perm.size(); ++q) if (sign[q] == sg && perm[q] == p) return true;
    return false;
  }
};
struct Op {
  std::vector<int> perm;
  int sign = 1;
  M3 R{};
};

double image_miss(const Mol& m, const std::vector<int>& p, const M3& R) {
  double w = 0;
  for (size_t i = 0; i < m.x.size(); ++i) w = std::max(w, norm(apply(R, m.x[i]) - m.x[size_t(p[i])]));
  return w;
}
double trace(const M3& R) { return R[0][0] + R[1][1] + R[2][2]; }

constexpr size_t kMaxOrder = 120;   // Ih, the largest finite point group a molecule can have

// Closes the generators into G. False when an operation misses by more than tol (with check) or the set outgrows Ih.
bool close_group(const Mol& m, const std::vector<Op>& gens, double tol, bool check, PermGroup& G) {
  G = PermGroup{};
  const size_t n = m.x.size();
  std::vector<int> id(n);
  for (size_t i = 0; i < n; ++i) id[i] = int(i);
  G.perm.push_back(id), G.sign.push_back(1), G.approx.push_back(identity()), G.fit.push_back(identity()), G.miss.push_back(0);
  for (size_t a = 0; a < G.perm.size(); ++a)
    for (const auto& s : gens) {
      std::vector<int> p(n);
      for (size_t i = 0; i < n; ++i) p[i] = G.perm[a][size_t(s.perm[i])];   // s, then a
      const int sg = G.sign[a] * s.sign;
      if (G.contains(p, sg)) continue;
      if (G.perm.size() >= kMaxOrder) return false;
      const M3 A = mul(G.approx[a], s.R);
      M3 F = fitted(m, p, A);
      double d = image_miss(m, p, A);
      if (det(F) * sg > 0) d = std::min(d, image_miss(m, p, F));
      else F = A;
      if (check && d > tol) return false;
      G.perm.push_back(std::move(p)), G.sign.push_back(sg), G.approx.push_back(A), G.fit.push_back(F), G.miss.push_back(d);
      G.worst = std::max(G.worst, d);
    }
  return true;
}

// The Schoenflies symbol of a closed group, read from the group itself: its order, the orders of its rotations, and the
// reflections and inversion among its improper operations (told apart by the trace, 1 against -3). Empty when the
// counts fit no point group, so a name never disagrees with the order.
std::string name_of(const PermGroup& G) {
  const size_t N = G.perm.size();
  int P = 0, maxord = 1, c2 = 0, refl = 0;
  bool inv = false, cyclic = false;
  std::vector<int> ord(N);
  for (size_t g = 0; g < N; ++g) ord[g] = element_order(G.perm[g], G.sign[g]);
  for (size_t g = 0; g < N; ++g) {
    if (G.sign[g] > 0) {
      ++P, maxord = std::max(maxord, ord[g]);
      if (ord[g] == 2) ++c2;
    } else if (ord[g] == 2) {
      if (trace(G.fit[g]) < -1) inv = true;
      else ++refl;
    }
  }
  for (size_t g = 0; g < N; ++g) if (G.sign[g] > 0 && ord[g] == P) cyclic = true;
  char kind = 0;
  int n = 1;
  if (P == 1) kind = 'C';
  else if (cyclic && c2 == (P % 2 == 0 ? 1 : 0)) kind = 'C', n = P;
  else if (P == 12 && maxord == 3 && c2 == 3) kind = 'T';
  else if (P == 24 && maxord == 4 && c2 == 9) kind = 'O';
  else if (P == 60 && maxord == 5 && c2 == 15) kind = 'I';
  else if (2 * maxord == P && c2 == maxord + (maxord % 2 == 0 ? 1 : 0)) kind = 'D', n = maxord;
  else return "";
  const std::string N_ = std::to_string(n);
  if (size_t(P) == N) {
    if (kind == 'C' || kind == 'D') return std::string(1, kind) + N_;
    return std::string(1, kind);
  }
  if (N != size_t(2 * P)) return "";
  switch (kind) {
    case 'C':
      if (n == 1) return refl == 1 && !inv ? "Cs" : refl == 0 && inv ? "Ci" : "";
      if (refl == 0 && inv == (n % 2 == 1)) return "S" + std::to_string(2 * n);
      if (refl == 1 && inv == (n % 2 == 0)) return "C" + N_ + "h";
      if (refl == n && !inv) return "C" + N_ + "v";
      return "";
    case 'D':
      if (refl == n + 1 && inv == (n % 2 == 0)) return "D" + N_ + "h";
      if (refl == n && inv == (n % 2 == 1)) return "D" + N_ + "d";
      return "";
    case 'T':
      if (inv && refl == 3) return "Th";
      if (!inv && refl == 6) return "Td";
      return "";
    case 'O': return inv && refl == 9 ? "Oh" : "";
    case 'I': return inv && refl == 15 ? "Ih" : "";
  }
  return "";
}

// the axis of a proper rotation: the eigenvector of its symmetric part with the largest eigenvalue (1, against cos θ)
Vec3 rotation_axis(const M3& R) {
  M3 S{};
  for (int r = 0; r < 3; ++r) for (int c = 0; c < 3; ++c) S[r][c] = 0.5 * (R[r][c] + R[c][r]);
  std::array<double, 3> w;
  M3 V;
  jacobi(S, w, V);
  int k = 0;
  for (int q = 1; q < 3; ++q) if (w[size_t(q)] > w[size_t(k)]) k = q;
  return canonical({V[0][k], V[1][k], V[2][k]});
}

// The symmetry elements of the group, one per distinct axis or plane (the highest order on each axis)
void elements_of(const PermGroup& G, PointGroupReport& r) {
  std::vector<SymmetryElement> rot, imp, mir;
  bool inv = false;
  auto merge = [](std::vector<SymmetryElement>& v, SymmetryElement e) {
    for (auto& x : v) if (std::fabs(dot(x.axis, e.axis)) > 0.98) { x.n = std::max(x.n, e.n); return; }
    v.push_back(std::move(e));
  };
  const size_t n = G.perm.empty() ? 0 : G.perm[0].size();
  for (size_t g = 0; g < G.perm.size(); ++g) {
    const int k = element_order(G.perm[g], G.sign[g]);
    if (k < 2) continue;
    if (G.sign[g] > 0) { merge(rot, {"C", k, rotation_axis(G.fit[g])}); continue; }
    M3 neg = G.fit[g];   // −(σ C_θ) is the rotation by θ + π about the same axis
    for (auto& row : neg) for (auto& v : row) v = -v;
    if (k == 2) {
      if (trace(G.fit[g]) < -1) inv = true;
      else merge(mir, {"sigma", 1, rotation_axis(neg)});
      continue;
    }
    int s = k;   // order 4m: S_4m; order 2 (mod 4): S_k when its k/2-th power is the inversion, else S_(k/2) with k/2 odd
    if (k % 4 != 0) {
      std::vector<int> q(n);
      for (size_t i = 0; i < n; ++i) q[i] = int(i);
      for (int t = 0; t < k / 2; ++t) {
        std::vector<int> nq(n);
        for (size_t i = 0; i < n; ++i) nq[i] = G.perm[g][size_t(q[i])];
        q.swap(nq);
      }
      for (size_t h = 0; h < G.perm.size(); ++h)
        if (G.sign[h] < 0 && G.perm[h] == q) { if (trace(G.fit[h]) >= -1) s = k / 2; break; }
    }
    merge(imp, {"S", s, rotation_axis(neg)});
  }
  std::stable_sort(rot.begin(), rot.end(), [](const SymmetryElement& a, const SymmetryElement& b) { return a.n > b.n; });
  std::stable_sort(imp.begin(), imp.end(), [](const SymmetryElement& a, const SymmetryElement& b) { return a.n > b.n; });
  for (auto& e : rot) r.elements.push_back(e);
  for (auto& e : imp) r.elements.push_back(e);
  for (auto& e : mir) r.elements.push_back(e);
  if (inv) r.elements.push_back({"i", 2, {0, 0, 0}});
}

// The operations to build the group from: each element found and the lower rotations on its axis (so that a C6 whose
// products fail can still give its C3 or C2), each with the order it must have, in the order the search prefers them
std::vector<Op> candidates(const Mol& m, const Found& f, double tol) {
  std::vector<Op> c;
  auto add = [&](std::vector<int> p, int sign, const M3& R) {
    for (const auto& o : c) if (o.sign == sign && o.perm == p) return;
    c.push_back({std::move(p), sign, R});
  };
  std::vector<std::pair<Vec3, int>> rots;
  for (const auto& [d, k] : f.rot) for (int q = k; q >= 2; --q) if (k % q == 0) rots.push_back({d, q});
  std::stable_sort(rots.begin(), rots.end(), [](const auto& a, const auto& b) { return a.second > b.second; });
  for (const auto& [d, q] : rots) {
    const M3 R = rotation(d, 2 * kPi / q);
    auto p = matches(m, R, tol, nullptr);
    if (cycles_of(p, q)) add(std::move(p), 1, R);
  }
  auto improper_op = [&](const M3& R, int order) {
    auto p = matches(m, R, tol, nullptr);
    if (!p.empty() && element_order(p, -1) == order) add(std::move(p), -1, R);
  };
  if (f.inv) improper_op(inversion(), 2);
  for (const auto& d : f.mirrors) improper_op(reflection(d), 2);
  for (const auto& [d, k] : f.imp) improper_op(improper(d, k), k % 2 == 0 ? k : 2 * k);
  return c;
}

// The group: every candidate if they close within the tolerance into a point group; otherwise the largest group found by
// adding candidates one at a time (from several starting ones), each kept only when the closure still holds and is a
// point group. Every operation of the result is within tol, and its name comes from the result itself.
struct Accepted {
  PermGroup G;
  std::string name = "C1";
  size_t used = 0;   // candidates inside G
};
Accepted accept(const Mol& m, const std::vector<Op>& cands, double tol) {
  Accepted best;
  {
    PermGroup G;
    if (close_group(m, cands, tol, true, G)) {
      const std::string nm = name_of(G);
      if (!nm.empty()) { best.G = std::move(G), best.name = nm, best.used = cands.size(); return best; }
    }
  }
  close_group(m, {}, tol, true, best.G);
  int seeds = 0;
  for (size_t s = 0; s < cands.size() && seeds < 24; ++s) {
    if (best.G.contains(cands[s].perm, cands[s].sign)) continue;   // a run from it would most likely find this group again
    ++seeds;
    std::vector<Op> gens;
    PermGroup cur;
    close_group(m, {}, tol, true, cur);
    std::string name = "C1";
    auto try_add = [&](const Op& o) {
      if (cur.contains(o.perm, o.sign)) return;
      gens.push_back(o);
      PermGroup G;
      std::string nm;
      if (close_group(m, gens, tol, true, G) && !(nm = name_of(G)).empty()) cur = std::move(G), name = nm;
      else gens.pop_back();
    };
    try_add(cands[s]);
    for (const auto& o : cands) try_add(o);
    if (cur.perm.size() > best.G.perm.size() || (cur.perm.size() == best.G.perm.size() && cur.worst < best.G.worst))
      best.G = std::move(cur), best.name = name;
  }
  for (const auto& o : cands) if (best.G.contains(o.perm, o.sign)) ++best.used;
  return best;
}

PointGroupReport report_of(const Mol& m, const Found& f, double tol) {
  PointGroupReport r;
  r.centre = m.centre;
  if (f.linear) {
    r.group = f.inv ? "D*h" : "C*v";
    r.order = 0;
    r.elements.push_back({"C", 0, f.linear_axis});
    if (f.inv) r.elements.push_back({"i", 2, {0, 0, 0}});
    double dev = 0;
    for (const auto& x : m.x) dev = std::max(dev, norm(x - f.linear_axis * dot(x, f.linear_axis)));
    r.max_deviation = dev;
    return r;
  }
  const std::vector<Op> cands = candidates(m, f, tol);
  const Accepted A = accept(m, cands, tol);
  r.group = A.name;
  r.order = int(A.G.perm.size());
  r.max_deviation = A.G.worst;
  elements_of(A.G, r);
  if (A.used < cands.size()) {   // say what was left out and why
    PermGroup all;
    char buf[400];
    if (!close_group(m, cands, tol, false, all))
      std::snprintf(buf, sizeof buf, "%zu of the %zu operations that hold one at a time are left out: together they generate more than %zu operations, which no point group has", cands.size() - A.used, cands.size(), kMaxOrder);
    else {
      const std::string nm = name_of(all);
      if (all.worst > tol)
        std::snprintf(buf, sizeof buf, "%zu of the %zu operations that hold one at a time are left out: together they close into %s of order %zu with a largest image miss of %.6f Å, above the tolerance %.3f Å", cands.size() - A.used, cands.size(), nm.empty() ? "a group" : nm.c_str(), all.perm.size(), all.worst, tol);
      else
        std::snprintf(buf, sizeof buf, "%zu of the %zu operations that hold one at a time are left out: together they close into a set of order %zu that is not a point group", cands.size() - A.used, cands.size(), all.perm.size());
    }
    r.notes.push_back(buf);
  }
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
    const PermGroup G = accept(m, candidates(m, f, tol), tol).G;   // every operation within tol
    std::vector<Vec3> y(m.x.size(), Vec3{0, 0, 0});
    for (size_t g = 0; g < G.perm.size(); ++g) {
      const M3& R = G.fit[g];
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
