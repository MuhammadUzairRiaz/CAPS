// Horn 1987 quaternion superposition (see superpose.hpp).
#include "caps/superpose.hpp"

#include <cmath>
#include <stdexcept>

namespace caps {

namespace {

// Largest eigenvalue's eigenvector of a symmetric 4×4 matrix (cyclic Jacobi).
std::array<double, 4> top_eigenvector(std::array<std::array<double, 4>, 4> a) {
  std::array<std::array<double, 4>, 4> v{};
  for (int i = 0; i < 4; ++i) v[i][i] = 1;
  for (int sweep = 0; sweep < 60; ++sweep) {
    double off = 0;
    for (int p = 0; p < 4; ++p)
      for (int q = p + 1; q < 4; ++q) off += a[p][q] * a[p][q];
    if (off < 1e-30) break;
    for (int p = 0; p < 4; ++p)
      for (int q = p + 1; q < 4; ++q) {
        if (std::abs(a[p][q]) < 1e-300) continue;
        const double theta = (a[q][q] - a[p][p]) / (2 * a[p][q]);
        const double t = (theta >= 0 ? 1.0 : -1.0) / (std::abs(theta) + std::sqrt(theta * theta + 1));
        const double c = 1 / std::sqrt(t * t + 1), s = t * c;
        for (int k = 0; k < 4; ++k) {
          const double akp = a[k][p], akq = a[k][q];
          a[k][p] = c * akp - s * akq;
          a[k][q] = s * akp + c * akq;
        }
        for (int k = 0; k < 4; ++k) {
          const double apk = a[p][k], aqk = a[q][k];
          a[p][k] = c * apk - s * aqk;
          a[q][k] = s * apk + c * aqk;
        }
        for (int k = 0; k < 4; ++k) {
          const double vkp = v[k][p], vkq = v[k][q];
          v[k][p] = c * vkp - s * vkq;
          v[k][q] = s * vkp + c * vkq;
        }
      }
  }
  int best = 0;
  for (int i = 1; i < 4; ++i)
    if (a[i][i] > a[best][best]) best = i;
  return {v[0][best], v[1][best], v[2][best], v[3][best]};
}

}  // namespace

Vec3 Superposition::apply(const Vec3& x) const {
  const Vec3 d = x - centre_mov;
  return {rot[0][0] * d[0] + rot[0][1] * d[1] + rot[0][2] * d[2] + centre_ref[0], rot[1][0] * d[0] + rot[1][1] * d[1] + rot[1][2] * d[2] + centre_ref[1],
          rot[2][0] * d[0] + rot[2][1] * d[1] + rot[2][2] * d[2] + centre_ref[2]};
}

Superposition superpose(const std::vector<Vec3>& ref, const std::vector<Vec3>& mov, const std::vector<double>& w) {
  if (ref.size() != mov.size()) throw std::invalid_argument("the two states have different atom counts");
  if (!w.empty() && w.size() != ref.size()) throw std::invalid_argument("one weight per atom");
  Superposition s;
  double wsum = 0;
  for (size_t i = 0; i < ref.size(); ++i) {
    const double wi = w.empty() ? 1.0 : w[i];
    if (wi <= 0) continue;
    s.centre_ref = s.centre_ref + ref[i] * wi;
    s.centre_mov = s.centre_mov + mov[i] * wi;
    wsum += wi;
    ++s.fitted;
  }
  if (s.fitted == 0) throw std::invalid_argument("no atoms to fit on");
  s.centre_ref = s.centre_ref * (1 / wsum);
  s.centre_mov = s.centre_mov * (1 / wsum);
  // correlation matrix S_ab = Σ w (mov − c_mov)_a (ref − c_ref)_b
  double S[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
  for (size_t i = 0; i < ref.size(); ++i) {
    const double wi = w.empty() ? 1.0 : w[i];
    if (wi <= 0) continue;
    const Vec3 a = mov[i] - s.centre_mov, b = ref[i] - s.centre_ref;
    for (int p = 0; p < 3; ++p)
      for (int q = 0; q < 3; ++q) S[p][q] += wi * a[p] * b[q];
  }
  if (s.fitted >= 3) {
    const double xx = S[0][0], xy = S[0][1], xz = S[0][2], yx = S[1][0], yy = S[1][1], yz = S[1][2], zx = S[2][0], zy = S[2][1], zz = S[2][2];
    const std::array<std::array<double, 4>, 4> N{{{xx + yy + zz, yz - zy, zx - xz, xy - yx},
                                                  {yz - zy, xx - yy - zz, xy + yx, zx + xz},
                                                  {zx - xz, xy + yx, -xx + yy - zz, yz + zy},
                                                  {xy - yx, zx + xz, yz + zy, -xx - yy + zz}}};
    const auto q = top_eigenvector(N);
    const double q0 = q[0], q1 = q[1], q2 = q[2], q3 = q[3];
    s.rot = {{{q0 * q0 + q1 * q1 - q2 * q2 - q3 * q3, 2 * (q1 * q2 - q0 * q3), 2 * (q1 * q3 + q0 * q2)},
              {2 * (q1 * q2 + q0 * q3), q0 * q0 - q1 * q1 + q2 * q2 - q3 * q3, 2 * (q2 * q3 - q0 * q1)},
              {2 * (q1 * q3 - q0 * q2), 2 * (q2 * q3 + q0 * q1), q0 * q0 - q1 * q1 - q2 * q2 + q3 * q3}}};
  }
  double e = 0;
  for (size_t i = 0; i < ref.size(); ++i) {
    const double wi = w.empty() ? 1.0 : w[i];
    if (wi <= 0) continue;
    const Vec3 d = s.apply(mov[i]) - ref[i];
    e += wi * dot(d, d);
  }
  s.rmsd = std::sqrt(e / wsum);
  return s;
}

double rmsd_after(const Superposition& fit, const std::vector<Vec3>& ref, const std::vector<Vec3>& mov, const std::vector<char>& mask) {
  double e = 0;
  size_t n = 0;
  for (size_t i = 0; i < ref.size() && i < mov.size(); ++i) {
    if (!mask.empty() && !mask[i]) continue;
    const Vec3 d = fit.apply(mov[i]) - ref[i];
    e += dot(d, d);
    ++n;
  }
  return n ? std::sqrt(e / double(n)) : 0.0;
}

}  // namespace caps
