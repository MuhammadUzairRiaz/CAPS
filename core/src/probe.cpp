#include "caps/probe.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace caps {

ProbeKind probe_kind(const std::string& n) {
  if (n == "point" || n == "centre") return ProbeKind::Point;
  if (n == "plane") return ProbeKind::Plane;
  if (n == "axis") return ProbeKind::Axis;
  if (n == "ellipsoid") return ProbeKind::Ellipsoid;
  throw std::invalid_argument("a probe is a point, plane, axis or ellipsoid (not '" + n + "')");
}

namespace {
// Jacobi eigen-decomposition of a symmetric 3×3 matrix: eigenvalues descending, eigenvectors as columns.
void eigen3(double a[3][3], double w[3], double v[3][3]) {
  for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) v[i][j] = i == j;
  for (int sweep = 0; sweep < 60; ++sweep) {
    const double off = a[0][1] * a[0][1] + a[0][2] * a[0][2] + a[1][2] * a[1][2];
    if (off < 1e-22) break;
    for (int p = 0; p < 2; ++p)
      for (int q = p + 1; q < 3; ++q) {
        if (std::abs(a[p][q]) < 1e-300) continue;
        const double th = 0.5 * std::atan2(2 * a[p][q], a[q][q] - a[p][p]);
        const double c = std::cos(th), s = std::sin(th);
        for (int k = 0; k < 3; ++k) {   // A ← Jᵀ A J
          const double akp = a[k][p], akq = a[k][q];
          a[k][p] = c * akp - s * akq, a[k][q] = s * akp + c * akq;
        }
        for (int k = 0; k < 3; ++k) {
          const double apk = a[p][k], aqk = a[q][k];
          a[p][k] = c * apk - s * aqk, a[q][k] = s * apk + c * aqk;
        }
        for (int k = 0; k < 3; ++k) {
          const double vkp = v[k][p], vkq = v[k][q];
          v[k][p] = c * vkp - s * vkq, v[k][q] = s * vkp + c * vkq;
        }
      }
  }
  int idx[3] = {0, 1, 2};
  std::sort(idx, idx + 3, [&](int x, int y) { return a[x][x] > a[y][y]; });
  double vv[3][3];
  for (int k = 0; k < 3; ++k) { w[k] = a[idx[k]][idx[k]]; for (int r = 0; r < 3; ++r) vv[r][k] = v[r][idx[k]]; }
  for (int r = 0; r < 3; ++r) for (int k = 0; k < 3; ++k) v[r][k] = vv[r][k];
}
}  // namespace

Probe make_probe(const System& s, const std::vector<size_t>& atoms, ProbeKind kind) {
  std::vector<size_t> list;
  for (size_t i : atoms) if (i < s.atoms.size()) list.push_back(i);
  if (list.empty()) throw std::invalid_argument("a probe needs atoms");
  std::vector<Vec3> p(list.size());
  std::vector<double> m(list.size());
  double M = 0;
  Vec3 c{0, 0, 0};
  for (size_t k = 0; k < list.size(); ++k) {
    const Vec3 d = s.atoms[list[k]].pos - s.atoms[list[0]].pos;
    p[k] = s.atoms[list[0]].pos + (s.cell.valid() ? s.cell.minimum_image(d) : d);
    m[k] = std::max(1e-6, s.mass_of(s.atoms[list[k]]));
    c = c + p[k] * m[k];
    M += m[k];
  }
  Probe pr;
  pr.kind = kind;
  pr.centre = c * (1.0 / M);
  double cov[3][3] = {};
  for (size_t k = 0; k < p.size(); ++k) {
    const Vec3 d = p[k] - pr.centre;
    for (int r = 0; r < 3; ++r) for (int q = 0; q < 3; ++q) cov[r][q] += m[k] * d[r] * d[q] / M;
  }
  double w[3], v[3][3];
  eigen3(cov, w, v);
  for (int k = 0; k < 3; ++k) {
    Vec3 e{v[0][k], v[1][k], v[2][k]};
    // one sign for each direction (its largest component positive): a plane's "above" does not flip from frame to frame
    int big = 0;
    for (int q = 1; q < 3; ++q) if (std::abs(e[q]) > std::abs(e[big])) big = q;
    if (e[big] < 0) e = e * -1.0;
    pr.axes[size_t(k)] = e;
  }
  pr.rms = std::sqrt(std::max(0.0, w[2]));
  // the inertia-shaped ellipsoid scaled the least to hold every atom
  double scale = 0;
  for (const auto& q : p) {
    const Vec3 d = q - pr.centre;
    double t = 0;
    for (int k = 0; k < 3; ++k) { const double x = dot(d, pr.axes[size_t(k)]); t += x * x / std::max(w[k], 1e-6); }
    scale = std::max(scale, std::sqrt(t));
  }
  for (int k = 0; k < 3; ++k) pr.semi[size_t(k)] = scale * std::sqrt(std::max(w[k], 1e-6));
  return pr;
}

double probe_measure(const System& s, const Probe& a, const Probe* b, const std::string& what) {
  if (what == "rms") return a.rms;
  if (what == "size") return a.kind == ProbeKind::Plane ? a.rms : a.semi[0];
  if (!b) throw std::invalid_argument(what + " needs a second probe");
  if (what == "distance") {
    Vec3 d = a.centre - b->centre;
    if (s.cell.valid()) d = s.cell.minimum_image(d);
    if (b->kind == ProbeKind::Plane) return dot(d, b->axes[2]);    // above the plane along its normal
    if (a.kind == ProbeKind::Plane) return -dot(d, a.axes[2]);
    return norm(d);
  }
  if (what == "angle") {
    const double c = std::clamp(std::abs(dot(a.direction(), b->direction())), 0.0, 1.0);
    const bool line_plane = (a.kind == ProbeKind::Plane) != (b->kind == ProbeKind::Plane);
    const double deg = std::acos(c) * 180.0 / M_PI;
    return line_plane ? 90.0 - deg : deg;   // a line against a plane: the angle to the plane, not to its normal
  }
  throw std::invalid_argument("a probe measures distance, angle, rms or size (not '" + what + "')");
}

}  // namespace caps
