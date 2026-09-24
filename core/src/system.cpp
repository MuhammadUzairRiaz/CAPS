#include "caps/system.hpp"

#include <cmath>
#include <numeric>
#include <unordered_map>

#include "caps/elements.hpp"

namespace caps {

double norm(const Vec3& a) { return std::sqrt(dot(a, a)); }

bool Cell::valid() const { return volume() > 1e-9; }

double Cell::volume() const { return std::fabs(dot(a, cross(b, c))); }

Vec3 Cell::to_fractional(const Vec3& r) const {
  // Solve r - origin = f0 a + f1 b + f2 c with the reciprocal vectors.
  const Vec3 d = r - origin;
  const double v = dot(a, cross(b, c));
  return {dot(d, cross(b, c)) / v, dot(d, cross(c, a)) / v, dot(d, cross(a, b)) / v};
}

Vec3 Cell::to_cartesian(const Vec3& f) const { return origin + a * f[0] + b * f[1] + c * f[2]; }

Vec3 Cell::minimum_image(Vec3 d) const {
  if (!valid()) return d;
  const double v = dot(a, cross(b, c));
  double f[3] = {dot(d, cross(b, c)) / v, dot(d, cross(c, a)) / v, dot(d, cross(a, b)) / v};
  for (int k = 0; k < 3; ++k)
    if (periodic[k]) f[k] -= std::round(f[k]);
  Vec3 r = a * f[0] + b * f[1] + c * f[2];
  // For strongly skewed cells the rounded image can miss the true nearest; check neighbours.
  double best = dot(r, r);
  for (int i = -1; i <= 1; ++i)
    for (int j = -1; j <= 1; ++j)
      for (int k = -1; k <= 1; ++k) {
        if ((i && !periodic[0]) || (j && !periodic[1]) || (k && !periodic[2])) continue;
        const Vec3 t = r + a * i + b * j + c * k;
        const double q = dot(t, t);
        if (q < best - 1e-12) { best = q; d = t; }
      }
  return best == dot(r, r) ? r : d;
}

Vec3 Cell::wrap(const Vec3& r) const {
  if (!valid()) return r;
  Vec3 f = to_fractional(r);
  for (int k = 0; k < 3; ++k)
    if (periodic[k]) f[k] -= std::floor(f[k]);
  return to_cartesian(f);
}

double System::mass_of(const Atom& at) const {
  for (const auto& t : types)
    if (t.type == at.type && t.type != 0 && t.mass > 0) return t.mass;
  return element(at.element).mass;
}

double System::total_mass() const {
  double m = 0;
  for (const auto& at : atoms) m += mass_of(at);
  return m;
}

double System::density() const {
  if (!cell.valid()) return 0.0;
  constexpr double kNA = 6.02214076e23;
  return total_mass() / kNA / (cell.volume() * 1e-24);
}

std::vector<std::vector<uint32_t>> System::neighbours() const {
  std::vector<std::vector<uint32_t>> nb(atoms.size());
  for (const auto& b : bonds) {
    nb[b.i].push_back(b.j);
    nb[b.j].push_back(b.i);
  }
  return nb;
}

std::vector<int> System::molecules(int* count) const {
  const size_t n = atoms.size();
  std::vector<int> m(n, -1);
  int k = 0;
  if (has_mol) {
    std::unordered_map<int64_t, int> seen;
    for (size_t i = 0; i < n; ++i) {
      auto [it, fresh] = seen.try_emplace(atoms[i].mol, k);
      if (fresh) ++k;
      m[i] = it->second;
    }
  } else {
    std::vector<uint32_t> parent(n);
    std::iota(parent.begin(), parent.end(), 0u);
    auto find = [&](uint32_t x) {
      while (parent[x] != x) x = parent[x] = parent[parent[x]];
      return x;
    };
    for (const auto& b : bonds) parent[find(b.i)] = find(b.j);
    std::vector<int> label(n, -1);
    for (size_t i = 0; i < n; ++i) {
      const uint32_t r = find(static_cast<uint32_t>(i));
      if (label[r] < 0) label[r] = k++;
      m[i] = label[r];
    }
  }
  if (count) *count = k;
  return m;
}

System Trajectory::frame(size_t k) const {
  System s = topology;
  if (k < positions.size()) {
    for (size_t i = 0; i < s.atoms.size() && i < positions[k].size(); ++i) s.atoms[i].pos = positions[k][i];
    if (k < cells.size()) s.cell = cells[k];
    if (k < timesteps.size()) s.timestep = timesteps[k];
  }
  return s;
}

}  // namespace caps
