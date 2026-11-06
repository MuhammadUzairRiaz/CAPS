#include "caps/cell_fit.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

#include "caps/elements.hpp"
#include "caps/field.hpp"
#include "caps/spacegroup.hpp"

namespace caps {

namespace {
constexpr double kTolerance = 0.45;   // Å over the covalent radii, as perceive_bonds

struct Frac {
  Vec3 bc, ca, ab;
  double v;
  explicit Frac(const Cell& c) : bc(cross(c.b, c.c)), ca(cross(c.c, c.a)), ab(cross(c.a, c.b)), v(dot(c.a, bc)) {}
  std::array<double, 3> of(const Vec3& d) const { return {dot(d, bc) / v, dot(d, ca) / v, dot(d, ab) / v}; }
};
}  // namespace

CellFit cell_fit(const System& s, const ForceField* field) {
  CellFit r;
  const Cell& c = s.cell;
  if (s.bonds.empty() || !c.valid() || !(c.periodic[0] || c.periodic[1] || c.periodic[2])) return r;
  const Frac fr(c);
  const Vec3 L[3] = {c.a, c.b, c.c};
  // the image the engines take: fractional separations rounded, axis by axis
  auto nearest = [&](Vec3 d) {
    const auto f = fr.of(d);
    for (int k = 0; k < 3; ++k)
      if (c.periodic[size_t(k)]) d = d - L[k] * std::round(f[size_t(k)]);
    return d;
  };
  std::vector<std::vector<std::pair<uint32_t, Vec3>>> adj(s.atoms.size());
  auto images = [&](auto&& visit) {
    for (int x = c.periodic[0] ? -1 : 0; x <= (c.periodic[0] ? 1 : 0); ++x)
      for (int y = c.periodic[1] ? -1 : 0; y <= (c.periodic[1] ? 1 : 0); ++y)
        for (int z = c.periodic[2] ? -1 : 0; z <= (c.periodic[2] ? 1 : 0); ++z)
          if (x || y || z) visit(x, y, z, L[0] * double(x) + L[1] * double(y) + L[2] * double(z));
  };
  auto mark_images = [&](int x, int y, int z) {
    if (x) r.axis[0] = true;
    if (y) r.axis[1] = true;
    if (z) r.axis[2] = true;
  };
  for (const auto& b : s.bonds) {
    if (b.i >= s.atoms.size() || b.j >= s.atoms.size() || b.i == b.j) continue;
    const Vec3 v = nearest(s.atoms[b.j].pos - s.atoms[b.i].pos);
    adj[b.i].push_back({b.j, v});
    adj[b.j].push_back({b.i, Vec3{0, 0, 0} - v});
    const int ei = s.atoms[b.i].element, ej = s.atoms[b.j].element;
    if (ei <= 0 || ej <= 0) continue;   // a coarse-grained bond has no covalent length to compare with
    const double lim = element(ei).covalent + element(ej).covalent + kTolerance;
    bool twice = false;
    images([&](int x, int y, int z, const Vec3& shift) {
      const Vec3 d = v + shift;
      if (dot(d, d) < lim * lim) twice = true, mark_images(x, y, z);
    });
    r.double_images += twice;
  }
  // an atom bonded to its own image (a chain with one atom per repeat)
  for (size_t i = 0; i < s.atoms.size(); ++i) {
    const int e = s.atoms[i].element;
    if (adj[i].empty() || e <= 1) continue;
    const double lim = 2 * element(e).covalent + kTolerance;
    bool self = false;
    images([&](int x, int y, int z, const Vec3& shift) {
      if (dot(shift, shift) < lim * lim) self = true, mark_images(x, y, z);
    });
    r.double_images += self;
  }
  // a path of bonds ends at one image of its last atom: the engines must read that same one
  auto check = [&](uint32_t a, uint32_t b, const Vec3& p) {
    if (s.atoms[a].element <= 0 || s.atoms[b].element <= 0) return;   // beads: a model's own exclusions, often 1-2 only
    if (a != b && field && a < field->excluded.size() && !std::binary_search(field->excluded[a].begin(), field->excluded[a].end(), b)) return;
    if (a == b) {   // back at the first atom: a real ring closes, a chain through the cell has moved on a period
      if (dot(p, p) < 0.01) return;
      const auto f = fr.of(p);
      for (int k = 0; k < 3; ++k) if (c.periodic[size_t(k)] && std::fabs(f[size_t(k)]) > 1e-3) r.axis[size_t(k)] = true;
      ++r.wrong_images;
      return;
    }
    const auto f = fr.of(p);
    bool wrong = false;
    for (int k = 0; k < 3; ++k)
      if (c.periodic[size_t(k)] && std::fabs(f[size_t(k)]) >= 0.5 - 1e-6) r.axis[size_t(k)] = true, wrong = true;
    r.wrong_images += wrong;
  };
  for (uint32_t j = 0; j < adj.size(); ++j) {
    const auto& nj = adj[j];
    for (size_t x = 0; x < nj.size(); ++x)
      for (size_t y = x + 1; y < nj.size(); ++y) check(nj[x].first, nj[y].first, nj[y].second - nj[x].second);   // 1-3
    for (const auto& [k, jk] : nj) {   // 1-4, each central bond once
      if (k < j) continue;
      for (const auto& [a, ja] : nj) {
        if (a == k) continue;
        for (const auto& [b, kb] : adj[k]) {
          if (b == j) continue;
          check(a, b, jk + kb - ja);
        }
      }
    }
  }
  r.fits = r.double_images == 0 && r.wrong_images == 0;
  return r;
}

std::array<int, 3> copies_to_fit(const System& s, size_t max_atoms) {
  std::array<int, 3> n{1, 1, 1};
  CellFit f = cell_fit(s);
  for (int round = 0; !f.fits && round < 64; ++round) {
    for (size_t k = 0; k < 3; ++k) n[k] += f.axis[k] ? 1 : 0;
    if (size_t(n[0]) * size_t(n[1]) * size_t(n[2]) * s.atoms.size() > max_atoms)
      throw std::invalid_argument("the copies this cell needs for its bonds exceed " + std::to_string(max_atoms) + " atoms");
    f = cell_fit(supercell(s, n[0], n[1], n[2]));
  }
  if (!f.fits) throw std::invalid_argument("no repetition of this cell gives its bonds a unique nearest image");
  return n;
}

std::string cell_fit_message(const System& s, const ForceField* field) {
  const CellFit f = cell_fit(s, field);
  if (f.fits) return "";
  std::string t = "this cell is too small for its bonds (";
  if (f.double_images) t += std::to_string(f.double_images) + " bonded to two images of one neighbour" + (f.wrong_images ? ", " : "");
  if (f.wrong_images) t += std::to_string(f.wrong_images) + " partners within three bonds read at the wrong image";
  t += "), so the force field would act on the wrong atoms: ";
  try {
    const auto n = copies_to_fit(s);
    t += "make a supercell of at least " + std::to_string(n[0]) + "×" + std::to_string(n[1]) + "×" + std::to_string(n[2]) + " first";
  } catch (const std::exception&) {
    t += "use a larger cell";
  }
  return t;
}

}  // namespace caps
