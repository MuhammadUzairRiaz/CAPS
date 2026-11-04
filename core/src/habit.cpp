// CAPS crystal habit: BFDH, attachment energies, surface energies, Wulff shapes and particles (see habit.hpp).
#include "caps/habit.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>

#include "caps/spacegroup.hpp"

namespace caps {

namespace {
constexpr double kJm2 = 4184.0 / 6.02214076e23 / 1e-20;   // J/m² per kcal/mol/Å²
using I3 = std::array<int, 3>;
using M3i = std::array<std::array<long, 3>, 3>;

struct Op {
  int R[3][3];
  double t[3];
};

// h·R (row vector times matrix)
I3 hR(const I3& h, const Op& op) {
  I3 r{0, 0, 0};
  for (int j = 0; j < 3; ++j)
    for (int i = 0; i < 3; ++i) r[size_t(j)] += h[size_t(i)] * op.R[i][j];
  return r;
}

bool absent(const I3& h, const std::vector<Op>& ops) {
  for (const auto& op : ops) {
    if (hR(h, op) != h) continue;
    double ph = 0;
    for (int i = 0; i < 3; ++i) ph += h[size_t(i)] * op.t[i];
    if (std::fabs(ph - std::round(ph)) > 1e-6) return true;
  }
  return false;
}

// reciprocal vectors without 2π: b_i · a_j = δ_ij
void reciprocal(const Cell& c, Vec3 r[3]) {
  const double V = dot(c.a, cross(c.b, c.c));
  r[0] = cross(c.b, c.c) * (1 / V), r[1] = cross(c.c, c.a) * (1 / V), r[2] = cross(c.a, c.b) * (1 / V);
}
Vec3 gvec(const Cell& c, const I3& h) {
  Vec3 r[3];
  reciprocal(c, r);
  return r[0] * double(h[0]) + r[1] * double(h[1]) + r[2] * double(h[2]);
}

// A bond with its lattice image: atom j at f_j + L is bonded to atom i. A short cell bonds an atom to two images of one
// neighbour (a chain along a 2.5 Å axis), which the minimum image would merge into one bond.
struct Edge {
  size_t i, j;
  I3 L;
};
std::vector<Edge> image_bonds(const System& s, std::vector<Vec3>& f) {
  const size_t N = s.atoms.size();
  f.resize(N);
  for (size_t i = 0; i < N; ++i) f[i] = s.cell.to_fractional(s.atoms[i].pos);
  std::set<std::tuple<size_t, size_t, int, int, int>> have;
  std::vector<Edge> e;
  auto add = [&](size_t i, size_t j, I3 L) {
    if (have.insert({i, j, L[0], L[1], L[2]}).second) e.push_back({i, j, L});
  };
  for (const auto& b : s.bonds) {
    double l0 = 1e300;
    for (int x = -2; x <= 2; ++x)
      for (int y = -2; y <= 2; ++y)
        for (int z = -2; z <= 2; ++z) {
          if (b.i == b.j && !x && !y && !z) continue;
          l0 = std::min(l0, norm(s.cell.to_cartesian(f[b.j] + Vec3{double(x), double(y), double(z)} - f[b.i])));
        }
    for (int x = -2; x <= 2; ++x)
      for (int y = -2; y <= 2; ++y)
        for (int z = -2; z <= 2; ++z) {
          if (b.i == b.j && !x && !y && !z) continue;
          if (norm(s.cell.to_cartesian(f[b.j] + Vec3{double(x), double(y), double(z)} - f[b.i])) > l0 + 0.05) continue;
          add(b.i, b.j, {x, y, z}), add(b.j, b.i, {-x, -y, -z});
        }
  }
  return e;
}

// A growth unit: a molecule (or an ion with its shell); its atoms' offsets from its centre of mass.
struct Unit {
  std::vector<size_t> atoms;
  std::vector<Vec3> rel;             // Cartesian, from the centre
  Vec3 com{0, 0, 0};                 // fractional, not wrapped
  std::vector<I3> periods;           // lattice periods of an endless unit (chain, layer, network)
};

struct Topology {
  std::vector<Vec3> f;               // fractional positions
  std::vector<Edge> edges;
  std::vector<std::vector<size_t>> out;   // per atom: its edges
  std::vector<Unit> units;
};

Topology topology_of(const System& s) {
  Topology t;
  t.edges = image_bonds(s, t.f);
  const size_t N = s.atoms.size();
  t.out.resize(N);
  for (size_t k = 0; k < t.edges.size(); ++k) t.out[t.edges[k].i].push_back(k);
  std::vector<Vec3> uf(N);
  std::vector<int> seen(N, -1);
  for (size_t a = 0; a < N; ++a) {
    if (seen[a] >= 0) continue;
    Unit u;
    std::vector<size_t> stack{a};
    seen[a] = int(t.units.size());
    uf[a] = t.f[a];
    std::set<I3> per;
    while (!stack.empty()) {
      const size_t i = stack.back();
      stack.pop_back();
      u.atoms.push_back(i);
      for (size_t k : t.out[i]) {
        const Edge& e = t.edges[k];
        const Vec3 want = uf[i] + (t.f[e.j] + Vec3{double(e.L[0]), double(e.L[1]), double(e.L[2])} - t.f[i]);
        if (seen[e.j] < 0) {
          seen[e.j] = seen[a];
          uf[e.j] = want;
          stack.push_back(e.j);
        } else {
          const Vec3 gap = uf[e.j] - want;
          I3 p{int(std::lround(gap[0])), int(std::lround(gap[1])), int(std::lround(gap[2]))};
          if (p != I3{0, 0, 0}) {
            if (p < I3{0, 0, 0}) p = {-p[0], -p[1], -p[2]};
            per.insert(p);
          }
        }
      }
    }
    double m = 0;
    Vec3 c{0, 0, 0};
    for (size_t i : u.atoms) {
      const double w = std::max(1e-3, s.mass_of(s.atoms[i]));
      c = c + uf[i] * w, m += w;
    }
    u.com = c * (1 / m);
    const Vec3 cc = s.cell.to_cartesian(u.com);
    for (size_t i : u.atoms) u.rel.push_back(s.cell.to_cartesian(uf[i]) - cc);
    u.periods.assign(per.begin(), per.end());
    t.units.push_back(std::move(u));
  }
  return t;
}

// The bonds of a structure built from copies of the cell's atoms (origin: each atom's index in the cell): every bond of the
// cell joins an atom to the copy of its partner at the bond's place, periodically in the structure's own cell.
void join_bonds(System& r, const std::vector<size_t>& origin, const System& cell, const Topology& t) {
  std::vector<std::vector<size_t>> of(cell.atoms.size());
  for (size_t k = 0; k < origin.size(); ++k) of[origin[k]].push_back(k);
  r.bonds.clear();
  for (size_t k = 0; k < origin.size(); ++k)
    for (size_t q : t.out[origin[k]]) {
      const Edge& e = t.edges[q];
      const Vec3 want = r.atoms[k].pos + cell.cell.to_cartesian(t.f[e.j] + Vec3{double(e.L[0]), double(e.L[1]), double(e.L[2])} - t.f[e.i]);
      for (size_t m : of[e.j]) {
        const Vec3 d = r.cell.valid() ? r.cell.minimum_image(r.atoms[m].pos - want) : r.atoms[m].pos - want;
        if (norm(d) < 0.1) {
          if (k < m) r.bonds.push_back({uint32_t(k), uint32_t(m), 1});
          break;
        }
      }
    }
  r.bonds_from_file = true;
}

// unimodular M with h·M = (0, 0, 1): columns 0, 1 span the plane's lattice, column 2 = w with h·w = 1
M3i plane_basis(const I3& h) {
  M3i M{};
  for (int i = 0; i < 3; ++i) M[size_t(i)][size_t(i)] = 1;
  long r[3] = {h[0], h[1], h[2]};
  auto colop = [&](int dst, int src, long q) {   // column dst -= q · column src
    for (int i = 0; i < 3; ++i) M[size_t(i)][size_t(dst)] -= q * M[size_t(i)][size_t(src)];
    r[dst] -= q * r[src];
  };
  auto swapc = [&](int a, int b) {
    for (int i = 0; i < 3; ++i) std::swap(M[size_t(i)][size_t(a)], M[size_t(i)][size_t(b)]);
    std::swap(r[a], r[b]);
  };
  // Euclid on the entries until only r[2] is non-zero
  for (int guard = 0; guard < 1000; ++guard) {
    int nz = 0, piv = -1;
    for (int c = 0; c < 3; ++c)
      if (r[c] != 0) {
        ++nz;
        if (piv < 0 || std::labs(r[c]) < std::labs(r[piv])) piv = c;
      }
    if (nz <= 1) {
      if (piv != 2) swapc(piv, 2);
      break;
    }
    for (int c = 0; c < 3; ++c)
      if (c != piv && r[c] != 0) colop(c, piv, r[c] / r[piv]);
  }
  if (r[2] < 0) {
    for (int i = 0; i < 3; ++i) M[size_t(i)][2] = -M[size_t(i)][2];
    r[2] = -r[2];
  }
  if (r[2] != 1) throw std::invalid_argument("plane_basis: the indices have a common divisor");
  return M;
}
Vec3 latvec(const Cell& c, long x, long y, long z) { return c.a * double(x) + c.b * double(y) + c.c * double(z); }

// the shortest in-plane basis (Lagrange–Gauss reduction)
void reduce2(Vec3& A, Vec3& B) {
  for (int it = 0; it < 100; ++it) {
    if (dot(A, A) > dot(B, B)) std::swap(A, B);
    const double q = std::round(dot(A, B) / dot(A, A));
    if (q == 0) break;
    B = B - A * q;
  }
}

// a cell's widths (perpendicular heights)
std::array<double, 3> widths(const Cell& c) {
  const double V = std::fabs(dot(c.a, cross(c.b, c.c)));
  return {V / norm(cross(c.b, c.c)), V / norm(cross(c.c, c.a)), V / norm(cross(c.a, c.b))};
}
int repeats(double width, double cutoff) { return std::max(1, int(std::ceil(2.02 * cutoff / width - 1e-9))); }

struct Built {
  System sys;
  std::vector<size_t> origin;
  size_t units = 0;      // growth units in it
  double area = 0;       // slab: its surface area, Å²
};

// the bulk repeated until each width is at least twice the cut-off
Built make_bulk(const System& cell, const Topology& t, double cutoff, size_t max_atoms) {
  const auto wd = widths(cell.cell);
  const int n0 = repeats(wd[0], cutoff), n1 = repeats(wd[1], cutoff), n2 = repeats(wd[2], cutoff);
  if (size_t(n0) * size_t(n1) * size_t(n2) * cell.atoms.size() > max_atoms)
    throw std::runtime_error("the bulk repeated to twice the cut-off would exceed " + std::to_string(max_atoms) + " atoms");
  Built b;
  b.sys.cell.a = cell.cell.a * double(n0), b.sys.cell.b = cell.cell.b * double(n1), b.sys.cell.c = cell.cell.c * double(n2);
  for (int i = 0; i < n0; ++i)
    for (int j = 0; j < n1; ++j)
      for (int k = 0; k < n2; ++k)
        for (size_t a = 0; a < cell.atoms.size(); ++a) {
          Atom at = cell.atoms[a];
          at.pos = cell.cell.to_cartesian(t.f[a] + Vec3{double(i), double(j), double(k)});
          at.image = {0, 0, 0};
          b.sys.atoms.push_back(at), b.origin.push_back(a);
        }
  b.units = t.units.size() * size_t(n0 * n1 * n2);
  b.sys.has_charges = cell.has_charges, b.sys.types = cell.types;
  join_bonds(b.sys, b.origin, cell, t);
  return b;
}

// the units whose centre lies at s ∈ [s0, s0 + width) (s = h·f, one lattice plane spacing per unit of s), whole, in a cell
// with the plane's lattice in-plane (repeated to twice the cut-off) and c along the normal (thickness width·d plus vacuum)
Built make_slab(const System& cell, const Topology& t, const I3& h, double s0, double width, double vacuum, double cutoff, size_t max_atoms) {
  const M3i M = plane_basis(h);
  Vec3 A = latvec(cell.cell, M[0][0], M[1][0], M[2][0]), B = latvec(cell.cell, M[0][1], M[1][1], M[2][1]);
  reduce2(A, B);
  const Vec3 w{double(M[0][2]), double(M[1][2]), double(M[2][2])};
  const Vec3 G = gvec(cell.cell, h);
  const double d = 1 / norm(G);
  const Vec3 n = G * d;
  if (dot(cross(A, B), n) < 0) std::swap(A, B);
  Cell one;
  one.a = A, one.b = B, one.c = n * (width * d + vacuum);
  const auto wd = widths(one);
  const int na = repeats(wd[0], cutoff), nb = repeats(wd[1], cutoff);
  Built sl;
  sl.sys.cell.a = A * double(na), sl.sys.cell.b = B * double(nb), sl.sys.cell.c = one.c;
  sl.area = norm(cross(A, B)) * na * nb;
  int64_t mol = 0;
  for (const auto& u : t.units) {
    const double su = double(h[0]) * u.com[0] + double(h[1]) * u.com[1] + double(h[2]) * u.com[2];
    for (long j = long(std::floor(s0 - su)) - 1; j <= long(std::ceil(s0 + width - su)) + 1; ++j) {
      const double s = su + double(j);
      if (s < s0 - 1e-9 || s >= s0 + width - 1e-9) continue;
      const Vec3 c = cell.cell.to_cartesian(u.com + w * double(j)) - n * (s0 * d);
      for (int ia = 0; ia < na; ++ia)
        for (int ib = 0; ib < nb; ++ib) {
          ++mol;
          for (size_t q = 0; q < u.atoms.size(); ++q) {
            Atom at = cell.atoms[u.atoms[q]];
            at.pos = c + u.rel[q] + A * double(ia) + B * double(ib);
            at.mol = mol;
            at.image = {0, 0, 0};
            sl.sys.atoms.push_back(at), sl.origin.push_back(u.atoms[q]);
          }
          ++sl.units;
        }
    }
  }
  if (sl.sys.atoms.size() > max_atoms) throw std::runtime_error("the slab would exceed " + std::to_string(max_atoms) + " atoms");
  sl.sys.has_mol = true;
  sl.sys.has_charges = cell.has_charges, sl.sys.types = cell.types;
  join_bonds(sl.sys, sl.origin, cell, t);
  return sl;
}

struct Energy {
  double E = 0;       // total
  double dipole = 0;  // along c, e·Å
  double charge = 0;
};
Energy energy_of(const System& s, const HabitOptions& o, bool slab) {
  const ForceField F = o.field_for(s);
  if (F.type_index.size() != s.atoms.size()) throw std::runtime_error("the force field does not describe the structure");
  Evaluator ev(F, o.energy);
  std::vector<double> x(3 * s.atoms.size()), f;
  for (size_t i = 0; i < s.atoms.size(); ++i)
    for (int k = 0; k < 3; ++k) x[3 * i + size_t(k)] = s.atoms[i].pos[size_t(k)];
  Energy e;
  e.E = ev.compute(x, s.cell, f).total();
  if (slab) {
    const Vec3 nz = s.cell.c * (1 / norm(s.cell.c));
    double z0 = 0;
    for (const auto& a : s.atoms) z0 += dot(a.pos, nz);
    z0 /= double(s.atoms.size());
    for (size_t i = 0; i < s.atoms.size(); ++i) e.dipole += F.charge[i] * (dot(s.atoms[i].pos, nz) - z0), e.charge += F.charge[i];
  }
  if (!std::isfinite(e.E)) throw std::runtime_error("the energy is not finite");
  return e;
}

std::string hkl_text(const I3& h) {
  char b[64];
  std::snprintf(b, sizeof b, "(%d %d %d)", h[0], h[1], h[2]);
  return b;
}
}  // namespace

WulffShape wulff_shape(const std::vector<Vec3>& normals, const std::vector<double>& distances) {
  WulffShape w;
  if (normals.size() != distances.size() || normals.empty()) return w;
  double hmax = 0;
  for (double h : distances) {
    if (!(h > 0)) return w;
    hmax = std::max(hmax, h);
  }
  const double L = 1e3 * hmax, eps = 1e-10 * hmax;
  // a large cube, then clipped by every plane in turn
  std::vector<std::vector<Vec3>> faces;
  std::vector<int> plane;
  for (int ax = 0; ax < 3; ++ax)
    for (int sg : {-1, 1}) {
      const int u = (ax + 1) % 3, v = (ax + 2) % 3;
      std::vector<Vec3> p;
      for (auto [a, b] : std::array<std::pair<int, int>, 4>{{{-1, -1}, {1, -1}, {1, 1}, {-1, 1}}}) {
        Vec3 x{0, 0, 0};
        x[size_t(ax)] = sg * L, x[size_t(u)] = a * L, x[size_t(v)] = b * L;
        p.push_back(x);
      }
      faces.push_back(p), plane.push_back(-1);
    }
  for (size_t pi = 0; pi < normals.size(); ++pi) {
    const Vec3 n = normals[pi];
    const double h = distances[pi];
    std::vector<std::vector<Vec3>> nf;
    std::vector<int> np;
    std::vector<Vec3> cut;
    for (size_t fi = 0; fi < faces.size(); ++fi) {
      const auto& P = faces[fi];
      std::vector<Vec3> Q;
      for (size_t k = 0; k < P.size(); ++k) {
        const Vec3 a = P[k], b = P[(k + 1) % P.size()];
        const double da = dot(n, a) - h, db = dot(n, b) - h;
        if (da <= eps) {
          Q.push_back(a);
          if (da > -eps) cut.push_back(a);
        }
        if ((da < -eps && db > eps) || (da > eps && db < -eps)) {
          const Vec3 x = a + (b - a) * (da / (da - db));
          Q.push_back(x), cut.push_back(x);
        }
      }
      if (Q.size() >= 3) nf.push_back(Q), np.push_back(plane[fi]);
    }
    // the new face: the cut points, merged and ordered around their centre (later planes clip it as a convex polygon)
    std::vector<Vec3> poly;
    for (const auto& x : cut) {
      bool dup = false;
      for (const auto& y : poly) dup = dup || norm(x - y) < 1e-9 * hmax;
      if (!dup) poly.push_back(x);
    }
    if (poly.size() >= 3) {
      Vec3 c{0, 0, 0};
      for (const auto& x : poly) c = c + x;
      c = c * (1.0 / double(poly.size()));
      Vec3 e1 = poly[0] - c;
      e1 = e1 * (1 / norm(e1));
      const Vec3 e2 = cross(n, e1);
      std::sort(poly.begin(), poly.end(), [&](const Vec3& a, const Vec3& b) {
        return std::atan2(dot(a - c, e2), dot(a - c, e1)) < std::atan2(dot(b - c, e2), dot(b - c, e1));
      });
      nf.push_back(poly), np.push_back(int(pi));
    }
    faces.swap(nf), plane.swap(np);
  }
  for (int p : plane)
    if (p < 0) return w;   // unbounded
  // merge the vertices, order each facet counter-clockwise seen from outside, drop the degenerate ones
  auto vid = [&](const Vec3& x) {
    for (size_t i = 0; i < w.vertices.size(); ++i)
      if (norm(w.vertices[i] - x) < 1e-7 * hmax) return int(i);
    w.vertices.push_back(x);
    return int(w.vertices.size() - 1);
  };
  for (size_t fi = 0; fi < faces.size(); ++fi) {
    std::vector<int> ids;
    for (const auto& x : faces[fi]) {
      const int id = vid(x);
      if (std::find(ids.begin(), ids.end(), id) == ids.end()) ids.push_back(id);
    }
    if (ids.size() < 3) continue;
    const Vec3 n = normals[size_t(plane[fi])];
    Vec3 c{0, 0, 0};
    for (int id : ids) c = c + w.vertices[size_t(id)];
    c = c * (1.0 / double(ids.size()));
    Vec3 e1 = w.vertices[size_t(ids[0])] - c;
    e1 = e1 * (1 / norm(e1));
    const Vec3 e2 = cross(n, e1);
    std::sort(ids.begin(), ids.end(), [&](int a, int b) {
      const Vec3 pa = w.vertices[size_t(a)] - c, pb = w.vertices[size_t(b)] - c;
      return std::atan2(dot(pa, e2), dot(pa, e1)) < std::atan2(dot(pb, e2), dot(pb, e1));
    });
    Vec3 s{0, 0, 0};
    for (size_t k = 0; k < ids.size(); ++k) s = s + cross(w.vertices[size_t(ids[k])], w.vertices[size_t(ids[(k + 1) % ids.size()])]);
    const double area = 0.5 * std::fabs(dot(s, n));
    if (area < 1e-12 * hmax * hmax) continue;
    w.facets.push_back(ids), w.facet_plane.push_back(plane[fi]), w.facet_area.push_back(area);
    w.area += area, w.volume += area * distances[size_t(plane[fi])] / 3;
  }
  // the widths across the facets' own normals: largest / smallest
  if (!w.facets.empty()) {
    double lo = 1e300, hi = 0;
    for (int pl : w.facet_plane) {
      const Vec3& n = normals[size_t(pl)];
      double a = 1e300, b = -1e300;
      for (const auto& v : w.vertices) a = std::min(a, dot(v, n)), b = std::max(b, dot(v, n));
      lo = std::min(lo, b - a), hi = std::max(hi, b - a);
    }
    w.aspect = hi / std::max(1e-12, lo);
  }
  return w;
}

void habit_planes(const System& cell, const HabitResult& r, const std::string& model, std::vector<Vec3>& normals, std::vector<double>& distances,
                  std::vector<int>* face) {
  normals.clear(), distances.clear();
  if (face) face->clear();
  for (size_t fi = 0; fi < r.faces.size(); ++fi) {
    const HabitFace& f = r.faces[fi];
    double h = 0;
    if (model == "bfdh") h = 1 / f.d_eff;
    else if (model == "growth") h = std::fabs(f.e_att);
    else if (model == "equilibrium") h = f.gamma;
    else throw std::invalid_argument("shape model '" + model + "': bfdh, growth or equilibrium");
    if (!std::isfinite(h) || h <= 0) continue;
    for (const auto& m : f.members) {
      const Vec3 G = gvec(cell.cell, m);
      normals.push_back(G * (1 / norm(G))), distances.push_back(h);
      if (face) face->push_back(int(fi));
    }
  }
}

HabitResult crystal_habit(const System& cell, const HabitOptions& o) {
  HabitResult r;
  if (!cell.cell.valid()) throw std::invalid_argument("habit: the crystal needs a periodic cell");
  if (cell.atoms.empty()) throw std::invalid_argument("habit: no atoms");
  if (o.max_index < 1 || o.max_index > 8) throw std::invalid_argument("habit: the largest index must be 1 to 8");
  auto stop = [&](const std::string& w, double f) {
    if (o.progress && !o.progress(w, f)) throw std::runtime_error("habit stopped");
  };
  stop("symmetry", 0);
  // ---- symmetry operations (point group and extinctions)
  std::vector<Op> ops;
  try {
    const SymmetryFound sf = find_symmetry(cell, o.symmetry_tolerance);
    if (const SpaceGroupSetting* g = find_space_group(sf.key); g && sf.number > 0) {
      r.space_group = sf.hm;
      for (const auto& so : hall_operations(g->hall)) {
        Op op;
        bool integral = true;
        for (int i = 0; i < 3; ++i)
          for (int j = 0; j < 3; ++j) {
            op.R[i][j] = int(std::lround(so.R[i][j]));
            integral = integral && std::fabs(so.R[i][j] - op.R[i][j]) < 1e-9;
          }
        for (int i = 0; i < 3; ++i) op.t[i] = so.t[i];
        if (integral) ops.push_back(op);
      }
    }
  } catch (const std::exception&) {
  }
  if (ops.empty()) {
    Op id{};
    for (int i = 0; i < 3; ++i) id.R[i][i] = 1;
    ops.push_back(id);
    r.space_group = "P 1";
    r.notes.push_back("no space group found: every face is its own form and no reflection is extinguished");
  }
  r.operations = int(ops.size());
  // ---- forms and effective spacings
  std::set<I3> done;
  const int M = o.max_index;
  for (int h = -M; h <= M; ++h)
    for (int k = -M; k <= M; ++k)
      for (int l = -M; l <= M; ++l) {
        const I3 hk{h, k, l};
        if (hk == I3{0, 0, 0} || done.count(hk)) continue;
        if (std::gcd(std::gcd(std::abs(h), std::abs(k)), std::abs(l)) != 1) continue;
        HabitFace f;
        std::set<I3> mem;
        for (const auto& op : ops) mem.insert(hR(hk, op));
        for (const auto& m : mem) done.insert(m);
        f.members.assign(mem.begin(), mem.end());
        // the representative: most non-negative indices, then the largest
        f.hkl = *std::max_element(f.members.begin(), f.members.end(), [](const I3& a, const I3& b) {
          auto key = [](const I3& x) { return std::make_tuple(int(x[0] >= 0) + int(x[1] >= 0) + int(x[2] >= 0), x[0], x[1], x[2]); };
          return key(a) < key(b);
        });
        f.d = 1 / norm(gvec(cell.cell, hk));
        f.order = 1;
        while (f.order < 12 && absent({hk[0] * f.order, hk[1] * f.order, hk[2] * f.order}, ops)) ++f.order;
        f.d_eff = f.d / f.order;
        if (f.d_eff < o.min_d) continue;
        r.faces.push_back(f);
      }
  std::stable_sort(r.faces.begin(), r.faces.end(), [](const HabitFace& a, const HabitFace& b) { return a.d_eff > b.d_eff + 1e-9; });
  for (size_t i = 0; i < r.faces.size(); ++i) r.faces[i].rank = int(i + 1);
  // ---- energies
  const Topology topo = topology_of(cell);
  const std::vector<Unit>& units = topo.units;
  r.units = int(units.size());
  if (o.attachment || o.surface) {
    if (!o.field_for) throw std::invalid_argument("habit: energies need a force field");
    stop("bulk", 0.02);
    const Built bb = make_bulk(cell, topo, o.energy.cutoff, o.max_atoms);
    const double e_bulk = energy_of(bb.sys, o, false).E / double(bb.units);
    r.e_latt = e_bulk;
    const double vac = o.vacuum > 0 ? o.vacuum : 2 * o.energy.cutoff + 10;
    const size_t nf = std::min(r.faces.size(), size_t(std::max(0, o.forms)));
    for (size_t fi = 0; fi < nf; ++fi) {
      HabitFace& f = r.faces[fi];
      stop("face " + hkl_text(f.hkl), 0.05 + 0.95 * double(fi) / double(nf));
      // a unit endless across the face would be cut
      bool cuts = false;
      for (const auto& u : units)
        for (const auto& p : u.periods)
          cuts = cuts || (f.hkl[0] * p[0] + f.hkl[1] * p[1] + f.hkl[2] * p[2]) != 0;
      if (cuts) {
        f.note = "the face cuts a covalent chain or network";
        continue;
      }
      // the cuts between units, within one slice period 1/order of s
      const double per = 1.0 / f.order;
      std::vector<double> sv;
      for (const auto& u : units) {
        double s = f.hkl[0] * u.com[0] + f.hkl[1] * u.com[1] + f.hkl[2] * u.com[2];
        s = s / per;
        s = (s - std::floor(s)) * per;
        sv.push_back(s);
      }
      std::sort(sv.begin(), sv.end());
      std::vector<double> uniq;
      for (double s : sv)
        if (uniq.empty() || s - uniq.back() > 1e-6) uniq.push_back(s);
      if (uniq.size() > 1 && uniq.front() + per - uniq.back() < 1e-6) uniq.pop_back();
      std::vector<double> cutsAt;
      for (size_t q = 0; q < uniq.size(); ++q) {
        const double a = uniq[q], b = q + 1 < uniq.size() ? uniq[q + 1] : uniq[0] + per;
        cutsAt.push_back(0.5 * (a + b));
      }
      double best = std::numeric_limits<double>::infinity();
      std::string why;
      for (double s0 : cutsAt) {
        try {
          const Built sl = make_slab(cell, topo, f.hkl, s0, per, vac, o.energy.cutoff, o.max_atoms);
          if (sl.units == 0) continue;
          const Energy e = energy_of(sl.sys, o, true);
          if (std::fabs(e.charge) > 1e-6 * double(sl.units)) { why = "every slice carries a net charge"; continue; }
          if (std::fabs(e.dipole) / sl.area > 1e-4) { why = "every slice carries a dipole along the normal (a polar face)"; continue; }
          const double e_slice = e.E / double(sl.units);
          const double e_att = e_bulk - e_slice;
          if (std::fabs(e_att) < std::fabs(best)) best = e_att, f.e_att = e_att, f.e_slice = e_slice, f.cut = s0;
        } catch (const std::runtime_error& e) {
          if (std::string(e.what()) == "habit stopped") throw;
          why = e.what();
        }
      }
      if (!std::isfinite(f.e_att)) {
        f.note = why.empty() ? "no slice could be built" : why;
        continue;
      }
      if (o.surface) {
        const int K = std::max(2, int(std::ceil(o.energy.cutoff / f.d_eff)) + 1);
        try {
          const Built sl = make_slab(cell, topo, f.hkl, f.cut, K * per, vac, o.energy.cutoff, o.max_atoms);
          const Energy e = energy_of(sl.sys, o, true);
          f.gamma = (e.E - double(sl.units) * e_bulk) / (2 * sl.area) * kJm2;
        } catch (const std::runtime_error& e) {
          if (std::string(e.what()) == "habit stopped") throw;
          f.note = e.what();
        }
      }
      if (!o.attachment) f.e_att = f.e_slice = std::numeric_limits<double>::quiet_NaN();
    }
    if (nf < r.faces.size())
      r.notes.push_back("energies for the " + std::to_string(nf) + " forms of largest d_eff; the others are left out of the growth and equilibrium shapes");
  }
  // ---- shapes and each face's share of their surfaces
  // a crystal of endless chains along one lattice direction: the faces that cut the chains have no energy, so the shape
  // is a prism open along them; caps across the chain axis close it and the shares are of its sides alone
  std::set<I3> chain_dirs;
  for (const auto& u : units)
    for (const auto& p : u.periods) chain_dirs.insert(p);
  const bool chains = chain_dirs.size() == 1;
  auto shape = [&](const std::string& model, WulffShape& w, std::vector<int>& face, double HabitFace::*share) {
    std::vector<Vec3> n;
    std::vector<double> d;
    habit_planes(cell, r, model, n, d, &face);
    if (n.empty()) return;
    w = wulff_shape(n, d);
    if (!w.ok() && chains && model != "bfdh") {
      const I3 p = *chain_dirs.begin();
      Vec3 ax = latvec(cell.cell, p[0], p[1], p[2]);
      ax = ax * (1 / norm(ax));
      const double hmax = *std::max_element(d.begin(), d.end());
      n.push_back(ax), d.push_back(3 * hmax), face.push_back(-1);
      n.push_back(ax * -1.0), d.push_back(3 * hmax), face.push_back(-1);
      w = wulff_shape(n, d);
      if (w.ok()) {
        r.notes.push_back("the " + model + " shape is a prism open along the chains (faces that cut them have no energy): its sides' shares and cross-section aspect are given");
        double lo = 1e300, hi = 0;   // the cross-section: widths across the sides alone
        for (size_t k = 0; k < w.facets.size(); ++k) {
          if (face[size_t(w.facet_plane[k])] < 0) continue;
          const Vec3& nn = n[size_t(w.facet_plane[k])];
          double a = 1e300, b = -1e300;
          for (const auto& v : w.vertices) a = std::min(a, dot(v, nn)), b = std::max(b, dot(v, nn));
          lo = std::min(lo, b - a), hi = std::max(hi, b - a);
        }
        w.aspect = hi / std::max(1e-12, lo);
      }
    }
    if (!w.ok()) {
      r.notes.push_back("the " + model + " faces do not enclose a shape (raise the largest index or the number of forms)");
      return;
    }
    double sides = 0;
    for (size_t k = 0; k < w.facets.size(); ++k)
      if (face[size_t(w.facet_plane[k])] >= 0) sides += w.facet_area[k];
    for (size_t k = 0; k < w.facets.size(); ++k)
      if (const int fi = face[size_t(w.facet_plane[k])]; fi >= 0) r.faces[size_t(fi)].*share += w.facet_area[k] / sides;
  };
  shape("bfdh", r.bfdh, r.bfdh_face, &HabitFace::area_bfdh);
  if (o.attachment) shape("growth", r.growth, r.growth_face, &HabitFace::area_growth);
  if (o.surface) shape("equilibrium", r.equilibrium, r.equilibrium_face, &HabitFace::area_equilibrium);
  stop("done", 1);
  return r;
}

System carve_particle(const System& cell, const std::vector<Vec3>& normals, const std::vector<double>& distances, double radius) {
  if (!cell.cell.valid()) throw std::invalid_argument("particle: the crystal needs a periodic cell");
  if (radius <= 0) throw std::invalid_argument("particle: the size must be positive");
  const WulffShape w = wulff_shape(normals, distances);
  if (!w.ok()) throw std::invalid_argument("particle: the faces do not enclose a shape");
  double rmax = 0;
  for (const auto& v : w.vertices) rmax = std::max(rmax, norm(v));
  const double scale = radius / rmax;
  const Topology topo = topology_of(cell);
  const std::vector<Unit>& units = topo.units;
  for (const auto& u : units)
    if (!u.periods.empty()) throw std::invalid_argument("particle: the crystal has an endless chain or network; carve it with the nano builder instead");
  // lattice translations whose cell can hold a unit within the radius
  Vec3 rc[3];
  reciprocal(cell.cell, rc);
  int n[3];
  for (int k = 0; k < 3; ++k) n[k] = int(std::ceil((radius + 10) * norm(rc[k]))) + 1;
  System out;
  out.has_mol = true;
  out.has_charges = cell.has_charges;
  out.types = cell.types;
  out.bonds_from_file = true;
  // the shape's centre at the cell's centre of mass of units
  Vec3 c0{0, 0, 0};
  for (const auto& u : units) c0 = c0 + cell.cell.to_cartesian(u.com);
  c0 = c0 * (1.0 / double(units.size()));
  int64_t mol = 0;
  for (int i = -n[0]; i <= n[0]; ++i)
    for (int j = -n[1]; j <= n[1]; ++j)
      for (int k = -n[2]; k <= n[2]; ++k)
        for (const auto& u : units) {
          const Vec3 c = cell.cell.to_cartesian(u.com + Vec3{double(i), double(j), double(k)}) - c0;
          bool in = true;
          for (size_t p = 0; p < normals.size() && in; ++p) in = dot(normals[p], c) <= distances[p] * scale + 1e-9;
          if (!in) continue;
          const size_t base = out.atoms.size();
          std::map<size_t, size_t> local;
          ++mol;
          for (size_t q = 0; q < u.atoms.size(); ++q) {
            Atom at = cell.atoms[u.atoms[q]];
            at.pos = c + u.rel[q];
            at.mol = mol;
            at.id = int64_t(out.atoms.size() + 1);
            at.image = {0, 0, 0};
            local[u.atoms[q]] = base + q;
            out.atoms.push_back(at);
          }
          for (size_t q : u.atoms)
            for (size_t k : topo.out[q])
              if (topo.edges[k].i < topo.edges[k].j && local.count(topo.edges[k].j))
                out.bonds.push_back({uint32_t(local[topo.edges[k].i]), uint32_t(local[topo.edges[k].j]), 1});
          if (out.atoms.size() > 15000) throw std::invalid_argument("particle: more than 15,000 atoms; choose a smaller size");
        }
  return out;
}

Json habit_json(const HabitResult& r) {
  auto val = [](double v) { return std::isfinite(v) ? Json(v) : Json(nullptr); };
  Json j = Json::object(), fs = Json::array();
  for (const auto& f : r.faces) {
    Json x = Json::object(), h = Json::array();
    for (int q : f.hkl) h.push_back(Json(q));
    x["hkl"] = h, x["faces"] = double(f.members.size()), x["d"] = f.d, x["d_eff"] = f.d_eff, x["order"] = double(f.order), x["rank"] = double(f.rank);
    x["e_att"] = val(f.e_att), x["e_slice"] = val(f.e_slice), x["gamma"] = val(f.gamma), x["note"] = f.note;
    x["area_bfdh"] = f.area_bfdh, x["area_growth"] = f.area_growth, x["area_equilibrium"] = f.area_equilibrium;
    fs.push_back(x);
  }
  j["space_group"] = r.space_group, j["operations"] = double(r.operations), j["units"] = double(r.units), j["e_latt"] = val(r.e_latt), j["faces"] = fs;
  Json shapes = Json::object();
  auto shp = [&](const char* name, const WulffShape& w, const std::vector<int>& face) {
    if (!w.ok()) return;
    Json x = Json::object(), vs = Json::array(), fc = Json::array(), ff = Json::array(), fa = Json::array();
    for (const auto& v : w.vertices) {
      Json p = Json::array();
      for (double c : v) p.push_back(c);
      vs.push_back(p);
    }
    for (size_t k = 0; k < w.facets.size(); ++k) {
      Json p = Json::array();
      for (int id : w.facets[k]) p.push_back(Json(id));
      fc.push_back(p);
      ff.push_back(Json(face[size_t(w.facet_plane[k])])), fa.push_back(w.facet_area[k]);
    }
    x["vertices"] = vs, x["facets"] = fc, x["facet_face"] = ff, x["facet_area"] = fa, x["area"] = w.area, x["volume"] = w.volume, x["aspect"] = w.aspect;
    shapes[name] = x;
  };
  shp("bfdh", r.bfdh, r.bfdh_face), shp("growth", r.growth, r.growth_face), shp("equilibrium", r.equilibrium, r.equilibrium_face);
  j["shapes"] = shapes;
  Json nt = Json::array();
  for (const auto& n : r.notes) nt.push_back(n);
  j["notes"] = nt;
  return j;
}

}  // namespace caps
