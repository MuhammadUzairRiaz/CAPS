// CAPS appearance: stereo labels, molecular surfaces, coordination polyhedra, ribbon paths (see caps/appearance.hpp).
#include "caps/appearance.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <unordered_map>

#include "caps/analysis.hpp"
#include "caps/config.hpp"
#include "caps/elements.hpp"
#include "caps/render.hpp"

namespace caps {

namespace {

Vec3 unitv(const Vec3& v) { const double n = norm(v); return n > 1e-12 ? v * (1.0 / n) : Vec3{0, 0, 1}; }

// ---------------------------------------------------------------- CIP

struct Neighbour { uint32_t atom; int order; };

struct Node {
  int z;             // atomic number (0: phantom)
  int atom;          // −1 for duplicates and phantoms
  int parent_atom;   // the atom it was reached from
  std::vector<int> path;   // atoms from the centre down to this node (for ring closures)
};

// Children of a node in the hierarchical digraph: neighbours other than the parent, a duplicate for an atom already on
// the path (ring closure), and duplicates for the extra bonds of double and triple bonds (both ends).
std::vector<Node> children(const System& s, const std::vector<std::vector<Neighbour>>& nb, const Node& n) {
  std::vector<Node> out;
  if (n.atom < 0) return out;
  for (const auto& e : nb[size_t(n.atom)]) {
    const int b = int(e.atom);
    const int dup = e.order == 2 ? 1 : e.order == 3 ? 2 : 0;
    const int zb = s.atoms[size_t(b)].element;
    for (int k = 0; k < dup; ++k) out.push_back({zb, -1, n.atom, {}});
    if (b == n.parent_atom) continue;
    if (std::find(n.path.begin(), n.path.end(), b) != n.path.end()) { out.push_back({zb, -1, n.atom, {}}); continue; }
    Node c{zb, b, n.atom, n.path};
    c.path.push_back(b);
    out.push_back(std::move(c));
  }
  std::sort(out.begin(), out.end(), [](const Node& a, const Node& b) { return a.z > b.z; });
  return out;
}

// −1, 0, +1: branch a ranks below, equal to, above branch b (rule 1a, sphere by sphere).
int compare_branches(const System& s, const std::vector<std::vector<Neighbour>>& nb, int centre, int ra, int rb) {
  const int za = s.atoms[size_t(ra)].element, zb = s.atoms[size_t(rb)].element;
  if (za != zb) return za > zb ? 1 : -1;
  std::vector<Node> la{{za, ra, centre, {centre, ra}}}, lb{{zb, rb, centre, {centre, rb}}};
  for (int depth = 0; depth < 24; ++depth) {
    std::vector<Node> na, nbv;
    std::vector<std::vector<int>> sa, sb;
    for (const auto& n : la) { auto c = children(s, nb, n); std::vector<int> z; for (const auto& x : c) z.push_back(x.z); sa.push_back(z); for (auto& x : c) na.push_back(std::move(x)); }
    for (const auto& n : lb) { auto c = children(s, nb, n); std::vector<int> z; for (const auto& x : c) z.push_back(x.z); sb.push_back(z); for (auto& x : c) nbv.push_back(std::move(x)); }
    // the sets of this sphere, in the order of their parents
    const size_t m = std::max(sa.size(), sb.size());
    for (size_t i = 0; i < m; ++i) {
      std::vector<int> A = i < sa.size() ? sa[i] : std::vector<int>{}, B = i < sb.size() ? sb[i] : std::vector<int>{};
      const size_t k = std::max(A.size(), B.size());
      A.resize(k, 0), B.resize(k, 0);
      if (A != B) return A > B ? 1 : -1;
    }
    if (na.empty() && nbv.empty()) return 0;
    la = std::move(na), lb = std::move(nbv);
    if (la.size() > 4000 || lb.size() > 4000) return 0;   // a large ring system: give up (no label)
  }
  return 0;
}

// ---------------------------------------------------------------- grid distance transform (for the excluded surface)

// 1D squared distance transform (Felzenszwalb & Huttenlocher 2012).
void edt1(const std::vector<double>& f, std::vector<double>& d, std::vector<int>& v, std::vector<double>& z) {
  const int n = int(f.size());
  int k = 0;
  v[0] = 0;
  z[0] = -std::numeric_limits<double>::infinity();
  z[1] = std::numeric_limits<double>::infinity();
  for (int q = 1; q < n; ++q) {
    double sp;
    for (;;) {
      sp = ((f[size_t(q)] + double(q) * q) - (f[size_t(v[size_t(k)])] + double(v[size_t(k)]) * v[size_t(k)])) / (2.0 * q - 2.0 * v[size_t(k)]);
      if (sp <= z[size_t(k)]) { if (--k < 0) { k = 0; break; } } else break;
    }
    ++k;
    v[size_t(k)] = q;
    z[size_t(k)] = sp;
    z[size_t(k) + 1] = std::numeric_limits<double>::infinity();
  }
  k = 0;
  for (int q = 0; q < n; ++q) {
    while (z[size_t(k) + 1] < q) ++k;
    const double dq = q - v[size_t(k)];
    d[size_t(q)] = dq * dq + f[size_t(v[size_t(k)])];
  }
}

}  // namespace

namespace {
// neighbours with bond orders for the CIP digraph (aromatic bonds: one Kekulé double bond each)
std::vector<std::vector<Neighbour>> cip_neighbours(const System& s) {
  std::vector<std::vector<Neighbour>> nb(s.atoms.size());
  for (const auto& b : s.bonds) {
    const int order = b.order == 2 ? 2 : b.order == 3 ? 3 : 1;
    nb[b.i].push_back({b.j, order}), nb[b.j].push_back({b.i, order});
  }
  for (const auto& b : s.bonds)
    if (b.order == 4) {
      auto dup = [&](uint32_t a, uint32_t other) {
        for (auto& e : nb[a]) if (e.atom == other && e.order == 1) { bool has2 = false; for (auto& f : nb[a]) has2 |= f.order == 2; if (!has2) e.order = 2; }
      };
      dup(b.i, b.j), dup(b.j, b.i);
    }
  return nb;
}
}  // namespace

std::vector<std::string> ez_labels(const System& s, std::vector<std::pair<uint32_t, uint32_t>>* bonds) {
  const size_t n = s.atoms.size();
  std::vector<std::string> out(n);
  const auto nb = cip_neighbours(s);
  // bonds in a ring: both ends still connected without the bond
  auto in_ring = [&](uint32_t a, uint32_t b) {
    std::vector<char> seen(n, 0);
    std::vector<uint32_t> st{a};
    seen[a] = 1;
    while (!st.empty()) {
      const uint32_t x = st.back();
      st.pop_back();
      for (const auto& e : nb[x]) {
        if ((x == a && e.atom == b) || (x == b && e.atom == a)) continue;
        if (e.atom == b) return true;
        if (!seen[e.atom]) seen[e.atom] = 1, st.push_back(e.atom);
      }
    }
    return false;
  };
  for (const auto& bd : s.bonds) {
    if (bd.order != 2) continue;
    const uint32_t a = bd.i, b = bd.j;
    if (nb[a].size() != 3 || nb[b].size() != 3) continue;   // both ends sp² with two substituents
    if (in_ring(a, b)) continue;
    auto top = [&](uint32_t end, uint32_t other) -> int {
      std::vector<int> r;
      for (const auto& e : nb[end]) if (e.atom != other) r.push_back(int(e.atom));
      if (r.size() != 2) return -1;
      const int k = compare_branches(s, nb, int(end), r[0], r[1]);
      return k == 0 ? -1 : k > 0 ? r[0] : r[1];
    };
    const int ta = top(a, b), tb = top(b, a);
    if (ta < 0 || tb < 0) continue;
    auto rel = [&](uint32_t from, int to) { Vec3 d = s.atoms[size_t(to)].pos - s.atoms[from].pos; return s.cell.valid() ? s.cell.minimum_image(d) : d; };
    const Vec3 ab = rel(a, int(b)), u = rel(a, ta), w = rel(b, tb);
    const Vec3 e = ab * (1 / norm(ab));
    const Vec3 up = u - e * dot(u, e), wp = w - e * dot(w, e);   // each group's offset across the double bond
    const std::string lab = dot(up, wp) > 0 ? "Z" : "E";          // the higher-ranked groups on one side: Z (cis)
    out[a] = out[b] = lab;
    if (bonds) bonds->push_back({a, b});
  }
  return out;
}

std::vector<std::string> stereo_labels(const System& s, const std::vector<char>& only) {
  const size_t n = s.atoms.size();
  std::vector<std::string> out(n);
  std::vector<std::vector<Neighbour>> nb(n);
  for (const auto& b : s.bonds) {
    const int order = b.order == 2 ? 2 : b.order == 3 ? 3 : 1;
    nb[b.i].push_back({b.j, order}), nb[b.j].push_back({b.i, order});
  }
  // aromatic atoms: one duplicate (a Kekulé double bond) towards the first aromatic neighbour
  for (const auto& b : s.bonds)
    if (b.order == 4) {
      auto dup = [&](uint32_t a, uint32_t other) {
        for (auto& e : nb[a]) if (e.atom == other && e.order == 1) { bool has2 = false; for (auto& f : nb[a]) has2 |= f.order == 2; if (!has2) e.order = 2; }
      };
      dup(b.i, b.j), dup(b.j, b.i);
    }
  for (size_t c = 0; c < n; ++c) {
    if (!only.empty() && (c >= only.size() || !only[c])) continue;
    if (nb[c].size() != 4 || s.atoms[c].element == 1) continue;
    bool multiple = false;
    for (const auto& e : nb[c]) multiple |= e.order != 1;
    if (multiple) continue;
    std::vector<int> r;
    for (const auto& e : nb[c]) r.push_back(int(e.atom));
    bool tie = false;
    std::sort(r.begin(), r.end(), [&](int a, int b) {
      const int k = compare_branches(s, nb, int(c), a, b);
      if (k == 0 && a != b) tie = true;
      return k > 0;
    });
    for (size_t i = 0; i + 1 < r.size() && !tie; ++i) tie = compare_branches(s, nb, int(c), r[i], r[i + 1]) == 0;
    if (tie) continue;
    Vec3 C = s.atoms[c].pos;
    auto rel = [&](int a) { Vec3 d = s.atoms[size_t(a)].pos - C; return s.cell.valid() ? s.cell.minimum_image(d) : d; };
    const double v = dot(rel(r[0]), cross(rel(r[1]), rel(r[2])));
    out[c] = v < 0 ? "R" : "S";
  }
  return out;
}

double Mesh::area() const {
  double a = 0;
  for (const auto& t : triangles) a += 0.5 * norm(cross(vertices[t[1]] - vertices[t[0]], vertices[t[2]] - vertices[t[0]]));
  return a;
}

Mesh surface_mesh(const System& s, const SurfaceOptions& o) {
  Mesh m;
  const size_t n = s.atoms.size();
  std::vector<size_t> use;
  for (size_t i = 0; i < n; ++i) if (o.atoms.empty() || (i < o.atoms.size() && o.atoms[i])) use.push_back(i);
  if (use.empty()) return m;
  const double probe = o.kind == SurfaceKind::VanDerWaals ? 0.0 : std::max(0.0, o.probe);
  auto R = [&](size_t i) { return element(s.atoms[i].element).vdw + probe; };
  Vec3 lo{1e300, 1e300, 1e300}, hi{-1e300, -1e300, -1e300};
  double rmax = 0;
  for (size_t i : use) {
    rmax = std::max(rmax, R(i));
    for (int k = 0; k < 3; ++k) lo[k] = std::min(lo[k], s.atoms[i].pos[k]), hi[k] = std::max(hi[k], s.atoms[i].pos[k]);
  }
  double h = std::max(0.2, o.spacing);
  const double margin = rmax + 2 * h;
  for (int k = 0; k < 3; ++k) lo[k] -= margin, hi[k] += margin;
  auto dims = [&](double hh, int* d) { for (int k = 0; k < 3; ++k) d[k] = int(std::ceil((hi[k] - lo[k]) / hh)) + 1; };
  int nd[3];
  dims(h, nd);
  while (double(nd[0]) * nd[1] * nd[2] > 3.0e7) { h *= 1.25; dims(h, nd); }
  const size_t N = size_t(nd[0]) * nd[1] * nd[2];
  auto idx = [&](int x, int y, int z) { return (size_t(z) * nd[1] + y) * nd[0] + x; };
  std::vector<float> f(N, float(margin));
  for (size_t i : use) {
    const Vec3 p = s.atoms[i].pos;
    const double r = R(i);
    const int x0 = std::max(0, int((p[0] - r - lo[0]) / h) - 1), x1 = std::min(nd[0] - 1, int((p[0] + r - lo[0]) / h) + 2);
    const int y0 = std::max(0, int((p[1] - r - lo[1]) / h) - 1), y1 = std::min(nd[1] - 1, int((p[1] + r - lo[1]) / h) + 2);
    const int z0 = std::max(0, int((p[2] - r - lo[2]) / h) - 1), z1 = std::min(nd[2] - 1, int((p[2] + r - lo[2]) / h) + 2);
    for (int z = z0; z <= z1; ++z)
      for (int y = y0; y <= y1; ++y)
        for (int x = x0; x <= x1; ++x) {
          const Vec3 g{lo[0] + x * h, lo[1] + y * h, lo[2] + z * h};
          const float v = float(norm(g - p) - r);
          float& t = f[idx(x, y, z)];
          if (v < t) t = v;
        }
  }
  if (o.kind == SurfaceKind::Excluded && probe > 0) {
    // the excluded surface: points farther than the probe from every accessible probe centre (f > 0) are inside
    std::vector<double> d(N);
    for (size_t k = 0; k < N; ++k) d[k] = f[k] > 0 ? 0.0 : 1e20;
    const int nmax = std::max({nd[0], nd[1], nd[2]});
    std::vector<double> line(static_cast<size_t>(nmax)), res(static_cast<size_t>(nmax)), zz(static_cast<size_t>(nmax) + 1);
    std::vector<int> v(static_cast<size_t>(nmax));
    for (int axis = 0; axis < 3; ++axis) {
      const int a1 = (axis + 1) % 3, a2 = (axis + 2) % 3;
      for (int i = 0; i < nd[a1]; ++i)
        for (int j = 0; j < nd[a2]; ++j) {
          line.resize(size_t(nd[axis])), res.resize(size_t(nd[axis]));
          for (int q = 0; q < nd[axis]; ++q) {
            int c[3]; c[axis] = q, c[a1] = i, c[a2] = j;
            line[size_t(q)] = d[idx(c[0], c[1], c[2])];
          }
          edt1(line, res, v, zz);
          for (int q = 0; q < nd[axis]; ++q) {
            int c[3]; c[axis] = q, c[a1] = i, c[a2] = j;
            d[idx(c[0], c[1], c[2])] = res[size_t(q)];
          }
        }
    }
    for (size_t k = 0; k < N; ++k) f[k] = float(probe - std::sqrt(d[k]) * h);   // < 0 inside the excluded volume
  }
  // marching tetrahedra: six per cube along the 0–7 diagonal; vertices welded on grid edges
  static const int tets[6][4] = {{0, 1, 3, 7}, {0, 3, 2, 7}, {0, 2, 6, 7}, {0, 6, 4, 7}, {0, 4, 5, 7}, {0, 5, 1, 7}};
  std::unordered_map<uint64_t, uint32_t> edge_vertex;
  auto vertex_on = [&](size_t a, size_t b, const Vec3& pa, const Vec3& pb, float fa, float fb) {
    const uint64_t key = uint64_t(std::min(a, b)) * uint64_t(N) + uint64_t(std::max(a, b));
    auto it = edge_vertex.find(key);
    if (it != edge_vertex.end()) return it->second;
    const double t = std::clamp(double(fa) / double(fa - fb), 0.0, 1.0);
    m.vertices.push_back(pa + (pb - pa) * t);
    m.normals.push_back({0, 0, 0});
    const uint32_t id = uint32_t(m.vertices.size() - 1);
    edge_vertex.emplace(key, id);
    return id;
  };
  auto tri = [&](uint32_t a, uint32_t b, uint32_t c, const Vec3& outward) {
    Vec3 nrm = cross(m.vertices[b] - m.vertices[a], m.vertices[c] - m.vertices[a]);
    if (dot(nrm, outward) < 0) { std::swap(b, c); nrm = nrm * -1.0; }
    if (norm(nrm) < 1e-14) return;
    m.triangles.push_back({a, b, c});
    for (uint32_t q : {a, b, c}) m.normals[q] = m.normals[q] + nrm;
  };
  for (int z = 0; z + 1 < nd[2]; ++z)
    for (int y = 0; y + 1 < nd[1]; ++y)
      for (int x = 0; x + 1 < nd[0]; ++x) {
        size_t ci[8];
        float cf[8];
        Vec3 cp[8];
        bool any_in = false, any_out = false;
        for (int c = 0; c < 8; ++c) {
          const int dx = c & 1, dy = (c >> 1) & 1, dz = (c >> 2) & 1;
          ci[c] = idx(x + dx, y + dy, z + dz);
          cf[c] = f[ci[c]];
          cp[c] = {lo[0] + (x + dx) * h, lo[1] + (y + dy) * h, lo[2] + (z + dz) * h};
          (cf[c] < 0 ? any_in : any_out) = true;
        }
        if (!any_in || !any_out) continue;
        for (const auto& t : tets) {
          std::vector<int> in, out;
          for (int q : t) (cf[q] < 0 ? in : out).push_back(q);
          if (in.empty() || out.empty()) continue;
          Vec3 cin{0, 0, 0}, cout{0, 0, 0};
          for (int q : in) cin = cin + cp[q];
          for (int q : out) cout = cout + cp[q];
          const Vec3 outward = cout * (1.0 / double(out.size())) - cin * (1.0 / double(in.size()));
          auto V = [&](int a, int b) { return vertex_on(ci[a], ci[b], cp[a], cp[b], cf[a], cf[b]); };
          if (in.size() == 1) tri(V(in[0], out[0]), V(in[0], out[1]), V(in[0], out[2]), outward);
          else if (in.size() == 3) tri(V(out[0], in[0]), V(out[0], in[1]), V(out[0], in[2]), outward);
          else {
            const uint32_t a = V(in[0], out[0]), b = V(in[0], out[1]), c = V(in[1], out[1]), d = V(in[1], out[0]);
            tri(a, b, c, outward);
            tri(a, c, d, outward);
          }
        }
      }
  for (auto& v : m.normals) v = unitv(v);
  // the nearest wrapped atom of every vertex (a bin grid over the atoms)
  m.nearest.assign(m.vertices.size(), -1);
  const double bin = std::max(3.0, rmax + 2 * h);
  int nb[3];
  for (int k = 0; k < 3; ++k) nb[k] = std::max(1, int((hi[k] - lo[k]) / bin) + 1);
  std::vector<std::vector<uint32_t>> bins(size_t(nb[0]) * nb[1] * nb[2]);
  auto bi = [&](const Vec3& p, int k) { return std::clamp(int((p[k] - lo[k]) / bin), 0, nb[k] - 1); };
  for (size_t i : use) bins[(size_t(bi(s.atoms[i].pos, 2)) * nb[1] + bi(s.atoms[i].pos, 1)) * nb[0] + bi(s.atoms[i].pos, 0)].push_back(uint32_t(i));
  for (size_t v = 0; v < m.vertices.size(); ++v) {
    const Vec3 p = m.vertices[v];
    double best = 1e300;
    const int b0 = bi(p, 0), b1 = bi(p, 1), b2 = bi(p, 2);
    for (int zz = std::max(0, b2 - 1); zz <= std::min(nb[2] - 1, b2 + 1); ++zz)
      for (int yy = std::max(0, b1 - 1); yy <= std::min(nb[1] - 1, b1 + 1); ++yy)
        for (int xx = std::max(0, b0 - 1); xx <= std::min(nb[0] - 1, b0 + 1); ++xx)
          for (uint32_t i : bins[(size_t(zz) * nb[1] + yy) * nb[0] + xx]) {
            const double d = norm(s.atoms[i].pos - p) - R(i);
            if (d < best) best = d, m.nearest[v] = int32_t(i);
          }
  }
  return m;
}

Mesh isosurface(const Cell& cell, const int n[3], const std::vector<double>& values, double level, bool periodic) {
  Mesh m;
  if (n[0] < 2 || n[1] < 2 || n[2] < 2 || values.size() != size_t(n[0]) * size_t(n[1]) * size_t(n[2])) return m;
  const size_t N = values.size();
  auto idx = [&](int i, int j, int k) { return (size_t(i) * size_t(n[1]) + size_t(j)) * size_t(n[2]) + size_t(k); };
  auto point = [&](int i, int j, int k) {
    return cell.origin + cell.a * ((i + 0.5) / n[0]) + cell.b * ((j + 0.5) / n[1]) + cell.c * ((k + 0.5) / n[2]);
  };
  // marching tetrahedra (six per cube along the 0–7 diagonal) on f = level − value: < 0 above the level ("inside")
  static const int tets[6][4] = {{0, 1, 3, 7}, {0, 3, 2, 7}, {0, 2, 6, 7}, {0, 6, 4, 7}, {0, 4, 5, 7}, {0, 5, 1, 7}};
  // vertices welded on grid edges; an edge is keyed by its two corners' unwrapped grid points (a wrapped corner keeps
  // its own position past the face, so the surface there is continuous with the image across it)
  std::unordered_map<uint64_t, uint32_t> edge_vertex;
  const int ext = periodic ? 1 : 0;
  auto key_of = [&](int i, int j, int k) { return (uint64_t(i) * uint64_t(n[1] + 1) + uint64_t(j)) * uint64_t(n[2] + 1) + uint64_t(k); };
  const uint64_t NK = uint64_t(n[0] + 1) * uint64_t(n[1] + 1) * uint64_t(n[2] + 1);
  auto tri = [&](uint32_t a, uint32_t b, uint32_t c, const Vec3& outward) {
    Vec3 nrm = cross(m.vertices[b] - m.vertices[a], m.vertices[c] - m.vertices[a]);
    if (dot(nrm, outward) < 0) { std::swap(b, c); nrm = nrm * -1.0; }
    if (norm(nrm) < 1e-14) return;
    m.triangles.push_back({a, b, c});
    for (uint32_t q : {a, b, c}) m.normals[q] = m.normals[q] + nrm;
  };
  (void)N;
  for (int x = 0; x + 1 < n[0] + ext; ++x)
    for (int y = 0; y + 1 < n[1] + ext; ++y)
      for (int z = 0; z + 1 < n[2] + ext; ++z) {
        uint64_t ck[8];
        double cf[8];
        Vec3 cp[8];
        bool any_in = false, any_out = false;
        for (int c = 0; c < 8; ++c) {
          const int i = x + (c & 1), j = y + ((c >> 1) & 1), k = z + ((c >> 2) & 1);
          ck[c] = key_of(i, j, k);
          cf[c] = level - values[idx(i % n[0], j % n[1], k % n[2])];
          if (cf[c] == 0) cf[c] = -1e-12;
          cp[c] = point(i, j, k);
          (cf[c] < 0 ? any_in : any_out) = true;
        }
        if (!any_in || !any_out) continue;
        auto V = [&](int a, int b) {
          const uint64_t key = std::min(ck[a], ck[b]) * NK + std::max(ck[a], ck[b]);
          auto it = edge_vertex.find(key);
          if (it != edge_vertex.end()) return it->second;
          const double t = std::clamp(cf[a] / (cf[a] - cf[b]), 0.0, 1.0);
          m.vertices.push_back(cp[a] + (cp[b] - cp[a]) * t);
          m.normals.push_back({0, 0, 0});
          const uint32_t id = uint32_t(m.vertices.size() - 1);
          edge_vertex.emplace(key, id);
          return id;
        };
        for (const auto& t : tets) {
          int in[4], out[4], ni = 0, no = 0;
          for (int q : t) (cf[q] < 0 ? in[ni++] : out[no++]) = q;
          if (ni == 0 || no == 0) continue;
          Vec3 cin{0, 0, 0}, cout{0, 0, 0};
          for (int q = 0; q < ni; ++q) cin = cin + cp[in[q]];
          for (int q = 0; q < no; ++q) cout = cout + cp[out[q]];
          const Vec3 outward = cout * (1.0 / no) - cin * (1.0 / ni);
          if (ni == 1) tri(V(in[0], out[0]), V(in[0], out[1]), V(in[0], out[2]), outward);
          else if (ni == 3) tri(V(out[0], in[0]), V(out[0], in[1]), V(out[0], in[2]), outward);
          else {
            const uint32_t a = V(in[0], out[0]), b = V(in[0], out[1]), c = V(in[1], out[1]), d = V(in[1], out[0]);
            tri(a, b, c, outward);
            tri(a, c, d, outward);
          }
        }
      }
  for (auto& v : m.normals) v = unitv(v);
  return m;
}

std::vector<double> surface_potential(const System& s, const Mesh& m, const std::vector<char>& atoms) {
  std::vector<double> phi(m.vertices.size(), 0.0);
  std::vector<size_t> charged;
  for (size_t i = 0; i < s.atoms.size(); ++i)
    if (std::fabs(s.atoms[i].charge) > 1e-9 && (atoms.empty() || (i < atoms.size() && atoms[i]))) charged.push_back(i);
  constexpr double kCoulomb = 332.0637;   // kcal Å / (mol e²)
  if (double(charged.size()) * double(m.vertices.size()) < 5e7) {
    for (size_t v = 0; v < m.vertices.size(); ++v) {
      double sum = 0;
      for (size_t i : charged) sum += s.atoms[i].charge / std::max(0.5, norm(m.vertices[v] - s.atoms[i].pos));
      phi[v] = kCoulomb * sum;
    }
    return phi;
  }
  // large systems: charges within 15 Å (bins of 15 Å)
  constexpr double kCut = 15.0;
  Vec3 lo{1e300, 1e300, 1e300};
  for (size_t i : charged) for (int k = 0; k < 3; ++k) lo[k] = std::min(lo[k], s.atoms[i].pos[k]);
  for (const auto& v : m.vertices) for (int k = 0; k < 3; ++k) lo[k] = std::min(lo[k], v[k]);
  std::unordered_map<int64_t, std::vector<size_t>> bins;
  auto key = [&](const Vec3& p, int dx, int dy, int dz) {
    const int64_t x = int64_t((p[0] - lo[0]) / kCut) + dx, y = int64_t((p[1] - lo[1]) / kCut) + dy, z = int64_t((p[2] - lo[2]) / kCut) + dz;
    return (x * 100003 + y) * 100003 + z;
  };
  for (size_t i : charged) bins[key(s.atoms[i].pos, 0, 0, 0)].push_back(i);
  for (size_t v = 0; v < m.vertices.size(); ++v) {
    double sum = 0;
    for (int dx = -1; dx <= 1; ++dx)
      for (int dy = -1; dy <= 1; ++dy)
        for (int dz = -1; dz <= 1; ++dz) {
          auto it = bins.find(key(m.vertices[v], dx, dy, dz));
          if (it == bins.end()) continue;
          for (size_t i : it->second) {
            const double r = norm(m.vertices[v] - s.atoms[i].pos);
            if (r < kCut) sum += s.atoms[i].charge / std::max(0.5, r);
          }
        }
    phi[v] = kCoulomb * sum;
  }
  return phi;
}

Mesh polyhedra(const System& s, const std::vector<char>& centres) {
  Mesh m;
  const size_t n = s.atoms.size();
  std::vector<std::vector<uint32_t>> nb(n);
  for (const auto& b : s.bonds) nb[b.i].push_back(b.j), nb[b.j].push_back(b.i);
  for (size_t c = 0; c < n; ++c) {
    if (!centres.empty() && (c >= centres.size() || !centres[c])) continue;
    std::vector<Vec3> p;
    for (uint32_t j : nb[c]) {
      if (s.atoms[j].element == 1) continue;
      Vec3 d = s.atoms[j].pos - s.atoms[c].pos;
      if (s.cell.valid()) d = s.cell.minimum_image(d);
      p.push_back(s.atoms[c].pos + d);
    }
    if (p.size() < 4 || p.size() > 8) continue;
    const Vec3 C = s.atoms[c].pos;
    const unsigned rgb = element_colour(s.atoms[c].element);
    for (size_t a = 0; a < p.size(); ++a)
      for (size_t b = a + 1; b < p.size(); ++b)
        for (size_t d = b + 1; d < p.size(); ++d) {
          Vec3 nrm = cross(p[b] - p[a], p[d] - p[a]);
          const double len = norm(nrm);
          if (len < 1e-8) continue;
          nrm = nrm * (1.0 / len);
          int above = 0, below = 0;
          for (size_t q = 0; q < p.size(); ++q) {
            if (q == a || q == b || q == d) continue;
            const double t = dot(p[q] - p[a], nrm);
            if (t > 1e-3) ++above; else if (t < -1e-3) ++below;
          }
          if (above && below) continue;
          // a hull face: outward, away from the centre
          const bool flip = dot(C - p[a], nrm) > 0;
          const uint32_t base = uint32_t(m.vertices.size());
          for (const Vec3& q : {p[a], p[b], p[d]}) { m.vertices.push_back(q); m.normals.push_back(flip ? nrm * -1.0 : nrm); m.colours.push_back(rgb); m.nearest.push_back(int32_t(c)); }
          m.triangles.push_back(flip ? std::array<uint32_t, 3>{base, base + 2, base + 1} : std::array<uint32_t, 3>{base, base + 1, base + 2});
        }
  }
  return m;
}

std::vector<std::vector<Vec3>> ribbon_paths(const System& s, const std::vector<char>& atoms) {
  std::vector<std::vector<Vec3>> out;
  for (const auto& chain : backbones(s, 4)) {
    // runs of the chain's flagged atoms, unwrapped along the chain
    std::vector<Vec3> pts;
    auto flush = [&] {
      if (pts.size() >= 2) {
        std::vector<Vec3> smooth;
        for (size_t i = 0; i + 1 < pts.size(); ++i) {
          const Vec3 &p0 = pts[i ? i - 1 : 0], &p1 = pts[i], &p2 = pts[i + 1], &p3 = pts[std::min(i + 2, pts.size() - 1)];
          const int steps = std::max(2, int(norm(p2 - p1) / 0.5));
          for (int k = 0; k < steps; ++k) {
            const double t = double(k) / steps, t2 = t * t, t3 = t2 * t;
            Vec3 q;
            for (int d = 0; d < 3; ++d)
              q[d] = 0.5 * (2 * p1[d] + (-p0[d] + p2[d]) * t + (2 * p0[d] - 5 * p1[d] + 4 * p2[d] - p3[d]) * t2 + (-p0[d] + 3 * p1[d] - 3 * p2[d] + p3[d]) * t3);
            smooth.push_back(q);
          }
        }
        smooth.push_back(pts.back());
        out.push_back(std::move(smooth));
      }
      pts.clear();
    };
    for (uint32_t a : chain) {
      if (!atoms.empty() && (a >= atoms.size() || !atoms[a])) { flush(); continue; }
      Vec3 p = s.atoms[a].pos;
      if (!pts.empty() && s.cell.valid()) p = pts.back() + s.cell.minimum_image(p - pts.back());
      pts.push_back(p);
    }
    flush();
  }
  return out;
}

}  // namespace caps
