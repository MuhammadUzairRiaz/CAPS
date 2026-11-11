// CAPS crystals and surfaces (see caps/crystal.hpp).
#include "caps/crystal.hpp"
#include "caps/spacegroup.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <functional>
#include <fstream>
#include <map>
#include <numeric>
#include <random>
#include <set>
#include <sstream>

#include "caps/analysis.hpp"
#include "caps/elements.hpp"

namespace caps {

namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr double kDeg = kPi / 180.0;

// ---------------------------------------------------------------- CIF tokens

struct CifData {
  std::string name;
  std::map<std::string, std::string> items;                                  // tag (lower case) → value
  std::vector<std::pair<std::vector<std::string>, std::vector<std::vector<std::string>>>> loops;   // tags, rows
};

std::vector<std::string> cif_tokens(const std::string& text) {
  std::vector<std::string> out;
  std::istringstream in(text);
  std::string line;
  while (std::getline(in, line)) {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (!line.empty() && line[0] == ';') {   // text field: up to the next line starting with ';'
      std::string field = line.substr(1);
      while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (!line.empty() && line[0] == ';') break;
        field += "\n" + line;
      }
      out.push_back(field);
      continue;
    }
    size_t i = 0;
    while (i < line.size()) {
      while (i < line.size() && std::isspace(static_cast<unsigned char>(line[i]))) ++i;
      if (i >= line.size()) break;
      if (line[i] == '#') break;
      if (line[i] == '\'' || line[i] == '"') {
        const char q = line[i];
        size_t j = i + 1;
        while (j < line.size() && !(line[j] == q && (j + 1 == line.size() || std::isspace(static_cast<unsigned char>(line[j + 1]))))) ++j;
        out.push_back(line.substr(i + 1, j - i - 1));
        i = j + 1;
        continue;
      }
      size_t j = i;
      while (j < line.size() && !std::isspace(static_cast<unsigned char>(line[j]))) ++j;
      out.push_back(line.substr(i, j - i));
      i = j;
    }
  }
  return out;
}

std::string lower(std::string s) {
  for (auto& c : s) c = char(std::tolower(static_cast<unsigned char>(c)));
  return s;
}

CifData cif_parse(const std::string& text) {
  const auto t = cif_tokens(text);
  CifData d;
  bool in_block = false;
  for (size_t i = 0; i < t.size();) {
    const std::string lt = lower(t[i]);
    if (lt.rfind("data_", 0) == 0) {
      if (in_block) break;   // the first data block only
      in_block = true;
      d.name = t[i].substr(5);
      ++i;
    } else if (lt == "loop_") {
      ++i;
      std::vector<std::string> tags;
      while (i < t.size() && !t[i].empty() && t[i][0] == '_') tags.push_back(lower(t[i++]));
      std::vector<std::string> vals;
      while (i < t.size() && !(t[i].size() && t[i][0] == '_') && lower(t[i]) != "loop_" && lower(t[i]).rfind("data_", 0) != 0) vals.push_back(t[i++]);
      std::vector<std::vector<std::string>> rows;
      if (!tags.empty())
        for (size_t r = 0; r + tags.size() <= vals.size(); r += tags.size()) rows.emplace_back(vals.begin() + long(r), vals.begin() + long(r + tags.size()));
      d.loops.push_back({tags, rows});
    } else if (!t[i].empty() && t[i][0] == '_') {
      if (i + 1 < t.size()) d.items[lt] = t[i + 1];
      i += 2;
    } else {
      ++i;
    }
  }
  return d;
}

double cif_number(const std::string& v, double fallback = 0) {
  if (v.empty() || v == "?" || v == ".") return fallback;
  std::string s = v.substr(0, v.find('('));
  try {
    return std::stod(s);
  } catch (...) {
    return fallback;
  }
}

// ---------------------------------------------------------------- symmetry operations

// Element from a CIF type symbol or label: "O2-" → O, "Si4+" → Si, "Fe1" → Fe, "OW1" → O.
int element_of_site(const std::string& sym) {
  std::string letters;
  for (char c : sym) {
    if (std::isalpha(static_cast<unsigned char>(c))) letters += c;
    else break;
  }
  if (letters.empty()) return 0;
  if (letters.size() >= 2) {
    std::string two = {char(std::toupper(static_cast<unsigned char>(letters[0]))), char(std::tolower(static_cast<unsigned char>(letters[1])))};
    if (int z = element_from_symbol(two)) return z;
  }
  return element_from_symbol(std::string(1, char(std::toupper(static_cast<unsigned char>(letters[0])))));
}

Vec3 wrap01(Vec3 f) {
  for (auto& x : f) {
    x -= std::floor(x);
    if (x >= 1 - 1e-9) x = 0;
  }
  return f;
}

// 3×3 inverse of the matrix with columns u, v, w (returns its rows).
std::array<Vec3, 3> inverse_cols(const Vec3& u, const Vec3& v, const Vec3& w) {
  const double det = dot(u, cross(v, w));
  if (std::fabs(det) < 1e-12) throw CrystalError("degenerate cell");
  return {cross(v, w) * (1 / det), cross(w, u) * (1 / det), cross(u, v) * (1 / det)};
}

Vec3 solve(const std::array<Vec3, 3>& inv, const Vec3& r) { return {dot(inv[0], r), dot(inv[1], r), dot(inv[2], r)}; }

std::string hill(const std::map<int, int>& count) {
  std::string f;
  auto part = [&](int z, int n) { f += std::string(element(z).symbol) + (n > 1 ? std::to_string(n) : ""); };
  const bool carbon = count.count(6) > 0;
  if (carbon) {
    part(6, count.at(6));
    if (count.count(1)) part(1, count.at(1));
  }
  std::vector<std::pair<std::string, int>> rest;
  for (const auto& [z, n] : count)
    if (!(carbon && (z == 6 || z == 1))) rest.push_back({element(z).symbol, z});
  std::sort(rest.begin(), rest.end());
  for (const auto& [sym, z] : rest) part(z, count.at(z));
  return f;
}

int gcd3(int a, int b, int c) { return std::gcd(std::gcd(std::abs(a), std::abs(b)), std::abs(c)); }

// The (hkl) plane of a bulk cell: the surface mesh s1, s2 (Lagrange-reduced, right-handed about n), the stacking
// vector s3 (one plane spacing along n), and every atom's coordinates in that basis, wrapped into [0, 1).
struct PlaneSetup {
  int h = 0, k = 0, l = 1;
  Vec3 s1, s2, s3, n;
  double d = 0, area = 0;
  std::vector<Vec3> f;   // per atom (α, β, γ)
};

PlaneSetup plane_setup(const System& bulk, int h, int k, int l) {
  if (!bulk.cell.valid()) throw CrystalError("the structure has no periodic cell to cleave");
  const int g = gcd3(h, k, l);
  if (g == 0) throw CrystalError("(000) is not a plane");
  PlaneSetup P;
  P.h = h / g, P.k = k / g, P.l = l / g;
  const Vec3 a = bulk.cell.a, b = bulk.cell.b, c = bulk.cell.c;
  const double V = std::fabs(dot(a, cross(b, c)));
  const double sgn = dot(a, cross(b, c)) > 0 ? 1 : -1;
  const Vec3 as = cross(b, c) * (sgn / V), bs = cross(c, a) * (sgn / V), cs = cross(a, b) * (sgn / V);
  const Vec3 G = as * P.h + bs * P.k + cs * P.l;
  P.d = 1 / norm(G);
  P.n = G * P.d;
  const double area_min = V / P.d;
  const int R = 10;
  std::vector<std::pair<double, Vec3>> in_plane;
  for (int u = -R; u <= R; ++u)
    for (int v = -R; v <= R; ++v)
      for (int w = -R; w <= R; ++w) {
        if ((u || v || w) && P.h * u + P.k * v + P.l * w == 0) {
          const Vec3 t = a * u + b * v + c * w;
          in_plane.push_back({norm(t), t});
        }
      }
  std::sort(in_plane.begin(), in_plane.end(), [](const auto& x, const auto& y) { return x.first < y.first; });
  const size_t m = std::min<size_t>(in_plane.size(), 300);
  double best = 1e300;
  for (size_t i = 0; i < m; ++i)
    for (size_t j = i + 1; j < m; ++j) {
      const double ar = norm(cross(in_plane[i].second, in_plane[j].second));
      if (std::fabs(ar - area_min) > 1e-4 * area_min) continue;
      const double score = in_plane[i].first + in_plane[j].first;
      if (score < best) { best = score; P.s1 = in_plane[i].second; P.s2 = in_plane[j].second; }
    }
  if (best > 1e299) throw CrystalError("no surface mesh found for this (hkl); try lower indices");
  // Lagrange reduction
  for (int it = 0; it < 100; ++it) {
    if (norm(P.s2) < norm(P.s1)) std::swap(P.s1, P.s2);
    const double mu = std::round(dot(P.s1, P.s2) / dot(P.s1, P.s1));
    if (mu == 0) break;
    P.s2 = P.s2 - P.s1 * mu;
  }
  // stacking vector: h u + k v + l w = 1, the smallest in-plane component
  double bi = 1e300;
  for (int u = -R; u <= R; ++u)
    for (int v = -R; v <= R; ++v)
      for (int w = -R; w <= R; ++w)
        if (P.h * u + P.k * v + P.l * w == 1) {
          const Vec3 t = a * u + b * v + c * w;
          const double ip = norm(t - P.n * dot(t, P.n));
          if (ip < bi - 1e-9) { bi = ip; P.s3 = t; }
        }
  if (bi > 1e299) throw CrystalError("no stacking vector for this (hkl)");
  if (dot(cross(P.s1, P.s2), P.n) < 0) std::swap(P.s1, P.s2);
  P.area = norm(cross(P.s1, P.s2));
  const auto inv = inverse_cols(P.s1, P.s2, P.s3);
  for (const auto& at : bulk.atoms) P.f.push_back(wrap01(solve(inv, at.pos - bulk.cell.origin)));
  return P;
}

// Bonds of a periodic structure with every image (a small cell bonds an atom to several images of another), from
// covalent radii + 0.45 Å as perceive_bonds; H–H pairs are not bonds.
struct ImageBond { uint32_t i, j; Vec3 d; };

bool metal(int z) {
  if (z == 3 || z == 4 || z == 11 || z == 12 || z == 13 || z == 19 || z == 20 || z == 37 || z == 38 || z == 55 || z == 56 || z == 87 || z == 88) return true;
  if ((z >= 21 && z <= 31) || (z >= 39 && z <= 50) || (z >= 57 && z <= 83) || z >= 89) return true;
  return false;
}

std::vector<ImageBond> image_bonds(const System& s) {
  std::vector<ImageBond> out;
  // in a crystal with non-metals, metal–metal contacts are not bonds (Ti–Ti in rutile); in a metal, only nearest
  // neighbours are (not the second shell of bcc iron)
  bool all_metal = !s.atoms.empty();
  for (const auto& at : s.atoms) all_metal = all_metal && metal(at.element);
  const Vec3 a = s.cell.a, b = s.cell.b, c = s.cell.c;
  const double V = std::fabs(dot(a, cross(b, c)));
  double rmax = 0;
  for (const auto& at : s.atoms) rmax = std::max(rmax, element(at.element).covalent);
  rmax = 2 * rmax + 0.8;
  const double w[3] = {V / norm(cross(b, c)), V / norm(cross(c, a)), V / norm(cross(a, b))};
  int n[3];
  for (int k = 0; k < 3; ++k) n[k] = s.cell.periodic[size_t(k)] ? int(std::ceil(rmax / w[k])) : 0;
  const size_t N = s.atoms.size();
  for (uint32_t i = 0; i < N; ++i)
    for (uint32_t j = i; j < N; ++j) {
      const int zi = s.atoms[i].element, zj = s.atoms[j].element;
      if (zi == 1 && zj == 1) continue;
      const double lim = element(zi).covalent + element(zj).covalent + (metal(zi) && metal(zj) ? 0.8 : 0.45);
      const Vec3 d0 = s.atoms[j].pos - s.atoms[i].pos;
      for (int p = -n[0]; p <= n[0]; ++p)
        for (int q = -n[1]; q <= n[1]; ++q)
          for (int r = -n[2]; r <= n[2]; ++r) {
            if (i == j && (p < 0 || (p == 0 && (q < 0 || (q == 0 && r <= 0))))) continue;   // each self-image pair once
            const Vec3 d = d0 + a * p + b * q + c * r;
            const double L = norm(d);
            if (L >= 0.4 && L <= lim) out.push_back({i, j, d});
          }
    }
  if (all_metal) {
    double dmin = 1e300;
    for (const auto& b : out) dmin = std::min(dmin, norm(b.d));
    out.erase(std::remove_if(out.begin(), out.end(), [&](const ImageBond& b) { return norm(b.d) > 1.15 * dmin; }), out.end());
  } else {
    out.erase(std::remove_if(out.begin(), out.end(), [&](const ImageBond& b) { return metal(s.atoms[b.i].element) && metal(s.atoms[b.j].element); }),
              out.end());
  }
  return out;
}


struct Plane { double z; std::vector<size_t> atoms; };

std::vector<Plane> atomic_planes(const PlaneSetup& P) {
  std::vector<std::pair<double, size_t>> zs;
  for (size_t i = 0; i < P.f.size(); ++i) zs.push_back({P.f[i][2] * P.d, i});
  std::sort(zs.begin(), zs.end());
  const double tol = 0.2;
  std::vector<Plane> planes;
  for (const auto& [z, i] : zs) {
    if (planes.empty() || z - planes.back().z > tol) planes.push_back({z, {i}});
    else planes.back().atoms.push_back(i);
  }
  if (planes.size() > 1 && planes.front().z + P.d - planes.back().z < tol) {   // the same plane across the layer boundary
    auto& f = planes.front();
    f.atoms.insert(f.atoms.end(), planes.back().atoms.begin(), planes.back().atoms.end());
    planes.pop_back();
  }
  return planes;
}

// Adds an atom with the element's type (created when new); returns its index.
uint32_t add_atom(System& s, int z, const Vec3& p, const std::string& name, int64_t mol) {
  int type = 0;
  for (const auto& t : s.types)
    if (t.label == element(z).symbol) type = t.type;
  if (!type) {
    TypeInfo ti;
    ti.type = int(s.types.size()) + 1;
    ti.mass = element(z).mass;
    ti.label = element(z).symbol;
    s.types.push_back(ti);
    type = ti.type;
  }
  Atom a;
  a.id = int64_t(s.atoms.size() + 1);
  a.mol = mol;
  a.element = z;
  a.type = type;
  a.name = name;
  a.pos = p;
  s.atoms.push_back(a);
  return uint32_t(s.atoms.size() - 1);
}

}  // namespace

// The structure's bonds (one per atom pair) from its image bonds.
std::vector<Bond> crystal_bonds(const System& s) {
  std::set<std::pair<uint32_t, uint32_t>> seen;
  std::vector<Bond> out;
  for (const auto& b : image_bonds(s))
    if (b.i != b.j && seen.insert({b.i, b.j}).second) out.push_back({b.i, b.j, 0});
  return out;
}

namespace {
// Makes every finite molecule of a periodic structure whole, in place, and leaves periodic networks as placed. Walking
// the image bonds (each with its bond vector d) from a root atom gives every atom of the component a position p_j =
// p_i + d. A finite molecule (one chain of a polymer melt, one unit of a molecular crystal, a guest in a framework)
// closes every ring at the same position, so the walk is consistent and p is the molecule whole. A component that
// meets one of its own atoms again one lattice vector away (a framework, a layer, an endless crystal chain), or holds
// an atom bonded to its own image, is periodic: no choice of images makes it whole, so its positions stay as read.
// Using the stored bond vectors rather than the minimum image also stays right in a cell narrower than twice a bond.
void make_finite_molecules_whole(System& s) {
  const auto ib = image_bonds(s);
  const size_t n = s.atoms.size();
  std::vector<std::vector<std::pair<uint32_t, Vec3>>> nb(n);
  std::vector<char> self_image(n, 0);
  for (const auto& b : ib) {
    if (b.i == b.j) { self_image[b.i] = 1; continue; }
    nb[b.i].push_back({b.j, b.d});
    nb[b.j].push_back({b.i, b.d * -1.0});
  }
  std::vector<Vec3> p(n);
  std::vector<char> seen(n, 0);
  std::vector<uint32_t> stack, members;
  for (uint32_t root = 0; root < n; ++root) {
    if (seen[root]) continue;
    seen[root] = 1;
    p[root] = s.atoms[root].pos;
    stack.assign(1, root);
    members.clear();
    bool periodic = false;
    while (!stack.empty()) {
      const uint32_t i = stack.back();
      stack.pop_back();
      members.push_back(i);
      periodic = periodic || self_image[i];
      for (const auto& [j, d] : nb[i]) {
        const Vec3 q = p[i] + d;
        if (!seen[j]) {
          seen[j] = 1, p[j] = q, stack.push_back(j);
        } else if (norm(q - p[j]) > 1e-3) {
          periodic = true;   // met again one lattice vector away
        }
      }
    }
    if (!periodic)
      for (uint32_t i : members) s.atoms[i].pos = p[i];
  }
}
}  // namespace

// Adds O–H to under-coordinated oxygens (S, Se), M–OH to cations bonded to oxygen in the bulk and H to under-coordinated
// covalent network atoms of `s`, on the free vertices of the atom's own bulk polyhedron fitted to its surviving bonds
// (outward_of(position) and room breaking ties). An O site two cations both lost is put back once, bridging them.
// Metals keep their bare surfaces. Returns the H and OH added.
std::pair<int, int> passivate_surface(System& s, const System& bulk, const std::function<Vec3(const Vec3&)>& outward_of) {
  int added_h = 0, added_oh = 0;
  // bulk coordination and partners per element
  std::vector<int> bcn(bulk.atoms.size(), 0);
  std::map<int, std::map<int, int>> partners;
  std::map<std::pair<int, int>, std::pair<double, int>> blen;
  for (const auto& b : image_bonds(bulk)) {
    ++bcn[b.i], ++bcn[b.j];
    const int zi = bulk.atoms[b.i].element, zj = bulk.atoms[b.j].element;
    partners[zi][zj]++, partners[zj][zi]++;
    const double r = norm(b.d);
    auto& bl = blen[{std::min(zi, zj), std::max(zi, zj)}];
    bl.first += r, bl.second++;
  }
  std::map<int, std::map<int, int>> cn_hist;
  for (size_t i = 0; i < bulk.atoms.size(); ++i) cn_hist[bulk.atoms[i].element][bcn[i]]++;
  std::map<int, int> want;
  for (const auto& [z, hgram] : cn_hist)
    want[z] = std::max_element(hgram.begin(), hgram.end(), [](const auto& x, const auto& y) { return x.second < y.second; })->first;
  auto main_partner = [&](int z) {
    const auto it = partners.find(z);
    if (it == partners.end() || it->second.empty()) return 0;
    return std::max_element(it->second.begin(), it->second.end(), [](const auto& x, const auto& y) { return x.second < y.second; })->first;
  };
  auto bond_length = [&](int a, int b) {
    const auto it = blen.find({std::min(a, b), std::max(a, b)});
    return it != blen.end() ? it->second.first / it->second.second : element(a).covalent + element(b).covalent;
  };
  // directions on a sphere
  std::vector<Vec3> dirs;
  const int nd = 800;
  for (int i = 0; i < nd; ++i) {
    const double y = 1 - 2 * (i + 0.5) / nd, r = std::sqrt(1 - y * y), phi = i * kPi * (3 - std::sqrt(5.0));
    dirs.push_back({r * std::cos(phi), r * std::sin(phi), y});
  }
  std::vector<std::vector<Vec3>> bonded(s.atoms.size());
  std::vector<std::vector<uint32_t>> bonded_to(s.atoms.size());
  for (const auto& b : image_bonds(s)) {
    bonded[b.i].push_back(b.d * (1 / norm(b.d)));
    bonded[b.j].push_back(b.d * (-1 / norm(b.d)));
    bonded_to[b.i].push_back(b.j), bonded_to[b.j].push_back(b.i);
  }
  // metals keep their bare surfaces: only oxygen (and S, Se), cations of oxides and covalent network atoms
  auto passivable = [&](int z) {
    if (z == 8 || z == 16 || z == 34) return true;
    if (main_partner(z) == 8) return true;
    return z == 5 || z == 6 || z == 7 || z == 14 || z == 15 || z == 32;
  };
  // The bulk polyhedra: the bond vectors of each element's atoms, one per distinct shape (an element on two kinds of
  // site, Fe in magnetite, has a tetrahedron and an octahedron). A surface atom's missing bonds are the vertices of its
  // own bulk polyhedron that its surviving bonds leave free, all placed together, so a cation that lost two or three
  // bonds still gets a tetrahedron or an octahedron.
  struct Poly { std::vector<Vec3> v; std::vector<int> pz; };
  std::map<int, std::vector<Poly>> polys;
  {
    std::vector<Poly> all(bulk.atoms.size());
    for (const auto& b : image_bonds(bulk)) {
      all[b.i].v.push_back(b.d), all[b.i].pz.push_back(bulk.atoms[b.j].element);
      all[b.j].v.push_back(b.d * -1.0), all[b.j].pz.push_back(bulk.atoms[b.i].element);
    }
    auto shape = [](const Poly& p) {
      std::vector<double> c;
      for (size_t a = 0; a < p.v.size(); ++a)
        for (size_t b = a + 1; b < p.v.size(); ++b) c.push_back(dot(p.v[a], p.v[b]) / (norm(p.v[a]) * norm(p.v[b])));
      std::sort(c.begin(), c.end());
      return c;
    };
    std::map<int, std::vector<std::vector<double>>> shapes;
    for (size_t i = 0; i < bulk.atoms.size(); ++i) {
      const int z = bulk.atoms[i].element;
      if (all[i].v.empty()) continue;
      auto& L = polys[z];
      if (L.size() >= 24) continue;
      const auto sh = shape(all[i]);
      bool seen = false;
      for (const auto& o : shapes[z]) {
        if (o.size() != sh.size()) continue;
        double m = 0;
        for (size_t k = 0; k < sh.size(); ++k) m = std::max(m, std::fabs(sh[k] - o[k]));
        if (m < 0.02) { seen = true; break; }
      }
      if (seen) continue;
      shapes[z].push_back(sh);
      L.push_back(all[i]);
    }
  }
  const bool periodic = s.cell.valid();
  auto sep = [&](const Vec3& a, const Vec3& b) {
    Vec3 d = a - b;
    if (periodic) d = s.cell.minimum_image(d);
    return norm(d);
  };
  auto is_chalcogen = [](int z) { return z == 8 || z == 16 || z == 34; };
  // the closest an added atom may come to another atom: H 1.2 Å from H and 1.5 Å from heavier atoms, a heavy atom
  // 2.2 Å from heavy atoms (shorter than any O-O edge of an oxide polyhedron)
  auto too_close = [](int za, int zb) { return za == 1 && zb == 1 ? 1.2 : za == 1 || zb == 1 ? 1.5 : 2.2; };
  // planned atoms (all placed before any is added): owner, element, position
  struct Plan { uint32_t owner; int z; Vec3 p; };
  constexpr double kShared = 1.2;
  std::vector<Plan> plans;
  // The atoms near an atom i (i and its partners left out) and those planned for other atoms, as separations from i:
  // kind 1 a cation that lost bonds too (a new O at bond length from it is the site the two shared, not a clash),
  // kind 2 an O planned by another cation (two O planned on one lattice site, within 1.2 Å since no two O come that
  // close, are one bridging O).
  struct Near { Vec3 r; int z; int kind; double r2_skip; };
  std::vector<char> deficient(s.atoms.size(), 0);   // set below: atoms short of a bulk polyhedron
  auto lost_bonds = [&](uint32_t k) { return k < deficient.size() && deficient[k]; };
  auto near_of = [&](uint32_t i, int nz, size_t n_atoms) {
    std::vector<Near> out;
    const Vec3 pos = s.atoms[i].pos;
    for (uint32_t k = 0; k < n_atoms; ++k) {
      if (k == i || std::find(bonded_to[i].begin(), bonded_to[i].end(), k) != bonded_to[i].end()) continue;
      Vec3 d = s.atoms[k].pos - pos;
      if (periodic) d = s.cell.minimum_image(d);
      if (norm(d) >= 5.0) continue;
      const int zk = s.atoms[k].element;
      const bool shared = nz == 8 && !is_chalcogen(zk) && main_partner(zk) == 8 && lost_bonds(k);
      const double skip = shared ? 0.8 * bond_length(zk, 8) : 0;
      out.push_back({d, zk, shared ? 1 : 0, skip * skip});
    }
    for (const Plan& pl : plans) {
      if (pl.owner == i) continue;
      Vec3 d = pl.p - pos;
      if (periodic) d = s.cell.minimum_image(d);
      if (norm(d) < 5.0) out.push_back({d, pl.z, nz == 8 && pl.z == 8 ? 2 : 0, 0});
    }
    return out;
  };
  // clashes of an atom of element nz at r from atom i with the atoms near it
  auto clashes = [&](int nz, const Vec3& r, const std::vector<Near>& near) {
    int n = 0;
    for (const Near& q : near) {
      const Vec3 d = q.r - r;
      const double d2 = dot(d, d);
      if (q.kind == 1 && d2 > q.r2_skip) continue;
      if (q.kind == 2 && d2 < kShared * kShared) continue;
      const double c = too_close(nz, q.z);
      if (d2 < c * c) ++n;
    }
    return n;
  };
  // the direction from `dirs` with the largest smallest angle to `taken` (or that angle closest to `target`), not below
  // `floor`, pointing outward and clear of other atoms (`clash` counts them)
  auto pick = [&](const std::vector<Vec3>& taken, const Vec3& outward, double target, double floor, const std::function<int(const Vec3&)>& clash) {
    Vec3 best = outward;
    double bs = -1e300;
    for (const auto& u : dirs) {
      double mina = 180;
      for (const auto& t : taken) mina = std::min(mina, std::acos(std::clamp(dot(u, t), -1.0, 1.0)) / kDeg);
      double sc = target > 0 ? -std::fabs(mina - target) * 2 : mina;
      if (mina < floor) sc -= 1000;
      sc += 40 * dot(outward, u);
      if (clash) sc -= 500 * clash(u);
      if (sc > bs) { bs = sc; best = u; }
    }
    return best;
  };
  // the length of an added bond: X-H on O, S, Se; the bulk M-O length on a cation; C-H style on covalent atoms
  auto added_length = [&](int z, int pz, double tlen) {
    if (is_chalcogen(z)) return z == 8 ? 0.97 : 1.34;
    if (main_partner(z) == 8) return pz == 8 && tlen > 0.5 ? tlen : bond_length(z, 8);
    return element(z).covalent + 0.31;
  };
  auto unit = [](const Vec3& v) { return v * (1 / norm(v)); };
  // a rotation taking the frame of (t1, t2) onto the frame of (b1, b2)
  struct Rot {
    Vec3 b1, b2, b3, t1, t2, t3;
    Vec3 operator()(const Vec3& x) const { return b1 * dot(t1, x) + b2 * dot(t2, x) + b3 * dot(t3, x); }
  };
  auto frame = [](const Vec3& u, const Vec3& v, Vec3& e1, Vec3& e2, Vec3& e3) {
    e1 = u * (1 / norm(u));
    const Vec3 w = v - e1 * dot(e1, v);
    const double l = norm(w);
    if (l < 1e-6) return false;
    e2 = w * (1 / l), e3 = cross(e1, e2);
    return true;
  };
  // Every assignment of an atom's surviving bonds (unit vectors b) to vertices of polyhedron P whose angles agree within
  // 25°, a bounded search; `leaf` gets each assignment with the rotation taking P onto the bonds and its worst angle
  // (one surviving bond: the polyhedron turned about it in 5° steps).
  using Leaf = std::function<void(const std::vector<int>&, const Rot&, double)>;
  auto fit = [&](const std::vector<Vec3>& b, const Poly& P, const Leaf& leaf) {
    const size_t n = P.v.size(), k = b.size();
    if (k == 0 || k > n) return;
    std::vector<Vec3> tv(n);
    for (size_t t = 0; t < n; ++t) tv[t] = unit(P.v[t]);
    std::vector<int> asg(k, -1);
    std::vector<char> used(n, 0);
    long budget = 5000;
    std::function<void(size_t)> dfs = [&](size_t q) {
      if (--budget < 0) return;
      if (q == k) {
        Rot R{};
        if (k == 1) {
          const Vec3& t1 = tv[size_t(asg[0])];
          size_t ref = n;
          double lo = 2;
          for (size_t t = 0; t < n; ++t)
            if (int(t) != asg[0] && std::fabs(dot(tv[t], t1)) < lo) lo = std::fabs(dot(tv[t], t1)), ref = t;
          if (ref == n) { leaf(asg, R, 0.0); return; }   // a lone bond: nothing left to place
          if (lo > 0.95 || !frame(t1, tv[ref], R.t1, R.t2, R.t3)) return;
          Vec3 e1, perp, e3;
          if (!frame(b[0], std::fabs(b[0][0]) < 0.9 ? Vec3{1, 0, 0} : Vec3{0, 1, 0}, e1, perp, e3)) return;
          for (int st = 0; st < 72; ++st) {
            const double ph = st * 5.0 * kDeg;
            const Vec3 e2 = perp * std::cos(ph) + e3 * std::sin(ph);
            R.b1 = e1, R.b2 = e2, R.b3 = cross(e1, e2);
            leaf(asg, R, 0.0);
          }
          return;
        }
        // the frame of the two most nearly perpendicular surviving bonds
        size_t pa = 0, pb = 1;
        double lo = 2;
        for (size_t x = 0; x < k; ++x)
          for (size_t y = x + 1; y < k; ++y)
            if (std::fabs(dot(b[x], b[y])) < lo) lo = std::fabs(dot(b[x], b[y])), pa = x, pb = y;
        if (!frame(tv[size_t(asg[pa])], tv[size_t(asg[pb])], R.t1, R.t2, R.t3) || !frame(b[pa], b[pb], R.b1, R.b2, R.b3)) return;
        double resid = 0;
        for (size_t x = 0; x < k; ++x)
          resid = std::max(resid, std::acos(std::clamp(dot(unit(R(tv[size_t(asg[x])])), b[x]), -1.0, 1.0)) / kDeg);
        if (resid < 20) leaf(asg, R, resid);
        return;
      }
      for (size_t t = 0; t < n; ++t) {
        if (used[t]) continue;
        bool ok = true;
        for (size_t x = 0; x < q && ok; ++x) {
          const double ab = std::acos(std::clamp(dot(b[q], b[x]), -1.0, 1.0)), at = std::acos(std::clamp(dot(tv[t], tv[size_t(asg[x])]), -1.0, 1.0));
          ok = std::fabs(ab - at) < 25 * kDeg;
        }
        if (!ok) continue;
        used[t] = 1, asg[q] = int(t);
        dfs(q + 1);
        used[t] = 0, asg[q] = -1;
      }
    };
    dfs(0);
  };
  // Which atoms lost bonds: fewer than the largest bulk polyhedron of their element, and not a whole smaller one (its
  // angles within 15°), or with no polyhedron at all fewer bonds than the element's usual coordination.
  const size_t n0 = s.atoms.size();
  for (size_t i = 0; i < n0; ++i) {
    const int z = s.atoms[i].element;
    if (!passivable(z)) continue;
    const size_t k = bonded[i].size();
    const auto pit = polys.find(z);
    if (pit == polys.end() || pit->second.empty()) {
      deficient[i] = int(k) < (want.count(z) ? want[z] : 0);
      continue;
    }
    size_t max_n = 0;
    for (const Poly& P : pit->second) max_n = std::max(max_n, P.v.size());
    if (k >= max_n) continue;
    bool whole = false;
    for (const Poly& P : pit->second) {
      if (P.v.size() != k || whole) continue;
      if (k == 0) { whole = true; break; }
      fit(bonded[i], P, [&](const std::vector<int>&, const Rot&, double resid) { if (resid < 15) whole = true; });
    }
    deficient[i] = !whole;
  }
  int greedy = 0;
  for (size_t i = 0; i < n0; ++i) {
    const int z = s.atoms[i].element;
    if (!deficient[i]) continue;
    const int nz = is_chalcogen(z) ? 1 : main_partner(z) == 8 ? 8 : 1;
    const Vec3 pos = s.atoms[i].pos;
    const Vec3 outward = outward_of(pos);
    const auto near = near_of(uint32_t(i), nz, n0);
    // each larger bulk polyhedron fitted to the surviving bonds: the best fit (smallest angle error, fewest clashes,
    // most outward) gives the new atoms
    const auto& b = bonded[i];
    const size_t k = b.size();
    struct Fit { double score = -1e300; std::vector<Vec3> p; } best;
    if (polys.count(z))
      for (const Poly& P : polys[z]) {
        if (P.v.size() <= k) continue;
        fit(b, P, [&](const std::vector<int>& asg, const Rot& R, double resid) {
          std::vector<char> used(P.v.size(), 0);
          for (int a : asg) used[size_t(a)] = 1;
          std::vector<Vec3> placed;
          int nclash = 0;
          double out = 0;
          for (size_t t = 0; t < P.v.size(); ++t) {
            if (used[t]) continue;
            const Vec3 u = unit(R(P.v[t]));
            const Vec3 r = u * added_length(z, P.pz[t], norm(P.v[t]));
            nclash += clashes(nz, r, near);
            out += dot(outward, u);
            placed.push_back(pos + r);
          }
          if (placed.empty()) return;
          const double sc = -resid - 60.0 * nclash + 10.0 * out / double(placed.size());
          if (sc > best.score) best.score = sc, best.p = std::move(placed);
        });
      }
    const int missing = (want.count(z) ? want[z] : 0) - int(k);
    if (best.p.empty() && missing <= 0) continue;
    std::vector<Vec3> placed;
    if (!best.p.empty()) {
      placed = best.p;
    } else {
      // no polyhedron fits (or none survives): one bond at a time, no closer than the polyhedron's own angles allow
      double target = 0, floor = 95;
      if (polys.count(z) && !polys[z].empty()) {
        const Poly* P = &polys[z].front();
        for (const Poly& Q : polys[z])
          if (int(Q.v.size()) == missing + int(k)) { P = &Q; break; }
        const auto& v = P->v;
        std::vector<double> nearest;
        for (size_t x = 0; x < v.size(); ++x) {
          double m = 180;
          for (size_t y = 0; y < v.size(); ++y)
            if (y != x) m = std::min(m, std::acos(std::clamp(dot(unit(v[x]), unit(v[y])), -1.0, 1.0)) / kDeg);
          nearest.push_back(m);
        }
        std::sort(nearest.begin(), nearest.end());
        target = nearest[nearest.size() / 2];
        floor = std::max(60.0, nearest.front() - 15);
      }
      const double len = added_length(z, 8, bond_length(z, 8));
      auto taken = b;
      for (int m = 0; m < missing; ++m) {
        const Vec3 u = pick(taken, outward, k + size_t(m) >= 2 ? target : 0, floor, [&](const Vec3& c) {
          const Vec3 p = pos + c * len;
          int n = clashes(nz, c * len, near);
          for (const Vec3& q : placed) n += sep(p, q) < too_close(nz, nz);
          return n;
        });
        taken.push_back(u);
        placed.push_back(pos + u * len);
      }
      ++greedy;
    }
    for (const Vec3& p : placed) plans.push_back({uint32_t(i), nz, p});
  }
  // the planned atoms added: an O planned by two or more cations on one site is the bulk's bridging O put back, bonded
  // to all of them; every new O with fewer cations than the bulk O has gets one H (M-OH at 115° to the M-O bond, a
  // bridging M-OH-M away from its cations), outward and clear of the others
  std::vector<uint32_t> new_o;
  std::map<uint32_t, std::vector<uint32_t>> o_owners;
  for (const Plan& pl : plans) {
    if (pl.z == 1) {
      const uint32_t hdx = add_atom(s, 1, pl.p, "H", 1);
      s.bonds.push_back({pl.owner, hdx, 1});
      ++added_h;
      continue;
    }
    long at = -1;
    for (uint32_t o : new_o)
      if (sep(pl.p, s.atoms[o].pos) < kShared) { at = long(o); break; }
    if (at < 0) {
      at = long(add_atom(s, 8, pl.p, "O", 1));
      new_o.push_back(uint32_t(at));
    } else {
      // the shared site: the mean of the positions its cations planned
      auto& a = s.atoms[size_t(at)];
      Vec3 d = pl.p - a.pos;
      if (periodic) d = s.cell.minimum_image(d);
      a.pos = a.pos + d * (1.0 / double(o_owners[uint32_t(at)].size() + 1));
    }
    o_owners[uint32_t(at)].push_back(pl.owner);
    s.bonds.push_back({pl.owner, uint32_t(at), 1});
  }
  int bridging = 0;
  const int o_cn = want.count(8) ? want[8] : 2;
  for (uint32_t o : new_o) {
    const auto& own = o_owners[o];
    if (own.size() > 1) ++bridging;
    if (int(own.size()) >= o_cn) continue;   // as many cations as in the bulk: a whole O
    const Vec3 po = s.atoms[o].pos;
    std::vector<Vec3> to_m;
    for (uint32_t m : own) {
      Vec3 d = s.atoms[m].pos - po;
      if (periodic) d = s.cell.minimum_image(d);
      to_m.push_back(unit(d));
    }
    std::vector<uint32_t> near;
    for (uint32_t q = 0; q < s.atoms.size(); ++q)
      if (q != o && std::find(own.begin(), own.end(), q) == own.end() && sep(po, s.atoms[q].pos) < 4.0) near.push_back(q);
    const Vec3 v = pick(to_m, outward_of(po), own.size() == 1 ? 115 : 0, 95, [&](const Vec3& c) {
      const Vec3 p = po + c * 0.97;
      int n = 0;
      for (uint32_t q : near) n += sep(p, s.atoms[q].pos) < too_close(1, s.atoms[q].element);
      return n;
    });
    const uint32_t hdx = add_atom(s, 1, po + v * 0.97, "H", 1);
    s.bonds.push_back({o, hdx, 1});
    ++added_oh;
  }
  if (bridging) s.notes.push_back(std::to_string(bridging) + " bridging O put back on bulk sites shared by surface cations that lost them");
  if (greedy) s.notes.push_back(std::to_string(greedy) + " surface atoms matched no bulk polyhedron: their missing bonds placed one at a time");
  // Each direction above was chosen from its own atom's bonds only: a hydrogen can land on a neighbour's (two silanols
  // on adjacent oxygens). Every added hydrogen closer than 1.7 Å to an atom it is not bonded to turns about its
  // parent's bond axis (length and angle kept) to the position farthest from the others.
  const auto nb = s.neighbours();
  auto closest = [&](uint32_t h, const Vec3& p) {
    double m = 1e300;
    for (uint32_t k = 0; k < s.atoms.size(); ++k) {
      if (k == h || (nb[h].size() == 1 && k == nb[h][0])) continue;
      Vec3 d = p - s.atoms[k].pos;
      if (s.cell.valid()) d = s.cell.minimum_image(d);
      m = std::min(m, norm(d));
    }
    return m;
  };
  for (uint32_t h = uint32_t(n0); h < s.atoms.size(); ++h) {
    if (s.atoms[h].element != 1 || nb[h].size() != 1) continue;
    if (closest(h, s.atoms[h].pos) >= 1.7) continue;
    const uint32_t par = nb[h][0];
    int anchor = -1;
    for (uint32_t q : nb[par]) if (q != h) { anchor = int(q); break; }
    if (anchor < 0) continue;
    Vec3 axis = s.atoms[par].pos - s.atoms[uint32_t(anchor)].pos;
    if (s.cell.valid()) axis = s.cell.minimum_image(axis);
    axis = axis * (1 / norm(axis));
    Vec3 arm = s.atoms[h].pos - s.atoms[par].pos;
    if (s.cell.valid()) arm = s.cell.minimum_image(arm);
    Vec3 best = s.atoms[h].pos;
    double bd = closest(h, best);
    for (int k = 1; k < 36; ++k) {
      const double t = k * 10.0 * kDeg, c = std::cos(t), sn = std::sin(t);
      const Vec3 r = arm * c + cross(axis, arm) * sn + axis * (dot(axis, arm) * (1 - c));
      const Vec3 p = s.atoms[par].pos + r;
      const double dd = closest(h, p);
      if (dd > bd) bd = dd, best = p;
    }
    s.atoms[h].pos = best;
  }
  return {added_h, added_oh};
}

// ---------------------------------------------------------------- CIF

System parse_cif(const std::string& text, const std::string& name, const CifDisorder& o, CifDisorderReport* report) {
  CifDisorderReport R;
  const CifData d = cif_parse(text);
  auto item = [&](const std::string& k) -> std::string {
    auto it = d.items.find(k);
    return it == d.items.end() ? std::string() : it->second;
  };
  const double A = cif_number(item("_cell_length_a")), B = cif_number(item("_cell_length_b")), C = cif_number(item("_cell_length_c"));
  const double al = cif_number(item("_cell_angle_alpha"), 90) * kDeg, be = cif_number(item("_cell_angle_beta"), 90) * kDeg,
               ga = cif_number(item("_cell_angle_gamma"), 90) * kDeg;
  if (A <= 0 || B <= 0 || C <= 0) throw CrystalError((name.empty() ? std::string("CIF") : name) + ": no cell lengths (_cell_length_a/b/c)");
  System s;
  s.source_format = "cif";
  s.title = !item("_chemical_name_mineral").empty() ? item("_chemical_name_mineral") : !item("_chemical_name_common").empty() ? item("_chemical_name_common") : d.name;
  s.cell.a = {A, 0, 0};
  s.cell.b = {B * std::cos(ga), B * std::sin(ga), 0};
  const double cx = C * std::cos(be), cy = C * (std::cos(al) - std::cos(be) * std::cos(ga)) / std::sin(ga);
  s.cell.c = {cx, cy, std::sqrt(std::max(0.0, C * C - cx * cx - cy * cy))};
  // symmetry
  std::vector<SymOp> ops;
  std::vector<std::vector<std::string>> sites;
  std::vector<std::string> site_tags;
  for (const auto& [tags, rows] : d.loops) {
    for (size_t t = 0; t < tags.size(); ++t)
      if (tags[t] == "_symmetry_equiv_pos_as_xyz" || tags[t] == "_space_group_symop_operation_xyz")
        for (const auto& r : rows) ops.push_back(parse_symop(r[t]));
    if (std::find(tags.begin(), tags.end(), "_atom_site_fract_x") != tags.end()) { site_tags = tags; sites = rows; }
  }
  if (ops.empty()) {
    for (const char* k : {"_symmetry_equiv_pos_as_xyz", "_space_group_symop_operation_xyz"})
      if (!item(k).empty()) ops.push_back(parse_symop(item(k)));
    if (ops.empty()) ops.push_back(parse_symop("x,y,z"));
    if (item("_symmetry_space_group_name_h-m").size() > 0 && lower(item("_symmetry_space_group_name_h-m")) != "p1" && lower(item("_symmetry_space_group_name_h-m")) != "p 1")
      s.notes.push_back("the file names space group " + item("_symmetry_space_group_name_h-m") + " but lists no symmetry operations; read as P1");
  }
  if (sites.empty()) throw CrystalError((name.empty() ? std::string("CIF") : name) + ": no atom sites with fractional coordinates");
  auto col = [&](const char* tag) {
    auto it = std::find(site_tags.begin(), site_tags.end(), tag);
    return it == site_tags.end() ? -1 : int(it - site_tags.begin());
  };
  const int cx_ = col("_atom_site_fract_x"), cy_ = col("_atom_site_fract_y"), cz_ = col("_atom_site_fract_z");
  if (cy_ < 0 || cz_ < 0)
    throw CrystalError((name.empty() ? std::string("CIF") : name) + ": the atom sites have no _atom_site_fract_" + (cy_ < 0 ? "y" : "z") + " column");
  const int csym = col("_atom_site_type_symbol"), clab = col("_atom_site_label"), cocc = col("_atom_site_occupancy");
  const int cass = col("_atom_site_disorder_assembly"), cgrp = col("_atom_site_disorder_group");
  const int cq = col("_atom_site_charge");   // partial charges (as CAPS writes them), copied to every image of a site
  if (o.mode != "major" && o.mode != "group" && o.mode != "random") throw CrystalError("disorder: major, group or random");
  // the sites
  struct Site { std::string label, assembly, group; int z = 0; Vec3 f0{0, 0, 0}; double occ = 1, q = 0; };
  std::vector<Site> S;
  auto blank = [](const std::string& v) { return v.empty() || v == "." || v == "?"; };
  for (const auto& r : sites) {
    Site st;
    st.occ = cocc >= 0 ? cif_number(r[size_t(cocc)], 1) : 1;
    st.z = csym >= 0 ? element_of_site(r[size_t(csym)]) : 0;
    if (!st.z && clab >= 0) st.z = element_of_site(r[size_t(clab)]);
    if (!st.z) throw CrystalError("site " + (clab >= 0 ? r[size_t(clab)] : std::string("?")) + ": unknown element");
    st.label = clab >= 0 ? r[size_t(clab)] : element(st.z).symbol;
    st.f0 = {cif_number(r[size_t(cx_)]), cif_number(r[size_t(cy_)]), cif_number(r[size_t(cz_)])};
    st.q = cq >= 0 && !blank(r[size_t(cq)]) ? cif_number(r[size_t(cq)], 0) : 0;
    if (cgrp >= 0 && !blank(r[size_t(cgrp)]) && r[size_t(cgrp)] != "0") {
      st.group = r[size_t(cgrp)];
      // a negative group (a fragment disordered across a symmetry element) is no alternative to the positive groups
      // of the same file: without a named assembly it gets one of its own, so neither choice nor draw drops the other
      st.assembly = cass >= 0 && !blank(r[size_t(cass)]) ? r[size_t(cass)] : st.group[0] == '-' ? "S" : "A";
    }
    S.push_back(st);
  }
  // the assemblies' groups and their occupancies (the largest of each group's sites)
  std::map<std::string, std::map<std::string, double>> asm_groups;
  for (const auto& st : S)
    if (!st.group.empty()) {
      double& g = asm_groups[st.assembly][st.group];
      g = std::max(g, st.occ);
    }
  R.assemblies = int(asm_groups.size());
  for (const auto& [a, gs] : asm_groups) R.groups += int(gs.size());
  std::map<std::string, std::string> chosen;   // major / group modes
  for (const auto& [a, gs] : asm_groups) {
    std::string best;
    double bo = -1;
    for (const auto& [g, oc] : gs)
      if (oc > bo + 1e-9) bo = oc, best = g;
    if (o.mode == "group" && gs.count(o.group)) best = o.group;
    else if (o.mode == "group") R.notes.push_back("assembly " + a + " has no group " + o.group + ": its major group " + best + " kept");
    chosen[a] = best;
  }
  // a negative disorder group: a fragment disordered across a symmetry element. The file lists one orientation and
  // the operations make the others, which overlap it, so of each overlapping set of copies only one is kept.
  auto across = [](const Site& st) { return !st.group.empty() && st.group[0] == '-'; };
  // every site through every operation: unwrapped (for which copy an atom belongs to), wrapped, one per position.
  // An image that falls on an earlier one of the same site (an atom on the symmetry element) is shared: for a site
  // across a symmetry element its operation and copy shift are kept as an alias, so the atom stays with either copy.
  struct Ex { size_t site; int op; Vec3 fw; std::array<int, 3> shift; std::vector<std::pair<int, std::array<int, 3>>> alias{}; };
  std::vector<Ex> X;
  for (size_t si = 0; si < S.size(); ++si) {
    const size_t first = X.size();
    for (size_t k = 0; k < ops.size(); ++k) {
      Vec3 f;
      for (int i = 0; i < 3; ++i) f[i] = ops[k].R[i][0] * S[si].f0[0] + ops[k].R[i][1] * S[si].f0[1] + ops[k].R[i][2] * S[si].f0[2] + ops[k].t[i];
      const Vec3 fw = wrap01(f);
      long dup = -1;
      for (size_t q = first; q < X.size() && dup < 0; ++q) {
        Vec3 df = fw - X[q].fw;
        for (auto& x : df) x -= std::round(x);
        if (norm(s.cell.to_cartesian(df) - s.cell.to_cartesian({0, 0, 0})) < 0.05) dup = long(q);
      }
      if (dup >= 0) {
        if (across(S[si])) {
          Ex& e = X[size_t(dup)];
          e.alias.push_back({int(k), {int(std::lround(f[0] - e.fw[0])), int(std::lround(f[1] - e.fw[1])), int(std::lround(f[2] - e.fw[2]))}});
        }
        continue;
      }
      X.push_back({si, int(k), fw, {int(std::lround(f[0] - fw[0])), int(std::lround(f[1] - fw[1])), int(std::lround(f[2] - fw[2]))}});
    }
  }
  // places: partially occupied atoms without a group within 0.8 Å of each other (and an occupied atom with them)
  std::vector<int> place(X.size(), -1);
  int nplaces = 0;
  auto dist = [&](const Vec3& a, const Vec3& b) {
    Vec3 df = a - b;
    for (auto& x : df) x -= std::round(x);
    return norm(s.cell.to_cartesian(df) - s.cell.to_cartesian({0, 0, 0}));
  };
  for (size_t i = 0; i < X.size(); ++i) {
    const Site& a = S[X[i].site];
    if (!a.group.empty() || a.occ >= 0.999 || place[i] >= 0) continue;
    place[i] = nplaces;
    std::vector<size_t> q{i};
    for (size_t h = 0; h < q.size(); ++h)
      for (size_t j = 0; j < X.size(); ++j)
        if (place[j] < 0 && S[X[j].site].group.empty() && S[X[j].site].occ < 0.999 && dist(X[q[h]].fw, X[j].fw) < 0.8) place[j] = nplaces, q.push_back(j);
    ++nplaces;
  }
  std::vector<std::vector<size_t>> members(static_cast<size_t>(nplaces));
  for (size_t i = 0; i < X.size(); ++i)
    if (place[i] >= 0) members[size_t(place[i])].push_back(i);
  R.places = nplaces;
  // the major occupant of each place
  std::vector<long> major_of(static_cast<size_t>(nplaces), -1);
  int split_kept = 0, mixed = 0;
  for (int p = 0; p < nplaces; ++p) {
    const auto& m = members[size_t(p)];
    double tot = 0, bo = -1;
    std::set<int> els;
    for (size_t i : m) {
      const Site& st = S[X[i].site];
      tot += st.occ, els.insert(st.z);
      if (st.occ > bo + 1e-9) bo = st.occ, major_of[size_t(p)] = long(i);
    }
    if (tot < 0.5 - 1e-9) major_of[size_t(p)] = -1;
    if (m.size() > 1) els.size() > 1 ? ++mixed : ++split_kept;
  }
  // the cells, the draws
  const bool random = o.mode == "random";
  const std::array<int, 3> N = o.supercell;   // every mode: major and group repeat the one chosen cell
  for (int k = 0; k < 3; ++k)
    if (N[size_t(k)] < 1) throw CrystalError("disorder: a supercell of at least 1 × 1 × 1");
  std::mt19937_64 rng(o.seed * 0x9E3779B97F4A7C15ull + 17);
  auto u01 = [&] { return double(rng() >> 11) * (1.0 / 9007199254740992.0); };
  std::map<std::tuple<std::string, int, int, int, int>, std::string> drawn;   // (assembly, op, copy cell) → group
  auto draw_group = [&](const std::string& a, int op, int c0, int c1, int c2) {
    const auto key = std::make_tuple(a, op, c0, c1, c2);
    auto it = drawn.find(key);
    if (it != drawn.end()) return it->second;
    double u = u01(), c = 0;
    std::string g;   // none: the assembly absent from this copy (occupancies adding up to less than 1)
    for (const auto& [name, oc] : asm_groups[a])
      if ((c += oc) > u) { g = name; break; }
    return drawn[key] = g;
  };
  auto cmod = [](int k, int n) { return ((k % n) + n) % n; };
  const int ncells = N[0] * N[1] * N[2];
  auto cell_index = [&](int c0, int c1, int c2) { return (cmod(c0, N[0]) * N[1] + cmod(c1, N[1])) * N[2] + cmod(c2, N[2]); };
  // copies across a symmetry element: one node per (assembly, group, operation) and copy cell. Two copies conflict
  // when atoms of the two come closer than 1.4 Å (shorter than any contact between molecules that are both there).
  struct Owner { int node; std::array<int, 3> shift; };
  std::vector<std::vector<Owner>> owners(X.size());
  std::map<std::tuple<std::string, std::string, int>, int> node_id;
  std::vector<std::pair<std::string, std::string>> node_ag;   // (assembly, group) of each node
  for (size_t i = 0; i < X.size(); ++i) {
    const Site& st = S[X[i].site];
    if (!across(st)) continue;
    auto own = [&](int op, const std::array<int, 3>& shift) {
      auto [it, added] = node_id.emplace(std::make_tuple(st.assembly, st.group, op), int(node_ag.size()));
      if (added) node_ag.push_back({st.assembly, st.group});
      owners[i].push_back({it->second, shift});
    };
    own(X[i].op, X[i].shift);
    for (const auto& [op, shift] : X[i].alias) own(op, shift);
  }
  std::vector<char> copy_on;   // per node × copy cell: 1 kept, 2 left out
  int across_left_out = 0;
  if (!node_ag.empty()) {
    // conflicts: node a in copy cell c overlaps node b in copy cell c + D
    std::vector<std::set<std::pair<int, std::array<int, 3>>>> conflict(node_ag.size());
    for (size_t i = 0; i < X.size(); ++i)
      for (size_t j = i + 1; j < X.size(); ++j) {
        if (owners[i].empty() || owners[j].empty() || S[X[i].site].assembly != S[X[j].site].assembly) continue;
        if (dist(X[i].fw, X[j].fw) >= 1.4) continue;
        for (const Owner& a : owners[i])
          for (const Owner& b : owners[j]) {
            if (a.node == b.node) continue;
            std::array<int, 3> D;
            for (int k = 0; k < 3; ++k)
              D[size_t(k)] = int(std::lround(X[i].fw[k] + a.shift[size_t(k)] - X[j].fw[k] - b.shift[size_t(k)]));
            conflict[size_t(a.node)].insert({b.node, D});
            conflict[size_t(b.node)].insert({a.node, {-D[0], -D[1], -D[2]}});
          }
      }
    copy_on.assign(node_ag.size() * size_t(ncells), 0);
    auto neighbours_of = [&](int v) {
      const int node = v / ncells, c = v % ncells, c0 = c / (N[1] * N[2]), c1 = (c / N[2]) % N[1], c2 = c % N[2];
      std::vector<int> out;
      for (const auto& [b, D] : conflict[size_t(node)]) {
        const int w = b * ncells + cell_index(c0 + D[0], c1 + D[1], c2 + D[2]);
        if (w != v && std::find(out.begin(), out.end(), w) == out.end()) out.push_back(w);
      }
      return out;
    };
    auto blocked = [&](int v) {
      for (int w : neighbours_of(v))
        if (copy_on[size_t(w)] == 1) return true;
      return false;
    };
    if (!random) {
      // the chosen group's copies in the order of the operations: a copy stays unless it overlaps one already kept
      // (the same choice in every copy cell; left-out copies counted per cell)
      for (const auto& [key, node] : node_id) {
        const auto& [a, g] = node_ag[size_t(node)];
        for (int c = 0; c < ncells; ++c) {
          const int v = node * ncells + c;
          const bool keep_it = chosen[a] == g && !blocked(v);
          copy_on[size_t(v)] = keep_it ? 1 : 2;
          if (c == 0 && chosen[a] == g && !keep_it) ++across_left_out;
        }
      }
    } else {
      // in a random order, each undecided copy and the free copies it overlaps: at most one of them, drawn by occupancy
      std::vector<int> order(copy_on.size());
      std::iota(order.begin(), order.end(), 0);
      for (size_t k = order.size(); k > 1; --k) std::swap(order[k - 1], order[std::min(k - 1, size_t(u01() * double(k)))]);
      for (int v : order) {
        if (copy_on[size_t(v)]) continue;
        if (blocked(v)) { copy_on[size_t(v)] = 2; continue; }
        std::vector<int> Q{v};
        for (int w : neighbours_of(v))
          if (!copy_on[size_t(w)] && !blocked(w)) Q.push_back(w);
        double u = u01(), c = 0;
        int pick = -1;
        for (int w : Q) {
          const auto& [a, g] = node_ag[size_t(w / ncells)];
          if ((c += asm_groups[a][g]) > u) { pick = w; break; }
        }
        for (int w : Q) copy_on[size_t(w)] = w == pick ? 1 : 2;
      }
    }
  }
  struct Kept { int z; Vec3 f; std::string label; double q; };
  std::vector<Kept> keep;
  for (int i0 = 0; i0 < N[0]; ++i0)
    for (int i1 = 0; i1 < N[1]; ++i1)
      for (int i2 = 0; i2 < N[2]; ++i2) {
        // each place's occupant in this cell
        std::vector<long> occ_of(static_cast<size_t>(nplaces), -1);
        for (int p = 0; p < nplaces; ++p) {
          if (!random) { occ_of[size_t(p)] = major_of[size_t(p)]; continue; }
          double u = u01(), c = 0;
          for (size_t i : members[size_t(p)])
            if ((c += S[X[i].site].occ) > u) { occ_of[size_t(p)] = long(i); break; }
        }
        for (size_t i = 0; i < X.size(); ++i) {
          const Site& st = S[X[i].site];
          bool on;
          if (across(st)) {
            on = false;
            for (const Owner& w : owners[i])
              if (copy_on[size_t(w.node) * size_t(ncells) + size_t(cell_index(i0 - w.shift[0], i1 - w.shift[1], i2 - w.shift[2]))] == 1) { on = true; break; }
          } else if (!st.group.empty()) {
            if (!random) on = chosen[st.assembly] == st.group;
            else {
              const int c0 = cmod(i0 - X[i].shift[0], N[0]), c1 = cmod(i1 - X[i].shift[1], N[1]), c2 = cmod(i2 - X[i].shift[2], N[2]);
              on = draw_group(st.assembly, X[i].op, c0, c1, c2) == st.group;
            }
          } else if (place[i] >= 0) {
            on = occ_of[size_t(place[i])] == long(i);
          } else {
            on = true;
          }
          for (const auto& [z, e] : std::map<int, double>{{st.z, st.occ}}) R.expected[z] += e;
          if (on) keep.push_back({st.z, {(X[i].fw[0] + i0) / N[0], (X[i].fw[1] + i1) / N[1], (X[i].fw[2] + i2) / N[2]}, st.label, st.q});
        }
      }
  s.cell.a = s.cell.a * double(N[0]), s.cell.b = s.cell.b * double(N[1]), s.cell.c = s.cell.c * double(N[2]);
  if (random) {
    R.notes.push_back("disorder filled at random in a " + std::to_string(N[0]) + " × " + std::to_string(N[1]) + " × " + std::to_string(N[2]) + " supercell (seed " +
                      std::to_string(o.seed) + ")");
  } else {
    if (R.assemblies) {
      std::string t;
      for (const auto& [a, g] : chosen) t += (t.empty() ? "" : ", ") + std::string("assembly ") + a + " group " + g;
      R.notes.push_back(std::to_string(R.assemblies) + " disorder assembl" + (R.assemblies == 1 ? "y" : "ies") + ": kept " + t);
    }
    if (across_left_out)
      R.notes.push_back(std::to_string(across_left_out) + " overlapping cop" + (across_left_out == 1 ? "y" : "ies") +
                        " of fragments disordered across a symmetry element left out: one orientation kept");
    int lone = 0;
    for (int p = 0; p < nplaces; ++p) lone += major_of[size_t(p)] < 0;
    if (lone) s.notes.push_back(std::to_string(lone) + " sites with occupancy below 0.5 left out");
    if (mixed) R.notes.push_back(std::to_string(mixed) + " mixed sites: the element of the largest occupancy kept");
    if (split_kept) R.notes.push_back(std::to_string(split_kept) + " split positions: the one of the largest occupancy kept");
    if (ncells > 1)
      R.notes.push_back("a " + std::to_string(N[0]) + " × " + std::to_string(N[1]) + " × " + std::to_string(N[2]) + " supercell of the ordered cell");
  }
  std::map<int, int> type_of;
  for (size_t i = 0; i < keep.size(); ++i) {
    Atom a;
    a.id = int64_t(i + 1);
    a.mol = 1;
    a.element = keep[i].z;
    a.name = keep[i].label;
    a.charge = keep[i].q;
    a.pos = s.cell.a * keep[i].f[0] + s.cell.b * keep[i].f[1] + s.cell.c * keep[i].f[2];
    auto it = type_of.find(keep[i].z);
    if (it == type_of.end()) {
      it = type_of.emplace(keep[i].z, int(type_of.size()) + 1).first;
      TypeInfo ti;
      ti.type = it->second;
      ti.mass = element(keep[i].z).mass;
      ti.label = element(keep[i].z).symbol;
      s.types.push_back(ti);
    }
    a.type = it->second;
    s.atoms.push_back(a);
    R.obtained[keep[i].z] += 1;
  }
  s.has_charges = cq >= 0;
  // a CIF numbers no molecules: they follow the perceived bonds (a network is one, a molecular crystal one per unit,
  // a polymer cell one per chain); mol = 1 above is only a placeholder
  s.has_mol = false;
  s.bonds = crystal_bonds(s);
  // A CIF holds fractional coordinates wrapped into the cell, so a finite molecule that crosses a cell face arrives cut
  // in pieces. Each finite molecule is made whole here along its image bonds (Rg, end-to-end vectors and C_n are only
  // defined on whole chains); a periodic network keeps its positions as placed, since no choice of images makes it
  // whole. This is per molecule, so guests in a framework are made whole while the framework stays as read.
  make_finite_molecules_whole(s);
  s.unwrapped = true;
  std::map<int, int> count;
  for (const auto& a : s.atoms) count[a.element]++;
  s.notes.insert(s.notes.begin(), std::to_string(ops.size()) + " symmetry operations · " + std::to_string(s.atoms.size()) + " atoms in the " +
                                      (ncells > 1 ? "supercell" : "cell") + " (" + hill(count) + ") · " + std::to_string(s.bonds.size()) + " bonds");
  for (const auto& n : R.notes) s.notes.push_back(n);
  // the space group stays with an ordered cell that its operations still map onto itself
  if (ops.size() > 1 && !random && ncells == 1 && R.groups == 0 && nplaces == 0) {
    Symmetry g;
    for (const char* k : {"_space_group_name_h-m_alt", "_symmetry_space_group_name_h-m"})
      if (g.name.empty() && !item(k).empty() && item(k) != "?") g.name = item(k);
    for (const char* k : {"_space_group_it_number", "_symmetry_int_tables_number"})
      if (g.number == 0 && !item(k).empty()) g.number = int(cif_number(item(k), 0));
    if (g.number == 0 && !g.name.empty())
      if (const auto* sg = find_space_group(g.name)) g.number = sg->number;   // the number from the symbol
    g.ops = ops;
    if (symmetry_holds(s, g, 0.05)) s.symmetry = std::move(g);
    else s.notes.push_back("the symmetry operations do not map the atoms read onto themselves: kept as P 1");
  }
  if (report) *report = R;
  return s;
}

System parse_cif(const std::string& text, const std::string& name) { return parse_cif(text, name, CifDisorder{}, nullptr); }

System read_cif(const std::string& path, const CifDisorder& o, CifDisorderReport* report) {
  std::ifstream in(path);
  if (!in) throw CrystalError("cannot open " + path);
  std::stringstream ss;
  ss << in.rdbuf();
  return parse_cif(ss.str(), path, o, report);
}

System read_cif(const std::string& path) { return read_cif(path, CifDisorder{}, nullptr); }

CifDisplacements cif_displacements(const std::string& text) {
  const CifData d = cif_parse(text);
  CifDisplacements out;
  constexpr double k8pi2 = 8 * 3.14159265358979323846 * 3.14159265358979323846;
  for (const auto& [tags, rows] : d.loops) {
    auto col = [&](const char* tag) {
      auto it = std::find(tags.begin(), tags.end(), tag);
      return it == tags.end() ? -1 : int(it - tags.begin());
    };
    const int lab = col("_atom_site_label"), ui = col("_atom_site_u_iso_or_equiv"), bi = col("_atom_site_b_iso_or_equiv");
    if (lab >= 0 && (ui >= 0 || bi >= 0))
      for (const auto& r : rows) {
        const std::string& v = r[size_t(ui >= 0 ? ui : bi)];
        if (v == "." || v == "?") continue;
        out.u_iso[r[size_t(lab)]] = ui >= 0 ? cif_number(v) : cif_number(v) / k8pi2;
      }
    const int al = col("_atom_site_aniso_label");
    if (al < 0) continue;
    const bool isB = col("_atom_site_aniso_b_11") >= 0 && col("_atom_site_aniso_u_11") < 0;
    const char* names[6] = {"11", "22", "33", "12", "13", "23"};
    int c[6];
    bool all = true;
    for (int k = 0; k < 6; ++k) {
      c[k] = col((std::string(isB ? "_atom_site_aniso_b_" : "_atom_site_aniso_u_") + names[k]).c_str());
      all = all && c[k] >= 0;
    }
    if (!all) continue;
    for (const auto& r : rows) {
      std::array<double, 6> u{};
      for (int k = 0; k < 6; ++k) u[size_t(k)] = cif_number(r[size_t(c[k])]) / (isB ? k8pi2 : 1.0);
      out.u_aniso[r[size_t(al)]] = u;
    }
  }
  return out;
}

std::string formula_of(const System& s, const std::vector<size_t>& atoms) {
  std::map<int, int> count;
  for (size_t i : atoms) count[s.atoms[i].element]++;
  return hill(count);
}

// ---------------------------------------------------------------- surfaces

std::vector<Termination> slab_terminations(const System& bulk, int h, int k, int l, double* dout) {
  const PlaneSetup P = plane_setup(bulk, h, k, l);
  if (dout) *dout = P.d;
  const auto planes = atomic_planes(P);
  const auto bonds = image_bonds(bulk);
  std::vector<Termination> out;
  const size_t np = planes.size();
  for (size_t i = 0; i < np; ++i) {
    const Plane& lo = planes[i];
    const Plane& hi = planes[(i + 1) % np];
    const double zhi = i + 1 < np ? hi.z : hi.z + P.d;
    Termination t;
    t.gap = zhi - lo.z;
    t.cut = std::fmod(lo.z + t.gap / 2, P.d);
    // bulk bonds crossing the cut (at any layer)
    for (const auto& b : bonds) {
      const double z1 = P.f[b.i][2] * P.d, z2 = z1 + dot(b.d, P.n);
      const double a = std::min(z1, z2), c = std::max(z1, z2);
      const long m0 = long(std::floor((a - t.cut) / P.d)) + 1, m1 = long(std::ceil((c - t.cut) / P.d)) - 1;
      for (long m = m0; m <= m1; ++m) {
        const double zc = t.cut + double(m) * P.d;
        if (zc > a + 1e-9 && zc < c - 1e-9) ++t.bonds_cut;
      }
    }
    t.bonds_per_nm2 = t.bonds_cut / (P.area / 100.0);
    t.top = formula_of(bulk, lo.atoms);
    t.bottom = formula_of(bulk, hi.atoms);
    char b[160];
    std::snprintf(b, sizeof b, "%s-terminated · %.2f Å gap · %.1f bonds/nm² cut", t.top.c_str(), t.gap, t.bonds_per_nm2);
    t.label = b;
    if (t.bottom != t.top) t.label += " · bottom " + t.bottom;
    out.push_back(t);
  }
  std::stable_sort(out.begin(), out.end(), [](const Termination& x, const Termination& y) {
    if (x.bonds_cut != y.bonds_cut) return x.bonds_cut < y.bonds_cut;
    return x.gap > y.gap;
  });
  // cuts related by the crystal's symmetry (screw axes, centring) look the same: keep one of each
  std::vector<Termination> uniq;
  for (const auto& t : out) {
    bool seen = false;
    for (const auto& u : uniq)
      seen = seen || (u.top == t.top && u.bottom == t.bottom && u.bonds_cut == t.bonds_cut && std::fabs(u.gap - t.gap) < 0.05);
    if (!seen) uniq.push_back(t);
  }
  return uniq;
}

System cleave(const System& bulk, const SlabOptions& o, SlabReport* rep) {
  SlabReport R;
  const PlaneSetup P = plane_setup(bulk, o.h, o.k, o.l);
  R.h = P.h, R.k = P.k, R.l = P.l, R.d = P.d;
  R.terminations = slab_terminations(bulk, o.h, o.k, o.l);
  if (o.termination < 0 || o.termination >= int(R.terminations.size()))
    throw CrystalError("termination " + std::to_string(o.termination + 1) + " of " + std::to_string(R.terminations.size()));
  const Termination& T = R.terminations[size_t(o.termination)];
  const int N = std::max(1, o.layers);
  // the frame: x along s1, z along the normal
  const Vec3 ex = P.s1 * (1 / norm(P.s1)), ez = P.n, ey = cross(ez, ex);
  double Ax = norm(P.s1), Bx = dot(P.s2, ex), By = dot(P.s2, ey);
  struct Site { int z; Vec3 p; std::string name; };
  std::vector<Site> sites;
  // whole molecules: each atom's unwrapped stacking coordinate within its molecule, and the molecule's centroid's
  std::vector<double> f3(bulk.atoms.size()), c3(bulk.atoms.size());
  for (size_t i = 0; i < bulk.atoms.size(); ++i) f3[i] = c3[i] = P.f[i][2];
  if (o.whole_molecules) {
    const auto ib = image_bonds(bulk);
    std::vector<std::vector<std::pair<uint32_t, Vec3>>> adj(bulk.atoms.size());
    for (const auto& b : ib) adj[b.i].push_back({b.j, b.d}), adj[b.j].push_back({b.i, b.d * -1.0});
    const auto inv = inverse_cols(P.s1, P.s2, P.s3);
    std::vector<char> seen(bulk.atoms.size(), 0);
    std::vector<Vec3> up(bulk.atoms.size());
    int endless = 0;
    for (size_t r = 0; r < bulk.atoms.size(); ++r) {
      if (seen[r]) continue;
      std::vector<uint32_t> mol{uint32_t(r)}, stack{uint32_t(r)};
      seen[r] = 1;
      up[r] = bulk.atoms[r].pos;
      bool loops = false;
      while (!stack.empty()) {
        const uint32_t a = stack.back();
        stack.pop_back();
        for (const auto& [b, d] : adj[a]) {
          const Vec3 pb = up[a] + d;
          if (!seen[b]) { seen[b] = 1, up[b] = pb, mol.push_back(b), stack.push_back(b); }
          else if (std::fabs(solve(inv, pb - up[b])[2]) > 0.5) loops = true;   // reaches its own image along the normal
        }
      }
      if (loops) { ++endless; continue; }   // endless along the normal: cut as the atoms are
      double c = 0;
      for (uint32_t a : mol) f3[a] = solve(inv, up[a] - bulk.cell.origin)[2], c += f3[a];
      c /= double(mol.size());
      // the molecule's centroid into [0, 1) and its atoms with it
      const double sh = -std::floor(c);
      for (uint32_t a : mol) f3[a] += sh, c3[a] = c + sh;
    }
    if (endless) R.notes.push_back(std::to_string(endless) + " molecule(s) endless along the surface normal: cut at the atomic planes");
    else R.notes.push_back("whole molecules: each placed by its centroid");
  }
  for (size_t i = 0; i < bulk.atoms.size(); ++i) {
    const Vec3 f{P.f[i][0], P.f[i][1], f3[i]};
    double g = c3[i] - T.cut / P.d;
    const double m = -std::floor(g);
    for (int L = 0; L < N; ++L) {
      const double gam = f[2] + m + L;
      const Vec3 r = P.s1 * f[0] + P.s2 * f[1] + P.s3 * gam;
      const double Z = gam * P.d - T.cut;
      sites.push_back({bulk.atoms[i].element, {dot(r, ex), dot(r, ey), Z}, bulk.atoms[i].name});
    }
  }
  auto wrap2 = [&](Vec3& p, double ax, double bx, double by) {   // into the mesh (ax, 0), (bx, by)
    double fb = p[1] / by;
    double fa = (p[0] - fb * bx) / ax;
    fa -= std::floor(fa), fb -= std::floor(fb);
    p[0] = fa * ax + fb * bx, p[1] = fb * by;
  };
  for (auto& s : sites) wrap2(s.p, Ax, Bx, By);
  // rectangular surface cell
  if (o.orthogonal && std::fabs(Bx) > 1e-6 * Ax) {
    int bi = 0, bj = 0;
    double be = 1e300;
    for (int j = 1; j <= 12 && be > 1e-9; ++j)
      for (int i = -12; i <= 12; ++i) {
        const double x = i * Ax + j * Bx, y = j * By;
        const double e = std::fabs(x) / y;
        if (e < be - 1e-12 || (std::fabs(e - be) < 1e-12 && j < bj)) { be = e; bi = i; bj = j; }
      }
    if (be > o.max_strain) {
      char b[160];
      std::snprintf(b, sizeof b, "no rectangular surface cell within %.1f %% shear (nearest %.2f %%); build without the orthogonal option", o.max_strain * 100, be * 100);
      throw CrystalError(b);
    }
    const double Bpx = bi * Ax + bj * Bx, Bpy = bj * By;
    std::vector<Site> out;
    for (const auto& s : sites)
      for (int q = 0; q < bj; ++q) {
        Site c = s;
        c.p = {s.p[0] + q * Bx, s.p[1] + q * By, s.p[2]};
        c.p[0] -= Bpx / Bpy * c.p[1];   // shear the lattice to 90°
        c.p[1] -= std::floor(c.p[1] / Bpy) * Bpy;
        c.p[0] -= std::floor(c.p[0] / Ax) * Ax;
        out.push_back(c);
      }
    sites = std::move(out);
    R.strain = Bpx / Bpy;
    if (be > 1e-9) {
      char b[120];
      std::snprintf(b, sizeof b, "surface cell sheared by %.2f %% to make it rectangular", be * 100);
      R.notes.push_back(b);
    }
    if (bj > 1) R.notes.push_back("rectangular surface cell = " + std::to_string(bj) + " surface meshes");
    Bx = 0, By = Bpy;
  }
  // surface supercell
  const int na = std::max(1, o.na), nb = std::max(1, o.nb);
  if (na > 1 || nb > 1) {
    std::vector<Site> out;
    for (int i = 0; i < na; ++i)
      for (int j = 0; j < nb; ++j)
        for (const auto& s : sites) out.push_back({s.z, {s.p[0] + i * Ax + j * Bx, s.p[1] + j * By, s.p[2]}, s.name});
    sites = std::move(out);
    Ax *= na, Bx *= nb, By *= nb;
  }
  double zmin = 1e300, zmax = -1e300;
  for (const auto& s : sites) zmin = std::min(zmin, s.p[2]), zmax = std::max(zmax, s.p[2]);
  const double vac = std::max(0.0, o.vacuum);
  const double Lz = vac > 0 ? zmax - zmin + vac : N * P.d;
  const double off = vac > 0 ? vac / 2 - zmin : 0;
  System s;
  s.source_format = "slab";
  s.title = (bulk.title.empty() ? std::string("crystal") : bulk.title) + " (" + std::to_string(P.h) + std::to_string(P.k) + std::to_string(P.l) + ") slab";
  s.cell.a = {Ax, 0, 0};
  s.cell.b = {Bx, By, 0};
  s.cell.c = {0, 0, Lz};
  std::map<int, int> type_of;
  auto add_atom = [&](int z, const Vec3& p, const std::string& name, int mol) {
    auto it = type_of.find(z);
    if (it == type_of.end()) {
      it = type_of.emplace(z, int(type_of.size()) + 1).first;
      TypeInfo ti;
      ti.type = it->second;
      ti.mass = element(z).mass;
      ti.label = element(z).symbol;
      s.types.push_back(ti);
    }
    Atom a;
    a.id = int64_t(s.atoms.size() + 1);
    a.mol = mol;
    a.element = z;
    a.type = it->second;
    a.name = name;
    a.pos = p;
    s.atoms.push_back(a);
    return uint32_t(s.atoms.size() - 1);
  };
  for (const auto& st : sites) add_atom(st.z, {st.p[0], st.p[1], st.p[2] + off}, st.name, 1);
  s.has_mol = true;
  s.unwrapped = true;
  s.bonds = crystal_bonds(s);
  R.thickness = zmax - zmin;
  R.a = Ax, R.b = std::hypot(Bx, By), R.gamma = std::atan2(By, Bx) / kDeg;

  if (o.passivate) {
    const double zmid = (zmin + zmax) / 2 + off;
    const auto [nh, noh] = passivate_surface(s, bulk, [&](const Vec3& p) { return Vec3{0, 0, p[2] > zmid ? 1.0 : -1.0}; });
    R.added_h = nh, R.added_oh = noh;
    if (R.added_h || R.added_oh)
      R.notes.push_back("passivated: " + std::to_string(R.added_oh) + " OH on cations, " + std::to_string(R.added_h) + " H on dangling bonds");
    // new atoms may stick out of the vacuum gap: keep the cell tall enough
    double top = -1e300, bot = 1e300;
    for (const auto& a : s.atoms) top = std::max(top, a.pos[2]), bot = std::min(bot, a.pos[2]);
    if (vac > 0) {
      const double shift = vac / 2 - bot;
      for (auto& a : s.atoms) a.pos[2] += shift;
      s.cell.c = {0, 0, top - bot + vac};
    }
  }
  std::vector<size_t> all(s.atoms.size());
  std::iota(all.begin(), all.end(), size_t(0));
  char b[200];
  std::snprintf(b, sizeof b, "(%d%d%d) slab · %s · %d layers of %.3f Å · %.2f × %.2f Å, γ %.1f° · %zu atoms (%s)", P.h, P.k, P.l, T.label.c_str(), N, P.d, R.a,
                R.b, R.gamma, s.atoms.size(), formula_of(s, all).c_str());
  R.notes.insert(R.notes.begin(), b);
  s.notes = R.notes;
  if (rep) *rep = R;
  return s;
}

}  // namespace caps
