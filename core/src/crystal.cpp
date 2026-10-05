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

// Adds O–H to under-coordinated oxygens (S, Se), M–OH to cations bonded to oxygen in the bulk and H to under-coordinated
// covalent network atoms of `s`, up to each element's coordination in `bulk`, pointing along outward_of(position).
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
  for (const auto& b : image_bonds(s)) {
    bonded[b.i].push_back(b.d * (1 / norm(b.d)));
    bonded[b.j].push_back(b.d * (-1 / norm(b.d)));
  }
  // metals keep their bare surfaces: only oxygen (and S, Se), cations of oxides and covalent network atoms
  auto passivable = [&](int z) {
    if (z == 8 || z == 16 || z == 34) return true;
    if (main_partner(z) == 8) return true;
    return z == 5 || z == 6 || z == 7 || z == 14 || z == 15 || z == 32;
  };
  auto pick = [&](const std::vector<Vec3>& taken, const Vec3& outward, double target) {
    Vec3 best = outward;
    double bs = -1e300;
    for (const auto& u : dirs) {
      double mina = 180;
      for (const auto& t : taken) mina = std::min(mina, std::acos(std::clamp(dot(u, t), -1.0, 1.0)) / kDeg);
      double sc = target > 0 ? -std::fabs(mina - target) * 2 : mina;
      if (mina < 95) sc -= 1000;
      sc += 40 * dot(outward, u);
      if (sc > bs) { bs = sc; best = u; }
    }
    return best;
  };
  const size_t n0 = s.atoms.size();
  for (size_t i = 0; i < n0; ++i) {
    const int z = s.atoms[i].element;
    const int missing = (want.count(z) ? want[z] : 0) - int(bonded[i].size());
    if (missing <= 0 || !passivable(z)) continue;
    const Vec3 outward = outward_of(s.atoms[i].pos);
    auto taken = bonded[i];
    for (int m = 0; m < missing; ++m) {
      const Vec3 u = pick(taken, outward, 0);
      taken.push_back(u);
      const int partner = main_partner(z);
      if (z == 8 || z == 16 || z == 34) {
        const uint32_t hdx = add_atom(s, 1, s.atoms[i].pos + u * (z == 8 ? 0.97 : 1.34), "H", 1);
        s.bonds.push_back({uint32_t(i), hdx, 1});
        ++added_h;
      } else if (partner == 8) {
        const Vec3 po = s.atoms[i].pos + u * bond_length(z, 8);
        const uint32_t odx = add_atom(s, 8, po, "O", 1);
        s.bonds.push_back({uint32_t(i), odx, 1});
        const Vec3 v = pick({u * -1.0}, outward, 115);
        const uint32_t hdx = add_atom(s, 1, po + v * 0.97, "H", 1);
        s.bonds.push_back({odx, hdx, 1});
        ++added_oh;
      } else {
        const uint32_t hdx = add_atom(s, 1, s.atoms[i].pos + u * (element(z).covalent + 0.31), "H", 1);
        s.bonds.push_back({uint32_t(i), hdx, 1});
        ++added_h;
      }
    }
  }
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

System parse_cif(const std::string& text, const std::string& name) {
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
  const int csym = col("_atom_site_type_symbol"), clab = col("_atom_site_label"), cocc = col("_atom_site_occupancy");
  int partial = 0;
  std::vector<Vec3> frac;
  std::vector<int> zs;
  std::vector<std::string> labels;
  for (const auto& r : sites) {
    const double occ = cocc >= 0 ? cif_number(r[size_t(cocc)], 1) : 1;
    if (occ < 0.5) { ++partial; continue; }
    int z = csym >= 0 ? element_of_site(r[size_t(csym)]) : 0;
    if (!z && clab >= 0) z = element_of_site(r[size_t(clab)]);
    if (!z) throw CrystalError("site " + (clab >= 0 ? r[size_t(clab)] : std::string("?")) + ": unknown element");
    const Vec3 f0{cif_number(r[size_t(cx_)]), cif_number(r[size_t(cy_)]), cif_number(r[size_t(cz_)])};
    for (const auto& op : ops) {
      Vec3 f;
      for (int i = 0; i < 3; ++i) f[i] = op.R[i][0] * f0[0] + op.R[i][1] * f0[1] + op.R[i][2] * f0[2] + op.t[i];
      f = wrap01(f);
      bool dup = false;
      for (size_t q = 0; q < frac.size() && !dup; ++q) {
        Vec3 df = f - frac[q];
        for (auto& x : df) x -= std::round(x);
        dup = norm(s.cell.to_cartesian(df) - s.cell.to_cartesian({0, 0, 0})) < 0.05;
      }
      if (dup) continue;
      frac.push_back(f);
      zs.push_back(z);
      labels.push_back(clab >= 0 ? r[size_t(clab)] : element(z).symbol);
    }
  }
  if (partial) s.notes.push_back(std::to_string(partial) + " sites with occupancy below 0.5 left out");
  std::map<int, int> type_of;
  for (size_t i = 0; i < frac.size(); ++i) {
    Atom a;
    a.id = int64_t(i + 1);
    a.mol = 1;
    a.element = zs[i];
    a.name = labels[i];
    a.pos = s.cell.a * frac[i][0] + s.cell.b * frac[i][1] + s.cell.c * frac[i][2];
    auto it = type_of.find(zs[i]);
    if (it == type_of.end()) {
      it = type_of.emplace(zs[i], int(type_of.size()) + 1).first;
      TypeInfo ti;
      ti.type = it->second;
      ti.mass = element(zs[i]).mass;
      ti.label = element(zs[i]).symbol;
      s.types.push_back(ti);
    }
    a.type = it->second;
    s.atoms.push_back(a);
  }
  s.has_mol = true;
  s.unwrapped = true;   // one periodic network: positions stay as placed (never "made whole")
  s.bonds = crystal_bonds(s);
  std::map<int, int> count;
  for (const auto& a : s.atoms) count[a.element]++;
  s.notes.insert(s.notes.begin(), std::to_string(ops.size()) + " symmetry operations · " + std::to_string(s.atoms.size()) + " atoms in the cell (" + hill(count) + ") · " +
                                      std::to_string(s.bonds.size()) + " bonds");
  return s;
}

System read_cif(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw CrystalError("cannot open " + path);
  std::stringstream ss;
  ss << in.rdbuf();
  return parse_cif(ss.str(), path);
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
