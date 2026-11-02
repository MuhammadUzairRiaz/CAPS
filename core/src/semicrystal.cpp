// CAPS semicrystalline cells (see semicrystal.hpp).
#include "caps/semicrystal.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <map>
#include <numeric>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>
#include <tuple>

#include "caps/crystal.hpp"
#include "caps/elements.hpp"
#include "caps/io.hpp"
#include "caps/json.hpp"
#include "caps/relax.hpp"
#include "caps/rng.hpp"
#include "caps/tags.hpp"
#include "caps/uff.hpp"

namespace caps {

namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr double kDeg = kPi / 180.0;
constexpr double kAmuPerA3 = 1.66053906660;      // g/cm³ for 1 g/mol per Å³

Vec3 unit(const Vec3& v) {
  const double n = norm(v);
  return n > 0 ? v * (1.0 / n) : Vec3{0, 0, 1};
}

// d bonded to c with |cd| = r, angle b-c-d = theta, torsion a-b-c-d = phi (radians; π = trans)
Vec3 place(const Vec3& a, const Vec3& b, const Vec3& c, double r, double theta, double phi) {
  const Vec3 bc = unit(c - b);
  Vec3 n = cross(b - a, bc);
  if (norm(n) < 1e-9) n = cross(bc, std::fabs(bc[0]) < 0.9 ? Vec3{1, 0, 0} : Vec3{0, 1, 0});
  n = unit(n);
  const Vec3 m = cross(n, bc);
  return c + bc * (-r * std::cos(theta)) + m * (r * std::sin(theta) * std::cos(phi)) + n * (r * std::sin(theta) * std::sin(phi));
}

double angle_deg(const Vec3& a, const Vec3& b, const Vec3& c) {
  const Vec3 u = unit(a - b), v = unit(c - b);
  return std::acos(std::clamp(dot(u, v), -1.0, 1.0)) / kDeg;
}

std::string formula_of(const System& s, const std::vector<uint32_t>& atoms) {
  std::map<int, int> c;
  for (uint32_t i : atoms) c[s.atoms[i].element]++;
  std::string f;
  for (int z : {6, 1}) if (c.count(z)) f += std::string(element(z).symbol) + (c[z] > 1 ? std::to_string(c[z]) : ""), c.erase(z);
  for (const auto& [z, k] : c) f += std::string(element(z).symbol) + (k > 1 ? std::to_string(k) : "");
  return f;
}

}  // namespace

// ---------------------------------------------------------------------------------------------------------------
// Unit cells

PolymerCrystalCell polymer_crystal_cell(const std::string& id, const std::string& data_dir) {
  PolymerCrystalCell pc;
  std::string low = id;
  for (auto& ch : low) ch = char(std::tolower((unsigned char)ch));
  // a file: a CIF or any structure with a cell
  if (std::filesystem::exists(id) && !std::filesystem::is_directory(id)) {
    pc.id = id;
    const bool cif = low.size() > 4 && low.substr(low.size() - 4) == ".cif";
    pc.cell = cif ? read_cif(id) : open_file(id).frame(0);
    if (!pc.cell.cell.valid()) throw std::invalid_argument(id + ": no unit cell");
    pc.name = pc.cell.title.empty() ? id : pc.cell.title;
    pc.source = id;
    return pc;
  }
  // an entry of CAPS's polymer crystal library (data/crystals/polymers/polymers.json): its id, or a short name
  if (low == "pe" || low == "polyethylene") low = "pe-orthorhombic";
  const std::filesystem::path dir = std::filesystem::path(data_dir) / "crystals" / "polymers";
  std::ifstream in(dir / "polymers.json");
  if (data_dir.empty() || !in)
    throw std::invalid_argument("'" + id + "' is not a file, and the polymer crystal library (data/crystals/polymers) was not found" +
                                (data_dir.empty() ? std::string() : " under " + data_dir));
  std::stringstream ss;
  ss << in.rdbuf();
  const Json lib = Json::parse(ss.str());
  std::string ids;
  for (const auto& e : lib["crystals"].items()) {
    ids += (ids.empty() ? "" : ", ") + e.text("id");
    if (e.text("id") != low) continue;
    pc.id = low;
    pc.name = e.text("name");
    pc.source = e.text("source") + (e.has("doi") ? ", doi:" + e.text("doi") : std::string()) + " (CAPS polymer crystal library, " + e.text("file") + ")";
    pc.cell = read_cif((dir / e.text("file")).string());
    pc.cell.title = pc.name;
    return pc;
  }
  throw std::invalid_argument("no polymer crystal '" + id + "' in the library (" + ids + ")");
}

std::vector<SemicrystallinePolymer> semicrystalline_polymers(const std::string& data_dir) {
  const std::filesystem::path dir = std::filesystem::path(data_dir) / "crystals" / "polymers";
  auto read = [](const std::filesystem::path& f) {
    std::ifstream in(f);
    if (!in) throw std::invalid_argument("cannot read " + f.string());
    std::stringstream ss;
    ss << in.rdbuf();
    return Json::parse(ss.str());
  };
  const Json lib = read(dir / "polymers.json"), def = read(dir / "semicrystalline.json");
  std::map<std::string, std::string> names;
  for (const auto& e : lib["crystals"].items()) names[e.text("id")] = e.text("name");
  std::vector<SemicrystallinePolymer> out;
  for (const auto& e : def["polymers"].items()) {
    SemicrystallinePolymer p;
    p.id = e.text("id"), p.name = names.count(p.id) ? names[p.id] : p.id;
    p.forcefield = e.text("forcefield"), p.density_source = e.text("density_source");
    if (e.has("repeats")) p.nx = int(e["repeats"][0].number()), p.ny = int(e["repeats"][1].number());
    p.crystal_nm = e.num("crystal_nm", p.crystal_nm), p.amorphous_nm = e.num("amorphous_nm", p.amorphous_nm);
    p.loops = e.num("loops", p.loops), p.bridges = e.num("bridges", p.bridges), p.tails = e.num("tails", p.tails);
    p.density = e.num("density", p.density);
    out.push_back(p);
  }
  return out;
}

bool semicrystalline_polymer(const std::string& id, const std::string& data_dir, SemicrystallinePolymer* out) {
  std::string low = id;
  for (auto& ch : low) ch = char(std::tolower((unsigned char)ch));
  if (low == "pe" || low == "polyethylene") low = "pe-orthorhombic";
  if (data_dir.empty() || !std::filesystem::exists(std::filesystem::path(data_dir) / "crystals" / "polymers" / "semicrystalline.json")) return false;
  for (const auto& p : semicrystalline_polymers(data_dir))
    if (p.id == low) {
      if (out) *out = p;
      return true;
    }
  return false;
}

std::vector<CrystalStem> find_stems(const System& s, double tol) {
  const size_t n = s.atoms.size();
  if (!s.cell.valid()) throw std::invalid_argument("find_stems: no unit cell");
  // bonds through the faces: (j, image shift) for each atom, by covalent radii + tolerance
  struct Link { uint32_t j; std::array<int, 3> sh; };
  std::vector<std::vector<Link>> nb(n);
  for (uint32_t i = 0; i < n; ++i)
    for (uint32_t j = 0; j < n; ++j)
      for (int a = -1; a <= 1; ++a)
        for (int b = -1; b <= 1; ++b)
          for (int c = -1; c <= 1; ++c) {
            if (i == j && a == 0 && b == 0 && c == 0) continue;
            const Vec3 d = s.atoms[j].pos + s.cell.a * a + s.cell.b * b + s.cell.c * c - s.atoms[i].pos;
            const double lim = element(s.atoms[i].element).covalent + element(s.atoms[j].element).covalent + tol;
            if (norm(d) < lim && norm(d) > 0.4) nb[i].push_back({j, {a, b, c}});
          }
  // components of the bond graph with each atom's image within its chain (positions unwrapped along the chain)
  std::vector<int> comp(n, -1);
  std::vector<std::array<int, 3>> img(n);
  std::vector<CrystalStem> out;
  for (uint32_t r = 0; r < n; ++r) {
    if (comp[r] >= 0) continue;
    const int id = int(out.size());
    std::vector<uint32_t> members{r}, stack{r};
    comp[r] = id;
    img[r] = {0, 0, 0};
    std::array<int, 3> trans{0, 0, 0};
    bool repeats = false;
    while (!stack.empty()) {
      const uint32_t u = stack.back();
      stack.pop_back();
      for (const auto& l : nb[u]) {
        const std::array<int, 3> w{img[u][0] + l.sh[0], img[u][1] + l.sh[1], img[u][2] + l.sh[2]};
        if (comp[l.j] < 0) {
          comp[l.j] = id;
          img[l.j] = w;
          members.push_back(l.j);
          stack.push_back(l.j);
        } else if (comp[l.j] == id && w != img[l.j]) {   // the chain meets itself one translation away
          const std::array<int, 3> t{w[0] - img[l.j][0], w[1] - img[l.j][1], w[2] - img[l.j][2]};
          const int nt = std::abs(t[0]) + std::abs(t[1]) + std::abs(t[2]), nr = std::abs(trans[0]) + std::abs(trans[1]) + std::abs(trans[2]);
          if (!repeats || nt < nr) trans = t, repeats = true;
        }
      }
    }
    CrystalStem st;
    st.atoms = members;
    st.formula = formula_of(s, members);
    // a chain repeats with one lattice translation: along a, b or c
    int ax = -1;
    if (repeats && std::abs(trans[0]) + std::abs(trans[1]) + std::abs(trans[2]) == 1) ax = trans[0] ? 0 : trans[1] ? 1 : 2;
    if (ax < 0) {
      out.push_back(st);   // a molecule, or a chain along no single cell edge: atoms only, no backbone
      continue;
    }
    st.lattice_axis = ax;
    const Vec3 edge = ax == 0 ? s.cell.a : ax == 1 ? s.cell.b : s.cell.c;
    st.axis = unit(edge);
    // backbone: heavy atoms that lie on the path from one cell's chain to the next (their removal splits the chain);
    // found as the heavy atoms with two or more heavy neighbours whose shortest path to their own +c image passes them
    std::vector<uint32_t> heavy;   // heavy atoms with two heavy neighbours or more (not a methyl, F or carbonyl O)
    for (uint32_t i : members) {
      if (s.atoms[i].element <= 1) continue;
      int hn = 0;
      for (const auto& l : nb[i]) hn += s.atoms[l.j].element > 1;
      if (hn >= 2) heavy.push_back(i);
    }
    // BFS over (atom, image) restricted to heavy atoms from the lowest heavy atom to its +c image
    auto unwrapped = [&](uint32_t i, const std::array<int, 3>& im) {
      return s.atoms[i].pos + s.cell.a * im[0] + s.cell.b * im[1] + s.cell.c * im[2];
    };
    uint32_t start = heavy.empty() ? members[0] : heavy[0];
    for (uint32_t i : heavy) if (dot(unwrapped(i, img[i]), st.axis) < dot(unwrapped(start, img[start]), st.axis)) start = i;
    using Key = std::pair<uint32_t, std::array<int, 3>>;
    std::map<Key, Key> prev;
    std::vector<Key> q{{start, img[start]}};
    prev[q[0]] = q[0];
    std::array<int, 3> gim = img[start];
    gim[size_t(ax)] += 1;
    const Key goal{start, gim};
    bool found = false;
    for (size_t h = 0; h < q.size() && !found; ++h) {
      const Key cur = q[h];
      for (const auto& l : nb[cur.first]) {
        if (std::find(heavy.begin(), heavy.end(), l.j) == heavy.end()) continue;
        const Key nx{l.j, {cur.second[0] + l.sh[0], cur.second[1] + l.sh[1], cur.second[2] + l.sh[2]}};
        if (prev.count(nx)) continue;
        if (std::abs(nx.second[size_t(ax)] - img[start][size_t(ax)]) > 2) continue;
        prev[nx] = cur;
        if (nx == goal) { found = true; break; }
        q.push_back(nx);
      }
    }
    if (found) {
      std::vector<Key> path;
      for (Key k = prev[goal]; !(k == q[0]); k = prev[k]) path.push_back(k);
      path.push_back(q[0]);
      std::reverse(path.begin(), path.end());
      for (const auto& k : path) st.backbone.push_back(k.first), st.image.push_back(k.second);
      // the backbone in order along +c, starting at the lowest atom (already: start is the lowest)
      // helical: the backbone's positions across the axis spread in two directions (a planar zig-zag in one)
      if (st.backbone.size() >= 3) {
        std::vector<Vec3> q;
        Vec3 m{0, 0, 0};
        for (const auto& k : path) {
          Vec3 p = unwrapped(k.first, k.second);
          p = p - st.axis * dot(p, st.axis);
          q.push_back(p);
          m = m + p;
        }
        m = m * (1.0 / double(q.size()));
        double c[3][3] = {{0}};
        for (const auto& p : q)
          for (int i2 = 0; i2 < 3; ++i2)
            for (int j2 = 0; j2 < 3; ++j2) c[i2][j2] += (p[i2] - m[i2]) * (p[j2] - m[j2]) / double(q.size());
        // the two in-plane eigenvalues: the smaller one is the trace minus the largest; the largest by power iteration
        Vec3 v{1, 0.3, 0.1};
        for (int it = 0; it < 60; ++it) {
          Vec3 w{c[0][0] * v[0] + c[0][1] * v[1] + c[0][2] * v[2], c[1][0] * v[0] + c[1][1] * v[1] + c[1][2] * v[2], c[2][0] * v[0] + c[2][1] * v[1] + c[2][2] * v[2]};
          if (norm(w) < 1e-12) break;
          v = unit(w);
        }
        const double l1 = dot(v, Vec3{c[0][0] * v[0] + c[0][1] * v[1] + c[0][2] * v[2], c[1][0] * v[0] + c[1][1] * v[1] + c[1][2] * v[2], c[2][0] * v[0] + c[2][1] * v[1] + c[2][2] * v[2]});
        const double l2 = c[0][0] + c[1][1] + c[2][2] - l1;
        st.helical = l2 > 0.04;   // more than 0.2 Å RMS across the zig-zag plane
      }
    }
    out.push_back(st);
  }
  return out;
}

// ---------------------------------------------------------------------------------------------------------------
// Builder

namespace {

// The lamellar cell: A along x and B in the xy plane (the crystal's a and b repeated), C along z (the stack normal).
struct Box {
  double ax, bx, by, lz;
  Vec3 mi(Vec3 d) const {
    d[2] -= lz * std::round(d[2] / lz);
    const double nb = std::round(d[1] / by);
    d[0] -= bx * nb, d[1] -= by * nb;
    d[0] -= ax * std::round(d[0] / ax);
    return d;
  }
  double lateral(Vec3 d) const {
    d[2] = 0;
    return norm(mi(d));
  }
  void frac(const Vec3& p, double f[3]) const {
    f[1] = p[1] / by;
    f[0] = (p[0] - bx * f[1]) / ax;
    f[2] = p[2] / lz;
  }
  double width(int k) const { return k == 0 ? ax * by / std::hypot(bx, by) : k == 1 ? by : lz; }   // between opposite faces
};

// Spatial hash of the heavy atoms placed so far (positions kept unwrapped by the caller), each with the backbone unit
// it belongs to (for the exclusions along a chain) and a dead flag (atoms taken back).
struct Hash {
  Box box;
  int n[3];
  std::vector<std::vector<uint32_t>> bins;
  std::vector<Vec3> pos;
  std::vector<int32_t> owner;
  std::vector<char> dead;
  Hash(const Box& b, double cell) : box(b) {
    for (int k = 0; k < 3; ++k) n[k] = std::max(1, int(b.width(k) / cell));
    bins.resize(size_t(n[0]) * n[1] * n[2]);
  }
  int bin(double f, int k) const {
    f -= std::floor(f);
    return std::min(n[k] - 1, int(f * n[k]));
  }
  uint32_t add(const Vec3& p, int32_t own) {
    const uint32_t id = uint32_t(pos.size());
    pos.push_back(p), owner.push_back(own), dead.push_back(0);
    double f[3];
    box.frac(p, f);
    bins[(size_t(bin(f[0], 0)) * n[1] + bin(f[1], 1)) * n[2] + bin(f[2], 2)].push_back(id);
    return id;
  }
  // f(id, distance) for every live point within the bins around p (the caller checks the distance)
  template <class F>
  void near(const Vec3& p, F f) const {
    double fr[3];
    box.frac(p, fr);
    const int b[3] = {bin(fr[0], 0), bin(fr[1], 1), bin(fr[2], 2)};
    std::vector<int> off[3];
    for (int k = 0; k < 3; ++k) {
      if (n[k] < 3) for (int q = 0; q < n[k]; ++q) off[k].push_back(q);
      else for (int d = -1; d <= 1; ++d) off[k].push_back(((b[k] + d) % n[k] + n[k]) % n[k]);
    }
    for (int x : off[0])
      for (int y : off[1])
        for (int z : off[2])
          for (uint32_t id : bins[(size_t(x) * n[1] + y) * n[2] + z])
            if (!dead[id]) f(id, norm(box.mi(pos[id] - p)));
  }
};

// torsion a-b-c-d with place()'s convention (π trans)
double dihedral(const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& d) {
  const Vec3 bc = unit(c - b);
  const Vec3 n = unit(cross(b - a, bc)), m = cross(n, bc);
  const Vec3 w = d - c;
  return std::atan2(dot(w, n), dot(w, m));
}

// The local frame of a backbone atom from its neighbours along the chain: e1 the bisector towards them, e2 the normal of
// their plane, e3 = e1 × e2. Walking the chain the other way (prev and next swapped) turns (x, y, z) into (x, −y, −z).
struct Frame {
  Vec3 o, e1, e2, e3;
};
Frame frame_of(const Vec3& prev, const Vec3& at, const Vec3& next) {
  const Vec3 u = unit(prev - at), v = unit(next - at);
  Frame f;
  f.o = at, f.e1 = unit(u + v), f.e2 = unit(cross(u, v)), f.e3 = cross(f.e1, f.e2);
  return f;
}
Vec3 to_frame(const Frame& f, const Vec3& p) {
  const Vec3 d = p - f.o;
  return {dot(d, f.e1), dot(d, f.e2), dot(d, f.e3)};
}
Vec3 from_frame(const Frame& f, const Vec3& q) { return f.o + f.e1 * q[0] + f.e2 * q[1] + f.e3 * q[2]; }
Vec3 walked_back(const Vec3& q) { return {q[0], -q[1], -q[2]}; }

int wrapi(int i, int p) { return ((i % p) + p) % p; }

// classes of a backbone bond for the torsions grown about it
enum TorsionClass { kSp3 = 0, kSp2Sp3 = 1, kConj = 2, kFixed = 3 };

struct Side {
  int el = 0;
  Vec3 q{0, 0, 0};   // frame coordinates (template) or position (a unit of a cell's stem)
  int parent = -1;   // the side atom it hangs from, −1 the backbone atom
};
// a bond the side groups make beyond their own tree: side a of the unit at sigma to side b (−1: the backbone atom) of the
// unit delta further along the chain (in the template's direction)
struct ExtraBond {
  int sigma, a, delta, b;
  bool operator<(const ExtraBond& o) const { return std::tie(sigma, a, delta, b) < std::tie(o.sigma, o.a, o.delta, o.b); }
};

// one backbone atom of a stem in the unit cell, with the atoms that hang from it (reached without passing another
// backbone atom), positions unwrapped along the chain
struct CellUnit {
  int el = 0, nbrs = 0;
  bool carbonyl = false;     // a carbon with a double-bonded O on it
  Vec3 pos{0, 0, 0};
  std::vector<Side> side;    // q = position (Å)
};
struct CellStem {
  std::vector<CellUnit> u;   // one period, upwards along c
  Vec3 period{0, 0, 0};
  std::vector<int> cls;      // torsion class of the bond from unit k to unit k + 1
  std::vector<ExtraBond> extra;   // sigma = the unit index here, a and b in this unit's raw side order
  // ways to read it against the template: sigma(k) = m + d k (mod P), and for each unit the template's side order
  struct Labeling {
    int m = 0, d = 1;
    double rms = 0;
    std::vector<std::vector<int>> perm;   // perm[k][t]: the raw side index matched to the template's side t
  };
  std::vector<Labeling> labels;
  Vec3 side_q(int k, int s) const {   // frame coordinates of side s of unit k, walked upwards
    const int n = int(u.size());
    const Vec3 prev = k > 0 ? u[size_t(k - 1)].pos : u.back().pos - period;
    const Vec3 next = k + 1 < n ? u[size_t(k + 1)].pos : u[0].pos + period;
    return to_frame(frame_of(prev, u[size_t(k)].pos, next), u[size_t(k)].side[size_t(s)].q);
  }
};

// The repeat unit the amorphous segments are grown from: the shortest period of the crystal's chain (one CH2 for PE,
// CH2-CF2 for PVDF, the whole repeat for nylon 6), with each backbone atom's bond, angle and torsion and its side atoms
// in the backbone atom's frame, read off the crystal (so the stereochemistry is the crystal's).
struct Template {
  int P = 1;
  std::vector<int> el, cls;
  std::vector<double> r, th, phi;   // bond sigma→sigma+1, angle at sigma, torsion sigma−1, sigma, sigma+1, sigma+2
  std::vector<std::vector<Side>> side;
  std::vector<ExtraBond> extra;
  double unit_mass = 0, reach = 0, r_mean = 0;
  bool polar = false;
  int next(int s, int d) const { return wrapi(s + d, P); }
  double bond(int s, int d) const { return d > 0 ? r[size_t(s)] : r[size_t(wrapi(s - 1, P))]; }
  // the torsion that places the atom after s (walking d): about the bond from the atom before s to s
  int tcls(int s, int d) const { return d > 0 ? cls[size_t(wrapi(s - 1, P))] : cls[size_t(s)]; }
  double tref(int s, int d) const { return d > 0 ? phi[size_t(wrapi(s - 1, P))] : phi[size_t(s)]; }
  // side positions of the unit at sigma s, its frame from prev and next as walked in direction d
  std::vector<Vec3> sides(int s, int d, const Vec3& prev, const Vec3& at, const Vec3& next) const {
    const Frame f = frame_of(prev, at, next);
    std::vector<Vec3> out;
    for (const Side& x : side[size_t(s)]) out.push_back(from_frame(f, d > 0 ? x.q : walked_back(x.q)));
    return out;
  }
};

// matches side groups by element and frame coordinates: the permutation (b's index for each of a's) and the largest
// deviation, or −1 when the groups differ
double match_sides(const std::vector<Side>& a, const std::vector<Vec3>& aq, const std::vector<Side>& b, const std::vector<Vec3>& bq,
                   std::vector<int>* perm) {
  if (a.size() != b.size()) return -1;
  std::vector<int> pm(a.size(), -1);
  std::vector<char> used(b.size(), 0);
  double worst = 0;
  for (size_t i = 0; i < a.size(); ++i) {
    int best = -1;
    double bd = 1e300;
    for (size_t j = 0; j < b.size(); ++j)
      if (!used[j] && b[j].el == a[i].el && norm(bq[j] - aq[i]) < bd) bd = norm(bq[j] - aq[i]), best = int(j);
    if (best < 0) return -1;
    used[size_t(best)] = 1, pm[i] = best;
    worst = std::max(worst, bd - (a[i].el == 1 ? 0.25 : 0.0));   // hydrogens placed by rule in some cells: a little more slack
  }
  if (perm) *perm = pm;
  return worst;
}

struct Link {
  uint32_t j;
  std::array<int, 3> sh;
};

std::vector<std::vector<Link>> cell_links(const System& s, double tol) {
  const size_t n = s.atoms.size();
  std::vector<std::vector<Link>> nb(n);
  for (uint32_t i = 0; i < n; ++i)
    for (uint32_t j = 0; j < n; ++j)
      for (int a = -1; a <= 1; ++a)
        for (int b = -1; b <= 1; ++b)
          for (int c = -1; c <= 1; ++c) {
            if (i == j && a == 0 && b == 0 && c == 0) continue;
            const Vec3 d = s.atoms[j].pos + s.cell.a * a + s.cell.b * b + s.cell.c * c - s.atoms[i].pos;
            const double lim = element(s.atoms[i].element).covalent + element(s.atoms[j].element).covalent + tol;
            if (norm(d) < lim && norm(d) > 0.4) nb[i].push_back({j, {a, b, c}});
          }
  return nb;
}

using Key = std::pair<uint32_t, std::array<int, 3>>;
std::array<int, 3> add3(const std::array<int, 3>& a, const std::array<int, 3>& b) { return {a[0] + b[0], a[1] + b[1], a[2] + b[2]}; }

// a stem of the unit cell as backbone units with their side atoms, and its bonds' torsion classes
CellStem cell_stem(const System& uc, const CrystalStem& st, const std::vector<std::vector<Link>>& links) {
  CellStem cs;
  const size_t nb = st.backbone.size();
  cs.period = uc.cell.c;
  auto at = [&](uint32_t i, const std::array<int, 3>& im) { return uc.atoms[i].pos + uc.cell.a * im[0] + uc.cell.b * im[1] + uc.cell.c * im[2]; };
  std::map<uint32_t, size_t> bbi;
  for (size_t k = 0; k < nb; ++k) bbi[st.backbone[k]] = k;
  cs.u.resize(nb);
  // side atoms: breadth first from every backbone atom at once, never through a backbone atom
  std::map<uint32_t, std::pair<size_t, int>> owner;   // atom -> unit, index in its side list
  std::map<uint32_t, std::array<int, 3>> img;
  std::vector<std::pair<uint32_t, std::array<int, 3>>> q;
  for (size_t k = 0; k < nb; ++k) {
    CellUnit& u = cs.u[k];
    u.el = uc.atoms[st.backbone[k]].element;
    u.pos = at(st.backbone[k], st.image[k]);
    u.nbrs = int(links[st.backbone[k]].size());
    q.push_back({st.backbone[k], st.image[k]});
  }
  for (size_t h = 0; h < q.size(); ++h) {
    const auto [i, im] = q[h];
    const bool bb = bbi.count(i) > 0;
    const size_t k = bb ? bbi[i] : owner[i].first;
    for (const Link& l : links[i]) {
      if (bbi.count(l.j) || owner.count(l.j)) continue;
      const std::array<int, 3> jm = add3(im, l.sh);
      CellUnit& u = cs.u[k];
      owner[l.j] = {k, int(u.side.size())};
      img[l.j] = jm;
      u.side.push_back({uc.atoms[l.j].element, at(l.j, jm), bb ? -1 : owner[i].second});
      if (bb && u.el == 6 && uc.atoms[l.j].element == 8 && links[l.j].size() == 1) u.carbonyl = true;   // C=O
      q.push_back({l.j, jm});
    }
  }
  // the global unit index (along the chain, unit k of the period shifted by t periods) of a backbone key
  auto unit_of = [&](uint32_t i, const std::array<int, 3>& im) {
    const size_t k = bbi[i];
    const std::array<int, 3> d{im[0] - st.image[k][0], im[1] - st.image[k][1], im[2] - st.image[k][2]};
    return int(k) + int(nb) * d[2];
  };
  std::set<ExtraBond> extra;
  for (const auto& [i, ow] : owner) {
    const int ki = int(ow.first);
    for (const Link& l : links[i]) {
      const std::array<int, 3> jm = add3(img[i], l.sh);
      if (bbi.count(l.j)) {
        const int g = unit_of(l.j, jm);
        if (g != ki || cs.u[size_t(ki)].side[size_t(ow.second)].parent != -1) extra.insert({ki, ow.second, g - ki, -1});
        continue;
      }
      const auto& oj = owner[l.j];
      const int g = int(oj.first) + int(nb) * (jm[2] - img[l.j][2]);
      const int delta = g - ki;
      const auto& si = cs.u[size_t(ki)].side[size_t(ow.second)];
      const auto& sj = cs.u[oj.first].side[size_t(oj.second)];
      if (delta == 0 && (si.parent == oj.second || sj.parent == ow.second)) continue;   // the side group's own tree
      if (delta > 0 || (delta == 0 && ow.second < oj.second)) extra.insert({ki, ow.second, delta, oj.second});
    }
  }
  cs.extra.assign(extra.begin(), extra.end());
  // torsion classes: a bond in a ring, an amide or ester C(=O)-N/O, or a C=C double bond is held; between two trigonal
  // atoms (aryl-carbonyl) two states; next to one trigonal atom three; sp3-sp3 the trans-gauche set
  auto trigonal = [&](int k) {
    const CellUnit& u = cs.u[size_t(wrapi(k, int(nb)))];
    const CellUnit& a = cs.u[size_t(wrapi(k - 1, int(nb)))];
    const CellUnit& b = cs.u[size_t(wrapi(k + 1, int(nb)))];
    if (u.el == 6) return u.nbrs == 3;
    if (u.el == 7) return u.nbrs == 3 && (a.carbonyl || b.carbonyl);
    if (u.el == 8) return a.carbonyl || b.carbonyl;
    return false;
  };
  for (size_t k = 0; k < nb; ++k) {
    const size_t k1 = (k + 1) % nb;
    const Key from{st.backbone[k], st.image[k]};
    std::array<int, 3> im1 = st.image[k1];
    if (k1 == 0) im1[2] += 1;
    const Key to{st.backbone[k1], im1};
    // ring: k+1 reached from k within eight bonds without the direct bond
    bool ring = false;
    {
      std::map<Key, int> dist{{from, 0}};
      std::vector<Key> bq{from};
      for (size_t h = 0; h < bq.size() && !ring; ++h) {
        const Key cur = bq[h];
        if (dist[cur] >= 8) continue;
        for (const Link& l : links[cur.first]) {
          const Key nx{l.j, add3(cur.second, l.sh)};
          if (cur == from && nx == to) continue;
          if (dist.count(nx)) continue;
          if (nx == to) { ring = true; break; }
          dist[nx] = dist[cur] + 1;
          bq.push_back(nx);
        }
      }
    }
    const CellUnit &u = cs.u[k], &v = cs.u[k1];
    const double len = norm(at(to.first, to.second) - u.pos);
    int c = kSp3;
    const bool tu = trigonal(int(k)), tv = trigonal(int(k1));
    if (ring) c = kFixed;
    else if ((u.carbonyl && (v.el == 7 || v.el == 8)) || (v.carbonyl && (u.el == 7 || u.el == 8))) c = kFixed;
    else if (u.el == 6 && v.el == 6 && u.nbrs == 3 && v.nbrs == 3 && len < 1.41) c = kFixed;
    else if (tu && tv) c = kConj;
    else if (tu || tv) c = kSp2Sp3;
    cs.cls.push_back(c);
  }
  return cs;
}

// the cell with its chains along c (a and b relabelled cyclically when they run along a or b), a along x, b in xy
System chains_along_c(const System& uc0, int axis) {
  const Vec3 e[3] = {uc0.cell.a, uc0.cell.b, uc0.cell.c};
  const int ia = (axis + 1) % 3, ib = (axis + 2) % 3, ic = axis;
  const double la = norm(e[ia]), lb = norm(e[ib]), lc = norm(e[ic]);
  const double ca = dot(e[ib], e[ic]) / (lb * lc), cb = dot(e[ia], e[ic]) / (la * lc), cg = dot(e[ia], e[ib]) / (la * lb);
  const double sg = std::sqrt(std::max(0.0, 1 - cg * cg));
  System s = uc0;
  s.cell.origin = {0, 0, 0};
  s.cell.a = {la, 0, 0};
  s.cell.b = {lb * cg, lb * sg, 0};
  const double cx = lc * cb, cy = lc * (ca - cb * cg) / sg;
  s.cell.c = {cx, cy, std::sqrt(std::max(0.0, lc * lc - cx * cx - cy * cy))};
  for (size_t i = 0; i < s.atoms.size(); ++i) {
    const Vec3 f = uc0.cell.to_fractional(uc0.atoms[i].pos - uc0.cell.origin);
    Vec3 g{f[ia], f[ib], f[ic]};
    for (int k = 0; k < 3; ++k) g[k] -= std::floor(g[k]);
    s.atoms[i].pos = s.cell.to_cartesian(g);
  }
  return s;
}

struct Seg {
  int type = 0;               // 0 loop, 1 bridge, 2 tail
  int a = -1, b = -1;         // stem ends (2 stem + 1 top, 2 stem bottom); b = −1 for a tail
  int n = 0;                  // backbone atoms to grow
  std::vector<int> beads;     // unit ids, in growth order
  int d = 1, q = 0;           // growth direction along the template; length ≡ q (mod P) for a loop or bridge
  Seg() = default;
  Seg(int t, int e1, int e2) : type(t), a(e1), b(e2) {}
};

// a backbone atom of the cell being built, with its side atoms (template order once its sigma is known)
struct URec {
  int el = 0, sigma = -1;
  bool crystal = false;
  Vec3 pos{0, 0, 0};
  std::vector<Vec3> side;
  std::vector<int> side_el;        // raw elements of a crystal unit before relabelling (template's otherwise)
  std::vector<uint32_t> hs;        // hash ids of its heavy side atoms
  int cs = -1, kc = -1, g = 0;     // crystal units: cell stem, unit index in the period, global index along the stem
  bool cap = false;
  Vec3 cap_pos{0, 0, 0};
};

double bond_to_h(int el) { return el == 7 ? 1.01 : el == 8 ? 0.96 : 1.09; }

}  // namespace

System build_semicrystalline(const SemicrystallineOptions& o, SemicrystallineReport* report) {
  auto cancelled = [&](const std::string& what, double f) { return o.progress && !o.progress(what, f); };
  SemicrystallineReport R;
  const PolymerCrystalCell pc = polymer_crystal_cell(o.polymer, o.data_dir);
  // --- the unit cell with its chains along c
  System uc = pc.cell;
  std::vector<CrystalStem> stems_uc = find_stems(uc);
  for (const auto& st : stems_uc) {
    if (st.backbone.empty()) throw std::invalid_argument(pc.name + ": a molecule, or a chain along no cell edge (" + st.formula + ")");
    if (st.lattice_axis != stems_uc[0].lattice_axis) throw std::invalid_argument(pc.name + ": chains along different cell edges");
  }
  if (stems_uc.empty()) throw std::invalid_argument(pc.name + ": no chains in the cell");
  if (stems_uc[0].lattice_axis != 2) R.notes.push_back(std::string("the cell's edges relabelled so that the chains run along c (they run along ") + char('a' + stems_uc[0].lattice_axis) + " in the file)");
  uc = chains_along_c(uc, stems_uc[0].lattice_axis);
  stems_uc = find_stems(uc);
  const auto links = cell_links(uc, 0.45);
  std::vector<CellStem> cst;
  for (const auto& st : stems_uc) {
    if (st.lattice_axis != 2 || st.backbone.empty()) throw std::invalid_argument(pc.name + ": a chain that does not repeat along c after relabelling");
    cst.push_back(cell_stem(uc, st, links));
  }
  const int nbb = int(cst[0].u.size());
  // --- the template: the shortest period of stem 0 along which the units repeat
  Template T;
  {
    const CellStem& s0 = cst[0];
    auto qs = [&](const CellStem& s, int k) {
      std::vector<Vec3> v;
      for (size_t i = 0; i < s.u[size_t(k)].side.size(); ++i) v.push_back(s.side_q(k, int(i)));
      return v;
    };
    auto geo = [&](const CellStem& s, int k, double& r, double& th, double& phi) {
      const int n = int(s.u.size());
      auto p = [&](int i) { const int w = int(std::floor(double(i) / n)); return s.u[size_t(wrapi(i, n))].pos + s.period * double(w); };
      r = norm(p(k + 1) - p(k));
      th = angle_deg(p(k - 1), p(k), p(k + 1)) * kDeg;
      phi = dihedral(p(k - 1), p(k), p(k + 1), p(k + 2));
    };
    int P = nbb;
    for (int p = 1; p < nbb; ++p) {
      if (nbb % p) continue;
      bool ok = true;
      for (int k = 0; k < nbb && ok; ++k) {
        const int k2 = (k + p) % nbb;
        const CellUnit &u = s0.u[size_t(k)], &v = s0.u[size_t(k2)];
        if (u.el != v.el || s0.cls[size_t(k)] != s0.cls[size_t(k2)]) { ok = false; break; }
        const double dv = match_sides(u.side, qs(s0, k), v.side, qs(s0, k2), nullptr);
        double r1, t1, f1, r2, t2, f2;
        geo(s0, k, r1, t1, f1), geo(s0, k2, r2, t2, f2);
        if (dv < 0 || dv > 0.3 || std::fabs(r1 - r2) > 0.06 || std::fabs(t1 - t2) > 4 * kDeg) ok = false;
      }
      if (ok) { P = p; break; }
    }
    T.P = P;
    for (int k = 0; k < P; ++k) {
      double r, th, phi;
      geo(s0, k, r, th, phi);
      T.el.push_back(s0.u[size_t(k)].el), T.cls.push_back(s0.cls[size_t(k)]);
      T.r.push_back(r), T.th.push_back(th), T.phi.push_back(phi);
      std::vector<Side> sd = s0.u[size_t(k)].side;
      const auto q = qs(s0, k);
      for (size_t i = 0; i < sd.size(); ++i) sd[i].q = q[i];
      T.side.push_back(sd);
      T.r_mean += r / P;
      T.reach += r * std::sin(th / 2) / P;
    }
    double mass = 0;
    for (const auto& u : s0.u) {
      mass += element(u.el).mass;
      for (const auto& x : u.side) mass += element(x.el).mass;
    }
    T.unit_mass = mass / nbb;
  }
  {   // the template's bonds and angles regrow the crystal's chain from three atoms, both ways along it (a check of the
      // bookkeeping; the period was found within 0.06 Å and 4° per unit, so the regrown chain may drift a little)
    const CellStem& s0 = cst[0];
    auto p = [&](int i) { return s0.u[size_t(wrapi(i, nbb))].pos + s0.period * std::floor(double(i) / nbb); };
    double worst = 0;
    for (int d : {1, -1}) {
      std::vector<Vec3> x{p(-2 * d), p(-d), p(0)};
      for (int i = 1; i <= 2 * nbb; ++i) {
        const int sc = wrapi(d * (i - 1), T.P);   // the head's sigma (stem 0 reads the template from m = 0)
        const size_t L = x.size();
        // the crystal's own torsion here (a helix's monomers may differ a little in it); bonds and angles the template's
        const double phi = dihedral(p(d * (i - 3)), p(d * (i - 2)), p(d * (i - 1)), p(d * i));
        x.push_back(place(x[L - 3], x[L - 2], x[L - 1], T.bond(sc, d), T.th[size_t(sc)], phi));
        worst = std::max(worst, norm(x.back() - p(d * i)));
      }
    }
    if (worst > 0.5) throw std::logic_error(pc.name + ": the repeat unit does not regrow its own chain (off by " + std::to_string(worst) + " Å)");
  }
  // --- every stem of the cell read against the template: its offsets and directions (both for a chain with no
  // direction, as PE or PVDF; one for nylon 6, iPP or PLLA, whose stems point up or down)
  for (size_t j = 0; j < cst.size(); ++j) {
    CellStem& s = cst[j];
    if (int(s.u.size()) != nbb) throw std::invalid_argument(pc.name + ": stems of different lengths per cell");
    for (int d : {1, -1}) {
      CellStem::Labeling best;
      best.rms = 1e300;
      for (int m = 0; m < T.P; ++m) {
        CellStem::Labeling L;
        L.m = m, L.d = d;
        bool ok = true;
        double worst = 0;
        for (int k = 0; k < nbb && ok; ++k) {
          const int sg = wrapi(m + d * k, T.P);
          if (s.u[size_t(k)].el != T.el[size_t(sg)]) { ok = false; break; }
          std::vector<Vec3> tq, uq;
          for (const Side& x : T.side[size_t(sg)]) tq.push_back(d > 0 ? x.q : walked_back(x.q));
          for (size_t i = 0; i < s.u[size_t(k)].side.size(); ++i) uq.push_back(s.side_q(k, int(i)));
          std::vector<int> pm;
          const double dv = match_sides(T.side[size_t(sg)], tq, s.u[size_t(k)].side, uq, &pm);
          if (dv < 0 || dv > 0.6) { ok = false; break; }
          worst = std::max(worst, dv);
          L.perm.push_back(pm);
        }
        if (ok && worst < best.rms) L.rms = worst, best = L;
      }
      if (best.rms < 1e299) s.labels.push_back(best);
    }
    if (s.labels.empty()) throw std::invalid_argument(pc.name + ": stem " + std::to_string(j + 1) + " (" + stems_uc[j].formula + ") does not match the others' repeat unit");
  }
  T.polar = cst[0].labels.size() == 1;
  for (const auto& s : cst)
    if ((s.labels.size() == 1) != T.polar) throw std::invalid_argument(pc.name + ": some stems have a direction and some not");
  {   // the side groups' extra bonds (rings across units) in template terms, from stem 0 read forwards
    const CellStem& s0 = cst[0];
    const CellStem::Labeling* L = nullptr;
    for (const auto& l : s0.labels) if (l.d == 1) L = &l;
    std::set<ExtraBond> ex;
    for (const ExtraBond& e : s0.extra) {
      const int k = e.sigma, k2 = wrapi(e.sigma + e.delta, nbb);
      auto tindex = [&](int kk, int raw) {
        if (raw < 0) return -1;
        const auto& pm = L->perm[size_t(kk)];
        return int(std::find(pm.begin(), pm.end(), raw) - pm.begin());
      };
      ex.insert({wrapi(L->m + k, T.P), tindex(k, e.a), e.delta, tindex(k2, e.b)});
    }
    T.extra.assign(ex.begin(), ex.end());
  }

  // --- the lamella's geometry
  const Vec3 ca = uc.cell.a, cb = uc.cell.b, cc = uc.cell.c;
  if (o.nx < 1 || o.ny < 1) throw std::invalid_argument("lateral repeats must be at least 1");
  const int nc = o.nc > 0 ? o.nc : std::max(2, int(std::lround(o.crystal_nm * 10 / cc[2])));
  const double La = o.amorphous_nm * 10;
  if (La < 5) throw std::invalid_argument("the amorphous layer must be at least 0.5 nm");
  const double Lc = nc * cc[2], Lz = Lc + La;
  const Box box{ca[0] * o.nx, cb[0] * o.ny, cb[1] * o.ny, Lz};
  const double area = box.ax * box.by;
  const double z0 = 0.5 * La;   // the crystal from z0 to z0 + Lc; the amorphous layer across z = 0
  R.nx = o.nx, R.ny = o.ny, R.nc = nc;
  R.lx = box.ax, R.ly = std::hypot(box.bx, box.by), R.lz = Lz;
  R.crystal_thickness = Lc, R.amorphous_thickness = La;
  R.cell_tilt_deg = std::acos(std::clamp(cc[2] / norm(cc), -1.0, 1.0)) / kDeg;
  R.repeat_period = T.P;
  R.polar = T.polar;
  R.polymer = pc.name;
  if (R.cell_tilt_deg > 0.5)
    R.notes.push_back("the crystal's c is " + std::to_string(R.cell_tilt_deg).substr(0, 4) + "° from the stack normal: the fold surfaces are the ab plane and the stems lean by that angle");

  // --- where each stem is cut: at its first backbone atom above the cell's floor, or a little higher where the chain is
  // flexible (the bonds about the cut sp3: CH2-CH2 in PET's glycol or nylon's methylenes), so the loops and bridges that
  // leave and enter there close without bending a ring, an amide or an ester; within c / (stems per cell), at least 4 Å,
  // so the fold surfaces stay flat
  std::vector<int> cut(cst.size(), 0);
  for (size_t j = 0; j < cst.size(); ++j) {
    const CellStem& s = cst[j];
    auto zof = [&](int g) { return (s.u[size_t(wrapi(g, nbb))].pos + s.period * std::floor(double(g) / nbb))[2]; };
    int g0 = -4 * nbb;
    while (zof(g0) < 0) ++g0;
    const double window = std::max(4.0, cc[2] / double(cst.size()));
    auto flex = [&](int g) {   // the bond from unit g − 1 to g
      const int c = s.cls[size_t(wrapi(g - 1, nbb))];
      return c == kSp3 ? 3 : c == kSp2Sp3 ? 2 : c == kConj ? 1 : 0;
    };
    int best = -1;
    for (int t = 0; t < nbb && zof(g0 + t) - zof(g0) <= window; ++t) {
      const int sc = 2 * flex(g0 + t) + flex(g0 + t - 1) + flex(g0 + t + 1) + flex(g0 + t - 2) + flex(g0 + t + 2);
      if (sc > best) best = sc, cut[j] = t;
    }
  }

  // --- the crystal: every stem cut to nc periods
  std::vector<URec> U;
  std::vector<std::vector<int>> stem;   // unit ids bottom → top
  std::vector<int> stem_cs;
  for (int ix = 0; ix < o.nx; ++ix)
    for (int iy = 0; iy < o.ny; ++iy)
      for (size_t j = 0; j < cst.size(); ++j) {
        const CellStem& s = cst[j];
        auto up = [&](int g) {   // global unit g along the stem: position and shift
          const int w = int(std::floor(double(g) / nbb));
          return std::make_pair(wrapi(g, nbb), s.period * double(w));
        };
        int g0 = -4 * nbb;
        while ((s.u[size_t(up(g0).first)].pos + up(g0).second)[2] < 0) ++g0;
        g0 += cut[j];
        std::vector<int> ids;
        for (int g = g0; g < g0 + nc * nbb; ++g) {
          const auto [k, sh] = up(g);
          const Vec3 off = sh + ca * double(ix) + cb * double(iy) + Vec3{0, 0, z0};
          URec r;
          r.el = s.u[size_t(k)].el, r.crystal = true, r.cs = int(j), r.kc = k, r.g = g;
          r.pos = s.u[size_t(k)].pos + off;
          for (const Side& x : s.u[size_t(k)].side) r.side.push_back(x.q + off), r.side_el.push_back(x.el);
          ids.push_back(int(U.size()));
          U.push_back(r);
        }
        stem.push_back(ids);
        stem_cs.push_back(int(j));
      }
  const int N = int(stem.size());
  R.stems = N;
  double ztop = -1e300, zbot = 1e300;
  for (const auto& st : stem) ztop = std::max(ztop, U[size_t(st.back())].pos[2]), zbot = std::min(zbot, U[size_t(st.front())].pos[2]);
  {
    double m = 0;
    for (const auto& u : U) {
      m += element(u.el).mass;
      for (int e : u.side_el) m += element(e).mass;
    }
    R.crystal_density = m / (area * Lc) * kAmuPerA3;
  }

  // --- the amorphous layer: backbone units at the target density
  const int M = int(std::lround(o.density * area * La / kAmuPerA3 / T.unit_mass));

  std::mt19937_64 rng(o.seed);
  // --- segment counts: S segments with fractions f use 2N ends: S (2 − f_tail) = 2N
  double fl = std::max(0.0, o.loops), fb = std::max(0.0, o.bridges), ft = std::max(0.0, o.tails);
  const double fs = fl + fb + ft;
  if (fs <= 0) throw std::invalid_argument("the loop, bridge and tail fractions are all zero");
  fl /= fs, fb /= fs, ft /= fs;
  const double S = 2.0 * N / (2.0 - ft);
  int nt = 2 * int(std::lround(ft * S / 2));
  int nbr = int(std::lround(fb * S));
  nt = std::max(nt, 2);   // chains need ends
  nbr = std::min(nbr, N - 1);

  // the direction of each stem's chain (polar chains): d = +1 when the template runs up the stem; an end's sign is the
  // direction the chain leaves the stem by (top: d, bottom: −d), and a loop or bridge joins two ends of opposite signs
  auto stem_d = [&](int i) { return T.polar ? cst[size_t(stem_cs[size_t(i)])].labels[0].d : 0; };
  auto sign_of = [&](int e) { return e % 2 ? stem_d(e / 2) : -stem_d(e / 2); };
  auto compat = [&](int e, int f) { return !T.polar || sign_of(e) != sign_of(f); };
  auto end_unit = [&](int e) { return e % 2 ? stem[size_t(e / 2)].back() : stem[size_t(e / 2)].front(); };
  auto lateral = [&](int e1, int e2) { return box.lateral(U[size_t(end_unit(e2))].pos - U[size_t(end_unit(e1))].pos); };

  // --- pair the stem ends: end e = 2 stem + (1 top | 0 bottom)
  std::vector<Seg> seg;
  std::vector<int> seg_of(2 * N, -1);
  std::vector<int> topE(N), botE(N);
  for (int i = 0; i < N; ++i) topE[i] = 2 * i + 1, botE[i] = 2 * i;
  shuffle(topE.begin(), topE.end(), rng);
  shuffle(botE.begin(), botE.end(), rng);
  std::vector<char> used(2 * N, 0);
  Normal<double> gauss;
  const double sig = std::max(5.0, 0.4 * La);
  int nbr_made = 0;
  for (int k = 0; k < N && nbr_made < nbr; ++k) {
    const int e = topE[k];
    const Vec3 p = U[size_t(end_unit(e))].pos;
    const Vec3 target{p[0] + sig * gauss(rng), p[1] + sig * gauss(rng), 0};
    int best = -1;
    double bd = 1e300;
    for (int f : botE) {
      if (used[f] || !compat(e, f)) continue;
      Vec3 d = box.mi(U[size_t(end_unit(f))].pos - target);
      d[2] = 0;
      if (norm(d) < bd) bd = norm(d), best = f;
    }
    if (best < 0) continue;
    used[e] = used[best] = 1;
    seg_of[e] = seg_of[best] = int(seg.size());
    seg.push_back(Seg(1, e, best));
    ++nbr_made;
  }
  if (nbr_made < nbr) R.notes.push_back(std::to_string(nbr_made) + " bridges of the " + std::to_string(nbr) + " asked (no free end of the other direction left)");
  if (o.reentry != "adjacent" && o.reentry != "random") throw std::invalid_argument("reentry: adjacent or random");
  // tails per face: as asked, raised until the rest pairs up (an even count; for chains with a direction as many
  // ends leaving up the template as down)
  int ntT = nt / 2, ntB = nt - ntT;
  auto face = [&](std::vector<int>& ends, int ntails) {
    std::vector<int> rest;
    for (int e : ends) if (!used[e]) rest.push_back(e);
    int plus = 0;
    for (int e : rest) plus += sign_of(e) > 0;
    const int minus = int(rest.size()) - plus;
    auto feasible = [&](int t) {
      if (t > int(rest.size())) return true;
      if (!T.polar) return (int(rest.size()) - t) % 2 == 0;
      return t >= std::abs(plus - minus) && (t - std::abs(plus - minus)) % 2 == 0;
    };
    while (!feasible(ntails)) ++ntails;
    ntails = std::min(ntails, int(rest.size()));
    // which ends become tails: for directed chains as many of each sign as balances the rest
    int tp = T.polar ? (ntails + plus - minus) / 2 : 0, tm = T.polar ? ntails - tp : 0;
    std::vector<int> tails, loopers;
    for (int e : rest) {
      bool tail = false;
      if (!T.polar) tail = int(tails.size()) < ntails;
      else if (sign_of(e) > 0 && tp > 0) tail = true, --tp;
      else if (sign_of(e) < 0 && tm > 0) tail = true, --tm;
      (tail ? tails : loopers).push_back(e);
    }
    for (int e : tails) {
      used[e] = 1;
      seg_of[e] = int(seg.size());
      seg.push_back(Seg(2, e, -1));
    }
    UniformReal<double> u01(0.0, 1.0);
    for (size_t k = 0; k < loopers.size(); ++k) {
      const int e = loopers[k];
      if (used[e]) continue;
      int partner = -1;
      std::vector<int> free;
      for (size_t q = k + 1; q < loopers.size(); ++q) if (!used[loopers[q]] && compat(e, loopers[q])) free.push_back(loopers[q]);
      if (free.empty()) throw std::logic_error("unpaired loop end");
      if (o.reentry == "random" || u01(rng) >= 0.7) {
        partner = free[UniformInt<size_t>(0, free.size() - 1)(rng)];
      } else {
        double bd = 1e300;
        for (int f : free) if (lateral(e, f) < bd) bd = lateral(e, f), partner = f;
      }
      used[e] = used[partner] = 1;
      seg_of[e] = seg_of[partner] = int(seg.size());
      seg.push_back(Seg(0, e, partner));
    }
    return ntails;
  };
  const int gotT = face(topE, ntT), gotB = face(botE, ntB);
  if (gotT + gotB > nt) R.notes.push_back(std::to_string(gotT + gotB - nt) + " more tails than asked, so that the rest of the stem ends pair up");
  if (gotT + gotB == 0) throw std::invalid_argument("chains need tails");

  // --- no chain may close on itself: a cycle (stems and segments with no tail) is joined to a chain with tails by
  // exchanging the partners of one of its segments and a segment of the same kind outside it
  auto partner_of = [&](int e) {
    const Seg& g = seg[size_t(seg_of[e])];
    return g.type == 2 ? -1 : (g.a == e ? g.b : g.a);
  };
  int merges = 0, broken = 0;
  for (int pass = 0; pass < 4 * N; ++pass) {
    std::vector<int> comp(N, -1);
    int ncomp = 0;
    std::vector<char> has_tail;
    for (int s0 = 0; s0 < N; ++s0) {
      if (comp[s0] >= 0) continue;
      bool tail = false;
      std::vector<int> st{s0};
      comp[s0] = ncomp;
      while (!st.empty()) {
        const int s1 = st.back();
        st.pop_back();
        for (int e : {2 * s1, 2 * s1 + 1}) {
          const int p = partner_of(e);
          if (p < 0) { tail = true; continue; }
          if (comp[p / 2] < 0) comp[p / 2] = ncomp, st.push_back(p / 2);
        }
      }
      has_tail.push_back(tail);
      ++ncomp;
    }
    int cyc = -1;
    for (int k = 0; k < ncomp; ++k) if (!has_tail[size_t(k)]) { cyc = k; break; }
    if (cyc < 0) break;
    int bx = -1, by = -1;
    double bd = 1e300;
    for (size_t x = 0; x < seg.size(); ++x) {
      if (seg[x].type == 2 || comp[seg[x].a / 2] != cyc) continue;
      for (size_t y = 0; y < seg.size(); ++y) {
        if (seg[y].type != seg[x].type || !has_tail[size_t(comp[seg[y].a / 2])]) continue;
        if (seg[x].type == 0 && (seg[x].a % 2) != (seg[y].a % 2)) continue;
        // the exchanged ends of the same sign, so both new segments still join opposite signs
        if (seg[x].type == 0 ? sign_of(seg[y].a) != sign_of(seg[x].b) : sign_of(seg[y].b) != sign_of(seg[x].b)) continue;
        const double d = lateral(seg[x].a, seg[y].a);
        if (d < bd) bd = d, bx = int(x), by = int(y);
      }
    }
    if (bx < 0) {   // nothing to exchange with: the cycle opens into two tails
      for (size_t x = 0; x < seg.size(); ++x)
        if (seg[x].type != 2 && comp[seg[x].a / 2] == cyc) {
          const int e1 = seg[x].a, e2 = seg[x].b;
          seg[x] = Seg(2, e1, -1);
          seg_of[e2] = int(seg.size());
          seg.push_back(Seg(2, e2, -1));
          ++broken;
          break;
        }
      continue;
    }
    Seg &X1 = seg[size_t(bx)], &Y1 = seg[size_t(by)];
    if (X1.type == 1) {   // bridges top → bottom: exchange the bottom ends
      std::swap(X1.b, Y1.b);
      seg_of[X1.b] = bx, seg_of[Y1.b] = by;
    } else {
      std::swap(X1.b, Y1.a);
      seg_of[X1.b] = bx, seg_of[Y1.a] = by;
    }
    ++merges;
  }
  if (merges) R.notes.push_back(std::to_string(merges) + " segment pairs exchanged partners so that no chain closes on itself");
  if (broken) R.notes.push_back(std::to_string(broken) + " closed rings opened into two tails each (no segment of their kind to exchange with)");

  // --- each stem's reading of the template, chain by chain (a stem of a chain with no direction takes whichever its
  // chain needs), then every crystal unit's sigma and side order
  std::vector<int> lab(N, -1);
  auto lab_with = [&](int i, int d) {
    const auto& L = cst[size_t(stem_cs[size_t(i)])].labels;
    for (size_t k = 0; k < L.size(); ++k) if (L[k].d == d) return int(k);
    throw std::logic_error("a stem without the direction its chain needs");
  };
  auto dir = [&](int i) { return cst[size_t(stem_cs[size_t(i)])].labels[size_t(lab[size_t(i)])].d; };
  for (size_t g0 = 0; g0 < seg.size(); ++g0) {
    if (seg[g0].type != 2) continue;
    int e = seg[g0].a;
    if (lab[size_t(e / 2)] >= 0) continue;   // this chain was walked from its other tail
    lab[size_t(e / 2)] = T.polar ? 0 : lab_with(e / 2, 1);
    for (;;) {
      // the chain runs through the stem from end e to e ^ 1 and leaves it with sign D
      const int other = e ^ 1;
      const int D = other % 2 ? dir(e / 2) : -dir(e / 2);
      const Seg& G = seg[size_t(seg_of[other])];
      if (G.type == 2) break;
      const int nxt = G.a == other ? G.b : G.a;
      // entering at the top the chain walks down the stem: −d = D; at the bottom d = D
      const int need = nxt % 2 ? -D : D;
      if (lab[size_t(nxt / 2)] >= 0) throw std::logic_error("a chain met a stem twice");
      lab[size_t(nxt / 2)] = lab_with(nxt / 2, need);
      e = nxt;
    }
  }
  for (int i = 0; i < N; ++i) {
    if (lab[size_t(i)] < 0) throw std::logic_error("a stem belongs to no chain");
    const auto& L = cst[size_t(stem_cs[size_t(i)])].labels[size_t(lab[size_t(i)])];
    for (int id : stem[size_t(i)]) {
      URec& u = U[size_t(id)];
      u.sigma = wrapi(L.m + L.d * u.g, T.P);
      const auto& pm = L.perm[size_t(u.kc)];
      std::vector<Vec3> sd(pm.size());
      for (size_t t = 0; t < pm.size(); ++t) sd[t] = u.side[size_t(pm[t])];
      u.side = sd;
      u.side_el.clear();
      for (const Side& x : T.side[size_t(u.sigma)]) u.side_el.push_back(x.el);
    }
  }

  // --- target end positions (unwrapped: the bottom face is reached in the image above, z + Lz)
  auto end_units = [&](int e, int& u1, int& u2, int& u3) {
    const auto& st = stem[size_t(e / 2)];
    if (e % 2) u1 = st[st.size() - 1], u2 = st[st.size() - 2], u3 = st[st.size() - 3];
    else u1 = st[0], u2 = st[1], u3 = st[2];
  };
  auto shifted = [&](int e) { return Vec3{0, 0, e % 2 ? 0.0 : Lz}; };
  auto image_near = [&](const Vec3& p, const Vec3& ref) {   // the lateral image of p nearest ref (z kept)
    Vec3 d = p - ref;
    const double dz = d[2];
    d = box.mi(d);
    d[2] = dz;
    return ref + d;
  };
  // grow every loop and bridge from its top end (bridges: top → bottom; bottom loops start at their first end)
  std::vector<int> nmin(seg.size(), 1);
  const double d13 = [&] {
    double s = 0;
    for (int k = 0; k < T.P; ++k) {
      const double r1 = T.r[size_t(k)], r2 = T.r[size_t(wrapi(k - 1, T.P))];
      s += std::sqrt(r1 * r1 + r2 * r2 - 2 * r1 * r2 * std::cos(T.th[size_t(k)])) / T.P;
    }
    return s;
  }();
  for (size_t g = 0; g < seg.size(); ++g) {
    Seg& G = seg[g];
    if (G.type == 1 && G.a % 2 == 0) std::swap(G.a, G.b);   // bridges grow from the top face
    int a1, a2, a3;
    end_units(G.a, a1, a2, a3);
    G.d = G.a % 2 ? dir(G.a / 2) : -dir(G.a / 2);
    if (G.type == 2) { nmin[g] = 1; continue; }
    int b1, b2, b3;
    end_units(G.b, b1, b2, b3);
    G.q = wrapi(G.d * (U[size_t(b1)].sigma - U[size_t(a1)].sigma) - 1, T.P);
    const Vec3 s1 = U[size_t(a1)].pos + shifted(G.a), e1 = image_near(U[size_t(b1)].pos + shifted(G.b), s1);
    const double chord = norm(e1 - s1);
    int n = G.type == 1 ? int(std::ceil(chord / (0.72 * T.reach))) : std::max(4, int(std::ceil(chord / (0.6 * T.reach))) + 2);
    if (T.P >= 6) n = std::max(n, G.type == 1 ? n + T.P : 2 * T.P);   // a long, stiff repeat (PET, nylons) folds in two at least
    n += wrapi(G.q - n, T.P);
    nmin[g] = n;
  }
  const int need = std::accumulate(nmin.begin(), nmin.end(), 0);
  if (need > M)
    throw std::invalid_argument("the amorphous layer holds " + std::to_string(M) + " backbone atoms at " + std::to_string(o.density).substr(0, 5) +
                                " g/cm³ but the loops, bridges and tails need at least " + std::to_string(need) +
                                ": a thicker amorphous layer, fewer bridges, or loops that re-enter nearby");
  {   // the rest shared out: exponential lengths (loops and bridges weight 1, tails 1.5: their free ends fill the middle of
      // the layer); loops and bridges in whole repeat units, so each closes in phase with the stem it enters
    std::vector<double> w(seg.size());
    UniformReal<double> u01(1e-12, 1.0);
    double sw = 0;
    for (size_t g = 0; g < seg.size(); ++g) sw += w[g] = -std::log(u01(rng)) * (seg[g].type == 2 ? 1.5 : 1.0);
    const int left = M - need;
    std::vector<int> extra(seg.size());
    int given = 0;
    std::vector<size_t> tails;
    for (size_t g = 0; g < seg.size(); ++g) {
      int x = int(std::floor(left * w[g] / sw));
      if (seg[g].type != 2) x -= x % T.P;
      else tails.push_back(g);
      given += extra[g] = x;
    }
    for (int k = given; k < left; ++k) extra[tails[UniformInt<size_t>(0, tails.size() - 1)(rng)]]++;
    for (size_t g = 0; g < seg.size(); ++g) seg[g].n = nmin[g] + extra[g];
    // a tail ends where the chain may be cut: after an sp3-sp3 bond (CH2-CH2, so the end is CH3 or OCH2CH3), or, in a
    // chain with none (PLLA), next to one trigonal atom; never in a ring, an amide or ester bond, or at a C=C
    int best_cls = -1;
    for (int k = 0; k < T.P; ++k)
      if (T.cls[size_t(k)] == kSp3 || T.cls[size_t(k)] == kSp2Sp3) best_cls = std::max(best_cls, T.cls[size_t(k)] == kSp3 ? 1 : 0);
    if (best_cls < 0) throw std::invalid_argument(pc.name + ": no backbone bond a chain end may be cut at (every one in a ring, an amide or ester, or double)");
    const int cut_cls = best_cls == 1 ? kSp3 : kSp2Sp3;
    std::vector<std::vector<int>> cuts(seg.size());   // each tail's cuttable lengths up to its share and a repeat beyond
    for (size_t g = 0; g < seg.size(); ++g) {
      Seg& G = seg[g];
      if (G.type != 2) continue;
      int a1, a2, a3;
      end_units(G.a, a1, a2, a3);
      // the last unit sigma_a + d n, cut from the next one: the bond (s, s + d)
      auto cuttable = [&](int n) {
        const int sl = wrapi(U[size_t(a1)].sigma + G.d * n, T.P);
        return T.cls[size_t(G.d > 0 ? sl : wrapi(sl - 1, T.P))] == cut_cls;
      };
      for (int n = 1; n <= std::max(1, G.n) + T.P; ++n) if (cuttable(n)) cuts[g].push_back(n);
      int best = cuts[g].back();   // the cuttable length nearest the tail's share
      for (int n : cuts[g]) if (std::abs(n - G.n) < std::abs(best - G.n)) best = n;
      G.n = best;
    }
    // back to the layer's budget: the longest tails shortened to their next cuttable length while that brings the total nearer
    for (;;) {
      int total = 0;
      for (const auto& G : seg) total += G.n;
      size_t lg = seg.size();
      for (size_t g = 0; g < seg.size(); ++g)
        if (seg[g].type == 2 && cuts[g].front() < seg[g].n && (lg == seg.size() || seg[g].n > seg[lg].n)) lg = g;
      if (lg == seg.size()) break;
      int lower = cuts[lg].front();
      for (int n : cuts[lg]) if (n < seg[lg].n) lower = n;
      if (std::abs(total - (seg[lg].n - lower) - M) >= std::abs(total - M)) break;
      seg[lg].n = lower;
    }
  }

  // --- growth: one backbone atom per segment per round, in random order; each atom placed brings the side atoms of the
  // one before it (their frame needs both neighbours)
  Hash hash(box, std::max(o.contact, 3.0));
  for (size_t i = 0; i < U.size(); ++i) {
    hash.add(U[i].pos, int32_t(i));
    for (size_t k = 0; k < U[i].side.size(); ++k)
      if (U[i].side_el[k] > 1) U[i].hs.push_back(hash.add(U[i].side[k], int32_t(i)));
  }
  const double zlo = ztop + 0.3, zhi = zbot + Lz - 0.3;
  struct Grow {
    std::vector<Vec3> path;     // the last atoms (stem end atoms first), unwrapped
    std::vector<int> ids;       // their units
    std::vector<uint32_t> hid;  // hash id of each bead
    Vec3 e1, e2;                // closure target and its stem neighbour (unwrapped image)
    int te1 = 0, te2 = 0, te3 = 0;
    int placed = 0;
  };
  std::vector<Grow> gs(seg.size());
  for (size_t g = 0; g < seg.size(); ++g) {
    int a1, a2, a3;
    end_units(seg[g].a, a1, a2, a3);
    const Vec3 sh = shifted(seg[g].a);
    gs[g].path = {U[size_t(a3)].pos + sh, U[size_t(a2)].pos + sh, U[size_t(a1)].pos + sh};
    gs[g].ids = {a3, a2, a1};
    if (seg[g].type != 2) {
      end_units(seg[g].b, gs[g].te1, gs[g].te2, gs[g].te3);
      const Vec3 e1 = U[size_t(gs[g].te1)].pos + shifted(seg[g].b), e2 = U[size_t(gs[g].te2)].pos + shifted(seg[g].b);
      const Vec3 im = image_near(e1, gs[g].path[2]) - e1;
      gs[g].e1 = e1 + im, gs[g].e2 = e2 + im;
    }
  }
  double worst = 1e300;
  int eased = 0, strained = 0;
  UniformReal<double> u01(0.0, 1.0);
  const double kC = 6.7;   // characteristic ratio of polyethylene (Flory), for the closure bias
  // the layer filled evenly: trials weighted down where the atoms grown so far already exceed their share (2 Å slices)
  const int nzb = std::max(1, int((zhi - zlo) / 2.0));
  std::vector<double> zfill(size_t(nzb), 0.0);
  const double zper = double(M) / nzb;
  double zsum = 0;
  auto zb = [&](double z) { return size_t(std::clamp(int((z - zlo) / (zhi - zlo) * nzb), 0, nzb - 1)); };
  std::vector<int> backtracks(seg.size(), 0);
  int backtracked = 0;
  // closest approach of p to the atoms placed, those of the excluded units left out, minus the limit
  auto margin_of = [&](const Vec3& p, const std::vector<int>& excl, double limit) {
    double m = 1e300;
    hash.near(p, [&](uint32_t id, double d) {
      if (d > limit + 1.0) return;
      for (int e : excl) if (e == hash.owner[id]) return;
      m = std::min(m, d - limit);
    });
    return m;
  };
  // the side atoms of unit u (sigma known), from its neighbours as walked in direction d, into the hash
  auto set_sides = [&](int u, int d, const Vec3& prev, const Vec3& at, const Vec3& next) {
    URec& x = U[size_t(u)];
    for (uint32_t h : x.hs) hash.dead[h] = 1;
    x.hs.clear();
    x.side = T.sides(x.sigma, d, prev, at, next);
    x.side_el.clear();
    for (const Side& s : T.side[size_t(x.sigma)]) x.side_el.push_back(s.el);
    for (size_t k = 0; k < x.side.size(); ++k)
      if (x.side_el[k] > 1) x.hs.push_back(hash.add(x.side[k], u));
  };
  auto clear_sides = [&](int u) {
    for (uint32_t h : U[size_t(u)].hs) hash.dead[h] = 1;
    U[size_t(u)].hs.clear();
    U[size_t(u)].side.clear();
  };
  // a torsion for the bond class: held bonds at the crystal's value, aryl-carbonyl at it or turned by 180°, others from
  // a trans-gauche set, with a share at the crystal's value where that is not trans (helices)
  auto torsion = [&](int cls, double ref) {
    if (cls == kFixed) return ref + 5.0 * kDeg * gauss(rng);
    if (cls == kConj) return (u01(rng) < 0.5 ? ref : ref + kPi) + 10.0 * kDeg * gauss(rng);
    if (std::fabs(ref) < 150 * kDeg && u01(rng) < 0.35) return ref + 15.0 * kDeg * gauss(rng);
    const double r = u01(rng);
    if (cls == kSp2Sp3) return ((r < 0.34 ? 180.0 : r < 0.67 ? 100.0 : -100.0) + 20.0 * gauss(rng)) * kDeg;
    return ((r < 0.55 ? 180.0 : r < 0.775 ? 65.0 : -65.0) + 15.0 * gauss(rng)) * kDeg;
  };
  std::vector<size_t> active;
  for (size_t g = 0; g < seg.size(); ++g) if (seg[g].n > 0) active.push_back(g);
  int64_t done_beads = 0;
  while (!active.empty()) {
    shuffle(active.begin(), active.end(), rng);
    std::vector<size_t> still;
    for (size_t g : active) {
      Seg& G = seg[g];
      Grow& W = gs[g];
      const bool closing = G.type != 2;
      const int k = G.n - W.placed - 1;   // atoms after this one
      const size_t L = W.path.size();
      const Vec3 A = W.path[L - 3], B = W.path[L - 2], C = W.path[L - 1];
      const int hc = W.ids[L - 1];
      const int sc = U[size_t(hc)].sigma, sp = T.next(sc, G.d);
      const double rb = T.bond(sc, G.d), th = T.th[size_t(sc)];
      std::vector<int> excl(W.ids.end() - 3, W.ids.end());
      Vec3 P;
      if (closing && k == 0) {
        // the last atom closes on the target: on the circle of points one bond from both the head and the target
        excl.push_back(W.te1), excl.push_back(W.te2);
        const Vec3 E = W.e1, dv = E - C;
        const double dd = norm(dv), r2 = T.bond(sp, G.d);
        const double thE = T.th[size_t(U[size_t(W.te1)].sigma)], thP = T.th[size_t(sp)];
        if (dd < rb + r2 - 0.02 && dd > std::fabs(rb - r2)) {
          const Vec3 ax = unit(dv);
          const double t0 = (dd * dd + rb * rb - r2 * r2) / (2 * dd);
          const double rho = std::sqrt(std::max(0.0, rb * rb - t0 * t0));
          const Vec3 m = C + ax * t0;
          const Vec3 uu = unit(cross(ax, std::fabs(ax[0]) < 0.9 ? Vec3{1, 0, 0} : Vec3{0, 1, 0})), vv = cross(ax, uu);
          double best = 1e300;
          for (int q = 0; q < 72; ++q) {
            const double t = q * 5 * kDeg;
            const Vec3 x = m + (uu * std::cos(t) + vv * std::sin(t)) * rho;
            const double t1 = (angle_deg(B, C, x) - th / kDeg), t2 = (angle_deg(x, E, W.e2) - thE / kDeg), t3 = (angle_deg(C, x, E) - thP / kDeg);
            const double mg = margin_of(x, excl, o.contact * 0.8);
            const double scr = (t1 * t1 + t2 * t2 + t3 * t3) / 100 + (mg < 0 ? 100 * mg * mg : 0) + ((x[2] < zlo || x[2] > zhi) ? 50 : 0);
            if (scr < best) best = scr, P = x;
          }
          if (best > 20) ++strained;
        } else {
          P = (C + E) * 0.5;
          ++strained;
        }
      } else {
        const bool dense = closing && k + 1 <= (T.P > 2 ? 5 : 3);   // a long repeat (stiff units): a wider scan further out
        const int ntr = dense ? 120 : std::max(4, o.trials);
        const int tc = T.tcls(sc, G.d);
        const double tr = T.tref(sc, G.d);
        std::vector<Vec3> cand(static_cast<size_t>(ntr));
        std::vector<double> dist(size_t(ntr), 0.0);
        for (int t = 0; t < ntr; ++t) {
          double phi;
          if (dense) {
            // near the target a held bond may bend ±20° (the relaxation brings it back) so a stiff unit can still close
            if (tc == kFixed) phi = tr + (t - ntr / 2) * (40.0 / ntr) * kDeg;
            else if (tc == kConj) phi = tr + (t % 2) * kPi + (t / 2 - ntr / 4) * (80.0 / ntr) * kDeg;
            else phi = (t * 3.0 + 1.5 * u01(rng)) * kDeg;
          } else phi = torsion(tc, tr);
          cand[size_t(t)] = place(A, B, C, rb, th, phi);
        }
        // hard limits: inside the amorphous layer, and (closing) the target still reachable with the bonds left
        std::vector<char> ok(size_t(ntr), 1);
        std::vector<double> wcl(size_t(ntr), 1.0);
        double d13here = d13;
        if (closing && k + 1 == 2) {
          const int s2 = T.next(sp, G.d);
          const double r1 = T.bond(sp, G.d), r2 = T.bond(s2, G.d), a2 = T.th[size_t(s2)];
          d13here = std::sqrt(r1 * r1 + r2 * r2 - 2 * r1 * r2 * std::cos(a2));
        }
        for (int t = 0; t < ntr; ++t) {
          const Vec3& p = cand[size_t(t)];
          if (p[2] < zlo || p[2] > zhi) ok[size_t(t)] = 0;
          if (closing) {
            const int kb = k + 1;   // bonds from this atom to the target
            const double d = norm(W.e1 - p);
            dist[size_t(t)] = d;
            if (d > kb * T.reach * (kb > 6 ? 0.8 : kb > 3 ? 0.9 : kb > 2 ? 0.97 : 1.0)) ok[size_t(t)] = 0;
            if (kb == 2 && (d < d13here - 0.33 || d > d13here + 0.31)) ok[size_t(t)] = 0;
            if (kb == 3 && d < 0.9 * d13) ok[size_t(t)] = 0;
            // ideal-chain closure probability, and a steep penalty once the target needs more than half the remaining reach
            const double f = d / (kb * T.reach);
            wcl[size_t(t)] = std::exp(-3 * d * d / (2.0 * kb * kC * T.r_mean * T.r_mean) - (f > 0.5 ? 25 * (f - 0.5) : 0.0));
          }
        }
        std::vector<int> ex2 = excl;
        if (closing && k + 1 <= 3) ex2.push_back(W.te1), ex2.push_back(W.te2);
        if (closing && k + 1 <= 2) ex2.push_back(W.te3);
        // the head's side atoms for each trial (they hang between the atom before and the trial)
        std::vector<std::vector<Vec3>> heavy(static_cast<size_t>(ntr));
        const auto& hside = T.side[size_t(sc)];
        bool has_heavy = false;
        for (const Side& x : hside) has_heavy |= x.el > 1;
        if (has_heavy)
          for (int t = 0; t < ntr; ++t) {
            const auto sd = T.sides(sc, G.d, B, C, cand[size_t(t)]);
            for (size_t i = 0; i < sd.size(); ++i) if (hside[i].el > 1) heavy[size_t(t)].push_back(sd[i]);
          }
        auto clear_of = [&](int t, double lim) {
          if (margin_of(cand[size_t(t)], ex2, lim) < 0) return false;
          for (const Vec3& h : heavy[size_t(t)]) if (margin_of(h, ex2, lim) < 0) return false;
          return true;
        };
        int pick = -1;
        static const double scales[] = {1.0, 0.9, 0.8, 0.7, 0.6, 0.5};
        for (int si = 0; si < 6 && pick < 0; ++si) {
          const double lim = o.contact * scales[si];
          std::vector<double> w(size_t(ntr), 0.0);
          double sw = 0;
          for (int t = 0; t < ntr; ++t) {
            if (!ok[size_t(t)]) continue;
            if (!clear_of(t, lim)) continue;
            const double fill = zfill[zb(cand[size_t(t)][2])] / std::max(zper * 0.25, zsum / nzb);
            if (si < 3 && fill > 1.15) continue;   // a full slice: room is looked for further out first
            // loops (their first half) and tails lean away from their face, so the interphase is not packed with short folds
            double lean = 0;
            if (G.type == 0 && 2 * W.placed < G.n) lean = 0.8;
            else if (G.type == 2) lean = 0.6;
            const double upw = (G.a % 2 ? 1.0 : -1.0) * (cand[size_t(t)][2] - C[2]);
            sw += w[size_t(t)] = wcl[size_t(t)] * std::exp(-8.0 * (fill - 1.0) + lean * upw);
          }
          if (sw <= 0) continue;
          double r = u01(rng) * sw;
          for (int t = 0; t < ntr; ++t) {
            r -= w[size_t(t)];
            if (w[size_t(t)] > 0 && r <= 0) { pick = t; break; }
          }
          if (pick < 0) for (int t = ntr - 1; t >= 0; --t) if (w[size_t(t)] > 0) { pick = t; break; }
          if (si > 0) ++eased;
        }
        if (pick < 0 && W.placed > 0 && backtracks[g] < (closing ? 240 : 80) && std::none_of(ok.begin(), ok.end(), [](char c) { return c != 0; })) {
          // a dead end (the target out of reach or the layer's faces in the way): take back a few atoms and grow them again
          const int back = std::min(W.placed, 2 + backtracks[g] / 8);
          for (int q = 0; q < back; ++q) {
            const int u = G.beads.back();
            hash.dead[W.hid.back()] = 1;
            zfill[zb(U[size_t(u)].pos[2])] -= 1;
            zsum -= 1;
            G.beads.pop_back();
            W.path.pop_back();
            W.ids.pop_back();
            W.hid.pop_back();
            clear_sides(W.ids.back());   // the atom before had its side atoms from the one taken back
            --W.placed;
          }
          ++backtracks[g];
          ++backtracked;
          still.push_back(g);
          continue;
        }
        if (pick < 0) {   // nothing clears even half the limit: the roomiest trial that keeps the hard limits, or the nearest to them
          double best = -1e300;
          for (int t = 0; t < ntr; ++t) {
            const double scr = margin_of(cand[size_t(t)], ex2, o.contact) - (ok[size_t(t)] ? 0 : 100 + std::fabs(dist[size_t(t)]));
            if (scr > best) best = scr, pick = t;
          }
          ++eased;
          if (!ok[size_t(pick)]) ++strained;
        }
        P = cand[size_t(pick)];
      }
      worst = std::min(worst, margin_of(P, excl, 0.0));
      const int id = int(U.size());
      URec nu;
      nu.el = T.el[size_t(sp)], nu.sigma = sp, nu.pos = P;
      U.push_back(nu);
      set_sides(hc, G.d, B, C, P);
      W.hid.push_back(hash.add(P, id));
      zfill[zb(P[2])] += 1;
      zsum += 1;
      G.beads.push_back(id);
      W.path.push_back(P);
      W.ids.push_back(id);
      ++W.placed;
      ++done_beads;
      if (closing && W.placed == G.n) {
        // closed: the last atom's side atoms and the target's, from their new neighbours
        set_sides(id, G.d, C, P, W.e1);
        set_sides(W.te1, G.d, P, W.e1, W.e2);
      } else if (!closing && W.placed == G.n) {
        // a free end: its side atoms from where the next atom would be, and a hydrogen there
        const int tc = T.tcls(sp, G.d);
        const Vec3 Q = place(B, C, P, T.bond(sp, G.d), T.th[size_t(sp)], tc == kSp3 || tc == kSp2Sp3 ? kPi : T.tref(sp, G.d));
        set_sides(id, G.d, C, P, Q);
        U[size_t(id)].cap = true;
        U[size_t(id)].cap_pos = P + unit(Q - P) * bond_to_h(U[size_t(id)].el);
      }
      if (W.placed < G.n) still.push_back(g);
    }
    active.swap(still);
    if (cancelled("growing the amorphous layer", double(done_beads) / std::max(1, M))) throw std::runtime_error("cancelled");
  }
  R.worst_contact = worst > 1e299 ? 0 : worst;
  if (backtracked) R.notes.push_back(std::to_string(backtracked) + " backtracks (atoms taken back and grown again at dead ends)");
  R.eased_steps = eased, R.strained_closures = strained;

  // --- chains: from each tail through stems and segments to the other tail (units in order, and the template's
  // direction along the walk)
  std::vector<std::vector<int>> chain;
  std::vector<int> chain_dir;
  std::vector<char> seen(seg.size(), 0);
  for (size_t g0 = 0; g0 < seg.size(); ++g0) {
    if (seg[g0].type != 2 || seen[g0]) continue;
    std::vector<int> ch(seg[g0].beads.rbegin(), seg[g0].beads.rend());
    seen[g0] = 1;
    int e = seg[g0].a;
    chain_dir.push_back(-seg[g0].d);
    for (;;) {
      const auto& st = stem[size_t(e / 2)];
      if (e % 2) ch.insert(ch.end(), st.rbegin(), st.rend());
      else ch.insert(ch.end(), st.begin(), st.end());
      const int other = e ^ 1;
      const int gi = seg_of[other];
      Seg& G = seg[size_t(gi)];
      seen[size_t(gi)] = 1;
      if (G.type == 2) {
        ch.insert(ch.end(), G.beads.begin(), G.beads.end());
        break;
      }
      if (G.a == other) {
        ch.insert(ch.end(), G.beads.begin(), G.beads.end());
        e = G.b;
      } else {
        ch.insert(ch.end(), G.beads.rbegin(), G.beads.rend());
        e = G.a;
      }
    }
    chain.push_back(ch);
  }
  for (size_t g = 0; g < seg.size(); ++g)
    if (!seen[g]) throw std::logic_error("a segment belongs to no chain (a ring left)");
  for (size_t ci = 0; ci < chain.size(); ++ci)
    for (size_t k = 1; k < chain[ci].size(); ++k)
      if (U[size_t(chain[ci][k])].sigma != T.next(U[size_t(chain[ci][k - 1])].sigma, chain_dir[ci]))
        throw std::logic_error("a chain out of step with its repeat unit");

  // --- each unit's class and segment as built: stems first (segment = stem), then the loops, bridges and tails
  std::vector<int8_t> ucls(U.size(), SemicrystallineReport::kNone);
  std::vector<int32_t> useg(U.size(), -1);
  for (int i = 0; i < N; ++i)
    for (int u : stem[size_t(i)]) ucls[size_t(u)] = SemicrystallineReport::kStem, useg[size_t(u)] = i;
  for (size_t g = 0; g < seg.size(); ++g)
    for (int u : seg[g].beads)
      ucls[size_t(u)] = seg[g].type == 0 ? SemicrystallineReport::kLoop : seg[g].type == 1 ? SemicrystallineReport::kBridge : SemicrystallineReport::kTail,
      useg[size_t(u)] = N + int32_t(g);

  // --- the System: backbone atoms in chain order, each followed by its side atoms
  System s;
  s.title = "semicrystalline " + pc.name + " lamellar stack";
  s.cell.a = {box.ax, 0, 0}, s.cell.b = {box.bx, box.by, 0}, s.cell.c = {0, 0, Lz};
  s.has_mol = true;
  s.tags = {{"crystal", "#3a7bd5"}, {"amorphous", "#e0913a"}};
  for (size_t ci = 0; ci < chain.size(); ++ci) {
    const auto& ch = chain[ci];
    const int D = chain_dir[ci];
    std::vector<uint32_t> first(ch.size()), base(ch.size());
    for (size_t k = 0; k < ch.size(); ++k) {
      const URec& u = U[size_t(ch[k])];
      if (u.side.size() != T.side[size_t(u.sigma)].size()) throw std::logic_error("a backbone atom without its side atoms");
      Atom at;
      at.element = u.el, at.name = element(u.el).symbol, at.mol = int64_t(ci) + 1;
      at.pos = s.cell.wrap(u.pos);
      at.tags = u.crystal ? 1u : 2u;
      base[k] = uint32_t(s.atoms.size());
      s.atoms.push_back(at);
      if (k) s.bonds.push_back({base[k - 1], base[k], 0});
      const auto& ts = T.side[size_t(u.sigma)];
      for (size_t i = 0; i < u.side.size(); ++i) {
        Atom sa = at;
        sa.element = ts[i].el, sa.name = element(ts[i].el).symbol, sa.pos = s.cell.wrap(u.side[i]);
        const uint32_t idx = uint32_t(s.atoms.size());
        s.atoms.push_back(sa);
        s.bonds.push_back({ts[i].parent < 0 ? base[k] : base[k] + 1 + uint32_t(ts[i].parent), idx, 0});
      }
      if (u.cap) {
        Atom h = at;
        h.element = 1, h.name = "H", h.pos = s.cell.wrap(u.cap_pos);
        s.bonds.push_back({base[k], uint32_t(s.atoms.size()), 0});
        s.atoms.push_back(h);
      }
      R.atom_class.resize(s.atoms.size(), ucls[size_t(ch[k])]);
      R.atom_segment.resize(s.atoms.size(), useg[size_t(ch[k])]);
    }
    // rings and other bonds between the side groups of different units
    for (size_t k = 0; k < ch.size(); ++k)
      for (const ExtraBond& e : T.extra) {
        if (e.sigma != U[size_t(ch[k])].sigma) continue;
        const long k2 = long(k) + long(D) * e.delta;
        if (k2 < 0 || k2 >= long(ch.size())) continue;
        const uint32_t i = base[k] + 1 + uint32_t(e.a), j = e.b < 0 ? base[size_t(k2)] : base[size_t(k2)] + 1 + uint32_t(e.b);
        s.bonds.push_back({i, j, 0});
      }
  }
  int64_t id = 1;
  for (auto& at : s.atoms) at.id = id++;
  s.bonds_from_file = true;

  // --- report: the builder's topology, the structure analysed
  SemicrystallineReport A = analyse_semicrystalline(s);
  R.chains = int(chain.size());
  R.atoms = o.united_atom ? size_t(std::count_if(s.atoms.begin(), s.atoms.end(), [](const Atom& q) { return q.element > 1; })) : s.atoms.size();
  for (const auto& G : seg) {
    const int units = int(G.beads.size());
    if (G.type == 0) {
      ++R.loops;
      (G.a % 2 ? R.loops_top : R.loops_bottom)++;
      R.loop_units.push_back(units);
    } else if (G.type == 1) {
      ++R.bridges;
      R.bridge_units.push_back(units);
    } else {
      ++R.tails;
      (G.a % 2 ? R.tails_top : R.tails_bottom)++;
      R.tail_units.push_back(units);
    }
  }
  for (const auto& ch : chain) R.chain_units.push_back(int(ch.size()));
  R.crystallinity_p2 = A.crystallinity_p2, R.crystallinity_mass = A.crystallinity_mass, R.p2_cut = A.p2_cut;
  R.tilt_deg = A.tilt_deg, R.tilt_sd = A.tilt_sd;
  R.z = A.z, R.rho = A.rho, R.p2 = A.p2;
  R.density = s.density();
  {
    double mass = 0;
    for (const auto& at : s.atoms) if (at.tags & 2u) mass += element(at.element).mass;
    R.amorphous_density = mass / (area * La) * kAmuPerA3;
  }
  R.notes.insert(R.notes.begin(), pc.source);
  if (A.loops != R.loops || A.bridges != R.bridges || A.tails != R.tails)
    R.notes.push_back("the order-parameter analysis of the built cell counts " + std::to_string(A.loops) + " loops, " + std::to_string(A.bridges) + " bridges, " +
                      std::to_string(A.tails) + " tails");
  if (R.strained_closures) R.notes.push_back(std::to_string(R.strained_closures) + " closures placed with strained angles (the relaxation takes them out)");
  if (report) *report = R;
  return s;
}

void relax_semicrystalline(System& s, std::shared_ptr<const ForceField> field, const EnergyOptions& e, double fmax,
                           std::function<bool(const std::string&, double)> progress) {
  const int ct = tag_index(s, "crystal");
  RelaxOptions r;
  r.field = field ? field : std::make_shared<const ForceField>(default_forcefield(s));
  r.energy = e;
  r.ftol = fmax;
  r.max_iterations = 3000;
  r.pushoff = true;
  if (ct >= 0) {
    r.fixed.assign(s.atoms.size(), 0);
    for (size_t i = 0; i < s.atoms.size(); ++i) r.fixed[i] = (s.atoms[i].tags >> ct) & 1u;
  }
  if (progress) r.progress = [&](const RelaxProgress& p) { return progress("push-off (crystal held): " + p.stage, double(p.stage_index) / std::max(1, p.stages)); };
  relax(s, r);
  r.fixed.clear();
  r.pushoff = false;
  r.max_iterations = 2000;
  if (progress) r.progress = [&](const RelaxProgress& p) { return progress("minimise: " + p.stage, double(p.stage_index) / std::max(1, p.stages)); };
  relax(s, r);
}

// ---------------------------------------------------------------------------------------------------------------
// Analysis

SemicrystallineReport analyse_semicrystalline(const System& s, double p2_cut, int min_stem, int bins) {
  SemicrystallineReport R;
  R.p2_cut = p2_cut;
  if (!s.cell.valid()) throw std::invalid_argument("analyse_semicrystalline: no cell");
  const double Lx = norm(s.cell.a), Ly = norm(s.cell.b), Lz = s.cell.c[2];
  R.lx = Lx, R.ly = Ly, R.lz = Lz;
  const double area = norm(cross(s.cell.a, s.cell.b));
  const size_t n = s.atoms.size();
  R.atoms = n;
  // chains: each molecule's backbone, the longest path through its heavy atoms (from the heavy atom farthest from any
  // one to the heavy atom farthest from that); side groups (H, F, CH3, C=O, the far side of a ring) hang off it
  std::vector<std::vector<uint32_t>> hn(n);
  for (const auto& bd : s.bonds)
    if (s.atoms[bd.i].element > 1 && s.atoms[bd.j].element > 1) hn[bd.i].push_back(bd.j), hn[bd.j].push_back(bd.i);
  std::vector<std::vector<uint32_t>> chains;
  std::vector<int> mark(n, -1), from(n, -1);
  auto farthest = [&](uint32_t r, int stamp) {
    std::vector<uint32_t> q{r};
    mark[r] = stamp, from[r] = -1;
    for (size_t h = 0; h < q.size(); ++h)
      for (uint32_t j : hn[q[h]])
        if (mark[j] != stamp) mark[j] = stamp, from[j] = int(q[h]), q.push_back(j);
    return q;
  };
  int stamp = 0;
  std::vector<char> done(n, 0);
  for (uint32_t i = 0; i < n; ++i) {
    if (s.atoms[i].element <= 1 || done[i]) continue;
    const auto comp = farthest(i, stamp++);
    for (uint32_t j : comp) done[j] = 1;
    if (comp.size() < 3) continue;
    const auto far = farthest(comp.back(), stamp++);
    std::vector<uint32_t> ch;
    for (int k = int(far.back()); k >= 0; k = from[size_t(k)]) ch.push_back(uint32_t(k));
    if (ch.size() >= 3) chains.push_back(ch);
  }
  R.chains = int(chains.size());
  // unwrapped positions along each chain
  std::vector<std::vector<Vec3>> up(chains.size());
  for (size_t c = 0; c < chains.size(); ++c) {
    up[c].push_back(s.atoms[chains[c][0]].pos);
    for (size_t k = 1; k < chains[c].size(); ++k) up[c].push_back(up[c][k - 1] + s.cell.minimum_image(s.atoms[chains[c][k]].pos - s.atoms[chains[c][k - 1]].pos));
  }
  // P2 of each backbone atom's chord against dir: the chord spans ±w atoms (w 1 for a zig-zag, whose neighbours' chord
  // lies along the stem; 2 to 5 for a helix: a whole turn of 2_1 or 3_1, about one of 7_2 or 10_3), averaged over ±3 atoms
  auto local_p2 = [&](const Vec3& dir, size_t wmax, std::vector<std::vector<double>>& out) {
    out.assign(chains.size(), {});
    for (size_t c = 0; c < chains.size(); ++c) {
      const size_t m = chains[c].size();
      std::vector<double> raw(m, std::nan(""));
      for (size_t k = 1; k + 1 < m; ++k) {
        const size_t w = std::min<size_t>({wmax, k, m - 1 - k});
        const double ct = dot(unit(up[c][k + w] - up[c][k - w]), dir);
        raw[k] = 1.5 * ct * ct - 0.5;
      }
      out[c].assign(m, 0.0);
      for (size_t k = 0; k < m; ++k) {
        double sum = 0;
        int cnt = 0;
        for (long q = long(k) - 3; q <= long(k) + 3; ++q)
          if (q >= 0 && q < long(m) && !std::isnan(raw[size_t(q)])) sum += raw[size_t(q)], ++cnt;
        out[c][k] = cnt ? sum / cnt : 0;
      }
    }
  };
  // the chords tried: the one that lines up the most atoms near perfectly (P2 > 0.95), the shorter on a near tie (it
  // ends the stems more sharply and smooths the amorphous layer less)
  auto ordered = [&](size_t w, std::vector<std::vector<double>>& p2) {
    // a first direction: the director of all the chords (the order tensor's leading axis; the stems, tilted or not,
    // outweigh the amorphous chords, which spread over every direction)
    double Q[3][3] = {{0}};
    for (size_t c = 0; c < chains.size(); ++c)
      for (size_t k = 1; k + 1 < chains[c].size(); ++k) {
        const size_t ww = std::min<size_t>({w, k, chains[c].size() - 1 - k});
        const Vec3 u = unit(up[c][k + ww] - up[c][k - ww]);
        for (int i = 0; i < 3; ++i)
          for (int j = 0; j < 3; ++j) Q[i][j] += u[i] * u[j];
      }
    Vec3 dir{0.05, 0.05, 1};
    for (int it = 0; it < 100; ++it) {
      const Vec3 v{Q[0][0] * dir[0] + Q[0][1] * dir[1] + Q[0][2] * dir[2], Q[1][0] * dir[0] + Q[1][1] * dir[1] + Q[1][2] * dir[2],
                   Q[2][0] * dir[0] + Q[2][1] * dir[1] + Q[2][2] * dir[2]};
      if (norm(v) < 1e-12) { dir = {0, 0, 1}; break; }
      dir = unit(v);
    }
    if (dir[2] < 0) dir = dir * -1.0;
    local_p2(dir, w, p2);
    Vec3 sum{0, 0, 0};   // the stems' mean direction from the crystalline chords, then P2 again against it
    for (size_t c = 0; c < chains.size(); ++c)
      for (size_t k = 1; k + 1 < chains[c].size(); ++k)
        if (p2[c][k] > p2_cut) {
          const size_t ww = std::min<size_t>({w, k, chains[c].size() - 1 - k});
          Vec3 u = unit(up[c][k + ww] - up[c][k - ww]);
          if (u[2] < 0) u = u * -1.0;
          sum = sum + u;
        }
    if (norm(sum) > 0) local_p2(unit(sum), w, p2);
    size_t sharp = 0;   // atoms in near-perfect order: only the chord that matches the helix's turn lines a stem up exactly
    for (const auto& v : p2) for (double x : v) sharp += x > 0.95;
    return sharp;
  };
  std::vector<std::vector<double>> p2, p2h;
  double best = 1.1 * double(ordered(1, p2));
  for (size_t w : {size_t(2), size_t(3), size_t(4), size_t(5)}) {
    const double n3 = double(ordered(w, p2h));
    if (n3 > best) best = n3, p2.swap(p2h);
  }
  // crystallinity, lamella centre (circular mean of the crystalline atoms' z)
  size_t total = 0, crys = 0;
  double cs = 0, sn = 0;
  for (size_t c = 0; c < chains.size(); ++c)
    for (size_t k = 0; k < chains[c].size(); ++k) {
      ++total;
      if (p2[c][k] > p2_cut) {
        ++crys;
        const double th = 2 * kPi * s.atoms[chains[c][k]].pos[2] / Lz;
        cs += std::cos(th), sn += std::sin(th);
      }
    }
  R.crystallinity_p2 = total ? double(crys) / double(total) : 0;
  const double zc = (crys ? std::atan2(sn, cs) / (2 * kPi) * Lz : 0.5 * Lz);
  auto side = [&](const Vec3& p) {   // +1 above the lamella's middle plane, −1 below
    double d = p[2] - zc;
    d -= Lz * std::round(d / Lz);
    return d >= 0 ? 1 : -1;
  };
  std::vector<double> tilts;
  R.atom_class.assign(n, SemicrystallineReport::kNone);
  R.atom_segment.assign(n, -1);
  int32_t nseg = 0;
  for (size_t c = 0; c < chains.size(); ++c) {
    const size_t m = chains[c].size();
    // stems: crystalline runs of at least min_stem atoms whose ends lie on both sides of the middle plane
    struct Run { size_t lo, hi; int f_lo, f_hi; };
    std::vector<Run> runs;
    for (size_t k = 0; k < m;) {
      if (p2[c][k] <= p2_cut) { ++k; continue; }
      size_t e = k;
      while (e + 1 < m && p2[c][e + 1] > p2_cut) ++e;
      if (e - k + 1 >= size_t(min_stem)) {
        const Vec3 span = up[c][e] - up[c][k];
        const int fl = side(s.atoms[chains[c][k]].pos), fh = side(s.atoms[chains[c][e]].pos);
        // the run crosses the middle plane: its ends on both sides, and it spans more than half the lamella
        if (fl != fh) {
          runs.push_back({k, e, fl, fh});
          tilts.push_back(std::acos(std::min(1.0, std::fabs(span[2]) / norm(span))) / kDeg);
        }
      }
      k = e + 1;
    }
    R.stems += int(runs.size());
    if (runs.empty()) continue;
    auto mark = [&](size_t lo, size_t hi, int8_t k) {   // chain positions lo..hi inclusive, one new segment
      for (size_t q = lo; q <= hi && q < m; ++q) R.atom_class[chains[c][q]] = k, R.atom_segment[chains[c][q]] = nseg;
      ++nseg;
    };
    if (runs.front().lo > 0) mark(0, runs.front().lo - 1, SemicrystallineReport::kTail);
    if (runs.back().hi + 1 < m) mark(runs.back().hi + 1, m - 1, SemicrystallineReport::kTail);
    for (size_t r = 0; r < runs.size(); ++r) {
      mark(runs[r].lo, runs[r].hi, SemicrystallineReport::kStem);
      if (r + 1 < runs.size() && runs[r + 1].lo > runs[r].hi + 1)
        mark(runs[r].hi + 1, runs[r + 1].lo - 1, runs[r].f_hi == runs[r + 1].f_lo ? SemicrystallineReport::kLoop : SemicrystallineReport::kBridge);
    }
    R.tails += 2;
    R.tail_units.push_back(int(runs.front().lo));
    R.tail_units.push_back(int(m - 1 - runs.back().hi));
    (runs.front().f_lo > 0 ? R.tails_top : R.tails_bottom)++;
    (runs.back().f_hi > 0 ? R.tails_top : R.tails_bottom)++;
    for (size_t r = 0; r + 1 < runs.size(); ++r) {
      const int units = int(runs[r + 1].lo - runs[r].hi - 1);
      if (runs[r].f_hi == runs[r + 1].f_lo) {
        ++R.loops;
        (runs[r].f_hi > 0 ? R.loops_top : R.loops_bottom)++;
        R.loop_units.push_back(units);
      } else {
        ++R.bridges;
        R.bridge_units.push_back(units);
      }
    }
    R.chain_units.push_back(int(m));
  }
  // every side atom as the backbone atom it hangs from (breadth first from the backbones, never through another one)
  {
    std::vector<std::vector<uint32_t>> nb(n);
    for (const auto& bd : s.bonds) nb[bd.i].push_back(bd.j), nb[bd.j].push_back(bd.i);
    std::vector<char> reached(n, 0);
    std::vector<uint32_t> q;
    for (const auto& ch : chains)
      for (uint32_t i : ch) reached[i] = 1, q.push_back(i);
    for (size_t h = 0; h < q.size(); ++h)
      for (uint32_t j : nb[q[h]])
        if (!reached[j]) {
          reached[j] = 1;
          R.atom_class[j] = R.atom_class[q[h]], R.atom_segment[j] = R.atom_segment[q[h]];
          q.push_back(j);
        }
  }
  if (!tilts.empty()) {
    double mt = 0, q = 0;
    for (double t : tilts) mt += t;
    mt /= double(tilts.size());
    for (double t : tilts) q += (t - mt) * (t - mt);
    R.tilt_deg = mt;
    R.tilt_sd = tilts.size() > 1 ? std::sqrt(q / double(tilts.size() - 1)) : 0;
  }
  // profiles along z
  bins = std::max(4, bins);
  const double dz = Lz / bins;
  R.z.resize(size_t(bins)), R.rho.assign(size_t(bins), 0.0), R.p2.assign(size_t(bins), 0.0);
  std::vector<int> pc(size_t(bins), 0);
  auto zbin = [&](double z) {
    double f = (z - s.cell.origin[2]) / Lz;
    f -= std::floor(f);
    return std::min(bins - 1, int(f * bins));
  };
  double mass_all = 0, mass_crys = 0;
  const int ct = tag_index(s, "crystal");
  for (const auto& at : s.atoms) {
    const double m = s.mass_of(at);
    R.rho[size_t(zbin(at.pos[2]))] += m;
    mass_all += m;
    if (ct >= 0 && ((at.tags >> ct) & 1u)) mass_crys += m;
  }
  for (size_t c = 0; c < chains.size(); ++c)
    for (size_t k = 0; k < chains[c].size(); ++k) {
      const int bz = zbin(s.atoms[chains[c][k]].pos[2]);
      R.p2[size_t(bz)] += p2[c][k];
      pc[size_t(bz)]++;
    }
  for (int k = 0; k < bins; ++k) {
    R.z[size_t(k)] = s.cell.origin[2] + (k + 0.5) * dz;
    R.rho[size_t(k)] *= kAmuPerA3 / (area * dz);
    R.p2[size_t(k)] = pc[size_t(k)] ? R.p2[size_t(k)] / pc[size_t(k)] : 0;
  }
  R.crystallinity_mass = ct >= 0 && mass_all > 0 ? mass_crys / mass_all : R.crystallinity_p2;
  R.density = s.density();
  {   // the amorphous layer's density: the bins whose chords are disordered (mean P2 below 0.2)
    double sr = 0;
    int nr = 0, nc = 0;
    for (int k = 0; k < bins; ++k) {
      if (R.p2[size_t(k)] < 0.2 && pc[size_t(k)]) sr += R.rho[size_t(k)], ++nr;
      if (R.p2[size_t(k)] > 0.5) ++nc;
    }
    R.amorphous_density = nr ? sr / nr : 0;
    R.crystal_thickness = nc * dz;
    R.amorphous_thickness = Lz - R.crystal_thickness;
  }
  return R;
}

std::string semicrystalline_json(const SemicrystallineReport& r) {
  Json j = Json::object();
  auto ints = [](const std::vector<int>& v) { Json a = Json::array(); for (int x : v) a.push_back(Json(double(x))); return a; };
  auto dbl = [](const std::vector<double>& v) { Json a = Json::array(); for (double x : v) a.push_back(Json(x)); return a; };
  j["format"] = Json("caps-semicrystalline/1");
  if (!r.polymer.empty()) j["polymer"] = Json(r.polymer);
  j["repeat_period"] = Json(double(r.repeat_period)), j["polar"] = Json(r.polar);
  j["cell_tilt_deg"] = Json(r.cell_tilt_deg), j["crystal_density"] = Json(r.crystal_density);
  j["nx"] = Json(double(r.nx)), j["ny"] = Json(double(r.ny)), j["nc"] = Json(double(r.nc));
  j["atoms"] = Json(double(r.atoms)), j["stems"] = Json(double(r.stems)), j["chains"] = Json(double(r.chains));
  j["loops"] = Json(double(r.loops)), j["bridges"] = Json(double(r.bridges)), j["tails"] = Json(double(r.tails));
  j["loops_top"] = Json(double(r.loops_top)), j["loops_bottom"] = Json(double(r.loops_bottom));
  j["tails_top"] = Json(double(r.tails_top)), j["tails_bottom"] = Json(double(r.tails_bottom));
  const double segs = r.loops + r.bridges + r.tails;
  j["fraction_loops"] = Json(segs ? r.loops / segs : 0.0), j["fraction_bridges"] = Json(segs ? r.bridges / segs : 0.0), j["fraction_tails"] = Json(segs ? r.tails / segs : 0.0);
  j["bridge_units"] = ints(r.bridge_units), j["loop_units"] = ints(r.loop_units), j["tail_units"] = ints(r.tail_units), j["chain_units"] = ints(r.chain_units);
  j["cell"] = dbl({r.lx, r.ly, r.lz});
  j["crystal_thickness"] = Json(r.crystal_thickness), j["amorphous_thickness"] = Json(r.amorphous_thickness);
  j["density"] = Json(r.density), j["amorphous_density"] = Json(r.amorphous_density);
  j["crystallinity_p2"] = Json(r.crystallinity_p2), j["crystallinity_mass"] = Json(r.crystallinity_mass), j["p2_cut"] = Json(r.p2_cut);
  j["tilt_deg"] = Json(r.tilt_deg), j["tilt_sd"] = Json(r.tilt_sd);
  j["z"] = dbl(r.z), j["rho"] = dbl(r.rho), j["p2"] = dbl(r.p2);
  j["worst_contact"] = Json(r.worst_contact), j["eased_steps"] = Json(double(r.eased_steps)), j["strained_closures"] = Json(double(r.strained_closures));
  Json notes = Json::array();
  for (const auto& x : r.notes) notes.push_back(Json(x));
  j["notes"] = notes;
  return j.dump(1);
}

std::string semicrystalline_summary(const SemicrystallineReport& r) {
  char b[1024];
  const double segs = std::max(1, r.loops + r.bridges + r.tails);
  auto mean = [](const std::vector<int>& v) { return v.empty() ? 0.0 : std::accumulate(v.begin(), v.end(), 0.0) / double(v.size()); };
  std::string out;
  std::snprintf(b, sizeof b, "cell %.2f × %.2f × %.2f Å · %d × %d × %d unit cells · lamella %.1f Å · amorphous %.1f Å\n", r.lx, r.ly, r.lz, r.nx, r.ny, r.nc,
                r.crystal_thickness, r.amorphous_thickness);
  out += b;
  std::snprintf(b, sizeof b, "%zu atoms · %d stems · %d chains (%.0f backbone atoms on average)\n", r.atoms, r.stems, r.chains, mean(r.chain_units));
  out += b;
  std::snprintf(b, sizeof b, "segments: %d loops (%.0f %%; top %d, bottom %d) · %d bridges (%.0f %%) · %d tails (%.0f %%)\n", r.loops, 100 * r.loops / segs,
                r.loops_top, r.loops_bottom, r.bridges, 100 * r.bridges / segs, r.tails, 100 * r.tails / segs);
  out += b;
  std::snprintf(b, sizeof b, "mean length (backbone atoms): loops %.1f · bridges %.1f · tails %.1f\n", mean(r.loop_units), mean(r.bridge_units), mean(r.tail_units));
  out += b;
  std::snprintf(b, sizeof b, "density %.3f g/cm³ (amorphous layer %.3f) · crystallinity %.3f by P2 > %.2f, %.3f by mass · stem tilt %.1f ± %.1f°\n", r.density,
                r.amorphous_density, r.crystallinity_p2, r.p2_cut, r.crystallinity_mass, r.tilt_deg, r.tilt_sd);
  out += b;
  return out;
}

std::string semicrystalline_lammps_run(const SemicrystallineEquilibration& q) {
  std::ostringstream o;
  const long n1 = std::lround(q.restrained_ps * 1000 / q.dt), n2 = std::lround(q.free_ps * 1000 / q.dt);
  o << "\n# ---- CAPS semicrystalline equilibration: minimise, NVT with the crystal tethered, then NPT (each axis on its own)\n"
    << "min_style cg\nminimize 1.0e-4 1.0e-6 5000 50000\nreset_timestep 0\n"
    << "timestep " << q.dt << "\n"
    << "velocity all create " << q.temperature << " 4928459 dist gaussian mom yes rot yes\n"
    << "thermo_style custom step temp press pe ke density lx ly lz\nthermo 1000\n"
    << "dump traj all custom 10000 semicrystalline.lammpstrj id mol type xu yu zu\n"
    << "fix tether crystal spring/self " << q.spring << "\n"
    << "fix stage1 all nvt temp " << q.temperature << " " << q.temperature << " 100.0\n"
    << "run " << n1 << "\n"
    << "unfix tether\nunfix stage1\n"
    << "fix stage2 all npt temp " << q.temperature << " " << q.temperature << " 100.0 x " << q.pressure << " " << q.pressure << " 1000.0 y " << q.pressure
    << " " << q.pressure << " 1000.0 z " << q.pressure << " " << q.pressure << " 1000.0 couple none\n"
    << "variable rho equal density\nfix rho all ave/time 100 100 10000 v_rho c_thermo_press file density.dat\n"
    << "run " << n2 << "\n"
    << "write_data semicrystalline_equilibrated.data pair ij\n";
  return o.str();
}

}  // namespace caps
