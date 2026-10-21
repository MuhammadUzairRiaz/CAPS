// CAPS semicrystalline cells (see semicrystal.hpp).
#include "caps/semicrystal.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <numeric>
#include <random>
#include <set>
#include <sstream>
#include <stdexcept>

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
constexpr double kCC = 1.53;                     // Å, C-C (Bunn 1939: 1.53)
constexpr double kTheta = 112.0 * kDeg;          // C-C-C (Bunn 1939: 112°)
constexpr double kCH = 1.09;                     // Å
constexpr double kHalfHCH = 0.5 * 109.4712206 * kDeg;
const double kReach = kCC * std::sin(kTheta / 2);   // projection of one bond on an all-trans chain's axis, 1.268 Å
constexpr double kMassCH2 = 14.027;              // g/mol
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

// Orthorhombic periodic box: minimum image and a spatial hash of points (positions kept unwrapped by the caller).
struct Box {
  double L[3];
  Vec3 mi(Vec3 d) const {
    for (int k = 0; k < 3; ++k) d[k] -= L[k] * std::round(d[k] / L[k]);
    return d;
  }
};

struct Hash {
  Box box;
  int n[3];
  double w[3];
  std::vector<std::vector<uint32_t>> bins;
  std::vector<Vec3> pos;
  Hash(const Box& b, double cell) : box(b) {
    for (int k = 0; k < 3; ++k) {
      n[k] = std::max(1, int(b.L[k] / cell));
      w[k] = b.L[k] / n[k];
    }
    bins.resize(size_t(n[0]) * n[1] * n[2]);
  }
  int bin(double x, int k) const {
    double f = x / box.L[k];
    f -= std::floor(f);
    return std::min(n[k] - 1, int(f * n[k]));
  }
  uint32_t add(const Vec3& p) {
    const uint32_t id = uint32_t(pos.size());
    pos.push_back(p);
    bins[(size_t(bin(p[0], 0)) * n[1] + bin(p[1], 1)) * n[2] + bin(p[2], 2)].push_back(id);
    return id;
  }
  // f(id, distance) for every point within the bins around p (the caller checks the distance)
  template <class F>
  void near(const Vec3& p, F f) const {
    const int b[3] = {bin(p[0], 0), bin(p[1], 1), bin(p[2], 2)};
    std::vector<int> off[3];
    for (int k = 0; k < 3; ++k) {
      if (n[k] < 3) for (int q = 0; q < n[k]; ++q) off[k].push_back(q);
      else for (int d = -1; d <= 1; ++d) off[k].push_back(((b[k] + d) % n[k] + n[k]) % n[k]);
    }
    for (int x : off[0])
      for (int y : off[1])
        for (int z : off[2])
          for (uint32_t id : bins[(size_t(x) * n[1] + y) * n[2] + z]) f(id, norm(box.mi(pos[id] - p)));
  }
};

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

PolymerCrystalCell polymer_crystal_cell(const std::string& id) {
  PolymerCrystalCell pc;
  std::string low = id;
  for (auto& ch : low) ch = char(std::tolower((unsigned char)ch));
  if (low == "pe" || low == "polyethylene" || low == "pe-orthorhombic") {
    // Bunn 1939, Pnam: C at (x, y, 1/4) etc. with x = 0.038, y = 0.065 (p. 491)
    pc.id = "pe";
    pc.name = "polyethylene, orthorhombic (Pnam)";
    pc.source = "C. W. Bunn, Trans. Faraday Soc. 35, 482 (1939): a 7.40, b 4.93, c 2.534 Å; C-C 1.53 Å, C-C-C 112°; hydrogens placed by CAPS (C-H 1.09 Å, H-C-H 109.47°)";
    System& s = pc.cell;
    s.title = pc.name;
    s.cell.a = {7.40, 0, 0};
    s.cell.b = {0, 4.93, 0};
    s.cell.c = {0, 0, 2.534};
    const double fx = 0.038, fy = 0.065;
    const Vec3 fc[4] = {{fx, 1 - fy, 0.25}, {1 - fx, fy, 0.75}, {0.5 - fx, 0.5 - fy, 0.75}, {0.5 + fx, 0.5 + fy, 0.25}};
    for (const auto& f : fc) {
      Atom a;
      a.element = 6, a.name = "C";
      a.pos = s.cell.to_cartesian(f);
      s.atoms.push_back(a);
    }
    // hydrogens: each carbon's chain neighbours are the other carbon of its chain at ±c/2 (minimum image along x, y)
    const Box box{{7.40, 4.93, 2.534}};
    for (int i = 0; i < 4; ++i) {
      const int j = i ^ 1;   // 0-1 the corner chain, 2-3 the centre chain
      const Vec3 c = s.atoms[i].pos;
      Vec3 d = box.mi(s.atoms[j].pos - c);
      d[2] = 0.5 * 2.534 * (d[2] >= 0 ? 1 : -1);
      const Vec3 a = c + d, b = c + Vec3{d[0], d[1], -d[2]};
      const Vec3 u = unit(unit(c - a) + unit(c - b)), n = unit(cross(a - c, b - c));
      for (double sg : {1.0, -1.0}) {
        Atom h;
        h.element = 1, h.name = "H";
        h.pos = s.cell.wrap(c + (u * std::cos(kHalfHCH) + n * (sg * std::sin(kHalfHCH))) * kCH);
        s.atoms.push_back(h);
      }
    }
    return pc;
  }
  // a file: a CIF or any structure with a cell (the polymer crystal library's entries)
  pc.id = id;
  const bool cif = low.size() > 4 && low.substr(low.size() - 4) == ".cif";
  pc.cell = cif ? read_cif(id) : open_file(id).frame(0);
  if (!pc.cell.cell.valid()) throw std::invalid_argument(id + ": no unit cell");
  pc.name = pc.cell.title.empty() ? id : pc.cell.title;
  pc.source = id;
  return pc;
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
          if (!repeats || std::abs(t[2]) < std::abs(trans[2]) || (t[2] != 0 && trans[2] == 0)) trans = t, repeats = true;
        }
      }
    }
    CrystalStem st;
    st.atoms = members;
    st.formula = formula_of(s, members);
    if (!repeats || trans[0] != 0 || trans[1] != 0 || std::abs(trans[2]) != 1) {
      out.push_back(st);   // a molecule, or a chain not along c: atoms only, no backbone
      continue;
    }
    st.axis = unit(s.cell.c);
    // backbone: heavy atoms that lie on the path from one cell's chain to the next (their removal splits the chain);
    // found as the heavy atoms with two or more heavy neighbours whose shortest path to their own +c image passes them
    std::vector<uint32_t> heavy;
    for (uint32_t i : members) if (s.atoms[i].element > 1) heavy.push_back(i);
    // BFS over (atom, image) restricted to heavy atoms from the lowest heavy atom to its +c image
    auto unwrapped = [&](uint32_t i, const std::array<int, 3>& im) {
      return s.atoms[i].pos + s.cell.a * im[0] + s.cell.b * im[1] + s.cell.c * im[2];
    };
    uint32_t start = heavy.empty() ? members[0] : heavy[0];
    for (uint32_t i : heavy) if (unwrapped(i, img[i])[2] < unwrapped(start, img[start])[2]) start = i;
    using Key = std::pair<uint32_t, std::array<int, 3>>;
    std::map<Key, Key> prev;
    std::vector<Key> q{{start, img[start]}};
    prev[q[0]] = q[0];
    const Key goal{start, {img[start][0], img[start][1], img[start][2] + 1}};
    bool found = false;
    for (size_t h = 0; h < q.size() && !found; ++h) {
      const Key cur = q[h];
      for (const auto& l : nb[cur.first]) {
        if (s.atoms[l.j].element <= 1) continue;
        const Key nx{l.j, {cur.second[0] + l.sh[0], cur.second[1] + l.sh[1], cur.second[2] + l.sh[2]}};
        if (prev.count(nx)) continue;
        if (std::abs(nx.second[2] - img[start][2]) > 2) continue;
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
      for (const auto& k : path) st.backbone.push_back(k.first);
      // the backbone in order along +c, starting at the lowest atom (already: start is the lowest)
      // helical: the backbone's projections on the ab plane span more than a line
      if (st.backbone.size() >= 3) {
        double area = 0;
        const Vec3 p0 = unwrapped(path[0].first, path[0].second);
        for (size_t k = 1; k + 1 < path.size(); ++k) {
          Vec3 u = unwrapped(path[k].first, path[k].second) - p0, v = unwrapped(path[k + 1].first, path[k + 1].second) - p0;
          u[2] = v[2] = 0;
          area += std::fabs(cross(u, v)[2]);
        }
        st.helical = area > 0.5;
      }
    }
    out.push_back(st);
  }
  return out;
}

// ---------------------------------------------------------------------------------------------------------------
// Builder

namespace {

struct Seg {
  int type = 0;               // 0 loop, 1 bridge, 2 tail
  int a = -1, b = -1;         // stem ends (2 stem + 1 top, 2 stem bottom); b = −1 for a tail
  int n = 0;                  // backbone atoms (CH2) to grow
  std::vector<uint32_t> beads;
};

}  // namespace

System build_semicrystalline(const SemicrystallineOptions& o, SemicrystallineReport* report) {
  auto cancelled = [&](const std::string& what, double f) { return o.progress && !o.progress(what, f); };
  SemicrystallineReport R;
  const PolymerCrystalCell pc = polymer_crystal_cell(o.polymer);
  const System& uc = pc.cell;
  const auto stems_uc = find_stems(uc);
  // polymethylene stems only for now: the amorphous segments are grown as CH2 chains
  std::vector<const CrystalStem*> chains;
  for (const auto& st : stems_uc) {
    if (st.backbone.empty()) throw std::invalid_argument(pc.name + ": a molecule or chain not along c (" + st.formula + ")");
    int c = 0, h = 0, other = 0;
    for (uint32_t i : st.atoms) (uc.atoms[i].element == 6 ? c : uc.atoms[i].element == 1 ? h : other)++;
    if (other || h != 2 * c || st.backbone.size() != size_t(c))
      throw std::invalid_argument(pc.name + ": stems of " + st.formula + " (" + (st.helical ? "helical" : "planar") +
                                  "); the amorphous layers are grown for polymethylene (CH2) stems only so far");
    chains.push_back(&st);
  }
  if (chains.empty()) throw std::invalid_argument(pc.name + ": no stems along c");
  const double a = norm(uc.cell.a), b = norm(uc.cell.b), c = norm(uc.cell.c);
  if (std::fabs(uc.cell.c[0]) > 1e-6 || std::fabs(uc.cell.c[1]) > 1e-6 || std::fabs(uc.cell.a[1]) > 1e-6 || std::fabs(uc.cell.b[0]) > 1e-6)
    throw std::invalid_argument(pc.name + ": the lamellar builder needs an orthorhombic cell with c along z (no tilt yet)");
  if (o.nx < 1 || o.ny < 1) throw std::invalid_argument("lateral repeats must be at least 1");
  const int nc = o.nc > 0 ? o.nc : std::max(2, int(std::lround(o.crystal_nm * 10 / c)));
  const double La = o.amorphous_nm * 10;
  if (La < 5) throw std::invalid_argument("the amorphous layer must be at least 0.5 nm");
  const double Lx = a * o.nx, Ly = b * o.ny, Lc = nc * c, Lz = Lc + La;
  const Box box{{Lx, Ly, Lz}};
  const double z0 = 0.5 * La;   // the crystal from z0 to z0 + Lc; the amorphous layer across z = 0
  R.nx = o.nx, R.ny = o.ny, R.nc = nc;
  R.lx = Lx, R.ly = Ly, R.lz = Lz;
  R.crystal_thickness = Lc, R.amorphous_thickness = La;

  // --- the crystal: stems bottom to top (backbone carbons; hydrogens are placed at the end for every carbon)
  std::vector<Vec3> X;                 // backbone carbons, unwrapped as grown
  std::vector<std::vector<uint32_t>> stem;   // carbon ids bottom → top
  for (int ix = 0; ix < o.nx; ++ix)
    for (int iy = 0; iy < o.ny; ++iy)
      for (const CrystalStem* st : chains) {
        std::vector<uint32_t> ids;
        for (int k = 0; k < nc; ++k)
          for (uint32_t i : st->backbone) {
            Vec3 p = uc.atoms[i].pos;
            // the backbone atom in the image that keeps the chain continuous (the stem's atoms near its first one)
            const Vec3 ref = uc.atoms[st->backbone[0]].pos;
            Vec3 d = p - ref;
            d[0] -= a * std::round(d[0] / a), d[1] -= b * std::round(d[1] / b);
            p = ref + d;
            p = p + uc.cell.a * ix + uc.cell.b * iy;
            p[2] = uc.atoms[i].pos[2] + z0 + k * c;
            ids.push_back(uint32_t(X.size()));
            X.push_back(p);
          }
        // order along +z
        std::sort(ids.begin(), ids.end(), [&](uint32_t u, uint32_t v) { return X[u][2] < X[v][2]; });
        stem.push_back(ids);
      }
  const int N = int(stem.size());
  const size_t ncrys = X.size();
  R.stems = N;
  double ztop = -1e300, zbot = 1e300;
  for (const auto& st : stem) ztop = std::max(ztop, X[st.back()][2]), zbot = std::min(zbot, X[st.front()][2]);

  // --- the amorphous layer: CH2 units at the target density
  const int M = int(std::lround(o.density * Lx * Ly * La / kAmuPerA3 / kMassCH2));
  {
    const size_t carbons = ncrys + size_t(M);
    const size_t atoms = o.united_atom ? carbons : 3 * carbons + 2;   // + the chain ends' third hydrogens (≈ tails)
    if (atoms > o.max_atoms)
      throw std::invalid_argument("the cell would hold about " + std::to_string(atoms) + " atoms (" + std::to_string(ncrys) + " crystal and " + std::to_string(M) +
                                  " amorphous CH2" + (o.united_atom ? "" : ", all atom") + "), above the limit of " + std::to_string(o.max_atoms) +
                                  ": fewer lateral repeats, thinner layers" + (o.united_atom ? "" : ", or a united-atom force field (oplsua2024)"));
  }

  std::mt19937_64 rng(o.seed);
  // --- segment counts: S segments with fractions f use 2N ends: S (2 − f_tail) = 2N
  double fl = std::max(0.0, o.loops), fb = std::max(0.0, o.bridges), ft = std::max(0.0, o.tails);
  const double fs = fl + fb + ft;
  if (fs <= 0) throw std::invalid_argument("the loop, bridge and tail fractions are all zero");
  fl /= fs, fb /= fs, ft /= fs;
  const double S = 2.0 * N / (2.0 - ft);
  int nt = 2 * int(std::lround(ft * S / 2));
  int nb = int(std::lround(fb * S));
  nt = std::max(nt, 2);   // chains need ends
  nb = std::min(nb, N - 1);
  int ntT = nt / 2;
  if ((N - nb - ntT) % 2 != 0) ++ntT;
  int ntB = nt - ntT;
  if ((N - nb - ntB) % 2 != 0) { --ntT, ++ntB; }   // ntT, ntB both ≡ N − nb (mod 2)
  if (ntT < 0 || ntB < 0 || (N - nb - ntT) < 0 || (N - nb - ntB) < 0) throw std::invalid_argument("these fractions leave no valid pairing of the stem ends");
  if (ntT == 0 && ntB == 0) throw std::invalid_argument("chains need tails");

  // --- pair the stem ends: end e = 2 stem + (1 top | 0 bottom)
  auto lateral = [&](int e1, int e2) {
    const Vec3 p = X[e1 % 2 ? stem[e1 / 2].back() : stem[e1 / 2].front()], q = X[e2 % 2 ? stem[e2 / 2].back() : stem[e2 / 2].front()];
    Vec3 d = box.mi(q - p);
    d[2] = 0;
    return norm(d);
  };
  std::vector<Seg> seg;
  std::vector<int> seg_of(2 * N, -1);
  std::vector<int> topE(N), botE(N);
  for (int i = 0; i < N; ++i) topE[i] = 2 * i + 1, botE[i] = 2 * i;
  shuffle(topE.begin(), topE.end(), rng);
  shuffle(botE.begin(), botE.end(), rng);
  std::vector<char> used(2 * N, 0);
  Normal<double> gauss;
  const double sig = std::max(5.0, 0.4 * La);
  for (int k = 0; k < nb; ++k) {
    const int e = topE[k];
    const Vec3 p = X[stem[e / 2].back()];
    const Vec3 target{p[0] + sig * gauss(rng), p[1] + sig * gauss(rng), 0};
    int best = -1;
    double bd = 1e300;
    for (int f : botE) {
      if (used[f]) continue;
      Vec3 d = box.mi(X[stem[f / 2].front()] - target);
      d[2] = 0;
      if (norm(d) < bd) bd = norm(d), best = f;
    }
    used[e] = used[best] = 1;
    seg_of[e] = seg_of[best] = int(seg.size());
    seg.push_back({1, e, best, 0, {}});
  }
  auto face = [&](std::vector<int>& ends, int ntails) {
    std::vector<int> rest;
    for (int e : ends) if (!used[e]) rest.push_back(e);
    for (int k = 0; k < ntails; ++k) {
      used[rest[k]] = 1;
      seg_of[rest[k]] = int(seg.size());
      seg.push_back({2, rest[k], -1, 0, {}});
    }
    std::vector<int> loopers(rest.begin() + ntails, rest.end());
    UniformReal<double> u01(0.0, 1.0);
    for (size_t k = 0; k < loopers.size(); ++k) {
      const int e = loopers[k];
      if (used[e]) continue;
      int partner = -1;
      std::vector<int> free;
      for (size_t q = k + 1; q < loopers.size(); ++q) if (!used[loopers[q]]) free.push_back(loopers[q]);
      if (free.empty()) throw std::logic_error("unpaired loop end");
      if (o.reentry == "random" || u01(rng) >= 0.7) {
        partner = free[UniformInt<size_t>(0, free.size() - 1)(rng)];
      } else {
        double bd = 1e300;
        for (int f : free) if (lateral(e, f) < bd) bd = lateral(e, f), partner = f;
      }
      used[e] = used[partner] = 1;
      seg_of[e] = seg_of[partner] = int(seg.size());
      seg.push_back({0, e, partner, 0, {}});
    }
  };
  if (o.reentry != "adjacent" && o.reentry != "random") throw std::invalid_argument("reentry: adjacent or random");
  face(topE, ntT);
  face(botE, ntB);

  // --- no chain may close on itself: a cycle (stems and segments with no tail) is joined to a chain with tails by
  // exchanging the partners of one of its segments and a segment of the same kind outside it
  auto partner_of = [&](int e) {
    const Seg& g = seg[seg_of[e]];
    return g.type == 2 ? -1 : (g.a == e ? g.b : g.a);
  };
  int merges = 0, broken = 0;
  for (int pass = 0; pass < 4 * N; ++pass) {
    // components over stems
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
    for (int k = 0; k < ncomp; ++k) if (!has_tail[k]) { cyc = k; break; }
    if (cyc < 0) break;
    // a segment x in the cycle and the nearest segment y of the same kind (same face for loops) in a chain with tails
    int bx = -1, by = -1;
    double bd = 1e300;
    for (size_t x = 0; x < seg.size(); ++x) {
      if (seg[x].type == 2 || comp[seg[x].a / 2] != cyc) continue;
      for (size_t y = 0; y < seg.size(); ++y) {
        if (seg[y].type != seg[x].type || !has_tail[comp[seg[y].a / 2]]) continue;
        if (seg[x].type == 0 && (seg[x].a % 2) != (seg[y].a % 2)) continue;
        const double d = lateral(seg[x].a, seg[y].a);
        if (d < bd) bd = d, bx = int(x), by = int(y);
      }
    }
    if (bx < 0) {   // nothing to exchange with: the cycle opens into two tails
      for (size_t x = 0; x < seg.size(); ++x)
        if (seg[x].type != 2 && comp[seg[x].a / 2] == cyc) {
          const int e1 = seg[x].a, e2 = seg[x].b;
          seg[x] = {2, e1, -1, 0, {}};
          seg_of[e2] = int(seg.size());
          seg.push_back({2, e2, -1, 0, {}});
          ++broken;
          break;
        }
      continue;
    }
    Seg &X1 = seg[bx], &Y1 = seg[by];
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

  // --- target end positions (unwrapped: the bottom face is reached in the image above, z + Lz)
  auto end_atoms = [&](int e, Vec3& s1, Vec3& s2, Vec3& s3) {
    const auto& st = stem[e / 2];
    const Vec3 up{0, 0, e % 2 ? 0.0 : Lz};
    if (e % 2) s1 = X[st[st.size() - 1]], s2 = X[st[st.size() - 2]], s3 = X[st[st.size() - 3]];
    else s1 = X[st[0]] + up, s2 = X[st[1]] + up, s3 = X[st[2]] + up;
  };
  auto image_near = [&](Vec3 p, const Vec3& ref) {
    p[0] += Lx * std::round((ref[0] - p[0]) / Lx);
    p[1] += Ly * std::round((ref[1] - p[1]) / Ly);
    return p;
  };
  // grow every loop and bridge from its top end (bridges: top → bottom; bottom loops start at their first end)
  std::vector<int> nmin(seg.size(), 1);
  for (size_t g = 0; g < seg.size(); ++g) {
    Seg& G = seg[g];
    if (G.type == 1 && G.a % 2 == 0) std::swap(G.a, G.b);   // bridges grow from the top face
    if (G.type == 2) { nmin[g] = 1; continue; }
    Vec3 s1, s2, s3, e1, e2, e3;
    end_atoms(G.a, s1, s2, s3);
    end_atoms(G.b, e1, e2, e3);
    e1 = image_near(e1, s1);
    const double chord = norm(e1 - s1);
    nmin[g] = G.type == 1 ? int(std::ceil(chord / (0.72 * kReach))) : std::max(4, int(std::ceil(chord / (0.6 * kReach))) + 2);
  }
  const int need = std::accumulate(nmin.begin(), nmin.end(), 0);
  if (need > M)
    throw std::invalid_argument("the amorphous layer holds " + std::to_string(M) + " CH2 at " + std::to_string(o.density).substr(0, 5) +
                                " g/cm³ but the loops, bridges and tails need at least " + std::to_string(need) +
                                ": a thicker amorphous layer, fewer bridges, or loops that re-enter nearby");
  {   // the rest shared out: exponential lengths (loops and bridges weight 1, tails 1.5: their free ends fill the middle of the layer)
    std::vector<double> w(seg.size());
    UniformReal<double> u01(1e-12, 1.0);
    double sw = 0;
    for (size_t g = 0; g < seg.size(); ++g) sw += w[g] = -std::log(u01(rng)) * (seg[g].type == 2 ? 1.5 : 1.0);
    int left = M - need;
    std::vector<int> extra(seg.size());
    int given = 0;
    for (size_t g = 0; g < seg.size(); ++g) given += extra[g] = int(std::floor(left * w[g] / sw));
    for (int k = given; k < left; ++k) extra[UniformInt<size_t>(0, seg.size() - 1)(rng)]++;
    for (size_t g = 0; g < seg.size(); ++g) seg[g].n = nmin[g] + extra[g];
  }

  // --- growth: one bead per segment per round, in random order
  Hash hash(box, std::max(o.contact, 3.0));
  for (const Vec3& p : X) hash.add(p);
  const double zlo = ztop + 0.3, zhi = zbot + Lz - 0.3;
  struct Grow {
    std::vector<Vec3> path;     // the last atoms (stem end atoms first)
    std::vector<uint32_t> ids;  // their ids (stem atoms too) for exclusions
    Vec3 e1, e2, e3;            // closure target and its stem neighbours (unwrapped image)
    uint32_t te1 = 0, te2 = 0, te3 = 0;
    int placed = 0;
  };
  std::vector<Grow> gs(seg.size());
  for (size_t g = 0; g < seg.size(); ++g) {
    Vec3 s1, s2, s3;
    end_atoms(seg[g].a, s1, s2, s3);
    const auto& st = stem[seg[g].a / 2];
    gs[g].path = {s3, s2, s1};
    gs[g].ids = seg[g].a % 2 ? std::vector<uint32_t>{st[st.size() - 3], st[st.size() - 2], st.back()} : std::vector<uint32_t>{st[2], st[1], st[0]};
    if (seg[g].type != 2) {
      end_atoms(seg[g].b, gs[g].e1, gs[g].e2, gs[g].e3);
      const Vec3 sh = image_near(gs[g].e1, s1) - gs[g].e1;
      gs[g].e1 = gs[g].e1 + sh, gs[g].e2 = gs[g].e2 + sh, gs[g].e3 = gs[g].e3 + sh;
      const auto& sb = stem[seg[g].b / 2];
      if (seg[g].b % 2) gs[g].te1 = sb.back(), gs[g].te2 = sb[sb.size() - 2], gs[g].te3 = sb[sb.size() - 3];
      else gs[g].te1 = sb[0], gs[g].te2 = sb[1], gs[g].te3 = sb[2];
    }
  }
  double worst = 1e300;
  int eased = 0, strained = 0;
  UniformReal<double> u01(0.0, 1.0);
  const double kC = 6.7;   // characteristic ratio of polyethylene (Flory), for the closure bias
  std::vector<char> dead(X.size(), 0);
  // the layer filled evenly: trials weighted down where the beads grown so far already exceed their share (2 Å slices)
  const int nzb = std::max(1, int((zhi - zlo) / 2.0));
  std::vector<double> zfill(size_t(nzb), 0.0);
  const double zper = double(M) / nzb;
  double zsum = 0;
  auto zb = [&](double z) { return size_t(std::clamp(int((z - zlo) / (zhi - zlo) * nzb), 0, nzb - 1)); };   // beads taken back by a backtrack
  std::vector<int> backtracks(seg.size(), 0);
  int backtracked = 0;
  auto margin_of = [&](const Vec3& p, const std::vector<uint32_t>& excl, double limit) {
    double m = 1e300;
    hash.near(p, [&](uint32_t id, double d) {
      if (d > limit + 1.0) return;
      if (dead[id]) return;
      for (uint32_t e : excl) if (e == id) return;
      m = std::min(m, d - limit);
    });
    return m;
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
      const int k = G.n - W.placed - 1;   // beads after this one
      const size_t L = W.path.size();
      const Vec3 A = W.path[L - 3], B = W.path[L - 2], C = W.path[L - 1];
      std::vector<uint32_t> excl(W.ids.end() - 3, W.ids.end());
      Vec3 P;
      if (closing && k == 0) {
        // the last bead closes on the target: on the circle of points 1.53 Å from both the head and the target
        excl.push_back(W.te1), excl.push_back(W.te2);
        const Vec3 E = W.e1, d = E - C;
        const double dd = norm(d);
        if (dd < 2 * kCC - 0.02) {
          const Vec3 m = (C + E) * 0.5, ax = unit(d);
          const double rho = std::sqrt(kCC * kCC - 0.25 * dd * dd);
          Vec3 u = unit(cross(ax, std::fabs(ax[0]) < 0.9 ? Vec3{1, 0, 0} : Vec3{0, 1, 0})), v = cross(ax, u);
          double best = 1e300;
          for (int q = 0; q < 72; ++q) {
            const double t = q * 5 * kDeg;
            const Vec3 x = m + (u * std::cos(t) + v * std::sin(t)) * rho;
            const double t1 = angle_deg(B, C, x) - 112, t2 = angle_deg(x, E, W.e2) - 112;
            const double mg = margin_of(x, excl, o.contact * 0.8);
            const double sc = (t1 * t1 + t2 * t2) / 100 + (mg < 0 ? 100 * mg * mg : 0) + ((x[2] < zlo || x[2] > zhi) ? 50 : 0);
            if (sc < best) best = sc, P = x;
          }
          if (best > 20) ++strained;
        } else {
          P = (C + E) * 0.5;
          ++strained;
        }
      } else {
        const bool dense = closing && k + 1 <= 3;
        const int ntr = dense ? 120 : std::max(4, o.trials);
        std::vector<Vec3> cand(static_cast<size_t>(ntr));
        std::vector<double> dist(size_t(ntr), 0.0);
        for (int t = 0; t < ntr; ++t) {
          double phi;
          if (dense) phi = (t * 3.0 + 1.5 * u01(rng)) * kDeg;
          else {
            const double r = u01(rng);
            const double c0 = r < 0.55 ? 180.0 : r < 0.775 ? 65.0 : -65.0;
            phi = (c0 + 15.0 * gauss(rng)) * kDeg;
          }
          cand[size_t(t)] = place(A, B, C, kCC, kTheta, phi);
        }
        // hard limits: inside the amorphous layer, and (closing) the target still reachable with the bonds left
        std::vector<char> ok(size_t(ntr), 1);
        std::vector<double> wcl(size_t(ntr), 1.0);
        for (int t = 0; t < ntr; ++t) {
          const Vec3& p = cand[size_t(t)];
          if (p[2] < zlo || p[2] > zhi) ok[size_t(t)] = 0;
          if (closing) {
            const int kb = k + 1;   // bonds from this bead to the target
            const double d = norm(W.e1 - p);
            dist[size_t(t)] = d;
            if (d > kb * kReach * (kb > 6 ? 0.8 : kb > 3 ? 0.9 : kb > 2 ? 0.97 : 1.0)) ok[size_t(t)] = 0;
            if (kb == 2 && (d < 2.2 || d > 2.85)) ok[size_t(t)] = 0;
            if (kb == 3 && d < 2.3) ok[size_t(t)] = 0;
            // ideal-chain closure probability, and a steep penalty once the target needs more than half the remaining reach
            const double f = d / (kb * kReach);
            wcl[size_t(t)] = std::exp(-3 * d * d / (2.0 * kb * kC * kCC * kCC) - (f > 0.5 ? 25 * (f - 0.5) : 0.0));
          }
        }
        std::vector<uint32_t> ex2 = excl;
        if (closing && k + 1 <= 3) ex2.push_back(W.te1), ex2.push_back(W.te2);
        if (closing && k + 1 <= 2) ex2.push_back(W.te3);
        int pick = -1;
        static const double scales[] = {1.0, 0.9, 0.8, 0.7, 0.6, 0.5};
        for (int si = 0; si < 6 && pick < 0; ++si) {
          const double lim = o.contact * scales[si];
          std::vector<double> w(size_t(ntr), 0.0);
          double sw = 0;
          for (int t = 0; t < ntr; ++t) {
            if (!ok[size_t(t)]) continue;
            if (margin_of(cand[size_t(t)], ex2, lim) < 0) continue;
            const double fill = zfill[zb(cand[size_t(t)][2])] / std::max(zper * 0.25, zsum / nzb);
            if (si < 3 && fill > 1.15) continue;   // a full slice: room is looked for further out first
            // loops (their first half) and tails lean away from their face, so the interphase is not packed with short folds
            double lean = 0;
            if (G.type == 0 && 2 * W.placed < G.n) lean = 0.8;
            else if (G.type == 2) lean = 0.6;
            const double up = (G.a % 2 ? 1.0 : -1.0) * (cand[size_t(t)][2] - C[2]);
            sw += w[size_t(t)] = wcl[size_t(t)] * std::exp(-8.0 * (fill - 1.0) + lean * up);
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
        if (pick < 0 && W.placed > 0 && backtracks[g] < 80 && std::none_of(ok.begin(), ok.end(), [](char c) { return c != 0; })) {
          // a dead end (the target out of reach or the layer's faces in the way): take back a few beads and grow them again
          const int back = std::min(W.placed, 2 + backtracks[g] / 8);
          for (int q = 0; q < back; ++q) {
            dead[G.beads.back()] = 1;
            zfill[zb(X[G.beads.back()][2])] -= 1;
            zsum -= 1;
            G.beads.pop_back();
            W.path.pop_back();
            W.ids.pop_back();
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
            const double sc = margin_of(cand[size_t(t)], ex2, o.contact) - (ok[size_t(t)] ? 0 : 100 + std::fabs(dist[size_t(t)]));
            if (sc > best) best = sc, pick = t;
          }
          ++eased;
          if (!ok[size_t(pick)]) ++strained;
        }
        P = cand[size_t(pick)];
      }
      const double mg = margin_of(P, excl, 0.0);
      worst = std::min(worst, mg);
      const uint32_t id = uint32_t(X.size());
      X.push_back(P);
      dead.push_back(0);
      hash.add(P);
      zfill[zb(P[2])] += 1;
      zsum += 1;
      G.beads.push_back(id);
      W.path.push_back(P);
      W.ids.push_back(id);
      ++W.placed;
      ++done_beads;
      if (W.placed < G.n) still.push_back(g);
    }
    active.swap(still);
    if (cancelled("growing the amorphous layer", double(done_beads) / std::max(1, M))) throw std::runtime_error("cancelled");
  }
  R.worst_contact = worst > 1e299 ? 0 : worst;
  if (backtracked) R.notes.push_back(std::to_string(backtracked) + " backtracks (beads taken back and grown again at dead ends)");
  R.eased_steps = eased, R.strained_closures = strained;

  // --- chains: from each tail through stems and segments to the other tail
  std::vector<std::vector<uint32_t>> chain;   // carbon ids in order
  std::vector<char> seen(seg.size(), 0);
  std::vector<std::vector<uint32_t>> segbeads_of_chain;
  for (size_t g0 = 0; g0 < seg.size(); ++g0) {
    if (seg[g0].type != 2 || seen[g0]) continue;
    std::vector<uint32_t> ch(seg[g0].beads.rbegin(), seg[g0].beads.rend());
    seen[g0] = 1;
    int e = seg[g0].a;
    for (;;) {
      const auto& st = stem[e / 2];
      if (e % 2) ch.insert(ch.end(), st.rbegin(), st.rend());
      else ch.insert(ch.end(), st.begin(), st.end());
      const int other = e ^ 1;
      const int gi = seg_of[other];
      Seg& G = seg[gi];
      seen[gi] = 1;
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

  // --- the System: carbons in chain order, each followed by its hydrogens
  System s;
  s.title = "semicrystalline " + pc.name + " lamellar stack";
  s.cell.a = {Lx, 0, 0}, s.cell.b = {0, Ly, 0}, s.cell.c = {0, 0, Lz};
  s.has_mol = true;
  s.tags = {{"crystal", "#3a7bd5"}, {"amorphous", "#e0913a"}};
  for (size_t ci = 0; ci < chain.size(); ++ci) {
    const auto& ch = chain[ci];
    const size_t m = ch.size();
    std::vector<uint32_t> cidx(m);
    auto mi = [&](const Vec3& p, const Vec3& ref) { return ref + box.mi(p - ref); };
    for (size_t k = 0; k < m; ++k) {
      const Vec3 cpos = X[ch[k]];
      Atom at;
      at.element = 6, at.name = "C", at.mol = int64_t(ci) + 1;
      at.pos = s.cell.wrap(cpos);
      at.tags = ch[k] < ncrys ? 1u : 2u;
      cidx[k] = uint32_t(s.atoms.size());
      s.atoms.push_back(at);
      if (k) s.bonds.push_back({cidx[k - 1], cidx[k], 1});
      std::vector<Vec3> hs;
      if (k > 0 && k + 1 < m) {
        const Vec3 pa = mi(X[ch[k - 1]], cpos), pb = mi(X[ch[k + 1]], cpos);
        const Vec3 u = unit(unit(cpos - pa) + unit(cpos - pb)), n = unit(cross(pa - cpos, pb - cpos));
        for (double sg : {1.0, -1.0}) hs.push_back(cpos + (u * std::cos(kHalfHCH) + n * (sg * std::sin(kHalfHCH))) * kCH);
      } else {   // a chain end: CH3, staggered
        const size_t k1 = k == 0 ? 1 : m - 2, k2 = k == 0 ? 2 : m - 3;
        const Vec3 p1 = mi(X[ch[k1]], cpos), p2 = mi(X[ch[std::min(k2, m - 1)]], cpos);
        for (double t : {60.0, 180.0, 300.0}) hs.push_back(place(p2, p1, cpos, kCH, 109.4712206 * kDeg, t * kDeg));
      }
      for (const Vec3& h : hs) {
        Atom ha;
        ha.element = 1, ha.name = "H", ha.mol = at.mol, ha.tags = at.tags;
        ha.pos = s.cell.wrap(h);
        s.bonds.push_back({cidx[k], uint32_t(s.atoms.size()), 1});
        s.atoms.push_back(ha);
      }
    }
  }
  int64_t id = 1;
  for (auto& at : s.atoms) at.id = id++;
  s.bonds_from_file = true;

  // --- report: the builder's topology, the structure analysed
  SemicrystallineReport A = analyse_semicrystalline(s);
  R.chains = int(chain.size());
  R.atoms = o.united_atom ? size_t(std::count_if(s.atoms.begin(), s.atoms.end(), [](const Atom& q) { return q.element == 6; })) : s.atoms.size();
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
    R.amorphous_density = mass / (Lx * Ly * La) * kAmuPerA3;
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
  // backbone carbons and their carbon neighbours
  std::vector<std::vector<uint32_t>> cn(n);
  for (const auto& bd : s.bonds)
    if (s.atoms[bd.i].element == 6 && s.atoms[bd.j].element == 6) cn[bd.i].push_back(bd.j), cn[bd.j].push_back(bd.i);
  std::vector<std::vector<uint32_t>> chains;
  std::vector<char> seen(n, 0);
  for (uint32_t i = 0; i < n; ++i) {
    if (s.atoms[i].element != 6 || seen[i] || cn[i].size() > 1) continue;
    std::vector<uint32_t> ch{i};
    seen[i] = 1;
    uint32_t prev = i, cur = i;
    for (;;) {
      uint32_t nx = UINT32_MAX;
      for (uint32_t j : cn[cur]) if (j != prev && !seen[j] && cn[j].size() <= 2) nx = j;
      if (nx == UINT32_MAX) break;
      seen[nx] = 1;
      ch.push_back(nx);
      prev = cur, cur = nx;
    }
    if (ch.size() >= 3) chains.push_back(ch);
  }
  R.chains = int(chains.size());
  // unwrapped positions along each chain
  std::vector<std::vector<Vec3>> up(chains.size());
  for (size_t c = 0; c < chains.size(); ++c) {
    up[c].push_back(s.atoms[chains[c][0]].pos);
    for (size_t k = 1; k < chains[c].size(); ++k) up[c].push_back(up[c][k - 1] + s.cell.minimum_image(s.atoms[chains[c][k]].pos - s.atoms[chains[c][k - 1]].pos));
  }
  auto local_p2 = [&](const Vec3& dir, std::vector<std::vector<double>>& out) {
    out.assign(chains.size(), {});
    for (size_t c = 0; c < chains.size(); ++c) {
      const size_t m = chains[c].size();
      std::vector<double> raw(m, std::nan(""));
      for (size_t k = 1; k + 1 < m; ++k) {
        const double ct = dot(unit(up[c][k + 1] - up[c][k - 1]), dir);
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
  std::vector<std::vector<double>> p2;
  Vec3 dir{0, 0, 1};
  local_p2(dir, p2);
  {   // the stems' mean direction from the crystalline chords, then P2 again against it
    Vec3 sum{0, 0, 0};
    for (size_t c = 0; c < chains.size(); ++c)
      for (size_t k = 1; k + 1 < chains[c].size(); ++k)
        if (p2[c][k] > p2_cut) {
          Vec3 u = unit(up[c][k + 1] - up[c][k - 1]);
          if (u[2] < 0) u = u * -1.0;
          sum = sum + u;
        }
    if (norm(sum) > 0) { dir = unit(sum); local_p2(dir, p2); }
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
  std::snprintf(b, sizeof b, "%zu atoms · %d stems · %d chains (%.0f CH2 on average)\n", r.atoms, r.stems, r.chains, mean(r.chain_units));
  out += b;
  std::snprintf(b, sizeof b, "segments: %d loops (%.0f %%; top %d, bottom %d) · %d bridges (%.0f %%) · %d tails (%.0f %%)\n", r.loops, 100 * r.loops / segs,
                r.loops_top, r.loops_bottom, r.bridges, 100 * r.bridges / segs, r.tails, 100 * r.tails / segs);
  out += b;
  std::snprintf(b, sizeof b, "mean length (CH2): loops %.1f · bridges %.1f · tails %.1f\n", mean(r.loop_units), mean(r.bridge_units), mean(r.tail_units));
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
    << "write_data semicrystalline_equilibrated.data\n";
  return o.str();
}

}  // namespace caps
