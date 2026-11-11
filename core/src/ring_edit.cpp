// Rings of 3 to 8 fused or spiro, an atom inserted into a bond, two structures joined (ring_edit.hpp).
#include "caps/ring_edit.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>

#include "caps/edit.hpp"
#include "caps/elements.hpp"

namespace caps {

namespace {

constexpr double kPi = 3.14159265358979323846;
const double kTetra = std::acos(-1.0 / 3.0);   // the tetrahedral angle, 109.47°
constexpr double kPucker = 0.25;                // Å: the start puckering of 6- to 8-rings (alternating out of the plane)

Vec3 unitv(const Vec3& v) { const double n = norm(v); return n > 1e-12 ? v * (1.0 / n) : Vec3{1, 0, 0}; }
Vec3 perpendicular(const Vec3& u) { const Vec3 t = std::fabs(u[0]) < 0.9 ? Vec3{1, 0, 0} : Vec3{0, 1, 0}; return unitv(cross(u, t)); }
Vec3 rotate(const Vec3& v, const Vec3& k, double t) {   // Rodrigues, k a unit vector
  return v * std::cos(t) + cross(k, v) * std::sin(t) + k * (dot(k, v) * (1 - std::cos(t)));
}

std::vector<std::vector<uint32_t>> neighbours(const System& s) {
  std::vector<std::vector<uint32_t>> nb(s.atoms.size());
  for (const auto& b : s.bonds) nb[b.i].push_back(b.j), nb[b.j].push_back(b.i);
  for (auto& x : nb) std::sort(x.begin(), x.end());
  return nb;
}
Vec3 rel(const System& s, uint32_t from, uint32_t to) {
  const Vec3 d = s.atoms[to].pos - s.atoms[from].pos;
  return s.cell.valid() ? s.cell.minimum_image(d) : d;
}
bool bonded(const System& s, uint32_t i, uint32_t j) {
  return std::any_of(s.bonds.begin(), s.bonds.end(), [&](const Bond& b) { return (b.i == i && b.j == j) || (b.i == j && b.j == i); });
}
double order_value(int order) { return order == kBondDative ? 0 : order == 2 ? 2 : order == 3 ? 3 : order == 4 ? 1.5 : 1; }
double order_sum(const System& s, uint32_t a) {
  double v = 0;
  for (const auto& b : s.bonds) if (b.i == a || b.j == a) v += order_value(b.order);
  return v;
}
int valence(const System& s, uint32_t a) {
  const double q = s.atoms[a].charge;
  const int c = int(std::lround(q));
  return default_valence(s.atoms[a].element, std::fabs(q - c) < 0.05 ? c : 0);
}
std::string label(const System& s, uint32_t a) { return std::string(element(s.atoms[a].element).symbol) + std::to_string(a + 1); }
double cov(int z) { return element(z).covalent; }

int type_for(System& s, int z) {
  const std::string sym = element(z).symbol;
  for (const auto& t : s.types) if (t.label == sym) return t.type;
  int next = 1;
  for (const auto& t : s.types) next = std::max(next, t.type + 1);
  TypeInfo ti;
  ti.type = next, ti.mass = element(z).mass, ti.label = sym;
  s.types.push_back(ti);
  return next;
}
uint32_t put(System& s, int z, const Vec3& p, uint32_t like) {
  Atom a;
  a.element = z;
  int64_t id = 0;
  for (const auto& x : s.atoms) id = std::max(id, x.id);
  a.id = id + 1;
  a.type = type_for(s, z);
  a.name = std::string(element(z).symbol) + std::to_string(s.atoms.size() + 1);
  a.pos = p;
  a.mol = s.atoms[like].mol, a.resname = s.atoms[like].resname, a.resid = s.atoms[like].resid, a.chain = s.atoms[like].chain;
  s.atoms.push_back(a);
  return uint32_t(s.atoms.size() - 1);
}

struct Member { int z = 6, exo = 0; };
std::vector<Member> parse_members(const std::string& text, int count) {
  std::string t = text;
  for (auto& c : t) if (c == ',' || c == ';') c = ' ';
  std::istringstream in(t);
  std::vector<Member> out;
  for (std::string w; in >> w;) {
    Member m;
    const auto eq = w.find('=');
    m.z = element_from_symbol(w.substr(0, eq));
    if (m.z <= 1) throw EditError("ring member '" + w + "': not a ring atom (an element symbol, optionally =O or =S)");
    if (eq != std::string::npos) {
      m.exo = element_from_symbol(w.substr(eq + 1));
      if (m.exo != 8 && m.exo != 16) throw EditError("ring member '" + w + "': the atom outside the ring is =O or =S");
    }
    out.push_back(m);
  }
  if (out.empty()) out.assign(size_t(count), Member{});
  if (int(out.size()) != count)
    throw EditError("the ring needs " + std::to_string(count) + " new atom(s), the members name " + std::to_string(out.size()));
  return out;
}

// the ring's new atoms placed, bonded in order from `first` round to `last`, their =O / =S added, hydrogens filled
std::vector<uint32_t> build_members(System& s, const std::vector<Member>& m, const std::vector<Vec3>& pos, uint32_t first, uint32_t last,
                                    int order, bool hydrogens, std::vector<uint32_t>* ring) {
  std::vector<uint32_t> made;
  for (size_t k = 0; k < m.size(); ++k) made.push_back(put(s, m[k].z, pos[k], first));
  uint32_t prev = first;
  for (uint32_t x : made) { s.bonds.push_back({prev, x, order}); prev = x; }
  s.bonds.push_back({prev, last, order});
  if (ring) ring->insert(ring->end(), made.begin(), made.end());
  std::vector<uint32_t> added = made;
  for (size_t k = 0; k < m.size(); ++k)
    if (m[k].exo) added.push_back(add_atom(s, int(made[k]), m[k].exo, 2));
  if (hydrogens) {
    std::vector<char> mask(s.atoms.size(), 0);
    for (uint32_t x : made) mask[x] = 1;
    const size_t before = s.atoms.size();
    add_hydrogens(s, mask);
    for (size_t x = before; x < s.atoms.size(); ++x) added.push_back(uint32_t(x));
  }
  return added;
}

// deletes the listed atoms and renumbers the result's lists
void drop(System& s, const std::vector<uint32_t>& gone, RingEditResult& r) {
  if (!s.velocities.empty() && s.velocities.size() < s.atoms.size()) s.velocities.resize(s.atoms.size(), Vec3{0, 0, 0});
  if (gone.empty()) return;
  std::vector<char> m(s.atoms.size(), 0);
  for (uint32_t g : gone) m[g] = 1;
  std::vector<uint32_t> shift(s.atoms.size(), 0);
  uint32_t c = 0;
  for (size_t i = 0; i < m.size(); ++i) { shift[i] = c; c += uint32_t(m[i]); }
  delete_atoms(s, m);
  for (auto* v : {&r.ring, &r.added}) for (auto& x : *v) x -= shift[x];
  r.removed = int(gone.size());
}

// UFF on the new atoms, the shared atoms and the shared atoms' hydrogens (the rest held)
void clean_region(System& s, const RingEditResult& r, const std::vector<uint32_t>& extra) {
  std::vector<char> mask(s.atoms.size(), 0);
  for (uint32_t x : r.added) mask[x] = 1;
  for (uint32_t x : extra) mask[x] = 1;
  const auto nb = neighbours(s);
  for (uint32_t x : extra) for (uint32_t q : nb[x]) if (s.atoms[q].element == 1) mask[q] = 1;
  clean_up(s, mask, 0.5);
}

// the terminal hydrogens of an atom (one bond each)
std::vector<uint32_t> terminal_h(const System& s, const std::vector<std::vector<uint32_t>>& nb, uint32_t a) {
  std::vector<uint32_t> h;
  for (uint32_t q : nb[a]) if (s.atoms[q].element == 1 && nb[q].size() == 1) h.push_back(q);
  return h;
}

void check_size(const RingSpec& sp) {
  if (sp.size < 3 || sp.size > 8) throw EditError("a ring of 3 to 8 atoms");
}

}  // namespace

RingEditResult fuse_ring(System& s, uint32_t i, uint32_t j, const RingSpec& sp) {
  check_size(sp);
  if (i >= s.atoms.size() || j >= s.atoms.size() || i == j) throw EditError("pick the two atoms of a bond");
  if (!bonded(s, i, j)) throw EditError("atoms " + label(s, i) + " and " + label(s, j) + " are not bonded");
  if (sp.aromatic && sp.size != 6) throw EditError("an aromatic fused ring has six atoms");
  const auto members = parse_members(sp.members, sp.size - 2);
  const auto nb = neighbours(s);
  const Vec3 bv = rel(s, i, j);
  const double L = norm(bv);
  const Vec3 e1 = unitv(bv);
  // the hydrogens that make way: one of each atom, on the same side
  const auto hi_all = terminal_h(s, nb, i), hj_all = terminal_h(s, nb, j);
  int hi = -1, hj = -1;
  if (!hi_all.empty() && !hj_all.empty()) {
    double best = -3;
    for (uint32_t a : hi_all)
      for (uint32_t b : hj_all) {
        const double c = dot(unitv(rel(s, i, a)), unitv(rel(s, j, b)));
        if (c > best) best = c, hi = int(a), hj = int(b);
      }
  } else if (!hi_all.empty()) {
    hi = int(hi_all[0]);
  } else if (!hj_all.empty()) {
    hj = int(hj_all[0]);
  }
  // an aromatic ring turns i-j into an aromatic bond and adds one more (1.5 each, half a bond of slack for the
  // aromatic counting); a plain ring adds one single bond
  int ij_order = 1;
  for (const auto& b : s.bonds) if ((b.i == i && b.j == j) || (b.i == j && b.j == i)) ij_order = b.order;
  for (auto [x, h] : {std::pair{i, hi}, std::pair{j, hj}}) {
    if (h >= 0) continue;
    const bool full = sp.aromatic ? order_sum(s, x) - order_value(ij_order) + 3.0 > valence(s, x) + 0.5 + 1e-9
                                  : order_sum(s, x) + 1 > valence(s, x) + 1e-9;
    if (full) throw EditError("atom " + label(s, x) + " has no hydrogen to give up and no free valence for the ring");
  }
  Vec3 u{0, 0, 0};
  if (hi >= 0) u = u + unitv(rel(s, i, uint32_t(hi)));
  if (hj >= 0) u = u + unitv(rel(s, j, uint32_t(hj)));
  if (hi < 0 && hj < 0) {   // the open side: away from the other neighbours
    for (uint32_t q : nb[i]) if (q != j) u = u - unitv(rel(s, i, q));
    for (uint32_t q : nb[j]) if (q != i) u = u - unitv(rel(s, j, q));
  }
  u = u - e1 * dot(u, e1);
  if (norm(u) < 1e-3) {
    if (hi >= 0 || hj >= 0) throw EditError("the hydrogens point along the bond: no side to fuse the ring on");
    u = perpendicular(e1);
  }
  u = unitv(u);
  const int n = sp.size;
  const double R = L / (2 * std::sin(kPi / n)), apo = L / (2 * std::tan(kPi / n));
  const Vec3 centre = s.atoms[i].pos + bv * 0.5 + u * apo, normal = unitv(cross(e1, u));
  std::vector<Vec3> pos;
  for (int k = 2; k < n; ++k) {   // vertex k of the polygon i (0), j (1), …
    const double t = (2 * k - 1) * kPi / n;
    Vec3 p = centre + (u * (-std::cos(t)) + e1 * std::sin(t)) * R;
    if (n >= 6 && !sp.aromatic) p = p + normal * (k % 2 ? kPucker : -kPucker);
    pos.push_back(p);
  }
  RingEditResult r;
  r.ring = {i, j};
  r.added = build_members(s, members, pos, j, i, sp.aromatic ? 4 : 1, sp.hydrogens, &r.ring);
  if (sp.aromatic)
    for (auto& b : s.bonds) if ((b.i == i && b.j == j) || (b.i == j && b.j == i)) b.order = 4;
  std::vector<uint32_t> gone;
  if (hi >= 0) gone.push_back(uint32_t(hi));
  if (hj >= 0) gone.push_back(uint32_t(hj));
  drop(s, gone, r);
  if (sp.clean) clean_region(s, r, {r.ring[0], r.ring[1]});
  r.notes.push_back(std::to_string(n) + "-membered ring fused onto " + label(s, r.ring[0]) + "–" + label(s, r.ring[1]));
  return r;
}

RingEditResult spiro_ring(System& s, uint32_t c, const RingSpec& sp) {
  check_size(sp);
  if (c >= s.atoms.size()) throw EditError("pick the atom the rings share");
  if (sp.aromatic) throw EditError("a spiro atom is saturated: the new ring cannot be aromatic");
  const auto members = parse_members(sp.members, sp.size - 1);
  const auto nb = neighbours(s);
  auto hs = terminal_h(s, nb, c);
  if (hs.size() > 2) hs.resize(2);
  std::vector<uint32_t> rest;
  for (uint32_t q : nb[c]) if (std::find(hs.begin(), hs.end(), q) == hs.end()) rest.push_back(q);
  double kept = 0;
  for (const auto& b : s.bonds)
    if ((b.i == c || b.j == c) && std::find(hs.begin(), hs.end(), b.i == c ? b.j : b.i) == hs.end()) kept += order_value(b.order);
  if (kept + 2 > valence(s, c) + 1e-9) throw EditError("atom " + label(s, c) + " has no room for two more ring bonds");
  Vec3 w1, w2;
  if (hs.size() == 2) {
    w1 = unitv(rel(s, c, hs[0])), w2 = unitv(rel(s, c, hs[1]));
  } else if (rest.size() == 2) {   // the two open tetrahedral directions
    const Vec3 d1 = unitv(rel(s, c, rest[0])), d2 = unitv(rel(s, c, rest[1]));
    const Vec3 b = unitv((d1 + d2) * -1.0), p = unitv(cross(d1, d2));
    const double half = kTetra / 2;
    w1 = b * std::cos(half) + p * std::sin(half), w2 = b * std::cos(half) - p * std::sin(half);
  } else {
    throw EditError("atom " + label(s, c) + " needs two hydrogens, or two neighbours and room for two bonds");
  }
  const int n = sp.size;
  double Lsum = 0;
  int prev = s.atoms[c].element;
  for (const auto& m : members) Lsum += cov(prev) + cov(m.z), prev = m.z;
  Lsum += cov(prev) + cov(s.atoms[c].element);
  const double L = Lsum / n, R = L / (2 * std::sin(kPi / n));
  Vec3 b = w1 + w2;
  b = norm(b) < 1e-6 ? perpendicular(w1) : unitv(b);
  const Vec3 q = unitv(w1 - b * dot(b, w1)), normal = unitv(cross(b, q)), centre = s.atoms[c].pos + b * R;
  std::vector<Vec3> pos;
  for (int k = 1; k < n; ++k) {
    const double t = 2 * kPi * k / n;
    Vec3 p = centre + (b * (-std::cos(t)) + q * std::sin(t)) * R;
    if (n >= 6) p = p + normal * (k % 2 ? kPucker : -kPucker);
    pos.push_back(p);
  }
  RingEditResult r;
  r.ring = {c};
  r.added = build_members(s, members, pos, c, c, 1, sp.hydrogens, &r.ring);
  drop(s, hs, r);
  if (sp.clean) clean_region(s, r, {r.ring[0]});
  r.notes.push_back(std::to_string(n) + "-membered spiro ring on " + label(s, r.ring[0]));
  return r;
}

RingEditResult insert_atom(System& s, uint32_t i, uint32_t j, int z, bool hydrogens, bool clean) {
  if (i >= s.atoms.size() || j >= s.atoms.size() || i == j) throw EditError("pick the two atoms of a bond");
  if (z <= 1) throw EditError("insert a heavy atom (not hydrogen)");
  if (!bonded(s, i, j)) throw EditError("atoms " + label(s, i) + " and " + label(s, j) + " are not bonded");
  int old = 0;
  for (const auto& b : s.bonds) if ((b.i == i && b.j == j) || (b.i == j && b.j == i)) old = b.order;
  RingEditResult r;
  if (old == 4) r.notes.push_back("the aromatic bond became two single bonds: the ring is no longer aromatic");
  remove_bond(s, i, j);
  const auto nb = neighbours(s);
  // the atoms reached from an atom with the bond gone; in a ring both sides are the same
  auto side = [&](uint32_t from) {
    std::vector<char> seen(s.atoms.size(), 0);
    std::vector<uint32_t> st{from}, out;
    seen[from] = 1;
    while (!st.empty()) {
      const uint32_t a = st.back();
      st.pop_back();
      out.push_back(a);
      for (uint32_t q : nb[a]) if (!seen[q]) seen[q] = 1, st.push_back(q);
    }
    return out;
  };
  const auto sj = side(j);
  const bool ring = std::find(sj.begin(), sj.end(), i) != sj.end();
  const double r1 = cov(s.atoms[i].element) + cov(z), r2 = cov(z) + cov(s.atoms[j].element);
  const Vec3 bv = rel(s, i, j), e1 = unitv(bv);
  double D = norm(bv);
  if (!ring) {   // the chain gives way: the smaller side moves out along the bond
    const double d = std::sqrt(r1 * r1 + r2 * r2 - 2 * r1 * r2 * std::cos(kTetra));
    const auto si = side(i);
    const bool move_j = sj.size() <= si.size();
    const Vec3 shift = e1 * ((d - D) * (move_j ? 1 : -1));
    for (uint32_t a : move_j ? sj : si) s.atoms[a].pos = s.atoms[a].pos + shift;
    if (!move_j) r.notes.push_back("the side of " + label(s, i) + " (" + std::to_string(si.size()) + " atoms) moved out along the bond");
    else r.notes.push_back("the side of " + label(s, j) + " (" + std::to_string(sj.size()) + " atoms) moved out along the bond");
    D = d;
  } else if (r1 + r2 < D) {
    throw EditError("the ring bond is longer than the two new bonds together");
  }
  // off the bond, away from the other neighbours of i and j
  const Vec3 pi = s.atoms[i].pos, mid = pi + e1 * (D / 2);
  Vec3 away{0, 0, 0};
  for (uint32_t a : {i, j})
    for (uint32_t q : nb[a]) away = away - (s.atoms[a].pos + rel(s, a, q) - mid);
  away = away - e1 * dot(away, e1);
  const Vec3 p = norm(away) > 1e-3 ? unitv(away) : perpendicular(e1);
  const double a = (r1 * r1 - r2 * r2 + D * D) / (2 * D), h = std::sqrt(std::max(0.0, r1 * r1 - a * a));
  const uint32_t x = put(s, z, pi + e1 * a + p * h, i);
  s.bonds.push_back({i, x, 1});
  s.bonds.push_back({x, j, 1});
  r.ring = {i, x, j};
  r.added = {x};
  if (hydrogens) {
    std::vector<char> mask(s.atoms.size(), 0);
    mask[x] = mask[i] = mask[j] = 1;
    const size_t before = s.atoms.size();
    add_hydrogens(s, mask);
    for (size_t k = before; k < s.atoms.size(); ++k) r.added.push_back(uint32_t(k));
  }
  drop(s, {}, r);
  if (clean) clean_region(s, r, {i, j});
  r.notes.insert(r.notes.begin(), std::string(element(z).symbol) + " inserted between " + label(s, i) + " and " + label(s, j));
  return r;
}

uint32_t default_leaving_atom(const System& b) {
  const auto nb = neighbours(b);
  for (uint32_t a = 0; a < b.atoms.size(); ++a)
    if (b.atoms[a].element != 1) {
      for (uint32_t q : nb[a]) if (b.atoms[q].element == 1 && nb[q].size() == 1) return q;
      break;
    }
  for (uint32_t a = 0; a < b.atoms.size(); ++a)
    if (b.atoms[a].element == 1 && nb[a].size() == 1) return a;
  for (uint32_t a = 0; a < b.atoms.size(); ++a)
    if (nb[a].size() == 1) return a;
  throw EditError("the piece has no terminal atom to leave");
}

RingEditResult join_pieces(System& a, const System& b, uint32_t la, uint32_t lb, const JoinOptions& o) {
  if (la >= a.atoms.size()) throw EditError("pick the leaving atom of the structure");
  if (lb >= b.atoms.size()) throw EditError("the piece has no atom " + std::to_string(lb + 1));
  if (o.order < 1 || o.order > 3) throw EditError("the new bond is single, double or triple");
  const auto na = neighbours(a), nbb = neighbours(b);
  if (na[la].size() != 1) throw EditError("the leaving atom " + label(a, la) + " must have exactly one bond");
  if (nbb[lb].size() != 1) throw EditError("the piece's leaving atom " + label(b, lb) + " must have exactly one bond");
  const uint32_t pa = na[la][0], pb = nbb[lb][0];
  const Vec3 da = unitv(rel(a, pa, la)), db = unitv(rel(b, pb, lb));
  const double L = cov(a.atoms[pa].element) + cov(b.atoms[pb].element) - (o.order == 2 ? 0.2 : o.order == 3 ? 0.34 : 0);   // add_atom's rule
  // turn the piece so its bond points along −da
  const Vec3 target = da * -1.0;
  Vec3 k = cross(db, target);
  const double ang = std::atan2(norm(k), dot(db, target));
  k = norm(k) > 1e-9 ? unitv(k) : perpendicular(db);
  const Vec3 at = a.atoms[pa].pos + da * L;
  std::vector<Vec3> local(b.atoms.size());   // the piece about pb, turned, before the torsion
  for (uint32_t q = 0; q < b.atoms.size(); ++q) local[q] = rotate(rel(b, pb, q), k, ang);
  auto placed = [&](uint32_t q, double phi) { return at + rotate(local[q], da, phi); };
  double phi = 0;
  int refa = -1, refb = -1;
  for (uint32_t q : na[pa]) if (q != la) { refa = int(q); break; }
  for (uint32_t q : nbb[pb]) if (q != lb) { refb = int(q); break; }
  auto dihedral = [&](double t) {   // refa–pa–pb–refb, degrees
    const Vec3 p0 = a.atoms[pa].pos + rel(a, pa, uint32_t(refa)), p1 = a.atoms[pa].pos, p2 = at, p3 = placed(uint32_t(refb), t);
    const Vec3 b1 = p1 - p0, b2 = p2 - p1, b3 = p3 - p2;
    const Vec3 n1 = cross(b1, b2), n2 = cross(b2, b3);
    return std::atan2(dot(cross(n1, n2), unitv(b2)), dot(n1, n2)) * 180 / kPi;
  };
  RingEditResult r;
  if (refa >= 0 && refb >= 0) {
    if (std::isfinite(o.torsion)) {
      phi = (o.torsion - dihedral(0)) * kPi / 180;   // turning the piece about the bond turns the dihedral one for one
      char t[160];
      std::snprintf(t, sizeof t, "dihedral %s–%s–(piece %s)–(piece %s) set to %.2f°", label(a, uint32_t(refa)).c_str(), label(a, pa).c_str(),
                    label(b, pb).c_str(), label(b, uint32_t(refb)).c_str(), o.torsion);
      r.notes.push_back(t);
    } else {   // the clearest turn: the largest smallest distance between the piece and the structure (bonding atoms aside)
      double reach = 0;
      for (uint32_t q = 0; q < b.atoms.size(); ++q) reach = std::max(reach, norm(local[q]));
      std::vector<Vec3> near;
      for (uint32_t q = 0; q < a.atoms.size(); ++q) {
        if (q == la || q == pa) continue;
        Vec3 d = a.atoms[q].pos - at;
        if (a.cell.valid()) d = a.cell.minimum_image(d);
        if (norm(d) < reach + 8) near.push_back(at + d);
      }
      auto clearance = [&](double t) {
        double m = std::numeric_limits<double>::infinity();
        for (uint32_t q = 0; q < b.atoms.size(); ++q) {
          if (q == lb || q == pb) continue;
          const Vec3 p = placed(q, t);
          for (const auto& x : near) m = std::min(m, norm(p - x));
        }
        return m;
      };
      if (!near.empty()) {
        double best = -1;
        for (int d = 0; d < 360; ++d) {
          const double c = clearance(d * kPi / 180);
          if (c > best + 1e-12) best = c, phi = d * kPi / 180;
        }
        double lo = phi - kPi / 180, hi = phi + kPi / 180;   // golden section on the bracket
        const double g = (std::sqrt(5.0) - 1) / 2;
        double x1 = hi - g * (hi - lo), x2 = lo + g * (hi - lo), f1 = clearance(x1), f2 = clearance(x2);
        for (int it = 0; it < 100 && hi - lo > 1e-11; ++it) {
          if (f1 < f2) lo = x1, x1 = x2, f1 = f2, x2 = lo + g * (hi - lo), f2 = clearance(x2);
          else hi = x2, x2 = x1, f2 = f1, x1 = hi - g * (hi - lo), f1 = clearance(x1);
        }
        phi = (lo + hi) / 2;
        char t[120];
        std::snprintf(t, sizeof t, "turned about the new bond for the most room: closest contact %.2f Å", clearance(phi));
        r.notes.push_back(t);
      }
    }
  }
  // the piece's atoms appended (its leaving atom skipped), types matched by label
  std::vector<int64_t> map(b.atoms.size(), -1);
  int64_t id = 0;
  for (const auto& x : a.atoms) id = std::max(id, x.id);
  for (uint32_t q = 0; q < b.atoms.size(); ++q) {
    if (q == lb) continue;
    Atom x = b.atoms[q];
    x.pos = placed(q, phi);
    x.id = ++id;
    x.mol = a.atoms[pa].mol, x.chain = a.atoms[pa].chain, x.tags = 0, x.image = {0, 0, 0};
    std::string lab = element(x.element).symbol;
    double mass = element(x.element).mass;
    for (const auto& t : b.types)
      if (t.type == x.type && !t.label.empty()) { lab = t.label; if (t.mass > 0) mass = t.mass; }
    int type = 0;
    for (const auto& t : a.types) if (t.label == lab) type = t.type;
    if (type == 0) {
      for (const auto& t : a.types) type = std::max(type, t.type);
      ++type;
      TypeInfo ti;
      ti.type = type, ti.mass = mass, ti.label = lab;
      a.types.push_back(ti);
    }
    x.type = type;
    map[q] = int64_t(a.atoms.size());
    r.added.push_back(uint32_t(a.atoms.size()));
    a.atoms.push_back(x);
  }
  for (const auto& bd : b.bonds)
    if (map[bd.i] >= 0 && map[bd.j] >= 0) a.bonds.push_back({uint32_t(map[bd.i]), uint32_t(map[bd.j]), bd.order});
  a.bonds.push_back({pa, uint32_t(map[pb]), o.order});
  if (b.has_charges) a.has_charges = true;
  if (b.cell.valid()) r.notes.push_back("the piece's cell is not kept: the joined structure keeps the first one's");
  r.ring = {pa, uint32_t(map[pb])};
  drop(a, {la}, r);
  if (o.clean) clean_region(a, r, {r.ring[0], r.ring[1]});
  r.notes.insert(r.notes.begin(), std::to_string(r.added.size()) + " atoms joined by the bond " + label(a, r.ring[0]) + "–" + label(a, r.ring[1]));
  return r;
}

}  // namespace caps
