// CAPS structure editing (see caps/edit.hpp).
#include "caps/edit.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <sstream>

#include "caps/analysis.hpp"
#include "caps/elements.hpp"
#include "caps/relax.hpp"
#include "caps/typing.hpp"
#include "caps/uff.hpp"

namespace caps {

namespace {

constexpr double kTet = 109.4712206 * 3.14159265358979323846 / 180.0;
Vec3 unitv(const Vec3& v) { const double n = norm(v); return n > 1e-12 ? v * (1.0 / n) : Vec3{1, 0, 0}; }
Vec3 perpendicular(const Vec3& u) { const Vec3 t = std::fabs(u[0]) < 0.9 ? Vec3{1, 0, 0} : Vec3{0, 1, 0}; return unitv(cross(u, t)); }

std::vector<std::vector<uint32_t>> neighbours(const System& s) {
  std::vector<std::vector<uint32_t>> nb(s.atoms.size());
  for (const auto& b : s.bonds) nb[b.i].push_back(b.j), nb[b.j].push_back(b.i);
  return nb;
}

double bond_order_sum(const System& s, uint32_t a) {
  double v = 0;
  for (const auto& b : s.bonds)
    if (b.i == a || b.j == a) v += b.order == 2 ? 2 : b.order == 3 ? 3 : b.order == 4 ? 1.5 : 1;
  return v;
}

bool has_multiple_bond(const System& s, uint32_t a) {
  for (const auto& b : s.bonds)
    if ((b.i == a || b.j == a) && (b.order == 2 || b.order == 3 || b.order == 4)) return true;
  return false;
}

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

Vec3 rel(const System& s, uint32_t from, uint32_t to) {
  Vec3 d = s.atoms[to].pos - s.atoms[from].pos;
  return s.cell.valid() ? s.cell.minimum_image(d) : d;
}

// The missing bond directions of an atom with `existing` unit bond vectors and `total` bonds in its hybridisation;
// `ref` (a unit vector) fixes the rotation when only one bond exists.
std::vector<Vec3> ideal_directions(const std::vector<Vec3>& existing, int total, const Vec3& ref) {
  std::vector<Vec3> out;
  const size_t k = existing.size();
  const int missing = total - int(k);
  if (missing <= 0) return out;
  if (k == 0) {
    const Vec3 a{1, 0, 0};
    out.push_back(a);
    if (total == 2) out.push_back(a * -1.0);
    if (total == 3) for (int q = 1; q < 3; ++q) out.push_back({std::cos(q * 2.0944), std::sin(q * 2.0944), 0});
    if (total >= 4) {
      out = {unitv({1, 1, 1}), unitv({1, -1, -1}), unitv({-1, 1, -1}), unitv({-1, -1, 1})};
    }
    out.resize(size_t(missing));
    return out;
  }
  if (total == 2) { out.push_back(unitv(existing[0] * -1.0)); return out; }   // sp: opposite
  if (k == 1) {
    const Vec3 u = existing[0];
    Vec3 p = ref - u * dot(ref, u);
    p = norm(p) > 1e-6 ? unitv(p) : perpendicular(u);
    const Vec3 q = cross(u, p);
    const double ang = total == 3 ? 2.0944 : kTet;   // 120° or 109.47° from the bond
    const int n = total - 1;
    for (int m = 0; m < missing; ++m) {
      // staggered around the bond: anti to the reference first
      const double t = 3.14159265358979323846 + 2 * 3.14159265358979323846 * m / n;
      out.push_back(unitv(u * std::cos(ang) + (p * std::cos(t) + q * std::sin(t)) * std::sin(ang)));
    }
    return out;
  }
  Vec3 sum{0, 0, 0};
  for (const auto& e : existing) sum = sum + e;
  const Vec3 b = unitv(sum * -1.0);
  if (k == 2 && total == 4 && missing == 2) {
    const Vec3 n = unitv(cross(existing[0], existing[1]));
    const double h = 0.5 * kTet;   // half the tetrahedral angle
    out.push_back(unitv(b * std::cos(h) + n * std::sin(h)));
    out.push_back(unitv(b * std::cos(h) - n * std::sin(h)));
    return out;
  }
  out.push_back(b);   // k = 2 in sp² (or one of sp³'s), k = 3 in sp³
  if (missing > 1) {
    const Vec3 n = norm(cross(existing[0], existing[1])) > 1e-6 ? unitv(cross(existing[0], existing[1])) : perpendicular(b);
    out.push_back(n);
  }
  out.resize(size_t(missing));
  return out;
}

double h_length(int z) { return z == 6 ? 1.09 : z == 7 ? 1.01 : z == 8 ? 0.96 : z == 16 ? 1.34 : element(z).covalent + 0.31; }

// Atoms reached from `start` without crossing `blocked` (the branch start hangs from blocked).
std::vector<uint32_t> branch(const std::vector<std::vector<uint32_t>>& nb, uint32_t blocked, uint32_t start, bool* ring) {
  std::vector<uint32_t> out{start}, stack{start};
  std::vector<char> seen(nb.size(), 0);
  seen[start] = seen[blocked] = 1;
  if (ring) *ring = false;
  while (!stack.empty()) {
    const uint32_t a = stack.back();
    stack.pop_back();
    for (uint32_t q : nb[a]) {
      if (q == blocked && a != start) { if (ring) *ring = true; continue; }
      if (seen[q]) continue;
      seen[q] = 1;
      out.push_back(q);
      stack.push_back(q);
    }
  }
  return out;
}

}  // namespace

int default_valence(int z, int charge) {
  switch (z) {
    case 1: return charge == 0 ? 1 : 0;
    case 5: return 3;
    case 6: return charge == 0 ? 4 : 3;
    case 7: return charge > 0 ? 4 : charge < 0 ? 2 : 3;
    case 8: return charge > 0 ? 3 : charge < 0 ? 1 : 2;
    case 9: case 17: case 35: case 53: return charge == 0 ? 1 : 0;
    case 14: return 4;
    case 15: return charge > 0 ? 4 : 3;
    case 16: return charge < 0 ? 1 : 2;
    default: return 0;
  }
}

void set_element(System& s, uint32_t atom, int z) {
  if (atom >= s.atoms.size()) throw EditError("no such atom");
  if (z <= 0) throw EditError("unknown element");
  auto& a = s.atoms[atom];
  a.element = z;
  a.type = type_for(s, z);
  a.name = std::string(element(z).symbol) + std::to_string(atom + 1);
}

uint32_t add_atom(System& s, int bonded_to, int z, int order, int geometry, Vec3 at) {
  if (z <= 0) throw EditError("unknown element");
  Atom a;
  a.element = z;
  a.id = s.atoms.empty() ? 1 : s.atoms.back().id + 1;
  a.type = type_for(s, z);
  a.name = std::string(element(z).symbol) + std::to_string(s.atoms.size() + 1);
  if (bonded_to < 0) {
    a.pos = at;
    a.mol = 0;
    for (const auto& x : s.atoms) a.mol = std::max(a.mol, x.mol);
    a.mol += 1;
    s.atoms.push_back(a);
    return uint32_t(s.atoms.size() - 1);
  }
  const uint32_t p = uint32_t(bonded_to);
  if (p >= s.atoms.size()) throw EditError("no such atom");
  const auto nb = neighbours(s);
  std::vector<Vec3> dirs;
  for (uint32_t q : nb[p]) dirs.push_back(unitv(rel(s, p, q)));
  int total = geometry == 1 ? 2 : geometry == 2 ? 3 : geometry == 3 ? 4 : 0;
  if (total == 0) total = order >= 3 || (has_multiple_bond(s, p) && std::any_of(s.bonds.begin(), s.bonds.end(), [&](const Bond& b) { return (b.i == p || b.j == p) && b.order == 3; })) ? 2
                        : order == 2 || has_multiple_bond(s, p) ? 3 : 4;
  total = std::max(total, int(dirs.size()) + 1);
  Vec3 ref = perpendicular(dirs.empty() ? Vec3{1, 0, 0} : dirs[0]);
  if (!nb[p].empty() && !nb[nb[p][0]].empty())   // a torsion reference: a neighbour of the first neighbour
    for (uint32_t r : nb[nb[p][0]]) if (r != p) { ref = unitv(rel(s, nb[p][0], r)); break; }
  const auto cand = ideal_directions(dirs, total, ref);
  const Vec3 d = cand.empty() ? perpendicular(dirs.empty() ? Vec3{1, 0, 0} : dirs[0]) : cand[0];
  const double len = z == 1 || s.atoms[p].element == 1 ? h_length(z == 1 ? s.atoms[p].element : z)
                                                        : element(z).covalent + element(s.atoms[p].element).covalent - (order == 2 ? 0.2 : order == 3 ? 0.34 : 0);
  a.pos = s.atoms[p].pos + d * len;
  a.mol = s.atoms[p].mol;
  a.resname = s.atoms[p].resname;
  a.resid = s.atoms[p].resid;
  s.atoms.push_back(a);
  const uint32_t id = uint32_t(s.atoms.size() - 1);
  s.bonds.push_back({p, id, order});
  return id;
}

void add_bond(System& s, uint32_t i, uint32_t j, int order) {
  if (i >= s.atoms.size() || j >= s.atoms.size() || i == j) throw EditError("a bond needs two different atoms");
  for (auto& b : s.bonds)
    if ((b.i == i && b.j == j) || (b.i == j && b.j == i)) { b.order = order; return; }   // an existing bond: its order
  s.bonds.push_back({i, j, order});
  // one molecule now: the second atom's molecule joins the first's
  const int64_t keep = s.atoms[i].mol, drop = s.atoms[j].mol;
  if (keep != drop && keep != 0)
    for (auto& a : s.atoms) if (a.mol == drop) a.mol = keep;
}

bool remove_bond(System& s, uint32_t i, uint32_t j) {
  const auto n0 = s.bonds.size();
  s.bonds.erase(std::remove_if(s.bonds.begin(), s.bonds.end(), [&](const Bond& b) { return (b.i == i && b.j == j) || (b.i == j && b.j == i); }), s.bonds.end());
  return s.bonds.size() != n0;
}

void delete_atoms(System& s, const std::vector<char>& remove) {
  const size_t n = s.atoms.size();
  std::vector<int64_t> map(n, -1);
  std::vector<Atom> kept;
  for (size_t i = 0; i < n; ++i)
    if (i >= remove.size() || !remove[i]) { map[i] = int64_t(kept.size()); kept.push_back(s.atoms[i]); }
  std::vector<Bond> bonds;
  for (const auto& b : s.bonds)
    if (map[b.i] >= 0 && map[b.j] >= 0) bonds.push_back({uint32_t(map[b.i]), uint32_t(map[b.j]), b.order});
  if (!s.velocities.empty() && s.velocities.size() == n) {
    std::vector<Vec3> v;
    for (size_t i = 0; i < n; ++i) if (map[i] >= 0) v.push_back(s.velocities[i]);
    s.velocities = v;
  }
  s.atoms = std::move(kept);
  s.bonds = std::move(bonds);
}

int add_hydrogens(System& s, const std::vector<char>& atoms) {
  const size_t n = s.atoms.size();
  int added = 0;
  for (size_t i = 0; i < n; ++i) {
    if (!atoms.empty() && (i >= atoms.size() || !atoms[i])) continue;
    const int z = s.atoms[i].element;
    if (z == 1) continue;
    const int charge = int(std::lround(s.atoms[i].charge));
    const int val = default_valence(z, std::fabs(s.atoms[i].charge - charge) < 0.05 ? charge : 0);
    const int missing = val - int(std::lround(bond_order_sum(s, uint32_t(i)) - 1e-9));
    if (val == 0 || missing <= 0) continue;
    const auto nb = neighbours(s);
    std::vector<Vec3> dirs;
    for (uint32_t q : nb[i]) dirs.push_back(unitv(rel(s, uint32_t(i), q)));
    const int bonds_now = int(nb[i].size());
    // hybridisation from the bonds it will have: σ bonds + lone pairs as for the neutral atom
    int total = bonds_now + missing;
    if (has_multiple_bond(s, uint32_t(i))) total = std::any_of(s.bonds.begin(), s.bonds.end(), [&](const Bond& b) { return (b.i == i || b.j == i) && b.order == 3; }) ? 2 : 3;
    else if ((z == 7 && charge == 0) || (z == 8 && charge == 0) || z == 16) total = 4;   // lone pairs fill the tetrahedron
    else total = std::max(total, 4);
    total = std::max(total, bonds_now + missing);
    Vec3 ref = dirs.empty() ? Vec3{0, 0, 1} : perpendicular(dirs[0]);
    if (!nb[i].empty())
      for (uint32_t r : nb[nb[i][0]]) if (r != i) { ref = unitv(rel(s, nb[i][0], r)); break; }
    const auto cand = ideal_directions(dirs, total, ref);
    for (int k = 0; k < missing && k < int(cand.size()); ++k) {
      Atom h;
      h.element = 1;
      h.type = type_for(s, 1);
      h.id = s.atoms.back().id + 1;
      h.mol = s.atoms[i].mol, h.resname = s.atoms[i].resname, h.resid = s.atoms[i].resid;
      h.name = "H" + std::to_string(s.atoms.size() + 1);
      h.pos = s.atoms[i].pos + cand[size_t(k)] * h_length(z);
      s.atoms.push_back(h);
      s.bonds.push_back({uint32_t(i), uint32_t(s.atoms.size() - 1), 1});
      ++added;
    }
  }
  return added;
}

namespace {

// Swaps the branches hanging from centre c at neighbours a and b: a 180° rotation of both about the bisector of their
// bonds, which inverts the centre's configuration.
void swap_branches(System& s, uint32_t c, uint32_t a, uint32_t b) {
  const auto nb = neighbours(s);
  bool ra = false, rb = false;
  const auto A = branch(nb, c, a, &ra), B = branch(nb, c, b, &rb);
  if (ra || rb) throw EditError("a swapped branch is in a ring with the centre");
  const Vec3 C = s.atoms[c].pos;
  const Vec3 axis = unitv(unitv(rel(s, c, a)) + unitv(rel(s, c, b)));
  for (const auto* br : {&A, &B})
    for (uint32_t i : *br) {
      Vec3 p = s.atoms[i].pos - C;
      if (s.cell.valid()) p = s.cell.minimum_image(p);
      s.atoms[i].pos = C + axis * (2 * dot(axis, p)) - p;
    }
}

}  // namespace

void invert_centre(System& s, uint32_t c) {
  if (c >= s.atoms.size()) throw EditError("no such atom");
  const auto nb = neighbours(s);
  if (nb[c].size() != 4) throw EditError("an inversion centre has four neighbours");
  // the two smallest branches that are not in a ring with the centre
  std::vector<std::pair<size_t, uint32_t>> br;
  for (uint32_t q : nb[c]) {
    bool ring = false;
    const auto b = branch(nb, c, q, &ring);
    if (!ring) br.emplace_back(b.size(), q);
  }
  if (br.size() < 2) throw EditError("the centre is in rings on three sides: swap by hand");
  std::sort(br.begin(), br.end());
  swap_branches(s, c, br[0].second, br[1].second);
}

TacticityReport tacticity(const System& s) {
  TacticityReport R;
  const auto nb = neighbours(s);
  for (const auto& chain : backbones(s, 4)) {
    std::set<uint32_t> on(chain.begin(), chain.end());
    TacticityChain T;
    for (size_t k = 1; k + 1 < chain.size(); ++k) {
      const uint32_t b = chain[k];
      if (nb[b].size() != 4) continue;
      std::vector<uint32_t> side;
      for (uint32_t q : nb[b]) if (!on.count(q)) side.push_back(q);
      if (side.size() != 2) continue;
      const bool h0 = s.atoms[side[0]].element == 1, h1 = s.atoms[side[1]].element == 1;
      if (h0 == h1) continue;   // CH2 or C(R)(R'): not a vinyl stereocentre here
      const uint32_t group = h0 ? side[1] : side[0];
      const double v = dot(cross(rel(s, b, chain[k - 1]), rel(s, b, chain[k + 1])), rel(s, b, group));
      T.centres.push_back(b);
      T.sign.push_back(v > 0 ? 1 : -1);
    }
    for (size_t k = 1; k < T.sign.size(); ++k) T.dyads += T.sign[k] == T.sign[k - 1] ? 'm' : 'r';
    for (char d : T.dyads) (d == 'm' ? R.m : R.r) += 1;
    for (size_t k = 1; k < T.dyads.size(); ++k) {
      const std::string t{T.dyads[k - 1], T.dyads[k]};
      if (t == "mm") ++R.mm; else if (t == "rr") ++R.rr; else ++R.mr;
    }
    if (!T.centres.empty()) R.chains.push_back(std::move(T));
  }
  if (R.m + R.r > 0) R.label = R.r == 0 ? "isotactic" : R.m == 0 ? "syndiotactic" : "atactic";
  return R;
}

int set_tacticity(System& s, bool iso) {
  const auto T = tacticity(s);
  int inverted = 0;
  for (const auto& ch : T.chains) {
    const std::set<uint32_t> chain_atoms = [&] {
      std::set<uint32_t> on;
      for (const auto& bb : backbones(s, 4)) for (uint32_t x : bb) on.insert(x);
      return on;
    }();
    for (size_t k = 1; k < ch.centres.size(); ++k) {
      const int want = iso ? ch.sign[0] : (k % 2 ? -ch.sign[0] : ch.sign[0]);
      if (ch.sign[k] == want) continue;
      // swap the centre's H and its side group; the backbone stays where it is
      const uint32_t b = ch.centres[k];
      const auto nb = neighbours(s);
      int h = -1, g = -1;
      for (uint32_t q : nb[b]) {
        if (chain_atoms.count(q)) continue;
        (s.atoms[q].element == 1 ? h : g) = int(q);
      }
      if (h < 0 || g < 0) continue;
      swap_branches(s, b, uint32_t(h), uint32_t(g));
      ++inverted;
    }
  }
  return inverted;
}

void clean_up(System& s, const std::vector<char>& atoms, double ftol) {
  auto ff = std::make_shared<ForceField>(assign_uff(s));
  RelaxOptions ro;
  ro.field = ff;
  ro.ftol = ftol;
  ro.max_iterations = 4000;
  ro.pushoff = true;
  if (!atoms.empty()) {
    ro.fixed.assign(s.atoms.size(), 1);
    for (size_t i = 0; i < s.atoms.size() && i < atoms.size(); ++i) ro.fixed[i] = atoms[i] ? 0 : 1;
  }
  relax(s, ro);
}

std::vector<char> select_smarts(const System& s, const std::string& pattern) {
  const Smarts sm(pattern);
  const Perception p = perceive(s);
  std::vector<char> out(s.atoms.size(), 0);
  for (size_t i = 0; i < s.atoms.size(); ++i) out[i] = sm.matches(s, p, uint32_t(i), {}) ? 1 : 0;
  return out;
}

std::vector<char> select_element(const System& s, const std::string& symbols) {
  std::set<int> zs;
  std::istringstream in(symbols);
  for (std::string t; in >> t;) {
    for (auto& ch : t) if (ch == ',') ch = ' ';
    std::istringstream in2(t);
    for (std::string u; in2 >> u;) {
      const int z = element_from_symbol(u);
      if (z <= 0) throw EditError("unknown element '" + u + "'");
      zs.insert(z);
    }
  }
  std::vector<char> out(s.atoms.size(), 0);
  for (size_t i = 0; i < s.atoms.size(); ++i) out[i] = zs.count(s.atoms[i].element) ? 1 : 0;
  return out;
}

std::vector<char> select_type(const System& s, const std::string& label) {
  std::set<int> types;
  for (const auto& t : s.types) if (t.label == label || std::to_string(t.type) == label) types.insert(t.type);
  if (types.empty()) throw EditError("no type '" + label + "' in this structure");
  std::vector<char> out(s.atoms.size(), 0);
  for (size_t i = 0; i < s.atoms.size(); ++i) out[i] = types.count(s.atoms[i].type) ? 1 : 0;
  return out;
}

std::vector<char> select_charge(const System& s, double lo, double hi) {
  std::vector<char> out(s.atoms.size(), 0);
  for (size_t i = 0; i < s.atoms.size(); ++i) out[i] = s.atoms[i].charge >= lo && s.atoms[i].charge <= hi ? 1 : 0;
  return out;
}

std::vector<char> select_within(const System& s, const std::vector<char>& from, double distance) {
  std::vector<char> out(s.atoms.size(), 0);
  std::vector<uint32_t> seeds;
  for (size_t i = 0; i < s.atoms.size() && i < from.size(); ++i) if (from[i]) seeds.push_back(uint32_t(i));
  if (seeds.empty()) return out;
  const double d2 = distance * distance;
  for (size_t i = 0; i < s.atoms.size(); ++i) {
    if (i < from.size() && from[i]) { out[i] = 1; continue; }
    for (uint32_t j : seeds) {
      const Vec3 d = rel(s, j, uint32_t(i));
      if (dot(d, d) <= d2) { out[i] = 1; break; }
    }
  }
  return out;
}

std::vector<char> select_grow(const System& s, const std::vector<char>& from, int steps) {
  std::vector<char> out(s.atoms.size(), 0);
  for (size_t i = 0; i < s.atoms.size() && i < from.size(); ++i) out[i] = from[i];
  for (int k = 0; k < steps; ++k) {
    std::vector<char> next = out;
    for (const auto& b : s.bonds) {
      if (out[b.i]) next[b.j] = 1;
      if (out[b.j]) next[b.i] = 1;
    }
    out = std::move(next);
  }
  return out;
}

}  // namespace caps
