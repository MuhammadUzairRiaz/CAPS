// CAPS structure editing (see caps/edit.hpp).
#include "caps/edit.hpp"

#include <algorithm>
#include <random>
#include <cmath>
#include <map>
#include <set>
#include <sstream>

#include "caps/analysis.hpp"
#include "caps/elements.hpp"
#include "caps/molecule.hpp"
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

namespace {
// Hydrogens an atom lacks: its usual valence (formal charge counted) less the bond orders it has (aromatic 1.5).
int missing_h(const System& s, uint32_t i, const std::vector<double>& order_sum) {
  const int z = s.atoms[i].element;
  if (z == 1) return 0;
  const int charge = int(std::lround(s.atoms[i].charge));
  const int val = default_valence(z, std::fabs(s.atoms[i].charge - charge) < 0.05 ? charge : 0);
  if (val == 0) return 0;
  return std::max(0, val - int(std::lround(order_sum[i] - 1e-9)));
}
std::vector<double> order_sums(const System& s) {
  std::vector<double> v(s.atoms.size(), 0.0);
  for (const auto& b : s.bonds) {
    const double o = b.order == 2 ? 2 : b.order == 3 ? 3 : b.order == 4 ? 1.5 : 1;
    v[b.i] += o, v[b.j] += o;
  }
  return v;
}
}  // namespace

std::vector<HydrogenPlanRow> hydrogen_plan(const System& s, const std::vector<char>& atoms) {
  const auto sums = order_sums(s);
  const auto nb = neighbours(s);
  std::vector<int> kind(s.atoms.size(), 3);   // 0 aromatic, 1 sp, 2 sp², 3 sp³
  for (const auto& b : s.bonds) {
    const int k = b.order == 4 ? 0 : b.order == 3 ? 1 : b.order == 2 ? 2 : 3;
    kind[b.i] = std::min(kind[b.i], k), kind[b.j] = std::min(kind[b.j], k);
  }
  std::map<std::tuple<int, int, int, bool>, HydrogenPlanRow> rows;   // (kind, element, heavy neighbours, all C)
  for (uint32_t i = 0; i < s.atoms.size(); ++i) {
    if (s.atoms[i].element == 1 || (!atoms.empty() && (i >= atoms.size() || !atoms[i]))) continue;
    int heavy = 0;
    bool all_c = true;
    for (uint32_t q : nb[i])
      if (s.atoms[q].element != 1) ++heavy, all_c = all_c && s.atoms[q].element == 6;
    auto& r = rows[{kind[i], s.atoms[i].element, heavy, all_c}];
    if (r.label.empty()) {
      static const char* hyb[] = {"aromatic", "sp", "sp²", "sp³"};
      const std::string sym = element(s.atoms[i].element).symbol;
      r.label = std::string(hyb[kind[i]]) + " " + sym + " with " + std::to_string(heavy) + (all_c && heavy ? " C" : " heavy") + (heavy == 1 ? " neighbour" : " neighbours");
    }
    ++r.atoms;
    r.hydrogens += missing_h(s, i, sums);
  }
  std::vector<HydrogenPlanRow> out;
  for (auto& [k, r] : rows) out.push_back(std::move(r));
  return out;
}

int add_hydrogens(System& s, const std::vector<char>& atoms) {
  const size_t n = s.atoms.size();
  int added = 0;
  const auto sums = order_sums(s);
  const auto nb = neighbours(s);   // each atom's own bonds: the H added to earlier atoms do not change them
  for (size_t i = 0; i < n; ++i) {
    if (!atoms.empty() && (i >= atoms.size() || !atoms[i])) continue;
    const int z = s.atoms[i].element;
    if (z == 1) continue;
    const int charge = int(std::lround(s.atoms[i].charge));
    const int missing = missing_h(s, uint32_t(i), sums);
    if (missing <= 0) continue;
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

namespace {

// "*" → a hydrogen stand-in (a written atom, so the attachment points keep their place in the atom order).
std::string star_as_hydrogen(const std::string& smiles) {
  std::string out;
  for (size_t k = 0; k < smiles.size(); ++k) {
    if (smiles[k] != '*') { out += smiles[k]; continue; }
    out += k > 0 && smiles[k - 1] == '[' ? "H" : "[H]";
  }
  return out;
}

struct Fragment3D {
  System s;
  std::vector<uint32_t> dummy, root;   // each attachment point and the atom it hangs on
};

const Fragment3D& fragment_3d(const std::string& smiles) {
  static std::map<std::string, Fragment3D> cache;
  auto it = cache.find(smiles);
  if (it != cache.end()) return it->second;
  const MolGraph g = parse_smiles(smiles);
  std::vector<uint32_t> stars;
  for (int i = 0; i < g.heavy; ++i) if (g.atoms[size_t(i)].element == 0) stars.push_back(uint32_t(i));
  if (stars.empty()) throw EditError("a fragment needs at least one * attachment point");
  BuildOptions bo;
  bo.forcefield = "uff";
  bo.seed = 5;
  Fragment3D f;
  f.s = build_molecule(star_as_hydrogen(smiles), bo).system;
  for (uint32_t d : stars) {
    int root = -1;
    for (const auto& b : f.s.bonds) {
      if (b.i == d) root = int(b.j);
      if (b.j == d) root = int(b.i);
    }
    if (root < 0) throw EditError("an attachment point must be bonded to one atom");
    f.dummy.push_back(d);
    f.root.push_back(uint32_t(root));
  }
  return cache.emplace(smiles, std::move(f)).first->second;
}

}  // namespace

std::vector<int> fragment_attach_atoms(const std::string& smiles) {
  const auto& f = fragment_3d(smiles);
  std::vector<int> out;
  for (uint32_t r : f.root) out.push_back(int(r));
  return out;
}

std::string silane_smiles(const std::string& name) {
  // each bonded through one silicon (one ethoxy condensed with a silanol); a bis-silane keeps its second silicon free
  if (name == "TESPT" || name == "Si69") return "*[Si](OCC)(OCC)CCCSSSSCCC[Si](OCC)(OCC)OCC";
  if (name == "TESPD" || name == "Si75") return "*[Si](OCC)(OCC)CCCSSCCC[Si](OCC)(OCC)OCC";
  if (name == "MPTES") return "*[Si](OCC)(OCC)CCCS";
  if (name == "APTES") return "*[Si](OCC)(OCC)CCCN";
  if (name == "VTES") return "*[Si](OCC)(OCC)C=C";
  if (name == "OCTEO") return "*[Si](OCC)(OCC)CCCCCCCC";
  throw EditError("unknown silane '" + name + "' (TESPT, TESPD, MPTES, APTES, VTES, OCTEO, or a SMILES with *)");
}

GraftReport graft_silanes(System& s, const GraftOptions& o) {
  GraftReport rep;
  if (o.smiles.find('*') == std::string::npos) throw EditError("the silane SMILES needs * where it bonds to the surface oxygen");
  const auto nb = neighbours(s);
  // surface silanols: an oxygen bonded to exactly one silicon and one hydrogen
  std::vector<uint32_t> sites;
  for (uint32_t o_ = 0; o_ < s.atoms.size(); ++o_) {
    if (s.atoms[o_].element != 8 || nb[o_].size() != 2) continue;
    int si = 0, h = 0;
    for (uint32_t q : nb[o_]) si += s.atoms[q].element == 14, h += s.atoms[q].element == 1;
    if (si == 1 && h == 1) sites.push_back(o_);
  }
  rep.silanols = sites.size();
  if (sites.empty()) throw EditError("no surface silanols (Si–O–H): build the silica passivated (Nanostructure or Surface builder)");
  const size_t want = o.count > 0 ? size_t(o.count) : size_t(std::lround(o.fraction * double(sites.size())));
  if (want == 0) throw EditError("nothing to graft: a count or a fraction above zero");
  // random order, then greedy: each accepted site at least min_spacing from those before it
  std::mt19937_64 rng(o.seed);
  std::shuffle(sites.begin(), sites.end(), rng);
  std::vector<Vec3> kept;
  std::vector<uint32_t> chosen;
  for (uint32_t site : sites) {
    if (chosen.size() >= want) break;
    bool ok = true;
    for (const auto& p : kept) {
      Vec3 d = s.atoms[site].pos - p;
      if (s.cell.valid()) d = s.cell.minimum_image(d);
      if (norm(d) < o.min_spacing) { ok = false; break; }
    }
    if (!ok) continue;
    kept.push_back(s.atoms[site].pos);
    chosen.push_back(site);
  }
  // attach from the highest index down: deleting each replaced hydrogen shifts only the indices above it
  std::sort(chosen.rbegin(), chosen.rend());
  for (uint32_t site : chosen) {
    rep.added_atoms += attach_fragment(s, site, o.smiles, 0, true).size();
    ++rep.grafted;
  }
  if (rep.grafted < want)
    rep.notes.push_back(std::to_string(want - rep.grafted) + " fewer than asked: silanols closer than " + std::to_string(o.min_spacing).substr(0, 4) +
                        " Å to a grafted one were passed over");
  rep.notes.push_back(o.name + " grafted on " + std::to_string(rep.grafted) + " of " + std::to_string(rep.silanols) + " silanols (" +
                      std::to_string(rep.added_atoms) + " atoms added; one ethanol released per graft); relax before dynamics");
  return rep;
}

std::vector<uint32_t> attach_fragment(System& s, uint32_t target, const std::string& smiles, int which, bool replace_h) {
  if (target >= s.atoms.size()) throw EditError("pick the atom to attach to");
  const auto& F = fragment_3d(smiles);
  if (which < 0 || size_t(which) >= F.dummy.size()) which = 0;
  const uint32_t D = F.dummy[size_t(which)], X = F.root[size_t(which)];
  const auto nb = neighbours(s);
  // where the fragment's attaching atom goes: in place of one of the target's hydrogens, or in its free direction
  int h = -1;
  if (replace_h)
    for (uint32_t q : nb[target]) if (s.atoms[q].element == 1) { h = int(q); break; }
  Vec3 dir;
  if (h >= 0) dir = unitv(rel(s, target, uint32_t(h)));
  else {
    std::vector<Vec3> dirs;
    for (uint32_t q : nb[target]) dirs.push_back(unitv(rel(s, target, q)));
    const auto cand = ideal_directions(dirs, std::max(4, int(dirs.size()) + 1), dirs.empty() ? Vec3{0, 0, 1} : perpendicular(dirs[0]));
    dir = cand.empty() ? Vec3{1, 0, 0} : cand[0];
  }
  const Vec3 T = s.atoms[target].pos;
  const double len = element(F.s.atoms[X].element).covalent + element(s.atoms[target].element).covalent;
  const Vec3 P = T + dir * len;
  // rotate the fragment so its attachment point (the dummy) points back at the target
  const Vec3 v = unitv(F.s.atoms[D].pos - F.s.atoms[X].pos), w = unitv(T - P);
  auto rotate = [](const Vec3& p, const Vec3& axis, double c, double sn) { return p * c + cross(axis, p) * sn + axis * (dot(axis, p) * (1 - c)); };
  const Vec3 axis0 = cross(v, w);
  const double s0 = norm(axis0), c0 = dot(v, w);
  auto align = [&](const Vec3& p) {
    if (s0 < 1e-9) return c0 > 0 ? p : rotate(p, perpendicular(v), -1, 0);
    return rotate(p, axis0 * (1.0 / s0), c0, s0);
  };
  // roll about the new bond: the orientation farthest from the structure
  double best = -1;
  std::vector<Vec3> placed;
  for (int k = 0; k < 24; ++k) {
    const double ang = 2 * 3.14159265358979323846 * k / 24;
    std::vector<Vec3> trial(F.s.atoms.size());
    double dmin = 1e300;
    for (size_t i = 0; i < F.s.atoms.size(); ++i) {
      const Vec3 local = align(F.s.atoms[i].pos - F.s.atoms[X].pos);
      trial[i] = P + rotate(local, w, std::cos(ang), std::sin(ang));
      if (i == D) continue;
      for (size_t j = 0; j < s.atoms.size(); ++j) {
        if (j == target || int(j) == h) continue;
        Vec3 d = trial[i] - s.atoms[j].pos;
        if (s.cell.valid()) d = s.cell.minimum_image(d);
        dmin = std::min(dmin, norm(d));
      }
    }
    if (dmin > best) best = dmin, placed = trial;
  }
  // into the structure: the fragment's atoms (not the dummy), the bond, the other attachment points as hydrogens
  std::vector<int64_t> map(F.s.atoms.size(), -1);
  std::vector<uint32_t> added;
  for (size_t i = 0; i < F.s.atoms.size(); ++i) {
    if (i == D) continue;
    Atom a = F.s.atoms[i];
    const bool other_dummy = std::find(F.dummy.begin(), F.dummy.end(), uint32_t(i)) != F.dummy.end();
    a.element = other_dummy ? 1 : a.element;
    a.pos = placed[i];
    if (other_dummy) {
      // a hydrogen at its own length along the dummy's bond
      uint32_t r = 0;
      for (const auto& b : F.s.bonds) { if (b.i == i) r = b.j; if (b.j == i) r = b.i; }
      a.pos = placed[r] + unitv(placed[i] - placed[r]) * h_length(F.s.atoms[r].element);
    }
    a.type = type_for(s, a.element);
    a.mol = s.atoms[target].mol;
    a.resname = s.atoms[target].resname;
    a.resid = s.atoms[target].resid;
    a.id = s.atoms.back().id + 1;
    a.charge = 0;
    a.name = std::string(element(a.element).symbol) + std::to_string(s.atoms.size() + 1);
    map[i] = int64_t(s.atoms.size());
    added.push_back(uint32_t(s.atoms.size()));
    s.atoms.push_back(a);
  }
  for (const auto& b : F.s.bonds)
    if (map[b.i] >= 0 && map[b.j] >= 0) s.bonds.push_back({uint32_t(map[b.i]), uint32_t(map[b.j]), b.order});
  s.bonds.push_back({target, uint32_t(map[X]), 1});
  if (h >= 0) {
    std::vector<char> del(s.atoms.size(), 0);
    del[size_t(h)] = 1;
    delete_atoms(s, del);
    for (auto& a : added) if (a > uint32_t(h)) --a;
  }
  return added;
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

std::vector<uint32_t> fuse_benzene(System& s, uint32_t i, uint32_t j) {
  if (i >= s.atoms.size() || j >= s.atoms.size() || i == j) throw EditError("pick the two atoms of a bond");
  if (!std::any_of(s.bonds.begin(), s.bonds.end(), [&](const Bond& b) { return (b.i == i && b.j == j) || (b.i == j && b.j == i); }))
    throw EditError("the two picked atoms are not bonded");
  const auto nb = neighbours(s);
  // the hydrogens of i and j that point to the same side
  int hi = -1, hj = -1;
  double best = -2;
  for (uint32_t a : nb[i]) {
    if (s.atoms[a].element != 1) continue;
    for (uint32_t b : nb[j]) {
      if (s.atoms[b].element != 1) continue;
      const double c = dot(unitv(rel(s, i, a)), unitv(rel(s, j, b)));
      if (c > best) best = c, hi = int(a), hj = int(b);
    }
  }
  if (hi < 0) throw EditError("both atoms need a hydrogen for the ring to replace");
  const Vec3 bv = rel(s, i, j);
  const double L = norm(bv);
  const Vec3 e1 = unitv(bv);
  Vec3 u = rel(s, i, uint32_t(hi)) * (1 / norm(rel(s, i, uint32_t(hi)))) + rel(s, j, uint32_t(hj)) * (1 / norm(rel(s, j, uint32_t(hj))));
  u = u - e1 * dot(u, e1);
  if (norm(u) < 1e-3) throw EditError("the hydrogens point along the bond: no side to fuse the ring on");
  u = unitv(u);
  const double h = std::sqrt(3.0) / 2 * L;
  const Vec3 o = s.atoms[i].pos;
  const Vec3 ring[4] = {o + e1 * (1.5 * L) + u * h, o + e1 * L + u * (2 * h), o + u * (2 * h), o + e1 * (-0.5 * L) + u * h};
  const Vec3 centre = o + e1 * (0.5 * L) + u * h;
  const uint32_t first = uint32_t(s.atoms.size());
  auto put = [&](int z, const Vec3& p) {
    Atom a;
    a.element = z;
    a.id = s.atoms.back().id + 1;
    a.type = type_for(s, z);
    a.name = std::string(element(z).symbol) + std::to_string(s.atoms.size() + 1);
    a.pos = p;
    a.mol = s.atoms[i].mol, a.resname = s.atoms[i].resname, a.resid = s.atoms[i].resid;
    s.atoms.push_back(a);
    return uint32_t(s.atoms.size() - 1);
  };
  uint32_t c[4];
  for (int k = 0; k < 4; ++k) c[k] = put(6, ring[k]);
  for (int k = 0; k < 4; ++k) {
    const uint32_t hh = put(1, ring[k] + unitv(ring[k] - centre) * 1.08);
    s.bonds.push_back({c[k], hh, 1});
  }
  s.bonds.push_back({j, c[0], 4});
  for (int k = 0; k < 3; ++k) s.bonds.push_back({c[k], c[k + 1], 4});
  s.bonds.push_back({c[3], i, 4});
  for (auto& b : s.bonds)
    if ((b.i == i && b.j == j) || (b.i == j && b.j == i)) b.order = 4;
  std::vector<char> gone(s.atoms.size(), 0);
  gone[size_t(hi)] = gone[size_t(hj)] = 1;
  delete_atoms(s, gone);
  // the new atoms moved down by the removed hydrogens that came before them (both did)
  std::vector<uint32_t> out;
  for (uint32_t k = first; k < first + 8; ++k) out.push_back(k - 2);
  return out;
}

int protonate_residues(System& s, double ph, std::vector<std::string>* notes) {
  const size_t n = s.atoms.size();
  const auto nb = neighbours(s);
  auto order = [&](uint32_t a, uint32_t b) {
    for (const auto& bd : s.bonds)
      if ((bd.i == a && bd.j == b) || (bd.i == b && bd.j == a)) return bd.order;
    return 0;
  };
  auto heavy_nb = [&](uint32_t a) {
    std::vector<uint32_t> v;
    for (uint32_t w : nb[a]) if (s.atoms[w].element != 1) v.push_back(w);
    return v;
  };
  auto upper = [](std::string r) { for (auto& c : r) c = char(std::toupper(static_cast<unsigned char>(c))); return r.substr(0, 3); };
  // the chains' first and last residues (termini), by molecule
  std::map<int64_t, std::pair<int64_t, int64_t>> span;
  for (const auto& a : s.atoms) {
    if (a.resid == 0) continue;
    auto it = span.find(a.mol);
    if (it == span.end()) span[a.mol] = {a.resid, a.resid};
    else it->second.first = std::min(it->second.first, a.resid), it->second.second = std::max(it->second.second, a.resid);
  }
  // what an earlier pH set on residue N, O and S atoms goes first
  for (auto& a : s.atoms)
    if (!a.resname.empty() && (a.element == 7 || a.element == 8 || a.element == 16) && std::fabs(a.charge - std::lround(a.charge)) < 1e-6) a.charge = 0;
  std::map<std::string, int> count;
  std::set<std::pair<int64_t, int64_t>> his_done;
  int sites = 0;
  auto set = [&](uint32_t a, int q, const std::string& what) {
    s.atoms[a].charge = q;
    ++count[what];
    ++sites;
  };
  // a carboxylate: its single-bonded oxygen (the other is C=O); none when both are double or both single with H
  auto carboxyl_o = [&](uint32_t c) -> int {
    int single = -1, dbl = 0;
    for (uint32_t w : heavy_nb(c))
      if (s.atoms[w].element == 8 && heavy_nb(w).size() == 1) (order(c, w) == 2 ? dbl : single) = order(c, w) == 2 ? dbl + 1 : int(w);
    return dbl == 1 ? single : -1;
  };
  for (uint32_t i = 0; i < n; ++i) {
    const Atom& a = s.atoms[i];
    if (a.element == 1 || a.resname.empty()) continue;
    const std::string r = upper(a.resname);
    const auto h = heavy_nb(i);
    const bool first = span.count(a.mol) && a.resid == span[a.mol].first, last = span.count(a.mol) && a.resid == span[a.mol].second;
    const bool backbone_n = a.element == 7 && (a.name == "N" || std::any_of(h.begin(), h.end(), [&](uint32_t w) { return s.atoms[w].name == "CA"; }));
    if (a.element == 7 && backbone_n && h.size() == 1 && first) {   // N-terminus
      if (ph < 8.0) set(i, 1, "N-terminus NH3+");
      continue;
    }
    if (a.element == 6 && (a.name == "C" || last) && carboxyl_o(i) >= 0 && std::any_of(h.begin(), h.end(), [&](uint32_t w) { return s.atoms[w].name == "CA"; })) {
      if (ph > 3.1) set(uint32_t(carboxyl_o(i)), -1, "C-terminus COO−");   // C-terminus
      continue;
    }
    if (backbone_n || a.name == "C" || a.name == "O" || a.name == "CA") continue;
    if ((r == "ASP" || r == "GLU") && a.element == 6 && carboxyl_o(i) >= 0) {
      if (ph > (r == "ASP" ? 3.9 : 4.3)) set(uint32_t(carboxyl_o(i)), -1, r == "ASP" ? "Asp−" : "Glu−");
    } else if (r == "LYS" && a.element == 7 && h.size() == 1) {
      if (ph < 10.5) set(i, 1, "Lys+");
    } else if (r == "ARG" && a.element == 7) {   // the guanidine N double-bonded to its central carbon
      for (uint32_t w : h)
        if (s.atoms[w].element == 6 && order(i, w) == 2 && std::count_if(nb[w].begin(), nb[w].end(), [&](uint32_t x) { return s.atoms[x].element == 7; }) == 3 && ph < 12.5)
          set(i, 1, "Arg+");
    } else if ((r == "HIS" || r == "HID" || r == "HIE" || r == "HIP") && a.element == 7 && h.size() == 2) {
      // His+: the ring N without its hydrogen (a double or aromatic bond, no H yet) takes the charge, once per residue
      const double sum = (order(i, h[0]) == 4 ? 1.5 : order(i, h[0])) + (order(i, h[1]) == 4 ? 1.5 : order(i, h[1]));
      const bool has_h = nb[i].size() > h.size();
      if (ph < 6.0 && sum >= 2.5 && !has_h && his_done.insert({a.mol, a.resid}).second) set(i, 1, "His+");
    } else if (r == "CYS" && a.element == 16 && h.size() == 1) {
      if (ph > 8.3) set(i, -1, "Cys−");
    } else if (r == "TYR" && a.element == 8 && h.size() == 1 && s.atoms[h[0]].element == 6) {
      bool aromatic = false;
      for (const auto& bd : s.bonds) aromatic = aromatic || ((bd.i == h[0] || bd.j == h[0]) && bd.order == 4);
      if (aromatic && ph > 10.1) set(i, -1, "Tyr−");
    }
  }
  if (notes) {
    std::string t;
    int q = 0;
    for (const auto& a : s.atoms) q += int(std::lround(a.charge));
    for (const auto& [k, c] : count) t += (t.empty() ? "" : ", ") + std::to_string(c) + " " + k;
    char b[96];
    std::snprintf(b, sizeof b, "pH %.1f (model pKa values): ", ph);
    notes->push_back(std::string(b) + (t.empty() ? "no charged residues" : t) + " · net formal charge " + std::to_string(q));
  }
  return sites;
}

int add_hydrogens_at_ph(System& s, double ph, const std::vector<char>& atoms, std::vector<std::string>* notes) {
  protonate_residues(s, ph, notes);
  int added = add_hydrogens(s, atoms);
  // histidine rings: which N carries H
  auto upper = [](std::string r) { for (auto& c : r) c = char(std::toupper(static_cast<unsigned char>(c))); return r.substr(0, 3); };
  std::map<std::pair<int64_t, int64_t>, std::vector<uint32_t>> ring_n;
  auto nb = neighbours(s);
  for (uint32_t i = 0; i < s.atoms.size(); ++i) {
    const auto& a = s.atoms[i];
    if (a.element != 7 || a.resname.empty()) continue;
    const std::string r = upper(a.resname);
    if (r != "HIS" && r != "HID" && r != "HIE" && r != "HIP") continue;
    int heavy = 0;
    bool backbone = a.name == "N";
    for (uint32_t w : nb[i]) heavy += s.atoms[w].element != 1, backbone = backbone || s.atoms[w].name == "CA" || s.atoms[w].name == "C";
    if (heavy == 2 && !backbone && !(atoms.size() == s.atoms.size() && !atoms[i])) ring_n[{a.mol, a.resid}].push_back(i);
  }
  for (auto& [key, ns] : ring_n) {
    if (ns.size() != 2) continue;
    auto has_h = [&](uint32_t i) { return std::any_of(nb[i].begin(), nb[i].end(), [&](uint32_t w) { return s.atoms[w].element == 1; }); };
    const bool charged = std::any_of(ns.begin(), ns.end(), [&](uint32_t i) { return s.atoms[i].charge > 0.5; });
    std::vector<uint32_t> want;
    if (charged) want = ns;
    else if (!has_h(ns[0]) && !has_h(ns[1])) {
      // HIE: the N farther (in bonds) from the residue's backbone CA, or from any non-ring heavy atom
      auto dist_to_ca = [&](uint32_t from) {
        std::map<uint32_t, int> d{{from, 0}};
        std::vector<uint32_t> q{from};
        for (size_t k = 0; k < q.size(); ++k)
          for (uint32_t w : nb[q[k]])
            if (d.emplace(w, d[q[k]] + 1).second) {
              if (s.atoms[w].name == "CA" && s.atoms[w].resid == key.second) return d[w];
              q.push_back(w);
            }
        return 0;
      };
      want.push_back(dist_to_ca(ns[0]) >= dist_to_ca(ns[1]) ? ns[0] : ns[1]);
    }
    for (uint32_t i : want)
      if (!has_h(i)) { add_atom(s, int(i), 1); ++added; nb = neighbours(s); }
  }
  return added;
}

int hydroxylate_phosphorus(System& s) {
  int changed = 0;
  const auto nb = neighbours(s);
  std::vector<std::pair<uint32_t, uint32_t>> ph;   // (P, H)
  for (uint32_t i = 0; i < s.atoms.size(); ++i)
    if (s.atoms[i].element == 15)
      for (uint32_t w : nb[i]) if (s.atoms[w].element == 1) ph.push_back({i, w});
  for (const auto& [p, h] : ph) {
    const Vec3 d = unitv(rel(s, p, h));
    Atom& o = s.atoms[h];
    o.element = 8;
    o.type = type_for(s, 8);
    o.name = "O" + std::to_string(h + 1);
    o.pos = s.atoms[p].pos + d * 1.61;
    add_atom(s, int(h), 1);
    ++changed;
  }
  return changed;
}

}  // namespace caps
