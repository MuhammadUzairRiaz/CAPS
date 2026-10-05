// CAPS structure editing (see caps/edit.hpp).
#include "caps/edit.hpp"

#include <algorithm>
#include <random>
#include <cmath>
#include <functional>
#include <map>
#include <set>
#include <sstream>

#include "caps/analysis.hpp"
#include "caps/appearance.hpp"
#include "caps/elements.hpp"
#include "caps/molecule.hpp"
#include "caps/relax.hpp"
#include "caps/torsion.hpp"
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
    if (b.i == a || b.j == a) v += b.order == kBondDative ? 0 : b.order == 2 ? 2 : b.order == 3 ? 3 : b.order == 4 ? 1.5 : 1;
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
    const double o = b.order == kBondDative ? 0 : b.order == 2 ? 2 : b.order == 3 ? 3 : b.order == 4 ? 1.5 : 1;   // coordinate: no valence
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

namespace {
Vec3 rotate_about(const Vec3& p, const Vec3& origin, const Vec3& axis, double rad) {
  const Vec3 r = p - origin, k = unitv(axis);
  const double c = std::cos(rad), sn = std::sin(rad);
  return origin + r * c + cross(k, r) * sn + k * (dot(k, r) * (1 - c));
}
std::vector<uint32_t> side_of(const System& s, uint32_t b, uint32_t c) {
  bool bonded = false;
  for (const auto& bd : s.bonds) bonded |= (bd.i == b && bd.j == c) || (bd.i == c && bd.j == b);
  if (!bonded) throw EditError("atoms " + std::to_string(b + 1) + " and " + std::to_string(c + 1) + " are not bonded");
  try {
    const auto m = moving_side(s, int(b), int(c));
    return m;
  } catch (const std::exception&) {
    throw EditError("the bond " + std::to_string(b + 1) + "–" + std::to_string(c + 1) + " is in a ring: its two sides cannot move apart");
  }
}
}  // namespace

void set_bond_length(System& s, uint32_t i, uint32_t j, double r) {
  if (i >= s.atoms.size() || j >= s.atoms.size() || i == j) throw EditError("pick two bonded atoms");
  if (!(r > 0.3 && r < 10)) throw EditError("a bond length between 0.3 and 10 Å");
  const auto mv = side_of(s, i, j);
  const Vec3 d = s.atoms[j].pos - s.atoms[i].pos;
  const double now = norm(d);
  if (now < 1e-9) throw EditError("the two atoms sit on each other");
  const Vec3 shift = d * ((r - now) / now);
  for (uint32_t a : mv) s.atoms[a].pos = s.atoms[a].pos + shift;
}

void set_bond_angle(System& s, uint32_t i, uint32_t j, uint32_t k, double theta_deg) {
  if (i >= s.atoms.size() || j >= s.atoms.size() || k >= s.atoms.size() || i == k || i == j || j == k) throw EditError("pick three atoms i–j–k");
  if (!(theta_deg > 1 && theta_deg < 179.5)) throw EditError("an angle between 1° and 179.5°");
  const auto mv = side_of(s, j, k);
  if (std::find(mv.begin(), mv.end(), i) != mv.end()) throw EditError("atom " + std::to_string(i + 1) + " is on the side that moves (a ring)");
  const Vec3 a = s.atoms[i].pos - s.atoms[j].pos, b = s.atoms[k].pos - s.atoms[j].pos;
  Vec3 n = cross(a, b);
  if (norm(n) < 1e-9) n = perpendicular(a);   // straight now: any plane through the bond
  const double now = std::acos(std::clamp(dot(a, b) / (norm(a) * norm(b)), -1.0, 1.0));
  const double rot = theta_deg * M_PI / 180 - now;
  for (uint32_t x : mv) s.atoms[x].pos = rotate_about(s.atoms[x].pos, s.atoms[j].pos, n, rot);
}

void set_torsion(System& s, uint32_t i, uint32_t j, uint32_t k, uint32_t l, double phi_deg) {
  if (std::max({i, j, k, l}) >= s.atoms.size()) throw EditError("pick four atoms i–j–k–l");
  const auto mv = side_of(s, j, k);
  if (std::find(mv.begin(), mv.end(), i) != mv.end()) throw EditError("atom " + std::to_string(i + 1) + " is on the side that moves (a ring)");
  set_dihedral(s, {int(i), int(j), int(k), int(l)}, phi_deg, mv);
}

void rotate_atoms(System& s, const std::vector<uint32_t>& atoms, const Vec3& axis, double degrees) {
  if (atoms.empty()) throw EditError("pick or select the atoms to rotate");
  if (norm(axis) < 1e-12) throw EditError("the rotation axis has no length");
  Vec3 c{0, 0, 0};
  for (uint32_t a : atoms) c = c + s.atoms[a].pos;
  c = c * (1.0 / double(atoms.size()));
  for (uint32_t a : atoms) s.atoms[a].pos = rotate_about(s.atoms[a].pos, c, axis, degrees * M_PI / 180);
}

std::vector<Vec3> coordination_directions(const std::string& g) {
  const double r3 = 1 / std::sqrt(3.0), c120 = -0.5, s120 = std::sqrt(3.0) / 2;
  if (g == "linear") return {{0, 0, 1}, {0, 0, -1}};
  if (g == "trigonal" || g == "trigonal_planar") return {{1, 0, 0}, {c120, s120, 0}, {c120, -s120, 0}};
  if (g == "tetrahedral") return {{r3, r3, r3}, {r3, -r3, -r3}, {-r3, r3, -r3}, {-r3, -r3, r3}};
  if (g == "square_planar") return {{1, 0, 0}, {0, 1, 0}, {-1, 0, 0}, {0, -1, 0}};
  if (g == "trigonal_bipyramidal") return {{0, 0, 1}, {0, 0, -1}, {1, 0, 0}, {c120, s120, 0}, {c120, -s120, 0}};
  if (g == "square_pyramidal") return {{0, 0, 1}, {1, 0, 0}, {0, 1, 0}, {-1, 0, 0}, {0, -1, 0}};
  if (g == "octahedral") return {{1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}};
  return {};
}

namespace {
// the rotation taking unit a to unit b (Rodrigues), applied to v
Vec3 turn(const Vec3& v, const Vec3& a, const Vec3& b) {
  Vec3 k = cross(a, b);
  const double sn = norm(k), cs = dot(a, b);
  if (sn < 1e-12) {
    if (cs > 0) return v;
    k = unitv(perpendicular(a));   // opposite: half a turn about any perpendicular
    return v * -1.0 + k * (2 * dot(k, v));
  }
  k = k * (1 / sn);
  return v * cs + cross(k, v) * sn + k * (dot(k, v) * (1 - cs));
}
}  // namespace

double set_coordination(System& s, uint32_t centre, const std::string& geometry) {
  if (centre >= s.atoms.size()) throw EditError("pick the centre atom");
  const std::vector<Vec3> ideal = coordination_directions(geometry);
  if (ideal.empty())
    throw EditError("geometry " + geometry + ": linear, trigonal, tetrahedral, square_planar, trigonal_bipyramidal, square_pyramidal or octahedral");
  std::vector<uint32_t> lig;
  for (const auto& b : s.bonds)
    if (b.i == centre) lig.push_back(b.j);
    else if (b.j == centre) lig.push_back(b.i);
  if (lig.empty()) throw EditError("atom " + std::to_string(centre + 1) + " has no bonded neighbours to place");
  if (lig.size() > ideal.size())
    throw EditError("atom " + std::to_string(centre + 1) + " has " + std::to_string(lig.size()) + " neighbours; " + geometry + " has " + std::to_string(ideal.size()) + " sites");
  const Vec3 c = s.atoms[centre].pos;
  std::vector<Vec3> u;
  for (uint32_t l : lig) {
    const Vec3 d = s.cell.valid() ? s.cell.minimum_image(s.atoms[l].pos - c) : s.atoms[l].pos - c;
    if (norm(d) < 1e-9) throw EditError("a neighbour sits on the centre");
    u.push_back(unitv(d));
  }
  // the ideal set turned so its site p lies on ligand 0 and site q in the plane of ligands 0 and 1; the rest assigned
  // greedily to the nearest free site; the turn and assignment that move the ligands least win
  double best = 1e300;
  std::vector<Vec3> best_dir;
  for (size_t p = 0; p < ideal.size(); ++p)
    for (size_t q = 0; q < ideal.size(); ++q) {
      if (q == p && ideal.size() > 1) continue;
      auto R = [&](const Vec3& v) {   // v in the turned frame
        Vec3 w = turn(v, ideal[p], u[0]);
        if (u.size() < 2) return w;
        // then about u[0] so that site q's image lies in the plane of u[0], u[1] on u[1]'s side
        const Vec3 qa = turn(ideal[q], ideal[p], u[0]);
        Vec3 a = qa - u[0] * dot(qa, u[0]), b = u[1] - u[0] * dot(u[1], u[0]);
        if (norm(a) < 1e-9 || norm(b) < 1e-9) return w;
        a = unitv(a), b = unitv(b);
        const double ang = std::atan2(dot(cross(a, b), u[0]), dot(a, b));
        return rotate_about(w, {0, 0, 0}, u[0], ang);
      };
      std::vector<Vec3> sites;
      for (const auto& v : ideal) sites.push_back(R(v));
      std::vector<char> used(sites.size(), 0);
      std::vector<Vec3> dir(u.size());
      double cost = 0;
      for (size_t k = 0; k < u.size(); ++k) {
        size_t bk = sites.size();
        double bd = -2;
        for (size_t m = 0; m < sites.size(); ++m)
          if (!used[m] && dot(sites[m], u[k]) > bd) bd = dot(sites[m], u[k]), bk = m;
        used[bk] = 1;
        dir[k] = sites[bk];
        cost += std::acos(std::clamp(bd, -1.0, 1.0));
      }
      if (cost < best) best = cost, best_dir = dir;
    }
  // each ligand turned about the centre onto its site, with its substituents (a chelate's ligand alone)
  double moved = 0;
  for (size_t k = 0; k < lig.size(); ++k) {
    auto side = side_of(s, centre, lig[k]);
    bool shared = std::find(side.begin(), side.end(), centre) != side.end();
    for (size_t m = 0; m < lig.size() && !shared; ++m) shared = m != k && std::find(side.begin(), side.end(), lig[m]) != side.end();
    if (shared) side = {lig[k]};
    moved = std::max(moved, std::acos(std::clamp(dot(u[k], best_dir[k]), -1.0, 1.0)) * 180 / M_PI);
    for (uint32_t a : side) {
      const Vec3 d = s.cell.valid() ? s.cell.minimum_image(s.atoms[a].pos - c) : s.atoms[a].pos - c;
      s.atoms[a].pos = c + turn(d, u[k], best_dir[k]);
    }
  }
  return moved;
}

void mirror_atoms(System& s, const std::vector<uint32_t>& atoms, const Vec3& normal) {
  if (atoms.empty()) throw EditError("pick or select the atoms to mirror");
  if (norm(normal) < 1e-12) throw EditError("the mirror plane's normal has no length");
  const Vec3 nn = unitv(normal);
  Vec3 c{0, 0, 0};
  for (uint32_t a : atoms) c = c + s.atoms[a].pos;
  c = c * (1.0 / double(atoms.size()));
  for (uint32_t a : atoms) s.atoms[a].pos = s.atoms[a].pos - nn * (2 * dot(s.atoms[a].pos - c, nn));
}

bool set_configuration(System& s, uint32_t centre, const std::string& rs) {
  if (rs != "R" && rs != "S") throw EditError("the configuration is R or S");
  if (centre >= s.atoms.size()) throw EditError("pick the stereocentre");
  std::vector<char> only(s.atoms.size(), 0);
  only[centre] = 1;
  const auto now = stereo_labels(s, only)[centre];
  if (now.empty()) return false;
  if (now != rs) invert_centre(s, centre);
  return true;
}

std::string thiolate_smiles(const std::string& name) {
  if (name == "C6" || name == "hexanethiolate") return "*SCCCCCC";
  if (name == "C12" || name == "dodecanethiolate") return "*SCCCCCCCCCCCC";
  if (name == "C18" || name == "octadecanethiolate") return "*SCCCCCCCCCCCCCCCCCC";
  if (name == "MPA") return "*SCCC(=O)O";
  if (name == "MUA") return "*SCCCCCCCCCCC(=O)O";
  if (name == "MHA") return "*SCCCCCCO";
  throw EditError("unknown thiolate '" + name + "' (C6, C12, C18, MPA, MUA, MHA, or a SMILES starting *S)");
}

ThiolateReport cap_thiolates(System& s, const ThiolateOptions& o) {
  ThiolateReport rep;
  if (o.smiles.find('*') == std::string::npos) throw EditError("the thiolate SMILES needs * where the sulfur binds the metal");
  auto metal = [](int z) { return z == 79 || z == 47 || z == 29 || z == 78 || z == 46; };
  std::vector<uint32_t> M;
  for (uint32_t i = 0; i < s.atoms.size(); ++i)
    if (metal(s.atoms[i].element)) M.push_back(i);
  if (M.size() < 4) throw EditError("no metal surface to cap: thiolates bind gold, silver, copper, platinum or palladium");
  auto rel = [&](const Vec3& a, const Vec3& b) { Vec3 d = b - a; return s.cell.valid() ? s.cell.minimum_image(d) : d; };
  // nearest-neighbour distance and the neighbour lists within 1.2 of it
  double dnn = 1e300;
  for (size_t a = 0; a < std::min<size_t>(M.size(), 200); ++a)
    for (size_t b = 0; b < M.size(); ++b)
      if (a != b) dnn = std::min(dnn, norm(rel(s.atoms[M[a]].pos, s.atoms[M[b]].pos)));
  const double cut = 1.2 * dnn;
  std::vector<std::vector<uint32_t>> nbm(M.size());
  for (size_t a = 0; a < M.size(); ++a)
    for (size_t b = a + 1; b < M.size(); ++b)
      if (norm(rel(s.atoms[M[a]].pos, s.atoms[M[b]].pos)) < cut) nbm[a].push_back(uint32_t(b)), nbm[b].push_back(uint32_t(a));
  size_t zmax = 0;
  for (const auto& l : nbm) zmax = std::max(zmax, l.size());
  // surface atoms: fewer neighbours than the bulk; their outward direction away from their neighbours' centre
  std::vector<char> surf(M.size(), 0);
  std::vector<Vec3> out(M.size(), Vec3{0, 0, 0});
  for (size_t a = 0; a < M.size(); ++a) {
    if (nbm[a].size() >= zmax) continue;
    Vec3 c{0, 0, 0};
    for (uint32_t b : nbm[a]) c = c + rel(s.atoms[M[a]].pos, s.atoms[M[b]].pos);
    const double l = norm(c);
    if (l < 1e-6) continue;
    surf[a] = 1;
    out[a] = c * (-1.0 / l);
    ++rep.surface_atoms;
  }
  if (!rep.surface_atoms) throw EditError("the metal has no surface atoms (a periodic bulk crystal?)");
  // candidate sites: each surface triangle (three mutually neighbouring surface atoms) facing out, else on top
  struct Site { Vec3 p, n; bool hollow; };
  std::vector<Site> sites;
  const double dms = 2.45;   // Å, metal–S
  for (size_t a = 0; a < M.size(); ++a) {
    if (!surf[a]) continue;
    for (uint32_t b : nbm[a]) {
      if (!surf[b] || b <= a) continue;
      for (uint32_t c : nbm[b]) {
        if (!surf[c] || c <= b || std::find(nbm[a].begin(), nbm[a].end(), c) == nbm[a].end()) continue;
        const Vec3 A = s.atoms[M[a]].pos, B = A + rel(A, s.atoms[M[b]].pos), C = A + rel(A, s.atoms[M[c]].pos);
        const Vec3 g = (A + B + C) * (1.0 / 3.0);
        Vec3 n = cross(B - A, C - A);
        const double ln = norm(n);
        if (ln < 1e-9) continue;
        n = n * (1 / ln);
        const Vec3 o3 = out[a] + out[b] + out[c];
        if (dot(n, o3) < 0) n = n * -1.0;
        if (dot(n, unitv(o3)) < 0.5) continue;   // an edge triangle standing across the surface
        // an atom below the hollow's centre (hcp site) or not (fcc) — either is a three-fold hollow
        const double R = norm(A - g), h = std::sqrt(std::max(0.25, dms * dms - R * R));
        sites.push_back({g + n * h, n, true});
      }
    }
  }
  const size_t hollows = sites.size();
  for (size_t a = 0; a < M.size(); ++a)
    if (surf[a]) sites.push_back({s.atoms[M[a]].pos + out[a] * dms, out[a], false});
  // greedy in random order, hollows first: each S at least the spacing from those before, and clear of the metal
  const double spacing = o.min_spacing > 0 ? o.min_spacing : std::sqrt(3.0) * dnn;
  std::mt19937_64 rng(o.seed);
  std::shuffle(sites.begin(), sites.begin() + std::ptrdiff_t(hollows), rng);
  std::shuffle(sites.begin() + std::ptrdiff_t(hollows), sites.end(), rng);
  std::vector<Site> chosen;
  for (const auto& st : sites) {
    bool ok = true;
    for (const auto& c : chosen)
      if (norm(rel(c.p, st.p)) < spacing * 0.999) { ok = false; break; }
    if (!ok) continue;
    for (uint32_t m : M)
      if (norm(rel(st.p, s.atoms[m].pos)) < dms - 0.05) { ok = false; break; }
    if (ok) chosen.push_back(st);
  }
  const size_t want = size_t(std::lround(std::clamp(o.fraction, 0.0, 1.0) * double(chosen.size())));
  if (want == 0) throw EditError("no site for a thiolate at this spacing");
  std::shuffle(chosen.begin(), chosen.end(), rng);
  chosen.resize(want);
  // each ligand: its S on the site, the tail's heavy-atom centre along the surface normal, rolled clear of the others
  const auto& F = fragment_3d(o.smiles);
  const uint32_t D = F.dummy[0], X = F.root[0];
  if (F.s.atoms[X].element != 16) throw EditError("the thiolate must bind through sulfur: write its SMILES as *S…");
  // the ligand straightened: every torsion along its longest heavy-atom path from the metal side anti (all-trans tail)
  System L = F.s;
  {
    const auto nbl = neighbours(L);
    std::vector<int> path{int(D), int(X)}, best_path;
    std::function<void(std::vector<int>&)> walk = [&](std::vector<int>& p) {
      bool ext = false;
      for (uint32_t q : nbl[size_t(p.back())]) {
        if (L.atoms[q].element <= 1 || std::find(p.begin(), p.end(), int(q)) != p.end()) continue;
        p.push_back(int(q));
        walk(p);
        p.pop_back();
        ext = true;
      }
      if (!ext && p.size() > best_path.size()) best_path = p;
    };
    walk(path);
    for (size_t k = 1; k + 2 < best_path.size(); ++k) {
      try {
        const auto mv = moving_side(L, best_path[k], best_path[k + 1]);
        set_dihedral(L, {best_path[k - 1], best_path[k], best_path[k + 1], best_path[k + 2]}, 180.0, mv);
      } catch (const std::exception&) {}   // a ring bond stays as built
    }
  }
  // the S–C bond along the surface normal (tilted up to 30° where that clears the neighbours), rolled about it
  int C1 = -1;
  for (const auto& b : L.bonds) {
    if (b.i == X && b.j != D) C1 = int(b.j);
    if (b.j == X && b.i != D) C1 = int(b.i);
  }
  const Vec3 v = C1 >= 0 ? unitv(L.atoms[size_t(C1)].pos - L.atoms[X].pos) : unitv(L.atoms[X].pos - L.atoms[D].pos);
  double reach = 0;
  for (const auto& a : L.atoms) reach = std::max(reach, norm(a.pos - L.atoms[X].pos));
  auto rotate = [](const Vec3& p, const Vec3& axis, double c, double sn) { return p * c + cross(axis, p) * sn + axis * (dot(axis, p) * (1 - c)); };
  const int mol = s.atoms.empty() ? 1 : std::max_element(s.atoms.begin(), s.atoms.end(), [](const Atom& a, const Atom& b) { return a.mol < b.mol; })->mol;
  int next_mol = mol + 1;
  constexpr double kPi = 3.14159265358979323846;
  for (const auto& st : chosen) {
    // the atoms a ligand here could touch
    std::vector<Vec3> near;
    for (const auto& a : s.atoms)
      if (norm(rel(st.p, a.pos)) < reach + 4) near.push_back(st.p + rel(st.p, a.pos));
    std::vector<Vec3> axes{st.n};
    const Vec3 t1 = perpendicular(st.n), t2 = cross(st.n, t1);
    for (double tilt : {15.0, 30.0})
      for (int az = 0; az < 6; ++az) {
        const double th = tilt * kPi / 180, ph = az * kPi / 3;
        axes.push_back(unitv(st.n * std::cos(th) + (t1 * std::cos(ph) + t2 * std::sin(ph)) * std::sin(th)));
      }
    double best = -1;
    std::vector<Vec3> placed;
    for (const Vec3& ax : axes) {
      const Vec3 axis0 = cross(v, ax);
      const double s0 = norm(axis0), c0 = dot(v, ax);
      auto align = [&](const Vec3& p) {
        if (s0 < 1e-9) return c0 > 0 ? p : rotate(p, perpendicular(v), -1, 0);
        return rotate(p, axis0 * (1.0 / s0), c0, s0);
      };
      for (int k = 0; k < 24; ++k) {
        const double ang = 2 * kPi * k / 24;
        std::vector<Vec3> trial(L.atoms.size());
        double dmin = 1e300;
        for (size_t i = 0; i < L.atoms.size(); ++i) {
          trial[i] = st.p + rotate(align(L.atoms[i].pos - L.atoms[X].pos), ax, std::cos(ang), std::sin(ang));
          if (i == D || i == X) continue;
          for (const auto& q : near) dmin = std::min(dmin, norm(trial[i] - q));
        }
        if (dmin > best) best = dmin, placed = trial;
      }
    }
    std::vector<int64_t> map(L.atoms.size(), -1);
    for (size_t i = 0; i < L.atoms.size(); ++i) {
      if (i == D) continue;
      Atom a = L.atoms[i];
      a.pos = placed[i];
      a.type = type_for(s, a.element);
      a.mol = next_mol;
      a.resname = "SR";
      a.resid = next_mol;
      a.id = s.atoms.empty() ? 1 : s.atoms.back().id + 1;
      a.charge = 0;
      a.name = std::string(element(a.element).symbol) + std::to_string(s.atoms.size() + 1);
      map[i] = int64_t(s.atoms.size());
      s.atoms.push_back(a);
      ++rep.added_atoms;
    }
    for (const auto& b : L.bonds)
      if (map[b.i] >= 0 && map[b.j] >= 0) s.bonds.push_back({uint32_t(map[b.i]), uint32_t(map[b.j]), b.order});
    ++next_mol;
    ++rep.ligands;
    (st.hollow ? rep.hollow : rep.on_top) += 1;
  }
  char b[320];
  std::snprintf(b, sizeof b, "%zu %s ligands on %zu surface metal atoms (%zu in three-fold hollows, %zu on top; S–S ≥ %.2f Å; %zu atoms added); "
                "no metal–S bonds written — relax with a force field that has metal–S terms before dynamics",
                rep.ligands, o.name.c_str(), rep.surface_atoms, rep.hollow, rep.on_top, spacing, rep.added_atoms);
  rep.notes.push_back(b);
  return rep;
}

std::vector<uint32_t> attach_fragment(System& s, uint32_t target, const std::string& smiles, int which, bool replace_h, const Vec3* direction) {
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
  else if (direction && norm(*direction) > 1e-9) dir = unitv(*direction);
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

std::pair<int, int> fix_hydrogens(System& s, const std::vector<char>& atoms) {
  // surplus first: a heavy atom with more bonds than its valence loses its extra hydrogens (the last bonded first)
  const auto sums = order_sums(s);
  const auto nb = neighbours(s);
  std::vector<char> drop(s.atoms.size(), 0);
  int removed = 0;
  for (uint32_t i = 0; i < s.atoms.size(); ++i) {
    if ((!atoms.empty() && (i >= atoms.size() || !atoms[i])) || s.atoms[i].element == 1) continue;
    const int charge = int(std::lround(s.atoms[i].charge));
    const int val = default_valence(s.atoms[i].element, std::fabs(s.atoms[i].charge - charge) < 0.05 ? charge : 0);
    if (val == 0) continue;
    int extra = int(std::lround(sums[i] - 1e-9)) - val;
    for (auto it = nb[i].rbegin(); it != nb[i].rend() && extra > 0; ++it)
      if (s.atoms[*it].element == 1 && !drop[*it] && nb[*it].size() == 1) { drop[*it] = 1; --extra; ++removed; }
  }
  std::vector<char> keep_marks = atoms;
  if (removed > 0) {
    // the marks follow the atoms that stay
    if (!atoms.empty()) {
      std::vector<char> m;
      for (size_t i = 0; i < s.atoms.size(); ++i) if (!drop[i]) m.push_back(i < atoms.size() ? atoms[i] : 0);
      keep_marks = m;
    }
    delete_atoms(s, drop);
  }
  const int added = add_hydrogens(s, keep_marks);
  return {added, removed};
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


int point_defects(System& s, const DefectOptions& o, std::vector<std::string>* notes) {
  if (o.from <= 0) throw std::invalid_argument("give the element the defects replace");
  if (o.to == o.from) throw std::invalid_argument("a substitution needs another element");
  std::vector<uint32_t> cand;
  for (uint32_t i = 0; i < s.atoms.size(); ++i)
    if (s.atoms[i].element == o.from && (o.region.empty() || (i < o.region.size() && o.region[i]))) cand.push_back(i);
  if (cand.empty()) throw std::invalid_argument(std::string("no ") + element(o.from).symbol + " atoms to pick");
  const size_t want = o.count > 0 ? size_t(o.count) : size_t(std::lround(o.fraction * double(cand.size())));
  std::mt19937_64 rng(o.seed);
  std::shuffle(cand.begin(), cand.end(), rng);
  std::vector<uint32_t> pick;
  for (uint32_t i : cand) {
    if (pick.size() >= want) break;
    bool ok = true;
    for (uint32_t p : pick) {
      const Vec3 d = s.cell.valid() ? s.cell.minimum_image(s.atoms[i].pos - s.atoms[p].pos) : s.atoms[i].pos - s.atoms[p].pos;
      if (norm(d) < o.min_spacing) { ok = false; break; }
    }
    if (ok) pick.push_back(i);
  }
  double q = 0;
  for (uint32_t i : pick) q += s.atoms[i].charge;
  const std::string sym = element(o.from).symbol;
  if (o.to == 0) {
    std::vector<char> rm(s.atoms.size(), 0);
    for (uint32_t i : pick) rm[i] = 1;
    delete_atoms(s, rm);
    if (notes) {
      notes->push_back(std::to_string(pick.size()) + " " + sym + " vacancies (of " + std::to_string(cand.size()) + " " + sym + ")");
      if (std::fabs(q) > 1e-6) notes->push_back("the removed atoms carried " + std::to_string(q).substr(0, 7) + " e: the structure is no longer neutral (compensate, or reassign the charges)");
    }
  } else {
    for (uint32_t i : pick) {
      s.atoms[i].element = o.to;
      s.atoms[i].name = element(o.to).symbol;
    }
    if (notes) notes->push_back(std::to_string(pick.size()) + " " + sym + " → " + element(o.to).symbol + " (of " + std::to_string(cand.size()) + " " + sym +
                                "); types and charges to be assigned again");
  }
  if (pick.size() < want && notes) notes->push_back("only " + std::to_string(pick.size()) + " of " + std::to_string(want) + " sites at " + std::to_string(o.min_spacing).substr(0, 4) + " Å spacing");
  return int(pick.size());
}

}  // namespace caps
