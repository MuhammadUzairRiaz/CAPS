// CAPS Exchange, topology (see topology_edit.hpp).
#include "caps/topology_edit.hpp"

#include <algorithm>

namespace caps {

namespace {

using Pair = TopologyEditor::Pair;

Pair pair_of(uint32_t i, uint32_t j) { return i < j ? Pair{i, j} : Pair{j, i}; }

std::array<uint32_t, 3> triple_of(uint32_t i, uint32_t j, uint32_t k) { return i < k ? std::array<uint32_t, 3>{i, j, k} : std::array<uint32_t, 3>{k, j, i}; }

std::array<uint32_t, 4> quad_of(uint32_t i, uint32_t j, uint32_t k, uint32_t l) {
  const std::array<uint32_t, 4> f{i, j, k, l}, r{l, k, j, i};
  return f < r ? f : r;
}

}  // namespace

TopologyEditor::Key2 TopologyEditor::key2(uint32_t i, uint32_t j) const {
  const int a = type(i), b = type(j);
  return a < b ? Key2{a, b} : Key2{b, a};
}

TopologyEditor::Key3 TopologyEditor::key3(uint32_t i, uint32_t j, uint32_t k) const {
  const int a = type(i), b = type(j), c = type(k);
  return a <= c ? Key3{a, b, c} : Key3{c, b, a};
}

TopologyEditor::Key4 TopologyEditor::key4(uint32_t i, uint32_t j, uint32_t k, uint32_t l) const {
  const Key4 f{type(i), type(j), type(k), type(l)}, r{type(l), type(k), type(j), type(i)};
  return f < r ? f : r;
}

TopologyEditor::ImproperKey TopologyEditor::improper_key(uint32_t c, const std::vector<std::vector<uint32_t>>& nb) const {
  std::vector<int> t;
  for (uint32_t a : nb[c]) t.push_back(type(a));
  std::sort(t.begin(), t.end());
  return {type(c), t};
}

bool TopologyEditor::bonded(uint32_t i, uint32_t j) const {
  return i < n_ && std::find(nb_[i].begin(), nb_[i].end(), j) != nb_[i].end();
}

TopologyEditor::BondParams TopologyEditor::bond_params(uint32_t i, uint32_t j) const {
  BondParams p;
  const Pair q = pair_of(i, j);
  for (const auto& b : ff_.bonds)
    if (pair_of(b.i, b.j) == q) { p.kind = 0, p.b = b; return p; }
  for (const auto& b : ff_.bonds_x)
    if (pair_of(b.i, b.j) == q) { p.kind = 1, p.x = b; return p; }
  for (const auto& b : ff_.bonds2)
    if (pair_of(b.i, b.j) == q) { p.kind = 2, p.c = b; return p; }
  return p;
}

void TopologyEditor::graph_pairs(const std::vector<std::vector<uint32_t>>& nb, std::set<Pair>& p12, std::set<Pair>& p13,
                                 std::set<Pair>& p14) const {
  p12.clear(), p13.clear(), p14.clear();
  for (uint32_t i = 0; i < nb.size(); ++i)
    for (uint32_t j : nb[i])
      if (i < j) p12.insert({i, j});
  for (uint32_t j = 0; j < nb.size(); ++j)
    for (size_t a = 0; a < nb[j].size(); ++a)
      for (size_t b = a + 1; b < nb[j].size(); ++b) {
        const Pair q = pair_of(nb[j][a], nb[j][b]);
        if (!p12.count(q)) p13.insert(q);
      }
  for (const auto& [j, k] : p12)
    for (uint32_t i : nb[j])
      for (uint32_t l : nb[k]) {
        if (i == k || l == j || i == l) continue;
        const Pair q = pair_of(i, l);
        if (!p12.count(q) && !p13.count(q)) p14.insert(q);
      }
}

TopologyEditor::TopologyEditor(System& s, ForceField& ff) : s_(s), ff_(ff), n_(s.atoms.size()) {
  if (ff.type_index.size() != n_) throw TopologyError("the force field is not for this structure's atoms");
  nb_.assign(n_, {});
  for (const auto& b : s.bonds) {
    if (b.i >= n_ || b.j >= n_ || b.i == b.j) continue;
    if (!bonded(b.i, b.j)) nb_[b.i].push_back(b.j), nb_[b.j].push_back(b.i);
  }
  for (auto& v : nb_) std::sort(v.begin(), v.end());
  auto check_bond = [&](uint32_t i, uint32_t j) {
    if (!bonded(i, j))
      throw TopologyError("the force field has a bond term " + std::to_string(i + 1) + "-" + std::to_string(j + 1) + " the structure does not have as a bond");
  };
  for (const auto& b : ff.bonds) check_bond(b.i, b.j);
  for (const auto& b : ff.bonds_x) check_bond(b.i, b.j);
  for (const auto& b : ff.bonds2) check_bond(b.i, b.j);

  // the exclusion rule: 1-3 and 1-4 partners excluded when most of them are in the force field's own list
  std::set<Pair> g12, g13, g14, ex;
  graph_pairs(nb_, g12, g13, g14);
  for (uint32_t i = 0; i < ff.excluded.size(); ++i)
    for (uint32_t j : ff.excluded[i])
      if (i < j) ex.insert({i, j});
  auto majority = [&](const std::set<Pair>& g, bool def) {
    if (g.empty()) return def;
    size_t in = 0;
    for (const auto& p : g) in += ex.count(p);
    return 2 * in > g.size();
  };
  x13_ = majority(g13, !ff.keep13);
  x14_ = majority(g14, !ff.pairs14.empty());
  for (const auto& p : ex)
    if (!g12.count(p) && !(x13_ && g13.count(p)) && !(x14_ && g14.count(p))) extra_excluded_.insert(p);
  // 1-4 pairs: rebuilt from the graph when they are its 1-4 pairs
  if (!ff.pairs14.empty()) {
    std::set<Pair> p;
    for (const auto& q : ff.pairs14) p.insert(pair_of(q[0], q[1]));
    p14_from_graph_ = p == g14;
    if (ff.per_pair14())
      for (size_t k = 0; k < ff.pairs14.size(); ++k) {
        const std::array<double, 2> sc{ff.pairs14_lj[k], k < ff.pairs14_coul.size() ? ff.pairs14_coul[k] : ff.coul14};
        p14_scale_[pair_of(ff.pairs14[k][0], ff.pairs14[k][1])] = sc;
        p14_type_scale_.emplace(key2(ff.pairs14[k][0], ff.pairs14[k][1]), sc);
      }
  }

  // the force field's terms as tables by type: each combination from the first set of atoms that has it
  for (const auto& p : g12) bond_seen_.insert(key2(p[0], p[1]));
  for (const auto& b : ff.bonds) { BondParams q; q.kind = 0, q.b = b; bond_tab_.emplace(key2(b.i, b.j), q); }
  for (const auto& b : ff.bonds_x) { BondParams q; q.kind = 1, q.x = b; bond_tab_.emplace(key2(b.i, b.j), q); }
  for (const auto& b : ff.bonds2) { BondParams q; q.kind = 2, q.c = b; bond_tab_.emplace(key2(b.i, b.j), q); }

  for (uint32_t j = 0; j < n_; ++j)
    for (size_t a = 0; a < nb_[j].size(); ++a)
      for (size_t b = a + 1; b < nb_[j].size(); ++b) angle_seen_.insert(key3(nb_[j][a], j, nb_[j][b]));
  std::map<std::array<uint32_t, 3>, AngleSet> by_triple;
  for (const auto& a : ff.angles) by_triple[triple_of(a.i, a.j, a.k)].h.push_back(a);
  for (const auto& a : ff.angles_x) by_triple[triple_of(a.i, a.j, a.k)].x.push_back(a);
  {
    std::map<Pair, std::vector<UreyBradley>> ub;
    for (const auto& u : ff.urey_bradley) ub[pair_of(u.i, u.k)].push_back(u);
    for (auto& [t, set] : by_triple) {
      const auto it = ub.find(pair_of(t[0], t[2]));
      if (it != ub.end()) set.ub = it->second;
    }
  }
  for (const auto& [t, set] : by_triple) angle_tab_.emplace(key3(t[0], t[1], t[2]), set);

  for (const auto& [j, k] : g12)
    for (uint32_t i : nb_[j])
      for (uint32_t l : nb_[k])
        if (i != k && l != j && i != l) dihedral_seen_.insert(key4(i, j, k, l));
  std::map<std::array<uint32_t, 4>, std::vector<TorsionTerm>> by_quad;
  for (const auto& d : ff.dihedrals) by_quad[quad_of(d.i, d.j, d.k, d.l)].push_back(d);
  for (const auto& [q, terms] : by_quad) dihedral_tab_.emplace(key4(q[0], q[1], q[2], q[3]), terms);

  // impropers: each term's centre is the atom bonded to the other three
  auto centre_of = [&](std::array<uint32_t, 4> a) -> int64_t {
    for (int c = 0; c < 4; ++c) {
      int k = 0;
      for (int o = 0; o < 4; ++o)
        if (o != c && bonded(a[c], a[o])) ++k;
      if (k == 3) return a[c];
    }
    return -1;
  };
  std::map<uint32_t, std::vector<ImproperPattern>> by_centre;
  auto pattern = [&](uint32_t c, std::array<uint32_t, 4> a) {
    std::array<int, 4> sl;
    for (int k = 0; k < 4; ++k) sl[k] = a[k] == c ? -1 : type(a[k]);
    return sl;
  };
  for (const auto& t : ff.impropers) {
    const int64_t c = centre_of({t.i, t.j, t.k, t.l});
    if (c >= 0) { ImproperPattern p{0, pattern(uint32_t(c), {t.i, t.j, t.k, t.l})}; p.t = t; by_centre[uint32_t(c)].push_back(p); }
  }
  for (const auto& t : ff.impropers_harmonic) {
    const int64_t c = centre_of({t.i, t.j, t.k, t.l});
    if (c >= 0) { ImproperPattern p{1, pattern(uint32_t(c), {t.i, t.j, t.k, t.l})}; p.h = t; by_centre[uint32_t(c)].push_back(p); }
  }
  for (const auto& t : ff.inversions) {
    ImproperPattern p{2, pattern(t.c, {t.c, t.a, t.b, t.d})};
    p.v = t;
    by_centre[t.c].push_back(p);
  }
  for (uint32_t c = 0; c < n_; ++c)
    if (nb_[c].size() == 3) improper_seen_.insert(improper_key(c, nb_));
  for (const auto& [c, pats] : by_centre)
    if (nb_[c].size() == 3) improper_tab_.emplace(improper_key(c, nb_), pats);
}

void TopologyEditor::change(const std::vector<Pair>& broken, const std::vector<Pair>& made) { change_impl(broken, made, {}); }

size_t TopologyEditor::learn(const System& frag, const ForceField& f, const std::vector<char>& complete) {
  const size_t m = frag.atoms.size();
  if (f.type_index.size() != m || complete.size() != m) throw TopologyError("learn: the force field and the marks must be for the fragment's atoms");
  std::map<std::string, int> ours;
  for (size_t k = 0; k < ff_.type_names.size(); ++k) ours.emplace(ff_.type_names[k], int(k));
  std::vector<int> ty(m, -1);   // each fragment atom's type here (−1: a type this force field does not have, or incomplete)
  for (size_t a = 0; a < m; ++a) {
    const int fi = f.type_index[a];
    if (!complete[a] || fi < 0 || size_t(fi) >= f.type_names.size()) continue;
    if (const auto it = ours.find(f.type_names[fi]); it != ours.end()) ty[a] = it->second;
  }
  std::vector<std::vector<uint32_t>> nb(m);
  for (const auto& b : frag.bonds)
    if (b.i < m && b.j < m && b.i != b.j) nb[b.i].push_back(b.j), nb[b.j].push_back(b.i);
  for (auto& v : nb) std::sort(v.begin(), v.end()), v.erase(std::unique(v.begin(), v.end()), v.end());
  auto ok = [&](std::initializer_list<uint32_t> a) {
    for (uint32_t x : a)
      if (x >= m || ty[x] < 0) return false;
    return true;
  };
  auto k2 = [&](uint32_t i, uint32_t j) { return ty[i] < ty[j] ? Key2{ty[i], ty[j]} : Key2{ty[j], ty[i]}; };
  auto k3 = [&](uint32_t i, uint32_t j, uint32_t k) { return ty[i] <= ty[k] ? Key3{ty[i], ty[j], ty[k]} : Key3{ty[k], ty[j], ty[i]}; };
  auto k4 = [&](uint32_t i, uint32_t j, uint32_t k, uint32_t l) {
    const Key4 a{ty[i], ty[j], ty[k], ty[l]}, b{ty[l], ty[k], ty[j], ty[i]};
    return a < b ? a : b;
  };
  size_t added = 0;
  auto bond = [&](uint32_t i, uint32_t j, const BondParams& q) {
    if (!ok({i, j})) return;
    bond_seen_.insert(k2(i, j));
    added += bond_tab_.emplace(k2(i, j), q).second;
  };
  for (const auto& b : f.bonds) { BondParams q; q.kind = 0, q.b = b; bond(b.i, b.j, q); }
  for (const auto& b : f.bonds_x) { BondParams q; q.kind = 1, q.x = b; bond(b.i, b.j, q); }
  for (const auto& b : f.bonds2) { BondParams q; q.kind = 2, q.c = b; bond(b.i, b.j, q); }
  std::map<std::array<uint32_t, 3>, AngleSet> by_triple;
  for (const auto& a : f.angles) by_triple[triple_of(a.i, a.j, a.k)].h.push_back(a);
  for (const auto& a : f.angles_x) by_triple[triple_of(a.i, a.j, a.k)].x.push_back(a);
  for (const auto& u : f.urey_bradley)
    for (auto& [tr, set] : by_triple)
      if (pair_of(tr[0], tr[2]) == pair_of(u.i, u.k)) set.ub.push_back(u);
  for (const auto& [tr, set] : by_triple)
    if (ok({tr[0], tr[1], tr[2]})) added += angle_tab_.emplace(k3(tr[0], tr[1], tr[2]), set).second;
  for (uint32_t j = 0; j < m; ++j)
    for (size_t a = 0; a < nb[j].size(); ++a)
      for (size_t b = a + 1; b < nb[j].size(); ++b)
        if (ok({nb[j][a], j, nb[j][b]})) angle_seen_.insert(k3(nb[j][a], j, nb[j][b]));
  std::map<std::array<uint32_t, 4>, std::vector<TorsionTerm>> by_quad;
  for (const auto& d : f.dihedrals) by_quad[quad_of(d.i, d.j, d.k, d.l)].push_back(d);
  for (const auto& [q, terms] : by_quad)
    if (ok({q[0], q[1], q[2], q[3]})) added += dihedral_tab_.emplace(k4(q[0], q[1], q[2], q[3]), terms).second;
  for (uint32_t j = 0; j < m; ++j)
    for (uint32_t k : nb[j]) {
      if (k < j) continue;
      for (uint32_t i : nb[j])
        for (uint32_t l : nb[k])
          if (i != k && l != j && i != l && ok({i, j, k, l})) dihedral_seen_.insert(k4(i, j, k, l));
    }
  return added;
}

void TopologyEditor::change(const std::vector<Pair>& broken, const std::vector<Pair>& made, const std::vector<Retype>& retype) {
  if (retype.empty()) return change_impl(broken, made, {});
  std::set<Pair> b, m, fresh;
  for (const auto& p : broken) b.insert(pair_of(p[0], p[1]));
  for (const auto& p : made) m.insert(pair_of(p[0], p[1]));
  for (const auto& r : retype) {
    if (r.atom >= n_) throw TopologyError("retype: an atom the structure does not have");
    if (r.type < 0 || size_t(r.type) >= ff_.type_names.size()) throw TopologyError("retype: a type the force field does not have");
    for (uint32_t q : nb_[r.atom])
      if (const Pair p = pair_of(r.atom, q); !b.count(p)) fresh.insert(p);
  }
  for (const auto& p : fresh) b.insert(p), m.insert(p);
  struct Old { int type; std::string name; double charge; };
  std::vector<Old> old;
  auto set = [&](uint32_t a, int type, const std::string& name, double charge) {
    ff_.type_index[a] = type;
    if (a < ff_.atom_type.size()) ff_.atom_type[a] = name;
    if (a < ff_.charge.size()) ff_.charge[a] = charge;
  };
  for (const auto& r : retype)
    old.push_back({ff_.type_index[r.atom], r.atom < ff_.atom_type.size() ? ff_.atom_type[r.atom] : std::string(), r.atom < ff_.charge.size() ? ff_.charge[r.atom] : 0.0});
  for (const auto& r : retype) set(r.atom, r.type, r.name, r.charge);
  try {
    change_impl(std::vector<Pair>(b.begin(), b.end()), std::vector<Pair>(m.begin(), m.end()), {}, fresh);
  } catch (...) {
    for (size_t k = retype.size(); k-- > 0;) set(retype[k].atom, old[k].type, old[k].name, old[k].charge);
    throw;
  }
}

void TopologyEditor::define_bond(const BondTerm& like) {
  if (like.i >= n_ || like.j >= n_) throw TopologyError("define_bond: its atoms must be atoms of the structure");
  BondParams q;
  q.kind = 0, q.b = like;
  bond_tab_[key2(like.i, like.j)] = q;
}

void TopologyEditor::define_bond(const BondX& like) {
  if (like.i >= n_ || like.j >= n_) throw TopologyError("define_bond: its atoms must be atoms of the structure");
  BondParams q;
  q.kind = 1, q.x = like;
  bond_tab_[key2(like.i, like.j)] = q;
}

void TopologyEditor::swap(uint32_t a1, uint32_t a2, uint32_t b1, uint32_t b2) {
  if (!bonded(a1, a2) || !bonded(b1, b2)) throw TopologyError("a swap needs the bonds a1-a2 and b1-b2");
  std::map<Pair, BondParams> hint;
  auto moved = [&](uint32_t keep, uint32_t from, uint32_t to) {
    BondParams p = bond_params(keep, from);
    auto relabel = [&](uint32_t& i, uint32_t& j) { (i == from ? i : j) = to; };
    if (p.kind == 0) relabel(p.b.i, p.b.j);
    if (p.kind == 1) relabel(p.x.i, p.x.j);
    if (p.kind == 2) relabel(p.c.i, p.c.j);
    hint[pair_of(keep, to)] = p;
  };
  moved(a1, a2, b2);
  moved(b1, b2, a2);
  change_impl({pair_of(a1, a2), pair_of(b1, b2)}, {pair_of(a1, b2), pair_of(b1, a2)}, hint);
}

void TopologyEditor::change_impl(const std::vector<Pair>& broken0, const std::vector<Pair>& made0, const std::map<Pair, BondParams>& hint,
                                 const std::set<Pair>& fresh) {
  // ---- checks: nothing is changed until every new term is known
  if (!ff_.angles2.empty() || !ff_.dihedrals2.empty() || !ff_.impropers2.empty() || !ff_.cbt.empty())
    throw TopologyError("class II cross terms and combined bending-torsion terms are not rebuilt by a bond change: assign the force field again to the new structure");
  if (!ff_.pairs14.empty() && !p14_from_graph_)
    throw TopologyError("the force field's 1-4 pairs are not the bond graph's: assign the force field again to the new structure");
  std::set<Pair> broken, made;
  std::set<uint32_t> touched;
  for (const auto& p0 : broken0) {
    const Pair p = pair_of(p0[0], p0[1]);
    if (p[1] >= n_ || !bonded(p[0], p[1])) throw TopologyError("no bond " + std::to_string(p[0] + 1) + "-" + std::to_string(p[1] + 1) + " to break");
    broken.insert(p), touched.insert(p[0]), touched.insert(p[1]);
  }
  for (const auto& p0 : made0) {
    const Pair p = pair_of(p0[0], p0[1]);
    if (p[0] == p[1] || p[1] >= n_) throw TopologyError("a bond needs two different atoms of the structure");
    if (bonded(p[0], p[1]) && !broken.count(p)) throw TopologyError("the bond " + std::to_string(p[0] + 1) + "-" + std::to_string(p[1] + 1) + " is there already");
    made.insert(p), touched.insert(p[0]), touched.insert(p[1]);
  }
  for (const auto& vs : ff_.vsites) {
    bool hit = touched.count(vs.site) > 0;
    for (uint32_t a : vs.from) hit = hit || touched.count(a) > 0;
    if (hit) throw TopologyError("a virtual site is built on an atom whose bonds change");
  }
  // the new graph
  auto nb = nb_;
  for (const auto& p : broken) {
    nb[p[0]].erase(std::find(nb[p[0]].begin(), nb[p[0]].end(), p[1]));
    nb[p[1]].erase(std::find(nb[p[1]].begin(), nb[p[1]].end(), p[0]));
  }
  for (const auto& p : made) nb[p[0]].push_back(p[1]), nb[p[1]].push_back(p[0]);
  for (uint32_t a : touched) std::sort(nb[a].begin(), nb[a].end());
  auto through = [&](uint32_t i, uint32_t j) { return broken.count(pair_of(i, j)) > 0; };
  auto name = [&](int t) { return t >= 0 && size_t(t) < ff_.type_names.size() ? ff_.type_names[t] : std::to_string(t); };

  // bonds
  auto bonds = ff_.bonds;
  auto bonds_x = ff_.bonds_x;
  auto bonds2 = ff_.bonds2;
  auto stash = stash_;
  for (const auto& p : broken) stash[p] = bond_params(p[0], p[1]);
  bonds.erase(std::remove_if(bonds.begin(), bonds.end(), [&](const BondTerm& b) { return through(b.i, b.j); }), bonds.end());
  bonds_x.erase(std::remove_if(bonds_x.begin(), bonds_x.end(), [&](const BondX& b) { return through(b.i, b.j); }), bonds_x.end());
  bonds2.erase(std::remove_if(bonds2.begin(), bonds2.end(), [&](const Class2Bond& b) { return through(b.i, b.j); }), bonds2.end());
  for (const auto& p : made) {
    BondParams q;
    bool known = false;
    if (const auto h = hint.find(p); h != hint.end()) q = h->second, known = true;
    else if (const auto st = stash.find(p); !fresh.count(p) && st != stash.end()) q = st->second, known = true;
    else if (const auto t = bond_tab_.find(key2(p[0], p[1])); t != bond_tab_.end()) q = t->second, known = true;
    else if (bond_seen_.count(key2(p[0], p[1]))) q.kind = -1, known = true;   // these types are bonded without a term
    if (!known)
      throw TopologyError("no bond parameters between types " + name(type(p[0])) + " and " + name(type(p[1])) + ": assign the force field again to the new structure");
    if (q.kind == 0) q.b.i = p[0], q.b.j = p[1], bonds.push_back(q.b);
    if (q.kind == 1) q.x.i = p[0], q.x.j = p[1], bonds_x.push_back(q.x);
    if (q.kind == 2) q.c.i = p[0], q.c.j = p[1], bonds2.push_back(q.c);
  }

  // angles through a broken bond go; every angle through a made bond comes from the table
  auto angles = ff_.angles;
  auto angles_x = ff_.angles_x;
  auto ub = ff_.urey_bradley;
  std::set<Pair> gone_ends;
  auto angle_gone = [&](uint32_t i, uint32_t j, uint32_t k) {
    if (!through(i, j) && !through(j, k)) return false;
    gone_ends.insert(pair_of(i, k));
    return true;
  };
  angles.erase(std::remove_if(angles.begin(), angles.end(), [&](const AngleTerm& a) { return angle_gone(a.i, a.j, a.k); }), angles.end());
  angles_x.erase(std::remove_if(angles_x.begin(), angles_x.end(), [&](const AngleX& a) { return angle_gone(a.i, a.j, a.k); }), angles_x.end());
  ub.erase(std::remove_if(ub.begin(), ub.end(), [&](const UreyBradley& u) { return gone_ends.count(pair_of(u.i, u.k)) > 0; }), ub.end());
  std::set<std::array<uint32_t, 3>> new_angles;
  for (const auto& p : made)
    for (int o = 0; o < 2; ++o) {
      const uint32_t a = p[o], b = p[1 - o];
      for (uint32_t x : nb[a])
        if (x != b) new_angles.insert(triple_of(x, a, b));
    }
  for (const auto& t : new_angles) {
    const auto it = angle_tab_.find(key3(t[0], t[1], t[2]));
    if (it == angle_tab_.end()) {
      if (angle_seen_.count(key3(t[0], t[1], t[2])) || (ff_.angles.empty() && ff_.angles_x.empty())) continue;   // no angle term for these
      throw TopologyError("no angle parameters for types " + name(type(t[0])) + "-" + name(type(t[1])) + "-" + name(type(t[2])) +
                          ": assign the force field again to the new structure");
    }
    for (AngleTerm a : it->second.h) a.i = t[0], a.j = t[1], a.k = t[2], angles.push_back(a);
    for (AngleX a : it->second.x) a.i = t[0], a.j = t[1], a.k = t[2], angles_x.push_back(a);
    for (UreyBradley u : it->second.ub) u.i = t[0], u.k = t[2], ub.push_back(u);
  }

  // proper torsions: likewise
  auto dihedrals = ff_.dihedrals;
  dihedrals.erase(std::remove_if(dihedrals.begin(), dihedrals.end(),
                                 [&](const TorsionTerm& d) { return through(d.i, d.j) || through(d.j, d.k) || through(d.k, d.l); }),
                  dihedrals.end());
  std::set<std::array<uint32_t, 4>> new_quads;
  for (const auto& p : made)
    for (int o = 0; o < 2; ++o) {
      const uint32_t a = p[o], b = p[1 - o];
      for (uint32_t i : nb[a])   // the made bond in the middle
        for (uint32_t l : nb[b])
          if (i != b && l != a && i != l) new_quads.insert(quad_of(i, a, b, l));
      for (uint32_t k : nb[b])   // the made bond at an end
        for (uint32_t l : nb[k])
          if (k != a && l != b && l != a) new_quads.insert(quad_of(a, b, k, l));
    }
  for (const auto& q : new_quads) {
    const auto it = dihedral_tab_.find(key4(q[0], q[1], q[2], q[3]));
    if (it == dihedral_tab_.end()) {
      if (dihedral_seen_.count(key4(q[0], q[1], q[2], q[3])) || ff_.dihedrals.empty()) continue;   // no torsion term for these
      throw TopologyError("no torsion parameters for types " + name(type(q[0])) + "-" + name(type(q[1])) + "-" + name(type(q[2])) + "-" +
                          name(type(q[3])) + ": assign the force field again to the new structure");
    }
    for (TorsionTerm d : it->second) d.i = q[0], d.j = q[1], d.k = q[2], d.l = q[3], dihedrals.push_back(d);
  }

  // impropers at the atoms whose neighbours changed: dropped and placed again from the table
  auto impropers = ff_.impropers;
  auto impropers_h = ff_.impropers_harmonic;
  auto inversions = ff_.inversions;
  auto at_touched = [&](std::array<uint32_t, 4> a) {
    for (uint32_t x : a)
      if (touched.count(x)) {
        int k = 0;   // x is the centre when it is bonded (before the change) to the other three
        for (uint32_t y : a)
          if (y != x && bonded(x, y)) ++k;
        if (k == 3) return true;
      }
    return false;
  };
  impropers.erase(std::remove_if(impropers.begin(), impropers.end(), [&](const TorsionTerm& t) { return at_touched({t.i, t.j, t.k, t.l}); }), impropers.end());
  impropers_h.erase(std::remove_if(impropers_h.begin(), impropers_h.end(), [&](const HarmonicTorsion& t) { return at_touched({t.i, t.j, t.k, t.l}); }),
                    impropers_h.end());
  inversions.erase(std::remove_if(inversions.begin(), inversions.end(), [&](const InversionTerm& t) { return touched.count(t.c) > 0; }), inversions.end());
  size_t missing = 0;
  const bool any_impropers = !ff_.impropers.empty() || !ff_.impropers_harmonic.empty() || !ff_.inversions.empty();
  for (uint32_t c : touched) {
    if (nb[c].size() != 3) continue;
    const auto key = improper_key(c, nb);
    const auto it = improper_tab_.find(key);
    if (it == improper_tab_.end()) {
      if (any_impropers && !improper_seen_.count(key)) ++missing;
      continue;
    }
    for (const auto& pat : it->second) {
      std::array<uint32_t, 4> a{};
      std::vector<char> used(3, 0);
      bool ok = true;
      for (int k = 0; k < 4 && ok; ++k) {
        if (pat.slot[k] < 0) { a[k] = c; continue; }
        ok = false;
        for (size_t m = 0; m < 3; ++m)
          if (!used[m] && type(nb[c][m]) == pat.slot[k]) { used[m] = 1, a[k] = nb[c][m], ok = true; break; }
      }
      if (!ok) continue;
      if (pat.kind == 0) { TorsionTerm t = pat.t; t.i = a[0], t.j = a[1], t.k = a[2], t.l = a[3]; impropers.push_back(t); }
      if (pat.kind == 1) { HarmonicTorsion t = pat.h; t.i = a[0], t.j = a[1], t.k = a[2], t.l = a[3]; impropers_h.push_back(t); }
      if (pat.kind == 2) { InversionTerm t = pat.v; t.c = a[0], t.a = a[1], t.b = a[2], t.d = a[3]; inversions.push_back(t); }
    }
  }

  // ---- commit
  ff_.bonds = std::move(bonds), ff_.bonds_x = std::move(bonds_x), ff_.bonds2 = std::move(bonds2);
  ff_.angles = std::move(angles), ff_.angles_x = std::move(angles_x), ff_.urey_bradley = std::move(ub);
  ff_.dihedrals = std::move(dihedrals);
  ff_.impropers = std::move(impropers), ff_.impropers_harmonic = std::move(impropers_h), ff_.inversions = std::move(inversions);
  if (!ff_.impropers_dlpoly.empty()) ff_.impropers_dlpoly_ok = false;   // the DL_POLY order is no longer known for the changed atoms
  stash_ = std::move(stash);
  missing_impropers_ += missing;
  nb_ = std::move(nb);
  s_.bonds.erase(std::remove_if(s_.bonds.begin(), s_.bonds.end(), [&](const Bond& b) { return broken.count(pair_of(b.i, b.j)) > 0; }), s_.bonds.end());
  for (const auto& p : made) s_.bonds.push_back({p[0], p[1], 1});

  std::set<Pair> g12, g13, g14;
  graph_pairs(nb_, g12, g13, g14);
  if (p14_from_graph_) {
    const bool per = ff_.per_pair14();
    ff_.pairs14.clear(), ff_.pairs14_lj.clear(), ff_.pairs14_coul.clear();
    for (const auto& p : g14) {
      ff_.pairs14.push_back({p[0], p[1]});
      if (!per) continue;
      std::array<double, 2> sc{ff_.lj14, ff_.coul14};
      if (const auto it = p14_scale_.find(p); it != p14_scale_.end()) sc = it->second;
      else if (const auto t = p14_type_scale_.find(key2(p[0], p[1])); t != p14_type_scale_.end()) sc = t->second;
      ff_.pairs14_lj.push_back(sc[0]), ff_.pairs14_coul.push_back(sc[1]);
    }
  }
  ff_.excluded.assign(n_, {});
  auto add = [&](const Pair& p) { ff_.excluded[p[0]].push_back(p[1]), ff_.excluded[p[1]].push_back(p[0]); };
  for (const auto& p : g12) add(p);
  if (x13_)
    for (const auto& p : g13) add(p);
  if (x14_)
    for (const auto& p : g14) add(p);
  for (const auto& p : extra_excluded_) add(p);
  for (auto& e : ff_.excluded) {
    std::sort(e.begin(), e.end());
    e.erase(std::unique(e.begin(), e.end()), e.end());
  }
  ++changes_;
}

}  // namespace caps
