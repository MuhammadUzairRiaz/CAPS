// A reaction scheme written as atom-mapped SMILES turned into a CAPS reaction template: see caps/react.hpp.
#include <algorithm>
#include <functional>
#include <map>
#include <set>
#include <sstream>

#include "caps/bond_react.hpp"
#include "caps/molecule.hpp"
#include "caps/react.hpp"

namespace caps {

namespace {

// The molecules of one side as one graph, hydrogens explicit.
struct Side {
  std::vector<int> element, map, comp, hcount;
  std::vector<char> aromatic, ring3;
  std::vector<std::set<int>> nb;
  std::vector<std::vector<int>> hs;   // each heavy atom's hydrogens
  int parts = 0;
  bool heavy(int i) const { return element[size_t(i)] != 1; }
  int degree(int i) const { return int(nb[size_t(i)].size()); }
};

Side side_of(const std::vector<std::string>& smiles) {
  Side s;
  for (const auto& sm : smiles) {
    MolGraph g = parse_smiles(sm);
    add_hydrogens(g);
    const int off = int(s.element.size());
    for (const auto& a : g.atoms) {
      s.element.push_back(a.element);
      s.map.push_back(a.map);
      s.aromatic.push_back(a.aromatic ? 1 : 0);
      s.comp.push_back(0);
    }
    s.nb.resize(s.element.size());
    for (const auto& b : g.bonds) {
      s.nb[size_t(off + b.a)].insert(off + b.b);
      s.nb[size_t(off + b.b)].insert(off + b.a);
    }
  }
  const size_t n = s.element.size();
  // components (a SMILES entry may hold several molecules)
  std::vector<int> comp(n, -1);
  for (size_t i = 0; i < n; ++i) {
    if (comp[i] >= 0) continue;
    std::vector<int> st{int(i)};
    comp[i] = s.parts;
    while (!st.empty()) {
      const int u = st.back();
      st.pop_back();
      for (int v : s.nb[size_t(u)]) if (comp[size_t(v)] < 0) { comp[size_t(v)] = s.parts; st.push_back(v); }
    }
    ++s.parts;
  }
  s.comp = comp;
  s.hcount.assign(n, 0);
  s.hs.assign(n, {});
  s.ring3.assign(n, 0);
  for (size_t i = 0; i < n; ++i)
    for (int v : s.nb[i])
      if (s.element[size_t(v)] == 1) { ++s.hcount[i]; s.hs[i].push_back(v); }
  for (size_t i = 0; i < n; ++i)
    for (int u : s.nb[i])
      for (int v : s.nb[i])
        if (u < v && s.nb[size_t(u)].count(v)) s.ring3[i] = 1;
  return s;
}

}  // namespace

ReactionTemplate template_from_scheme(const std::vector<std::string>& reactant_smiles, const std::vector<std::string>& product_smiles,
                                      const std::string& name, std::vector<std::string>* notes) {
  if (reactant_smiles.empty() || product_smiles.empty()) throw ReactError("a scheme needs reactants and products");
  const Side R = side_of(reactant_smiles), P = side_of(product_smiles);
  const int nr = int(R.element.size()), np = int(P.element.size());
  std::vector<int> to(size_t(nr), -1), from(size_t(np), -1);
  auto pair = [&](int r, int p) { to[size_t(r)] = p; from[size_t(p)] = r; };
  // 1. map numbers
  std::map<int, int> pmap;
  for (int p = 0; p < np; ++p) if (P.map[size_t(p)] > 0) pmap[P.map[size_t(p)]] = p;
  std::set<int> used_maps;
  for (int r = 0; r < nr; ++r) {
    const int m = R.map[size_t(r)];
    if (m <= 0) continue;
    if (!used_maps.insert(m).second) throw ReactError("map number " + std::to_string(m) + " is used twice among the reactants");
    auto it = pmap.find(m);
    if (it == pmap.end()) continue;
    if (P.element[size_t(it->second)] != R.element[size_t(r)]) throw ReactError("map " + std::to_string(m) + " changes element");
    pair(r, it->second);
  }
  // 2. the unmapped heavy atoms: grown from paired neighbours (same element, number of bonds and hydrogens), then the rest by
  //    element, bonds and hydrogens, preferring the candidate whose neighbours correspond
  auto same_kind = [&](int r, int p) {
    return R.element[size_t(r)] == P.element[size_t(p)] && R.degree(r) == P.degree(p) && R.hcount[size_t(r)] == P.hcount[size_t(p)] &&
           R.aromatic[size_t(r)] == P.aromatic[size_t(p)];
  };
  auto score = [&](int r, int p) {
    int k = 0;
    for (int v : R.nb[size_t(r)]) if (to[size_t(v)] >= 0 && P.nb[size_t(p)].count(to[size_t(v)])) ++k;
    return k;
  };
  for (bool grew = true; grew;) {
    grew = false;
    for (int r = 0; r < nr; ++r) {
      if (!R.heavy(r) || to[size_t(r)] >= 0) continue;
      int best = -1, bs = 0;
      bool tie = false;
      for (int p = 0; p < np; ++p) {
        if (!P.heavy(p) || from[size_t(p)] >= 0 || P.map[size_t(p)] > 0 || !same_kind(r, p)) continue;
        const int sc = score(r, p);
        if (sc > bs) best = p, bs = sc, tie = false;
        else if (sc == bs && sc > 0) tie = true;
      }
      if (best >= 0 && bs > 0 && !tie) { pair(r, best); grew = true; }
    }
    if (grew) continue;
    // nothing anchored: the first unmapped atom with a unique candidate of its kind (symmetric groups are equivalent)
    for (int r = 0; r < nr && !grew; ++r) {
      if (!R.heavy(r) || to[size_t(r)] >= 0 || R.map[size_t(r)] > 0) continue;
      int best = -1, bs = -1;
      for (int p = 0; p < np; ++p) {
        if (!P.heavy(p) || from[size_t(p)] >= 0 || P.map[size_t(p)] > 0 || !same_kind(r, p)) continue;
        const int sc = score(r, p);
        if (sc > bs) best = p, bs = sc;
      }
      if (best >= 0) { pair(r, best); grew = true; }
    }
  }
  // a mapped atom that is missing from the products leaves; an unmapped product atom of its element (water written "O")
  // is where it went
  for (int r = 0; r < nr; ++r) {
    if (!R.heavy(r) || to[size_t(r)] >= 0) continue;
    for (int p = 0; p < np; ++p)
      if (P.heavy(p) && from[size_t(p)] < 0 && P.element[size_t(p)] == R.element[size_t(r)] && P.map[size_t(p)] == 0) { pair(r, p); break; }
  }
  for (int p = 0; p < np; ++p)
    if (P.heavy(p) && from[size_t(p)] < 0)
      throw ReactError("product atom " + std::to_string(p + 1) + " (element " + std::to_string(P.element[size_t(p)]) + ") has no reactant atom: check the scheme's balance");
  // 3. bonds between heavy atoms that form and break
  std::vector<std::pair<int, int>> form, brk;
  for (int r = 0; r < nr; ++r) {
    if (!R.heavy(r) || to[size_t(r)] < 0) continue;
    for (int r2 = r + 1; r2 < nr; ++r2) {
      if (!R.heavy(r2) || to[size_t(r2)] < 0) continue;
      const bool before = R.nb[size_t(r)].count(r2) > 0, after = P.nb[size_t(to[size_t(r)])].count(to[size_t(r2)]) > 0;
      if (before && !after) brk.push_back({r, r2});
      if (!before && after) form.push_back({r, r2});
    }
  }
  std::vector<int> deleted;
  for (int r = 0; r < nr; ++r) if (R.heavy(r) && to[size_t(r)] < 0) deleted.push_back(r);
  for (auto [a, b] : form)
    if (R.map[size_t(a)] == 0 && R.map[size_t(b)] == 0) throw ReactError("a bond forms between two unmapped atoms: tag the reacting atoms with map numbers");
  for (auto [a, b] : brk)
    if (R.map[size_t(a)] == 0 && R.map[size_t(b)] == 0) throw ReactError("a bond breaks between two unmapped atoms: tag the reacting atoms with map numbers");
  {   // two molecules meet per step
    std::set<int> comps;
    for (auto [a, b] : form) comps.insert(R.comp[size_t(a)]), comps.insert(R.comp[size_t(b)]);
    for (auto [a, b] : brk) comps.insert(R.comp[size_t(a)]), comps.insert(R.comp[size_t(b)]);
    if (comps.size() > 2)
      throw ReactError(std::to_string(comps.size()) + " molecules react at once: CAPS runs one bond-forming encounter of two at a time — write it as steps");
  }
  // 4. hydrogens: atoms that lose them give them to atoms that gain them (a move); the rest leave with the leaving atoms
  std::vector<int> donors, acceptors;   // one entry per hydrogen
  for (int r = 0; r < nr; ++r) {
    if (!R.heavy(r) || to[size_t(r)] < 0) continue;
    const int d = R.hcount[size_t(r)] - P.hcount[size_t(to[size_t(r)])];
    for (int k = 0; k < d; ++k) donors.push_back(r);
    for (int k = 0; k < -d; ++k) acceptors.push_back(r);
  }
  std::vector<std::pair<int, int>> moves;   // (hydrogen, new partner)
  std::map<int, int> h_used;
  auto take_h = [&](int heavy) {
    const auto& hs = R.hs[size_t(heavy)];
    const int k = h_used[heavy]++;
    if (k >= int(hs.size())) throw ReactError("internal: atom has too few hydrogens to give");
    return hs[size_t(k)];
  };
  // a donor that joins the acceptor's product molecule pairs first (a proton moving within the new link)
  std::vector<char> donor_used(donors.size(), 0);
  for (int acc : acceptors) {
    size_t pick = donors.size();
    for (size_t k = 0; k < donors.size() && pick == donors.size(); ++k)
      if (!donor_used[k] && P.comp[size_t(to[size_t(donors[k])])] == P.comp[size_t(to[size_t(acc)])]) pick = k;
    for (size_t k = 0; k < donors.size() && pick == donors.size(); ++k)
      if (!donor_used[k]) pick = k;
    if (pick == donors.size()) throw ReactError("an atom gains a hydrogen no atom gives: check the scheme's hydrogens");
    donor_used[pick] = 1;
    moves.push_back({take_h(donors[pick]), acc});
  }
  std::vector<int> leaving_h;   // hydrogens that leave with a leaving atom (HCl, HBr) or alone
  for (size_t k = 0; k < donors.size(); ++k)
    if (!donor_used[k]) leaving_h.push_back(take_h(donors[k]));
  // 5. byproducts: small product molecules apart from the largest — water, alcohols, CO2 — made of reactant atoms that stay
  std::vector<int> psize(size_t(P.parts), 0);
  for (int p = 0; p < np; ++p) ++psize[size_t(P.comp[size_t(p)])];
  const int largest = int(std::max_element(psize.begin(), psize.end()) - psize.begin());
  std::set<int> byproduct;
  for (int r = 0; r < nr; ++r) {
    if (!R.heavy(r) || to[size_t(r)] < 0) continue;
    const int c = P.comp[size_t(to[size_t(r)])];
    if (P.parts > 1 && c != largest && psize[size_t(c)] <= 12) byproduct.insert(r);
  }
  // a leaving heavy atom (Cl, Br) and leftover hydrogens leave together as one byproduct molecule (HCl) when they can
  std::vector<std::pair<int, int>> bp_bonds;   // bonds formed inside the byproduct
  if (!deleted.empty()) {
    for (int h : leaving_h) bp_bonds.push_back({h, deleted.front()});
    for (int d : deleted) byproduct.insert(d);
    for (int h : leaving_h) byproduct.insert(h);
    leaving_h.clear();
  }
  if (!leaving_h.empty()) {   // hydrogens alone: H2
    if (leaving_h.size() % 2) throw ReactError("an odd hydrogen leaves on its own: check the scheme");
    for (size_t k = 0; k + 1 < leaving_h.size(); k += 2) bp_bonds.push_back({leaving_h[k], leaving_h[k + 1]});
    for (int h : leaving_h) byproduct.insert(h);
  }
  // a byproduct's own hydrogens go with it
  for (int b : std::vector<int>(byproduct.begin(), byproduct.end()))
    if (R.heavy(b))
      for (int h : R.hs[size_t(b)]) {
        bool moved = false;
        for (auto [mh, to_] : moves) moved = moved || mh == h;
        if (!moved) byproduct.insert(h);
      }
  for (auto [h, acc] : moves) if (byproduct.count(acc)) byproduct.insert(h);

  // 6. the pattern: every atom whose bonds or hydrogens change, the hydrogens that move or leave, the byproduct, and the heavy
  //    atoms bonded to the changing ones
  std::set<int> changed;
  for (auto [a, b] : form) changed.insert(a), changed.insert(b);
  for (auto [a, b] : brk) changed.insert(a), changed.insert(b);
  for (int d : donors) changed.insert(d);
  for (int a : acceptors) changed.insert(a);
  for (int b : byproduct) if (R.heavy(b)) changed.insert(b);
  if (changed.empty()) throw ReactError("nothing changes between the reactants and the products");
  std::set<int> pattern(changed.begin(), changed.end());
  for (int c : changed)
    for (int v : R.nb[size_t(c)]) if (R.heavy(v)) pattern.insert(v);
  for (auto [h, a] : moves) pattern.insert(h);
  for (int b : byproduct) pattern.insert(b);
  // initiators: the first bond formed between two reactant molecules, else the first formed
  std::pair<int, int> init{-1, -1};
  for (auto [a, b] : form) if (R.comp[size_t(a)] != R.comp[size_t(b)]) { init = {a, b}; break; }
  if (init.first < 0 && !form.empty()) init = form.front();
  if (init.first < 0) throw ReactError("no bond forms: CAPS templates start from a new bond");
  // every pattern atom must hang from an initiator through pattern bonds: only two molecules can meet in one step
  {
    std::set<int> seen{init.first, init.second};
    std::vector<int> st{init.first, init.second};
    while (!st.empty()) {
      const int u = st.back();
      st.pop_back();
      for (int v : R.nb[size_t(u)]) if (pattern.count(v) && seen.insert(v).second) st.push_back(v);
    }
    for (int a : pattern)
      if (!seen.count(a)) {
        std::set<int> comps;
        for (int x : pattern) comps.insert(R.comp[size_t(x)]);
        throw ReactError(std::to_string(comps.size()) + " molecules react at once: CAPS runs one bond-forming encounter of two at a time — write it as steps");
      }
  }
  // order: from the initiators outwards, each atom bonded to one before it
  std::vector<int> order;
  {
    std::set<int> seen;
    for (int s0 : {init.first, init.second}) {
      if (seen.count(s0)) continue;
      std::vector<int> q{s0};
      seen.insert(s0);
      for (size_t h = 0; h < q.size(); ++h)
        for (int v : R.nb[size_t(q[h])]) if (pattern.count(v) && seen.insert(v).second) q.push_back(v);
      order.insert(order.end(), q.begin(), q.end());
    }
  }
  std::map<int, int> mapno;
  for (size_t k = 0; k < order.size(); ++k) mapno[order[k]] = int(k + 1);
  ReactionTemplate t;
  t.name = name;
  for (int a : order) {
    TemplateAtom x;
    x.map = mapno[a];
    x.element = R.element[size_t(a)];
    if (R.heavy(a)) {
      x.degree = R.degree(a);
      if (changed.count(a) || byproduct.count(a)) x.h_min = x.h_max = R.hcount[size_t(a)];
      x.ring3 = R.ring3[size_t(a)] != 0;
      if (x.element == 6) (R.aromatic[size_t(a)] ? x.aromatic : x.not_aromatic) = true;
    }
    for (int v : R.nb[size_t(a)]) if (mapno.count(v) && mapno[v] < x.map) x.bonded.push_back(mapno[v]);
    t.atoms.push_back(x);
  }
  for (auto& a : t.atoms)
    for (int m : std::vector<int>(a.bonded))
      for (auto& b : t.atoms) if (b.map == m) b.bonded.push_back(a.map);
  t.init_a = mapno[init.first];
  t.init_b = mapno[init.second];
  t.capture = 5.0;
  // two molecules meet (min_path 0: different molecules only); within one molecule (a ring closing): any two atoms of it
  t.min_path = R.comp[size_t(init.first)] == R.comp[size_t(init.second)] ? 1 : 0;
  for (auto [a, b] : form) t.form.push_back({mapno[a], mapno[b]});
  for (auto [a, b] : brk) t.brk.push_back({mapno[a], mapno[b]});
  for (auto [h, a] : moves) t.move.push_back({mapno[h], mapno[a]});
  for (auto [a, b] : bp_bonds) t.form.push_back({mapno[a], mapno[b]});
  for (int b : byproduct) t.byproduct.push_back(mapno[b]);
  std::sort(t.byproduct.begin(), t.byproduct.end());
  t.sites = {t.init_a};
  if (notes) {
    notes->push_back(std::to_string(form.size()) + " bonds formed, " + std::to_string(brk.size()) + " broken, " + std::to_string(moves.size()) +
                     " hydrogens moved" + (byproduct.empty() ? "" : ", " + std::to_string(byproduct.size()) + " atoms leave as a byproduct"));
    notes->push_back("pattern: the changing atoms by element, bonds and hydrogens, their neighbours by element and bonds — the model compound's other atoms stand for the chain");
  }
  // through text, so the parser checks it as any template
  return parse_templates(template_text(t)).front();
}

}  // namespace caps
