// The reaction-template editor's view of a template: drawings of the pattern before and after, changes and checks.
#include <algorithm>
#include <cstdio>
#include <functional>
#include <map>
#include <set>

#include "caps/elements.hpp"
#include "caps/json.hpp"
#include "caps/molecule.hpp"
#include "caps/react.hpp"

namespace caps {

namespace {
// A byproduct leaves the network: its bonds to the other atoms are gone in the product, its bonds among themselves stay.
template <class Edges>
void cut_byproduct(const ReactionTemplate& t, Edges& post) {
  auto in = [&](int m) { return std::find(t.byproduct.begin(), t.byproduct.end(), m) != t.byproduct.end(); };
  for (auto it = post.begin(); it != post.end();) it = in(it->first) != in(it->second) ? post.erase(it) : std::next(it);
}
}  // namespace


namespace {

using Edge = std::pair<int, int>;   // map numbers, low first
Edge edge(int a, int b) { return {std::min(a, b), std::max(a, b)}; }

Json drawing(const std::vector<int>& maps, const std::map<int, int>& element, const std::set<Edge>& bonds, const std::set<int>& reacting) {
  MolGraph g;
  std::map<int, int> index;
  for (int m : maps) {
    MolAtom a;
    a.element = element.at(m);
    a.bracket = true;
    a.hcount = 0;
    a.map = m;
    index[m] = int(g.atoms.size());
    g.atoms.push_back(a);
  }
  for (const auto& [a, b] : bonds) {
    if (!index.count(a) || !index.count(b)) continue;
    MolBond mb;
    mb.a = index[a];
    mb.b = index[b];
    g.bonds.push_back(mb);
    g.atoms[size_t(mb.a)].order.push_back(mb.b);
    g.atoms[size_t(mb.b)].order.push_back(mb.a);
  }
  g.heavy = int(g.atoms.size());
  const auto p = g.atoms.empty() ? std::vector<Vec3>{} : depict(g);
  Json atoms = Json::array();
  for (size_t k = 0; k < g.atoms.size(); ++k) {
    Json o = Json::object();
    o["map"] = maps[k];
    o["symbol"] = std::string(caps::element(g.atoms[k].element).symbol);
    o["x"] = p[k][0];
    o["y"] = p[k][1];
    o["reacting"] = reacting.count(maps[k]) > 0;
    atoms.push_back(std::move(o));
  }
  Json bl = Json::array();
  for (const auto& [a, b] : bonds) {
    Json e = Json::array();
    e.push_back(a);
    e.push_back(b);
    bl.push_back(std::move(e));
  }
  Json out = Json::object();
  out["atoms"] = std::move(atoms);
  out["bonds"] = std::move(bl);
  return out;
}

// One side of a reaction SMARTS: each connected part by depth-first walk with ring closures, parts joined by '.'.
// Pattern bonds carry no order (the template is connectivity), so each is '~' (any bond).
std::string smarts_side(const std::vector<int>& maps, const std::set<Edge>& bonds, const std::function<std::string(int)>& atom,
                        const std::set<Edge>& single = {}) {
  auto sym = [&](int a, int b) { return single.count(edge(a, b)) ? std::string("-") : std::string("~"); };
  std::map<int, std::vector<int>> adj;
  for (int m : maps) adj[m];
  for (const auto& [a, b] : bonds)
    if (adj.count(a) && adj.count(b)) adj[a].push_back(b), adj[b].push_back(a);
  for (auto& [m, v] : adj) std::sort(v.begin(), v.end());
  std::set<int> seen;
  std::string out;
  for (int start : maps) {
    if (seen.count(start)) continue;
    // ring closures: the bonds a spanning tree from start does not use
    std::map<int, int> parent{{start, 0}};
    std::vector<int> order, stack{start};
    std::set<int> tree;
    while (!stack.empty()) {
      const int u = stack.back();
      stack.pop_back();
      if (tree.count(u)) continue;
      tree.insert(u);
      order.push_back(u);
      for (auto it = adj[u].rbegin(); it != adj[u].rend(); ++it)
        if (!tree.count(*it)) parent[*it] = u, stack.push_back(*it);
    }
    std::set<Edge> tree_edges;
    for (const auto& [c, p] : parent) if (p) tree_edges.insert(edge(c, p));
    std::map<int, std::vector<std::pair<int, std::string>>> closures;   // atom → ring numbers (and bond symbol) opened or closed there
    int ring = 0;
    for (const auto& e : bonds)
      if (tree.count(e.first) && tree.count(e.second) && !tree_edges.count(e)) {
        ++ring;
        closures[e.first].push_back({ring, sym(e.first, e.second)}), closures[e.second].push_back({ring, sym(e.first, e.second)});
      }
    std::function<void(int, int)> walk = [&](int u, int from) {
      seen.insert(u);
      out += atom(u);
      for (const auto& [r, bs] : closures[u]) out += bs + (r < 10 ? std::to_string(r) : "%" + std::to_string(r));
      std::vector<int> kids;
      for (int w : adj[u]) if (w != from && parent.count(w) && parent[w] == u && !seen.count(w)) kids.push_back(w);
      for (size_t k = 0; k < kids.size(); ++k) {
        const bool branch = k + 1 < kids.size();
        if (branch) out += "(";
        out += sym(u, kids[k]);
        walk(kids[k], u);
        if (branch) out += ")";
      }
    };
    if (!out.empty()) out += ".";
    walk(start, 0);
  }
  return out;
}

}  // namespace

std::string reaction_smarts(const ReactionTemplate& t) {
  std::map<int, const TemplateAtom*> by;
  std::vector<int> maps;
  std::set<Edge> pre;
  for (const auto& a : t.atoms) {
    by[a.map] = &a;
    maps.push_back(a.map);
    for (int b : a.bonded) pre.insert(edge(a.map, b));
  }
  // the query side: element, connections, hydrogens, ring and aromaticity constraints, the map number
  auto query = [&](int m) {
    const TemplateAtom& a = *by.at(m);
    std::string q = "[#" + std::to_string(a.element);
    if (a.degree >= 0) q += ";X" + std::to_string(a.degree);
    if (a.h_min >= 0 || a.h_max >= 0) {
      const int lo = std::max(0, a.h_min), hi = a.h_max >= 0 ? a.h_max : 4;
      if (lo == hi) q += ";H" + std::to_string(lo);
      else if (lo >= 1 && hi == 4 && a.h_max < 0) q += lo == 1 ? ";!H0" : ";!H0;!H1" + std::string(lo >= 3 ? ";!H2" : "");
      else {
        q += ";";
        for (int h = lo; h <= hi; ++h) q += (h > lo ? "," : "") + std::string("H") + std::to_string(h);
      }
    }
    if (a.ring3) q += ";r3";
    if (a.not_aromatic) q += ";A";
    if (a.aromatic) q += ";a";
    return q + ":" + std::to_string(m) + "]";
  };
  std::set<Edge> post = pre;
  for (const auto& [a, b] : t.form) post.insert(edge(a, b));
  for (const auto& [a, b] : t.brk) post.erase(edge(a, b));
  for (const auto& [a, b] : t.move) {
    for (auto it = post.begin(); it != post.end();) it = (it->first == a || it->second == a) ? post.erase(it) : std::next(it);
    post.insert(edge(a, b));
  }
  std::vector<int> post_maps;
  for (int m : maps)
    if (std::find(t.remove.begin(), t.remove.end(), m) == t.remove.end()) post_maps.push_back(m);
  for (int m : t.remove)
    for (auto it = post.begin(); it != post.end();) it = (it->first == m || it->second == m) ? post.erase(it) : std::next(it);
  cut_byproduct(t, post);
  // the product side: element and map; a deleted atom is absent there (a mapped reactant atom missing from the products
  // is removed, as RDKit and Daylight read reaction SMARTS)
  auto product = [&](int m) { return "[#" + std::to_string(by.at(m)->element) + ":" + std::to_string(m) + "]"; };
  // a formed bond whose two atoms keep their number of bonds in the pattern (a substitution: an H or a partner lost, this
  // bond gained) is single; other bonds stay '~' (CAPS perceives the orders after the reaction)
  std::map<int, int> deg_pre, deg_post;
  for (const auto& [a, b] : pre) ++deg_pre[a], ++deg_pre[b];
  for (const auto& [a, b] : post) ++deg_post[a], ++deg_post[b];
  std::set<Edge> single;
  for (const auto& [a, b] : t.form)
    if (deg_post[a] == deg_pre[a] && deg_post[b] == deg_pre[b]) single.insert(edge(a, b));
  for (const auto& [a, b] : t.move)
    if (deg_post[a] == deg_pre[a] && deg_post[b] == deg_pre[b]) single.insert(edge(a, b));
  const std::string right = smarts_side(post_maps, post, product, single);
  return smarts_side(maps, pre, query) + ">>" + right;
}

std::string template_view(const ReactionTemplate& t) {
  std::map<int, int> element;
  std::vector<int> maps;
  std::set<Edge> pre;
  for (const auto& a : t.atoms) {
    element[a.map] = a.element;
    maps.push_back(a.map);
    for (int b : a.bonded) pre.insert(edge(a.map, b));
  }
  auto sym = [&](int m) { return (element.count(m) ? std::string(caps::element(element[m]).symbol) : std::string("?")) + std::to_string(m); };
  std::set<int> reacting{t.init_a, t.init_b};
  std::set<Edge> post = pre;
  Json changes = Json::array();
  auto change = [&](const char* kind, const std::string& text) {
    Json o = Json::object();
    o["kind"] = kind;
    o["text"] = text;
    changes.push_back(std::move(o));
  };
  std::map<int, int> delta;   // bonds gained − lost, per atom
  for (const auto& [a, b] : t.form) { post.insert(edge(a, b)); reacting.insert(a); reacting.insert(b); ++delta[a]; ++delta[b]; change("formed", sym(a) + " – " + sym(b)); }
  for (const auto& [a, b] : t.brk) { post.erase(edge(a, b)); reacting.insert(a); reacting.insert(b); --delta[a]; --delta[b]; change("broken", sym(a) + " – " + sym(b)); }
  for (const auto& [a, b] : t.move) {
    // the atom leaves its partners and bonds to the second atom
    std::vector<Edge> gone;
    for (const auto& e : post) if (e.first == a || e.second == a) gone.push_back(e);
    std::string from;
    for (const auto& e : gone) {
      post.erase(e);
      const int other = e.first == a ? e.second : e.first;
      --delta[other];
      --delta[a];
      from += (from.empty() ? "" : ", ") + sym(other);
    }
    post.insert(edge(a, b));
    ++delta[a];
    ++delta[b];
    reacting.insert(a);
    reacting.insert(b);
    change("moved", sym(a) + " from " + (from.empty() ? std::string("—") : from) + " to " + sym(b));
  }
  std::vector<int> post_maps;
  for (int m : maps) if (std::find(t.remove.begin(), t.remove.end(), m) == t.remove.end()) post_maps.push_back(m);
  for (int m : t.remove) {
    for (auto it = post.begin(); it != post.end();) it = (it->first == m || it->second == m) ? post.erase(it) : std::next(it);
    reacting.insert(m);
    change("deleted", sym(m));
  }
  if (!t.byproduct.empty()) {
    cut_byproduct(t, post);
    std::string bp;
    for (int m : t.byproduct) { bp += (bp.empty() ? "" : " ") + sym(m); reacting.insert(m); }
    change("byproduct", bp + " leave as a molecule (kept or removed by the run)");
  }

  Json checks = Json::array();
  auto check = [&](bool ok, const std::string& text) {
    Json o = Json::object();
    o["ok"] = ok;
    o["text"] = text;
    checks.push_back(std::move(o));
  };
  const std::set<int> known(maps.begin(), maps.end());
  check(known.count(t.init_a) && known.count(t.init_b) && t.init_a != t.init_b,
        "Initiators " + sym(t.init_a) + " and " + sym(t.init_b) + " are two atoms of the pattern");
  bool refs = true;
  for (const auto& v : {t.form, t.brk, t.move})
    for (const auto& [a, b] : v) refs = refs && known.count(a) && known.count(b);
  for (int m : t.remove) refs = refs && known.count(m);
  check(refs, "Every change refers to atoms of the pattern");
  bool broken_exist = true, formed_new = true;
  for (const auto& [a, b] : t.brk) broken_exist = broken_exist && pre.count(edge(a, b));
  for (const auto& [a, b] : t.form) formed_new = formed_new && !pre.count(edge(a, b));
  check(broken_exist, "Broken bonds exist before the reaction");
  check(formed_new, "Formed bonds are new");
  check(!changes.items().empty(), "The template changes something");
  const int kept = int(post_maps.size()), all = int(maps.size());
  check(true, "Mapped atoms kept on both sides: " + std::to_string(kept) + " of " + std::to_string(all) + (kept < all ? " (the rest are deleted)" : ""));
  std::string changed;
  for (const auto& [m, d] : delta) if (d != 0 && std::find(t.remove.begin(), t.remove.end(), m) == t.remove.end()) changed += (changed.empty() ? "" : ", ") + sym(m) + (d > 0 ? " +" : " ") + std::to_string(d);
  check(changed.empty(), changed.empty() ? "Every atom keeps its number of bonds" : "Bond counts change: " + changed + " (check the hydrogens)");
  check(t.capture > 0 && t.capture < 15, "Capture distance " + std::to_string(t.capture).substr(0, 4) + " Å");

  Json j = Json::object();
  j["name"] = t.name;
  j["pre"] = drawing(maps, element, pre, reacting);
  j["post"] = drawing(post_maps, element, post, reacting);
  j["changes"] = std::move(changes);
  j["checks"] = std::move(checks);
  Json init = Json::array();
  init.push_back(t.init_a);
  init.push_back(t.init_b);
  j["initiators"] = std::move(init);
  j["capture"] = t.capture;
  j["probability"] = t.probability;
  j["min_path"] = t.min_path;
  j["smarts"] = reaction_smarts(t);
  return j.dump(0);
}

}  // namespace caps
