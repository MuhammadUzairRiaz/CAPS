// The reaction-template editor's view of a template: drawings of the pattern before and after, changes and checks.
#include <algorithm>
#include <cstdio>
#include <map>
#include <set>

#include "caps/elements.hpp"
#include "caps/json.hpp"
#include "caps/molecule.hpp"
#include "caps/react.hpp"

namespace caps {

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

}  // namespace

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
  return j.dump(0);
}

}  // namespace caps
