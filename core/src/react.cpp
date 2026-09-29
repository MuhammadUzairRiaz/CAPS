// CAPS React: template matching, bond edits, Polymatic-style cycles and network analysis.
#include "caps/uff.hpp"
#include "caps/react.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <functional>
#include <limits>
#include <map>
#include <numeric>
#include <random>
#include <set>
#include <sstream>

#include "caps/analysis.hpp"
#include "caps/dynamics.hpp"
#include "caps/elements.hpp"
#include "caps/relax.hpp"
#include "cell_list.hpp"

namespace caps {

// ---------------------------------------------------------------------------------------------------------------
// Templates

namespace {

const std::map<std::string, std::string>& builtins() {
  static const std::map<std::string, std::string> m = {
      {"cc_crosslink",
       "reaction cc_crosslink   # C–C crosslink between saturated carbons, H2 released (radiation crosslinking)\n"
       "atom 1 C degree=4 H>=1 not_aromatic\n"
       "atom 2 C degree=4 H>=1 not_aromatic\n"
       "atom 3 H bonded 1\n"
       "atom 4 H bonded 2\n"
       "initiators 3 4   # the leaving hydrogens: close H pairs mean C–H bonds pointing at each other\n"
       "capture 3.0\n"
       "probability 1.0\n"
       "min_path 6\n"
       "form 1 2\n"
       "delete 3 4\n"
       "sites 1\n"},
      {"sulfur_allylic",
       "reaction sulfur_allylic   # accelerated sulfur cure (disulfide-donor model): an allylic C–H of a diene rubber and an\n"
       "                          # S–H end of H–S–S–H (or a growing C–S–S–H) form C–S, H2 released; two in turn give C–S–S–C\n"
       "atom 1 C degree=4 H>=1 not_aromatic   # the allylic carbon\n"
       "atom 2 C degree=3 not_aromatic bonded 1   # its double-bond neighbour\n"
       "atom 3 H bonded 1\n"
       "atom 4 S degree=2 H=1   # an S–H end of the sulfur donor\n"
       "atom 5 H bonded 4\n"
       "initiators 1 4\n"
       "capture 5.0\n"
       "probability 1.0\n"
       "min_path 0\n"
       "form 1 4\n"
       "delete 3 5\n"
       "sites 4   # the donor's S–H ends\n"},
      {"polysulfide_allylic",
       "reaction polysulfide_allylic   # a silane's polysulfide (TESPT, TESPD) couples to a diene rubber: an S–S bond opens, one\n"
       "                               # sulfur bonds to an allylic carbon, that carbon's hydrogen moves to the other (S–H)\n"
       "atom 1 C degree=4 H>=1 not_aromatic   # the allylic carbon\n"
       "atom 2 C degree=3 not_aromatic bonded 1   # its double-bond neighbour\n"
       "atom 3 H bonded 1\n"
       "atom 4 S degree=2\n"
       "atom 5 S degree=2 bonded 4\n"
       "initiators 1 4\n"
       "capture 5.0\n"
       "probability 1.0\n"
       "min_path 0\n"
       "form 1 4\n"
       "break 4 5\n"
       "move 3 5\n"
       "sites 4 5   # one S–S bond\n"},
      {"peroxide_allylic",
       "reaction peroxide_allylic   # peroxide cure of a diene rubber: C–C crosslink between allylic carbons, H2 released\n"
       "atom 1 C degree=4 H>=1 not_aromatic\n"
       "atom 2 C degree=3 not_aromatic bonded 1\n"
       "atom 3 H bonded 1\n"
       "atom 4 C degree=4 H>=1 not_aromatic\n"
       "atom 5 C degree=3 not_aromatic bonded 4\n"
       "atom 6 H bonded 4\n"
       "initiators 3 6   # the leaving hydrogens\n"
       "capture 3.0\n"
       "probability 1.0\n"
       "min_path 6\n"
       "form 1 4\n"
       "delete 3 6\n"
       "sites 1\n"},
      {"epoxy_amine_primary",
       "reaction epoxy_amine_primary   # epoxide CH2 + primary amine → β-hydroxy secondary amine\n"
       "atom 1 C ring3 H=2\n"
       "atom 2 O ring3 bonded 1\n"
       "atom 3 C ring3 bonded 1 2\n"
       "atom 4 N H=2\n"
       "atom 5 H bonded 4\n"
       "initiators 4 1\n"
       "capture 4.5\n"
       "probability 1.0\n"
       "min_path 0\n"
       "form 4 1\n"
       "break 1 2\n"
       "move 5 2\n"
       "sites 1 2 3   # one epoxide ring\n"},
      {"epoxy_amine_secondary",
       "reaction epoxy_amine_secondary   # epoxide CH2 + secondary amine → β-hydroxy tertiary amine\n"
       "atom 1 C ring3 H=2\n"
       "atom 2 O ring3 bonded 1\n"
       "atom 3 C ring3 bonded 1 2\n"
       "atom 4 N H=1 degree=3\n"
       "atom 5 H bonded 4\n"
       "initiators 4 1\n"
       "capture 4.5\n"
       "probability 1.0\n"
       "min_path 0\n"
       "form 4 1\n"
       "break 1 2\n"
       "move 5 2\n"
       "sites 1 2 3   # one epoxide ring\n"},
  };
  return m;
}

}  // namespace

std::vector<std::string> builtin_template_names() {
  std::vector<std::string> n;
  for (const auto& [k, v] : builtins()) n.push_back(k);
  return n;
}

std::string builtin_template(const std::string& name) {
  auto it = builtins().find(name);
  if (it == builtins().end()) throw ReactError("no built-in template '" + name + "'");
  return it->second;
}

std::vector<ReactionTemplate> parse_templates(const std::string& text) {
  std::vector<ReactionTemplate> out;
  std::istringstream in(text);
  std::string line;
  int lineno = 0;
  ReactionTemplate* t = nullptr;
  auto bad = [&](const std::string& why) { return ReactError("template line " + std::to_string(lineno) + ": " + why); };
  auto find = [&](int map) -> TemplateAtom* {
    for (auto& a : t->atoms)
      if (a.map == map) return &a;
    return nullptr;
  };
  while (std::getline(in, line)) {
    ++lineno;
    const std::string raw = line;
    if (auto h = line.find('#'); h != std::string::npos) line = line.substr(0, h);
    std::istringstream ls(line);
    std::vector<std::string> w;
    for (std::string x; ls >> x;) w.push_back(x);
    if (w.empty()) continue;
    auto integer = [&](size_t i) {
      if (i >= w.size()) throw bad("a number is missing after '" + w[0] + "'");
      try { return std::stoi(w[i]); } catch (...) { throw bad("'" + w[i] + "' is not a whole number"); }
    };
    auto real = [&](size_t i) {
      if (i >= w.size()) throw bad("a number is missing after '" + w[0] + "'");
      try { return std::stod(w[i]); } catch (...) { throw bad("'" + w[i] + "' is not a number"); }
    };
    if (w[0] == "reaction") {
      out.emplace_back();
      t = &out.back();
      t->name = w.size() > 1 ? w[1] : "reaction " + std::to_string(out.size());
      t->text = raw + "\n";
      continue;
    }
    if (!t) throw bad("statements must follow a 'reaction' line");
    t->text += raw + "\n";
    if (w[0] == "atom") {
      TemplateAtom a;
      a.map = integer(1);
      if (a.map <= 0 || find(a.map)) throw bad("map numbers must be positive and unique");
      if (w.size() < 3) throw bad("atom needs an element");
      a.element = element_from_symbol(w[2]);
      if (a.element <= 0) throw bad("unknown element '" + w[2] + "'");
      for (size_t i = 3; i < w.size(); ++i) {
        const std::string& k = w[i];
        if (k == "ring3") a.ring3 = true;
        else if (k == "not_aromatic") a.not_aromatic = true;
        else if (k.rfind("H>=", 0) == 0) a.h_min = std::stoi(k.substr(3));
        else if (k.rfind("H<=", 0) == 0) a.h_max = std::stoi(k.substr(3));
        else if (k.rfind("H=", 0) == 0) a.h_min = a.h_max = std::stoi(k.substr(2));
        else if (k.rfind("degree=", 0) == 0) a.degree = std::stoi(k.substr(7));
        else if (k == "bonded") {
          for (++i; i < w.size(); ++i) {
            const int m = integer(i);
            TemplateAtom* o = find(m);
            if (!o) throw bad("atom " + std::to_string(a.map) + " is bonded to " + std::to_string(m) + ", which is not defined above");
            a.bonded.push_back(m);
          }
        } else throw bad("does not understand '" + k + "'");
      }
      t->atoms.push_back(a);
      for (int m : a.bonded) find(m)->bonded.push_back(a.map);
    } else if (w[0] == "initiators") {
      t->init_a = integer(1);
      t->init_b = integer(2);
    } else if (w[0] == "capture") t->capture = real(1);
    else if (w[0] == "probability") t->probability = real(1);
    else if (w[0] == "min_path") t->min_path = integer(1);
    else if (w[0] == "charges") {   // after the reaction: keep (the atoms' charges, a deleted atom's to its partner) or forcefield
      if (w.size() < 2 || (w[1] != "keep" && w[1] != "forcefield")) throw bad("charges keep | forcefield");
      t->keep_charges = w[1] == "keep";
    }
    else if (w[0] == "sites") { for (size_t i = 1; i < w.size(); ++i) t->sites.push_back(integer(i)); }
    else if (w[0] == "form") t->form.push_back({integer(1), integer(2)});
    else if (w[0] == "break") t->brk.push_back({integer(1), integer(2)});
    else if (w[0] == "move") t->move.push_back({integer(1), integer(2)});
    else if (w[0] == "delete") { for (size_t i = 1; i < w.size(); ++i) t->remove.push_back(integer(i)); }
    else throw bad("does not understand '" + w[0] + "'");
  }
  for (auto& r : out) {
    t = &r;
    auto need = [&](int m, const std::string& what) {
      if (!find(m)) throw ReactError("reaction " + r.name + ": " + what + " refers to atom " + std::to_string(m) + ", which is not defined");
    };
    if (!r.init_a || !r.init_b) throw ReactError("reaction " + r.name + ": 'initiators a b' is missing");
    need(r.init_a, "initiators");
    need(r.init_b, "initiators");
    if (r.init_a == r.init_b) throw ReactError("reaction " + r.name + ": the initiators must be two different atoms");
    if (r.sites.empty()) r.sites.push_back(r.init_a);
    for (int m : r.sites) need(m, "sites");
    for (auto [a, b] : r.form) { need(a, "form"); need(b, "form"); }
    for (auto [a, b] : r.brk) {
      need(a, "break"); need(b, "break");
      if (std::find(find(a)->bonded.begin(), find(a)->bonded.end(), b) == find(a)->bonded.end())
        throw ReactError("reaction " + r.name + ": break " + std::to_string(a) + " " + std::to_string(b) + " is not a bond of the pattern");
    }
    for (auto [a, b] : r.move) { need(a, "move"); need(b, "move"); }
    for (int a : r.remove) need(a, "delete");
    if (r.capture <= 0 || r.probability < 0 || r.probability > 1) throw ReactError("reaction " + r.name + ": capture must be > 0, probability in [0, 1]");
    // every atom must connect to an initiator through pattern bonds
    std::set<int> reach{r.init_a, r.init_b};
    for (bool grew = true; grew;) {
      grew = false;
      for (const auto& a : r.atoms)
        if (!reach.count(a.map))
          for (int m : a.bonded)
            if (reach.count(m)) { reach.insert(a.map); grew = true; break; }
    }
    if (reach.size() != r.atoms.size()) throw ReactError("reaction " + r.name + ": every pattern atom must be bonded (through the pattern) to an initiator");
  }
  if (out.empty()) throw ReactError("no 'reaction' found");
  return out;
}

// ---------------------------------------------------------------------------------------------------------------
// Matching

namespace {

struct Chem {
  std::vector<std::vector<uint32_t>> nb;
  std::vector<int> hcount;
  std::vector<char> ring3, aromatic;
  std::vector<int> mol;
};

Chem chem_of(const System& s) {
  Chem c;
  const size_t n = s.atoms.size();
  c.nb = s.neighbours();
  c.hcount.assign(n, 0);
  c.ring3.assign(n, 0);
  c.aromatic.assign(n, 0);
  for (size_t i = 0; i < n; ++i)
    for (uint32_t j : c.nb[i]) c.hcount[i] += s.atoms[j].element == 1;
  for (size_t i = 0; i < n; ++i)
    for (size_t x = 0; x < c.nb[i].size(); ++x)
      for (size_t y = x + 1; y < c.nb[i].size(); ++y) {
        const uint32_t u = c.nb[i][x], v = c.nb[i][y];
        if (std::find(c.nb[u].begin(), c.nb[u].end(), v) != c.nb[u].end()) c.ring3[i] = c.ring3[u] = c.ring3[v] = 1;
      }
  // aromatic: three-connected carbons in a six-ring of three-connected carbons
  std::vector<char> tri(n);
  for (size_t i = 0; i < n; ++i) tri[i] = s.atoms[i].element == 6 && c.nb[i].size() == 3;
  for (size_t i = 0; i < n; ++i) {
    if (!tri[i] || c.aromatic[i]) continue;
    std::vector<uint32_t> path{uint32_t(i)};
    std::function<bool()> dfs = [&]() -> bool {
      const uint32_t u = path.back();
      for (uint32_t v : c.nb[u]) {
        if (!tri[v]) continue;
        if (path.size() == 6) { if (v == i) return true; continue; }
        if (std::find(path.begin(), path.end(), v) != path.end()) continue;
        path.push_back(v);
        if (dfs()) return true;
        path.pop_back();
      }
      return false;
    };
    if (dfs())
      for (uint32_t v : path) c.aromatic[v] = 1;
  }
  c.mol = s.molecules();
  return c;
}

bool fits(const TemplateAtom& t, const System& s, const Chem& c, uint32_t i) {
  if (s.atoms[i].element != t.element) return false;
  if (t.h_min >= 0 && c.hcount[i] < t.h_min) return false;
  if (t.h_max >= 0 && c.hcount[i] > t.h_max) return false;
  if (t.degree >= 0 && int(c.nb[i].size()) != t.degree) return false;
  if (t.ring3 && !c.ring3[i]) return false;
  if (t.not_aromatic && c.aromatic[i]) return false;
  return true;
}

bool bonded(const Chem& c, uint32_t a, uint32_t b) { return std::find(c.nb[a].begin(), c.nb[a].end(), b) != c.nb[a].end(); }

// Complete the pattern from the two initiators by backtracking over neighbours.
bool complete(const ReactionTemplate& t, const System& s, const Chem& c, const std::vector<int>& order, std::vector<int64_t>& asg,
              std::vector<char>& used, size_t k) {
  if (k == order.size()) return true;
  const int ti = order[k];
  const TemplateAtom& ta = t.atoms[ti];
  auto index_of = [&](int map) {
    for (size_t q = 0; q < t.atoms.size(); ++q)
      if (t.atoms[q].map == map) return int(q);
    return -1;
  };
  // an already assigned partner to grow from
  int anchor = -1;
  for (int m : ta.bonded) {
    const int q = index_of(m);
    if (asg[q] >= 0) { anchor = q; break; }
  }
  if (anchor < 0) return false;
  for (uint32_t cand : c.nb[asg[anchor]]) {
    if (used[cand] || !fits(ta, s, c, cand)) continue;
    bool ok = true;
    for (int m : ta.bonded) {
      const int q = index_of(m);
      if (asg[q] >= 0 && !bonded(c, uint32_t(asg[q]), cand)) { ok = false; break; }
    }
    if (!ok) continue;
    asg[ti] = cand;
    used[cand] = 1;
    if (complete(t, s, c, order, asg, used, k + 1)) return true;
    asg[ti] = -1;
    used[cand] = 0;
  }
  return false;
}

std::vector<int> growth_order(const ReactionTemplate& t) {
  std::vector<int> order;
  std::set<int> have{t.init_a, t.init_b};
  std::vector<char> done(t.atoms.size(), 0);
  for (size_t q = 0; q < t.atoms.size(); ++q)
    if (t.atoms[q].map == t.init_a || t.atoms[q].map == t.init_b) done[q] = 1;
  for (bool grew = true; grew;) {
    grew = false;
    for (size_t q = 0; q < t.atoms.size(); ++q) {
      if (done[q]) continue;
      for (int m : t.atoms[q].bonded)
        if (have.count(m)) {
          order.push_back(int(q));
          have.insert(t.atoms[q].map);
          done[q] = 1;
          grew = true;
          break;
        }
    }
  }
  return order;
}

// Are a and b within `limit` − 1 bonds of each other?
bool within_path(const Chem& c, uint32_t a, uint32_t b, int limit) {
  if (limit <= 1) return false;
  std::vector<std::pair<uint32_t, int>> q{{a, 0}};
  std::set<uint32_t> seen{a};
  for (size_t h = 0; h < q.size(); ++h) {
    auto [u, d] = q[h];
    if (u == b) return true;
    if (d + 1 >= limit) continue;
    for (uint32_t v : c.nb[u])
      if (seen.insert(v).second) q.push_back({v, d + 1});
  }
  return false;
}

}  // namespace

std::vector<Match> find_matches(const System& s, const ReactionTemplate& t, int reaction_index) {
  const Chem c = chem_of(s);
  const size_t n = s.atoms.size();
  auto idx = [&](int map) {
    for (size_t q = 0; q < t.atoms.size(); ++q)
      if (t.atoms[q].map == map) return int(q);
    return -1;
  };
  const int ia = idx(t.init_a), ib = idx(t.init_b);
  std::vector<char> ca(n, 0), cb(n, 0);
  for (size_t i = 0; i < n; ++i) {
    ca[i] = fits(t.atoms[ia], s, c, uint32_t(i));
    cb[i] = fits(t.atoms[ib], s, c, uint32_t(i));
  }
  const std::vector<int> order = growth_order(t);
  Grid g(s, t.capture);
  const double cap2 = t.capture * t.capture;
  std::vector<Match> out;
  for (uint32_t a = 0; a < n; ++a) {
    if (!ca[a]) continue;
    const int x = g.bin(g.frac[a], 0), y = g.bin(g.frac[a], 1), z = g.bin(g.frac[a], 2);
    g.for_neighbour_bins(x, y, z, [&](const std::vector<uint32_t>& bin) {
      for (uint32_t b : bin) {
        if (b == a || !cb[b]) continue;
        const Vec3 d = g.sep(a, b);
        const double r2 = dot(d, d);
        if (r2 > cap2) continue;
        if (c.mol[a] == c.mol[b] && (t.min_path <= 0 || within_path(c, a, b, t.min_path))) continue;
        std::vector<int64_t> asg(t.atoms.size(), -1);
        std::vector<char> used(n, 0);
        asg[ia] = a;
        asg[ib] = b;
        used[a] = used[b] = 1;
        if (!complete(t, s, c, order, asg, used, 0)) continue;
        Match m;
        m.reaction = reaction_index;
        m.distance = std::sqrt(r2);
        for (auto v : asg) m.atoms.push_back(uint32_t(v));
        out.push_back(std::move(m));
      }
    });
  }
  std::sort(out.begin(), out.end(), [](const Match& p, const Match& q) {
    return p.distance != q.distance ? p.distance < q.distance : p.atoms < q.atoms;
  });
  return out;
}

int count_sites(const System& s, const ReactionTemplate& t) {
  const Chem c = chem_of(s);
  auto idx = [&](int map) {
    for (size_t q = 0; q < t.atoms.size(); ++q)
      if (t.atoms[q].map == map) return int(q);
    return -1;
  };
  // the site atoms, grown from the first through pattern bonds among themselves
  std::vector<int> order{idx(t.sites[0])};
  std::set<int> in{t.sites[0]};
  for (bool grew = true; grew;) {
    grew = false;
    for (int m : t.sites) {
      if (in.count(m)) continue;
      for (int b : t.atoms[idx(m)].bonded)
        if (in.count(b)) { order.push_back(idx(m)); in.insert(m); grew = true; break; }
    }
  }
  std::set<std::vector<uint32_t>> groups;
  const size_t n = s.atoms.size();
  std::vector<int64_t> asg(t.atoms.size(), -1);
  std::vector<char> used(n, 0);
  std::function<void(size_t)> grow = [&](size_t k) {
    if (k == order.size()) {
      std::vector<uint32_t> g;
      for (int q : order) g.push_back(uint32_t(asg[q]));
      std::sort(g.begin(), g.end());
      groups.insert(g);
      return;
    }
    const TemplateAtom& ta = t.atoms[order[k]];
    int anchor = -1;
    for (int b : ta.bonded)
      if (in.count(b) && asg[idx(b)] >= 0) { anchor = idx(b); break; }
    for (uint32_t cand : c.nb[asg[anchor]]) {
      if (used[cand] || !fits(ta, s, c, cand)) continue;
      bool ok = true;
      for (int b : ta.bonded)
        if (in.count(b) && asg[idx(b)] >= 0 && !bonded(c, uint32_t(asg[idx(b)]), cand)) ok = false;
      if (!ok) continue;
      asg[order[k]] = cand;
      used[cand] = 1;
      grow(k + 1);
      asg[order[k]] = -1;
      used[cand] = 0;
    }
  };
  for (uint32_t i = 0; i < n; ++i) {
    if (!fits(t.atoms[order[0]], s, c, i)) continue;
    asg[order[0]] = i;
    used[i] = 1;
    grow(1);
    asg[order[0]] = -1;
    used[i] = 0;
  }
  return int(groups.size());
}

int apply_matches(System& s, const std::vector<ReactionTemplate>& templates, const std::vector<Match>& matches) {
  const size_t n = s.atoms.size();
  std::set<std::pair<uint32_t, uint32_t>> bonds;
  for (const auto& b : s.bonds) bonds.insert({std::min(b.i, b.j), std::max(b.i, b.j)});
  auto key = [](uint32_t a, uint32_t b) { return std::make_pair(std::min(a, b), std::max(a, b)); };
  std::vector<char> used(n, 0), dead(n, 0);
  int applied = 0;
  // charges kept through the reaction (templates with "charges keep", on a structure that has charges): a deleted atom's
  // charge goes to its bonded partner that stays, so every molecule's net charge is conserved
  bool keep = s.has_charges;
  std::vector<double> q(n);
  for (size_t i = 0; i < n; ++i) q[i] = s.atoms[i].charge;
  std::vector<std::pair<uint32_t, std::vector<uint32_t>>> gone;   // deleted atom, its partners before the reaction
  for (const auto& m : matches) {
    bool free = true;
    for (uint32_t a : m.atoms) free = free && !used[a];
    if (!free) continue;
    const ReactionTemplate& t = templates[m.reaction];
    auto at = [&](int map) {
      for (size_t q = 0; q < t.atoms.size(); ++q)
        if (t.atoms[q].map == map) return m.atoms[q];
      throw ReactError("internal: map " + std::to_string(map) + " not in the match");
    };
    for (uint32_t a : m.atoms) used[a] = 1;
    keep = keep && t.keep_charges;
    for (int a : t.remove) {
      std::vector<uint32_t> partners;
      for (const auto& e : bonds)
        if (e.first == at(a) || e.second == at(a)) partners.push_back(e.first == at(a) ? e.second : e.first);
      gone.push_back({at(a), partners});
    }
    for (auto [a, b] : t.brk) bonds.erase(key(at(a), at(b)));
    for (auto [a, b] : t.form) bonds.insert(key(at(a), at(b)));
    for (auto [h, x] : t.move) {
      const uint32_t hi = at(h), xi = at(x);
      for (auto it = bonds.begin(); it != bonds.end();) it = (it->first == hi || it->second == hi) ? bonds.erase(it) : std::next(it);
      bonds.insert(key(hi, xi));
      // put the moved atom at a bond length from its new partner, on the side it came from
      Vec3 d = s.cell.valid() ? s.cell.minimum_image(s.atoms[hi].pos - s.atoms[xi].pos) : s.atoms[hi].pos - s.atoms[xi].pos;
      double len = norm(d);
      if (len < 1e-6) { d = {1, 0, 0}; len = 1; }
      const double r0 = element(s.atoms[hi].element).covalent + element(s.atoms[xi].element).covalent;
      s.atoms[hi].pos = s.atoms[xi].pos + d * (r0 / len);
    }
    for (int a : t.remove) dead[at(a)] = 1;
    ++applied;
  }
  if (keep)
    for (const auto& [d, partners] : gone) {
      uint32_t to = d;
      for (uint32_t p : partners) if (!dead[p]) { to = p; break; }
      if (to == d) { keep = false; break; }   // nothing stays bonded to it: the charge has nowhere to go
      q[to] += q[d];
    }
  if (keep) for (size_t i = 0; i < n; ++i) s.atoms[i].charge = q[i];
  // compact: drop deleted atoms and their bonds
  std::vector<int64_t> remap(n, -1);
  std::vector<Atom> atoms;
  for (size_t i = 0; i < n; ++i)
    if (!dead[i]) { remap[i] = int64_t(atoms.size()); atoms.push_back(s.atoms[i]); }
  s.bonds.clear();
  for (auto [a, b] : bonds)
    if (remap[a] >= 0 && remap[b] >= 0) s.bonds.push_back({uint32_t(remap[a]), uint32_t(remap[b])});
  s.atoms = std::move(atoms);
  for (size_t i = 0; i < s.atoms.size(); ++i) s.atoms[i].id = int64_t(i + 1);
  s.velocities.clear();
  s.has_charges = keep && applied > 0;   // else charges are recomputed for the new chemistry
  s.bonds_from_file = true;
  // molecules follow the new bonds
  s.has_mol = false;
  const auto mol = s.molecules();
  for (size_t i = 0; i < s.atoms.size(); ++i) s.atoms[i].mol = mol[i] + 1;
  s.has_mol = true;
  return applied;
}

// ---------------------------------------------------------------------------------------------------------------
// Networks

ClusterStats cluster_stats(const System& s) {
  ClusterStats c;
  System t = s;
  t.has_mol = false;
  int nm = 0;
  const auto mol = t.molecules(&nm);
  std::vector<double> mass(nm, 0.0);
  for (size_t i = 0; i < s.atoms.size(); ++i) mass[mol[i]] += s.mass_of(s.atoms[i]);
  c.clusters = nm;
  double tot = 0, sq = 0, big = 0;
  for (double m : mass) { tot += m; sq += m * m; big = std::max(big, m); }
  if (tot <= 0) return c;
  c.largest_fraction = big / tot;
  c.mw = sq / tot;
  const double rt = tot - big, rs = sq - big * big;
  c.reduced_mw = rt > 0 ? rs / rt : 0;
  return c;
}

double flory_stockmayer(double r, double fa, double fb) {
  const double d = r * (fa - 1) * (fb - 1);
  return d > 0 ? 1.0 / std::sqrt(d) : std::numeric_limits<double>::infinity();
}

void react(System& s, const ReactOptions& o, ReactReport* rep_out) {
  const auto t0 = std::chrono::steady_clock::now();
  ReactReport rep;
  if (o.templates.empty()) throw ReactError("no reaction templates");
  std::mt19937_64 rng(o.seed);
  std::uniform_real_distribution<double> uni(0, 1);
  for (const auto& t : o.templates) rep.initial_sites = std::max(rep.initial_sites, count_sites(s, t));
  if (rep.initial_sites == 0) throw ReactError("no atom matches the counted site of any template; check the templates against the structure");
  if (s.cell.valid() && !s.unwrapped) make_molecules_whole(s);

  auto relax_now = [&](CycleRow& row) {
    if (!o.relax) return;
    RelaxOptions r;
    r.pushoff = false;
    r.ftol = o.relax_ftol;
    r.max_iterations = o.relax_iterations;
    r.energy = o.energy;
    RelaxReport rr;
    try {
      relax(s, r, &rr);
    } catch (const FieldError& e) {
      throw ReactError(std::string("cannot relax after the reaction: ") + e.what() +
                       ". Run with relaxation off (topology only) or add parameters for these atoms.");
    }
    row.energy = rr.final.total();
    row.max_force = rr.fmax_final;
    // names follow the new types
    const ForceField ff = default_forcefield(s);
    for (size_t i = 0; i < s.atoms.size(); ++i) s.atoms[i].name = ff.atom_type[i];
  };

  if (o.during_md && o.md_ps <= 0) throw ReactError("reactions during MD need the check interval (md_ps > 0)");
  if (o.during_md) {
    char b[240];
    std::snprintf(b, sizeof b, "REACTER-style (Gissinger et al. 2017): continuous NVT at %.0f K, reactions checked every %.3g ps; reacted sites "
                  "stabilised by a local capped-force minimisation (the rest held), no global minimisation", o.temperature, o.md_ps);
    rep.notes.push_back(b);
  }
  // REACTER-style: the reacted sites (atom ids, the atoms within two bonds and within 5 Å of those) settle by a local
  // minimisation, everything else held
  auto stabilise = [&](const std::set<int64_t>& ids, CycleRow& row) {
    std::vector<char> site(s.atoms.size(), 0);
    for (size_t i = 0; i < s.atoms.size(); ++i) site[i] = ids.count(s.atoms[i].id) ? 1 : 0;
    const auto nb = s.neighbours();
    for (int ring = 0; ring < 2; ++ring) {   // and the atoms within two bonds
      std::vector<char> grow = site;
      for (size_t i = 0; i < s.atoms.size(); ++i)
        if (site[i]) for (uint32_t w : nb[i]) grow[w] = 1;
      site = grow;
    }
    {   // and every atom within 5 Å of those, whose contacts the new geometry changes
      std::vector<Vec3> core;
      for (size_t i = 0; i < s.atoms.size(); ++i) if (site[i]) core.push_back(s.atoms[i].pos);
      for (size_t i = 0; i < s.atoms.size(); ++i) {
        if (site[i]) continue;
        for (const auto& c : core)
          if (norm(s.cell.valid() ? s.cell.minimum_image(s.atoms[i].pos - c) : s.atoms[i].pos - c) < 5.0) { site[i] = 2; break; }
      }
      for (auto& x : site) x = x ? 1 : 0;
    }
    RelaxOptions r;
    r.pushoff = true;   // new bonds start stretched: capped forces first
    r.ftol = std::min(o.relax_ftol, 1.0);   // tight: dynamics goes on from here
    r.max_iterations = std::max(o.relax_iterations, 3000);
    r.energy = o.energy;
    r.fixed.assign(s.atoms.size(), 1);
    for (size_t i = 0; i < s.atoms.size(); ++i) r.fixed[i] = site[i] ? 0 : 1;
    const auto v = s.velocities;
    RelaxReport rr;
    try {
      relax(s, r, &rr);
    } catch (const FieldError& e) {
      throw ReactError(std::string("cannot stabilise the reacted sites: ") + e.what());
    }
    if (v.size() == s.atoms.size()) {   // the held atoms keep their velocities; the site starts from rest
      s.velocities = v;
      for (size_t i = 0; i < s.atoms.size(); ++i) if (site[i]) s.velocities[i] = {0, 0, 0};
    }
    row.energy = rr.final.total();
    row.max_force = rr.fmax_final;
  };

  int stall = 0;
  for (int cycle = 1; cycle <= o.max_cycles; ++cycle) {
    // the checkpoint: this cycle failing gives back the structure as the last one left it
    const System keep = o.keep_on_failure && !rep.cycles.empty() ? s : System{};
    const int reactions_before = rep.reactions, stall_before = stall;
    int applied = 0;
    try {
    std::vector<Match> all;
    for (size_t k = 0; k < o.templates.size(); ++k) {
      auto m = find_matches(s, o.templates[k], int(k));
      all.insert(all.end(), m.begin(), m.end());
    }
    std::stable_sort(all.begin(), all.end(), [](const Match& a, const Match& b) { return a.distance < b.distance; });
    // probability, then at most max_per_cycle non-overlapping, closest first
    std::vector<Match> chosen;
    std::set<uint32_t> busy;
    for (const auto& m : all) {
      if (int(chosen.size()) >= o.max_per_cycle) break;
      if (uni(rng) > o.templates[m.reaction].probability) continue;
      bool free = true;
      for (uint32_t a : m.atoms) free = free && !busy.count(a);
      if (!free) continue;
      for (uint32_t a : m.atoms) busy.insert(a);
      chosen.push_back(m);
    }
    std::set<int64_t> site_ids;
    for (const auto& m : chosen)
      for (uint32_t a : m.atoms) site_ids.insert(s.atoms[a].id);
    applied = chosen.empty() ? 0 : apply_matches(s, o.templates, chosen);
    CycleRow row;
    row.cycle = cycle;
    row.reactions = applied;
    rep.reactions += applied;
    row.total = rep.reactions;
    row.conversion = double(rep.reactions) / rep.initial_sites;
    if (applied > 0) {
      stall = 0;
      if (o.during_md) stabilise(site_ids, row);
      else relax_now(row);
    } else {
      ++stall;
    }
    if (o.md_ps > 0 && (o.during_md || applied > 0 || stall <= o.stall_cycles)) {
      DynamicsOptions d;
      d.steps = std::max<int64_t>(1, std::llround(o.md_ps * 1000));
      d.temperature = o.temperature;
      d.new_velocities = s.velocities.size() != s.atoms.size();
      d.seed = o.seed + uint64_t(cycle);
      d.thermo_every = int(d.steps);
      d.frame_every = 0;
      d.energy = o.energy;
      DynamicsReport dr;
      try {
        run_dynamics(s, d, &dr);
      } catch (const FieldError& e) {
        throw ReactError(std::string("cannot run dynamics between cycles: ") + e.what());
      }
      if (!dr.thermo.empty()) row.energy = dr.thermo.back().potential;
    }
    row.clusters = cluster_stats(s);
    row.atoms = int(s.atoms.size());
    rep.cycles.push_back(row);
    if (o.frame) o.frame(s, cycle);   // a frame that cannot be kept fails the cycle too
    } catch (const std::exception& e) {
      if (!rep.cycles.empty() && rep.cycles.back().cycle == cycle) rep.cycles.pop_back();
      if (!o.keep_on_failure || rep.cycles.empty()) throw;
      s = keep;
      rep.reactions = reactions_before;
      stall = stall_before;
      rep.failed_cycle = cycle;
      rep.failure = e.what();
      rep.notes.push_back("failed at cycle " + std::to_string(cycle) + ": " + rep.failure + " — the structure after cycle " + std::to_string(cycle - 1) +
                          " is kept (restart from it)");
      break;
    }
    const CycleRow& row = rep.cycles.back();
    if (o.progress && !o.progress(row)) throw ReactError("reaction run cancelled");
    if (row.conversion >= o.target_conversion) {
      rep.notes.push_back("target conversion reached");
      break;
    }
    if (applied == 0 && !o.during_md && (o.md_ps <= 0 || stall > o.stall_cycles)) {
      rep.notes.push_back(o.md_ps > 0 ? "no reactive pairs within the capture distance after dynamics; stopped"
                                      : "no more reactive pairs within the capture distance; stopped (enable dynamics between cycles to let "
                                        "groups diffuse)");
      break;
    }
    if (cycle == o.max_cycles) rep.notes.push_back("cycle limit reached");
  }

  // gel point: the reduced weight-average mass (without the largest cluster) peaks there
  double peak = 0;
  int at = -1;
  for (size_t k = 0; k < rep.cycles.size(); ++k)
    if (rep.cycles[k].clusters.reduced_mw > peak) { peak = rep.cycles[k].clusters.reduced_mw; at = int(k); }
  if (at >= 0 && at + 1 < int(rep.cycles.size()) && rep.cycles.back().clusters.reduced_mw < 0.8 * peak &&
      rep.cycles.back().clusters.largest_fraction > 0.5)
    rep.gel_conversion = rep.cycles[at].conversion;
  rep.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  char b[256];
  const auto& last = rep.cycles.empty() ? CycleRow{} : rep.cycles.back();
  std::snprintf(b, sizeof b, "%d reactions in %zu cycles (%d counted sites) · conversion %.3f · largest cluster %.1f%% of the mass · %.1f s",
                rep.reactions, rep.cycles.size(), rep.initial_sites, last.conversion, 100 * last.clusters.largest_fraction, rep.seconds);
  rep.notes.insert(rep.notes.begin(), b);
  if (rep.gel_conversion >= 0) {
    std::snprintf(b, sizeof b, "gel point from cluster analysis at conversion %.3f (peak of the reduced weight-average mass)", rep.gel_conversion);
    rep.notes.push_back(b);
  }
  s.unwrapped = true;
  if (rep_out) *rep_out = std::move(rep);
}

}  // namespace caps
