// CAPS React: template matching, bond edits, Polymatic-style cycles and network analysis.
#include "caps/rng.hpp"
#include "caps/uff.hpp"
#include "caps/react.hpp"

#include <tuple>
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
       "reaction cc_crosslink   # C–C crosslink of two CH2 sites (PE, EPDM, polyisoprene CH2): CH2–CH2 → CH–CH; each CH2 once\n"
       "                        # two H leave (H2 in radiation crosslinking; the peroxide's alcohols in a peroxide cure — same network)\n"
       "atom 1 C degree=4 H=2 not_aromatic\n"
       "atom 2 C degree=4 H=2 not_aromatic\n"
       "atom 3 H bonded 1\n"
       "atom 4 H bonded 2\n"
       "initiators 3 4   # the leaving hydrogens: close H pairs mean C–H bonds pointing at each other\n"
       "capture 3.0\n"
       "probability 1.0\n"
       "min_path 6\n"
       "form 1 2\n"
       "form 3 4   # the two hydrogens leave together as H2\n"
       "byproduct 3 4\n"
       "sites 1\n"},
      {"cc_crosslink_any",
       "reaction cc_crosslink_any   # C–C crosslink of any two sp3 C–H carbons (CH3, CH2 or CH: PP's tertiary CH too); two H leave\n"
       "atom 1 C degree=4 H>=1 not_aromatic\n"
       "atom 2 C degree=4 H>=1 not_aromatic\n"
       "atom 3 H bonded 1\n"
       "atom 4 H bonded 2\n"
       "initiators 3 4\n"
       "capture 3.0\n"
       "probability 1.0\n"
       "min_path 6\n"
       "form 1 2\n"
       "form 3 4\n"
       "byproduct 3 4\n"
       "sites 1\n"},
      {"sulfur_allylic",
       "reaction sulfur_allylic   # accelerated sulfur cure (disulfide-donor model): an allylic C–H of a diene rubber and an\n"
       "                          # S–H end of H–S–S–H (or a growing C–S–S–H) form C–S, H2 released; two in turn give C–S–S–C\n"
       "atom 1 C degree=4 H>=1 not_aromatic   # the allylic carbon\n"
       "atom 2 C degree=3 not_aromatic bonded 1   # its double-bond neighbour\n"
       "atom 3 H bonded 1\n"
       "atom 4 S degree=2 H>=1  # an S–H end of the sulfur donor (H2S: a monosulfide C–S–C in two steps)\n"
       "atom 5 H bonded 4\n"
       "initiators 1 4\n"
       "capture 5.0\n"
       "probability 1.0\n"
       "min_path 0\n"
       "form 1 4\n"
       "form 3 5   # the two hydrogens leave as H2 (H2S in a real accelerated cure; the mass balance is the same for the network)\n"
       "byproduct 3 5\n"
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
       "form 3 6   # written as H2: in a real cure the peroxide's radicals take the hydrogens (cumyl alcohol from DCP)\n"
       "byproduct 3 6\n"
       "sites 1\n"},
      {"enr_acid_ester",
       "reaction enr_acid_ester   # epoxidised natural rubber (trisubstituted epoxide) + carboxylic acid → β-hydroxy ester: the acid\n"
       "                          # oxygen bonds to the tertiary ring carbon, the ring opens, the acid H goes to the epoxide O (an OH on\n"
       "                          # the secondary carbon). A ring opening: no water. PBS or maleic-acid COOH ends, a half-ester's COOH\n"
       "atom 1 C ring3 H=0 degree=4 not_aromatic   # the tertiary epoxide carbon\n"
       "atom 2 O ring3 bonded 1\n"
       "atom 3 C ring3 H=1 bonded 1 2              # the secondary epoxide carbon\n"
       "atom 4 C degree=3 not_aromatic             # the acid carbon\n"
       "atom 5 O degree=1 bonded 4                 # its C=O\n"
       "atom 6 O H=1 degree=2 bonded 4             # the acid OH\n"
       "atom 7 H bonded 6\n"
       "initiators 6 1\n"
       "capture 5.0\n"
       "probability 1.0\n"
       "min_path 0\n"
       "form 6 1\n"
       "break 1 2\n"
       "move 7 2\n"
       "sites 1 2 3   # one epoxide ring\n"},
      {"ester_condensation",
       "reaction ester_condensation   # carboxylic acid + alcohol → ester + water (Fischer esterification, a condensation): the acid's\n"
       "                              # OH and the alcohol's H leave as H2O. PBS ends (COOH + HO–CH2), the OH an epoxide opening left\n"
       "atom 1 C degree=3 not_aromatic   # the acid carbon\n"
       "atom 2 O degree=1 bonded 1       # its C=O\n"
       "atom 3 O H=1 degree=2 bonded 1   # the acid OH: leaves in the water\n"
       "atom 4 H bonded 3\n"
       "atom 5 O H=1 degree=2            # the alcohol oxygen: becomes the ester oxygen\n"
       "atom 6 C degree=4 bonded 5       # its sp3 carbon (not another acid)\n"
       "atom 7 H bonded 5\n"
       "initiators 5 1\n"
       "capture 5.0\n"
       "probability 1.0\n"
       "min_path 0\n"
       "form 5 1\n"
       "break 1 3\n"
       "break 5 7\n"
       "form 3 7   # H2O: the acid's O with its H and the alcohol's H\n"
       "byproduct 3 4 7\n"
       "sites 1 2 3   # one COOH\n"},
      {"anhydride_alcohol",
       "reaction anhydride_alcohol   # cyclic anhydride (maleic anhydride, MAH) + alcohol → half-ester with a free COOH: the alcohol O\n"
       "                             # bonds to one carbonyl carbon, the ring C–O–C opens, the alcohol H goes to the bridging O. No water.\n"
       "                             # The new COOH can open an epoxide in turn (enr_acid_ester): ENR–MAH–ENR, ENR–MAH–PBS\n"
       "atom 1 C degree=3 not_aromatic   # a carbonyl carbon of the anhydride\n"
       "atom 2 O degree=1 bonded 1       # its C=O\n"
       "atom 3 O H=0 degree=2 bonded 1   # the bridging oxygen\n"
       "atom 4 C degree=3 bonded 3       # the other carbonyl carbon\n"
       "atom 5 O degree=1 bonded 4       # its C=O\n"
       "atom 6 O H=1 degree=2            # the alcohol oxygen\n"
       "atom 7 H bonded 6\n"
       "atom 8 C degree=4 bonded 6       # its sp3 carbon: an alcohol, not a carboxylic acid's OH\n"
       "initiators 6 1\n"
       "capture 5.0\n"
       "probability 1.0\n"
       "min_path 0\n"
       "form 6 1\n"
       "break 1 3\n"
       "move 7 3\n"
       "sites 1 2 3 4 5   # one anhydride group\n"},
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
        else if (k == "aromatic") a.aromatic = true;
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
    else if (w[0] == "byproduct") { for (size_t i = 1; i < w.size(); ++i) t->byproduct.push_back(integer(i)); }
    else if (w[0] == "weight") t->weight = real(1);
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
    for (int a : r.byproduct) {
      need(a, "byproduct");
      if (std::find(r.remove.begin(), r.remove.end(), a) != r.remove.end()) throw ReactError("reaction " + r.name + ": atom " + std::to_string(a) + " is both deleted and a byproduct");
      if (a == r.init_a || a == r.init_b) {
        // an initiator may leave (the H of an H-abstraction); it then must not also stay in the network through a formed bond
        for (auto [x, y] : r.form)
          if ((x == a || y == a) && std::find(r.byproduct.begin(), r.byproduct.end(), x == a ? y : x) == r.byproduct.end())
            throw ReactError("reaction " + r.name + ": byproduct atom " + std::to_string(a) + " forms a bond to the network");
      }
    }
    if (r.weight <= 0) throw ReactError("reaction " + r.name + ": weight must be > 0");
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
  if (t.aromatic && !c.aromatic[i]) return false;
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

std::vector<Match> find_matches(const System& s, const ReactionTemplate& t, int reaction_index, const std::function<bool(uint32_t, uint32_t)>& allow,
                                double capture) {
  const double cap = capture > 0 ? capture : t.capture;
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
  Grid g(s, cap);
  const double cap2 = cap * cap;
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
        if (allow && !allow(a, b)) continue;
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

int count_sites(const System& s, const ReactionTemplate& t) { return int(site_groups(s, t).size()); }

std::vector<ChainSites> chain_sites(const System& s0, const std::vector<ReactionTemplate>& templates, const std::vector<int64_t>& chains) {
  System s = s0;
  s.has_mol = false;
  std::vector<int64_t> mol(s.atoms.size());
  if (chains.size() == s.atoms.size()) {
    for (size_t i = 0; i < mol.size(); ++i) mol[i] = chains[i] - 1;   // the chains a run started from (byproducts negative: left out below)
  } else {
    const auto m = s.molecules();
    for (size_t i = 0; i < mol.size(); ++i) mol[i] = m[i];
  }
  std::map<int64_t, ChainSites> by;
  for (size_t i = 0; i < s.atoms.size(); ++i) {
    if (mol[i] + 1 <= 0) continue;   // a byproduct of a run
    auto& c = by[mol[i] + 1];
    c.chain = mol[i] + 1;
    ++c.atoms;
    c.mass += s.mass_of(s.atoms[i]);
  }
  int largest = 0;
  for (const auto& [m, c] : by) largest = std::max(largest, c.atoms);
  std::map<int64_t, std::set<int64_t>> units;
  for (size_t i = 0; i < s.atoms.size(); ++i) if (s.atoms[i].resid > 0 && mol[i] + 1 > 0) units[mol[i] + 1].insert(s.atoms[i].resid);
  // each chain's sites: the groups of every template whose site atoms lie on it (a group counted once over templates)
  std::map<int64_t, std::set<std::vector<uint32_t>>> groups;
  for (const auto& t : templates)
    for (const auto& g : site_groups(s, t)) if (mol[g[0]] + 1 > 0) groups[mol[g[0]] + 1].insert(g);
  std::vector<ChainSites> r;
  for (auto& [m, c] : by) {
    if (c.atoms < 30 || c.atoms * 5 < largest) continue;
    c.sites = int(groups[m].size());
    c.units = int(units[m].size());
    r.push_back(c);
  }
  return r;
}

std::vector<std::vector<uint32_t>> site_groups(const System& s, const ReactionTemplate& t) {
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
  return {groups.begin(), groups.end()};
}

int apply_matches(System& s, const std::vector<ReactionTemplate>& templates, const std::vector<Match>& matches, bool keep_byproducts,
                  int* byproducts, std::vector<int64_t>* tag, std::vector<int64_t>* carry) {
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
    keep = keep && t.keep_charges && !(keep_byproducts && !t.byproduct.empty());   // a kept H2 / H2O takes no share: recomputed
    std::vector<int> leaving = t.remove;
    if (!keep_byproducts) leaving.insert(leaving.end(), t.byproduct.begin(), t.byproduct.end());
    for (int a : leaving) {
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
    if (keep_byproducts && !t.byproduct.empty()) {
      // the byproduct leaves the network as its own molecule: bonds from it to any other atom break, bonds among its atoms stay
      std::set<uint32_t> bp;
      for (int a : t.byproduct) bp.insert(at(a));
      for (auto it = bonds.begin(); it != bonds.end();)
        it = (bp.count(it->first) != bp.count(it->second)) ? bonds.erase(it) : std::next(it);
      if (tag && tag->size() == n) {   // a molecule of its own, belonging to no chain
        int64_t low = 0;
        for (int64_t v : *tag) low = std::min(low, v);
        const int64_t id = low - 1;   // negative: not a chain of the start
        for (uint32_t a : bp) (*tag)[a] = id;
      }
    }
    if (byproducts && !t.byproduct.empty()) ++*byproducts;
    for (int a : leaving) dead[at(a)] = 1;
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
  if (tag && tag->size() == n) {
    std::vector<int64_t> kept;
    for (size_t i = 0; i < n; ++i) if (!dead[i]) kept.push_back((*tag)[i]);
    *tag = std::move(kept);
  }
  if (carry && carry->size() == n) {
    std::vector<int64_t> kept;
    for (size_t i = 0; i < n; ++i) if (!dead[i]) kept.push_back((*carry)[i]);
    *carry = std::move(kept);
  }
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

namespace {

constexpr double kAvogadro = 6.02214076e23;

// Which chain each atom belongs to for the between-chains rule and the crosslink count. tag: the molecule each atom started
// in (negative: a byproduct molecule); polymer: those molecules that are chains. A chain atom belongs to its chain; an atom
// of a small molecule (a curative, a crosslinker) belongs to every chain it has bonded to so far, reached through small-
// molecule atoms only (extra: bonds chosen earlier in the same cycle). Empty: a small molecule bonded to no chain.
std::set<int64_t> owners(uint32_t a, const std::vector<std::vector<uint32_t>>& nb, const std::vector<int64_t>& tag, const std::set<int64_t>& polymer,
                         const std::multimap<uint32_t, uint32_t>& extra) {
  if (polymer.count(tag[a])) return {tag[a]};
  std::set<int64_t> own;
  std::vector<uint32_t> q{a};
  std::set<uint32_t> seen{a};
  auto visit = [&](uint32_t v) {
    if (polymer.count(tag[v])) { own.insert(tag[v]); return; }
    if (seen.insert(v).second) q.push_back(v);
  };
  for (size_t h = 0; h < q.size() && q.size() < 5000; ++h) {
    const uint32_t u = q[h];
    for (uint32_t v : nb[u]) visit(v);
    for (auto [it, e] = extra.equal_range(u); it != e; ++it) visit(it->second);
  }
  return own;
}

bool meet(const std::set<int64_t>& a, const std::set<int64_t>& b) {
  for (int64_t x : a) if (b.count(x)) return true;
  return false;
}

}  // namespace

void react(System& s, const ReactOptions& o, ReactReport* rep_out) {
  const auto t0 = std::chrono::steady_clock::now();
  ReactReport rep;
  if (o.templates.empty()) throw ReactError("no reaction templates");
  std::mt19937_64 rng(o.seed);
  caps::UniformReal<double> uni(0, 1);
  for (const auto& t : o.templates) rep.initial_sites = std::max(rep.initial_sites, count_sites(s, t));
  if (rep.initial_sites == 0) throw ReactError("no atom matches the counted site of any template; check the templates against the structure");
  if (s.cell.valid() && !s.unwrapped) make_molecules_whole(s);
  rep.field = o.field_name.empty() ? (o.retype ? "the assigned force field" : "built-in default (GAFF for C and H, UFF otherwise)") : o.field_name;

  // the chains of the start: molecules of at least 30 atoms and a fifth of the largest (curatives, crosslinkers, solvent
  // and fillers' small molecules are not chains)
  std::vector<int64_t> carry = o.carry.size() == s.atoms.size() ? o.carry : std::vector<int64_t>{};
  std::vector<int64_t> tag(s.atoms.size());
  std::set<int64_t> polymer;
  {
    if (o.chains.size() == s.atoms.size()) {
      tag = o.chains;   // the chains an earlier run started from
      rep.notes.push_back("chains as the earlier reaction run left them (the molecules it started from)");
    } else {
      System t = s;
      t.has_mol = false;
      const auto mol = t.molecules();
      for (size_t i = 0; i < s.atoms.size(); ++i) tag[i] = mol[i] + 1;
    }
    std::map<int64_t, int> count;
    std::map<int64_t, double> mass;
    for (size_t i = 0; i < s.atoms.size(); ++i) { ++count[tag[i]]; mass[tag[i]] += s.mass_of(s.atoms[i]); }
    int largest = 0;
    for (const auto& [m, k] : count) if (m > 0) largest = std::max(largest, k);
    for (const auto& [m, k] : count)
      if (m > 0 && k >= 30 && k * 5 >= largest) { polymer.insert(m); rep.chain_mass += mass[m]; }
    rep.chains = int(polymer.size());
    std::set<std::pair<int64_t, int64_t>> units;
    for (size_t i = 0; i < s.atoms.size(); ++i)
      if (polymer.count(tag[i]) && s.atoms[i].resid > 0) units.insert({tag[i], s.atoms[i].resid});
    rep.monomers = int(units.size());
    rep.volume = s.cell.valid() ? s.cell.volume() : 0;
  }
  const double avogadro_volume = rep.volume * 1e-30 * kAvogadro;   // links → mol/m³: divide by this
  // where each small molecule (a crosslinker) attached to a chain: its start id → (chain, repeat unit), for the link records
  std::map<int64_t, std::vector<std::pair<int64_t, int64_t>>> attached;
  std::map<int64_t, std::string> formula;   // each small molecule's formula at the start
  {
    std::map<int64_t, std::map<int, int>> el;
    for (size_t i = 0; i < s.atoms.size(); ++i) if (!polymer.count(tag[i])) el[tag[i]][s.atoms[i].element]++;
    for (const auto& [m, e] : el) {
      std::string f;
      for (int z : {6, 1}) if (e.count(z)) f += std::string(element(z).symbol) + (e.at(z) > 1 ? std::to_string(e.at(z)) : "");
      for (const auto& [z, k] : e) if (z != 6 && z != 1) f += std::string(element(z).symbol) + (k > 1 ? std::to_string(k) : "");
      formula[m] = f;
    }
  }
  if (o.target != ReactTarget::Conversion) {
    if (o.target_value <= 0) throw ReactError("the crosslink target must be > 0");
    if (rep.chains == 0) throw ReactError("a crosslink target needs chains (molecules of 30 atoms or more); this structure has none");
    double n = 0;
    switch (o.target) {
      case ReactTarget::Crosslinks: n = o.target_value; break;
      case ReactTarget::PerChain: n = o.target_value * rep.chains / 2.0; break;   // each link joins two chains
      case ReactTarget::Density:
        if (rep.volume <= 0) throw ReactError("a crosslink density target needs a periodic cell");
        n = o.target_value * avogadro_volume;
        break;
      case ReactTarget::Mc: n = rep.chain_mass / (2.0 * o.target_value); break;   // strands = 2 × links
      case ReactTarget::DegreePercent:
        if (rep.monomers == 0) throw ReactError("a degree-of-crosslinking target needs the chains' repeat units (residue numbers, as Grow writes them)");
        n = o.target_value / 100.0 * rep.monomers / 2.0;
        break;
      default: break;
    }
    rep.target_crosslinks = int(std::llround(n));
    if (rep.target_crosslinks < 1) {
      char b[200];
      std::snprintf(b, sizeof b, "the target is %.2f links in this cell (fewer than one): build a larger cell or ask for more", n);
      throw ReactError(b);
    }
  }

  // the force field for the structure as it is now
  auto field_now = [&](int cycle) -> std::shared_ptr<const ForceField> {
    if (!o.retype) return nullptr;
    try {
      return o.retype(s);
    } catch (const std::exception& e) {
      throw ReactError(rep.field + " cannot describe the structure after cycle " + std::to_string(cycle) + ": " + e.what() +
                       " — add the missing parameters (Force field › Fill gaps from a file), choose a force field that has them for the reaction, or react without retyping (topology only; type the product afterwards)");
    }
  };
  auto name_types = [&](const std::shared_ptr<const ForceField>& ff) {
    const ForceField def = ff ? ForceField{} : default_forcefield(s);
    const ForceField& f = ff ? *ff : def;
    for (size_t i = 0; i < s.atoms.size() && i < f.atom_type.size(); ++i) s.atoms[i].name = f.atom_type[i];
  };

  auto relax_now = [&](CycleRow& row, const std::shared_ptr<const ForceField>& ff) {
    if (!o.relax) return;
    RelaxOptions r;
    r.pushoff = false;
    r.ftol = o.relax_ftol;
    r.max_iterations = o.relax_iterations;
    r.energy = o.energy;
    r.field = ff;
    RelaxReport rr;
    try {
      relax(s, r, &rr);
    } catch (const FieldError& e) {
      throw ReactError(std::string("cannot relax after the reaction: ") + e.what() +
                       ". Run with relaxation off (topology only) or add parameters for these atoms.");
    }
    row.energy = rr.final.total();
    row.max_force = rr.fmax_final;
    name_types(ff);
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
  auto stabilise = [&](const std::set<int64_t>& ids, CycleRow& row, const std::shared_ptr<const ForceField>& ff) {
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
    r.field = ff;
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
    name_types(ff);
  };

  std::vector<char> linked(s.atoms.size(), 0);   // atoms of the bonds that joined two chains
  std::map<int64_t, int> site_used;                // each chain's sites that have reacted (sites_per_chain)
  if (o.sites_per_chain > 0) {
    char b[300];
    std::snprintf(b, sizeof b, "at most %d reactive sites per chain react", o.sites_per_chain);
    rep.notes.push_back(b);
    // a link between chains uses one site on each of its two chains: at most chains × sites / 2 links
    const int most = rep.chains * o.sites_per_chain / 2;
    if (rep.target_crosslinks > most) {
      std::snprintf(b, sizeof b, "the target (%d links) needs more sites than %d per chain allow: at most %d links (%d chains × %d sites ÷ 2) — the run stops there",
                    rep.target_crosslinks, o.sites_per_chain, most, rep.chains, o.sites_per_chain);
      rep.notes.push_back(b);
    }
  }
  // capture distances: the templates' own, raised together when auto capture finds no pair
  double extra = 0;
  auto capture_of = [&](size_t k) { return o.templates[k].capture + extra; };
  double top_capture = 0;
  for (const auto& t : o.templates) top_capture = std::max(top_capture, t.capture);
  if (o.auto_capture && o.capture_max < top_capture) throw ReactError("auto capture: the largest capture must be at least the templates' own");
  auto weight_of = [&](size_t k) { return k < o.weights.size() && o.weights[k] > 0 ? o.weights[k] : o.templates[k].weight; };

  int stall = 0;
  for (int cycle = 1; cycle <= o.max_cycles; ++cycle) {
    // the checkpoint: this cycle failing gives back the structure as the last one left it
    const System keep = o.keep_on_failure && !rep.cycles.empty() ? s : System{};
    const std::vector<int64_t> keep_tag = tag;
    const std::vector<int64_t> keep_carry = carry;
    const std::map<int64_t, int> keep_site_used = site_used;
    const std::vector<char> keep_linked = linked;
    const int reactions_before = rep.reactions, stall_before = stall, links_before = rep.crosslinks, intra_before = rep.intrachain,
              byproducts_before = rep.byproducts;
    const double extra_before = extra;
    int applied = 0;
    try {
    const auto nb = s.neighbours();
    const std::multimap<uint32_t, uint32_t> none;
    auto allow = [&](uint32_t a, uint32_t b) {
      return !o.between_chains || !meet(owners(a, nb, tag, polymer, none), owners(b, nb, tag, polymer, none));
    };
    std::vector<std::vector<Match>> per(o.templates.size());
    size_t found = 0;
    for (;;) {
      found = 0;
      for (size_t k = 0; k < o.templates.size(); ++k) {
        per[k] = find_matches(s, o.templates[k], int(k), allow, capture_of(k));
        found += per[k].size();
      }
      if (found || !o.auto_capture || top_capture + extra + o.capture_step > o.capture_max + 1e-9) break;
      extra += o.capture_step;   // no pair: look further
    }
    if (extra > extra_before) {
      char b[160];
      std::snprintf(b, sizeof b, "cycle %d: no pair within %.2f Å — auto capture raised to %.2f Å", cycle, top_capture + extra_before, top_capture + extra);
      rep.notes.push_back(b);
    }
    // choose: each pick re-checks the chain rule with the bonds already chosen this cycle and counts the links it makes
    std::vector<Match> chosen;
    std::set<uint32_t> busy;
    std::multimap<uint32_t, uint32_t> formed;
    int links = 0, intra = 0;
    std::vector<std::pair<uint32_t, uint32_t>> link_bonds;   // formed bonds of the links chosen this cycle
    std::vector<LinkRecord> pending;                          // this cycle's link records
    std::vector<std::pair<int64_t, std::pair<int64_t, int64_t>>> pending_attach;
    auto take = [&](const Match& m) {
      if (uni(rng) > o.templates[m.reaction].probability) return false;
      for (uint32_t a : m.atoms) if (busy.count(a)) return false;
      const ReactionTemplate& t = o.templates[m.reaction];
      auto at = [&](int map) { for (size_t q = 0; q < t.atoms.size(); ++q) if (t.atoms[q].map == map) return m.atoms[q]; return m.atoms[0]; };
      const auto oa = owners(at(t.init_a), nb, tag, polymer, formed), ob = owners(at(t.init_b), nb, tag, polymer, formed);
      const bool same = meet(oa, ob);
      if (o.between_chains && same) return false;
      const bool link = !oa.empty() && !ob.empty() && !same;
      if (link && rep.target_crosslinks > 0 && rep.crosslinks + links >= rep.target_crosslinks) return false;   // the target is met
      // each chain's own sites: at most sites_per_chain react (the template's site atoms on a chain)
      std::map<int64_t, int> uses;
      if (o.sites_per_chain > 0) {
        // every atom that forms a bond and stays in the network uses one of its chain's sites (both carbons of a C–C link)
        std::set<int> ends;
        for (auto [x, y] : t.form) ends.insert(x), ends.insert(y);
        for (int mp : t.byproduct) ends.erase(mp);
        for (int mp : t.remove) ends.erase(mp);
        for (int mp : ends) { const uint32_t a = at(mp); if (polymer.count(tag[a])) ++uses[tag[a]]; }
        for (const auto& [c, k] : uses) if (site_used[c] + k > o.sites_per_chain) return false;
      }
      for (uint32_t a : m.atoms) busy.insert(a);
      for (const auto& [c, k] : uses) site_used[c] += k;
      for (auto [x, y] : t.form) { formed.insert({at(x), at(y)}); formed.insert({at(y), at(x)}); }
      links += link;
      intra += same && !oa.empty();
      if (link) for (auto [x, y] : t.form) link_bonds.push_back({at(x), at(y)});
      {
        // the record: each side's chain and unit — a chain atom's own, or where the crosslinker it sits on attached
        const uint32_t ia = at(t.init_a), ib = at(t.init_b);
        auto side = [&](uint32_t i, const std::set<int64_t>& own, int64_t other_chain) -> std::tuple<int64_t, int64_t, int64_t> {
          if (polymer.count(tag[i])) return {tag[i], s.atoms[i].resid, 0};
          for (const auto& [c, u] : attached[tag[i]])
            if (own.count(c) && c != other_chain) return {c, u, tag[i]};
          return {own.empty() ? 0 : *own.begin(), 0, tag[i]};
        };
        if (link) {
          const int64_t cb = polymer.count(tag[ib]) ? tag[ib] : (ob.empty() ? 0 : *ob.begin());
          auto [ca, ua, va] = side(ia, oa, cb);
          auto [cb2, ub, vb] = side(ib, ob, ca);
          LinkRecord lr;
          lr.cycle = int(rep.cycles.size()) + 1;
          lr.reaction = t.name;
          lr.chain_a = ca, lr.unit_a = ua, lr.chain_b = cb2, lr.unit_b = ub;
          lr.via = va ? va : vb;
          if (lr.via) lr.via_name = formula[lr.via];
          pending.push_back(lr);
        } else if (polymer.count(tag[ia]) != polymer.count(tag[ib])) {
          // a crosslinker meets a chain: remember where (a later link through it names this unit)
          const uint32_t chain_atom = polymer.count(tag[ia]) ? ia : ib, small = chain_atom == ia ? ib : ia;
          pending_attach.push_back({tag[small], {tag[chain_atom], s.atoms[chain_atom].resid}});
        }
      }
      chosen.push_back(m);
      return true;
    };
    if (o.selection == 1 && o.templates.size() > 1) {
      std::vector<size_t> next(o.templates.size(), 0);
      while (int(chosen.size()) < o.max_per_cycle) {
        double total = 0;
        for (size_t k = 0; k < per.size(); ++k) if (next[k] < per[k].size()) total += weight_of(k);
        if (total <= 0) break;
        double r = uni(rng) * total;
        size_t k = 0;
        for (; k < per.size(); ++k) {
          if (next[k] >= per[k].size()) continue;
          r -= weight_of(k);
          if (r <= 0) break;
        }
        if (k == per.size()) for (k = per.size(); k-- > 0;) if (next[k] < per[k].size()) break;
        while (next[k] < per[k].size() && !take(per[k][next[k]])) ++next[k];
        if (next[k] < per[k].size()) ++next[k];
      }
    } else {
      std::vector<Match> all;
      for (auto& v : per) all.insert(all.end(), v.begin(), v.end());
      std::stable_sort(all.begin(), all.end(), [](const Match& a, const Match& b) { return a.distance < b.distance; });
      for (const auto& m : all) {
        if (int(chosen.size()) >= o.max_per_cycle) break;
        take(m);
      }
    }
    std::set<int64_t> site_ids;
    for (const auto& m : chosen)
      for (uint32_t a : m.atoms) site_ids.insert(s.atoms[a].id);
    if (!chosen.empty()) {
      // the linked atoms follow the compaction: deleted atoms (and byproducts, unless kept) leave, the rest keep their order
      std::vector<char> gone(linked.size(), 0);
      for (const auto& m : chosen) {
        const ReactionTemplate& t = o.templates[m.reaction];
        auto at = [&](int map) { for (size_t q = 0; q < t.atoms.size(); ++q) if (t.atoms[q].map == map) return m.atoms[q]; return m.atoms[0]; };
        for (int a : t.remove) gone[at(a)] = 1;
        if (!o.keep_byproducts) for (int a : t.byproduct) gone[at(a)] = 1;
      }
      for (auto [x, y] : link_bonds) linked[x] = linked[y] = 1;
      std::vector<char> next;
      for (size_t i = 0; i < linked.size(); ++i) if (!gone[i]) next.push_back(linked[i]);
      linked = std::move(next);
    }
    applied = chosen.empty() ? 0 : apply_matches(s, o.templates, chosen, o.keep_byproducts, &rep.byproducts, &tag, &carry);
    // links join chains across the cell walls: every molecule (now the network) whole again, so each cycle's frame and the
    // live view show the joined chains side by side and bonded atoms carry consistent images
    if (applied > 0 && s.cell.valid()) make_molecules_whole(s);
    if (applied) {
      for (auto& lr : pending) rep.links.push_back(lr);
      for (const auto& [m, cu] : pending_attach) attached[m].push_back(cu);
    }
    CycleRow row;
    row.cycle = cycle;
    row.reactions = applied;
    row.capture = top_capture + extra;
    rep.reactions += applied;
    rep.crosslinks += applied ? links : 0;
    rep.intrachain += applied ? intra : 0;
    row.total = rep.reactions;
    row.crosslinks = rep.crosslinks;
    row.target = rep.target_crosslinks;
    if (avogadro_volume > 0) row.density = rep.crosslinks / avogadro_volume;
    if (rep.monomers > 0) row.degree = 200.0 * rep.crosslinks / rep.monomers;
    row.conversion = double(rep.reactions) / rep.initial_sites;
    std::shared_ptr<const ForceField> ff;
    if (applied > 0) {
      stall = 0;
      ff = field_now(cycle);
      if (o.during_md) stabilise(site_ids, row, ff);
      else relax_now(row, ff);
    } else {
      ++stall;
    }
    if (o.md_ps > 0 && (o.during_md || applied > 0 || stall <= o.stall_cycles)) {
      if (!ff) ff = field_now(cycle);
      DynamicsOptions d;
      d.steps = std::max<int64_t>(1, std::llround(o.md_ps * 1000));
      d.temperature = o.temperature;
      d.new_velocities = s.velocities.size() != s.atoms.size();
      d.seed = o.seed + uint64_t(cycle);
      d.thermo_every = int(d.steps);
      d.frame_every = 0;
      d.energy = o.energy;
      d.field = ff;
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
    if (o.live) o.live(s, row, tag, linked);
    } catch (const std::exception& e) {
      if (!rep.cycles.empty() && rep.cycles.back().cycle == cycle) rep.cycles.pop_back();
      if (!o.keep_on_failure || rep.cycles.empty()) throw;
      s = keep;
      tag = keep_tag;
      carry = keep_carry;
      site_used = keep_site_used;
      linked = keep_linked;
      rep.reactions = reactions_before;
      rep.crosslinks = links_before;
      rep.intrachain = intra_before;
      rep.byproducts = byproducts_before;
      rep.links.resize(size_t(links_before));
      extra = extra_before;
      stall = stall_before;
      rep.failed_cycle = cycle;
      rep.failure = e.what();
      rep.notes.push_back("failed at cycle " + std::to_string(cycle) + ": " + rep.failure + " — the structure after cycle " + std::to_string(cycle - 1) +
                          " is kept (restart from it)");
      break;
    }
    const CycleRow& row = rep.cycles.back();
    if (o.progress && !o.progress(row)) throw ReactError("reaction run cancelled");
    if (rep.target_crosslinks > 0 && rep.crosslinks >= rep.target_crosslinks) {
      rep.notes.push_back("crosslink target reached");
      break;
    }
    if (o.target == ReactTarget::Conversion && row.conversion >= o.target_conversion) {
      rep.notes.push_back("target conversion reached");
      break;
    }
    if (applied == 0 && !o.during_md && (o.md_ps <= 0 || stall > o.stall_cycles)) {
      rep.notes.push_back(o.md_ps > 0 ? "no reactive pairs within the capture distance after dynamics; stopped"
                                      : std::string("no more reactive pairs within the capture distance; stopped (") +
                                            (o.auto_capture ? "auto capture reached its largest distance; " : "") +
                                            "enable dynamics between cycles to let groups diffuse)");
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
  // the network
  if (avogadro_volume > 0) rep.density = rep.crosslinks / avogadro_volume;
  if (rep.chains > 0) rep.per_chain = 2.0 * rep.crosslinks / rep.chains;
  if (rep.crosslinks > 0) rep.mc = rep.chain_mass / (2.0 * rep.crosslinks);
  if (rep.monomers > 0) rep.degree = 200.0 * rep.crosslinks / rep.monomers;
  rep.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  char b[320];
  const auto& last = rep.cycles.empty() ? CycleRow{} : rep.cycles.back();
  std::snprintf(b, sizeof b, "%d reactions in %zu cycles (%d counted sites) · conversion %.3f · largest cluster %.1f%% of the mass · %.1f s",
                rep.reactions, rep.cycles.size(), rep.initial_sites, last.conversion, 100 * last.clusters.largest_fraction, rep.seconds);
  rep.notes.insert(rep.notes.begin(), b);
  if (rep.chains > 0) {
    std::snprintf(b, sizeof b, "%d links between chains (%d chains; %.2f per chain) · ν = %.4g mol/m³ · Mc ≈ %.4g g/mol", rep.crosslinks, rep.chains,
                  rep.per_chain, rep.density, rep.mc);
    rep.notes.insert(rep.notes.begin() + 1, b);
    if (rep.monomers > 0) {
      std::snprintf(b, sizeof b, "degree of crosslinking DC = 2 × links / monomers = %.2f %% (%d repeat units)", rep.degree, rep.monomers);
      rep.notes.insert(rep.notes.begin() + 2, b);
    }
    if (rep.intrachain) rep.notes.push_back(std::to_string(rep.intrachain) + " links closed within one chain (loops; not counted as crosslinks)");
  }
  if (rep.target_crosslinks > 0) {
    std::snprintf(b, sizeof b, "target %d links, achieved %d (%.0f%%)", rep.target_crosslinks, rep.crosslinks, 100.0 * rep.crosslinks / rep.target_crosslinks);
    rep.notes.insert(rep.notes.begin() + 1, b);
  }
  for (const auto& lr : rep.links) {
    char b[240];
    std::snprintf(b, sizeof b, "link %zu (cycle %d, %s): chain %lld unit %lld — %s — chain %lld unit %lld", size_t(&lr - rep.links.data()) + 1, lr.cycle,
                  lr.reaction.c_str(), static_cast<long long>(lr.chain_a), static_cast<long long>(lr.unit_a),
                  lr.via ? (lr.via_name + " #" + std::to_string(lr.via)).c_str() : "direct bond", static_cast<long long>(lr.chain_b), static_cast<long long>(lr.unit_b));
    rep.notes.push_back(b);
  }
  if (rep.byproducts)
    rep.notes.push_back(std::to_string(rep.byproducts) + (o.keep_byproducts ? " byproduct molecules kept in the cell" : " byproduct molecules removed"));
  rep.notes.push_back("force field during the run: " + rep.field);
  rep.chains_after = tag;
  rep.carry_after = carry;
  // new bonds join molecules across the cell: every molecule whole again, so bonded atoms carry consistent image flags
  if (s.cell.valid()) make_molecules_whole(s);
  s.unwrapped = true;
  if (rep_out) *rep_out = std::move(rep);
}

}  // namespace caps
