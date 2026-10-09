// CAPS chemistry-aware coarse-grained mapping (see caps/cg_rules.hpp).
#include "caps/cg_rules.hpp"

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <cmath>
#include <cstdio>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>

#include "caps/elements.hpp"
#include "caps/resolution.hpp"
#include "caps/typing.hpp"

namespace caps {

namespace {

uint64_t mix(uint64_t h, uint64_t v) {   // FNV-1a over the bytes of v, then a final avalanche
  for (int k = 0; k < 8; ++k) { h ^= (v >> (8 * k)) & 0xff; h *= 1099511628211ULL; }
  return h;
}
uint64_t avalanche(uint64_t x) {
  x ^= x >> 33; x *= 0xff51afd7ed558ccdULL; x ^= x >> 33; x *= 0xc4ceb9fe1a85ec53ULL; x ^= x >> 33;
  return x;
}

std::string hex(uint64_t v) {
  char b[17];
  std::snprintf(b, sizeof b, "%016llx", static_cast<unsigned long long>(v));
  return b;
}

// Hill formula: C, H, then the rest alphabetically (alphabetical throughout without carbon)
std::string hill(const std::map<std::string, int>& counts) {
  std::string out;
  auto put = [&](const std::string& s, int n) { out += s; if (n > 1) out += std::to_string(n); };
  const bool carbon = counts.count("C") > 0;
  if (carbon) {
    put("C", counts.at("C"));
    if (counts.count("H")) put("H", counts.at("H"));
  }
  for (const auto& [s, n] : counts)
    if (!carbon || (s != "C" && s != "H")) put(s, n);
  return out;
}

}  // namespace

std::string cg_key(const std::vector<std::string>& seq) {
  std::vector<std::string> r(seq.rbegin(), seq.rend());
  const auto& use = r < seq ? r : seq;
  std::string out;
  for (size_t k = 0; k < use.size(); ++k) out += (k ? "-" : "") + use[k];
  return out;
}

CgTypes cg_types_from_json(const Json& j) {
  CgTypes t;
  auto list = [&](const char* k, std::vector<std::string>& v) {
    if (!j.has(k)) return;
    for (const auto& x : j[k].items()) v.push_back(x.str());
  };
  list("beads", t.beads);
  list("bonds", t.bonds);
  list("angles", t.angles);
  list("dihedrals", t.dihedrals);
  // keys as given or reversed: written canonically
  for (auto* v : {&t.bonds, &t.angles, &t.dihedrals})
    for (auto& s : *v) {
      std::vector<std::string> seq;
      size_t a = 0;
      for (size_t b; (b = s.find('-', a)) != std::string::npos; a = b + 1) seq.push_back(s.substr(a, b - a));
      seq.push_back(s.substr(a));
      s = cg_key(seq);
    }
  return t;
}

Json cg_types_json(const CgTypes& t) {
  Json j = Json::object();
  auto arr = [](const std::vector<std::string>& v) {
    Json a = Json::array();
    for (const auto& s : v) a.push_back(s);
    return a;
  };
  j["beads"] = arr(t.beads);
  j["bonds"] = arr(t.bonds);
  j["angles"] = arr(t.angles);
  j["dihedrals"] = arr(t.dihedrals);
  return j;
}

CgTypes cg_types_union(const std::vector<CgTypes>& lists) {
  std::set<std::string> b, bo, an, di;
  for (const auto& t : lists) {
    b.insert(t.beads.begin(), t.beads.end());
    bo.insert(t.bonds.begin(), t.bonds.end());
    an.insert(t.angles.begin(), t.angles.end());
    di.insert(t.dihedrals.begin(), t.dihedrals.end());
  }
  return {{b.begin(), b.end()}, {bo.begin(), bo.end()}, {an.begin(), an.end()}, {di.begin(), di.end()}};
}

std::vector<int> CgMapping::bead_of() const {
  std::vector<int> out(size_t(atoms), -1);
  for (size_t b = 0; b < bead_atoms.size(); ++b)
    for (uint32_t a : bead_atoms[b]) out[a] = int(b);
  return out;
}

namespace {

// Bonds, angles, dihedrals and chain order from the beads (bead_atoms, bead_mol set).
void finish_topology(CgMapping& m, const System& aa) {
  const auto of = m.bead_of();
  std::set<std::pair<int, int>> bb;
  for (const auto& b : aa.bonds) {
    int x = of[b.i], y = of[b.j];
    if (x == y || x < 0 || y < 0) continue;
    if (x > y) std::swap(x, y);
    bb.insert({x, y});
  }
  m.bonds.assign(bb.begin(), bb.end());
  const size_t nb = m.beads();
  std::vector<std::vector<int>> nbr(nb);
  for (const auto& [x, y] : m.bonds) nbr[size_t(x)].push_back(y), nbr[size_t(y)].push_back(x);
  for (auto& v : nbr) std::sort(v.begin(), v.end());
  m.angles.clear();
  m.dihedrals.clear();
  for (size_t j = 0; j < nb; ++j)
    for (size_t a = 0; a < nbr[j].size(); ++a)
      for (size_t c = a + 1; c < nbr[j].size(); ++c) m.angles.push_back({nbr[j][a], int(j), nbr[j][c]});
  for (const auto& [j, k] : m.bonds)
    for (int i : nbr[size_t(j)])
      if (i != k)
        for (int l : nbr[size_t(k)])
          if (l != j && l != i) m.dihedrals.push_back({i, j, k, l});
  // chains: per molecule, end to end
  std::map<int, std::vector<int>> by_mol;
  for (size_t b = 0; b < nb; ++b) by_mol[m.bead_mol[b]].push_back(int(b));
  auto first_atom = [&](int b) { return *std::min_element(m.bead_atoms[size_t(b)].begin(), m.bead_atoms[size_t(b)].end()); };
  m.chains.clear();
  m.branched = false;
  for (auto& [mol, beads] : by_mol) {
    int start = -1;
    bool branched = false;
    for (int b : beads) {
      if (nbr[size_t(b)].size() > 2) branched = true;
      if (nbr[size_t(b)].size() <= 1 && (start < 0 || first_atom(b) < first_atom(start))) start = b;
    }
    if (start < 0) {   // a ring (or one bead): the lowest atom first
      start = beads.front();
      for (int b : beads) if (first_atom(b) < first_atom(start)) start = b;
    }
    m.branched = m.branched || branched;
    std::vector<int> order{start};
    std::set<int> seen{start};
    for (size_t h = 0; h < order.size(); ++h)   // breadth first: a linear chain comes out end to end
      for (int w : nbr[size_t(order[h])])
        if (seen.insert(w).second) order.push_back(w);
    if (order.size() != beads.size())   // a molecule id holding several bonded pieces: each piece after the first
      for (int b : beads)
        if (!seen.count(b)) {
          order.push_back(b);
          seen.insert(b);
          for (size_t h = order.size() - 1; h < order.size(); ++h)
            for (int w : nbr[size_t(order[h])])
              if (seen.insert(w).second) order.push_back(w);
        }
    m.chains.push_back(std::move(order));
  }
  if (m.branched) m.notes.push_back("branched molecules: their beads are listed breadth first from an end");
}

// Formula, aromaticity and a canonical key of each fragment; fragments of one class share the key.
void classify(CgMapping& m, const System& aa, const Perception& p, const std::vector<int>& of) {
  const size_t nb = m.beads();
  m.bead_class.assign(nb, "");
  for (size_t b = 0; b < nb; ++b) {
    const auto& at = m.bead_atoms[b];
    std::map<std::string, int> counts;
    bool arom = false;
    std::vector<uint32_t> heavy;
    for (uint32_t a : at) {
      counts[element(aa.atoms[a].element).symbol]++;
      if (aa.atoms[a].element != 1) heavy.push_back(a);
      arom = arom || p.aromatic[a];
    }
    // Weisfeiler–Lehman over the heavy atoms: element, aromatic, hydrogens, bonds leaving the bead
    std::map<uint32_t, uint64_t> lab;
    for (uint32_t a : heavy) {
      int out = 0;
      for (uint32_t w : p.nb[a]) out += of[w] != int(b) && aa.atoms[w].element != 1;
      lab[a] = mix(mix(mix(mix(1469598103934665603ULL, uint64_t(aa.atoms[a].element)), p.aromatic[a] ? 1 : 0), uint64_t(p.hcount[a])), uint64_t(out));
    }
    for (int round = 0; round < 4; ++round) {
      std::map<uint32_t, uint64_t> next;
      for (uint32_t a : heavy) {
        std::vector<uint64_t> ns;
        for (size_t k = 0; k < p.nb[a].size(); ++k) {
          const uint32_t w = p.nb[a][k];
          if (lab.count(w)) ns.push_back(mix(lab[w], uint64_t(p.order[a][k])));
        }
        std::sort(ns.begin(), ns.end());
        uint64_t h = mix(1469598103934665603ULL, lab[a]);
        for (uint64_t v : ns) h = mix(h, v);
        next[a] = avalanche(h);
      }
      lab = std::move(next);
    }
    std::vector<uint64_t> all;
    for (const auto& [a, l] : lab) all.push_back(l);
    std::sort(all.begin(), all.end());
    uint64_t key = mix(1469598103934665603ULL, uint64_t(all.size()));
    for (uint64_t v : all) key = mix(key, v);
    const std::string k = hex(avalanche(key));
    m.bead_class[b] = k;
    auto it = std::find_if(m.classes.begin(), m.classes.end(), [&](const CgFragmentClass& c) { return c.key == k; });
    if (it == m.classes.end()) m.classes.push_back({k, hill(counts), arom, "", 1});
    else ++it->count;
  }
}

// The first naming rule that embeds in one fragment of each class names the class.
void name_classes(CgMapping& m, const System& aa, const Perception& p, const CgRules& r) {
  std::vector<Smarts> sms;
  for (const auto& nr : r.names) sms.emplace_back(nr.smarts);
  std::vector<char> outside(aa.atoms.size(), 1);
  int auto_n = 0;
  std::vector<std::string> unnamed;
  for (auto& c : m.classes) {
    if (!c.name.empty()) continue;
    size_t rep = 0;
    while (m.bead_class[rep] != c.key) ++rep;
    const auto& at = m.bead_atoms[rep];
    for (uint32_t a : at) outside[a] = 0;
    for (size_t k = 0; k < sms.size() && c.name.empty(); ++k)
      for (uint32_t a : at) {
        if (!sms[k].embeddings(aa, p, a, &outside).empty()) { c.name = r.names[k].name; break; }
      }
    for (uint32_t a : at) outside[a] = 1;
    if (c.name.empty()) {
      if (r.strict && !r.names.empty()) {
        unnamed.push_back(c.formula + (c.aromatic ? " (aromatic)" : "") + " × " + std::to_string(c.count) + " [class " + c.key.substr(0, 8) + "]");
        continue;
      }
      c.name = "F" + std::to_string(++auto_n);
      m.notes.push_back("fragment " + c.formula + " (class " + c.key.substr(0, 8) + ", " + std::to_string(c.count) + " beads) named " + c.name + " by its class");
    }
  }
  if (!unnamed.empty()) {
    std::string msg = std::to_string(unnamed.size()) + " fragment class" + (unnamed.size() > 1 ? "es" : "") + " no naming rule matches: ";
    for (size_t k = 0; k < unnamed.size(); ++k) msg += (k ? "; " : "") + unnamed[k];
    throw std::invalid_argument(msg + " (add a rule for each, or set strict off to name them by class)");
  }
  m.bead_kind.assign(m.beads(), "");
  for (size_t b = 0; b < m.beads(); ++b)
    for (const auto& c : m.classes)
      if (c.key == m.bead_class[b]) m.bead_kind[b] = c.name;
}

void anchors(CgMapping& m, const System& aa, const Perception& p) {
  m.bead_anchor.assign(m.beads(), -1);
  if (m.position.rfind("atom:", 0) != 0) return;
  const Smarts sm(m.position.substr(5));
  std::vector<char> outside(aa.atoms.size(), 1);
  int missing = 0;
  for (size_t b = 0; b < m.beads(); ++b) {
    for (uint32_t a : m.bead_atoms[b]) outside[a] = 0;
    for (uint32_t a : m.bead_atoms[b]) {
      const auto e = sm.embeddings(aa, p, a, &outside);
      if (!e.empty()) { m.bead_anchor[b] = int(e.front().front()); break; }
    }
    for (uint32_t a : m.bead_atoms[b]) outside[a] = 1;
    if (m.bead_anchor[b] < 0) ++missing;
  }
  if (missing) m.notes.push_back(std::to_string(missing) + " beads have no atom matching " + m.position.substr(5) + ": they sit at their centre of mass");
}

}  // namespace

CgMapping cg_mapping(const System& aa, const CgRules& r) {
  const size_t n = aa.atoms.size();
  if (n == 0) throw std::invalid_argument("no atoms to map");
  const int forms = int(!r.cut.empty()) + int(!r.beads.empty()) + int(!r.explicit_bead.empty());
  if (forms != 1) throw std::invalid_argument("give one mapping form: cut (bond SMARTS), beads (fragment SMARTS) or explicit (atom → bead)");
  for (const auto& a : aa.atoms)
    if (a.element <= 0) throw std::invalid_argument("atoms without an element (the masses of the data file name none): rules need elements");
  if (r.position != "com" && r.position != "cog" && r.position.rfind("atom:", 0) != 0)
    throw std::invalid_argument("bead position: com, cog or atom:<SMARTS>");
  const Perception p = perceive(aa);
  int nmol = 0;
  const auto mol = aa.molecules(&nmol);
  CgMapping m;
  m.atoms = int(n);
  m.position = r.position;
  std::vector<int> of(n, -1);
  std::vector<std::string> cover_kind;   // the exact-cover form: each bead's fragment name
  if (!r.cut.empty()) {
    std::set<std::pair<uint32_t, uint32_t>> cut;
    for (const auto& text : r.cut) {
      const Smarts sm(text);
      if (sm.size() < 2) throw std::invalid_argument("a cut pattern needs at least two atoms (the bond between its first and last): " + text);
      size_t found = 0;
      for (uint32_t a = 0; a < n; ++a)
        for (const auto& e : sm.embeddings(aa, p, a)) {
          const auto& nb = p.nb[e.front()];
          if (std::find(nb.begin(), nb.end(), e.back()) == nb.end()) continue;   // the first and last pattern atoms must be bonded
          if (cut.insert(std::minmax(e.front(), e.back())).second) ++found;
        }
      m.notes.push_back("cut " + text + ": " + std::to_string(found) + " bonds");
    }
    for (uint32_t s = 0; s < n; ++s) {
      if (of[s] >= 0) continue;
      const int b = int(m.bead_atoms.size());
      std::vector<uint32_t> q{s};
      of[s] = b;
      for (size_t h = 0; h < q.size(); ++h)
        for (uint32_t w : p.nb[q[h]])
          if (of[w] < 0 && !cut.count(std::minmax(q[h], w))) of[w] = b, q.push_back(w);
      std::sort(q.begin(), q.end());
      m.bead_atoms.push_back(std::move(q));
    }
  } else if (!r.beads.empty()) {
    std::vector<BeadRule> rules;
    for (const auto& b : r.beads) rules.push_back({{b.name}, b.smarts, ""});
    BeadMapReport rep;
    const System cg = map_to_beads(aa, rules, {}, {}, &rep);
    if (!rep.uncovered.empty())
      throw std::invalid_argument(std::to_string(rep.uncovered.size()) + " atoms are covered by no bead fragment (first: atom " + std::to_string(rep.uncovered.front() + 1) +
                                  ", " + element(aa.atoms[rep.uncovered.front()].element).symbol + "): add a fragment for them");
    m.bead_atoms.assign(cg.atoms.size(), {});
    for (uint32_t a = 0; a < n; ++a) of[a] = rep.site_of[a], m.bead_atoms[size_t(of[a])].push_back(a);
    for (const auto& b : cg.atoms) cover_kind.push_back(b.name);
  } else {
    if (r.explicit_bead.size() != n) throw std::invalid_argument("the explicit mapping has " + std::to_string(r.explicit_bead.size()) + " entries for " + std::to_string(n) + " atoms");
    int nb = 0;
    for (int b : r.explicit_bead) {
      if (b < 0) throw std::invalid_argument("every atom needs a bead in the explicit mapping");
      nb = std::max(nb, b + 1);
    }
    m.bead_atoms.assign(size_t(nb), {});
    for (uint32_t a = 0; a < n; ++a) of[a] = r.explicit_bead[a], m.bead_atoms[size_t(of[a])].push_back(a);
    for (size_t b = 0; b < m.bead_atoms.size(); ++b)
      if (m.bead_atoms[b].empty()) throw std::invalid_argument("bead " + std::to_string(b) + " of the explicit mapping has no atoms");
  }
  const size_t nb = m.bead_atoms.size();
  m.bead_mol.assign(nb, 0);
  m.bead_mass.assign(nb, 0.0);
  for (size_t b = 0; b < nb; ++b) {
    m.bead_mol[b] = mol[m.bead_atoms[b].front()];
    for (uint32_t a : m.bead_atoms[b]) {
      if (mol[a] != m.bead_mol[b]) throw std::invalid_argument("bead " + std::to_string(b + 1) + " spans two molecules");
      m.bead_mass[b] += aa.mass_of(aa.atoms[a]);
    }
    m.mass_cg += m.bead_mass[b];
  }
  for (const auto& a : aa.atoms) m.mass_aa += aa.mass_of(a);
  classify(m, aa, p, of);
  if (!r.beads.empty()) {   // the exact-cover form names its beads itself
    m.bead_kind = cover_kind;
    for (auto& c : m.classes)
      for (size_t b = 0; b < nb; ++b)
        if (m.bead_class[b] == c.key) { c.name = cover_kind[b]; break; }
  } else {
    name_classes(m, aa, p, r);
  }
  anchors(m, aa, p);
  finish_topology(m, aa);
  return m;
}

CgRules cg_rules_from_json(const Json& j) {
  CgRules r;
  auto named = [](const Json& x) {
    std::vector<CgNameRule> out;
    if (x.is_object())
      for (const auto& [k, v] : x.members()) out.push_back({k, v.str()});
    else
      for (const auto& e : x.items()) out.push_back({e["name"].str(), e["smarts"].str()});
    return out;
  };
  if (j.has("cut"))
    for (const auto& c : j["cut"].items()) r.cut.push_back(c.str());
  if (j.has("beads")) r.beads = named(j["beads"]);
  if (j.has("names")) r.names = named(j["names"]);
  if (j.has("explicit"))
    for (const auto& c : j["explicit"].items()) r.explicit_bead.push_back(int(c.number()));
  r.position = j.text("position", "com");
  if (j.has("strict")) r.strict = j["strict"].boolean();
  return r;
}

Json cg_mapping_json(const CgMapping& m, const System& aa) {
  Json j = Json::object();
  j["format"] = "caps-cg-map";
  j["version"] = 1;
  j["position"] = m.position;
  j["atoms"] = m.atoms;
  auto id = [&](uint32_t a) { return aa.atoms[a].id > 0 ? double(aa.atoms[a].id) : double(a + 1); };
  Json beads = Json::array();
  for (size_t b = 0; b < m.beads(); ++b) {
    Json e = Json::object();
    e["kind"] = m.bead_kind[b];
    e["class"] = m.bead_class[b];
    e["mol"] = m.bead_mol[b];
    e["mass"] = m.bead_mass[b];
    Json at = Json::array();
    for (uint32_t a : m.bead_atoms[b]) at.push_back(id(a));
    e["atoms"] = at;
    if (!m.bead_anchor.empty() && m.bead_anchor[b] >= 0) e["anchor"] = id(uint32_t(m.bead_anchor[b]));
    beads.push_back(e);
  }
  j["beads"] = beads;
  Json chains = Json::array();
  for (const auto& c : m.chains) {
    Json x = Json::array();
    for (int b : c) x.push_back(b);
    chains.push_back(x);
  }
  j["chains"] = chains;
  Json bonds = Json::array();
  for (const auto& [x, y] : m.bonds) {
    Json e = Json::array();
    e.push_back(x);
    e.push_back(y);
    bonds.push_back(e);
  }
  j["bonds"] = bonds;
  Json cls = Json::array();
  for (const auto& c : m.classes) {
    Json e = Json::object();
    e["key"] = c.key;
    e["formula"] = c.formula;
    e["aromatic"] = c.aromatic;
    e["name"] = c.name;
    e["count"] = c.count;
    cls.push_back(e);
  }
  j["classes"] = cls;
  return j;
}

CgMapping cg_mapping_from_json(const Json& j, const System& aa) {
  if (j.text("format") != "caps-cg-map") throw std::invalid_argument("not a CAPS mapping (format caps-cg-map)");
  CgMapping m;
  m.atoms = int(aa.atoms.size());
  if (int(j.num("atoms", 0)) != m.atoms)
    throw std::invalid_argument("the mapping is for " + std::to_string(int(j.num("atoms", 0))) + " atoms, the structure has " + std::to_string(m.atoms));
  m.position = j.text("position", "com");
  std::map<int64_t, uint32_t> index;
  for (uint32_t a = 0; a < aa.atoms.size(); ++a) index[aa.atoms[a].id > 0 ? aa.atoms[a].id : int64_t(a + 1)] = a;
  auto at = [&](double id) {
    auto it = index.find(int64_t(id));
    if (it == index.end()) throw std::invalid_argument("the mapping names atom " + std::to_string(int64_t(id)) + ", which the structure does not have");
    return it->second;
  };
  for (const auto& e : j["beads"].items()) {
    std::vector<uint32_t> v;
    for (const auto& x : e["atoms"].items()) v.push_back(at(x.number()));
    double mass = 0;
    for (uint32_t a : v) mass += aa.mass_of(aa.atoms[a]);
    m.bead_atoms.push_back(v);
    m.bead_kind.push_back(e["kind"].str());
    m.bead_class.push_back(e.text("class"));
    m.bead_mol.push_back(int(e.num("mol", 0)));
    m.bead_mass.push_back(mass);
    m.bead_anchor.push_back(e.has("anchor") ? int(at(e["anchor"].number())) : -1);
    m.mass_cg += mass;
  }
  for (const auto& a : aa.atoms) m.mass_aa += aa.mass_of(a);
  std::vector<int> seen(aa.atoms.size(), 0);
  for (const auto& b : m.bead_atoms)
    for (uint32_t a : b) ++seen[a];
  for (size_t a = 0; a < seen.size(); ++a)
    if (seen[a] != 1) throw std::invalid_argument("atom " + std::to_string(a + 1) + " is in " + std::to_string(seen[a]) + " beads of the mapping (each atom belongs to one)");
  if (j.has("classes"))
    for (const auto& c : j["classes"].items()) m.classes.push_back({c.text("key"), c.text("formula"), c.has("aromatic") && c["aromatic"].boolean(), c.text("name"), int(c.num("count", 0))});
  finish_topology(m, aa);
  return m;
}

std::vector<Vec3> cg_positions(const CgMapping& m, const System& aa, const std::vector<Vec3>& pos, const Cell& cell) {
  if (pos.size() != size_t(m.atoms)) throw std::invalid_argument("the frame has " + std::to_string(pos.size()) + " atoms, the mapping " + std::to_string(m.atoms));
  const bool pbc = cell.valid();
  const bool com = m.position == "com";
  std::vector<Vec3> out(m.beads());
  for (size_t b = 0; b < m.beads(); ++b) {
    if (!m.bead_anchor.empty() && m.bead_anchor[b] >= 0) { out[b] = pos[size_t(m.bead_anchor[b])]; continue; }
    const auto& at = m.bead_atoms[b];
    const Vec3 ref = pos[at.front()];
    Vec3 c{0, 0, 0};
    double w = 0;
    for (uint32_t a : at) {
      const double wa = com ? aa.mass_of(aa.atoms[a]) : 1.0;
      const Vec3 p = pbc ? ref + cell.minimum_image(pos[a] - ref) : pos[a];
      c = c + p * wa;
      w += wa;
    }
    out[b] = c * (1.0 / w);
  }
  return out;
}

CgTypes cg_types_of(const CgMapping& m) {
  std::set<std::string> b, bo, an, di;
  for (const auto& k : m.bead_kind) b.insert(k);
  for (const auto& [x, y] : m.bonds) bo.insert(cg_key({m.bead_kind[size_t(x)], m.bead_kind[size_t(y)]}));
  for (const auto& t : m.angles) an.insert(cg_key({m.bead_kind[size_t(t[0])], m.bead_kind[size_t(t[1])], m.bead_kind[size_t(t[2])]}));
  for (const auto& t : m.dihedrals) di.insert(cg_key({m.bead_kind[size_t(t[0])], m.bead_kind[size_t(t[1])], m.bead_kind[size_t(t[2])], m.bead_kind[size_t(t[3])]}));
  return {{b.begin(), b.end()}, {bo.begin(), bo.end()}, {an.begin(), an.end()}, {di.begin(), di.end()}};
}

System cg_structure(const CgMapping& m, const System& aa, const std::vector<Vec3>& beads, const Cell& cell, const CgTypes& types_in) {
  const CgTypes own = cg_types_of(m);
  const CgTypes& t = types_in.empty() ? own : types_in;
  auto check = [](const std::vector<std::string>& have, const std::vector<std::string>& list, const char* what) {
    std::vector<std::string> miss;
    for (const auto& h : have)
      if (std::find(list.begin(), list.end(), h) == list.end()) miss.push_back(h);
    if (!miss.empty()) {
      std::string s;
      for (size_t k = 0; k < miss.size(); ++k) s += (k ? ", " : "") + miss[k];
      throw std::invalid_argument(std::string("the shared type list has no ") + what + " " + s);
    }
  };
  check(own.beads, t.beads, "bead type");
  check(own.bonds, t.bonds, "bond type");
  check(own.angles, t.angles, "angle type");
  check(own.dihedrals, t.dihedrals, "dihedral type");
  System s;
  s.title = "coarse-grained " + (aa.title.empty() ? std::string("structure") : aa.title);
  s.cell = cell;
  s.has_mol = true;
  // a type's mass: the most common bead mass of that kind (chain-end beads carry their end group's atoms too)
  std::map<std::string, std::map<long, int>> masses;
  for (size_t b = 0; b < m.beads(); ++b) masses[m.bead_kind[b]][std::lround(m.bead_mass[b] * 1000)]++;
  for (size_t k = 0; k < t.beads.size(); ++k) {
    TypeInfo ti;
    ti.type = int(k + 1);
    ti.label = t.beads[k];
    if (masses.count(t.beads[k])) {
      const auto& h = masses[t.beads[k]];
      ti.mass = double(std::max_element(h.begin(), h.end(), [](auto& x, auto& y) { return x.second < y.second; })->first) / 1000.0;
    }
    s.types.push_back(ti);
  }
  for (size_t b = 0; b < m.beads(); ++b) {
    Atom a;
    a.id = int64_t(b + 1);
    a.mol = m.bead_mol[b] + 1;
    a.name = m.bead_kind[b];
    a.type = int(std::find(t.beads.begin(), t.beads.end(), m.bead_kind[b]) - t.beads.begin()) + 1;
    a.element = 0;
    a.pos = beads[b];
    s.atoms.push_back(a);
  }
  for (const auto& [x, y] : m.bonds) s.bonds.push_back({uint32_t(x), uint32_t(y), 1});
  return s;
}

CgMappingSummary cg_mapping_summary(const CgMapping& m, int sequences) {
  CgMappingSummary s;
  for (const auto& k : m.bead_kind) s.kinds[k]++;
  for (const auto& [k, c] : s.kinds) s.fraction[k] = double(c) / double(std::max<size_t>(1, m.beads()));
  s.mass_error = std::fabs(m.mass_cg - m.mass_aa);
  s.chains = int(m.chains.size());
  s.min_beads = 1 << 30;
  for (const auto& c : m.chains) s.min_beads = std::min(s.min_beads, int(c.size())), s.max_beads = std::max(s.max_beads, int(c.size()));
  if (m.chains.empty()) s.min_beads = 0;
  for (int k = 0; k < sequences && k < int(m.chains.size()); ++k) {
    std::string q;
    for (int b : m.chains[size_t(k)]) q += (q.empty() ? "" : " ") + m.bead_kind[size_t(b)];
    s.sequences.push_back(q);
  }
  return s;
}

void write_cg_lammps_data(const CgMapping& m, const System& cg, const CgTypes& types_in, const std::string& path) {
  const CgTypes own = cg_types_of(m);
  const CgTypes& t = types_in.empty() ? own : types_in;
  auto num = [](const std::vector<std::string>& list, const std::string& k) {
    const auto it = std::find(list.begin(), list.end(), k);
    if (it == list.end()) throw std::invalid_argument("the shared type list has no type " + k);
    return int(it - list.begin()) + 1;
  };
  auto kind = [&](int b) { return m.bead_kind[size_t(b)]; };
  std::ofstream f(path);
  if (!f) throw std::runtime_error("cannot write " + path);
  char line[256];
  f << "# CAPS coarse-grained beads (" << m.beads() << " beads, position: " << m.position << ")\n\n";
  f << cg.atoms.size() << " atoms\n" << t.beads.size() << " atom types\n";
  f << m.bonds.size() << " bonds\n" << std::max<size_t>(1, t.bonds.size()) << " bond types\n";
  f << m.angles.size() << " angles\n" << std::max<size_t>(1, t.angles.size()) << " angle types\n";
  f << m.dihedrals.size() << " dihedrals\n" << std::max<size_t>(1, t.dihedrals.size()) << " dihedral types\n\n";
  const Cell& c = cg.cell;
  if (c.valid()) {
    std::snprintf(line, sizeof line, "%.8f %.8f xlo xhi\n%.8f %.8f ylo yhi\n%.8f %.8f zlo zhi\n", c.origin[0], c.origin[0] + c.a[0], c.origin[1],
                  c.origin[1] + c.b[1], c.origin[2], c.origin[2] + c.c[2]);
    f << line;
    if (std::fabs(c.b[0]) > 1e-12 || std::fabs(c.c[0]) > 1e-12 || std::fabs(c.c[1]) > 1e-12) {
      std::snprintf(line, sizeof line, "%.8f %.8f %.8f xy xz yz\n", c.b[0], c.c[0], c.c[1]);
      f << line;
    }
  } else {
    double lo = 1e300, hi = -1e300;
    for (const auto& a : cg.atoms) for (int k = 0; k < 3; ++k) lo = std::min(lo, a.pos[size_t(k)]), hi = std::max(hi, a.pos[size_t(k)]);
    std::snprintf(line, sizeof line, "%.8f %.8f xlo xhi\n%.8f %.8f ylo yhi\n%.8f %.8f zlo zhi\n", lo - 10, hi + 10, lo - 10, hi + 10, lo - 10, hi + 10);
    f << line;
  }
  f << "\nMasses\n\n";
  for (size_t k = 0; k < t.beads.size(); ++k) {
    double mass = 0;
    for (const auto& ti : cg.types) if (ti.label == t.beads[k]) mass = ti.mass;
    std::snprintf(line, sizeof line, "%zu %.6f  # %s\n", k + 1, mass > 0 ? mass : 1.0, t.beads[k].c_str());
    f << line;
  }
  f << "\nAtoms  # full\n\n";
  for (size_t b = 0; b < cg.atoms.size(); ++b) {
    const auto& a = cg.atoms[b];
    std::snprintf(line, sizeof line, "%zu %lld %d 0.0 %.6f %.6f %.6f  # %s\n", b + 1, static_cast<long long>(a.mol), num(t.beads, a.name), a.pos[0], a.pos[1], a.pos[2], a.name.c_str());
    f << line;
  }
  auto section = [&](const char* title, const std::vector<std::string>& list, size_t count, auto&& term) {
    if (!count) return;
    f << "\n" << title << "\n\n";
    for (size_t k = 0; k < count; ++k) {
      const auto [key, ids] = term(k);
      f << k + 1 << " " << num(list, key) << ids << "\n";
    }
  };
  section("Bonds", t.bonds, m.bonds.size(), [&](size_t k) {
    const auto& [x, y] = m.bonds[k];
    return std::make_pair(cg_key({kind(x), kind(y)}), " " + std::to_string(x + 1) + " " + std::to_string(y + 1));
  });
  section("Angles", t.angles, m.angles.size(), [&](size_t k) {
    const auto& q = m.angles[k];
    return std::make_pair(cg_key({kind(q[0]), kind(q[1]), kind(q[2])}), " " + std::to_string(q[0] + 1) + " " + std::to_string(q[1] + 1) + " " + std::to_string(q[2] + 1));
  });
  section("Dihedrals", t.dihedrals, m.dihedrals.size(), [&](size_t k) {
    const auto& q = m.dihedrals[k];
    return std::make_pair(cg_key({kind(q[0]), kind(q[1]), kind(q[2]), kind(q[3])}),
                          " " + std::to_string(q[0] + 1) + " " + std::to_string(q[1] + 1) + " " + std::to_string(q[2] + 1) + " " + std::to_string(q[3] + 1));
  });
  // the type names, for the tables that go with them
  f << "\n# bond types:";
  for (size_t k = 0; k < t.bonds.size(); ++k) f << " " << k + 1 << "=" << t.bonds[k];
  f << "\n# angle types:";
  for (size_t k = 0; k < t.angles.size(); ++k) f << " " << k + 1 << "=" << t.angles[k];
  f << "\n# dihedral types:";
  for (size_t k = 0; k < t.dihedrals.size(); ++k) f << " " << k + 1 << "=" << t.dihedrals[k];
  f << "\n";
}

namespace {

// fast field splitting of one dump line
inline const char* next_field(const char* p, const char*& end) {
  while (*p == ' ' || *p == '\t') ++p;
  end = p;
  while (*end && *end != ' ' && *end != '\t' && *end != '\n' && *end != '\r') ++end;
  return p;
}

}  // namespace

CgTrajectoryReport map_lammps_dump(const CgMapping& m, const System& aa, const std::string& dump, const std::string& out, const CgTrajectoryOptions& o) {
  const auto t0 = std::chrono::steady_clock::now();
  CgTrajectoryReport rep;
  std::ifstream in(dump, std::ios::binary);
  if (!in) throw std::runtime_error("cannot open " + dump);
  in.seekg(0, std::ios::end);
  const double size = double(std::max<std::streamoff>(1, in.tellg()));
  in.seekg(0);
  std::ofstream w(out);
  if (!w) throw std::runtime_error("cannot write " + out);
  const CgTypes own = cg_types_of(m);
  const CgTypes& t = o.types.empty() ? own : o.types;
  std::vector<int> btype(m.beads());
  for (size_t b = 0; b < m.beads(); ++b) {
    const auto it = std::find(t.beads.begin(), t.beads.end(), m.bead_kind[b]);
    if (it == t.beads.end()) throw std::invalid_argument("the shared type list has no bead type " + m.bead_kind[b]);
    btype[b] = int(it - t.beads.begin()) + 1;
  }
  // atom id → index
  int64_t maxid = 0;
  for (size_t a = 0; a < aa.atoms.size(); ++a) maxid = std::max(maxid, aa.atoms[a].id > 0 ? aa.atoms[a].id : int64_t(a + 1));
  std::vector<int> index(size_t(maxid + 1), -1);
  for (size_t a = 0; a < aa.atoms.size(); ++a) index[size_t(aa.atoms[a].id > 0 ? aa.atoms[a].id : int64_t(a + 1))] = int(a);
  // chain walk for whole chains: each bead's predecessor along the bonds
  std::vector<int> parent(m.beads(), -1);
  {
    std::vector<std::vector<int>> nbr(m.beads());
    for (const auto& [x, y] : m.bonds) nbr[size_t(x)].push_back(y), nbr[size_t(y)].push_back(x);
    std::vector<char> seen(m.beads(), 0);
    for (const auto& ch : m.chains)
      for (int s : ch) {
        if (seen[size_t(s)]) continue;
        seen[size_t(s)] = 1;
        std::vector<int> q{s};
        for (size_t h = 0; h < q.size(); ++h)
          for (int v : nbr[size_t(q[h])])
            if (!seen[size_t(v)]) seen[size_t(v)] = 1, parent[size_t(v)] = q[h], q.push_back(v);
      }
  }
  std::vector<int> order;   // beads parents first
  for (const auto& ch : m.chains) {
    std::vector<int> q;
    for (int b : ch) if (parent[size_t(b)] < 0) q.push_back(b);
    std::vector<std::vector<int>> kids(m.beads());
    for (int b : ch) if (parent[size_t(b)] >= 0) kids[size_t(parent[size_t(b)])].push_back(b);
    for (size_t h = 0; h < q.size(); ++h) for (int k : kids[size_t(q[h])]) q.push_back(k);
    order.insert(order.end(), q.begin(), q.end());
  }
  std::vector<Vec3> pos(aa.atoms.size());
  std::vector<char> got(aa.atoms.size(), 0);
  std::string line;
  size_t frame = 0;
  char buf[512];
  while (std::getline(in, line)) {
    if (line.rfind("ITEM: TIMESTEP", 0) != 0) continue;
    std::string ts;
    std::getline(in, ts);
    std::getline(in, line);   // ITEM: NUMBER OF ATOMS
    std::getline(in, line);
    const size_t natoms = size_t(std::strtoull(line.c_str(), nullptr, 10));
    std::string boxhead;
    std::getline(in, boxhead);   // ITEM: BOX BOUNDS …
    std::string bl[3];
    double lo[3], hi[3], tilt[3] = {0, 0, 0};
    const bool tri = boxhead.find("xy") != std::string::npos;
    for (int k = 0; k < 3; ++k) {
      std::getline(in, bl[k]);
      char* e = nullptr;
      lo[k] = std::strtod(bl[k].c_str(), &e);
      hi[k] = std::strtod(e, &e);
      if (tri) tilt[k] = std::strtod(e, &e);
    }
    std::getline(in, line);   // ITEM: ATOMS …
    std::vector<std::string> cols;
    {
      std::istringstream ss(line.substr(std::min(line.size(), size_t(11))));
      for (std::string c; ss >> c;) cols.push_back(c);
    }
    const bool want = frame >= o.first && frame <= o.last && (frame - o.first) % std::max<size_t>(1, o.stride) == 0;
    if (!want) {
      for (size_t k = 0; k < natoms; ++k) std::getline(in, line);
      ++frame;
      ++rep.frames_read;
      if (frame > o.last) break;
      continue;
    }
    auto col = [&](const char* name) { for (size_t k = 0; k < cols.size(); ++k) if (cols[k] == name) return int(k); return -1; };
    const int cid = col("id");
    int cx = col("xu"), cy = col("yu"), cz = col("zu");
    int mode = 0;   // 0 unwrapped, 1 wrapped (+ images), 2 scaled, 3 scaled unwrapped
    if (cx < 0) { cx = col("x"), cy = col("y"), cz = col("z"); mode = 1; }
    if (cx < 0) { cx = col("xsu"), cy = col("ysu"), cz = col("zsu"); mode = 3; }
    if (cx < 0) { cx = col("xs"), cy = col("ys"), cz = col("zs"); mode = 2; }
    if (cid < 0 || cx < 0 || cy < 0 || cz < 0) throw std::invalid_argument("the dump needs id and xu yu zu, x y z, xs ys zs or xsu ysu zsu columns");
    const int ix = col("ix"), iy = col("iy"), iz = col("iz");
    const bool images = ix >= 0 && iy >= 0 && iz >= 0;
    // the cell (LAMMPS restricted triclinic: bounds include the tilts)
    double xlo = lo[0], xhi = hi[0], ylo = lo[1], yhi = hi[1];
    if (tri) {
      xlo -= std::min({0.0, tilt[0], tilt[1], tilt[0] + tilt[1]});
      xhi -= std::max({0.0, tilt[0], tilt[1], tilt[0] + tilt[1]});
      ylo -= std::min(0.0, tilt[2]);
      yhi -= std::max(0.0, tilt[2]);
    }
    Cell cell;
    cell.origin = {xlo, ylo, lo[2]};
    cell.a = {xhi - xlo, 0, 0};
    cell.b = {tilt[0], yhi - ylo, 0};
    cell.c = {tilt[1], tilt[2], hi[2] - lo[2]};
    std::fill(got.begin(), got.end(), 0);
    const int ncol = int(cols.size());
    std::vector<double> v(static_cast<size_t>(ncol));
    for (size_t k = 0; k < natoms; ++k) {
      if (!std::getline(in, line)) throw std::invalid_argument("the dump ends inside frame " + std::to_string(frame + 1));
      const char* p = line.c_str();
      for (int c = 0; c < ncol; ++c) {
        const char* e;
        p = next_field(p, e);
        v[size_t(c)] = std::strtod(p, nullptr);
        p = e;
      }
      const int64_t id = int64_t(v[size_t(cid)]);
      if (id <= 0 || id > maxid || index[size_t(id)] < 0) throw std::invalid_argument("the dump has atom id " + std::to_string(id) + ", which the structure does not");
      const size_t a = size_t(index[size_t(id)]);
      Vec3 r{v[size_t(cx)], v[size_t(cy)], v[size_t(cz)]};
      if (mode >= 2) r = cell.to_cartesian(r) + cell.origin;
      if (mode == 1 && images) r = r + cell.a * v[size_t(ix)] + cell.b * v[size_t(iy)] + cell.c * v[size_t(iz)];
      pos[a] = r;
      got[a] = 1;
    }
    if (natoms != aa.atoms.size() || std::find(got.begin(), got.end(), 0) != got.end())
      throw std::invalid_argument("frame " + std::to_string(frame + 1) + " has " + std::to_string(natoms) + " atoms; the structure has " + std::to_string(aa.atoms.size()));
    const bool unwrapped = mode == 0 || mode == 3 || images;
    rep.unwrapped_input = unwrapped;
    rep.images = images;
    std::vector<Vec3> beads = cg_positions(m, aa, pos, cell);
    if (!unwrapped)
      for (int b : order)
        if (parent[size_t(b)] >= 0) beads[size_t(b)] = beads[size_t(parent[size_t(b)])] + cell.minimum_image(beads[size_t(b)] - beads[size_t(parent[size_t(b)])]);
    w << "ITEM: TIMESTEP\n" << ts << "\nITEM: NUMBER OF ATOMS\n" << m.beads() << "\n" << boxhead << "\n" << bl[0] << "\n" << bl[1] << "\n" << bl[2] << "\nITEM: ATOMS id mol type xu yu zu\n";
    for (size_t b = 0; b < m.beads(); ++b) {
      const int len = std::snprintf(buf, sizeof buf, "%zu %d %d %.4f %.4f %.4f\n", b + 1, m.bead_mol[b] + 1, btype[b], beads[b][0], beads[b][1], beads[b][2]);
      w.write(buf, len);
    }
    ++rep.frames_written;
    ++rep.frames_read;
    ++frame;
    if (o.progress && !o.progress(double(in.tellg()) / size, rep.frames_written)) { rep.stopped = true; break; }
    if (frame > o.last) break;
  }
  if (!rep.frames_written) rep.notes.push_back("no frame of the dump was selected");
  if (!rep.unwrapped_input && rep.frames_written) rep.notes.push_back("wrapped coordinates without image flags: chains made whole along their bead bonds");
  rep.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  return rep;
}

size_t for_each_dump_frame(const std::string& path, size_t natoms,
                           const std::function<bool(size_t, int64_t, const std::vector<Vec3>&, const Cell&, bool)>& f, size_t stride) {
  std::ifstream in(path, std::ios::binary);
  if (!in) throw std::runtime_error("cannot open " + path);
  std::vector<Vec3> pos(natoms);
  std::vector<char> got(natoms);
  std::string line;
  size_t frame = 0, passed = 0;
  stride = std::max<size_t>(1, stride);
  while (std::getline(in, line)) {
    if (line.rfind("ITEM: TIMESTEP", 0) != 0) continue;
    std::getline(in, line);
    const int64_t ts = std::strtoll(line.c_str(), nullptr, 10);
    std::getline(in, line);
    std::getline(in, line);
    const size_t n = size_t(std::strtoull(line.c_str(), nullptr, 10));
    std::string boxhead;
    std::getline(in, boxhead);
    double lo[3], hi[3], tilt[3] = {0, 0, 0};
    const bool tri = boxhead.find("xy") != std::string::npos;
    for (int k = 0; k < 3; ++k) {
      std::getline(in, line);
      char* e = nullptr;
      lo[k] = std::strtod(line.c_str(), &e);
      hi[k] = std::strtod(e, &e);
      if (tri) tilt[k] = std::strtod(e, &e);
    }
    std::getline(in, line);
    std::vector<std::string> cols;
    {
      std::istringstream ss(line.substr(std::min(line.size(), size_t(11))));
      for (std::string c; ss >> c;) cols.push_back(c);
    }
    const bool want = frame % stride == 0;
    if (!want) {
      for (size_t k = 0; k < n; ++k) std::getline(in, line);
      ++frame;
      continue;
    }
    if (n != natoms) throw std::invalid_argument(path + ": frame " + std::to_string(frame + 1) + " has " + std::to_string(n) + " atoms, not " + std::to_string(natoms));
    auto col = [&](const char* name) { for (size_t k = 0; k < cols.size(); ++k) if (cols[k] == name) return int(k); return -1; };
    const int cid = col("id");
    int cx = col("xu"), cy = col("yu"), cz = col("zu"), mode = 0;
    if (cx < 0) { cx = col("x"), cy = col("y"), cz = col("z"); mode = 1; }
    if (cx < 0) { cx = col("xsu"), cy = col("ysu"), cz = col("zsu"); mode = 3; }
    if (cx < 0) { cx = col("xs"), cy = col("ys"), cz = col("zs"); mode = 2; }
    if (cid < 0 || cx < 0 || cy < 0 || cz < 0) throw std::invalid_argument(path + ": the dump needs id and xu yu zu, x y z, xs ys zs or xsu ysu zsu");
    const int ix = col("ix"), iy = col("iy"), iz = col("iz");
    const bool images = ix >= 0 && iy >= 0 && iz >= 0;
    double xlo = lo[0], xhi = hi[0], ylo = lo[1], yhi = hi[1];
    if (tri) {
      xlo -= std::min({0.0, tilt[0], tilt[1], tilt[0] + tilt[1]});
      xhi -= std::max({0.0, tilt[0], tilt[1], tilt[0] + tilt[1]});
      ylo -= std::min(0.0, tilt[2]);
      yhi -= std::max(0.0, tilt[2]);
    }
    Cell cell;
    cell.origin = {xlo, ylo, lo[2]};
    cell.a = {xhi - xlo, 0, 0};
    cell.b = {tilt[0], yhi - ylo, 0};
    cell.c = {tilt[1], tilt[2], hi[2] - lo[2]};
    std::fill(got.begin(), got.end(), 0);
    std::vector<double> v(cols.size());
    for (size_t k = 0; k < n; ++k) {
      if (!std::getline(in, line)) throw std::invalid_argument(path + ": the dump ends inside frame " + std::to_string(frame + 1));
      const char* p = line.c_str();
      for (size_t c = 0; c < cols.size(); ++c) {
        const char* e;
        p = next_field(p, e);
        v[c] = std::strtod(p, nullptr);
        p = e;
      }
      const int64_t id = int64_t(v[size_t(cid)]);
      if (id < 1 || size_t(id) > natoms) throw std::invalid_argument(path + ": atom id " + std::to_string(id) + " outside 1 … " + std::to_string(natoms));
      Vec3 r{v[size_t(cx)], v[size_t(cy)], v[size_t(cz)]};
      if (mode >= 2) r = cell.to_cartesian(r) + cell.origin;
      if (mode == 1 && images) r = r + cell.a * v[size_t(ix)] + cell.b * v[size_t(iy)] + cell.c * v[size_t(iz)];
      pos[size_t(id - 1)] = r;
      got[size_t(id - 1)] = 1;
    }
    if (std::find(got.begin(), got.end(), 0) != got.end()) throw std::invalid_argument(path + ": frame " + std::to_string(frame + 1) + " misses atoms");
    ++passed;
    const bool wrapped = !(mode == 0 || mode == 3 || images);
    if (!f(frame, ts, pos, cell, wrapped)) break;
    ++frame;
  }
  return passed;
}

}  // namespace caps
