// Martini 3 proteins, as martinize2 builds them (vermouth-martinize, force field martini3001). The model is data:
// data/martini/martini3-protein.json (bench/ff/convert_vermouth_martini3.py) holds the residue blocks, modifications,
// links and mappings of vermouth's files. This file applies them with vermouth's rules: residues mapped onto beads
// (mass-weighted centres of the mapped atoms), modifications for termini and protonation states, the links matched as
// vermouth's DoLinks matches them (node-induced subgraph isomorphism of the link graph, with the node attributes, the
// residue order, non-edges, patterns and the molecule's meta), their interactions added or replaced and removed in file
// order, parameters measured from the beads (the side-chain fix dihedrals), and the rubber-band elastic network.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <functional>
#include <map>
#include <memory>
#include <numeric>
#include <queue>
#include <set>
#include <sstream>
#include <stdexcept>

#include "caps/analysis.hpp"
#include "caps/elements.hpp"
#include "caps/json.hpp"
#include "caps/martini_protein.hpp"

namespace caps {

namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr double kDeg = kPi / 180;
constexpr double kKJ = 4.184;

const Json& model3(const std::string& path) {
  static std::map<std::string, std::unique_ptr<Json>> cache;
  auto& slot = cache[path];
  if (!slot) {
    std::ifstream in(path);
    if (!in) throw std::invalid_argument("cannot open " + path);
    std::stringstream ss;
    ss << in.rdbuf();
    slot = std::make_unique<Json>(Json::parse(ss.str()));
    if (slot->text("format") != "caps-vermouth-model") throw std::invalid_argument(path + ": not a CAPS vermouth model");
  }
  return *slot;
}

using Attrs = std::map<std::string, Json>;

bool same(const Json& a, const Json& b) {
  if (a.is_number() && b.is_number()) return a.number() == b.number();
  return a.dump(-1) == b.dump(-1);
}

// vermouth's attributes_match: every template attribute equals the node's (an absent attribute is null), or is a choice
// the node's value is in
bool attr_match(const Attrs& node, const std::string& key, const Json& tv) {
  auto it = node.find(key);
  const Json* mv = it == node.end() ? nullptr : &it->second;
  if (tv.is_object() && tv.has("$choice")) {
    if (!mv) return false;
    for (const auto& c : tv["$choice"].items())
      if (same(*mv, c)) return true;
    return false;
  }
  if (!mv) return tv.is_null();
  return same(*mv, tv);
}

bool attrs_match(const Attrs& node, const Json& tmpl, std::initializer_list<const char*> ignore = {}) {
  if (!tmpl.is_object()) return true;
  for (const auto& [k, v] : tmpl.members()) {
    bool skip = false;
    for (const char* x : ignore) skip = skip || k == x;
    if (!skip && !attr_match(node, k, v)) return false;
  }
  return true;
}

// a link node or pattern atom against a molecule node (order, replace and modifications are not matched here)
bool node_match(const Attrs& node, const Json& tmpl) { return attrs_match(node, tmpl, {"order", "replace", "modifications"}); }

struct Inter {
  std::vector<int> atoms;
  std::vector<std::string> params;
  Json meta = Json::object();
};

struct Mol {
  std::vector<Attrs> nodes;
  std::vector<Vec3> pos;          // Å
  std::vector<std::set<int>> adj;
  std::map<std::string, std::vector<Inter>> inter;
  std::vector<std::string> inter_order;   // interaction types in the order they first appear
  std::vector<int> comp;          // connected component (vermouth molecule)
  Attrs meta;
  void add_edge(int a, int b) {
    if (a == b) return;
    adj[size_t(a)].insert(b);
    adj[size_t(b)].insert(a);
  }
  std::vector<Inter>& list(const std::string& t) {
    if (!inter.count(t)) inter_order.push_back(t);
    return inter[t];
  }
};

int resid_of(const Attrs& a) { return int(a.at("resid").number()); }

// vermouth's match_order: the residue order a link asks for ("+"/"-" as numbers, ">" / "<" / "*" series)
struct Order { int type; int value; };   // type 0 number, 1 '><', 2 '*'
Order interpret(const Json& o) {
  if (o.is_number()) return {0, int(o.number())};
  const std::string s = o.str();
  if (s.empty()) throw std::invalid_argument("empty link order");
  if (s[0] == '>') return {1, int(s.size())};
  if (s[0] == '<') return {1, -int(s.size())};
  return {2, int(s.size())};
}
int sgn(int x) { return (x > 0) - (x < 0); }
bool match_order(const Json& o1, int r1, const Json& o2, int r2) {
  const Order a = interpret(o1), b = interpret(o2);
  if (a.type == 0) {
    if (b.type == 0) return (b.value - a.value) == (r2 - r1);
    if (a.value == 0) {
      if (b.type == 1 && sgn(r2 - r1) != sgn(b.value)) return false;
      if (b.type == 2 && r1 == r2) return false;
    }
  } else if (a.type == 1) {
    if (b.type == 0 && b.value == 0 && sgn(r1 - r2) != sgn(a.value)) return false;
    if (b.type == 1 && sgn(r2 - r1) != sgn(b.value - a.value)) return false;
  } else {
    if (b.type == 0 && b.value == 0 && r1 == r2) return false;
    if (b.type == 2 && ((a.value == b.value) != (r1 == r2))) return false;
  }
  return true;
}

double dihedral_deg(const Vec3& p0, const Vec3& p1, const Vec3& p2, const Vec3& p3) {
  const Vec3 ab = p1 - p0, bc = p2 - p1, cd = p3 - p2;
  const Vec3 n1 = cross(ab, bc), n2 = cross(bc, cd);
  const double y = dot(n1, cd) * norm(bc), x = dot(n1, n2);
  return std::atan2(y, x) / kDeg;
}

// a parameter the link measures from the beads (vermouth's LinkParameterEffector), formatted as vermouth formats it
std::string effector(const Json& e, const Mol& m, const std::map<std::string, int>& match) {
  const std::string name = e["$effector"].str();
  std::vector<Vec3> p;
  for (const auto& k : e["keys"].items()) p.push_back(m.pos[size_t(match.at(k.str()))] * 0.1);   // nm
  double v = 0;
  if (name == "dihphase" || name == "dihedral") {
    if (p.size() != 4) throw std::invalid_argument(name + " needs four atoms");
    v = dihedral_deg(p[0], p[1], p[2], p[3]);
    if (name == "dihphase") {
      v -= 180;
      if (v > 180) v -= 360;
      if (v < -180) v += 360;
    }
  } else if (name == "angle") {
    const Vec3 a = p[0] - p[1], c = p[2] - p[1];
    v = std::acos(std::clamp(dot(a, c) / (norm(a) * norm(c)), -1.0, 1.0)) / kDeg;
  } else if (name == "distance") {
    v = norm(p[1] - p[0]);
  } else {
    throw std::invalid_argument("parameter effector " + name + " is not known");
  }
  if (e.has("format") && !e["format"].is_null()) {
    std::string f = e["format"].str();   // Python's format spec: ".01f", ".3f" ...
    char buf[64];
    std::snprintf(buf, sizeof buf, ("%" + f).c_str(), v);
    return buf;
  }
  char buf[64];
  std::snprintf(buf, sizeof buf, "%.17g", v);
  return buf;
}

std::vector<std::string> params_of(const Json& ps, const Mol& m, const std::map<std::string, int>& match) {
  std::vector<std::string> out;
  for (const auto& p : ps.items()) out.push_back(p.is_object() ? effector(p, m, match) : p.str());
  return out;
}

// every match of a link in the molecule (link key → molecule node), vermouth's match_link
void match_link(const Mol& m, const Json& link, const std::function<void(const std::map<std::string, int>&)>& found) {
  if (link.has("molmeta") && !attrs_match(m.meta, link["molmeta"])) return;
  const Json& N = link["nodes"];
  std::vector<std::string> keys;
  for (const auto& [k, v] : N.members()) keys.push_back(k);
  const size_t nk = keys.size();
  if (!nk) return;
  std::map<std::string, size_t> kidx;
  for (size_t i = 0; i < nk; ++i) kidx[keys[i]] = i;
  std::vector<std::vector<char>> ladj(nk, std::vector<char>(nk, 0));
  for (const auto& e : link["edges"].items()) {
    const size_t a = kidx.at(e[0].str()), b = kidx.at(e[1].str());
    ladj[a][b] = ladj[b][a] = 1;
  }
  // search order: breadth first from the first node of each component, so later nodes have a matched neighbour
  std::vector<size_t> order;
  std::vector<int> anchor(nk, -1);
  std::vector<char> seen(nk, 0);
  for (size_t s = 0; s < nk; ++s) {
    if (seen[s]) continue;
    std::queue<size_t> q;
    q.push(s);
    seen[s] = 1;
    while (!q.empty()) {
      const size_t u = q.front();
      q.pop();
      order.push_back(u);
      for (size_t v = 0; v < nk; ++v)
        if (ladj[u][v] && !seen[v]) seen[v] = 1, anchor[v] = int(u), q.push(v);
    }
  }
  std::vector<const Json*> tmpl(nk);
  for (size_t i = 0; i < nk; ++i) tmpl[i] = &N[keys[i]];
  const size_t n = m.nodes.size();
  std::vector<int> img(nk, -1);
  std::vector<char> used(n, 0);
  auto check_rest = [&]() -> bool {
    std::map<std::string, int> match;
    for (size_t i = 0; i < nk; ++i) match[keys[i]] = img[i];
    // one molecule (vermouth applies links per molecule)
    for (size_t i = 1; i < nk; ++i)
      if (m.comp[size_t(img[i])] != m.comp[size_t(img[0])]) return false;
    // non-edges: no neighbour of the node that is the described atom of the residue at that order
    for (const auto& ne : link["non_edges"].items()) {
      const std::string from = ne[0].str();
      if (!match.count(from)) continue;
      const int f = match.at(from);
      const int fres = resid_of(m.nodes[size_t(f)]);
      const Json& to = ne[1];
      const int to_order = to.has("order") && to["order"].is_number() ? int(to["order"].number()) : 0;
      for (int nb : m.adj[size_t(f)])
        if (resid_of(m.nodes[size_t(nb)]) == fres + to_order && node_match(m.nodes[size_t(nb)], to)) return false;
    }
    // patterns: at least one must match when there are any
    const auto& pats = link["patterns"].items();
    if (!pats.empty()) {
      bool any = false;
      for (const auto& p : pats) {
        bool ok = true;
        for (const auto& at : p.items()) {
          const std::string k = at[0].str();
          if (!match.count(k)) { ok = false; break; }
          if (!node_match(m.nodes[size_t(match.at(k))], at[1])) { ok = false; break; }
        }
        if (ok) { any = true; break; }
      }
      if (!any) return false;
    }
    // residue orders: nodes with the same order share a residue; the orders agree pairwise
    std::vector<std::pair<Json, int>> om;
    for (size_t i = 0; i < nk; ++i) {
      if (!tmpl[i]->has("order")) continue;
      const Json& o = (*tmpl[i])["order"];
      const int r = resid_of(m.nodes[size_t(img[i])]);
      bool have = false;
      for (const auto& [oo, rr] : om)
        if (same(oo, o)) {
          have = true;
          if (rr != r) return false;
        }
      if (!have) om.push_back({o, r});
    }
    for (size_t a = 0; a < om.size(); ++a)
      for (size_t b = a + 1; b < om.size(); ++b)
        if (!match_order(om[a].first, om[a].second, om[b].first, om[b].second)) return false;
    found(match);
    return true;
  };
  std::function<void(size_t)> rec = [&](size_t depth) {
    if (depth == nk) {
      check_rest();
      return;
    }
    const size_t li = order[depth];
    std::vector<int> cand;
    if (anchor[li] >= 0) cand.assign(m.adj[size_t(img[size_t(anchor[li])])].begin(), m.adj[size_t(img[size_t(anchor[li])])].end());
    else {
      cand.resize(n);
      std::iota(cand.begin(), cand.end(), 0);
    }
    for (int c : cand) {
      if (used[size_t(c)] || !node_match(m.nodes[size_t(c)], *tmpl[li])) continue;
      // induced: edges to every matched node exactly as in the link
      bool ok = true;
      for (size_t d = 0; d < depth && ok; ++d) {
        const size_t lj = order[d];
        const bool me = m.adj[size_t(c)].count(img[lj]) > 0;
        ok = me == bool(ladj[li][lj]);
      }
      if (!ok) continue;
      img[li] = c;
      used[size_t(c)] = 1;
      rec(depth + 1);
      used[size_t(c)] = 0;
      img[li] = -1;
    }
  };
  rec(0);
}

int version_of(const Json& meta) { return meta.has("version") ? int(meta["version"].number()) : 0; }

// vermouth's add_or_replace_interaction: the same atoms in the same order and the same version replace
void add_or_replace(Mol& m, const std::string& type, Inter it) {
  auto& L = m.list(type);
  for (auto& x : L)
    if (x.atoms == it.atoms && version_of(x.meta) == version_of(it.meta)) {
      x = std::move(it);
      return;
    }
  L.push_back(std::move(it));
}

// vermouth's remove_matching_interaction: the first interaction with these atoms (and parameters, atom attributes and
// meta when the template gives them)
void remove_matching(Mol& m, const std::string& type, const Inter& t, const Json& atom_attrs) {
  if (!m.inter.count(type)) return;
  auto& L = m.inter[type];
  for (size_t k = 0; k < L.size(); ++k) {
    const Inter& x = L[k];
    if (x.atoms != t.atoms) continue;
    if (!t.params.empty() && x.params != t.params) continue;
    bool ok = true;
    for (size_t a = 0; a < x.atoms.size() && ok && a < atom_attrs.items().size(); ++a)
      ok = attrs_match(m.nodes[size_t(x.atoms[a])], atom_attrs[a]);
    if (!ok) continue;
    Attrs mm;
    for (const auto& [kk, vv] : x.meta.members()) mm[kk] = vv;
    if (!attrs_match(mm, t.meta)) continue;
    L.erase(L.begin() + long(k));
    return;
  }
}

void apply_links(Mol& m, const Json& links) {
  for (const auto& link : links.items()) {
    std::vector<std::map<std::string, int>> matches;
    match_link(m, link, [&](const std::map<std::string, int>& mt) { matches.push_back(mt); });
    for (const auto& mt : matches) {
      for (const auto& [k, v] : link["nodes"].members())
        if (v.has("replace"))
          for (const auto& [rk, rv] : v["replace"].members()) m.nodes[size_t(mt.at(k))][rk] = rv;
      for (const auto& [type, list] : link["removed"].members())
        for (const auto& it : list.items()) {
          Inter t;
          for (const auto& a : it["atoms"].items()) t.atoms.push_back(mt.at(a.str()));
          t.params = params_of(it["params"], m, mt);
          t.meta = it["meta"];
          remove_matching(m, type, t, it.has("atom_attrs") ? it["atom_attrs"] : Json::array());
        }
      for (const auto& [type, list] : link["interactions"].members())
        for (const auto& it : list.items()) {
          Inter t;
          for (const auto& a : it["atoms"].items()) t.atoms.push_back(mt.at(a.str()));
          t.params = params_of(it["params"], m, mt);
          t.meta = it["meta"];
          add_or_replace(m, type, std::move(t));
        }
    }
  }
}

double to_num(const Json& v) { return v.is_number() ? v.number() : std::stod(v.str()); }

// vermouth's atom masses (processors/attach_mass.py), whole numbers; other elements 30
double vermouth_mass(int z) {
  switch (z) {
    case 1: return 1;
    case 6: return 12;
    case 7: return 14;
    case 8: return 16;
    case 15: return 31;
    case 16: return 32;
    default: return 30;
  }
}

}  // namespace

System martini3_protein(const System& aa_in, const Martini3Options& opt, const std::string& data_path, MartiniProteinReport* rep_out) {
  MartiniProteinReport rep;
  const Json& M = model3(data_path);
  const Json& B = M["blocks"];
  const Json& MAP = M["mapping"];
  const Json& MODS = M["modifications"];
  const Json& MH = M["mod_h"];
  // the protein alone: residues with a Martini 3 block (waters, ions and other molecules are not converted)
  auto canon = [&](const std::string& n) -> std::string {
    static const std::map<std::string, std::string> alias = {{"HID", "HSD"}, {"HIE", "HSE"}, {"HIP", "HSP"}, {"CYX", "CYS"}, {"CYM", "CYS"},
                                                             {"ASH", "ASPP"}, {"GLH", "GLUP"}, {"LYN", "LSN"}, {"HSE", "HSE"}};
    auto it = alias.find(n);
    return it == alias.end() ? n : it->second;
  };
  System aa;
  {
    std::vector<int> keep(aa_in.atoms.size(), -1);
    for (size_t i = 0; i < aa_in.atoms.size(); ++i)
      if (MAP.has(canon(aa_in.atoms[i].resname))) {
        keep[i] = int(aa.atoms.size());
        aa.atoms.push_back(aa_in.atoms[i]);
      }
    for (const auto& b : aa_in.bonds)
      if (keep[b.i] >= 0 && keep[b.j] >= 0) aa.bonds.push_back({uint32_t(keep[b.i]), uint32_t(keep[b.j]), b.order});
    aa.cell = aa_in.cell;
    aa.title = aa_in.title;
    aa.has_mol = false;   // residues are split into chains by their bonds
    // martinize2 takes bonds from the file and from distances (-bonds-from both): a PDB whose CONECT records list only
    // its disulfides still has its peptide bonds
    std::set<std::pair<uint32_t, uint32_t>> have;
    for (const auto& b : aa.bonds) have.insert({std::min(b.i, b.j), std::max(b.i, b.j)});
    for (const auto& b : perceive_bonds(aa))
      if (have.insert({std::min(b.i, b.j), std::max(b.i, b.j)}).second) aa.bonds.push_back(b);
    if (aa.atoms.size() < aa_in.atoms.size())
      rep.notes.push_back(std::to_string(aa_in.atoms.size() - aa.atoms.size()) + " atoms not in amino acids (water, ions, other molecules) left out");
  }
  auto res = protein_residues(aa);
  const size_t nres = res.size();
  if (!nres) throw std::invalid_argument("no amino-acid residues");
  // secondary structure: DSSP, given letters (one letter for every residue), or none ("-")
  bool have_ss = opt.ss != "-";
  std::string letters;
  if (have_ss) {
    letters = opt.ss.empty() ? dssp(aa) : opt.ss.size() == 1 ? std::string(nres, opt.ss[0]) : opt.ss;
    if (letters.size() != nres)
      throw std::invalid_argument("the secondary structure has " + std::to_string(letters.size()) + " letters for " + std::to_string(nres) + " residues");
  }
  std::string cg = have_ss ? dssp_to_martini(letters, data_path) : std::string(nres, ' ');
  rep.dssp = letters;
  rep.cg_ss = have_ss ? cg : "";
  const auto nb = aa.neighbours();
  auto atom_named = [&](const ProteinResidue& r, const std::string& name) -> int {
    for (uint32_t a : r.atoms)
      if (aa.atoms[a].name == name) return int(a);
    return -1;
  };
  auto nh = [&](int a) {
    int h = 0;
    if (a >= 0) for (uint32_t v : nb[size_t(a)]) h += aa.atoms[v].element == 1;
    return h;
  };
  std::vector<char> nter(nres), cter(nres);
  for (size_t k = 0; k < nres; ++k) {
    nter[k] = k == 0 || res[k - 1].chain != res[k].chain;
    cter[k] = k + 1 == nres || res[k + 1].chain != res[k].chain;
  }
  // blocks and modifications: protonation from the hydrogens present, termini charged unless neutral ones are asked for
  // (or the structure's hydrogens say so)
  std::vector<std::string> block(nres);
  std::vector<std::vector<std::string>> mods(nres);
  for (size_t k = 0; k < nres; ++k) {
    const std::string n = canon(res[k].name);
    bool anyh = false;
    for (uint32_t a : res[k].atoms) anyh = anyh || aa.atoms[a].element == 1;
    block[k] = n;
    auto hon = [&](const char* x) { return nh(atom_named(res[k], x)) > 0; };
    if (n == "HIS" && anyh) {
      const bool d = hon("ND1"), e = hon("NE2");
      mods[k].push_back(d && e ? "HIS-HP" : d ? "HIS-HD" : "HIS-HE");
    } else if (n == "ASP" && anyh) {
      if (hon("OD1")) mods[k].push_back("ASP-HD1");
      else if (hon("OD2")) mods[k].push_back("ASP-HD2");
    } else if (n == "GLU" && anyh) {
      if (hon("OE1")) mods[k].push_back("GLU-HE1");
      else if (hon("OE2")) mods[k].push_back("GLU-HE2");
    } else if (n == "LYS" && anyh) {
      const int h = nh(atom_named(res[k], "NZ"));
      if (h == 2) mods[k].push_back("LYS-LSN");
      else if (h == 3) mods[k].push_back("LYS-HZ3");
    }
    if (!B.has(block[k]) || !MAP.has(block[k])) throw std::invalid_argument("residue " + res[k].name + std::to_string(res[k].resid) + ": Martini 3 has no " + block[k]);
    const int oxt = atom_named(res[k], "OXT") >= 0 ? atom_named(res[k], "OXT") : atom_named(res[k], "OT2");
    if (nter[k]) mods[k].push_back(opt.neutral_termini || (anyh && nh(res[k].n) == (n == "PRO" ? 1 : 2)) ? "NH2-ter" : "N-ter");
    if (cter[k]) mods[k].push_back(opt.neutral_termini || nh(oxt) > 0 || (res[k].o >= 0 && nh(res[k].o) > 0) ? "COOH-ter" : "C-ter");
  }
  std::vector<int> res_idx(aa.atoms.size(), -1);
  for (size_t k = 0; k < nres; ++k)
    for (uint32_t a : res[k].atoms) res_idx[a] = int(k);
  auto res_of = [&](uint32_t a) { return res_idx[a]; };
  // beads: the block's atoms, at the mass-weighted centre of their mapped atoms (vermouth's DoAverageBead)
  Mol m;
  std::vector<std::map<std::string, int>> bead(nres);
  std::map<uint32_t, std::vector<int>> atom_beads;   // all-atom atom → beads it is part of (weight 0 included)
  // vermouth numbers residues 1 … n per molecule ("-resid mol"); molecules are found after the edges, so number per
  // residue here and renumber below
  for (size_t k = 0; k < nres; ++k) {
    const Json& blk = B[block[k]];
    const Json& mp = MAP[block[k]];
    std::map<std::string, std::pair<Vec3, double>> acc;
    auto add = [&](uint32_t a, const std::string& b, double w) {
      auto& s = acc[b];
      const double wm = w * vermouth_mass(aa.atoms[a].element);
      s.first = s.first + aa.atoms[a].pos * wm;
      s.second += wm;
    };
    std::vector<std::pair<uint32_t, std::string>> members;
    auto input_name = [&](const std::string& map_name) {
      if (map_name == "CD" && block[k] == "ILE" && atom_named(res[k], "CD") < 0) return std::string("CD1");
      if (map_name == "O" && atom_named(res[k], "O") < 0) return std::string("OT1");
      return map_name;
    };
    for (const auto& h : mp["heavy"].items()) {
      const std::string an = h[0].str();
      const int a = atom_named(res[k], input_name(an));
      if (a < 0) continue;   // missing in the structure: no position (vermouth's repaired atoms)
      for (const auto& e : h[1].items()) add(uint32_t(a), e[0].str(), e[1].number()), members.push_back({uint32_t(a), e[0].str()});
      // hydrogens on this atom weigh as the map's hydrogens on it do; hydrogens on the backbone N the residue's map
      // does not list (a terminal proline) join the backbone bead, as the N-terminal modification maps them
      const Json* hw = mp["h"].has(an) ? &mp["h"][an] : nullptr;
      // a modification mapping that lists this atom's hydrogens puts them in its bead with weight 1
      std::string mod_bead;
      for (const auto& md : mods[k])
        if (MH.has(md) && MH[md].has(an)) mod_bead = MH[md][an].str();
      for (uint32_t v : nb[size_t(a)]) {
        if (aa.atoms[v].element != 1) continue;
        if (!mod_bead.empty()) {
          add(v, mod_bead, 1.0), members.push_back({v, mod_bead});
        } else if (hw) {
          for (const auto& e : hw->items()) add(v, e[0].str(), e[1].number()), members.push_back({v, e[0].str()});
        } else if (aa.atoms[size_t(a)].element != 6) {
          // a hydrogen the reference residue does not have (a terminal NH3+, a protonated acid, a neutral lysine's
          // missing one): the modification mappings put it in the bead of its heavy atom
          const std::string bb = h[1].items().empty() ? "BB" : h[1][0][0].str();
          add(v, bb, 1.0), members.push_back({v, bb});
        }
      }
    }
    // the C-terminal oxygen (and its hydrogen), whatever its name (OXT, OT2, O01 ...): an oxygen on the backbone C
    // besides O joins the backbone bead (the C-ter / COOH-ter modification mappings)
    if (res[k].c >= 0)
      for (uint32_t a : nb[size_t(res[k].c)]) {
        if (aa.atoms[a].element != 8 || int(a) == res[k].o || res_of(a) != int(k)) continue;
        add(a, "BB", 1.0), members.push_back({a, "BB"});
        for (uint32_t v : nb[a])
          if (aa.atoms[v].element == 1) add(v, "BB", 1.0), members.push_back({v, "BB"});
      }
    for (const auto& at : blk["atoms"].items()) {
      const std::string an = at["atomname"].str();
      Attrs A;
      for (const auto& [kk, vv] : at.members()) A[kk] = vv;
      A["resid"] = Json(int(k + 1));
      A["chain"] = Json(res[k].chain);
      A["stash_resid"] = Json(double(res[k].resid));
      if (!A.count("charge")) A["charge"] = Json(0.0);
      if (have_ss) A["cgsecstruct"] = Json(std::string(1, cg[k])), A["aasecstruct"] = Json(std::string(1, letters[k]));
      for (const auto& [a, b] : opt.idr)
        if (res[k].resid >= a && res[k].resid <= b) {
          A["cgidr"] = Json(true);
          if (have_ss && cg[k] != 'C') A["cgsecstruct"] = Json("C"), cg[k] = 'C';
        }
      auto it = acc.find(an);
      if (it == acc.end() || std::fabs(it->second.second) < 1e-7)
        throw std::invalid_argument("residue " + res[k].name + std::to_string(res[k].resid) + ": no atoms for bead " + an);
      bead[k][an] = int(m.nodes.size());
      m.nodes.push_back(std::move(A));
      m.pos.push_back(it->second.first * (1 / it->second.second));
    }
    for (const auto& [a, b] : members)
      if (bead[k].count(b)) atom_beads[a].push_back(bead[k][b]);
  }
  const size_t n = m.nodes.size();
  m.adj.assign(n, {});
  // edges: within residues from the blocks, between residues where mapped atoms are bonded (the peptide bond, a
  // disulfide)
  for (size_t k = 0; k < nres; ++k)
    for (const auto& e : B[block[k]]["edges"].items()) m.add_edge(bead[k].at(e[0].str()), bead[k].at(e[1].str()));
  std::vector<int> res_of_atom(aa.atoms.size(), -1);
  for (size_t k = 0; k < nres; ++k)
    for (uint32_t a : res[k].atoms) res_of_atom[a] = int(k);
  for (const auto& b : aa.bonds) {
    if (res_of_atom[b.i] == res_of_atom[b.j] || res_of_atom[b.i] < 0 || res_of_atom[b.j] < 0) continue;
    if (!atom_beads.count(b.i) || !atom_beads.count(b.j)) continue;
    const bool ss_bond = aa.atoms[b.i].element == 16 && aa.atoms[b.j].element == 16;
    if (ss_bond && !opt.disulfides) continue;
    for (int x : atom_beads[b.i])
      for (int y : atom_beads[b.j]) m.add_edge(x, y);
    if (ss_bond) ++rep.disulfides;
  }
  // molecules (connected beads) and residues numbered 1 … n in each (vermouth's "-resid mol")
  m.comp.assign(n, -1);
  int ncomp = 0;
  for (size_t s = 0; s < n; ++s) {
    if (m.comp[s] >= 0) continue;
    std::queue<size_t> q;
    q.push(s);
    m.comp[s] = ncomp;
    while (!q.empty()) {
      const size_t u = q.front();
      q.pop();
      for (int v : m.adj[u])
        if (m.comp[size_t(v)] < 0) m.comp[size_t(v)] = ncomp, q.push(size_t(v));
    }
    ++ncomp;
  }
  {
    std::vector<int> next(size_t(ncomp), 0);
    for (size_t k = 0; k < nres; ++k) {
      const int c = m.comp[size_t(bead[k].begin()->second)];
      const int r = ++next[size_t(c)];
      for (const auto& [nm, i] : bead[k]) m.nodes[size_t(i)]["resid"] = Json(r);
    }
  }
  // the blocks' interactions
  for (size_t k = 0; k < nres; ++k)
    for (const auto& [type, list] : B[block[k]]["interactions"].members())
      for (const auto& it : list.items()) {
        Inter t;
        for (const auto& a : it["atoms"].items()) t.atoms.push_back(bead[k].at(a.str()));
        for (const auto& p : it["params"].items()) t.params.push_back(p.str());
        t.meta = it["meta"];
        m.list(type).push_back(std::move(t));
      }
  // modifications: the attributes they replace
  for (size_t k = 0; k < nres; ++k)
    for (const auto& mod : mods[k]) {
      if (!MODS.has(mod)) throw std::invalid_argument("modification " + mod + " is not in the model");
      for (const auto& [key, at] : MODS[mod]["atoms"].members()) {
        const std::string an = at.text("atomname", key);
        if (!bead[k].count(an)) continue;
        Attrs& node = m.nodes[size_t(bead[k][an])];
        if (!attrs_match(node, at, {"order", "replace", "modifications", "PTM_atom"})) continue;
        if (at.has("replace"))
          for (const auto& [rk, rv] : at["replace"].members()) node[rk] = rk == "charge" ? Json(to_num(rv)) : rv;
      }
    }
  m.meta["scfix"] = Json(opt.scfix);
  m.meta["extdih"] = Json(opt.extdih);
  m.meta["idr"] = Json(!opt.idr.empty());
  apply_links(m, M["links"]);
  // the rubber-band elastic network (vermouth's apply_rubber_band) between backbone beads
  int n_rubber = 0;
  if (opt.elastic) {
    std::vector<int> sel;
    for (size_t i = 0; i < n; ++i)
      if (m.nodes[i].at("atomname").str() == M["variables"].text("bb_atomname", "BB")) sel.push_back(int(i));
    // residues: beads of the same chain / resid / resname; separation along the residue graph
    std::map<std::tuple<int, int, std::string>, int> rid;
    std::vector<int> rnode(n);
    for (size_t i = 0; i < n; ++i) {
      auto key = std::make_tuple(int(m.nodes[i].at("chain").number()), resid_of(m.nodes[i]), m.nodes[i].at("resname").str());
      auto it = rid.emplace(key, int(rid.size())).first;
      rnode[i] = it->second;
    }
    const size_t nr = rid.size();
    std::vector<std::set<int>> radj(nr);
    for (size_t i = 0; i < n; ++i)
      for (int j : m.adj[i])
        if (rnode[i] != rnode[size_t(j)]) radj[size_t(rnode[i])].insert(rnode[size_t(j)]);
    const int rmd = opt.ermd >= 0 ? opt.ermd : 2;
    auto within = [&](int ra, int rb) {   // residue-graph distance ≤ rmd
      if (ra == rb) return true;
      std::vector<int> d(nr, -1);
      std::queue<int> q;
      q.push(ra);
      d[size_t(ra)] = 0;
      while (!q.empty()) {
        const int u = q.front();
        q.pop();
        if (d[size_t(u)] >= rmd) continue;
        for (int v : radj[size_t(u)])
          if (d[size_t(v)] < 0) {
            d[size_t(v)] = d[size_t(u)] + 1;
            if (v == rb) return true;
            q.push(v);
          }
      }
      return false;
    };
    for (size_t a = 0; a < sel.size(); ++a)
      for (size_t b = a + 1; b < sel.size(); ++b) {
        const int i = sel[a], j = sel[b];
        const double d = norm(m.pos[size_t(i)] - m.pos[size_t(j)]) * 0.1;   // nm
        double kf = std::exp(-opt.ea * std::pow(d - opt.es, opt.ep)) * opt.ef;
        if (kf < opt.em) kf = 0;
        if (kf > opt.ef) kf = opt.ef;
        if (d > opt.eu || d < opt.el) kf = 0;
        if (!(kf > opt.em)) continue;
        bool domain = true;
        if (opt.eunit == "molecule") domain = m.comp[size_t(i)] == m.comp[size_t(j)];
        else if (opt.eunit == "chain") domain = m.nodes[size_t(i)].at("chain").number() == m.nodes[size_t(j)].at("chain").number();
        else if (opt.eunit == "all") domain = true;
        if (!domain) continue;
        if (within(rnode[size_t(i)], rnode[size_t(j)])) continue;
        Inter t;
        t.atoms = {i, j};
        char len[32], kk[32];
        std::snprintf(len, sizeof len, "%.5f", std::round(d * 1e5) / 1e5);
        std::snprintf(kk, sizeof kk, "%.10g", kf);
        t.params = {M["variables"].text("elastic_network_bond_type", "1"), len, kk};
        t.meta = Json::object();
        t.meta["group"] = Json("Rubber band");
        m.list("bonds").push_back(std::move(t));
        ++n_rubber;
      }
    rep.notes.push_back(std::to_string(n_rubber) + " elastic bonds (" + std::to_string(opt.ef) + " kJ/mol/nm², " + std::to_string(opt.el) + "–" +
                        std::to_string(opt.eu) + " nm)");
  }
  // the beads and their explicit topology (CAPS units; constraints as stiff bonds, the FLEXIBLE variants left out)
  System out;
  out.title = aa_in.title;
  out.cell = aa_in.cell;
  out.source_format = "caps-martini3-protein";
  auto topo = std::make_shared<ExplicitTopology>();
  topo->source = "Martini 3 protein";
  topo->natoms = n;
  topo->masses.assign(n, std::nan(""));
  for (size_t i = 0; i < n; ++i) {
    const Attrs& A = m.nodes[i];
    Atom at;
    at.id = int64_t(i + 1);
    at.mol = m.comp[i] + 1;
    at.resname = A.at("resname").str();
    at.resid = int64_t(A.at("stash_resid").number());
    at.element = 0;
    at.name = A.at("atype").str();
    at.charge = to_num(A.at("charge"));
    at.pos = m.pos[i];
    out.atoms.push_back(at);
    if (A.count("mass")) topo->masses[i] = to_num(A.at("mass"));
  }
  if (std::all_of(topo->masses.begin(), topo->masses.end(), [](double x) { return std::isnan(x); })) topo->masses.clear();
  // one type entry per bead type, with Martini 3's mass for its size (regular 72, small 54, tiny 36)
  std::map<std::string, int> tid;
  for (auto& a : out.atoms) {
    auto [it, fresh] = tid.emplace(a.name, int(tid.size()) + 1);
    a.type = it->second;
    if (fresh) {
      TypeInfo t;
      t.type = it->second;
      t.label = a.name;
      t.mass = a.name[0] == 'T' ? 36.0 : a.name[0] == 'S' ? 54.0 : 72.0;
      out.types.push_back(t);
    }
  }
  out.has_charges = true;
  out.has_mol = true;
  auto flexible_only = [](const Json& meta) { return meta.text("ifdef") == "FLEXIBLE"; };
  auto group = [](const Json& meta, const std::string& def) { return meta.text("group", def); };
  const double stiff = opt.constraint_kj;
  auto P = [](const Inter& t, size_t k) { return k < t.params.size() ? std::stod(t.params[k]) : 0.0; };
  for (const auto& type : m.inter_order) {
    for (const auto& t : m.inter[type]) {
      if (flexible_only(t.meta)) continue;
      auto u = [&](size_t k) { return uint32_t(t.atoms[k]); };
      if (type == "bonds") {
        if (t.params.empty() || t.params[0] != "1") throw std::invalid_argument("bond function " + (t.params.empty() ? "?" : t.params[0]) + " is not handled");
        topo->bonds.push_back({u(0), u(1), P(t, 2) / (2 * kKJ * 100), P(t, 1) * 10, group(t.meta, "bond")});
        out.bonds.push_back({u(0), u(1), 1});
      } else if (type == "constraints") {
        topo->bonds.push_back({u(0), u(1), stiff / (2 * kKJ * 100), P(t, 1) * 10, "constraint: " + group(t.meta, "constraint")});
        out.bonds.push_back({u(0), u(1), 1});
      } else if (type == "angles") {
        const int f = int(P(t, 0));
        const int form = f == 1 ? 0 : f == 2 ? 1 : f == 10 ? 5 : -1;
        if (form < 0) throw std::invalid_argument("angle function " + t.params[0] + " is not handled");
        topo->angles.push_back({u(0), u(1), u(2), form, P(t, 2) / (2 * kKJ), P(t, 1) * kDeg, group(t.meta, "angle")});
      } else if (type == "dihedrals") {
        const int f = int(P(t, 0));
        if (f != 1 && f != 9 && f != 4) throw std::invalid_argument("dihedral function " + t.params[0] + " is not handled");
        topo->dihedrals.push_back({u(0), u(1), u(2), u(3), f, P(t, 2) / kKJ, P(t, 1) * kDeg, int(P(t, 3)), group(t.meta, "dihedral")});
      } else if (type == "impropers") {
        topo->dihedrals.push_back({u(0), u(1), u(2), u(3), 2, P(t, 2) / (2 * kKJ), P(t, 1) * kDeg, 0, group(t.meta, "improper")});
      } else if (type == "exclusions") {
        for (size_t k = 1; k < t.atoms.size(); ++k) topo->exclusions.push_back({u(0), u(k)});
      } else if (type == "virtual_sitesn") {
        ExplicitTopology::VSite v;
        v.site = u(0);
        for (size_t k = 1; k < t.atoms.size(); ++k) v.from.push_back(u(k));
        const int f = int(P(t, 0));
        if (f == 1) v.w.assign(v.from.size(), 1.0 / double(v.from.size()));   // centre of geometry
        else if (f != 2) throw std::invalid_argument("virtual_sitesn function " + t.params[0] + " is not handled");
        topo->vsites.push_back(v);   // function 2: the centre of mass
      } else {
        throw std::invalid_argument("interactions of type " + type + " are not handled");
      }
    }
  }
  out.topology = topo;
  out.bonds_from_file = true;
  rep.residues = int(nres);
  rep.beads = int(n);
  rep.chains = res.back().chain + 1;
  for (size_t k = 0; k < nres; ++k) rep.residue_names.push_back(m.nodes[size_t(bead[k].at("BB"))].at("resname").str());
  if (rep_out) *rep_out = rep;
  return out;
}

bool is_martini3_model(const std::string& path) {
  std::ifstream in(path);
  std::string head(512, '\0');
  in.read(&head[0], long(head.size()));
  return head.find("\"caps-vermouth-model\"") != std::string::npos;
}

std::string martini3_itp(const System& s, double constraint_kj) {
  if (!s.topology) throw std::invalid_argument("the structure has no explicit topology");
  const ExplicitTopology& T = *s.topology;
  const double kc = constraint_kj / (2 * kKJ * 100);
  std::ostringstream o;
  char b[256];
  o << "[ moleculetype ]\nprotein 1\n\n[ atoms ]\n";
  for (size_t i = 0; i < s.atoms.size(); ++i) {
    const auto& a = s.atoms[i];
    std::snprintf(b, sizeof b, "%4zu %-5s %4lld %-4s %-4s %4zu %6.2f", i + 1, a.name.c_str(), static_cast<long long>(a.resid), a.resname.c_str(), "B",
                  i + 1, a.charge);
    o << b;
    if (!T.masses.empty() && !std::isnan(T.masses[i])) {
      std::snprintf(b, sizeof b, " %.1f", T.masses[i]);
      o << b;
    }
    o << "\n";
  }
  auto is_constraint = [&](const ExplicitTopology::Bond& x) { return x.group.rfind("constraint", 0) == 0 && std::fabs(x.k - kc) < 1e-9; };
  o << "\n[ bonds ]\n";
  for (const auto& x : T.bonds)
    if (!is_constraint(x)) {
      std::snprintf(b, sizeof b, "%4u %4u 1 %.5f %.4f ; %s\n", x.i + 1, x.j + 1, x.r0 / 10, x.k * 2 * kKJ * 100, x.group.c_str());
      o << b;
    }
  o << "\n[ constraints ]\n";
  for (const auto& x : T.bonds)
    if (is_constraint(x)) {
      std::snprintf(b, sizeof b, "%4u %4u 1 %.5f ; %s\n", x.i + 1, x.j + 1, x.r0 / 10, x.group.c_str());
      o << b;
    }
  o << "\n[ angles ]\n";
  for (const auto& x : T.angles) {
    std::snprintf(b, sizeof b, "%4u %4u %4u %d %.2f %.2f ; %s\n", x.i + 1, x.j + 1, x.k + 1, x.form == 0 ? 1 : x.form == 1 ? 2 : 10, x.theta0 / kDeg,
                  x.kt * 2 * kKJ, x.group.c_str());
    o << b;
  }
  o << "\n[ dihedrals ]\n";
  for (const auto& x : T.dihedrals) {
    if (x.form == 2) std::snprintf(b, sizeof b, "%4u %4u %4u %4u 2 %.2f %.2f ; %s\n", x.i + 1, x.j + 1, x.k + 1, x.l + 1, x.phi0 / kDeg, x.kd * 2 * kKJ, x.group.c_str());
    else std::snprintf(b, sizeof b, "%4u %4u %4u %4u %d %.2f %.2f %d ; %s\n", x.i + 1, x.j + 1, x.k + 1, x.l + 1, x.form, x.phi0 / kDeg, x.kd * kKJ, x.n,
                       x.group.c_str());
    o << b;
  }
  if (!T.vsites.empty()) {
    o << "\n[ virtual_sitesn ]\n";
    for (const auto& v : T.vsites) {
      o << v.site + 1 << (v.w.empty() ? " 2" : " 1");
      for (uint32_t a : v.from) o << " " << a + 1;
      o << "\n";
    }
  }
  if (!T.exclusions.empty()) {
    o << "\n[ exclusions ]\n";
    for (const auto& e : T.exclusions) o << e.first + 1 << " " << e.second + 1 << "\n";
  }
  return o.str();
}

}  // namespace caps
