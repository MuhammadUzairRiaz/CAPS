// Martini 3 molecules: vermouth-martinize's martini3001 small molecules (data/martini/martini3-small-molecules.json,
// bench/ff/convert_vermouth_martini3_small.py) mapped from all-atom structures by graph, molecules given as GROMACS
// terms (data/martini/martini3-molecules.json), and their topology as a GROMACS .itp.
#include <tuple>
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
#include "caps/martini3_molecules.hpp"
#include "caps/resolution.hpp"

namespace caps {

namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr double kDeg = kPi / 180;
constexpr double kKJ = 4.184;

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

// One GROMACS-style term of a block or link in CAPS units: bonds (function 1), constraints (stiff bonds), angles (1, 2,
// 10), dihedrals (1, 9, 4), impropers (2), exclusions, virtual sites (virtual_sitesn 1 / 2, linear virtual_sites2 /
// virtual_sites3). Terms only for GROMACS's FLEXIBLE variant are left out (their constraints are used).
void add_term(ExplicitTopology& topo, System& out, const std::string& type, const std::vector<int>& atoms, const std::vector<std::string>& params,
              const Json& meta, double stiff) {
  if (!meta.text("ifdef").empty()) return;   // nothing is defined: FLEXIBLE, TI, POSRE_… variants are left out
  auto group = [&](const std::string& def) { return meta.text("group", def); };
  auto P = [&](size_t k) { return k < params.size() ? std::stod(params[k]) : 0.0; };
  auto u = [&](size_t k) { return uint32_t(atoms.at(k)); };
  if (type == "bonds") {
    if (params.empty() || params[0] != "1") throw std::invalid_argument("bond function " + (params.empty() ? "?" : params[0]) + " is not handled");
    topo.bonds.push_back({u(0), u(1), P(2) / (2 * kKJ * 100), P(1) * 10, group("bond")});
    out.bonds.push_back({u(0), u(1), 1});
  } else if (type == "constraints") {
    topo.bonds.push_back({u(0), u(1), stiff / (2 * kKJ * 100), P(1) * 10, "constraint: " + group("constraint")});
    out.bonds.push_back({u(0), u(1), 1});
  } else if (type == "angles") {
    const int f = int(P(0));
    const int form = f == 1 ? 0 : f == 2 ? 1 : f == 10 ? 5 : -1;
    if (form < 0) throw std::invalid_argument("angle function " + params[0] + " is not handled");
    topo.angles.push_back({u(0), u(1), u(2), form, P(2) / (2 * kKJ), P(1) * kDeg, group("angle")});
  } else if (type == "dihedrals" && int(P(0)) == 11) {   // combined bending–torsion: k a0 … a4
    ExplicitTopology::Dihedral d{u(0), u(1), u(2), u(3), 11, 0, 0, 0, group("dihedral")};
    for (int n = 0; n < 5; ++n) d.c[size_t(n)] = P(1) * P(size_t(2 + n)) / kKJ;
    topo.dihedrals.push_back(d);
  } else if (type == "pairs") {   // an explicit LJ pair: function 1, σ, ε
    if (int(P(0)) != 1 || params.size() < 3) throw std::invalid_argument("pairs without their own σ, ε are not handled");
    topo.pairs.push_back({u(0), u(1), P(2) / kKJ, P(1) * 10});
  } else if (type == "dihedrals") {
    const int f = int(P(0));
    if (f != 1 && f != 9 && f != 4) throw std::invalid_argument("dihedral function " + params[0] + " is not handled");
    topo.dihedrals.push_back({u(0), u(1), u(2), u(3), f, P(2) / kKJ, P(1) * kDeg, int(P(3)), group("dihedral")});
  } else if (type == "impropers") {
    topo.dihedrals.push_back({u(0), u(1), u(2), u(3), 2, P(2) / (2 * kKJ), P(1) * kDeg, 0, group("improper")});
  } else if (type == "exclusions") {
    for (size_t k = 1; k < atoms.size(); ++k) topo.exclusions.push_back({u(0), u(k)});
  } else if (type == "virtual_sitesn") {
    ExplicitTopology::VSite v;
    v.site = u(0);
    for (size_t k = 1; k < atoms.size(); ++k) v.from.push_back(u(k));
    const int f = int(P(0));
    if (f == 1) v.w.assign(v.from.size(), 1.0 / double(v.from.size()));   // centre of geometry
    else if (f != 2) throw std::invalid_argument("virtual_sitesn function " + params[0] + " is not handled");
    topo.vsites.push_back(v);   // function 2: the centre of mass
  } else if (type == "virtual_sites2" || type == "virtual_sites3") {
    const bool three = type == "virtual_sites3";
    if (int(P(0)) != 1) throw std::invalid_argument(type + " function " + params[0] + " is not handled (only linear ones)");
    ExplicitTopology::VSite v;
    v.site = u(0);
    if (three) v.from = {u(1), u(2), u(3)}, v.w = {1 - P(1) - P(2), P(1), P(2)};
    else v.from = {u(1), u(2)}, v.w = {1 - P(1), P(1)};
    topo.vsites.push_back(v);
  } else {
    throw std::invalid_argument("interactions of type " + type + " are not handled");
  }
}



const Json& small_model(const std::string& path) {
  static std::map<std::string, std::unique_ptr<Json>> cache;
  auto& slot = cache[path];
  if (!slot) {
    std::ifstream in(path);
    if (!in) throw std::invalid_argument("cannot open " + path);
    std::stringstream ss;
    ss << in.rdbuf();
    slot = std::make_unique<Json>(Json::parse(ss.str()));
    if (slot->text("format") != "caps-martini3-small-molecules") throw std::invalid_argument(path + ": not a Martini 3 small-molecule file");
  }
  return *slot;
}

double size_mass(const std::string& type) { return type[0] == 'T' ? 36.0 : type[0] == 'S' ? 54.0 : 72.0; }

// append b's atoms, bonds, types and explicit topology to a (b's molecules numbered after a's)
void append_system(System& a, const System& b) {
  const uint32_t off = uint32_t(a.atoms.size());
  int64_t mol0 = 0;
  for (const auto& x : a.atoms) mol0 = std::max(mol0, x.mol);
  auto ta = a.topology ? std::make_shared<ExplicitTopology>(*a.topology) : std::make_shared<ExplicitTopology>();
  const ExplicitTopology empty;
  const ExplicitTopology& tb = b.topology ? *b.topology : empty;
  if (!ta->masses.empty() || !tb.masses.empty()) {
    ta->masses.resize(off, std::nan(""));
    for (size_t i = 0; i < b.atoms.size(); ++i) ta->masses.push_back(i < tb.masses.size() ? tb.masses[i] : std::nan(""));
  }
  std::map<std::string, int> tid;
  for (const auto& t : a.types) tid[t.label] = t.type;
  for (auto at : b.atoms) {
    at.id = int64_t(a.atoms.size() + 1);
    at.mol += mol0;
    auto [it, fresh] = tid.emplace(at.name, int(tid.size()) + 1);
    at.type = it->second;
    if (fresh) a.types.push_back({it->second, size_mass(at.name), at.name});
    a.atoms.push_back(at);
  }
  for (auto x : b.bonds) a.bonds.push_back({x.i + off, x.j + off, x.order});
  for (auto x : tb.bonds) x.i += off, x.j += off, ta->bonds.push_back(x);
  for (auto x : tb.angles) x.i += off, x.j += off, x.k += off, ta->angles.push_back(x);
  for (auto x : tb.dihedrals) x.i += off, x.j += off, x.k += off, x.l += off, ta->dihedrals.push_back(x);
  for (auto x : tb.exclusions) ta->exclusions.push_back({x.first + off, x.second + off});
  for (auto x : tb.pairs) x.i += off, x.j += off, ta->pairs.push_back(x);
  for (auto v : tb.vsites) {
    v.site += off;
    for (auto& f : v.from) f += off;
    ta->vsites.push_back(v);
  }
  ta->natoms = a.atoms.size();
  ta->source = ta->source.empty() ? tb.source : (tb.source.empty() || tb.source == ta->source ? ta->source : ta->source + " and " + tb.source);
  a.topology = ta;
  a.has_charges = true;
  a.has_mol = true;
  a.bonds_from_file = true;
}

// one isomorphism of a reference residue's heavy atoms onto a molecule's (elements, heavy-atom bonds; hydrogen counts
// too when the molecule has hydrogens): ref heavy index → molecule atom, or empty
std::vector<int> match_heavy(const std::vector<int>& rz, const std::vector<std::vector<int>>& radj, const std::vector<int>& rh,
                             const std::vector<int>& mz, const std::vector<std::vector<int>>& madj, const std::vector<int>& mh, bool use_h) {
  const size_t n = rz.size();
  if (n != mz.size()) return {};
  // search order: breadth first from the atom of the rarest element, so each later atom has a matched neighbour
  std::vector<size_t> order;
  std::vector<int> anchor(n, -1);
  std::vector<char> seen(n, 0);
  std::map<int, int> count;
  for (int z : rz) ++count[z];
  std::vector<size_t> starts(n);
  std::iota(starts.begin(), starts.end(), 0);
  std::stable_sort(starts.begin(), starts.end(), [&](size_t a, size_t b) { return count[rz[a]] < count[rz[b]]; });
  for (size_t s0 : starts) {
    if (seen[s0]) continue;
    std::queue<size_t> q;
    q.push(s0);
    seen[s0] = 1;
    while (!q.empty()) {
      const size_t u = q.front();
      q.pop();
      order.push_back(u);
      for (int v : radj[u])
        if (!seen[size_t(v)]) seen[size_t(v)] = 1, anchor[size_t(v)] = int(u), q.push(size_t(v));
    }
  }
  std::vector<int> img(n, -1), used(n, 0);
  std::vector<std::set<int>> mset(n);
  for (size_t i = 0; i < n; ++i) mset[i] = std::set<int>(madj[i].begin(), madj[i].end());
  std::function<bool(size_t)> rec = [&](size_t d) -> bool {
    if (d == n) return true;
    const size_t r = order[d];
    std::vector<int> cand;
    if (anchor[r] >= 0) cand.assign(madj[size_t(img[size_t(anchor[r])])].begin(), madj[size_t(img[size_t(anchor[r])])].end());
    else {
      cand.resize(n);
      std::iota(cand.begin(), cand.end(), 0);
    }
    for (int c : cand) {
      if (used[size_t(c)] || mz[size_t(c)] != rz[r] || madj[size_t(c)].size() != radj[r].size() || (use_h && mh[size_t(c)] != rh[r])) continue;
      bool ok = true;
      for (size_t e = 0; e < d && ok; ++e) {
        const size_t r2 = order[e];
        const bool re = std::find(radj[r].begin(), radj[r].end(), int(r2)) != radj[r].end();
        ok = re == (mset[size_t(c)].count(img[r2]) > 0);
      }
      if (!ok) continue;
      img[r] = c, used[size_t(c)] = 1;
      if (rec(d + 1)) return true;
      img[r] = -1, used[size_t(c)] = 0;
    }
    return false;
  };
  return rec(0) ? img : std::vector<int>{};
}

}  // namespace

System martini3_small_molecules(const System& aa, const std::string& data_path, double constraint_kj, std::vector<std::string>* unmatched,
                                std::vector<std::string>* notes, bool geometric) {
  const Json& M = small_model(data_path)["molecules"];
  System out;
  out.title = aa.title;
  out.cell = aa.cell;
  out.source_format = "caps-martini3-small-molecules";
  if (aa.atoms.empty()) return out;
  int ncomp = 0;
  System tmp = aa;
  tmp.has_mol = false;
  const auto comp = tmp.molecules(&ncomp);
  const auto nb = aa.neighbours();
  std::vector<std::vector<uint32_t>> members(static_cast<size_t>(ncomp));
  for (uint32_t i = 0; i < aa.atoms.size(); ++i) members[size_t(comp[i])].push_back(i);
  std::map<std::string, int> found;
  for (const auto& atoms : members) {
    // the molecule's heavy atoms, their heavy neighbours and hydrogen counts
    std::vector<uint32_t> heavy;
    std::map<uint32_t, int> hidx;
    for (uint32_t a : atoms)
      if (aa.atoms[a].element > 1) hidx[a] = int(heavy.size()), heavy.push_back(a);
    std::vector<int> mz, mh;
    std::vector<std::vector<int>> madj(heavy.size());
    bool has_h = false;
    for (size_t k = 0; k < heavy.size(); ++k) {
      mz.push_back(aa.atoms[heavy[k]].element);
      int h = 0;
      for (uint32_t v : nb[heavy[k]]) {
        if (aa.atoms[v].element == 1) ++h;
        else if (hidx.count(v)) madj[k].push_back(hidx[v]);
      }
      mh.push_back(h);
      has_h = has_h || h > 0;
    }
    std::string hit;
    std::vector<int> img;
    std::vector<int> rheavy;
    const Json* ref = nullptr;
    for (const auto& [name, mol] : M.members()) {
      const auto& ra = mol["aa"]["atoms"].items();
      std::vector<int> rz, rh;
      rheavy.clear();
      std::map<int, int> rk;
      for (size_t i = 0; i < ra.size(); ++i) {
        const int z = element_from_symbol(ra[i]["element"].str());
        if (z > 1) rk[int(i)] = int(rz.size()), rz.push_back(z), rheavy.push_back(int(i));
      }
      if (rz.size() != mz.size()) continue;
      {
        std::vector<int> a = rz, b = mz;
        std::sort(a.begin(), a.end()), std::sort(b.begin(), b.end());
        if (a != b) continue;
      }
      std::vector<std::vector<int>> radj(rz.size());
      rh.assign(rz.size(), 0);
      for (const auto& bd : mol["aa"]["bonds"].items()) {
        const int i = int(bd[0].number()), j = int(bd[1].number());
        const bool hi = !rk.count(i), hj = !rk.count(j);
        if (!hi && !hj) radj[size_t(rk[i])].push_back(rk[j]), radj[size_t(rk[j])].push_back(rk[i]);
        else if (hi && !hj) ++rh[size_t(rk[j])];
        else if (hj && !hi) ++rh[size_t(rk[i])];
      }
      img = match_heavy(rz, radj, rh, mz, madj, mh, has_h);
      if (!img.empty()) { hit = name; ref = &mol; break; }
    }
    if (!ref) {
      if (unmatched) {   // its formula, Hill order
        std::map<std::string, int> f;
        for (uint32_t a : atoms) ++f[element(aa.atoms[a].element).symbol];
        std::string formula;
        auto put = [&](const std::string& e) { if (f.count(e)) formula += e + (f[e] > 1 ? std::to_string(f[e]) : ""), f.erase(e); };
        if (f.count("C")) put("C"), put("H");
        for (const auto& [e, c] : f) formula += e + (c > 1 ? std::to_string(c) : "");
        unmatched->push_back(formula);
      }
      continue;
    }
    ++found[hit];
    // reference atom index → the molecule's atom: heavy atoms by the isomorphism, hydrogens in order on their atom
    const auto& ra = (*ref)["aa"]["atoms"].items();
    std::vector<int> to(ra.size(), -1);
    for (size_t k = 0; k < rheavy.size(); ++k) to[size_t(rheavy[k])] = int(heavy[size_t(img[k])]);
    {
      std::vector<std::vector<int>> rhs(ra.size());
      for (const auto& bd : (*ref)["aa"]["bonds"].items()) {
        const int i = int(bd[0].number()), j = int(bd[1].number());
        if (ra[size_t(i)]["element"].str() == "H" && ra[size_t(j)]["element"].str() != "H") rhs[size_t(j)].push_back(i);
        if (ra[size_t(j)]["element"].str() == "H" && ra[size_t(i)]["element"].str() != "H") rhs[size_t(i)].push_back(j);
      }
      for (int rhv : rheavy) {
        std::vector<uint32_t> mhs;
        for (uint32_t v : nb[size_t(to[size_t(rhv)])])
          if (aa.atoms[v].element == 1) mhs.push_back(v);
        for (size_t q = 0; q < rhs[size_t(rhv)].size() && q < mhs.size(); ++q) to[size_t(rhs[size_t(rhv)][q])] = int(mhs[q]);
      }
    }
    // beads at vermouth's mass-weighted centres of their mapped atoms
    const Json& blk = (*ref)["block"];
    std::map<std::string, std::pair<Vec3, double>> acc;
    const Vec3 origin = aa.atoms[size_t(to[size_t(rheavy[0])])].pos;
    for (const auto& e : (*ref)["map"].items()) {
      const int ai = to[size_t(e[0].number())];
      if (ai < 0) continue;
      const Vec3 d = aa.cell.valid() ? aa.cell.minimum_image(aa.atoms[size_t(ai)].pos - origin) : aa.atoms[size_t(ai)].pos - origin;
      for (const auto& bw : e[1].items()) {
        const double wm = bw[1].number() * (geometric ? 1.0 : vermouth_mass(aa.atoms[size_t(ai)].element));
        auto& s = acc[bw[0].str()];
        s.first = s.first + (origin + d) * wm;
        s.second += wm;
      }
    }
    System one;
    auto topo = std::make_shared<ExplicitTopology>();
    topo->source = "Martini 3 small molecules";
    std::map<std::string, int> idx;
    bool any_mass = false;
    for (const auto& at : blk["atoms"].items()) {
      Atom a;
      a.id = int64_t(one.atoms.size() + 1);
      a.mol = 1;
      a.name = at["atype"].str();
      a.resname = hit;
      a.resid = 1;
      a.charge = at.has("charge") ? to_num(at["charge"]) : 0.0;
      auto it = acc.find(at["atomname"].str());
      if (it != acc.end() && std::fabs(it->second.second) > 1e-7) a.pos = it->second.first * (1 / it->second.second);
      idx[at["atomname"].str()] = int(one.atoms.size());
      one.atoms.push_back(a);
      topo->masses.push_back(at.has("mass") ? to_num(at["mass"]) : std::nan(""));
      any_mass = any_mass || at.has("mass");
    }
    for (const auto& [type, list] : blk["interactions"].members())
      for (const auto& it : list.items()) {
        std::vector<int> ids;
        for (const auto& x : it["atoms"].items()) ids.push_back(idx.at(x.str()));
        std::vector<std::string> ps;
        for (const auto& x : it["params"].items()) ps.push_back(x.str());
        add_term(*topo, one, type, ids, ps, it["meta"], constraint_kj);
      }
    // virtual sites where their beads put them
    for (const auto& v : topo->vsites) {
      Vec3 c{0, 0, 0};
      double wt = 0;
      for (size_t k = 0; k < v.from.size(); ++k) {
        const double mk = !std::isnan(topo->masses[v.from[k]]) ? topo->masses[v.from[k]] : size_mass(one.atoms[v.from[k]].name);
        const double w = v.w.empty() ? mk : v.w[k];
        c = c + one.atoms[v.from[k]].pos * w;
        wt += w;
      }
      if (wt != 0) one.atoms[v.site].pos = c * (1 / wt);
    }
    if (!any_mass) topo->masses.clear();
    topo->natoms = one.atoms.size();
    one.topology = topo;
    append_system(out, one);
  }
  if (notes && !found.empty()) {
    std::string l;
    for (const auto& [n, c] : found) l += (l.empty() ? "" : ", ") + std::to_string(c) + " " + n;
    notes->push_back("small molecules mapped onto Martini 3 beads: " + l);
  }
  return out;
}

System martini3_all_atom(const System& aa_in, const std::string& small_path, double constraint_kj, std::vector<std::string>* notes_out,
                         bool geometric, bool refuse_unmatched) {
  std::vector<std::string> notes;
  // bonds from the file and from distances
  System aa = aa_in;
  {
    std::set<std::pair<uint32_t, uint32_t>> have;
    for (const auto& b : aa.bonds) have.insert({std::min(b.i, b.j), std::max(b.i, b.j)});
    for (const auto& b : perceive_bonds(aa))
      if (have.insert({std::min(b.i, b.j), std::max(b.i, b.j)}).second) aa.bonds.push_back(b);
  }
  // water and ions are not mapped (a W bead is four waters); the other molecules are Martini 3's small molecules by graph
  System rest;
  rest.cell = aa.cell;
  std::vector<int> keep(aa.atoms.size(), -1);
  int water = 0;
  static const std::set<std::string> waters = {"HOH", "WAT", "SOL", "TIP3", "TIP4", "SPC", "H2O", "T3P", "T4P"};
  const auto anb = aa.neighbours();
  std::vector<char> is_water(aa.atoms.size(), 0);   // an O with two H and nothing else, and those H
  for (size_t i = 0; i < aa.atoms.size(); ++i) {
    if (aa.atoms[i].element != 8) continue;
    int h = 0, x = 0;
    for (uint32_t v : anb[i]) (aa.atoms[v].element == 1 ? h : x)++;
    if ((h == 2 && x == 0) || (waters.count(aa.atoms[i].resname) && h + x <= 2)) {
      is_water[i] = 1;
      for (uint32_t v : anb[i]) is_water[v] = 1;
    }
  }
  int ions = 0;
  for (size_t i = 0; i < aa.atoms.size(); ++i) {
    const auto& a = aa.atoms[i];
    if (is_water[i] || waters.count(a.resname)) {
      water += a.element == 8;
      continue;
    }
    if (anb[i].empty() && a.element > 1) {   // a lone atom: an ion (Martini's ions are hydrated beads, not mapped)
      ++ions;
      continue;
    }
    keep[i] = int(rest.atoms.size());
    rest.atoms.push_back(a);
  }
  for (const auto& b : aa.bonds)
    if (keep[b.i] >= 0 && keep[b.j] >= 0) rest.bonds.push_back({uint32_t(keep[b.i]), uint32_t(keep[b.j]), b.order});
  System out;
  if (!rest.atoms.empty()) {
    std::vector<std::string> un;
    out = martini3_small_molecules(rest, small_path, constraint_kj, &un, &notes, geometric);
    if (!un.empty()) {
      std::map<std::string, int> c;
      for (const auto& f : un) ++c[f];
      std::string l;
      for (const auto& [f, k] : c) l += (l.empty() ? "" : ", ") + std::to_string(k) + " × " + f;
      // typing refuses rather than drop part of the structure; an explicit mapping (caps martini) leaves them out
      if (refuse_unmatched)
        throw std::invalid_argument(std::to_string(un.size()) + " molecules match no Martini 3 small molecule (" + l +
                                    "): Martini 3 has no beads for them (build them from its templates or bead SMILES)");
      notes.push_back(std::to_string(un.size()) + " molecules match no Martini 3 small molecule and are left out (" + l + ")");
    }
  }
  if (ions) notes.push_back(std::to_string(ions) + " ions left out (Martini 3's ions are hydrated beads: add them to the coarse-grained structure)");
  if (water) notes.push_back(std::to_string(water) + " water molecules left out (a Martini W bead is four waters: solvate the coarse-grained structure)");
  if (out.atoms.empty()) throw std::invalid_argument("no Martini 3 small molecules in the structure");
  out.title = aa.title;
  out.cell = aa.cell;
  if (notes_out) *notes_out = notes;
  return out;
}

System gromacs_molecule(const Json& mol, const std::function<double(const std::string&)>& type_mass, double constraint_kj, uint64_t seed) {
  BeadMolecule m;
  std::vector<double> len;
  const auto& atoms = mol["atoms"].items();
  for (const auto& a : atoms) {
    m.type.push_back(a["type"].str());
    m.charge.push_back(a["charge"].number());
  }
  auto topo = std::make_shared<ExplicitTopology>();
  topo->source = "Martini 3 molecule";
  const size_t n = atoms.size();
  topo->masses.assign(n, std::nan(""));
  bool any_mass = false;
  for (size_t i = 0; i < n; ++i)
    if (atoms[i].has("mass") && !atoms[i]["mass"].is_null()) topo->masses[i] = atoms[i]["mass"].number(), any_mass = true;
  if (!any_mass) topo->masses.clear();
  auto u = [](const Json& v) { return uint32_t(v.number()); };
  for (const auto& b : mol["bonds"].items()) {
    if (int(b[2].number()) != 1) throw std::invalid_argument("bond function " + std::to_string(int(b[2].number())) + " is not handled");
    topo->bonds.push_back({u(b[0]), u(b[1]), b[4].number() / (2 * kKJ * 100), b[3].number() * 10, "bond"});
    m.bonds.push_back({int(b[0].number()), int(b[1].number())});
    len.push_back(b[3].number() * 10);
  }
  for (const auto& c : mol["constraints"].items()) {
    topo->bonds.push_back({u(c[0]), u(c[1]), constraint_kj / (2 * kKJ * 100), c[3].number() * 10, "constraint"});
    m.bonds.push_back({int(c[0].number()), int(c[1].number())});
    len.push_back(c[3].number() * 10);
  }
  for (const auto& a : mol["angles"].items()) {
    const int f = int(a[3].number());
    const double th = a[4].number() * kDeg, k = a[5].number();
    if (f == 1) topo->angles.push_back({u(a[0]), u(a[1]), u(a[2]), 0, k / (2 * kKJ), th, "angle"});
    else if (f == 2) topo->angles.push_back({u(a[0]), u(a[1]), u(a[2]), 1, k / (2 * kKJ), th, "angle"});
    else if (f == 10) topo->angles.push_back({u(a[0]), u(a[1]), u(a[2]), 5, k / (2 * kKJ), th, "angle"});
    else throw std::invalid_argument("angle function " + std::to_string(f) + " is not handled");
  }
  for (const auto& d : mol["dihedrals"].items()) {
    const int f = int(d[4].number());
    const uint32_t a = u(d[0]), b = u(d[1]), c = u(d[2]), e = u(d[3]);
    if (f == 1 || f == 9 || f == 4) {
      const int mult = d.size() > 7 ? int(d[7].number()) : 1;
      topo->dihedrals.push_back({a, b, c, e, f == 4 ? 4 : 1, d[6].number() / kKJ, d[5].number() * kDeg, mult, "dihedral"});
    } else if (f == 2) {
      topo->dihedrals.push_back({a, b, c, e, 2, d[6].number() / (2 * kKJ), d[5].number() * kDeg, 0, "improper"});
    } else {
      throw std::invalid_argument("dihedral function " + std::to_string(f) + " is not handled");
    }
  }
  for (const auto& x : mol["exclusions"].items())
    for (size_t k = 1; k < x.size(); ++k) topo->exclusions.push_back({u(x[0]), u(x[k])});
  for (const auto& v : mol["vsites"].items()) {
    ExplicitTopology::VSite vs;
    vs.site = u(v["site"]);
    for (const auto& a : v["from"].items()) vs.from.push_back(u(a));
    if (v.has("weights") && !v["weights"].is_null())
      for (const auto& w : v["weights"].items()) vs.w.push_back(w.number());
    topo->vsites.push_back(vs);
  }
  BeadBuildOptions o;
  o.seed = seed;
  o.mass = type_mass;
  System s = build_bead_graph(m, len, o);
  // virtual sites where their atoms put them (weights, or the centre of mass)
  for (const auto& v : topo->vsites) {
    Vec3 c{0, 0, 0};
    double wt = 0;
    for (size_t k = 0; k < v.from.size(); ++k) {
      const double w = !v.w.empty() ? v.w[k] : (!topo->masses.empty() && !std::isnan(topo->masses[v.from[k]]) ? topo->masses[v.from[k]] : type_mass(m.type[v.from[k]]));
      c = c + s.atoms[v.from[k]].pos * w;
      wt += w;
    }
    if (wt != 0) s.atoms[v.site].pos = c * (1 / wt);
  }
  topo->natoms = n;
  s.topology = topo;
  return s;
}

std::string martini3_itp(const System& s, double constraint_kj) {
  if (!s.topology) throw std::invalid_argument("the structure has no explicit topology");
  const ExplicitTopology& T = *s.topology;
  const double kc = constraint_kj / (2 * kKJ * 100);
  std::ostringstream o;
  char b[256];
  o << "[ moleculetype ]\nmolecule 1\n\n[ atoms ]\n";
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
