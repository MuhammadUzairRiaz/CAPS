// Building blocks (GROMOS's residue topologies with their termini): residues recognised in a structure by their atoms
// and bonds take the block's types, charges and bonded terms, as pdb2gmx applies aminoacids.rtp and the .tdb patches.
#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include "caps/elements.hpp"
#include "caps/ffdef.hpp"
#include "caps/json.hpp"

namespace caps {

namespace {

constexpr double kKJ = 4.184;
constexpr double kDeg = 3.14159265358979323846 / 180;

struct TAtom {
  std::string name, type;
  double q = 0;
  int z = 0;
  std::vector<int> nb;       // bonded atoms within the block
  std::vector<int> links;    // elements of the atoms of neighbouring residues it bonds to (-C, +N ...)
};

struct Block {
  std::string key, resname;
  std::vector<TAtom> atoms;
  std::map<std::string, int> idx;
  const Json* j = nullptr;
  int deg(int t) const { return int(atoms[size_t(t)].nb.size() + atoms[size_t(t)].links.size()); }
};

int link_element(const std::string& name) {   // "+N" -> N, "-O3'" -> O, "-C" -> C
  const std::string n = name.substr(1);
  for (const char* two : {"Cl", "CL", "Fe", "FE", "Zn", "ZN", "Mg", "MG", "Ca", "CA"})
    if (n == two) return element_from_symbol(two);
  return element_from_symbol(n.substr(0, 1));
}

}  // namespace

std::string apply_residue_templates(System& s, const FFDef& ff) {
  const Json& R = *ff.residue_templates;
  const Json& codes = R["codes"];
  std::map<std::string, int> zof;
  for (const auto& t : ff.types) zof[t.name] = t.element;
  // the blocks, larger first (a terminal variant before the residue it patches)
  std::vector<Block> blocks;
  for (const auto& [key, rj] : R["residues"].members()) {
    Block b;
    b.key = key;
    b.resname = key.substr(0, key.find(':'));
    b.j = &rj;
    for (const auto& a : rj["atoms"].items()) {
      TAtom t;
      t.name = a[0].str();
      t.type = a[1].str();
      t.q = a[2].number();
      auto it = zof.find(t.type);
      if (it == zof.end()) throw FFError(ff.name + ": block " + key + " names type " + t.type + ", not in the force field");
      t.z = it->second;
      b.idx[t.name] = int(b.atoms.size());
      b.atoms.push_back(t);
    }
    for (const auto& bd : rj["bonds"].items()) {
      const std::string x = bd[0].str(), y = bd[1].str();
      const bool lx = x[0] == '+' || x[0] == '-', ly = y[0] == '+' || y[0] == '-';
      if (lx && ly) continue;
      if (lx || ly) {
        const std::string own = lx ? y : x, other = lx ? x : y;
        if (b.idx.count(own)) b.atoms[size_t(b.idx[own])].links.push_back(link_element(other));
        continue;
      }
      if (!b.idx.count(x) || !b.idx.count(y)) throw FFError(ff.name + ": block " + key + " bonds an atom it lacks (" + x + "-" + y + ")");
      b.atoms[size_t(b.idx[x])].nb.push_back(b.idx[y]);
      b.atoms[size_t(b.idx[y])].nb.push_back(b.idx[x]);
    }
    // a "+N" written in the residue before means a "-C" link here: every link atom counts it from its own block
    blocks.push_back(std::move(b));
  }
  // the link a block's first atom takes from the residue before it is listed there ("C +N"): give each block's own atom
  // that link too, from the residue's own terms that reach back ("-C N ..."): an atom bonded across by an angle's middle
  for (auto& b : blocks) {
    std::set<std::pair<int, std::string>> have;
    for (const auto& bd : (*b.j)["bonds"].items()) {
      const std::string x = bd[0].str(), y = bd[1].str();
      if (x[0] == '+' || x[0] == '-') have.insert({b.idx.count(y) ? b.idx[y] : -1, x});
      if (y[0] == '+' || y[0] == '-') have.insert({b.idx.count(x) ? b.idx[x] : -1, y});
    }
    for (const auto& an : (*b.j)["angles"].items()) {
      // "-C N CA": N (the middle) bonds to -C
      const std::string a0 = an[0].str(), a1 = an[1].str(), a2 = an[2].str();
      for (const auto& [end, mid] : {std::pair{a0, a1}, std::pair{a2, a1}})
        if ((end[0] == '-' || end[0] == '+') && b.idx.count(mid) && !have.count({b.idx[mid], end})) {
          // only when no other atom of the block takes that link (CA -C N: N's -C)
          bool other = false;
          for (const auto& [ia, nm] : have) other = other || nm == end;
          if (other) continue;
          have.insert({b.idx[mid], end});
          b.atoms[size_t(b.idx[mid])].links.push_back(link_element(end));
        }
    }
  }
  std::stable_sort(blocks.begin(), blocks.end(), [](const Block& a, const Block& b) { return a.atoms.size() > b.atoms.size(); });

  const size_t n = s.atoms.size();
  std::vector<std::vector<uint32_t>> nb(n);
  for (const auto& b : s.bonds) nb[b.i].push_back(b.j), nb[b.j].push_back(b.i);
  std::vector<int> owner(n, -1), slot(n, -1);   // residue index, block atom index
  struct Placed { const Block* b; std::vector<uint32_t> at; };
  std::vector<Placed> placed;

  for (uint32_t a0 = 0; a0 < n; ++a0) {
    if (owner[a0] >= 0) continue;
    // the residue holding a0: the first block (largest first) matching with some block atom on a0
    bool found = false;
    for (const auto& B : blocks) {
      if (found) break;
      for (size_t t0 = 0; t0 < B.atoms.size() && !found; ++t0) {
        if (B.atoms[t0].z != s.atoms[a0].element || B.deg(int(t0)) != int(nb[a0].size())) continue;
        // breadth-first order of the block from t0 (a block is connected)
        std::vector<int> order{int(t0)}, parent{-1};
        std::vector<char> seen(B.atoms.size(), 0);
        seen[t0] = 1;
        for (size_t k = 0; k < order.size(); ++k)
          for (int u : B.atoms[size_t(order[k])].nb)
            if (!seen[size_t(u)]) seen[size_t(u)] = 1, order.push_back(u), parent.push_back(order[k]);
        if (order.size() != B.atoms.size()) continue;
        std::vector<int64_t> map(B.atoms.size(), -1);
        std::set<uint32_t> used;
        std::function<bool(size_t)> extend = [&](size_t k) -> bool {
          if (k == order.size()) return true;
          const int t = order[k];
          std::vector<uint32_t> cand;
          if (k == 0) cand = {a0};
          else
            for (uint32_t x : nb[size_t(map[size_t(parent[k])])]) cand.push_back(x);
          // atoms already named as the block names them first (a structure from pdb2gmx, a PDB)
          std::stable_sort(cand.begin(), cand.end(), [&](uint32_t x, uint32_t y) { return (s.atoms[x].name == B.atoms[size_t(t)].name) > (s.atoms[y].name == B.atoms[size_t(t)].name); });
          for (uint32_t x : cand) {
            if (owner[x] >= 0 || used.count(x) || s.atoms[x].element != B.atoms[size_t(t)].z || int(nb[x].size()) != B.deg(t)) continue;
            bool ok = true;
            for (size_t m = 0; m < k && ok; ++m) {   // bonded in the block exactly where bonded in the structure
              const int u = order[m];
              const bool inb = std::find(B.atoms[size_t(t)].nb.begin(), B.atoms[size_t(t)].nb.end(), u) != B.atoms[size_t(t)].nb.end();
              const bool ins = std::find(nb[x].begin(), nb[x].end(), uint32_t(map[size_t(u)])) != nb[x].end();
              ok = inb == ins;
            }
            if (!ok) continue;
            map[size_t(t)] = x;
            used.insert(x);
            if (extend(k + 1)) return true;
            used.erase(x);
            map[size_t(t)] = -1;
          }
          return false;
        };
        if (!extend(0)) continue;
        // the bonds leaving the residue: exactly the block's links, to atoms of the elements they name
        bool ok = true;
        for (size_t t = 0; t < B.atoms.size() && ok; ++t) {
          std::vector<int> out;
          for (uint32_t x : nb[size_t(map[t])])
            if (!used.count(x)) out.push_back(s.atoms[x].element);
          std::vector<int> want = B.atoms[t].links;
          std::sort(out.begin(), out.end());
          std::sort(want.begin(), want.end());
          ok = out == want;
        }
        if (!ok) continue;
        Placed p{&B, {}};
        for (size_t t = 0; t < B.atoms.size(); ++t) {
          p.at.push_back(uint32_t(map[t]));
          owner[size_t(map[t])] = int(placed.size());
          slot[size_t(map[t])] = int(t);
        }
        placed.push_back(std::move(p));
        found = true;
      }
    }
  }
  // rigid solvent models (water, methanol, DMSO, chloroform, CCl4): their blocks list the constraints that hold the
  // geometry as bonds (water's H–H, chloroform's H–Cl but no C–H), so a whole molecule matches such a block when the
  // elements agree and every chemical bond is at most two block bonds apart; the molecule then takes the block's bonds,
  // as GROMACS reads them
  std::vector<std::pair<uint32_t, uint32_t>> rigid_bonds;
  {
    // molecules by their bonds (a structure's molecule ids may group several)
    std::vector<int> comp(n, -1);
    std::vector<std::vector<uint32_t>> members;
    for (uint32_t i = 0; i < n; ++i) {
      if (comp[i] >= 0) continue;
      members.emplace_back();
      std::vector<uint32_t> todo{i};
      comp[i] = int(members.size() - 1);
      while (!todo.empty()) {
        const uint32_t x = todo.back();
        todo.pop_back();
        members.back().push_back(x);
        for (uint32_t y : nb[x])
          if (comp[y] < 0) comp[y] = comp[i], todo.push_back(y);
      }
      std::sort(members.back().begin(), members.back().end());
    }
    for (const auto& m : members) {
      if (m.empty() || m.size() > 8 || std::any_of(m.begin(), m.end(), [&](uint32_t i) { return owner[i] >= 0; })) continue;
      for (const auto& B : blocks) {
        if (B.atoms.size() != m.size() || std::any_of(B.atoms.begin(), B.atoms.end(), [](const TAtom& t) { return !t.links.empty(); })) continue;
        const size_t k = m.size();
        std::vector<std::vector<int>> d(k, std::vector<int>(k, 99));   // block-graph distances
        for (size_t a = 0; a < k; ++a) {
          d[a][a] = 0;
          for (int b : B.atoms[a].nb) d[a][size_t(b)] = 1;
        }
        for (size_t c = 0; c < k; ++c)
          for (size_t a = 0; a < k; ++a)
            for (size_t b = 0; b < k; ++b) d[a][b] = std::min(d[a][b], d[a][c] + d[c][b]);
        std::vector<int> perm(k);
        for (size_t a = 0; a < k; ++a) perm[a] = int(a);
        bool ok = false;
        do {
          ok = true;
          for (size_t a = 0; a < k && ok; ++a) ok = B.atoms[size_t(perm[a])].z == s.atoms[m[a]].element;
          for (size_t a = 0; a < k && ok; ++a)
            for (uint32_t y : nb[m[a]]) {
              const size_t b = size_t(std::find(m.begin(), m.end(), y) - m.begin());
              ok = ok && d[size_t(perm[a])][size_t(perm[b])] <= 2;
            }
        } while (!ok && std::next_permutation(perm.begin(), perm.end()));
        if (!ok) continue;
        Placed p{&B, std::vector<uint32_t>(k)};
        for (size_t a = 0; a < k; ++a) p.at[size_t(perm[a])] = m[a];
        for (size_t t = 0; t < k; ++t) owner[p.at[t]] = int(placed.size()), slot[p.at[t]] = int(t);
        for (size_t t = 0; t < k; ++t)
          for (int u : B.atoms[t].nb)
            if (size_t(u) > t) rigid_bonds.push_back({p.at[t], p.at[size_t(u)]});
        placed.push_back(std::move(p));
        break;
      }
    }
    if (!rigid_bonds.empty()) {   // the block's bonds in place of the molecule's
      std::set<uint32_t> in;
      for (const auto& [a, b] : rigid_bonds) in.insert(a), in.insert(b);
      std::vector<Bond> kept;
      for (const auto& b : s.bonds)
        if (!in.count(b.i) && !in.count(b.j)) kept.push_back(b);
      for (const auto& [a, b] : rigid_bonds) kept.push_back({a, b, 1});
      s.bonds = std::move(kept);
      s.bonds_from_file = true;
      for (auto& v : nb) v.clear();
      for (const auto& b : s.bonds) nb[b.i].push_back(b.j), nb[b.j].push_back(b.i);
    }
  }
  std::vector<uint32_t> loose;
  for (uint32_t i = 0; i < n; ++i)
    if (owner[i] < 0) loose.push_back(i);
  if (!loose.empty()) {
    std::string l;
    for (size_t k = 0; k < loose.size() && k < 8; ++k)
      l += (k ? ", " : "") + std::to_string(loose[k] + 1) + " (" + element(s.atoms[loose[k]].element).symbol + ")";
    throw FFError(ff.name + " is a building-block force field: " + std::to_string(loose.size()) + " atoms are in no residue it defines (" + l +
                  (loose.size() > 8 ? ", …" : "") + "); its blocks are the residues of its topology database (amino acids, nucleotides, "
                  "sugars, lipids, solvents)");
  }
  // neighbours: the residue across each '-' and '+' link
  std::vector<int> prev(placed.size(), -1), next(placed.size(), -1);
  for (size_t r = 0; r < placed.size(); ++r)
    for (uint32_t x : placed[r].at)
      for (uint32_t y : nb[x])
        if (owner[y] != int(r)) {
          // which way: the block's own terms say whether an atom of this element is reached by - or +
          const Block& B = *placed[r].b;
          const std::string& own = B.atoms[size_t(slot[x])].name;
          char dir = 0;
          for (const char* kind : {"bonds", "angles", "dihedrals", "impropers"})
            for (const auto& t : (*B.j)[kind].items()) {
              bool has_own = false;
              char d = 0;
              for (size_t q = 0; q + 1 < t.size(); ++q) {
                const std::string nm = t[q].str();
                has_own = has_own || nm == own;
                if ((nm[0] == '-' || nm[0] == '+') && link_element(nm) == s.atoms[y].element) d = nm[0];
              }
              if (has_own && d && std::string(kind) == "bonds") dir = d;
            }
          if (!dir) {   // the link listed by the other residue ("C +N" there): the other way round
            const Block& O = *placed[size_t(owner[y])].b;
            const std::string& theirs = O.atoms[size_t(slot[y])].name;
            for (const auto& t : (*O.j)["bonds"].items()) {
              const std::string p0 = t[0].str(), p1 = t[1].str();
              if ((p0 == theirs && (p1[0] == '+' || p1[0] == '-')) || (p1 == theirs && (p0[0] == '+' || p0[0] == '-')))
                dir = (p0 == theirs ? p1[0] : p0[0]) == '+' ? '-' : '+';
            }
          }
          if (dir == '-') prev[r] = owner[y];
          else if (dir == '+') next[r] = owner[y];
        }
  auto resolve = [&](size_t r, const std::string& nm) -> int64_t {
    size_t rr = r;
    std::string name = nm;
    if (nm[0] == '-' || nm[0] == '+') {
      const int o = nm[0] == '-' ? prev[r] : next[r];
      if (o < 0) return -1;
      rr = size_t(o);
      name = nm.substr(1);
    }
    const auto it = placed[rr].b->idx.find(name);
    return it == placed[rr].b->idx.end() ? -1 : int64_t(placed[rr].at[size_t(it->second)]);
  };
  auto topo = std::make_shared<ExplicitTopology>();
  topo->natoms = n;
  topo->types.resize(n);
  topo->source = ff.name + " building blocks";
  std::map<std::string, int> count;
  int resid = 0;
  for (size_t r = 0; r < placed.size(); ++r) {
    const Block& B = *placed[r].b;
    ++count[B.key];
    ++resid;
    for (size_t t = 0; t < B.atoms.size(); ++t) {
      Atom& a = s.atoms[placed[r].at[t]];
      a.name = B.atoms[t].name;
      a.resname = B.resname;
      a.resid = resid;
      a.charge = B.atoms[t].q;
      topo->types[placed[r].at[t]] = B.atoms[t].type;
    }
    auto term = [&](const Json& t, size_t k, std::vector<uint32_t>& at, std::string& code) {
      at.clear();
      for (size_t q = 0; q < k; ++q) {
        const int64_t x = resolve(r, t[q].str());
        if (x < 0) return false;
        at.push_back(uint32_t(x));
      }
      code = t[k].str();
      return true;
    };
    std::vector<uint32_t> at;
    std::string code;
    for (const auto& t : (*B.j)["bonds"].items()) {
      if (!term(t, 2, at, code)) continue;   // a link to a residue that is not there (the chain's end)
      const auto& c = codes[code];
      topo->bonds.push_back({at[0], at[1], c[1].number() / 4 / kKJ / 1e4, c[0].number() * 10, code, 2});
    }
    for (const auto& t : (*B.j)["angles"].items()) {
      if (!term(t, 3, at, code)) continue;
      const auto& c = codes[code];
      topo->angles.push_back({at[0], at[1], at[2], 1, c[1].number() / 2 / kKJ, c[0].number() * kDeg, code});
    }
    for (const auto& t : (*B.j)["dihedrals"].items()) {
      if (!term(t, 4, at, code)) continue;
      const auto& c = codes[code];
      topo->dihedrals.push_back({at[0], at[1], at[2], at[3], 1, c[1].number() / kKJ, c[0].number() * kDeg, int(c[2].number()), code});
    }
    for (const auto& t : (*B.j)["impropers"].items()) {
      if (!term(t, 4, at, code)) continue;
      const auto& c = codes[code];
      topo->dihedrals.push_back({at[0], at[1], at[2], at[3], 2, c[1].number() / 2 / kKJ, c[0].number() * kDeg, 0, code});
    }
  }
  // every structure bond is a block's bond
  {
    std::set<std::pair<uint32_t, uint32_t>> tb;
    for (const auto& b : topo->bonds) tb.insert({std::min(b.i, b.j), std::max(b.i, b.j)});
    std::vector<std::string> miss;
    for (const auto& b : s.bonds)
      if (!tb.count({std::min(b.i, b.j), std::max(b.i, b.j)}))
        miss.push_back(s.atoms[b.i].resname + std::to_string(s.atoms[b.i].resid) + ":" + s.atoms[b.i].name + "-" + s.atoms[b.j].resname +
                       std::to_string(s.atoms[b.j].resid) + ":" + s.atoms[b.j].name);
    if (!miss.empty())
      throw FFError(ff.name + ": " + std::to_string(miss.size()) + " bonds between residues no block links (" + miss[0] +
                    (miss.size() > 1 ? ", …" : "") + ")");
  }
  // the 1-4 pairs (three bonds apart; nrexcl 3), less the blocks' own exclusions, which are excluded outright
  const int nrexcl = int(R.num("nrexcl", 3));
  std::set<std::pair<uint32_t, uint32_t>> extra;
  for (size_t r = 0; r < placed.size(); ++r)
    if (placed[r].b->j->has("exclusions"))
      for (const auto& e : (*placed[r].b->j)["exclusions"].items()) {
        const int64_t x = resolve(r, e[0].str()), y = resolve(r, e[1].str());
        if (x >= 0 && y >= 0) extra.insert({uint32_t(std::min(x, y)), uint32_t(std::max(x, y))});
      }
  for (uint32_t i = 0; i < n; ++i) {
    std::map<uint32_t, int> d{{i, 0}};
    std::vector<uint32_t> front{i};
    for (int step = 1; step <= nrexcl; ++step) {
      std::vector<uint32_t> nxt;
      for (uint32_t x : front)
        for (uint32_t y : nb[x])
          if (!d.count(y)) d[y] = step, nxt.push_back(y);
      front = std::move(nxt);
    }
    for (const auto& [j, k] : d) {   // 1-3 listed (an explicit topology has no angles to find them by), 1-4 with their own terms
      if (j <= i) continue;
      if (k == 2) topo->exclusions.push_back({i, j});
      if (k == nrexcl && !extra.count({i, j})) topo->pairs14.push_back({i, j});
    }
  }
  for (const auto& e : extra) topo->exclusions.push_back({e.first, e.second});
  s.topology = topo;
  s.has_charges = true;
  std::string list;
  int shown = 0;
  for (const auto& [k, c] : count) {
    if (shown++ < 12) list += (list.empty() ? "" : ", ") + k + (c > 1 ? " ×" + std::to_string(c) : "");
  }
  if (shown > 12) list += ", …";
  return std::to_string(placed.size()) + " residues recognised as " + ff.name + " building blocks (" + list + "): their types, charges and " +
         std::to_string(topo->bonds.size()) + " bonds, " + std::to_string(topo->angles.size()) + " angles, " +
         std::to_string(topo->dihedrals.size()) + " dihedrals and impropers as the blocks list them; " + std::to_string(topo->pairs14.size()) +
         " 1-4 pairs";
}

}  // namespace caps
