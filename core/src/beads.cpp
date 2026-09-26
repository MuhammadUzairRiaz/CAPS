// Coarse-grained molecules from bead SMILES (the template notation of the MARTINI and SDK sources), and the mapping
// of all-atom structures onto beads by fragment rules (the SDK typing files' "beads").
#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <map>
#include <random>
#include <set>
#include <stdexcept>

#include "caps/resolution.hpp"
#include "caps/typing.hpp"

namespace caps {

BeadMolecule parse_bead_smiles(const std::string& text) {
  BeadMolecule m;
  std::vector<int> stack;
  std::map<int, int> ring;   // ring-closure number → bead waiting for it
  int prev = -1;
  size_t p = 0;
  auto fail = [&](const std::string& why) { throw std::invalid_argument("bead SMILES, position " + std::to_string(p + 1) + ": " + why); };
  while (p < text.size()) {
    const char c = text[p];
    if (std::isspace(static_cast<unsigned char>(c))) { ++p; continue; }
    if (c == '[') {
      const size_t e = text.find(']', p);
      if (e == std::string::npos) fail("a bead has no closing ]");
      std::string tok = text.substr(p + 1, e - p - 1);
      double q = 0;
      // a charge suffix: the last + or − followed by digits only (Qa-1, Qd+1, Qa1-1, Qd+2)
      const size_t sgn = tok.find_last_of("+-");
      if (sgn != std::string::npos && sgn > 0) {
        const std::string d = tok.substr(sgn + 1);
        if (std::all_of(d.begin(), d.end(), [](char x) { return std::isdigit(static_cast<unsigned char>(x)); })) {
          q = d.empty() ? 1.0 : std::stod(d);
          if (tok[sgn] == '-') q = -q;
          tok = tok.substr(0, sgn);
        }
      }
      if (tok.empty()) fail("an empty bead");
      const int id = int(m.type.size());
      m.type.push_back(tok);
      m.charge.push_back(q);
      if (prev >= 0) m.bonds.push_back({prev, id});
      prev = id;
      p = e + 1;
    } else if (c == '(') {
      if (prev < 0) fail("a branch before any bead");
      stack.push_back(prev);
      ++p;
    } else if (c == ')') {
      if (stack.empty()) fail("an unmatched )");
      prev = stack.back();
      stack.pop_back();
      ++p;
    } else if (std::isdigit(static_cast<unsigned char>(c)) || c == '%') {
      int n = 0;
      if (c == '%') {
        if (p + 2 >= text.size() || !std::isdigit(static_cast<unsigned char>(text[p + 1])) || !std::isdigit(static_cast<unsigned char>(text[p + 2])))
          fail("% needs two digits");
        n = (text[p + 1] - '0') * 10 + (text[p + 2] - '0');
        // three-digit closures (%100) as EMC writes past 99
        if (p + 3 < text.size() && std::isdigit(static_cast<unsigned char>(text[p + 3]))) n = n * 10 + (text[p + 3] - '0'), ++p;
        p += 3;
      } else {
        n = c - '0';
        ++p;
      }
      if (prev < 0) fail("a ring closure before any bead");
      if (auto it = ring.find(n); it != ring.end()) {
        m.bonds.push_back({it->second, prev});
        ring.erase(it);
      } else {
        ring[n] = prev;
      }
    } else if (c == '-' || c == '.') {
      if (c == '.') prev = -1;
      ++p;
    } else {
      fail(std::string("unexpected '") + c + "'");
    }
  }
  if (!stack.empty()) throw std::invalid_argument("bead SMILES: an unclosed (");
  if (!ring.empty()) throw std::invalid_argument("bead SMILES: ring closure " + std::to_string(ring.begin()->first) + " is not closed");
  if (m.type.empty()) throw std::invalid_argument("bead SMILES: no beads");
  return m;
}

System build_beads(const std::string& text, const BeadBuildOptions& o) {
  const BeadMolecule m = parse_bead_smiles(text);
  const int n = int(m.type.size());
  std::vector<std::vector<int>> nb(static_cast<size_t>(n));
  for (const auto& [a, b] : m.bonds) nb[size_t(a)].push_back(b), nb[size_t(b)].push_back(a);
  std::map<std::pair<int, int>, double> len;   // each bond's target length
  for (const auto& [a, b] : m.bonds) {
    const double l = o.bond_length ? o.bond_length(m.type[size_t(a)], m.type[size_t(b)]) : 0.0;
    len[{std::min(a, b), std::max(a, b)}] = l > 0 ? l : o.default_bond;
  }
  // a bond longer than another path between its beads (Cooke–Deserno's head-to-tail bond across its lipid) is taken at
  // that path's length: the beads in a straight line
  for (auto& [ab, l] : len) {
    std::vector<double> d(size_t(n), 1e30);
    d[size_t(ab.first)] = 0;
    std::vector<char> done(size_t(n), 0);
    for (int it = 0; it < n; ++it) {
      int u = -1;
      for (int v = 0; v < n; ++v)
        if (!done[size_t(v)] && (u < 0 || d[size_t(v)] < d[size_t(u)])) u = v;
      if (u < 0 || d[size_t(u)] >= 1e30) break;
      done[size_t(u)] = 1;
      for (const auto& [e, le] : len) {
        if (e == ab) continue;
        const int w = e.first == u ? e.second : e.second == u ? e.first : -1;
        if (w >= 0 && d[size_t(u)] + le < d[size_t(w)]) d[size_t(w)] = d[size_t(u)] + le;
      }
    }
    if (d[size_t(ab.second)] < l) l = d[size_t(ab.second)];
  }
  auto b0 = [&](int a, int b) {
    auto it = len.find({std::min(a, b), std::max(a, b)});
    return it != len.end() ? it->second : o.default_bond;
  };
  // graph distances (for the spring model: bonded at b0, others kept apart)
  std::vector<std::vector<int>> gd(size_t(n), std::vector<int>(size_t(n), 1 << 20));
  for (int s = 0; s < n; ++s) {
    auto& d = gd[size_t(s)];
    d[size_t(s)] = 0;
    std::vector<int> q{s};
    for (size_t h = 0; h < q.size(); ++h)
      for (int v : nb[size_t(q[h])])
        if (d[size_t(v)] > d[size_t(q[h])] + 1) d[size_t(v)] = d[size_t(q[h])] + 1, q.push_back(v);
  }
  // placement: breadth first, each bead at its parent's bond length in a random direction away from the others
  std::mt19937_64 rng(o.seed);
  std::normal_distribution<double> nd(0, 1);
  std::vector<Vec3> x(size_t(n), Vec3{0, 0, 0});
  std::vector<char> placed(size_t(n), 0);
  for (int root = 0; root < n; ++root) {
    if (placed[size_t(root)]) continue;
    placed[size_t(root)] = 1;
    std::vector<int> q{root};
    for (size_t h = 0; h < q.size(); ++h) {
      const int u = q[h];
      for (int v : nb[size_t(u)]) {
        if (placed[size_t(v)]) continue;
        Vec3 best{0, 0, 0};
        double bestd = -1;
        for (int tr = 0; tr < 40; ++tr) {
          Vec3 dir{nd(rng), nd(rng), nd(rng)};
          dir = dir * (1 / std::max(1e-9, norm(dir)));
          const Vec3 cand = x[size_t(u)] + dir * b0(u, v);
          double dmin = 1e30;
          for (int w = 0; w < n; ++w)
            if (placed[size_t(w)] && w != u) dmin = std::min(dmin, norm(cand - x[size_t(w)]));
          if (dmin > bestd) bestd = dmin, best = cand;
        }
        x[size_t(v)] = best;
        placed[size_t(v)] = 1;
        q.push_back(v);
      }
    }
  }
  // spring model: bonds at b0; beads two bonds apart at least 1.4 b0 (rings excepted: they are bonded), the rest
  // at least 0.9 of the longer bond; steepest descent with a capped step
  double bmax = o.default_bond;
  for (const auto& [a, b] : m.bonds) bmax = std::max(bmax, b0(a, b));
  for (int it = 0; it < 3000; ++it) {
    std::vector<Vec3> g(size_t(n), Vec3{0, 0, 0});
    for (const auto& [a, b] : m.bonds) {
      const Vec3 d = x[size_t(b)] - x[size_t(a)];
      const double r = std::max(1e-6, norm(d)), k = 2 * (r - b0(a, b)) / r;
      g[size_t(b)] = g[size_t(b)] + d * k;
      g[size_t(a)] = g[size_t(a)] - d * k;
    }
    for (int a = 0; a < n; ++a)
      for (int b = a + 1; b < n; ++b) {
        const int gdist = gd[size_t(a)][size_t(b)];
        if (gdist <= 1) continue;
        const double dmin = gdist == 2 ? 1.4 * std::min(b0(a, b), bmax) : 0.9 * bmax;
        const Vec3 d = x[size_t(b)] - x[size_t(a)];
        const double r = std::max(1e-6, norm(d));
        if (r >= dmin) continue;
        const double k = 2 * (r - dmin) / r;
        g[size_t(b)] = g[size_t(b)] + d * k;
        g[size_t(a)] = g[size_t(a)] - d * k;
      }
    double gmax = 0;
    for (const auto& v : g) gmax = std::max(gmax, norm(v));
    if (gmax < 1e-4) break;
    const double step = std::min(0.1, 0.2 / gmax);
    for (int a = 0; a < n; ++a) x[size_t(a)] = x[size_t(a)] - g[size_t(a)] * step;
  }
  // a small jitter: the spring model leaves chains exactly straight, where an angle's force needs 1/sin θ (LAMMPS clamps
  // sin θ at 0.001 and so computes a different force there)
  std::normal_distribution<double> jit(0, 0.05);
  for (auto& p : x) p = p + Vec3{jit(rng), jit(rng), jit(rng)};
  System s;
  s.title = "beads";
  s.source_format = "caps-beads";
  std::map<std::string, int> tindex;
  for (int i = 0; i < n; ++i) {
    Atom a;
    a.id = i + 1;
    a.mol = 1;
    a.name = m.type[size_t(i)];
    a.element = 0;   // a bead, not an element
    a.charge = m.charge[size_t(i)];
    a.pos = x[size_t(i)];
    auto [it, fresh] = tindex.emplace(a.name, int(tindex.size()) + 1);
    a.type = it->second;
    if (fresh) {
      TypeInfo t;
      t.type = it->second;
      t.label = a.name;
      t.mass = o.mass ? o.mass(a.name) : 0.0;
      if (t.mass <= 0) t.mass = o.default_mass;
      s.types.push_back(t);
    }
    s.atoms.push_back(a);
  }
  for (const auto& [a, b] : m.bonds) s.bonds.push_back({uint32_t(a), uint32_t(b), 1});
  s.bonds_from_file = true;
  s.has_charges = true;   // the bead SMILES gives every bead's charge (zero where it writes none)
  s.has_mol = true;
  return s;
}

System map_to_beads(const System& s, const std::vector<BeadRule>& rules, const std::vector<BeadGroup>& groups,
                    const std::function<double(const std::string&, const std::string&)>& bond_length, BeadMapReport* rep_out) {
  BeadMapReport rep;
  const size_t n = s.atoms.size();
  const Perception p = perceive(s);
  auto heavy = [&](uint32_t a) { return s.atoms[a].element != 1; };
  auto sep = [&](const Vec3& a, const Vec3& b) { return s.cell.valid() ? s.cell.minimum_image(b - a) : b - a; };
  // molecules (connected components)
  std::vector<int> comp(n, -1);
  int ncomp = 0;
  for (uint32_t a = 0; a < n; ++a) {
    if (comp[a] >= 0) continue;
    std::vector<uint32_t> q{a};
    comp[a] = ncomp;
    for (size_t h = 0; h < q.size(); ++h)
      for (uint32_t v : p.nb[q[h]])
        if (comp[v] < 0) comp[v] = ncomp, q.push_back(v);
    ++ncomp;
  }
  std::vector<std::vector<uint32_t>> members(static_cast<size_t>(ncomp));
  for (uint32_t a = 0; a < n; ++a) members[size_t(comp[a])].push_back(a);

  struct Frag { int rule; std::vector<uint32_t> heavy; };   // rule < 0: a group (−1 − group index)
  std::vector<Frag> chosen;
  std::vector<char> taken(n, 0);

  // groups of small molecules (SDK's W = three waters): each molecule one heavy atom matching the pattern
  for (size_t g = 0; g < groups.size(); ++g) {
    const Smarts sm(groups[g].molecule);
    std::vector<uint32_t> mol;   // the heavy atom of each such molecule
    for (int c = 0; c < ncomp; ++c) {
      std::vector<uint32_t> hv;
      for (uint32_t a : members[size_t(c)])
        if (heavy(a)) hv.push_back(a);
      if (hv.size() == 1 && !taken[hv[0]] && sm.matches(s, p, hv[0], {})) mol.push_back(hv[0]);
    }
    const int k = std::max(1, groups[g].count);
    // sweep in space (a 4 Å grid, row by row) so each group gathers neighbours, not whatever came first in the file
    {
      auto key = [&](uint32_t a) {
        Vec3 r = s.atoms[a].pos;
        if (s.cell.valid()) r = s.cell.wrap(r);
        return std::array<long, 3>{long(std::floor(r[2] / 4.0)), long(std::floor(r[1] / 4.0)), long(std::floor(r[0] / 4.0))};
      };
      std::stable_sort(mol.begin(), mol.end(), [&](uint32_t a, uint32_t b) { return key(a) < key(b); });
    }
    std::vector<char> used(mol.size(), 0);
    int left = 0;
    for (size_t i = 0; i < mol.size(); ++i) {
      if (used[i]) continue;
      std::vector<std::pair<double, size_t>> d;
      for (size_t j = 0; j < mol.size(); ++j)
        if (!used[j] && j != i) d.push_back({norm(sep(s.atoms[mol[i]].pos, s.atoms[mol[j]].pos)), j});
      if (int(d.size()) < k - 1) { left += int(mol.size() - std::count(used.begin(), used.end(), 1)); break; }
      std::partial_sort(d.begin(), d.begin() + (k - 1), d.end());
      Frag f{-1 - int(g), {mol[i]}};
      used[i] = 1;
      for (int t = 0; t < k - 1; ++t) f.heavy.push_back(mol[d[size_t(t)].second]), used[d[size_t(t)].second] = 1;
      for (uint32_t a : f.heavy) taken[a] = 1;
      chosen.push_back(std::move(f));
    }
    if (left) rep.notes.push_back(std::to_string(left) + " " + groups[g].type + " molecules left over (" + std::to_string(k) + " per bead)");
  }

  // fragment candidates, first rule first; the same atom set from a later rule is dropped
  std::vector<Smarts> sms;
  for (const auto& r : rules) sms.emplace_back(r.smarts);
  std::vector<Frag> cand;
  std::set<std::vector<uint32_t>> seen;
  for (size_t r = 0; r < rules.size(); ++r)
    for (uint32_t a = 0; a < n; ++a) {
      if (!heavy(a) || taken[a]) continue;
      for (auto e : sms[r].embeddings(s, p, a, &taken)) {
        bool ok = true;
        for (uint32_t x : e) ok = ok && heavy(x);
        if (!ok) continue;
        std::sort(e.begin(), e.end());
        if (seen.insert(e).second) cand.push_back({int(r), e});
      }
    }
  std::vector<std::vector<int>> of(n);   // candidates holding each atom
  for (size_t c = 0; c < cand.size(); ++c)
    for (uint32_t a : cand[c].heavy) of[a].push_back(int(c));

  // exact cover per molecule: the uncovered atom with the fewest live candidates first (chain ends decide early)
  std::vector<int> blocked(cand.size(), 0), alive(n, 0);
  for (uint32_t a = 0; a < n; ++a) alive[a] = int(of[a].size());
  std::vector<char> covered(n, 0);
  auto choose = [&](int c, int dir) {
    for (uint32_t a : cand[size_t(c)].heavy) {
      covered[a] = dir > 0;
      for (int c2 : of[a]) {
        if (dir > 0) {
          if (blocked[size_t(c2)]++ == 0)
            for (uint32_t b : cand[size_t(c2)].heavy) --alive[b];
        } else {
          if (--blocked[size_t(c2)] == 0)
            for (uint32_t b : cand[size_t(c2)].heavy) ++alive[b];
        }
      }
    }
  };
  int mapped_mols = 0;
  for (int c = 0; c < ncomp; ++c) {
    std::vector<uint32_t> hv;
    for (uint32_t a : members[size_t(c)])
      if (heavy(a) && !taken[a]) hv.push_back(a);
    if (hv.empty()) continue;
    struct Level { uint32_t atom; size_t next; int cand; };
    std::vector<Level> st;
    long budget = 2000000;
    bool ok = false;
    for (;;) {
      if (--budget < 0) break;
      int best = -1, bestn = 1 << 30;
      for (uint32_t a : hv)
        if (!covered[a] && alive[a] < bestn) bestn = alive[a], best = int(a);
      if (best < 0) { ok = true; break; }
      bool advanced = false;
      if (bestn > 0) {
        for (size_t k = 0; k < of[size_t(best)].size(); ++k) {
          const int cc = of[size_t(best)][k];
          if (blocked[size_t(cc)]) continue;
          choose(cc, +1);
          st.push_back({uint32_t(best), k + 1, cc});
          advanced = true;
          break;
        }
      }
      if (advanced) continue;
      // backtrack: the next live candidate of the deepest level that has one
      bool resumed = false;
      while (!st.empty() && !resumed) {
        Level lv = st.back();
        st.pop_back();
        choose(lv.cand, -1);
        for (size_t k = lv.next; k < of[lv.atom].size(); ++k) {
          const int cc = of[lv.atom][k];
          if (blocked[size_t(cc)]) continue;
          choose(cc, +1);
          st.push_back({lv.atom, k + 1, cc});
          resumed = true;
          break;
        }
      }
      if (!resumed) break;
    }
    if (ok) {
      for (const auto& lv : st) chosen.push_back(cand[size_t(lv.cand)]);
      ++mapped_mols;
    } else {
      for (const auto& lv : st) choose(lv.cand, -1);
      for (uint32_t a : members[size_t(c)]) rep.uncovered.push_back(a);
    }
  }

  // beads: fragment heavy atoms plus their hydrogens, at the centre of mass
  System out;
  out.title = s.title;
  out.cell = s.cell;
  out.source_format = s.source_format;
  std::vector<int> bead_of(n, -1);
  std::vector<std::string> btype;
  std::vector<int> brule;
  for (const auto& f : chosen) {
    std::vector<uint32_t> at = f.heavy;
    for (uint32_t a : f.heavy)
      for (uint32_t v : p.nb[a])
        if (s.atoms[v].element == 1) at.push_back(v);
    const Vec3 ref = s.atoms[at[0]].pos;
    Vec3 c{0, 0, 0};
    double m = 0, q = 0;
    for (uint32_t a : at) {
      const double w = s.mass_of(s.atoms[a]);
      c = c + (ref + sep(ref, s.atoms[a].pos)) * w;
      m += w;
      q += s.atoms[a].charge;
      bead_of[a] = int(out.atoms.size());
    }
    Atom b;
    b.id = int64_t(out.atoms.size() + 1);
    b.mol = s.atoms[at[0]].mol;
    b.resname = s.atoms[at[0]].resname;
    b.resid = s.atoms[at[0]].resid;
    b.element = 0;
    b.charge = q;
    b.pos = c * (1 / m);
    const std::string t = f.rule >= 0 ? rules[size_t(f.rule)].types.front() : groups[size_t(-1 - f.rule)].type;
    b.name = t;
    btype.push_back(t);
    brule.push_back(f.rule);
    out.atoms.push_back(b);
    TypeInfo ti;
    ti.mass = m;
    out.types.push_back(ti);   // one entry per bead for now: merged by name below
  }
  // atoms not covered stay as they were
  std::vector<int> keep_of(n, -1);
  for (uint32_t a = 0; a < n; ++a) {
    if (bead_of[a] >= 0) continue;
    keep_of[a] = int(out.atoms.size());
    Atom x = s.atoms[a];
    x.id = int64_t(out.atoms.size() + 1);
    out.atoms.push_back(x);
    out.types.push_back(TypeInfo{});
  }
  auto site = [&](uint32_t a) { return bead_of[a] >= 0 ? bead_of[a] : keep_of[a]; };
  std::set<std::pair<int, int>> bonds;
  for (const auto& b : s.bonds) {
    const int u = site(b.i), v = site(b.j);
    if (u != v) bonds.insert({std::min(u, v), std::max(u, v)});
  }
  for (const auto& [u, v] : bonds) out.bonds.push_back({uint32_t(u), uint32_t(v), 1});
  // alternatives (SDK's EST1 / EST2 on one GL): the assignment whose distances best match the force field's bond lengths
  std::vector<std::vector<uint32_t>> bnb(out.atoms.size());
  for (const auto& b : out.bonds) bnb[b.i].push_back(b.j), bnb[b.j].push_back(b.i);
  int resolved = 0, unresolved = 0;
  for (size_t r = 0; r < rules.size(); ++r) {
    const auto& alts = rules[r].types;
    if (alts.size() < 2) continue;
    std::map<uint32_t, std::vector<uint32_t>> by_anchor;   // a neighbouring bead of another rule → the alternatives on it
    for (uint32_t b = 0; b < btype.size(); ++b)
      if (brule[b] == int(r))
        for (uint32_t nbb : bnb[b])
          if (nbb < btype.size() && brule[nbb] != int(r)) by_anchor[nbb].push_back(b);
    std::set<uint32_t> done;
    for (auto& [anc, list] : by_anchor) {
      if (list.size() != alts.size()) continue;
      std::vector<size_t> perm(alts.size());
      for (size_t k = 0; k < perm.size(); ++k) perm[k] = k;
      double bestc = 1e300;
      std::vector<size_t> bestp = perm;
      bool known = false;
      do {
        double c = 0;
        for (size_t k = 0; k < list.size(); ++k) {
          const double r0 = bond_length ? bond_length(alts[perm[k]], btype[anc]) : 0.0;
          if (r0 > 0) known = true;
          const double d = norm(sep(out.atoms[anc].pos, out.atoms[list[k]].pos));
          c += r0 > 0 ? (d - r0) * (d - r0) : 0.0;
        }
        if (c < bestc) bestc = c, bestp = perm;
      } while (std::next_permutation(perm.begin(), perm.end()));
      for (size_t k = 0; k < list.size(); ++k) {
        btype[list[k]] = alts[bestp[k]];
        out.atoms[list[k]].name = alts[bestp[k]];
        done.insert(list[k]);
      }
      known ? ++resolved : ++unresolved;
    }
    for (uint32_t b = 0; b < btype.size(); ++b)
      if (brule[b] == int(r) && !done.count(b)) ++unresolved;
  }
  if (resolved) rep.notes.push_back(std::to_string(resolved) + " sets of alternative beads told apart by the force field's bond lengths");
  if (unresolved) rep.notes.push_back(std::to_string(unresolved) + " alternative beads could not be told apart (the first type used)");
  // one type per name, its mass the first bead's
  std::map<std::string, int> tid;
  std::vector<TypeInfo> types;
  for (size_t i = 0; i < out.atoms.size(); ++i) {
    auto& a = out.atoms[i];
    const bool bead = i < btype.size();
    const std::string key = bead ? a.name : "atom:" + a.name + ":" + std::to_string(a.type);
    auto [it, fresh] = tid.emplace(key, int(tid.size()) + 1);
    if (fresh) {
      TypeInfo t;
      t.type = it->second;
      t.label = a.name;
      t.mass = bead ? out.types[i].mass : 0.0;   // atoms: mass from the element
      types.push_back(t);
    }
    a.type = it->second;
    if (bead) rep.by_type[a.name]++;
  }
  out.types = types;
  out.bonds_from_file = true;
  out.has_charges = s.has_charges;
  out.has_mol = s.has_mol;
  rep.beads = int(btype.size());
  rep.molecules_mapped = mapped_mols;
  if (!rep.uncovered.empty())
    rep.notes.push_back(std::to_string(rep.uncovered.size()) + " atoms in molecules the fragments cannot cover are left as atoms");
  if (rep_out) *rep_out = std::move(rep);
  return out;
}

}  // namespace caps
