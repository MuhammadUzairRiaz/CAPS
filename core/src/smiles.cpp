// SMILES parsing, hydrogens, valence checks and formula (see molecule.hpp).
#include <algorithm>
#include <cctype>
#include <cmath>
#include <functional>
#include <map>
#include <set>
#include <sstream>

#include "caps/elements.hpp"
#include "caps/molecule.hpp"

namespace caps {
namespace {

// Organic subset normal valences (OpenSMILES 3.5).
std::vector<int> normal_valences(int z) {
  switch (z) {
    case 5: return {3};
    case 6: return {4};
    case 7: return {3, 5};
    case 8: return {2};
    case 15: return {3, 5};
    case 16: return {2, 4, 6};
    case 9: case 17: case 35: case 53: return {1};
    default: return {};
  }
}

bool organic(int z) { return !normal_valences(z).empty(); }

int bond_valence(int order) { return order == 4 ? 1 : order; }

struct Ring {
  int atom, order, dir, slot;
  size_t pos;
};

}  // namespace

MolGraph parse_smiles(const std::string& text) {
  MolGraph g;
  g.smiles = text;
  std::string s = text;
  while (!s.empty() && std::isspace(static_cast<unsigned char>(s.back()))) s.pop_back();
  size_t i = 0;
  while (i < s.size() && std::isspace(static_cast<unsigned char>(s[i]))) ++i;
  if (i == s.size()) throw SmilesError("empty SMILES", 0);

  int prev = -1, pending = 0, pdir = 0;
  std::vector<int> branches;
  std::map<int, Ring> rings;
  g.parts = 1;
  bool expect_atom = true;   // after '.', '(' or a bond an atom must follow

  auto add_bond = [&](int a, int b, int order, int dir, size_t pos) {
    for (const auto& x : g.bonds)
      if ((x.a == a && x.b == b) || (x.a == b && x.b == a)) throw SmilesError("two bonds between the same atoms", pos);
    if (a == b) throw SmilesError("a ring closure bonds an atom to itself", pos);
    if (order == 0) order = g.atoms[a].aromatic && g.atoms[b].aromatic ? 4 : 1;
    g.bonds.push_back({a, b, order, dir});
  };

  auto add_atom = [&](MolAtom a, size_t pos) {
    const int k = int(g.atoms.size());
    if (prev >= 0) {
      a.order.insert(a.order.begin(), prev);
      g.atoms[prev].order.push_back(k);
    }
    g.atoms.push_back(std::move(a));
    if (prev >= 0) add_bond(prev, k, pending, pdir, pos);
    pending = 0;
    pdir = 0;
    prev = k;
    expect_atom = false;
  };

  while (i < s.size()) {
    const char c = s[i];
    const size_t at = i;
    if (c == '(') {
      if (prev < 0) throw SmilesError("a branch must follow an atom", at);
      branches.push_back(prev);
      expect_atom = true;
      ++i;
    } else if (c == ')') {
      if (branches.empty()) throw SmilesError("')' without '('", at);
      if (expect_atom || pending) throw SmilesError("empty branch or a bond with no atom", at);
      prev = branches.back();
      branches.pop_back();
      ++i;
    } else if (c == '-' || c == '=' || c == '#' || c == '$' || c == ':' || c == '/' || c == '\\') {
      if (prev < 0) throw SmilesError(std::string("bond '") + c + "' with no atom before it", at);
      if (pending) throw SmilesError("two bond symbols in a row", at);
      pending = c == '=' ? 2 : c == '#' || c == '$' ? 3 : c == ':' ? 4 : 1;
      pdir = c == '/' ? 1 : c == '\\' ? -1 : 0;
      expect_atom = true;
      ++i;
    } else if (c == '.') {
      if (prev < 0 || pending) throw SmilesError("'.' must separate two parts", at);
      prev = -1;
      ++g.parts;
      expect_atom = true;
      ++i;
    } else if (std::isdigit(static_cast<unsigned char>(c)) || c == '%') {
      if (prev < 0) throw SmilesError("a ring-closure digit must follow an atom", at);
      int num;
      if (c == '%') {
        if (i + 2 >= s.size() + 0 || !std::isdigit(static_cast<unsigned char>(s[i + 1])) || !std::isdigit(static_cast<unsigned char>(s[i + 2])))
          throw SmilesError("'%' must be followed by two digits", at);
        num = (s[i + 1] - '0') * 10 + (s[i + 2] - '0');
        i += 3;
      } else {
        num = c - '0';
        ++i;
      }
      auto it = rings.find(num);
      if (it == rings.end()) {
        rings[num] = {prev, pending, pdir, int(g.atoms[prev].order.size()), at};
        g.atoms[prev].order.push_back(-1);
      } else {
        Ring r = it->second;
        rings.erase(it);
        if (r.order && pending && r.order != pending) throw SmilesError("ring closure " + std::to_string(num) + " has two different bond orders", at);
        const int order = pending ? pending : r.order;
        // direction: written at the opening (r.atom → prev) or at the closure (prev → r.atom)
        int dir = r.dir ? r.dir : 0;
        int a = r.atom, b = prev;
        if (!r.dir && pdir) { dir = pdir; a = prev; b = r.atom; }
        add_bond(a, b, order, dir, at);
        g.atoms[r.atom].order[size_t(r.slot)] = prev;
        g.atoms[prev].order.push_back(r.atom);
        pending = 0;
        pdir = 0;
      }
    } else if (c == '[') {
      const size_t end = s.find(']', i);
      if (end == std::string::npos) throw SmilesError("'[' without ']'", at);
      const std::string b = s.substr(i + 1, end - i - 1);
      MolAtom a;
      a.bracket = true;
      a.hcount = 0;
      size_t k = 0;
      while (k < b.size() && std::isdigit(static_cast<unsigned char>(b[k]))) a.isotope = a.isotope * 10 + (b[k++] - '0');
      if (k >= b.size()) throw SmilesError("bracket atom without an element", at);
      // element: aromatic lower case (b c n o p s se as te) or a symbol
      std::string sym;
      if (std::islower(static_cast<unsigned char>(b[k]))) {
        if (b.compare(k, 2, "se") == 0 || b.compare(k, 2, "as") == 0 || b.compare(k, 2, "te") == 0) sym = b.substr(k, 2), k += 2;
        else sym = b.substr(k, 1), ++k;
        a.aromatic = true;
        sym[0] = char(std::toupper(static_cast<unsigned char>(sym[0])));
      } else if (b[k] == '*') {
        sym = "*";
        ++k;
      } else {
        sym = b.substr(k, 1);
        ++k;
        if (k < b.size() && std::islower(static_cast<unsigned char>(b[k]))) {
          const std::string two = sym + b[k];
          if (element_from_symbol(two) > 0) sym = two, ++k;
        }
      }
      a.element = sym == "*" ? 0 : element_from_symbol(sym);
      if (sym != "*" && a.element == 0) throw SmilesError("unknown element '" + sym + "'", at + 1);
      if (k < b.size() && b[k] == '@') {
        ++k;
        a.chiral = 1;
        if (k < b.size() && b[k] == '@') { a.chiral = 2; ++k; }
        else if (b.compare(k, 2, "TH") == 0 && k + 2 < b.size()) { a.chiral = b[k + 2] == '2' ? 2 : 1; k += 3; }
        else if (k + 1 < b.size() && std::isupper(static_cast<unsigned char>(b[k])) && std::isupper(static_cast<unsigned char>(b[k + 1]))) {
          a.chiral = 0;   // allene, square-planar, bipyramidal and octahedral classes: not used
          k += 2;
          while (k < b.size() && std::isdigit(static_cast<unsigned char>(b[k]))) ++k;
        }
      }
      if (k < b.size() && b[k] == 'H') {
        ++k;
        a.hcount = 1;
        if (k < b.size() && std::isdigit(static_cast<unsigned char>(b[k]))) a.hcount = b[k++] - '0';
      }
      if (k < b.size() && (b[k] == '+' || b[k] == '-')) {
        const char sign = b[k++];
        int q = 1;
        if (k < b.size() && std::isdigit(static_cast<unsigned char>(b[k]))) {
          q = 0;
          while (k < b.size() && std::isdigit(static_cast<unsigned char>(b[k]))) q = q * 10 + (b[k++] - '0');
        } else
          while (k < b.size() && b[k] == sign) ++q, ++k;
        a.charge = sign == '+' ? q : -q;
      }
      if (k < b.size() && b[k] == ':') {
        ++k;
        while (k < b.size() && std::isdigit(static_cast<unsigned char>(b[k]))) a.map = a.map * 10 + (b[k++] - '0');
      }
      if (k != b.size()) throw SmilesError("cannot read bracket atom [" + b + "]", at);
      if (a.element == 1 && a.hcount > 0) throw SmilesError("[H] cannot carry hydrogens", at);
      if (a.chiral && a.hcount == 1) a.order.push_back(-2);   // its own hydrogen follows the atom before it
      add_atom(std::move(a), at);
      i = end + 1;
    } else if (std::isalpha(static_cast<unsigned char>(c)) || c == '*') {
      MolAtom a;
      if (c == '*') { a.element = 0; ++i; }
      else if (s.compare(i, 2, "Cl") == 0) { a.element = 17; i += 2; }
      else if (s.compare(i, 2, "Br") == 0) { a.element = 35; i += 2; }
      else {
        const char u = char(std::toupper(static_cast<unsigned char>(c)));
        const std::string organic_set = "BCNOPSFI";
        if (organic_set.find(u) == std::string::npos) {
          if (c == 'H') throw SmilesError("hydrogen outside brackets: write [H] or leave hydrogens implicit", at);
          throw SmilesError(std::string("'") + c + "' is not in the organic subset; write it in brackets, e.g. [" + c + "]", at);
        }
        if (std::islower(static_cast<unsigned char>(c))) {
          if (u == 'F' || u == 'I') throw SmilesError(std::string("'") + c + "' cannot be aromatic", at);
          a.aromatic = true;
        }
        a.element = element_from_symbol(std::string(1, u));
        ++i;
      }
      add_atom(std::move(a), at);
    } else if (std::isspace(static_cast<unsigned char>(c))) {
      break;   // a name or comment after the SMILES
    } else {
      throw SmilesError(std::string("unexpected '") + c + "'", at);
    }
  }
  if (!branches.empty()) throw SmilesError("'(' without ')'", s.size());
  if (!rings.empty()) throw SmilesError("ring bond " + std::to_string(rings.begin()->first) + " is not closed", rings.begin()->second.pos);
  if (pending || expect_atom) throw SmilesError("the SMILES ends with a bond or a branch", s.size());
  g.heavy = int(g.atoms.size());
  return g;
}

void add_hydrogens(MolGraph& g) {
  const int n0 = int(g.atoms.size());
  std::vector<int> bsum(size_t(n0), 0), nb(size_t(n0), 0), narom(size_t(n0), 0);
  for (const auto& b : g.bonds) {
    bsum[size_t(b.a)] += bond_valence(b.order);
    bsum[size_t(b.b)] += bond_valence(b.order);
    ++nb[size_t(b.a)];
    ++nb[size_t(b.b)];
    if (b.order == 4) ++narom[size_t(b.a)], ++narom[size_t(b.b)];
  }
  for (int i = 0; i < n0; ++i) {
    MolAtom& a = g.atoms[size_t(i)];
    int h = 0;
    if (a.hcount >= 0) h = a.hcount;
    else if (organic(a.element)) {
      if (a.aromatic) {
        // one π electron from c, n (pyridine type), p, b; o and s give a lone pair and take no hydrogen
        if (a.element == 8 || a.element == 16) h = 0;
        else {
          const int target = a.element == 6 ? 4 : 3;
          h = std::max(0, target - (bsum[size_t(i)] + 1));
        }
      } else {
        const int v = bsum[size_t(i)];
        for (int nv : normal_valences(a.element))
          if (nv >= v) { h = nv - v; break; }
      }
    }
    for (int k = 0; k < h; ++k) {
      MolAtom H;
      H.element = 1;
      H.hcount = 0;
      H.order.push_back(i);
      const int hi = int(g.atoms.size());
      g.atoms.push_back(H);
      MolAtom& p = g.atoms[size_t(i)];   // (the vector may have grown)
      g.bonds.push_back({i, hi, 1, 0});
      auto own = std::find(p.order.begin(), p.order.end(), -2);
      if (own != p.order.end()) *own = hi;   // a bracket's own hydrogen keeps its written place
      else p.order.push_back(hi);
    }
    g.atoms[size_t(i)].hcount = h;
  }
}

std::vector<std::string> valence_problems(const MolGraph& g) {
  std::vector<std::string> out;
  std::vector<int> v(g.atoms.size(), 0);
  std::vector<int> ar(g.atoms.size(), 0);
  for (const auto& b : g.bonds) {
    v[size_t(b.a)] += bond_valence(b.order);
    v[size_t(b.b)] += bond_valence(b.order);
    if (b.order == 4) ++ar[size_t(b.a)], ++ar[size_t(b.b)];
  }
  for (size_t i = 0; i < g.atoms.size(); ++i) {
    const MolAtom& a = g.atoms[i];
    if (a.element == 0) continue;
    const auto nv = normal_valences(a.element);
    if (a.aromatic && ar[i] == 0 && !a.bracket)
      out.push_back("atom " + std::to_string(i + 1) + " (" + element(a.element).symbol + ") is aromatic but not in an aromatic ring");
    if (nv.empty()) continue;
    int used = v[i] + (a.aromatic && a.element != 8 && a.element != 16 && ar[i] ? 1 : 0);
    int maxv = nv.back();
    if (a.charge) {
      // N+, P+, O+, S+ gain a bond; C+, C- , B- lose or gain one; anions of O, S, halogens lose one
      if (a.element == 7 || a.element == 15 || a.element == 8 || a.element == 16) maxv += a.charge > 0 ? a.charge : a.charge;
      else if (a.element == 6) maxv = 3;
      else if (a.element == 5 && a.charge < 0) maxv = 4;
      else maxv -= std::abs(a.charge);
    }
    if (used > std::max(maxv, 0))
      out.push_back("atom " + std::to_string(i + 1) + " (" + element(a.element).symbol + ") has valence " + std::to_string(used) +
                    ", more than " + std::to_string(std::max(maxv, 0)));
  }
  return out;
}

MolInfo molecule_info(const MolGraph& g) {
  MolInfo m;
  std::map<int, int> count;
  for (const auto& a : g.atoms) {
    ++count[a.element];
    m.mass += element(a.element).mass;
    m.charge += a.charge;
    if (a.element != 1) ++m.heavy;
    if (a.chiral) ++m.stereocentres;
  }
  m.atoms = int(g.atoms.size());
  m.bonds = int(g.bonds.size());
  // rings: bonds − atoms + connected components
  std::vector<int> comp(g.atoms.size(), -1);
  int ncomp = 0;
  std::vector<std::vector<int>> adj(g.atoms.size());
  for (const auto& b : g.bonds) adj[size_t(b.a)].push_back(b.b), adj[size_t(b.b)].push_back(b.a);
  for (size_t s = 0; s < g.atoms.size(); ++s) {
    if (comp[s] >= 0) continue;
    std::vector<int> st{int(s)};
    comp[s] = ncomp;
    while (!st.empty()) {
      const int u = st.back();
      st.pop_back();
      for (int w : adj[size_t(u)])
        if (comp[size_t(w)] < 0) comp[size_t(w)] = ncomp, st.push_back(w);
    }
    ++ncomp;
  }
  m.rings = m.bonds - m.atoms + ncomp;
  // E/Z: double bonds with a directional single bond on each end
  for (const auto& b : g.bonds) {
    if (b.order != 2) continue;
    bool ea = false, eb = false;
    for (const auto& x : g.bonds) {
      if (!x.dir || &x == &b) continue;
      if (x.a == b.a || x.b == b.a) ea = true;
      if (x.a == b.b || x.b == b.b) eb = true;
    }
    m.stereo_bonds += ea && eb;
  }
  // Hill order: C, H, then alphabetical (no carbon: all alphabetical)
  auto sym = [](int z) { return std::string(z ? element(z).symbol : "*"); };
  std::vector<std::pair<std::string, int>> rest;
  std::ostringstream f;
  const bool carbon = count.count(6) > 0;
  for (const auto& [z, k] : count) {
    if (carbon && (z == 6 || z == 1)) continue;
    rest.push_back({sym(z), k});
  }
  std::sort(rest.begin(), rest.end());
  auto put = [&](const std::string& e, int k) { f << e; if (k > 1) f << k; };
  if (carbon) {
    put("C", count[6]);
    if (count.count(1)) put("H", count[1]);
  }
  for (const auto& [e, k] : rest) put(e, k);
  if (m.charge) f << (m.charge > 0 ? "+" : "-") << (std::abs(m.charge) > 1 ? std::to_string(std::abs(m.charge)) : "");
  m.formula = f.str();
  m.problems = valence_problems(g);
  return m;
}

std::string write_smiles(const MolGraph& g) {
  const int n = g.heavy > 0 ? std::min<int>(g.heavy, int(g.atoms.size())) : int(g.atoms.size());
  std::vector<std::vector<std::pair<int, int>>> adj(static_cast<size_t>(n));   // (neighbour, bond index)
  for (size_t k = 0; k < g.bonds.size(); ++k) {
    const auto& b = g.bonds[k];
    if (b.a >= n || b.b >= n) continue;
    adj[size_t(b.a)].push_back({b.b, int(k)});
    adj[size_t(b.b)].push_back({b.a, int(k)});
  }
  // pass 1: spanning tree and ring closures
  std::vector<int> seen(static_cast<size_t>(n), 0);
  std::vector<std::vector<std::pair<int, int>>> kids(static_cast<size_t>(n));
  struct Closure { int anc, desc, bond, number = 0; };
  std::vector<Closure> rings;
  std::vector<char> used(g.bonds.size(), 0);
  std::function<void(int)> dfs = [&](int u) {
    seen[size_t(u)] = 1;
    for (auto [v, bi] : adj[size_t(u)]) {
      if (used[size_t(bi)]) continue;
      used[size_t(bi)] = 1;
      if (!seen[size_t(v)]) {
        kids[size_t(u)].push_back({v, bi});
        dfs(v);
      } else
        rings.push_back({v, u, bi});   // v is an ancestor still open: the ring bond opens at v, closes at u
    }
  };
  std::vector<int> roots;
  for (int i = 0; i < n; ++i)
    if (!seen[size_t(i)]) roots.push_back(i), dfs(i);

  auto sym = [&](const MolBond& b, int from) {
    const bool ar = g.atoms[size_t(b.a)].aromatic && g.atoms[size_t(b.b)].aromatic;
    if (b.order == 2) return std::string("=");
    if (b.order == 3) return std::string("#");
    if (b.order == 4) return ar ? std::string() : std::string(":");
    if (b.dir) {
      const int d = from == b.a ? b.dir : -b.dir;
      return std::string(d > 0 ? "/" : "\\");
    }
    return ar ? std::string("-") : std::string();
  };
  std::vector<int> free_numbers;
  int next_number = 1;
  std::string out;
  std::function<void(int, int)> write = [&](int u, int from_bond) {
    const MolAtom& a = g.atoms[size_t(u)];
    // neighbours in the order this writing gives them
    std::vector<int> w;
    if (from_bond >= 0) w.push_back(g.bonds[size_t(from_bond)].a == u ? g.bonds[size_t(from_bond)].b : g.bonds[size_t(from_bond)].a);
    const bool own_h = a.bracket && a.hcount > 0 && std::find(a.order.begin(), a.order.end(), -2) != a.order.end();
    if (own_h) w.push_back(-2);
    std::string digits;
    for (auto& r : rings)
      if (r.desc == u) {
        digits += (r.number < 10 ? std::to_string(r.number) : "%" + std::to_string(r.number));
        free_numbers.push_back(r.number);
        w.push_back(r.anc);
      }
    for (auto& r : rings)
      if (r.anc == u) {
        if (!free_numbers.empty()) {
          std::sort(free_numbers.begin(), free_numbers.end());
          r.number = free_numbers.front();
          free_numbers.erase(free_numbers.begin());
        } else
          r.number = next_number++;
        digits += sym(g.bonds[size_t(r.bond)], u) + (r.number < 10 ? std::to_string(r.number) : "%" + std::to_string(r.number));
        w.push_back(r.desc);
      }
    for (auto [v, bi] : kids[size_t(u)]) w.push_back(v);
    // chirality: the written order is a permutation of the parsed one; an odd permutation swaps @ and @@
    int chiral = 0;
    if (a.chiral) {
      std::vector<int> o = a.order;
      std::vector<int> ws = w, os = o;
      std::sort(ws.begin(), ws.end());
      std::sort(os.begin(), os.end());
      if (ws == os) {
        int inv = 0;
        std::vector<int> idx;
        for (int x : w) idx.push_back(int(std::find(o.begin(), o.end(), x) - o.begin()));
        for (size_t p = 0; p < idx.size(); ++p)
          for (size_t q = p + 1; q < idx.size(); ++q) inv += idx[p] > idx[q];
        chiral = inv % 2 ? 3 - a.chiral : a.chiral;
      }
    }
    std::string e = a.element ? element(a.element).symbol : "*";
    if (a.aromatic) for (auto& c : e) c = char(std::tolower(static_cast<unsigned char>(c)));
    const bool bracket = a.bracket || (a.element && !organic(a.element)) || a.charge || a.isotope || a.map || chiral || a.element == 1;
    if (bracket) {
      std::string t = "[";
      if (a.isotope) t += std::to_string(a.isotope);
      t += e;
      if (chiral) t += chiral == 1 ? "@" : "@@";
      if (a.hcount > 0) t += a.hcount == 1 ? "H" : "H" + std::to_string(a.hcount);
      if (a.charge) t += (a.charge > 0 ? "+" : "-") + (std::abs(a.charge) > 1 ? std::to_string(std::abs(a.charge)) : "");
      if (a.map) t += ":" + std::to_string(a.map);
      e = t + "]";
    }
    out += e + digits;
    const auto& ks = kids[size_t(u)];
    for (size_t c = 0; c < ks.size(); ++c) {
      const auto [v, bi] = ks[c];
      const bool last = c + 1 == ks.size();
      if (!last) out += "(";
      out += sym(g.bonds[size_t(bi)], u);
      write(v, bi);
      if (!last) out += ")";
    }
  };
  for (size_t r = 0; r < roots.size(); ++r) {
    if (r) out += ".";
    write(roots[r], -1);
  }
  return out;
}

System molecule_system(const MolGraph& g, const std::vector<Vec3>& pos) {
  System s;
  s.title = g.smiles;
  std::map<int, int> seen;
  for (size_t i = 0; i < g.atoms.size(); ++i) {
    Atom a;
    a.id = int64_t(i + 1);
    a.mol = 1;
    a.element = g.atoms[i].element;
    a.charge = g.atoms[i].charge;
    a.name = std::string(a.element ? element(a.element).symbol : "X") + std::to_string(++seen[a.element]);
    a.resname = "MOL";
    a.pos = i < pos.size() ? pos[i] : Vec3{0, 0, 0};
    s.atoms.push_back(a);
  }
  for (const auto& b : g.bonds) s.bonds.push_back({uint32_t(b.a), uint32_t(b.b), b.order});
  return s;
}

}  // namespace caps
