#include "caps/query.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <functional>
#include <set>
#include <stdexcept>

#include "caps/edit.hpp"
#include "caps/elements.hpp"
#include "caps/typing.hpp"

namespace caps {

namespace {

struct Tok {
  enum Kind { Word, Number, String, LParen, RParen, End } kind = End;
  std::string text;
};

std::vector<Tok> lex(const std::string& q) {
  std::vector<Tok> out;
  size_t i = 0;
  while (i < q.size()) {
    const char c = q[i];
    if (std::isspace((unsigned char)c) || c == ',') { ++i; continue; }
    if (c == '(') { out.push_back({Tok::LParen, "("}); ++i; continue; }
    if (c == ')') { out.push_back({Tok::RParen, ")"}); ++i; continue; }
    if (c == '"' || c == '\'') {
      const size_t j = q.find(c, i + 1);
      if (j == std::string::npos) throw std::invalid_argument("unclosed quote in the query");
      out.push_back({Tok::String, q.substr(i + 1, j - i - 1)});
      i = j + 1;
      continue;
    }
    size_t j = i;
    while (j < q.size() && !std::isspace((unsigned char)q[j]) && q[j] != '(' && q[j] != ')' && q[j] != ',' && q[j] != '"') ++j;
    std::string w = q.substr(i, j - i);
    const bool num = !w.empty() && (std::isdigit((unsigned char)w[0]) || ((w[0] == '-' || w[0] == '.') && w.size() > 1 && std::isdigit((unsigned char)w[1])));
    out.push_back({num ? Tok::Number : Tok::Word, w});
    i = j;
  }
  out.push_back({Tok::End, ""});
  return out;
}

std::string lower(std::string s) { for (auto& c : s) c = char(std::tolower((unsigned char)c)); return s; }

bool keyword(const std::string& w) {
  static const std::set<std::string> k = {"and", "or", "not", "smarts", "element", "type", "chain", "molecule", "index", "ring", "stereo",
                                          "within", "of", "sel", "selection", "all", "none", "x", "y", "z", "charge", "mass", "bonds", "hybrid", "tag"};
  return k.count(lower(w)) > 0;
}

// 1-4, 7 → {1,2,3,4,7} (tokens already split at commas)
std::vector<long> ranges(const std::vector<Tok>& t, size_t& p) {
  std::vector<long> v;
  while (t[p].kind == Tok::Number) {
    const std::string& w = t[p].text;
    const auto dash = w.find('-', 1);
    if (dash != std::string::npos) {
      const long a = std::stol(w.substr(0, dash)), b = std::stol(w.substr(dash + 1));
      for (long x = std::min(a, b); x <= std::max(a, b); ++x) v.push_back(x);
    } else {
      v.push_back(std::stol(w));
    }
    ++p;
  }
  if (v.empty()) throw std::invalid_argument("expected numbers or ranges such as 1-4 after '" + t[p > 0 ? p - 1 : 0].text + "'");
  return v;
}

struct Parser {
  const System& s;
  const std::vector<char>& current;
  std::vector<Tok> t;
  size_t p = 0;
  const Perception* per = nullptr;
  std::vector<char> cip;
  bool cip_done = false;
  std::unique_ptr<Perception> own;

  const Perception& perception() {
    if (!per) { own = std::make_unique<Perception>(perceive(s)); per = own.get(); }
    return *per;
  }
  std::vector<char> none() const { return std::vector<char>(s.atoms.size(), 0); }

  std::vector<char> expr() {
    auto a = term();
    while (lower(t[p].text) == "or" && t[p].kind == Tok::Word) {
      ++p;
      const auto b = term();
      for (size_t i = 0; i < a.size(); ++i) a[i] = a[i] || b[i];
    }
    return a;
  }
  std::vector<char> term() {
    auto a = factor();
    while (lower(t[p].text) == "and" && t[p].kind == Tok::Word) {
      ++p;
      const auto b = factor();
      for (size_t i = 0; i < a.size(); ++i) a[i] = a[i] && b[i];
    }
    return a;
  }
  std::vector<char> factor() {
    if (t[p].kind == Tok::Word && lower(t[p].text) == "not") {
      ++p;
      auto a = factor();
      for (auto& x : a) x = !x;
      return a;
    }
    if (t[p].kind == Tok::LParen) {
      ++p;
      auto a = expr();
      if (t[p].kind != Tok::RParen) throw std::invalid_argument("missing ')' in the query");
      ++p;
      return a;
    }
    return primary();
  }
  std::vector<char> primary() {
    if (t[p].kind != Tok::Word) throw std::invalid_argument(t[p].kind == Tok::End ? "the query ends too early" : "unexpected '" + t[p].text + "'");
    const std::string w = lower(t[p].text);
    ++p;
    const size_t n = s.atoms.size();
    auto words = [&] {
      std::string out;
      while ((t[p].kind == Tok::Word && !keyword(t[p].text)) || t[p].kind == Tok::String) out += (out.empty() ? "" : " ") + t[p++].text;
      if (out.empty()) throw std::invalid_argument("'" + w + "' needs a value");
      return out;
    };
    if (w == "smarts") {
      if (t[p].kind != Tok::String && t[p].kind != Tok::Word) throw std::invalid_argument("smarts needs a pattern, e.g. smarts \"c1ccccc1\"");
      return select_smarts(s, t[p++].text);
    }
    if (w == "element") return select_element(s, words());
    if (w == "type") {
      auto out = none();
      for (const auto& label : [&] { std::vector<std::string> v; std::string x = words(), cur; for (char c : x + " ") { if (c == ' ') { if (!cur.empty()) v.push_back(cur); cur.clear(); } else cur += c; } return v; }()) {
        const auto m = select_type(s, label);
        for (size_t i = 0; i < n; ++i) out[i] = out[i] || m[i];
      }
      return out;
    }
    if (w == "chain" || w == "molecule") {
      const auto want = ranges(t, p);
      const auto mol = s.molecules();
      std::set<long> ws(want.begin(), want.end());
      auto out = none();
      for (size_t i = 0; i < n; ++i) out[i] = ws.count(long(mol[i]) + 1) ? 1 : 0;
      return out;
    }
    if (w == "index") {
      auto out = none();
      for (long k : ranges(t, p)) if (k >= 1 && size_t(k) <= n) out[size_t(k - 1)] = 1;
      return out;
    }
    if (w == "ring") {
      const auto& rings = perception().rings;
      auto out = none();
      for (long k : ranges(t, p)) {
        if (k < 1 || size_t(k) > rings.size()) throw std::invalid_argument("there is no ring " + std::to_string(k) + " (" + std::to_string(rings.size()) + " rings)");
        for (uint32_t a : rings[size_t(k - 1)]) out[a] = 1;
      }
      return out;
    }
    if (w == "stereo") {
      if (!cip_done) { cip = cip_labels(s); cip_done = true; }
      const std::string v = t[p].kind == Tok::End ? "*" : t[p++].text;
      auto out = none();
      for (size_t i = 0; i < n; ++i) {
        if (!cip[i]) continue;
        if (v == "*" || (v == "R" || v == "r" ? cip[i] == 'R' : v == "S" || v == "s" ? cip[i] == 'S' : false)) out[i] = 1;
        if (v != "*" && v != "R" && v != "r" && v != "S" && v != "s") throw std::invalid_argument("stereo takes R, S or *");
      }
      return out;
    }
    if (w == "within") {
      if (t[p].kind != Tok::Number) throw std::invalid_argument("within needs a distance, e.g. within 5.0 of sel");
      const double d = std::stod(t[p++].text);
      if (lower(t[p].text) != "of") throw std::invalid_argument("within " + std::to_string(d) + " needs 'of', e.g. within 5 of sel");
      ++p;
      const auto from = factor();
      return select_within(s, from, d);
    }
    if (w == "x" || w == "y" || w == "z" || w == "charge" || w == "mass" || w == "bonds") {   // a range A..B, or one value
      if (t[p].kind != Tok::Number) throw std::invalid_argument(w + " needs a range, e.g. " + w + " 7..11");
      const std::string r = t[p++].text;
      double lo, hi;
      try {
        const auto dd = r.find("..", 1);
        lo = std::stod(r.substr(0, dd));
        hi = dd == std::string::npos ? lo : std::stod(r.substr(dd + 2));
      } catch (...) { throw std::invalid_argument("'" + r + "' is not a range (A..B)"); }
      if (hi < lo) std::swap(lo, hi);
      const auto col = atom_column(s, w);
      const double eps = w == "bonds" ? 0.5 : 1e-9;
      auto out = none();
      for (size_t i = 0; i < n; ++i) out[i] = col[i] >= lo - eps && col[i] <= hi + eps;
      return out;
    }
    if (w == "hybrid") {
      const std::string v = lower(t[p].kind == Tok::End ? "" : t[p++].text);
      const int want = v == "sp" ? 1 : v == "sp2" ? 2 : v == "sp3" ? 3 : -1;
      if (want < 0) throw std::invalid_argument("hybrid takes sp, sp2 or sp3");
      const auto col = atom_column(s, "hybrid");
      auto out = none();
      for (size_t i = 0; i < n; ++i) out[i] = int(col[i]) == want;
      return out;
    }
    if (w == "tag") {   // tag NAME: the atoms carrying that tag
      const std::string name = t[p++].text;
      int k = -1;
      for (size_t q = 0; q < s.tags.size(); ++q) if (s.tags[q].name == name) k = int(q);
      if (k < 0) throw std::invalid_argument("no tag '" + name + "' in this structure");
      auto out = none();
      for (size_t i = 0; i < n; ++i) out[i] = (s.atoms[i].tags >> k) & 1u;
      return out;
    }
    if (w == "sel" || w == "selection") {
      auto out = none();
      for (size_t i = 0; i < n && i < current.size(); ++i) out[i] = current[i];
      return out;
    }
    if (w == "all") return std::vector<char>(n, 1);
    if (w == "none") return none();
    throw std::invalid_argument("unknown word '" + t[p - 1].text + "' (smarts, element, type, chain, index, ring, stereo, within, tag, sel, and, or, not)");
  }
};

}  // namespace

QueryResult select_query(const System& s, const std::string& query, const std::vector<char>& current) {
  Parser ps{s, current, lex(query)};
  QueryResult r;
  r.atoms = ps.expr();
  if (ps.t[ps.p].kind != Tok::End) throw std::invalid_argument("unexpected '" + ps.t[ps.p].text + "' in the query");
  const auto& rings = ps.perception().rings;
  for (const auto& ring : rings)
    if (std::all_of(ring.begin(), ring.end(), [&](uint32_t a) { return r.atoms[a] != 0; })) ++r.rings;
  return r;
}

std::vector<char> cip_labels(const System& s) {
  const size_t n = s.atoms.size();
  std::vector<char> out(n, 0);
  const Perception per = perceive(s);
  const auto& nb = per.nb;
  // one sphere of a branch: the atomic numbers (duplicates for multiple bonds) reached from `front`, sorted high first
  struct Node { uint32_t atom; uint32_t from; bool duplicate; int z; };
  auto zof = [&](uint32_t a) { return s.atoms[a].element; };
  // the children of a node: its other neighbours, a duplicate for each extra bond order; an atom already reached in this
  // branch (a ring closure) and a duplicate carry three phantom substituents (Z 0), which have none
  std::vector<std::vector<char>> seen(4, std::vector<char>(n, 0));   // reused: only the touched entries are cleared
  std::vector<std::vector<uint32_t>> touched(4);
  auto children = [&](const Node& x, std::vector<Node>& kids, int branch) {
    if (x.duplicate) {
      if (x.z > 0) for (int k = 0; k < 3; ++k) kids.push_back({x.atom, x.atom, true, 0});
      return;
    }
    for (size_t k = 0; k < nb[x.atom].size(); ++k) {
      const uint32_t b = nb[x.atom][k];
      if (b == x.from) {
        // a multiple bond back to where we came from still adds its duplicates
        for (int e = 1; e < per.order[x.atom][k]; ++e) kids.push_back({b, x.atom, true, zof(b)});
        continue;
      }
      const bool again = seen[size_t(branch)][b] != 0;
      if (!again) { seen[size_t(branch)][b] = 1; touched[size_t(branch)].push_back(b); }
      kids.push_back({b, x.atom, again, zof(b)});
      for (int e = 1; e < per.order[x.atom][k]; ++e) kids.push_back({b, x.atom, true, zof(b)});
    }
  };
  for (uint32_t c = 0; c < n; ++c) {
    if (nb[c].size() != 4 || per.aromatic[c]) continue;
    const int z = zof(c);
    if (z != 6 && z != 7 && z != 14 && z != 15 && z != 16) continue;
    // rank the four branches sphere by sphere
    std::vector<std::vector<Node>> front(4);
    for (int k = 0; k < 4; ++k) {
      for (uint32_t a : touched[size_t(k)]) seen[size_t(k)][a] = 0;
      touched[size_t(k)] = {c, nb[c][size_t(k)]};
      front[size_t(k)] = {{nb[c][size_t(k)], c, false, zof(nb[c][size_t(k)])}};
      seen[size_t(k)][c] = seen[size_t(k)][nb[c][size_t(k)]] = 1;
    }
    std::vector<int> rank(4, -1);
    std::vector<std::vector<int>> key(4);
    bool resolved = false, exhausted = false;
    for (int sphere = 0; sphere < 16 && !resolved; ++sphere) {
      for (int k = 0; k < 4; ++k) {
        std::vector<int> zs;
        for (const auto& x : front[size_t(k)]) zs.push_back(x.z);
        std::sort(zs.rbegin(), zs.rend());
        key[size_t(k)].push_back(-1);   // sphere separator
        key[size_t(k)].insert(key[size_t(k)].end(), zs.begin(), zs.end());
      }
      std::vector<int> order = {0, 1, 2, 3};
      std::sort(order.begin(), order.end(), [&](int a, int b) { return key[size_t(a)] > key[size_t(b)]; });
      bool distinct = true;
      for (int k = 0; k + 1 < 4; ++k) if (key[size_t(order[size_t(k)])] == key[size_t(order[size_t(k + 1)])]) distinct = false;
      if (distinct) {
        for (int k = 0; k < 4; ++k) rank[size_t(order[size_t(k)])] = k;   // 0 highest
        resolved = true;
        break;
      }
      // the next sphere, explored from the atoms of this one in their order
      for (int k = 0; k < 4; ++k) {
        std::vector<Node> next;
        auto cur = front[size_t(k)];
        std::stable_sort(cur.begin(), cur.end(), [](const Node& a, const Node& b) { return a.z > b.z; });
        for (const auto& x : cur) children(x, next, k);
        front[size_t(k)] = std::move(next);
      }
      // every branch fully explored and still a tie: two identical substituents, not a stereocentre
      if (std::all_of(front.begin(), front.end(), [](const std::vector<Node>& f) {
            return std::all_of(f.begin(), f.end(), [](const Node& x) { return x.duplicate && x.z == 0; }) || f.empty();
          })) { exhausted = true; break; }
    }
    // four different neighbours in the first sphere, or distinct after exploration: a stereocentre
    const bool h2 = std::count_if(nb[c].begin(), nb[c].end(), [&](uint32_t b) { return zof(b) == 1; }) >= 2;
    if (h2) continue;
    if (!resolved) { if (!exhausted) out[c] = '*'; continue; }
    Vec3 v[4];
    for (int k = 0; k < 4; ++k) v[rank[size_t(k)]] = s.atoms[nb[c][size_t(k)]].pos - s.atoms[c].pos;
    if (s.cell.valid()) for (auto& x : v) x = s.cell.minimum_image(x);
    // looking with the lowest priority (v[3]) away: a → b → c clockwise is R
    const double chir = dot(v[0], cross(v[1], v[2]));
    out[c] = chir < 0 ? 'R' : 'S';
  }
  return out;
}

std::vector<double> atom_column(const System& s, const std::string& name) {
  const size_t n = s.atoms.size();
  std::vector<double> v(n, 0.0);
  if (name == "x" || name == "y" || name == "z") {
    const int k = name == "x" ? 0 : name == "y" ? 1 : 2;
    for (size_t i = 0; i < n; ++i) v[i] = s.atoms[i].pos[k];
  } else if (name == "charge") {
    for (size_t i = 0; i < n; ++i) v[i] = s.atoms[i].charge;
  } else if (name == "mass") {
    for (size_t i = 0; i < n; ++i) v[i] = s.mass_of(s.atoms[i]);
  } else if (name == "molecule") {
    const auto mol = s.molecules();
    for (size_t i = 0; i < n; ++i) v[i] = s.has_mol ? double(s.atoms[i].mol) : double(mol[i] + 1);
  } else if (name == "bonds" || name == "hybrid") {
    std::vector<int> deg(n, 0), doubles(n, 0), triples(n, 0), aromatic(n, 0);
    for (const auto& b : s.bonds) {
      if (b.i >= n || b.j >= n) continue;
      ++deg[b.i], ++deg[b.j];
      if (b.order == 2) ++doubles[b.i], ++doubles[b.j];
      if (b.order == 3) ++triples[b.i], ++triples[b.j];
      if (b.order == 4) ++aromatic[b.i], ++aromatic[b.j];
    }
    for (size_t i = 0; i < n; ++i) {
      if (name == "bonds") { v[i] = deg[i]; continue; }
      const int z = s.atoms[i].element, d = deg[i];
      int h = 0;
      if (triples[i] > 0 || doubles[i] >= 2) h = 1;
      else if (doubles[i] > 0 || aromatic[i] > 0) h = 2;
      else if (z == 6 || z == 14 || z == 32) h = d == 4 ? 3 : d == 3 ? 2 : d == 2 ? 1 : 0;   // carbon group: 4 σ bonds sp3
      else if (z == 7 || z == 15) h = d == 3 || d == 4 ? 3 : d == 2 ? 2 : d == 1 ? 1 : 0;
      else if (z == 8 || z == 16) h = d == 2 ? 3 : d == 1 ? 2 : 0;
      else if (z == 5) h = d == 3 ? 2 : d == 4 ? 3 : 0;
      v[i] = h;
    }
  } else {
    throw std::invalid_argument("no atom column '" + name + "' (x, y, z, charge, mass, bonds, hybrid, molecule)");
  }
  return v;
}

}  // namespace caps
