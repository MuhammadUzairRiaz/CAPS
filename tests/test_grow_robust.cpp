// The chain grower on hard requests: dense melts, many short chains, single units, a long single chain, rigid ring
// backbones and a cross-section of the polymer library at 1.1 g/cm³, every cell checked for its atoms, chains, units,
// density, bonds and closest contacts. bench/perf/grow_sweep.py runs the whole library on a larger grid.
#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "caps/elements.hpp"
#include "caps/json.hpp"
#include "caps/polymer.hpp"

using namespace caps;

namespace {

ChainSpec units_spec(const std::vector<std::string>& smiles, int dp, Sequence seq = Sequence::Homopolymer) {
  ChainSpec c;
  for (size_t k = 0; k < smiles.size(); ++k) c.units.push_back({std::string("U") + char('A' + int(k)), smiles[k]});
  c.dp = dp;
  c.sequence = seq;
  return c;
}

// Closest distance between atoms more than three bonds apart (minimum image, orthorhombic cell), by a cell list.
double closest_nonbonded(const System& s) {
  const size_t n = s.atoms.size();
  std::vector<std::vector<uint32_t>> nb(n);
  for (const auto& b : s.bonds) nb[b.i].push_back(b.j), nb[b.j].push_back(b.i);
  const double L[3] = {s.cell.a[0], s.cell.b[1], s.cell.c[2]};
  int nc[3];
  for (int k = 0; k < 3; ++k) nc[k] = std::max(1, int(L[k] / 3.0));
  auto bin = [&](const Vec3& p, int k) { return ((int(std::floor(p[k] / (L[k] / nc[k]))) % nc[k]) + nc[k]) % nc[k]; };
  std::map<std::array<int, 3>, std::vector<uint32_t>> bins;
  for (uint32_t i = 0; i < n; ++i) bins[{bin(s.atoms[i].pos, 0), bin(s.atoms[i].pos, 1), bin(s.atoms[i].pos, 2)}].push_back(i);
  double best = std::min({L[0], L[1], L[2]});
  std::vector<int> mark(n, -1);
  for (uint32_t i = 0; i < n; ++i) {
    std::vector<uint32_t> fr{i};
    mark[i] = int(i);
    for (int d = 0; d < 3; ++d) {
      std::vector<uint32_t> nx;
      for (uint32_t u : fr)
        for (uint32_t w : nb[u])
          if (mark[w] != int(i)) mark[w] = int(i), nx.push_back(w);
      fr = std::move(nx);
    }
    const int bx = bin(s.atoms[i].pos, 0), by = bin(s.atoms[i].pos, 1), bz = bin(s.atoms[i].pos, 2);
    std::set<std::array<int, 3>> seen;
    for (int dx = -1; dx <= 1; ++dx)
      for (int dy = -1; dy <= 1; ++dy)
        for (int dz = -1; dz <= 1; ++dz) {
          const std::array<int, 3> key{(bx + dx + nc[0]) % nc[0], (by + dy + nc[1]) % nc[1], (bz + dz + nc[2]) % nc[2]};
          if (!seen.insert(key).second) continue;
          auto it = bins.find(key);
          if (it == bins.end()) continue;
          for (uint32_t j : it->second) {
            if (j <= i || mark[j] == int(i)) continue;
            Vec3 v = s.atoms[j].pos - s.atoms[i].pos;
            for (int k = 0; k < 3; ++k) v[k] -= L[k] * std::round(v[k] / L[k]);
            best = std::min(best, norm(v));
          }
        }
  }
  return best;
}

// Everything a grown cell must be: the chains and units asked for, the atoms of its units plus two end caps per chain,
// the density, every bond a bond (unwrapped, none across the cell), no two atoms more than three bonds apart closer
// than `contact` Å (eased limits go down to 0.4 of the full ones: H–H 0.78 Å, four bonds apart 85 % of that).
void expect_cell(const System& s, const ChainSpec& c, int chains, double density, double contact, const std::string& what) {
  std::map<int64_t, std::vector<size_t>> mols;
  for (size_t i = 0; i < s.atoms.size(); ++i) mols[s.atoms[i].mol].push_back(i);
  EXPECT_EQ(int(mols.size()), chains) << what;
  std::map<std::string, int> unit_atoms;
  for (const auto& u : c.units) unit_atoms[u.name] = repeat_unit_info(u.smiles).atoms;
  for (const auto& [m, idx] : mols) {
    std::map<int64_t, std::string> res;
    for (size_t i : idx) res[s.atoms[i].resid] = s.atoms[i].resname;
    EXPECT_EQ(int(res.size()), c.dp) << what << " molecule " << m;
    int want = 2;
    for (const auto& [r, name] : res) want += unit_atoms.count(name) ? unit_atoms[name] : 0;
    EXPECT_EQ(int(idx.size()), want) << what << " molecule " << m;
  }
  EXPECT_NEAR(s.density(), density, 0.01 * density) << what;
  int bad = 0;
  for (const auto& b : s.bonds) {
    const double d = norm(s.atoms[b.i].pos - s.atoms[b.j].pos);
    bad += d < 0.8 || d > 2.4;
  }
  EXPECT_EQ(bad, 0) << what << ": bonds outside 0.8–2.4 Å";
  EXPECT_GE(closest_nonbonded(s), contact) << what;
}

// Worst (distance − limit) over pairs more than three bonds apart, at the grower's own full limits: 0.88 × the sum of
// the van der Waals radii less 0.08 Å per hydrogen (H–H 1.95, C–H 2.47, C–C 2.99 Å), 85 % for pairs four bonds apart,
// 90 % where an end-cap hydrogen (the last two atoms of each chain) is one of the two. Minimum image, orthorhombic cell.
double worst_margin_at_full_limits(const System& s) {
  const size_t n = s.atoms.size();
  std::vector<char> cap(n, 0);
  for (size_t i = 1; i < n; ++i)
    if (i + 1 == n || s.atoms[i + 1].mol != s.atoms[i].mol) cap[i] = cap[i - 1] = 1;
  std::vector<std::vector<uint32_t>> nb(n);
  for (const auto& b : s.bonds) nb[b.i].push_back(b.j), nb[b.j].push_back(b.i);
  const double L[3] = {s.cell.a[0], s.cell.b[1], s.cell.c[2]};
  double worst = 1e9;
  std::vector<int> d(n, -1);
  for (uint32_t i = 0; i < n; ++i) {
    std::fill(d.begin(), d.end(), -1);
    d[i] = 0;
    std::vector<uint32_t> q{i};
    for (size_t h = 0; h < q.size(); ++h)
      if (d[q[h]] < 4)
        for (uint32_t w : nb[q[h]])
          if (d[w] < 0) d[w] = d[q[h]] + 1, q.push_back(w);
    for (uint32_t j = i + 1; j < n; ++j) {
      if (d[j] >= 0 && d[j] <= 3) continue;
      Vec3 r = s.atoms[j].pos - s.atoms[i].pos;
      for (int k = 0; k < 3; ++k) r[k] -= L[k] * std::round(r[k] / L[k]);
      const int zi = s.atoms[i].element, zj = s.atoms[j].element;
      const double lim = 0.88 * (element(zi).vdw + element(zj).vdw) - 0.08 * ((zi == 1) + (zj == 1));
      worst = std::min(worst, norm(r) - lim * (d[j] == 4 ? 0.85 : 1.0) * (cap[i] || cap[j] ? 0.9 : 1.0));
    }
  }
  return worst;
}

}  // namespace

// A PE melt above its real density: chains stuck in the crowd ease their own contact limits (not below 0.4) and the
// cell is complete; at fixed limits the same request keeps the full limits everywhere
TEST(GrowRobust, DenseMeltWithEasedLimits) {
  const ChainSpec c = units_spec({"*CC*"}, 60);
  GrowOptions o;
  o.chains = 8;
  o.density = 1.0;
  o.auto_scale = true;
  GrowReport rep;
  const System s = grow_chains(c, o, &rep);
  expect_cell(s, c, 8, 1.0, 0.6, "PE 8 x 60 at 1.0");
  EXPECT_GT(rep.eased_steps, 0);
  EXPECT_GE(rep.min_scale, 0.4 - 1e-9);
  EXPECT_GE(rep.worst_margin, o.accept - 1e-9);   // against the limits each unit was placed at

  // several seeds: one growth path can pass or fail on rounding alone (each platform's maths library takes its own)
  for (uint64_t seed = 1; seed <= 10; ++seed) {
    GrowOptions f;
    f.chains = 4;
    f.density = 0.4;
    f.seed = seed;
    GrowReport fr;
    const System fixed = grow_chains(units_spec({"*CC*"}, 30), f, &fr);
    EXPECT_EQ(fr.eased_steps, 0) << "seed " << seed;
    EXPECT_DOUBLE_EQ(fr.min_scale, 1.0) << "seed " << seed;
    EXPECT_GE(fr.worst_margin, f.accept - 1e-9) << "seed " << seed;
    // the cell as written, end caps included, against the full limits: within the accepted 0.05 Å
    EXPECT_GE(worst_margin_at_full_limits(fixed), f.accept - 1e-9) << "seed " << seed;
  }
}

// Many short chains, one unit, one long chain: the start points and the first units are the hard part
TEST(GrowRobust, ShortChainsSingleUnitsAndOneLongChain) {
  struct Case { int chains, dp; double density; };
  for (const Case& k : {Case{200, 2, 0.9}, Case{50, 5, 1.1}, Case{1, 1, 0.3}, Case{3, 1, 1.1}, Case{1, 150, 1.0}}) {
    const ChainSpec c = units_spec({"*CC(*)C"}, k.dp);
    GrowOptions o;
    o.chains = k.chains;
    o.density = k.density;
    o.auto_scale = true;
    const System s = grow_chains(c, o);
    expect_cell(s, c, k.chains, k.density, 0.6, "PP " + std::to_string(k.chains) + " x " + std::to_string(k.dp) + " at " + std::to_string(k.density));
  }
}

// Ethylidene norbornene (2,3-linked bicyclic units): 12 chains of 10 units at 0.3 g/cm³ took over ten minutes; the
// probe chain finds the scale the rigid units allow and the stuck chains ease below it
TEST(GrowRobust, RingBackboneCell) {
  const ChainSpec c = units_spec({"[*][C@H]1[C@@H]2CC(=CC)[C@@H](C2)[C@H]1[*]"}, 10);
  GrowOptions o;
  o.chains = 12;
  o.density = 0.3;
  o.seed = 7;
  o.auto_scale = true;
  GrowReport rep;
  const System s = grow_chains(c, o, &rep);
  expect_cell(s, c, 12, 0.3, 0.6, "ENB 12 x 10");
  EXPECT_GE(rep.min_scale, 0.4 - 1e-9);
}

// A cross-section of the library, homopolymers and copolymers, as dense cells at 1.1 g/cm³ by best-of-k trials and at
// 0.9 g/cm³ by Rosenbluth growth
TEST(GrowRobust, LibraryCrossSectionDense) {
  std::ifstream f(std::string(CAPS_SOURCE_DIR) + "/data/polymers/library.json");
  std::stringstream ss;
  ss << f.rdbuf();
  const Json lib = Json::parse(ss.str());
  std::map<std::string, std::string> smiles;
  for (const auto& p : lib["polymers"].items()) smiles[p["id"].str()] = p["smiles"].str();
  struct Case { std::vector<std::string> ids; Sequence seq; int method; };
  const std::vector<Case> cases = {
      {{"P002"}, Sequence::Homopolymer, 0}, {{"P013"}, Sequence::Homopolymer, 0}, {{"P019"}, Sequence::Homopolymer, 0},
      {{"P056"}, Sequence::Homopolymer, 0}, {{"P202"}, Sequence::Homopolymer, 0}, {{"P201"}, Sequence::Homopolymer, 0},
      {{"P017"}, Sequence::Homopolymer, 0}, {{"P045"}, Sequence::Homopolymer, 0}, {{"P072"}, Sequence::Homopolymer, 0},
      {{"P084"}, Sequence::Homopolymer, 0}, {{"P013", "P005"}, Sequence::Random, 0}, {{"P013", "P019"}, Sequence::Block, 0},
      {{"P001"}, Sequence::Homopolymer, 1}, {{"P013"}, Sequence::Homopolymer, 2}, {{"P019"}, Sequence::Homopolymer, 1}};
  for (const auto& k : cases) {
    std::vector<std::string> u;
    std::string what;
    for (const auto& id : k.ids) u.push_back(smiles.at(id)), what += id + " ";
    ChainSpec c = units_spec(u, 8, k.seq);
    if (k.seq == Sequence::Random) c.weights = {0.25, 0.75};
    if (k.seq == Sequence::Block) c.blocks = {4, 4};
    GrowOptions o;
    o.chains = 10;
    o.density = k.method == 0 ? 1.1 : 0.9;
    o.method = k.method;
    o.seed = 5;
    o.auto_scale = true;
    what += "method " + std::to_string(k.method);
    System s;
    try {
      s = grow_chains(c, o);
    } catch (const std::exception& e) {
      ADD_FAILURE() << what << ": " << e.what();
      continue;
    }
    expect_cell(s, c, 10, o.density, 0.6, what);
  }
}

// The Rosenbluth trials scored in parallel: the same draws, the same weights, the same cell as one at a time
TEST(GrowRobust, RosenbluthParallelEqualsSerial) {
  for (int method : {1, 2}) {
    std::vector<Vec3> pos[2];
    for (int t = 0; t < 2; ++t) {
      GrowOptions o;
      o.chains = 3;
      o.density = 0.4;
      o.method = method;
      o.seed = 11;
      o.threads = t == 0 ? 1 : 4;
      const System s = grow_chains(units_spec({"*CC(*)c1ccccc1"}, 10), o);
      for (const auto& a : s.atoms) pos[t].push_back(a.pos);
    }
    ASSERT_EQ(pos[0].size(), pos[1].size());
    for (size_t i = 0; i < pos[0].size(); ++i) EXPECT_EQ(pos[0][i], pos[1][i]) << "method " << method << " atom " << i;
  }
}
