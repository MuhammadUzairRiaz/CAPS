#include <gtest/gtest.h>

#include <cmath>
#include <map>
#include <string>
#include <vector>

#include "caps/polymer.hpp"

using namespace caps;

namespace {
ChainSpec spec(std::vector<std::string> smiles, Sequence s = Sequence::Homopolymer, int dp = 10) {
  ChainSpec c;
  for (auto& x : smiles) c.units.push_back({x, x});
  c.sequence = s;
  c.dp = dp;
  return c;
}
}  // namespace

TEST(Polymer, RepeatUnits) {
  const UnitInfo ps = repeat_unit_info("*CC(*)c1ccccc1");
  EXPECT_EQ(ps.formula, "C8H8");
  EXPECT_NEAR(ps.mass, 104.15, 0.01);
  EXPECT_EQ(ps.stereocentres, 1);
  EXPECT_EQ(repeat_unit_info("*CC*").formula, "C2H4");
  EXPECT_EQ(repeat_unit_info("*CC*").stereocentres, 0);
  EXPECT_EQ(repeat_unit_info("[*]OCCOC(=O)c1ccc(C(=O)[*])cc1").formula, "C10H8O4");   // PET
  EXPECT_THROW(repeat_unit_info("CC*"), std::exception);
  EXPECT_THROW(repeat_unit_info("*C=C*C"), std::exception);
}

TEST(Polymer, Sequences) {
  auto s = chain_sequence(spec({"*CC*", "*CC(*)C", "*CCO*"}, Sequence::Alternating, 7), 1);
  EXPECT_EQ(s, (std::vector<int>{0, 1, 2, 0, 1, 2, 0}));
  ChainSpec b = spec({"*CC*", "*CC(*)C"}, Sequence::Block, 8);
  b.blocks = {3, 2};
  EXPECT_EQ(chain_sequence(b, 1), (std::vector<int>{0, 0, 0, 1, 1, 0, 0, 0}));
  ChainSpec p = spec({"*CC*", "*CC(*)C"}, Sequence::Pattern, 6);
  p.pattern = "AAB";
  EXPECT_EQ(chain_sequence(p, 1), (std::vector<int>{0, 0, 1, 0, 0, 1}));
  ChainSpec r = spec({"*CC*", "*CC(*)C", "*CCO*"}, Sequence::Random, 3000);
  r.weights = {0.6, 0.3, 0.1};
  std::map<int, int> n;
  for (int k : chain_sequence(r, 5)) ++n[k];
  EXPECT_NEAR(n[0] / 3000.0, 0.6, 0.04);
  EXPECT_NEAR(n[2] / 3000.0, 0.1, 0.03);
  auto g = chain_sequence(spec({"*CC*", "*CC(*)C"}, Sequence::Gradient, 200), 3);
  int a0 = 0, a1 = 0;
  for (int i = 0; i < 50; ++i) a0 += g[size_t(i)] == 0, a1 += g[size_t(150 + i)] == 0;
  EXPECT_GT(a0, a1 + 20);
}

TEST(Polymer, ChainGraph) {
  const ChainSpec c = spec({"*CC(*)c1ccccc1"}, Sequence::Homopolymer, 5);
  const MolGraph g = chain_graph(c, chain_sequence(c, 1));
  EXPECT_EQ(molecule_info(g).formula, "C40H42");   // five units and two hydrogen ends
  const ChainSpec d = spec({"*O[Si](C)(C)*"}, Sequence::Homopolymer, 3);   // bracket ends still get their caps
  EXPECT_EQ(molecule_info(chain_graph(d, chain_sequence(d, 1))).formula, "C6H20O3Si3");
}

TEST(Polymer, GrowsCellsWithRoom) {
  GrowOptions o;
  o.chains = 3;
  o.density = 0.3;
  o.seed = 4;
  o.contact_scale = 0.7;   // methacrylate backbones (quaternary carbons) grow at reduced contact limits, then Relax
  const ChainSpec c = spec({"*CC(*)c1ccccc1", "*CC(*)(C)C(=O)OC"}, Sequence::Alternating, 10);
  GrowReport rep;
  const System s = grow_chains(c, o, &rep);
  EXPECT_EQ(s.atoms.size(), size_t(3 * (5 * 16 + 5 * 15 + 2)));
  EXPECT_NEAR(s.density(), 0.3, 0.005);
  // no atoms of different chains closer than 1.9 Å (minimum image)
  const double L = s.cell.a[0];
  double dmin = 1e9;
  for (size_t i = 0; i < s.atoms.size(); ++i)
    for (size_t j = i + 1; j < s.atoms.size(); ++j) {
      if (s.atoms[i].mol == s.atoms[j].mol) continue;
      Vec3 d = s.atoms[i].pos - s.atoms[j].pos;
      for (int k = 0; k < 3; ++k) d[k] -= L * std::round(d[k] / L);
      dmin = std::min(dmin, norm(d));
    }
  EXPECT_GT(dmin, 0.7 * 2.0 - 0.1);   // the scaled H–H limit, less the accepted overlap
}

TEST(Polymer, Tacticity) {
  // sign of each backbone CH centre, neighbours taken in chain order: previous CH2, next CH2, ring carbon, H
  auto signs = [](Tacticity t) {
    ChainSpec c = spec({"*CC(*)c1ccccc1"}, Sequence::Homopolymer, 12);
    c.tacticity = t;
    GrowOptions o;
    o.chains = 1;
    o.density = 0.05;
    const System s = grow_chains(c, o);
    std::vector<std::vector<int>> nb(s.atoms.size());
    for (const auto& b : s.bonds) nb[b.i].push_back(int(b.j)), nb[b.j].push_back(int(b.i));
    std::vector<int> out;
    // walk the chain: units are 16 atoms, head CH2 first then the CH (tree order)
    for (int u = 0; u + 1 < 12; ++u) {
      const int ch2 = 16 * u, ch = 16 * u + 1, next = 16 * (u + 1);
      int ring = -1, h = -1;
      for (int w : nb[size_t(ch)]) {
        if (w == ch2 || w == next) continue;
        (s.atoms[size_t(w)].element == 1 ? h : ring) = w;
      }
      const Vec3 c0 = s.atoms[size_t(ch)].pos;
      const double v = dot(s.atoms[size_t(next)].pos - c0, cross(s.atoms[size_t(ring)].pos - c0, s.atoms[size_t(h)].pos - c0)) +
                       0 * dot(s.atoms[size_t(ch2)].pos, s.atoms[size_t(ch2)].pos);
      out.push_back(v > 0 ? 1 : -1);
    }
    return out;
  };
  const auto iso = signs(Tacticity::Isotactic), syn = signs(Tacticity::Syndiotactic);
  for (size_t k = 1; k < iso.size(); ++k) EXPECT_EQ(iso[k], iso[0]) << "isotactic unit " << k;
  for (size_t k = 1; k < syn.size(); ++k) EXPECT_EQ(syn[k], -syn[k - 1]) << "syndiotactic unit " << k;
}

TEST(Polymer, RubbersGrowAtFullContactLimits) {
  // cis-1,4-polyisoprene and a random SBR at 0.5 g/cm³
  for (auto units : {std::vector<std::string>{"[*]C/C=C(C)\\C[*]"}, std::vector<std::string>{"[*]C/C=C\\C[*]", "*CC(*)c1ccccc1"}}) {
    ChainSpec c = spec(units, units.size() > 1 ? Sequence::Random : Sequence::Homopolymer, 20);
    c.weights = {0.86, 0.14};
    GrowOptions o;
    o.chains = 5;
    o.density = 0.5;
    const System s = grow_chains(c, o);
    EXPECT_NEAR(s.density(), 0.5, 0.01);
  }
}

TEST(Polymer, BlendsOfTwoRubbers) {
  BlendComponent nr, br;
  nr.spec.units = {{"NR", "[*]C/C=C(C)\\C[*]"}};
  nr.spec.dp = 10;
  nr.weight = 0.7;
  br.spec.units = {{"BR", "[*]C/C=C\\C[*]"}};
  br.spec.dp = 10;
  br.weight = 0.3;
  BlendOptions o;
  o.chains = 6;
  BlendReport r;
  const System s = grow_blend({nr, br}, o, &r);
  ASSERT_EQ(r.chains.size(), 2u);
  EXPECT_EQ(r.chains[0], 6);
  EXPECT_NEAR(r.weight_fraction[0], 0.7, 0.08);   // whole chains: close to the asked fractions
  EXPECT_EQ(r.molecules[1].first, r.molecules[0].second + 1);
  int64_t top = 0;
  for (const auto& a : s.atoms) top = std::max(top, a.mol);
  EXPECT_EQ(top, r.chains[0] + r.chains[1]);
  // two slabs: each component in its own half of the cell along z
  o.morphology = BlendMorphology::Slabs;
  const System t = grow_blend({nr, br}, o, &r);
  const double half = norm(t.cell.c) / 2;
  for (const auto& a : t.atoms) {
    if (a.element == 1) continue;
    if (a.mol <= r.molecules[0].second) EXPECT_LT(a.pos[2], half);
    else EXPECT_GT(a.pos[2], half);
  }
}

TEST(Polymer, BlendDroplet) {
  BlendComponent nr, ps;
  nr.spec.units = {{"NR", "[*]C/C=C(C)\\C[*]"}};
  nr.spec.dp = 8;
  nr.weight = 0.8;
  ps.spec.units = {{"PS", "*CC(*)c1ccccc1"}};
  ps.spec.dp = 8;
  ps.weight = 0.2;
  BlendOptions o;
  o.chains = 6;
  o.morphology = BlendMorphology::Droplet;
  BlendReport r;
  const System s = grow_blend({nr, ps}, o, &r);
  // the minor component (PS) sits inside the droplet at the centre, the rubber outside it
  const Vec3 c = (s.cell.a + s.cell.b + s.cell.c) * 0.5;
  double rin = 0, rout = 1e9;
  for (const auto& a : s.atoms) {
    if (a.element == 1) continue;
    const double d = norm(s.cell.minimum_image(a.pos - c));
    if (a.mol >= r.molecules[1].first && a.mol <= r.molecules[1].second) rin = std::max(rin, d);
    else rout = std::min(rout, d);
  }
  EXPECT_LT(rin, rout + 0.5);
}
