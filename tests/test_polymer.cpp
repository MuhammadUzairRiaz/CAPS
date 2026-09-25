#include <gtest/gtest.h>

#include <cmath>
#include <map>
#include <string>
#include <vector>

#include "caps/analysis.hpp"
#include "caps/edit.hpp"
#include "caps/polymer.hpp"
#include "caps/query.hpp"

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

namespace {
// every carbon with valence four (aromatic bonds 1.5), every hydrogen with one bond, one molecule per star or comb, and no two atoms more than
// three bonds apart closer than the scaled H–H limit
void expect_sound(const System& s, int molecules, double scale) {
  std::vector<std::vector<uint32_t>> nb(s.atoms.size());
  std::vector<double> valence(s.atoms.size(), 0.0);
  for (const auto& b : s.bonds) {
    nb[b.i].push_back(b.j), nb[b.j].push_back(b.i);
    const double v = b.order == 4 ? 1.5 : b.order == 2 ? 2 : b.order == 3 ? 3 : 1;
    valence[b.i] += v, valence[b.j] += v;
  }
  for (size_t i = 0; i < s.atoms.size(); ++i) EXPECT_NEAR(valence[i], s.atoms[i].element == 6 ? 4 : 1, 1e-9) << "atom " << i;
  int nm = 0;
  s.molecules(&nm);
  EXPECT_EQ(nm, molecules);
  const double L = s.cell.a[0];
  double dmin = 1e9;
  for (size_t i = 0; i < s.atoms.size(); ++i) {
    std::map<uint32_t, int> d{{uint32_t(i), 0}};
    std::vector<uint32_t> q{uint32_t(i)};
    for (size_t k = 0; k < q.size(); ++k)
      if (d[q[k]] < 3)
        for (uint32_t w : nb[q[k]])
          if (d.emplace(w, d[q[k]] + 1).second) q.push_back(w);
    for (size_t j = i + 1; j < s.atoms.size(); ++j) {
      if (d.count(uint32_t(j))) continue;
      Vec3 v = s.atoms[i].pos - s.atoms[j].pos;
      for (int k = 0; k < 3; ++k) v[k] -= L * std::round(v[k] / L);
      dmin = std::min(dmin, norm(v));
    }
  }
  EXPECT_GT(dmin, scale * 2.0 - 0.2);
}
}  // namespace

TEST(Polymer, StarsOfFourArmsOnOneCarbon) {
  GrowOptions o;
  o.chains = 3;
  o.density = 0.3;
  o.seed = 2;
  o.auto_scale = true;   // branch points are crowded: the contact scale steps down when needed, as in the Studio
  ChainSpec c = spec({"*CC*"}, Sequence::Homopolymer, 10);
  c.architecture = Architecture::Star;
  c.arms = 4;
  GrowReport rep;
  const System s = grow_chains(c, o, &rep);
  // 4 arms × 10 C₂H₄ units, a tail cap on each arm, and the core's head cap and two hydrogens given to arms
  EXPECT_EQ(s.atoms.size(), size_t(3 * (4 * 10 * 6 + 4 - 2)));
  expect_sound(s, 3, 0.6);
  // each core carbon is bonded to four carbons
  std::vector<int> carbons(s.atoms.size(), 0);
  for (const auto& b : s.bonds)
    if (s.atoms[b.i].element == 6 && s.atoms[b.j].element == 6) ++carbons[b.i], ++carbons[b.j];
  EXPECT_EQ(std::count(carbons.begin(), carbons.end(), 4), 3);
  EXPECT_NE(rep.notes.front().find("3 stars of 4 arms"), std::string::npos) << rep.notes.front();
}

TEST(Polymer, CombsAndRandomBranches) {
  GrowOptions o;
  o.chains = 2;
  o.density = 0.3;
  o.seed = 5;
  o.auto_scale = true;
  ChainSpec c = spec({"*CC(*)c1ccccc1"}, Sequence::Homopolymer, 12);
  c.architecture = Architecture::Comb;
  c.spacing = 4;
  c.arm_dp = 3;
  const System s = grow_chains(c, o);
  // side chains on units 4, 8, 12: 12 + 9 styrene units per comb, caps on both backbone ends and each side chain's
  // tail, one backbone hydrogen per side chain given up
  EXPECT_EQ(s.atoms.size(), size_t(2 * (21 * 16 + 2 + 3 - 3)));
  expect_sound(s, 2, 0.6);
  c.architecture = Architecture::Branched;
  c.branch_probability = 0.3;
  GrowReport rep;
  const System b = grow_chains(c, o, &rep);
  expect_sound(b, 2, 0.6);
  EXPECT_NE(rep.notes.front().find("branched chains"), std::string::npos) << rep.notes.front();
}

// Region shapes: a slab with vacuum (heavy atoms inside the film heights) and chains around a cylinder along z (none
// inside its radius); the cell edges follow from the density of the space the chains may use
TEST(Polymer, SlabAndCylinderRegions) {
  ChainSpec c = spec({"*CC=CC*"}, Sequence::Homopolymer, 15);
  for (int shape = 0; shape < 3; ++shape) {
    GrowOptions o;
    o.chains = 6;
    o.density = 0.6;
    o.seed = 3;
    o.auto_scale = true;
    if (shape == 0) o.slab_thickness = 20, o.slab_vacuum = 30;
    if (shape == 1) o.cylinder_radius = 11;
    if (shape == 2) o.cylinder_radius = 8, o.cylinder_outside = true;
    const System s = grow_chains(c, o);
    const double Lx = s.cell.a[0], Lz = s.cell.c[2];
    for (const auto& a : s.atoms) {
      if (a.element == 1) continue;   // end caps are placed after growth
      const double z = a.pos[2] - Lz * std::floor(a.pos[2] / Lz);
      double dx = a.pos[0] - Lx / 2, dy = a.pos[1] - Lx / 2;
      dx -= Lx * std::round(dx / Lx), dy -= Lx * std::round(dy / Lx);
      const double r = std::hypot(dx, dy);
      if (shape == 0) { EXPECT_GE(z, 15.0 - 0.1); EXPECT_LE(z, 35.0 + 0.1); }
      if (shape == 1) EXPECT_LE(r, 11.0 + 0.1);
      if (shape == 2) EXPECT_GE(r, 8.0 - 0.1);
    }
    if (shape == 0) EXPECT_NEAR(Lz, 50.0, 1e-9);
    // the density over the space the chains may use is the target
    double mass = 0;
    for (const auto& a : s.atoms) mass += s.mass_of(a);
    const double room = shape == 0 ? Lx * s.cell.b[1] * 20 : shape == 1 ? M_PI * 121 * Lz : (Lx * Lx - M_PI * 64) * Lz;
    EXPECT_NEAR(mass / (0.602214076 * room), 0.6, 0.02) << "shape " << shape;
  }
}

// Growth methods: keeping the roomiest trial stretches polyethylene (trans-rich); Rosenbluth selection with butane
// torsions at 450 K gives the trans share of a Boltzmann chain (about 0.57 for three states) and reports ln W
TEST(Polymer, RosenbluthGrowthIsLessStretched) {
  ChainSpec c = spec({"*CC*"}, Sequence::Homopolymer, 40);
  double trans[3];
  for (int method = 0; method < 3; ++method) {
    GrowOptions o;
    o.chains = 6;
    o.density = 0.4;
    o.seed = 5;
    o.method = method;
    o.auto_scale = true;
    GrowReport r;
    System s = grow_chains(c, o, &r);
    int nt = 0, n = 0;
    for (const auto& b : backbones(s))
      for (size_t k = 0; k + 3 < b.size(); ++k, ++n) nt += std::fabs(measure(s, {b[k], b[k + 1], b[k + 2], b[k + 3]})) > 120;
    trans[method] = double(nt) / n;
    if (method > 0) {
      EXPECT_TRUE(std::isfinite(r.ln_rosenbluth) && r.ln_rosenbluth < 0) << r.ln_rosenbluth;
      EXPECT_NE(r.notes.back().find("Rosenbluth growth"), std::string::npos);
    }
  }
  EXPECT_GT(trans[0], 0.68);
  for (int m : {1, 2}) {
    EXPECT_LT(trans[m], trans[0] - 0.08) << "method " << m;
    EXPECT_GT(trans[m], 0.40) << "method " << m;
  }
}

// A DNA strand from nucleotide units (head on O5′, tail on P): every sugar keeps its D configuration (C4′ R, C1′ R, C3′ S)
// when configurations are kept; the 3′ P–H cap becomes a phosphate P–OH
TEST(Polymer, DnaStrandKeepsDSugars) {
  const char* bases[4] = {"n2cnc3c(N)ncnc32", "N2C=CC(N)=NC2=O", "N2C=NC3=C2N=C(N)NC3=O", "N2C=C(C)C(=O)NC2=O"};
  ChainSpec c;
  for (int b = 0; b < 4; ++b)
    c.units.push_back({std::string(1, "ACGT"[b]), std::string("*OC[C@H]1O[C@@H](") + bases[b] + ")C[C@@H]1OP(=O)([O-])*"});
  c.sequence = Sequence::Pattern;
  c.pattern = "ACGTTGCA";
  c.dp = 8;
  c.keep_configuration = true;
  GrowOptions o;
  o.chains = 1;
  o.density = 0.02;
  o.seed = 4;
  o.auto_scale = true;
  System s = grow_chains(c, o);
  EXPECT_EQ(hydroxylate_phosphorus(s), 1);
  const auto cip = cip_labels(s);
  int r = 0, sl = 0;
  for (size_t i = 0; i < cip.size(); ++i)   // the sugar carbons (phosphorus reads as a centre too: its =O and O− differ)
    if (s.atoms[i].element == 6) r += cip[i] == 'R', sl += cip[i] == 'S';
  EXPECT_EQ(r, 16);   // C4′ and C1′ of eight sugars
  EXPECT_EQ(sl, 8);   // C3′
}
