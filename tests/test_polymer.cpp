#include <cstdio>
#include <gtest/gtest.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "caps/cg_map.hpp"
#include "caps/relax.hpp"
#include "caps/elements.hpp"
#include "caps/analysis.hpp"
#include "caps/edit.hpp"
#include "caps/polymer.hpp"
#include "caps/molecule.hpp"
#include "caps/molinfo.hpp"
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
  // exact composition: PBSA BS:BA 80:20 at DP 25 holds 20 and 5 in every chain, in a random order
  ChainSpec x = spec({"[*]OCCCCOC(=O)CCC(=O)[*]", "[*]OCCCCOC(=O)CCCCC(=O)[*]"}, Sequence::Shuffled, 25);
  x.weights = {0.8, 0.2};
  std::set<std::vector<int>> orders;
  for (unsigned seed = 1; seed <= 5; ++seed) {
    const auto q = chain_sequence(x, seed);
    EXPECT_EQ(std::count(q.begin(), q.end(), 0), 20);
    EXPECT_EQ(std::count(q.begin(), q.end(), 1), 5);
    orders.insert(q);
  }
  EXPECT_GT(orders.size(), 1u);
  ChainSpec y = spec({"*CC*", "*CC(*)C", "*CCO*"}, Sequence::Shuffled, 10);   // 1/3 each of 10: largest remainders 4, 3, 3
  const auto qy = chain_sequence(y, 2);
  EXPECT_EQ(std::count(qy.begin(), qy.end(), 0) + std::count(qy.begin(), qy.end(), 1) + std::count(qy.begin(), qy.end(), 2), 10);
  EXPECT_EQ(std::count(qy.begin(), qy.end(), 0), 4);
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

// A dendrimer: a star core of 3 arms, every end splitting in two for 2 generations (3 + 6 + 12 segments)
TEST(Polymer, DendrimerGenerations) {
  GrowOptions o;
  o.chains = 2;
  o.density = 0.2;
  o.seed = 4;
  o.auto_scale = true;
  ChainSpec c = spec({"*CC*"}, Sequence::Homopolymer, 4);
  c.architecture = Architecture::Dendrimer;
  c.arms = 3;
  c.arm_dp = 3;
  c.generations = 2;
  GrowReport rep;
  const System s = grow_chains(c, o, &rep);
  // core: 3 arms × 4 C₂H₄ with 2 hydrogens net; 18 branches × 3 C₂H₄, each a tail cap for the hydrogen it replaced
  EXPECT_EQ(s.atoms.size(), size_t(2 * (3 * 4 * 6 + 2 + 18 * 3 * 6)));
  expect_sound(s, 2, 0.6);
  std::vector<int> heavy(s.atoms.size(), 0), all(s.atoms.size(), 0);
  for (const auto& b : s.bonds) {
    ++all[b.i], ++all[b.j];
    if (s.atoms[b.i].element == 6 && s.atoms[b.j].element == 6) ++heavy[b.i], ++heavy[b.j];
  }
  for (size_t i = 0; i < s.atoms.size(); ++i) EXPECT_EQ(all[i], s.atoms[i].element == 6 ? 4 : 1) << i;
  // chain ends: 12 outer branch ends per dendrimer, each a CH₃
  int ends = 0;
  for (size_t i = 0; i < s.atoms.size(); ++i) ends += s.atoms[i].element == 6 && heavy[i] == 1;
  EXPECT_EQ(ends, 2 * 12);
  EXPECT_NE(rep.notes.front().find("dendrimers of generation 2"), std::string::npos) << rep.notes.front();
}

// Head and tail on one ring (addition norbornenes): flagged until their configuration is written; exo,exo grows
TEST(Polymer, RingBackboneStereo) {
  const UnitInfo open = repeat_unit_info("[*]C1C2CCC(C2)C1[*]");
  EXPECT_TRUE(open.ring_backbone);
  EXPECT_TRUE(open.ring_stereo_open);
  const UnitInfo exo = repeat_unit_info("[*][C@H]1[C@@H]2CC[C@@H](C2)[C@H]1[*]");
  EXPECT_TRUE(exo.ring_backbone);
  EXPECT_FALSE(exo.ring_stereo_open);
  EXPECT_FALSE(repeat_unit_info("*CC(*)c1ccccc1").ring_backbone);
  EXPECT_FALSE(repeat_unit_info("[*]C1CCC(CC[*])C1").ring_stereo_open);   // ROMP: the tail is off the ring
  GrowOptions o;
  o.chains = 1;
  o.density = 0.05;
  o.seed = 1;
  o.auto_scale = true;
  const System s = grow_chains(spec({"[*][C@H]1[C@@H]2CC[C@@H](C2)[C@H]1[*]"}, Sequence::Homopolymer, 8), o);
  EXPECT_EQ(s.atoms.size(), size_t(8 * 17 + 2));   // C₇H₁₀ units, two end caps
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

// End groups: methyl heads and tert-butyl tails on polyethylene chains — each chain gains 1 + 4 carbons, loses no
// hydrogen count beyond the two caps replaced, and stays one molecule
TEST(Polymer, EndGroupsReplaceTheCaps) {
  GrowOptions o;
  o.chains = 2;
  o.density = 0.3;
  o.seed = 3;
  ChainSpec c = spec({"*CC*"}, Sequence::Homopolymer, 5);
  const System plain = grow_chains(c, o);
  c.head_cap = "methyl";
  c.tail_cap = "tert-butyl";
  GrowReport rep;
  const System capped = grow_chains(c, o, &rep);
  auto count = [](const System& s, int z) { int n = 0; for (const auto& a : s.atoms) n += a.element == z; return n; };
  EXPECT_EQ(count(capped, 6), count(plain, 6) + 2 * 5);
  EXPECT_EQ(count(capped, 1), count(plain, 1) - 4 + 2 * (3 + 9));
  int nm = 0;
  capped.molecules(&nm);
  EXPECT_EQ(nm, 2);
  bool said = false;
  for (const auto& n : rep.notes) said |= n.find("chain ends capped") != std::string::npos;
  EXPECT_TRUE(said);
  EXPECT_EQ(chain_end_smiles("hydrogen"), "");
  ChainSpec bad = c;
  bad.head_cap = "nonsense";
  EXPECT_THROW(grow_chains(bad, o), std::invalid_argument);
}

TEST(Polymer, HeadToHeadLinkageReversesEverySecondUnit) {
  GrowOptions o;
  o.chains = 2;
  o.density = 0.3;
  o.seed = 5;
  ChainSpec c = spec({"*CC(*)Cl"}, Sequence::Homopolymer, 6);
  // bonds between two chlorinated carbons: none head-to-tail, three per chain head-to-head (units 1–2, 3–4, 5–6)
  auto hh_bonds = [](const System& s) {
    std::vector<int> cl(s.atoms.size(), 0);
    for (const auto& b : s.bonds) {
      if (s.atoms[b.i].element == 17) cl[b.j] = 1;
      if (s.atoms[b.j].element == 17) cl[b.i] = 1;
    }
    int n = 0;
    for (const auto& b : s.bonds) n += cl[b.i] && cl[b.j];
    return n;
  };
  EXPECT_EQ(hh_bonds(grow_chains(c, o)), 0);
  c.linkage = Linkage::HeadToHead;
  GrowReport rep;
  const System hh = grow_chains(c, o, &rep);
  EXPECT_EQ(hh_bonds(hh), 2 * 3);
  bool said = false;
  for (const auto& n : rep.notes) said |= n.find("head-to-head linkage: 6 of 12") != std::string::npos;
  EXPECT_TRUE(said);
  const auto inv = chain_inversions(c, 6, 1);
  const std::string smi = chain_graph(c, {0, 0, 0, 0, 0, 0}, inv).smiles;
  EXPECT_NE(smi, chain_graph(c, {0, 0, 0, 0, 0, 0}).smiles) << smi;
  EXPECT_TRUE(smi.find("C(Cl)C(CC") != std::string::npos) << smi;
  // random: a share of units reversed, the first never
  c.linkage = Linkage::Random;
  c.inversion = 0.5;
  const auto r = chain_inversions(c, 400, 9);
  int n = 0;
  for (char x : r) n += x;
  EXPECT_EQ(r[0], 0);
  EXPECT_NEAR(n / 400.0, 0.5, 0.08);
  EXPECT_EQ(linkage_from_string("head-to-head"), Linkage::HeadToHead);
  EXPECT_THROW(linkage_from_string("sideways"), std::invalid_argument);
}

TEST(Polymer, OrientedGrowthAlignsTheBackbone) {
  GrowOptions o;
  o.chains = 4;
  o.density = 0.3;
  o.seed = 2;
  const ChainSpec c = spec({"*CC*"}, Sequence::Homopolymer, 30);
  GrowReport iso, best, ros;
  grow_chains(c, o, &iso);
  EXPECT_LT(std::fabs(iso.orientation), 0.2);   // isotropic: ⟨P₂⟩ about 0 (z for the report)
  o.orient_axis = {0, 0, 1};
  o.orient_strength = 4;
  grow_chains(c, o, &best);
  EXPECT_GT(best.orientation, 0.5) << "best-of-k with the field";
  o.method = 1;
  grow_chains(c, o, &ros);
  EXPECT_GT(ros.orientation, 0.4) << "Rosenbluth with the field";
  bool said = false;
  for (const auto& n : ros.notes) said |= n.find("oriented growth") != std::string::npos;
  EXPECT_TRUE(said);
  std::printf("orientation: isotropic %.3f · best-of-k %.3f · Rosenbluth %.3f\n", iso.orientation, best.orientation, ros.orientation);
}

TEST(Polymer, DeeperLookAheadGrows) {
  GrowOptions o;
  o.chains = 6;
  o.density = 0.6;
  o.seed = 4;
  const ChainSpec c = spec({"*CC(*)c1ccccc1"}, Sequence::Homopolymer, 12);
  // look-ahead 3 grows the same cell within the acceptance margin; whether it backtracks less depends on the path
  // (one growth is chaotic: it differs between math libraries in the last bits, and the counts went both ways on
  // different machines), so the counts are reported, not compared
  int b1 = 0, b3 = 0;
  for (uint64_t seed : {4, 5}) {
    GrowReport r1, r3;
    o.seed = seed;
    o.lookahead = 1;
    const System a = grow_chains(c, o, &r1);
    o.lookahead = 3;
    const System b = grow_chains(c, o, &r3);
    EXPECT_EQ(a.atoms.size(), b.atoms.size());
    EXPECT_GE(r3.worst_margin, o.accept);
    b1 += r1.backtracks, b3 += r3.backtracks;
  }
  std::printf("look-ahead 1: %d backtracks · look-ahead 3: %d backtracks (two seeds)\n", b1, b3);
}

TEST(Polymer, LogNormalAndHistogramChainLengths) {
  // log-normal: the sample's Nn and Đ come back (continuous draws, rounded to integers)
  const auto n = draw_chain_lengths("log-normal", 60, 1.4, 40000, 7);
  double s1 = 0, s2 = 0;
  for (int x : n) s1 += x, s2 += double(x) * x;
  const double nn = s1 / n.size(), nw = s2 / s1;
  EXPECT_NEAR(nn, 60, 0.6);
  EXPECT_NEAR(nw / nn, 1.4, 0.02);
  // its density integrates to one
  double area = 0;
  for (double x = 0.5; x < 2000; x += 0.5) area += chain_length_pdf("log-normal", 60, 1.4, x) * 0.5;
  EXPECT_NEAR(area, 1.0, 1e-3);
  // a histogram: lengths drawn in proportion to their weights
  const auto h = draw_chain_lengths(std::vector<std::pair<int, double>>{{10, 1}, {20, 3}}, 20000, 3);
  const double f20 = double(std::count(h.begin(), h.end(), 20)) / h.size();
  EXPECT_NEAR(f20, 0.75, 0.01);
  EXPECT_THROW(draw_chain_lengths(std::vector<std::pair<int, double>>{{1, 1}}, 5, 1), std::invalid_argument);
}

TEST(Polymer, ParallelTrialsGiveTheSameCell) {
  GrowOptions o;
  o.chains = 6;
  o.density = 0.5;
  o.seed = 9;
  const ChainSpec c = spec({"*CC(*)c1ccccc1"}, Sequence::Homopolymer, 12);
  o.threads = 1;
  const auto t0 = std::chrono::steady_clock::now();
  const System a = grow_chains(c, o);
  const auto t1 = std::chrono::steady_clock::now();
  o.threads = 4;
  const System b = grow_chains(c, o);
  const auto t2 = std::chrono::steady_clock::now();
  ASSERT_EQ(a.atoms.size(), b.atoms.size());
  double d = 0;
  for (size_t i = 0; i < a.atoms.size(); ++i) d = std::max(d, norm(a.atoms[i].pos - b.atoms[i].pos));
  EXPECT_EQ(d, 0.0);
  std::printf("serial %.2f s · 4 threads %.2f s\n", std::chrono::duration<double>(t1 - t0).count(), std::chrono::duration<double>(t2 - t1).count());
}

// Structure-based coarse-graining (cg_map.hpp): polyethylene, three backbone carbons per bead — masses conserved, each
// bond type's r0 the mean of its mapped lengths and k = k_B T / (2 var), 1-2 and 1-3 bead pairs excluded, a repulsive
// size where the non-bonded g(r) reaches 1/e
TEST(Polymer, StructureBasedCoarseGraining) {
  ChainSpec c;
  c.units = {{"ethylene", "*CC*"}};
  c.dp = 30;
  GrowOptions o;
  o.chains = 8;
  o.density = 0.3;
  o.auto_scale = true;
  o.seed = 5;
  System aa = grow_chains(c, o);
  RelaxOptions ro;
  ro.target_density = 0.85;
  ro.ftol = 2.0;
  relax(aa, ro);
  make_molecules_whole(aa);
  CgMapOptions m;
  m.scheme = "backbone_n";
  m.per_bead = 3;
  m.temperature = 450;
  const CgMapResult r = cg_map(aa, m);
  double maa = 0, mcg = 0;
  for (const auto& a : aa.atoms) maa += element(a.element).mass;
  for (double x : r.ff->mass) mcg += x;
  EXPECT_NEAR(mcg, maa, 1e-6 * maa);
  EXPECT_EQ(r.beads.atoms.size(), 8u * 20u);   // 60 backbone carbons per chain, 3 per bead
  ASSERT_EQ(r.bonds.size(), 1u);
  // r0 and k from the mapped lengths, independently
  std::vector<double> len;
  for (const auto& b : r.beads.bonds) len.push_back(norm(r.cells.back().minimum_image(r.frames.back()[b.j] - r.frames.back()[b.i])));
  double mean = 0, var = 0;
  for (double x : len) mean += x / double(len.size());
  for (double x : len) var += (x - mean) * (x - mean) / double(len.size() - 1);
  EXPECT_NEAR(r.bonds[0].r0, mean, 1e-9);
  EXPECT_NEAR(r.bonds[0].k, 0.0019872067 * 450 / (2 * var), 1e-6 * r.bonds[0].k);
  EXPECT_GT(r.bonds[0].r0, 3.0);
  EXPECT_LT(r.bonds[0].r0, 4.5);
  EXPECT_GT(r.sigma, 3.0);
  EXPECT_LT(r.sigma, 5.5);
  EXPECT_NEAR(r.cut, r.sigma * std::pow(2.0, 1.0 / 6), 1e-9);
  // 1-2 and 1-3 bead pairs excluded, 1-4 not
  const auto nb = r.beads.neighbours();
  const uint32_t a = r.beads.bonds.front().i;
  for (uint32_t w : nb[a]) {
    EXPECT_TRUE(std::binary_search(r.ff->excluded[a].begin(), r.ff->excluded[a].end(), w));
    for (uint32_t v : nb[w])
      if (v != a) EXPECT_TRUE(std::binary_search(r.ff->excluded[a].begin(), r.ff->excluded[a].end(), v));
  }
  EXPECT_THROW(cg_map(aa, CgMapOptions{"unknown", 3, 300}), std::invalid_argument);
}

TEST(Polymer, RepeatUnitFromPickedAtoms) {
  auto molecule = [](const std::string& smiles) { BuildOptions o; o.forcefield = "uff"; return build_molecule(smiles, o).system; };
  // ethylbenzene: C1 the CH3, C2 the CH2 (then the ring); head on the methyl, tail on the benzylic carbon → styrene's unit
  const System eb = molecule("CCc1ccccc1");
  const std::string smi = repeat_unit_smiles(eb, 0, 1);
  ASSERT_EQ(smi.front(), '*');
  EXPECT_EQ(std::count(smi.begin(), smi.end(), '*'), 2);
  const UnitInfo u = repeat_unit_info(smi);
  EXPECT_EQ(u.formula, "C8H8");
  EXPECT_EQ(u.head_element, "C");
  // hydrogens picked directly: the same unit
  uint32_t h0 = 0, h1 = 0;
  for (const auto& b : eb.bonds) {
    if (b.i == 0 && eb.atoms[b.j].element == 1 && !h0) h0 = b.j;
    if (b.j == 0 && eb.atoms[b.i].element == 1 && !h0) h0 = b.i;
    if (b.i == 1 && eb.atoms[b.j].element == 1 && !h1) h1 = b.j;
    if (b.j == 1 && eb.atoms[b.i].element == 1 && !h1) h1 = b.i;
  }
  EXPECT_EQ(repeat_unit_info(repeat_unit_smiles(eb, h0, h1)).formula, "C8H8");
  // a quaternary carbon has no hydrogen to give
  const System neo = molecule("CC(C)(C)C");
  EXPECT_THROW(repeat_unit_smiles(neo, 0, 1), std::invalid_argument);
  EXPECT_THROW(repeat_unit_smiles(eb, 0, 0), std::invalid_argument);
}

TEST(Polymer, DyadPatternOnAGrownChain) {
  // an atactic chain with its dyads given: the tacticity analysis reads the same pattern back
  ChainSpec c = spec({"*CC(*)c1ccccc1"}, Sequence::Homopolymer, 13);
  c.dyads = "mmr";
  GrowOptions o;
  o.chains = 1;
  o.density = 0.05;
  o.seed = 3;
  o.auto_scale = true;
  const System s = grow_chains(c, o);
  const TacticityReport t = tacticity(s);
  ASSERT_EQ(t.chains.size(), 1u);
  // every dyad between the stereocentres found (a chain end's carbon is not one) follows the pattern
  const std::string& d = t.chains[0].dyads;
  ASSERT_GE(d.size(), 9u);
  for (size_t k = 0; k < d.size(); ++k) EXPECT_EQ(d[k], "mmr"[k % 3]) << d;
}
