#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <map>

#include "caps/pack.hpp"
#include "caps/properties.hpp"
#include "caps/react.hpp"
#include "caps/io.hpp"
#include "caps/json.hpp"
#include "caps/bond_react.hpp"
#include "caps/uff.hpp"
#include <filesystem>
#include <fstream>
#include <sstream>
#include "caps/polymer.hpp"
#include "caps/molecule.hpp"
#include "caps/grow.hpp"
#include "caps/relax.hpp"

using namespace caps;

namespace {

System molecule(const std::vector<std::pair<int, Vec3>>& atoms, const std::vector<Bond>& bonds) {
  System m;
  for (const auto& [el, p] : atoms) {
    Atom a;
    a.element = el;
    a.pos = p;
    a.id = int64_t(m.atoms.size() + 1);
    m.atoms.push_back(a);
  }
  m.bonds = bonds;
  return m;
}

System ethylene_oxide() {
  return molecule({{6, {-0.735, 0, 0}}, {6, {0.735, 0, 0}}, {8, {0, 1.227, 0}},
                   {1, {-1.25, -0.5, 0.9}}, {1, {-1.25, -0.5, -0.9}}, {1, {1.25, -0.5, 0.9}}, {1, {1.25, -0.5, -0.9}}},
                  {{0, 1}, {0, 2}, {1, 2}, {0, 3}, {0, 4}, {1, 5}, {1, 6}});
}

System methylamine() {
  return molecule({{6, {0, 0, 0}}, {7, {1.47, 0, 0}}, {1, {-0.36, 1.03, 0}}, {1, {-0.36, -0.51, 0.89}}, {1, {-0.36, -0.51, -0.89}},
                   {1, {1.81, -0.47, 0.82}}, {1, {1.81, -0.47, -0.82}}},
                  {{0, 1}, {0, 2}, {0, 3}, {0, 4}, {1, 5}, {1, 6}});
}

System packed_epoxy_amine(int eo, int ma, double box) {
  Region b;
  b.kind = Region::InsideBox;
  b.a = {0, 0, 0};
  b.b = {box, box, box};
  PackOptions o;
  o.cell.a = {box, 0, 0};
  o.cell.b = {0, box, 0};
  o.cell.c = {0, 0, box};
  o.seed = 2;
  return pack({{"EO", ethylene_oxide(), eo, {b}}, {"MA", methylamine(), ma, {b}}}, o);
}

std::map<int, int> element_count(const System& s) {
  std::map<int, int> c;
  for (const auto& a : s.atoms) c[a.element]++;
  return c;
}

}  // namespace

TEST(React, TemplatesParseAndValidate) {
  for (const auto& n : builtin_template_names()) {
    const auto t = parse_templates(builtin_template(n));
    ASSERT_EQ(t.size(), 1u) << n;
    EXPECT_EQ(t[0].name, n);
  }
  const auto e = parse_templates(builtin_template("epoxy_amine_primary"))[0];
  EXPECT_EQ(e.atoms.size(), 5u);
  EXPECT_EQ(e.init_a, 4);
  EXPECT_EQ(e.form.size(), 1u);
  EXPECT_EQ(e.move[0], (std::pair<int, int>{5, 2}));
  EXPECT_THROW(parse_templates("atom 1 C"), ReactError);                                             // no reaction line
  EXPECT_THROW(parse_templates("reaction r\natom 1 C\natom 2 C\nform 1 2"), ReactError);             // no initiators
  EXPECT_THROW(parse_templates("reaction r\natom 1 C\natom 2 C bonded 3\ninitiators 1 2"), ReactError);   // bonded to undefined
  EXPECT_THROW(parse_templates("reaction r\natom 1 C\natom 2 C\natom 3 H\ninitiators 1 2"), ReactError);  // 3 not connected
  EXPECT_THROW(parse_templates("reaction r\natom 1 C\natom 2 C\ninitiators 1 2\nbreak 1 2"), ReactError); // break a non-bond
  EXPECT_THROW(parse_templates("reaction r\natom 1 Xx\natom 2 C\ninitiators 1 2"), ReactError);
}

TEST(React, EpoxyAmineTopology) {
  System s = packed_epoxy_amine(20, 10, 16);
  const auto before = element_count(s);
  const size_t nb = s.bonds.size();
  ReactOptions o;
  o.templates = parse_templates(builtin_template("epoxy_amine_primary") + builtin_template("epoxy_amine_secondary"));
  for (auto& t : o.templates) t.capture = 8.0;
  o.relax = false;   // no O / N parameters yet: topology only
  o.max_per_cycle = 2;
  o.max_cycles = 40;
  ReactReport r;
  react(s, o, &r);
  EXPECT_EQ(r.initial_sites, 20);   // one epoxide ring per ethylene oxide (both carbons are CH2)
  EXPECT_GT(r.reactions, 5);
  EXPECT_EQ(element_count(s), before);        // no atoms lost
  EXPECT_EQ(s.bonds.size(), nb);              // +N–C −C–O, N–H → O–H
  // every nitrogen has at most 2 hydrogens; oxygens that opened carry one H; no N–H–O left over
  const auto nbr = s.neighbours();
  int opened = 0;
  for (size_t i = 0; i < s.atoms.size(); ++i) {
    if (s.atoms[i].element == 8) {
      int h = 0, c = 0;
      for (uint32_t j : nbr[i]) (s.atoms[j].element == 1 ? h : c)++;
      if (h == 1) { ++opened; EXPECT_EQ(c, 1); }
      else EXPECT_EQ(c, 2);   // still an epoxide
    }
    if (s.atoms[i].element == 1) EXPECT_EQ(nbr[i].size(), 1u);
  }
  EXPECT_EQ(opened, r.reactions);
  // conversion bookkeeping and clusters
  EXPECT_NEAR(r.cycles.back().conversion, r.reactions / 20.0, 1e-12);
  EXPECT_LT(r.cycles.back().clusters.clusters, 30);
}

TEST(React, CrosslinkRelaxesAndRemovesHydrogen) {
  GrowOptions g;
  g.chains = 4;
  g.dp = 5;
  g.density = 0.4;
  g.seed = 3;
  System s = grow(g);
  RelaxOptions rl;
  rl.target_density = 0.95;
  rl.ftol = 1.0;
  relax(s, rl);
  const size_t atoms = s.atoms.size(), bonds = s.bonds.size();
  ReactOptions o;
  o.templates = parse_templates(builtin_template("cc_crosslink"));
  o.max_per_cycle = 2;
  o.max_cycles = 3;
  ReactReport r;
  react(s, o, &r);
  ASSERT_GT(r.reactions, 0);
  EXPECT_EQ(s.atoms.size(), atoms - 2 * r.reactions);   // H2 leaves per crosslink
  EXPECT_EQ(s.bonds.size(), bonds - r.reactions);        // +C–C, −2 C–H
  EXPECT_TRUE(std::isfinite(r.cycles.back().energy));
  EXPECT_LE(r.cycles.back().clusters.clusters, 4);
  // still typed: every carbon has four neighbours or is aromatic
  const ForceField ff = assign_gaff(s);
  EXPECT_EQ(ff.atom_type.size(), s.atoms.size());
  EXPECT_GT(r.cycles.front().max_force, 0);   // the force after each cycle's relaxation
}

// FailedJob: a cycle that fails leaves the structure as the last completed cycle did, and the report says which one
TEST(React, FailedCycleKeepsTheCompletedOnes) {
  GrowOptions g;
  g.chains = 4;
  g.dp = 5;
  g.density = 0.4;
  g.seed = 3;
  System s = grow(g);
  RelaxOptions rl;
  rl.target_density = 0.95;
  rl.ftol = 1.0;
  relax(s, rl);
  ReactOptions o;
  o.templates = parse_templates(builtin_template("cc_crosslink"));
  o.max_per_cycle = 1;
  o.max_cycles = 4;
  System after1;
  o.frame = [&](const System& x, int cycle) {
    if (cycle == 1) after1 = x;
    if (cycle == 2) throw std::runtime_error("disk full");
  };
  ReactReport r;
  react(s, o, &r);
  EXPECT_EQ(r.failed_cycle, 2);
  EXPECT_NE(r.failure.find("disk full"), std::string::npos);
  ASSERT_EQ(r.cycles.size(), 1u);
  EXPECT_EQ(r.reactions, r.cycles[0].total);
  EXPECT_EQ(s.atoms.size(), after1.atoms.size());
  EXPECT_EQ(s.bonds.size(), after1.bonds.size());
  // off: the failure is thrown, as before; and a failure in the first cycle always is
  o.keep_on_failure = false;
  System t = s;
  EXPECT_THROW(react(t, o), std::runtime_error);
  o.keep_on_failure = true;
  o.frame = [](const System&, int) { throw std::runtime_error("at once"); };
  EXPECT_THROW(react(t, o), std::runtime_error);
}

TEST(React, ClustersAndFloryStockmayer) {
  EXPECT_NEAR(flory_stockmayer(1, 2, 4), 1 / std::sqrt(3.0), 1e-12);   // DGEBA + IPDA
  EXPECT_NEAR(flory_stockmayer(1, 3, 3), 0.5, 1e-12);
  // three carbons bonded in a chain plus one isolated carbon
  System s = molecule({{6, {0, 0, 0}}, {6, {1.5, 0, 0}}, {6, {3, 0, 0}}, {6, {10, 0, 0}}}, {{0, 1}, {1, 2}});
  const auto c = cluster_stats(s);
  EXPECT_EQ(c.clusters, 2);
  EXPECT_NEAR(c.largest_fraction, 0.75, 1e-12);
  const double m = 12.011;
  EXPECT_NEAR(c.mw, (9 * m * m + m * m) / (4 * m), 1e-9);
  EXPECT_NEAR(c.reduced_mw, m, 1e-9);
}

TEST(React, SameSeedSameNetwork) {
  System a = packed_epoxy_amine(12, 6, 14), b = a;
  ReactOptions o;
  o.templates = parse_templates(builtin_template("epoxy_amine_primary"));
  o.templates[0].capture = 7;
  o.templates[0].probability = 0.5;
  o.relax = false;
  o.max_cycles = 5;
  o.seed = 4;
  react(a, o);
  react(b, o);
  ASSERT_EQ(a.bonds.size(), b.bonds.size());
  for (size_t k = 0; k < a.bonds.size(); ++k) EXPECT_EQ(a.bonds[k].i, b.bonds[k].i);
}

TEST(React, SulfurCureOfNaturalRubber) {
  ChainSpec spec;
  spec.units = {{"cis-1,4-isoprene", "[*]C/C=C(C)\\C[*]"}};
  spec.dp = 10;
  GrowOptions g;
  g.chains = 4;
  g.density = 0.6;
  System s = grow_chains(spec, g);
  BuildResult donor = build_molecule("SS");   // H–S–S–H
  PackReport pr;
  s = insert_molecules(s, donor.system, 8, PackOptions{}, &pr);
  ASSERT_EQ(s.atoms.size(), size_t(4 * (10 * 13 + 2) + 8 * 4));
  int64_t top = 0;
  for (const auto& a : s.atoms) top = std::max(top, a.mol);
  EXPECT_EQ(top, 4 + 8);
  ReactOptions r;
  r.templates = parse_templates(builtin_template("sulfur_allylic"));
  r.relax = false;
  r.max_cycles = 20;
  ReactReport rep;
  react(s, r, &rep);
  EXPECT_GT(rep.reactions, 4);
  // every S keeps two neighbours; some donors bridge two chains (C–S–S–C)
  std::vector<std::vector<uint32_t>> nb(s.atoms.size());
  for (const auto& b : s.bonds) nb[b.i].push_back(b.j), nb[b.j].push_back(b.i);
  int bridges = 0;
  for (uint32_t i = 0; i < s.atoms.size(); ++i) {
    if (s.atoms[i].element != 16) continue;
    EXPECT_EQ(nb[i].size(), 2u);
    for (uint32_t j : nb[i])
      if (s.atoms[j].element == 16 && i < j) {
        bool ci = false, cj = false;
        for (uint32_t k : nb[i]) ci = ci || s.atoms[k].element == 6;
        for (uint32_t k : nb[j]) cj = cj || s.atoms[k].element == 6;
        bridges += ci && cj;
      }
  }
  EXPECT_GT(bridges, 0);
  // the crosslink analysis counts the same bridges
  Trajectory t;
  t.topology = s;
  std::vector<Vec3> p;
  for (const auto& a : s.atoms) p.push_back(a.pos);
  t.positions.push_back(p);
  t.cells.push_back(s.cell);
  t.timesteps.push_back(0);
  const auto props = analyze(t, {"crosslinks"}, AnalyzeOptions{});
  ASSERT_EQ(props.size(), 1u);
  EXPECT_EQ(int(props[0].extra.at("sulfur bridges")), bridges);
  EXPECT_GT(props[0].value, 0.0);
}

// REACTER-style: one continuous NVT run with reaction checks every 0.05 ps; reacted sites settle locally, the run goes on
// through checks without reactions (no stall stop), and the cell keeps its velocities
TEST(React, DuringMdChecksAtIntervals) {
  GrowOptions g;
  g.chains = 4;
  g.dp = 5;
  g.density = 0.4;
  g.seed = 3;
  System s = grow(g);
  RelaxOptions rl;
  rl.target_density = 0.95;
  rl.ftol = 1.0;
  relax(s, rl);
  const size_t atoms = s.atoms.size();
  ReactOptions o;
  o.templates = parse_templates(builtin_template("cc_crosslink"));
  o.max_per_cycle = 2;
  o.max_cycles = 5;
  o.during_md = true;
  o.md_ps = 0.05;
  o.temperature = 400;
  ReactReport r;
  react(s, o, &r);
  EXPECT_EQ(r.cycles.size(), 5u);   // every check runs, reacting or not
  ASSERT_GT(r.reactions, 0);
  EXPECT_EQ(s.atoms.size(), atoms - 2 * r.reactions);
  EXPECT_EQ(s.velocities.size(), s.atoms.size());
  for (const auto& c : r.cycles) EXPECT_TRUE(std::isfinite(c.energy));
  EXPECT_TRUE(std::any_of(r.notes.begin(), r.notes.end(), [](const std::string& n) { return n.find("REACTER-style") != std::string::npos; }));
  o.md_ps = 0;
  EXPECT_THROW(react(s, o), ReactError);
}

// Silane coupling to rubber: the tetrasulfide of diethyl tetrasulfide (TESPT's polysulfide without its silyl groups)
// opens S–S bonds onto allylic carbons of natural rubber; every sulfur keeps two neighbours, each moved H is an S–H
TEST(React, PolysulfideCouplesToNaturalRubber) {
  ChainSpec spec;
  spec.units = {{"cis-1,4-isoprene", "[*]C/C=C(C)\\C[*]"}};
  spec.dp = 10;
  GrowOptions g;
  g.chains = 4;
  g.density = 0.6;
  System s = grow_chains(spec, g);
  BuildResult donor = build_molecule("CCSSSSCC");
  s = insert_molecules(s, donor.system, 6, PackOptions{});
  auto count_ss = [](const System& t) {
    int n = 0;
    for (const auto& b : t.bonds) n += t.atoms[b.i].element == 16 && t.atoms[b.j].element == 16;
    return n;
  };
  const int ss0 = count_ss(s);
  ReactOptions r;
  r.templates = parse_templates(builtin_template("polysulfide_allylic"));
  r.relax = false;
  r.max_cycles = 20;
  ReactReport rep;
  react(s, r, &rep);
  ASSERT_GT(rep.reactions, 2);
  EXPECT_EQ(count_ss(s), ss0 - rep.reactions);   // one S–S opened per reaction
  std::vector<std::vector<uint32_t>> nb(s.atoms.size());
  for (const auto& b : s.bonds) nb[b.i].push_back(b.j), nb[b.j].push_back(b.i);
  int sh = 0, sc_rubber = 0;
  for (uint32_t i = 0; i < s.atoms.size(); ++i) {
    if (s.atoms[i].element != 16) continue;
    EXPECT_EQ(nb[i].size(), 2u);
    for (uint32_t j : nb[i]) sh += s.atoms[j].element == 1;
    // a rubber carbon: allylic, bonded to a carbon with three neighbours (the C=C), unlike the donors' ethyl carbons
    for (uint32_t j : nb[i]) {
      if (s.atoms[j].element != 6) continue;
      bool allylic = false;
      for (uint32_t k : nb[j]) allylic = allylic || (s.atoms[k].element == 6 && nb[k].size() == 3);
      sc_rubber += allylic;
    }
  }
  EXPECT_EQ(sh, rep.reactions);
  EXPECT_EQ(sc_rubber, rep.reactions);
}

// Reaction SMARTS (design/boards/ReactionTemplate "Saved as atom-mapped reaction SMARTS"): the query atoms with their
// constraints, products by map number, deleted atoms absent, substitution bonds single. Checked with an external
// reaction-SMARTS toolkit (2026): run on ethane + ethane / propylene oxide + methylamine, they give butane and
// 1-(methylamino)propan-2-ol.
TEST(React, TemplatesAsReactionSmarts) {
  const auto cc = parse_templates(builtin_template("cc_crosslink"));
  EXPECT_EQ(reaction_smarts(cc[0]), "[#6;X4;H2;A:1]~[#1:3].[#6;X4;H2;A:2]~[#1:4]>>[#6:1]-[#6:2].[#1:3]-[#1:4]");
  const auto ep = parse_templates(builtin_template("epoxy_amine_primary"));
  EXPECT_EQ(reaction_smarts(ep[0]), "[#6;H2;r3:1]~1~[#8;r3:2]~[#6;r3:3]~1.[#7;H2:4]~[#1:5]>>[#6:1](~[#6:3]~[#8:2]-[#1:5])-[#7:4]");
  // the editor's view carries it
  EXPECT_NE(template_view(cc[0]).find("\"smarts\""), std::string::npos);
}

// "charges keep": the atoms keep their charges and a deleted hydrogen's joins its carbon, so the net charge is conserved
// exactly; without it the charges are dropped for the force field to recompute
TEST(React, ChargesKeptThroughTheReaction) {
  GrowOptions g;
  g.chains = 4;
  g.dp = 5;
  g.density = 0.4;
  g.seed = 3;
  System s = grow(g);
  RelaxOptions rl;
  rl.target_density = 0.95;
  rl.ftol = 1.0;
  relax(s, rl);
  const ForceField ff = assign_gaff(s);
  for (size_t i = 0; i < s.atoms.size(); ++i) s.atoms[i].charge = ff.charge[i] + (i % 7 == 0 ? 0.01 : 0.0);   // a net charge to follow
  s.has_charges = true;
  double q0 = 0;
  for (const auto& a : s.atoms) q0 += a.charge;
  auto t = parse_templates(builtin_template("cc_crosslink") + "charges keep\n");
  ASSERT_TRUE(t[0].keep_charges);
  System k = s;
  const auto m = find_matches(k, t[0]);
  ASSERT_FALSE(m.empty());
  ASSERT_EQ(apply_matches(k, t, {m.front()}), 1);
  double q1 = 0;
  for (const auto& a : k.atoms) q1 += a.charge;
  EXPECT_TRUE(k.has_charges);
  EXPECT_EQ(k.atoms.size(), s.atoms.size() - 2);
  EXPECT_NEAR(q1, q0, 1e-12);
  // the default: recomputed for the new chemistry
  auto d = parse_templates(builtin_template("cc_crosslink"));
  System r = s;
  apply_matches(r, d, {find_matches(r, d[0]).front()});
  EXPECT_FALSE(r.has_charges);
  EXPECT_THROW(parse_templates(builtin_template("cc_crosslink") + "charges maybe\n"), ReactError);
}

namespace {

System pe_cell(int chains, int dp, uint64_t seed) {
  GrowOptions g;
  g.chains = chains;
  g.dp = dp;
  g.density = 0.4;
  g.seed = seed;
  System s = grow(g);
  RelaxOptions rl;
  rl.target_density = 0.9;
  rl.ftol = 2.0;
  relax(s, rl);
  return s;
}

System packed(const std::vector<std::pair<std::string, int>>& mols, double box, uint64_t seed = 4) {
  Region b;
  b.kind = Region::InsideBox;
  b.a = {0, 0, 0};
  b.b = {box, box, box};
  PackOptions o;
  o.cell.a = {box, 0, 0};
  o.cell.b = {0, box, 0};
  o.cell.c = {0, 0, box};
  o.seed = seed;
  std::vector<PackItem> c;
  for (const auto& [smiles, n] : mols) c.push_back({smiles, build_molecule(smiles).system, n, {b}});
  return pack(c, o);
}

// the molecules of s as element formulas ("C4H8O2" …), counted
std::map<std::string, int> formulas(const System& s) {
  System t = s;
  t.has_mol = false;
  int nm = 0;
  const auto mol = t.molecules(&nm);
  std::vector<std::map<int, int>> el(static_cast<size_t>(nm));
  for (size_t i = 0; i < s.atoms.size(); ++i) el[size_t(mol[i])][s.atoms[i].element]++;
  std::map<std::string, int> out;
  for (const auto& e : el) {
    std::string f;
    for (int z : {6, 1, 8, 7, 16})
      if (e.count(z)) f += std::string(z == 6 ? "C" : z == 1 ? "H" : z == 8 ? "O" : z == 7 ? "N" : "S") + (e.at(z) > 1 ? std::to_string(e.at(z)) : "");
    out[f]++;
  }
  return out;
}

}  // namespace

// Ten chains (the default polymer of Grow) crosslinked only between chains, to a target of 1.2 links per chain: the run stops at the target
// (6 links), counts no loop within a chain, and every C–C link joins two chains of the start
TEST(React, BetweenChainsToACrosslinkTarget) {
  System s = pe_cell(10, 12, 5);
  std::vector<int64_t> start(s.atoms.size());
  {
    System t = s;
    t.has_mol = false;
    const auto m = t.molecules();
    for (size_t i = 0; i < s.atoms.size(); ++i) start[i] = m[i];
  }
  // carbons keep their order (only hydrogens leave): the k-th carbon of the product is the k-th carbon of the start
  std::vector<int64_t> carbon_chain;
  for (size_t i = 0; i < s.atoms.size(); ++i) if (s.atoms[i].element == 6) carbon_chain.push_back(start[i]);
  int within_before = 0;
  for (const auto& b : s.bonds)
    if (s.atoms[b.i].element == 6 && s.atoms[b.j].element == 6) within_before += start[b.i] == start[b.j];
  ReactOptions o;
  o.templates = parse_templates(builtin_template("cc_crosslink"));
  o.between_chains = true;
  o.target = ReactTarget::PerChain;
  o.target_value = 1.2;
  o.max_per_cycle = 2;
  o.max_cycles = 40;
  o.auto_capture = true;
  o.relax = false;
  ReactReport r;
  react(s, o, &r);
  EXPECT_EQ(r.chains, 10);
  EXPECT_EQ(r.target_crosslinks, 6);
  EXPECT_EQ(r.crosslinks, 6);
  EXPECT_EQ(r.intrachain, 0);
  EXPECT_NEAR(r.per_chain, 1.2, 1e-12);
  EXPECT_GT(r.density, 0);
  EXPECT_NEAR(r.mc, r.chain_mass / 12.0, 1e-9);
  std::vector<int> carbon_index(s.atoms.size(), -1);
  int k = 0;
  for (size_t i = 0; i < s.atoms.size(); ++i) if (s.atoms[i].element == 6) carbon_index[i] = k++;
  int cc_between = 0, cc_within = 0;
  for (const auto& b : s.bonds)
    if (carbon_index[b.i] >= 0 && carbon_index[b.j] >= 0)
      (carbon_chain[size_t(carbon_index[b.i])] == carbon_chain[size_t(carbon_index[b.j])] ? cc_within : cc_between)++;
  EXPECT_EQ(cc_between, 6);
  EXPECT_EQ(cc_within, within_before);   // each chain's own C-C bonds, untouched
}

// Byproducts kept: each C–C crosslink leaves an H2 molecule in the cell; removed, the atoms go
TEST(React, ByproductsKeptOrRemoved) {
  for (bool keep : {true, false}) {
    System s = pe_cell(4, 8, 7);
    const size_t atoms = s.atoms.size();
    ReactOptions o;
    o.templates = parse_templates(builtin_template("cc_crosslink"));
    o.keep_byproducts = keep;
    o.max_cycles = 3;
    o.relax = false;
    ReactReport r;
    react(s, o, &r);
    ASSERT_GT(r.reactions, 0);
    EXPECT_EQ(r.byproducts, r.reactions);
    const auto f = formulas(s);
    if (keep) {
      EXPECT_EQ(s.atoms.size(), atoms);
      EXPECT_EQ(f.count("H2") ? f.at("H2") : 0, r.reactions);
    } else {
      EXPECT_EQ(s.atoms.size(), atoms - 2 * r.reactions);
      EXPECT_EQ(f.count("H2"), 0u);
    }
  }
}

// The force field given for the run types the structure after every cycle that reacted (the user's assignment in the
// Studio), and relaxes it
TEST(React, RunsWithTheGivenForceField) {
  System s = pe_cell(4, 8, 9);
  ReactOptions o;
  o.templates = parse_templates(builtin_template("cc_crosslink"));
  o.max_cycles = 3;
  int calls = 0;
  o.retype = [&](const System& x) {
    ++calls;
    return std::make_shared<const ForceField>(assign_gaff(x));
  };
  o.field_name = "GAFF (test)";
  ReactReport r;
  react(s, o, &r);
  ASSERT_GT(r.reactions, 0);
  int reacted = 0;
  for (const auto& c : r.cycles) reacted += c.reactions > 0;
  EXPECT_EQ(calls, reacted);
  EXPECT_EQ(r.field, "GAFF (test)");
  EXPECT_TRUE(std::isfinite(r.cycles.back().energy));
}

// Auto capture: a template whose capture is too short for any pair widens until pairs are found
TEST(React, AutoCaptureWidens) {
  System s = pe_cell(4, 8, 11);
  ReactOptions o;
  o.templates = parse_templates(builtin_template("cc_crosslink"));
  o.templates[0].capture = 0.5;
  o.max_cycles = 2;
  o.relax = false;
  ReactReport r;
  react(s, o, &r);
  EXPECT_EQ(r.reactions, 0);
  o.auto_capture = true;
  o.capture_max = 4.0;
  s = pe_cell(4, 8, 11);
  react(s, o, &r);
  EXPECT_GT(r.reactions, 0);
  EXPECT_GT(r.cycles.front().capture, 0.5);
  EXPECT_LE(r.cycles.back().capture, 4.0 + 1e-9);
}

// The ENR / PBS / MAH chemistry on model compounds: an ENR-type trisubstituted epoxide opened by acetic acid gives the
// β-hydroxy ester (no water); acetic acid + ethanol condense to ethyl acetate and water (kept); maleic anhydride opened by
// methanol gives the half-ester acid (no water)
TEST(React, EpoxideAcidEsterAndAnhydrideChemistry) {
  {
    System s = packed({{"CC1(C)OC1C", 6}, {"CC(=O)O", 6}}, 16);
    ReactOptions o;
    o.templates = parse_templates(builtin_template("enr_acid_ester"));
    o.relax = false;
    o.auto_capture = true;
    o.max_cycles = 10;
    ReactReport r;
    react(s, o, &r);
    ASSERT_GT(r.reactions, 0);
    const auto f = formulas(s);
    EXPECT_EQ(f.at("C7H14O3"), r.reactions);   // 2-methyl-2,3-epoxybutane + acetic acid, one molecule, nothing lost
    EXPECT_EQ(f.count("H2O"), 0u);
    EXPECT_EQ(r.byproducts, 0);
  }
  {
    System s = packed({{"CC(=O)O", 6}, {"CCO", 6}}, 14);
    ReactOptions o;
    o.templates = parse_templates(builtin_template("ester_condensation"));
    o.relax = false;
    o.keep_byproducts = true;
    o.auto_capture = true;
    o.max_cycles = 10;
    ReactReport r;
    react(s, o, &r);
    ASSERT_GT(r.reactions, 0);
    const auto f = formulas(s);
    EXPECT_EQ(f.at("C4H8O2"), r.reactions);   // ethyl acetate
    EXPECT_EQ(f.at("H2O"), r.reactions);
  }
  {
    System s = packed({{"O=C1OC(=O)C=C1", 6}, {"CO", 6}}, 14);
    ReactOptions o;
    o.templates = parse_templates(builtin_template("anhydride_alcohol"));
    o.relax = false;
    o.auto_capture = true;
    o.max_cycles = 10;
    ReactReport r;
    react(s, o, &r);
    ASSERT_GT(r.reactions, 0);
    const auto f = formulas(s);
    EXPECT_EQ(f.at("C5H6O4"), r.reactions);   // monomethyl maleate: one ester, one COOH
    EXPECT_EQ(f.count("H2O"), 0u);
  }
}

// fix bond/react: templates cut from a cell and typed with its force field, then read back as a CAPS template that
// crosslinks another cell; the data file keeps only the cell's atoms while holding every type the reaction creates
TEST(React, BondReactExportAndImport) {
  System s = pe_cell(6, 10, 13);
  const auto dir = (std::filesystem::temp_directory_path() / "caps_bond_react_test").string();
  std::filesystem::remove_all(dir);
  BondReactOptions bo;
  bo.max_variants = 2;
  auto gaff = [](const System& x) { return std::make_shared<const ForceField>(assign_gaff(x)); };
  const auto rep = write_bond_react(s, parse_templates(builtin_template("cc_crosslink")), gaff, dir, "rx", bo);
  ASSERT_GE(rep.variants.size(), 1u);
  EXPECT_GT(rep.candidates, 0);
  EXPECT_LE(rep.covered, rep.candidates);
  // the data file: the cell's atoms only
  {
    std::ifstream in(dir + "/rx.data");
    std::string line;
    std::getline(in, line);
    std::getline(in, line);
    std::getline(in, line);
    EXPECT_EQ(line, std::to_string(s.atoms.size()) + " atoms");
  }
  const std::string v = dir + "/rx_" + rep.variants[0].name;
  std::ifstream map(v + "_map.txt");
  std::stringstream ms;
  ms << map.rdbuf();
  EXPECT_NE(ms.str().find("InitiatorIDs"), std::string::npos);
  EXPECT_NE(ms.str().find("2 deleteIDs"), std::string::npos);   // the two hydrogens (H2 removed)
  std::vector<std::string> notes;
  const ReactionTemplate t = read_bond_react(v + "_pre.mol", v + "_post.mol", v + "_map.txt", dir + "/rx.data", "back", 3.0, &notes);
  EXPECT_EQ(t.form.size(), 1u);
  EXPECT_EQ(t.remove.size(), 2u);
  System s2 = pe_cell(6, 10, 17);
  ReactOptions o;
  o.templates = {t};
  o.relax = false;
  o.max_cycles = 3;
  ReactReport r;
  react(s2, o, &r);
  EXPECT_GT(r.reactions, 0);
  map.close();   // Windows will not remove a file still open
  std::filesystem::remove_all(dir);
}

// fix bond/react with the cure options: types split by component (a renamed copy: the energy unchanged; each group's types
// numbered in the given order, a LAMMPS group), a virtual-cure survey, molecule ids kept by molmap (templates with a
// Molecules section) and crosslink-density targets in the input
TEST(React, BondReactCureOptions) {
  System s = pe_cell(6, 10, 13);
  auto gaff = [](const System& x) { return std::make_shared<const ForceField>(assign_gaff(x)); };
  // split types: same energy
  {
    const ForceField ff = assign_gaff(s);
    std::vector<int> g(s.atoms.size());
    for (size_t i = 0; i < s.atoms.size(); ++i) g[i] = s.atoms[i].mol <= 3 ? 1 : 0;
    const ForceField sp = split_types_by_group(ff, g, {"B", "A"});
    EXPECT_EQ(sp.type_names.size(), 2 * ff.type_names.size());
    EXPECT_EQ(sp.type_index[0], int(ff.type_names.size()) + ff.type_index[0]);   // molecule 1 is in A, numbered after B's
    EnergyOptions eo;
    std::vector<double> x, f1, f2;
    for (const auto& a : s.atoms) x.insert(x.end(), a.pos.begin(), a.pos.end());
    Evaluator e1(ff, eo), e2(sp, eo);
    const auto t1 = e1.compute(x, s.cell, f1), t2 = e2.compute(x, s.cell, f2);
    EXPECT_NEAR(t1.vdw, t2.vdw, 1e-9);
    EXPECT_NEAR(t1.bond + t1.angle + t1.dihedral, t2.bond + t2.angle + t2.dihedral, 1e-9);
  }
  const auto dir = (std::filesystem::temp_directory_path() / "caps_bond_react_cure").string();
  std::filesystem::remove_all(dir);
  BondReactOptions bo;
  bo.max_variants = 2;
  bo.between_chains = true;
  bo.survey_after = {0.5};
  bo.survey_relax = false;
  bo.mol_ids = "molmap";
  bo.targets = {0.5, 1.0};
  bo.limiting = 2;
  bo.check_every = 500;
  bo.max_steps = 5000;
  for (const auto& a : s.atoms) bo.type_group.push_back(a.mol <= 3 ? 0 : 1);
  bo.type_group_names = {"first", "second"};
  const auto rep = write_bond_react(s, parse_templates(builtin_template("cc_crosslink")), gaff, dir, "cure", bo);
  ASSERT_GE(rep.variants.size(), 1u);
  ASSERT_EQ(rep.frames.size(), 1u);
  EXPECT_GT(rep.frames[0].reactions, 0);
  ASSERT_FALSE(rep.steps.empty());
  EXPECT_EQ(rep.steps[0].covered, std::min(rep.steps[0].candidates, rep.steps[0].covered));
  EXPECT_FALSE(rep.link_reactions.empty());   // a C–C link between two chains
  auto text = [](const std::string& path) {
    std::ifstream in(path);
    std::stringstream ss;
    ss << in.rdbuf();
    return ss.str();
  };
  const std::string in = text(dir + "/cure.in");
  EXPECT_NE(in.find("reset_mol_ids molmap"), std::string::npos);
  EXPECT_NE(in.find("group           first          type"), std::string::npos);
  EXPECT_NE(in.find("group           second         type"), std::string::npos);
  EXPECT_NE(in.find("variable        tgt index 50 100"), std::string::npos);
  EXPECT_NE(in.find("write_data      cure_XL${tgt}.data"), std::string::npos);
  EXPECT_NE(in.find("crosslink_progress.dat"), std::string::npos);
  EXPECT_NE(text(dir + "/cure_" + rep.variants[0].name + "_pre.mol").find("Molecules"), std::string::npos);
  std::filesystem::remove_all(dir);
}

// Polybutadiene vulcanised with trisulfide donors (H–S3–H) between chains to a degree of crosslinking DC = 2 links / monomers
// of 10 % (the definition of Vasilev et al. 2021 and Alamfard et al. 2023): C–S3–C bridges form, one link per donor that
// reached two chains, and the run stops at the target
TEST(React, SulfurBridgesToADegreeOfCrosslinking) {
  ChainSpec spec;
  spec.units = {{"cis-1,4-butadiene", "[*]C/C=C\\C[*]"}};
  spec.dp = 20;
  GrowOptions g;
  g.chains = 6;
  g.density = 0.6;
  g.seed = 21;
  System s = grow_chains(spec, g);
  BuildResult donor = build_molecule("SSS");
  PackReport pr;
  s = insert_molecules(s, donor.system, 12, PackOptions{}, &pr);
  ReactOptions r;
  r.templates = parse_templates(builtin_template("sulfur_allylic"));
  r.relax = false;
  r.between_chains = true;
  r.auto_capture = true;
  r.capture_max = 9;
  r.target = ReactTarget::DegreePercent;
  r.target_value = 10;
  r.max_cycles = 60;
  r.max_per_cycle = 2;
  ReactReport rep;
  react(s, r, &rep);
  EXPECT_EQ(rep.monomers, 120);
  EXPECT_EQ(rep.target_crosslinks, 6);   // 10 % × 120 / 2
  EXPECT_LE(rep.crosslinks, 6);
  EXPECT_GT(rep.crosslinks, 0);
  EXPECT_NEAR(rep.degree, 200.0 * rep.crosslinks / 120, 1e-9);
  // bridges: a sulfur chain with a carbon at both ends
  std::vector<std::vector<uint32_t>> nb(s.atoms.size());
  for (const auto& b : s.bonds) nb[b.i].push_back(b.j), nb[b.j].push_back(b.i);
  int ends_on_carbon = 0;
  for (uint32_t i = 0; i < s.atoms.size(); ++i)
    if (s.atoms[i].element == 16)
      for (uint32_t j : nb[i]) ends_on_carbon += s.atoms[j].element == 6;
  EXPECT_GE(ends_on_carbon, 2 * rep.crosslinks);   // each link: both ends of its S3 on carbon
}

// Crosslinking goes on once every chain is in one network: "different molecules" (a template's min_path 0) are the
// molecules the run started from, not the connected pieces now — so a donor hanging from chain A still meets chain B after
// the network percolates. 16 cis-PB chains of 200 units at 1.0 g/cm³ with 227 H2S donors reach 197 links between chains
// (it stopped at about 17 before); the network written and read back (its chains no longer recorded) has its 16 chains
// recovered, its pendant S–H ends go on reacting, and its sites are counted on the original chains and repeat units
TEST(React, CrosslinkingContinuesAfterPercolation) {
  ChainSpec spec;
  spec.units = {{"cis-1,4-butadiene", "[*]C/C=C\\C[*]"}};
  spec.dp = 200;
  GrowOptions g;
  g.chains = 16;
  g.density = 1.0;
  g.seed = 1;
  g.auto_scale = true;   // as caps.polymer grows: contact limits stepped down only where a chain cannot be placed
  System s = grow_chains(spec, g);
  PackOptions po;
  po.tolerance = 2.0;
  po.seed = 1;
  PackReport pr;
  s = insert_molecules(s, build_molecule("S").system, 227, po, &pr);
  const auto templates = parse_templates(builtin_template("sulfur_allylic"));
  ReactOptions r;
  r.templates = templates;
  r.relax = false;
  r.between_chains = true;
  r.target = ReactTarget::Crosslinks;
  r.target_value = 197;
  r.auto_capture = true;
  r.capture_max = 8;
  r.capture_step = 0.5;
  r.max_cycles = 600;
  r.max_per_cycle = 20;
  r.seed = 1;
  ReactReport rep;
  react(s, r, &rep);
  EXPECT_EQ(rep.chains, 16);
  EXPECT_EQ(rep.monomers, 3200);
  EXPECT_GE(rep.crosslinks, 190);
  EXPECT_EQ(rep.intrachain, 0);
  EXPECT_NEAR(rep.degree, 200.0 * rep.crosslinks / 3200, 1e-9);
  EXPECT_GT(cluster_stats(s).largest_fraction, 0.99);   // one network
  // written and read back: the chains it started from recorded on the atoms ("# chain N")
  const auto path = (std::filesystem::temp_directory_path() / "caps_percolated.data").string();
  write_lammps_data(s, path);
  System back = read_lammps_data(path);
  std::filesystem::remove(path);
  {
    std::set<int64_t> chains;
    for (const auto& a : back.atoms) chains.insert(a.chain);
    EXPECT_FALSE(chains.count(0));
    std::string note;
    const auto oc = original_chains(back, templates, &note);
    EXPECT_NE(note.find("records"), std::string::npos);
    EXPECT_EQ(std::set<int64_t>(oc.begin(), oc.end()).size(), chains.size());
  }
  // without the record (a file another program wrote): the chains recovered from the C–S bonds the template forms
  for (auto& a : back.atoms) a.chain = 0;
  const auto cs = chain_sites(back, templates);
  ASSERT_EQ(cs.size(), 16u);
  int units = 0, pendants = 0;
  double mass = 0;
  for (const auto& c : cs) units += c.units, pendants += c.sites, mass += c.mass;
  EXPECT_EQ(units, 3200);
  EXPECT_NEAR(mass / units, 54.09, 0.2);   // C4H6, less the hydrogens the links took
  ASSERT_GT(pendants, 0);
  ReactOptions more = r;
  more.target = ReactTarget::Conversion;
  more.target_conversion = 1.0;
  more.between_chains = false;
  more.max_cycles = 50;
  ReactReport rep2;
  react(back, more, &rep2);
  EXPECT_EQ(rep2.chains, 16);
  EXPECT_TRUE(std::any_of(rep2.notes.begin(), rep2.notes.end(), [](const std::string& n) { return n.find("16 chains recovered") != std::string::npos; }));
  for (const auto& a : back.atoms) EXPECT_GT(a.chain, 0);   // the run records them in turn
  EXPECT_GE(rep2.reactions, pendants * 9 / 10);   // most pendant S–H ends react
}

// Two chains already joined by one link take a second when links between chains are asked for: the pendant donor on one
// chain bonds to the other although the two are one molecule now
TEST(React, SecondLinkBetweenJoinedChains) {
  ChainSpec spec;
  spec.units = {{"cis-1,4-butadiene", "[*]C/C=C\\C[*]"}};
  spec.dp = 30;
  GrowOptions g;
  g.chains = 2;
  g.density = 0.6;
  g.seed = 4;
  System s = grow_chains(spec, g);
  PackReport pr;
  s = insert_molecules(s, build_molecule("S").system, 6, PackOptions{}, &pr);
  ReactOptions r;
  r.templates = parse_templates(builtin_template("sulfur_allylic"));
  r.relax = false;
  r.between_chains = true;
  r.target = ReactTarget::Crosslinks;
  r.target_value = 2;
  r.auto_capture = true;
  r.capture_max = 12;
  r.max_cycles = 80;
  r.max_per_cycle = 1;
  ReactReport rep;
  react(s, r, &rep);
  EXPECT_EQ(rep.chains, 2);
  EXPECT_EQ(rep.crosslinks, 2);
  EXPECT_EQ(rep.intrachain, 0);
}

// The reaction library (data/reactions): every scheme of two molecules becomes a template, and run on its own model
// compounds it makes the scheme's products — the main product's formula, and the byproducts' when they are kept
TEST(React, ReactionLibraryRunsOnItsModelCompounds) {
  std::ifstream in(std::string(CAPS_SOURCE_DIR) + "/data/reactions/library.json");
  std::stringstream ss;
  ss << in.rdbuf();
  const Json lib = Json::parse(ss.str());
  int converted = 0, refused = 0, ran = 0;
  std::vector<std::string> failures;
  for (const auto& e : lib["reactions"].items()) {
    std::vector<std::string> rs, ps;
    for (const auto& r : e["reactants"].items()) rs.push_back(r.text("smiles", ""));
    for (const auto& p : e["products"].items()) ps.push_back(p.text("smiles", ""));
    const std::string id = e.text("id", "");
    ReactionTemplate t;
    try {
      t = template_from_scheme(rs, ps, id);
      ++converted;
    } catch (const ReactError& x) {
      ++refused;
      if (rs.size() < 3) failures.push_back(id + ": " + x.what());   // only three-molecule schemes may be refused
      continue;
    }
    if (rs.size() > 2) continue;
    // the model compounds packed together, reacted (topology only), the products looked for by formula
    std::vector<std::pair<std::string, int>> mols;
    for (const auto& r : rs) mols.push_back({r, rs.size() == 1 ? 8 : 5});
    System s;
    try { s = packed(mols, 16, 7); } catch (const std::exception& x) { failures.push_back(id + ": packing: " + x.what()); continue; }
    ReactOptions o;
    o.templates = {t};
    o.relax = false;
    o.keep_byproducts = true;
    o.auto_capture = true;
    o.capture_max = 10;
    o.max_cycles = 6;
    o.max_per_cycle = 1;
    // the formulas seen after every cycle: a product may react on later (a hemiacetal's OH is an alcohol too, a
    // methylol bridge links again), so the scheme's product must appear at some point, not survive to the end
    std::map<std::string, int> f;
    o.frame = [&](const System& sys, int) { for (const auto& [h, c] : formulas(sys)) f[h] = std::max(f[h], c); };
    ReactReport rep;
    try { react(s, o, &rep); } catch (const std::exception& x) { failures.push_back(id + ": react: " + x.what()); continue; }
    if (rep.reactions == 0) { failures.push_back(id + ": no reaction"); continue; }
    for (const auto& [h, c] : formulas(s)) f[h] = std::max(f[h], c);
    for (const auto& p : ps) {
      const auto want = formulas(build_molecule(p).system);
      for (const auto& [formula, k] : want)
        if (!f.count(formula)) {
          std::string have;
          for (const auto& [h, c] : f) have += " " + h + "×" + std::to_string(c);
          failures.push_back(id + ": no " + formula + " among the products (" + std::to_string(rep.reactions) + " reactions:" + have + ")");
        }
    }
    ++ran;
  }
  for (const auto& x : failures) ADD_FAILURE() << x;
  std::printf("reaction library: %d converted, %d refused (three molecules at once), %d run on their model compounds\n", converted, refused, ran);
  EXPECT_GE(converted, 45);
  EXPECT_LE(refused, 4);
  EXPECT_GE(ran, 40);
}

// A grown cell keeps its repeat units through a LAMMPS data file (a "# res N NAME" comment LAMMPS ignores), so a cure run
// on the file still counts DC and places each link on its unit
TEST(React, RepeatUnitsSurviveADataFile) {
  ChainSpec spec;
  spec.units = {{"cis-1,4-butadiene", "[*]C/C=C\\C[*]"}};
  spec.dp = 6;
  GrowOptions g;
  g.chains = 2;
  g.density = 0.5;
  System s = grow_chains(spec, g);
  ASSERT_GT(s.atoms[0].resid, 0);
  const auto path = (std::filesystem::temp_directory_path() / "caps_res_test.data").string();
  write_lammps_data(s, path);
  const System back = read_lammps_data(path);
  ASSERT_EQ(back.atoms.size(), s.atoms.size());
  for (size_t i = 0; i < s.atoms.size(); ++i) {
    EXPECT_EQ(back.atoms[i].resid, s.atoms[i].resid);
    EXPECT_EQ(back.atoms[i].resname, s.atoms[i].resname);
  }
  std::filesystem::remove(path);
}

// A centre pushed exactly flat (an sp3 centre with an angle at 180°, as an opened ring or a moved hydrogen can leave one)
// sits on a saddle where the angle's force is zero; relax() moves it off the line and minimises again, so the angle
// comes back bent and dynamics does not blow up from it
TEST(React, RelaxLeavesAFlatSaddle) {
  // water laid out straight, H–O–H at 180°: every force lies along the line, so nothing but a kick bends it
  System s = build_molecule("O").system;
  ASSERT_EQ(s.atoms.size(), 3u);
  s.atoms[0].pos = {0, 0, 0};
  s.atoms[1].pos = {0.96, 0, 0};
  s.atoms[2].pos = {-0.96, 0, 0};
  RelaxOptions o;
  o.pushoff = false;
  o.ftol = 0.05;
  RelaxReport rep;
  relax(s, o, &rep);
  const Vec3 u = s.atoms[1].pos - s.atoms[0].pos, v = s.atoms[2].pos - s.atoms[0].pos;   // the oxygen is the vertex
  const double theta = std::acos(dot(u, v) / (norm(u) * norm(v))) * 180 / 3.14159265358979;
  EXPECT_LT(theta, 130.0);
  EXPECT_GT(theta, 100.0);
  bool noted = false;
  for (const auto& n : rep.notes) noted = noted || n.find("stuck at 180") != std::string::npos;
  EXPECT_TRUE(noted);
}

// Reactive sites per chain: counted for the templates, and capped during a run — every atom forming a bond on a chain uses
// one of its sites (both carbons of a C–C link), so with 2 per chain no chain reacts more than twice
TEST(React, SitesPerChainAreCountedAndCapped) {
  caps::System s = caps::open_file(std::string(CAPS_SAMPLES) + "/ps_melt.data").frame(0);
  const auto cc = caps::parse_templates(caps::builtin_template("cc_crosslink"));
  const auto per = caps::chain_sites(s, cc);
  ASSERT_EQ(per.size(), 10u);
  for (const auto& c : per) {
    int ch2 = 0;
    // the chain's CH2 carbons, counted directly
    caps::System t = s;
    t.has_mol = false;
    const auto mol = t.molecules();
    std::vector<int> h(s.atoms.size(), 0);
    for (const auto& b : s.bonds) { if (s.atoms[b.i].element == 1) ++h[b.j]; if (s.atoms[b.j].element == 1) ++h[b.i]; }
    for (size_t i = 0; i < s.atoms.size(); ++i)
      if (mol[i] + 1 == c.chain && s.atoms[i].element == 6 && h[i] == 2) {
        int deg = 0;
        for (const auto& b : s.bonds) deg += (b.i == i) + (b.j == i);
        ch2 += deg == 4;
      }
    EXPECT_EQ(c.sites, ch2);
  }
  caps::ReactOptions o;
  o.templates = cc;
  o.relax = false;
  o.max_cycles = 20;
  o.max_per_cycle = 4;
  o.between_chains = true;
  o.auto_capture = true;
  o.capture_max = 8;
  o.sites_per_chain = 2;
  o.target = caps::ReactTarget::Crosslinks;
  o.target_value = 30;   // more than 10 chains × 2 ÷ 2 = 10 allow
  caps::ReactReport rep;
  caps::System r = s;
  caps::react(r, o, &rep);
  EXPECT_LE(rep.crosslinks, 10);
  EXPECT_GT(rep.crosslinks, 0);
  EXPECT_TRUE(std::any_of(rep.notes.begin(), rep.notes.end(), [](const std::string& n) { return n.find("at most 10 links") != std::string::npos; }));
  // no chain lost more than 2 CH2 sites
  const auto after = caps::chain_sites(r, cc);
  int lost = 0;
  for (const auto& c : per) lost += c.sites;
  for (const auto& c : after) lost -= c.sites;
  EXPECT_EQ(lost, 2 * rep.crosslinks);
}

// ---- fix bond/react sets that LAMMPS runs as written (the PBS/ENR-50/MAH report): templates obey LAMMPS's rule that every
// atom is linked to an initiator without passing an edge atom, no template atom is bonded to nothing, a hydrogen that
// changes partner carries a distance constraint, every reaction has stabilize_steps, and the Rmax raise compares the
// absolute Rmax with its limit
namespace {
struct MolFile { int atoms = 0; std::vector<std::pair<int, int>> bonds; std::vector<std::string> type_names; };
MolFile read_mol_file(const std::string& path) {
  std::ifstream f(path);
  MolFile m;
  std::string line, sec;
  while (std::getline(f, line)) {
    const std::string t = line.substr(0, line.find('#'));
    std::istringstream ss(t);
    std::vector<std::string> w;
    for (std::string x; ss >> x;) w.push_back(x);
    if (w.empty()) continue;
    if (w.size() == 2 && w[1] == "atoms") { m.atoms = std::stoi(w[0]); continue; }
    if (w.size() == 1 && std::isalpha(static_cast<unsigned char>(w[0][0]))) { sec = w[0]; continue; }
    if (sec == "Bonds" && w.size() >= 4) m.bonds.push_back({std::stoi(w[2]), std::stoi(w[3])});
    if (sec == "Types" && w.size() >= 2) {
      const auto c = line.find('#');
      std::string name = c == std::string::npos ? "" : line.substr(c + 1);
      name.erase(0, name.find_first_not_of(' '));
      m.type_names.push_back(name);
    }
  }
  return m;
}
struct MapFile { std::vector<int> init, edge; std::vector<std::vector<std::string>> constraints; int equivalences = 0; };
MapFile read_map_file(const std::string& path) {
  std::ifstream f(path);
  MapFile m;
  std::string line, sec;
  while (std::getline(f, line)) {
    const std::string t = line.substr(0, line.find('#'));
    std::istringstream ss(t);
    std::vector<std::string> w;
    for (std::string x; ss >> x;) w.push_back(x);
    if (w.empty()) continue;
    if (w.size() == 1 && std::isalpha(static_cast<unsigned char>(w[0][0]))) { sec = w[0]; continue; }
    if (w.size() == 2 && std::isalpha(static_cast<unsigned char>(w[1][0]))) continue;   // a count line
    if (sec == "InitiatorIDs") m.init.push_back(std::stoi(w[0]));
    if (sec == "EdgeIDs") m.edge.push_back(std::stoi(w[0]));
    if (sec == "Equivalences") ++m.equivalences;
    if (sec == "Constraints") m.constraints.push_back(w);
  }
  return m;
}
// every template of a set: the LAMMPS rules, and the hydrogen-transfer constraints
void check_bond_react_set(const std::string& dir, const std::string& stem, const BondReactReport& rep, double hmax) {
  ASSERT_FALSE(rep.variants.empty());
  for (const auto& v : rep.variants) {
    const std::string base = dir + "/" + stem + "_" + v.name;
    const MolFile pre = read_mol_file(base + "_pre.mol"), post = read_mol_file(base + "_post.mol");
    const MapFile map = read_map_file(base + "_map.txt");
    EXPECT_EQ(pre.atoms, post.atoms) << v.name;
    EXPECT_EQ(map.equivalences, pre.atoms) << v.name;
    std::vector<std::set<int>> nb(size_t(pre.atoms) + 1), nbp(size_t(pre.atoms) + 1);
    for (auto [a, b] : pre.bonds) nb[size_t(a)].insert(b), nb[size_t(b)].insert(a);
    for (auto [a, b] : post.bonds) nbp[size_t(a)].insert(b), nbp[size_t(b)].insert(a);
    // reached from an initiator without passing an edge atom
    const std::set<int> edge(map.edge.begin(), map.edge.end());
    std::set<int> seen(map.init.begin(), map.init.end());
    std::vector<int> q(map.init.begin(), map.init.end());
    for (size_t h = 0; h < q.size(); ++h) {
      if (edge.count(q[h])) continue;
      for (int w : nb[size_t(q[h])]) if (seen.insert(w).second) q.push_back(w);
    }
    EXPECT_EQ(int(seen.size()), pre.atoms) << v.name << ": atoms beyond an edge";
    for (int a = 1; a <= pre.atoms && pre.atoms > 1; ++a) EXPECT_FALSE(nb[size_t(a)].empty()) << v.name << ": atom " << a << " bonded to nothing";
    // a hydrogen whose partner changes: one constraint to its new partner
    int moving = 0;
    for (int a = 1; a <= pre.atoms; ++a) {
      if (pre.type_names[size_t(a - 1)].rfind("H", 0) != 0 || nb[size_t(a)] == nbp[size_t(a)] || nbp[size_t(a)].empty()) continue;
      ++moving;
      int partner = 0;
      for (int w : nbp[size_t(a)]) if (!nb[size_t(a)].count(w)) partner = w;
      int found = 0;
      for (const auto& c : map.constraints)
        if (c.size() == 5 && c[0] == "distance" && std::stoi(c[1]) == a && std::stoi(c[2]) == partner && std::stod(c[4]) == hmax) ++found;
      EXPECT_EQ(found, 1) << v.name << ": hydrogen " << a << " → " << partner;
    }
    EXPECT_EQ(int(map.constraints.size()), moving) << v.name;
  }
  std::ifstream in(dir + "/" + stem + ".in");
  std::stringstream ss;
  ss << in.rdbuf();
  std::istringstream lines(ss.str());
  int reacts = 0;
  for (std::string l; std::getline(lines, l);)
    if (l.rfind("  react ", 0) == 0) { ++reacts; EXPECT_NE(l.find("stabilize_steps 200"), std::string::npos) << l; }
  EXPECT_EQ(reacts, int(rep.variants.size()));
}
System enr_chains(int chains, int dp, uint64_t seed) {
  ChainSpec spec;
  spec.units = {{"ENR", "[*]CC1(C)OC1C[*]"}, {"NR", "[*]C/C=C(C)\\C[*]"}};
  spec.sequence = Sequence::Alternating;
  spec.dp = dp;
  GrowOptions g;
  g.chains = chains;
  g.density = 0.3;
  g.seed = seed;
  g.auto_scale = true;
  return grow_chains(spec, g);
}
}  // namespace

TEST(React, BondReactMaleicAcidTwoStepsRunsInLammps) {
  System s = enr_chains(4, 6, 31);
  PackReport pr;
  s = insert_molecules(s, build_molecule("OC(=O)/C=C\\C(=O)O").system, 6, PackOptions{}, &pr);
  auto uff = [](const System& x) { return std::make_shared<const ForceField>(assign_uff(x)); };
  const auto dir = (std::filesystem::temp_directory_path() / "caps_br_maleic").string();
  std::filesystem::remove_all(dir);
  BondReactOptions bo;
  bo.between_chains = true;
  bo.survey_after = {0.25, 0.5, 0.8};
  bo.survey_relax = false;
  bo.survey_capture = 12;
  bo.mol_ids = "molmap";
  bo.targets = {0.5};
  bo.limiting = 6;
  bo.stall_chunks = 20;
  bo.rmax = 3.5, bo.rmax_step = 0.25, bo.rmax_limit = 4.0;
  const auto rep = write_bond_react(s, parse_templates(builtin_template("enr_acid_ester")), uff, dir, "ma", bo);
  check_bond_react_set(dir, "ma", rep, 3.5);
  if (std::getenv("CAPS_TEST_VERBOSE")) { for (const auto& n : rep.notes) std::printf("note: %s\n", n.c_str()); for (const auto& v : rep.variants) std::printf("%s %s %d atoms\n", v.name.c_str(), v.step.c_str(), v.pre_atoms); }
  EXPECT_GT(rep.h_transfers, 0);   // the acid H moves to the epoxide O
  // the Rmax raise compares the absolute Rmax with its limit and logs it
  std::ifstream in(dir + "/ma.in");
  std::stringstream ss;
  ss << in.rdbuf();
  EXPECT_NE(ss.str().find("$(3.5+v_rmax_add) < 4"), std::string::npos);
  EXPECT_EQ(ss.str().find("${rmax_add} < 4"), std::string::npos);
  EXPECT_NE(ss.str().find("Rmax now"), std::string::npos);
  // with LAMMPS on the PATH (a build with REACTER and reset_mol_ids molmap): a dense liquid of the model compounds (the ENR
  // epoxide unit and maleic acid) exported the same way runs 2000 steps with reactions and no atom lost (CAPS_TEST_LMP=0
  // skips it)
  const char* skip = std::getenv("CAPS_TEST_LMP");
  if (!(skip && std::string(skip) == "0") && std::system("command -v lmp > /dev/null 2>&1") == 0) {
    const System liq = packed({{"CC1(C)OC1C", 12}, {"OC(=O)/C=C\\C(=O)O", 6}}, 15, 7);
    const auto dl = (std::filesystem::temp_directory_path() / "caps_br_liquid").string();
    std::filesystem::remove_all(dl);
    BondReactOptions lo;
    lo.between_chains = true;
    lo.mol_ids = "molmap";
    lo.nevery = 20;
    lo.steps = 2000;
    lo.temperature = 400;
    const auto lrep = write_bond_react(liq, parse_templates(builtin_template("enr_acid_ester")), uff, dl, "liq", lo);
    check_bond_react_set(dl, "liq", lrep, 3.5);
    const int rc = std::system(("cd \"" + dl + "\" && lmp -in liq.in > check.log 2>&1").c_str());
    std::ifstream lg(dl + "/check.log");
    std::stringstream ls;
    ls << lg.rdbuf();
    const std::string log = ls.str();
    if (log.find("Unrecognized fix style 'bond/react'") != std::string::npos || (log.find("ERROR") != std::string::npos && log.find("molmap") != std::string::npos)) {
      GTEST_LOG_(INFO) << "this LAMMPS lacks REACTER or reset_mol_ids molmap: not run";
    } else {
      EXPECT_EQ(rc, 0) << log.substr(log.size() > 2000 ? log.size() - 2000 : 0);
      EXPECT_EQ(log.find("ERROR"), std::string::npos);
      EXPECT_EQ(log.find("Bond atoms"), std::string::npos);
      // at least one reaction: the sum of the f_rxns columns on the last thermo line
      std::istringstream ll(log);
      std::vector<std::string> head, last;
      for (std::string l; std::getline(ll, l);) {
        std::istringstream w(l);
        std::vector<std::string> t;
        for (std::string x; w >> x;) t.push_back(x);
        if (!t.empty() && t[0] == "Step") head = t;
        else if (!head.empty() && t.size() == head.size() && std::isdigit(static_cast<unsigned char>(t[0][0]))) last = t;
      }
      double reactions = 0;
      for (size_t c = 0; c < head.size() && c < last.size(); ++c)
        if (head[c].rfind("f_rxns[", 0) == 0) reactions += std::stod(last[c]);
      EXPECT_GT(reactions, 0) << "no reaction in 2000 steps";
    }
    if (!std::getenv("CAPS_TEST_KEEP")) std::filesystem::remove_all(dl);
  }
  if (!std::getenv("CAPS_TEST_KEEP")) std::filesystem::remove_all(dir);
}

TEST(React, BondReactAnhydrideAlcoholAndEpoxide) {
  System s = enr_chains(3, 6, 33);
  PackReport pr;
  s = insert_molecules(s, build_molecule("OCCCCO").system, 4, PackOptions{}, &pr);   // butanediol: PBS's hydroxyl ends
  s = insert_molecules(s, build_molecule("O=C1OC(=O)C=C1").system, 6, PackOptions{}, &pr);
  auto uff = [](const System& x) { return std::make_shared<const ForceField>(assign_uff(x)); };
  const auto dir = (std::filesystem::temp_directory_path() / "caps_br_mah").string();
  std::filesystem::remove_all(dir);
  BondReactOptions bo;
  bo.between_chains = true;
  bo.survey_after = {0.25, 0.5, 0.8};
  bo.survey_relax = false;
  bo.survey_capture = 12;
  bo.mol_ids = "molmap";
  const auto t = parse_templates(builtin_template("anhydride_alcohol") + "\n" + builtin_template("enr_acid_ester"));
  const auto rep = write_bond_react(s, t, uff, dir, "mah", bo);
  check_bond_react_set(dir, "mah", rep, 3.5);
  if (std::getenv("CAPS_TEST_VERBOSE")) { for (const auto& n : rep.notes) std::printf("note: %s\n", n.c_str()); for (const auto& v : rep.variants) std::printf("%s %s %d atoms\n", v.name.c_str(), v.step.c_str(), v.pre_atoms); }
  std::filesystem::remove_all(dir);
}
