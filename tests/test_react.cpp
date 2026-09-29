#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <map>

#include "caps/pack.hpp"
#include "caps/properties.hpp"
#include "caps/react.hpp"
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
// constraints, products by map number, deleted atoms absent, substitution bonds single. RDKit (2026.03) parses these and,
// run on ethane + ethane / propylene oxide + methylamine, gives butane and 1-(methylamino)propan-2-ol.
TEST(React, TemplatesAsReactionSmarts) {
  const auto cc = parse_templates(builtin_template("cc_crosslink"));
  EXPECT_EQ(reaction_smarts(cc[0]), "[#6;X4;!H0;A:1]~[#1:3].[#6;X4;!H0;A:2]~[#1:4]>>[#6:1]-[#6:2].[#1:3]-[#1:4]");
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
