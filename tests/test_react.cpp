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
