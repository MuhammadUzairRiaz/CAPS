#include <gtest/gtest.h>

#include <cmath>
#include <set>

#include "caps/graft_chains.hpp"
#include "caps/polymer.hpp"

using namespace caps;

namespace {
int hydrogens(const System& s, const std::vector<std::vector<uint32_t>>& nb, uint32_t i) {
  int h = 0;
  for (uint32_t j : nb[i]) h += s.atoms[j].element == 1;
  return h;
}
System pe_chain(int dp) {
  ChainSpec pe;
  pe.units = {{"PE", "[*]CC[*]"}};
  pe.dp = dp;
  GrowOptions g;
  g.chains = 1, g.density = 0.05, g.seed = 3;
  return grow_chains(pe, g);
}
}  // namespace

// Grafts grown with the backbone: polyisoprene units with the methine carbon marked [CH:3] carry PMMA grafts there
// only, at most graft_max per chain; each graft's first atom is bonded to a marked carbon that has given up its hydrogen
TEST(Graft, GraftsGrowOnTheMarkedAtom) {
  ChainSpec spec;
  spec.units = {{"isoprene", "[*]CC(C)=[CH:3]C[*]"}};
  spec.dp = 20;
  spec.architecture = Architecture::Graft;
  spec.graft_units = {{"MMA", "[*]CC(C)(C(=O)OC)[*]"}};
  spec.graft_dp = 4;
  spec.graft_probability = 0.3;
  spec.graft_max = 3;
  GrowOptions g;
  g.chains = 2, g.density = 0.1, g.seed = 5;
  g.auto_scale = true;
  GrowReport rep;
  const System s = grow_chains(spec, g, &rep);
  const int grafts = int(rep.sequences.size()) - 2;
  ASSERT_GT(grafts, 0);
  EXPECT_LE(grafts, 2 * spec.graft_max);
  const auto nb = s.neighbours();
  int links = 0;
  for (const auto& b : s.bonds) {
    const auto &a = s.atoms[b.i], &c = s.atoms[b.j];
    if ((a.resname == "ISO") == (c.resname == "ISO")) continue;
    const uint32_t bb = a.resname == "ISO" ? b.i : b.j;
    ++links;
    EXPECT_EQ(s.atoms[bb].element, 6);
    EXPECT_EQ(hydrogens(s, nb, bb), 0) << "the backbone atom of a graft kept a hydrogen";
    // the marked carbon: its heavy neighbours are three carbons (=C(CH3), the next CH2 and the graft)
    int carbons = 0;
    for (uint32_t j : nb[bb]) carbons += s.atoms[j].element == 6;
    EXPECT_EQ(carbons, 3);
  }
  EXPECT_EQ(links, grafts);
  std::printf("%d grafts of PMMA on 2 polyisoprene chains (%zu atoms)\n", grafts, s.atoms.size());
}

// The atoms equivalent to a CH2 of polyethylene are its 58 CH2 carbons (the two CH3 ends differ)
TEST(Graft, EquivalentSites) {
  const System pe = pe_chain(30);
  const auto nb = pe.neighbours();
  uint32_t ch2 = 0;
  for (uint32_t i = 0; i < pe.atoms.size(); ++i)
    if (pe.atoms[i].element == 6 && hydrogens(pe, nb, i) == 2) { ch2 = i; break; }
  EXPECT_EQ(equivalent_atoms(pe, ch2).size(), 58u);
}

// Grafting onto an existing chain: four grafts at least 5 Å apart on equivalent CH2 sites; the host's atoms stay where
// they were, each site loses one hydrogen and gains one graft, and the whole is one molecule; a host without a cell
// is grafted in a box and given back without one
TEST(Graft, OntoAnExistingChain) {
  const System pe = pe_chain(30);
  const auto nb = pe.neighbours();
  uint32_t ch2 = 0;
  for (uint32_t i = 0; i < pe.atoms.size(); ++i)
    if (pe.atoms[i].element == 6 && hydrogens(pe, nb, i) == 2) { ch2 = i; break; }
  ChainSpec g;
  g.units = {{"MMA", "[*]CC(C)(C(=O)OC)[*]"}};
  g.dp = 3;
  ChainGraftOptions o;
  o.sites = {ch2};
  o.equivalent = true;
  o.max_grafts = 4;
  o.min_spacing = 5.0;
  o.grow.seed = 2;
  for (bool cell : {true, false}) {
    System host = pe;
    if (!cell) host.cell = Cell{};
    ChainGraftReport rep;
    const System s = graft_chains(host, g, o, &rep);
    EXPECT_EQ(rep.grafted, 4);
    EXPECT_EQ(rep.candidates, 58);
    EXPECT_EQ(s.cell.valid(), cell);
    // the host minus four hydrogens first, at their places (a periodic host's atoms may come back wrapped)
    auto dist = [&](const Vec3& a, const Vec3& b) { return norm(host.cell.valid() ? host.cell.minimum_image(a - b) : a - b); };
    size_t k = 0, dropped = 0;
    for (size_t i = 0; i < host.atoms.size(); ++i) {
      ASSERT_LT(k, s.atoms.size());
      if (dist(s.atoms[k].pos, host.atoms[i].pos) > 1e-6) {
        EXPECT_EQ(host.atoms[i].element, 1) << "host atom " << i << " moved";
        ++dropped;
        continue;
      }
      ++k;
    }
    EXPECT_EQ(dropped, 4u);
    std::set<uint32_t> sites(rep.sites.begin(), rep.sites.end());
    EXPECT_EQ(sites.size(), 4u);
    for (uint32_t a : sites) {
      for (uint32_t b : sites)
        if (a < b) EXPECT_GE(norm(host.cell.valid() ? host.cell.minimum_image(host.atoms[a].pos - host.atoms[b].pos) : host.atoms[a].pos - host.atoms[b].pos), 5.0);
    }
    int nmol = 0;
    (void)s.molecules(&nmol);
    EXPECT_EQ(nmol, 1);
    for (const auto& a : s.atoms) EXPECT_EQ(a.mol, 1);
    std::printf("graft onto PE (%s): %d grafts, %zu atoms\n", cell ? "periodic" : "no cell", rep.grafted, s.atoms.size());
  }
}
