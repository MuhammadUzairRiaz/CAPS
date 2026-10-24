#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "caps/exchange_moves.hpp"
#include "caps/grow.hpp"
#include "caps/molecule.hpp"
#include "caps/pack.hpp"
#include "caps/polymer.hpp"
#include "caps/react.hpp"
#include "caps/uff.hpp"

using namespace caps;

namespace {

std::vector<double> coords(const System& s) {
  std::vector<double> x;
  for (const auto& a : s.atoms) x.insert(x.end(), {a.pos[0], a.pos[1], a.pos[2]});
  return x;
}

// Chains of the given units grown into a cell at a low density (room for the curatives inserted next).
System chains(const std::vector<RepeatUnit>& units, Sequence seq, int n, int dp, uint64_t seed) {
  ChainSpec spec;
  spec.units = units;
  spec.sequence = seq;
  spec.dp = dp;
  GrowOptions g;
  g.chains = n;
  g.density = 0.5;
  g.seed = seed;
  g.auto_scale = true;
  return grow_chains(spec, g);
}

// Molecules from SMILES packed into a cubic cell.
System packed(const std::vector<std::pair<std::string, int>>& mols, double box, uint64_t seed) {
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

System react_with(System s, const std::string& templ) {
  ReactOptions r;
  r.templates = parse_templates(builtin_template(templ));
  r.relax = false;
  r.auto_capture = true;
  r.max_cycles = 30;
  ReactReport rep;
  react(s, r, &rep);
  EXPECT_GT(rep.reactions, 0) << templ;
  return s;
}

System cured(System s, const std::string& curative, int count, const std::string& templ) {
  s = insert_molecules(s, build_molecule(curative).system, count, PackOptions{});
  ReactOptions r;
  r.templates = parse_templates(builtin_template(templ));
  r.relax = false;
  r.max_cycles = 30;
  ReactReport rep;
  react(s, r, &rep);
  EXPECT_GT(rep.reactions, 0) << templ;
  return s;
}

struct Run {
  ExchangeStats stats;
  size_t candidates = 0;
};

// Exchange moves at the network's positions, then checks: every hydrogen holds one bond, every atom keeps its valence,
// and the run's types are a fresh typing's of the structure as it now is (types followed roles).
Run exchange(System& s, const std::string& templ, int tries, double temperature, uint64_t seed) {
  ForceField ff = assign_uff(s);
  ExchangeOptions xo;
  xo.kind = ExchangeKind::Template;
  xo.templates = parse_templates(builtin_template(templ));
  xo.relax_steps = 150;
  xo.seed = seed;
  xo.parameterize = [](const System& f) { return assign_uff(f); };   // new combinations of types learned
  std::vector<int> valence(s.atoms.size(), 0);
  for (const auto& b : s.bonds) ++valence[b.i], ++valence[b.j];
  ExchangeEngine e(s, ff, xo, EnergyOptions{});
  Run r;
  auto t = xo.templates[0];
  if (t.min_path <= 0) t.min_path = 4;
  r.candidates = find_matches(s, t).size();
  auto x = coords(s);
  for (int k = 0; k < tries; ++k) e.run(x, s.cell, temperature, nullptr);
  r.stats = e.stats();
  std::vector<int> now(s.atoms.size(), 0);
  for (const auto& b : s.bonds) ++now[b.i], ++now[b.j];
  EXPECT_EQ(now, valence);
  for (size_t i = 0; i < s.atoms.size(); ++i) s.atoms[i].pos = {x[3 * i], x[3 * i + 1], x[3 * i + 2]};
  const ForceField fresh = assign_uff(s);
  EXPECT_EQ(ff.atom_type, fresh.atom_type);
  EXPECT_EQ(ff.bonds.size(), fresh.bonds.size());
  EXPECT_EQ(ff.angles.size(), fresh.angles.size());
  EXPECT_EQ(ff.dihedrals.size(), fresh.dihedrals.size());
  EXPECT_GT(r.candidates, 0u) << templ;
  EXPECT_EQ(r.stats.unparameterised, 0) << templ << ": " << r.stats.note;   // every new combination of types learned
  EXPECT_EQ(r.stats.relaxed, r.stats.attempts - r.stats.no_candidate) << templ;
  return r;
}

}  // namespace

// 1. ENR-50 cured with 4,4'-dithiodibutyric acid (Imbernon et al. 2015): β-hydroxy esters from the epoxide opening, and
// disulfides in the crosslinks: both exchange.
TEST(RubberVitrimers, EnrDiacid) {
  System s = packed({{"CC1(C)OC1C", 12}, {"OC(=O)CCCSSCCCC(=O)O", 6}}, 17, 3);   // ENR's epoxidised unit as its model compound
  s = react_with(s, "enr_acid_ester");
  System t = s;
  EXPECT_GT(exchange(s, "transesterification_exchange", 30, 450, 1).stats.accepted, 0);
  exchange(t, "polysulfide_exchange", 30, 450, 2);   // the crosslinks' disulfides are found; at these positions they stay
}

// 3. Natural rubber cured with sulfur donors (H–S4–H: polysulfide bridges): S–S metathesis keeps every sulfur.
TEST(RubberVitrimers, SulfurCuredNr) {
  System s = chains({{"cis-1,4-isoprene", "[*]C/C=C(C)\\C[*]"}}, Sequence::Alternating, 4, 10, 5);
  s = cured(s, "SSSS", 8, "sulfur_allylic");
  const auto r = exchange(s, "polysulfide_exchange", 30, 450, 3);
  EXPECT_GT(r.stats.accepted, 0);
  EXPECT_GT(r.stats.learned, 0);   // the first C–S–S–C torsion: no disulfide bridge before the exchange
}

namespace {

// Bonds from a carbonyl carbon (an ester's acyl group) to an oxygen of the first `n` atoms (one component's molecules,
// packed first).
int acyl_on(const System& s, uint32_t n) {
  std::vector<std::vector<uint32_t>> nb(s.atoms.size());
  for (const auto& b : s.bonds) nb[b.i].push_back(b.j), nb[b.j].push_back(b.i);
  int k = 0;
  for (uint32_t o = 0; o < n; ++o) {
    if (s.atoms[o].element != 8) continue;
    for (uint32_t c : nb[o]) {
      if (c < n || s.atoms[c].element != 6) continue;
      bool carbonyl = false;
      for (uint32_t q : nb[c]) carbonyl = carbonyl || (s.atoms[q].element == 8 && nb[q].size() == 1);
      k += carbonyl;
    }
  }
  return k;
}

}  // namespace

// 2. ENR on cellulose: the β-hydroxy esters of an acid-cured ENR pass their acyl groups to the sugar's hydroxyls
// (glucose stands for the cellulose surface), the bond between the rubber and the fibre.
TEST(RubberVitrimers, EnrOnCellulose) {
  System s = packed({{"OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O", 6}, {"CC1(C)OC1C", 8}, {"OC(=O)CCCC(=O)O", 4}}, 17, 7);
  const uint32_t sugar = 6 * 24;   // the glucose atoms come first
  s = react_with(s, "enr_acid_ester");
  EXPECT_EQ(acyl_on(s, sugar), 0);
  EXPECT_GT(exchange(s, "transesterification_exchange", 60, 450, 4).stats.accepted, 0);
  EXPECT_GT(acyl_on(s, sugar), 0);
}

// 4. ENR/PBS through maleic acid: PBS's own ester links exchange with ENR's β-hydroxyls (grafts), next to the
// ENR–maleate–ENR bridges.
TEST(RubberVitrimers, EnrPbsMaleate) {
  System s = packed({{"OCCCCOC(=O)CCC(=O)OCCCCO", 4}, {"CC1(C)OC1C", 8}, {"OC(=O)/C=C\\C(=O)O", 3}}, 17, 9);
  s = react_with(s, "enr_acid_ester");
  EXPECT_GT(exchange(s, "transesterification_exchange", 60, 450, 6).stats.accepted, 0);
}
