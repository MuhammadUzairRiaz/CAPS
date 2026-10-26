#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <memory>

#include "caps/exchange_protocols.hpp"
#include "caps/kremer_grest.hpp"

using namespace caps;

namespace {

// Straight bead–spring chains along x on a grid in y and z (σ = 1 Å, ε = 1 kcal/mol, m = 1 g/mol), their ends 2.97 Å
// apart through the boundary so that no bond and no capture reaches across it.
System chains(int side, int beads, double spacing = 1.15) {
  System s;
  s.cell.a = {beads * 0.97 + 2.0, 0, 0};
  s.cell.b = {0, side * spacing, 0};
  s.cell.c = {0, 0, side * spacing};
  for (int cy = 0; cy < side; ++cy)
    for (int cz = 0; cz < side; ++cz) {
      const uint32_t first = uint32_t(s.atoms.size());
      for (int k = 0; k < beads; ++k) {
        Atom a;
        a.mol = 1 + cy * side + cz;
        a.type = 1;
        a.pos = {1.0 + 0.97 * k, (cy + 0.5) * spacing, (cz + 0.5) * spacing};
        s.atoms.push_back(a);
        if (k) s.bonds.push_back({first + uint32_t(k) - 1, first + uint32_t(k), 1});
      }
    }
  s.has_mol = true;
  return s;
}

}  // namespace

// A cut takes the bonds crossing the plane, a gap moves the upper side and lengthens the cell, a void takes atoms with
// their bonds, and fragments are cuts at mid-planes.
TEST(ExchangeProtocols, Builders) {
  System s = chains(4, 12);
  EXPECT_EQ(cut_bonds_across(s, 0, 0.5), 16u);   // one bond per chain crosses the mid-plane
  EXPECT_EQ(cut_bonds_across(s, 0, 0.5), 0u);
  const double lx = s.cell.a[0];
  open_gap(s, 0, 0.5, 3.0);
  EXPECT_NEAR(s.cell.a[0], lx + 3.0, 1e-12);
  EXPECT_NEAR(s.atoms[11].pos[0], 1.0 + 0.97 * 11 + 3.0, 1e-12);   // above the plane: moved
  EXPECT_NEAR(s.atoms[0].pos[0], 1.0, 1e-12);                      // below: not
  System f = chains(4, 12);
  EXPECT_EQ(fragment(f, {0, 1}), 16u);   // across x; nothing crosses y (the chains run along x)
  System v = chains(4, 12);
  const size_t before = v.bonds.size();
  EXPECT_EQ(remove_sphere(v, v.atoms[5].pos, 0.5), 1u);
  EXPECT_EQ(v.atoms.size(), size_t(4 * 4 * 12 - 1));
  EXPECT_EQ(v.bonds.size(), before - 2);
  for (const auto& b : v.bonds) EXPECT_LT(std::max(b.i, b.j), v.atoms.size());
  EXPECT_THROW(cut_bonds_across(v, 3, 0.5), std::invalid_argument);
}

// A cut healed by bind moves (a bead holds two bonds, so the cut ends have one free; a binding energy of many k_B T):
// bonds bridging the two sides come back.
TEST(ExchangeProtocols, HealRunRebondsACut) {
  System s = chains(4, 12);
  const auto side = sides_from_positions(s, 0, 0.5);
  ASSERT_EQ(cut_bonds_across(s, 0, 0.5), 16u);
  KgOptions ko;
  ProtocolOptions o;
  o.field = std::make_shared<ForceField>(kremer_grest_forcefield(s, ko));
  o.temperature = 503.2, o.dt = 0.25, o.tau_t = 25;
  o.moves.kind = ExchangeKind::Bind;
  for (uint32_t i = 0; i < s.atoms.size(); ++i) o.moves.group_a.push_back(i);
  o.moves.valence_a = 2;
  o.moves.binding = 30;
  o.moves.capture = 1.3;
  o.moves.attempts = 20;
  EXPECT_EQ(interface_bonds(s, side, 0, 0.5, 1.0).bridges, 0u);
  const auto r = heal_run(s, o, side, 0, 0.5, 0.5, 2);
  ASSERT_EQ(r.bridges.size(), 2u);
  EXPECT_NEAR(r.t_ps.back(), 0.5, 1e-12);
  EXPECT_GT(r.bridges.back(), 0);
  EXPECT_GT(r.stats.accepted, 0);
}

// A step strain stretches the cell by 1 + ε along the axis and keeps the volume, and the stretched chains pull; a Tv scan
// gives one specific volume per temperature.
TEST(ExchangeProtocols, StepStrainAndTvScan) {
  System s = chains(4, 12);
  KgOptions ko;
  ProtocolOptions o;
  o.field = std::make_shared<ForceField>(kremer_grest_forcefield(s, ko));
  o.temperature = 503.2, o.dt = 0.25, o.tau_t = 25;
  o.exchange = false;
  const double v0 = s.cell.volume(), lx = s.cell.a[0];
  const auto r = step_strain(s, o, 0, 0.05, 0.25, 50);
  EXPECT_NEAR(s.cell.a[0], lx * 1.05, 1e-9);
  EXPECT_NEAR(s.cell.volume(), v0, 1e-9 * v0);
  ASSERT_FALSE(r.stress.empty());
  EXPECT_EQ(r.stress.size(), r.t_ps.size());
  EXPECT_GT(r.stress.front(), 0);
  const auto t = tv_scan(s, o, 400, 600, 100, 0.25);
  ASSERT_EQ(t.T.size(), 3u);
  EXPECT_DOUBLE_EQ(t.T.back(), 600);
  for (double v : t.specific_volume) EXPECT_GT(v, 0);
}

// An error in a worker thread comes back to the caller: a FENE bond stretched past R0, evaluated on four threads, throws
// FieldError (it used to end the process).
TEST(ExchangeProtocols, WorkerErrorsComeBackToTheCaller) {
  System s = chains(4, 12);
  KgOptions ko;
  const ForceField ff = kremer_grest_forcefield(s, ko);
  s.atoms[1].pos[0] += 2.0;   // its bonds to atoms 0 and 2 past R0 = 1.5 Å
  std::vector<double> x, f;
  for (const auto& a : s.atoms) x.insert(x.end(), {a.pos[0], a.pos[1], a.pos[2]});
  for (int threads : {1, 4}) {
    EnergyOptions eo;
    eo.threads = threads;
    Evaluator ev(ff, eo);
    EXPECT_THROW(ev.compute(x, s.cell, f), FieldError) << threads << " threads";
  }
}
