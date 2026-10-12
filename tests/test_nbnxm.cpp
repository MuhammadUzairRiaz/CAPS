#include <gtest/gtest.h>

#include <cmath>
#include <random>
#include <string>
#include <vector>

#include "caps/dynamics.hpp"
#include "caps/field.hpp"
#include "caps/grow.hpp"
#include "caps/relax.hpp"

using namespace caps;

namespace {

// A small relaxed polystyrene cell (3 chains of DP 4, 0.9 g/cm³): narrower than twice the cut-offs below, so atoms meet
// their own periodic images
const System& ps_cell() {
  static const System s = [] {
    GrowOptions g;
    g.chains = 3;
    g.dp = 4;
    g.density = 0.4;
    g.seed = 3;
    System c = grow(g);
    RelaxOptions r;
    r.target_density = 0.9;
    r.ftol = 1.0;
    relax(c, r);
    return c;
  }();
  return s;
}

std::vector<double> positions(const System& s) {
  std::vector<double> x;
  for (const auto& a : s.atoms) x.insert(x.end(), a.pos.begin(), a.pos.end());
  return x;
}

// mixed precision against double: pair energies, forces and virial
void compare(const ForceField& ff, EnergyOptions eo, const System& s, const std::string& what) {
  eo.mixed = false;
  Evaluator ed(ff, eo);
  eo.mixed = true;
  Evaluator em(ff, eo);
  std::vector<double> x = positions(s);
  std::mt19937 rng(5);
  std::normal_distribution<double> g(0.0, 0.05);
  for (int round = 0; round < 3; ++round) {   // the start, then moved twice (the lists kept, pruned again, rebuilt)
    std::vector<double> fd, fm, fm2;
    const EnergyTerms d = ed.compute(x, s.cell, fd), m = em.compute(x, s.cell, fm);
    EXPECT_NEAR(m.vdw, d.vdw, 2e-6 * std::fabs(d.vdw) + 1e-4) << what << " round " << round;
    EXPECT_NEAR(m.coulomb, d.coulomb, 2e-6 * std::fabs(d.coulomb) + 1e-4) << what << " round " << round;
    EXPECT_NEAR(m.bond, d.bond, 1e-9 * std::fabs(d.bond) + 1e-9) << what;
    double rms = 0, err = 0;
    for (size_t k = 0; k < fd.size(); ++k) rms += fd[k] * fd[k], err += (fm[k] - fd[k]) * (fm[k] - fd[k]);
    // single-precision pair forces: a relaxed cell's net forces are small differences of large pair forces, about 5e-6 rms
    EXPECT_LT(std::sqrt(err / rms), 1e-5) << what << " round " << round;
    for (int c = 0; c < 6; ++c) EXPECT_NEAR(m.w[c], d.w[c], 2e-5 * std::fabs(d.virial) + 1e-3) << what << " w" << c;
    // forces and virial without the pair energies are the same
    em.energies(false);
    const EnergyTerms m2 = em.compute(x, s.cell, fm2);
    em.energies(true);
    for (size_t k = 0; k < fm.size(); ++k) ASSERT_EQ(fm2[k], fm[k]) << what << " component " << k;
    EXPECT_EQ(m2.virial, m.virial) << what;
    for (double& v : x) v += g(rng) * (round + 1);
  }
}

}  // namespace

TEST(Nbnxm, MixedPrecisionMatchesDoubleForEveryPairForm) {
  const System& s = ps_cell();
  ASSERT_TRUE(s.cell.valid());
  ForceField base = assign_gaff(s);
  struct Case { const char* form; const char* mixing; };
  for (const Case c : {Case{"lj12-6", "geometric"}, Case{"lj12-6", "arithmetic"}, Case{"lj9-6", "sixthpower"}}) {
    ForceField ff = base;
    ff.pair_form = c.form;
    ff.mixing = c.mixing;
    for (const bool pme : {true, false})
      for (const bool tail : {true, false}) {
        EnergyOptions eo;
        eo.cutoff = 8.0;
        eo.tail = tail;
        eo.threads = 3;
        eo.electrostatics = pme ? EnergyOptions::Electrostatics::PME : EnergyOptions::Electrostatics::DSF;
        compare(ff, eo, s, std::string(c.form) + " " + c.mixing + (pme ? " PME" : " DSF") + (tail ? " tail" : " shifted"));
      }
  }
  // no charges, one thread
  EnergyOptions eo;
  eo.cutoff = 9.0;
  eo.coulomb = false;
  eo.threads = 1;
  compare(base, eo, s, "no Coulomb");
}

TEST(Nbnxm, MixedPrecisionDynamicsConservesEnergy) {
  // NVE, 1 fs, 2000 steps: the conserved quantity drifts no more in mixed precision than in double
  auto drift = [](bool mixed) {
    System s = ps_cell();
    DynamicsOptions o;
    o.steps = 2000;
    o.dt = 1.0;
    o.temperature = 300;
    o.thermostat = Thermostat::None;
    o.thermo_every = 100;
    o.seed = 4;
    o.energy.cutoff = 8.0;
    o.energy.electrostatics = EnergyOptions::Electrostatics::PME;
    o.energy.mixed = mixed;
    o.energy.threads = 2;
    DynamicsReport r;
    run_dynamics(s, o, &r);
    double lo = 1e300, hi = -1e300;
    for (const auto& row : r.thermo) lo = std::min(lo, row.conserved), hi = std::max(hi, row.conserved);
    return hi - lo;
  };
  const double dd = drift(false), dm = drift(true);
  EXPECT_LT(dm, 1.5 * dd + 0.05) << "double " << dd << " mixed " << dm;
}
