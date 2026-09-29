#include <gtest/gtest.h>

#include <cmath>
#include <cstdio>
#include <vector>

#include "caps/analysis.hpp"
#include "caps/cbmc.hpp"
#include "caps/field.hpp"
#include "caps/molecule.hpp"
#include "caps/polymer.hpp"
#include "caps/provenance.hpp"
#include "caps/recipe.hpp"
#include "caps/uff.hpp"
#include "caps/yaml.hpp"

using namespace caps;

namespace {
constexpr double kPi = 3.14159265358979323846;

double dihedral(const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& d) {
  const Vec3 b1 = b - a, b2 = c - b, b3 = d - c;
  const Vec3 m = cross(b1, b2), n = cross(b2, b3);
  return std::atan2(norm(b2) * dot(b1, n), dot(m, n));
}

Vec3 turn(const Vec3& p, const Vec3& c, const Vec3& axis, double t) {
  const Vec3 u = axis * (1 / norm(axis)), v = p - c;
  return c + v * std::cos(t) + cross(u, v) * std::sin(t) + u * (dot(u, v) * (1 - std::cos(t)));
}

double evaluate(Evaluator& ev, const System& s) {
  std::vector<double> x(3 * s.atoms.size()), f;
  for (size_t i = 0; i < s.atoms.size(); ++i)
    for (int d = 0; d < 3; ++d) x[3 * i + size_t(d)] = s.atoms[i].pos[size_t(d)];
  const auto t = ev.compute(x, s.cell, f);
  return t.total();
}
}  // namespace

// Butane in vacuum: the trans population CBMC samples equals the Boltzmann one, integrated over the central torsion and
// both methyl rotations on a grid with the Evaluator's energy (van der Waals and torsions, no electrostatics)
TEST(Cbmc, ButaneTorsionsFollowBoltzmann) {
  System s = build_molecule("CCCC").system;
  const ForceField ff = default_forcefield(s);
  std::vector<uint32_t> C;
  for (uint32_t i = 0; i < s.atoms.size(); ++i)
    if (s.atoms[i].element == 6) C.push_back(i);
  ASSERT_EQ(C.size(), 4u);
  const auto nb = s.neighbours();
  auto hydrogens = [&](uint32_t c) { std::vector<uint32_t> h; for (uint32_t w : nb[c]) if (s.atoms[w].element == 1) h.push_back(w); return h; };
  const auto h1 = hydrogens(C[0]), h3 = hydrogens(C[2]), h4 = hydrogens(C[3]);
  // the exact trans fraction at 300 K
  EnergyOptions eo;
  eo.coulomb = false;
  eo.tail = false;
  eo.threads = 1;
  Evaluator ev(ff, eo);
  const double beta = 1 / (0.0019872043 * 300);
  double z = 0, zt = 0, emin = 1e300;
  std::vector<std::array<double, 2>> grid;
  for (int a = 0; a < 72; ++a)
    for (int b = 0; b < 12; ++b)
      for (int c = 0; c < 12; ++c) {
        System t = s;
        auto& P = t.atoms;
        for (uint32_t h : h1) P[h].pos = turn(P[h].pos, P[C[0]].pos, P[C[0]].pos - P[C[1]].pos, 2 * kPi * b / 12);
        std::vector<uint32_t> far = h4;
        far.insert(far.end(), h3.begin(), h3.end());
        far.push_back(C[3]);
        for (uint32_t i : far) P[i].pos = turn(P[i].pos, P[C[2]].pos, P[C[2]].pos - P[C[1]].pos, 2 * kPi * a / 72);
        for (uint32_t h : h4) P[h].pos = turn(P[h].pos, P[C[3]].pos, P[C[3]].pos - P[C[2]].pos, 2 * kPi * c / 12);
        const double e = evaluate(ev, t), phi = dihedral(P[C[0]].pos, P[C[1]].pos, P[C[2]].pos, P[C[3]].pos);
        grid.push_back({e, phi});
        emin = std::min(emin, e);
      }
  for (const auto& [e, phi] : grid) {
    const double w = std::exp(-beta * (e - emin));
    z += w;
    if (std::fabs(phi) > 2 * kPi / 3) zt += w;
  }
  const double exact = zt / z;
  // CBMC
  CbmcOptions o;
  o.coulomb = false;
  o.moves = 4;
  o.trials = 6;
  o.max_torsions = 3;
  int trans = 0, samples = 0, acc = 0, att = 0;
  for (int it = 0; it < 6000; ++it) {
    o.seed = uint64_t(it) + 1;
    CbmcReport r;
    cbmc_regrow(s, ff, o, &r);
    acc += r.accepted, att += r.attempted;
    trans += std::fabs(dihedral(s.atoms[C[0]].pos, s.atoms[C[1]].pos, s.atoms[C[2]].pos, s.atoms[C[3]].pos)) > 2 * kPi / 3;
    ++samples;
  }
  const double sampled = double(trans) / samples;
  std::printf("butane: exact trans %.4f, CBMC %.4f, acceptance %.3f\n", exact, sampled, double(acc) / att);
  EXPECT_GT(exact, 0.5);
  EXPECT_LT(exact, 0.9);
  EXPECT_NEAR(sampled, exact, 0.025) << "exact " << exact << " sampled " << sampled << " acceptance " << double(acc) / att;
  EXPECT_GT(double(acc) / att, 0.3);
}

// A polyethylene melt: bond lengths and angles are kept, the trial-energy model's change over the accepted moves is the
// Evaluator's energy change (van der Waals shifted at the same cut-off, torsions; no electrostatics), and a cell with no
// chain end is refused
TEST(Cbmc, MeltKeepsGeometryAndCountsEnergy) {
  ChainSpec c;
  c.units.push_back({"*CC*", "*CC*"});
  c.dp = 20;
  GrowOptions g;
  g.chains = 6;
  g.density = 0.6;
  g.seed = 3;
  System s = grow_chains(c, g);
  const ForceField ff = default_forcefield(s);
  CbmcOptions o;
  o.coulomb = false;
  o.moves = 400;
  o.cutoff = 9;
  o.seed = 7;
  EnergyOptions eo;
  eo.coulomb = false;
  eo.tail = false;
  eo.cutoff = 9;
  eo.threads = 1;
  make_molecules_whole(s);
  const System before = s;
  Evaluator ev(ff, eo);
  const double e0 = evaluate(ev, s);
  CbmcReport r;
  cbmc_regrow(s, ff, o, &r);
  const double e1 = evaluate(ev, s);
  std::printf("melt: %d/%d accepted, dE Evaluator %.6f, CBMC %.6f, R2 %.1f -> %.1f, %.2f s\n", r.accepted, r.attempted, e1 - e0, r.energy_change, r.r2_before, r.r2_after, r.seconds);
  EXPECT_EQ(r.chains, 12);
  EXPECT_GT(r.accepted, 20);
  EXPECT_LT(r.accepted, r.attempted);
  EXPECT_NEAR(e1 - e0, r.energy_change, 1e-6 * std::max(1.0, std::fabs(e1 - e0))) << "Evaluator " << e1 - e0 << " CBMC " << r.energy_change;
  double worst = 0;
  for (const auto& b : s.bonds) {
    const double l0 = norm(before.cell.minimum_image(before.atoms[b.j].pos - before.atoms[b.i].pos));
    const double l1 = norm(s.cell.minimum_image(s.atoms[b.j].pos - s.atoms[b.i].pos));
    worst = std::max(worst, std::fabs(l1 - l0));
  }
  EXPECT_LT(worst, 1e-9);
  System ring = build_molecule("c1ccccc1").system;
  EXPECT_THROW(cbmc_regrow(ring, default_forcefield(ring), o), std::invalid_argument);
}

// A recipe stage: grown, relaxed, regrown; the manifest records the step and the methods text cites Siepmann & Frenkel
TEST(Cbmc, RecipeStage) {
  RecipeOptions o;
  const auto res = run_recipe(yaml_parse("build: {polymer: {smiles: \"*CC*\", dp: 10, chains: 4}}\ntype: {forcefield: default}\ngrow: {density: 0.5, seed: 2}\n"
                                         "relax: {fmax: 5}\ncbmc: {moves: 100, trials: 6, seed: 4}\n"), o);
  bool step = false;
  for (const auto& st : res.manifest.steps) step |= st.engine == "cbmc.regrow";
  EXPECT_TRUE(step);
  EXPECT_NE(methods_text(res.manifest).find("configurational-bias Monte Carlo"), std::string::npos);
}

// Recipes build crystals, surfaces, nanostructures and solvent boxes
TEST(Recipe, BuildCrystalSurfaceNanoSolvate) {
  RecipeOptions o;
  const std::string cif = std::string(CAPS_SOURCE_DIR) + "/data/crystals/";
  const auto rut = run_recipe(yaml_parse("build: {crystal: {group: \"P 42/m n m\", cell: [4.594, 4.594, 2.959], sites: \"Ti1 Ti 0 0 0; O1 O 0.3048 0.3048 0\", "
                                         "supercell: [2, 2, 3]}}\n"), o);
  EXPECT_EQ(rut.system.atoms.size(), 72u);
  const auto slab = run_recipe(yaml_parse("build: {surface: {cif: \"" + cif + "alpha-quartz.cif\", hkl: [0, 0, 1], layers: 3, vacuum: 15}}\n"), o);
  EXPECT_GT(slab.system.atoms.size(), 20u);
  const auto tube = run_recipe(yaml_parse("build: {nano: {kind: tube, n: 5, m: 5, length: 10}}\n"), o);
  EXPECT_GT(tube.system.atoms.size(), 50u);
  const auto water = run_recipe(yaml_parse("build: {solvate: {solvent: water, water_model: \"SPC/E\", edge: 16, ions: none, seed: 2}}\n"), o);
  EXPECT_NEAR(water.system.density(), 1.0, 0.03);   // SPC/E water at its density
  EXPECT_THROW(run_recipe(yaml_parse("build: {crystal: {group: \"P 1\"}}\n"), o), RecipeError);
}
