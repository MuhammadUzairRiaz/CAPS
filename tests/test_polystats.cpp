// Polymer statistics: chain-length draws, the terminal copolymer model, stereo sequences, Flory–Huggins, Ewald sizes.
#include <gtest/gtest.h>

#include <cmath>
#include <numeric>

#include "caps/analysis.hpp"
#include "caps/entangle.hpp"
#include "caps/io.hpp"
#include "caps/kspace.hpp"
#include "caps/polymer.hpp"
#include "caps/polystats.hpp"

using namespace caps;

TEST(Polystats, SchulzZimmDrawsMatchTheTargetForLargeSamples) {
  const auto L = draw_chain_lengths("schulz-zimm", 40, 1.1, 20000, 7);
  double s1 = 0, s2 = 0;
  for (int n : L) s1 += n, s2 += double(n) * n;
  const double nn = s1 / double(L.size()), pdi = s2 / s1 / nn;
  EXPECT_NEAR(nn, 40, 0.3);
  EXPECT_NEAR(pdi, 1.1, 0.01);
  EXPECT_EQ(L, draw_chain_lengths("schulz-zimm", 40, 1.1, 20000, 7));   // deterministic for a seed
  const auto F = draw_chain_lengths("flory", 50, 2, 20000, 3);
  s1 = s2 = 0;
  for (int n : F) s1 += n, s2 += double(n) * n;
  EXPECT_NEAR(s2 / s1 / (s1 / double(F.size())), 2 - 1.0 / 50, 0.05);
  // the Gamma density integrates to 1
  double area = 0;
  for (double x = 0.05; x < 200; x += 0.1) area += chain_length_pdf("schulz-zimm", 40, 1.1, x) * 0.1;
  EXPECT_NEAR(area, 1, 1e-3);
}

TEST(Polystats, TerminalModelMatchesMayoLewis) {
  const auto m = copolymer_terminal(0.52, 0.46, 0.5);   // styrene / MMA
  EXPECT_NEAR(m.F1, 0.510, 5e-4);
  EXPECT_NEAR(m.azeotrope, 0.529, 5e-4);
  EXPECT_NEAR(m.run_a, 1.52, 5e-3);
  EXPECT_NEAR(m.run_b, 1.46, 5e-3);
  EXPECT_LT(copolymer_terminal(2, 0.5, 0.5).azeotrope, 0);   // one ratio above 1, one below: no azeotrope
  // long chains drawn with the terminal model reach the model's composition
  ChainSpec spec;
  spec.units.resize(2);
  spec.sequence = Sequence::Terminal;
  spec.dp = 40000, spec.r1 = 0.52, spec.r2 = 0.46, spec.weights = {0.5, 0.5};
  const auto seq = chain_sequence(spec, 12);
  EXPECT_NEAR(double(std::count(seq.begin(), seq.end(), 0)) / double(seq.size()), m.F1, 0.01);
}

TEST(Polystats, StereoModelsAndFits) {
  const auto b = stereo_bernoulli(0.5);
  EXPECT_NEAR(b.mm, 0.25, 1e-12);
  EXPECT_NEAR(b.mr, 0.5, 1e-12);
  EXPECT_NEAR(std::accumulate(b.pentads.begin(), b.pentads.end(), 0.0), 1, 1e-12);
  EXPECT_NEAR(b.pentads[0], 0.0625, 1e-12);   // mmmm = Pm⁴
  const auto iso = stereo_bernoulli(0.8);
  double rms = 1;
  EXPECT_NEAR(fit_bernoulli(iso.pentads, &rms).pm, 0.8, 1e-6);
  EXPECT_LT(rms, 1e-8);
  const auto mk = stereo_markov(0.2, 0.6);
  const auto fit = fit_markov(mk.pentads, &rms);
  EXPECT_NEAR(fit.p_mr, 0.2, 1e-5);
  EXPECT_NEAR(fit.p_rm, 0.6, 1e-5);
  const auto c = count_stereo("mmrmr");
  EXPECT_EQ(c.m, 3);
  EXPECT_EQ(c.mm, 1);
  EXPECT_EQ(c.mr, 3);
  EXPECT_EQ(c.pentad_total, 2);
}

TEST(Polystats, FloryHugginsCriticalPointSpinodalAndBinodal) {
  const auto c = blend_critical(100, 200);
  EXPECT_NEAR(c.chi_c, 0.01457, 1e-5);
  EXPECT_NEAR(c.phi_c, 0.5858, 1e-4);
  double lo, hi;
  ASSERT_TRUE(blend_spinodal(100, 200, 0.03, lo, hi));
  EXPECT_NEAR(lo, 0.1857, 1e-4);
  EXPECT_NEAR(hi, 0.8977, 1e-4);
  ASSERT_TRUE(blend_binodal(100, 200, 0.03, lo, hi));
  EXPECT_NEAR(lo, 0.0385, 1e-4);
  EXPECT_NEAR(hi, 0.9932, 1e-4);
  EXPECT_FALSE(blend_binodal(100, 200, 0.01, lo, hi));   // below χc: one phase
  // symmetric blend: binodal symmetric about ½
  ASSERT_TRUE(blend_binodal(50, 50, 0.05, lo, hi));
  EXPECT_NEAR(lo + hi, 1, 1e-9);
}

TEST(Polystats, EwaldAndSolventEstimates) {
  EXPECT_NEAR(ewald_beta(12, 1e-5), 0.2603, 1e-4);
  EXPECT_EQ(pme_mesh_size(45.3, 1.2), 40);
  EXPECT_EQ(pme_mesh_size(16.0, 1.2), 14);   // 2·7
  EXPECT_NEAR(hildebrand_chi(106.3, 18.2, 18.6, 298.15), 0.347, 1e-3);
}

// ---------------------------------------------------------------- model resolution, hydrogens from geometry

#include "caps/edit.hpp"
#include "caps/import.hpp"
#include "caps/molecule.hpp"
#include "caps/resolution.hpp"

TEST(Resolution, EicosaneKeepsItsMassAtEveryResolution) {
  BuildOptions b;
  b.forcefield = "uff";
  const System s = build_molecule("CCCCCCCCCCCCCCCCCCCC", b).system;
  ResolutionReport aa = all_atom_summary(s), ua, cg;
  const System u = united_atom(s, &ua);
  const System c = coarse_grain(s, 5, &cg);
  EXPECT_EQ(aa.sites, 62);
  EXPECT_EQ(int(u.atoms.size()), 20);
  EXPECT_EQ(int(c.atoms.size()), 4);
  EXPECT_NEAR(aa.mass, 282.556, 1e-3);
  EXPECT_NEAR(ua.mass, aa.mass, 1e-9);
  EXPECT_NEAR(cg.mass, aa.mass, 1e-9);
  EXPECT_EQ(int(u.bonds.size()), 19);
  EXPECT_EQ(int(c.bonds.size()), 3);
}

TEST(Resolution, HeavyAtomsGetTheirHydrogensFromTheGeometry) {
  BuildOptions b;
  b.forcefield = "uff";
  for (const auto& [smiles, h] : std::vector<std::pair<std::string, int>>{{"Cc1ccccc1", 8}, {"CC(=O)Nc1ccc(O)cc1", 9}, {"CN1C=NC2=C1C(=O)N(C)C(=O)N2C", 10}}) {
    System s = build_molecule(smiles, b).system;
    std::vector<char> heavy(s.atoms.size(), 0);
    for (size_t i = 0; i < s.atoms.size(); ++i) heavy[i] = s.atoms[i].element == 1;
    delete_atoms(s, heavy);
    for (auto& bd : s.bonds) bd.order = 1;   // as a PDB of heavy atoms reads
    orders_from_geometry(s);
    EXPECT_EQ(add_hydrogens(s), h) << smiles;
  }
}

// Chain A runs from (0, 1, 0) to (20, 1, 0) but loops behind chain B, which stands along z at x = 10, y = 0. The
// primitive path of A must stay hooked behind B; without B it pulls straight to 20 Å.
static System hooked_pair(bool with_b) {
  System s;
  auto add = [&](Vec3 p, int64_t mol) { Atom a; a.element = 6; a.mol = mol; a.id = int64_t(s.atoms.size()) + 1; a.pos = p; s.atoms.push_back(a); };
  for (int k = 0; k <= 40; ++k) {
    const double t = k / 40.0, x = 20 * t;
    add({x, 1 - 4 * std::sin(M_PI * t), 0}, 1);
    if (k) s.bonds.push_back({uint32_t(k - 1), uint32_t(k), 1});
  }
  if (with_b)
    for (int k = 0; k <= 20; ++k) {
      add({10, 0, -10.0 + k}, 2);
      if (k) s.bonds.push_back({uint32_t(s.atoms.size() - 2), uint32_t(s.atoms.size() - 1), 1});
    }
  s.has_mol = true;
  return s;
}

TEST(Entanglement, APathStaysHookedBehindAnotherChain) {
  for (bool with_b : {false, true}) {
    const System s = hooked_pair(with_b);
    std::vector<std::vector<uint32_t>> bb(1);
    for (uint32_t i = 0; i <= 40; ++i) bb[0].push_back(i);
    if (with_b) { bb.emplace_back(); for (uint32_t i = 41; i <= 61; ++i) bb[1].push_back(i); }
    PrimitivePathOptions o;
    o.sigma = 1.0;
    const PrimitivePaths p = primitive_paths(s, bb, o);
    ASSERT_EQ(p.lpp.size(), bb.size());
    if (!with_b) {
      EXPECT_NEAR(p.lpp[0], 20.0, 0.05);
    } else {
      // both paths bend where they meet (equal tension): A cannot pull straight, and B is pushed out of line
      EXPECT_GT(p.lpp[0], 20.1);
      EXPECT_GT(p.lpp[1], 20.02);
      double closest = 1e9;
      for (const auto& a : p.paths[0])
        for (const auto& b : p.paths[1]) closest = std::min(closest, norm(a - b));
      EXPECT_GT(closest, 0.9);   // the beads keep about σ apart
      // A's bead nearest x = 10 is still on the far side (y < 0) of B
      double ymid = 1e9, best = 1e9;
      for (const auto& q : p.paths[0]) if (std::abs(q[0] - 10) < best) best = std::abs(q[0] - 10), ymid = q[1];
      EXPECT_LT(ymid, 0.0);
      const EntanglementEstimate e = entanglement_estimate(p);
      EXPECT_EQ(e.chains, 2);
      EXPECT_GT(e.ne_mscoil, 0.0);
    }
  }
}

// samples/kg_melt.data: 20 Kremer–Grest chains of 100 beads (pushed off and run 200 τ in LAMMPS). The same primitive-path
// protocol run in LAMMPS (FENE K = 30, R₀ = 1.5 without the LJ part, WCA between chains only, neigh_modify exclude
// molecule/intra, ends fixed, Langevin T = 0.001, 200 000 steps of 0.006 τ) gives ⟨L_pp⟩ = 21.225 σ, N_e (modified
// S-coil) = 68.1 and N_e (classical) = 42.2.
TEST(Entanglement, KremerGrestMeltMatchesTheLammpsPrimitivePaths) {
  System s = read_lammps_data(std::string(CAPS_SAMPLES) + "/kg_melt.data");
  for (const auto& a : s.atoms) ASSERT_EQ(a.element, 6);   // mass-1 beads are read as beads, not hydrogens
  if (!s.unwrapped) make_molecules_whole(s);
  const auto bb = backbones(s);
  ASSERT_EQ(bb.size(), 20u);
  const PrimitivePaths p = primitive_paths(s, bb);
  EXPECT_TRUE(p.converged);
  const EntanglementEstimate e = entanglement_estimate(p);
  EXPECT_NEAR(e.lpp, 21.225, 0.01 * 21.225);
  EXPECT_NEAR(e.ne_mscoil, 68.1, 0.05 * 68.1);
  EXPECT_NEAR(e.ne_coil, 42.2, 0.03 * 42.2);
}

// Backmapping: beads of eicosane (5 backbone atoms each) turned and moved rigidly carry their atoms rigidly — exact for
// beads with two bonded neighbours; every bond inside a bead keeps its length; the wrong bead count is refused
TEST(Resolution, BackmapFollowsRigidlyMovedBeads) {
  BuildOptions b;
  b.forcefield = "uff";
  const System aa = build_molecule("CCCCCCCCCCCCCCCCCCCC", b).system;
  ResolutionReport rr;
  System beads = coarse_grain(aa, 5, &rr);
  ASSERT_EQ(beads.atoms.size(), 4u);
  // rotate 40° about (1, 2, 3) and shift
  const Vec3 ax = Vec3{1, 2, 3} * (1 / std::sqrt(14.0));
  const double th = 40 * M_PI / 180, c = std::cos(th), s = std::sin(th);
  auto rot = [&](const Vec3& v) { return v * c + cross(ax, v) * s + ax * (dot(ax, v) * (1 - c)); };
  for (auto& a : beads.atoms) a.pos = rot(a.pos) + Vec3{5, -3, 2};
  BackmapReport rep;
  const System back = backmap(aa, beads, 5, &rep);
  ASSERT_EQ(back.atoms.size(), aa.atoms.size());
  const auto nbb = beads.neighbours();
  for (size_t i = 0; i < aa.atoms.size(); ++i) {
    if (nbb[size_t(rr.site_of[i])].size() < 2) continue;   // end beads may roll about their one bond
    const Vec3 want = rot(aa.atoms[i].pos) + Vec3{5, -3, 2};
    EXPECT_NEAR(norm(back.atoms[i].pos - want), 0.0, 1e-6) << "atom " << i;
  }
  for (const auto& bd : aa.bonds)
    if (rr.site_of[bd.i] == rr.site_of[bd.j])
      EXPECT_NEAR(norm(back.atoms[bd.j].pos - back.atoms[bd.i].pos), norm(aa.atoms[bd.j].pos - aa.atoms[bd.i].pos), 1e-9);
  EXPECT_GT(rep.rms_turn, 25.0);   // the two middle beads turn 40°, the ends the smallest rotation that aligns their bond
  EXPECT_LE(rep.rms_turn, 40.0 + 1e-6);
  System three = beads;
  three.atoms.pop_back();
  EXPECT_THROW(backmap(aa, three, 5), std::invalid_argument);
}

// RIS by generator matrices: with equal weights the chain rotates freely (C₂ = 1 + cos θ, C∞ = (1 + cos θ)/(1 − cos θ));
// polyethylene at 413 K (Flory's parameters) gives C∞ ≈ 6.9, near the measured 6.7 ± 0.3
TEST(Polystats, RisGeneratorMatrices) {
  RisModel fr;
  fr.e_sigma = 0, fr.e_omega = 0;
  const auto a = ris_cn(fr, 300, 3000);
  const double c = std::cos((180 - 112) * M_PI / 180);
  EXPECT_NEAR(a[0], 1.0, 1e-12);
  EXPECT_NEAR(a[1], 1 + c, 1e-9);
  EXPECT_NEAR(a.back(), (1 + c) / (1 - c), 0.01);
  const auto pe = ris_cn(RisModel{}, 413, 3000);
  EXPECT_GT(pe.back(), 6.5);
  EXPECT_LT(pe.back(), 7.2);
  EXPECT_GT(ris_cn(RisModel{}, 300, 3000).back(), pe.back());   // colder: more trans, stiffer
}
