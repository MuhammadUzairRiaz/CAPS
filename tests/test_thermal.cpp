// Semicrystalline lamellar cells and MD thermal conductivity.
#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <numeric>
#include <random>
#include <set>
#include <string>
#include <cstdlib>

#include "caps/field.hpp"
#include "caps/semicrystal.hpp"
#include "caps/tags.hpp"
#include "caps/thermal.hpp"
#include "caps/uff.hpp"

using namespace caps;

namespace {

// Lennard-Jones argon (ε 0.2381 kcal/mol, σ 3.405 Å: the LAMMPS manual's Green–Kubo example) on an fcc lattice of
// n³ cells with lattice constant a, or a liquid start (the same lattice, melted by the caller).
System argon_fcc(int n, double a) {
  System s;
  const double f[4][3] = {{0, 0, 0}, {0.5, 0.5, 0}, {0.5, 0, 0.5}, {0, 0.5, 0.5}};
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      for (int k = 0; k < n; ++k)
        for (auto& b : f) {
          Atom at;
          at.element = 18, at.name = "Ar";
          at.pos = {(i + b[0]) * a, (j + b[1]) * a, (k + b[2]) * a};
          s.atoms.push_back(at);
        }
  s.cell.a = {n * a, 0, 0}, s.cell.b = {0, n * a, 0}, s.cell.c = {0, 0, n * a};
  return s;
}

std::shared_ptr<ForceField> argon_ff(size_t n) {
  auto ff = std::make_shared<ForceField>();
  ff->name = "argon";
  ff->type_names = {"Ar"};
  ff->lj = {{0.2381, 3.405}};
  ff->atom_type.assign(n, "Ar");
  ff->why.assign(n, "");
  ff->type_index.assign(n, 0);
  ff->charge.assign(n, 0.0);
  ff->mass.assign(n, 39.948);
  ff->excluded.assign(n, {});
  return ff;
}

SemicrystallineOptions small_stack() {
  SemicrystallineOptions o;
  o.nx = 3, o.ny = 4, o.nc = 8, o.amorphous_nm = 2.5;
  o.loops = 0.5, o.bridges = 0.2, o.tails = 0.3;
  o.seed = 7;
  o.data_dir = std::string(CAPS_SOURCE_DIR) + "/data";
  return o;
}

}  // namespace

TEST(Semicrystalline, PeCellFromBunn) {
  const PolymerCrystalCell pc = polymer_crystal_cell("pe", std::string(CAPS_SOURCE_DIR) + "/data");
  EXPECT_EQ(pc.id, "pe-orthorhombic");
  EXPECT_NE(pc.source.find("Bunn"), std::string::npos);
  ASSERT_EQ(pc.cell.atoms.size(), 12u);
  EXPECT_NEAR(pc.cell.cell.volume(), 7.40 * 4.93 * 2.534, 1e-9);
  const auto st = find_stems(pc.cell);
  ASSERT_EQ(st.size(), 2u);   // the corner and centre chains
  for (const auto& s : st) {
    EXPECT_EQ(s.formula, "C2H4");
    EXPECT_EQ(s.backbone.size(), 2u);
    EXPECT_FALSE(s.helical);
  }
  // C-C 1.53 Å and the zig-zag angle 112° (Bunn 1939, p. 490)
  const Cell& c = pc.cell.cell;
  const Vec3 c0 = pc.cell.atoms[0].pos, c1 = pc.cell.atoms[1].pos;
  const double d = norm(c.minimum_image(c1 - c0));
  EXPECT_NEAR(d, 1.53, 0.01);
  // density of the crystal: 4 CH2 in a b c
  EXPECT_NEAR(pc.cell.density(), 4 * 14.027 / 6.02214076e23 / (7.40 * 4.93 * 2.534e-24), 1e-3);
}

TEST(Semicrystalline, FindsHelicalStemsInANonPeCell) {
  // a synthetic 3_1 helix of carbons with methyl branches (isotactic-PP-like geometry, not crystal data): one stem,
  // backbone 6 atoms per repeat, helical
  System s;
  const double c = 6.5, r = 0.8;
  s.cell.a = {10, 0, 0}, s.cell.b = {0, 10, 0}, s.cell.c = {0, 0, c};
  std::vector<Vec3> bb;
  for (int k = 0; k < 6; ++k) {
    const double t = 2 * M_PI * k / 3.0;
    bb.push_back({5 + r * std::cos(t), 5 + r * std::sin(t), (k + 0.25) * c / 6});
  }
  for (const auto& p : bb) {
    Atom a;
    a.element = 6, a.pos = p;
    s.atoms.push_back(a);
  }
  for (int k = 1; k < 6; k += 2) {   // a methyl on every second backbone carbon, outward
    Vec3 out = bb[size_t(k)] - Vec3{5, 5, bb[size_t(k)][2]};
    out = out * (1.53 / norm(out));
    Atom a;
    a.element = 6, a.pos = bb[size_t(k)] + out;
    s.atoms.push_back(a);
  }
  const auto st = find_stems(s);
  ASSERT_EQ(st.size(), 1u);
  EXPECT_EQ(st[0].backbone.size(), 6u);
  EXPECT_TRUE(st[0].helical);
  EXPECT_EQ(st[0].atoms.size(), 9u);
}

TEST(Semicrystalline, TopologyIsValid) {
  SemicrystallineReport rep;
  const SemicrystallineOptions o = small_stack();
  const System s = build_semicrystalline(o, &rep);
  // every chain one connected path of carbons with two ends; tails two per chain
  const int nchains = rep.chains;
  EXPECT_EQ(rep.tails, 2 * nchains);
  EXPECT_EQ(rep.stems, 2 * o.nx * o.ny);
  std::map<int64_t, int> ends;
  std::vector<int> cdeg(s.atoms.size(), 0);
  for (const auto& b : s.bonds)
    if (s.atoms[b.i].element == 6 && s.atoms[b.j].element == 6) {
      cdeg[b.i]++, cdeg[b.j]++;
      EXPECT_EQ(s.atoms[b.i].mol, s.atoms[b.j].mol);
    }
  for (size_t i = 0; i < s.atoms.size(); ++i) {
    if (s.atoms[i].element != 6) continue;
    EXPECT_GE(cdeg[i], 1);
    EXPECT_LE(cdeg[i], 2);
    if (cdeg[i] == 1) ends[s.atoms[i].mol]++;
  }
  EXPECT_EQ(int(ends.size()), nchains);
  for (const auto& [m, k] : ends) EXPECT_EQ(k, 2) << "molecule " << m;
  // each carbon has 4 bonds (CH2 inside, CH3 at the ends)
  std::vector<int> deg(s.atoms.size(), 0);
  for (const auto& b : s.bonds) deg[b.i]++, deg[b.j]++;
  for (size_t i = 0; i < s.atoms.size(); ++i) EXPECT_EQ(deg[i], s.atoms[i].element == 6 ? 4 : 1);
  // segment counts within one or two of the targets (fractions per segment), loops + bridges + tails use every end
  const double S = rep.loops + rep.bridges + rep.tails;
  EXPECT_EQ(2 * rep.loops + 2 * rep.bridges + rep.tails, 2 * rep.stems);
  EXPECT_NEAR(rep.bridges / S, 0.2, 2.0 / S);
  EXPECT_NEAR(rep.tails / S, 0.3, 2.5 / S);
  // the order-parameter analysis of the cell finds the builder's topology
  const SemicrystallineReport a = analyse_semicrystalline(s);
  EXPECT_EQ(a.loops, rep.loops);
  EXPECT_EQ(a.bridges, rep.bridges);
  EXPECT_EQ(a.tails, rep.tails);
  EXPECT_EQ(a.stems, rep.stems);
  // groups: crystal and amorphous cover every atom once
  const int ct = tag_index(s, "crystal"), am = tag_index(s, "amorphous");
  ASSERT_GE(ct, 0);
  ASSERT_GE(am, 0);
  size_t nc = 0, na = 0;
  for (const auto& at : s.atoms) {
    const bool c = (at.tags >> ct) & 1u, m = (at.tags >> am) & 1u;
    EXPECT_NE(c, m);
    nc += c, na += m;
  }
  EXPECT_EQ(nc, size_t(3 * 2 * o.nx * o.ny * 2 * o.nc));   // CH2 × stems × carbons per stem
  // the amorphous layer at its density (CH2 count rounded; end hydrogens)
  EXPECT_NEAR(rep.amorphous_density, o.density, 0.02);
  // no clashes: non-bonded carbons (more than 3 bonds apart) at least 2 Å apart as grown
  EXPECT_GT(rep.worst_contact, 2.0 - o.contact);
}

TEST(Semicrystalline, CrystalMatchesTheUnitCell) {
  SemicrystallineReport rep;
  const SemicrystallineOptions o = small_stack();
  const System s = build_semicrystalline(o, &rep);
  EXPECT_NEAR(rep.lx, 7.40 * o.nx, 1e-9);
  EXPECT_NEAR(rep.ly, 4.93 * o.ny, 1e-9);
  EXPECT_NEAR(rep.crystal_thickness, 2.534 * o.nc, 1e-9);
  // crystal carbons sit on the Bunn sites: fractional x, y in the unit cell at ±0.038 / ±0.065 about a chain axis
  const int ct = tag_index(s, "crystal");
  int checked = 0;
  for (const auto& at : s.atoms) {
    if (at.element != 6 || !((at.tags >> ct) & 1u)) continue;
    double fx = at.pos[0] / 7.40, fy = at.pos[1] / 4.93;
    fx -= std::floor(fx), fy -= std::floor(fy);
    // distance to the nearest site of the four in the cell (with images)
    const double sites[4][2] = {{0.038, 0.935}, {0.962, 0.065}, {0.462, 0.435}, {0.538, 0.565}};
    double best = 1;
    for (auto& q : sites)
      for (int dx = -1; dx <= 1; ++dx)
        for (int dy = -1; dy <= 1; ++dy) best = std::min(best, std::hypot(fx - q[0] - dx, fy - q[1] - dy));
    EXPECT_LT(best, 1e-6);
    ++checked;
  }
  EXPECT_EQ(checked, 2 * o.nx * o.ny * 2 * o.nc);
  EXPECT_LT(rep.tilt_deg, 3.0);   // straight stems (the runs' ends reach a little into the interphase)
  EXPECT_GT(rep.crystallinity_p2, rep.crystallinity_mass - 0.02);
}

TEST(Semicrystalline, AtomBudget) {
  SemicrystallineOptions o = small_stack();
  o.nx = 10, o.ny = 10, o.nc = 30, o.amorphous_nm = 6;
  EXPECT_THROW(build_semicrystalline(o), std::invalid_argument);
  o.loops = 0, o.bridges = 1, o.tails = 0;
  o.nx = 2, o.ny = 2, o.nc = 4, o.amorphous_nm = 0.6;   // bridges cannot cross at this density in so thin a layer
  EXPECT_THROW(build_semicrystalline(o), std::invalid_argument);
}

TEST(HeatFlux, PerAtomTermsSumToTheEvaluator) {
  // a small lamellar cell with GAFF (bonds, angles, torsions, LJ with 1-4 pairs, Gasteiger charges with DSF)
  SemicrystallineOptions o = small_stack();
  o.nx = 3, o.ny = 4, o.nc = 6, o.amorphous_nm = 2.5;
  System s = build_semicrystalline(o);
  const ForceField ff = default_forcefield(s);
  EnergyOptions e;
  e.cutoff = 7.0;
  e.tail = false;
  e.threads = 1;
  std::vector<double> x(3 * s.atoms.size());
  for (size_t i = 0; i < s.atoms.size(); ++i)
    for (int k = 0; k < 3; ++k) x[3 * i + k] = s.atoms[i].pos[k];
  Evaluator ev(ff, e);
  std::vector<double> f;
  const EnergyTerms et = ev.compute(x, s.cell, f);
  HeatFlux hf(ff, e);
  std::vector<double> pe, fh;
  std::vector<std::array<double, 9>> w;
  hf.per_atom(x, s.cell, pe, w, &fh);
  const double sum = std::accumulate(pe.begin(), pe.end(), 0.0) + hf.constant_energy();
  EXPECT_NEAR(sum, et.total(), 1e-6 * std::fabs(et.total()) + 1e-6);
  // (the unrelaxed cell has close H-H contacts with huge forces: relative agreement)
  for (size_t k = 0; k < f.size(); ++k) EXPECT_NEAR(fh[k], f[k], 1e-9 * std::fabs(f[k]) + 1e-5) << k;
  // the centroid virials sum to the virial tensor (its symmetric part: CAPS keeps W symmetrised)
  double W[9] = {0};
  for (const auto& a : w)
    for (int k = 0; k < 9; ++k) W[k] += a[size_t(k)];
  const double tol = 1e-6 * std::fabs(et.virial) + 1e-6;
  EXPECT_NEAR(W[0], et.w[0], tol);
  EXPECT_NEAR(W[4], et.w[1], tol);
  EXPECT_NEAR(W[8], et.w[2], tol);
  EXPECT_NEAR(0.5 * (W[1] + W[3]), et.w[3], tol);
  EXPECT_NEAR(0.5 * (W[2] + W[6]), et.w[4], tol);
  EXPECT_NEAR(0.5 * (W[5] + W[7]), et.w[5], tol);
}

TEST(HeatFlux, GreenKuboOfAKnownProcess) {
  // AR(1) samples q_{t+1} = a q_t + noise: autocorrelation σ² a^k, trapezoid integral Δt σ² (1 + a) / (2 (1 − a))
  std::mt19937_64 g(5);
  std::normal_distribution<double> nd(0, 1);
  const double a = 0.9, sig2 = 1 / (1 - a * a);
  const size_t n = 400000;
  std::vector<Vec3> q(n);
  Vec3 cur{0, 0, 0};
  for (size_t t = 0; t < n; ++t) {
    for (int k = 0; k < 3; ++k) cur[k] = a * cur[k] + nd(g);
    q[t] = cur;
  }
  const double dt = 2.0, V = 1000.0, T = 300.0;
  const ThermalResult r = green_kubo_kappa(q, dt, V, T, 0.4, 5);   // 200 lags, the correlation long gone
  const double expect = dt * sig2 * (1 + a) / (2 * (1 - a)) / (V * 0.0019872043 * T * T) * kKappaSI;
  for (int k = 0; k < 3; ++k) EXPECT_NEAR(r.kappa[size_t(k)], expect, 0.05 * expect);
  EXPECT_TRUE(r.plateau);
}

TEST(HeatFlux, NemdOfALinearProfile) {
  // slabs of 2 Å, T rising 1 K per slab from the cold slab to the hot one and back: 0.5 K/Å; 10 kcal/mol over 1 ns
  // through 400 Å² (both directions): κ = (10 / (2·400·1e6)) / 0.5 in (kcal/mol)/(Å fs K)
  std::vector<double> T;
  for (int k = 0; k < 20; ++k) T.push_back(300 + (k <= 10 ? k : 20 - k));
  const ThermalResult r = nemd_kappa(T, 2.0, 10.0, 1e6, 400.0);
  EXPECT_NEAR(r.kappa_mean, 10.0 / (2 * 400 * 1e6) / 0.5 * kKappaSI, 1e-9);
  EXPECT_NEAR(r.gradient, 0.5e10, 1e-3);
}

TEST(HeatFlux, TwoDimensionalSheetUsesAnEffectiveThickness) {
  // an argon monolayer (triangular, 3.8 Å) in xy with 40 Å of vacuum: in-plane κ only, over A·t
  System s;
  const int nx = 8, ny = 8;
  const double d = 3.8;
  for (int i = 0; i < nx; ++i)
    for (int j = 0; j < ny; ++j) {
      Atom at;
      at.element = 18;
      at.pos = {(i + 0.5 * (j % 2)) * d, j * d * std::sqrt(3.0) / 2, 20.0};
      s.atoms.push_back(at);
    }
  s.cell.a = {nx * d, 0, 0}, s.cell.b = {0, ny * d * std::sqrt(3.0) / 2, 0}, s.cell.c = {0, 0, 40};
  EXPECT_NEAR(sheet_thickness(s), 2 * 1.88, 1e-9);   // Bondi argon 1.88 Å on both sides
  ThermalOptions o;
  o.field = argon_ff(s.atoms.size());
  o.energy.cutoff = 8.5;
  o.energy.coulomb = false;
  o.dimensionality = 2;
  o.thickness = 3.4;
  o.temperature = 30;
  o.equilibrate_ps = 2;
  o.ps = 15;
  o.dt = 5;
  o.sample_every = 2;
  o.corr_ps = 2;
  const ThermalResult r = thermal_conductivity(s, o);
  EXPECT_EQ(r.dimensionality, 2);
  EXPECT_EQ(r.kappa[2], 0.0);
  EXPECT_NEAR(r.kappa_mean, 0.5 * (r.kappa[0] + r.kappa[1]), 1e-12);
  EXPECT_NEAR(r.sheet_conductance, r.kappa_mean * 3.4e-10, 1e-20);
  // the vacuum must keep the sheet from its image; a flux across the sheet is refused
  o.method = "nemd";
  o.axis = 2;
  EXPECT_THROW(thermal_conductivity(s, o), std::invalid_argument);
}

TEST(HeatFlux, LennardJonesGreenKuboAgreesWithNemd) {
  // a short check on a small Lennard-Jones liquid (ρ* 0.6, T* 1.35, rc 2.5σ: the LAMMPS examples/KAPPA state, κ* about
  // 3.3–3.4 by NEMD there): Green–Kubo after 1500 steps of NVT (the fcc start melts) with 3000 of NVE, then NEMD from its end, 4500 steps; agreement within the noise of so short
  // runs (a factor of 1.6), and κ* in a sane range
  const double sigma = 3.405, eps = 0.2381, kB = 0.0019872043;
  // CAPS_LJ_SEED picks another start (repeat the check over seeds); CAPS_LJ_STATE=triple runs T* 0.722, ρ* 0.8442 (κ*
  // 6.99 ± 0.03 extrapolated to large N by Viscardy, Servantie & Gaspard, J. Chem. Phys. 126, 184513 (2007); 6.93 at N 500)
  const bool triple = std::getenv("CAPS_LJ_STATE") && std::string(std::getenv("CAPS_LJ_STATE")) == "triple";
  const double rstar = triple ? 0.8442 : 0.6, tstar = triple ? 0.722 : 1.35;
  const double rho = rstar / (sigma * sigma * sigma);
  const int n = 5;   // 500 atoms
  const double a = std::cbrt(4 / rho);
  System s = argon_fcc(n, a);
  const double unit = 1.380649e-23 / (sigma * sigma * 1e-20) * std::sqrt(eps * 4184 / 6.02214076e23 / (39.948e-3 / 6.02214076e23));
  ThermalOptions o;
  o.field = argon_ff(s.atoms.size());
  o.energy.cutoff = 2.5 * sigma;
  o.energy.coulomb = false;
  o.temperature = tstar * eps / kB;
  o.dt = 5;   // 0.0023 τ
  o.equilibrate_ps = 7.5;   // the fcc start melts and settles
  o.ps = 15;
  o.sample_every = 2;
  o.corr_ps = std::getenv("CAPS_LJ_CORR") ? std::stod(std::getenv("CAPS_LJ_CORR")) : 3;
  o.seed = std::getenv("CAPS_LJ_SEED") ? std::stoull(std::getenv("CAPS_LJ_SEED")) : 3;
  System g = s;
  const ThermalResult gk = thermal_conductivity(g, o);
  o.method = "nemd";
  o.axis = 2;
  o.slabs = 10;
  o.swap_every = std::getenv("CAPS_LJ_SWAP") ? std::stoi(std::getenv("CAPS_LJ_SWAP")) : 20;
  o.equilibrate_ps = 0;   // from the liquid the Green–Kubo run ends with
  o.ps = 22.5;
  System m = g;
  const ThermalResult ne = thermal_conductivity(m, o);
  const double kg = gk.kappa_mean / unit, kn = ne.kappa_mean / unit;
  if (std::getenv("CAPS_LJ_SEED")) {
    for (size_t l = 0; l < gk.t_ps.size(); l += gk.t_ps.size() / 12)
      std::printf("  t %.2f ps running κ* %.2f acf %.3f\n", gk.t_ps[l], (gk.running[0][l] + gk.running[1][l] + gk.running[2][l]) / 3 / unit, (gk.acf[0][l] + gk.acf[1][l] + gk.acf[2][l]) / 3);
  }
  std::printf("T %.1f K LJ ρ* %.4g T* %.4g seed %llu: κ* Green–Kubo %.2f ± %.2f, NEMD %.2f (LAMMPS examples/KAPPA: 3.3–3.9)\n", gk.temperature, rstar, tstar, (unsigned long long)o.seed, kg,
              gk.kappa_mean_error / unit, kn);
  const double ref = triple ? 6.93 : 3.4;
  EXPECT_GT(kg, 0.45 * ref);
  EXPECT_LT(kg, 1.8 * ref);
  EXPECT_GT(kn, 0.45 * ref);
  EXPECT_LT(kn, 1.8 * ref);
  if (!std::getenv("CAPS_LJ_SEED")) EXPECT_LT(std::max(kg, kn) / std::min(kg, kn), 1.6);
}
TEST(Semicrystalline, StemsOfTheLibrarysOtherCrystals) {
  // the polymer crystal library's cells (data/crystals/polymers): stems found along the chain axis, helices told apart
  const std::string dir = std::string(CAPS_SOURCE_DIR) + "/data";
  struct Want { const char* id; size_t stems, backbone; bool helical; int axis; };
  const Want want[] = {{"pe-orthorhombic", 2, 2, false, 2}, {"ipp-alpha", 4, 6, true, 2}, {"peo", 4, 21, true, 2}, {"pom-orthorhombic", 2, 4, true, 2},
                       {"pvdf-beta", 2, 2, false, 2}, {"nylon66-alpha", 1, 14, false, 2}, {"natural-rubber", 4, 8, true, 2}};
  for (const auto& w : want) {
    const auto st = find_stems(polymer_crystal_cell(w.id, dir).cell);
    ASSERT_EQ(st.size(), w.stems) << w.id;
    for (const auto& s : st) {
      EXPECT_EQ(s.backbone.size(), w.backbone) << w.id;
      EXPECT_EQ(s.helical, w.helical) << w.id;
      EXPECT_EQ(s.lattice_axis, w.axis) << w.id;
    }
  }
  // the builder grows polymethylene amorphous layers only so far: a helical iPP stack is refused with the reason
  SemicrystallineOptions o = small_stack();
  o.polymer = "ipp-alpha";
  try {
    build_semicrystalline(o);
    ADD_FAILURE() << "iPP built";
  } catch (const std::invalid_argument& e) {
    EXPECT_NE(std::string(e.what()).find("helical"), std::string::npos) << e.what();
  }
}
