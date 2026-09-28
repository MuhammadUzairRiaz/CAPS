#include <gtest/gtest.h>

#include <cmath>
#include <random>
#include <string>
#include <vector>

#include "caps/elements.hpp"
#include "caps/properties.hpp"
#include "caps/uff.hpp"

using namespace caps;

namespace {

constexpr double kPi = 3.14159265358979323846;

Cell cube(double l) {
  Cell c;
  c.a = {l, 0, 0};
  c.b = {0, l, 0};
  c.c = {0, 0, l};
  return c;
}

// A one-frame trajectory of carbon atoms at the given positions.
Trajectory carbons(const std::vector<Vec3>& pos, double box) {
  Trajectory t;
  for (const auto& p : pos) {
    Atom a;
    a.element = 6;
    a.pos = p;
    t.topology.atoms.push_back(a);
  }
  t.topology.cell = cube(box);
  t.positions.push_back(pos);
  t.cells.push_back(t.topology.cell);
  t.timesteps.push_back(0);
  return t;
}

Trajectory simple_cubic(double a, int n) {
  std::vector<Vec3> pos;
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      for (int k = 0; k < n; ++k) pos.push_back({i * a, j * a, k * a});
  return carbons(pos, n * a);
}

const Property& get(const std::vector<Property>& ps, const std::string& id) {
  for (const auto& p : ps)
    if (p.id == id) return p;
  throw std::runtime_error("no " + id);
}

}  // namespace

TEST(Properties, ScatteringTables) {
  EXPECT_NEAR(xray_f(6, 0.0), 6.0, 0.01);     // f(0) = Z
  EXPECT_NEAR(xray_f(8, 0.0), 8.0, 0.01);
  EXPECT_LT(xray_f(6, 4.0), 3.0);             // form factors fall with q
  EXPECT_DOUBLE_EQ(neutron_b(1), -3.739);
  EXPECT_THROW(neutron_b(92), std::invalid_argument);
}

TEST(Properties, DensityOfKnownCell) {
  const Trajectory t = simple_cubic(3.0, 4);   // 64 C in 12 Å
  AnalyzeOptions o;
  const auto ps = analyze(t, {"density"}, o);
  EXPECT_NEAR(get(ps, "density").value, 64 * element(6).mass / 6.02214076e23 / (12.0 * 12.0 * 12.0 * 1e-24), 1e-4);
}

// Simple cubic crystal: Bragg peaks at the reciprocal lattice of the crystal only; the first at 2π/a.
TEST(Properties, StructureFactorBraggPeak) {
  const double a = 3.0;
  const Trajectory t = simple_cubic(a, 8);
  AnalyzeOptions o;
  o.qmax = 2.5;   // below the second shell (2π√2/a), which is the higher one
  const auto ps_sq = analyze(t, {"sq"}, o);
  const auto& sq = get(ps_sq, "sq");
  EXPECT_NEAR(sq.value, 2 * kPi / a, 0.03);
  // between Bragg shells, S = 0 exactly (every other lattice vector of the box cancels)
  const auto& s = sq.series[0];
  for (size_t k = 0; k < s.x.size(); ++k)
    if (s.x[k] < 1.9) EXPECT_NEAR(s.y[k], 0.0, 1e-9) << "q = " << s.x[k];
}

// Ideal gas: S(q) = 1 on average at every q.
TEST(Properties, StructureFactorIdealGas) {
  std::mt19937 rng(3);
  std::uniform_real_distribution<double> u(0, 20);
  std::vector<Vec3> pos(800);
  for (auto& p : pos) p = {u(rng), u(rng), u(rng)};
  const Trajectory t = carbons(pos, 20);
  AnalyzeOptions o;
  o.qmax = 3.5;
  o.q_direct = 3.5;
  const auto ps = analyze(t, {"sq"}, o);
  const auto& s = get(ps, "sq").series[0];
  double mean = 0;
  size_t n = 0;
  for (size_t k = 0; k < s.x.size(); ++k)
    if (s.x[k] > 1.0) { mean += s.y[k]; ++n; }
  ASSERT_GT(n, 20u);
  EXPECT_NEAR(mean / n, 1.0, 0.08);
}

// One isolated atom: the occupied fraction is its sphere's volume over the cell's.
TEST(Properties, FreeVolumeOfOneSphere) {
  const Trajectory t = carbons({{7.3, 6.9, 7.1}}, 14.0);
  AnalyzeOptions o;
  o.grid = 0.2;
  const auto ps_fv = analyze(t, {"ffv"}, o);
  const auto& fv = get(ps_fv, "ffv");
  const double r = element(6).vdw;
  EXPECT_NEAR(fv.extra.at("van der Waals occupied fraction"), 4.0 / 3.0 * kPi * r * r * r / (14.0 * 14.0 * 14.0), 2e-4);
}

// Simple cubic lattice of spacing a: the largest cavity is the cube centre, diameter 2(a√3/2 − r).
TEST(Properties, LargestPoreOfLattice) {
  const double a = 6.0;
  const Trajectory t = simple_cubic(a, 4);
  AnalyzeOptions o;
  o.grid = 0.3;
  const auto ps_psd = analyze(t, {"psd"}, o);
  const auto& psd = get(ps_psd, "psd");
  EXPECT_NEAR(psd.extra.at("largest diameter (Å)"), 2 * (a * std::sqrt(3.0) / 2 - element(6).vdw), 1e-3);
}

// Independent atoms on Gaussian random walks plus a common drift: the drift is removed and D recovered (Einstein).
TEST(Properties, DiffusionOfRandomWalk) {
  const double D = 0.5;           // Å²/ps
  const double dt = 0.5;          // ps between frames
  const int n = 400, frames = 200;
  std::mt19937 rng(11);
  std::normal_distribution<double> step(0.0, std::sqrt(2 * D * dt));
  std::vector<Vec3> pos(n);
  std::uniform_real_distribution<double> u(0, 30);
  for (auto& p : pos) p = {u(rng), u(rng), u(rng)};
  Trajectory t = carbons(pos, 30);
  t.topology.unwrapped = true;
  for (int f = 1; f < frames; ++f) {
    for (auto& p : pos) p = p + Vec3{step(rng), step(rng), step(rng)};
    std::vector<Vec3> shifted = pos;
    for (auto& p : shifted) p = p + Vec3{0.3 * f, 0, 0};   // drift of the whole system
    t.positions.push_back(shifted);
    t.cells.push_back(t.topology.cell);
    t.timesteps.push_back(f);
  }
  AnalyzeOptions o;
  o.frame_ps = dt;
  const auto ps = analyze(t, {"msd", "diffusion"}, o);
  const auto& d = get(ps, "diffusion");
  EXPECT_NEAR(d.value, D * 10.0, 0.1 * D * 10.0);     // 1 Å²/ps = 10 × 10⁻⁵ cm²/s
  EXPECT_NEAR(d.extra.at("log-log slope β"), 1.0, 0.1);
  EXPECT_TRUE(d.notes.empty());
}

TEST(Properties, UnknownIdAndJson) {
  const Trajectory t = simple_cubic(3.0, 3);
  AnalyzeOptions o;
  EXPECT_THROW(analyze(t, {"nope"}, o), std::invalid_argument);
  const std::string j = properties_json(analyze(t, {"density", "msd"}, o));
  EXPECT_NE(j.find("\"id\":\"density\""), std::string::npos);
  EXPECT_NE(j.find("null"), std::string::npos);   // MSD of one frame: no value
}

TEST(Properties, OrientationOfAlignedAndRandomChains) {
  // nine all-trans zigzag chains along z (a crude PE crystal): S = 1, Herman f along z = 1
  System s;
  s.cell.a = {15, 0, 0};
  s.cell.b = {0, 15, 0};
  s.cell.c = {0, 0, 25.4};
  s.unwrapped = true;
  s.has_mol = true;
  for (int cx = 0; cx < 3; ++cx)
    for (int cy = 0; cy < 3; ++cy) {
      const uint32_t first = uint32_t(s.atoms.size());
      for (int k = 0; k < 20; ++k) {
        Atom a;
        a.element = 6;
        a.mol = cx * 3 + cy + 1;
        a.pos = {2.5 + cx * 5 + (k % 2 ? 0.42 : -0.42), 2.5 + cy * 5, 0.5 + k * 1.27};
        s.atoms.push_back(a);
        if (k > 0) s.bonds.push_back({first + uint32_t(k - 1), first + uint32_t(k), 1});
      }
    }
  Trajectory t;
  t.topology = s;
  std::vector<Vec3> p;
  for (const auto& a : s.atoms) p.push_back(a.pos);
  t.positions.push_back(p);
  t.cells.push_back(s.cell);
  t.timesteps.push_back(0);
  const auto props = analyze(t, {"orientation"}, AnalyzeOptions{});
  ASSERT_EQ(props.size(), 1u);
  EXPECT_NEAR(props[0].value, 1.0, 1e-9);
  EXPECT_NEAR(props[0].extra.at("Herman f along z"), 1.0, 1e-9);
  EXPECT_NEAR(std::fabs(props[0].extra.at("director z")), 1.0, 1e-9);
  EXPECT_GT(props[0].extra.at("local crystallinity (fraction)"), 0.5);
}

// X-ray form factors of every element to Cf: f(0) = Z (the neutral atom's electrons). Electron scattering factors (Peng
// 1996) against the X-ray ones through the Mott–Bethe relation f_e(s) = 0.023934 (Z − f_x(s)) / s² (Å, s = sin θ/λ in Å⁻¹),
// two independent tables agreeing at 0.3 ≤ s ≤ 1 Å⁻¹ within a few per cent.
TEST(Properties, FormFactorsXrayAndElectron) {
  for (int z = 1; z <= 98; ++z)
    if (z != 69) EXPECT_NEAR(xray_f(z, 0), double(z), 0.02 * z + 0.01) << z;
  EXPECT_THROW(xray_f(69, 0), std::invalid_argument);   // thulium's row in the source sums to 63.85: left out
  EXPECT_NEAR(electron_f(6, 0), 2.5088, 1e-4);   // Σ a_i for carbon
  for (int z : {6, 8, 14, 26, 29, 30, 79})
    for (double s : {0.3, 0.5, 1.0}) {
      const double q = 4 * M_PI * s, mb = 0.023934 * (z - xray_f(z, q)) / (s * s);
      EXPECT_NEAR(electron_f(z, q), mb, 0.05 * mb) << z << " s " << s;
    }
}

// Free-volume radii: Bondi's, UFF's x/2, or half the assigned force field's Lennard-Jones minimum; a larger radius leaves
// less free volume.
TEST(Properties, FreeVolumeRadii) {
  System s;
  s.cell.a = {10, 0, 0}, s.cell.b = {0, 10, 0}, s.cell.c = {0, 0, 10};
  Atom c;
  c.element = 6, c.pos = {5, 5, 5};
  s.atoms.push_back(c);
  const auto bondi = free_volume_radii(s, "bondi"), uff = free_volume_radii(s, "uff");
  double x = 0, d = 0;
  ASSERT_TRUE(uff_vdw(6, x, d));
  EXPECT_NEAR(bondi[0], 1.70, 1e-9);
  EXPECT_NEAR(uff[0], 0.5 * x, 1e-12);
  ForceField ff;
  ff.type_index = {0};
  ff.lj = {{0.1, 3.4}};
  EXPECT_NEAR(free_volume_radii(s, "forcefield", &ff)[0], 0.5 * 3.4 * std::pow(2.0, 1.0 / 6), 1e-12);
  ff.pair_form = "lj9-6";
  EXPECT_NEAR(free_volume_radii(s, "forcefield", &ff)[0], 1.7, 1e-12);
  EXPECT_THROW(free_volume_radii(s, "forcefield"), std::invalid_argument);
  EXPECT_THROW(free_volume_radii(s, "connolly"), std::invalid_argument);
  const Trajectory t = carbons({{5, 5, 5}}, 10);
  AnalyzeOptions o;
  const double fb = analyze(t, {"ffv"}, o)[0].extra.at("van der Waals occupied fraction");
  o.radii = "uff";
  const double fu = analyze(t, {"ffv"}, o)[0].extra.at("van der Waals occupied fraction");
  EXPECT_NEAR(fb / fu, std::pow(1.70 / (0.5 * x), 3), 0.08 * fb / fu);   // sphere volumes on a 0.4 Å grid
}
