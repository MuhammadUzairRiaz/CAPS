#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstdio>

#include "caps/normal_modes.hpp"
#include "caps/phonons.hpp"
#include "caps/uff.hpp"

using namespace caps;

namespace {
constexpr double kPi = 3.14159265358979323846;
constexpr double kEps = 0.2385, kSig = 3.405, kA = 5.26, kRc = 8.5;
// primitive fcc argon, one atom
System argon() {
  System s;
  s.cell.a = {0, kA / 2, kA / 2}, s.cell.b = {kA / 2, 0, kA / 2}, s.cell.c = {kA / 2, kA / 2, 0};
  Atom at;
  at.element = 18;
  at.pos = {0, 0, 0};
  s.atoms.push_back(at);
  return s;
}
ForceField lj(const System& s) {
  ForceField ff = default_forcefield(s);
  for (auto& t : ff.lj) t = {kEps, kSig};
  ff.charge.assign(s.atoms.size(), 0.0);
  return ff;
}
// the analytic dynamical matrix of a monatomic Bravais lattice with pair potential φ (cut at rc): 3 eigen-wavenumbers
std::vector<double> analytic(const Vec3& k, double mass) {
  std::vector<double> D(9, 0.0);
  const System s = argon();
  for (int i = -6; i <= 6; ++i)
    for (int j = -6; j <= 6; ++j)
      for (int l = -6; l <= 6; ++l) {
        if (!i && !j && !l) continue;
        const Vec3 R = s.cell.a * double(i) + s.cell.b * double(j) + s.cell.c * double(l);
        const double r = norm(R);
        if (r >= kRc) continue;
        const double sr6 = std::pow(kSig / r, 6);
        const double d1 = 4 * kEps * (-12 * sr6 * sr6 + 6 * sr6) / r, d2 = 4 * kEps * (156 * sr6 * sr6 - 42 * sr6) / (r * r);
        const double c = 1 - std::cos(dot(k, R));
        for (int a = 0; a < 3; ++a)
          for (int b = 0; b < 3; ++b) {
            const double ra = R[size_t(a)] / r, rb = R[size_t(b)] / r;
            D[size_t(3 * a + b)] += (d2 * ra * rb + d1 / r * ((a == b) - ra * rb)) * c / mass;
          }
      }
  std::vector<double> w;
  symmetric_eigen(D, 3, w);
  for (auto& x : w) x = (x < 0 ? -1 : 1) * 108.5913 * std::sqrt(std::fabs(x));
  return w;
}
}  // namespace

// fcc argon with a cut Lennard-Jones pair: every band on the path equals the analytic lattice sum, the DOS holds three
// modes per cell, C_v reaches 3k at high T, and the sound speeds are the analytic slopes
TEST(Phonons, LennardJonesArgonIsExact) {
  const System s = argon();
  PhononOptions o;
  o.field_for = lj;
  o.energy.cutoff = kRc;
  o.energy.coulomb = false;
  o.energy.tail = false;
  o.path_points = 60;
  o.mesh = {10, 10, 10};
  o.temperatures = {50, 3000};
  const auto r = phonons(s, o);
  const double m = lj(s).mass[0];
  ASSERT_FALSE(r.bands.empty());
  // the path: Γ X S Y Γ Z R of the primitive cell's reciprocal axes; check every tenth point
  const double vol = s.cell.volume();
  const Vec3 B0 = cross(s.cell.b, s.cell.c) * (2 * kPi / vol), B1 = cross(s.cell.c, s.cell.a) * (2 * kPi / vol), B2 = cross(s.cell.a, s.cell.b) * (2 * kPi / vol);
  const std::vector<std::array<double, 3>> pts = {{0, 0, 0}, {0.5, 0, 0}, {0.5, 0.5, 0}, {0, 0.5, 0}, {0, 0, 0}, {0, 0, 0.5}, {0.5, 0.5, 0.5}};
  double worst = 0;
  for (size_t i = 0; i < r.path_x.size(); i += 7) {
    // locate the point on the path by its x
    double xs = 0;
    Vec3 k{0, 0, 0};
    for (size_t sgi = 0; sgi + 1 < pts.size(); ++sgi) {
      const Vec3 ka = B0 * pts[sgi][0] + B1 * pts[sgi][1] + B2 * pts[sgi][2], kb = B0 * pts[sgi + 1][0] + B1 * pts[sgi + 1][1] + B2 * pts[sgi + 1][2];
      const double L = norm(kb - ka);
      if (r.path_x[i] <= xs + L + 1e-12) { k = ka + (kb - ka) * ((r.path_x[i] - xs) / L); break; }
      xs += L;
    }
    const auto w = analytic(k, m);
    for (int b = 0; b < 3; ++b) worst = std::max(worst, std::fabs(r.bands[i][size_t(b)] - w[size_t(b)]) / std::max(5.0, w[size_t(b)]));
  }
  std::printf("argon: supercell %d×%d×%d, max |ν − ν_analytic| / ν %.2e, top of the bands %.2f cm⁻¹\n", r.supercell[0], r.supercell[1], r.supercell[2], worst,
              analytic(B0 * 0.5 + B1 * 0.5 + B2 * 0.5, m)[2]);
  EXPECT_LT(worst, 2e-4);
  // DOS: three modes per cell
  double area = 0;
  for (size_t i = 0; i < r.dos.size(); ++i) area += r.dos[i] * (i + 1 < r.dos_w.size() ? r.dos_w[i + 1] - r.dos_w[i] : 0.0);
  EXPECT_NEAR(area, 3.0, 0.01);
  // C_v → 3k per cell (cal/mol/K) at 3000 K; the 50 K value lies below it
  ASSERT_EQ(r.Cv.size(), 2u);
  EXPECT_NEAR(r.Cv[1], 3 * 1.9872043 * (1 - 1.0 / 1000), 0.01 * 3 * 1.9872043);
  EXPECT_LT(r.Cv[0], r.Cv[1]);
  // sound speeds along a*: the analytic slope at the same small k
  const Vec3 k = B0 * 0.01;
  const auto w = analytic(k, m);
  for (int b = 0; b < 3; ++b) {
    const double v = 2 * kPi * 2.99792458e10 * w[size_t(b)] / (norm(k) * 1e10);
    EXPECT_NEAR(r.sound[0][size_t(b)], v, 1e-3 * v);
  }
  std::printf("argon sound speeds along a*: %.0f %.0f %.0f m/s\n", r.sound[0][0], r.sound[0][1], r.sound[0][2]);
}

// The relaxation in the supercell's copies: the conventional fcc cell with one atom pushed 0.05 Å returns to the
// lattice, and its Γ-point modes are then the primitive ones at Γ and at the three X points folded onto it
TEST(Phonons, RelaxesInTheCopiesAndFoldsX) {
  System s;
  s.cell.a = {kA, 0, 0}, s.cell.b = {0, kA, 0}, s.cell.c = {0, 0, kA};
  for (const Vec3 f : std::vector<Vec3>{{0, 0, 0}, {0.5, 0.5, 0}, {0.5, 0, 0.5}, {0, 0.5, 0.5}}) {
    Atom at;
    at.element = 18;
    at.pos = f * kA;
    s.atoms.push_back(at);
  }
  s.atoms[1].pos = s.atoms[1].pos + Vec3{0.05, 0, 0};
  PhononOptions o;
  o.field_for = lj;
  o.energy.cutoff = kRc;
  o.energy.coulomb = false;
  o.energy.tail = false;
  o.path = {{0, 0, 0}, {0.5, 0, 0}};
  o.path_points = 4;
  o.mesh = {2, 2, 2};
  o.relax_ftol = 1e-5;
  const auto r = phonons(s, o);
  EXPECT_LT(r.max_force, 1e-5);
  EXPECT_NEAR(r.relax_shift, 0.05 * 0.75, 0.01);   // the pushed atom comes back, its neighbours give way a little
  // Γ of the conventional cell
  const double m = lj(s).mass[0];
  std::vector<double> want;
  for (const Vec3 k : std::vector<Vec3>{{0, 0, 0}, {2 * kPi / kA, 0, 0}, {0, 2 * kPi / kA, 0}, {0, 0, 2 * kPi / kA}})
    for (double w : analytic(k, m)) want.push_back(w);
  std::sort(want.begin(), want.end());
  ASSERT_EQ(r.bands.front().size(), 12u);
  for (size_t i = 0; i < 12; ++i) EXPECT_NEAR(r.bands.front()[i], want[i], 2e-3 * std::max(5.0, want[i])) << i;
  std::printf("conventional fcc Γ: %.2f … %.2f cm⁻¹ (X points of the primitive cell %.2f … %.2f) · relaxed %d steps, %.3f Å\n", r.bands.front()[3],
              r.bands.front()[11], want[3], want[11], r.relax_steps, r.relax_shift);
}
