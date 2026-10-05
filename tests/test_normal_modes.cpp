#include <gtest/gtest.h>

#include <cmath>
#include <random>

#include "caps/normal_modes.hpp"
#include "caps/uff.hpp"

using namespace caps;

// The eigen-solver: A v = λ v for every pair, eigenvectors orthonormal, eigenvalues ascending.
TEST(NormalModes, SymmetricEigen) {
  const size_t n = 37;
  std::mt19937 rng(5);
  std::uniform_real_distribution<double> u(-1, 1);
  std::vector<double> A(n * n);
  for (size_t i = 0; i < n; ++i)
    for (size_t j = i; j < n; ++j) A[i * n + j] = A[j * n + i] = u(rng);
  auto V = A;
  std::vector<double> w;
  symmetric_eigen(V, n, w);
  for (size_t k = 0; k + 1 < n; ++k) EXPECT_LE(w[k], w[k + 1]);
  for (size_t k = 0; k < n; ++k) {
    for (size_t i = 0; i < n; ++i) {
      double av = 0;
      for (size_t j = 0; j < n; ++j) av += A[i * n + j] * V[j * n + k];
      EXPECT_NEAR(av, w[k] * V[i * n + k], 1e-10);
    }
    for (size_t l = 0; l < n; ++l) {
      double d = 0;
      for (size_t i = 0; i < n; ++i) d += V[i * n + k] * V[i * n + l];
      EXPECT_NEAR(d, k == l ? 1.0 : 0.0, 1e-10);
    }
  }
}

// A harmonic diatomic, E = k (r − r0)²: one vibration at ν̃ = (1/2πc) √(2k/μ), five rigid motions projected out,
// and its zero-point energy hcν̃/2.
TEST(NormalModes, DiatomicMatchesTheSpring) {
  System s;
  Atom a;
  a.element = 7;
  a.pos = {0, 0, 0};
  s.atoms.push_back(a);
  a.pos = {1.1, 0, 0};
  s.atoms.push_back(a);
  s.bonds = {{0, 1}};
  auto ff = std::make_shared<ForceField>(default_forcefield(s));
  ff->bonds = {{0, 1, 500.0, 1.1}};
  ff->charge = {0, 0};
  NormalModesOptions o;
  o.field = ff;
  const auto r = normal_modes(s, o);
  ASSERT_EQ(r.wavenumber.size(), 1u);
  EXPECT_EQ(r.projected, 5);
  const double mu = ff->mass[0] * ff->mass[1] / (ff->mass[0] + ff->mass[1]);
  const double expect = std::sqrt(2 * 500.0 / mu * 4.184e-4) / (2 * M_PI * 2.99792458e-5);
  EXPECT_NEAR(r.wavenumber[0], expect, 1e-3 * expect);
  EXPECT_NEAR(r.reduced_mass[0], ff->mass[0], 1e-6 * mu);   // 1/Σ|q_i/√m_i|² of the unit mass-weighted mode
  EXPECT_NEAR(r.zpe, 0.5 * expect * 2.859144e-3, 1e-3 * r.zpe);
  // the animation: the largest displacement is the amplitude, a quarter period in
  const auto fr = mode_frames(s, r, 0, 0.2, 8);
  EXPECT_NEAR(std::fabs(fr[2][0][0] - 0.0), 0.2, 1e-9);
}

// Water by the default force field: three vibrations, all real, the bend lowest; six rigid motions projected.
TEST(NormalModes, WaterHasThreeRealModes) {
  System s;
  Atom o, h1, h2;
  o.element = 8, h1.element = 1, h2.element = 1;
  o.pos = {0, 0, 0}, h1.pos = {0.9572, 0, 0}, h2.pos = {-0.2399872, 0.9266272, 0};
  s.atoms = {o, h1, h2};
  s.bonds = {{0, 1}, {0, 2}};
  NormalModesOptions op;
  const auto r = normal_modes(s, op);
  EXPECT_EQ(r.projected, 6);
  ASSERT_EQ(r.wavenumber.size(), 3u);
  // at the force field's own minimum only when the geometry is relaxed: here the check is that the bend is the lowest
  // and the two stretches lie far above it
  EXPECT_LT(r.wavenumber[0], 0.6 * r.wavenumber[1]);
}
