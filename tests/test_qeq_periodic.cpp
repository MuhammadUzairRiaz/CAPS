#include <gtest/gtest.h>

#include <cmath>

#include "caps/crystal.hpp"
#include "caps/qeq.hpp"

using namespace caps;

namespace {
constexpr double kCoul = 14.399645, kBohr = 0.529177210903, kHartree = 27.211386245988;

System water() {
  System s;
  Atom o, h1, h2;
  o.element = 8, o.pos = {0, 0, 0};
  h1.element = 1, h1.pos = {0.9572, 0, 0};
  h2.element = 1, h2.pos = {0.9572 * std::cos(104.52 * M_PI / 180), 0.9572 * std::sin(104.52 * M_PI / 180), 0};
  s.atoms = {o, h1, h2};
  return s;
}

System quartz333() {
  const System x = read_cif(std::string(CAPS_SOURCE_DIR) + "/data/crystals/alpha-quartz.cif");
  System s;
  s.cell = x.cell;
  s.cell.a = x.cell.a * 3.0, s.cell.b = x.cell.b * 3.0, s.cell.c = x.cell.c * 3.0;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      for (int k = 0; k < 3; ++k)
        for (auto a : x.atoms) {
          a.pos = a.pos + x.cell.a * double(i) + x.cell.b * double(j) + x.cell.c * double(k);
          s.atoms.push_back(a);
        }
  return s;
}
}  // namespace

// The 1s–1s Coulomb integral against its closed form (Roothaan): J(R) = [1 − (1 + 11x/8 + 3x²/4 + x³/6) e^(−2x)] / R,
// x = ζR, and J(0) = 5ζ/8; unequal shells the same either way round; 1/R once the densities part
TEST(QEqPeriodic, SlaterIntegrals) {
  const double z = 1.3;
  EXPECT_NEAR(slater_coulomb_integral(1, z, 1, z, 0.0), 5 * z / 8, 1e-12);
  for (double R : {0.3, 1.0, 2.2, 4.5}) {
    const double x = z * R;
    const double exact = (1 - (1 + 11 * x / 8 + 0.75 * x * x + x * x * x / 6) * std::exp(-2 * x)) / R;
    EXPECT_NEAR(slater_coulomb_integral(1, z, 1, z, R), exact, 1e-11) << R;
  }
  EXPECT_NEAR(slater_coulomb_integral(2, 0.9, 3, 1.4, 1.7), slater_coulomb_integral(3, 1.4, 2, 0.9, 1.7), 1e-11);
  EXPECT_NEAR(slater_coulomb_integral(4, 1.1, 2, 1.6, 0.0), slater_coulomb_integral(2, 1.6, 4, 1.1, 0.0), 1e-11);
  EXPECT_NEAR(slater_coulomb_integral(2, 1.6, 3, 1.2, 40.0), 1 / 40.0, 1e-13);
}

// The Ewald sum is exact: rock salt's Madelung energy, −M e²/(4πε0 r0) per ion pair with M = 1.747564594633
TEST(QEqPeriodic, MadelungOfRockSalt) {
  System s;
  const double a = 5.64;
  s.cell.a = {2 * a, 0, 0}, s.cell.b = {0, 2 * a, 0}, s.cell.c = {0, 0, 2 * a};
  std::vector<double> q;
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j)
      for (int k = 0; k < 4; ++k) {
        Atom x;
        const bool na = (i + j + k) % 2 == 0;
        x.element = na ? 11 : 17, x.pos = {i * a / 2, j * a / 2, k * a / 2};
        s.atoms.push_back(x);
        q.push_back(na ? 1.0 : -1.0);
      }
  QEqOptions o;
  o.kernel = "point", o.ewald = true, o.ewald_accuracy = 1e-12, o.cutoff = a;
  const double E = qeq_coulomb_energy(s, q, o);
  EXPECT_NEAR(E / 32, -1.747564594633 * kCoul / (a / 2), 1e-8);
}

// A water molecule against the equalisation written out by hand with the same Slater integrals: three equations, the
// total charge, and hydrogen's charge-dependent hardness iterated to self-consistency
TEST(QEqPeriodic, WaterByHand) {
  const System s = water();
  double chi[3], J[3], R[3], z[3];
  int n[3];
  for (int i = 0; i < 3; ++i) {
    ASSERT_TRUE(qeq_parameters(s.atoms[size_t(i)].element, chi[i], J[i], R[i]));
    n[i] = principal_quantum_number(s.atoms[size_t(i)].element);
    z[i] = 0.4913 * (2 * n[i] + 1) / (2 * R[i] / kBohr);
  }
  double K[3][3];
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      if (i != j) K[i][j] = kHartree * slater_coulomb_integral(n[i], z[i], n[j], z[j], norm(s.atoms[size_t(j)].pos - s.atoms[size_t(i)].pos) / kBohr);
  // symmetric water: q_O = −2 q_H; χ_O + J_O q_O + 2 K_OH q_H = χ_H + J_H(q_H) q_H + K_OH q_O + K_HH q_H
  auto solve = [&](bool hyd) {
    double qh = 0;
    for (int it = 0; it < 2000; ++it) {
      const double JH = J[1] * (1 + (hyd ? qh / z[1] : 0));
      // χ_O − 2 J_O qh + 2 K_OH qh = χ_H + JH qh − 2 K_OH qh + K_HH qh
      const double next = (chi[0] - chi[1]) / (JH + 2 * J[0] - 4 * K[0][1] + K[1][2]);
      if (std::fabs(next - qh) < 1e-14) return next;
      qh += 0.5 * (next - qh);
    }
    return qh;
  };
  QEqOptions o;
  o.kernel = "slater", o.cutoff = 1000, o.tolerance = 1e-12, o.scf_tolerance = 1e-11, o.scf_iterations = 500;
  const auto q0 = qeq_charges(s, o);
  const double h0 = solve(false);
  EXPECT_NEAR(q0[1], h0, 2e-6);
  EXPECT_NEAR(q0[2], h0, 2e-6);
  EXPECT_NEAR(q0[0], -2 * h0, 4e-6);
  o.hydrogen = true;
  QEqReport rep;
  const auto q1 = qeq_charges(s, o, &rep);
  const double h1 = solve(true);
  EXPECT_NEAR(q1[1], h1, 2e-6);
  EXPECT_NEAR(q1[0], -2 * h1, 4e-6);
  EXPECT_GT(std::fabs(h1 - h0), 1e-3);   // the correction matters
  EXPECT_GT(rep.scf, 1);
}

// Periodic quartz: Ewald charges are the same at two real-space cut-offs (the Ewald split moves; silicon's diffuse 3s
// density differs from a point by 1·10⁻⁷ of 1/r at 12 Å), every Si alike and every
// O alike, neutral; a molecule alone in a large box is the molecule without the box
TEST(QEqPeriodic, QuartzAndTheBox) {
  const System s = quartz333();
  QEqOptions o;
  o.kernel = "slater", o.ewald = true, o.ewald_accuracy = 1e-9, o.tolerance = 1e-11;
  o.cutoff = 14.0;
  const auto a = qeq_charges(s, o);
  o.cutoff = 12.0;
  const auto b = qeq_charges(s, o);
  double sum = 0, si = 0, ox = 0;
  for (size_t i = 0; i < a.size(); ++i) {
    EXPECT_NEAR(a[i], b[i], 1e-5);
    sum += a[i];
    (s.atoms[i].element == 14 ? si : ox) = a[i];
  }
  EXPECT_NEAR(sum, 0, 1e-8);
  for (size_t i = 0; i < a.size(); ++i) EXPECT_NEAR(a[i], s.atoms[i].element == 14 ? si : ox, 1e-6);
  EXPECT_GT(si, 0);
  EXPECT_NEAR(si, -2 * ox, 1e-6);

  System w = water();
  QEqOptions m;
  m.kernel = "slater", m.cutoff = 1000, m.tolerance = 1e-12;
  const auto free = qeq_charges(w, m);
  for (auto& at : w.atoms) at.pos = at.pos + Vec3{20, 20, 20};
  w.cell.a = {40, 0, 0}, w.cell.b = {0, 40, 0}, w.cell.c = {0, 0, 40};
  m.ewald = true, m.cutoff = 12, m.ewald_accuracy = 1e-10;
  const auto boxed = qeq_charges(w, m);
  for (size_t i = 0; i < 3; ++i) EXPECT_NEAR(boxed[i], free[i], 2e-4);
}

// Bounds: oxygen held at −0.5 e (free: −0.63 e) leaves silicon +1.0 e (neutral), every atom exactly there
TEST(QEqPeriodic, ChargeBounds) {
  const System s = quartz333();
  QEqOptions o;
  o.kernel = "slater", o.ewald = true, o.cutoff = 10, o.tolerance = 1e-11;
  const auto free = qeq_charges(s, o);
  for (size_t i = 0; i < s.atoms.size(); ++i)
    if (s.atoms[i].element == 8) ASSERT_LT(free[i], -0.5);
  o.bounds[8] = {-0.5, 0.0};
  QEqReport rep;
  const auto q = qeq_charges(s, o, &rep);
  for (size_t i = 0; i < s.atoms.size(); ++i) EXPECT_NEAR(q[i], s.atoms[i].element == 8 ? -0.5 : 1.0, 1e-9);
  EXPECT_EQ(rep.bounded, 162);
}
