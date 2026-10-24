#include <gtest/gtest.h>

#include <cmath>
#include <memory>
#include <random>
#include <vector>

#include "caps/dynamics.hpp"
#include "caps/field.hpp"

using namespace caps;

namespace {

Cell box(double l) {
  Cell c;
  c.a = {l, 0, 0};
  c.b = {0, l, 0};
  c.c = {0, 0, l};
  return c;
}

// Reactive sites only: site classes as given, no other interaction (ε = 0 for the types), one pair type 1–2.
ForceField sites(const std::vector<int>& cls, const SwapModel::Pair& p) {
  ForceField ff;
  const size_t n = cls.size();
  ff.name = "swap test";
  ff.atom_type.assign(n, "S");
  ff.type_index.assign(n, 0);
  ff.type_names = {"S"};
  ff.lj = {{0.0, 1.0}};
  ff.charge.assign(n, 0.0);
  ff.mass.assign(n, 12.0);
  ff.excluded.assign(n, {});
  ff.swap3.site = cls;
  ff.swap3.pairs = {p};
  return ff;
}

EnergyTerms eval(const ForceField& ff, const std::vector<double>& x, const Cell& c, std::vector<double>& f) {
  EnergyOptions eo;
  eo.coulomb = false;
  eo.tail = false;
  eo.cutoff = 3.0;
  Evaluator ev(ff, eo);
  return ev.compute(x, c, f);
}

double energy(const ForceField& ff, const std::vector<double>& x, const Cell& c) {
  std::vector<double> f;
  return eval(ff, x, c, f).total();
}

// The model written out directly (all pairs by the minimum image): the reference for the evaluator's sums.
double reference(const std::vector<int>& cls, const SwapModel::Pair& p, const std::vector<double>& x, const Cell& c) {
  const double rmin = p.sigma * std::pow(2.0, 1.0 / (p.n - p.m));
  auto v = [&](double r) { return 4 * p.eps * (std::pow(p.sigma / r, p.n) - std::pow(p.sigma / r, p.m)); };
  const double vs = p.shift ? v(p.rc) : 0;
  auto hat = [&](double r) { return r <= rmin ? 1.0 : -(v(r) - vs) / p.eps; };
  const size_t n = cls.size();
  auto dist = [&](size_t i, size_t j) {
    Vec3 d{x[3 * j] - x[3 * i], x[3 * j + 1] - x[3 * i + 1], x[3 * j + 2] - x[3 * i + 2]};
    return norm(c.minimum_image(d));
  };
  auto reacts = [&](size_t i, size_t j) { return cls[i] && cls[j] && ((cls[i] == p.a && cls[j] == p.b) || (cls[i] == p.b && cls[j] == p.a)); };
  double e = 0;
  for (size_t i = 0; i < n; ++i)
    for (size_t j = i + 1; j < n; ++j)
      if (reacts(i, j) && dist(i, j) < p.rc) e += v(dist(i, j)) - vs;
  for (size_t i = 0; i < n; ++i)
    for (size_t j = 0; j < n; ++j)
      for (size_t k = j + 1; k < n; ++k) {
        if (j == i || k == i || !reacts(i, j) || !reacts(i, k)) continue;
        const double rj = dist(i, j), rk = dist(i, k);
        if (rj < p.rc && rk < p.rc) e += p.lambda * p.eps * hat(rj) * hat(rk);
      }
  return e;
}

}  // namespace

// A bound pair sits at −ε; a third site joining it at the minimum makes −2ε + λε: the swap passes over ε(λ − 1).
TEST(SwapModel, BarrierIsEpsilonTimesLambdaMinusOne) {
  SwapModel::Pair p;
  p.eps = 10, p.sigma = 1, p.rc = 2.5, p.lambda = 1.2;
  const double rmin = std::pow(2.0, 0.1);
  const Cell c = box(30);
  // B (class 2) at the centre, A (class 1) at the minimum, the other A far away
  ForceField ff = sites({1, 2, 1}, p);
  std::vector<double> bound{10, 10, 10, 10 + rmin, 10, 10, 20, 20, 20};
  EXPECT_NEAR(energy(ff, bound, c), -p.eps, 1e-9);
  // both A at the minimum of B, on opposite sides (A–A do not interact)
  std::vector<double> triplet{10 - rmin, 10, 10, 10, 10, 10, 10 + rmin, 10, 10};
  EXPECT_NEAR(energy(ff, triplet, c), -2 * p.eps + p.lambda * p.eps, 1e-9);
  EXPECT_NEAR(energy(ff, triplet, c) - energy(ff, bound, c), p.eps * (p.lambda - 1), 1e-9);
  // λ = 1: no barrier, the triplet costs nothing
  p.lambda = 1;
  ForceField flat = sites({1, 2, 1}, p);
  EXPECT_NEAR(energy(flat, triplet, c), energy(flat, bound, c), 1e-9);
}

// Forces are the energy's gradient; the virial tensor is its response to strain (stretch and shear); the evaluator's sums
// are the model's (a direct reference over all pairs and triplets, 200 sites binned by the cell list).
TEST(SwapModel, ForcesVirialAndSumsMatchTheEnergy) {
  SwapModel::Pair p;
  p.eps = 5, p.sigma = 1, p.rc = 2.5, p.lambda = 1.5;
  std::mt19937_64 rng(4);
  std::uniform_real_distribution<double> u(0, 1);
  for (bool shift : {false, true}) {
    p.shift = shift;
    // a small cluster in a box: every site within reach of several others, none closer than 0.95
    const Cell c = box(9);
    std::vector<int> cls;
    std::vector<double> x;
    while (cls.size() < 10) {
      const Vec3 r{3 + 3 * u(rng), 3 + 3 * u(rng), 3 + 3 * u(rng)};
      bool ok = true;
      for (size_t k = 0; k < cls.size(); ++k) ok = ok && norm(r - Vec3{x[3 * k], x[3 * k + 1], x[3 * k + 2]}) > 0.95;
      if (!ok) continue;
      cls.push_back(cls.size() % 2 ? 2 : 1);
      x.insert(x.end(), {r[0], r[1], r[2]});
    }
    ForceField ff = sites(cls, p);
    std::vector<double> f;
    const EnergyTerms t = eval(ff, x, c, f);
    EXPECT_NEAR(t.total(), reference(cls, p, x, c), 1e-9 * (1 + std::fabs(t.total())));
    const double h = 1e-6;
    for (size_t k = 0; k < x.size(); ++k) {
      auto xp = x, xm = x;
      xp[k] += h, xm[k] -= h;
      const double g = (energy(ff, xp, c) - energy(ff, xm, c)) / (2 * h);
      EXPECT_NEAR(f[k], -g, 1e-5 * (1 + std::fabs(g))) << "coordinate " << k;
    }
    // affine strains: x along x stretched by (1 + δ); x moved by δ y (shear)
    auto strained = [&](double dl, bool shear, Cell& cs) {
      std::vector<double> xs = x;
      cs = c;
      for (size_t i = 0; i < cls.size(); ++i) xs[3 * i] = shear ? x[3 * i] + dl * x[3 * i + 1] : x[3 * i] * (1 + dl);
      if (shear) cs.b[0] += dl * c.b[1];
      else cs.a[0] *= 1 + dl;
      return xs;
    };
    for (bool shear : {false, true}) {
      Cell cp, cm;
      const auto xp = strained(1e-6, shear, cp), xm = strained(-1e-6, shear, cm);
      const double dEdl = (energy(ff, xp, cp) - energy(ff, xm, cm)) / 2e-6;
      EXPECT_NEAR(shear ? t.w[3] : t.w[0], -dEdl, 1e-5 * (1 + std::fabs(dEdl))) << (shear ? "xy" : "xx");
    }
    EXPECT_NEAR(t.virial, t.w[0] + t.w[1] + t.w[2], 1e-9 * (1 + std::fabs(t.virial)));
  }
  // many sites, binned: the same sums as the direct reference
  const Cell c = box(15);
  std::vector<int> cls;
  std::vector<double> x;
  while (cls.size() < 200) {
    const Vec3 r{15 * u(rng), 15 * u(rng), 15 * u(rng)};
    bool ok = true;
    for (size_t k = 0; k < cls.size() && ok; ++k) ok = norm(c.minimum_image(r - Vec3{x[3 * k], x[3 * k + 1], x[3 * k + 2]})) > 0.95;
    if (!ok) continue;
    cls.push_back(1 + int(cls.size() % 2));
    x.insert(x.end(), {r[0], r[1], r[2]});
  }
  p.shift = false;
  ForceField ff = sites(cls, p);
  const double e = energy(ff, x, c);
  EXPECT_NEAR(e, reference(cls, p, x, c), 1e-9 * (1 + std::fabs(e)));
  EXPECT_LT(e, 0);
}

// Excluded (bonded) pairs take no part; a cut-off over half the cell is refused.
TEST(SwapModel, ExclusionsAndCellSize) {
  SwapModel::Pair p;
  p.eps = 10, p.sigma = 1, p.rc = 2.5;
  const double rmin = std::pow(2.0, 0.1);
  ForceField ff = sites({1, 2}, p);
  std::vector<double> x{10, 10, 10, 10 + rmin, 10, 10};
  EXPECT_NEAR(energy(ff, x, box(30)), -10, 1e-9);
  ff.excluded = {{1}, {0}};
  EXPECT_NEAR(energy(ff, x, box(30)), 0, 1e-12);
  std::vector<double> f;
  EXPECT_THROW(eval(sites({1, 2}, p), x, box(4), f), FieldError);
}

// NVE dynamics of free A and B sites that pair and swap (λ = 1: no barrier), with a WCA core between like sites: the
// pair shifted to zero at the cut-off conserves the energy, triplets included.
TEST(SwapModel, NveConservesEnergy) {
  SwapModel::Pair p;
  p.eps = 3, p.sigma = 1, p.rc = 2.5, p.lambda = 1, p.shift = true;
  System s;
  const double L = 12;
  s.cell = box(L);
  std::vector<int> cls;
  for (int ix = 0; ix < 6; ++ix)
    for (int iy = 0; iy < 6; ++iy)
      for (int iz = 0; iz < 6; ++iz) {
        Atom a;
        a.element = 6;
        a.pos = {(ix + 0.5) * L / 6, (iy + 0.5) * L / 6, (iz + 0.5) * L / 6};
        s.atoms.push_back(a);
        cls.push_back((ix + iy + iz) % 2 ? 2 : 1);
      }
  ForceField ff = sites(cls, p);
  // like sites repel (WCA, ε = 1 kcal/mol, σ = 1.6 Å); unlike sites only through the swap model
  ff.type_names = {"A", "B"};
  ff.lj = {{1.0, 1.6}, {1.0, 1.6}};
  ff.atom_type.clear();
  for (size_t i = 0; i < cls.size(); ++i) ff.type_index[i] = cls[i] - 1, ff.atom_type.push_back(cls[i] == 1 ? "A" : "B");
  ff.pair_override[{0, 1}] = {0.0, 1.0};
  ff.lj_shift = true;
  ff.cutoff = 1.6 * std::pow(2.0, 1.0 / 6);
  DynamicsOptions d;
  d.field = std::make_shared<ForceField>(ff);
  d.thermostat = Thermostat::None;
  d.temperature = 300;
  d.dt = 0.5;
  d.steps = 4000;
  d.thermo_every = 50;
  d.energy.coulomb = false;
  d.energy.tail = false;
  DynamicsReport r;
  run_dynamics(s, d, &r);
  double drift = 0;
  for (const auto& row : r.thermo) drift = std::max(drift, std::fabs(row.conserved - r.thermo[0].conserved));
  const double kin = r.thermo[0].kinetic;
  EXPECT_LT(drift, 0.01 * kin) << "drift " << drift << " kcal/mol, kinetic energy " << kin;
  EXPECT_LT(r.thermo.back().potential, r.thermo[0].potential);   // pairs formed
}
