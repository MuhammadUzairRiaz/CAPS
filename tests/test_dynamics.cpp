#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <map>
#include <numeric>

#include "caps/dynamics.hpp"
#include "caps/grow.hpp"
#include "caps/io.hpp"
#include "caps/relax.hpp"
#include "caps/uff.hpp"

using namespace caps;

namespace {

// A small relaxed polystyrene cell (3 chains of DP 4, 0.9 g/cm³), built once.
const System& relaxed_cell() {
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

double mean(const std::vector<ThermoRow>& rows, size_t from, double ThermoRow::*field) {
  double s = 0;
  size_t n = 0;
  for (size_t k = from; k < rows.size(); ++k, ++n) s += rows[k].*field;
  return n ? s / n : 0;
}

}  // namespace

TEST(Dynamics, InitialVelocitiesHaveTargetTemperatureAndNoDrift) {
  System s = relaxed_cell();
  DynamicsOptions o;
  o.steps = 0;
  o.temperature = 350;
  DynamicsReport r;
  run_dynamics(s, o, &r);
  ASSERT_EQ(r.thermo.size(), 1u);
  EXPECT_NEAR(r.thermo[0].temperature, 350.0, 1e-9);
  ASSERT_EQ(s.velocities.size(), s.atoms.size());
  const ForceField ff = assign_gaff(s);
  double p[3] = {0, 0, 0};
  for (size_t i = 0; i < s.atoms.size(); ++i)
    for (int k = 0; k < 3; ++k) p[k] += ff.mass[i] * s.velocities[i][k];
  EXPECT_LT(std::fabs(p[0]) + std::fabs(p[1]) + std::fabs(p[2]), 1e-9);
}

// Velocity Verlet conserves energy: after the first 100 fs the total stays within a small band and does not drift.
TEST(Dynamics, NveConservesEnergy) {
  System s = relaxed_cell();
  DynamicsOptions o;
  o.thermostat = Thermostat::None;
  o.dt = 0.5;
  o.steps = 2200;
  o.thermo_every = 20;
  o.temperature = 300;
  DynamicsReport r;
  run_dynamics(s, o, &r);
  const auto& t = r.thermo;
  const size_t from = 10;   // skip 100 fs
  double lo = 1e300, hi = -1e300;
  for (size_t k = from; k < t.size(); ++k) { lo = std::min(lo, t[k].total); hi = std::max(hi, t[k].total); }
  const double kin = mean(t, from, &ThermoRow::kinetic);
  EXPECT_LT(hi - lo, 0.01 * kin) << "band " << hi - lo << " kcal/mol, kinetic " << kin;
  // drift: first vs last fifth
  const size_t q = (t.size() - from) / 5;
  double a = 0, b = 0;
  for (size_t k = 0; k < q; ++k) { a += t[from + k].total; b += t[t.size() - 1 - k].total; }
  EXPECT_LT(std::fabs(b - a) / q, 0.005 * kin);
}

TEST(Dynamics, ThermostatsHoldTheTemperature) {
  for (Thermostat th : {Thermostat::Bussi, Thermostat::Langevin}) {
    System s = relaxed_cell();
    DynamicsOptions o;
    o.thermostat = th;
    o.temperature = 400;
    o.tau_t = 50;
    o.steps = 4000;
    o.thermo_every = 10;
    o.seed = 11;
    DynamicsReport r;
    run_dynamics(s, o, &r);
    const double tm = mean(r.thermo, 50, &ThermoRow::temperature);
    EXPECT_NEAR(tm, 400.0, 16.0) << to_string(th);
  }
}

// The conserved quantity of the Bussi thermostat (total energy plus heat given to the bath) does not drift.
TEST(Dynamics, BussiConservedQuantity) {
  System s = relaxed_cell();
  DynamicsOptions o;
  o.thermostat = Thermostat::Bussi;
  o.dt = 0.5;
  o.steps = 2000;
  o.thermo_every = 20;
  DynamicsReport r;
  run_dynamics(s, o, &r);
  const auto& t = r.thermo;
  double lo = 1e300, hi = -1e300;
  for (size_t k = 10; k < t.size(); ++k) { lo = std::min(lo, t[k].conserved); hi = std::max(hi, t[k].conserved); }
  EXPECT_LT(hi - lo, 0.01 * mean(t, 10, &ThermoRow::kinetic));
}

TEST(Dynamics, BarostatRespondsToTargetPressure) {
  double rho[2];
  int k = 0;
  for (double p0 : {1.0, 20000.0}) {
    System s = relaxed_cell();
    DynamicsOptions o;
    o.barostat = Barostat::Berendsen;
    o.pressure = p0;
    o.tau_p = 200;
    o.steps = 2000;
    o.seed = 4;
    DynamicsReport r;
    run_dynamics(s, o, &r);
    rho[k++] = s.density();
    EXPECT_NE(s.cell.volume(), relaxed_cell().cell.volume());
  }
  EXPECT_GT(rho[1], rho[0] + 0.02);
}

TEST(Dynamics, SameSeedSameTrajectory) {
  System a = relaxed_cell(), b = relaxed_cell();
  DynamicsOptions o;
  o.thermostat = Thermostat::Langevin;
  o.barostat = Barostat::CRescale;
  o.steps = 300;
  o.seed = 9;
  run_dynamics(a, o);
  run_dynamics(b, o);
  for (size_t i = 0; i < a.atoms.size(); ++i)
    for (int k = 0; k < 3; ++k) ASSERT_EQ(a.atoms[i].pos[k], b.atoms[i].pos[k]);
}

TEST(Dynamics, VelocitiesAndTrajectoryRoundTrip) {
  System s = relaxed_cell();
  DynamicsOptions o;
  o.steps = 100;
  o.frame_every = 50;
  Trajectory t;
  t.topology = s;
  o.frame = [&](const std::vector<double>& x, const Cell& c, int64_t step) {
    std::vector<Vec3> p(x.size() / 3);
    for (size_t i = 0; i < p.size(); ++i) p[i] = {x[3 * i], x[3 * i + 1], x[3 * i + 2]};
    t.positions.push_back(p);
    t.cells.push_back(c);
    t.timesteps.push_back(step);
  };
  run_dynamics(s, o);
  EXPECT_EQ(t.frames(), 3u);   // steps 0, 50, 100
  const auto dir = std::filesystem::temp_directory_path();
  const std::string data = (dir / "caps_md_vel.data").string(), dump = (dir / "caps_md.lammpstrj").string();
  write_lammps_data_ff(s, assign_gaff(s), EnergyOptions{}, data);
  const System back = read_lammps_data(data);
  ASSERT_EQ(back.velocities.size(), s.velocities.size());
  for (size_t i = 0; i < s.atoms.size(); ++i)
    for (int k = 0; k < 3; ++k) EXPECT_NEAR(back.velocities[i][k], s.velocities[i][k], 1e-8);
  write_lammps_dump(t, dump);
  const Trajectory tb = read_lammps_dump(dump, &back);
  ASSERT_EQ(tb.frames(), 3u);
  EXPECT_EQ(tb.timesteps[2], 100);
  EXPECT_NEAR(tb.positions[2][7][1], t.positions[2][7][1], 1e-4);
}

TEST(Dynamics, RejectsBadSettings) {
  System s = relaxed_cell();
  DynamicsOptions o;
  o.dt = 10;
  EXPECT_THROW(run_dynamics(s, o), std::invalid_argument);
  o.dt = 1;
  o.thermostat = Thermostat::None;
  o.barostat = Barostat::CRescale;
  EXPECT_THROW(run_dynamics(s, o), std::invalid_argument);
}

// NPH: the Berendsen barostat with no thermostat (stochastic cell rescaling needs one and is refused): the cell follows
// the target pressure while the temperature is left to the dynamics
TEST(Dynamics, NphScalesTheCellWithoutAThermostat) {
  System s = relaxed_cell();
  DynamicsOptions o;
  o.thermostat = Thermostat::None;
  o.barostat = Barostat::Berendsen;
  o.pressure = 20000.0;
  o.tau_p = 200;
  o.steps = 2000;
  o.seed = 4;
  DynamicsReport r;
  run_dynamics(s, o, &r);
  EXPECT_GT(s.density(), relaxed_cell().density() + 0.01);
  const double t = mean(r.thermo, r.thermo.size() / 2, &ThermoRow::temperature);
  EXPECT_GT(t, 100);
  EXPECT_LT(t, 900);
  System c = relaxed_cell();
  o.barostat = Barostat::CRescale;
  EXPECT_THROW(run_dynamics(c, o), std::invalid_argument);
}

// r-RESPA: bonded forces at 0.5 fs inside 2 fs steps conserve energy as plain 0.5 fs steps do (and far better than
// plain 2 fs steps, which the C–H stretches do not allow); respa = 1 is the ordinary integrator
TEST(Dynamics, RespaConservesEnergyWithLongOuterSteps) {
  auto band = [](int respa, double dt, int steps) {
    System s = relaxed_cell();
    DynamicsOptions o;
    o.thermostat = Thermostat::None;
    o.dt = dt;
    o.respa = respa;
    o.steps = steps;
    o.thermo_every = 5;
    o.seed = 7;
    o.new_velocities = true;
    DynamicsReport r;
    run_dynamics(s, o, &r);
    double lo = 1e300, hi = -1e300;
    for (size_t k = 4; k < r.thermo.size(); ++k) lo = std::min(lo, r.thermo[k].total), hi = std::max(hi, r.thermo[k].total);
    return std::make_pair(hi - lo, mean(r.thermo, 4, &ThermoRow::kinetic));
  };
  const auto [plain_fine, kin] = band(1, 0.5, 1200);
  const auto [respa, kin2] = band(4, 2.0, 300);
  const auto [plain_coarse, kin3] = band(1, 2.0, 300);
  EXPECT_LT(respa, 0.02 * kin2) << "r-RESPA band " << respa;
  EXPECT_LT(respa, 0.5 * plain_coarse) << "r-RESPA " << respa << " vs plain 2 fs " << plain_coarse;
  EXPECT_LT(plain_fine, 0.01 * kin);
  (void)kin3;
  System c = relaxed_cell();
  DynamicsOptions o;
  o.respa = 2;
  o.thermostat = Thermostat::Langevin;
  EXPECT_THROW(run_dynamics(c, o), std::invalid_argument);
}

// SHAKE/RATTLE: bonds to hydrogen stay at the force field's r0 to the tolerance, with no relative velocity along them;
// 2 fs steps conserve energy better than plain 1 fs steps (plain 2 fs steps do not); each constraint takes one degree of freedom
TEST(Dynamics, BondConstraintsHoldAndConserveEnergy) {
  const System& c0 = relaxed_cell();
  const ForceField ff = assign_gaff(c0);
  auto run = [&](ConstraintMode mode, double dt, int steps, System& s) {
    s = c0;
    DynamicsOptions o;
    o.thermostat = Thermostat::None;
    o.dt = dt;
    o.constraints = mode;
    o.steps = steps;
    o.thermo_every = 5;
    o.seed = 7;
    o.new_velocities = true;
    DynamicsReport r;
    run_dynamics(s, o, &r);
    double lo = 1e300, hi = -1e300;
    for (size_t k = 4; k < r.thermo.size(); ++k) lo = std::min(lo, r.thermo[k].total), hi = std::max(hi, r.thermo[k].total);
    return std::make_pair(hi - lo, mean(r.thermo, 4, &ThermoRow::kinetic));
  };
  System s;
  const auto [band, kin] = run(ConstraintMode::HBonds, 2.0, 300, s);
  System u;
  const auto [plain1, kin1] = run(ConstraintMode::None, 1.0, 600, u);
  const auto [plain2, kin2] = run(ConstraintMode::None, 2.0, 300, u);
  EXPECT_LT(band, 0.03 * kin) << "SHAKE 2 fs band " << band;
  EXPECT_LT(band, plain1) << "SHAKE 2 fs " << band << " vs plain 1 fs " << plain1;
  EXPECT_LT(band, 0.3 * plain2) << "SHAKE 2 fs " << band << " vs plain 2 fs " << plain2;
  (void)kin1, (void)kin2;
  // every C–H at r0, no velocity along it
  std::map<std::pair<uint32_t, uint32_t>, double> r0;
  for (const auto& b : ff.bonds) r0[{std::min(b.i, b.j), std::max(b.i, b.j)}] = b.r0;
  size_t nh = 0;
  double worst = 0, worst_v = 0;
  for (const auto& b : s.bonds) {
    if (s.atoms[b.i].element != 1 && s.atoms[b.j].element != 1) continue;
    ++nh;
    const Vec3 d = s.cell.minimum_image(s.atoms[b.i].pos - s.atoms[b.j].pos);
    worst = std::max(worst, std::fabs(norm(d) - r0[{std::min(b.i, b.j), std::max(b.i, b.j)}]));
    worst_v = std::max(worst_v, std::fabs(dot(d, s.velocities[b.i] - s.velocities[b.j])) / norm(d));
  }
  ASSERT_GT(nh, 0u);
  EXPECT_LT(worst, 1e-6);
  EXPECT_LT(worst_v, 1e-7);
  // degrees of freedom: 3N − 3 − constraints
  System t = c0;
  DynamicsOptions o;
  o.steps = 0;
  o.temperature = 300;
  o.constraints = ConstraintMode::HBonds;
  DynamicsReport r;
  run_dynamics(t, o, &r);
  const double k = kinetic_energy([&] {
    std::vector<double> v;
    for (const auto& q : t.velocities) v.insert(v.end(), {q[0], q[1], q[2]});
    return v;
  }(), ff.mass);
  EXPECT_NEAR(2 * k / ((3.0 * t.atoms.size() - 3 - nh) * 0.0019872041), 300.0, 1e-6);
  o.respa = 2;
  EXPECT_THROW(run_dynamics(t, o), std::invalid_argument);
}

// Nosé–Hoover chains and MTK pressure coupling: the temperature settles at the target, the extended system's energy
// does not drift (0.5 fs steps), and a high target pressure packs the cell denser than 1 atm does. Step-by-step
// agreement with LAMMPS's fix nvt / npt iso is checked by scripts/check_lammps.sh.
TEST(Dynamics, NoseHooverChainsAndMtkBarostat) {
  {
    System s = relaxed_cell();
    DynamicsOptions o;
    o.thermostat = Thermostat::NoseHoover;
    o.temperature = 400;
    o.tau_t = 50;
    o.steps = 4000;
    o.thermo_every = 10;
    o.seed = 11;
    DynamicsReport r;
    run_dynamics(s, o, &r);
    EXPECT_NEAR(mean(r.thermo, 100, &ThermoRow::temperature), 400.0, 20.0);
  }
  {
    System s = relaxed_cell();
    DynamicsOptions o;
    o.thermostat = Thermostat::NoseHoover;
    o.barostat = Barostat::MTK;
    o.dt = 0.5;
    o.steps = 2000;
    o.thermo_every = 20;
    DynamicsReport r;
    run_dynamics(s, o, &r);
    double lo = 1e300, hi = -1e300;
    for (size_t k = 10; k < r.thermo.size(); ++k) lo = std::min(lo, r.thermo[k].conserved), hi = std::max(hi, r.thermo[k].conserved);
    EXPECT_LT(hi - lo, 0.02 * mean(r.thermo, 10, &ThermoRow::kinetic));
  }
  double rho[2];
  int k = 0;
  for (double p0 : {1.0, 20000.0}) {
    System s = relaxed_cell();
    DynamicsOptions o;
    o.thermostat = Thermostat::NoseHoover;
    o.barostat = Barostat::MTK;
    o.pressure = p0;
    o.tau_p = 200;
    o.steps = 2000;
    DynamicsReport r;
    run_dynamics(s, o, &r);
    rho[k++] = s.density();
  }
  EXPECT_GT(rho[1], rho[0] + 0.02);
  // MTK is isotropic and runs without constraints
  System s = relaxed_cell();
  DynamicsOptions bad;
  bad.barostat = Barostat::MTK;
  bad.thermostat = Thermostat::Bussi;
  EXPECT_THROW(run_dynamics(s, bad), std::invalid_argument);
}

// LINCS solves the same constraint equations as SHAKE: the two trajectories agree to the tolerance, energy is
// conserved as well, and every constrained bond holds its length
TEST(Dynamics, LincsMatchesShake) {
  const System& c0 = relaxed_cell();
  auto run = [&](ConstraintAlgorithm alg, DynamicsReport& r) {
    System s = c0;
    DynamicsOptions o;
    o.thermostat = Thermostat::None;
    o.dt = 2.0;
    o.constraints = ConstraintMode::HBonds;
    o.constraint_algorithm = alg;
    o.steps = 200;
    o.thermo_every = 5;
    o.seed = 7;
    o.new_velocities = true;
    run_dynamics(s, o, &r);
    return s;
  };
  DynamicsReport rs, rl;
  const System a = run(ConstraintAlgorithm::Shake, rs), b = run(ConstraintAlgorithm::Lincs, rl);
  double worst = 0;
  for (size_t i = 0; i < a.atoms.size(); ++i) worst = std::max(worst, norm(a.atoms[i].pos - b.atoms[i].pos));
  EXPECT_LT(worst, 1e-4);
  EXPECT_NEAR(rl.thermo.back().total, rs.thermo.back().total, 1e-3 * std::fabs(rs.thermo.back().kinetic));
  EXPECT_TRUE(std::any_of(rl.notes.begin(), rl.notes.end(), [](const std::string& n) { return n.find("LINCS") != std::string::npos; }));
}

// An electric field: two +1 e argon atoms at rest, 15 Å apart (beyond the cut-off, Coulomb off), nothing else acting, moves ½ a t² along E with
// a = q E · 23.0605 kcal/mol/Å per e·V/Å over its mass; a neutral atom does not move.
TEST(Dynamics, ElectricFieldAcceleratesACharge) {
  for (double q : {1.0, 0.0}) {
    System s;
    Atom a;
    a.element = 18;
    a.pos = {10, 10, 10};
    s.atoms.push_back(a);
    a.pos = {25, 10, 10};
    s.atoms.push_back(a);
    s.cell.a = {40, 0, 0}, s.cell.b = {0, 40, 0}, s.cell.c = {0, 0, 40};
    s.velocities = {Vec3{0, 0, 0}, Vec3{0, 0, 0}};
    auto ff = std::make_shared<ForceField>(default_forcefield(s));
    ff->charge = {q, q};
    DynamicsOptions o;
    o.field = ff;
    o.thermostat = Thermostat::None;
    o.steps = 100;
    o.dt = 1;
    o.frame_every = 0;
    o.energy.coulomb = false;
    o.energy.cutoff = 10;
    o.efield = {0, 0, 0.1};
    DynamicsReport rep;
    run_dynamics(s, o, &rep);
    const double acc = q * 0.1 * 23.060548 * 4.184e-4 / 39.948;   // Å/fs²
    for (const auto& at : s.atoms) EXPECT_NEAR(at.pos[2] - 10, 0.5 * acc * 100 * 100, 1e-6);
    EXPECT_NEAR(s.atoms[0].pos[0] - 10, 0, 1e-9);
  }
}

// SLLOD shear of the Lennard-Jones liquid near its triple point (ρ* = 0.844, T* = 0.722, argon units ε/k_B = 119.8 K,
// σ = 3.405 Å): at γ̇* = 1 the viscosity η* = −⟨P_xy⟩/γ̇ in units of √(mε)/σ² is shear-thinned to about 2.0–2.4 (Evans &
// Morriss; zero-shear ≈ 3.2). 256 atoms, 2,000 steps to melt and 3,000 sheared: the range allows for the noise. The cell
// tilts at γ̇ L_y and is flipped back.
TEST(Dynamics, SllodShearViscosityOfTheLjLiquid) {
  const double sigma = 3.405, eps = 119.8 * 0.0019872043, rho = 0.844;
  const double a = std::cbrt(4 / rho) * sigma;
  System s;
  for (int i = 0; i < 4; ++i)
    for (int j = 0; j < 4; ++j)
      for (int k = 0; k < 4; ++k)
        for (const auto& b : std::vector<Vec3>{{0, 0, 0}, {0.5, 0.5, 0}, {0.5, 0, 0.5}, {0, 0.5, 0.5}}) {
          Atom at;
          at.element = 18;
          at.mol = int64_t(s.atoms.size() + 1);
          at.pos = {(i + b[0]) * a, (j + b[1]) * a, (k + b[2]) * a};
          s.atoms.push_back(at);
        }
  s.cell.a = {4 * a, 0, 0}, s.cell.b = {0, 4 * a, 0}, s.cell.c = {0, 0, 4 * a};
  auto ff = std::make_shared<ForceField>(default_forcefield(s));
  for (auto& t : ff->lj) t = {eps, sigma};
  std::fill(ff->charge.begin(), ff->charge.end(), 0.0);
  DynamicsOptions o;
  o.field = ff;
  o.temperature = 0.722 * 119.8;
  o.thermostat = Thermostat::Bussi;
  o.tau_t = 100;
  o.dt = 5;
  o.energy.cutoff = 2.5 * sigma;
  o.energy.coulomb = false;
  o.energy.tail = false;
  o.frame_every = 0;
  o.steps = 2000;
  run_dynamics(s, o);
  const double tau = sigma * 1e-10 * std::sqrt(39.948e-3 / 6.02214076e23 / (119.8 * 1.380649e-23));   // s
  o.shear_rate = 1.0 / tau * 1e-12;   // γ̇* = 1, in 1/ps
  o.tau_t = 10;                        // the viscous heat (~80 K/ps here) leaves T above the target by about heat rate × τ
  o.steps = 3000;
  o.thermo_every = 10;
  DynamicsReport rep;
  run_dynamics(s, o, &rep);
  double pxy = 0;
  size_t np = 0;
  for (size_t k = rep.thermo.size() / 4; k < rep.thermo.size(); ++k) pxy += rep.thermo[k].p[3], ++np;
  pxy /= double(np);
  const double eta = -pxy * 101325 / (o.shear_rate * 1e12);                     // Pa·s
  const double eta_star = eta / (std::sqrt(39.948e-3 / 6.02214076e23 * 119.8 * 1.380649e-23) / std::pow(sigma * 1e-10, 2));
  EXPECT_GT(eta_star, 1.5);
  EXPECT_LT(eta_star, 3.0);
  // the tilt: γ̇ t L_y, wrapped into ±L_x/2
  const double L = 4 * a, tilt = std::fmod(o.shear_rate * 1e-3 * o.steps * o.dt * L + 0.5 * L, L) - 0.5 * L;
  EXPECT_NEAR(s.cell.b[0], tilt, 1e-6 * L);
  double tm = 0;
  for (size_t k = rep.thermo.size() / 4; k < rep.thermo.size(); ++k) tm += rep.thermo[k].temperature / double(np);
  EXPECT_NEAR(tm, o.temperature, 0.03 * o.temperature);   // the thermostat takes the viscous heat out of the peculiar motion
  std::printf("SLLOD LJ liquid: eta* = %.3f at gamma* = 1 (<T> %.1f K)\n", eta_star, tm);
}
