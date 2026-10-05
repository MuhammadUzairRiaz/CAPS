#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <map>
#include <numeric>
#include <random>

#include "caps/dynamics.hpp"
#include "caps/grow.hpp"
#include "caps/io.hpp"
#include "caps/relax.hpp"
#include "caps/ibi.hpp"
#include "caps/free_energy.hpp"
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

// Atoms held along some axes only (fixed bits 2 x, 4 y, 8 z): in MD their z stays exactly while x and y move, the
// temperature counts one coordinate fewer each; in a minimisation the held x stays.
TEST(Dynamics, PartlyHeldAtoms) {
  System s = relaxed_cell();
  const size_t n = s.atoms.size(), k = 20;
  DynamicsOptions o;
  o.fixed.assign(n, 0);
  for (size_t i = 0; i < k; ++i) o.fixed[i] = 8;
  o.steps = 300;
  o.frame_every = 0;
  o.temperature = 300;
  System t = s;
  DynamicsReport rep;
  run_dynamics(t, o, &rep);
  double moved_xy = 0;
  for (size_t i = 0; i < k; ++i) {
    EXPECT_DOUBLE_EQ(t.atoms[i].pos[2], s.atoms[i].pos[2]);
    moved_xy += std::fabs(t.atoms[i].pos[0] - s.atoms[i].pos[0]) + std::fabs(t.atoms[i].pos[1] - s.atoms[i].pos[1]);
  }
  EXPECT_GT(moved_xy, 0.1);
  bool noted = false;
  for (const auto& note : rep.notes) noted = noted || note.find("20 single coordinates held") != std::string::npos;
  EXPECT_TRUE(noted);
  RelaxOptions r;
  r.fixed.assign(n, 0);
  for (size_t i = 0; i < k; ++i) r.fixed[i] = 2;
  r.ftol = 2.0;
  System u = t;
  relax(u, r);
  for (size_t i = 0; i < k; ++i) EXPECT_DOUBLE_EQ(u.atoms[i].pos[0], t.atoms[i].pos[0]);
  EXPECT_TRUE(holds_axis(1, 0) && holds_axis(8, 2) && !holds_axis(8, 0) && holds_all(14) && !holds_all(12));
}

// A tabulated pair: the Lennard-Jones potential written into a table (0.005 Å grid) gives the analytic energy and
// forces to the interpolation error.
TEST(Dynamics, TabulatedPairMatchesLennardJones) {
  System s = relaxed_cell();
  ForceField ff = default_forcefield(s);
  std::fill(ff.charge.begin(), ff.charge.end(), 0.0);
  EnergyOptions e;
  e.coulomb = false, e.tail = false, e.cutoff = 9;
  std::vector<double> x, fa, ft;
  for (const auto& a : s.atoms) x.insert(x.end(), a.pos.begin(), a.pos.end());
  Evaluator ea(ff, e);
  const auto Ea = ea.compute(x, s.cell, fa);
  // every type pair as a table of its mixed LJ, shifted to zero at the cut-off as the plain pairs are
  ForceField ft_ff = ff;
  const int nt = int(ff.type_names.size());
  for (int a = 0; a < nt; ++a)
    for (int b = a; b < nt; ++b) {
      const PairType p = mixed_pair(ff, a, b);
      TabulatedPair tb;
      tb.r0 = 1.0, tb.dr = 0.005;
      auto V = [&](double r) { const double q = std::pow(p.sigma / r, 6); return 4 * p.eps * (q * q - q); };
      for (double r = tb.r0; r <= 9.0 + 1e-9; r += tb.dr) {
        const double q = std::pow(p.sigma / r, 6);
        tb.e.push_back(V(r) - V(9.0));
        tb.f.push_back(24 * p.eps * (2 * q * q - q) / r);
      }
      ft_ff.tables.push_back(tb);
      ft_ff.pair_func[{a, b}] = {kPairTable, double(ft_ff.tables.size() - 1), 0, 0};
    }
  Evaluator et(ft_ff, e);
  const auto Et = et.compute(x, s.cell, ft);
  EXPECT_NEAR(Et.vdw, Ea.vdw, 1e-3 * std::fabs(Ea.vdw) + 0.05);
  double worst = 0;
  for (size_t k = 0; k < fa.size(); ++k) worst = std::max(worst, std::fabs(ft[k] - fa[k]));
  EXPECT_LT(worst, 0.05);
}

// Iterative Boltzmann inversion: the Lennard-Jones liquid's own g(r) as the target, started from −k_B T ln g (which
// over-structures a dense liquid); the iterations bring the model's g(r) back toward the target.
TEST(Dynamics, IbiRecoversTheLjLiquidStructure) {
  const double sigma = 3.405, eps = 119.8 * 0.0019872043, rho = 0.844, a = std::cbrt(4 / rho) * sigma;
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
  s.has_mol = true;
  s.cell.a = {4 * a, 0, 0}, s.cell.b = {0, 4 * a, 0}, s.cell.c = {0, 0, 4 * a};
  auto ff = std::make_shared<ForceField>(default_forcefield(s));
  for (auto& t : ff->lj) t = {eps, sigma};
  std::fill(ff->charge.begin(), ff->charge.end(), 0.0);
  // the target: the LJ liquid's g(r) at T* = 0.722
  DynamicsOptions d;
  d.field = ff;
  d.temperature = 0.722 * 119.8, d.dt = 5, d.tau_t = 100;
  d.energy.cutoff = 2.5 * sigma, d.energy.coulomb = false, d.energy.tail = false;
  d.steps = 2000, d.frame_every = 0;
  run_dynamics(s, d);
  std::vector<std::vector<Vec3>> frames;
  std::vector<Cell> cells;
  d.steps = 2000, d.frame_every = 20;
  d.frame = [&](const std::vector<double>& x, const Cell& c, int64_t step) {
    if (step == 0) return;
    std::vector<Vec3> p(x.size() / 3);
    for (size_t i = 0; i < p.size(); ++i) p[i] = {x[3 * i], x[3 * i + 1], x[3 * i + 2]};
    frames.push_back(p), cells.push_back(c);
  };
  run_dynamics(s, d);
  const double dr = 0.1;
  const size_t nbin = size_t(2.5 * sigma / dr);
  const auto g = nonbonded_gr(s, frames, cells, dr, nbin);
  std::vector<double> r(nbin);
  for (size_t k = 0; k < nbin; ++k) r[k] = (double(k) + 0.5) * dr;
  IbiOptions o;
  o.temperature = 0.722 * 119.8;
  o.iterations = 4;
  o.run_ps = 8, o.first_equilibrate_ps = 4, o.dt = 5;
  o.frame_every = 20;
  const auto res = run_ibi(s, *ff, r, g, o);
  ASSERT_EQ(res.history.size(), 5u);
  EXPECT_LT(res.history.back().residual, res.history.front().residual);
  EXPECT_LT(res.history.back().residual, 0.02);
  std::printf("IBI: residual %.4f -> %.4f\n", res.history.front().residual, res.history.back().residual);
}

// Alchemical pairs: at λ = 1 the soft-core Lennard-Jones and scaled Coulomb are the plain pairs; the ∂U/∂λ the kernel
// reports is the numerical derivative of the energy in λ (LJ at λ_lj = 0.5, Coulomb at λ_coul = 0.3).
TEST(Dynamics, AlchemicalPairsAndTheirDerivatives) {
  System s = relaxed_cell();
  ForceField ff = default_forcefield(s);
  for (size_t i = 0; i < ff.charge.size(); ++i) ff.charge[i] = (i % 2 ? 0.1 : -0.1);   // some charges
  std::vector<char> al(s.atoms.size(), 0);
  const auto mol = s.molecules();
  for (size_t i = 0; i < al.size(); ++i) al[i] = mol[i] == 0;   // the first chain is the solute
  std::vector<double> x, f;
  for (const auto& a : s.atoms) x.insert(x.end(), a.pos.begin(), a.pos.end());
  EnergyOptions e;
  e.tail = false;
  Evaluator plain(ff, e);
  const auto E0 = plain.compute(x, s.cell, f);
  EnergyOptions ea = e;
  ea.alchemical = &al;
  auto at = [&](double ll, double lc) {
    EnergyOptions o = ea;
    o.lambda_lj = ll, o.lambda_coul = lc;
    Evaluator ev(ff, o);
    std::vector<double> g;
    return ev.compute(x, s.cell, g);
  };
  const auto E1 = at(1, 1);
  EXPECT_NEAR(E1.vdw, E0.vdw, 1e-8 * std::fabs(E0.vdw) + 1e-8);
  EXPECT_NEAR(E1.coulomb, E0.coulomb, 1e-8 * std::fabs(E0.coulomb) + 1e-8);
  const double h = 1e-5;
  const auto Em = at(0.5, 0.3);
  const double num_lj = (at(0.5 + h, 0.3).vdw - at(0.5 - h, 0.3).vdw) / (2 * h);
  const double num_c = (at(0.5, 0.3 + h).coulomb - at(0.5, 0.3 - h).coulomb) / (2 * h);
  EXPECT_NEAR(Em.dudl_lj, num_lj, 1e-4 * std::fabs(num_lj) + 1e-6);
  EXPECT_NEAR(Em.dudl_coul, num_c, 1e-4 * std::fabs(num_c) + 1e-6);
  // λ = 0: the solute does not see the rest (its pairs with it carry nothing)
  const auto E00 = at(0, 0);
  EXPECT_NE(E00.vdw, E0.vdw);
  // the soft-core forces are −∇U: a solute atom moved by ±h along x
  EnergyOptions o = ea;
  o.lambda_lj = 0.4, o.lambda_coul = 0.6;
  Evaluator ev(ff, o);
  std::vector<double> g;
  ev.compute(x, s.cell, g);
  size_t k = 0;
  while (!al[k]) ++k;
  auto U = [&](double dx) {
    auto y = x;
    y[3 * k] += dx;
    Evaluator e2(ff, o);
    std::vector<double> gg;
    return e2.compute(y, s.cell, gg).total();
  };
  const double fnum = -(U(1e-5) - U(-1e-5)) / 2e-5;
  EXPECT_NEAR(g[3 * k], fnum, 1e-4 * std::fabs(fnum) + 1e-4);
}

// Solvation free energy by TI against Widom's test-particle insertion: the excess chemical potential of a Lennard-Jones
// fluid (ρ* = 0.5, T* = 2, the same shifted potential at 2.5 σ) — decoupling one particle by soft-core TI gives μ_ex.
TEST(Dynamics, SolvationTiMatchesWidomInsertion) {
  const double sigma = 3.405, eps = 119.8 * 0.0019872043, T = 2.0 * 119.8, rho = 0.5;
  const double a = std::cbrt(4 / rho) * sigma, L = 4 * a, rc = 2.5 * sigma;
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
  s.has_mol = true;
  s.cell.a = {L, 0, 0}, s.cell.b = {0, L, 0}, s.cell.c = {0, 0, L};
  auto ff = std::make_shared<ForceField>(default_forcefield(s));
  for (auto& t : ff->lj) t = {eps, sigma};
  std::fill(ff->charge.begin(), ff->charge.end(), 0.0);
  // Widom: random insertions into frames of the fluid
  DynamicsOptions d;
  d.field = ff, d.temperature = T, d.dt = 5, d.tau_t = 100;
  d.energy.cutoff = rc, d.energy.coulomb = false, d.energy.tail = false;
  d.steps = 2000, d.frame_every = 0;
  System w = s;
  run_dynamics(w, d);
  std::vector<std::vector<Vec3>> frames;
  d.steps = 6000, d.frame_every = 30;
  d.frame = [&](const std::vector<double>& x, const Cell&, int64_t step) {
    if (step == 0) return;
    std::vector<Vec3> p(x.size() / 3);
    for (size_t i = 0; i < p.size(); ++i) p[i] = {x[3 * i], x[3 * i + 1], x[3 * i + 2]};
    frames.push_back(p);
  };
  run_dynamics(w, d);
  const double s6 = std::pow(sigma, 6), vc = 4 * eps * (s6 * s6 / std::pow(rc, 12) - s6 / std::pow(rc, 6)), kT = 0.0019872043 * T;
  std::mt19937 rng(3);
  std::uniform_real_distribution<double> U(0, L);
  double wsum = 0;
  long nins = 0;
  for (const auto& f : frames)
    for (int t = 0; t < 400; ++t) {
      const Vec3 p{U(rng), U(rng), U(rng)};
      double du = 0;
      for (const auto& q : f) {
        Vec3 dd = q - p;
        for (int c = 0; c < 3; ++c) dd[c] -= L * std::round(dd[c] / L);
        const double r2 = dd[0] * dd[0] + dd[1] * dd[1] + dd[2] * dd[2];
        if (r2 < rc * rc) { const double x6 = s6 / (r2 * r2 * r2); du += 4 * eps * (x6 * x6 - x6) - vc; }
      }
      wsum += std::exp(-du / kT), ++nins;
    }
  const double mu_widom = -kT * std::log(wsum / double(nins));
  // TI: one particle decoupled
  SolvationOptions o;
  o.field = ff;
  o.energy.cutoff = rc, o.energy.coulomb = false;
  o.solute.assign(s.atoms.size(), 0);
  o.solute[0] = 1;
  o.temperature = T, o.dt = 5, o.tau_t = 100;
  o.ps = 10, o.equilibrate_ps = 2;
  const auto r = solvation_free_energy(w, o);
  std::printf("mu_ex: Widom %.3f kcal/mol, TI %.3f +- %.3f (kT %.3f)\n", mu_widom, r.dg, r.dg_err, kT);
  EXPECT_NEAR(r.dg, mu_widom, std::max(0.3 * kT, 3 * r.dg_err));   // within three standard errors of TI (short windows)
}
