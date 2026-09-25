#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <numeric>

#include "caps/dynamics.hpp"
#include "caps/grow.hpp"
#include "caps/io.hpp"
#include "caps/relax.hpp"

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
