// CAPS Dynamics: molecular dynamics with the Field evaluator (units real: fs, Å, kcal/mol, K, atm).
#pragma once
#include <cstdint>
#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "caps/field.hpp"
#include "caps/system.hpp"

namespace caps {

enum class Thermostat { None, Bussi, Langevin };
enum class Barostat { None, CRescale, Berendsen };

Thermostat thermostat_from_string(const std::string& s);   // "none" | "bussi" | "langevin"
Barostat barostat_from_string(const std::string& s);       // "none" | "crescale" | "berendsen"
const char* to_string(Thermostat t);
const char* to_string(Barostat b);

struct ThermoRow {
  int64_t step = 0;
  double time_ps = 0;
  double temperature = 0;      // K, 3N − 3 degrees of freedom
  double potential = 0, kinetic = 0, total = 0;   // kcal/mol
  double conserved = 0;        // total plus the energy the thermostat and barostat exchanged (drift check)
  double pressure = 0;         // atm, virial + kinetic
  double target_temperature = 0;   // K, the thermostat target at this step
  double volume = 0, density = 0;
  double p[6] = {0, 0, 0, 0, 0, 0};   // pressure tensor, atm (kinetic + virial): xx yy zz xy xz yz; stress = −p
  double lx = 0, ly = 0, lz = 0;      // cell edge lengths, Å
  double pull_force = 0;              // spring force along the pull direction, kcal/mol/Å
  double pull_disp = 0;               // displacement of the pulled group's centre along it, Å
};

struct DynamicsOptions {
  // The force field to use (from CAPS Field); null: CAPS's built-in GAFF typing of C and H.
  std::shared_ptr<const ForceField> field;
  double dt = 1.0;                  // fs
  // r-RESPA (Tuckerman, Berne & Martyna, J. Chem. Phys. 97, 1990 (1992)): > 1 splits each step: non-bonded forces at dt,
  // bonded forces at dt / respa. NVE or the Bussi thermostat; barostats act on the outer step.
  int respa = 1;
  int64_t steps = 10000;
  double temperature = 300.0;       // K, thermostat target and initial velocities
  double temperature_end = -1.0;    // K; ≥ 0 ramps the thermostat target linearly to this over the run
  Thermostat thermostat = Thermostat::Bussi;
  double tau_t = 100.0;             // fs; Bussi relaxation time, Langevin 1/γ
  Barostat barostat = Barostat::None;
  double pressure = 1.0;            // atm
  double tau_p = 1000.0;            // fs
  double compressibility = 4.5e-5;  // atm⁻¹, isothermal (sets the barostat's response, not the result)
  int barostat_every = 10;          // steps between cell updates
  // Pressure coupling: false: isotropic (one scale for the cell). true: each Cartesian axis with couple_axis[k] set is
  // scaled on its own from the diagonal pressure component P_kk (Berendsen only), the others are left alone — for
  // uniaxial tests, where the lateral axes follow the target pressure while the pulled axis is deformed.
  bool anisotropic = false;
  bool couple_axis[3] = {true, true, true};
  // Deformation at a constant engineering strain rate (as LAMMPS fix deform erate, remap x): the cell length along
  // deform_axis (0 x, 1 y, 2 z) is L0 (1 + deform_rate t), with t in ps from the start of this run.
  int deform_axis = -1;
  double deform_rate = 0.0;         // 1/ps
  bool new_velocities = false;      // draw Maxwell–Boltzmann velocities even when the system has some
  // Per atom: 1 = held in place (a substrate under a film): no velocity, no force, not counted in the temperature.
  std::vector<char> fixed;
  // Steered pulling (pull-out, debonding): a spring of stiffness pull_k on the centre of mass of the atoms in pull_group,
  // its anchor moving along pull_dir at pull_rate from where the centre starts; the spring force is shared over the
  // group by mass. The force and the centre's displacement along pull_dir go into each thermo row.
  std::vector<char> pull_group;
  Vec3 pull_dir{1, 0, 0};
  double pull_k = 10.0;             // kcal/mol/Å²
  double pull_rate = 0.0;           // Å/ps
  uint64_t seed = 1;
  int thermo_every = 100;           // steps between thermo rows
  int64_t step_offset = 0;          // added to reported steps and times (runs chained into a protocol)
  int frame_every = 1000;           // steps between recorded frames (0 = none)
  EnergyOptions energy;
  // Called with each thermo row; return false to cancel (throws DynamicsCancelled).
  std::function<bool(const ThermoRow&)> progress;
  // Called after every force evaluation of the trajectory (each MD step) with the energy terms (virial tensor
  // included), positions and cell: for on-the-fly averages such as stress fluctuations.
  std::function<void(const EnergyTerms&, const std::vector<double>&, const Cell&, int64_t)> each_step;
  // Called for each recorded frame with positions (3N, unwrapped), cell and step.
  std::function<void(const std::vector<double>&, const Cell&, int64_t)> frame;
};

struct DynamicsReport {
  std::vector<ThermoRow> thermo;
  int64_t steps = 0;
  double seconds = 0, ns_per_day = 0;
  int list_builds = 0;
  std::vector<std::string> notes;
};

struct DynamicsCancelled : std::runtime_error {
  DynamicsCancelled() : std::runtime_error("dynamics cancelled") {}
};

// Runs MD on `s` in place: positions (unwrapped), cell, velocities and charges are updated. Velocities are taken from
// the system when present (unless new_velocities), otherwise drawn from Maxwell–Boltzmann at the target temperature.
void run_dynamics(System& s, const DynamicsOptions& o, DynamicsReport* report = nullptr);

// Kinetic energy (kcal/mol) of velocities v (3N, Å/fs) with masses m (g/mol).
double kinetic_energy(const std::vector<double>& v, const std::vector<double>& m);

}  // namespace caps
