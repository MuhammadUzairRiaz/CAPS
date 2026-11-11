// CAPS melting point by the Z method.
//
//  Runs       one NVE run per starting temperature T0, each from the perfect crystal (the structure as given, at rest)
//             with Maxwell–Boltzmann velocities at T0, at the crystal's volume. Its total energy E is fixed by T0; its
//             temperature T is the mean over the second half of the run (a harmonic crystal settles at T0 / 2).
//  Order      the crystal's own Bragg order ρ(G) = |Σ_j exp(i G·r_j)| / N over the star of the reciprocal vectors G
//             (d ≥ 2 Å) where it is largest in the perfect cell, as a fraction of that start value: near 1 − G²⟨u²⟩/3 in
//             a crystal, ~ N^(−1/2) in a liquid (Verlet's translational order parameter). A run whose order over its last
//             fifth stays below `melt_order` has melted.
//  Z          T(E) rises along the crystal branch up to the critical superheating point (E_LS, T_LS); the next energy
//             that melts drops onto the liquid at the melting temperature: T_m is the temperature of the lowest-energy
//             melted run (Belonoshko, Skorodumova, Rosengren & Johansson, Phys. Rev. B 73, 012201 (2006)). Bisection
//             runs between the last crystal and the first melted start narrow the step. The waiting time for melting
//             grows sharply near E_LS: short runs overestimate the superheating point but not the bottom of the Z
//             unless the melted run nearest it is missing.
#pragma once

#include <array>
#include <cstdint>
#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "caps/field.hpp"
#include "caps/json.hpp"
#include "caps/system.hpp"

namespace caps {

struct ZMethodOptions {
  std::shared_ptr<const ForceField> field;   // null: CAPS's built-in typing
  EnergyOptions energy;
  std::vector<double> temperatures;   // starting temperatures T0 (K), any order
  int64_t steps = 10000;
  double dt = 1.0;                    // fs
  int sample_every = 20;              // steps between samples of T, E, P and the order
  int bisect = 0;                     // runs added between the last crystal and the first melted start
  double melt_order = 0.3;            // order fraction below which a run has melted
  uint64_t seed = 1;
  std::function<bool(const std::string&, double)> progress;
};

struct ZPoint {
  double t0 = 0;                      // K, starting
  double energy = 0;                  // kcal/mol per atom, total (kinetic + potential), mean over the run
  double temperature = 0, temperature_sd = 0;   // K, second half
  double pressure = 0;                // atm, second half
  double order = 0;                   // last fifth, fraction of the start value
  double msd = 0;                     // Å², from the lattice sites at the end
  double drift = 0;                   // kcal/mol per atom, last minus first total energy
  bool melted = false;
  std::vector<double> time_ps, t_trace, order_trace;
};

struct ZResult {
  std::vector<ZPoint> points;         // by increasing energy
  double t_melt = 0;                  // K (0: no Z found, i.e. no melted run colder than the hottest crystal)
  double t_superheat = 0, e_superheat = 0;   // the hottest crystal run below the first melted energy
  double e_melt = 0;                  // energy of the lowest melted run (0: no Z found)
  std::vector<std::array<double, 3>> g;   // the reciprocal vectors of the order parameter, 1/Å
  double order0 = 0;                  // ρ(G) of the perfect cell
  std::vector<std::string> notes;
};

ZResult z_method(const System& crystal, const ZMethodOptions& o);
Json z_method_json(const ZResult& r);

}  // namespace caps
