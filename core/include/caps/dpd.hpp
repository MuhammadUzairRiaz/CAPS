// CAPS DPD: dissipative particle dynamics of coarse-grained beads (Hoogerbrugge & Koelman 1992; Groot & Warren,
// J. Chem. Phys. 107, 4423 (1997)) — the mesoscale of polymer blends, block copolymers and solutions: microdomains
// (lamellae, cylinders, spheres), phase separation, the morphology a formulation takes.
//
// Reduced units: r_c = 1, bead mass 1, kT = 1. Beads interact within r_c by a soft conservative force a_ij (1 − r) r̂, a
// dissipative force −γ w(r)² (r̂·v_ij) r̂ and a random force σ w(r) θ_ij r̂ / √dt with w(r) = 1 − r and σ² = 2γkT (the
// fluctuation–dissipation pair: the pair thermostat conserves momentum). Bonded beads are joined by springs −C (r − r₀).
// Groot–Warren's integrator (their modified velocity Verlet, λ = 0.65).
//
// Chemistry enters through a_ij: like beads a_ii = 75 kT / ρ (25 at ρ = 3, water's compressibility), unlike beads
// a_ij = a_ii + χ_ij / 0.286 at ρ = 3 (Groot & Warren's mapping; at another density give a_ij directly).
//
// Reported: kT and the pressure (virial) with their errors, the segregation order parameter ψ (0 mixed, 1 separated:
// the mean |φ_A − φ_B| / (φ_A + φ_B) over cells of about 2 r_c), the structure factor S(q) of the first type's
// density with its peak q* and domain spacing 2π/q*, and frames for the view (bead types as elements, r_c = rc_angstrom).
#pragma once

#include <cstdint>
#include <functional>
#include <map>
#include <string>
#include <vector>

#include "caps/system.hpp"

namespace caps {

struct DpdSpecies {
  std::string name;       // "PS-b-PMMA", "solvent"
  std::string sequence;   // one letter per bead, A–Z: "AAAAABBBBB" (a diblock), "A" (a monomer)
  int count = 0;          // molecules
};

struct DpdOptions {
  std::vector<DpdSpecies> species;
  double density = 3.0;                           // beads per r_c³ (sets the box)
  std::map<std::string, double> chi;              // "AB" → χ_AB (a_ij from Groot–Warren's mapping at ρ = 3)
  std::map<std::string, double> a;                // "AB" → a_ij directly (wins over χ); "AA" → a_ii
  double gamma = 4.5, dt = 0.04;
  double bond_k = 4.0, bond_r0 = 0.0;             // springs −C (r − r₀)
  long steps = 20000, equilibration = 5000;       // the first `equilibration` steps are not sampled
  int frame_every = 500;
  double rc_angstrom = 6.46;                      // r_c in Å for the frames (Groot & Rabone: 3 waters per bead)
  uint64_t seed = 1;
  std::function<bool(long step, double kT)> progress;   // about every 1 %; false stops
};

struct DpdReport {
  int beads = 0, molecules = 0;
  double box = 0;                                  // r_c
  double kT = 0, kT_error = 0, pressure = 0, pressure_error = 0;
  double order = 0;                                // ψ at the end
  std::vector<double> q, sq;                       // S(q) of type A, radially averaged
  double q_peak = 0, spacing = 0;                  // 2π/q* (r_c)
  std::vector<std::pair<double, double>> order_series;   // (step, ψ)
  std::string types;                               // bead letters in use
  Trajectory frames;                               // beads as atoms (A → C, B → O, C → N, D → S, others → P), Å
  double seconds = 0;
  std::vector<std::string> notes;
};

DpdReport run_dpd(const DpdOptions& o);

// a_ij for types i, j under the options (like beads 75/ρ, unlike from χ at ρ = 3, or given)
double dpd_repulsion(const DpdOptions& o, char i, char j);

}  // namespace caps
