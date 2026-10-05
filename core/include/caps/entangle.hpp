// CAPS entanglements: primitive-path analysis (PPA) after Everaers, Sukumaran, Grest, Svaneborg, Sivasubramanian &
// Kremer, Science 303, 823 (2004). The backbone atoms of every chain become beads; chain ends are held fixed, bonds pull
// with zero rest length, and beads of different chains repel within a diameter σ, so chains cannot pass through each
// other. Interactions within a chain are switched off. Minimising the energy shrinks each backbone to its primitive path;
// the contour lengths give the entanglement length:
//
//   classical coil      N_e = N_b ⟨R²⟩ / ⟨L_pp⟩²                       (Everaers et al. 2004)
//   modified S-coil     N_e = N_b / (⟨L_pp²⟩ / ⟨R²⟩ − 1)               (Hoy, Foteinopoulou & Kröger, PRE 80, 031803, 2009)
//   tube step           a_pp = ⟨R²⟩ / ⟨L_pp⟩
//
// N_b is backbone bonds per chain. This is not Z1 (Kröger's geometric shortest-path reduction) and it counts no kinks.
// The paths lose self-entanglements (a chain passes through itself and its own periodic images).
#pragma once
#include <cstdint>
#include <functional>
#include <string>
#include <vector>

#include "caps/system.hpp"

namespace caps {

struct PrimitivePathOptions {
  double sigma = 0;           // bead diameter σ, Å; 0: the mean backbone bond / 0.97 (the Kremer–Grest bond)
  int max_steps = 200000;     // FIRE steps at most
  double dt = 0.01;           // initial FIRE time step, τ (at most 10×)
  bool per_chain = true;      // FIRE per chain (a chain that overshoots does not reset the others)
  double force_tol = 1e-4;    // stop when the largest force is below this × k σ
  double repulsion = 1;       // ε of the WCA repulsion between chains (Everaers et al.: 1, with K = 30 ε/σ²)
  int record_every = 50;       // steps between the recorded mean L_pp
  std::function<bool(double fraction)> progress;   // return false to stop early
};

struct PrimitivePaths {
  std::vector<std::vector<Vec3>> paths;     // per chain: bead positions (unwrapped, ends as in the frame)
  std::vector<double> lpp, r2, bonds;       // per chain: contour length, squared end-to-end distance (Å²), backbone bonds
  std::vector<double> trace_step, trace_lpp;   // mean L_pp during the minimisation
  double sigma = 0;
  int steps = 0;
  bool converged = false;
  bool stopped = false;
};

// frame must hold whole molecules (make_molecules_whole); backbones from backbones(frame).
PrimitivePaths primitive_paths(const System& frame, const std::vector<std::vector<uint32_t>>& backbones,
                               const PrimitivePathOptions& o = {});

struct EntanglementEstimate {
  double nb = 0;              // mean backbone bonds per chain
  double r2 = 0, lpp = 0, lpp2 = 0;
  double ne_coil = 0, ne_mscoil = 0;   // backbone bonds per entanglement (0 when the paths are straight)
  double a_pp = 0;            // Å
  double z = 0;               // entanglements per chain (N_b / N_e, modified S-coil)
  int chains = 0;
};
EntanglementEstimate entanglement_estimate(const PrimitivePaths& p);

}  // namespace caps
