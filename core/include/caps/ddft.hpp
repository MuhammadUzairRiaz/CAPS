// CAPS density-field dynamics: the mesoscale morphology of polymer blends and block copolymers as fields on a grid —
// dynamic mean-field density functional theory of Gaussian chains (Fraaije, J. Chem. Phys. 99, 9202 (1993); Fraaije et
// al., J. Chem. Phys. 106, 4260 (1997); Maurits & Fraaije, J. Chem. Phys. 107, 5879 (1997)).
//
// Each bead type I is a volume-fraction field θ_I(r) on a periodic grid (spacing h = 1, the unit of length). Molecules are
// Gaussian chains of beads (a sequence such as A8B8), bond length a: the bond kernel exp(−q²a²/6) in Fourier space.
//   Chain functional  an external potential U_I (kT) per type gives the densities by the chain propagators
//                     q_1 = w_1, q_s = w_s (K ⋆ q_{s−1}), q†_N = w_N, q†_s = w_s (K ⋆ q†_{s+1}), w_s = e^(−U_{type(s)}),
//                     θ_s ∝ q_s q†_s / w_s, normalised to the molecule's share of the beads.
//   Free energy       F = F_ideal[θ] + Σ_{I<J} χ_IJ ∫ θ_I (K ⋆ θ_J) + (κ_H / 2) ∫ (Σ_I θ_I − 1)² (kT per bead volume),
//                     μ_I = −U_I + Σ_J χ_IJ (K ⋆ θ_J) + κ_H (Σ θ − 1): U_I found for the present θ by inverting the chain
//                     functional (U_I ← U_I + λ ln(θ_I[U] / θ_I)).
//   Dynamics          ∂θ_I/∂τ = ∇·(θ_I ∇μ_I) + ∇·ζ_I, τ = M kT t / h² (local mobility), finite volumes on the cell faces
//                     (mass conserved exactly), conserved Gaussian noise ⟨ζζ⟩ = (2 θ / Ω) δ (Ω the noise parameter: beads
//                     per cell, larger is quieter; 0 none).
// A grid of n_z = 1 is a two-dimensional system.
#pragma once

#include <array>
#include <cstdint>
#include <functional>
#include <map>
#include <string>
#include <vector>

namespace caps {

struct DdftMolecule {
  std::string name;
  std::string sequence;   // bead types as letters, runs allowed: "A8B8", "AAAABBBB", "A10"
  double fraction = 1;    // volume fraction of all beads (normalised over the molecules)
};

struct DdftOptions {
  std::vector<DdftMolecule> molecules;
  std::map<std::string, double> chi;   // "AB": χ between bead types A and B (others 0)
  std::array<int, 3> grid{32, 32, 32};
  double bond = 1.1543;                // a, grid spacings
  double kappa = 10;                   // Helfand compressibility, kT
  double noise = 75;                   // Ω (0: none)
  double dt = 0.02;                    // Δτ
  int steps = 2000;
  int frame_every = 100;
  int inverse_iterations = 60;
  double inverse_tol = 1e-4;           // largest |θ[U] − θ| allowed
  uint64_t seed = 1;
  int threads = 0;
  std::function<bool(double)> progress;   // fraction; false stops
};

struct DdftResult {
  std::vector<std::string> types;      // bead types, in order
  std::array<int, 3> grid{0, 0, 0};
  std::vector<double> times;           // τ of each frame
  std::vector<std::vector<float>> frames;   // per frame: θ of every type, type-major, then x fastest (i + nx (j + ny k))
  std::vector<std::vector<double>> order;   // per frame, per type: ⟨(θ − θ̄)²⟩
  std::vector<double> interaction;     // per frame: Σ χ_IJ ⟨θ_I (K ⋆ θ_J)⟩ (kT per bead volume)
  std::vector<double> compressibility; // per frame: (κ_H / 2) ⟨(Σθ − 1)²⟩
  std::vector<double> dominant_spacing;   // per frame: 2π / q of the strongest |θ̂|² of type 0, grid spacings
  int worst_inverse = 0;               // most inverse iterations a step needed
  double seconds = 0;
  std::vector<std::string> notes;
};

DdftResult ddft_run(const DdftOptions& o);

// The ideal-chain functional alone (for tests and analysis): the volume fractions of each type for potentials U (type-major,
// as the frames), molecules and bond length as in o.
std::vector<double> ddft_chain_densities(const DdftOptions& o, const std::vector<double>& U);

}  // namespace caps
