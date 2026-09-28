// CAPS coarse-grained melts (design/boards/CoarseGrained): Kremer–Grest bead-spring chains (J. Chem. Phys. 92, 5057
// (1990)) — M chains of N beads at a reduced density ρσ³ in a cubic box L = (M N / ρ)^{1/3} σ — built as random walks
// with a bond of 0.97 σ that never folds straight back (bond angles above 60°; a bending stiffness k_θ makes straighter
// walks more likely), and written for LAMMPS in reduced units with the standard preparation: a soft-core push-off whose
// strength is ramped up (Auhl et al. 2003), then FENE bonds (K = 30 ε/σ², R₀ = 1.5 σ) and the WCA repulsion (LJ cut at
// 2^{1/6} σ and shifted), Langevin thermostat at T = 1.
//
// In the structure, σ = 1 Å numerically (positions and box in σ); beads carry element 6 for drawing and mass 1.
#pragma once
#include <cstdint>
#include <string>

#include "caps/field.hpp"
#include "caps/system.hpp"

namespace caps {

struct KgOptions {
  int chains = 50;            // M
  int beads = 100;            // N
  double density = 0.85;      // ρσ³
  double k_theta = 0.0;       // bending stiffness, ε (0: flexible)
  double bond = 0.97;         // σ, the FENE–WCA bond length at T = 1
  uint64_t seed = 1;
  // Real units (all three > 0): the reduced model mapped onto a polymer by the user's choice of σ, of the temperature that
  // sets ε = k_B T and of the bead mass m — every length × σ, energy × ε, mass × m, time × τ = σ √(m/ε) — written in LAMMPS
  // units real. CAPS supplies no mapping of its own: the values come from the user (e.g. Everaers et al. 2020's Kuhn-scale
  // mapping of a commodity polymer).
  double sigma = 0;           // Å
  double temperature = 0;     // K
  double bead_mass = 0;       // g/mol
};

// The mapped units of a Kremer–Grest model: ε (kcal/mol), τ (fs), and the mass density (g/cm³) at ρσ³; zero when unmapped.
struct KgUnits { double eps = 0, tau_fs = 0, density = 0; };
KgUnits kg_units(const KgOptions& o);

struct KgReport {
  double box = 0;             // σ
  double closest = 0;         // closest non-bonded pair, σ (overlaps are expected before the push-off)
  double mean_r2 = 0;         // ⟨R_ee²⟩ / (N − 1) of the built chains, σ² per bond
};

System kremer_grest(const KgOptions& o, KgReport* rep = nullptr);

// The Kremer–Grest force field of a melt, for CAPS's own runs and the Field page: one bead type, FENE bonds with their WCA
// core (K = 30 ε/σ², R₀ = 1.5 σ), the WCA pair (LJ cut at 2^{1/6} σ and shifted), the cosine bending k_θ ε (1 + cos θ) when
// k_θ > 0, and special_bonds fene (1-2 pairs out, 1-3 and 1-4 in full). Mapped (sigma, temperature, bead_mass all > 0):
// in Å, kcal/mol and g/mol for a structure in Å; otherwise reduced units read as σ = 1 Å, ε = 1 kcal/mol, m = 1 g/mol
// (then T* = 1 is 503.2 K and τ = 48.9 fs in CAPS's units).
ForceField kremer_grest_forcefield(const System& s, const KgOptions& o);

// LAMMPS data (atom_style angle; bond type 1; angle type 1 when k_θ > 0) and an input deck: push-off with pair soft, then
// the Kremer–Grest force field; `stem` gets .data and .in appended.
void write_kg_lammps(const System& s, const KgOptions& o, const std::string& stem, double pushoff_steps = 20000, double run_steps = 100000);

}  // namespace caps
