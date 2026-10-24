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
#include <vector>

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

// Bond exchange in LAMMPS (fix bond/swap of the MC package; Sides et al., J. Polym. Sci. B 42, 199 (2004)) for a
// Kremer–Grest system that is already equilibrated: the structure's own bonds and positions, FENE + WCA (+ cosine bending
// on the angles of the bond graph), special_bonds fene (the 0 1 1 the fix needs), Langevin at T = 1 and the swaps. The
// atoms that may swap get atom type 2 (the same bead) and form the group "exch". Molecule IDs decide which swaps LAMMPS
// allows, so they are written for one of its two uses:
//   chains (stickers = false): every atom may swap and the molecule ID is the bead's place along its chain (1 … N, or
//     1 … N/2 … 1 with from_both_ends), so swaps keep chain lengths; CAPS's ExchangeEngine takes the same numbers as its
//     labels (chain_position_labels)
//   stickers: the group's atoms share molecule ID 1 (any two stickers may swap), the others 2; bonds between two stickers
//     carry no angle (the manual asks for that) and the bonds are stored on both atoms (newton on off)
// The fix's two outputs (swaps done, swaps tried) go into every thermo line.
struct BondSwapExport {
  std::vector<uint32_t> group;     // the atoms that may swap (empty: every atom)
  bool stickers = false;
  bool from_both_ends = false;
  int every = 50;                  // steps between tries (Nevery)
  double fraction = 0.5;           // fraction of the group's atoms considered at each try
  double cutoff = 1.3;             // σ: every one of the four distances must be shorter
  uint64_t seed = 1;
  double run_steps = 100000;
};
void write_kg_bond_swap(const System& s, const KgOptions& o, const BondSwapExport& e, const std::string& stem);

}  // namespace caps
