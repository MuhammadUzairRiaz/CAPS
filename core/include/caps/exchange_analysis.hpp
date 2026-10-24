// CAPS Exchange, analysis: what is measured on networks whose bonds exchange (vitrimers, dynamic covalent and reversible
// networks) and on their healing. Stage E0 of docs/studies/vitrimers/caps_module_proposal.md: these functions work on
// data from any run (CAPS's own, LAMMPS fix bond/react or bond/swap, an experiment's table), before CAPS's own MD can
// exchange bonds.
//
//   fit_arrhenius        ln τ = ln τ0 + Ea / (R T), least squares in 1/T with standard errors; Ea in kJ/mol, or in the
//                        temperature's own unit (ε for T in ε/k_B, reduced Lennard-Jones runs)
//   crossing_temperature where two relaxation times measured at the same temperatures cross (the exchange time
//                        overtaking the segmental time: Branham-Ferrari et al., Nat. Phys. 2026, 10.1038/s41567-026-03468-0)
//   stress_relaxation    G(t) = V/(10 k_B T) Σ_αβ ⟨P°_αβ(0) P°_αβ(t)⟩ from pressure-tensor samples (the traceless symmetric
//                        tensor, as the Green–Kubo viscosity), correlated by a multi-tau correlator (Ramírez et al., J. Chem.
//                        Phys. 133, 154103, 2010) so lags span many decades at a fixed cost; G'(ω), G''(ω) by the exact
//                        transform of the piecewise-linear G(t); a KWW fit of the decay after the fast vibrational part
//   bond_memory          the fraction of exchangeable sites still bonded to their time-zero partner, all time origins; a
//                        KWW fit, the 1/e time, and the two-step form f = ½[1 + e^(−t/τ)] exp[−(t/τ_d)^β] of Branham-Ferrari
//                        et al. (Eq. 3), where bonds come back to old partners before diffusion carries them apart
//   interface_bonds      bonds joining atoms that started on different sides of a cut, and the number of bonds crossing
//                        each plane along an axis, B(z) (Ciarella & Ellenbroek, Coatings 9, 114, 2019; Stukalin et al.,
//                        Macromolecules 46, 7525, 2013)
//   healing_efficiency   a healed or welded sample's stress–strain curve against the intact one's: toughness (area) ratio
//                        (Chen et al., Polymer 284, 126284, 2023), peak-stress ratio, modulus ratio (Yuan et al.,
//                        Materials 19, 3702, 2026) and the stress at a set strain (Wick et al., J. Polym. Sci. 63, 9, 2025).
//                        The definitions disagree, so each comes back named.
#pragma once
#include <array>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

#include "caps/properties.hpp"
#include "caps/system.hpp"

namespace caps {

// ---------------------------------------------------------------- Arrhenius

struct ArrheniusFit {
  bool ok = false;
  bool reduced = false;          // Ea in the temperature's unit (k_B = 1), not kJ/mol
  double ea = std::numeric_limits<double>::quiet_NaN(), ea_err = std::numeric_limits<double>::quiet_NaN();
  double ln_tau0 = std::numeric_limits<double>::quiet_NaN(), ln_tau0_err = std::numeric_limits<double>::quiet_NaN();
  double r2 = std::numeric_limits<double>::quiet_NaN();
  size_t points = 0;             // the points fitted (inside t_min..t_max, τ > 0)
  std::vector<double> T, tau;    // those points
  std::string note;
};
// τ against T (K, or ε/k_B with reduced). Points outside [t_min, t_max] (0: no bound) or with τ ≤ 0 are left out; fewer
// than two points, or one temperature, give ok = false and a note.
ArrheniusFit fit_arrhenius(const std::vector<double>& T, const std::vector<double>& tau, bool reduced = false, double t_min = 0,
                           double t_max = 0);
// The temperature where ln τ_a − ln τ_b changes sign between neighbouring temperatures (linear in 1/T between them); NaN
// when it never does. T, tau_a and tau_b have one entry per temperature, in any order.
double crossing_temperature(const std::vector<double>& T, const std::vector<double>& tau_a, const std::vector<double>& tau_b);
std::vector<Property> arrhenius_properties(const ArrheniusFit& f, const std::string& what = "τ");

// ---------------------------------------------------------------- stress relaxation

// Multi-tau correlator: lags j·m^k·dt for level k (level 0: j = 0 … p−1; above: j = p/m … p−1), each level fed the
// averages of m samples of the level below. Several channels at once; result() gives Σ_channels ⟨x(0) x(t)⟩.
class MultiTauCorrelator {
 public:
  MultiTauCorrelator(size_t channels, int p = 16, int m = 2, int levels = 0);   // levels 0: 40
  void add(const std::vector<double>& x);
  // lags (in samples) and the correlation at each lag that has at least one product
  void result(std::vector<double>& lag, std::vector<double>& corr) const;
 private:
  void add_at(int level, const std::vector<double>& x);
  size_t nc_;
  int p_, m_, nlev_;
  std::vector<std::vector<std::vector<double>>> buf_;   // [level][slot][channel], slot 0 the newest
  std::vector<int> filled_;                             // samples held per level (≤ p)
  std::vector<std::vector<double>> sum_, count_;        // [level][lag]
  std::vector<std::vector<double>> acc_;                // [level][channel]
  std::vector<int> nacc_;
};

struct StressRelaxationOptions {
  double dt_fs = 1.0;            // time between samples
  double volume = 0;             // Å³
  double temperature = 300;      // K
  int p = 16, m = 2;             // correlator
  bool subtract_mean = true;     // correlate fluctuations about each component's mean (a solid keeps a mean stress)
  double fit_from_ps = 1.0;      // the KWW fit starts here (after the vibrational drop)
  int omega_points = 60;         // G'(ω), G''(ω) on a log grid from 1/t_max to 1/t_1
  double g_inf = 0;              // MPa, G beyond the last lag (0: a liquid; the equilibrium modulus of a solid network)
  double max_lag_fraction = 0.1; // lags up to this fraction of the run (the longest lags of a multi-tau correlator rest on
                                 // a handful of products)
  double integrate_ps = 0;       // the viscosity integral runs to here; 0: to the first lag where G falls to zero or below
};
struct StressRelaxation {
  std::vector<double> t_ps, g;                 // MPa
  std::vector<double> omega, storage, loss;    // rad/ps, MPa
  double g0 = 0;                               // G(0), MPa
  double viscosity = 0;                        // ∫G dt, mPa·s (the Green–Kubo integral, to integrated_to)
  double integrated_to = 0;                    // ps
  double kww_tau = std::numeric_limits<double>::quiet_NaN(), kww_beta = std::numeric_limits<double>::quiet_NaN(),
         kww_mean = std::numeric_limits<double>::quiet_NaN();   // ps; mean = τ Γ(1/β)/β
  size_t samples = 0;
  std::vector<std::string> notes;
};
// p: pressure tensor samples, atm, in the order xx yy zz xy xz yz.
StressRelaxation stress_relaxation(const std::vector<std::array<double, 6>>& p, const StressRelaxationOptions& o);
// G*(ω) = iω ∫ G(t) e^(−iωt) dt with G linear between the given points and g_inf after the last one; t ascending.
void dynamic_moduli(const std::vector<double>& t, const std::vector<double>& g, const std::vector<double>& omega, double g_inf,
                    std::vector<double>& storage, std::vector<double>& loss);
std::vector<Property> stress_relaxation_properties(const StressRelaxation& r);

// KWW fit y/y0 = exp(−(t/τ)^β) by the straight line ln(−ln C) against ln t, over the points with lo < C < hi (C = y/y0, y0
// the first point at or after t_from). Returns {τ, β, ⟨τ⟩ = τ Γ(1/β)/β}; NaNs when there are too few points.
std::array<double, 3> fit_kww(const std::vector<double>& t, const std::vector<double>& y, double t_from = 0, double lo = 0.05,
                              double hi = 0.95);

// ---------------------------------------------------------------- bond memory

struct TwoStepFit {
  bool ok = false;
  double tau = std::numeric_limits<double>::quiet_NaN();       // the first, exponential step (returns to the old partner)
  double tau_d = std::numeric_limits<double>::quiet_NaN();     // the stretched decay (partners carried apart)
  double beta = std::numeric_limits<double>::quiet_NaN();
  double rms = std::numeric_limits<double>::quiet_NaN();
};
// f(t) = ½[1 + e^(−t/τ)] exp[−(t/τ_d)^β] by least squares (Nelder–Mead in ln τ, ln τ_d, β; 0.1 ≤ β ≤ 2).
TwoStepFit fit_two_step(const std::vector<double>& t, const std::vector<double>& f);

struct BondMemory {
  std::vector<double> t, c;      // c(t): bonded at the origin and to the same partner at t, over bonded at the origin
  double tau_e = std::numeric_limits<double>::quiet_NaN();      // c = 1/e (linear between lags)
  double kww_tau = std::numeric_limits<double>::quiet_NaN(), kww_beta = std::numeric_limits<double>::quiet_NaN(),
         kww_mean = std::numeric_limits<double>::quiet_NaN();
  TwoStepFit two_step;
  size_t sites = 0, frames = 0;
  size_t changes = 0;            // partner changes between consecutive frames, all sites
  double change_time = std::numeric_limits<double>::quiet_NaN();  // site-frames bonded / changes × frame time: the mean
                                 // time between partner changes seen at this frame spacing (events between frames are missed)
  std::vector<std::string> notes;
};
// partners[f][i]: the partner of exchangeable site i in frame f (any id), −1 when free; every frame lists the same sites.
// frame_dt: the time between frames (any unit; the results carry it).
BondMemory bond_memory(const std::vector<std::vector<int64_t>>& partners, double frame_dt);
std::vector<Property> bond_memory_properties(const BondMemory& r, const std::string& time_unit = "ps");

// ---------------------------------------------------------------- healing

struct InterfaceBonds {
  size_t bonds = 0;              // bonds in the structure
  size_t bridges = 0;            // bonds between an atom of side 0 and one of side 1
  size_t crossing_plane = 0;     // bonds crossing the plane at `plane` (fraction of the cell along the axis)
  std::vector<double> z, crossings;   // B(z): bonds crossing the plane at each bin edge (Å along the axis)
};
// side[i] ∈ {0, 1} for each atom (where it started), or −1 to leave the atom out; bonds by the minimum image.
InterfaceBonds interface_bonds(const System& s, const std::vector<int8_t>& side, int axis = 2, double plane = 0.5, double bin = 1.0);
// The side of each atom from its position: 0 below `plane` (fraction of the cell along the axis), 1 above.
std::vector<int8_t> sides_from_positions(const System& s, int axis = 2, double plane = 0.5);
std::vector<Property> interface_properties(const InterfaceBonds& r, int axis);

struct StressStrain {
  std::vector<double> strain, stress;   // strain as a fraction; stress in any unit, the same for both curves
};
struct HealingEfficiency {
  double toughness = std::numeric_limits<double>::quiet_NaN();   // ∫σ dε healed / intact, %
  double peak = std::numeric_limits<double>::quiet_NaN();        // σ_max healed / intact, %
  double modulus = std::numeric_limits<double>::quiet_NaN();     // initial slope healed / intact, %
  double stress_at = std::numeric_limits<double>::quiet_NaN();   // σ(ε*) healed / intact, % (NaN when no ε* given)
  double up_to = 0;              // strain the areas run to (0: each curve to its own end)
  double modulus_strain = 0, at_strain = 0;
  std::vector<std::string> notes;
};
// modulus_strain: the initial slope is a straight-line fit (with intercept) of the points with 0 ≤ ε ≤ modulus_strain;
// at_strain > 0: the stress there (linear between points, averaged over ±at_window); up_to > 0: both areas to that strain.
HealingEfficiency healing_efficiency(const StressStrain& healed, const StressStrain& intact, double modulus_strain = 0.03,
                                     double at_strain = 0, double up_to = 0, double at_window = 0);
std::vector<Property> healing_properties(const HealingEfficiency& h);

}  // namespace caps
