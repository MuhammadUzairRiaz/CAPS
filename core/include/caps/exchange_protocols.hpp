// CAPS Exchange, protocols (stage E6 of docs/studies/vitrimers/caps_module_proposal.md): the structures a healing test
// starts from, and the runs that measure exchange. Every run takes the exchange moves of caps/exchange_moves.hpp (none
// when ProtocolOptions::exchange is off) and hands back the force field the moves left, so runs chain.
//
//   builders  cut_bonds_across (an ideal cut: Stukalin et al. 2013, Ciarella & Ellenbroek 2019), open_gap (two faces
//             apart, a weld: Yang et al. 2016; a gap to close: Chen et al. 2023), remove_sphere (a void or a crack tip),
//             fragment (the network cut at mid-planes into pieces: recycling, Wick et al. 2025)
//   step_strain  a volume-keeping uniaxial step strain, then a hold with exchange: σ(t) and its KWW fit (Yang et al. 2016)
//   tv_scan      NPT at each temperature with exchange: specific volume and its hinge (Tv from the expansion, Perego &
//                Khabaz 2020), and the exchange rate per temperature with its Arrhenius fit
//   heal_run     MD with exchange in chunks: bonds bridging the two sides of a cut and crossing its plane after each
//                (weld and recycle are heal runs on a gap or on fragments)
//   creep        run_creep (caps/mechanics.hpp) takes the same moves through CreepOptions::topology
#pragma once
#include <cstdint>
#include <memory>
#include <vector>

#include "caps/dynamics.hpp"
#include "caps/exchange_analysis.hpp"
#include "caps/exchange_moves.hpp"
#include "caps/mechanics.hpp"

namespace caps {

// ---------------------------------------------------------------- builders

// The builders drop the bonded terms and force field a file carried (they no longer match), so a run retypes.
size_t cut_bonds_across(System& s, int axis, double plane);        // bonds crossing the plane removed; returns how many
void open_gap(System& s, int axis, double plane, double gap);      // atoms above the plane moved by gap Å, the cell lengthened by it
size_t remove_sphere(System& s, const Vec3& centre, double radius);   // atoms within radius removed with their bonds; returns how many
size_t fragment(System& s, const std::vector<int>& axes);          // the bonds crossing the mid-plane of each axis removed

// ---------------------------------------------------------------- protocols

struct ProtocolOptions {
  std::shared_ptr<const ForceField> field;   // the structure's force field (the protocol replaces it as exchange changes it)
  EnergyOptions energy;
  double temperature = 300, dt = 1.0, tau_t = 100, tau_p = 1000, pressure = 1.0;
  bool npt = false;                          // heal runs: hold the pressure (a gap closes)
  bool exchange = true;
  ExchangeOptions moves;
  int64_t exchange_every = 10;
  uint64_t seed = 1;
};

struct StepStrainResult {
  std::vector<double> t_ps, stress;          // MPa, σ_axis − (σ_other1 + σ_other2)/2
  double kww_tau = 0, kww_beta = 0, kww_mean = 0;   // ps (NaN when the decay is too short to fit)
  ExchangeStats stats;
};
StepStrainResult step_strain(System& s, ProtocolOptions& o, int axis, double strain, double hold_ps, int sample_every = 10);

struct TvScanResult {
  std::vector<double> T, specific_volume, accepted_per_ps;
  BilinearFit fit;                           // hinge of the specific volume: Tv (or Tg) where the expansion changes
  ArrheniusFit rate;                         // of the time between accepted moves
};
TvScanResult tv_scan(System& s, ProtocolOptions& o, double t_from, double t_to, double t_step, double ps_per_step);

struct HealResult {
  std::vector<double> t_ps, bridges, crossing;   // after each chunk
  ExchangeStats stats;
};
// side: per atom, where it started (sides_from_positions on the structure before the cut); plane as interface_bonds.
HealResult heal_run(System& s, ProtocolOptions& o, const std::vector<int8_t>& side, int axis, double plane, double ps, int chunks);

}  // namespace caps
