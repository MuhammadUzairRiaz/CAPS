// CAPS coarse-grained polymer builder (the coarse-graining workflow's "Build CG melt" stage): chains of repeat units, each
// unit a short bead sequence (BS = B S, BA = B A, BT = B T), drawn with real sequence statistics, as random walks whose bond
// lengths, angles and dihedrals are drawn from the inverted distributions (caps/cg_bonded.hpp) — the chains start with the
// local structure of the all-atom model and are pushed apart afterwards (Auhl, Everaers, Grest, Kremer & Plimpton,
// J. Chem. Phys. 119, 12718 (2003)).
//
//   sequence   bernoulli  each unit drawn by the composition
//              markov     first order: the next unit by a transition matrix (rows: the unit before)
//              block      blocks of given lengths, cycling through the units
//              gradient   the first unit's share falls linearly along each chain from 1 to 0 (others share the rest)
//              alternating, or a pattern ("AAB" over the units in order)
//   lengths    every chain dp units, or drawn (schulz-zimm, flory, poisson, log-normal) with a number-average dp and Đ
//   placement  chain starts uniform in the cell, each walk's first bond in a random direction; the cell is cubic, sized
//              by the density; the chains' own overlaps are left for the push-off
//
// The report gives each chain's R_ee, the mean ⟨R_ee²⟩^½, the box edge L and L/⟨R_ee²⟩^½ (a warning when the box is smaller
// than the chains), the units drawn and the internal distances ⟨R²(n)⟩/n of the walks (compare with the all-atom model's:
// cg_internal_distances on its mapped frames).
#pragma once
#include <map>
#include <string>
#include <vector>

#include "caps/cg_bonded.hpp"
#include "caps/cg_rules.hpp"

namespace caps {

struct CgUnit { std::string name; std::vector<std::string> beads; };

struct CgBuildOptions {
  std::vector<CgUnit> units;
  std::map<std::string, double> composition;                          // unit → share (normalised)
  std::string sequence = "bernoulli";                                 // bernoulli | markov | block | gradient | alternating | pattern
  std::map<std::string, std::map<std::string, double>> markov;        // markov: unit before → (next unit → probability)
  std::vector<int> blocks;                                             // block: lengths
  std::string pattern;                                                 // pattern: letters A, B, … for the units in order
  int dp = 100, chains = 10;
  std::string lengths = "monodisperse";                                // or a distribution of caps/polymer.hpp
  double pdi = 1.0;
  double density = 1.2;                                                // g/cm³
  uint64_t seed = 1;
};

struct CgBuildResult {
  System beads;                      // the melt (element 0, named and typed by kind; molecule per chain; cubic cell)
  CgMapping topology;                // bead kinds, molecules, bonds, angles, dihedrals, chains (no atoms)
  std::vector<double> ree;           // per chain, Å
  double ree_rms = 0, contour = 0;   // ⟨R_ee²⟩^½ and the mean contour length (Σ bond lengths), Å
  double box = 0;                    // cubic edge, Å
  std::map<std::string, int> units_drawn;
  std::vector<int> chain_dp;
  std::vector<double> internal;      // ⟨R²(n)⟩/n (Å²) for n = 1, 2, … beads apart
  std::vector<std::string> sequences;   // the first chains as unit names
  std::vector<std::string> notes;
};

// masses: per bead kind (g/mol), e.g. from the maps
CgBuildResult build_cg_polymer(const CgBuildOptions& o, const CgBondedResult& bonded, const std::map<std::string, double>& masses);

// ⟨R²(n)⟩/n (Å²) over chains of a topology and frames: chains made whole along their bonds (minimum image), n up to the
// longest chain − 1; counts give how many pairs went into each n.
std::vector<double> cg_internal_distances(const CgTopology& t, const std::vector<std::vector<Vec3>>& frames, const std::vector<Cell>& cells,
                                          std::vector<double>* counts = nullptr);

// The bead structure and topology as files: STEM.cg.data (LAMMPS, types numbered by `types`) and STEM.map.json (a mapping
// without atoms, for cgfit and the analyses).
std::vector<std::string> write_cg_build(const CgBuildResult& r, const CgTypes& types, const std::string& stem);

// LAMMPS equilibration of a built melt (in.cg_equil): soft push-off (pair_style soft, its prefactor ramped from 0 to
// `soft_max` kcal/mol over `pushoff_steps`, Langevin at the anneal temperature, bonded terms on), then the model's pairs with
// steps still limited (a third as long), a hot anneal (`anneal_T` K, `anneal_steps`), cooling to `T` and NPT at `P` atm; dumps for the checks.
struct CgEquilOptions {
  double T = 300, P = 1, anneal_T = 500, dt = 10, soft_max = 60;
  int64_t pushoff_steps = 300000, anneal_steps = 500000, cool_steps = 200000, npt_steps = 500000, dump_every = 10000;
  int exclude = 3;
};
std::string cg_equil_deck(const CgEquilOptions& o);

}  // namespace caps
