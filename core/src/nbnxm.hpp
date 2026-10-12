// Cluster pair list and kernel for plain pair terms in mixed precision (the scheme of GROMACS's nbnxm: Páll & Hess,
// Comput. Phys. Commun. 184, 2641 (2013)). Atoms are put in a grid of columns, sorted along each column and cut into
// clusters of four; the list holds pairs of clusters, each with its periodic image and a 16-bit mask of the atom pairs
// that interact (exclusions, padding, the half list, the list radius). The kernel takes four j atoms per SIMD register
// against each i atom in turn, in single precision; energies, virial and forces leave it in double.
#pragma once

#include "caps/system.hpp"
#include <cstdint>
#include <vector>

namespace caps {

class ThreadPool;

namespace nbnxm {

constexpr uint32_t kPad = 0xffffffffu;   // a slot with no atom

struct Entry {
  uint32_t cj;      // j cluster
  uint16_t code;    // periodic image of the j cluster
  uint16_t mask;    // bit 4a + l: i atom a with j atom l interact
  uint16_t ex;      // the same bits for excluded pairs (1-2, 1-3, 1-4 at their nearest image) when the kernel takes
                    // their Ewald / DSF correction (List::build with p14), else 0
  uint16_t ex14;    // those of them whose Coulomb keeps a fraction (1-4 pairs)
};

struct List {
  // build inputs
  double rlist = 0;
  // slots: cluster c holds slots 4c … 4c + 3; slot s holds atom perm[s] (kPad for padding)
  size_t nclusters = 0;
  std::vector<uint32_t> perm, iperm;
  std::vector<double> off;        // 3 per slot: lattice offset to the binned position, Cartesian, for cells
  std::vector<int32_t> ioff;      // 3 per slot: the same as integer lattice coefficients
  std::vector<int32_t> qint;      // 3 per image code
  std::vector<double> qvec;       // 3 per image code, Cartesian, for cells
  std::vector<size_t> first;      // per i cluster: its entries [first[c], first[c + 1])
  std::vector<Entry> entries;
  // Work in chunks of i clusters [tb[t], te[t]) with about equal entries (4 per worker), and dynamic pruning (as GROMACS):
  // chunk t's entries with only the atom pairs inside a shorter radius at the last prune, i cluster tb[t] + k at
  // [pf[t][k], pf[t][k + 1])
  std::vector<size_t> tb, te;
  std::vector<std::vector<Entry>> pe;
  std::vector<std::vector<uint32_t>> pf;
  // per step, by slot: positions relative to their cluster's origin (3 per cluster, double), so single precision keeps
  // its digits for the separations; an entry's shift is the image plus the difference of the origins
  std::vector<double> org;
  std::vector<float> X, Y, Z, Q, P12, P6;
  std::vector<uint32_t> T;        // type index (0 for padding)
  // When the coefficients are not products of per-type factors: for i type t and j cluster c the four j atoms' c12
  // and c6 side by side, CT12[4 (t·nclusters + c) + l] (a row loads them at once instead of gathering four); empty
  // when that would take more than 64 MB (the kernel then gathers)
  std::vector<float> CT12, CT6;
  Cell cells;

  // the excluded partners (sorted), types, the type pairs that never interact (empty: none), charges, and per-type
  // factors whose products are the pair coefficients (empty: the kernel looks them up per type pair)
  void build(const std::vector<double>& x, const Cell& cell_in, double rlist, const std::vector<std::vector<uint32_t>>& excluded,
             const std::vector<int>& type_index, size_t ntypes, const std::vector<char>& skip_type, const std::vector<double>& q,
             const std::vector<double>& f12, const std::vector<double>& f6, ThreadPool& pool, const float* c12 = nullptr,
             const float* c6 = nullptr, const std::vector<std::vector<uint32_t>>* p14 = nullptr);
  // p14 (per atom, sorted: the excluded partners whose Coulomb keeps a fraction): the entries then carry the excluded
  // pairs too (ex, ex14), for the kernel's corrections; without it excluded pairs are left out
  void rescale(const Cell& cell);                       // offsets and images for a deformed cell
  void positions(const std::vector<double>& x, ThreadPool& pool);   // X, Y, Z for this step
  void prune(double radius, ThreadPool& pool);          // pe, pf from the entries at the positions X, Y, Z
};

struct Params {
  float rc2 = 0, rc = 0;
  int lj = 1;                     // 1: 12-6, 2: class II 9-6
  bool comb = false;              // the coefficients are products of per-type factors (P12, P6 per slot)
  const float* c12 = nullptr;     // per type pair, when not comb
  const float* c6 = nullptr;
  size_t nt = 0;
  float shift = 0;                // 1: energies shifted to zero at the cut-off
  float rcn = 0, rc6 = 0;         // r_c⁻¹² (or r_c⁻⁹ for 9-6) and r_c⁻⁶
  int coul = 0;                   // 0 none, 1 Ewald real space (PME), 2 damped shifted force
  float kappa = 0, kc = 0;        // κ and 2κ/√π
  float e0 = 0, f0 = 0;           // DSF shift constants
  float kcoul = 0;                // Coulomb constant
  bool energies = true;           // false: forces and virial only (the energies left at 0)
  float f14 = 0;                  // the fraction of Coulomb a 1-4 pair keeps (entries' ex14 lanes)
};

// Forces added to fx, fy, fz (by slot, float), energies added, and per image code the summed forces on the j atoms
// taken through that image (fsh, 3 per code), for the virial: Σ d ⊗ f over the pairs = Σ x ⊗ F over the slots (F the
// pair forces) + Σ S ⊗ fsh over the images (pair_virial, shift_virial)
struct Sums {
  double evdw = 0, ecoul = 0;
  std::vector<double> fsh;
};
// Σ x ⊗ F over slots [sb, se) (F by slot, 3 per slot, double; x the listed positions), added to W (row-major 3 × 3)
void pair_virial(const List& L, const double* F, size_t sb, size_t se, double W[9]);
// Σ S ⊗ fsh over the images, added to W
void shift_virial(const List& L, const std::vector<double>& fsh, double W[9]);
// chunk t's i clusters, from its pruned entries, with the widest kernel this processor runs
void kernel(const List& L, const Params& p, int t, float* fx, float* fy, float* fz, Sums& out);
const char* instruction_set();   // the one kernel() uses: "AVX2", "SSE2", "NEON" or "scalar"

// the eight-lane kernel (nbnxm_avx2.cpp, built with AVX2 and FMA on x86-64): false when not built in
bool kernel_avx2(const List& L, const Params& p, int t, float* fx, float* fy, float* fz, Sums& out);
bool avx2_built();

}  // namespace nbnxm
}  // namespace caps
