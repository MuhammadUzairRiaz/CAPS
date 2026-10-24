#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "caps/dynamics.hpp"
#include "caps/exchange_moves.hpp"
#include "caps/field.hpp"
#include "caps/kremer_grest.hpp"
#include "caps/molecule.hpp"

using namespace caps;

namespace {

constexpr double kB = 0.0019872043;

Cell box(double x, double y, double z) {
  Cell c;
  c.a = {x, 0, 0};
  c.b = {0, y, 0};
  c.c = {0, 0, z};
  return c;
}

std::vector<double> coords(const System& s) {
  std::vector<double> x;
  for (const auto& a : s.atoms) x.insert(x.end(), {a.pos[0], a.pos[1], a.pos[2]});
  return x;
}

double whole_energy(const ForceField& ff, const System& s, const EnergyOptions& eo) {
  Evaluator ev(ff, eo);
  std::vector<double> f;
  return ev.compute(coords(s), s.cell, f).total();
}

// Atoms at fixed places with harmonic bonds (k (r − 1)², kcal/mol/Å²) and no pair terms: every state's energy is known.
struct Frozen {
  System s;
  ForceField ff;
};
Frozen frozen(const std::vector<Vec3>& pos, const std::vector<std::array<uint32_t, 2>>& bonds, double k) {
  Frozen f;
  f.s.cell = box(30, 30, 30);
  for (const auto& p : pos) {
    Atom a;
    a.element = 6;
    a.pos = p + Vec3{10, 10, 10};
    f.s.atoms.push_back(a);
  }
  for (const auto& b : bonds) f.s.bonds.push_back({b[0], b[1], 1});
  const size_t n = pos.size();
  f.ff.name = "frozen";
  f.ff.atom_type.assign(n, "X");
  f.ff.type_index.assign(n, 0);
  f.ff.type_names = {"X"};
  f.ff.lj = {{0.0, 1.0}};
  f.ff.charge.assign(n, 0.0);
  f.ff.mass.assign(n, 12.0);
  f.ff.excluded.assign(n, {});
  for (const auto& b : bonds) {
    f.ff.bonds.push_back({b[0], b[1], k, 1.0});
    f.ff.excluded[b[0]].push_back(b[1]), f.ff.excluded[b[1]].push_back(b[0]);
  }
  return f;
}

// The frequency of each topology over many moves at frozen positions (key: the sorted bond list).
std::map<std::vector<std::array<uint32_t, 2>>, double> sample(ExchangeEngine& e, System& s, double T, int moves) {
  std::map<std::vector<std::array<uint32_t, 2>>, double> h;
  const auto x = coords(s);
  for (int m = 0; m < moves; ++m) {
    e.run(x, s.cell, T, nullptr);
    std::vector<std::array<uint32_t, 2>> key;
    for (const auto& b : s.bonds) key.push_back({std::min(b.i, b.j), std::max(b.i, b.j)});
    std::sort(key.begin(), key.end());
    h[key] += 1.0 / moves;
  }
  return h;
}

double bond_energy(const System& s, const std::vector<std::array<uint32_t, 2>>& bonds, double k) {
  double u = 0;
  for (const auto& b : bonds) {
    const double r = norm(s.atoms[b[1]].pos - s.atoms[b[0]].pos);
    u += k * (r - 1) * (r - 1);
  }
  return u;
}

System two_butanes() {
  MolGraph g = parse_smiles("CCCC");
  add_hydrogens(g);
  const System a = molecule_system(g, embed(g));
  System s = a;
  const uint32_t n = uint32_t(a.atoms.size());
  for (auto& at : s.atoms) at.mol = 1, at.pos = at.pos + Vec3{10, 10, 10};
  for (auto at : a.atoms) at.mol = 2, at.pos = at.pos + Vec3{13, 10.5, 10}, s.atoms.push_back(at);
  for (auto b : a.bonds) b.i += n, b.j += n, s.bonds.push_back(b);
  s.has_mol = true;
  s.cell = box(30, 30, 30);
  return s;
}

System kg_grid(int side, int beads, double spacing = 1.15) {
  System s;
  s.cell = box(beads * 0.97 + 0.2, side * spacing, side * spacing);
  for (int cy = 0; cy < side; ++cy)
    for (int cz = 0; cz < side; ++cz) {
      const uint32_t first = uint32_t(s.atoms.size());
      for (int k = 0; k < beads; ++k) {
        Atom a;
        a.mol = 1 + cy * side + cz;
        a.type = 1;
        a.pos = {0.1 + 0.97 * k, (cy + 0.5) * spacing, (cz + 0.5) * spacing};
        s.atoms.push_back(a);
        if (k) s.bonds.push_back({first + uint32_t(k) - 1, first + uint32_t(k), 1});
      }
    }
  s.has_mol = true;
  return s;
}

}  // namespace

// ΔU from the atoms within three bonds of the change equals the whole system's, with pairwise (DSF) and with Ewald (PME)
// electrostatics: the k-space sums and every pair outside the cluster are the same before and after.
TEST(ExchangeMoves, LocalEnergyEqualsTheWholeSystems) {
  for (auto el : {EnergyOptions::Electrostatics::DSF, EnergyOptions::Electrostatics::PME}) {
    System s = two_butanes();
    const uint32_t n = uint32_t(s.atoms.size() / 2);
    ForceField ff = assign_gaff(s);
    EnergyOptions eo;
    eo.electrostatics = el;
    ExchangeOptions xo;
    xo.group_a = {0, 1, n, n + 1};
    ExchangeEngine local(s, ff, xo, eo);
    const auto x = coords(s);
    auto whole = [&] { return whole_energy(ff, s, eo); };
    const double u0 = whole();
    // swap the two butanes' methyls, and break one C–C bond alone
    const double d1 = local.delta({{0, 1}, {n, n + 1}}, {{1, n}, {0, n + 1}}, x, s.cell, nullptr);
    const double d2 = local.delta({{1, 2}}, {}, x, s.cell, nullptr);
    EXPECT_NEAR(whole(), u0, 1e-9 * std::fabs(u0));   // delta leaves nothing changed
    ExchangeOptions full = xo;
    full.local = false;
    ExchangeEngine global(s, ff, full, eo);
    EXPECT_NEAR(d1, global.delta({{0, 1}, {n, n + 1}}, {{1, n}, {0, n + 1}}, x, s.cell, whole), 1e-7 * (1 + std::fabs(d1)));
    EXPECT_NEAR(d2, global.delta({{1, 2}}, {}, x, s.cell, whole), 1e-7 * (1 + std::fabs(d2)));
    EXPECT_GT(std::fabs(d1), 1.0);   // the methyls are far from their new partners: a real energy change
  }
  // bead–spring chains with bending (1-3 and 1-4 pairs in full)
  System s = kg_grid(3, 8);
  KgOptions ko;
  ko.k_theta = 1.5;
  ForceField ff = kremer_grest_forcefield(s, ko);
  EnergyOptions eo;
  ExchangeOptions xo;
  for (uint32_t i = 0; i < s.atoms.size(); ++i) xo.group_a.push_back(i);
  ExchangeEngine local(s, ff, xo, eo);
  ExchangeOptions full = xo;
  full.local = false;
  ExchangeEngine global(s, ff, full, eo);
  const auto x = coords(s);
  auto whole = [&] { return whole_energy(ff, s, eo); };
  const double d = local.delta({{3, 4}, {11, 12}}, {{3, 11}, {4, 12}}, x, s.cell, nullptr);
  EXPECT_NEAR(d, global.delta({{3, 4}, {11, 12}}, {{3, 11}, {4, 12}}, x, s.cell, whole), 1e-8 * (1 + std::fabs(d)));
}

// Swaps among four atoms at fixed places visit the three pairings with Boltzmann weights; a barrier slows the moves but
// leaves the weights alone.
TEST(ExchangeMoves, SwapSamplesBoltzmannWeights) {
  const std::vector<Vec3> pos{{0, 0, 0}, {1.0, 0, 0}, {1.1, 0.9, 0}, {0.1, 1.05, 0}};
  const double k = 5, T = 300;
  for (double barrier : {0.0, 0.3}) {
    Frozen f = frozen(pos, {{0, 1}, {2, 3}}, k);
    ExchangeOptions xo;
    xo.group_a = {0, 1, 2, 3};
    xo.capture = 1.6;
    xo.barrier = barrier;
    xo.seed = 3;
    ExchangeEngine e(f.s, f.ff, xo, EnergyOptions{});
    const auto h = sample(e, f.s, T, 60000);
    const std::vector<std::vector<std::array<uint32_t, 2>>> states{{{0, 1}, {2, 3}}, {{0, 3}, {1, 2}}, {{0, 2}, {1, 3}}};
    double z = 0;
    for (const auto& st : states) z += std::exp(-bond_energy(f.s, st, k) / (kB * T));
    for (const auto& st : states) {
      const double p = std::exp(-bond_energy(f.s, st, k) / (kB * T)) / z;
      const auto it = h.find(st);
      EXPECT_NEAR(it == h.end() ? 0.0 : it->second, p, 0.015) << "barrier " << barrier;
    }
    if (barrier > 0) EXPECT_LT(e.stats().accepted, e.stats().attempts);
  }
}

// Two pivots share three residues: with one valence per residue the pivots never hold the same one, with two they may;
// either way every arrangement comes with its Boltzmann weight.
TEST(ExchangeMoves, ShiftRespectsValencesAndSamplesBoltzmannWeights) {
  // pivots 0 and 4, residues 1, 2, 3 at different distances from each
  const std::vector<Vec3> pos{{0, 0, 0}, {1.0, 0, 0}, {0, 1.1, 0}, {0.3, 0.3, 1.0}, {1.0, 1.0, 0.3}};
  const double k = 5, T = 300;
  for (int valence : {1, 2}) {
    Frozen f = frozen(pos, {{0, 1}, {2, 4}}, k);
    ExchangeOptions xo;
    xo.kind = ExchangeKind::Shift;
    xo.group_a = {0, 4};
    xo.group_b = {1, 2, 3};
    xo.valence_b = valence;
    xo.capture = 1.6;
    xo.seed = 5;
    ExchangeEngine e(f.s, f.ff, xo, EnergyOptions{});
    const auto h = sample(e, f.s, T, 60000);
    double z = 0;
    std::map<std::vector<std::array<uint32_t, 2>>, double> w;
    for (uint32_t rp : {1u, 2u, 3u})
      for (uint32_t rq : {1u, 2u, 3u}) {
        if (rp == rq && valence == 1) continue;
        std::vector<std::array<uint32_t, 2>> st{{0, rp}, {std::min(rq, 4u), std::max(rq, 4u)}};
        std::sort(st.begin(), st.end());
        w[st] = std::exp(-bond_energy(f.s, st, k) / (kB * T));
        z += w[st];
      }
    for (const auto& [st, wt] : h) EXPECT_TRUE(w.count(st)) << "a state the valences forbid was visited";
    for (const auto& [st, wt] : w) {
      const auto it = h.find(st);
      EXPECT_NEAR(it == h.end() ? 0.0 : it->second, wt / z, 0.015) << "valence " << valence;
    }
  }
}

// Binding and breaking among two A and two B sites at fixed places with no bond energy: the populations of 0, 1 and 2
// bonds follow Z = 1 + 4K + 2K² (Wilcoxson et al. 2026, their four-reaction check).
TEST(ExchangeMoves, BindSamplesThePartitionFunction) {
  const std::vector<Vec3> pos{{0, 0, 0}, {1.0, 0, 0}, {0, 1.0, 0}, {1.0, 1.0, 0}};
  const double T = 300, K = 1.5;
  Frozen f = frozen(pos, {}, 0);
  ExchangeOptions xo;
  xo.kind = ExchangeKind::Bind;
  xo.group_a = {0, 3};
  xo.group_b = {1, 2};
  xo.capture = 1.6;
  xo.binding = kB * T * std::log(K);
  xo.seed = 7;
  ExchangeEngine e(f.s, f.ff, xo, EnergyOptions{});
  e.editor().define_bond(BondTerm{0, 1, 0.0, 1.0});
  const auto h = sample(e, f.s, T, 200000);
  double p[3] = {0, 0, 0};
  for (const auto& [st, wt] : h) p[st.size()] += wt;
  const double z = 1 + 4 * K + 2 * K * K;
  EXPECT_NEAR(p[0], 1 / z, 0.01);
  EXPECT_NEAR(p[1], 4 * K / z, 0.01);
  EXPECT_NEAR(p[2], 2 * K * K / z, 0.01);
  for (const auto& [st, wt] : h)
    for (const auto& b : st) EXPECT_TRUE((b[0] == 0 || b[0] == 3) != (b[1] == 0 || b[1] == 3)) << "only A–B bonds";
}

// Exchange during CAPS's own MD: swaps happen on a melt, the bond count stays and no bead holds more than two bonds.
TEST(ExchangeMoves, SwapsDuringDynamics) {
  System s = kg_grid(4, 10);
  KgOptions ko;
  DynamicsOptions d;
  d.field = std::make_shared<ForceField>(kremer_grest_forcefield(s, ko));
  d.temperature = 503.2;
  d.dt = 0.25;
  d.steps = 1500;
  d.tau_t = 25;
  run_dynamics(s, d);   // melted first
  d.steps = 1500;
  d.topology_every = 10;
  ExchangeOptions xo;
  for (uint32_t i = 0; i < s.atoms.size(); ++i) xo.group_a.push_back(i);
  xo.capture = 1.3;
  xo.attempts = 20;
  auto st = std::make_shared<ExchangeStats>();
  d.topology = exchange_callback(xo, d.energy, st);
  DynamicsReport r;
  run_dynamics(s, d, &r);
  EXPECT_GT(st->accepted, 0);
  EXPECT_EQ(r.topology_changes > 0, st->accepted > 0);
  EXPECT_TRUE(st->local);
  EXPECT_EQ(s.bonds.size(), size_t(16 * 9));
  std::vector<int> deg(s.atoms.size(), 0);
  for (const auto& b : s.bonds) ++deg[b.i], ++deg[b.j];
  for (int k : deg) EXPECT_LE(k, 2);
}

// Labels as LAMMPS molecule IDs: with each bead's place along its chain, a1 and b1 share theirs and so do a2 and b2, so
// every swap keeps chain lengths; without labels the same moves make chains of other lengths.
TEST(ExchangeMoves, ChainLabelsKeepChainLengths) {
  // two chains of four beads, side by side 0.6 apart, harmonic bonds and no pair terms
  std::vector<Vec3> pos;
  for (double y : {0.0, 0.6})
    for (int k = 0; k < 4; ++k) pos.push_back({double(k), y, 0});
  for (bool labelled : {true, false}) {
    Frozen f = frozen(pos, {{0, 1}, {1, 2}, {2, 3}, {4, 5}, {5, 6}, {6, 7}}, 5);
    ExchangeOptions xo;
    xo.group_a = {0, 1, 2, 3, 4, 5, 6, 7};
    if (labelled) xo.labels = chain_position_labels(f.s);
    xo.capture = 1.5;
    xo.seed = 11;
    ExchangeEngine e(f.s, f.ff, xo, EnergyOptions{});
    const auto h = sample(e, f.s, 300, 4000);
    bool other_lengths = false;
    for (const auto& [st, wt] : h) {
      std::vector<std::vector<uint32_t>> nb(8);
      for (const auto& b : st) nb[b[0]].push_back(b[1]), nb[b[1]].push_back(b[0]);
      std::vector<char> seen(8, 0);
      for (uint32_t a = 0; a < 8; ++a) {
        if (seen[a]) continue;
        std::vector<uint32_t> stack{a};
        seen[a] = 1;
        int size = 0;
        while (!stack.empty()) {
          const uint32_t c = stack.back();
          stack.pop_back();
          ++size;
          for (uint32_t q : nb[c])
            if (!seen[q]) seen[q] = 1, stack.push_back(q);
        }
        if (size != 4) other_lengths = true;
      }
    }
    EXPECT_GT(e.stats().accepted, 0) << (labelled ? "labelled" : "unlabelled");
    EXPECT_EQ(other_lengths, !labelled);
  }
}

TEST(ExchangeMoves, ChainPositionLabels) {
  const System s = kg_grid(1, 6);
  EXPECT_EQ(chain_position_labels(s), (std::vector<int64_t>{1, 2, 3, 4, 5, 6}));
  EXPECT_EQ(chain_position_labels(s, true), (std::vector<int64_t>{1, 2, 3, 3, 2, 1}));
}

// The LAMMPS deck: the fix with its arguments, the 0 1 1 special bonds, molecule IDs for either use, newton off for bonds
// with stickers.
TEST(ExchangeMoves, BondSwapDeckForLammps) {
  System s = kg_grid(2, 6);
  KgOptions ko;
  ko.k_theta = 1.5;
  const auto dir = std::filesystem::temp_directory_path();
  auto read = [](const std::string& p) {
    std::ifstream in(p);
    return std::string((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
  };
  BondSwapExport e;
  const std::string chain = (dir / "caps_bondswap_chain").string();
  write_kg_bond_swap(s, ko, e, chain);
  const std::string in = read(chain + ".in"), data = read(chain + ".data");
  for (const char* key : {"fix swap exch bond/swap 50 0.5 1.3", "special_bonds fene", "group exch type 2", "f_swap[1] f_swap[2]", "atom_style angle"})
    EXPECT_NE(in.find(key), std::string::npos) << key;
  EXPECT_EQ(in.find("newton on off"), std::string::npos);
  EXPECT_NE(data.find("\n1 1 2 "), std::string::npos);   // atom 1: molecule ID 1 (first on its chain), type 2
  EXPECT_NE(data.find("\n2 2 2 "), std::string::npos);   // atom 2: molecule ID 2
  EXPECT_NE(data.find("16 angles"), std::string::npos);   // four chains of six beads: four angles each
  // stickers: the chain ends, two of them bonded across the chains (atoms 6 and 7)
  s.bonds.push_back({5, 6, 1});
  e.stickers = true;
  e.group = {0, 5, 6, 11};
  const std::string st = (dir / "caps_bondswap_sticker").string();
  write_kg_bond_swap(s, ko, e, st);
  const std::string in2 = read(st + ".in"), data2 = read(st + ".data");
  EXPECT_NE(in2.find("newton on off"), std::string::npos);
  EXPECT_NE(data2.find("\n6 1 2 "), std::string::npos);   // sticker: molecule ID 1, type 2
  EXPECT_NE(data2.find("\n2 2 1 "), std::string::npos);   // not a sticker: molecule ID 2, type 1
  EXPECT_NE(data2.find("16 angles"), std::string::npos);   // none across the sticker–sticker bond
}
