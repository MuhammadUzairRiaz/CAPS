#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <memory>
#include <random>
#include <vector>

#include "caps/dynamics.hpp"
#include "caps/field.hpp"
#include "caps/kremer_grest.hpp"
#include "caps/molecule.hpp"
#include "caps/topology_edit.hpp"

using namespace caps;

namespace {

Cell box(double x, double y, double z) {
  Cell c;
  c.a = {x, 0, 0};
  c.b = {0, y, 0};
  c.c = {0, 0, z};
  return c;
}

// Two copies of a molecule, the second shifted along x, in a large cell.
System two_copies(const std::string& smiles, double shift) {
  MolGraph g = parse_smiles(smiles);
  add_hydrogens(g);
  const System a = molecule_system(g, embed(g));
  System s = a;
  const uint32_t n = uint32_t(a.atoms.size());
  for (auto& at : s.atoms) at.mol = 1, at.pos = at.pos + Vec3{10, 10, 10};
  for (auto at : a.atoms) {
    at.mol = 2, at.pos = at.pos + Vec3{10 + shift, 10, 10};
    s.atoms.push_back(at);
  }
  for (auto b : a.bonds) b.i += n, b.j += n, s.bonds.push_back(b);
  s.has_mol = true;
  s.cell = box(40, 40, 40);
  return s;
}

EnergyTerms energy(const ForceField& ff, const System& s, Evaluator* ev = nullptr) {
  std::vector<double> x, f;
  for (const auto& a : s.atoms) x.insert(x.end(), {a.pos[0], a.pos[1], a.pos[2]});
  if (ev) return ev->compute(x, s.cell, f);
  Evaluator e(ff, EnergyOptions{});
  return e.compute(x, s.cell, f);
}

void expect_same(const EnergyTerms& a, const EnergyTerms& b) {
  auto near = [](double p, double q) { return std::fabs(p - q) <= 1e-8 * (1 + std::fabs(q)); };
  EXPECT_TRUE(near(a.bond, b.bond)) << a.bond << " vs " << b.bond;
  EXPECT_TRUE(near(a.angle, b.angle)) << a.angle << " vs " << b.angle;
  EXPECT_TRUE(near(a.dihedral, b.dihedral)) << a.dihedral << " vs " << b.dihedral;
  EXPECT_TRUE(near(a.improper, b.improper)) << a.improper << " vs " << b.improper;
  EXPECT_TRUE(near(a.vdw, b.vdw)) << a.vdw << " vs " << b.vdw;
  EXPECT_TRUE(near(a.coulomb, b.coulomb)) << a.coulomb << " vs " << b.coulomb;
}

// Straight Kremer–Grest chains on a square grid (bonds 0.97 σ along x, chains `spacing` apart): no overlaps to push off.
System kg_grid(int side, int beads, double spacing = 1.15) {
  System s;
  const double lx = beads * 0.97 + 0.2;
  s.cell = box(lx, side * spacing, side * spacing);
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

// A swap of groups between two copies of a molecule leaves two copies of it: every term the editor rebuilds (bonds,
// angles, torsions, impropers, 1-4 pairs, exclusions) must be what the force field gives the new structure from scratch.
TEST(Topology, SwapMatchesAFreshAssignmentGaff) {
  for (const char* smiles : {"CCCC", "Cc1ccccc1"}) {
    System s = two_copies(smiles, 8.0);
    const uint32_t n = uint32_t(s.atoms.size() / 2);
    const ForceField ff0 = assign_gaff(s);
    for (size_t i = 0; i < s.atoms.size(); ++i) s.atoms[i].charge = ff0.charge[i];
    s.has_charges = true;   // the same charges for every assignment
    ForceField ff = assign_gaff(s);
    const EnergyTerms e0 = energy(ff, s);
    Evaluator ev(ff, EnergyOptions{});
    energy(ff, s, &ev);
    TopologyEditor ed(s, ff);
    EXPECT_TRUE(ed.excludes_13());
    EXPECT_TRUE(ed.excludes_14());
    // atom 1 (CH2, or the ring carbon) takes the other copy's atom 0 (the methyl), and the other way round
    ed.swap(1, 0, n + 1, n);
    EXPECT_TRUE(ed.bonded(1, n));
    EXPECT_TRUE(ed.bonded(n + 1, 0));
    EXPECT_FALSE(ed.bonded(1, 0));
    const ForceField fresh = assign_gaff(s);
    EXPECT_EQ(ff.bonds.size(), fresh.bonds.size()) << smiles;
    EXPECT_EQ(ff.angles.size(), fresh.angles.size()) << smiles;
    EXPECT_EQ(ff.dihedrals.size(), fresh.dihedrals.size()) << smiles;
    EXPECT_EQ(ff.impropers.size(), fresh.impropers.size()) << smiles;
    EXPECT_EQ(ff.pairs14.size(), fresh.pairs14.size()) << smiles;
    EXPECT_EQ(ff.excluded, fresh.excluded) << smiles;
    EXPECT_EQ(ed.missing_impropers(), 0u);
    const EnergyTerms e1 = energy(fresh, s);
    expect_same(energy(ff, s), e1);
    // the evaluator built before the change gives the same after a reload
    ev.reload();
    expect_same(energy(ff, s, &ev), e1);
    // and back: the original energy exactly
    ed.swap(1, n, n + 1, 0);
    expect_same(energy(ff, s), e0);
    EXPECT_EQ(ed.changes(), 2u);
  }
}

// Bead–spring chains with bending: FENE bonds and cosine angles follow a swap; 1-3 and 1-4 pairs stay in (special_bonds fene).
TEST(Topology, SwapMatchesAFreshAssignmentKremerGrest) {
  System s = kg_grid(2, 8, 1.0);
  KgOptions o;
  o.k_theta = 1.5;
  ForceField ff = kremer_grest_forcefield(s, o);
  TopologyEditor ed(s, ff);
  EXPECT_FALSE(ed.excludes_13());
  EXPECT_FALSE(ed.excludes_14());
  // chain 1 beads 3–4 and chain 2 beads 3–4 (side by side): 3 joins the other chain's 4
  ed.swap(3, 4, 8 + 3, 8 + 4);
  const ForceField fresh = kremer_grest_forcefield(s, o);
  EXPECT_EQ(ff.bonds_x.size(), fresh.bonds_x.size());
  EXPECT_EQ(ff.angles_x.size(), fresh.angles_x.size());
  EXPECT_EQ(ff.excluded, fresh.excluded);
  expect_same(energy(ff, s), energy(fresh, s));
}

// A change the force field cannot describe leaves everything as it was.
TEST(Topology, ImpossibleChangeLeavesStructureAndFieldAlone) {
  System s = two_copies("CCCC", 8.0);
  ForceField ff = assign_gaff(s);
  const auto bonds = s.bonds;
  const auto terms = ff.angles.size();
  TopologyEditor ed(s, ff);
  EXPECT_THROW(ed.change({}, {{0, 1}}), TopologyError);          // already bonded
  EXPECT_THROW(ed.change({{0, 2}}, {}), TopologyError);          // no bond to break
  // a hydrogen bonded to a second carbon: H-C bonds exist, but an H with two bonds makes angles GAFF never had (C-H-C)
  EXPECT_THROW(ed.change({}, {{4, 15}}), TopologyError);
  EXPECT_EQ(s.bonds.size(), bonds.size());
  EXPECT_EQ(ff.angles.size(), terms);
  EXPECT_EQ(ed.changes(), 0u);
  // breaking a bond and making it again restores it, its own parameters included
  const EnergyTerms e0 = energy(ff, s);
  ed.change({{0, 1}}, {});
  EXPECT_FALSE(ed.bonded(0, 1));
  ed.change({}, {{0, 1}});
  expect_same(energy(ff, s), e0);
}

// NVE dynamics of a melted bead–spring system with Metropolis bond swaps every 10 steps (new bonds shorter than 1.25 σ):
// the energy the swaps add is kept out of the conserved quantity, so it moves from row to row no more than without
// swaps, though the swaps moved far more energy than that; a rejected swap that is undone leaves the energy exactly as it
// was; the bonds stay two per bead at most.
TEST(Topology, SwapsDuringNveKeepTheConservedQuantity) {
  // the grid melted first (Bussi thermostat), so the runs below start from the same disordered state
  System melt = kg_grid(6, 10);
  KgOptions o;
  auto ff = std::make_shared<ForceField>(kremer_grest_forcefield(melt, o));
  DynamicsOptions d;
  d.field = ff;
  d.temperature = 503.2;   // T* = 1 (σ = 1 Å, ε = 1 kcal/mol, m = 1 g/mol)
  d.dt = 0.25;             // about 0.005 τ
  d.steps = 2000;
  d.tau_t = 25;
  d.seed = 5;
  run_dynamics(melt, d);
  auto run = [&](bool swaps, DynamicsReport& rep, size_t& accepted, double& moved) {
    System s = melt;
    DynamicsOptions n = d;
    n.thermostat = Thermostat::None;
    n.steps = 2000;
    n.thermo_every = 10;
    std::unique_ptr<TopologyEditor> ed;
    std::mt19937_64 rng(9);
    std::uniform_real_distribution<double> u01(0, 1);
    accepted = 0, moved = 0;
    if (swaps) {
      n.topology_every = 10;
      n.topology = [&](TopologyContext& c) {
        if (!ed) ed = std::make_unique<TopologyEditor>(c.system, c.field);
        auto dist = [&](uint32_t i, uint32_t j) {
          Vec3 r{c.x[3 * j] - c.x[3 * i], c.x[3 * j + 1] - c.x[3 * i + 1], c.x[3 * j + 2] - c.x[3 * i + 2]};
          return norm(c.cell.minimum_image(r));
        };
        bool any = false;
        double u = c.potential();
        for (int t = 0; t < 40; ++t) {
          const auto bs = c.system.bonds;
          const auto& a = bs[size_t(u01(rng) * bs.size()) % bs.size()];
          const uint32_t a1 = u01(rng) < 0.5 ? a.i : a.j, a2 = a1 == a.i ? a.j : a.i;
          int64_t pick = -1;
          const size_t start = size_t(u01(rng) * bs.size());
          for (size_t k = 0; k < bs.size() && pick < 0; ++k) {
            const auto& b = bs[(k + start) % bs.size()];
            for (int o2 = 0; o2 < 2 && pick < 0; ++o2) {
              const uint32_t b1 = o2 ? b.j : b.i, b2 = o2 ? b.i : b.j;
              if (b1 == a1 || b1 == a2 || b2 == a1 || b2 == a2) continue;
              if (ed->bonded(a1, b2) || ed->bonded(b1, a2)) continue;
              if (dist(a1, b2) < 1.25 && dist(b1, a2) < 1.25) pick = int64_t(b1) << 32 | b2;
            }
          }
          if (pick < 0) continue;
          const uint32_t b1 = uint32_t(pick >> 32), b2 = uint32_t(pick & 0xffffffff);
          // Metropolis; a rejected trial is undone and must leave the energy exactly as it was
          ed->swap(a1, a2, b1, b2);
          const double u1 = c.potential();
          if (u1 <= u || u01(rng) < std::exp(-(u1 - u) / (0.0019872067 * c.temperature))) {
            moved += std::fabs(u1 - u);
            u = u1, any = true, ++accepted;
          } else {
            ed->swap(a1, b2, b1, a2);
            EXPECT_NEAR(c.potential(), u, 1e-9 * (1 + std::fabs(u)));
          }
        }
        return any;
      };
    }
    run_dynamics(s, n, &rep);
    for (size_t i = 0; i < s.atoms.size(); ++i) {
      size_t deg = 0;
      for (const auto& b : s.bonds) deg += b.i == i || b.j == i;
      EXPECT_LE(deg, 2u);
    }
    EXPECT_EQ(s.bonds.size(), size_t(36 * 9));
  };
  // the largest change of the conserved quantity between consecutive rows (10 steps, one swap call apart)
  auto jitter = [](const DynamicsReport& r) {
    double m = 0;
    for (size_t k = 1; k < r.thermo.size(); ++k) m = std::max(m, std::fabs(r.thermo[k].conserved - r.thermo[k - 1].conserved));
    return m;
  };
  DynamicsReport ref, ex;
  size_t acc0 = 0, acc = 0;
  double moved0 = 0, moved = 0;
  run(false, ref, acc0, moved0);
  run(true, ex, acc, moved);
  EXPECT_GT(acc, 10u);
  EXPECT_GT(ex.topology_changes, 5);
  const double j0 = jitter(ref), j1 = jitter(ex);
  EXPECT_LT(j1, 3 * j0 + 1e-6) << "row-to-row change of the conserved quantity without swaps " << j0 << ", with " << j1;
  EXPECT_GT(moved, 10 * j0) << "the swaps must move more energy than the integrator's own jitter for the check to mean anything";
}

// Bond constraints and bond exchange are alternatives.
TEST(Topology, NoExchangeWithConstraints) {
  System s = two_copies("CCCC", 8.0);
  DynamicsOptions d;
  d.field = std::make_shared<ForceField>(assign_gaff(s));
  d.constraints = ConstraintMode::HBonds;
  d.steps = 1;
  d.topology_every = 1;
  d.topology = [](TopologyContext&) { return false; };
  EXPECT_THROW(run_dynamics(s, d), std::invalid_argument);
}

// The general route: the callback changes the structure's bonds and puts in the force field assigned again to it (any
// chemistry, new types included); the run goes on with it.
TEST(Topology, CallbackCanReassignTheForceField) {
  System s = two_copies("CCCC", 6.0);
  const ForceField ff0 = assign_gaff(s);
  for (size_t i = 0; i < s.atoms.size(); ++i) s.atoms[i].charge = ff0.charge[i];
  s.has_charges = true;
  const uint32_t n = uint32_t(s.atoms.size() / 2);
  DynamicsOptions d;
  d.field = std::make_shared<ForceField>(assign_gaff(s));
  d.thermostat = Thermostat::None;
  d.steps = 20;
  d.dt = 0.5;
  d.thermo_every = 5;
  d.topology_every = 10;
  int calls = 0;
  d.topology = [&](TopologyContext& c) {
    if (calls++) return false;
    auto& b = c.system.bonds;
    auto drop = [&](uint32_t i, uint32_t j) {
      b.erase(std::remove_if(b.begin(), b.end(), [&](const Bond& x) { return (x.i == i && x.j == j) || (x.i == j && x.j == i); }), b.end());
    };
    drop(1, 0), drop(n + 1, n);
    b.push_back({1, n, 1}), b.push_back({n + 1, 0, 1});
    c.field = assign_gaff(c.system);
    return true;
  };
  DynamicsReport r;
  run_dynamics(s, d, &r);
  EXPECT_EQ(r.topology_changes, 1);
  EXPECT_NE(r.topology_energy, 0.0);
  bool swapped = false;
  for (const auto& b : s.bonds) swapped |= (b.i == 1 && b.j == n) || (b.i == n && b.j == 1);
  EXPECT_TRUE(swapped);
  for (const auto& row : r.thermo) EXPECT_TRUE(std::isfinite(row.conserved));
}
