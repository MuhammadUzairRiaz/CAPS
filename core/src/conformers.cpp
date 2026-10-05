#include "caps/conformers.hpp"

#include <algorithm>
#include <cmath>
#include <random>
#include <stdexcept>

#include "caps/dynamics.hpp"
#include "caps/relax.hpp"
#include "caps/superpose.hpp"
#include "caps/torsion.hpp"
#include "caps/uff.hpp"

namespace caps {

std::vector<std::array<uint32_t, 4>> rotatable_bonds(const System& s) {
  const auto nb = s.neighbours();
  std::vector<std::array<uint32_t, 4>> out;
  auto heavy_other = [&](uint32_t at, uint32_t not_this) -> int64_t {
    for (uint32_t k : nb[at])
      if (k != not_this && s.atoms[k].element > 1) return k;
    return -1;
  };
  for (const auto& b : s.bonds) {
    if (b.order > 1) continue;   // double, triple, aromatic, amide, dative
    const int64_t a = heavy_other(b.i, b.j), d = heavy_other(b.j, b.i);
    if (a < 0 || d < 0) continue;
    try {
      moving_side(s, int(b.i), int(b.j));   // throws for a bond in a ring
    } catch (const std::exception&) {
      continue;
    }
    out.push_back({uint32_t(a), b.i, b.j, uint32_t(d)});
  }
  return out;
}

ConformerSearchResult conformer_search(const System& s0, const ConformerSearchOptions& o) {
  ConformerSearchResult r;
  if (s0.atoms.size() < 4) throw std::invalid_argument("a conformer search needs a molecule of at least four atoms");
  System s = s0;
  s.cell = Cell{};   // molecules in vacuum
  s.velocities.clear();
  const ForceField ff = o.field ? *o.field : default_forcefield(s);
  auto ffp = std::make_shared<const ForceField>(ff);
  r.field = ff.name;
  r.method = o.method;
  r.rotor_atoms = rotatable_bonds(s);
  r.rotors = int(r.rotor_atoms.size());
  std::vector<std::vector<uint32_t>> moving;
  for (const auto& d : r.rotor_atoms) moving.push_back(moving_side(s, int(d[1]), int(d[2])));
  RelaxOptions ro;
  ro.field = ffp;
  ro.energy = o.energy;
  ro.ftol = o.ftol;
  ro.max_iterations = 5000;
  std::mt19937_64 rng(o.seed);
  std::normal_distribution<double> jitter(0.0, 15.0);
  std::uniform_int_distribution<int> pick(0, 2);
  struct Min { std::vector<Vec3> pos; double e; };
  std::vector<Min> mins;
  auto minimise = [&](System& t) {
    RelaxReport rep;
    relax(t, ro, &rep);
    std::vector<Vec3> p;
    for (const auto& a : t.atoms) p.push_back(a.pos);
    mins.push_back({std::move(p), rep.final.total()});
    r.minima.push_back(mins.back().e);
  };
  if (o.method == "anneal") {
    // a hot NVT run, snapshots every anneal_ps quenched
    System hot = s;
    {
      System first = s;
      minimise(first);
    }
    DynamicsOptions d;
    d.field = ffp;
    d.energy = o.energy;
    d.temperature = o.anneal_temperature;
    d.thermostat = Thermostat::Bussi;
    d.tau_t = 100;
    d.dt = 0.5;
    d.seed = o.seed;
    d.new_velocities = true;
    d.frame_every = std::max(1, int(std::lround(o.anneal_ps * 1000 / d.dt)));
    d.steps = int64_t(d.frame_every) * std::max(1, o.trials - 1);
    std::vector<std::vector<double>> snaps;
    d.frame = [&](const std::vector<double>& x, const Cell&, int64_t step) { if (step > 0) snaps.push_back(x); };
    d.progress = [&](const ThermoRow& row) { return !o.progress || o.progress(0.5 * double(row.step) / double(d.steps)); };
    run_dynamics(hot, d);
    for (size_t k = 0; k < snaps.size(); ++k) {
      if (o.progress && !o.progress(0.5 + 0.5 * double(k) / snaps.size())) throw std::runtime_error("conformer search cancelled");
      System t = s;
      for (size_t i = 0; i < t.atoms.size(); ++i) t.atoms[i].pos = {snaps[k][3 * i], snaps[k][3 * i + 1], snaps[k][3 * i + 2]};
      minimise(t);
    }
    r.notes.push_back(std::to_string(snaps.size()) + " snapshots of NVT at " + std::to_string(int(o.anneal_temperature)) + " K, " +
                      std::to_string(o.anneal_ps).substr(0, 4) + " ps apart, each minimised");
  } else {
    for (int t = 0; t < std::max(1, o.trials); ++t) {
      if (o.progress && !o.progress(double(t) / std::max(1, o.trials))) throw std::runtime_error("conformer search cancelled");
      System trial = s;
      if (t > 0)
        for (size_t k = 0; k < r.rotor_atoms.size(); ++k) {
          const double phi = 60.0 + 120.0 * pick(rng) + jitter(rng);
          set_dihedral(trial, {int(r.rotor_atoms[k][0]), int(r.rotor_atoms[k][1]), int(r.rotor_atoms[k][2]), int(r.rotor_atoms[k][3])},
                       phi > 180 ? phi - 360 : phi, moving[k]);
        }
      minimise(trial);
      if (r.rotors == 0) break;   // rigid: one minimum
    }
    if (r.rotors == 0) r.notes.push_back("no rotatable bonds: one conformer (rings are not sampled)");
  }
  // clusters, lowest first
  std::vector<size_t> ord(mins.size());
  for (size_t k = 0; k < ord.size(); ++k) ord[k] = k;
  std::sort(ord.begin(), ord.end(), [&](size_t a, size_t b) { return mins[a].e < mins[b].e; });
  std::vector<double> w(s.atoms.size());
  for (size_t i = 0; i < w.size(); ++i) w[i] = s.atoms[i].element > 1 ? 1.0 : 0.0;
  size_t heavy = 0;
  for (double x : w) heavy += x > 0;
  if (heavy < 3) std::fill(w.begin(), w.end(), 1.0);
  const double emin = mins.empty() ? 0 : mins[ord[0]].e;
  for (size_t k : ord) {
    if (mins[k].e - emin > o.window) break;
    bool joined = false;
    for (auto& c : r.conformers) {
      const auto fit = superpose(c.pos, mins[k].pos, w);
      if (fit.rmsd < o.rmsd) { ++c.found; joined = true; break; }
    }
    if (!joined) {
      ConformerHit h;
      h.pos = mins[k].pos;
      h.energy = mins[k].e;
      h.relative = mins[k].e - emin;
      h.found = 1;
      r.conformers.push_back(std::move(h));
    }
  }
  constexpr double kR = 0.0019872043;
  double z = 0;
  for (const auto& c : r.conformers) z += std::exp(-c.relative / (kR * o.temperature));
  for (auto& c : r.conformers) c.population = std::exp(-c.relative / (kR * o.temperature)) / z;
  r.notes.push_back("populations from the minimised energies alone (Boltzmann at " + std::to_string(int(std::lround(o.temperature))) +
                    " K): no vibrational entropy, no solvent");
  if (!r.conformers.empty() && r.conformers[0].found == 1 && mins.size() > 3)
    r.notes.push_back("the lowest conformer was found once: more trials may find lower ones");
  return r;
}

}  // namespace caps
