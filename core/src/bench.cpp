// CAPS Bench: the built-in validation suite (see bench.hpp).
#include "caps/bench.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <functional>
#include <map>
#include <memory>
#include <mutex>
#include <numeric>
#include <random>
#include <sstream>
#include <thread>

#include "caps/cg_bonded.hpp"
#include "caps/cg_nonbonded.hpp"
#include "caps/cg_rules.hpp"
#include "caps/config.hpp"
#include "caps/dynamics.hpp"
#include "caps/field.hpp"
#include "caps/io.hpp"
#include "caps/properties.hpp"
#include "caps/json.hpp"
#include "caps/kspace.hpp"
#include "caps/molecule.hpp"
#include "caps/pack.hpp"
#include "caps/relax.hpp"
#include "caps/polymer.hpp"
#include "caps/recipe.hpp"
#include "caps/render.hpp"
#include "caps/typing.hpp"

namespace caps {
namespace {

using Clock = std::chrono::steady_clock;
double since(Clock::time_point t) { return std::chrono::duration<double>(Clock::now() - t).count(); }

std::string fmt(const char* f, double v) {
  char b[64];
  std::snprintf(b, sizeof b, f, v);
  return b;
}
std::string num(double v, int digits = 3) { return fmt(("%." + std::to_string(digits) + "g").c_str(), v); }
std::string thousands(long long v) {
  std::string s = std::to_string(v), out;
  for (size_t i = 0; i < s.size(); ++i) {
    if (i && (s.size() - i) % 3 == 0) out += ",";
    out += s[i];
  }
  return out;
}
std::string mean_sd(const std::vector<double>& v, const char* f) {
  if (v.empty()) return "—";
  const double m = std::accumulate(v.begin(), v.end(), 0.0) / double(v.size());
  if (v.size() < 2) return fmt(f, m);
  double s = 0;
  for (double x : v) s += (x - m) * (x - m);
  return fmt(f, m) + " ± " + fmt(f, std::sqrt(s / double(v.size() - 1)));
}

struct Meta {
  const char* id;
  const char* title;
  const char* scope;
  std::vector<std::string> columns;
  const char* not_run;   // non-null: the built-in suite does not run this table, and why
};

const std::vector<Meta>& metas() {
  static const std::vector<Meta> m = {
      {"T1", "Forces and virial", "analytic vs finite differences of the energy", {"Check", "Atoms", "Max |ΔF| (kcal/mol/Å)", "RMS F", "Status"}, nullptr},
      {"T2", "Packing", "CAPS Pack on water boxes", {"System", "Molecules", "Atoms", "Wall (s)", "d min (Å)", "Success"}, nullptr},
      {"T3", "Amorphous cell builders", "vs Materials Studio, Polymatic, PySIMM, RadonPy", {}, "needs the other builders installed; compare their cells with CAPS Grow outside the Studio"},
      {"T4", "MD throughput", "ns/day, NVT with Bussi", {"System", "Atoms", "Threads", "ns/day", "Status"}, nullptr},
      {"T5", "NVE energy conservation", "drift, kT/ns/atom", {"System", "dt (fs)", "Length (ps)", "Drift (kT/ns/atom)", "RMS fluct. (kcal/mol)", "Status"}, nullptr},
      {"T6", "Properties vs experiment", "your equilibrated cells against data/reference/polymers.json", {"Material", "Cell", "Frames", "Density (g/cm³)", "Experiment", "Status"}, nullptr},
      {"T7", "Chain statistics", "C∞ of the same cells vs literature", {"Material", "Cell", "Chains", "C∞", "Literature", "Status"}, nullptr},
      {"T8", "Parallel scaling", "force evaluation, threads", {"System", "Threads", "ms / evaluation", "Speed-up", "Efficiency"}, nullptr},
      {"T9", "Reproducibility", "the same run twice", {"Run", "Steps", "Max |Δx| (Å)", "Status"}, nullptr},
      {"T10", "Builder feature coverage", "against vendor documentation", {}, "a written comparison, not a measurement"},
      {"T11", "Rendering", "CPU renderer, 1920 × 1080", {"System", "Atoms", "Supersampling", "ms / frame"}, nullptr},
      {"T12", "Molecule builder", "stereochemistry, round trips, embedding", {"Check", "Cases", "Correct", "Status"}, nullptr},
      {"T13", "Coarse-grained mapping", "ester-cut beads of grown PBS, PBSA, PBAT", {"System", "Beads", "Beads / unit", "Composition", "Mass error (g/mol)", "Status"}, nullptr},
      {"T14", "Bonded Boltzmann inversion", "an ideal chain drawn from known potentials, inverted", {"Term", "Samples", "Largest |ΔU| (kT)", "Status"}, nullptr},
      {"T15", "IBI self-consistency", "a two-type Lennard-Jones fluid's g(r) as the target, from −kT ln g", {"Iteration", "RDF error (%)", "|ΔP| (atm)", "Status"}, nullptr},
      {"T16", "Entanglement length of a Kremer–Grest melt", "PPA, N = 500, ρ = 0.85, k_θ = 0: N_e ≈ 85 ± 10 beads (Everaers 2004, Hoy 2009)", {},
       "needs an equilibrated Kremer–Grest melt (millions of τ): build it with caps.build.kremer_grest, equilibrate on a cluster, then caps ppa --method caps|lammps"},
      {"T17", "Strain hardening of a Kremer–Grest glass", "T = 0.2: G_R rises with N and saturates for N ≳ 10 N_e (Hoy & Robbins 2006)", {},
       "needs equilibrated glasses of several chain lengths and their tension runs: caps mech decks, then caps mech analyze on each"},
      {"T18", "API regressions", "polymer(forcefield=…) then export; one set of charges by both routes", {"Check", "Result", "Status"}, nullptr},
  };
  return m;
}

const Meta* meta(const std::string& id) {
  for (const auto& m : metas())
    if (m.id == id) return &m;
  return nullptr;
}

// ---- shared structures: the PS melt from the samples, relaxed once per process
std::mutex g_cache_mutex;
std::map<std::string, System> g_cache;

System load_sample(const BenchOptions& o) {
  const std::string path = o.samples + "/ps_melt.data";
  System s = read_lammps_data(path);
  if (s.atoms.empty()) throw std::runtime_error("cannot read " + path);
  return s;
}

System relaxed_sample(const BenchOptions& o) {
  std::lock_guard<std::mutex> l(g_cache_mutex);
  auto it = g_cache.find(o.samples);
  if (it != g_cache.end()) return it->second;
  System s = load_sample(o);
  RelaxOptions r;
  r.ftol = 1.0;
  r.max_iterations = 3000;
  r.target_density = 0.9;   // a melt, not the grown cell's 0.39 g/cm³
  relax(s, r);
  s.velocities.clear();
  g_cache[o.samples] = s;
  return s;
}

System replicate(const System& s, int n) {
  System r = s;
  r.atoms.clear();
  r.bonds.clear();
  r.velocities.clear();
  const size_t na = s.atoms.size();
  int64_t mols = 0;
  for (const auto& a : s.atoms) mols = std::max(mols, a.mol);
  int k = 0;
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      for (int l = 0; l < n; ++l, ++k) {
        const Vec3 shift = s.cell.a * double(i) + s.cell.b * double(j) + s.cell.c * double(l);
        for (const auto& a : s.atoms) {
          Atom b = a;
          b.pos = a.pos + shift;
          b.id = int64_t(r.atoms.size() + 1);
          b.mol = a.mol + mols * k;
          r.atoms.push_back(b);
        }
        for (const auto& b : s.bonds) r.bonds.push_back({uint32_t(b.i + na * size_t(k)), uint32_t(b.j + na * size_t(k)), b.order});
      }
  r.cell.a = s.cell.a * double(n);
  r.cell.b = s.cell.b * double(n);
  r.cell.c = s.cell.c * double(n);
  return r;
}

std::vector<double> flat(const System& s) {
  std::vector<double> x;
  x.reserve(3 * s.atoms.size());
  for (const auto& a : s.atoms) x.insert(x.end(), {a.pos[0], a.pos[1], a.pos[2]});
  return x;
}

bool tick(const BenchOptions& o, const std::string& id, const std::string& what, double f) {
  return !o.progress || o.progress(id, what, f);
}

struct Cancelled : std::runtime_error {
  Cancelled() : std::runtime_error("bench cancelled") {}
};

void step(const BenchOptions& o, const std::string& id, const std::string& what, double f) {
  if (!tick(o, id, what, f)) throw Cancelled();
}

// ---- T1: forces and virial against finite differences of the energy (the pair set frozen at the start, so the
// energy is a smooth function and the difference quotient converges)
void t1(BenchTable& t, const BenchOptions& o) {
  System s = relaxed_sample(o);
  const ForceField ff = assign_gaff(s);
  const size_t n = s.atoms.size();
  for (int pass = 0; pass < 2; ++pass) {
    step(o, "T1", pass == 0 ? "forces, all terms" : "forces, Coulomb off", 0.1 + 0.4 * pass);
    EnergyOptions eo;
    eo.coulomb = pass == 0;
    Evaluator ev(ff, eo);
    std::vector<double> x = flat(s), f, g;
    ev.compute(x, s.cell, f);
    ev.freeze_pairs(x, s.cell);
    double rms = 0;
    for (double v : f) rms += v * v;
    rms = std::sqrt(rms / double(f.size()));
    const double h = 1e-4;
    double worst = 0;
    const size_t stride = std::max<size_t>(1, n / 16);
    for (size_t i = 0; i < n; i += stride)
      for (int c = 0; c < 3; ++c) {
        const size_t k = 3 * i + size_t(c);
        const double x0 = x[k];
        x[k] = x0 + h;
        const double ep = ev.compute(x, s.cell, g).total();
        x[k] = x0 - h;
        const double em = ev.compute(x, s.cell, g).total();
        x[k] = x0;
        worst = std::max(worst, std::fabs(-(ep - em) / (2 * h) - f[k]));
      }
    const bool ok = worst < 1e-3;
    t.rows.push_back({{pass == 0 ? "Forces · bonded, LJ + tail, DSF Coulomb" : "Forces · bonded, LJ + tail", thousands(long(n)), num(worst, 2), fmt("%.1f", rms),
                       ok ? "pass" : "fail"},
                      ok ? "pass" : "fail"});
  }
  // virial: dE/dλ under affine scaling r → λr of positions and cell equals −Σ r·f (bonded, pair and tail terms)
  step(o, "T1", "virial", 0.9);
  EnergyOptions eo;
  eo.tail = false;   // the tail pressure is the homogeneous-fluid formula (as LAMMPS), not −dE_tail/dλ at a fixed pair set
  Evaluator ev(ff, eo);
  std::vector<double> x = flat(s), f;
  const EnergyTerms e0 = ev.compute(x, s.cell, f);
  ev.freeze_pairs(x, s.cell);
  const double eps = 1e-5;
  auto at = [&](double lam) {
    std::vector<double> y = x, g;
    for (double& v : y) v *= lam;
    Cell c = s.cell;
    c.a = c.a * lam, c.b = c.b * lam, c.c = c.c * lam, c.origin = c.origin * lam;
    return ev.compute(y, c, g).total();
  };
  const double dEdl = (at(1 + eps) - at(1 - eps)) / (2 * eps);
  const double diff = std::fabs(dEdl + e0.virial);
  const bool ok = diff < 1e-4 * std::max(1.0, std::fabs(e0.virial));
  t.rows.push_back({{"Virial · Σ r·f vs dE/dλ (no tail)", thousands(long(n)), num(diff, 2) + " kcal/mol", num(std::fabs(e0.virial), 4) + " kcal/mol", ok ? "pass" : "fail"},
                    ok ? "pass" : "fail"});
  // PME's reciprocal part against the plain Ewald sum (both in Fortran), with the default grid
  step(o, "T1", "PME vs Ewald", 0.95);
  {
    std::vector<double> xs = flat(s);
    const double beta = ewald_beta(10.0, 1e-5);
    const PmeGrid g = pme_grid(s.cell, beta, 1.0, 5);
    std::vector<double> fp(xs.size(), 0.0), fe(xs.size(), 0.0);
    double vp[6], ve[6];
    const int km = int(std::ceil(0.55 * std::max({norm(s.cell.a), norm(s.cell.b), norm(s.cell.c)})));
    const int kmax[3] = {km, km, km};
    const double ep = pme_reciprocal(xs, ff.charge, s.cell, g, fp, vp), ee = ewald_reciprocal(xs, ff.charge, s.cell, beta, kmax, fe, ve);
    double rms = 0, err = 0;
    for (size_t k = 0; k < xs.size(); ++k) rms += fe[k] * fe[k], err += (fp[k] - fe[k]) * (fp[k] - fe[k]);
    const double rel = std::sqrt(err / std::max(rms, 1e-30));
    const bool okp = rel < 1e-3 && std::fabs(ep - ee) < 1e-3 * std::max(1.0, std::fabs(ee));
    t.rows.push_back({{"PME reciprocal vs Ewald sum (" + std::to_string(g.k[0]) + "³ grid, order 5)", thousands(long(n)), num(rel, 2) + " (relative)",
                       fmt("%.4f", ee) + " kcal/mol", okp ? "pass" : "fail"},
                      okp ? "pass" : "fail"});
  }
  t.note = "Central differences (h = 1e-4 Å, λ = 1 ± 1e-5) of the total energy with the pair list frozen; pass: |ΔF| < 1e-3 kcal/mol/Å, virial within 1e-4 relative. "
           "The LJ tail pressure is the homogeneous-fluid correction, which counts pairs crossing the cut-off and so is not the derivative "
           "at a frozen pair list; it and the full virial tensor are compared with LAMMPS in bench/ff/check_data_lammps.py.";
}

// ---- T2: Pack on water boxes
void t2(BenchTable& t, const BenchOptions& o) {
  Trajectory w = read_pdb(o.samples + "/water.pdb");
  System mol = w.topology;
  const std::vector<int> counts = o.quick ? std::vector<int>{300, 1000} : std::vector<int>{1000, 5000, 20000};
  for (size_t k = 0; k < counts.size(); ++k) {
    const int nw = counts[k];
    // edge for 1.0 g/cm³: V = N · 18.015 / 0.6022 Å³
    const double edge = std::cbrt(nw * 18.015 / 0.60221);
    std::vector<double> times, dmins;
    int ok = 0, runs = std::max(1, nw >= 20000 ? 1 : o.repeats);
    for (int r = 0; r < runs; ++r) {
      step(o, "T2", "water × " + std::to_string(nw), (double(k) + double(r) / runs) / double(counts.size()));
      PackItem it;
      it.name = "water";
      it.molecule = mol;
      it.count = nw;
      Region box;
      box.a = {0, 0, 0};
      box.b = {edge, edge, edge};
      it.regions.push_back(box);
      PackOptions po;
      po.tolerance = 2.0;
      po.cell.a = {edge, 0, 0};
      po.cell.b = {0, edge, 0};
      po.cell.c = {0, 0, edge};
      po.seed = uint64_t(r + 1);
      PackReport rep;
      const auto t0 = Clock::now();
      try {
        pack({it}, po, &rep);
      } catch (const std::exception&) {
      }
      times.push_back(since(t0));
      dmins.push_back(rep.dmin);
      ok += rep.success;
    }
    t.rows.push_back({{"Water box, 1.0 g/cm³", thousands(nw), thousands(3L * nw), mean_sd(times, "%.2f"), mean_sd(dmins, "%.3f"),
                       std::to_string(ok) + " / " + std::to_string(runs)},
                      ok == runs ? "pass" : "fail"});
  }
  t.note = "Tolerance 2.0 Å between molecules, periodic cell; a run passes when no pair is closer than the tolerance. "
           "Against packmol: bench/pack/run_bench.py (needs packmol).";
}

// ---- T4: throughput
void t4(BenchTable& t, const BenchOptions& o) {
  const System base = relaxed_sample(o);
  const int threads = max_threads() > 0 ? max_threads() : std::min(16, int(std::max(1u, std::thread::hardware_concurrency())));
  for (int big = 0; big < 2; ++big) {
    System s = big ? replicate(base, 2) : base;
    std::vector<double> nsd;
    for (int r = 0; r < std::max(1, o.repeats); ++r) {
      step(o, "T4", big ? "10k atoms" : "1.3k atoms", (big + double(r) / o.repeats) / 2);
      System run = s;
      DynamicsOptions d;
      d.steps = o.quick ? 300 : 2000;
      d.thermo_every = int(d.steps);
      d.frame_every = 0;
      d.seed = uint64_t(r + 1);
      DynamicsReport rep;
      run_dynamics(run, d, &rep);
      nsd.push_back(rep.ns_per_day);
    }
    t.rows.push_back({{"PS melt" + std::string(big ? " × 8" : ""), thousands(long(s.atoms.size())), std::to_string(threads), mean_sd(nsd, "%.1f"), "info"}, "info"});
  }
  {
    // the same melt with particle-mesh Ewald
    std::vector<double> nsd;
    for (int r = 0; r < std::max(1, o.repeats); ++r) {
      step(o, "T4", "1.3k atoms, PME", 0.9);
      System run = base;
      DynamicsOptions d;
      d.steps = o.quick ? 300 : 2000;
      d.thermo_every = int(d.steps);
      d.frame_every = 0;
      d.seed = uint64_t(r + 1);
      d.energy.electrostatics = EnergyOptions::Electrostatics::PME;
      DynamicsReport rep;
      run_dynamics(run, d, &rep);
      nsd.push_back(rep.ns_per_day);
    }
    t.rows.push_back({{"PS melt, PME (1 Å grid, order 5)", thousands(long(base.atoms.size())), std::to_string(threads), mean_sd(nsd, "%.1f"), "info"}, "info"});
  }
  t.note = "Built-in GAFF, 10 Å cut-off, DSF Coulomb (PME where named), dt 1 fs, Bussi thermostat; mean ± sd over the repeats.";
}

// ---- T5: NVE drift
void t5(BenchTable& t, const BenchOptions& o) {
  System s = relaxed_sample(o);
  const double kT = 0.0019872043 * 300.0;
  // thermalise at 300 K first
  step(o, "T5", "thermalising", 0.05);
  {
    DynamicsOptions d;
    d.steps = o.quick ? 500 : 2000;
    d.frame_every = 0;
    d.thermo_every = int(d.steps);
    run_dynamics(s, d);
  }
  const double dts[] = {1.0, 0.5};
  for (int k = 0; k < 2; ++k) {
    step(o, "T5", "NVE dt " + fmt("%.1f", dts[k]) + " fs", 0.2 + 0.4 * k);
    System run = s;
    DynamicsOptions d;
    d.thermostat = Thermostat::None;
    d.dt = dts[k];
    const double ps = o.quick ? 2.0 : 10.0;
    d.steps = int64_t(ps * 1000.0 / d.dt);
    d.thermo_every = 10;
    d.frame_every = 0;
    DynamicsReport rep;
    run_dynamics(run, d, &rep);
    // linear fit of the conserved energy against time
    double sx = 0, sy = 0, sxx = 0, sxy = 0;
    const double m = double(rep.thermo.size());
    for (const auto& r : rep.thermo) sx += r.time_ps, sy += r.conserved, sxx += r.time_ps * r.time_ps, sxy += r.time_ps * r.conserved;
    const double slope = (m * sxy - sx * sy) / std::max(1e-30, m * sxx - sx * sx);   // kcal/mol/ps
    const double icpt = (sy - slope * sx) / m;
    double res = 0;
    for (const auto& r : rep.thermo) res += std::pow(r.conserved - (icpt + slope * r.time_ps), 2);
    const double drift = slope * 1000.0 / double(run.atoms.size()) / kT;   // kT/ns/atom
    // GROMACS's default verlet-buffer-tolerance: 0.005 kJ/mol/ps per atom = 2.0 kT/ns/atom at 300 K
    const double limit = 0.005 / 4.184 * 1000.0 / kT;
    const bool ok = std::fabs(drift) < limit;
    t.rows.push_back({{"PS melt, 300 K", fmt("%.1f", dts[k]), fmt("%.0f", ps), fmt("%+.4f", drift), fmt("%.3f", std::sqrt(res / std::max(1.0, m))), ok ? "pass" : "fail"},
                      ok ? "pass" : "fail"});
  }
  t.note = "Velocity Verlet without a thermostat after thermalising with Bussi; drift is the slope of a straight-line fit to the total energy. "
           "Pass: below GROMACS's default tolerance (verlet-buffer-tolerance 0.005 kJ/mol/ps per atom, 2.0 kT/ns/atom at 300 K).";
}

// ---- T8: thread scaling
void t8(BenchTable& t, const BenchOptions& o) {
  const System s = replicate(relaxed_sample(o), 2);
  const ForceField ff = assign_gaff(s);
  const std::vector<double> x = flat(s);
  const int hw = int(std::max(1u, std::thread::hardware_concurrency()));
  std::vector<int> ts{1, 2, 4, 8, 16};
  ts.erase(std::remove_if(ts.begin(), ts.end(), [&](int v) { return v > hw; }), ts.end());
  if (std::find(ts.begin(), ts.end(), std::min(hw, 16)) == ts.end()) ts.push_back(std::min(hw, 16));
  double t1ms = 0;
  for (size_t k = 0; k < ts.size(); ++k) {
    step(o, "T8", std::to_string(ts[k]) + " threads", double(k) / double(ts.size()));
    EnergyOptions eo;
    eo.threads = ts[k];
    Evaluator ev(ff, eo);
    std::vector<double> f;
    ev.compute(x, s.cell, f);
    ev.compute(x, s.cell, f);
    const int reps = o.quick ? 5 : 20;
    const auto t0 = Clock::now();
    for (int r = 0; r < reps; ++r) ev.compute(x, s.cell, f);
    const double ms = since(t0) * 1000.0 / reps;
    if (k == 0) t1ms = ms;
    t.rows.push_back({{"PS melt × 8", std::to_string(ts[k]), fmt("%.2f", ms), fmt("%.2f×", t1ms / ms), fmt("%.0f %%", 100.0 * t1ms / ms / ts[k])}, "info"});
  }
  t.note = "Energy and forces for the whole system (pair list reused), after two warm-up evaluations.";
}

// ---- T9: reproducibility
void t9(BenchTable& t, const BenchOptions& o) {
  const System s0 = relaxed_sample(o);
  auto run = [&](int threads) {
    System s = s0;
    DynamicsOptions d;
    d.steps = o.quick ? 200 : 1000;
    d.seed = 7;
    d.frame_every = 0;
    d.thermo_every = int(d.steps);
    d.energy.threads = threads;
    run_dynamics(s, d);
    return s;
  };
  auto maxdiff = [](const System& a, const System& b) {
    double m = 0;
    for (size_t i = 0; i < a.atoms.size(); ++i) m = std::max(m, norm(a.atoms[i].pos - b.atoms[i].pos));
    return m;
  };
  const int hw = std::min(16, int(std::max(1u, std::thread::hardware_concurrency())));
  step(o, "T9", "first run", 0.1);
  const System a = run(hw);
  step(o, "T9", "second run", 0.4);
  const System b = run(hw);
  const double d1 = maxdiff(a, b);
  const int64_t steps = o.quick ? 200 : 1000;
  t.rows.push_back({{"Same seed, same threads (" + std::to_string(hw) + ")", std::to_string(steps), d1 == 0 ? "0 (bit-identical)" : num(d1, 2), d1 == 0 ? "pass" : "fail"},
                    d1 == 0 ? "pass" : "fail"});
  step(o, "T9", "one thread", 0.7);
  const System c = run(1);
  const double d2 = maxdiff(a, c);
  t.rows.push_back({{"Same seed, 1 vs " + std::to_string(hw) + " threads", std::to_string(steps), d2 == 0 ? "0 (bit-identical)" : num(d2, 2), "info"}, "info"});
  t.note = "NVT with Bussi, seed 7. Each worker owns a fixed share of the pairs and results are reduced in worker order, so a run repeats exactly "
           "for a given thread count; another thread count sums in another order and differs in the last bits, which dynamics amplifies.";
}

// ---- T6 / T7: the user's equilibrated cells against the reference ranges
struct RefCell { std::string id, name, file; Json values; };

std::vector<RefCell> reference_cells(const BenchOptions& o, std::vector<std::string>& missing) {
  std::vector<RefCell> r;
  if (o.reference.empty()) throw std::runtime_error("no reference file (data/reference/polymers.json)");
  std::ifstream in(o.reference);
  if (!in) throw std::runtime_error("cannot read " + o.reference);
  std::stringstream ss;
  ss << in.rdbuf();
  const Json ref = Json::parse(ss.str());
  for (const auto& m : ref["materials"].items()) {
    const std::string id = m.text("id", "");
    const std::string f = o.cells.empty() ? "" : o.cells + "/" + id + ".data";
    if (!f.empty() && std::filesystem::exists(f)) r.push_back({id, m.text("name", id), f, m["values"]});
    else missing.push_back(m.text("name", id) + " (" + id + ".data)");
  }
  return r;
}

Trajectory cell_frames(const RefCell& c) {
  const std::string stem = c.file.substr(0, c.file.size() - 5);
  for (const char* ext : {".lammpstrj", ".dump", ".dcd", ".xtc"})
    if (std::filesystem::exists(stem + ext)) return open_file(stem + ext, c.file);
  return open_file(c.file);
}

void t67(BenchTable& t, const BenchOptions& o, bool chains) {
  if (o.reference.empty() || !std::filesystem::exists(o.reference)) {
    t.note = "No reference ranges (data/reference/polymers.json) given to the suite: nothing to compare with.";
    throw Cancelled();   // reported as not run
  }
  std::vector<std::string> missing;
  const auto cells = reference_cells(o, missing);
  const std::string prop = chains ? "cn" : "density";
  size_t k = 0;
  for (const auto& c : cells) {
    step(o, t.id, c.name, double(k++) / std::max<size_t>(1, cells.size()));
    if (!c.values.has(prop)) {
      t.rows.push_back({{c.name, std::filesystem::path(c.file).filename().string(), "—", "—", "no reference value", "info"}, "info"});
      continue;
    }
    const Json& rv = c.values[prop];
    const double lo = rv.num("lo", 0), hi = rv.num("hi", 0);
    const Trajectory tr = cell_frames(c);
    AnalyzeOptions ao;
    ao.blocks = 5;
    const auto ps = analyze(tr, {prop}, ao);
    const Property& p = ps.front();
    const bool ok = std::isfinite(p.value) && p.value >= lo && p.value <= hi;
    std::string counted = std::to_string(tr.frames());
    if (chains) {
      int nm = 0;
      tr.topology.molecules(&nm);
      counted = std::to_string(nm);
    }
    const std::string val = std::isfinite(p.value) ? num(p.value, 4) + (std::isfinite(p.error) && p.error > 0 ? " ± " + num(p.error, 2) : "") : "—";
    t.rows.push_back({{c.name, std::filesystem::path(c.file).filename().string(), counted, val, num(lo, 4) + " – " + num(hi, 4) + " (" + rv.text("source", "") + ")",
                       std::isfinite(p.value) ? (ok ? "pass" : "fail") : "info"},
                      std::isfinite(p.value) ? (ok ? "pass" : "fail") : "info", c.file});
  }
  for (const auto& m : missing) t.rows.push_back({{m, "no cell", "—", "—", "—", "info"}, "info"});
  if (cells.empty()) {
    t.rows.clear();
    t.status = "not run";
    t.note = "No equilibrated cells in " + (o.cells.empty() ? std::string("the cells folder") : o.cells) + ": put a relaxed, equilibrated cell per material there, named by its "
             "reference id (" + [&] { std::string l; for (const auto& m : missing) l += (l.empty() ? "" : ", ") + m; return l; }() +
             "), with its trajectory of the same name for frame averages; CAPS then measures " + (chains ? "C∞ (Analyze cn)" : "the density") +
             " and compares it with the range in data/reference/polymers.json. Materials without reference values there are not judged (add them with their source).";
    throw Cancelled();   // reported as not run
  }
  t.note = chains ? "C∞ from the internal distances ⟨R²(n)⟩/(n⟨b²⟩) extrapolated in 1/n (Analyze cn), over the cell's chains and frames; pass when inside the literature range."
                  : "Density averaged over the cell's frames (block errors, Analyze density); pass when inside the experimental range. The cells are yours: their equilibration is what is judged.";
}

// ---- T11: rendering
void t11(BenchTable& t, const BenchOptions& o) {
  const System base = relaxed_sample(o);
  int k = 0;
  for (int big = 0; big < 2; ++big) {
    const System s = big ? replicate(base, 2) : base;
    for (int ss : {1, 2}) {
      step(o, "T11", "render", k++ / 4.0);
      Renderer r;
      RenderOptions ro;
      ro.width = 1920;
      ro.height = 1080;
      ro.supersample = ss;
      Camera cam;
      r.render(s, cam, ro);
      std::vector<double> ms;
      for (int rep = 0; rep < std::max(1, o.repeats); ++rep) {
        const auto t0 = Clock::now();
        r.render(s, cam, ro);
        ms.push_back(since(t0) * 1000.0);
      }
      t.rows.push_back({{big ? "PS melt × 8" : "PS melt", thousands(long(s.atoms.size())), std::to_string(ss) + "×", mean_sd(ms, "%.0f")}, "info"});
    }
  }
  t.note = "Ball and stick with outlines and depth cue, one thread for the rasteriser; mean ± sd over the repeats.";
}

// ---- T12: molecule builder
void t12(BenchTable& t, const BenchOptions& o) {
  // (S) at the alpha carbon of alanine, found by elements
  auto alanine_S = [](const BuildResult& r) {
    const auto& g = r.graph;
    const auto& p = r.conformers.front().pos;
    auto nb = [&](int i) {
      std::vector<int> v;
      for (const auto& b : g.bonds) {
        if (b.a == i) v.push_back(b.b);
        if (b.b == i) v.push_back(b.a);
      }
      return v;
    };
    int N = -1, CA = -1, CO = -1, CM = -1, H = -1;
    for (size_t i = 0; i < g.atoms.size(); ++i)
      if (g.atoms[i].element == 7) N = int(i);
    for (int v : nb(N))
      if (g.atoms[size_t(v)].element == 6) CA = v;
    for (int v : nb(CA)) {
      const int z = g.atoms[size_t(v)].element;
      if (z == 1) H = v;
      if (z != 6) continue;
      int ox = 0;
      for (int w : nb(v)) ox += g.atoms[size_t(w)].element == 8;
      (ox ? CO : CM) = v;
    }
    return dot(cross(p[size_t(CO)] - p[size_t(N)], p[size_t(CM)] - p[size_t(N)]), p[size_t(CA)] - p[size_t(H)]) > 0;
  };
  BuildOptions bo;
  if (!o.forcefields.empty()) bo.forcefield = o.forcefields + "/gaff-amber25.json";
  struct Case { const char* smiles; bool S; };
  const Case ala[] = {{"N[C@@H](C)C(=O)O", true}, {"C[C@@H](C(=O)O)N", true}, {"OC(=O)[C@@H](N)C", true}, {"N[C@H](C)C(=O)O", false},
                      {"C[C@H](C(=O)O)N", false}, {"[C@@H](N)(C)C(=O)O", false}};
  int ok = 0, n = 0;
  step(o, "T12", "configurations", 0.1);
  for (const auto& c : ala) {
    ++n;
    try {
      ok += alanine_S(build_molecule(c.smiles, bo)) == c.S;
    } catch (const std::exception&) {
    }
  }
  t.rows.push_back({{"Tetrahedral configuration (alanine, 6 spellings)", std::to_string(n), std::to_string(ok), ok == n ? "pass" : "fail"}, ok == n ? "pass" : "fail"});

  step(o, "T12", "E/Z", 0.3);
  struct EZ { const char* smiles; int a, b, c, d; bool trans; };
  const EZ ez[] = {{"F/C=C/F", 0, 1, 2, 3, true}, {"F/C=C\\F", 0, 1, 2, 3, false}, {"C(\\F)=C/F", 1, 0, 2, 3, true}, {"C(/F)=C/F", 1, 0, 2, 3, false},
                   {"C/C=C/C", 0, 1, 2, 3, true}, {"C/C=C\\C", 0, 1, 2, 3, false}};
  ok = n = 0;
  for (const auto& c : ez) {
    ++n;
    try {
      const auto r = build_molecule(c.smiles, bo);
      const auto& p = r.conformers.front().pos;
      const Vec3 b1 = p[size_t(c.b)] - p[size_t(c.a)], b2 = p[size_t(c.c)] - p[size_t(c.b)], b3 = p[size_t(c.d)] - p[size_t(c.c)];
      const Vec3 n1 = cross(b1, b2), n2 = cross(b2, b3);
      const double phi = std::fabs(std::atan2(dot(cross(n1, b2 * (1.0 / norm(b2))), n2), dot(n1, n2))) * 180 / 3.14159265358979;
      ok += c.trans ? phi > 170 : phi < 10;
    } catch (const std::exception&) {
    }
  }
  t.rows.push_back({{"Double-bond E/Z", std::to_string(n), std::to_string(ok), ok == n ? "pass" : "fail"}, ok == n ? "pass" : "fail"});

  step(o, "T12", "round trips", 0.5);
  const char* set[] = {"CC(=O)Oc1ccccc1C(=O)O", "C1CC2CCC1C2", "c1ccc2cc3ccccc3cc2c1", "[NH4+].[Cl-]", "F/C=C/F", "N[C@@H](C)C(=O)O",
                       "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O", "C[C@]12CC[C@H]3[C@@H](CCC4=CC(=O)CC[C@@]34C)[C@@H]1CC[C@@H]2O",
                       "CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "O=C(O)CCCCCCC/C=C\\CCCCCCCC"};
  ok = n = 0;
  for (const char* smi : set) {
    ++n;
    try {
      MolGraph g = parse_smiles(smi);
      const std::string w = write_smiles(g);
      MolGraph h = parse_smiles(w);
      const bool same_text = write_smiles(h) == w;
      add_hydrogens(g);
      add_hydrogens(h);
      const MolInfo a = molecule_info(g), b = molecule_info(h);
      ok += same_text && a.formula == b.formula && a.stereocentres == b.stereocentres && a.stereo_bonds == b.stereo_bonds;
    } catch (const std::exception&) {
    }
  }
  t.rows.push_back({{"SMILES written and read back", std::to_string(n), std::to_string(ok), ok == n ? "pass" : "fail"}, ok == n ? "pass" : "fail"});

  step(o, "T12", "embedding", 0.7);
  ok = n = 0;
  int kept = 0;
  for (const char* smi : set) {
    ++n;
    try {
      const auto r = build_molecule(smi, bo);
      ++ok;
      const auto c = chirality_check(r.graph, r.conformers.front().pos);
      kept += std::all_of(c.begin(), c.end(), [](int v) { return v > 0; });
    } catch (const std::exception&) {
    }
  }
  t.rows.push_back({{std::string("3D embedding") + (bo.forcefield.empty() ? "" : " + GAFF2 clean-up"), std::to_string(n), std::to_string(ok), ok == n ? "pass" : "fail"},
                    ok == n ? "pass" : "fail"});
  t.rows.push_back({{"Stereocentres kept after the clean-up", std::to_string(ok), std::to_string(kept), kept == ok ? "pass" : "fail"}, kept == ok ? "pass" : "fail"});
  t.note = "Configurations are checked geometrically from the built coordinates (no CIP code involved), dihedrals within 10° of 0 or 180°.";
}


// ---- T13: chemistry-aware mapping of grown copolyesters (ester cut: diol B, diacid S, A, T)
CgRules ester_cut_rules() {
  CgRules r;
  r.cut = {"[CX3](=O)-[OX2;!H1]"};
  r.names = {{"T", "[CX3](=O)c1ccc(cc1)[CX3]=O"}, {"A", "[CX3](=O)[CH2][CH2][CH2][CH2][CX3]=O"}, {"S", "[CX3](=O)[CH2][CH2][CX3]=O"}, {"B", "[OX2][CH2][CH2][CH2][CH2][OX2]"}};
  return r;
}
void t13(BenchTable& t, const BenchOptions& o) {
  struct Case { const char* name; std::vector<std::string> units; Sequence seq; std::string pattern; std::vector<double> w; };
  const std::vector<Case> cases = {{"PBS", {"[*]OCCCCOC(=O)CCC(=O)[*]"}, Sequence::Homopolymer, "", {}},
                                   {"PBSA 80/20", {"[*]OCCCCOC(=O)CCC(=O)[*]", "[*]OCCCCOC(=O)CCCCC(=O)[*]"}, Sequence::Pattern, "AAAAB", {}},
                                   {"PBAT 56/44", {"[*]OCCCCOC(=O)CCCCC(=O)[*]", "[*]OCCCCOC(=O)c1ccc(cc1)C(=O)[*]"}, Sequence::Random, "", {0.56, 0.44}}};
  int k = 0;
  for (const auto& c : cases) {
    step(o, "T13", c.name, double(k++) / cases.size());
    ChainSpec spec;
    for (const auto& u : c.units) spec.units.push_back({u, u});
    spec.sequence = c.seq, spec.pattern = c.pattern, spec.weights = c.w, spec.dp = o.quick ? 10 : 25, spec.tail_cap = "hydroxyl";
    GrowOptions g;
    g.chains = o.quick ? 4 : 10, g.density = 0.1, g.seed = 11, g.auto_scale = true;
    const System aa = grow_chains(spec, g);
    const CgMapping m = cg_mapping(aa, ester_cut_rules());
    const auto sum = cg_mapping_summary(m);
    const double per_unit = double(m.beads()) / double(spec.dp * g.chains);
    std::string comp;
    for (const auto& [kk, f] : sum.fraction) { char b[40]; std::snprintf(b, sizeof b, "%s%s %.1f%%", comp.empty() ? "" : " ", kk.c_str(), 100 * f); comp += b; }
    // the composition the chains were grown with: units counted from the structure (aromatic rings, adipate beads)
    bool comp_ok = true;
    if (std::string(c.name).rfind("PBSA", 0) == 0) comp_ok = sum.kinds.at("S") * 1 == 4 * sum.kinds.at("A");
    if (std::string(c.name).rfind("PBAT", 0) == 0) {
      const Perception p = perceive(aa);
      int ring = 0;
      for (size_t i = 0; i < aa.atoms.size(); ++i) ring += p.aromatic[i] ? 1 : 0;
      comp_ok = sum.kinds.count("T") && sum.kinds.at("T") * 6 == ring;
    }
    const bool ok = std::fabs(per_unit - 2) < 1e-9 && sum.mass_error < 1e-6 && comp_ok;
    char me[32];
    std::snprintf(me, sizeof me, "%.1e", sum.mass_error);
    t.rows.push_back({{c.name, std::to_string(m.beads()), std::to_string(per_unit).substr(0, 4), comp, me, ok ? "pass" : "fail"}, ok ? "pass" : "fail"});
  }
  t.note = "Grown with CAPS Grow (acid ends), mapped by the ester-cut rules (data/cg/mapping_rules.json): two beads per repeat unit, mass conserved, PBSA's AAAAB pattern 80/20 S/A exactly, PBAT's T beads one per terephthalate ring.";
}

// ---- T14: Boltzmann inversion of an ideal chain drawn exactly from known potentials
void t14(BenchTable& t, const BenchOptions& o) {
  constexpr double kB = 0.0019872067, T = 300, kT = kB * T, kPi = 3.14159265358979323846;
  auto Ub = [](double r) { return 3.0 * (r - 5) * (r - 5) + 0.6 * (r - 5) * (r - 5) * (r - 5); };
  auto Ua = [](double a) { return 0.0012 * (a - 130) * (a - 130); };
  auto Ud = [](double p) { const double x = p * kPi / 180; return 0.8 * (1 + std::cos(x)) + 0.35 * (1 - std::cos(2 * x)) + 0.15 * std::sin(x); };
  std::mt19937_64 g(11);
  auto unit = [&]() { return double(g() >> 11) * (1.0 / 9007199254740992.0); };
  struct Draw { std::vector<double> x, c; };
  auto make = [&](const std::function<double(double)>& w, double lo, double hi) {
    Draw d;
    const int n = 200000;
    double c = 0;
    for (int i = 0; i < n; ++i) { const double x = lo + (hi - lo) * (i + 0.5) / n; c += w(x); d.x.push_back(x), d.c.push_back(c); }
    for (auto& v : d.c) v /= c;
    return d;
  };
  auto draw = [&](const Draw& d) {
    const size_t i = std::min(d.x.size() - 1, size_t(std::lower_bound(d.c.begin(), d.c.end(), unit()) - d.c.begin()));
    return d.x[i] + (unit() - 0.5) * (d.x[1] - d.x[0]);
  };
  const Draw rb = make([&](double r) { return r * r * std::exp(-Ub(r) / kT); }, 2, 8), ta = make([&](double a) { return std::sin(a * kPi / 180) * std::exp(-Ua(a) / kT); }, 0, 180),
             pd = make([&](double p) { return std::exp(-Ud(p) / kT); }, -180, 180);
  auto place = [](const Vec3& a, const Vec3& b, const Vec3& c, double r, double th, double ph) {
    const double tt = th * kPi / 180, pp = ph * kPi / 180;
    Vec3 bc = c - b;
    bc = bc * (1.0 / norm(bc));
    Vec3 n = cross(b - a, bc);
    n = n * (1.0 / norm(n));
    const Vec3 m = cross(n, bc);
    return c + bc * (-r * std::cos(tt)) + m * (r * std::sin(tt) * std::cos(pp)) + n * (r * std::sin(tt) * std::sin(pp));
  };
  CgTopology top;
  const int mols = 400, frames = o.quick ? 200 : 1000;
  for (int m = 0; m < mols; ++m) {
    const int b = 4 * m;
    for (int k = 0; k < 4; ++k) top.kind.push_back(k % 2 ? "B" : "A"), top.mol.push_back(m);
    top.bonds.push_back({b, b + 1}), top.bonds.push_back({b + 1, b + 2}), top.bonds.push_back({b + 2, b + 3});
    top.angles.push_back({b, b + 1, b + 2}), top.angles.push_back({b + 1, b + 2, b + 3});
    top.dihedrals.push_back({b, b + 1, b + 2, b + 3});
  }
  const CgTypes types{{"A", "B"}, {"A-B"}, {"A-B-A", "B-A-B"}, {"A-B-A-B"}};
  CgBondedOptions bo;
  CgBondedAccumulator acc(types, bo);
  const int si = acc.add_system(top);
  for (int f = 0; f < frames; ++f) {
    if (f % 50 == 0) step(o, "T14", "drawing the chains", double(f) / frames);
    std::vector<Vec3> pos;
    for (int m = 0; m < mols; ++m) {
      const double r1 = draw(rb), r2 = draw(rb), r3 = draw(rb), t1 = draw(ta), t2 = draw(ta), ph = draw(pd);
      const Vec3 a{0, 0, 0}, b{r1, 0, 0};
      const Vec3 c = b + Vec3{-r2 * std::cos(t1 * kPi / 180), r2 * std::sin(t1 * kPi / 180), 0};
      const Vec3 d = place(a, b, c, r3, t2, ph);
      const Vec3 off{double(m % 20) * 30.0, double(m / 20) * 30.0, 0};
      pos.push_back(a + off), pos.push_back(b + off), pos.push_back(c + off), pos.push_back(d + off);
    }
    acc.add_frame(si, pos, Cell{});
  }
  const CgBondedResult r = invert_bonded(acc);
  auto worst = [&](const CgBondedTable& tb, const std::function<double(double)>& U) {
    double umin = 1e300;
    for (double x : tb.x) umin = std::min(umin, U(x));
    std::vector<std::pair<double, double>> v;
    for (size_t i = 0; i < tb.x.size(); ++i)
      if (U(tb.x[i]) - umin < 2.5 * kT && tb.x[i] >= tb.lo && tb.x[i] <= tb.hi) v.push_back({tb.U[i], U(tb.x[i])});
    double mean = 0;
    for (const auto& [a, b] : v) mean += a - b;
    mean /= double(v.size());
    double w = 0;
    for (const auto& [a, b] : v) w = std::max(w, std::fabs(a - b - mean) / kT);
    return w;
  };
  auto row = [&](const std::string& name, const CgBondedTable& tb, const std::function<double(double)>& U) {
    const double w = worst(tb, U);
    char b[32];
    std::snprintf(b, sizeof b, "%.3f", w);
    t.rows.push_back({{name, std::to_string(tb.count), b, w < 0.1 ? "pass" : "fail"}, w < 0.1 ? "pass" : "fail"});
  };
  row("bond A-B", r.bonds[0], Ub);
  row("angle A-B-A", r.angles[0], Ua);
  row("angle B-A-B", r.angles[1], Ua);
  row("dihedral A-B-A-B", r.dihedrals[0], Ud);
  t.note = "Internal coordinates drawn independently from r² e^(−U/kT), sin θ e^(−U/kT) and e^(−U/kT) (an ideal chain without non-bonded terms), inverted by caps cgfit bonded; compared where U < 2.5 kT above its minimum. Pass: within 0.1 kT.";
}

// ---- T15: IBI self-consistency on a two-type Lennard-Jones fluid with CAPS's engine
void t15(BenchTable& t, const BenchOptions& o) {
  const CgTypes types{{"A", "B"}, {}, {}, {}};
  CgTopology top;
  System s;
  const int n = o.quick ? 500 : 1000;
  const double L = std::cbrt(double(n) * 42.0 * 42.0 * 42.0 / 500.0);
  s.cell.a = {L, 0, 0}, s.cell.b = {0, L, 0}, s.cell.c = {0, 0, L};
  const int m = int(std::ceil(std::cbrt(double(n))));
  for (int i = 0; i < n; ++i) {
    Atom a;
    a.id = i + 1, a.mol = i + 1, a.name = i % 2 ? "B" : "A", a.type = i % 2 + 1;
    a.pos = {(i % m + 0.5) * L / m, ((i / m) % m + 0.5) * L / m, (i / (m * m) + 0.5) * L / m};
    s.atoms.push_back(a);
    top.kind.push_back(a.name), top.mol.push_back(i);
  }
  s.types = {{1, 50, "A"}, {2, 60, "B"}};
  const std::map<std::string, double> mass = {{"A", 50}, {"B", 60}};
  CgPairSet grid;
  grid.r0 = 1.0, grid.dr = 0.05, grid.rc = 11.0, grid.temperature = 300;
  const CgPairSet ref = tables_of_fits({{"A-A", "lj126", 0.45, 4.2}, {"A-B", "lj126", 0.50, 4.4}, {"B-B", "lj126", 0.55, 4.6}}, grid);
  CgEngineOptions eo;
  eo.steps = o.quick ? 2000 : 10000, eo.equilibrate = o.quick ? 1000 : 4000, eo.frame_every = 50, eo.seed = 7;
  step(o, "T15", "the reference fluid", 0);
  const CgEngineRun r0 = run_cg_engine(s, cg_forcefield(top, types, mass, ref), eo);
  CgPairOptions po;
  po.rmax = 11.0;
  CgRdfAccumulator acc(types, po);
  const int si = acc.add_system(top);
  for (size_t f = 0; f < r0.frames.size(); ++f) acc.add_frame(si, r0.frames[f], r0.cells[f]);
  const CgTargets tg = cg_targets(acc, 300, {r0.pressure});
  CgIbiOptions io;
  io.rc = 11.0, io.alpha = 0.3;
  CgPairSet p = ibi_start(tg, io);
  System cur = r0.last;
  const int iters = o.quick ? 4 : 20;
  CgEngineOptions it = eo;
  it.new_velocities = false, it.steps = o.quick ? 1000 : 5000, it.equilibrate = o.quick ? 500 : 1000;
  double err = 0, dp = 0;
  for (int k = 0; k < iters; ++k) {
    step(o, "T15", "IBI iteration " + std::to_string(k + 1), double(k + 1) / (iters + 1));
    const CgEngineRun run = run_cg_engine(cur, cg_forcefield(top, types, mass, p), it);
    CgRdfAccumulator a2(types, po);
    const int s2 = a2.add_system(top);
    for (size_t f = 0; f < run.frames.size(); ++f) a2.add_frame(s2, run.frames[f], run.cells[f]);
    CgIbiStepReport rep;
    p = ibi_step(p, tg, {a2.rdf(s2)}, {run.pressure}, io, &rep);
    err = 100 * std::sqrt(rep.residual), dp = std::fabs(run.pressure - r0.pressure);
    char e[32], d[32];
    std::snprintf(e, sizeof e, "%.2f", err), std::snprintf(d, sizeof d, "%.0f", dp);
    const bool last = k + 1 == iters;
    const bool ok = err < 1.0 && dp < 50;
    t.rows.push_back({{std::to_string(k + 1), e, d, last ? (ok ? "pass" : o.quick ? "info" : "fail") : "info"}, last ? (ok ? "pass" : o.quick ? "info" : "fail") : "info"});
    cur = run.last;
  }
  t.note = "Target: the per-pair g(r) of a fluid with known Lennard-Jones pairs (" + std::to_string(n) + " beads, 300 K, CAPS's engine); IBI from −kT ln g with the pressure ramp; RDF error = √(∫(g − g_t)² dr / ∫ g_t² dr), the largest over the pairs. Pass: below 1 % with |ΔP| < 50 atm at the last iteration" +
           std::string(o.quick ? " (quick: 4 short iterations, informative only)" : "") + ".";
}

// ---- T18: API regressions (polymer with a force field, then export; one set of charges)
void t18(BenchTable& t, const BenchOptions& o) {
  step(o, "T18", "recipe", 0);
  Json r = Json::parse(R"({"recipe": 1, "name": "pbs", "build": {"polymer": {"units": ["[*]OCCCCOC(=O)CCC(=O)[*]"], "dp": 3, "chains": 1, "tacticity": "atactic", "sequence": "homopolymer"}},
                         "grow": {"density": 0.1, "seed": 1}, "type": {"forcefield": "opls2005"}})");
  RecipeOptions ro;
  ro.forcefield_dir = o.forcefields;
  const auto res = run_recipe(r, ro);
  const bool has_ff = res.field && res.field->charge.size() == res.system.atoms.size();
  t.rows.push_back({{"polymer(forcefield = opls2005) keeps its force field", has_ff ? res.forcefield : "none", has_ff ? "pass" : "fail"}, has_ff ? "pass" : "fail"});
  double dq = has_ff ? 0 : 1;
  if (has_ff)
    for (size_t i = 0; i < res.system.atoms.size(); ++i) dq = std::max(dq, std::fabs(res.system.atoms[i].charge - res.field->charge[i]));
  char b[48];
  std::snprintf(b, sizeof b, "max |Δq| %.1e e", dq);
  t.rows.push_back({{"the structure carries the force field's charges", b, dq < 1e-9 ? "pass" : "fail"}, dq < 1e-9 ? "pass" : "fail"});
  std::string label;
  for (const auto& st : res.manifest.steps)
    if (st.engine == "field.assign") label = st.summary;
  const bool named = label.find("bond increments") != std::string::npos;
  t.rows.push_back({{"OPLS 2005 charges named as bond increments", label, named ? "pass" : "fail"}, named ? "pass" : "fail"});
  t.note = "The issues reported with the PBS/PBSA/PBAT cells: export after polymer(forcefield=…) asked for a force field; the recipe and Field described one set of charges in two ways.";
}
}  // namespace

std::vector<std::string> bench_ids() {
  std::vector<std::string> v;
  for (const auto& m : metas()) v.push_back(m.id);
  return v;
}

BenchTable bench_describe(const std::string& id) {
  const Meta* m = meta(id);
  if (!m) throw std::runtime_error("no bench table " + id);
  BenchTable t;
  t.id = m->id;
  t.title = m->title;
  t.scope = m->scope;
  t.columns = m->columns;
  t.status = "not run";
  if (m->not_run) t.note = std::string("Not part of the built-in suite: ") + m->not_run + ".";
  return t;
}

BenchTable run_bench(const std::string& id, const BenchOptions& o) {
  BenchTable t = bench_describe(id);
  const Meta* m = meta(id);
  if (m->not_run) return t;
  const auto t0 = Clock::now();
  try {
    if (id == "T1") t1(t, o);
    else if (id == "T2") t2(t, o);
    else if (id == "T4") t4(t, o);
    else if (id == "T5") t5(t, o);
    else if (id == "T6") t67(t, o, false);
    else if (id == "T7") t67(t, o, true);
    else if (id == "T8") t8(t, o);
    else if (id == "T9") t9(t, o);
    else if (id == "T11") t11(t, o);
    else if (id == "T12") t12(t, o);
    else if (id == "T13") t13(t, o);
    else if (id == "T14") t14(t, o);
    else if (id == "T15") t15(t, o);
    else if (id == "T18") t18(t, o);
  } catch (const Cancelled&) {
    t.rows.clear();
    t.status = "not run";
    if (t.note.empty()) t.note = "cancelled";
    return t;
  } catch (const std::exception& e) {
    t.status = "fail";
    t.note = std::string("could not run: ") + e.what();
    t.seconds = since(t0);
    return t;
  }
  t.seconds = since(t0);
  // the structure each row ran on, for Open run: the PS melt sample, or the water box for packing (T6/T7 name their cells)
  const bool own = id == "T12" || id == "T13" || id == "T14" || id == "T15" || id == "T18";   // tables that build their own systems
  const std::string sample = id == "T2" ? o.samples + "/water.pdb" : own ? "" : o.samples + "/ps_melt.data";
  for (auto& r : t.rows)
    if (r.file.empty() && !sample.empty() && !o.samples.empty() && std::filesystem::exists(sample)) r.file = sample;
  bool any_fail = false, any_pass = false;
  for (const auto& r : t.rows) any_fail |= r.status == "fail", any_pass |= r.status == "pass";
  t.status = any_fail ? "fail" : any_pass ? "pass" : "info";
  return t;
}

std::string bench_markdown(const std::vector<BenchTable>& tables) {
  std::ostringstream o;
  o << "# CAPS benchmark tables\n\n| Table | Scope | Rows | Status |\n|---|---|---|---|\n";
  for (const auto& t : tables) o << "| " << t.id << " " << t.title << " | " << t.scope << " | " << t.rows.size() << " | " << t.status << " |\n";
  for (const auto& t : tables) {
    o << "\n## " << t.id << " · " << t.title << "\n\n";
    if (!t.rows.empty()) {
      o << "|";
      for (const auto& c : t.columns) o << " " << c << " |";
      o << "\n|";
      for (size_t k = 0; k < t.columns.size(); ++k) o << "---|";
      o << "\n";
      for (const auto& r : t.rows) {
        o << "|";
        for (const auto& c : r.cells) o << " " << c << " |";
        o << "\n";
      }
      o << "\n";
    }
    if (!t.note.empty()) o << t.note << "\n";
  }
  return o.str();
}

std::string bench_csv(const BenchTable& t) {
  auto q = [](const std::string& s) {
    if (s.find_first_of(",\"\n") == std::string::npos) return s;
    std::string r = "\"";
    for (char c : s) r += c == '"' ? std::string("\"\"") : std::string(1, c);
    return r + "\"";
  };
  std::ostringstream o;
  for (size_t k = 0; k < t.columns.size(); ++k) o << (k ? "," : "") << q(t.columns[k]);
  o << "\n";
  for (const auto& r : t.rows) {
    for (size_t k = 0; k < r.cells.size(); ++k) o << (k ? "," : "") << q(r.cells[k]);
    o << "\n";
  }
  return o.str();
}

std::string bench_latex(const std::vector<BenchTable>& tables) {
  auto esc = [](const std::string& s) {
    std::string r;
    for (size_t i = 0; i < s.size(); ++i) {
      const char c = s[i];
      if (c == '%' || c == '&' || c == '_' || c == '#') r += std::string("\\") + c;
      else if (s.compare(i, 2, "±") == 0) { r += "$\\pm$"; ++i; }
      else if (s.compare(i, 2, "×") == 0) { r += "$\\times$"; ++i; }
      else if (s.compare(i, 2, "\xC3\x85") == 0) { r += "\\AA{}"; ++i; }
      else r += c;
    }
    return r;
  };
  std::ostringstream o;
  for (const auto& t : tables) {
    if (t.rows.empty()) continue;
    o << "% " << t.id << " " << t.title << "\n\\begin{table}\n\\centering\n\\caption{" << esc(t.title) << " (" << esc(t.scope) << ")}\n\\begin{tabular}{";
    for (size_t k = 0; k < t.columns.size(); ++k) o << (k ? "r" : "l");
    o << "}\n\\hline\n";
    for (size_t k = 0; k < t.columns.size(); ++k) o << (k ? " & " : "") << esc(t.columns[k]);
    o << " \\\\\n\\hline\n";
    for (const auto& r : t.rows) {
      for (size_t k = 0; k < r.cells.size(); ++k) o << (k ? " & " : "") << esc(r.cells[k]);
      o << " \\\\\n";
    }
    o << "\\hline\n\\end{tabular}\n\\end{table}\n\n";
  }
  return o.str();
}

}  // namespace caps
