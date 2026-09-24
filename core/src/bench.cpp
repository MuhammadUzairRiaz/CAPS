// CAPS Bench: the built-in validation suite (see bench.hpp).
#include "caps/bench.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <map>
#include <memory>
#include <mutex>
#include <numeric>
#include <sstream>
#include <thread>

#include "caps/config.hpp"
#include "caps/dynamics.hpp"
#include "caps/field.hpp"
#include "caps/io.hpp"
#include "caps/molecule.hpp"
#include "caps/pack.hpp"
#include "caps/relax.hpp"
#include "caps/render.hpp"

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
      {"T6", "Properties vs experiment", "PE PP PS PMMA PET PC PA6 PDMS", {}, "needs equilibrated cells (hours per polymer); Analyze compares each cell with the experimental ranges"},
      {"T7", "Chain statistics", "C∞ vs literature", {}, "needs long equilibrated melts; Analyze reports C∞ from the internal distances"},
      {"T8", "Parallel scaling", "force evaluation, threads", {"System", "Threads", "ms / evaluation", "Speed-up", "Efficiency"}, nullptr},
      {"T9", "Reproducibility", "the same run twice", {"Run", "Steps", "Max |Δx| (Å)", "Status"}, nullptr},
      {"T10", "Builder feature coverage", "against vendor documentation", {}, "a written comparison, not a measurement"},
      {"T11", "Rendering", "CPU renderer, 1920 × 1080", {"System", "Atoms", "Supersampling", "ms / frame"}, nullptr},
      {"T12", "Molecule builder", "stereochemistry, round trips, embedding", {"Check", "Cases", "Correct", "Status"}, nullptr},
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
  t.note = "Built-in GAFF, 10 Å cut-off, DSF Coulomb, dt 1 fs, Bussi thermostat; mean ± sd over the repeats.";
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
  if (!o.forcefields.empty()) bo.forcefield = o.forcefields + "/gaff-amber25-dlfield.json";
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
    else if (id == "T8") t8(t, o);
    else if (id == "T9") t9(t, o);
    else if (id == "T11") t11(t, o);
    else if (id == "T12") t12(t, o);
  } catch (const Cancelled&) {
    t.rows.clear();
    t.status = "not run";
    t.note = "cancelled";
    return t;
  } catch (const std::exception& e) {
    t.status = "fail";
    t.note = std::string("could not run: ") + e.what();
    t.seconds = since(t0);
    return t;
  }
  t.seconds = since(t0);
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
