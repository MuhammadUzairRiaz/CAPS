// CAPS Equilibrate: protocols as chained Dynamics stages, block-average convergence checks.
#include "caps/equilibrate.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <sstream>

#include "caps/analysis.hpp"
#include "caps/field.hpp"

namespace caps {

namespace {

constexpr double kBarToAtm = 1.0 / 1.01325;

Stage nvt(const std::string& label, double ps, double t) {
  Stage s;
  s.label = label;
  s.ensemble = Ensemble::NVT;
  s.ps = ps;
  s.t_start = t;
  return s;
}

Stage npt(const std::string& label, double ps, double t, double p) {
  Stage s = nvt(label, ps, t);
  s.ensemble = Ensemble::NPT;
  s.pressure = p;
  return s;
}

std::string fmt(const char* f, double a) {
  char b[64];
  std::snprintf(b, sizeof b, f, a);
  return b;
}

}  // namespace

std::vector<Stage> larsen21(const ProtocolParams& p) {
  const double k = p.time_scale, tm = p.t_max, tf = p.t_final, pm = p.p_max;
  // Larsen et al. 2011, Table 1: step, ensemble, conditions, duration (ps).
  std::vector<Stage> s = {
      nvt("1 · heat", 50, tm),         nvt("2 · cool", 50, tf),         npt("3 · compress 0.02 Pmax", 50, tf, 0.02 * pm),
      nvt("4 · heat", 50, tm),         nvt("5 · cool", 100, tf),        npt("6 · compress 0.6 Pmax", 50, tf, 0.6 * pm),
      nvt("7 · heat", 50, tm),         nvt("8 · cool", 100, tf),        npt("9 · compress Pmax", 50, tf, pm),
      nvt("10 · heat", 50, tm),        nvt("11 · cool", 100, tf),       npt("12 · decompress 0.5 Pmax", 5, tf, 0.5 * pm),
      nvt("13 · heat", 5, tm),         nvt("14 · cool", 10, tf),        npt("15 · decompress 0.1 Pmax", 5, tf, 0.1 * pm),
      nvt("16 · heat", 5, tm),         nvt("17 · cool", 10, tf),        npt("18 · decompress 0.01 Pmax", 5, tf, 0.01 * pm),
      nvt("19 · heat", 5, tm),         nvt("20 · cool", 10, tf),        npt("21 · final", 800, tf, p.p_final)};
  for (auto& st : s) st.ps *= k;
  return s;
}

std::vector<Stage> annealing(const ProtocolParams& p) {
  std::vector<Stage> s;
  const double k = p.time_scale;
  s.push_back(npt("settle", p.hold_ps * k, p.t_low, p.p_final));
  for (int c = 1; c <= p.cycles; ++c) {
    const std::string n = std::to_string(c);
    Stage up = npt("cycle " + n + " · heat", p.ramp_ps * k, p.t_low, p.p_final);
    up.t_end = p.t_high;
    s.push_back(up);
    s.push_back(npt("cycle " + n + " · hold", p.hold_ps * k, p.t_high, p.p_final));
    Stage down = npt("cycle " + n + " · cool", p.ramp_ps * k, p.t_high, p.p_final);
    down.t_end = p.t_low;
    s.push_back(down);
    s.push_back(npt("cycle " + n + " · hold", p.hold_ps * k, p.t_low, p.p_final));
  }
  return s;
}

std::vector<Stage> pushoff(const ProtocolParams& p) {
  std::vector<Stage> s;
  const double k = p.time_scale;
  for (double cap : p.caps) {
    Stage st = nvt("push-off · cap " + fmt("%g", cap), p.cap_ps * k, p.t_final);
    st.force_cap = cap;
    s.push_back(st);
  }
  s.push_back(nvt("uncapped", 20 * k, p.t_final));
  s.push_back(npt("NPT", 50 * k, p.t_final, p.p_final));
  return s;
}

std::vector<Stage> protocol_by_name(const std::string& name, const ProtocolParams& p) {
  if (name == "larsen21" || name == "21-step") return larsen21(p);
  if (name == "annealing" || name == "anneal") return annealing(p);
  if (name == "pushoff" || name == "push-off") return pushoff(p);
  throw std::invalid_argument("unknown protocol '" + name + "' (larsen21, annealing, pushoff)");
}

std::vector<Stage> parse_protocol(const std::string& text) {
  std::vector<Stage> out;
  std::istringstream in(text);
  std::string line;
  int lineno = 0;
  while (std::getline(in, line)) {
    ++lineno;
    std::string label;
    if (auto h = line.find('#'); h != std::string::npos) {
      label = line.substr(h + 1);
      line = line.substr(0, h);
      while (!label.empty() && label.front() == ' ') label.erase(0, 1);
      while (!label.empty() && (label.back() == ' ' || label.back() == '\r')) label.pop_back();
    }
    std::istringstream ls(line);
    std::vector<std::string> t;
    for (std::string w; ls >> w;) {
      std::transform(w.begin(), w.end(), w.begin(), [](unsigned char c) { return char(std::tolower(c)); });
      t.push_back(w);
    }
    if (t.empty()) continue;
    auto bad = [&](const std::string& why) { return std::invalid_argument("protocol line " + std::to_string(lineno) + ": " + why); };
    auto num = [&](size_t i) {
      if (i >= t.size()) throw bad("a number is missing after '" + t[i - 1] + "'");
      try { return std::stod(t[i]); } catch (...) { throw bad("'" + t[i] + "' is not a number"); }
    };
    Stage s;
    if (t[0] == "nvt") s.ensemble = Ensemble::NVT;
    else if (t[0] == "npt") s.ensemble = Ensemble::NPT;
    else if (t[0] == "nve") s.ensemble = Ensemble::NVE;
    else throw bad("starts with '" + t[0] + "'; expected nvt, npt or nve");
    bool have_time = false, have_p = false;
    for (size_t i = 1; i < t.size(); ++i) {
      if (i + 1 < t.size() && (t[i + 1] == "ps" || t[i + 1] == "ns" || t[i + 1] == "fs")) {
        const double v = num(i);
        s.ps = t[i + 1] == "ns" ? v * 1000 : t[i + 1] == "fs" ? v / 1000 : v;
        have_time = true;
        ++i;
      } else if (t[i] == "t") {
        s.t_start = num(++i);
        if (i + 1 < t.size() && (t[i + 1] == "to" || t[i + 1] == "->")) {
          i += 2;
          s.t_end = num(i);
        }
      } else if (t[i] == "p") {
        s.pressure = num(++i);
        have_p = true;
        if (i + 1 < t.size() && (t[i + 1] == "bar" || t[i + 1] == "atm")) {
          if (t[i + 1] == "bar") s.pressure *= kBarToAtm;
          ++i;
        }
      } else if (t[i] == "cap") {
        s.force_cap = num(++i);
      } else {
        throw bad("does not understand '" + t[i] + "'");
      }
    }
    if (!have_time || s.ps <= 0) throw bad("needs a duration, e.g. '50 ps'");
    if (s.ensemble == Ensemble::NPT && !have_p) throw bad("npt needs a pressure, e.g. 'P 1 atm'");
    s.label = label.empty() ? std::to_string(out.size() + 1) + " · " + t[0] : label;
    out.push_back(s);
  }
  if (out.empty()) throw std::invalid_argument("the protocol has no stages");
  return out;
}

std::string protocol_text(const std::vector<Stage>& stages) {
  std::string out;
  char b[200];
  for (const auto& s : stages) {
    const char* e = s.ensemble == Ensemble::NPT ? "npt" : s.ensemble == Ensemble::NVE ? "nve" : "nvt";
    std::snprintf(b, sizeof b, "%s %g ps", e, s.ps);
    std::string l = b;
    if (s.ensemble != Ensemble::NVE) {
      std::snprintf(b, sizeof b, " T %g", s.t_start);
      l += b;
      if (s.t_end >= 0) { std::snprintf(b, sizeof b, " to %g", s.t_end); l += b; }
    }
    if (s.ensemble == Ensemble::NPT) {
      std::snprintf(b, sizeof b, "%.2f", s.pressure);   // 0.01 atm, trailing zeros dropped
      std::string pv = b;
      while (pv.back() == '0') pv.pop_back();
      if (pv.back() == '.') pv.pop_back();
      l += " P " + pv + " atm";
    }
    if (s.force_cap > 0) { std::snprintf(b, sizeof b, " cap %g", s.force_cap); l += b; }
    while (l.size() < 30) l += ' ';
    l += ' ';
    out += l + "# " + s.label + "\n";
  }
  return out;
}

double protocol_ps(const std::vector<Stage>& stages) {
  double t = 0;
  for (const auto& s : stages) t += s.ps;
  return t;
}

void equilibrate(System& s, const EquilibrateOptions& o, EquilibrateReport* rep_out) {
  EquilibrateReport rep;
  if (o.stages.empty()) throw std::invalid_argument("the protocol has no stages");
  const double dt = o.md.dt;
  const auto t0 = std::chrono::steady_clock::now();
  int64_t offset = 0;
  const int nstages = int(o.stages.size());
  const int total_stages = nstages + (o.until_converged ? o.max_blocks : 0);
  const double natoms = double(s.atoms.size());

  auto run_stage = [&](const Stage& st, int index, const std::string& label, uint64_t seed, int64_t skip = 0) {
    DynamicsOptions d = o.md;
    const int64_t full = std::max<int64_t>(1, std::llround(st.ps * 1000.0 / dt));
    skip = std::clamp<int64_t>(skip, 0, full - 1);
    d.steps = full - skip;
    // a stage continued after a stop: the ramp from where it was
    d.temperature = st.t_end >= 0 ? st.t_start + (st.t_end - st.t_start) * double(skip) / double(full) : st.t_start;
    d.temperature_end = st.t_end;
    offset += skip;
    if (o.checkpoint)
      d.checkpoint = [&, index, skip](const std::vector<double>& x, const std::vector<double>& v, const Cell& c, int64_t step) {
        o.checkpoint(x, v, c, index - 1, skip + (step - d.step_offset));
      };
    d.thermostat = st.ensemble == Ensemble::NVE ? Thermostat::None
                                                : (o.md.thermostat == Thermostat::None ? Thermostat::Bussi : o.md.thermostat);
    d.barostat = st.ensemble == Ensemble::NPT ? (o.md.barostat == Barostat::None ? Barostat::CRescale : o.md.barostat) : Barostat::None;
    d.pressure = st.pressure;
    d.energy.force_cap = st.force_cap;
    d.seed = seed;
    d.step_offset = offset;
    d.thermo_every = std::max(1, int(std::lround(o.thermo_ps * 1000.0 / dt)));
    d.frame_every = o.frame ? std::max(1, int(std::lround(o.frame_ps * 1000.0 / dt))) : 0;
    d.new_velocities = false;
    const size_t first = rep.thermo.size();
    d.progress = [&](const ThermoRow& r) {
      // the first row of a stage repeats the last row of the previous one
      if (!rep.thermo.empty() && rep.thermo.back().step == r.step) return true;
      rep.thermo.push_back(r);
      return !o.progress || o.progress(index, total_stages, label, r);
    };
    bool first_frame = true;
    d.frame = [&](const std::vector<double>& x, const Cell& c, int64_t step) {
      // the start of a stage is the end of the previous one; record it once
      if (first_frame && offset > 0) { first_frame = false; return; }
      first_frame = false;
      if (o.frame) o.frame(x, c, step);
    };
    DynamicsReport dr;
    run_dynamics(s, d, &dr);
    offset += d.steps;
    StageSummary sum;
    sum.label = label;
    sum.ps = d.steps * dt / 1000.0;
    const size_t a = first + (rep.thermo.size() - first) / 2, n = rep.thermo.size() - a;
    for (size_t k = a; k < rep.thermo.size(); ++k) {
      sum.temperature += rep.thermo[k].temperature / n;
      sum.pressure += rep.thermo[k].pressure / n;
      sum.density += rep.thermo[k].density / n;
      sum.potential += rep.thermo[k].potential / n;
    }
    sum.density_end = s.density();
    rep.stages.push_back(sum);
    return std::make_pair(first, rep.thermo.size());
  };

  const int first_stage = std::clamp(o.start_stage, 0, nstages);
  if (first_stage > 0 || o.start_step > 0) {   // continued after a stop: the stages done before count in the time axis
    for (int k = 0; k < first_stage; ++k) offset += std::max<int64_t>(1, std::llround(o.stages[k].ps * 1000.0 / dt));
    rep.notes.push_back("continued from a checkpoint at stage " + std::to_string(first_stage + 1) + ", step " + std::to_string(o.start_step) +
                        " of it: the stages before are not in this report");
  }
  for (int k = first_stage; k < nstages; ++k) run_stage(o.stages[k], k + 1, o.stages[k].label, o.md.seed + 1000003ull * k, k == first_stage ? o.start_step : 0);

  if (o.until_converged) {
    const Stage& last = o.stages.back();
    Stage prod = npt("production", o.block_ps, last.t_end >= 0 ? last.t_end : last.t_start,
                     last.ensemble == Ensemble::NPT ? last.pressure : 1.0);
    ConvergenceCheck cd{"density", {}, 0, o.tol_density, false};
    ConvergenceCheck ce{"potential energy per atom", {}, 0, o.tol_energy, false};
    ConvergenceCheck cr{"mean Rg", {}, 0, o.tol_rg, false};
    ConvergenceCheck ci{"internal distances", {}, 0, o.tol_internal, false};
    std::vector<std::vector<double>> profiles;   // ⟨R²(n)⟩/(n⟨b²⟩) of each block, by n
    bool chains = true;
    for (int b = 1; b <= o.max_blocks; ++b) {
      const auto [from, to] = run_stage(prod, nstages + b, "production block " + std::to_string(b), o.md.seed + 7919ull * b);
      double dm = 0, em = 0;
      for (size_t k = from; k < to; ++k) { dm += rep.thermo[k].density; em += rep.thermo[k].potential / natoms; }
      cd.blocks.push_back(dm / double(to - from));
      ce.blocks.push_back(em / double(to - from));
      const auto shapes = molecule_shapes(s);
      double rg = 0;
      for (const auto& m : shapes) rg += m.rg;
      cr.blocks.push_back(shapes.empty() ? 0 : rg / shapes.size());
      {
        System whole = s;
        if (!whole.unwrapped) make_molecules_whole(whole);
        const auto id = internal_distances(whole);
        chains = id.chains > 0 && id.n.size() >= 3;
        std::vector<double> prof(id.n.empty() ? 0 : size_t(id.n.back()) + 1, 0.0);
        for (size_t k = 0; k < id.n.size(); ++k) prof[size_t(id.n[k])] = id.ratio[k];
        profiles.push_back(prof);
        double mean = 0;
        for (double r : id.ratio) mean += r;
        ci.blocks.push_back(id.ratio.empty() ? 0 : mean / double(id.ratio.size()));
      }
      rep.blocks = b;
      if (b < 2) continue;
      // the last two changes between consecutive blocks must both be within tolerance
      auto check = [&](ConvergenceCheck& c, bool relative) {
        auto change = [&](size_t i) {
          const double d = c.blocks[i] - c.blocks[i - 1];
          return relative ? std::fabs(d) / std::max(1e-12, std::fabs(c.blocks[i - 1])) : std::fabs(d);
        };
        const size_t m = c.blocks.size();
        c.change = change(m - 1);
        c.ok = m >= 3 && c.change <= c.tolerance && change(m - 2) <= c.tolerance;
      };
      check(cd, true);
      check(ce, false);
      check(cr, true);
      // the internal distances: the largest relative change at any n between the last two blocks, and from the target
      if (chains && profiles.size() >= 2) {
        auto dev = [&](const std::vector<double>& a, const std::vector<double>& ref) {
          double m = 0;
          for (size_t k = 2; k < std::min(a.size(), ref.size()); ++k)
            if (ref[k] > 0 && a[k] > 0) m = std::max(m, std::fabs(a[k] / ref[k] - 1));
          return m;
        };
        const auto& now = profiles.back();
        ci.change = dev(now, profiles[profiles.size() - 2]);
        if (!o.internal_target.empty()) ci.change = std::max(ci.change, dev(now, o.internal_target));
        ci.ok = profiles.size() >= 3 && ci.change <= ci.tolerance;
      } else {
        ci.ok = true;   // no chains: nothing to check
      }
      if (b >= o.min_blocks && cd.ok && ce.ok && cr.ok && ci.ok) { rep.converged = true; break; }
    }
    rep.checks = {cd, ce, cr};
    if (chains) rep.checks.push_back(ci);
    char b[200];
    std::snprintf(b, sizeof b, "%s after %d production blocks of %g ps (density ±%.1f%%, energy ±%.3g kcal/mol per atom, Rg ±%.0f%%, internal distances ±%.0f%%%s)",
                  rep.converged ? "converged" : "not converged", rep.blocks, o.block_ps, 100 * o.tol_density, o.tol_energy, 100 * o.tol_rg,
                  100 * o.tol_internal, o.internal_target.empty() ? "" : " and of the target curve");
    rep.notes.push_back(b);
    rep.notes.push_back("the checks catch drift in density, energy and chain size over the run; long chains relax far more slowly than any "
                        "MD run, so a passed Rg check is necessary, not sufficient");
  }

  rep.ps = offset * dt / 1000.0;
  rep.seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
  char b[200];
  std::snprintf(b, sizeof b, "%zu stages · %.1f ps in %.1f s · final density %.4f g/cm³", rep.stages.size(), rep.ps, rep.seconds, s.density());
  rep.notes.insert(rep.notes.begin(), b);
  if (rep_out) *rep_out = std::move(rep);
}

}  // namespace caps
