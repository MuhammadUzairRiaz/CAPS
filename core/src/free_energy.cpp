#include "caps/free_energy.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

#include "caps/dynamics.hpp"
#include "caps/uff.hpp"

namespace caps {

SolvationResult solvation_free_energy(const System& s0, const SolvationOptions& o) {
  const size_t n = s0.atoms.size();
  if (o.solute.size() != n || std::none_of(o.solute.begin(), o.solute.end(), [](char c) { return c != 0; }))
    throw std::invalid_argument("mark the solute's atoms");
  if (std::all_of(o.solute.begin(), o.solute.end(), [](char c) { return c != 0; })) throw std::invalid_argument("the solute is the whole structure: no solvent");
  if (!s0.cell.valid()) throw std::invalid_argument("a solvation free energy needs a periodic box of solvent");
  const ForceField ff = o.field ? *o.field : default_forcefield(s0);
  bool charged = false;
  for (size_t i = 0; i < n; ++i) charged = charged || (o.solute[i] && std::fabs(ff.charge[i]) > 1e-12);
  SolvationResult R;
  System s = s0;
  DynamicsOptions d;
  d.field = std::make_shared<const ForceField>(ff);
  d.energy = o.energy;
  d.energy.tail = false;
  if (d.energy.electrostatics == EnergyOptions::Electrostatics::PME) d.energy.electrostatics = EnergyOptions::Electrostatics::DSF;
  d.energy.alchemical = &o.solute;
  d.dt = o.dt;
  d.temperature = o.temperature;
  d.thermostat = Thermostat::Bussi;
  d.tau_t = o.tau_t;
  d.seed = o.seed;
  d.frame_every = 0;
  d.thermo_every = 1000;
  d.new_velocities = s.velocities.size() != n;
  const int64_t neq = std::llround(o.equilibrate_ps * 1000 / o.dt), nrun = std::max<int64_t>(10, std::llround(o.ps * 1000 / o.dt));
  // the windows in order: electrostatics off first (when the solute has charges), then the Lennard-Jones
  struct W { std::string leg; double lc, ll; };
  std::vector<W> plan;
  if (charged) for (double l : o.coul_windows) plan.push_back({"coulomb", l, 1.0});
  for (double l : o.lj_windows) plan.push_back({"lj", 0.0, l});
  if (!charged) for (auto& w : plan) w.lc = 0.0;
  for (size_t k = 0; k < plan.size(); ++k) {
    d.energy.lambda_coul = plan[k].lc;
    d.energy.lambda_lj = plan[k].ll;
    const std::string stage = plan[k].leg + " λ = " + std::to_string(plan[k].leg == "coulomb" ? plan[k].lc : plan[k].ll).substr(0, 5);
    if (neq > 0) {
      d.steps = neq;
      d.each_step = nullptr;
      d.progress = [&](const ThermoRow&) { return !o.progress || o.progress(stage + " (settling)", double(k) / plan.size()); };
      run_dynamics(s, d);
      d.new_velocities = false;
    }
    std::vector<double> samples;
    d.steps = nrun;
    d.each_step = [&](const EnergyTerms& e, const std::vector<double>&, const Cell&, int64_t) {
      samples.push_back(plan[k].leg == "coulomb" ? e.dudl_coul : e.dudl_lj);
    };
    d.progress = [&](const ThermoRow& r) { return !o.progress || o.progress(stage, (double(k) + double(r.step) / double(nrun)) / plan.size()); };
    run_dynamics(s, d);
    d.new_velocities = false;
    // block mean and error
    const int nb = std::clamp(o.blocks, 2, int(samples.size()));
    std::vector<double> bm;
    for (int b = 0; b < nb; ++b) {
      const size_t a = samples.size() * size_t(b) / nb, e = samples.size() * size_t(b + 1) / nb;
      bm.push_back(std::accumulate(samples.begin() + long(a), samples.begin() + long(e), 0.0) / double(e - a));
    }
    const double m = std::accumulate(bm.begin(), bm.end(), 0.0) / nb;
    double v = 0;
    for (double x : bm) v += (x - m) * (x - m);
    R.windows.push_back({plan[k].leg, plan[k].leg == "coulomb" ? plan[k].lc : plan[k].ll, m, std::sqrt(v / (nb - 1) / nb)});
  }
  // trapezoid rule over each leg, λ ascending; errors added in quadrature with the trapezoid weights
  auto integrate = [&](const std::string& leg, double& err2) {
    std::vector<TiWindow> w;
    for (const auto& x : R.windows) if (x.leg == leg) w.push_back(x);
    std::sort(w.begin(), w.end(), [](const TiWindow& a, const TiWindow& b) { return a.lambda < b.lambda; });
    double sum = 0;
    for (size_t k = 0; k + 1 < w.size(); ++k) {
      const double h = w[k + 1].lambda - w[k].lambda;
      sum += 0.5 * h * (w[k].dudl + w[k + 1].dudl);
      err2 += 0.25 * h * h * (w[k].error * w[k].error + w[k + 1].error * w[k + 1].error);
    }
    return sum;
  };
  double e2 = 0;
  R.dg_coul = charged ? integrate("coulomb", e2) : 0.0;
  R.dg_lj = integrate("lj", e2);
  R.dg = R.dg_coul + R.dg_lj;
  R.dg_err = std::sqrt(e2);
  char b[200];
  std::snprintf(b, sizeof b, "TI over %zu windows of %.1f ps at %.0f K: ΔG_solv = %.3f ± %.3f kcal/mol (electrostatics %.3f, Lennard-Jones %.3f)", plan.size(), o.ps,
                o.temperature, R.dg, R.dg_err, R.dg_coul, R.dg_lj);
  R.notes.push_back(b);
  if (!charged) R.notes.push_back("the solute has no charges: the Lennard-Jones leg alone");
  R.notes.push_back("pairwise electrostatics (DSF) and no tail correction, as the soft-core pairs need; compare like with like");
  return R;
}

}  // namespace caps
