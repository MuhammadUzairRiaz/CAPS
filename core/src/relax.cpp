// CAPS Relax: steepest descent, Polak–Ribière conjugate gradient, L-BFGS and FIRE minimisers; capped-force push-off;
// affine compression to a target density; isotropic box relaxation; LAMMPS data export with the force field.
//
// L-BFGS: Liu and Nocedal, Math. Program. 45, 503 (1989), two-loop recursion, m = 10.
// FIRE: Bitzek, Koskinen, Gähler, Moseler and Gumbsch, Phys. Rev. Lett. 97, 170201 (2006), standard parameters.
// Push-off: Auhl, Everaers, Grest, Kremer and Plimpton, J. Chem. Phys. 119, 12718 (2003), here as capped-force
// minimisation stages rather than capped-force dynamics.
#include "caps/uff.hpp"
#include "caps/relax.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <deque>
#include <fstream>
#include <map>
#include <array>
#include <tuple>

#include "caps/elements.hpp"

namespace caps {

Minimiser minimiser_from_string(const std::string& s) {
  if (s == "sd" || s == "steepest") return Minimiser::SteepestDescent;
  if (s == "cg") return Minimiser::ConjugateGradient;
  if (s == "lbfgs" || s == "l-bfgs") return Minimiser::LBFGS;
  if (s == "fire") return Minimiser::FIRE;
  throw std::invalid_argument("unknown minimiser '" + s + "' (sd, cg, lbfgs, fire)");
}

const char* to_string(Minimiser m) {
  switch (m) {
    case Minimiser::SteepestDescent: return "steepest descent";
    case Minimiser::ConjugateGradient: return "Polak–Ribière CG";
    case Minimiser::LBFGS: return "L-BFGS";
    case Minimiser::FIRE: return "FIRE";
  }
  return "?";
}

double max_force(const std::vector<double>& f) {
  double m = 0;
  for (size_t i = 0; i + 2 < f.size(); i += 3) m = std::max(m, f[i] * f[i] + f[i + 1] * f[i + 1] + f[i + 2] * f[i + 2]);
  return std::sqrt(m);
}

namespace {

double dotv(const std::vector<double>& a, const std::vector<double>& b) {
  double s = 0;
  for (size_t i = 0; i < a.size(); ++i) s += a[i] * b[i];
  return s;
}

double max_disp(const std::vector<double>& d) { return max_force(d); }

}  // namespace

RelaxStage minimise(Evaluator& ev, std::vector<double>& x, const Cell& cell, const RelaxOptions& o, const std::string& name,
                    int* evaluations) {
  RelaxStage st;
  st.name = name;
  const size_t n3 = x.size();
  std::vector<double> f, fn, xn(n3), d(n3);
  int nev = 0;
  auto eval = [&](const std::vector<double>& p, std::vector<double>& out) {
    ++nev;
    double en = ev.compute(p, cell, out).total();
    for (const auto& rs : o.restraints) {   // k (r − r0)²
      if (3 * size_t(std::max(rs.i, rs.j)) + 2 >= p.size()) continue;
      Vec3 dv{p[3 * rs.j] - p[3 * rs.i], p[3 * rs.j + 1] - p[3 * rs.i + 1], p[3 * rs.j + 2] - p[3 * rs.i + 2]};
      if (cell.valid()) dv = cell.minimum_image(dv);
      const double r = norm(dv);
      if (r < 1e-9) continue;
      en += rs.k * (r - rs.r0) * (r - rs.r0);
      const double g = 2 * rs.k * (r - rs.r0) / r;   // −dE/dx_j = −g d
      for (int c = 0; c < 3; ++c) out[3 * rs.j + c] -= g * dv[c], out[3 * rs.i + c] += g * dv[c];
    }
    for (size_t i = 0; i < o.fixed.size() && 3 * i + 2 < out.size(); ++i)
      if (o.fixed[i]) out[3 * i] = out[3 * i + 1] = out[3 * i + 2] = 0;   // held atoms feel no force and do not move
    return en;
  };

  double e = eval(x, f);
  double fmax = max_force(f);
  RelaxProgress pr;
  pr.stage = name;
  auto report = [&](int it) {
    if (!o.progress) return;
    pr.iteration = it;
    pr.energy = e;
    pr.fmax = fmax;
    if (!o.progress(pr)) throw RelaxCancelled();
  };

  int it = 0;
  st.stopped_by = "iterations";
  if (fmax < o.ftol) st.stopped_by = "force";

  if (o.method == Minimiser::FIRE) {
    const int nmin = 5;
    const double finc = 1.1, fdec = 0.5, astart = 0.1, falpha = 0.99, dtmax = 0.5, maxstep = 0.1;
    double dt = 0.05, alpha = astart;
    int npos = 0;
    std::vector<double> v(n3, 0.0);
    for (; it < o.max_iterations && st.stopped_by != "force"; ++it) {
      const double p = dotv(f, v);
      if (p > 0) {
        const double vn = std::sqrt(dotv(v, v)), fnn = std::sqrt(dotv(f, f));
        for (size_t k = 0; k < n3; ++k) v[k] = (1 - alpha) * v[k] + alpha * vn * f[k] / std::max(fnn, 1e-300);
        if (++npos > nmin) { dt = std::min(dt * finc, dtmax); alpha *= falpha; }
      } else {
        std::fill(v.begin(), v.end(), 0.0);
        dt *= fdec;
        alpha = astart;
        npos = 0;
      }
      for (size_t k = 0; k < n3; ++k) { v[k] += dt * f[k]; d[k] = dt * v[k]; }
      const double md = max_disp(d);
      if (md > maxstep) for (auto& q : d) q *= maxstep / md;
      for (size_t k = 0; k < n3; ++k) x[k] += d[k];
      e = eval(x, f);
      fmax = max_force(f);
      if (fmax < o.ftol) st.stopped_by = "force";
      if ((it & 15) == 0) report(it);
    }
  } else {
    // Line-search methods. Backtracking (Armijo) with the first trial step limited to a largest atomic move.
    const int m = 10;
    std::deque<std::vector<double>> S, Y;
    std::deque<double> rho;
    std::vector<double> fprev, xprev, dprev(n3, 0.0);
    double step = 0.1;        // Å, largest atomic displacement of the first trial
    const double step_max = o.method == Minimiser::SteepestDescent ? 0.3 : 0.5;
    for (; it < o.max_iterations && st.stopped_by != "force"; ++it) {
      // direction
      if (o.method == Minimiser::SteepestDescent) {
        d = f;
      } else if (o.method == Minimiser::ConjugateGradient) {
        if (fprev.empty()) d = f;
        else {
          double num = 0, den = dotv(fprev, fprev);
          for (size_t k = 0; k < n3; ++k) num += f[k] * (f[k] - fprev[k]);
          const double beta = den > 0 ? std::max(0.0, num / den) : 0.0;
          for (size_t k = 0; k < n3; ++k) d[k] = f[k] + beta * dprev[k];
        }
      } else {
        // L-BFGS two-loop on the gradient g = −f
        std::vector<double> q(n3);
        for (size_t k = 0; k < n3; ++k) q[k] = -f[k];
        std::vector<double> al(S.size());
        for (int j = static_cast<int>(S.size()) - 1; j >= 0; --j) {
          al[j] = rho[j] * dotv(S[j], q);
          for (size_t k = 0; k < n3; ++k) q[k] -= al[j] * Y[j][k];
        }
        double gamma = 1.0 / 300.0;   // initial inverse-Hessian scale (Å² mol/kcal), about a stiff bond
        if (!S.empty()) gamma = dotv(S.back(), Y.back()) / dotv(Y.back(), Y.back());
        for (auto& v : q) v *= gamma;
        for (size_t j = 0; j < S.size(); ++j) {
          const double b = rho[j] * dotv(Y[j], q);
          for (size_t k = 0; k < n3; ++k) q[k] += S[j][k] * (al[j] - b);
        }
        for (size_t k = 0; k < n3; ++k) d[k] = -q[k];
      }
      double slope = dotv(f, d);   // −g·d, positive for a descent direction
      if (slope <= 0) {
        d = f;
        slope = dotv(f, f);
        S.clear(); Y.clear(); rho.clear();
      }
      // trial step
      const double md = max_disp(d);
      double a;
      if (o.method == Minimiser::LBFGS && !S.empty()) a = std::min(1.0, step_max / md);
      else a = step / md;
      bool ok = false;
      double en = 0;
      for (int ls = 0; ls < 30; ++ls) {
        for (size_t k = 0; k < n3; ++k) xn[k] = x[k] + a * d[k];
        en = eval(xn, fn);
        if (en <= e - 1e-4 * a * slope) { ok = true; break; }
        a *= 0.5;
      }
      if (!ok) {
        if (!S.empty() || !fprev.empty()) {   // start over from steepest descent once
          S.clear(); Y.clear(); rho.clear(); fprev.clear();
          step = 0.02;
          continue;
        }
        st.stopped_by = "line search";
        break;
      }
      const double moved = a * md;
      step = std::clamp(moved * 1.5, 1e-4, step_max);
      if (o.method == Minimiser::LBFGS) {
        std::vector<double> s(n3), y(n3);
        for (size_t k = 0; k < n3; ++k) { s[k] = xn[k] - x[k]; y[k] = f[k] - fn[k]; }
        const double sy = dotv(s, y);
        if (sy > 1e-10) {
          S.push_back(std::move(s)); Y.push_back(std::move(y)); rho.push_back(1.0 / sy);
          if (static_cast<int>(S.size()) > m) { S.pop_front(); Y.pop_front(); rho.pop_front(); }
        }
      }
      fprev = f;
      dprev = d;
      const double de = e - en;
      x.swap(xn);
      f.swap(fn);
      e = en;
      fmax = max_force(f);
      if (fmax < o.ftol) st.stopped_by = "force";
      else if (de >= 0 && de < o.etol * std::max(1.0, std::fabs(e))) { st.stopped_by = "energy"; ++it; break; }
      if ((it & 7) == 0) report(it);
    }
  }
  report(it);
  st.iterations = it;
  st.energy = e;
  st.fmax = fmax;
  if (evaluations) *evaluations += nev;
  return st;
}

namespace {

std::vector<double> flat(const System& s) {
  std::vector<double> x(3 * s.atoms.size());
  for (size_t i = 0; i < s.atoms.size(); ++i)
    for (int k = 0; k < 3; ++k) x[3 * i + k] = s.atoms[i].pos[k];
  return x;
}

// Scale the cell and all positions about the cell origin (affine, keeps molecules whole).
void scale_affine(std::vector<double>& x, Cell& c, double s) {
  for (size_t i = 0; i < x.size(); i += 3)
    for (int k = 0; k < 3; ++k) x[i + k] = c.origin[k] + s * (x[i + k] - c.origin[k]);
  c.a = c.a * s;
  c.b = c.b * s;
  c.c = c.c * s;
}

}  // namespace

void relax(System& s, const RelaxOptions& o, RelaxReport* rep_out) {
  RelaxReport rep;
  if (s.atoms.empty()) throw FieldError("nothing to relax: no atoms");
  if (o.field && o.field->atom_type.size() != s.atoms.size())
    throw FieldError("the assigned force field is for " + std::to_string(o.field->atom_type.size()) + " atoms, the structure has " +
                     std::to_string(s.atoms.size()));
  const ForceField ff = o.field ? *o.field : default_forcefield(s);
  rep.field = ff.name;
  for (const auto& n : ff.notes) rep.notes.push_back(n);

  const double mass = s.total_mass();
  auto density_of = [&](const Cell& c) {
    constexpr double kNA = 6.02214076e23;
    return c.valid() ? mass / kNA / (c.volume() * 1e-24) : 0.0;
  };

  Cell cell = s.cell;
  std::vector<double> x = flat(s), f;
  if (cell.valid()) {
    const double v = cell.volume();
    const double w = std::min({v / norm(cross(cell.b, cell.c)), v / norm(cross(cell.c, cell.a)), v / norm(cross(cell.a, cell.b))});
    if (o.target_density > 0) {
      const double wt = w * std::cbrt(density_of(cell) / o.target_density);
      if (wt < 2 * o.energy.cutoff)
        rep.notes.push_back("the cell will be narrower than twice the cut-off; periodic images are included in the pair list");
    }
  } else if (o.target_density > 0 || o.relax_box) {
    throw FieldError("compression and box relaxation need a periodic cell");
  }

  Evaluator ev(ff, o.energy);
  rep.initial = ev.compute(x, cell, f);
  rep.fmax_initial = max_force(f);
  rep.density_initial = density_of(cell);

  // Plan the stages so progress can show "k of n".
  int stages = 1;
  if (o.pushoff) stages += static_cast<int>(o.pushoff_caps.size());
  std::vector<double> dens;
  if (o.target_density > 0 && rep.density_initial > 0) {
    double d = rep.density_initial;
    const double step = std::max(0.005, o.compress_step);
    if (o.target_density > d) {
      while (d < o.target_density * (1 - 1e-9)) { d = std::min(o.target_density, d * (1 + step)); dens.push_back(d); }
    } else {
      while (d > o.target_density * (1 + 1e-9)) { d = std::max(o.target_density, d / (1 + step)); dens.push_back(d); }
    }
  }
  stages += static_cast<int>(dens.size());
  if (o.relax_box) stages += 1;

  int stage_index = 0;
  auto run = [&](const std::string& name, const EnergyOptions& eo, RelaxOptions mo) {
    ev.set_options(eo);
    ++stage_index;
    auto user = o.progress;
    const double rho = density_of(cell);
    mo.progress = [&, user, rho, name](const RelaxProgress& p) {
      if (!user) return true;
      RelaxProgress q = p;
      q.stage = name;
      q.stage_index = stage_index;
      q.stages = stages;
      q.density = rho;
      return user(q);
    };
    RelaxStage st = minimise(ev, x, cell, mo, name, &rep.evaluations);
    st.density = rho;
    const EnergyTerms t = ev.compute(x, cell, f);
    ++rep.evaluations;
    st.pressure = cell.valid() ? pressure_atm(t.virial, cell.volume()) : 0.0;
    rep.iterations += st.iterations;
    rep.stages.push_back(st);
    if (o.snapshot) o.snapshot(x, cell, name);
    return st;
  };

  auto soft = [&](double cap) {
    EnergyOptions eo = o.energy;
    eo.force_cap = cap;
    eo.coulomb = false;
    return eo;
  };
  RelaxOptions quick = o;
  quick.method = Minimiser::FIRE;
  quick.ftol = std::max(o.ftol, 5.0);
  quick.max_iterations = std::min(o.max_iterations, 400);

  if (o.pushoff)
    for (double cap : o.pushoff_caps) {
      char nm[64];
      std::snprintf(nm, sizeof nm, "push-off, force cap %g", cap);
      run(nm, soft(cap), quick);
    }

  // Compression: scale, push off the new close contacts gently, then minimise loosely.
  RelaxOptions loose = o;
  loose.ftol = std::max(o.ftol, 2.0);
  loose.max_iterations = std::min(o.max_iterations, 1500);
  for (double d : dens) {
    scale_affine(x, cell, std::cbrt(density_of(cell) / d));
    char nm[64];
    std::snprintf(nm, sizeof nm, "compress to %.3f g/cm³", d);
    ev.set_options(soft(100));
    minimise(ev, x, cell, quick, nm, &rep.evaluations);
    run(nm, o.energy, loose);
  }

  RelaxStage last = run("minimise", o.energy, o);

  if (o.relax_box) {
    // Secant iteration on ln V towards the target pressure, minimising positions at each volume.
    double bulk = 3.0e4;   // atm, first guess (about 3 GPa)
    double lnv_prev = 0, p_prev = 0;
    bool have_prev = false;
    double max_step = 0.03;   // in ln V per cycle
    int cycles = 0;
    RelaxStage st = last;
    for (; cycles < o.box_cycles; ++cycles) {
      const double p = st.pressure;
      if (std::fabs(p - o.pressure) < o.pressure_tol) break;
      const double lnv = std::log(cell.volume());
      if (have_prev && std::fabs(lnv - lnv_prev) > 1e-9) {
        const double b = -(p - p_prev) / (lnv - lnv_prev);
        if (b > 1e3 && b < 5e6) bulk = b;   // up to 500 GPa: ceramics (MgO 160 GPa) and diamond (443 GPa) too
        if ((p - o.pressure) * (p_prev - o.pressure) < 0) max_step *= 0.5;   // crossed the target: smaller steps
      }
      lnv_prev = lnv;
      p_prev = p;
      have_prev = true;
      const double dlnv = std::clamp((p - o.pressure) / bulk, -max_step, max_step);
      scale_affine(x, cell, std::exp(dlnv / 3));
      ev.set_options(o.energy);
      RelaxOptions mo = o;
      auto user = o.progress;
      const double rho = density_of(cell);
      mo.progress = [&, user, rho](const RelaxProgress& q0) {
        if (!user) return true;
        RelaxProgress q = q0;
        q.stage = "box";
        q.stage_index = stages;
        q.stages = stages;
        q.density = rho;
        q.pressure = p_prev;
        return user(q);
      };
      st = minimise(ev, x, cell, mo, "box", &rep.evaluations);
      const EnergyTerms t = ev.compute(x, cell, f);
      ++rep.evaluations;
      st.pressure = pressure_atm(t.virial, cell.volume());
      st.density = rho;
      rep.iterations += st.iterations;
    }
    char nm[96];
    std::snprintf(nm, sizeof nm, "box relaxed to %.0f atm in %d cycles", st.pressure, cycles);
    st.name = nm;
    rep.stages.push_back(st);
    if (o.snapshot) o.snapshot(x, cell, "box");
    if (std::fabs(st.pressure - o.pressure) >= o.pressure_tol)
      rep.notes.push_back("box pressure did not reach the tolerance in " + std::to_string(o.box_cycles) + " cycles");
    rep.notes.push_back("box relaxation is 0 K mechanical equilibrium; the density at a temperature needs Dynamics (NPT)");
    last = st;
  }

  rep.final = ev.compute(x, cell, f);
  rep.fmax_final = max_force(f);
  rep.density_final = density_of(cell);
  rep.pressure_final = cell.valid() ? pressure_atm(rep.final.virial, cell.volume()) : 0.0;
  rep.list_builds = ev.list_builds();
  rep.converged = last.stopped_by == "force" || last.stopped_by == "energy";
  if (!rep.converged) rep.notes.push_back("the final minimisation stopped by " + last.stopped_by + " before the force tolerance");
  for (size_t k = 0; k < o.restraints.size() && k < 8; ++k) {   // where each restrained pair ended up
    const auto& rs = o.restraints[k];
    if (3 * size_t(std::max(rs.i, rs.j)) + 2 >= x.size()) continue;
    Vec3 dv{x[3 * rs.j] - x[3 * rs.i], x[3 * rs.j + 1] - x[3 * rs.i + 1], x[3 * rs.j + 2] - x[3 * rs.i + 2]};
    if (cell.valid()) dv = cell.minimum_image(dv);
    char b[160];
    std::snprintf(b, sizeof b, "restraint %u–%u: target %.3f Å, k %.1f kcal/mol/Å², final %.3f Å", rs.i + 1, rs.j + 1, rs.r0, rs.k, norm(dv));
    rep.notes.push_back(b);
  }
  if (o.restraints.size() > 8) rep.notes.push_back(std::to_string(o.restraints.size() - 8) + " more restraints");

  for (size_t i = 0; i < s.atoms.size(); ++i) {
    for (int k = 0; k < 3; ++k) s.atoms[i].pos[k] = x[3 * i + k];
    s.atoms[i].charge = ff.charge[i];
  }
  s.has_charges = true;
  s.cell = cell;
  s.unwrapped = true;

  char buf[256];
  std::snprintf(buf, sizeof buf, "energy %.1f → %.1f kcal/mol, largest force %.3g → %.3g kcal/mol/Å, %d iterations (%s)",
                rep.initial.total(), rep.final.total(), rep.fmax_initial, rep.fmax_final, rep.iterations, to_string(o.method));
  rep.notes.insert(rep.notes.begin(), buf);
  if (cell.valid()) {
    std::snprintf(buf, sizeof buf, "density %.4f → %.4f g/cm³, pressure %.0f atm (0 K virial)", rep.density_initial, rep.density_final,
                  rep.pressure_final);
    rep.notes.insert(rep.notes.begin() + 1, buf);
  }
  if (rep_out) *rep_out = std::move(rep);
}

// LAMMPS data with the force field: lammps_data.cpp

}  // namespace caps
