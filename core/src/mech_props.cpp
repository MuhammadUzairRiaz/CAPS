// Mechanics results as Analyze properties.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>

#include "caps/mechanics.hpp"

namespace caps {

namespace {

constexpr double NaN = std::numeric_limits<double>::quiet_NaN();

Property prop(const std::string& id, const std::string& name, const std::string& unit, const std::string& method, double v, double e = NaN) {
  Property p;
  p.id = id;
  p.name = name;
  p.unit = unit;
  p.method = method;
  p.value = v;
  p.error = e;
  return p;
}

}  // namespace

std::vector<Property> elastic_properties(const ElasticResult& r, const std::string& suffix) {
  const bool fl = suffix == "_fluct";
  const std::string how = fl ? " (fluctuations)" : " (static strain)";
  std::vector<Property> out;
  const double diag = (r.C[0][0] + r.C[1][1] + r.C[2][2]) / 3;
  double derr = NaN;
  if (!std::isnan(r.err[0][0])) derr = std::sqrt(r.err[0][0] * r.err[0][0] + r.err[1][1] * r.err[1][1] + r.err[2][2] * r.err[2][2]) / 3;
  Property c = prop("cij" + suffix, "Elastic constants" + how, "GPa", r.method, diag, derr);
  for (int I = 0; I < 6; ++I)
    for (int J = I; J < 6; ++J) {
      char k[48];
      std::snprintf(k, sizeof k, "C%d%d (GPa)", I + 1, J + 1);
      c.extra[k] = r.C[I][J];
      if (!std::isnan(r.err[I][J]) && (I == J || (I < 3 && J < 3))) {
        std::snprintf(k, sizeof k, "C%d%d error (GPa)", I + 1, J + 1);
        c.extra[k] = r.err[I][J];
      }
    }
  const char* sn[6] = {"xx", "yy", "zz", "yz", "xz", "xy"};
  for (int v = 0; v < 6; ++v) c.extra[std::string("prestress σ") + sn[v] + " (GPa)"] = r.prestress[v];
  if (fl)
    for (int I = 0; I < 3; ++I) {
      c.extra["Born C" + std::to_string(I + 1) + std::to_string(I + 1) + " (GPa)"] = r.born[I][I];
      c.extra["fluctuation C" + std::to_string(I + 1) + std::to_string(I + 1) + " (GPa)"] = r.fluct[I][I];
    }
  c.extra["configurations"] = r.configurations;
  c.notes = r.notes;
  out.push_back(c);
  const std::string avg = "Hill average (mean of Voigt and Reuss) of the elastic constants" + how;
  Property E = prop("youngs" + suffix, "Young's modulus" + how, "GPa", avg, r.E_hill);
  E.extra["Voigt bound K (GPa)"] = r.K_voigt;
  E.extra["Reuss bound K (GPa)"] = r.K_reuss;
  out.push_back(E);
  out.push_back(prop("bulk" + suffix, "Bulk modulus" + how, "GPa", avg, r.K_hill));
  Property G = prop("shear" + suffix, "Shear modulus" + how, "GPa", avg, r.G_hill);
  G.extra["Voigt bound (GPa)"] = r.G_voigt;
  G.extra["Reuss bound (GPa)"] = r.G_reuss;
  out.push_back(G);
  Property nu = prop("poisson" + suffix, "Poisson's ratio" + how, "", avg, r.nu_hill);
  nu.extra["Lamé λ (GPa)"] = r.lambda_hill;
  out.push_back(nu);
  return out;
}

std::vector<Property> tensile_properties(const TensileResult& r) {
  std::vector<Property> out;
  Property m = prop("tensile_modulus", "Tensile modulus", "GPa", r.method, r.modulus, r.modulus_err);
  if (!std::isnan(r.poisson)) {
    m.extra["Poisson's ratio"] = r.poisson;
    if (!std::isnan(r.poisson_err)) m.extra["Poisson's ratio error"] = r.poisson_err;
  }
  m.extra["peak stress (MPa)"] = r.peak_stress;
  m.extra["strain at peak"] = r.peak_strain;
  m.notes = r.notes;
  Series raw{"stress–strain", "strain", "stress (MPa)", {}, {}}, sm{"stress–strain, smoothed ±0.5 %", "strain", "stress (MPa)", {}, {}};
  Series lat{"lateral strain", "axial strain", "mean lateral strain", {}, {}};
  for (size_t i = 0; i < r.curve.size(); ++i) {
    raw.x.push_back(r.curve[i].strain);
    raw.y.push_back(r.curve[i].stress);
    sm.x.push_back(r.curve[i].strain);
    sm.y.push_back(i < r.smooth.size() ? r.smooth[i] : r.curve[i].stress);
    lat.x.push_back(r.curve[i].strain);
    lat.y.push_back(0.5 * (r.curve[i].lateral1 + r.curve[i].lateral2));
  }
  m.series = {sm, raw, lat};
  out.push_back(m);
  Property y = prop("yield", "Yield stress (0.2 % offset)", "MPa", "first point of the smoothed curve below E (ε − 0.002)",
                    r.yield_strain > 0 ? r.yield_stress : NaN);
  if (r.yield_strain > 0) y.extra["yield strain"] = r.yield_strain;
  else y.notes.push_back("not reached within the deformation");
  out.push_back(y);
  return out;
}

std::vector<Property> pull_properties(const PullResult& r, bool normal) {
  Property p = prop(normal ? "pull_normal" : "pull_shear", normal ? "Peak normal stress (debonding)" : "Interfacial shear strength", "MPa", r.method, r.strength);
  p.extra["peak force (kcal/mol/Å)"] = r.peak_force;
  p.extra["displacement at peak (Å)"] = r.peak_displacement;
  p.extra["work (mJ/m²)"] = r.work;
  p.extra["interfaces"] = r.interfaces;
  p.extra["area per interface (Å²)"] = r.area;
  p.notes = r.notes;
  Series sm{"force, smoothed ±0.5 Å", "displacement (Å)", "force (kcal/mol/Å)", {}, {}}, raw{"force", "displacement (Å)", "force (kcal/mol/Å)", {}, {}};
  for (size_t i = 0; i < r.curve.size(); ++i) {
    raw.x.push_back(r.curve[i].displacement);
    raw.y.push_back(r.curve[i].force);
    sm.x.push_back(r.curve[i].displacement);
    sm.y.push_back(i < r.smooth.size() ? r.smooth[i] : r.curve[i].force);
  }
  p.series = {sm, raw};
  return {p};
}

std::vector<Property> cooling_properties(const CoolingResult& r) {
  std::vector<Property> out;
  Property tg = prop("tg", "Glass transition Tg", "K", r.method, r.fit.ok ? r.fit.tg : NaN, r.fit.ok ? r.fit.tg_err : NaN);
  if (r.fit.ok) {
    tg.extra["expansion below Tg (1/K)"] = r.fit.alpha_low;
    tg.extra["expansion above Tg (1/K)"] = r.fit.alpha_high;
    tg.extra["specific volume at Tg (cm³/g)"] = r.fit.value_at_tg;
  }
  tg.notes = r.notes;
  Series v{"specific volume", "T (K)", "specific volume (cm³/g)", {}, {}}, fit{"two-line fit", "T (K)", "specific volume (cm³/g)", {}, {}};
  Series d{"density", "T (K)", "density (g/cm³)", {}, {}};
  for (const auto& p : r.points) {
    v.x.push_back(p.temperature);
    v.y.push_back(p.specific_volume);
    d.x.push_back(p.temperature);
    d.y.push_back(p.density);
  }
  if (r.fit.ok && !r.points.empty()) {
    double lo = r.points.front().temperature, hi = lo;
    for (const auto& p : r.points) { lo = std::min(lo, p.temperature); hi = std::max(hi, p.temperature); }
    for (int k = 0; k <= 100; ++k) {
      const double T = lo + (hi - lo) * k / 100;
      fit.x.push_back(T);
      fit.y.push_back(r.fit.value_at_tg + r.fit.slope_low * std::min(T - r.fit.tg, 0.0) + r.fit.slope_high * std::max(T - r.fit.tg, 0.0));
    }
  }
  tg.series = {v, fit, d};
  out.push_back(tg);
  return out;
}

}  // namespace caps
