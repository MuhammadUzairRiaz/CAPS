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
  // the compliance and what it gives along the cell axes: E_i = 1/S_ii, ν_ij = −S_ij/S_ii (strain along j over strain
  // along i, stretched along i), G_yz = 1/S44 …; the universal anisotropy index A^U = 5 G_V/G_R + K_V/K_R − 6 (0: isotropic;
  // Ranganathan & Ostoja-Starzewski 2008)
  const bool hasS = r.S[0][0] != 0;
  if (hasS) {
    for (int I = 0; I < 6; ++I)
      for (int J = I; J < 6; ++J) {
        char k[48];
        std::snprintf(k, sizeof k, "S%d%d (1/GPa)", I + 1, J + 1);
        c.extra[k] = r.S[I][J];
      }
    const char* ax[3] = {"x", "y", "z"};
    for (int i = 0; i < 3; ++i) c.extra[std::string("E") + ax[i] + " = 1/S" + std::to_string(i + 1) + std::to_string(i + 1) + " (GPa)"] = 1 / r.S[i][i];
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        if (i != j) c.extra[std::string("ν") + ax[i] + ax[j] + " = −S" + std::to_string(i + 1) + std::to_string(j + 1) + "/S" + std::to_string(i + 1) + std::to_string(i + 1)] = -r.S[i][j] / r.S[i][i];
    const char* sh[3] = {"yz", "xz", "xy"};
    for (int i = 0; i < 3; ++i) c.extra[std::string("G") + sh[i] + " = 1/S" + std::to_string(i + 4) + std::to_string(i + 4) + " (GPa)"] = 1 / r.S[i + 3][i + 3];
    if (r.G_reuss > 0 && r.K_reuss > 0) c.extra["anisotropy index A^U"] = 5 * r.G_voigt / r.G_reuss + r.K_voigt / r.K_reuss - 6;
  }
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
  nu.extra["Lamé μ = G (GPa)"] = r.G_hill;
  out.push_back(nu);
  // sound speeds of the isotropic (Hill) solid: v_L = √((K + 4G/3)/ρ), v_T = √(G/ρ), and the Debye mean
  // v_m = [(2/v_T³ + 1/v_L³)/3]^(−1/3) (Anderson 1963); GPa over g/cm³ is 10⁶ m²/s²
  if (r.density > 0 && r.K_hill > 0 && r.G_hill > 0) {
    const double vl = 1000 * std::sqrt((r.K_hill + 4 * r.G_hill / 3) / r.density), vt = 1000 * std::sqrt(r.G_hill / r.density);
    Property v = prop("sound" + suffix, "Sound speeds" + how, "m/s", "longitudinal v_L = √((K + 4G/3)/ρ) from the Hill K and G" + how, vl);
    v.extra["v_L longitudinal (m/s)"] = vl;
    v.extra["v_T transverse (m/s)"] = vt;
    v.extra["v_m Debye mean (m/s)"] = std::pow((2 / (vt * vt * vt) + 1 / (vl * vl * vl)) / 3, -1.0 / 3);
    v.extra["ρ (g/cm³)"] = r.density;
    v.extra["P-wave modulus M = K + 4G/3 (GPa)"] = r.K_hill + 4 * r.G_hill / 3;
    out.push_back(v);
  }
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
  const bool energy = r.property == 1;
  const std::string what = energy ? "potential energy per atom" : "specific volume", unit = energy ? "kcal/mol" : "cm³/g";
  if (r.fit.ok) {
    if (!energy) {
      tg.extra["expansion below Tg (1/K)"] = r.fit.alpha_low;
      tg.extra["expansion above Tg (1/K)"] = r.fit.alpha_high;
    } else {
      tg.extra["heat capacity (potential part) below Tg (kcal/mol/K per atom)"] = r.fit.slope_low;
      tg.extra["heat capacity (potential part) above Tg (kcal/mol/K per atom)"] = r.fit.slope_high;
    }
    tg.extra[what + " at Tg (" + unit + ")"] = r.fit.value_at_tg;
    tg.extra["slope below Tg (" + unit + "/K)"] = r.fit.slope_low;
    tg.extra["slope above Tg (" + unit + "/K)"] = r.fit.slope_high;
    if (!r.points.empty()) tg.extra["fit residual rms (" + unit + ")"] = std::sqrt(r.fit.rss / double(r.points.size()));
  }
  tg.notes = r.notes;
  Series v{what, "T (K)", what + " (" + unit + ")", {}, {}}, fit{"fit", "T (K)", what + " (" + unit + ")", {}, {}};
  Series d{"density", "T (K)", "density (g/cm³)", {}, {}};
  for (size_t k = 0; k < r.points.size(); ++k) {
    const auto& p = r.points[k];
    v.x.push_back(p.temperature);
    v.y.push_back(k < r.fitted.size() ? r.fitted[k] : p.specific_volume);
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
