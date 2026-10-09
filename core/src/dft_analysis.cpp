// Analysis of finished VASP runs (caps/dft_analysis.hpp).
#include "caps/dft_analysis.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <deque>
#include <filesystem>
#include <fstream>
#include <map>
#include <numeric>
#include <regex>
#include <set>
#include <sstream>
#include <stdexcept>

#include "caps/elements.hpp"
#include "caps/slab2d.hpp"
#include "caps/vasp_out.hpp"

namespace caps {

namespace fs = std::filesystem;

namespace {

std::string slurp(const std::string& p) {
  std::ifstream in(p, std::ios::binary);
  if (!in) return "";
  std::stringstream ss;
  ss << in.rdbuf();
  return ss.str();
}
void spit(const std::string& p, const std::string& t) { std::ofstream f(p, std::ios::binary); f << t; }
std::string join(const std::string& a, const std::string& b) { return (fs::path(a) / b).string(); }
std::string fmt(const char* f, double v) { char b[96]; std::snprintf(b, sizeof b, f, v); return b; }
bool exists(const std::string& p) { std::error_code ec; return fs::exists(p, ec); }

Vec3 mic_ab(const Cell& c, Vec3 d) {
  const Vec3 f = c.to_fractional(d + c.origin);
  const double r0 = std::round(f[0]), r1 = std::round(f[1]);
  Vec3 best = d;
  double bl = 1e300;
  for (int i = -1; i <= 1; ++i)
    for (int j = -1; j <= 1; ++j) {
      const Vec3 e = d - c.a * (r0 + i) - c.b * (r1 + j);
      if (dot(e, e) < bl) bl = dot(e, e), best = e;
    }
  return best;
}

std::pair<std::vector<size_t>, std::vector<size_t>> split_molecule(const System& at, const System& slab) {
  std::map<int, int> sc, seen;
  for (const auto& a : slab.atoms) ++sc[a.element];
  std::vector<size_t> mol, sl;
  for (size_t i = 0; i < at.atoms.size(); ++i) (++seen[at.atoms[i].element] > sc[at.atoms[i].element] ? mol : sl).push_back(i);
  return {mol, sl};
}

double top_plane(const System& at, const std::vector<size_t>& slab) {
  double zx = -1e300;
  for (size_t i : slab) if (at.atoms[i].element != 1) zx = std::max(zx, at.atoms[i].pos[2]);
  double sum = 0;
  int n = 0;
  for (size_t i : slab) if (at.atoms[i].element != 1 && at.atoms[i].pos[2] > zx - 1.0) sum += at.atoms[i].pos[2], ++n;
  return n ? sum / n : zx;
}

// the molecule made whole along bonds < 1.9 Å
std::vector<Vec3> whole(const System& at, const std::vector<size_t>& idx) {
  std::vector<Vec3> out(idx.size());
  if (idx.empty()) return out;
  std::vector<char> done(idx.size(), 0);
  out[0] = at.atoms[idx[0]].pos;
  done[0] = 1;
  std::deque<size_t> q{0};
  while (!q.empty()) {
    const size_t i = q.front();
    q.pop_front();
    for (size_t j = 0; j < idx.size(); ++j) {
      if (done[j]) continue;
      const Vec3 D = mic_ab(at.cell, at.atoms[idx[j]].pos - at.atoms[idx[i]].pos);
      if (norm(D) < 1.9) out[j] = out[i] + D, done[j] = 1, q.push_back(j);
    }
  }
  for (size_t j = 0; j < idx.size(); ++j) if (!done[j]) out[j] = at.atoms[idx[j]].pos;
  return out;
}

std::string structure_in(const std::string& d) {
  for (const char* f : {"CONTCAR", "POSCAR.in", "POSCAR"}) {
    const std::string p = join(d, f);
    std::error_code ec;
    if (exists(p) && fs::file_size(p, ec) > 0) return p;
  }
  throw std::runtime_error("no structure in " + d);
}

std::string find_slab_poscar(const std::string& dir) {
  fs::path d = fs::absolute(dir);
  for (int k = 0; k < 3; ++k) {
    const fs::path p = d.parent_path() / "slab" / "POSCAR";
    if (exists(p.string())) return p.string();
    d = d.parent_path();
  }
  throw std::runtime_error("cannot find slab/POSCAR next to " + dir);
}

std::map<std::string, std::string> incar_settings(const std::string& path) {
  static const std::set<std::string> keys = {"ENCUT", "IVDW", "ISMEAR", "SIGMA", "PREC", "LDIPOL", "IDIPOL", "LREAL", "GGA"};
  std::map<std::string, std::string> out;
  std::stringstream ss(slurp(path));
  for (std::string l; std::getline(ss, l);) {
    l = l.substr(0, l.find('#'));
    const auto e = l.find('=');
    if (e == std::string::npos) continue;
    std::string k = l.substr(0, e), v = l.substr(e + 1);
    auto trim = [](std::string& s) { s.erase(0, s.find_first_not_of(" \t")); s.erase(s.find_last_not_of(" \t\r") + 1); };
    trim(k), trim(v);
    if (keys.count(k)) out[k] = v;
  }
  return out;
}

double linfit_slope(const std::vector<double>& x, const std::vector<double>& y) {
  const double n = double(x.size());
  double mx = 0, my = 0;
  for (size_t i = 0; i < x.size(); ++i) mx += x[i] / n, my += y[i] / n;
  double sxy = 0, sxx = 0;
  for (size_t i = 0; i < x.size(); ++i) sxy += (x[i] - mx) * (y[i] - my), sxx += (x[i] - mx) * (x[i] - mx);
  return sxx > 0 ? sxy / sxx : 0;
}

std::vector<std::string> subdirs(const std::string& d, const std::string& prefix) {
  std::vector<std::string> out;
  std::error_code ec;
  for (const auto& e : fs::directory_iterator(d, ec))
    if (e.is_directory() && e.path().filename().string().rfind(prefix, 0) == 0) out.push_back(e.path().string());
  std::sort(out.begin(), out.end());
  return out;
}

}  // namespace

// ---------------------------------------------------------------- binding energies
Json binding_energies(const std::string& root) {
  Json r = Json::object();
  auto energy = [](const std::string& folder, bool* ok, std::string* st) {
    for (const auto& [s, need] : std::vector<std::pair<std::string, std::string>>{{"04_static", "EDIFF is reached"}, {"03_relax", "reached required accuracy"}}) {
      const auto o = read_outcar(join(join(folder, s), "OUTCAR"));
      if (!o.e0.empty()) { *ok = s == "04_static" ? o.ediff_reached : o.relaxed; *st = s; return o.e0.back(); }
    }
    *ok = false;
    *st = "";
    return std::nan("");
  };
  bool oks, okm;
  std::string sts, stm;
  const double es = energy(join(root, "slab"), &oks, &sts), em = energy(join(root, "molecule"), &okm, &stm);
  if (!std::isfinite(es) || !std::isfinite(em)) { r["error"] = "slab or molecule has no finished OUTCAR yet"; return r; }
  r["E_slab"] = es, r["E_slab_stage"] = sts, r["E_slab_converged"] = oks;
  r["E_molecule"] = em, r["E_molecule_stage"] = stm, r["E_molecule_converged"] = okm;
  const auto ref = incar_settings(join(join(root, "slab"), "INCAR.static"));
  Json warn = Json::array();
  std::vector<std::string> dirs = subdirs(root, "complex_");
  for (const auto& d : std::vector<std::string>{join(root, "molecule")}) dirs.insert(dirs.begin(), d);
  for (const auto& d : dirs)
    if (incar_settings(join(d, "INCAR.static")) != ref) warn.push_back(fs::path(d).filename().string() + " uses different settings than slab/ (ENCUT/IVDW/ISMEAR/SIGMA/PREC/IDIPOL/LREAL/GGA)");
  r["settings_warnings"] = warn;
  struct Row { std::string name; double ec, eb; bool final_; };
  std::vector<Row> rows;
  for (const auto& d : subdirs(root, "complex_")) {
    bool ok;
    std::string st;
    const double ec = energy(d, &ok, &st);
    if (!std::isfinite(ec)) continue;
    rows.push_back({fs::path(d).filename().string(), ec, ec - es - em, ok && oks && okm && st == "04_static" && sts == "04_static" && stm == "04_static"});
  }
  std::sort(rows.begin(), rows.end(), [](const Row& a, const Row& b) { return a.eb < b.eb; });
  Json rj = Json::array();
  std::string csv = "complex,E_complex_eV,E_bind_eV,E_bind_kJmol,final\n";
  for (const auto& x : rows) {
    Json o = Json::object();
    o["complex"] = x.name, o["E_complex"] = x.ec, o["E_bind"] = x.eb, o["E_bind_kJmol"] = x.eb * 96.485, o["final"] = x.final_;
    rj.push_back(o);
    csv += x.name + "," + fmt("%.5f", x.ec) + "," + fmt("%.5f", x.eb) + "," + fmt("%.2f", x.eb * 96.485) + "," + (x.final_ ? "True" : "False") + "\n";
  }
  spit(join(root, "binding_energies.csv"), csv);
  r["complexes"] = rj;
  if (!rows.empty()) {
    r["best"] = rows[0].name, r["best_E_bind"] = rows[0].eb;
    r["spread_best3"] = rows[std::min<size_t>(2, rows.size() - 1)].eb - rows[0].eb;
    Json close = Json::array();
    for (size_t k = 1; k < rows.size(); ++k) if (rows[k].eb - rows[0].eb < 0.05) close.push_back(rows[k].name);
    r["within_0.05_eV"] = close;
    if (rows[0].eb > 0) r["note"] = "E_bind > 0 means no binding: check the contact distance, dispersion (IVDW) and convergence";
  }
  return r;
}

// ---------------------------------------------------------------- geometry
Json adsorption_geometry_of(const System& at, const System& slab_poscar) {
  const auto [mol, slab] = split_molecule(at, slab_poscar);
  Json r = Json::object();
  if (mol.empty()) { r["error"] = "no molecule atoms (the structure has the slab's composition)"; return r; }
  const auto m = whole(at, mol);
  const double plane = top_plane(at, slab);
  double hmin = 1e300, M = 0;
  Vec3 com{0, 0, 0};
  std::vector<size_t> heavy;
  for (size_t k = 0; k < mol.size(); ++k) {
    const int z = at.atoms[mol[k]].element;
    if (z != 1) heavy.push_back(k), hmin = std::min(hmin, m[k][2] - plane);
    com = com + m[k] * element(z).mass;
    M += element(z).mass;
  }
  com = com * (1.0 / M);
  r["h_min"] = hmin - 0.0;
  r["h_com"] = com[2] - plane;
  // the anchor: the nitrile N (an N bonded to C under 1.25 Å), else the lowest heavy atom
  size_t kA = heavy.empty() ? 0 : heavy[0];
  bool nitrile = false;
  for (size_t k = 0; k < mol.size() && !nitrile; ++k)
    if (at.atoms[mol[k]].element == 7)
      for (size_t q = 0; q < mol.size(); ++q) if (at.atoms[mol[q]].element == 6 && norm(m[q] - m[k]) < 1.25) { kA = k; nitrile = true; break; }
  if (!nitrile) for (size_t k : heavy) if (m[k][2] < m[kA][2]) kA = k;
  double dA = 1e300;
  size_t partner = 0;
  for (size_t j : slab) { const double d = norm(mic_ab(at.cell, at.atoms[j].pos - m[kA])); if (d < dA) dA = d, partner = j; }
  r["anchor"] = std::string(element(at.atoms[mol[kA]].element).symbol) + (nitrile ? " (nitrile)" : " (lowest heavy atom)");
  r["d_anchor"] = dA;
  r["d_anchor_partner"] = std::string(element(at.atoms[partner].element).symbol);
  // O–H···X hydrogen bonds from surface OH to the anchor
  int nh = 0;
  for (size_t h : slab) {
    if (at.atoms[h].element != 1) continue;
    const Vec3 vX = mic_ab(at.cell, m[kA] - at.atoms[h].pos);
    if (norm(vX) >= 2.5) continue;
    double best = 1e300;
    Vec3 vO{};
    for (size_t o : slab) if (at.atoms[o].element == 8) { const Vec3 v = mic_ab(at.cell, at.atoms[o].pos - at.atoms[h].pos); if (norm(v) < best) best = norm(v), vO = v; }
    if (best > 1.3) continue;
    const double c = dot(vO, vX) / (norm(vO) * norm(vX));
    if (std::acos(std::clamp(c, -1.0, 1.0)) * 180 / M_PI > 120) ++nh;
  }
  r["Hbond_OH_anchor"] = nh;
  int nch = 0, nc = 0;
  for (size_t k = 0; k < mol.size(); ++k) {
    double dacc = 1e300, dany = 1e300;
    for (size_t j : slab) {
      const double d = norm(mic_ab(at.cell, at.atoms[j].pos - m[k]));
      dany = std::min(dany, d);
      if (at.atoms[j].element == 8 || at.atoms[j].element == 9) dacc = std::min(dacc, d);
    }
    if (at.atoms[mol[k]].element == 1 && dacc < 2.7) ++nch;
    nc += dany < 3.5;
  }
  r["CH_X_contacts"] = nch;
  r["n_contacts"] = nc;
  // tilt of the long axis of the heavy atoms (largest principal axis of their spread) against the surface plane
  if (heavy.size() >= 2) {
    Vec3 c{0, 0, 0};
    for (size_t k : heavy) c = c + m[k];
    c = c * (1.0 / double(heavy.size()));
    double S[3][3] = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};
    for (size_t k : heavy) { const Vec3 d = m[k] - c; for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b) S[a][b] += d[size_t(a)] * d[size_t(b)]; }
    // power iteration for the largest axis
    Vec3 v{1, 0.3, 0.2};
    for (int it = 0; it < 200; ++it) {
      Vec3 w{0, 0, 0};
      for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b) w[size_t(a)] += S[a][b] * v[size_t(b)];
      v = w * (1.0 / norm(w));
    }
    r["tilt_deg"] = std::asin(std::fabs(v[2])) * 180 / M_PI;
  }
  return r;
}

Json adsorption_geometry(const std::vector<std::string>& dirs) {
  Json rows = Json::array();
  std::string csv = "complex,h_min,h_com,d_anchor,partner,Hbonds,CH_X,n_contacts,tilt_deg\n";
  std::string parent;
  for (const auto& d : dirs) {
    const std::string stage = fs::is_directory(join(d, "03_relax")) ? join(d, "03_relax") : d;
    const System at = read_vasp_poscar(structure_in(stage));
    Json r = adsorption_geometry_of(at, read_vasp_poscar(find_slab_poscar(d)));
    r["complex"] = fs::path(d).filename().string();
    spit(join(d, "geometry.json"), r.dump(1));
    rows.push_back(r);
    csv += r.text("complex") + "," + fmt("%.3f", r.num("h_min", 0)) + "," + fmt("%.3f", r.num("h_com", 0)) + "," + fmt("%.3f", r.num("d_anchor", 0)) + "," +
           r.text("d_anchor_partner") + "," + fmt("%.0f", r.num("Hbond_OH_anchor", 0)) + "," + fmt("%.0f", r.num("CH_X_contacts", 0)) + "," +
           fmt("%.0f", r.num("n_contacts", 0)) + "," + fmt("%.1f", r.num("tilt_deg", 0)) + "\n";
    parent = fs::absolute(d).parent_path().string();
  }
  if (!parent.empty()) spit(join(parent, "geometry_table.csv"), csv);
  Json out = Json::object();
  out["complexes"] = rows;
  return out;
}

// ---------------------------------------------------------------- work function
Json work_function(const std::string& d, const std::string& ref) {
  Json r = Json::object();
  r["folder"] = d;
  if (!exists(join(d, "LOCPOT"))) { r["error"] = "no LOCPOT (static run without LVHAR?)"; return r; }
  const auto g = read_vasp_grid(join(d, "LOCPOT"));
  const auto o = read_outcar(join(d, "OUTCAR"));
  if (o.efermi.empty()) { r["error"] = "no E-fermi in OUTCAR"; return r; }
  const double ef = o.efermi.back();
  const int nz = g.n[2];
  std::vector<double> prof(size_t(nz), 0.0);
  for (int k = 0; k < nz; ++k) {
    double s = 0;
    for (int j = 0; j < g.n[1]; ++j) for (int i = 0; i < g.n[0]; ++i) s += g.at(i, j, k);
    prof[size_t(k)] = s / double(g.n[0] * g.n[1]);
  }
  const double c = g.system.cell.c[2];
  double zmin = 1e300, zmax = -1e300;
  for (const auto& a : g.system.atoms) { const double z = a.pos[2] - std::floor(a.pos[2] / c) * c; zmin = std::min(zmin, z), zmax = std::max(zmax, z); }
  std::vector<double> top, bot;
  for (int k = 0; k < nz; ++k) {
    const double z = c * k / nz;
    if (z > zmax + 4.0 && z < c - 1.0) top.push_back(prof[size_t(k)]);
    if (z > 1.0 && z < zmin - 4.0) bot.push_back(prof[size_t(k)]);
  }
  auto median = [](std::vector<double> v) { std::sort(v.begin(), v.end()); return v.empty() ? std::nan("") : v[v.size() / 2] * (v.size() % 2 ? 1 : 0.5) + (v.size() % 2 ? 0 : 0.5 * v[v.size() / 2 - 1]); };
  r["E_fermi"] = ef;
  if (!top.empty()) {
    r["V_vac_top"] = median(top), r["phi_top"] = median(top) - ef;
    r["plateau_flatness_top"] = *std::max_element(top.begin(), top.end()) - *std::min_element(top.begin(), top.end());
  }
  if (!bot.empty()) r["V_vac_bottom"] = median(bot), r["phi_bottom"] = median(bot) - ef;
  Json pj = Json::array();
  for (int k = 0; k < nz; ++k) { Json p = Json::array(); p.push_back(c * k / nz); p.push_back(prof[size_t(k)] - ef); pj.push_back(p); }
  r["profile"] = pj;
  if (!ref.empty()) {
    const Json rr = work_function(ref);
    if (rr.has("phi_top") && r.has("phi_top")) r["delta_phi"] = r["phi_top"].number() - rr["phi_top"].number();
  }
  if (o.dipole.size()) r["dipole"] = o.dipole.back();
  Json slim = r;
  slim["profile"] = Json::array();
  spit(join(d, "workfunction.json"), slim.dump(1));
  return r;
}

// ---------------------------------------------------------------- DOS
Json dos_analysis(const std::string& d, const std::string& slab_poscar) {
  Json r = Json::object();
  r["folder"] = d;
  const Doscar ds = read_doscar(join(d, "DOSCAR"));
  const System at = read_vasp_poscar(structure_in(d));
  std::vector<size_t> mol, slab;
  const bool free_mol = fs::absolute(d).string().find("/molecule/") != std::string::npos;
  if (free_mol) for (size_t i = 0; i < at.atoms.size(); ++i) mol.push_back(i);
  else {
    std::string sp = slab_poscar;
    if (sp.empty()) try { sp = find_slab_poscar(d); } catch (...) {}
    if (!sp.empty()) std::tie(mol, slab) = split_molecule(at, read_vasp_poscar(sp));
    else for (size_t i = 0; i < at.atoms.size(); ++i) slab.push_back(i);
  }
  const size_t ne = ds.energy.size();
  auto group = [&](const std::vector<size_t>& idx) {
    std::vector<double> g(ne, 0.0);
    for (size_t i : idx) if (i < ds.site.size()) for (size_t e = 0; e < ne && e < ds.site[i].size(); ++e) g[e] += ds.site[i][e];
    return g;
  };
  std::map<std::string, std::vector<double>> groups;
  groups["total"] = ds.total;
  // the surface metal: the most common non-C/N/O/F/H element of the slab
  std::map<int, int> cnt;
  for (size_t i : slab) { const int z = at.atoms[i].element; if (z != 1 && z != 6 && z != 7 && z != 8 && z != 9) ++cnt[z]; }
  int metal = 0, mc = 0;
  for (const auto& [z, c] : cnt) if (c > mc) mc = c, metal = z;
  if (metal) { std::vector<size_t> mi; for (size_t i : slab) if (at.atoms[i].element == metal) mi.push_back(i); groups[element(metal).symbol] = group(mi); }
  if (!slab.empty()) {
    const double plane = top_plane(at, slab);
    std::vector<size_t> term;
    for (size_t i : slab) if ((at.atoms[i].element == 8 || at.atoms[i].element == 9) && at.atoms[i].pos[2] > plane - 0.5) term.push_back(i);
    groups["top_termination"] = group(term);
  }
  if (!mol.empty()) {
    groups["molecule"] = group(mol);
    std::vector<size_t> n;
    for (size_t i : mol) if (at.atoms[i].element == 7) n.push_back(i);
    if (!n.empty()) groups["anchor_N"] = group(n);
  }
  size_t i0 = 0;
  for (size_t e = 0; e < ne; ++e) if (std::fabs(ds.energy[e]) < std::fabs(ds.energy[i0])) i0 = e;
  r["E_fermi"] = ds.efermi;
  r["DOS_at_EF_total"] = ds.total.empty() ? 0.0 : ds.total[i0];
  auto peak = [&](const std::vector<double>& y, int sign) -> double {
    double ymax = 0;
    for (size_t e = 0; e < ne; ++e) if (sign < 0 ? ds.energy[e] < -0.05 : ds.energy[e] > 0.05) ymax = std::max(ymax, y[e]);
    if (ymax <= 0) return std::nan("");
    double best = std::nan("");
    for (size_t e = 1; e + 1 < ne; ++e) {
      if (!(sign < 0 ? ds.energy[e] < -0.05 : ds.energy[e] > 0.05)) continue;
      if (y[e] > 0.2 * ymax && y[e] >= y[e - 1] && y[e] >= y[e + 1]) {
        if (sign < 0) best = std::isnan(best) ? ds.energy[e] : std::max(best, ds.energy[e]);
        else if (std::isnan(best)) best = ds.energy[e];
      }
    }
    return best;
  };
  if (groups.count("molecule")) {
    const double h = peak(groups["molecule"], -1), l = peak(groups["molecule"], +1);
    if (std::isfinite(h)) r["mol_HOMO_like"] = h;
    if (std::isfinite(l)) r["mol_LUMO_like"] = l;
    if (groups.count("anchor_N")) { const double n = peak(groups["anchor_N"], -1); if (std::isfinite(n)) r["N_peak_below_EF"] = n; }
  }
  Json gj = Json::object();
  std::string csv = "E_minus_EF";
  for (const auto& [k, v] : groups) csv += "," + k;
  csv += "\n";
  for (size_t e = 0; e < ne; ++e) {
    csv += fmt("%.5f", ds.energy[e]);
    for (const auto& [k, v] : groups) csv += "," + fmt("%.5f", v[e]);
    csv += "\n";
  }
  for (const auto& [k, v] : groups) {
    Json a = Json::array();
    for (size_t e = 0; e < ne; ++e) a.push_back(v[e]);
    gj[k] = a;
  }
  Json ej = Json::array();
  for (double e : ds.energy) ej.push_back(e);
  Json slim = r;
  spit(join(d, "dos.json"), slim.dump(1));
  spit(join(d, "dos_groups.csv"), csv);
  r["energy"] = ej;
  r["groups"] = gj;
  return r;
}

// ---------------------------------------------------------------- charge-density difference
Json charge_density_difference(const std::string& cx, const std::string& slab_poscar) {
  const std::string pc = join(join(join(cx, "charge"), "04_static"), "CHGCAR"), ps = join(join(join(cx, "cdd_slab"), "04_static"), "CHGCAR"),
                    pm = join(join(join(cx, "cdd_mol"), "04_static"), "CHGCAR");
  for (const auto& p : {pc, ps, pm}) if (!exists(p)) throw std::runtime_error("missing " + p);
  const auto c = read_vasp_grid(pc), s = read_vasp_grid(ps), m = read_vasp_grid(pm);
  if (c.n != s.n || c.n != m.n) throw std::runtime_error("FFT grids differ: settings or cell not identical");
  VaspGrid diff = c;
  for (size_t i = 0; i < diff.values.size(); ++i) diff.values[i] = c.values[i] - s.values[i] - m.values[i];
  write_chgcar(diff, join(cx, "CHGCAR_diff"));
  const Cell& cell = c.system.cell;
  const double vol = std::fabs(dot(cell.a, cross(cell.b, cell.c)));
  const double area = norm(cross(cell.a, cell.b));
  const int nz = c.n[2];
  const double dz = cell.c[2] / nz;
  std::vector<double> prof(static_cast<size_t>(nz)), cum(static_cast<size_t>(nz));
  double run = 0, mx = 0;
  for (int k = 0; k < nz; ++k) {
    double sum = 0;
    for (int j = 0; j < c.n[1]; ++j) for (int i = 0; i < c.n[0]; ++i) { const double v = diff.at(i, j, k) / vol; sum += v; mx = std::max(mx, std::fabs(v)); }
    prof[size_t(k)] = sum / double(c.n[0] * c.n[1]) * area;   // e/Å
    run += prof[size_t(k)] * dz;
    cum[size_t(k)] = run;
  }
  const System at = c.system;
  const auto [mol, slab] = split_molecule(at, read_vasp_poscar(slab_poscar));
  const double plane = top_plane(at, slab);
  double zmol = 1e300;
  for (size_t i : mol) if (at.atoms[i].element != 1) zmol = std::min(zmol, at.atoms[i].pos[2]);
  const double zmid = 0.5 * (plane + zmol);
  const int k = int(std::lround(zmid / dz)) % nz;
  Json r = Json::object();
  r["complex"] = cx, r["z_surface"] = plane, r["z_molecule_min"] = zmol, r["z_mid"] = zmid, r["dQ_at_mid"] = cum[size_t((k + nz) % nz)], r["max_abs_drho"] = mx;
  spit(join(cx, "cdd.json"), r.dump(1));
  std::string csv = "z_A,drho_e_per_A,dQ_cumulative_e\n";
  Json pj = Json::array();
  for (int q = 0; q < nz; ++q) {
    csv += fmt("%.5f", q * dz) + "," + fmt("%.6e", prof[size_t(q)]) + "," + fmt("%.6e", cum[size_t(q)]) + "\n";
    Json p = Json::array();
    p.push_back(q * dz); p.push_back(prof[size_t(q)]); p.push_back(cum[size_t(q)]);
    pj.push_back(p);
  }
  spit(join(cx, "cdd_profile.csv"), csv);
  r["profile"] = pj;
  r["note"] = "dQ(z) integrated from the bottom of the cell: dQ at the mid-plane > 0 means the slab side GAINED electrons (molecule donates); must agree in sign with Bader";
  return r;
}

// ---------------------------------------------------------------- Bader
Json bader_transfer(const std::string& slab_p, const std::string& cx_p, const std::string& acf, const std::string& potcar) {
  const auto sn = poscar_species(slab_p), cn = poscar_species(cx_p);
  std::map<std::string, int> slab_count;
  for (const auto& [s, c] : sn) slab_count[s] += c;
  // ZVAL: from the complex's POTCAR (in species order), else the valences of the reference PAW sets
  std::map<std::string, double> zval = {{"Ti", 10.0}, {"C", 4.0}, {"N", 5.0}, {"O", 6.0}, {"F", 7.0}, {"H", 1.0}};
  std::string zsource = "Ti_pv/C/N/O/F/H valences (no POTCAR given)";
  if (!potcar.empty()) {
    std::vector<double> zs;
    const std::string t = slurp(potcar);
    static const std::regex re(R"(ZVAL\s*=\s*([\d.]+))");
    for (auto it = std::sregex_iterator(t.begin(), t.end(), re); it != std::sregex_iterator(); ++it) zs.push_back(std::stod((*it)[1]));
    if (zs.size() == cn.size()) { zval.clear(); for (size_t k = 0; k < cn.size(); ++k) zval[cn[k].first] = zs[k]; zsource = "POTCAR"; }
  }
  const auto q = read_acf(acf);
  int total = 0;
  for (const auto& [s, c] : cn) total += c;
  if (int(q.size()) != total) throw std::runtime_error("ACF.dat has " + std::to_string(q.size()) + " atoms, the POSCAR " + std::to_string(total));
  double mol = 0, sl = 0;
  size_t k = 0;
  for (const auto& [s, c] : cn) {
    if (!zval.count(s)) throw std::runtime_error("no ZVAL for " + s + ": give the POTCAR");
    for (int i = 0; i < c; ++i, ++k) ((i < slab_count[s]) ? sl : mol) += q[k] - zval[s];
  }
  Json r = Json::object();
  r["molecule_gained"] = mol, r["slab_gained"] = sl, r["total"] = mol + sl, r["zval_source"] = zsource;
  if (std::fabs(mol + sl) > 0.05) r["note"] = "the total should be ~0 for a neutral system: the Bader analysis is not converged (finer grid, AECCAR0+2 reference)";
  return r;
}

// ---------------------------------------------------------------- frequencies
Json nitrile_shift(const std::string& cx) {
  const std::string molcase = join(fs::absolute(cx).parent_path().string(), "molecule");
  const auto fc = read_outcar(join(join(join(cx, "freq"), "04_static"), "OUTCAR")), fm = read_outcar(join(join(join(molcase, "freq"), "04_static"), "OUTCAR"));
  Json r = Json::object();
  if (fc.freq_real.empty() || fm.freq_real.empty()) { r["error"] = "frequency runs not finished (need <complex>/freq and molecule/freq)"; return r; }
  auto cn = [](const std::vector<double>& f) { double b = std::nan(""); for (double x : f) if (x > 2000 && x < 2500 && (std::isnan(b) || x > b)) b = x; return b; };
  const double a = cn(fc.freq_real), f = cn(fm.freq_real);
  if (std::isfinite(a)) r["nu_CN_adsorbed"] = a;
  if (std::isfinite(f)) r["nu_CN_free"] = f;
  if (std::isfinite(a) && std::isfinite(f)) r["shift"] = a - f;
  Json ic = Json::array(), im = Json::array();
  for (double x : fc.freq_imag) ic.push_back(x);
  for (double x : fm.freq_imag) im.push_back(x);
  r["imaginary_complex"] = ic, r["imaginary_free"] = im;
  if (!fc.freq_imag.empty() || !fm.freq_imag.empty()) r["note"] = "imaginary modes in a PARTIAL Hessian are common frozen-neighbour artefacts; only the C≡N stretch is interpreted";
  spit(join(cx, "freq.json"), r.dump(1));
  return r;
}

// ---------------------------------------------------------------- AIMD
Json md_analysis(const std::string& dir, double equil_ps, double dt_fs) {
  std::vector<System> frames;
  std::vector<IonicStep> en;
  for (const auto& s : subdirs(dir, "seg_")) {
    if (!exists(join(s, "XDATCAR")) || !exists(join(s, "OSZICAR"))) continue;
    auto fr = read_xdatcar(join(s, "XDATCAR"));
    auto e = read_oszicar(join(s, "OSZICAR"));
    const size_t n = std::min(fr.size(), e.size());
    frames.insert(frames.end(), fr.begin(), fr.begin() + long(n));
    en.insert(en.end(), e.begin(), e.begin() + long(n));
  }
  if (frames.size() < 10) throw std::runtime_error("fewer than 10 MD frames found in " + dir);
  const System& f0 = frames[0];
  const size_t natoms = f0.atoms.size();
  // the molecule: the fragment of C/N/H atoms (bonds < 1.75 Å) that holds an N, else the largest C/H fragment
  std::vector<size_t> cand;
  for (size_t i = 0; i < natoms; ++i) { const int z = f0.atoms[i].element; if (z == 6 || z == 7 || z == 1) cand.push_back(i); }
  std::vector<int> comp(natoms, -1);
  std::vector<std::vector<size_t>> frags;
  for (size_t c : cand) {
    if (comp[c] >= 0) continue;
    std::vector<size_t> fr{c};
    comp[c] = int(frags.size());
    for (size_t h = 0; h < fr.size(); ++h)
      for (size_t d : cand)
        if (comp[d] < 0 && norm(mic_ab(f0.cell, f0.atoms[d].pos - f0.atoms[fr[h]].pos)) < 1.75) comp[d] = int(frags.size()), fr.push_back(d);
    frags.push_back(fr);
  }
  std::vector<size_t> mol;
  for (const auto& fr : frags) {
    const bool hasN = std::any_of(fr.begin(), fr.end(), [&](size_t i) { return f0.atoms[i].element == 7; });
    if ((hasN && (mol.empty() || !std::any_of(mol.begin(), mol.end(), [&](size_t i) { return f0.atoms[i].element == 7; }))) || (!hasN && fr.size() > mol.size() && mol.empty())) mol = fr;
  }
  std::sort(mol.begin(), mol.end());
  if (mol.size() < 2) throw std::runtime_error("no molecule found among the C/N/H atoms");
  std::set<size_t> ms(mol.begin(), mol.end());
  std::vector<size_t> slab, slab_heavy;
  for (size_t i = 0; i < natoms; ++i) if (!ms.count(i)) { slab.push_back(i); if (f0.atoms[i].element != 1) slab_heavy.push_back(i); }
  size_t iA = mol[0];
  for (size_t i : mol) if (f0.atoms[i].element == 7) { iA = i; break; }
  std::vector<Vec3> ref = whole(f0, mol);
  struct B { size_t p, q; double L; };
  std::vector<B> bonds;
  for (size_t p = 0; p < mol.size(); ++p)
    for (size_t q = p + 1; q < mol.size(); ++q) {
      const double L = norm(ref[p] - ref[q]);
      if (L < 1.2 * (element(f0.atoms[mol[p]].element).covalent + element(f0.atoms[mol[q]].element).covalent)) bonds.push_back({p, q, L});
    }
  struct Row { double t, T, E, F, h, hA, dmin; int cont, broken; };
  std::vector<Row> R;
  for (size_t k = 0; k < frames.size(); ++k) {
    const System& at = frames[k];
    std::vector<Vec3> m(mol.size());
    for (size_t p = 0; p < mol.size(); ++p) m[p] = ref[p] + mic_ab(at.cell, at.atoms[mol[p]].pos - ref[p]);
    ref = m;
    double top = -1e300;
    for (size_t i : slab_heavy) top = std::max(top, at.atoms[i].pos[2]);
    double ps = 0;
    int pn = 0;
    for (size_t i : slab_heavy) if (at.atoms[i].pos[2] > top - 1.0) ps += at.atoms[i].pos[2], ++pn;
    const double plane = pn ? ps / pn : top;
    double M = 0;
    Vec3 com{0, 0, 0};
    for (size_t p = 0; p < mol.size(); ++p) { const double w = element(at.atoms[mol[p]].element).mass; com = com + m[p] * w; M += w; }
    com = com * (1.0 / M);
    const size_t pA = size_t(std::find(mol.begin(), mol.end(), iA) - mol.begin());
    double dmin = 1e300;
    int cont = 0;
    for (size_t p = 0; p < mol.size(); ++p) {
      double dany = 1e300;
      for (size_t j : slab) {
        const double d = norm(mic_ab(at.cell, at.atoms[j].pos - m[p]));
        dany = std::min(dany, d);
        if (at.atoms[mol[p]].element != 1 && at.atoms[j].element != 1) dmin = std::min(dmin, d);
      }
      cont += dany < 3.5;
    }
    int broken = 0;
    for (const auto& b : bonds) broken += norm(m[b.p] - m[b.q]) > 1.5 * b.L;
    R.push_back({double(k) * dt_fs / 1000.0, en[k].temperature, en[k].etot, en[k].f, com[2] - plane, m[pA][2] - plane, dmin, cont, broken});
  }
  std::vector<const Row*> prod;
  for (const auto& r : R) if (r.t >= equil_ps) prod.push_back(&r);
  if (prod.size() < 10) { prod.clear(); for (size_t k = R.size() / 2; k < R.size(); ++k) prod.push_back(&R[k]); }
  double Tm = 0, Ts = 0, hp = 0, hs = 0, hAp = 0, dm = 0, dmax = 0, cm = 0;
  int bmax = 0;
  std::vector<double> tx, ey;
  for (const Row* r : prod) Tm += r->T, hp += r->h, hAp += r->hA, dm += r->dmin, cm += r->cont, dmax = std::max(dmax, r->dmin), bmax = std::max(bmax, r->broken), tx.push_back(r->t), ey.push_back(r->E);
  const double np = double(prod.size());
  Tm /= np, hp /= np, hAp /= np, dm /= np, cm /= np;
  for (const Row* r : prod) Ts += (r->T - Tm) * (r->T - Tm), hs += (r->h - hp) * (r->h - hp);
  Ts = std::sqrt(Ts / np), hs = std::sqrt(hs / np);
  const double Texp = Tm * std::sqrt(2.0 / (3.0 * double(natoms)));
  const double drift = linfit_slope(tx, ey) / double(natoms) * 1000;   // meV/atom/ps
  const double h0 = R[0].h;
  const bool ok = bmax == 0 && hp <= h0 + 1.5 && dmax < 4.5;
  Json r = Json::object();
  r["frames"] = int(R.size()), r["length_ps"] = R.back().t, r["production_from_ps"] = prod.front()->t, r["atoms"] = int(natoms), r["molecule_atoms"] = int(mol.size());
  r["T_mean"] = Tm, r["T_std"] = Ts, r["T_std_expected"] = Texp, r["drift_meV_atom_ps"] = drift;
  r["height_start"] = h0, r["height_mean"] = hp, r["height_std"] = hs, r["anchor_height_start"] = R[0].hA, r["anchor_height_mean"] = hAp;
  r["contact_mean"] = dm, r["contact_max"] = dmax, r["contacts_mean"] = cm, r["broken_bonds_max"] = bmax;
  r["verdict"] = ok ? "STAYS ADSORBED" : "DESORBS / UNSTABLE";
  std::string csv = "t_ps,T_K,E_eV,F_eV,height_A,anchor_height_A,dmin_heavy_A,contacts,broken_bonds\n";
  Json ts = Json::array();
  for (const auto& x : R) {
    csv += fmt("%.4f", x.t) + "," + fmt("%.2f", x.T) + "," + fmt("%.5f", x.E) + "," + fmt("%.5f", x.F) + "," + fmt("%.4f", x.h) + "," + fmt("%.4f", x.hA) + "," +
           fmt("%.4f", x.dmin) + "," + std::to_string(x.cont) + "," + std::to_string(x.broken) + "\n";
    Json row = Json::array();
    for (double v : {x.t, x.T, x.E, x.F, x.h, x.hA, x.dmin, double(x.cont), double(x.broken)}) row.push_back(v);
    ts.push_back(row);
  }
  spit(join(dir, "md_timeseries.csv"), csv);
  std::string txt = "AIMD analysis of " + dir + "\n";
  txt += "temperature      : " + fmt("%.1f", Tm) + " +/- " + fmt("%.1f", Ts) + " K   (expected fluctuation: +/- " + fmt("%.1f", Texp) + " K)\n";
  txt += "energy drift (E) : " + fmt("%+.3f", drift) + " meV/atom/ps\n";
  txt += "COM height       : start " + fmt("%.2f", h0) + " A, production mean " + fmt("%.2f", hp) + " +/- " + fmt("%.2f", hs) + " A\n";
  txt += "closest contact  : production mean " + fmt("%.2f", dm) + " A, max " + fmt("%.2f", dmax) + " A\n";
  txt += "broken bonds     : max " + std::to_string(bmax) + "\n";
  txt += "VERDICT          : " + r.text("verdict") + "\n";
  spit(join(dir, "md_summary.txt"), txt);
  r["timeseries"] = ts;
  return r;
}

// ---------------------------------------------------------------- summary
Json summarize_results(const std::string& main) {
  Json rows = Json::array();
  const std::vector<std::string> cases = subdirs(join(main, "structures"), "");
  const std::vector<std::string> sets = subdirs(join(main, "adsorption"), "");
  std::string csv = "case,a_relaxed_A,X_bond_A,O_H_A,E0_slab_eV,best_complex,E_bind_eV,E_bind_kJmol,spread_best3_eV,bader_molecule_e,h_min_A,d_anchor_A,phi_slab_eV,delta_phi_eV,DOS_EF,cdd_dQ_mid_e,nu_shift_cm1,aimd_verdict\n";
  auto num = [](const Json& j, const char* k) { return j.has(k) && j[k].is_number() ? fmt("%.4f", j[k].number()) : std::string("n/a"); };
  auto load = [](const std::string& p) { const std::string t = slurp(p); return t.empty() ? Json::object() : Json::parse(t); };
  for (const auto& c : cases) {
    const std::string name = fs::path(c).filename().string();
    Json r = Json::object();
    r["case"] = name;
    // lattice constant per 1 × 1 cell: |a| of the relaxed cell over the in-plane repeat of the POSCAR's metal layer
    for (const char* st : {"02_cell2", "03_relax"}) {
      const std::string p = join(join(c, st), "CONTCAR");
      if (!exists(p)) continue;
      const System s = read_vasp_poscar(p);
      // repeat: atoms of the most common element per layer → n² cells
      std::map<int, int> cnt;
      for (const auto& a : s.atoms) ++cnt[a.element];
      const auto sites = surface_sites(s);
      int top = 0;
      for (const auto& x : sites) top += x.kind == "top" && x.face == "top";
      const double n = std::sqrt(std::max(1, top));
      r["a"] = norm(s.cell.a) / std::round(n);
      // mean surface–X and O–H bonds
      double sx = 0, soh = 0;
      int nx = 0, noh = 0;
      for (const auto& p2 : neighbours_pbc(s, 2.6)) {
        const int zi = s.atoms[p2.i].element, zj = s.atoms[p2.j].element;
        if ((zi == 8 || zi == 9) && zj != 1 && zj != 6 && zj != 7 && zj != 8 && zj != 9 && p2.r < 2.6) sx += p2.r, ++nx;
        if (zi == 8 && zj == 1 && p2.r < 1.2) soh += p2.r, ++noh;
      }
      if (nx) r["X_bond"] = sx / nx;
      if (noh) r["O_H"] = soh / noh;
      break;
    }
    const auto os = read_outcar(join(join(c, "04_static"), "OUTCAR"));
    if (!os.e0.empty()) r["E0_slab"] = os.e0.back();
    // the adsorption set of this case: the longest case name it starts with
    for (const auto& s : sets) {
      const std::string sn = fs::path(s).filename().string();
      if (sn.rfind(name + "_", 0) != 0) continue;
      bool longer = false;
      for (const auto& c2 : cases) { const std::string n2 = fs::path(c2).filename().string(); longer |= n2.size() > name.size() && sn.rfind(n2 + "_", 0) == 0; }
      if (longer) continue;
      const Json b = binding_energies(s);
      if (!b.has("best")) continue;
      const std::string best = b.text("best");
      r["set"] = sn, r["best_complex"] = best, r["E_bind"] = b["best_E_bind"].number(), r["E_bind_kJmol"] = b["best_E_bind"].number() * 96.485, r["spread_best3"] = b["spread_best3"].number();
      const std::string bd = join(s, best);
      const std::string chg = join(join(bd, "charge"), "04_static");
      try {
        if (exists(join(chg, "ACF.dat"))) r["bader_molecule"] = bader_transfer(join(join(s, "slab"), "POSCAR"), join(chg, "POSCAR.in"), join(chg, "ACF.dat"), join(bd, "POTCAR"))["molecule_gained"].number();
      } catch (...) {}
      const Json g = load(join(bd, "geometry.json"));
      if (g.has("h_min")) r["h_min"] = g["h_min"].number();
      if (g.has("d_anchor")) r["d_anchor"] = g["d_anchor"].number();
      const Json ws = load(join(join(join(s, "slab"), "04_static"), "workfunction.json")), wc = load(join(chg, "workfunction.json"));
      if (ws.has("phi_top")) r["phi_slab"] = ws["phi_top"].number();
      if (ws.has("phi_top") && wc.has("phi_top")) r["delta_phi"] = wc["phi_top"].number() - ws["phi_top"].number();
      const Json dj = load(join(join(bd, "04_static"), "dos.json"));
      if (dj.has("DOS_at_EF_total")) r["DOS_EF"] = dj["DOS_at_EF_total"].number();
      const Json cj = load(join(bd, "cdd.json"));
      if (cj.has("dQ_at_mid")) r["cdd_dQ_mid"] = cj["dQ_at_mid"].number();
      const Json fj = load(join(bd, "freq.json"));
      if (fj.has("shift")) r["nu_shift"] = fj["shift"].number();
    }
    const std::string md = slurp(join(join(join(main, "aimd"), name + "_best"), "md_summary.txt"));
    const auto vp = md.find("VERDICT");
    if (vp != std::string::npos) { std::string v = md.substr(md.find(':', vp) + 1); v = v.substr(0, v.find('\n')); v.erase(0, v.find_first_not_of(' ')); r["aimd_verdict"] = v; }
    rows.push_back(r);
    csv += name + "," + num(r, "a") + "," + num(r, "X_bond") + "," + num(r, "O_H") + "," + num(r, "E0_slab") + "," + r.text("best_complex", "n/a") + "," + num(r, "E_bind") + "," +
           num(r, "E_bind_kJmol") + "," + num(r, "spread_best3") + "," + num(r, "bader_molecule") + "," + num(r, "h_min") + "," + num(r, "d_anchor") + "," + num(r, "phi_slab") + "," +
           num(r, "delta_phi") + "," + num(r, "DOS_EF") + "," + num(r, "cdd_dQ_mid") + "," + num(r, "nu_shift") + "," + r.text("aimd_verdict", "n/a") + "\n";
  }
  spit(join(main, "results_summary.csv"), csv);
  Json out = Json::object();
  out["rows"] = rows;
  std::vector<std::pair<double, std::string>> rank;
  for (const auto& r : rows.items()) if (r.has("E_bind")) rank.push_back({r["E_bind"].number(), r.text("case")});
  std::sort(rank.begin(), rank.end());
  std::string order;
  for (const auto& [e, n] : rank) order += (order.empty() ? "" : " > ") + n;
  out["ranking"] = order;
  std::string txt = "binding strength ranking (strongest first): " + order + "\n";
  spit(join(main, "results_summary.txt"), csv + "\n" + txt);
  return out;
}

}  // namespace caps
