// Reading VASP output and the checks built on it (caps/vasp_out.hpp).
#include "caps/vasp_out.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <map>
#include <regex>
#include <set>
#include <sstream>
#include <stdexcept>

#include "caps/elements.hpp"
#include "caps/slab2d.hpp"

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
std::vector<std::string> lines_of(const std::string& t) {
  std::vector<std::string> out;
  std::stringstream ss(t);
  for (std::string l; std::getline(ss, l);) { if (!l.empty() && l.back() == '\r') l.pop_back(); out.push_back(l); }
  return out;
}
std::vector<std::string> split(const std::string& s) {
  std::stringstream ss(s);
  std::vector<std::string> t;
  for (std::string w; ss >> w;) t.push_back(w);
  return t;
}
double to_d(const std::string& s) { try { return std::stod(s); } catch (...) { return std::nan(""); } }
std::string fmt(const char* f, double v) { char b[96]; std::snprintf(b, sizeof b, f, v); return b; }
double age_seconds(const std::string& p) {
  std::error_code ec;
  const auto t = fs::last_write_time(p, ec);
  if (ec) return 0;
  return std::chrono::duration<double>(fs::file_time_type::clock::now() - t).count();
}
bool exists(const std::string& p) { std::error_code ec; return fs::exists(p, ec); }
std::string join(const std::string& a, const std::string& b) { return (fs::path(a) / b).string(); }

// in-plane minimum image (a, b periodic)
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

}  // namespace

// ---------------------------------------------------------------- readers
Outcar read_outcar(const std::string& path) {
  Outcar o;
  const std::string t = slurp(path);
  if (t.empty()) return o;
  o.exists = true;
  o.finished = t.find("General timing and accounting") != std::string::npos;
  o.relaxed = t.find("reached required accuracy") != std::string::npos;
  o.ediff_reached = t.find("EDIFF is reached") != std::string::npos;
  o.scf_failed = t.find("was not achieved") != std::string::npos;
  const auto L = lines_of(t);
  std::set<std::string> warn;
  for (size_t k = 0; k < L.size(); ++k) {
    const std::string& l = L[k];
    size_t p;
    if ((p = l.find("energy(sigma->0) =")) != std::string::npos) o.e0.push_back(to_d(l.substr(p + 18)));
    else if (l.find("TOTAL-FORCE (eV/Angst)") != std::string::npos && o.nions > 0) {
      std::vector<Vec3> f;
      double fm = 0;
      for (int a = 0; a < o.nions && k + 2 + size_t(a) < L.size(); ++a) {
        const auto w = split(L[k + 2 + size_t(a)]);
        if (w.size() < 6) break;
        const Vec3 v{to_d(w[3]), to_d(w[4]), to_d(w[5])};
        f.push_back(v);
        fm = std::max(fm, norm(v));
      }
      if (int(f.size()) == o.nions) o.forces.push_back(f), o.fmax.push_back(fm);
    } else if (l.find("  in kB") == 0 || l.find("in kB ") != std::string::npos) {
      const auto w = split(l);
      if (w.size() >= 8 && w[0] == "in" && w[1] == "kB") o.stress.push_back({to_d(w[2]), to_d(w[3]), to_d(w[4]), to_d(w[5]), to_d(w[6]), to_d(w[7])});
    } else if ((p = l.find("magnetization")) != std::string::npos && l.find("number of electron") != std::string::npos) {
      const auto w = split(l.substr(p + 13));
      if (!w.empty()) o.magnetization.push_back(to_d(w[0]));
    } else if ((p = l.find("E-fermi :")) != std::string::npos) {
      const auto w = split(l.substr(p + 9));
      if (!w.empty()) o.efermi.push_back(to_d(w[0]));
    } else if ((p = l.find("dipolmoment")) != std::string::npos) {
      const auto w = split(l.substr(p + 11));
      if (w.size() >= 3) o.dipole.push_back(to_d(w[2]));
    } else if (l.find("WARNING") != std::string::npos) {
      std::string s = l;
      s.erase(0, s.find_first_not_of(" |"));
      if (warn.size() < 20) warn.insert(s.substr(0, 120));
    } else if ((p = l.find("NIONS =")) != std::string::npos) o.nions = int(to_d(split(l.substr(p + 7)).at(0)));
    else if ((p = l.find("NELM   =")) != std::string::npos) { const auto w = split(l.substr(p + 8)); if (!w.empty()) o.nelm = int(to_d(w[0].substr(0, w[0].find(';')))); }
    else if ((p = l.find("NSW    =")) != std::string::npos) { const auto w = split(l.substr(p + 8)); if (!w.empty()) o.nsw = int(to_d(w[0])); }
    else if ((p = l.find("IBRION =")) != std::string::npos) { const auto w = split(l.substr(p + 8)); if (!w.empty()) o.ibrion = int(to_d(w[0])); }
    else if ((p = l.find("ENCUT  =")) != std::string::npos && o.encut == 0) { const auto w = split(l.substr(p + 8)); if (!w.empty()) o.encut = to_d(w[0]); }
    else if ((p = l.find("Elapsed time (sec):")) != std::string::npos) o.seconds = to_d(l.substr(p + 19));
    else if (l.find(" f  =") != std::string::npos && l.find("cm-1") != std::string::npos) {
      const auto w = split(l.substr(0, l.find("cm-1")));
      if (!w.empty()) o.freq_real.push_back(to_d(w.back()));
    } else if (l.find(" f/i=") != std::string::npos && l.find("cm-1") != std::string::npos) {
      const auto w = split(l.substr(0, l.find("cm-1")));
      if (!w.empty()) o.freq_imag.push_back(to_d(w.back()));
    }
  }
  o.warnings.assign(warn.begin(), warn.end());
  return o;
}

std::vector<IonicStep> read_oszicar(const std::string& path) {
  std::vector<IonicStep> out;
  int scf = 0;
  static const std::regex md(R"(T=\s*([-+\d.Ee]+)\s+E=\s*([-+\d.Ee]+)\s+F=\s*([-+\d.Ee]+)\s+E0=\s*([-+\d.Ee]+))");
  for (const auto& l : lines_of(slurp(path))) {
    if (l.rfind("DAV", 0) == 0 || l.rfind("RMM", 0) == 0) { ++scf; continue; }
    const size_t pf = l.find(" F= ");
    if (pf == std::string::npos) continue;
    IonicStep s;
    std::smatch m;
    if (std::regex_search(l, m, md)) s.temperature = to_d(m[1]), s.etot = to_d(m[2]), s.f = to_d(m[3]), s.e0 = to_d(m[4]);
    else {
      s.f = to_d(split(l.substr(pf + 4)).at(0));
      const size_t pe = l.find("E0=");
      s.e0 = pe == std::string::npos ? s.f : to_d(split(l.substr(pe + 3)).at(0));
    }
    s.scf = scf;
    scf = 0;
    out.push_back(s);
  }
  return out;
}

std::vector<IonicStep> read_ionic_steps(const std::string& dir, std::string* used) {
  std::vector<IonicStep> best;
  std::string b;
  for (const char* f : {"OSZICAR", "vasp.out"}) {
    const std::string p = join(dir, f);
    if (!exists(p)) continue;
    auto s = read_oszicar(p);
    if (s.size() > best.size() || b.empty()) best = s, b = p;
  }
  if (used) *used = b;
  return best;
}

std::vector<System> read_xdatcar(const std::string& path) {
  const auto L = lines_of(slurp(path));
  if (L.size() < 8) throw std::runtime_error(path + ": not an XDATCAR");
  System base;
  const double sc = to_d(split(L[1]).at(0));
  Vec3 v[3];
  for (int r = 0; r < 3; ++r) { const auto w = split(L[size_t(2 + r)]); v[r] = Vec3{to_d(w[0]), to_d(w[1]), to_d(w[2])} * sc; }
  base.cell.a = v[0], base.cell.b = v[1], base.cell.c = v[2];
  const auto names = split(L[5]);
  std::vector<int> counts;
  for (const auto& w : split(L[6])) counts.push_back(int(to_d(w)));
  std::vector<int> el;
  for (size_t s = 0; s < names.size() && s < counts.size(); ++s) for (int i = 0; i < counts[s]; ++i) el.push_back(element_from_symbol(names[s].substr(0, names[s].find('_'))));
  std::vector<System> frames;
  for (size_t k = 7; k < L.size();) {
    if (L[k].find("configuration") == std::string::npos) { ++k; continue; }
    System f = base;
    for (size_t i = 0; i < el.size() && k + 1 + i < L.size(); ++i) {
      const auto w = split(L[k + 1 + i]);
      Atom a;
      a.element = el[i];
      a.id = int64_t(i + 1);
      a.pos = base.cell.to_cartesian({to_d(w[0]), to_d(w[1]), to_d(w[2])});
      f.atoms.push_back(a);
    }
    frames.push_back(std::move(f));
    k += 1 + el.size();
  }
  return frames;
}

Doscar read_doscar(const std::string& path) {
  const auto L = lines_of(slurp(path));
  if (L.size() < 7) throw std::runtime_error(path + ": not a DOSCAR");
  Doscar d;
  const int nat = int(to_d(split(L[0]).at(0)));
  const auto h = split(L[5]);
  const int nedos = int(to_d(h.at(2)));
  d.efermi = to_d(h.at(3));
  size_t k = 6;
  for (int e = 0; e < nedos && k < L.size(); ++e, ++k) {
    const auto w = split(L[k]);
    d.energy.push_back(to_d(w[0]) - d.efermi);
    d.total.push_back(w.size() >= 5 ? to_d(w[1]) + to_d(w[2]) : to_d(w[1]));
  }
  for (int a = 0; a < nat && k < L.size(); ++a) {
    ++k;   // the atom's header
    std::vector<double> s;
    for (int e = 0; e < nedos && k < L.size(); ++e, ++k) {
      const auto w = split(L[k]);
      double sum = 0;
      for (size_t q = 1; q < w.size(); ++q) sum += to_d(w[q]);
      s.push_back(sum);
    }
    d.site.push_back(s);
  }
  return d;
}

VaspGrid read_vasp_grid(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("cannot open " + path);
  VaspGrid g;
  std::string l;
  std::vector<std::string> head;
  for (int k = 0; k < 7 && std::getline(in, l); ++k) head.push_back(l);
  if (head.size() < 7) throw std::runtime_error(path + ": not a VASP grid file");
  const double sc = to_d(split(head[1]).at(0));
  Vec3 v[3];
  for (int r = 0; r < 3; ++r) { const auto w = split(head[size_t(2 + r)]); v[r] = Vec3{to_d(w[0]), to_d(w[1]), to_d(w[2])} * sc; }
  g.system.cell.a = v[0], g.system.cell.b = v[1], g.system.cell.c = v[2];
  const auto names = split(head[5]);
  std::vector<int> counts;
  for (const auto& w : split(head[6])) counts.push_back(int(to_d(w)));
  std::getline(in, l);   // Direct
  for (size_t s = 0; s < names.size(); ++s)
    for (int i = 0; i < counts[s]; ++i) {
      std::getline(in, l);
      const auto w = split(l);
      Atom a;
      a.element = element_from_symbol(names[s].substr(0, names[s].find('_')));
      a.id = int64_t(g.system.atoms.size() + 1);
      a.pos = g.system.cell.to_cartesian({to_d(w[0]), to_d(w[1]), to_d(w[2])});
      g.system.atoms.push_back(a);
    }
  while (std::getline(in, l) && split(l).size() < 3) {}
  const auto nw = split(l);
  g.n = {int(to_d(nw[0])), int(to_d(nw[1])), int(to_d(nw[2]))};
  const size_t total = size_t(g.n[0]) * size_t(g.n[1]) * size_t(g.n[2]);
  g.values.reserve(total);
  double x;
  while (g.values.size() < total && in >> x) g.values.push_back(x);
  if (g.values.size() != total) throw std::runtime_error(path + ": the grid ends early");
  return g;
}

void write_chgcar(const VaspGrid& g, const std::string& path) {
  // atoms grouped as given (species contiguous), Direct coordinates, five values a line
  std::ofstream f(path);
  if (!f) throw std::runtime_error("cannot write " + path);
  char b[200];
  f << "CAPS grid\n1.0\n";
  for (const Vec3& v : {g.system.cell.a, g.system.cell.b, g.system.cell.c}) { std::snprintf(b, sizeof b, " %12.6f %12.6f %12.6f\n", v[0], v[1], v[2]); f << b; }
  std::vector<std::pair<int, int>> runs;
  for (const auto& a : g.system.atoms) { if (!runs.empty() && runs.back().first == a.element) ++runs.back().second; else runs.push_back({a.element, 1}); }
  for (const auto& r : runs) f << "   " << element(r.first).symbol;
  f << "\n";
  for (const auto& r : runs) f << "   " << r.second;
  f << "\nDirect\n";
  for (const auto& a : g.system.atoms) { const Vec3 fr = g.system.cell.to_fractional(a.pos); std::snprintf(b, sizeof b, " %10.6f %10.6f %10.6f\n", fr[0], fr[1], fr[2]); f << b; }
  f << "\n" << g.n[0] << " " << g.n[1] << " " << g.n[2] << "\n";
  for (size_t i = 0; i < g.values.size(); ++i) {
    std::snprintf(b, sizeof b, " %17.10E", g.values[i]);
    f << b << ((i + 1) % 5 == 0 ? "\n" : "");
  }
  f << "\n";
}

std::vector<double> read_acf(const std::string& path) {
  std::vector<double> q;
  for (const auto& l : lines_of(slurp(path))) {
    const auto w = split(l);
    if (w.size() >= 5 && std::all_of(w[0].begin(), w[0].end(), ::isdigit)) q.push_back(to_d(w[4]));
  }
  return q;
}

std::vector<std::pair<std::string, int>> poscar_species(const std::string& path) {
  const auto L = lines_of(slurp(path));
  if (L.size() < 7) throw std::runtime_error(path + ": not a POSCAR");
  const auto n = split(L[5]), c = split(L[6]);
  std::vector<std::pair<std::string, int>> out;
  for (size_t k = 0; k < n.size() && k < c.size(); ++k) out.push_back({n[k], int(to_d(c[k]))});
  return out;
}

int steps_to_target(const std::vector<double>& fm, double target) {
  if (fm.size() < 4 || fm.back() <= target) return fm.empty() || fm.back() > target ? -1 : 0;
  const size_t k = std::min<size_t>(fm.size(), 8);
  double mx = 0, my = 0;
  std::vector<double> x, y;
  for (size_t i = fm.size() - k; i < fm.size(); ++i) x.push_back(double(i)), y.push_back(std::log(std::max(fm[i], 1e-6)));
  for (size_t i = 0; i < k; ++i) mx += x[i] / double(k), my += y[i] / double(k);
  double sxy = 0, sxx = 0;
  for (size_t i = 0; i < k; ++i) sxy += (x[i] - mx) * (y[i] - my), sxx += (x[i] - mx) * (x[i] - mx);
  const double slope = sxy / std::max(sxx, 1e-9);
  if (slope >= -1e-3) return -1;
  return std::max(1, int((std::log(target) - y.back()) / slope));
}

// ---------------------------------------------------------------- checks
namespace {

const std::vector<std::string> kStages = {"01_cell", "02_cell2", "03_relax", "04_static"};

double lattice_a(const std::string& poscar) {
  const auto L = lines_of(slurp(poscar));
  if (L.size() < 3) return std::nan("");
  const double s = to_d(split(L[1]).at(0));
  const auto w = split(L[2]);
  return s * std::sqrt(to_d(w[0]) * to_d(w[0]) + to_d(w[1]) * to_d(w[1]) + to_d(w[2]) * to_d(w[2]));
}

Json check_stage(const std::string& d, bool complex, const std::string& data_dir) {
  Json r = Json::object();
  r["folder"] = d;
  const Outcar o = read_outcar(join(d, "OUTCAR"));
  if (!o.exists) { r["status"] = "no OUTCAR yet"; r["ok"] = false; return r; }
  const std::string base = fs::path(d).filename().string();
  const bool is_static = base == "04_static" || o.nsw == 0;
  const auto steps = read_ionic_steps(d);
  r["finished"] = o.finished;
  r["converged"] = is_static ? (o.ediff_reached ? "single point" : "single point, EDIFF not reached") : (o.relaxed ? "yes" : "no");
  r["scf_failed"] = o.scf_failed;
  r["ionic_steps"] = int(steps.size());
  if (!o.e0.empty()) r["E0"] = o.e0.back();
  if (!o.fmax.empty()) r["max_force"] = o.fmax.back();
  if (!o.magnetization.empty()) r["magnetization"] = o.magnetization.back();
  Json w = Json::array();
  for (size_t k = 0; k < o.warnings.size() && k < 3; ++k) w.push_back(o.warnings[k]);
  r["warnings"] = w;
  bool ok = o.finished && !o.scf_failed && (o.relaxed || is_static);
  const std::string pin = join(d, "POSCAR.in"), cont = join(d, "CONTCAR");
  if (base.find("cell") != std::string::npos && exists(pin) && exists(cont)) {
    const double a0 = lattice_a(pin), a1 = lattice_a(cont);
    r["a_before"] = a0, r["a_after"] = a1;
    if (base == "02_cell2" && std::fabs(a1 - a0) > 0.002) {
      r["note"] = "a still changes in the restart: run 02_cell2 once more (Pulay stress)";
      ok = false;
    }
  }
  if (!complex && exists(cont) && fs::file_size(cont) > 0) {
    try {
      std::map<std::string, std::string> bad;
      const System s = read_vasp_poscar(cont, &bad);
      ValidateOptions vo;
      vo.bad_labels = bad;
      vo.data_dir = data_dir;
      const auto v = validate_2d(s, vo);
      r["structure"] = v.status;
      Json e = Json::array();
      for (const auto& f : v.findings) if (f.level == "ERROR") e.push_back(f.text);
      r["structure_errors"] = e;
      ok = ok && v.status == "PASS";
    } catch (const std::exception& e) { r["structure"] = std::string("skipped (") + e.what() + ")"; }
  }
  r["ok"] = ok;
  r["status"] = ok ? "OK" : "LOOK AT THIS RUN";
  return r;
}

struct LogState { std::vector<std::string> done, started; std::string ws; bool not_converged = false; std::string file; };
LogState read_job_log(const std::string& case_dir) {
  LogState s;
  std::string newest;
  fs::file_time_type nt{};
  std::error_code ec;
  for (const auto& e : fs::directory_iterator(case_dir, ec)) {
    const std::string n = e.path().filename().string();
    if (n.rfind("slurm-", 0) == 0 && n.size() > 10 && n.substr(n.size() - 4) == ".out") {
      const auto t = fs::last_write_time(e.path(), ec);
      if (newest.empty() || t > nt) newest = e.path().string(), nt = t;
    }
  }
  if (newest.empty()) return s;
  s.file = newest;
  static const std::regex done(R"(^== (\S+) done)"), start(R"(^== (\S+)\s+start)"), ws(R"(^workspace: (\S+))");
  for (const auto& l : lines_of(slurp(newest))) {
    std::smatch m;
    if (std::regex_search(l, m, done)) s.done.push_back(m[1]);
    else if (std::regex_search(l, m, start)) s.started.push_back(m[1]);
    else if (std::regex_search(l, m, ws)) s.ws = m[1];
    if (l.find("NOT converged") != std::string::npos) s.not_converged = true;
  }
  return s;
}

std::vector<std::string> job_stages(const std::string& case_dir) {
  static const std::regex st("^STAGES=\"(.*)\"");
  for (const auto& l : lines_of(slurp(join(case_dir, "job.slurm")))) {
    std::smatch m;
    if (std::regex_search(l, m, st)) return split(m[1]);
  }
  return kStages;
}

double last_e0(const std::string& d) {
  const auto o = read_outcar(join(d, "OUTCAR"));
  return o.e0.empty() ? std::nan("") : o.e0.back();
}

// molecule atoms = those beyond the slab's composition per element (slab atoms first in each element block)
std::pair<std::vector<size_t>, std::vector<size_t>> split_molecule(const System& at, const std::string& slab_poscar) {
  std::map<std::string, int> sc;
  for (const auto& [n, c] : poscar_species(slab_poscar)) sc[n] += c;
  std::map<std::string, int> seen;
  std::vector<size_t> mol, slab;
  for (size_t i = 0; i < at.atoms.size(); ++i) {
    const std::string e = element(at.atoms[i].element).symbol;
    (++seen[e] > sc[e] ? mol : slab).push_back(i);
  }
  return {mol, slab};
}

double top_plane(const System& at, const std::vector<size_t>& slab) {
  double zx = -1e300;
  for (size_t i : slab) if (at.atoms[i].element != 1) zx = std::max(zx, at.atoms[i].pos[2]);
  double sum = 0;
  int n = 0;
  for (size_t i : slab) if (at.atoms[i].element != 1 && at.atoms[i].pos[2] > zx - 1.0) sum += at.atoms[i].pos[2], ++n;
  return n ? sum / n : zx;
}

std::set<std::pair<size_t, size_t>> mol_bonds(const System& at, const std::vector<size_t>& mol) {
  std::set<std::pair<size_t, size_t>> b;
  for (size_t a = 0; a < mol.size(); ++a)
    for (size_t c = a + 1; c < mol.size(); ++c) {
      const double d = norm(mic_ab(at.cell, at.atoms[mol[c]].pos - at.atoms[mol[a]].pos));
      if (d < 1.2 * (element(at.atoms[mol[a]].element).covalent + element(at.atoms[mol[c]].element).covalent)) b.insert({a, c});
    }
  return b;
}

}  // namespace

Json check_runs(const std::vector<std::string>& dirs, bool complex, const std::string& data_dir) {
  Json out = Json::array();
  for (const auto& d : dirs) {
    bool any = false;
    for (const auto& s : kStages)
      if (fs::is_directory(join(d, s))) { out.push_back(check_stage(join(d, s), complex, data_dir)); any = true; }
    if (!any) out.push_back(check_stage(d, complex, data_dir));
  }
  return out;
}

Json run_progress(const std::string& case_dir) {
  Json r = Json::object();
  r["case"] = case_dir;
  const auto lg = read_job_log(case_dir);
  const auto stages = job_stages(case_dir);
  Json sj = Json::array();
  for (const auto& s : stages) sj.push_back(s);
  r["stages"] = sj;
  if (lg.file.empty()) { r["state"] = "not started (no slurm-*.out)"; return r; }
  r["log"] = fs::path(lg.file).filename().string();
  Json dj = Json::array();
  for (const auto& s : lg.done) dj.push_back(s);
  r["done"] = dj;
  if (lg.not_converged) r["note"] = "a stage stopped NOT converged: submit job.slurm again (it resumes)";
  std::string running;
  for (const auto& s : lg.started) if (std::find(lg.done.begin(), lg.done.end(), s) == lg.done.end()) running = s;
  if (running.empty()) { r["state"] = lg.done.size() == stages.size() ? "all stages done" : "nothing running"; return r; }
  r["state"] = "running";
  r["running"] = running;
  const std::string d = !lg.ws.empty() && fs::is_directory(join(lg.ws, running)) ? join(lg.ws, running) : join(case_dir, running);
  r["live"] = d;
  const auto steps = read_ionic_steps(d);
  r["ionic_steps"] = int(steps.size());
  Json last = Json::array();
  for (size_t k = steps.size() > 3 ? steps.size() - 3 : 0; k < steps.size(); ++k) last.push_back(steps[k].e0);
  r["last_E0"] = last;
  const auto o = read_outcar(join(d, "OUTCAR"));
  if (!o.fmax.empty()) r["max_force"] = o.fmax.back(), r["first_force"] = o.fmax.front();
  Json fj = Json::array();
  for (double f : o.fmax) fj.push_back(f);
  r["fmax"] = fj;
  const std::string start = exists(join(d, "POSCAR")) ? join(d, "POSCAR") : "";
  if (!start.empty() && !steps.empty()) {
    const double per = age_seconds(start) / (double(steps.size()) + 0.5);
    r["seconds_per_step"] = per;
    if (running != "04_static") {
      const int left = steps_to_target(o.fmax);
      r["steps_left"] = left;
      if (left > 0) r["eta_seconds"] = left * per;
    }
  }
  if (exists(join(d, "vasp.out"))) r["seconds_since_output"] = age_seconds(join(d, "vasp.out"));
  return r;
}

Json health_check(const std::string& root) {
  Json r = Json::object();
  auto ref = [&](const std::string& c, std::string* st) {
    for (const char* s : {"04_static", "03_relax"}) { const double e = last_e0(join(c, s)); if (std::isfinite(e)) { *st = s; return e; } }
    return std::nan("");
  };
  std::string ss, sm;
  const double es = ref(join(root, "slab"), &ss), em = ref(join(root, "molecule"), &sm);
  if (!std::isfinite(es) || !std::isfinite(em)) { r["error"] = "slab or molecule has no energy yet: run the health check after they finished"; return r; }
  const double eref = es + em;
  r["E_ref"] = eref, r["E_slab"] = es, r["E_slab_stage"] = ss, r["E_molecule"] = em, r["E_molecule_stage"] = sm;
  const std::string slab_poscar = join(join(root, "slab"), "POSCAR");
  const System mol0 = read_vasp_poscar(join(join(root, "molecule"), "POSCAR"));
  std::vector<size_t> all(mol0.atoms.size());
  for (size_t i = 0; i < all.size(); ++i) all[i] = i;
  const auto ref_bonds = mol_bonds(mol0, all);
  Json rows = Json::array();
  int bad = 0;
  std::vector<std::string> cases;
  std::error_code ec;
  for (const auto& e : fs::directory_iterator(root, ec)) if (e.is_directory() && e.path().filename().string().rfind("complex_", 0) == 0) cases.push_back(e.path().string());
  std::sort(cases.begin(), cases.end());
  for (const auto& c : cases) {
    Json row = Json::object();
    row["complex"] = fs::path(c).filename().string();
    Json flags = Json::array();
    auto flag = [&](const std::string& kind, const std::string& text, const std::string& action) {
      Json f = Json::object();
      f["kind"] = kind, f["text"] = text, f["action"] = action;
      flags.push_back(f);
    };
    // the stage running or that ran last
    std::string d;
    bool running = false;
    const auto lg = read_job_log(c);
    for (const auto& s : lg.started)
      if (std::find(lg.done.begin(), lg.done.end(), s) == lg.done.end() && !lg.ws.empty() && fs::is_directory(join(lg.ws, s))) d = join(lg.ws, s), running = true;
    if (d.empty()) for (auto it = kStages.rbegin(); it != kStages.rend(); ++it) if (fs::is_directory(join(c, *it))) { d = join(c, *it); break; }
    if (d.empty()) { row["state"] = "not started"; row["flags"] = flags; rows.push_back(row); continue; }
    row["stage"] = fs::path(d).filename().string();
    row["running"] = running;
    const auto steps = read_ionic_steps(d);
    row["steps"] = int(steps.size());
    const auto o = read_outcar(join(d, "OUTCAR"));
    const int nelm = o.nelm > 0 ? o.nelm : 120;
    if (!steps.empty()) {
      const double de = steps.back().e0 - eref;
      row["dE"] = de;
      if (std::fabs(de) > 3.0) flag("ENERGY", fmt("%+.2f", de) + " eV from slab+molecule", "a physisorbed molecule binds by at most ~2–3 eV: check the geometry and the dipole settings; reset the stage if it ran away");
      double jump = 0;
      for (size_t k = 1; k < steps.size(); ++k) jump = std::max(jump, std::fabs(steps[k].e0 - steps[k - 1].e0));
      if (jump > 2.0) flag("JUMP", fmt("%.2f", jump) + " eV between steps", "electronic or geometric runaway: reset the stage and restart from the last sane geometry");
      int hit = 0;
      for (const auto& s : steps) hit += s.scf >= nelm;
      if (hit) flag("SCF", "hit NELM=" + std::to_string(nelm) + " in " + std::to_string(hit) + " step(s)", "the electronic loop did not converge: ALGO = All or a smaller mixing, then restart");
    }
    if (!o.dipole.empty()) {
      row["dipole"] = o.dipole.back();
      if (std::fabs(o.dipole.back()) > 1.0) flag("DIPOLE", fmt("%+.2f", o.dipole.back()) + " e*A", "a neutral molecule on a metal cannot polarise that much: check IDIPOL/DIPOL (no LDIPOL) and the geometry");
    }
    try {
      const std::string cont = join(d, "CONTCAR");
      const System at = read_vasp_poscar(exists(cont) && fs::file_size(cont) > 0 ? cont : join(c, "POSCAR"));
      const auto [mol, slab] = split_molecule(at, slab_poscar);
      double dmin = 1e300;
      for (size_t i : mol)
        for (size_t j = 0; j < at.atoms.size(); ++j)
          if (j != i) dmin = std::min(dmin, norm(mic_ab(at.cell, at.atoms[j].pos - at.atoms[i].pos)));
      if (dmin < 0.8) flag("GEOM", "atoms " + fmt("%.2f", dmin) + " A apart", "atoms collapsed: reset the stage");
      const auto now = mol_bonds(at, mol);
      int lost = 0;
      for (const auto& b : ref_bonds) lost += !now.count(b);
      if (lost) flag("GEOM", std::to_string(lost) + " molecule bond(s) broken", "the molecule reacted or fell apart: look at it in the 3D view before using any number");
      const double plane = top_plane(at, slab);
      double hmin = 1e300;
      for (size_t i : mol) if (at.atoms[i].element != 1) hmin = std::min(hmin, at.atoms[i].pos[2] - plane);
      row["h_min"] = hmin;
      if (hmin > 6.0) flag("GEOM", "molecule desorbed (h_min " + fmt("%.1f", hmin) + " A)", "no binding at this site: the complex is not a minimum; keep it out of E_bind");
      if (hmin < 0.0) flag("GEOM", "molecule below the termination plane (h_min " + fmt("%.1f", hmin) + " A)", "the molecule entered the slab: reset the stage");
    } catch (const std::exception& e) { flag("GEOM", std::string("not checked (") + e.what() + ")", ""); }
    if (running && exists(join(d, "vasp.out"))) {
      const double age = age_seconds(join(d, "vasp.out"));
      if (age > 1200) flag("HANG", "no output for " + fmt("%.0f", age / 60) + " min", "VASP hung: the job's watchdog restarts the stage; else cancel and resubmit (it resumes)");
    }
    row["flags"] = flags;
    bad += flags.size() > 0;
    rows.push_back(row);
  }
  r["complexes"] = rows;
  r["flagged"] = bad;
  return r;
}

}  // namespace caps
