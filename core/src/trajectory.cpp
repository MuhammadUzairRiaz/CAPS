// CAPS trajectory player data: LAMMPS logs, per-frame series, chain ends, smoothing (see caps/trajectory.hpp).
#include "caps/trajectory.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>

#include "caps/analysis.hpp"

namespace caps {

namespace {

bool is_number(const std::string& t, double& v) {
  if (t.empty()) return false;
  char* end = nullptr;
  v = std::strtod(t.c_str(), &end);
  return end && *end == '\0';
}

std::vector<std::string> words(const std::string& line) {
  std::istringstream in(line);
  std::vector<std::string> w;
  for (std::string x; in >> x;) w.push_back(x);
  return w;
}

}  // namespace

ThermoLog parse_lammps_log(const std::string& text) {
  ThermoLog L;
  std::istringstream in(text);
  std::string line;
  std::vector<std::string> header;   // the current one-line table's columns
  bool in_table = false;
  std::map<std::string, size_t> col_of;
  auto column = [&](const std::string& name) {
    auto it = col_of.find(name);
    if (it != col_of.end()) return it->second;
    col_of[name] = L.columns.size();
    L.columns.push_back(name);
    for (auto& r : L.rows) r.push_back(std::numeric_limits<double>::quiet_NaN());
    return L.columns.size() - 1;
  };
  column("Step");
  auto new_row = [&] { L.rows.emplace_back(L.columns.size(), std::numeric_limits<double>::quiet_NaN()); return L.rows.size() - 1; };
  bool multi_run_open = false;
  while (std::getline(in, line)) {
    if (!line.empty() && line.back() == '\r') line.pop_back();
    const auto w = words(line);
    if (w.empty()) continue;
    // multi style: "---------------- Step N ----- CPU = … ----------"
    if (w.size() >= 3 && w[0].rfind("----", 0) == 0 && w[1] == "Step") {
      double step;
      if (!is_number(w[2], step)) continue;
      if (!multi_run_open) { L.run_starts.push_back(L.rows.size()); multi_run_open = true; }
      const size_t r = new_row();
      L.rows[r][0] = step;
      // "Name = value" pairs on the following lines, until a blank line or the next block
      std::streampos mark = in.tellg();
      while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        const auto v = words(line);
        if (v.empty() || v[0].rfind("----", 0) == 0 || v[0] == "Loop") { in.seekg(mark); break; }
        bool any = false;
        for (size_t k = 0; k + 2 < v.size(); k += 3) {
          double x;
          if (v[k + 1] != "=" || !is_number(v[k + 2], x)) break;
          const size_t c = column(v[k]);
          L.rows[r][c] = x;
          any = true;
        }
        if (!any) { in.seekg(mark); break; }
        mark = in.tellg();
      }
      continue;
    }
    if (w[0] == "Loop" && w.size() > 1 && w[1] == "time") { in_table = false; multi_run_open = false; continue; }
    // one-line style: a header starting with Step, then numeric rows of the same width
    if (w[0] == "Step" && w.size() > 1) {
      bool numeric_header = false;
      double x;
      for (const auto& h : w) numeric_header |= is_number(h, x);
      if (numeric_header) continue;
      header = w;
      in_table = true;
      L.run_starts.push_back(L.rows.size());
      for (const auto& h : header) column(h);
      continue;
    }
    if (in_table) {
      if (w.size() != header.size()) { in_table = false; continue; }
      std::vector<double> vals(w.size());
      bool ok = true;
      for (size_t k = 0; k < w.size() && ok; ++k) ok = is_number(w[k], vals[k]);
      if (!ok) { in_table = false; continue; }
      const size_t r = new_row();
      for (size_t k = 0; k < header.size(); ++k) L.rows[r][col_of[header[k]]] = vals[k];
    }
  }
  if (L.rows.empty()) throw std::runtime_error("no thermo output in the log");
  return L;
}

ThermoLog read_lammps_log(const std::string& path) {
  std::ifstream f(path);
  if (!f) throw std::runtime_error("cannot read " + path);
  std::stringstream ss;
  ss << f.rdbuf();
  return parse_lammps_log(ss.str());
}

ChainEnds chain_ends(const System& s, int molecule) {
  ChainEnds e;
  const auto mol = s.molecules();
  if (s.atoms.empty()) return e;
  // molecule ids as the Studio shows them (1-based over molecules())
  std::map<int, int> size;
  for (int m : mol) ++size[m + 1];
  if (molecule <= 0 || !size.count(molecule)) {
    int best = 0;
    for (const auto& [m, n] : size) if (n > best) best = n, molecule = m;
  }
  e.molecule = molecule;
  for (const auto& bb : backbones(s, 2)) {
    if (bb.empty() || mol[bb.front()] + 1 != molecule) continue;
    e.first = int(bb.front()), e.last = int(bb.back());
    return e;
  }
  // no backbone (a small or branched molecule): its two farthest atoms
  double best = -1;
  std::vector<size_t> atoms;
  for (size_t i = 0; i < s.atoms.size(); ++i) if (mol[i] + 1 == molecule) atoms.push_back(i);
  for (size_t a = 0; a < atoms.size() && atoms.size() < 3000; ++a)
    for (size_t b = a + 1; b < atoms.size(); ++b) {
      const double d = norm(s.atoms[atoms[a]].pos - s.atoms[atoms[b]].pos);
      if (d > best) best = d, e.first = int(atoms[a]), e.last = int(atoms[b]);
    }
  return e;
}

DataTable trajectory_series(const Trajectory& t, const TrajectorySeriesOptions& o) {
  DataTable T;
  T.name = "trajectory";
  T.title = "Per-frame series";
  const size_t nf = t.frames();
  if (nf == 0) throw std::runtime_error("no frames");
  System first = t.frame(0);
  if (!first.unwrapped) make_molecules_whole(first);
  const ChainEnds ends = chain_ends(first, o.molecule);
  const bool cell = first.cell.valid();
  T.columns = {"Frame", "Timestep", "Time (ps)"};
  if (cell) T.columns.push_back("Density (g/cm³)"), T.columns.push_back("Volume (Å³)");
  T.columns.push_back("Rg (Å)");
  T.columns.push_back("Ree (Å)");
  std::vector<size_t> log_cols;
  if (o.log)
    for (size_t c = 1; c < o.log->columns.size(); ++c) {
      T.columns.push_back(o.log->columns[c]);
      log_cols.push_back(c);
    }
  // the log's steps, sorted with their rows, for the nearest-step lookup
  std::vector<std::pair<double, size_t>> steps;
  double interval = 0;
  if (o.log) {
    for (size_t r = 0; r < o.log->rows.size(); ++r) if (!std::isnan(o.log->rows[r][0])) steps.emplace_back(o.log->rows[r][0], r);
    std::sort(steps.begin(), steps.end());
    for (size_t k = 1; k < steps.size(); ++k) if (steps[k].first > steps[k - 1].first) { interval = steps[k].first - steps[k - 1].first; break; }
  }
  const int stride = std::max(1, o.stride);
  const int total = int((nf + size_t(stride) - 1) / size_t(stride));
  int done = 0;
  for (size_t k = 0; k < nf; k += size_t(stride)) {
    System f = t.frame(k);
    if (!f.unwrapped) make_molecules_whole(f);
    std::vector<double> row;
    const double step = k < t.timesteps.size() ? double(t.timesteps[k]) : double(k);
    row.push_back(double(k));
    row.push_back(step);
    row.push_back(step * o.dt_fs * 1e-3);
    if (cell) row.push_back(f.density()), row.push_back(f.cell.volume());
    // Rg of the chain (mass-weighted), Ree of its backbone ends
    double m = 0;
    Vec3 com{0, 0, 0};
    const auto mol = f.molecules();
    for (size_t i = 0; i < f.atoms.size(); ++i)
      if (mol[i] + 1 == ends.molecule) { const double w = f.mass_of(f.atoms[i]); com = com + f.atoms[i].pos * w; m += w; }
    double rg = std::numeric_limits<double>::quiet_NaN();
    if (m > 0) {
      com = com * (1.0 / m);
      double s2 = 0;
      for (size_t i = 0; i < f.atoms.size(); ++i)
        if (mol[i] + 1 == ends.molecule) { const Vec3 d = f.atoms[i].pos - com; s2 += f.mass_of(f.atoms[i]) * dot(d, d); }
      rg = std::sqrt(s2 / m);
    }
    row.push_back(rg);
    row.push_back(ends.first >= 0 && size_t(ends.last) < f.atoms.size() ? norm(f.atoms[size_t(ends.last)].pos - f.atoms[size_t(ends.first)].pos)
                                                                          : std::numeric_limits<double>::quiet_NaN());
    // the log at this step
    if (!log_cols.empty()) {
      const auto it = std::lower_bound(steps.begin(), steps.end(), std::make_pair(step, size_t(0)));
      size_t best = SIZE_MAX;
      double gap = 1e300;
      for (auto jt : {it, it == steps.begin() ? it : std::prev(it)})
        if (jt != steps.end() && std::fabs(jt->first - step) < gap) gap = std::fabs(jt->first - step), best = jt->second;
      const bool near = best != SIZE_MAX && (gap == 0 || (interval > 0 && gap <= interval / 2));
      for (size_t c : log_cols) row.push_back(near ? o.log->rows[best][c] : std::numeric_limits<double>::quiet_NaN());
    }
    T.rows.push_back(std::move(row));
    if (o.progress && !o.progress(++done, total)) break;
  }
  return T;
}

std::vector<Vec3> smoothed_positions(const Trajectory& t, size_t k, int window) {
  std::vector<Vec3> p = t.positions.at(k);
  if (window <= 1 || t.frames() < 2) return p;
  const int h = window / 2;
  const size_t lo = k >= size_t(h) ? k - size_t(h) : 0, hi = std::min(t.frames() - 1, k + size_t(h));
  const Cell& cell = k < t.cells.size() ? t.cells[k] : t.topology.cell;
  std::vector<Vec3> sum(p.size(), Vec3{0, 0, 0});
  int n = 0;
  for (size_t f = lo; f <= hi; ++f) {
    const auto& q = t.positions[f];
    if (q.size() != p.size()) continue;
    for (size_t i = 0; i < p.size(); ++i) {
      Vec3 d = q[i] - p[i];
      if (cell.valid()) d = cell.minimum_image(d);
      sum[i] = sum[i] + d;
    }
    ++n;
  }
  if (n > 0)
    for (size_t i = 0; i < p.size(); ++i) p[i] = p[i] + sum[i] * (1.0 / n);
  return p;
}

}  // namespace caps
