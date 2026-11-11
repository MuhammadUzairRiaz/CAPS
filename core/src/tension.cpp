// Surface and interfacial tension from the pressure tensor (tension.hpp).
#include "caps/tension.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace caps {

namespace {
constexpr double kAtmAngstromToMNm = 101325.0 * 1e-10 * 1e3;   // 1 atm·Å in mN/m
std::string lower(std::string s) {
  for (auto& c : s) c = char(std::tolower(static_cast<unsigned char>(c)));
  return s;
}
}  // namespace

TensionResult surface_tension(const std::vector<TensionSample>& s, int interfaces, double skip, int blocks) {
  if (interfaces < 1) throw std::invalid_argument("at least one interface");
  if (skip < 0 || skip >= 1) throw std::invalid_argument("skip a share between 0 and 1");
  TensionResult r;
  const size_t first = size_t(std::floor(skip * double(s.size())));
  for (size_t k = first; k < s.size(); ++k) {
    const auto& x = s[k];
    if (x.ln <= 0) throw std::invalid_argument("a sample without the cell length along the normal");
    r.series.push_back(x.ln / interfaces * (x.pn - 0.5 * (x.pt1 + x.pt2)) * kAtmAngstromToMNm);
    r.pn += x.pn, r.pt += 0.5 * (x.pt1 + x.pt2), r.ln += x.ln;
  }
  const size_t n = r.series.size();
  if (n == 0) throw std::invalid_argument("no samples after the skipped share");
  r.samples = int(n);
  r.pn /= double(n), r.pt /= double(n), r.ln /= double(n);
  for (double g : r.series) r.gamma += g;
  r.gamma /= double(n);
  r.blocks = std::max(0, std::min(blocks, int(n / 2)));
  if (r.blocks >= 2) {
    const size_t bl = n / size_t(r.blocks);
    std::vector<double> means;
    for (int b = 0; b < r.blocks; ++b) {
      double m = 0;
      for (size_t k = size_t(b) * bl; k < size_t(b + 1) * bl; ++k) m += r.series[k];
      means.push_back(m / double(bl));
    }
    double mm = 0, v = 0;
    for (double m : means) mm += m;
    mm /= double(means.size());
    for (double m : means) v += (m - mm) * (m - mm);
    r.error = std::sqrt(v / double(means.size() - 1) / double(means.size()));
  } else {
    r.blocks = 0;
    r.notes.push_back("too few samples for an error estimate");
  }
  return r;
}

std::vector<TensionSample> tension_samples(const std::vector<ThermoRow>& rows, int axis) {
  if (axis < 0 || axis > 2) throw std::invalid_argument("axis x, y or z");
  std::vector<TensionSample> out;
  for (const auto& t : rows) {
    TensionSample x;
    x.pn = t.p[axis], x.pt1 = t.p[(axis + 1) % 3], x.pt2 = t.p[(axis + 2) % 3];
    x.ln = axis == 0 ? t.lx : axis == 1 ? t.ly : t.lz;
    out.push_back(x);
  }
  return out;
}

std::vector<TensionSample> tension_samples(const ThermoLog& log, int axis, double cell_length) {
  if (axis < 0 || axis > 2) throw std::invalid_argument("axis x, y or z");
  const char* ax = "xyz";
  auto col = [&](const std::string& name) {
    for (size_t k = 0; k < log.columns.size(); ++k) if (lower(log.columns[k]) == name) return int(k);
    return -1;
  };
  const std::string a0(1, ax[axis]), a1(1, ax[(axis + 1) % 3]), a2(1, ax[(axis + 2) % 3]);
  const int pn = col("p" + a0 + a0), p1 = col("p" + a1 + a1), p2 = col("p" + a2 + a2);
  if (pn < 0 || p1 < 0 || p2 < 0) throw std::invalid_argument("the log needs the columns Pxx, Pyy and Pzz (thermo_style custom … pxx pyy pzz)");
  const int len = col("l" + a0), hi = col(a0 + "hi"), lo = col(a0 + "lo");
  if (len < 0 && (hi < 0 || lo < 0) && cell_length <= 0)
    throw std::invalid_argument("the log has no L" + a0 + " (or " + a0 + "hi and " + a0 + "lo): give the cell length");
  std::vector<TensionSample> out;
  for (const auto& row : log.rows) {
    TensionSample x;
    x.pn = row[size_t(pn)], x.pt1 = row[size_t(p1)], x.pt2 = row[size_t(p2)];
    // a run in the same log that did not print these columns leaves NaN in its rows: skip them
    if (!std::isfinite(x.pn) || !std::isfinite(x.pt1) || !std::isfinite(x.pt2)) continue;
    // the cell length: L column, else hi - lo, else the given length, taking the first one this row has
    x.ln = NAN;
    if (len >= 0) x.ln = row[size_t(len)];
    if (!std::isfinite(x.ln) && hi >= 0 && lo >= 0) x.ln = row[size_t(hi)] - row[size_t(lo)];
    if (!std::isfinite(x.ln) && cell_length > 0) x.ln = cell_length;
    if (!std::isfinite(x.ln)) continue;
    out.push_back(x);
  }
  return out;
}

std::vector<TensionSample> tension_samples_from_file(const std::string& path, int axis, double cell_length) {
  std::ifstream in(path);
  if (!in) throw std::invalid_argument("cannot read " + path);
  std::stringstream ss;
  ss << in.rdbuf();
  const std::string text = ss.str();
  const auto nl = text.find('\n');
  const std::string head = text.substr(0, nl);
  if (head.find("pxx_atm") == std::string::npos) return tension_samples(parse_lammps_log(text), axis, cell_length);
  // CAPS's thermo CSV: the header names the columns; the stage column is text
  ThermoLog log;
  {
    std::string cur;
    for (char c : head + ",") {
      if (c == ',') { std::string n = cur; if (n.size() > 4 && n.substr(n.size() - 4) == "_atm") n = n.substr(0, n.size() - 4);
                      if (n.size() > 2 && n.substr(n.size() - 2) == "_A") n = n.substr(0, n.size() - 2);
                      log.columns.push_back(n); cur.clear(); }
      else if (c != '\r') cur += c;
    }
  }
  std::istringstream rows(text.substr(nl == std::string::npos ? text.size() : nl + 1));
  for (std::string line; std::getline(rows, line);) {
    if (line.empty()) continue;
    std::vector<double> v;
    std::string cur;
    for (char c : line + ",") {
      if (c == ',') { char* e = nullptr; const double x = std::strtod(cur.c_str(), &e); v.push_back(e && e != cur.c_str() ? x : NAN); cur.clear(); }
      else if (c != '\r') cur += c;
    }
    if (v.size() >= log.columns.size()) log.rows.push_back(v);
  }
  return tension_samples(log, axis, cell_length);
}

}  // namespace caps
