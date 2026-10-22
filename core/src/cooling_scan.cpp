// Cooling scans that ran elsewhere (caps tg on a cluster) read back and fitted: one scan, or replicas pooled.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <map>
#include <regex>
#include <sstream>
#include <stdexcept>

#include "caps/json.hpp"
#include "caps/mechanics.hpp"

namespace caps {

namespace {

namespace fs = std::filesystem;
constexpr double NaN = std::numeric_limits<double>::quiet_NaN();

std::string slurp(const fs::path& p) {
  std::ifstream f(p, std::ios::binary);
  std::stringstream ss;
  ss << f.rdbuf();
  return ss.str();
}

std::string lower(std::string s) {
  for (char& c : s) c = char(std::tolower(static_cast<unsigned char>(c)));
  return s;
}

std::string g(double v, int digits) {
  char b[32];
  std::snprintf(b, sizeof b, "%.*g", digits, v);
  return b;
}

CoolingPoint from_density(double T, double d) {
  CoolingPoint p;
  p.temperature = T, p.density = d, p.specific_volume = d > 0 ? 1 / d : 0;
  return p;
}
CoolingPoint from_volume(double T, double v) {
  CoolingPoint p;
  p.temperature = T, p.specific_volume = v, p.density = v > 0 ? 1 / v : 0;
  return p;
}

// caps tg's checkpoint of the points: {"points": [{T, density, density_err, specific_volume, potential}]}
bool points_file(const Json& j, std::vector<CoolingPoint>& pts) {
  if (!j.is_object() || !j.has("points") || !j["points"].is_array()) return false;
  for (const auto& q : j["points"].items()) {
    CoolingPoint p;
    p.temperature = q.num("T", 0), p.density = q.num("density", 0), p.density_err = q.num("density_err", 0);
    p.specific_volume = q.num("specific_volume", p.density > 0 ? 1 / p.density : 0), p.potential = q.num("potential", 0);
    if (p.density <= 0 && p.specific_volume > 0) p.density = 1 / p.specific_volume;
    pts.push_back(p);
  }
  return !pts.empty();
}

// a properties file (caps tg --json, a Studio replica's out/STEM.properties.json): the tg card's density curve
bool properties_file(const Json& j, std::vector<CoolingPoint>& pts) {
  if (!j.is_array()) return false;
  for (const auto& p : j.items()) {
    if (!p.is_object() || p.text("id") != "tg" || !p.has("series")) continue;
    const Json* dens = nullptr;
    const Json* vol = nullptr;
    for (const auto& s : p["series"].items()) {
      const std::string label = s.text("label");
      if (label == "density") dens = &s;
      else if (label == "specific volume") vol = &s;
    }
    const Json* use = dens ? dens : vol;
    if (!use) continue;
    const auto& x = (*use)["x"].items();
    const auto& y = (*use)["y"].items();
    for (size_t k = 0; k < std::min(x.size(), y.size()); ++k)
      if (x[k].is_number() && y[k].is_number()) pts.push_back(dens ? from_density(x[k].number(), y[k].number()) : from_volume(x[k].number(), y[k].number()));
    return !pts.empty();
  }
  return false;
}

// a table: T and density or specific volume (comma, semicolon, tab or space separated; # comments). kind: what the
// header names ("density", "specific volume", "fit", "potential" or "" when it names nothing)
bool table_file(const fs::path& path, bool volume, std::vector<CoolingPoint>& pts, std::string& kind) {
  std::ifstream f(path);
  if (!f) return false;
  std::vector<std::pair<double, double>> rows;
  for (std::string line; std::getline(f, line);) {
    if (line.empty()) continue;
    const std::string l = lower(line);
    if (line[0] == '#') {   // caps tg --csv: "# Glass transition Tg: <curve>"
      const auto colon = l.rfind(':');
      std::string what = colon == std::string::npos ? "" : l.substr(colon + 1);
      what.erase(0, what.find_first_not_of(" \t"));
      while (!what.empty() && std::isspace(static_cast<unsigned char>(what.back()))) what.pop_back();
      if (what == "fit") kind = "fit";
      continue;
    }
    for (char& c : line) if (c == ',' || c == ';' || c == '\t') c = ' ';
    std::istringstream ls(line);
    double a, b;
    if (ls >> a >> b) {
      rows.emplace_back(a, b);
      continue;
    }
    if (kind != "fit") {   // a header row
      if (l.find("potential") != std::string::npos) kind = "potential";
      else if (l.find("specific volume") != std::string::npos || l.find("volume") != std::string::npos) kind = "specific volume";
      else if (l.find("density") != std::string::npos) kind = "density";
    }
  }
  if (rows.empty() || kind == "fit" || kind == "potential") return false;
  // no header: a density lies near 1 g/cm³ as a specific volume does; --volume says which, else density
  const bool vol = volume || kind == "specific volume";
  for (const auto& [T, y] : rows) pts.push_back(vol ? from_volume(T, y) : from_density(T, y));
  return true;
}

bool csv_folder(const fs::path& dir, bool volume, std::vector<CoolingPoint>& pts, std::string& source) {
  std::vector<fs::path> files;
  std::error_code ec;
  for (const auto& e : fs::directory_iterator(dir, ec))
    if (e.is_regular_file() && lower(e.path().extension().string()) == ".csv") files.push_back(e.path());
  std::sort(files.begin(), files.end());
  std::vector<CoolingPoint> best;
  std::string best_kind;
  for (const auto& f : files) {
    std::vector<CoolingPoint> p;
    std::string kind;
    if (!table_file(f, volume, p, kind)) continue;
    // the density curve first (the specific volume shown may be the property fitted), then specific volume
    const int rank = kind == "density" ? 2 : kind == "specific volume" ? 1 : 0;
    const int have = best_kind == "density" ? 2 : best_kind == "specific volume" ? 1 : best.empty() ? -1 : 0;
    if (rank > have) best = p, best_kind = kind, source = f.string();
  }
  pts = best;
  return !pts.empty();
}

// the run's own log: the one whose printed points are these
void match_log(CoolingScan& sc, const fs::path& dir) {
  std::vector<fs::path> dirs = {dir};
  if (dir.has_parent_path() && dir.parent_path() != dir) {
    dirs.push_back(dir.parent_path());
    if (dir.parent_path().has_parent_path() && dir.parent_path().parent_path() != dir.parent_path()) dirs.push_back(dir.parent_path().parent_path());
  }
  static const std::regex point(R"(^\s*([0-9.eE+-]+) K\s+([0-9.eE+-]+) ± ([0-9.eE+-]+) g/cm³)");
  static const std::regex sched(R"(stepwise cooling ([0-9.eE+-]+) → ([0-9.eE+-]+) K in ([0-9.eE+-]+) K steps of ([0-9.eE+-]+) ps)");
  static const std::regex tgline(R"(^Glass transition Tg: ([0-9.eE+-]+))");
  for (const auto& d : dirs) {
    std::vector<fs::path> logs;
    std::error_code ec;
    for (const auto& e : fs::directory_iterator(d, ec)) {
      if (!e.is_regular_file()) continue;
      const std::string ext = lower(e.path().extension().string());
      if ((ext == ".log" || ext == ".out" || ext == ".txt") && e.file_size(ec) < (5u << 20)) logs.push_back(e.path());
    }
    std::sort(logs.begin(), logs.end());
    for (const auto& lp : logs) {
      const std::string text = slurp(lp);
      if (text.find("Glass transition Tg") == std::string::npos || text.find("stepwise cooling") == std::string::npos) continue;
      std::vector<std::pair<double, double>> printed;
      std::string method;
      double tg = NaN;
      std::istringstream in(text);
      for (std::string line; std::getline(in, line);) {
        std::smatch m;
        if (std::regex_search(line, m, point)) printed.emplace_back(std::stod(m[1]), std::stod(m[2]));
        else if (std::regex_search(line, m, tgline)) tg = std::stod(m[1]);
        else if (const auto at = line.find("method: stepwise cooling"); at != std::string::npos) method = line.substr(at + 8);
      }
      if (printed.size() != sc.points.size() || method.empty()) continue;
      bool same = true;
      for (size_t k = 0; k < printed.size() && same; ++k)
        same = std::fabs(printed[k].first - sc.points[k].temperature) < 1e-6 && std::fabs(printed[k].second - sc.points[k].density) < 6e-6;
      if (!same) continue;
      sc.log = lp.string();
      sc.logged_tg = tg;
      const auto cut = method.find("; Tg =");
      sc.schedule = cut == std::string::npos ? method : method.substr(0, cut);
      std::smatch m;
      if (std::regex_search(method, m, sched)) sc.t_step = std::stod(m[3]), sc.hold_ps = std::stod(m[4]);
      return;
    }
  }
}

// the atom count of a LAMMPS data file in the folder (the potential energy is per cell; the fit takes it per atom)
size_t data_atoms(const fs::path& dir) {
  std::vector<fs::path> files;
  std::error_code ec;
  for (const auto& e : fs::directory_iterator(dir, ec))
    if (e.is_regular_file() && lower(e.path().extension().string()) == ".data") files.push_back(e.path());
  std::sort(files.begin(), files.end());
  static const std::regex atoms(R"(^\s*(\d+)\s+atoms\s*$)");
  for (const auto& f : files) {
    std::ifstream in(f);
    std::string line;
    for (int k = 0; k < 40 && std::getline(in, line); ++k) {
      std::smatch m;
      if (std::regex_search(line, m, atoms)) return size_t(std::stoull(m[1]));
    }
  }
  return 0;
}

bool read_file(const fs::path& p, bool volume, CoolingScan& sc) {
  if (lower(p.extension().string()) == ".json") {
    Json j;
    try { j = Json::parse(slurp(p)); } catch (const std::exception&) { return false; }
    if (points_file(j, sc.points) || properties_file(j, sc.points)) { sc.source = p.string(); return true; }
    return false;
  }
  std::string kind;
  if (table_file(p, volume, sc.points, kind)) { sc.source = p.string(); return true; }
  return false;
}

}  // namespace

CoolingScan read_cooling_scan(const std::string& path, bool volume) {
  CoolingScan sc;
  std::error_code ec;
  const fs::path p = fs::absolute(fs::path(path), ec).lexically_normal();
  if (!fs::exists(p, ec)) throw std::runtime_error("no such file or folder: " + path);
  fs::path dir = p;
  if (fs::is_directory(p, ec)) {
    // the run's folder: its points file, else its csv folder, else a properties file (here or in out/), else its own tables
    auto first = [&](const fs::path& d, const std::string& tail) {
      std::vector<fs::path> f;
      std::error_code e2;
      for (const auto& e : fs::directory_iterator(d, e2)) {
        const std::string n = e.path().filename().string();
        if (e.is_regular_file() && n.size() > tail.size() && n.compare(n.size() - tail.size(), tail.size(), tail) == 0) f.push_back(e.path());
      }
      std::sort(f.begin(), f.end());
      return f.empty() ? fs::path() : f.front();
    };
    bool ok = false;
    if (const auto f = first(p, ".tg-points.json"); !f.empty()) ok = read_file(f, volume, sc);
    if (!ok && fs::is_directory(p / "csv", ec)) ok = csv_folder(p / "csv", volume, sc.points, sc.source);
    for (const fs::path d : {p, p / "out"})
      if (!ok && fs::is_directory(d, ec))
        if (const auto f = first(d, ".properties.json"); !f.empty()) ok = read_file(f, volume, sc);
    if (!ok) ok = csv_folder(p, volume, sc.points, sc.source);
    if (!ok) throw std::runtime_error(path + ": no cooling scan here (OUT.tg-points.json, a csv folder of caps tg, or a properties file)");
  } else {
    if (!read_file(p, volume, sc)) throw std::runtime_error(path + ": no cooling points in it (T against density or specific volume)");
    dir = p.parent_path().empty() ? fs::path(".") : p.parent_path();
  }
  if (dir.filename() == "csv") dir = dir.parent_path();
  if (sc.points.size() < 2) throw std::runtime_error(path + ": fewer than two temperatures");
  match_log(sc, dir);
  sc.atoms = data_atoms(dir);
  return sc;
}

std::vector<Property> fit_cooling_scans(const std::vector<CoolingScan>& scans, int property, int fit, double glassy_max, double rubbery_min) {
  if (scans.empty()) throw std::invalid_argument("no cooling scans to fit");
  const bool energy = property == 1;
  const std::string what = energy ? "potential energy per atom" : "specific volume", unit = energy ? "kcal/mol" : "cm³/g";
  const std::string how = fit == 1 ? "the crossing of straight lines through the glassy and the rubbery range" : "the hinge of a continuous two-line fit";
  std::vector<CoolingResult> runs;
  std::vector<Property> out;
  for (size_t k = 0; k < scans.size(); ++k) {
    const auto& sc = scans[k];
    CoolingResult r;
    r.points = sc.points;
    r.property = property;
    std::vector<double> T, y;
    bool any_energy = false;
    for (const auto& p : sc.points) any_energy = any_energy || p.potential != 0;
    if (energy && (sc.atoms == 0 || !any_energy))
      throw std::invalid_argument(fs::path(sc.source).filename().string() + ": Tg from the potential energy needs the run's points file (it carries the energy) and its structure beside it (the atom count)");
    for (const auto& p : sc.points) { T.push_back(p.temperature); y.push_back(energy ? p.potential / double(sc.atoms) : p.specific_volume); }
    r.fitted = y;
    r.fit = fit == 1 ? fit_two_ranges(T, y, glassy_max, rubbery_min) : fit_bilinear(T, y);
    const std::string name = fs::path(sc.source).filename().string();
    r.method = (sc.schedule.empty() ? "a finished cooling scan read back from " + name + " (" + std::to_string(sc.points.size()) + " temperatures)"
                                    : sc.schedule + " (read back from " + name + ")") +
               "; Tg = " + how + " of " + what + ", error by bootstrap";
    if (sc.t_step > 0 && sc.hold_ps > 0)
      r.notes.push_back("effective cooling rate " + g(sc.t_step / sc.hold_ps * 1e12, 3) +
                        " K/s: simulated Tg sits above the experimental value, roughly 3 K per decade of rate (Williams–Landel–Ferry)");
    else
      r.notes.push_back("the cooling rate is not known: no log of the run beside the points");
    if (!r.fit.note.empty()) r.notes.push_back(r.fit.note);
    if (std::isfinite(sc.logged_tg) && r.fit.ok && fit == 0 && !energy)
      r.notes.push_back("the run itself reported Tg " + g(sc.logged_tg, 6) + " K");
    auto props = cooling_properties(r);
    if (scans.size() > 1) {
      props[0].id = "tg_r" + std::to_string(k + 1);
      props[0].name = "Glass transition Tg, replica " + std::to_string(k + 1) + " of " + std::to_string(scans.size());
    }
    out.push_back(props[0]);
    runs.push_back(std::move(r));
  }
  if (scans.size() == 1) return out;

  // the replicas pooled at the temperatures of the first: the property's mean ± SD, Tg's mean ± SD, the slopes' means
  auto sd = [](const std::vector<double>& v) {
    if (v.size() < 2) return NaN;
    double m = 0;
    for (double x : v) m += x;
    m /= double(v.size());
    double s = 0;
    for (double x : v) s += (x - m) * (x - m);
    return std::sqrt(s / double(v.size() - 1));
  };
  auto mean = [](const std::vector<double>& v) {
    double m = 0;
    for (double x : v) m += x;
    return v.empty() ? NaN : m / double(v.size());
  };
  Property tg;
  tg.id = "tg";
  tg.name = "Glass transition Tg";
  tg.unit = "K";
  std::vector<double> tgs;
  std::vector<size_t> ok;
  for (size_t k = 0; k < runs.size(); ++k)
    if (runs[k].fit.ok) tgs.push_back(runs[k].fit.tg), ok.push_back(k);
  tg.value = mean(tgs);
  tg.error = tgs.size() > 1 ? sd(tgs) : tgs.size() == 1 ? runs[ok[0]].fit.tg_err : NaN;
  tg.method = "mean ± SD of " + std::to_string(runs.size()) + " replicas, each: " + runs[0].method;
  tg.extra["replicas"] = double(runs.size());
  for (size_t k = 0; k < runs.size(); ++k) tg.extra["Tg replica " + std::to_string(k + 1) + " (K)"] = runs[k].fit.ok ? runs[k].fit.tg : NaN;
  // the extras every replica has, averaged
  if (!ok.empty()) {
    for (const auto& [key, v0] : out[ok[0]].extra) {
      std::vector<double> v;
      for (size_t k : ok)
        if (auto it = out[k].extra.find(key); it != out[k].extra.end() && std::isfinite(it->second)) v.push_back(it->second);
      if (v.size() == ok.size()) tg.extra[key] = mean(v);
    }
  }
  tg.notes = {runs[0].notes.front()};   // the cooling rate (or that it is not known)
  for (size_t k = 0; k < runs.size(); ++k)
    if (!runs[k].fit.note.empty()) tg.notes.push_back("replica " + std::to_string(k + 1) + ": " + runs[k].fit.note);
  if (ok.size() < runs.size()) tg.notes.push_back(std::to_string(runs.size() - ok.size()) + " replica(s) without a break are left out of Tg");
  Series m{what + " (mean of replicas)", "T (K)", what + " (" + unit + ")", {}, {}}, s{what + " (SD of replicas)", "T (K)", what + " (" + unit + ")", {}, {}};
  Series fitline{"fit", "T (K)", what + " (" + unit + ")", {}, {}}, d{"density", "T (K)", "density (g/cm³)", {}, {}};
  bool missing = false;
  for (size_t i = 0; i < runs[0].points.size(); ++i) {
    const double T = runs[0].points[i].temperature;
    std::vector<double> v, rho;
    for (const auto& r : runs)
      for (size_t j = 0; j < r.points.size(); ++j)
        if (std::fabs(r.points[j].temperature - T) < 1e-6) { v.push_back(r.fitted[j]); rho.push_back(r.points[j].density); break; }
    missing = missing || v.size() < runs.size();
    m.x.push_back(T), m.y.push_back(mean(v));
    s.x.push_back(T), s.y.push_back(sd(v));
    d.x.push_back(T), d.y.push_back(mean(rho));
  }
  if (missing) tg.notes.push_back("not every replica has every temperature: each mean is over the replicas that have it");
  if (!ok.empty()) {
    std::vector<double> at, lo, hi;
    for (size_t k : ok) at.push_back(runs[k].fit.value_at_tg), lo.push_back(runs[k].fit.slope_low), hi.push_back(runs[k].fit.slope_high);
    const double a = mean(at), bl = mean(lo), bh = mean(hi);
    const double tmin = *std::min_element(m.x.begin(), m.x.end()), tmax = *std::max_element(m.x.begin(), m.x.end());
    for (int k = 0; k <= 100; ++k) {
      const double T = tmin + (tmax - tmin) * k / 100;
      fitline.x.push_back(T), fitline.y.push_back(a + bl * std::min(T - tg.value, 0.0) + bh * std::max(T - tg.value, 0.0));
    }
  }
  tg.series = {m, fitline, d, s};
  out.push_back(tg);
  return out;
}

}  // namespace caps
