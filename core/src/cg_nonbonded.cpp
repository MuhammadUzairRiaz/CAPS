// CAPS non-bonded coarse-grained potentials (see caps/cg_nonbonded.hpp).
#include "caps/cg_nonbonded.hpp"

#include "caps/dynamics.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <limits>
#include <set>
#include <stdexcept>

namespace caps {

namespace {
constexpr double kB = 0.0019872067;   // kcal/(mol K)
constexpr double kPi = 3.14159265358979323846;

std::vector<double> smooth_open(const std::vector<double>& y, double s) {
  if (s <= 0) return y;
  const int w = int(std::ceil(3 * s)), n = int(y.size());
  std::vector<double> out(y.size());
  for (int i = 0; i < n; ++i) {
    double a = 0, ws = 0;
    for (int j = -w; j <= w; ++j) {
      const int q = i + j;
      if (q < 0 || q >= n) continue;
      const double k = std::exp(-0.5 * (j / s) * (j / s));
      a += k * y[size_t(q)], ws += k;
    }
    out[size_t(i)] = a / ws;
  }
  return out;
}

double interp(const std::vector<double>& x, const std::vector<double>& y, double t) {
  if (t <= x.front()) return y.front();
  if (t >= x.back()) return y.back();
  const size_t k = size_t(std::upper_bound(x.begin(), x.end(), t) - x.begin());
  const double w = (t - x[k - 1]) / (x[k] - x[k - 1]);
  return y[k - 1] * (1 - w) + y[k] * w;
}

void forces(CgPairTable& p, double dr) {
  const size_t n = p.U.size();
  p.F.assign(n, 0.0);
  for (size_t k = 0; k < n; ++k) {
    const size_t a = k == 0 ? 0 : k - 1, b = std::min(n - 1, k + 1);
    p.F[k] = -(p.U[b] - p.U[a]) / (double(b - a) * dr);
  }
}

double half_width(const Cell& c) {
  const double V = std::fabs(c.volume());
  return 0.5 * std::min({V / norm(cross(c.b, c.c)), V / norm(cross(c.c, c.a)), V / norm(cross(c.a, c.b))});
}
}  // namespace

std::vector<std::string> cg_pair_keys(const CgTypes& t) {
  std::vector<std::string> out;
  for (size_t i = 0; i < t.beads.size(); ++i)
    for (size_t j = i; j < t.beads.size(); ++j) out.push_back(cg_key({t.beads[i], t.beads[j]}));
  return out;
}

// ------------------------------------------------------------------------------------------------ g(r)
CgRdfAccumulator::CgRdfAccumulator(const CgTypes& types, const CgPairOptions& o) : types_(types), o_(o) {
  if (o.dr <= 0 || o.rmax <= o.dr) throw std::invalid_argument("g(r): give a positive bin and range");
  keys_ = cg_pair_keys(types);
}

int CgRdfAccumulator::add_system(const CgTopology& t, double weight, const std::string& name) {
  Sys s;
  s.t = t;
  s.weight = weight;
  s.name = name;
  s.ntype.assign(types_.beads.size(), 0.0);
  for (const auto& k : t.kind) {
    const auto it = std::find(types_.beads.begin(), types_.beads.end(), k);
    if (it == types_.beads.end()) throw std::invalid_argument("the shared type list has no bead type " + k);
    s.tix.push_back(int(it - types_.beads.begin()));
    s.ntype[size_t(s.tix.back())] += 1;
  }
  // exclusions: beads fewer than `exclude` bonds apart
  std::vector<std::vector<int>> nb(t.beads());
  for (const auto& [a, b] : t.bonds) nb[size_t(a)].push_back(b), nb[size_t(b)].push_back(a);
  s.excl.assign(t.beads(), {});
  for (size_t i = 0; i < t.beads(); ++i) {
    std::vector<int> frontier{int(i)};
    std::set<int> seen{int(i)};
    for (int depth = 1; depth < o_.exclude; ++depth) {
      std::vector<int> next;
      for (int u : frontier)
        for (int v : nb[size_t(u)])
          if (seen.insert(v).second) next.push_back(v);
      frontier = next;
    }
    for (int v : seen) if (v != int(i)) s.excl[i].push_back(v);
    std::sort(s.excl[i].begin(), s.excl[i].end());
  }
  const size_t nbin = size_t(std::ceil(o_.rmax / o_.dr));
  s.h.assign(keys_.size(), std::vector<double>(nbin, 0.0));
  sys_.push_back(std::move(s));
  return int(sys_.size()) - 1;
}

void CgRdfAccumulator::add_frame(int si, const std::vector<Vec3>& p, const Cell& cell) {
  Sys& s = sys_.at(size_t(si));
  if (!cell.valid()) throw std::invalid_argument("g(r) needs a periodic cell");
  if (p.size() != s.t.beads()) throw std::invalid_argument("a frame of " + s.name + " has another number of beads than its map");
  const double rc = std::min(o_.rmax, 0.999 * half_width(cell));
  if (rmax_used_ == 0 || rc < rmax_used_) rmax_used_ = rc;
  const size_t n = p.size(), ntp = types_.beads.size();
  // pair index by type pair
  std::vector<int> pix(ntp * ntp);
  {
    size_t k = 0;
    for (size_t i = 0; i < ntp; ++i)
      for (size_t j = i; j < ntp; ++j) pix[i * ntp + j] = pix[j * ntp + i] = int(k++);
  }
  // cell list on fractional coordinates
  std::vector<Vec3> f(n);
  for (size_t i = 0; i < n; ++i) {
    Vec3 u = cell.to_fractional(p[i] - cell.origin);
    for (int d = 0; d < 3; ++d) u[size_t(d)] -= std::floor(u[size_t(d)]);
    f[i] = u;
  }
  const double V = std::fabs(cell.volume());
  const double w[3] = {V / norm(cross(cell.b, cell.c)), V / norm(cross(cell.c, cell.a)), V / norm(cross(cell.a, cell.b))};
  int nc[3];
  for (int d = 0; d < 3; ++d) nc[d] = std::max(1, int(std::floor(w[d] / rc)));
  const bool brute = nc[0] < 3 || nc[1] < 3 || nc[2] < 3;
  const double rc2 = rc * rc, inv = 1.0 / o_.dr;
  auto add = [&](size_t i, size_t j) {
    const Vec3 d = cell.minimum_image(p[j] - p[i]);
    const double r2 = dot(d, d);
    if (r2 >= rc2) return;
    const auto& ex = s.excl[i];
    if (std::binary_search(ex.begin(), ex.end(), int(j))) return;
    const size_t b = size_t(std::sqrt(r2) * inv);
    if (b < s.h[0].size()) s.h[size_t(pix[size_t(s.tix[i]) * ntp + size_t(s.tix[j])])][b] += 1;
  };
  if (brute) {
    for (size_t i = 0; i < n; ++i)
      for (size_t j = i + 1; j < n; ++j) add(i, j);
  } else {
    std::vector<std::vector<uint32_t>> bins(size_t(nc[0] * nc[1] * nc[2]));
    auto cid = [&](int x, int y, int z) { return size_t(((x + nc[0]) % nc[0]) + nc[0] * (((y + nc[1]) % nc[1]) + nc[1] * ((z + nc[2]) % nc[2]))); };
    std::vector<std::array<int, 3>> cof(n);
    for (size_t i = 0; i < n; ++i) {
      cof[i] = {std::min(nc[0] - 1, int(f[i][0] * nc[0])), std::min(nc[1] - 1, int(f[i][1] * nc[1])), std::min(nc[2] - 1, int(f[i][2] * nc[2]))};
      bins[cid(cof[i][0], cof[i][1], cof[i][2])].push_back(uint32_t(i));
    }
    for (size_t i = 0; i < n; ++i)
      for (int dx = -1; dx <= 1; ++dx)
        for (int dy = -1; dy <= 1; ++dy)
          for (int dz = -1; dz <= 1; ++dz)
            for (uint32_t j : bins[cid(cof[i][0] + dx, cof[i][1] + dy, cof[i][2] + dz)])
              if (j > i) add(i, j);
  }
  s.volume += V;
  ++s.frames;
}

std::vector<double> CgRdfAccumulator::r() const {
  const size_t nbin = size_t(std::ceil(o_.rmax / o_.dr));
  std::vector<double> out(nbin);
  for (size_t k = 0; k < nbin; ++k) out[k] = (double(k) + 0.5) * o_.dr;
  return out;
}

std::vector<CgRdf> CgRdfAccumulator::rdf(int si) const {
  const Sys& s = sys_.at(size_t(si));
  std::vector<CgRdf> out;
  const size_t ntp = types_.beads.size();
  size_t k = 0;
  const size_t nbin = size_t(std::ceil(o_.rmax / o_.dr));
  const size_t nuse = std::min(nbin, size_t(rmax_used_ / o_.dr));
  for (size_t i = 0; i < ntp; ++i)
    for (size_t j = i; j < ntp; ++j, ++k) {
      CgRdf R;
      R.key = keys_[k];
      R.frames = s.frames;
      R.ideal_pairs = i == j ? s.ntype[i] * (s.ntype[i] - 1) / 2 : s.ntype[i] * s.ntype[j];
      R.counts = s.h[k];
      R.g.assign(nbin, 0.0);
      if (s.frames > 0 && R.ideal_pairs > 0) {
        const double Vm = s.volume / double(s.frames);
        for (size_t b = 0; b < nuse; ++b) {
          const double r0 = double(b) * o_.dr, r1 = r0 + o_.dr, shell = 4.0 / 3.0 * kPi * (r1 * r1 * r1 - r0 * r0 * r0);
          R.g[b] = s.h[k][b] / double(s.frames) / (R.ideal_pairs * shell / Vm);
        }
      }
      out.push_back(std::move(R));
    }
  return out;
}

// ------------------------------------------------------------------------------------------------ targets, tables (JSON)
namespace {
Json arr(const std::vector<double>& v) {
  Json a = Json::array();
  for (double x : v) a.push_back(x);
  return a;
}
std::vector<double> vec(const Json& a) {
  std::vector<double> v;
  for (const auto& x : a.items()) v.push_back(x.number());
  return v;
}
}  // namespace

CgTargets cg_targets(const CgRdfAccumulator& acc, double temperature, const std::vector<double>& pressure) {
  CgTargets t;
  t.r = acc.r();
  t.temperature = temperature;
  for (int s = 0; s < acc.systems(); ++s) {
    t.systems.push_back(acc.name(s));
    t.weights.push_back(acc.weight(s));
    t.rdf.push_back(acc.rdf(s));
    t.pressure.push_back(size_t(s) < pressure.size() ? pressure[size_t(s)] : 1.0);
  }
  return t;
}

Json cg_targets_json(const CgTargets& t) {
  Json j = Json::object();
  j["format"] = "caps-cg-targets";
  j["version"] = 1;
  j["temperature"] = t.temperature;
  j["r"] = arr(t.r);
  Json sys = Json::array();
  for (size_t s = 0; s < t.systems.size(); ++s) {
    Json e = Json::object();
    e["name"] = t.systems[s];
    e["weight"] = t.weights[s];
    e["pressure"] = t.pressure[s];
    Json pairs = Json::array();
    for (const auto& R : t.rdf[s]) {
      Json p = Json::object();
      p["key"] = R.key;
      p["frames"] = double(R.frames);
      p["ideal_pairs"] = R.ideal_pairs;
      p["g"] = arr(R.g);
      p["counts"] = arr(R.counts);
      pairs.push_back(p);
    }
    e["pairs"] = pairs;
    sys.push_back(e);
  }
  j["systems"] = sys;
  return j;
}

CgTargets cg_targets_from_json(const Json& j) {
  if (j.text("format") != "caps-cg-targets") throw std::invalid_argument("not CAPS g(r) targets (format caps-cg-targets)");
  CgTargets t;
  t.temperature = j.num("temperature", 300);
  t.r = vec(j["r"]);
  for (const auto& e : j["systems"].items()) {
    t.systems.push_back(e.text("name"));
    t.weights.push_back(e.num("weight", 1));
    t.pressure.push_back(e.num("pressure", 1));
    std::vector<CgRdf> v;
    for (const auto& p : e["pairs"].items()) {
      CgRdf R;
      R.key = p["key"].str();
      R.frames = long(p.num("frames", 0));
      R.ideal_pairs = p.num("ideal_pairs", 0);
      R.g = vec(p["g"]);
      R.counts = vec(p["counts"]);
      v.push_back(std::move(R));
    }
    t.rdf.push_back(std::move(v));
  }
  return t;
}

Json cg_pairs_json(const CgPairSet& p) {
  Json j = Json::object();
  j["format"] = "caps-cg-pairs";
  j["version"] = 1;
  j["r0"] = p.r0;
  j["dr"] = p.dr;
  j["rc"] = p.rc;
  j["temperature"] = p.temperature;
  j["iteration"] = p.iteration;
  j["source"] = p.source;
  Json a = Json::array();
  for (const auto& t : p.pairs) {
    Json e = Json::object();
    e["key"] = t.key;
    e["U"] = arr(t.U);
    a.push_back(e);
  }
  j["pairs"] = a;
  return j;
}

CgPairSet cg_pairs_from_json(const Json& j) {
  if (j.text("format") != "caps-cg-pairs") throw std::invalid_argument("not a CAPS pair set (format caps-cg-pairs)");
  CgPairSet p;
  p.r0 = j.num("r0", 0.5);
  p.dr = j.num("dr", 0.05);
  p.rc = j.num("rc", 15);
  p.temperature = j.num("temperature", 300);
  p.iteration = int(j.num("iteration", 0));
  p.source = j.text("source");
  for (const auto& e : j["pairs"].items()) {
    CgPairTable t;
    t.key = e["key"].str();
    t.U = vec(e["U"]);
    forces(t, p.dr);
    p.pairs.push_back(std::move(t));
  }
  return p;
}

// ------------------------------------------------------------------------------------------------ IBI
namespace {

// the potential on the targets' grid from g (−kT ln g where resolved), a wall below, zero at rc
std::vector<double> potential_of(const std::vector<double>& r, const std::vector<double>& g, double kT, double rc, size_t* first) {
  const size_t n = r.size();
  size_t k0 = 0;
  while (k0 < n && g[k0] < 0.02) ++k0;
  if (k0 + 5 >= n || r[k0] >= rc) throw std::invalid_argument("a target g(r) is not resolved inside the cut-off");
  std::vector<double> U(n, 0.0);
  for (size_t k = k0; k < n; ++k) U[k] = -kT * std::log(std::max(g[k], 1e-6));
  const double dr = r[1] - r[0];
  const double slope = std::max(5 * kT, (U[k0] - U[k0 + 1]) / dr);
  for (size_t k = 0; k < k0; ++k) U[k] = U[k0] + slope * double(k0 - k) * dr;
  if (first) *first = k0;
  return U;
}

// onto the table grid, a straight wall below the target grid's first point, shifted to zero at rc
CgPairTable to_table(const std::string& key, const std::vector<double>& r, const std::vector<double>& U, const CgPairSet& p) {
  CgPairTable t;
  t.key = key;
  const size_t n = size_t(std::lround((p.rc - p.r0) / p.dr)) + 1;
  t.U.resize(n);
  const double dr = r[1] - r[0];
  const double slope = (U[0] - U[1]) / dr;
  for (size_t k = 0; k < n; ++k) {
    const double x = p.r0 + double(k) * p.dr;
    t.U[k] = x < r.front() ? U.front() + slope * (r.front() - x) : interp(r, U, x);
  }
  const double end = t.U.back();
  for (auto& u : t.U) u -= end;
  forces(t, p.dr);
  return t;
}

// the table back on the targets' grid
std::vector<double> on_grid(const CgPairTable& t, const CgPairSet& p, const std::vector<double>& r) {
  std::vector<double> x(t.U.size());
  for (size_t k = 0; k < x.size(); ++k) x[k] = p.r0 + double(k) * p.dr;
  std::vector<double> out(r.size());
  for (size_t k = 0; k < r.size(); ++k) out[k] = r[k] > p.rc ? 0.0 : interp(x, t.U, r[k]);
  return out;
}

size_t pair_index(const std::vector<CgRdf>& v, const std::string& key) {
  for (size_t k = 0; k < v.size(); ++k) if (v[k].key == key) return k;
  return SIZE_MAX;
}
}  // namespace

CgPairSet ibi_start(const CgTargets& t, const CgIbiOptions& o) {
  if (t.rdf.empty()) throw std::invalid_argument("no targets");
  CgPairSet p;
  p.temperature = t.temperature;
  p.dr = o.table_dr;
  p.rc = std::min(o.rc, t.r.back());
  p.r0 = std::max(0.5, p.dr);
  p.source = "IBI";
  const double kT = kB * t.temperature;
  for (const auto& R0 : t.rdf.front()) {
    const std::string& key = R0.key;
    // pooled over the systems that have the pair, weighted by their counts there and their weight
    std::vector<double> g(t.r.size(), 0.0), w(t.r.size(), 0.0);
    for (size_t s = 0; s < t.rdf.size(); ++s) {
      const size_t k = pair_index(t.rdf[s], key);
      if (k == SIZE_MAX || t.rdf[s][k].ideal_pairs <= 0) continue;
      // statistics per bin: the ideal count (smooth in r) — counts would weigh against the low g a pair has near contact
      for (size_t b = 0; b < t.r.size(); ++b) {
        const double ws = t.weights[s] * t.rdf[s][k].ideal_pairs;
        g[b] += ws * t.rdf[s][k].g[b], w[b] += ws;
      }
    }
    bool any = false;
    for (size_t b = 0; b < g.size(); ++b) if (w[b] > 0) g[b] /= w[b], any = true;
    if (!any) continue;   // a pair no system has: no table (said by the caller)
    const auto U = potential_of(t.r, g, kT, p.rc, nullptr);
    p.pairs.push_back(to_table(key, t.r, U, p));
  }
  return p;
}

CgIbiStepReport ibi_compare(const CgTargets& t, const std::vector<std::vector<CgRdf>>& cur, const std::vector<double>& pressure) {
  CgIbiStepReport rep;
  rep.pressure = pressure;
  for (const auto& R0 : t.rdf.front()) {
    CgIbiPairReport pr;
    pr.key = R0.key;
    double num = 0, den = 0, wsum = 0;
    for (size_t s = 0; s < t.rdf.size() && s < cur.size(); ++s) {
      const size_t k = pair_index(t.rdf[s], R0.key), kc = pair_index(cur[s], R0.key);
      if (k == SIZE_MAX || kc == SIZE_MAX || t.rdf[s][k].ideal_pairs <= 0) continue;
      double n = 0, d = 0;
      for (size_t b = 0; b < t.r.size(); ++b) {
        const double e = cur[s][kc].g[b] - t.rdf[s][k].g[b];
        n += e * e, d += t.rdf[s][k].g[b] * t.rdf[s][k].g[b];
        pr.max_dev = std::max(pr.max_dev, std::fabs(e));
      }
      num += t.weights[s] * n, den += t.weights[s] * d, wsum += t.weights[s];
    }
    if (wsum <= 0) continue;
    pr.residual = den > 0 ? num / den : 0;
    rep.residual = std::max(rep.residual, pr.residual);
    rep.pairs.push_back(pr);
  }
  return rep;
}

CgPairSet ibi_step(const CgPairSet& cur, const CgTargets& t, const std::vector<std::vector<CgRdf>>& current, const std::vector<double>& pressure,
                   const CgIbiOptions& o, CgIbiStepReport* rep_out) {
  if (current.size() != t.rdf.size()) throw std::invalid_argument("one CG g(r) set per target system");
  const double kT = kB * t.temperature;
  CgIbiStepReport rep = ibi_compare(t, current, pressure);
  CgPairSet next = cur;
  ++next.iteration;
  // the pressure correction: the weighted mean error (atm → bar)
  double A = 0;
  if (o.pressure_correction && pressure.size() == t.rdf.size()) {
    double dp = 0, ws = 0;
    for (size_t s = 0; s < pressure.size(); ++s) dp += t.weights[s] * (pressure[s] - t.pressure[s]), ws += t.weights[s];
    dp = ws > 0 ? dp / ws * 1.01325 : 0;
    A = -o.ramp * kT * (dp > 0 ? 1 : dp < 0 ? -1 : 0) * std::min(1.0, 0.0003 * std::fabs(dp));
  }
  rep.ramp_A = A;
  for (auto& tab : next.pairs) {
    std::vector<double> num(t.r.size(), 0.0), den(t.r.size(), 0.0);
    for (size_t s = 0; s < t.rdf.size(); ++s) {
      const size_t k = pair_index(t.rdf[s], tab.key), kc = pair_index(current[s], tab.key);
      if (k == SIZE_MAX || kc == SIZE_MAX || t.rdf[s][k].ideal_pairs <= 0) continue;
      const double ws = t.weights[s] * t.rdf[s][k].ideal_pairs;
      for (size_t b = 0; b < t.r.size(); ++b) {
        const double gt = t.rdf[s][k].g[b], gc = current[s][kc].g[b];
        if (gt > 1e-3 && gc > 1e-3) num[b] += ws * std::log(gc / gt), den[b] += ws;
      }
    }
    std::vector<double> dU(t.r.size(), 0.0);
    size_t first = SIZE_MAX;
    for (size_t b = 0; b < t.r.size(); ++b)
      if (den[b] > 0) { dU[b] = o.alpha * kT * num[b] / den[b]; if (first == SIZE_MAX) first = b; }
    if (first == SIZE_MAX) continue;
    // below the first resolved bin the wall moves with it; smoothed; then the ramp
    for (size_t b = 0; b < first; ++b) dU[b] = dU[first];
    dU = smooth_open(dU, o.smooth);
    std::vector<double> U = on_grid(tab, cur, t.r);
    for (size_t b = 0; b < t.r.size(); ++b)
      if (t.r[b] <= next.rc) U[b] += dU[b] + A * (1 - t.r[b] / next.rc);
    // the wall below the targets' grid follows from the first two points (to_table)
    tab = to_table(tab.key, t.r, U, next);
  }
  if (rep_out) *rep_out = rep;
  return next;
}

std::vector<std::string> write_pairs(const CgPairSet& p, const CgTypes& types, const std::string& dir) {
  namespace fs = std::filesystem;
  fs::create_directories(dir);
  std::vector<std::string> files;
  char line[256];
  {
    const std::string path = (fs::path(dir) / "pairs.table").string();
    std::ofstream f(path);
    if (!f) throw std::runtime_error("cannot write " + path);
    f << "# CAPS coarse-grained pair potentials (" << p.source << ", iteration " << p.iteration << ", " << p.temperature << " K): r (Å), E (kcal/mol), −dE/dr\n";
    for (const auto& t : p.pairs) {
      std::snprintf(line, sizeof line, "\n%s\nN %zu R %.6f %.6f\n\n", t.key.c_str(), t.U.size(), p.r0, p.r0 + double(t.U.size() - 1) * p.dr);
      f << line;
      for (size_t k = 0; k < t.U.size(); ++k) {
        std::snprintf(line, sizeof line, "%zu %.6f %.8g %.8g\n", k + 1, p.r0 + double(k) * p.dr, t.U[k], t.F[k]);
        f << line;
      }
    }
    files.push_back(path);
  }
  {
    const std::string path = (fs::path(dir) / "pair.in").string();
    std::ofstream f(path);
    f << "# CAPS coarse-grained pair potentials: include after read_data; type numbers follow the shared type list\n";
    f << "# variable PAIR is the folder of pairs.table (e.g. -var PAIR " << fs::path(dir).filename().string() << ")\n";
    f << "pair_style table linear " << p.points() << "\n";
    for (size_t i = 0; i < types.beads.size(); ++i)
      for (size_t j = i; j < types.beads.size(); ++j) {
        const std::string key = cg_key({types.beads[i], types.beads[j]});
        const bool have = std::any_of(p.pairs.begin(), p.pairs.end(), [&](const CgPairTable& t) { return t.key == key; });
        if (!have) { f << "# no table for " << key << " (no system has both types): give it a potential before a run\n"; continue; }
        std::snprintf(line, sizeof line, "pair_coeff %zu %zu ${PAIR}/pairs.table %s %.4f\n", i + 1, j + 1, key.c_str(), p.rc);
        f << line;
      }
    files.push_back(path);
  }
  {
    const std::string path = (fs::path(dir) / "pairs.json").string();
    std::ofstream f(path);
    f << cg_pairs_json(p).dump(0) << "\n";
    files.push_back(path);
  }
  return files;
}

// ------------------------------------------------------------------------------------------------ analytic fits
double pair_energy(const CgPairFit& f, double r) {
  if (f.form == "lj126") { const double x = std::pow(f.sigma / r, 6); return 4 * f.epsilon * (x * x - x); }
  if (f.form == "lj96") { const double x = std::pow(f.sigma / r, 3); return f.epsilon * (2 * x * x * x - 3 * x * x); }
  if (f.form == "morse") { const double e = std::exp(-f.a * (r - f.sigma)); return f.epsilon * (e * e - 2 * e); }
  if (f.form == "mie") {
    const double C = (f.n / (f.n - f.m)) * std::pow(f.n / f.m, f.m / (f.n - f.m));
    return C * f.epsilon * (std::pow(f.sigma / r, f.n) - std::pow(f.sigma / r, f.m));
  }
  throw std::invalid_argument("pair form: lj126, lj96, morse or mie");
}

std::vector<CgPairFit> fit_pairs(const CgPairSet& p, const std::string& form, double repulsion) {
  if (form != "lj126" && form != "lj96" && form != "morse" && form != "mie") throw std::invalid_argument("pair form: lj126, lj96, morse or mie");
  const double kT = kB * p.temperature;
  std::vector<CgPairFit> out;
  for (const auto& t : p.pairs) {
    std::vector<double> x, y;
    for (size_t k = 0; k < t.U.size(); ++k) {
      const double r = p.r0 + double(k) * p.dr;
      if (t.U[k] < repulsion * kT) x.push_back(r), y.push_back(t.U[k]);
    }
    // the repulsion side starts where U first drops below the limit: drop points before it
    while (!x.empty() && x.size() > 1 && x[1] - x[0] > 1.5 * p.dr) x.erase(x.begin()), y.erase(y.begin());
    CgPairFit f;
    f.key = t.key;
    f.form = form;
    if (x.size() < 8) { f.rms = std::numeric_limits<double>::quiet_NaN(); out.push_back(f); continue; }
    f.lo = x.front(), f.hi = x.back();
    const size_t imin = size_t(std::min_element(y.begin(), y.end()) - y.begin());
    const double depth = std::max(0.05 * kT, -y[imin]), rmin = x[imin];
    // parameters (log scale for the positive ones)
    std::vector<double> q;
    if (form == "lj126") q = {std::log(depth), std::log(rmin / std::pow(2.0, 1.0 / 6))};
    else if (form == "lj96") q = {std::log(depth), std::log(rmin)};
    else if (form == "morse") q = {std::log(depth), std::log(rmin), std::log(1.0)};
    else q = {std::log(depth), std::log(rmin / std::pow(2.0, 1.0 / 6)), std::log(12.0 - 6.0)};
    const double rc = p.rc;
    auto make = [&](const std::vector<double>& v) {
      CgPairFit g = f;
      g.epsilon = std::exp(v[0]), g.sigma = std::exp(v[1]);
      if (form == "morse") g.a = std::exp(v[2]);
      if (form == "mie") g.n = 6.0 + std::exp(v[2]), g.m = 6.0;
      return g;
    };
    auto resid = [&](const std::vector<double>& v) {
      const CgPairFit g = make(v);
      const double shift = pair_energy(g, rc);   // the table is zero at rc: so is the shifted analytic form
      std::vector<double> e(x.size());
      for (size_t i = 0; i < x.size(); ++i) e[i] = pair_energy(g, x[i]) - shift - y[i];
      return e;
    };
    auto cost = [&](const std::vector<double>& v) { double c = 0; for (double e : resid(v)) c += e * e; return c; };
    // Levenberg–Marquardt with a numerical Jacobian
    double lambda = 1e-2, c0 = cost(q);
    const size_t np = q.size();
    for (int it = 0; it < 300; ++it) {
      const auto e0 = resid(q);
      std::vector<std::vector<double>> J(x.size(), std::vector<double>(np));
      for (size_t a = 0; a < np; ++a) {
        auto qp = q, qm = q;
        qp[a] += 1e-6, qm[a] -= 1e-6;
        const auto ep = resid(qp), em = resid(qm);
        for (size_t i = 0; i < x.size(); ++i) J[i][a] = (ep[i] - em[i]) / 2e-6;
      }
      std::vector<std::vector<double>> H(np, std::vector<double>(np, 0.0));
      std::vector<double> gvec(np, 0.0);
      for (size_t i = 0; i < x.size(); ++i)
        for (size_t a = 0; a < np; ++a) {
          gvec[a] += J[i][a] * e0[i];
          for (size_t b = 0; b < np; ++b) H[a][b] += J[i][a] * J[i][b];
        }
      bool improved = false;
      for (int tries = 0; tries < 12 && !improved; ++tries) {
        auto M = H;
        for (size_t a = 0; a < np; ++a) M[a][a] *= 1 + lambda;
        // solve M d = −g (Gauss–Jordan, tiny system)
        std::vector<double> d(np);
        auto A = M;
        std::vector<double> b(np);
        for (size_t a = 0; a < np; ++a) b[a] = -gvec[a];
        bool singular = false;
        for (size_t c = 0; c < np && !singular; ++c) {
          size_t piv = c;
          for (size_t r = c + 1; r < np; ++r) if (std::fabs(A[r][c]) > std::fabs(A[piv][c])) piv = r;
          if (std::fabs(A[piv][c]) < 1e-300) { singular = true; break; }
          std::swap(A[c], A[piv]), std::swap(b[c], b[piv]);
          for (size_t r = 0; r < np; ++r) {
            if (r == c) continue;
            const double m = A[r][c] / A[c][c];
            for (size_t k = c; k < np; ++k) A[r][k] -= m * A[c][k];
            b[r] -= m * b[c];
          }
        }
        if (singular) { lambda *= 10; continue; }
        for (size_t a = 0; a < np; ++a) d[a] = b[a] / A[a][a];
        auto qn = q;
        for (size_t a = 0; a < np; ++a) qn[a] += std::clamp(d[a], -1.0, 1.0);
        const double cn = cost(qn);
        if (std::isfinite(cn) && cn < c0) { q = qn, improved = (c0 - cn) > 1e-14 * (1 + c0), c0 = cn, lambda = std::max(1e-9, lambda / 3); if (!improved) break; }
        else lambda *= 10;
      }
      if (!improved) break;
    }
    f = make(q);
    f.lo = x.front(), f.hi = x.back();
    f.rms = std::sqrt(c0 / double(x.size()));
    out.push_back(f);
  }
  return out;
}

CgPairSet tables_of_fits(const std::vector<CgPairFit>& fits, const CgPairSet& grid) {
  CgPairSet p = grid;
  p.pairs.clear();
  p.source = fits.empty() ? "fit" : fits.front().form + " fit";
  for (const auto& f : fits) {
    if (!std::isfinite(f.rms)) continue;
    CgPairTable t;
    t.key = f.key;
    const size_t n = grid.points() ? grid.points() : size_t(std::lround((grid.rc - grid.r0) / grid.dr)) + 1;
    const double shift = pair_energy(f, grid.rc);
    for (size_t k = 0; k < n; ++k) t.U.push_back(pair_energy(f, grid.r0 + double(k) * grid.dr) - shift);
    forces(t, grid.dr);
    p.pairs.push_back(std::move(t));
  }
  return p;
}

std::string lammps_pair_lines(const std::vector<CgPairFit>& fits, const CgTypes& types, double rc) {
  if (fits.empty()) return "";
  const std::string form = fits.front().form;
  char b[256];
  std::string out = "# CAPS coarse-grained pairs fitted by " + form + " (shifted to zero at the cut-off, as the tables they came from)\n";
  const char* style = form == "lj126" ? "lj/cut" : form == "lj96" ? "lj/class2" : form == "morse" ? "morse" : "mie/cut";
  std::snprintf(b, sizeof b, "pair_style %s %.4f\npair_modify shift yes\n", style, rc);
  out += b;
  for (size_t i = 0; i < types.beads.size(); ++i)
    for (size_t j = i; j < types.beads.size(); ++j) {
      const std::string key = cg_key({types.beads[i], types.beads[j]});
      const auto it = std::find_if(fits.begin(), fits.end(), [&](const CgPairFit& f) { return f.key == key && std::isfinite(f.rms); });
      if (it == fits.end()) { out += "# no fit for " + key + "\n"; continue; }
      if (form == "morse") std::snprintf(b, sizeof b, "pair_coeff %zu %zu %.6f %.6f %.6f  # %s, rms %.3f kcal/mol\n", i + 1, j + 1, it->epsilon, it->a, it->sigma, key.c_str(), it->rms);
      else if (form == "mie") std::snprintf(b, sizeof b, "pair_coeff %zu %zu %.6f %.6f %.4f %.4f  # %s, rms %.3f kcal/mol\n", i + 1, j + 1, it->epsilon, it->sigma, it->n, it->m, key.c_str(), it->rms);
      else std::snprintf(b, sizeof b, "pair_coeff %zu %zu %.6f %.6f  # %s, rms %.3f kcal/mol\n", i + 1, j + 1, it->epsilon, it->sigma, key.c_str(), it->rms);
      out += b;
    }
  return out;
}

// ------------------------------------------------------------------------------------------------ calibration
CgCalibrationStep calibration_step(const std::vector<CgCalibrationPoint>& h, double rho_t, double tg_t, double rho_tol, double tg_tol) {
  if (h.empty()) throw std::invalid_argument("calibration needs at least one run's density (and T_g)");
  CgCalibrationStep st;
  const auto& last = h.back();
  st.s_sigma = last.s_sigma, st.s_eps = last.s_eps;
  // density: log ρ against log s_σ (ρ ∝ s_σ⁻³ for a start)
  if (rho_t > 0 && last.density > 0) {
    if (std::fabs(last.density - rho_t) <= rho_tol * rho_t) st.density_done = true;
    else {
      double slope = -3;
      for (size_t k = h.size() - 1; k-- > 0;)
        if (h[k].density > 0 && std::fabs(std::log(h[k].s_sigma / last.s_sigma)) > 1e-6) {
          const double s = (std::log(last.density) - std::log(h[k].density)) / (std::log(last.s_sigma) - std::log(h[k].s_sigma));
          if (s < -0.3) slope = s;   // keep the physical sign; a flat or wrong-signed secant falls back to −3
          break;
        }
      const double step = (std::log(rho_t) - std::log(last.density)) / slope;
      st.s_sigma = last.s_sigma * std::exp(std::clamp(step, std::log(0.8), std::log(1.25)));
    }
  }
  // T_g: log T_g against log s_ε (T_g ∝ s_ε for a start), with σ held where the density was reached
  if (tg_t > 0 && last.tg > 0) {
    if (std::fabs(last.tg - tg_t) <= tg_tol) st.tg_done = true;
    else {
      double slope = 1;
      for (size_t k = h.size() - 1; k-- > 0;)
        if (h[k].tg > 0 && std::fabs(std::log(h[k].s_eps / last.s_eps)) > 1e-6) {
          const double s = (std::log(last.tg) - std::log(h[k].tg)) / (std::log(last.s_eps) - std::log(h[k].s_eps));
          if (s > 0.2) slope = s;
          break;
        }
      const double step = (std::log(tg_t) - std::log(last.tg)) / slope;
      st.s_eps = last.s_eps * std::exp(std::clamp(step, std::log(0.7), std::log(1.4)));
    }
  }
  char b[256];
  std::snprintf(b, sizeof b, "density %s (%.4g → target %.4g g/cm³) · T_g %s (%.4g → target %.4g K) · next s_σ %.4f, s_ε %.4f", st.density_done ? "reached" : "not yet",
                last.density, rho_t, st.tg_done ? "reached" : last.tg > 0 ? "not yet" : "not measured", last.tg, tg_t, st.s_sigma, st.s_eps);
  st.note = b;
  return st;
}

std::vector<CgPairFit> scale_fits(const std::vector<CgPairFit>& fits, double s_sigma, double s_eps) {
  std::vector<CgPairFit> out = fits;
  for (auto& f : out) {
    f.sigma *= s_sigma;
    f.epsilon *= s_eps;
    if (f.form == "morse") f.a /= s_sigma;
  }
  return out;
}

// ------------------------------------------------------------------------------------------------ CAPS's engine
ForceField cg_forcefield(const CgTopology& t, const CgTypes& types, const std::map<std::string, double>& masses, const CgPairSet& pairs,
                         const CgBondedResult* bonded, int exclude) {
  const size_t n = t.beads();
  ForceField F;
  F.name = "coarse-grained beads · pairs " + pairs.source + (bonded ? " · bonded springs at the tables' wells" : "");
  F.type_names = types.beads;
  F.lj.assign(types.beads.size(), PairType{0, 1});
  F.pair_form = "lj12-6";
  F.mixing = "arithmetic";
  F.cutoff = pairs.rc;
  F.lj14 = 1.0, F.coul14 = 1.0;
  F.type_index.resize(n);
  F.atom_type = t.kind;
  F.why.assign(n, "coarse-grained bead");
  F.charge.assign(n, 0.0);
  F.mass.resize(n);
  for (size_t b = 0; b < n; ++b) {
    const auto it = std::find(types.beads.begin(), types.beads.end(), t.kind[b]);
    if (it == types.beads.end()) throw std::invalid_argument("the shared type list has no bead type " + t.kind[b]);
    F.type_index[b] = int(it - types.beads.begin());
    const auto m = masses.find(t.kind[b]);
    F.mass[b] = m != masses.end() ? m->second : 100.0;
  }
  // pairs: one table per type pair
  for (size_t i = 0; i < types.beads.size(); ++i)
    for (size_t j = i; j < types.beads.size(); ++j) {
      const std::string key = cg_key({types.beads[i], types.beads[j]});
      const auto it = std::find_if(pairs.pairs.begin(), pairs.pairs.end(), [&](const CgPairTable& p) { return p.key == key; });
      if (it == pairs.pairs.end()) continue;
      TabulatedPair tb;
      tb.r0 = pairs.r0, tb.dr = pairs.dr, tb.e = it->U, tb.f = it->F;
      F.pair_func[{int(i), int(j)}] = PairFunc{kPairTable, double(F.tables.size()), 0, 0};
      F.tables.push_back(std::move(tb));
    }
  // bonded: harmonic at the wells
  if (!t.bonds.empty() && !bonded) throw std::invalid_argument("bonded beads need bonded potentials (cgfit bonded)");
  auto find = [](const std::vector<CgBondedTable>& v, const std::string& k) -> const CgBondedTable* {
    for (const auto& x : v) if (x.key == k && x.sampled) return &x;
    return nullptr;
  };
  for (const auto& [a, b] : t.bonds) {
    const auto* T = find(bonded->bonds, cg_key({t.kind[size_t(a)], t.kind[size_t(b)]}));
    if (!T) throw std::invalid_argument("no bonded table for bond " + cg_key({t.kind[size_t(a)], t.kind[size_t(b)]}));
    F.bonds.push_back({uint32_t(a), uint32_t(b), std::max(T->k_harmonic, 1.0), T->x0});
  }
  for (const auto& q : t.angles) {
    const auto* T = bonded ? find(bonded->angles, cg_key({t.kind[size_t(q[0])], t.kind[size_t(q[1])], t.kind[size_t(q[2])]})) : nullptr;
    if (!T) continue;
    const double k = T->k_harmonic * (180.0 / kPi) * (180.0 / kPi);   // per deg² → per rad²
    F.angles.push_back({uint32_t(q[0]), uint32_t(q[1]), uint32_t(q[2]), k, T->x0 * kPi / 180.0});
  }
  if (!t.dihedrals.empty()) F.notes.push_back("CAPS's engine runs these beads without their dihedral tables (LAMMPS uses them all)");
  // exclusions
  std::vector<std::vector<int>> nb(n);
  for (const auto& [a, b] : t.bonds) nb[size_t(a)].push_back(b), nb[size_t(b)].push_back(a);
  F.excluded.assign(n, {});
  for (size_t i = 0; i < n; ++i) {
    std::vector<int> frontier{int(i)};
    std::set<int> seen{int(i)};
    for (int d = 1; d < exclude; ++d) {
      std::vector<int> next;
      for (int u : frontier) for (int v : nb[size_t(u)]) if (seen.insert(v).second) next.push_back(v);
      frontier = next;
    }
    for (int v : seen) if (v != int(i)) F.excluded[i].push_back(uint32_t(v));
    std::sort(F.excluded[i].begin(), F.excluded[i].end());
  }
  F.native_timestep = 5;
  return F;
}

CgEngineRun run_cg_engine(const System& beads, const ForceField& ff, const CgEngineOptions& o) {
  if (!beads.cell.valid()) throw std::invalid_argument("CG runs need a periodic cell");
  CgEngineRun out;
  System s = beads;
  DynamicsOptions d;
  d.field = std::make_shared<ForceField>(ff);
  d.energy.cutoff = ff.cutoff;
  d.energy.coulomb = false;
  d.energy.tail = false;
  d.dt = o.dt;
  d.temperature = o.temperature;
  d.thermostat = Thermostat::Bussi;
  d.tau_t = 100 * o.dt;
  d.seed = o.seed;
  d.new_velocities = o.new_velocities;
  if (o.npt) d.barostat = Barostat::Berendsen, d.pressure = o.pressure;
  if (o.equilibrate > 0) {
    d.steps = o.equilibrate;
    run_dynamics(s, d);
    d.new_velocities = false;
  }
  d.steps = o.steps;
  d.frame_every = int(std::max<int64_t>(1, o.frame_every));
  d.thermo_every = d.frame_every;
  double ps = 0, rho = 0;
  long np = 0;
  d.frame = [&](const std::vector<double>& x, const Cell& c, int64_t step) {
    if (step == 0) return;
    std::vector<Vec3> p(x.size() / 3);
    for (size_t i = 0; i < p.size(); ++i) p[i] = {x[3 * i], x[3 * i + 1], x[3 * i + 2]};
    out.frames.push_back(std::move(p));
    out.cells.push_back(c);
  };
  d.progress = [&](const ThermoRow& r) { if (r.step > 0) ps += r.pressure, rho += r.density, ++np; return true; };
  run_dynamics(s, d);
  out.pressure = np ? ps / double(np) : 0;
  out.density = np ? rho / double(np) : 0;
  out.last = s;
  return out;
}

}  // namespace caps
