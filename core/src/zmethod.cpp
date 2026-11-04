// CAPS Z-method melting point (see zmethod.hpp).
#include "caps/zmethod.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <stdexcept>

#include "caps/dynamics.hpp"

namespace caps {

namespace {
constexpr double kPi = 3.14159265358979323846;

// ρ(G) = |Σ exp(i G·r)| / N, averaged over the star
double bragg_order(const std::vector<std::array<double, 3>>& g, const double* x, size_t n) {
  double sum = 0;
  for (const auto& G : g) {
    double c = 0, s = 0;
    for (size_t i = 0; i < n; ++i) {
      const double ph = G[0] * x[3 * i] + G[1] * x[3 * i + 1] + G[2] * x[3 * i + 2];
      c += std::cos(ph), s += std::sin(ph);
    }
    sum += std::hypot(c, s) / double(n);
  }
  return g.empty() ? 0 : sum / double(g.size());
}

// the star of reciprocal vectors of the cell (d ≥ 2 Å) where the perfect cell's order is largest
std::vector<std::array<double, 3>> order_star(const System& s, double& rho0) {
  const Cell& c = s.cell;
  const double V = c.volume();
  const Vec3 b[3] = {cross(c.b, c.c) * (2 * kPi / V), cross(c.c, c.a) * (2 * kPi / V), cross(c.a, c.b) * (2 * kPi / V)};
  const double gmax = 2 * kPi / 2.0;
  int hm[3];
  for (int k = 0; k < 3; ++k) {
    const Vec3 e = k == 0 ? c.a : k == 1 ? c.b : c.c;
    hm[k] = int(std::floor(gmax * norm(e) / (2 * kPi))) + 1;   // |h| ≤ G_max |a| / 2π
  }
  const size_t n = s.atoms.size();
  std::vector<double> x(3 * n);
  for (size_t i = 0; i < n; ++i)
    for (int k = 0; k < 3; ++k) x[3 * i + size_t(k)] = s.atoms[i].pos[size_t(k)];
  struct Cand { std::array<double, 3> g; double len, rho; };
  std::vector<Cand> cands;
  for (int h = -hm[0]; h <= hm[0]; ++h)
    for (int k = -hm[1]; k <= hm[1]; ++k)
      for (int l = -hm[2]; l <= hm[2]; ++l) {
        if (h == 0 && k == 0 && l == 0) continue;
        const Vec3 G = b[0] * double(h) + b[1] * double(k) + b[2] * double(l);
        const double len = norm(G);
        if (len > gmax) continue;
        // one of each ± pair
        if (h < 0 || (h == 0 && (k < 0 || (k == 0 && l < 0)))) continue;
        const std::array<double, 3> g{G[0], G[1], G[2]};
        cands.push_back({g, len, bragg_order({g}, x.data(), n)});
      }
  if (cands.empty()) throw std::invalid_argument("Z method: the cell is too small for a reciprocal vector with d ≥ 2 Å");
  double best = 0;
  for (const auto& q : cands) best = std::max(best, q.rho);
  if (best < 0.2) throw std::invalid_argument("Z method: the structure has no strong Bragg order (is it a crystal?)");
  // the shortest G that comes within 10 % of the strongest, and every G of that length that does too
  double len = 1e300;
  for (const auto& q : cands)
    if (q.rho >= 0.9 * best) len = std::min(len, q.len);
  std::vector<std::array<double, 3>> star;
  double sum = 0;
  for (const auto& q : cands)
    if (q.rho >= 0.9 * best && std::fabs(q.len - len) < 1e-6 * len) star.push_back(q.g), sum += q.rho;
  rho0 = sum / double(star.size());
  return star;
}
}  // namespace

ZResult z_method(const System& crystal, const ZMethodOptions& o) {
  ZResult r;
  if (!crystal.cell.valid()) throw std::invalid_argument("Z method: a periodic crystal");
  if (o.temperatures.empty()) throw std::invalid_argument("Z method: starting temperatures");
  const size_t n = crystal.atoms.size();
  {
    const Cell& c = crystal.cell;
    const double V = c.volume();
    const double w[3] = {V / norm(cross(c.b, c.c)), V / norm(cross(c.c, c.a)), V / norm(cross(c.a, c.b))};
    const double need = 2 * o.energy.cutoff;
    for (int k = 0; k < 3; ++k)
      if (w[k] < need) {
        int m[3];
        for (int q = 0; q < 3; ++q) m[q] = int(std::ceil(need / w[q] - 1e-9));
        throw std::invalid_argument("Z method: the cell is " + std::to_string(w[k]).substr(0, 5) + " Å wide, under twice the cut-off (" + std::to_string(need).substr(0, 5) +
                                    " Å): make a " + std::to_string(m[0]) + " × " + std::to_string(m[1]) + " × " + std::to_string(m[2]) + " supercell first");
      }
  }
  r.g = order_star(crystal, r.order0);
  std::vector<double> t0s = o.temperatures;
  std::sort(t0s.begin(), t0s.end());
  const int total_runs = int(t0s.size()) + std::max(0, o.bisect);
  int done = 0;
  auto run = [&](double t0) {
    ZPoint p;
    p.t0 = t0;
    System s = crystal;
    for (auto& a : s.atoms) a.image = {0, 0, 0};
    DynamicsOptions m;
    m.field = o.field;
    m.energy = o.energy;
    m.dt = o.dt;
    m.steps = o.steps;
    m.temperature = t0;
    m.thermostat = Thermostat::None;
    m.barostat = Barostat::None;
    m.new_velocities = true;
    m.seed = o.seed + uint64_t(1000 * t0);
    m.thermo_every = std::max(1, o.sample_every);
    m.frame_every = 0;
    std::vector<double> order;
    m.after_step = [&](const std::vector<double>& x, std::vector<double>&, const Cell&, int64_t step) {
      if (step % m.thermo_every == 0) order.push_back(bragg_order(r.g, x.data(), n) / r.order0);
    };
    const int idx = done;
    m.progress = [&](const ThermoRow& row) {
      if (!o.progress) return true;
      char b[96];
      std::snprintf(b, sizeof b, "T0 %.0f K: %.0f K", t0, row.temperature);
      return o.progress(b, (double(idx) + double(row.step) / double(std::max<int64_t>(1, o.steps))) / double(total_runs));
    };
    DynamicsReport rep;
    run_dynamics(s, m, &rep);
    const auto& th = rep.thermo;
    if (th.size() < 4) throw std::invalid_argument("Z method: too few samples (steps / sample_every)");
    double e = 0, t = 0, t2 = 0, pr = 0;
    size_t nh = 0;
    for (const auto& row : th) e += row.total;
    p.energy = e / double(th.size()) / double(n);
    p.drift = (th.back().total - th.front().total) / double(n);
    for (size_t k = th.size() / 2; k < th.size(); ++k) t += th[k].temperature, t2 += th[k].temperature * th[k].temperature, pr += th[k].pressure, ++nh;
    p.temperature = t / double(nh);
    p.temperature_sd = std::sqrt(std::max(0.0, t2 / double(nh) - p.temperature * p.temperature));
    p.pressure = pr / double(nh);
    double ord = 0;
    size_t no = 0;
    for (size_t k = order.size() - std::max<size_t>(1, order.size() / 5); k < order.size(); ++k) ord += order[k], ++no;
    p.order = no ? ord / double(no) : 0;
    p.melted = p.order < o.melt_order;
    double msd = 0;
    for (size_t i = 0; i < n; ++i) {
      const Vec3 d = s.atoms[i].pos - crystal.atoms[i].pos;
      msd += dot(d, d);
    }
    p.msd = msd / double(n);
    for (size_t k = 0; k < th.size(); ++k) p.time_ps.push_back(th[k].time_ps), p.t_trace.push_back(th[k].temperature), p.order_trace.push_back(k < order.size() ? order[k] : order.back());
    ++done;
    return p;
  };
  for (double t0 : t0s) r.points.push_back(run(t0));
  // bisection of the starting temperature between the last crystal and the first melted run
  for (int b = 0; b < o.bisect; ++b) {
    double lo = -1, hi = -1;
    std::sort(r.points.begin(), r.points.end(), [](const ZPoint& a, const ZPoint& c) { return a.t0 < c.t0; });
    for (const auto& p : r.points) {
      if (p.melted) { hi = p.t0; break; }
      lo = p.t0;
    }
    if (lo < 0 || hi < 0) break;
    r.points.push_back(run(0.5 * (lo + hi)));
  }
  std::sort(r.points.begin(), r.points.end(), [](const ZPoint& a, const ZPoint& c) { return a.energy < c.energy; });
  const ZPoint* first_melt = nullptr;
  for (const auto& p : r.points)
    if (p.melted) { first_melt = &p; break; }
  if (!first_melt) {
    r.notes.push_back("no run melted: start hotter or run longer");
  } else {
    r.t_melt = first_melt->temperature, r.e_melt = first_melt->energy;
    for (const auto& p : r.points)
      if (!p.melted && p.energy < first_melt->energy && p.temperature > r.t_superheat) r.t_superheat = p.temperature, r.e_superheat = p.energy;
    if (r.t_superheat == 0) r.notes.push_back("every run melted: start colder for the crystal branch");
    else if (r.t_melt >= r.t_superheat) r.notes.push_back("the melted run is not colder than the last crystal: no Z yet (a smaller energy step, or longer runs)");
    bool late = false;
    for (const auto& p : r.points)
      if (!p.melted && p.energy > first_melt->energy) late = true;
    if (late) r.notes.push_back("a crystal run above the first melted energy: the runs are too short for melting to start reliably");
    r.notes.push_back("superheating limit T_LS " + std::to_string(int(std::lround(r.t_superheat))) + " K, melting T_m " + std::to_string(int(std::lround(r.t_melt))) + " K at the crystal's volume");
  }
  return r;
}

Json z_method_json(const ZResult& r) {
  Json j = Json::object(), pts = Json::array(), g = Json::array(), nt = Json::array();
  auto arr = [](const std::vector<double>& v) { Json a = Json::array(); for (double x : v) a.push_back(x); return a; };
  for (const auto& p : r.points) {
    Json q = Json::object();
    q["t0"] = p.t0, q["energy"] = p.energy, q["temperature"] = p.temperature, q["temperature_sd"] = p.temperature_sd, q["pressure"] = p.pressure;
    q["order"] = p.order, q["msd"] = p.msd, q["drift"] = p.drift, q["melted"] = p.melted;
    q["time_ps"] = arr(p.time_ps), q["t_trace"] = arr(p.t_trace), q["order_trace"] = arr(p.order_trace);
    pts.push_back(q);
  }
  for (const auto& v : r.g) g.push_back(arr({v[0], v[1], v[2]}));
  for (const auto& x : r.notes) nt.push_back(x);
  j["points"] = pts, j["t_melt"] = r.t_melt, j["t_superheat"] = r.t_superheat, j["e_superheat"] = r.e_superheat, j["e_melt"] = r.e_melt;
  j["g"] = g, j["order0"] = r.order0, j["notes"] = nt;
  return j;
}

}  // namespace caps
