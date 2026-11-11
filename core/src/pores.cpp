// CAPS pore network (see caps/pores.hpp).
#include "caps/pores.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <deque>
#include <stdexcept>
#include <unordered_map>

#include "caps/domains.hpp"
#include "caps/elements.hpp"
#include "caps/voids.hpp"

namespace caps {

namespace {

size_t gidx(const int n[3], int i, int j, int k) { return (size_t(i) * size_t(n[1]) + size_t(j)) * size_t(n[2]) + size_t(k); }

// Along which axes some connected set of grid points with value > r reaches its own periodic image; `best` gets the
// largest value in a cluster percolating along `axis` (−1: any axis).
std::array<bool, 3> percolation(const GridField& g, double r, int axis, double* best) {
  const size_t N = g.values.size();
  std::vector<int> lab(N, -1);
  std::vector<std::array<int, 3>> shift(N);
  std::array<bool, 3> any{false, false, false};
  if (best) *best = 0;
  int L = 0;
  for (size_t s0 = 0; s0 < N; ++s0) {
    if (lab[s0] >= 0 || !(g.values[s0] > r)) continue;
    std::array<bool, 3> perc{false, false, false};
    double top = 0;
    std::deque<size_t> q{s0};
    lab[s0] = L, shift[s0] = {0, 0, 0};
    while (!q.empty()) {
      const size_t cur = q.front();
      q.pop_front();
      top = std::max(top, g.values[cur]);
      const int i = int(cur / (size_t(g.n[1]) * size_t(g.n[2]))), j = int((cur / size_t(g.n[2])) % size_t(g.n[1])), k = int(cur % size_t(g.n[2]));
      const int d[6][3] = {{1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}};
      for (const auto& dd : d) {
        int x[3] = {i + dd[0], j + dd[1], k + dd[2]};
        std::array<int, 3> sh = shift[cur];
        for (int a = 0; a < 3; ++a) {
          if (x[a] < 0) x[a] += g.n[a], sh[size_t(a)] -= 1;
          else if (x[a] >= g.n[a]) x[a] -= g.n[a], sh[size_t(a)] += 1;
        }
        const size_t nb = gidx(g.n, x[0], x[1], x[2]);
        if (!(g.values[nb] > r)) continue;
        if (lab[nb] < 0) {
          lab[nb] = L, shift[nb] = sh;
          q.push_back(nb);
        } else {
          for (int a = 0; a < 3; ++a)
            if (shift[nb][size_t(a)] != sh[size_t(a)]) perc[size_t(a)] = true;
        }
      }
    }
    for (int a = 0; a < 3; ++a) any[size_t(a)] = any[size_t(a)] || perc[size_t(a)];
    const bool counts = axis < 0 ? (perc[0] || perc[1] || perc[2]) : perc[size_t(axis)];
    if (best && counts) *best = std::max(*best, top);
    ++L;
  }
  return any;
}

}  // namespace

PoreNetworkReport pore_network(const System& s, const PoreNetworkOptions& o) {
  if (!s.cell.valid()) throw std::invalid_argument("pores need a periodic cell");
  if (s.atoms.empty()) throw std::invalid_argument("pores: no atoms");
  if (!(o.probe >= 0)) throw std::invalid_argument("pores: the probe radius cannot be negative");
  PoreNetworkReport R;
  const Cell& c = s.cell;
  const size_t na = s.atoms.size();
  std::vector<double> rad(na);
  for (size_t i = 0; i < na; ++i) rad[i] = i < o.radii.size() ? o.radii[i] : element(s.atoms[i].element).vdw > 0 ? element(s.atoms[i].element).vdw : 1.7;
  // the distance field, reaching far enough for the largest included sphere
  VoidOptions vo;
  vo.grid = std::max(0.05, o.grid);
  vo.radii = rad;
  const double V = c.volume();
  const double widths[3] = {V / norm(cross(c.b, c.c)), V / norm(cross(c.c, c.a)), V / norm(cross(c.a, c.b))};
  vo.reach = 5.0;   // beyond this from every surface the distance is found exactly below
  R.field = void_distance_field(s, vo);
  {
    GridField& f = R.field;
    for (int i = 0; i < f.n[0]; ++i)
      for (int j = 0; j < f.n[1]; ++j)
        for (int k = 0; k < f.n[2]; ++k) {
          double& v = f.values[gidx(f.n, i, j, k)];
          if (v < vo.reach - 1e-5) continue;
          const Vec3 p = c.to_cartesian({(i + 0.5) / f.n[0], (j + 0.5) / f.n[1], (k + 0.5) / f.n[2]});
          double best = 1e300;
          for (size_t a = 0; a < na; ++a) best = std::min(best, norm(c.minimum_image(p - s.atoms[a].pos)) - rad[a]);
          v = best;
        }
  }
  const GridField& g = R.field;
  const size_t N = g.values.size();
  const double vox = V / double(N);
  double mass = 0;
  for (const auto& a : s.atoms) mass += element(a.element).mass;
  R.density = mass / 0.602214076 / V;   // g/mol per Å³ → g/cm³
  const double per_gram = 0.602214076 / mass;   // cm³/g per Å³ of the cell: N_A·10⁻²⁴ / M
  double fmax = -1e300;
  size_t geo = 0;
  for (double v : g.values) fmax = std::max(fmax, v), geo += v > 0;
  R.lcd = 2 * std::max(0.0, fmax);
  R.geometric_fraction = double(geo) / double(N);
  // clusters of probe-centre space
  const DomainResult D = label_domains(g, o.probe, true, 1);
  std::vector<char> channel_of(D.domains.size(), 0);
  std::vector<double> top(D.domains.size(), 0.0);
  for (size_t q = 0; q < N; ++q)
    if (D.label[q] >= 0) top[size_t(D.label[q])] = std::max(top[size_t(D.label[q])], g.values[q]);
  for (size_t k = 0; k < D.domains.size(); ++k) {
    PoreCluster pc;
    pc.voxels = D.domains[k].voxels;
    pc.volume = D.domains[k].volume;
    pc.percolates = D.domains[k].percolates;
    pc.channel = !D.domains[k].enclosed;
    pc.largest = 2 * top[k];
    channel_of[k] = pc.channel;
    (pc.channel ? R.av_fraction : R.nav_fraction) += pc.volume / V;
    R.clusters.push_back(pc);
  }
  std::array<bool, 3> dims{false, false, false};
  for (const auto& pc : R.clusters)
    if (pc.channel)
      for (int a = 0; a < 3; ++a) dims[size_t(a)] = dims[size_t(a)] || pc.percolates[size_t(a)];
  R.channel_dims = int(dims[0]) + int(dims[1]) + int(dims[2]);
  R.av = R.av_fraction * V * per_gram;
  R.nav = R.nav_fraction * V * per_gram;
  // pore-limiting diameter per axis: bisection on the radius at which the free space stops percolating
  for (int a = 0; a < 3; ++a) {
    double lo = -1e-12, hi = fmax;
    if (!percolation(g, 0.0, a, nullptr)[size_t(a)]) continue;   // not even a point passes along this axis
    lo = 0.0;
    for (int it = 0; it < 40 && hi - lo > 1e-7; ++it) {
      const double mid = 0.5 * (lo + hi);
      (percolation(g, mid, a, nullptr)[size_t(a)] ? lo : hi) = mid;
    }
    R.pld_axis[size_t(a)] = 2 * lo;
  }
  int best_axis = -1;
  for (int a = 0; a < 3; ++a)
    if (R.pld_axis[size_t(a)] > R.pld) R.pld = R.pld_axis[size_t(a)], best_axis = a;
  if (best_axis >= 0) {
    double d = 0;
    percolation(g, 0.5 * R.pld - 1e-6, best_axis, &d);
    R.dif = 2 * d;
  }
  // surface: Fibonacci points on every atom's sphere of radius R + probe, against the neighbours' spheres (a supercell
  // view by explicit images: the cell may be narrower than the spheres)
  const int M = std::max(50, o.samples);
  std::vector<Vec3> unit(static_cast<size_t>(M));
  const double golden = M_PI * (3 - std::sqrt(5.0));
  for (int k = 0; k < M; ++k) {
    const double z = 1 - (2.0 * k + 1) / M, r = std::sqrt(std::max(0.0, 1 - z * z)), phi = golden * k;
    unit[size_t(k)] = {r * std::cos(phi), r * std::sin(phi), z};
  }
  double rmax = 0;
  for (double x : rad) rmax = std::max(rmax, x + o.probe);
  int span[3];
  for (int a = 0; a < 3; ++a) span[a] = int(std::ceil(2 * rmax / widths[a]));
  // every centre wrapped into the home cell first: unwrapped or made-whole frames can put atoms several cells away,
  // and the ±span images below only reach the neighbours of atoms that sit in the cell
  std::vector<Vec3> pos(na);
  for (size_t i = 0; i < na; ++i) {
    Vec3 f = c.to_fractional(s.atoms[i].pos);
    for (int a = 0; a < 3; ++a) f[size_t(a)] -= std::floor(f[size_t(a)]);
    pos[i] = c.to_cartesian(f);
  }
  struct Img { Vec3 p; double r; };
  std::vector<std::vector<size_t>> near(na);
  std::vector<Img> imgs;
  for (size_t j = 0; j < na; ++j)
    for (int u = -span[0]; u <= span[0]; ++u)
      for (int v = -span[1]; v <= span[1]; ++v)
        for (int w = -span[2]; w <= span[2]; ++w) imgs.push_back({pos[j] + c.a * double(u) + c.b * double(v) + c.c * double(w), rad[j] + o.probe});
  // neighbours through a hash of bins 2 rmax wide over the images
  const double bin = 2 * rmax;
  auto key = [&](const Vec3& p, int dx, int dy, int dz) {
    const int64_t x = int64_t(std::floor(p[0] / bin)) + dx, y = int64_t(std::floor(p[1] / bin)) + dy, z = int64_t(std::floor(p[2] / bin)) + dz;
    return (x * 73856093) ^ (y * 19349663) ^ (z * 83492791);
  };
  std::unordered_map<int64_t, std::vector<size_t>> hash;
  for (size_t m = 0; m < imgs.size(); ++m) hash[key(imgs[m].p, 0, 0, 0)].push_back(m);
  for (size_t i = 0; i < na; ++i)
    for (int dx = -1; dx <= 1; ++dx)
      for (int dy = -1; dy <= 1; ++dy)
        for (int dz = -1; dz <= 1; ++dz) {
          const auto it = hash.find(key(pos[i], dx, dy, dz));
          if (it == hash.end()) continue;
          for (size_t m : it->second) {
            const double cut = rad[i] + o.probe + imgs[m].r;
            const Vec3 d = imgs[m].p - pos[i];
            if (dot(d, d) < cut * cut && dot(d, d) > 1e-12) near[i].push_back(m);
          }
        }
  for (auto& v : near) {   // hash collisions can list an image twice
    std::sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(), v.end()), v.end());
  }
  auto label_at = [&](const Vec3& p) {   // the cluster a probe centre at p opens into: the best labelled grid point around it
    Vec3 f = c.to_fractional(p);
    int b[3];
    for (int a = 0; a < 3; ++a) b[a] = int(std::floor((f[size_t(a)] - std::floor(f[size_t(a)])) * g.n[a]));
    int lab = -1;
    double best = -1e300;
    for (int x = -1; x <= 1; ++x)
      for (int y = -1; y <= 1; ++y)
        for (int z = -1; z <= 1; ++z) {
          const size_t q = gidx(g.n, ((b[0] + x) % g.n[0] + g.n[0]) % g.n[0], ((b[1] + y) % g.n[1] + g.n[1]) % g.n[1], ((b[2] + z) % g.n[2] + g.n[2]) % g.n[2]);
          if (D.label[q] >= 0 && g.values[q] > best) best = g.values[q], lab = D.label[q];
        }
    return lab;
  };
  size_t unassigned = 0;
  for (size_t i = 0; i < na; ++i) {
    const double Ri = rad[i] + o.probe;
    const double dA = 4 * M_PI * Ri * Ri / M;
    for (const auto& u : unit) {
      const Vec3 p = pos[i] + u * Ri;
      bool covered = false;
      for (size_t m : near[i]) {
        const Vec3 d = p - imgs[m].p;
        if (dot(d, d) < imgs[m].r * imgs[m].r) { covered = true; break; }
      }
      if (covered) continue;
      const int lab = label_at(p);
      if (lab < 0) { ++unassigned; continue; }
      (channel_of[size_t(lab)] ? R.asa : R.nasa) += dA;
    }
  }
  R.asa_volume = R.asa / V * 1e4, R.nasa_volume = R.nasa / V * 1e4;
  R.asa_mass = R.asa * 1e4 * per_gram, R.nasa_mass = R.nasa * 1e4 * per_gram;
  char b[240];
  std::snprintf(b, sizeof b, "probe %.2f Å, grid %d × %d × %d, %d surface points per atom; %zu channel and pocket cluster(s), channels along %d axis(es)", o.probe,
                g.n[0], g.n[1], g.n[2], M, R.clusters.size(), R.channel_dims);
  R.notes.push_back(b);
  if (unassigned) {
    std::snprintf(b, sizeof b, "%zu uncovered surface point(s) fell between grid points of no cluster (a finer grid places them)", unassigned);
    R.notes.push_back(b);
  }
  if (R.channel_dims == 0) R.notes.push_back("no channel: the probe cannot pass through the cell along any axis");
  return R;
}

Json pore_json(const PoreNetworkReport& r) {
  Json j = Json::object();
  j["lcd"] = r.lcd, j["pld"] = r.pld, j["dif"] = r.dif;
  Json pa = Json::array();
  for (double x : r.pld_axis) pa.push_back(x);
  j["pld_axis"] = pa;
  j["channel_dims"] = double(r.channel_dims), j["density"] = r.density;
  j["av_fraction"] = r.av_fraction, j["nav_fraction"] = r.nav_fraction, j["av"] = r.av, j["nav"] = r.nav, j["geometric_fraction"] = r.geometric_fraction;
  j["asa"] = r.asa, j["nasa"] = r.nasa, j["asa_volume"] = r.asa_volume, j["nasa_volume"] = r.nasa_volume, j["asa_mass"] = r.asa_mass, j["nasa_mass"] = r.nasa_mass;
  Json cl = Json::array();
  for (const auto& c : r.clusters) {
    Json x = Json::object();
    x["volume"] = c.volume, x["voxels"] = double(c.voxels), x["channel"] = c.channel, x["largest"] = c.largest;
    Json p = Json::array();
    for (bool b : c.percolates) p.push_back(b);
    x["percolates"] = p;
    cl.push_back(x);
  }
  j["clusters"] = cl;
  Json g = Json::array();
  for (int k = 0; k < 3; ++k) g.push_back(double(r.field.n[k]));
  j["grid"] = g;
  Json n = Json::array();
  for (const auto& x : r.notes) n.push_back(x);
  j["notes"] = n;
  return j;
}

}  // namespace caps
