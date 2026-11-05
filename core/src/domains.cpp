// CAPS domains of a density field (see domains.hpp).
#include "caps/domains.hpp"

#include <algorithm>
#include <cmath>
#include <deque>
#include <stdexcept>

#include "caps/analysis.hpp"
#include "caps/query.hpp"

namespace caps {

namespace {
constexpr double kPi = 3.14159265358979323846;
size_t idx(const int n[3], int i, int j, int k) { return (size_t(i) * size_t(n[1]) + size_t(j)) * size_t(n[2]) + size_t(k); }
int wrapi(int x, int n) { return ((x % n) + n) % n; }
}  // namespace

GridField density_field(const System& s, const std::vector<char>& atoms, double grid, double sigma) {
  if (!s.cell.valid()) throw std::invalid_argument("a density field needs a periodic cell");
  if (!(grid > 0) || !(sigma > 0)) throw std::invalid_argument("density field: a positive grid and width");
  GridField g;
  g.name = "mass density", g.unit = "g/cm³", g.cell = s.cell;
  const Cell& c = s.cell;
  const double V = c.volume();
  const Vec3 A[3] = {c.a, c.b, c.c};
  double width[3];   // the cell's width perpendicular to each pair of the other edges
  for (int k = 0; k < 3; ++k) width[k] = V / norm(cross(A[(k + 1) % 3], A[(k + 2) % 3]));
  for (int k = 0; k < 3; ++k) g.n[k] = std::max(4, int(std::ceil(norm(A[k]) / grid)));
  const size_t N = size_t(g.n[0]) * size_t(g.n[1]) * size_t(g.n[2]);
  g.values.assign(N, 0.0);
  const double vox = V / double(N);
  int reach[3];
  for (int k = 0; k < 3; ++k) reach[k] = int(std::ceil(3 * sigma / (width[k] / g.n[k]))) + 1;
  for (size_t a = 0; a < s.atoms.size(); ++a) {
    if (!atoms.empty() && (a >= atoms.size() || !atoms[a])) continue;
    const double m = s.mass_of(s.atoms[a]) * 1.66053906660;   // g/mol per Å³ → g/cm³
    Vec3 f = c.to_fractional(s.atoms[a].pos);
    int c0[3];
    for (int k = 0; k < 3; ++k) c0[k] = int(std::floor((f[size_t(k)] - std::floor(f[size_t(k)])) * g.n[k]));
    std::vector<std::pair<size_t, double>> ker;
    double sum = 0;
    for (int di = -reach[0]; di <= reach[0]; ++di)
      for (int dj = -reach[1]; dj <= reach[1]; ++dj)
        for (int dk = -reach[2]; dk <= reach[2]; ++dk) {
          const int i = wrapi(c0[0] + di, g.n[0]), j = wrapi(c0[1] + dj, g.n[1]), k = wrapi(c0[2] + dk, g.n[2]);
          const Vec3 p = c.origin + c.a * ((i + 0.5) / g.n[0]) + c.b * ((j + 0.5) / g.n[1]) + c.c * ((k + 0.5) / g.n[2]);
          const double d2 = dot(c.minimum_image(p - s.atoms[a].pos), c.minimum_image(p - s.atoms[a].pos));
          if (d2 > 9 * sigma * sigma) continue;
          const double w = std::exp(-d2 / (2 * sigma * sigma));
          ker.push_back({idx(g.n, i, j, k), w});
          sum += w;
        }
    if (sum <= 0) continue;
    for (const auto& [q, w] : ker) g.values[q] += m * w / (sum * vox);
  }
  return g;
}

DomainResult label_domains(const GridField& g, double threshold, bool above, int min_voxels) {
  DomainResult R;
  const size_t N = size_t(g.n[0]) * size_t(g.n[1]) * size_t(g.n[2]);
  if (g.values.size() != N || N == 0) throw std::invalid_argument("domains: an empty or malformed grid");
  R.field = g;
  R.threshold = threshold;
  auto in = [&](size_t q) { return above ? g.values[q] > threshold : g.values[q] < threshold; };
  std::vector<int> lab(N, -1);
  std::vector<std::array<int, 3>> shift(N, {0, 0, 0});
  const Cell& c = g.cell;
  const double vox = c.volume() / double(N);
  struct Raw { std::vector<size_t> vox; std::array<bool, 3> perc{false, false, false}; };
  std::vector<Raw> raw;
  for (int i0 = 0; i0 < g.n[0]; ++i0)
    for (int j0 = 0; j0 < g.n[1]; ++j0)
      for (int k0 = 0; k0 < g.n[2]; ++k0) {
        const size_t s0 = idx(g.n, i0, j0, k0);
        if (lab[s0] >= 0 || !in(s0)) continue;
        const int L = int(raw.size());
        raw.push_back({});
        Raw& r = raw.back();
        std::deque<std::array<int, 3>> q{{i0, j0, k0}};
        lab[s0] = L, shift[s0] = {0, 0, 0};
        while (!q.empty()) {
          const auto [i, j, k] = q.front();
          q.pop_front();
          const size_t cur = idx(g.n, i, j, k);
          r.vox.push_back(cur);
          const int d[6][3] = {{1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1}};
          for (const auto& dd : d) {
            int x[3] = {i + dd[0], j + dd[1], k + dd[2]};
            std::array<int, 3> sh = shift[cur];
            for (int a = 0; a < 3; ++a) {
              if (x[a] < 0) x[a] += g.n[a], sh[size_t(a)] -= 1;
              else if (x[a] >= g.n[a]) x[a] -= g.n[a], sh[size_t(a)] += 1;
            }
            const size_t nb = idx(g.n, x[0], x[1], x[2]);
            if (!in(nb)) continue;
            if (lab[nb] < 0) {
              lab[nb] = L, shift[nb] = sh;
              q.push_back({x[0], x[1], x[2]});
            } else if (lab[nb] == L) {
              for (int a = 0; a < 3; ++a)
                if (shift[nb][size_t(a)] != sh[size_t(a)]) r.perc[size_t(a)] = true;
            }
          }
        }
      }
  // the regions' numbers
  std::vector<std::pair<size_t, int>> order;
  for (size_t L = 0; L < raw.size(); ++L)
    if (int(raw[L].vox.size()) >= min_voxels) order.push_back({raw[L].vox.size(), int(L)});
  std::sort(order.begin(), order.end(), [](const auto& a, const auto& b) { return a.first > b.first; });
  std::vector<int> newl(raw.size(), -1);
  for (size_t k = 0; k < order.size(); ++k) newl[size_t(order[k].second)] = int(k);
  R.label.assign(N, -1);
  for (size_t q = 0; q < N; ++q)
    if (lab[q] >= 0) R.label[q] = newl[size_t(lab[q])];
  size_t total = 0;
  for (const auto& [sz, L] : order) {
    const Raw& r = raw[size_t(L)];
    Domain D;
    D.voxels = r.vox.size();
    D.volume = double(D.voxels) * vox;
    D.equivalent_radius = std::cbrt(3 * D.volume / (4 * kPi));
    D.percolates = r.perc;
    D.enclosed = !(r.perc[0] || r.perc[1] || r.perc[2]);
    // unwrapped voxel centres: the centroid and the gyration tensor
    std::vector<Vec3> pts;
    pts.reserve(r.vox.size());
    for (size_t q : r.vox) {
      const int i = int(q / (size_t(g.n[1]) * size_t(g.n[2]))), j = int((q / size_t(g.n[2])) % size_t(g.n[1])), k = int(q % size_t(g.n[2]));
      const auto& sh = shift[q];
      pts.push_back(c.origin + c.a * ((i + 0.5) / g.n[0] + sh[0]) + c.b * ((j + 0.5) / g.n[1] + sh[1]) + c.c * ((k + 0.5) / g.n[2] + sh[2]));
    }
    Vec3 m{0, 0, 0};
    for (const auto& p : pts) m = m + p;
    m = m * (1.0 / double(pts.size()));
    double G[3][3] = {};
    for (const auto& p : pts) {
      const Vec3 d = p - m;
      for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b) G[a][b] += d[size_t(a)] * d[size_t(b)] / double(pts.size());
    }
    double w[3], Vv[3][3];
    symmetric_eigen3(G, w, Vv);
    D.axes = {std::sqrt(5 * std::max(0.0, w[2])), std::sqrt(5 * std::max(0.0, w[1])), std::sqrt(5 * std::max(0.0, w[0]))};
    D.centroid = c.wrap(m);
    R.domains.push_back(D);
    total += D.voxels;
  }
  R.fraction = double(total) / double(N);
  // nearest enclosed neighbours
  for (size_t a = 0; a < R.domains.size(); ++a) {
    if (!R.domains[a].enclosed) continue;
    double best = 0;
    for (size_t b = 0; b < R.domains.size(); ++b) {
      if (a == b || !R.domains[b].enclosed) continue;
      const double d = norm(c.minimum_image(R.domains[b].centroid - R.domains[a].centroid));
      if (best == 0 || d < best) best = d;
    }
    R.domains[a].nearest = best;
  }
  // the voxel values
  double lo = 1e300, hi = -1e300;
  for (double v : g.values) lo = std::min(lo, v), hi = std::max(hi, v);
  const int nb = 50;
  const double bw = hi > lo ? (hi - lo) / nb : 1.0;
  R.hist_y.assign(nb, 0.0);
  for (int k = 0; k < nb; ++k) R.hist_x.push_back(lo + (k + 0.5) * bw);
  for (double v : g.values) R.hist_y[size_t(std::clamp(int((v - lo) / bw), 0, nb - 1))] += 1.0 / double(N);
  int perc = 0;
  for (const auto& D : R.domains) perc += !D.enclosed;
  R.notes.push_back(std::to_string(R.domains.size()) + " domains (" + std::to_string(perc) + " percolating) covering " + std::to_string(int(std::lround(100 * R.fraction))) + " % of the cell");
  return R;
}

DomainResult density_domains(const System& s, const DomainOptions& o) {
  std::vector<char> atoms;
  if (!o.group.empty()) {
    atoms = select_query(s, o.group).atoms;
    if (std::none_of(atoms.begin(), atoms.end(), [](char x) { return x != 0; })) throw std::invalid_argument("domains: the group selects no atoms");
  }
  const GridField g = density_field(s, atoms, o.grid, o.sigma);
  double mean = 0;
  for (double v : g.values) mean += v;
  mean /= double(g.values.size());
  DomainResult R = label_domains(g, o.threshold > 0 ? o.threshold : mean, o.above, o.min_voxels);
  if (!(o.threshold > 0)) R.notes.push_back("threshold: the field's mean, " + std::to_string(mean).substr(0, 6) + " " + g.unit);
  return R;
}

std::vector<double> field_slice(const GridField& g, const Vec3& point, const Vec3& u, const Vec3& v, double width, double height, int nx, int ny) {
  if (nx < 1 || ny < 1) throw std::invalid_argument("slice: at least one sample each way");
  const Vec3 uu = u * (1.0 / norm(u)), vv = v * (1.0 / norm(v));
  std::vector<double> out;
  out.reserve(size_t(nx) * size_t(ny));
  auto at = [&](int i, int j, int k) { return g.values[idx(g.n, wrapi(i, g.n[0]), wrapi(j, g.n[1]), wrapi(k, g.n[2]))]; };
  for (int y = 0; y < ny; ++y)
    for (int x = 0; x < nx; ++x) {
      const double sx = nx > 1 ? width * (double(x) / (nx - 1) - 0.5) : 0, sy = ny > 1 ? height * (double(y) / (ny - 1) - 0.5) : 0;
      const Vec3 p = point + uu * sx + vv * sy;
      const Vec3 f = g.cell.to_fractional(p);
      double t[3];
      int i0[3];
      for (int k = 0; k < 3; ++k) {
        const double q = f[size_t(k)] * g.n[k] - 0.5;
        i0[k] = int(std::floor(q));
        t[k] = q - i0[k];
      }
      double val = 0;
      for (int a = 0; a < 2; ++a)
        for (int b = 0; b < 2; ++b)
          for (int cc = 0; cc < 2; ++cc)
            val += (a ? t[0] : 1 - t[0]) * (b ? t[1] : 1 - t[1]) * (cc ? t[2] : 1 - t[2]) * at(i0[0] + a, i0[1] + b, i0[2] + cc);
      out.push_back(val);
    }
  return out;
}

}  // namespace caps
