// CAPS substitutional disorder (see sqs.hpp).
#include "caps/sqs.hpp"

#include "caps/io.hpp"

#include <algorithm>
#include <cmath>
#include <map>
#include <random>
#include <set>
#include <stdexcept>

namespace caps {

namespace {
std::vector<uint32_t> sublattice(const System& s, int z) {
  std::vector<uint32_t> v;
  for (uint32_t i = 0; i < s.atoms.size(); ++i)
    if (s.atoms[i].element == z) v.push_back(i);
  if (v.empty()) throw std::invalid_argument("disorder: no atoms of element " + std::to_string(z));
  return v;
}

// the cell's symmetry operations as permutations of the sublattice: integer rotations R of the fractional basis that keep
// the metric, with a translation that maps every atom onto an atom of its element
std::vector<std::vector<int>> site_permutations(const System& s, const std::vector<uint32_t>& sub, double tol) {
  const size_t N = s.atoms.size();
  std::vector<Vec3> f(N);
  for (size_t i = 0; i < N; ++i) {
    f[i] = s.cell.to_fractional(s.atoms[i].pos);
    for (auto& x : f[i]) x -= std::floor(x);
  }
  const Vec3 A[3] = {s.cell.a, s.cell.b, s.cell.c};
  double G[3][3];
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) G[i][j] = dot(A[i], A[j]);
  std::vector<int> subidx(N, -1);
  for (size_t k = 0; k < sub.size(); ++k) subidx[sub[k]] = int(k);
  std::vector<std::vector<int>> perms;
  std::set<std::vector<int>> seen;
  int R[3][3];
  for (int code = 0; code < 19683; ++code) {
    int c = code;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) R[i][j] = c % 3 - 1, c /= 3;
    const int det = R[0][0] * (R[1][1] * R[2][2] - R[1][2] * R[2][1]) - R[0][1] * (R[1][0] * R[2][2] - R[1][2] * R[2][0]) + R[0][2] * (R[1][0] * R[2][1] - R[1][1] * R[2][0]);
    if (det != 1 && det != -1) continue;
    bool metric = true;   // Rᵀ G R = G
    for (int i = 0; i < 3 && metric; ++i)
      for (int j = 0; j < 3 && metric; ++j) {
        double v = 0;
        for (int a = 0; a < 3; ++a)
          for (int b = 0; b < 3; ++b) v += R[a][i] * G[a][b] * R[b][j];
        metric = std::fabs(v - G[i][j]) < 1e-3 * (1 + std::fabs(G[i][j]));
      }
    if (!metric) continue;
    auto rot = [&](const Vec3& x) { return Vec3{R[0][0] * x[0] + R[0][1] * x[1] + R[0][2] * x[2], R[1][0] * x[0] + R[1][1] * x[1] + R[1][2] * x[2], R[2][0] * x[0] + R[2][1] * x[1] + R[2][2] * x[2]}; };
    const Vec3 r0 = rot(f[sub[0]]);
    for (uint32_t tgt : sub) {
      const Vec3 t = f[tgt] - r0;
      std::vector<int> perm(sub.size(), -1);
      bool ok = true;
      for (size_t i = 0; i < N && ok; ++i) {
        Vec3 y = rot(f[i]) + t;
        bool hit = false;
        for (size_t j = 0; j < N; ++j) {
          if (s.atoms[j].element != s.atoms[i].element) continue;
          Vec3 d = y - f[j];
          for (auto& x : d) x -= std::round(x);
          if (norm(s.cell.to_cartesian(d)) < tol) {
            if (subidx[i] >= 0) perm[size_t(subidx[i])] = subidx[j];
            hit = true;
            break;
          }
        }
        ok = hit;
      }
      if (ok && seen.insert(perm).second) perms.push_back(perm);
    }
  }
  return perms;
}
}  // namespace

SqsResult special_quasirandom(const System& s, const SqsOptions& o) {
  SqsResult r;
  if (!s.cell.valid()) throw std::invalid_argument("SQS: a periodic cell");
  if (o.to <= 0) throw std::invalid_argument("SQS: the substituting element");
  const auto sub = sublattice(s, o.from);
  const size_t n = sub.size();
  const int k = int(std::lround(o.fraction * double(n)));
  if (k <= 0 || k >= int(n)) throw std::invalid_argument("SQS: the fraction leaves no disorder on " + std::to_string(n) + " sites");
  r.sites = int(n), r.substituted = k;
  const double x = double(k) / double(n);
  r.target = (1 - 2 * x) * (1 - 2 * x);
  // shells: distinct minimum-image distances between sublattice sites
  std::vector<double> dist;
  std::vector<std::vector<double>> D(n, std::vector<double>(n, 0));
  for (size_t i = 0; i < n; ++i)
    for (size_t j = i + 1; j < n; ++j) {
      D[i][j] = D[j][i] = norm(s.cell.minimum_image(s.atoms[sub[j]].pos - s.atoms[sub[i]].pos));
      dist.push_back(D[i][j]);
    }
  std::sort(dist.begin(), dist.end());
  for (double d : dist)
    if (r.shell_distance.empty() || d - r.shell_distance.back() > 1e-3) r.shell_distance.push_back(d);
  if (r.shell_distance.size() > size_t(o.shells)) r.shell_distance.resize(size_t(std::max(1, o.shells)));
  const size_t M = r.shell_distance.size();
  // per site, per shell: its neighbours
  std::vector<std::vector<std::vector<int>>> nb(n, std::vector<std::vector<int>>(M));
  std::vector<double> npair(M, 0);
  for (size_t i = 0; i < n; ++i)
    for (size_t j = 0; j < n; ++j) {
      if (i == j) continue;
      for (size_t m = 0; m < M; ++m)
        if (std::fabs(D[i][j] - r.shell_distance[m]) < 1e-3) {
          nb[i][m].push_back(int(j));
          if (i < j) npair[m] += 1;
        }
    }
  std::mt19937_64 rng(o.seed);
  std::uniform_real_distribution<double> U(0, 1);
  auto sums = [&](const std::vector<int>& sg) {
    std::vector<double> S(M, 0);
    for (size_t i = 0; i < n; ++i)
      for (size_t m = 0; m < M; ++m)
        for (int j : nb[i][m]) if (size_t(j) > i) S[m] += sg[i] * sg[size_t(j)];
    return S;
  };
  auto objective = [&](const std::vector<double>& S) {
    double f = 0;
    for (size_t m = 0; m < M; ++m) f += std::fabs(S[m] / npair[m] - r.target) / double(m + 1);
    return f;
  };
  std::vector<int> best;
  double bestf = 1e300;
  for (int rs = 0; rs < std::max(1, o.restarts); ++rs) {
    std::vector<int> sg(n, 1);
    std::vector<size_t> idx(n);
    for (size_t i = 0; i < n; ++i) idx[i] = i;
    std::shuffle(idx.begin(), idx.end(), rng);
    for (int q = 0; q < k; ++q) sg[idx[size_t(q)]] = -1;
    std::vector<double> S = sums(sg);
    double f = objective(S);
    for (int it = 0; it < o.steps && f > 1e-12; ++it) {
      const double T = 0.1 * std::pow(1e-4, double(it) / std::max(1, o.steps - 1));
      size_t a = size_t(U(rng) * n) % n, b = size_t(U(rng) * n) % n;
      if (sg[a] == sg[b]) continue;
      std::vector<double> S2 = S;
      for (size_t m = 0; m < M; ++m) {
        for (int j : nb[a][m]) if (size_t(j) != b) S2[m] += -2.0 * sg[a] * sg[size_t(j)];
        for (int j : nb[b][m]) if (size_t(j) != a) S2[m] += -2.0 * sg[b] * sg[size_t(j)];
      }
      const double f2 = objective(S2);
      if (f2 <= f || U(rng) < std::exp(-(f2 - f) / T)) {
        std::swap(sg[a], sg[b]), S = S2, f = f2;
        if (f < bestf) bestf = f, best = sg;
      }
    }
    if (f < bestf) bestf = f, best = sg;
  }
  r.objective = bestf;
  const std::vector<double> S = sums(best);
  for (size_t m = 0; m < M; ++m) r.correlation.push_back(S[m] / npair[m]);
  r.structure = s;
  for (size_t i = 0; i < n; ++i)
    if (best[i] < 0) {
      Atom& at = r.structure.atoms[sub[i]];
      at.element = o.to;
      at.name.clear();
    }
  r.notes.push_back(std::to_string(k) + " of " + std::to_string(n) + " sites substituted; " + std::to_string(M) + " pair shells matched to (1 − 2x)² = " +
                    std::to_string(r.target));
  return r;
}

EnumResult enumerate_occupations(const System& s, const EnumOptions& o) {
  EnumResult r;
  if (!s.cell.valid()) throw std::invalid_argument("enumerate: a periodic cell");
  r.sites = sublattice(s, o.from);
  const size_t n = r.sites.size();
  const int k = o.count;
  if (k <= 0 || k > int(n)) throw std::invalid_argument("enumerate: the count must be between 1 and the sites");
  double total = 1;
  for (int i = 0; i < k; ++i) total = total * double(n - size_t(i)) / double(i + 1);
  if (total > double(o.max_configurations)) throw std::invalid_argument("enumerate: C(n, k) = " + std::to_string(size_t(total)) + " occupations is too many");
  r.total = size_t(std::llround(total));
  const auto perms = site_permutations(s, r.sites, o.tolerance);
  r.operations = int(perms.size());
  std::map<std::vector<int>, size_t> orbit;
  std::vector<int> comb(static_cast<size_t>(k), 0);
  for (int i = 0; i < k; ++i) comb[size_t(i)] = i;
  while (true) {
    std::vector<int> best = comb;
    for (const auto& p : perms) {
      std::vector<int> img(static_cast<size_t>(k), 0);
      for (int i = 0; i < k; ++i) img[size_t(i)] = p[size_t(comb[size_t(i)])];
      std::sort(img.begin(), img.end());
      if (img < best) best = img;
    }
    ++orbit[best];
    int i = k - 1;
    while (i >= 0 && comb[size_t(i)] == int(n) - k + i) --i;
    if (i < 0) break;
    ++comb[size_t(i)];
    for (int j = i + 1; j < k; ++j) comb[size_t(j)] = comb[size_t(j - 1)] + 1;
  }
  for (const auto& [c, m] : orbit) r.configurations.push_back(c), r.multiplicity.push_back(m);
  r.notes.push_back(std::to_string(r.total) + " occupations, " + std::to_string(r.configurations.size()) + " distinct under " + std::to_string(r.operations) + " operations");
  return r;
}

}  // namespace caps

namespace caps {

System occupied(const System& cell, const EnumResult& r, size_t which, int to) {
  if (which >= r.configurations.size()) throw std::out_of_range("enumerate: no configuration " + std::to_string(which + 1));
  System s = cell;
  for (int k : r.configurations[which]) {
    Atom& a = s.atoms[r.sites[size_t(k)]];
    a.element = to;
    a.name.clear();
  }
  return s;
}

void write_disordered(const System& s, const std::string& path) {
  auto ends = [&](const std::string& x) { return path.size() >= x.size() && path.compare(path.size() - x.size(), x.size(), x) == 0; };
  const std::string base = path.substr(path.find_last_of("/\\") == std::string::npos ? 0 : path.find_last_of("/\\") + 1);
  if (ends(".cif")) write_cif(s, path);
  else if (ends(".vasp") || base.rfind("POSCAR", 0) == 0) write_poscar(s, path);
  else write_pdb(s, path);
}

Json sqs_json(const SqsResult& r) {
  Json j = Json::object(), d = Json::array(), c = Json::array(), n = Json::array();
  for (double x : r.shell_distance) d.push_back(x);
  for (double x : r.correlation) c.push_back(x);
  for (const auto& x : r.notes) n.push_back(x);
  j["shell_distance"] = d, j["correlation"] = c, j["target"] = r.target, j["objective"] = r.objective;
  j["sites"] = double(r.sites), j["substituted"] = double(r.substituted), j["notes"] = n;
  return j;
}

Json enumeration_json(const EnumResult& r) {
  Json j = Json::object(), cf = Json::array(), m = Json::array(), st = Json::array(), n = Json::array();
  for (const auto& c : r.configurations) {
    Json a = Json::array();
    for (int k : c) a.push_back(double(k));
    cf.push_back(a);
  }
  for (size_t x : r.multiplicity) m.push_back(double(x));
  for (uint32_t x : r.sites) st.push_back(double(x));
  for (const auto& x : r.notes) n.push_back(x);
  j["configurations"] = cf, j["multiplicity"] = m, j["sites"] = st, j["operations"] = double(r.operations), j["total"] = double(r.total), j["notes"] = n;
  return j;
}

}  // namespace caps
