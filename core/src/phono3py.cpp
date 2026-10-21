// Lattice thermal conductivity with VASP and phono3py: plan, inputs and step scripts, checks and analysis (caps/phono3py.hpp).
#include "caps/phono3py.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

#include "caps/elements.hpp"
#include "caps/slab2d.hpp"
#include "caps/vasp_out.hpp"
#include "caps/vasp_set.hpp"

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
void spit(const std::string& p, const std::string& t) {
  std::ofstream f(p, std::ios::binary);
  if (!f) throw std::runtime_error("cannot write " + p);
  f << t;
}
void make_exec(const std::string& p) {
  std::error_code ec;
  fs::permissions(p, fs::perms::owner_exec | fs::perms::group_exec | fs::perms::others_exec, fs::perm_options::add, ec);
}
std::string join(const std::string& a, const std::string& b) { return (fs::path(a) / b).string(); }
std::string fmt(const char* f, double v) { char b[96]; std::snprintf(b, sizeof b, f, v); return b; }
std::string dims(const std::array<int, 3>& d) { return std::to_string(d[0]) + " " + std::to_string(d[1]) + " " + std::to_string(d[2]); }
Json arr3(const std::array<int, 3>& d) { Json a = Json::array(); for (int x : d) a.push_back(x); return a; }
std::vector<std::string> species_of(const System& s) {
  std::vector<std::string> out;
  for (const auto& a : s.atoms) { const std::string e = element(a.element).symbol; if (out.empty() || out.back() != e) out.push_back(e); }
  return out;
}
Json read_json_file(const std::string& path) {
  const std::string t = slurp(path);
  if (t.empty()) throw std::runtime_error("cannot read " + path);
  return Json::parse(t);
}
double heights(const System& s, int k) {
  const Vec3 v[3] = {s.cell.a, s.cell.b, s.cell.c};
  const Vec3 n = cross(v[(k + 1) % 3], v[(k + 2) % 3]);
  return std::fabs(dot(v[k], n)) / norm(n);
}
using IMat = std::array<std::array<int, 3>, 3>;
using IVec = std::array<int, 3>;
IVec mul(const IMat& R, const IVec& d) {
  IVec o{0, 0, 0};
  for (int i = 0; i < 3; ++i) for (int k = 0; k < 3; ++k) o[size_t(i)] += R[size_t(i)][size_t(k)] * d[size_t(k)];
  return o;
}
long det3(const IVec& a, const IVec& b, const IVec& c) {
  return long(a[0]) * (b[1] * c[2] - b[2] * c[1]) - long(a[1]) * (b[0] * c[2] - b[2] * c[0]) + long(a[2]) * (b[0] * c[1] - b[1] * c[0]);
}
bool covalent(int zi, int zj, double r) { return r < 1.2 * (element(zi).covalent + element(zj).covalent) && r > 0.4; }

}  // namespace

// ---------------------------------------------------------------- structure
std::vector<int> chain_axes(const System& s) {
  const size_t n = s.atoms.size();
  std::vector<std::vector<std::pair<size_t, IVec>>> adj(n);
  for (const auto& nb : neighbours_pbc(s, 2.6)) {
    if (!covalent(s.atoms[nb.i].element, s.atoms[nb.j].element, nb.r)) continue;
    // the lattice image of j that bond reaches: f_i + d − f_j
    const Vec3 fi = s.cell.to_fractional(s.atoms[nb.i].pos), fj = s.cell.to_fractional(s.atoms[nb.j].pos);
    const Vec3 fd = s.cell.to_fractional(s.atoms[nb.i].pos + nb.d) - fi;
    IVec im;
    for (int k = 0; k < 3; ++k) im[size_t(k)] = int(std::lround(fi[size_t(k)] + fd[size_t(k)] - fj[size_t(k)]));
    adj[nb.i].push_back({nb.j, im});
  }
  std::vector<IVec> off(n);
  std::vector<char> seen(n, 0);
  std::set<int> axes;
  for (size_t s0 = 0; s0 < n; ++s0) {
    if (seen[s0]) continue;
    seen[s0] = 1, off[s0] = {0, 0, 0};
    std::vector<size_t> st{s0};
    while (!st.empty()) {
      const size_t i = st.back();
      st.pop_back();
      for (const auto& [j, im] : adj[i]) {
        const IVec want{off[i][0] + im[0], off[i][1] + im[1], off[i][2] + im[2]};
        if (!seen[j]) { seen[j] = 1, off[j] = want, st.push_back(j); continue; }
        for (int k = 0; k < 3; ++k) if (want[size_t(k)] != off[j][size_t(k)]) axes.insert(k);   // a cycle through the boundary
      }
    }
  }
  return {axes.begin(), axes.end()};
}

std::array<int, 3> supercell_for_length(const System& s, double length, double slack) {
  std::array<int, 3> n{};
  for (int k = 0; k < 3; ++k) n[size_t(k)] = std::max(1, int(std::ceil(length * (1 - slack) / heights(s, k) - 1e-9)));
  return n;
}

std::vector<SymOpFrac> crystal_operations(const System& s, double tol) {
  const Vec3 L[3] = {s.cell.a, s.cell.b, s.cell.c};
  double G[3][3], gmax = 0;
  for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) G[i][j] = dot(L[i], L[j]), gmax = std::max(gmax, std::fabs(G[i][j]));
  std::vector<Vec3> f;
  for (const auto& a : s.atoms) { Vec3 x = s.cell.to_fractional(a.pos); for (auto& v : x) v -= std::floor(v); f.push_back(x); }
  auto cart = [&](const Vec3& d) { return L[0] * d[0] + L[1] * d[1] + L[2] * d[2]; };
  auto close = [&](Vec3 d) { for (auto& v : d) v -= std::round(v); return norm(cart(d)) < tol; };
  std::vector<SymOpFrac> out;
  if (s.atoms.empty()) return out;
  // the reference: an atom of the rarest element (fewest translations to try)
  std::map<int, int> count;
  for (const auto& a : s.atoms) ++count[a.element];
  size_t ref = 0;
  for (size_t i = 0; i < s.atoms.size(); ++i) if (count[s.atoms[i].element] < count[s.atoms[ref].element]) ref = i;
  IMat R;
  for (int code = 0; code < 19683; ++code) {
    int c = code;
    for (int i = 0; i < 3; ++i) for (int j = 0; j < 3; ++j) R[size_t(i)][size_t(j)] = c % 3 - 1, c /= 3;
    const long d = det3(R[0], R[1], R[2]);   // rows; |det| is the same as for columns
    if (d != 1 && d != -1) continue;
    bool metric = true;   // Rᵀ G R = G
    for (int i = 0; i < 3 && metric; ++i)
      for (int j = 0; j < 3 && metric; ++j) {
        double x = 0;
        for (int k = 0; k < 3; ++k) for (int l = 0; l < 3; ++l) x += R[size_t(k)][size_t(i)] * G[k][l] * R[size_t(l)][size_t(j)];
        metric = std::fabs(x - G[i][j]) < 1e-3 * gmax;
      }
    if (!metric) continue;
    auto rot = [&](const Vec3& x) { Vec3 y{0, 0, 0}; for (int i = 0; i < 3; ++i) for (int k = 0; k < 3; ++k) y[size_t(i)] += R[size_t(i)][size_t(k)] * x[size_t(k)]; return y; };
    const Vec3 r0 = rot(f[ref]);
    for (size_t j = 0; j < s.atoms.size(); ++j) {
      if (s.atoms[j].element != s.atoms[ref].element) continue;
      Vec3 t = f[j] - r0;
      for (auto& v : t) v -= std::floor(v + 1e-9);
      bool ok = true;
      for (size_t i = 0; i < s.atoms.size() && ok; ++i) {
        const Vec3 y = rot(f[i]) + t;
        bool hit = false;
        for (size_t k = 0; k < s.atoms.size() && !hit; ++k) hit = s.atoms[k].element == s.atoms[i].element && close(y - f[k]);
        ok = hit;
      }
      if (!ok) continue;
      bool dup = false;
      for (const auto& o : out) if (o.R == R && close(o.t - t)) dup = true;
      if (!dup) out.push_back({R, t});
    }
  }
  return out;
}

// ---------------------------------------------------------------- phono3py's displacement count
namespace {

struct Super {
  std::vector<Vec3> f;              // fractional in the supercell
  std::vector<int> z;
  Vec3 L[3];
  std::vector<std::pair<IMat, Vec3>> ops;
  double tol = 1e-3;
  std::unordered_map<long long, std::vector<int>> grid;
  int g = 64;
  long long key(const Vec3& x) const {
    long long k = 0;
    for (int i = 0; i < 3; ++i) { double v = x[size_t(i)] - std::floor(x[size_t(i)]); int c = int(v * g) % g; k = k * g + c; }
    return k;
  }
  void index() { for (size_t i = 0; i < f.size(); ++i) grid[key(f[i])].push_back(int(i)); }
  int find(Vec3 x, int elem) const {
    for (auto& v : x) v -= std::floor(v);
    // the neighbouring grid cells too
    for (int dx = -1; dx <= 1; ++dx)
      for (int dy = -1; dy <= 1; ++dy)
        for (int dz = -1; dz <= 1; ++dz) {
          const Vec3 y{x[0] + double(dx) / g, x[1] + double(dy) / g, x[2] + double(dz) / g};
          const auto it = grid.find(key(y));
          if (it == grid.end()) continue;
          for (int k : it->second) {
            if (z[size_t(k)] != elem) continue;
            Vec3 d = x - f[size_t(k)];
            for (auto& v : d) v -= std::round(v);
            if (norm(L[0] * d[0] + L[1] * d[1] + L[2] * d[2]) < tol) return k;
          }
        }
    return -1;
  }
  Vec3 apply(const std::pair<IMat, Vec3>& op, const Vec3& x) const {
    Vec3 y = op.second;
    for (int i = 0; i < 3; ++i) for (int k = 0; k < 3; ++k) y[size_t(i)] += op.first[size_t(i)][size_t(k)] * x[size_t(k)];
    return y;
  }
};

Super make_super(const System& unit, const std::vector<SymOpFrac>& uops, const std::array<int, 3>& d, double tol) {
  Super S;
  S.tol = tol;
  S.L[0] = unit.cell.a * d[0], S.L[1] = unit.cell.b * d[1], S.L[2] = unit.cell.c * d[2];
  for (int i = 0; i < d[0]; ++i)
    for (int j = 0; j < d[1]; ++j)
      for (int k = 0; k < d[2]; ++k)
        for (const auto& a : unit.atoms) {
          Vec3 f = unit.cell.to_fractional(a.pos);
          for (auto& v : f) v -= std::floor(v);
          S.f.push_back({(f[0] + i) / d[0], (f[1] + j) / d[1], (f[2] + k) / d[2]});
          S.z.push_back(a.element);
        }
  S.index();
  // supercell operations: R_s = D⁻¹ R D (integer, else not a symmetry of the supercell lattice), t_s = D⁻¹ (t + n)
  for (const auto& op : uops) {
    IMat Rs;
    bool integer = true;
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) {
        const int num = op.R[size_t(i)][size_t(j)] * d[size_t(j)];
        if (num % d[size_t(i)] != 0) integer = false;
        Rs[size_t(i)][size_t(j)] = num / d[size_t(i)];
      }
    if (!integer) continue;
    for (int i = 0; i < d[0]; ++i)
      for (int j = 0; j < d[1]; ++j)
        for (int k = 0; k < d[2]; ++k) S.ops.push_back({Rs, Vec3{(op.t[0] + i) / d[0], (op.t[1] + j) / d[1], (op.t[2] + k) / d[2]}});
  }
  return S;
}

// phonopy's least displacements for a site symmetry (rotations in lattice coordinates): one direction whose images span
// 3D, else two (the first, one of its images and a second direction), else the three axes; diagonal candidates or axes
const std::vector<IVec>& axis_dirs() { static const std::vector<IVec> d = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}; return d; }
const std::vector<IVec>& diag_dirs() {
  static const std::vector<IVec> d = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {1, 1, 0}, {1, 0, 1}, {0, 1, 1}, {1, -1, 0}, {1, 0, -1}, {0, 1, -1},
                                      {1, 1, 1}, {1, 1, -1}, {1, -1, 1}, {-1, 1, 1}};
  return d;
}
std::vector<IVec> least_dirs(const std::vector<IMat>& site, bool diagonal) {
  const auto& dirs = diagonal ? diag_dirs() : axis_dirs();
  std::vector<IVec> chosen;
  for (const auto& d : dirs) {
    std::vector<IVec> im;
    for (const auto& R : site) im.push_back(mul(R, d));
    for (size_t i = 0; i < im.size() && chosen.empty(); ++i)
      for (size_t j = i + 1; j < im.size() && chosen.empty(); ++j)
        if (det3(d, im[i], im[j]) != 0) chosen = {d};
    if (!chosen.empty()) break;
  }
  if (chosen.empty())
    for (const auto& d1 : dirs) {
      for (const auto& R : site) {
        const IVec r = mul(R, d1);
        for (const auto& d2 : dirs) if (det3(d1, r, d2) != 0) { chosen = {d1, d2}; break; }
        if (!chosen.empty()) break;
      }
      if (!chosen.empty()) break;
    }
  if (chosen.empty()) chosen = axis_dirs();
  std::vector<IVec> out;
  for (const auto& d : chosen) {
    out.push_back(d);
    bool minus_equivalent = false;
    for (const auto& R : site) { const IVec r = mul(R, d); if (r[0] == -d[0] && r[1] == -d[1] && r[2] == -d[2]) minus_equivalent = true; }
    if (!minus_equivalent) out.push_back({-d[0], -d[1], -d[2]});
  }
  return out;
}
int least_count(const std::vector<IMat>& site, bool diagonal) { return int(least_dirs(site, diagonal).size()); }

// the independent atoms of the supercell and each one's site symmetry (affine operations that fix it)
std::vector<int> independent(const Super& S) {
  std::vector<char> done(S.f.size(), 0);
  std::vector<int> out;
  for (size_t i = 0; i < S.f.size(); ++i) {
    if (done[i]) continue;
    out.push_back(int(i));
    for (const auto& op : S.ops) { const int k = S.find(S.apply(op, S.f[i]), S.z[i]); if (k >= 0) done[size_t(k)] = 1; }
  }
  return out;
}
std::vector<std::pair<IMat, Vec3>> site_ops(const Super& S, int i) {
  std::vector<std::pair<IMat, Vec3>> out;
  for (const auto& op : S.ops) if (S.find(S.apply(op, S.f[size_t(i)]), S.z[size_t(i)]) == i) out.push_back(op);
  return out;
}

}  // namespace

DisplacementCount count_displacements(const System& unit, const std::array<int, 3>& d3, const std::array<int, 3>& d2, double tol, double cutoff) {
  DisplacementCount c;
  const auto uops = crystal_operations(unit, tol);
  c.operations = int(uops.size());
  // fc2: least displacements (axes, +/- auto) of the independent atoms of the fc2 supercell
  {
    const Super S = make_super(unit, uops, d2, tol);
    for (int i : independent(S)) {
      std::vector<IMat> rs;
      for (const auto& op : site_ops(S, i)) rs.push_back(op.first);
      c.fc2 += least_count(rs, false);
    }
  }
  const Super S = make_super(unit, uops, d3, tol);
  const auto ind = independent(S);
  c.independent_atoms = int(ind.size());
  for (int i : ind) {
    const auto site = site_ops(S, i);
    std::vector<IMat> rs;
    for (const auto& op : site) rs.push_back(op.first);
    for (const auto& d1 : least_dirs(rs, false)) {
      ++c.first;
      // the site symmetry that keeps this displacement
      std::vector<std::pair<IMat, Vec3>> red;
      for (const auto& op : site) { const IVec r = mul(op.first, d1); if (r == d1) red.push_back(op); }
      // orbits of every supercell atom under it; per representative its bond symmetry and least displacements
      std::vector<char> done(S.f.size(), 0);
      for (size_t j = 0; j < S.f.size(); ++j) {
        if (done[j]) continue;
        std::vector<IMat> bond;
        for (const auto& op : red) {
          const int k = S.find(S.apply(op, S.f[j]), S.z[j]);
          if (k >= 0) done[size_t(k)] = 1;
          if (k == int(j)) bond.push_back(op.first);
        }
        if (cutoff > 0) {
          Vec3 df = S.f[j] - S.f[size_t(i)];
          for (auto& v : df) v -= std::round(v);
          if (norm(S.L[0] * df[0] + S.L[1] * df[1] + S.L[2] * df[2]) > cutoff) continue;
        }
        c.pairs += least_count(bond, true);
      }
    }
  }
  return c;
}

// ---------------------------------------------------------------- the plan
int valence_electrons(const System& s, const std::string& data_dir) {
  std::map<std::string, int> zv;
  try {
    const Json j = read_json_file(data_dir + "/dft/enmax.json");
    if (j.has("zval")) for (const auto& [k, v] : j["zval"].members()) zv[k] = int(v.number());
  } catch (...) {}
  int n = 0;
  for (const auto& a : s.atoms) {
    const auto it = zv.find(element(a.element).symbol);
    n += it != zv.end() ? it->second : a.element;   // unknown: all electrons (an upper bound)
  }
  return n;
}

double single_point_core_hours(int electrons) { const double x = electrons / 288.0; return 4.33 * x * x; }

Phono3pyPlan plan_phono3py(const System& unit, const Phono3pyOptions& o, const std::string& mode) {
  Phono3pyPlan p;
  p.label = o.label;
  p.atoms = int(unit.atoms.size());
  p.dimensionality = resolve_dimensionality(unit, o.dimensionality, &p.notes);
  const bool two_d = p.dimensionality == "2d";
  p.dim_fc3 = o.dim_fc3[0] > 0 ? o.dim_fc3 : supercell_for_length(unit, o.fc3_length);
  p.dim_fc2 = o.dim_fc2[0] > 0 ? o.dim_fc2 : supercell_for_length(unit, o.fc2_length);
  if (two_d) {
    if (o.dim_fc3[0] > 0 && o.dim_fc3[2] != 1) throw std::invalid_argument("2D: --dim must be 1 along c (the vacuum)");
    if (o.dim_fc2[0] > 0 && o.dim_fc2[2] != 1) throw std::invalid_argument("2D: --dim-fc2 must be 1 along c (the vacuum)");
    p.dim_fc3[2] = p.dim_fc2[2] = 1;
    const Vec3 n = cross(unit.cell.a, unit.cell.b);
    p.cell_height = std::fabs(dot(unit.cell.c, n)) / norm(n);
    p.thickness = o.thickness;
  }
  const int n3 = p.dim_fc3[0] * p.dim_fc3[1] * p.dim_fc3[2], n2 = p.dim_fc2[0] * p.dim_fc2[1] * p.dim_fc2[2];
  p.atoms_fc3 = p.atoms * n3, p.atoms_fc2 = p.atoms * n2;
  if (p.atoms_fc3 > 15000 || p.atoms_fc2 > 15000) throw std::invalid_argument("supercells above 15,000 atoms: choose smaller --dim / --dim-fc2");
  p.full = count_displacements(unit, p.dim_fc3, p.dim_fc2, 1e-3, mode == "random" ? 0.0 : o.cutoff_pair);
  if (mode != "random" && o.cutoff_pair > 0) p.notes.push_back("pairs farther apart than " + fmt("%g", o.cutoff_pair) + " Å left out (--cutoff-pair)");
  p.mode = mode;
  if (mode == "random") {
    p.rd = o.rd > 0 ? o.rd : std::max(50, (p.full.fc3() + 9) / 10);
    if (o.rd_auto) p.notes.push_back("phono3py --rd auto chooses the number of random supercells; the estimate assumes a tenth of the full set");
    p.fc3_supercells = p.rd;
  } else {
    p.fc3_supercells = p.full.fc3();
  }
  // supercell k-meshes (the guide's density), 1 along c for a sheet
  System s3 = unit, s2 = unit;
  s3.cell.a = unit.cell.a * p.dim_fc3[0], s3.cell.b = unit.cell.b * p.dim_fc3[1], s3.cell.c = unit.cell.c * p.dim_fc3[2];
  s2.cell.a = unit.cell.a * p.dim_fc2[0], s2.cell.b = unit.cell.b * p.dim_fc2[1], s2.cell.c = unit.cell.c * p.dim_fc2[2];
  p.kpoints_fc3 = kmesh_bulk(s3, o.kdens, two_d), p.kpoints_fc2 = kmesh_bulk(s2, o.kdens, two_d);
  // the chain axis and the q-point mesh: 74 Å across and 405 Å along the chain (the guide's PE mesh, 160 15 10, as lengths)
  if (o.chain_axis >= -1) p.chain_axis = o.chain_axis;
  else {
    const auto ax = chain_axes(unit);
    p.chain_axis = ax.size() == 1 ? ax[0] : -1;
  }
  const double len[3] = {norm(unit.cell.a), norm(unit.cell.b), norm(unit.cell.c)};
  for (int k = 0; k < 3; ++k)
    p.mesh[size_t(k)] = o.mesh[0] > 0 ? o.mesh[size_t(k)] : std::max(1, int(std::lround((k == p.chain_axis ? 405.0 : 74.0) / len[k])));
  if (two_d) p.mesh[2] = 1;
  p.electrons = valence_electrons(unit, o.data_dir);
  p.core_h_fc3 = single_point_core_hours(p.electrons * n3);
  p.core_h_fc2 = single_point_core_hours(p.electrons * n2);
  p.core_h_total = p.fc3_supercells * p.core_h_fc3 + p.full.fc2 * p.core_h_fc2;
  if (mode == "random")
    p.notes.push_back("random displacements: " + std::to_string(p.rd) + " fc3 supercells (phono3py --rd, fitted with symfc) instead of the full " +
                      std::to_string(p.full.fc3()));
  return p;
}

Json Phono3pyPlan::json() const {
  Json j = Json::object();
  j["label"] = label, j["atoms"] = atoms, j["electrons"] = electrons, j["dimensionality"] = dimensionality;
  j["dim_fc3"] = arr3(dim_fc3), j["dim_fc2"] = arr3(dim_fc2), j["atoms_fc3"] = atoms_fc3, j["atoms_fc2"] = atoms_fc2;
  Json f = Json::object();
  f["first"] = full.first, f["pairs"] = full.pairs, f["fc3"] = full.fc3(), f["fc2"] = full.fc2, f["independent_atoms"] = full.independent_atoms, f["operations"] = full.operations;
  j["full_set"] = f;
  j["mode"] = mode, j["fc3_supercells"] = fc3_supercells, j["fc2_supercells"] = full.fc2, j["rd"] = rd;
  j["kpoints_fc3"] = arr3(kpoints_fc3), j["kpoints_fc2"] = arr3(kpoints_fc2), j["mesh"] = arr3(mesh), j["chain_axis"] = chain_axis;
  j["encut"] = encut;
  j["core_hours_per_fc3_supercell"] = core_h_fc3, j["core_hours_per_fc2_supercell"] = core_h_fc2, j["core_hours_total"] = core_h_total;
  j["ranks_fc3"] = ranks_fc3, j["ranks_fc2"] = ranks_fc2;
  if (dimensionality == "2d") j["cell_height"] = cell_height, j["thickness"] = thickness;
  Json n = Json::array();
  for (const auto& x : notes) n.push_back(x);
  j["notes"] = n;
  return j;
}

std::string Phono3pyPlan::text() const {
  std::ostringstream o;
  const char* ax = chain_axis == 0 ? "a" : chain_axis == 1 ? "b" : chain_axis == 2 ? "c" : "none";
  o << label << " (" << dimensionality << "): " << atoms << " atoms per cell, chain axis " << ax << "\n";
  o << "  fc3 supercell " << dims(dim_fc3) << " (" << atoms_fc3 << " atoms), fc2 supercell " << dims(dim_fc2) << " (" << atoms_fc2 << " atoms)\n";
  o << "  full finite-difference set: " << full.fc3() << " fc3 supercells (" << full.first << " single + " << full.pairs << " pairs), " << full.fc2
    << " fc2 supercells; " << full.operations << " symmetry operations, " << full.independent_atoms << " independent atoms\n";
  o << "  this set (" << mode << "): " << fc3_supercells << " fc3 + " << full.fc2 << " fc2 single points\n";
  o << "  k-mesh fc3 " << dims(kpoints_fc3) << ", fc2 " << dims(kpoints_fc2) << "; q-mesh for kappa " << dims(mesh) << "\n";
  o << "  estimated cost: " << fmt("%.2f", core_h_fc3) << " core-h per fc3 supercell, " << fmt("%.1f", core_h_fc2) << " per fc2 supercell; total about "
    << fmt("%.0f", core_h_total) << " core-hours (+ the kappa step)\n";
  for (const auto& n : notes) o << "  note: " << n << "\n";
  return o.str();
}

// ---------------------------------------------------------------- the step scripts
namespace {

std::string head(const ClusterProfile& p, const std::string& name, int ntasks, int cpus, const std::string& time, const std::string& array = "") {
  std::string h = "#!/bin/bash\n#SBATCH --job-name=" + name.substr(0, 40) + "\n";
  if (!p.account.empty()) h += "#SBATCH --account=" + p.account + "\n";
  h += "#SBATCH --nodes=1\n#SBATCH --ntasks=" + std::to_string(ntasks) + "\n#SBATCH --cpus-per-task=" + std::to_string(cpus) + "\n";
  if (p.mem_per_cpu > 0) h += "#SBATCH --mem-per-cpu=" + std::to_string(p.mem_per_cpu) + "\n";
  h += "#SBATCH --time=" + time + "\n";
  if (!array.empty()) h += "#SBATCH --array=" + array + "\n";
  h += "#SBATCH --output=" + std::string(array.empty() ? "slurm-%j.out" : "logs/slurm-%A_%a.out") + "\n";
  return h;
}
std::string hms(double hours) {
  const long s = std::max(900L, long(hours * 3600));
  char b[32];
  std::snprintf(b, sizeof b, "%02ld:%02ld:00", s / 3600, (s % 3600) / 60);
  return b;
}

std::string env_script(const ClusterProfile& p) {
  std::string t =
      "#!/bin/bash\n"
      "## How the step scripts find phono3py and VASP: edit this file once for your cluster (sourced by every step).\n"
      "## phono3py >= 4 splits the setup (phono3py-init: -d, --cf3, --cf2) from the calculation (phono3py); 3.x does\n"
      "## both with phono3py and computes with phono3py-load (the guide's commands). Both are handled here.\n";
  if (!p.modules.empty()) t += "module purge > /dev/null 2>&1\nmodule load " + p.modules + "\n";
  t += "export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-1}\n"
       "if [ -n \"${PHONO3PY_CONDA:-}\" ] && [ -f \"$PHONO3PY_CONDA\" ]; then source \"$PHONO3PY_CONDA\" ${PHONO3PY_ENV:-}\n"
       "elif [ -f \"$HOME/miniconda3/bin/activate\" ]; then source \"$HOME/miniconda3/bin/activate\" ${PHONO3PY_ENV:-}; fi   # guide ch. 6\n"
       "if command -v phono3py-init > /dev/null 2>&1; then INIT=phono3py-init; RUN=phono3py     # phono3py >= 4\n"
       "else INIT=phono3py; RUN=phono3py-load; fi                                                # phono3py 3.x\n"
       "VASP=\"" + p.vasp + "\"\n";
  return t;
}

std::string task_script(const ClusterProfile& p) {
  std::string t =
      "#!/bin/bash\n"
      "## One displaced supercell (sourced by run_step2.sh, run_step3.sh and run_dryrun.sh): VASP in a scratch workspace,\n"
      "## then vasprun.xml checked (complete, with forces, at least MIN_XML bytes) and the run repeated up to two times if\n"
      "## not (the guide's restart fix); WAVECAR and the other large files are deleted, only vasprun.xml, OSZICAR and OUT\n"
      "## are kept (50 GB home quota). DONE / FAILED mark the folder.\n"
      "MIN_XML=${MIN_XML:-20000}\n"
      "vasprun_ok() {\n"
      "  [ -s \"$1\" ] || return 1\n"
      "  [ \"$(stat -c %s \"$1\" 2> /dev/null || stat -f %z \"$1\")\" -ge \"$MIN_XML\" ] || return 1\n"
      "  tail -c 400 \"$1\" | grep -q \"</modeling>\" && grep -q '<varray name=\"forces\"' \"$1\"\n"
      "}\n"
      "kpoints_of() { case $1 in disp-FC2-*) echo KPOINTS-fc2;; *) echo KPOINTS-fc3;; esac; }\n"
      "run_disp() {   # $1 = disp-00001 or disp-FC2-00001\n"
      "  local d=$1 set=$SET\n"
      "  if [ -f \"$set/$d/DONE\" ]; then echo \"$d already done\"; return 0; fi\n"
      "  rm -f \"$set/$d/FAILED\"\n"
      "  local ws=\"\" WS=ph3-${SLURM_JOB_ID:-0}-${d}\n";
  if (!p.ws_allocate.empty()) t += "  " + p.ws_allocate + " > /dev/null 2>&1 && ws=$(" + p.ws_find + " 2> /dev/null)\n";
  t += "  if [ -z \"$ws\" ] || [ ! -d \"$ws\" ]; then ws=$set/scratch/$d; mkdir -p \"$ws\"; fi\n"
       "  cp \"$set/$d/POSCAR\" \"$ws/POSCAR\"\n"
       "  cp \"$set/INCAR\" \"$set/POTCAR\" \"$ws/\"\n"
       "  cp \"$set/$(kpoints_of \"$d\")\" \"$ws/KPOINTS\"\n"
       "  cd \"$ws\" || return 1\n"
       "  local ok=0\n"
       "  for attempt in 1 2 3; do\n"
       "    ${CAPS_SRUN:-srun} $VASP > OUT\n"
       "    if vasprun_ok vasprun.xml && grep -qs \"EDIFF is reached\" OUTCAR; then ok=1; break; fi\n"
       "    echo \"$d: vasprun.xml incomplete or SCF not converged (attempt $attempt)\"\n"
       "  done\n"
       "  rm -f WAVECAR CHG CHGCAR XDATCAR PROCAR DOSCAR EIGENVAL vaspout.h5 PCDAT REPORT IBZKPT\n"
       "  cp vasprun.xml OSZICAR OUT \"$set/$d/\" 2> /dev/null\n"
       "  cd \"$set\"\n"
       "  if [ $ok = 1 ]; then touch \"$d/DONE\"; else touch \"$d/FAILED\"; fi\n";
  if (!p.ws_release.empty()) t += "  case $ws in \"$set\"/scratch/*) rm -rf \"$ws\";; *) rm -rf \"${ws:?}\"/*; " + p.ws_release + " > /dev/null 2>&1;; esac\n";
  else t += "  rm -rf \"${ws:?}\"\n";
  t += "  [ $ok = 1 ]\n"
       "}\n";
  return t;
}

const char* kAnalyzePy = R"PY(#!/usr/bin/env python3
"""Extracts a phono3py kappa-m*.hdf5 for CAPS (caps phono3py analyze) and plots it (the guide's ch. 8 plots).

usage: python3 analyze_kappa.py kappa-mXXX.hdf5 [--out kappa_extract] [--ts 300 ...] [--plot]

Writes into --out: kappa_summary.json (temperatures, the kappa tensor xx yy zz yz xz xy per temperature, the mesh and the
weight sum) and modes_T<T>.csv per temperature: frequency (THz), weight, mode kappa xx yy zz yz xz xy (W/m/K, summed over
the irreducible q-points and bands and divided by the weight sum they give kappa), gamma (THz), group velocity (THz*A).
Needs h5py and numpy (as phono3py itself); --plot needs matplotlib.
"""
import argparse, json, os, sys

import h5py
import numpy as np


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("hdf5")
    ap.add_argument("--out", default="kappa_extract")
    ap.add_argument("--ts", type=float, nargs="*", default=None)
    ap.add_argument("--plot", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    with h5py.File(a.hdf5, "r") as f:
        T = np.array(f["temperature"][:], dtype=float)
        kappa = np.array(f["kappa"][:], dtype=float)
        freq = np.array(f["frequency"][:], dtype=float)
        w = np.array(f["weight"][:], dtype=float)
        mk = np.array(f["mode_kappa"][:], dtype=float)
        g = np.array(f["gamma"][:], dtype=float)
        gv = np.array(f["group_velocity"][:], dtype=float)
        mesh = [int(x) for x in f["mesh"][:]] if "mesh" in f else []
    ws = float(w.sum())
    pick = range(len(T)) if not a.ts else [int(np.argmin(abs(T - t))) for t in a.ts]
    summary = {"file": os.path.abspath(a.hdf5), "mesh": mesh, "weight_sum": ws, "temperatures": [float(T[i]) for i in pick],
               "kappa": [[float(x) for x in kappa[i]] for i in pick],
               "check": "sum(mode_kappa) / sum(weight) = kappa", "modes": []}
    for i in pick:
        name = "modes_T%g.csv" % T[i]
        with open(os.path.join(a.out, name), "w") as out:
            out.write("frequency_THz,weight,k_xx,k_yy,k_zz,k_yz,k_xz,k_xy,gamma_THz,v_x,v_y,v_z\n")
            nq, nb = freq.shape
            for q in range(nq):
                for b in range(nb):
                    k = mk[i, q, b] / ws
                    v = gv[q, b]
                    out.write("%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g,%.8g\n" % (
                        freq[q, b], w[q], k[0], k[1], k[2], k[3], k[4], k[5], g[i, q, b], v[0], v[1], v[2]))
        summary["modes"].append(name)
    json.dump(summary, open(os.path.join(a.out, "kappa_summary.json"), "w"), indent=1)
    print("wrote", a.out, "(%d temperatures)" % len(summary["temperatures"]))
    if a.plot:
        plot(a.out, summary)


def plot(out, s):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    k = np.array(s["kappa"])
    fig, ax = plt.subplots()
    for c, lab in enumerate(["xx", "yy", "zz"]):
        ax.plot(s["temperatures"], k[:, c], "o-", label=lab)
    ax.set_xlabel("T (K)"); ax.set_ylabel("kappa (W/m/K)"); ax.legend(); fig.savefig(os.path.join(out, "kappa_T.png"), dpi=150)
    for name in s["modes"]:
        d = np.loadtxt(os.path.join(out, name), delimiter=",", skiprows=1)
        if d.ndim == 1:
            d = d[None, :]
        f, kk, g, v = d[:, 0], d[:, 2:5], d[:, 8], d[:, 9:12]
        tau = np.where(g > 0, 1.0 / (4 * np.pi * np.where(g > 0, g, 1)), 0.0)
        mfp = np.linalg.norm(v, axis=1) * tau / 10.0   # nm
        o = np.argsort(f)
        fig, ax = plt.subplots(1, 3, figsize=(14, 4))
        for c, lab in enumerate(["xx", "yy", "zz"]):
            ax[0].plot(f[o], np.cumsum(kk[o, c]), label=lab)
        ax[0].set_xlabel("frequency (THz)"); ax[0].set_ylabel("cumulative kappa (W/m/K)"); ax[0].legend()
        ax[1].scatter(f, tau, s=2); ax[1].set_yscale("log"); ax[1].set_xlabel("frequency (THz)"); ax[1].set_ylabel("lifetime (ps)")
        ax[2].scatter(f, mfp, s=2); ax[2].set_yscale("log"); ax[2].set_xlabel("frequency (THz)"); ax[2].set_ylabel("mean free path (nm)")
        fig.tight_layout(); fig.savefig(os.path.join(out, name.replace(".csv", ".png")), dpi=150)


if __name__ == "__main__":
    sys.exit(main())
)PY";

}  // namespace

Phono3pyResult write_phono3py_set(const System& unit0, const std::string& dir, const Phono3pyOptions& o) {
  Phono3pyResult res;
  const System unit = group_by_species(unit0, {"C", "N", "O", "F", "S", "Si", "Ti", "H"});
  // the mode: auto → full unless the full fc3 set is large; both → full/ and random/
  std::vector<std::string> modes;
  Phono3pyPlan probe = plan_phono3py(unit, o, "full");
  if (o.mode == "full" || o.mode == "random") modes = {o.mode};
  else if (o.mode == "both") modes = {"full", "random"};
  else if (o.mode == "auto") modes = {probe.full.fc3() > o.random_above ? "random" : "full"};
  else throw std::invalid_argument("--mode is full, random, both or auto");
  // ENCUT: relax4's (chain.json of the relaxation this cell came from), else about 2 × ENMAX
  const auto species = species_of(unit);
  double enmax = 0;
  std::string enc_note;
  if (!o.pp_dir.empty()) {
    const auto p = assemble_potcar(species, o.pp_dir, potcar_map(o.data_dir));
    for (double e : p.enmax) enmax = std::max(enmax, e);
  } else {
    const auto t = enmax_table(o.data_dir);
    for (const auto& el : species) { const auto it = t.find(el); if (it != t.end()) enmax = std::max(enmax, it->second); }
  }
  double encut = o.encut;
  if (encut <= 0) {
    if (enmax <= 0) throw std::invalid_argument("no ENMAX for the elements: give --pp-dir or --encut");
    encut = 10.0 * std::ceil(2.0 * enmax / 10.0 - 1e-9);
    enc_note = "ENCUT " + fmt("%g", encut) + " eV = about 2 × ENMAX (relax1/relax3/relax4's cutoff of caps vasp-relax); --encut for another";
  }
  ClusterProfile prof = cluster_profile(o.data_dir, o.profile);
  if (!o.account.empty()) prof.account = o.account;
  for (const auto& mode : modes) {
    const std::string d = modes.size() > 1 ? join(dir, mode) : dir;
    fs::create_directories(join(d, "logs"));
    Phono3pyPlan p = plan_phono3py(unit, o, mode);
    p.encut = encut;
    if (!enc_note.empty()) p.notes.push_back(enc_note);
    const bool two_d = p.dimensionality == "2d";
    const ClusterRule r3 = parallel_rule(prof, p.atoms_fc3, false), r2 = parallel_rule(prof, p.atoms_fc2, false);
    p.ranks_fc3 = r3.ranks, p.ranks_fc2 = r2.ranks;
    // dipole for a polar sheet
    double dipz = -1;
    if (two_d && (o.dipole == "on" || (o.dipole == "auto" && polar_along_c(unit)))) {
      double zmin = 1e300, zmax = -1e300;
      for (const auto& a : unit.atoms) { const double f = unit.cell.to_fractional(a.pos)[2]; zmin = std::min(zmin, f), zmax = std::max(zmax, f); }
      dipz = 0.5 * (zmin + zmax);
      p.notes.push_back("polar sheet: dipole correction (IDIPOL 3, LDIPOL) in every single point");
    }
    // files
    write_vasp_poscar(unit, join(d, "POSCAR-unitcell"), formula_ordered(unit, species) + " CAPS " + o.label + " unit cell");
    std::vector<IncarLine> L;
    auto add = [&](const std::string& k, const std::string& v, const std::string& why) { L.push_back({k, v, why}); };
    add("ISTART", "0", "new run");
    add("ICHARG", "2", "charge from atoms");
    add("ENCUT", fmt("%g", encut), "the relaxation's cutoff: forces of the displaced supercells consistent with the relaxed cell");
    add("PREC", "Accurate", "accurate forces (guide ch. 3.4)");
    add("ALGO", "Normal", "robust blocked Davidson");
    add("GGA", "PE", "PBE");
    if (o.ivdw) add("IVDW", std::to_string(o.ivdw), "D3-BJ dispersion, as in the relaxation");
    add("LASPH", ".TRUE.", "non-spherical PAW contributions");
    add("EDIFF", "1E-8", "tight SCF: finite-difference forces of 0.03 Å displacements");
    add("NELM", "200", "electronic steps");
    add("ISMEAR", "0", "Gaussian smearing (insulator)");
    add("SIGMA", "0.05", "smearing width");
    add("IBRION", "-1", "single point: no ionic update (guide ch. 3.4)");
    add("ISIF", "0", "no stress (forces only)");
    add("NSW", "0", "single point");
    add("LREAL", ".FALSE.", "exact projection: accurate small forces");
    add("NCORE", std::to_string(r3.ncore), "cores per orbital (cluster profile)");
    add("KPAR", std::to_string(r3.kpar), "k-point groups (cluster profile)");
    add("LWAVE", ".FALSE.", "no WAVECAR (quota)");
    add("LCHARG", ".FALSE.", "no CHGCAR (quota)");
    if (dipz >= 0) {
      add("IDIPOL", "3", "polar sheet: dipole correction along z");
      add("LDIPOL", ".TRUE.", "corrected potential: the forces carry it");
      add("DIPOL", "0.5 0.5 " + fmt("%.4f", dipz), "centre of the atoms (fractional)");
    }
    spit(join(d, "INCAR"), incar_text(o.label + " phono3py single point", L));
    spit(join(d, "KPOINTS-fc3"), "CAPS fc3 supercell " + dims(p.dim_fc3) + "\n0\nGamma\n  " + dims(p.kpoints_fc3) + "\n  0  0  0\n");
    spit(join(d, "KPOINTS-fc2"), "CAPS fc2 supercell " + dims(p.dim_fc2) + "\n0\nGamma\n  " + dims(p.kpoints_fc2) + "\n  0  0  0\n");
    const auto map = potcar_map(o.data_dir);
    if (!o.pp_dir.empty()) assemble_potcar(species, o.pp_dir, map, join(d, "POTCAR"));
    spit(join(d, "make_potcar.sh"), make_potcar_script(map));
    make_exec(join(d, "make_potcar.sh"));
    spit(join(d, "phono3py_env.sh"), env_script(prof));
    spit(join(d, "vasp_task.sh"), task_script(prof));
    spit(join(d, "analyze_kappa.py"), kAnalyzePy);
    const std::string lim = o.array_limit > 0 ? "%" + std::to_string(o.array_limit) : "";
    const std::string disp_args = "--dim " + dims(p.dim_fc3) + " --dim-fc2 " + dims(p.dim_fc2) + " --pa auto -c POSCAR-unitcell" +
                                  (mode == "random" ? std::string(" --rd ") + (o.rd_auto ? "auto" : std::to_string(p.rd)) + " --random-seed " + std::to_string(o.random_seed)
                                                     : (o.cutoff_pair > 0 ? " --cutoff-pair " + fmt("%g", o.cutoff_pair) : ""));
    const std::string common = "set -o pipefail\ncd \"${SLURM_SUBMIT_DIR:-$(dirname \"$0\")}\"\nSET=$(pwd)\nsource ./phono3py_env.sh\n";
    // step 1: displacements and their folders; the array lines of steps 2 and 3 set to the real counts
    std::string s1 = head(prof, "p3_step1", 1, 4, "01:00:00") +
                     "## Step 1 (guide 7.2): the displaced supercells (phono3py -d) and one folder per supercell.\n"
                     "## CAPS estimated " + std::to_string(p.fc3_supercells) + " fc3 and " + std::to_string(p.full.fc2) + " fc2 supercells.\n" + common +
                     "if [ -f fc3_dirs.txt ]; then echo \"step 1 already done ($(wc -l < fc3_dirs.txt) fc3, $(wc -l < fc2_dirs.txt) fc2)\"; exit 0; fi\n"
                     "$INIT -d " + disp_args + " > out_step1.out 2>&1 || { cat out_step1.out; exit 1; }\n"
                     ": > fc3_dirs.txt; : > fc2_dirs.txt\n"
                     "for f in $(ls POSCAR-[0-9]* 2> /dev/null | sort); do d=disp-${f#POSCAR-}; mkdir -p \"$d\"; mv \"$f\" \"$d/POSCAR\"; echo \"$d\" >> fc3_dirs.txt; done\n"
                     "for f in $(ls POSCAR_FC2-[0-9]* 2> /dev/null | sort); do d=disp-FC2-${f#POSCAR_FC2-}; mkdir -p \"$d\"; mv \"$f\" \"$d/POSCAR\"; echo \"$d\" >> fc2_dirs.txt; done\n"
                     "n3=$(wc -l < fc3_dirs.txt | tr -d ' '); n2=$(wc -l < fc2_dirs.txt | tr -d ' ')\n"
                     "[ \"$n3\" -gt 0 ] && [ \"$n2\" -gt 0 ] || { echo \"phono3py wrote no supercells (out_step1.out)\"; exit 1; }\n"
                     "sed -i.bak \"s/^#SBATCH --array=.*/#SBATCH --array=1-$n3" + lim + "/\" run_step2.sh && rm -f run_step2.sh.bak\n"
                     "sed -i.bak \"s/^#SBATCH --array=.*/#SBATCH --array=1-$n2" + lim + "/\" run_step3.sh && rm -f run_step3.sh.bak\n"
                     "echo \"fc3: $n3 supercells (CAPS estimated " + std::to_string(p.fc3_supercells) + "), fc2: $n2 (estimated " + std::to_string(p.full.fc2) + ")\"\n"
                     "echo \"next: sbatch run_dryrun.sh (a few supercells through the whole chain), then sbatch run_step2.sh and sbatch run_step3.sh\"\n";
    spit(join(d, "run_step1.sh"), s1);
    // steps 2 and 3: the arrays
    const double t3 = std::max(0.25, 2.0 * p.core_h_fc3 / r3.ranks), t2 = std::max(0.25, 2.0 * p.core_h_fc2 / r2.ranks);
    auto arr = [&](const std::string& name, const std::string& list, const ClusterRule& r, double t, int n, const std::string& what) {
      return head(prof, name, r.ranks, 1, hms(t), "1-" + std::to_string(n) + lim) + "## " + what +
             " Array task N runs line N of " + list + " (run_step1.sh sets the array to the real count;\n"
             "## run_dryrun.sh sets --time from the measured time). Rerun failed ones with sbatch --array=LIST " + name + ".sh.\n" + common +
             "source ./vasp_task.sh\nd=$(sed -n \"${SLURM_ARRAY_TASK_ID}p\" " + list + ")\n[ -n \"$d\" ] || { echo \"no line ${SLURM_ARRAY_TASK_ID} in " + list +
             "\"; exit 1; }\nrun_disp \"$d\"\n";
    };
    spit(join(d, "run_step2.sh"), arr("run_step2", "fc3_dirs.txt", r3, t3, p.fc3_supercells, "Step 2 (guide 7.3): the fc3 supercells."));
    spit(join(d, "run_step3.sh"), arr("run_step3", "fc2_dirs.txt", r2, t2, p.full.fc2, "Step 3 (guide 7.3): the fc2 supercells (submit after step 2; both may run together)."));
    // step 4: check + quota
    std::string s4 = head(prof, "p3_step4", 1, 1, "01:00:00") +
                     "## Step 4 (guide 7.4): every folder checked (DONE, vasprun.xml complete), large files deleted; prints the\n"
                     "## array lists to rerun. May run while steps 2 and 3 are still running (before the nightly quota check).\n" + common +
                     "source ./vasp_task.sh\n"
                     "check() {   # $1 list file, $2 step script\n"
                     "  local n=0 bad=\"\" miss=0 i=0\n"
                     "  while read -r d; do\n"
                     "    i=$((i + 1))\n"
                     "    rm -f \"$d\"/WAVECAR \"$d\"/CHG \"$d\"/CHGCAR \"$d\"/XDATCAR \"$d\"/vaspout.h5\n"
                     "    if [ -f \"$d/DONE\" ] && vasprun_ok \"$d/vasprun.xml\"; then n=$((n + 1))\n"
                     "    elif [ -f \"$d/FAILED\" ] || [ -f \"$d/vasprun.xml\" ]; then rm -f \"$d/DONE\"; bad=\"$bad${bad:+,}$i\"\n"
                     "    else miss=$((miss + 1)); fi\n"
                     "  done < \"$1\"\n"
                     "  echo \"$1: $n done, $miss not run yet, failed: ${bad:-none}\"\n"
                     "  [ -n \"$bad\" ] && echo \"  rerun: sbatch --array=$bad $2\"\n"
                     "  [ -z \"$bad\" ] && [ $miss = 0 ]\n"
                     "}\n"
                     "ok=0\ncheck fc3_dirs.txt run_step2.sh || ok=1\ncheck fc2_dirs.txt run_step3.sh || ok=1\n"
                     "du -sh . 2> /dev/null\n"
                     "[ $ok = 0 ] && echo \"all supercells done: sbatch run_step5.sh\"\nexit $ok\n";
    spit(join(d, "run_step4.sh"), s4);
    // step 5: forces
    std::string s5 = head(prof, "p3_step5", 1, 4, "02:00:00") +
                     "## Step 5 (guide 7.5): FORCES_FC3 and FORCES_FC2 from the vasprun.xml files, in the order of the displacements.\n" + common +
                     "bash ./run_step4.sh > out_step4_check.out || { cat out_step4_check.out; echo \"not every supercell is done: see run_step4.sh\"; exit 1; }\n"
                     "$INIT --cf3 $(sed 's|$|/vasprun.xml|' fc3_dirs.txt) > out_step5.out 2>&1 || { cat out_step5.out; exit 1; }\n"
                     "$INIT --cf2 $(sed 's|$|/vasprun.xml|' fc2_dirs.txt) >> out_step5.out 2>&1 || { cat out_step5.out; exit 1; }\n"
                     "[ -s FORCES_FC3 ] && [ -s FORCES_FC2 ] || { echo \"FORCES_FC3 / FORCES_FC2 missing (out_step5.out)\"; exit 1; }\n"
                     "echo \"FORCES_FC3 and FORCES_FC2 written: sbatch run_step6.sh\"\n";
    spit(join(d, "run_step5.sh"), s5);
    // step 6: force constants
    const int cpn = prof.cores_per_node > 0 ? prof.cores_per_node : 32;
    std::string s6 = head(prof, "p3_step6", 1, cpn, "12:00:00") +
                     "## Step 6 (guide 7.5): fc2.hdf5 and fc3.hdf5 " + std::string(mode == "random" ? "fitted by symfc (random displacements need it)" : "by finite differences") + ".\n" + common +
                     "$RUN" + (mode == "random" ? " --symfc" : "") + " -v > out_step6.out 2>&1 || { tail -30 out_step6.out; exit 1; }\n"
                     "[ -s fc3.hdf5 ] && [ -s fc2.hdf5 ] || { echo \"fc3.hdf5 / fc2.hdf5 missing (out_step6.out)\"; exit 1; }\n"
                     "echo \"force constants written: sbatch run_step7.sh\"\n";
    spit(join(d, "run_step6.sh"), s6);
    // step 7: kappa
    std::string ts;
    for (double t : o.temperatures) ts += (ts.empty() ? "" : " ") + fmt("%g", t);
    std::string s7 = head(prof, "p3_step7", 1, cpn, "48:00:00") +
                     "## Step 7 (guide 7.6): kappa on the " + dims(p.mesh) + " q-mesh, " + (o.method == "lbte" ? "direct solution (--lbte)" : "relaxation-time approximation (--br)") +
                     ", every temperature in one run.\n" + common +
                     "$RUN --fc3 --fc2 --mesh " + dims(p.mesh) + " --" + (o.method == "lbte" ? "lbte" : "br") + " --ts " + ts + " > out_run7.out 2>&1 || { tail -30 out_run7.out; exit 1; }\n"
                     "k=$(ls -t kappa-m*.hdf5 2> /dev/null | head -1)\n"
                     "[ -n \"$k\" ] || { echo \"no kappa-m*.hdf5 (out_run7.out)\"; exit 1; }\n"
                     "tail -20 out_run7.out\n"
                     "python3 analyze_kappa.py \"$k\" --out kappa_extract --plot || echo \"analyze_kappa.py needs h5py, numpy and matplotlib: run it where they are\"\n"
                     "echo \"then: caps phono3py analyze $SET\"\n";
    spit(join(d, "run_step7.sh"), s7);
    // dry run: step 1, then a few supercells through the same task code, timed; the arrays' time limits from it
    const int nd = std::max(1, o.dry_run);
    std::string dr = head(prof, "p3_dryrun", std::max(r3.ranks, r2.ranks), 1, hms(std::max(1.0, 1.5 * (nd * p.core_h_fc3 / r3.ranks + p.core_h_fc2 / r2.ranks)))) +
                     "## The dry run (before any array): step 1, then the first " + std::to_string(nd) + " fc3 supercells and the first fc2 supercell through\n"
                     "## exactly the code of steps 2 and 3, timed; vasprun.xml checked for forces. Projects the total core-hours from\n"
                     "## the measured times and sets --time of run_step2.sh / run_step3.sh to 1.5 x the measured time.\n" + common +
                     "bash ./run_step1.sh || exit 1\nsource ./vasp_task.sh\n"
                     "fmt_time() { local s=$1; [ $s -lt 900 ] && s=900; printf '%02d:%02d:00' $((s / 3600)) $(((s % 3600) / 60)); }\n"
                     "t3=0; k=0\n"
                     "for d in $(head -" + std::to_string(nd) + " fc3_dirs.txt); do\n"
                     "  t0=$(date +%s); run_disp \"$d\" || { echo \"dry run: $d failed (see $d/OUT)\"; exit 1; }; t=$(( $(date +%s) - t0 )); t3=$((t3 > t ? t3 : t)); k=$((k + 1))\n"
                     "  echo \"$d: $t s\"\n"
                     "done\n"
                     "d=$(head -1 fc2_dirs.txt); t0=$(date +%s); run_disp \"$d\" || { echo \"dry run: $d failed (see $d/OUT)\"; exit 1; }; t2=$(( $(date +%s) - t0 ))\n"
                     "echo \"$d: $t2 s\"\n"
                     "n3=$(wc -l < fc3_dirs.txt | tr -d ' '); n2=$(wc -l < fc2_dirs.txt | tr -d ' ')\n"
                     "R3=" + std::to_string(r3.ranks) + "; R2=" + std::to_string(r2.ranks) + "\n"
                     "awk -v n3=$n3 -v n2=$n2 -v t3=$t3 -v t2=$t2 -v r3=$R3 -v r2=$R2 -v e=" + fmt("%.0f", p.core_h_total) + " 'BEGIN {\n"
                     "  c = (n3 * t3 * r3 + n2 * t2 * r2) / 3600\n"
                     "  printf \"measured: %d s per fc3 supercell (%d ranks), %d s per fc2 supercell (%d ranks)\\n\", t3, r3, t2, r2\n"
                     "  printf \"projected: %d fc3 + %d fc2 supercells = %.0f core-hours (CAPS estimated %.0f)\\n\", n3, n2, c, e }' | tee dryrun_report.txt\n"
                     "sed -i.bak \"s/^#SBATCH --time=.*/#SBATCH --time=$(fmt_time $((t3 * 3 / 2)))/\" run_step2.sh && rm -f run_step2.sh.bak\n"
                     "sed -i.bak \"s/^#SBATCH --time=.*/#SBATCH --time=$(fmt_time $((t2 * 3 / 2)))/\" run_step3.sh && rm -f run_step3.sh.bak\n"
                     "echo \"dry run passed: the arrays' time limits are set; sbatch run_step2.sh, then sbatch run_step3.sh\" | tee -a dryrun_report.txt\n";
    spit(join(d, "run_dryrun.sh"), dr);
    // quota cleanup
    std::string cl = "#!/bin/bash\n"
                     "## Quota cleanup (guide 4.1 and 7.3): deletes WAVECAR, CHG, CHGCAR, XDATCAR, vaspout.h5 and the scratch folders of\n"
                     "## this set; with --outcar also every OUTCAR (vasprun.xml keeps the forces). Shows what it would delete unless --yes.\n"
                     "cd \"$(dirname \"$0\")\" || exit 1\nyes=0; outcar=0\nfor a in \"$@\"; do case $a in --yes) yes=1;; --outcar) outcar=1;; esac; done\n"
                     "big() { find . -type f \\( -name WAVECAR -o -name CHG -o -name CHGCAR -o -name XDATCAR -o -name vaspout.h5 $( [ $outcar = 1 ] && echo -o -name OUTCAR ) \\) \"$@\"; }\n"
                     "du -sh . 2> /dev/null\n"
                     "if [ $yes = 1 ]; then big -delete; rm -rf scratch; du -sh . 2> /dev/null\n"
                     "else echo \"$(big | wc -l | tr -d ' ') files would be deleted\"; big -exec du -ch {} + 2> /dev/null | tail -1; echo \"(dry run: ./cleanup.sh --yes to delete)\"; fi\n";
    spit(join(d, "cleanup.sh"), cl);
    for (const char* f : {"run_step1.sh", "run_step2.sh", "run_step3.sh", "run_step4.sh", "run_step5.sh", "run_step6.sh", "run_step7.sh", "run_dryrun.sh", "cleanup.sh", "analyze_kappa.py"})
      make_exec(join(d, f));
    spit(join(d, "plan.json"), p.json().dump(1));
    spit(join(d, "README.txt"),
         "CAPS phono3py set (" + mode + "): " + o.label + "\n\n" + p.text() +
         "\nOrder (TU Dresden group VASP guide, ch. 7):\n"
         "  ./make_potcar.sh (PP_DIR=...) unless POTCAR is here\n"
         "  sbatch run_dryrun.sh     step 1 + a few supercells timed; sets the arrays' time limits\n"
         "  sbatch run_step2.sh      fc3 array (" + std::to_string(p.fc3_supercells) + " estimated)\n"
         "  sbatch run_step3.sh      fc2 array (" + std::to_string(p.full.fc2) + ")\n"
         "  sbatch run_step4.sh      check + delete large files (also while 2 and 3 run); rerun the failed ones it lists\n"
         "  sbatch run_step5.sh      FORCES_FC3, FORCES_FC2\n"
         "  sbatch run_step6.sh      fc3.hdf5, fc2.hdf5\n"
         "  sbatch run_step7.sh      kappa-m*.hdf5 (+ analyze_kappa.py)\n"
         "  ./cleanup.sh [--yes]     quota cleanup at any time\n"
         "  caps phono3py check DIR / caps phono3py analyze DIR on your machine\n"
         "phono3py_env.sh says how phono3py is found (conda, version 3 or 4).\n");
    res.dirs.push_back(d);
    res.plans.push_back(p);
  }
  return res;
}

// ---------------------------------------------------------------- check
Json phono3py_check(const std::string& dir) {
  Json r = Json::object();
  auto one = [&](const std::string& list, const std::string& script) {
    Json j = Json::object();
    std::ifstream in(join(dir, list));
    if (!in) { j["error"] = list + " missing: run run_step1.sh first"; return j; }
    int done = 0, failed = 0, missing = 0, i = 0;
    std::string spec;
    for (std::string d; std::getline(in, d);) {
      if (d.empty()) continue;
      ++i;
      const std::string p = join(dir, d);
      const std::string xml = slurp(join(p, "vasprun.xml"));
      const bool complete = xml.size() > 200 && xml.find("</modeling>", xml.size() - 400) != std::string::npos && xml.find("name=\"forces\"") != std::string::npos;
      if (fs::exists(join(p, "DONE")) && complete) ++done;
      else if (fs::exists(join(p, "FAILED")) || !xml.empty()) ++failed, spec += (spec.empty() ? "" : ",") + std::to_string(i);
      else ++missing;
    }
    j["total"] = i, j["done"] = done, j["failed"] = failed, j["not_run"] = missing;
    if (!spec.empty()) j["rerun"] = "sbatch --array=" + spec + " " + script;
    return j;
  };
  r["fc3"] = one("fc3_dirs.txt", "run_step2.sh");
  r["fc2"] = one("fc2_dirs.txt", "run_step3.sh");
  const bool ready = !r["fc3"].has("error") && !r["fc2"].has("error") && r["fc3"].num("done", 0) == r["fc3"].num("total", -1) && r["fc2"].num("done", 0) == r["fc2"].num("total", -1);
  r["ready_for_step5"] = ready;
  for (const char* f : {"FORCES_FC3", "FORCES_FC2", "fc3.hdf5", "fc2.hdf5"}) r[f] = fs::exists(join(dir, f));
  return r;
}

// ---------------------------------------------------------------- analysis
Json analyze_kappa_extract(const std::string& ex, const std::string& out, int chain_axis, const std::string& dimensionality, double cell_height, double thickness) {
  const Json sum = read_json_file(join(ex, "kappa_summary.json"));
  fs::create_directories(out);
  Json res = Json::object();
  const bool two_d = dimensionality == "2d";
  res["dimensionality"] = dimensionality;
  res["chain_axis"] = chain_axis;
  Json rows = Json::array();
  const auto& temps = sum["temperatures"].items();
  const auto& kap = sum["kappa"].items();
  const auto& modes = sum["modes"].items();
  const double scale = two_d && thickness > 0 ? cell_height / thickness : 1.0;
  if (two_d) {
    res["cell_height"] = cell_height;
    if (thickness > 0) res["thickness"] = thickness, res["kappa_scale"] = scale;
    else res["note"] = "2D: give --thickness (Å; e.g. the interlayer distance of the bulk parent) for kappa of the sheet itself; the sheet conductance needs none";
  }
  std::ofstream tab(join(out, "kappa_tensor.csv"));
  tab << "T_K,k_xx,k_yy,k_zz,k_yz,k_xz,k_xy" << (two_d ? ",sheet_G_xx_nW_K,sheet_G_yy_nW_K" : "") << "\n";
  for (size_t t = 0; t < temps.size(); ++t) {
    Json row = Json::object();
    const double T = temps[t].number();
    row["T"] = T;
    std::array<double, 6> k{};
    for (size_t c = 0; c < 6; ++c) k[c] = kap[t].items()[c].number();
    Json kt = Json::array();
    for (double x : k) kt.push_back(x);
    row["kappa_cell"] = kt;   // phono3py's: per volume of the cell (with its vacuum for a sheet)
    tab << T;
    for (double x : k) tab << "," << x;
    if (two_d) {
      // sheet conductance G = kappa_cell × h (h: the cell height along c, vacuum included), W/K → nW/K
      const double gx = k[0] * cell_height * 1e-10 * 1e9, gy = k[1] * cell_height * 1e-10 * 1e9;
      row["sheet_conductance_nW_per_K"] = Json::array();
      row["sheet_conductance_nW_per_K"].push_back(gx), row["sheet_conductance_nW_per_K"].push_back(gy);
      tab << "," << gx << "," << gy;
      if (thickness > 0) {
        Json ks = Json::array();
        ks.push_back(k[0] * scale), ks.push_back(k[1] * scale), ks.push_back(k[5] * scale);
        row["kappa_sheet_xx_yy_xy"] = ks;
      }
    }
    tab << "\n";
    if (chain_axis >= 0 && chain_axis <= 2) {
      const double along = k[size_t(chain_axis)];
      double across = 0;
      for (int c = 0; c < 3; ++c) if (c != chain_axis) across += k[size_t(c)] / 2;
      row["along_chain"] = along, row["transverse_mean"] = across, row["anisotropy"] = across > 0 ? along / across : 0.0;
    }
    // the modes at this temperature
    if (t < modes.size()) {
      std::ifstream in(join(ex, modes[t].str()));
      std::string line;
      std::getline(in, line);
      struct M { double f, w, k[6], g, v[3], tau, mfp; };
      std::vector<M> ms;
      while (std::getline(in, line)) {
        M m{};
        std::stringstream ss(line);
        std::string x;
        std::vector<double> v;
        while (std::getline(ss, x, ',')) v.push_back(std::stod(x));
        if (v.size() < 12) continue;
        m.f = v[0], m.w = v[1];
        for (int c = 0; c < 6; ++c) m.k[c] = v[size_t(2 + c)];
        m.g = v[8];
        for (int c = 0; c < 3; ++c) m.v[c] = v[size_t(9 + c)];
        m.tau = m.g > 0 ? 1.0 / (4 * M_PI * m.g) : 0.0;                                     // ps
        m.mfp = std::sqrt(m.v[0] * m.v[0] + m.v[1] * m.v[1] + m.v[2] * m.v[2]) * m.tau * 0.1;   // THz·Å × ps = Å → nm
        ms.push_back(m);
      }
      double sum_k[3] = {0, 0, 0};
      for (const auto& m : ms) for (int c = 0; c < 3; ++c) sum_k[c] += m.k[c];
      Json chk = Json::array();
      for (int c = 0; c < 3; ++c) chk.push_back(sum_k[c]);
      row["sum_of_modes"] = chk;   // equals kappa_cell xx yy zz
      const std::string tag = fmt("%g", T);
      // cumulative against frequency and against mean free path
      auto cumulative = [&](const char* file, auto key, const char* col) {
        std::vector<size_t> o(ms.size());
        std::iota(o.begin(), o.end(), size_t(0));
        std::sort(o.begin(), o.end(), [&](size_t a, size_t b) { return key(ms[a]) < key(ms[b]); });
        std::ofstream c(join(out, std::string(file) + "_T" + tag + ".csv"));
        c << col << ",cum_xx,cum_yy,cum_zz\n";
        double acc[3] = {0, 0, 0};
        Json half = Json::array();
        bool hit[3] = {false, false, false};
        double at[3] = {0, 0, 0};
        for (size_t i : o) {
          for (int k2 = 0; k2 < 3; ++k2) {
            acc[k2] += ms[i].k[k2];
            if (!hit[k2] && sum_k[k2] > 0 && acc[k2] >= 0.5 * sum_k[k2]) hit[k2] = true, at[k2] = key(ms[i]);
          }
          c << key(ms[i]) << "," << acc[0] << "," << acc[1] << "," << acc[2] << "\n";
        }
        for (double x : at) half.push_back(x);
        return half;
      };
      row["half_kappa_below_frequency_THz"] = cumulative("cumulative_frequency", [](const M& m) { return m.f; }, "frequency_THz");
      row["half_kappa_below_mfp_nm"] = cumulative("cumulative_mfp", [](const M& m) { return m.mfp; }, "mfp_nm");
      // lifetimes and group velocities: kappa-weighted means, and the per-mode table
      std::ofstream lt(join(out, "modes_T" + tag + ".csv"));
      lt << "frequency_THz,lifetime_ps,group_velocity_km_s,mfp_nm,k_xx,k_yy,k_zz\n";
      double wk = 0, wt = 0, wv = 0, vmax = 0;
      for (const auto& m : ms) {
        const double v = std::sqrt(m.v[0] * m.v[0] + m.v[1] * m.v[1] + m.v[2] * m.v[2]) * 0.1;   // THz·Å = 100 m/s = 0.1 km/s
        const double kk = m.k[0] + m.k[1] + m.k[2];
        wk += kk, wt += kk * m.tau, wv += kk * v, vmax = std::max(vmax, v);
        lt << m.f << "," << m.tau << "," << v << "," << m.mfp << "," << m.k[0] << "," << m.k[1] << "," << m.k[2] << "\n";
      }
      row["modes"] = int(ms.size());
      row["kappa_weighted_lifetime_ps"] = wk > 0 ? wt / wk : 0.0;
      row["kappa_weighted_group_velocity_km_s"] = wk > 0 ? wv / wk : 0.0;
      row["max_group_velocity_km_s"] = vmax;
    }
    rows.push_back(row);
  }
  res["temperatures"] = rows;
  res["out"] = out;
  spit(join(out, "kappa_analysis.json"), res.dump(1));
  return res;
}

}  // namespace caps
