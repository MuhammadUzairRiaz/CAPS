// CAPS atomic displacement parameters from a trajectory (see adp.hpp).
#include "caps/adp.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <stdexcept>

#include "caps/analysis.hpp"
#include "caps/elements.hpp"
#include "caps/io.hpp"

namespace caps {

namespace {
constexpr double kPi = 3.14159265358979323846;

// the 3 × 3 of the cell vectors as columns, and its inverse
void cell_matrix(const Cell& c, double A[3][3]) {
  for (int i = 0; i < 3; ++i) A[i][0] = c.a[size_t(i)], A[i][1] = c.b[size_t(i)], A[i][2] = c.c[size_t(i)];
}
void invert3(const double M[3][3], double R[3][3]) {
  const double det = M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0]) + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]);
  if (std::fabs(det) < 1e-300) throw std::invalid_argument("a singular cell");
  R[0][0] = (M[1][1] * M[2][2] - M[1][2] * M[2][1]) / det, R[0][1] = (M[0][2] * M[2][1] - M[0][1] * M[2][2]) / det, R[0][2] = (M[0][1] * M[1][2] - M[0][2] * M[1][1]) / det;
  R[1][0] = (M[1][2] * M[2][0] - M[1][0] * M[2][2]) / det, R[1][1] = (M[0][0] * M[2][2] - M[0][2] * M[2][0]) / det, R[1][2] = (M[0][2] * M[1][0] - M[0][0] * M[1][2]) / det;
  R[2][0] = (M[1][0] * M[2][1] - M[1][1] * M[2][0]) / det, R[2][1] = (M[0][1] * M[2][0] - M[0][0] * M[2][1]) / det, R[2][2] = (M[0][0] * M[1][1] - M[0][1] * M[1][0]) / det;
}
void full(const Sym6& u, double U[3][3]) {
  U[0][0] = u[0], U[1][1] = u[1], U[2][2] = u[2], U[0][1] = U[1][0] = u[3], U[0][2] = U[2][0] = u[4], U[1][2] = U[2][1] = u[5];
}
Sym6 packed(const double U[3][3]) { return {U[0][0], U[1][1], U[2][2], 0.5 * (U[0][1] + U[1][0]), 0.5 * (U[0][2] + U[2][0]), 0.5 * (U[1][2] + U[2][1])}; }
// reciprocal lengths a*, b*, c*
std::array<double, 3> recip_lengths(const Cell& c) {
  const double V = c.volume();
  return {norm(cross(c.b, c.c)) / V, norm(cross(c.c, c.a)) / V, norm(cross(c.a, c.b)) / V};
}
// χ²₃ CDF
double chi2_3(double x) { return x <= 0 ? 0 : std::erf(std::sqrt(x / 2)) - std::sqrt(2 * x / kPi) * std::exp(-x / 2); }
}  // namespace

double ellipsoid_scale(double p) {
  if (!(p > 0 && p < 1)) throw std::invalid_argument("ellipsoid probability between 0 and 1");
  double lo = 0, hi = 100;
  for (int k = 0; k < 200; ++k) {
    const double mid = 0.5 * (lo + hi);
    (chi2_3(mid) < p ? lo : hi) = mid;
  }
  return std::sqrt(0.5 * (lo + hi));
}

void finish_adp(AtomAdp& a, double c) {
  double U[3][3], w[3], V[3][3];
  full(a.U, U);
  symmetric_eigen3(U, w, V);
  a.u_iso = (a.U[0] + a.U[1] + a.U[2]) / 3;
  a.b = 8 * kPi * kPi * a.u_iso;
  for (int k = 0; k < 3; ++k) {   // largest first
    const int src = 2 - k;
    a.axes[size_t(k)] = c * std::sqrt(std::max(0.0, w[src]));
    a.dirs[size_t(k)] = {V[0][src], V[1][src], V[2][src]};
  }
  // λmin relative to λmax: a smallest eigenvalue at round-off level (a rank-deficient tensor from few frames) is zero
  a.anisotropy = w[2] > 0 && w[0] > 1e-8 * w[2] ? w[2] / w[0] : std::numeric_limits<double>::infinity();
}

Sym6 u_cif_from_cartesian(const Sym6& u, const Cell& cell) {
  double A[3][3], Ai[3][3], U[3][3], R[3][3];
  cell_matrix(cell, A);
  invert3(A, Ai);
  full(u, U);
  const auto n = recip_lengths(cell);
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) {
      double s = 0;
      for (int k = 0; k < 3; ++k)
        for (int l = 0; l < 3; ++l) s += Ai[i][k] * U[k][l] * Ai[j][l];
      R[i][j] = s / (n[size_t(i)] * n[size_t(j)]);
    }
  return packed(R);
}

Sym6 u_cartesian_from_cif(const Sym6& u, const Cell& cell) {
  double A[3][3], U[3][3], R[3][3];
  cell_matrix(cell, A);
  full(u, U);
  const auto n = recip_lengths(cell);
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) {
      double s = 0;
      for (int k = 0; k < 3; ++k)
        for (int l = 0; l < 3; ++l) s += A[i][k] * n[size_t(k)] * U[k][l] * n[size_t(l)] * A[j][l];
      R[i][j] = s;
    }
  return packed(R);
}

double u_eq_from_cif(const Sym6& u, const Cell& cell) {
  double U[3][3];
  full(u, U);
  const auto n = recip_lengths(cell);
  const Vec3 a[3] = {cell.a, cell.b, cell.c};
  double s = 0;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j) s += U[i][j] * n[size_t(i)] * n[size_t(j)] * dot(a[i], a[j]);
  return s / 3;
}

AdpResult displacement_tensors(const Trajectory& t, const std::vector<size_t>& frames, const AdpOptions& o) {
  AdpResult R;
  const size_t n = t.topology.atoms.size();
  if (frames.size() < 2) throw std::invalid_argument("displacement tensors: two frames or more");
  if (n == 0) throw std::invalid_argument("displacement tensors: no atoms");
  for (int k = 0; k < 3; ++k)
    if (o.supercell[size_t(k)] < 1) throw std::invalid_argument("displacement tensors: each supercell repeat must be 1 or more");
  R.frames = int(frames.size());
  R.scale = ellipsoid_scale(o.probability);
  std::vector<double> m(n);
  double mt = 0;
  for (size_t i = 0; i < n; ++i) m[i] = t.topology.mass_of(t.topology.atoms[i]), mt += m[i];
  // positions followed continuously, the centre's drift removed
  std::vector<Vec3> cur(t.positions[frames[0]].begin(), t.positions[frames[0]].end());
  std::vector<std::vector<Vec3>> X;
  X.reserve(frames.size());
  for (size_t f = 0; f < frames.size(); ++f) {
    const auto& P = t.positions[frames[f]];
    if (P.size() != n) throw std::invalid_argument("displacement tensors: a frame with another number of atoms");
    const Cell& c = frames[f] < t.cells.size() ? t.cells[frames[f]] : t.topology.cell;
    if (f > 0)
      for (size_t i = 0; i < n; ++i) {
        Vec3 d = P[i] - cur[i];
        if (c.valid()) d = c.minimum_image(d);
        cur[i] = cur[i] + d;
      }
    std::vector<Vec3> x = cur;
    if (o.remove_drift) {
      Vec3 com{0, 0, 0};
      for (size_t i = 0; i < n; ++i) com = com + x[i] * (m[i] / mt);
      for (auto& v : x) v = v - com;
    }
    X.push_back(std::move(x));
  }
  // the centre of the first frame back in, so the mean positions sit where the atoms are
  Vec3 com0{0, 0, 0};
  for (size_t i = 0; i < n; ++i) com0 = com0 + t.positions[frames[0]][i] * (m[i] / mt);
  if (!o.remove_drift) com0 = {0, 0, 0};
  const double nf = double(frames.size());
  R.atoms.resize(n);
  for (size_t i = 0; i < n; ++i) {
    Vec3 mean{0, 0, 0};
    for (const auto& x : X) mean = mean + x[i];
    mean = mean * (1.0 / nf);
    Sym6 U{0, 0, 0, 0, 0, 0};
    for (const auto& x : X) {
      const Vec3 d = x[i] - mean;
      U[0] += d[0] * d[0], U[1] += d[1] * d[1], U[2] += d[2] * d[2], U[3] += d[0] * d[1], U[4] += d[0] * d[2], U[5] += d[1] * d[2];
    }
    for (auto& u : U) u /= nf;
    AtomAdp& a = R.atoms[i];
    a.mean = mean + com0;
    a.U = U;
    finish_adp(a, R.scale);
  }
  // sites of the unit cell
  const Cell& cell = t.topology.cell;
  if (o.supercell[0] < 1 || o.supercell[1] < 1 || o.supercell[2] < 1)
    throw std::invalid_argument("displacement tensors: the supercell repeats must be whole numbers of at least 1");
  const bool sc = o.supercell[0] > 1 || o.supercell[1] > 1 || o.supercell[2] > 1;
  if (sc) {
    if (!cell.valid()) throw std::invalid_argument("displacement tensors: sites need the supercell's cell");
    Cell u = cell;
    u.a = cell.a * (1.0 / o.supercell[0]), u.b = cell.b * (1.0 / o.supercell[1]), u.c = cell.c * (1.0 / o.supercell[2]);
    R.unit_cell = u;
    // each site of the unit cell has one image per repeat, so the repeats must divide the atoms
    const size_t reps = size_t(o.supercell[0]) * size_t(o.supercell[1]) * size_t(o.supercell[2]);
    const std::string sc_name = std::to_string(o.supercell[0]) + " x " + std::to_string(o.supercell[1]) + " x " + std::to_string(o.supercell[2]);
    if (n % reps != 0)
      throw std::invalid_argument("displacement tensors: " + std::to_string(n) + " atoms cannot split into " + sc_name + " = " + std::to_string(reps) +
                                  " images of a unit cell; the repeats must divide the trajectory's cell");
    std::vector<Vec3> fold;   // fractional in the unit cell
    for (size_t i = 0; i < n; ++i) {
      Vec3 f = u.to_fractional(R.atoms[i].mean);
      for (auto& x : f) x -= std::floor(x);
      int s = -1;
      for (size_t k = 0; k < fold.size() && s < 0; ++k) {
        Vec3 df = f - fold[k];
        for (auto& x : df) x -= std::round(x);
        if (norm(u.to_cartesian(df) - u.to_cartesian({0, 0, 0})) < 0.5 && t.topology.atoms[R.site_atoms[k][0]].element == t.topology.atoms[i].element) s = int(k);
      }
      if (s < 0) s = int(fold.size()), fold.push_back(f), R.site_atoms.push_back({});
      R.site_atoms[size_t(s)].push_back(uint32_t(i));
    }
    // every site must collect exactly one atom per image; anything else means the repeats are not the trajectory's period
    bool even = fold.size() * reps == n;
    for (const auto& sa : R.site_atoms) even = even && sa.size() == reps;
    if (!even) {
      std::map<size_t, size_t> counts;   // atoms per site -> number of sites
      for (const auto& sa : R.site_atoms) ++counts[sa.size()];
      std::string found;
      for (const auto& [per, ns] : counts) found += (found.empty() ? "" : ", ") + std::to_string(ns) + " of " + std::to_string(per) + (per == 1 ? " atom" : " atoms");
      throw std::invalid_argument("displacement tensors: the supercell repeats " + sc_name + " do not match the trajectory's periodicity: each site should hold " +
                                  std::to_string(reps) + " atoms (one per image), but the atoms fold onto " + std::to_string(fold.size()) +
                                  " sites (" + found + ")");
    }
    for (size_t k = 0; k < fold.size(); ++k) {
      AtomAdp a;
      const double ns = double(R.site_atoms[k].size());
      // the mean over all images of the site, each taken by its minimum image about the site's first atom so that
      // images near a cell face do not split
      Vec3 acc{0, 0, 0};
      for (uint32_t i : R.site_atoms[k]) {
        Vec3 df = u.to_fractional(R.atoms[i].mean);
        for (int q = 0; q < 3; ++q) df[size_t(q)] -= fold[k][size_t(q)], df[size_t(q)] -= std::round(df[size_t(q)]);
        acc = acc + df * (1.0 / ns);
        for (int q = 0; q < 6; ++q) a.U[size_t(q)] += R.atoms[i].U[size_t(q)] / ns;
      }
      Vec3 fm = fold[k] + acc;
      for (auto& x : fm) x -= std::floor(x);
      a.mean = u.to_cartesian(fm);
      finish_adp(a, R.scale);
      R.sites.push_back(a);
    }
    R.notes.push_back(std::to_string(fold.size()) + " sites in the unit cell, " + std::to_string(reps) + " atoms each");
  } else {
    R.unit_cell = cell;
  }
  if (frames.size() < 50) R.notes.push_back("only " + std::to_string(frames.size()) + " frames: the tensors are rough");
  return R;
}

void write_cif_adp(const std::vector<int>& elements, const std::vector<AtomAdp>& rows, const Cell& cell, const std::string& path, const std::string& title) {
  if (!cell.valid()) throw std::runtime_error("a CIF needs a periodic cell");
  if (elements.size() != rows.size()) throw std::invalid_argument("write_cif_adp: an element per row");
  std::ofstream out(path);
  if (!out) throw std::runtime_error("cannot write " + path);
  auto ang = [](const Vec3& u, const Vec3& v) { return std::acos(std::clamp(dot(u, v) / (norm(u) * norm(v)), -1.0, 1.0)) * 180 / kPi; };
  const std::string name = cif_block_name(title);   // printable ASCII, at most 75 characters, never empty
  char b[256];
  out << "# written by CAPS: mean positions and displacement parameters from a trajectory\ndata_" << name << "\n"
      << "_space_group_name_H-M_alt 'P 1'\n_space_group_IT_number 1\n_symmetry_space_group_name_H-M 'P 1'\n_symmetry_Int_Tables_number 1\n";
  std::snprintf(b, sizeof b, "_cell_length_a %.6f\n_cell_length_b %.6f\n_cell_length_c %.6f\n_cell_angle_alpha %.6f\n_cell_angle_beta %.6f\n_cell_angle_gamma %.6f\n",
                norm(cell.a), norm(cell.b), norm(cell.c), ang(cell.b, cell.c), ang(cell.a, cell.c), ang(cell.a, cell.b));
  out << b << "loop_\n_space_group_symop_operation_xyz\n'x, y, z'\n";
  out << "loop_\n_atom_site_label\n_atom_site_type_symbol\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n_atom_site_U_iso_or_equiv\n_atom_site_adp_type\n";
  std::map<int, int> count;
  std::vector<std::string> labels;
  for (size_t i = 0; i < rows.size(); ++i) {
    Vec3 f = cell.to_fractional(rows[i].mean);
    for (auto& x : f) x -= std::floor(x);
    const char* sym = element(elements[i]).symbol;
    labels.push_back(std::string(sym) + std::to_string(++count[elements[i]]));
    std::snprintf(b, sizeof b, "%s %s %.6f %.6f %.6f %.6f Uani\n", labels.back().c_str(), sym, f[0], f[1], f[2], rows[i].u_iso);
    out << b;
  }
  out << "loop_\n_atom_site_aniso_label\n_atom_site_aniso_U_11\n_atom_site_aniso_U_22\n_atom_site_aniso_U_33\n_atom_site_aniso_U_12\n_atom_site_aniso_U_13\n_atom_site_aniso_U_23\n";
  for (size_t i = 0; i < rows.size(); ++i) {
    const Sym6 u = u_cif_from_cartesian(rows[i].U, cell);
    std::snprintf(b, sizeof b, "%s %.6f %.6f %.6f %.6f %.6f %.6f\n", labels[i].c_str(), u[0], u[1], u[2], u[3], u[4], u[5]);
    out << b;
  }
}

}  // namespace caps
