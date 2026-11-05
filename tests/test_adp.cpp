#include <gtest/gtest.h>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <sstream>

#include "caps/adp.hpp"
#include "caps/crystal.hpp"

using namespace caps;

namespace {
// R: a rotation about (1, 1, 1)/√3 by 40°
std::array<Vec3, 3> rotated_axes() {
  const double t = 40 * M_PI / 180, c = std::cos(t), s = std::sin(t);
  const Vec3 k{1 / std::sqrt(3.0), 1 / std::sqrt(3.0), 1 / std::sqrt(3.0)};
  auto rot = [&](const Vec3& v) { return v * c + cross(k, v) * s + k * (dot(k, v) * (1 - c)); };
  return {rot({1, 0, 0}), rot({0, 1, 0}), rot({0, 0, 1})};
}
// two atoms; atom 0 visits r0 ± s_k e_k (six frames), atom 1 stays; atom 0 sits at the cell's edge so its images wrap
Trajectory six_frames(const Vec3& shift_per_frame) {
  Trajectory t;
  Atom a;
  a.element = 18, a.pos = {19.98, 5, 5};
  Atom b;
  b.element = 18, b.pos = {10, 10, 10};
  t.topology.atoms = {a, b};
  t.topology.cell.a = {20, 0, 0}, t.topology.cell.b = {0, 20, 0}, t.topology.cell.c = {0, 0, 20};
  const auto e = rotated_axes();
  const double s[3] = {0.12, 0.2, 0.3};
  int f = 0;
  for (int k = 0; k < 3; ++k)
    for (int sg : {1, -1}) {
      Vec3 p0 = a.pos + e[size_t(k)] * (sg * s[k]) + shift_per_frame * double(f);
      Vec3 p1 = b.pos + shift_per_frame * double(f);
      for (int q = 0; q < 3; ++q) p0[size_t(q)] -= 20 * std::floor(p0[size_t(q)] / 20), p1[size_t(q)] -= 20 * std::floor(p1[size_t(q)] / 20);
      t.positions.push_back({p0, p1});
      t.cells.push_back(t.topology.cell);
      t.timesteps.push_back(f++);
    }
  return t;
}
}  // namespace

// U = (1/6) Σ (±s_k e_k)(±s_k e_k)ᵀ = R diag(s²/3) Rᵀ exactly, though the atom wraps across the cell edge; the ellipsoid
// axes are c s_k/√3 along e_k; B = 8π² U_iso
TEST(Adp, ExactTensorAcrossTheCellEdge) {
  const Trajectory t = six_frames({0, 0, 0});
  AdpOptions o;
  o.remove_drift = false;
  const auto r = displacement_tensors(t, {0, 1, 2, 3, 4, 5}, o);
  const auto e = rotated_axes();
  const double s[3] = {0.12, 0.2, 0.3};
  double U[3][3] = {};
  for (int k = 0; k < 3; ++k)
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j) U[i][j] += s[k] * s[k] / 3 * e[size_t(k)][size_t(i)] * e[size_t(k)][size_t(j)];
  const Sym6 want{U[0][0], U[1][1], U[2][2], U[0][1], U[0][2], U[1][2]};
  for (int q = 0; q < 6; ++q) EXPECT_NEAR(r.atoms[0].U[size_t(q)], want[size_t(q)], 1e-12) << q;
  EXPECT_NEAR(r.atoms[0].u_iso, (0.12 * 0.12 + 0.04 + 0.09) / 9, 1e-12);
  EXPECT_NEAR(r.atoms[0].b, 8 * M_PI * M_PI * r.atoms[0].u_iso, 1e-12);
  for (int k = 0; k < 3; ++k) {
    EXPECT_NEAR(r.atoms[0].axes[size_t(k)], r.scale * s[2 - k] / std::sqrt(3.0), 1e-9);
    EXPECT_NEAR(std::fabs(dot(r.atoms[0].dirs[size_t(k)], e[size_t(2 - k)])), 1.0, 1e-9);
  }
  EXPECT_NEAR(norm(t.topology.cell.minimum_image(r.atoms[0].mean - Vec3{19.98, 5, 5})), 0.0, 1e-12);
  for (double u : r.atoms[1].U) EXPECT_NEAR(u, 0.0, 1e-24);
}

// The whole cell drifting frame by frame leaves the tensors as they are without the drift
TEST(Adp, DriftRemoved) {
  const auto a = displacement_tensors(six_frames({0, 0, 0}), {0, 1, 2, 3, 4, 5});
  const auto b = displacement_tensors(six_frames({0.7, -0.3, 1.1}), {0, 1, 2, 3, 4, 5});
  for (size_t i = 0; i < 2; ++i)
    for (int q = 0; q < 6; ++q) EXPECT_NEAR(a.atoms[i].U[size_t(q)], b.atoms[i].U[size_t(q)], 1e-12);
}

// CIF basis: a round trip in a triclinic cell; an isotropic U gives U_eq = U (Fischer & Tillmanns's sum) and U^ij of
// an orthogonal cell are the Cartesian components; the 50 % and 99 % probability scales of the trivariate normal
TEST(Adp, CifBasisAndScales) {
  Cell c;
  c.a = {6.1, 0, 0}, c.b = {1.3, 7.2, 0}, c.c = {-0.9, 1.1, 8.4};
  const Sym6 U{0.031, 0.022, 0.045, 0.004, -0.006, 0.009};
  const Sym6 back = u_cartesian_from_cif(u_cif_from_cartesian(U, c), c);
  for (int q = 0; q < 6; ++q) EXPECT_NEAR(back[size_t(q)], U[size_t(q)], 1e-14);
  const Sym6 iso{0.02, 0.02, 0.02, 0, 0, 0};
  EXPECT_NEAR(u_eq_from_cif(u_cif_from_cartesian(iso, c), c), 0.02, 1e-14);
  EXPECT_NEAR(u_eq_from_cif(u_cif_from_cartesian(U, c), c), (U[0] + U[1] + U[2]) / 3, 1e-14);
  Cell o;
  o.a = {5, 0, 0}, o.b = {0, 6, 0}, o.c = {0, 0, 7};
  const Sym6 uo = u_cif_from_cartesian(U, o);
  for (int q = 0; q < 6; ++q) EXPECT_NEAR(uo[size_t(q)], U[size_t(q)], 1e-14);
  EXPECT_NEAR(ellipsoid_scale(0.5), 1.5382, 1e-4);
  EXPECT_NEAR(ellipsoid_scale(0.99), 3.3682, 1e-4);
}

// Sites of a 2 × 2 × 2 supercell of a one-atom cubic cell: one site whose U is the mean of the eight atoms'; a CIF
// written with the site's U_ij reads back the same tensor
TEST(Adp, SupercellSitesAndCifRoundTrip) {
  Trajectory t;
  for (int i = 0; i < 2; ++i)
    for (int j = 0; j < 2; ++j)
      for (int k = 0; k < 2; ++k) {
        Atom a;
        a.element = 29, a.pos = {3.0 * i + 0.4, 3.0 * j + 0.4, 3.0 * k + 0.4};
        t.topology.atoms.push_back(a);
      }
  t.topology.cell.a = {6, 0, 0}, t.topology.cell.b = {0, 6, 0}, t.topology.cell.c = {0, 0, 6};
  // atom n moves ±(0.05 (n + 1)) along x over two frames, and the other way along y
  for (int f = 0; f < 2; ++f) {
    std::vector<Vec3> p;
    for (size_t n = 0; n < 8; ++n) {
      const double d = (f ? -1 : 1) * 0.05 * double(n + 1);
      p.push_back(t.topology.atoms[n].pos + Vec3{d, -0.5 * d, 0});
    }
    t.positions.push_back(p);
    t.cells.push_back(t.topology.cell);
    t.timesteps.push_back(f);
  }
  AdpOptions o;
  o.supercell = {2, 2, 2};
  o.remove_drift = false;
  const auto r = displacement_tensors(t, {0, 1}, o);
  ASSERT_EQ(r.sites.size(), 1u);
  double uxx = 0;
  for (size_t n = 0; n < 8; ++n) uxx += std::pow(0.05 * double(n + 1), 2) / 8;
  EXPECT_NEAR(r.sites[0].U[0], uxx, 1e-12);
  EXPECT_NEAR(r.sites[0].U[1], uxx / 4, 1e-12);
  EXPECT_NEAR(r.sites[0].U[3], -uxx / 2, 1e-12);
  const auto path = (std::filesystem::temp_directory_path() / "caps_adp.cif").string();
  write_cif_adp({29}, r.sites, r.unit_cell, path, "test");
  std::ifstream in(path);
  std::stringstream ss;
  ss << in.rdbuf();
  const auto d = cif_displacements(ss.str());
  ASSERT_EQ(d.u_aniso.count("Cu1"), 1u);
  const Sym6 back = u_cartesian_from_cif(d.u_aniso.at("Cu1"), r.unit_cell);
  for (int q = 0; q < 6; ++q) EXPECT_NEAR(back[size_t(q)], r.sites[0].U[size_t(q)], 1e-6);
  EXPECT_NEAR(d.u_iso.at("Cu1"), r.sites[0].u_iso, 1e-6);
}
