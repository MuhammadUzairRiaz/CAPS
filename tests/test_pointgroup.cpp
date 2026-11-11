#include <gtest/gtest.h>

#include <cmath>

#include "caps/pointgroup.hpp"

using namespace caps;

namespace {
System mol(const std::vector<std::pair<int, Vec3>>& atoms) {
  System s;
  for (const auto& [z, p] : atoms) {
    Atom a;
    a.element = z, a.pos = p;
    s.atoms.push_back(a);
  }
  return s;
}
const double t = 1.09 / std::sqrt(3.0);
System methane() { return mol({{6, {0, 0, 0}}, {1, {t, t, t}}, {1, {t, -t, -t}}, {1, {-t, t, -t}}, {1, {-t, -t, t}}}); }
System ethane(bool staggered) {
  std::vector<std::pair<int, Vec3>> a{{6, {0, 0, -0.765}}, {6, {0, 0, 0.765}}};
  for (int k = 0; k < 3; ++k) {
    const double p = 2 * M_PI * k / 3, q = p + (staggered ? M_PI / 3 : 0);
    a.push_back({1, {1.02 * std::cos(p), 1.02 * std::sin(p), -1.16}});
    a.push_back({1, {1.02 * std::cos(q), 1.02 * std::sin(q), 1.16}});
  }
  return mol(a);
}
}  // namespace

// The groups of molecules with exact geometries, and their orders
TEST(PointGroup, KnownMolecules) {
  struct Case { const char* name; System s; std::string group; int order; };
  std::vector<std::pair<int, Vec3>> benz, sf6{{16, {0, 0, 0}}}, icos;
  for (int k = 0; k < 6; ++k) {
    const double p = M_PI * k / 3;
    benz.push_back({6, {1.39 * std::cos(p), 1.39 * std::sin(p), 0}});
    benz.push_back({1, {2.47 * std::cos(p), 2.47 * std::sin(p), 0}});
  }
  for (int k = 0; k < 3; ++k)
    for (double sg : {-1.56, 1.56}) { Vec3 v{0, 0, 0}; v[size_t(k)] = sg; sf6.push_back({9, v}); }
  const double phi = (1 + std::sqrt(5.0)) / 2;
  for (double a : {-1.0, 1.0})
    for (double b : {-phi, phi}) {
      icos.push_back({5, {0, a, b}});
      icos.push_back({5, {a, b, 0}});
      icos.push_back({5, {b, 0, a}});
    }
  const std::vector<Case> cases = {
      {"water", mol({{8, {0, 0, 0.1173}}, {1, {0, 0.7572, -0.4692}}, {1, {0, -0.7572, -0.4692}}}), "C2v", 4},
      {"ammonia", mol({{7, {0, 0, 0.1}}, {1, {0.94, 0, -0.27}}, {1, {0.94 * std::cos(2 * M_PI / 3), 0.94 * std::sin(2 * M_PI / 3), -0.27}},
                       {1, {0.94 * std::cos(4 * M_PI / 3), 0.94 * std::sin(4 * M_PI / 3), -0.27}}}), "C3v", 6},
      {"methane", methane(), "Td", 24},
      {"benzene", mol(benz), "D6h", 24},
      {"SF6", mol(sf6), "Oh", 48},
      {"B12 icosahedron", mol(icos), "Ih", 120},
      {"staggered ethane", ethane(true), "D3d", 12},
      {"eclipsed ethane", ethane(false), "D3h", 12},
      {"H2O2", mol({{8, {0.7, 0, 0}}, {8, {-0.7, 0, 0}}, {1, {0.9, 0.9, 0.2}}, {1, {-0.9, -0.9, 0.2}}}), "C2", 2},
      {"CHFClBr", mol({{6, {0, 0, 0}}, {1, {0.6, 0.6, 0.6}}, {9, {0.8, -0.8, -0.8}}, {17, {-1, 1, -1}}, {35, {-1.1, -1.1, 1.1}}}), "C1", 1},
  };
  for (const auto& c : cases) {
    const auto r = point_group(c.s, 0.05);
    EXPECT_EQ(r.group, c.group) << c.name;
    EXPECT_EQ(r.order, c.order) << c.name;
    EXPECT_LT(r.max_deviation, 1e-9) << c.name;
  }
  EXPECT_EQ(point_group(mol({{8, {0, 0, -1.16}}, {6, {0, 0, 0}}, {8, {0, 0, 1.16}}})).group, "D*h");
  EXPECT_EQ(point_group(mol({{1, {0, 0, -1.06}}, {6, {0, 0, 0}}, {7, {0, 0, 1.15}}})).group, "C*v");
}

// Methane with every atom moved by up to 0.03 Å is found Td at 0.1 Å, and symmetrised to an exact tetrahedron: the
// group's images meet to 10⁻⁹ Å, the four C–H bonds equal, every H–C–H angle the tetrahedral one
TEST(PointGroup, SymmetrizeNoisyMethane) {
  System s = methane();
  const Vec3 noise[5] = {{0.01, -0.02, 0.015}, {0.03, 0, -0.01}, {-0.02, 0.025, 0.01}, {0.015, -0.01, -0.03}, {-0.01, -0.02, 0.02}};
  for (int k = 0; k < 5; ++k) s.atoms[size_t(k)].pos = s.atoms[size_t(k)].pos + noise[k];
  EXPECT_EQ(point_group(s, 0.01).group, "C1");   // too noisy for a tight tolerance
  const auto r = symmetrize_molecule(s, 0.1);
  EXPECT_EQ(r.group, "Td");
  EXPECT_LT(r.max_deviation, 1e-9);
  const Vec3 c = s.atoms[0].pos;
  const double d1 = norm(s.atoms[1].pos - c);
  for (int k = 2; k < 5; ++k) EXPECT_NEAR(norm(s.atoms[size_t(k)].pos - c), d1, 1e-9);
  for (int i = 1; i < 5; ++i)
    for (int j = i + 1; j < 5; ++j) {
      const Vec3 a = s.atoms[size_t(i)].pos - c, b = s.atoms[size_t(j)].pos - c;
      EXPECT_NEAR(std::acos(dot(a, b) / (norm(a) * norm(b))) * 180 / M_PI, std::acos(-1.0 / 3.0) * 180 / M_PI, 1e-7);
    }
}

// A planar molecule: benzene with every atom moved by up to 0.04 Å (out of the plane too) symmetrised to D6h: flat,
// a regular hexagon, the group's images meeting to 10⁻⁹ Å
TEST(PointGroup, SymmetrizeNoisyBenzene) {
  std::vector<std::pair<int, Vec3>> a;
  for (int k = 0; k < 6; ++k) {
    const double p = M_PI * k / 3;
    a.push_back({6, {1.39 * std::cos(p) + 0.03 * std::sin(3.1 * k), 1.39 * std::sin(p) - 0.02 * std::cos(1.7 * k), 0.04 * std::sin(2.3 * k)}});
    a.push_back({1, {2.47 * std::cos(p) - 0.03 * std::cos(2.9 * k), 2.47 * std::sin(p) + 0.025 * std::sin(1.3 * k), -0.03 * std::cos(0.7 * k)}});
  }
  System s = mol(a);
  const auto r = symmetrize_molecule(s, 0.15);
  EXPECT_EQ(r.group, "D6h");
  EXPECT_EQ(r.order, 24);
  EXPECT_LT(r.max_deviation, 1e-9);
  const double cc = norm(s.atoms[0].pos - s.atoms[2].pos);
  for (int k = 0; k < 6; ++k) EXPECT_NEAR(norm(s.atoms[size_t(2 * k)].pos - s.atoms[size_t((2 * k + 2) % 12)].pos), cc, 1e-9);
}

// Near-symmetric geometries: the group named, its order and its largest miss all come from one set of operations, every
// one of them within the tolerance. Noisy SF6 is not Oh at 0.1 Å (Oh misses by 0.14 Å) but a subgroup of matching order;
// noisy methane and ethane close into C3v and D3d, not into a group larger than their names.
TEST(PointGroup, NoisyGeometriesNameMatchesOrder) {
  auto order_of = [](const std::string& g) {
    if (g == "C1") return 1;
    if (g == "Cs" || g == "Ci") return 2;
    if (g == "T") return 12;
    if (g == "Td" || g == "Th" || g == "O") return 24;
    if (g == "Oh") return 48;
    if (g == "I") return 60;
    if (g == "Ih") return 120;
    const int n = std::stoi(g.substr(1));
    if (g[0] == 'S') return n;
    if (g[0] == 'C') return g.size() > 1 && (g.back() == 'v' || g.back() == 'h') ? 2 * n : n;
    return g.back() == 'h' || g.back() == 'd' ? 4 * n : 2 * n;   // Dn, Dnh, Dnd
  };
  struct Case { System s; double tol; std::string group; };
  const std::vector<Case> cases = {
      {mol({{16, {0.035326, 0.023879, -0.030708}}, {9, {1.568931, 0.003024, 0.045016}}, {9, {-1.532467, -0.025447, 0.057621}},
            {9, {-0.045832, 1.550175, 0.030857}}, {9, {-0.041762, -1.561324, -0.055295}}, {9, {0.020186, 0.031749, 1.568763}},
            {9, {0.045057, -0.022350, -1.536565}}}), 0.1, "D2"},
      {mol({{6, {0.019270, 0.049173, 0.033876}}, {1, {0.660017, 0.627364, 0.591423}}, {1, {-0.595304, -0.650098, 0.666099}},
            {1, {-0.573401, 0.617501, -0.641834}}, {1, {0.683616, -0.603024, -0.669600}}}), 0.1, "C3v"},
      {mol({{6, {-0.016813, -0.010222, 0.756432}}, {6, {0.039903, 0.007134, -0.776143}}, {1, {1.014244, -0.017988, 1.123861}},
            {1, {0.478137, 0.910120, -1.177150}}, {1, {-0.475153, 0.863292, 1.141258}}, {1, {-1.019123, -0.024812, -1.170132}},
            {1, {-0.473507, -0.852605, 1.184957}}, {1, {0.520472, -0.850272, -1.124744}}}), 0.2, "D3d"},
  };
  for (const auto& c : cases) {
    const auto r = point_group(c.s, c.tol);
    EXPECT_EQ(r.group, c.group);
    EXPECT_EQ(r.order, order_of(r.group)) << r.group;
    EXPECT_LE(r.max_deviation, c.tol) << r.group;
  }
  // benzene: C7 never passes for a six-ring at a loose tolerance
  std::vector<std::pair<int, Vec3>> benz;
  for (int k = 0; k < 6; ++k) {
    const double p = M_PI * k / 3;
    benz.push_back({6, {1.39 * std::cos(p) + 0.02 * std::sin(2.1 * k), 1.39 * std::sin(p), 0.03 * std::cos(1.3 * k)}});
    benz.push_back({1, {2.47 * std::cos(p), 2.47 * std::sin(p) - 0.03 * std::sin(0.7 * k), 0}});
  }
  for (double tol : {0.02, 0.05, 0.1, 0.5}) {
    const auto r = point_group(mol(benz), tol);
    EXPECT_EQ(r.order, order_of(r.group)) << r.group << " at " << tol;
    EXPECT_LE(r.max_deviation, tol) << r.group << " at " << tol;
    for (const auto& e : r.elements) EXPECT_LE(e.n, 6) << r.group << " at " << tol;
  }
}
