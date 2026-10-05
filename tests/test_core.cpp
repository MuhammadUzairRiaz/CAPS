#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <random>
#include <regex>
#include <fstream>
#include <set>

#include "caps/analysis.hpp"
#include "caps/appearance.hpp"
#include "caps/checks.hpp"
#include "caps/elements.hpp"
#include "caps/io.hpp"
#include "caps/render.hpp"
#include "caps/pipeline.hpp"
#include "caps/bundle.hpp"

using namespace caps;

static const std::string S = CAPS_SAMPLES;

namespace {
std::string tmp(const std::string& name) { return (std::filesystem::temp_directory_path() / name).string(); }
void write(const std::string& p, const std::string& text) { std::ofstream(p) << text; }
}  // namespace

TEST(Elements, SymbolsMassesNames) {
  EXPECT_EQ(element_from_symbol("c"), 6);
  EXPECT_EQ(element_from_symbol("Cl"), 17);
  EXPECT_EQ(element_from_mass(12.011), 6);
  EXPECT_EQ(element_from_mass(1.008), 1);
  EXPECT_EQ(element_from_mass(50.0), 0);   // nothing within 0.1
  EXPECT_EQ(element_from_name("CA"), 6);   // alpha carbon, not calcium
  EXPECT_EQ(element_from_name("Ca"), 20);
  EXPECT_EQ(element_from_name("HW1"), 1);
  EXPECT_EQ(element_from_name("NA"), 11);  // GROMACS sodium ion
  EXPECT_EQ(element_from_name("1HB"), 1);
}

TEST(Cell, MinimumImageAndWrapOrthorhombic) {
  Cell c;
  c.a = {10, 0, 0}; c.b = {0, 20, 0}; c.c = {0, 0, 30};
  const Vec3 d = c.minimum_image({9, -11, 16});
  EXPECT_NEAR(d[0], -1, 1e-12);
  EXPECT_NEAR(d[1], 9, 1e-12);
  EXPECT_NEAR(d[2], -14, 1e-12);
  const Vec3 w = c.wrap({-1, 25, 61});
  EXPECT_NEAR(w[0], 9, 1e-12);
  EXPECT_NEAR(w[1], 5, 1e-12);
  EXPECT_NEAR(w[2], 1, 1e-12);
  EXPECT_NEAR(c.volume(), 6000, 1e-9);
}

TEST(Cell, MinimumImageTriclinic) {
  Cell c;
  c.a = {10, 0, 0}; c.b = {6, 10, 0}; c.c = {0, 0, 10};
  // Brute force over images must agree with the cell's answer.
  const Vec3 d0{7.3, 8.1, -4.2};
  const Vec3 d = c.minimum_image(d0);
  double best = 1e300;
  for (int i = -2; i <= 2; ++i)
    for (int j = -2; j <= 2; ++j)
      for (int k = -2; k <= 2; ++k) best = std::min(best, norm(d0 + c.a * i + c.b * j + c.c * k));
  EXPECT_NEAR(norm(d), best, 1e-9);
}

TEST(LammpsData, SampleCountsTypesAndGeometry) {
  System s = read_lammps_data(S + "/ps_melt.data");
  EXPECT_EQ(s.atoms.size(), 1300u);
  EXPECT_EQ(s.bonds.size(), 1370u);
  EXPECT_TRUE(s.bonds_from_file);
  EXPECT_TRUE(s.has_mol);
  EXPECT_TRUE(s.unwrapped);   // image flags were applied
  int nm = 0;
  s.molecules(&nm);
  EXPECT_EQ(nm, 10);
  int c = 0, h = 0;
  for (const auto& a : s.atoms) { c += a.element == 6; h += a.element == 1; }
  EXPECT_EQ(c, 640);
  EXPECT_EQ(h, 660);
  // Builder geometry: every bond is C–C 1.39–1.53 Å or C–H 1.08–1.09 Å once unwrapped.
  for (const auto& b : s.bonds) {
    const double r = norm(s.atoms[b.i].pos - s.atoms[b.j].pos);
    EXPECT_GT(r, 1.07);
    EXPECT_LT(r, 1.54);
  }
  EXPECT_NEAR(s.cell.volume(), 33.0 * 33.0 * 33.0, 1e-6);
  EXPECT_NEAR(s.density(), 0.386, 0.001);   // loosely packed built cell (matches the design canvas)
  double q = 0;
  for (const auto& a : s.atoms) q += a.charge;
  EXPECT_NEAR(q, 0.0, 1300 * 5e-7);   // charges stored to 6 decimals: rounding bound
}

TEST(LammpsData, MissingAtomsInBondsAreReportedNotFatal) {
  const std::string p = tmp("caps_bad.data");
  write(p, "bad\n\n2 atoms\n2 bonds\n1 atom types\n\n0 10 xlo xhi\n0 10 ylo yhi\n0 10 zlo zhi\n\nMasses\n\n1 12.011\n\n"
           "Atoms # full\n\n1 1 1 0.0 1 1 1\n2 1 1 0.0 2.5 1 1\n\nBonds\n\n1 1 1 2\n2 1 2 3\n");
  System s = read_lammps_data(p);
  EXPECT_EQ(s.bonds.size(), 1u);
  bool noted = false;
  for (const auto& n : s.notes) noted |= n.find("not defined") != std::string::npos;
  EXPECT_TRUE(noted);
}

TEST(LammpsDump, FramesAndTopologyJoin) {
  System top = read_lammps_data(S + "/ps_melt.data");
  Trajectory t = read_lammps_dump(S + "/ps_melt.lammpstrj", &top);
  EXPECT_EQ(t.frames(), 3u);
  EXPECT_EQ(t.timesteps[2], 2000);
  EXPECT_EQ(t.topology.bonds.size(), 1370u);
  // The sample shifts frame k by 0.5·k Å in x.
  EXPECT_NEAR(t.positions[2][10][0] - t.positions[0][10][0], 1.0, 1e-4);
}

TEST(LammpsDump, ScaledCoordinatesAndTriclinicBounds) {
  const std::string p = tmp("caps_tric.lammpstrj");
  // Box with xy = 2: bounds written as the bounding box (xlo_bound = xlo + min(0,xy,...)).
  write(p, "ITEM: TIMESTEP\n5\nITEM: NUMBER OF ATOMS\n1\nITEM: BOX BOUNDS xy xz yz pp pp pp\n0 12 2\n0 10 0\n0 10 0\n"
           "ITEM: ATOMS id type xs ys zs\n1 1 0.5 0.5 0.5\n");
  Trajectory t = read_lammps_dump(p);
  const Cell& c = t.cells[0];
  EXPECT_NEAR(norm(c.a), 10, 1e-12);
  EXPECT_NEAR(c.b[0], 2, 1e-12);
  EXPECT_NEAR(t.positions[0][0][0], 6, 1e-12);   // 0.5·10 + 0.5·2
  EXPECT_NEAR(t.positions[0][0][1], 5, 1e-12);
}

TEST(Gro, UnitsResiduesAndBox) {
  Trajectory t = read_gro(S + "/ps_melt.gro");
  const System& s = t.topology;
  EXPECT_EQ(s.atoms.size(), 1300u);
  EXPECT_NEAR(norm(s.cell.a), 33.0, 1e-3);   // nm -> Å
  System d = read_lammps_data(S + "/ps_melt.data");
  EXPECT_NEAR(s.atoms[0].pos[0], d.atoms[0].pos[0], 0.006);   // .gro keeps 3 decimals in nm
  EXPECT_EQ(s.atoms[0].element, 6);
}

TEST(Gro, MoleculesFromPerceivedBonds) {
  Trajectory t = open_file(S + "/ps_melt.gro");
  int nm = 0;
  t.topology.molecules(&nm);
  EXPECT_EQ(nm, 10);
  EXPECT_EQ(t.topology.bonds.size(), 1370u);
}

TEST(Bonds, PerceivedMatchTopologyForBuilderGeometry) {
  System s = read_lammps_data(S + "/ps_melt.data");
  System w = s;
  for (auto& a : w.atoms) a.pos = w.cell.wrap(a.pos);   // perception must work across the periodic faces
  auto b = perceive_bonds(w);
  auto key = [](Bond x) { return std::make_pair(std::min(x.i, x.j), std::max(x.i, x.j)); };
  std::set<std::pair<uint32_t, uint32_t>> A, B;
  for (auto x : s.bonds) A.insert(key(x));
  for (auto x : b) B.insert(key(x));
  EXPECT_EQ(A, B);
}

TEST(Analysis, WholeMoleculesAndShapes) {
  System s = read_lammps_data(S + "/ps_melt.data");
  System w = s;
  for (auto& a : w.atoms) a.pos = w.cell.wrap(a.pos);
  make_molecules_whole(w);
  for (const auto& b : w.bonds) EXPECT_LT(norm(w.atoms[b.i].pos - w.atoms[b.j].pos), 1.6);
  const auto a = molecule_shapes(s), b = molecule_shapes(w);
  ASSERT_EQ(a.size(), 10u);
  for (size_t k = 0; k < a.size(); ++k) {
    EXPECT_NEAR(a[k].rg, b[k].rg, 1e-6);   // shape does not depend on which image is used
    EXPECT_GE(a[k].kappa2, 0.0);
    EXPECT_LE(a[k].kappa2, 1.0);
  }
}

TEST(Analysis, EigenOfKnownMatrix) {
  const double A[3][3] = {{2, 1, 0}, {1, 2, 0}, {0, 0, 5}};
  double w[3], V[3][3];
  symmetric_eigen3(A, w, V);
  EXPECT_NEAR(w[0], 1, 1e-10);
  EXPECT_NEAR(w[1], 3, 1e-10);
  EXPECT_NEAR(w[2], 5, 1e-10);
}

TEST(Analysis, IdealGasRdfIsNearOne) {
  // Uniform random points: g(r) ≈ 1 beyond the first bins.
  System s;
  s.cell.a = {20, 0, 0}; s.cell.b = {0, 20, 0}; s.cell.c = {0, 0, 20};
  uint64_t x = 88172645463325252ull;
  auto rnd = [&] { x ^= x << 13; x ^= x >> 7; x ^= x << 17; return (x >> 11) * (1.0 / 9007199254740992.0); };
  for (int i = 0; i < 2000; ++i) { Atom a; a.element = 18; a.pos = {20 * rnd(), 20 * rnd(), 20 * rnd()}; s.atoms.push_back(a); }
  const auto g = rdf(s, 0, 0, 8, 1.0, false);
  for (size_t k = 3; k < g.size(); ++k) EXPECT_NEAR(g[k].second, 1.0, 0.08);
}

TEST(Render, BackgroundsAndAlpha) {
  System s = read_lammps_data(S + "/ps_melt.data");
  Renderer r;
  Camera cam;
  RenderOptions o;
  o.width = 200; o.height = 150;
  o.background = Background::Transparent;
  Image t = r.render(s, cam, o);
  ASSERT_EQ(t.rgba.size(), 200u * 150 * 4);
  EXPECT_EQ(t.rgba[3], 0);   // corner pixel fully transparent
  int opaque = 0;
  for (size_t k = 3; k < t.rgba.size(); k += 4) opaque += t.rgba[k] == 255;
  EXPECT_GT(opaque, 2000);   // the structure covers part of the frame
  o.background = Background::White;
  Image w = r.render(s, cam, o);
  EXPECT_EQ(w.rgba[0], 255); EXPECT_EQ(w.rgba[1], 255); EXPECT_EQ(w.rgba[2], 255); EXPECT_EQ(w.rgba[3], 255);
  o.background = Background::Dark;
  Image d = r.render(s, cam, o);
  EXPECT_EQ(d.rgba[0], 0x0F); EXPECT_EQ(d.rgba[3], 255);
  // Picking returns an atom somewhere near the centre.
  int hits = 0;
  for (int y = 50; y < 100; ++y)
    for (int x = 70; x < 130; ++x) hits += r.pick(x, y) >= 0;
  EXPECT_GT(hits, 0);
}

TEST(Render, PngIsValid) {
  Image img;
  img.width = 3; img.height = 2;
  img.rgba.assign(3 * 2 * 4, 128);
  const auto png = encode_png(img);
  ASSERT_GT(png.size(), 40u);
  EXPECT_EQ(png[0], 0x89);
  EXPECT_EQ(png[1], 'P');
  EXPECT_EQ(std::string(png.end() - 8, png.end() - 4), "IEND");
}

TEST(Render, SvgHasNoBackgroundWhenTransparent) {
  System s = read_lammps_data(S + "/ps_melt.data");
  RenderOptions o;
  o.width = 300; o.height = 200;
  o.background = Background::Transparent;
  const std::string t = render_svg(s, Camera{}, o);
  EXPECT_EQ(t.find("<rect"), std::string::npos);
  o.background = Background::White;
  EXPECT_NE(render_svg(s, Camera{}, o).find("<rect width=\"100%\" height=\"100%\" fill=\"#FFFFFF\""), std::string::npos);
}

TEST(Io, DetectAndRoundTrip) {
  EXPECT_EQ(detect_format(S + "/ps_melt.data"), "lammps-data");
  EXPECT_EQ(detect_format(S + "/ps_melt.lammpstrj"), "lammps-dump");
  EXPECT_EQ(detect_format(S + "/ps_melt.gro"), "gro");
  System s = read_lammps_data(S + "/ps_melt.data");
  const std::string p = tmp("caps_rt.data");
  write_lammps_data(s, p);
  System r = read_lammps_data(p);
  ASSERT_EQ(r.atoms.size(), s.atoms.size());
  ASSERT_EQ(r.bonds.size(), s.bonds.size());
  for (size_t i = 0; i < s.atoms.size(); ++i)
    for (int k = 0; k < 3; ++k) EXPECT_NEAR(r.atoms[i].pos[k], s.atoms[i].pos[k], 1e-5);
  const std::string x = tmp("caps_rt.xyz");
  write_xyz(s, x);
  Trajectory tx = read_xyz(x);
  EXPECT_EQ(tx.topology.atoms.size(), 1300u);
  EXPECT_NEAR(norm(tx.topology.cell.a), 33.0, 1e-9);
}

TEST(Analysis, MeasureDistanceAngleDihedral) {
  System s;
  s.cell.a = {10, 0, 0}; s.cell.b = {0, 10, 0}; s.cell.c = {0, 0, 10};
  auto add = [&](Vec3 p) { Atom a; a.element = 6; a.pos = p; s.atoms.push_back(a); };
  add({9.5, 1, 1}); add({0.5, 1, 1}); add({0.5, 2, 1}); add({0.5, 2, 2});   // first pair straddles the x face
  EXPECT_NEAR(measure(s, {0, 1}), 1.0, 1e-12);
  EXPECT_NEAR(measure(s, {0, 1, 2}), 90.0, 1e-9);
  EXPECT_NEAR(std::fabs(measure(s, {0, 1, 2, 3})), 90.0, 1e-9);
  {   // IUPAC sign: A +x, B origin, C +z, D +y from C — the near bond turns clockwise onto the far one: +90°
    System q;
    for (Vec3 r : {Vec3{1, 0, 0}, Vec3{0, 0, 0}, Vec3{0, 0, 1}, Vec3{0, 1, 1}}) { Atom a; a.element = 6; a.pos = r; q.atoms.push_back(a); }
    EXPECT_NEAR(measure(q, {0, 1, 2, 3}), 90.0, 1e-9);
    EXPECT_NEAR(measure(q, {3, 2, 1, 0}), 90.0, 1e-9);   // the same dihedral read backwards keeps its sign
  }
  // trans chain: dihedral 180
  s.atoms[3].pos = {1.5, 2, 1};
  EXPECT_NEAR(std::fabs(measure(s, {0, 1, 2, 3})), 180.0, 1e-9);
}

TEST(Analysis, RdfCellListMatchesAllPairs) {
  System s = read_lammps_data(S + "/ps_melt.data");
  const auto fast = rdf(s, 6, 6, 12.0, 0.5, true);    // 12 < 33/2: cell list path
  // Brute force reference
  const auto mol = s.molecules();
  std::vector<double> h(fast.size(), 0);
  size_t na = 0;
  for (size_t i = 0; i < s.atoms.size(); ++i) {
    if (s.atoms[i].element != 6) continue;
    ++na;
    for (size_t j = 0; j < s.atoms.size(); ++j) {
      if (j == i || s.atoms[j].element != 6 || mol[i] == mol[j]) continue;
      const double r = norm(s.cell.minimum_image(s.atoms[j].pos - s.atoms[i].pos));
      if (r < 12.0) h[size_t(r / 0.5)] += 1;
    }
  }
  const double rho = na / s.cell.volume();
  for (size_t k = 0; k < fast.size(); ++k) {
    const double r0 = k * 0.5, r1 = r0 + 0.5;
    EXPECT_NEAR(fast[k].second, h[k] / (na * rho * 4.0 / 3.0 * M_PI * (r1 * r1 * r1 - r0 * r0 * r0)), 1e-12);
  }
}

TEST(Pdb, RoundTripWithConectAndCell) {
  System s = read_lammps_data(S + "/ps_melt.data");
  const std::string p = tmp("caps_rt.pdb");
  write_pdb(s, p);
  EXPECT_EQ(detect_format(p), "pdb");
  Trajectory t = open_file(p);
  const System& r = t.topology;
  ASSERT_EQ(r.atoms.size(), 1300u);
  EXPECT_EQ(r.bonds.size(), 1370u);
  EXPECT_TRUE(r.bonds_from_file);
  EXPECT_NEAR(norm(r.cell.a), 33.0, 1e-3);
  EXPECT_NEAR(r.cell.volume(), 33.0 * 33.0 * 33.0, 1.0);
  for (size_t i = 0; i < 1300; i += 97) {
    EXPECT_EQ(r.atoms[i].element, s.atoms[i].element);
    for (int k = 0; k < 3; ++k) EXPECT_NEAR(r.atoms[i].pos[k], s.atoms[i].pos[k], 6e-4);   // 3 decimals
  }
  int nm = 0;
  r.molecules(&nm);
  EXPECT_EQ(nm, 10);
}

TEST(Pdb, TriclinicCryst1) {
  const std::string p = tmp("caps_tric.pdb");
  write(p, "CRYST1   10.000   12.000   14.000  80.00  95.00 110.00 P 1           1\n"
           "ATOM      1  O   HOH A   1       1.000   2.000   3.000  1.00  0.00           O\nEND\n");
  Trajectory t = read_pdb(p);
  const Cell& c = t.topology.cell;
  EXPECT_NEAR(norm(c.a), 10, 1e-9);
  EXPECT_NEAR(norm(c.b), 12, 1e-9);
  EXPECT_NEAR(norm(c.c), 14, 1e-9);
  EXPECT_NEAR(std::acos(dot(c.b, c.c) / (12 * 14)) * 180 / M_PI, 80, 1e-9);
  EXPECT_NEAR(std::acos(dot(c.a, c.b) / (10 * 12)) * 180 / M_PI, 110, 1e-9);
  EXPECT_EQ(t.topology.atoms[0].element, 8);
}

// ---------------------------------------------------------------- Grow

#include "caps/grow.hpp"

namespace {

// Non-bonded pairs more than three bonds apart (or in different molecules), minimum image.
double worst_contact_margin(const System& s) {
  const auto nb = s.neighbours();
  auto lim = [](int a, int b) { return (a == 6 && b == 6) ? 3.0 : (a == 1 && b == 1) ? 2.0 : 2.45; };
  double worst = 9.0;
  for (size_t i = 0; i < s.atoms.size(); ++i) {
    // atoms within three bonds of i
    std::set<uint32_t> near{uint32_t(i)};
    std::vector<uint32_t> frontier{uint32_t(i)};
    for (int hop = 0; hop < 3; ++hop) {
      std::vector<uint32_t> next;
      for (auto x : frontier)
        for (auto y : nb[x])
          if (near.insert(y).second) next.push_back(y);
      frontier = next;
    }
    for (size_t j = i + 1; j < s.atoms.size(); ++j) {
      if (near.count(uint32_t(j))) continue;
      const double d = norm(s.cell.minimum_image(s.atoms[j].pos - s.atoms[i].pos));
      worst = std::min(worst, d - lim(s.atoms[i].element, s.atoms[j].element));
    }
  }
  return worst;
}

}  // namespace

TEST(Grow, SingleChainGeometryIsExact) {
  GrowOptions o;
  o.chains = 1; o.dp = 6; o.box = 60; o.seed = 3;
  System s = grow(o);
  ASSERT_EQ(s.atoms.size(), size_t(16 * 6 + 2));
  EXPECT_EQ(s.bonds.size(), s.atoms.size() - 1 + 6);   // a tree plus one ring closure per unit
  for (const auto& b : s.bonds) {
    const double r = norm(s.atoms[b.j].pos - s.atoms[b.i].pos);
    bool ok = false;
    for (double ref : {1.53, 1.51, 1.39, 1.08, 1.09}) ok |= std::fabs(r - ref) < 1e-6;
    EXPECT_TRUE(ok) << "bond length " << r;
  }
  for (int k = 1; k + 1 < 12; ++k) EXPECT_NEAR(measure(s, {uint32_t(k - 1), uint32_t(k), uint32_t(k + 1)}), 114.0, 1e-6);
  EXPECT_GE(worst_contact_margin(s), -0.05 - 1e-9);
}

TEST(Grow, CellMeetsDensityContactsAndCharge) {
  GrowOptions o;
  o.chains = 10; o.dp = 8; o.density = 0.4; o.seed = 1;
  GrowReport rep;
  System s = grow(o, &rep);
  EXPECT_EQ(s.atoms.size(), 1300u);
  EXPECT_EQ(s.bonds.size(), 1370u);
  EXPECT_NEAR(s.density(), 0.4, 1e-9);
  int nm = 0;
  s.molecules(&nm);
  EXPECT_EQ(nm, 10);
  const double w = worst_contact_margin(s);
  EXPECT_GE(w, -0.05 - 1e-9);                 // independent check, periodic images included
  EXPECT_NEAR(w, rep.worst_margin, 0.3);      // the report is in the same range (it sees the look-ahead atoms too)
  double q = 0;
  for (const auto& a : s.atoms) q += a.charge;
  EXPECT_NEAR(q, 0.0, 1e-9);                  // Gasteiger–Marsili conserves charge
  // Types: 1 c3, 2 ca, 3 hc, 4 ha — 16 c3 + 48 ca per DP-8 chain.
  int counts[5] = {0, 0, 0, 0, 0};
  for (const auto& a : s.atoms) counts[a.type]++;
  EXPECT_EQ(counts[1], 160);
  EXPECT_EQ(counts[2], 480);
  EXPECT_EQ(counts[3], 10 * (8 + 16 + 2));    // one H per CH, two per CH2, two end caps
  EXPECT_EQ(counts[4], 400);
}

TEST(Grow, SameSeedSameCell) {
  GrowOptions o;
  o.chains = 3; o.dp = 5; o.box = 30; o.seed = 42;
  System a = grow(o), b = grow(o);
  ASSERT_EQ(a.atoms.size(), b.atoms.size());
  for (size_t i = 0; i < a.atoms.size(); ++i)
    for (int k = 0; k < 3; ++k) EXPECT_EQ(a.atoms[i].pos[k], b.atoms[i].pos[k]);
}

TEST(Grow, TacticityControlsSideOfRings) {
  for (auto t : {Tacticity::Isotactic, Tacticity::Syndiotactic}) {
    GrowOptions o;
    o.chains = 1; o.dp = 10; o.box = 60; o.seed = 7; o.tacticity = t;
    System s = grow(o);
    // Ring side at each CH (odd backbone atom): sign of (prev − C) × (next − C) · (ipso − C); ipso atoms follow the backbone.
    std::vector<int> sides;
    uint32_t x = 20;   // first substituent index after 20 backbone atoms
    for (uint32_t k = 0; k < 20; ++k) {
      if (k % 2 == 1) {
        if (k + 1 < 20) {
          const Vec3 c = s.atoms[k].pos;
          const double v = dot(cross(s.atoms[k - 1].pos - c, s.atoms[k + 1].pos - c), s.atoms[x].pos - c);
          sides.push_back(v > 0 ? 1 : -1);
        }
        x += 12;
      } else {
        x += 2;
      }
    }
    for (size_t k = 1; k < sides.size(); ++k) {
      if (t == Tacticity::Isotactic) EXPECT_EQ(sides[k], sides[0]);
      else EXPECT_EQ(sides[k], -sides[k - 1]);
    }
  }
}

TEST(Grow, RejectsImpossibleRequestsClearly) {
  GrowOptions o;
  o.chains = 1; o.dp = 8; o.box = 8;
  EXPECT_THROW(grow(o), std::invalid_argument);
  o.box = 0; o.density = 0;
  EXPECT_THROW(grow(o), std::invalid_argument);
  EXPECT_NEAR(box_for_density([] { GrowOptions g; g.chains = 10; g.dp = 8; g.density = 0.386; return g; }()), 33.0, 0.01);
}

TEST(Render, ViewScaleForScaleBars) {
  System s;
  s.cell.a = {20, 0, 0}; s.cell.b = {0, 20, 0}; s.cell.c = {0, 0, 20};
  Atom a; a.element = 6; a.pos = {10, 10, 10};
  s.atoms.push_back(a);
  Camera cam; cam.yaw = 0; cam.pitch = 0; cam.zoom = 1; cam.perspective = false;
  RenderOptions o; o.width = 800; o.height = 800;
  EXPECT_NEAR(view_scale(s, cam, o), 800 * 0.45 / 11.0, 1e-9);   // half the cell + 1 Å pad
  RenderOptions big = o; big.width = 1600; big.height = 1600;
  EXPECT_NEAR(view_scale(s, cam, big), 2 * view_scale(s, cam, o), 1e-9);
  cam.zoom = 1.5;
  EXPECT_NEAR(view_scale(s, cam, o), 1.5 * 800 * 0.45 / 11.0, 1e-9);
}

TEST(Render, AmbientAccessibility) {
  // a 5 × 5 × 5 block of atoms 1.5 Å apart (as dense as bonded matter), plus one lone atom far away
  System s;
  for (int x = 0; x < 5; ++x)
    for (int y = 0; y < 5; ++y)
      for (int z = 0; z < 5; ++z) { Atom a; a.element = 6; a.pos = {1.5 * x, 1.5 * y, 1.5 * z}; s.atoms.push_back(a); }
  Atom lone; lone.element = 6; lone.pos = {60, 60, 60};
  s.atoms.push_back(lone);
  const std::vector<double> r(s.atoms.size(), 0.8);
  const std::vector<char> show(s.atoms.size(), 1);
  const auto acc = ambient_accessibility(s, r, show);
  const size_t centre = (2 * 5 + 2) * 5 + 2, corner = 0;
  EXPECT_FLOAT_EQ(acc.back(), 1.f);
  EXPECT_LT(acc[centre], 0.1f);
  EXPECT_GT(acc[corner], acc[centre] + 0.3f);
  // hidden atoms do not occlude
  std::vector<char> only_centre(s.atoms.size(), 0);
  only_centre[centre] = 1;
  EXPECT_FLOAT_EQ(ambient_accessibility(s, r, only_centre)[centre], 1.f);
  // and the image darkens
  Renderer R;
  RenderOptions o; o.width = 120; o.height = 120; o.depth_cue = false; o.outlines = false; o.show_cell = false;
  Camera cam;
  const Image plain = R.render(s, cam, o);
  o.ambient_occlusion = true;
  const Image ao = R.render(s, cam, o);
  long sp = 0, sa = 0;
  for (size_t k = 0; k < plain.rgba.size(); k += 4) { sp += plain.rgba[k]; sa += ao.rgba[k]; }
  EXPECT_LT(sa, sp);
}

TEST(Pipeline, StepsOnPolystyrene) {
  const Trajectory t = open_file(std::string(CAPS_SOURCE_DIR) + "/samples/ps_melt.lammpstrj", std::string(CAPS_SOURCE_DIR) + "/samples/ps_melt.data");
  const System f = t.frame(0);
  size_t nh = 0;
  for (const auto& a : f.atoms) nh += a.element == 1;
  auto run = [&](const std::string& json) { return run_pipeline(f, pipeline_from_json(Json::parse(json)), 0, 0); };

  // select hydrogens (listed top to bottom, run bottom to top), then delete them
  auto st = run(R"([{"type":"delete_selected"},{"type":"select_expression","expression":"Element == \"H\""}])");
  EXPECT_EQ(st.steps[1].summary, std::to_string(nh) + " of " + std::to_string(f.atoms.size()) + " selected");
  EXPECT_EQ(st.system.atoms.size(), f.atoms.size() - nh);
  EXPECT_EQ(st.attribute("Particles"), double(f.atoms.size() - nh));
  for (const auto& b : st.system.bonds) EXPECT_TRUE(st.system.atoms[b.i].element != 1 && st.system.atoms[b.j].element != 1);
  EXPECT_EQ(st.origin.size(), st.system.atoms.size());
  EXPECT_EQ(f.atoms[size_t(st.origin[5])].id, st.system.atoms[5].id);

  // clusters by bonds are the ten chains; C–H coordination is 66 H on 64 C per chain (C8H8 units, two end caps)
  st = run(R"([{"type":"coordination","cutoff":1.25,"element_a":6,"element_b":1},{"type":"cluster","mode":"bonds"}])");
  EXPECT_EQ(st.attribute("ClusterAnalysis.cluster_count"), 10.0);
  EXPECT_EQ(st.attribute("ClusterAnalysis.largest_size"), 130.0);
  EXPECT_NEAR(st.attribute("CoordinationAnalysis.mean"), 66.0 / 64.0, 1e-12);
  ASSERT_EQ(st.tables.size(), 2u);
  EXPECT_EQ(st.tables[0].name, "clusters");
  EXPECT_EQ(st.tables[1].name, "rdf");

  // colour coding by molecule is categorical with ten entries; by z it is continuous
  st = run(R"([{"type":"colour_coding","property":"Molecule"}])");
  EXPECT_FALSE(st.legend.continuous);
  EXPECT_EQ(st.legend.entries.size(), 10u);
  EXPECT_NE(st.colour[0], kNoColour);
  st = run(R"([{"type":"colour_coding","property":"Position.Z","map":"viridis"}])");
  EXPECT_TRUE(st.legend.continuous);
  EXPECT_LT(st.legend.lo, st.legend.hi);

  // slice keeps a 12 Å slab through the centre; replicate doubles the particles and the cell
  st = run(R"([{"type":"slice","normal":[0,0,1],"width":12}])");
  EXPECT_GT(st.system.atoms.size(), 0u);
  EXPECT_LT(st.system.atoms.size(), f.atoms.size());
  for (const auto& a : st.system.atoms) EXPECT_LE(std::fabs(a.pos[2] - (f.cell.origin[2] + f.cell.c[2] / 2)), 6.0 + 1e-9);
  st = run(R"([{"type":"replicate","nx":2,"ny":1,"nz":1}])");
  EXPECT_EQ(st.system.atoms.size(), 2 * f.atoms.size());
  EXPECT_EQ(st.system.bonds.size(), 2 * f.bonds.size());
  EXPECT_NEAR(st.system.cell.a[0], 2 * f.cell.a[0], 1e-9);
  EXPECT_NEAR(st.attribute("Density"), f.density(), 1e-9);

  // density profile averages to the cell density; histogram counts every particle
  st = run(R"([{"type":"histogram","property":"Charge","bins":20},{"type":"binning","axis":2,"bins":10,"reduction":"density"}])");
  double mean = 0;
  for (const auto& r : st.tables[0].rows) mean += r[1];
  EXPECT_NEAR(mean / 10, f.density(), 1e-6);
  double count = 0;
  for (const auto& r : st.tables[1].rows) count += r[1];
  EXPECT_EQ(count, double(f.atoms.size()));

  // computed properties feed later steps; errors stay on their step
  st = run(R"([{"type":"select_expression","expression":"Twice > 0 && Type != 3"},{"type":"compute_property","name":"Twice","expression":"2 * Charge"},{"type":"select_expression","expression":"Position.Q > 1"}])");
  EXPECT_EQ(st.steps[2].level, "error");
  EXPECT_NE(st.steps[2].summary.find("unknown property Position.Q"), std::string::npos);
  EXPECT_EQ(st.steps[1].level, "ok");
  EXPECT_EQ(st.props.at("Twice")[3], 2 * f.atoms[3].charge);
  EXPECT_GT(st.selected_count(), 0u);
  EXPECT_THROW(evaluate_expression(st, "1 +"), std::invalid_argument);
  EXPECT_THROW(evaluate_expression(st, "sqrt(2, 3)"), std::invalid_argument);
  const auto v = evaluate_expression(st, "max(Index, 3) % 2 == 1 || !(Index >= 0)");
  EXPECT_EQ(v[0], 1.0);   // max(0, 3) = 3, odd
  EXPECT_EQ(v[4], 0.0);

  // result JSON and particle table
  const Json r = pipeline_result_json(st);
  EXPECT_TRUE(r["attributes"].size() >= 10);
  const Json pj = particles_json(st, "Molecule == 1", 0, 5);
  EXPECT_EQ(pj["total"].number(), 130.0);
  EXPECT_EQ(pj["rows"].size(), 5u);
  EXPECT_EQ(pipeline_to_json(pipeline_from_json(Json::parse(R"([{"type":"wrap","enabled":false}])")))["steps"][size_t(0)]["enabled"].boolean(), false);
}

TEST(Pipeline, TopologyShapeAndFrames) {
  const Trajectory t = open_file(std::string(CAPS_SOURCE_DIR) + "/samples/ps_melt.lammpstrj", std::string(CAPS_SOURCE_DIR) + "/samples/ps_melt.data");
  auto run = [&](const std::string& json, int fr) { return run_pipeline(t.frame(size_t(fr)), pipeline_from_json(Json::parse(json)), fr, 0, &t); };
  // bonds perceived from distances match the data file's; every angle and dihedral is counted once
  auto st = run(R"([{"type":"topology"},{"type":"create_bonds","mode":"perceive","replace":true}])", 0);
  EXPECT_EQ(st.system.bonds.size(), t.topology.bonds.size());
  size_t angles = 0;
  for (const auto& nb : t.topology.neighbours()) angles += nb.size() * (nb.size() - 1) / 2;
  double counted = 0;
  for (const auto& r : st.tables[1].rows) counted += r[1];
  EXPECT_EQ(st.tables[1].name, "angles");
  EXPECT_EQ(counted, double(angles));
  EXPECT_NEAR(st.attribute("Topology.mean_bond"), 1.26, 0.05);
  // molecule shape: one row per chain, the per-atom Rg equals its chain's
  st = run(R"([{"type":"molecule_shape"}])", 0);
  ASSERT_EQ(st.tables[0].rows.size(), 10u);
  EXPECT_NEAR(st.props.at("MoleculeRg")[0], st.tables[0].rows[0][3], 1e-12);
  // displacements: nothing moves against itself; the sample's frames are shifted rigidly
  st = run(R"([{"type":"displacements","reference":"frame","frame":1}])", 1);
  EXPECT_NEAR(st.attribute("Displacements.msd"), 0.0, 1e-18);
  st = run(R"([{"type":"displacements","reference":"first"}])", 2);
  EXPECT_GT(st.attribute("Displacements.msd"), 0.0);
  EXPECT_EQ(st.props.at("Displacement").size(), t.topology.atoms.size());
  // smoothing over all three frames puts atom 0 at its mean position
  st = run(R"([{"type":"smooth","window":5}])", 1);
  Vec3 mean{0, 0, 0};
  for (size_t k = 0; k < 3; ++k) mean = mean + t.positions[k][0];
  mean = mean * (1.0 / 3);
  EXPECT_NEAR(norm(st.system.atoms[0].pos - mean), 0.0, 1e-9);
  EXPECT_EQ(st.attribute("Smoothed"), 1.0);
  // trailing: frames 0–1 only; per-atom properties made below are averaged too (the stored x over the window)
  st = run(R"([{"type":"smooth","window":5,"kind":"trailing","properties":true,"positions":false,"mark":false},
               {"type":"compute_property","name":"X","expression":"Position.X"}])", 1);
  const double x01 = 0.5 * (t.positions[0][0][0] + t.positions[1][0][0]);
  EXPECT_NEAR(st.props.at("X")[0], x01, 1e-9);
  EXPECT_NEAR(st.system.atoms[0].pos[0], t.positions[1][0][0], 1e-9);   // positions left alone
  EXPECT_EQ(st.attribute("Smoothed", 0), 0.0);
  EXPECT_NE(st.steps[0].summary.find("1 per-atom property averaged"), std::string::npos) << st.steps[0].summary;
  // manual selection: indices of the source frame, as text with ranges or an array; add and subtract
  st = run(R"([{"type":"manual_selection","atoms":"0 5 10-14"}])", 0);
  EXPECT_EQ(st.selected_count(), 7u);
  EXPECT_TRUE(st.selected[12] && !st.selected[6]);
  st = run(R"([{"type":"manual_selection","atoms":[11, 12],"mode":"subtract"},{"type":"manual_selection","atoms":"10-14"}])", 0);
  EXPECT_EQ(st.selected_count(), 3u);
  st = run(R"([{"type":"manual_selection","atoms":"0 x"}])", 0);
  EXPECT_EQ(st.steps[0].level, "error");
  // neighbour terms: Σ 1 over bonded neighbours is the degree; Σ Distance over them the summed bond lengths; the
  // neighbour's own properties are read (Σ Element over bonded neighbours of a CH₂ carbon: 6 + 6 + 1 + 1)
  st = run(R"([{"type":"compute_property","name":"Deg","expression":"0","neighbours":true,"neighbour_mode":"bonds","neighbour_expression":"1"},
               {"type":"compute_property","name":"Z","expression":"0","neighbours":true,"neighbour_mode":"bonds","neighbour_expression":"Element"},
               {"type":"compute_property","name":"L","expression":"0","neighbours":true,"neighbour_mode":"bonds","neighbour_expression":"Distance"}])", 0);
  {
    const auto nb = st.system.neighbours();
    for (size_t i : {size_t(0), size_t(5), size_t(40)}) {
      EXPECT_EQ(st.props.at("L").size(), st.system.atoms.size());
      double zsum = 0, lsum = 0;
      for (uint32_t j : nb[i]) {
        zsum += st.system.atoms[j].element;
        Vec3 d = st.system.atoms[j].pos - st.system.atoms[i].pos;
        lsum += norm(st.system.cell.minimum_image(d));
      }
      EXPECT_NEAR(st.props.at("Z")[i], zsum, 1e-9);
      EXPECT_NEAR(st.props.at("L")[i], lsum, 1e-9);
    }
  }
  // within a cutoff: every atom closer than 1.2 Å to a hydrogen-bearing carbon is one of its hydrogens
  st = run(R"([{"type":"compute_property","name":"H","expression":"0","neighbours":true,"cutoff":1.2,"neighbour_expression":"Element == 1"}])", 0);
  {
    const auto nb = st.system.neighbours();
    int hc = 0;
    for (uint32_t j : nb[0]) hc += st.system.atoms[j].element == 1;
    EXPECT_EQ(st.props.at("H")[0], double(hc));
  }
  // Gaussian weights favour the frame itself: between it and the plain mean
  st = run(R"([{"type":"smooth","window":5,"kind":"gaussian"}])", 1);
  EXPECT_LT(norm(st.system.atoms[0].pos - t.positions[1][0]), norm(mean - t.positions[1][0]) + 1e-12);
  // clusters of whole molecules joined through heavy atoms within 6 Å: the melt percolates
  st = run(R"([{"type":"cluster","mode":"cutoff","cutoff":6,"heavy_only":true,"unit":"molecules"}])", 0);
  EXPECT_EQ(st.attribute("ClusterAnalysis.cluster_count"), 1.0);
  // intermolecular coordination leaves out the bonded neighbours: within 1.8 Å only a few close contacts remain
  st = run(R"([{"type":"coordination","cutoff":1.8,"inter_only":true}])", 0);
  const double inter = st.attribute("CoordinationAnalysis.mean");
  st = run(R"([{"type":"coordination","cutoff":1.8}])", 0);
  EXPECT_GT(st.attribute("CoordinationAnalysis.mean"), 2.0);
  EXPECT_LT(inter, 0.02);
}

TEST(Pipeline, VectorsAndPaths) {
  const Trajectory t = open_file(std::string(CAPS_SOURCE_DIR) + "/samples/ps_melt.lammpstrj", std::string(CAPS_SOURCE_DIR) + "/samples/ps_melt.data");
  System whole = t.frame(0);
  if (!whole.unwrapped) make_molecules_whole(whole);
  auto run = [&](const std::string& json) { return run_pipeline(t.frame(0), pipeline_from_json(Json::parse(json)), 0, 0, &t); };
  // end-to-end arrows: one per chain, ⟨R²⟩ as the chain statistics give it
  auto st = run(R"([{"type":"vectors","property":"end_to_end"}])");
  EXPECT_EQ(st.segments.size(), 10u);
  EXPECT_TRUE(st.segments[0].arrow);
  EXPECT_NEAR(st.attribute("Vectors.mean_ree2"), internal_distances(whole).r2_end, 1e-6);
  // displacement arrows need the displacements step below them
  st = run(R"([{"type":"vectors","property":"displacement"}])");
  EXPECT_EQ(st.steps[0].level, "error");
  st = run(R"([{"type":"vectors","property":"displacement"},{"type":"displacements"}])");
  EXPECT_EQ(st.steps[0].level, "ok");
  // paths of the ten chain centres over the three frames: two segments each
  st = run(R"([{"type":"trajectory_lines"}])");
  EXPECT_EQ(st.segments.size(), 20u);
  EXPECT_FALSE(st.segments[0].arrow);
  // the renderer draws them
  Renderer R;
  RenderOptions o;
  o.width = 160; o.height = 120; o.show_cell = false;
  o.segments = run(R"([{"type":"vectors","property":"end_to_end","radius":0.6}])").segments;
  const Image with = R.render(t.frame(0), Camera{}, o);
  o.segments.clear();
  const Image without = R.render(t.frame(0), Camera{}, o);
  EXPECT_NE(with.rgba, without.rgba);
}

TEST(Pipeline, GridFields) {
  const Trajectory t = open_file(std::string(CAPS_SOURCE_DIR) + "/samples/ps_melt.lammpstrj", std::string(CAPS_SOURCE_DIR) + "/samples/ps_melt.data");
  const System f = t.frame(0);
  auto run = [&](const std::string& json) { return run_pipeline(f, pipeline_from_json(Json::parse(json)), 0, 0, &t); };
  // Voronoi cells fill the box; radical weights move volume from hydrogens (on the chains' surface, facing the empty
  // space of this loose cell) towards the larger carbons
  auto carbon = [&](const PipelineState& x) {
    double v = 0;
    for (size_t i = 0; i < f.atoms.size(); ++i) if (f.atoms[i].element == 6) v += x.props.at("AtomicVolume")[i];
    return v;
  };
  auto st = run(R"([{"type":"voronoi","method":"radical","grid":0.6}])");
  EXPECT_NEAR(st.attribute("Voronoi.sum"), f.cell.volume(), 1e-6 * f.cell.volume());
  const double radical = carbon(st);
  st = run(R"([{"type":"voronoi","method":"grid","grid":0.6}])");
  EXPECT_NEAR(st.attribute("Voronoi.sum"), f.cell.volume(), 1e-6 * f.cell.volume());
  EXPECT_GT(radical, carbon(st));
  // accessibility falls as the probe grows
  st = run(R"([{"type":"voids","probe":1.4,"grid":0.7,"show":false}])");
  const auto& sweep = st.tables[1].rows;
  for (size_t k = 1; k < sweep.size(); ++k) EXPECT_LE(sweep[k][1], sweep[k - 1][1]);
  EXPECT_GT(st.attribute("Voids.accessible_fraction"), 0.0);
  EXPECT_TRUE(st.segments.empty());
  // the smoothed density holds the cell's mass exactly; its profile averages to it
  st = run(R"([{"type":"density_field","grid":0.9,"sigma":1.5}])");
  EXPECT_NEAR(st.attribute("DensityField.mean"), f.density(), 1e-9);
  double m = 0;
  for (const auto& r : st.tables[0].rows) m += r[1];
  EXPECT_NEAR(m / st.tables[0].rows.size(), f.density(), 1e-9);
  EXPECT_FALSE(st.segments.empty());
  EXPECT_TRUE(st.meshes.empty());
  // the isosurface at the mean density: a closed surface between the dense and the sparse half of the cell
  st = run(R"([{"type":"density_field","grid":0.9,"sigma":1.5,"isosurface":true,"slice":false}])");
  ASSERT_EQ(st.meshes.size(), 1u);
  EXPECT_GT(st.meshes[0].mesh->triangles.size(), 100u);
  EXPECT_NEAR(st.attribute("DensityField.iso_level"), f.density(), 1e-9);
  const double share = st.attribute("DensityField.iso_volume_fraction");
  EXPECT_GT(share, 0.05);
  EXPECT_LT(share, 0.95);
  EXPECT_GT(st.attribute("DensityField.iso_area"), 0.0);
  EXPECT_TRUE(st.segments.empty());
  // no cell: an error on the step, not a crash
  System open = f;
  open.cell = Cell{};
  const auto e = run_pipeline(open, pipeline_from_json(Json::parse(R"([{"type":"voids"}])")));
  EXPECT_EQ(e.steps[0].level, "error");
}

TEST(Pipeline, MsdRecoversDiffusionAndScatter) {
  // each chain of the sample random-walks rigidly with D = 0.01 Å² per timestep for 400 frames
  const Trajectory t0 = open_file(std::string(CAPS_SOURCE_DIR) + "/samples/ps_melt.lammpstrj", std::string(CAPS_SOURCE_DIR) + "/samples/ps_melt.data");
  Trajectory t;
  t.topology = t0.topology;
  t.topology.unwrapped = true;
  const double D = 0.01;
  std::mt19937 rng(7);
  std::normal_distribution<double> step(0.0, std::sqrt(2 * D));
  const auto mol = t0.topology.molecules();
  std::vector<Vec3> shift(10, Vec3{0, 0, 0});
  for (int f = 0; f < 400; ++f) {
    std::vector<Vec3> pos = t0.positions[0];
    for (size_t i = 0; i < pos.size(); ++i) pos[i] = pos[i] + shift[size_t(mol[i])];
    t.positions.push_back(pos);
    t.cells.push_back(t0.cells[0]);
    t.timesteps.push_back(f);
    for (auto& s : shift) s = s + Vec3{step(rng), step(rng), step(rng)};
  }
  const auto st = run_pipeline(t.frame(0), pipeline_from_json(Json::parse(R"([{"type":"msd","timestep_fs":2}])")), 0, 0, &t);
  EXPECT_NEAR(st.attribute("MSD.D_centres"), D, 0.35 * D);
  EXPECT_NEAR(st.attribute("MSD.D_centres_cm2s"), st.attribute("MSD.D_centres") / 2 * 0.1, 1e-15);
  // rigid motion: every atom moves with its chain, so the atom and centre curves agree
  for (const auto& r : st.tables[0].rows) EXPECT_NEAR(r[1], r[2], 1e-6 * std::max(1.0, r[2]));
  // the same walk with every chain carried along 0.3 Å per frame: the drift swamps the MSD unless removed; the
  // per-molecule table has one curve per chain
  Trajectory td = t;
  for (size_t f = 0; f < td.positions.size(); ++f)
    for (auto& r : td.positions[f]) r = r + Vec3{0.3 * double(f), 0, 0};
  const auto sd = run_pipeline(td.frame(0), pipeline_from_json(Json::parse(R"([{"type":"msd"}])")), 0, 0, &td);
  const auto sr = run_pipeline(td.frame(0), pipeline_from_json(Json::parse(R"([{"type":"msd","remove_drift":true,"per_molecule":true}])")), 0, 0, &td);
  EXPECT_GT(sd.attribute("MSD.D_centres"), 20 * D);
  EXPECT_NEAR(sr.attribute("MSD.D_centres"), D * 0.9, 0.35 * D);   // the system's own walk (1/N of it) goes with the drift
  ASSERT_EQ(sr.tables.size(), 2u);
  EXPECT_EQ(sr.tables[1].columns.size(), 11u);
  // scatter: a property against itself is perfectly correlated; the table is drawn as points
  const auto sc = run_pipeline(t.frame(0), pipeline_from_json(Json::parse(R"([{"type":"scatter","x":"Position.X","y":"Position.X"}])")));
  EXPECT_NEAR(sc.attribute("Scatter.pearson_r"), 1.0, 1e-9);
  EXPECT_TRUE(sc.tables[0].points);
}

// Cluster unwrap: the melt's chains, clustered by bonds and unwrapped, have every bond short in the shown positions
TEST(Pipeline, ClusterUnwrapMakesChainsWhole) {
  const Trajectory t = open_file(std::string(CAPS_SOURCE_DIR) + "/samples/ps_melt.lammpstrj", std::string(CAPS_SOURCE_DIR) + "/samples/ps_melt.data");
  System w = t.frame(0);
  for (auto& a : w.atoms) a.pos = w.cell.wrap(a.pos);   // folded into the cell: chains cut at the faces
  auto longest = [](const System& s) {
    double m = 0;
    for (const auto& b : s.bonds) m = std::max(m, norm(s.atoms[b.j].pos - s.atoms[b.i].pos));
    return m;
  };
  ASSERT_GT(longest(w), 5.0);
  const auto st = run_pipeline(w, pipeline_from_json(Json::parse(R"([{"type":"cluster","mode":"bonds","unwrap":true}])")));
  EXPECT_LT(longest(st.system), 1.8);
  EXPECT_EQ(st.tables[0].rows.size(), 10u);
}

TEST(Bundle, Sha256ZipAndReproduce) {
  EXPECT_EQ(sha256_hex("abc"), "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad");
  EXPECT_EQ(sha256_hex(""), "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855");
  EXPECT_EQ(sha256_hex(std::string(1000, 'a')), "41edece42d63e8d9bf515a9ba6932e1c20cbc9f5a5d134645adb5db1b9737ea3");
  const auto dir = std::filesystem::temp_directory_path();
  const std::string zip = (dir / "caps_test.zip").string();
  write_zip(zip, {{"a.txt", "hello"}, {"data/b.csv", std::string(5000, 'x') + "end"}});
  const auto back = read_zip(zip);
  EXPECT_EQ(back.at("a.txt"), "hello");
  EXPECT_EQ(back.at("data/b.csv").size(), 5003u);
  // a bundle with its input and pipeline reproduces its data files
  const std::string in = std::string(CAPS_SOURCE_DIR) + "/samples/ps_melt.lammpstrj", top = std::string(CAPS_SOURCE_DIR) + "/samples/ps_melt.data";
  const Trajectory t = open_file(in, top);
  const Pipeline p = pipeline_from_json(Json::parse(R"([{"type":"molecule_shape"},{"type":"coordination","cutoff":5,"bins":50}])"));
  BundleOptions o;
  o.name = "test";
  o.input = in;
  o.topology = top;
  o.include_input = true;
  o.width = 160;
  o.height = 90;
  const auto files = bundle_files(t, p, o, Camera{}, RenderOptions{});
  std::vector<std::string> names;
  for (const auto& f : files) names.push_back(f.name);
  for (const char* want : {"figure.png", "figure.svg", "data/rdf.csv", "data/molecules.csv", "pipeline.json", "provenance.json", "README.txt"})
    EXPECT_NE(std::find(names.begin(), names.end(), want), names.end()) << want;
  const std::string bundle = (dir / "caps_test.caps-bundle.zip").string();
  write_bundle(bundle, files);
  std::vector<std::string> report;
  EXPECT_TRUE(reproduce_bundle(bundle, report));
  EXPECT_EQ(report.size(), 2u);
  // tampering with a data file is caught
  auto z = read_zip(bundle);
  z["provenance.json"] = std::regex_replace(z["provenance.json"], std::regex("(\"data/rdf.csv\": \")[0-9a-f]"), "$1z");
  std::vector<std::pair<std::string, std::string>> items(z.begin(), z.end());
  write_zip(bundle, items);
  report.clear();
  EXPECT_FALSE(reproduce_bundle(bundle, report));
  std::filesystem::remove(zip);
  std::filesystem::remove(bundle);
}

TEST(Io, InspectBeforeOpening) {
  const std::string dump = std::string(CAPS_SOURCE_DIR) + "/samples/ps_melt.lammpstrj", data = std::string(CAPS_SOURCE_DIR) + "/samples/ps_melt.data";
  const auto r = inspect_file(dump, data);
  EXPECT_EQ(r.format, "lammps-dump");
  EXPECT_EQ(r.atoms, 1300u);
  EXPECT_EQ(r.frames, 3u);
  ASSERT_FALSE(r.columns.empty());
  EXPECT_EQ(r.columns[0].maps_to, "Particle Identifier");
  EXPECT_EQ(r.types.size(), 4u);
  EXPECT_NE(r.bonds_from.find("1370 bonds"), std::string::npos);
  EXPECT_EQ(r.head[0], "ITEM: TIMESTEP");
  const auto g = inspect_file(std::string(CAPS_SOURCE_DIR) + "/samples/ps_melt.gro");
  EXPECT_EQ(g.format, "gro");
  EXPECT_EQ(g.atoms, 1300u);
  EXPECT_EQ(g.frames, 1u);
  EXPECT_FALSE(g.types.empty());
}

TEST(Pipeline, YamlRoundTrip) {
  const Pipeline p = pipeline_from_json(Json::parse(R"([
    {"type":"colour_coding","property":"DistanceToCOM","map":"viridis"},
    {"type":"select_expression","expression":"Type == 2 && Position.Z > 13","enabled":false},
    {"type":"slice","normal":[0,0,1],"width":12.5,"invert":true},
    {"type":"unwrap"}])"));
  const std::string y = pipeline_to_yaml(p, "PS melt · structure report", "PS_melt.lammpstrj", "PS_melt.data");
  EXPECT_NE(y.find("  - unwrap: {}"), std::string::npos);
  EXPECT_LT(y.find("unwrap"), y.find("colour_coding"));   // written in the order they run
  std::string name, file, topo;
  const Pipeline q = pipeline_from_yaml(y, &name, &file, &topo);
  EXPECT_EQ(name, "PS melt · structure report");
  EXPECT_EQ(file, "PS_melt.lammpstrj");
  EXPECT_EQ(topo, "PS_melt.data");
  EXPECT_EQ(pipeline_to_json(q).dump(0), pipeline_to_json(p).dump(0));
  EXPECT_THROW(pipeline_from_yaml("steps:\n  - wrap: {}\n"), std::invalid_argument);
  EXPECT_THROW(pipeline_from_yaml("caps_pipeline: 1\nsteps:\n  - wrap: {a: [1, 2}\n"), std::invalid_argument);
}

TEST(Pipeline, BranchesShareTheTrunk) {
  const Trajectory t = open_file(std::string(CAPS_SOURCE_DIR) + "/samples/ps_melt.lammpstrj", std::string(CAPS_SOURCE_DIR) + "/samples/ps_melt.data");
  Pipeline p = pipeline_from_json(Json::parse(R"([
    {"type":"colour_coding","property":"Molecule","branch":"Figure"},
    {"type":"compute_property","name":"X","expression":"Position.X","branch":"Numbers"},
    {"type":"coordination","cutoff":3.0,"rmax":8}])"));
  EXPECT_EQ(pipeline_trunk(p), 2u);
  p.branch = "Numbers";
  auto st = run_pipeline(t.frame(0), p, 0, 0, &t);
  EXPECT_TRUE(st.props.count("X"));
  EXPECT_EQ(st.steps[0].level, "off");
  EXPECT_NE(st.steps[0].summary.find("branch Figure"), std::string::npos);
  EXPECT_GT(st.attribute("CoordinationAnalysis.mean"), 0.0);
  p.branch = "";
  st = run_pipeline(t.frame(0), p, 0, 0, &t);
  EXPECT_FALSE(st.props.count("X"));
  EXPECT_EQ(st.steps[1].level, "off");
  // the trunk run once, each branch on a copy of it: the same result as the whole run
  p.branch = "Numbers";
  PipelineState trunk = pipeline_begin(t.frame(0), p, 0, 0, &t);
  pipeline_run_steps(trunk, p, 3, 2);
  PipelineState a = trunk;
  pipeline_run_steps(a, p, 2, 0);
  pipeline_finish(a);
  const auto whole = run_pipeline(t.frame(0), p, 0, 0, &t);
  ASSERT_EQ(a.attributes.size(), whole.attributes.size());
  for (size_t k = 0; k < a.attributes.size(); ++k) EXPECT_EQ(a.attributes[k], whole.attributes[k]);
  EXPECT_EQ(a.props.at("X"), whole.props.at("X"));
  EXPECT_EQ(pipeline_from_json(pipeline_to_json(p)).branch, "Numbers");
  EXPECT_EQ(pipeline_from_yaml(pipeline_to_yaml(p)).branch, "Numbers");
  EXPECT_EQ(pipeline_to_json(pipeline_from_yaml(pipeline_to_yaml(p))).dump(0), pipeline_to_json(p).dump(0));
}

TEST(Pipeline, OutputsBlock) {
  // the outputs: block survives YAML and JSON, and writes the files it names
  const std::string y = "caps_pipeline: 1\nsteps:\n  - coordination: {cutoff: 3.0, rmax: 8}\noutputs:\n  - table: rdf -> rdf.csv\n"
                        "  - plot: rdf -> plots/rdf.svg\n  - attributes: -> attributes.csv\n  - render: {file: view.png, size: [320, 240]}\n"
                        "  - table: missing -> missing.csv\n";
  const Pipeline p = pipeline_from_yaml(y);
  ASSERT_EQ(p.outputs.size(), 5u);
  EXPECT_EQ(p.outputs[1].kind, "plot");
  EXPECT_EQ(p.outputs[1].what, "rdf");
  EXPECT_EQ(p.outputs[1].path, "plots/rdf.svg");
  EXPECT_EQ(p.outputs[3].width, 320);
  const Pipeline q = pipeline_from_yaml(pipeline_to_yaml(p));
  EXPECT_EQ(pipeline_to_json(q).dump(0), pipeline_to_json(p).dump(0));
  EXPECT_EQ(pipeline_to_json(pipeline_from_json(pipeline_to_json(p))).dump(0), pipeline_to_json(p).dump(0));
  const Trajectory t = open_file(std::string(CAPS_SOURCE_DIR) + "/samples/ps_melt.lammpstrj", std::string(CAPS_SOURCE_DIR) + "/samples/ps_melt.data");
  const auto st = run_pipeline(t.frame(0), p, 0, 0, &t);
  const auto dir = std::filesystem::temp_directory_path() / "caps_outputs_test";
  std::filesystem::remove_all(dir);
  const auto log = write_pipeline_outputs(st, p, dir.string());
  ASSERT_EQ(log.size(), 5u);
  EXPECT_EQ(log[4].rfind("not written: missing.csv", 0), 0u) << log[4];
  EXPECT_TRUE(std::filesystem::file_size(dir / "rdf.csv") > 100);
  std::ifstream svg(dir / "plots" / "rdf.svg");
  const std::string sv((std::istreambuf_iterator<char>(svg)), {});
  EXPECT_NE(sv.find("<path d=\"M"), std::string::npos);
  EXPECT_TRUE(std::filesystem::file_size(dir / "view.png") > 1000);
  std::filesystem::remove_all(dir);
  EXPECT_THROW(pipeline_from_yaml("caps_pipeline: 1\noutputs:\n  - table: rdf\n"), std::invalid_argument);
}

TEST(Pipeline, PythonStep) {
#ifdef _WIN32
  const char* python = "python";
#else
  const char* python = "python3";
#endif
  if (std::system((std::string(python) + " --version > " + (std::filesystem::temp_directory_path() / "caps_py.txt").string() + " 2>&1").c_str()) != 0)
    GTEST_SKIP() << "no Python";
  const auto dir = std::filesystem::temp_directory_path();
  const std::string script = (dir / "caps_test_step.py").string();
  {
    std::ofstream f(script);
    f << "from caps.pipeline import step\n\n"
         "@step(name=\"Heavy atoms\")\n"
         "def modify(frame, data):\n"
         "    el = data.particles[\"Element\"]\n"
         "    heavy = [1 if e != \"H\" else 0 for e in el]\n"
         "    data.attributes[\"HeavyCount\"] = sum(heavy)\n"
         "    data.particles[\"Heavy\"] = heavy\n"
         "    data.tables[\"per_molecule\"] = {1: 64, 2: 64}\n"
         "    data.selection = heavy\n";
  }
  const Trajectory t = open_file(std::string(CAPS_SOURCE_DIR) + "/samples/ps_melt.lammpstrj", std::string(CAPS_SOURCE_DIR) + "/samples/ps_melt.data");
  Json step = Json::object();
  step["type"] = "python";
  step["file"] = script;
  step["path"] = std::string(CAPS_SOURCE_DIR) + "/data/python";
  Json arr = Json::array();
  arr.push_back(step);
  const auto st = run_pipeline(t.frame(0), pipeline_from_json(arr), 0, 0, &t);
  ASSERT_EQ(st.steps[0].level, "ok") << st.steps[0].summary;
  EXPECT_EQ(st.steps[0].title, "Heavy atoms");
  EXPECT_EQ(st.attribute("HeavyCount"), 640.0);
  EXPECT_EQ(st.props.at("Heavy")[0], 1.0);
  EXPECT_EQ(st.selected_count(), 640u);
  ASSERT_EQ(st.tables.size(), 1u);
  EXPECT_EQ(st.tables[0].rows.size(), 2u);
  // a failing script reports its error on the step
  { std::ofstream f(script); f << "from caps.pipeline import step\n@step()\ndef modify(frame, data):\n    raise ValueError('bad input')\n"; }
  const auto e = run_pipeline(t.frame(0), pipeline_from_json(arr), 0, 0, &t);
  EXPECT_EQ(e.steps[0].level, "error");
  EXPECT_NE(e.steps[0].summary.find("ValueError: bad input"), std::string::npos) << e.steps[0].summary;
  std::filesystem::remove(script);
  // the step typed in (code instead of a file)
  Json typed = Json::object();
  typed["type"] = "python";
  typed["code"] = std::string("from caps.pipeline import step\n\n@step(name=\"Typed\")\ndef modify(frame, data):\n    data.attributes[\"N\"] = len(data.particles[\"Element\"])\n");
  typed["path"] = std::string(CAPS_SOURCE_DIR) + "/data/python";
  Json arr2 = Json::array();
  arr2.push_back(typed);
  const auto c = run_pipeline(t.frame(0), pipeline_from_json(arr2), 0, 0, &t);
  ASSERT_EQ(c.steps[0].level, "ok") << c.steps[0].summary;
  EXPECT_EQ(c.steps[0].title, "Typed");
  EXPECT_EQ(c.attribute("N"), 1300.0);
  // what it prints is its console; a failure keeps the traceback there
  typed["code"] = std::string("from caps.pipeline import step\n\n@step()\ndef modify(frame, data):\n    print('hello from the step')\n    raise KeyError('Nope')\n");
  Json arr3 = Json::array();
  arr3.push_back(typed);
  const auto f = run_pipeline(t.frame(0), pipeline_from_json(arr3), 0, 0, &t);
  EXPECT_EQ(f.steps[0].level, "error");
  EXPECT_NE(f.steps[0].output.find("hello from the step"), std::string::npos) << f.steps[0].output;
  EXPECT_NE(f.steps[0].output.find("KeyError"), std::string::npos) << f.steps[0].output;
}

TEST(Io, FileWithoutAtomsIsAnError) {
  const std::string path = (std::filesystem::temp_directory_path() / "caps_test_garbage.data").string();
  { std::ofstream f(path); f << "garbage\n"; }
  EXPECT_THROW(open_file(path), ReadError);
  std::filesystem::remove(path);
}

TEST(Io, GroWriterRoundTrip) {
  const Trajectory t = open_file(std::string(CAPS_SOURCE_DIR) + "/samples/ps_melt.lammpstrj", std::string(CAPS_SOURCE_DIR) + "/samples/ps_melt.data");
  const System s = t.frame(0);
  const std::string path = (std::filesystem::temp_directory_path() / "caps_test_roundtrip.gro").string();
  write_gro(s, path);
  const Trajectory g = read_gro(path);
  ASSERT_EQ(g.topology.atoms.size(), s.atoms.size());
  for (size_t i = 0; i < s.atoms.size(); i += 97)
    for (int k = 0; k < 3; ++k) EXPECT_NEAR(g.positions[0][i][k], s.atoms[i].pos[k], 0.006);   // 0.001 nm precision
  EXPECT_NEAR(g.topology.cell.a[0], s.cell.a[0], 1e-3);
  std::filesystem::remove(path);
}

TEST(Io, StagedOpen) {
  const std::string dump = std::string(CAPS_SOURCE_DIR) + "/samples/ps_melt.lammpstrj", data = std::string(CAPS_SOURCE_DIR) + "/samples/ps_melt.data";
  const Trajectory full = open_file(dump, data);
  ASSERT_GE(full.frames(), 3u);
  std::vector<int> stages;
  double last = 0;
  OpenProgress p;
  p.report = [&](int st, double f, const std::string&) { if (stages.empty() || stages.back() != st) stages.push_back(st); if (st == 3) last = f; return true; };
  const Trajectory all = open_file(dump, data, p);
  EXPECT_EQ(all.frames(), full.frames());
  EXPECT_EQ(stages, (std::vector<int>{0, 1, 2, 3}));
  EXPECT_DOUBLE_EQ(last, 1.0);
  // frame 0 alone, with the same bonds
  p.max_frames = 1;
  const Trajectory first = open_file(dump, data, p);
  EXPECT_EQ(first.frames(), 1u);
  EXPECT_EQ(first.topology.bonds.size(), full.topology.bonds.size());
  EXPECT_EQ(first.positions[0][7][2], full.positions[0][7][2]);
  // stopped after the second frame: the frames read are kept, with a note
  p.max_frames = 0;
  int calls = 0;
  p.report = [&](int st, double, const std::string&) { return st != 3 || ++calls < 2; };
  const Trajectory part = open_file(dump, data, p);
  EXPECT_EQ(part.frames(), 2u);
  EXPECT_NE(std::find_if(part.topology.notes.begin(), part.topology.notes.end(), [](const std::string& n) { return n.find("stopped") != std::string::npos; }),
            part.topology.notes.end());
}

TEST(FileChecks, SampleAndBrokenStructures) {
  const Trajectory t = open_file(std::string(CAPS_SOURCE_DIR) + "/samples/ps_melt.lammpstrj", std::string(CAPS_SOURCE_DIR) + "/samples/ps_melt.data");
  const auto c = file_checks(t);
  auto has = [&](const std::string& level, const std::string& title) {
    return std::any_of(c.begin(), c.end(), [&](const FileCheck& x) { return x.level == level && x.title.find(title) != std::string::npos; });
  };
  EXPECT_TRUE(has("pass", "Atom counts agree"));
  EXPECT_TRUE(has("pass", "Bonds are plausible"));
  EXPECT_TRUE(has("pass", "Charges are neutral"));
  EXPECT_TRUE(has("note", "outside the box"));
  // two atoms on top of each other, a stretched bond, a net charge
  Trajectory b;
  System& s = b.topology;
  s.cell.a = {20, 0, 0};
  s.cell.b = {0, 20, 0};
  s.cell.c = {0, 0, 20};
  for (int k = 0; k < 4; ++k) {
    Atom a;
    a.id = k + 1;
    a.element = 6;
    a.charge = 0.25;
    a.pos = {5.0 + (k == 3 ? 0.3 : 3.0 * k), 5, 5};
    s.atoms.push_back(a);
  }
  s.has_charges = true;
  s.bonds = {{0, 1, 1}};   // 3.0 Å: stretched
  std::vector<Vec3> p;
  for (const auto& a : s.atoms) p.push_back(a.pos);
  b.positions.push_back(p);
  b.cells.push_back(s.cell);
  b.timesteps.push_back(0);
  const auto d = file_checks(b);
  auto dhas = [&](const std::string& level, const std::string& title) {
    return std::any_of(d.begin(), d.end(), [&](const FileCheck& x) { return x.level == level && x.title.find(title) != std::string::npos; });
  };
  EXPECT_TRUE(dhas("warn", "stretched"));
  EXPECT_TRUE(dhas("warn", "closer than 0.7"));
  EXPECT_TRUE(dhas("warn", "Net charge"));
}

// Horn superposition: a rotated, translated copy fits back exactly; noise shows as the RMSD; weights pick the atoms
#include "caps/superpose.hpp"
TEST(Superpose, RecoversRotationAndTranslation) {
  std::vector<caps::Vec3> ref;
  uint64_t seed = 12345;
  auto rnd = [&] { seed = seed * 6364136223846793005ULL + 1442695040888963407ULL; return double(seed >> 11) / double(1ULL << 53) * 10 - 5; };
  for (int i = 0; i < 40; ++i) ref.push_back({rnd(), rnd(), rnd()});
  // rotation about (1, 2, 3) by 1.1 rad, then a shift
  const double th = 1.1, l = std::sqrt(14.0), ux = 1 / l, uy = 2 / l, uz = 3 / l, c = std::cos(th), s = std::sin(th);
  const double R[3][3] = {{c + ux * ux * (1 - c), ux * uy * (1 - c) - uz * s, ux * uz * (1 - c) + uy * s},
                          {uy * ux * (1 - c) + uz * s, c + uy * uy * (1 - c), uy * uz * (1 - c) - ux * s},
                          {uz * ux * (1 - c) - uy * s, uz * uy * (1 - c) + ux * s, c + uz * uz * (1 - c)}};
  std::vector<caps::Vec3> mov;
  for (const auto& p : ref)
    mov.push_back({R[0][0] * p[0] + R[0][1] * p[1] + R[0][2] * p[2] + 7, R[1][0] * p[0] + R[1][1] * p[1] + R[1][2] * p[2] - 3,
                   R[2][0] * p[0] + R[2][1] * p[1] + R[2][2] * p[2] + 11});
  const auto fit = caps::superpose(ref, mov);
  EXPECT_LT(fit.rmsd, 1e-9);
  EXPECT_EQ(fit.fitted, 40u);
  for (size_t i = 0; i < ref.size(); ++i) EXPECT_LT(caps::norm(fit.apply(mov[i]) - ref[i]), 1e-9);
  // one atom moved by 3 Å: fitting on the others leaves exactly that shift on it
  mov[5] = mov[5] + caps::Vec3{3, 0, 0};
  std::vector<double> w(ref.size(), 1.0);
  w[5] = 0;
  const auto fit2 = caps::superpose(ref, mov, w);
  EXPECT_LT(fit2.rmsd, 1e-9);
  EXPECT_NEAR(caps::norm(fit2.apply(mov[5]) - ref[5]), 3.0, 1e-9);
  EXPECT_THROW(caps::superpose(ref, std::vector<caps::Vec3>(3)), std::invalid_argument);
}

// Wrap by molecule keeps every molecule whole with its centre of mass in the cell; unwrap nojump follows an atom that
// crosses the boundary each frame; unwrap along bonds leaves no split bond
TEST(Pipeline, WrapModesAndUnwrapMethods) {
  const Trajectory t = open_file(std::string(CAPS_SOURCE_DIR) + "/samples/ps_melt.lammpstrj", std::string(CAPS_SOURCE_DIR) + "/samples/ps_melt.data");
  const System f = t.frame(0);
  auto run = [&](const System& s, const std::string& json, int frame = 0, const Trajectory* tr = nullptr) {
    return run_pipeline(s, pipeline_from_json(Json::parse(json)), frame, 0, tr);
  };
  auto st = run(f, R"([{"type":"wrap","mode":"molecules"}])");
  EXPECT_EQ(st.attribute("Wrap.bonds_across_faces"), 0.0);
  int nm = 0;
  const auto mol = st.system.molecules(&nm);
  std::vector<Vec3> com(size_t(nm), Vec3{0, 0, 0});
  std::vector<double> mass(size_t(nm), 0);
  for (size_t i = 0; i < mol.size(); ++i) {
    const double m = st.system.mass_of(st.system.atoms[i]);
    com[size_t(mol[i])] = com[size_t(mol[i])] + st.system.atoms[i].pos * m, mass[size_t(mol[i])] += m;
  }
  for (int m = 0; m < nm; ++m) {
    const Vec3 fr = st.system.cell.to_fractional(com[size_t(m)] * (1.0 / mass[size_t(m)]));
    for (int k = 0; k < 3; ++k) EXPECT_TRUE(fr[k] >= -1e-9 && fr[k] < 1 + 1e-9) << m;
  }
  // atoms folded one by one cut bonds; unwrapping along bonds joins them again
  st = run(f, R"([{"type":"unwrap"},{"type":"wrap"}])");
  EXPECT_EQ(st.attribute("Unwrap.bonds_split"), 0.0);
  EXPECT_GT(st.attribute("Unwrap.atoms_moved"), 0.0);
  ASSERT_FALSE(st.tables.empty());
  EXPECT_EQ(st.tables.back().name, "images");

  // one argon atom stepping +4 Å a frame through a 10 Å box, written wrapped
  const auto dir = std::filesystem::temp_directory_path() / "caps_nojump";
  std::filesystem::create_directories(dir);
  const auto gro = (dir / "hop.gro").string();
  {
    std::ofstream o(gro);
    for (int k = 0; k < 5; ++k) {
      const double x = std::fmod(1.0 + 4.0 * k, 10.0) / 10.0;
      char line[128];
      std::snprintf(line, sizeof line, "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n", 1, "AR", "AR", 1, x, 0.5, 0.5);
      o << "hop t= " << k << "\n    1\n" << line << "   1.00000   1.00000   1.00000\n";
    }
  }
  const Trajectory hop = open_file(gro);
  ASSERT_EQ(hop.frames(), 5u);
  const auto nj = run(hop.frame(4), R"([{"type":"unwrap","method":"nojump"}])", 4, &hop);
  EXPECT_NEAR(nj.system.atoms[0].pos[0], 1.0 + 16.0, 1e-6);
  const auto plain = run(hop.frame(4), R"([{"type":"unwrap","method":"bonds"}])", 4, &hop);
  EXPECT_NEAR(plain.system.atoms[0].pos[0], 7.0, 1e-6);
}

// Freeze property: frame-0 values on a later frame, matched by identifier, including a property made by a step below
TEST(Pipeline, FreezePropertyAtAReferenceFrame) {
  const Trajectory t = open_file(std::string(CAPS_SOURCE_DIR) + "/samples/ps_melt.lammpstrj", std::string(CAPS_SOURCE_DIR) + "/samples/ps_melt.data");
  ASSERT_GE(t.frames(), 2u);
  const int last = int(t.frames()) - 1;
  const System f0 = t.frame(0), fl = t.frame(size_t(last));
  const auto st = run_pipeline(fl, pipeline_from_json(Json::parse(
      R"([{"type":"freeze_property","property":"Custom"},{"type":"compute_property","name":"Custom","expression":"Position.X + 2"}])")), last, 0, &t);
  ASSERT_EQ(st.steps[0].level, "ok") << st.steps[0].summary;
  const auto& fr = st.props.at("Custom frozen");
  const auto& now = st.props.at("Custom");
  bool moved = false;
  for (size_t i = 0; i < fl.atoms.size(); ++i) {
    EXPECT_NEAR(fr[i], f0.atoms[i].pos[0] + 2, 1e-9) << i;
    EXPECT_NEAR(now[i], fl.atoms[i].pos[0] + 2, 1e-9) << i;
    moved |= std::fabs(fr[i] - now[i]) > 1e-3;
  }
  EXPECT_TRUE(moved);
}

// Translucent atoms blend over what is behind them (the background here); a radius override draws bigger atoms; the
// transparency step writes the property the view draws
TEST(Render, TransparencyAndRadiusPerAtom) {
  System s;
  Atom a;
  a.element = 6, a.id = 1, a.pos = {0, 0, 0};
  s.atoms.push_back(a);
  a.id = 2, a.pos = {4, 0, 0};
  s.atoms.push_back(a);
  Renderer R;
  RenderOptions o;
  o.width = 120, o.height = 80, o.supersample = 1, o.outlines = false, o.depth_cue = false;
  Camera cam;
  auto px = [](const Image& im, int x, int y) { const size_t k = (size_t(y) * im.width + x) * 4; return std::array<int, 3>{im.rgba[k], im.rgba[k + 1], im.rgba[k + 2]}; };
  const Image opaque = R.render(s, cam, o);
  const auto bg = px(opaque, 1, 1);
  // find the two atoms' centres: columns where the middle row differs from the background
  std::vector<int> cols;
  for (int x = 0; x < opaque.width; ++x) if (px(opaque, x, opaque.height / 2) != bg) cols.push_back(x);
  ASSERT_FALSE(cols.empty());
  const int left = cols.front() + 2, row = opaque.height / 2;
  o.transparency = {0.5f, 0.f};
  const Image half = R.render(s, cam, o);
  o.transparency = {1.f, 0.f};
  const Image gone = R.render(s, cam, o);
  for (int c = 0; c < 3; ++c) {
    const int lo = std::min(bg[size_t(c)], px(opaque, left, row)[size_t(c)]), hi = std::max(bg[size_t(c)], px(opaque, left, row)[size_t(c)]);
    EXPECT_GE(px(half, left, row)[size_t(c)], lo - 1);
    EXPECT_LE(px(half, left, row)[size_t(c)], hi + 1);
    EXPECT_NEAR(px(gone, left, row)[size_t(c)], bg[size_t(c)], 1);
  }
  EXPECT_NE(px(half, left, row), px(opaque, left, row));
  o.transparency.clear();
  auto covered = [&](const Image& im) { int k = 0; for (int y = 0; y < im.height; ++y) for (int x = 0; x < im.width; ++x) k += px(im, x, y) != bg; return k; };
  o.radius = {1.5f, 0.f};
  const Image big = R.render(s, cam, o);
  EXPECT_GT(covered(big), covered(opaque));
  // the step: the selection translucent
  const auto st = run_pipeline(s, pipeline_from_json(Json::parse(R"([{"type":"transparency","value":0.6},{"type":"select_expression","expression":"Identifier == 1"}])")), 0, 0);
  EXPECT_NEAR(st.props.at("Transparency")[0], 0.6, 1e-12);
  EXPECT_EQ(st.props.at("Transparency")[1], 0.0);
}

// A dump's other columns survive: velocities into each frame, every numeric column by its name in the atoms' id order
// (the file lists them unsorted), and the magnitudes |f| and |v| from their components
TEST(Io, DumpKeepsVelocitiesAndColumns) {
  const auto path = (std::filesystem::temp_directory_path() / "caps_cols.lammpstrj").string();
  {
    std::ofstream o(path);
    for (int step : {0, 100}) {
      o << "ITEM: TIMESTEP\n" << step << "\nITEM: NUMBER OF ATOMS\n3\nITEM: BOX BOUNDS pp pp pp\n0 10\n0 10\n0 10\n";
      o << "ITEM: ATOMS id type x y z vx vy vz fx fy fz c_pe\n";
      o << "3 1 3 3 3 0.003 0 0 0 0 4 -3." << step << "\n";   // listed first, id 3
      o << "1 1 1 1 1 0.001 0 0 3 4 0 -1\n";
      o << "2 1 2 2 2 0 0.002 0 0 0 0 -2\n";
    }
  }
  const Trajectory t = read_lammps_dump(path);
  std::filesystem::remove(path);
  ASSERT_EQ(t.frames(), 2u);
  ASSERT_EQ(t.velocities.size(), 2u);
  EXPECT_NEAR(t.velocities[0][0][0], 0.001, 1e-12);   // atom id 1
  EXPECT_NEAR(t.velocities[0][2][0], 0.003, 1e-12);   // atom id 3
  EXPECT_NEAR(t.frame(1).velocities[1][1], 0.002, 1e-12);
  ASSERT_TRUE(t.columns.count("c_pe") && t.columns.count("|f|") && t.columns.count("|v|"));
  EXPECT_FLOAT_EQ(t.columns.at("c_pe")[0][0], -1.0f);
  EXPECT_FLOAT_EQ(t.columns.at("c_pe")[1][2], -3.1f);   // frame 2, id 3: "-3.100"
  EXPECT_FLOAT_EQ(t.columns.at("|f|")[0][0], 5.0f);     // (3, 4, 0)
  EXPECT_FLOAT_EQ(t.columns.at("|f|")[0][2], 4.0f);
}

// atoms stored wrapped one by one (chains split at the walls): the data file's image flags still join every bond —
// LAMMPS reads no "Inconsistent image flags"
TEST(LammpsData, ImageFlagsConsistentForWrappedAtoms) {
  System s = read_lammps_data(S + "/ps_melt.data");
  ASSERT_TRUE(s.cell.valid());
  size_t split = 0;
  for (auto& a : s.atoms) a.pos = s.cell.wrap(a.pos);
  s.unwrapped = false;
  const double half = 0.5 * std::min({norm(s.cell.a), norm(s.cell.b), norm(s.cell.c)});
  for (const auto& b : s.bonds) split += norm(s.atoms[b.j].pos - s.atoms[b.i].pos) > half;
  ASSERT_GT(split, 0u);   // the test makes bonds that cross the walls
  const std::string p = tmp("caps_imageflags.data");
  write_lammps_data(s, p);
  System r = read_lammps_data(p);   // positions unwrapped with the file's image flags, as LAMMPS does
  ASSERT_EQ(r.bonds.size(), s.bonds.size());
  for (const auto& b : r.bonds) EXPECT_LT(norm(r.atoms[b.j].pos - r.atoms[b.i].pos), 2.0);
  // every written position stays the same atom: its wrapped place unchanged
  for (size_t i = 0; i < s.atoms.size(); ++i) {
    const Vec3 d = s.cell.minimum_image(r.atoms[i].pos - s.atoms[i].pos);
    EXPECT_LT(norm(d), 1e-5);
  }
}
