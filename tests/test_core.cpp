#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <set>

#include "caps/analysis.hpp"
#include "caps/checks.hpp"
#include "caps/elements.hpp"
#include "caps/io.hpp"
#include "caps/render.hpp"
#include "caps/pipeline.hpp"

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
  EXPECT_EQ(st.steps[1].summary, std::to_string(nh) + " selected");
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
