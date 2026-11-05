#include <gtest/gtest.h>

#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>

#include "caps/io.hpp"
#include "caps/measured.hpp"
#include "caps/refine.hpp"

using namespace caps;

namespace {
std::string tmp(const std::string& name) { return (std::filesystem::temp_directory_path() / name).string(); }
void put(const std::string& path, const std::string& text) { std::ofstream(path) << text; }
const std::string kJcampHead =
    "##TITLE=test spectrum\n##JCAMP-DX=4.24\n##DATA TYPE=INFRARED SPECTRUM\n##XUNITS=1/CM\n##YUNITS=ABSORBANCE\n"
    "##XFACTOR=1.0\n##YFACTOR=0.5\n##FIRSTX=4000\n##LASTX=3994\n##NPOINTS=7\n##XYDATA=(X++(Y..Y))\n";
}  // namespace

// JCAMP-DX: the ordinates 10 12 15 15 15 14 −3 written plainly (AFFN), packed (PAC) and compressed (SQZ, DIF, DUP with the
// Y check values that repeat each DIF line's last ordinate) all decode to the same curve, × YFACTOR, from FIRSTX to LASTX
TEST(Measured, JcampForms) {
  const std::vector<double> want{10, 12, 15, 15, 15, 14, -3};
  for (const std::string body : {std::string("4000 10 12 15 15 15 14 -3\n"), std::string("4000+10+12+15+15+15+14-3\n"),
                                 std::string("4000A0KL%T\n3996A5jj7\n3994c\n")}) {
    const MeasuredCurve c = parse_measured(kJcampHead + body + "##END=\n", "x.jdx");
    ASSERT_EQ(c.y.size(), want.size()) << body;
    for (size_t k = 0; k < want.size(); ++k) {
      EXPECT_DOUBLE_EQ(c.y[k], 0.5 * want[k]) << body;
      EXPECT_DOUBLE_EQ(c.x[k], 4000.0 - double(k)) << body;
    }
    EXPECT_EQ(c.xunits, "1/CM");
    EXPECT_TRUE(c.notes.empty()) << body;
  }
  // a damaged check value is reported
  const MeasuredCurve bad = parse_measured(kJcampHead + "4000A0KL%T\n3996A4jj7\n3994c\n##END=\n", "x.jdx");
  ASSERT_FALSE(bad.notes.empty());
  EXPECT_NE(bad.notes[0].find("check"), std::string::npos);
  // (XY..XY) pairs
  const MeasuredCurve pts = parse_measured("##TITLE=p\n##JCAMP-DX=4.24\n##XFACTOR=2\n##YFACTOR=10\n##XYPOINTS=(XY..XY)\n1.0, 2.0; 3.0, 4.0\n##END=\n");
  ASSERT_EQ(pts.x.size(), 2u);
  EXPECT_DOUBLE_EQ(pts.x[1], 6.0);
  EXPECT_DOUBLE_EQ(pts.y[1], 40.0);
}

// XRDML: 2θ from start to end over the points (step = (end − start)/(n − 1)), the Omega axis ignored; listed positions
// and the 2.x counts element; read_pattern takes the file directly
TEST(Measured, Xrdml) {
  const std::string v15 =
      "<?xml version=\"1.0\"?><xrdMeasurements xmlns=\"http://www.xrdml.com/XRDMeasurement/1.5\"><xrdMeasurement><scan><dataPoints>"
      "<positions axis=\"Omega\" unit=\"deg\"><startPosition>5</startPosition><endPosition>5.18</endPosition></positions>"
      "<positions axis=\"2Theta\" unit=\"deg\"><startPosition>10.00</startPosition><endPosition>10.18</endPosition></positions>"
      "<commonCountingTime unit=\"seconds\">1.5</commonCountingTime><intensities unit=\"counts\">1 2 3 4 5 6 7 8 9 10</intensities>"
      "</dataPoints></scan></xrdMeasurement></xrdMeasurements>";
  const MeasuredCurve c = parse_measured(v15);
  ASSERT_EQ(c.x.size(), 10u);
  for (size_t k = 0; k < 10; ++k) {
    EXPECT_NEAR(c.x[k], 10.0 + 0.02 * double(k), 1e-12);
    EXPECT_DOUBLE_EQ(c.y[k], double(k + 1));
  }
  const std::string p = tmp("caps_measured.xrdml");
  put(p, v15);
  EXPECT_EQ(read_pattern(p).two_theta.size(), 10u);
  const std::string v21 =
      "<xrdMeasurements xmlns=\"http://www.xrdml.com/XRDMeasurement/2.1\"><xrdMeasurement><scan><dataPoints>"
      "<positions axis=\"2Theta\" unit=\"deg\"><listPositions>20 20.5 21.5</listPositions></positions>"
      "<counts unit=\"counts\">7 8 9</counts></dataPoints></scan><scan><dataPoints></dataPoints></scan></xrdMeasurement></xrdMeasurements>";
  const MeasuredCurve d = parse_measured(v21);
  ASSERT_EQ(d.x.size(), 3u);
  EXPECT_DOUBLE_EQ(d.x[2], 21.5);
  EXPECT_DOUBLE_EQ(d.y[2], 9);
}

// GSAS: a CONS bank in STD records (10 × (I2, F6.0); σ = √(Y/NCTR)) and in FXYE (centidegrees, intensity, σ)
TEST(Measured, Gsas) {
  std::string t = "test pattern\nBANK 1 12 2 CONST 1000.00 2.00 0 0 STD\n";
  char b[16];
  for (int r = 0; r < 2; ++r) {
    for (int f = 0; f < 10 && r * 10 + f < 12; ++f) {
      const int i = r * 10 + f;
      std::snprintf(b, sizeof b, "%2d%6d", i == 3 ? 4 : 1, 100 + i);
      t += b;
    }
    t += "\n";
  }
  const MeasuredCurve c = parse_measured(t, "x.gsa");
  ASSERT_EQ(c.y.size(), 12u);
  for (int i = 0; i < 12; ++i) {
    EXPECT_NEAR(c.x[size_t(i)], 10.0 + 0.02 * i, 1e-12);
    EXPECT_DOUBLE_EQ(c.y[size_t(i)], 100 + i);
    EXPECT_NEAR(c.sigma[size_t(i)], std::sqrt((100.0 + i) / (i == 3 ? 4 : 1)), 1e-12);
  }
  const MeasuredCurve e = parse_measured("t\nBANK 1 3 3 CONS 1000 2 0 0 FXYE\n1000.0 50.0 7.0\n1002.0 60.0 8.0\n1004.0 70.0 9.0\n");
  ASSERT_EQ(e.y.size(), 3u);
  EXPECT_DOUBLE_EQ(e.x[1], 10.02);
  EXPECT_DOUBLE_EQ(e.sigma[2], 9.0);
}

// UXD (single-column counts from _START by _STEPSIZE, the first range of two; two-column data) and Rigaku *TYPE text (a
// group's values from *START by *STEP)
TEST(Measured, UxdAndRigaku) {
  const MeasuredCurve u = parse_measured(";comment\n_FILEVERSION=1\n_SAMPLE='quartz'\n_DRIVE='COUPLED'\n_START=5.000\n_STEPSIZE=0.020\n_COUNTS\n 10 20 30 40\n 50 60\n"
                                         "_DRIVE='COUPLED'\n_START=50\n_STEPSIZE=0.02\n_COUNTS\n1 2\n");
  ASSERT_EQ(u.y.size(), 6u);
  EXPECT_NEAR(u.x[5], 5.10, 1e-12);
  EXPECT_DOUBLE_EQ(u.y[5], 60);
  EXPECT_EQ(u.title, "'quartz'");
  ASSERT_FALSE(u.notes.empty());
  const MeasuredCurve two = parse_measured("_FILEVERSION=1\n_DRIVE='COUPLED'\n_2THETACOUNTS\n 20.0 5\n 20.5 6\n");
  ASSERT_EQ(two.x.size(), 2u);
  EXPECT_DOUBLE_EQ(two.x[1], 20.5);
  const MeasuredCurve g = parse_measured("*TYPE = Raw\n*GROUP_COUNT = 1\n*BEGIN\n*START = 10.00\n*STOP = 10.06\n*STEP = 0.02\n*COUNT = 4\n 1048, 1162, 1108, 1163\n*END\n*EOF\n");
  ASSERT_EQ(g.y.size(), 4u);
  EXPECT_NEAR(g.x[3], 10.06, 1e-12);
  EXPECT_DOUBLE_EQ(g.y[1], 1162);
  EXPECT_TRUE(g.notes.empty());
}

// Gaussian cube: a grid CAPS writes reads back (values to the figures written, the cell and origin, the atoms); a file
// with an orbital list and Å voxels (negative counts) keeps the first value at each point
TEST(Measured, CubeRoundTripAndOrbitals) {
  GridField g;
  g.name = "field";
  g.cell.a = {6, 0, 0}, g.cell.b = {1, 8, 0}, g.cell.c = {0, 0.5, 10};
  g.cell.origin = {0.3, -0.2, 0.1};
  g.n[0] = 3, g.n[1] = 4, g.n[2] = 5;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 4; ++j)
      for (int k = 0; k < 5; ++k) g.values.push_back(std::sin(i + 2.0 * j) + 0.1 * k);
  System atoms;
  Atom a;
  a.element = 8, a.pos = {1.0, 2.0, 3.0};
  atoms.atoms = {a};
  const std::string p = tmp("caps_measured.cube");
  write_grid(g, atoms, p);
  System back;
  const GridField r = read_cube(p, &back);
  ASSERT_EQ(r.values.size(), g.values.size());
  for (size_t q = 0; q < g.values.size(); ++q) EXPECT_NEAR(r.values[q], g.values[q], 1e-5 * std::max(1.0, std::fabs(g.values[q])));
  EXPECT_NEAR(norm(r.cell.a - g.cell.a), 0, 1e-5);
  EXPECT_NEAR(norm(r.cell.b - g.cell.b), 0, 1e-5);
  EXPECT_NEAR(norm(r.cell.c - g.cell.c), 0, 1e-5);
  EXPECT_NEAR(norm(r.cell.origin - g.cell.origin), 0, 1e-5);
  ASSERT_EQ(back.atoms.size(), 1u);
  EXPECT_EQ(back.atoms[0].element, 8);
  EXPECT_NEAR(norm(back.atoms[0].pos - a.pos), 0, 1e-5);

  const std::string mo = "title\nsecond\n   -1    0.0 0.0 0.0\n   -2    1.0 0.0 0.0\n    2    0.0 1.0 0.0\n    1    0.0 0.0 1.0\n"
                         "    8    0.0 0.5 0.5 0.5\n    2   11   12\n 1.0 10.0 2.0 20.0\n 3.0 30.0 4.0 40.0\n";
  System at;
  const GridField m = parse_cube(mo, &at);
  ASSERT_EQ(m.values.size(), 4u);
  EXPECT_EQ(m.values, (std::vector<double>{1, 2, 3, 4}));
  EXPECT_NEAR(norm(m.cell.a - Vec3{2, 0, 0}), 0, 1e-12);   // Å: two voxels of 1 Å
  EXPECT_NEAR(at.atoms[0].pos[0], 0.5, 1e-12);
}

// An .arc archive: two .car frames, one trajectory
TEST(Measured, ArcFrames) {
  const std::string p = tmp("caps_measured.arc");
  put(p,
      "!BIOSYM archive 3\nPBC=OFF\nframe one\n!DATE Thu Nov 05 2026\n"
      "O1       0.000000000    0.000000000    0.000000000 XXXX 1      o*      O  -0.820\n"
      "H1       0.957200000    0.000000000    0.000000000 XXXX 1      h*      H   0.410\nend\nend\n"
      "frame two\n!DATE Thu Nov 05 2026\n"
      "O1       0.100000000    0.000000000    0.000000000 XXXX 1      o*      O  -0.820\n"
      "H1       1.000000000    0.000000000    0.000000000 XXXX 1      h*      H   0.410\nend\nend\n");
  const Trajectory t = open_file(p);
  ASSERT_EQ(t.frames(), 2u);
  EXPECT_EQ(t.topology.atoms.size(), 2u);
  EXPECT_DOUBLE_EQ(t.frame(1).atoms[0].pos[0], 0.1);
  EXPECT_DOUBLE_EQ(t.frame(1).atoms[1].pos[0], 1.0);
  EXPECT_DOUBLE_EQ(t.topology.atoms[0].charge, -0.82);
}
