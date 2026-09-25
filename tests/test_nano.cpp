#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>

#include "caps/crystal.hpp"
#include "caps/mechanics.hpp"
#include "caps/nano.hpp"

using namespace caps;

namespace {
const std::string kCrystals = std::string(CAPS_SOURCE_DIR) + "/data/crystals/";
std::vector<int> degrees(const System& s) {
  std::vector<int> d(s.atoms.size(), 0);
  for (const auto& b : s.bonds) ++d[b.i], ++d[b.j];
  return d;
}
}  // namespace

TEST(Nano, NanotubeGeometryAndBonds) {
  const auto g = nanotube_geometry(10, 10);
  EXPECT_NEAR(g[0], 13.56, 0.01);   // (10,10) diameter
  EXPECT_NEAR(g[1], 30.0, 1e-9);
  EXPECT_NEAR(g[2], 2.46, 0.01);
  EXPECT_NEAR(nanotube_geometry(10, 0)[0], 7.83, 0.01);
  NanotubeOptions o;
  o.n = 8, o.m = 4;
  NanoReport r;
  const System t = nanotube(o, &r);
  EXPECT_EQ(r.atoms_per_period, 112);
  for (int d : degrees(t)) EXPECT_EQ(d, 3);   // periodic: every carbon sp2
  for (const auto& b : t.bonds) EXPECT_NEAR(norm(t.cell.minimum_image(t.atoms[b.j].pos - t.atoms[b.i].pos)), 1.42, 0.03);
  o.periodic = false;
  const System f = nanotube(o, &r);
  EXPECT_GT(r.capped, 0);
  for (size_t i = 0; i < f.atoms.size(); ++i)
    if (f.atoms[i].element == 6) EXPECT_EQ(degrees(f)[i], 3);
}

TEST(Nano, GrapheneSheetAndFlake) {
  SheetOptions o;
  o.layers = 2;
  const System s = graphene_sheet(o);
  for (int d : degrees(s)) EXPECT_EQ(d, 3);
  EXPECT_EQ(s.bonds.size() * 2, s.atoms.size() * 3);
  o.layers = 1;
  o.periodic = false;
  NanoReport r;
  const System f = graphene_sheet(o, &r);
  const auto d = degrees(f);
  for (size_t i = 0; i < f.atoms.size(); ++i) EXPECT_EQ(d[i], f.atoms[i].element == 6 ? 3 : 1);
}

TEST(Nano, ParticlesCutFromCrystals) {
  ParticleOptions o;
  o.shape = ParticleShape::Cuboctahedron;
  o.radius = 10;
  const System cu = nanoparticle(read_cif(kCrystals + "copper.cif"), o);
  EXPECT_EQ(cu.atoms.size(), 147u);   // the magic-number cuboctahedron
  o.shape = ParticleShape::Sphere;
  o.radius = 9;
  o.passivate = true;
  NanoReport r;
  const System si = nanoparticle(read_cif(kCrystals + "alpha-quartz.cif"), o, &r);
  EXPECT_GT(r.added_oh + r.added_h, 0);
  const auto d = degrees(si);
  for (size_t i = 0; i < si.atoms.size(); ++i) {
    if (si.atoms[i].element == 14) EXPECT_EQ(d[i], 4);
    if (si.atoms[i].element == 1) EXPECT_EQ(d[i], 1);
  }
}

TEST(Nano, FillerInARubberMatrix) {
  ParticleOptions po;
  po.radius = 7;
  po.passivate = true;
  const System f = nanoparticle(read_cif(kCrystals + "alpha-quartz.cif"), po);
  ChainSpec spec;
  spec.units = {{"cis-1,4-isoprene", "*C/C=C(/C)C*"}};
  spec.dp = 10;
  FillerMatrixOptions fo;
  fo.chains = 6;
  fo.density = 0.6;
  FillerReport fr;
  const System s = embed_filler(f, spec, fo, &fr);
  EXPECT_EQ(s.atoms.size(), f.atoms.size() + 6u * (10 * 13 + 2));
  EXPECT_GT(fr.filler_mass_fraction, 0.1);
  EXPECT_LT(fr.filler_volume_fraction, 0.5);
  size_t held = 0;
  for (const auto& a : s.atoms) held += a.mol == 1;
  EXPECT_EQ(held, f.atoms.size());
  // a periodic nanotube keeps its axis: the cell is one tube period long
  NanotubeOptions to;
  to.n = to.m = 5;
  to.length = 12;
  const System tube = nanotube(to);
  fo.keep_axis = {false, false, true};
  fo.chains = 4;
  spec.dp = 6;
  const System c = embed_filler(tube, spec, fo, &fr);
  EXPECT_NEAR(norm(c.cell.c), norm(tube.cell.c), 1e-9);
}

TEST(Nano, SilicaFibreInRubberAndPullOut) {
  ParticleOptions po;
  po.shape = ParticleShape::Fibre;
  po.radius = 6;
  po.length = 16;
  po.passivate = true;
  const System f = nanoparticle(read_cif(kCrystals + "alpha-quartz.cif"), po);
  EXPECT_NEAR(norm(f.cell.c), 3 * 5.4052, 1e-3);   // whole cells along z, periodic
  const auto d = degrees(f);
  for (size_t i = 0; i < f.atoms.size(); ++i) {
    if (f.atoms[i].element == 14) EXPECT_EQ(d[i], 4);
    if (f.atoms[i].element == 8) EXPECT_EQ(d[i], 2);
  }
  ChainSpec spec;
  spec.units = {{"cis-1,4-isoprene", "*C/C=C(/C)C*"}};
  spec.dp = 6;
  FillerMatrixOptions fo;
  fo.chains = 5;
  fo.density = 0.6;
  fo.keep_axis = {false, false, true};
  System s = embed_filler(f, spec, fo);
  EXPECT_NEAR(norm(s.cell.c), norm(f.cell.c), 1e-9);
  PullOptions pp;
  pp.axis = 2;
  pp.distance = 1.5;
  pp.rate = 10;
  pp.equilibrate_ps = 0.1;
  pp.dt = 0.5;
  const PullResult r = run_pull(s, pp);
  const double side = 2 * 3.14159265358979 * norm(f.cell.c);   // 2πRL per Å of radius; R reaches the surface hydroxyls
  EXPECT_GT(r.area, 5 * side);
  EXPECT_LT(r.area, 9.5 * side);
  EXPECT_EQ(r.interfaces, 1);
  EXPECT_GT(r.curve.size(), 3u);
}
