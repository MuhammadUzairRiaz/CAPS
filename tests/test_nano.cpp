#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <numeric>

#include "caps/crystal.hpp"
#include "caps/mechanics.hpp"
#include "caps/edit.hpp"
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

// Multi-walled tubes: (5,5)@(10,10)@(15,15) has 20 + 40 + 60 atoms per period and walls 3.39 Å apart; zigzag walls
// step by (9,0); chiral walls are refused (their periods along the axis differ)
TEST(Nano, MultiWalledTubes) {
  NanotubeOptions o;
  o.n = 5, o.m = 5, o.walls = 3, o.length = 24.6;
  NanoReport r;
  const System s = nanotube(o, &r);
  EXPECT_EQ(r.atoms_per_period, 120);
  EXPECT_EQ(int(s.atoms.size()), 120 * 10);
  EXPECT_NEAR(r.diameter, 2 * 15 * std::sqrt(3.0) * 1.42 * std::sqrt(3.0) / (2 * M_PI), 1e-6);
  EXPECT_NE(r.notes.front().find("(5,5)@(10,10)@(15,15) 3-walled"), std::string::npos) << r.notes.front();
  EXPECT_NE(r.notes.front().find("3.39 Å apart"), std::string::npos) << r.notes.front();
  o.n = 9, o.m = 0, o.walls = 2;
  EXPECT_EQ((nanotube(o, &r), r.atoms_per_period), 36 + 72);
  o.n = 8, o.m = 4;
  EXPECT_THROW(nanotube(o), std::invalid_argument);
}

// Truncated octahedron and icosahedron: at the same circumscribed radius a copper particle holds the atoms of their
// volumes relative to the sphere (32/√5³ ÷ 4π/3 = 0.683 and 0.6055) once the particle is large enough (30 Å) for the
// lattice's layering at the faces to matter little (at 18 Å every shape, the cube and octahedron too, is 5–10 % off)
TEST(Nano, TruncatedOctahedronAndIcosahedron) {
  const System bulk = read_cif(kCrystals + "copper.cif");
  auto count = [&](ParticleShape sh) {
    ParticleOptions o;
    o.shape = sh;
    o.radius = 30;
    const System p = nanoparticle(bulk, o);
    Vec3 c{0, 0, 0};
    for (const auto& a : p.atoms) c = c + a.pos;
    c = c * (1.0 / double(p.atoms.size()));
    double rmax = 0;
    for (const auto& a : p.atoms) rmax = std::max(rmax, norm(a.pos - c));
    EXPECT_LT(rmax, 30.0 + 1.5) << to_string(sh);
    return double(p.atoms.size());
  };
  const double sphere = count(ParticleShape::Sphere);
  EXPECT_NEAR(count(ParticleShape::TruncatedOctahedron) / sphere, 32 / std::pow(5.0, 1.5) / (4 * M_PI / 3), 0.06);
  EXPECT_NEAR(count(ParticleShape::Icosahedron) / sphere, 0.6055, 0.06);
  EXPECT_EQ(particle_shape_from_string("icosahedron"), ParticleShape::Icosahedron);
}

// TESPT grafted on a passivated silica particle: each graft trades a silanol H for the silane bonded through its Si
// (64 atoms: net +63), the sites at least the spacing apart, no atom closer than 0.9 Å to another
TEST(Nano, SilaneGraftOnSilica) {
  ParticleOptions po;
  po.radius = 12;
  po.passivate = true;
  System p = nanoparticle(read_cif(kCrystals + "alpha-quartz.cif"), po);
  const size_t n0 = p.atoms.size();
  {
    double d0 = 1e9;
    size_t bi = 0, bj = 0;
    for (size_t i = 0; i < n0; ++i)
      for (size_t j = i + 1; j < n0; ++j) if (norm(p.atoms[i].pos - p.atoms[j].pos) < d0) d0 = norm(p.atoms[i].pos - p.atoms[j].pos), bi = i, bj = j;
    const auto nb0 = p.neighbours();
    EXPECT_GT(d0, 0.9) << "passivation left atoms " << bi << " (Z " << p.atoms[bi].element << ", " << nb0[bi].size() << " bonds) and " << bj << " " << d0 << " Å apart";
  }
  GraftOptions g;
  g.count = 4;
  g.min_spacing = 6;
  const GraftReport r = graft_silanes(p, g);
  ASSERT_EQ(r.grafted, 4u);
  EXPECT_GT(r.silanols, 4u);
  EXPECT_EQ(p.atoms.size(), n0 + 4 * 63);
  size_t sulfur = 0;
  for (const auto& a : p.atoms) sulfur += a.element == 16;
  EXPECT_EQ(sulfur, 4u * 4u);
  double dmin = 1e9;
  for (size_t i = 0; i < p.atoms.size(); ++i)
    for (size_t j = i + 1; j < p.atoms.size(); ++j) dmin = std::min(dmin, norm(p.atoms[i].pos - p.atoms[j].pos));
  EXPECT_GT(dmin, 0.9);
  EXPECT_THROW(graft_silanes(p, GraftOptions{"CCC", "bad"}), EditError);
}

// A gold sphere capped with dodecanethiolates: the 429-atom particle of the Nanoparticle board, each S 2.45 Å from its
// nearest gold, S atoms at least √3 · d(Au–Au) apart, tails all-trans and clear of the gold and of each other
TEST(Nano, GoldParticleCappedWithThiolates) {
  const System bulk = read_cif(kCrystals + "gold.cif");
  ParticleOptions o;
  o.radius = 12;
  o.thiolate = "C12";
  const System p = nanoparticle(bulk, o);
  std::vector<Vec3> au, sulfur;
  for (const auto& a : p.atoms) {
    if (a.element == 79) au.push_back(a.pos);
    if (a.element == 16) sulfur.push_back(a.pos);
  }
  EXPECT_EQ(au.size(), 429u);
  ASSERT_GT(sulfur.size(), 30u);
  EXPECT_EQ(p.atoms.size(), au.size() + sulfur.size() * (1 + 12 + 25));   // S, C12H25
  const double dnn = 4.0782 / std::sqrt(2.0);
  for (size_t i = 0; i < sulfur.size(); ++i) {
    double dmin = 1e9;
    for (const auto& g : au) dmin = std::min(dmin, norm(sulfur[i] - g));
    EXPECT_NEAR(dmin, 2.45, 0.02);
    for (size_t j = i + 1; j < sulfur.size(); ++j) EXPECT_GE(norm(sulfur[i] - sulfur[j]), std::sqrt(3.0) * dnn - 1e-6);
  }
  // no ligand atom closer than 3 Å to the gold (besides its S), none of two ligands closer than 2.5 Å
  for (const auto& a : p.atoms) {
    if (a.element == 79 || a.element == 16) continue;
    for (const auto& g : au) EXPECT_GT(norm(a.pos - g), 3.0);
  }
  for (size_t i = 0; i < p.atoms.size(); ++i)
    for (size_t j = i + 1; j < p.atoms.size(); ++j)
      if (p.atoms[i].element != 79 && p.atoms[j].element != 79 && p.atoms[i].mol != p.atoms[j].mol)
        EXPECT_GT(norm(p.atoms[i].pos - p.atoms[j].pos), 2.5);
  // every ligand inside the cell with vacuum around it
  for (const auto& a : p.atoms)
    for (int k = 0; k < 3; ++k) { EXPECT_GT(a.pos[k], 5.0); EXPECT_LT(a.pos[k], p.cell.a[0] - 5.0); }
}

// Boron nitride: a periodic h-BN sheet has as many B as N, every bond B–N at 1.446 Å and three per atom; a (10,10) BN tube
// has the rolled diameter a√(n² + nm + m²)/π with a = 2.504 Å and only B–N bonds; a flake's edges get B–H 1.19, N–H 1.01 Å
TEST(Nano, BoronNitrideSheetAndTube) {
  SheetOptions so;
  so.material = "h-BN";
  so.layers = 2;
  const System sh = graphene_sheet(so);
  int nb = 0, nn = 0;
  for (const auto& a : sh.atoms) nb += a.element == 5, nn += a.element == 7;
  EXPECT_EQ(nb, nn);
  EXPECT_GT(nb, 0);
  std::vector<int> deg(sh.atoms.size(), 0);
  for (const auto& b : sh.bonds) {
    const int zi = sh.atoms[b.i].element, zj = sh.atoms[b.j].element;
    EXPECT_EQ(zi + zj, 12) << "bond " << zi << "–" << zj;
    EXPECT_NEAR(norm(sh.cell.minimum_image(sh.atoms[b.j].pos - sh.atoms[b.i].pos)), 2.504 / std::sqrt(3.0), 1e-6);
    ++deg[b.i], ++deg[b.j];
  }
  for (int d : deg) EXPECT_EQ(d, 3);
  // AA′: the second layer's atom above each B is an N 3.328 Å up
  NanotubeOptions to;
  to.material = "h-BN";
  to.n = 10, to.m = 10;
  NanoReport r;
  const System t = nanotube(to, &r);
  EXPECT_NEAR(r.diameter, 2.504 * std::sqrt(300.0) / 3.14159265358979323846, 1e-6);
  for (const auto& b : t.bonds) EXPECT_EQ(t.atoms[b.i].element + t.atoms[b.j].element, 12);
  SheetOptions fo = so;
  fo.layers = 1;
  fo.periodic = false;
  const System fl = graphene_sheet(fo);
  int caps = 0;
  for (const auto& b : fl.bonds) {
    const int zi = fl.atoms[b.i].element, zj = fl.atoms[b.j].element;
    if (zi != 1 && zj != 1) continue;
    const int heavy = zi == 1 ? zj : zi;
    EXPECT_NEAR(norm(fl.atoms[b.j].pos - fl.atoms[b.i].pos), heavy == 5 ? 1.19 : 1.01, 1e-6);
    ++caps;
  }
  EXPECT_GT(caps, 0);
  EXPECT_THROW(graphene_sheet(SheetOptions{"MoS2"}), std::invalid_argument);
}

TEST(Nano, MoreParticleShapesFollowTheirVolumes) {
  // gold (fcc, a 4.0782 Å): 4 atoms per a³; each shape's atom count ≈ density × its volume (the surface's granularity
  // allows ~15 % at these sizes)
  const System au = read_cif(kCrystals + "gold.cif");
  const double a = 4.0782, rho = 4 / (a * a * a), R = 14, H = 24, t = 0.5;
  struct Case { const char* name; double volume; };
  // the tetrahedron's faces are fcc {111} planes (spacing a/√3): with its inradius R/3 midway between two planes the
  // atoms fill exactly the geometric volume (on a plane, whole layers tip the count by up to half a spacing per face)
  const double Rt = 3 * 4.5 * a / std::sqrt(3.0);
  const double tet_edge = Rt * std::sqrt(8.0 / 3.0);
  const Case cases[] = {
      {"rod", M_PI * R * R * H},
      {"cone", M_PI * R * R * H / 3},
      {"frustum", M_PI * H / 3 * (R * R + R * R * t + R * R * t * t)},
      {"tetrahedron", tet_edge * tet_edge * tet_edge / (6 * std::sqrt(2.0))},
      {"pyramid", 2 * R * R * H / 3},
      {"hemisphere", 2.0 / 3 * M_PI * R * R * R},
  };
  for (const auto& c : cases) {
    ParticleOptions o;
    o.shape = particle_shape_from_string(c.name);
    o.radius = std::string(c.name) == "tetrahedron" ? Rt : R;
    o.height = std::string(c.name) == "tetrahedron" || std::string(c.name) == "hemisphere" ? 0 : H;
    o.top_ratio = t;
    const System p = nanoparticle(au, o);
    const double expect = rho * c.volume;
    EXPECT_NEAR(double(p.atoms.size()), expect, 0.15 * expect) << c.name;
    EXPECT_STREQ(to_string(o.shape), c.name);
  }
}

// Each wall its own (n, m): (5,5)@(10,10)@(26,0), 3.39 and 3.40 Å apart — the zigzag outer wall's period (√3 a) differs
// from the armchair walls' (a), so the periodic tube is made long enough for both within 2 % — and ropes: 7 tubes about one, and the periodic triangular
// lattice, axes 2R + 3.4 Å apart, each tube its own molecule.
TEST(Nano, WallChiralitiesAndRopes) {
  NanotubeOptions o;
  o.wall_chiralities = {{5, 5}, {10, 10}, {26, 0}};
  o.length = 60;
  NanoReport r;
  const System t = nanotube(o, &r);
  const double ar = 1.42 * std::sqrt(3.0);
  const auto atoms_per = [&](int n, int m) { const int d = std::gcd(2 * m + n, 2 * n + m); return 4 * (n * n + n * m + m * m) / d; };
  EXPECT_EQ(r.atoms_per_period, atoms_per(5, 5) + atoms_per(10, 10) + atoms_per(26, 0));
  // every atom on one of the three radii
  for (const auto& a : t.atoms) {
    const double rr = std::hypot(a.pos[0] - t.cell.a[0] / 2, a.pos[1] - t.cell.b[1] / 2);
    const double r5 = ar * std::sqrt(75.0) / (2 * M_PI), r10 = ar * std::sqrt(300.0) / (2 * M_PI), r26 = ar * 26 / (2 * M_PI);
    EXPECT_LT(std::min({std::fabs(rr - r5), std::fabs(rr - r10), std::fabs(rr - r26)}), 1e-6);
  }
  // too close walls are refused, and walls whose periods cannot share a short length
  o.wall_chiralities = {{5, 5}, {6, 6}};
  EXPECT_THROW(nanotube(o), std::invalid_argument);
  o.wall_chiralities = {{5, 5}, {10, 10}, {26, 0}};
  o.length = 12;
  EXPECT_THROW(nanotube(o), std::invalid_argument);
  // a rope of seven (5,5) tubes, each its own molecule, neighbours 2R + 3.4 Å apart
  NanotubeOptions b;
  b.n = 5, b.m = 5, b.length = 10, b.bundle = 7;
  NanoReport br;
  const System rope = nanotube(b, &br);
  int64_t top = 0;
  for (const auto& a : rope.atoms) top = std::max(top, a.mol);
  EXPECT_EQ(top, 7);
  EXPECT_EQ(rope.atoms.size() % 7, 0u);
  // the periodic lattice: two tubes per cell, the cell a = D, b = D√3
  b.bundle_lattice = true;
  const System lat = nanotube(b);
  const double D = 2 * ar * std::sqrt(75.0) / (2 * M_PI) + 3.4;
  EXPECT_NEAR(lat.cell.a[0], D, 1e-9);
  EXPECT_NEAR(lat.cell.b[1], D * std::sqrt(3.0), 1e-9);
  int64_t ltop = 0;
  for (const auto& a : lat.atoms) ltop = std::max(ltop, a.mol);
  EXPECT_EQ(ltop, 2);
}

// Point defects in an aluminium nanoparticle: 10 % vacancies (that many atoms fewer), then copper on 5 Al sites at least
// 5 Å apart (the count kept, the spacing held).
TEST(Nano, VacanciesAndDoping) {
  const System al = read_cif(std::string(CAPS_SOURCE_DIR) + "/data/crystals/aluminium.cif");
  ParticleOptions po;
  po.radius = 10;
  System p = nanoparticle(al, po);
  const size_t n0 = p.atoms.size();
  DefectOptions v;
  v.from = 13;
  v.fraction = 0.1;
  std::vector<std::string> notes;
  const int nv = point_defects(p, v, &notes);
  EXPECT_EQ(nv, int(std::lround(0.1 * double(n0))));
  EXPECT_EQ(p.atoms.size(), n0 - size_t(nv));
  DefectOptions d;
  d.from = 13, d.to = 29, d.count = 5, d.min_spacing = 5.0, d.seed = 3;
  EXPECT_EQ(point_defects(p, d), 5);
  std::vector<Vec3> cu;
  for (const auto& a : p.atoms) if (a.element == 29) cu.push_back(a.pos);
  ASSERT_EQ(cu.size(), 5u);
  for (size_t i = 0; i < cu.size(); ++i)
    for (size_t j = i + 1; j < cu.size(); ++j) EXPECT_GE(norm(p.cell.valid() ? p.cell.minimum_image(cu[i] - cu[j]) : cu[i] - cu[j]), 5.0 - 1e-9);
  d.to = 13;
  d.from = 13;
  EXPECT_THROW(point_defects(p, d), std::invalid_argument);
}
