#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <map>

#include "caps/crystal.hpp"
#include "caps/dynamics.hpp"
#include "caps/mechanics.hpp"
#include "caps/properties.hpp"
#include "caps/elements.hpp"
#include "caps/io.hpp"
#include "caps/polymer.hpp"
#include "caps/relax.hpp"
#include "caps/uff.hpp"

using namespace caps;

namespace {

const std::string kCrystals = std::string(CAPS_SOURCE_DIR) + "/data/crystals/";

std::vector<int> coordination(const System& s) {
  std::vector<int> cn(s.atoms.size(), 0);
  for (const auto& b : s.bonds) ++cn[b.i], ++cn[b.j];
  return cn;
}

}  // namespace

TEST(Crystal, CifSymmetryExpansion) {
  const System q = read_cif(kCrystals + "alpha-quartz.cif");
  ASSERT_EQ(q.atoms.size(), 9u);   // Si3O6 from two sites and six operations
  EXPECT_NEAR(q.density(), 2.649, 0.005);
  ASSERT_EQ(q.bonds.size(), 12u);
  for (const auto& b : q.bonds) {
    const double r = norm(q.cell.minimum_image(q.atoms[b.j].pos - q.atoms[b.i].pos));
    EXPECT_NEAR(r, 1.61, 0.01);
  }
  const auto cn = coordination(q);
  for (size_t i = 0; i < q.atoms.size(); ++i) EXPECT_EQ(cn[i], q.atoms[i].element == 14 ? 4 : 2);
  const System nacl = read_cif(kCrystals + "rock-salt.cif");   // face centring from four operations
  EXPECT_EQ(nacl.atoms.size(), 8u);
  const System g = open_file(kCrystals + "graphite.cif").frame(0);   // CIF through the general reader
  EXPECT_EQ(g.atoms.size(), 4u);
  EXPECT_THROW(parse_cif("data_x\n_cell_length_a 3\n"), CrystalError);
  const System p1 = parse_cif("data_x\n_cell_length_a 4\n_cell_length_b 4\n_cell_length_c 4\nloop_\n_atom_site_label\n_atom_site_fract_x\n_atom_site_fract_y\n"
                              "_atom_site_fract_z\nAr1 0 0 0\n");
  EXPECT_EQ(p1.atoms.size(), 1u);
  EXPECT_EQ(p1.atoms[0].element, 18);
}

TEST(Crystal, TerminationsCountBrokenBonds) {
  const System q = read_cif(kCrystals + "alpha-quartz.cif");
  double d = 0;
  const auto tq = slab_terminations(q, 0, 0, 1, &d);
  EXPECT_NEAR(d, 5.4052, 1e-4);
  EXPECT_NEAR(tq.front().bonds_per_nm2, 9.6, 0.1);   // the hydroxylated quartz (001) silanol density
  const auto tr = slab_terminations(read_cif(kCrystals + "rutile.cif"), 1, 1, 0);
  EXPECT_NEAR(tr.front().bonds_per_nm2, 10.4, 0.1);   // rutile (110): the stoichiometric bridging-O surface
  const auto td = slab_terminations(read_cif(kCrystals + "diamond.cif"), 1, 1, 1);
  EXPECT_NEAR(td.front().bonds_per_nm2, 18.2, 0.2);   // one dangling bond per (111) surface atom
  const auto tg = slab_terminations(read_cif(kCrystals + "graphite.cif"), 0, 0, 1);
  EXPECT_EQ(tg.front().bonds_cut, 0);                 // graphite cleaves between sheets
  EXPECT_THROW(slab_terminations(q, 0, 0, 0), CrystalError);
}

TEST(Crystal, CleaveOrthogonalSlabAndPassivate) {
  const System q = read_cif(kCrystals + "alpha-quartz.cif");
  SlabOptions o;
  o.layers = 2;
  o.na = 2;
  SlabReport r;
  const System s = cleave(q, o, &r);
  EXPECT_NEAR(s.cell.b[0], 0.0, 1e-9);   // rectangular
  EXPECT_NEAR(r.gamma, 90.0, 1e-6);
  EXPECT_EQ(s.atoms.size(), 9u * 2 * 2 * 2);   // 9 per mesh × 2 layers × 2 meshes (rectangular) × 2 (na)
  std::map<int, int> count;
  for (const auto& a : s.atoms) count[a.element]++;
  EXPECT_EQ(count[8], 2 * count[14]);   // whole layers: stoichiometric
  o.passivate = true;
  const System p = cleave(q, o, &r);
  EXPECT_GT(r.added_oh + r.added_h, 0);
  const auto cn = coordination(p);
  for (size_t i = 0; i < p.atoms.size(); ++i) {
    if (p.atoms[i].element == 14) EXPECT_EQ(cn[i], 4) << "Si " << i;
    if (p.atoms[i].element == 8) EXPECT_EQ(cn[i], 2) << "O " << i;
    if (p.atoms[i].element == 1) EXPECT_EQ(cn[i], 1) << "H " << i;
  }
  // graphite sheets 3.355 Å apart
  SlabOptions og;
  og.layers = 2;
  const System g = cleave(read_cif(kCrystals + "graphite.cif"), og);
  std::vector<double> zs;
  for (const auto& a : g.atoms) zs.push_back(a.pos[2]);
  std::sort(zs.begin(), zs.end());
  zs.erase(std::unique(zs.begin(), zs.end(), [](double x, double y) { return std::fabs(x - y) < 0.1; }), zs.end());
  ASSERT_EQ(zs.size(), 4u);
  EXPECT_NEAR(zs[1] - zs[0], 3.3555, 1e-3);
}

TEST(Crystal, RubberFilmOnSilicaAndRelaxWithTheSurfaceHeld) {
  const System q = read_cif(kCrystals + "alpha-quartz.cif");
  SlabOptions so;
  so.layers = 1;
  so.na = 2, so.nb = 1;
  so.passivate = true;
  const System slab = cleave(q, so);
  ChainSpec spec;
  spec.units = {{"cis-1,4-isoprene", "*C/C=C(/C)C*"}};
  spec.dp = 6;
  InterfaceOptions io;
  io.film = 14;
  io.density = 0.6;
  GrowReport rep;
  const System s = build_interface(slab, spec, io, &rep);
  ASSERT_GT(s.atoms.size(), slab.atoms.size());
  size_t nsub = 0;
  double top = -1e300, zlo = 1e300;
  for (const auto& a : s.atoms) {
    if (a.mol == 1) ++nsub, top = std::max(top, a.pos[2]);
  }
  EXPECT_EQ(nsub, slab.atoms.size());
  for (const auto& a : s.atoms)
    if (a.mol != 1 && a.element != 1) zlo = std::min(zlo, a.pos[2]);
  EXPECT_GT(zlo, top);   // the film sits above the surface
  // UFF on silica + rubber, the surface held in place
  RelaxOptions ro;
  ro.field = std::make_shared<ForceField>(assign_uff(s));
  ro.fixed.assign(s.atoms.size(), 0);
  for (size_t i = 0; i < s.atoms.size(); ++i) ro.fixed[i] = s.atoms[i].mol == 1;
  ro.ftol = 2.0;
  System r = s;
  relax(r, ro);
  double moved = 0;
  for (size_t i = 0; i < s.atoms.size(); ++i)
    if (ro.fixed[i]) moved = std::max(moved, norm(r.atoms[i].pos - s.atoms[i].pos));
  EXPECT_LT(moved, 1e-12);
}

TEST(Crystal, DynamicsHoldsTheSurface) {
  SlabOptions so;
  so.layers = 1;
  so.na = 2;
  so.passivate = true;
  const System slab = cleave(read_cif(kCrystals + "alpha-quartz.cif"), so);
  ChainSpec spec;
  spec.units = {{"butadiene", "*C/C=C\\C*"}};
  spec.dp = 5;
  InterfaceOptions io;
  io.film = 12;
  io.density = 0.5;
  System s = build_interface(slab, spec, io);
  DynamicsOptions d;
  d.field = std::make_shared<ForceField>(assign_uff(s));
  d.fixed.assign(s.atoms.size(), 0);
  size_t nh = 0;
  for (size_t i = 0; i < s.atoms.size(); ++i) nh += (d.fixed[i] = s.atoms[i].mol == 1);
  d.steps = 50;
  d.dt = 0.5;
  d.thermostat = Thermostat::Langevin;
  d.frame_every = 0;
  const System s0 = s;
  DynamicsReport rep;
  run_dynamics(s, d, &rep);
  double moved_held = 0, moved_free = 0;
  for (size_t i = 0; i < s.atoms.size(); ++i) {
    const double dd = norm(s.atoms[i].pos - s0.atoms[i].pos);
    if (d.fixed[i]) moved_held = std::max(moved_held, dd);
    else moved_free = std::max(moved_free, dd);
  }
  EXPECT_EQ(moved_held, 0.0);
  EXPECT_GT(moved_free, 0.01);
  EXPECT_GT(nh, 0u);
  EXPECT_NEAR(rep.thermo.front().temperature, 300.0, 1.0);   // drawn over the free atoms only
}

TEST(Crystal, InterfaceProfileAndAdhesion) {
  SlabOptions so;
  so.layers = 2;
  so.na = 3, so.nb = 2;
  so.passivate = true;
  const System slab = cleave(read_cif(kCrystals + "alpha-quartz.cif"), so);
  ChainSpec spec;
  spec.units = {{"cis-1,4-isoprene", "*C/C=C(/C)C*"}};
  spec.dp = 8;
  InterfaceOptions io;
  io.film = 18;
  io.density = 0.7;
  io.vacuum = 15;   // one interface
  const System s = build_interface(slab, spec, io);
  Trajectory t;
  t.topology = s;
  std::vector<Vec3> p;
  for (const auto& a : s.atoms) p.push_back(a.pos);
  t.positions.push_back(p);
  t.cells.push_back(s.cell);
  t.timesteps.push_back(0);
  const ForceField ff = assign_uff(s);
  AnalyzeOptions o;
  o.ff = &ff;
  const auto props = analyze(t, {"zprofile", "adhesion"}, o);
  ASSERT_EQ(props.size(), 2u);
  EXPECT_EQ(props[0].series.size(), 3u);
  EXPECT_NEAR(props[0].value, 0.7, 0.25);   // the film's own density
  EXPECT_GT(props[1].value, 5.0);           // the film sticks: positive work of adhesion, mJ/m²
  EXPECT_LT(props[1].value, 500.0);
  EXPECT_EQ(props[1].extra.at("interfaces"), 1.0);
}

TEST(Crystal, PullOutFromTheSurface) {
  SlabOptions so;
  so.layers = 1;
  so.na = 3, so.nb = 2;
  so.passivate = true;
  const System slab = cleave(read_cif(kCrystals + "alpha-quartz.cif"), so);
  ChainSpec spec;
  spec.units = {{"cis-1,4-isoprene", "*C/C=C(/C)C*"}};
  spec.dp = 6;
  InterfaceOptions io;
  io.film = 12;
  io.density = 0.6;
  System s = build_interface(slab, spec, io);
  const System s0 = s;
  PullOptions po;
  po.distance = 2.0;
  po.rate = 10.0;
  po.equilibrate_ps = 0.1;
  po.dt = 0.5;
  po.sample_every = 10;
  const PullResult r = run_pull(s, po);
  ASSERT_GT(r.curve.size(), 5u);
  EXPECT_GT(r.peak_force, 0.0);
  EXPECT_GT(r.curve.back().displacement, 0.0);   // the film moved along +x
  EXPECT_EQ(r.interfaces, 2);
  for (size_t i = 0; i < s.atoms.size(); ++i)
    if (s.atoms[i].mol == 1) EXPECT_EQ(norm(s.atoms[i].pos - s0.atoms[i].pos), 0.0);
}
