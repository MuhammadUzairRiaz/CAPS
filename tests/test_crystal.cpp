#include <gtest/gtest.h>
#include "caps/lattice.hpp"
#include "caps/spacegroup.hpp"

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
  EXPECT_EQ(props[0].series.size(), 3u + 3u);   // film, surface, all atoms + the surface's Si, O and H
  EXPECT_NEAR(props[0].value, 0.7, 0.25);   // the film's own density
  // the interface: the film reaches half its plateau a few Å above the surface; the adsorbed layer is at least that thick
  const double gap = props[0].extra.at("gap to the surface (Å)");
  EXPECT_GT(gap, 0.5);
  EXPECT_LT(gap, 8.0);
  EXPECT_GE(props[0].extra.at("adsorbed layer thickness (Å)"), gap);
  EXPECT_NEAR(props[0].extra.at("film reaches half its plateau at z (Å)") - props[0].extra.at("surface top (Å)"), gap, 1e-9);
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

namespace {
caps::System scrambled(const caps::System& s) {
  // a unimodular change of basis (det +1): the same lattice in an ugly cell
  const caps::Mat3 u{{{1, 1, 0}, {0, 1, 1}, {1, 1, 1}}};   // det = 1
  return caps::transform_cell(s, u);
}
caps::System fcc_copper() {
  caps::CrystalSpec c;
  c.space_group = "F m -3 m";
  c.a = c.b = c.c = 3.615;
  c.sites = {{"Cu1", 29, {0, 0, 0}}};
  return caps::build_crystal(c);
}
}  // namespace

TEST(Lattice, NiggliAndConventionalCellOfCopper) {
  const caps::System conv = fcc_copper();
  ASSERT_EQ(conv.atoms.size(), 4u);
  const caps::System prim = caps::primitive_cell(conv, 'F');
  ASSERT_EQ(prim.atoms.size(), 1u);
  const caps::System ugly = scrambled(prim);
  ASSERT_EQ(ugly.atoms.size(), 1u);
  caps::NiggliResult nr;
  const caps::System red = caps::niggli_cell(ugly, &nr);
  const double ap = 3.615 / std::sqrt(2.0);
  EXPECT_NEAR(nr.a, ap, 1e-6);
  EXPECT_NEAR(nr.b, ap, 1e-6);
  EXPECT_NEAR(nr.c, ap, 1e-6);
  // the Niggli cell of FCC: all angles 60° (type I, +++)
  EXPECT_NEAR(nr.alpha, 60, 1e-6);
  EXPECT_NEAR(nr.beta, 60, 1e-6);
  EXPECT_NEAR(nr.gamma, 60, 1e-6);
  caps::ConventionalResult cr;
  const caps::System back = caps::conventional_cell(ugly, 0.1, &cr);
  EXPECT_EQ(cr.number, 225) << cr.hm;
  EXPECT_EQ(back.atoms.size(), 4u);
  EXPECT_NEAR(caps::norm(back.cell.a), 3.615, 1e-6);
  EXPECT_NEAR(back.cell.volume(), 3.615 * 3.615 * 3.615, 1e-4);
}

TEST(Lattice, ConventionalCellsOfRutileAndWurtzite) {
  caps::CrystalSpec r;
  r.space_group = "P 42/m n m";
  r.a = r.b = 4.594, r.c = 2.959;
  r.sites = {{"Ti1", 22, {0, 0, 0}}, {"O1", 8, {0.3048, 0.3048, 0}}};
  caps::ConventionalResult cr;
  const caps::System rut = caps::conventional_cell(scrambled(caps::build_crystal(r)), 0.1, &cr);
  EXPECT_EQ(cr.number, 136) << cr.hm;
  EXPECT_EQ(rut.atoms.size(), 6u);
  caps::CrystalSpec w;
  w.space_group = "P 63 m c";
  w.a = w.b = 3.25, w.c = 5.207, w.gamma = 120;
  w.sites = {{"Zn1", 30, {1.0 / 3, 2.0 / 3, 0}}, {"O1", 8, {1.0 / 3, 2.0 / 3, 0.382}}};
  const caps::System wz = caps::conventional_cell(scrambled(caps::build_crystal(w)), 0.1, &cr);
  EXPECT_EQ(cr.number, 186) << cr.hm;
  EXPECT_EQ(wz.atoms.size(), 4u);
  // a supercell comes back primitive (rutile's cell is primitive: 6 atoms), then to the conventional cell
  int k = 0;
  const caps::System prim = caps::find_primitive_cell(caps::supercell(caps::build_crystal(r), 2, 1, 3), 0.1, &k);
  EXPECT_EQ(k, 6);
  EXPECT_EQ(prim.atoms.size(), 6u);
  EXPECT_EQ(caps::find_primitive_cell(fcc_copper()).atoms.size(), 1u);   // F-centred: 4 lattice points
  caps::conventional_cell(caps::supercell(caps::build_crystal(r), 2, 2, 1), 0.1, &cr);
  EXPECT_EQ(cr.number, 136) << cr.hm;
  // a matrix that does not map the lattice onto itself is refused
  EXPECT_THROW(caps::transform_cell(caps::build_crystal(r), caps::Mat3{{{0.5, 0, 0}, {0, 1, 0}, {0, 0, 1}}}), std::invalid_argument);
}

TEST(Lattice, VacuumSlabAndNanowire) {
  const caps::System cu = caps::supercell(fcc_copper(), 2, 2, 3);
  caps::SlabResult sr;
  const caps::System slab = caps::vacuum_slab(cu, 15, true, &sr);
  EXPECT_EQ(slab.atoms.size(), cu.atoms.size());
  EXPECT_NEAR(sr.thickness, 3 * 3.615 - 3.615 / 2, 1e-6);   // six (002) layers, 1.8075 Å apart
  EXPECT_NEAR(caps::norm(slab.cell.c), sr.thickness + 15, 1e-6);
  double lo = 1e9, hi = -1e9;
  for (const auto& a : slab.atoms) lo = std::min(lo, a.pos[2]), hi = std::max(hi, a.pos[2]);
  EXPECT_NEAR(lo, 7.5, 1e-6);   // centred: 7.5 Å of vacuum each side
  EXPECT_NEAR(caps::norm(slab.cell.c) - hi, 7.5, 1e-6);
  caps::WireOptions o;
  o.uvw = {0, 0, 1};
  o.radius = 6;
  o.repeats = 2;
  caps::WireResult wr;
  const caps::System wire = caps::nanowire(fcc_copper(), o, &wr);
  EXPECT_NEAR(wr.period, 3.615, 1e-9);
  EXPECT_NEAR(caps::norm(wire.cell.c), 2 * 3.615, 1e-9);
  // the atom count follows the density: 4 atoms per a³ in the cylinder, within the edge's granularity
  const double expect = 4 / std::pow(3.615, 3) * M_PI * 36 * 2 * 3.615;
  EXPECT_NEAR(double(wire.atoms.size()), expect, 0.15 * expect);
  for (const auto& a : wire.atoms) {
    const double dx = a.pos[0] - caps::norm(wire.cell.a) / 2, dy = a.pos[1] - caps::norm(wire.cell.b) / 2;
    EXPECT_LE(std::sqrt(dx * dx + dy * dy), 6 + 1e-9);
  }
}
