#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>

#include "caps/provenance.hpp"

using namespace caps;

namespace {

Manifest sample(uint64_t seed) {
  Manifest m;
  m.inputs = {{"PS_atactic_DP40.caps", "9f2c"}, {"gaff2.json", "71be"}};
  ProvStep g;
  g.engine = "grow.trials";
  g.params = {{"chains", "20"}, {"DP", "40"}};
  g.rng = "mt19937-64 · seed " + std::to_string(seed);
  g.cites = {"matsumoto1998"};
  ProvStep r;
  r.engine = "relax.lbfgs";
  r.params = {{"|F|max", "0.02 kcal/mol/Å"}};
  r.cites = {"liu1989", "auhl2003"};
  r.approximations = {{"van der Waals", "cut-off 12 Å"}, {"Electrostatics", "off"}};
  ProvStep md;
  md.engine = "dynamics.npt";
  md.cites = {"swope1982", "bussi2007", "bernetti2020", "essmann1995"};
  md.approximations = {{"van der Waals", "cut-off 12 Å + tail correction"}, {"Electrostatics", "SPME · relative tolerance 1e-05"}, {"Precision", "double · reproducible with 10 threads"}};
  m.steps = {g, r, md};
  return m;
}

}  // namespace

TEST(Provenance, JsonRoundTripAndSidecar) {
  const Manifest m = sample(20260923);
  const Manifest back = manifest_from_json(Json::parse(manifest_json(m).dump(2)));
  ASSERT_EQ(back.steps.size(), 3u);
  EXPECT_EQ(back.steps[0].rng, "mt19937-64 · seed 20260923");
  EXPECT_EQ(back.steps[1].cites, (std::vector<std::string>{"liu1989", "auhl2003"}));
  EXPECT_EQ(back.inputs, m.inputs);

  const auto data = (std::filesystem::temp_directory_path() / "caps_prov_cell.data").string();
  write_manifest(m, data);
  EXPECT_TRUE(std::filesystem::exists(data + ".provenance.json"));
  const auto r = read_manifest(data);
  ASSERT_TRUE(r.has_value());
  EXPECT_EQ(r->steps[2].engine, "dynamics.npt");
  EXPECT_FALSE(read_manifest(data + ".missing").has_value());
}

TEST(Provenance, ApproximationsTakeTheLatestStep) {
  const auto a = approximations(sample(1));
  ASSERT_GE(a.size(), 3u);
  EXPECT_EQ(a[0].first, "van der Waals");
  EXPECT_EQ(a[0].second, "cut-off 12 Å + tail correction");
  EXPECT_EQ(a[1].second, "SPME · relative tolerance 1e-05");
  EXPECT_EQ(a[2].first, "Precision");
}

TEST(Provenance, CompareFindsOnlyTheSeed) {
  const auto d = compare(sample(20260923), sample(20260924));
  ASSERT_EQ(d.rows.size(), 1u);
  EXPECT_EQ(d.rows[0].engine, "grow.trials");
  EXPECT_EQ(d.rows[0].key, "rng");
  EXPECT_EQ(d.differing_steps, 1);
  EXPECT_TRUE(d.same_inputs);
  EXPECT_TRUE(d.notes.empty());

  Manifest b = sample(20260923);
  b.steps.pop_back();
  b.inputs[1].second = "0000";
  const auto e = compare(sample(20260923), b);
  EXPECT_FALSE(e.same_inputs);
  EXPECT_EQ(e.notes.size(), 2u);   // inputs differ; dynamics only in the first
}

TEST(Provenance, BibtexForEveryCitedMethod) {
  const std::string bib = bibtex(all_cites(sample(1)));
  for (const char* k : {"matsumoto1998", "liu1989", "auhl2003", "swope1982", "bussi2007", "bernetti2020", "essmann1995"})
    EXPECT_NE(bib.find(std::string("{") + k + ","), std::string::npos) << k;
  EXPECT_NE(bib.find("doi = {10.1063/1.2408420}"), std::string::npos);
  // every key the C API cites is in the table
  for (const char* k : {"abascal2005", "auhl2003", "berendsen1984", "berendsen1987", "bernetti2020", "bitzek2006", "bussi2007", "cordero2008",
                        "engh1991", "essmann1995", "fennell2006", "gasteiger1980", "hall1981", "jorgensen1983", "larsen2011", "liu1989",
                        "martinez2009", "matsumoto1998", "parsons2005", "polak1969", "rappe1991", "rappe1992", "swope1982", "wang2004"})
    EXPECT_TRUE(known_citation(k)) << k;
}

// Neutron contrast by deuteration (design/boards/Scattering): PS has aliphatic backbone H and aromatic ring H.
#include "caps/io.hpp"
#include "caps/properties.hpp"

TEST(Scattering, DeuterationPatternsSplitTheHydrogens) {
  const Trajectory t = open_file(std::string(CAPS_SAMPLES) + "/ps_melt.data");
  size_t nh = 0;
  for (const auto& a : t.topology.atoms) nh += a.element == 1;
  auto count = [&](int p) { size_t n = 0; for (char c : deuterated_hydrogens(t.topology, p)) n += c; return n; };
  EXPECT_EQ(count(1), nh);
  EXPECT_EQ(count(2) + count(3), nh);   // every PS hydrogen is on an aliphatic or an aromatic carbon
  EXPECT_EQ(count(3) % 5, 0u);          // five ring H on every phenyl
  EXPECT_GT(count(3), count(2));        // C8H8: five ring H, three backbone H per unit (plus end groups)
  EXPECT_EQ(count(4), 0u);

  AnalyzeOptions o;
  o.q_direct = 2.0;
  o.qmax = 6.0;
  const auto h = analyze(t, {"neutron"}, o);
  o.deuterate = 1;
  const auto d = analyze(t, {"neutron"}, o);
  ASSERT_EQ(h.size(), 1u);
  ASSERT_EQ(d.size(), 1u);
  EXPECT_EQ(d[0].extra.at("deuterated hydrogens"), double(nh));
  EXPECT_NE(d[0].notes.front().find("scatter as ²H"), std::string::npos);
  // deuterium changes the contrast: the curves differ
  ASSERT_FALSE(h[0].series.empty());
  double diff = 0;
  for (size_t k = 0; k < std::min(h[0].series[0].y.size(), d[0].series[0].y.size()); ++k) diff += std::fabs(h[0].series[0].y[k] - d[0].series[0].y[k]);
  EXPECT_GT(diff, 1.0);
}

#include "caps/voids.hpp"

TEST(Voids, OneAtomLeavesTheFarCornerEmpty) {
  System s;
  s.cell.a = {20, 0, 0};
  s.cell.b = {0, 20, 0};
  s.cell.c = {0, 0, 20};
  Atom a;
  a.element = 6;
  s.atoms.push_back(a);
  VoidOptions o;
  o.grid = 0.5;
  o.reach = 25;
  o.max_count = 3;
  const auto r = largest_voids(s, o);
  ASSERT_FALSE(r.spheres.empty());
  const double far = std::sqrt(300.0) - 1.7;   // the body centre, minus Bondi C
  EXPECT_NEAR(r.largest, far, 0.5);
  for (int k = 0; k < 3; ++k) EXPECT_NEAR(r.spheres[0].centre[k], 10.0, 0.6);
  const double vdw = 4.0 / 3.0 * 3.14159265 * 1.7 * 1.7 * 1.7 / 8000.0;
  EXPECT_NEAR(r.accessible_point, 1 - vdw, 2e-3);
  EXPECT_LT(r.accessible_probe, r.accessible_point);
  const Mesh m = void_mesh({r.spheres[0]}, 2);
  EXPECT_EQ(m.vertices.size(), 162u);
  EXPECT_EQ(m.triangles.size(), 320u);
  const auto pdb = (std::filesystem::temp_directory_path() / "caps_voids.pdb").string();
  write_voids_pdb(s, r.spheres, pdb);
  std::ifstream f(pdb);
  std::string all((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
  EXPECT_NE(all.find("CRYST1   20.000   20.000   20.000  90.00  90.00  90.00"), std::string::npos);
  EXPECT_NE(all.find("HETATM    1  VO  VOI"), std::string::npos);
}

#include "caps/crystal.hpp"
#include "caps/molecule.hpp"
#include "caps/nano.hpp"

TEST(Pore, SlitKeepsTheFluidBetweenTheWalls) {
  BuildOptions bo;
  bo.forcefield = "uff";
  System methane = build_molecule("C", bo).system;
  methane.title = "methane";
  PoreOptions o;
  o.width = 10.0;
  o.fluid = &methane;
  o.count = 12;
  PoreReport r;
  const System s = build_pore(o, &r);
  ASSERT_GT(r.wall_atoms, 100);
  EXPECT_EQ(int(s.atoms.size()), r.wall_atoms + 12 * 5);
  EXPECT_NEAR(s.cell.c[2], 10.0 + 3.35, 1e-9);   // one sheet per wall: the walls meet through the boundary
  for (size_t i = size_t(r.wall_atoms); i < s.atoms.size(); ++i) {
    EXPECT_GT(s.atoms[i].pos[2], 0.5);
    EXPECT_LT(s.atoms[i].pos[2], 9.5);
    EXPECT_GT(s.atoms[i].mol, 1);
  }
  for (int i = 0; i < r.wall_atoms; ++i) EXPECT_EQ(s.atoms[size_t(i)].mol, 1);
  EXPECT_GE(r.dmin, 1.9);
  EXPECT_GT(r.fluid_density, 0);
}

TEST(Pore, ChannelInQuartzIsEmptyAlongItsAxis) {
  const System quartz = read_cif(std::string(CAPS_SOURCE_DIR) + "/data/crystals/alpha-quartz.cif");
  BuildOptions bo;
  bo.forcefield = "uff";
  System water = build_molecule("O", bo).system;
  PoreOptions o;
  o.kind = PoreKind::Cylinder;
  o.crystal = &quartz;
  o.width = 14.0;
  o.wall = 5.0;
  o.length = 15.0;
  o.fluid = &water;
  o.count = 8;
  PoreReport r;
  const System s = build_pore(o, &r);
  const Vec3 axis = s.cell.origin + (s.cell.a + s.cell.b) * 0.5;
  for (int i = 0; i < r.wall_atoms; ++i) {
    Vec3 d = s.cell.minimum_image(s.atoms[size_t(i)].pos - axis);
    d[2] = 0;
    EXPECT_GE(norm(d), 7.0 - 1e-6);
  }
  for (size_t i = size_t(r.wall_atoms); i < s.atoms.size(); ++i) {
    Vec3 d = s.cell.minimum_image(s.atoms[i].pos - axis);
    d[2] = 0;
    EXPECT_LT(norm(d), 7.0);
  }
  EXPECT_EQ(int(s.atoms.size()), r.wall_atoms + 8 * 3);
}

#include "caps/grow.hpp"

TEST(Gasteiger, HeteroatomsGetTheirOwnParameters) {
  BuildOptions bo;
  bo.forcefield = "uff";
  const System ethanol = build_molecule("CCO", bo).system;
  const auto q = gasteiger_ch(ethanol, {});
  double qo = 0, qho = 0, total = 0;
  for (size_t i = 0; i < ethanol.atoms.size(); ++i) {
    total += q[i];
    if (ethanol.atoms[i].element == 8) qo = q[i];
  }
  for (const auto& b : ethanol.bonds) {
    const auto& ai = ethanol.atoms[b.i];
    const auto& aj = ethanol.atoms[b.j];
    if (ai.element == 8 && aj.element == 1) qho = q[b.j];
    if (aj.element == 8 && ai.element == 1) qho = q[b.i];
  }
  EXPECT_NEAR(total, 0.0, 1e-9);
  EXPECT_LT(qo, -0.35);    // PEOE ethanol O ≈ −0.39
  EXPECT_GT(qho, 0.18);    // hydroxyl H ≈ +0.21
  const System nitrile = build_molecule("CC#N", bo).system;   // acrylonitrile units in NBR carry C≡N
  const auto qn = gasteiger_ch(nitrile, {});
  for (size_t i = 0; i < nitrile.atoms.size(); ++i)
    if (nitrile.atoms[i].element == 7) EXPECT_NEAR(qn[i], -0.197, 0.01);   // PEOE acetonitrile N
  const System chloro = build_molecule("C=CCl", bo).system;   // chloroprene-like vinyl chloride
  const auto qc = gasteiger_ch(chloro, {});
  for (size_t i = 0; i < chloro.atoms.size(); ++i)
    if (chloro.atoms[i].element == 17) EXPECT_LT(qc[i], -0.05);
  const System silane = build_molecule("[SiH4]", bo).system;
  EXPECT_THROW(gasteiger_ch(silane, {}), std::invalid_argument);
}

TEST(Provenance, CitationTextReadsTheEntries) {
  EXPECT_EQ(citation_text("bussi2007"), "Bussi, G., Donadio, D., Parrinello, M., \"Canonical sampling through velocity rescaling\", J. Chem. Phys. 126, 014101 (2007). doi:10.1063/1.2408420");
  EXPECT_NE(citation_text("cordero2008").find("Gómez"), std::string::npos);
  EXPECT_NE(citation_text("einstein1905").find("Über"), std::string::npos);
  EXPECT_NE(citation_text("prince2004").find("(ed.)"), std::string::npos);
  EXPECT_EQ(citation_text("nobody2099"), "");
}

TEST(Provenance, MethodsTextNumbersItsReferences) {
  Manifest m;
  ProvStep r;
  r.engine = "relax.lbfgs";
  r.params = {{"minimiser", "L-BFGS"}, {"|F|max", "0.5 kcal/mol/Å"}, {"push-off", "on"}};
  r.cites = {"liu1989", "auhl2003", "fennell2006"};
  r.approximations = {{"van der Waals", "cut-off 12 Å + tail correction"}, {"Electrostatics", "damped shifted force · cut-off 12 Å"}};
  ProvStep d;
  d.engine = "dynamics.npt";
  d.params = {{"length", "1000 ps · 1000000 steps of 1 fs"}, {"temperature", "300 K"}, {"thermostat", "Bussi velocity rescaling · τ 100 fs"},
              {"barostat", "stochastic cell rescaling · 1 atm · τ 1000 fs"}};
  d.cites = {"swope1982", "bussi2007", "bernetti2020", "fennell2006"};
  m.steps = {r, d};
  std::vector<std::string> refs;
  const std::string t = methods_text(m, &refs, {m});
  EXPECT_NE(t.find("minimised with L-BFGS [1]"), std::string::npos) << t;
  EXPECT_NE(t.find("push-off stages [2]"), std::string::npos) << t;
  EXPECT_NE(t.find("(velocity Verlet [3], 1 fs time step) was run for 1000 ps at 300 K"), std::string::npos) << t;
  EXPECT_NE(t.find("truncated at 12 Å with analytic tail corrections"), std::string::npos) << t;
  EXPECT_NE(t.find("damped shifted force method [6] (cut-off 12 Å)"), std::string::npos) << t;
  EXPECT_NE(t.find("2 independent replicas"), std::string::npos) << t;
  ASSERT_EQ(refs.size(), 6u);
  EXPECT_EQ(refs[0].rfind("Liu, D. C.", 0), 0u);
}

#include "caps/kremer_grest.hpp"

TEST(KremerGrest, ChainsBoxAndDeck) {
  KgOptions o;
  o.chains = 10;
  o.beads = 30;
  o.k_theta = 1.5;
  KgReport r;
  const System s = kremer_grest(o, &r);
  EXPECT_EQ(s.atoms.size(), 300u);
  EXPECT_EQ(s.bonds.size(), 290u);
  EXPECT_NEAR(r.box, std::cbrt(300 / 0.85), 1e-9);
  for (const auto& b : s.bonds) {
    const Vec3 d = s.atoms[b.j].pos - s.atoms[b.i].pos;
    EXPECT_NEAR(std::sqrt(dot(d, d)), 0.97, 1e-9);
  }
  for (size_t i = 1; i + 1 < s.atoms.size(); ++i) {
    if (s.atoms[i - 1].mol != s.atoms[i].mol || s.atoms[i + 1].mol != s.atoms[i].mol) continue;
    const Vec3 a = s.atoms[i].pos - s.atoms[i - 1].pos, b = s.atoms[i + 1].pos - s.atoms[i].pos;
    EXPECT_GE(dot(a, b) / (0.97 * 0.97), -0.5 - 1e-9);   // no folding back past 120° between bonds
  }
  EXPECT_GT(r.mean_r2, 1.0);   // stiffer than a freely jointed chain
  const auto stem = (std::filesystem::temp_directory_path() / "caps_kg").string();
  write_kg_lammps(s, o, stem, 100, 100);
  std::ifstream in(stem + ".in");
  const std::string deck((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
  for (const char* key : {"units lj", "bond_coeff 1 30.0 1.5 1.0 1.0", "special_bonds fene", "pair_style soft", "pair_modify shift yes", "angle_style cosine"})
    EXPECT_NE(deck.find(key), std::string::npos) << key;
  std::ifstream dat(stem + ".data");
  const std::string data((std::istreambuf_iterator<char>(dat)), std::istreambuf_iterator<char>());
  EXPECT_NE(data.find("300 atoms"), std::string::npos);
  EXPECT_NE(data.find("280 angles"), std::string::npos);
}
