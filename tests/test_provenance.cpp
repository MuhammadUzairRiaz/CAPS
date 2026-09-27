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

#include "caps/react.hpp"

TEST(ReactTemplate, ViewOfTheEpoxyAmineTemplate) {
  const auto t = parse_templates(builtin_template("epoxy_amine_primary"));
  ASSERT_FALSE(t.empty());
  const Json v = Json::parse(template_view(t[0]));
  EXPECT_GE(v["pre"]["atoms"].size(), 4u);
  EXPECT_GE(v["changes"].size(), 2u);
  bool formed = false, broken = false;
  for (const auto& c : v["changes"].items()) { formed |= c["kind"].str() == "formed"; broken |= c["kind"].str() == "broken"; }
  EXPECT_TRUE(formed);
  EXPECT_TRUE(broken);
  for (const auto& c : v["checks"].items()) EXPECT_TRUE(c["ok"].boolean()) << c["text"].str();
  // a broken template: breaking a bond that is not there
  auto bad = t[0];
  bad.brk.push_back({bad.atoms.front().map, bad.atoms.back().map});
  const Json w = Json::parse(template_view(bad));
  bool flagged = false;
  for (const auto& c : w["checks"].items()) flagged |= !c["ok"].boolean();
  EXPECT_TRUE(flagged);
}

#include "caps/yaml.hpp"

TEST(Yaml, TheRecipeSubset) {
  const Json j = yaml_parse(R"(recipe: 1   # a CAPS recipe
name: ps_cell
build:
  polymer:
    smiles: "*CC(*)c1ccccc1"
    dp: 40
    chains: 20
    tacticity: atactic
type: { forcefield: gaff2, charges: types }
grow: { density: 0.50 }
relax: { method: lbfgs, fmax: 0.02 }
export: [lammps, gromacs]
steps:
  - md: { ps: 10, temperature: 300 }
  - analyze:
      properties: [density, rg]
units: '[*:1]CC[*:2]'
)");
  EXPECT_EQ(j["recipe"].number(), 1);
  EXPECT_EQ(j["build"]["polymer"]["smiles"].str(), "*CC(*)c1ccccc1");
  EXPECT_EQ(j["build"]["polymer"]["dp"].number(), 40);
  EXPECT_EQ(j["type"]["forcefield"].str(), "gaff2");
  EXPECT_DOUBLE_EQ(j["grow"]["density"].number(), 0.5);
  EXPECT_EQ(j["export"].size(), 2u);
  EXPECT_EQ(j["export"][1].str(), "gromacs");
  EXPECT_EQ(j["steps"].size(), 2u);
  EXPECT_EQ(j["steps"][0]["md"]["ps"].number(), 10);
  EXPECT_EQ(j["steps"][1]["analyze"]["properties"][1].str(), "rg");
  EXPECT_EQ(j["units"].str(), "[*:1]CC[*:2]");
  EXPECT_THROW(yaml_parse("a: 1\n  b: 2\n"), std::invalid_argument);
  EXPECT_THROW(yaml_parse("a: 1\na: 2\n"), std::invalid_argument);
  EXPECT_THROW(yaml_parse("a: {b: 1\n"), std::invalid_argument);
}

#include "caps/recipe.hpp"
#include "caps/yaml.hpp"

TEST(Recipe, BuildsGrowsRelaxesAndExportsWithProvenance) {
  const auto dir = std::filesystem::temp_directory_path() / "caps_recipe_test";
  std::filesystem::remove_all(dir);
  const caps::Json r = caps::yaml_parse(
      "recipe: 1\nname: pe\nbuild:\n  polymer: { smiles: \"*CC*\", dp: 5, chains: 3 }\ntype: { forcefield: default }\n"
      "grow: { density: 0.5, seed: 2 }\nrelax: { fmax: 5 }\nanalyze: { properties: [density] }\nexport: [lammps, pdb]\n");
  EXPECT_EQ(caps::recipe_stages(r), (std::vector<std::string>{"build", "type", "grow", "relax", "analyze", "export"}));
  caps::RecipeOptions o;
  o.out_dir = dir.string();
  std::vector<std::string> seen;
  o.progress = [&](const caps::RecipeEvent& e) { if (e.status == "done") seen.push_back(e.name); };
  const auto res = caps::run_recipe(r, o);
  EXPECT_EQ(seen.size(), 6u);
  EXPECT_EQ(res.system.atoms.size(), 3u * (5 * 6 + 2));
  ASSERT_FALSE(res.properties.empty());
  EXPECT_NEAR(res.properties[0].value, 0.5, 0.05);
  EXPECT_TRUE(std::filesystem::exists(dir / "pe.data"));
  EXPECT_TRUE(std::filesystem::exists(dir / "pe.in"));
  EXPECT_TRUE(std::filesystem::exists(dir / "pe.pdb"));
  const auto m = caps::read_manifest((dir / "pe.data").string());
  ASSERT_TRUE(m.has_value());
  EXPECT_NE(caps::methods_text(*m).find("L-BFGS"), std::string::npos);
  // the same seed gives the same cell
  const auto again = caps::run_recipe(r, o);
  EXPECT_DOUBLE_EQ(again.system.atoms.back().pos[0], res.system.atoms.back().pos[0]);
  std::filesystem::remove_all(dir);
}

// A sulfur cure in a recipe: H–S–S–H donors packed into a natural-rubber cell, cured, and the network typed again with
// PCFF (every crosslink atom typed, no missing term); an unknown template is an input error
TEST(Recipe, SulfurCureTypedWithPcff) {
  const caps::Json r = caps::yaml_parse(
      "recipe: 1\nname: nr\nbuild:\n  polymer: { smiles: \"[*]C/C=C(C)\\\\C[*]\", dp: 10, chains: 3 }\ntype: { forcefield: " + std::string(CAPS_SOURCE_DIR) + "/data/forcefields/pcff-frc.json }\n"
      "grow: { density: 0.5, seed: 2 }\nreact:\n  insert: { smiles: SS, count: 6 }\n  templates: [sulfur_allylic]\n  relax: false\n  seed: 3\n");
  EXPECT_EQ(caps::recipe_stages(r), (std::vector<std::string>{"build", "type", "grow", "react"}));
  caps::RecipeOptions o;
  std::string react_summary;
  o.progress = [&](const caps::RecipeEvent& e) { if (e.status == "done" && e.name == "react") react_summary = e.detail; };
  const auto res = caps::run_recipe(r, o);
  EXPECT_NE(react_summary.find("typed again: PCFF"), std::string::npos) << react_summary;
  int ss = 0, cs = 0;
  for (const auto& b : res.system.bonds) {
    const int ei = res.system.atoms[b.i].element, ej = res.system.atoms[b.j].element;
    ss += ei == 16 && ej == 16;
    cs += (ei == 16 && ej == 6) || (ei == 6 && ej == 16);
  }
  EXPECT_EQ(ss, 6);    // every donor keeps its S–S bond
  EXPECT_GT(cs, 0);    // and some are bonded to the rubber
  EXPECT_THROW(caps::run_recipe(caps::yaml_parse("build: {polymer: {smiles: \"*CC*\", dp: 3, chains: 1}}\ngrow: {density: 0.3}\nreact: {templates: [nope]}\n"), o),
               caps::RecipeError);
}

// type.fill_from: OPLS-AA 2024 lacks the CM–CT–CT–CM torsion of polyisoprene; OPLS 2005's fills it (and only where
// 2024 has none), the provenance says how many terms; without it the recipe stops with missing parameters (exit 3)
TEST(Recipe, FillGapsFromAnotherForceField) {
  const std::string lib = std::string(CAPS_SOURCE_DIR) + "/data/forcefields/";
  auto yaml = [&](const std::string& fill) {
    return caps::yaml_parse("build:\n  polymer: { smiles: \"[*]C/C=C(C)\\\\C[*]\", dp: 6, chains: 1 }\ntype: { forcefield: " + lib + "oplsaa2024-moltemplate.json" + fill + " }\n");
  };
  caps::RecipeOptions o;
  int code = 0;
  try { caps::run_recipe(yaml(""), o); } catch (const caps::RecipeError& e) { code = e.code; }
  EXPECT_EQ(code, 3);
  const auto res = caps::run_recipe(yaml(", fill_from: " + lib + "opls2005.json"), o);
  bool said = false;
  for (const auto& st : res.manifest.steps)
    if (st.engine == "field.assign")
      for (const auto& [k, v] : st.params) said |= k == "gaps filled from" && v.find("OPLS 2005") != std::string::npos && v.find(" terms") != std::string::npos;
  EXPECT_TRUE(said);
}

TEST(Recipe, ExitCodes) {
  caps::RecipeOptions o;
  auto code = [&](const std::string& y) {
    try { caps::run_recipe(caps::yaml_parse(y), o); } catch (const caps::RecipeError& e) { return e.code; }
    return 0;
  };
  EXPECT_EQ(code("build: {molecule: CCO}\nbogus: 1\n"), 2);
  EXPECT_EQ(code("type: {forcefield: uff}\n"), 2);
  EXPECT_EQ(code("build: {molecule: CCO}\ngrow: {density: 0.5}\n"), 2);
  EXPECT_EQ(code("build: {polymer: {smiles: \"*CC*\", dp: 2, chains: 1}}\ntype: {forcefield: nothing-here.json}\n"), 2);
  EXPECT_EQ(code("build: {molecule: CCO}\ntype: {forcefield: uff}\nrelax: {fmax: 2}\n"), 0);
}

#include "caps/colourvision.hpp"

TEST(ColourVision, MachadoSimulationAndDeltaE) {
  using caps::Vision;
  // greys are unchanged (each Machado row sums to 1); black–white is ΔE 100
  for (Vision v : {Vision::Protan, Vision::Deutan, Vision::Tritan}) {
    const unsigned g = caps::simulate_vision(0x808080, v);
    EXPECT_NEAR(double((g >> 16) & 255), 128.0, 1.0);
    EXPECT_NEAR(double(g & 255), 128.0, 1.0);
  }
  EXPECT_NEAR(caps::delta_e76(0x000000, 0xFFFFFF), 100.0, 1e-4);
  EXPECT_EQ(caps::simulate_vision(0xE35049, Vision::Normal), 0xE35049u);
  // a red and a green of similar lightness stay far apart normally and collapse for deuteranopes
  const unsigned red = 0xC0504D, green = 0x6E9A3C;
  EXPECT_GT(caps::delta_e76(red, green), 40);
  EXPECT_LT(caps::delta_e76(caps::simulate_vision(red, Vision::Deutan), caps::simulate_vision(green, Vision::Deutan)), 15);
  const auto pairs = caps::confusable_pairs({{"test", {"red", "green", "blue"}, {red, green, 0x2271DB}}}, 15);
  ASSERT_FALSE(pairs.empty());
  EXPECT_EQ(pairs[0].a, 0);
  EXPECT_EQ(pairs[0].b, 1);
  // an image is simulated in place, alpha kept
  caps::Image img;
  img.width = 1; img.height = 1;
  img.rgba = {0xC0, 0x50, 0x4D, 77};
  caps::simulate_vision(img, Vision::Protan);
  EXPECT_EQ(img.rgba[3], 77);
  EXPECT_EQ((unsigned(img.rgba[0]) << 16) | (unsigned(img.rgba[1]) << 8) | img.rgba[2], caps::simulate_vision(red, Vision::Protan));
}

#include "caps/molecule.hpp"
#include "caps/query.hpp"
#include "caps/polymer.hpp"

TEST(Query, CipLabelsOfAlanine) {
  caps::BuildOptions b;
  b.forcefield = "uff";
  const auto l = caps::build_molecule("N[C@@H](C)C(=O)O", b);   // L-alanine: S
  const auto d = caps::build_molecule("N[C@H](C)C(=O)O", b);    // D-alanine: R
  const auto cl = caps::cip_labels(l.system), cd = caps::cip_labels(d.system);
  EXPECT_EQ(std::count(cl.begin(), cl.end(), 'S'), 1);
  EXPECT_EQ(std::count(cl.begin(), cl.end(), 'R'), 0);
  EXPECT_EQ(std::count(cd.begin(), cd.end(), 'R'), 1);
  // the methyl carbon (three H) and the carboxyl carbon are not stereocentres
  EXPECT_EQ(std::count_if(cl.begin(), cl.end(), [](char c) { return c != 0; }), 1);
  EXPECT_EQ(caps::select_query(l.system, "stereo S").atoms[1], 1);
}

TEST(Query, GrammarOnAPolystyreneChain) {
  caps::ChainSpec spec;
  spec.units = {{"A", "*CC(*)c1ccccc1"}};
  spec.dp = 10;
  caps::GrowOptions g;
  g.chains = 1;
  g.density = 0.05;
  const caps::System s = caps::grow_chains(spec, g, nullptr);
  auto count = [&](const std::string& q) { const auto r = caps::select_query(s, q); return std::count(r.atoms.begin(), r.atoms.end(), 1); };
  const auto rings = caps::select_query(s, "smarts \"c1ccccc1\"");
  EXPECT_EQ(std::count(rings.atoms.begin(), rings.atoms.end(), 1), 60);
  EXPECT_EQ(rings.rings, 10);
  EXPECT_EQ(count("element H"), 8 * 10 + 2);
  // the backbone CH of every unit but the last, whose H cap leaves two hydrogens on it
  EXPECT_EQ(count("stereo *"), 9);
  EXPECT_EQ(count("ring 1"), 6);
  EXPECT_EQ(count("smarts \"c1ccccc1\" and chain 1"), 60);
  EXPECT_EQ(count("not element H"), long(s.atoms.size()) - 82);
  EXPECT_EQ(count("(element C or element H) and not smarts \"c\""), long(s.atoms.size()) - 60);
  const long near = count("within 5 of ring 1 and not element H");
  EXPECT_GT(near, 6);
  EXPECT_LT(near, count("within 5 of ring 1"));
  EXPECT_EQ(count("index 1-5"), 5);
  std::vector<char> sel(s.atoms.size(), 0);
  sel[3] = 1;
  const auto only = caps::select_query(s, "sel", sel);
  EXPECT_EQ(std::count(only.atoms.begin(), only.atoms.end(), 1), 1);
  EXPECT_THROW(caps::select_query(s, "within 5 sel"), std::invalid_argument);
  EXPECT_THROW(caps::select_query(s, "element C and"), std::invalid_argument);
  EXPECT_THROW(caps::select_query(s, "colour red"), std::invalid_argument);
}

#include "caps/charges.hpp"

TEST(Charges, GasteigerReportByGroupAndChgFiles) {
  caps::ChainSpec spec;
  spec.units = {{"A", "*CC(*)c1ccccc1"}};
  spec.dp = 8;
  caps::GrowOptions g;
  g.chains = 1;
  g.density = 0.05;
  const caps::System s = caps::grow_chains(spec, g, nullptr);
  const auto r = caps::compute_charges(s, "gasteiger");
  EXPECT_NEAR(r.net, 0.0, 1e-9);
  EXPECT_GT(r.max_abs, 0.03);
  std::map<std::string, int> n;
  for (const auto& gr : r.groups) n[gr.name] = gr.n;
  EXPECT_EQ(n["C aromatic"], 48);
  EXPECT_EQ(n["H on aromatic C"], 40);
  EXPECT_EQ(n["C sp³"], 16);
  EXPECT_EQ(n["H on sp³ C"], 26);
  double counted = 0;
  for (double c : r.counts) counted += c;
  EXPECT_EQ(counted, double(s.atoms.size()));
  // a .chg file: one per line, or "index charge"
  const auto path = (std::filesystem::temp_directory_path() / "caps_test.chg").string();
  { std::ofstream f(path); f << "# test\n"; for (size_t i = 0; i < s.atoms.size(); ++i) f << (i + 1) << " " << (i % 2 ? 0.1 : -0.1) << "\n"; }
  const auto q = caps::read_charge_file(path, s.atoms.size());
  EXPECT_DOUBLE_EQ(q[0], -0.1);
  EXPECT_DOUBLE_EQ(q[1], 0.1);
  { std::ofstream f(path); f << "0.5\n"; }
  EXPECT_THROW(caps::read_charge_file(path, s.atoms.size()), std::runtime_error);
  std::filesystem::remove(path);
}

#include "caps/spacegroup.hpp"
#include "caps/pipeline.hpp"
#include "caps/io.hpp"

TEST(Orientation, PolyethyleneCrystalIsPerfectlyOrdered) {
  // Bunn's orthorhombic PE, 3 × 4 × 8 cells: 24 chains bonded through the cell; the estimator must give S = 1 exactly
  caps::CrystalSpec spec;
  spec.space_group = "Pnam";
  spec.a = 7.40, spec.b = 4.93, spec.c = 2.534;
  spec.sites = {{"C1", 6, {0.0380, 0.0650, 0.25}}, {"H1", 1, {0.1848, 0.0466, 0.25}}, {"H2", 1, {0.0068, 0.2811, 0.25}}};
  spec.supercell = {3, 4, 8};
  const caps::System s = caps::build_crystal(spec);
  caps::Trajectory t;
  t.topology = s;
  std::vector<caps::Vec3> p;
  for (const auto& a : s.atoms) p.push_back(a.pos);
  t.positions.push_back(p);
  t.cells.push_back(s.cell);
  t.timesteps.push_back(0);
  const auto props = caps::analyze(t, {"orientation"}, caps::AnalyzeOptions{});
  ASSERT_EQ(props.size(), 1u);
  EXPECT_NEAR(props[0].value, 1.0, 1e-9);
  EXPECT_EQ(props[0].extra.at("chord vectors per frame"), 336.0);
  EXPECT_NEAR(std::fabs(props[0].extra.at("director z")), 1.0, 1e-9);
  EXPECT_NEAR(props[0].extra.at("local crystallinity (fraction)"), 1.0, 1e-9);
}

// The pipeline step gives the same order per atom: every PE atom crystalline, P₂ = 1 along z; the amorphous melt is
// nearly isotropic; an affine strain scales positions and cell about the centre
TEST(Orientation, PipelineStepPerAtomAndAffineStrain) {
  caps::CrystalSpec spec;
  spec.space_group = "Pnam";
  spec.a = 7.40, spec.b = 4.93, spec.c = 2.534;
  spec.sites = {{"C1", 6, {0.0380, 0.0650, 0.25}}, {"H1", 1, {0.1848, 0.0466, 0.25}}, {"H2", 1, {0.0068, 0.2811, 0.25}}};
  spec.supercell = {3, 4, 8};
  const caps::System s = caps::build_crystal(spec);
  auto run = [](const caps::System& f, const std::string& json) { return caps::run_pipeline(f, caps::pipeline_from_json(caps::Json::parse(json)), 0, 0); };
  auto st = run(s, R"([{"type":"orientation","axis":"z"}])");
  EXPECT_NEAR(st.attribute("Orientation.S"), 1.0, 1e-9);
  EXPECT_NEAR(st.attribute("Orientation.P2_axis"), 1.0, 1e-9);
  EXPECT_NEAR(st.attribute("Crystallinity.fraction"), 1.0, 1e-9);
  for (size_t i = 0; i < s.atoms.size(); ++i) {
    EXPECT_NEAR(st.props.at("Orientation")[i], 1.0, 1e-9) << i;   // hydrogens take their carbon's value
    EXPECT_EQ(st.props.at("Crystalline")[i], 1.0) << i;
  }
  const caps::Trajectory t = caps::open_file(std::string(CAPS_SOURCE_DIR) + "/samples/ps_melt.lammpstrj", std::string(CAPS_SOURCE_DIR) + "/samples/ps_melt.data");
  const caps::System melt = t.frame(0);
  st = run(melt, R"([{"type":"orientation"}])");
  EXPECT_LT(st.attribute("Orientation.S"), 0.3);
  EXPECT_LT(st.attribute("Crystallinity.fraction"), 0.2);
  st = run(melt, R"([{"type":"affine_transform","strain":[0.1,0,0]}])");
  EXPECT_NEAR(st.attribute("AffineTransformation.volume_ratio"), 1.1, 1e-12);
  EXPECT_NEAR(st.system.cell.a[0], 1.1 * melt.cell.a[0], 1e-9);
  const caps::Vec3 c0 = melt.cell.origin + (melt.cell.a + melt.cell.b + melt.cell.c) * 0.5, c1 = st.system.cell.origin + (st.system.cell.a + st.system.cell.b + st.system.cell.c) * 0.5;
  EXPECT_NEAR(caps::norm(c1 - c0), 0.0, 1e-9);
  EXPECT_NEAR(st.system.atoms[0].pos[0] - c0[0], 1.1 * (melt.atoms[0].pos[0] - c0[0]), 1e-9);
}

// Adaptive CNA and centrosymmetry: FCC gold, BCC iron, HCP magnesium each whole; a 13-atom icosahedron's centre;
// centrosymmetry zero in FCC, not at a vacancy's neighbours
TEST(Structure, CommonNeighbourAnalysisAndCentrosymmetry) {
  auto crystal = [](const std::string& sg, double a, double c, std::vector<caps::CrystalSite> sites, std::array<int, 3> sc) {
    caps::CrystalSpec spec;
    spec.space_group = sg;
    spec.a = a, spec.b = a, spec.c = c;
    if (sg == "P63/mmc") spec.gamma = 120;
    spec.sites = std::move(sites);
    spec.supercell = sc;
    return caps::build_crystal(spec);
  };
  auto run = [](const caps::System& f, const std::string& json) { return caps::run_pipeline(f, caps::pipeline_from_json(caps::Json::parse(json)), 0, 0); };
  const auto au = crystal("Fm-3m", 4.078, 4.078, {{"Au", 79, {0, 0, 0}}}, {4, 4, 4});
  auto st = run(au, R"([{"type":"centrosymmetry"},{"type":"cna"}])");
  EXPECT_EQ(st.attribute("CommonNeighborAnalysis.counts.FCC"), double(au.atoms.size()));
  EXPECT_LT(st.attribute("Centrosymmetry.max"), 1e-9);
  const auto fe = crystal("Im-3m", 2.8665, 2.8665, {{"Fe", 26, {0, 0, 0}}}, {5, 5, 5});
  st = run(fe, R"([{"type":"cna"}])");
  EXPECT_EQ(st.attribute("CommonNeighborAnalysis.counts.BCC"), double(fe.atoms.size()));
  const auto mg = crystal("P63/mmc", 3.209, 5.211, {{"Mg", 12, {1.0 / 3, 2.0 / 3, 0.25}}}, {5, 5, 3});
  st = run(mg, R"([{"type":"cna"}])");
  EXPECT_EQ(st.attribute("CommonNeighborAnalysis.counts.HCP"), double(mg.atoms.size()));
  // a vacancy: its twelve neighbours lose their centrosymmetry and their FCC signature
  caps::System vac = au;
  vac.atoms.erase(vac.atoms.begin() + 100);
  st = run(vac, R"([{"type":"centrosymmetry"},{"type":"cna"}])");
  EXPECT_EQ(st.attribute("CommonNeighborAnalysis.counts.FCC"), double(vac.atoms.size() - 12));
  int high = 0;
  for (double v : st.props.at("Centrosymmetry")) high += v > 1.0;
  EXPECT_EQ(high, 12);
  // a Mackay icosahedron of 13 atoms (no cell): the centre is icosahedral
  caps::System ico;
  const double g = (1 + std::sqrt(5.0)) / 2, r = 2.88 / std::sqrt(1 + g * g);
  std::vector<caps::Vec3> v{{0, 0, 0}};
  for (int s1 : {-1, 1})
    for (int s2 : {-1, 1}) v.push_back({0, s1 * r, s2 * g * r}), v.push_back({s1 * r, s2 * g * r, 0}), v.push_back({s2 * g * r, 0, s1 * r});
  for (const auto& x : v) { caps::Atom a; a.element = 79; a.pos = x; a.id = int64_t(ico.atoms.size() + 1); ico.atoms.push_back(a); }
  st = run(ico, R"([{"type":"cna"}])");
  EXPECT_EQ(st.props.at("Structure Type")[0], 4.0);
}

// Polyhedral template matching: FCC gold, BCC iron, HCP magnesium, the icosahedron's centre, each with RMSD ≈ 0; the
// nearest-neighbour distance recovered; a rotated crystal stays FCC without strain; a simple shear γ leaves γ/2
TEST(Structure, PolyhedralTemplateMatching) {
  auto crystal = [](const std::string& sg, double a, double c, std::vector<caps::CrystalSite> sites, std::array<int, 3> sc) {
    caps::CrystalSpec spec;
    spec.space_group = sg;
    spec.a = a, spec.b = a, spec.c = c;
    if (sg == "P63/mmc") spec.gamma = 120;
    spec.sites = std::move(sites);
    spec.supercell = sc;
    return caps::build_crystal(spec);
  };
  auto run = [](const caps::System& f, const std::string& json) { return caps::run_pipeline(f, caps::pipeline_from_json(caps::Json::parse(json)), 0, 0); };
  const auto au = crystal("Fm-3m", 4.078, 4.078, {{"Au", 79, {0, 0, 0}}}, {4, 4, 4});
  auto st = run(au, R"([{"type":"ptm"}])");
  EXPECT_EQ(st.attribute("PolyhedralTemplateMatching.counts.FCC"), double(au.atoms.size()));
  for (size_t i = 0; i < au.atoms.size(); ++i) {
    EXPECT_LT(st.props.at("RMSD")[i], 1e-6);
    EXPECT_NEAR(st.props.at("Interatomic Distance")[i], 4.078 / std::sqrt(2.0), 1e-6);
    EXPECT_LT(st.props.at("Shear Strain")[i], 1e-6);
  }
  const auto fe = crystal("Im-3m", 2.8665, 2.8665, {{"Fe", 26, {0, 0, 0}}}, {5, 5, 5});
  st = run(fe, R"([{"type":"ptm"}])");
  EXPECT_EQ(st.attribute("PolyhedralTemplateMatching.counts.BCC"), double(fe.atoms.size()));
  EXPECT_NEAR(st.props.at("Interatomic Distance")[0], 2.8665 * std::sqrt(3.0) / 2, 1e-6);
  const auto mg = crystal("P63/mmc", 3.209, 5.211, {{"Mg", 12, {1.0 / 3, 2.0 / 3, 0.25}}}, {5, 5, 3});
  st = run(mg, R"([{"type":"ptm"}])");
  EXPECT_EQ(st.attribute("PolyhedralTemplateMatching.counts.HCP"), double(mg.atoms.size()));
  // rotated by 30° about [1 2 3] (cell and all): still FCC, no strain
  {
    caps::System r = au;
    const caps::Vec3 ax = caps::Vec3{1, 2, 3} * (1 / std::sqrt(14.0));
    const double c = std::cos(0.5236), sn = std::sin(0.5236);
    auto rot = [&](const caps::Vec3& p) { return p * c + caps::cross(ax, p) * sn + ax * (caps::dot(ax, p) * (1 - c)); };
    for (auto& a : r.atoms) a.pos = rot(a.pos);
    r.cell.a = rot(r.cell.a), r.cell.b = rot(r.cell.b), r.cell.c = rot(r.cell.c), r.cell.origin = rot(r.cell.origin);
    st = run(r, R"([{"type":"ptm"}])");
    EXPECT_EQ(st.attribute("PolyhedralTemplateMatching.counts.FCC"), double(r.atoms.size()));
    for (double v : st.props.at("Shear Strain")) EXPECT_LT(v, 1e-6);
  }
  // simple shear γ = 0.02 (x += γ y): the shear strain measure is γ/2 to first order
  {
    st = run(au, R"([{"type":"ptm"},{"type":"affine_transform","matrix":[1,0.02,0,0,1,0,0,0,1]}])");
    EXPECT_EQ(st.attribute("PolyhedralTemplateMatching.counts.FCC"), double(au.atoms.size()));
    EXPECT_NEAR(st.props.at("Shear Strain")[7], 0.01, 5e-4);
  }
  // a Mackay icosahedron of 13 atoms: the centre is icosahedral
  caps::System ico;
  const double g = (1 + std::sqrt(5.0)) / 2, r = 2.88 / std::sqrt(1 + g * g);
  std::vector<caps::Vec3> v{{0, 0, 0}};
  for (int s1 : {-1, 1})
    for (int s2 : {-1, 1}) v.push_back({0, s1 * r, s2 * g * r}), v.push_back({s1 * r, s2 * g * r, 0}), v.push_back({s2 * g * r, 0, s1 * r});
  for (const auto& x : v) { caps::Atom a; a.element = 79; a.pos = x; a.id = int64_t(ico.atoms.size() + 1); ico.atoms.push_back(a); }
  st = run(ico, R"([{"type":"ptm"}])");
  EXPECT_EQ(st.props.at("Structure Type")[0], 4.0);
  EXPECT_LT(st.props.at("RMSD")[0], 1e-6);
  // Combine datasets: a second file's particles appended
  const std::string f2 = (std::filesystem::temp_directory_path() / "caps_combine.xyz").string();
  caps::write_xyz(ico, f2);
  st = run(au, R"([{"type":"combine","path":")" + f2 + R"("}])");
  EXPECT_EQ(st.system.atoms.size(), au.atoms.size() + 13);
  EXPECT_EQ(st.attribute("CombineDatasets.added"), 13.0);
}

TEST(Recipe, CheckedWithoutRunning) {
  const auto ok = caps::check_recipe(caps::yaml_parse(
      "recipe: 1\nbuild: {polymer: {smiles: \"*CC(*)c1ccccc1\", dp: 40, chains: 20}}\ntype: {forcefield: gaff2}\ngrow: {density: 0.5}\n"
      "relax: {method: lbfgs, fmax: 0.02}\nequilibrate: {protocol: larsen21, t_max: 600}\nexport: [lammps, gromacs]\n"));
  EXPECT_EQ(ok.code, 0) << ok.error;
  ASSERT_EQ(ok.stages.size(), 6u);
  EXPECT_NE(ok.stages[0].summary.find("C8H8 · DP 40 × 20 chains"), std::string::npos);
  EXPECT_EQ(ok.protocol, "larsen21");
  EXPECT_EQ(ok.schedule.size(), 21u);
  const auto bad = caps::check_recipe(caps::yaml_parse("build: {polymer: {smiles: \"CC\"}}\nexport: [xyz, tiff]\n"));
  EXPECT_EQ(bad.code, 2);
  EXPECT_FALSE(bad.stages[0].ok);
  EXPECT_FALSE(bad.stages[1].ok);
  EXPECT_NE(bad.error.find("build"), std::string::npos);
}

#include "caps/molinfo.hpp"

TEST(MoleculeInfo, StyreneAndButane) {
  caps::BuildOptions b;
  b.forcefield = "uff";
  const auto sty = caps::build_molecule("C=Cc1ccccc1", b);
  const auto i = caps::molecule_info(sty.system, 0);
  EXPECT_EQ(i.formula, "C8H8");
  EXPECT_NEAR(i.mass, 104.152, 0.01);
  EXPECT_NEAR(i.monoisotopic, 104.0626, 1e-3);
  EXPECT_NEAR(i.dbe, 5.0, 1e-12);
  EXPECT_EQ(i.atoms, 16);
  EXPECT_EQ(i.rings, 1);
  EXPECT_EQ(i.rotatable, 1);
  EXPECT_LT(std::fabs(i.inertia_defect), 6.0);   // (near) planar
  EXPECT_LE(i.inertia[0], i.inertia[1]);
  const auto g = caps::parse_smiles(i.smiles);   // the written SMILES reads back as the same molecule
  caps::MolGraph gh = g;
  caps::add_hydrogens(gh);
  EXPECT_EQ(gh.atoms.size(), 16u) << i.smiles;
  const auto but = caps::molecule_info(caps::build_molecule("CCCC", b).system, 0);
  EXPECT_EQ(but.rotatable, 1);
  EXPECT_EQ(caps::molecule_info(caps::build_molecule("CCO", b).system, 0).rotatable, 0);
}

TEST(Sasa, IsolatedAtomAndTwoTouching) {
  caps::System s;
  caps::Atom a;
  a.element = 18;   // argon, Bondi 1.88 Å
  a.pos = {0, 0, 0};
  s.atoms.push_back(a);
  const auto one = caps::sasa(s, 1.4, 400);
  EXPECT_NEAR(one.total, 4 * M_PI * 3.28 * 3.28, 1e-9);
  a.pos = {3.0, 0, 0};
  s.atoms.push_back(a);
  const auto two = caps::sasa(s, 1.4, 400);
  EXPECT_LT(two.total, 2 * one.total);
  EXPECT_NEAR(two.area[0], two.area[1], 0.02 * two.area[0]);
  // the lens cut from each sphere: 2π R h with h = R − d/2
  const double R = 3.28, h = R - 1.5, exact = 2 * (4 * M_PI * R * R - 2 * M_PI * R * h);
  EXPECT_NEAR(two.total, exact, 0.02 * exact);
}
