#include <gtest/gtest.h>
#include <filesystem>

#include <cmath>
#include <string>
#include <vector>

#include "caps/molecule.hpp"

using namespace caps;

namespace {

const std::string kGaff = std::string(CAPS_SOURCE_DIR) + "/data/forcefields/gaff-amber25.json";

double dihedral(const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& d) {
  const Vec3 b1 = b - a, b2 = c - b, b3 = d - c;
  const Vec3 n1 = cross(b1, b2), n2 = cross(b2, b3);
  const Vec3 m = cross(n1, b2 * (1.0 / norm(b2)));
  return std::atan2(dot(m, n2), dot(n1, n2)) * 180 / 3.14159265358979323846;
}

int count(const MolGraph& g, int z) {
  int k = 0;
  for (const auto& a : g.atoms) k += a.element == z;
  return k;
}

}  // namespace

TEST(Smiles, ParsesTheGrammar) {
  MolGraph g = parse_smiles("CC(=O)Oc1ccccc1C(=O)O");   // aspirin
  EXPECT_EQ(g.heavy, 13);
  add_hydrogens(g);
  EXPECT_EQ(count(g, 1), 8);
  EXPECT_EQ(molecule_info(g).formula, "C9H8O4");
  EXPECT_NEAR(molecule_info(g).mass, 180.16, 0.01);
  EXPECT_EQ(molecule_info(g).rings, 1);

  g = parse_smiles("[NH4+].[Cl-]");
  EXPECT_EQ(g.parts, 2);
  add_hydrogens(g);
  EXPECT_EQ(molecule_info(g).formula, "ClH4N");   // Hill order without carbon: alphabetical
  EXPECT_EQ(molecule_info(g).charge, 0);

  g = parse_smiles("c1ccc2c(c1)[nH]c1ccccc12");   // carbazole: fused aromatic rings, [nH]
  add_hydrogens(g);
  EXPECT_EQ(molecule_info(g).formula, "C12H9N");
  EXPECT_EQ(molecule_info(g).rings, 3);

  g = parse_smiles("C%10CCCCC%10");
  add_hydrogens(g);
  EXPECT_EQ(molecule_info(g).formula, "C6H12");

  g = parse_smiles("c1ccsc1");   // thiophene: aromatic s takes no hydrogen
  add_hydrogens(g);
  EXPECT_EQ(molecule_info(g).formula, "C4H4S");

  g = parse_smiles("[2H]C([2H])([2H])O");
  EXPECT_EQ(g.atoms[0].isotope, 2);

  // a bond symbol before a ring-closure digit belongs to the ring bond, at the opening or the closing
  g = parse_smiles("C=1CCCCC1");
  add_hydrogens(g);
  EXPECT_EQ(molecule_info(g).formula, "C6H10");
  g = parse_smiles("[CH]1[C]2=[C]([NH][CH]=1)[CH]=[CH][CH]=[CH]2");   // indole, Kekulé, the ring bond double at its closure
  add_hydrogens(g);
  EXPECT_EQ(molecule_info(g).formula, "C8H7N");
}

TEST(Smiles, RejectsWhatIsNotSmiles) {
  for (const char* bad : {"C(C", "CC)", "C1CC", "C==C", "Q", "CH4", "[C", "C.", "C(=)C", "[Xy]"}) EXPECT_THROW(parse_smiles(bad), SmilesError) << bad;
  MolGraph g = parse_smiles("C(C)(C)(C)(C)C");   // five bonds on a carbon
  add_hydrogens(g);
  EXPECT_FALSE(valence_problems(g).empty());
}

TEST(Smiles, ChiralNeighbourOrder) {
  MolGraph g = parse_smiles("N[C@@H](C)C(=O)O");
  add_hydrogens(g);
  const MolAtom& c = g.atoms[1];
  ASSERT_EQ(c.order.size(), 4u);
  EXPECT_EQ(c.order[0], 0);                     // N, the atom before
  EXPECT_EQ(g.atoms[size_t(c.order[1])].element, 1);   // its own hydrogen comes next
  EXPECT_EQ(c.order[2], 2);                     // methyl
  EXPECT_EQ(c.order[3], 3);                     // carboxyl
}

TEST(Embed, LAlanineIsS) {
  // N[C@@H](C)C(=O)O is L-alanine, (S): with H pointing away, N → COOH → CH3 runs anticlockwise
  const BuildResult r = build_molecule("N[C@@H](C)C(=O)O", {});
  const auto& p = r.conformers.front().pos;
  const int H = r.graph.atoms[1].order[1];
  const Vec3 view = p[1] - p[size_t(H)];   // from the side opposite the hydrogen
  const Vec3 n = cross(p[3] - p[0], p[2] - p[0]);   // N → C(OOH) → CH3
  EXPECT_GT(dot(n, view), 0);
  for (int s : chirality_check(r.graph, p)) EXPECT_EQ(s, 1);

  const BuildResult d = build_molecule("N[C@H](C)C(=O)O", {});
  const auto& q = d.conformers.front().pos;
  const int Hd = d.graph.atoms[1].order[1];
  EXPECT_LT(dot(cross(q[3] - q[0], q[2] - q[0]), q[1] - q[size_t(Hd)]), 0);
}

TEST(Embed, DoubleBondGeometry) {
  // F–C=C–F dihedral; the atom indices of the four atoms in SMILES order
  auto torsion = [](const std::string& smi, int a, int b, int c, int d) {
    const BuildResult r = build_molecule(smi, {});
    const auto& p = r.conformers.front().pos;
    return std::fabs(dihedral(p[size_t(a)], p[size_t(b)], p[size_t(c)], p[size_t(d)]));
  };
  EXPECT_GT(torsion("F/C=C/F", 0, 1, 2, 3), 175);    // E
  EXPECT_LT(torsion("F/C=C\\F", 0, 1, 2, 3), 5);    // Z
  EXPECT_GT(torsion("C(\\F)=C/F", 1, 0, 2, 3), 175); // the same E, written from the other side
  EXPECT_LT(torsion("C(/F)=C/F", 1, 0, 2, 3), 5);    // Z
  EXPECT_LT(torsion("CC=CC", 0, 1, 2, 3) * (180 - torsion("CC=CC", 0, 1, 2, 3)), 5 * 180);   // unspecified: planar either way
}

TEST(Embed, AromaticRingIsFlatAndBondsAreSensible) {
  const BuildResult r = build_molecule("c1ccccc1", {});
  const auto& p = r.conformers.front().pos;
  for (const auto& b : r.graph.bonds) {
    const double d = norm(p[size_t(b.a)] - p[size_t(b.b)]);
    const bool ch = r.graph.atoms[size_t(b.a)].element == 1 || r.graph.atoms[size_t(b.b)].element == 1;
    EXPECT_NEAR(d, ch ? 1.07 : 1.40, 0.05);
  }
  EXPECT_LT(std::fabs(dihedral(p[0], p[1], p[2], p[3])), 5);
}

TEST(Embed, GaffCleanUpGivesAChairCyclohexane) {
  BuildOptions o;
  o.forcefield = kGaff;
  o.conformers = 4;
  const BuildResult r = build_molecule("C1CCCCC1", o);
  ASSERT_TRUE(r.conformers.front().minimised) << (r.notes.empty() ? "" : r.notes.front());
  const auto& p = r.conformers.front().pos;
  for (int k = 0; k < 6; ++k) {
    const double t = dihedral(p[size_t(k)], p[size_t((k + 1) % 6)], p[size_t((k + 2) % 6)], p[size_t((k + 3) % 6)]);
    EXPECT_NEAR(std::fabs(t), 55, 8) << "ring torsion " << k;
  }
  for (size_t k = 1; k < r.conformers.size(); ++k) EXPECT_GE(r.conformers[k].energy, r.conformers[0].energy);
}

TEST(Embed, DisconnectedPartsDoNotOverlap) {
  const BuildResult r = build_molecule("O.O.O", {});
  const auto& p = r.conformers.front().pos;
  EXPECT_GT(norm(p[0] - p[1]), 2.3);
  EXPECT_GT(norm(p[0] - p[2]), 2.3);
}

namespace {
// (S) at the alpha carbon of alanine built from `smi`, found by elements: with H away, N → COOH → CH3 anticlockwise
bool alanine_is_S(const std::string& smi) {
  const BuildResult r = build_molecule(smi, {});
  const auto& g = r.graph;
  const auto& p = r.conformers.front().pos;
  auto nbrs = [&](int i) {
    std::vector<int> v;
    for (const auto& b : g.bonds) {
      if (b.a == i) v.push_back(b.b);
      if (b.b == i) v.push_back(b.a);
    }
    return v;
  };
  int N = -1, CA = -1, CO = -1, CM = -1, H = -1;
  for (size_t i = 0; i < g.atoms.size(); ++i) if (g.atoms[i].element == 7) N = int(i);
  for (int v : nbrs(N)) if (g.atoms[size_t(v)].element == 6) CA = v;
  for (int v : nbrs(CA)) {
    const int z = g.atoms[size_t(v)].element;
    if (z == 1) H = v;
    if (z != 6) continue;
    int o = 0;
    for (int w : nbrs(v)) o += g.atoms[size_t(w)].element == 8;
    (o ? CO : CM) = v;
  }
  const Vec3 n = cross(p[size_t(CO)] - p[size_t(N)], p[size_t(CM)] - p[size_t(N)]);
  return dot(n, p[size_t(CA)] - p[size_t(H)]) > 0;
}
}  // namespace

TEST(Smiles, WriterRoundTrips) {
  for (const char* smi : {"CC(=O)Oc1ccccc1C(=O)O", "C1CC2CCC1C2", "c1ccc2cc3ccccc3cc2c1", "[NH4+].[Cl-]", "F/C=C/F", "F/C=C\\F",
                          "N[C@@H](C)C(=O)O", "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O", "C%10CCCCC%10", "[2H]C([2H])([2H])O", "N#Cc1ccccc1"}) {
    MolGraph g = parse_smiles(smi);
    const std::string w = write_smiles(g);
    MolGraph h = parse_smiles(w);
    EXPECT_EQ(write_smiles(h), w) << smi << " → " << w;   // idempotent
    add_hydrogens(g);
    add_hydrogens(h);
    EXPECT_EQ(molecule_info(g).formula, molecule_info(h).formula) << smi << " → " << w;
    EXPECT_EQ(molecule_info(g).stereocentres, molecule_info(h).stereocentres) << smi << " → " << w;
    EXPECT_EQ(molecule_info(g).stereo_bonds, molecule_info(h).stereo_bonds) << smi << " → " << w;
  }
}

TEST(Smiles, WriterKeepsTheConfiguration) {
  EXPECT_TRUE(alanine_is_S("N[C@@H](C)C(=O)O"));
  EXPECT_TRUE(alanine_is_S("C[C@@H](C(=O)O)N"));   // L-alanine as PubChem writes it
  EXPECT_FALSE(alanine_is_S("N[C@H](C)C(=O)O"));
  EXPECT_TRUE(alanine_is_S(write_smiles(parse_smiles("N[C@@H](C)C(=O)O"))));
  // start the writing from another atom: move the carboxyl oxygen to the front
  MolGraph g = parse_smiles("OC(=O)[C@@H](C)N");
  EXPECT_EQ(alanine_is_S("OC(=O)[C@@H](C)N"), alanine_is_S(write_smiles(g)));
}

TEST(Depict, RingsAreRegularAndChainsZigZag) {
  const MolGraph g = parse_smiles("c1ccc2ccccc2c1CCCC");
  const auto p = depict(g);
  ASSERT_EQ(p.size(), g.atoms.size());
  for (const auto& b : g.bonds) EXPECT_NEAR(norm(p[size_t(b.a)] - p[size_t(b.b)]), 1.0, 0.05);
  // naphthalene: every ring 1-3 distance √3
  EXPECT_NEAR(norm(p[0] - p[2]), std::sqrt(3.0), 0.08);
  // no two atoms on top of each other
  for (size_t i = 0; i < p.size(); ++i)
    for (size_t j = i + 1; j < p.size(); ++j) EXPECT_GT(norm(p[i] - p[j]), 0.8) << i << " " << j;
  for (const auto& v : p) EXPECT_EQ(v[2], 0.0);
}

#include "caps/bench.hpp"

TEST(Bench, QuickTablesPassAndExport) {
  BenchOptions o;
  o.samples = std::string(CAPS_SOURCE_DIR) + "/samples";
  o.forcefields = std::string(CAPS_SOURCE_DIR) + "/data/forcefields";
  o.quick = true;
  o.repeats = 1;
  std::vector<BenchTable> ts;
  for (const char* id : {"T1", "T9", "T12"}) {
    ts.push_back(run_bench(id, o));
    EXPECT_EQ(ts.back().status, "pass") << id << ": " << ts.back().note;
  }
  EXPECT_EQ(run_bench("T3", o).status, "not run");
  EXPECT_NE(bench_markdown(ts).find("| T12 Molecule builder |"), std::string::npos);
  EXPECT_NE(bench_latex(ts).find("\\begin{tabular}"), std::string::npos);
  EXPECT_NE(bench_csv(ts[0]).find("Check,Atoms"), std::string::npos);
}

// Stereo the SMILES leaves open is named in the notes; written stereo is not
TEST(Molecule, UnspecifiedStereoIsReported) {
  BuildOptions bo;
  bo.forcefield = "uff";
  auto has = [](const BuildResult& r, const char* what) {
    return std::any_of(r.notes.begin(), r.notes.end(), [&](const std::string& n) { return n.find(what) != std::string::npos; });
  };
  EXPECT_TRUE(has(build_molecule("CC(N)O", bo), "1 stereocentre not specified"));
  EXPECT_FALSE(has(build_molecule("C[C@H](N)O", bo), "not specified"));
  EXPECT_TRUE(has(build_molecule("CC=CC", bo), "1 double bond without"));
  EXPECT_FALSE(has(build_molecule("C/C=C/C", bo), "double bond without"));
  EXPECT_FALSE(has(build_molecule("CC(C)=CC", bo), "double bond without"));   // not stereogenic
}

// Bench T6 / T7: the user's cells named by reference id measured against data/reference/polymers.json; no cells, not run.
TEST(Bench, PropertiesAndChainsAgainstReference) {
  namespace fs = std::filesystem;
  const std::string root = CAPS_SOURCE_DIR;
  BenchOptions o;
  o.samples = root + "/samples";
  o.reference = root + "/data/reference/polymers.json";
  o.cells = (fs::temp_directory_path() / "caps_bench_cells_none").string();
  EXPECT_EQ(run_bench("T6", o).status, "not run");
  const fs::path cells = fs::temp_directory_path() / "caps_bench_cells";
  fs::create_directories(cells);
  fs::copy_file(root + "/samples/ps_melt.data", cells / "ps-atactic.data", fs::copy_options::overwrite_existing);
  o.cells = cells.string();
  const BenchTable t6 = run_bench("T6", o), t7 = run_bench("T7", o);
  ASSERT_GE(t6.rows.size(), 2u);
  EXPECT_EQ(t6.rows[0].cells[0], "Polystyrene, atactic");
  EXPECT_EQ(t6.rows[0].status, "fail");   // the sample is a small unequilibrated cell far below 1.04 g/cm³
  EXPECT_EQ(t6.rows[1].cells[1], "no cell");
  EXPECT_EQ(t7.rows[0].cells[4].rfind("9.5 – 10", 0), 0u);
  fs::remove_all(cells);
}

// Rotor search: every conformer of octane with it is at least as low as without; heavy atoms only when asked.
TEST(Molecule, RotorSearchAndHydrogensAsWritten) {
  BuildOptions o;
  o.forcefield = "uff";
  o.conformers = 4;
  o.seed = 3;
  const BuildResult plain = build_molecule("CCCCCCCC", o);
  o.rotor_search = true;
  const BuildResult rot = build_molecule("CCCCCCCC", o);
  EXPECT_LE(rot.conformers.front().energy, plain.conformers.front().energy + 1e-6);
  // all-trans: the end-to-end C1–C8 distance of the lowest conformer near the extended chain's (≈ 8.8 Å)
  const double ee = norm(rot.system.atoms[7].pos - rot.system.atoms[0].pos);
  EXPECT_GT(ee, 8.3);
  EXPECT_NE(rot.method.find("rotor search"), std::string::npos);
  BuildOptions ua;
  ua.implicit_hydrogens = false;
  const BuildResult h = build_molecule("CCCC", ua);
  EXPECT_EQ(h.system.atoms.size(), 4u);
}

// Analogs from R groups: a benzene core with two attachment points, three substituents on one and two on the other —
// six whole molecules, each built (one molecule, the right formula), hydrogen leaving the core atom's own H.
TEST(Molecule, AnalogsFromRGroups) {
  const auto an = enumerate_analogs("c1cc([*:1])ccc1[*:2]", {{1, {"H", "Cl", "C(=O)O"}}, {2, {"C", "*OC"}}});
  ASSERT_EQ(an.size(), 6u);
  std::set<std::string> formulas;
  for (const auto& a : an) {
    BuildOptions b;
    const auto r = build_molecule(a.smiles, b);
    int nm = 0;
    r.system.molecules(&nm);
    EXPECT_EQ(nm, 1) << a.smiles;
    formulas.insert(r.info.formula);
  }
  EXPECT_EQ(formulas.size(), 6u);
  EXPECT_TRUE(formulas.count("C7H8"));      // R1 = H, R2 = methyl: toluene
  EXPECT_TRUE(formulas.count("C8H8O3"));    // R1 = COOH, R2 = OMe: anisic acid
  EXPECT_EQ(an[0].name, "R1=H, R2=C");
  EXPECT_THROW(enumerate_analogs("c1ccccc1", {{1, {"C"}}}), std::invalid_argument);
}
