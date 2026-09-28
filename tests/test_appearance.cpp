#include <gtest/gtest.h>

#include <cmath>

#include "caps/appearance.hpp"
#include "caps/crystal.hpp"
#include "caps/molecule.hpp"
#include "caps/peptide.hpp"
#include "caps/render.hpp"
#include "caps/polymer.hpp"

using namespace caps;

namespace {

const std::string kCrystals = std::string(CAPS_SOURCE_DIR) + "/data/crystals/";

System molecule(const std::string& smiles) {
  BuildOptions o;
  o.forcefield = "uff";
  return build_molecule(smiles, o).system;
}

std::string label_of_first_centre(const System& s) {
  const auto l = stereo_labels(s);
  for (const auto& x : l) if (!x.empty()) return x;
  return "";
}

}  // namespace

TEST(Appearance, CipLabelsOfAminoAcids) {
  EXPECT_EQ(label_of_first_centre(molecule("N[C@@H](C)C(=O)O")), "S");    // L-alanine
  EXPECT_EQ(label_of_first_centre(molecule("N[C@H](C)C(=O)O")), "R");     // D-alanine
  EXPECT_EQ(label_of_first_centre(molecule("N[C@@H](CS)C(=O)O")), "R");   // L-cysteine: S outranks the carboxyl
  EXPECT_EQ(label_of_first_centre(molecule("C[C@@H](O)CC")), "R");        // (R)-butan-2-ol (PubChem CC[C@@H](C)O)
  EXPECT_EQ(label_of_first_centre(molecule("CC(C)CC")), "");              // no centre
  // L-threonine is (2S,3R)
  const System thr = molecule("C[C@H]([C@@H](C(=O)O)N)O");
  const auto l = stereo_labels(thr);
  int r = 0, s = 0;
  for (const auto& x : l) r += x == "R", s += x == "S";
  EXPECT_EQ(r, 1);
  EXPECT_EQ(s, 1);
}

TEST(Appearance, AccessibleSurfaceOfOneAtom) {
  System s;
  Atom a;
  a.element = 18;   // argon, Bondi 1.88 Å
  s.atoms.push_back(a);
  SurfaceOptions o;
  o.spacing = 0.3;
  const Mesh m = surface_mesh(s, o);
  const double r = 1.88 + 1.4;
  EXPECT_NEAR(m.area(), 4 * M_PI * r * r, 0.04 * 4 * M_PI * r * r);
  for (size_t v = 0; v < m.vertices.size(); v += 97) EXPECT_NEAR(norm(m.vertices[v]), r, 0.1);
  o.kind = SurfaceKind::VanDerWaals;
  EXPECT_NEAR(surface_mesh(s, o).area(), 4 * M_PI * 1.88 * 1.88, 0.05 * 4 * M_PI * 1.88 * 1.88);
}

TEST(Appearance, ExcludedSurfaceLiesBetween) {
  PeptideOptions po;
  po.sequence = "GAG";
  po.cleanup = false;
  const System pep = build_peptide(po);
  SurfaceOptions o;
  o.spacing = 0.4;
  const double sas = surface_mesh(pep, o).area();
  o.kind = SurfaceKind::VanDerWaals;
  const double vdw = surface_mesh(pep, o).area();
  o.kind = SurfaceKind::Excluded;
  const Mesh ses = surface_mesh(pep, o);
  EXPECT_LT(ses.area(), sas);
  EXPECT_LT(ses.area(), vdw);   // the excluded surface smooths the crevices of the van der Waals one
  EXPECT_GT(ses.area(), 0.3 * vdw);
  // the potential of a zwitterion: positive near NH3+, negative near COO−
  const auto phi = surface_potential(pep, ses);
  double lo = 1e9, hi = -1e9;
  for (double x : phi) lo = std::min(lo, x), hi = std::max(hi, x);
  EXPECT_LT(lo, -10);
  EXPECT_GT(hi, 10);
}

TEST(Appearance, TetrahedraOfQuartz) {
  const System q = read_cif(kCrystals + "alpha-quartz.cif");
  std::vector<char> si(q.atoms.size(), 0);
  int nsi = 0;
  for (size_t i = 0; i < q.atoms.size(); ++i) if (q.atoms[i].element == 14) si[i] = 1, ++nsi;
  const Mesh m = polyhedra(q, si);
  EXPECT_EQ(m.triangles.size(), size_t(4 * nsi));   // one SiO4 tetrahedron each
}

TEST(Appearance, RendersMixedStylesAndMeshes) {
  PeptideOptions po;
  po.sequence = "AEAAAKA";
  po.structure = "HHHHHHH";
  po.cleanup = false;
  const System pep = build_peptide(po);
  RenderOptions ro;
  ro.width = 160, ro.height = 120;
  ro.atom_style.assign(pep.atoms.size(), uint8_t(Style::Wireframe));
  for (size_t i = 0; i < pep.atoms.size() / 2; ++i) ro.atom_style[i] = uint8_t(Style::SpaceFilling);
  SurfaceOptions so;
  const Mesh m = surface_mesh(pep, so);
  ro.meshes.push_back({&m, 0x8FB8D8, 0.5f});
  ro.colour_by = ColourBy::Property;
  ro.ramp = Ramp::BlueOrange;
  ro.symmetric = true;
  for (const auto& a : pep.atoms) ro.property.push_back(a.charge);
  Renderer r;
  const Image img = r.render(pep, Camera{}, ro);
  ASSERT_EQ(img.rgba.size(), size_t(160 * 120 * 4));
  size_t lit = 0;
  for (size_t k = 0; k < img.rgba.size(); k += 4) lit += img.rgba[k] + img.rgba[k + 1] + img.rgba[k + 2] > 90;
  EXPECT_GT(lit, 1000u);
  EXPECT_FALSE(ribbon_paths(pep, {}).empty());
}

// E/Z of double bonds from the 3D geometry and CIP ranks: trans- and cis-2-butene, and every backbone double bond of a
// grown cis-1,4-polyisoprene chain Z (natural rubber)
TEST(Appearance, EzLabelsOfDoubleBonds) {
  BuildOptions bo;
  bo.forcefield = "uff";
  const System e = build_molecule("C/C=C/C", bo).system;
  const System z = build_molecule("C/C=C\\C", bo).system;
  std::vector<std::pair<uint32_t, uint32_t>> be, bz;
  const auto le = ez_labels(e, &be), lz = ez_labels(z, &bz);
  ASSERT_EQ(be.size(), 1u);
  ASSERT_EQ(bz.size(), 1u);
  EXPECT_EQ(le[be[0].first], "E");
  EXPECT_EQ(lz[bz[0].first], "Z");
  // isobutylene has two methyls on one end: no E/Z
  EXPECT_TRUE(ez_labels(build_molecule("CC(C)=CC", bo).system).at(1).empty());
  ChainSpec spec;
  spec.units = {RepeatUnit{"isoprene", "[*]C/C=C(C)\\C[*]"}};
  spec.dp = 4;
  GrowOptions g;
  g.chains = 1;
  g.dp = 0;
  g.density = 0.05;
  g.seed = 1;
  const System nr = grow_chains(spec, g);
  std::vector<std::pair<uint32_t, uint32_t>> bn;
  const auto ln = ez_labels(nr, &bn);
  ASSERT_EQ(bn.size(), 3u);   // the tail unit's double bond carries two methyls at its capped end: not stereogenic
  for (const auto& [a, b] : bn) EXPECT_EQ(ln[a], "Z");
}

TEST(Appearance, IsosurfaceOfASphereAndAPeriodicSlab) {
  // a field falling off from the cell centre: the isosurface at 6 Å is a sphere of area 4πr²
  Cell c;
  c.a = {30, 0, 0}, c.b = {0, 30, 0}, c.c = {0, 0, 30};
  const int n[3] = {60, 60, 60};
  std::vector<double> v(size_t(n[0]) * n[1] * n[2]);
  for (int i = 0; i < n[0]; ++i)
    for (int j = 0; j < n[1]; ++j)
      for (int k = 0; k < n[2]; ++k) {
        const Vec3 p{30.0 * (i + 0.5) / n[0] - 15, 30.0 * (j + 0.5) / n[1] - 15, 30.0 * (k + 0.5) / n[2] - 15};
        v[(size_t(i) * n[1] + j) * n[2] + k] = 10 - norm(p);
      }
  const Mesh m = isosurface(c, n, v, 4.0);
  EXPECT_NEAR(m.area(), 4 * M_PI * 36, 0.02 * 4 * M_PI * 36);
  for (size_t q = 0; q < m.vertices.size(); q += 97) EXPECT_NEAR(norm(m.vertices[q] - Vec3{15, 15, 15}), 6.0, 0.05);
  // outward normals point down the field (away from the centre)
  EXPECT_GT(dot(m.normals[0], m.vertices[0] - Vec3{15, 15, 15}), 0);
  // a slab across the periodic cell: two faces, each the cell's cross-section, closed through the faces
  for (int i = 0; i < n[0]; ++i)
    for (int j = 0; j < n[1]; ++j)
      for (int k = 0; k < n[2]; ++k) v[(size_t(i) * n[1] + j) * n[2] + k] = std::fabs(30.0 * (k + 0.5) / n[2] - 15) < 5 ? 1.0 : 0.0;
  const Mesh slab = isosurface(c, n, v, 0.5);
  EXPECT_NEAR(slab.area(), 2 * 30 * 30, 1.0);
}
