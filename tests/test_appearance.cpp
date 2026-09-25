#include <gtest/gtest.h>

#include <cmath>

#include "caps/appearance.hpp"
#include "caps/crystal.hpp"
#include "caps/molecule.hpp"
#include "caps/peptide.hpp"
#include "caps/render.hpp"

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
