#include <gtest/gtest.h>

#include <cmath>

#include "caps/appearance.hpp"
#include "caps/edit.hpp"
#include "caps/elements.hpp"
#include "caps/json.hpp"
#include <fstream>
#include <sstream>
#include "caps/molecule.hpp"
#include "caps/peptide.hpp"
#include "caps/polymer.hpp"
#include "caps/torsion.hpp"
#include "caps/typing.hpp"

using namespace caps;

namespace {

System molecule(const std::string& smiles) {
  BuildOptions o;
  o.forcefield = "uff";
  return build_molecule(smiles, o).system;
}

double angle_deg(const Vec3& a, const Vec3& b, const Vec3& c) {
  const Vec3 u = a - b, v = c - b;
  return std::acos(dot(u, v) / (norm(u) * norm(v))) * 180 / M_PI;
}

System polystyrene(Tacticity t) {
  ChainSpec spec;
  spec.units.push_back({"styrene", "[*]CC([*])c1ccccc1"});
  spec.dp = 12;
  spec.tacticity = t;
  GrowOptions o;
  o.chains = 1;
  o.density = 0.3;
  o.seed = 3;
  return grow_chains(spec, o);
}

}  // namespace

TEST(Edit, HydrogensFillValencesTetrahedrally) {
  System s;
  for (int z : {6, 6, 8}) { Atom a; a.element = z; s.atoms.push_back(a); }
  s.atoms[1].pos = {1.52, 0, 0};
  s.atoms[2].pos = {2.0, 1.35, 0};
  s.bonds = {{0, 1, 1}, {1, 2, 1}};
  EXPECT_EQ(add_hydrogens(s), 6);   // CH3 CH2 OH
  ASSERT_EQ(s.atoms.size(), 9u);
  int per[3] = {0, 0, 0};
  for (const auto& b : s.bonds)
    for (uint32_t k : {b.i, b.j}) if (k < 3 && s.atoms[b.i + b.j - k].element == 1) ++per[k];
  EXPECT_EQ(per[0], 3);
  EXPECT_EQ(per[1], 2);
  EXPECT_EQ(per[2], 1);
  // every angle at the methyl carbon near tetrahedral
  std::vector<uint32_t> around;
  for (const auto& b : s.bonds) { if (b.i == 0) around.push_back(b.j); if (b.j == 0) around.push_back(b.i); }
  for (size_t a = 0; a < around.size(); ++a)
    for (size_t b = a + 1; b < around.size(); ++b)
      EXPECT_NEAR(angle_deg(s.atoms[around[a]].pos, s.atoms[0].pos, s.atoms[around[b]].pos), 109.47, 3.0);
  EXPECT_EQ(add_hydrogens(s), 0);   // nothing left to fill
  EXPECT_EQ(default_valence(7, 1), 4);
}

TEST(Edit, AtomsBondsElementsAndDeletion) {
  System s = molecule("C");
  const size_t n0 = s.atoms.size();   // CH4
  // replace one H by a carbon: ethane's skeleton, then fill it
  uint32_t h = 1;
  set_element(s, h, 6);
  EXPECT_EQ(s.atoms[h].element, 6);
  EXPECT_EQ(add_hydrogens(s), 3);
  EXPECT_EQ(s.atoms.size(), n0 + 3);
  const uint32_t o = add_atom(s, int(h), 8, 1);
  EXPECT_NEAR(norm(s.atoms[o].pos - s.atoms[h].pos), element(8).covalent + element(6).covalent, 1e-9);
  add_bond(s, 0, o, 1);
  EXPECT_TRUE(remove_bond(s, 0, o));
  std::vector<char> del(s.atoms.size(), 0);
  del[o] = 1;
  const size_t nb = s.bonds.size();
  delete_atoms(s, del);
  EXPECT_EQ(s.atoms.size(), n0 + 3);
  EXPECT_EQ(s.bonds.size(), nb - 1);
}

TEST(Edit, InvertingTheCentreOfAlanine) {
  System ala = molecule("N[C@@H](C)C(=O)O");   // L = S
  uint32_t ca = 1;
  auto lab = stereo_labels(ala);
  ASSERT_EQ(lab[ca], "S");
  invert_centre(ala, ca);
  lab = stereo_labels(ala);
  EXPECT_EQ(lab[ca], "R");
}

TEST(Edit, TacticityOfGrownPolystyrene) {
  EXPECT_EQ(tacticity(polystyrene(Tacticity::Isotactic)).label, "isotactic");
  EXPECT_EQ(tacticity(polystyrene(Tacticity::Syndiotactic)).label, "syndiotactic");
  System at = polystyrene(Tacticity::Atactic);
  const auto T = tacticity(at);
  ASSERT_EQ(T.chains.size(), 1u);
  EXPECT_GE(T.chains[0].centres.size(), 10u);
  EXPECT_EQ(T.m + T.r, int(T.chains[0].centres.size()) - 1);
  set_tacticity(at, true);
  EXPECT_EQ(tacticity(at).label, "isotactic");
  set_tacticity(at, false);
  EXPECT_EQ(tacticity(at).label, "syndiotactic");
}

TEST(Edit, Selections) {
  const System tol = molecule("Cc1ccccc1");
  const auto ring = select_smarts(tol, "c1ccccc1");
  EXPECT_EQ(std::count(ring.begin(), ring.end(), 1), 6);
  const auto c = select_element(tol, "C");
  EXPECT_EQ(std::count(c.begin(), c.end(), 1), 7);
  auto grown = select_grow(tol, ring, 1);
  EXPECT_EQ(std::count(grown.begin(), grown.end(), 1), 6 + 5 + 1);   // the ring, its 5 H and the methyl C
  const auto near = select_within(tol, ring, 1.2);
  EXPECT_EQ(std::count(near.begin(), near.end(), 1), 6 + 5);          // ring H within 1.2 Å, the methyl C (1.5 Å) not
  EXPECT_THROW(select_element(tol, "Xx"), EditError);
}

TEST(Edit, AttachingFragments) {
  System s = molecule("C");   // methane
  const auto added = attach_fragment(s, 0, "[CH3]*");
  EXPECT_EQ(s.atoms.size(), 8u);   // ethane
  EXPECT_EQ(added.size(), 4u);     // C and three H (the target's H replaced)
  int cc = 0;
  for (const auto& b : s.bonds)
    if (s.atoms[b.i].element == 6 && s.atoms[b.j].element == 6) { ++cc; EXPECT_NEAR(norm(s.atoms[b.i].pos - s.atoms[b.j].pos), 1.52, 0.05); }
  EXPECT_EQ(cc, 1);
  System bz = molecule("c1ccccc1");
  attach_fragment(bz, 0, "*C(=O)O*");   // benzoic acid: the second point becomes the acid H
  EXPECT_EQ(bz.atoms.size(), 15u);
  int o = 0;
  for (const auto& a : bz.atoms) o += a.element == 8;
  EXPECT_EQ(o, 2);
  // nothing too close outside bonds
  for (size_t i = 0; i < bz.atoms.size(); ++i)
    for (size_t j = i + 1; j < bz.atoms.size(); ++j) {
      bool bonded = false;
      for (const auto& b : bz.bonds) bonded |= (b.i == i && b.j == j) || (b.i == j && b.j == i);
      if (!bonded) EXPECT_GT(norm(bz.atoms[i].pos - bz.atoms[j].pos), 0.9) << i << " " << j;
    }
  EXPECT_EQ(fragment_attach_atoms("*C(=O)O*").size(), 2u);
  EXPECT_THROW(attach_fragment(bz, 0, "CCO"), EditError);   // no attachment point
}

TEST(Edit, EveryLibraryFragmentBuilds) {
  // data/fragments/catalogue.json: each fragment embeds (with its * points as hydrogens) and, when it has * points,
  // reports the atoms they hang on
  std::ifstream in(std::string(CAPS_SOURCE_DIR) + "/data/fragments/catalogue.json");
  ASSERT_TRUE(in.good());
  std::stringstream ss;
  ss << in.rdbuf();
  const Json cat = Json::parse(ss.str());
  int n = 0, failed = 0;
  for (const auto& f : cat["fragments"].items()) {
    const std::string smi = f.text("smiles");
    try {
      if (smi.find('*') != std::string::npos) EXPECT_FALSE(fragment_attach_atoms(smi).empty()) << smi;
      else {
        BuildOptions o;
        o.forcefield = "uff";
        EXPECT_GT(build_molecule(smi, o).system.atoms.size(), 0u) << smi;
      }
    } catch (const std::exception& e) {
      ADD_FAILURE() << f.text("name") << " " << smi << ": " << e.what();
      ++failed;
    }
    ++n;
  }
  EXPECT_GE(n, 100);
  EXPECT_EQ(failed, 0);
}

// Fusing a benzene ring onto a C–C bond of benzene gives naphthalene: C10H8, eleven aromatic ring bonds, nothing missing
TEST(Edit, FuseBenzeneMakesNaphthalene) {
  BuildOptions b;
  b.forcefield = "uff";
  System s = build_molecule("c1ccccc1", b).system;
  uint32_t i = 0, j = 0;
  for (const auto& bd : s.bonds)
    if (s.atoms[bd.i].element == 6 && s.atoms[bd.j].element == 6) { i = bd.i, j = bd.j; break; }
  const auto added = fuse_benzene(s, i, j);
  EXPECT_EQ(added.size(), 8u);
  int c = 0, h = 0, arom = 0;
  for (const auto& a : s.atoms) c += a.element == 6, h += a.element == 1;
  for (const auto& bd : s.bonds) arom += bd.order == 4;
  EXPECT_EQ(c, 10);
  EXPECT_EQ(h, 8);
  EXPECT_EQ(arom, 11);
  EXPECT_EQ(add_hydrogens(s), 0);
  // the new ring is a hexagon: every new C–C 1.39 ± 0.05 Å from its neighbours
  for (const auto& bd : s.bonds)
    if (s.atoms[bd.i].element == 6 && s.atoms[bd.j].element == 6) EXPECT_NEAR(norm(s.atoms[bd.j].pos - s.atoms[bd.i].pos), 1.40, 0.06);
  EXPECT_THROW(fuse_benzene(s, 0, 17), EditError);
}

// Protonation by pH gives back the peptide builder's own hydrogens: a peptide built at pH 7 and at pH 2, stripped to
// its heavy atoms and protonated again at the same pH, has the same number of hydrogens and the same net charge
TEST(Edit, ProtonationByPhMatchesThePeptideBuilder) {
  for (double ph : {7.0, 2.0, 11.0}) {
    PeptideOptions p;
    p.sequence = "KDEHRYCA";
    p.ph = ph;
    p.cleanup = false;
    // the builder's terminal forms are set, not titrated: those of the model pKa values (8.0, 3.1) at this pH
    p.n_term = ph < 8.0 ? "NH3+" : "NH2";
    p.c_term = ph > 3.1 ? "COO-" : "COOH";
    System s = build_peptide(p);
    int h0 = 0;
    double q0 = 0;
    for (const auto& a : s.atoms) h0 += a.element == 1, q0 += a.charge;
    // heavy atoms only, charges gone (as a PDB of heavy atoms reads)
    std::vector<char> hyd(s.atoms.size(), 0);
    for (size_t i = 0; i < s.atoms.size(); ++i) hyd[i] = s.atoms[i].element == 1;
    delete_atoms(s, hyd);
    for (auto& a : s.atoms) a.charge = 0;
    std::vector<std::string> notes;
    add_hydrogens_at_ph(s, ph, {}, &notes);
    int h1 = 0, q1 = 0;
    for (const auto& a : s.atoms) h1 += a.element == 1, q1 += int(std::lround(a.charge));
    EXPECT_EQ(h1, h0) << "pH " << ph << " · " << notes.front();
    EXPECT_NEAR(q1, std::lround(q0), 0.01) << "pH " << ph << " · " << notes.front();
  }
}

// Exact geometry edits: butane's C1–C2 bond, C1–C2–C3 angle and dihedral set by moving one side, the rest untouched;
// a rigid rotation keeps every distance; a mirror image inverts the stereocentre; R/S set by name; a ring bond refused
TEST(Edit, ExactGeometryRotateMirrorAndConfiguration) {
  BuildOptions bo;
  bo.forcefield = "uff";
  System s = build_molecule("CCCC", bo).system;
  auto dist = [&](uint32_t a, uint32_t b) { return norm(s.atoms[a].pos - s.atoms[b].pos); };
  auto angle = [&](uint32_t a, uint32_t b, uint32_t c) {
    const Vec3 u = s.atoms[a].pos - s.atoms[b].pos, v = s.atoms[c].pos - s.atoms[b].pos;
    return std::acos(dot(u, v) / (norm(u) * norm(v))) * 180 / M_PI;
  };
  const double c34 = dist(2, 3);
  set_bond_length(s, 1, 0, 1.60);   // C2–C1: C1's side (with its hydrogens) moves
  EXPECT_NEAR(dist(0, 1), 1.60, 1e-9);
  EXPECT_NEAR(dist(2, 3), c34, 1e-9);
  set_bond_angle(s, 0, 1, 2, 120.0);
  EXPECT_NEAR(angle(0, 1, 2), 120.0, 1e-7);
  set_torsion(s, 0, 1, 2, 3, 60.0);
  EXPECT_NEAR(dihedral_angle(s, {0, 1, 2, 3}), 60.0, 1e-6);
  EXPECT_NEAR(dist(0, 1), 1.60, 1e-9);
  // rotation of the whole molecule about an arbitrary axis: every distance kept
  std::vector<uint32_t> all(s.atoms.size());
  for (uint32_t k = 0; k < all.size(); ++k) all[k] = k;
  const double d03 = dist(0, 3);
  rotate_atoms(s, all, {1, 2, 3}, 37.0);
  EXPECT_NEAR(dist(0, 3), d03, 1e-9);
  // a stereocentre: mirrored it turns over; set_configuration puts it back
  System c = build_molecule("C[C@H](N)O", bo).system;
  std::vector<char> only(c.atoms.size(), 0);
  only[1] = 1;
  const std::string before = stereo_labels(c, only)[1];
  ASSERT_FALSE(before.empty());
  std::vector<uint32_t> every(c.atoms.size());
  for (uint32_t k = 0; k < every.size(); ++k) every[k] = k;
  mirror_atoms(c, every, {0, 0, 1});
  const std::string mirrored = stereo_labels(c, only)[1];
  EXPECT_NE(mirrored, before);
  EXPECT_TRUE(set_configuration(c, 1, before));
  EXPECT_EQ(stereo_labels(c, only)[1], before);
  EXPECT_FALSE(set_configuration(c, 0, "R"));   // the methyl carbon is no stereocentre
  // a bond in a ring cannot be stretched by moving one side
  System ring = build_molecule("C1CCCCC1", bo).system;
  EXPECT_THROW(set_bond_length(ring, 0, 1, 1.7), EditError);
}

// Coordination geometries: an octahedral Fe made exact from a distorted start (cis 90°, trans 180°, bond lengths kept, a
// hydroxo ligand's hydrogen carried along); square planar from a tetrahedral start; a dative bond in no atom's valence.
TEST(Edit, CoordinationGeometryAndDativeBonds) {
  System s;
  auto add = [&](int z, Vec3 p) { Atom a; a.element = z; a.pos = p; s.atoms.push_back(a); return uint32_t(s.atoms.size() - 1); };
  const uint32_t fe = add(26, {0, 0, 0});
  const Vec3 start[6] = {{2.3, 0.3, 0.1}, {-2.2, 0.2, -0.4}, {0.4, 2.1, 0.3}, {-0.2, -2.4, 0.5}, {0.3, -0.2, 2.2}, {0.5, 0.4, -2.0}};
  std::vector<uint32_t> lig;
  for (const auto& p : start) { lig.push_back(add(17, p)); s.bonds.push_back({fe, lig.back(), 1}); }
  // one ligand a hydroxo: O on the Fe, its H bonded to it
  s.atoms[lig[0]].element = 8;
  const uint32_t h = add(1, s.atoms[lig[0]].pos + Vec3{0.6, 0.7, 0});
  s.bonds.push_back({lig[0], h, 1});
  std::vector<double> len;
  for (uint32_t l : lig) len.push_back(norm(s.atoms[l].pos - s.atoms[fe].pos));
  const double oh = norm(s.atoms[h].pos - s.atoms[lig[0]].pos);
  set_coordination(s, fe, "octahedral");
  for (size_t a = 0; a < lig.size(); ++a) {
    EXPECT_NEAR(norm(s.atoms[lig[a]].pos - s.atoms[fe].pos), len[a], 1e-9);
    for (size_t b = a + 1; b < lig.size(); ++b) {
      const Vec3 u = s.atoms[lig[a]].pos - s.atoms[fe].pos, v = s.atoms[lig[b]].pos - s.atoms[fe].pos;
      const double ang = std::acos(std::clamp(dot(u, v) / (norm(u) * norm(v)), -1.0, 1.0)) * 180 / M_PI;
      EXPECT_TRUE(std::fabs(ang - 90) < 1e-6 || std::fabs(ang - 180) < 1e-6) << a << " " << b << " " << ang;
    }
  }
  EXPECT_NEAR(norm(s.atoms[h].pos - s.atoms[lig[0]].pos), oh, 1e-9);   // the hydroxo's H came with its O
  EXPECT_THROW(set_coordination(s, fe, "tetrahedral"), EditError);   // six neighbours, four sites
  // square planar from a tetrahedral start
  System p;
  p.atoms.push_back(s.atoms[fe]);
  const double r3 = 1 / std::sqrt(3.0);
  for (Vec3 d : {Vec3{r3, r3, r3}, Vec3{r3, -r3, -r3}, Vec3{-r3, r3, -r3}, Vec3{-r3, -r3, r3}}) {
    Atom a; a.element = 17; a.pos = d * 2.3; p.atoms.push_back(a);
    p.bonds.push_back({0, uint32_t(p.atoms.size() - 1), 1});
  }
  set_coordination(p, 0, "square_planar");
  double sum = 0;
  for (uint32_t a = 1; a <= 4; ++a)
    for (uint32_t b = a + 1; b <= 4; ++b) {
      const Vec3 u = p.atoms[a].pos - p.atoms[0].pos, v = p.atoms[b].pos - p.atoms[0].pos;
      sum += std::acos(std::clamp(dot(u, v) / (norm(u) * norm(v)), -1.0, 1.0)) * 180 / M_PI;
    }
  EXPECT_NEAR(sum, 4 * 90 + 2 * 180, 1e-6);
  // a water on a metal: dative, no formal charge on O; a plain single bond would make it +1
  System w;
  auto addw = [&](int z, Vec3 q) { Atom a; a.element = z; a.pos = q; w.atoms.push_back(a); };
  addw(29, {0, 0, 0}), addw(8, {2.0, 0, 0}), addw(1, {2.6, 0.75, 0}), addw(1, {2.6, -0.75, 0});
  w.bonds = {{0, 1, kBondDative}, {1, 2, 1}, {1, 3, 1}};
  EXPECT_EQ(perceive(w).charge[1], 0);
  EXPECT_EQ(add_hydrogens(w), 0);
  w.bonds[0].order = 1;
  EXPECT_EQ(perceive(w).charge[1], 1);
}
