#include <gtest/gtest.h>

#include <cmath>
#include <map>
#include <random>

#include "caps/elements.hpp"
#include "caps/field.hpp"
#include "caps/molecule.hpp"
#include "caps/uff.hpp"

using namespace caps;

namespace {

System from_smiles(const std::string& smi, bool embedded = true) {
  MolGraph g = parse_smiles(smi);
  add_hydrogens(g);
  return molecule_system(g, embedded ? embed(g) : std::vector<Vec3>(g.atoms.size()));
}

std::map<std::string, int> label_counts(const std::string& smi) {
  std::map<std::string, int> c;
  for (const auto& t : uff_types(from_smiles(smi, false))) c[t]++;
  return c;
}

double bond_length(const System& s, uint32_t i, uint32_t j) { return norm(s.atoms[i].pos - s.atoms[j].pos); }

}  // namespace

TEST(Uff, TableCoversEveryElementToLawrencium) {
  EXPECT_EQ(uff_label_count(), 127);
  EXPECT_TRUE(uff_has_label("C_R"));
  EXPECT_TRUE(uff_has_label("Si3"));
  EXPECT_TRUE(is_uff("uff"));
  EXPECT_TRUE(is_uff("UFF"));
  EXPECT_TRUE(is_uff("/some/where/uff.json"));
  EXPECT_FALSE(is_uff("gaff2-moltemplate.json"));
  // every element as a dichloride takes a label of its own element
  for (int z = 3; z <= 103; ++z) {
    System s;
    Atom m;
    m.element = z;
    s.atoms.push_back(m);
    for (int k = 0; k < 2; ++k) {
      Atom cl;
      cl.element = 17;
      cl.pos = {k ? 2.3 : -2.3, 0.1 * k, 0};
      s.atoms.push_back(cl);
      s.bonds.push_back({0, uint32_t(k + 1), 1});
    }
    const ForceField ff = assign_uff(s);
    const std::string sym = element(z).symbol;
    EXPECT_EQ(ff.atom_type[0].rfind(sym.size() == 1 ? sym + "_" : sym, 0), 0u) << sym << " typed " << ff.atom_type[0];
    EXPECT_EQ(ff.bonds.size(), 2u);
  }
}

TEST(Uff, TypesFollowHybridisationAndConjugation) {
  auto c = label_counts("CC");
  EXPECT_EQ(c["C_3"], 2);
  EXPECT_EQ(c["H_"], 6);
  EXPECT_EQ(label_counts("c1ccccc1")["C_R"], 6);
  EXPECT_EQ(label_counts("C=C")["C_2"], 2);
  EXPECT_EQ(label_counts("CC#N")["C_1"], 1);
  EXPECT_EQ(label_counts("CC#N")["N_1"], 1);
  c = label_counts("CC(=O)C");   // acetone: an isolated carbonyl
  EXPECT_EQ(c["C_2"], 1);
  EXPECT_EQ(c["O_2"], 1);
  c = label_counts("CC(=O)OC");   // methyl acetate: the ester is conjugated
  EXPECT_EQ(c["C_R"], 1);
  EXPECT_EQ(c["O_R"], 2);
  EXPECT_EQ(label_counts("CC(=O)NC")["N_R"], 1);   // amide N
  EXPECT_EQ(label_counts("CN(C)C")["N_3"], 1);
  EXPECT_EQ(label_counts("c1ccncc1")["N_R"], 1);
  EXPECT_EQ(label_counts("CS(=O)C")["S_3+4"], 1);
  EXPECT_EQ(label_counts("CS(=O)(=O)C")["S_3+6"], 1);
  EXPECT_EQ(label_counts("CSC")["S_3+2"], 1);
  EXPECT_EQ(label_counts("c1ccsc1")["S_R"], 1);
  EXPECT_EQ(label_counts("COP(=O)(OC)OC")["P_3+5"], 1);
  EXPECT_EQ(label_counts("CP(C)C")["P_3+3"], 1);
  c = label_counts("C[Si](C)(C)O[Si](C)(C)C");   // a siloxane (silicone rubber)
  EXPECT_EQ(c["Si3"], 2);
  EXPECT_EQ(c["O_3"], 1);
  EXPECT_EQ(label_counts("ClC(Cl)Cl")["Cl"], 3);
  EXPECT_EQ(label_counts("CC(=C)C=C")["C_R"], 4);   // isoprene: conjugated diene
  EXPECT_EQ(label_counts("CC(C)=CC")["C_2"], 2);    // a polyisoprene unit: isolated double bond
  EXPECT_EQ(label_counts("F[S](F)(F)(F)(F)F")["S_3+6"], 1);
}

TEST(Uff, ForcesMatchFiniteDifferences) {
  // every UFF term form: sp3 and sp2 / sp / trigonal angles, a linear centre, sp2 and P inversions, an octahedral
  // centre (SF6), the torsion rules, and full 1-4 van der Waals
  System s = from_smiles("CC#CC(=O)Oc1ccc(cc1)P(C)C.F[S](F)(F)(F)(F)F");
  const ForceField ff = assign_uff(s);
  std::map<int, int> forms;
  for (const auto& a : ff.angles_x) forms[a.form]++;
  EXPECT_GT(forms[3], 0);
  EXPECT_GT(forms[11], 0);
  EXPECT_GT(forms[13], 0);
  EXPECT_GT(forms[14], 0);
  int inv1 = 0, inv2 = 0;
  for (const auto& v : ff.inversions) (v.form == 1 ? inv1 : inv2)++;
  EXPECT_GT(inv1, 0);
  EXPECT_EQ(inv2, 1);
  EXPECT_TRUE(ff.pairs14.empty());

  std::mt19937 rng(3);
  std::normal_distribution<double> nd(0, 0.08);
  std::vector<double> x;
  for (const auto& a : s.atoms) for (int k = 0; k < 3; ++k) x.push_back(a.pos[k] + nd(rng));
  EnergyOptions eo;
  eo.coulomb = false;
  eo.tail = false;
  eo.cutoff = 30;
  Evaluator ev(ff, eo);
  std::vector<double> f;
  const Cell none;
  ev.compute(x, none, f);
  const double h = 1e-5;
  double worst = 0;
  for (size_t i = 0; i < x.size(); ++i) {
    auto xp = x, xm = x;
    xp[i] += h;
    xm[i] -= h;
    std::vector<double> ft;
    const double ep = ev.compute(xp, none, ft).total(), em = ev.compute(xm, none, ft).total();
    const double fd = -(ep - em) / (2 * h);
    worst = std::max(worst, std::fabs(fd - f[i]) / std::max(1.0, std::fabs(fd)));
  }
  EXPECT_LT(worst, 1e-4);
}

TEST(Uff, CleanUpGivesTextbookGeometry) {
  BuildOptions o;
  o.forcefield = "uff";
  // butane, a siloxane, an ester and dimethyl sulfone: every conformer minimises and bond lengths are UFF's
  for (const std::string smi : {"CCCC", "C[Si](C)(C)O[Si](C)(C)C", "CC(=O)OC", "CS(=O)(=O)C", "c1ccccc1Br"}) {
    const BuildResult r = build_molecule(smi, o);
    ASSERT_TRUE(r.conformers.front().minimised) << smi;
    EXPECT_NE(r.method.find("UFF"), std::string::npos);
    for (const auto& b : r.system.bonds) {
      const int zi = r.system.atoms[b.i].element, zj = r.system.atoms[b.j].element;
      const double d = bond_length(r.system, b.i, b.j);
      if (zi == 6 && zj == 6 && b.order == 1) EXPECT_NEAR(d, 1.51, 0.04) << smi;
      if ((zi == 14 && zj == 8) || (zi == 8 && zj == 14)) EXPECT_NEAR(d, 1.73, 0.03) << smi;   // UFF's Si–O (long, as published)
      if ((zi == 1 || zj == 1) && (zi == 6 || zj == 6)) EXPECT_NEAR(d, 1.10, 0.04) << smi;
      if ((zi == 35 || zj == 35)) EXPECT_NEAR(d, 1.90, 0.08) << smi;
    }
  }
  // a planar benzene ring after clean-up
  const BuildResult bz = build_molecule("c1ccccc1", o);
  const auto& p = bz.system.atoms;
  const Vec3 n = cross(p[1].pos - p[0].pos, p[2].pos - p[0].pos);
  for (int k = 3; k < 6; ++k) EXPECT_NEAR(dot(p[k].pos - p[0].pos, n) / norm(n), 0.0, 0.02);
}

TEST(Uff, CleansUpMetalComplexesAndMainGroupHypervalence) {
  BuildOptions o;
  o.forcefield = "uff";
  // cisplatin (square planar Pt), a tetrahedral zinc complex, SF6 and a phosphine: all minimise
  for (const std::string smi : {"N[Pt](N)(Cl)Cl", "Cl[Zn](Cl)(Cl)Cl", "F[S](F)(F)(F)(F)F", "CP(C)C"}) {
    BuildResult r;
    ASSERT_NO_THROW(r = build_molecule(smi, o)) << smi;
    EXPECT_TRUE(r.conformers.front().minimised) << smi;
  }
}
