#include <gtest/gtest.h>

#include <string>
#include <vector>

#include "caps/elements.hpp"
#include "caps/ffdef.hpp"
#include "caps/typing.hpp"

using namespace caps;

namespace {

// A molecule from heavy atoms ("C", "N", ...), their hydrogen counts and heavy-atom bonds; hydrogens are appended
// after the heavy atoms. No bond orders: perception must find them.
System mol(const std::vector<std::string>& heavy, const std::vector<int>& nh, const std::vector<std::pair<int, int>>& bonds) {
  System s;
  for (const auto& e : heavy) {
    Atom a;
    a.element = element_from_symbol(e);
    s.atoms.push_back(a);
  }
  for (auto [i, j] : bonds) s.bonds.push_back({uint32_t(i), uint32_t(j), 0});
  for (size_t i = 0; i < heavy.size(); ++i)
    for (int k = 0; k < nh[i]; ++k) {
      Atom h;
      h.element = 1;
      s.atoms.push_back(h);
      s.bonds.push_back({uint32_t(i), uint32_t(s.atoms.size() - 1), 0});
    }
  return s;
}

System benzene() { return mol({"C", "C", "C", "C", "C", "C"}, {1, 1, 1, 1, 1, 1}, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 0}}); }

bool match(const std::string& sm, const System& s, uint32_t atom) {
  const Perception p = perceive(s);
  return Smarts(sm).matches(s, p, atom, {});
}

}  // namespace

TEST(Perception, AromaticRingsAndRingSizes) {
  const System b = benzene();
  const Perception p = perceive(b);
  for (int i = 0; i < 6; ++i) {
    EXPECT_TRUE(p.aromatic[i]);
    EXPECT_EQ(p.smallest_ring[i], 6);
    EXPECT_EQ(p.ring_count[i], 1);
  }
  EXPECT_FALSE(p.aromatic[6]);
  // cyclohexane: a ring, not aromatic
  const System ch = mol({"C", "C", "C", "C", "C", "C"}, {2, 2, 2, 2, 2, 2}, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 0}});
  const Perception q = perceive(ch);
  EXPECT_FALSE(q.aromatic[0]);
  EXPECT_EQ(q.smallest_ring[0], 6);
  // naphthalene: two rings, the fusion atoms in both, all aromatic
  const System n = mol({"C", "C", "C", "C", "C", "C", "C", "C", "C", "C"}, {1, 1, 1, 1, 0, 1, 1, 1, 1, 0},
                       {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 9}, {9, 0}, {4, 5}, {5, 6}, {6, 7}, {7, 8}, {8, 9}});
  const Perception r = perceive(n);
  EXPECT_EQ(r.rings.size(), 2u);
  for (int i = 0; i < 10; ++i) EXPECT_TRUE(r.aromatic[i]) << i;
  EXPECT_EQ(r.ring_count[4], 2);
}

TEST(Perception, HeteroaromaticsAndLonePairs) {
  // pyrrole (N-H gives two π electrons), furan, pyridine
  const System pyrrole = mol({"N", "C", "C", "C", "C"}, {1, 1, 1, 1, 1}, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 0}});
  EXPECT_TRUE(perceive(pyrrole).aromatic[0]);
  const System furan = mol({"O", "C", "C", "C", "C"}, {0, 1, 1, 1, 1}, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 0}});
  EXPECT_TRUE(perceive(furan).aromatic[0]);
  const System pyridine = mol({"N", "C", "C", "C", "C", "C"}, {0, 1, 1, 1, 1, 1}, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 0}});
  EXPECT_TRUE(perceive(pyridine).aromatic[0]);
  // cyclopentadiene: the CH2 breaks the ring
  const System cp = mol({"C", "C", "C", "C", "C"}, {2, 1, 1, 1, 1}, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 0}});
  EXPECT_FALSE(perceive(cp).aromatic[1]);
}

TEST(Perception, BondOrdersAndCharges) {
  // acetonitrile: C#N
  const System acn = mol({"C", "C", "N"}, {3, 0, 0}, {{0, 1}, {1, 2}});
  EXPECT_EQ(perceive(acn).bond_order(1, 2), 3);
  // acetate: one C=O, one C-O(−)
  const System ac = mol({"C", "C", "O", "O"}, {3, 0, 0, 0}, {{0, 1}, {1, 2}, {1, 3}});
  const Perception pa = perceive(ac);
  EXPECT_EQ(pa.bond_order(1, 2) + pa.bond_order(1, 3), 3);
  EXPECT_EQ(pa.charge[2] + pa.charge[3], -1);
  // nitromethane: N+(=O)O−
  const System nm = mol({"C", "N", "O", "O"}, {3, 0, 0, 0}, {{0, 1}, {1, 2}, {1, 3}});
  const Perception pn = perceive(nm);
  EXPECT_EQ(pn.charge[1], 1);
  EXPECT_EQ(pn.charge[2] + pn.charge[3], -1);
  // methylammonium: four-connected N+
  const System ma = mol({"C", "N"}, {3, 3}, {{0, 1}});
  EXPECT_EQ(perceive(ma).charge[1], 1);
  // imidazolium (both N protonated): C=N+ and an aromatic ring
  const System im = mol({"N", "C", "N", "C", "C"}, {1, 1, 1, 1, 1}, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 0}});
  const Perception pi = perceive(im);
  EXPECT_EQ(pi.charge[0] + pi.charge[2], 1);
  for (int i = 0; i < 5; ++i) EXPECT_TRUE(pi.aromatic[i]) << i;
  // bond orders given in the file are kept (mol2): propene written with its double bond
  System pr = mol({"C", "C", "C"}, {2, 1, 3}, {{0, 1}, {1, 2}});
  pr.bonds[0].order = 2;
  for (size_t k = 1; k < pr.bonds.size(); ++k) pr.bonds[k].order = 1;
  EXPECT_EQ(perceive(pr).bond_order(0, 1), 2);
}

TEST(Smarts, Primitives) {
  const System b = benzene();
  EXPECT_TRUE(match("c", b, 0));
  EXPECT_FALSE(match("C", b, 0));
  EXPECT_TRUE(match("[cH1]", b, 0));
  EXPECT_TRUE(match("[#6;X3;R1;r6]", b, 0));
  EXPECT_TRUE(match("c1ccccc1", b, 0));
  EXPECT_FALSE(match("c1cccc1", b, 0));
  EXPECT_TRUE(match("[H]c", b, 6));
  EXPECT_TRUE(match("[H][c;$(c1ccccc1)]", b, 6));
  EXPECT_TRUE(match("[!N;!O]", b, 0));
  EXPECT_TRUE(match("[N,c]", b, 0));
  EXPECT_TRUE(match("c:c", b, 0));
  EXPECT_FALSE(match("c=c", b, 0));   // aromatic bonds are ':' only
  EXPECT_TRUE(match("c~c", b, 0));
  EXPECT_TRUE(match("c@c", b, 0));
  // ethanol: C(H3)-C(H2)-O-H
  const System et = mol({"C", "C", "O"}, {3, 2, 1}, {{0, 1}, {1, 2}});
  EXPECT_TRUE(match("[CX4H3][CX4H2][OX2H1]", et, 0));
  EXPECT_TRUE(match("[OX2H1]", et, 2));
  EXPECT_TRUE(match("[CH2](C)O", et, 1));
  EXPECT_TRUE(match("[C;D4;v4;+0]", et, 0));
  EXPECT_FALSE(match("[C;R]", et, 0));
  EXPECT_TRUE(match("[C;R0]", et, 0));
  // acetone carbonyl
  const System ac = mol({"C", "C", "O", "C"}, {3, 0, 0, 3}, {{0, 1}, {1, 2}, {1, 3}});
  EXPECT_TRUE(match("[CX3](=O)([#6])[#6]", ac, 1));
  EXPECT_TRUE(match("O=C", ac, 2));
  EXPECT_FALSE(match("O-C", ac, 2));
}

TEST(Smarts, Errors) {
  for (const char* bad : {"", "[C", "C(", "C)", "C1CC", "[Q]", "C.C", "%x"}) EXPECT_THROW(Smarts{bad}, FFError) << bad;
}

TEST(Typing, PrecedenceOverridesAndTypeReferences) {
  FFDef ff;
  ff.name = "test";
  for (const char* t : {"c", "c3", "hc", "ho", "oh", "cx"}) {
    FFType x;
    x.name = t;
    ff.types.push_back(x);
  }
  ff.typing = {{"c", "[CX4]", "", {}, 0},
               {"c3", "[CX4H3]", "", {}, 1},
               {"oh", "[OX2H1]", "", {}, 0},
               {"hc", "[H]C", "", {}, 0},
               {"ho", "[H]O", "", {}, 0},
               // a type defined through another type (foyer's %type): C next to a hydroxyl O
               {"cx", "[C;$(C[%oh])]", "", {"c"}, 0}};
  const System et = mol({"C", "C", "O"}, {3, 2, 1}, {{0, 1}, {1, 2}});
  const TypingResult r = assign_types(et, ff);
  EXPECT_EQ(r.types[0], "c3");   // priority 1 over c
  EXPECT_EQ(r.types[1], "cx");   // overrides c, found once oh was assigned
  EXPECT_EQ(r.types[2], "oh");
  for (size_t i = 3; i < et.atoms.size() - 1; ++i) EXPECT_EQ(r.types[i], "hc");
  EXPECT_EQ(r.types.back(), "ho");
  EXPECT_EQ(r.untyped, 0);
  // an atom no rule matches is reported
  System w = et;
  Atom na;
  na.element = 11;
  w.atoms.push_back(na);
  EXPECT_EQ(assign_types(w, ff).untyped, 1);
}

TEST(Typing, ShippedRulesTypeCommonMolecules) {
  for (const char* lib : {"pcff", "cvff"}) {
    FFDef ff = load_forcefield(std::string(CAPS_SOURCE_DIR) + "/data/forcefields/" + lib + "-dlfield.json");
    load_typing(ff, std::string(CAPS_SOURCE_DIR) + "/data/typing/" + lib + "-dlfield.typing.json");
    const TypingResult r = assign_types(benzene(), ff);
    EXPECT_EQ(r.types[0], "cp") << lib;
    EXPECT_EQ(r.untyped, 0) << lib;
    const System ac = mol({"C", "C", "O", "O"}, {3, 0, 0, 1}, {{0, 1}, {1, 2}, {1, 3}});   // acetic acid
    const TypingResult q = assign_types(ac, ff);
    EXPECT_EQ(q.untyped, 0) << lib;
    EXPECT_EQ(q.types[1], std::string(lib) == "pcff" ? "c_1" : "c'");
    EXPECT_EQ(q.types[3], std::string(lib) == "pcff" ? "o_2" : "oh");
  }
}

namespace {

TypingResult type_with(const char* ff_name, const System& s) {
  FFDef ff = load_forcefield(std::string(CAPS_SOURCE_DIR) + "/data/forcefields/" + ff_name + ".json");
  if (ff.typing.empty()) load_typing(ff, std::string(CAPS_SOURCE_DIR) + "/data/typing/" + ff_name + ".typing.json");
  return assign_types(s, ff);
}

System butadiene() { return mol({"C", "C", "C", "C"}, {2, 1, 1, 2}, {{0, 1}, {1, 2}, {2, 3}}); }
System furan() { return mol({"O", "C", "C", "C", "C"}, {0, 1, 1, 1, 1}, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 0}}); }

}  // namespace

TEST(Perception, CationsFoundBySearch) {
  // benzamidinium: the ring stays aromatic and one amidine N carries the charge
  const System b = mol({"C", "C", "C", "C", "C", "C", "C", "N", "N"}, {0, 1, 1, 1, 1, 1, 0, 2, 2},
                       {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 0}, {0, 6}, {6, 7}, {6, 8}});
  const Perception p = perceive(b);
  for (int i = 0; i < 6; ++i) EXPECT_TRUE(p.aromatic[i]) << i;
  EXPECT_EQ(p.charge[7] + p.charge[8], 1);
  // methyl azide: R-N=N+=N- or R-N(-)-N+#N
  const System az = mol({"C", "N", "N", "N"}, {3, 0, 0, 0}, {{0, 1}, {1, 2}, {2, 3}});
  const Perception pa = perceive(az);
  EXPECT_EQ(pa.charge[2], 1);
  EXPECT_EQ(pa.charge[1] + pa.charge[3], -1);
}

TEST(Typing, GaffFollowsAntechamber) {
  // antechamber: furan os / cc cd cd cc / h4 on the carbons next to O, ha on the others; butadiene c2 ce ce c2
  const TypingResult f = type_with("gaff-amber16-dlfield", furan());
  EXPECT_EQ(f.types[0], "os");
  EXPECT_EQ(f.types[1], "cc");
  EXPECT_EQ(f.types[2], "cd");
  EXPECT_EQ(f.types[3], "cd");
  EXPECT_EQ(f.types[4], "cc");
  EXPECT_EQ(f.types[5], "h4");
  EXPECT_EQ(f.types[6], "ha");
  const TypingResult b = type_with("gaff-amber16-dlfield", butadiene());
  EXPECT_EQ(b.types[0], "c2");
  EXPECT_EQ(b.types[1], "ce");
  EXPECT_EQ(b.types[2], "ce");
  EXPECT_EQ(b.types[3], "c2");
  // GAFF2 amine N by hydrogen count: methylamine n8
  const System ma = mol({"C", "N"}, {3, 2}, {{0, 1}});
  EXPECT_EQ(type_with("gaff-amber25-dlfield", ma).types[1], "n8");
  EXPECT_EQ(type_with("gaff-amber16-dlfield", ma).types[1], "n3");
}

TEST(Typing, CgenffAndOpls) {
  // CGenFF: conjugated diene CG2DC3 CG2DC2 CG2DC1 CG2DC3 (same type across a double bond, switched across the single)
  const TypingResult b = type_with("cgenff-dlfield", butadiene());
  EXPECT_EQ(b.types[0], "CG2DC3");
  EXPECT_EQ(b.types[1], "CG2DC2");
  EXPECT_EQ(b.types[2], "CG2DC1");
  EXPECT_EQ(b.types[3], "CG2DC3");
  const TypingResult f = type_with("cgenff-dlfield", furan());
  EXPECT_EQ(f.types[0], "OG2R50");
  EXPECT_EQ(f.types[1], "CG2R51");
  // OPLS-AA (DL_FIELD names): methyl acetate CT CO4 O OES CT
  const System ma = mol({"C", "C", "O", "O", "C"}, {3, 0, 0, 0, 3}, {{0, 1}, {1, 2}, {1, 3}, {3, 4}});
  const TypingResult o = type_with("opls2005-dlfield", ma);
  EXPECT_EQ(o.types[0], "CT");
  EXPECT_EQ(o.types[1], "CO4");
  EXPECT_EQ(o.types[2], "O");
  EXPECT_EQ(o.types[3], "OES");
  EXPECT_EQ(o.types[4], "CT");
}
