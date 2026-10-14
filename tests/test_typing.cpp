#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
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

// @ / @@ from the positions: CHFClBr with F above the carbon and Cl, Br, H anticlockwise below it, seen from F
TEST(Typing, SmartsChiralityFromPositions) {
  System s = mol({"C", "F", "Cl", "Br"}, {1, 0, 0, 0}, {{0, 1}, {0, 2}, {0, 3}});
  const double h = -1.0 / 3, c = std::cos(2 * M_PI / 3), sn = std::sin(2 * M_PI / 3);
  s.atoms[1].pos = {0, 0, 1};
  s.atoms[2].pos = {1, 0, h};
  s.atoms[3].pos = {c, sn, h};
  s.atoms[4].pos = {c, -sn, h};   // H
  const Perception p = perceive(s);
  EXPECT_TRUE(Smarts("[C@](F)(Cl)(Br)[H]").matches(s, p, 0, {}));
  EXPECT_FALSE(Smarts("[C@@](F)(Cl)(Br)[H]").matches(s, p, 0, {}));
  EXPECT_TRUE(Smarts("[C@@](F)(Br)(Cl)[H]").matches(s, p, 0, {}));   // two neighbours swapped: the other sense
  EXPECT_TRUE(Smarts("[C@@H](F)(Cl)Br").matches(s, p, 0, {}));        // the bracket's H first: an odd permutation
  EXPECT_TRUE(Smarts("[C@](F)(Cl)Br").matches(s, p, 0, {}));          // the unnamed fourth neighbour last
  EXPECT_TRUE(Smarts("F[C@](Cl)(Br)[H]").matches(s, p, 1, {}));       // reached from F: F is still first
  EXPECT_FALSE(Smarts("F[C@@](Cl)(Br)[H]").matches(s, p, 1, {}));
  EXPECT_TRUE(Smarts("[C](F)(Cl)(Br)[H]").matches(s, p, 0, {}));
  for (auto& a : s.atoms) a.pos = {0, 0, 0};   // no positions: no handedness, a chiral pattern does not match
  EXPECT_FALSE(Smarts("[C@](F)(Cl)(Br)[H]").matches(s, perceive(s), 0, {}));
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
  for (const char* lib : {"pcff-frc", "cvff-frc"}) {
    const FFDef ff = load_forcefield(std::string(CAPS_SOURCE_DIR) + "/data/forcefields/" + lib + ".json");   // its own typing
    const TypingResult r = assign_types(benzene(), ff);
    EXPECT_EQ(r.types[0], "cp") << lib;
    EXPECT_EQ(r.untyped, 0) << lib;
    const System ac = mol({"C", "C", "O", "O"}, {3, 0, 0, 1}, {{0, 1}, {1, 2}, {1, 3}});   // acetic acid
    const TypingResult q = assign_types(ac, ff);
    EXPECT_EQ(q.untyped, 0) << lib;
    EXPECT_EQ(q.types[1], std::string(lib) == "pcff-frc" ? "c_1" : "c'");
    EXPECT_EQ(q.types[3], std::string(lib) == "pcff-frc" ? "o_2" : "oh");
  }
}

// PCFF's c=1 (pcff.frc: "next to end doubly bonded carbon") as DL_FIELD's typing gives it: a C=C carbon single-bonded to
// another non-aromatic C=C carbon; an aromatic ring or a C=O next to it does not make it c=1
TEST(Typing, PcffConjugatedAlkeneCarbonsAsDlField) {
  const FFDef ff = load_forcefield(std::string(CAPS_SOURCE_DIR) + "/data/forcefields/pcff-frc.json");
  const System bd = mol({"C", "C", "C", "C"}, {2, 1, 1, 2}, {{0, 1}, {1, 2}, {2, 3}});   // buta-1,3-diene
  const TypingResult b = assign_types(bd, ff);
  EXPECT_EQ(b.types[0], "c=");
  EXPECT_EQ(b.types[1], "c=1");
  EXPECT_EQ(b.types[2], "c=1");
  const System ht = mol({"C", "C", "C", "C", "C", "C"}, {2, 1, 1, 1, 1, 2}, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}});   // hexa-1,3,5-triene
  const TypingResult h = assign_types(ht, ff);
  EXPECT_EQ(h.types[2], "c=1");
  EXPECT_EQ(h.types[3], "c=1");
  const System st = mol({"C", "C", "C", "C", "C", "C", "C", "C"}, {2, 1, 0, 1, 1, 1, 1, 1},
                        {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 6}, {6, 7}, {7, 2}});   // styrene
  const TypingResult y = assign_types(st, ff);
  EXPECT_EQ(y.types[0], "c=");
  EXPECT_EQ(y.types[1], "c=");
  EXPECT_EQ(y.types[2], "cp");
  const System ac = mol({"C", "C", "C", "O"}, {2, 1, 1, 0}, {{0, 1}, {1, 2}, {2, 3}});   // acrolein
  EXPECT_EQ(assign_types(ac, ff).types[1], "c=");
}

// the automatic bends of cvff.frc / pcff.frc (#quadratic_angle *_auto) by their numbered wildcards: the entry of the end
// atom with the lowest number (*2 H before *7 / *9 C), the generic * X * last. CVFF: propanal's C-C(=O)-H takes
// *2 c'_ h_ (110, 55); PCFF: styrene's ring C-CH=-H takes *2 c=_ h_ (120, 37.5), not * c=_ * (120, 60)
TEST(Typing, FrcAutomaticBendsByNumberedWildcard) {
  const System pr = mol({"C", "C", "C", "O"}, {3, 2, 1, 0}, {{0, 1}, {1, 2}, {2, 3}});
  const System st = mol({"C", "C", "C", "C", "C", "C", "C", "C"}, {2, 1, 0, 1, 1, 1, 1, 1},
                        {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 6}, {6, 7}, {7, 2}});
  struct Case { const char* lib; System s; uint32_t centre, heavy_end; double theta0, k; };
  for (const Case& c : {Case{"cvff-frc", pr, 2, 1, 110.0, 55.0}, Case{"pcff-frc", st, 1, 2, 120.0, 37.5}}) {
    const FFDef def = load_forcefield(std::string(CAPS_SOURCE_DIR) + "/data/forcefields/" + c.lib + ".json");
    System s = c.s;
    const TypingResult tr = assign_types(s, def);
    ASSERT_EQ(tr.untyped, 0) << c.lib;
    ParamReport rep;
    const ForceField f = parameterize(s, def, tr.types, "types", &rep, false);
    bool found = false;
    const auto check = [&](uint32_t i, uint32_t j, uint32_t k, double theta0, double kt) {
      if (j != c.centre || !((i == c.heavy_end && s.atoms[k].element == 1) || (k == c.heavy_end && s.atoms[i].element == 1))) return;
      found = true;
      EXPECT_NEAR(theta0 * 180 / 3.14159265358979323846, c.theta0, 1e-6) << c.lib;
      EXPECT_NEAR(kt, c.k, 1e-6) << c.lib;
    };
    for (const auto& a : f.angles) check(a.i, a.j, a.k, a.theta0, a.kt);
    for (const auto& a : f.angles2) check(a.i, a.j, a.k, a.theta0, a.k2);
    EXPECT_TRUE(found) << c.lib;
  }
}

namespace {

TypingResult type_with(const char* ff_name, const System& s) {
  FFDef ff = load_forcefield(std::string(CAPS_SOURCE_DIR) + "/data/forcefields/" + ff_name + ".json");
  if (ff.typing.empty()) load_typing(ff, std::string(CAPS_SOURCE_DIR) + "/data/typing/" + ff_name + ".typing.json");
  return assign_types(s, ff);
}

// a library built on this machine from files CAPS does not distribute (CHARMM's toppar: bench/ff/open/build_charmm.sh)
bool built_here(const char* ff_name) {
  return std::filesystem::exists(std::string(CAPS_SOURCE_DIR) + "/data/forcefields/" + ff_name + ".json");
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
  const TypingResult f = type_with("gaff-1.81-ambertools", furan());
  EXPECT_EQ(f.types[0], "os");
  EXPECT_EQ(f.types[1], "cc");
  EXPECT_EQ(f.types[2], "cd");
  EXPECT_EQ(f.types[3], "cd");
  EXPECT_EQ(f.types[4], "cc");
  EXPECT_EQ(f.types[5], "h4");
  EXPECT_EQ(f.types[6], "ha");
  const TypingResult b = type_with("gaff-1.81-ambertools", butadiene());
  EXPECT_EQ(b.types[0], "c2");
  EXPECT_EQ(b.types[1], "ce");
  EXPECT_EQ(b.types[2], "ce");
  EXPECT_EQ(b.types[3], "c2");
  // GAFF2 amine N by hydrogen count: methylamine n8
  const System ma = mol({"C", "N"}, {3, 2}, {{0, 1}});
  EXPECT_EQ(type_with("gaff2-2.2.20-ambertools", ma).types[1], "n8");
  EXPECT_EQ(type_with("gaff-1.81-ambertools", ma).types[1], "n3");
}

TEST(Typing, GaffRingsAndPolarHydrogensAsAntechamber) {
  // checked against antechamber 22.0 (AmberTools 23.6), bench/typing/compare_antechamber.py:
  // naphthalene: all ten carbons ca (a double bond into the fused ring is not exocyclic)
  const System nap = mol({"C", "C", "C", "C", "C", "C", "C", "C", "C", "C"}, {1, 1, 1, 0, 1, 1, 1, 1, 0, 1},
                         {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 6}, {6, 7}, {7, 8}, {8, 3}, {8, 9}, {9, 0}});
  for (const char* ff : {"gaff-1.81-ambertools", "gaff2-2.2.20-ambertools"}) {
    const TypingResult t = type_with(ff, nap);
    for (int i = 0; i < 10; ++i) EXPECT_EQ(t.types[i], "ca") << ff << " atom " << i;
  }
  // cyclooctatetraene: a planar conjugated ring, cc / cd (not the chain types ce / cf)
  const System cot = mol({"C", "C", "C", "C", "C", "C", "C", "C"}, {1, 1, 1, 1, 1, 1, 1, 1},
                         {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 6}, {6, 7}, {7, 0}});
  for (const std::string& t : type_with("gaff2-2.2.20-ambertools", cot).types)
    if (t[0] == 'c') EXPECT_TRUE(t == "cc" || t == "cd") << t;
  // acetaldehyde: the formyl H is h4 (the O of C=O counts as electron-withdrawing); DMF's formyl H is h5 (N and O)
  const System ald = mol({"C", "C", "O"}, {3, 1, 0}, {{0, 1}, {1, 2}});
  EXPECT_EQ(type_with("gaff2-2.2.20-ambertools", ald).types[6], "h4");
  const System dmf = mol({"C", "N", "C", "C", "O"}, {1, 0, 3, 3, 0}, {{0, 1}, {1, 2}, {1, 3}, {0, 4}});
  EXPECT_EQ(type_with("gaff2-2.2.20-ambertools", dmf).types[5], "h5");
  // guanidinium: all three N alike, nv (GAFF2) / nh (GAFF), the formal =NH2+ included
  const System gua = mol({"C", "N", "N", "N"}, {0, 2, 2, 2}, {{0, 1}, {0, 2}, {0, 3}});
  for (int i = 1; i <= 3; ++i) {
    EXPECT_EQ(type_with("gaff2-2.2.20-ambertools", gua).types[i], "nv") << i;
    EXPECT_EQ(type_with("gaff-1.81-ambertools", gua).types[i], "nh") << i;
  }
}

TEST(Typing, OplsAa2024AsTheGromacsTopologies) {
  // checked against the GROMACS OPLS-AA topologies of Caleman et al. 2012 (bench/typing/compare_oplsaa_top.py);
  // the type name starts with the OPLS number
  const auto num = [](const std::string& t) { return t.substr(0, t.find('_')); };
  // dichloromethane: Cl 151, C 152, H 153 (formerly left untyped: the table has no CH2Cl2 numbers)
  const TypingResult dcm = type_with("oplsaa2024-moltemplate", mol({"C", "Cl", "Cl"}, {2, 0, 0}, {{0, 1}, {0, 2}}));
  EXPECT_EQ(dcm.untyped, 0);
  EXPECT_EQ(num(dcm.types[0]), "152");
  EXPECT_EQ(num(dcm.types[1]), "151");
  EXPECT_EQ(num(dcm.types[3]), "153");
  // diethyl sulfide: CH2 on S 210, CH3 135 (the sulfide C pattern asked the S for three neighbours)
  const TypingResult des = type_with("oplsaa2024-moltemplate", mol({"C", "C", "S", "C", "C"}, {3, 2, 0, 2, 3}, {{0, 1}, {1, 2}, {2, 3}, {3, 4}}));
  EXPECT_EQ(num(des.types[1]), "210");
  EXPECT_EQ(num(des.types[2]), "202");
  EXPECT_EQ(num(des.types[0]), "135");
  // pyrimidine: N 530, C2 531, C4 532, C5 533
  const TypingResult pym = type_with("oplsaa2024-moltemplate", mol({"N", "C", "N", "C", "C", "C"}, {0, 1, 0, 1, 1, 1},
                                                                   {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 0}}));
  EXPECT_EQ(num(pym.types[0]), "530");
  EXPECT_EQ(num(pym.types[1]), "531");
  EXPECT_EQ(num(pym.types[3]), "532");
  EXPECT_EQ(num(pym.types[4]), "533");
  // glycerol: the triol types O 171, CH2 173, CH 174
  const TypingResult gly = type_with("oplsaa2024-moltemplate", mol({"O", "C", "C", "O", "C", "O"}, {1, 2, 1, 1, 2, 1},
                                                                   {{0, 1}, {1, 2}, {2, 3}, {2, 4}, {4, 5}}));
  EXPECT_EQ(num(gly.types[0]), "171");
  EXPECT_EQ(num(gly.types[1]), "173");
  EXPECT_EQ(num(gly.types[2]), "174");
}

TEST(Typing, MoleculeClassConditions) {
  // a rule keyed by its molecule (in_molecule / not_in_molecule, as the CHARMM family keys an sp3 CH3): a sulfoxide's
  // methyls C_L3, an ether's CAE3, any other C3
  const auto dir = std::filesystem::temp_directory_path();
  const std::string rules = (dir / "caps_molecule_class.typing.json").string(), ffp = (dir / "caps_molecule_class.json").string();
  std::ofstream(rules) << R"({"format": "caps-typing", "version": 1, "rules": [
 {"type": "C3", "smarts": "[CX4H3]", "priority": 1},
 {"type": "C_L3", "smarts": "[CX4H3]", "priority": 2, "in_molecule": ["[SX3]=[OX1]"]},
 {"type": "CAE3", "smarts": "[CX4H3]", "priority": 2, "in_molecule": ["[OX2]([CX4])[CX4]"], "not_in_molecule": ["[#7]"]},
 {"type": "X", "smarts": "[!#6]", "priority": 0}, {"type": "X", "smarts": "[#6;!X4]", "priority": 0}, {"type": "X", "smarts": "[CX4;!H3]", "priority": 0}]})";
  std::ofstream(ffp) << R"({"format": "caps-forcefield", "format_version": 1, "name": "molecule classes", "typing": "caps_molecule_class.typing.json",
 "atom_types": [{"name": "C3", "element": "C"}, {"name": "C_L3", "element": "C"}, {"name": "CAE3", "element": "C"}, {"name": "X", "element": "X"}]})";
  const FFDef ff = load_forcefield(ffp);
  const System dmso = mol({"C", "S", "O", "C"}, {3, 0, 0, 3}, {{0, 1}, {1, 2}, {1, 3}});
  EXPECT_EQ(assign_types(dmso, ff).types[0], "C_L3");
  const System dme = mol({"C", "O", "C"}, {3, 0, 3}, {{0, 1}, {1, 2}});
  EXPECT_EQ(assign_types(dme, ff).types[0], "CAE3");
  // an amino ether: not_in_molecule keeps it out of CAE3
  const System ame = mol({"C", "O", "C", "N"}, {3, 0, 2, 2}, {{0, 1}, {1, 2}, {2, 3}});
  EXPECT_EQ(assign_types(ame, ff).types[0], "C3");
  // both together in one structure: each molecule keeps its own
  System two = dmso;
  const uint32_t off = uint32_t(two.atoms.size());
  for (const auto& a : dme.atoms) two.atoms.push_back(a);
  for (const auto& b : dme.bonds) two.bonds.push_back({b.i + off, b.j + off, b.order});
  const TypingResult t2 = assign_types(two, ff);
  EXPECT_EQ(t2.types[0], "C_L3");
  EXPECT_EQ(t2.types[off], "CAE3");
  std::filesystem::remove(rules);
  std::filesystem::remove(ffp);
}

TEST(Typing, CgenffConjugationPhaseAsTheRtf) {
  if (!built_here("cgenff-toppar")) GTEST_SKIP() << "CGenFF is built locally from CHARMM's toppar (bench/ff/open/build_charmm.sh)";
  // top_all36_cgenff.rtf, RESI MOBU (1-methoxybutadiene): the enol ether C is CG2D1O, outside the alternation; the C
  // double-bonded to it anchors the phase (CG2DC1), then CG2DC2 across the single bond, the terminal =CH2 CG2DC3
  const System mobu = mol({"C", "O", "C", "C", "C", "C"}, {3, 0, 1, 1, 1, 2}, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}});
  const TypingResult t = type_with("cgenff-toppar", mobu);
  EXPECT_EQ(t.types[2], "CG2D1O");
  EXPECT_EQ(t.types[3], "CG2DC1");
  EXPECT_EQ(t.types[4], "CG2DC2");
  EXPECT_EQ(t.types[5], "CG2DC3");
  // RESI PRAL (acrolein): a conjugated C alone in its system is CG2DC1
  const TypingResult a = type_with("cgenff-toppar", mol({"C", "C", "C", "O"}, {2, 1, 1, 0}, {{0, 1}, {1, 2}, {2, 3}}));
  EXPECT_EQ(a.types[1], "CG2DC1");
}

TEST(Typing, Cgenff2a5FromThePaperSI) {
  // CGenFF v. 2a5 as published with Vanommeslaeghe et al., J. Comput. Chem. 31, 671 (2010): in the repository
  const std::string lib = std::string(CAPS_SOURCE_DIR) + "/data/forcefields/cgenff-2a5-jcc2010.json";
  const FFDef def = load_forcefield(lib);
  EXPECT_EQ(def.types.size(), 140u);
  EXPECT_NE(def.source.find("doi:10.1002/jcc.21367"), std::string::npos);
  EXPECT_NE(def.source.find("CGenFF v. 2a5"), std::string::npos);
  {
    std::ifstream in(lib);
    const std::string text((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    EXPECT_NE(text.find("parameter values as published in the Supporting Information of Vanommeslaeghe et al. 2010; CAPS's conversion"), std::string::npos);
    EXPECT_EQ(text.find("DL_FIELD"), std::string::npos);
  }
  // 2a5's own typing where later releases added types: DMSO's S is SG3O1 (no SG3O3), a terminal alkyne C CG1T1 (no
  // CG1T2), phosphate P PG1 outside a pyrophosphate, PG2 in one, PG3 for its protonated P (RESI PPI2; no PG0)
  const TypingResult dmso = type_with("cgenff-2a5-jcc2010", mol({"C", "S", "O", "C"}, {3, 0, 0, 3}, {{0, 1}, {1, 2}, {1, 3}}));
  EXPECT_EQ(dmso.types[1], "SG3O1");
  EXPECT_EQ(dmso.types[2], "OG2P1");
  const TypingResult propyne = type_with("cgenff-2a5-jcc2010", mol({"C", "C", "C"}, {1, 0, 3}, {{0, 1}, {1, 2}}));
  EXPECT_EQ(propyne.types[0], "CG1T1");
  // CH3-O-P(O2)-O-P(O2)-OH
  const TypingResult ppi2 = type_with("cgenff-2a5-jcc2010", mol({"C", "O", "P", "O", "O", "O", "P", "O", "O", "O"}, {3, 0, 0, 0, 0, 0, 0, 1, 0, 0},
                                                                {{0, 1}, {1, 2}, {2, 3}, {2, 4}, {2, 5}, {3, 6}, {6, 7}, {6, 8}, {6, 9}}));
  EXPECT_EQ(ppi2.types[2], "PG2");
  EXPECT_EQ(ppi2.types[3], "OG303");
  EXPECT_EQ(ppi2.types[6], "PG3");
  const TypingResult mp = type_with("cgenff-2a5-jcc2010", mol({"C", "O", "P", "O", "O", "O"}, {3, 0, 0, 0, 0, 0}, {{0, 1}, {1, 2}, {2, 3}, {2, 4}, {2, 5}}));
  EXPECT_EQ(mp.types[2], "PG1");
  // toluene: every term from the file, Urey-Bradley on the ring angles, the methyl C's own 1-4 Lennard-Jones
  const System tol = mol({"C", "C", "C", "C", "C", "C", "C"}, {3, 0, 1, 1, 1, 1, 1}, {{0, 1}, {1, 2}, {2, 3}, {3, 4}, {4, 5}, {5, 6}, {6, 1}});
  FFDef d2 = def;
  if (d2.typing.empty()) load_typing(d2, std::string(CAPS_SOURCE_DIR) + "/data/typing/cgenff-2a5-jcc2010.typing.json");
  const TypingResult tt = assign_types(tol, d2);
  EXPECT_EQ(tt.untyped, 0);
  EXPECT_EQ(tt.types[0], "CG331");
  EXPECT_EQ(tt.types[1], "CG2R61");
  ParamReport rep;
  const ForceField ff = parameterize(tol, d2, tt.types, "gasteiger", &rep);
  EXPECT_TRUE(rep.missing.empty());
  EXPECT_EQ(rep.estimated_terms, 0);
  EXPECT_FALSE(ff.urey_bradley.empty());
  EXPECT_FALSE(ff.lj14_types.empty());
  // impropers only where the topology's IMPR records put them: none on toluene (the file's CG2R61 X X HGR61, which
  // no residue lists, used to go on every ring CH), the amide C and N of N-methylacetamide as RESI NMA lists them
  EXPECT_TRUE(def.improper_placement);
  EXPECT_TRUE(ff.impropers.empty());
  EXPECT_TRUE(ff.impropers_harmonic.empty());
  const System nma = mol({"C", "N", "C", "O", "C"}, {3, 1, 0, 0, 3}, {{0, 1}, {1, 2}, {2, 3}, {2, 4}});
  const TypingResult tn = assign_types(nma, d2);
  ASSERT_EQ(tn.untyped, 0);
  ParamReport rn;
  const ForceField fn = parameterize(nma, d2, tn.types, "gasteiger", &rn);
  ASSERT_EQ(fn.impropers_harmonic.size(), 2u);
  EXPECT_EQ(fn.impropers_harmonic[0].i, 1u);   // N  C  CAY(the N methyl)  HN: NG2S1 X X HGP1, centre first
  EXPECT_EQ(fn.impropers_harmonic[1].i, 2u);   // C  CH3  N  O: CG2O1 X X OG2D1
  EXPECT_DOUBLE_EQ(fn.impropers_harmonic[1].k2, 120.0);
  // the placement survives writing the library out and reading it back
  const std::string tmp = (std::filesystem::temp_directory_path() / "caps_cgenff2a5_roundtrip.json").string();
  save_forcefield(def, tmp);
  FFDef back = load_forcefield(tmp);
  std::filesystem::remove(tmp);
  EXPECT_EQ(back.improper_sites.size(), def.improper_sites.size());
  EXPECT_EQ(back.impropers_unlisted, def.impropers_unlisted);
  size_t exc = 0, exc_back = 0;
  for (const auto& st : def.improper_sites) exc += st.exceptions.size();
  for (const auto& st : back.improper_sites) exc_back += st.exceptions.size();
  EXPECT_EQ(exc_back, exc);
  EXPECT_GT(exc, 0u);
}

TEST(Typing, CgenffAndOpls) {
  if (!built_here("cgenff-toppar")) GTEST_SKIP() << "CGenFF is built locally from CHARMM's toppar (bench/ff/open/build_charmm.sh)";
  // CGenFF: conjugated diene CG2DC3 CG2DC2 CG2DC1 CG2DC3 (same type across a double bond, switched across the single)
  const TypingResult b = type_with("cgenff-toppar", butadiene());
  EXPECT_EQ(b.types[0], "CG2DC3");
  EXPECT_EQ(b.types[1], "CG2DC2");
  EXPECT_EQ(b.types[2], "CG2DC1");
  EXPECT_EQ(b.types[3], "CG2DC3");
  const TypingResult f = type_with("cgenff-toppar", furan());
  EXPECT_EQ(f.types[0], "OG2R50");
  EXPECT_EQ(f.types[1], "CG2R51");
  // OPLS-AA 2024: methyl acetate CH3 135, ester C 465, =O 466, -O- 467, O-CH3 468
  const System ma = mol({"C", "C", "O", "O", "C"}, {3, 0, 0, 0, 3}, {{0, 1}, {1, 2}, {1, 3}, {3, 4}});
  const TypingResult o = type_with("oplsaa2024", ma);
  EXPECT_EQ(o.types[0], "opls_135");
  EXPECT_EQ(o.types[1], "opls_465");
  EXPECT_EQ(o.types[2], "opls_466");
  EXPECT_EQ(o.types[3], "opls_467");
  EXPECT_EQ(o.types[4], "opls_468");
}

#include "caps/example_typing.hpp"
#include "caps/io.hpp"
#include "caps/polymer.hpp"

// Typing by example: a polystyrene trimer typed with OPLS-AA 2024's rules teaches the types; the whole melt (ten
// chains of another length, other conformations) gets exactly the types the rules give it; a type changed by hand on
// one body-unit atom reaches every atom of that environment
TEST(Typing, ByExampleFromATrimer) {
  FFDef ff = load_forcefield(std::string(CAPS_SOURCE_DIR) + "/data/forcefields/oplsaa2024-moltemplate.json");
  load_typing(ff, std::string(CAPS_SOURCE_DIR) + "/data/typing/oplsaa2024-moltemplate.typing.json");
  ChainSpec c;
  RepeatUnit u;
  u.smiles = "*CC(*)c1ccccc1";
  u.name = "styrene";
  c.units = {u};
  c.dp = 3;
  GrowOptions g;
  g.chains = 1;
  g.density = 0.05;
  g.seed = 1;
  g.auto_scale = true;
  const System tri = grow_chains(c, g);
  const TypingResult tt = assign_types(tri, ff);
  ASSERT_EQ(tt.untyped, 0);
  const ExampleTypes learned = learn_types(tri, tt.types);
  EXPECT_TRUE(learned.conflicts.empty());
  const Trajectory melt = open_file(std::string(CAPS_SOURCE_DIR) + "/samples/ps_melt.data");
  const System m = melt.frame(0);
  const TypingResult mt = assign_types(m, ff);
  const ExampleMatch em = apply_types(m, learned);
  EXPECT_EQ(em.unmatched, 0u);
  size_t same = 0;
  for (size_t i = 0; i < m.atoms.size(); ++i) same += em.types[i] == mt.types[i];
  EXPECT_EQ(same, m.atoms.size());
  // by hand: the body unit's backbone CH carbon (a unit with neighbours on both sides) given another type
  std::vector<std::string> edited = tt.types;
  uint32_t ch = 0;
  for (uint32_t i = 0; i < tri.atoms.size(); ++i)
    if (tri.atoms[i].resid == 2 && tri.atoms[i].element == 6 && tt.types[i].rfind("515_", 0) == 0) ch = i;
  ASSERT_NE(ch, 0u);
  const auto eq = equivalent_atoms(tri, ch, learned.radius);
  for (uint32_t i : eq) edited[i] = "CUSTOM";
  const ExampleMatch em2 = apply_types(m, learn_types(tri, edited));
  size_t custom = 0;
  for (const auto& t : em2.types) custom += t == "CUSTOM";
  EXPECT_EQ(custom, 10u * 6);   // every backbone CH with a unit on both sides: six in each 8-unit chain
}

#include "caps/molecule.hpp"

// OPLS 2005's own charges: charge keys (the typing file's charge_rules) and its bond charge increments give the ester
// charges of PMMA's side group (C=O 0.51, O= −0.43, O −0.33, OCH3 0.16 with H 0.03), as the reference assigns them
TEST(Typing, Opls2005ChargesFromBondIncrements) {
  if (!built_here("opls2005")) GTEST_SKIP() << "OPLS 2005 has no open source CAPS can build it from; a local copy is needed";
  FFDef ff = load_forcefield(std::string(CAPS_SOURCE_DIR) + "/data/forcefields/opls2005.json");
  load_typing(ff, std::string(CAPS_SOURCE_DIR) + "/data/typing/opls2005.typing.json");
  BuildOptions bo;
  bo.forcefield = "uff";
  const System s = build_molecule("COC(=O)C(C)(C)C", bo).system;   // methyl pivalate
  const TypingResult tr = assign_types(s, ff);
  ASSERT_EQ(tr.untyped, 0);
  ParamReport rep;
  const ForceField f = parameterize(s, ff, tr.types, "types", &rep, false);
  EXPECT_TRUE(rep.missing.empty()) << rep.missing.front();
  double net = 0;
  for (size_t i = 0; i < s.atoms.size(); ++i) net += f.charge[i];
  EXPECT_NEAR(net, 0.0, 1e-9);
  auto q_of = [&](const std::string& type, int element) {
    for (size_t i = 0; i < s.atoms.size(); ++i)
      if (tr.types[i] == type && s.atoms[i].element == element) return f.charge[i];
    return std::nan("");
  };
  EXPECT_NEAR(q_of("CO4", 6), 0.51, 1e-9);
  EXPECT_NEAR(q_of("O", 8), -0.43, 1e-9);
  EXPECT_NEAR(q_of("OES", 8), -0.33, 1e-9);
  EXPECT_NEAR(f.charge[0], 0.16, 1e-9);   // the methoxy C (charge key 181)
  EXPECT_EQ(rep.charge_keys[0], "181");
}

// OPLS-AA 2024 with OPLS 2005's charges: 2-ethylpyridine's fixed OPLS-AA 2024 charges do not balance (the ring C1 at
// +0.473 expects a substituent worth an H); OPLS 2005 types the same molecule and its bond increments give charges that
// add up to zero, used with OPLS-AA 2024's own types
TEST(Typing, CompanionChargesFromOpls2005) {
  if (!built_here("opls2005")) GTEST_SKIP() << "OPLS 2005 has no open source CAPS can build it from; a local copy is needed";
  const std::string lib = std::string(CAPS_SOURCE_DIR) + "/data/forcefields/";
  FFDef ff = load_forcefield(lib + "oplsaa2024-moltemplate.json");
  load_typing(ff, std::string(CAPS_SOURCE_DIR) + "/data/typing/oplsaa2024-moltemplate.typing.json");
  ASSERT_EQ(ff.charge_increments_from, "opls2005.json");
  BuildOptions bo;
  bo.forcefield = "uff";
  const System s = build_molecule("CC(C)c1ccccn1", bo).system;   // 2-isopropylpyridine
  const TypingResult tr = assign_types(s, ff);
  ASSERT_EQ(tr.untyped, 0);
  EXPECT_EQ(std::count_if(tr.types.begin(), tr.types.end(), [](const std::string& t) { return t.rfind("520_", 0) == 0; }), 1);
  std::string note;
  const std::vector<double> q = companion_charges(s, ff, lib + "oplsaa2024-moltemplate.json", &note);
  ASSERT_EQ(q.size(), s.atoms.size());
  double net = 0;
  for (double v : q) net += v;
  EXPECT_NEAR(net, 0.0, 1e-9);
  EXPECT_NE(note.find("OPLS 2005"), std::string::npos) << note;
  // the same charges OPLS 2005 gives on its own
  FFDef o5 = load_forcefield(lib + "opls2005.json");
  load_typing(o5, std::string(CAPS_SOURCE_DIR) + "/data/typing/opls2005.typing.json");
  const TypingResult t5 = assign_types(s, o5);
  ParamReport r5;
  const ForceField f5 = parameterize(s, o5, t5.types, "types", &r5, false);
  for (size_t i = 0; i < q.size(); ++i) EXPECT_NEAR(q[i], f5.charge[i], 1e-12);
}

#include "caps/molecule.hpp"

// OPLS-AA 2024 on maleic anhydride (MAH), maleic acid, the maleate esters and the other anhydrides: every atom typed,
// every term found (CAPS's QM fits for the maleate torsions, the CM-C-OH bend and the MAH ring torsion) and each
// molecule neutral (the anhydride's own caps_anh_* charges; the ester types 465-467 leave it at -0.17 e). An opened
// anhydride, the half-ester, takes the ordinary ester and acid types again
TEST(Typing, OplsAa2024MaleateComplete) {
  const FFDef def = load_forcefield(std::string(CAPS_SOURCE_DIR) + "/data/forcefields/oplsaa2024.json");
  for (const std::string smi : {"O=C1C=CC(=O)O1", "OC(=O)/C=C\\C(=O)O", "COC(=O)/C=C\\C(=O)O", "COC(=O)/C=C\\C(=O)OC",
                                "O=C1CCC(=O)O1", "CC(=O)OC(C)=O", "CCC(C)OC(=O)/C=C\\C(=O)OCCCCO"}) {
    System s = build_molecule(smi).system;
    std::string ch = "types";
    prepare_for_forcefield(s, def, ch);
    const TypingResult tr = assign_types(s, def);
    EXPECT_EQ(tr.untyped, 0) << smi;
    ParamReport rep;
    const ForceField f = parameterize(s, def, tr.types, ch, &rep, false);
    EXPECT_TRUE(rep.missing.empty()) << smi << ": " << (rep.missing.empty() ? "" : rep.missing[0]);
    double q = 0;
    for (double x : f.charge) q += x;
    EXPECT_NEAR(q, 0.0, 1e-9) << smi;
    const bool anhydride = smi.find("O1") != std::string::npos || smi == "CC(=O)OC(C)=O";
    const bool anh_types = std::any_of(tr.types.begin(), tr.types.end(), [](const std::string& t) { return t.rfind("caps_anh", 0) == 0; });
    EXPECT_EQ(anh_types, anhydride) << smi;   // the half-ester and the acid: ester and acid types only
  }
}
