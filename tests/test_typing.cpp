#include <algorithm>
#include <cmath>
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
    FFDef ff = load_forcefield(std::string(CAPS_SOURCE_DIR) + "/data/forcefields/" + lib + ".json");
    load_typing(ff, std::string(CAPS_SOURCE_DIR) + "/data/typing/" + lib + ".typing.json");
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
  const TypingResult f = type_with("gaff-amber16", furan());
  EXPECT_EQ(f.types[0], "os");
  EXPECT_EQ(f.types[1], "cc");
  EXPECT_EQ(f.types[2], "cd");
  EXPECT_EQ(f.types[3], "cd");
  EXPECT_EQ(f.types[4], "cc");
  EXPECT_EQ(f.types[5], "h4");
  EXPECT_EQ(f.types[6], "ha");
  const TypingResult b = type_with("gaff-amber16", butadiene());
  EXPECT_EQ(b.types[0], "c2");
  EXPECT_EQ(b.types[1], "ce");
  EXPECT_EQ(b.types[2], "ce");
  EXPECT_EQ(b.types[3], "c2");
  // GAFF2 amine N by hydrogen count: methylamine n8
  const System ma = mol({"C", "N"}, {3, 2}, {{0, 1}});
  EXPECT_EQ(type_with("gaff-amber25", ma).types[1], "n8");
  EXPECT_EQ(type_with("gaff-amber16", ma).types[1], "n3");
}

TEST(Typing, CgenffAndOpls) {
  // CGenFF: conjugated diene CG2DC3 CG2DC2 CG2DC1 CG2DC3 (same type across a double bond, switched across the single)
  const TypingResult b = type_with("cgenff", butadiene());
  EXPECT_EQ(b.types[0], "CG2DC3");
  EXPECT_EQ(b.types[1], "CG2DC2");
  EXPECT_EQ(b.types[2], "CG2DC1");
  EXPECT_EQ(b.types[3], "CG2DC3");
  const TypingResult f = type_with("cgenff", furan());
  EXPECT_EQ(f.types[0], "OG2R50");
  EXPECT_EQ(f.types[1], "CG2R51");
  // OPLS-AA (DL_FIELD names): methyl acetate CT CO4 O OES CT
  const System ma = mol({"C", "C", "O", "O", "C"}, {3, 0, 0, 0, 3}, {{0, 1}, {1, 2}, {1, 3}, {3, 4}});
  const TypingResult o = type_with("opls2005", ma);
  EXPECT_EQ(o.types[0], "CT");
  EXPECT_EQ(o.types[1], "CO4");
  EXPECT_EQ(o.types[2], "O");
  EXPECT_EQ(o.types[3], "OES");
  EXPECT_EQ(o.types[4], "CT");
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
