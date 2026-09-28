#include <gtest/gtest.h>

#include <cmath>

#include "caps/analysis.hpp"
#include "caps/ffdef.hpp"
#include "caps/ffmerge.hpp"
#include "caps/io.hpp"
#include "caps/typing.hpp"

using namespace caps;

namespace {

const std::string kFF = std::string(CAPS_SOURCE_DIR) + "/data/forcefields/";

// the part of s made of the atoms given (bonds among them), and its force field
System part_of(const System& s, const std::vector<uint32_t>& atoms) {
  System p;
  p.cell = s.cell;
  std::vector<int> at(s.atoms.size(), -1);
  for (size_t k = 0; k < atoms.size(); ++k) { at[atoms[k]] = int(k); p.atoms.push_back(s.atoms[atoms[k]]); }
  for (const auto& b : s.bonds) if (at[b.i] >= 0 && at[b.j] >= 0) p.bonds.push_back({uint32_t(at[b.i]), uint32_t(at[b.j]), b.order});
  p.bonds_from_file = true;
  return p;
}

ForceField typed(const System& s, const FFDef& def) {
  const TypingResult tr = assign_types(s, def);
  EXPECT_EQ(tr.untyped, 0);
  ParamReport rep;
  return parameterize(s, def, tr.types, "gasteiger", &rep, false);
}

EnergyTerms energy(const ForceField& ff, const System& s) {
  EnergyOptions o;
  Evaluator ev(ff, o);
  std::vector<double> x, f;
  for (const auto& a : s.atoms) x.insert(x.end(), a.pos.begin(), a.pos.end());
  return ev.compute(x, s.cell, f);
}

}  // namespace

TEST(FFMerge, SameForceFieldInTwoPartsIsExact) {
  // the polystyrene melt split into chains 1-5 and 6-10, both GAFF2, merged with GAFF's own cross rule (ε geometric,
  // σ arithmetic): every term equals the single assignment's
  System s = open_file(std::string(CAPS_SOURCE_DIR) + "/samples/ps_melt.data").frame(0);
  make_molecules_whole(s);
  const FFDef def = load_forcefield(kFF + "gaff-amber25.json");
  const auto mol = s.molecules();
  std::vector<uint32_t> a, b;
  for (size_t i = 0; i < s.atoms.size(); ++i) (mol[i] < 5 ? a : b).push_back(uint32_t(i));
  // one charge model for both: Gasteiger–Marsili per chain equals per cell (it is per molecule)
  const ForceField whole = typed(s, def);
  const System sa = part_of(s, a), sb = part_of(s, b);
  const ForceField fa = typed(sa, def), fb = typed(sb, def);
  std::vector<std::string> notes;
  const ForceField m = merge_forcefields(s.atoms.size(), {{&fa, a, "A"}, {&fb, b, "B"}}, MergeOptions{}, &notes);
  const EnergyTerms e1 = energy(whole, s), e2 = energy(m, s);
  EXPECT_NEAR(e2.bond, e1.bond, 1e-9);
  EXPECT_NEAR(e2.angle, e1.angle, 1e-9);
  EXPECT_NEAR(e2.dihedral, e1.dihedral, 1e-9);
  EXPECT_NEAR(e2.improper, e1.improper, 1e-9);
  EXPECT_NEAR(e2.vdw, e1.vdw, 1e-8);
  EXPECT_NEAR(e2.coulomb, e1.coulomb, 1e-8);
  EXPECT_FALSE(notes.empty());
  EXPECT_EQ(m.type_names.size(), fa.type_names.size() + fb.type_names.size());   // B's names tagged
}

TEST(FFMerge, DifferentFamiliesAreRefusedOrMergedAsAsked) {
  System s = open_file(std::string(CAPS_SOURCE_DIR) + "/samples/ps_melt.data").frame(0);
  make_molecules_whole(s);
  const auto mol = s.molecules();
  std::vector<uint32_t> a, b;
  for (size_t i = 0; i < s.atoms.size(); ++i) (mol[i] < 5 ? a : b).push_back(uint32_t(i));
  const FFDef gaff = load_forcefield(kFF + "gaff-amber25.json"), opls = load_forcefield(kFF + "opls2005.json");
  const ForceField fa = typed(part_of(s, a), gaff), fb = typed(part_of(s, b), opls);
  // GAFF 1-4 Coulomb 1/1.2, OPLS 1/2: refused unless asked
  EXPECT_THROW(merge_forcefields(s.atoms.size(), {{&fa, a, "GAFF"}, {&fb, b, "OPLS"}}, MergeOptions{}), FieldError);
  MergeOptions o;
  o.scaling14 = "first";
  o.eps_rule = "arithmetic";
  std::vector<std::string> notes;
  const ForceField m = merge_forcefields(s.atoms.size(), {{&fa, a, "GAFF"}, {&fb, b, "OPLS"}}, o, &notes);
  // a cross pair: ε arithmetic, σ arithmetic of the two sites
  const int ta = m.type_index[a[0]], tb = m.type_index[b[0]];
  const PairType x = mixed_pair(m, ta, tb);
  EXPECT_NEAR(x.eps, 0.5 * (m.lj[size_t(ta)].eps + m.lj[size_t(tb)].eps), 1e-12);
  EXPECT_NEAR(x.sigma, 0.5 * (m.lj[size_t(ta)].sigma + m.lj[size_t(tb)].sigma), 1e-12);
  // OPLS's own unlike pairs keep OPLS's geometric σ
  if (fb.type_names.size() > 1) {
    const int t0 = m.type_index[b[0]];
    int t1 = t0;
    for (uint32_t i : b) if (m.type_index[i] != t0) { t1 = m.type_index[i]; break; }
    const PairType y = mixed_pair(m, t0, t1);
    EXPECT_NEAR(y.sigma, std::sqrt(m.lj[size_t(t0)].sigma * m.lj[size_t(t1)].sigma), 1e-12);
  }
  EXPECT_TRUE(std::any_of(notes.begin(), notes.end(), [](const std::string& n) { return n.find("1-4 scaling") != std::string::npos; }));
  EXPECT_TRUE(std::isfinite(energy(m, s).total()));
  // an atom left out
  std::vector<uint32_t> short_b(b.begin(), b.end() - 1);
  EXPECT_THROW(merge_forcefields(s.atoms.size(), {{&fa, a, "GAFF"}, {&fb, short_b, "OPLS"}}, o), FieldError);
}
