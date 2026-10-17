// Completion of an assignment by CAPS ML (caps/ffcomplete.hpp): the source a fitted rule carries follows its terms into
// the assigned force field, its saved form and the LAMMPS and GROMACS files. The provenance tests use a hand-made overlay
// with an ML source, so they run in every build; FfCompleteMl.* fit with the model and need a build with libtorch.
#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>
#include <algorithm>
#include <regex>
#include <set>
#include <sstream>

#include "caps/ffcomplete.hpp"
#include "caps/ffdef.hpp"
#include "caps/ffio.hpp"
#include "caps/molecule.hpp"
#include "caps/relax.hpp"
#include "caps/typing.hpp"

using namespace caps;

namespace {

std::string slurp(const std::string& p) {
  std::ifstream f(p);
  std::stringstream ss;
  ss << f.rdbuf();
  return ss.str();
}

std::filesystem::path temp_dir(const std::string& name) {
  auto d = std::filesystem::temp_directory_path() / ("caps_test_" + name);
  std::filesystem::create_directories(d);
  return d;
}

// 4-phenyl-1-butene: OPLS-AA 2024 has no CM-CT-CT-CA torsion (the one SBR misses about its styrene-butadiene junction)
struct Case {
  FFDef def;
  System s;
  std::vector<std::string> types;
};

Case phenylbutene() {
  Case c;
  c.def = load_forcefield(std::string(CAPS_SOURCE_DIR) + "/data/forcefields/oplsaa2024.json");
  c.s = build_molecule("C=CCCc1ccccc1").system;
  const TypingResult tr = assign_types(c.s, c.def);
  EXPECT_EQ(tr.untyped, 0);
  c.types = tr.types;
  return c;
}

// an overlay as CAPS ML writes one: a rule for each missing dihedral, with its source (values made up: no model here)
FFDef fake_ml_overlay(const ParamReport& rep) {
  FFDef ov;
  ov.name = "CAPS ML (test)";
  for (const auto& m : rep.missing) {
    std::istringstream in(m);
    std::string kind, w;
    in >> kind;
    if (kind != "dihedral") continue;
    FFRule r;
    while (in >> w) {
      const auto p = w.find('(');
      r.match.push_back(p == std::string::npos ? w : w.substr(p + 1, w.size() - p - 2));
    }
    r.name = "fake";
    r.style = "opls";
    r.params = {0.5, -0.25, 0.125, 0.0};
    r.comment = "CAPS ff complete (CAPS MACE v1, caps-mace-v1.model): test";
    r.source = ml_source("data/ml/caps-mace-v1.model", 0.314);
    ov.dihedrals.push_back(r);
  }
  return ov;
}

}  // namespace

TEST(FfComplete, OptionValuesAndSourceText) {
  EXPECT_TRUE(parse_complete_option("ml"));
  EXPECT_FALSE(parse_complete_option(""));
  EXPECT_FALSE(parse_complete_option("none"));
  EXPECT_THROW(parse_complete_option("qm"), std::invalid_argument);
  EXPECT_EQ(ml_source("/x/data/ml/caps-mace-v1.model", 0.314), "CAPS ML (CAPS MACE v1, wB97M-D3(BJ)/def2-TZVPPD reference), RMS 0.31 kcal/mol");
  EXPECT_TRUE(is_ml_source(ml_source("other.model", 1.0)));
}

// The source of an overlay's rule reaches every term it gives, survives the saved force field and the overlay file
// (and an older overlay without "source" gets it from its comment)
TEST(FfComplete, SourceFollowsTheTerms) {
  Case c = phenylbutene();
  ParamReport rep;
  (void)parameterize(c.s, c.def, c.types, "types", &rep, true);
  ASSERT_FALSE(rep.missing.empty()) << "OPLS-AA 2024 now has CM-CT-CT-CA: pick another missing torsion for this test";
  const FFDef ov = fake_ml_overlay(rep);
  ASSERT_FALSE(ov.dihedrals.empty());
  FFDef def = c.def;
  merge_forcefield(def, ov);
  ParamReport rep2;
  const ForceField f = parameterize(c.s, def, c.types, "types", &rep2, false);
  EXPECT_TRUE(rep2.missing.empty());
  ASSERT_EQ(f.term_sources.size(), ov.dihedrals.size());
  EXPECT_EQ(f.term_sources[0].kind, "dihedral");
  EXPECT_TRUE(is_ml_source(f.term_sources[0].source));
  EXPECT_FALSE(f.term_sources[0].terms.empty());
  const TermSourceIndex idx(f);
  auto t = f.term_sources[0].terms[0];
  EXPECT_NE(idx.find("dihedral", t), nullptr);
  std::reverse(t.begin(), t.end());
  EXPECT_NE(idx.find("dihedral", t), nullptr);   // either direction
  EXPECT_EQ(idx.find("angle", {t[0], t[1], t[2]}), nullptr);
  // the assigned force field written whole and read back
  const ForceField back = forcefield_from_json(forcefield_to_json(f));
  ASSERT_EQ(back.term_sources.size(), f.term_sources.size());
  EXPECT_EQ(back.term_sources[0].source, f.term_sources[0].source);
  EXPECT_EQ(back.term_sources[0].terms, f.term_sources[0].terms);
  // the overlay file keeps the source; an older one (comment only) gets it from the comment
  const auto dir = temp_dir("ffcomplete_src");
  save_forcefield(ov, (dir / "ov.json").string());
  const FFDef ov2 = load_forcefield((dir / "ov.json").string());
  ASSERT_FALSE(ov2.dihedrals.empty());
  EXPECT_EQ(ov2.dihedrals[0].source, ov.dihedrals[0].source);
  FFDef old = ov;   // as caps ff complete wrote overlays before rules carried "source": the comment says it
  for (auto& r : old.dihedrals) r.source.clear();
  save_forcefield(old, (dir / "old.json").string());
  const FFDef o2 = load_forcefield((dir / "old.json").string());
  ASSERT_FALSE(o2.dihedrals.empty());
  EXPECT_EQ(o2.dihedrals[0].source, "CAPS ML (CAPS MACE v1, wB97M-D3(BJ)/def2-TZVPPD reference)");   // no RMS in this comment
  std::filesystem::remove_all(dir);
}

// The LAMMPS data file and input mark CAPS ML's terms (each coefficient line, and a header block), the GROMACS .top and
// .itp too (each line of the term, and a header note); LAMMPS's six decimals stay
TEST(FfComplete, ExportersMarkMlTerms) {
  Case c = phenylbutene();
  ParamReport rep;
  (void)parameterize(c.s, c.def, c.types, "types", &rep, true);
  ASSERT_FALSE(rep.missing.empty());
  FFDef def = c.def;
  merge_forcefield(def, fake_ml_overlay(rep));
  const ForceField f = parameterize(c.s, def, c.types, "types", nullptr, false);
  const auto dir = temp_dir("ffcomplete_export");
  EnergyOptions e;
  System s = c.s;
  for (size_t i = 0; i < s.atoms.size(); ++i) s.atoms[i].charge = f.charge[i];
  s.has_charges = true;
  LammpsStyle ls;
  ls.native = true;   // the force field's own styles (caps ff apply's default): OPLS torsions as dihedral_style opls
  write_lammps_data_ff(s, f, e, (dir / "m.data").string(), false, ls);
  write_lammps_input(s, f, e, "m.data", (dir / "m.in").string(), 0, true, LammpsRun{}, ls, nullptr);
  const std::string data = slurp((dir / "m.data").string());
  // header block before the counts
  const auto counts = data.find(" atoms\n");
  const auto head = data.find("# 1 term");
  ASSERT_NE(head, std::string::npos) << data.substr(0, 600);
  EXPECT_LT(head, counts);
  EXPECT_NE(data.find("CAPS ML (CAPS MACE v1, wB97M-D3(BJ)/def2-TZVPPD reference), RMS 0.31 kcal/mol"), std::string::npos);
  // the dihedral coefficient line of the fitted class: six decimals, then "# types  CAPS ML: classes, RMS ..."
  const std::string cls = f.term_sources[0].classes;
  const std::regex line(R"((\d+)( -?\d+\.\d{6}){4}  # \S+  CAPS ML: )" + std::regex_replace(cls, std::regex(R"([-])"), R"(\-)") + R"(, RMS 0\.31 kcal/mol)");
  const auto sec = data.find("Dihedral Coeffs");
  ASSERT_NE(sec, std::string::npos);
  const std::string block = data.substr(sec, data.find("\n\n", data.find("\n\n", sec) + 2) - sec);
  EXPECT_TRUE(std::regex_search(block, line)) << block;
  int marked = 0;
  for (std::istringstream in(block); in.good();) {
    std::string l;
    std::getline(in, l);
    marked += l.find("CAPS ML") != std::string::npos;
  }
  EXPECT_EQ(marked, 1);   // only the fitted class's type
  EXPECT_NE(slurp((dir / "m.in").string()).find("# 1 term"), std::string::npos);
  // GROMACS: the fitted dihedral's lines say so, the headers count it
  write_gromacs(s, f, e, (dir / "g").string());
  const std::string top = slurp((dir / "g.top").string()), itp = slurp((dir / "g.itp").string());
  EXPECT_NE(top.find("; 1 term of 1 parameter set fitted by CAPS ML"), std::string::npos) << top.substr(0, 400);
  EXPECT_NE(itp.find("; 1 term of 1 parameter set fitted by CAPS ML"), std::string::npos);
  const auto dsec = itp.find("[ dihedrals ]");
  ASSERT_NE(dsec, std::string::npos);
  int glines = 0;
  for (std::istringstream in(itp.substr(dsec)); in.good();) {
    std::string l;
    std::getline(in, l);
    if (l.find("; CAPS ML: " + cls + ", RMS 0.31 kcal/mol") != std::string::npos) ++glines;
  }
  EXPECT_EQ(glines, 3);   // V1, V2, V3 of the one fitted torsion (V4 is 0: no line)
  std::filesystem::remove_all(dir);
}

// Every coefficient GROMACS reads from CAPS's .top and .itp has six decimals (or six significant digits under 0.001), as
// the LAMMPS files have them
TEST(Gromacs, SixDecimalsEverywhere) {
  Case c = phenylbutene();
  FFDef def = c.def;
  ParamReport rep;
  (void)parameterize(c.s, def, c.types, "types", &rep, true);
  merge_forcefield(def, fake_ml_overlay(rep));
  const ForceField f = parameterize(c.s, def, c.types, "types", nullptr, false);
  const auto dir = temp_dir("gmx_six");
  EnergyOptions e;
  write_gromacs(c.s, f, e, (dir / "g").string());
  const std::regex six(R"(^-?\d+\.\d{6}$|^-?\d\.\d{6}e[-+]\d+$)");
  int checked = 0;
  for (const std::string file : {"g.top", "g.itp"}) {
    std::string sec;
    for (std::istringstream in(slurp((dir / file).string())); in.good();) {
      std::string l;
      std::getline(in, l);
      if (const auto p = l.find(';'); p != std::string::npos) l = l.substr(0, p);
      if (l.rfind("#", 0) == 0) continue;   // #include
      std::smatch m;
      if (std::regex_search(l, m, std::regex(R"(\[\s*(\S+)\s*\])"))) { sec = m[1]; continue; }
      std::istringstream ws(l);
      std::vector<std::string> w;
      for (std::string x; ws >> x;) w.push_back(x);
      if (w.empty()) continue;
      // the numbers with a fraction: every one a coefficient, charge or mass
      static const std::set<std::string> secs = {"defaults", "atomtypes", "nonbond_params", "pairs", "bonds", "angles", "dihedrals", "atoms"};
      if (!secs.count(sec)) continue;
      for (const auto& x : w)
        if (x.find('.') != std::string::npos) {
          EXPECT_TRUE(std::regex_match(x, six)) << file << " [" << sec << "] " << x << " in: " << l;
          ++checked;
        }
    }
  }
  EXPECT_GT(checked, 200);
  // charges at six decimals keep the molecule neutral
  double q = 0;
  const std::string itp = slurp((dir / "g.itp").string());
  bool in_atoms = false;
  for (std::istringstream in(itp); in.good();) {
    std::string l;
    std::getline(in, l);
    if (l.find("[ atoms ]") != std::string::npos) { in_atoms = true; continue; }
    if (l.rfind("[", 0) == 0 || (in_atoms && l.empty())) { if (in_atoms && l.empty()) break; continue; }
    if (!in_atoms || l[0] == ';') continue;
    std::istringstream ws(l);
    std::vector<std::string> w;
    for (std::string x; ws >> x;) w.push_back(x);
    if (w.size() >= 8) q += std::stod(w[6]);
  }
  EXPECT_NEAR(q, 0.0, 1e-9);
  std::filesystem::remove_all(dir);
}

// A build without CAPS's ML engine says so when completion is asked for, and fits nothing
TEST(FfComplete, WithoutTheEngineSaysSo) {
  if (ml_engine_available()) GTEST_SKIP() << "this build has the ML engine";
  Case c = phenylbutene();
  FFDef def = c.def;
  try {
    (void)fill_missing_with_ml(c.s, def, c.types, CompleteOptions{});
    FAIL() << "no error";
  } catch (const std::runtime_error& e) {
    EXPECT_NE(std::string(e.what()).find("built without its ML engine"), std::string::npos) << e.what();
  }
  EXPECT_EQ(def.dihedrals.size(), c.def.dihedrals.size());
}

#ifdef CAPS_WITH_TORCH
// With the model: the missing CM-CT-CT-CA torsion of 4-phenyl-1-butene is fitted, the assignment is complete and the
// term carries CAPS ML's source into the force field
TEST(FfCompleteMl, FillsTheMissingTorsion) {
  Case c = phenylbutene();
  FFDef def = c.def;
  CompleteOptions co;
  co.model = std::string(CAPS_SOURCE_DIR) + "/data/ml/caps-mace-v1.model";
  co.step = 30;   // a coarser scan: the test checks the flow, not the profile
  const MlFill f = fill_missing_with_ml(c.s, def, c.types, co);
  EXPECT_FALSE(f.report.missing_before.empty());
  EXPECT_GE(f.fitted, 1) << (f.left.empty() ? "" : f.left[0]);
  EXPECT_TRUE(f.report.complete);
  ParamReport rep;
  const ForceField ff = parameterize(c.s, def, c.types, "types", &rep, false);
  EXPECT_TRUE(rep.missing.empty());
  ASSERT_FALSE(ff.term_sources.empty());
  EXPECT_EQ(ff.term_sources[0].source.rfind("CAPS ML (CAPS MACE v1, wB97M-D3(BJ)/def2-TZVPPD reference), RMS ", 0), 0u) << ff.term_sources[0].source;
  EXPECT_FALSE(find_ml_model({std::string(CAPS_SOURCE_DIR) + "/data/forcefields/oplsaa2024.json"}).empty());
}
#endif
