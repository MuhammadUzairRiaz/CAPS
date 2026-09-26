#include <gtest/gtest.h>

#include <cmath>
#include <fstream>
#include <map>
#include <regex>
#include <set>
#include <sstream>
#include <string>
#include <tuple>

#include "caps/ffdef.hpp"
#include "caps/field.hpp"
#include "caps/io.hpp"
#include "caps/martini_protein.hpp"
#include "caps/peptide.hpp"
#include "caps/typing.hpp"

using namespace caps;

namespace {

const std::string kData = std::string(CAPS_SOURCE_DIR) + "/data/martini/martini22-protein.json";
const std::string kRef = std::string(CAPS_SOURCE_DIR) + "/tests/data/vermouth/";

// the sections of a GROMACS .itp as sorted tuples (indices from 1, parameters rounded as martinize writes them)
std::map<std::string, std::set<std::string>> itp_terms(const std::string& text) {
  std::map<std::string, std::set<std::string>> out;
  std::istringstream in(text);
  std::string sec, line;
  while (std::getline(in, line)) {
    line = line.substr(0, line.find(';'));
    std::smatch m;
    if (std::regex_search(line, m, std::regex(R"(\[\s*(\S+)\s*\])"))) { sec = m[1]; continue; }
    std::istringstream ls(line);
    std::vector<std::string> w;
    for (std::string x; ls >> x;) w.push_back(x);
    if (w.empty() || sec == "moleculetype") continue;
    auto num = [&](size_t k, int digits) { std::ostringstream o; o.setf(std::ios::fixed); o.precision(digits); o << std::stod(w[k]); return o.str(); };
    std::string key;
    if (sec == "atoms") key = w[0] + " " + w[1] + " " + num(6, 2);
    else if (sec == "bonds") key = std::to_string(std::min(std::stoi(w[0]), std::stoi(w[1]))) + "-" + std::to_string(std::max(std::stoi(w[0]), std::stoi(w[1]))) + " " + num(3, 4) + " " + num(4, 0);
    else if (sec == "constraints") key = std::to_string(std::min(std::stoi(w[0]), std::stoi(w[1]))) + "-" + std::to_string(std::max(std::stoi(w[0]), std::stoi(w[1]))) + " " + num(3, 4);
    else if (sec == "angles") {
      int i = std::stoi(w[0]), k = std::stoi(w[2]);
      if (i > k) std::swap(i, k);
      key = std::to_string(i) + "-" + w[1] + "-" + std::to_string(k) + " " + num(4, 1) + " " + num(5, 1);
    } else if (sec == "dihedrals") key = w[0] + "-" + w[1] + "-" + w[2] + "-" + w[3] + " " + w[4] + " " + num(5, 1) + " " + num(6, 1);
    else continue;
    out[sec].insert(key);
  }
  return out;
}

std::string slurp(const std::string& p) {
  std::ifstream f(p);
  std::stringstream s;
  s << f.rdbuf();
  return s.str();
}

}  // namespace

// DSSP against DSSP 2.0 (CMBI) on vermouth's test proteins: a β-sheet with a turn and bends, and the Trp-cage (α and 3₁₀)
TEST(Dssp, MatchesDssp2OnReferenceProteins) {
  EXPECT_EQ(dssp(open_file(kRef + "1ico_aa.pdb").frame(0)), "CEEEEEETTEEEEEECCCCCCTTCEEEEC");
  EXPECT_EQ(dssp(open_file(kRef + "trpcage_aa.pdb").frame(0)), "CHHHHHHHTTGGGGTCCCCC");
  // vermouth's conversion to Martini's codes: helix ends and short helices
  EXPECT_EQ(dssp_to_martini("CHHHHHHHTTGGGGTCCCCC", kData), "C1113222TT3333TCCCCC");
  EXPECT_EQ(dssp_to_martini("CHHHHHHHHC", kData), "C11112222C");
}

// The Martini 2.2 protein of 1ICO is martinize2's, term by term (bead types, charges, bonds, constraints — here stiff bonds,
// written back as constraints —, angles, dihedrals, the disulfide) and bead by bead
TEST(MartiniProtein, MatchesMartinize2) {
  MartiniProteinReport rep;
  const System cg = martini22_protein(open_file(kRef + "1ico_aa.pdb").frame(0), "", kData, &rep);
  EXPECT_EQ(rep.beads, 67);
  EXPECT_EQ(rep.disulfides, 1);
  const auto mine = itp_terms(martini_itp(cg, kData)), ref = itp_terms(slurp(kRef + "1ico_martini22.itp"));
  for (const char* sec : {"atoms", "bonds", "constraints", "angles", "dihedrals"}) {
    ASSERT_TRUE(ref.count(sec)) << sec;
    EXPECT_EQ(mine.at(sec), ref.at(sec)) << sec;
  }
  const System pos = open_file(kRef + "1ico_martini22_cg.pdb").frame(0);
  ASSERT_EQ(pos.atoms.size(), cg.atoms.size());
  double worst = 0;
  for (size_t k = 0; k < cg.atoms.size(); ++k) worst = std::max(worst, norm(cg.atoms[k].pos - pos.atoms[k].pos));
  EXPECT_LT(worst, 0.005);
  // Martini 2.2's own bead types and pairs (martini_v2.2.itp): every bead typed and every pair found, AC1 / AC2 included
  FFDef def = load_forcefield(std::string(CAPS_SOURCE_DIR) + "/data/forcefields/martini22-proteins.json");
  System s = open_file(kRef + "1ico_aa.pdb").frame(0);
  std::string ch = "auto";
  prepare_for_forcefield(s, def, ch);
  const TypingResult tr = assign_types(s, def);
  EXPECT_EQ(tr.untyped, 0);
  ParamReport pr;
  EXPECT_NO_THROW(parameterize(s, def, tr.types, ch, &pr, false));
  EXPECT_TRUE(pr.missing.empty());
}

// Helices by the .ff's rules on CAPS's own helical peptide: helix-end bead types, helix constraints, 96° / 700 BBB angles in
// the helix core, −120° / 400 backbone dihedrals; and the force field types, parameterises and evaluates it
TEST(MartiniProtein, HelixRules) {
  PeptideOptions po;
  po.sequence = "AEAAAKEAAAKEAAAKA";
  po.structure = std::string(po.sequence.size(), 'H');
  const System aa = build_peptide(po);
  MartiniProteinReport rep;
  const System cg = martini22_protein(aa, "", kData, &rep);
  EXPECT_EQ(rep.cg_ss.substr(1, 4), "1111");
  const auto t = itp_terms(martini_itp(cg, kData));
  int helix_constraints = 0, dihedrals = 0, core_angles = 0;
  for (const auto& c : t.at("constraints")) helix_constraints += c.find(" 0.3100") != std::string::npos;
  for (const auto& d : t.at("dihedrals")) dihedrals += d.find(" 1 -120.0 400.0") != std::string::npos;
  for (const auto& a : t.at("angles")) core_angles += a.find(" 96.0 700.0") != std::string::npos;
  EXPECT_GT(helix_constraints, 8);
  EXPECT_GT(dihedrals, 8);
  EXPECT_GT(core_angles, 5);
  FFDef def = load_forcefield(std::string(CAPS_SOURCE_DIR) + "/data/forcefields/martini22-proteins.json");
  System s = aa;
  std::string ch = "auto";
  prepare_for_forcefield(s, def, ch);
  ASSERT_TRUE(s.topology);
  const TypingResult tr = assign_types(s, def);
  EXPECT_EQ(tr.untyped, 0);
  ParamReport pr;
  const ForceField ff = parameterize(s, def, tr.types, ch, &pr, false);
  EXPECT_EQ(ff.dihedrals.size(), size_t(dihedrals));
  std::vector<double> x, f;
  for (const auto& a : s.atoms) x.insert(x.end(), {a.pos[0], a.pos[1], a.pos[2]});
  Evaluator ev(ff, EnergyOptions{});
  EXPECT_TRUE(std::isfinite(ev.compute(x, s.cell, f).total()));
}
