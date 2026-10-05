#include <algorithm>
#include <gtest/gtest.h>

#include <cmath>
#include <map>

#include "caps/molecule.hpp"
#include "caps/peptide.hpp"

using namespace caps;

namespace {

// IUPAC dihedral a-b-c-d in degrees.
double dihedral(const Vec3& a, const Vec3& b, const Vec3& c, const Vec3& d) {
  const Vec3 b0 = a - b, b1 = (c - b) * (1.0 / norm(c - b)), b2 = d - c;
  const Vec3 v = b0 - b1 * dot(b0, b1), w = b2 - b1 * dot(b2, b1);
  return std::atan2(dot(cross(b1, v), w), dot(v, w)) * 180.0 / 3.14159265358979323846;
}

// Atoms of a residue by name, residues counted from the N atoms in order.
std::vector<std::map<std::string, Vec3>> residues(const System& s) {
  std::vector<std::map<std::string, Vec3>> out;
  for (const auto& a : s.atoms) {
    if (a.resname == "ACE" || a.resname == "NME") continue;
    if (a.name == "N") out.emplace_back();
    if (!out.empty() && !out.back().count(a.name)) out.back()[a.name] = a.pos;
  }
  return out;
}

double signed_volume(const Vec3& centre, const Vec3& a, const Vec3& b, const Vec3& c) { return dot(a - centre, cross(b - centre, c - centre)); }

}  // namespace

TEST(Peptide, HelixKeepsItsTorsions) {
  PeptideOptions o;
  o.sequence = "AEAAAKEAAAKA";
  o.structure = std::string(12, 'H');
  o.cleanup = false;
  PeptideReport rep;
  const System s = build_peptide(o, &rep);
  EXPECT_EQ(rep.residues, 12);
  EXPECT_EQ(rep.charge, 0);   // NH3+ and COO− with two Glu− and two Lys+
  const auto r = residues(s);
  ASSERT_EQ(r.size(), 12u);
  for (size_t i = 1; i + 1 < r.size(); ++i) {
    EXPECT_NEAR(dihedral(r[i - 1].at("C"), r[i].at("N"), r[i].at("CA"), r[i].at("C")), -57, 0.5) << i;
    EXPECT_NEAR(dihedral(r[i].at("N"), r[i].at("CA"), r[i].at("C"), r[i + 1].at("N")), -47, 0.5) << i;
    EXPECT_NEAR(std::fabs(dihedral(r[i].at("CA"), r[i].at("C"), r[i + 1].at("N"), r[i + 1].at("CA"))), 180, 0.5) << i;
  }
  // an α-helix rises 1.5 Å per residue: CA 1 to CA 11 about 15 Å
  EXPECT_NEAR(norm(r[10].at("CA") - r[0].at("CA")), 15.0, 1.0);
  // i → i+4 backbone hydrogen bonds: O(i) to N(i+4) about 2.9–3.1 Å
  EXPECT_NEAR(norm(r[2].at("O") - r[6].at("N")), 3.0, 0.35);
}

TEST(Peptide, EveryResidueIsL) {
  PeptideOptions o;
  o.sequence = "ARNDCQEGHILKMFPSTWYV";
  o.structure = "PPPPPPPPPPPPPPPPPPPP";
  o.cleanup = false;
  PeptideReport rep;
  const System s = build_peptide(o, &rep);
  EXPECT_EQ(rep.residues, 20);
  // the same handedness at every CA: N, C, CB (Gly has none)
  const auto r = residues(s);
  ASSERT_EQ(r.size(), 20u);
  double ref = 0;
  for (size_t i = 0; i < r.size(); ++i) {
    if (!r[i].count("CB")) continue;
    const double v = signed_volume(r[i].at("CA"), r[i].at("N"), r[i].at("C"), r[i].at("CB"));
    if (ref == 0) ref = v;
    EXPECT_GT(v * ref, 0) << "residue " << i;
  }
  // no two atoms on top of each other
  size_t close = 0;
  for (size_t i = 0; i < s.atoms.size(); ++i)
    for (size_t j = i + 1; j < s.atoms.size(); ++j) close += norm(s.atoms[i].pos - s.atoms[j].pos) < 0.7;
  EXPECT_EQ(close, 0u);
}

TEST(Peptide, SecondCentresOfIleAndThrMatchPubChem) {
  // the side chains as written in the builder against PubChem's isomeric SMILES of L-Ile (2S,3S) and L-Thr (2S,3R)
  auto cb_volume = [](const std::string& smiles, bool thr) {
    MolGraph g = parse_smiles(smiles);
    add_hydrogens(g);
    const auto p = embed(g);
    int n = -1, ca = -1;
    for (int i = 0; i < g.heavy; ++i) if (g.atoms[size_t(i)].element == 7) n = i;
    std::vector<std::vector<int>> nb(g.atoms.size());
    for (const auto& b : g.bonds) nb[size_t(b.a)].push_back(b.b), nb[size_t(b.b)].push_back(b.a);
    ca = nb[size_t(n)][0];
    int cb = -1;
    for (int q : nb[size_t(ca)]) {
      if (g.atoms[size_t(q)].element != 6) continue;
      bool carbonyl = false;
      for (int x : nb[size_t(q)]) carbonyl |= g.atoms[size_t(x)].element == 8 && nb[size_t(x)].size() == 1 && !thr ? true : (g.atoms[size_t(x)].element == 8 && std::count_if(nb[size_t(q)].begin(), nb[size_t(q)].end(), [&](int y) { return g.atoms[size_t(y)].element == 8; }) == 2);
      if (!carbonyl) cb = q;
    }
    // CB's heavy neighbours other than CA: Ile CG1 (ethyl: has a carbon neighbour) and CG2 (methyl); Thr OG1 and CG2
    int big = -1, methyl = -1;
    for (int q : nb[size_t(cb)]) {
      if (q == ca || g.atoms[size_t(q)].element == 1) continue;
      const bool is_methyl = g.atoms[size_t(q)].element == 6 && std::count_if(nb[size_t(q)].begin(), nb[size_t(q)].end(), [&](int y) { return g.atoms[size_t(y)].element != 1; }) == 1;
      if (is_methyl && methyl < 0) methyl = q; else big = q;
    }
    return signed_volume(p[size_t(cb)], p[size_t(ca)], p[size_t(big)], p[size_t(methyl)]);
  };
  EXPECT_GT(cb_volume("N[C@@H]([C@@H](C)CC)C(=O)O", false) * cb_volume("CC[C@H](C)[C@@H](C(=O)O)N", false), 0);
  EXPECT_GT(cb_volume("N[C@@H]([C@H](O)C)C(=O)O", true) * cb_volume("C[C@H]([C@@H](C(=O)O)N)O", true), 0);
}

TEST(Peptide, TerminiChargesAndCaps) {
  PeptideOptions o;
  o.sequence = "GKDG";
  o.cleanup = false;
  PeptideReport rep;
  build_peptide(o, &rep);
  EXPECT_EQ(rep.charge, 0);   // NH3+ … COO−, K+ D−
  o.ph = 12.8;                 // Lys neutral, Asp−: NH3+ COO− D−
  build_peptide(o, &rep);
  EXPECT_EQ(rep.charge, -1);   // the chosen NH3+ stays
  o.neutral = true;
  o.n_term = "ACE";
  o.c_term = "NME";
  const System capped = build_peptide(o, &rep);
  EXPECT_EQ(rep.charge, 0);
  int ace = 0, nme = 0;
  for (const auto& a : capped.atoms) ace += a.resname == "ACE", nme += a.resname == "NME";
  EXPECT_EQ(ace, 6);   // CH3 C O + 3 H
  EXPECT_EQ(nme, 6);   // N CH3 + H + 3 H
  EXPECT_THROW({ o.sequence = "GXB"; build_peptide(o); }, std::invalid_argument);
  EXPECT_EQ(parse_fasta(">sp|P1|test\nACDE\nFGH\n>second\nKLM\n"), "ACDEFGH");
  EXPECT_EQ(residue_name('W'), "TRP");
}

TEST(Peptide, CleanUpKeepsTheHelix) {
  PeptideOptions o;
  o.sequence = "AAAAAAAAAA";
  o.structure = "HHHHHHHHHH";
  o.n_term = "ACE";
  o.c_term = "NME";
  PeptideReport rep;
  const System s = build_peptide(o, &rep);
  bool cleaned = false;
  for (const auto& n : rep.notes) cleaned |= n.rfind("UFF clean-up: E", 0) == 0;
  EXPECT_TRUE(cleaned) << (rep.notes.empty() ? "" : rep.notes.front());
  const auto r = residues(s);
  ASSERT_EQ(r.size(), 10u);
  for (size_t i = 2; i + 2 < r.size(); ++i)
    EXPECT_NEAR(dihedral(r[i - 1].at("C"), r[i].at("N"), r[i].at("CA"), r[i].at("C")), -57, 25) << i;
}
