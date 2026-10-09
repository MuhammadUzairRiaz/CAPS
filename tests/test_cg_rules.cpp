// Chemistry-aware coarse-grained mapping (caps/cg_rules.hpp): polyesters cut at their ester bonds into diol and diacid
// beads (PBS, PBSA, PBAT share B), mass conserved, chains end to end, strict naming, the mapping file read back, one type
// list for several systems, the exact-cover and anchor forms, and the streaming dump mapper on a shuffled, wrapped dump.
#include <gtest/gtest.h>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <sstream>

#include "caps/cg_rules.hpp"
#include "caps/elements.hpp"
#include "caps/polymer.hpp"
#include "caps/typing.hpp"

using namespace caps;

namespace {
const std::string kData = std::string(CAPS_SOURCE_DIR) + "/data";
const std::string kPBS = "[*]OCCCCOC(=O)CCC(=O)[*]";
const std::string kBA = "[*]OCCCCOC(=O)CCCCC(=O)[*]";
const std::string kBT = "[*]OCCCCOC(=O)c1ccc(cc1)C(=O)[*]";

System grow(std::vector<std::string> units, Sequence seq, int dp, int chains, const std::string& pattern = "", std::vector<double> w = {},
            const std::string& tail = "hydroxyl") {
  ChainSpec c;
  for (auto& u : units) c.units.push_back({u, u});
  c.sequence = seq;
  c.dp = dp;
  c.pattern = pattern;
  c.weights = w;
  c.tail_cap = tail;
  GrowOptions g;
  g.chains = chains;
  g.density = 0.1;
  g.auto_scale = true;
  g.seed = 7;
  return grow_chains(c, g);
}

CgRules ester_cut() {
  std::ifstream f(kData + "/cg/mapping_rules.json");
  std::stringstream ss;
  ss << f.rdbuf();
  const Json lib = Json::parse(ss.str());
  for (const auto& p : lib["presets"].items())
    if (p.text("id") == "ester-cut") return cg_rules_from_json(p);
  throw std::runtime_error("no ester-cut preset");
}

int count(const CgMapping& m, const std::string& k) { return int(std::count(m.bead_kind.begin(), m.bead_kind.end(), k)); }

bool bonded(const CgMapping& m, int x, int y) {
  return std::find(m.bonds.begin(), m.bonds.end(), std::make_pair(std::min(x, y), std::max(x, y))) != m.bonds.end();
}
}  // namespace

TEST(CgRules, PbsTwoBeadsPerUnitMassConserved) {
  const System aa = grow({kPBS}, Sequence::Homopolymer, 6, 2);
  const CgMapping m = cg_mapping(aa, ester_cut());
  EXPECT_EQ(m.beads(), 24u);
  EXPECT_EQ(count(m, "B"), 12);
  EXPECT_EQ(count(m, "S"), 12);
  EXPECT_LT(std::fabs(m.mass_cg - m.mass_aa), 1e-6);
  ASSERT_EQ(m.chains.size(), 2u);
  for (const auto& ch : m.chains) {
    ASSERT_EQ(ch.size(), 12u);
    for (size_t k = 0; k + 1 < ch.size(); ++k) {
      EXPECT_TRUE(bonded(m, ch[k], ch[k + 1])) << "chain order";
      EXPECT_NE(m.bead_kind[size_t(ch[k])], m.bead_kind[size_t(ch[k + 1])]) << "B and S alternate";
    }
  }
  // end groups stay in their bead: the acid end's S bead holds one more O and H than an inner S
  const CgTypes t = cg_types_of(m);
  EXPECT_EQ(t.bonds, (std::vector<std::string>{"B-S"}));
  EXPECT_EQ(t.angles, (std::vector<std::string>{"B-S-B", "S-B-S"}));
  EXPECT_EQ(t.dihedrals, (std::vector<std::string>{"B-S-B-S"}));
  double smin = 1e9, smax = 0;
  for (size_t b = 0; b < m.beads(); ++b)
    if (m.bead_kind[b] == "S") smin = std::min(smin, m.bead_mass[b]), smax = std::max(smax, m.bead_mass[b]);
  EXPECT_NEAR(smin, 84.07, 0.02);                  // –C(=O)CH2CH2C(=O)–
  EXPECT_NEAR(smax - smin, 15.999 + 1.008, 0.02);  // –C(=O)CH2CH2C(=O)OH at the acid end
}

TEST(CgRules, CopolyesterCompositions) {
  const System pbsa = grow({kPBS, kBA}, Sequence::Pattern, 10, 2, "AAAAB");
  const CgMapping a = cg_mapping(pbsa, ester_cut());
  EXPECT_EQ(count(a, "B"), 20);
  EXPECT_EQ(count(a, "S"), 16);   // 80 / 20
  EXPECT_EQ(count(a, "A"), 4);
  const System pbat = grow({kBA, kBT}, Sequence::Random, 30, 2, "", {0.56, 0.44});
  const CgMapping t = cg_mapping(pbat, ester_cut());
  int ring_c = 0;
  const Perception p = perceive(pbat);
  for (size_t i = 0; i < pbat.atoms.size(); ++i) ring_c += p.aromatic[i] ? 1 : 0;
  EXPECT_EQ(count(t, "T") * 6, ring_c);   // one terephthalate bead per ring
  EXPECT_EQ(count(t, "A") + count(t, "T"), 60);
  EXPECT_EQ(count(t, "B"), 60);
  EXPECT_LT(std::fabs(t.mass_cg - t.mass_aa), 1e-6);
  // one list for the three: identical numbering
  const CgTypes all = cg_types_union({cg_types_of(cg_mapping(grow({kPBS}, Sequence::Homopolymer, 4, 1), ester_cut())), cg_types_of(a), cg_types_of(t)});
  EXPECT_EQ(all.beads, (std::vector<std::string>{"A", "B", "S", "T"}));
  EXPECT_EQ(all.bonds, (std::vector<std::string>{"A-B", "B-S", "B-T"}));
  std::vector<Vec3> pos;
  for (const auto& at : pbsa.atoms) pos.push_back(at.pos);
  const System cg = cg_structure(a, pbsa, cg_positions(a, pbsa, pos, pbsa.cell), pbsa.cell, all);
  for (const auto& b : cg.atoms) EXPECT_EQ(all.beads[size_t(b.type - 1)], b.name);
}

TEST(CgRules, OtherPolyestersAndStrictNaming) {
  const CgMapping pcl = cg_mapping(grow({"[*]OCCCCCC(=O)[*]"}, Sequence::Homopolymer, 8, 1), ester_cut());
  EXPECT_EQ(count(pcl, "C"), 8);
  const CgMapping pla = cg_mapping(grow({"[*]OC(C)C(=O)[*]"}, Sequence::Homopolymer, 8, 1), ester_cut());
  EXPECT_EQ(count(pla, "L"), 8);
  // naming rules that leave the terephthalate unnamed: reported, not merged
  CgRules r = ester_cut();
  r.names = {{"B", "[OX2][CH2][CH2][CH2][CH2][OX2]"}, {"A", "[CX3](=O)[CH2][CH2][CH2][CH2][CX3]=O"}};
  const System pbat = grow({kBA, kBT}, Sequence::Alternating, 4, 1);
  try {
    cg_mapping(pbat, r);
    FAIL() << "an unnamed fragment must stop the mapping";
  } catch (const std::invalid_argument& e) {
    EXPECT_NE(std::string(e.what()).find("C8H4O2"), std::string::npos) << e.what();
  }
  r.strict = false;
  const CgMapping loose = cg_mapping(pbat, r);
  EXPECT_EQ(count(loose, "F1") + count(loose, "F2"), 2);   // the two terephthalates (one with the acid end) by class
}

TEST(CgRules, MapFileRoundTripAndOtherForms) {
  const System aa = grow({kPBS}, Sequence::Homopolymer, 5, 2);
  const CgMapping m = cg_mapping(aa, ester_cut());
  const CgMapping back = cg_mapping_from_json(Json::parse(cg_mapping_json(m, aa).dump(1)), aa);
  EXPECT_EQ(back.bead_atoms, m.bead_atoms);
  EXPECT_EQ(back.bead_kind, m.bead_kind);
  EXPECT_EQ(back.chains, m.chains);
  EXPECT_EQ(back.bonds, m.bonds);
  EXPECT_EQ(back.dihedrals.size(), m.dihedrals.size());
  // the explicit form gives the same beads
  CgRules ex;
  ex.explicit_bead = m.bead_of();
  ex.strict = false;
  EXPECT_EQ(cg_mapping(aa, ex).bead_atoms, m.bead_atoms);
  // exact cover by fragment SMARTS (no acid end: every heavy atom in one fragment)
  const System plain = grow({kPBS}, Sequence::Homopolymer, 5, 1, "", {}, "");
  CgRules cover;
  cover.beads = {{"B", "[OX2][CH2][CH2][CH2][CH2][OX2]"}, {"S", "[CX3](=O)[CH2][CH2][CX3]=O"}};
  const CgMapping c = cg_mapping(plain, cover);
  EXPECT_EQ(count(c, "B"), 5);
  EXPECT_EQ(count(c, "S"), 5);
  // a bead on one atom
  CgRules on = ester_cut();
  on.position = "atom:[CX3]=O";
  const CgMapping a = cg_mapping(aa, on);
  std::vector<Vec3> pos;
  for (const auto& at : aa.atoms) pos.push_back(at.pos);
  const auto beads = cg_positions(a, aa, pos, aa.cell);
  for (size_t b = 0; b < a.beads(); ++b)
    if (a.bead_kind[b] == "S") {
      ASSERT_GE(a.bead_anchor[b], 0);
      EXPECT_EQ(aa.atoms[size_t(a.bead_anchor[b])].element, 6);
      EXPECT_LT(norm(beads[b] - pos[size_t(a.bead_anchor[b])]), 1e-12);
    }
}

TEST(CgRules, StreamingDumpShuffledAndWrapped) {
  const System aa = grow({kPBS}, Sequence::Homopolymer, 8, 3);
  const CgMapping m = cg_mapping(aa, ester_cut());
  const auto dir = std::filesystem::temp_directory_path() / "caps_cg_rules_test";
  std::filesystem::create_directories(dir);
  const std::string dump = (dir / "aa.lammpstrj").string(), out = (dir / "cg.lammpstrj").string();
  const Cell& c = aa.cell;
  ASSERT_TRUE(c.valid());
  {
    std::ofstream f(dump);
    for (int fr = 0; fr < 3; ++fr) {
      f << "ITEM: TIMESTEP\n" << fr * 100 << "\nITEM: NUMBER OF ATOMS\n" << aa.atoms.size() << "\nITEM: BOX BOUNDS pp pp pp\n";
      for (int k = 0; k < 3; ++k) f << c.origin[size_t(k)] << " " << c.origin[size_t(k)] + (k == 0 ? c.a[0] : k == 1 ? c.b[1] : c.c[2]) << "\n";
      f << "ITEM: ATOMS type z x y id\n";   // any column order, atoms backwards, wrapped, no image flags
      for (size_t i = aa.atoms.size(); i-- > 0;) {
        const Vec3 w = c.wrap(aa.atoms[i].pos);
        f << 1 << " " << w[2] << " " << w[0] << " " << w[1] << " " << aa.atoms[i].id << "\n";
      }
    }
  }
  CgTrajectoryOptions o;
  o.stride = 2;   // frames 1 and 3
  const auto rep = map_lammps_dump(m, aa, dump, out, o);
  EXPECT_EQ(rep.frames_read, 3u);
  EXPECT_EQ(rep.frames_written, 2u);
  EXPECT_FALSE(rep.unwrapped_input);
  std::vector<Vec3> pos;
  for (const auto& at : aa.atoms) pos.push_back(at.pos);
  const auto want = cg_positions(m, aa, pos, c);
  std::ifstream f(out);
  std::string line;
  std::vector<Vec3> got;
  while (std::getline(f, line))
    if (line.rfind("ITEM: ATOMS", 0) == 0) {
      for (size_t b = 0; b < m.beads(); ++b) {
        std::getline(f, line);
        std::istringstream s(line);
        int id, mol, type;
        Vec3 r;
        s >> id >> mol >> type >> r[0] >> r[1] >> r[2];
        got.push_back(r);
      }
      break;
    }
  ASSERT_EQ(got.size(), m.beads());
  for (size_t b = 0; b < m.beads(); ++b) EXPECT_LT(norm(c.minimum_image(got[b] - want[b])), 2e-3);   // the same beads (4 decimals)
  for (const auto& [x, y] : m.bonds) EXPECT_LT(norm(got[size_t(y)] - got[size_t(x)]), 10.0) << "chains whole along their bonds";
  std::filesystem::remove_all(dir);
}
