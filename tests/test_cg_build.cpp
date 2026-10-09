// Coarse-grained polymer builder (caps/cg_build.hpp): the walks follow the bonded distributions (a freely rotating chain's
// C∞ when the dihedrals are uniform), compositions and sequence statistics, chain lengths, the box from the density, the
// files, the equilibration deck.
#include <gtest/gtest.h>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <functional>

#include "caps/cg_build.hpp"

using namespace caps;

namespace {
constexpr double kB = 0.0019872067, kT = kB * 300, kPi = 3.14159265358979323846;

CgBondedTable table(int kind, const std::string& key, double lo, double hi, double bin, const std::function<double(double)>& U) {
  CgBondedTable t;
  t.kind = kind, t.key = key, t.sampled = true, t.count = 100000;
  for (double x = lo + bin / 2; x < hi; x += bin) t.xh.push_back(x), t.P.push_back(std::exp(-U(x) / kT));
  return t;
}

// A–B units: bonds A-B around 5 Å, angles around 120°, no dihedral tables (uniform)
CgBondedResult tables() {
  CgBondedResult r;
  r.bonds.push_back(table(0, "A-B", 0, 20, 0.02, [](double x) { return 3.0 * (x - 5) * (x - 5); }));
  for (const char* k : {"A-B-A", "B-A-B"}) r.angles.push_back(table(1, k, 0, 180, 1, [](double t) { return 0.002 * (t - 120) * (t - 120); }));
  return r;
}
const std::map<std::string, double> kMass = {{"A", 100}, {"B", 80}, {"C", 120}};
}  // namespace

TEST(CgBuild, WalksFollowTheDistributions) {
  CgBuildOptions o;
  o.units = {{"AB", {"A", "B"}}};
  o.dp = 300, o.chains = 30, o.density = 1.0, o.seed = 3;
  const CgBuildResult r = build_cg_polymer(o, tables(), kMass);
  ASSERT_EQ(r.beads.atoms.size(), 30u * 600u);
  // the box: total mass / density
  const double mass = 30 * 300 * 180.0;
  EXPECT_NEAR(r.box, std::cbrt(mass / 6.02214076e23 / 1.0 * 1e24), 1e-6);
  // bond lengths and angles as drawn
  double rs = 0, r2 = 0, cs = 0;
  long nb = 0, na = 0;
  for (const auto& ch : r.topology.chains)
    for (size_t i = 0; i + 1 < ch.size(); ++i) {
      const Vec3 d = r.beads.atoms[size_t(ch[i + 1])].pos - r.beads.atoms[size_t(ch[i])].pos;
      rs += norm(d), r2 += dot(d, d), ++nb;
      if (i + 2 < ch.size()) {
        const Vec3 e = r.beads.atoms[size_t(ch[i + 2])].pos - r.beads.atoms[size_t(ch[i + 1])].pos;
        cs += dot(d * -1.0, e) / (norm(d) * norm(e)), ++na;
      }
    }
  // the expected means of the distributions drawn from: r² e^(−U/kT) and sin θ e^(−U/kT) (their Jacobians move them off
  // the potentials' minima: ⟨r⟩ ≈ 5.04 Å, ⟨cos θ⟩ ≈ −0.47)
  double zr = 0, mr = 0, zc = 0, mc = 0;
  for (double x = 0.001; x < 20; x += 0.001) { const double w = x * x * std::exp(-3.0 * (x - 5) * (x - 5) / kT); zr += w, mr += w * x; }
  for (double t = 0.05; t < 180; t += 0.1) { const double w = std::sin(t * kPi / 180) * std::exp(-0.002 * (t - 120) * (t - 120) / kT); zc += w, mc += w * std::cos(t * kPi / 180); }
  EXPECT_NEAR(rs / nb, mr / zr, 0.01);
  const double b2 = r2 / nb, c = cs / na;
  EXPECT_NEAR(c, mc / zc, 0.01);
  // freely rotating chain: ⟨R²(n)⟩/n → ⟨b²⟩ (1 − ⟨cos θ⟩)/(1 + ⟨cos θ⟩) for long n
  const double cinf = b2 * (1 - c) / (1 + c);
  const double at = r.internal[400];
  EXPECT_NEAR(at / cinf, 1.0, 0.1) << "⟨R²(n)⟩/n at n = 401: " << at << " vs " << cinf;
  EXPECT_FALSE(r.notes.empty());
}

TEST(CgBuild, SequencesLengthsAndFiles) {
  CgBuildOptions o;
  o.units = {{"AB", {"A", "B"}}, {"AC", {"A", "C"}}};
  CgBondedResult t = tables();
  t.bonds.push_back(table(0, "A-C", 0, 20, 0.02, [](double x) { return 3.0 * (x - 5.5) * (x - 5.5); }));
  for (const char* k : {"A-C-A", "C-A-C", "B-A-C"}) t.angles.push_back(table(1, k, 0, 180, 1, [](double a) { return 0.002 * (a - 120) * (a - 120); }));
  o.dp = 200, o.chains = 20, o.density = 1.1;
  // Bernoulli 56 / 44
  o.composition = {{"AB", 0.56}, {"AC", 0.44}};
  auto r = build_cg_polymer(o, t, kMass);
  const double f = double(r.units_drawn["AB"]) / 4000.0;
  EXPECT_NEAR(f, 0.56, 0.03);
  // Markov with long blocks: few AB–AC junctions
  o.sequence = "markov";
  o.markov = {{"AB", {{"AB", 0.95}, {"AC", 0.05}}}, {"AC", {{"AB", 0.05}, {"AC", 0.95}}}};
  r = build_cg_polymer(o, t, kMass);
  int junctions = 0;
  for (const auto& ch : r.topology.chains)
    for (size_t i = 3; i < ch.size(); i += 2)
      junctions += r.beads.atoms[size_t(ch[i])].name != r.beads.atoms[size_t(ch[i - 2])].name;
  EXPECT_LT(junctions, int(20 * 199 * 0.1)) << "blocky";
  // gradient: the first unit early, the second late
  o.sequence = "gradient";
  r = build_cg_polymer(o, t, kMass);
  int early = 0, late = 0;
  for (const auto& ch : r.topology.chains) {
    early += r.beads.atoms[size_t(ch[1])].name == "B";
    late += r.beads.atoms[size_t(ch[ch.size() - 1])].name == "B";
  }
  EXPECT_GT(early, late);
  // Schulz–Zimm lengths with Đ = 2
  o.sequence = "bernoulli";
  o.lengths = "schulz-zimm", o.pdi = 2.0, o.chains = 400, o.dp = 50;
  r = build_cg_polymer(o, t, kMass);
  double sn = 0, sw = 0;
  for (int n : r.chain_dp) sn += n, sw += double(n) * n;
  EXPECT_NEAR((sw / sn) / (sn / 400), 2.0, 0.3);
  // a small box: the warning
  o.lengths = "monodisperse", o.chains = 2, o.dp = 400;
  r = build_cg_polymer(o, t, kMass);
  EXPECT_TRUE(std::any_of(r.notes.begin(), r.notes.end(), [](const std::string& n) { return n.find("box is smaller") != std::string::npos; }));
  // files: the map reads back as a topology, the data names its types
  const auto dir = std::filesystem::temp_directory_path() / "caps_cg_build_test";
  std::filesystem::create_directories(dir);
  const CgTypes types{{"A", "B", "C"}, {"A-B", "A-C"}, {"A-B-A", "A-C-A", "B-A-B", "B-A-C", "C-A-C"}, {}};
  o.chains = 3, o.dp = 10;
  r = build_cg_polymer(o, t, kMass);
  std::vector<std::string> keys;
  for (const auto& q : r.topology.dihedrals) keys.push_back(cg_key({r.topology.bead_kind[size_t(q[0])], r.topology.bead_kind[size_t(q[1])], r.topology.bead_kind[size_t(q[2])], r.topology.bead_kind[size_t(q[3])]}));
  CgTypes with = types;
  std::sort(keys.begin(), keys.end());
  keys.erase(std::unique(keys.begin(), keys.end()), keys.end());
  with.dihedrals = keys;
  const auto files = write_cg_build(r, with, (dir / "melt").string());
  ASSERT_EQ(files.size(), 2u);
  std::ifstream mf(files[1]);
  std::string js((std::istreambuf_iterator<char>(mf)), std::istreambuf_iterator<char>());
  const CgTopology back = cg_topology_from_map(Json::parse(js));
  EXPECT_EQ(back.beads(), r.beads.atoms.size());
  EXPECT_EQ(back.angles.size(), r.topology.angles.size());
  const std::string deck = cg_equil_deck({});
  EXPECT_NE(deck.find("fix push all adapt 1 pair soft a * * v_A"), std::string::npos);
  EXPECT_NE(deck.find("include ${PAIR}/pair.in"), std::string::npos);
  std::filesystem::remove_all(dir);
}
