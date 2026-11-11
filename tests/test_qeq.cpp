#include <gtest/gtest.h>

#include <cmath>
#include <map>
#include <numeric>

#include "caps/charges.hpp"
#include "caps/crystal.hpp"
#include "caps/molecule.hpp"
#include "caps/qeq.hpp"

using namespace caps;

namespace {
const std::string kCrystals = std::string(CAPS_SOURCE_DIR) + "/data/crystals/";
}

TEST(QEq, MoleculesAreNeutralWithSensibleSigns) {
  for (const std::string smi : {"O", "CO", "CC(=O)OC", "C[Si](C)(C)O"}) {
    const BuildResult b = build_molecule(smi);
    QEqReport r;
    const auto q = qeq_charges(b.system, QEqOptions{}, &r);
    EXPECT_NEAR(std::accumulate(q.begin(), q.end(), 0.0), 0.0, 1e-9) << smi;
    EXPECT_LT(r.residual, 1e-6);
    for (size_t i = 0; i < q.size(); ++i) {
      if (b.system.atoms[i].element == 8) EXPECT_LT(q[i], -0.2) << smi;    // oxygen draws electrons
      if (b.system.atoms[i].element == 14) EXPECT_GT(q[i], 0.3) << smi;   // silicon gives them
    }
  }
  const BuildResult w = build_molecule("O");
  const auto qw = qeq_charges(w.system);
  EXPECT_GT(qw[0], -1.2);   // water O between −0.5 and −1.2 e (QEq gives about −0.8)
  EXPECT_LT(qw[0], -0.5);
}

TEST(QEq, PeriodicCrystals) {
  const System nacl = read_cif(kCrystals + "rock-salt.cif");
  const auto q = qeq_charges(nacl);
  for (size_t i = 0; i < q.size(); ++i) EXPECT_EQ(q[i] > 0, nacl.atoms[i].element == 11);
  QEqOptions o;
  o.total_charge = 0;
  // a bulk block of quartz large enough for the 10 Å cut-off (the unit cell is 4.9 Å)
  SlabOptions so;
  so.layers = 4;
  so.na = 5, so.nb = 3;
  so.vacuum = 0;
  const System quartz = cleave(read_cif(kCrystals + "alpha-quartz.cif"), so);
  QEqReport qr;
  const auto qq = qeq_charges(quartz, o, &qr);
  EXPECT_NEAR(qr.cutoff, 10.0, 1e-9);
  std::map<int, double> mean;
  std::map<int, int> count;
  for (size_t i = 0; i < qq.size(); ++i) mean[quartz.atoms[i].element] += qq[i], count[quartz.atoms[i].element]++;
  mean[14] /= count[14], mean[8] /= count[8];
  EXPECT_NEAR(mean[14], -2 * mean[8], 1e-6);   // SiO2 neutral: q(Si) = −2 q(O)
  EXPECT_GT(mean[14], 0.8);   // QEq silica: Si about +1.3 e
}

// A cell narrower than twice the cut-off keeps the whole cut-off: quartz's unit cell (narrowest width 4.25 Å) takes the
// Ewald sum by default (Si about +1.3 e, as the 3×3×3 supercell); without Ewald the tapered kernel runs over every image
// within 10 Å, so the unit cell and the supercell agree, and a molecule in a 12 Å box matches the molecule alone
TEST(QEq, NarrowCellsKeepTheCutoff) {
  const System x = read_cif(kCrystals + "alpha-quartz.cif");
  System big;
  big.cell = x.cell;
  big.cell.a = x.cell.a * 3.0, big.cell.b = x.cell.b * 3.0, big.cell.c = x.cell.c * 3.0;
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      for (int k = 0; k < 3; ++k)
        for (auto a : x.atoms) {
          a.pos = a.pos + x.cell.a * double(i) + x.cell.b * double(j) + x.cell.c * double(k);
          big.atoms.push_back(a);
        }
  QEqReport r;
  const auto q = qeq_charges(x, QEqOptions{}, &r);
  EXPECT_NEAR(r.cutoff, 10.0, 1e-9);
  EXPECT_GT(r.alpha, 0.0);   // Ewald on by default for a cell periodic along a, b and c
  for (size_t i = 0; i < q.size(); ++i)
    if (x.atoms[i].element == 14) EXPECT_NEAR(q[i], 1.26, 0.03);
  QEqOptions o;
  o.ewald = false;
  QEqReport r1;
  const auto q1 = qeq_charges(x, o, &r1);
  const auto q3 = qeq_charges(big, o);
  EXPECT_NEAR(r1.cutoff, 10.0, 1e-9);
  EXPECT_EQ(r1.alpha, 0.0);
  for (size_t i = 0; i < q3.size(); ++i) EXPECT_NEAR(q3[i], q1[i % x.atoms.size()], 1e-6);
  for (size_t i = 0; i < q1.size(); ++i)
    if (x.atoms[i].element == 14) EXPECT_GT(q1[i], 0.8);   // a 2.1 Å minimum-image cut-off gave +0.38 e

  System m = build_molecule("C[Si](C)(C)O[Si](C)(C)C").system;
  m.cell.a = {0, 0, 0}, m.cell.b = {0, 0, 0}, m.cell.c = {0, 0, 0};
  const auto alone = qeq_charges(m);
  Vec3 c{0, 0, 0};
  for (const auto& a : m.atoms) c = c + a.pos;
  c = c * (1.0 / double(m.atoms.size()));
  for (auto& a : m.atoms) a.pos = a.pos - c + Vec3{6, 6, 6};
  m.cell.origin = {0, 0, 0}, m.cell.a = {12, 0, 0}, m.cell.b = {0, 12, 0}, m.cell.c = {0, 0, 12};
  m.cell.periodic = {true, true, true};
  const auto boxed = qeq_charges(m, o);
  for (size_t i = 0; i < boxed.size(); ++i) EXPECT_NEAR(boxed[i], alone[i], 0.02);

  // the default Ewald stays off for a cell too large to sum it quickly (7×7×7 quartz: 3087 atoms × about 11000
  // k-vectors), where the tapered image sum stands in with a note; Ewald asked for explicitly is kept on
  System q7;
  q7.cell = x.cell;
  q7.cell.a = x.cell.a * 7.0, q7.cell.b = x.cell.b * 7.0, q7.cell.c = x.cell.c * 7.0;
  for (int i = 0; i < 7; ++i)
    for (int j = 0; j < 7; ++j)
      for (int k = 0; k < 7; ++k)
        for (auto a : x.atoms) {
          a.pos = a.pos + x.cell.a * double(i) + x.cell.b * double(j) + x.cell.c * double(k);
          q7.atoms.push_back(a);
        }
  QEqReport r7;
  const auto qd = qeq_charges(q7, QEqOptions{}, &r7);
  EXPECT_EQ(r7.alpha, 0.0);
  bool noted = false;
  for (const auto& note : r7.notes) noted = noted || note.find("Ewald left off by default") != std::string::npos;
  EXPECT_TRUE(noted);
  for (size_t i = 0; i < qd.size(); ++i) EXPECT_NEAR(qd[i], q1[i % x.atoms.size()], 1e-6);
  Json js = Json::parse("{\"ewald\": true}");
  EXPECT_FALSE(qeq_options_from_json(js).ewald_auto);
  EXPECT_TRUE(qeq_options_from_json(Json::parse("{}")).ewald_auto);
}

// Charge adjustments: methanol's three methyl hydrogens equivalent (one class), averaged to one charge; scaled × 0.8;
// neutralised to a chosen total, evenly or in proportion to |q|.
TEST(Charges, AverageScaleNeutralise) {
  BuildOptions b;
  b.forcefield = "uff";
  const System s = build_molecule("CO", b).system;
  const auto cls = equivalent_atoms(s);
  std::map<int, int> per;
  for (size_t i = 0; i < s.atoms.size(); ++i) per[cls[i]]++;
  int three = 0;
  for (const auto& [c, k] : per) three += k == 3;
  EXPECT_EQ(three, 1);   // the methyl hydrogens
  auto rep = compute_charges(s, "gasteiger");
  std::vector<double> q = rep.q;
  ChargeAdjust a;
  a.average = true;
  adjust_charges(s, q, a, 0);
  for (size_t i = 0; i < s.atoms.size(); ++i)
    for (size_t j = 0; j < s.atoms.size(); ++j)
      if (cls[i] == cls[j]) EXPECT_NEAR(q[i], q[j], 1e-12);
  EXPECT_NEAR(std::accumulate(q.begin(), q.end(), 0.0), std::accumulate(rep.q.begin(), rep.q.end(), 0.0), 1e-12);   // the total kept
  auto q2 = q;
  ChargeAdjust sc;
  sc.scale = 0.8;
  adjust_charges(s, q2, sc, 0);
  for (size_t i = 0; i < q.size(); ++i) EXPECT_NEAR(q2[i], 0.8 * q[i], 1e-12);
  for (const char* mode : {"even", "proportional"}) {
    auto q3 = q;
    ChargeAdjust nz;
    nz.neutralise = mode;
    nz.target = 0.25;
    adjust_charges(s, q3, nz, 0);
    EXPECT_NEAR(std::accumulate(q3.begin(), q3.end(), 0.0), 0.25, 1e-12);
  }
}

// Ewald charges must not depend on the real-space cut-off. The Ohno kernel's difference from 1/r falls as −a²/(2r³)
// with a per element pair, so summed over a 3D lattice it grew like ln r_c (quartz Si +3.0 e at 6 Å, +5.2 e at 20 Å);
// switched off over a fixed 5–7 Å it settles. The Slater default settles by itself, at Si +1.26 e, and the image sum
// without Ewald comes to the same charges
TEST(QEq, ChargesIndependentOfTheCutoff) {
  const System x = read_cif(kCrystals + "alpha-quartz.cif");
  auto si_charge = [&](const char* kernel, double rc, bool ewald) {
    QEqOptions o;
    o.kernel = kernel, o.cutoff = rc, o.ewald = ewald, o.ewald_auto = false, o.tolerance = 1e-10;
    const auto q = qeq_charges(x, o);
    for (size_t i = 0; i < q.size(); ++i)
      if (x.atoms[i].element == 14) return q[i];
    return 0.0;
  };
  EXPECT_EQ(QEqOptions{}.kernel, "slater");
  for (const char* k : {"slater", "ohno"}) {
    const double s10 = si_charge(k, 10, true);
    EXPECT_NEAR(si_charge(k, 14, true), s10, 0.01) << k;
    EXPECT_NEAR(si_charge(k, 20, true), s10, 0.01) << k;
    EXPECT_NEAR(si_charge(k, 20, false), s10, 0.2) << k;
  }
  EXPECT_NEAR(si_charge("slater", 10, true), 1.2553, 0.005);
}
