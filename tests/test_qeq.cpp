#include <gtest/gtest.h>

#include <cmath>
#include <map>
#include <numeric>

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
