// CAPS pair model: the non-bonded energy of one atom pair with a force field's own parameters, for Monte Carlo codes
// that move a few atoms at a time (CBMC regrowth, the adsorption locator, test-particle insertion).
//   van der Waals: Lennard-Jones 12-6, or class II 9-6, with the force field's mixing rule and explicit pairs, shifted to
//   zero at the cut-off (1-4 pairs scaled by lj14, or with the separate 1-4 parameters mixed by the same rule);
//   electrostatics: damped shifted force (Fennell & Gezelter 2006), 1-4 pairs scaled by coul14, over the dielectric.
// The same forms as the Evaluator with tail corrections off and DSF electrostatics.
#pragma once

#include <cstdint>
#include <vector>

#include "caps/field.hpp"

namespace caps {

class PairModel {
 public:
  PairModel(const ForceField& ff, double cutoff, bool coulomb = true, double dsf_alpha = 0.2);
  double cutoff() const { return rc_; }
  bool coulomb() const { return coul_; }
  // atoms a and b (the force field's atom indices) at squared distance r2 (Å²); p14: a 1-4 pair (scaled). Beyond the
  // cut-off: 0; closer than 10⁻³ Å: +∞.
  double energy(uint32_t a, uint32_t b, double r2, bool p14 = false) const;

 private:
  const ForceField& ff_;
  double rc_, rc2_, alpha_, erfc_c_, dsf_f_, qscale_;
  bool coul_, lj96_;
  size_t nt_;
  std::vector<double> eps_, sig_, shift_, eps14_, sig14_, shift14_;
  double lj(double r2, double e, double s) const;
};

}  // namespace caps
