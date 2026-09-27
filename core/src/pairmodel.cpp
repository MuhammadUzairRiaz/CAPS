// CAPS pair model (see pairmodel.hpp).
#include "caps/pairmodel.hpp"

#include <cmath>
#include <limits>

namespace caps {

namespace {
constexpr double kCoulomb = 332.06371;   // kcal·Å/(mol·e²), as the Evaluator
constexpr double kPi = 3.14159265358979323846;
}  // namespace

double PairModel::lj(double r2, double e, double s) const {
  if (lj96_) {
    const double t = s / std::sqrt(r2), t3 = t * t * t;
    return e * (2 * t3 * t3 * t3 - 3 * t3 * t3);
  }
  const double q = s * s / r2, q3 = q * q * q;
  return 4 * e * (q3 * q3 - q3);
}

PairModel::PairModel(const ForceField& ff, double cutoff, bool coulomb, double dsf_alpha)
    : ff_(ff), rc_(cutoff), rc2_(cutoff * cutoff), alpha_(dsf_alpha), lj96_(ff.pair_form == "lj9-6"), nt_(ff.lj.size()) {
  erfc_c_ = std::erfc(alpha_ * rc_) / rc_;
  dsf_f_ = erfc_c_ / rc_ + 2 * alpha_ / std::sqrt(kPi) * std::exp(-alpha_ * alpha_ * rc2_) / rc_;
  qscale_ = kCoulomb / (ff.dielectric > 0 ? ff.dielectric : 1.0);
  coul_ = coulomb;
  bool any = false;
  for (double q : ff.charge) any |= q != 0;
  coul_ = coul_ && any;
  eps_.resize(nt_ * nt_), sig_.resize(nt_ * nt_), shift_.resize(nt_ * nt_);
  for (size_t a = 0; a < nt_; ++a)
    for (size_t b = 0; b < nt_; ++b) {
      const bool off = ff.excluded_type_pairs.count({int(std::min(a, b)), int(std::max(a, b))}) > 0;
      const PairType p = off ? PairType{0, 1} : mixed_pair(ff, int(a), int(b));
      eps_[a * nt_ + b] = p.eps, sig_[a * nt_ + b] = p.sigma;
      shift_[a * nt_ + b] = p.eps != 0 ? lj(rc2_, p.eps, p.sigma) : 0;
    }
  if (!ff.lj14_types.empty()) {
    eps14_.resize(nt_ * nt_), sig14_.resize(nt_ * nt_), shift14_.resize(nt_ * nt_);
    for (size_t a = 0; a < nt_; ++a)
      for (size_t b = 0; b < nt_; ++b) {
        const auto& A = ff.lj14_types[a];
        const auto& B = ff.lj14_types[b];
        const double sg = ff.mixing == "geometric" ? std::sqrt(A.sigma * B.sigma) : 0.5 * (A.sigma + B.sigma);
        eps14_[a * nt_ + b] = std::sqrt(A.eps * B.eps), sig14_[a * nt_ + b] = sg;
        shift14_[a * nt_ + b] = lj(rc2_, eps14_[a * nt_ + b], sg);
      }
  }
}

double PairModel::energy(uint32_t a, uint32_t b, double r2, bool p14) const {
  if (r2 >= rc2_) return 0.0;
  if (r2 < 1e-6) return std::numeric_limits<double>::infinity();
  const size_t tp = size_t(ff_.type_index[a]) * nt_ + size_t(ff_.type_index[b]);
  const double fl = p14 ? ff_.lj14 : 1.0, fq = p14 ? ff_.coul14 : 1.0;
  double e = 0;
  if (p14 && !eps14_.empty()) e = fl * (lj(r2, eps14_[tp], sig14_[tp]) - shift14_[tp]);
  else if (eps_[tp] != 0) e = fl * (lj(r2, eps_[tp], sig_[tp]) - shift_[tp]);
  if (coul_ && fq != 0) {
    const double qq = ff_.charge[a] * ff_.charge[b];
    if (qq != 0) {
      const double r = std::sqrt(r2);
      e += fq * qscale_ * qq * (std::erfc(alpha_ * r) / r - erfc_c_ + dsf_f_ * (r - rc_));
    }
  }
  return e;
}

}  // namespace caps
