// Optimal superposition of two conformations of the same atoms: Horn's closed form with unit quaternions (J. Opt. Soc.
// Am. A 4, 629 (1987)), weighted, and the per-atom shifts that remain.
#pragma once
#include <array>
#include <vector>

#include "caps/system.hpp"

namespace caps {

struct Superposition {
  std::array<std::array<double, 3>, 3> rot{{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}};   // moving → reference: x' = R (x − c_mov) + c_ref
  Vec3 centre_ref{0, 0, 0}, centre_mov{0, 0, 0};
  double rmsd = 0;   // over the fitted atoms (weighted), Å
  size_t fitted = 0;
  Vec3 apply(const Vec3& x) const;
};

// Fits mov onto ref over the atoms with a positive weight (w empty: every atom, weight 1). Throws std::invalid_argument
// when the sizes differ or fewer than three atoms are fitted (one or two fit trivially: then only the centres are
// matched and rot stays the identity).
Superposition superpose(const std::vector<Vec3>& ref, const std::vector<Vec3>& mov, const std::vector<double>& w = {});

// RMSD over the atoms with mask[i] after applying fit to mov (mask empty: every atom).
double rmsd_after(const Superposition& fit, const std::vector<Vec3>& ref, const std::vector<Vec3>& mov, const std::vector<char>& mask = {});

}  // namespace caps
