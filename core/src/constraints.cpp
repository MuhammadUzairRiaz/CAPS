// CAPS bond constraints: SHAKE for positions, RATTLE for velocities (see constraints.hpp).
#include "caps/constraints.hpp"

#include <algorithm>
#include <cmath>
#include <tuple>
#include <cstdio>
#include <map>

#include "caps/elements.hpp"

namespace caps {

namespace {
constexpr double kAcc = 4.184e-4;   // (kcal/mol/Å)/(g/mol) → Å/fs²
}

ConstraintMode constraints_from_string(const std::string& s) {
  if (s == "none" || s.empty()) return ConstraintMode::None;
  if (s == "h-bonds" || s == "hbonds" || s == "h") return ConstraintMode::HBonds;
  if (s == "all-bonds" || s == "bonds" || s == "all") return ConstraintMode::AllBonds;
  throw std::invalid_argument("unknown constraints '" + s + "' (none, h-bonds, all-bonds)");
}

const char* to_string(ConstraintMode m) {
  switch (m) {
    case ConstraintMode::None: return "none";
    case ConstraintMode::HBonds: return "bonds to hydrogen, rigid water";
    case ConstraintMode::AllBonds: return "all bonds";
  }
  return "?";
}

ConstraintAlgorithm constraint_algorithm_from_string(const std::string& s) {
  if (s == "shake" || s == "rattle" || s.empty()) return ConstraintAlgorithm::Shake;
  if (s == "lincs") return ConstraintAlgorithm::Lincs;
  throw std::invalid_argument("unknown constraint solver '" + s + "' (shake, lincs)");
}

ConstraintSet make_constraints(const System& s, const ForceField& ff, ConstraintMode mode, const std::vector<char>& held) {
  ConstraintSet out;
  if (mode == ConstraintMode::None) return out;
  const size_t n = s.atoms.size();
  auto key = [](uint32_t a, uint32_t b) { return a < b ? std::make_pair(a, b) : std::make_pair(b, a); };
  // equilibrium lengths from the force field's bond terms
  std::map<std::pair<uint32_t, uint32_t>, double> r0;
  for (const auto& b : ff.bonds) r0[key(b.i, b.j)] = b.r0;
  for (const auto& b : ff.bonds2) r0[key(b.i, b.j)] = b.r0;
  for (const auto& b : ff.bonds_x) {
    if (b.form == 1) r0[key(b.i, b.j)] = b.c;        // Morse
    else if (b.form == 0) r0[key(b.i, b.j)] = b.b;   // GROMOS quartic
  }
  std::map<std::tuple<uint32_t, uint32_t, uint32_t>, double> th0;   // (end, vertex, end), ends ordered
  auto akey = [](uint32_t a, uint32_t j, uint32_t b) { return a < b ? std::make_tuple(a, j, b) : std::make_tuple(b, j, a); };
  for (const auto& a : ff.angles) th0[akey(a.i, a.j, a.k)] = a.theta0;
  for (const auto& a : ff.angles2) th0[akey(a.i, a.j, a.k)] = a.theta0;
  for (const auto& a : ff.angles_x) th0[akey(a.i, a.j, a.k)] = a.b;

  std::vector<char> skip(n, 0);
  for (size_t i = 0; i < n && i < held.size(); ++i) skip[i] = held[i];
  for (const auto& vs : ff.vsites) skip[vs.site] = 1;
  const auto nb = s.neighbours();
  auto length = [&](uint32_t i, uint32_t j, size_t& own) {
    auto it = r0.find(key(i, j));
    if (it != r0.end() && it->second > 0.3) return it->second;
    ++own;
    return norm(s.cell.valid() ? s.cell.minimum_image(s.atoms[j].pos - s.atoms[i].pos) : s.atoms[j].pos - s.atoms[i].pos);
  };
  size_t own = 0;
  std::vector<char> water(n, 0);
  if (mode == ConstraintMode::HBonds)
    for (uint32_t o = 0; o < n; ++o)
      if (s.atoms[o].element == 8 && nb[o].size() == 2 && s.atoms[nb[o][0]].element == 1 && s.atoms[nb[o][1]].element == 1 && !skip[o] &&
          !skip[nb[o][0]] && !skip[nb[o][1]])
        water[o] = 1;
  for (const auto& b : s.bonds) {
    if (skip[b.i] || skip[b.j] || s.atoms[b.i].element == 0 || s.atoms[b.j].element == 0) continue;
    const bool h = s.atoms[b.i].element == 1 || s.atoms[b.j].element == 1;
    if (mode == ConstraintMode::HBonds && !h) continue;
    out.c.push_back({b.i, b.j, length(b.i, b.j, own)});
    ++out.bonds;
  }
  for (uint32_t o = 0; o < n; ++o) {
    if (!water[o]) continue;
    const uint32_t a = nb[o][0], b = nb[o][1];
    const double da = length(o, a, own), db = length(o, b, own);
    double th;
    auto it = th0.find(akey(a, o, b));
    if (it != th0.end()) th = it->second;
    else {
      const Vec3 u = s.atoms[a].pos - s.atoms[o].pos, v = s.atoms[b].pos - s.atoms[o].pos;
      th = std::acos(std::clamp(dot(u, v) / (norm(u) * norm(v)), -1.0, 1.0));
    }
    out.c.push_back({a, b, std::sqrt(da * da + db * db - 2 * da * db * std::cos(th))});
    ++out.waters;
  }
  char buf[200];
  std::snprintf(buf, sizeof buf, "constraints: %zu bond%s%s held at fixed length", out.bonds, out.bonds == 1 ? "" : "s",
                out.waters ? (", " + std::to_string(out.waters) + " rigid water" + (out.waters == 1 ? "" : "s")).c_str() : "");
  out.notes.push_back(buf);
  if (own) out.notes.push_back(std::to_string(own) + " constrained length(s) taken from the structure (no bond length in the force field)");
  return out;
}

ConstraintSolver::ConstraintSolver(ConstraintSet set, const std::vector<double>& m, const Cell& cell, double tol, int max_iter, ConstraintAlgorithm algorithm)
    : set_(std::move(set)), inv_m_(m.size()), tol_(tol), max_iter_(max_iter), alg_(algorithm) {
  (void)cell;
  for (size_t i = 0; i < m.size(); ++i) inv_m_[i] = m[i] > 0 ? 1.0 / m[i] : 0.0;
  ref_.resize(set_.c.size());
  if (alg_ == ConstraintAlgorithm::Lincs) {
    std::vector<std::vector<uint32_t>> on(m.size());
    for (uint32_t k = 0; k < set_.c.size(); ++k) on[set_.c[k].i].push_back(k), on[set_.c[k].j].push_back(k);
    auto S = [&](uint32_t k) { const auto& c = set_.c[k]; return 1.0 / std::sqrt(inv_m_[c.i] + inv_m_[c.j]); };
    for (size_t a = 0; a < on.size(); ++a)
      for (uint32_t k : on[a])
        for (uint32_t l : on[a]) {
          if (k == l) continue;
          const double sk = set_.c[k].i == a ? 1.0 : -1.0, sl = set_.c[l].i == a ? 1.0 : -1.0;
          couple_.push_back({k, l, -S(k) * S(l) * sk * sl * inv_m_[a]});
        }
  }
}

// LINCS on the positions (the reference vectors give the old directions B)
void ConstraintSolver::lincs(std::vector<double>& x, std::vector<double>* v, double h) {
  const size_t n = set_.c.size();
  std::vector<Vec3> B(n);
  std::vector<double> S(n), a_kl(couple_.size());
  for (size_t k = 0; k < n; ++k) {
    const double l = norm(ref_[k]);
    if (l < 1e-12) throw FieldError("two constrained atoms sit on each other");
    B[k] = ref_[k] * (1 / l);
    S[k] = 1.0 / std::sqrt(inv_m_[set_.c[k].i] + inv_m_[set_.c[k].j]);
  }
  for (size_t q = 0; q < couple_.size(); ++q) a_kl[q] = couple_[q].coef * dot(B[couple_[q].k], B[couple_[q].l]);
  auto sep = [&](size_t k) {
    const auto& c = set_.c[k];
    const Vec3& r = ref_[k];
    return Vec3{r[0] + (x[3 * c.i] - x0_[3 * c.i]) - (x[3 * c.j] - x0_[3 * c.j]), r[1] + (x[3 * c.i + 1] - x0_[3 * c.i + 1]) - (x[3 * c.j + 1] - x0_[3 * c.j + 1]),
                r[2] + (x[3 * c.i + 2] - x0_[3 * c.i + 2]) - (x[3 * c.j + 2] - x0_[3 * c.j + 2])};
  };
  // (I − A)⁻¹ rhs ≈ Σ Aⁿ rhs, n ≤ 4 (Hess 1997's expansion order)
  std::vector<double> sol(n), t(n), t2(n);
  auto solve_apply = [&](const std::vector<double>& rhs) {
    sol = rhs;
    t = rhs;
    for (int order = 0; order < 4; ++order) {
      std::fill(t2.begin(), t2.end(), 0.0);
      for (size_t q = 0; q < couple_.size(); ++q) t2[couple_[q].k] += a_kl[q] * t[couple_[q].l];
      for (size_t k = 0; k < n; ++k) sol[k] += t2[k];
      t.swap(t2);
    }
    for (size_t k = 0; k < n; ++k) {
      const auto& c = set_.c[k];
      const double g = S[k] * sol[k];
      for (int d = 0; d < 3; ++d) {
        x[3 * c.i + d] -= inv_m_[c.i] * g * B[k][d];
        x[3 * c.j + d] += inv_m_[c.j] * g * B[k][d];
      }
    }
  };
  std::vector<double> before;
  if (v) before = x;
  std::vector<double> rhs(n);
  for (size_t k = 0; k < n; ++k) rhs[k] = S[k] * (dot(B[k], sep(k)) - set_.c[k].d);
  solve_apply(rhs);
  // the rotational lengthening: aim the projection along B at p = √(d² − ⊥²), ⊥ the separation's part across B (Hess's
  // √(2d² − l²) when the projection is already d), until every length holds
  bool done = false;
  for (int it = 0; it < max_iter_ && !done; ++it) {
    done = true;
    for (size_t k = 0; k < n; ++k) {
      const Vec3 s = sep(k);
      const double d = set_.c[k].d, l2 = dot(s, s), along = dot(B[k], s);
      if (std::abs(d * d - l2) > 2 * tol_ * d * d) done = false;
      if (along < 0) throw FieldError("a constrained bond turned by more than 90° in one step; lower the time step or relax first");
      const double perp2 = std::max(0.0, l2 - along * along);
      rhs[k] = S[k] * (along - std::sqrt(std::max(0.0, d * d - perp2)));
    }
    if (!done) solve_apply(rhs);
  }
  if (!done) throw FieldError("the bond constraints did not converge (LINCS); lower the time step or relax first");
  if (v && h > 0)
    for (size_t k = 0; k < x.size(); ++k) (*v)[k] += (x[k] - before[k]) / h;
}

void ConstraintSolver::reference(const std::vector<double>& x, const Cell& cell) {
  x0_ = x;
  for (size_t k = 0; k < set_.c.size(); ++k) {
    const auto& c = set_.c[k];
    const Vec3 d{x[3 * c.i] - x[3 * c.j], x[3 * c.i + 1] - x[3 * c.j + 1], x[3 * c.i + 2] - x[3 * c.j + 2]};
    ref_[k] = cell.valid() ? cell.minimum_image(d) : d;
  }
}

void ConstraintSolver::shake(std::vector<double>& x, std::vector<double>* v, double h) {
  if (set_.c.empty()) return;
  if (alg_ == ConstraintAlgorithm::Lincs) { lincs(x, v, h); return; }
  // separations follow from the reference vector and the two atoms' displacements, so no image is looked up per sweep
  auto sep = [&](size_t k) {
    const auto& c = set_.c[k];
    const Vec3& r = ref_[k];
    return Vec3{r[0] + (x[3 * c.i] - x0_[3 * c.i]) - (x[3 * c.j] - x0_[3 * c.j]), r[1] + (x[3 * c.i + 1] - x0_[3 * c.i + 1]) - (x[3 * c.j + 1] - x0_[3 * c.j + 1]),
                r[2] + (x[3 * c.i + 2] - x0_[3 * c.i + 2]) - (x[3 * c.j + 2] - x0_[3 * c.j + 2])};
  };
  std::vector<double> before;
  if (v) before = x;
  bool done = false;
  for (int it = 0; it < max_iter_ && !done; ++it) {
    done = true;
    for (size_t k = 0; k < set_.c.size(); ++k) {
      const auto& c = set_.c[k];
      const Vec3 s = sep(k);
      const double d2 = c.d * c.d, diff = d2 - dot(s, s);
      if (std::abs(diff) <= 2 * tol_ * d2) continue;
      done = false;
      const Vec3& r = ref_[k];
      const double sr = dot(s, r);
      if (sr < 1e-6 * d2) throw FieldError("a constrained bond turned by more than 90° in one step; lower the time step or relax first");
      const double g = diff / (2 * (inv_m_[c.i] + inv_m_[c.j]) * sr);
      for (int a = 0; a < 3; ++a) {
        x[3 * c.i + a] += g * inv_m_[c.i] * r[a];
        x[3 * c.j + a] -= g * inv_m_[c.j] * r[a];
      }
    }
  }
  if (!done) throw FieldError("the bond constraints did not converge (SHAKE); lower the time step or relax first");
  if (v && h > 0)
    for (size_t k = 0; k < x.size(); ++k) (*v)[k] += (x[k] - before[k]) / h;
}

void ConstraintSolver::rattle(const std::vector<double>& x, std::vector<double>& v, const Cell& cell, double h) {
  for (double& q : w_) q = 0;
  if (set_.c.empty()) return;
  std::vector<Vec3> r(set_.c.size());
  std::vector<double> lam(set_.c.size(), 0.0);
  for (size_t k = 0; k < set_.c.size(); ++k) {
    const auto& c = set_.c[k];
    const Vec3 d{x[3 * c.i] - x[3 * c.j], x[3 * c.i + 1] - x[3 * c.j + 1], x[3 * c.i + 2] - x[3 * c.j + 2]};
    r[k] = cell.valid() ? cell.minimum_image(d) : d;
  }
  bool done = false;
  for (int it = 0; it < max_iter_ && !done; ++it) {
    done = true;
    for (size_t k = 0; k < set_.c.size(); ++k) {
      const auto& c = set_.c[k];
      const Vec3 dv{v[3 * c.i] - v[3 * c.j], v[3 * c.i + 1] - v[3 * c.j + 1], v[3 * c.i + 2] - v[3 * c.j + 2]};
      const double d2 = c.d * c.d, rv = dot(r[k], dv);
      if (std::abs(rv) <= tol_ * d2) continue;   // tol in 1/fs: relative bond-length rate
      done = false;
      const double g = -rv / ((inv_m_[c.i] + inv_m_[c.j]) * d2);
      lam[k] += g;
      for (int a = 0; a < 3; ++a) {
        v[3 * c.i + a] += g * inv_m_[c.i] * r[k][a];
        v[3 * c.j + a] -= g * inv_m_[c.j] * r[k][a];
      }
    }
  }
  if (!done) throw FieldError("the bond constraints did not converge (RATTLE); lower the time step");
  if (h <= 0) return;
  // the impulse g r on atom i over a half kick of h/2 is the force 2 g r / (h kAcc); virial Σ r_ij · F_i
  static const int ia[6] = {0, 1, 2, 0, 0, 1}, ib[6] = {0, 1, 2, 1, 2, 2};
  for (size_t k = 0; k < set_.c.size(); ++k) {
    const double s = 2 * lam[k] / (h * kAcc);
    for (int q = 0; q < 6; ++q) w_[q] += s * r[k][ia[q]] * r[k][ib[q]];
  }
}

}  // namespace caps
