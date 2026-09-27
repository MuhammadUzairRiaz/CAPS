// CAPS bond constraints for dynamics: SHAKE (Ryckaert, Ciccotti & Berendsen, J. Comput. Phys. 23, 327 (1977)) for the
// positions and RATTLE (Andersen, J. Comput. Phys. 52, 24 (1983)) for the velocities, with velocity Verlet.
#pragma once
#include <cstdint>
#include <string>
#include <vector>

#include "caps/field.hpp"
#include "caps/system.hpp"

namespace caps {

enum class ConstraintMode { None, HBonds, AllBonds };
// How the positions are brought back onto the constraints: SHAKE (iterative, one constraint at a time) or LINCS (Hess,
// Bekker, Berendsen & Fraaije, J. Comput. Chem. 18, 1463 (1997): a matrix expansion of order 4 over coupled constraints,
// then corrections for the lengthening by rotation until every constraint holds to the tolerance). Both solve the same
// equations along the old bond directions; velocities are corrected by RATTLE in both.
enum class ConstraintAlgorithm { Shake, Lincs };

ConstraintMode constraints_from_string(const std::string& s);   // "none" | "h-bonds" | "all-bonds"
const char* to_string(ConstraintMode m);
ConstraintAlgorithm constraint_algorithm_from_string(const std::string& s);   // "shake" | "lincs"

struct DistanceConstraint {
  uint32_t i, j;
  double d;   // Å
};

struct ConstraintSet {
  std::vector<DistanceConstraint> c;
  size_t bonds = 0, waters = 0;   // bond constraints; rigid waters (their H···H distance)
  std::vector<std::string> notes;
};

// The constraints of a mode: HBonds holds every bond to a hydrogen at the force field's r0 and makes water rigid (an
// oxygen bonded to exactly two hydrogens: their H···H distance from r0 and θ0, as LAMMPS fix shake's "a" option);
// AllBonds every bond. Bonds without a length in the force field keep the length they have in s. Virtual sites and
// held atoms are never constrained.
ConstraintSet make_constraints(const System& s, const ForceField& ff, ConstraintMode mode, const std::vector<char>& held);

// The solver, for one run: positions x (3N, Å), velocities v (3N, Å/fs), masses m (g/mol).
class ConstraintSolver {
 public:
  ConstraintSolver(ConstraintSet set, const std::vector<double>& m, const Cell& cell, double tol = 1e-8, int max_iter = 1000,
                   ConstraintAlgorithm algorithm = ConstraintAlgorithm::Shake);

  size_t size() const { return set_.c.size(); }
  const ConstraintSet& set() const { return set_; }

  // Before a drift: the constrained vectors (minimum image) at the current positions.
  void reference(const std::vector<double>& x, const Cell& cell);
  // After a drift of length h: moves x back onto the constraints along the reference vectors (SHAKE) and, when v is
  // given, adds the same correction divided by h to the velocities. Throws FieldError when it does not converge.
  void shake(std::vector<double>& x, std::vector<double>* v, double h);
  // Removes the relative velocity along each constraint (RATTLE's second half). With h > 0 the impulse is read as a
  // force over a half kick of h/2 and its virial (kcal/mol) is kept for virial()/tensor().
  void rattle(const std::vector<double>& x, std::vector<double>& v, const Cell& cell, double h = 0);

  double virial() const { return w_[0] + w_[1] + w_[2]; }
  const double* tensor() const { return w_; }   // xx yy zz xy xz yz

 private:
  ConstraintSet set_;
  std::vector<double> inv_m_;
  std::vector<Vec3> ref_;
  std::vector<double> x0_;
  double tol_;
  int max_iter_;
  ConstraintAlgorithm alg_;
  // LINCS: for each pair of constraints sharing an atom, S_k S_l s_ka s_la / m_a (the coupling without the directions)
  struct Coupling { uint32_t k, l; double coef; };
  std::vector<Coupling> couple_;
  void lincs(std::vector<double>& x, std::vector<double>* v, double h);
  double w_[6] = {0, 0, 0, 0, 0, 0};
};

}  // namespace caps
