#pragma once
// CAPS's machine-learned potential (data/ml/caps-mace-v1.model, a compiled MACE model) evaluated in the core through
// libtorch: energies and forces of a molecule or fragment, with no Python. Built only with -DCAPS_WITH_TORCH=ON (see
// the top-level CMakeLists); otherwise available() is false and Potential's constructor says how to build it.
#include <array>
#include <memory>
#include <string>
#include <vector>

namespace caps::ml {

struct Result {
  double energy = 0;                          // eV
  std::vector<std::array<double, 3>> forces;  // eV/Å, one per atom
};

// whether this build has libtorch
bool available();

// libtorch's own intra-op threads (call once, before the first evaluation; 1 when evaluations run in parallel threads)
void set_threads(int n);

class Potential {
 public:
  explicit Potential(const std::string& model_path);   // throws std::runtime_error (no libtorch, no file, not a model)
  ~Potential();
  Potential(const Potential&) = delete;
  Potential& operator=(const Potential&) = delete;
  // a non-periodic molecule: atomic numbers and positions (Å). Throws when an element is outside the model's.
  Result evaluate(const std::vector<int>& z, const std::vector<std::array<double, 3>>& pos) const;
  double cutoff() const;                       // Å
  const std::vector<int>& elements() const;    // the atomic numbers the model was trained on
  bool covers(const std::vector<int>& z) const;
 private:
  struct Impl;
  std::unique_ptr<Impl> p_;
};

// An internal coordinate held at a value while the rest relaxes: a bond (i, j), an angle (i, j, k, j the apex) or a
// dihedral (i, j, k, l), in atom indices (0-based); n is 2, 3 or 4.
struct Hold {
  int n = 0;
  std::array<int, 4> idx{0, 0, 0, 0};
  double target = 0;   // Å or radians; set from the starting geometry when hold_start is used
};

struct RelaxResult {
  double energy = 0;                          // eV, the potential's (the restraints' energy left out)
  std::vector<std::array<double, 3>> pos;     // Å
  std::vector<double> held;                   // each Hold's value at the end (Å or radians)
  bool converged = false;
  int steps = 0;
};

// the value of a Hold's coordinate in a geometry (Å or radians; a dihedral in (-pi, pi])
double hold_value(const Hold& h, const std::vector<std::array<double, 3>>& pos);

// L-BFGS minimisation on the potential (ASE's settings) with each Hold kept exactly at its target (the force projected
// off the held coordinates, the geometry brought back onto them after each step, as ASE's FixInternals does); stops when
// no atom's projected force exceeds fmax eV/Å, or after max_steps. k is unused (kept for the signature).
RelaxResult relax(const Potential& pot, const std::vector<int>& z, std::vector<std::array<double, 3>> pos,
                  const std::vector<Hold>& holds, double fmax = 0.01, int max_steps = 400, double k = 100.0);

}  // namespace caps::ml
