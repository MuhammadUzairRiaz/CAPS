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

}  // namespace caps::ml
