// CAPS's machine-learned potential through libtorch (see caps/mlpot.hpp). The compiled MACE model takes one dict of
// tensors: positions, the one-hot elements (node_attrs, in the model's own element order), the directed edges within its
// cutoff (both directions, no self-edges), their periodic shifts (zero here: one molecule), batch / ptr for a single
// graph, and the training labels and weights, which inference does not use but the compiled graph expects. The input
// was checked against MACE's own Python calculator: the same energy and forces to 1e-12 eV.
#include "caps/mlpot.hpp"

#include <stdexcept>

#ifdef CAPS_WITH_TORCH
#include <torch/script.h>

#include <cmath>
#include <map>
#include <mutex>
#endif

namespace caps::ml {

#ifndef CAPS_WITH_TORCH

bool available() { return false; }
struct Potential::Impl {};
Potential::Potential(const std::string&) {
  throw std::runtime_error("this CAPS was built without libtorch: configure with -DCAPS_WITH_TORCH=ON "
                           "-DTorch_DIR=<torch>/share/cmake/Torch to use the machine-learned potential");
}
Potential::~Potential() = default;
Result Potential::evaluate(const std::vector<int>&, const std::vector<std::array<double, 3>>&) const { return {}; }
double Potential::cutoff() const { return 0; }
const std::vector<int>& Potential::elements() const { static const std::vector<int> none; return none; }
bool Potential::covers(const std::vector<int>&) const { return false; }

#else

bool available() { return true; }

struct Potential::Impl {
  mutable torch::jit::Module module;
  mutable std::mutex lock;    // one evaluation at a time: a TorchScript module is not reentrant
  double r_max = 0;
  std::vector<int> z;         // the model's elements, in its node_attrs order
  std::map<int, int> index;   // atomic number -> column of node_attrs
};

Potential::Potential(const std::string& path) : p_(std::make_unique<Impl>()) {
  try {
    p_->module = torch::jit::load(path, torch::kCPU);
  } catch (const std::exception& e) {
    throw std::runtime_error("cannot load the machine-learned potential " + path + ": " + e.what());
  }
  p_->module.eval();
  if (!p_->module.hasattr("r_max") || !p_->module.hasattr("atomic_numbers"))
    throw std::runtime_error(path + " is not a compiled MACE model (no r_max / atomic_numbers)");
  p_->r_max = p_->module.attr("r_max").toTensor().item<double>();
  auto an = p_->module.attr("atomic_numbers").toTensor().to(torch::kLong);
  for (int64_t k = 0; k < an.size(0); ++k) {
    p_->z.push_back(static_cast<int>(an[k].item<int64_t>()));
    p_->index[p_->z.back()] = static_cast<int>(k);
  }
}

Potential::~Potential() = default;
double Potential::cutoff() const { return p_->r_max; }
const std::vector<int>& Potential::elements() const { return p_->z; }

bool Potential::covers(const std::vector<int>& z) const {
  for (int e : z)
    if (!p_->index.count(e)) return false;
  return true;
}

Result Potential::evaluate(const std::vector<int>& z, const std::vector<std::array<double, 3>>& pos) const {
  const int64_t n = static_cast<int64_t>(z.size());
  if (n == 0 || pos.size() != z.size()) throw std::runtime_error("ml::Potential: no atoms, or positions and elements differ in number");
  if (!covers(z)) throw std::runtime_error("the machine-learned potential does not cover an element of this structure");
  const auto f64 = torch::TensorOptions().dtype(torch::kFloat64);
  const auto i64 = torch::TensorOptions().dtype(torch::kLong);
  const int64_t ne = static_cast<int64_t>(p_->z.size());
  auto positions = torch::empty({n, 3}, f64);
  auto node_attrs = torch::zeros({n, ne}, f64);
  auto pa = positions.accessor<double, 2>();
  auto na = node_attrs.accessor<double, 2>();
  for (int64_t i = 0; i < n; ++i) {
    for (int k = 0; k < 3; ++k) pa[i][k] = pos[i][k];
    na[i][p_->index.at(z[i])] = 1.0;
  }
  // directed edges within the cutoff (all pairs: fragments are small)
  std::vector<int64_t> snd, rcv;
  const double rc2 = p_->r_max * p_->r_max;
  for (int64_t i = 0; i < n; ++i)
    for (int64_t j = 0; j < n; ++j) {
      if (i == j) continue;
      double d2 = 0;
      for (int k = 0; k < 3; ++k) d2 += (pos[i][k] - pos[j][k]) * (pos[i][k] - pos[j][k]);
      if (d2 < rc2) {
        snd.push_back(i);
        rcv.push_back(j);
      }
    }
  const int64_t e = static_cast<int64_t>(snd.size());
  auto edge_index = torch::empty({2, e}, i64);
  auto ea = edge_index.accessor<int64_t, 2>();
  for (int64_t k = 0; k < e; ++k) {
    ea[0][k] = snd[k];
    ea[1][k] = rcv[k];
  }
  c10::Dict<std::string, torch::Tensor> d;
  d.insert("positions", positions);
  d.insert("node_attrs", node_attrs);
  d.insert("edge_index", edge_index);
  d.insert("shifts", torch::zeros({e, 3}, f64));
  d.insert("unit_shifts", torch::zeros({e, 3}, f64));
  d.insert("cell", torch::zeros({3, 3}, f64));
  d.insert("rcell", torch::zeros({3, 3}, f64));
  d.insert("volume", torch::zeros({1}, f64));
  d.insert("pbc", torch::zeros({1, 3}, torch::TensorOptions().dtype(torch::kBool)));
  d.insert("batch", torch::zeros({n}, i64));
  d.insert("ptr", torch::tensor({int64_t(0), n}, i64));
  d.insert("head", torch::zeros({1}, i64));
  // labels and weights: unused in inference, expected by the compiled graph
  d.insert("energy", torch::zeros({1}, f64));
  d.insert("forces", torch::zeros({n, 3}, f64));
  d.insert("stress", torch::zeros({1, 3, 3}, f64));
  d.insert("virials", torch::zeros({1, 3, 3}, f64));
  d.insert("dipole", torch::zeros({1, 3}, f64));
  d.insert("polarizability", torch::zeros({1, 3, 3}, f64));
  d.insert("charges", torch::zeros({n}, f64));
  d.insert("density_coefficients", torch::zeros({n, 1}, f64));
  d.insert("total_charge", torch::zeros({1}, f64));
  d.insert("total_spin", torch::ones({1}, f64));
  d.insert("elec_temp", torch::zeros({1}, f64));
  d.insert("fermi_level", torch::zeros({1}, f64));
  d.insert("external_field", torch::zeros({1, 3}, f64));
  for (const char* w : {"weight", "energy_weight", "forces_weight", "stress_weight", "virials_weight", "charges_weight"})
    d.insert(w, torch::ones({1}, f64));
  d.insert("dipole_weight", torch::ones({1, 3}, f64));
  d.insert("polarizability_weight", torch::ones({1, 3, 3}, f64));

  std::lock_guard<std::mutex> g(p_->lock);
  std::vector<torch::jit::IValue> args{d, false, true};   // training = false, compute_force = true
  auto out = p_->module.forward(args).toGenericDict();
  Result r;
  r.energy = out.at("energy").toTensor().item<double>();
  auto f = out.at("forces").toTensor().to(torch::kFloat64).contiguous();
  auto fa = f.accessor<double, 2>();
  r.forces.resize(n);
  for (int64_t i = 0; i < n; ++i) r.forces[i] = {fa[i][0], fa[i][1], fa[i][2]};
  return r;
}

#endif

}  // namespace caps::ml
