// Internal: PME's reciprocal part run on a thread pool (the evaluator's). See kspace.hpp.
#pragma once
#include <vector>

#include "caps/kspace.hpp"
#include "parallel.hpp"

namespace caps {

double pme_reciprocal(const std::vector<double>& x, const std::vector<double>& q, const Cell& cell, const PmeGrid& g, std::vector<double>& f,
                      double vir[6], ThreadPool& pool);

}  // namespace caps
