// CAPS k-space electrostatics: smooth particle-mesh Ewald and the plain Ewald sum (Fortran kernels in
// core/fortran/caps_kspace.f90). The real-space part, the self energy and the corrections for excluded pairs live in
// the force-field evaluator (field.cpp), which calls these for the reciprocal part.
#pragma once
#include <complex>
#include <vector>

#include "caps/system.hpp"

namespace caps {

constexpr double kCoulombConstant = 332.06371;   // kcal·Å/(mol·e²), as LAMMPS units real

// Ewald coefficient β (1/Å) with erfc(β rc) = rtol, as GROMACS's ewald-rtol.
double ewald_beta(double cutoff, double rtol);

struct PmeGrid {
  int k[3] = {0, 0, 0};   // grid points along a, b, c (products of 2, 3, 5 and 7)
  int order = 4;          // B-spline order (4: cubic, as GROMACS's default)
  double beta = 0;
};
// Grid for a cell: the smallest 2-3-5-7 sizes with spacing at most `spacing` Å along each edge.
PmeGrid pme_grid(const Cell& cell, double beta, double spacing = 1.2, int order = 4);

// Reciprocal-space energy (kcal/mol); forces are added to f (3N); vir receives Σ r·f (xx yy zz xy xz yz).
double pme_reciprocal(const std::vector<double>& x, const std::vector<double>& q, const Cell& cell, const PmeGrid& g, std::vector<double>& f,
                      double vir[6]);
double ewald_reciprocal(const std::vector<double>& x, const std::vector<double>& q, const Cell& cell, double beta, const int kmax[3],
                        std::vector<double>& f, double vir[6]);

// The Fortran FFT (for tests): in-place 3-D transform of data[i + k1 (j + k2 k)], X = Σ x exp(sign 2πi m·k/K).
void fft3d(std::vector<std::complex<double>>& data, int k1, int k2, int k3, int sign);
int fft_good_size(int n);

}  // namespace caps
