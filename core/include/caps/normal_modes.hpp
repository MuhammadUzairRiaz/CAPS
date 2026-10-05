// CAPS normal modes: the harmonic vibrations of a structure at a minimum of its force field.
//
//  Hessian      central differences of the forces, ±h on each coordinate of the moving atoms (the pair list frozen at
//               the minimum, so the second derivatives are smooth), symmetrised, mass-weighted: H_ij / √(m_i m_j).
//  Projection   the rigid translations (and, without a periodic cell, the rotations) projected out of the mass-weighted
//               Hessian, so they come back as exact zeros and are dropped (Miller, Handy & Adams 1980).
//  Frequencies  ν̃ = √λ / (2π c) with λ in kcal/mol/Å²/(g/mol): ν̃ = 108.59 √λ cm⁻¹; imaginary modes (λ < 0: not at a
//               minimum) are reported as negative wavenumbers.
//  Intensities  fixed partial charges: dμ/dQ_k = Σ_i q_i e_ik / √m_i, the IR intensity ∝ |dμ/dQ_k|² (relative, the
//               polarisation of the charges ignored).
//  Thermo       quantum harmonic oscillators at T: zero-point energy, vibrational energy, entropy and heat capacity.
#pragma once
#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "caps/field.hpp"
#include "caps/system.hpp"

namespace caps {

struct NormalModesOptions {
  std::shared_ptr<const ForceField> field;   // null: CAPS's default force field
  EnergyOptions energy;
  double step = 0.005;                       // Å, finite-difference displacement
  std::vector<char> moving;                  // per atom: 1 vibrates (empty: all); the others are held
  double temperature = 298.15;               // K, for the vibrational thermodynamics
  size_t max_atoms = 600;                    // the dense Hessian's limit (3N ≤ 1800: seconds)
  std::function<bool(double)> progress;      // fraction done → false cancels
};

struct NormalModesResult {
  std::vector<double> wavenumber;            // cm⁻¹, ascending; negative: imaginary
  std::vector<std::vector<double>> mode;     // Cartesian displacements (3N, every atom; held ones zero), unit length
  std::vector<double> ir;                    // relative IR intensity (largest 1), 0 without charges
  std::vector<double> reduced_mass;          // g/mol, 1/Σ|l_i|² of the mass-weighted unit mode l = q/√m (the quantum-chemistry
                                             // convention: a homonuclear diatomic gives the atomic mass)
  int projected = 0;                         // rigid motions removed (3 or 6, 5 for a linear molecule; 0 with held atoms)
  int imaginary = 0;                         // modes with ν̃ < −10 cm⁻¹
  double max_force = 0;                      // kcal/mol/Å at the structure: should be small (a minimum)
  double zpe = 0, e_vib = 0, s_vib = 0, cv_vib = 0;   // kcal/mol, kcal/mol (incl. ZPE), cal/mol/K, cal/mol/K
  std::string field;
  std::vector<std::string> notes;
};

NormalModesResult normal_modes(const System& s, const NormalModesOptions& o);

// Frames of mode k (index into result.mode) for the trajectory player: one period in `frames` steps, the largest atom
// displacement `amplitude` Å.
std::vector<std::vector<Vec3>> mode_frames(const System& s, const NormalModesResult& r, size_t k, double amplitude, int frames);

// Eigen-decomposition of a symmetric n × n matrix (row-major a, overwritten): eigenvalues ascending in w, eigenvectors
// as the columns of a (Householder tridiagonalisation, then the implicit QL algorithm with Wilkinson shifts).
void symmetric_eigen(std::vector<double>& a, size_t n, std::vector<double>& w);

}  // namespace caps
