// CAPS dielectric, piezoelectric and optical response of a crystal from its force field's second derivatives (the lattice
// dynamics of the relaxed crystal at k = 0; Gale & Rohl, Mol. Simul. 29, 291 (2003)).
//
//  Force constants  K = ∂²E/∂x∂x of the cell's particles (cores and shells), the crystal's periodic copies moving
//                   together (central differences of the copies' mean forces); Ewald sums with conducting boundaries give
//                   the k → 0 transverse response.
//  Static           ε⁰_αβ = δ_αβ + (1/ε₀V) Σ_ij q_i q_j (K⁺)_iα,jβ over every particle (the uniform translation projected out).
//  High frequency   ε^∞ the same over the shells alone with the cores held (1 for rigid ions); refractive indices √ of its
//                   principal values.
//  Piezoelectric    proper constants of the relaxed ions, e_αj = (1/V) Σ_i q_i (K⁺Θ)_iα,j, Θ_iα,j = ∂F_iα/∂ε_j the forces of
//                   a homogeneous strain at fixed fractional coordinates (C/m²); d = e S with the relaxed-ion compliance
//                   (pC/N).
//  Frequency        the shells condensed out (K_eff = K_cc − K_cs K_ss⁻¹ K_sc, Z_eff = q_c − K_cs K_ss⁻¹ q_s), the Γ modes of
//                   K_eff/√(m m) and their oscillator strengths: ε(ω) = ε^∞ + Σ_m S_m ω_m² / (ω_m² − ω² − iγω), and the
//                   normal-incidence reflectivity R = |(√ε − 1)/(√ε + 1)|² along the principal axes of ε⁰ (γ a damping chosen by
//                   the user). ε(0) = ε⁰.
// Lyddane–Sachs–Teller (Phys. Rev. 59, 673 (1941)) for a diatomic cubic crystal: ε⁰/ε^∞ = ω_LO²/ω_TO².
#pragma once

#include <array>
#include <functional>
#include <string>
#include <vector>

#include "caps/field.hpp"
#include "caps/system.hpp"

namespace caps {

using Mat3x3 = std::array<std::array<double, 3>, 3>;

struct DielectricOptions {
  std::function<ForceField(const System&)> field_for;   // required
  EnergyOptions energy;
  std::vector<char> shell;           // per atom of the cell: a shell
  bool relax = true;                 // relax the atoms first (Newton steps on the force constants); false: the cores as
                                     // given, their shells relaxed
  bool relax_cell = true;            // and the cell (at zero stress) before them
  double relax_ftol = 1e-4;          // kcal/mol/Å, the largest gradient of the cell and elastic-constant relaxations
  int relax_iterations = 1000;       // and their step budget each
  double step = 0.002;               // Å, the force constants' displacement
  double strain = 1e-4;              // the internal-strain tensor's strain
  bool piezo = true;
  double damping = 10;               // cm⁻¹, γ of the oscillators
  double w_max = 0;                  // cm⁻¹, the spectrum's upper end (0: 1.2 × the highest mode)
  int points = 400;
  std::function<bool(const std::string&, double)> progress;
};

struct DielectricResult {
  Mat3x3 eps0{}, epsinf{};
  std::array<double, 3> eps0_principal{}, epsinf_principal{}, refractive{};
  std::array<std::array<double, 6>, 3> e{};   // C/m², Voigt j = xx yy zz yz xz xy
  std::array<std::array<double, 6>, 3> d{};   // pC/N
  bool piezo = false;
  std::array<std::array<double, 6>, 6> C{};   // GPa, relaxed ions (for d)
  std::vector<double> modes;                  // cm⁻¹, the Γ modes of the cores (the three translations dropped)
  std::vector<double> w, reflectivity[3], eps_re[3], eps_im[3];   // spectrum along the principal axes of ε⁰
  std::array<Vec3, 3> axes{};                 // principal axes of ε⁰
  std::array<double, 6> cell{};               // a b c α β γ of the (relaxed) cell
  Cell cell_vectors;
  std::vector<Vec3> positions;                // the cell's particles where the response was taken
  double max_force = 0;                       // kcal/mol/Å there
  double volume = 0;
  bool shells = false;
  std::vector<std::string> notes;
};

DielectricResult crystal_dielectric(const System& cell, const DielectricOptions& o);

}  // namespace caps
