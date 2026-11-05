#pragma once
// Cell relaxation with chosen lattice parameters held (monoclinic PE and iPP, a fibre axis kept, a layer spacing kept):
// the variables are the free ones of a, b, c, α, β, γ and the atoms in fractional coordinates; the cell keeps its
// orientation (h = R₀ h_std(a, b, c, α, β, γ), R₀ fixed at the start, h_std with a along x and b in the xy plane).
// The enthalpy H = E + P V is minimised by L-BFGS (Liu & Nocedal, Math. Program. 45, 503 (1989)) with its exact
// gradient: atoms ∂E/∂s = hᵀ ∂E/∂x; cell ∂H/∂h = (P V I − W) h⁻ᵀ at fixed fractional coordinates, W = Σ r ⊗ f the
// virial tensor (the cell derivative of Parrinello & Rahman, J. Appl. Phys. 52, 7182 (1981)), carried to each lattice
// parameter through ∂h/∂p. A cell narrower than twice the cut-off is evaluated as a supercell of copies that move
// together, so every neighbour within the cut-off is seen.
#include <array>
#include <functional>
#include <string>
#include <vector>

#include "caps/field.hpp"
#include "caps/system.hpp"

namespace caps {

// a, b, c (Å), α, β, γ (degrees)
std::array<double, 6> lattice_parameters(const Cell& c);
// The cell of those parameters in the standard orientation (a along x, b in the xy plane, c with positive z).
Cell standard_cell(const std::array<double, 6>& p);

struct CellRelaxOptions {
  std::array<bool, 6> hold{};       // a b c α β γ kept as they are
  double pressure = 0.0;            // atm, hydrostatic
  double ftol = 0.01;               // kcal/mol/Å: atoms; the cell variables are scaled to the same units (Å for angles: rad × V^⅓)
  int max_steps = 3000;
  bool move_atoms = true;
  std::vector<char> fixed;          // per atom: kept at its fractional position
  std::function<bool(int step, double enthalpy, double gmax)> progress;   // false cancels (throws std::runtime_error)
};
struct CellRelaxReport {
  std::array<double, 6> before{}, after{};   // a b c α β γ
  std::array<double, 6> stress{};            // GPa, Voigt xx yy zz yz xz xy, tensile positive, at the end
  std::array<double, 6> gradient{};          // ∂H/∂p at the end: kcal/mol/Å for lengths, kcal/mol/rad for angles
  double energy_before = 0, energy_after = 0, enthalpy_after = 0;   // kcal/mol (the cell given, copies divided out)
  double volume_before = 0, volume_after = 0;                        // Å³
  double pressure_after = 0;                                         // atm, −tr σ / 3
  std::array<int, 3> copies{1, 1, 1};
  int steps = 0;
  bool converged = false;
  std::vector<std::string> notes;
};
// field_for: the force field of any supercell of s (the copies). Positions, cell and nothing else of s change.
CellRelaxReport relax_cell(System& s, const std::function<ForceField(const System&)>& field_for, const EnergyOptions& eo,
                           const CellRelaxOptions& o = {});

}  // namespace caps
