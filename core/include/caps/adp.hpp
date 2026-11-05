// CAPS atomic displacement parameters from a trajectory.
//
//  Tensors    U = ⟨Δr Δr^T⟩ about each atom's mean position over the frames (positions followed continuously across
//             the cell: each step by its minimum image; the system's mass-weighted centre drift removed frame by
//             frame), U_iso = tr U / 3 and B = 8π² U_iso (Trueblood et al., Acta Cryst. A52, 770 (1996)).
//  Sites      for a supercell of a crystal (`supercell` repeats of the unit cell along a, b, c), the atoms whose mean
//             positions fold onto one place of the unit cell form a site; its U is their mean.
//  Ellipsoids the semi-axes c·√λ_k along U's eigenvectors, c from the probability p of the trivariate normal inside:
//             P(χ²₃ ≤ c²) = p (c = 1.5382 for 50 %).
//  CIF        U_cif = N⁻¹ A⁻¹ U A⁻ᵀ N⁻¹ with A the cell vectors as columns and N = diag(a*, b*, c*) (Grosse-Kunstleve &
//             Adams, J. Appl. Cryst. 35, 477 (2002)), written as _atom_site_aniso_U_11 … _23; U_eq from U_cif by
//             U_eq = ⅓ Σ_ij U^ij a*_i a*_j (a_i · a_j).
#pragma once

#include <array>
#include <string>
#include <vector>

#include "caps/system.hpp"

namespace caps {

using Sym6 = std::array<double, 6>;   // xx yy zz xy xz yz

struct AdpOptions {
  std::array<int, 3> supercell{1, 1, 1};
  double probability = 0.5;            // ellipsoids
  bool remove_drift = true;
};

struct AtomAdp {
  Vec3 mean{0, 0, 0};                  // Å
  Sym6 U{0, 0, 0, 0, 0, 0};            // Å², Cartesian
  double u_iso = 0, b = 0;             // Å²
  Vec3 axes{0, 0, 0};                  // ellipsoid semi-axes at the probability, Å, largest first
  std::array<Vec3, 3> dirs{};          // their directions
  double anisotropy = 1;               // largest / smallest eigenvalue
};

struct AdpResult {
  std::vector<AtomAdp> atoms;          // per atom of the topology
  std::vector<AtomAdp> sites;          // per site of the unit cell (supercell given), else empty
  std::vector<std::vector<uint32_t>> site_atoms;
  Cell unit_cell;                      // the cell the sites live in
  int frames = 0;
  double scale = 0;                    // c for the probability
  std::vector<std::string> notes;
};

AdpResult displacement_tensors(const Trajectory& t, const std::vector<size_t>& frames, const AdpOptions& o = {});

// c with P(χ²₃ ≤ c²) = p
double ellipsoid_scale(double p);
// a tensor's derived numbers (U_iso, B, axes at scale c, directions, anisotropy) from U
void finish_adp(AtomAdp& a, double c);
// Cartesian ↔ CIF (U^ij on the normalised reciprocal axes) for a cell
Sym6 u_cif_from_cartesian(const Sym6& U, const Cell& cell);
Sym6 u_cartesian_from_cif(const Sym6& Ucif, const Cell& cell);
double u_eq_from_cif(const Sym6& Ucif, const Cell& cell);

// A CIF of the mean positions (sites when there are, else atoms) in `cell` with their U_ij. labels: one per row.
void write_cif_adp(const std::vector<int>& elements, const std::vector<AtomAdp>& rows, const Cell& cell, const std::string& path,
                   const std::string& title = "caps");

}  // namespace caps
