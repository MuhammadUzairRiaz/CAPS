// CAPS crystal habit: which faces a crystal shows and the shape it takes.
//
//  BFDH         faces ranked by their effective interplanar spacing (Donnay & Harker, Am. Mineral. 22, 446 (1937)): d_hkl of
//               the first order n(hkl) the crystal's symmetry does not extinguish (lattice centring, screw axes, glide
//               planes: an operation with h·R = h and h·t ∉ Z), growth rate ∝ 1/d_eff.
//  Attachment   E_att = E_latt − E_slice per growth unit (Hartman & Perdok, Acta Cryst. 8, 49 (1955); Hartman & Bennema,
//  energy       J. Cryst. Growth 49, 145 (1980)): the slice is one effective spacing of whole growth units (molecules, ions
//               with their shells) and E_att is the energy of a unit in the bulk less its energy in that slice alone (both
//               unrelaxed, so intramolecular energy cancels); of the cuts between units the one with the lowest |E_att| is
//               taken; growth rate ∝ |E_att|.
//  Surface      γ = (E_slab − N e_bulk) / 2A of an unrelaxed slab of several slices (at least the cut-off thick), J/m².
//  energy
//  Wulff shape  the intersection of the half-spaces n·x ≤ h over every face of the forms (Wulff, Z. Kristallogr. 34, 449
//               (1901)), with h = 1/d_eff (BFDH), |E_att| (growth) or γ (equilibrium); facets, areas and volume.
//  Particle     a crystal particle carved from the bulk in that shape: every growth unit whose centre lies inside.
// Faces are grouped into forms by the crystal's point group (from its space group; P1 when none is found). A slice that
// carries a dipole along its normal (Tasker type 3) or a net charge, or that would cut a covalent network or chain, has
// no energy: the face stays in BFDH only.
#pragma once

#include <array>
#include <functional>
#include <limits>
#include <string>
#include <vector>

#include "caps/field.hpp"
#include "caps/json.hpp"
#include "caps/system.hpp"

namespace caps {

struct HabitOptions {
  int max_index = 3;               // |h|, |k|, |l| ≤ this
  double min_d = 0;                // Å: faces with a smaller effective spacing are left out
  int forms = 12;                  // the energies are computed for this many forms, largest d_eff first
  bool attachment = false;         // attachment energies (growth shape)
  bool surface = false;            // surface energies (equilibrium shape)
  std::function<ForceField(const System&)> field_for;   // the force field of any structure built from the crystal (energies)
  EnergyOptions energy;
  double vacuum = 0;               // Å above each slab (0: twice the cut-off plus 10 Å)
  double symmetry_tolerance = 0.1; // Å
  size_t max_atoms = 15000;
  std::function<bool(const std::string&, double)> progress;   // false stops
};

struct HabitFace {
  std::array<int, 3> hkl{0, 0, 1};                  // the form's representative
  std::vector<std::array<int, 3>> members;          // every face of the form
  double d = 0, d_eff = 0;                          // Å
  int order = 1;                                    // d_eff = d / order
  int rank = 0;                                     // BFDH rank, 1 = largest d_eff
  double e_att = std::numeric_limits<double>::quiet_NaN();     // kcal/mol per growth unit
  double e_slice = std::numeric_limits<double>::quiet_NaN();   // kcal/mol per growth unit
  double gamma = std::numeric_limits<double>::quiet_NaN();     // J/m²
  double cut = 0;                                   // the slice's start, as a fraction of d along the normal
  std::string note;                                 // why an energy is missing
  double area_bfdh = 0, area_growth = 0, area_equilibrium = 0;   // fraction of the shape's surface
};

struct WulffShape {
  std::vector<Vec3> vertices;
  std::vector<std::vector<int>> facets;             // vertex indices, counter-clockwise seen from outside
  std::vector<int> facet_plane;                     // per facet: the plane it lies in (index into the planes given)
  std::vector<double> facet_area;
  double area = 0, volume = 0;
  double aspect = 0;                                // largest / smallest width across the facets' normals
  bool ok() const { return !facets.empty(); }
};

struct HabitResult {
  std::string space_group;
  int operations = 1;
  std::vector<HabitFace> faces;                     // BFDH order
  int units = 0;                                    // growth units per cell
  double e_latt = std::numeric_limits<double>::quiet_NaN();   // bulk energy per growth unit (with intramolecular), kcal/mol
  WulffShape bfdh, growth, equilibrium;
  std::vector<int> bfdh_face, growth_face, equilibrium_face;   // per shape plane: the face (form) it belongs to
  std::vector<std::string> notes;
};

HabitResult crystal_habit(const System& cell, const HabitOptions& o);

// The convex shape bounded by n_i·x ≤ h_i (unit normals n_i, h_i > 0). Empty when the planes do not bound a shape.
WulffShape wulff_shape(const std::vector<Vec3>& normals, const std::vector<double>& distances);

// A particle of the crystal in a shape (planes n_i·x ≤ h_i scaled so the farthest vertex is `radius` Å from the centre):
// every growth unit of the bulk whose centre of mass lies inside, kept whole. No cell.
System carve_particle(const System& cell, const std::vector<Vec3>& normals, const std::vector<double>& distances, double radius);

// The planes of a result's shape ("bfdh", "growth", "equilibrium"): unit Cartesian normals, distances, the face of each.
void habit_planes(const System& cell, const HabitResult& r, const std::string& model, std::vector<Vec3>& normals, std::vector<double>& distances,
                  std::vector<int>* face = nullptr);

// The result as JSON: space_group, operations, units, e_latt, faces [{hkl, faces, d, d_eff, order, rank, e_att, e_slice,
// gamma, note, area_bfdh, area_growth, area_equilibrium}], shapes {bfdh, growth, equilibrium: {vertices, facets,
// facet_face, facet_area, area, volume, aspect}}, notes (null for a missing number).
Json habit_json(const HabitResult& r);

}  // namespace caps
