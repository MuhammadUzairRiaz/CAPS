// CAPS crystals and surfaces: CIF files, slabs cleaved along (hkl), and the surfaces fibre–polymer interfaces start from.
//
//  read_cif / parse_cif  CIF 1.1: cell lengths and angles, the symmetry operations (_symmetry_equiv_pos_as_xyz or
//                        _space_group_symop_operation_xyz; P1 when absent), atom sites (type symbol or label, fractional
//                        coordinates, occupancy ≥ 0.5 kept). Sites are expanded by symmetry into the full cell, duplicates
//                        within 0.01 Å merged, and bonds perceived from covalent radii across the periodic cell.
//  slab_terminations     The distinct cuts of a (hkl) plane: every gap between atomic planes along the surface normal,
//                        with the atoms of the planes either side and the bonds each cut breaks per nm².
//  cleave                A slab of whole layers (stoichiometric by construction) cut at the chosen termination: the
//                        surface mesh from the two shortest in-plane lattice vectors (Lagrange-reduced), x along the first,
//                        z along the surface normal; optionally made orthogonal (a rectangular supercell, sheared to 90°
//                        when the nearest one is within the strain limit), repeated na × nb, with vacuum above and below,
//                        and with the dangling bonds passivated (O–H on oxygens, M–OH on cations bonded to oxygen, H on
//                        others) up to the bulk coordination of each element.
#pragma once
#include <functional>
#include <stdexcept>
#include <string>
#include <vector>

#include "caps/system.hpp"

namespace caps {

struct CrystalError : std::runtime_error {
  using std::runtime_error::runtime_error;
};

System parse_cif(const std::string& text, const std::string& name = "");
System read_cif(const std::string& path);

struct Termination {
  double cut = 0;            // cut height within one layer, 0 ≤ cut < d (Å above the lowest plane)
  double gap = 0;            // distance between the planes either side of the cut, Å
  int bonds_cut = 0;         // bulk bonds the cut breaks, per surface cell
  double bonds_per_nm2 = 0;
  std::string top, bottom;   // formula of the atomic plane at the top surface and at the bottom surface
  std::string label;         // "O-terminated · 1.18 Å gap · 6.5 bonds/nm²"
};

struct SlabOptions {
  int h = 0, k = 0, l = 1;
  int layers = 3;            // repeats of the plane spacing d along the normal
  int termination = 0;       // index into slab_terminations (fewest bonds cut first)
  double vacuum = 15.0;      // Å, split above and below the slab
  bool orthogonal = true;    // rectangular surface cell (needed by the film grower)
  double max_strain = 0.02;  // shear allowed to make the surface cell rectangular
  int na = 1, nb = 1;        // surface supercell
  bool passivate = false;    // O–H / M–OH / H on dangling bonds (hydroxylated oxide surfaces)
};

struct SlabReport {
  int h = 0, k = 0, l = 1;
  double d = 0;                          // plane spacing, Å
  double a = 0, b = 0, gamma = 90;       // surface cell, Å and degrees
  double thickness = 0;                  // slab thickness (lowest to highest atom), Å
  double strain = 0;                     // shear applied to make the cell rectangular
  int added_h = 0, added_oh = 0;
  std::vector<Termination> terminations;
  std::vector<std::string> notes;
};

// The plane spacing d of (hkl) (reduced by their common divisor) and its terminations, fewest bonds cut first.
std::vector<Termination> slab_terminations(const System& bulk, int h, int k, int l, double* d = nullptr);
System cleave(const System& bulk, const SlabOptions& o, SlabReport* rep = nullptr);

// Adds O–H to under-coordinated oxygens (S, Se), M–OH to cations bonded to oxygen in the bulk and H to under-coordinated
// covalent network atoms (B, C, N, Si, P, Ge) of `s`, up to each element's coordination in `bulk`, pointing along
// outward_of(position) where it can. Metals keep bare surfaces. Returns (H added, OH added).
std::pair<int, int> passivate_surface(System& s, const System& bulk, const std::function<Vec3(const Vec3&)>& outward_of);

// The bonds of a periodic structure from covalent radii, counting every periodic image (no metal–metal bonds in
// ionic crystals; nearest neighbours in metals).
std::vector<Bond> crystal_bonds(const System& s);

// Hill-order formula of a set of atoms ("SiO2", "C6H12").
std::string formula_of(const System& s, const std::vector<size_t>& atoms);

}  // namespace caps
