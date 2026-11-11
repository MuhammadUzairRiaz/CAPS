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
//                        others) up to the bulk coordination of each element. Molecular crystals: whole_molecules puts
//                        each molecule by its centroid, so the surfaces are made of whole molecules (no broken bonds).
#pragma once
#include <functional>
#include <stdexcept>
#include <string>
#include <array>
#include <cstdint>
#include <map>
#include <vector>

#include "caps/system.hpp"

namespace caps {

struct CrystalError : std::runtime_error {
  using std::runtime_error::runtime_error;
};

System parse_cif(const std::string& text, const std::string& name = "");
System read_cif(const std::string& path);

// Disorder in a CIF (Hall, Allen & Brown, Acta Cryst. A47, 655 (1991); the core dictionary's _atom_site_occupancy,
// _atom_site_disorder_assembly and _atom_site_disorder_group). Sites of a disorder group ("." and 0: none) belong to
// an assembly (one unnamed assembly when not given); partially occupied sites without a group that lie within 0.8 Å
// of each other are alternatives for one place (split positions, mixed elements on one site), the rest of the
// occupancy a vacancy.
//   major   each assembly keeps its group of the largest occupancy; each place keeps its occupant of the largest
//           occupancy when the place's total is at least 0.5 (the default; a lone site below 0.5 is left out)
//   group   each assembly keeps the named group (its major group when it has none of that name); places as major
//   random  a supercell filled at random to the occupancies: one group drawn for each copy of an assembly (each
//           symmetry operation and cell, so a disordered molecule stays whole), one occupant or a vacancy for each
//           place, from `seed`
struct CifDisorder {
  std::string mode = "major";              // major | group | random
  std::string group;                       // mode group: the group label kept
  std::array<int, 3> supercell{1, 1, 1};   // copies of the cell (every mode; random draws each copy)
  uint64_t seed = 1;
};
struct CifDisorderReport {
  int assemblies = 0, groups = 0;          // disorder assemblies and groups in the file
  int places = 0;                          // places with alternatives or a vacancy (per cell)
  std::map<int, double> expected;          // element → atoms the occupancies give (whole result)
  std::map<int, int> obtained;             // element → atoms written
  std::vector<std::string> notes;
};
System parse_cif(const std::string& text, const std::string& name, const CifDisorder& o, CifDisorderReport* report = nullptr);
System read_cif(const std::string& path, const CifDisorder& o, CifDisorderReport* report = nullptr);

// Displacement parameters a CIF gives, by site label: U_iso (Å², from _atom_site_U_iso_or_equiv or B/8π²) and the
// anisotropic U^ij (_atom_site_aniso_U_11 … _23 or B^ij/8π², on the normalised reciprocal axes: adp.hpp converts them).
struct CifDisplacements {
  std::map<std::string, double> u_iso;
  std::map<std::string, std::array<double, 6>> u_aniso;   // 11 22 33 12 13 23
};
CifDisplacements cif_displacements(const std::string& text);

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
  // molecular crystals (polymer crystals, cellulose, organic solids): each molecule goes into the layer of its centroid,
  // whole, instead of being cut at the atomic planes; a molecule endless along the normal is cut as the atoms are
  bool whole_molecules = false;
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
// covalent network atoms (B, C, N, Si, P, Ge) of `s`, on the free vertices of the atom's bulk polyhedron (from `bulk`)
// fitted to its surviving bonds, outward_of(position) breaking ties. An O site lost by two surface cations is put back
// once, bridging them. Metals keep bare surfaces. Returns (H added, OH added).
std::pair<int, int> passivate_surface(System& s, const System& bulk, const std::function<Vec3(const Vec3&)>& outward_of);

// The bonds of a periodic structure from covalent radii, counting every periodic image (no metal–metal bonds in
// ionic crystals; nearest neighbours in metals).
std::vector<Bond> crystal_bonds(const System& s);

// Hill-order formula of a set of atoms ("SiO2", "C6H12").
std::string formula_of(const System& s, const std::vector<size_t>& atoms);

}  // namespace caps
