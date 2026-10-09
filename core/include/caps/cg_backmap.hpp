// CAPS backmapping of chemistry-specific coarse-grained melts to all atoms (the coarse-graining workflow's last stage):
// fragment per bead, by the mapping of caps/cg_rules.hpp.
//
//   library      from a reference all-atom cell (a LAMMPS data file with its force field: types, charges, coefficients) and
//                its map.json: each bead class (kind and place in the chain: head, inner, tail — an inner B and a B with the
//                HO– end are different fragments) with up to `conformers` instances (their atoms in one order) and, for
//                every bond, angle, dihedral and improper of the reference chains, a template keyed by the classes of the
//                consecutive beads it spans and each atom's role in its bead — so the terms across cut bonds come back
//                with their own types
//   placement    each target bead gets a fragment of its class (an instance drawn at random), its centre of mass on the
//                bead, turned so that the directions to its two nearest chain neighbours match the reference instance's
//                (Horn's fit of the bead and its two neighbours); chains whose ends do not match the reference's are tried
//                reversed
//   topology     every template instantiated wherever its class sequence occurs: the cut bonds restored, angles and
//                dihedrals across them, charges and types per atom role
//   relaxation   a staged LAMMPS deck: bonded terms only (pairs off), a soft-core push-off, the full force field (the
//                reference input's pair, kspace and special_bonds lines), a short NPT
//   check        bond lengths and angles against the harmonic force field's r₀ and θ₀ (per type: mean and largest deviation)
#pragma once
#include <array>
#include <map>
#include <string>
#include <vector>

#include "caps/cg_bonded.hpp"
#include "caps/cg_rules.hpp"
#include "caps/system.hpp"

namespace caps {

// A LAMMPS data file of atom_style full, kept as it is: coefficient sections verbatim, terms with their type numbers.
struct LammpsFull {
  std::string title;
  Cell cell;
  std::vector<std::string> masses;                                  // "type mass  # label" lines
  std::vector<std::pair<std::string, std::vector<std::string>>> coeffs;   // ("Pair Coeffs  # lj/cut/coul/long", lines) in file order
  int atom_types = 0, bond_types = 0, angle_types = 0, dihedral_types = 0, improper_types = 0;
  struct Atom { int64_t id = 0, mol = 0; int type = 0; double q = 0; Vec3 x{0, 0, 0}; };
  std::vector<Atom> atoms;                                          // positions unwrapped (image flags applied)
  struct Term { int type = 0; std::vector<int64_t> ids; };
  std::vector<Term> bonds, angles, dihedrals, impropers;
};
LammpsFull read_lammps_full(const std::string& path);
// pair_coeffs: false leaves Pair Coeffs out (written to an include file instead, so a deck can switch pair styles)
void write_lammps_full(const LammpsFull& d, const std::string& path, bool pair_coeffs = true);

struct BackmapLibrary {
  struct Fragment { std::vector<int> type; std::vector<double> q; std::vector<Vec3> rel; std::array<Vec3, 2> dirs; };   // rel: from the centre of mass
  struct ClassInfo { std::string key, kind, role; std::vector<Fragment> instances; int atoms = 0; double mass = 0; };
  std::vector<ClassInfo> classes;
  std::map<std::pair<std::string, std::string>, int> by_kind_role;   // (kind, head|inner|tail|single) → class index
  struct Template { int kind = 0; int type = 0; std::vector<int> classes; std::vector<std::pair<int, int>> atoms; };   // atoms: (bead offset, role)
  std::vector<Template> templates;                                  // kind 0 bond, 1 angle, 2 dihedral, 3 improper
  LammpsFull reference;                                             // its masses and coefficients (atoms cleared)
  std::vector<std::string> notes;
};
BackmapLibrary backmap_library(const LammpsFull& aa, const Json& map_json, int conformers = 20);

struct BackmapResult {
  LammpsFull data;
  int beads = 0, terms_unmatched = 0;
  double cut_min = 0, cut_mean = 0, cut_max = 0;   // the restored cut bonds' lengths before relaxation (Å)
  double charge = 0;                                // net charge
  int reversed_chains = 0;
  std::vector<std::string> notes;
};
BackmapResult backmap_fragments(const BackmapLibrary& lib, const CgTopology& t, const std::vector<Vec3>& beads, const Cell& cell, uint64_t seed = 1);

// The staged relaxation deck; style lines (pair_style, kspace_style, special_bonds, bond_style …) from the reference input.
struct BackmapDeckOptions { double T = 300, P = 1; int64_t soft_steps = 20000, npt_steps = 20000; double dt_soft = 0.5, dt = 1.0; double soft_max = 30; };
std::string backmap_deck(const std::vector<std::string>& style_lines, const BackmapDeckOptions& o);
std::vector<std::string> style_lines_of(const std::string& lammps_input);

// Bond lengths and angles against harmonic coefficients: per type the count, mean and largest |Δ| (Å or degrees).
struct BackmapCheckRow { std::string kind; int type = 0; int count = 0; double ref = 0, mean_dev = 0, max_dev = 0; };
std::vector<BackmapCheckRow> backmap_check(const LammpsFull& d);

}  // namespace caps
