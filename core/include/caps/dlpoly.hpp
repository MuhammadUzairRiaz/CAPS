// CAPS DL_POLY export: FIELD, CONFIG and a generic CONTROL for DL_POLY 4 (units kcal), in the conventions the DL_POLY
// force-field tools write them: harmonic terms as ½k (k twice the force field's k), class II bonds and angles as
// "quar" (2K2, 3K3, 4K4), torsions as "cos" A [1 + cos(mφ − δ)] (OPLS triples as "cos3"), the 1-4 scale factors on one
// dihedral term per 1-4 pair (a zero-amplitude term when the force field has no torsion there), inversions with the
// centre first, and every van der Waals type pair written out ("lj", "nm" 9-6 for class II, "buck").
// DL_POLY has no class II cross terms, no separate 1-4 Lennard-Jones parameters and no combined bending–torsion: those
// are left out and named in the notes (its energy is then the force field's diagonal part only).
#pragma once

#include <string>
#include <vector>

#include "caps/field.hpp"
#include "caps/system.hpp"

namespace caps {

struct DlpolyOptions {
  std::string title = "CAPS";
  double cutoff = 10.0;          // Å (CONTROL rcut)
  double temperature = 300;      // K (CONTROL)
  double timestep_fs = 1.0;
  long steps = 10000;
  bool group_molecules = true;   // consecutive identical molecules become one molecular type (nummols N)
};

// Writes FIELD, CONFIG and CONTROL into dir (created when missing) and returns notes (what was left out, how atoms
// were ordered). Throws std::invalid_argument for terms DL_POLY cannot express at all (e.g. SDK pairs, virtual sites).
std::vector<std::string> write_dlpoly(const System& s, const ForceField& ff, const std::string& dir, const DlpolyOptions& o = {});

}  // namespace caps
