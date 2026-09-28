// CAPS literature many-body potentials for a group of atoms (a crystal or filler in a polymer composite): Tersoff,
// Stillinger–Weber, Vashishta, EAM … as LAMMPS reads them from their own files. CAPS never writes or invents their
// parameters: the user gives the file (LAMMPS's potentials folder, the NIST Interatomic Potentials Repository, a paper's
// supplement) and CAPS checks it, maps the group's elements onto it and writes it next to the LAMMPS inputs.
//
// The group becomes one atom type per element, named by its symbol, with the standard atomic mass and no charge (the
// potential carries all the interactions inside the group). Its pairs with the rest of the system are Lennard-Jones:
// UFF's parameters for the element (Rappé et al. 1992), mixed with the other group's by the cross rule chosen, or given
// explicitly — the usual treatment of a Tersoff / AIREBO filler in a polymer (the potential inside, van der Waals across).
//
// LAMMPS reads Tersoff, Stillinger–Weber, Vashishta, Gao–Weber and EAM files in metal units (eV) and converts them to the
// real units (kcal/mol) CAPS writes, when the file says its units on its first line ("UNITS: metal"); a file that does not
// say is read as it stands, so its units must be given (the copy written with the inputs then says them). AIREBO, REBO,
// MEAM, BOP and COMB are read in metal units only: they are refused, with the reason, since CAPS's inputs are in real units.
#pragma once

#include <string>
#include <vector>

#include "caps/field.hpp"
#include "caps/system.hpp"

namespace caps {

struct ManyBodySpec {
  std::string style;   // LAMMPS pair style
  std::string file;    // the potential file
  std::string units;   // metal | real; "" to take the file's own "UNITS:" tag
  // the Lennard-Jones form of the cross pairs, the other groups' own: lj12-6 (σ = x / 2^(1/6)) or lj9-6 (class II, PCFF /
  // COMPASS: σ = the minimum x), UFF's well depth D and minimum x either way
  std::string pair_form = "lj12-6";
};

// The styles CAPS writes (the file format of each is checked), and whether a style is one LAMMPS reads in metal units only.
const std::vector<std::string>& manybody_styles();

// Checks the file against the style and the group's elements (every element triplet, or every element, has its entry)
// and returns the group's force field (see above). Throws FieldError with the reason.
ForceField manybody_part(const System& group, const ManyBodySpec& spec, std::vector<std::string>* notes = nullptr);

// The name the potential file is written under beside the inputs: its own, or NAME-UNITS.EXT when CAPS adds the units
// (never the user's file changed).
std::string manybody_file_name(const ManyBodyFile& mb);

// Writes the potential file into dir under manybody_file_name, adding its units to the first line when the file did not say
// them; returns the file name written (no directory). Throws when it cannot be read or written.
std::string write_manybody_file(const ManyBodyFile& mb, const std::string& dir);

}  // namespace caps
