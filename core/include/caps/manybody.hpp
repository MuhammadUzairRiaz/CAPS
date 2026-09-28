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
// say is read as it stands, so its units must be given (the copy written with the inputs then says them). AIREBO, AIREBO-M
// and REBO (carbon and hydrogen: nanotubes, graphene) are read by LAMMPS in metal units only: a system with one is written
// in metal units (eV, ps, bar), every parameter of the force field converted, or — asked for — in real units with a copy
// of the file CAPS converts (A, B and the ε's; checked to reproduce LAMMPS's metal-unit energy and forces). MEAM (metal
// units only) maps each element to a library entry chosen by the user (by default the first of that atomic number). BOP
// and COMB (charge equilibration) are refused with the reason.
#pragma once

#include <map>
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
  std::string args;    // the pair style's arguments; "" the usual ones (airebo, airebo/morse: 3.0 1 1 — LJ cut-off 3σ, LJ and torsion on)
  // MEAM: file is the element library (library.meam), file2 the alloy parameter file ("" none); entries maps an element to
  // the library entry it takes ('SiS', 'Ni4' …) — by default the first entry of that atomic number
  // the order counts: the parameter file names elements by their index in this list (SiC.meam: 1 Si, 2 C), so every
  // entry it refers to is extracted, in this order, even for an element the group does not have
  std::string file2;
  std::vector<std::pair<std::string, std::string>> entries;
};

// A MEAM library's entries, as LAMMPS reads them (the first of a repeated name counts): name, lattice, atomic number, mass.
struct MeamEntry { std::string name, lattice; int z = 0; double mass = 0; };
std::vector<MeamEntry> meam_library(const std::string& path);

// The styles CAPS writes (the file format of each is checked), and whether a style is one LAMMPS reads in metal units only.
const std::vector<std::string>& manybody_styles();

// Checks the file against the style and the group's elements (every element triplet, or every element, has its entry)
// and returns the group's force field (see above). Throws FieldError with the reason.
ForceField manybody_part(const System& group, const ManyBodySpec& spec, std::vector<std::string>* notes = nullptr);

// The name the potential file is written under beside the inputs: its own, or NAME-UNITS.EXT when CAPS adds the units
// (never the user's file changed).
std::string manybody_file_name(const ManyBodyFile& mb, const std::string& units = "");
// MEAM's parameter file beside the inputs (its own name)
std::string manybody_file2_name(const ManyBodyFile& mb);

// Writes the potential file into dir under manybody_file_name, adding its units to the first line when the file did not say
// them; returns the file name written (no directory). Throws when it cannot be read or written.
// units: the inputs' ("real" or "metal"). AIREBO / AIREBO-M / REBO in real units are written converted by CAPS (LAMMPS
// cannot convert them): the energy parameters A, B, ε (LJ, torsion, Morse) times 23.060549, the dimensionless splines as
// they are, the first line saying UNITS: real. The other styles are written as given (LAMMPS converts them itself).
std::string write_manybody_file(const ManyBodyFile& mb, const std::string& dir, const std::string& units = "");
// Whether CAPS writes a real-units copy of this style's file (the Brenner family).
bool manybody_caps_converts(const std::string& style);

}  // namespace caps
