# Literature many-body potentials

Potential files for the crystal or filler group of a composite (Field · by group · Literature potential): Tersoff,
Stillinger–Weber, Vashishta, EAM, AIREBO, AIREBO-M, REBO and MEAM parameters as LAMMPS reads them. `catalogue.json`
lists each file with its pair style, the elements it covers, what it is used for and the paper to cite.

The files in `lammps/` are copied unchanged from the `potentials` folder of the LAMMPS distribution
(https://www.lammps.org, GNU General Public License version 2); each keeps its header with the date, the units, the
contributor and the citation. CAPS does not refit them: a system that uses one gets the file written beside its `.data`
and `.in` files. How the units are handled depends on the style:

- Tersoff, Stillinger–Weber, Vashishta and EAM: the file is written unchanged, and LAMMPS converts its metal units (eV)
  to the real units (kcal/mol) of the CAPS inputs.
- AIREBO, AIREBO-M and REBO: LAMMPS reads these in metal units only, so CAPS writes the inputs in metal units (eV, ps,
  bar) by default. If real units are asked for, CAPS writes its own real-unit copy of the file beside the inputs, with
  the energy parameters (A, B and the epsilons) multiplied by 23.060549; the copy in `lammps/` stays as it is.
- MEAM: LAMMPS reads these in metal units only and CAPS has no conversion for them, so the inputs are written in metal
  units, and asking for real units is refused.

Other files (the NIST Interatomic Potentials Repository, a paper's supplement) can be chosen directly in Field; a file
that does not state its units on its first line needs them given there.
