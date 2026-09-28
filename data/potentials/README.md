# Literature many-body potentials

Potential files for the crystal or filler group of a composite (Field · by group · Literature potential): Tersoff,
Stillinger–Weber, Vashishta and EAM parameters as LAMMPS reads them. `catalogue.json` lists each file with its pair
style, the elements it covers, what it is used for and the paper to cite.

The files in `lammps/` are copied unchanged from the `potentials` folder of the LAMMPS distribution
(https://www.lammps.org, GNU General Public License version 2); each keeps its header with the date, the units, the
contributor and the citation. CAPS does not modify or refit them: a system that uses one gets the file written beside
its `.data` and `.in` files, and LAMMPS converts its metal units (eV) to the real units (kcal/mol) of the CAPS inputs.

Other files (the NIST Interatomic Potentials Repository, a paper's supplement) can be chosen directly in Field; a file
that does not state its units on its first line needs them given there.
