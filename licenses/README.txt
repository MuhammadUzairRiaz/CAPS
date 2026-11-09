Reference data built into CAPS and where it comes from

UFF parameters (core/src/uff_params.inc)
  Rappé, Casewit, Colwell, Goddard & Skiff, J. Am. Chem. Soc. 114, 10024 (1992), Tables I and III;
  Rappé & Goddard, J. Phys. Chem. 95, 3358 (1991), Table I; the other elements' GMP electronegativities
  as UFF implementations share them (ref. 13 of the UFF paper).

Element names and mass numbers (core/src/element_extra.inc)
  IUPAC / CIAAW: Standard Atomic Weights and Isotopic Compositions of the Elements (2024);
  NIST: Atomic Weights and Isotopic Compositions; NUBASE2020 (Kondev et al., Chin. Phys. C 45, 030001 (2021)).

Pauling electronegativities and electron configurations (core/src/element_extra.inc)
  Blue Obelisk Data Repository (MIT licence), as distributed with VTK 9.6.1: see BlueObelisk-MIT.txt.

Libraries shipped inside the installers built with CAPS's ML potential (caps ml, caps ff complete)
  libtorch / PyTorch 2.14.1 (BSD-3-Clause, with the third-party notices of its LICENSE file): see PyTorch-BSD-3.txt.
  CAPS's model data/ml/caps-mace-v1.model is CAPS's own, trained by CAPS (MACE architecture, Batatia et al., MIT licence).

Force-field and potential data (each library's source and licence are in data/forcefields/catalogue.json)
  LAMMPS potential files in data/potentials/lammps, unchanged with their headers: GNU GPL version 2, see LAMMPS-GPL-2.0.txt.
  GROMACS OPLS-AA (data/forcefields/oplsaa-gromacs.json): GNU LGPL 2.1, see LGPL-2.1.txt.
  CL&P and its deep-eutectic set: MIT, see CLandP-MIT.txt.
  TraPPE-UA as distributed with foyer (data/forcefields/trappe-ua-foyer.json): MIT, see foyer-MIT.txt.
  PCFF, CVFF, COMPASS: the .frc files as distributed with LAMMPS (tools/msi2lmp/frc_files).
  INTERFACE force field 1.5: files of the Heinz group. Martini 3: files of the Martini developers (cgmartini.nl).

Libraries in the Studio and the build
  Avalonia, .NET: MIT. zlib: zlib licence. IBM Plex fonts: SIL Open Font License, see studio/CapsStudio/Assets/Fonts/LICENSE.txt.
