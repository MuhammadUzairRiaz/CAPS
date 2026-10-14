#!/bin/bash
# CAPS's CGenFF library, built on this machine from the MacKerell lab's toppar release (no CHARMM file or data
# converted from it is kept in the repository: the generated file is gitignored). The biomolecular CHARMM sets
# (CHARMM36 proteins / lipids / nucleic acids / carbohydrates, CHARMM22, CHARMM19) were removed on 2026-10-14: CAPS
# focuses on polymers, crystals and materials.
#   bench/ff/open/build_charmm.sh TOPPAR_DIR [OUT_DIR]
# TOPPAR_DIR: the unpacked toppar_c36_jul17.tgz (mackerell.umaryland.edu/charmm_ff.shtml).
set -e
T=${1:?toppar directory}; O=${2:-data/forcefields}
C=bench/ff/convert_charmm_prm.py
W=$T/toppar_water_ions.str
python3 $C $O/cgenff-toppar.json --name "CGenFF (toppar)" --prm $T/par_all36_cgenff.prm --prm $W --rtf $T/top_all36_cgenff.rtf \
  --reference "K. Vanommeslaeghe et al., J. Comput. Chem. 31, 671 (2010)" --version "toppar_c36_jul17 (CGenFF 3.0.1 / 4.0)" --typing ../typing/cgenff-toppar.typing.json --impropers-from-rtf
