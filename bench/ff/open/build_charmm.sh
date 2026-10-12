#!/bin/bash
# CAPS's CHARMM libraries, built on this machine from the MacKerell lab's toppar release (no CHARMM file or data
# converted from it is kept in the repository: the generated files are gitignored).
#   bench/ff/open/build_charmm.sh TOPPAR_DIR [OUT_DIR]
# TOPPAR_DIR: the unpacked toppar_c36_jul17.tgz (mackerell.umaryland.edu/charmm_ff.shtml).
set -e
T=${1:?toppar directory}; O=${2:-data/forcefields}
C=bench/ff/convert_charmm_prm.py
W=$T/toppar_water_ions.str
streams() { for f in $T/stream/$1/*.str; do printf -- "--prm %s " "$f"; done; }
REF="A. D. MacKerell Jr. et al., J. Phys. Chem. B 102, 3586 (1998)"
python3 $C $O/charmm36-prot-toppar.json --name "CHARMM36 proteins (toppar)" --prm $T/par_all36_prot.prm --prm $W $(streams prot) --rtf $T/top_all36_prot.rtf \
  --reference "R. B. Best et al., J. Chem. Theory Comput. 8, 3257 (2012)" --reference "$REF" --version "toppar_c36_jul17" --typing ../typing/charmm36-prot-toppar.typing.json
python3 $C $O/charmm36-lipid-toppar.json --name "CHARMM36 lipids (toppar)" --prm $T/par_all36_lipid.prm --prm $W $(streams lipid) --rtf $T/top_all36_lipid.rtf \
  --reference "J. B. Klauda et al., J. Phys. Chem. B 114, 7830 (2010)" --version "toppar_c36_jul17" --typing ../typing/charmm36-lipid-toppar.typing.json
python3 $C $O/charmm36-nucl-toppar.json --name "CHARMM36 nucleic acids (toppar)" --prm $T/par_all36_na.prm --prm $W $(streams na) --rtf $T/top_all36_na.rtf \
  --reference "K. Hart et al., J. Chem. Theory Comput. 8, 348 (2012)" --version "toppar_c36_jul17" --typing ../typing/charmm36-nucl-toppar.typing.json
python3 $C $O/charmm36-carb-toppar.json --name "CHARMM36 carbohydrates (toppar)" --prm $T/par_all36_carb.prm --prm $W $(streams carb) --rtf $T/top_all36_carb.rtf \
  --reference "O. Guvench et al., J. Chem. Theory Comput. 7, 3162 (2011)" --version "toppar_c36_jul17" --typing ../typing/charmm36-carb-toppar.typing.json
python3 $C $O/cgenff-toppar.json --name "CGenFF (toppar)" --prm $T/par_all36_cgenff.prm --prm $W --rtf $T/top_all36_cgenff.rtf \
  --reference "K. Vanommeslaeghe et al., J. Comput. Chem. 31, 671 (2010)" --version "toppar_c36_jul17 (CGenFF 3.0.1 / 4.0)" --typing ../typing/cgenff-toppar.typing.json
python3 $C $O/charmm22-prot-toppar.json --name "CHARMM22 proteins (toppar)" --prm $T/par_all22_prot.prm --prm $W --rtf $T/top_all22_prot.rtf \
  --reference "$REF" --version "toppar_c36_jul17" --typing ../typing/charmm22-prot-toppar.typing.json
python3 $C $O/charmm19-toppar.json --name "CHARMM19 (toppar)" --prm $T/param19.inp --rtf $T/toph19.inp --one-torsion-per-bond \
  --reference "E. Neria, S. Fischer, M. Karplus, J. Chem. Phys. 105, 1902 (1996)" --version "param19 / toph19" --typing ../typing/charmm19-toppar.typing.json
