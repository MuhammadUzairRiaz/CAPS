#!/bin/bash
# CAPS's AMBER-family libraries from AmberTools' own parameter files (dat/leap/parm). GAFF / GAFF2 / ff03 are committed
# (the files are public domain); GLYCAM is built locally only.
#   bench/ff/open/build_amber.sh PARM_DIR [OUT_DIR]
# PARM_DIR holds gaff.dat, gaff211.dat, gaff2.dat, parm99.dat, frcmod.ff03 (AmberTools, ambermd.org; its LICENSE places them in
# the public domain) and, for the sugar library, GLYCAM_06h.dat (Woods group, its own copyright notice: built here, never shipped).
set -e
P=${1:?AmberTools dat/leap/parm directory}; O=${2:-data/forcefields}
C=bench/ff/convert_amber_dat.py
GAFF="J. Wang, R. M. Wolf, J. W. Caldwell, P. A. Kollman, D. A. Case, J. Comput. Chem. 25, 1157 (2004)"
python3 $C $O/gaff-1.81-ambertools.json --typing ../typing/gaff-1.81-ambertools.typing.json --name "GAFF 1.81 (AmberTools)" --dat $P/gaff.dat --reference "$GAFF" --version "1.81 (May 2017)"
python3 $C $O/gaff2-2.11-ambertools.json --typing ../typing/gaff2-2.11-ambertools.typing.json --name "GAFF2 2.11 (AmberTools)" --dat $P/gaff211.dat --reference "$GAFF" --version "2.11 (May 2016)"
python3 $C $O/gaff2-2.2.20-ambertools.json --typing ../typing/gaff2-2.2.20-ambertools.typing.json --name "GAFF2 2.2.20 (AmberTools)" --dat $P/gaff2.dat --reference "$GAFF" --version "2.2.20 (March 2021)"
python3 $C $O/amber-ff03-ambertools.json --typing ../typing/amber-ff03-ambertools.typing.json --name "AMBER ff03 (parm99 + frcmod.ff03, AmberTools)" --dat $P/parm99.dat --frcmod $P/frcmod.ff03 \
  --reference "Y. Duan et al., J. Comput. Chem. 24, 1999 (2003)" --reference "J. Wang, P. Cieplak, P. A. Kollman, J. Comput. Chem. 21, 1049 (2000)"
if [ -f $P/GLYCAM_06h.dat ]; then
  # GLYCAM sets SCEE = SCNB = 1 on its torsions; CAPS applies one 1-4 scaling per structure, so GLYCAM is a library of its own
  python3 $C $O/glycam06h-ambertools.json --name "GLYCAM06h (AmberTools)" --dat $P/GLYCAM_06h.dat --scale14 1 1 \
    --reference "K. N. Kirschner et al., J. Comput. Chem. 29, 622 (2008)" --version "GLYCAM_06h"
fi
