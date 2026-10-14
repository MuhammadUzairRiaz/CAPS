#!/bin/bash
# CAPS's GAFF / GAFF2 libraries from AmberTools' own parameter files (dat/leap/parm), which are public domain.
# (The biomolecular AMBER sets, ff03 / ff14SB and GLYCAM, were removed on 2026-10-14: CAPS focuses on polymers, crystals
# and materials.)
#   bench/ff/open/build_amber.sh PARM_DIR [OUT_DIR]
# PARM_DIR holds gaff.dat, gaff211.dat and gaff2.dat (AmberTools, ambermd.org; its LICENSE places them in the public domain).
set -e
P=${1:?AmberTools dat/leap/parm directory}; O=${2:-data/forcefields}
C=bench/ff/convert_amber_dat.py
GAFF="J. Wang, R. M. Wolf, J. W. Caldwell, P. A. Kollman, D. A. Case, J. Comput. Chem. 25, 1157 (2004)"
python3 $C $O/gaff-1.81-ambertools.json --typing ../typing/gaff-1.81-ambertools.typing.json --name "GAFF 1.81 (AmberTools)" --dat $P/gaff.dat --reference "$GAFF" --version "1.81 (May 2017)"
python3 $C $O/gaff2-2.11-ambertools.json --typing ../typing/gaff2-2.11-ambertools.typing.json --name "GAFF2 2.11 (AmberTools)" --dat $P/gaff211.dat --reference "$GAFF" --version "2.11 (May 2016)"
python3 $C $O/gaff2-2.2.20-ambertools.json --typing ../typing/gaff2-2.2.20-ambertools.typing.json --name "GAFF2 2.2.20 (AmberTools)" --dat $P/gaff2.dat --reference "$GAFF" --version "2.2.20 (March 2021)"
