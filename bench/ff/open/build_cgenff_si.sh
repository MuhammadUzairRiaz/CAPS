#!/bin/bash
# CGenFF v. 2a5, the version published with the CGenFF paper, from the paper's Supporting Information:
# K. Vanommeslaeghe et al., J. Comput. Chem. 31, 671 (2010), doi:10.1002/jcc.21367 (Wiley's jcc_21367_sm_suppinfo.tgz).
#   bench/ff/open/build_cgenff_si.sh SUPP_DIR [OUT_DIR]
# SUPP_DIR: the unpacked SI (supp/, with par_all36_cgenff.prm and top_all36_cgenff.rtf). Writes cgenff-2a5-jcc2010.json;
# the typing rules come from bench/typing/make_cgenff_rules.py (data/typing/cgenff-2a5-jcc2010.typing.json).
set -e
S=${1:?SI directory}; O=${2:-data/forcefields}
python3 bench/ff/convert_charmm_prm.py $O/cgenff-2a5-jcc2010.json --name "CGenFF 2a5 (J. Comput. Chem. 2010)" \
  --prm $S/par_all36_cgenff.prm --rtf $S/top_all36_cgenff.rtf \
  --source "Supporting Information of K. Vanommeslaeghe et al., J. Comput. Chem. 31, 671 (2010), doi:10.1002/jcc.21367 (jcc_21367_sm_suppinfo.tgz): par_all36_cgenff.prm and top_all36_cgenff.rtf, CGenFF v. 2a5" \
  --reference "K. Vanommeslaeghe, E. Hatcher, C. Acharya, S. Kundu, S. Zhong, J. Shim, E. Darian, O. Guvench, P. Lopes, I. Vorobyov, A. D. MacKerell Jr., J. Comput. Chem. 31, 671 (2010), doi:10.1002/jcc.21367" \
  --version "2a5 (J. Comput. Chem. 2010 Supporting Information)" --typing ../typing/cgenff-2a5-jcc2010.typing.json --impropers-from-rtf \
  --note "parameter values as published in the Supporting Information of Vanommeslaeghe et al. 2010; CAPS's conversion" \
  --note "the file's NONBONDED line asks for an energy switch of the Lennard-Jones from 10 to 12 A (vswitch) and shifted electrostatics; CAPS uses CHARMM's force switch over the same range (vfswitch, as later CGenFF releases ask for) and long-range electrostatics" \
  --note "no water or ions: the Supporting Information's files define none (TIP3P comes with the CHARMM36 libraries)" \
  --note "atom types added in later CGenFF releases (thiocarbonyl, sulfone, sulfoxide, thiophosphate, epoxide, cyclobutyl, beta-lactam N, bipyrrole and N-aryl ring types, the terminal alkyne CG1T2) are not in 2a5; its typing rules leave them out"
