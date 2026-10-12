#!/bin/bash
# CAPS's typing against every DL_FIELD library that has a CAPS force field (bench/ff/validate_typing.py): whole
# molecules and residues with their open valences capped; DL_FIELD's own re-typing applied where it was run
# (validate_family.py work directories in $R/dlf_<ff>). One line per library; details in $R/typing_<ff>_dl.log.
#   DLFIELD=~/project/dl_f_4.13 bench/typing/sweep_dlfield.sh [OUTDIR]
cd "$(dirname "$0")/../.."
R=${1:-docs/paper/runs/checks_1012}
CH="water_spce?=data/typing/charmm-spc.typing.json;water_tip3p=data/typing/charmm-tip3p.typing.json;water_tip3p_pme=data/typing/charmm-tip3p-pme.typing.json"
AM="water_tip3p=data/typing/amber-tip3p.typing.json;water_tip3p_pme=data/typing/amber-tip3p-pme.typing.json"
CV="water_tip3p=data/typing/cvff-tip3p.typing.json;water_spce?=data/typing/cvff-spc.typing.json"
OP="water_tip3p_c=data/typing/opls2005-tip3p-c.typing.json;water_spce?=data/typing/opls2005-spc.typing.json"
while read lib ff ov; do
  [ -z "$lib" ] && continue
  case "$ov" in charmm) o="$CH";; amber) o="$AM";; cvff) o="$CV";; opls) o="$OP";; *) o="";; esac
  d=${ff%-frc}   # the .frc entries take the DL_FIELD run of their DL-derived twin
  dl=""; [ -d "$R/dlf_$d" ] && dl="--dl $R/dlf_$d"
  ex='tip4p|tip5p|TIP4P|TIP5P'
  [ "$ff" = opls2005 ] && ex="$ex|popc"   # DL_FIELD cannot read POPC's five-character atom names: no output to compare
  python3 bench/ff/validate_typing.py "$lib" "data/forcefields/$ff.json" "data/typing/$ff.typing.json" $dl --show 40 --all --cap-open \
    --exclude "$ex" ${o:+--overlay "$o"} > "$R/typing_${ff}_dl.log" 2>&1
  head -1 "$R/typing_${ff}_dl.log" | cut -c1-200
done <<'LIST'
PCFF pcff-frc
PCFF pcff
COMPASS compass
COMPASS compass-frc
CVFF cvff cvff
OPLS2005 opls2005 opls
AMBER amber amber
AMBER16_gaff gaff-amber16
AMBER25_gaff gaff-amber25
CHARMM charmm charmm
CHARMM19 charmm19
CHARMM22_prot charmm22-prot charmm
CHARMM36_prot charmm36-prot charmm
CHARMM36_lipid charmm36-lipid charmm
CHARMM36_nucl charmm36-nucl charmm
CHARMM36_carb charmm36-carb charmm
CHARMM36_cgenff cgenff charmm
DREIDING dreiding
GROMOS_G54A7 gromos-54a7
OPLS_AAM opls-aam
OPLS_CL_P opls-clp
OPLS_DES opls-des
TRAPPE_EH trappe-eh
TRAPPE_UA trappe-ua
MISC_FF misc
INORGANIC_binary_halides inorganic-binary-halides
INORGANIC_binary_misc inorganic-binary-misc
INORGANIC_binary_oxides inorganic-binary-oxides
INORGANIC_clay inorganic-clay
INORGANIC_glass inorganic-glass
INORGANIC_ternary_oxides inorganic-ternary-oxides
INORGANIC_zeolite inorganic-zeolite
INORGANIC_zeolite_Hill_Sauer inorganic-zeolite-hill-sauer
LIST
