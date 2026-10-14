#!/bin/bash
# caps relax on a LAMMPS data file caps ff apply wrote with OPLS-AA 2024, without --ff: the file's own coefficients,
# charges, 1-4 weights and cut-off (its starting energy terms are caps ff apply's, to the six decimals the file holds), not the built-in GAFF / UFF; a file whose
# force field CAPS cannot rebuild (here: without CAPS's record of its 1-4 weights) is refused with a pointer to --ff.
# usage: test_cli_data_ff.sh CAPS SOURCE_DIR
set -u
CAPS=$1
ROOT=$2
d=$(mktemp -d)
trap 'rm -rf "$d"' EXIT
cd "$d" || exit 1
fail() { echo "FAIL: $*"; exit 1; }
# two "energy (kcal/mol): bond …" lines agree term by term to 2e-5 relative: the data file holds every coefficient at six
# decimals (a mixed pair's 0.105924501… as 0.105925), so the rebuilt force field differs from the one that wrote it there
same_terms() {
  awk -v a="$1" -v b="$2" 'BEGIN {
    na = split(a, x, " "); nb = split(b, y, " ")
    if (na != nb) exit 1
    for (i = 1; i <= na; ++i)
      if (x[i] ~ /^-?[0-9.]+$/) { d = x[i] - y[i]; if (d < 0) d = -d; s = x[i] < 0 ? -x[i] : x[i]; if (s < 1) s = 1; if (d > 2e-5 * s) exit 1 }
  }'
}

"$CAPS" build 'CC#N' -o acn.pdb > /dev/null || fail "caps build"
printf 'tolerance 2.0\nfiletype pdb\noutput box.pdb\npbc 0 0 0 16 16 16\nstructure acn.pdb\n  number 40\n  inside box 0 0 0 16 16 16\nend structure\n' > box.inp
"$CAPS" pack box.inp -o box.pdb --quiet > /dev/null || fail "caps pack"
apply=$("$CAPS" ff apply box.pdb --ff "$ROOT/data/forcefields/oplsaa2024.json" --cutoff 7.5 -o box.data | grep '^energy') || fail "caps ff apply"

out=$("$CAPS" relax box.data -o relaxed.data --cutoff 7.5 --iterations 20 2> /dev/null)
[ $? -le 3 ] || fail "caps relax: $out"   # 3: not converged in 20 iterations
echo "$out" | grep -q "^force field: OPLS-AA 2024 (Jorgensen) (the file's own" || fail "relax did not use the file's force field: $out"
start=$(echo "$out" | grep '^starting energy')
same_terms "${apply#energy (kcal/mol): }" "${start#starting energy (kcal/mol): }" || fail "energy terms differ
  caps ff apply: $apply
  caps relax:    $start"
echo "caps ff apply: $apply"
echo "caps relax:    $start"

# without --cutoff: the force field's own (OPLS-AA 2024's 12 Å, recorded in the data file), as caps ff apply uses it
apply12=$("$CAPS" ff apply box.pdb --ff "$ROOT/data/forcefields/oplsaa2024.json" | grep '^energy')
out12=$("$CAPS" relax box.data -o relaxed12.data --iterations 1 2> /dev/null)
echo "$out12" | grep -q "^cut-off: 12 Å (the force field's own" || fail "relax did not take the force field's cut-off: $out12"
start12=$(echo "$out12" | grep '^starting energy')
same_terms "${apply12#energy (kcal/mol): }" "${start12#starting energy (kcal/mol): }" || fail "energy terms at the force field's cut-off differ
  caps ff apply: $apply12
  caps relax:    $start12"
"$CAPS" md box.data -o md.data --steps 2 --thermo 1 2> /dev/null | grep -q "^cut-off: 12 Å (the force field's own" || fail "md did not take the force field's cut-off"
"$CAPS" relax box.data -o x.data --cutoff 9 --iterations 1 2> /dev/null | grep -q "^cut-off" && fail "--cutoff given, yet a cut-off line"
echo "cut-off 12 Å: $start12"

# --overlay A,B: both overlays merged, the later on top (the nitrile N of oplsaa2024-lj-refit.json)
F="$ROOT/data/forcefields"
"$CAPS" ff apply box.pdb --ff "$F/oplsaa2024.json" --overlay "$F/oplsaa2024-rubber-refit.json,$F/oplsaa2024-lj-refit.json" -o ov.data > /dev/null ||
  fail "caps ff apply with two overlays"
grep -q "^[0-9]* [0-9]* 0.250000 3.200000  # opls_753 opls_753" ov.data || fail "the second overlay's nitrile N is not in the file"

# the relaxed file carries the force field on
"$CAPS" relax relaxed.data -o again.data --cutoff 7.5 --iterations 5 2> /dev/null | grep -q "^force field: OPLS-AA 2024 (Jorgensen) (the file's own" ||
  fail "the relaxed file lost its force field"

# without CAPS's record of its 1-4 weights the coefficients are not used, and not replaced by GAFF / UFF either
grep -v '^# caps-ff' box.data > other.data
err=$("$CAPS" relax other.data -o x.data 2>&1) && fail "a file without its 1-4 weights was relaxed: $err"
echo "$err" | grep -q "Give --ff" || fail "no pointer to --ff: $err"
echo "$err" | grep -qi "GAFF\|UFF" && fail "fell back to GAFF / UFF: $err"
echo "refused: $err"
exit 0
