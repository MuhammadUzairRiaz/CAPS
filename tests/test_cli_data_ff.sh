#!/bin/bash
# caps relax on a LAMMPS data file caps ff apply wrote with OPLS-AA 2024, without --ff: the file's own coefficients,
# charges and 1-4 weights (its starting energy terms are caps ff apply's), not the built-in GAFF / UFF; a file whose
# force field CAPS cannot rebuild (here: without CAPS's record of its 1-4 weights) is refused with a pointer to --ff.
# usage: test_cli_data_ff.sh CAPS SOURCE_DIR
set -u
CAPS=$1
ROOT=$2
d=$(mktemp -d)
trap 'rm -rf "$d"' EXIT
cd "$d" || exit 1
fail() { echo "FAIL: $*"; exit 1; }

"$CAPS" build 'CC#N' -o acn.pdb > /dev/null || fail "caps build"
printf 'tolerance 2.0\nfiletype pdb\noutput box.pdb\npbc 0 0 0 16 16 16\nstructure acn.pdb\n  number 40\n  inside box 0 0 0 16 16 16\nend structure\n' > box.inp
"$CAPS" pack box.inp -o box.pdb --quiet > /dev/null || fail "caps pack"
apply=$("$CAPS" ff apply box.pdb --ff "$ROOT/data/forcefields/oplsaa2024.json" --cutoff 7.5 -o box.data | grep '^energy') || fail "caps ff apply"

out=$("$CAPS" relax box.data -o relaxed.data --cutoff 7.5 --iterations 20 2> /dev/null)
[ $? -le 3 ] || fail "caps relax: $out"   # 3: not converged in 20 iterations
echo "$out" | grep -q "^force field: OPLS-AA 2024 (Jorgensen) (the file's own" || fail "relax did not use the file's force field: $out"
start=$(echo "$out" | grep '^starting energy')
[ "${apply#energy (kcal/mol): }" = "${start#starting energy (kcal/mol): }" ] || fail "energy terms differ
  caps ff apply: $apply
  caps relax:    $start"
echo "caps ff apply: $apply"
echo "caps relax:    $start"

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
