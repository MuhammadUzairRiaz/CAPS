#!/usr/bin/env bash
# Cross-check CAPS Field against LAMMPS on one relaxed cell: energy terms, pressure and per-atom forces.
# Needs an `lmp` executable with the MOLECULE and EXTRA-MOLECULE / EXTRA-PAIR packages (fourier, cvff, lj/cut/coul/dsf).
set -euo pipefail
cd "$(dirname "$0")/.."
LMP=${LMP:-lmp}
command -v "$LMP" >/dev/null || { echo "no LAMMPS executable ($LMP); set LMP=/path/to/lmp"; exit 2; }
out=$(mktemp -d)
build/cli/caps grow -o "$out/grown.data" --chains 6 --dp 6 --density 0.4 --seed 4 >/dev/null
build/cli/caps relax "$out/grown.data" -o "$out/cell.data" --density 0.9 --quiet >/dev/null
build/cli/caps field "$out/cell.data" --forces "$out/caps_f.txt" > "$out/caps.txt"
cat > "$out/in.check" <<IN
units real
atom_style full
pair_style lj/cut/coul/dsf 0.2 10.0
pair_modify mix arithmetic tail yes
bond_style harmonic
angle_style harmonic
dihedral_style fourier
improper_style cvff
special_bonds amber
read_data $out/cell.data
neighbor 1.0 bin
thermo_style custom ebond eangle edihed eimp evdwl ecoul press
thermo_modify format float %.6f
run 0
write_dump all custom $out/lmp_f.dump id fx fy fz modify sort id format float %.8f
IN
"$LMP" -in "$out/in.check" -log "$out/log.lammps" -screen none
python3 - "$out" <<'PY'
import re, sys
d = sys.argv[1]
log = open(f"{d}/log.lammps").read().splitlines()
i = next(k for k, l in enumerate(log) if l.split()[:1] == ["E_bond"])
lmp = dict(zip(["bond", "angle", "dihedral", "improper", "vdW", "Coulomb", "pressure"], map(float, log[i + 1].split())))
caps = open(f"{d}/caps.txt").read()
e = {k: float(v) for k, v in re.findall(r"(bond|angle|dihedral|improper|vdW|Coulomb) (-?[\d.]+)", caps)}
e["pressure"] = float(re.search(r"pressure (-?\d+) atm", caps).group(1))
print(f"{'term':10s} {'CAPS':>14s} {'LAMMPS':>14s}")
for k in ["bond", "angle", "dihedral", "improper", "vdW", "Coulomb", "pressure"]:
    print(f"{k:10s} {e[k]:14.2f} {lmp[k]:14.2f}")
cf = {int(l.split()[0]): list(map(float, l.split()[1:])) for l in open(f"{d}/caps_f.txt")}
worst = big = 0.0
for l in open(f"{d}/lmp_f.dump").read().splitlines()[9:]:
    i, *f = l.split()
    f = list(map(float, f))
    worst = max(worst, sum((a - b) ** 2 for a, b in zip(f, cf[int(i)])) ** 0.5)
    big = max(big, sum(a * a for a in f) ** 0.5)
print(f"largest per-atom force difference {worst:.2e} kcal/mol/Å (largest force {big:.3f})")
print("CAPS reports the pressure rounded to 1 atm.")
sys.exit(0 if worst < 1e-3 else 1)
PY

# NVE trajectory: same start (positions and velocities), 200 steps of 1 fs in both codes.
build/cli/caps md "$out/cell.data" -o "$out/start.data" --steps 0 --temp 300 --seed 7 --quiet >/dev/null
build/cli/caps md "$out/start.data" -o "$out/end.data" --steps 200 --thermostat none --dump "$out/caps.lammpstrj" --every 100 --quiet >/dev/null
sed -e "s#read_data .*#read_data $out/start.data#" -e '/thermo_style/,$d' "$out/in.check" > "$out/in.nve"
cat >> "$out/in.nve" <<IN
neigh_modify every 1 delay 0 check yes
timestep 1.0
fix 1 all nve
dump d all custom 100 $out/lmp.lammpstrj id xu yu zu
dump_modify d sort id format float %.6f
run 200
IN
"$LMP" -in "$out/in.nve" -log none -screen none
echo "NVE, 200 steps of 1 fs from the same positions and velocities:"
python3 scripts/compare_dumps.py "$out/caps.lammpstrj" "$out/lmp.lammpstrj" 1e-3

# Bond constraints: bonds to hydrogen held by SHAKE/RATTLE (fix rattle in LAMMPS), 200 steps of 2 fs from a start that
# already meets them (caps projects positions and velocities onto the constraints before its first step).
build/cli/caps md "$out/cell.data" -o "$out/start_c.data" --steps 0 --temp 300 --seed 7 --constraints h-bonds --quiet >/dev/null
build/cli/caps md "$out/start_c.data" -o "$out/end_c.data" --steps 200 --dt 2 --thermostat none --constraints h-bonds --dump "$out/caps_c.lammpstrj" --every 100 --quiet >/dev/null
sed -e "s#read_data .*#read_data $out/start_c.data#" -e '/thermo_style/,$d' "$out/in.check" > "$out/in.rattle"
cat >> "$out/in.rattle" <<IN
neigh_modify every 1 delay 0 check yes
timestep 2.0
fix 1 all nve
fix 2 all rattle 1e-10 500 0 m 1.008
dump d all custom 100 $out/lmp_c.lammpstrj id xu yu zu
dump_modify d sort id format float %.6f
run 200
IN
"$LMP" -in "$out/in.rattle" -log none -screen none
echo "NVE with bonds to hydrogen constrained (fix rattle), 200 steps of 2 fs from the same start:"
python3 scripts/compare_dumps.py "$out/caps_c.lammpstrj" "$out/lmp_c.lammpstrj" 1e-3

# Nosé–Hoover chains (fix nvt) and MTK pressure coupling (fix npt iso): 400 steps of 1 fs from the same start
build/cli/caps md "$out/start.data" -o "$out/end_nvt.data" --steps 400 --thermostat nose-hoover --tau-t 100 --dump "$out/caps_nvt.lammpstrj" --every 200 --quiet >/dev/null
build/cli/caps md "$out/start.data" -o "$out/end_npt.data" --steps 400 --thermostat nose-hoover --tau-t 100 --barostat mtk --pressure 1 --tau-p 1000 \
  --dump "$out/caps_npt.lammpstrj" --every 200 --quiet >/dev/null
for kind in nvt npt; do
  fix="fix 1 all nvt temp 300 300 100.0"
  [ "$kind" = npt ] && fix="fix 1 all npt temp 300 300 100.0 iso 1 1 1000.0"
  sed -e "s#read_data .*#read_data $out/start.data#" -e '/thermo_style/,$d' "$out/in.check" > "$out/in.$kind"
  cat >> "$out/in.$kind" <<IN
neigh_modify every 1 delay 0 check yes
timestep 1.0
$fix
dump d all custom 200 $out/lmp_$kind.lammpstrj id xu yu zu
dump_modify d sort id format float %.6f
run 400
IN
  "$LMP" -in "$out/in.$kind" -log none -screen none
  echo "Nosé–Hoover $kind vs LAMMPS fix $kind, 400 steps of 1 fs from the same start:"
  python3 scripts/compare_dumps.py "$out/caps_$kind.lammpstrj" "$out/lmp_$kind.lammpstrj" 1e-3
done
