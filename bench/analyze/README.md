# Analyze checks

`check_analyze.py DATA DUMP [--frame-ps 0.5]` compares `caps analyze` with independent implementations on the same
files: MDAnalysis (g(r), Rg, atom MSD), numpy (MSD with drift removed, D, the direct reciprocal-lattice S(q)), Monte
Carlo point insertion (free fraction) and an analytic lattice (largest pore). It needs a Python with MDAnalysis and
scipy.

```bash
PY=/path/to/python-with-mdanalysis
./build/cli/caps md out/eq_full.data -o out/ps500.data --dump out/ps500.lammpstrj --steps 100000 --dt 1 --temp 500 --every 500 --new-velocities --seed 7
$PY bench/analyze/check_analyze.py out/eq_full.data out/ps500.lammpstrj --frame-ps 0.5
```

2026-09-24: 10 of 10 checks pass (see the table in the top-level README).
