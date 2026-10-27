#!/usr/bin/env bash
# Regenerates the Studio screenshots of the website (docs/manual/img/*.jpg): every picture is the real Studio,
# rendered off-screen at 1440 x 900 by CapsStudio --screenshot, after the module has run on a real system.
#
#   tools/docs/make_screenshots.sh              every picture
#   tools/docs/make_screenshots.sh 07_relax ... only those
#
# Needs the core (build/, cmake Release) and the Studio (dotnet build -c Release in studio/CapsStudio).
# Runs stay short (about 5,000 MD steps at most, no cell over 15,000 atoms). PNGs and work files go to $SHOT_WORK
# (default: a temporary folder); the JPEGs replace the ones in docs/manual/img.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
DLL="$ROOT/studio/CapsStudio/bin/Release/net10.0/CapsStudio.dll"
CAPS="$ROOT/build/cli/caps"
IMG="$ROOT/docs/manual/img"
WORK="${SHOT_WORK:-$(mktemp -d -t caps-shots)}"
mkdir -p "$WORK"

# shot NAME FILE ARGS...: the Studio opened on FILE ("" for none; TOPO=DATA for a trajectory's topology), driven by ARGS,
# saved as NAME.jpg
shot() {
  local name="$1" file="$2"; shift 2
  echo "== $name"
  # a fresh settings folder (jobs, recent files, built molecules) for each picture; the user's own are never touched
  rm -rf "$WORK/state-$name"; mkdir -p "$WORK/state-$name"
  CAPS_SETTINGS="$WORK/state-$name/settings.json" CAPS_RECENT_DIR="$WORK/state-$name/recent" \
  dotnet "$DLL" --screenshot "$WORK/$name.png" "$file" "${TOPO:-}" "$@"
  sips -s format jpeg -s formatOptions 85 "$WORK/$name.png" --out "$IMG/$name.jpg" >/dev/null
}

# Real runs of atactic polystyrene (10 chains x 30 units, 4,820 atoms, OPLS-AA 2024) on an HPC cluster, kept
# outside the repository in $PS_RUNS (default $HOME/project/CAPS-builds/ps-barnard-runs):
#   equil.csv, ps-equil.*.log, ps_eq.data  caps equilibrate ps_min.data --protocol larsen21 --dt 0.5 --tfinal 300 --tmax 600
#                                          --until-converged --block 20 --log equil.csv -o ps_eq.data (its printed output: the .log)
#   ps_300K.data, ps_300K.lammpstrj        2 ns NPT at 300 K and 1 atm from ps_eq.data, a frame every 10 ps
#   tg600_1 .. tg600_3, *.log              caps tg 600 -> 300 K, three seeds (each folder: ps_tg.tg-points.json, csv/, ps_tg.data)
# The pictures copy what they need to $DEMO and run from there, so no home-folder path shows in them.
PS_RUNS="${PS_RUNS:-$HOME/project/CAPS-builds/ps-barnard-runs}"
DEMO=/tmp/caps-demo
ps_runs() {
  if [ ! -d "$PS_RUNS" ]; then echo "$1 needs the polystyrene cluster runs in \$PS_RUNS ($PS_RUNS); skipped" >&2; return 1; fi
  mkdir -p "$DEMO/$2"
}

# the polystyrene cell most pages start from: 10 atactic chains of 30 units grown at 0.90 g/cm3, OPLS-AA 2024 assigned
PS="growcell=P013:oplsaa2024:10:30:0.9"

# a project folder of real runs: three recipes (build, type, grow, relax, 1 ps of MD, analysis), each with the
# provenance CAPS writes beside its cell
project_folder() {
  local P="$WORK/project"
  rm -rf "$P"; mkdir -p "$P"
  cat > "$WORK/pe_cell.yaml" <<'Y'
recipe: 1
name: pe_cell
build:
  polymer: { smiles: "[*]CC[*]", dp: 40, chains: 8, tacticity: atactic }
type: { forcefield: oplsaa2024 }
grow: { trials: 16, density: 0.60, seed: 1 }
relax: { method: lbfgs, fmax: 0.5 }
md: { ps: 1, temperature: 300, ensemble: nvt }
analyze: { properties: [density, rg] }
export: [lammps, gromacs]
Y
  cat > "$WORK/nr_cell.yaml" <<'Y'
recipe: 1
name: nr_cell
build:
  polymer: { smiles: "[*]C/C=C(C)\\C[*]", dp: 20, chains: 8 }
type: { forcefield: gaff2 }
grow: { trials: 16, density: 0.60, seed: 1 }
relax: { method: lbfgs, fmax: 0.5 }
md: { ps: 1, temperature: 300, ensemble: nvt }
analyze: { properties: [density, rg] }
export: [lammps]
Y
  for r in "$ROOT/samples/ps_cell.yaml" "$WORK/pe_cell.yaml" "$WORK/nr_cell.yaml"; do "$CAPS" run "$r" --out "$P" >/dev/null; done
  echo "$P"
}

make() {
  case "$1" in
    01_start) shot 01_start "" "wait=500" ;;
    02_studio) shot 02_studio "" "$PS" "relax=1" "module=8" "tab=0" "wait=1500" ;;
    03_molecule) shot 03_molecule "" "molecule=C1CCC(CC1)NSc1nc2ccccc2s1" "wait=3000" ;;
    # one atactic polystyrene chain (DP 8) previewed from the library
    04_polymer) shot 04_polymer "" "polydp=8" "polymer=P013" "wait=1000" ;;
    # PBSA (butylene succinate : butylene adipate 80:20 by moles) at DP 25, its composition table
    04b_composition) shot 04b_composition "" "polydp=25" "polymer=C112" "wait=1000" "pagescroll=end" ;;
    05_grow) shot 05_grow "" "$PS" "growwrap=1" "wait=1500" ;;
    06_field) shot 06_field "" "$PS" "module=7" "coverage=1" "fieldpick=0" "wait=1000" ;;
    07_relax) shot 07_relax "" "$PS" "relax=1" "wait=1000" ;;
    # 2,000 steps (2 ps) NPT at 300 K and 1 atm on the minimised cell, a frame every 250 steps: the curves and the report
    08_dynamics) shot 08_dynamics "" "$PS" "relax=1" "ensemble=2" "mdsteps=2000" "mdevery=250" "md=1" "wait=1500" "pagescroll=end" ;;
    # the cluster equilibration opened as a finished run (thermo log, printed output, structure), scrolled to the criteria
    09_equilibrate) ps_runs 09_equilibrate ps-equil || return 0
      cp "$PS_RUNS"/equil.csv "$PS_RUNS"/ps-equil.*.log "$PS_RUNS"/ps_eq.data "$DEMO/ps-equil/"
      (cd "$DEMO" && shot 09_equilibrate "" "wait=1500" "eqload=$DEMO/ps-equil" "pagescroll=1030") ;;
    # Analyze on the 2 ns NPT run at 300 K (201 frames): what a glass this short supports, against experiment; no Tg (that
    # would run a cooling here) and no D (the MSD is sub-diffusive)
    18_analyze) ps_runs 18_analyze ps-300K || return 0
      cp "$PS_RUNS"/ps_300K.data "$PS_RUNS"/ps_300K.lammpstrj "$DEMO/ps-300K/"
      (cd "$DEMO" && TOPO="$DEMO/ps-300K/ps_300K.data" shot 18_analyze "$DEMO/ps-300K/ps_300K.lammpstrj" "wait=1500" \
         "props=density,rdf,rg,cn,ced,delta,msd" "pair=1" "ref=ps-atactic" "analyze=1" "curvefor=Mean-square") ;;
    # sulfur S8 and the accelerator CBS packed into the free space of a natural-rubber cell
    10_pack) shot 10_pack "" "growcell=P056:gaff2-2.2.20-ambertools:8:20:0.6" "packinto=1" \
               "packmol=S1SSSSSSS1*12*sulfur_S8+C1CCC(CC1)NSc1nc2ccccc2s1*4*CBS" "pack=1" "wait=1500" ;;
    # sulfur cure: 20 H-S-S-H donors inserted into a natural-rubber cell (8 chains x 20 units, GAFF2), linked to
    # allylic carbons until 8 links join chains; 0.5 ps MD after each cycle, at most 8 cycles (4,000 steps)
    16_react) shot 16_react "" "growcell=P056:gaff2-2.2.20-ambertools:8:20:0.85" "rxset=2" "rxinsert=SS*20" "field=1" \
               "rxtarget=1:8" "rxcycles=8" "rxmd=0.5" "react=1" "wait=1500" ;;
    # natural rubber / cis-polybutadiene 70:30 by weight, DP 20, grown together into one cell
    11_blend) shot 11_blend "" "blend=1" "blendbuild=1" "wrap=1" "tab=3" "wait=1500" ;;
    # a poly(ethylene oxide) chain (DP 20) in TIP4P/2005 water with 0.15 mol/L NaCl
    21_solvation) shot 21_solvation "" "growcell=P014:gaff2-2.2.20-ambertools:1:20:0.1" "solvation=1" "solvate=1" "tab=3" "layer=1:eye" "wait=1500" ;;
    # a 30 x 30 A graphene sheet with carboxyl groups grafted on 5 % of its sites
    12_nano_sheet) shot 12_nano_sheet "" "nano=0" "sheet=30,30,1" "graft=carboxyl:0.05" "colour=0" "tab=3" "wait=1000" ;;
    # a (10,10) carbon nanotube with 10 cis-polyisoprene chains grown around it
    13_nano_matrix) shot 13_nano_matrix "" "nano=1" "nanomatrix=1" "nanobuild=1" "wrap=1" "wait=1500" ;;
    # a hydroxylated alpha-quartz (001) slab with a cis-polyisoprene film grown on it
    15_surface) shot 15_surface "" "surface=alpha-quartz" "surfbuild=1" "wrap=1" "wait=1500" ;;
    # the cluster cooling scans opened as finished runs, three seeds pooled (mean +/- SD); TG_RUNS names the run folders in
    # $PS_RUNS (default the 600 -> 300 K rerun; TG_RUNS="tg1 tg2 tg3" for the earlier 500 -> 250 K scans)
    22_glass) ps_runs 22_glass ps-tg || return 0
      local runs="${TG_RUNS:-tg600_1 tg600_2 tg600_3}" paths="" r
      for r in $runs; do
        if [ ! -d "$PS_RUNS/$r" ]; then echo "22_glass: no $PS_RUNS/$r yet; skipped" >&2; return 0; fi
        cp -R "$PS_RUNS/$r" "$DEMO/ps-tg/"
        paths="${paths:+$paths,}$DEMO/ps-tg/$r"
      done
      cp "$PS_RUNS"/*.log "$DEMO/ps-tg/" 2>/dev/null || true   # the runs' printed output: their schedules
      set -- $runs
      (cd "$DEMO" && shot 22_glass "$DEMO/ps-tg/$1/ps_tg.data" "wait=2500" "tgload=$paths") ;;
    # the project page on that folder, the polystyrene cell picked: its results and the methods text from its provenance
    20_project) local P; P="$(project_folder)"; shot 20_project "$P/ps_cell.data" "project=-" "wait=1500" ;;
    # Jobs after a grown, minimised and simulated cell: the last run (2,000 steps NPT) selected with its log and provenance
    23_jobs) shot 23_jobs "" "$PS" "relax=1" "ensemble=2" "mdsteps=2000" "md=1" "jobs=1" "wait=1500" ;;
    # DFT workbench: a 5 x 5 Ti3C2 sheet with OH on the fcc sites of both faces (one per surface Ti, Ti3C2(OH)2), built and
    # validated, seen from above at a tilt
    27_dft_sheets) shot 27_dft_sheets "" "sheetcell=5x5" "dft=sheets" "wait=1000" "pageview=0.25,0.75,1.7" "pagescroll=45" ;;
    # orthorhombic polyethylene (Pnam), a 2 x 3 x 2 supercell
    14_crystal) shot 14_crystal "" "crystal=1" "crystalsuper=2x3x2" "wait=1000" ;;
    # a Kremer-Grest melt of 50 chains x 100 beads, previewed in the builder
    17_cg) shot 17_cg "" "cg=1" "wait=1000" ;;
    # Export for the minimised OPLS-AA 2024 polystyrene cell: LAMMPS and GROMACS files, the checks
    19_export) shot 19_export "" "$PS" "relax=1" "exportpage=1" "wait=1000" ;;
    *) echo "no picture called $1" >&2; return 1 ;;
  esac
}

if [ $# -eq 0 ]; then set -- 01_start 02_studio 03_molecule 04_polymer 04b_composition 05_grow 06_field 07_relax 08_dynamics 09_equilibrate 10_pack 11_blend 12_nano_sheet 13_nano_matrix 14_crystal 15_surface 16_react 17_cg 18_analyze 19_export 20_project 21_solvation 22_glass 23_jobs 27_dft_sheets; fi
for n in "$@"; do make "$n"; done
echo "PNGs and work files in $WORK"
