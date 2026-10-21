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

# shot NAME FILE ARGS...: the Studio opened on FILE ("" for none), driven by ARGS, saved as NAME.jpg
shot() {
  local name="$1" file="$2"; shift 2
  echo "== $name"
  # a fresh settings folder (jobs, recent files, built molecules) for each picture; the user's own are never touched
  rm -rf "$WORK/state-$name"; mkdir -p "$WORK/state-$name"
  CAPS_SETTINGS="$WORK/state-$name/settings.json" CAPS_RECENT_DIR="$WORK/state-$name/recent" \
  dotnet "$DLL" --screenshot "$WORK/$name.png" "$file" "" "$@"
  sips -s format jpeg -s formatOptions 85 "$WORK/$name.png" --out "$IMG/$name.jpg" >/dev/null
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
    # the 21-step protocol at 0.002 of its published length (3.1 ps), then two 0.5 ps production blocks
    09_equilibrate) shot 09_equilibrate "" "$PS" "relax=1" "eqprotocol=0" "eqscale=0.002" "equntil=1" "eqblock=0.5" "eqmax=2" "eq=1" "wait=1500" ;;
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
    # a short cooling scan of a polyethylene cell (6 chains x 30 units, OPLS-AA 2024): 400 to 200 K in 50 K steps (five,
    # the fewest a two-line fit takes), 1 ps at each after 0.1 ps at 400 K, tau_P 100 fs, one replica (5,100 steps):
    # it shows the page at work, not a converged Tg
    22_glass) shot 22_glass "" "growcell=P001:oplsaa2024:6:30:0.85" "relax=1" "analysis=glass:run:400:200:50:1:0.1:1:100" "wait=1500" ;;
    # the project page on that folder, the polystyrene cell picked: its results and the methods text from its provenance
    20_project) local P; P="$(project_folder)"; shot 20_project "$P/ps_cell.data" "project=-" "wait=1500" ;;
    # Jobs after a grown, minimised and simulated cell: the last run (2,000 steps NPT) selected with its log and provenance
    23_jobs) shot 23_jobs "" "$PS" "relax=1" "ensemble=2" "mdsteps=2000" "md=1" "jobs=1" "wait=1500" ;;
    # DFT workbench: a Ti3C2 sheet with OH on the fcc sites of its top face, built and validated
    27_dft_sheets) shot 27_dft_sheets "" "dft=sheets" "wait=1000" ;;
    # orthorhombic polyethylene (Pnam), a 2 x 3 x 2 supercell
    14_crystal) shot 14_crystal "" "crystal=1" "crystalsuper=2x3x2" "wait=1000" ;;
    # a Kremer-Grest melt of 50 chains x 100 beads, previewed in the builder
    17_cg) shot 17_cg "" "cg=1" "wait=1000" ;;
    # Export for the minimised OPLS-AA 2024 polystyrene cell: LAMMPS and GROMACS files, the checks
    19_export) shot 19_export "" "$PS" "relax=1" "exportpage=1" "wait=1000" ;;
    *) echo "no picture called $1" >&2; return 1 ;;
  esac
}

if [ $# -eq 0 ]; then set -- 01_start 02_studio 03_molecule 04_polymer 04b_composition 05_grow 06_field 07_relax 08_dynamics 09_equilibrate 10_pack 11_blend 12_nano_sheet 13_nano_matrix 14_crystal 15_surface 16_react 17_cg 19_export 20_project 21_solvation 22_glass 23_jobs 27_dft_sheets; fi
for n in "$@"; do make "$n"; done
echo "PNGs and work files in $WORK"
