// caps — command-line front end over the same core the Studio uses.
#include <cstdlib>
#include <cctype>
#include <chrono>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <iostream>
#include <atomic>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

#include "caps/amber.hpp"
#include "caps/config.hpp"
#include "caps/analysis.hpp"
#include "caps/cbmc.hpp"
#include "caps/dlpoly.hpp"
#include "caps/bench.hpp"
#include "caps/dynamics.hpp"
#include "caps/edit.hpp"
#include "caps/elements.hpp"
#include "caps/ffdef.hpp"
#include "caps/ffio.hpp"
#include "caps/typing.hpp"
#include "caps/uff.hpp"
#include "caps/checks.hpp"
#include "caps/pipeline.hpp"
#include "caps/bundle.hpp"
#include "caps/crystal.hpp"
#include "caps/cg_commands.hpp"
#include "caps/cluster_job.hpp"
#include "caps/dft_commands.hpp"
#include "caps/spacegroup.hpp"
#include "caps/solvate.hpp"
#include "caps/nano.hpp"
#include "caps/properties.hpp"
#include "caps/equilibrate.hpp"
#include "caps/exchange_analysis.hpp"
#include "caps/exchange_moves.hpp"
#include "caps/exchange_library.hpp"
#include "caps/bond_react.hpp"
#include "caps/kremer_grest.hpp"
#include "caps/grow.hpp"
#include "caps/io.hpp"
#include "caps/martini3_molecules.hpp"
#include "caps/live.hpp"
#include "caps/mechanics.hpp"
#include "caps/molecule.hpp"
#include "caps/pack.hpp"
#include "caps/polymer.hpp"
#include "caps/react.hpp"
#include "caps/relax.hpp"
#include "caps/render.hpp"
#include "caps/ffcomplete.hpp"
#include "caps/mlpot.hpp"
#include "caps/update.hpp"
#include "caps/provenance.hpp"
#include "caps/recipe.hpp"
#include "caps/semicrystal.hpp"
#include "caps/tags.hpp"
#include "caps/thermal.hpp"
#include "caps/yaml.hpp"
#ifdef _WIN32
#include <io.h>
// the one Win32 call needed (no <windows.h>: its min/max macros)
extern "C" __declspec(dllimport) unsigned long __stdcall GetModuleFileNameA(void* module, char* name, unsigned long size);
#define isatty _isatty
#define fileno _fileno
#else
#include <unistd.h>
#endif
#ifdef __APPLE__
#include <mach-o/dyld.h>
#endif

using namespace caps;

namespace {

// Where the shipped data (data/, samples/) is: $CAPS_HOME, beside this program (Windows and Linux installs: data/ next to
// caps; the macOS app: Contents/Resources beside Contents/MacOS), the working directory, the source tree it was built from.
std::filesystem::path executable_dir() {
  namespace fs = std::filesystem;
  std::string p;
#if defined(_WIN32)
  char buf[4096];
  const unsigned long n = GetModuleFileNameA(nullptr, buf, sizeof buf);
  if (n > 0 && n < sizeof buf) p.assign(buf, n);
#elif defined(__APPLE__)
  char buf[4096];
  uint32_t size = sizeof buf;
  if (_NSGetExecutablePath(buf, &size) == 0) p = buf;
#else
  std::error_code link_ec;
  p = fs::read_symlink("/proc/self/exe", link_ec).string();
#endif
  if (p.empty()) return {};
  std::error_code ec;
  const fs::path c = fs::weakly_canonical(fs::path(p), ec);
  return (ec ? fs::path(p) : c).parent_path();
}

const std::vector<std::string>& caps_roots() {
  static const std::vector<std::string> roots = [] {
    std::vector<std::string> r;
    if (const char* h = std::getenv("CAPS_HOME")) r.push_back(h);
    const std::filesystem::path exe = executable_dir();
    if (!exe.empty()) {
      r.push_back(exe.string());
      r.push_back((exe.parent_path() / "Resources").string());
      r.push_back(exe.parent_path().string());
    }
    r.push_back(".");
    r.push_back(CAPS_SOURCE_ROOT);
    return r;
  }();
  return roots;
}

// A library id from a neutral short name its catalogue entry lists under "aliases" (spce, tip3p, martini …): the first
// entry naming it whose file is installed, so the name follows whichever library provides the model
std::string catalogue_alias(const std::string& root, const std::string& lower) {
  std::ifstream in(root + "/data/forcefields/catalogue.json");
  if (!in) return "";
  std::stringstream ss;
  ss << in.rdbuf();
  try {
    const Json cat = Json::parse(ss.str());
    for (const auto& e : cat["forcefields"].items()) {
      if (!e.has("aliases") || !e.has("file") || !e["file"].is_string()) continue;
      for (const auto& a : e["aliases"].items()) {
        std::string x = a.str();
        for (auto& c : x) c = char(std::tolower((unsigned char)c));
        if (x == lower && std::filesystem::exists(root + "/data/forcefields/" + e["file"].str())) return e["file"].str();
      }
    }
  } catch (const std::exception&) {}
  return "";
}

// A force field by path, or by its library name (pcff, compass, gaff2, oplsaa2024 …) from data/forcefields.
std::string ff_path(const std::string& s0) {
  if (s0.empty() || std::filesystem::exists(s0)) return s0;
  // the short names of the libraries rebuilt from their own sources (as caps.core._forcefield_path in Python)
  static const std::map<std::string, std::string> alias = {
      {"gaff2", "gaff2-2.2.20-ambertools"}, {"gaff", "gaff-1.81-ambertools"}, {"opls", "oplsaa2024"}, {"pcff", "pcff-frc"},
      {"opls-ua", "oplsua2024"}, {"oplsua", "oplsua2024"}, {"compass", "compass-frc"}, {"cvff", "cvff-frc"},
      {"dreiding", "dreiding1990"}, {"cgenff", "cgenff-toppar"},
      {"trappe-ua", "trappe-ua-umn"}, {"trappe", "trappe-ua-umn"},
      {"clayff", "clayff-cygan2004"}, {"clayff2004", "clayff-cygan2004"}, {"trappe-eh", "trappe-eh-umn"},
      // water models and the CG lipid model by their usual names (the libraries built from their papers)
      {"spc", "spc-berendsen1981"}, {"spce", "spce-berendsen1987"}, {"tip3p", "tip3p-jorgensen1983"},
      {"tip3p-charmm", "tip3p-charmm-mackerell1998"}, {"tip4p", "tip4p-jorgensen1983"}, {"tip4p2005", "tip4p2005-abascal2005"},
      {"tip5p", "tip5p-mahoney2000"}, {"opc", "opc-izadi2014"}, {"opc3", "opc3-izadi2016"}, {"mw", "mw-molinero2009"},
      {"cooke-deserno", "cooke-deserno2005"}};
  std::string lower = s0;
  for (auto& c : lower) c = char(std::tolower((unsigned char)c));
  auto al = alias.find(lower);
  // a short name whose first library is built only on the user's machine: the published one when it is absent
  // (cgenff: the local toppar build, else CGenFF 2a5 from the 2010 paper's SI, in the repository)
  static const std::map<std::string, std::string> fallback = {{"cgenff", "cgenff-2a5-jcc2010"}};
  auto fb = fallback.find(lower);
  std::vector<std::string> tries = {s0};
  if (al != alias.end()) tries.push_back(al->second);
  if (fb != fallback.end()) tries.push_back(fb->second);
  for (const std::string& s : tries)
    for (const auto& root : caps_roots())
      for (const std::string& name : {s, s + ".json"})
        if (std::filesystem::exists(root + "/data/forcefields/" + name)) {
          if (s != s0 && fb != fallback.end()) {   // say which one a short name with a fallback took (once)
            static std::set<std::string> said;
            if (said.insert(lower).second)
              std::fprintf(stderr, "force field %s: %s%s\n", s0.c_str(), s.c_str(),
                           (s == fb->second && al != alias.end()
                                ? " (" + al->second + " is not built here: bench/ff/open/build_charmm.sh builds it)"
                                : std::string("")).c_str());
          }
          return root + "/data/forcefields/" + name;
        }
  for (const auto& root : caps_roots())
    if (const std::string f = catalogue_alias(root, lower); !f.empty()) return root + "/data/forcefields/" + f;
  return s0;
}

// Charges worth keeping: a file whose charge column is all zero (a builder's placeholder) has none.
bool file_charges(const System& s) {
  if (!s.has_charges) return false;
  for (const auto& a : s.atoms)
    if (std::fabs(a.charge) > 1e-9) return true;
  return false;
}

int usage() {
  std::cerr << "caps 0.1.0\n"
               "usage:\n"
               "  caps --version | --check-update   this build; whether a newer CAPS release is out (exit 0 no, 1 yes, 2 cannot tell)\n"
               "  caps info    FILE [--topology DATA]\n"
               "  caps render  FILE -o OUT.png|OUT.svg [--topology DATA] [--frame N] [--size WxH]\n"
               "               [--bg dark|white|transparent] [--colour molecule|element|type|distance]\n"
               "               [--style ball|space|sticks|noh|backbone] [--yaw DEG] [--pitch DEG] [--zoom Z] [--no-cell] [--wrap]\n"
               "  caps shape   FILE [--topology DATA]          per-molecule Rg and shape\n"
               "  caps rdf     FILE [--topology DATA] [--rmax 12] [--dr 0.2] [--pair C-C] [--inter]\n"
               "  caps analyze FILE [--topology DATA] [--props density,rdf,sq,xray,neutron,rg,ree,cn,persistence,msd,diffusion,\n"
               "               relaxation,ced,delta,ffv,psd,crosslinks,entanglements] [--first N --last N --stride N] [--frame-ps X | --timestep-fs 1]\n"
               "               [--pair C-C --inter] [--qmax 25 --dq 0.02 --qdirect 4] [--probe 0] [--grid 0.4] [--ff FF.json] [--json OUT] [--csv DIR]\n"
               "               diffusion: [--fit 0.2:0.5] [--viscosity ETA (mPa·s) --temp T (Yeh–Hummer finite-size correction)]\n"
               "               interfaces (zprofile, adhesion, interaction, orientation): [--surface 1-3,7 (molecule ids; 1)] [--axis x|y|z (z)] [--zbin 0.5]\n"
               "  caps elastic FILE [--topology DATA] [--method strain|fluct|fluct-run] [--configs N] [--strain 1e-4] [--temp T] [--ps 100] [--ff FF.json] [--json OUT]\n"
               "  caps tensile DATA -o OUT.data [--axis x] [--rate 1e-3] [--strain 0.2] [--temp 300] [--fixed-lateral] [--precision mixed|double] [--ff FF.json [--overlay OV.json] [--charges MODE] [--complete ml]] [--csv DIR]\n"
               "  caps tg DATA -o OUT.data [--from 500 --to 200 --step 20 --ps 100] [--precision mixed|double] [--ff FF.json [--overlay OV.json] [--charges MODE] [--complete ml]] [--csv DIR]   |   caps tg --fit RUN [RUN ...] [--volume]\n"
               "  caps exchange arrhenius TABLE [TABLE2] [--reduced] [--from T --to T]   |   caps exchange modulus PTENSOR --cell-volume A3 --temp K --dt FS [--fit PS] [--g-inf MPa] [--corr-ps PS]\n"
               "               | caps exchange memory PARTNERS --frame-ps PS   |   caps exchange bridges START NOW [--axis z] [--plane 0.5] [--bin 1]\n"
               "               | caps exchange heal HEALED INTACT [--modulus-strain 0.03] [--at-strain E] [--up-to E]   [--json OUT] [--csv DIR]   (bond exchange and healing)\n"
               "               | caps exchange list   |   caps exchange show CHEMISTRY   (exchange chemistries and their barriers, data/exchange)\n"
               "               | caps exchange run DATA -o OUT [--chemistry ID] --kind swap|shift|bind [--sites SEL] [--partners SEL] [--labels chain|ends] [--capture 1.3]\n"
               "                 [--barrier E] [--binding E] [--every 10] [--attempts 1] [--steps 10000] [--temp 300] [--dt 1] [--kg [--k-theta K] | --ff FF.json]\n"
               "               | caps exchange bond-react DATA -o STEM --chemistry ID [--temp 300] [--t-ref T] [--steps N] [--ff FF.json]   (LAMMPS fix bond/react)\n"
               "               | caps exchange bond-swap DATA -o STEM [--sites SEL] [--stickers | --both-ends] [--every 50] [--fraction 0.5] [--cutoff 1.3] [--seed 1] [--steps N] [--k-theta K]\n"
               "  caps frames  FILE OUT.lammpstrj|dcd|xyz|pdb|gro|trr [--topology DATA] [--first N] [--last N] [--stride N]   (read thinned)\n"
               "  caps convert FILE OUT.data|xyz|pdb|mol2|car|gro|sdf|cif|vasp|POSCAR [--topology DATA]   (.car: with its .mdf)\n"
               "  caps edit    FILE --ops 'OP ARGS; OP ARGS …' | --ops-file OPS.txt -o OUT   edits in order, atoms numbered from 1:\n"
               "               element SEL Sym · delete SEL · bond I J [order] · unbond I J · addh [SEL] · attach I SMILES ·\n"
               "               length I J Å · angle I J K ° · torsion I J K L ° · invert I · config I R|S · rotate SEL x,y,z ° ·\n"
               "               mirror SEL nx,ny,nz · move SEL dx,dy,dz · clean [SEL] · tacticity iso|syndio ·\n"
               "               crystals: supercell na nb nc · primitive · niggli · conventional · redefine m11 … m33 · vacuum Å ·\n"
               "               nanowire u v w radius [repeats] [cylinder|hexagonal|square] [vacuum] ·\n"
               "               SEL: 3,5-9 · all · element:C,N · smarts:PATTERN · type:LABEL (numbers after a delete shift)\n"
               "  caps provenance FILE [--json | --bibtex | --methods] [--compare OTHER]   the steps that produced FILE (FILE.provenance.json)\n"
               "  caps bench   [T1 T2 … | --all] [--repeats 3] [--quick] [--out DIR] [--samples DIR]   the built-in validation suite\n"
               "  caps build   SMILES -o OUT.mol2|OUT.pdb|OUT.xyz|OUT.data [--conformers 1] [--seed 1] [--ff FF.json] [--all]\n"
               "               a 3D molecule from SMILES; --ff cleans each conformer up with that force field (with typing rules)\n"
               "  caps check   FILE [--topology DATA] [--report OUT.md]   file checks (counts, bonds, contacts, charges, cell)\n"
               "  caps pipeline FILE [--topology DATA] --steps STEPS.json|STEPS.yaml|'[…]' [--frame N] [--table NAME] [--particles EXPR] [--out DIR] [--branch NAME]\n"
               "                                   visualize pipeline on one frame: step status, attributes, a table as CSV\n"
               "  DFT surfaces & adsorption (caps <command> --help for every option, its default and why):\n"
               "  caps sheet [--preset Ti3C2 | --from FILE --formula Ti3C2 --remove Al] -o OUT     a 2D sheet\n"
               "  caps terminate SHEET --top O:0.5,OH:0.25,F:0.25 [--supercell 3x3] [--vasp-set DIR] -o OUT\n"
               "  caps validate FILE… [--expect Ti3C2O2]       2D slab checks (exit 1 on FAIL)\n"
               "  caps adsorb-dft SLAB [--from-relaxed] [--supercell 5x5] [--smiles S] --out ROOT\n"
               "  caps vasp-set STRUCTURE --out DIR [--pp-dir PAW] [--profile slurm-workspace]\n"
               "  caps vasp-relax CRYSTAL --out DIR [--dimensionality 3d|2d] [--encuts 900,1350,900,900] [--account A]   relax1 … relax4 (guide ch. 5)\n"
               "  caps phono3py setup|estimate CRYSTAL [--out DIR] [--mode full|random|both] [--dimensionality 3d|2d] · phono3py check|analyze DIR   thermal conductivity\n"
               "  caps vasp-conv setup|collect CASE · vasp-scan make|fit · vasp-derived charge|cdd|freq|aimd\n"
               "  caps vasp-jobs submit|update|reset|cleanup|store · vasp-check · vasp-progress · vasp-health · vasp-bind\n"
               "  caps vasp-analyze geom|wf|dos|cdd|bader|freq|md|summary …\n"
               "  Cluster jobs (caps job help): CAPS's own runs as batch jobs on any SLURM or PBS cluster, or a workstation\n"
               "  caps job profile [--preset slurm|slurm-workspace|pbs|workstation] [--set key=value,…]   your host profile\n"
               "  caps job new --title T --kind md [--input FILES] [--cpus N] [--time HH:MM:SS] [--submit] -- caps md …\n"
               "  caps job status|tail [-f]|poll|collect|cancel|resume DIR · caps job list · caps job new --array LIST.txt\n"
               "  Coarse-graining (caps <command> --help for every option):\n"
               "  caps cgmap STRUCTURE… --preset ester-cut -o DIR [--dump DUMP]    chemistry-aware coarse-grained mapping\n"
               "  caps cgfit bonded|refine|targets|ibi-start|ibi-step|ibi-run|fit|tg|calibrate|sample-chain …   CG potentials\n"
               "  caps cgbuild --units BS=B+S --dp 100 --chains 50 --bonded DIR --maps MAP -o DIR   a CG melt of repeat units\n"
               "  caps ppa MAP FRAMES [--method caps|z1|lammps] -o DIR      entanglements: primitive paths, N_e estimators\n"
               "  caps mech decks|analyze …       tension decks and their analysis · caps cgdyn MAP DUMP -o DIR  chain dynamics\n"
               "  caps backmap REF_MAP REF_DATA --cg CG_MAP --frame CG_FRAME -o DIR   CG beads back to all atoms\n"
               "  caps bundle  FILE [--topology DATA] --steps S.json [-o OUT.caps-bundle.zip] [--include-input] [--frame N]\n"
               "                                   a figure with its data, pipeline, provenance and hashes (and the input)\n"
               "  caps reproduce BUNDLE.caps-bundle.zip   rebuild a bundle's data from its input and pipeline, compare sha256\n"
               "  caps run     RECIPE.yaml|json [--seed N] [--threads N] [--out DIR] [--json] [--log thermo.csv] [--dump traj.lammpstrj]   build → type → grow → relax → equilibrate →\n"
               "               md → analyze → export from one file (exit 0 ok · 2 input · 3 missing params · 4 failed run)\n"
               "  caps run     PIPELINE.yaml|json [--input 'runs/*/X.lammpstrj'] [--frame first|last] [--csv OUT] [--out DIR: the outputs: block] [--branch NAME]\n"
               "                                   a saved pipeline over many inputs: one row of attributes per input\n"
               "  caps crystal --group 'P 42/m n m' --cell a,b,c[,α,β,γ] --sites 'Ti1 Ti 0 0 0; O1 O 0.3048 0.3048 0' -o OUT\n"
               "               (OUT.cif in P 1, OUT.vasp or POSCAR for VASP, .data, .pdb, .xyz …)\n"
               "               [--supercell 2,2,2] [--primitive] [--symmetrize] [--tolerance 0.01]   a crystal from a space group\n"
               "  caps crystal CRYSTAL.cif --find-symmetry [--tolerance 0.1]   its space group and asymmetric unit\n"
               "  caps crystal --groups [QUERY]    the 530 space-group settings (key, number, Hermann–Mauguin, Hall)\n"
               "  caps solvate [SOLUTE] -o OUT [--edge 30 | --box a,b,c | --padding 10] [--solvent water|toluene|…] [--model TIP4P/2005]\n"
               "               [--salt NaCl --conc 0.15 | --neutralise | --ions 3,6 | --no-ions] [--molecules N] [--tolerance 2] [--solvents]\n"
               "                                   solvent and ions packed around a solute (CAPS Pack)\n"
               "  caps surface CRYSTAL.cif -o OUT.data|mol2|pdb|xyz [--hkl 0,0,1] [--layers 3] [--termination 1] [--vacuum 15]\n"
               "               [--supercell 2,2] [--no-orthogonal] [--max-strain 2] [--passivate] [--list]   a slab (terminations listed)\n"
               "  caps interface CRYSTAL.cif|SLAB -o OUT --units SMILES[,…] [surface options] [--film 30] [--film-density 0.9]\n"
               "               [--chains N] [--dp 10] [--gap 1] [--vacuum 0] [--sequence …] [--ff FF]   a polymer film on a surface\n"
               "  caps nano    tube [--n 10 --m 10 --length 25 --finite] | sheet [--lx 20 --ly 20 --layers 1 --flake] |\n"
               "               particle CRYSTAL.cif [--shape sphere|cube|octahedron|cuboctahedron|truncated-octahedron|icosahedron|rod|cone|frustum|tetrahedron|pyramid|hemisphere|fibre\n               --radius 12 --height 24 --top-ratio 0.5 --length 20 --passivate]\n"
               "               [--units SMILES --chains 10 --dp 20 --density 0.9]   -o OUT   fillers, alone or in a polymer matrix\n"
               "  caps martini FILE -o CG.data [--itp CG.itp] [--centre geometric|mass]   Martini 3 small molecules mapped from\n"
               "               an all-atom structure (vermouth's mappings); water, ions and other molecules left out\n"
               "  caps pore    slit [--width 10 --layers 1 --lx 26 --ly 22 --vacuum] | cylinder CRYSTAL.cif [--width 14 --wall 6 --length 20 --passivate] |\n"
               "               framework CRYSTAL.cif [--supercell 2,2,2]   [--fluid SMILES --count N --tolerance 2 --seed 1] -o OUT   a fluid in a pore\n"
               "  caps pull    FILE [--normal | --axis x|y|z] [--distance 10] [--rate 5] [--spring 10] [--temp 300] [--surface 1] [--csv OUT]\n"
               "               pull-out / debonding of a film from a held surface: interfacial shear strength, work of separation\n"
               "  caps blend   --components SMILES1,SMILES2 [--weights 0.5,0.5] [--chains 8] [--dp 20] [--density 0.5] [--slabs | --droplet] -o OUT\n"
               "  caps semicrystalline -o OUT.data [--polymer pe|pe-orthorhombic|CELL.cif] [--repeats 4,6] [--crystal 5 (nm) | --crystal-cells 20]\n"
               "               [--amorphous 4 (nm)] [--loops 0.5 --bridges 0.2 --tails 0.3] [--reentry adjacent|random] [--density 0.855]\n"
               "               [--seed 1] [--ff oplsaa2024|oplsua2024] [--no-relax] [--report R.json] [--lammps-input IN.lmp] [--gromacs STEM]\n"
               "               a lamellar stack: crystalline stems, loops, bridges and tails in the amorphous layer (groups crystal, amorphous)\n"
               "  caps semicrystalline FILE --analyse [--report R.json]   loops, bridges, tails, P2 crystallinity, profiles of any frame\n"
               "  caps grow    -o OUT.data|OUT.pdb|OUT.xyz [--chains 10] [--dp 8] [--density 0.5 | --box 33]\n"
               "               [--tacticity atactic|isotactic|syndiotactic] [--seed 1] [--trans] [--scale 1.0]\n"
               "               [--units '*CC(*)c1ccccc1,*CC(*)(C)C(=O)OC' --sequence homopolymer|alternating|block|random|shuffled|gradient|pattern\n"
               "                --weights 0.7,0.3 --blocks 20,20 --pattern AAB --ff FF.json --method trials|rosenbluth|rosenbluth_lj\n"
               "                --threads N --lookahead 1]   any repeat units (else polystyrene); the contact scale is automatic\n"
               "                (lowered only where a chain is stuck) unless --scale fixes it\n"
               "  caps field   FILE [--topology DATA] [--forces OUT.txt]   the structure's force field (below): types, terms, energy (and per-atom forces)\n"
               "  caps relax   FILE -o OUT.data|OUT.pdb|OUT.xyz [--method lbfgs|cg|sd|fire] [--ftol 0.5] [--iterations 5000]\n"
               "               [--density 1.05] [--step 0.06] [--box-relax] [--pressure 1] [--no-pushoff] [--cutoff Å]\n"
               "               [--no-coulomb] [--ff FF.json [--overlay OV.json[,MORE.json]] [--charges auto|types|gasteiger|keep|qeq|increments] [--complete ml]] [--quiet]\n"
               "               (.data output carries the force field for LAMMPS)\n"
               "               without --ff (relax, md, equilibrate, field): the force field a data file CAPS wrote holds (its coefficients,\n"
               "               charges and 1-4 weights; refused when CAPS cannot rebuild it from the file), else GAFF (C and H) or UFF;\n"
               "               without --cutoff (relax, md, equilibrate): the force field's own cut-off, else 10 Å\n"
               "  caps cbmc    FILE -o OUT.data [--moves 1000] [--trials 8] [--max-torsions 4] [--temp 300] [--cutoff 9] [--seed 1]\n"
               "               [--ff FF] [--no-coulomb]   configurational-bias regrowth of chain ends (Siepmann & Frenkel)\n"
               "  caps thermal FILE --dimensionality 3d|2d [--method gk|nemd] [--temp 300] [--eq-ps 10] [--ps 100] [--dt 1] [--every 5]\n"
               "               [--corr-ps 5] [--axis x|y|z] [--layers 20] [--swap-every 100] [--transient 0.3] [--thickness 3.35] [--seed 1]\n"
               "               [--ff FF] [--json OUT.json] [--lammps-input IN.lmp]   thermal conductivity by MD: Green-Kubo (heat-flux\n"
               "               autocorrelation, centroid form for bonded terms) or Muller-Plathe NEMD; 2d: a sheet with vacuum along z,\n"
               "               in-plane only, kappa over an effective thickness and the sheet conductance\n"
               "  caps md      FILE -o OUT.data [--steps 10000] [--dt 1] [--temp 300] [--thermostat bussi|langevin|nose-hoover|none]\n"
               "               [--tau-t 100] [--barostat none|crescale|berendsen|mtk] [--pressure 1] [--tau-p 1000]\n"
               "               [--constraints none|h-bonds|all-bonds] [--constraint-solver shake|lincs] [--hmr 3.024] [--seed 1] [--new-velocities] [--thermo 100] [--dump TRAJ.lammpstrj --every 1000]\n"
               "               [--log thermo.csv] [--cutoff Å] [--skin 1.5] [--threads N] [--no-coulomb] [--precision mixed|double] [--ff FF.json [--overlay OV.json] [--charges MODE] [--complete ml]] [--quiet]\n"
               "               [--checkpoint-every N [--checkpoint OUT.restart.data]]   the full state every N steps\n"
               "               [--progress-file progress.jsonl]   one JSON line per thermo row (default in a SLURM job)\n"
               "               caps md OUT.restart.data --resume --steps TOTAL -o OUT   continues from a checkpoint\n"
               "  caps equilibrate FILE -o OUT.data [--protocol larsen21|annealing|pushoff|PROTOCOL.txt] [--print-protocol]\n"
               "               [--tfinal 300] [--tmax 600] [--pfinal 1] [--pmax 49346] [--scale 1] (atm, K; --scale shortens every stage)\n"
               "               [--cycles 3] [--tlow 300] [--thigh 600] [--ramp 50] [--hold 50]   (annealing, ps)\n"
               "               [--until-converged] [--block 20] [--max-blocks 20] [--dt 1] [--thermostat bussi|langevin|nose-hoover] [--tau-t 100]\n"
               "               [--barostat crescale|berendsen|mtk] [--tau-p 1000] [--seed 1] [--dump TRAJ.lammpstrj --every-ps 10] [--log thermo.csv]\n"
               "               [--precision mixed|double] [--constraints h-bonds] [--hmr 3.024] [--cutoff Å] [--ff FF.json [--overlay OV.json] [--charges MODE] [--complete ml]] [--quiet]\n"
               "               (mixed: pair terms in single precision\n"
               "               with SIMD, as GROMACS; --hmr: hydrogens at that mass, taken from their bonded atoms, for steps past 2 fs)\n"
               "  caps chains  FILE [--topology DATA]          backbones and mean-square internal distances\n"
               "  caps pack    INPUT.inp [-o OUT] [--threads N] [--quiet]   packmol-style input (structure … end structure)\n"
               "  caps contacts FILE [--tol 2.0] [--no-pbc] [--molecule-size N]   smallest distance between atoms of different molecules\n"
               "  caps react   FILE -o OUT [--template NAME|FILE.txt ...] [--capture Å] [--per-cycle 10] [--cycles 100] [--target 1.0]\n"
               "               [--no-relax] [--relax-iterations 500] [--md-ps 0] [--during-md] [--temp 300] [--seed 1] [--cutoff Å] [--no-coulomb]\n"
               "               [--ff FF.json [--overlay OV.json] [--charges MODE] [--complete ml]] [--between-chains] [--crosslinks links|per-chain|density|mc|degree:VALUE]\n"
               "               [--sites-per-chain N] [--weights 2,1] [--auto-capture [--capture-max 8] [--capture-step 0.5]] [--keep-byproducts]\n"
               "               [--insert SMILES --count N] [--remove-unreacted] [--summary NET.json] [--fa 2 --fb 4 --ratio 1] [--list-templates] [--quiet]\n"
               "  caps ml eval FILE [--model M.model] [--json OUT.json]   energy (eV) and forces (eV/A) of a molecule from CAPS's\n"
               "               machine-learned potential (data/ml/caps-mace-v1.model; a build with -DCAPS_WITH_TORCH=ON)\n"
               "  caps ml relax JOBS.json -o OUT.json [--model M.model]   constrained relaxations (bonds / angles / dihedrals held)\n"
               "  caps ff complete FILE --ff FF.json [--overlay OV.json] -o OVERLAY.json [--report R.json] [--model M.model] [--step 15]\n"
               "               the bonded terms FF lacks for FILE, found and fitted with CAPS's machine-learned potential\n"
               "  caps ff import-lt FILE.lt -o FF.json          convert a moltemplate force field into the CAPS format\n"
               "  caps ff info FF.json                           types, rules, styles, references\n"
               "  caps ff overlays [FF_ID]                       the library's overlays (parameter refits) for --overlay\n"
               "  caps ff type FILE --ff FF.json [--typing RULES.json] [-o TYPES.txt] [--explain]   assign atom types from SMARTS rules\n"
               "  caps ff apply FILE --ff FF.json [-o OUT.data [--lammps-input OUT.in [--lammps-run check|minimize|nvt|npt --temp 300 --press 1 --dt FS (default: the force field's, Martini 20, else 0.5) --steps N] [--moltemplate SYSTEM.lt]]] [--gromacs STEM] [--dlpoly DIR] [--overlay USER.json[,MORE.json]] [--types TYPES.txt] [--charges auto|keep|types|gasteiger]\n"
               "               [--complete ml [--model M.model]]   (terms the force field lacks fitted by CAPS ML, marked as such in every file)\n"
               "               [--lammps-style native|exact] [--hybrid] [--kspace auto|pppm|ewald|dsf|cut] [--kspace-accuracy 1e-4] [--lammps-cutoff Å] [--units auto|real|metal]\n"
               "               [--list] [-o OUT.data]   parameters for a structure whose atoms carry type names (or TYPES.txt)\n";
  return 2;
}

// Files matching a pattern with * and ? in any path segment (cells/seed*/PS_melt.data), sorted.
bool wild_match(const char* p, const char* t) {
  if (!*p) return !*t;
  if (*p == '*') return wild_match(p + 1, t) || (*t && wild_match(p, t + 1));
  return *t && (*p == '?' || *p == *t) && wild_match(p + 1, t + 1);
}
std::vector<std::string> glob_files(const std::string& pattern) {
  namespace fs = std::filesystem;
  if (pattern.find_first_of("*?") == std::string::npos) return fs::exists(pattern) ? std::vector<std::string>{pattern} : std::vector<std::string>{};
  const fs::path full = fs::absolute(pattern);
  std::vector<fs::path> current = {full.root_path()};
  std::vector<std::string> parts;
  for (const auto& part : full.relative_path()) parts.push_back(part.string());
  for (size_t k = 0; k < parts.size(); ++k) {
    std::vector<fs::path> next;
    const bool last = k + 1 == parts.size();
    for (const auto& dir : current) {
      std::error_code ec;
      if (parts[k].find_first_of("*?") == std::string::npos) {
        const fs::path p = dir / parts[k];
        if (last ? fs::is_regular_file(p, ec) : fs::is_directory(p, ec)) next.push_back(p);
        continue;
      }
      for (const auto& e : fs::directory_iterator(dir, ec))
        if ((last ? e.is_regular_file() : e.is_directory()) && wild_match(parts[k].c_str(), e.path().filename().string().c_str())) next.push_back(e.path());
    }
    std::sort(next.begin(), next.end());
    current = std::move(next);
  }
  std::vector<std::string> out;
  for (const auto& p : current) out.push_back(p.string());
  return out;
}

// Every option any command reads (a typo such as --cutof must not run with the default): keep in step with the
// options the commands look up.
const std::set<std::string>& known_options() {
  static const std::set<std::string> k = {
    "--all", "--allow-missing", "--auto-capture", "--auto-scale", "--between-chains", "--capture-max", "--capture-step", "--crosslinks",
    "--during-md", "--keep-byproducts", "--relax-iterations", "--remove-unreacted", "--sites-per-chain", "--summary", "--axis", "--barostat", "--beads", "--bg", "--bibtex",
    "--block", "--blocks", "--born", "--branch", "--born-every", "--born-strain", "--box", "--box-relax",
    "--capture", "--cell", "--centre", "--chains", "--charges", "--colour", "--comfortable", "--compare", "--complete",
    "--components", "--conc", "--constraint-solver", "--constraints", "--configs", "--conformers", "--count", "--csv", "--cutoff", "--cycles",
    "--density", "--deterministic", "--distance", "--dp", "--dq", "--dr", "--droplet", "--dt", "--dump",
    "--edge", "--eq-ps", "--equilibrate", "--dlpoly",
    "--escalate", "--every", "--every-ps", "--ewald-rtol", "--exclude-mol", "--explain",
    "--fa", "--fb", "--ff", "--film", "--film-density", "--find-symmetry", "--finite", "--first",
    "--amber", "--digits", "--dsf-alpha", "--fit", "--fix-mol", "--fixed-lateral", "--flake", "--fluid", "--forcefields", "--forces", "--frame",
    "--frame-ps", "--from", "--ftol", "--gap", "--grid", "--gromacs", "--group", "--groups", "--hkl", "--hmr",
    "--hold", "--hybrid", "--include-input", "--input", "--insert", "--inter", "--ions", "--itp", "--iterations",
    "--json", "--kspace", "--kspace-accuracy", "--lammps-cutoff", "--units", "--lammps-input", "--lammps-run", "--moltemplate",
    "--lammps-style", "--last", "--layers", "--length", "--list", "--list-templates", "--log", "--lx", "--ly", "--m",
    "--lookahead", "--martini", "--max-blocks", "--max-torsions", "--moves", "--max-strain", "--md-ps", "--method", "--methods", "--model", "--molecule-size",
    "--molecules", "--n", "--names", "--neutralise", "--new-velocities", "--no-cell", "--wrap",
    "--no-coulomb", "--no-ions", "--no-orthogonal", "--no-pbc", "--no-pushoff", "--no-relax",
    "--no-tail", "--normal", "--out", "--overlay", "--padding", "--pair", "--particles",
    "--ops", "--ops-file", "--checkpoint-every", "--checkpoint", "--resume", "--resume-at", "--bench-json", "--resume-points", "--resume-curve", "--height", "--top-ratio", "--passivate", "--pattern", "--per-cycle", "--perspective", "--pfinal", "--pitch", "--pmax", "--pme",
    "--pme-order", "--pme-spacing", "--precision", "--press", "--pressure", "--primitive", "--print-protocol", "--probe", "--progress-file",
    "--props", "--protocol", "--ps", "--qdirect", "--qmax", "--quick", "--quiet", "--radius", "--ramp", "--rate",
    "--ratio", "--batch", "--repeats", "--report", "--rmax", "--salt", "--samples", "--scale", "--seed", "--sequence",
    "--shape", "--sites", "--size", "--skin", "--slabs", "--solvent", "--solvents", "--spring", "--step",
    "--steps", "--strain", "--stride", "--style", "--supercell", "--surface",
    "--symmetrize", "--table", "--tacticity", "--target", "--tau-p", "--tau-t", "--temp", "--template",
    "--termination", "--tfinal", "--thermo", "--thermostat", "--thigh", "--threads", "--timestep-fs", "--tlow",
    "--tmax", "--to", "--tol", "--tolerance", "--topology", "--trans", "--trials", "--types", "--typing", "--units",
    "--polymer", "--crystal", "--crystal-cells", "--amorphous", "--loops", "--bridges", "--tails", "--reentry", "--analyse",
    "--corr-ps", "--long", "--swap-every", "--transient", "--dimensionality", "--thickness",
    "--until-converged", "--vacuum", "--viscosity", "--volume", "--cell-volume", "--reduced", "--kind", "--chemistry", "--t-ref", "--barrier", "--binding", "--attempts", "--partners", "--labels", "--kg", "--k-theta", "--stickers", "--both-ends", "--fraction", "--g-inf", "--plane", "--bin", "--modulus-strain", "--at-strain", "--up-to", "--wall", "--weights", "--width", "--yaw", "--zbin", "--zoom"};
  return k;
}

// The known option closest to a mistyped one (edit distance), for the error message.
std::string closest_option(const std::string& a) {
  std::string best;
  size_t bd = 99;
  for (const auto& k : known_options()) {
    std::vector<size_t> d(k.size() + 1);
    for (size_t j = 0; j <= k.size(); ++j) d[j] = j;
    for (size_t i = 1; i <= a.size(); ++i) {
      size_t prev = d[0];
      d[0] = i;
      for (size_t j = 1; j <= k.size(); ++j) {
        const size_t t = d[j];
        d[j] = std::min({d[j] + 1, d[j - 1] + 1, prev + (a[i - 1] == k[j - 1] ? 0 : 1)});
        prev = t;
      }
    }
    if (d[k.size()] < bd) bd = d[k.size()], best = k;
  }
  return bd <= 3 ? best : std::string();
}

// Options without a value (the rest take the next word)
bool is_cli_switch(const std::string& a) {
  return a == "--no-cell" || a == "--wrap" || a == "--inter" || a == "--perspective" || a == "--trans" || a == "--escalate" ||
         a == "--box-relax" || a == "--no-pushoff" || a == "--no-coulomb" || a == "--quiet" || a == "--new-velocities" ||
         a == "--until-converged" || a == "--print-protocol" || a == "--no-pbc" ||
         a == "--no-relax" || a == "--list-templates" || a == "--list" || a == "--allow-missing" || a == "--no-tail" || a == "--explain" || a == "--names" || a == "--fixed-lateral" || a == "--volume" || a == "--quick" || a == "--all" || a == "--pme" || a == "--no-orthogonal" || a == "--passivate" || a == "--auto-scale" || a == "--finite" || a == "--flake" || a == "--normal" || a == "--slabs" || a == "--droplet" || a == "--include-input" || a == "--primitive" || a == "--symmetrize" || a == "--find-symmetry" || a == "--groups" || a == "--neutralise" || a == "--no-ions" || a == "--solvents" || a == "--bibtex" || a == "--json" || a == "--deterministic" || a == "--vacuum" || a == "--methods" ||
         a == "--hybrid" || a == "--resume" || a == "--between-chains" || a == "--keep-byproducts" || a == "--auto-capture" ||
         a == "--during-md" || a == "--remove-unreacted" || a == "--analyse" || a == "--long" || a == "--reduced" || a == "--kg" || a == "--stickers" || a == "--both-ends";
}

std::map<std::string, std::string> parse(int argc, char** argv, int from, std::vector<std::string>& pos) {
  std::map<std::string, std::string> o;
  for (int i = from; i < argc; ++i) {
    std::string a = argv[i];
    if (a.rfind("--", 0) == 0 && !known_options().count(a)) {
      const std::string near = closest_option(a);
      throw std::invalid_argument("unknown option " + a + (near.empty() ? "" : " (did you mean " + near + "?)"));
    }
    if (a.rfind("--", 0) == 0 || a == "-o") {
      // --json is a switch (provenance, run, job) except when a .json file name follows (analyze, elastic: --json OUT.json)
      const auto next_json = [&] { const std::string n = i + 1 < argc ? argv[i + 1] : ""; return n.size() > 5 && n.compare(n.size() - 5, 5, ".json") == 0; };
      const bool flag = (is_cli_switch(a) && !(a == "--json" && next_json())) ||
                        (a == "--types" && (i + 1 >= argc || std::string(argv[i + 1]).rfind("--", 0) == 0));
      o[a] = flag ? "1" : (i + 1 < argc ? argv[++i] : "");
    } else {
      pos.push_back(a);
    }
  }
  return o;
}

System load(const std::string& path, const std::map<std::string, std::string>& o) {
  auto it = o.find("--topology");
  // --frame N: only that frame is kept (a dump's earlier frames are passed over unread, reading stops after it)
  OpenProgress sel;
  if (auto fr = o.find("--frame"); fr != o.end()) sel.frames.first = sel.frames.last = std::stoul(fr->second);
  Trajectory t;
  try {
    t = open_file(path, it == o.end() ? "" : it->second, sel);
  } catch (const std::exception& e) {
    if (!sel.frames.all() && std::string(e.what()).find("no frames in the selection") != std::string::npos)
      throw std::runtime_error("frame " + std::to_string(sel.frames.first) + " out of range: " + e.what());
    throw;
  }
  if (!sel.frames.all() && t.frames_read == 0 && sel.frames.first >= t.frames())   // a one-frame file: nothing to select
    throw std::runtime_error("frame " + std::to_string(sel.frames.first) + " out of range (file has " + std::to_string(t.frames()) + ")");
  System s = t.frame(0);
  if (!s.unwrapped) make_molecules_whole(s);
  s.notes.insert(s.notes.begin(), sel.frames.all() ? std::to_string(t.frames()) + " frame(s)" : "frame " + std::to_string(sel.frames.first) + " of the file");
  return s;
}


// --complete ml (caps ff apply, relax, md, equilibrate): the bonded terms the force field lacks for the structure, fitted
// by CAPS ML and merged onto it before the assignment goes on (caps/ffcomplete.hpp). Off unless asked, so a result stays
// what the force field alone gives. Prints what was fitted and what is still missing (with the reason).
bool completion_asked(std::map<std::string, std::string>& o) { return o.count("--complete") && parse_complete_option(o["--complete"]); }

// --model M.model, else CAPS's model found as the core finds it (find_ml_model: beside the force field's library, CAPS_HOME,
// the program), else under one of the CAPS roots; throws (resolve_ml_model) when there is none. caps ml, caps ff complete
// and --complete ml all take their model from here.
std::string cli_ml_model(std::map<std::string, std::string>& o, const std::string& ffpath = "") {
  std::string m = o.count("--model") ? o["--model"] : find_ml_model({ffpath});
  if (m.empty())
    for (const auto& r : caps_roots())
      if (!r.empty() && std::filesystem::exists(r + "/data/ml/caps-mace-v1.model")) { m = r + "/data/ml/caps-mace-v1.model"; break; }
  return resolve_ml_model(m, {ffpath});
}

MlFill cli_complete_ml(const System& s, FFDef& def, const std::vector<std::string>& types, std::map<std::string, std::string>& o,
                       const std::string& ffpath, bool quiet) {
  if (!ml_engine_available()) throw std::runtime_error(ml_engine_missing_message());
  CompleteOptions co;
  co.model = cli_ml_model(o, ffpath);
  if (o.count("--step")) co.step = std::stoi(o["--step"]);
  auto log = [&](const std::string& m) { if (!quiet) { std::printf("  CAPS ML: %s\n", m.c_str()); std::fflush(stdout); } };
  MlFill f = fill_missing_with_ml(s, def, types, co, log);
  if (!quiet) {
    if (f.report.missing_before.empty()) std::printf("CAPS ML: nothing missing in %s for this structure\n", def.name.c_str());
    else std::printf("CAPS ML: %d parameter set(s) fitted where %s has none; %zu still missing\n", f.fitted, def.name.c_str(), f.report.missing_after.size());
    for (const auto& t : f.fitted_terms) std::printf("  fitted %s\n", t.c_str());
    for (const auto& t : f.left) std::printf("  left: %s\n", t.c_str());
  }
  return f;
}

// The provenance step of a completion: the model, each fitted term with its source, what was left
ProvStep ml_prov_step(const MlFill& f, const std::string& model) {
  ProvStep p;
  p.engine = "field.complete_ml";
  p.summary = std::to_string(f.fitted) + " parameter set(s) fitted by CAPS ML where the force field has none";
  p.params.push_back({"model", model});
  for (const auto& t : f.fitted_terms) {
    const size_t dot = t.find(" · ");
    p.params.push_back({dot == std::string::npos ? t : t.substr(0, dot), dot == std::string::npos ? "" : t.substr(dot + 4)});
  }
  for (const auto& t : f.left) p.params.push_back({"still missing", t});
  p.approximations.push_back({"Estimated parameters", f.fitted ? std::to_string(f.fitted) + " bonded parameter set(s) from CAPS ML (caps ff complete)" : "none"});
  p.time = now_iso();
  return p;
}

// The force field for commands that take --ff FF.json [--typing RULES] [--charges MODE]: typed by the force field's
// rules (or the atom names in the file); without --ff the built-in GAFF of C and H.
ForceField cli_forcefield(System& s0, std::map<std::string, std::string>& o, bool quiet = false) {
  // a force field written whole for this structure (the Studio's structure.ff.json, caps/ffio.hpp): used as it is
  if (o.count("--ff") && std::filesystem::exists(o["--ff"])) {
    std::ifstream ff(o["--ff"]);
    std::stringstream ss;
    ss << ff.rdbuf();
    if (ss.str().find("caps-assigned-forcefield") != std::string::npos) {
      ForceField F = forcefield_from_json(ss.str());
      if (F.type_index.size() != s0.atoms.size())
        throw std::invalid_argument(o["--ff"] + " is for " + std::to_string(F.type_index.size()) + " atoms; the structure has " + std::to_string(s0.atoms.size()));
      for (size_t i = 0; i < s0.atoms.size(); ++i) s0.atoms[i].charge = F.charge[i];
      s0.has_charges = true;
      if (!quiet) std::printf("force field: %s (read whole from %s)\n", F.name.c_str(), o["--ff"].c_str());
      return F;
    }
  }
  if (!o.count("--ff") && s0.forcefield && s0.forcefield->charge.size() == s0.atoms.size()) {
    if (!quiet) std::printf("force field: %s (the file's own coefficients and charges; give --ff for another)\n", s0.forcefield->name.c_str());
    return *s0.forcefield;
  }
  if (!o.count("--ff")) {
    // coefficients in the file that CAPS cannot rebuild: refused, not replaced by GAFF or UFF
    if (!s0.forcefield_error.empty() && s0.forcefield_atoms == s0.atoms.size()) throw std::runtime_error(s0.forcefield_error);
    const bool ch = std::all_of(s0.atoms.begin(), s0.atoms.end(), [](const Atom& a) { return a.element == 1 || a.element == 6; });
    if (ch) {
      try {
        ForceField F = assign_gaff(s0);
        if (!quiet) std::printf("force field: the built-in GAFF (C and H; the structure has no force field of its own); give --ff for another\n");
        return F;
      } catch (const FieldError&) {   // the built-in GAFF knows sp3 and aromatic carbon only (not C=C of rubbers)
      }
    }
    if (!quiet) std::printf("force field: UFF (the structure has no force field of its own%s); give --ff for another\n", ch ? "; C=C is beyond the built-in GAFF" : "");
    return assign_uff(s0);
  }
  if (is_uff(o["--ff"])) {
    UffOptions uo;
    uo.keep_charges = o.count("--charges") && o["--charges"] == "keep";
    uo.qeq = o.count("--charges") && o["--charges"] == "qeq";
    if (!quiet) std::printf("force field: UFF (every element; %s)\n", uo.qeq ? "QEq charges" : uo.keep_charges ? "charges from the file" : "no charges");
    return assign_uff(s0, uo);
  }
  FFDef def = load_forcefield(ff_path(o["--ff"]));
  if (o.count("--overlay")) merge_overlays(def, o["--overlay"]);   // as caps ff apply: relax, md, equilibrate
  if (o.count("--typing")) load_typing(def, o["--typing"]);
  if (needs_prepare(def)) {   // united atom, shells, ionic bonds: the structure the force field describes
    std::string ch = o.count("--charges") ? o["--charges"] : (file_charges(s0) ? "keep" : "auto");
    const std::string note = prepare_for_forcefield(s0, def, ch);
    if (!note.empty() && !quiet) std::printf("%s\n", note.c_str());
    if (ch == "keep") o["--charges"] = "keep";
  }
  std::vector<std::string> types;
  if (!def.typing.empty()) {
    const TypingResult tr = assign_types(s0, def);
    if (tr.untyped) throw std::runtime_error(untyped_message(def, s0, tr.untyped, &tr.types));
    types = tr.types;
  } else {
    for (const auto& a : s0.atoms) types.push_back(a.name);
  }
  if (completion_asked(o)) (void)cli_complete_ml(s0, def, types, o, ff_path(o["--ff"]), quiet);
  ParamReport rep;
  if (!o.count("--charges") && file_charges(s0) && !quiet)
    std::fprintf(stderr, "charges: the file's own charges are kept, not %s's (--charges auto gives the force field's own)\n", def.name.c_str());
  // --charges auto | increments as caps ff apply and the Studio take them (parameterize_charges), the rest as parameterize
  std::vector<std::string> said;
  ForceField ff = parameterize_charges(s0, def, types, o.count("--charges") ? o["--charges"] : (file_charges(s0) ? "keep" : "types"),
                                       ff_path(o["--ff"]), &rep, false, &said);
  if (!quiet) for (const auto& l : said) std::printf("%s\n", l.c_str());
  if (!rep.missing.empty()) throw std::runtime_error(std::to_string(rep.missing.size()) + " parameters missing in " + def.name + " (caps ff apply lists them)");
  if (!quiet) std::printf("force field: %s\n", def.name.c_str());
  return ff;
}

// The cut-off of relax, md and equilibrate: --cutoff, else the force field's own (its library's "cutoff", a data file's
// record of it), else 10 Å. A coarse-grained model's own cut-off (ForceField::cutoff) applies whatever this says.
void run_cutoff(EnergyOptions& e, const ForceField& ff, std::map<std::string, std::string>& o, bool quiet) {
  if (o.count("--cutoff")) {
    e.cutoff = std::stod(o["--cutoff"]);
    return;
  }
  if (ff.cutoff > 0) return;
  e.cutoff = ff.native_cutoff > 0 ? ff.native_cutoff : 10.0;
  if (!quiet)
    std::printf("cut-off: %g Å (%s; --cutoff for another)\n", e.cutoff, ff.native_cutoff > 0 ? "the force field's own" : "the default, the force field gives none");
}

// one job of caps ml relax: {"elements", "positions", "hold", "fmax", "steps"} -> {"energy_ev", "positions", "held"
// (A / degrees), "converged", "steps"} or {"error"}
Json relax_job(const caps::ml::Potential& pot, const Json& job) {
  Json r = Json::object();
  try {
    std::vector<int> z;
    std::vector<std::array<double, 3>> xyz;
    for (const Json& e : job["elements"].items()) z.push_back(element_from_symbol(e.str()));
    for (const Json& p : job["positions"].items()) xyz.push_back({p[0].number(), p[1].number(), p[2].number()});
    std::vector<caps::ml::Hold> holds;
    if (job.has("hold"))
      for (const Json& h : job["hold"].items()) {
        caps::ml::Hold hd;
        hd.n = static_cast<int>(h.size());
        if (hd.n < 2 || hd.n > 4) throw std::runtime_error("a held coordinate has 2, 3 or 4 atoms");
        for (int k = 0; k < hd.n; ++k) hd.idx[k] = static_cast<int>(h[k].number());
        hd.target = caps::ml::hold_value(hd, xyz);
        holds.push_back(hd);
      }
    const double fmax = job.has("fmax") ? job["fmax"].number() : 0.01;
    const int steps = job.has("steps") ? static_cast<int>(job["steps"].number()) : 400;
    const auto res = caps::ml::relax(pot, z, xyz, holds, fmax, steps);
    r["energy_ev"] = res.energy;
    Json ps = Json::array();
    for (const auto& p : res.pos) {
      Json v = Json::array();
      for (double c : p) v.push_back(Json(c));
      ps.push_back(v);
    }
    r["positions"] = ps;
    Json hv = Json::array();
    for (size_t k = 0; k < holds.size(); ++k)
      hv.push_back(Json(holds[k].n == 2 ? res.held[k] : res.held[k] * 180.0 / M_PI));
    r["held"] = hv;
    r["converged"] = Json(res.converged);
    r["steps"] = Json(static_cast<double>(res.steps));
  } catch (const std::exception& e) {
    r["error"] = std::string(e.what());
  }
  return r;
}

// caps run RECIPE.yaml: every stage on one line as it runs ([k/n] stage · detail · state, a bar while it works), or one
// JSON object per event with --json. Exit 0 done, 2 the recipe or an input is wrong, 3 missing parameters, 4 a run failed.
int cli_recipe(const Json& r, const std::string& file, std::map<std::string, std::string>& o) {
  RecipeOptions ro;
  ro.base_dir = std::filesystem::absolute(file).parent_path().string();
  {
    std::ifstream rf(file);
    ro.sha256 = sha256_hex(std::string((std::istreambuf_iterator<char>(rf)), std::istreambuf_iterator<char>()));
  }
  ro.out_dir = o.count("--out") ? o["--out"] : ".";
  if (o.count("--seed")) ro.seed = std::stoll(o["--seed"]);
  if (o.count("--threads")) ro.threads = std::stoi(o["--threads"]);
  for (const std::string& root : caps_roots())
    if (!root.empty() && std::filesystem::exists(root + "/data/forcefields/catalogue.json")) { ro.forcefield_dir = root + "/data/forcefields"; break; }
  const bool json = o.count("--json") > 0, tty = isatty(fileno(stdout));
  auto esc = [](const std::string& x) { return Json(x).dump(0); };
  // live files (a cluster job): thermo rows, progress lines, frames of the stages that ask for them
  LiveOutput live("run");
  const bool as_job = !default_progress_file().empty();
  const bool resuming = o.count("--resume") > 0;   // the files of the run it continues grow on
  if (o.count("--log") || as_job) live.open_thermo(o.count("--log") ? o["--log"] : "thermo.csv", resuming);
  if (const std::string pf = o.count("--progress-file") ? o["--progress-file"] : default_progress_file(); !pf.empty()) live.open_progress(pf, resuming);
  ro.live = &live;
  ro.frames_path = o.count("--dump") ? o["--dump"] : (std::filesystem::path(ro.out_dir) / "traj.lammpstrj").string();
  // stage checkpoints in a job (or with --resume): a run stopped near its time limit goes on from its last finished stage
  install_stop_signals();
  ro.resume = o.count("--resume") > 0;
  if (as_job || ro.resume || o.count("--checkpoint")) ro.checkpoint_dir = o.count("--checkpoint") ? o["--checkpoint"] : ro.out_dir;
  bool open_line = false;
  auto last_line = std::chrono::steady_clock::now() - std::chrono::hours(1);
  ro.progress = [&](const RecipeEvent& e) {
    if (json) {
      std::printf("{\"stage\":%d,\"stages\":%d,\"name\":%s,\"status\":%s,\"fraction\":%.3f,\"detail\":%s}\n", e.stage, e.stages, esc(e.name).c_str(), esc(e.status).c_str(), e.fraction,
                  esc(e.detail).c_str());
      std::fflush(stdout);
      return;
    }
    char head[64];
    std::snprintf(head, sizeof head, "[%d/%d] %-12s", e.stage, e.stages, e.name.c_str());
    if (e.status == "running") {
      if (!tty) {   // a log file (a batch job): a line every 30 s, so tail -f shows the run moving
        const auto now = std::chrono::steady_clock::now();
        if (now - last_line < std::chrono::seconds(30)) return;
        last_line = now;
        std::printf("%s %3.0f%%  %s\n", head, 100 * e.fraction, e.detail.c_str());
        std::fflush(stdout);
        return;
      }
      const int w = 20, f = int(std::lround(std::clamp(e.fraction, 0.0, 1.0) * w));
      std::string bar = std::string(size_t(f), '#') + std::string(size_t(w - f), '.');
      std::printf("\r%s %s %3.0f%%  %-50.50s", head, bar.c_str(), 100 * e.fraction, e.detail.c_str());
      std::fflush(stdout);
      open_line = true;
      return;
    }
    if (open_line) std::printf("\r%*s\r", 120, "");
    open_line = false;
    size_t cols = 0;   // display width (UTF-8 continuation bytes take no column)
    for (unsigned char c : e.detail) cols += (c & 0xC0) != 0x80;
    std::printf("%s %s%*s %s\n", head, e.detail.c_str(), int(cols < 64 ? 64 - cols : 0), "", e.status.c_str());
    std::fflush(stdout);
  };
  try {
    const auto res = run_recipe(r, ro);
    std::filesystem::remove("resume.txt");
    live.done(true);
    if (json) {
      Json out = Json::object();
      out["status"] = Json("done");
      Json files = Json::array();
      for (const auto& f : res.files) files.push_back(Json(f));
      out["files"] = files;
      out["atoms"] = Json(double(res.system.atoms.size()));
      std::printf("%s\n", out.dump(0).c_str());
    } else {
      for (const auto& p : res.properties)
        if (std::isfinite(p.value)) std::printf("%s: %.6g%s%s\n", p.name.c_str(), p.value, p.unit.empty() ? "" : " ", p.unit.c_str());
      for (const auto& f : res.files) std::printf("wrote %s\n", f.c_str());
    }
    return 0;
  } catch (const RecipeError& e) {
    if (stop_requested() && !ro.checkpoint_dir.empty()) {   // stopped from outside: go on from the last finished stage
      std::string again = "caps run " + file + " --resume --out " + ro.out_dir + (ro.checkpoint_dir != ro.out_dir ? " --checkpoint " + ro.checkpoint_dir : "");
      if (o.count("--seed")) again += " --seed " + o["--seed"];
      if (o.count("--log")) again += " --log " + o["--log"];
      std::ofstream("resume.txt") << again << "\n";
      std::printf("\nstopped (a stop was asked for): finished stages are in %s/recipe.state.json; continue with\n  %s\n", ro.checkpoint_dir.c_str(), again.c_str());
      live.done(false, "interrupted; resume: " + again);
      return 75;
    }
    live.done(false, e.what());
    if (json) std::printf("{\"status\":\"failed\",\"exit\":%d,\"error\":%s}\n", e.code, esc(e.what()).c_str());
    else std::fprintf(stderr, "caps run: %s (exit %d)\n", e.what(), e.code);
    return e.code;
  }
}

// Prints properties; --json OUT and --csv DIR (one CSV per curve).
void cli_report(const std::vector<Property>& props, std::map<std::string, std::string>& o) {
  for (const auto& p : props) {
    std::printf("\n%s", p.name.c_str());
    if (std::isfinite(p.value)) {
      std::printf(": %.6g", p.value);
      if (std::isfinite(p.error)) std::printf(" ± %.2g", p.error);
      if (!p.unit.empty()) std::printf(" %s", p.unit.c_str());
    }
    std::printf("\n  method: %s\n", p.method.c_str());
    for (const auto& [k, v] : p.extra) std::printf("  %s: %.6g\n", k.c_str(), v);
    for (const auto& n : p.notes) std::printf("  note: %s\n", n.c_str());
  }
  if (o.count("--json")) {
    std::ofstream f(o["--json"]);
    f << properties_json(props) << "\n";
    std::printf("\nwrote %s\n", o["--json"].c_str());
  }
  if (o.count("--csv")) {
    std::filesystem::create_directories(o["--csv"]);
    for (const auto& p : props)
      for (size_t k = 0; k < p.series.size(); ++k) {
        const auto path = o["--csv"] + "/" + p.id + (p.series.size() > 1 ? "_" + std::to_string(k + 1) : "") + ".csv";
        std::ofstream f(path);
        f << "# " << p.name << ": " << p.series[k].label << "\n" << p.series[k].x_label << "," << p.series[k].y_label << "\n";
        for (size_t i = 0; i < p.series[k].x.size(); ++i) f << p.series[k].x[i] << "," << p.series[k].y[i] << "\n";
      }
    std::printf("wrote curves to %s\n", o["--csv"].c_str());
  }
}

// Atoms by a selection: all, type:LABEL (a type label or number), or 1-based numbers and ranges (1-10,15).
std::vector<uint32_t> cli_selection(const System& s, const std::string& sel) {
  std::vector<uint32_t> out;
  if (sel == "all") {
    for (uint32_t i = 0; i < s.atoms.size(); ++i) out.push_back(i);
    return out;
  }
  if (sel.rfind("type:", 0) == 0) {
    const auto m = select_type(s, sel.substr(5));
    for (uint32_t i = 0; i < m.size(); ++i)
      if (m[i]) out.push_back(i);
    if (out.empty()) throw std::invalid_argument("no atom has type " + sel.substr(5));
    return out;
  }
  std::stringstream ss(sel);
  std::string part;
  while (std::getline(ss, part, ',')) {
    const size_t dash = part.find('-');
    const long a = std::stol(part.substr(0, dash)), b = dash == std::string::npos ? a : std::stol(part.substr(dash + 1));
    if (a < 1 || b < a || size_t(b) > s.atoms.size()) throw std::invalid_argument("selection " + part + ": atoms are numbered 1 to " + std::to_string(s.atoms.size()));
    for (long k = a; k <= b; ++k) out.push_back(uint32_t(k - 1));
  }
  return out;
}

// A table of numbers: one row per line, columns split by commas, semicolons, tabs or spaces; comment lines (#) and lines
// that are not all numbers (a header) are skipped.
std::vector<std::vector<double>> read_number_table(const std::string& path) {
  std::ifstream f(path);
  if (!f) throw std::runtime_error("cannot open " + path);
  std::vector<std::vector<double>> rows;
  std::string line;
  while (std::getline(f, line)) {
    if (line.empty() || line[0] == '#') continue;
    for (char& c : line)
      if (c == ',' || c == ';' || c == '\t' || c == '\r') c = ' ';
    std::istringstream in(line);
    std::vector<double> r;
    std::string w;
    bool numbers = true;
    while (in >> w) {
      char* end = nullptr;
      const double v = std::strtod(w.c_str(), &end);
      if (end == w.c_str() || *end) { numbers = false; break; }
      r.push_back(v);
    }
    if (numbers && !r.empty()) rows.push_back(std::move(r));
  }
  if (rows.empty()) throw std::runtime_error(path + ": no rows of numbers");
  return rows;
}

int axis_of(const std::string& a) {
  if (a == "x" || a == "0") return 0;
  if (a == "y" || a == "1") return 1;
  if (a == "z" || a == "2") return 2;
  if (a == "xyz" || a == "all" || a == "3") return 3;   // tensile: the three directions averaged
  throw std::invalid_argument("axis must be x, y, z (or xyz for a tensile test averaged over the three)");
}

// --pme [--ewald-rtol 1e-5] [--pme-spacing 1.0] [--pme-order 5]: particle-mesh Ewald instead of damped shifted force
void electrostatics(EnergyOptions& e, std::map<std::string, std::string>& o) {
  if (o.count("--pme")) e.electrostatics = EnergyOptions::Electrostatics::PME;
  if (o.count("--ewald-rtol")) e.ewald_rtol = std::stod(o["--ewald-rtol"]);
  if (o.count("--pme-spacing")) e.pme_spacing = std::stod(o["--pme-spacing"]);
  if (o.count("--pme-order")) e.pme_order = std::stoi(o["--pme-order"]);
  if (o.count("--dsf-alpha")) e.dsf_alpha = std::stod(o["--dsf-alpha"]);
}

// --precision mixed|double (MD, equilibration, Tg scans and tensile pulls): plain pair terms in single precision with SIMD (the default), or all in
// double
void precision(EnergyOptions& e, std::map<std::string, std::string>& o) {
  e.mixed = true;
  if (!o.count("--precision")) return;
  if (o["--precision"] == "double") e.mixed = false;
  else if (o["--precision"] != "mixed") throw std::runtime_error("--precision takes mixed or double");
}

void save_structure(const System& s, const ForceField& ff, const EnergyOptions& e, const std::string& out) {
  auto ends = [&](const char* x) { return out.size() > 4 && out.substr(out.size() - 4) == x; };
  if (ends(".pdb")) write_pdb(s, out);
  else if (ends(".gro")) write_gro(s, out);   // same atoms and order as a GROMACS topology read with --topology
  else if (ends(".xyz")) write_xyz(s, out);
  else if (const std::string why = write_lammps_data_or_structure(s, ff, e, out); !why.empty()) std::fprintf(stderr, "%s: %s\n", out.c_str(), why.c_str());
}

}  // namespace

// POSCAR, CONTCAR (and POSCAR_…) by the file's name, as VASP names them
static bool poscar_name(const std::string& p) {
  std::string n = std::filesystem::path(p).filename().string();
  for (char& c : n) c = char(std::toupper(static_cast<unsigned char>(c)));
  return n.rfind("POSCAR", 0) == 0 || n.rfind("CONTCAR", 0) == 0;
}

// A structure written by its file extension: .pdb .xyz .mol2 .car .gro .sdf/.mol .cif, .vasp or a POSCAR/CONTCAR name,
// else a LAMMPS data file.
static void write_structure_file(const System& s, const std::string& out) {
  auto ends = [&](const char* e) { const std::string x(e); return out.size() > x.size() && out.compare(out.size() - x.size(), x.size(), x) == 0; };
  if (ends(".pdb")) write_pdb(s, out);
  else if (ends(".xyz")) write_xyz(s, out);
  else if (ends(".mol2")) write_mol2(s, out);
  else if (ends(".car")) write_car(s, out);
  else if (ends(".gro")) write_gro(s, out);
  else if (ends(".sdf") || ends(".mol")) write_sdf(s, out);
  else if (ends(".cif")) write_cif(s, out);
  else if (ends(".vasp") || ends(".poscar") || poscar_name(out)) write_poscar(s, out);
  else write_lammps_data(s, out);
}

// The DFT surface & adsorption workbench (caps/dft_commands.hpp): the same commands as the Studio's DFT pages
static int dft_main(const std::string& cmd, int argc, char** argv) {
  std::vector<std::string> pos;
  std::vector<std::pair<std::string, std::string>> flags;
  bool json = false;
  for (int i = 2; i < argc; ++i) {
    const std::string a = argv[i];
    if (a == "--help" || a == "-h") { std::printf("%s", dft_help(cmd).c_str()); return 0; }
    if (a == "--json") { json = true; continue; }
    if (a.rfind("--", 0) == 0 || a == "-o") {
      if (dft_is_switch(cmd, a)) flags.push_back({a, "true"});
      else if (i + 1 < argc) flags.push_back({a, argv[++i]});
      else { std::fprintf(stderr, "caps %s: %s needs a value\n", cmd.c_str(), a.c_str()); return 2; }
    } else pos.push_back(a);
  }
  std::string data;
  for (const std::string& root : caps_roots())
    if (!root.empty() && std::filesystem::exists(root + "/data/sheets/sheets.json")) { data = root + "/data"; break; }
  if (data.empty()) { std::fprintf(stderr, "caps: data/sheets not found (set CAPS_HOME)\n"); return 2; }
  try {
    const Json r = dft_run(cmd, dft_args_from_cli(pos, flags), data);
    if (json) std::printf("%s\n", r.dump(1).c_str());
    else if (r.has("text")) std::printf("%s", r.text("text").c_str());
    else if (r.has("report")) std::printf("%s\nwrote %s\n%s", r.text("report").c_str(), r.text("wrote").c_str(), r.has("vasp_set") ? ("VASP set: " + r["vasp_set"].text("dir") + "\n").c_str() : "");
    else std::printf("%s\n", r.dump(1).c_str());
    return r.has("ok") && !r["ok"].boolean() ? 1 : 0;
  } catch (const std::exception& e) {
    std::fprintf(stderr, "caps %s: %s\n", cmd.c_str(), e.what());
    return 2;
  }
}

// caps job …: CAPS's own engines as batch jobs on a cluster or workstation (caps/cluster_job.hpp). Run on the cluster
// itself (or by the Studio over ssh): the queue commands (sbatch, squeue, sacct, scancel; qsub, qstat, qdel) are called
// here, the rest is the core's.
namespace {

std::string iso_now_cli() {
  const std::time_t t = std::time(nullptr);
  char b[32];
  std::strftime(b, sizeof b, "%Y-%m-%dT%H:%M:%SZ", std::gmtime(&t));
  return b;
}

std::string run_capture(const std::string& cmd) {
#ifdef _WIN32
  FILE* p = _popen(cmd.c_str(), "r");
#else
  FILE* p = popen(cmd.c_str(), "r");
#endif
  if (!p) return "";
  std::string out;
  char buf[4096];
  size_t n;
  while ((n = std::fread(buf, 1, sizeof buf, p)) > 0) out.append(buf, n);
#ifdef _WIN32
  _pclose(p);
#else
  pclose(p);
#endif
  return out;
}

std::string shq(const std::string& s) {   // a single-quoted shell word
  std::string o = "'";
  for (char c : s) o += c == '\'' ? std::string("'\\''") : std::string(1, c);
  return o + "'";
}

std::string first_line(const std::string& s) { return s.substr(0, s.find('\n')); }

int job_usage() {
  std::fprintf(stderr,
               "usage (run on the cluster, or by the Studio over ssh):\n"
               "  caps job profile [--preset slurm|slurm-workspace|pbs|workstation] [--set key=value,…] [--force]   the host profile\n"
               "                   (~/CAPS/host.json: scheduler, account, partition, modules, root, scratch, cpus, mem_per_cpu, time …)\n"
               "  caps job new --title T --kind md [--input FILE,…] [--cpus N] [--mem 2G] [--time 24:00:00] [--partition P]\n"
               "               [--keep-workspace] [--host-profile FILE] [--submit] -- caps md structure.data -o out.data …\n"
               "  caps job new --kind md --array LIST.txt [--input FILE,…] [--submit]   one task per line: TITLE caps …\n"
               "  caps job submit DIR          queue a job folder made by caps job new\n"
               "  caps job status [DIR] [--json]   state, node, scratch, step reached, the latest thermo row, ETA, the log's end\n"
               "  caps job tail [DIR] [-f]     the live log (run.log in the scratch folder while it runs)\n"
               "  caps job poll DIR [--log-offset N --progress-offset N --err-offset N]   everything new since a look, as JSON\n"
               "  caps job collect [DIR]       copy the scratch folder to out/ (a job stopped by its time limit)\n"
               "  caps job resume DIR [--submit] [--time …]   a job stopped near its time limit, again from its checkpoint\n"
               "  caps job list [--json]       every job folder under the root (the Studio's Find my jobs)\n"
               "  caps job cancel DIR          scancel / qdel / kill\n"
               "  caps job scaling --input STRUCTURE [--threads 1,8,32,full] [--steps 2000] [--submit]   threads vs ns/day on one node\n"
               "  caps job collect-bench DIR   the scaling check's table (ns/day, speedup, efficiency, a suggested thread count)\n"
               "  caps job template            the built-in job script ({placeholders})\n");
  return 2;
}

struct JobArgs {
  std::string sub;
  std::vector<std::string> pos;
  std::map<std::string, std::string> opt;
  std::vector<std::string> inputs;
  std::string command;
};

JobArgs parse_job_args(int argc, char** argv) {
  static const std::set<std::string> switches = {"--submit", "--keep-workspace", "--json", "-f", "--force", "--no-keep-workspace"};
  JobArgs a;
  a.sub = argc > 2 ? argv[2] : "";
  for (int i = 3; i < argc; ++i) {
    const std::string s = argv[i];
    if (s == "--") {   // the command: everything after, as one line (each word quoted only when it needs it)
      for (int k = i + 1; k < argc; ++k) {
        std::string w = argv[k];
        const bool plain = !w.empty() && w.find_first_of(" \t'\"$`\\*?;&|<>()") == std::string::npos;
        a.command += (a.command.empty() ? "" : " ") + (plain ? w : shq(w));
      }
      break;
    }
    if (s.rfind("-", 0) == 0 && s.size() > 1) {
      if (switches.count(s)) { a.opt[s] = "1"; continue; }
      if (i + 1 >= argc) throw std::invalid_argument(s + " needs a value");
      const std::string v = argv[++i];
      if (s == "--input") {
        std::stringstream ss(v);
        std::string x;
        while (std::getline(ss, x, ',')) if (!x.empty()) a.inputs.push_back(x);
      } else {
        a.opt[s] = v;
      }
      continue;
    }
    a.pos.push_back(s);
  }
  return a;
}

HostProfile job_profile(const JobArgs& a) {
  const std::string path = a.opt.count("--host-profile") ? a.opt.at("--host-profile") : default_host_profile_path();
  HostProfile h = load_host_profile(path);
  if (a.opt.count("--root")) h.root = a.opt.at("--root");
  if (a.opt.count("--scheduler")) h.scheduler = a.opt.at("--scheduler");
  return h;
}

std::string job_dir_arg(const JobArgs& a) {
  std::string d = a.pos.empty() ? std::string(".") : a.pos[0];
  if (!std::filesystem::exists(std::filesystem::path(d) / "caps-job.json"))
    throw std::invalid_argument(d + " is not a CAPS job folder (no caps-job.json)");
  return std::filesystem::absolute(d).lexically_normal().string();
}

// queue the script; the scheduler's id ("" when it refused)
std::string submit_script(const HostProfile& h, const JobFolder& f, int tasks, std::string* message) {
  const std::string dir = shq(f.dir), script = shq(f.script);
  std::string out;
  if (h.scheduler == "slurm") {
    out = run_capture("cd " + dir + " && sbatch --parsable " + script + " 2>&1");
    const std::string id = first_line(out).substr(0, first_line(out).find(';'));
    if (!id.empty() && std::all_of(id.begin(), id.end(), [](char c) { return std::isdigit(static_cast<unsigned char>(c)); })) return id;
  } else if (h.scheduler == "pbs") {
    out = run_capture("cd " + dir + " && qsub " + script + " 2>&1");
    const std::string id = first_line(out);
    if (!id.empty() && std::isdigit(static_cast<unsigned char>(id[0]))) return id;
  } else {   // a workstation: in the background, tasks one after another
    const std::string run = tasks > 0 ? "for i in $(seq 0 " + std::to_string(tasks - 1) + "); do CAPS_TASK=$i bash " + script + "; done"
                                      : "bash " + script;
    out = run_capture("cd " + dir + " && nohup bash -c " + shq(run) + " > job.out 2>&1 < /dev/null & echo $!");
    const std::string pid = first_line(out);
    if (!pid.empty() && std::isdigit(static_cast<unsigned char>(pid[0]))) return "pid:" + pid;
  }
  if (message) *message = out;
  return "";
}

// before the scheduler sees the script: state submitted (the script's own writes come later and win)
void mark_submitting(const HostProfile& h, const JobFolder& f) {
  Json c = Json::object();
  c["state"] = std::string("submitted");
  c["scheduler"] = h.scheduler;
  c["submitted"] = std::string(iso_now_cli());
  update_job_state(f.dir, c);
  for (const auto& d : f.task_dirs) update_job_state(d, c);
}
// after: the scheduler's id, unless the script already wrote its own (a job that started at once)
void mark_submitted(const HostProfile& h, const JobFolder& f, const std::string& id) {
  auto set_id = [&](const std::string& dir, const std::string& jid, const std::string& task) {
    const Json st = read_job_state(dir);
    if (st.has("slurm_job") && st["slurm_job"].is_string() && !st["slurm_job"].str().empty()) return;
    Json c = Json::object();
    c["slurm_job"] = jid;
    if (!task.empty()) c["array_task"] = task;
    update_job_state(dir, c);
  };
  set_id(f.dir, id, "");
  for (size_t k = 0; k < f.task_dirs.size(); ++k)
    set_id(f.task_dirs[k], h.scheduler == "slurm" ? id + "_" + std::to_string(k) : h.scheduler == "pbs" ? id.substr(0, id.find('[')) + "[" + std::to_string(k) + "]" : id,
           std::to_string(k));
}

// the queue's line for a job ("%T|%r|%M|%l|%N"; empty when it has left the queue) and the accounting after it
std::string queue_line(const std::string& scheduler, const std::string& id) {
  if (id.empty()) return "";
  if (id.rfind("pid:", 0) == 0) {
    const std::string pid = id.substr(4);
    return run_capture("kill -0 " + pid + " 2>/dev/null && echo 'RUNNING|None|||'$(hostname)");
  }
  if (scheduler == "pbs") {
    const std::string s = run_capture("qstat " + shq(id) + " 2>/dev/null | tail -n 1");
    std::stringstream ss(s);
    std::vector<std::string> f;
    std::string x;
    while (ss >> x) f.push_back(x);
    if (f.size() < 5 || f[0].find(id.substr(0, id.find('.'))) == std::string::npos) return "";
    const std::string st = f[4] == "Q" || f[4] == "H" || f[4] == "W" ? "PENDING" : f[4] == "R" ? "RUNNING" : f[4] == "E" ? "COMPLETING" : "";
    return st.empty() ? "" : st + "|None|" + f[3] + "||";
  }
  return first_line(run_capture("squeue -h -j " + shq(id) + " -o '%T|%r|%M|%l|%N' 2>/dev/null"));
}
std::string acct_text(const std::string& scheduler, const std::string& id) {
  if (scheduler != "slurm" || id.empty() || id.rfind("pid:", 0) == 0) return "";
  return run_capture("sacct -j " + shq(id) + " -o State,ExitCode,Elapsed,MaxRSS -P -n 2>/dev/null");
}

Json poll_dir(const HostProfile& h, const std::string& dir, const PollRequest& base) {
  const Json st = read_job_state(dir);
  const std::string id = st.has("slurm_job") && st["slurm_job"].is_string() ? st["slurm_job"].str() : "";
  const std::string sched = st.has("scheduler") && st["scheduler"].is_string() ? st["scheduler"].str() : h.scheduler;
  const std::string q = queue_line(sched, id);
  const std::string acct = q.empty() ? acct_text(sched, id) : "";
  PollRequest p = base;
  p.dir = dir;
  return job_poll(p, q, acct);
}

std::string fmt_duration(double s) {
  if (s < 0) return "unknown";
  char b[64];
  const long t = long(s + 0.5);
  if (t >= 86400) std::snprintf(b, sizeof b, "%ldd %ldh %02ldm", t / 86400, (t % 86400) / 3600, (t % 3600) / 60);
  else std::snprintf(b, sizeof b, "%ldh %02ldm %02lds", t / 3600, (t % 3600) / 60, t % 60);
  return b;
}

int job_main(int argc, char** argv) {
  JobArgs a;
  try {
    a = parse_job_args(argc, argv);
  } catch (const std::exception& e) {
    std::fprintf(stderr, "caps job: %s\n", e.what());
    return 2;
  }
  try {
    if (a.sub.empty() || a.sub == "--help" || a.sub == "-h" || a.sub == "help") return job_usage();
    if (a.sub == "template") {
      const HostProfile h = job_profile(a);
      std::printf("%s", h.job_template.empty() ? job_template_for(h.scheduler).c_str() : h.job_template.c_str());
      return 0;
    }
    if (a.sub == "profile") {
      const std::string path = a.opt.count("--host-profile") ? a.opt.at("--host-profile") : default_host_profile_path();
      HostProfile h = load_host_profile(path);
      bool changed = false;
      if (a.opt.count("--preset")) {
        if (std::filesystem::exists(path) && !a.opt.count("--force")) throw std::invalid_argument(path + " exists: --force to replace it with the preset");
        bool found = false;
        for (const auto& [k, p] : host_presets()) if (k == a.opt.at("--preset")) h = p, found = true;
        if (!found) throw std::invalid_argument("presets: slurm, slurm-workspace, pbs, workstation");
        changed = true;
      }
      if (a.opt.count("--set")) {
        Json j = host_profile_json(h);
        std::stringstream ss(a.opt.at("--set"));
        std::string kv;
        while (std::getline(ss, kv, ',')) {
          const auto eq = kv.find('=');
          if (eq == std::string::npos) throw std::invalid_argument("--set key=value,…");
          const std::string k = kv.substr(0, eq), v = kv.substr(eq + 1);
          if (!j.has(k)) throw std::invalid_argument("no profile key " + k);
          if (j[k].is_number()) j[k] = std::stod(v);
          else if (j[k].kind() == Json::Bool) j[k] = v == "1" || v == "true" || v == "yes";
          else j[k] = v;
        }
        h = host_profile_from_json(j);
        changed = true;
      }
      if (changed) save_host_profile(h, path);
      std::printf("%s\n%s", path.c_str(), host_profile_json(h).dump(2).c_str());
      std::printf("\n");
      return 0;
    }
    if (a.sub == "new") {
      const HostProfile h = job_profile(a);
      JobRequest r;
      r.title = a.opt.count("--title") ? a.opt.at("--title") : "";
      r.kind = a.opt.count("--kind") ? a.opt.at("--kind") : "run";
      r.command = a.command;
      r.inputs = a.inputs;
      if (a.opt.count("--cpus")) r.cpus = std::stoi(a.opt.at("--cpus"));
      if (a.opt.count("--mem")) r.mem = a.opt.at("--mem");
      if (a.opt.count("--time")) r.time = a.opt.at("--time");
      if (a.opt.count("--partition")) r.partition = a.opt.at("--partition");
      if (a.opt.count("--caps")) r.caps_bin = a.opt.at("--caps");
      if (a.opt.count("--keep-workspace")) r.keep_workspace = 1;
      if (a.opt.count("--no-keep-workspace")) r.keep_workspace = 0;
      JobFolder f;
      int tasks = 0;
      if (a.opt.count("--array")) {
        std::ifstream in(a.opt.at("--array"));
        if (!in) throw std::invalid_argument("cannot read " + a.opt.at("--array"));
        std::stringstream ss;
        ss << in.rdbuf();
        const auto list = read_array_list(ss.str());
        tasks = int(list.size());
        f = make_array_job(h, r, list);
      } else {
        if (r.command.empty()) throw std::invalid_argument("give the command after -- (caps job new --title T --kind md -- caps md …)");
        if (r.title.empty()) throw std::invalid_argument("--title (the structure's or study's name)");
        f = make_job(h, r);
      }
      Json out = Json::object();
      out["dir"] = f.dir, out["job"] = f.job, out["script"] = f.script;
      Json td = Json::array();
      for (const auto& d : f.task_dirs) td.push_back(d);
      out["tasks"] = td;
      if (a.opt.count("--submit")) {
        std::string msg;
        mark_submitting(h, f);
        const std::string id = submit_script(h, f, tasks, &msg);
        if (id.empty()) {
          Json back = Json::object();
          back["state"] = std::string("created");
          update_job_state(f.dir, back);
          out["submitted"] = false, out["error"] = msg;
          std::printf("%s\n", a.opt.count("--json") ? out.dump(0).c_str() : ("made " + f.dir + "; the scheduler refused it:\n" + msg).c_str());
          return 1;
        }
        mark_submitted(h, f, id);
        out["submitted"] = true, out["id"] = id;
      }
      if (a.opt.count("--json")) std::printf("%s\n", out.dump(0).c_str());
      else {
        std::printf("%s\n", f.dir.c_str());
        for (const auto& d : f.task_dirs) std::printf("  task %s\n", d.c_str());
        if (out.has("id")) std::printf("submitted: %s\n", out["id"].str().c_str());
        else std::printf("next: caps job submit %s\n", f.dir.c_str());
      }
      return 0;
    }
    if (a.sub == "resume") {   // a job stopped near its time limit: the same title and kind again, from its checkpoint
      const HostProfile h = job_profile(a);
      JobRequest ov;
      if (a.opt.count("--cpus")) ov.cpus = std::stoi(a.opt.at("--cpus"));
      if (a.opt.count("--mem")) ov.mem = a.opt.at("--mem");
      if (a.opt.count("--time")) ov.time = a.opt.at("--time");
      if (a.opt.count("--partition")) ov.partition = a.opt.at("--partition");
      ov.kind = "";
      JobFolder f = make_resume_job(h, job_dir_arg(a), ov);
      Json out = Json::object();
      out["dir"] = f.dir, out["job"] = f.job, out["script"] = f.script;
      if (a.opt.count("--submit")) {
        std::string msg;
        mark_submitting(h, f);
        const std::string id = submit_script(h, f, 0, &msg);
        if (id.empty()) { out["submitted"] = false, out["error"] = msg; }
        else { mark_submitted(h, f, id); out["submitted"] = true, out["id"] = id; }
      }
      if (a.opt.count("--json")) std::printf("%s\n", out.dump(0).c_str());
      else std::printf("%s\n%s\n", f.dir.c_str(), out.has("id") ? ("submitted: " + out["id"].str()).c_str() : ("next: caps job submit " + f.dir).c_str());
      return out.has("submitted") && !out["submitted"].boolean() ? 1 : 0;
    }
    if (a.sub == "scaling") {   // the same short MD at several thread counts, one array task each, on one node shape
      HostProfile h = job_profile(a);
      if (a.inputs.empty()) throw std::invalid_argument("--input STRUCTURE (a cell of the size you will run)");
      std::vector<std::string> tlist;
      if (a.opt.count("--threads")) {
        std::stringstream ss(a.opt.at("--threads"));
        std::string x;
        while (std::getline(ss, x, ',')) if (!x.empty()) tlist.push_back(x);
      } else {
        for (int t : scaling_threads(h.node_cores)) tlist.push_back(std::to_string(t));
      }
      int most = 1;
      for (const auto& t : tlist) if (t != "full") most = std::max(most, std::stoi(t));
      const int cores = h.node_cores > 0 ? h.node_cores : most;
      const std::string steps = a.opt.count("--steps") ? a.opt.at("--steps") : "2000";
      const std::string input = std::filesystem::path(a.inputs[0]).filename().string();
      std::vector<ArrayTask> tasks;
      for (const auto& t : tlist) {
        const std::string th = t == "full" ? std::to_string(cores) : t;
        tasks.push_back({"scaling", "caps md " + input + " -o bench.data --steps " + steps + " --threads " + th + " --thermo 1000 --quiet --bench-json bench.json"});
      }
      JobRequest r;
      r.kind = "bench";
      r.inputs = a.inputs;
      r.cpus = std::max(most, cores);   // every task gets the whole node: the measurement is not crowded
      if (a.opt.count("--time")) r.time = a.opt.at("--time");
      else r.time = "01:00:00";
      JobFolder f = make_array_job(h, r, tasks);
      Json out = Json::object();
      out["dir"] = f.dir, out["job"] = f.job;
      Json td = Json::array();
      for (const auto& d : f.task_dirs) td.push_back(d);
      out["tasks"] = td;
      if (a.opt.count("--submit")) {
        std::string msg;
        mark_submitting(h, f);
        const std::string id = submit_script(h, f, int(tasks.size()), &msg);
        if (id.empty()) out["submitted"] = false, out["error"] = msg;
        else mark_submitted(h, f, id), out["submitted"] = true, out["id"] = id;
      }
      if (a.opt.count("--json")) std::printf("%s\n", out.dump(0).c_str());
      else std::printf("%s\n%s · then: caps job collect-bench %s\n", f.dir.c_str(), out.has("id") ? ("submitted: " + out["id"].str()).c_str() : "not submitted", f.dir.c_str());
      return 0;
    }
    if (a.sub == "collect-bench") {   // the scaling check's table from each task's bench.json
      const std::string dir = a.pos.empty() ? std::string(".") : a.pos[0];
      std::ifstream tl(std::filesystem::path(dir) / "tasks.txt");
      if (!tl) throw std::invalid_argument(dir + " is not a scaling job (no tasks.txt)");
      std::vector<Json> bench;
      int waiting = 0;
      for (std::string line; std::getline(tl, line);) {
        if (line.empty()) continue;
        const Json st = read_job_state(line);
        const auto b = std::filesystem::path(job_live_dir(line, st)) / "bench.json";
        if (!std::filesystem::exists(b)) { ++waiting; continue; }
        std::ifstream bf(b);
        std::stringstream ss;
        ss << bf.rdbuf();
        try { bench.push_back(Json::parse(ss.str())); } catch (const std::exception&) { ++waiting; }
      }
      Json t = scaling_table(bench);
      t["waiting"] = double(waiting);
      std::ofstream(std::filesystem::path(dir) / "scaling.json") << t.dump(1) << "\n";
      if (a.opt.count("--json")) std::printf("%s\n", t.dump(0).c_str());
      else std::printf("%s%s", scaling_text(t).c_str(), waiting ? (std::to_string(waiting) + " task(s) not finished yet\n").c_str() : "");
      return 0;
    }
    if (a.sub == "submit") {
      const HostProfile h = job_profile(a);
      JobFolder f;
      f.dir = job_dir_arg(a);
      f.script = (std::filesystem::path(f.dir) / "job.sh").string();
      int tasks = 0;
      if (std::ifstream tl(std::filesystem::path(f.dir) / "tasks.txt"); tl) {
        std::string line;
        while (std::getline(tl, line)) if (!line.empty()) f.task_dirs.push_back(line), ++tasks;
      }
      std::string msg;
      mark_submitting(h, f);
      const std::string id = submit_script(h, f, tasks, &msg);
      if (id.empty()) { std::fprintf(stderr, "caps job: the scheduler refused it:\n%s", msg.c_str()); return 1; }
      mark_submitted(h, f, id);
      std::printf("submitted: %s\n", id.c_str());
      return 0;
    }
    if (a.sub == "poll") {
      const HostProfile h = job_profile(a);
      PollRequest p;
      if (a.opt.count("--log-offset")) p.log_offset = std::stoll(a.opt.at("--log-offset"));
      if (a.opt.count("--progress-offset")) p.progress_offset = std::stoll(a.opt.at("--progress-offset"));
      if (a.opt.count("--err-offset")) p.err_offset = std::stoll(a.opt.at("--err-offset"));
      std::printf("%s\n", poll_dir(h, job_dir_arg(a), p).dump(0).c_str());
      return 0;
    }
    if (a.sub == "status") {
      const HostProfile h = job_profile(a);
      PollRequest p;
      p.max_bytes = 4000;   // the log's end and the newest progress lines
      const std::string dir = job_dir_arg(a);
      const Json j = poll_dir(h, dir, p);
      if (a.opt.count("--json")) { std::printf("%s\n", j.dump(1).c_str()); return 0; }
      const Json& st = j["job_state"];
      auto sv = [](const Json& o, const char* k) { return o.has(k) && o[k].is_string() ? o[k].str() : std::string(); };
      std::printf("job      %s / %s   (%s)\n", sv(st, "title").c_str(), sv(st, "job").c_str(), dir.c_str());
      std::printf("state    %s", j["state"].str().c_str());
      const Json& q = j["queue"];
      if (q["found"].boolean()) {
        std::printf("  (%s%s%s)", q["state"].str().c_str(), q["reason"].str().empty() ? "" : ", ", q["reason"].str().c_str());
        std::printf("\nnode     %s\nelapsed  %s of %s", q["node"].str().c_str(), q["elapsed"].str().c_str(), q["limit"].str().c_str());
      } else if (j.has("acct")) {
        std::printf("  (sacct %s, exit %s, %s, MaxRSS %s)", j["acct"]["state"].str().c_str(), j["acct"]["exit_code"].str().c_str(), j["acct"]["elapsed"].str().c_str(), j["acct"]["max_rss"].str().c_str());
      }
      std::printf("\nfiles    %s\n", j["live_dir"].str().c_str());
      const auto& pl = j["progress"].items();
      for (auto it = pl.rbegin(); it != pl.rend(); ++it)
        if (it->has("step")) {
          const Json& r = *it;
          std::printf("reached  %s · step %.0f · %.3f ps · T %.1f K · P %.0f atm · ρ %.4f g/cm³ · Epot %.1f\n", sv(r, "stage").c_str(), r["step"].number(),
                      r.has("time_ps") ? r["time_ps"].number() : 0.0, r.has("T") ? r["T"].number() : 0.0, r.has("P") ? r["P"].number() : 0.0,
                      r.has("rho") ? r["rho"].number() : 0.0, r.has("pe") ? r["pe"].number() : 0.0);
          break;
        }
      if (j.has("fraction")) std::printf("done     %.1f %%%s\n", 100 * j["fraction"].number(), j.has("eta_s") ? ("   ETA " + fmt_duration(j["eta_s"].number())).c_str() : "");
      std::string log = j["log"].str();
      std::vector<std::string> lines;
      std::stringstream ls(log);
      std::string line;
      while (std::getline(ls, line)) lines.push_back(line);
      if (!lines.empty()) {
        std::printf("log      (the last lines of %s)\n", "run.log");
        for (size_t k = lines.size() > 6 ? lines.size() - 6 : 0; k < lines.size(); ++k) std::printf("  %s\n", lines[k].c_str());
      }
      return 0;
    }
    if (a.sub == "tail") {
      const std::string dir = job_dir_arg(a);
      const Json st = read_job_state(dir);
      const std::string log = (std::filesystem::path(job_live_dir(dir, st)) / (st.has("log") && st["log"].is_string() ? st["log"].str() : "run.log")).string();
      if (!std::filesystem::exists(log)) { std::fprintf(stderr, "caps job: no log yet (%s)\n", log.c_str()); return 1; }
      std::fprintf(stderr, "%s\n", log.c_str());
      return std::system(("tail -n 40 " + std::string(a.opt.count("-f") ? "-F " : "") + shq(log)).c_str()) == 0 ? 0 : 1;
    }
    if (a.sub == "collect") {
      const std::string dir = job_dir_arg(a);
      const Json st = read_job_state(dir);
      const std::string scr = st.has("scratch") && st["scratch"].is_string() ? st["scratch"].str() : "";
      const auto out = std::filesystem::path(dir) / "out";
      if (scr.empty() || !std::filesystem::is_directory(scr) || std::filesystem::equivalent(scr, out)) {
        std::printf("nothing to collect: the scratch folder is gone or is out/ itself\n");
        return 0;
      }
      std::filesystem::create_directories(out);
      std::filesystem::copy(scr, out, std::filesystem::copy_options::recursive | std::filesystem::copy_options::overwrite_existing);
      std::printf("copied %s → %s\n", scr.c_str(), out.string().c_str());
      return 0;
    }
    if (a.sub == "list") {
      const HostProfile h = job_profile(a);
      const auto root = std::filesystem::path(expand_home(h.root));
      Json all = Json::array();
      std::error_code ec;
      if (std::filesystem::is_directory(root, ec))
        for (const auto& t : std::filesystem::directory_iterator(root, ec)) {
          if (!t.is_directory(ec) || t.path().filename() == "_arrays" || t.path().filename() == "bin" || t.path().filename() == "src") continue;
          for (const auto& jd : std::filesystem::directory_iterator(t.path(), ec)) {
            if (!std::filesystem::exists(jd.path() / "caps-job.json")) continue;
            Json st = read_job_state(jd.path().string());
            st["dir"] = jd.path().string();
            all.push_back(st);
          }
        }
      if (a.opt.count("--json")) { std::printf("%s\n", all.dump(0).c_str()); return 0; }
      for (const auto& st : all.items()) {
        auto sv = [&](const char* k) { return st.has(k) && st[k].is_string() ? st[k].str() : std::string(); };
        std::printf("%-12s %-10s %-24s %s\n", sv("state").c_str(), sv("slurm_job").c_str(), (sv("title") + "/" + sv("job")).c_str(), sv("dir").c_str());
      }
      return 0;
    }
    if (a.sub == "cancel") {
      const HostProfile h = job_profile(a);
      const std::string dir = job_dir_arg(a);
      const Json st = read_job_state(dir);
      const std::string id = st.has("slurm_job") && st["slurm_job"].is_string() ? st["slurm_job"].str() : "";
      if (id.empty()) throw std::invalid_argument("the job was not submitted");
      const std::string sched = st.has("scheduler") && st["scheduler"].is_string() ? st["scheduler"].str() : h.scheduler;
      std::ofstream(std::filesystem::path(dir) / "cancel-requested") << iso_now_cli() << "\n";   // the script records cancelled, not timeout
      const std::string cmd = id.rfind("pid:", 0) == 0 ? "kill " + id.substr(4) : sched == "pbs" ? "qdel " + shq(id) : "scancel " + shq(id);
      const int rc = std::system((cmd + " 2>&1").c_str());
      Json c = Json::object();
      c["state"] = std::string("cancelled");
      update_job_state(dir, c);
      std::printf("%s: %s\n", cmd.c_str(), rc == 0 ? "done" : "the scheduler reported a problem");
      return rc == 0 ? 0 : 1;
    }
    return job_usage();
  } catch (const std::exception& e) {
    std::fprintf(stderr, "caps job: %s\n", e.what());
    return 2;
  }
}

}  // namespace

// The coarse-graining workflow (caps/cg_commands.hpp): the same commands as the Studio's Coarse-grain page
static int cg_main(const std::string& cmd, int argc, char** argv) {
  std::vector<std::string> pos;
  std::vector<std::pair<std::string, std::string>> flags;
  bool json = false;
  for (int i = 2; i < argc; ++i) {
    const std::string a = argv[i];
    if (a == "--help" || a == "-h") { std::printf("%s", cg_help(cmd).c_str()); return 0; }
    if (a == "--json") { json = true; continue; }
    // options: --name, or one dash and a letter (-o, -T, -P); a value starting with - is consumed below
    if (a.rfind("--", 0) == 0 || (a.size() >= 2 && a[0] == '-' && std::isalpha(static_cast<unsigned char>(a[1])))) {
      if (!cg_is_option(cmd, a)) { std::fprintf(stderr, "caps %s: unknown option %s (caps %s --help lists them)\n", cmd.c_str(), a.c_str(), cmd.c_str()); return 2; }
      if (cg_is_switch(cmd, a)) flags.push_back({a, "true"});
      else if (i + 1 < argc) flags.push_back({a, argv[++i]});
      else { std::fprintf(stderr, "caps %s: %s needs a value\n", cmd.c_str(), a.c_str()); return 2; }
    } else pos.push_back(a);
  }
  std::string data;
  for (const std::string& root : caps_roots())
    if (!root.empty() && std::filesystem::exists(root + "/data/cg/mapping_rules.json")) { data = root + "/data"; break; }
  if (data.empty()) { std::fprintf(stderr, "caps: data/cg not found (set CAPS_HOME)\n"); return 2; }
  try {
    const Json r = cg_run(cmd, dft_args_from_cli(pos, flags), data);
    if (json) std::printf("%s\n", r.dump(1).c_str());
    else if (r.has("text")) std::printf("%s", r.text("text").c_str());
    else std::printf("%s\n", r.dump(1).c_str());
    return r.has("ok") && !r["ok"].boolean() ? 1 : 0;
  } catch (const std::exception& e) {
    std::fprintf(stderr, "caps %s: %s\n", cmd.c_str(), e.what());
    return 2;
  }
}

// The command that continues a stopped run (resume.txt): the same options, the checkpoint as the input, the extra ones.
std::string resume_line(const std::string& head, const std::string& input, const std::map<std::string, std::string>& o,
                        const std::map<std::string, std::string>& extra) {
  auto word = [](const std::string& w) {
    if (!w.empty() && w.find_first_of(" \t'\"$`\\*?;&|<>()") == std::string::npos) return w;
    std::string q = "'";
    for (char c : w) q += c == '\'' ? std::string("'\\''") : std::string(1, c);
    return q + "'";
  };
  std::string line = head + " " + word(input);
  for (const auto& [k, v] : o) {
    if (extra.count(k) || k == "--resume" || k == "--resume-at" || k == "--print-protocol") continue;
    line += " " + k;
    if (!is_cli_switch(k)) line += " " + word(v);
  }
  for (const auto& [k, v] : extra) line += " " + k + (v.empty() ? "" : " " + word(v));
  return line;
}

// caps --check-update: the manifest of the latest release (one GET, nothing sent), compared with this build.
// Exit 0: this build is the latest or newer; 1: a newer release is out; 2: cannot tell (offline, turned off, no manifest).
int check_update_main() {
  const auto local = this_build();
  std::printf("this build:     CAPS %s · commit %s · %s\n", local.version.c_str(), local.commit.c_str(), local.date.c_str());
  if (update_checks_disabled()) {
    std::printf("update check:   off (CAPS_NO_UPDATE_CHECK is set)\n");
    return 2;
  }
  ReleaseManifest m;
  try {
    m = parse_release_manifest(fetch_release_manifest(update_feed()));
  } catch (const std::exception& e) {
    std::fprintf(stderr, "caps: could not check for updates: %s\n", e.what());
    return 2;
  }
  std::printf("latest release: CAPS %s · commit %s · %s%s%s\n", m.version.c_str(), m.commit.substr(0, 10).c_str(), m.date.c_str(),
              m.title.empty() ? "" : " · ", m.title.c_str());
  const auto v = decide_update(local, m);
  switch (v) {
    case UpdateVerdict::Update:
      std::printf("update:         a newer CAPS is out · %s\n", download_url(m, platform_key()).c_str());
      for (const auto& h : m.highlights) std::printf("  - %s\n", h.title.c_str());
      return 1;
    case UpdateVerdict::Current: std::printf("update:         none, this is the latest CAPS\n"); return 0;
    case UpdateVerdict::LocalNewer: std::printf("update:         none, this build is newer than the latest release\n"); return 0;
    default: std::printf("update:         cannot tell (this build's date is unknown)\n"); return 2;
  }
}

int main(int argc, char** argv) {
  if (argc == 2 && (std::string(argv[1]) == "--version" || std::string(argv[1]) == "version")) {
    // "caps 0.1.0 (commit a1b2c3d4e5)": the Studio compares it with its own before running jobs on a host
    std::printf("caps %s (commit %s)\n", version_string(), commit_string());
    return 0;
  }
  if (argc == 2 && std::string(argv[1]) == "--check-update") return check_update_main();
  if (argc >= 2 && std::string(argv[1]) == "job") return job_main(argc, argv);
  if (argc >= 2 && is_dft_command(argv[1])) return dft_main(argv[1], argc, argv);
  if (argc >= 2 && is_cg_command(argv[1])) return cg_main(argv[1], argc, argv);
  if (argc < 3) return usage();
  const std::string cmd = argv[1];
  std::vector<std::string> pos;
  std::map<std::string, std::string> o;
  try {
    o = parse(argc, argv, 2, pos);
  } catch (const std::exception& e) {
    std::fprintf(stderr, "caps: %s\n", e.what());
    return 2;
  }
  if (cmd == "grow") {
    try {
      GrowOptions g;
      if (o.count("--chains")) g.chains = std::stoi(o["--chains"]);
      if (o.count("--dp")) g.dp = std::stoi(o["--dp"]);
      if (o.count("--box")) g.box = std::stod(o["--box"]);
      g.density = o.count("--density") ? std::stod(o["--density"]) : 0.5;
      if (o.count("--tacticity")) g.tacticity = tacticity_from_string(o["--tacticity"]);
      if (o.count("--seed")) g.seed = std::stoull(o["--seed"]);
      g.curve = !o.count("--trans");
      if (o.count("--scale")) g.contact_scale = std::stod(o["--scale"]);
      if (o.count("--trials")) g.trials = std::stoi(o["--trials"]);
      if (o.count("--comfortable")) g.comfortable = std::stod(o["--comfortable"]);
      g.escalate = o.count("--escalate");
      g.auto_scale = o.count("--auto-scale") || !o.count("--scale");   // automatic unless --scale fixes it
      if (!o.count("-o")) return usage();
      GrowReport rep;
      System s;
      if (o.count("--units")) {
        // any repeat units (SMILES with two attachment points), comma separated
        auto split = [](const std::string& t) {
          std::vector<std::string> v;
          std::string cur;
          for (char c : t) {
            if (c == ',') { v.push_back(cur); cur.clear(); } else cur += c;
          }
          if (!cur.empty()) v.push_back(cur);
          return v;
        };
        ChainSpec spec;
        for (const auto& u : split(o["--units"])) spec.units.push_back({u, u});
        spec.dp = g.dp;
        spec.tacticity = g.tacticity;
        if (o.count("--sequence")) spec.sequence = sequence_from_string(o["--sequence"]);
        if (o.count("--weights")) for (const auto& w : split(o["--weights"])) spec.weights.push_back(std::stod(w));
        if (o.count("--blocks")) for (const auto& b : split(o["--blocks"])) spec.blocks.push_back(std::stoi(b));
        if (o.count("--pattern")) spec.pattern = o["--pattern"];
        if (o.count("--ff")) spec.forcefield = o["--ff"];
        if (o.count("--method")) {
          const std::string m = o["--method"];
          if (m != "trials" && m != "rosenbluth" && m != "rosenbluth_lj") throw std::invalid_argument("--method trials|rosenbluth|rosenbluth_lj");
          g.method = m == "rosenbluth" ? 1 : m == "rosenbluth_lj" ? 2 : 0;
        }
        if (o.count("--threads")) g.threads = std::stoi(o["--threads"]);
        if (o.count("--lookahead")) g.lookahead = std::clamp(std::stoi(o["--lookahead"]), 1, 4);
        s = grow_chains(spec, g, &rep);
      } else
        s = grow(g, &rep);
      const std::string out = o["-o"];
      auto ends = [&](const char* e) { return out.size() > 4 && out.substr(out.size() - 4) == e; };
      if (ends(".pdb")) write_pdb(s, out);
      else if (ends(".xyz")) write_xyz(s, out);
      else write_lammps_data(s, out);
      for (const auto& n : rep.notes) std::printf("%s\n", n.c_str());
      std::printf("%zu atoms · %zu bonds · %d restarts · %d backtracks · closest contact margin %.2f Å · wrote %s\n", s.atoms.size(),
                  s.bonds.size(), rep.restarts, rep.backtracks, rep.worst_margin, out.c_str());
      return 0;
    } catch (const std::exception& e) {
      std::fprintf(stderr, "caps grow: %s\n", e.what());
      return 1;
    }
  }
  if (cmd == "semicrystalline") {
    try {
      auto split = [](const std::string& t) {
        std::vector<std::string> v;
        std::stringstream ss(t);
        for (std::string x; std::getline(ss, x, ',');) v.push_back(x);
        return v;
      };
      auto write_report = [&](const SemicrystallineReport& r, const std::string& path) {
        std::ofstream f(path);
        if (!f) throw std::runtime_error("cannot write " + path);
        f << semicrystalline_json(r) << "\n";
        std::printf("wrote %s\n", path.c_str());
      };
      if (o.count("--analyse")) {
        if (pos.empty()) return usage();
        System s = open_file(pos[0]).frame(0);
        read_tags(s, pos[0]);
        const SemicrystallineReport r = analyse_semicrystalline(s);
        std::printf("%s", semicrystalline_summary(r).c_str());
        if (o.count("--report")) write_report(r, o["--report"]);
        return 0;
      }
      if (!o.count("-o")) return usage();
      SemicrystallineOptions so;
      if (o.count("--polymer")) so.polymer = o["--polymer"];
      for (const auto& root : caps_roots())   // the polymer crystal library
        if (std::filesystem::exists(root + "/data/crystals/polymers/polymers.json")) { so.data_dir = root + "/data"; break; }
      if (o.count("--repeats")) {
        const auto v = split(o["--repeats"]);
        so.nx = std::stoi(v.at(0)), so.ny = v.size() > 1 ? std::stoi(v[1]) : so.nx;
      }
      if (o.count("--crystal")) so.crystal_nm = std::stod(o["--crystal"]);
      if (o.count("--crystal-cells")) so.nc = std::stoi(o["--crystal-cells"]);
      if (o.count("--amorphous")) so.amorphous_nm = std::stod(o["--amorphous"]);
      if (o.count("--loops")) so.loops = std::stod(o["--loops"]);
      if (o.count("--bridges")) so.bridges = std::stod(o["--bridges"]);
      if (o.count("--tails")) so.tails = std::stod(o["--tails"]);
      if (o.count("--reentry")) so.reentry = o["--reentry"];
      if (o.count("--density")) so.density = std::stod(o["--density"]);
      if (o.count("--seed")) so.seed = std::stoull(o["--seed"]);
      const bool quiet = o.count("--quiet");
      if (o.count("--ff") && !std::filesystem::exists(o["--ff"] + ".json") && !is_uff(o["--ff"])) {
        try { so.united_atom = load_forcefield(ff_path(o["--ff"])).united_atom; } catch (const std::exception&) {}
      }
      so.progress = [&](const std::string& what, double f) {
        if (!quiet) std::fprintf(stderr, "\r%-40s %3.0f %%", what.c_str(), 100 * std::min(1.0, f));
        return true;
      };
      SemicrystallineReport rep;
      System s = build_semicrystalline(so, &rep);
      if (!quiet) std::fprintf(stderr, "\r%-60s\r", "");
      ForceField ff;
      const bool typed = o.count("--ff");
      if (typed) ff = cli_forcefield(s, o, quiet);
      EnergyOptions e;
      e.coulomb = !o.count("--no-coulomb");
      electrostatics(e, o);
      if (typed) run_cutoff(e, ff, o, quiet);
      else ff = default_forcefield(s);
      const auto field = std::make_shared<const ForceField>(ff);
      if (!o.count("--no-relax")) {
        relax_semicrystalline(s, field, e, 1.0, [&](const std::string& what, double f) {
          if (!quiet) std::fprintf(stderr, "\r%-60s %3.0f %%", what.substr(0, 60).c_str(), 100 * std::min(1.0, f));
          return true;
        });
        if (!quiet) std::fprintf(stderr, "\r%-70s\r", "");
        // the structure after relaxation: order parameter, tilt and profiles again (the topology is the builder's)
        const SemicrystallineReport a = analyse_semicrystalline(s);
        rep.crystallinity_p2 = a.crystallinity_p2, rep.tilt_deg = a.tilt_deg, rep.tilt_sd = a.tilt_sd;
        rep.z = a.z, rep.rho = a.rho, rep.p2 = a.p2;
        rep.notes.push_back("relaxed: push-off with the crystal held, then everything minimised (" + ff.name + ")");
      }
      rep.atoms = s.atoms.size();
      const std::string out = o["-o"];
      save_structure(s, ff, e, out);
      write_tags(s, out);
      for (const auto& n : rep.notes) std::printf("%s\n", n.c_str());
      std::printf("%s", semicrystalline_summary(rep).c_str());
      std::printf("wrote %s (groups crystal, amorphous: %s)\n", out.c_str(), tags_sidecar(out).c_str());
      write_report(rep, o.count("--report") ? o["--report"] : std::filesystem::path(out).replace_extension(".semicrystalline.json").string());
      if (o.count("--lammps-input")) {
        const std::string in = o["--lammps-input"];
        const std::string data = std::filesystem::path(out).filename().string();
        LammpsRun run;
        run.kind = LammpsRun::Kind::None;
        write_lammps_input(s, ff, e, data, in, 0, false, run);
        std::ofstream f(in, std::ios::app);
        f << semicrystalline_lammps_run(SemicrystallineEquilibration{});
        std::printf("wrote %s (minimise, 200 ps NVT with the crystal tethered, 2 ns NPT)\n", in.c_str());
      }
      if (o.count("--gromacs")) {
        const std::string stem = o["--gromacs"];
        for (const auto& n : write_gromacs(s, ff, e, stem)) std::printf("gromacs: %s\n", n.c_str());
        std::ofstream ndx(stem + ".ndx");
        for (size_t k = 0; k < s.tags.size(); ++k) {
          ndx << "[ " << tag_group_name(s.tags[k].name) << " ]\n";
          int col = 0;
          for (size_t i = 0; i < s.atoms.size(); ++i)
            if ((s.atoms[i].tags >> k) & 1u) ndx << (i + 1) << (++col % 15 ? " " : "\n");
          ndx << "\n";
        }
        std::ofstream mdp(stem + "_npt.mdp");
        mdp << "; CAPS semicrystalline stack: NPT at 300 K, 1 bar, each box axis on its own (the lamellae keep their shape)\n"
            << "integrator = md\ndt = 0.001\nnsteps = 2000000\nnstxout-compressed = 10000\nnstenergy = 1000\n"
            << "tcoupl = v-rescale\ntc-grps = System\ntau-t = 0.1\nref-t = 300\n"
            << "pcoupl = c-rescale\npcoupltype = anisotropic\ntau-p = 2.0\nref-p = 1 1 1 0 0 0\ncompressibility = 4.5e-5 4.5e-5 4.5e-5 0 0 0\n"
            << "gen-vel = yes\ngen-temp = 300\n"
            << "; non-bonded settings as " << stem << ".mdp (CAPS's matching cut-offs):\n" << gromacs_mdp(s, ff, e);
        std::printf("wrote %s.top, %s.gro, %s.ndx (crystal, amorphous) and %s_npt.mdp\n", stem.c_str(), stem.c_str(), stem.c_str(), stem.c_str());
      }
      return 0;
    } catch (const std::exception& e) {
      std::fprintf(stderr, "caps semicrystalline: %s\n", e.what());
      return 1;
    }
  }
  if (cmd == "solvate") {   // before the generic open: the solute is optional
    try {
      if (o.count("--solvents")) {
        for (const auto& sv : solvent_library()) std::printf("%-12s %-20s %6.3f g/cm³  %s\n", sv.id.c_str(), sv.name.c_str(), sv.density, sv.use.c_str());
        std::printf("salts:");
        for (const auto& sl : salt_library()) std::printf(" %s", sl.id.c_str());
        std::printf("\n");
        return 0;
      }
      if (!o.count("-o")) return usage();
      SolvateOptions so;
      if (o.count("--edge")) so.shape = 0, so.edge = std::stod(o["--edge"]);
      if (o.count("--box")) {
        std::vector<double> v;
        std::stringstream ss(o["--box"]);
        for (std::string x; std::getline(ss, x, ',');) v.push_back(std::stod(x));
        if (v.size() != 3) throw std::invalid_argument("--box a,b,c");
        so.shape = 1, so.edges = {v[0], v[1], v[2]};
      }
      if (o.count("--padding")) so.shape = 2, so.padding = std::stod(o["--padding"]);
      if (o.count("--solvent")) so.solvent = o["--solvent"];
      if (o.count("--model")) so.water_model = o["--model"];
      if (o.count("--density")) so.density = std::stod(o["--density"]);
      if (o.count("--molecules")) so.molecules = std::stoi(o["--molecules"]);
      if (o.count("--salt")) so.salt = o["--salt"];
      if (o.count("--conc")) so.concentration = std::stod(o["--conc"]);
      if (o.count("--neutralise")) so.ion_mode = 1;
      if (o.count("--no-ions")) so.ion_mode = 0;
      if (o.count("--ions")) {
        so.ion_mode = 3;
        const auto c = o["--ions"].find(',');
        if (c == std::string::npos) throw std::invalid_argument("--ions CATIONS,ANIONS");
        so.cations = std::stoi(o["--ions"].substr(0, c)), so.anions = std::stoi(o["--ions"].substr(c + 1));
      }
      if (o.count("--tolerance")) so.tolerance = std::stod(o["--tolerance"]);
      if (o.count("--seed")) so.seed = std::stoull(o["--seed"]);
      System solute;
      const bool has = !pos.empty();
      if (has) solute = load(pos[0], o);
      so.progress = [](const PackProgress& p) {
        std::fprintf(stderr, "\r%-12s loop %d/%d  bad %d  dmin %.2f Å   ", p.stage.c_str(), p.loop, p.loops, p.bad, p.dmin);
        return true;
      };
      SolvateReport rep;
      const System s = solvate(has ? &solute : nullptr, so, &rep);
      std::fprintf(stderr, "\n");
      for (const auto& n : s.notes) std::printf("%s\n", n.c_str());
      std::printf("free volume %.0f of %.0f Å³ · min. distance %.2f Å · %.1f s\n", rep.plan.free_volume, rep.plan.box_volume, rep.pack.dmin, rep.pack.seconds);
      const std::string out = o["-o"];
      auto ends = [&](const char* e) { return out.size() > 4 && out.substr(out.size() - 4) == e; };
      if (ends(".pdb")) write_pdb(s, out);
      else if (ends(".xyz")) write_xyz(s, out);
      else if (ends("mol2")) write_mol2(s, out);
      else if (ends(".car")) write_car(s, out);
      else if (ends(".gro")) write_gro(s, out);
      else write_lammps_data(s, out);
      std::printf("%zu atoms · wrote %s\n", s.atoms.size(), out.c_str());
      return 0;
    } catch (const std::exception& e) {
      std::fprintf(stderr, "caps solvate: %s\n", e.what());
      return 1;
    }
  }
  if (cmd == "provenance") {
    if (pos.empty()) return usage();
    const auto m = read_manifest(pos[0]);
    if (!m) { std::fprintf(stderr, "caps provenance: no %s\n", sidecar_path(pos[0]).c_str()); return 1; }
    if (o.count("--json")) { std::printf("%s\n", manifest_json(*m).dump(2).c_str()); return 0; }
    if (o.count("--bibtex")) { std::printf("%s", bibtex(all_cites(*m)).c_str()); return 0; }
    if (o.count("--methods")) {
      std::vector<std::string> refs;
      std::printf("%s\n\n", methods_text(*m, &refs).c_str());
      for (size_t k = 0; k < refs.size(); ++k) std::printf("[%zu] %s\n", k + 1, refs[k].c_str());
      return 0;
    }
    if (auto c = o.find("--compare"); c != o.end()) {
      const auto other = read_manifest(c->second);
      if (!other) { std::fprintf(stderr, "caps provenance: no %s\n", sidecar_path(c->second).c_str()); return 1; }
      const auto d = compare(*m, *other);
      for (const auto& r : d.rows) std::printf("step %d %-22s %-18s %s  vs  %s\n", r.step + 1, r.engine.c_str(), r.key.c_str(), r.a.c_str(), r.b.c_str());
      for (const auto& n : d.notes) std::printf("%s\n", n.c_str());
      if (d.rows.empty() && d.notes.empty()) std::printf("identical: same inputs, steps, parameters and seeds\n");
      return 0;
    }
    std::printf("%s · %zu steps%s\n", m->generator.c_str(), m->steps.size(), m->deterministic ? " · deterministic" : "");
    for (size_t k = 0; k < m->steps.size(); ++k) {
      const auto& st = m->steps[k];
      std::printf("%2zu  %-22s %s\n", k + 1, st.engine.c_str(), st.summary.c_str());
      for (const auto& [key, v] : st.params) std::printf("      %-18s %s\n", key.c_str(), v.c_str());
      if (!st.rng.empty()) std::printf("      %-18s %s\n", "rng", st.rng.c_str());
      if (!st.cites.empty()) {
        std::string c;
        for (const auto& x : st.cites) c += (c.empty() ? "" : " · ") + x;
        std::printf("      %-18s %s\n", "cites", c.c_str());
      }
    }
    for (const auto& [k, v] : approximations(*m)) std::printf("approximation  %-20s %s\n", k.c_str(), v.c_str());
    for (const auto& [n, h] : m->inputs) std::printf("input          %-20s %s\n", n.c_str(), h.c_str());
    return 0;
  }
  if (cmd == "crystal") {
    try {
      if (o.count("--groups")) {
        const std::string q = pos.empty() ? "" : pos[0];
        const auto* hit = q.empty() ? nullptr : find_space_group(q);
        for (const auto& sg : space_group_settings()) {
          if (hit && sg.number != hit->number) continue;
          std::printf("%-8s %3d  %-18s %-16s %s\n", sg.key.c_str(), sg.number, sg.hm.c_str(), sg.hall.c_str(), crystal_system(sg.number).c_str());
        }
        if (!q.empty() && !hit) throw std::invalid_argument("no space group matches '" + q + "'");
        return 0;
      }
      const double tol = o.count("--tolerance") ? std::stod(o["--tolerance"]) : -1;
      if (o.count("--find-symmetry")) {
        if (pos.empty()) return usage();
        const auto ext = std::filesystem::path(pos[0]).extension().string();
        const System s = ext == ".cif" ? read_cif(pos[0]) : load(pos[0], o);
        const SymmetryFound f = find_symmetry(s, tol > 0 ? tol : 0.1);
        std::printf("%s (No. %d, %s) · %d operations · setting %s\n", f.hm.c_str(), f.number, f.system.c_str(), f.operations, f.key.c_str());
        std::printf("origin shift %.4f %.4f %.4f\n", f.origin[0], f.origin[1], f.origin[2]);
        std::printf("%-8s %-3s %9s %9s %9s\n", "label", "el", "x", "y", "z");
        for (const auto& site : f.sites)
          std::printf("%-8s %-3s %9.5f %9.5f %9.5f\n", site.label.c_str(), element(site.element).symbol, site.frac[0], site.frac[1], site.frac[2]);
        return 0;
      }
      if (!o.count("-o")) return usage();
      CrystalSpec spec;
      if (o.count("--group")) spec.space_group = o["--group"];
      if (o.count("--cell")) {
        std::vector<double> v;
        std::stringstream ss(o["--cell"]);
        for (std::string x; std::getline(ss, x, ',');) v.push_back(std::stod(x));
        if (v.size() != 3 && v.size() != 6) throw std::invalid_argument("--cell a,b,c or a,b,c,alpha,beta,gamma");
        spec.a = v[0], spec.b = v[1], spec.c = v[2];
        if (v.size() == 6) spec.alpha = v[3], spec.beta = v[4], spec.gamma = v[5];
      }
      if (o.count("--sites")) {
        std::stringstream ss(o["--sites"]);
        for (std::string line; std::getline(ss, line, ';');) {
          std::istringstream ls(line);
          CrystalSite site;
          std::string el;
          if (!(ls >> site.label >> el >> site.frac[0] >> site.frac[1] >> site.frac[2])) {
            if (line.find_first_not_of(" \t") == std::string::npos) continue;
            throw std::invalid_argument("a site is 'LABEL EL x y z', not '" + line + "'");
          }
          site.element = element_from_symbol(el);
          if (site.element <= 0) throw std::invalid_argument("unknown element '" + el + "'");
          spec.sites.push_back(site);
        }
      }
      if (o.count("--supercell")) {
        std::stringstream ss(o["--supercell"]);
        size_t k = 0;
        for (std::string x; std::getline(ss, x, ',') && k < 3;) spec.supercell[k++] = std::max(1, std::stoi(x));
      }
      if (tol > 0) spec.tolerance = tol;
      if (spec.sites.empty()) throw std::invalid_argument("no sites: --sites 'LABEL EL x y z; …'");
      if (o.count("--symmetrize")) {
        int moved = 0;
        spec = symmetrize_sites(spec, 0.3, &moved);
        std::printf("%d site(s) moved onto special positions\n", moved);
      }
      CrystalReport rep;
      System s;
      if (o.count("--primitive")) {
        const auto sc = spec.supercell;
        spec.supercell = {1, 1, 1};
        s = build_crystal(spec, &rep);
        const auto* sg = find_space_group(spec.space_group);
        const bool rhomb_axes = sg->key.size() > 2 && sg->key.substr(sg->key.size() - 2) == ":r";
        s = primitive_cell(s, rhomb_axes ? 'P' : sg->hm[0]);
        std::printf("primitive cell: %zu atoms\n", s.atoms.size());
        if (sc[0] * sc[1] * sc[2] > 1) s = supercell(s, sc[0], sc[1], sc[2]);
        s.bonds = crystal_bonds(s);
      } else {
        s = build_crystal(spec, &rep);
      }
      for (const auto& n : s.notes) std::printf("%s\n", n.c_str());
      for (size_t k = 0; k < spec.sites.size(); ++k) std::printf("  %-6s %-2s multiplicity %d\n", spec.sites[k].label.c_str(), element(spec.sites[k].element).symbol, rep.multiplicity[k]);
      for (const auto& n : rep.notes) std::printf("note: %s\n", n.c_str());
      const std::string out = o["-o"];
      write_structure_file(s, out);
      std::printf("%zu atoms · %zu bonds · %.4f g/cm³ · wrote %s\n", s.atoms.size(), s.bonds.size(), s.density(), out.c_str());
      return 0;
    } catch (const std::exception& e) {
      std::fprintf(stderr, "caps crystal: %s\n", e.what());
      return 1;
    }
  }
  if (cmd == "pore") {
    try {
      if (pos.empty() || !o.count("-o")) return usage();
      PoreOptions po;
      const std::string kind = pos[0];
      po.kind = kind == "cylinder" ? PoreKind::Cylinder : kind == "framework" ? PoreKind::Framework : PoreKind::Slit;
      System crystal, fluid;
      if (po.kind != PoreKind::Slit) {
        if (pos.size() < 2) { std::fprintf(stderr, "caps pore: %s needs a crystal (CIF)\n", kind.c_str()); return 1; }
        crystal = read_cif(pos[1]);
        po.crystal = &crystal;
      }
      if (o.count("--width")) po.width = std::stod(o["--width"]);
      if (o.count("--layers")) po.layers = std::stoi(o["--layers"]);
      if (o.count("--lx")) po.lx = std::stod(o["--lx"]);
      if (o.count("--ly")) po.ly = std::stod(o["--ly"]);
      po.vacuum = o.count("--vacuum") > 0;
      if (o.count("--wall")) po.wall = std::stod(o["--wall"]);
      if (o.count("--length")) po.length = std::stod(o["--length"]);
      po.passivate = o.count("--passivate") > 0;
      if (o.count("--supercell")) {
        std::stringstream ss(o["--supercell"]);
        std::string part;
        for (int k = 0; k < 3 && std::getline(ss, part, ','); ++k) po.repeat[k] = std::stoi(part);
      }
      if (o.count("--tolerance")) po.tolerance = std::stod(o["--tolerance"]);
      if (o.count("--seed")) po.seed = std::stoull(o["--seed"]);
      if (o.count("--fluid") && o.count("--count")) {
        BuildOptions bo;
        bo.forcefield = "uff";
        fluid = build_molecule(o["--fluid"], bo).system;
        fluid.title = o["--fluid"];
        po.fluid = &fluid;
        po.count = std::stoi(o["--count"]);
      }
      PoreReport r;
      const System s = build_pore(po, &r);
      const std::string out = o["-o"];
      auto ends = [&](const char* e) { const std::string x = e; return out.size() >= x.size() && out.compare(out.size() - x.size(), x.size(), x) == 0; };
      if (ends(".pdb")) write_pdb(s, out);
      else if (ends(".xyz")) write_xyz(s, out);
      else if (ends("mol2")) write_mol2(s, out);
      else if (ends(".car")) write_car(s, out);
      else write_lammps_data(s, out);
      for (const auto& n : r.notes) std::printf("%s\n", n.c_str());
      std::printf("wrote %s · %zu atoms · walls are molecule 1 (hold them: --fix-mol 1)\n", out.c_str(), s.atoms.size());
      return 0;
    } catch (const std::exception& e) {
      std::fprintf(stderr, "caps pore: %s\n", e.what());
      return 1;
    }
  }
  if (cmd == "nano") {
    // caps nano tube|sheet|particle [CIF] … [--embed --units SMILES --chains N --dp N --density ρ] -o OUT
    try {
      if (pos.empty() || !o.count("-o")) return usage();
      const std::string kind = pos[0];
      NanoReport nr;
      System f;
      std::array<bool, 3> keep{false, false, false};
      if (kind == "tube") {
        NanotubeOptions t;
        if (o.count("--n")) t.n = std::stoi(o["--n"]);
        if (o.count("--m")) t.m = std::stoi(o["--m"]);
        if (o.count("--length")) t.length = std::stod(o["--length"]);
        t.periodic = !o.count("--finite");
        f = nanotube(t, &nr);
        keep = {false, false, t.periodic};
      } else if (kind == "sheet") {
        SheetOptions sh;
        if (o.count("--lx")) sh.lx = std::stod(o["--lx"]);
        if (o.count("--ly")) sh.ly = std::stod(o["--ly"]);
        if (o.count("--layers")) sh.layers = std::stoi(o["--layers"]);
        sh.periodic = !o.count("--flake");
        f = graphene_sheet(sh, &nr);
        keep = {sh.periodic, sh.periodic, false};
      } else if (kind == "particle") {
        if (pos.size() < 2) throw std::invalid_argument("caps nano particle CRYSTAL.cif …");
        ParticleOptions po;
        if (o.count("--shape")) po.shape = particle_shape_from_string(o["--shape"]);
        if (o.count("--radius")) po.radius = std::stod(o["--radius"]);
        po.on_atom = !(o.count("--centre") && o["--centre"] == "cell");
        po.passivate = o.count("--passivate");
        if (o.count("--length")) po.length = std::stod(o["--length"]);
        if (o.count("--height")) po.height = std::stod(o["--height"]);
        if (o.count("--top-ratio")) po.top_ratio = std::stod(o["--top-ratio"]);
        f = nanoparticle(read_cif(pos[1]), po, &nr);
        keep = {false, false, po.shape == ParticleShape::Fibre};
      } else {
        throw std::invalid_argument("kind must be tube, sheet or particle");
      }
      for (const auto& n : nr.notes) std::printf("%s\n", n.c_str());
      System s = f;
      if (o.count("--units")) {   // embed in a polymer matrix
        ChainSpec spec;
        std::stringstream ss(o["--units"]);
        for (std::string u; std::getline(ss, u, ',');) spec.units.push_back({u, u});
        spec.dp = o.count("--dp") ? std::stoi(o["--dp"]) : 20;
        if (o.count("--ff")) spec.forcefield = o["--ff"];
        FillerMatrixOptions fo;
        fo.chains = o.count("--chains") ? std::stoi(o["--chains"]) : 10;
        if (o.count("--density")) fo.density = std::stod(o["--density"]);
        if (o.count("--seed")) fo.grow.seed = std::stoull(o["--seed"]);
        fo.keep_axis = keep;
        FillerReport fr;
        s = embed_filler(f, spec, fo, &fr);
        for (const auto& n : fr.notes) std::printf("%s\n", n.c_str());
      }
      const std::string out = o["-o"];
      auto ends = [&](const char* e) { return out.size() > 4 && out.substr(out.size() - 4) == e; };
      if (ends(".pdb")) write_pdb(s, out);
      else if (ends(".xyz")) write_xyz(s, out);
      else if (ends("mol2")) write_mol2(s, out);
      else if (ends(".car")) write_car(s, out);
      else write_lammps_data(s, out);
      std::printf("%zu atoms · %zu bonds · wrote %s\n", s.atoms.size(), s.bonds.size(), out.c_str());
      return 0;
    } catch (const std::exception& e) {
      std::fprintf(stderr, "caps nano: %s\n", e.what());
      return 1;
    }
  }
  if (cmd == "blend") {
    // caps blend --components SMILES1,SMILES2[,…] [--weights 0.5,0.5] [--chains 8] [--dp 20] [--density 0.5] [--slabs | --droplet] -o OUT
    try {
      if (!o.count("-o") || !o.count("--components")) return usage();
      auto list = [](const std::string& t) {
        std::vector<std::string> v;
        std::stringstream ss(t);
        for (std::string x; std::getline(ss, x, ',');) v.push_back(x);
        return v;
      };
      std::vector<BlendComponent> comps;
      const auto smi = list(o["--components"]);
      const auto w = o.count("--weights") ? list(o["--weights"]) : std::vector<std::string>{};
      for (size_t k = 0; k < smi.size(); ++k) {
        BlendComponent c;
        c.spec.units.push_back({smi[k], smi[k]});
        c.spec.dp = o.count("--dp") ? std::stoi(o["--dp"]) : 20;
        if (o.count("--ff")) c.spec.forcefield = o["--ff"];
        c.weight = k < w.size() ? std::stod(w[k]) : 1.0;
        comps.push_back(c);
      }
      BlendOptions bo;
      if (o.count("--chains")) bo.chains = std::stoi(o["--chains"]);
      if (o.count("--density")) bo.density = std::stod(o["--density"]);
      if (o.count("--seed")) bo.grow.seed = std::stoull(o["--seed"]);
      if (o.count("--slabs")) bo.morphology = BlendMorphology::Slabs;
      if (o.count("--droplet")) bo.morphology = BlendMorphology::Droplet;
      BlendReport br;
      const System s = grow_blend(comps, bo, &br);
      for (const auto& n : br.notes) std::printf("%s\n", n.c_str());
      const std::string out = o["-o"];
      auto ends = [&](const char* e) { return out.size() > 4 && out.substr(out.size() - 4) == e; };
      if (ends(".pdb")) write_pdb(s, out);
      else if (ends(".xyz")) write_xyz(s, out);
      else if (ends("mol2")) write_mol2(s, out);
      else if (ends(".car")) write_car(s, out);
      else write_lammps_data(s, out);
      std::printf("%zu atoms · %zu bonds · wrote %s\n", s.atoms.size(), s.bonds.size(), out.c_str());
      return 0;
    } catch (const std::exception& e) {
      std::fprintf(stderr, "caps blend: %s\n", e.what());
      return 1;
    }
  }
  if (cmd == "interface") {
    // a polymer film grown onto a slab: a slab file, or a crystal (CIF) cleaved with the surface options
    try {
      if (pos.empty() || !o.count("-o") || !o.count("--units")) return usage();
      auto list = [](const std::string& t) {
        std::vector<std::string> v;
        std::stringstream ss(t);
        for (std::string x; std::getline(ss, x, ',');) v.push_back(x);
        return v;
      };
      System slab;
      const std::string src = pos[0];
      if (src.size() > 4 && (src.substr(src.size() - 4) == ".cif" || src.substr(src.size() - 4) == ".CIF")) {
        SlabOptions so;
        if (o.count("--hkl")) {
          const auto v = list(o["--hkl"]);
          so.h = std::stoi(v.at(0)), so.k = std::stoi(v.at(1)), so.l = std::stoi(v.at(2));
        }
        if (o.count("--layers")) so.layers = std::stoi(o["--layers"]);
        if (o.count("--termination")) so.termination = std::stoi(o["--termination"]) - 1;
        if (o.count("--supercell")) {
          const auto v = list(o["--supercell"]);
          so.na = std::stoi(v.at(0)), so.nb = v.size() > 1 ? std::stoi(v[1]) : so.na;
        }
        so.passivate = o.count("--passivate");
        so.vacuum = 10;   // free surfaces (the interface builder sets the final cell)
        SlabReport sr;
        slab = cleave(read_cif(src), so, &sr);
        for (const auto& n : sr.notes) std::printf("%s\n", n.c_str());
      } else {
        slab = open_file(src).frame(0);
      }
      ChainSpec spec;
      for (const auto& u : list(o["--units"])) spec.units.push_back({u, u});
      spec.dp = o.count("--dp") ? std::stoi(o["--dp"]) : 10;
      if (o.count("--tacticity")) spec.tacticity = tacticity_from_string(o["--tacticity"]);
      if (o.count("--sequence")) spec.sequence = sequence_from_string(o["--sequence"]);
      if (o.count("--weights")) for (const auto& w : list(o["--weights"])) spec.weights.push_back(std::stod(w));
      if (o.count("--blocks")) for (const auto& b : list(o["--blocks"])) spec.blocks.push_back(std::stoi(b));
      if (o.count("--pattern")) spec.pattern = o["--pattern"];
      if (o.count("--ff")) spec.forcefield = o["--ff"];
      InterfaceOptions io;
      if (o.count("--film")) io.film = std::stod(o["--film"]);
      if (o.count("--film-density")) io.density = std::stod(o["--film-density"]);
      if (o.count("--chains")) io.chains = std::stoi(o["--chains"]);
      if (o.count("--gap")) io.gap = std::stod(o["--gap"]);
      if (o.count("--vacuum")) io.vacuum = std::stod(o["--vacuum"]);
      if (o.count("--seed")) io.grow.seed = std::stoull(o["--seed"]);
      if (o.count("--scale")) io.grow.contact_scale = std::stod(o["--scale"]);
      GrowReport rep;
      const System s = build_interface(slab, spec, io, &rep);
      for (const auto& n : rep.notes) std::printf("%s\n", n.c_str());
      const std::string out = o["-o"];
      auto ends = [&](const char* e) { return out.size() > 4 && out.substr(out.size() - 4) == e; };
      if (ends(".pdb")) write_pdb(s, out);
      else if (ends(".xyz")) write_xyz(s, out);
      else if (ends("mol2")) write_mol2(s, out);
      else if (ends(".car")) write_car(s, out);
      else write_lammps_data(s, out);
      std::printf("%zu atoms · %zu bonds · wrote %s\n", s.atoms.size(), s.bonds.size(), out.c_str());
      return 0;
    } catch (const std::exception& e) {
      std::fprintf(stderr, "caps interface: %s\n", e.what());
      return 1;
    }
  }
  if (cmd == "surface") {
    // a slab cleaved from a crystal (CIF), with its terminations
    try {
      if (pos.empty()) return usage();
      const System bulk = read_cif(pos[0]);
      for (const auto& n : bulk.notes) std::printf("%s\n", n.c_str());
      SlabOptions so;
      auto ints = [](const std::string& t) {
        std::vector<int> v;
        std::stringstream ss(t);
        for (std::string x; std::getline(ss, x, ',');) v.push_back(std::stoi(x));
        return v;
      };
      if (o.count("--hkl")) {
        const auto v = ints(o["--hkl"]);
        if (v.size() != 3) throw std::invalid_argument("--hkl takes h,k,l");
        so.h = v[0], so.k = v[1], so.l = v[2];
      }
      if (o.count("--layers")) so.layers = std::stoi(o["--layers"]);
      if (o.count("--termination")) so.termination = std::stoi(o["--termination"]) - 1;
      if (o.count("--vacuum")) so.vacuum = std::stod(o["--vacuum"]);
      if (o.count("--supercell")) {
        const auto v = ints(o["--supercell"]);
        so.na = v.at(0), so.nb = v.size() > 1 ? v[1] : v[0];
      }
      if (o.count("--max-strain")) so.max_strain = std::stod(o["--max-strain"]) / 100;
      so.orthogonal = !o.count("--no-orthogonal");
      so.passivate = o.count("--passivate");
      double d = 0;
      const auto terms = slab_terminations(bulk, so.h, so.k, so.l, &d);
      std::printf("(%d%d%d): plane spacing %.4f Å · %zu terminations\n", so.h, so.k, so.l, d, terms.size());
      for (size_t i = 0; i < terms.size(); ++i) std::printf("  %zu  %s\n", i + 1, terms[i].label.c_str());
      if (o.count("--list")) return 0;
      if (!o.count("-o")) return usage();
      SlabReport rep;
      const System slab = cleave(bulk, so, &rep);
      for (const auto& n : rep.notes) std::printf("%s\n", n.c_str());
      const std::string out = o["-o"];
      auto ends = [&](const char* e) { return out.size() > 4 && out.substr(out.size() - 4) == e; };
      if (ends(".pdb")) write_pdb(slab, out);
      else if (ends(".xyz")) write_xyz(slab, out);
      else if (ends("mol2")) write_mol2(slab, out);
      else if (ends(".car")) write_car(slab, out);
      else write_lammps_data(slab, out);
      std::printf("wrote %s\n", out.c_str());
      return 0;
    } catch (const std::exception& e) {
      std::fprintf(stderr, "caps surface: %s\n", e.what());
      return 1;
    }
  }
  if (cmd == "run") {   // a saved pipeline over many inputs: caps run PIPE.yaml --input 'runs/*/X.lammpstrj' [--frame first|last] [--csv OUT] [--out DIR]
    if (pos.empty()) return usage();
    std::ifstream pf(pos[0]);
    if (!pf) throw std::runtime_error("cannot read " + pos[0]);
    const std::string text((std::istreambuf_iterator<char>(pf)), std::istreambuf_iterator<char>());
    if (text.find("caps_pipeline") == std::string::npos) {   // a recipe: caps run RECIPE.yaml [--seed N] [--threads N] [--out DIR] [--json]
      Json r;
      try {
        const auto t0 = text.find_first_not_of(" \t\r\n");
        r = t0 != std::string::npos && text[t0] == '{' ? Json::parse(text) : yaml_parse(text);
      } catch (const std::exception& e) {
        std::fprintf(stderr, "caps run: %s: %s\n", pos[0].c_str(), e.what());
        return 2;
      }
      if (r.is_object() && (r.has("recipe") || r.has("build"))) return cli_recipe(r, pos[0], o);
    }
    std::string name, file, topo;
    Pipeline pl = text.find("caps_pipeline") != std::string::npos ? pipeline_from_yaml(text, &name, &file, &topo) : pipeline_from_json(Json::parse(text));
    if (o.count("--branch")) pl.branch = o["--branch"];   // the branch to run (steps of the others are skipped)
    const std::string pattern = o.count("--input") ? o["--input"] : file;
    const auto inputs = glob_files(pattern);
    if (inputs.empty()) throw std::runtime_error("no files match " + pattern);
    const bool last = !o.count("--frame") || o["--frame"] != "first";
    std::vector<std::string> keys;
    std::vector<std::pair<std::string, std::map<std::string, double>>> rows;
    // the pipeline's outputs: under --out (default outputs/), one folder per input when there are several
    const std::string out_root = o.count("--out") ? o["--out"] : "outputs";
    size_t idx = 0;
    for (const auto& in : inputs) {
      ++idx;
      try {
        std::string tp = o.count("--topology") ? o["--topology"] : "";
        if (tp.empty() && (in.size() > 10 && (in.rfind(".lammpstrj") == in.size() - 10 || in.rfind(".dump") == in.size() - 5))) {
          const auto d = std::filesystem::path(in).replace_extension(".data");
          if (std::filesystem::exists(d)) tp = d.string();
        }
        const Trajectory t = open_file(in, tp);
        const int fr = last ? int(t.frames()) - 1 : 0;
        System f0 = t.frame(size_t(fr));
        if (!f0.unwrapped && f0.cell.valid()) make_molecules_whole(f0);
        const auto st = run_pipeline(f0, pl, fr, t.timesteps.empty() ? 0 : t.timesteps[size_t(fr)], &t);
        std::map<std::string, double> m;
        for (const auto& [k, v] : st.attributes) {
          m[k] = v;
          if (std::find(keys.begin(), keys.end(), k) == keys.end()) keys.push_back(k);
        }
        rows.push_back({in, m});
        if (!pl.outputs.empty()) {
          std::string sub = std::filesystem::path(in).parent_path().filename().string();
          if (sub.empty() || sub == ".") sub = std::filesystem::path(in).stem().string();
          const std::string dir = inputs.size() == 1 ? out_root : out_root + "/" + std::to_string(idx) + "_" + sub;
          for (const auto& line : write_pipeline_outputs(st, pl, dir)) std::fprintf(stderr, "        %s\n", line.c_str());
        }
        std::fprintf(stderr, "done    %s\n", in.c_str());
      } catch (const std::exception& e) {
        rows.push_back({in, {}});
        std::fprintf(stderr, "failed  %s · %s\n", in.c_str(), e.what());
      }
    }
    std::string csv = "input,state";
    for (const auto& k : keys) csv += "," + k;
    csv += "\n";
    char b[40];
    for (const auto& [in, m] : rows) {
      csv += in + (m.empty() ? ",failed" : ",done");
      for (const auto& k : keys) {
        auto it = m.find(k);
        if (it == m.end()) csv += ",";
        else { std::snprintf(b, sizeof b, ",%.10g", it->second); csv += b; }
      }
      csv += "\n";
    }
    if (o.count("--csv")) {
      std::ofstream out(o["--csv"]);
      out << csv;
      std::fprintf(stderr, "wrote %s\n", o["--csv"].c_str());
    } else {
      std::fputs(csv.c_str(), stdout);
    }
    return 0;
  }
  if (cmd == "reproduce") {   // a figure bundle: rebuild its data from its input and pipeline, compare the hashes
    std::vector<std::string> report;
    const bool ok = reproduce_bundle(pos[0], report);
    for (const auto& r : report) std::printf("%s\n", r.c_str());
    std::printf("%s\n", ok ? "reproduced: every data file matches" : "not reproduced");
    return ok ? 0 : 2;
  }
  if (cmd == "bench") {
    try {
      BenchOptions b;
      // the shipped data (caps_roots)
      auto find = [](const std::string& rel, const std::string& probe) {
        for (const auto& r : caps_roots())
          if (std::filesystem::exists(r + "/" + rel + "/" + probe)) return r + "/" + rel;
        return rel;
      };
      b.samples = o.count("--samples") ? o["--samples"] : find("samples", "ps_melt.data");
      b.forcefields = o.count("--forcefields") ? o["--forcefields"] : find("data/forcefields", "catalogue.json");
      if (o.count("--repeats")) b.repeats = std::stoi(o["--repeats"]);
      b.quick = o.count("--quick") > 0;
      std::vector<std::string> ids = pos;
      if (ids.empty() || o.count("--all")) ids = bench_ids();
      b.progress = [](const std::string& id, const std::string& what, double f) {
        std::fprintf(stderr, "\r%-4s %-40s %3.0f %%", id.c_str(), what.c_str(), 100 * f);
        return true;
      };
      std::vector<BenchTable> tables;
      for (const auto& id : ids) {
        tables.push_back(run_bench(id, b));
        const auto& t = tables.back();
        std::fprintf(stderr, "\r%-60s\r", "");
        std::printf("%-4s %-28s %-8s %zu rows  %.1f s\n", t.id.c_str(), t.title.c_str(), t.status.c_str(), t.rows.size(), t.seconds);
      }
      const std::string md = bench_markdown(tables);
      if (o.count("--out")) {
        const std::string dir = o["--out"];
        std::filesystem::create_directories(dir);
        std::ofstream(dir + "/results.md") << md;
        std::ofstream(dir + "/results.tex") << bench_latex(tables);
        for (const auto& t : tables)
          if (!t.rows.empty()) std::ofstream(dir + "/" + t.id + ".csv") << bench_csv(t);
        std::printf("wrote %s/results.md, results.tex and one CSV per table\n", dir.c_str());
      } else
        std::printf("\n%s", md.c_str());
      bool fail = false;
      for (const auto& t : tables) fail |= t.status == "fail";
      return fail ? 1 : 0;
    } catch (const std::exception& e) {
      std::fprintf(stderr, "caps bench: %s\n", e.what());
      return 1;
    }
  }
  if (cmd == "build" && (o.count("--beads") || o.count("--template"))) {
    // coarse-grained: bead SMILES, or a force field's bead template by name
    try {
      if (!o.count("-o")) return usage();
      if (o.count("--template") && !o.count("--ff")) throw std::runtime_error("--template needs --ff (the force field whose templates to use)");
      FFDef def;
      if (o.count("--ff")) def = load_forcefield(ff_path(o["--ff"]));
      if (o.count("--template") && o["--template"] == "list") {
        for (const auto& [k, v] : bead_template_list(def)) std::printf("%-20s %s\n", k.c_str(), v.size() > 90 ? (v.substr(0, 87) + "...").c_str() : v.c_str());
        return 0;
      }
      const std::string text = o.count("--template") ? o["--template"] : o["--beads"];
      if (o.count("--template") && !has_bead_template(def, text)) throw std::runtime_error(def.name + " has no bead template '" + text + "' (--template list)");
      const System s = build_bead_molecule(text, def, o.count("--seed") ? std::stoull(o["--seed"]) : 1);
      const std::string out = o["-o"];
      auto ends = [&](const char* e) { return out.size() > 4 && out.substr(out.size() - std::strlen(e)) == e; };
      if (ends(".pdb")) write_pdb(s, out);
      else if (ends(".xyz")) write_xyz(s, out);
      else write_lammps_data(s, out);
      double q = 0;
      for (const auto& a : s.atoms) q += a.charge;
      std::printf("%zu beads · %zu bonds · charge %+g e%s\n", s.atoms.size(), s.bonds.size(), q, def.name.empty() ? "" : (" · " + def.name).c_str());
      std::printf("wrote %s (relax it with the force field: caps relax %s --ff …)\n", out.c_str(), out.c_str());
      return 0;
    } catch (const std::exception& e) {
      std::fprintf(stderr, "caps build: %s\n", e.what());
      return 1;
    }
  }
  if (cmd == "build") {
    try {
      if (pos.empty() || !o.count("-o")) return usage();
      BuildOptions b;
      if (o.count("--conformers")) b.conformers = std::stoi(o["--conformers"]);
      if (o.count("--seed")) b.seed = std::stoull(o["--seed"]);
      if (o.count("--ff")) b.forcefield = o["--ff"];
      const BuildResult r = build_molecule(pos[0], b);
      const std::string out = o["-o"];
      auto ends = [&](const char* e) { return out.size() > 4 && out.substr(out.size() - std::strlen(e)) == e; };
      auto write = [&](const System& s, const std::string& path) {
        if (ends(".pdb")) write_pdb(s, path);
        else if (ends(".xyz")) write_xyz(s, path);
        else if (ends(".mol2")) write_mol2(s, path);
        else write_lammps_data(s, path);
      };
      if (o.count("--all") && r.conformers.size() > 1) {
        const auto dot = out.find_last_of('.');
        for (size_t k = 0; k < r.conformers.size(); ++k)
          write(molecule_system(r.graph, r.conformers[k].pos), out.substr(0, dot) + "_" + std::to_string(k + 1) + out.substr(dot));
      } else
        write(r.system, out);
      std::printf("%s · %.2f g/mol · %d atoms (%d heavy) · %d rings · %d stereocentres\n", r.info.formula.c_str(), r.info.mass, r.info.atoms,
                  r.info.heavy, r.info.rings, r.info.stereocentres + r.info.stereo_bonds);
      std::printf("%s\n", r.method.c_str());
      const double e0 = r.conformers.front().energy;
      for (size_t k = 0; k < r.conformers.size(); ++k)
        std::printf("  conformer %zu  %s\n", k + 1,
                    r.conformers[k].minimised ? (std::to_string(r.conformers[k].energy - e0).substr(0, 6) + " kcal/mol").c_str() : "not minimised");
      for (const auto& n : r.notes) std::printf("note: %s\n", n.c_str());
      std::printf("wrote %s\n", out.c_str());
      return 0;
    } catch (const std::exception& e) {
      std::fprintf(stderr, "caps build: %s\n", e.what());
      return 1;
    }
  }
  if (cmd == "martini") {
    // Martini 3 small molecules: every molecule of an all-atom structure that is one of the model's molecules becomes its
    // beads with the explicit topology; water and ions are left out, other molecules too (listed)
    try {
      if (pos.empty() || !o.count("-o")) return usage();
      std::string mdir;
      for (const std::string& root : caps_roots())
        if (!root.empty() && std::filesystem::exists(root + "/data/martini/martini3-small-molecules.json")) { mdir = root + "/data/martini"; break; }
      if (mdir.empty()) throw std::runtime_error("data/martini/martini3-small-molecules.json not found (set CAPS_HOME)");
      if (o.count("--martini") && o["--martini"] != "3") throw std::runtime_error("only Martini 3 small molecules are mapped (--martini 3)");
      const System aa = load(pos[0], o);
      const bool geometric = !(o.count("--centre") && o["--centre"] == "mass");   // geometry (default) or mass
      std::vector<std::string> notes;
      const System cg = martini3_all_atom(aa, mdir + "/martini3-small-molecules.json", 1e6, &notes, geometric, false);
      write_lammps_data(cg, o["-o"]);
      if (o.count("--itp")) {
        std::ofstream f(o["--itp"]);
        f << "; Martini 3 molecules written by CAPS\n" << martini3_itp(cg);
      }
      for (const auto& n : notes) std::printf("%s\n", n.c_str());
      std::printf("%zu atoms -> %zu beads; %zu bonds, %zu angles, %zu dihedrals, %zu virtual sites\nwrote %s\n", aa.atoms.size(), cg.atoms.size(),
                  cg.topology->bonds.size(), cg.topology->angles.size(), cg.topology->dihedrals.size(), cg.topology->vsites.size(), o["-o"].c_str());
      return 0;
    } catch (const std::exception& e) {
      std::fprintf(stderr, "caps martini: %s\n", e.what());
      return 1;
    }
  }
  if (cmd == "pack") {
    if (pos.empty()) return usage();
    try {
      PackOptions po;
      std::string output;
      auto items = read_packmol_input(pos[0], po, &output);
      if (o.count("-o")) output = o["-o"];
      if (output.empty()) throw PackError("no output file: give 'output' in the input or -o");
      if (o.count("--threads")) po.threads = std::stoi(o["--threads"]);
      const bool quiet = o.count("--quiet");
      po.progress = [&](const PackProgress& p) {
        if (!quiet && p.stage == "optimisation" && p.bad >= 0 && p.penalty >= 0 && p.loop > 0 && p.dmin == 0 && p.bad > 0)
          std::fprintf(stderr, "  round %d: penalty %.3g, %d molecules in violation\n", p.loop, p.penalty, p.bad);
        return true;
      };
      PackReport rep;
      System s;
      try {
        s = pack(items, po, &rep);
      } catch (const PackError& e) {
        for (const auto& n : rep.notes) std::printf("%s\n", n.c_str());
        std::fprintf(stderr, "caps pack: %s\n", e.what());
        return 4;
      }
      auto ends = [&](const char* e) { const std::string x = e; return output.size() >= x.size() && output.compare(output.size() - x.size(), x.size(), x) == 0; };
      if (ends(".pdb")) write_pdb(s, output);
      else if (ends(".xyz")) write_xyz(s, output);
      else write_lammps_data(s, output);
      for (const auto& n : rep.notes) std::printf("%s\n", n.c_str());
      std::printf("wrote %s\n", output.c_str());
      return 0;
    } catch (const std::exception& e) {
      std::fprintf(stderr, "caps pack: %s\n", e.what());
      return 2;
    }
  }
  if (cmd == "ml") {
    // caps ml eval FILE: CAPS's machine-learned potential on one molecule (non-periodic), natively through libtorch
    try {
      if (pos.size() < 2 || (pos[0] != "eval" && pos[0] != "relax")) return usage();
      const std::string model = cli_ml_model(o);
      if (pos[0] == "relax") {
        // caps ml relax JOBS.json -o OUT.json: a batch of constrained relaxations in one process (the model loaded once).
        // JOBS: {"jobs": [{"elements": ["C", ...], "positions": [[x, y, z], ...], "hold": [[i, j], [i, j, k], [i, j, k, l]],
        // "fmax": 0.01, "steps": 400}]}, atoms 0-based; each held coordinate stays at its value in the given positions.
        // OUT: {"model": ..., "results": [{"energy_ev", "positions", "held" (A / degrees), "converged", "steps"} | {"error"}]}
        if (!o.count("-o")) return usage();
        std::ifstream jf(pos[1]);
        const Json jobs = Json::parse(std::string((std::istreambuf_iterator<char>(jf)), std::istreambuf_iterator<char>()));
        const auto& list = jobs["jobs"].items();
        // the jobs are independent: --threads N workers (default: the machine's cores, at most one per job), each with
        // its own copy of the model and libtorch held to one thread, so the workers do not compete
        int workers = o.count("--threads") ? std::stoi(o["--threads"]) : static_cast<int>(std::thread::hardware_concurrency());
        workers = std::max(1, std::min<int>(workers, static_cast<int>(list.size())));
        caps::ml::set_threads(workers > 1 ? 1 : 0);
        std::vector<std::unique_ptr<caps::ml::Potential>> pots;
        for (int w = 0; w < workers; ++w) pots.push_back(std::make_unique<caps::ml::Potential>(model));
        std::vector<Json> rs(list.size());
        std::atomic<size_t> next{0};
        auto work = [&](int w) {
          for (size_t i = next++; i < list.size(); i = next++) rs[i] = relax_job(*pots[w], list[i]);
        };
        std::vector<std::thread> pool;
        for (int w = 1; w < workers; ++w) pool.emplace_back(work, w);
        work(0);
        for (auto& t : pool) t.join();
        Json results = Json::array();
        for (auto& r : rs) results.push_back(r);
        Json out = Json::object();
        out["model"] = model;
        out["results"] = results;
        std::ofstream(o["-o"]) << out.dump_exact(0) << "\n";
        return 0;
      }
      System s = load(pos[1], o);
      std::vector<int> z;
      std::vector<std::array<double, 3>> xyz;
      for (const auto& a : s.atoms) {
        z.push_back(a.element);
        xyz.push_back(a.pos);
      }
      caps::ml::Potential pot(model);
      if (!pot.covers(z)) {
        std::fprintf(stderr, "caps ml: the model covers only");
        for (int e : pot.elements()) std::fprintf(stderr, " %s", element(e).symbol);
        std::fprintf(stderr, "; this structure has others\n");
        return 3;
      }
      if (o.count("--repeats")) {
        // caps ml eval FILE --repeats N [--batch K] [--threads T]: timing, N calls of K copies of the molecule each
        const int reps = std::stoi(o["--repeats"]), k = o.count("--batch") ? std::stoi(o["--batch"]) : 1;
        if (o.count("--threads")) caps::ml::set_threads(std::stoi(o["--threads"]));
        std::vector<std::vector<int>> zb(k, z);
        std::vector<std::vector<std::array<double, 3>>> pb(k, xyz);
        (void)pot.evaluate_batch(zb, pb);   // warm-up (the first call compiles the graph)
        const auto t0 = std::chrono::steady_clock::now();
        double check = 0;
        for (int r = 0; r < reps; ++r) check += pot.evaluate_batch(zb, pb).at(0).energy;
        const double dt = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
        std::printf("batch %d: %.2f ms per call, %.2f ms per molecule (energy %.10f eV)\n", k, 1e3 * dt / reps, 1e3 * dt / (reps * k), check / reps);
        return 0;
      }
      const auto r = pot.evaluate(z, xyz);
      double fmax = 0, frms = 0;
      for (const auto& f : r.forces) {
        const double m = std::sqrt(f[0] * f[0] + f[1] * f[1] + f[2] * f[2]);
        fmax = std::max(fmax, m);
        frms += m * m;
      }
      frms = std::sqrt(frms / std::max<size_t>(1, r.forces.size()));
      std::printf("model       %s (cutoff %.1f A)\natoms       %zu\nenergy      %.10f eV\nforce       rms %.6f, max %.6f eV/A\n",
                  model.c_str(), pot.cutoff(), z.size(), r.energy, frms, fmax);
      if (o.count("--json")) {
        Json j = Json::object();
        j["model"] = model;
        j["energy_ev"] = r.energy;
        Json fs = Json::array();
        for (const auto& f : r.forces) {
          Json v = Json::array();
          for (double c : f) v.push_back(Json(c));
          fs.push_back(v);
        }
        j["forces_ev_per_a"] = fs;
        std::ofstream(o["--json"]) << j.dump(1) << "\n";
      }
      return 0;
    } catch (const std::exception& e) {
      std::fprintf(stderr, "caps ml: %s\n", e.what());
      return 1;
    }
  }
  if (cmd == "ff") {
    try {
      if (pos.empty()) return usage();
      const std::string sub = pos[0];
      if (sub == "complete") {
        // caps ff complete FILE --ff FF.json -o OVERLAY.json [--report REPORT.json] [--model M.model] [--step 15]: the
        // bonded terms the force field lacks for FILE, found and fitted with CAPS's machine-learned potential
        if (pos.size() < 2 || !o.count("--ff") || !o.count("-o")) return usage();
        if (!ml_engine_available()) throw std::runtime_error(ml_engine_missing_message());
        const std::string model = cli_ml_model(o, ff_path(o["--ff"]));
        System s = load(pos[1], o);
        FFDef ff = load_forcefield(ff_path(o["--ff"]));
        if (o.count("--overlay")) merge_overlays(ff, o["--overlay"]);   // complete the force field as it will run
        CompleteOptions co;
        co.model = model;
        if (o.count("--step")) co.step = std::stoi(o["--step"]);
        const auto rep = complete_forcefield(s, ff, co, [](const std::string& m) { std::printf("%s\n", m.c_str()); std::fflush(stdout); });
        save_forcefield(rep.overlay, o["-o"]);
        int fitted = 0;
        for (const auto& t : rep.terms) fitted += t.status == "fitted";
        std::printf("%s: %zu missing before, %d fitted, %zu missing after; complete %s\nwrote %s\n", ff.name.c_str(),
                    rep.missing_before.size(), fitted, rep.missing_after.size(), rep.complete ? "yes" : "no", o["-o"].c_str());
        for (const auto& t : rep.terms)
          if (t.status != "fitted") {
            std::string c;
            for (size_t k = 0; k < t.classes.size(); ++k) c += (k ? "-" : "") + t.classes[k];
            std::printf("  %s %s %s: %s\n", t.status.c_str(), t.kind.c_str(), c.c_str(), t.why.c_str());
          }
        if (o.count("--report")) std::ofstream(o["--report"]) << complete_report_json(rep, pos[1], ff.name, model) << "\n";
        return rep.complete ? 0 : 3;
      }
      if (sub == "overlays") {
        // caps ff overlays [FF_ID]: the overlay files of the library (parameter refits for --overlay; the Studio's Field
        // page lists the same, caps_forcefield_overlays)
        std::string dir;
        for (const auto& r : caps_roots())
          if (!r.empty() && std::filesystem::exists(r + "/data/forcefields/catalogue.json")) { dir = r + "/data/forcefields"; break; }
        if (dir.empty()) throw std::runtime_error("data/forcefields not found");
        const std::string base = pos.size() > 1 ? pos[1] : "";
        int n = 0;
        for (const auto& ov : library_overlays(dir)) {
          if (!base.empty() && ov.base != base) continue;
          std::printf("%-34s %-12s %s (pairs %d, bonds %d, angles %d, dihedrals %d, impropers %d)\n", ov.file.c_str(),
                      ov.base.empty() ? "-" : ov.base.c_str(), ov.name.c_str(), ov.pairs, ov.bonds, ov.angles, ov.dihedrals, ov.impropers);
          ++n;
        }
        if (!n) std::printf("no overlays%s\n", base.empty() ? "" : (" for " + base).c_str());
        return 0;
      }
      if (sub == "import-lt") {
        if (pos.size() < 2 || !o.count("-o")) return usage();
        FFDef ff = import_moltemplate(pos[1]);
        save_forcefield(ff, o["-o"]);
        std::printf("%s: %zu types, %zu pair, %zu bond, %zu angle, %zu dihedral, %zu improper rules · wrote %s\n", ff.name.c_str(), ff.types.size(),
                    ff.pairs.size(), ff.bonds.size(), ff.angles.size(), ff.dihedrals.size(), ff.impropers.size(), o["-o"].c_str());
        for (const auto& n : ff.notes) std::printf("  %s\n", n.c_str());
        return 0;
      }
      if (sub == "info") {
        if (pos.size() < 2) return usage();
        const FFDef ff = load_forcefield(ff_path(pos[1]));
        std::printf("%s %s\nsource %s\nstyles pair %s · bond %s · angle %s · dihedral %s · improper %s\nmixing %s · special lj %g %g %g · coul %g %g %g · "
                    "cut-off %g · impropers %s\n",
                    ff.name.c_str(), ff.version.c_str(), ff.source.c_str(), ff.pair_style.c_str(), ff.bond_style.c_str(), ff.angle_style.c_str(),
                    ff.dihedral_style.c_str(), ff.improper_style.c_str(), ff.mixing.c_str(), ff.special_lj[0], ff.special_lj[1], ff.special_lj[2],
                    ff.special_coul[0], ff.special_coul[1], ff.special_coul[2], ff.cutoff, ff.improper_order.c_str());
        std::printf("%zu types · rules: %zu pair, %zu bond, %zu angle, %zu dihedral, %zu improper, %zu bond increment\n", ff.types.size(),
                    ff.pairs.size(), ff.bonds.size(), ff.angles.size(), ff.dihedrals.size(), ff.impropers.size(), ff.bond_increments.size());
        for (const auto& r : ff.references) std::printf("reference: %s\n", r.c_str());
        for (const auto& n : ff.notes) std::printf("note: %s\n", n.c_str());
        if (o.count("--types"))
          for (const auto& t : ff.types)
            std::printf("  %-28s %-2s %8.4f %8s  %s\n", t.name.c_str(), t.element ? element(t.element).symbol : "?", t.mass,
                        std::isnan(t.charge) ? "" : std::to_string(t.charge).substr(0, 7).c_str(), t.description.c_str());
        return 0;
      }
      if (sub == "type") {
        // caps ff type FILE --ff FF.json [--typing RULES.json] [-o TYPES.txt] [--explain]
        if (pos.size() < 2 || !o.count("--ff")) return usage();
        System s = load(pos[1], o);
        FFDef ff = load_forcefield(ff_path(o["--ff"]));
        if (o.count("--typing"))
          for (std::stringstream ts(o["--typing"]); ts.good();) {   // several files: comma-separated, later ones on top
            std::string f;
            std::getline(ts, f, ',');
            if (!f.empty()) load_typing(ff, f);
          }
        {
          std::string ch = "types";   // typing only: no charges
          if (const std::string ua = prepare_for_forcefield(s, ff, ch); !ua.empty()) std::printf("note: %s\n", ua.c_str());
        }
        const TypingResult r = assign_types(s, ff);
        for (const auto& n : r.notes) std::printf("note: %s\n", n.c_str());
        const bool explain = o.count("--explain");
        for (size_t i = 0; i < s.atoms.size(); ++i) {
          if (!explain) continue;
          std::string c;
          for (const auto& x : r.candidates[i]) c += (c.empty() ? "" : " ") + x;
          std::printf("%5zu %-3s %-6s %-8s %s   [%s]\n", i + 1, element(s.atoms[i].element).symbol, s.atoms[i].name.c_str(),
                      r.types[i].empty() ? "?" : r.types[i].c_str(), r.why[i].c_str(), c.c_str());
        }
        std::map<std::string, int> count;
        for (const auto& t : r.types) ++count[t.empty() ? "(untyped)" : t];
        std::printf("%zu atoms typed with %s: ", s.atoms.size(), ff.name.c_str());
        for (const auto& [t, k] : count) std::printf("%s×%d ", t.c_str(), k);
        std::printf("\n%d untyped, %d ambiguous\n", r.untyped, r.ambiguous);
        if (o.count("-o")) {
          std::ofstream f(o["-o"]);
          if (!f) throw std::runtime_error("cannot write " + o["-o"]);
          for (const auto& t : r.types) f << (t.empty() ? "?" : t) << "\n";
          std::printf("wrote %s\n", o["-o"].c_str());
        }
        return r.untyped ? 3 : 0;
      }
      if (sub == "apply") {
        if (pos.size() < 2 || !o.count("--ff")) return usage();
        System s = load(pos[1], o);
        const bool uff = is_uff(o["--ff"]);
        ForceField f;
        ParamReport rep;
        double cutoff = 10.0;
        const bool own = o["--ff"] == "file";   // the force field the file carries (an AMBER prmtop)
        MlFill ml_fill;                          // --complete ml: what CAPS ML fitted
        std::string ml_model;
        bool ml_ran = false;
        if ((own || uff) && completion_asked(o)) throw std::runtime_error("--complete ml fills a library force field's gaps (--ff NAME or FF.json), not UFF's or a file's own");
        if (own) {
          if (!s.forcefield || s.forcefield->charge.size() != s.atoms.size())
            throw std::runtime_error(!s.forcefield_error.empty() ? s.forcefield_error
                                                                 : "--ff file: " + pos[1] + " carries no force field (open an AMBER prmtop with its coordinates, or a LAMMPS data file CAPS wrote with one)");
          f = *s.forcefield;
          for (const auto& n : f.notes) std::printf("%s\n", n.c_str());
        } else if (uff) {
          UffOptions uo;
          uo.keep_charges = o.count("--charges") && o["--charges"] == "keep";
          uo.qeq = o.count("--charges") && o["--charges"] == "qeq";
          f = assign_uff(s, uo);
          for (const auto& n : f.notes) std::printf("%s\n", n.c_str());
        } else {
          FFDef ff = load_forcefield(ff_path(o["--ff"]));
          if (o.count("--overlay")) merge_overlays(ff, o["--overlay"]);
          if (o.count("--typing"))
            for (std::stringstream ts(o["--typing"]); ts.good();) {   // several files: comma-separated, later ones on top
              std::string f;
              std::getline(ts, f, ',');
              if (!f.empty()) load_typing(ff, f);
            }
          {   // a united-atom force field: hydrogens on carbon fold into their carbons (charges computed first, then summed)
            std::string ch = o.count("--charges") ? o["--charges"] : "auto";
            if (const std::string ua = prepare_for_forcefield(s, ff, ch); !ua.empty()) {
              std::printf("%s\n", ua.c_str());
              if (ch == "keep") o["--charges"] = "keep";
            }
          }
          std::vector<std::string> types;
          // without --types, atoms are typed by the force field's rules when it has them (else the file's atom names)
          const bool auto_type = !ff.typing.empty() && !o.count("--types") && !o.count("--names");
          if (auto_type) {
            TypingResult r = assign_types(s, ff);
            for (const auto& n : r.notes) std::printf("note: %s\n", n.c_str());
            // a Materials Studio .car carries its force-field types: kept where the force field has them (IFF's inorganic
            // types, which no rule assigns), the rules type the rest
            if (s.source_format == "car") {
              std::set<std::string> ffnames;
              for (const auto& t : ff.types) ffnames.insert(t.name);
              int kept = 0;
              for (size_t i = 0; i < s.atoms.size() && i < r.types.size(); ++i)
                if (ffnames.count(s.atoms[i].name)) {
                  if (r.types[i].empty()) --r.untyped;
                  r.types[i] = s.atoms[i].name;
                  ++kept;
                }
              std::printf("note: %d atoms keep their type from the .car file\n", kept);
            }
            if (r.untyped) {
              std::string l;
              for (size_t i = 0; i < s.atoms.size(); ++i)
                if (r.types[i].empty()) l += " " + std::to_string(i + 1) + element(s.atoms[i].element).symbol;
              throw std::runtime_error(std::to_string(r.untyped) + " atoms match no typing rule:" + l + " (give them with --types)");
            }
            std::printf("typed %zu atoms automatically (%d ambiguous; caps ff type --explain shows why)\n", s.atoms.size(), r.ambiguous);
            types = r.types;
          }
          if (auto_type) {
          } else if (o.count("--types")) {
            // one line per atom, or "index type" lines (1-based) to change only some atoms
            for (const auto& a : s.atoms) types.push_back(a.name);
            std::ifstream tf(o["--types"]);
            if (!tf) throw std::runtime_error("cannot open " + o["--types"]);
            std::string line;
            size_t k = 0;
            while (std::getline(tf, line)) {
              if (auto h = line.find('#'); h != std::string::npos) line = line.substr(0, h);
              std::istringstream ls(line);
              std::vector<std::string> w;
              for (std::string x; ls >> x;) w.push_back(x);
              if (w.empty()) continue;
              if (w.size() >= 2) {
                const size_t i = std::stoul(w[0]);
                if (i < 1 || i > types.size()) throw std::runtime_error("atom " + w[0] + " out of range in " + o["--types"]);
                types[i - 1] = w[1];
              } else {
                if (k >= types.size()) throw std::runtime_error("more types than atoms in " + o["--types"]);
                types[k++] = w[0];
              }
            }
          } else {
            for (const auto& a : s.atoms) types.push_back(a.name);
          }
          if (completion_asked(o)) {
            ml_fill = cli_complete_ml(s, ff, types, o, ff_path(o["--ff"]), false);
            ml_model = ml_fill.model;
            ml_ran = true;
          }
          std::string charges = o.count("--charges") ? o["--charges"] : (file_charges(s) ? "keep" : "auto");
          if (!o.count("--charges") && charges == "keep")
            std::fprintf(stderr, "charges: the file's own charges are kept, not %s's (--charges auto gives the force field's own)\n", ff.name.c_str());
          ParamReport rep;
          // auto: the force field's own charges, else Gasteiger–Marsili, else QEq; the companion's increments when the
          // per-type charges miss the formal charge (the same core function as relax, md, react and the Studio)
          std::vector<std::string> said;
          f = parameterize_charges(s, ff, types, charges, ff_path(o["--ff"]), &rep, o.count("--allow-missing"), &said);
          for (const auto& l : said) std::printf("%s\n", l.c_str());
          cutoff = ff.cutoff;
          for (const auto& n : f.notes) std::printf("%s\n", n.c_str());
        }
        for (const auto& n : rep.notes) std::printf("note: %s\n", n.c_str());
        if (!rep.missing.empty()) {
          std::printf("missing parameters (%zu):\n", rep.missing.size());
          for (const auto& m : rep.missing) std::printf("  %s\n", m.c_str());
        }
        if (o.count("--list")) {
          for (const auto& b : f.bonds) std::printf("bond %u %u  %g %g\n", b.i + 1, b.j + 1, b.k, b.r0);
          for (const auto& a : f.angles) std::printf("angle %u %u %u  %g %g\n", a.i + 1, a.j + 1, a.k + 1, a.kt, a.theta0 * 180 / M_PI);
          for (const auto& d : f.dihedrals) std::printf("dihedral %u %u %u %u  %g %d %g\n", d.i + 1, d.j + 1, d.k + 1, d.l + 1, d.v, d.n, d.delta * 180 / M_PI);
          for (const auto& d : f.impropers) std::printf("improper %u %u %u %u  cvff %g %d %g\n", d.i + 1, d.j + 1, d.k + 1, d.l + 1, d.v, d.n, d.delta * 180 / M_PI);
          for (const auto& d : f.impropers_harmonic) std::printf("improper %u %u %u %u  harmonic %g %g\n", d.i + 1, d.j + 1, d.k + 1, d.l + 1, d.k2, d.chi0 * 180 / M_PI);
          // class II: every coefficient in LAMMPS order, angles in degrees
          const double R = 180 / M_PI;
          for (const auto& b : f.bonds2) std::printf("bond2 %u %u  %.10g %.10g %.10g %.10g\n", b.i + 1, b.j + 1, b.r0, b.k2, b.k3, b.k4);
          for (const auto& a : f.angles2)
            std::printf("angle2 %u %u %u  %.10g %.10g %.10g %.10g  %.10g %.10g %.10g  %.10g %.10g %.10g %.10g\n", a.i + 1, a.j + 1, a.k + 1, a.theta0 * R, a.k2,
                        a.k3, a.k4, a.bb_m, a.bb_r1, a.bb_r2, a.ba_n1, a.ba_n2, a.ba_r1, a.ba_r2);
          for (const auto& d : f.dihedrals2)
            std::printf("dihedral2 %u %u %u %u  %.10g %.10g %.10g %.10g %.10g %.10g  %.10g %.10g %.10g %.10g  %.10g %.10g %.10g %.10g %.10g %.10g %.10g %.10g  "
                        "%.10g %.10g %.10g %.10g %.10g %.10g %.10g %.10g  %.10g %.10g %.10g  %.10g %.10g %.10g\n",
                        d.i + 1, d.j + 1, d.k + 1, d.l + 1, d.k1, d.phi1 * R, d.k2, d.phi2 * R, d.k3, d.phi3 * R, d.mbt[0], d.mbt[1], d.mbt[2], d.mbt_r2,
                        d.ebt_b[0], d.ebt_b[1], d.ebt_b[2], d.ebt_c[0], d.ebt_c[1], d.ebt_c[2], d.ebt_r1, d.ebt_r3, d.at_d[0], d.at_d[1], d.at_d[2],
                        d.at_e[0], d.at_e[1], d.at_e[2], d.at_theta1 * R, d.at_theta2 * R, d.aat_m, d.aat_theta1 * R, d.aat_theta2 * R, d.bb13_n,
                        d.bb13_r1, d.bb13_r3);
          for (const auto& d : f.impropers2)
            std::printf("improper2 %u %u %u %u  %.10g %.10g  %.10g %.10g %.10g %.10g %.10g %.10g\n", d.i + 1, d.j + 1, d.k + 1, d.l + 1, d.kchi,
                        d.chi0 * R, d.m1, d.m2, d.m3, d.theta1 * R, d.theta2 * R, d.theta3 * R);
          for (const auto& d : f.inversions)
            std::printf("inversion %u %u %u %u  %.10g %.10g %s\n", d.c + 1, d.a + 1, d.b + 1, d.d + 1, d.kw, d.w0 * R,
                        d.form == 1 ? "planar" : d.form == 2 ? "fourier" : "harmonic");
          for (const auto& b : f.bonds_x)
            std::printf("bondx %u %u  %s %.10g %.10g %.10g\n", b.i + 1, b.j + 1, b.form == 1 ? "morse" : "gromos", b.a, b.b, b.c);
          for (const auto& a : f.angles_x)
            std::printf("anglex %u %u %u  %s %.10g %.10g\n", a.i + 1, a.j + 1, a.k + 1,
                        a.form == 2 ? "cosine" : a.form == 3 ? "fourier" : a.form > 10 ? ("periodic/" + std::to_string(a.form - 10)).c_str() : "cosine/squared", a.a, a.b * R);
          for (const auto& u : f.urey_bradley) std::printf("ub %u %u  %.10g %.10g\n", u.i + 1, u.k + 1, u.kub, u.r0);
          for (const auto& [ab, pf] : f.pair_func)
            std::printf("pairfunc %s %s  %s %.10g %.10g %.10g\n", f.type_names[ab.first].c_str(), f.type_names[ab.second].c_str(),
                        pf.form == 1 ? "buck" : pf.form == kPairRydberg ? "rydberg" : "morse", pf.a, pf.b, pf.c);
          {
            ForceField f14;
            f14.mixing = f.mixing;
            f14.lj = f.lj14_types.empty() ? f.lj : f.lj14_types;
            for (const auto& p14 : f.pairs14) {
              const PairType pt = f.lj14_types.empty() ? mixed_pair(f, f.type_index[p14[0]], f.type_index[p14[1]])
                                                        : mixed_pair(f14, f.type_index[p14[0]], f.type_index[p14[1]]);
              std::printf("pair14 %u %u  %.10g %.10g  %.10g %.10g\n", p14[0] + 1, p14[1] + 1, pt.eps * f.lj14, pt.sigma, f.lj14, f.coul14);
            }
          }
          for (size_t i = 0; i < f.charge.size(); ++i) std::printf("charge %zu %.10g\n", i + 1, f.charge[i]);
          for (size_t a2 = 0; a2 < f.type_names.size(); ++a2)
            for (size_t b2 = a2; b2 < f.type_names.size(); ++b2) {
              const PairType pt = mixed_pair(f, int(a2), int(b2));
              std::printf("pair %s %s  %.10g %.10g  %s\n", f.type_names[a2].c_str(), f.type_names[b2].c_str(), pt.eps, pt.sigma, f.pair_form.c_str());
            }
        }
        EnergyOptions eo;
        eo.cutoff = o.count("--cutoff") ? std::stod(o["--cutoff"]) : cutoff;
        if (o.count("--no-tail")) eo.tail = false;
        electrostatics(eo, o);
        Evaluator ev(f, eo);
        std::vector<double> x, g;
        for (const auto& a : s.atoms) x.insert(x.end(), a.pos.begin(), a.pos.end());
        const EnergyTerms e = ev.compute(x, s.cell, g);
        std::printf("energy (kcal/mol): bond %.6f  angle %.6f  dihedral %.6f  improper %.6f  vdW %.6f  Coulomb %.6f  total %.6f\n", e.bond, e.angle,
                    e.dihedral, e.improper, e.vdw, e.coulomb, e.total());
        std::printf("virial tensor (kcal/mol): xx %.8g  yy %.8g  zz %.8g  xy %.8g  xz %.8g  yz %.8g\n", e.w[0], e.w[1], e.w[2], e.w[3], e.w[4], e.w[5]);
        if (o.count("--born") && s.cell.valid()) {
          // Born matrix ∂²U/∂η² (kcal/mol, LAMMPS compute born/matrix order); --born lj|bonded|all selects the terms
          // (lj: Lennard-Jones pairs only, charges off; bonded: bonds, angles, torsions, impropers only)
          ForceField fb = f;
          const std::string which = o["--born"];
          if (which == "lj") {
            std::fill(fb.charge.begin(), fb.charge.end(), 0.0);
            fb.bonds.clear(); fb.angles.clear(); fb.dihedrals.clear(); fb.impropers.clear(); fb.impropers_harmonic.clear();
            fb.inversions.clear(); fb.bonds_x.clear(); fb.angles_x.clear(); fb.urey_bradley.clear();
            fb.bonds2.clear(); fb.angles2.clear(); fb.dihedrals2.clear(); fb.impropers2.clear();
          } else if (which == "bonded") {
            for (auto& p : fb.lj) p.eps = 0;
            for (auto& p : fb.lj14_types) p.eps = 0;
            fb.pair_func.clear();
            fb.pair_override.clear();
            std::fill(fb.charge.begin(), fb.charge.end(), 0.0);
          }
          Evaluator eb(fb, eo);
          const Mat6 B = born_matrix(eb, x, s.cell, o.count("--born-strain") ? std::stod(o["--born-strain"]) : 1e-5);
          const double kc = s.cell.volume() / (4184.0 / 6.02214076e23 / 1e-30 / 1e9);   // GPa → kcal/mol
          static const int order[21][2] = {{0, 0}, {1, 1}, {2, 2}, {3, 3}, {4, 4}, {5, 5}, {0, 1}, {0, 2}, {0, 3}, {0, 4}, {0, 5},
                                           {1, 2}, {1, 3}, {1, 4}, {1, 5}, {2, 3}, {2, 4}, {2, 5}, {3, 4}, {3, 5}, {4, 5}};
          std::printf("born matrix (kcal/mol, %s):", which.c_str());
          for (const auto& ij : order) std::printf(" %.10g", B[ij[0]][ij[1]] * kc);
          std::printf("\n");
        }
        if (o.count("--forces")) {
          std::ofstream fo(o["--forces"]);
          char b[128];
          for (size_t i = 0; i < s.atoms.size(); ++i) {
            std::snprintf(b, sizeof b, "%lld %.8f %.8f %.8f\n", static_cast<long long>(s.atoms[i].id), g[3 * i], g[3 * i + 1], g[3 * i + 2]);
            fo << b;
          }
        }
        if (o.count("-o")) {
          for (size_t i = 0; i < s.atoms.size(); ++i) { s.atoms[i].charge = f.charge[i]; s.atoms[i].name = f.atom_type[i]; }
          s.has_charges = true;
          // with an input script beside it, the pair coefficients go in the script and the data file stays plain
          const bool with_in = o.count("--lammps-input") > 0;
          // the force field's own LAMMPS styles (default), or exactly CAPS's energy (--lammps-style exact: the parity
          // benches); --hybrid writes every style as hybrid; --kspace pppm|ewald|dsf|cut; --lammps-cutoff Å
          LammpsStyle ls;
          ls.native = !(o.count("--lammps-style") && o["--lammps-style"] == "exact");
          ls.full_precision = !ls.native;
          ls.hybrid = o.count("--hybrid") > 0;
          if (o.count("--kspace")) ls.coulomb = o["--kspace"];
          if (o.count("--kspace-accuracy")) ls.kspace_accuracy = std::stod(o["--kspace-accuracy"]);
          if (o.count("--lammps-cutoff")) ls.cutoff = std::stod(o["--lammps-cutoff"]);
          if (o.count("--units")) ls.units = o["--units"];
          write_lammps_data_ff(s, f, eo, o["-o"], !with_in, ls);
          std::printf("wrote %s\n", o["-o"].c_str());
          if (with_in) {   // the LAMMPS commands that reproduce this energy with the data file
            namespace fs = std::filesystem;
            const fs::path dp(o["-o"]), ip(o["--lammps-input"]);
            // read_data as the script will see it: the file name when both sit in one folder
            const std::string rel = fs::absolute(dp).parent_path() == fs::absolute(ip).parent_path() ? dp.filename().string() : o["-o"];
            LammpsRun run;   // --lammps-run check (default: a single point with every term) | minimize | nvt | npt | none
            const std::string rk = o.count("--lammps-run") ? o["--lammps-run"] : "check";
            run.kind = rk == "minimize" ? LammpsRun::Kind::Minimize : rk == "nvt" ? LammpsRun::Kind::NVT : rk == "npt" ? LammpsRun::Kind::NPT
                     : rk == "none" ? LammpsRun::Kind::None : LammpsRun::Kind::Check;
            if (o.count("--temp")) run.temperature = std::stod(o["--temp"]);
            if (o.count("--press")) run.pressure = std::stod(o["--press"]);
            if (o.count("--steps")) run.steps = std::stoll(o["--steps"]);
            if (o.count("--dt")) run.dt = std::stod(o["--dt"]);
            std::vector<std::string> lnotes;
            write_lammps_input(s, f, eo, rel, o["--lammps-input"], o.count("--fix-mol") ? std::stoll(o["--fix-mol"]) : 0, true, run, ls, &lnotes);
            for (const auto& n : lnotes) std::printf("lammps: %s\n", n.c_str());
            std::printf("wrote %s (%s)\n", o["--lammps-input"].c_str(), ls.native ? (ls.hybrid ? "the force field's own styles, hybrid" : "the force field's own styles") : "CAPS-exact styles");
            if (o.count("--moltemplate")) {   // the same as a moltemplate system
              std::ofstream lt(o["--moltemplate"]);
              lt << lammps_to_moltemplate(o["-o"], o["--lammps-input"], f.name);
              std::printf("wrote %s (moltemplate)\n", o["--moltemplate"].c_str());
            }
          }
        }
        if (o.count("--gromacs")) {   // STEM.top, STEM.gro and STEM.mdp with the same force field
          for (size_t i = 0; i < s.atoms.size(); ++i) s.atoms[i].charge = f.charge[i];
          for (const auto& n : write_gromacs(s, f, eo, o["--gromacs"])) std::printf("gromacs: %s\n", n.c_str());
          std::printf("wrote %s.top, %s.gro and %s.mdp\n", o["--gromacs"].c_str(), o["--gromacs"].c_str(), o["--gromacs"].c_str());
        }
        if (o.count("--amber")) {   // STEM.prmtop and STEM.inpcrd with the same force field (AMBER, OpenMM, ParmEd)
          for (const auto& n : write_amber(s, f, o["--amber"])) std::printf("amber: %s\n", n.c_str());
          std::printf("wrote %s.prmtop and %s.inpcrd\n", o["--amber"].c_str(), o["--amber"].c_str());
        }
        if (o.count("--dlpoly")) {   // DIR/FIELD, CONFIG and CONTROL
          DlpolyOptions dop;
          dop.title = std::filesystem::path(pos[1]).stem().string();
          if (o.count("--cutoff")) dop.cutoff = std::stod(o["--cutoff"]);
          for (const auto& n : write_dlpoly(s, f, o["--dlpoly"], dop)) std::printf("dlpoly: %s\n", n.c_str());
          std::printf("wrote %s/FIELD, CONFIG and CONTROL\n", o["--dlpoly"].c_str());
        }
        if (ml_ran) {   // the provenance beside the files: the assignment and each term CAPS ML fitted, with its source
          Manifest m;
          ProvStep a;
          a.engine = "field.assign";
          a.summary = f.name + " · " + std::to_string(rep.missing.size()) + " parameters missing after completion";
          a.params = {{"force field", std::filesystem::path(ff_path(o["--ff"])).filename().string()}, {"completion", "ml"}};
          a.time = now_iso();
          m.steps.push_back(a);
          m.steps.push_back(ml_prov_step(ml_fill, ml_model));
          for (const std::string& out : {o.count("-o") ? o["-o"] : std::string(), o.count("--gromacs") ? o["--gromacs"] + ".top" : std::string()})
            if (!out.empty()) {
              write_manifest(m, out);
              std::printf("wrote %s\n", sidecar_path(out).c_str());
            }
        }
        return rep.missing.empty() ? 0 : 3;
      }
      return usage();
    } catch (const std::exception& e) {
      std::fprintf(stderr, "caps ff: %s\n", e.what());
      return 2;
    }
  }
  if (cmd == "elastic") {
    // caps elastic FILE [--topology DATA] [--method strain|fluct] [--ff FF.json ...] [--first --last --stride] [--configs N]
    //              [--strain E] [--ftol F] [--temp T] [--cutoff RC] [--json OUT] [--csv DIR]
    try {
      if (pos.empty()) return usage();
      Trajectory t = open_file(pos[0], o.count("--topology") ? o["--topology"] : "");
      AnalyzeOptions ao;
      if (o.count("--first")) ao.first = std::stol(o["--first"]);
      if (o.count("--last")) ao.last = std::stol(o["--last"]);
      if (o.count("--stride")) ao.stride = std::stol(o["--stride"]);
      auto fr = analysis_frames(t, ao);
      System s0 = t.frame(fr.empty() ? 0 : fr[0]);
      if (!s0.unwrapped) make_molecules_whole(s0);
      auto ff = std::make_shared<ForceField>(cli_forcefield(s0, o));
      EnergyOptions eo;
      if (o.count("--cutoff")) eo.cutoff = std::stod(o["--cutoff"]);
      electrostatics(eo, o);
      if (o.count("--threads")) eo.threads = std::stoi(o["--threads"]);
      const std::string method = o.count("--method") ? o["--method"] : (t.frames() > 1 && o.count("--temp") ? "fluct" : "strain");
      std::vector<Property> props;
      if (method == "fluct-run") {
        if (!o.count("--temp")) throw std::invalid_argument("--temp: the temperature of the NVT run");
        System s = t.frame(fr.empty() ? 0 : fr.back());
        if (!s.unwrapped) make_molecules_whole(s);
        FluctuationRunOptions fo;
        fo.field = ff;
        fo.energy = eo;
        fo.temperature = std::stod(o["--temp"]);
        if (o.count("--ps")) fo.ps = std::stod(o["--ps"]);
        if (o.count("--equilibrate")) fo.equilibrate_ps = std::stod(o["--equilibrate"]);
        if (o.count("--born-every")) fo.born_every = std::stoi(o["--born-every"]);
        if (o.count("--dt")) fo.dt = std::stod(o["--dt"]);
        if (o.count("--seed")) fo.seed = std::stoull(o["--seed"]);
        fo.new_velocities = o.count("--new-velocities") > 0;
        if (!o.count("--quiet"))
          fo.progress = [](const std::string& w, double f) { std::fprintf(stderr, "\r%-40s %3.0f %%", w.c_str(), 100 * f); return true; };
        props = elastic_properties(fluctuation_run(s, fo), "_fluct");
        if (!o.count("--quiet")) std::fprintf(stderr, "\n");
      } else if (method == "fluct") {
        if (!o.count("--temp")) throw std::invalid_argument("--temp: the temperature of the NVT run");
        FluctuationOptions fo;
        fo.ff = ff.get();
        fo.energy = eo;
        fo.temperature = std::stod(o["--temp"]);
        if (o.count("--strain")) fo.strain = std::stod(o["--strain"]);
        if (o.count("--blocks")) fo.blocks = std::stoi(o["--blocks"]);
        props = elastic_properties(fluctuation_elastic(t, fr, fo), "_fluct");
      } else {
        // static: up to --configs frames spread over the selection
        const size_t nc = std::min<size_t>(fr.size(), o.count("--configs") ? std::stoul(o["--configs"]) : 1);
        std::vector<System> cs;
        for (size_t k = 0; k < nc; ++k) {
          System s = t.frame(fr[nc == 1 ? fr.size() - 1 : k * (fr.size() - 1) / (nc - 1)]);
          if (!s.unwrapped) make_molecules_whole(s);
          cs.push_back(std::move(s));
        }
        StaticElasticOptions so;
        so.field = ff;
        so.energy = eo;
        if (o.count("--strain")) so.strain = std::stod(o["--strain"]);
        if (o.count("--ftol")) so.ftol = std::stod(o["--ftol"]);
        if (!o.count("--quiet"))
          so.progress = [](const std::string& w, double f) { std::fprintf(stderr, "\r%-60s %3.0f %%", w.c_str(), 100 * f); return true; };
        props = elastic_properties(static_elastic(cs, so), "");
        if (!o.count("--quiet")) std::fprintf(stderr, "\n");
      }
      cli_report(props, o);
      return 0;
    } catch (const std::exception& e) {
      std::fprintf(stderr, "caps elastic: %s\n", e.what());
      return 1;
    }
  }
  if (cmd == "tensile") {
    // caps tensile DATA -o OUT.data [--axis x|y|z|xyz] [--rate 1e-3 (1/ps)] [--strain 0.2] [--temp 300] [--dt 1] [--fixed-lateral]
    //              [--pressure 1] [--tau-p 1000] [--fit 0.02] [--seed 1] [--new-velocities] [--ff ...] [--dump T.lammpstrj --every N]
    try {
      if (pos.empty() || !o.count("-o")) return usage();
      System s = load(pos[0], o);
      auto ff = std::make_shared<ForceField>(cli_forcefield(s, o));
      TensileOptions to;
      to.field = ff;
      if (o.count("--cutoff")) to.energy.cutoff = std::stod(o["--cutoff"]);
      electrostatics(to.energy, o);
      precision(to.energy, o);
      if (o.count("--threads")) to.energy.threads = std::stoi(o["--threads"]);
      if (o.count("--axis")) to.axis = axis_of(o["--axis"]);
      if (o.count("--rate")) to.rate = std::stod(o["--rate"]);
      if (o.count("--strain")) to.max_strain = std::stod(o["--strain"]);
      if (o.count("--temp")) to.temperature = std::stod(o["--temp"]);
      if (o.count("--dt")) to.dt = std::stod(o["--dt"]);
      if (o.count("--pressure")) to.pressure = std::stod(o["--pressure"]);
      if (o.count("--tau-p")) to.tau_p = std::stod(o["--tau-p"]);
      if (o.count("--fit")) to.fit_strain = std::stod(o["--fit"]);
      if (o.count("--equilibrate")) to.equilibrate_ps = std::stod(o["--equilibrate"]);
      if (o.count("--seed")) to.seed = std::stoull(o["--seed"]);
      to.new_velocities = o.count("--new-velocities") || s.velocities.size() != s.atoms.size();
      to.lateral_pressure = !o.count("--fixed-lateral");
      Trajectory traj;
      if (o.count("--dump")) {
        traj.topology = s;
        to.frame_every = o.count("--every") ? std::stoi(o["--every"]) : 1000;
        to.frame = [&](const std::vector<double>& x, const Cell& c, int64_t step) {
          std::vector<Vec3> p(x.size() / 3);
          for (size_t i = 0; i < p.size(); ++i) p[i] = {x[3 * i], x[3 * i + 1], x[3 * i + 2]};
          traj.positions.push_back(std::move(p));
          traj.cells.push_back(c);
          traj.timesteps.push_back(step);
        };
      }
      const bool quiet = o.count("--quiet");
      size_t shown = 0;
      to.progress = [&](const TensilePoint& p) {
        if (!quiet && shown++ % 20 == 0) std::printf("strain %7.4f  stress %9.2f MPa  lateral %8.4f %8.4f  T %6.1f K\n", p.strain, p.stress, p.lateral1, p.lateral2, p.temperature);
        return true;
      };
      // checkpoints: the strained state with the curve so far and L0, when a stop is asked for (a cluster job near its time
      // limit: exit 75 and resume.txt) and every --checkpoint-every steps
      install_stop_signals();
      TensileResume resume;
      const std::string t_cp = std::filesystem::path(o["-o"]).replace_extension(".restart.data").string();
      const std::string t_curve = std::filesystem::path(o["-o"]).replace_extension(".tensile-curve.json").string();
      if (o.count("--resume-curve")) {
        if (to.axis == 3) throw std::invalid_argument("--resume-curve: a pull averaged over three axes starts again (one axis at a time continues)");
        std::ifstream cf(o["--resume-curve"]);
        std::stringstream ss;
        ss << cf.rdbuf();
        const Json j = Json::parse(ss.str());
        for (int k = 0; k < 3; ++k) resume.L0[k] = j["L0"][size_t(k)].number();
        for (const auto& q : j["curve"].items()) {
          TensilePoint p;
          p.strain = q.num("strain", 0), p.stress = q.num("stress", 0), p.lateral1 = q.num("lateral1", 0), p.lateral2 = q.num("lateral2", 0);
          p.temperature = q.num("T", 0), p.time_ps = q.num("time_ps", 0);
          resume.curve.push_back(p);
        }
        to.resume = &resume;
      }
      if (o.count("--checkpoint-every")) to.checkpoint_every = std::stoll(o["--checkpoint-every"]);
      if (to.axis != 3)
        to.checkpoint = [&](const std::vector<double>& x, const std::vector<double>& v, const Cell& c, const std::vector<TensilePoint>& curve, const double L0[3]) {
          System k = s;
          k.cell = c;
          k.velocities.resize(k.atoms.size());
          for (size_t i = 0; i < k.atoms.size(); ++i) k.atoms[i].pos = {x[3 * i], x[3 * i + 1], x[3 * i + 2]}, k.velocities[i] = {v[3 * i], v[3 * i + 1], v[3 * i + 2]};
          k.unwrapped = true;
          save_structure(k, *ff, to.energy, t_cp);
          Json j = Json::object(), l0 = Json::array(), arr = Json::array();
          for (int a = 0; a < 3; ++a) l0.push_back(L0[a]);
          for (const auto& p : curve) {
            Json q = Json::object();
            q["strain"] = p.strain, q["stress"] = p.stress, q["lateral1"] = p.lateral1, q["lateral2"] = p.lateral2, q["T"] = p.temperature, q["time_ps"] = p.time_ps;
            arr.push_back(q);
          }
          j["L0"] = l0, j["curve"] = arr;
          std::ofstream(t_curve) << j.dump(0) << "\n";
          std::ofstream("resume.txt") << resume_line("caps tensile", t_cp, o, {{"--resume-curve", t_curve}}) << "\n";
        };
      TensileResult r;
      try {
        r = run_tensile(s, to);
      } catch (const DynamicsInterrupted& ex) {
        if (!std::filesystem::exists("resume.txt")) std::ofstream("resume.txt") << resume_line("caps tensile", pos[0], o, {}) << "\n";
        std::ifstream rf("resume.txt");
        std::string again;
        std::getline(rf, again);
        if (o.count("--dump")) write_lammps_dump(traj, o["--dump"]);
        std::printf("stopped at step %lld (a stop was asked for); continue with\n  %s\n", static_cast<long long>(ex.step), again.c_str());
        return 75;
      }
      std::filesystem::remove("resume.txt");   // finished: nothing to continue
      save_structure(s, *ff, to.energy, o["-o"]);
      if (o.count("--dump")) write_lammps_dump(traj, o["--dump"]);
      cli_report(tensile_properties(r), o);
      return 0;
    } catch (const std::exception& e) {
      std::fprintf(stderr, "caps tensile: %s\n", e.what());
      return 1;
    }
  }
  if (cmd == "pull") {
    // caps pull FILE [--normal] [--distance 10] [--rate 5] [--spring 10] [--temp 300] [--eq-ps 5] [--surface 1] [--ff FF] [--csv OUT.csv]
    try {
      if (pos.empty()) return usage();
      System s = open_file(pos[0], o.count("--topology") ? o["--topology"] : "").frame(0);
      PullOptions po;
      po.normal = o.count("--normal");
      if (o.count("--axis")) po.axis = axis_of(o["--axis"]);
      if (o.count("--distance")) po.distance = std::stod(o["--distance"]);
      if (o.count("--rate")) po.rate = std::stod(o["--rate"]);
      if (o.count("--spring")) po.spring = std::stod(o["--spring"]);
      if (o.count("--temp")) po.temperature = std::stod(o["--temp"]);
      if (o.count("--eq-ps")) po.equilibrate_ps = std::stod(o["--eq-ps"]);
      if (o.count("--surface")) po.surface_mol = std::stoll(o["--surface"]);
      if (o.count("--seed")) po.seed = std::stoull(o["--seed"]);
      po.field = std::make_shared<ForceField>(cli_forcefield(s, o));
      electrostatics(po.energy, o);
      const PullResult r = run_pull(s, po);
      std::printf("%s\n", r.method.c_str());
      std::printf("peak force %.3f kcal/mol/Å at %.2f Å · %s %.3f MPa · work %.2f mJ/m²\n", r.peak_force, r.peak_displacement,
                  po.normal ? "peak normal stress" : "interfacial shear strength", r.strength, r.work);
      for (const auto& n : r.notes) std::printf("note: %s\n", n.c_str());
      if (o.count("--csv")) {
        std::ofstream f(o["--csv"]);
        f << "time_ps,displacement_A,force_kcal_mol_A,force_smoothed,temperature_K\n";
        for (size_t i = 0; i < r.curve.size(); ++i)
          f << r.curve[i].time_ps << "," << r.curve[i].displacement << "," << r.curve[i].force << "," << r.smooth[i] << "," << r.curve[i].temperature << "\n";
        std::printf("wrote %s\n", o["--csv"].c_str());
      }
      return 0;
    } catch (const std::exception& e) {
      std::fprintf(stderr, "caps pull: %s\n", e.what());
      return 1;
    }
  }
  if (cmd == "exchange") {
    // caps exchange arrhenius TABLE [TABLE2] [--reduced] [--from T --to T]
    //     Ea from relaxation times (columns T τ); a second table of times at the same temperatures: where the two cross
    // caps exchange modulus PTENSOR --cell-volume Å³ --temp K --dt FS [--fit PS] [--g-inf MPa] [--corr-ps PS]
    //     G(t), G'(ω), G''(ω) from pressure-tensor samples, atm, the last six columns xx yy zz xy xz yz (LAMMPS
    //     fix ave/time of pxx pyy pzz pxy pxz pyz in units real)
    // caps exchange memory PARTNERS --frame-ps PS
    //     bond lifetime: a row per frame, a column per exchangeable site, the partner's id (−1: free)
    // caps exchange bridges START NOW [--axis z] [--plane 0.5] [--bin 1]
    //     bonds across a cut: sides from START's positions, bonds and positions from NOW (the same atoms)
    // caps exchange heal HEALED INTACT [--modulus-strain 0.03] [--at-strain E] [--up-to E]
    //     healing efficiency from two stress–strain tables (strain as a fraction, stress)
    try {
      if (pos.empty()) return usage();
      const std::string sub = pos[0];
      std::vector<Property> props;
      auto column = [](const std::vector<std::vector<double>>& t, size_t c, const std::string& path) {
        std::vector<double> v;
        for (const auto& r : t) {
          if (r.size() <= c) throw std::invalid_argument(path + ": a row has fewer than " + std::to_string(c + 1) + " columns");
          v.push_back(r[c]);
        }
        return v;
      };
      if (sub == "arrhenius") {
        if (pos.size() < 2) return usage();
        const auto t1 = read_number_table(pos[1]);
        const auto T = column(t1, 0, pos[1]), tau = column(t1, 1, pos[1]);
        const bool reduced = o.count("--reduced") > 0;
        const double tmin = o.count("--from") ? std::stod(o["--from"]) : 0, tmax = o.count("--to") ? std::stod(o["--to"]) : 0;
        props = arrhenius_properties(fit_arrhenius(T, tau, reduced, tmin, tmax), "τ");
        if (pos.size() > 2) {
          const auto t2 = read_number_table(pos[2]);
          const auto T2 = column(t2, 0, pos[2]), tau2 = column(t2, 1, pos[2]);
          auto p2 = arrhenius_properties(fit_arrhenius(T2, tau2, reduced, tmin, tmax), "τ₂");
          p2[0].id = "arrhenius_2";
          if (T2 != T) throw std::invalid_argument("the two tables must list the same temperatures in the same order");
          const double tx = crossing_temperature(T, tau, tau2);
          if (std::isfinite(tx)) props[0].extra["temperature where τ and τ₂ cross"] = tx;
          else props[0].notes.push_back("τ and τ₂ do not cross over these temperatures");
          props.push_back(p2[0]);
        }
      } else if (sub == "modulus") {
        if (pos.size() < 2) return usage();
        if (!o.count("--cell-volume") || !o.count("--temp") || !o.count("--dt"))
          throw std::invalid_argument("modulus needs --cell-volume (Å³), --temp (K) and --dt (fs between samples)");
        const auto t = read_number_table(pos[1]);
        std::vector<std::array<double, 6>> p;
        for (const auto& r : t) {
          if (r.size() < 6) throw std::invalid_argument(pos[1] + ": a row has fewer than six columns (xx yy zz xy xz yz)");
          const size_t k = r.size() - 6;
          p.push_back({r[k], r[k + 1], r[k + 2], r[k + 3], r[k + 4], r[k + 5]});
        }
        StressRelaxationOptions so;
        so.volume = std::stod(o["--cell-volume"]), so.temperature = std::stod(o["--temp"]), so.dt_fs = std::stod(o["--dt"]);
        if (o.count("--fit")) so.fit_from_ps = std::stod(o["--fit"]);
        if (o.count("--g-inf")) so.g_inf = std::stod(o["--g-inf"]);
        if (o.count("--corr-ps")) so.integrate_ps = std::stod(o["--corr-ps"]);
        props = stress_relaxation_properties(stress_relaxation(p, so));
      } else if (sub == "memory") {
        if (pos.size() < 2) return usage();
        if (!o.count("--frame-ps")) throw std::invalid_argument("memory needs --frame-ps (time between the rows)");
        std::vector<std::vector<int64_t>> partners;
        for (const auto& r : read_number_table(pos[1])) {
          std::vector<int64_t> f;
          for (double v : r) f.push_back(v < 0 ? -1 : int64_t(std::llround(v)));
          partners.push_back(std::move(f));
        }
        props = bond_memory_properties(bond_memory(partners, std::stod(o["--frame-ps"])), "ps");
      } else if (sub == "bridges") {
        if (pos.size() < 3) return usage();
        const System start = load(pos[1], o), now = load(pos[2], o);
        if (start.atoms.size() != now.atoms.size()) throw std::invalid_argument("START and NOW must hold the same atoms");
        const int axis = o.count("--axis") ? axis_of(o["--axis"]) : 2;
        const double plane = o.count("--plane") ? std::stod(o["--plane"]) : 0.5, bin = o.count("--bin") ? std::stod(o["--bin"]) : 1.0;
        props = interface_properties(interface_bonds(now, sides_from_positions(start, axis, plane), axis, plane, bin), axis);
      } else if (sub == "heal") {
        if (pos.size() < 3) return usage();
        auto curve = [&](const std::string& path) {
          const auto t = read_number_table(path);
          return StressStrain{column(t, 0, path), column(t, 1, path)};
        };
        const double ms = o.count("--modulus-strain") ? std::stod(o["--modulus-strain"]) : 0.03,
                     at = o.count("--at-strain") ? std::stod(o["--at-strain"]) : 0, up = o.count("--up-to") ? std::stod(o["--up-to"]) : 0;
        props = healing_properties(healing_efficiency(curve(pos[1]), curve(pos[2]), ms, at, up));
      } else {
        auto exchange_lib = [] {
          for (const std::string& root : caps_roots())
            if (!root.empty() && std::filesystem::exists(root + "/data/exchange/library.json")) return load_exchange_library(root + "/data/exchange/library.json");
          throw std::runtime_error("the exchange library (data/exchange/library.json) was not found");
        };
        if (sub == "list") {
          const auto lib = exchange_lib();
          std::printf("%-26s %-13s %-6s %-22s %s\n", "id", "mechanism", "move", "scheme", "default barrier (kcal/mol) and sources");
          for (const auto& c : lib) {
            char bar[64] = "none";
            if (std::isfinite(c.default_kcal_per_mol)) std::snprintf(bar, sizeof bar, "%.4g", c.default_kcal_per_mol);
            std::printf("%-26s %-13s %-6s %-22s %s, %zu value%s from the literature\n", c.id.c_str(), c.mechanism.substr(0, 13).c_str(), c.move.c_str(),
                        c.scheme.empty() ? "(none yet)" : c.scheme.c_str(), bar, c.barriers.size(), c.barriers.size() == 1 ? "" : "s");
          }
          return 0;
        }
        if (sub == "show") {
          if (pos.size() < 2) return usage();
          const auto lib = exchange_lib();
          const auto& c = find_exchange_chemistry(lib, pos[1]);
          std::printf("%s (%s)\n  family: %s\n  mechanism: %s; move: %s; scheme: %s\n", c.name.c_str(), c.id.c_str(), c.family.c_str(), c.mechanism.c_str(),
                      c.move.c_str(), c.scheme.empty() ? "none yet" : c.scheme.c_str());
          for (const auto& b : c.barriers)
            std::printf("  %s %s: %s; %s%s%s [%s, doi:%s]%s%s\n", b.value.c_str(), b.unit.c_str(), b.kind.c_str(), b.method.c_str(), b.conditions.empty() ? "" : "; ",
                        b.conditions.c_str(), b.status.c_str(), b.source.c_str(), b.note.empty() ? "" : "\n    note: ", b.note.c_str());
          if (std::isfinite(c.default_kcal_per_mol)) std::printf("  default barrier: %.4g kcal/mol (%s)\n", c.default_kcal_per_mol, c.default_why.c_str());
          else std::printf("  default barrier: none (no molecular barrier read yet)\n");
          if (std::isfinite(c.capture_angstrom)) std::printf("  capture: %s %.4g A [doi:%s]\n", c.capture_pair.c_str(), c.capture_angstrom, c.capture_source.c_str());
          for (const auto& x : c.systems) std::printf("  system: %s\n", x.c_str());
          for (const auto& x : c.to_read) std::printf("  to read: %s\n", x.c_str());
          if (!c.notes.empty()) std::printf("  note: %s\n", c.notes.c_str());
          return 0;
        }
        // a chemistry's reaction template, from its scheme in data/reactions/library.json
        auto chemistry_template = [&](const ExchangeChemistry& c) {
          if (c.scheme.empty()) throw std::invalid_argument(c.id + " has no reaction scheme yet (see caps exchange show " + c.id + ")");
          for (const std::string& root : caps_roots()) {
            const std::string path = root + "/data/reactions/library.json";
            if (root.empty() || !std::filesystem::exists(path)) continue;
            std::ifstream in(path);
            std::stringstream ss;
            ss << in.rdbuf();
            const Json lib = Json::parse(ss.str());
            for (const auto& e : lib["reactions"].items()) {
              if (e.text("id") != c.scheme) continue;
              std::vector<std::string> rs, ps;
              for (const auto& r : e["reactants"].items()) rs.push_back(r.text("smiles"));
              for (const auto& p : e["products"].items()) ps.push_back(p.text("smiles"));
              return template_from_scheme(rs, ps, c.scheme);
            }
          }
          throw std::runtime_error("scheme " + c.scheme + " not found in data/reactions/library.json");
        };
        if (sub == "bond-react") {
          // LAMMPS fix bond/react for an exchange chemistry, its prob scaled by the Arrhenius factor of its barrier
          if (pos.size() < 2 || !o.count("-o") || !o.count("--chemistry")) return usage();
          System s = load(pos[1], o);
          const auto lib = exchange_lib();
          const auto& c = find_exchange_chemistry(lib, o["--chemistry"]);
          const ReactionTemplate t = chemistry_template(c);
          BondReactOptions bo;
          if (o.count("--temp")) bo.temperature = std::stod(o["--temp"]);
          if (o.count("--steps")) bo.steps = std::stoll(o["--steps"]);
          if (o.count("--seed")) bo.seed = std::stoull(o["--seed"]);
          if (o.count("--t-ref")) {
            bo.t_ref = std::stod(o["--t-ref"]);
            if (!std::isfinite(c.default_kcal_per_mol)) throw std::invalid_argument(c.id + " has no default barrier: --t-ref has nothing to scale");
            bo.barrier = c.default_kcal_per_mol;
          }
          auto field = [&](const System& x) {
            System y = x;
            return std::make_shared<const ForceField>(cli_forcefield(y, o, true));
          };
          const std::filesystem::path stem(o["-o"]);
          const auto rep = write_bond_react(s, {t}, field, stem.has_parent_path() ? stem.parent_path().string() : ".", stem.filename().string(), bo);
          std::printf("%s: %zu template variant%s; prob scaled by %.4g (barrier %.4g kcal/mol, %g K against %g K)\n", c.id.c_str(), rep.variants.size(),
                      rep.variants.size() == 1 ? "" : "s",
                      bo.barrier > 0 ? arrhenius_probability(bo.barrier, bo.temperature, bo.t_ref) : 1.0, bo.barrier, bo.temperature, bo.t_ref);
          for (const auto& f : rep.files) std::printf("wrote %s\n", f.c_str());
          return 0;
        }
        if (sub == "run") {
          // CAPS's own MD with exchange moves (caps/exchange_moves.hpp)
          if (pos.size() < 2 || !o.count("-o")) return usage();
          System s = load(pos[1], o);
          KgOptions ko;
          if (o.count("--k-theta")) ko.k_theta = std::stod(o["--k-theta"]);
          const bool kg = o.count("--kg") > 0;
          ForceField ff = kg ? kremer_grest_forcefield(s, ko) : cli_forcefield(s, o);
          ExchangeOptions xo;
          if (o.count("--chemistry")) {
            const auto lib = exchange_lib();
            xo = exchange_options_for(find_exchange_chemistry(lib, o["--chemistry"]));
            std::printf("chemistry %s: %s moves, barrier %.4g kcal/mol, capture %.4g A\n", o["--chemistry"].c_str(), to_string(xo.kind), xo.barrier, xo.capture);
          }
          if (o.count("--kind") || !o.count("--chemistry")) xo.kind = exchange_kind_from_string(o.count("--kind") ? o["--kind"] : "swap");
          if (xo.kind == ExchangeKind::Template) {
            if (!o.count("--chemistry")) throw std::invalid_argument("--kind template takes its template from --chemistry");
            const auto lib = exchange_lib();
            xo.templates = {chemistry_template(find_exchange_chemistry(lib, o["--chemistry"]))};
          }
          if (xo.kind != ExchangeKind::Template || o.count("--sites")) xo.group_a = cli_selection(s, o.count("--sites") ? o["--sites"] : "all");
          if (o.count("--partners")) xo.group_b = cli_selection(s, o["--partners"]);
          if (o.count("--labels")) {
            if (o["--labels"] != "chain" && o["--labels"] != "ends") throw std::invalid_argument("--labels takes chain or ends");
            xo.labels = chain_position_labels(s, o["--labels"] == "ends");
          }
          if (o.count("--capture")) xo.capture = std::stod(o["--capture"]);
          if (o.count("--barrier")) xo.barrier = std::stod(o["--barrier"]);
          if (o.count("--binding")) xo.binding = std::stod(o["--binding"]);
          if (o.count("--attempts")) xo.attempts = std::stoi(o["--attempts"]);
          if (o.count("--seed")) xo.seed = std::stoull(o["--seed"]);
          DynamicsOptions d;
          d.field = std::make_shared<ForceField>(ff);
          d.steps = o.count("--steps") ? std::stoll(o["--steps"]) : 10000;
          d.temperature = o.count("--temp") ? std::stod(o["--temp"]) : 300;
          if (o.count("--dt")) d.dt = std::stod(o["--dt"]);
          if (o.count("--seed")) d.seed = std::stoull(o["--seed"]);
          if (o.count("--cutoff")) d.energy.cutoff = std::stod(o["--cutoff"]);
          electrostatics(d.energy, o);
          precision(d.energy, o);
          if (o.count("--threads")) d.energy.threads = std::stoi(o["--threads"]);
          d.topology_every = o.count("--every") ? std::stoll(o["--every"]) : 10;
          auto st = std::make_shared<ExchangeStats>();
          d.topology = exchange_callback(xo, d.energy, st);
          DynamicsReport r;
          run_dynamics(s, d, &r);
          std::printf("%s moves every %lld steps: %lld tried, %lld accepted, %lld without a partner in reach, %lld the force field could not describe; "
                      "%.6g kcal/mol added by the accepted moves%s\n",
                      to_string(xo.kind), (long long)d.topology_every, (long long)st->attempts, (long long)st->accepted, (long long)st->no_candidate,
                      (long long)st->unparameterised, st->energy, st->local ? " (local energies)" : " (whole-system energies)");
          if (!st->note.empty()) std::printf("note: %s\n", st->note.c_str());
          for (const auto& n : r.notes) std::printf("%s\n", n.c_str());
          const ForceField after = kg ? kremer_grest_forcefield(s, ko) : cli_forcefield(s, o, true);
          save_structure(s, after, d.energy, o["-o"]);
          std::printf("wrote %s\n", o["-o"].c_str());
          return 0;
        }
        if (sub == "bond-swap") {
          // LAMMPS fix bond/swap for a Kremer–Grest system (write_kg_bond_swap)
          if (pos.size() < 2 || !o.count("-o")) return usage();
          const System s = load(pos[1], o);
          KgOptions ko;
          if (o.count("--k-theta")) ko.k_theta = std::stod(o["--k-theta"]);
          BondSwapExport e;
          if (o.count("--sites")) e.group = cli_selection(s, o["--sites"]);
          e.stickers = o.count("--stickers") > 0;
          e.from_both_ends = o.count("--both-ends") > 0;
          if (o.count("--every")) e.every = std::stoi(o["--every"]);
          if (o.count("--fraction")) e.fraction = std::stod(o["--fraction"]);
          if (o.count("--cutoff")) e.cutoff = std::stod(o["--cutoff"]);
          if (o.count("--seed")) e.seed = std::stoull(o["--seed"]);
          if (o.count("--steps")) e.run_steps = std::stod(o["--steps"]);
          write_kg_bond_swap(s, ko, e, o["-o"]);
          std::printf("wrote %s.data and %s.in (fix bond/swap, %s mode)\n", o["-o"].c_str(), o["-o"].c_str(), e.stickers ? "sticker" : "chain");
          return 0;
        }
        throw std::invalid_argument("caps exchange arrhenius | modulus | memory | bridges | heal | run | bond-swap");
      }
      cli_report(props, o);
      return 0;
    } catch (const std::exception& e) {
      std::fprintf(stderr, "caps exchange: %s\n", e.what());
      return 1;
    }
  }
  if (cmd == "tg") {
    // caps tg DATA -o OUT.data [--from 500 --to 200 --step 20 --ps 100] [--dt 1] [--pressure 1] [--barostat crescale|berendsen]
    //         [--seed 1] [--new-velocities] [--ff ...] [--json] [--csv DIR]
    // caps tg --fit RUN [RUN ...]    two-line fit of finished scans (a run folder, OUT.tg-points.json, a csv folder or a table
    //                               of T and density or specific volume); several are replicas: mean ± SD
    try {
      if (o.count("--fit")) {
        // one finished scan or several replicas (their folders, OUT.tg-points.json, csv folders or tables): each fitted, and
        // replicas pooled as mean ± SD
        std::vector<std::string> paths = {o["--fit"]};
        paths.insert(paths.end(), pos.begin(), pos.end());
        std::vector<CoolingScan> scans;
        for (const auto& p : paths) scans.push_back(read_cooling_scan(p, o.count("--volume") > 0));
        const auto props = fit_cooling_scans(scans);
        for (size_t k = 0; k < scans.size(); ++k) {
          const auto& sc = scans[k];
          std::printf("%s: %zu temperatures %g to %g K%s%s\n", sc.source.c_str(), sc.points.size(), sc.points.front().temperature, sc.points.back().temperature,
                      sc.log.empty() ? "" : (", schedule from " + sc.log).c_str(), sc.atoms ? (", " + std::to_string(sc.atoms) + " atoms").c_str() : "");
        }
        cli_report(props, o);
        return 0;
      }
      if (pos.empty() || !o.count("-o")) return usage();
      System s = load(pos[0], o);
      auto ff = std::make_shared<ForceField>(cli_forcefield(s, o));
      CoolingOptions co;
      co.field = ff;
      if (o.count("--cutoff")) co.energy.cutoff = std::stod(o["--cutoff"]);
      electrostatics(co.energy, o);
      precision(co.energy, o);
      if (o.count("--threads")) co.energy.threads = std::stoi(o["--threads"]);
      if (o.count("--from")) co.t_start = std::stod(o["--from"]);
      if (o.count("--to")) co.t_end = std::stod(o["--to"]);
      if (o.count("--step")) co.t_step = std::stod(o["--step"]);
      if (o.count("--ps")) co.ps_per_step = std::stod(o["--ps"]);
      if (o.count("--equilibrate")) co.equilibrate_ps = std::stod(o["--equilibrate"]);
      if (o.count("--dt")) co.dt = std::stod(o["--dt"]);
      if (o.count("--pressure")) co.pressure = std::stod(o["--pressure"]);
      if (o.count("--barostat")) co.barostat = barostat_from_string(o["--barostat"]);
      if (o.count("--seed")) co.seed = std::stoull(o["--seed"]);
      co.new_velocities = o.count("--new-velocities") || s.velocities.size() != s.atoms.size();
      const bool quiet = o.count("--quiet");
      // checkpoints after each temperature (a cluster job near its time limit stops with exit 75 and resume.txt)
      install_stop_signals();
      const std::string tg_cp = std::filesystem::path(o["-o"]).replace_extension(".restart.data").string();
      const std::string tg_state = std::filesystem::path(o["-o"]).replace_extension(".tg-points.json").string();
      if (o.count("--resume-points")) {
        std::ifstream pf(o["--resume-points"]);
        std::stringstream ss;
        ss << pf.rdbuf();
        for (const auto& q : Json::parse(ss.str())["points"].items()) {
          CoolingPoint p;
          p.temperature = q.num("T", 0), p.density = q.num("density", 0), p.density_err = q.num("density_err", 0);
          p.specific_volume = q.num("specific_volume", 0), p.potential = q.num("potential", 0);
          co.done.push_back(p);
        }
      }
      co.after_point = [&](const System& now, const std::vector<CoolingPoint>& pts) {
        save_structure(now, *ff, co.energy, tg_cp);
        Json j = Json::object(), arr = Json::array();
        for (const auto& p : pts) {
          Json q = Json::object();
          q["T"] = p.temperature, q["density"] = p.density, q["density_err"] = p.density_err, q["specific_volume"] = p.specific_volume, q["potential"] = p.potential;
          arr.push_back(q);
        }
        j["points"] = arr;
        std::ofstream(tg_state) << j.dump(1) << "\n";
        std::ofstream("resume.txt") << resume_line("caps tg", tg_cp, o, {{"--resume-points", tg_state}}) << "\n";
      };
      int last = -1;
      co.progress = [&](const ThermoRow& r, int k, int n) {
        if (!quiet && k != last) {
          if (k < 0) std::printf("equilibrating at %.1f K\n", r.target_temperature);
          else std::printf("step %d of %d: %.1f K\n", k + 1, n, r.target_temperature);
          last = k;
        }
        return true;
      };
      CoolingResult r;
      try {
        r = run_cooling(s, co);
      } catch (const DynamicsInterrupted& ex) {
        if (!std::filesystem::exists("resume.txt")) std::ofstream("resume.txt") << resume_line("caps tg", pos[0], o, {}) << "\n";   // before the first temperature
        std::ifstream rf("resume.txt");
        std::string again;
        std::getline(rf, again);
        std::printf("stopped at step %lld (a stop was asked for); continue with\n  %s\n", static_cast<long long>(ex.step), again.c_str());
        return 75;
      }
      std::filesystem::remove("resume.txt");   // finished: nothing to continue
      if (!quiet)
        for (const auto& p : r.points) std::printf("%8.1f K  %.5f ± %.5f g/cm³\n", p.temperature, p.density, p.density_err);
      save_structure(s, *ff, co.energy, o["-o"]);
      cli_report(cooling_properties(r), o);
      return 0;
    } catch (const std::exception& e) {
      std::fprintf(stderr, "caps tg: %s\n", e.what());
      return 1;
    }
  }
  if (cmd == "frames") {
    // caps frames FILE OUT.lammpstrj|dcd|xyz|pdb|gro|trr [--topology DATA] [--first N] [--last N] [--stride N] [--timestep-fs 1]
    // A long trajectory thinned as it is read: only the frames kept are held and written.
    if (pos.size() < 2) return usage();
    OpenProgress rd;
    if (o.count("--first")) rd.frames.first = std::stoul(o["--first"]);
    if (o.count("--last")) rd.frames.last = std::stoul(o["--last"]);
    if (o.count("--stride")) rd.frames.stride = std::max<size_t>(1, std::stoul(o["--stride"]));
    const Trajectory t = open_file(pos[0], o.count("--topology") ? o["--topology"] : "", rd);
    write_trajectory(t, pos[1], o.count("--timestep-fs") ? std::stod(o["--timestep-fs"]) : 1.0);
    std::printf("wrote %s: %zu frames (%zu read)\n", pos[1].c_str(), t.frames(), t.frames_read ? t.frames_read : t.frames());
    return 0;
  }
  if (cmd == "analyze") {
    // caps analyze FILE [--topology DATA] --props density,rdf,... [--first --last --stride] [--frame-ps | --timestep-fs]
    //              [--pair C-C] [--inter] [--rmax] [--dr] [--qmax] [--dq] [--qdirect] [--probe] [--grid] [--ff FF.json [--typing R] [--charges]]
    //              [--json OUT.json] [--csv DIR]
    try {
      if (pos.empty()) return usage();
      // only the frames analysed are read (a dump's others passed over unread): the analysis then takes them all
      OpenProgress rd;
      if (o.count("--first")) rd.frames.first = size_t(std::max(0L, std::stol(o["--first"])));
      if (o.count("--last") && std::stol(o["--last"]) >= 0) rd.frames.last = size_t(std::stol(o["--last"]));
      if (o.count("--stride")) rd.frames.stride = size_t(std::max(1L, std::stol(o["--stride"])));
      Trajectory t = open_file(pos[0], o.count("--topology") ? o["--topology"] : "", rd);
      if (t.frames_read == 0 && !rd.frames.all()) t.select(rd.frames);   // a file read whole (one frame): the same frames
      AnalyzeOptions ao;
      if (o.count("--frame-ps")) ao.frame_ps = std::stod(o["--frame-ps"]) * double(rd.frames.stride);   // between the frames kept
      if (o.count("--timestep-fs")) ao.timestep_fs = std::stod(o["--timestep-fs"]);
      if (o.count("--blocks")) ao.blocks = std::stoi(o["--blocks"]);
      if (o.count("--rmax")) ao.rdf_rmax = std::stod(o["--rmax"]);
      if (o.count("--dr")) ao.rdf_dr = std::stod(o["--dr"]);
      if (o.count("--qmax")) ao.qmax = std::stod(o["--qmax"]);
      if (o.count("--dq")) ao.dq = std::stod(o["--dq"]);
      if (o.count("--qdirect")) ao.q_direct = std::stod(o["--qdirect"]);
      if (o.count("--probe")) ao.probe = std::stod(o["--probe"]);
      if (o.count("--grid")) ao.grid = std::stod(o["--grid"]);
      if (o.count("--fit")) {   // --fit 0.2:0.5
        const auto f = o["--fit"];
        ao.fit_from = std::stod(f.substr(0, f.find(':')));
        ao.fit_to = std::stod(f.substr(f.find(':') + 1));
      }
      ao.inter_only = o.count("--inter") > 0;
      if (o.count("--pair")) {
        const auto pr = o["--pair"];
        const auto dash = pr.find('-');
        ao.elem_a = pr.substr(0, dash) == "all" ? 0 : element_from_symbol(pr.substr(0, dash));
        const std::string b = dash == std::string::npos ? pr : pr.substr(dash + 1);
        ao.elem_b = b == "all" ? 0 : element_from_symbol(b);
      }
      std::vector<std::string> ids;
      {
        std::stringstream ss(o.count("--props") ? o["--props"] : "density,rdf,rg,ree,cn,persistence,msd,diffusion,relaxation,ffv");
        for (std::string w; std::getline(ss, w, ',');)
          if (!w.empty()) ids.push_back(w);
      }
      // cohesive energy and fluctuation elastic constants: the force field of the structure (frame 0, molecules whole)
      ForceField ff;
      auto has = [&](const char* k) { return std::find(ids.begin(), ids.end(), k) != ids.end(); };
      if (o.count("--zbin")) ao.zbin = std::stod(o["--zbin"]);
      if (o.count("--exclude-mol")) ao.exclude_mol = std::stoll(o["--exclude-mol"]);
      if (o.count("--axis")) {
        const std::string a = o["--axis"];
        if (a == "a" || a == "x") ao.axis = 0;
        else if (a == "b" || a == "y") ao.axis = 1;
        else if (a == "c" || a == "z") ao.axis = 2;
        else throw std::invalid_argument("--axis: x, y or z (a, b or c)");
      }
      if (o.count("--surface")) {
        std::string t = o["--surface"];
        for (auto& c : t) if (c == ',') c = ' ';
        std::istringstream is(t);
        for (std::string w; is >> w;) {
          const auto dash = w.find('-', 1);
          const long long a = std::stoll(w.substr(0, dash)), b = dash == std::string::npos ? a : std::stoll(w.substr(dash + 1));
          if (b < a) throw std::invalid_argument("--surface: bad range " + w);
          for (long long m = a; m <= b; ++m) ao.surface_mols.push_back(m);
        }
      }
      if (has("ced") || has("delta") || has("cij_fluct") || has("adhesion")) {
        System s0 = t.frame(0);
        if (!s0.unwrapped) make_molecules_whole(s0);
        ff = cli_forcefield(s0, o);
        ao.ff = &ff;
        // --cutoff, else the force field's own cut-off (OPLS-AA 2024: 12 Å), as relax and md take it
        run_cutoff(ao.energy, ff, o, false);
      }
      if (o.count("--temp")) ao.temperature = std::stod(o["--temp"]);
      if (o.count("--viscosity")) ao.viscosity = std::stod(o["--viscosity"]);   // mPa·s: Yeh–Hummer correction of D (with --temp)
      const auto props = analyze(t, ids, ao);
      const auto fr = analysis_frames(t, ao);
      std::printf("%s · %zu frames used of %zu · %zu atoms\n", pos[0].c_str(), fr.size(), t.frames(), t.topology.atoms.size());
      cli_report(props, o);
      return 0;
    } catch (const std::exception& e) {
      std::fprintf(stderr, "caps analyze: %s\n", e.what());
      return 1;
    }
  }
  if (cmd == "react" && o.count("--list-templates")) {
    for (const auto& n : builtin_template_names()) std::printf("%s", builtin_template(n).c_str());
    return 0;
  }
  if (pos.empty()) return usage();
  try {
    System s = load(pos[0], o);
    if (cmd == "react") {
      ReactOptions r;
      std::string text;
      // --template may be given several times: collect them from argv
      for (int i = 2; i + 1 < argc; ++i)
        if (std::string(argv[i]) == "--template") {
          const std::string t = argv[i + 1];
          if (std::ifstream f(t); f && t.find('.') != std::string::npos) { std::stringstream ss; ss << f.rdbuf(); text += ss.str() + "\n"; }
          else text += builtin_template(t) + "\n";
        }
      if (text.empty()) text = builtin_template("cc_crosslink");
      r.templates = parse_templates(text);
      if (o.count("--capture")) for (auto& t : r.templates) t.capture = std::stod(o["--capture"]);
      if (o.count("--per-cycle")) r.max_per_cycle = std::stoi(o["--per-cycle"]);
      if (o.count("--cycles")) r.max_cycles = std::stoi(o["--cycles"]);
      if (o.count("--target")) r.target_conversion = std::stod(o["--target"]);
      if (o.count("--md-ps")) r.md_ps = std::stod(o["--md-ps"]);
      if (o.count("--temp")) r.temperature = std::stod(o["--temp"]);
      if (o.count("--seed")) r.seed = std::stoull(o["--seed"]);
      r.relax = !o.count("--no-relax");
      if (o.count("--relax-iterations")) r.relax_iterations = std::stoi(o["--relax-iterations"]);
      // the options of the Studio's React page and caps.react (the same ReactOptions; caps_react in the C API)
      r.during_md = o.count("--during-md") > 0;
      r.between_chains = o.count("--between-chains") > 0;
      r.keep_byproducts = o.count("--keep-byproducts") > 0;
      r.remove_unreacted = o.count("--remove-unreacted") > 0;
      if (o.count("--sites-per-chain")) r.sites_per_chain = std::max(0, std::stoi(o["--sites-per-chain"]));
      r.auto_capture = o.count("--auto-capture") > 0;
      if (o.count("--capture-max")) r.capture_max = std::stod(o["--capture-max"]);
      if (o.count("--capture-step")) r.capture_step = std::stod(o["--capture-step"]);
      if (o.count("--weights")) {   // one relative rate per template, in order: each pick chooses a template by weight
        for (std::stringstream ss(o["--weights"]); ss.good();) {
          std::string w;
          std::getline(ss, w, ',');
          if (!w.empty()) r.weights.push_back(std::stod(w));
        }
        if (r.weights.size() != r.templates.size())
          throw std::invalid_argument("--weights: " + std::to_string(r.weights.size()) + " given for " + std::to_string(r.templates.size()) + " templates");
        for (double w : r.weights) if (w <= 0) throw std::invalid_argument("--weights must be > 0");
        r.selection = 1;
      }
      if (o.count("--crosslinks")) {   // KIND:VALUE — links, per-chain, density (mol/m³), mc (g/mol), degree (%)
        const std::string v = o["--crosslinks"];
        const auto colon = v.find(':');
        if (colon == std::string::npos) throw std::invalid_argument("--crosslinks KIND:VALUE (links, per-chain, density, mc or degree)");
        std::string kind = v.substr(0, colon);
        for (auto& ch : kind) if (ch == '_') ch = '-';
        if (kind == "links" || kind == "crosslinks") r.target = ReactTarget::Crosslinks;
        else if (kind == "per-chain") r.target = ReactTarget::PerChain;
        else if (kind == "density") r.target = ReactTarget::Density;
        else if (kind == "mc") r.target = ReactTarget::Mc;
        else if (kind == "degree") r.target = ReactTarget::DegreePercent;
        else throw std::invalid_argument("--crosslinks: kind " + kind + " (links, per-chain, density, mc or degree)");
        r.target_value = std::stod(v.substr(colon + 1));
      }
      r.energy.coulomb = !o.count("--no-coulomb");
      if (o.count("--cutoff")) r.energy.cutoff = std::stod(o["--cutoff"]);
      if (!o.count("-o")) return usage();
      if (o.count("--insert")) {   // curatives first: e.g. --insert SS --count 40 (H–S–S–H sulfur donors for sulfur_allylic)
        BuildOptions bo;
        bo.forcefield = "uff";
        const BuildResult br = build_molecule(o["--insert"], bo);
        PackOptions po;
        if (o.count("--tolerance")) po.tolerance = std::stod(o["--tolerance"]);
        PackReport pr;
        System guest = br.system;
        guest.title = o["--insert"];
        s = insert_molecules(s, guest, o.count("--count") ? std::stoi(o["--count"]) : 10, po, &pr);
        for (const auto& n : pr.notes) std::printf("%s\n", n.c_str());
      }
      const bool quiet = o.count("--quiet");
      // --ff [--overlay --typing --charges]: that force field types and relaxes every state of the network (as the Studio's
      // assigned force field; a complete assignment is needed from the start), else the built-in default (GAFF / UFF)
      std::shared_ptr<const ForceField> field_now;
      if (o.count("--ff")) {
        System s0 = s;
        field_now = std::make_shared<ForceField>(cli_forcefield(s0, o, quiet));
        const auto opts = o;
        r.retype = [opts](const System& x) {
          auto oo = opts;
          System y = x;
          return std::shared_ptr<const ForceField>(std::make_shared<ForceField>(cli_forcefield(y, oo, true)));
        };
        r.field_name = field_now->name;
        run_cutoff(r.energy, *field_now, o, quiet);   // the force field's own cut-off without --cutoff, as relax and md
      }
      electrostatics(r.energy, o);
      if (!quiet) std::printf("%6s %6s %6s %11s %9s %11s %12s %13s\n", "cycle", "react", "total", "conversion", "clusters", "largest %", "reduced Mw", "E/kcal·mol⁻¹");
      r.progress = [&](const CycleRow& c) {
        if (!quiet)
          std::printf("%6d %6d %6d %11.4f %9d %11.1f %12.1f %13.1f\n", c.cycle, c.reactions, c.total, c.conversion, c.clusters.clusters,
                      100 * c.clusters.largest_fraction, c.clusters.reduced_mw, c.energy);
        std::fflush(stdout);
        return true;
      };
      ReactReport rep;
      react(s, r, &rep);
      const std::string out = o["-o"];
      auto ends = [&](const char* e) { return out.size() > 4 && out.substr(out.size() - 4) == e; };
      if (ends(".pdb")) write_pdb(s, out);
      else if (ends(".xyz")) write_xyz(s, out);
      else if (r.retype) write_lammps_data_ff(s, *r.retype(s), r.energy, out);   // the network with its force field re-assigned
      else if (r.relax) write_lammps_data_ff(s, default_forcefield(s), r.energy, out);
      else write_lammps_data(s, out);
      for (const auto& n : rep.notes) std::printf("%s\n", n.c_str());
      if (o.count("--summary")) std::ofstream(o["--summary"]) << react_summary_json(rep, s) << "\n";
      if (const SulfurCensus sc = sulfur_census(s); sc.sulfur_atoms && !quiet) {
        std::string line = "sulfur:";
        for (const auto& [kind, m] : {std::pair{"bridges", &sc.bridges}, {"cyclic", &sc.cyclic}, {"pendant", &sc.pendant}, {"free", &sc.free}}) {
          if (m->empty()) continue;
          line += std::string(" ") + kind;
          for (const auto& [rank, n] : *m) line += " S" + std::to_string(rank) + "×" + std::to_string(n);
          line += ";";
        }
        std::printf("%s mean bridge rank %.3f (%d S atoms)\n", line.c_str(), sc.mean_bridge_rank, sc.sulfur_atoms);
      }
      if (o.count("--fa") && o.count("--fb")) {
        const double ratio = o.count("--ratio") ? std::stod(o["--ratio"]) : 1.0;
        const double ac = flory_stockmayer(ratio, std::stod(o["--fa"]), std::stod(o["--fb"]));
        std::printf("Flory–Stockmayer: α_c = 1/√(r (fA−1)(fB−1)) = %.3f (r %.2f, fA %s, fB %s)%s\n", ac, ratio, o["--fa"].c_str(), o["--fb"].c_str(),
                    rep.gel_conversion >= 0 ? (" · simulated " + std::to_string(rep.gel_conversion).substr(0, 5)).c_str() : " · simulated gel point not reached");
      }
      std::printf("wrote %s\n", out.c_str());
      return 0;
    }
    if (cmd == "contacts") {
      const double tol = o.count("--tol") ? std::stod(o["--tol"]) : 2.0;
      const bool per = !o.count("--no-pbc") && s.cell.valid();
      if (o.count("--molecule-size")) {   // molecules as consecutive blocks (for files whose bonds cannot be trusted)
        const int k = std::stoi(o["--molecule-size"]);
        if (k <= 0 || s.atoms.size() % size_t(k) != 0) throw std::runtime_error("--molecule-size must divide the atom count");
        for (size_t i = 0; i < s.atoms.size(); ++i) s.atoms[i].mol = int64_t(i / size_t(k) + 1);
        s.has_mol = true;
      } else if (s.atoms.size() > 99999 && s.source_format == "pdb") {
        std::fprintf(stderr, "note: PDB files over 99 999 atoms cannot number CONECT records reliably; molecules may be wrong "
                             "(use --molecule-size N for one-species files)\n");
      }
      int nm = 0;
      s.molecules(&nm);
      const auto [dmin, close] = intermolecular_contacts(s, tol, per);
      std::printf("%zu atoms, %d molecules, %s · smallest intermolecular distance %.4f Å · %d pairs closer than %.2f Å\n", s.atoms.size(), nm,
                  per ? "periodic (minimum image)" : "not periodic", dmin, close, tol);
      return close == 0 ? 0 : 1;
    }
    if (cmd == "bundle") {   // write a figure bundle: caps bundle FILE --steps S.json -o OUT.caps-bundle.zip [--include-input]
      const std::string top = o.count("--topology") ? o["--topology"] : "";
      const Trajectory t = open_file(pos[0], top);
      std::string text = o.count("--steps") ? o["--steps"] : "[]";
      if (!text.empty() && text[0] != '[' && text[0] != '{') {
        std::ifstream f(text);
        text.assign(std::istreambuf_iterator<char>(f), {});
      }
      BundleOptions bo;
      bo.input = pos[0];
      bo.topology = top;
      bo.include_input = o.count("--include-input") > 0;
      bo.frame = o.count("--frame") ? std::stoi(o["--frame"]) : 0;
      const std::string out = o.count("-o") ? o["-o"] : "figure.caps-bundle.zip";
      bo.name = std::filesystem::path(out).stem().stem().string();
      const auto files = bundle_files(t, pipeline_from_json(Json::parse(text)), bo, Camera{}, RenderOptions{});
      write_bundle(out, files);
      for (const auto& f : files) std::printf("%-40s %8zu bytes  %s\n", f.name.c_str(), f.bytes.size(), sha256_hex(f.bytes).substr(0, 12).c_str());
      std::printf("wrote %s\n", out.c_str());
      return 0;
    }
    if (cmd == "pipeline") {   // visualize pipeline: steps from a JSON file (or --steps '[…]'), on one frame
      const Trajectory t = open_file(pos[0], o.count("--topology") ? o["--topology"] : "");
      std::string text = o.count("--steps") ? o["--steps"] : "[]";
      if (!text.empty() && text[0] != '[' && text[0] != '{') {
        std::ifstream f(text);
        if (!f) throw std::runtime_error("cannot read " + text);
        text.assign(std::istreambuf_iterator<char>(f), {});
      }
      const int fr = o.count("--frame") ? std::stoi(o["--frame"]) : 0;
      if (fr < 0 || size_t(fr) >= t.frames()) throw std::runtime_error("frame out of range");
      const System f0 = t.frame(size_t(fr));
      Pipeline pl = text.rfind("caps_pipeline", 0) == 0 || text.find("\ncaps_pipeline") != std::string::npos ? pipeline_from_yaml(text) : pipeline_from_json(Json::parse(text));
      if (o.count("--branch")) pl.branch = o["--branch"];
      const auto st = run_pipeline(f0, pl, fr, t.timesteps.empty() ? 0 : t.timesteps[size_t(fr)], &t);
      for (size_t k = st.steps.size(); k-- > 0;)
        std::printf("%-8s %-22s %s\n", st.steps[k].level.c_str(), st.steps[k].title.c_str(), st.steps[k].summary.c_str());
      for (const auto& [k, v] : st.attributes) std::printf("  %-34s %.6g\n", k.c_str(), v);
      if (o.count("--table")) {
        for (const auto& tb : st.tables) {
          if (tb.name != o["--table"]) continue;
          for (size_t c = 0; c < tb.columns.size(); ++c) std::printf("%s%s", c ? "," : "", tb.columns[c].c_str());
          std::printf("\n");
          for (const auto& r : tb.rows) {
            for (size_t c = 0; c < r.size(); ++c) std::printf("%s%.6g", c ? "," : "", r[c]);
            std::printf("\n");
          }
        }
      }
      if (o.count("--particles")) std::printf("%s\n", particles_json(st, o["--particles"], 0, 50).dump(2).c_str());
      if (!pl.outputs.empty() || o.count("--out"))   // the pipeline's outputs block
        for (const auto& line : write_pipeline_outputs(st, pl, o.count("--out") ? o["--out"] : ".")) std::printf("%s\n", line.c_str());
      int errors = 0;
      for (const auto& s2 : st.steps) errors += s2.level == "error";
      return errors ? 2 : 0;
    }
    if (cmd == "check") {   // file checks: what was found, what was done, what to change
      const Trajectory t = open_file(pos[0], o.count("--topology") ? o["--topology"] : "");
      const auto checks = file_checks(t);
      int bad = 0;
      for (const auto& c : checks) {
        std::printf("%-5s  %s\n       %s\n", c.level.c_str(), c.title.c_str(), c.detail.c_str());
        bad += c.level == "error";
      }
      if (o.count("--report")) {
        std::ofstream f(o["--report"]);
        f << file_checks_text(checks, pos[0]);
        std::printf("wrote %s\n", o["--report"].c_str());
      }
      return bad ? 2 : 0;
    }
    if (cmd == "info") {
      int nm = 0;
      s.molecules(&nm);
      double q = 0;
      for (const auto& a : s.atoms) q += a.charge;
      std::map<std::string, int> counts;
      for (const auto& a : s.atoms) counts[element(a.element).symbol]++;
      std::printf("file        %s\nformat      %s\natoms       %zu\nbonds       %zu (%s)\nmolecules   %d\n", pos[0].c_str(), s.source_format.c_str(), s.atoms.size(),
                  s.bonds.size(), s.bonds_from_file ? "from file" : "perceived", nm);
      std::printf("elements   ");
      for (const auto& [e, c] : counts) std::printf(" %s %d", e.c_str(), c);
      std::printf("\n");
      if (s.cell.valid())
        std::printf("cell        %.4f %.4f %.4f Å · volume %.1f Å³\ndensity     %.4f g/cm³\n", norm(s.cell.a), norm(s.cell.b), norm(s.cell.c), s.cell.volume(), s.density());
      std::printf("mass        %.3f g/mol\ncharge      %+.3e e\n", s.total_mass(), q);
      for (const auto& n : s.notes) std::printf("note        %s\n", n.c_str());
      return 0;
    }
    if (cmd == "render") {
      if (!o.count("-o")) return usage();
      const std::string out = o["-o"];
      Camera cam;
      RenderOptions opt;
      if (o.count("--yaw")) cam.yaw = std::stod(o["--yaw"]) * M_PI / 180;
      if (o.count("--pitch")) cam.pitch = std::stod(o["--pitch"]) * M_PI / 180;
      if (o.count("--zoom")) cam.zoom = std::stod(o["--zoom"]);
      cam.perspective = o.count("--perspective");
      if (o.count("--size")) std::sscanf(o["--size"].c_str(), "%dx%d", &opt.width, &opt.height);
      const std::string bg = o.count("--bg") ? o["--bg"] : "dark";
      opt.background = bg == "white" ? Background::White : bg == "transparent" ? Background::Transparent : Background::Dark;
      const std::string c = o.count("--colour") ? o["--colour"] : "molecule";
      opt.colour_by = c == "element" ? ColourBy::Element : c == "type" ? ColourBy::Type : c == "distance" ? ColourBy::Property : ColourBy::Molecule;
      if (opt.colour_by == ColourBy::Property) {
        const auto sh = molecule_shapes(s);
        const auto mol = s.molecules();
        for (size_t i = 0; i < s.atoms.size(); ++i) opt.property.push_back(norm(s.atoms[i].pos - sh[mol[i]].com));
      }
      const std::string st = o.count("--style") ? o["--style"] : "ball";
      opt.style = st == "space" ? Style::SpaceFilling : st == "sticks" ? Style::Sticks : st == "noh" ? Style::NoHydrogens : st == "backbone" ? Style::Backbone : Style::BallAndStick;
      opt.show_cell = !o.count("--no-cell");
      // --wrap: atoms folded into the cell as in the simulation box; bonds across the boundary are not drawn; each chain
      // keeps its colour (taken from the whole molecules, before the bonds are cut)
      if (o.count("--wrap") && s.cell.valid()) {
        if (opt.colour_by == ColourBy::Molecule) {
          const auto mol = s.molecules();
          opt.colours.resize(s.atoms.size());
          for (size_t i = 0; i < s.atoms.size(); ++i) {
            unsigned c = molecule_colour(mol[i]);
            if (s.atoms[i].element == 1) {   // hydrogens lighter, as the renderer draws them
              unsigned r = c >> 16 & 255, g = c >> 8 & 255, b = c & 255;
              r += unsigned(0.55 * (255 - r)), g += unsigned(0.55 * (255 - g)), b += unsigned(0.55 * (255 - b));
              c = r << 16 | g << 8 | b;
            }
            opt.colours[i] = c;
          }
        }
        for (auto& a : s.atoms) a.pos = s.cell.wrap(a.pos);
        std::erase_if(s.bonds, [&](const Bond& b) { return norm(s.atoms[b.i].pos - s.atoms[b.j].pos) > 4.0; });
      }
      if (out.size() > 4 && out.substr(out.size() - 4) == ".svg") {
        std::ofstream f(out);
        f << render_svg(s, cam, opt);
      } else {
        Renderer r;
        write_png(r.render(s, cam, opt), out);
      }
      std::printf("wrote %s (%dx%d, %s background)\n", out.c_str(), opt.width, opt.height, bg.c_str());
      return 0;
    }
    if (cmd == "shape") {
      std::printf("mol  atoms      mass(g/mol)   Rg(Å)   kappa2\n");
      for (const auto& m : molecule_shapes(s)) std::printf("%4d %6d %14.3f %8.3f %8.4f\n", m.molecule + 1, m.atoms, m.mass, m.rg, m.kappa2);
      return 0;
    }
    if (cmd == "field") {
      const ForceField ff = default_forcefield(s);
      std::printf("force field %s\n", ff.name.c_str());
      for (const auto& n : ff.notes) std::printf("  %s\n", n.c_str());
      std::map<std::string, std::string> why;
      for (size_t i = 0; i < ff.atom_type.size(); ++i) why[ff.atom_type[i]] = ff.why[i];
      for (const auto& [t, w] : why) std::printf("  %-3s %s\n", t.c_str(), w.c_str());
      EnergyOptions eo;
      if (o.count("--cutoff")) eo.cutoff = std::stod(o["--cutoff"]);
      if (o.count("--no-tail")) eo.tail = false;
      electrostatics(eo, o);
      Evaluator ev(ff, eo);
      std::vector<double> x, f;
      for (const auto& a : s.atoms) x.insert(x.end(), a.pos.begin(), a.pos.end());
      const EnergyTerms e = ev.compute(x, s.cell, f);
      const int dg = o.count("--digits") ? std::clamp(std::stoi(o["--digits"]), 0, 10) : 2;
      std::printf("energy (kcal/mol): bond %.*f  angle %.*f  dihedral %.*f  improper %.*f  vdW %.*f  Coulomb %.*f  total %.*f\n", dg, e.bond, dg, e.angle,
                  dg, e.dihedral, dg, e.improper, dg, e.vdw, dg, e.coulomb, dg, e.total());
      std::printf("largest force %.3f kcal/mol/Å", max_force(f));
      if (s.cell.valid()) std::printf(" · pressure %.0f atm (0 K virial)", pressure_atm(e.virial, s.cell.volume()));
      std::printf("\n");
      if (o.count("--forces")) {
        std::ofstream fo(o["--forces"]);
        char b[128];
        for (size_t i = 0; i < s.atoms.size(); ++i) {
          std::snprintf(b, sizeof b, "%zu %.8f %.8f %.8f\n", i + 1, f[3 * i], f[3 * i + 1], f[3 * i + 2]);
          fo << b;
        }
      }
      return 0;
    }
    if (cmd == "relax") {
      if (!o.count("-o")) return usage();
      RelaxOptions r;
      if (o.count("--method")) r.method = minimiser_from_string(o["--method"]);
      if (o.count("--ftol")) r.ftol = std::stod(o["--ftol"]);
      if (o.count("--iterations")) r.max_iterations = std::stoi(o["--iterations"]);
      if (o.count("--density")) r.target_density = std::stod(o["--density"]);
      if (o.count("--step")) r.compress_step = std::stod(o["--step"]);
      if (o.count("--pressure")) r.pressure = std::stod(o["--pressure"]);
      r.relax_box = o.count("--box-relax");
      r.pushoff = !o.count("--no-pushoff");
      r.energy.coulomb = !o.count("--no-coulomb");
      electrostatics(r.energy, o);
      r.field = std::make_shared<ForceField>(cli_forcefield(s, o, o.count("--quiet") > 0));
      run_cutoff(r.energy, *r.field, o, o.count("--quiet") > 0);
      if (o.count("--fix-mol")) {   // hold one molecule in place (the substrate of an interface is molecule 1)
        const int64_t m = std::stoll(o["--fix-mol"]);
        r.fixed.assign(s.atoms.size(), 0);
        size_t nf = 0;
        for (size_t i = 0; i < s.atoms.size(); ++i)
          if (s.atoms[i].mol == m) r.fixed[i] = 1, ++nf;
        std::printf("holding %zu atoms of molecule %lld in place\n", nf, static_cast<long long>(m));
      }
      const bool quiet = o.count("--quiet");
      std::string last;
      r.progress = [&](const RelaxProgress& p) {
        if (!quiet && p.stage != last) {
          if (!last.empty()) std::fprintf(stderr, "\n");
          last = p.stage;
        }
        if (!quiet)
          std::fprintf(stderr, "\r  [%d/%d] %-34s it %5d  E %12.1f  |F|max %9.3f", p.stage_index, p.stages, p.stage.c_str(), p.iteration,
                       p.energy, p.fmax);
        return true;
      };
      RelaxReport rep;
      relax(s, r, &rep);
      if (!quiet) std::fprintf(stderr, "\n");
      const std::string out = o["-o"];
      auto ends = [&](const char* e) { return out.size() > 4 && out.substr(out.size() - 4) == e; };
      if (ends(".pdb")) write_pdb(s, out);
      else if (ends(".gro")) write_gro(s, out);   // same atoms and order as a GROMACS topology read with --topology
      else if (ends(".xyz")) write_xyz(s, out);
      else if (ends("mol2")) write_mol2(s, out);
      else if (ends(".car")) write_car(s, out);
      else if (const std::string why = write_lammps_data_or_structure(s, r.field ? *r.field : default_forcefield(s), r.energy, out); !why.empty()) std::fprintf(stderr, "%s: %s\n", out.c_str(), why.c_str());
      std::printf("%s\n", rep.field.c_str());
      for (const auto& n : rep.notes) std::printf("%s\n", n.c_str());
      std::printf("starting energy (kcal/mol): bond %.6f  angle %.6f  dihedral %.6f  improper %.6f  vdW %.6f  Coulomb %.6f  total %.6f\n", rep.initial.bond,
                  rep.initial.angle, rep.initial.dihedral, rep.initial.improper, rep.initial.vdw, rep.initial.coulomb, rep.initial.total());
      for (const auto& st : rep.stages)
        std::printf("  %-34s %5d it  E %12.1f  |F|max %8.3f  (%s)\n", st.name.c_str(), st.iterations, st.energy, st.fmax, st.stopped_by.c_str());
      std::printf("%d iterations, %d energy and force evaluations\n", rep.iterations, rep.evaluations);
      std::printf("wrote %s\n", out.c_str());
      return rep.converged ? 0 : 3;
    }
    if (cmd == "cbmc") {
      if (!o.count("-o")) return usage();
      CbmcOptions c;
      if (o.count("--moves")) c.moves = std::stoi(o["--moves"]);
      if (o.count("--trials")) c.trials = std::stoi(o["--trials"]);
      if (o.count("--max-torsions")) c.max_torsions = std::stoi(o["--max-torsions"]);
      if (o.count("--temp")) c.temperature = std::stod(o["--temp"]);
      if (o.count("--cutoff")) c.cutoff = std::stod(o["--cutoff"]);
      if (o.count("--seed")) c.seed = std::stoull(o["--seed"]);
      c.coulomb = !o.count("--no-coulomb");
      const ForceField ff = o.count("--ff") ? cli_forcefield(s, o) : default_forcefield(s);
      CbmcReport rep;
      cbmc_regrow(s, ff, c, &rep);
      const std::string out = o["-o"];
      auto ends = [&](const char* e) { return out.size() > 4 && out.substr(out.size() - 4) == e; };
      if (ends(".pdb")) write_pdb(s, out);
      else if (ends(".gro")) write_gro(s, out);
      else if (ends(".xyz")) write_xyz(s, out);
      else if (ends("mol2")) write_mol2(s, out);
      else if (ends(".car")) write_car(s, out);
      else {
        EnergyOptions eo;
        eo.coulomb = c.coulomb;
        if (const std::string why = write_lammps_data_or_structure(s, ff, eo, out); !why.empty()) std::fprintf(stderr, "%s: %s\n", out.c_str(), why.c_str());
      }
      for (const auto& n : rep.notes) std::printf("%s\n", n.c_str());
      std::printf("energy change %.3f kcal/mol · %.1f s\nwrote %s\n", rep.energy_change, rep.seconds, out.c_str());
      return 0;
    }
    if (cmd == "thermal") {
      if (!o.count("--dimensionality")) throw std::invalid_argument("--dimensionality 3d (bulk) or 2d (a sheet with vacuum along z): say which; CAPS does not guess it");
      ThermalOptions to;
      const std::string dim = o["--dimensionality"];
      if (dim != "2d" && dim != "3d" && dim != "2" && dim != "3") throw std::invalid_argument("--dimensionality 2d|3d");
      to.dimensionality = dim[0] == '2' ? 2 : 3;
      if (o.count("--thickness")) to.thickness = std::stod(o["--thickness"]);
      if (o.count("--method")) to.method = o["--method"];
      if (o.count("--temp")) to.temperature = std::stod(o["--temp"]);
      if (o.count("--eq-ps")) to.equilibrate_ps = std::stod(o["--eq-ps"]);
      if (o.count("--ps")) to.ps = std::stod(o["--ps"]);
      if (o.count("--dt")) to.dt = std::stod(o["--dt"]);
      if (o.count("--every")) to.sample_every = std::stoi(o["--every"]);
      if (o.count("--corr-ps")) to.corr_ps = std::stod(o["--corr-ps"]);
      if (o.count("--axis")) {
        const std::string a = o["--axis"];
        to.axis = a == "x" ? 0 : a == "y" ? 1 : a == "z" ? 2 : throw std::invalid_argument("--axis x|y|z");
      } else if (to.dimensionality == 2) to.axis = 0;
      if (o.count("--layers")) to.slabs = std::stoi(o["--layers"]);
      if (o.count("--swap-every")) to.swap_every = std::stoi(o["--swap-every"]);
      if (o.count("--transient")) to.transient = std::stod(o["--transient"]);
      if (o.count("--seed")) to.seed = std::stoull(o["--seed"]);
      if (o.count("--threads")) to.energy.threads = std::stoi(o["--threads"]);
      to.energy.coulomb = !o.count("--no-coulomb");
      electrostatics(to.energy, o);
      const bool quiet = o.count("--quiet");
      to.field = std::make_shared<ForceField>(cli_forcefield(s, o, quiet));
      run_cutoff(to.energy, *to.field, o, quiet);
      if (o.count("--lammps-input")) {   // the long run for the cluster: data file beside the input, κ computed by LAMMPS
        const std::string in = o["--lammps-input"];
        const std::string data = std::filesystem::path(in).replace_extension(".data").string();
        save_structure(s, *to.field, to.energy, data);
        LammpsRun run;
        run.kind = LammpsRun::Kind::None;
        write_lammps_input(s, *to.field, to.energy, std::filesystem::path(data).filename().string(), in, 0, false, run);
        ThermalLammps q;
        q.method = to.method, q.temperature = to.temperature, q.dt = to.dt, q.equilibrate_ps = std::max(to.equilibrate_ps, 100.0);
        q.ps = std::max(to.ps, 2000.0), q.sample_every = to.sample_every, q.corr_ps = to.corr_ps, q.axis = to.axis, q.slabs = to.slabs;
        q.swap_every = to.swap_every, q.dimensionality = to.dimensionality;
        q.thickness = to.dimensionality == 2 ? (to.thickness > 0 ? to.thickness : sheet_thickness(s)) : 0;
        q.kspace = to.energy.electrostatics == EnergyOptions::Electrostatics::PME;
        std::ofstream f(in, std::ios::app);
        f << thermal_lammps_run(q);
        std::printf("wrote %s and %s (%s, %.0f ps equilibration, %.0f ps sampled)\n", in.c_str(), data.c_str(), to.method.c_str(), q.equilibrate_ps, q.ps);
        if (!o.count("--ps")) return 0;   // only the inputs unless a local run is asked for too
      }
      if (std::llround(to.ps * 1000 / to.dt) + std::llround(to.equilibrate_ps * 1000 / to.dt) > 5000 && !o.count("--long"))
        std::fprintf(stderr, "note: %lld steps; long runs belong on a cluster (--lammps-input writes the inputs)\n",
                     std::llround((to.ps + to.equilibrate_ps) * 1000 / to.dt));
      to.progress = [&](const std::string& w, double f) {
        if (!quiet) std::fprintf(stderr, "\r%-45s %3.0f %%", w.c_str(), 100 * f);
        return true;
      };
      const ThermalResult r = thermal_conductivity(s, to);
      if (!quiet) std::fprintf(stderr, "\r%-60s\r", "");
      for (const auto& n : r.notes) std::printf("%s\n", n.c_str());
      std::printf("%s", thermal_summary(r).c_str());
      if (o.count("--json")) {
        Json j = Json::object();
        j["method"] = Json(r.method), j["dimensionality"] = Json(double(r.dimensionality)), j["thickness"] = Json(r.thickness);
        Json k = Json::array(), e = Json::array();
        for (int a = 0; a < 3; ++a) k.push_back(Json(r.kappa[size_t(a)])), e.push_back(Json(r.error[size_t(a)]));
        j["kappa"] = k, j["error"] = e, j["kappa_mean"] = Json(r.kappa_mean), j["kappa_mean_error"] = Json(r.kappa_mean_error);
        j["sheet_conductance"] = Json(r.sheet_conductance), j["temperature"] = Json(r.temperature), j["plateau"] = Json(r.plateau);
        j["flux"] = Json(r.flux), j["gradient"] = Json(r.gradient), j["exchanges"] = Json(double(r.exchanges));
        Json t = Json::array(), x = Json::array(), tt = Json::array();
        for (double v : r.t_ps) t.push_back(Json(v));
        for (double v : r.slab_x) x.push_back(Json(v));
        for (double v : r.slab_t) tt.push_back(Json(v));
        j["t_ps"] = t, j["slab_x"] = x, j["slab_t"] = tt;
        for (int a = 0; a < 3; ++a) {
          Json run = Json::array();
          for (double v : r.running[size_t(a)]) run.push_back(Json(v));
          j[std::string("running_") + char('x' + a)] = run;
        }
        std::ofstream(o["--json"]) << j.dump(1) << "\n";
        std::printf("wrote %s\n", o["--json"].c_str());
      }
      return 0;
    }
    if (cmd == "md") {
      if (!o.count("-o")) return usage();
      const auto md_clock = std::chrono::steady_clock::now();   // setup + run, for --bench-json
      DynamicsOptions d;
      if (o.count("--steps")) d.steps = std::stoll(o["--steps"]);
      if (o.count("--dt")) d.dt = std::stod(o["--dt"]);
      if (o.count("--temp")) d.temperature = std::stod(o["--temp"]);
      if (o.count("--thermostat")) d.thermostat = thermostat_from_string(o["--thermostat"]);
      if (o.count("--tau-t")) d.tau_t = std::stod(o["--tau-t"]);
      if (o.count("--barostat")) d.barostat = barostat_from_string(o["--barostat"]);
      if (o.count("--pressure")) d.pressure = std::stod(o["--pressure"]);
      if (o.count("--tau-p")) d.tau_p = std::stod(o["--tau-p"]);
      if (o.count("--constraints")) d.constraints = constraints_from_string(o["--constraints"]);
      if (o.count("--constraint-solver")) d.constraint_algorithm = constraint_algorithm_from_string(o["--constraint-solver"]);
      if (o.count("--hmr")) d.hmr = std::stod(o["--hmr"]);
      if (o.count("--seed")) d.seed = std::stoull(o["--seed"]);
      if (o.count("--thermo")) d.thermo_every = std::stoi(o["--thermo"]);
      if (o.count("--every")) d.frame_every = std::stoi(o["--every"]);
      if (o.count("--skin")) d.energy.skin = std::stod(o["--skin"]);
      if (o.count("--threads")) d.energy.threads = std::stoi(o["--threads"]);
      d.energy.coulomb = !o.count("--no-coulomb");
      electrostatics(d.energy, o);
      precision(d.energy, o);
      d.new_velocities = o.count("--new-velocities");
      d.field = std::make_shared<ForceField>(cli_forcefield(s, o, o.count("--quiet") > 0));
      run_cutoff(d.energy, *d.field, o, o.count("--quiet") > 0);
      if (o.count("--fix-mol")) {   // hold one molecule in place (the substrate of an interface is molecule 1)
        const int64_t mm = std::stoll(o["--fix-mol"]);
        d.fixed.assign(s.atoms.size(), 0);
        for (size_t i = 0; i < s.atoms.size(); ++i) d.fixed[i] = s.atoms[i].mol == mm;
      }
      const bool quiet = o.count("--quiet");
      // resume: the input is a checkpoint written by --checkpoint-every; --steps is the run's total length
      if (o.count("--resume")) {
        const auto at = s.title.find("checkpoint step ");
        if (at == std::string::npos) throw std::invalid_argument("--resume needs a checkpoint written by caps md --checkpoint-every (its title says the step)");
        const int64_t done = std::stoll(s.title.substr(at + 16));
        if (s.velocities.size() != s.atoms.size()) throw std::invalid_argument("the checkpoint carries no velocities");
        if (done >= d.steps) throw std::invalid_argument("the checkpoint is at step " + std::to_string(done) + ": --steps (the total) must be larger");
        d.step_offset = done;
        d.steps -= done;
        d.new_velocities = false;
        std::printf("resuming at step %lld: %lld steps to go\n", static_cast<long long>(done), static_cast<long long>(d.steps));
      }
      // checkpoints: the full state every N steps, written whole (to a temporary file, then renamed); also when a stop is
      // asked for (SIGUSR1 / SIGTERM: a cluster job near its time limit), so the run can continue with --resume
      install_stop_signals();
      const std::string cp = o.count("--checkpoint") ? o["--checkpoint"] : std::filesystem::path(o["-o"]).replace_extension(".restart.data").string();
      {
        if (o.count("--checkpoint-every")) d.checkpoint_every = std::stoll(o["--checkpoint-every"]);
        const int64_t total = d.step_offset + d.steps;
        d.checkpoint = [&, cp, total](const std::vector<double>& x, const std::vector<double>& v, const Cell& c, int64_t step) {
          System k = s;
          k.cell = c;
          k.velocities.resize(k.atoms.size());
          for (size_t i = 0; i < k.atoms.size(); ++i) {
            k.atoms[i].pos = {x[3 * i], x[3 * i + 1], x[3 * i + 2]};
            k.velocities[i] = {v[3 * i], v[3 * i + 1], v[3 * i + 2]};
          }
          k.unwrapped = true;
          k.title = "checkpoint step " + std::to_string(step) + " of " + std::to_string(total);
          const std::string tmp = cp + ".tmp";
          const std::string why = write_lammps_data_or_structure(k, d.field ? *d.field : default_forcefield(k), d.energy, tmp);
          if (!why.empty()) write_lammps_data(k, tmp);
          std::filesystem::rename(tmp, cp);
          if (!quiet) std::printf("checkpoint: step %lld → %s\n", static_cast<long long>(step), cp.c_str());
        };
      }
      // live output: the dump, the thermo CSV and the progress lines are written while the run goes (a cluster job is
      // followed by reading them; a run stopped by a time limit leaves what it had done)
      LiveOutput live("md");
      const bool resumed = o.count("--resume") > 0;   // the files of the run it continues grow on
      if (o.count("--dump")) {
        live.open_frames(o["--dump"], s, resumed);
        d.frame = [&](const std::vector<double>& x, const Cell& c, int64_t step) { live.frame(x, c, step); };
      } else {
        d.frame_every = 0;
      }
      if (o.count("--log")) live.open_thermo(o["--log"], resumed);
      if (const std::string pf = o.count("--progress-file") ? o["--progress-file"] : default_progress_file(); !pf.empty()) live.open_progress(pf, resumed);
      const double total_steps = double(d.step_offset + d.steps);
      if (!quiet) std::printf("%10s %9s %9s %13s %13s %13s %10s %8s\n", "step", "time/ps", "T/K", "Epot", "Etotal", "conserved", "P/atm", "ρ/g·cm⁻³");
      d.progress = [&](const ThermoRow& r) {
        if (!quiet) {
          std::printf("%10lld %9.3f %9.2f %13.3f %13.3f %13.3f %10.1f %8.4f\n", static_cast<long long>(r.step), r.time_ps, r.temperature, r.potential,
                      r.total, r.conserved, r.pressure, r.density);
          std::fflush(stdout);
        }
        live.thermo(r, "md");
        live.progress("md", 0, 0, total_steps > 0 ? double(r.step) / total_steps : -1, &r);
        return true;
      };
      DynamicsReport rep;
      try {
        run_dynamics(s, d, &rep);
      } catch (const DynamicsInterrupted& e) {
        // the checkpoint holds this step; the command that continues the run is written beside it
        std::string again = "caps md " + cp + " --resume --steps " + std::to_string(d.step_offset + d.steps) + " -o " + o["-o"];
        for (const char* k : {"--dt", "--temp", "--thermostat", "--tau-t", "--barostat", "--pressure", "--tau-p", "--constraints", "--constraint-solver", "--thermo",
                              "--every", "--cutoff", "--skin", "--threads", "--checkpoint-every", "--log", "--dump"})
          if (o.count(k)) again += std::string(" ") + k + " " + o[k];
        std::ofstream("resume.txt") << again << "\n";
        std::printf("stopped at step %lld (a stop was asked for): the state is in %s; continue with\n  %s\n", static_cast<long long>(e.step), cp.c_str(), again.c_str());
        live.done(false, "interrupted at step " + std::to_string(e.step) + "; resume: " + again);
        return 75;
      }
      std::filesystem::remove("resume.txt");   // finished: nothing to continue
      live.done(true);
      if (o.count("--bench-json")) {   // a scaling check's measurement (caps job scaling)
        const double wall = std::chrono::duration<double>(std::chrono::steady_clock::now() - md_clock).count();
        Json b = Json::object();
        b["version"] = std::string(version_string()), b["commit"] = std::string(commit_string());
        const char* node = std::getenv("SLURMD_NODENAME");
        if (!node) node = std::getenv("HOSTNAME");
        b["node"] = std::string(node ? node : "");
        b["threads"] = double(rep.threads), b["atoms"] = double(s.atoms.size()), b["steps"] = double(rep.steps);
        b["wall_s"] = wall, b["loop_s"] = rep.seconds, b["setup_s"] = std::max(0.0, wall - rep.seconds);
        b["steps_per_s"] = rep.seconds > 0 ? double(rep.steps) / rep.seconds : 0.0;
        b["ns_per_day"] = rep.ns_per_day;
        std::ofstream(o["--bench-json"]) << b.dump(1) << "\n";
      }
      const std::string out = o["-o"];
      auto ends = [&](const char* e) { return out.size() > 4 && out.substr(out.size() - 4) == e; };
      if (ends(".pdb")) write_pdb(s, out);
      else if (ends(".gro")) write_gro(s, out);   // same atoms and order as a GROMACS topology read with --topology
      else if (ends(".xyz")) write_xyz(s, out);
      else if (ends("mol2")) write_mol2(s, out);
      else if (ends(".car")) write_car(s, out);
      else if (const std::string why = write_lammps_data_or_structure(s, d.field ? *d.field : default_forcefield(s), d.energy, out); !why.empty()) std::fprintf(stderr, "%s: %s\n", out.c_str(), why.c_str());
      for (const auto& n : rep.notes) std::printf("%s\n", n.c_str());
      std::printf("wrote %s%s\n", out.c_str(), o.count("--dump") ? (" and " + o["--dump"]).c_str() : "");
      return 0;
    }
    if (cmd == "chains") {
      const InternalDistances d = internal_distances(s);
      std::printf("%d backbones · ⟨b²⟩ %.4f Å² · ⟨R²⟩ end to end %.1f Å²\n# n  <R2(n)>/(n b2)\n", d.chains, d.b2, d.r2_end);
      for (size_t k = 0; k < d.n.size(); ++k) std::printf("%d %.4f\n", d.n[k], d.ratio[k]);
      return 0;
    }
    if (cmd == "equilibrate") {
      install_stop_signals();
      ProtocolParams pp;
      if (o.count("--tfinal")) pp.t_final = pp.t_low = std::stod(o["--tfinal"]);
      if (o.count("--tmax")) pp.t_max = pp.t_high = std::stod(o["--tmax"]);
      if (o.count("--pfinal")) pp.p_final = std::stod(o["--pfinal"]);
      if (o.count("--pmax")) pp.p_max = std::stod(o["--pmax"]);
      if (o.count("--scale")) pp.time_scale = std::stod(o["--scale"]);
      if (o.count("--cycles")) pp.cycles = std::stoi(o["--cycles"]);
      if (o.count("--tlow")) pp.t_low = std::stod(o["--tlow"]);
      if (o.count("--thigh")) pp.t_high = std::stod(o["--thigh"]);
      if (o.count("--ramp")) pp.ramp_ps = std::stod(o["--ramp"]);
      if (o.count("--hold")) pp.hold_ps = std::stod(o["--hold"]);
      const std::string proto = o.count("--protocol") ? o["--protocol"] : "larsen21";
      EquilibrateOptions e;
      if (std::ifstream pf(proto); pf && proto.find('.') != std::string::npos) {
        std::stringstream ss;
        ss << pf.rdbuf();
        e.stages = parse_protocol(ss.str());
      } else {
        e.stages = protocol_by_name(proto, pp);
      }
      if (o.count("--print-protocol")) {
        std::printf("%s# total %.1f ps\n", protocol_text(e.stages).c_str(), protocol_ps(e.stages));
        return 0;
      }
      if (!o.count("-o")) return usage();
      if (o.count("--dt")) e.md.dt = std::stod(o["--dt"]);
      if (o.count("--thermostat")) e.md.thermostat = thermostat_from_string(o["--thermostat"]);
      if (o.count("--barostat")) e.md.barostat = barostat_from_string(o["--barostat"]);
      if (o.count("--tau-t")) e.md.tau_t = std::stod(o["--tau-t"]);   // as the Studio's Equilibrate (caps_equil_opts.tau_t / tau_p)
      if (o.count("--tau-p")) e.md.tau_p = std::stod(o["--tau-p"]);
      if (o.count("--seed")) e.md.seed = std::stoull(o["--seed"]);
      electrostatics(e.md.energy, o);
      precision(e.md.energy, o);
      if (o.count("--hmr")) e.md.hmr = std::stod(o["--hmr"]);
      if (o.count("--constraints")) e.md.constraints = constraints_from_string(o["--constraints"]);
      if (o.count("--threads")) e.md.energy.threads = std::stoi(o["--threads"]);
      // the force field: --ff (a library name, a file, or one written whole by CAPS), else the built-in default
      e.md.field = std::make_shared<ForceField>(cli_forcefield(s, o, o.count("--quiet") > 0));
      run_cutoff(e.md.energy, *e.md.field, o, o.count("--quiet") > 0);
      e.until_converged = o.count("--until-converged");
      // continued after a stop: --resume-at STAGE:STEP (written into resume.txt by the stopped run)
      if (o.count("--resume-at")) {
        const std::string ra = o["--resume-at"];
        const auto colon = ra.find(':');
        e.start_stage = std::stoi(ra.substr(0, colon)) - 1;
        e.start_step = colon == std::string::npos ? 0 : std::stoll(ra.substr(colon + 1));
        if (s.velocities.size() != s.atoms.size()) throw std::invalid_argument("--resume-at needs the checkpoint (it carries the velocities)");
      }
      const std::string eq_cp = o.count("--checkpoint") ? o["--checkpoint"] : std::filesystem::path(o["-o"]).replace_extension(".restart.data").string();
      e.checkpoint = [&](const std::vector<double>& x, const std::vector<double>& v, const Cell& c, int stage, int64_t at) {
        System k = s;
        k.cell = c;
        k.velocities.resize(k.atoms.size());
        for (size_t i = 0; i < k.atoms.size(); ++i) k.atoms[i].pos = {x[3 * i], x[3 * i + 1], x[3 * i + 2]}, k.velocities[i] = {v[3 * i], v[3 * i + 1], v[3 * i + 2]};
        k.unwrapped = true;
        k.title = "checkpoint equilibrate stage " + std::to_string(stage + 1) + " step " + std::to_string(at);
        const std::string tmp = eq_cp + ".tmp";
        const std::string why = write_lammps_data_or_structure(k, e.md.field ? *e.md.field : default_forcefield(k), e.md.energy, tmp);
        if (!why.empty()) write_lammps_data(k, tmp);
        std::filesystem::rename(tmp, eq_cp);
        std::ofstream("resume.txt") << resume_line("caps equilibrate", eq_cp, o, {{"--resume-at", std::to_string(stage + 1) + ":" + std::to_string(at)}}) << "\n";
      };
      if (o.count("--block")) e.block_ps = std::stod(o["--block"]);
      if (o.count("--max-blocks")) e.max_blocks = std::stoi(o["--max-blocks"]);
      if (o.count("--every-ps")) e.frame_ps = std::stod(o["--every-ps"]);
      const bool quiet = o.count("--quiet");
      LiveOutput live("equilibrate");
      const bool resumed = o.count("--resume-at") > 0;   // the files of the run it continues grow on
      if (o.count("--dump")) {
        live.open_frames(o["--dump"], s, resumed);
        e.frame = [&](const std::vector<double>& x, const Cell& c, int64_t step) { live.frame(x, c, step); };
      }
      if (o.count("--log")) live.open_thermo(o["--log"], resumed);
      if (const std::string pf = o.count("--progress-file") ? o["--progress-file"] : default_progress_file(); !pf.empty()) live.open_progress(pf, resumed);
      const double total_ps = protocol_ps(e.stages);
      std::printf("%s# total %.1f ps\n", protocol_text(e.stages).c_str(), total_ps);
      int last_stage = -1;
      double last_print = -1e9;
      e.progress = [&](int st, int n, const std::string& label, const ThermoRow& r) {
        live.thermo(r, label);
        live.progress(label, st, n, total_ps > 0 && !e.until_converged ? r.time_ps / total_ps : -1, &r);
        if (quiet) return true;
        if (st != last_stage || r.time_ps - last_print >= 5.0) {
          std::printf("[%2d/%d] %-28s %9.2f ps  T %7.1f K  P %9.0f atm  ρ %.4f  Epot %11.1f\n", st, n, label.c_str(), r.time_ps, r.temperature,
                      r.pressure, r.density, r.potential);
          std::fflush(stdout);
          last_stage = st;
          last_print = r.time_ps;
        }
        return true;
      };
      EquilibrateReport rep;
      try {
        equilibrate(s, e, &rep);
      } catch (const DynamicsInterrupted& ex) {
        std::ifstream rf("resume.txt");
        std::string again;
        std::getline(rf, again);
        std::printf("stopped at step %lld (a stop was asked for): the state is in %s; continue with\n  %s\n", static_cast<long long>(ex.step), eq_cp.c_str(), again.c_str());
        live.done(false, "interrupted; resume: " + again);
        return 75;
      }
      std::filesystem::remove("resume.txt");
      live.done(true);
      const std::string out = o["-o"];
      auto ends = [&](const char* x) { return out.size() > 4 && out.substr(out.size() - 4) == x; };
      if (ends(".pdb")) write_pdb(s, out);
      else if (ends(".gro")) write_gro(s, out);   // same atoms and order as a GROMACS topology read with --topology
      else if (ends(".xyz")) write_xyz(s, out);
      else if (const std::string why = write_lammps_data_or_structure(s, e.md.field ? *e.md.field : default_forcefield(s), e.md.energy, out); !why.empty()) std::fprintf(stderr, "%s: %s\n", out.c_str(), why.c_str());
      std::printf("\n%-30s %8s %8s %10s %9s\n", "stage", "ps", "T/K", "P/atm", "ρ/g·cm⁻³");
      for (const auto& st : rep.stages)
        std::printf("%-30s %8.1f %8.1f %10.0f %9.4f\n", st.label.c_str(), st.ps, st.temperature, st.pressure, st.density);
      for (const auto& c : rep.checks) {
        std::printf("check %-26s %s  change %.4g (tolerance %.4g)  blocks:", c.quantity.c_str(), c.ok ? "ok  " : "open", c.change, c.tolerance);
        for (double v : c.blocks) std::printf(" %.4f", v);
        std::printf("\n");
      }
      for (const auto& n : rep.notes) std::printf("%s\n", n.c_str());
      std::printf("wrote %s\n", out.c_str());
      return e.until_converged && !rep.converged ? 3 : 0;
    }
    if (cmd == "rdf") {
      const double rmax = o.count("--rmax") ? std::stod(o["--rmax"]) : 12.0, dr = o.count("--dr") ? std::stod(o["--dr"]) : 0.2;
      int ea = 0, eb = 0;
      if (o.count("--pair")) {
        const auto p = o["--pair"];
        const auto dash = p.find('-');
        ea = element_from_symbol(p.substr(0, dash));
        eb = element_from_symbol(dash == std::string::npos ? p : p.substr(dash + 1));
      }
      std::printf("# r(Å) g(r)\n");
      for (auto [r, g] : rdf(s, ea, eb, rmax, dr, o.count("--inter"))) std::printf("%.3f %.5f\n", r, g);
      return 0;
    }
    if (cmd == "edit") {
      std::string script = o.count("--ops") ? o["--ops"] : "";
      if (o.count("--ops-file")) {
        std::ifstream in(o["--ops-file"]);
        if (!in) throw std::runtime_error("cannot read " + o["--ops-file"]);
        for (std::string l; std::getline(in, l);) script += "\n" + l;   // comments: edit_script
      }
      if (script.find_first_not_of("; \t\r\n") == std::string::npos || !o.count("-o"))
        throw std::invalid_argument("caps edit FILE --ops 'OP ARGS; …' (or --ops-file) -o OUT");
      const auto report = edit_script(s, script);
      for (const auto& line : report) std::printf("%s\n", line.c_str());
      write_structure_file(s, o["-o"]);
      std::printf("%zu atoms · %zu bonds · wrote %s\n", s.atoms.size(), s.bonds.size(), o["-o"].c_str());
      return 0;
    }
    if (cmd == "convert") {
      if (pos.size() < 2) return usage();
      const std::string out = pos[1];
      write_structure_file(s, out);   // every writer by the extension (.car with its .mdf, POSCAR by name)
      std::printf("wrote %s\n", out.c_str());
      return 0;
    }
    return usage();
  } catch (const std::exception& e) {
    std::fprintf(stderr, "caps: %s\n", e.what());
    return 1;
  }
}
