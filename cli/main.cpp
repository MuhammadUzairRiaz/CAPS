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
#if !defined(_WIN32)
#include <sys/utsname.h>
#endif

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
#include "caps/exchange_protocols.hpp"
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
#include "caps/powder.hpp"
#include "caps/phonons.hpp"
#include "caps/crystalfit.hpp"
#include "caps/habit.hpp"
#include "caps/qspr.hpp"
#include "caps/refine.hpp"
#include "caps/dielectric.hpp"
#include "caps/powder_index.hpp"
#include "caps/sqs.hpp"
#include "caps/zmethod.hpp"
#include "caps/polystats.hpp"
#include "caps/restraints.hpp"
#include "caps/graft_chains.hpp"
#include "caps/adp.hpp"
#include "caps/domains.hpp"
#include "caps/ddft.hpp"
#include "caps/dendrimer.hpp"
#include "caps/miller.hpp"
#include "caps/cell_fit.hpp"
#include "caps/cell_relax.hpp"
#include "caps/table_tools.hpp"
#include "caps/pointgroup.hpp"
#include "caps/tension.hpp"
#include "caps/qha.hpp"
#include "caps/dpd.hpp"
#include "caps/qeq.hpp"
#include "caps/pores.hpp"
#include "caps/measured.hpp"
#include "caps/reaction_path.hpp"
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

// Whether a path exists, false when it cannot be looked at: the last root is the source folder of the build, which an
// installed copy on another computer may not be allowed to read (status() would throw "Permission denied")
bool exists_here(const std::string& p) {
  std::error_code ec;
  return std::filesystem::exists(p, ec);
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
        if (x == lower && exists_here(root + "/data/forcefields/" + e["file"].str())) return e["file"].str();
      }
    }
  } catch (const std::exception&) {}
  return "";
}

// A force field by path, or by its library name (pcff, compass, gaff2, oplsaa2024 …) from data/forcefields.
std::string ff_path(const std::string& s0) {
  if (s0.empty() || exists_here(s0)) return s0;
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
        if (exists_here(root + "/data/forcefields/" + name)) {
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

// A crystal by path, or from the library in data/crystals by its catalogue id (pe-orthorhombic, silicon), its file
// (polymers/pe-orthorhombic.cif) or that file's name (pe-orthorhombic.cif, periclase.cif), as caps phono3py and vasp-relax
// take it. An existing file always wins; a name found nowhere is refused naming both places looked in.
std::string crystal_path(const std::string& in) {
  if (in.empty() || exists_here(in)) return in;
  for (const auto& root : caps_roots()) {
    std::ifstream f(root + "/data/crystals/catalogue.json");
    if (root.empty() || !f) continue;
    std::stringstream ss;
    ss << f.rdbuf();
    try {
      const Json cat = Json::parse(ss.str());
      for (const auto& e : cat["crystals"].items()) {
        const std::string file = e.text("file");
        if (file.empty()) continue;
        if (e.text("id") == in || file == in || std::filesystem::path(file).filename().string() == in) {
          const std::string p = root + "/data/crystals/" + file;
          if (exists_here(p)) return p;
        }
      }
    } catch (const std::exception&) {}
  }
  throw std::invalid_argument(in + ": no such file, nor a crystal of data/crystals/catalogue.json");
}

// The library in data/forcefields whose "name" is this one (a CAPS data file records its force field by that name in its
// "# caps-ff name" line), or "" when none is. Only the head of each file is read: the name is one of its first keys.
std::string library_of_field_name(const std::string& name) {
  namespace fs = std::filesystem;
  if (name.empty()) return "";
  for (const auto& root : caps_roots()) {
    const fs::path dir = fs::path(root) / "data" / "forcefields";
    std::error_code ec;
    if (!fs::is_directory(dir, ec)) continue;
    std::vector<fs::path> files;
    for (const auto& e : fs::directory_iterator(dir, ec))
      if (e.path().extension() == ".json" && e.path().filename() != "catalogue.json") files.push_back(e.path());
    std::sort(files.begin(), files.end());
    for (const auto& f : files) {
      std::ifstream in(f);
      std::string head(4096, '\0');
      in.read(head.data(), std::streamsize(head.size()));
      head.resize(size_t(in.gcount()));
      const size_t k = head.find("\n \"name\": \"");
      if (k == std::string::npos) continue;
      std::string v;
      for (size_t i = k + 11; i < head.size() && head[i] != '"'; ++i) {
        if (head[i] == '\\' && i + 1 < head.size()) ++i;
        v += head[i];
      }
      if (v != name) continue;
      try {
        if (load_forcefield(f.string()).name == name) return f.string();
      } catch (const std::exception&) {
      }
    }
  }
  return "";
}

// Charges worth keeping: a file whose charge column is all zero (a builder's placeholder) has none.
bool file_charges(const System& s) {
  if (!s.has_charges) return false;
  for (const auto& a : s.atoms)
    if (std::fabs(a.charge) > 1e-9) return true;
  return false;
}

// --charges keep asked of a structure with no charges (a .pdb, .xyz or .gro has no charge column; a builder writes zeros):
// every atom would get 0 and the Coulomb energy 0, so say so rather than write it silently
void warn_keep_without_charges(const std::map<std::string, std::string>& o, const System& s) {
  const auto c = o.find("--charges");
  if (c != o.end() && c->second == "keep" && !file_charges(s))
    std::fprintf(stderr, "charges: --charges keep, but the structure has no charges (none in its format, or all zero); every atom gets "
                         "charge 0 and Coulomb is 0 (--charges auto gives the force field's own)\n");
}

const char* usage_text() {
  return "caps 0.1.0\n"
               "usage:\n"
               "  caps --version | --check-update   this build; whether a newer CAPS release is out (exit 0 no, 1 yes, 2 cannot tell)\n"
               "  caps report  [--email]                a problem report for support@caps-studio.org (version, system, the last failed command;\n"
               "                                        a wrong argument is not kept)\n"
               "  caps info    FILE [--topology DATA]\n"
               "  caps render  FILE -o OUT.png|OUT.svg [--topology DATA] [--frame N] [--size WxH]\n"
               "               [--bg dark|white|transparent] [--colour molecule|element|type|distance]\n"
               "               [--style ball|space|sticks|noh|backbone] [--yaw DEG] [--pitch DEG] [--zoom Z] [--no-cell] [--wrap]\n"
               "  caps shape   FILE [--topology DATA]          per-molecule Rg and shape\n"
               "  caps rdf     FILE [--topology DATA] [--rmax 12] [--dr 0.2] [--pair C-C] [--inter]   g(r) of a periodic cell (all images within\n"
               "               --rmax count, so it may pass half the cell); --pair takes element symbols or all\n"
               "  caps analyze FILE [--topology DATA] [--props density,rdf,sq,xray,neutron,rg,ree,cn,persistence,msd,diffusion,\n"
               "               relaxation,ced,delta,ffv,psd,crosslinks,entanglements] [--first N --last N --stride N] [--frame-ps X | --timestep-fs 1]\n"
               "               [--pair C-C --inter] [--qmax 25 --dq 0.02 --qdirect 4] [--probe 0] [--grid 0.4] [--ff FF.json] [--json OUT] [--csv DIR]\n"
               "               diffusion: [--fit 0.2:0.5] [--viscosity ETA (mPa·s) --temp T (Yeh–Hummer finite-size correction)]\n"
               "               interfaces (zprofile, adhesion, interaction, orientation): [--surface 1-3,7 (molecule ids; 1)] [--axis x|y|z (z)] [--zbin 0.5]\n"
               "               ir (infrared spectrum from the dipole): --temp T [--ir-max-ps 2] [--ir-max-wavenumber 4000]\n"
               "               adp (displacement parameters U, B; a profile of U_iso along --axis in --zbin bins)\n"
               "               pattern: --pattern '[#6][S:1][S:2][#6]' (2, 3 or 4 points by atom-map number; shared numbers: a centroid) [--pattern-bin W]\n"
               "               hbonds: [--hb-donors 'N O F'] [--hb-acceptors 'N O F S'] [--hb-criterion da|ha] [--hb-da 3.5] [--hb-ha Å]\n"
               "                [--hb-angle °] (da: ∠H–D···A at most, 30; ha: ∠D–H···A at least, 110) [--hb-between 'molecule 1 | molecule 2']\n"
               "  caps elastic FILE [--topology DATA] [--method strain|fluct|fluct-run] [--configs N] [--strain 1e-4] [--temp T] [--ps 100] [--ff FF.json] [--json OUT]\n"
               "  caps tensile DATA -o OUT.data [--axis x|y|z|xyz] [--rate 1e-3] [--strain 0.2] [--temp 300] [--dt 1] [--equilibrate 20] [--pressure 1] [--tau-p 1000] [--fixed-lateral]\n"
               "               [--fit 0.02] [--seed 1] [--new-velocities] [--dump T.lammpstrj [--every 1000]] [--checkpoint-every N] [--resume-curve OUT.tensile-curve.json]\n"
               "               [--precision mixed|double] [--ff FF.json [--overlay OV.json] [--charges MODE] [--complete ml]] [--csv DIR]\n"
               "               runs --equilibrate ps of NPT (each axis at --pressure atm; 20 ps is 20000 steps at 1 fs) before the pull, so the box starts\n"
               "               stress-free; --equilibrate 0 skips it for a cell already equilibrated (--eq-ps is the same option)\n"
               "  caps tg DATA -o OUT.data [--from 500 --to 200 --step 20 --ps 100] [--precision mixed|double] [--ff FF.json [--overlay OV.json] [--charges MODE] [--complete ml]] [--csv DIR]   |   caps tg --fit RUN [RUN ...] [--volume]\n"
               "  caps exchange arrhenius TABLE [TABLE2] [--reduced] [--from T --to T]   |   caps exchange modulus PTENSOR --cell-volume A3 --temp K --dt FS [--fit PS] [--g-inf MPa] [--corr-ps PS]\n"
               "               | caps exchange memory PARTNERS --frame-ps PS   |   caps exchange bridges START NOW [--axis z] [--plane 0.5] [--bin 1]\n"
               "               | caps exchange heal HEALED INTACT [--modulus-strain 0.03] [--at-strain E] [--up-to E]   [--json OUT] [--csv DIR]   (bond exchange and healing)\n"
               "               | caps exchange list   |   caps exchange show CHEMISTRY   (exchange chemistries and their barriers, data/exchange)\n"
               "               | caps exchange run DATA -o OUT [--chemistry ID] --kind swap|shift|bind [--sites SEL] [--partners SEL] [--labels chain|ends] [--capture 1.3]\n"
               "                 [--barrier E] [--binding E] [--valence N[,M]] [--relax 1000 [--relax-depth 3]] [--every 10] [--attempts 1] [--steps 10000] [--temp 300] [--dt FS (default: the force field's own step, 0.49 with --kg; else 1)] [--kg [--k-theta K] | --ff FF.json]\n"
               "               | caps exchange bond-react DATA -o STEM --chemistry ID [--temp 300] [--t-ref T] [--steps N] [--ff FF.json]   (LAMMPS fix bond/react)\n"
               "               | caps exchange bond-swap DATA -o STEM [--sites SEL] [--stickers | --both-ends] [--every 50] [--fraction 0.5] [--cutoff 1.3] [--seed 1] [--steps N] [--k-theta K]\n"
               "               | caps exchange cut|gap|crack|fragment DATA -o OUT [--axis z] [--plane 0.5] [--gap A] [--radius A [--centre x,y,z]] [--axes xyz]   (healing-test starts)\n"
               "               | caps exchange step-strain DATA --strain E --ps PS [--sample 10] | tv DATA --from T --to T --step T --ps PS [--pressure 1]\n"
               "               | caps exchange heal-run DATA --ps PS [--start START] [--axis z] [--plane 0.5] [--chunks 10] [--npt]\n"
               "                 (with the moves and field options of run, or --no-exchange; [-o OUT] [--json OUT] [--csv DIR])\n"
               "  caps frames  FILE OUT.lammpstrj|dcd|xyz|pdb|gro|trr [--topology DATA] [--first N] [--last N] [--stride N]   (read thinned)\n"
               "  caps convert FILE OUT.data|lmp|xyz|extxyz|pdb|mol2|car|gro|sdf|cif|vasp|POSCAR [--topology DATA] [--supercell auto|na,nb,nc]   (.car: with its .mdf)\n"
               "               a crystal cell too small for its bonds is refused for a format with bonds: --supercell auto writes the copies that hold them\n"
               "  caps table   CSV [--sort 'A,-B'] [--hist COLUMN [--rule fd|scott|sturges|N]] [--sdf OUT.sdf --structure COLUMN]\n"
               "               a table's rows ordered by columns (- descending), a histogram of a column, an SD file of the rows\n"
               "  caps pointgroup FILE [--tolerance 0.1] [--symmetrize -o OUT]   the molecule's point group and symmetry elements;\n"
               "               --symmetrize makes it exactly symmetric\n"
               "  caps qha     CRYSTAL [--mode volume|axes] [--lock c] [--range 0.02] [--points 7] [--temps 0,100,200,298.15,400] [--pressure 0]\n"
               "               [--no-relax-cell] [--ff FF.json] [--cutoff 10] [--supercell n,n,n] [--mesh m,m,m] [--csv OUT]\n"
               "               quasi-harmonic free energy: phonons at strained cells, the cell, α, B_T, C_v, C_p and S at each temperature\n"
               "               (--points: cells on the volume grid, at least 4, volume mode only; the axes mode uses a fixed grid of the\n"
               "               centre, 4 cells per free axis, 4 per pair of free axes and, with a, b and c all free, 4 along the diagonal: 29)\n"
               "  caps tension THERMO.csv|log.lammps [--axis z] [--interfaces 2] [--skip 0.2] [--blocks 5] [--cell Lz]\n"
               "               surface or interfacial tension of a slab cell from its pressure tensor (Kirkwood–Buff), mN/m\n"
               "  caps planes  FILE --hkl h,k,l [--shift 0] [--max N] [--json OUT]   lattice planes (hkl) of the cell: d-spacing,\n"
               "               normal, each plane cut to the cell, and the view angles that look down the normal\n"
               "               (--shift in units of d, 0 through the origin; planes on the cell's faces are not listed)\n"
               "  caps edit    FILE --ops 'OP ARGS; OP ARGS …' | --ops-file OPS.txt -o OUT   edits in order, atoms numbered from 1:\n"
               "               element SEL Sym · delete SEL · bond I J [order] · unbond I J · addh [SEL] · attach I SMILES ·\n"
               "               length I J Å · angle I J K ° · torsion I J K L ° · invert I · config I R|S · rotate SEL x,y,z ° ·\n"
               "               mirror SEL nx,ny,nz · move SEL dx,dy,dz · clean [SEL] · tacticity iso|syndio ·\n"
               "               rings: ring I J N [MEMBERS] [aromatic] · spiro I N [MEMBERS] · insert I J Sym ·\n"
               "               join I FILE [J|auto] [torsion°]   (MEMBERS: new ring atoms, e.g. O or N,C=O,C,C) ·\n"
               "               crystals: supercell na nb nc · primitive · niggli · conventional · redefine m11 … m33 · vacuum Å ·\n"
               "               nanowire u v w radius [repeats] [cylinder|hexagonal|square] [vacuum] ·\n"
               "               SEL: 3,5-9 · all · element:C,N · smarts:PATTERN · type:LABEL (numbers after a delete shift)\n"
               "  caps provenance FILE [--json | --bibtex | --methods] [--compare OTHER]   the steps that produced FILE (FILE.provenance.json)\n"
               "  caps bench   [T1 T2 … | --all] [--repeats 3] [--quick] [--out DIR] [--samples DIR] [--cells DIR] [--reference polymers.json]\n"
               "               the built-in validation suite; T6/T7 measure the equilibrated cells in --cells (else CAPS_BENCH_CELLS,\n"
               "               else ~/.caps/bench/cells), one per material named by its reference id, against data/reference/polymers.json\n"
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
               "  caps ppa MAP FRAMES [--method caps|z1|both|lammps] -o DIR entanglements: primitive paths, N_e estimators\n"
               "  caps mech decks|analyze …       tension decks and their analysis · caps cgdyn MAP DUMP -o DIR  chain dynamics\n"
               "  caps backmap REF_MAP REF_DATA --cg CG_MAP --frame CG_FRAME -o DIR   CG beads back to all atoms\n"
               "  caps bundle  FILE [--topology DATA] --steps S.json [-o OUT.caps-bundle.zip] [--include-input] [--frame N]\n"
               "                                   a figure with its data, pipeline, provenance and hashes (and the input)\n"
               "  caps reproduce BUNDLE.caps-bundle.zip   check every shipped file's sha256, rebuild the data from its input and pipeline, compare\n"
               "  caps run     RECIPE.yaml|json [--seed N] [--threads N] [--out DIR] [--json] [--log thermo.csv] [--dump traj.lammpstrj]   build → type → grow → relax → equilibrate →\n"
               "               md → analyze → export from one file (exit 0 ok · 2 input · 3 missing params · 4 failed run); an unknown\n"
               "               key in any stage is refused (exit 2); --dump: frames every frame_ps (or frame_every steps) of an md or\n"
               "               equilibrate stage, about 100 per stage that sets neither\n"
               "  caps run     PIPELINE.yaml|json [--input 'runs/*/X.lammpstrj'] [--topology DATA] [--frame first|last] [--csv OUT] [--out DIR: the outputs: block] [--branch NAME]\n"
               "                                   a saved pipeline over many inputs: one row of attributes per input; the topology is\n"
               "                                   --topology, else the pipeline's source: topology:, else STEM.data next to each dump\n"
               "  caps crystal --group 'P 42/m n m' --cell a,b,c[,α,β,γ] --sites 'Ti1 Ti 0 0 0; O1 O 0.3048 0.3048 0' -o OUT\n"
               "               (OUT.cif keeps the space group and its asymmetric unit, P 1 once caps symmetry clears it;\n"
               "               OUT.vasp or POSCAR for VASP, .data, .pdb, .xyz … hold every atom)\n"
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
               "               particle CRYSTAL.cif [--shape sphere|cube|octahedron|cuboctahedron|truncated-octahedron|icosahedron|rod|cone|frustum|tetrahedron|pyramid|hemisphere|fibre\n               --radius 12 --height 24 --top-ratio 0.5 --length 20 --passivate --whole-molecules --hkl 0,0,1]\n"
               "               [--units SMILES --chains 10 --dp 20 --density 0.9]   -o OUT   fillers, alone or in a polymer matrix\n"
               "  caps martini FILE -o CG.data [--itp CG.itp] [--centre geometric|mass]   Martini 3 small molecules mapped from\n"
               "               an all-atom structure (vermouth's mappings); water, ions and other molecules left out\n"
               "  caps pore    slit [--width 10 --layers 1 --lx 26 --ly 22 --vacuum] | cylinder CRYSTAL.cif [--width 14 --wall 6 --length 20 --passivate] |\n"
               "               framework CRYSTAL.cif [--supercell 2,2,2]   [--fluid SMILES --count N --tolerance 2 --seed 1] -o OUT   a fluid in a pore\n"
               "  caps pull    FILE [--normal | --axis x|y|z] [--distance 10] [--rate 5] [--spring 10] [--temp 300] [--surface 1] [--csv OUT]\n"
               "               pull-out / debonding of a film from a held surface: interfacial shear strength, work of separation\n"
               "  caps blend   --components SMILES1,SMILES2 [--weights 0.5,0.5] [--chains 8] [--dp 20] [--density 0.5] [--slabs | --droplet] -o OUT\n"
               "  caps semicrystalline -o OUT.data [--polymer pe|pvdf-beta|ipp-alpha|nylon6-alpha|pet|...|CELL.cif | --polymer list]\n"
               "               [--repeats 4,6] [--crystal 7.6 (nm) | --crystal-cells 30] [--amorphous 7.6 (nm)] [--loops 0.35 --bridges 0.1 --tails 0.55]\n"
               "               [--reentry adjacent|random] [--density 0.852] [--seed 1] [--ff auto|oplsaa2024|pcff|compass|oplsua2024]\n"
               "               [--no-relax: skip the CAPS relax, tens of minutes on one thread for a library cell of ~14,000 atoms]\n"
               "               [--report R.json] [--lammps-input IN.lmp] [--gromacs STEM]   a lamellar stack of any library crystal: stems, loops,\n"
               "               bridges and tails grown from the crystal's repeat unit (groups crystal, amorphous); defaults and force field per polymer\n"
               "  caps semicrystalline FILE --analyse [--report R.json] [--labels L.txt]   loops, bridges, tails, P2 crystallinity, profiles\n"
               "               of any frame; --labels: each atom's id, class (stem, loop, bridge, tail, none) and segment (side atoms as their backbone atom)\n"
               "  caps grow    -o OUT.data|OUT.pdb|OUT.xyz|OUT.mol2|OUT.gro|OUT.cif [--chains 10] [--dp 8] [--density 0.5 | --box 33]\n"
               "               copolymers: --sequence terminal|markov|penultimate --weights FEED [--r1 R --r2 R | --reactivity 'r11,r12;r21,r22'\n"
               "                | --transition 'p11,p12;p21,p22' | --penultimate r1,r1p,r2,r2p] [--drift none|chains|along --conversion 0.8]\n"
               "                [--sequence-report OUT.json]   composition, dyads, triads and run lengths, model against built\n"
               "               architecture: --architecture linear|star|comb|branched|dendrimer|graft [--arms 4] [--arm-dp 5] [--spacing 4]\n"
               "                [--branch-probability 0.1] [--generations 2]; graft: --graft-units 'SMILES,…' [--graft-dp 10] [--graft-on K]\n"
               "                [--graft-probability 0.1] [--graft-max N] [--graft-sequence KIND --graft-weights W]   (graft site: the atom marked [CH:3])\n"
               "               [--tacticity atactic|isotactic|syndiotactic] [--seed 1] [--trans] [--scale 1.0]\n"
               "               [--units '*CC(*)c1ccccc1,*CC(*)(C)C(=O)OC' --sequence homopolymer|alternating|block|random|shuffled|gradient|pattern\n"
               "                --weights 0.7,0.3 --blocks 20,20 --pattern AAB --ff FF.json --method trials|rosenbluth|rosenbluth_lj\n"
               "                --threads N --lookahead 1]   any repeat units (else polystyrene); several units and no --sequence:\n"
               "                alternating (random with --weights); gradient: A falls from the head to the tail, its share the first\n"
               "                weight, the other units in their weights' ratio; the contact scale is automatic\n"
               "                (lowered only where a chain is stuck) unless --scale fixes it\n"
               "  caps field   FILE [--topology DATA] [--ff FF.json [--overlay OV.json] [--charges MODE]] [--cutoff Å] [--forces OUT.txt]   the structure's force field (below): types, terms, energy (and per-atom forces)\n"
               "  caps relax   FILE -o OUT.data|OUT.pdb|OUT.xyz [--method lbfgs|cg|sd|fire] [--ftol 0.5] [--iterations 5000]\n"
               "               [--density 1.05] [--step 0.06] [--box-relax] [--pressure 1] [--no-pushoff] [--cutoff Å]\n"
               "               [--no-coulomb] [--ff FF.json [--overlay OV.json[,MORE.json]] [--charges auto|types|gasteiger|keep|qeq|increments] [--complete ml]] [--quiet]\n"
               "               (.data output carries the force field for LAMMPS; --quiet hides only the live progress line)\n"
               "               exit 3: OUT was written but the last minimisation stopped on --iterations before --ftol (not converged)\n"
               "               --restraints 'distance 1 2 3.0 20; angle 1 2 3 109.5; torsion 1 2 3 4 60:90 50; position 1-40 5 axes z radius 0.5'\n"
               "               or a file of them (relax, md, equilibrate; ff apply writes them into the LAMMPS and GROMACS files): ½ k Δ², flat between\n"
               "               LO:HI, atoms from 1, k kcal/mol/Å² or /rad² (defaults 10 and 50), 'ramp KEND' to change k over a run\n"
               "               without --ff (relax, md, equilibrate, field): the force field a data file CAPS wrote holds (its coefficients,\n"
               "               charges and 1-4 weights; refused when CAPS cannot rebuild it from the file), else GAFF (C and H) or UFF;\n"
               "               without --cutoff (relax, md, equilibrate): the force field's own cut-off, else 10 Å\n"
               "  caps cbmc    FILE -o OUT.data [--moves 1000] [--trials 8] [--max-torsions 4] [--temp 300] [--cutoff 9] [--seed 1]\n"
               "               [--ff FF] [--no-coulomb]   configurational-bias regrowth of chain ends (Siepmann & Frenkel)\n"
               "  caps fields --molecules 'A8B8:1; A8:0.5' [--chi AB=0.9] [--grid 32,32,32 (nz 1: 2D)] [--bond 1.1543] [--kappa 10] [--noise 75] [--dt 0.02]\n"
               "               [--steps 2000] [--every 100] [--seed 1] [--vtk OUT.vtk] [--csv OUT.csv]   density-field dynamics of Gaussian chains (blends, block copolymers)\n"
               "  caps cellrelax FILE -o OUT [--hold c,beta] [--pressure 0] [--ftol 0.01] [--steps 3000] [--cutoff 10] [--ff FF.json]\n"
               "               the cell and atoms relaxed with chosen lattice parameters held (a b c alpha beta gamma)\n"
               "  caps path A B -o PATH.xyz [--no-match] [--method idpp|lst|linear] [--images 7] [--neb [--ff FF.json [--overlay OV.json] [--charges auto|types|gasteiger|keep|qeq|increments]]\n"
               "           [--cutoff 10] [--no-coulomb] [--relax-ends] [--spring 5] [--fmax 0.05] [--no-climb] [--steps 2000]]\n"
               "           atoms of B matched to A, the images between them, and the climbing-image NEB barrier with the force field\n"
               "  caps data FILE [-o OUT.csv|.npy (a curve) | OUT.cube|.vtk|.npy (a .cube grid)] [--json OUT]   a measured curve (columns, .xrdml, GSAS, .uxd, Rigaku, JCAMP-DX) or a .cube grid\n"
               "  caps symmetry FILE [--find [--tolerance 0.1]] [--clear] [--cell a,b,c,alpha,beta,gamma] [-o OUT.cif]\n"
               "             the space group kept with the structure: found or read, the cell fitted to its crystal system, the CIF as its asymmetric unit\n"
               "  caps pores FILE [--probe 1.4] [--grid 0.25] [--samples 2000] [--radii bondi|uff] [--json OUT] [--cube FIELD.cube]\n"
               "             channels and pockets, LCD, PLD per axis, accessible volume and surface (per cell, m²/cm³, m²/g)\n"
               "  caps qeq FILE [--kernel ohno|slater|point] [--ewald|--no-ewald] [--accuracy 1e-7] [--lambda 0.4913] [--hydrogen] [--bounds 'O:-1.2:0,Si:0:2.4']\n"
               "           [--cutoff 10] [--total 0] [--chg OUT.chg] [-o OUT.data|mol2|car]   charge equilibration (Slater shielding; Ewald by default for cells\n"
               "           periodic along a, b and c unless the cell is too large for it to be quick; --ewald always, --no-ewald never:\n"
               "           the tapered kernel summed over the images within the cut-off); charges that did not converge exit 1 and are not\n"
               "           written unless [--allow-unconverged]\n"
               "  caps dendrimer -o OUT [--preset pamam|ppi|bismpa --generation G] [--core SMILES] [--unit SMILES] [--end SMILES]\n"
               "               [--focal SMILES] [--shells 2] [--growth perfect|sma|flory --units 30 --seed 1] [--no-relax]\n"
               "               AB_n units with sites [*:1] (A) and [*:2] (B); prints units per shell and the degree of branching\n"
               "  caps dpd --molecules 'AAAABBBB:300@shell:in; C:9000' [--regions 'core droplet r=3 prio=2; shell shell r=5 in=3 prio=1'] [--box 16]\n"
               "           [--chi AB=25] [--density 3] [--steps 20000] [--equilibration STEPS/4] [--every STEPS/40] [--seed 1] [--spec SPEC.json] [-o OUT.xyz] [--json OUT]\n"
               "           dissipative particle dynamics; a sequence is bead letters A-Z with optional counts (A4B4 = AAAABBBB), branches in parentheses;\n"
               "           regions (droplet, rod, slab, shell, tube; c=fx,fy,fz r= in= t= axis= prio=) compose the start;\n"
               "           frames write bead types as elements A=C, B=O, C=N, D=S, others P (select them with 'element O' in domains, analyze)\n"
               "  caps index DATA.xy [--wavelength 1.54056] [--tolerance 0.03] [--systems cubic,tetragonal,hexagonal,orthorhombic,monoclinic] [--max-length 20]\n"
               "             [--lines 20] [--unindexed 0] [--max-boxes 20000000] [--json OUT]   cells that index the pattern's lines (successive\n"
               "             dichotomy), by M20; a cell that indexes every line with M20 100 or more ends the search of lower symmetries\n"
               "  caps solve DATA.xy MOLECULE --cell a,b,c,al,be,ga --space-group 'P 21/c' [--wavelength Å] [--steps 6000] [--restarts 2] [--seed 1]\n"
               "             [--uvw U,V,W] [--xy X,Y] [--zero DEG] -o CRYSTAL.cif   a molecular crystal found by simulated annealing against the\n"
               "             pattern (direct space); --steps counts all moves, shared by the restarts; the peak profile is fitted to the data\n"
               "             (Pawley, in the cell) first unless --uvw, --xy or --zero give it\n"
               "  caps zmelt CRYSTAL --temps 100,200,300,400 [--ff FF] [--steps 10000] [--dt 1] [--bisect 3] [--cutoff 10] [--sample-every 20]\n"
               "              [--csv OUT] [--json OUT]   melting point by the Z method: NVE runs from the perfect crystal (Belonoshko et al.)\n"
               "  caps domains FILE [--group QUERY] [--grid 0.8] [--sigma 1.5] [--threshold ρ | the mean] [--below] [--min-voxels 1]\n"
               "              [--csv OUT] [--cube OUT.cube] [--slice 'px,py,pz; ux,uy,uz; vx,vy,vz; w,h; nx,ny' [--slice-csv OUT]]\n"
               "              domains of a group's density field: volumes, radii, ellipsoid axes, percolation, spacing; slices on any plane\n"
               "  caps adp TRAJ [--topology DATA] [--first N] [--last N] [--stride N] [--supercell 2,2,2] [--probability 0.5] [--no-drift]\n"
               "              [--csv OUT] [--cif OUT.cif]   atomic displacement tensors U_ij, U_iso and B per atom (per site of a supercell)\n"
               "  caps cif FILE.cif [--disorder major|group|random] [--group LABEL] [--supercell 4,4,4] [--seed 1] [-o OUT] [--json OUT]\n"
               "              disorder groups, mixed and split sites, vacancies: the major group, a named group, or a supercell filled at random\n"
               "  caps graft FILE --units 'SMILES,…' --sites 12,40 [--equivalent] [--probability 1] [--per-site 1] [--max N] [--min-spacing Å]\n"
               "              [--dp 10] [--sequence KIND --weights W] [--padding Å] [--seed 1] -o OUT   chains grown from chosen atoms (from 1)\n"
               "  caps copolymer --feed 0.3,0.3,0.4 (--r1 R --r2 R | --reactivity 'r11,r12,r13;…' | --transition '…' | --penultimate r1,r1p,r2,r2p)\n"
               "               [--conversion 0.9 [--points 21]] [--csv OUT.csv] [--json OUT]   copolymer model: composition, dyads, triads, runs, drift\n"
               "  caps sqs CRYSTAL --from Cu --to Au [--fraction 0.5] [--supercell 2,2,2] [--shells 3] [--steps 20000] [--restarts 4] [--seed 1] -o OUT.cif\n"
               "              special quasirandom structure: the random alloy's pair correlations on the first shells (Zunger et al.)\n"
               "  caps enumerate CRYSTAL --from Cu --to Au [--count 1] [--supercell 2,2,2] [--max-configurations 200000] [--write DIR]\n"
               "              the symmetry-distinct occupations of a sublattice by COUNT substituents, with orbit sizes\n"
               "  caps dielectric CRYSTAL [--ff FF.json] [--cutoff 10] [--no-relax] [--fixed-cell] [--no-piezo] [--damping 10] [--csv STEM] [--json OUT]\n"
               "              dielectric tensors (static, high-frequency), refractive indices, piezoelectric e and d, reflectivity\n"
               "  caps refine DATA.xy CRYSTAL [CRYSTAL …] [--pawley] [--anode Cu | --wavelength Å] [--neutron] [--background 6] [--two-theta 10,90]\n"
               "              [--refine U,W,Y,zero] [--b-iso] [--preferred h,k,l] [--iterations 50] [--csv OUT] [--json OUT]   Rietveld or Pawley fit of a measured pattern\n"
               "  caps qspr descriptors [FILE ...] [--smiles 'CCO;CC(C)O'] [--csv OUT]   molecular descriptors (graph, shape, charged surface areas)\n"
               "  caps qspr model TABLE.csv --y COL [--x A,B,…] [--method mlr|pls|gfa|nn] [--stepwise --squares --products] [--components N] [--hidden 4]\n"
               "                  [--folds K] [--randomise N] [--seed 1] [--json OUT]   a model of y with cross-validation, y-randomisation and outliers\n"
               "  caps qspr cluster TABLE.csv [--x …] [--method ward|average|complete|single|kmeans] [--k 3] [--no-scale]   clusters of the rows\n"
               "  caps qspr lda TABLE.csv --class COL [--x …]   Fisher's discriminant analysis of the rows' classes\n"
               "  caps habit CRYSTAL [--ff FF.json] [--max-index 3] [--min-d Å] [--attachment] [--surface-energy] [--forms 12] [--particle RADIUS_Å --model bfdh|growth|equilibrium --out P.pdb]\n"
               "             [--shape-obj SHAPE.obj] [--json OUT]   crystal habit: BFDH faces, attachment and surface energies, Wulff shapes, a particle in that shape\n"
               "  caps fit SPEC.json [--ff FF.json] [--relaxed | --conventional] [--iterations 50] [--out FITTED.json] [--json REPORT.json]\n"
               "               force-field parameters fitted to crystal structures, lattice energies and elastic constants (the spec lists them)\n"
               "  caps phonons CRYSTAL [--ff FF.json] [--cutoff Å] [--supercell n1,n2,n3] [--mesh m1,m2,m3] [--path 'G 0 0 0; X 0.5 0 0; …'] [--displacement 0.01]\n"
               "               [--sigma 5 (cm⁻¹)] [--temps 100,300] [--csv STEM] [--json OUT]   force-field phonons: dispersion, DOS, harmonic thermodynamics, sound speeds\n"
               "               (CRYSTAL for phonons, qha and dielectric: a file, or a library crystal by id or file name: silicon, pe-orthorhombic.cif)\n"
               "  caps powder FILE [--topology DATA] [--anode Cu|Mo|Co|Fe|Cr|Ag | --wavelength Å] [--no-kalpha2] [--neutron] [--two-theta 5:90] [--step 0.02]\n"
               "               [--profile tch|pseudo-voigt|gaussian|lorentzian] [--uvw U,V,W] [--xy X,Y] [--size Å | --size-abc La,Lb,Lc] [--strain ε] [--biso B]\n"
               "               [--preferred h,k,l --march r] [--zero °] [--monochromator 2θ°] [--polarised] [--first N --last N --stride N] [--peaks 20] [--csv OUT] [--json OUT]\n"
               "               the powder pattern of the cell (its reciprocal lattice, frames averaged): reflections with hkl, d, 2θ, multiplicity, |F|, I\n"
               "  caps thermal FILE --dimensionality 3d|2d [--method gk|nemd] [--temp 300] [--eq-ps 10] [--ps 100] [--dt 1] [--every 5]\n"
               "               [--corr-ps 5] [--axis x|y|z] [--layers 20] [--swap-every 100] [--transient 0.3] [--thickness 3.35] [--seed 1]\n"
               "               [--ff FF] [--json OUT.json] [--lammps-input IN.lmp [--lammps-eq-ps 100] [--lammps-ps 2000]]\n"
               "               thermal conductivity by MD: Green-Kubo (heat-flux autocorrelation, centroid form for bonded terms)\n"
               "               or Muller-Plathe NEMD; 2d: a sheet with vacuum along z,\n"
               "               in-plane only, kappa over an effective thickness and the sheet conductance; NEMD takes any force field\n"
               "               and an even --layers of at least 10; Green-Kubo here takes harmonic bonds and angles, Fourier torsions\n"
               "               and impropers and 12-6 LJ with DSF Coulomb (for other terms use nemd or --lammps-input); --eq-ps/--ps are the\n"
               "               local run, the LAMMPS input's lengths are --lammps-eq-ps/--lammps-ps (with --lammps-input, --ps also runs here)\n"
               "  caps md      FILE -o OUT.data [--steps 10000] [--dt 1] [--temp 300] [--thermostat bussi|langevin|nose-hoover|none]\n"
               "               [--tau-t 100] [--barostat none|crescale|berendsen|mtk] [--pressure 1] [--tau-p 1000]\n"
               "               [--constraints none|h-bonds|all-bonds] [--constraint-solver shake|lincs] [--hmr 3.024] [--seed 1] [--new-velocities] [--thermo 100] [--dump TRAJ.lammpstrj --every 1000]\n"
               "               [--log thermo.csv] [--cutoff Å] [--skin 1.5] [--threads N] [--no-coulomb] [--precision mixed|double] [--ff FF.json [--overlay OV.json] [--charges MODE] [--complete ml]] [--quiet]\n"
               "               [--checkpoint-every N [--checkpoint OUT.restart.data]]   the full state every N steps\n"
               "               [--progress-file progress.jsonl]   one JSON line per thermo row (default in a SLURM job)\n"
               "               caps md OUT.restart.data --resume --steps TOTAL -o OUT   continues from a checkpoint with the first run's settings\n"
               "               (kept beside it in OUT.restart.settings.json; --dt, --constraints and --hmr cannot change on a resume)\n"
               "  caps equilibrate FILE -o OUT.data [--protocol larsen21|annealing|pushoff|PROTOCOL.txt] [--print-protocol]\n"
               "               [--tfinal 300] [--tmax 600] [--pfinal 1] [--pmax 49346] [--scale 1] (atm, K; --scale shortens every stage)\n"
               "               [--cycles 3] [--tlow 300] [--thigh 600] [--ramp 50] [--hold 50]   (annealing, ps)\n"
               "               [--until-converged] [--block 20] [--max-blocks 20] [--dt 1] [--thermostat bussi|langevin|nose-hoover] [--tau-t 100]\n"
               "               (--until-converged exits 3 when --max-blocks run out before every check passes; OUT is still written)\n"
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
               "               [--insert SMILES --count N] [--remove-unreacted] [--remove-pendant] [--summary NET.json] [--fa 2 --fb 4 --ratio 1] [--list-templates] [--quiet]\n"
               "  caps ml eval FILE [--model M.model] [--json OUT.json]   energy (eV) and forces (eV/A) of a molecule from CAPS's\n"
               "               machine-learned potential (data/ml/caps-mace-v1.model; a build with -DCAPS_WITH_TORCH=ON)\n"
               "  caps ml relax JOBS.json -o OUT.json [--model M.model]   constrained relaxations (bonds / angles / dihedrals held)\n"
               "  caps ff complete FILE --ff FF.json [--overlay OV.json] -o OVERLAY.json [--report R.json] [--model M.model] [--step 15]\n"
               "               the bonded terms FF lacks for FILE, found and fitted with CAPS's machine-learned potential\n"
               "  caps ff import-lt FILE.lt -o FF.json          convert a moltemplate force field into the CAPS format\n"
               "  caps ff info FF.json|NAME|uff                  types, rules, styles, references\n"
               "  caps ff overlays [FF_ID]                       the library's overlays (parameter refits) for --overlay\n"
               "  caps ff type FILE --ff FF.json|NAME|uff [--typing RULES.json] [-o TYPES.txt] [--explain]   assign atom types from SMARTS rules (UFF: its own typer)\n"
               "  caps ff apply FILE --ff FF.json [-o OUT.data [--lammps-input OUT.in [--lammps-run check|minimize|nvt|npt --temp 300 --press 1 --dt FS (default: the force field's, Martini 20, else 0.5) --steps N] [--moltemplate SYSTEM.lt]]] [--gromacs STEM] [--dlpoly DIR] [--overlay USER.json[,MORE.json]] [--types TYPES.txt] [--charges auto|keep|types|gasteiger]\n"
               "               [--complete ml [--model M.model]]   (terms the force field lacks fitted by CAPS ML, marked as such in every file)\n"
               "               [--lammps-style native|exact] [--hybrid] [--kspace auto|pppm|ewald|dsf|cut] [--kspace-accuracy 1e-4] [--lammps-cutoff Å] [--units auto|real|metal]\n"
               "               [--supercell auto|NA,NB,NC]   copies of the cell first (a --types file names the cell's atoms); a cell too small for its bonds is refused without them (ff type types it on its copies)\n"
               "               [--list] [-o OUT.data]   parameters for a structure whose atoms carry type names (or TYPES.txt)\n";
}

std::string usage_of(const std::string& cmd);
std::string g_usage_cmd;   // the command main() is running: usage() shows its entries, not all 250 lines

// the usage on stderr, exit 2: a command line CAPS could not run (the running command's entries when it has some)
int usage() {
  const std::string u = g_usage_cmd.empty() ? std::string() : usage_of(g_usage_cmd);
  if (u.empty()) std::cerr << usage_text();
  else std::cerr << "caps " << g_usage_cmd << ": a required input or option is missing\nusage:\n" << u << "(caps --help lists every command)\n";
  return 2;
}

// One command's entries of the usage text (a line '  caps CMD …' and its indented continuation lines, cut as the
// website's command reference cuts them); empty when the usage has no entry for it.
std::string usage_of(const std::string& cmd) {
  std::string out, line;
  bool in = false;
  std::istringstream u(usage_text());
  while (std::getline(u, line)) {
    if (line.rfind("  caps ", 0) == 0 && line.size() > 7 && line[7] != '-') {
      size_t e = 7;
      while (e < line.size() && line[e] != ' ') ++e;
      in = line.compare(7, e - 7, cmd) == 0 && e - 7 == cmd.size();
    } else if (line.rfind("    ", 0) != 0) {
      in = false;
    }
    if (in) out += line.substr(2) + "\n";
  }
  return out;
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
    "--during-md", "--keep-byproducts", "--relax-iterations", "--remove-unreacted", "--remove-pendant", "--sites-per-chain", "--summary", "--axis", "--barostat", "--beads", "--bg", "--bibtex",
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
    "--json", "--kspace", "--kspace-accuracy", "--lammps-cutoff", "--lammps-ps", "--lammps-eq-ps", "--units", "--lammps-input", "--lammps-run", "--moltemplate",
    "--lammps-style", "--last", "--layers", "--length", "--list", "--list-templates", "--log", "--lx", "--ly", "--m",
    "--labels", "--email", "--lookahead", "--martini", "--max-blocks", "--max-boxes", "--max-torsions", "--moves", "--max-strain", "--md-ps", "--method", "--methods", "--model", "--molecule-size",
    "--molecules", "--n", "--names", "--neutralise", "--new-velocities", "--no-cell", "--wrap",
    "--no-coulomb", "--no-ions", "--no-orthogonal", "--no-pbc", "--no-pushoff", "--no-relax",
    "--no-tail", "--normal", "--out", "--overlay", "--padding", "--pair", "--particles",
    "--ops", "--ops-file", "--checkpoint-every", "--checkpoint", "--resume", "--resume-at", "--bench-json", "--resume-points", "--resume-curve", "--height", "--top-ratio", "--passivate", "--pattern", "--per-cycle", "--perspective", "--pfinal", "--pitch", "--pmax", "--pme",
    "--pme-order", "--pme-spacing", "--precision", "--press", "--pressure", "--primitive", "--print-protocol", "--probe", "--progress-file",
    "--props", "--protocol", "--ps", "--qdirect", "--qmax", "--quick", "--quiet", "--radius", "--ramp", "--rate",
    "--ratio", "--batch", "--repeats", "--cells", "--reference", "--report", "--rmax", "--salt", "--samples", "--scale", "--seed", "--sequence",
    "--shape", "--sites", "--size", "--skin", "--slabs", "--solvent", "--solvents", "--spring", "--step",
    "--steps", "--strain", "--stride", "--style", "--supercell", "--surface",
    "--symmetrize", "--table", "--tacticity", "--target", "--tau-p", "--tau-t", "--temp", "--template",
    "--termination", "--tfinal", "--thermo", "--thermostat", "--thigh", "--threads", "--timestep-fs", "--tlow",
    "--tmax", "--to", "--tol", "--tolerance", "--topology", "--trans", "--trials", "--types", "--typing", "--units",
    "--polymer", "--crystal", "--crystal-cells", "--amorphous", "--loops", "--bridges", "--tails", "--reentry", "--analyse",
    "--corr-ps", "--long", "--swap-every", "--transient", "--dimensionality", "--thickness",
    "--until-converged", "--vacuum", "--viscosity", "--volume", "--cell-volume", "--reduced", "--kind", "--chemistry", "--t-ref", "--barrier", "--binding", "--attempts", "--partners", "--labels", "--kg", "--k-theta", "--stickers", "--both-ends", "--fraction", "--shells", "--restraints", "--restarts", "--max-configurations", "--write", "--bisect", "--sample-every", "--r1", "--r2", "--reactivity", "--transition", "--penultimate", "--drift", "--conversion", "--sequence-report", "--feed", "--points", "--architecture", "--arms", "--arm-dp", "--spacing", "--branch-probability", "--generations", "--graft-units", "--graft-sequence", "--graft-weights", "--graft-dp", "--graft-on", "--graft-probability", "--graft-max", "--equivalent", "--per-site", "--disorder", "--probability", "--no-drift", "--cif", "--hb-donors", "--hb-acceptors", "--hb-criterion", "--hb-da", "--hb-angle", "--hb-ha", "--hb-between", "--pattern-bin", "--threshold", "--below", "--min-voxels", "--slice", "--slice-csv", "--cube", "--whole-molecules", "--min-spacing", "--max", "--probability", "--axes", "--sample", "--no-exchange", "--start", "--npt", "--anode", "--wavelength", "--no-kalpha2", "--neutron", "--two-theta", "--profile", "--uvw", "--xy", "--eta", "--size-abc", "--scherrer-k", "--biso", "--preferred", "--march", "--zero", "--monochromator", "--polarised", "--peaks", "--mesh", "--path", "--displacement", "--sigma", "--temps", "--tolerance", "--systems", "--max-length", "--unindexed", "--lines", "--space-group", "--cell", "--fixed-cell", "--no-piezo", "--damping", "--pawley", "--b-iso", "--background", "--refine", "--smiles", "--y", "--x", "--stepwise", "--squares", "--products", "--folds", "--randomise", "--k", "--class", "--no-scale", "--hidden", "--ir-max-ps", "--ir-max-wavenumber", "--relaxed", "--conventional", "--max-index", "--forms", "--attachment", "--surface-energy", "--particle", "--shape-obj", "--min-d", "--chi", "--molecules", "--images", "--neb", "--fmax", "--no-climb", "--relax-ends", "--no-match", "--find", "--clear", "--samples", "--radii", "--kernel", "--lambda", "--total", "--bounds", "--hydrogen", "--ewald", "--no-ewald", "--accuracy", "--chg", "--allow-unconverged", "--spec", "--regions", "--gamma", "--equilibration", "--rc-angstrom", "--grid", "--bond", "--kappa", "--noise", "--vtk", "--chunks", "--valence", "--relax", "--relax-depth", "--mode", "--range", "--lock", "--no-relax-cell", "--slab", "--interfaces", "--skip", "--blocks", "--cell", "--shift", "--hold", "--hist", "--rule", "--sort", "--sdf", "--structure", "--preset", "--unit", "--core", "--end", "--focal", "--generation", "--growth", "--g-inf", "--plane", "--bin", "--modulus-strain", "--at-strain", "--up-to", "--wall", "--weights", "--width", "--yaw", "--zbin", "--zoom"};
  return k;
}

// The candidate closest to a mistyped word (edit distance), for an error message. Only a plausible typo is offered:
// at most a third of the word's letters (leading dashes not counted, at least one) changed, and of two equally near
// candidates the one of the same length; none when that still leaves a tie.
template <class Words>
std::string closest_word(const std::string& a, const Words& words) {
  size_t letters = a.size();
  for (size_t k = 0; k < a.size() && a[k] == '-'; ++k) --letters;
  const size_t limit = std::max<size_t>(1, letters / 3);
  std::string best;
  size_t bd = 99;
  bool tie = false;
  for (const std::string& k : words) {
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
    const size_t dk = d[k.size()];
    if (dk < bd) bd = dk, best = k, tie = false;
    else if (dk == bd) {
      const bool same_k = k.size() == a.size(), same_best = best.size() == a.size();
      if (same_k && !same_best) best = k, tie = false;
      else if (same_k == same_best) tie = true;
    }
  }
  return bd <= limit && !tie ? best : std::string();
}

// The known option closest to a mistyped one, for the error message.
std::string closest_option(const std::string& a) { return closest_word(a, known_options()); }

// Options without a value (the rest take the next word)
bool is_cli_switch(const std::string& a) {
  return a == "--no-cell" || a == "--wrap" || a == "--inter" || a == "--perspective" || a == "--trans" || a == "--escalate" ||
         a == "--box-relax" || a == "--no-pushoff" || a == "--no-coulomb" || a == "--quiet" || a == "--new-velocities" ||
         a == "--until-converged" || a == "--print-protocol" || a == "--no-pbc" ||
         a == "--no-relax" || a == "--no-relax-cell" || a == "--list-templates" || a == "--list" || a == "--allow-missing" || a == "--no-tail" || a == "--explain" || a == "--names" || a == "--fixed-lateral" || a == "--volume" || a == "--quick" || a == "--all" || a == "--pme" || a == "--no-orthogonal" || a == "--passivate" || a == "--auto-scale" || a == "--finite" || a == "--flake" || a == "--normal" || a == "--slabs" || a == "--droplet" || a == "--include-input" || a == "--primitive" || a == "--symmetrize" || a == "--find-symmetry" || a == "--groups" || a == "--neutralise" || a == "--no-ions" || a == "--solvents" || a == "--bibtex" || a == "--json" || a == "--deterministic" || a == "--vacuum" || a == "--methods" ||
         a == "--hybrid" || a == "--resume" || a == "--between-chains" || a == "--keep-byproducts" || a == "--auto-capture" ||
         a == "--during-md" || a == "--remove-unreacted" || a == "--remove-pendant" || a == "--analyse" || a == "--long" || a == "--reduced" || a == "--kg" || a == "--stickers" || a == "--both-ends" || a == "--no-exchange" || a == "--npt" ||
         a == "--no-kalpha2" || a == "--neutron" || a == "--polarised" || a == "--relaxed" || a == "--conventional" || a == "--attachment" || a == "--surface-energy" || a == "--stepwise" || a == "--squares" || a == "--products" || a == "--no-scale" || a == "--pawley" || a == "--b-iso" || a == "--fixed-cell" || a == "--no-piezo" || a == "--equivalent" || a == "--no-drift" || a == "--below" || a == "--whole-molecules" || a == "--ewald" || a == "--no-ewald" || a == "--hydrogen" || a == "--allow-unconverged" || a == "--find" || a == "--clear" || a == "--neb" || a == "--no-climb" || a == "--relax-ends" || a == "--no-match";
}

// A switch for this command? --json is one for run and provenance (print JSON) and names the output file for the
// rest (analyze, elastic, data ...: --json OUT); --vacuum is one for pore (a slit in vacuum) and a length for the
// rest (surface, interface, habit: --vacuum 15)
bool is_cli_switch_for(const std::string& cmd, const std::string& a) {
  if (a == "--json") return cmd == "run" || cmd == "provenance";
  if (a == "--vacuum") return cmd == "pore";
  return is_cli_switch(a);
}

std::map<std::string, std::string> parse(int argc, char** argv, int from, std::vector<std::string>& pos, const std::string& cmd) {
  std::map<std::string, std::string> o;
  for (int i = from; i < argc; ++i) {
    std::string a = argv[i];
    if (a.rfind("--", 0) == 0 && !known_options().count(a)) {
      const std::string near = closest_option(a);
      throw std::invalid_argument("unknown option " + a + (near.empty() ? " (caps " + cmd + " --help shows the command's options)" : " (did you mean " + near + "?)"));
    }
    if (a.rfind("--", 0) == 0 || a == "-o") {
      if (a == "--json" && !is_cli_switch_for(cmd, a) && (i + 1 >= argc || argv[i + 1][0] == '-'))
        throw std::invalid_argument("--json needs a file name here (--json OUT.json)");
      const bool flag = is_cli_switch_for(cmd, a) ||
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
      if (!r.empty() && exists_here(r + "/data/ml/caps-mace-v1.model")) { m = r + "/data/ml/caps-mace-v1.model"; break; }
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
// copolymer options (grow, interface): --r1 --r2, --reactivity 'r11,r12,r13;r21,…' (rows by ';'), --transition (the same
// layout), --penultimate r1,r1',r2,r2', --drift none|chains|along, --conversion X
std::vector<double> cli_matrix(const std::string& t, const char* what) {
  std::vector<double> v;
  size_t rows = 0, cols = 0;
  std::stringstream rs(t);
  for (std::string row; std::getline(rs, row, ';');) {
    size_t c = 0;
    std::stringstream cs(row);
    for (std::string x; std::getline(cs, x, ',');) v.push_back(std::stod(x)), ++c;
    if (rows == 0) cols = c;
    else if (c != cols) throw std::invalid_argument(std::string(what) + ": every row needs " + std::to_string(cols) + " numbers");
    ++rows;
  }
  if (rows != cols || rows == 0) throw std::invalid_argument(std::string(what) + ": a square matrix, rows split by ';'");
  return v;
}
// architecture options (grow): --architecture linear|star|comb|branched|dendrimer|graft, --arms, --arm-dp, --spacing,
// --branch-probability, --generations; grafts: --graft-units 'SMILES,…' --graft-sequence --graft-weights --graft-dp
// --graft-on K (backbone unit, from 1) --graft-probability --graft-max
void architecture_options(ChainSpec& spec, std::map<std::string, std::string>& o) {
  auto list = [](const std::string& t) {
    std::vector<std::string> v;
    std::stringstream ss(t);
    for (std::string x; std::getline(ss, x, ',');) v.push_back(x);
    return v;
  };
  if (o.count("--architecture")) spec.architecture = architecture_from_string(o["--architecture"]);
  if (o.count("--arms")) spec.arms = std::stoi(o["--arms"]);
  if (o.count("--arm-dp")) spec.arm_dp = std::stoi(o["--arm-dp"]);
  if (o.count("--spacing")) spec.spacing = std::stoi(o["--spacing"]);
  if (o.count("--branch-probability")) spec.branch_probability = std::stod(o["--branch-probability"]);
  if (o.count("--generations")) spec.generations = std::stoi(o["--generations"]);
  if (o.count("--graft-units")) for (const auto& u : list(o["--graft-units"])) spec.graft_units.push_back({u, u});
  if (o.count("--graft-sequence")) spec.graft_sequence = sequence_from_string(o["--graft-sequence"]);
  if (o.count("--graft-weights")) for (const auto& w : list(o["--graft-weights"])) spec.graft_weights.push_back(std::stod(w));
  if (o.count("--graft-dp")) spec.graft_dp = std::stoi(o["--graft-dp"]);
  if (o.count("--graft-on")) spec.graft_on = std::stoi(o["--graft-on"]) - 1;
  if (o.count("--graft-probability")) spec.graft_probability = std::stod(o["--graft-probability"]);
  if (o.count("--graft-max")) spec.graft_max = std::stoi(o["--graft-max"]);
  if (spec.architecture == Architecture::Graft && spec.graft_units.empty()) throw std::invalid_argument("--architecture graft needs --graft-units");
}

void copolymer_options(ChainSpec& spec, std::map<std::string, std::string>& o) {
  if (o.count("--r1")) spec.r1 = std::stod(o["--r1"]);
  if (o.count("--r2")) spec.r2 = std::stod(o["--r2"]);
  if (o.count("--reactivity")) spec.reactivity = cli_matrix(o["--reactivity"], "--reactivity");
  if (o.count("--transition")) spec.transition = cli_matrix(o["--transition"], "--transition");
  if (o.count("--penultimate")) {
    std::stringstream ss(o["--penultimate"]);
    for (std::string x; std::getline(ss, x, ',');) spec.penultimate.push_back(std::stod(x));
    if (spec.penultimate.size() != 4) throw std::invalid_argument("--penultimate r1,r1',r2,r2'");
  }
  if (o.count("--drift")) spec.drift = o["--drift"];
  if (o.count("--conversion")) spec.conversion = std::stod(o["--conversion"]);
  const size_t n = spec.units.size();
  if (!spec.reactivity.empty() && spec.reactivity.size() != n * n) throw std::invalid_argument("--reactivity: " + std::to_string(n) + " × " + std::to_string(n) + " for the units given");
  if (!spec.transition.empty() && spec.transition.size() != n * n) throw std::invalid_argument("--transition: " + std::to_string(n) + " × " + std::to_string(n) + " for the units given");
}
// the sequence report of a grown copolymer: printed, and whole to --sequence-report FILE.json
void print_sequence_report(const ChainSpec& spec, const GrowReport& rep, std::map<std::string, std::string>& o) {
  if (rep.sequences.empty()) return;
  const Json r = sequence_report(spec, rep.sequences);
  const int n = int(r["units"].number());
  // one unit: nothing to tabulate (composition 1, run length DP), but a --sequence-report file is still written
  if (spec.units.size() >= 2) {
    auto row = [&](const char* name, const Json& a) {
      std::printf("  %-22s", name);
      for (const auto& x : a.items()) std::printf(" %9.4f", x.is_number() ? x.number() : std::nan(""));
      std::printf("\n");
    };
    std::printf("sequence (%s), units", r["kind"].str().c_str());
    for (int k = 0; k < n; ++k) std::printf(" %9c", char('A' + k));
    std::printf("\n");
    if (r.has("model")) row("model composition", r["model"]["composition"]);
    if (r.has("drift")) row("model cumulative", r["drift"].items().back()["cumulative"]);
    row("built composition", r["built"]["composition"]);
    if (r.has("model")) row("model run length", r["model"]["run_length"]);
    row("built run length", r["built"]["run_length"]);
    if (r.has("model")) {
      double worst = 0;
      const auto& md = r["model"]["dyads"].items();
      const auto& bd = r["built"]["dyads"].items();
      for (size_t k = 0; k < md.size() && k < bd.size(); ++k) worst = std::max(worst, std::fabs(md[k].number() - bd[k].number()));
      std::printf("  largest dyad difference built − model %.4f over %.0f units\n", worst, r["built"]["count"].number());
    }
  }
  if (o.count("--sequence-report") && o["--sequence-report"] != "1") {
    std::ofstream f(o["--sequence-report"]);
    if (!f) throw std::runtime_error("cannot write " + o["--sequence-report"]);
    f << r.dump(1) << "\n";
    std::printf("wrote %s\n", o["--sequence-report"].c_str());
  }
}

// --restraints 'distance 1 2 3.0 20; position 1-10 5' or a file of them (relax, md, equilibrate, ff apply's exports)
std::vector<Restraint> cli_restraints(const System& s, std::map<std::string, std::string>& o, bool quiet) {
  if (!o.count("--restraints")) return {};
  std::string text = o["--restraints"];
  if (std::ifstream f(text); f && text.find(' ') == std::string::npos) {
    std::stringstream ss;
    ss << f.rdbuf();
    text = ss.str();
  }
  auto rs = parse_restraints(text, s);
  if (!quiet) {
    std::vector<double> x;
    for (const auto& a : s.atoms) x.insert(x.end(), a.pos.begin(), a.pos.end());
    std::vector<double> f(x.size(), 0.0);
    std::printf("%zu restraint(s), %.4f kcal/mol at the start\n", rs.size(), add_restraint_forces(rs, x, s.cell, f).energy);
  }
  return rs;
}

// One row of a CSV or tab-separated table (RFC 4180): a field in double quotes is one value whatever separators it holds,
// the quotes are not part of it, and "" inside it is a literal quote. The separator is a tab when one appears outside
// quotes, else a comma (sep 0), or the one given.
std::vector<std::string> split_csv_record(const std::string& line, char sep = 0) {
  if (!sep) {
    sep = ',';
    bool q = false;
    for (char c : line) {
      if (c == '"') q = !q;
      else if (c == '\t' && !q) {
        sep = '\t';
        break;
      }
    }
  }
  std::vector<std::string> v;
  std::string cur;
  bool q = false;
  for (size_t i = 0; i < line.size(); ++i) {
    const char c = line[i];
    if (c == '"') {
      if (q && i + 1 < line.size() && line[i + 1] == '"') cur += '"', ++i;
      else q = !q;
    } else if (c == sep && !q) {
      v.push_back(cur), cur.clear();
    } else if (c != '\r') {
      cur += c;
    }
  }
  v.push_back(cur);
  return v;
}

// A pipeline's steps as caps pipeline, caps run and caps bundle read them: JSON when the text starts with '[' or '{',
// the YAML a Studio pipeline is saved as otherwise. Either way a mistake is the user's (invalid_argument) and names
// the format it was read as.
Pipeline cli_pipeline(const std::string& text, std::string* name = nullptr, std::string* file = nullptr, std::string* topology = nullptr) {
  const auto t0 = text.find_first_not_of(" \t\r\n");
  if (t0 == std::string::npos) throw std::invalid_argument("the pipeline is empty: give a JSON list of steps or a pipeline YAML");
  if (text[t0] == '[' || text[t0] == '{') {
    try {
      return pipeline_from_json(Json::parse(text));
    } catch (const JsonError& e) {
      throw std::invalid_argument(std::string("pipeline JSON: ") + e.what());
    }
  }
  if (text.rfind("caps_pipeline", 0) != 0 && text.find("\ncaps_pipeline") == std::string::npos)
    throw std::invalid_argument("not a CAPS pipeline: a pipeline YAML starts with the line \"caps_pipeline: 1\" and lists its steps under "
                                "\"steps:\", one a line, as \"  - wrap: {}\" or \"  - cluster: {mode: cutoff, cutoff: 3.3}\" (a JSON list of steps works too)");
  return pipeline_from_yaml(text, name, file, topology);
}

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
      // A crystal cell too small for its bonds lists one bond per pair of atoms (bonds.cpp), so a polyethylene carbon
      // bonded to two images of one neighbour looks three-bonded, an alkene GAFF cannot type. The crystal commands
      // (phonons, qha, dielectric, cellrelax) type the copies that hold every bond (field_for), so GAFF is judged on
      // those copies here too and the force field named is the one they use.
      std::array<int, 3> n{1, 1, 1};
      if (s0.cell.valid() && !s0.bonds.empty() && !cell_fit(s0).fits) {
        try {
          n = copies_to_fit(s0);
        } catch (const std::invalid_argument&) {   // no repetition within 15000 atoms: judged on the cell as it is
        }
      }
      const bool copies = n != std::array<int, 3>{1, 1, 1};
      try {
        if (!copies) {
          ForceField F = assign_gaff(s0);
          if (!quiet) std::printf("force field: the built-in GAFF (C and H; the structure has no force field of its own); give --ff for another\n");
          return F;
        }
        (void)assign_gaff(supercell(s0, n[0], n[1], n[2]));   // throws FieldError when the copies are beyond it too
        if (!quiet)
          std::printf("force field: the built-in GAFF (C and H; the structure has no force field of its own), typed on the %d×%d×%d copies this cell "
                      "needs for its bonds; give --ff for another\n", n[0], n[1], n[2]);
        try {
          return assign_gaff(s0);
        } catch (const FieldError&) {   // the cell alone misses bonds: only its copies can be computed on (dynamics and relax refuse it)
          ForceField F = assign_uff(s0);
          F.name = "UFF on the cell alone (too small for its bonds; GAFF on its " + std::to_string(n[0]) + "×" + std::to_string(n[1]) + "×" +
                   std::to_string(n[2]) + " copies)";
          return F;
        }
      } catch (const FieldError&) {   // the built-in GAFF knows sp3 and aromatic carbon only (not C=C of rubbers)
      }
    }
    if (!quiet) std::printf("force field: UFF (the structure has no force field of its own%s); give --ff for another\n", ch ? "; C=C is beyond the built-in GAFF" : "");
    return assign_uff(s0);
  }
  if (!quiet) warn_keep_without_charges(o, s0);
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
  // --charges auto | increments as caps ff apply and the Studio take them (parameterize_charges), the rest as parameterize.
  // No --charges and no charges in the file: auto (the force field's own, else Gasteiger-Marsili, else QEq), the default of
  // caps ff apply, the recipe runner and the Studio, so a force field without per-type charges (GAFF, DREIDING) still runs
  std::vector<std::string> said;
  ForceField ff = parameterize_charges(s0, def, types, o.count("--charges") ? o["--charges"] : (file_charges(s0) ? "keep" : "auto"),
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
    if (z.size() != xyz.size())
      throw std::runtime_error(std::to_string(z.size()) + " elements but " + std::to_string(xyz.size()) + " positions");
    std::vector<caps::ml::Hold> holds;
    if (job.has("hold"))
      for (const Json& h : job["hold"].items()) {
        caps::ml::Hold hd;
        hd.n = static_cast<int>(h.size());
        if (hd.n < 2 || hd.n > 4) throw std::runtime_error("a held coordinate has 2, 3 or 4 atoms");
        for (int k = 0; k < hd.n; ++k) {
          const double a = h[k].number();
          if (!(a >= 0 && a < double(xyz.size()))) {
            char b[64];
            std::snprintf(b, sizeof b, "%g", a);
            throw std::runtime_error("held atom " + std::string(b) + " is not one of the " + std::to_string(xyz.size()) + " atoms (atoms count from 0)");
          }
          hd.idx[k] = static_cast<int>(a);
        }
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
    if (!root.empty() && exists_here(root + "/data/forcefields/catalogue.json")) { ro.forcefield_dir = root + "/data/forcefields"; break; }
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
  ro.frames_requested = o.count("--dump") > 0;   // md / equilibrate stages without frame_ps then write about 100 frames
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
    if (ro.frames_requested && !live.frames_open())
      std::fprintf(stderr, "caps run: --dump wrote no frames (no md or equilibrate stage ran, or each set frame_ps / frame_every to 0)\n");
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
    for (const auto& [k, v] : p.extra)
      if (std::isfinite(v)) std::printf("  %s: %.6g\n", k.c_str(), v);
      else std::printf("  %s: n/a\n", k.c_str());   // not defined (see the notes); never a bare "nan"
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

// Atoms by a selection, SEL as caps edit reads it (core parse_selection): all, element:C,N, smarts:PATTERN, type:LABEL
// (a type label or number), or 1-based numbers and ranges (1-10,15).
std::vector<uint32_t> cli_selection(const System& s, const std::string& sel) {
  std::vector<char> m;
  try {
    m = parse_selection(s, sel);
  } catch (const std::exception& e) {
    throw std::invalid_argument("selection '" + sel + "': " + e.what() + " (use all, element:SYM[,SYM], smarts:PATTERN, type:LABEL or atom numbers 1-" +
                                std::to_string(s.atoms.size()) + " such as 3,5-9)");
  }
  std::vector<uint32_t> out;
  for (uint32_t i = 0; i < m.size(); ++i)
    if (m[i]) out.push_back(i);
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

// --pme [--ewald-rtol 1e-5] [--pme-spacing 1.0] [--pme-order 5] [--slab 3]: particle-mesh Ewald instead of damped shifted force;
// --slab: a film with vacuum along z, the reciprocal cell lengthened and the dipole correction added (implies --pme)
void electrostatics(EnergyOptions& e, std::map<std::string, std::string>& o) {
  if (o.count("--pme")) e.electrostatics = EnergyOptions::Electrostatics::PME;
  if (o.count("--ewald-rtol")) e.ewald_rtol = std::stod(o["--ewald-rtol"]);
  if (o.count("--pme-spacing")) e.pme_spacing = std::stod(o["--pme-spacing"]);
  if (o.count("--pme-order")) e.pme_order = std::stoi(o["--pme-order"]);
  if (o.count("--dsf-alpha")) e.dsf_alpha = std::stod(o["--dsf-alpha"]);
  if (o.count("--slab")) {   // PME slab correction for a film with vacuum along z
    e.slab = std::stod(o["--slab"]);
    if (e.slab <= 1) throw std::invalid_argument("--slab takes the lengthening factor along z, above 1 (3 is usual)");
    e.electrostatics = EnergyOptions::Electrostatics::PME;
  }
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

// A structure written by its file extension: .pdb .xyz/.extxyz .mol2 .car .gro .sdf/.mol .cif, .vasp or a POSCAR/CONTCAR
// name, a LAMMPS data file for .data/.lmp/.lammps or a name with no extension. Any other extension is refused (a
// .lammpstrj or .top holding a data file would not read back as one).
static void write_structure_file(const System& s, const std::string& out);
// the same name check before a long build, so a run is not lost to an output name write_structure_file refuses
static void require_structure_name(const std::string& out) {
  std::string low = out;
  for (char& c : low) c = char(std::tolower(static_cast<unsigned char>(c)));
  auto ends = [&](const char* e) { const std::string x(e); return low.size() > x.size() && low.compare(low.size() - x.size(), x.size(), x) == 0; };
  for (const char* e : {".pdb", ".xyz", ".extxyz", ".mol2", ".car", ".gro", ".sdf", ".mol", ".cif", ".vasp", ".poscar", ".data", ".lmp", ".lammps"})
    if (ends(e)) return;
  if (poscar_name(out) || !std::filesystem::path(out).has_extension()) return;
  throw std::invalid_argument("output format '" + std::filesystem::path(out).extension().string() +
                              "': write .data/.lmp, .xyz/.extxyz, .pdb, .mol2, .car, .gro, .sdf/.mol, .cif, .vasp or a "
                              "POSCAR/CONTCAR name (caps frames writes trajectories such as .lammpstrj)");
}

static void write_structure_file(const System& s, const std::string& out) {
  std::string low = out;
  for (char& c : low) c = char(std::tolower(static_cast<unsigned char>(c)));
  auto ends = [&](const char* e) { const std::string x(e); return low.size() > x.size() && low.compare(low.size() - x.size(), x.size(), x) == 0; };
  if (ends(".pdb")) write_pdb(s, out);
  else if (ends(".xyz") || ends(".extxyz")) write_xyz(s, out);
  else if (ends(".mol2")) write_mol2(s, out);
  else if (ends(".car")) write_car(s, out);
  else if (ends(".gro")) write_gro(s, out);
  else if (ends(".sdf") || ends(".mol")) write_sdf(s, out);
  else if (ends(".cif")) write_cif(s, out);
  else if (ends(".vasp") || ends(".poscar") || poscar_name(out)) write_poscar(s, out);
  else if (ends(".data") || ends(".lmp") || ends(".lammps") || !std::filesystem::path(out).has_extension()) write_lammps_data(s, out);
  else
    throw std::invalid_argument("output format '" + std::filesystem::path(out).extension().string() +
                                "': write .data/.lmp, .xyz/.extxyz, .pdb, .mol2, .car, .gro, .sdf/.mol, .cif, .vasp or a "
                                "POSCAR/CONTCAR name (caps frames writes trajectories such as .lammpstrj)");
}

// write_structure_file, but a LAMMPS data name carries the force field (Coeffs sections, charges, core-shell types and
// masses) as caps relax writes it
static void write_structure_file(const System& s, const ForceField& ff, const EnergyOptions& e, const std::string& out) {
  std::string low = out;
  for (char& c : low) c = char(std::tolower(static_cast<unsigned char>(c)));
  auto ends = [&](const char* x) { const std::string y(x); return low.size() > y.size() && low.compare(low.size() - y.size(), y.size(), y) == 0; };
  const bool data = (ends(".data") || ends(".lmp") || ends(".lammps") || !std::filesystem::path(out).has_extension()) &&
                    !ends(".poscar") && !poscar_name(out);
  if (!data) return write_structure_file(s, out);
  if (const std::string why = write_lammps_data_or_structure(s, ff, e, out); !why.empty()) std::fprintf(stderr, "%s: %s\n", out.c_str(), why.c_str());
}

// True for a name write_structure_file writes as a LAMMPS data file (it carries the force field, core-shell types too)
static bool lammps_data_name(const std::string& out) {
  std::string low = out;
  for (char& c : low) c = char(std::tolower(static_cast<unsigned char>(c)));
  auto ends = [&](const char* x) { const std::string y(x); return low.size() > y.size() && low.compare(low.size() - y.size(), y.size(), y) == 0; };
  return (ends(".data") || ends(".lmp") || ends(".lammps") || !std::filesystem::path(out).has_extension()) && !ends(".poscar") && !poscar_name(out);
}

// The chemistry a force field read from a file is bound to, atom by atom: elements, charges and the bond list. A geometry
// edit (move, rotate, a torsion, the cell) leaves it as it was, so the file's force field still describes the structure;
// an edit of the elements, the bonds, the charges or the atom count does not.
struct ChemistryKey {
  std::vector<int> z;
  std::vector<double> q;
  std::vector<std::pair<uint32_t, uint32_t>> bonds;
  bool operator==(const ChemistryKey& k) const { return z == k.z && q == k.q && bonds == k.bonds; }
};
static ChemistryKey chemistry_key(const System& s) {
  ChemistryKey k;
  for (const auto& a : s.atoms) k.z.push_back(a.element), k.q.push_back(a.charge);
  for (const auto& b : s.bonds) k.bonds.emplace_back(std::min(b.i, b.j), std::max(b.i, b.j));
  std::sort(k.bonds.begin(), k.bonds.end());
  return k;
}

// caps edit / caps convert: the result as a file. A LAMMPS data output keeps the force field the input carried (its
// coefficients, read back) when the chemistry is unchanged (same_chemistry); otherwise it is the structure only, and a
// note says the force field was not carried over, so the loss is never silent.
static void write_edited_structure(const System& s, bool had_ff, bool same_chemistry, const std::string& from, const std::string& out,
                                   std::map<std::string, std::string>& o, const std::string& why_dropped) {
  if (!lammps_data_name(out)) return write_structure_file(s, out);   // a structure format: coordinates and bonds only
  if (had_ff && same_chemistry && s.forcefield && s.forcefield->charge.size() == s.atoms.size()) {
    EnergyOptions e;
    electrostatics(e, o);
    run_cutoff(e, *s.forcefield, o, true);
    return write_structure_file(s, *s.forcefield, e, out);
  }
  write_structure_file(s, out);
  if (had_ff)
    std::fprintf(stderr, "note: %s is the structure only: the force field of %s (%s) was not carried over because %s; give it again with caps ff apply %s --ff FF\n",
                 out.c_str(), from.c_str(), s.forcefield ? s.forcefield->name.c_str() : "its coefficients", why_dropped.c_str(), out.c_str());
  else if (!s.forcefield_error.empty() && s.forcefield_atoms > 0)
    std::fprintf(stderr, "note: %s is the structure only: the coefficients of %s were not read back (%s)\n", out.c_str(), from.c_str(), s.forcefield_error.c_str());
}

// True for a name write_structure_file writes in a structure format other than LAMMPS data (.cif, .pdb, POSCAR, ...)
static bool structure_format_name(const std::string& out) {
  std::string low = out;
  for (char& c : low) c = char(std::tolower(static_cast<unsigned char>(c)));
  auto ends = [&](const char* x) { const std::string y(x); return low.size() > y.size() && low.compare(low.size() - y.size(), y.size(), y) == 0; };
  for (const char* e : {".pdb", ".xyz", ".extxyz", ".mol2", ".car", ".gro", ".sdf", ".mol", ".cif", ".vasp", ".poscar"})
    if (ends(e)) return true;
  return poscar_name(out);
}

// The structure as it was read, at the positions and cell of `prepared`: the copy cli_forcefield prepared for a shell
// model (a shell appended per core) or a water model (massless sites inserted) holds atoms that are not atoms of the
// material, and a CIF, PDB or XYZ file must not list them. The atoms are matched by order (shells come after every core;
// sites have element 0); anything else (united atoms, beads) has no such match and `prepared` is returned as it is.
static System as_read(const System& in, const System& prepared) {
  const size_t n = in.atoms.size();
  if (prepared.atoms.size() <= n) return prepared;
  std::vector<size_t> keep;
  bool order = true;
  for (size_t i = 0; i < n && order; ++i) order = prepared.atoms[i].element == in.atoms[i].element;
  if (order) {
    for (size_t i = 0; i < n; ++i) keep.push_back(i);
  } else {
    for (size_t i = 0; i < prepared.atoms.size(); ++i) if (prepared.atoms[i].element != 0) keep.push_back(i);
    if (keep.size() != n) return prepared;
    for (size_t k = 0; k < n; ++k) if (prepared.atoms[keep[k]].element != in.atoms[k].element) return prepared;
  }
  System out = in;
  out.cell = prepared.cell;
  for (size_t k = 0; k < n; ++k) out.atoms[k].pos = prepared.atoms[keep[k]].pos;
  if (prepared.velocities.size() == prepared.atoms.size()) {
    out.velocities.resize(n);
    for (size_t k = 0; k < n; ++k) out.velocities[k] = prepared.velocities[keep[k]];
  }
  return out;
}

void keep_cli_error(int argc, char** argv, const std::exception& e);

// The DFT surface & adsorption workbench (caps/dft_commands.hpp): the same commands as the Studio's DFT pages
static int dft_main(const std::string& cmd, int argc, char** argv) {
  std::vector<std::string> pos;
  std::vector<std::pair<std::string, std::string>> flags;
  bool json = false;
  for (int i = 2; i < argc; ++i) {
    const std::string a = argv[i];
    if (a == "--help" || a == "-h") { std::printf("%s", dft_help(cmd).c_str()); return 0; }
    if (a == "--json") { json = true; continue; }
    const std::vector<std::string> own = dft_option_names(cmd);
    // a single dash before a one-letter option the command declares (-a), as --help prints it; a negative number stays
    // a positional
    const bool short_opt = a.size() == 2 && a[0] == '-' && std::isalpha(static_cast<unsigned char>(a[1])) &&
                           std::find(own.begin(), own.end(), a) != own.end();
    const bool stray_short = a.size() == 2 && a[0] == '-' && std::isalpha(static_cast<unsigned char>(a[1])) && a != "-o" && !short_opt;
    if (a.rfind("--", 0) == 0 || a == "-o" || short_opt || stray_short) {
      if (stray_short || (a != "-o" && !dft_known_option(a))) {
        std::vector<std::string> words = own;
        words.insert(words.end(), {"--json", "--help", "-o"});
        const std::string near = closest_word(a, words);
        std::fprintf(stderr, "caps %s: unknown option %s%s; caps %s --help lists the options\n", cmd.c_str(), a.c_str(),
                     near.empty() ? "" : (" (did you mean " + near + "?)").c_str(), cmd.c_str());
        return 2;
      }
      if (dft_is_switch(cmd, a)) flags.push_back({a, "true"});
      else if (i + 1 < argc) flags.push_back({a, argv[++i]});
      else { std::fprintf(stderr, "caps %s: %s needs a value\n", cmd.c_str(), a.c_str()); return 2; }
    } else pos.push_back(a);
  }
  std::string data;
  for (const std::string& root : caps_roots())
    if (!root.empty() && exists_here(root + "/data/sheets/sheets.json")) { data = root + "/data"; break; }
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
    keep_cli_error(argc, argv, e);
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
// the scheduler refused the script: the job and its tasks were never queued, so they are created again, not submitted
void mark_refused(const JobFolder& f) {
  Json c = Json::object();
  c["state"] = std::string("created");
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
    keep_cli_error(argc, argv, e);
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
          // checked here, before host_profile_from_json (which quietly mends an old or hand-edited file): a value
          // outside the allowed set is refused, not replaced
          if (k == "scheduler" && v != "slurm" && v != "pbs" && v != "direct")
            throw std::invalid_argument("scheduler must be slurm, pbs or direct (got '" + v + "')");
          if (k == "scratch" && v != "workspace" && v != "env" && v != "none")
            throw std::invalid_argument("scratch must be workspace, env or none (got '" + v + "')");
          if (j[k].is_number()) {
            // port, cpus, ws_days, node_cores: whole numbers, the whole value read
            std::size_t used = 0;
            long n = 0;
            try { n = std::stol(v, &used); } catch (const std::exception&) { used = 0; }
            if (v.empty() || used != v.size()) throw std::invalid_argument(k + " needs a whole number (got '" + v + "')");
            if (k == "cpus" && n < 1) throw std::invalid_argument("cpus must be at least 1 (got " + v + ")");
            if (k == "port" && (n < 1 || n > 65535)) throw std::invalid_argument("port must be 1 to 65535 (got " + v + ")");
            if (n < 0 || n > 1000000) throw std::invalid_argument(k + " must be 0 to 1000000 (got " + v + ")");
            j[k] = double(n);
          } else if (j[k].kind() == Json::Bool) {
            if (v == "1" || v == "true" || v == "yes") j[k] = true;
            else if (v == "0" || v == "false" || v == "no") j[k] = false;
            else throw std::invalid_argument(k + " takes true, false, yes, no, 1 or 0 (got '" + v + "')");
          } else j[k] = v;
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
          mark_refused(f);
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
        if (id.empty()) { mark_refused(f); out["submitted"] = false, out["error"] = msg; }
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
        tasks.push_back({"scaling", "caps md " + shq(input) + " -o bench.data --steps " + shq(steps) + " --threads " + th + " --thermo 1000 --quiet --bench-json bench.json"});
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
        if (id.empty()) mark_refused(f), out["submitted"] = false, out["error"] = msg;
        else mark_submitted(h, f, id), out["submitted"] = true, out["id"] = id;
      }
      if (a.opt.count("--json")) std::printf("%s\n", out.dump(0).c_str());
      else std::printf("%s\n%s · then: caps job collect-bench %s\n", f.dir.c_str(), out.has("id") ? ("submitted: " + out["id"].str()).c_str() : "not submitted", f.dir.c_str());
      return out.has("submitted") && !out["submitted"].boolean() ? 1 : 0;
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
      if (id.empty()) {
        mark_refused(f);
        std::fprintf(stderr, "caps job: the scheduler refused it:\n%s", msg.c_str());
        return 1;
      }
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
      // out/ is not there yet when the job was stopped before it copied back, and equivalent() throws on a missing path
      std::error_code same_ec;
      if (scr.empty() || !std::filesystem::is_directory(scr) || (std::filesystem::exists(out) && std::filesystem::equivalent(scr, out, same_ec))) {
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
    keep_cli_error(argc, argv, e);
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
    if (!root.empty() && exists_here(root + "/data/cg/mapping_rules.json")) { data = root + "/data"; break; }
  if (data.empty()) { std::fprintf(stderr, "caps: data/cg not found (set CAPS_HOME)\n"); return 2; }
  try {
    const Json r = cg_run(cmd, dft_args_from_cli(pos, flags), data);
    if (json) std::printf("%s\n", r.dump(1).c_str());
    else if (r.has("text")) std::printf("%s", r.text("text").c_str());
    else std::printf("%s\n", r.dump(1).c_str());
    return r.has("ok") && !r["ok"].boolean() ? 1 : 0;
  } catch (const std::exception& e) {
    std::fprintf(stderr, "caps %s: %s\n", cmd.c_str(), e.what());
    keep_cli_error(argc, argv, e);
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
  const std::string cmd = head.substr(head.find_last_of(' ') + 1);   // "caps tensile" -> tensile
  for (const auto& [k, v] : o) {
    if (extra.count(k) || k == "--resume" || k == "--resume-at" || k == "--print-protocol") continue;
    line += " " + k;
    if (!is_cli_switch_for(cmd, k)) line += " " + word(v);
  }
  for (const auto& [k, v] : extra) line += " " + k + (v.empty() ? "" : " " + word(v));
  return line;
}

// --pair A-B of rdf and analyze: element symbols, "all" for any element; an unknown symbol is an error (0 would be read
// as "any element" and give the all-atom g(r) under the pair's name)
std::pair<int, int> rdf_pair_elements(const std::string& p) {
  const auto dash = p.find('-');
  auto one = [&](const std::string& w) {
    if (w == "all" || w == "*") return 0;
    const int z = element_from_symbol(w);
    if (z == 0) throw std::invalid_argument("--pair " + p + ": unknown element '" + w + "' (element symbols such as C-H, or all)");
    return z;
  };
  return {one(p.substr(0, dash)), one(dash == std::string::npos ? p : p.substr(dash + 1))};
}

// caps md's run settings beside its checkpoint (run.restart.data -> run.restart.settings.json): every option of the run
// (time step, thermostat, barostat, constraints, HMR, restraints, cut-off, electrostatics, precision, seed ...) except
// those of one invocation, so caps md CHECKPOINT --resume continues with them when they are not given again
std::string md_settings_path(const std::string& checkpoint) {
  return std::filesystem::path(checkpoint).replace_extension(".settings.json").string();
}

std::string md_settings_text(const std::map<std::string, std::string>& o) {
  static const std::set<std::string> skip = {"-o", "--steps", "--resume", "--topology", "--frame", "--new-velocities", "--bench-json",
                                             "--quiet", "--threads", "--progress-file"};
  Json j = Json::object(), opts = Json::object();
  j["caps_md_settings"] = 1;
  for (const auto& [k, v] : o)
    if (!skip.count(k)) opts[k] = v;
  j["options"] = opts;
  return j.dump(1) + "\n";
}

// --resume: the stored settings fill every option not given again. A change of time step, constraints or HMR breaks
// the trajectory, so it is refused; other changes are said. A checkpoint with no settings file (an older CAPS) runs
// with the options given now, and says which.
void md_resume_settings(const std::string& checkpoint, const System& s, std::map<std::string, std::string>& o) {
  const DynamicsOptions def;
  const std::string sp = md_settings_path(checkpoint);
  std::ifstream f(sp);
  if (!f) {
    std::fprintf(stderr,
                 "note: no %s beside the checkpoint, so the first run's settings are unknown; the run continues with the options given "
                 "now (--dt %s, --thermostat %s, --barostat %s, --constraints %s%s): give them again if the first run used others\n",
                 sp.c_str(), o.count("--dt") ? o["--dt"].c_str() : "1", o.count("--thermostat") ? o["--thermostat"].c_str() : "bussi",
                 o.count("--barostat") ? o["--barostat"].c_str() : "none", o.count("--constraints") ? o["--constraints"].c_str() : "none",
                 o.count("--hmr") ? (", --hmr " + o["--hmr"]).c_str() : "");
    return;
  }
  std::stringstream ss;
  ss << f.rdbuf();
  Json j;
  try {
    j = Json::parse(ss.str());
  } catch (const JsonError& e) {
    throw std::invalid_argument(sp + ": " + e.what());
  }
  if (!j.is_object() || !j.has("options") || !j["options"].is_object()) throw std::invalid_argument(sp + " is not a caps md settings file");
  std::map<std::string, std::string> stored;
  for (const auto& [k, v] : j["options"].members()) stored[k] = v.is_string() ? v.str() : v.dump(0);
  // the checkpoint carries the force field it ran with (its coefficients and charges): used as it is, not retyped
  const bool own_ff = s.forcefield && s.forcefield->charge.size() == s.atoms.size();
  static const std::set<std::string> ff_keys = {"--ff", "--overlay", "--charges", "--typing", "--complete"};
  auto same = [](const std::string& a, const std::string& b) {
    if (a == b) return true;
    try {
      size_t ia = 0, ib = 0;
      const double x = std::stod(a, &ia), y = std::stod(b, &ib);
      return ia == a.size() && ib == b.size() && x == y;
    } catch (const std::exception&) {
      return false;
    }
  };
  // the settings a resumed run cannot change; a value not stored is the first run's default
  auto as_run = [&](const std::string& k, const std::string& v) {
    if (k == "--constraints") return std::to_string(int(constraints_from_string(v)));
    if (k == "--constraint-solver") return std::to_string(int(constraint_algorithm_from_string(v)));
    return v;
  };
  const std::map<std::string, std::string> fixed_default = {{"--dt", std::to_string(def.dt)},
                                                            {"--constraints", std::to_string(int(def.constraints))},
                                                            {"--constraint-solver", std::to_string(int(def.constraint_algorithm))},
                                                            {"--hmr", std::to_string(def.hmr)}};
  for (const auto& [k, dv] : fixed_default) {
    if (!o.count(k)) continue;
    const std::string was = stored.count(k) ? as_run(k, stored[k]) : dv;
    if (!same(as_run(k, o[k]), was))
      throw std::invalid_argument(k + " " + o[k] + " is not the checkpoint's run (" + k + " " + (stored.count(k) ? stored[k] : "not given") +
                                  "): a resumed run keeps it; leave out --resume to start a new run from this state");
  }
  std::string restored, changed;
  for (const auto& [k, v] : stored) {
    if (own_ff && ff_keys.count(k)) continue;
    if (!o.count(k)) {
      o[k] = v;
      restored += " " + k + (is_cli_switch_for("md", k) ? "" : " " + v);
    } else if (!same(o[k], v)) {
      changed += " " + k + " " + o[k] + " (was " + v + ")";
    }
  }
  if (!restored.empty()) std::printf("the checkpoint's run settings:%s\n", restored.c_str());
  if (!changed.empty()) std::printf("changed from the checkpoint's run:%s\n", changed.c_str());
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

// ---- problem reports: the newest failed command in the CAPS logs folder, beside the Studio's crash.log (macOS
// ~/Library/Application Support/CAPS/logs, Linux $XDG_DATA_HOME or ~/.local/share/CAPS/logs, Windows
// %LOCALAPPDATA%\CAPS\logs); caps report prints it with the version and system for support@caps-studio.org
constexpr const char* kSupportEmail = "support@caps-studio.org";

std::filesystem::path caps_log_dir() {
  auto env = [](const char* k) { const char* v = std::getenv(k); return std::string(v ? v : ""); };
#if defined(_WIN32)
  return std::filesystem::path(env("LOCALAPPDATA")) / "CAPS" / "logs";
#elif defined(__APPLE__)
  return std::filesystem::path(env("HOME")) / "Library" / "Application Support" / "CAPS" / "logs";
#else
  const auto xdg = env("XDG_DATA_HOME");
  return (xdg.empty() ? std::filesystem::path(env("HOME")) / ".local" / "share" : std::filesystem::path(xdg)) / "CAPS" / "logs";
#endif
}

void log_cli_error(int argc, char** argv, const std::string& what) {
  try {
    const auto dir = caps_log_dir();
    std::filesystem::create_directories(dir);
    std::ofstream f(dir / "cli-error.log", std::ios::trunc);
    const auto now = std::time(nullptr);
    char when[32];
    std::strftime(when, sizeof when, "%Y-%m-%d %H:%M:%S", std::localtime(&now));
    f << "---- " << when << " · caps " << version_string() << " (commit " << commit_string() << ")\n$";
    for (int i = 0; i < argc; ++i) f << ' ' << argv[i];
    f << "\ncaps: " << what << "\n";
  } catch (...) {}
}

// a command that failed: a wrong argument is the user's to fix; anything else may be a bug in CAPS, so it is kept
// for caps report (every command's own catch calls this, so the log always holds the newest failure)
void keep_cli_error(int argc, char** argv, const std::exception& e) {
  if (!dynamic_cast<const std::invalid_argument*>(&e)) log_cli_error(argc, argv, e.what());
}

std::string system_description() {
#if defined(_WIN32)
  return "Windows";
#else
  struct utsname u {};
  if (uname(&u) != 0) return "unknown";
  return std::string(u.sysname) + " " + u.release + " " + u.machine;
#endif
}

std::string url_escape(const std::string& s) {
  std::string o;
  char b[4];
  for (unsigned char c : s) {
    if (std::isalnum(c) || c == '-' || c == '_' || c == '.' || c == '~') o += char(c);
    else { std::snprintf(b, sizeof b, "%%%02X", c); o += b; }
  }
  return o;
}

int report_main(int argc, char** argv) {
  const bool email = argc >= 3 && std::string(argv[2]) == "--email";
  std::ostringstream r;
  r << "What I was doing when it happened:\n\n\n---- CAPS (filled in by caps report)\n"
    << "caps " << version_string() << " · commit " << commit_string() << "\nSystem " << system_description() << "\n";
  std::ifstream last(caps_log_dir() / "cli-error.log");
  if (last) r << "\n---- the last command that failed\n" << std::string(std::istreambuf_iterator<char>(last), {});
  const std::string text = r.str();
  std::printf("%s\nSend this to %s (caps report --email opens it in your mail program).\n", text.c_str(), kSupportEmail);
  if (!email) return 0;
  const std::string url = std::string("mailto:") + kSupportEmail + "?subject=" + url_escape(std::string("CAPS problem report (caps ") +
                          version_string() + ")") + "&body=" + url_escape(text.substr(0, 1600));
#if defined(_WIN32)
  const std::string open = "start \"\" \"" + url + "\"";
#elif defined(__APPLE__)
  const std::string open = "open '" + url + "'";
#else
  const std::string open = "xdg-open '" + url + "' >/dev/null 2>&1";
#endif
  if (std::system(open.c_str()) != 0) {
    std::fprintf(stderr, "caps: no mail program opened; copy the report above into an email to %s\n", kSupportEmail);
    return 1;
  }
  return 0;
}

// The commands main() runs itself (the DFT, coarse-graining and job commands have their own tables): a name not here
// is refused before any file is read. Keep in step with the `cmd == "…"` branches below.
const std::set<std::string>& main_commands() {
  static const std::set<std::string> k = {
    "adp", "analyze", "bench", "blend", "build", "bundle", "cbmc", "cellrelax", "chains", "check", "cif", "contacts",
    "convert", "copolymer", "crystal", "data", "dendrimer", "dielectric", "domains", "dpd", "edit", "elastic",
    "enumerate", "equilibrate", "exchange", "ff", "field", "fields", "fit", "frames", "graft", "grow", "habit", "index",
    "info", "interface", "martini", "md", "ml", "nano", "pack", "path", "phonons", "pipeline", "planes", "pointgroup",
    "pore", "pores", "powder", "provenance", "pull", "qeq", "qha", "qspr", "rdf", "react", "refine", "relax", "render",
    "reproduce", "run", "semicrystalline", "shape", "solvate", "solve", "sqs", "surface", "symmetry", "table", "tensile",
    "tension", "tg", "thermal", "zmelt"};
  return k;
}

// caps CMD --help (or caps help CMD): the command's own help where it has one, else its entries of the usage text
int command_help(const std::string& cmd, char* argv0) {
  if (cmd == "job") {
    char help[] = "--help", job[] = "job";
    char* a[] = {argv0, job, help, nullptr};
    return job_main(3, a);
  }
  if (is_dft_command(cmd)) { std::printf("%s", dft_help(cmd).c_str()); return 0; }
  if (is_cg_command(cmd)) { std::printf("%s", cg_help(cmd).c_str()); return 0; }
  const std::string u = usage_of(cmd);
  if (u.empty()) { std::printf("%s", usage_text()); return 0; }
  std::printf("%s(caps --help lists every command)\n", u.c_str());
  return 0;
}

int main(int argc, char** argv) {
  if (argc >= 2 && std::string(argv[1]) == "report") return report_main(argc, argv);
  if (argc >= 2) {
    const std::string a1 = argv[1];
    if (a1 == "--help" || a1 == "-h" || a1 == "help") {
      if (argc >= 3) return command_help(argv[2], argv[0]);
      std::printf("%s", usage_text());
      return 0;
    }
    if ((a1 == "--version" || a1 == "version" || a1 == "--check-update") && argc > 2) {
      std::fprintf(stderr, "caps: %s takes no arguments\n", a1.c_str());
      return 2;
    }
  }
  if (argc == 2 && (std::string(argv[1]) == "--version" || std::string(argv[1]) == "version")) {
    // "caps 0.1.0 (commit a1b2c3d4e5)": the Studio compares it with its own before running jobs on a host
    std::printf("caps %s (commit %s)\n", version_string(), commit_string());
    return 0;
  }
  if (argc == 2 && std::string(argv[1]) == "--check-update") return check_update_main();
  if (argc >= 2 && std::string(argv[1]) == "job") return job_main(argc, argv);
  if (argc >= 2 && is_dft_command(argv[1])) return dft_main(argv[1], argc, argv);
  if (argc >= 2 && is_cg_command(argv[1])) return cg_main(argv[1], argc, argv);
  if (argc >= 2 && !main_commands().count(argv[1])) {
    // an unknown command is said to be one before any file is opened (a missing file is not blamed on CAPS)
    const std::string c = argv[1];
    std::vector<std::string> all(main_commands().begin(), main_commands().end());
    all.insert(all.end(), dft_commands().begin(), dft_commands().end());
    all.insert(all.end(), cg_commands().begin(), cg_commands().end());
    all.insert(all.end(), {"job", "report"});
    std::string near = closest_word(c, all);
    if (near.empty() && c.size() >= 3)   // the start of a longer name: adsorb, adsorb-dft
      for (const auto& k : all)
        if (k.rfind(c + "-", 0) == 0) { near = k; break; }
    std::fprintf(stderr, "caps: unknown command %s%s; caps --help lists the commands\n", c.c_str(),
                 near.empty() ? "" : (" (did you mean " + near + "?)").c_str());
    return 2;
  }
  for (int i = 2; i < argc; ++i)
    if (std::strcmp(argv[i], "--help") == 0 || std::strcmp(argv[i], "-h") == 0) return command_help(argv[1], argv[0]);
  if (argc >= 2) g_usage_cmd = argv[1];
  if (argc < 3) return usage();
  const std::string cmd = argv[1];
  std::vector<std::string> pos;
  std::map<std::string, std::string> o;
  try {
    o = parse(argc, argv, 2, pos, cmd);
  } catch (const std::exception& e) {
    std::fprintf(stderr, "caps: %s\n", e.what());
    return 2;
  }
  if (cmd == "grow") {
    try {
      // a number option read whole: "--dp abc" or "--dp 3x" names the option rather than passing or failing as "stoi"
      auto num_opt = [&](const char* key, bool whole) {
        const std::string& v = o[key];
        size_t pos = 0;
        double x = 0;
        try { x = whole ? double(std::stoll(v, &pos)) : std::stod(v, &pos); } catch (const std::exception&) { pos = 0; }
        if (pos == 0 || pos != v.size()) throw std::invalid_argument(std::string(key) + ": expected " + (whole ? "a whole number" : "a number") + ", got '" + v + "'");
        return x;
      };
      auto int_opt = [&](const char* key) { return int(std::clamp(num_opt(key, true), -2147483648.0, 2147483647.0)); };
      GrowOptions g;
      if (o.count("--chains")) g.chains = int_opt("--chains");
      if (o.count("--dp")) g.dp = int_opt("--dp");
      if (o.count("--box")) g.box = num_opt("--box", false);
      g.density = o.count("--density") ? num_opt("--density", false) : 0.5;
      if (o.count("--tacticity")) g.tacticity = tacticity_from_string(o["--tacticity"]);
      if (o.count("--seed")) g.seed = std::stoull(o["--seed"]);
      g.curve = !o.count("--trans");
      if (o.count("--scale")) g.contact_scale = num_opt("--scale", false);
      if (o.count("--trials")) g.trials = int_opt("--trials");
      if (o.count("--comfortable")) g.comfortable = num_opt("--comfortable", false);
      g.escalate = o.count("--escalate");
      g.auto_scale = o.count("--auto-scale") || !o.count("--scale");   // automatic unless --scale fixes it
      if (!o.count("-o")) return usage();
      // the core would build DP 1 or one chain for these and still report the number asked for
      if (o.count("--chains") && g.chains < 1) throw std::invalid_argument("--chains must be at least 1 (got " + o["--chains"] + ")");
      if (o.count("--dp") && g.dp < 1) throw std::invalid_argument("--dp must be at least 1 (got " + o["--dp"] + ")");
      require_structure_name(o["-o"]);
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
        // several units and no sequence: all of them, as caps graft does (homopolymer would keep unit A only)
        if (spec.units.size() > 1 && !o.count("--sequence")) spec.sequence = spec.weights.empty() ? Sequence::Alternating : Sequence::Random;
        copolymer_options(spec, o);
        architecture_options(spec, o);
        if (o.count("--blocks")) for (const auto& b : split(o["--blocks"])) spec.blocks.push_back(std::stoi(b));
        if (o.count("--pattern")) spec.pattern = o["--pattern"];
        if (o.count("--ff")) spec.forcefield = ff_path(o["--ff"]);
        if (o.count("--method")) {
          const std::string m = o["--method"];
          if (m != "trials" && m != "rosenbluth" && m != "rosenbluth_lj") throw std::invalid_argument("--method trials|rosenbluth|rosenbluth_lj");
          g.method = m == "rosenbluth" ? 1 : m == "rosenbluth_lj" ? 2 : 0;
        }
        if (o.count("--threads")) g.threads = int_opt("--threads");
        if (o.count("--lookahead")) g.lookahead = std::clamp(int_opt("--lookahead"), 1, 4);
        s = grow_chains(spec, g, &rep);
        print_sequence_report(spec, rep, o);
      } else {
        if (o.count("--sequence-report"))
          throw std::invalid_argument("--sequence-report needs --units (the built-in chain has no repeat-unit sequence)");
        s = grow(g, &rep);
      }
      const std::string out = o["-o"];
      write_structure_file(s, out);   // .data, .pdb, .xyz, .mol2, .gro, .car, .sdf, .cif … by the name's extension
      for (const auto& n : rep.notes) std::printf("%s\n", n.c_str());
      char margin[64];
      if (std::isfinite(rep.worst_margin)) std::snprintf(margin, sizeof margin, "%.2f Å", rep.worst_margin);
      else std::snprintf(margin, sizeof margin, "n/a (no non-bonded pairs)");
      std::printf("%zu atoms · %zu bonds · %d restarts · %d backtracks · closest contact margin %s · wrote %s\n", s.atoms.size(),
                  s.bonds.size(), rep.restarts, rep.backtracks, margin, out.c_str());
      return 0;
    } catch (const std::exception& e) {
      std::fprintf(stderr, "caps grow: %s\n", e.what());
      keep_cli_error(argc, argv, e);
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
        if (o.count("--labels")) {
          std::ofstream f(o["--labels"]);
          if (!f) throw std::runtime_error("cannot write " + o["--labels"]);
          static const char* names[] = {"none", "stem", "loop", "bridge", "tail"};
          size_t count[5] = {0, 0, 0, 0, 0};
          f << "# CAPS semicrystalline labels of " << pos[0] << ": atom id (as in the file), class, segment (0-based, -1 none); side atoms as their backbone atom\n";
          for (size_t i = 0; i < s.atoms.size(); ++i) {
            const int k = r.atom_class[i];
            ++count[k];
            f << (s.atoms[i].id ? s.atoms[i].id : int64_t(i + 1)) << " " << names[k] << " " << r.atom_segment[i] << "\n";
          }
          std::printf("wrote %s (stem %zu, loop %zu, bridge %zu, tail %zu, none %zu atoms)\n", o["--labels"].c_str(), count[1], count[2], count[3], count[4], count[0]);
        }
        return 0;
      }
      SemicrystallineOptions so;
      if (o.count("--polymer")) so.polymer = o["--polymer"];
      for (const auto& root : caps_roots())   // the polymer crystal library
        if (exists_here(root + "/data/crystals/polymers/polymers.json")) { so.data_dir = root + "/data"; break; }
      if (so.polymer == "list") {   // the crystals the builder takes, with their defaults
        for (const auto& p : semicrystalline_polymers(so.data_dir))
          std::printf("%-18s %-48s %-11s %d,%d  crystal %.1f nm, amorphous %.1f nm, %.3f g/cm³\n", p.id.c_str(), p.name.c_str(), p.forcefield.c_str(), p.nx, p.ny,
                      p.crystal_nm, p.amorphous_nm, p.density);
        return 0;
      }
      if (!o.count("-o")) return usage();
      // the crystal's defaults (semicrystalline.json), each overridden by its option
      SemicrystallinePolymer def;
      const bool known = semicrystalline_polymer(so.polymer, so.data_dir, &def);
      if (known) {
        so.nx = def.nx, so.ny = def.ny, so.crystal_nm = def.crystal_nm, so.amorphous_nm = def.amorphous_nm;
        so.loops = def.loops, so.bridges = def.bridges, so.tails = def.tails, so.density = def.density;
      }
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
      // --ff auto (the default for a library crystal): the force field semicrystalline.json names for it
      if ((!o.count("--ff") && known) || (o.count("--ff") && o["--ff"] == "auto")) {
        if (!known) throw std::invalid_argument("--ff auto: " + so.polymer + " is not in the library's list (caps semicrystalline --polymer list); give --ff");
        o["--ff"] = def.forcefield;
        if (!quiet) std::printf("force field for %s: %s (data/crystals/polymers/semicrystalline.json; --ff for another)\n", def.id.c_str(), def.forcefield.c_str());
      }
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
        // the push-off and two minimisations cost about 1.5 min per 1000 atoms on one thread (1844 atoms of nylon 6 with
        // PCFF: 159 s), so a library-size cell takes tens of minutes: said before it starts, with the way round it
        if (!quiet && s.atoms.size() > 3000)
          std::fprintf(stderr, "relaxing %zu atoms in CAPS (%s): roughly %.0f min on one thread, less with more threads; --no-relax skips it%s\n",
                       s.atoms.size(), ff.name.c_str(), std::max(1.0, 1.5 * double(s.atoms.size()) / 1000.0),
                       o.count("--lammps-input") ? " (the --lammps-input deck minimises the cell itself)" : "");
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
      keep_cli_error(argc, argv, e);
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
      keep_cli_error(argc, argv, e);
      return 1;
    }
  }
  if (cmd == "provenance") {
    if (pos.empty()) return usage();
    std::string why;   // missing, not JSON, or not a caps-manifest/1.0 file
    const auto m = read_manifest(pos[0], &why);
    if (!m) {
      std::fprintf(stderr, "caps provenance: %s\n", why.c_str());
      log_cli_error(argc, argv, why);   // kept for caps report, as every failed command is
      return 1;
    }
    if (o.count("--json")) { std::printf("%s\n", manifest_json(*m).dump(2).c_str()); return 0; }
    if (o.count("--bibtex")) {
      const std::string b = bibtex(all_cites(*m));
      std::printf("%s", b.c_str());
      if (b.empty()) {
        std::string engines;
        for (const auto& st : m->steps) engines += (engines.empty() ? "" : ", ") + st.engine;
        std::fprintf(stderr, "caps provenance: no citations recorded for %s\n", engines.empty() ? "a manifest with no steps" : ("these steps (" + engines + ")").c_str());
      }
      return 0;
    }
    if (o.count("--methods")) {
      std::vector<std::string> refs;
      const std::string text = methods_text(*m, &refs);
      if (text.empty()) { std::fprintf(stderr, "caps provenance: the manifest has no steps to describe\n"); return 0; }
      std::printf("%s\n\n", text.c_str());
      for (size_t k = 0; k < refs.size(); ++k) std::printf("[%zu] %s\n", k + 1, refs[k].c_str());
      return 0;
    }
    if (auto c = o.find("--compare"); c != o.end()) {
      const auto other = read_manifest(c->second, &why);
      if (!other) {
        std::fprintf(stderr, "caps provenance: %s\n", why.c_str());
        log_cli_error(argc, argv, why);
        return 1;
      }
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
        if (f.supercell[0] * f.supercell[1] * f.supercell[2] > 1) {
          const auto& c = f.cell;
          std::printf("the file is a %d x %d x %d supercell of a = %.4f, b = %.4f, c = %.4f Å; the sites below are in that cell\n", f.supercell[0], f.supercell[1], f.supercell[2], norm(c.a), norm(c.b), norm(c.c));
        }
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
      keep_cli_error(argc, argv, e);
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
      keep_cli_error(argc, argv, e);
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
        po.whole_molecules = o.count("--whole-molecules") > 0;
        if (o.count("--hkl")) {
          std::string v = o["--hkl"];
          for (auto& ch : v) if (ch == ' ') ch = ',';
          if (std::sscanf(v.c_str(), "%d,%d,%d", &po.hkl[0], &po.hkl[1], &po.hkl[2]) != 3) throw std::invalid_argument("--hkl h,k,l");
        }
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
        if (o.count("--ff")) spec.forcefield = ff_path(o["--ff"]);
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
      keep_cli_error(argc, argv, e);
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
      // one weight per component, each a number of 0 or more (a missing one is not taken as 1, a negative one not as 0)
      if (o.count("--weights") && w.size() != smi.size())
        throw std::invalid_argument("--weights needs one value per component (" + std::to_string(smi.size()) + " components, " +
                                    std::to_string(w.size()) + " weights)");
      std::vector<double> wt(smi.size(), 1.0);
      for (size_t k = 0; k < w.size(); ++k) {
        size_t used = 0;
        try { wt[k] = std::stod(w[k], &used); } catch (const std::exception&) { used = 0; }
        if (used == 0 || used != w[k].size() || !std::isfinite(wt[k]) || wt[k] < 0)
          throw std::invalid_argument("--weights: '" + w[k] + "' is not a weight of 0 or more");
      }
      if (!w.empty() && *std::max_element(wt.begin(), wt.end()) <= 0) throw std::invalid_argument("--weights: at least one weight must be above 0");
      for (size_t k = 0; k < smi.size(); ++k) {
        BlendComponent c;
        c.spec.units.push_back({smi[k], smi[k]});
        c.spec.dp = o.count("--dp") ? std::stoi(o["--dp"]) : 20;
        if (o.count("--ff")) c.spec.forcefield = ff_path(o["--ff"]);
        c.weight = wt[k];
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
      keep_cli_error(argc, argv, e);
      return 1;
    }
  }
  if (cmd == "interface") {
    // a polymer film grown onto a slab: a slab file, or a crystal (CIF) cleaved with the surface options
    try {
      if (pos.empty() || !o.count("-o") || !o.count("--units")) return usage();
      // a word left over (a value after an option that takes none) is refused, not dropped
      if (pos.size() > 1) throw std::invalid_argument("unexpected argument " + pos[1] + " (caps interface takes one SLAB or CRYSTAL.cif)");
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
      if (spec.units.size() > 1 && !o.count("--sequence")) spec.sequence = spec.weights.empty() ? Sequence::Alternating : Sequence::Random;   // all the units given
      copolymer_options(spec, o);
      if (o.count("--blocks")) for (const auto& b : list(o["--blocks"])) spec.blocks.push_back(std::stoi(b));
      if (o.count("--pattern")) spec.pattern = o["--pattern"];
      if (o.count("--ff")) spec.forcefield = ff_path(o["--ff"]);
      InterfaceOptions io;
      if (o.count("--film")) io.film = std::stod(o["--film"]);
      if (o.count("--film-density")) io.density = std::stod(o["--film-density"]), io.density_given = true;
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
      keep_cli_error(argc, argv, e);
      return 1;
    }
  }
  if (cmd == "surface") {
    // a slab cleaved from a crystal (CIF), with its terminations
    try {
      if (pos.empty()) return usage();
      // a word left over (a value after an option that takes none) is refused, not dropped
      if (pos.size() > 1) throw std::invalid_argument("unexpected argument " + pos[1] + " (caps surface takes one CRYSTAL.cif)");
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
      if (o.count("--list") && !o.count("-o")) return 0;   // the list alone; with -o the slab is written as well
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
      keep_cli_error(argc, argv, e);
      return 1;
    }
  }
  if (cmd == "run") {   // a saved pipeline over many inputs: caps run PIPE.yaml --input 'runs/*/X.lammpstrj' [--frame first|last] [--csv OUT] [--out DIR]
    if (pos.empty()) return usage();
    try {
      std::ifstream pf(pos[0]);
      if (!pf) throw std::runtime_error("cannot read " + pos[0]);
      const std::string text((std::istreambuf_iterator<char>(pf)), std::istreambuf_iterator<char>());
      if (text.find("caps_pipeline") == std::string::npos) {   // a recipe: caps run RECIPE.yaml [--seed N] [--threads N] [--out DIR] [--json]
        Json r;
        try {
          const auto t0 = text.find_first_not_of(" \t\r\n");
          r = t0 != std::string::npos && (text[t0] == '{' || text[t0] == '[') ? Json::parse(text) : yaml_parse(text);   // '[': a JSON list of steps
        } catch (const std::exception& e) {
          std::fprintf(stderr, "caps run: %s: %s\n", pos[0].c_str(), e.what());
          keep_cli_error(argc, argv, e);
          return 2;
        }
        if (r.is_object() && (r.has("recipe") || r.has("build"))) return cli_recipe(r, pos[0], o);
      }
      std::string name, file, topo;
      Pipeline pl = cli_pipeline(text, &name, &file, &topo);
      if (o.count("--branch")) pl.branch = o["--branch"];   // the branch to run (steps of the others are skipped)
      // relative paths in the pipeline's source: block are read next to the pipeline file first, then from here
      const std::filesystem::path pipe_dir = std::filesystem::path(pos[0]).parent_path();
      const std::string pattern = o.count("--input") ? o["--input"] : file;
      if (pattern.empty()) throw std::runtime_error(pos[0] + " names no input: give --input FILE or a source: file: line");
      auto inputs = glob_files(pattern);
      if (inputs.empty() && !o.count("--input") && !pipe_dir.empty() && std::filesystem::path(pattern).is_relative())
        inputs = glob_files((pipe_dir / pattern).string());
      if (inputs.empty()) throw std::runtime_error("no files match " + pattern);
      // the topology: --topology, else the pipeline's source: topology:, else STEM.data next to each trajectory
      std::string topology = o.count("--topology") ? o["--topology"] : "";
      if (topology.empty() && !topo.empty()) {
        topology = topo;
        if (std::filesystem::path(topo).is_relative() && !pipe_dir.empty() && std::filesystem::exists(pipe_dir / topo)) topology = (pipe_dir / topo).string();
        if (!std::filesystem::exists(topology)) throw std::runtime_error("the topology " + topo + " named in " + pos[0] + " does not exist");
      }
      const bool last = !o.count("--frame") || o["--frame"] != "first";
      std::vector<std::string> keys;
      std::vector<std::pair<std::string, std::map<std::string, double>>> rows;
      // the pipeline's outputs: under --out (default outputs/), one folder per input when there are several
      const std::string out_root = o.count("--out") ? o["--out"] : "outputs";
      size_t idx = 0;
      for (const auto& in : inputs) {
        ++idx;
        try {
          std::string tp = topology;
          const std::string ext = std::filesystem::path(in).extension().string();
          const bool dump = ext == ".lammpstrj" || ext == ".dump";
          if (tp.empty() && dump) {
            const auto d = std::filesystem::path(in).replace_extension(".data");
            if (std::filesystem::exists(d)) tp = d.string();
          }
          const Trajectory t = open_file(in, tp);
          const int fr = last ? int(t.frames()) - 1 : 0;
          System f0 = t.frame(size_t(fr));
          if (tp.empty() && dump && !f0.atoms.empty() && f0.total_mass() <= 0)
            std::fprintf(stderr, "note    %s has no topology: no masses or bonds, so Mass and Density are 0 and bonds are guessed (give --topology DATA)\n", in.c_str());
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
    } catch (const std::exception& e) {
      std::fprintf(stderr, "caps run: %s\n", e.what());
      keep_cli_error(argc, argv, e);
      return 2;
    }
  }
  if (cmd == "reproduce") {   // a figure bundle: rebuild its data from its input and pipeline, compare the hashes
    if (pos.empty()) return usage();
    try {
      std::vector<std::string> report;
      const bool ok = reproduce_bundle(pos[0], report);
      for (const auto& r : report) std::printf("%s\n", r.c_str());
      std::printf("%s\n", ok ? "reproduced: every data file matches" : "not reproduced");
      return ok ? 0 : 2;
    } catch (const std::exception& e) {
      std::fprintf(stderr, "caps reproduce: %s\n", e.what());
      keep_cli_error(argc, argv, e);
      return 2;
    }
  }
  if (cmd == "bench") {
    try {
      BenchOptions b;
      // the shipped data (caps_roots)
      auto find = [](const std::string& rel, const std::string& probe) {
        for (const auto& r : caps_roots())
          if (exists_here(r + "/" + rel + "/" + probe)) return r + "/" + rel;
        return rel;
      };
      b.samples = o.count("--samples") ? o["--samples"] : find("samples", "ps_melt.data");
      b.forcefields = o.count("--forcefields") ? o["--forcefields"] : find("data/forcefields", "catalogue.json");
      if (o.count("--repeats")) b.repeats = std::stoi(o["--repeats"]);
      b.quick = o.count("--quick") > 0;
      // T6 / T7, as the C API sets them: the reference ranges beside the force-field library (data/reference/
      // polymers.json), the equilibrated cells from --cells, else CAPS_BENCH_CELLS, else ~/.caps/bench/cells
      b.reference = o.count("--reference") ? o["--reference"]
                                           : (std::filesystem::path(b.forcefields).parent_path() / "reference" / "polymers.json").string();
      if (o.count("--cells")) b.cells = o["--cells"];
      else if (const char* c = std::getenv("CAPS_BENCH_CELLS"); c && *c) b.cells = c;
      else if (const char* h = std::getenv("HOME"); h && *h) b.cells = std::string(h) + "/.caps/bench/cells";
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
      keep_cli_error(argc, argv, e);
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
      keep_cli_error(argc, argv, e);
      return 1;
    }
  }
  if (cmd == "build") {
    try {
      if (pos.empty() || !o.count("-o")) return usage();
      BuildOptions b;
      if (o.count("--conformers")) b.conformers = std::stoi(o["--conformers"]);
      if (o.count("--seed")) b.seed = std::stoull(o["--seed"]);
      if (o.count("--ff")) {
        // a library name (oplsaa2024, gaff2, pcff …) or a path, as every other command takes it; a name that is
        // neither stops here rather than cleaning up with UFF in its place
        b.forcefield = ff_path(o["--ff"]);
        if (!b.forcefield.empty() && !is_uff(b.forcefield) && !exists_here(b.forcefield))
          throw std::runtime_error("--ff " + o["--ff"] + ": no such force field (a library name in data/forcefields, a path, or uff)");
      }
      const BuildResult r = build_molecule(pos[0], b);
      const std::string out = o["-o"];
      auto ends = [&](const char* e) { return out.size() > 4 && out.substr(out.size() - std::strlen(e)) == e; };
      auto write = [&](const System& s, const std::string& path) {
        if (ends(".pdb")) write_pdb(s, path);
        else if (ends(".xyz")) write_xyz(s, path);
        else if (ends(".mol2")) write_mol2(s, path);
        else write_lammps_data(s, path);
      };
      std::vector<std::string> written;   // the files actually written, for the report
      if (o.count("--all") && r.conformers.size() > 1) {
        // confs.pdb -> confs_1.pdb ... (the number goes on the file's stem, in the same folder; a name may have no extension)
        const std::filesystem::path op(out);
        const std::string stem = op.stem().string(), ext = op.extension().string();
        for (size_t k = 0; k < r.conformers.size(); ++k) {
          const std::string p = (op.parent_path() / (stem + "_" + std::to_string(k + 1) + ext)).string();
          write(molecule_system(r.graph, r.conformers[k].pos), p);
          written.push_back(p);
        }
      } else {
        write(r.system, out);
        written.push_back(out);
      }
      // stereo: every centre (and stereogenic double bond) in the built structure, with how many the SMILES left open
      const int centres = r.info.stereocentres + r.open_centres, dbonds = r.info.stereo_bonds + r.open_bonds;
      std::string stereo = std::to_string(centres) + (centres == 1 ? " stereocentre" : " stereocentres");
      if (r.open_centres) stereo += " (" + std::to_string(r.open_centres) + " unspecified)";
      if (dbonds) {
        stereo += " · " + std::to_string(dbonds) + (dbonds == 1 ? " stereo double bond" : " stereo double bonds");
        if (r.open_bonds) stereo += " (" + std::to_string(r.open_bonds) + " unspecified)";
      }
      std::printf("%s · %.2f g/mol · %d atoms (%d heavy) · %d rings · %s\n", r.info.formula.c_str(), r.info.mass, r.info.atoms,
                  r.info.heavy, r.info.rings, stereo.c_str());
      std::printf("%s\n", r.method.c_str());
      const double e0 = r.conformers.front().energy;
      for (size_t k = 0; k < r.conformers.size(); ++k)
        std::printf("  conformer %zu  %s\n", k + 1,
                    r.conformers[k].minimised ? (std::to_string(r.conformers[k].energy - e0).substr(0, 6) + " kcal/mol").c_str() : "not minimised");
      for (const auto& n : r.notes) std::printf("note: %s\n", n.c_str());
      if (written.size() > 1)
        std::printf("wrote %s … %s (%zu conformers)\n", written.front().c_str(), written.back().c_str(), written.size());
      else
        std::printf("wrote %s\n", written.front().c_str());
      return 0;
    } catch (const std::exception& e) {
      std::fprintf(stderr, "caps build: %s\n", e.what());
      keep_cli_error(argc, argv, e);
      return 1;
    }
  }
  if (cmd == "martini") {
    // Martini 3 small molecules: every molecule of an all-atom structure that is one of the model's molecules becomes its
    // beads with the explicit topology; water and ions are left out, other molecules too (listed)
    try {
      if (pos.empty()) return usage();
      if (!o.count("-o") || o["-o"].empty()) {
        std::fprintf(stderr, "caps martini: no output file: give -o CG.data\n");
        return 2;
      }
      if (o.count("--centre") && o["--centre"] != "geometric" && o["--centre"] != "mass")
        throw std::invalid_argument("--centre " + o["--centre"] + ": geometric or mass");
      std::string mdir;
      for (const std::string& root : caps_roots())
        if (!root.empty() && exists_here(root + "/data/martini/martini3-small-molecules.json")) { mdir = root + "/data/martini"; break; }
      if (mdir.empty()) throw std::runtime_error("data/martini/martini3-small-molecules.json not found (set CAPS_HOME)");
      if (o.count("--martini") && o["--martini"] != "3") throw std::runtime_error("only Martini 3 small molecules are mapped (--martini 3)");
      const System aa = load(pos[0], o);
      const bool geometric = !(o.count("--centre") && o["--centre"] == "mass");   // geometry (default) or mass
      std::vector<std::string> notes;
      const System cg = martini3_all_atom(aa, mdir + "/martini3-small-molecules.json", 1e6, &notes, geometric, false);
      write_structure_file(cg, o["-o"]);   // by its extension: LAMMPS data for .data, .gro, .pdb ... in their own format
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
      keep_cli_error(argc, argv, e);
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
      keep_cli_error(argc, argv, e);
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
        if (!jf) throw std::runtime_error("cannot open " + pos[1]);
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
      keep_cli_error(argc, argv, e);
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
          if (!r.empty() && exists_here(r + "/data/forcefields/catalogue.json")) { dir = r + "/data/forcefields"; break; }
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
        // UFF is built in (no library file): its description, as the C API and the Studio list it
        const FFDef ff = is_uff(pos[1]) ? uff_definition() : load_forcefield(ff_path(pos[1]));
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
        const bool uff = is_uff(o["--ff"]);
        if (uff && o.count("--typing"))
          throw std::invalid_argument("--typing: UFF types atoms from element, hybridisation and oxidation state, not from SMARTS rules");
        FFDef ff;
        TypingResult r;
        if (uff) {   // the built-in UFF typer (as ff apply --ff uff and the C API type)
          ff.name = "UFF";
          if (cell_fit(s).fits) {
            r.types = uff_types(s, &r.why);
          } else {   // a cell too small for its bonds: typed on the copies that fit, each atom given its first copy's type
            std::array<int, 3> n{1, 1, 1};
            try { n = copies_to_fit(s); } catch (const std::invalid_argument&) { throw std::invalid_argument(cell_fit_message(s)); }
            r.types = uff_types(supercell(s, n[0], n[1], n[2]), &r.why);
            r.types.resize(s.atoms.size());
            r.why.resize(s.atoms.size());
            r.notes.push_back("this cell is too small for its bonds (its bond list would give an atom the wrong neighbours), so it was typed as " +
                              std::to_string(n[0]) + "×" + std::to_string(n[1]) + "×" + std::to_string(n[2]) +
                              " copies and each atom given the type of its first copy");
          }
          r.rule.assign(s.atoms.size(), -1);
          r.candidates.assign(s.atoms.size(), {});
          for (size_t i = 0; i < r.types.size(); ++i)
            if (r.types[i].empty()) ++r.untyped;
            else r.candidates[i] = {r.types[i]};
        } else {
          ff = load_forcefield(ff_path(o["--ff"]));
          if (o.count("--typing"))
            for (std::stringstream ts(o["--typing"]); ts.good();) {   // several files: comma-separated, later ones on top
              std::string f;
              std::getline(ts, f, ',');
              if (!f.empty()) load_typing(ff, f);
            }
          std::string ch = "types";   // typing only: no charges
          if (const std::string ua = prepare_for_forcefield(s, ff, ch); !ua.empty()) std::printf("note: %s\n", ua.c_str());
          r = assign_types_fitted(s, ff);   // a cell too small for its bonds: typed on the copies that fit (cell_fit.hpp)
        }
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
        size_t cell_atoms = 0;   // with --supercell copies: the cell's own atoms (a --types file names them, each copy takes them)
        if (o.count("--supercell")) {   // auto: the copies this cell needs for its bonds (cell_fit.hpp); or na,nb,nc
          if (!s.cell.valid()) throw std::invalid_argument("--supercell: " + pos[1] + " has no cell");
          std::array<int, 3> n{1, 1, 1};
          std::string v = o["--supercell"];
          if (v == "auto") {
            n = copies_to_fit(s);
          } else {
            for (auto& c : v) if (c == 'x' || c == 'X') c = ',';
            if (std::sscanf(v.c_str(), "%d,%d,%d", &n[0], &n[1], &n[2]) != 3 || n[0] < 1 || n[1] < 1 || n[2] < 1)
              throw std::invalid_argument("--supercell: auto, or na,nb,nc (whole numbers of at least 1)");
          }
          if (n != std::array<int, 3>{1, 1, 1}) {
            if (o["--ff"] == "file") throw std::invalid_argument("--supercell: a file's own force field describes its own atoms only; make the copies with --ff NAME");
            cell_atoms = s.atoms.size();
            s = supercell(s, n[0], n[1], n[2]);
            std::printf("note: %d×%d×%d copies of the cell, %zu atoms\n", n[0], n[1], n[2], s.atoms.size());
          }
        }
        // a cell too small for its bonds gives a file whose bonds, angles and exclusions miss their image partners (a
        // polyethylene cell typed as an alkene): refused, not written with a warning
        if (const std::string small = cell_fit_message(s, o["--ff"] == "file" && s.forcefield ? s.forcefield.get() : nullptr); !small.empty())
          throw std::runtime_error(small + (o["--ff"] == "file" ? "" : " (caps ff apply --supercell auto makes those copies and writes them)"));
        const bool uff = is_uff(o["--ff"]);
        ForceField f;
        ParamReport rep;
        double cutoff = 10.0;
        const bool own = o["--ff"] == "file";   // the force field the file carries (an AMBER prmtop)
        MlFill ml_fill;                          // --complete ml: what CAPS ML fitted
        std::string ml_model;
        bool ml_ran = false;
        if ((own || uff) && completion_asked(o)) throw std::runtime_error("--complete ml fills a library force field's gaps (--ff NAME or FF.json), not UFF's or a file's own");
        if (!own) warn_keep_without_charges(o, s);
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
            if (r.untyped) throw std::runtime_error(untyped_message(ff, s, r.untyped, &r.types) + " (give them with --types)");
            std::printf("typed %zu atoms automatically (%d ambiguous; caps ff type --explain shows why)\n", s.atoms.size(), r.ambiguous);
            types = r.types;
          }
          if (auto_type) {
          } else if (o.count("--types")) {
            // one line per atom, or "index type" lines (1-based) to change only some atoms
            for (const auto& a : s.atoms) types.push_back(a.name);
            std::ifstream tf(o["--types"]);
            if (!tf) throw std::runtime_error("cannot open " + o["--types"]);
            const size_t named = cell_atoms ? cell_atoms : types.size();   // with --supercell: the cell's atoms, as caps ff type writes them
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
                if (i < 1 || i > named) throw std::runtime_error("atom " + w[0] + " out of range in " + o["--types"]);
                types[i - 1] = w[1];
              } else {
                if (k >= named) throw std::runtime_error("more types than atoms in " + o["--types"]);
                types[k++] = w[0];
              }
            }
            // one-per-line entries must cover every atom: the rest would keep their file names as types
            if (k > 0 && k < named)
              throw std::runtime_error(o["--types"] + " has " + std::to_string(k) + " types for " + std::to_string(named) +
                                       " atoms (give one per atom, or 'index type' lines to change only some)");
            for (size_t i = named; i < types.size(); ++i) types[i] = types[i % named];   // each copy: its cell atom's type
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
          ls.restraints = cli_restraints(s, o, true);
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
          for (const auto& n : write_gromacs(s, f, eo, o["--gromacs"], cli_restraints(s, o, true))) std::printf("gromacs: %s\n", n.c_str());
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
      keep_cli_error(argc, argv, e);
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
      const size_t n_read = s0.atoms.size();
      const std::map<std::string, std::string> o_read = o;   // the options before cli_forcefield records --charges keep
      auto ff = std::make_shared<ForceField>(cli_forcefield(s0, o));
      // A shell model (or another field that changes the atoms: shells added, united atoms) describes the structure
      // cli_forcefield prepared, not the frames as read: every frame handed to the solver gets the same preparation.
      auto prepared = [&](System& s) {
        if (s.atoms.size() == ff->atom_type.size()) return;
        std::map<std::string, std::string> oq = o_read;
        (void)cli_forcefield(s, oq, true);
        if (s.atoms.size() != ff->atom_type.size())
          throw std::invalid_argument("the force field is for " + std::to_string(ff->atom_type.size()) + " atoms; the frame has " +
                                      std::to_string(s.atoms.size()) + " after preparation");
      };
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
        prepared(s);
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
        if (ff->atom_type.size() != n_read)
          throw std::invalid_argument("the fluctuation method reads the frames as they are, and " + std::to_string(n_read) + " atoms per frame do not match the " +
                                      std::to_string(ff->atom_type.size()) + " the force field describes (a shell model adds shells the frames do not hold); "
                                      "use --method strain or --method fluct-run");
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
          prepared(s);
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
      keep_cli_error(argc, argv, e);
      return 1;
    }
  }
  if (cmd == "tensile") {
    // caps tensile DATA -o OUT.data [--axis x|y|z|xyz] [--rate 1e-3 (1/ps)] [--strain 0.2] [--temp 300] [--dt 1] [--fixed-lateral]
    //              [--equilibrate 20 (ps of NPT before the pull; 0 none; --eq-ps the same)] [--pressure 1] [--tau-p 1000] [--fit 0.02]
    //              [--seed 1] [--new-velocities] [--ff ...] [--dump T.lammpstrj --every N] [--checkpoint-every N] [--resume-curve FILE]
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
      if (o.count("--eq-ps")) to.equilibrate_ps = std::stod(o["--eq-ps"]);   // the name caps pull uses
      if (o.count("--equilibrate")) to.equilibrate_ps = std::stod(o["--equilibrate"]);
      if (to.equilibrate_ps < 0) throw std::invalid_argument("--equilibrate: ps of NPT before the pull, 0 or more");
      if (to.dt <= 0) throw std::invalid_argument("--dt: a time step in fs, above 0");
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
      if (!quiet && to.equilibrate_ps > 0 && !to.resume) {   // the prelude prints nothing while it runs: say what it is
        const long long neq = std::llround(to.equilibrate_ps * 1000 / to.dt);
        std::printf("NPT %g ps (%lld steps of %g fs%s) before the pull, each axis at %g atm; --equilibrate 0 skips it\n", to.equilibrate_ps, neq, to.dt,
                    to.axis == 3 ? ", once per axis" : "", to.pressure);
        std::fflush(stdout);
      }
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
      keep_cli_error(argc, argv, e);
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
      keep_cli_error(argc, argv, e);
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
            if (!root.empty() && exists_here(root + "/data/exchange/library.json")) return load_exchange_library(root + "/data/exchange/library.json");
          throw std::runtime_error("the exchange library (data/exchange/library.json) was not found");
        };
        if (sub == "list") {
          const auto lib = exchange_lib();
          // columns as wide as their longest entry: a mechanism such as "associative or dissociative" is printed whole
          int iw = 2, mw = 9, vw = 4, sw = 6;
          for (const auto& c : lib) {
            iw = std::max<int>(iw, int(c.id.size()));
            mw = std::max<int>(mw, int(c.mechanism.size()));
            vw = std::max<int>(vw, int(c.move.size()));
            sw = std::max<int>(sw, int(c.scheme.empty() ? 10 : c.scheme.size()));
          }
          std::printf("%-*s %-*s %-*s %-*s %s\n", iw, "id", mw, "mechanism", vw, "move", sw, "scheme", "default barrier (kcal/mol) and sources");
          for (const auto& c : lib) {
            char bar[64] = "none";
            if (std::isfinite(c.default_kcal_per_mol)) std::snprintf(bar, sizeof bar, "%.4g", c.default_kcal_per_mol);
            std::printf("%-*s %-*s %-*s %-*s %s, %zu value%s from the literature\n", iw, c.id.c_str(), mw, c.mechanism.c_str(), vw, c.move.c_str(), sw,
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
        // a chemistry's reaction template: its built-in network template, else from its scheme in data/reactions/library.json
        auto chemistry_template = [&](const ExchangeChemistry& c) {
          if (!c.template_name.empty()) return parse_templates(builtin_template(c.template_name)).at(0);
          if (c.scheme.empty()) throw std::invalid_argument(c.id + " has no reaction scheme yet (see caps exchange show " + c.id + ")");
          for (const std::string& root : caps_roots()) {
            const std::string path = root + "/data/reactions/library.json";
            if (root.empty() || !exists_here(path)) continue;
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
        KgOptions ko;
        if (o.count("--k-theta")) ko.k_theta = std::stod(o["--k-theta"]);
        const bool kg = o.count("--kg") > 0;
        const int axis = o.count("--axis") ? axis_of(o["--axis"]) : 2;
        const double plane = o.count("--plane") ? std::stod(o["--plane"]) : 0.5;
        if (axis > 2) throw std::invalid_argument("--axis takes x, y or z");
        const char* ax = axis == 0 ? "x" : axis == 1 ? "y" : "z";
        if (sub == "cut" || sub == "gap" || sub == "crack" || sub == "fragment") {
          // the structures a healing test starts from (caps/exchange_protocols.hpp)
          if (pos.size() < 2 || !o.count("-o")) return usage();
          System s = load(pos[1], o);
          if (sub == "cut" || sub == "gap") {
            std::printf("%zu bonds cut across the plane at %.4g of the cell along %s\n", cut_bonds_across(s, axis, plane), plane, ax);
            if (sub == "gap") {
              if (!o.count("--gap")) throw std::invalid_argument("gap needs --gap (Å between the two faces)");
              open_gap(s, axis, plane, std::stod(o["--gap"]));
              std::printf("the upper side moved %s A along %s, the cell lengthened by as much\n", o["--gap"].c_str(), ax);
            }
          } else if (sub == "crack") {
            if (!o.count("--radius")) throw std::invalid_argument("crack needs --radius (Å)");
            Vec3 c = s.cell.valid() ? s.cell.to_cartesian({0.5, 0.5, 0.5}) : Vec3{0, 0, 0};
            if (o.count("--centre")) {
              std::stringstream ss(o["--centre"]);
              std::string v;
              for (int k = 0; k < 3; ++k) {
                if (!std::getline(ss, v, ',')) throw std::invalid_argument("--centre takes x,y,z in Å");
                c[k] = std::stod(v);
              }
            }
            const size_t n = remove_sphere(s, c, std::stod(o["--radius"]));
            std::printf("%zu atoms removed within %s A of (%.4g, %.4g, %.4g)\n", n, o["--radius"].c_str(), c[0], c[1], c[2]);
          } else {
            std::vector<int> axes;
            for (char a : o.count("--axes") ? o["--axes"] : std::string("xyz")) axes.push_back(axis_of(std::string(1, a)));
            std::printf("%zu bonds cut across the mid-planes\n", fragment(s, axes));
          }
          write_structure_file(s, o["-o"]);
          std::printf("wrote %s\n", o["-o"].c_str());
          return 0;
        }
        // the exchange moves of run, step-strain, tv and heal-run
        auto exchange_moves = [&](const System& s) {
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
          // relaxed moves: all-atom templates by default (a new bond made at the capture distance), else when asked
          xo.relax_steps = o.count("--relax") ? std::stoi(o["--relax"]) : xo.kind == ExchangeKind::Template ? 1000 : 0;
          if (o.count("--relax-depth")) xo.relax_depth = std::stoi(o["--relax-depth"]);
          // new combinations of types: the run's own assignment types each change's neighbourhood, with --overlay and
          // --complete ml filling the terms the force field lacks. One definition serves the whole run, so what CAPS ML
          // fits joins it and is fitted once.
          if (!kg && o.count("--ff")) {
            std::shared_ptr<FFDef> def;
            if (!is_uff(o["--ff"])) {
              def = std::make_shared<FFDef>(load_forcefield(ff_path(o["--ff"])));
              if (o.count("--overlay")) merge_overlays(*def, o["--overlay"]);
              if (o.count("--typing")) load_typing(*def, o["--typing"]);
              if (needs_prepare(*def) || def->typing.empty()) def = nullptr;   // united atoms, shells, types by name: as the run assigned
            }
            if (!def)
              xo.parameterize = [&o](const System& f) {
                System y = f;
                return cli_forcefield(y, o, true);
              };
            else {
              auto given_up = std::make_shared<std::set<std::string>>();   // gaps CAPS ML could not fit: not tried again
              xo.parameterize = [&o, def, given_up](const System& f) {
                System y = f;
                const TypingResult tr = assign_types(y, *def);
                if (tr.untyped) throw std::runtime_error(untyped_message(*def, y, tr.untyped, &tr.types));
                auto assign = [&](ParamReport& rep) {
                  std::vector<std::string> said;
                  return parameterize_charges(y, *def, tr.types, file_charges(y) ? "keep" : "types", ff_path(o["--ff"]), &rep, true, &said);
                };
                ParamReport rep;
                ForceField ff = assign(rep);
                if (rep.missing.empty()) return ff;
                const bool fresh = std::any_of(rep.missing.begin(), rep.missing.end(), [&](const std::string& m) { return !given_up->count(m); });
                if (fresh && completion_asked(o)) {
                  // CAPS ML scans from a clean geometry of the same bond graph: the run's positions hold the new bond at
                  // the capture distance among crowded neighbours, and the scan's relaxations would clash
                  System z = y;
                  try {
                    z = embed_system(y);
                  } catch (const std::exception&) {
                    z = y;   // no embedding: the run's own geometry
                  }
                  const MlFill fill = cli_complete_ml(z, *def, tr.types, o, ff_path(o["--ff"]), true);
                  for (const auto& term : fill.fitted_terms) std::printf("CAPS ML, for an exchange: fitted %s\n", term.c_str());
                  for (const auto& term : fill.left) std::printf("CAPS ML, for an exchange: could not fit %s\n", term.c_str());
                  std::fflush(stdout);
                  rep = ParamReport{};
                  ff = assign(rep);
                  if (rep.missing.empty()) return ff;
                }
                for (const auto& m : rep.missing) given_up->insert(m);
                throw std::runtime_error(def->name + " has no parameters for " + rep.missing.front() +
                                         (rep.missing.size() > 1 ? " and " + std::to_string(rep.missing.size() - 1) + " more" : std::string()) +
                                         (completion_asked(o) ? "" : " (--complete ml fits them with CAPS ML)"));
              };
            }
          }
          if (o.count("--valence")) {   // bonds a site can hold: N for both groups, or N,M for A and B (a bead in a chain holds two)
            const std::string v = o["--valence"];
            const size_t c = v.find(',');
            xo.valence_a = std::stoi(v.substr(0, c));
            xo.valence_b = c == std::string::npos ? xo.valence_a : std::stoi(v.substr(c + 1));
          }
          if (o.count("--attempts")) xo.attempts = std::stoi(o["--attempts"]);
          if (o.count("--seed")) xo.seed = std::stoull(o["--seed"]);
          return xo;
        };
        auto energy_options = [&](EnergyOptions& e) {
          if (o.count("--cutoff")) e.cutoff = std::stod(o["--cutoff"]);
          electrostatics(e, o);
          precision(e, o);
          if (o.count("--threads")) e.threads = std::stoi(o["--threads"]);
        };
        const int64_t every = o.count("--every") ? std::stoll(o["--every"]) : 10;
        auto print_stats = [&](ExchangeKind kind, const ExchangeStats& st) {
          std::printf("%s moves every %lld steps: %lld tried, %lld accepted, %lld without a partner in reach, %lld the force field could not describe; "
                      "%.6g kcal/mol added by the accepted moves%s\n",
                      to_string(kind), (long long)every, (long long)st.attempts, (long long)st.accepted, (long long)st.no_candidate,
                      (long long)st.unparameterised, st.energy, st.local ? " (local energies)" : " (whole-system energies)");
          if (!st.note.empty()) std::printf("note: %s\n", st.note.c_str());
        };
        // --dt, else the force field's own step (the Kremer-Grest model's 0.01 tau, a coarse-grained model's), else the default
        auto time_step = [&](const ForceField& f, double fallback) {
          double dt = fallback;
          if (o.count("--dt")) dt = std::stod(o["--dt"]);
          else if (f.native_timestep > 0) {
            dt = f.native_timestep;
            if (!o.count("--quiet")) std::printf("time step: %.4g fs (the force field's own; --dt for another)\n", dt);
          }
          if (!(dt > 0) || !std::isfinite(dt)) throw std::invalid_argument("--dt: a time step in fs, above 0");
          return dt;
        };
        auto save_after = [&](System& s, const EnergyOptions& e) {
          if (!o.count("-o")) return;
          const ForceField after = kg ? kremer_grest_forcefield(s, ko) : cli_forcefield(s, o, true);
          save_structure(s, after, e, o["-o"]);
          std::printf("wrote %s\n", o["-o"].c_str());
        };
        if (sub == "run") {
          // CAPS's own MD with exchange moves (caps/exchange_moves.hpp)
          if (pos.size() < 2 || !o.count("-o")) return usage();
          System s = load(pos[1], o);
          ForceField ff = kg ? kremer_grest_forcefield(s, ko) : cli_forcefield(s, o);
          const ExchangeOptions xo = exchange_moves(s);
          DynamicsOptions d;
          d.field = std::make_shared<ForceField>(ff);
          d.steps = o.count("--steps") ? std::stoll(o["--steps"]) : 10000;
          d.temperature = o.count("--temp") ? std::stod(o["--temp"]) : 300;
          d.dt = time_step(ff, d.dt);
          if (o.count("--seed")) d.seed = std::stoull(o["--seed"]);
          energy_options(d.energy);
          d.topology_every = every;
          auto st = std::make_shared<ExchangeStats>();
          d.topology = exchange_callback(xo, d.energy, st);
          DynamicsReport r;
          run_dynamics(s, d, &r);
          print_stats(xo.kind, *st);
          for (const auto& n : r.notes) std::printf("%s\n", n.c_str());
          save_after(s, d.energy);
          return 0;
        }
        if (sub == "step-strain" || sub == "tv" || sub == "heal-run") {
          // the protocols (caps/exchange_protocols.hpp): MD with the exchange moves of run (none with --no-exchange)
          if (pos.size() < 2) return usage();
          System s = load(pos[1], o);
          ProtocolOptions po;
          po.field = std::make_shared<const ForceField>(kg ? kremer_grest_forcefield(s, ko) : cli_forcefield(s, o));
          energy_options(po.energy);
          po.temperature = o.count("--temp") ? std::stod(o["--temp"]) : 300;
          po.dt = time_step(*po.field, po.dt);
          if (o.count("--pressure")) po.pressure = std::stod(o["--pressure"]);
          if (o.count("--seed")) po.seed = std::stoull(o["--seed"]);
          po.exchange = !o.count("--no-exchange");
          if (po.exchange) po.moves = exchange_moves(s);
          po.exchange_every = every;
          auto need = [&](const char* flag) {
            if (!o.count(flag)) throw std::invalid_argument(sub + " needs " + flag);
            return std::stod(o[flag]);
          };
          ExchangeStats stats;
          const double kNaN = std::numeric_limits<double>::quiet_NaN();
          if (sub == "step-strain") {
            const double strain = need("--strain"), ps = need("--ps");
            const auto r = step_strain(s, po, axis, strain, ps, o.count("--sample") ? std::stoi(o["--sample"]) : 10);
            stats = r.stats;
            Property p{"step_strain", "Stress relaxation after a step strain", "ps", "", kNaN, kNaN, {}, {}, {}};
            p.method = std::string("volume-keeping step strain ") + o["--strain"] + " along " + ax + ", then MD at " + std::to_string(int(po.temperature)) +
                       " K; σ = σ_" + ax + " − (σ_other + σ_other)/2 from the pressure tensor; value: the KWW mean time τ Γ(1/β)/β";
            p.value = r.kww_mean;
            p.extra["KWW tau (ps)"] = r.kww_tau;
            p.extra["KWW beta"] = r.kww_beta;
            if (!r.stress.empty()) p.extra["stress at the start (MPa)"] = r.stress.front(), p.extra["stress at the end (MPa)"] = r.stress.back();
            p.series.push_back(Series{"σ(t)", "t (ps)", "stress (MPa)", r.t_ps, r.stress});
            if (!std::isfinite(r.kww_mean)) p.notes.push_back("the stress did not decay enough to fit; hold longer");
            props.push_back(p);
          } else if (sub == "tv") {
            const double from = need("--from"), to = need("--to"), step = need("--step"), ps = need("--ps");
            const auto r = tv_scan(s, po, from, to, step, ps);
            Property p{"tv", "Specific volume against temperature", "K", "", kNaN, kNaN, {}, {}, {}};
            p.method = "NPT (C-rescale) at each temperature with exchange; density from the second half of each; value: the hinge of a two-line fit "
                       "(the topology freezing temperature Tv when exchange sets it, else Tg)";
            if (r.fit.ok) {
              p.value = r.fit.tg, p.error = r.fit.tg_err;
              p.extra["expansion below (1/K)"] = r.fit.alpha_low, p.extra["expansion above (1/K)"] = r.fit.alpha_high;
            } else p.notes.push_back(r.fit.note.empty() ? "no hinge found" : r.fit.note);
            p.series.push_back(Series{"v(T)", "T (K)", "specific volume (cm³/g)", r.T, r.specific_volume});
            p.series.push_back(Series{"exchange rate", "T (K)", "accepted moves per ps", r.T, r.accepted_per_ps});
            props.push_back(p);
            if (r.rate.ok) {
              auto a = arrhenius_properties(r.rate, "time between accepted moves");
              props.insert(props.end(), a.begin(), a.end());
            }
          } else {
            const double ps = need("--ps");
            const System start = o.count("--start") ? load(o["--start"], o) : s;
            if (start.atoms.size() != s.atoms.size()) throw std::invalid_argument("--start must hold the same atoms");
            po.npt = o.count("--npt") > 0;
            const auto r = heal_run(s, po, sides_from_positions(start, axis, plane), axis, plane, ps, o.count("--chunks") ? std::stoi(o["--chunks"]) : 10);
            stats = r.stats;
            Property p{"heal_run", "Bonds bridging a cut during exchange", "bonds", "", kNaN, kNaN, {}, {}, {}};
            char where[64];
            std::snprintf(where, sizeof where, "%.4g along %s", plane, ax);
            p.method = std::string("MD with exchange") + (po.npt ? " at constant pressure" : "") + "; bridges: bonds joining atoms that started on the two sides of the plane at " +
                       where + " (sides from " + (o.count("--start") ? "--start" : "the input") +
                       "), across the cut plane; bonds joining the sides through the periodic cell face are not counted";
            p.value = r.bridges.empty() ? kNaN : r.bridges.back();
            p.series.push_back(Series{"bridges", "t (ps)", "bonds", r.t_ps, r.bridges});
            p.series.push_back(Series{"crossing the plane", "t (ps)", "bonds", r.t_ps, r.crossing});
            props.push_back(p);
          }
          if (po.exchange && sub != "tv") print_stats(po.moves.kind, stats);
          save_after(s, po.energy);
          cli_report(props, o);
          return 0;
        }
        if (sub == "bond-swap") {
          // LAMMPS fix bond/swap for a Kremer–Grest system (write_kg_bond_swap)
          if (pos.size() < 2 || !o.count("-o")) return usage();
          const System s = load(pos[1], o);
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
        throw std::invalid_argument("caps exchange arrhenius | modulus | memory | bridges | heal | list | show | run | bond-react | bond-swap | cut | gap | crack | fragment | step-strain | tv | heal-run");
      }
      cli_report(props, o);
      return 0;
    } catch (const std::exception& e) {
      std::fprintf(stderr, "caps exchange: %s\n", e.what());
      keep_cli_error(argc, argv, e);
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
        j["atoms"] = double(now.atoms.size());   // tg --fit takes the potential energy per atom
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
      keep_cli_error(argc, argv, e);
      return 1;
    }
  }
  if (cmd == "fields") {
    // caps fields --molecules 'A4B4:1; A8:0.5' [--chi AB=0.9,AC=0.5] [--grid 32,32,32] [--bond 1.1543] [--kappa 10] [--noise 75]
    //             [--dt 0.02] [--steps 2000] [--every 100] [--seed 1] [--vtk OUT.vtk] [--csv OUT.csv] [--json OUT.json]
    try {
      if (!o.count("--molecules")) throw std::invalid_argument("--molecules 'SEQUENCE:fraction; …' (e.g. 'A8B8:1')");
      DdftOptions d;
      {
        std::stringstream ss(o["--molecules"]);
        std::string item;
        while (std::getline(ss, item, ';')) {
          item.erase(0, item.find_first_not_of(" \t"));
          item.erase(item.find_last_not_of(" \t") + 1);
          if (item.empty()) continue;
          const auto c = item.find(':');
          const std::string seq = item.substr(0, c);
          d.molecules.push_back({seq, seq, c == std::string::npos ? 1.0 : std::stod(item.substr(c + 1))});
        }
      }
      if (o.count("--chi")) {
        std::stringstream ss(o["--chi"]);
        std::string item;
        while (std::getline(ss, item, ',')) {
          const auto e = item.find('=');
          if (e == std::string::npos) throw std::invalid_argument("--chi AB=0.9,AC=0.5");
          d.chi[item.substr(0, e)] = std::stod(item.substr(e + 1));
        }
      }
      if (o.count("--grid")) {
        std::string s = o["--grid"];
        for (auto& ch : s) if (ch == 'x' || ch == ':') ch = ',';
        std::stringstream ss(s);
        std::string x;
        for (int k = 0; k < 3 && std::getline(ss, x, ','); ++k) d.grid[size_t(k)] = std::stoi(x);
      }
      if (o.count("--bond")) d.bond = std::stod(o["--bond"]);
      if (o.count("--kappa")) d.kappa = std::stod(o["--kappa"]);
      if (o.count("--noise")) d.noise = std::stod(o["--noise"]);
      if (o.count("--dt")) d.dt = std::stod(o["--dt"]);
      if (o.count("--steps")) d.steps = std::stoi(o["--steps"]);
      if (o.count("--every")) d.frame_every = std::stoi(o["--every"]);
      if (o.count("--seed")) d.seed = std::stoull(o["--seed"]);
      const bool quiet = o.count("--quiet");
      if (!quiet) d.progress = [&](double f) { std::fprintf(stderr, "\rdensity fields %3.0f %%", 100 * f); return true; };
      const DdftResult r = ddft_run(d);
      if (!quiet) std::fprintf(stderr, "\r%-30s\r", "");
      for (const auto& n : r.notes) std::printf("%s\n", n.c_str());
      std::printf("%10s", "τ");
      for (const auto& ty : r.types) std::printf("  order(%s)", ty.c_str());
      std::printf("  %12s %10s\n", "χ energy", "spacing");
      for (size_t i = 0; i < r.times.size(); ++i) {
        std::printf("%10.3f", r.times[i]);
        for (double v : r.order[i]) std::printf("  %10.4g", v);
        std::printf("  %12.5g %10.3f\n", r.interaction[i], r.dominant_spacing[i]);
      }
      if (o.count("--csv")) {
        std::ofstream f(o["--csv"]);
        f << "tau";
        for (const auto& ty : r.types) f << ",order_" << ty;
        f << ",interaction,compressibility,spacing\n";
        for (size_t i = 0; i < r.times.size(); ++i) {
          f << r.times[i];
          for (double v : r.order[i]) f << "," << v;
          f << "," << r.interaction[i] << "," << r.compressibility[i] << "," << r.dominant_spacing[i] << "\n";
        }
        std::printf("wrote %s\n", o["--csv"].c_str());
      }
      if (o.count("--vtk") && !r.frames.empty()) {   // legacy VTK structured points, one scalar field per bead type (ParaView)
        std::ofstream f(o["--vtk"]);
        const size_t M = size_t(r.grid[0]) * size_t(r.grid[1]) * size_t(r.grid[2]);
        f << "# vtk DataFile Version 3.0\nCAPS density fields at tau " << r.times.back() << "\nASCII\nDATASET STRUCTURED_POINTS\n";
        f << "DIMENSIONS " << r.grid[0] << " " << r.grid[1] << " " << r.grid[2] << "\nORIGIN 0 0 0\nSPACING 1 1 1\nPOINT_DATA " << M << "\n";
        for (size_t t = 0; t < r.types.size(); ++t) {
          f << "SCALARS theta_" << r.types[t] << " float 1\nLOOKUP_TABLE default\n";
          for (size_t c = 0; c < M; ++c) f << r.frames.back()[t * M + c] << (c % 8 == 7 ? "\n" : " ");
          f << "\n";
        }
        std::printf("wrote %s\n", o["--vtk"].c_str());
      }
      return 0;
    } catch (const std::exception& e) {
      std::fprintf(stderr, "caps fields: %s\n", e.what());
      keep_cli_error(argc, argv, e);
      return 1;
    }
  }
  if (cmd == "tension") {
    try {
      if (pos.empty()) return usage();
      const std::string ax = o.count("--axis") ? o["--axis"] : "z";
      const int axis = ax == "x" ? 0 : ax == "y" ? 1 : ax == "z" ? 2 : -1;
      if (axis < 0) throw std::invalid_argument("--axis x, y or z");
      const auto samples = tension_samples_from_file(pos[0], axis, o.count("--cell") ? std::stod(o["--cell"]) : 0.0);
      const auto t = surface_tension(samples, o.count("--interfaces") ? std::stoi(o["--interfaces"]) : 2, o.count("--skip") ? std::stod(o["--skip"]) : 0.2,
                                     o.count("--blocks") ? std::stoi(o["--blocks"]) : 5);
      std::printf("γ = %.3f ± %.3f mN/m · %d samples in %d blocks · ⟨P_N⟩ %.2f atm, ⟨P_T⟩ %.2f atm, L %.3f Å (normal %s)\n", t.gamma, t.error, t.samples,
                  t.blocks, t.pn, t.pt, t.ln, ax.c_str());
      for (const auto& n : t.notes) std::printf("%s\n", n.c_str());
      return 0;
    } catch (const std::exception& e) {
      std::fprintf(stderr, "caps tension: %s\n", e.what());
      keep_cli_error(argc, argv, e);
      return 1;
    }
  }
  if (cmd == "table") {
    try {
      if (pos.empty()) return usage();
      std::ifstream in(pos[0]);
      if (!in) throw std::invalid_argument("cannot read " + pos[0]);
      std::vector<std::vector<std::string>> cells;
      for (std::string line; std::getline(in, line);) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        std::vector<std::string> f = split_csv_record(line);   // quoted fields ("Nylon 6,6") are one value, without their quotes
        for (auto& x : f) {
          const auto a = x.find_first_not_of(' '), b = x.find_last_not_of(' ');
          x = a == std::string::npos ? "" : x.substr(a, b - a + 1);
        }
        cells.push_back(f);
      }
      if (cells.size() < 2) throw std::invalid_argument(pos[0] + ": a header row and at least one row");
      const auto header = cells[0];
      cells.erase(cells.begin());
      auto column = [&](const std::string& name) {
        for (size_t k = 0; k < header.size(); ++k) if (header[k] == name) return k;
        throw std::invalid_argument("no column '" + name + "' (" + std::to_string(header.size()) + " columns in " + pos[0] + ")");
      };
      auto number = [](const std::string& t) {
        char* e = nullptr;
        const double v = std::strtod(t.c_str(), &e);
        return e && *e == 0 && !t.empty() ? v : NAN;
      };
      std::vector<size_t> order(cells.size());
      for (size_t r = 0; r < order.size(); ++r) order[r] = r;
      if (o.count("--sort")) {
        std::vector<std::vector<double>> keys;
        std::vector<bool> desc;
        std::string spec = o["--sort"], cur;
        for (char c : spec + ",") {
          if (c != ',') { cur += c; continue; }
          if (cur.empty()) continue;
          const bool d = cur[0] == '-';
          const size_t k = column(d ? cur.substr(1) : cur);
          std::vector<double> v;
          for (const auto& row : cells) v.push_back(k < row.size() ? number(row[k]) : NAN);
          keys.push_back(v), desc.push_back(d);
          cur.clear();
        }
        order = sort_rows(keys, desc);
        for (size_t k = 0; k < header.size(); ++k) std::printf("%s%s", k ? "\t" : "", header[k].c_str());
        std::printf("\n");
        for (size_t r : order) {
          for (size_t k = 0; k < header.size(); ++k) std::printf("%s%s", k ? "\t" : "", k < cells[r].size() ? cells[r][k].c_str() : "");
          std::printf("\n");
        }
      }
      if (o.count("--hist")) {
        const size_t k = column(o["--hist"]);
        std::vector<double> v;
        for (const auto& row : cells) v.push_back(k < row.size() ? number(row[k]) : NAN);
        std::string rule = o.count("--rule") ? o["--rule"] : "fd";
        int bins = 0;
        if (!rule.empty() && std::isdigit(static_cast<unsigned char>(rule[0]))) bins = std::stoi(rule), rule = "bins";
        const auto h = histogram(v, rule, bins);
        std::printf("# %s: %d values (%d left out), %s, width %.6g\n# from        to          count\n", o["--hist"].c_str(), h.n, h.left_out, h.rule.c_str(), h.width);
        for (size_t b = 0; b < h.counts.size(); ++b) std::printf("%-12.6g %-12.6g %d\n", h.edges[b], h.edges[b + 1], h.counts[b]);
      }
      if (o.count("--sdf")) {
        if (!o.count("--structure")) throw std::invalid_argument("--sdf needs --structure COLUMN (the column of structure files)");
        const size_t ks = column(o["--structure"]);
        const auto base = std::filesystem::path(pos[0]).parent_path();
        std::vector<SdfRow> rows;
        for (size_t r : order) {
          SdfRow row;
          const auto& c = cells[r];
          std::string f = ks < c.size() ? c[ks] : "";
          if (!f.empty() && std::filesystem::path(f).is_relative()) f = (base / f).string();
          row.structure = f, row.title = std::filesystem::path(f).stem().string();
          for (size_t k = 0; k < header.size(); ++k) if (k != ks) row.data.push_back({header[k], k < c.size() ? c[k] : ""});
          rows.push_back(row);
        }
        std::vector<std::string> notes;
        const int n = write_sdf_table(rows, o["--sdf"], &notes);
        for (const auto& x : notes) std::printf("%s\n", x.c_str());
        std::printf("wrote %s (%d records)\n", o["--sdf"].c_str(), n);
      }
      return 0;
    } catch (const std::exception& e) {
      std::fprintf(stderr, "caps table: %s\n", e.what());
      keep_cli_error(argc, argv, e);
      return 1;
    }
  }
  if (cmd == "cellrelax") {
    try {
      if (pos.empty() || !o.count("-o")) return usage();
      Trajectory t = open_file(pos[0], o.count("--topology") ? o["--topology"] : "");
      System s = t.frames() ? t.frame(t.frames() - 1) : t.topology;
      CellRelaxOptions co;
      if (o.count("--hold")) {
        std::string h = o["--hold"];
        for (auto& ch : h) if (ch == ',') ch = ' ';
        std::istringstream in(h);
        static const char* names[6] = {"a", "b", "c", "alpha", "beta", "gamma"};
        for (std::string w; in >> w;) {
          bool found = false;
          for (int k = 0; k < 6; ++k) if (w == names[k]) co.hold[size_t(k)] = true, found = true;
          if (!found) throw std::invalid_argument("--hold takes a, b, c, alpha, beta, gamma (comma-separated), not '" + w + "'");
        }
      }
      if (o.count("--pressure")) co.pressure = std::stod(o["--pressure"]);
      if (o.count("--ftol")) co.ftol = std::stod(o["--ftol"]);
      if (o.count("--steps")) co.max_steps = std::stoi(o["--steps"]);
      const bool quiet = o.count("--quiet");
      EnergyOptions eo;
      eo.cutoff = o.count("--cutoff") ? std::stod(o["--cutoff"]) : 10.0;
      electrostatics(eo, o);
      auto field_for = [&](const System& sup) { System x = sup; return cli_forcefield(x, o, true); };
      const System s_in = s;
      (void)cli_forcefield(s, o, quiet);   // the structure as the force field describes it (shells, united atoms, charges); says which field
      const auto r = relax_cell(s, field_for, eo, co);
      for (const auto& n : r.notes) std::printf("%s\n", n.c_str());
      const char* nm[6] = {"a", "b", "c", "alpha", "beta", "gamma"};
      std::printf("parameter      before        after   %s\n", "");
      for (int k = 0; k < 6; ++k)
        std::printf("%-8s %12.6f %12.6f%s\n", nm[k], r.before[size_t(k)], r.after[size_t(k)], co.hold[size_t(k)] ? "   held" : "");
      std::printf("volume %.4f → %.4f Å³ · energy %.6f → %.6f kcal/mol · pressure %.2f atm · %d steps%s\n", r.volume_before, r.volume_after,
                  r.energy_before, r.energy_after, r.pressure_after, r.steps, r.converged ? "" : " (not converged)");
      std::printf("stress GPa (xx yy zz yz xz xy): %.5f %.5f %.5f %.5f %.5f %.5f\n", r.stress[0], r.stress[1], r.stress[2], r.stress[3], r.stress[4], r.stress[5]);
      // a .data file keeps the field it was relaxed with (shells as core-shell types); any other format the atoms as read
      if (lammps_data_name(o["-o"])) write_structure_file(s, field_for(s), eo, o["-o"]);
      else write_structure_file(as_read(s_in, s), o["-o"]);
      std::printf("wrote %s\n", o["-o"].c_str());
      return r.converged ? 0 : 2;
    } catch (const std::exception& e) {
      std::fprintf(stderr, "caps cellrelax: %s\n", e.what());
      keep_cli_error(argc, argv, e);
      return 1;
    }
  }
  if (cmd == "path") {
    try {
      if (pos.size() < 2 || !o.count("-o")) return usage();
      Trajectory ta = open_file(pos[0], o.count("--topology") ? o["--topology"] : ""), tb = open_file(pos[1]);
      System a = ta.frames() ? ta.frame(ta.frames() - 1) : ta.topology;
      const System b0 = tb.frames() ? tb.frame(tb.frames() - 1) : tb.topology;
      AtomMatch m;
      if (o.count("--no-match")) {   // B's file order kept (a permutation matching would undo, as a methyl turned by 120°)
        if (b0.atoms.size() != a.atoms.size()) throw std::invalid_argument("A and B have different atom counts");
        for (size_t i = 0; i < a.atoms.size(); ++i) m.map.push_back(int(i));
      } else {
        m = match_atoms(a, b0);
      }
      for (const auto& n : m.notes) std::printf("%s\n", n.c_str());
      System b = reorder_atoms(b0, m.map, a);
      if (o.count("--relax-ends")) {   // both ends minimised with the force field first
        auto ffp = std::make_shared<ForceField>(cli_forcefield(a, o, true));
        RelaxOptions ro;
        ro.field = ffp, ro.ftol = 0.05, ro.max_iterations = 20000;
        relax(a, ro);
        relax(b, ro);
        std::printf("ends minimised with the force field to 0.05 kcal/mol/Å\n");
      }
      PathOptions po;
      if (o.count("--method")) po.method = o["--method"];
      if (o.count("--images")) po.images = std::stoi(o["--images"]);
      PathReport pr;
      auto path = interpolate_path(a, b, po, &pr);
      for (const auto& n : pr.notes) std::printf("%s\n", n.c_str());
      std::vector<double> energy;
      if (o.count("--neb")) {
        const bool quiet = o.count("--quiet");
        const ForceField ff = cli_forcefield(a, o, quiet);
        EnergyOptions eo;
        if (o.count("--cutoff")) eo.cutoff = std::stod(o["--cutoff"]);
        run_cutoff(eo, ff, o, quiet);
        electrostatics(eo, o);
        eo.coulomb = !o.count("--no-coulomb");
        NebOptions no;
        if (o.count("--spring")) no.spring = std::stod(o["--spring"]);
        if (o.count("--fmax")) no.fmax = std::stod(o["--fmax"]);
        if (o.count("--steps")) no.max_steps = std::stoi(o["--steps"]);
        no.climb = !o.count("--no-climb");
        const NebResult r = neb_forcefield(a, path, ff, eo, no);
        for (const auto& n : r.notes) std::printf("%s\n", n.c_str());
        std::printf("image   arc Å     E kcal/mol   ΔE from A\n");
        for (size_t k = 0; k < r.images.size(); ++k)
          std::printf("%5zu %8.3f %14.4f %11.4f%s\n", k, r.arc[k], r.energy[k], r.energy[k] - r.energy[0], int(k) == r.climbing ? "   climbing" : "");
        path = r.images;
        energy = r.energy;
      } else {
        std::printf("image   shortest distance Å\n");
        for (size_t k = 0; k < pr.min_distance.size(); ++k) std::printf("%5zu %10.4f\n", k, pr.min_distance[k]);
      }
      Trajectory out;
      out.topology = a;
      for (size_t k = 0; k < path.size(); ++k) {
        out.positions.push_back(path[k]);
        out.cells.push_back(a.cell);
        out.timesteps.push_back(int64_t(k));
      }
      write_trajectory(out, o["-o"]);
      std::printf("wrote %s (%zu images)\n", o["-o"].c_str(), path.size());
      return 0;
    } catch (const std::exception& e) {
      std::fprintf(stderr, "caps path: %s\n", e.what());
      keep_cli_error(argc, argv, e);
      return 1;
    }
  }
  if (cmd == "data") {
    try {
      if (pos.empty()) return usage();
      const std::string& in = pos[0];
      auto ends = [](const std::string& x, const char* e) { const std::string t(e); return x.size() >= t.size() && x.compare(x.size() - t.size(), t.size(), t) == 0; };
      if (ends(in, ".cube") || ends(in, ".CUBE")) {
        System atoms;
        const GridField g = read_cube(in, &atoms);
        double lo = 1e300, hi = -1e300, sum = 0;
        for (double v : g.values) lo = std::min(lo, v), hi = std::max(hi, v), sum += v;
        std::printf("Gaussian cube \"%s\": %d × %d × %d points, %zu atoms, volume %.3f Å³\nvalues %.6g … %.6g, mean %.6g\n", g.name.c_str(), g.n[0], g.n[1], g.n[2],
                    atoms.atoms.size(), g.cell.volume(), lo, hi, sum / double(g.values.size()));
        if (o.count("-o")) {
          write_grid(g, atoms, o["-o"]);
          std::printf("wrote %s\n", o["-o"].c_str());
        }
        return 0;
      }
      const MeasuredCurve c = read_measured(in);
      if (c.x.empty()) {   // a header with no readable data lines after it
        for (const auto& n : c.notes) std::printf("%s\n", n.c_str());
        throw std::runtime_error(in + ": no data points could be read" + (c.format.empty() ? "" : " (" + c.format + ")"));
      }
      std::printf("%s%s%s: %zu points, x %.6g … %.6g%s%s, y %s%s\n", c.format.c_str(), c.title.empty() ? "" : " · ", c.title.c_str(), c.x.size(), c.x.front(), c.x.back(),
                  c.xunits.empty() ? "" : " ", c.xunits.c_str(), c.yunits.empty() ? "as written" : c.yunits.c_str(), c.sigma.empty() ? "" : " with σ");
      for (const auto& n : c.notes) std::printf("%s\n", n.c_str());
      if (o.count("-o") && (ends(o["-o"], ".npy") || ends(o["-o"], ".NPY"))) {
        // NumPy .npy v1.0: float64, shape (N, 2) of x, y, or (N, 3) with sigma, row by row (as write_grid writes a grid)
        const size_t cols = c.sigma.empty() ? 2 : 3;
        std::ofstream f(o["-o"], std::ios::binary);
        if (!f) throw std::runtime_error("cannot write " + o["-o"]);
        std::string h = "{'descr': '<f8', 'fortran_order': False, 'shape': (" + std::to_string(c.x.size()) + ", " + std::to_string(cols) + "), }";
        const size_t total = 10 + h.size() + 1;
        h += std::string((64 - total % 64) % 64, ' ') + "\n";
        const unsigned char magic[8] = {0x93, 'N', 'U', 'M', 'P', 'Y', 1, 0};
        f.write(reinterpret_cast<const char*>(magic), 8);
        const uint16_t hl = uint16_t(h.size());
        const unsigned char l2[2] = {uint8_t(hl & 0xFF), uint8_t(hl >> 8)};
        f.write(reinterpret_cast<const char*>(l2), 2);
        f.write(h.data(), std::streamsize(h.size()));
        std::vector<double> rows;
        rows.reserve(c.x.size() * cols);
        for (size_t k = 0; k < c.x.size(); ++k) {
          rows.push_back(c.x[k]), rows.push_back(c.y[k]);
          if (cols == 3) rows.push_back(c.sigma[k]);
        }
        f.write(reinterpret_cast<const char*>(rows.data()), std::streamsize(rows.size() * sizeof(double)));
        if (!f) throw std::runtime_error("cannot write " + o["-o"]);
        std::printf("wrote %s\n", o["-o"].c_str());
      } else if (o.count("-o")) {
        {
          std::string low = std::filesystem::path(o["-o"]).extension().string();
          for (char& ch : low) ch = char(std::tolower(static_cast<unsigned char>(ch)));
          if (!low.empty() && low != ".csv" && low != ".txt" && low != ".dat")
            throw std::invalid_argument("-o " + o["-o"] + ": write a curve as .csv (or .txt, .dat) or .npy; .cube and .vtk are for a .cube grid");
        }
        std::ofstream f(o["-o"]);
        if (!f) throw std::runtime_error("cannot write " + o["-o"]);
        f << "# " << c.format << (c.title.empty() ? "" : ": " + c.title) << "\n# x" << (c.xunits.empty() ? "" : " (" + c.xunits + ")") << ", y" << (c.sigma.empty() ? "" : ", sigma") << "\n";
        char b[96];
        for (size_t k = 0; k < c.x.size(); ++k) {
          if (c.sigma.empty()) std::snprintf(b, sizeof b, "%.8g,%.8g\n", c.x[k], c.y[k]);
          else std::snprintf(b, sizeof b, "%.8g,%.8g,%.8g\n", c.x[k], c.y[k], c.sigma[k]);
          f << b;
        }
        std::printf("wrote %s\n", o["-o"].c_str());
      }
      if (o.count("--json")) {
        std::ofstream(o["--json"]) << measured_json(c).dump();
        std::printf("wrote %s\n", o["--json"].c_str());
      }
      return 0;
    } catch (const std::exception& e) {
      std::fprintf(stderr, "caps data: %s\n", e.what());
      keep_cli_error(argc, argv, e);
      return 1;
    }
  }
  if (cmd == "symmetry") {
    try {
      if (pos.empty()) return usage();
      Trajectory t = open_file(pos[0], o.count("--topology") ? o["--topology"] : "");
      System s = t.frames() ? t.frame(t.frames() - 1) : t.topology;
      if (o.count("--clear")) s.symmetry = {};
      if (o.count("--find")) s.symmetry = symmetry_of(s, o.count("--tolerance") ? std::stod(o["--tolerance"]) : 0.1);
      if (o.count("--cell")) {
        std::vector<double> v;
        std::stringstream ss(o["--cell"]);
        for (std::string x; std::getline(ss, x, ',');) v.push_back(std::stod(x));
        if (v.size() != 6) throw std::invalid_argument("--cell a,b,c,alpha,beta,gamma");
        set_cell_symmetric(s, v[0], v[1], v[2], v[3], v[4], v[5]);
        const double la = norm(s.cell.a), lb = norm(s.cell.b), lc = norm(s.cell.c);
        auto ang = [](const Vec3& x, const Vec3& y) { return std::acos(std::clamp(dot(x, y) / (norm(x) * norm(y)), -1.0, 1.0)) * 180 / M_PI; };
        std::printf("cell %.6f %.6f %.6f  %.4f %.4f %.4f\n", la, lb, lc, ang(s.cell.b, s.cell.c), ang(s.cell.a, s.cell.c), ang(s.cell.a, s.cell.b));
      }
      const auto& g = s.symmetry;
      if (g.ops.empty()) {
        std::printf("no space group kept (P 1); --find looks for one\n");
      } else {
        double worst = 0;
        const bool holds = symmetry_holds(s, g, 0.05, &worst);
        std::printf("%s (No. %d, %s) · %zu operations · %s (worst image %.2e Å)\n", g.name.c_str(), g.number, g.number > 0 ? crystal_system(g.number).c_str() : "?",
                    g.ops.size(), holds ? "holds" : "broken: written as P 1", worst);
        if (holds) std::printf("asymmetric unit: %zu of %zu atoms\n", asymmetric_unit(s, g, 0.05).size(), s.atoms.size());
      }
      if (o.count("-o")) {
        const std::string& out = o["-o"];
        if (out.size() < 4 || out.compare(out.size() - 4, 4, ".cif") != 0) throw std::invalid_argument("-o takes a .cif (other formats carry no space group)");
        write_cif(s, out);
        std::printf("wrote %s\n", out.c_str());
      }
      return 0;
    } catch (const std::exception& e) {
      std::fprintf(stderr, "caps symmetry: %s\n", e.what());
      keep_cli_error(argc, argv, e);
      return 1;
    }
  }
  if (cmd == "pores") {
    try {
      if (pos.empty()) return usage();
      Trajectory t = open_file(pos[0], o.count("--topology") ? o["--topology"] : "");
      const System s = t.frames() ? t.frame(t.frames() - 1) : t.topology;
      PoreNetworkOptions po;
      if (o.count("--probe")) po.probe = std::stod(o["--probe"]);
      if (o.count("--grid")) po.grid = std::stod(o["--grid"]);
      if (o.count("--samples")) po.samples = std::stoi(o["--samples"]);
      if (o.count("--radii") && o["--radii"] != "bondi") po.radii = free_volume_radii(s, o["--radii"]);
      const PoreNetworkReport r = pore_network(s, po);
      for (const auto& n : r.notes) std::printf("%s\n", n.c_str());
      std::printf("density %.4f g/cm³\n", r.density);
      std::printf("LCD %.3f Å   PLD %.3f Å (a %.3f, b %.3f, c %.3f)   Dif %.3f Å\n", r.lcd, r.pld, r.pld_axis[0], r.pld_axis[1], r.pld_axis[2], r.dif);
      std::printf("probe-centre volume: channels %.4f of the cell (%.4f cm³/g), pockets %.4f (%.4f cm³/g); geometric %.4f\n", r.av_fraction, r.av,
                  r.nav_fraction, r.nav, r.geometric_fraction);
      std::printf("surface: accessible %.2f Å² (%.2f m²/cm³, %.2f m²/g), non-accessible %.2f Å² (%.2f m²/cm³, %.2f m²/g)\n", r.asa, r.asa_volume,
                  r.asa_mass, r.nasa, r.nasa_volume, r.nasa_mass);
      std::printf("cluster   volume Å³   largest Å   kind      percolates\n");
      for (size_t k = 0; k < r.clusters.size() && k < 30; ++k) {
        const auto& c = r.clusters[k];
        std::printf("%6zu  %10.2f  %9.3f   %-8s  %c%c%c\n", k + 1, c.volume, c.largest, c.channel ? "channel" : "pocket", c.percolates[0] ? 'a' : '-',
                    c.percolates[1] ? 'b' : '-', c.percolates[2] ? 'c' : '-');
      }
      if (r.clusters.size() > 30) std::printf("… %zu more\n", r.clusters.size() - 30);
      if (o.count("--json")) {
        std::ofstream(o["--json"]) << pore_json(r).dump();
        std::printf("wrote %s\n", o["--json"].c_str());
      }
      if (o.count("--cube")) {
        write_grid(r.field, s, o["--cube"]);
        std::printf("wrote %s\n", o["--cube"].c_str());
      }
      return 0;
    } catch (const std::exception& e) {
      std::fprintf(stderr, "caps pores: %s\n", e.what());
      keep_cli_error(argc, argv, e);
      return 1;
    }
  }
  if (cmd == "qeq") {
    try {
      if (pos.empty()) return usage();
      if (o.count("-o")) {   // checked before the solve: a PDB has no partial-charge column (its columns 79-80 hold an integer formal charge)
        const std::string ext = std::filesystem::path(o["-o"]).extension().string();
        if (ext == ".pdb" || ext == ".ent")
          throw std::invalid_argument("-o " + o["-o"] + ": PDB has no column for partial charges, so the QEq charges would be lost; write .mol2, .car or "
                                      ".data, or the charges alone with --chg");
        if (ext != ".data" && ext != ".mol2" && ext != ".car") throw std::invalid_argument("-o takes .data, .mol2 or .car (the charges alone: --chg)");
      }
      Trajectory t = open_file(pos[0], o.count("--topology") ? o["--topology"] : "");
      System s = t.frames() ? t.frame(t.frames() - 1) : t.topology;
      QEqOptions qo;
      if (o.count("--kernel")) qo.kernel = o["--kernel"];
      if (o.count("--ewald") && o.count("--no-ewald")) throw std::invalid_argument("--ewald and --no-ewald together");
      qo.ewald = o.count("--no-ewald") == 0;   // Ewald by default; used only for a cell periodic along a, b and c
      qo.ewald_auto = !o.count("--ewald") && !o.count("--no-ewald");   // by default not for a cell too large to sum quickly
      if (o.count("--accuracy")) qo.ewald_accuracy = std::stod(o["--accuracy"]);
      if (o.count("--lambda")) qo.lambda = std::stod(o["--lambda"]);
      qo.hydrogen = o.count("--hydrogen") > 0;
      if (o.count("--bounds")) qo.bounds = qeq_bounds_from_text(o["--bounds"]);
      if (o.count("--cutoff")) qo.cutoff = std::stod(o["--cutoff"]);
      if (o.count("--total")) qo.total_charge = std::stod(o["--total"]);
      QEqReport rep;
      const auto q = qeq_charges(s, qo, &rep);
      for (const auto& n : rep.notes) std::printf("%s\n", n.c_str());
      std::printf("equalised electronegativity %.4f eV\n", rep.mu);
      std::map<int, std::vector<double>> by;
      for (size_t i = 0; i < q.size(); ++i) by[s.atoms[i].element].push_back(q[i]);
      std::printf("element  atoms     mean q        min        max\n");
      for (const auto& [z, v] : by) {
        double m = 0;
        for (double x : v) m += x;
        std::printf("%-7s %6zu  %9.4f  %9.4f  %9.4f\n", element(z).symbol, v.size(), m / double(v.size()), *std::min_element(v.begin(), v.end()),
                    *std::max_element(v.begin(), v.end()));
      }
      const bool allow = o.count("--allow-unconverged") > 0;
      if (!rep.converged && !allow) {
        std::fprintf(stderr, "caps qeq: the charges did not converge (see the notes above); nothing written (--allow-unconverged writes them anyway)\n");
        return 1;
      }
      if (!rep.converged) std::fprintf(stderr, "caps qeq: the charges did not converge; written as asked (--allow-unconverged)\n");
      if (o.count("--chg")) {
        std::ofstream f(o["--chg"]);
        f << "# CAPS QEq charges (e), one per atom in order\n";
        char b[40];
        for (double x : q) {
          std::snprintf(b, sizeof b, "%.6f\n", x);
          f << b;
        }
        std::printf("wrote %s\n", o["--chg"].c_str());
      }
      if (o.count("-o")) {
        for (size_t i = 0; i < q.size(); ++i) s.atoms[i].charge = q[i];
        s.has_charges = true;
        const std::string& out = o["-o"];
        auto ends = [&](const char* x) { const std::string e(x); return out.size() >= e.size() && out.compare(out.size() - e.size(), e.size(), e) == 0; };
        if (ends(".data")) write_lammps_data(s, out);
        else if (ends(".mol2")) write_mol2(s, out);
        else if (ends(".car")) write_car(s, out);
        else throw std::invalid_argument("-o takes .data, .mol2 or .car (the charges alone: --chg)");
        std::printf("wrote %s\n", out.c_str());
      }
      return 0;
    } catch (const std::exception& e) {
      std::fprintf(stderr, "caps qeq: %s\n", e.what());
      keep_cli_error(argc, argv, e);
      return 1;
    }
  }
  if (cmd == "dendrimer") {
    try {
      if (!o.count("-o")) return usage();
      DendrimerSpec sp;
      if (o.count("--preset")) sp.preset = o["--preset"];
      if (o.count("--core")) sp.core = o["--core"];
      if (o.count("--unit")) sp.unit = o["--unit"];
      if (o.count("--end")) sp.end = o["--end"];
      if (o.count("--focal")) sp.focal = o["--focal"];
      if (o.count("--generation")) sp.generation = std::stoi(o["--generation"]);
      if (o.count("--shells")) sp.shells = std::stoi(o["--shells"]);
      if (o.count("--growth")) sp.growth = o["--growth"];
      if (o.count("--units")) sp.units = std::stoi(o["--units"]);
      if (o.count("--seed")) sp.seed = std::stoull(o["--seed"]);
      sp.relax = !o.count("--no-relax");
      DendrimerReport r;
      const System s = build_dendrimer(sp, &r);
      for (const auto& n : r.notes) std::printf("%s\n", n.c_str());
      std::printf("%d AB%d units in %d shells · %s · %.2f g/mol · %zu atoms%s\n", r.units, r.n, r.shells, r.formula.c_str(), r.mass, s.atoms.size(),
                  r.ends ? (" · " + std::to_string(r.ends) + " end groups").c_str() : "");
      std::printf("shell  units\n");
      for (size_t k = 0; k < r.counts.per_shell.size(); ++k) std::printf("%5zu %6d\n", k + 1, r.counts.per_shell[k]);
      std::printf("B sites used  units\n");
      for (size_t k = 0; k < r.counts.by_used.size(); ++k) std::printf("%12zu %6d\n", k, r.counts.by_used[k]);
      if (r.n == 2)
        std::printf("D %d · L %d · T %d · DB (Frey 2D/(2D+L)) %.4f · DB (Fréchet (D+T)/(D+T+L)) %.4f · large-N limit of this growth %.4f\n", r.counts.dendritic,
                    r.counts.linear, r.counts.terminal, r.counts.db_frey, r.counts.db_frechet, r.db_limit);
      write_structure_file(s, o["-o"]);
      std::printf("wrote %s\n", o["-o"].c_str());
      return 0;
    } catch (const std::exception& e) {
      std::fprintf(stderr, "caps dendrimer: %s\n", e.what());
      keep_cli_error(argc, argv, e);
      return 1;
    }
  }
  if (cmd == "dpd") {
    // caps dpd --molecules 'SEQ:COUNT[@REGION[:in|out]]; …' [--regions 'NAME SHAPE c=fx,fy,fz r=R in=R0 t=T axis=z prio=P; …']
    //          [--box L] [--chi AB=25,…] [--density 3] [--gamma 4.5] [--dt 0.04] [--steps N] [--equilibration N] [--every N]
    //          [--rc-angstrom 6.46] [--seed 1] [--spec SPEC.json (the C API's JSON; options given here win)] [-o OUT] [--json OUT]
    try {
      Json j = Json::object();
      if (o.count("--spec")) {
        std::ifstream f(o["--spec"]);
        if (!f) throw std::invalid_argument("cannot read " + o["--spec"]);
        std::stringstream b;
        b << f.rdbuf();
        j = Json::parse(b.str());
      }
      auto trim = [](std::string x) {
        x.erase(0, x.find_first_not_of(" \t"));
        x.erase(x.find_last_not_of(" \t") + 1);
        return x;
      };
      if (o.count("--regions")) {
        Json rs = Json::array();
        std::stringstream ss(o["--regions"]);
        for (std::string item; std::getline(ss, item, ';');) {
          item = trim(item);
          if (item.empty()) continue;
          std::stringstream ws(item);
          std::string name, shape, kv;
          ws >> name >> shape;
          if (shape.empty()) throw std::invalid_argument("--regions 'NAME SHAPE key=value …'");
          Json r = Json::object();
          r["name"] = name, r["shape"] = shape;
          while (ws >> kv) {
            const auto e = kv.find('=');
            if (e == std::string::npos) throw std::invalid_argument("--regions: " + kv + " is not key=value");
            const std::string k = kv.substr(0, e), v = kv.substr(e + 1);
            if (k == "c") {
              Json c = Json::array();
              std::stringstream cs(v);
              for (std::string x; std::getline(cs, x, ',');) c.push_back(std::stod(x));
              if (c.items().size() != 3) throw std::invalid_argument("--regions: c=fx,fy,fz");
              r["centre"] = c;
            } else if (k == "r") r["radius"] = std::stod(v);
            else if (k == "in") r["inner"] = std::stod(v);
            else if (k == "t") r["thickness"] = std::stod(v);
            else if (k == "axis") r["axis"] = v;
            else if (k == "prio") r["priority"] = std::stod(v);
            else throw std::invalid_argument("--regions: unknown key " + k + " (c, r, in, t, axis, prio)");
          }
          rs.push_back(r);
        }
        j["regions"] = rs;
      }
      if (o.count("--molecules")) {
        Json ms = Json::array();
        std::stringstream ss(o["--molecules"]);
        for (std::string item; std::getline(ss, item, ';');) {
          item = trim(item);
          if (item.empty()) continue;
          Json m = Json::object();
          const auto at = item.find('@');
          std::string head = item.substr(0, at);
          if (at != std::string::npos) {
            std::string where = item.substr(at + 1);
            const auto c = where.find(':');
            m["region"] = trim(where.substr(0, c));
            if (c != std::string::npos) m["orient"] = trim(where.substr(c + 1));
          }
          const auto c = head.find(':');
          if (c == std::string::npos) throw std::invalid_argument("--molecules 'SEQUENCE:COUNT[@REGION[:in|out]]; …'");
          m["name"] = trim(head.substr(0, c)), m["sequence"] = trim(head.substr(0, c)), m["count"] = std::stod(head.substr(c + 1));
          ms.push_back(m);
        }
        j["species"] = ms;
      }
      if (o.count("--chi")) {
        Json chi = Json::object();
        std::stringstream ss(o["--chi"]);
        for (std::string item; std::getline(ss, item, ',');) {
          const auto e = item.find('=');
          if (e == std::string::npos) throw std::invalid_argument("--chi AB=25,AC=5");
          chi[trim(item.substr(0, e))] = std::stod(item.substr(e + 1));
        }
        j["chi"] = chi;
      }
      const std::pair<const char*, const char*> nums[] = {{"--box", "box"}, {"--density", "density"}, {"--gamma", "gamma"}, {"--dt", "dt"}, {"--steps", "steps"},
                                                          {"--equilibration", "equilibration"}, {"--every", "frame_every"}, {"--rc-angstrom", "rc_angstrom"}, {"--seed", "seed"}};
      for (const auto& [flag, key] : nums)
        if (o.count(flag)) j[key] = std::stod(o[flag]);
      if (!j.has("species")) throw std::invalid_argument("--molecules 'SEQUENCE:COUNT; …' or --spec SPEC.json");
      DpdOptions d = dpd_options_from_json(j);
      const bool quiet = o.count("--quiet");
      if (!quiet) d.progress = [&](long step, double kT) { std::fprintf(stderr, "\rDPD step %ld  kT %.3f   ", step, kT); return true; };
      const DpdReport r = run_dpd(d);
      if (!quiet) std::fprintf(stderr, "\r%-40s\r", "");
      for (const auto& n : r.notes) std::printf("%s\n", n.c_str());
      std::printf("%d beads in %d molecules, box %.4g r_c\n", r.beads, r.molecules, r.box);
      if (d.steps > 0)
        std::printf("kT %.4f ± %.4f   pressure %.4f ± %.4f   ψ %.4f   spacing %.4g r_c\n", r.kT, r.kT_error, r.pressure, r.pressure_error, r.order, r.spacing);
      if (o.count("-o")) {
        write_trajectory(r.frames, o["-o"]);
        std::printf("wrote %s (%zu frames)\n", o["-o"].c_str(), r.frames.frames());
      }
      if (o.count("--json")) {
        Json out = Json::object();
        out["beads"] = double(r.beads), out["molecules"] = double(r.molecules), out["box"] = r.box, out["kT"] = r.kT, out["pressure"] = r.pressure;
        out["order"] = r.order, out["spacing"] = r.spacing, out["q_peak"] = r.q_peak;
        Json notes = Json::array();
        for (const auto& n : r.notes) notes.push_back(n);
        out["notes"] = notes;
        std::ofstream(o["--json"]) << out.dump();
        std::printf("wrote %s\n", o["--json"].c_str());
      }
      return 0;
    } catch (const std::exception& e) {
      std::fprintf(stderr, "caps dpd: %s\n", e.what());
      keep_cli_error(argc, argv, e);
      return 1;
    }
  }
  if (cmd == "qha") {
    try {
      if (pos.empty()) return usage();
      Trajectory t = open_file(crystal_path(pos[0]), o.count("--topology") ? o["--topology"] : "");   // a path or a library crystal
      System s0 = t.frames() ? t.frame(0) : t.topology;
      const bool quiet = o.count("--quiet");
      const ForceField ff0 = cli_forcefield(s0, o, quiet);
      QhaOptions q;
      if (o.count("--cutoff")) q.phonon.energy.cutoff = std::stod(o["--cutoff"]);
      run_cutoff(q.phonon.energy, ff0, o, quiet);
      electrostatics(q.phonon.energy, o);
      q.phonon.energy.coulomb = !o.count("--no-coulomb");
      q.phonon.field_for = [&](const System& sup) { System c = sup; return cli_forcefield(c, o, true); };
      auto ints = [&](const std::string& key) {
        std::array<int, 3> v{0, 0, 0};
        std::string w = o[key];
        for (auto& ch : w) if (ch == ',' || ch == 'x') ch = ' ';
        std::istringstream in(w);
        if (!(in >> v[0] >> v[1] >> v[2])) throw std::invalid_argument(key + ": three numbers, comma-separated");
        for (int x : v)
          if (x < 0) throw std::invalid_argument(key + ": whole numbers of 0 (chosen automatically) or more");
        return v;
      };
      if (o.count("--supercell")) q.phonon.supercell = ints("--supercell");
      if (o.count("--mesh")) q.phonon.mesh = ints("--mesh");
      if (o.count("--mode")) q.mode = o["--mode"];
      if (o.count("--range")) q.range = std::stod(o["--range"]);
      if (o.count("--points")) q.points = std::stoi(o["--points"]);
      if (o.count("--pressure")) q.pressure = std::stod(o["--pressure"]);
      q.relax_cell = !o.count("--no-relax-cell");
      if (o.count("--temps")) {
        q.temperatures.clear();
        std::stringstream ss(o["--temps"]);
        std::string x;
        while (std::getline(ss, x, ',')) q.temperatures.push_back(std::stod(x));
      }
      if (o.count("--lock")) {
        std::stringstream ss(o["--lock"]);
        std::string x;
        while (std::getline(ss, x, ',')) {
          if (x != "a" && x != "b" && x != "c") throw std::invalid_argument("--lock takes a, b, c (comma-separated)");
          q.axes[size_t(x[0] - 'a')] = false;
        }
      }
      if (!quiet) q.progress = [&](const std::string& what, double f) { std::fprintf(stderr, "\r%3.0f %% %-60s", 100 * f, what.c_str()); return true; };
      const QhaResult r = quasi_harmonic(s0, q);
      if (!quiet) std::fprintf(stderr, "\n");
      for (const auto& n : r.notes) std::printf("%s\n", n.c_str());
      std::printf("%zu strained cells · %s · V0 %.3f Å³ · fit rms %.2e kcal/mol\n", r.samples.size(), r.field.c_str(), r.volume0, r.fit_rms);
      std::printf("%8s %9s %9s %9s %10s %10s %10s %10s %10s %8s %9s %9s %9s\n", "T/K", "a/Å", "b/Å", "c/Å", "V/Å3", "aV/1e-6K", "aa/1e-6K", "ab/1e-6K", "ac/1e-6K", "B/GPa", "Cv", "Cp", "S");
      for (const auto& w : r.rows)
        std::printf("%8.2f %9.4f %9.4f %9.4f %10.3f %10.2f %10.2f %10.2f %10.2f %8.3f %9.3f %9.3f %9.3f%s%s\n", w.T, w.length[0], w.length[1], w.length[2], w.volume, 1e6 * w.alpha_V,
                    1e6 * w.alpha[0], 1e6 * w.alpha[1], 1e6 * w.alpha[2], w.B, w.Cv, w.Cp, w.S, w.extrapolated ? "  extrapolated" : "", w.unstable ? "  unstable" : "");
      std::printf("Cv, Cp, S in cal/mol/K per cell\n");
      if (o.count("--csv")) {
        std::ofstream f(o["--csv"]);
        f << "T_K,a_A,b_A,c_A,V_A3,alpha_V_per_K,alpha_a_per_K,alpha_b_per_K,alpha_c_per_K,B_GPa,Cv,Cp,S,G_kcal_mol,extrapolated,unstable\n";
        char b[400];
        for (const auto& w : r.rows) {
          std::snprintf(b, sizeof b, "%.3f,%.6f,%.6f,%.6f,%.5f,%.6e,%.6e,%.6e,%.6e,%.5f,%.5f,%.5f,%.5f,%.6f,%d,%d\n", w.T, w.length[0], w.length[1], w.length[2], w.volume, w.alpha_V,
                        w.alpha[0], w.alpha[1], w.alpha[2], w.B, w.Cv, w.Cp, w.S, w.G, int(w.extrapolated), int(w.unstable));
          f << b;
        }
        std::printf("wrote %s\n", o["--csv"].c_str());
      }
      return 0;
    } catch (const std::exception& e) {
      std::fprintf(stderr, "caps qha: %s\n", e.what());
      keep_cli_error(argc, argv, e);
      return 1;
    }
  }
  if (cmd == "phonons") {
    // caps phonons CRYSTAL [--ff FF.json ...] [--cutoff Å] [--supercell n1,n2,n3] [--mesh m1,m2,m3] [--path 'G 0 0 0; X 0.5 0 0; ...']
    //              [--displacement 0.01] [--sigma 5] [--temps 100,300] [--csv STEM] [--json OUT]
    try {
      if (pos.empty()) return usage();
      Trajectory t = open_file(crystal_path(pos[0]), o.count("--topology") ? o["--topology"] : "");   // a path or a library crystal
      System s0 = t.frames() ? t.frame(0) : t.topology;
      const bool quiet = o.count("--quiet");
      const ForceField ff0 = cli_forcefield(s0, o, quiet);
      PhononOptions po;
      if (o.count("--cutoff")) po.energy.cutoff = std::stod(o["--cutoff"]);
      run_cutoff(po.energy, ff0, o, quiet);
      electrostatics(po.energy, o);
      po.energy.coulomb = !o.count("--no-coulomb");
      po.field_for = [&](const System& sup) { System c = sup; return cli_forcefield(c, o, true); };
      auto ints = [&](const std::string& key) {
        std::array<int, 3> v{0, 0, 0};
        std::string s = o[key];
        for (auto& ch : s) if (ch == ':' || ch == 'x') ch = ',';
        size_t at = 0;
        for (int k = 0; k < 3; ++k) {
          const size_t c = s.find(',', at);
          v[size_t(k)] = std::stoi(s.substr(at, c == std::string::npos ? std::string::npos : c - at));
          if (c == std::string::npos && k < 2) throw std::invalid_argument(key + ": three numbers, comma-separated");
          at = c + 1;
        }
        for (int x : v)
          if (x < 0) throw std::invalid_argument(key + ": whole numbers of 0 (chosen automatically) or more");
        return v;
      };
      if (o.count("--supercell")) po.supercell = ints("--supercell");
      if (o.count("--mesh")) po.mesh = ints("--mesh");
      if (o.count("--displacement")) {
        po.step = std::stod(o["--displacement"]);
        if (!(po.step > 0 && po.step <= 0.5)) throw std::invalid_argument("--displacement: a length above 0 and at most 0.5 Å, e.g. 0.01");
      }
      if (o.count("--sigma")) po.dos_sigma = std::stod(o["--sigma"]);
      if (o.count("--temps")) {
        po.temperatures.clear();
        std::stringstream ss(o["--temps"]);
        std::string x;
        while (std::getline(ss, x, ',')) po.temperatures.push_back(std::stod(x));
      }
      if (o.count("--path")) {   // 'LABEL f1 f2 f3; LABEL f1 f2 f3; …'
        std::stringstream ss(o["--path"]);
        std::string seg;
        while (std::getline(ss, seg, ';')) {
          std::stringstream one(seg);
          std::string lab;
          std::array<double, 3> f{};
          if (!(one >> lab >> f[0] >> f[1] >> f[2])) throw std::invalid_argument("--path: 'LABEL f1 f2 f3; …' in fractions of the reciprocal axes");
          po.labels.push_back(lab == "G" ? "Γ" : lab), po.path.push_back(f);
        }
      }
      if (!quiet)
        po.progress = [&](const std::string& w, double f) { std::fprintf(stderr, "\r%-30s %3.0f %%", w.c_str(), 100 * f); return true; };
      const PhononResult r = phonons(s0, po);
      if (!quiet) std::fprintf(stderr, "\r%-40s\r", "");
      for (const auto& n : r.notes) std::printf("%s\n", n.c_str());
      std::printf("lowest wavenumber %.1f cm⁻¹ · highest %.1f cm⁻¹ · largest force at the reference %.3g kcal/mol/Å\n", r.min_wavenumber,
                  r.dos_w.empty() ? 0.0 : r.dos_w.back(), r.max_force);
      const char* ax[3] = {"a*", "b*", "c*"};
      for (int d = 0; d < 3; ++d)
        std::printf("sound speeds along %s: %.0f %.0f %.0f m/s\n", ax[d], r.sound[size_t(d)][0], r.sound[size_t(d)][1], r.sound[size_t(d)][2]);
      std::printf("%8s %12s %12s %12s %12s %12s   (per cell)\n", "T (K)", "ZPE", "U", "F", "S", "Cv");
      for (size_t i = 0; i < r.T.size(); ++i)
        std::printf("%8.2f %12.4f %12.4f %12.4f %12.4f %12.4f\n", r.T[i], r.zpe[i], r.U[i], r.F[i], r.S[i], r.Cv[i]);
      std::printf("%8s %12s %12s %12s %12s %12s\n", "", "kcal/mol", "kcal/mol", "kcal/mol", "cal/mol/K", "cal/mol/K");
      if (o.count("--csv")) {
        const std::string stem = o["--csv"];
        std::ofstream b(stem + "_bands.csv");
        b << "k_path_inv_angstrom";
        for (size_t j = 0; j < (r.bands.empty() ? 0 : r.bands[0].size()); ++j) b << ",band" << j + 1 << "_cm-1";
        b << "\n";
        for (size_t i = 0; i < r.bands.size(); ++i) {
          b << r.path_x[i];
          for (double v : r.bands[i]) b << "," << v;
          b << "\n";
        }
        std::ofstream d(stem + "_dos.csv");
        d << "wavenumber_cm-1,dos";
        for (const auto& l : r.pdos_labels) d << "," << l;
        d << "\n";
        for (size_t i = 0; i < r.dos_w.size(); ++i) {
          d << r.dos_w[i] << "," << r.dos[i];
          for (const auto& p : r.pdos) d << "," << p[i];
          d << "\n";
        }
        std::printf("wrote %s_bands.csv and %s_dos.csv\n", stem.c_str(), stem.c_str());
      }
      if (o.count("--json")) {
        Json j = Json::object();
        auto arr = [](const std::vector<double>& v) { Json a = Json::array(); for (double x : v) a.push_back(Json(x)); return a; };
        Json bands = Json::array(), ticks = Json::array(), snd = Json::array(), pd = Json::array();
        for (const auto& b : r.bands) bands.push_back(arr(b));
        for (const auto& l : r.tick_labels) ticks.push_back(Json(l));
        for (const auto& row : r.sound) snd.push_back(arr({row[0], row[1], row[2]}));
        for (const auto& p : r.pdos) pd.push_back(arr(p));
        j["path_x"] = arr(r.path_x), j["bands"] = bands, j["tick_x"] = arr(r.tick_x), j["tick_labels"] = ticks;
        j["dos_w"] = arr(r.dos_w), j["dos"] = arr(r.dos), j["pdos"] = pd, j["T"] = arr(r.T), j["zpe"] = arr(r.zpe), j["U"] = arr(r.U), j["F"] = arr(r.F);
        j["S"] = arr(r.S), j["Cv"] = arr(r.Cv), j["sound"] = snd, j["min_wavenumber"] = Json(r.min_wavenumber);
        std::ofstream(o["--json"]) << j.dump(1) << "\n";
        std::printf("wrote %s\n", o["--json"].c_str());
      }
      return 0;
    } catch (const std::exception& e) {
      std::fprintf(stderr, "caps phonons: %s\n", e.what());
      keep_cli_error(argc, argv, e);
      return 1;
    }
  }
  if (cmd == "index") {
    try {
      if (pos.empty()) return usage();
      const MeasuredPattern d = read_pattern(pos[0]);
      IndexOptions io;
      if (o.count("--wavelength")) io.wavelength = std::stod(o["--wavelength"]);
      if (o.count("--tolerance")) io.tolerance = std::stod(o["--tolerance"]);
      if (o.count("--max-length")) io.max_length = std::stod(o["--max-length"]);
      if (o.count("--lines")) io.lines = std::stoi(o["--lines"]);
      if (o.count("--unindexed")) io.unindexed = std::stoi(o["--unindexed"]);
      if (o.count("--max-boxes")) io.max_boxes = std::stol(o["--max-boxes"]);
      if (!o.count("--quiet")) io.progress = [&](const std::string& w, double f) { std::fprintf(stderr, "\r%-20s %3.0f %%", w.c_str(), 100 * f); return true; };
      if (o.count("--systems")) {
        io.systems.clear();
        std::stringstream ss(o["--systems"]);
        std::string w;
        while (std::getline(ss, w, ',')) io.systems.push_back(w);
      }
      PeakOptions po;
      po.lambda1 = io.wavelength;
      const auto pk = find_peaks(d, po);
      std::vector<double> tt;
      for (const auto& p : pk) tt.push_back(p.two_theta);
      std::printf("%zu peaks:", tt.size());
      for (size_t k = 0; k < tt.size() && k < 30; ++k) std::printf(" %.3f", tt[k]);
      std::printf("\n");
      std::fflush(stdout);
      const IndexResult r = index_pattern(tt, io);
      if (!o.count("--quiet")) std::fprintf(stderr, "\r%-30s\r", "");
      std::printf("%-13s %9s %9s %9s %8s %8s %8s %9s %8s %8s %s\n", "system", "a", "b", "c", "alpha", "beta", "gamma", "V", "M20", "F_N", "indexed");
      for (const auto& s2 : r.solutions)
        std::printf("%-13s %9.4f %9.4f %9.4f %8.3f %8.3f %8.3f %9.2f %8.1f %8.1f %d/%zu\n", s2.system.c_str(), s2.cell[0], s2.cell[1], s2.cell[2], s2.cell[3], s2.cell[4], s2.cell[5],
                    s2.volume, s2.M20, s2.FN, s2.indexed, r.lines.size());
      if (!r.solutions.empty()) {
        const auto& b = r.solutions[0];
        std::printf("the best cell's lines (2θ obs, calc, hkl):\n");
        for (size_t i = 0; i < b.hkl.size(); ++i)
          std::printf("  %8.3f %8.3f  %d %d %d\n", b.two_theta_obs[i], b.two_theta_calc[i], b.hkl[i][0], b.hkl[i][1], b.hkl[i][2]);
      }
      for (const auto& n : r.notes) std::printf("%s\n", n.c_str());
      std::printf("%.1f s\n", r.seconds);
      if (o.count("--json")) {
        auto arr = [](const auto& v) { Json a = Json::array(); for (double x : v) a.push_back(Json(x)); return a; };
        Json j = Json::object(), sols = Json::array(), nt = Json::array();
        for (const auto& s2 : r.solutions) {
          Json x = Json::object(), hk = Json::array();
          for (const auto& h : s2.hkl) {
            Json t3 = Json::array();
            for (int q : h) t3.push_back(Json(q));
            hk.push_back(t3);
          }
          x["system"] = s2.system, x["cell"] = arr(s2.cell), x["volume"] = s2.volume, x["M20"] = s2.M20, x["FN"] = s2.FN, x["indexed"] = double(s2.indexed);
          x["hkl"] = hk, x["two_theta_obs"] = arr(s2.two_theta_obs), x["two_theta_calc"] = arr(s2.two_theta_calc);
          sols.push_back(x);
        }
        for (const auto& n : r.notes) nt.push_back(Json(n));
        j["lines"] = arr(r.lines), j["solutions"] = sols, j["notes"] = nt, j["seconds"] = r.seconds;
        std::ofstream(o["--json"]) << j.dump(1) << "\n";
        std::printf("wrote %s\n", o["--json"].c_str());
      }
      return r.solutions.empty() ? 2 : 0;
    } catch (const std::exception& e) {
      std::fprintf(stderr, "caps index: %s\n", e.what());
      keep_cli_error(argc, argv, e);
      return 1;
    }
  }
  if (cmd == "solve") {
    try {
      if (pos.size() < 2 || !o.count("--cell")) return usage();
      const MeasuredPattern d = read_pattern(pos[0]);
      Trajectory t = open_file(pos[1], "");
      System mol = t.frames() ? t.frame(0) : t.topology;
      std::array<double, 6> cp{};
      {
        std::string v = o["--cell"];
        for (auto& c : v) if (c == ',') c = ' ';
        std::istringstream is(v);
        for (auto& x : cp) if (!(is >> x)) throw std::invalid_argument("--cell a,b,c,alpha,beta,gamma");
      }
      SolveOptions so;
      if (o.count("--space-group")) so.space_group = o["--space-group"];
      if (o.count("--wavelength")) so.refine.wavelength = std::stod(o["--wavelength"]);
      if (o.count("--steps")) so.steps = std::stoi(o["--steps"]);
      if (o.count("--restarts")) so.restarts = std::stoi(o["--restarts"]);
      if (o.count("--seed")) so.seed = std::stoull(o["--seed"]);
      {
        // a known instrument profile, held instead of the Pawley fit of it
        auto nums = [&](const std::string& key, size_t n) {
          std::string v = o[key];
          for (auto& c : v) if (c == ',') c = ' ';
          std::istringstream is(v);
          std::vector<double> x(n);
          for (auto& q : x) if (!(is >> q)) throw std::invalid_argument(key + ": " + std::to_string(n) + " numbers, comma-separated");
          return x;
        };
        if (o.count("--uvw")) { const auto v = nums("--uvw", 3); so.refine.U = v[0], so.refine.V = v[1], so.refine.W = v[2], so.fit_profile = false; }
        if (o.count("--xy")) { const auto v = nums("--xy", 2); so.refine.X = v[0], so.refine.Y = v[1], so.fit_profile = false; }
        if (o.count("--zero")) so.refine.zero = std::stod(o["--zero"]), so.fit_profile = false;
      }
      if (!o.count("--quiet")) so.progress = [&](const std::string& w, double f) { std::fprintf(stderr, "\r%-20s %3.0f %%", w.c_str(), 100 * f); return true; };
      const SolveResult r = solve_structure(d, mol, cp, so);
      if (!o.count("--quiet")) std::fprintf(stderr, "\r%-30s\r", "");
      std::printf("Rwp %.4f after %d moves (%d accepted, %d torsions)\n", r.rwp, r.steps, r.accepted, r.torsions);
      for (const auto& n : r.notes) std::printf("%s\n", n.c_str());
      const std::string out = o.count("-o") ? o["-o"] : (o.count("--out") ? o["--out"] : "solved.cif");
      write_structure_file(r.crystal, out);
      std::printf("wrote %s (%zu atoms)\n", out.c_str(), r.crystal.atoms.size());
      return 0;
    } catch (const std::exception& e) {
      std::fprintf(stderr, "caps solve: %s\n", e.what());
      keep_cli_error(argc, argv, e);
      return 1;
    }
  }
  if (cmd == "domains") {
    try {
      if (pos.empty()) return usage();
      Trajectory t = open_file(pos[0], o.count("--topology") ? o["--topology"] : "");
      const System s = t.frames() ? t.frame(t.frames() - 1) : t.topology;
      DomainOptions dop;
      if (o.count("--group")) dop.group = o["--group"];
      if (o.count("--grid")) dop.grid = std::stod(o["--grid"]);
      if (o.count("--sigma")) dop.sigma = std::stod(o["--sigma"]);
      if (o.count("--threshold")) dop.threshold = std::stod(o["--threshold"]);
      dop.above = !o.count("--below");
      if (o.count("--min-voxels")) dop.min_voxels = std::stoi(o["--min-voxels"]);
      const DomainResult r = density_domains(s, dop);
      for (const auto& n : r.notes) std::printf("%s\n", n.c_str());
      std::printf("grid %d × %d × %d, level %.4g %s\n", r.field.n[0], r.field.n[1], r.field.n[2], r.threshold, r.field.unit.c_str());
      std::printf("domain      volume Å³   R_eq Å   axes Å                    percolates  nearest Å\n");
      for (size_t k = 0; k < r.domains.size() && k < 40; ++k) {
        const auto& d = r.domains[k];
        std::printf("%6zu  %12.1f  %7.2f   %6.2f %6.2f %6.2f      %c%c%c        %8.2f\n", k + 1, d.volume, d.equivalent_radius, d.axes[0], d.axes[1], d.axes[2],
                    d.percolates[0] ? 'a' : '-', d.percolates[1] ? 'b' : '-', d.percolates[2] ? 'c' : '-', d.nearest);
      }
      if (r.domains.size() > 40) std::printf("… %zu more\n", r.domains.size() - 40);
      if (o.count("--csv")) {
        std::ofstream f(o["--csv"]);
        f << "domain,voxels,volume_A3,equivalent_radius_A,axis1_A,axis2_A,axis3_A,cx,cy,cz,percolates_a,percolates_b,percolates_c,nearest_A\n";
        for (size_t k = 0; k < r.domains.size(); ++k) {
          const auto& d = r.domains[k];
          f << k + 1 << "," << d.voxels << "," << d.volume << "," << d.equivalent_radius << "," << d.axes[0] << "," << d.axes[1] << "," << d.axes[2] << "," << d.centroid[0] << ","
            << d.centroid[1] << "," << d.centroid[2] << "," << d.percolates[0] << "," << d.percolates[1] << "," << d.percolates[2] << "," << d.nearest << "\n";
        }
        std::printf("wrote %s\n", o["--csv"].c_str());
      }
      if (o.count("--cube")) {
        write_grid(r.field, s, o["--cube"]);
        std::printf("wrote %s\n", o["--cube"].c_str());
      }
      if (o.count("--slice")) {   // 'px,py,pz; ux,uy,uz; vx,vy,vz; width,height; nx,ny' → CSV of x, y, value
        std::vector<std::vector<double>> parts;
        std::stringstream ss(o["--slice"]);
        for (std::string part; std::getline(ss, part, ';');) {
          std::vector<double> v;
          std::stringstream ps(part);
          for (std::string x; std::getline(ps, x, ',');) v.push_back(std::stod(x));
          parts.push_back(v);
        }
        if (parts.size() != 5 || parts[0].size() != 3 || parts[1].size() != 3 || parts[2].size() != 3 || parts[3].size() != 2 || parts[4].size() != 2)
          throw std::invalid_argument("--slice 'px,py,pz; ux,uy,uz; vx,vy,vz; width,height; nx,ny'");
        const int nx = int(parts[4][0]), ny = int(parts[4][1]);
        const auto vals = field_slice(r.field, {parts[0][0], parts[0][1], parts[0][2]}, {parts[1][0], parts[1][1], parts[1][2]}, {parts[2][0], parts[2][1], parts[2][2]},
                                      parts[3][0], parts[3][1], nx, ny);
        const std::string path = o.count("--slice-csv") ? o["--slice-csv"] : std::string("slice.csv");
        std::ofstream f(path);
        f << "u_A,v_A,value\n";
        for (int y = 0; y < ny; ++y)
          for (int x = 0; x < nx; ++x)
            f << (nx > 1 ? parts[3][0] * (double(x) / (nx - 1) - 0.5) : 0.0) << "," << (ny > 1 ? parts[3][1] * (double(y) / (ny - 1) - 0.5) : 0.0) << "," << vals[size_t(y * nx + x)] << "\n";
        std::printf("wrote %s (%d × %d samples)\n", path.c_str(), nx, ny);
      }
      return 0;
    } catch (const std::exception& e) {
      std::fprintf(stderr, "caps domains: %s\n", e.what());
      keep_cli_error(argc, argv, e);
      return 1;
    }
  }
  if (cmd == "adp") {
    try {
      if (pos.empty()) return usage();
      Trajectory t = open_file(pos[0], o.count("--topology") ? o["--topology"] : "");
      AnalyzeOptions ao;
      if (o.count("--first")) ao.first = std::stol(o["--first"]);
      if (o.count("--last")) ao.last = std::stol(o["--last"]);
      if (o.count("--stride")) ao.stride = std::stol(o["--stride"]);
      const auto fr = analysis_frames(t, ao);
      AdpOptions po;
      if (o.count("--supercell")) {
        std::string v = o["--supercell"];
        for (auto& c : v) if (c == 'x' || c == 'X') c = ',';
        if (std::sscanf(v.c_str(), "%d,%d,%d", &po.supercell[0], &po.supercell[1], &po.supercell[2]) != 3) throw std::invalid_argument("--supercell: na,nb,nc");
      }
      if (o.count("--probability")) po.probability = std::stod(o["--probability"]);
      po.remove_drift = !o.count("--no-drift");
      const AdpResult r = displacement_tensors(t, fr, po);
      const auto& rows = r.sites.empty() ? r.atoms : r.sites;
      std::vector<int> els;
      for (size_t k = 0; k < rows.size(); ++k) els.push_back(t.topology.atoms[r.sites.empty() ? k : r.site_atoms[k][0]].element);
      std::map<int, std::pair<double, int>> per;
      for (size_t k = 0; k < rows.size(); ++k) per[els[k]].first += rows[k].b, per[els[k]].second += 1;
      std::printf("%d frames · %zu %s · ellipsoids at %.0f %% (scale %.4f)\n", r.frames, rows.size(), r.sites.empty() ? "atoms" : "sites", 100 * po.probability, r.scale);
      for (const auto& [z, e] : per) std::printf("  %-3s mean B %.3f Å² (U_iso %.5f Å²) over %d\n", element(z).symbol, e.first / e.second, e.first / e.second / (8 * M_PI * M_PI), e.second);
      for (const auto& n : r.notes) std::printf("%s\n", n.c_str());
      if (o.count("--csv")) {
        std::ofstream f(o["--csv"]);
        f << (r.sites.empty() ? "atom" : "site") << ",element,x,y,z,U11,U22,U33,U12,U13,U23,U_iso,B,axis1,axis2,axis3,anisotropy\n";
        for (size_t k = 0; k < rows.size(); ++k) {
          const auto& a = rows[k];
          f << k + 1 << "," << element(els[k]).symbol << "," << a.mean[0] << "," << a.mean[1] << "," << a.mean[2];
          for (double u : a.U) f << "," << u;
          f << "," << a.u_iso << "," << a.b << "," << a.axes[0] << "," << a.axes[1] << "," << a.axes[2] << "," << a.anisotropy << "\n";
        }
        std::printf("wrote %s\n", o["--csv"].c_str());
      }
      if (o.count("--cif")) {
        write_cif_adp(els, rows, r.unit_cell, o["--cif"], t.topology.title.empty() ? "caps" : t.topology.title);
        std::printf("wrote %s (U_ij on the CIF's axes)\n", o["--cif"].c_str());
      }
      return 0;
    } catch (const std::exception& e) {
      std::fprintf(stderr, "caps adp: %s\n", e.what());
      keep_cli_error(argc, argv, e);
      return 1;
    }
  }
  if (cmd == "cif") {
    try {
      if (pos.empty()) return usage();
      CifDisorder dop;
      if (o.count("--disorder")) dop.mode = o["--disorder"];
      if (o.count("--group")) dop.group = o["--group"], dop.mode = o.count("--disorder") ? dop.mode : "group";
      if (o.count("--supercell")) {
        std::string v = o["--supercell"];
        for (auto& c : v) if (c == 'x' || c == 'X') c = ',';
        if (std::sscanf(v.c_str(), "%d,%d,%d", &dop.supercell[0], &dop.supercell[1], &dop.supercell[2]) != 3) throw std::invalid_argument("--supercell: na,nb,nc");
      }
      if (o.count("--seed")) dop.seed = std::stoull(o["--seed"]);
      CifDisorderReport rep;
      const System s = read_cif(pos[0], dop, &rep);
      if (s.atoms.size() > 15000) throw std::invalid_argument(std::to_string(s.atoms.size()) + " atoms (at most 15000): a smaller supercell");
      for (const auto& n : s.notes) std::printf("%s\n", n.c_str());
      std::printf("%d disorder assembl%s, %d groups, %d places with alternatives or vacancies (per cell)\n", rep.assemblies, rep.assemblies == 1 ? "y" : "ies", rep.groups, rep.places);
      std::printf("element   occupancies give   written\n");
      for (const auto& [z, e] : rep.expected) std::printf("%-8s %17.2f %9d\n", element(z).symbol, e, rep.obtained.count(z) ? rep.obtained.at(z) : 0);
      if (o.count("-o")) {
        const std::string out = o["-o"];
        auto ends = [&](const char* e) { return out.size() > 4 && out.substr(out.size() - 4) == e; };
        if (ends(".cif")) write_cif(s, out);
        else if (ends(".pdb")) write_pdb(s, out);
        else if (ends(".xyz")) write_xyz(s, out);
        else write_lammps_data(s, out);
        std::printf("wrote %s\n", out.c_str());
      }
      if (o.count("--json") && o["--json"] != "1") {
        Json j = Json::object(), ex = Json::object(), ob = Json::object();
        for (const auto& [z, e] : rep.expected) ex[element(z).symbol] = e;
        for (const auto& [z, n] : rep.obtained) ob[element(z).symbol] = double(n);
        j["assemblies"] = double(rep.assemblies), j["groups"] = double(rep.groups), j["places"] = double(rep.places), j["expected"] = ex, j["obtained"] = ob;
        std::ofstream(o["--json"]) << j.dump(1) << "\n";
        std::printf("wrote %s\n", o["--json"].c_str());
      }
      return 0;
    } catch (const std::exception& e) {
      std::fprintf(stderr, "caps cif: %s\n", e.what());
      keep_cli_error(argc, argv, e);
      return 1;
    }
  }
  if (cmd == "graft") {
    try {
      if (pos.empty() || !o.count("--units") || !o.count("--sites") || !o.count("-o")) return usage();
      Trajectory t = open_file(pos[0], o.count("--topology") ? o["--topology"] : "");
      const System host = t.frames() ? t.frame(t.frames() - 1) : t.topology;
      auto list = [](const std::string& x) {
        std::vector<std::string> v;
        std::stringstream ss(x);
        for (std::string y; std::getline(ss, y, ',');) v.push_back(y);
        return v;
      };
      ChainSpec spec;
      for (const auto& u : list(o["--units"])) spec.units.push_back({u, u});
      spec.dp = o.count("--dp") ? std::stoi(o["--dp"]) : 10;
      if (spec.dp < 1) throw std::invalid_argument("--dp must be at least 1 (got " + o["--dp"] + ")");
      require_structure_name(o["-o"]);
      if (o.count("--sequence")) spec.sequence = sequence_from_string(o["--sequence"]);
      if (o.count("--weights")) for (const auto& w : list(o["--weights"])) spec.weights.push_back(std::stod(w));
      if (spec.units.size() > 1 && spec.sequence == Sequence::Homopolymer) spec.sequence = spec.weights.empty() ? Sequence::Alternating : Sequence::Random;
      if (o.count("--ff")) spec.forcefield = ff_path(o["--ff"]);
      ChainGraftOptions go;
      for (const auto& x : list(o["--sites"])) {
        const long a = std::stol(x);
        if (a < 1) throw std::invalid_argument("--sites: atoms from 1");
        go.sites.push_back(uint32_t(a - 1));
      }
      go.equivalent = o.count("--equivalent") > 0;
      if (o.count("--probability")) go.probability = std::stod(o["--probability"]);
      if (o.count("--per-site")) go.per_site = std::stoi(o["--per-site"]);
      if (o.count("--max")) go.max_grafts = std::stoi(o["--max"]);
      if (o.count("--min-spacing")) go.min_spacing = std::stod(o["--min-spacing"]);
      if (o.count("--padding")) go.padding = std::stod(o["--padding"]);
      if (o.count("--seed")) go.grow.seed = std::stoull(o["--seed"]);
      ChainGraftReport rep;
      const System s = graft_chains(host, spec, go, &rep);
      if (s.atoms.size() > 15000) throw std::invalid_argument(std::to_string(s.atoms.size()) + " atoms (at most 15000): fewer or shorter grafts");
      const std::string out = o["-o"];
      write_structure_file(s, out);   // the format the name asks for (.mol2, .gro, .cif …), LAMMPS data for .data
      for (const auto& n : rep.notes) std::printf("%s\n", n.c_str());
      std::printf("%d grafts of %d units on %d of %d candidate sites; %zu atoms · wrote %s\n", rep.grafted, spec.dp,
                  int(std::set<uint32_t>(rep.sites.begin(), rep.sites.end()).size()), rep.candidates, s.atoms.size(), out.c_str());
      return 0;
    } catch (const std::exception& e) {
      std::fprintf(stderr, "caps graft: %s\n", e.what());
      keep_cli_error(argc, argv, e);
      return 1;
    }
  }
  if (cmd == "copolymer") {
    try {
      if (!o.count("--feed")) return usage();
      std::vector<double> feed;
      {
        std::stringstream ss(o["--feed"]);
        for (std::string x; std::getline(ss, x, ',');) feed.push_back(std::stod(x));
      }
      const int n = int(feed.size());
      if (n < 2) throw std::invalid_argument("--feed: two or more units");
      ChainSpec spec;
      spec.units.resize(size_t(n));
      spec.weights = feed;
      copolymer_options(spec, o);
      spec.sequence = !spec.transition.empty() ? Sequence::Markov : !spec.penultimate.empty() ? Sequence::Penultimate : Sequence::Terminal;
      if (spec.sequence == Sequence::Terminal && spec.reactivity.empty() && n > 2) throw std::invalid_argument("three or more units need --reactivity");
      if (spec.sequence == Sequence::Penultimate && n != 2) throw std::invalid_argument("--penultimate is a model of two units: give --feed two fractions");
      SequenceModel m;
      std::vector<double> r = spec.reactivity;
      if (r.empty() && n == 2) r = {1, spec.r1, spec.r2, 1};
      if (spec.sequence == Sequence::Markov) m = transition_model(spec.transition, n);
      else if (spec.sequence == Sequence::Penultimate) m = penultimate_model(feed[0] / (feed[0] + feed[1]), spec.penultimate[0], spec.penultimate[1], spec.penultimate[2], spec.penultimate[3]);
      else m = terminal_model(feed, r);
      std::printf("%s model, %d units\nunit  composition  run length\n", to_string(spec.sequence), n);
      for (int i = 0; i < n; ++i) std::printf("%4c  %11.6f  %10.4f\n", char('A' + i), m.composition[size_t(i)], m.run_length[size_t(i)]);
      std::printf("dyads (row: first unit)\n");
      for (int i = 0; i < n; ++i) {
        std::printf("  %c ", char('A' + i));
        for (int j = 0; j < n; ++j) std::printf(" %9.6f", m.dyads[size_t(i * n + j)]);
        std::printf("\n");
      }
      std::vector<DriftPoint> drift;
      if (o.count("--conversion")) {
        if (spec.sequence == Sequence::Markov) throw std::invalid_argument("--conversion: a transition matrix has no feed to drift");
        drift = composition_drift(feed, r, std::stod(o["--conversion"]), o.count("--points") ? std::stoi(o["--points"]) : 21, spec.penultimate);
        std::printf("drift: conversion, feed, instantaneous and cumulative copolymer (unit A%s)\n", n > 2 ? ", then the others in the CSV" : "");
        for (const auto& p : drift) std::printf("  %6.3f  %9.6f  %9.6f  %9.6f\n", p.conversion, p.feed[0], p.instantaneous[0], p.cumulative[0]);
      }
      if (o.count("--csv")) {
        std::ofstream f(o["--csv"]);
        if (drift.empty()) {
          f << "unit,composition,run_length\n";
          for (int i = 0; i < n; ++i) f << char('A' + i) << "," << m.composition[size_t(i)] << "," << m.run_length[size_t(i)] << "\n";
        } else {
          f << "conversion";
          for (const char* w : {"feed", "instantaneous", "cumulative"})
            for (int i = 0; i < n; ++i) f << "," << w << "_" << char('A' + i);
          f << "\n";
          for (const auto& p : drift) {
            f << p.conversion;
            for (const auto* v : {&p.feed, &p.instantaneous, &p.cumulative})
              for (double x : *v) f << "," << x;
            f << "\n";
          }
        }
        std::printf("wrote %s\n", o["--csv"].c_str());
      }
      if (o.count("--json") && o["--json"] != "1") {
        // the model and the drift table above, as the C API's caps_copolymer_model gives them:
        // no "built" block (this command builds no chains), and the drift from the same --points
        const Json full = sequence_report(spec, {});
        Json r = Json::object();
        r["kind"] = full["kind"], r["units"] = full["units"];
        if (full.has("model")) r["model"] = full["model"];
        if (!drift.empty()) {
          auto arr = [](const std::vector<double>& v) {
            Json a = Json::array();
            for (double x : v) a.push_back(Json(x));
            return a;
          };
          Json d = Json::array();
          for (const auto& p : drift) {
            Json q = Json::object();
            q["conversion"] = p.conversion, q["feed"] = arr(p.feed), q["instantaneous"] = arr(p.instantaneous), q["cumulative"] = arr(p.cumulative);
            d.push_back(q);
          }
          r["drift"] = d;
        }
        std::ofstream(o["--json"]) << r.dump(1) << "\n";
        std::printf("wrote %s\n", o["--json"].c_str());
      }
      return 0;
    } catch (const std::exception& e) {
      std::fprintf(stderr, "caps copolymer: %s\n", e.what());
      keep_cli_error(argc, argv, e);
      return 1;
    }
  }
  if (cmd == "zmelt") {
    try {
      if (pos.empty()) return usage();
      if (!o.count("--temps")) throw std::invalid_argument("--temps T1,T2,... is required (the NVE starting temperatures, K)");
      Trajectory t = open_file(pos[0], o.count("--topology") ? o["--topology"] : "");
      System s0 = t.frames() ? t.frame(0) : t.topology;
      if (s0.atoms.size() > 15000) throw std::invalid_argument("the cell has " + std::to_string(s0.atoms.size()) + " atoms (at most 15000)");
      const bool quiet = o.count("--quiet");
      ZMethodOptions zo;
      zo.field = std::make_shared<ForceField>(cli_forcefield(s0, o, quiet));
      zo.energy.cutoff = o.count("--cutoff") ? std::stod(o["--cutoff"]) : 10.0;
      zo.energy.tail = false;   // a constant at the crystal's volume; off as in the Studio, the pair shifted at the cut-off
      electrostatics(zo.energy, o);
      {
        std::stringstream ss(o["--temps"]);
        std::string x;
        while (std::getline(ss, x, ',')) zo.temperatures.push_back(std::stod(x));
      }
      if (o.count("--steps")) zo.steps = std::stoll(o["--steps"]);
      if (o.count("--dt")) zo.dt = std::stod(o["--dt"]);
      if (o.count("--bisect")) zo.bisect = std::stoi(o["--bisect"]);
      if (o.count("--sample-every")) zo.sample_every = std::stoi(o["--sample-every"]);
      if (o.count("--seed")) zo.seed = std::stoull(o["--seed"]);
      if (!quiet) zo.progress = [&](const std::string& w, double f) { std::fprintf(stderr, "\r%-36s %3.0f %%", w.c_str(), 100 * f); return true; };
      const ZResult r = z_method(s0, zo);
      if (!quiet) std::fprintf(stderr, "\r%-44s\r", "");
      std::printf("order parameter on %zu reciprocal vectors, |G| %.4f 1/Å, perfect cell %.4f\n", r.g.size(), std::sqrt(r.g[0][0] * r.g[0][0] + r.g[0][1] * r.g[0][1] + r.g[0][2] * r.g[0][2]), r.order0);
      std::printf("%8s %14s %9s %8s %10s %7s %8s  state\n", "T0 (K)", "E (kcal/mol/at)", "T (K)", "±", "P (atm)", "order", "MSD Å²");
      for (const auto& p : r.points)
        std::printf("%8.1f %14.6f %9.2f %8.2f %10.0f %7.3f %8.2f  %s\n", p.t0, p.energy, p.temperature, p.temperature_sd, p.pressure, p.order, p.msd, p.melted ? "melted" : "crystal");
      for (const auto& n : r.notes) std::printf("%s\n", n.c_str());
      if (o.count("--csv")) {
        std::ofstream f(o["--csv"]);
        f << "t0_K,energy_kcal_mol_atom,temperature_K,temperature_sd_K,pressure_atm,order,msd_A2,melted\n";
        for (const auto& p : r.points) f << p.t0 << "," << p.energy << "," << p.temperature << "," << p.temperature_sd << "," << p.pressure << "," << p.order << "," << p.msd << "," << (p.melted ? 1 : 0) << "\n";
        std::printf("wrote %s\n", o["--csv"].c_str());
      }
      if (o.count("--json") && o["--json"] != "1") {
        std::ofstream(o["--json"]) << z_method_json(r).dump(1) << "\n";
        std::printf("wrote %s\n", o["--json"].c_str());
      }
      return 0;
    } catch (const std::exception& e) {
      std::fprintf(stderr, "caps zmelt: %s\n", e.what());
      keep_cli_error(argc, argv, e);
      return 1;
    }
  }
  if (cmd == "sqs" || cmd == "enumerate") {
    try {
      if (pos.empty() || !o.count("--from") || !o.count("--to")) return usage();
      Trajectory t = open_file(pos[0], o.count("--topology") ? o["--topology"] : "");
      System s0 = t.frames() ? t.frame(0) : t.topology;
      if (o.count("--supercell")) {
        std::string v = o["--supercell"];
        for (auto& c : v) if (c == 'x' || c == 'X') c = ',';
        int n[3] = {1, 1, 1};
        if (std::sscanf(v.c_str(), "%d,%d,%d", &n[0], &n[1], &n[2]) != 3 || n[0] < 1 || n[1] < 1 || n[2] < 1) throw std::invalid_argument("--supercell: na,nb,nc");
        s0 = supercell(s0, n[0], n[1], n[2]);
      }
      if (s0.atoms.size() > 15000) throw std::invalid_argument("the cell has " + std::to_string(s0.atoms.size()) + " atoms (at most 15000)");
      const int from = element_from_symbol(o["--from"]), to = element_from_symbol(o["--to"]);
      if (!from || !to) throw std::invalid_argument("--from and --to: element symbols");
      if (cmd == "sqs") {
        if (!o.count("-o")) return usage();
        SqsOptions so;
        so.from = from, so.to = to;
        if (o.count("--fraction")) so.fraction = std::stod(o["--fraction"]);
        if (o.count("--shells")) so.shells = std::stoi(o["--shells"]);
        if (o.count("--steps")) so.steps = std::stoi(o["--steps"]);
        if (o.count("--restarts")) so.restarts = std::stoi(o["--restarts"]);
        if (o.count("--seed")) so.seed = std::stoull(o["--seed"]);
        const SqsResult r = special_quasirandom(s0, so);
        std::printf("%d of %d %s sites take %s (x = %.4f); target (1 - 2x)^2 = %.4f\n", r.substituted, r.sites, element(from).symbol, element(to).symbol,
                    double(r.substituted) / r.sites, r.target);
        std::printf("shell  distance (Å)  correlation\n");
        for (size_t m = 0; m < r.correlation.size(); ++m) std::printf("%5zu  %12.5f  %11.5f\n", m + 1, r.shell_distance[m], r.correlation[m]);
        std::printf("objective %.6f\n", r.objective);
        write_disordered(r.structure, o["-o"]);
        std::printf("wrote %s\n", o["-o"].c_str());
      } else {
        EnumOptions eo;
        eo.from = from, eo.to = to;
        if (o.count("--count")) eo.count = std::stoi(o["--count"]);
        if (o.count("--max-configurations")) eo.max_configurations = size_t(std::stoull(o["--max-configurations"]));
        if (o.count("--tolerance")) eo.tolerance = std::stod(o["--tolerance"]);
        const EnumResult r = enumerate_occupations(s0, eo);
        std::printf("%zu occupations of %zu %s sites by %d %s: %zu symmetry-distinct (%d operations)\n", r.total, r.sites.size(), element(from).symbol, eo.count,
                    element(to).symbol, r.configurations.size(), r.operations);
        std::printf("config  multiplicity  sites\n");
        for (size_t k = 0; k < r.configurations.size(); ++k) {
          std::printf("%6zu  %12zu ", k + 1, r.multiplicity[k]);
          for (int x : r.configurations[k]) std::printf(" %u", r.sites[size_t(x)] + 1);
          std::printf("\n");
        }
        if (o.count("--write")) {
          std::filesystem::create_directories(o["--write"]);
          for (size_t k = 0; k < r.configurations.size(); ++k)
            write_cif(occupied(s0, r, k, to), (std::filesystem::path(o["--write"]) / ("config_" + std::to_string(k + 1) + ".cif")).string());
          std::printf("wrote %zu CIF files to %s\n", r.configurations.size(), o["--write"].c_str());
        }
      }
      return 0;
    } catch (const std::exception& e) {
      std::fprintf(stderr, "caps %s: %s\n", cmd.c_str(), e.what());
      keep_cli_error(argc, argv, e);
      return 1;
    }
  }
  if (cmd == "dielectric") {
    try {
      if (pos.empty()) return usage();
      Trajectory t = open_file(crystal_path(pos[0]), o.count("--topology") ? o["--topology"] : "");   // a path or a library crystal
      System s0 = t.frames() ? t.frame(0) : t.topology;
      const bool quiet = o.count("--quiet");
      (void)cli_forcefield(s0, o, quiet);   // shells added for a shell model
      DielectricOptions dop;
      dop.energy.cutoff = o.count("--cutoff") ? std::stod(o["--cutoff"]) : 10.0;
      dop.energy.tail = false;
      dop.energy.electrostatics = EnergyOptions::Electrostatics::PME;
      electrostatics(dop.energy, o);
      dop.relax = !o.count("--no-relax");
      dop.relax_cell = !o.count("--fixed-cell");
      dop.piezo = !o.count("--no-piezo");
      if (o.count("--damping")) dop.damping = std::stod(o["--damping"]);
      if (o.count("--ff") && !is_uff(o["--ff"])) {
        try {
          const FFDef def = load_forcefield(ff_path(o["--ff"]));
          std::set<std::string> sh;
          for (const auto& [c, x] : def.shells) sh.insert(x);
          for (const auto& a : s0.atoms) dop.shell.push_back(sh.count(a.name) ? 1 : 0);
        } catch (const std::exception&) {
        }
      }
      dop.field_for = [o](const System& x) mutable { System c = x; return cli_forcefield(c, o, true); };
      if (!quiet) dop.progress = [&](const std::string& w, double f) { std::fprintf(stderr, "\r%-30s %3.0f %%", w.c_str(), 100 * f); return true; };
      const DielectricResult r = crystal_dielectric(s0, dop);
      if (!quiet) std::fprintf(stderr, "\r%-40s\r", "");
      auto m3 = [](const char* n, const Mat3x3& m) {
        std::printf("%s\n", n);
        for (const auto& row : m) std::printf("  %10.4f %10.4f %10.4f\n", row[0], row[1], row[2]);
      };
      std::printf("cell a b c %.5f %.5f %.5f Å, angles %.3f %.3f %.3f° (%s), largest force %.1e kcal/mol/Å\n", r.cell[0], r.cell[1], r.cell[2], r.cell[3], r.cell[4], r.cell[5],
                  dop.relax ? (dop.relax_cell ? "relaxed" : "atoms relaxed in the cell") : (r.shells ? "cores as given, shells relaxed" : "as given"), r.max_force);
      m3("static dielectric tensor ε⁰", r.eps0);
      m3("high-frequency dielectric tensor ε^∞", r.epsinf);
      std::printf("principal values: ε⁰ %.4f %.4f %.4f · ε^∞ %.4f %.4f %.4f", r.eps0_principal[0], r.eps0_principal[1], r.eps0_principal[2], r.epsinf_principal[0],
                  r.epsinf_principal[1], r.epsinf_principal[2]);
      if (r.shells) std::printf(" · refractive indices %.4f %.4f %.4f", r.refractive[0], r.refractive[1], r.refractive[2]);
      std::printf("\nΓ modes (cm⁻¹):");
      for (double w : r.modes) std::printf(" %.1f", w);
      std::printf("\n");
      if (r.piezo) {
        std::printf("piezoelectric e (C/m²), columns xx yy zz yz xz xy\n");
        for (const auto& row : r.e) std::printf("  %9.5f %9.5f %9.5f %9.5f %9.5f %9.5f\n", row[0], row[1], row[2], row[3], row[4], row[5]);
        std::printf("piezoelectric d (pC/N)\n");
        for (const auto& row : r.d) std::printf("  %9.4f %9.4f %9.4f %9.4f %9.4f %9.4f\n", row[0], row[1], row[2], row[3], row[4], row[5]);
      }
      for (const auto& n : r.notes) std::printf("%s\n", n.c_str());
      if (o.count("--csv")) {
        std::ofstream f(o["--csv"] + "_reflectivity.csv");
        f << "wavenumber_cm-1,R1,R2,R3,eps_re1,eps_im1,eps_re2,eps_im2,eps_re3,eps_im3\n";
        for (size_t i = 0; i < r.w.size(); ++i) {
          f << r.w[i];
          for (int k = 0; k < 3; ++k) f << "," << r.reflectivity[k][i];
          for (int k = 0; k < 3; ++k) f << "," << r.eps_re[k][i] << "," << r.eps_im[k][i];
          f << "\n";
        }
        std::printf("wrote %s_reflectivity.csv\n", o["--csv"].c_str());
      }
      if (o.count("--json") && o["--json"] != "1") {
        Json j = Json::object();
        auto arr = [](const auto& v) { Json a = Json::array(); for (double x : v) a.push_back(Json(x)); return a; };
        Json e0 = Json::array(), ei = Json::array(), ee = Json::array(), dd = Json::array();
        for (int k = 0; k < 3; ++k) e0.push_back(arr(r.eps0[size_t(k)])), ei.push_back(arr(r.epsinf[size_t(k)])), ee.push_back(arr(r.e[size_t(k)])), dd.push_back(arr(r.d[size_t(k)]));
        j["eps0"] = e0, j["epsinf"] = ei, j["e"] = ee, j["d"] = dd, j["modes"] = arr(r.modes), j["refractive"] = arr(r.refractive);
        std::ofstream(o["--json"]) << j.dump(1) << "\n";
        std::printf("wrote %s\n", o["--json"].c_str());
      }
      return 0;
    } catch (const std::exception& e) {
      std::fprintf(stderr, "caps dielectric: %s\n", e.what());
      keep_cli_error(argc, argv, e);
      return 1;
    }
  }
  if (cmd == "refine") {
    try {
      if (pos.size() < 2) return usage();
      const MeasuredPattern data = read_pattern(pos[0]);
      std::vector<RefinePhase> phases;
      for (size_t k = 1; k < pos.size(); ++k) {
        RefinePhase p;
        p.name = std::filesystem::path(pos[k]).stem().string();
        Trajectory t = open_file(pos[k], "");
        p.cell = t.frames() ? t.frame(0) : t.topology;
        p.b_iso = o.count("--b-iso");
        if (o.count("--preferred")) {
          std::string v = o["--preferred"];
          for (auto& c : v) if (c == ',') c = ' ';
          std::istringstream is(v);
          is >> p.preferred[0] >> p.preferred[1] >> p.preferred[2];
          p.march = true;
        }
        phases.push_back(p);
      }
      RefineOptions ro;
      ro.mode = o.count("--pawley") ? "pawley" : "rietveld";
      if (o.count("--neutron")) ro.radiation = "neutron";
      if (o.count("--anode")) ro.anode = o["--anode"];
      if (o.count("--wavelength")) ro.wavelength = std::stod(o["--wavelength"]);
      if (o.count("--background")) ro.background = std::stoi(o["--background"]);
      if (o.count("--iterations")) ro.iterations = std::stoi(o["--iterations"]);
      if (o.count("--two-theta")) {
        std::string v = o["--two-theta"];
        for (auto& c : v) if (c == ',' || c == ':') c = ' ';
        std::istringstream is(v);
        is >> ro.two_theta_min >> ro.two_theta_max;
      }
      if (o.count("--refine")) {
        ro.refine_U = ro.refine_V = ro.refine_W = ro.refine_X = ro.refine_Y = ro.refine_zero = false;
        std::stringstream ss(o["--refine"]);
        std::string w;
        while (std::getline(ss, w, ',')) {
          if (w == "U") ro.refine_U = true;
          else if (w == "V") ro.refine_V = true;
          else if (w == "W") ro.refine_W = true;
          else if (w == "X") ro.refine_X = true;
          else if (w == "Y") ro.refine_Y = true;
          else if (w == "zero") ro.refine_zero = true;
          else throw std::invalid_argument("--refine: U V W X Y zero");
        }
      }
      const RefineResult r = refine_pattern(data, phases, ro);
      std::printf("%s · Rp %.4f · Rwp %.4f · Rexp %.4f · χ² %.3f · %d iterations, %d parameters\n", ro.mode.c_str(), r.rp, r.rwp, r.rexp, r.chi2, r.iterations, r.parameters);
      for (const auto& p : r.phases) {
        std::printf("%s: a b c %.5f %.5f %.5f Å (± %.1e %.1e %.1e), angles %.3f %.3f %.3f°, V %.3f Å³", p.name.c_str(), p.lattice[0], p.lattice[1], p.lattice[2],
                    p.lattice_esd[0], p.lattice_esd[1], p.lattice_esd[2], p.lattice[3], p.lattice[4], p.lattice[5], p.volume);
        if (ro.mode == "rietveld") std::printf(", scale %.5g, weight fraction %.4f, B %.3f Å²", p.scale, p.weight_fraction, p.b);
        std::printf("\n");
      }
      std::printf("profile: U %.5g V %.5g W %.5g X %.5g Y %.5g zero %.5f°", r.U, r.V, r.W, r.X, r.Y, r.zero);
      if (r.size > 0) std::printf(" · crystallite size %.0f Å (Scherrer, from Y)", r.size);
      if (r.strain > 0) std::printf(" · microstrain %.2e (from X)", r.strain);
      std::printf("\n");
      for (const auto& n : r.notes) std::printf("%s\n", n.c_str());
      if (o.count("--csv")) {
        std::ofstream f(o["--csv"]);
        f << "two_theta,y_obs,y_calc,difference,background\n";
        for (size_t i = 0; i < r.two_theta.size(); ++i)
          f << r.two_theta[i] << "," << r.y_obs[i] << "," << r.y_calc[i] << "," << r.y_obs[i] - r.y_calc[i] << "," << r.background[i] << "\n";
        std::printf("wrote %s\n", o["--csv"].c_str());
      }
      if (o.count("--json") && o["--json"] != "1") {
        std::ofstream(o["--json"]) << refine_json(r).dump(1) << "\n";
        std::printf("wrote %s\n", o["--json"].c_str());
      }
      return 0;
    } catch (const std::exception& e) {
      std::fprintf(stderr, "caps refine: %s\n", e.what());
      keep_cli_error(argc, argv, e);
      return 1;
    }
  }
  if (cmd == "qspr") {
    try {
      if (pos.empty()) return usage();
      const std::string op = pos[0];
      // a CSV table: header names, then rows; numeric cells as numbers, the rest kept as text
      struct Table {
        std::vector<std::string> names;
        std::vector<std::vector<std::string>> cells;
      };
      auto read_csv = [](const std::string& path) {
        std::ifstream in(path);
        if (!in) throw std::runtime_error("cannot read " + path);
        Table t;
        std::string line;
        auto split = [](const std::string& l) { return split_csv_record(l, ','); };   // the reader caps table uses
        if (!std::getline(in, line)) throw std::runtime_error(path + " is empty");
        t.names = split(line);
        while (std::getline(in, line))
          if (!line.empty()) t.cells.push_back(split(line));
        return t;
      };
      auto column = [](const Table& t, const std::string& name) {
        for (size_t c = 0; c < t.names.size(); ++c)
          if (t.names[c] == name) return c;
        throw std::invalid_argument("no column '" + name + "'");
      };
      auto numeric = [](const std::string& x, double& v) {
        try {
          size_t k = 0;
          v = std::stod(x, &k);
          return k == x.size();
        } catch (...) {
          return false;
        }
      };
      // the x columns: the given ones, else every numeric column other than the excluded
      auto xcols = [&](const Table& t, const std::vector<size_t>& exclude) {
        std::vector<size_t> c;
        if (o.count("--x")) {
          std::stringstream ss(o["--x"]);
          std::string w;
          while (std::getline(ss, w, ',')) c.push_back(column(t, w));
          return c;
        }
        for (size_t k = 0; k < t.names.size(); ++k) {
          if (std::find(exclude.begin(), exclude.end(), k) != exclude.end()) continue;
          bool num = !t.cells.empty();
          double v;
          for (const auto& r : t.cells) num = num && k < r.size() && numeric(r[k], v);
          if (num) c.push_back(k);
        }
        return c;
      };
      auto value = [&](const Table& t, size_t r, size_t c) {
        double v;
        if (c >= t.cells[r].size() || !numeric(t.cells[r][c], v)) throw std::invalid_argument("row " + std::to_string(r + 1) + ", column " + t.names[c] + ": not a number");
        return v;
      };
      if (op == "descriptors") {
        std::vector<std::pair<std::string, System>> mols;
        for (size_t k = 1; k < pos.size(); ++k) {
          Trajectory t = open_file(pos[k], "");
          System s = t.frames() ? t.frame(0) : t.topology;
          s.cell = Cell{};
          mols.push_back({pos[k], s});
        }
        if (o.count("--smiles")) {
          std::stringstream ss(o["--smiles"]);
          std::string w;
          while (std::getline(ss, w, ';'))
            if (!w.empty()) mols.push_back({w, build_molecule(w).system});
        }
        if (mols.empty()) throw std::invalid_argument("give structure files or --smiles");
        std::vector<std::vector<Descriptor>> all;
        for (auto& [name, s] : mols) {
          std::string used;
          all.push_back(molecular_descriptors(s, &used));
          if (!o.count("--quiet")) std::printf("%s: %zu atoms, charges %s\n", name.c_str(), s.atoms.size(), used.c_str());
        }
        auto num = [](double v) { char b[40]; std::snprintf(b, sizeof b, "%.6g", v); return std::isfinite(v) ? std::string(b) : std::string(""); };
        if (o.count("--csv")) {
          std::ofstream f(o["--csv"]);
          f << "structure";
          for (const auto& d : all[0]) f << "," << d.name;
          f << "\n";
          for (size_t m = 0; m < all.size(); ++m) {
            f << "\"" << mols[m].first << "\"";
            for (const auto& d : all[m]) f << "," << num(d.value);
            f << "\n";
          }
          std::printf("wrote %s (%zu structures, %zu descriptors)\n", o["--csv"].c_str(), all.size(), all[0].size());
        } else {
          for (size_t k = 0; k < all[0].size(); ++k) {
            std::printf("%-20s", all[0][k].name.c_str());
            for (const auto& a : all) std::printf(" %14s", num(a[k].value).c_str());
            std::printf("  %s%s%s\n", all[0][k].unit.c_str(), all[0][k].unit.empty() ? "" : "  ", all[0][k].what.c_str());
          }
        }
        std::printf("AlogP and molar refractivity are not given: their atom-type tables are not in CAPS\n");
        return 0;
      }
      if (pos.size() < 2) return usage();
      const Table t = read_csv(pos[1]);
      if (t.cells.size() < 3) throw std::invalid_argument("three rows or more");
      if (op == "model") {
        if (!o.count("--y")) throw std::invalid_argument("--y COL");
        const size_t yc = column(t, o["--y"]);
        const auto xc = xcols(t, {yc});
        if (xc.empty()) throw std::invalid_argument("no numeric x columns");
        std::vector<double> y;
        std::vector<std::vector<double>> x(xc.size());
        std::vector<std::string> names;
        for (size_t c : xc) names.push_back(t.names[c]);
        for (size_t r = 0; r < t.cells.size(); ++r) {
          y.push_back(value(t, r, yc));
          for (size_t k = 0; k < xc.size(); ++k) x[k].push_back(value(t, r, xc[k]));
        }
        ModelOptions mo;
        if (o.count("--method")) mo.method = o["--method"];
        mo.stepwise = o.count("--stepwise"), mo.squares = o.count("--squares"), mo.products = o.count("--products");
        if (o.count("--components")) mo.components = std::stoi(o["--components"]);
        if (o.count("--folds")) mo.folds = std::stoi(o["--folds"]);
        if (o.count("--randomise")) mo.randomise = std::stoi(o["--randomise"]);
        if (o.count("--hidden")) mo.nn_hidden = std::stoi(o["--hidden"]);
        if (o.count("--seed")) mo.seed = std::stoull(o["--seed"]);
        const ModelResult r = fit_model(y, x, names, mo);
        std::printf("%s of %s on %zu columns, %zu rows\n%s\n", mo.method.c_str(), o["--y"].c_str(), names.size(), y.size(), r.equation.c_str());
        if (!r.terms.empty() && mo.method != "nn") {
          std::printf("%-28s %14s %12s\n", "term", "coefficient", "std. error");
          for (const auto& tm : r.terms) std::printf("%-28s %14.6g %12s\n", tm.name.c_str(), tm.coef, std::isfinite(tm.se) ? std::to_string(tm.se).c_str() : "-");
        }
        std::printf("R² %.4f · adjusted R² %s · RMSE %.4g", r.r2, std::isfinite(r.r2_adj) ? std::to_string(r.r2_adj).c_str() : "-", r.rmse);
        if (std::isfinite(r.q2)) std::printf(" · q² %.4f · RMSE(CV) %.4g", r.q2, r.rmse_cv);
        if (std::isfinite(r.lof)) std::printf(" · LOF %.4g", r.lof);
        std::printf("\n");
        if (!r.q2_by_components.empty()) {
          std::printf("q² by PLS components:");
          for (double v : r.q2_by_components) std::printf(" %.4f", v);
          std::printf("\n");
        }
        for (const auto& [e, l] : r.gfa_models) std::printf("  GFA model (LOF %.4g): %s\n", l, e.c_str());
        if (!r.cooks.empty()) {
          std::vector<size_t> idx(r.cooks.size());
          std::iota(idx.begin(), idx.end(), 0);
          std::sort(idx.begin(), idx.end(), [&](size_t a, size_t b) { return r.cooks[a] > r.cooks[b]; });
          std::printf("most influential rows (Cook's distance, leverage, standardised residual):");
          for (size_t k = 0; k < std::min<size_t>(3, idx.size()); ++k)
            std::printf(" %zu (%.3g, %.3g, %.3g)", idx[k] + 1, r.cooks[idx[k]], r.leverage[idx[k]], r.std_residual[idx[k]]);
          std::printf("\n");
        }
        for (const auto& n : r.notes) std::printf("%s\n", n.c_str());
        if (o.count("--json") && o["--json"] != "1") {
          std::ofstream(o["--json"]) << model_json(r).dump(1) << "\n";
          std::printf("wrote %s\n", o["--json"].c_str());
        }
        return 0;
      }
      if (op == "cluster" || op == "lda") {
        std::vector<size_t> ex;
        if (op == "lda") {
          if (!o.count("--class")) throw std::invalid_argument("--class COL");
          ex.push_back(column(t, o["--class"]));
        }
        const auto xc = xcols(t, ex);
        if (xc.empty()) throw std::invalid_argument("no numeric columns");
        std::vector<std::vector<double>> rows;
        for (size_t r = 0; r < t.cells.size(); ++r) {
          std::vector<double> v;
          for (size_t c : xc) v.push_back(value(t, r, c));
          rows.push_back(v);
        }
        auto label = [&](size_t r) { return t.cells[r].empty() ? std::to_string(r + 1) : t.cells[r][0]; };
        if (op == "cluster") {
          const std::string m = o.count("--method") ? o["--method"] : "ward";
          const int k = o.count("--k") ? std::stoi(o["--k"]) : 3;
          const RowClusters c = cluster_rows(rows, m, k, !o.count("--no-scale"));
          std::printf("%s clustering of %zu rows on %zu columns into %d (within sum of squares %.4g, %s columns)\n", m.c_str(), rows.size(), xc.size(), k, c.within_ss,
                      o.count("--no-scale") ? "raw" : "autoscaled");
          for (int q = 0; q < k; ++q) {
            std::printf("cluster %d:", q + 1);
            for (size_t r = 0; r < rows.size(); ++r) if (c.label[r] == q) std::printf(" %s", label(r).c_str());
            std::printf("\n");
          }
          return 0;
        }
        std::vector<std::string> cls;
        for (size_t r = 0; r < t.cells.size(); ++r) {
          if (ex[0] >= t.cells[r].size()) throw std::invalid_argument("row " + std::to_string(r + 1) + " has no " + t.names[ex[0]] + " (the class column)");
          cls.push_back(t.cells[r][ex[0]]);
        }
        const DiscriminantResult d = discriminant(rows, cls);
        std::printf("discriminant analysis of %zu rows, %zu classes, %zu columns: %.1f %% right on the rows, %.1f %% left out one at a time\n", rows.size(), d.classes.size(),
                    xc.size(), 100 * d.accuracy, 100 * d.accuracy_loo);
        for (size_t f = 0; f < d.directions.size(); ++f) {
          std::printf("function %zu (eigenvalue %.4g):", f + 1, d.eigenvalues[f]);
          for (size_t k = 0; k < xc.size(); ++k) std::printf(" %+.4g·%s", d.directions[f][k], t.names[xc[k]].c_str());
          std::printf("\n");
        }
        for (size_t r = 0; r < rows.size(); ++r)
          if (d.predicted_loo[r] != int(std::find(d.classes.begin(), d.classes.end(), cls[r]) - d.classes.begin()))
            std::printf("row %s (%s) taken for %s when left out\n", label(r).c_str(), cls[r].c_str(), d.classes[size_t(d.predicted_loo[r])].c_str());
        return 0;
      }
      throw std::invalid_argument("qspr: descriptors, model, cluster or lda");
    } catch (const std::exception& e) {
      std::fprintf(stderr, "caps qspr: %s\n", e.what());
      keep_cli_error(argc, argv, e);
      return 1;
    }
  }
  if (cmd == "habit") {
    try {
      if (pos.empty()) return usage();
      if (pos.size() > 1) throw std::invalid_argument("unexpected argument " + pos[1] + " (caps habit takes one CRYSTAL)");
      Trajectory t = open_file(pos[0], o.count("--topology") ? o["--topology"] : "");
      System s0 = t.frames() ? t.frame(0) : t.topology;
      const bool quiet = o.count("--quiet");
      HabitOptions ho;
      if (o.count("--max-index")) ho.max_index = std::stoi(o["--max-index"]);
      if (o.count("--min-d")) ho.min_d = std::stod(o["--min-d"]);
      if (o.count("--forms")) ho.forms = std::stoi(o["--forms"]);
      if (o.count("--vacuum")) ho.vacuum = std::stod(o["--vacuum"]);
      ho.attachment = o.count("--attachment");
      ho.surface = o.count("--surface-energy");
      // the structure as the force field describes it (shells, ionic bonds): with --ff alone too, so a force field that makes
      // ions the growth units shapes the BFDH analysis and the particle as well
      // the force field's charges and the names of the shells it adds (at their core's place, bonded to it) are kept for the
      // particle: its net charge is reported and the shells are folded into their cores before it is written
      const size_t atoms_read = s0.atoms.size();
      ForceField hff;
      std::set<std::string> shell_names;
      if (ho.attachment || ho.surface || o.count("--ff")) {
        hff = cli_forcefield(s0, o, quiet);
        for (const auto& b : s0.bonds) {
          const size_t lo = std::min(b.i, b.j), hi = std::max(b.i, b.j);
          if (hi >= atoms_read && lo < atoms_read && norm(s0.atoms[hi].pos - s0.atoms[lo].pos) < 1e-6) shell_names.insert(s0.atoms[hi].name);
        }
      }
      if (ho.attachment || ho.surface) {
        ho.energy.cutoff = o.count("--cutoff") ? std::stod(o["--cutoff"]) : 10.0;
        ho.energy.tail = false;
        ho.energy.electrostatics = EnergyOptions::Electrostatics::PME;   // slabs of ionic crystals need Ewald sums
        electrostatics(ho.energy, o);
        ho.energy.coulomb = !o.count("--no-coulomb");
        ho.field_for = [o](const System& x) mutable { System c = x; return cli_forcefield(c, o, true); };
      }
      if (!quiet)
        ho.progress = [&](const std::string& w, double f) { std::fprintf(stderr, "\r%-30s %3.0f %%", w.c_str(), 100 * f); return true; };
      const HabitResult r = crystal_habit(s0, ho);
      if (!quiet) std::fprintf(stderr, "\r%-40s\r", "");
      std::printf("%s (%d operations) · %d growth units per cell", r.space_group.c_str(), r.operations, r.units);
      if (std::isfinite(r.e_latt)) std::printf(" · bulk energy %.4f kcal/mol per unit", r.e_latt);
      std::printf("\n%5s %-10s %5s %9s %9s %12s %12s %10s %8s %8s %8s\n", "rank", "form", "faces", "d (Å)", "d_eff", "E_att", "E_slice", "γ (J/m²)", "BFDH", "growth", "equil.");
      auto num = [](double v, const char* f) { char b[32]; if (std::isfinite(v)) std::snprintf(b, sizeof b, f, v); else std::snprintf(b, sizeof b, "-"); return std::string(b); };
      for (const auto& f : r.faces) {
        if (f.rank > 20 && f.area_bfdh == 0 && f.area_growth == 0 && f.area_equilibrium == 0) continue;
        char hk[32];
        std::snprintf(hk, sizeof hk, "{%d %d %d}", f.hkl[0], f.hkl[1], f.hkl[2]);
        std::printf("%5d %-10s %5zu %9.4f %9.4f %12s %12s %10s %7.1f%% %7.1f%% %7.1f%%%s%s\n", f.rank, hk, f.members.size(), f.d, f.d_eff, num(f.e_att, "%.4f").c_str(),
                    num(f.e_slice, "%.4f").c_str(), num(f.gamma, "%.4f").c_str(), 100 * f.area_bfdh, 100 * f.area_growth, 100 * f.area_equilibrium,
                    f.note.empty() ? "" : "  ", f.note.c_str());
      }
      auto shape_line = [](const char* n, const WulffShape& w) {
        if (w.ok()) std::printf("%s shape: %zu facets, aspect ratio %.2f\n", n, w.facets.size(), w.aspect);
      };
      shape_line("BFDH", r.bfdh), shape_line("growth (attachment energy)", r.growth), shape_line("equilibrium (surface energy)", r.equilibrium);
      for (const auto& n : r.notes) std::printf("%s\n", n.c_str());
      const std::string model = o.count("--model") ? o["--model"] : (r.growth.ok() ? "growth" : "bfdh");
      if (o.count("--shape-obj")) {
        std::vector<Vec3> n;
        std::vector<double> d;
        habit_planes(s0, r, model, n, d);
        const WulffShape w = wulff_shape(n, d);
        if (!w.ok()) throw std::runtime_error("the " + model + " shape is not available");
        std::ofstream obj(o["--shape-obj"]);
        obj << "# CAPS " << model << " shape of " << pos[0] << "\n";
        for (const auto& v : w.vertices) obj << "v " << v[0] << " " << v[1] << " " << v[2] << "\n";
        for (const auto& f : w.facets) {
          obj << "f";
          for (int id : f) obj << " " << id + 1;
          obj << "\n";
        }
        std::printf("wrote %s (%s shape)\n", o["--shape-obj"].c_str(), model.c_str());
      }
      if (o.count("--particle")) {
        std::vector<Vec3> n;
        std::vector<double> d;
        habit_planes(s0, r, model, n, d);
        System sc = s0;   // the force field's charges on the atoms, so the carved particle carries them
        const bool charged = hff.charge.size() == sc.atoms.size() && !sc.atoms.empty();
        if (charged) {
          for (size_t i = 0; i < sc.atoms.size(); ++i) sc.atoms[i].charge = hff.charge[i];
          sc.has_charges = true;
        }
        System p = carve_particle(sc, n, d, std::stod(o["--particle"]));
        double q = 0;
        for (const auto& a : p.atoms) q += a.charge;
        // shells are the force field's particles, not atoms: each folds into its core (its charge with it) and the
        // structure file holds one atom per ion
        size_t shells = 0;
        if (!shell_names.empty()) {
          std::vector<char> rm(p.atoms.size(), 0);
          for (size_t i = 0; i < p.atoms.size(); ++i) rm[i] = shell_names.count(p.atoms[i].name) ? 1 : 0;
          for (const auto& b : p.bonds) {
            if (rm[b.i] && !rm[b.j]) p.atoms[b.j].charge += p.atoms[b.i].charge;
            else if (rm[b.j] && !rm[b.i]) p.atoms[b.i].charge += p.atoms[b.j].charge;
          }
          shells = size_t(std::count(rm.begin(), rm.end(), char(1)));
          if (shells) delete_atoms(p, rm);
          for (size_t i = 0; i < p.atoms.size(); ++i) p.atoms[i].id = int64_t(i + 1);
        }
        const std::string out = o.count("--out") ? o["--out"] : "particle.pdb";
        write_structure_file(p, out);
        std::map<std::string, int> by_element;
        for (const auto& a : p.atoms) by_element[a.element > 0 ? element(a.element).symbol : a.name]++;
        std::string formula;
        for (const auto& [e, c] : by_element) formula += (formula.empty() ? "" : " ") + e + std::to_string(c);
        std::printf("wrote %s: a %s-shaped particle of %zu atoms (%s), radius %s Å", out.c_str(), model.c_str(), p.atoms.size(), formula.c_str(), o["--particle"].c_str());
        if (shells) std::printf(", the %zu shells folded into their cores", shells);
        std::printf("\n");
        if (charged && std::abs(q) > 1e-3)
          std::printf("note: the particle's net charge is %+.3f e (the force field's charges): the cut kept or dropped single ions at its "
                      "faces, so it is not stoichiometric; neutralise it (remove surface ions) before a run that needs a neutral particle\n", q);
        else if (charged)
          std::printf("net charge %+.3f e\n", q);
      }
      if (o.count("--json") && o["--json"] != "1") {
        const Json j = habit_json(r);
        std::ofstream(o["--json"]) << j.dump(1) << "\n";
        std::printf("wrote %s\n", o["--json"].c_str());
      }
      return 0;
    } catch (const std::exception& e) {
      std::fprintf(stderr, "caps habit: %s\n", e.what());
      keep_cli_error(argc, argv, e);
      return 1;
    }
  }
  if (cmd == "fit") {
    // caps fit SPEC.json [--ff FF.json] [--relaxed | --conventional] [--iterations N] [--out FITTED.json] [--json REPORT.json]
    // SPEC: {"forcefield": "catlow1976-mgo.json", "structures": [{"file": "MgO.cif", "name": "MgO", "formula_units": 4}],
    //        "parameters": [{"kind": "pair", "types": ["Mg", "O"], "name": "A"}], "targets": [{"structure": 0, "kind": "C11", "value": 297}],
    //        "relaxed": true, "cutoff": 10, "electrostatics": "pme", "weights": {...}}
    try {
      if (pos.empty()) return usage();
      std::ifstream in(pos[0]);
      if (!in) throw std::runtime_error("cannot read " + pos[0]);
      std::stringstream ss;
      ss << in.rdbuf();
      const Json spec = Json::parse(ss.str());
      const std::filesystem::path base = std::filesystem::path(pos[0]).parent_path();
      auto near = [&](const std::string& f) {   // a file named in the spec: beside the spec unless found as given
        if (f.empty() || std::filesystem::exists(f)) return f;
        const auto q = base / f;
        return std::filesystem::exists(q) ? q.string() : f;
      };
      const bool quiet = o.count("--quiet");
      if (!o.count("--ff") && spec.has("forcefield")) o["--ff"] = near(spec["forcefield"].str());
      if (!o.count("--ff")) throw std::invalid_argument("fit: give the force field (--ff, or \"forcefield\" in the spec)");
      FitOptions fo;
      fit_options_from_json(spec, fo);
      if (o.count("--relaxed")) fo.relaxed = true;
      if (o.count("--conventional")) fo.relaxed = false;
      if (o.count("--iterations")) fo.iterations = std::stoi(o["--iterations"]);
      fo.energy.cutoff = o.count("--cutoff") ? std::stod(o["--cutoff"]) : spec.num("cutoff", 10.0);
      fo.energy.tail = false;
      // ionic crystals need Ewald sums: PME unless the spec says dsf
      fo.energy.electrostatics = spec.text("electrostatics", "pme") == "dsf" ? EnergyOptions::Electrostatics::DSF : EnergyOptions::Electrostatics::PME;
      electrostatics(fo.energy, o);
      fo.energy.coulomb = !o.count("--no-coulomb");
      FFDef def = load_forcefield(ff_path(o["--ff"]));
      std::set<std::string> shell_types;
      for (const auto& [core, sh] : def.shells) shell_types.insert(sh);
      if (!spec.has("structures") || !spec["structures"].is_array() || spec["structures"].items().empty())
        throw std::invalid_argument("fit: the spec lists no structures");
      for (const auto& js : spec["structures"].items()) {
        FitStructure fs;
        const std::string file = near(js.text("file"));
        Trajectory t = open_file(file, "");
        System s0 = t.frames() ? t.frame(0) : t.topology;
        if (!quiet) std::printf("%s: ", file.c_str());
        (void)cli_forcefield(s0, o, quiet);   // the structure as the force field describes it (shells added)
        fs.name = js.text("name", std::filesystem::path(file).stem().string());
        fs.cell = s0;
        fs.formula_units = js.num("formula_units", 1);
        if (js.has("fit_structure")) fs.fit_structure = js["fit_structure"].boolean();
        if (js.has("supercell") && js["supercell"].is_array() && js["supercell"].items().size() == 3)
          for (size_t q = 0; q < 3; ++q) fs.supercell[q] = int(js["supercell"].items()[q].number());
        for (const auto& a : s0.atoms) fs.shell.push_back(shell_types.count(a.name) ? 1 : 0);
        fs.field_for = [o](const System& sup) mutable { System c = sup; return cli_forcefield(c, o, true); };
        fo.structures.push_back(fs);
      }
      if (!quiet)
        fo.progress = [&](const std::string& w, double f) { std::fprintf(stderr, "\r%-30s %3.0f %%", w.c_str(), 100 * f); return true; };
      const FitResult r = fit_to_crystals(fo);
      if (!quiet) std::fprintf(stderr, "\r%-40s\r", "");
      std::printf("%s fitting: S %.6g -> %.6g after %d iterations (%s), %d evaluations, %.1f s\n", fo.relaxed ? "relaxed" : "conventional", r.ss_initial, r.ss_final,
                  r.iterations, r.stopped_by.c_str(), r.evaluations, r.seconds);
      std::printf("%-34s %14s %14s %12s\n", "parameter", "start", "fitted", "std. error");
      for (size_t k = 0; k < r.names.size(); ++k)
        std::printf("%-34s %14.6f %14.6f %12s\n", r.names[k].c_str(), r.initial[k], r.final[k],
                    std::isfinite(r.error[k]) ? std::to_string(r.error[k]).c_str() : "-");
      if (!r.observations.empty()) {
        std::printf("%-22s %12s %12s %12s  %s\n", "observable", "target", "start", "fitted", "unit");
        for (const auto& ob : r.observations)
          std::printf("%-22s %12.5f %12.5f %12.5f  %s\n", (r.structures[ob.structure].name + " " + ob.what).c_str(), ob.target, ob.initial, ob.final, ob.unit.c_str());
      }
      for (const auto& st : r.structures) {
        std::printf("%s (supercell %d x %d x %d, %zu atoms): a b c %.5f %.5f %.5f Å, angles %.3f %.3f %.3f° (observed %.5f %.5f %.5f, %.3f %.3f %.3f)\n", st.name.c_str(),
                    st.supercell[0], st.supercell[1], st.supercell[2], st.supercell_atoms, st.final[0], st.final[1], st.final[2], st.final[3], st.final[4], st.final[5],
                    st.observed[0], st.observed[1], st.observed[2], st.observed[3], st.observed[4], st.observed[5]);
        if (fo.relaxed) std::printf("  RMS shift of the atoms from observed %.4f Å (start %.4f)\n", st.shift_final, st.shift_initial);
        else std::printf("  largest force at the observed structure %.4g kcal/mol/Å (start %.4g); stress xx %.4g GPa\n", st.force_final, st.force_initial, st.stress_final[0]);
        std::printf("  lattice energy %.4f kcal/mol per formula unit\n", st.energy_final);
      }
      for (const auto& n : r.notes) std::printf("%s\n", n.c_str());
      if (o.count("--out")) {
        FFDef out = def;
        for (const auto& n : apply_fitted_terms(out, r.terms)) std::printf("%s\n", n.c_str());
        out.name += " (fitted)";
        save_forcefield(out, o["--out"]);
        std::printf("wrote %s\n", o["--out"].c_str());
      }
      if (o.count("--json") && o["--json"] != "1") {
        std::ofstream(o["--json"]) << fit_result_json(r).dump(1) << "\n";
        std::printf("wrote %s\n", o["--json"].c_str());
      }
      return 0;
    } catch (const std::exception& e) {
      std::fprintf(stderr, "caps fit: %s\n", e.what());
      keep_cli_error(argc, argv, e);
      return 1;
    }
  }
  if (cmd == "powder") {
    // caps powder FILE [--topology DATA] [--first N --last N --stride N] [--anode Cu | --wavelength Å] [--no-kalpha2] [--neutron]
    //             [--two-theta 5:90] [--step 0.02] [--profile tch|pseudo-voigt|gaussian|lorentzian] [--uvw U,V,W] [--xy X,Y] [--eta 0.5]
    //             [--size Å | --size-abc La,Lb,Lc] [--strain ε] [--biso B] [--preferred h,k,l --march r] [--zero °] [--monochromator 2θ°]
    //             [--polarised] [--peaks N] [--csv OUT.csv] [--json OUT.json]
    try {
      if (pos.empty()) return usage();
      OpenProgress rd;
      if (o.count("--first")) rd.frames.first = size_t(std::max(0L, std::stol(o["--first"])));
      if (o.count("--last") && std::stol(o["--last"]) >= 0) rd.frames.last = size_t(std::stol(o["--last"]));
      if (o.count("--stride")) rd.frames.stride = size_t(std::max(1L, std::stol(o["--stride"])));
      if (rd.frames.first > rd.frames.last) throw std::invalid_argument("--first " + std::to_string(rd.frames.first) + " is after --last " + std::to_string(rd.frames.last) + ": the selection is empty");
      Trajectory t = open_file(pos[0], o.count("--topology") ? o["--topology"] : "", rd);
      if (t.frames_read == 0 && !rd.frames.all()) t.select(rd.frames);
      PowderOptions po;
      auto pair = [&](const std::string& key, double& a, double& b) {
        const std::string v = o[key];
        const auto c = v.find_first_of(":,");
        if (c == std::string::npos) throw std::invalid_argument(key + " A:B");
        a = std::stod(v.substr(0, c)), b = std::stod(v.substr(c + 1));
      };
      auto list = [&](const std::string& key, size_t n) {
        std::vector<double> v;
        std::string s = o[key];
        for (auto& ch : s) if (ch == ':') ch = ',';
        size_t at = 0;
        while (at <= s.size()) {
          const size_t c = s.find(',', at);
          v.push_back(std::stod(s.substr(at, c == std::string::npos ? std::string::npos : c - at)));
          if (c == std::string::npos) break;
          at = c + 1;
        }
        if (v.size() != n) throw std::invalid_argument(key + ": " + std::to_string(n) + " numbers, comma-separated");
        return v;
      };
      if (o.count("--anode")) po.anode = o["--anode"];
      if (o.count("--wavelength")) po.wavelength = std::stod(o["--wavelength"]);
      if (o.count("--no-kalpha2")) po.kalpha2 = false;
      if (o.count("--neutron")) po.radiation = "neutron";
      if (o.count("--two-theta")) pair("--two-theta", po.two_theta_min, po.two_theta_max);
      if (o.count("--step")) po.step = std::stod(o["--step"]);
      if (o.count("--profile")) po.profile = o["--profile"];
      if (o.count("--uvw")) { const auto v = list("--uvw", 3); po.U = v[0], po.V = v[1], po.W = v[2]; }
      if (o.count("--xy")) { const auto v = list("--xy", 2); po.X = v[0], po.Y = v[1]; }
      if (o.count("--eta")) po.eta = std::stod(o["--eta"]);
      if (o.count("--size")) po.size = std::stod(o["--size"]);
      if (o.count("--size-abc")) { const auto v = list("--size-abc", 3); po.size_abc = {v[0], v[1], v[2]}; }
      if (o.count("--scherrer-k")) po.scherrer_k = std::stod(o["--scherrer-k"]);
      if (o.count("--strain")) po.strain = std::stod(o["--strain"]);
      if (o.count("--biso")) po.b_iso = std::stod(o["--biso"]);
      if (o.count("--preferred")) { const auto v = list("--preferred", 3); po.preferred = {int(v[0]), int(v[1]), int(v[2])}; }
      if (o.count("--march")) po.march_r = std::stod(o["--march"]);
      if (o.count("--zero")) po.zero_shift = std::stod(o["--zero"]);
      if (o.count("--monochromator")) po.monochromator_2theta = std::stod(o["--monochromator"]);
      if (o.count("--polarised")) po.polarised = true;
      if (o.count("--threads")) po.threads = std::stoi(o["--threads"]);
      const int show = o.count("--peaks") ? std::stoi(o["--peaks"]) : 20;
      std::vector<size_t> fr;
      for (size_t k = 0; k < std::max<size_t>(1, t.frames()); ++k) fr.push_back(k);
      const PowderReport r = powder_pattern(t, fr, po);
      for (const auto& n : r.notes) std::printf("%s\n", n.c_str());
      std::printf("%5s %4s %4s %4s %9s %9s %5s %10s %7s %7s\n", "#", "h", "k", "l", "d (Å)", "2θ (°)", "m", "|F|", "I", "FWHM");
      for (size_t i = 0; i < r.peaks.size() && int(i) < show; ++i) {
        const auto& p = r.peaks[i];
        std::printf("%5zu %4d %4d %4d %9.4f %9.3f %5d %10.3f %7.2f %7.3f\n", i + 1, p.h, p.k, p.l, p.d, p.two_theta, p.multiplicity, p.f, p.intensity, p.fwhm);
      }
      if (o.count("--csv")) {
        std::ofstream f(o["--csv"]);
        f << "two_theta_deg,q_inv_angstrom,intensity\n";
        for (size_t i = 0; i < r.two_theta.size(); ++i) f << r.two_theta[i] << "," << r.q[i] << "," << r.intensity[i] << "\n";
        const std::string pk = std::filesystem::path(o["--csv"]).replace_extension("").string() + "_peaks.csv";
        std::ofstream g(pk);
        g << "h,k,l,d_angstrom,two_theta_deg,q_inv_angstrom,multiplicity,F,intensity,fwhm_deg\n";
        for (const auto& p : r.peaks)
          g << p.h << "," << p.k << "," << p.l << "," << p.d << "," << p.two_theta << "," << p.q << "," << p.multiplicity << "," << p.f << "," << p.intensity << "," << p.fwhm << "\n";
        std::printf("wrote %s and %s\n", o["--csv"].c_str(), pk.c_str());
      }
      if (o.count("--json")) {
        Json j = Json::object(), pk = Json::array(), x = Json::array(), y = Json::array(), lam = Json::array(), nt = Json::array();
        for (const auto& p : r.peaks) {
          Json e = Json::object();
          e["h"] = Json(double(p.h)), e["k"] = Json(double(p.k)), e["l"] = Json(double(p.l)), e["d"] = Json(p.d), e["two_theta"] = Json(p.two_theta);
          e["q"] = Json(p.q), e["multiplicity"] = Json(double(p.multiplicity)), e["F"] = Json(p.f), e["intensity"] = Json(p.intensity), e["fwhm"] = Json(p.fwhm);
          pk.push_back(e);
        }
        for (size_t i = 0; i < r.two_theta.size(); ++i) x.push_back(Json(r.two_theta[i])), y.push_back(Json(r.intensity[i]));
        for (double v : r.wavelengths) lam.push_back(Json(v));
        for (const auto& n : r.notes) nt.push_back(Json(n));
        j["peaks"] = pk, j["two_theta"] = x, j["intensity"] = y, j["wavelengths"] = lam, j["notes"] = nt, j["frames"] = Json(double(r.frames));
        std::ofstream(o["--json"]) << j.dump(1) << "\n";
        std::printf("wrote %s\n", o["--json"].c_str());
      }
      return 0;
    } catch (const std::exception& e) {
      std::fprintf(stderr, "caps powder: %s\n", e.what());
      keep_cli_error(argc, argv, e);
      return 1;
    }
  }
  if (cmd == "frames") {
    // caps frames FILE OUT.lammpstrj|dcd|xyz|pdb|gro|trr [--topology DATA] [--first N] [--last N] [--stride N] [--timestep-fs 1]
    // A long trajectory thinned as it is read: only the frames kept are held and written.
    if (pos.size() < 2) return usage();
    try {
      OpenProgress rd;
      if (o.count("--first")) rd.frames.first = std::stoul(o["--first"]);
      if (o.count("--last")) rd.frames.last = std::stoul(o["--last"]);
      if (o.count("--stride")) rd.frames.stride = std::max<size_t>(1, std::stoul(o["--stride"]));
      if (rd.frames.first > rd.frames.last) throw std::invalid_argument("--first " + std::to_string(rd.frames.first) + " is after --last " + std::to_string(rd.frames.last) + ": the selection is empty");
      const Trajectory t = open_file(pos[0], o.count("--topology") ? o["--topology"] : "", rd);
      write_trajectory(t, pos[1], o.count("--timestep-fs") ? std::stod(o["--timestep-fs"]) : 1.0);
      std::printf("wrote %s: %zu frames (%zu read)\n", pos[1].c_str(), t.frames(), t.frames_read ? t.frames_read : t.frames());
      return 0;
    } catch (const std::exception& e) {
      std::fprintf(stderr, "caps frames: %s\n", e.what());
      keep_cli_error(argc, argv, e);
      return 1;
    }
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
      if (rd.frames.first > rd.frames.last) throw std::invalid_argument("--first " + std::to_string(rd.frames.first) + " is after --last " + std::to_string(rd.frames.last) + ": the selection is empty");
      Trajectory t = open_file(pos[0], o.count("--topology") ? o["--topology"] : "", rd);
      if (t.frames_read == 0 && !rd.frames.all()) t.select(rd.frames);   // a file read whole (one frame): the same frames
      AnalyzeOptions ao;
      if (o.count("--frame-ps")) ao.frame_ps = std::stod(o["--frame-ps"]) * double(rd.frames.stride);   // between the frames kept
      if (o.count("--timestep-fs")) ao.timestep_fs = std::stod(o["--timestep-fs"]);
      if (o.count("--blocks")) ao.blocks = std::stoi(o["--blocks"]);
      if (o.count("--rmax")) ao.rdf_rmax = std::stod(o["--rmax"]);
      if (o.count("--dr")) ao.rdf_dr = std::stod(o["--dr"]);
      if (!(ao.rdf_dr > 0)) throw std::invalid_argument("--dr: the g(r) bin width in Å, more than 0 (got " + o["--dr"] + ")");
      if (!(ao.rdf_rmax >= 0)) throw std::invalid_argument("--rmax: the g(r) range in Å, 0 or more (0: from the cell; got " + o["--rmax"] + ")");
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
        const auto [za, zb] = rdf_pair_elements(pr);
        ao.elem_a = za, ao.elem_b = zb;
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
      if (o.count("--ir-max-ps")) ao.ir_max_ps = std::stod(o["--ir-max-ps"]);
      if (o.count("--pattern")) ao.pattern = o["--pattern"];
      if (o.count("--pattern-bin")) ao.pattern_bin = std::stod(o["--pattern-bin"]);
      if (o.count("--hb-donors")) ao.hb_donors = o["--hb-donors"];
      if (o.count("--hb-acceptors")) ao.hb_acceptors = o["--hb-acceptors"];
      if (o.count("--hb-criterion")) ao.hb_criterion = o["--hb-criterion"];
      if (o.count("--hb-da")) ao.hb_da = std::stod(o["--hb-da"]);
      if (o.count("--hb-angle")) ao.hb_hda = ao.hb_dha = std::stod(o["--hb-angle"]);
      if (o.count("--hb-ha")) ao.hb_ha = std::stod(o["--hb-ha"]);
      if (o.count("--hb-between")) {   // 'QUERY A | QUERY B'
        const std::string v = o["--hb-between"];
        const auto bar = v.find('|');
        if (bar == std::string::npos) throw std::invalid_argument("--hb-between 'QUERY | QUERY'");
        ao.hb_group_a = v.substr(0, bar), ao.hb_group_b = v.substr(bar + 1);
      }
      if (o.count("--ir-max-wavenumber")) ao.ir_max_wavenumber = std::stod(o["--ir-max-wavenumber"]);
      const auto props = analyze(t, ids, ao);
      const auto fr = analysis_frames(t, ao);
      std::printf("%s · %zu frames used of %zu · %zu atoms\n", pos[0].c_str(), fr.size(), t.frames(), t.topology.atoms.size());
      cli_report(props, o);
      return 0;
    } catch (const std::exception& e) {
      std::fprintf(stderr, "caps analyze: %s\n", e.what());
      keep_cli_error(argc, argv, e);
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
      r.remove_pendant = o.count("--remove-pendant") > 0;
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
      // no --ff, and the file carries a force field (a CAPS data file names it): every state of the network is typed with
      // that field's library, as md, relax and tensile keep a file's own field, not the built-in GAFF or UFF in its place
      // (a file that carries the built-in GAFF stays on the built-in default, which types C and H with it)
      if (!o.count("--ff") && s.forcefield && !s.forcefield->name.empty() && s.forcefield->name.rfind("GAFF 1.81 (hydrocarbon subset", 0) != 0) {
        const bool uff = s.forcefield->name.rfind("UFF", 0) == 0;
        const std::string lib = uff ? std::string("uff") : library_of_field_name(s.forcefield->name);
        if (!lib.empty()) {
          o["--ff"] = lib;
          if (!o.count("--charges") && !uff) o["--charges"] = "auto";   // the field's own charges, for the new chemistry too
          if (!quiet) std::printf("force field: %s, the file's own, typed again after each reaction (give --ff for another)\n", s.forcefield->name.c_str());
        } else if (r.relax) {   // relaxing each state needs a field for the new bonds: not GAFF or UFF in the file's place unasked
          throw std::runtime_error("the file's force field (" + s.forcefield->name + ") is not a library CAPS can type the reacted network with; "
                                   "give --ff (oplsaa2024, pcff, gaff2 ...), or --no-relax for the topology only");
        }
      }
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
        if (!quiet) {
          char e[32] = "-";   // 0: no energy evaluated (relaxation off, or nothing reacted yet)
          if (c.energy != 0) std::snprintf(e, sizeof e, "%.1f", c.energy);
          std::printf("%6d %6d %6d %11.4f %9d %11.1f %12.1f %13s\n", c.cycle, c.reactions, c.total, c.conversion, c.clusters.clusters,
                      100 * c.clusters.largest_fraction, c.clusters.reduced_mw, e);
        }
        std::fflush(stdout);
        return true;
      };
      ReactReport rep;
      react(s, r, &rep);
      // the built-in default is GAFF only while every carbon is one it types (no C=C): the note names the one the network got
      std::shared_ptr<const ForceField> builtin;
      if (!r.retype && r.relax) {
        builtin = std::make_shared<ForceField>(default_forcefield(s));
        for (auto& n : rep.notes)
          if (n.rfind("force field during the run: built-in default", 0) == 0) n += "; the network written: " + builtin->name;
      }
      const std::string out = o["-o"];
      auto ends = [&](const char* e) { return out.size() > 4 && out.substr(out.size() - 4) == e; };
      if (ends(".pdb")) write_pdb(s, out);
      else if (ends(".xyz")) write_xyz(s, out);
      else if (r.retype) write_lammps_data_ff(s, *r.retype(s), r.energy, out);   // the network with its force field re-assigned
      else if (builtin) write_lammps_data_ff(s, *builtin, r.energy, out);
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
      char dtext[64] = "only one molecule: no intermolecular distance";
      if (std::isfinite(dmin)) std::snprintf(dtext, sizeof dtext, "smallest intermolecular distance %.4f Å", dmin);
      std::printf("%zu atoms, %d molecules, %s · %s · %d pairs closer than %.2f Å\n", s.atoms.size(), nm,
                  per ? "periodic (minimum image)" : "not periodic", dtext, close, tol);
      return close == 0 ? 0 : 1;
    }
    if (cmd == "bundle") {   // write a figure bundle: caps bundle FILE --steps S.json -o OUT.caps-bundle.zip [--include-input]
      const std::string top = o.count("--topology") ? o["--topology"] : "";
      const Trajectory t = open_file(pos[0], top);
      std::string text = o.count("--steps") ? o["--steps"] : "[]";
      if (!text.empty() && text[0] != '[' && text[0] != '{') {
        std::ifstream f(text);
        if (!f) throw std::runtime_error("cannot read " + text);
        text.assign(std::istreambuf_iterator<char>(f), {});
      }
      BundleOptions bo;
      bo.input = pos[0];
      bo.topology = top;
      bo.include_input = o.count("--include-input") > 0;
      bo.frame = o.count("--frame") ? std::stoi(o["--frame"]) : 0;
      const std::string out = o.count("-o") ? o["-o"] : "figure.caps-bundle.zip";
      bo.name = std::filesystem::path(out).stem().stem().string();
      const auto files = bundle_files(t, cli_pipeline(text), bo, Camera{}, RenderOptions{});
      write_bundle(out, files);
      for (const auto& f : files) std::printf("%-40s %8zu bytes  %s\n", f.name.c_str(), f.bytes.size(), sha256_hex(f.bytes).substr(0, 12).c_str());
      std::printf("wrote %s\n", out.c_str());
      return 0;
    }
    if (cmd == "pipeline") {   // visualize pipeline: steps from a JSON file (or --steps '[…]'), on one frame
      // the steps first: a pipeline YAML's source: topology: is the topology when --topology is not given
      std::string text = o.count("--steps") ? o["--steps"] : "[]";
      std::filesystem::path steps_dir;
      if (!text.empty() && text[0] != '[' && text[0] != '{') {
        std::ifstream f(text);
        if (!f) throw std::runtime_error("cannot read " + text);
        steps_dir = std::filesystem::path(text).parent_path();
        text.assign(std::istreambuf_iterator<char>(f), {});
      }
      std::string topo;
      Pipeline pl = cli_pipeline(text, nullptr, nullptr, &topo);
      std::string tp = o.count("--topology") ? o["--topology"] : "";
      if (tp.empty() && !topo.empty()) {
        tp = topo;
        if (std::filesystem::path(topo).is_relative() && !steps_dir.empty() && std::filesystem::exists(steps_dir / topo)) tp = (steps_dir / topo).string();
        if (!std::filesystem::exists(tp)) throw std::runtime_error("the topology " + topo + " named in " + o["--steps"] + " does not exist");
      }
      const Trajectory t = open_file(pos[0], tp);
      const int fr = o.count("--frame") ? std::stoi(o["--frame"]) : 0;
      if (fr < 0 || size_t(fr) >= t.frames()) throw std::runtime_error("frame out of range");
      const System f0 = t.frame(size_t(fr));
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
      // a crystal cell too small for its bonds holds one bond per pair of atoms (bonds.cpp): the count is short and the
      // chains look cut, so the crystal's own count (from the copies that hold every bond once) and the reason are given
      const std::string small = s.cell.valid() ? cell_fit_message(s) : "";
      std::string per_cell;
      if (!small.empty() && !s.bonds_from_file) {
        try {
          const auto n = copies_to_fit(s);
          const size_t k = size_t(n[0]) * size_t(n[1]) * size_t(n[2]);
          if (k * s.atoms.size() <= 4000)
            per_cell = "; the crystal has " + std::to_string(supercell(s, n[0], n[1], n[2]).bonds.size() / k) + " per cell, see the note";
        } catch (const std::exception&) {
        }
      }
      std::printf("file        %s\nformat      %s\natoms       %zu\nbonds       %zu (%s%s)\nmolecules   %d\n", pos[0].c_str(), s.source_format.c_str(), s.atoms.size(),
                  s.bonds.size(), s.bonds_from_file ? "from file" : "perceived", per_cell.c_str(), nm);
      std::printf("elements   ");
      for (const auto& [e, c] : counts) std::printf(" %s %d", e.c_str(), c);
      std::printf("\n");
      if (s.cell.valid())
        std::printf("cell        %.4f %.4f %.4f Å · volume %.1f Å³\ndensity     %.4f g/cm³\n", norm(s.cell.a), norm(s.cell.b), norm(s.cell.c), s.cell.volume(), s.density());
      std::printf("mass        %.3f g/mol\ncharge      %+.3e e\n", s.total_mass(), q);
      for (const auto& n : s.notes) std::printf("note        %s\n", n.c_str());
      if (!small.empty()) std::printf("note        %s\n", small.c_str());
      return 0;
    }
    if (cmd == "render") {
      if (!o.count("-o")) return usage();
      const std::string out = o["-o"];
      // every value checked before anything is drawn: an unknown one is refused, not drawn as the default
      auto pick = [&](const char* key, const char* def, std::initializer_list<const char*> ok) {
        const std::string v = o.count(key) ? o[key] : def;
        std::string list;
        for (const char* k : ok) {
          if (v == k) return v;
          list += (list.empty() ? "" : "|") + std::string(k);
        }
        throw std::invalid_argument(std::string(key) + ": unknown value '" + v + "' (" + list + ")");
      };
      std::string ext = std::filesystem::path(out).extension().string();
      for (char& ch : ext) ch = char(std::tolower(static_cast<unsigned char>(ch)));
      if (ext != ".png" && ext != ".svg") throw std::invalid_argument("-o: the picture is written as .png or .svg (got " + out + ")");
      const bool svg = ext == ".svg";
      Camera cam;
      RenderOptions opt;
      if (o.count("--yaw")) cam.yaw = std::stod(o["--yaw"]) * M_PI / 180;
      if (o.count("--pitch")) cam.pitch = std::stod(o["--pitch"]) * M_PI / 180;
      if (o.count("--zoom")) cam.zoom = std::stod(o["--zoom"]);
      cam.perspective = o.count("--perspective");
      if (o.count("--size")) {
        int w = 0, h = 0;
        char tail = 0;
        if (std::sscanf(o["--size"].c_str(), "%dx%d%c", &w, &h, &tail) != 2 || w <= 0 || h <= 0)
          throw std::invalid_argument("--size takes WIDTHxHEIGHT in pixels, e.g. 1280x800 (got " + o["--size"] + ")");
        // the renderer holds about 80 bytes per output pixel (2x2 supersampled colour, depth and atom id): 8192 a side and
        // 8K-poster area at most keeps it near 3 GB
        constexpr int side_max = 8192;
        constexpr long long area_max = 8192LL * 4608;
        if (w > side_max || h > side_max || (long long)w * h > area_max)
          throw std::invalid_argument("--size " + o["--size"] + " is too large: at most " + std::to_string(side_max) + " pixels a side and " +
                                      std::to_string(area_max / 1000000) + " megapixels in all (e.g. 7680x4320)");
        opt.width = w, opt.height = h;
      }
      const std::string bg = pick("--bg", "dark", {"dark", "white", "transparent"});
      const std::string colour_v = pick("--colour", "molecule", {"molecule", "element", "type", "distance"});
      const std::string style_v = pick("--style", "ball", {"ball", "space", "sticks", "noh", "backbone"});
      opt.background = bg == "white" ? Background::White : bg == "transparent" ? Background::Transparent : Background::Dark;
      const std::string& c = colour_v;
      opt.colour_by = c == "element" ? ColourBy::Element : c == "type" ? ColourBy::Type : c == "distance" ? ColourBy::Property : ColourBy::Molecule;
      if (opt.colour_by == ColourBy::Property) {
        const auto sh = molecule_shapes(s);
        const auto mol = s.molecules();
        for (size_t i = 0; i < s.atoms.size(); ++i) opt.property.push_back(norm(s.atoms[i].pos - sh[mol[i]].com));
      }
      const std::string& st = style_v;
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
      if (svg) {
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
      // --ff as relax, md and equilibrate take it (the file's own coefficients otherwise): a data file whose pair
      // coefficients sit in the LAMMPS input is then typed afresh, as its refusal without --ff advises
      const ForceField ff = o.count("--ff") ? cli_forcefield(s, o, true) : default_forcefield(s);
      std::printf("force field %s\n", ff.name.c_str());
      for (const auto& n : ff.notes) std::printf("  %s\n", n.c_str());
      std::map<std::string, std::string> why;
      for (size_t i = 0; i < ff.atom_type.size(); ++i) why[ff.atom_type[i]] = i < ff.why.size() ? ff.why[i] : "";
      for (const auto& [t, w] : why) std::printf("  %-3s %s\n", t.c_str(), w.c_str());
      EnergyOptions eo;
      if (o.count("--ff")) run_cutoff(eo, ff, o, true);   // the force field's own cut-off, as relax uses
      else if (o.count("--cutoff")) eo.cutoff = std::stod(o["--cutoff"]);
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
      const System s_in = s;   // shells or sites cli_forcefield adds are not written to a structure format (as_read)
      r.field = std::make_shared<ForceField>(cli_forcefield(s, o, o.count("--quiet") > 0));
      run_cutoff(r.energy, *r.field, o, o.count("--quiet") > 0);
      r.restraint_terms = cli_restraints(s, o, o.count("--quiet") > 0);
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
      if (ends(".pdb")) write_pdb(as_read(s_in, s), out);
      else if (ends(".gro")) write_gro(s, out);   // same atoms and order as a GROMACS topology read with --topology
      else if (ends(".xyz")) write_xyz(as_read(s_in, s), out);
      else if (ends("mol2")) write_mol2(as_read(s_in, s), out);
      else if (ends(".car")) write_car(as_read(s_in, s), out);
      else if (!lammps_data_name(out) && structure_format_name(out)) write_structure_file(as_read(s_in, s), out);
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
        // The cluster run has its own lengths (default 100 ps equilibration, 2000 ps sampled), as the C API's
        // lammps_equilibrate_ps / lammps_ps: --eq-ps / --ps are the local run's and are never raised here.
        q.method = to.method, q.temperature = to.temperature, q.dt = to.dt;
        q.equilibrate_ps = o.count("--lammps-eq-ps") ? std::stod(o["--lammps-eq-ps"]) : 100.0;
        q.ps = o.count("--lammps-ps") ? std::stod(o["--lammps-ps"]) : 2000.0;
        if (!(q.equilibrate_ps >= 0)) throw std::invalid_argument("--lammps-eq-ps: ps of equilibration in the LAMMPS input, 0 or more");
        if (!(q.ps > 0)) throw std::invalid_argument("--lammps-ps: ps sampled in the LAMMPS input, more than 0");
        q.sample_every = to.sample_every, q.corr_ps = to.corr_ps, q.axis = to.axis, q.slabs = to.slabs;
        q.swap_every = to.swap_every, q.dimensionality = to.dimensionality;
        q.thickness = to.dimensionality == 2 ? (to.thickness > 0 ? to.thickness : sheet_thickness(s)) : 0;
        q.kspace = to.energy.electrostatics == EnergyOptions::Electrostatics::PME;
        std::ofstream f(in, std::ios::app);
        f << thermal_lammps_run(q);
        std::printf("wrote %s and %s (%s, %g ps equilibration, %g ps sampled; --lammps-eq-ps / --lammps-ps set these)\n", in.c_str(), data.c_str(),
                    to.method.c_str(), q.equilibrate_ps, q.ps);
        if (!o.count("--ps")) return 0;   // only the inputs unless a local run is asked for too
        if (!quiet) std::printf("--ps given: running %g ps (after %g ps equilibration) here as well\n", to.ps, to.equilibrate_ps);
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
      // --resume: the first run's settings (beside the checkpoint) for every option not given again, before any is read
      if (o.count("--resume") && s.title.find("checkpoint step ") != std::string::npos) md_resume_settings(pos[0], s, o);
      DynamicsOptions d;
      if (o.count("--steps")) d.steps = std::stoll(o["--steps"]);
      if (d.steps < 0) throw std::invalid_argument("--steps must be 0 or more (got " + o["--steps"] + ")");
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
      d.restraints = cli_restraints(s, o, o.count("--quiet") > 0);
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
        const std::string settings = md_settings_text(o);   // the run's options, for --resume
        d.checkpoint = [&, cp, total, settings](const std::vector<double>& x, const std::vector<double>& v, const Cell& c, int64_t step) {
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
          const std::string sp = md_settings_path(cp);
          std::ofstream(sp + ".tmp") << settings;
          std::filesystem::rename(sp + ".tmp", sp);
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
        // every option of this run again (the physics: electrostatics, HMR, held molecule, restraints, precision ...), the
        // words quoted; the checkpoint is the input, so the original input's --topology and --frame do not apply to it
        auto ro = o;
        ro.erase("--topology");
        ro.erase("--frame");
        const std::string again = resume_line("caps md", cp, ro, {{"--resume", ""}, {"--steps", std::to_string(d.step_offset + d.steps)}});
        std::ofstream("resume.txt") << again << "\n";
        std::printf("stopped at step %lld (a stop was asked for): the state is in %s; continue with\n  %s\n", static_cast<long long>(e.step), cp.c_str(), again.c_str());
        live.done(false, "interrupted at step " + std::to_string(e.step) + "; resume: " + again);
        return 75;
      }
      std::filesystem::remove("resume.txt");   // finished: nothing to continue
      // the last step too (the core writes one only before the end), so the checkpoint is the finished run, not the last
      // multiple of N before it, and a longer --resume continues from the end instead of repeating the steps after it
      if (d.checkpoint_every > 0 && d.steps > 0 && s.velocities.size() == s.atoms.size()) {
        std::vector<double> fx(3 * s.atoms.size()), fv(3 * s.atoms.size());
        for (size_t i = 0; i < s.atoms.size(); ++i)
          for (int k = 0; k < 3; ++k) fx[3 * i + size_t(k)] = s.atoms[i].pos[k], fv[3 * i + size_t(k)] = s.velocities[i][k];
        d.checkpoint(fx, fv, s.cell, d.step_offset + d.steps);
      }
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
      if (s.total_mass() <= 0 && s.bonds.empty() && std::none_of(s.atoms.begin(), s.atoms.end(), [](const Atom& a) { return a.element > 0; }))
        throw std::runtime_error(pos[0] + " has no elements or masses (a LAMMPS dump without element or mass columns), so its backbones are unknown; "
                                 "give a topology with --topology DATA");
      const InternalDistances d = internal_distances(s);
      if (d.periodic > 0)
        std::printf("# %d of the backbones are infinite chains bonded through the cell, each cut once at a bond: their end-to-end distance is about one period, not a coil size\n", d.periodic);
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
      if (!(pp.time_scale > 0)) throw std::invalid_argument("--scale must be above 0 (got " + o["--scale"] + "); it multiplies every stage's length");
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
      e.md.restraints = cli_restraints(s, o, o.count("--quiet") > 0);
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
      if (!(dr > 0)) throw std::invalid_argument("--dr: the bin width in Å, more than 0 (got " + o["--dr"] + ")");
      if (!(rmax > 0)) throw std::invalid_argument("--rmax: the range in Å, more than 0 (got " + o["--rmax"] + ")");
      // g(r) is the pair density over the bulk density: without a cell there is no density, and every value would be 0
      if (!s.cell.valid())
        throw std::invalid_argument(pos[0] + " has no periodic cell, and g(r) is normalised by the number density of the cell; "
                                    "give a file with a cell (a LAMMPS data file, a PDB with CRYST1, a CIF ...)");
      int ea = 0, eb = 0;
      if (o.count("--pair")) {
        const auto [za, zb] = rdf_pair_elements(o["--pair"]);
        ea = za, eb = zb;
        auto present = [&](int z) { return z == 0 || std::any_of(s.atoms.begin(), s.atoms.end(), [z](const Atom& a) { return a.element == z; }); };
        for (int z : {ea, eb})
          if (!present(z)) throw std::invalid_argument("--pair " + o["--pair"] + ": the system has no " + std::string(element(z).symbol) + " atoms");
      }
      std::printf("# r(Å) g(r)\n");
      for (auto [r, g] : rdf(s, ea, eb, rmax, dr, o.count("--inter"))) std::printf("%.3f %.5f\n", r, g);
      return 0;
    }
    if (cmd == "pointgroup") {
      const double tol = o.count("--tolerance") ? std::stod(o["--tolerance"]) : o.count("--tol") ? std::stod(o["--tol"]) : 0.1;
      PointGroupReport r = point_group(s, tol);
      if (o.count("--symmetrize")) {
        if (!o.count("-o")) throw std::invalid_argument("--symmetrize needs -o OUT");
        const double before = r.max_deviation;
        r = symmetrize_molecule(s, tol);
        std::printf("symmetrised: the images met to %.4f Å, now to %.2e Å\n", before, r.max_deviation);
      }
      std::printf("point group %s · order %d · largest image miss %.6f Å (tolerance %.3f Å)\n", r.group.c_str(), r.order, r.max_deviation, tol);
      for (const auto& e : r.elements) {
        const std::string k = e.kind == "sigma" ? "σ" : e.kind == "i" ? "i" : e.kind + (e.n > 0 ? std::to_string(e.n) : std::string("∞"));
        if (e.kind == "i") std::printf("  i\n");
        else std::printf("  %-4s %s (%.5f, %.5f, %.5f)\n", k.c_str(), e.kind == "sigma" ? "normal" : "axis", e.axis[0], e.axis[1], e.axis[2]);
      }
      for (const auto& n : r.notes) std::printf("%s\n", n.c_str());
      if (o.count("--symmetrize")) {
        write_structure_file(s, o["-o"]);
        std::printf("wrote %s\n", o["-o"].c_str());
      }
      return 0;
    }
    if (cmd == "planes") {
      if (!o.count("--hkl")) throw std::invalid_argument("caps planes FILE --hkl h,k,l [--shift 0] [--max N] [--json OUT]");
      std::vector<int> hk;
      {
        std::string t = o["--hkl"];
        for (auto& ch : t) if (ch == ',') ch = ' ';
        std::istringstream in(t);
        for (int x; in >> x;) hk.push_back(x);
      }
      if (hk.size() != 3) throw std::invalid_argument("--hkl takes three integers, e.g. 1,1,-1");
      const double shift = o.count("--shift") ? std::stod(o["--shift"]) : 0.0;   // 0: through the origin, as the Studio and caps_miller
      const auto mp = miller_planes(s.cell, hk[0], hk[1], hk[2], shift, o.count("--max") ? std::stoi(o["--max"]) : 0);
      std::printf("%s  d = %.6f Å  normal (%.6f, %.6f, %.6f)\n", mp.label().c_str(), mp.d, mp.normal[0], mp.normal[1], mp.normal[2]);
      std::printf("view down the normal: yaw %.4f°, pitch %.4f°\n", mp.yaw * 180 / M_PI, mp.pitch * 180 / M_PI);
      std::printf("plane  offset Å  corners  area Å²\n");
      Json polys = Json::array();
      for (size_t q = 0; q < mp.polygons.size(); ++q) {
        const auto& poly = mp.polygons[q];
        Vec3 a2{0, 0, 0};
        for (size_t k = 0; k < poly.size(); ++k) a2 = a2 + cross(poly[k], poly[(k + 1) % poly.size()]);
        std::printf("%5zu %9.4f %8zu %8.3f\n", q + 1, mp.offsets[q], poly.size(), norm(a2) / 2);
        Json pj = Json::array();
        for (const auto& v : poly) { Json vj = Json::array(); for (int k = 0; k < 3; ++k) vj.push_back(v[size_t(k)]); pj.push_back(vj); }
        polys.push_back(pj);
      }
      if (mp.polygons.empty())   // (h00) at shift 0: the only planes are the cell's faces, which are not listed
        std::printf("no plane of %s crosses the inside of the cell at shift %g (planes on the cell's faces are not listed); "
                    "try --shift 0.5\n", mp.label().c_str(), shift);
      if (o.count("--json")) {
        Json j = Json::object();
        j["label"] = mp.label(), j["d"] = mp.d, j["yaw"] = mp.yaw, j["pitch"] = mp.pitch, j["polygons"] = polys;
        Json off = Json::array(), nrm = Json::array();
        for (double x : mp.offsets) off.push_back(x);
        for (int k = 0; k < 3; ++k) nrm.push_back(mp.normal[size_t(k)]);
        j["offsets"] = off, j["normal"] = nrm;
        std::ofstream f(o["--json"]);
        f << j.dump(1) << "\n";
        std::printf("wrote %s\n", o["--json"].c_str());
      }
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
      const bool had_ff = s.forcefield && s.forcefield->charge.size() == s.atoms.size();
      const ChemistryKey before = chemistry_key(s);
      const auto report = edit_script(s, script);
      for (const auto& line : report) std::printf("%s\n", line.c_str());
      // a geometry edit keeps the file's force field; one that changes atoms, elements, bonds or charges does not
      write_edited_structure(s, had_ff, chemistry_key(s) == before, pos[0], o["-o"], o, "the edit changed the atoms, elements, bonds or charges");
      std::printf("%zu atoms · %zu bonds · wrote %s\n", s.atoms.size(), s.bonds.size(), o["-o"].c_str());
      return 0;
    }
    if (cmd == "convert") {
      if (pos.size() < 2) return usage();
      const std::string out = pos[1];
      const bool had_ff = s.forcefield && s.forcefield->charge.size() == s.atoms.size();
      const ChemistryKey before = chemistry_key(s);
      // --supercell auto (the copies this cell needs for its bonds) or na,nb,nc; 1,1,1 writes the cell as it is
      bool asked = false;
      if (o.count("--supercell") && s.cell.valid()) {
        std::array<int, 3> n{1, 1, 1};
        std::string v = o["--supercell"];
        if (v == "auto") {
          n = copies_to_fit(s);
        } else {
          for (auto& c : v) if (c == 'x' || c == 'X' || c == ':') c = ',';
          if (std::sscanf(v.c_str(), "%d,%d,%d", &n[0], &n[1], &n[2]) != 3 || n[0] < 1 || n[1] < 1 || n[2] < 1)
            throw std::invalid_argument("--supercell: auto, or na,nb,nc (whole numbers of at least 1)");
        }
        asked = true;
        if (n != std::array<int, 3>{1, 1, 1}) {
          s = supercell(s, n[0], n[1], n[2]);
          std::printf("%d×%d×%d copies of the cell, %zu atoms, %zu bonds\n", n[0], n[1], n[2], s.atoms.size(), s.bonds.size());
        }
      } else if (o.count("--supercell")) {
        throw std::invalid_argument("--supercell: " + pos[0] + " has no cell");
      }
      // A cell too small for its bonds holds one bond per pair of atoms (bonds.cpp), so a polyethylene cell keeps 10 of its
      // 12 bonds: a format that carries bonds (LAMMPS data, mol2, PDB CONECT, SDF, CAR/MDF) would describe broken chains.
      // Coordinates alone (xyz, gro, cif, POSCAR) are right as they are.
      std::string low = out;
      for (char& c : low) c = char(std::tolower(static_cast<unsigned char>(c)));
      auto ends = [&](const char* e) { const std::string x(e); return low.size() > x.size() && low.compare(low.size() - x.size(), x.size(), x) == 0; };
      // the writers that carry bonds, as write_structure_file picks them (an unknown extension is left to its own refusal)
      const bool coordinates_only = ends(".xyz") || ends(".extxyz") || ends(".gro") || ends(".cif") || ends(".vasp") || ends(".poscar") || poscar_name(out);
      const bool with_bonds = !coordinates_only && (ends(".pdb") || ends(".mol2") || ends(".car") || ends(".sdf") || ends(".mol") || ends(".data") ||
                                                    ends(".lmp") || ends(".lammps") || !std::filesystem::path(out).has_extension());
      if (with_bonds && s.cell.valid()) {
        if (const std::string small = cell_fit_message(s); !small.empty()) {
          if (!asked) throw std::invalid_argument(small + " (caps convert " + pos[0] + " " + out + " --supercell auto writes those copies; --supercell 1,1,1 writes this cell anyway)");
          std::fprintf(stderr, "warning: %s; the bonds written here do not describe the crystal\n", small.c_str());
        }
      }
      // every writer by the extension (.car with its .mdf, POSCAR by name); a LAMMPS data file keeps the input's force
      // field unless --supercell made copies (the field is bound to the atoms it was read with)
      write_edited_structure(s, had_ff, chemistry_key(s) == before, pos[0], out, o, "--supercell made copies of the cell");
      std::printf("wrote %s\n", out.c_str());
      return 0;
    }
    return usage();
  } catch (const std::exception& e) {
    std::fprintf(stderr, "caps: %s\n", e.what());
    // a wrong argument is the user's to fix; anything else may be a bug in CAPS, so it is kept for caps report
    if (!dynamic_cast<const std::invalid_argument*>(&e)) {
      log_cli_error(argc, argv, e.what());
      std::fprintf(stderr, "caps: if this looks like a bug in CAPS, `caps report --email` sends the details to %s\n", kSupportEmail);
    }
    return 1;
  }
}
