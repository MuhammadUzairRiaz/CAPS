// caps — command-line front end over the same core the Studio uses.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "caps/analysis.hpp"
#include "caps/bench.hpp"
#include "caps/dynamics.hpp"
#include "caps/elements.hpp"
#include "caps/ffdef.hpp"
#include "caps/typing.hpp"
#include "caps/uff.hpp"
#include "caps/crystal.hpp"
#include "caps/nano.hpp"
#include "caps/properties.hpp"
#include "caps/equilibrate.hpp"
#include "caps/grow.hpp"
#include "caps/io.hpp"
#include "caps/mechanics.hpp"
#include "caps/molecule.hpp"
#include "caps/pack.hpp"
#include "caps/polymer.hpp"
#include "caps/react.hpp"
#include "caps/relax.hpp"
#include "caps/render.hpp"

using namespace caps;

namespace {

int usage() {
  std::cerr << "caps 0.1.0\n"
               "usage:\n"
               "  caps info    FILE [--topology DATA]\n"
               "  caps render  FILE -o OUT.png|OUT.svg [--topology DATA] [--frame N] [--size WxH]\n"
               "               [--bg dark|white|transparent] [--colour molecule|element|type|distance]\n"
               "               [--style ball|space|sticks|noh|backbone] [--yaw DEG] [--pitch DEG] [--zoom Z] [--no-cell]\n"
               "  caps shape   FILE [--topology DATA]          per-molecule Rg and shape\n"
               "  caps rdf     FILE [--topology DATA] [--rmax 12] [--dr 0.2] [--pair C-C] [--inter]\n"
               "  caps analyze FILE [--topology DATA] [--props density,rdf,sq,xray,neutron,rg,ree,cn,persistence,msd,diffusion,\n"
               "               relaxation,ced,delta,ffv,psd] [--first N --last N --stride N] [--frame-ps X | --timestep-fs 1]\n"
               "               [--pair C-C --inter] [--qmax 25 --dq 0.02 --qdirect 4] [--probe 0] [--grid 0.4] [--ff FF.json] [--json OUT] [--csv DIR]\n"
               "  caps elastic FILE [--topology DATA] [--method strain|fluct|fluct-run] [--configs N] [--strain 1e-4] [--temp T] [--ps 100] [--ff FF.json] [--json OUT]\n"
               "  caps tensile DATA -o OUT.data [--axis x] [--rate 1e-3] [--strain 0.2] [--temp 300] [--fixed-lateral] [--ff FF.json] [--csv DIR]\n"
               "  caps tg DATA -o OUT.data [--from 500 --to 200 --step 20 --ps 100] [--ff FF.json] [--csv DIR]   |   caps tg --fit TABLE.csv\n"
               "  caps convert FILE OUT.data|OUT.xyz|OUT.pdb [--topology DATA]\n"
               "  caps bench   [T1 T2 … | --all] [--repeats 3] [--quick] [--out DIR] [--samples DIR]   the built-in validation suite\n"
               "  caps build   SMILES -o OUT.mol2|OUT.pdb|OUT.xyz|OUT.data [--conformers 1] [--seed 1] [--ff FF.json] [--all]\n"
               "               a 3D molecule from SMILES; --ff cleans each conformer up with that force field (with typing rules)\n"
               "  caps surface CRYSTAL.cif -o OUT.data|mol2|pdb|xyz [--hkl 0,0,1] [--layers 3] [--termination 1] [--vacuum 15]\n"
               "               [--supercell 2,2] [--no-orthogonal] [--max-strain 2] [--passivate] [--list]   a slab (terminations listed)\n"
               "  caps interface CRYSTAL.cif|SLAB -o OUT --units SMILES[,…] [surface options] [--film 30] [--film-density 0.9]\n"
               "               [--chains N] [--dp 10] [--gap 1] [--vacuum 0] [--sequence …] [--ff FF]   a polymer film on a surface\n"
               "  caps nano    tube [--n 10 --m 10 --length 25 --finite] | sheet [--lx 20 --ly 20 --layers 1 --flake] |\n"
               "               particle CRYSTAL.cif [--shape sphere|cube|octahedron|cuboctahedron|fibre --radius 12 --length 20 --passivate]\n"
               "               [--units SMILES --chains 10 --dp 20 --density 0.9]   -o OUT   fillers, alone or in a polymer matrix\n"
               "  caps pull    FILE [--normal | --axis x|y|z] [--distance 10] [--rate 5] [--spring 10] [--temp 300] [--surface 1] [--csv OUT]\n"
               "               pull-out / debonding of a film from a held surface: interfacial shear strength, work of separation\n"
               "  caps blend   --components SMILES1,SMILES2 [--weights 0.5,0.5] [--chains 8] [--dp 20] [--density 0.5] [--slabs] -o OUT\n"
               "  caps grow    -o OUT.data|OUT.pdb|OUT.xyz [--chains 10] [--dp 8] [--density 0.5 | --box 33]\n"
               "               [--tacticity atactic|isotactic|syndiotactic] [--seed 1] [--trans] [--scale 1.0]\n"
               "               [--units '*CC(*)c1ccccc1,*CC(*)(C)C(=O)OC' --sequence homopolymer|alternating|block|random|gradient|pattern\n"
               "                --weights 0.7,0.3 --blocks 20,20 --pattern AAB --ff FF.json] [--auto-scale]   any repeat units (else polystyrene)\n"
               "  caps field   FILE [--topology DATA] [--forces OUT.txt]   GAFF types, terms, energy (and per-atom forces)\n"
               "  caps relax   FILE -o OUT.data|OUT.pdb|OUT.xyz [--method lbfgs|cg|sd|fire] [--ftol 0.5] [--iterations 5000]\n"
               "               [--density 1.05] [--step 0.06] [--box-relax] [--pressure 1] [--no-pushoff] [--cutoff 10]\n"
               "               [--no-coulomb] [--quiet]   (.data output carries the force field for LAMMPS)\n"
               "  caps md      FILE -o OUT.data [--steps 10000] [--dt 1] [--temp 300] [--thermostat bussi|langevin|none]\n"
               "               [--tau-t 100] [--barostat none|crescale|berendsen] [--pressure 1] [--tau-p 1000]\n"
               "               [--seed 1] [--new-velocities] [--thermo 100] [--dump TRAJ.lammpstrj --every 1000]\n"
               "               [--log thermo.csv] [--cutoff 10] [--skin 1.5] [--threads N] [--no-coulomb] [--quiet]\n"
               "  caps equilibrate FILE -o OUT.data [--protocol larsen21|annealing|pushoff|PROTOCOL.txt] [--print-protocol]\n"
               "               [--tfinal 300] [--tmax 600] [--pfinal 1] [--pmax 49346] [--scale 1] (atm, K; --scale shortens every stage)\n"
               "               [--cycles 3] [--tlow 300] [--thigh 600] [--ramp 50] [--hold 50]   (annealing, ps)\n"
               "               [--until-converged] [--block 20] [--max-blocks 20] [--dt 1] [--thermostat bussi|langevin]\n"
               "               [--barostat crescale|berendsen] [--seed 1] [--dump TRAJ.lammpstrj --every-ps 10] [--log thermo.csv] [--quiet]\n"
               "  caps chains  FILE [--topology DATA]          backbones and mean-square internal distances\n"
               "  caps pack    INPUT.inp [-o OUT] [--threads N] [--quiet]   packmol-style input (structure … end structure)\n"
               "  caps contacts FILE [--tol 2.0] [--no-pbc] [--molecule-size N]   smallest distance between atoms of different molecules\n"
               "  caps react   FILE -o OUT [--template NAME|FILE.txt ...] [--capture Å] [--per-cycle 10] [--cycles 100] [--target 1.0]\n"
               "               [--no-relax] [--md-ps 0] [--temp 300] [--seed 1] [--fa 2 --fb 4 --ratio 1] [--list-templates] [--quiet]\n"
               "  caps ff import-lt FILE.lt -o FF.json          convert a moltemplate force field into the CAPS format\n"
               "  caps ff import-dlf LIB/NAME.par -o FF.json    convert a DL_FIELD library (.par + .sf + .bci)\n"
               "  caps ff info FF.json                           types, rules, styles, references\n"
               "  caps ff type FILE --ff FF.json [--typing RULES.json] [-o TYPES.txt] [--explain]   assign atom types from SMARTS rules\n"
               "  caps ff apply FILE --ff FF.json [-o OUT.data [--lammps-input OUT.in]] [--overlay USER.json] [--types TYPES.txt] [--charges keep|types|gasteiger]\n"
               "               [--list] [-o OUT.data]   parameters for a structure whose atoms carry type names (or TYPES.txt)\n";
  return 2;
}

std::map<std::string, std::string> parse(int argc, char** argv, int from, std::vector<std::string>& pos) {
  std::map<std::string, std::string> o;
  for (int i = from; i < argc; ++i) {
    std::string a = argv[i];
    if (a.rfind("--", 0) == 0 || a == "-o") {
      const bool flag = a == "--no-cell" || a == "--inter" || a == "--perspective" || a == "--trans" || a == "--escalate" ||
                        a == "--box-relax" || a == "--no-pushoff" || a == "--no-coulomb" || a == "--quiet" || a == "--new-velocities" ||
                        a == "--until-converged" || a == "--print-protocol" || a == "--no-pbc" ||
                        a == "--no-relax" || a == "--list-templates" || a == "--list" || a == "--allow-missing" || a == "--no-tail" || a == "--explain" || a == "--names" || a == "--fixed-lateral" || a == "--volume" || a == "--quick" || a == "--all" || a == "--pme" || a == "--no-orthogonal" || a == "--passivate" || a == "--auto-scale" || a == "--finite" || a == "--flake" || a == "--normal" || a == "--slabs" ||
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
  Trajectory t = open_file(path, it == o.end() ? "" : it->second);
  size_t f = 0;
  if (auto fr = o.find("--frame"); fr != o.end()) f = std::stoul(fr->second);
  if (f >= t.frames()) throw std::runtime_error("frame " + std::to_string(f) + " out of range (file has " + std::to_string(t.frames()) + ")");
  System s = t.frame(f);
  if (!s.unwrapped) make_molecules_whole(s);
  s.notes.insert(s.notes.begin(), std::to_string(t.frames()) + " frame(s)");
  return s;
}

// The force field for commands that take --ff FF.json [--typing RULES] [--charges MODE]: typed by the force field's
// rules (or the atom names in the file); without --ff the built-in GAFF of C and H.
ForceField cli_forcefield(const System& s0, std::map<std::string, std::string>& o, bool quiet = false) {
  if (!o.count("--ff")) {
    const bool ch = std::all_of(s0.atoms.begin(), s0.atoms.end(), [](const Atom& a) { return a.element == 1 || a.element == 6; });
    if (!ch) {
      if (!quiet) std::printf("force field: UFF (elements beyond C and H); give --ff for another\n");
      return assign_uff(s0);
    }
    if (!quiet) std::printf("force field: the built-in GAFF (C and H); give --ff for another\n");
    return assign_gaff(s0);
  }
  if (is_uff(o["--ff"])) {
    UffOptions uo;
    uo.keep_charges = o.count("--charges") && o["--charges"] == "keep";
    uo.qeq = o.count("--charges") && o["--charges"] == "qeq";
    if (!quiet) std::printf("force field: UFF (every element; %s)\n", uo.qeq ? "QEq charges" : uo.keep_charges ? "charges from the file" : "no charges");
    return assign_uff(s0, uo);
  }
  FFDef def = load_forcefield(o["--ff"]);
  if (o.count("--typing")) load_typing(def, o["--typing"]);
  std::vector<std::string> types;
  if (!def.typing.empty()) {
    const TypingResult tr = assign_types(s0, def);
    if (tr.untyped) throw std::runtime_error(std::to_string(tr.untyped) + " atoms match no typing rule of " + def.name);
    types = tr.types;
  } else {
    for (const auto& a : s0.atoms) types.push_back(a.name);
  }
  ParamReport rep;
  ForceField ff = parameterize(s0, def, types, o.count("--charges") ? o["--charges"] : (s0.has_charges ? "keep" : "types"), &rep, false);
  if (!rep.missing.empty()) throw std::runtime_error(std::to_string(rep.missing.size()) + " parameters missing in " + def.name + " (caps ff apply lists them)");
  if (!quiet) std::printf("force field: %s\n", def.name.c_str());
  return ff;
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

int axis_of(const std::string& a) {
  if (a == "x" || a == "0") return 0;
  if (a == "y" || a == "1") return 1;
  if (a == "z" || a == "2") return 2;
  throw std::invalid_argument("axis must be x, y or z");
}

// --pme [--ewald-rtol 1e-5] [--pme-spacing 1.2] [--pme-order 4]: particle-mesh Ewald instead of damped shifted force
void electrostatics(EnergyOptions& e, std::map<std::string, std::string>& o) {
  if (o.count("--pme")) e.electrostatics = EnergyOptions::Electrostatics::PME;
  if (o.count("--ewald-rtol")) e.ewald_rtol = std::stod(o["--ewald-rtol"]);
  if (o.count("--pme-spacing")) e.pme_spacing = std::stod(o["--pme-spacing"]);
  if (o.count("--pme-order")) e.pme_order = std::stoi(o["--pme-order"]);
}

void save_structure(const System& s, const ForceField& ff, const EnergyOptions& e, const std::string& out) {
  auto ends = [&](const char* x) { return out.size() > 4 && out.substr(out.size() - 4) == x; };
  if (ends(".pdb")) write_pdb(s, out);
  else if (ends(".xyz")) write_xyz(s, out);
  else write_lammps_data_ff(s, ff, e, out);
}

}  // namespace

int main(int argc, char** argv) {
  if (argc < 3) return usage();
  const std::string cmd = argv[1];
  std::vector<std::string> pos;
  auto o = parse(argc, argv, 2, pos);
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
      g.auto_scale = o.count("--auto-scale");
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
        s = grow_chains(spec, g, &rep);
      } else
        s = grow(g, &rep);
      const std::string out = o["-o"];
      auto ends = [&](const char* e) { return out.size() > 4 && out.substr(out.size() - 4) == e; };
      if (ends(".pdb")) write_pdb(s, out);
      else if (ends(".xyz")) write_xyz(s, out);
      else write_lammps_data(s, out);
      for (const auto& n : rep.notes) std::printf("%s\n", n.c_str());
      std::printf("%zu atoms · %zu bonds · wrote %s\n", s.atoms.size(), s.bonds.size(), out.c_str());
      return 0;
    } catch (const std::exception& e) {
      std::fprintf(stderr, "caps grow: %s\n", e.what());
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
      else write_lammps_data(s, out);
      std::printf("%zu atoms · %zu bonds · wrote %s\n", s.atoms.size(), s.bonds.size(), out.c_str());
      return 0;
    } catch (const std::exception& e) {
      std::fprintf(stderr, "caps nano: %s\n", e.what());
      return 1;
    }
  }
  if (cmd == "blend") {
    // caps blend --components SMILES1,SMILES2[,…] [--weights 0.5,0.5] [--chains 8] [--dp 20] [--density 0.5] [--slabs] -o OUT
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
      BlendReport br;
      const System s = grow_blend(comps, bo, &br);
      for (const auto& n : br.notes) std::printf("%s\n", n.c_str());
      const std::string out = o["-o"];
      auto ends = [&](const char* e) { return out.size() > 4 && out.substr(out.size() - 4) == e; };
      if (ends(".pdb")) write_pdb(s, out);
      else if (ends(".xyz")) write_xyz(s, out);
      else if (ends("mol2")) write_mol2(s, out);
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
      else write_lammps_data(slab, out);
      std::printf("wrote %s\n", out.c_str());
      return 0;
    } catch (const std::exception& e) {
      std::fprintf(stderr, "caps surface: %s\n", e.what());
      return 1;
    }
  }
  if (cmd == "bench") {
    try {
      BenchOptions b;
      // the shipped data: $CAPS_HOME, the working directory, or the source tree this was built from
      auto find = [](const std::string& rel, const std::string& probe) {
        std::vector<std::string> roots;
        if (const char* h = std::getenv("CAPS_HOME")) roots.push_back(h);
        roots.push_back(".");
        roots.push_back(CAPS_SOURCE_ROOT);
        for (const auto& r : roots)
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
  if (cmd == "ff") {
    try {
      if (pos.empty()) return usage();
      const std::string sub = pos[0];
      if (sub == "import-dlf") {
        // caps ff import-dlf LIB/PCFF.par [--sf LIB/PCFF.sf] [--bci LIB/PCFF.bci] -o FF.json (sf / bci found next to the .par by default)
        if (pos.size() < 2 || !o.count("-o")) return usage();
        const std::string par = pos[1], stem = par.substr(0, par.size() - 4);
        auto exists = [](const std::string& f) { return std::ifstream(f).good(); };
        const std::string sf = o.count("--sf") ? o["--sf"] : exists(stem + ".sf") ? stem + ".sf" : "";
        const std::string bci = o.count("--bci") ? o["--bci"] : exists(stem + ".bci") ? stem + ".bci" : "";
        FFDef ff = import_dlfield(par, sf, bci);
        save_forcefield(ff, o["-o"]);
        std::printf("%s: %zu types, %zu pair, %zu bond, %zu angle, %zu dihedral, %zu improper rules, %zu bond increments · wrote %s\n", ff.name.c_str(),
                    ff.types.size(), ff.pairs.size(), ff.bonds.size(), ff.angles.size(), ff.dihedrals.size(), ff.impropers.size(),
                    ff.bond_increments.size(), o["-o"].c_str());
        for (const auto& n : ff.notes) std::printf("  %s\n", n.c_str());
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
        const FFDef ff = load_forcefield(pos[1]);
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
        FFDef ff = load_forcefield(o["--ff"]);
        if (o.count("--typing"))
          for (std::stringstream ts(o["--typing"]); ts.good();) {   // several files: comma-separated, later ones on top
            std::string f;
            std::getline(ts, f, ',');
            if (!f.empty()) load_typing(ff, f);
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
        if (uff) {
          UffOptions uo;
          uo.keep_charges = o.count("--charges") && o["--charges"] == "keep";
          uo.qeq = o.count("--charges") && o["--charges"] == "qeq";
          f = assign_uff(s, uo);
          for (const auto& n : f.notes) std::printf("%s\n", n.c_str());
        } else {
          FFDef ff = load_forcefield(o["--ff"]);
          if (o.count("--overlay")) merge_forcefield(ff, load_forcefield(o["--overlay"]));
          if (o.count("--typing"))
            for (std::stringstream ts(o["--typing"]); ts.good();) {   // several files: comma-separated, later ones on top
              std::string f;
              std::getline(ts, f, ',');
              if (!f.empty()) load_typing(ff, f);
            }
          std::vector<std::string> types;
          // without --types, atoms are typed by the force field's rules when it has them (else the file's atom names)
          const bool auto_type = !ff.typing.empty() && !o.count("--types") && !o.count("--names");
          if (auto_type) {
            const TypingResult r = assign_types(s, ff);
            for (const auto& n : r.notes) std::printf("note: %s\n", n.c_str());
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
          const std::string charges = o.count("--charges") ? o["--charges"] : (s.has_charges ? "keep" : "types");
          ParamReport rep;
          f = parameterize(s, ff, types, charges, &rep, o.count("--allow-missing"));
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
                        pf.form == 1 ? "buck" : "morse", pf.a, pf.b, pf.c);
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
          write_lammps_data_ff(s, f, eo, o["-o"]);
          std::printf("wrote %s\n", o["-o"].c_str());
          if (o.count("--lammps-input")) {   // the LAMMPS commands that reproduce this energy with the data file
            write_lammps_input(s, f, eo, o["-o"], o["--lammps-input"], o.count("--fix-mol") ? std::stoll(o["--fix-mol"]) : 0);
            std::printf("wrote %s\n", o["--lammps-input"].c_str());
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
    // caps tensile DATA -o OUT.data [--axis x] [--rate 1e-3 (1/ps)] [--strain 0.2] [--temp 300] [--dt 1] [--fixed-lateral]
    //              [--pressure 1] [--tau-p 1000] [--fit 0.02] [--seed 1] [--new-velocities] [--ff ...] [--dump T.lammpstrj --every N]
    try {
      if (pos.empty() || !o.count("-o")) return usage();
      System s = load(pos[0], o);
      auto ff = std::make_shared<ForceField>(cli_forcefield(s, o));
      TensileOptions to;
      to.field = ff;
      if (o.count("--cutoff")) to.energy.cutoff = std::stod(o["--cutoff"]);
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
      const TensileResult r = run_tensile(s, to);
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
  if (cmd == "tg") {
    // caps tg DATA -o OUT.data [--from 500 --to 200 --step 20 --ps 100] [--dt 1] [--pressure 1] [--barostat crescale|berendsen]
    //         [--seed 1] [--new-velocities] [--ff ...] [--json] [--csv DIR]
    // caps tg --fit TABLE.csv       two-line fit of an existing table (T, density or specific volume)
    try {
      if (o.count("--fit")) {
        std::ifstream f(o["--fit"]);
        if (!f) throw std::runtime_error("cannot read " + o["--fit"]);
        std::vector<double> T, y;
        for (std::string line; std::getline(f, line);) {
          if (line.empty() || line[0] == '#') continue;
          for (char& c : line) if (c == ',' || c == ';' || c == '\t') c = ' ';
          std::istringstream ls(line);
          double a, b;
          if (ls >> a >> b) { T.push_back(a); y.push_back(b); }
        }
        const bool density = !y.empty() && y[0] > 0.5 && y[0] < 3 && !o.count("--volume");
        CoolingResult r;
        for (size_t i = 0; i < T.size(); ++i) {
          CoolingPoint p;
          p.temperature = T[i];
          p.density = density ? y[i] : 1 / y[i];
          p.specific_volume = density ? 1 / y[i] : y[i];
          r.points.push_back(p);
        }
        std::vector<double> v;
        for (const auto& p : r.points) v.push_back(p.specific_volume);
        r.fit = fit_bilinear(T, v);
        r.method = "two-line fit (free hinge) of specific volume against T from " + o["--fit"] + " (" + std::to_string(T.size()) + " points" +
                   (density ? ", densities inverted" : "") + "), error by bootstrap";
        if (!r.fit.note.empty()) r.notes.push_back(r.fit.note);
        cli_report(cooling_properties(r), o);
        return 0;
      }
      if (pos.empty() || !o.count("-o")) return usage();
      System s = load(pos[0], o);
      auto ff = std::make_shared<ForceField>(cli_forcefield(s, o));
      CoolingOptions co;
      co.field = ff;
      if (o.count("--cutoff")) co.energy.cutoff = std::stod(o["--cutoff"]);
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
      int last = -1;
      co.progress = [&](const ThermoRow& r, int k, int n) {
        if (!quiet && k != last) {
          if (k < 0) std::printf("equilibrating at %.1f K\n", r.target_temperature);
          else std::printf("step %d of %d: %.1f K\n", k + 1, n, r.target_temperature);
          last = k;
        }
        return true;
      };
      const CoolingResult r = run_cooling(s, co);
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
  if (cmd == "analyze") {
    // caps analyze FILE [--topology DATA] --props density,rdf,... [--first --last --stride] [--frame-ps | --timestep-fs]
    //              [--pair C-C] [--inter] [--rmax] [--dr] [--qmax] [--dq] [--qdirect] [--probe] [--grid] [--ff FF.json [--typing R] [--charges]]
    //              [--json OUT.json] [--csv DIR]
    try {
      if (pos.empty()) return usage();
      Trajectory t = open_file(pos[0], o.count("--topology") ? o["--topology"] : "");
      AnalyzeOptions ao;
      if (o.count("--first")) ao.first = std::stol(o["--first"]);
      if (o.count("--last")) ao.last = std::stol(o["--last"]);
      if (o.count("--stride")) ao.stride = std::stol(o["--stride"]);
      if (o.count("--frame-ps")) ao.frame_ps = std::stod(o["--frame-ps"]);
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
      if (has("ced") || has("delta") || has("cij_fluct") || has("adhesion")) {
        System s0 = t.frame(0);
        if (!s0.unwrapped) make_molecules_whole(s0);
        ff = cli_forcefield(s0, o);
        ao.ff = &ff;
        if (o.count("--cutoff")) ao.energy.cutoff = std::stod(o["--cutoff"]);
      }
      if (o.count("--temp")) ao.temperature = std::stod(o["--temp"]);
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
      else if (r.relax) write_lammps_data_ff(s, default_forcefield(s), r.energy, out);
      else write_lammps_data(s, out);
      for (const auto& n : rep.notes) std::printf("%s\n", n.c_str());
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
      Evaluator ev(ff, eo);
      std::vector<double> x, f;
      for (const auto& a : s.atoms) x.insert(x.end(), a.pos.begin(), a.pos.end());
      const EnergyTerms e = ev.compute(x, s.cell, f);
      std::printf("energy (kcal/mol): bond %.2f  angle %.2f  dihedral %.2f  improper %.2f  vdW %.2f  Coulomb %.2f  total %.2f\n", e.bond, e.angle,
                  e.dihedral, e.improper, e.vdw, e.coulomb, e.total());
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
      if (o.count("--cutoff")) r.energy.cutoff = std::stod(o["--cutoff"]);
      r.energy.coulomb = !o.count("--no-coulomb");
      electrostatics(r.energy, o);
      if (o.count("--ff")) r.field = std::make_shared<ForceField>(cli_forcefield(s, o));
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
      else if (ends(".xyz")) write_xyz(s, out);
      else if (ends("mol2")) write_mol2(s, out);
      else write_lammps_data_ff(s, r.field ? *r.field : default_forcefield(s), r.energy, out);
      std::printf("%s\n", rep.field.c_str());
      for (const auto& n : rep.notes) std::printf("%s\n", n.c_str());
      for (const auto& st : rep.stages)
        std::printf("  %-34s %5d it  E %12.1f  |F|max %8.3f  (%s)\n", st.name.c_str(), st.iterations, st.energy, st.fmax, st.stopped_by.c_str());
      std::printf("wrote %s\n", out.c_str());
      return rep.converged ? 0 : 3;
    }
    if (cmd == "md") {
      if (!o.count("-o")) return usage();
      DynamicsOptions d;
      if (o.count("--steps")) d.steps = std::stoll(o["--steps"]);
      if (o.count("--dt")) d.dt = std::stod(o["--dt"]);
      if (o.count("--temp")) d.temperature = std::stod(o["--temp"]);
      if (o.count("--thermostat")) d.thermostat = thermostat_from_string(o["--thermostat"]);
      if (o.count("--tau-t")) d.tau_t = std::stod(o["--tau-t"]);
      if (o.count("--barostat")) d.barostat = barostat_from_string(o["--barostat"]);
      if (o.count("--pressure")) d.pressure = std::stod(o["--pressure"]);
      if (o.count("--tau-p")) d.tau_p = std::stod(o["--tau-p"]);
      if (o.count("--seed")) d.seed = std::stoull(o["--seed"]);
      if (o.count("--thermo")) d.thermo_every = std::stoi(o["--thermo"]);
      if (o.count("--every")) d.frame_every = std::stoi(o["--every"]);
      if (o.count("--cutoff")) d.energy.cutoff = std::stod(o["--cutoff"]);
      if (o.count("--skin")) d.energy.skin = std::stod(o["--skin"]);
      if (o.count("--threads")) d.energy.threads = std::stoi(o["--threads"]);
      d.energy.coulomb = !o.count("--no-coulomb");
      electrostatics(d.energy, o);
      d.new_velocities = o.count("--new-velocities");
      if (o.count("--ff")) d.field = std::make_shared<ForceField>(cli_forcefield(s, o));
      if (o.count("--fix-mol")) {   // hold one molecule in place (the substrate of an interface is molecule 1)
        const int64_t mm = std::stoll(o["--fix-mol"]);
        d.fixed.assign(s.atoms.size(), 0);
        for (size_t i = 0; i < s.atoms.size(); ++i) d.fixed[i] = s.atoms[i].mol == mm;
      }
      const bool quiet = o.count("--quiet");
      Trajectory traj;
      traj.topology = s;
      if (o.count("--dump"))
        d.frame = [&](const std::vector<double>& x, const Cell& c, int64_t step) {
          std::vector<Vec3> p(x.size() / 3);
          for (size_t i = 0; i < p.size(); ++i) p[i] = {x[3 * i], x[3 * i + 1], x[3 * i + 2]};
          traj.positions.push_back(std::move(p));
          traj.cells.push_back(c);
          traj.timesteps.push_back(step);
        };
      else
        d.frame_every = 0;
      if (!quiet) std::printf("%10s %9s %9s %13s %13s %13s %10s %8s\n", "step", "time/ps", "T/K", "Epot", "Etotal", "conserved", "P/atm", "ρ/g·cm⁻³");
      d.progress = [&](const ThermoRow& r) {
        if (!quiet)
          std::printf("%10lld %9.3f %9.2f %13.3f %13.3f %13.3f %10.1f %8.4f\n", static_cast<long long>(r.step), r.time_ps, r.temperature, r.potential,
                      r.total, r.conserved, r.pressure, r.density);
        return true;
      };
      DynamicsReport rep;
      run_dynamics(s, d, &rep);
      const std::string out = o["-o"];
      auto ends = [&](const char* e) { return out.size() > 4 && out.substr(out.size() - 4) == e; };
      if (ends(".pdb")) write_pdb(s, out);
      else if (ends(".xyz")) write_xyz(s, out);
      else if (ends("mol2")) write_mol2(s, out);
      else write_lammps_data_ff(s, d.field ? *d.field : default_forcefield(s), d.energy, out);
      if (o.count("--dump")) write_lammps_dump(traj, o["--dump"]);
      if (o.count("--log")) {
        std::ofstream lg(o["--log"]);
        lg << "step,time_ps,temperature_K,potential,kinetic,total,conserved,pressure_atm,volume_A3,density_g_cm3\n";
        char b[256];
        for (const auto& r : rep.thermo) {
          std::snprintf(b, sizeof b, "%lld,%.4f,%.3f,%.4f,%.4f,%.4f,%.4f,%.2f,%.2f,%.5f\n", static_cast<long long>(r.step), r.time_ps, r.temperature,
                        r.potential, r.kinetic, r.total, r.conserved, r.pressure, r.volume, r.density);
          lg << b;
        }
      }
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
      if (o.count("--seed")) e.md.seed = std::stoull(o["--seed"]);
      if (o.count("--cutoff")) e.md.energy.cutoff = std::stod(o["--cutoff"]);
      e.until_converged = o.count("--until-converged");
      if (o.count("--block")) e.block_ps = std::stod(o["--block"]);
      if (o.count("--max-blocks")) e.max_blocks = std::stoi(o["--max-blocks"]);
      if (o.count("--every-ps")) e.frame_ps = std::stod(o["--every-ps"]);
      const bool quiet = o.count("--quiet");
      Trajectory traj;
      traj.topology = s;
      if (o.count("--dump"))
        e.frame = [&](const std::vector<double>& x, const Cell& c, int64_t step) {
          std::vector<Vec3> p(x.size() / 3);
          for (size_t i = 0; i < p.size(); ++i) p[i] = {x[3 * i], x[3 * i + 1], x[3 * i + 2]};
          traj.positions.push_back(std::move(p));
          traj.cells.push_back(c);
          traj.timesteps.push_back(step);
        };
      std::printf("%s# total %.1f ps\n", protocol_text(e.stages).c_str(), protocol_ps(e.stages) * 1.0);
      int last_stage = -1;
      double last_print = -1e9;
      e.progress = [&](int st, int n, const std::string& label, const ThermoRow& r) {
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
      equilibrate(s, e, &rep);
      const std::string out = o["-o"];
      auto ends = [&](const char* x) { return out.size() > 4 && out.substr(out.size() - 4) == x; };
      if (ends(".pdb")) write_pdb(s, out);
      else if (ends(".xyz")) write_xyz(s, out);
      else write_lammps_data_ff(s, default_forcefield(s), e.md.energy, out);
      if (o.count("--dump")) write_lammps_dump(traj, o["--dump"]);
      if (o.count("--log")) {
        std::ofstream lg(o["--log"]);
        lg << "step,time_ps,target_K,temperature_K,potential,kinetic,total,pressure_atm,volume_A3,density_g_cm3\n";
        char b[256];
        for (const auto& r : rep.thermo) {
          std::snprintf(b, sizeof b, "%lld,%.4f,%.2f,%.3f,%.4f,%.4f,%.4f,%.2f,%.2f,%.5f\n", static_cast<long long>(r.step), r.time_ps, r.target_temperature,
                        r.temperature, r.potential, r.kinetic, r.total, r.pressure, r.volume, r.density);
          lg << b;
        }
      }
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
    if (cmd == "convert") {
      if (pos.size() < 2) return usage();
      const std::string out = pos[1];
      if (out.size() > 4 && out.substr(out.size() - 4) == ".xyz") write_xyz(s, out);
      else if (out.size() > 4 && out.substr(out.size() - 4) == ".pdb") write_pdb(s, out);
      else if (out.size() > 5 && out.substr(out.size() - 5) == ".mol2") write_mol2(s, out);
      else write_lammps_data(s, out);
      std::printf("wrote %s\n", out.c_str());
      return 0;
    }
    return usage();
  } catch (const std::exception& e) {
    std::fprintf(stderr, "caps: %s\n", e.what());
    return 1;
  }
}
