// The DFT surface & adsorption workbench as commands (caps/dft_commands.hpp).
#include "caps/dft_commands.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>

#include "caps/adsorb_dft.hpp"
#include "caps/crystal.hpp"
#include "caps/elements.hpp"
#include "caps/dft_analysis.hpp"
#include "caps/io.hpp"
#include "caps/molecule.hpp"
#include "caps/provenance.hpp"
#include "caps/slab2d.hpp"
#include "caps/trajectory.hpp"
#include "caps/vasp_out.hpp"
#include "caps/vasp_set.hpp"

namespace caps {

namespace fs = std::filesystem;

namespace {

// ---------------------------------------------------------------- options
struct Opt { const char* name; const char* def; const char* why; };
struct Cmd { const char* name; const char* usage; const char* about; std::vector<Opt> opts; std::vector<const char*> examples; };

const std::vector<Cmd>& table() {
  static const std::vector<Cmd> t = {
      {"sheet", "caps sheet [--preset NAME | --from FILE] -o OUT",
       "One 2D sheet: a preset from data/sheets/sheets.json (layer lists: element, stacking position, height), or ONE complete layer cut out of a bulk, stacked or MAX-type structure.",
       {{"preset", "Ti3C2", "a layer list from data/sheets/sheets.json (Ti3C2: Materials Project mp-1094034)"},
        {"a", "", "in-plane lattice constant (Å); heights scale with it — a start for the DFT cell relaxation"},
        {"from", "", "a CIF / POSCAR / any structure to cut one layer out of"},
        {"formula", "", "the layer's composition up to a multiple (Ti3C2), so incomplete layers are skipped"},
        {"remove", "", "elements taken out first (a MAX phase's A element: Al)"},
        {"index", "0", "which complete layer, from the bottom"},
        {"gap", "2.0", "Å along z that separates two layers"},
        {"o", "sheet.vasp", "the output (POSCAR, CIF, PDB, XYZ, LAMMPS data … by extension)"}},
       {"caps sheet --preset Ti3C2 -o Ti3C2.vasp", "caps sheet --from Ti3AlC2.cif --remove Al --formula Ti3C2 -o sheet.vasp"}},
      {"terminate", "caps terminate SHEET [--top O:0.5,OH:0.25,F:0.25] -o OUT",
       "Terminations on both faces of a sheet: fcc/hcp hollows, top or bridge sites; height from the bond length h = √(d² − r²); mixed faces by fractions with a seed; the bottom face the same as the top unless --janus.",
       {{"top", "O", "species and fractions on the top face (O, OH, F, Cl, Br, I, S, Se, Te, NH, H; data/sheets/terminations.json)"},
        {"bottom", "", "the bottom face with --janus"},
        {"janus", "false", "different faces (a net dipole: the dipole correction is switched on)"},
        {"site", "fcc", "fcc (above the third layer), hcp (above the second), top or bridge"},
        {"site_bottom", "", "the bottom face's site with --janus"},
        {"supercell", "1x1", "in-plane repeat; a mixed face needs several sites (3x3)"},
        {"vacuum", "20", "Å between periodic images, centred along z (< 15 Å warns, < 10 Å fails)"},
        {"seed", "1", "the random arrangement of a mixed face (NumPy default_rng: the same on every platform)"},
        {"bond", "", "overrides, e.g. O=2.05,OH=2.2,OH.tail=0.97 (Å)"},
        {"order", "Ti,C,N,O,F,H", "species order of the POSCAR (the POTCAR order follows it)"},
        {"surface", "", "the outer layer's element (default: the most common heavy element there)"},
        {"o", "slab.vasp", "the output"},
        {"vasp_set", "", "also write a VASP set (01_cell … 04_static) into this folder"}},
       {"caps terminate Ti3C2.vasp --top OH -o Ti3C2_OH.vasp", "caps terminate Ti3C2.vasp --top O:0.5,OH:0.25,F:0.25 --supercell 3x3 --seed 1 -o mixed.vasp --vasp-set structures/Ti3C2_mixed"}},
      {"sites", "caps sites SLAB [--element Ti]",
       "The surface sites of both faces of a slab: three-fold hollows classified fcc (above the third layer) or hcp (above the second), top and bridge, with the in-plane distance r a termination's height comes from (h = √(d² − r²)).",
       {{"element", "", "the outer layer's element (default: the most common heavy element there)"},
        {"species", "", "a termination (O, OH, F …): the height it takes on each site"}}, {"caps sites Ti3C2.vasp --species O --json"}},
      {"validate", "caps validate FILE [FILE …] [--expect Ti3C2O2]",
       "PASS/FAIL report of a 2D slab: labels, composition × N, one complete sheet and its layers, fragments, floating atoms, vacuum, close contacts, isolated atoms, termination distances and coordination, O–H bonds, site classes, top/bottom asymmetry. Exit code 1 on FAIL.",
       {{"expect", "", "expected formula per cell, e.g. Ti3C2O2 (× N)"}, {"core", "", "the bare layer's formula, e.g. Ti3C2"}},
       {"caps validate structures/*/POSCAR --core Ti3C2", "caps validate CONTCAR --expect Ti3C2O2H2 --json"}},
      {"adsorb-dft", "caps adsorb-dft SLAB [--smiles S | --molecule FILE] --out ROOT",
       "A molecule on a slab for DFT: anchors by SMARTS, anchor-down / parallel / upright × azimuths (symmetry-equivalent ones marked), lowered to the contact limits, image distance checked; the trio slab/ molecule/ complex_* written with one settings set, then audited independently.",
       {{"smiles", "C/C=C/CCC(C#N)C/C=C/C", "the molecule (stereo kept), embedded and cleaned up with UFF"},
        {"molecule", "", "a structure file instead of SMILES"},
        {"from_relaxed", "false", "the slab is a relaxed CONTCAR: repeated by --supercell and re-centred"},
        {"supercell", "1x1", "repeat of SLAB (with --from-relaxed: of the relaxed cell)"},
        {"modes", "", "anchor:nitrile,anchor:vinyl,parallel,upright (default: every anchor found + parallel)"},
        {"azimuths", "0,90,180,270", "degrees about z"},
        {"skip_equivalent", "false", "leave out azimuths the outer layer's symmetry makes equivalent"},
        {"prescreen", "0", "also the N lowest configurations of the force-field Adsorption Locator (UFF, simulated annealing) as complexes"},
        {"prescreen_steps", "20000", "Monte Carlo steps per annealing cycle"}, {"prescreen_cycles", "3", "annealing cycles"},
        {"dmin", "2.3", "Å, closest molecule–slab contact (any atoms)"},
        {"dheavy", "3.0", "Å, closest heavy-atom contact"},
        {"vacuum", "20", "Å"},
        {"force", "false", "write complexes that fail their checks"},
        {"out", "adsorption", "the set's folder"},
        {"pp_dir", "", "your licensed PAW directory (POTCARs written into every folder)"},
        {"profile", "slurm-workspace", "cluster profile (data/dft/clusters.json)"},
        {"kdens", "45", "Å: N_i = ceil(kdens / |a_i|)"}},
       {"caps adsorb-dft structures/Ti3C2_OH/03_relax/CONTCAR --from-relaxed --supercell 5x5 --out adsorption/Ti3C2_OH_NBR"}},
      {"vasp-set", "caps vasp-set STRUCTURE --out DIR",
       "One VASP case folder: POSCAR (species grouped), INCAR.cell / .relax / .static with a reason per tag, KPOINTS, job.slurm (stages, resume, backup, watchdog, storage), make_potcar.sh, POTCAR from your PAW directory, the validation report.",
       {{"out", "", "the case folder"}, {"label", "", "SYSTEM and the job name (default: the folder name)"},
        {"stages", "slab", "slab (01_cell 02_cell2 03_relax 04_static), fixed (03_relax 04_static) or a list"},
        {"static_out", "pot", "04_static outputs: none (energy + DOS), pot (+ LOCPOT), all (+ CHGCAR, AECCAR)"},
        {"encut", "", "eV (default 1.3 × the largest ENMAX of the POTCARs; 520 without them)"},
        {"kdens", "45", "Å"}, {"ivdw", "12", "D3-BJ"}, {"ispin", "1", "non-magnetic"},
        {"dipole", "auto", "auto (when the faces differ), on, off — energy-only, no LDIPOL"},
        {"gamma", "false", "Γ only (a molecule in a box)"},
        {"pp_dir", "", "your licensed PAW directory"}, {"profile", "slurm-workspace", "cluster profile"}},
       {"caps vasp-set Ti3C2_O.vasp --out structures/Ti3C2_O --pp-dir $PP_DIR"}},
      {"vasp-conv", "caps vasp-conv setup|collect CASE",
       "ENCUT and k-mesh convergence single points from a case folder (conv_encut_*, conv_k_*), and the choice: the cheapest setting for which it AND every more expensive one pass (ENCUT: in-plane stress ≤ 1 kB and max force ≤ 0.01 eV/Å; k: E0 ≤ 1 meV/atom and stress ≤ 1 kB).",
       {{"encuts", "400,450,520,600,700,800", "eV"}, {"kmeshes", "6,9,12,15,18", "N × N × 1"}},
       {"caps vasp-conv setup structures/Ti3C2_O", "caps vasp-conv collect structures/Ti3C2_O --json"}},
      {"vasp-scan", "caps vasp-scan make CASE LO HI STEP | fit CASE",
       "In-plane lattice-constant scan (fixed-cell relaxations, ascan_<a>) and its fit: the parabola's minimum and Birch–Murnaghan (3rd order) on the cell volume. A cross-check of 01_cell, or its replacement on VASP < 6.4.",
       {}, {"caps vasp-scan make structures/Ti3C2_O 2.98 3.18 0.04", "caps vasp-scan fit structures/Ti3C2_O"}},
      {"vasp-derived", "caps vasp-derived charge|cdd|freq|aimd CASE",
       "Sets made from a finished run: charge (the complex's single point with CHGCAR, AECCAR, LOCPOT), cdd (slab part and molecule part at the frozen complex geometry: same cell, k-mesh, ENCUT → same FFT grid), freq (partial Hessian, IBRION 5, the nitrile group freed), aimd (CSVR NVT, H as deuterium, 1 fs, segments of NSW steps).",
       {{"pp_dir", "", "PAW directory (cdd, aimd)"}, {"free", "", "freq: 0-based atoms to free (default: the nitrile N, C, α-C and its H)"},
        {"temperature", "300", "aimd, K"}, {"nsw", "1000", "aimd steps (fs) per segment"}, {"segments", "6", "aimd: 1 equilibration + 5 production"},
        {"out", "", "aimd folder"}, {"profile", "slurm-workspace", "cluster profile"}},
       {"caps vasp-derived charge adsorption/Ti3C2_OH_NBR/complex_nitrile_down_az0", "caps vasp-derived aimd complex_X/03_relax/CONTCAR --out aimd/Ti3C2_OH_best"}},
      {"vasp-jobs", "caps vasp-jobs submit|update|reset|cleanup|store …",
       "Job tools: submit (folders whose last stage is not DONE and that are not queued; --run calls sbatch), update (rewrite job.slurm from the profile, keeping name, ranks, time, stages, hang), reset CASE STAGE (renames the real stage folder <stage>_bad), cleanup ROOT (only finished stages; --keep-best), store ROOT PATH (stage folders to long-term storage with links, writes .store). Dry runs unless --yes.",
       {{"yes", "false", "really do it (cleanup, store)"}, {"keep_best", "false", "cleanup: also the charge files of all but the best complex of each set"},
        {"run", "false", "submit: call sbatch"}, {"profile", "slurm-workspace", "update: the profile to write"}},
       {"caps vasp-jobs submit adsorption/Ti3C2_OH_NBR/*", "caps vasp-jobs cleanup . --keep-best --yes"}},
      {"vasp-check", "caps vasp-check CASE|STAGE …", "Finished? converged? SCF failures, ionic steps, final E0, max force, magnetisation, warnings, the lattice change of the cell stages (02_cell2 < 0.002 Å), the CONTCAR re-validated.",
       {{"complex", "false", "an adsorption complex (no slab re-validation)"}}, {"caps vasp-check structures/Ti3C2_O"}},
      {"vasp-progress", "caps vasp-progress CASE …", "Stages done and running, live files in the workspace, steps, last energies, max force against 0.01 eV/Å, time per step, ETA from a log(fmax) fit.", {}, {"caps vasp-progress structures/*"}},
      {"vasp-health", "caps vasp-health SET", "Every complex of an adsorption set against E_ref = E(slab) + E(molecule): ENERGY, JUMP, DIPOLE, SCF, GEOM, HANG — each with what to do. Exit code 1 when anything is flagged.", {}, {"caps vasp-health adsorption/Ti3C2_OH_NBR"}},
      {"vasp-bind", "caps vasp-bind SET", "E_bind = E0(complex) − E0(slab) − E0(molecule) from 04_static, eV and kJ/mol, ranked; settings compared across the trio; binding_energies.csv.", {}, {"caps vasp-bind adsorption/Ti3C2_OH_NBR"}},
      {"vasp-analyze", "caps vasp-analyze geom|wf|dos|cdd|bader|freq|md|summary …",
       "geom COMPLEX… (heights, anchor distance, H-bonds, contacts, tilt) · wf STAGE [--reference STAGE] (work function, Δφ) · dos STAGE (grouped DOS, DOS(E_F), HOMO/LUMO-like peaks) · cdd COMPLEX (Δρ, CHGCAR_diff, ΔQ) · bader SLAB_POSCAR COMPLEX_POSCAR ACF.dat [POTCAR] · freq COMPLEX (ν(C≡N) shift) · md AIMD (stability verdict) · summary MAIN (results_summary.csv).",
       {{"reference", "", "wf: the slab-alone stage for Δφ"}, {"equil_ps", "1.0", "md: equilibration cut"}, {"dt_fs", "1.0", "md: time step"}},
       {"caps vasp-analyze geom adsorption/Ti3C2_OH_NBR/complex_*", "caps vasp-analyze summary ~/mxene"}}};
  return t;
}

const Cmd& cmd_of(const std::string& c) {
  for (const auto& x : table()) if (c == x.name) return x;
  throw std::invalid_argument("unknown command " + c);
}

std::string S(const Json& a, const std::string& k, const std::string& def = "") {
  if (!a.has(k)) return def;
  const Json& v = a[k];
  if (v.kind() == Json::String) return v.str();
  if (v.is_number()) { char b[64]; std::snprintf(b, sizeof b, "%.10g", v.number()); return b; }
  if (v.kind() == Json::Bool) return v.boolean() ? "true" : "false";
  return def;
}
double N(const Json& a, const std::string& k, double def) {
  if (!a.has(k)) return def;
  if (a[k].is_number()) return a[k].number();
  const std::string s = S(a, k);
  if (s.empty()) return def;
  try { return std::stod(s); } catch (...) { throw std::invalid_argument("--" + k + " needs a number, not '" + s + "'"); }
}
bool B(const Json& a, const std::string& k, bool def = false) {
  if (!a.has(k)) return def;
  if (a[k].kind() == Json::Bool) return a[k].boolean();
  const std::string s = S(a, k);
  return s.empty() || s == "true" || s == "1" || s == "yes" || s == "on";
}
std::vector<std::string> L(const Json& a, const std::string& k) {
  std::vector<std::string> out;
  if (!a.has(k)) return out;
  if (a[k].kind() == Json::Array) { for (const auto& x : a[k].items()) out.push_back(x.kind() == Json::String ? x.str() : S(Json(Json::object()), "", "")); return out; }
  std::stringstream ss(S(a, k));
  for (std::string t; std::getline(ss, t, ',');) { t.erase(0, t.find_first_not_of(' ')); t.erase(t.find_last_not_of(' ') + 1); if (!t.empty()) out.push_back(t); }
  return out;
}
std::vector<std::string> inputs(const Json& a) {
  std::vector<std::string> out;
  if (a.has("inputs")) for (const auto& x : a["inputs"].items()) out.push_back(x.str());
  return out;
}
std::pair<int, int> supercell(const Json& a, const char* k = "supercell") {
  const std::string s = S(a, k, "1x1");
  const auto x = s.find_first_of("x×, ");
  if (x == std::string::npos) { const int n = std::stoi(s); return {n, n}; }
  return {std::stoi(s.substr(0, x)), std::stoi(s.substr(x + (s.compare(x, 2, "×") == 0 ? 2 : 1)))};
}

bool poscar_like(const std::string& p) {
  const std::string n = fs::path(p).filename().string();
  auto ends = [&](const char* e) { const std::string x(e); return n.size() >= x.size() && n.compare(n.size() - x.size(), x.size(), x) == 0; };
  return n.rfind("POSCAR", 0) == 0 || n.rfind("CONTCAR", 0) == 0 || ends(".vasp") || ends(".poscar");
}
System read_any(const std::string& p, std::map<std::string, std::string>* bad = nullptr, std::string* comment = nullptr) {
  if (poscar_like(p)) return read_vasp_poscar(p, bad, comment);
  const std::string n = fs::path(p).extension().string();
  if (n == ".cif") return read_cif(p);
  const Trajectory t = open_file(p);
  return t.frame(0);
}
void write_any(const System& s, const std::string& p, const std::string& comment) {
  const std::string e = fs::path(p).extension().string();
  if (poscar_like(p)) write_vasp_poscar(s, p, comment);
  else if (e == ".cif") write_cif(s, p);
  else if (e == ".pdb") write_pdb(s, p);
  else if (e == ".xyz") write_xyz(s, p);
  else if (e == ".mol2") write_mol2(s, p);
  else if (e == ".gro") write_gro(s, p);
  else write_lammps_data(s, p);
}

void provenance(const std::string& path, const std::string& engine, const std::string& summary, const KeyValues& params, const std::string& rng = "") {
  Manifest m;
  ProvStep st;
  st.engine = engine;
  st.summary = summary;
  st.params = params;
  st.rng = rng;
  st.time = now_iso();
  m.steps.push_back(st);
  try { write_manifest(m, path); } catch (...) {}
}

std::vector<std::string> order_of(const Json& a) {
  auto o = L(a, "order");
  return o.empty() ? std::vector<std::string>{"Ti", "C", "N", "O", "F", "H"} : o;
}

VaspSetOptions set_options(const Json& a, const std::string& data_dir) {
  VaspSetOptions o;
  o.data_dir = data_dir;
  o.pp_dir = S(a, "pp_dir");
  o.profile = S(a, "profile", "slurm-workspace");
  o.settings.kdens = N(a, "kdens", 45);
  o.settings.encut = N(a, "encut", 0);
  o.settings.ivdw = int(N(a, "ivdw", 12));
  o.settings.ispin = int(N(a, "ispin", 1));
  o.settings.static_out = S(a, "static_out", "pot");
  const std::string dip = S(a, "dipole", "auto");
  o.settings.dipole = dip == "on" || dip == "true" ? 1 : dip == "off" || dip == "false" ? 0 : -1;
  o.gamma_only = B(a, "gamma");
  return o;
}

Json set_result(const VaspSetResult& r) {
  Json j = Json::object();
  j["dir"] = r.dir;
  Json f = Json::array();
  for (const auto& x : r.files) f.push_back(x);
  j["files"] = f;
  j["encut"] = r.encut;
  Json k = Json::array();
  for (int x : r.kpoints) k.push_back(x);
  j["kpoints"] = k;
  j["ranks"] = r.parallel.ranks, j["kpar"] = r.parallel.kpar, j["ncore"] = r.parallel.ncore, j["time"] = r.parallel.time;
  Json n = Json::array();
  for (const auto& x : r.notes) n.push_back(x);
  j["notes"] = n;
  return j;
}

Json plan_json(const std::vector<ToolAction>& plan) {
  Json a = Json::array();
  long long tot = 0;
  for (const auto& p : plan) { Json x = Json::object(); x["path"] = p.path, x["what"] = p.what, x["bytes"] = double(p.bytes); a.push_back(x); tot += p.bytes; }
  Json r = Json::object();
  r["items"] = a;
  r["bytes"] = double(tot);
  return r;
}

}  // namespace

const std::vector<std::string>& dft_commands() {
  static std::vector<std::string> c;
  if (c.empty()) for (const auto& x : table()) c.push_back(x.name);
  return c;
}
bool dft_is_switch(const std::string& c, std::string name) {
  while (!name.empty() && name[0] == '-') name.erase(0, 1);
  std::replace(name.begin(), name.end(), '-', '_');
  if (name == "json" || name == "help" || name == "force" || name == "yes" || name == "run") return true;
  for (const auto& x : table())
    if (c == x.name) for (const auto& o : x.opts) if (name == o.name) return std::string(o.def) == "false";
  return false;
}

bool is_dft_command(const std::string& c) { return std::find(dft_commands().begin(), dft_commands().end(), c) != dft_commands().end(); }

std::string dft_help(const std::string& c) {
  const Cmd& x = cmd_of(c);
  std::ostringstream o;
  o << "usage: " << x.usage << "\n\n" << x.about << "\n";
  if (!x.opts.empty()) {
    o << "\noptions (default · why):\n";
    for (const auto& p : x.opts) {
      std::string n = p.name;
      std::replace(n.begin(), n.end(), '_', '-');
      o << "  " << (n.size() == 1 ? "-" : "--") << n << (*p.def ? std::string("  [") + p.def + "]" : "") << "\n      " << p.why << "\n";
    }
  }
  o << "  --json\n      machine-readable output\n";
  if (!x.examples.empty()) { o << "\nexamples:\n"; for (const auto* e : x.examples) o << "  " << e << "\n"; }
  return o.str();
}

Json dft_args_from_cli(const std::vector<std::string>& pos, const std::vector<std::pair<std::string, std::string>>& flags) {
  Json a = Json::object();
  Json in = Json::array();
  for (const auto& p : pos) in.push_back(p);
  a["inputs"] = in;
  for (auto [k, v] : flags) {
    while (!k.empty() && k[0] == '-') k.erase(0, 1);
    std::replace(k.begin(), k.end(), '-', '_');
    if (k == "output") k = "o";
    a[k] = v;
  }
  return a;
}

std::string dft_command_line(const std::string& c, const Json& a) {
  std::string s = "caps " + c;
  auto q = [](const std::string& x) { return x.find_first_of(" \t'\"$*?") == std::string::npos && !x.empty() ? x : "'" + x + "'"; };
  if (a.has("inputs")) for (const auto& x : a["inputs"].items()) s += " " + q(x.str());
  for (const auto& [k, v] : a.members()) {
    if (k == "inputs" || k == "json") continue;
    std::string n = k;
    std::replace(n.begin(), n.end(), '_', '-');
    if (v.kind() == Json::Bool) { if (v.boolean()) s += " --" + n; continue; }
    std::string val = v.kind() == Json::String ? v.str() : v.dump(0);
    if (v.kind() == Json::Array) { val.clear(); for (const auto& x : v.items()) val += (val.empty() ? "" : ",") + (x.kind() == Json::String ? x.str() : x.dump(0)); }
    s += (n.size() == 1 ? " -" : " --") + n + " " + q(val);
  }
  return s;
}

Json dft_run(const std::string& c, const Json& a, const std::string& data_dir) {
  cmd_of(c);
  const auto in = inputs(a);
  Json r = Json::object();
  r["ok"] = true;
  r["command"] = dft_command_line(c, a);
  if (c == "sheet") {
    System s;
    std::vector<std::string> notes;
    if (!S(a, "from").empty() || (!in.empty() && S(a, "preset").empty())) {
      const std::string src = S(a, "from", in.empty() ? "" : in[0]);
      IsolateOptions io;
      io.formula = S(a, "formula");
      io.remove = L(a, "remove");
      io.index = int(N(a, "index", 0));
      io.gap = N(a, "gap", 2.0);
      IsolateReport rep;
      s = isolate_layer(read_any(src), io, &rep);
      r["layers_along_z"] = rep.layers, r["complete"] = rep.complete;
      notes = rep.notes;
    } else {
      s = sheet_from_layers(sheet_preset(data_dir, S(a, "preset", "Ti3C2"), N(a, "a", 0), &notes));
    }
    s = group_by_species(s, order_of(a));
    const std::string out = S(a, "o", "sheet.vasp");
    write_any(s, out, formula_ordered(s) + " CAPS sheet");
    provenance(out, "dft.sheet", "2D sheet " + formula_ordered(s), {{"preset", S(a, "preset")}, {"from", S(a, "from")}, {"a", S(a, "a")}});
    r["wrote"] = out, r["formula"] = formula_ordered(s), r["atoms"] = int(s.atoms.size()), r["a"] = norm(s.cell.a);
    Json n = Json::array();
    for (const auto& x : notes) n.push_back(x);
    r["notes"] = n;
    return r;
  }
  if (c == "terminate") {
    if (in.empty()) throw std::invalid_argument("terminate: give the sheet (a file)");
    const System sheet = read_any(in[0]);
    TerminateOptions o;
    o.top.fractions = parse_fractions(S(a, "top", "O"));
    o.top.site = S(a, "site", "fcc");
    o.janus = B(a, "janus");
    o.bottom.fractions = parse_fractions(S(a, "bottom"));
    o.bottom.site = S(a, "site_bottom", o.top.site);
    std::tie(o.na, o.nb) = supercell(a);
    o.vacuum = N(a, "vacuum", 20);
    o.seed = uint64_t(N(a, "seed", 1));
    o.order = order_of(a);
    o.surface_element = S(a, "surface");
    o.data_dir = data_dir;
    for (const auto& kv : L(a, "bond")) { const auto e = kv.find('='); if (e != std::string::npos) o.bond[kv.substr(0, e)] = std::stod(kv.substr(e + 1)); }
    TerminateReport rep;
    const System s = terminate_slab(sheet, o, &rep);
    ValidateOptions vo;
    vo.data_dir = data_dir;
    const auto v = validate_2d(s, vo);
    const std::string out = S(a, "o", "slab.vasp");
    const std::string label = fs::path(out).stem().string();
    write_any(s, out, formula_ordered(s, o.order) + " CAPS " + label);
    provenance(out, "dft.terminate", "terminated " + formula_ordered(s, o.order), {{"top", S(a, "top", "O")}, {"site", o.top.site}, {"supercell", std::to_string(o.na) + "x" + std::to_string(o.nb)}, {"vacuum", S(a, "vacuum", "20")}},
               "numpy default_rng (PCG64) · seed " + std::to_string(o.seed));
    r["wrote"] = out, r["formula"] = formula_ordered(s, o.order), r["atoms"] = int(s.atoms.size()), r["r_site"] = rep.r;
    r["validation"] = v.json();
    r["report"] = v.text(label);
    if (!S(a, "vasp_set").empty()) {
      if (v.status == "FAIL" && !B(a, "force")) throw std::runtime_error("validation failed; nothing written (--force to override)\n" + v.text(label));
      VaspSetOptions so = set_options(a, data_dir);
      so.label = fs::path(S(a, "vasp_set")).filename().string();
      so.stages = stages_slab();
      so.asymmetric = v.info.has("asymmetric") && v.info["asymmetric"].boolean();
      so.report_text = v.text(so.label);
      so.report_json = v.json();
      r["vasp_set"] = set_result(write_vasp_set(s, S(a, "vasp_set"), so));
    }
    return r;
  }
  if (c == "sites") {
    if (in.empty()) throw std::invalid_argument("sites: give the slab");
    const System s = read_any(in[0]);
    Json a2 = Json::array();
    std::map<std::string, int> count;
    // with --species: the height a termination of that kind takes on each site, h = √(d² − r²) (top: d)
    double bond = 0;
    const std::string sp = S(a, "species");
    std::string surf = S(a, "element");
    if (surf.empty()) {
      std::map<int, int> c;
      double zx = -1e300;
      for (const auto& at : s.atoms) if (at.element != 1) zx = std::max(zx, at.pos[2]);
      for (const auto& at : s.atoms) if (at.element != 1 && at.pos[2] > zx - 0.2) ++c[at.element];
      int best = 0, bn = 0;
      for (const auto& [z, n] : c) if (n > bn) bn = n, best = z;
      surf = best ? std::string(element(best).symbol) : "";
    }
    if (!sp.empty()) for (const auto& k : termination_library(data_dir, surf)) if (k.name == sp) bond = k.bond, r["bond"] = k.bond, r["bond_source"] = k.source;
    r["surface_element"] = surf;
    for (const auto& x : surface_sites(s, surf)) {
      Json j = Json::object();
      Json p = Json::array();
      for (double v : x.pos) p.push_back(v);
      j["pos"] = p, j["kind"] = x.kind, j["face"] = x.face, j["beneath"] = x.beneath, j["r"] = x.r;
      if (bond > 0) j["h"] = x.kind == "top" ? bond : (bond > x.r ? std::sqrt(bond * bond - x.r * x.r) : -1.0);
      a2.push_back(j);
      ++count[x.face + " " + x.kind];
    }
    r["sites"] = a2;
    std::string t;
    for (const auto& [k, n] : count) t += k + ": " + std::to_string(n) + "\n";
    r["text"] = t;
    return r;
  }
  if (c == "validate") {
    Json rows = Json::array();
    int bad = 0;
    std::string text;
    for (const auto& p : in) {
      std::map<std::string, std::string> labels;
      std::string comment;
      const System s = read_any(p, &labels, &comment);
      ValidateOptions vo;
      vo.expect = S(a, "expect");
      vo.core = S(a, "core");
      vo.comment = comment;
      vo.bad_labels = labels;
      vo.data_dir = data_dir;
      const auto v = validate_2d(s, vo);
      Json j = v.json();
      j["file"] = p;
      rows.push_back(j);
      text += v.text(fs::path(p).filename().string()) + "\n";
      bad += v.status == "FAIL";
    }
    r["files"] = rows, r["text"] = text, r["failed"] = bad;
    r["ok"] = bad == 0;
    return r;
  }
  if (c == "adsorb-dft") {
    if (in.empty()) throw std::invalid_argument("adsorb-dft: give the slab (a POSCAR / CONTCAR)");
    System slab;
    const auto [na, nb] = supercell(a);
    if (B(a, "from_relaxed")) slab = slab_from_relaxed(read_any(in[0]), na, nb, N(a, "vacuum", 20), order_of(a));
    else {
      slab = read_any(in[0]);
      if (na > 1 || nb > 1) slab = slab_from_relaxed(slab, na, nb, N(a, "vacuum", 20), order_of(a));
    }
    System mol;
    if (!S(a, "molecule").empty()) mol = read_any(S(a, "molecule"));
    else {
      BuildOptions bo;
      bo.forcefield = "uff";
      bo.seed = uint64_t(N(a, "seed", 7));
      mol = build_molecule(S(a, "smiles", "C/C=C/CCC(C#N)C/C=C/C"), bo).system;
    }
    const auto anchors = find_anchors(mol, anchor_library(data_dir));
    AdsorbSetOptions ao;
    ao.modes = L(a, "modes");
    if (a.has("azimuths")) { ao.azimuths.clear(); for (const auto& x : L(a, "azimuths")) ao.azimuths.push_back(std::stod(x)); }
    ao.dmin = N(a, "dmin", 2.3);
    ao.dheavy = N(a, "dheavy", 3.0);
    ao.vacuum = N(a, "vacuum", 20);
    ao.skip_equivalent = B(a, "skip_equivalent");
    ao.order = order_of(a);
    auto set = build_adsorption_set(slab, mol, anchors, ao);
    if (N(a, "prescreen", 0) > 0) {
      // the force-field Adsorption Locator's lowest configurations as further complexes
      PrescreenOptions po;
      po.keep = int(N(a, "prescreen", 0));
      po.cycles = int(N(a, "prescreen_cycles", 3));
      po.steps = int(N(a, "prescreen_steps", 20000));
      po.seed = uint64_t(N(a, "seed", 1));
      std::vector<double> e;
      auto more = prescreen_complexes(set.slab, mol, ao, po, &e);
      Json pe = Json::array();
      for (size_t k = 0; k < more.size(); ++k) { pe.push_back(e[k]); set.complexes.push_back(std::move(more[k])); }
      r["prescreen_energies_kcal"] = pe;
      r["prescreen_note"] = std::to_string(e.size()) + " distinct minimum" + (e.size() == 1 ? "" : "s") + " of the " + std::to_string(po.keep) + " asked (configurations within 0.5 Å RMS of a kept one are the same minimum)";
    }
    Json cx = Json::array();
    for (const auto& c2 : set.complexes) {
      Json j = Json::object();
      j["name"] = c2.name, j["mode"] = c2.mode, j["azimuth"] = c2.azimuth, j["any"] = c2.any, j["heavy"] = c2.heavy, j["image"] = c2.image, j["status"] = c2.status;
      j["equivalent_to"] = c2.equivalent_to, j["suggested_supercell"] = c2.suggested_supercell;
      Json is = Json::array();
      for (const auto& f : c2.issues) { Json x = Json::object(); x["level"] = f.level, x["text"] = f.text; is.push_back(x); }
      j["issues"] = is;
      cx.push_back(j);
    }
    r["complexes"] = cx;
    Json an = Json::object();
    for (const auto& [k, v] : anchors) { Json at = Json::array(); for (int i : v) at.push_back(i); an[k] = at; }
    r["anchors"] = an;
    r["molecule_atoms"] = int(mol.atoms.size());
    Json sym = Json::array();
    for (double d : set.symmetry_rotations) sym.push_back(d);
    r["symmetry_rotations"] = sym;
    Json notes = Json::array();
    for (const auto& n : set.notes) notes.push_back(n);
    r["notes"] = notes;
    const std::string root = S(a, "out");
    if (!root.empty()) {
      VaspSetOptions so = set_options(a, data_dir);
      const auto written = write_adsorption_set(set, root, so, B(a, "force"));
      Json w = Json::array();
      for (const auto& x : written) w.push_back(x.dir);
      r["written"] = w;
      // the independent audit of what was written
      std::vector<std::pair<std::string, System>> cxs;
      for (const auto& c2 : set.complexes) {
        const std::string p = (fs::path(root) / c2.name / "POSCAR").string();
        if (fs::exists(p)) cxs.push_back({c2.name, read_vasp_poscar(p)});
      }
      Json au = Json::array();
      bool all_ok = true;
      for (const auto& row : audit_complexes(read_vasp_poscar((fs::path(root) / "slab" / "POSCAR").string()), read_vasp_poscar((fs::path(root) / "molecule" / "POSCAR").string()), cxs)) {
        Json x = Json::object();
        x["name"] = row.name, x["n_mol"] = row.n_mol, x["any"] = row.any, x["heavy"] = row.heavy, x["image"] = row.image, x["intact"] = row.intact, x["ok"] = row.ok;
        au.push_back(x);
        all_ok = all_ok && row.ok;
      }
      r["audit"] = au;
      r["ok"] = all_ok;
      provenance((fs::path(root) / "slab" / "POSCAR").string(), "dft.adsorb", "adsorption set: " + std::to_string(set.complexes.size()) + " complexes",
                 {{"smiles", S(a, "smiles", "C/C=C/CCC(C#N)C/C=C/C")}, {"dmin", S(a, "dmin", "2.3")}, {"dheavy", S(a, "dheavy", "3.0")}});
    }
    return r;
  }
  if (c == "vasp-set") {
    if (in.empty()) throw std::invalid_argument("vasp-set: give the structure");
    std::map<std::string, std::string> labels;
    std::string comment;
    const System s0 = read_any(in[0], &labels, &comment);
    const System s = group_by_species(s0, order_of(a));
    VaspSetOptions o = set_options(a, data_dir);
    const std::string out = S(a, "out");
    if (out.empty()) throw std::invalid_argument("vasp-set: --out DIR");
    o.label = S(a, "label", fs::path(out).filename().string());
    const std::string st = S(a, "stages", "slab");
    o.stages = st == "slab" ? stages_slab() : st == "fixed" ? stages_fixed() : L(a, "stages");
    ValidateOptions vo;
    vo.data_dir = data_dir;
    vo.bad_labels = labels;
    const auto v = validate_2d(s, vo);
    o.asymmetric = v.info.has("asymmetric") && v.info["asymmetric"].boolean();
    o.report_text = v.text(o.label);
    o.report_json = v.json();
    r["set"] = set_result(write_vasp_set(s, out, o));
    r["validation"] = v.status;
    provenance((fs::path(out) / "POSCAR").string(), "dft.vasp_set", "VASP set " + o.label, {{"stages", st}, {"kdens", S(a, "kdens", "45")}, {"profile", o.profile}});
    return r;
  }
  if (c == "vasp-conv") {
    if (in.size() < 2) throw std::invalid_argument("vasp-conv setup|collect CASE");
    if (in[0] == "setup") {
      std::vector<int> e, k;
      for (const auto& x : L(a, "encuts")) e.push_back(std::stoi(x));
      for (const auto& x : L(a, "kmeshes")) k.push_back(std::stoi(x));
      const auto made = e.empty() && k.empty() ? write_convergence_set(in[1]) : write_convergence_set(in[1], e.empty() ? std::vector<int>{400, 450, 520, 600, 700, 800} : e, k.empty() ? std::vector<int>{6, 9, 12, 15, 18} : k);
      Json m = Json::array();
      for (const auto& d : made) m.push_back(d);
      r["made"] = m;
      return r;
    }
    const auto res = collect_convergence(in[1]);
    auto rows = [](const std::vector<ConvRow>& v) {
      Json a2 = Json::array();
      for (const auto& x : v) { Json j = Json::object(); j["value"] = x.value, j["E_atom"] = x.e_atom, j["stress"] = x.stress, j["fmax"] = x.fmax, j["dE_meV"] = x.de_mev, j["dstress"] = x.ds, j["dforce"] = x.df, j["pass"] = x.pass; a2.push_back(j); }
      return a2;
    };
    r["encut"] = rows(res.encut), r["k"] = rows(res.k), r["encut_choice"] = res.encut_choice, r["k_choice"] = res.k_choice;
    return r;
  }
  if (c == "vasp-scan") {
    if (in.size() >= 5 && in[0] == "make") {
      const auto made = write_lattice_scan(in[1], std::stod(in[2]), std::stod(in[3]), std::stod(in[4]));
      Json m = Json::array();
      for (const auto& d : made) m.push_back(d);
      r["made"] = m;
      return r;
    }
    if (in.size() >= 2 && in[0] == "fit") {
      const auto f = fit_lattice_scan(in[1]);
      Json p = Json::array();
      for (const auto& [x, e] : f.points) { Json q = Json::array(); q.push_back(x); q.push_back(e); p.push_back(q); }
      r["points"] = p, r["a_parabola"] = f.a_parabola, r["inside"] = f.inside;
      if (f.a_bm > 0) r["a_birch_murnaghan"] = f.a_bm, r["B0_GPa"] = f.b0_gpa;
      return r;
    }
    throw std::invalid_argument("vasp-scan make CASE LO HI STEP | fit CASE");
  }
  if (c == "vasp-derived") {
    if (in.size() < 2) throw std::invalid_argument("vasp-derived charge|cdd|freq|aimd CASE");
    const std::string kind = in[0], cd = in[1];
    if (kind == "charge") r["made"] = write_charge_set(cd);
    else if (kind == "cdd") {
      fs::path d = fs::absolute(cd);
      const std::string slab = (d.parent_path() / "slab" / "POSCAR").string();
      const auto made = write_cdd_set(cd, read_vasp_poscar(slab), S(a, "pp_dir"), potcar_map(data_dir));
      Json m = Json::array();
      for (const auto& x : made) m.push_back(x);
      r["made"] = m;
    } else if (kind == "freq") {
      std::vector<int> free;
      for (const auto& x : L(a, "free")) free.push_back(std::stoi(x));
      r["made"] = write_freq_set(cd, free);
      // and the free molecule's, for the shift
      const fs::path mol = fs::absolute(cd).parent_path() / "molecule";
      if (fs::exists(mol / "INCAR.static")) try { r["made_molecule"] = write_freq_set(mol.string(), {}); } catch (const std::exception& e) { r["molecule_note"] = e.what(); }
    } else if (kind == "aimd") {
      if (S(a, "out").empty()) throw std::invalid_argument("vasp-derived aimd STRUCTURE --out DIR");
      VaspSetOptions o = set_options(a, data_dir);
      r["made"] = write_aimd_set(group_by_species(read_any(cd), order_of(a)), S(a, "out"), o, N(a, "temperature", 300), int(N(a, "nsw", 1000)), int(N(a, "segments", 6)));
    } else throw std::invalid_argument("vasp-derived charge|cdd|freq|aimd");
    return r;
  }
  if (c == "vasp-jobs") {
    if (in.empty()) throw std::invalid_argument("vasp-jobs submit|update|reset|cleanup|store …");
    const std::string act = in[0];
    const std::vector<std::string> rest(in.begin() + 1, in.end());
    if (act == "submit") {
      std::vector<std::string> queued;
#ifndef _WIN32
      if (FILE* p = popen("squeue -u \"$USER\" -h -o %Z 2>/dev/null", "r")) {
        char b[4096];
        while (std::fgets(b, sizeof b, p)) { std::string l(b); while (!l.empty() && (l.back() == '\n' || l.back() == '\r')) l.pop_back(); if (!l.empty()) queued.push_back(l); }
        pclose(p);
      }
#endif
      const auto go = submittable(rest, queued);
      Json s = Json::array();
      for (const auto& d : go) {
        Json j = Json::object();
        j["dir"] = d;
        if (B(a, "run")) {
          const std::string cmd = "cd '" + d + "' && { [ -f POTCAR ] || ./make_potcar.sh; } && sbatch job.slurm";
          j["exit"] = std::system(cmd.c_str());
        }
        s.push_back(j);
      }
      r["submit"] = s;
      r["skipped"] = int(rest.size() - go.size());
      return r;
    }
    if (act == "update") { r["updated"] = update_job_scripts(rest, cluster_profile(data_dir, S(a, "profile", "slurm-workspace"))); return r; }
    if (act == "reset") { if (rest.size() < 2) throw std::invalid_argument("vasp-jobs reset CASE STAGE"); r["moved_to"] = reset_stage(rest[0], rest[1]); return r; }
    if (act == "cleanup") {
      const auto plan = cleanup_plan(rest.empty() ? "." : rest[0], B(a, "keep_best"));
      r["plan"] = plan_json(plan);
      if (B(a, "yes")) apply_plan(plan, "cleanup"), r["done"] = true;
      return r;
    }
    if (act == "store") {
      if (rest.size() < 2) throw std::invalid_argument("vasp-jobs store ROOT STORAGE_PATH");
      const auto plan = store_plan(rest[0], fs::absolute(rest[1]).string());
      r["plan"] = plan_json(plan);
      if (B(a, "yes")) apply_plan(plan, "store", rest[0], fs::absolute(rest[1]).string()), r["done"] = true;
      return r;
    }
    throw std::invalid_argument("vasp-jobs submit|update|reset|cleanup|store");
  }
  if (c == "vasp-check") { r["runs"] = check_runs(in, B(a, "complex"), data_dir); bool ok = true; for (const auto& x : r["runs"].items()) ok = ok && x.has("ok") && x["ok"].boolean(); r["ok"] = ok; return r; }
  if (c == "vasp-progress") { Json p = Json::array(); for (const auto& d : in) p.push_back(run_progress(d)); r["cases"] = p; return r; }
  if (c == "vasp-health") { if (in.empty()) throw std::invalid_argument("vasp-health SET"); r["health"] = health_check(in[0]); r["ok"] = !r["health"].has("error") && r["health"].num("flagged", 0) == 0; return r; }
  if (c == "vasp-bind") { if (in.empty()) throw std::invalid_argument("vasp-bind SET"); r["binding"] = binding_energies(in[0]); r["ok"] = !r["binding"].has("error"); return r; }
  if (c == "vasp-analyze") {
    if (in.empty()) throw std::invalid_argument("vasp-analyze geom|wf|dos|cdd|bader|freq|md|summary …");
    const std::string k = in[0];
    const std::vector<std::string> rest(in.begin() + 1, in.end());
    auto need = [&](size_t n) { if (rest.size() < n) throw std::invalid_argument("vasp-analyze " + k + ": missing arguments (caps vasp-analyze --help)"); };
    if (k == "geom") r["result"] = adsorption_geometry(rest);
    else if (k == "wf") { need(1); r["result"] = work_function(rest[0], S(a, "reference")); }
    else if (k == "dos") { need(1); r["result"] = dos_analysis(rest[0], S(a, "slab_poscar")); }
    else if (k == "cdd") { need(1); r["result"] = charge_density_difference(rest[0], (fs::absolute(rest[0]).parent_path() / "slab" / "POSCAR").string()); }
    else if (k == "bader") { need(3); r["result"] = bader_transfer(rest[0], rest[1], rest[2], rest.size() > 3 ? rest[3] : ""); }
    else if (k == "freq") { need(1); r["result"] = nitrile_shift(rest[0]); }
    else if (k == "md") { need(1); r["result"] = md_analysis(rest[0], N(a, "equil_ps", 1.0), N(a, "dt_fs", 1.0)); }
    else if (k == "summary") r["result"] = summarize_results(rest.empty() ? "." : rest[0]);
    else throw std::invalid_argument("vasp-analyze geom|wf|dos|cdd|bader|freq|md|summary");
    if (r["result"].has("error")) r["ok"] = false;
    return r;
  }
  throw std::invalid_argument("unknown command " + c);
}

}  // namespace caps
