// The coarse-graining workflow as commands (caps/cg_commands.hpp).
#include "caps/cg_commands.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <limits>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>

#include "caps/bundle.hpp"
#include "caps/cg_bonded.hpp"
#include "caps/cg_nonbonded.hpp"
#include "caps/mechanics.hpp"
#include "caps/trajectory.hpp"
#include "caps/cg_rules.hpp"
#include "caps/dft_commands.hpp"
#include "caps/io.hpp"
#include "caps/provenance.hpp"

namespace caps {

namespace fs = std::filesystem;

namespace {

struct Opt { const char* name; const char* def; const char* why; };
struct Cmd { const char* name; const char* usage; const char* about; std::vector<Opt> opts; std::vector<const char*> examples; };

const std::vector<Cmd>& table() {
  static const std::vector<Cmd> t = {
      {"cgmap", "caps cgmap STRUCTURE [STRUCTURE …] (--preset ID | --rules FILE | --cut SMARTS --names N=SMARTS,… | --map FILE) [-o DIR] [--dump DUMP]",
       "Chemistry-aware coarse-grained mapping: bonds cut by SMARTS (the fragments are the beads), fragment SMARTS, or an atom → bead "
       "list; beads named by rules (unknown fragments reported, never merged), at the centre of mass, centre of geometry or one atom. "
       "Several structures share one type list (bead, bond, angle, dihedral types numbered alike, so tables are shared). Writes per "
       "structure STEM.cg.data and STEM.map.json; with --dump maps a LAMMPS dump frame by frame (one frame in memory).",
       {{"preset", "", "a rule set from data/cg/mapping_rules.json (ester-cut: polyesters, B S A T beads)"},
        {"rules", "", "a rule file: {\"cut\": [SMARTS …] | \"beads\": {NAME: SMARTS} | \"explicit\": [bead per atom], \"names\": {NAME: SMARTS}, \"position\", \"strict\"}"},
        {"cut", "", "bond SMARTS, comma-separated: the bond between each pattern's first and last atoms is cut"},
        {"names", "", "naming rules NAME=SMARTS, comma-separated, first match wins"},
        {"map", "", "an existing map.json (the same atoms): used as it is"},
        {"position", "com", "com (centre of mass), cog (centre of geometry) or atom:SMARTS (the first atom of its first match in the bead)"},
        {"no_strict", "false", "name fragments no rule matches by their class (F1, F2 …) instead of stopping"},
        {"types", "", "a shared type list (JSON: beads, bonds, angles, dihedrals); default: the union over the structures given, written to DIR/types.json"},
        {"o", ".", "the output folder"},
        {"dump", "", "a LAMMPS dump of the (single) structure: id plus xu yu zu, x y z (ix iy iz), or xs ys zs"},
        {"out_dump", "", "the mapped dump (default DIR/STEM.cg.lammpstrj)"},
        {"stride", "1", "every n-th frame of the dump"},
        {"first", "1", "the first frame kept (1-based)"},
        {"last", "", "the last frame kept (1-based; default: the end)"}},
       {"caps cgmap PBS/system.data PBSA/system.data PBAT/system.data --preset ester-cut -o cg",
        "caps cgmap system.data --map cg/system.map.json --types cg/types.json --dump hold.lammpstrj --stride 5 -o cg"}},
      {"cgfit", "caps cgfit bonded|refine|targets|ibi-start|ibi-step|ibi-run|fit|tg|calibrate …  (see below)",
       "Bonded coarse-grained potentials by tabulated Boltzmann inversion over many frames of many systems mapped alike: each MAP "
       "(STEM.map.json from cgmap) followed by its frames (CG dumps STEM.cg.lammpstrj or the single frame STEM.cg.data); each "
       "system's histograms normalised, then weighted. U(r) = −kT ln[P/r²], U(θ) = −kT ln[P/sin θ], U(φ) = −kT ln P, inverted only "
       "where P exceeds a share of its maximum, with slope-continuous walls outside. refine: one bonded IBI step from a CG run "
       "with the tables (U ← U + α kT ln(P_CG/P_target)), for the shift the non-bonded 1–3 and 1–4 terms cause.\n"
       "Non-bonded: targets (per-pair g(r) per system from MAP FRAMES …); ibi-start (−kT ln g tables, LAMMPS decks and run_ibi.sh, "
       "the loop for a cluster); ibi-step (one joint update from each system's CG dump and log: MAP DUMP [LOG] …); ibi-run (the loop "
       "with CAPS's engine, small systems: MAP DATA …); fit (LJ 12-6, LJ 9-6, Morse or Mie fitted to the tables); tg (T_g from a "
       "cooling scan's T and density columns); calibrate (the next σ and ε scale factors from the runs so far, with decks for them).",
       {{"types", "", "the shared type list (types.json from cgmap); default: the union of the maps given"},
        {"T", "300", "temperature of the inversion (K): the AA run's"},
        {"weights", "", "one weight per system, comma-separated (default: equal; each system is normalised, so frame counts do not weigh)"},
        {"stride", "1", "every n-th frame of each dump"},
        {"threshold", "0.05", "invert where the smoothed distribution exceeds this share of its maximum; walls outside"},
        {"smooth", "1", "Gaussian smoothing of the histograms, σ in bins (1: the round trip on an ideal chain stays within 0.06 kT)"},
        {"bin_bond", "0.02", "bond histogram bin (Å)"},
        {"bin_angle", "1", "angle histogram bin (degrees)"},
        {"bin_dihedral", "5", "dihedral histogram bin (degrees)"},
        {"bond_wall", "20", "curvature of the bond walls, kcal/mol/Å² (at least the well's)"},
        {"angle_wall", "0.005", "curvature of the angle walls, kcal/mol/deg²"},
        {"aat", "1,170,180", "dihedral_style table/cut switch K, θ1, θ2: the dihedral turns off as a neighbouring angle goes from θ1 to θ2"},
        {"tables", "", "refine: the bonded.json of the tables the CG run used"},
        {"alpha", "0.5", "refine: damping of the update"},
        {"rmax", "15", "targets: g(r) range (Å; at most half the box)"},
        {"dr", "0.05", "targets: g(r) bin (Å)"},
        {"exclude", "3", "targets: bead pairs fewer than this many bonds apart left out (3: 1-2 and 1-3, special_bonds lj 0 0 1)"},
        {"pressure", "1", "targets: each system's AA pressure (atm, comma-separated): the IBI pressure target"},
        {"targets", "", "ibi-*: the targets.json"},
        {"pairs", "", "ibi-step / fit: the pairs.json of the current tables"},
        {"bonded", "", "ibi-start / ibi-run: the bonded folder (bonded.in, bonded.json)"},
        {"rc", "15", "ibi-start: the tables' cut-off (Å)"},
        {"ramp", "0.1", "ibi-*: the pressure correction's size, A = −ramp kT sign(ΔP) min(1, 0.0003|ΔP|/bar) (0: off)"},
        {"iterations", "3", "ibi-run: iterations with CAPS's engine"},
        {"steps", "2000", "ibi-run: steps per iteration (CAPS's engine; keep runs short on a laptop) / the decks' production steps"},
        {"dt", "10", "time step of the runs (fs)"},
        {"form", "lj126", "fit: lj126, lj96, morse or mie"},
        {"repulsion", "4", "fit: the fitted range starts where U falls below this many kT"},
        {"fits", "", "calibrate: the fits.json to scale"},
        {"history", "", "calibrate: calibration.json (the runs so far; created when missing)"},
        {"s_sigma", "1", "calibrate: the σ scale factor of the run being reported"},
        {"s_eps", "1", "calibrate: the ε scale factor of the run being reported"},
        {"density", "", "calibrate: that run's density (g/cm³)"},
        {"tg", "", "calibrate: that run's T_g (K), e.g. from caps cgfit tg"},
        {"target_density", "", "calibrate: the AA (or experimental) density, g/cm³"},
        {"target_tg", "", "calibrate: the AA (or experimental) T_g, K"},
        {"t_hi", "500", "decks: the cooling scan's start (K)"},
        {"t_lo", "150", "decks: the cooling scan's end (K)"},
        {"t_step", "25", "decks: the cooling scan's step (K)"},
        {"o", "bonded", "the output folder"}},
       {"caps cgfit bonded cg/PBS.map.json cg/PBS.cg.lammpstrj cg/PBSA.map.json cg/PBSA.cg.lammpstrj cg/PBAT.map.json cg/PBAT.cg.lammpstrj --types cg/types.json -T 300 -o bonded",
        "caps cgfit refine cg/PBS.map.json run1/cg.lammpstrj --tables bonded/bonded.json --types cg/types.json -o bonded_it2",
        "caps cgfit targets cg/PBS.map.json cg/PBS.cg.lammpstrj cg/PBAT.map.json cg/PBAT.cg.lammpstrj --types cg/types.json -T 300 -o ibi",
        "caps cgfit ibi-start --targets ibi/targets.json --bonded bonded --types cg/types.json -o ibi   (then: bash ibi/run_ibi.sh 20)",
        "caps cgfit ibi-step cg/PBS.map.json ibi/it000/PBS.lammpstrj ibi/it000/PBS.log … --targets ibi/targets.json --pairs ibi/it000/pairs.json -o ibi/it001",
        "caps cgfit fit --pairs ibi/it020/pairs.json --form lj126 --types cg/types.json -o lj",
        "caps cgfit calibrate --fits lj/fits.json --history calib/calibration.json --s-sigma 1 --s-eps 1 --density 1.31 --tg 365 --target-density 1.25 --target-tg 241 -o calib/step1"}},
  };
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
// a comma-separated list (SMARTS hold commas inside brackets: split only outside [ ])
std::vector<std::string> L(const Json& a, const std::string& k) {
  std::vector<std::string> out;
  if (!a.has(k)) return out;
  if (a[k].kind() == Json::Array) { for (const auto& x : a[k].items()) out.push_back(x.str()); return out; }
  const std::string s = S(a, k);
  std::string cur;
  int depth = 0;
  for (char ch : s) {
    if (ch == '[') ++depth;
    if (ch == ']') --depth;
    if (ch == ',' && depth == 0) { if (!cur.empty()) out.push_back(cur); cur.clear(); continue; }
    cur += ch;
  }
  if (!cur.empty()) out.push_back(cur);
  for (auto& x : out) { x.erase(0, x.find_first_not_of(' ')); x.erase(x.find_last_not_of(' ') + 1); }
  return out;
}
std::vector<std::string> inputs(const Json& a) {
  std::vector<std::string> out;
  if (a.has("inputs")) for (const auto& x : a["inputs"].items()) out.push_back(x.str());
  return out;
}
Json read_json_file(const std::string& p) {
  std::ifstream f(p);
  if (!f) throw std::invalid_argument("cannot open " + p);
  std::stringstream ss;
  ss << f.rdbuf();
  return Json::parse(ss.str());
}
void write_text(const std::string& p, const std::string& text) {
  std::ofstream f(p);
  if (!f) throw std::runtime_error("cannot write " + p);
  f << text;
}
std::string stem_of(const std::string& p) {
  fs::path x(p);
  std::string s = x.stem().string();
  // system.data in PBS/DP-25-40: name it by its folders so several systems do not overwrite each other
  if (s == "system" || s == "equil" || s == "data") {
    const auto parent = x.parent_path();
    if (!parent.empty()) {
      const std::string a = parent.filename().string(), b = parent.parent_path().filename().string();
      s = (b.empty() ? "" : b + "_") + a;
    }
  }
  // file names that are easy to type: letters, digits, . _ - (folder names such as "PBSA 11.27.12 AM" hold odd spaces)
  std::string out;
  for (unsigned char ch : s) {
    const bool ok = std::isalnum(ch) || ch == '.' || ch == '_' || ch == '-';
    if (ok) out += char(ch);
    else if (ch < 0x80 || (ch & 0xC0) == 0xC0) { if (out.empty() || out.back() != '_') out += '_'; }   // one _ per space or non-ASCII character
  }
  while (!out.empty() && out.back() == '_') out.pop_back();
  return out.empty() ? "system" : out;
}
std::string sha_of(const std::string& p) {
  std::ifstream f(p, std::ios::binary);
  std::stringstream ss;
  ss << f.rdbuf();
  return sha256_hex(ss.str());
}

void provenance(const std::string& path, const std::string& engine, const std::string& summary, const KeyValues& params,
                const KeyValues& ins, const std::vector<std::string>& cites, const KeyValues& approx = {}) {
  Manifest m;
  m.inputs = ins;
  ProvStep st;
  st.engine = engine;
  st.summary = summary;
  st.params = params;
  st.cites = cites;
  st.approximations = approx;
  st.time = now_iso();
  m.steps.push_back(st);
  try { write_manifest(m, path); } catch (...) {}
}

CgRules rules_of(const Json& a, const std::string& data_dir, std::string* label) {
  CgRules r;
  if (!S(a, "preset").empty()) {
    const Json lib = read_json_file((fs::path(data_dir) / "cg" / "mapping_rules.json").string());
    bool found = false;
    for (const auto& p : lib["presets"].items())
      if (p.text("id") == S(a, "preset")) { r = cg_rules_from_json(p); found = true; }
    if (!found) {
      std::string ids;
      for (const auto& p : lib["presets"].items()) ids += (ids.empty() ? "" : ", ") + p.text("id");
      throw std::invalid_argument("no mapping preset " + S(a, "preset") + " (there are: " + ids + ")");
    }
    *label = "preset " + S(a, "preset");
  } else if (!S(a, "rules").empty()) {
    r = cg_rules_from_json(read_json_file(S(a, "rules")));
    *label = "rules " + fs::path(S(a, "rules")).filename().string();
  }
  if (a.has("cut")) r.cut = L(a, "cut"), r.beads.clear(), r.explicit_bead.clear();
  if (a.has("names")) {
    r.names.clear();
    for (const auto& x : L(a, "names")) {
      const auto e = x.find('=');
      if (e == std::string::npos) throw std::invalid_argument("--names takes NAME=SMARTS, not " + x);
      r.names.push_back({x.substr(0, e), x.substr(e + 1)});
    }
  }
  if (label->empty()) *label = "rules given on the command line";
  if (a.has("position")) r.position = S(a, "position");
  if (B(a, "no_strict")) r.strict = false;
  return r;
}


// ---------------------------------------------------------------- coarse-graining helpers
struct CgSysIn { std::string map; std::vector<std::string> frames, logs; Json mj; CgTopology t; };

std::vector<CgSysIn> cg_systems(const std::vector<std::string>& in, size_t from) {
  std::vector<CgSysIn> sys;
  for (size_t k = from; k < in.size(); ++k) {
    const std::string& p = in[k];
    const std::string name = fs::path(p).filename().string();
    if (p.size() >= 9 && p.compare(p.size() - 9, 9, ".map.json") == 0) {
      CgSysIn x;
      x.map = p;
      x.mj = read_json_file(p);
      x.t = cg_topology_from_map(x.mj);
      sys.push_back(std::move(x));
    } else {
      if (sys.empty()) throw std::invalid_argument("give a mapping (STEM.map.json) before its frames: " + p);
      if (fs::path(p).extension() == ".log" || name.rfind("log.lammps", 0) == 0) sys.back().logs.push_back(p);
      else sys.back().frames.push_back(p);
    }
  }
  return sys;
}

CgTypes cg_types_for(const Json& a, const std::vector<CgSysIn>& sys) {
  if (!S(a, "types").empty()) return cg_types_from_json(read_json_file(S(a, "types")));
  std::vector<CgTypes> lists;
  for (const auto& x : sys) {
    std::set<std::string> b, bo, an, di;
    for (const auto& k : x.t.kind) b.insert(k);
    for (const auto& [i, j] : x.t.bonds) bo.insert(cg_key({x.t.kind[size_t(i)], x.t.kind[size_t(j)]}));
    for (const auto& q : x.t.angles) an.insert(cg_key({x.t.kind[size_t(q[0])], x.t.kind[size_t(q[1])], x.t.kind[size_t(q[2])]}));
    for (const auto& q : x.t.dihedrals) di.insert(cg_key({x.t.kind[size_t(q[0])], x.t.kind[size_t(q[1])], x.t.kind[size_t(q[2])], x.t.kind[size_t(q[3])]}));
    lists.push_back({{b.begin(), b.end()}, {bo.begin(), bo.end()}, {an.begin(), an.end()}, {di.begin(), di.end()}});
  }
  return cg_types_union(lists);
}

// each kind's most common bead mass (chain ends carry their end groups' atoms too)
std::map<std::string, double> cg_masses(const std::vector<CgSysIn>& sys) {
  std::map<std::string, std::map<long, int>> h;
  for (const auto& x : sys)
    for (const auto& b : x.mj["beads"].items()) h[b["kind"].str()][std::lround(b.num("mass", 0) * 1000)]++;
  std::map<std::string, double> out;
  for (const auto& [k, m] : h) out[k] = double(std::max_element(m.begin(), m.end(), [](auto& x, auto& y) { return x.second < y.second; })->first) / 1000.0;
  return out;
}

// every frame of a system: CG dumps streamed, or STEM.cg.data
template <class F>
size_t cg_frames(const CgSysIn& x, size_t stride, F&& f) {
  size_t n = 0;
  for (const auto& p : x.frames) {
    if (fs::path(p).extension() == ".data") {
      const System d = read_lammps_data(p);
      if (d.atoms.size() != x.t.beads()) throw std::invalid_argument(p + " has " + std::to_string(d.atoms.size()) + " beads, its map " + std::to_string(x.t.beads()));
      std::vector<Vec3> q(d.atoms.size());
      for (size_t b = 0; b < d.atoms.size(); ++b) q[size_t(d.atoms[b].id > 0 ? d.atoms[b].id - 1 : int64_t(b))] = d.atoms[b].pos;
      f(q, d.cell);
      ++n;
    } else {
      n += for_each_dump_frame(p, x.t.beads(), [&](size_t, int64_t, const std::vector<Vec3>& q, const Cell& c, bool) { f(q, c); return true; }, stride);
    }
  }
  if (!n) throw std::invalid_argument("no frames for " + x.map + ": give its CG dump(s) or STEM.cg.data after it");
  return n;
}

// the mean of a LAMMPS log's thermo column over the second half of its last run
double log_mean(const std::string& path, const char* column) {
  const ThermoLog L = read_lammps_log(path);
  size_t c = SIZE_MAX;
  for (size_t k = 0; k < L.columns.size(); ++k) if (L.columns[k] == column) c = k;
  if (c == SIZE_MAX || L.rows.empty()) throw std::invalid_argument(path + ": no " + std::string(column) + " column in its thermo output");
  const size_t start = L.run_starts.empty() ? 0 : L.run_starts.back();
  const size_t n = L.rows.size() - start, from = start + n / 2;
  double sum = 0;
  for (size_t k = from; k < L.rows.size(); ++k) sum += L.rows[k][c];
  return sum / double(L.rows.size() - from);
}

std::string stem_of_map(const std::string& map) {
  std::string n = fs::path(map).filename().string();
  return n.size() > 9 ? n.substr(0, n.size() - 9) : n;
}

// LAMMPS inputs for the CG model: in.cg_run (NVT or NPT, a dump for g(r) and the next IBI step) and in.cg_tg (a stepwise
// cooling scan writing T and density to tg.dat for caps cgfit tg)
void write_cg_decks(const std::string& dir, int exclude, double T, double dt, int64_t steps, double t_hi, double t_lo, double t_step, std::vector<std::string>& files) {
  const char* special = exclude <= 2 ? "0 1 1" : exclude == 3 ? "0 0 1" : "0 0 0";
  char b[512];
  {
    const std::string p = (fs::path(dir) / "in.cg_run").string();
    std::ofstream f(p);
    f << "# CAPS coarse-grained run (LAMMPS; dihedral_style table/cut needs the EXTRA-MOLECULE package)\n"
         "#   lmp -in in.cg_run -var DATA system.cg.data -var BONDED ../bonded -var PAIR it000 -var OUT it000/system [-var ENS npt]\n";
    std::snprintf(b, sizeof b, "variable T index %g\nvariable P index 1.0\nvariable DT index %g\nvariable STEPS index %lld\nvariable ENS index nvt\nvariable DUMP index 1000\nvariable SEED index 4928459\n",
                  T, dt, static_cast<long long>(steps));
    f << b;
    f << "units real\natom_style full\nboundary p p p\n";
    f << "special_bonds lj " << special << "\n";
    f << "read_data ${DATA}\ninclude ${BONDED}/bonded.in\ninclude ${PAIR}/pair.in\n"
         "neighbor 3.0 bin\nneigh_modify delay 0 every 1 check yes\ntimestep ${DT}\n"
         "velocity all create ${T} ${SEED} mom yes rot yes dist gaussian\n"
         "if \"${ENS} == npt\" then \"fix md all npt temp ${T} ${T} $(100.0*dt) iso ${P} ${P} $(1000.0*dt)\" else \"fix md all nvt temp ${T} ${T} $(100.0*dt)\"\n"
         "thermo_style custom step temp press density pe ke etotal\nthermo 1000\n"
         "dump cg all custom ${DUMP} ${OUT}.lammpstrj id mol type xu yu zu\ndump_modify cg sort id\n"
         "log ${OUT}.log\nrun ${STEPS}\nwrite_data ${OUT}.data\n";
    files.push_back(p);
  }
  {
    const std::string p = (fs::path(dir) / "in.cg_tg").string();
    std::ofstream f(p);
    f << "# CAPS coarse-grained T_g scan (LAMMPS): NPT at each temperature, the mean density written to ${OUT}.tg.dat\n"
         "#   lmp -in in.cg_tg -var DATA system.cg.data -var BONDED ../bonded -var PAIR lj -var OUT tg/system;  caps cgfit tg tg/system.tg.dat\n";
    std::snprintf(b, sizeof b, "variable DT index %g\nvariable EQ index 100000\nvariable SAMPLE index 100000\nvariable P index 1.0\nvariable SEED index 4928459\n", dt);
    f << b;
    f << "units real\natom_style full\nboundary p p p\nspecial_bonds lj " << special << "\n"
         "read_data ${DATA}\ninclude ${BONDED}/bonded.in\ninclude ${PAIR}/pair.in\nneighbor 3.0 bin\nneigh_modify delay 0 every 1 check yes\ntimestep ${DT}\n";
    std::snprintf(b, sizeof b, "velocity all create %g ${SEED} mom yes rot yes dist gaussian\nthermo_style custom step temp press density pe\nthermo 5000\n"
                  "variable rho equal density\nprint \"# T (K)  density (g/cm3)\" file ${OUT}.tg.dat\n", t_hi);
    f << b;
    for (double t = t_hi; t >= t_lo - 1e-9; t -= t_step) {
      std::snprintf(b, sizeof b, "fix md all npt temp %g %g $(100.0*dt) iso ${P} ${P} $(1000.0*dt)\nrun ${EQ}\nfix avg all ave/time 10 $(v_SAMPLE/10) ${SAMPLE} v_rho\nrun ${SAMPLE}\n"
                    "print \"%g $(f_avg)\" append ${OUT}.tg.dat\nunfix avg\nunfix md\n", t, t, t);
      f << b;
    }
    files.push_back(p);
  }
}

// the non-bonded modes of cgfit
void cgfit_nonbonded(const std::string& mode, const Json& a, const std::vector<std::string>& in, Json& r, Json& files) {
  std::ostringstream t;
  char b[512];
  const double T = N(a, "T", 300);
  if (mode == "targets") {
    const auto sys = cg_systems(in, 1);
    if (sys.empty()) throw std::invalid_argument("targets: give each system's map.json and its frames");
    const CgTypes types = cg_types_for(a, sys);
    CgPairOptions po;
    po.rmax = N(a, "rmax", 15), po.dr = N(a, "dr", 0.05), po.exclude = int(N(a, "exclude", 3));
    CgRdfAccumulator acc(types, po);
    std::vector<double> press;
    const auto pl = L(a, "pressure");
    KeyValues ins;
    for (size_t k = 0; k < sys.size(); ++k) {
      const int si = acc.add_system(sys[k].t, 1, stem_of_map(sys[k].map));
      const size_t n = cg_frames(sys[k], size_t(std::max(1.0, N(a, "stride", 1))), [&](const std::vector<Vec3>& q, const Cell& c) { acc.add_frame(si, q, c); });
      press.push_back(k < pl.size() ? std::stod(pl[k]) : pl.size() == 1 ? std::stod(pl[0]) : 1.0);
      std::snprintf(b, sizeof b, "%s: %zu frames\n", stem_of_map(sys[k].map).c_str(), n);
      t << b;
      ins.push_back({fs::path(sys[k].map).filename().string(), sha_of(sys[k].map)});
    }
    const std::string dir = S(a, "o", "ibi");
    fs::create_directories(dir);
    CgTargets tg = cg_targets(acc, T, press);
    const std::string tp = (fs::path(dir) / "targets.json").string();
    write_text(tp, cg_targets_json(tg).dump(0) + "\n");
    write_text((fs::path(dir) / "types.json").string(), cg_types_json(types).dump(1) + "\n");
    files.push_back(tp);
    files.push_back((fs::path(dir) / "types.json").string());
    // a CSV per system for plots: r and every pair's g
    for (size_t s2 = 0; s2 < tg.systems.size(); ++s2) {
      const std::string cp = (fs::path(dir) / (tg.systems[s2] + ".gr.csv")).string();
      std::ofstream f(cp);
      f << "r";
      for (const auto& R : tg.rdf[s2]) f << "," << R.key;
      f << "\n";
      for (size_t k = 0; k < tg.r.size(); ++k) {
        f << tg.r[k];
        for (const auto& R : tg.rdf[s2]) f << "," << R.g[k];
        f << "\n";
      }
      files.push_back(cp);
    }
    // per pair: where g first reaches 1 (contact) and its first peak, pooled
    for (const auto& R0 : tg.rdf.front()) {
      double contact = 0, peak = 0, gpk = 0;
      for (size_t s2 = 0; s2 < tg.rdf.size(); ++s2)
        for (const auto& R : tg.rdf[s2])
          if (R.key == R0.key && R.ideal_pairs > 0)
            for (size_t k = 0; k < tg.r.size(); ++k) {
              if (contact == 0 && R.g[k] >= 1) contact = tg.r[k];
              if (R.g[k] > gpk) gpk = R.g[k], peak = tg.r[k];
            }
      if (gpk > 0) { std::snprintf(b, sizeof b, "  %-5s g = 1 at %.2f Å · first peak %.2f at %.2f Å\n", R0.key.c_str(), contact, gpk, peak); t << b; }
      else t << "  " << R0.key << ": no pairs in any system\n";
    }
    t << "wrote " << tp << "\n";
    provenance(tp, "cg.targets", "per-pair g(r) of " + std::to_string(sys.size()) + " systems", {{"rmax", S(a, "rmax", "15")}, {"dr", S(a, "dr", "0.05")}, {"exclude", S(a, "exclude", "3")}}, ins, {"reith2003"});
  } else if (mode == "ibi-start") {
    if (S(a, "targets").empty()) throw std::invalid_argument("ibi-start needs --targets (targets.json)");
    const CgTargets tg = cg_targets_from_json(read_json_file(S(a, "targets")));
    const std::string dir = S(a, "o", fs::path(S(a, "targets")).parent_path().string());
    CgIbiOptions io;
    io.rc = N(a, "rc", 15);
    io.table_dr = N(a, "table_dr", 0.05);
    const CgPairSet p = ibi_start(tg, io);
    CgTypes types;
    if (!S(a, "types").empty()) types = cg_types_from_json(read_json_file(S(a, "types")));
    else types = cg_types_from_json(read_json_file((fs::path(S(a, "targets")).parent_path() / "types.json").string()));
    const std::string it0 = (fs::path(dir) / "it000").string();
    for (const auto& f : write_pairs(p, types, it0)) files.push_back(f);
    std::vector<std::string> deckfiles;
    write_cg_decks(dir, int(N(a, "exclude", 3)), tg.temperature, N(a, "dt", 10), int64_t(N(a, "steps", 500000)), N(a, "t_hi", 500), N(a, "t_lo", 150), N(a, "t_step", 25), deckfiles);
    for (const auto& f : deckfiles) files.push_back(f);
    // the loop for a cluster (bash; LAMMPS as $CAPS_LMP, CAPS as $CAPS_BIN)
    const std::string bonded = S(a, "bonded", "../bonded");
    const std::string sh = (fs::path(dir) / "run_ibi.sh").string();
    {
      std::ofstream f(sh);
      f << "#!/bin/bash\n# CAPS joint IBI loop: each iteration runs every system with the current tables (LAMMPS), then caps cgfit ibi-step.\n"
           "#   bash run_ibi.sh N [START]   (from this folder; LAMMPS: $CAPS_LMP, default lmp; CAPS: $CAPS_BIN, default caps)\n"
           "# Each system needs its starting CG structure as DATA_<name> below (STEM.cg.data from caps cgmap, or an equilibrated one).\n"
           "set -euo pipefail\nLMP=${CAPS_LMP:-lmp}\nCAPS=${CAPS_BIN:-caps}\nN=${1:-20}\nSTART=${2:-0}\n";
      f << "BONDED=" << bonded << "\nTARGETS=targets.json\nTYPES=types.json\n";
      // shell variables by number (system names hold '-' and '.'), each system named in a comment; override with SYS1_DATA=… bash run_ibi.sh
      for (size_t k = 0; k < tg.systems.size(); ++k) {
        const std::string& nme = tg.systems[k];
        const std::string v = "SYS" + std::to_string(k + 1);
        f << "# system " << k + 1 << ": " << nme << "\n" << v << "_DATA=${" << v << "_DATA:-../cg/" << nme << ".cg.data}\n" << v << "_MAP=${" << v << "_MAP:-../cg/" << nme << ".map.json}\n";
      }
      f << "for ((i=START; i<N; i++)); do\n  it=$(printf 'it%03d' $i); nx=$(printf 'it%03d' $((i+1))); prev=$(printf 'it%03d' $((i-1)))\n  args=()\n";
      for (size_t k = 0; k < tg.systems.size(); ++k) {
        const std::string& nme = tg.systems[k];
        const std::string v = "SYS" + std::to_string(k + 1);
        f << "  data=$" << v << "_DATA; if [ $i -gt 0 ] && [ -f \"$prev/" << nme << ".data\" ]; then data=\"$prev/" << nme << ".data\"; fi\n"
             "  \"$LMP\" -in in.cg_run -var DATA \"$data\" -var BONDED \"$BONDED\" -var PAIR $it -var OUT \"$it/" << nme << "\" -var ENS nvt\n"
             "  args+=(\"$" << v << "_MAP\" \"$it/" << nme << ".lammpstrj\" \"$it/" << nme << ".log\")\n";
      }
      f << "  \"$CAPS\" cgfit ibi-step \"${args[@]}\" --targets $TARGETS --pairs $it/pairs.json --types $TYPES -o $nx | tee $nx.report.txt\n"
           "done\n";
      files.push_back(sh);
    }
    std::snprintf(b, sizeof b, "IBI start: %zu pair tables (−kT ln g, cut-off %.1f Å) in %s\n", p.pairs.size(), p.rc, it0.c_str());
    t << b;
    for (const auto& R0 : tg.rdf.front())
      if (std::none_of(p.pairs.begin(), p.pairs.end(), [&](const CgPairTable& q) { return q.key == R0.key; })) t << "  no table for " << R0.key << " (no system has both types)\n";
    t << "LAMMPS: in.cg_run (NVT/NPT run with a dump for the next step), in.cg_tg (cooling scan); loop: bash " << sh << " 20\n";
  } else if (mode == "ibi-step" || mode == "ibi-run") {
    if (S(a, "targets").empty()) throw std::invalid_argument(mode + " needs --targets (targets.json)");
    const CgTargets tg = cg_targets_from_json(read_json_file(S(a, "targets")));
    const auto sys = cg_systems(in, 1);
    if (sys.size() != tg.systems.size()) throw std::invalid_argument(mode + ": give the " + std::to_string(tg.systems.size()) + " systems of the targets, in their order (" + std::to_string(sys.size()) + " given)");
    CgTypes types = cg_types_for(a, sys);
    CgPairOptions po;
    po.rmax = tg.r.back() + 0.5 * (tg.r[1] - tg.r[0]), po.dr = tg.r[1] - tg.r[0], po.exclude = int(N(a, "exclude", 3));
    CgIbiOptions io;
    io.alpha = N(a, "alpha", 0.2);
    io.ramp = N(a, "ramp", 0.1);
    io.pressure_correction = io.ramp > 0;
    const std::string dir = S(a, "o", "ibi");
    if (mode == "ibi-step") {
      if (S(a, "pairs").empty()) throw std::invalid_argument("ibi-step needs --pairs (the pairs.json the runs used)");
      const CgPairSet cur = cg_pairs_from_json(read_json_file(S(a, "pairs")));
      CgRdfAccumulator acc(types, po);
      std::vector<double> press;
      for (size_t k = 0; k < sys.size(); ++k) {
        const int si = acc.add_system(sys[k].t, 1, stem_of_map(sys[k].map));
        cg_frames(sys[k], size_t(std::max(1.0, N(a, "stride", 1))), [&](const std::vector<Vec3>& q, const Cell& c) { acc.add_frame(si, q, c); });
        if (sys[k].logs.empty()) {
          if (io.pressure_correction) throw std::invalid_argument(stem_of_map(sys[k].map) + ": give its LAMMPS log after its dump (for the pressure), or --ramp 0");
          press.push_back(tg.pressure[k]);
        } else press.push_back(log_mean(sys[k].logs.back(), "Press"));
      }
      std::vector<std::vector<CgRdf>> cur_rdf;
      for (size_t k = 0; k < sys.size(); ++k) cur_rdf.push_back(acc.rdf(int(k)));
      CgIbiStepReport rep;
      const CgPairSet next = ibi_step(cur, tg, cur_rdf, press, io, &rep);
      for (const auto& f : write_pairs(next, types, dir)) files.push_back(f);
      std::snprintf(b, sizeof b, "IBI iteration %d → %d · largest residual %.4f · pressure ramp A = %.4f kcal/mol\n", cur.iteration, next.iteration, rep.residual, rep.ramp_A);
      t << b;
      for (size_t k = 0; k < sys.size(); ++k) { std::snprintf(b, sizeof b, "  %s: P = %.0f atm (target %.0f)\n", tg.systems[k].c_str(), press[k], tg.pressure[k]); t << b; }
      Json pr = Json::array();
      for (const auto& x : rep.pairs) {
        std::snprintf(b, sizeof b, "  %-5s residual %.4f · largest |Δg| %.3f\n", x.key.c_str(), x.residual, x.max_dev);
        t << b;
        Json e = Json::object();
        e["key"] = x.key, e["residual"] = x.residual, e["max_dev"] = x.max_dev;
        pr.push_back(e);
      }
      r["pairs"] = pr;
      r["residual"] = rep.residual;
      Json pj = Json::array();
      for (double x : press) pj.push_back(x);
      r["pressure"] = pj;
    } else {
      // the loop with CAPS's engine: each system from its starting CG structure (MAP DATA)
      CgPairSet p = !S(a, "pairs").empty() ? cg_pairs_from_json(read_json_file(S(a, "pairs"))) : ibi_start(tg, io);
      std::unique_ptr<CgBondedResult> bonded;
      if (!S(a, "bonded").empty()) bonded = std::make_unique<CgBondedResult>(bonded_from_json(read_json_file((fs::path(S(a, "bonded")) / "bonded.json").string())));
      const auto masses = cg_masses(sys);
      std::vector<System> start;
      for (const auto& x : sys) {
        if (x.frames.empty()) throw std::invalid_argument("ibi-run: give each system's starting CG structure (STEM.cg.data) after its map");
        System d = read_lammps_data(x.frames.front());
        if (d.atoms.size() != x.t.beads()) throw std::invalid_argument(x.frames.front() + " does not match its map");
        start.push_back(d);
      }
      CgEngineOptions eo;
      eo.temperature = tg.temperature;
      eo.dt = N(a, "dt", 10);
      eo.steps = int64_t(N(a, "steps", 2000));
      eo.equilibrate = eo.steps / 4;
      eo.frame_every = std::max<int64_t>(1, eo.steps / 40);
      const int iters = int(N(a, "iterations", 3));
      for (int it = 0; it < iters; ++it) {
        std::vector<std::vector<CgRdf>> cur_rdf;
        std::vector<double> press;
        for (size_t k = 0; k < sys.size(); ++k) {
          const ForceField ff = cg_forcefield(sys[k].t, types, masses, p, bonded.get(), po.exclude);
          const CgEngineRun run = run_cg_engine(start[k], ff, eo);
          CgRdfAccumulator acc(types, po);
          const int si = acc.add_system(sys[k].t);
          for (size_t f = 0; f < run.frames.size(); ++f) acc.add_frame(si, run.frames[f], run.cells[f]);
          cur_rdf.push_back(acc.rdf(si));
          press.push_back(run.pressure);
          start[k].atoms = run.last.atoms, start[k].cell = run.last.cell;
          eo.new_velocities = false;
        }
        CgIbiStepReport rep;
        p = ibi_step(p, tg, cur_rdf, press, io, &rep);
        std::snprintf(b, sizeof b, "iteration %d: largest residual %.4f · ramp %.4f kcal/mol · P", p.iteration, rep.residual, rep.ramp_A);
        t << b;
        for (double x : press) { std::snprintf(b, sizeof b, " %.0f", x); t << b; }
        t << " atm\n";
        const std::string itdir = (fs::path(dir) / ("it" + std::string(p.iteration < 10 ? "00" : p.iteration < 100 ? "0" : "") + std::to_string(p.iteration))).string();
        for (const auto& f : write_pairs(p, types, itdir)) files.push_back(f);
      }
      if (bonded && !sys.front().t.dihedrals.empty()) t << "note: CAPS's engine runs without the dihedral tables (LAMMPS uses them): use the LAMMPS loop for the production model\n";
    }
  } else if (mode == "fit") {
    if (S(a, "pairs").empty()) throw std::invalid_argument("fit needs --pairs (pairs.json)");
    const CgPairSet p = cg_pairs_from_json(read_json_file(S(a, "pairs")));
    const std::string form = S(a, "form", "lj126");
    const auto fits = fit_pairs(p, form, N(a, "repulsion", 4));
    CgTypes types;
    if (!S(a, "types").empty()) types = cg_types_from_json(read_json_file(S(a, "types")));
    else { std::set<std::string> bs; for (const auto& x : p.pairs) { const auto d = x.key.find('-'); bs.insert(x.key.substr(0, d)), bs.insert(x.key.substr(d + 1)); } types.beads.assign(bs.begin(), bs.end()); }
    const std::string dir = S(a, "o", form);
    fs::create_directories(dir);
    write_text((fs::path(dir) / "pair.in").string(), lammps_pair_lines(fits, types, p.rc));
    files.push_back((fs::path(dir) / "pair.in").string());
    Json fj = Json::object();
    fj["format"] = "caps-cg-fits";
    fj["form"] = form;
    fj["rc"] = p.rc;
    fj["temperature"] = p.temperature;
    Json arr = Json::array();
    for (const auto& f : fits) {
      Json e = Json::object();
      e["key"] = f.key, e["epsilon"] = f.epsilon, e["sigma"] = f.sigma, e["a"] = f.a, e["n"] = f.n, e["m"] = f.m;
      e["rms"] = std::isfinite(f.rms) ? Json(f.rms) : Json();
      e["lo"] = f.lo, e["hi"] = f.hi;
      arr.push_back(e);
      if (!std::isfinite(f.rms)) { t << "  " << f.key << ": too few points to fit\n"; continue; }
      if (form == "morse") std::snprintf(b, sizeof b, "  %-5s D0 %.4f kcal/mol · α %.4f /Å · r0 %.3f Å · rms %.3f kcal/mol over %.2f–%.2f Å\n", f.key.c_str(), f.epsilon, f.a, f.sigma, f.rms, f.lo, f.hi);
      else if (form == "mie") std::snprintf(b, sizeof b, "  %-5s ε %.4f kcal/mol · σ %.3f Å · n %.2f, m 6 · rms %.3f kcal/mol over %.2f–%.2f Å\n", f.key.c_str(), f.epsilon, f.sigma, f.n, f.rms, f.lo, f.hi);
      else std::snprintf(b, sizeof b, "  %-5s ε %.4f kcal/mol · σ %.3f Å · rms %.3f kcal/mol over %.2f–%.2f Å\n", f.key.c_str(), f.epsilon, f.sigma, f.rms, f.lo, f.hi);
      t << b;
    }
    fj["fits"] = arr;
    write_text((fs::path(dir) / "fits.json").string(), fj.dump(1) + "\n");
    files.push_back((fs::path(dir) / "fits.json").string());
    for (const auto& f : write_pairs(tables_of_fits(fits, p), types, (fs::path(dir) / "tables").string())) files.push_back(f);
    t << form << " fitted to " << fits.size() << " pairs: " << dir << "/pair.in (analytic), " << dir << "/tables (the same as tables)\n";
  } else if (mode == "tg") {
    if (in.size() < 2) throw std::invalid_argument("tg: give the cooling scan's file (T and density per line)");
    std::ifstream f(in[1]);
    if (!f) throw std::invalid_argument("cannot open " + in[1]);
    std::vector<double> Tv, rho;
    for (std::string line; std::getline(f, line);) {
      if (line.empty() || line[0] == '#') continue;
      std::istringstream ss(line);
      double x, y;
      if (ss >> x >> y) Tv.push_back(x), rho.push_back(1.0 / y);   // specific volume
    }
    if (Tv.size() < 5) throw std::invalid_argument("tg: at least five temperatures are needed");
    const BilinearFit fit = fit_bilinear(Tv, rho);
    if (!fit.ok) throw std::invalid_argument("tg: no kink in the specific volume (" + fit.note + ")");
    std::snprintf(b, sizeof b, "T_g = %.1f ± %.1f K (bilinear fit of the specific volume; α %.2e → %.2e /K below → above)\n", fit.tg, fit.tg_err, fit.alpha_low, fit.alpha_high);
    t << b;
    r["tg"] = fit.tg;
    r["tg_err"] = fit.tg_err;
  } else if (mode == "calibrate") {
    if (S(a, "fits").empty()) throw std::invalid_argument("calibrate needs --fits (fits.json)");
    const Json fj = read_json_file(S(a, "fits"));
    std::vector<CgPairFit> fits;
    for (const auto& e : fj["fits"].items()) {
      CgPairFit f;
      f.key = e["key"].str(), f.form = fj["form"].str(), f.epsilon = e.num("epsilon", 0), f.sigma = e.num("sigma", 0), f.a = e.num("a", 0), f.n = e.num("n", 12), f.m = e.num("m", 6);
      f.rms = e["rms"].is_null() ? std::numeric_limits<double>::quiet_NaN() : e["rms"].number();
      fits.push_back(f);
    }
    const std::string hist = S(a, "history", "calibration.json");
    Json h = fs::exists(hist) ? read_json_file(hist) : Json::object();
    if (!h.has("runs")) h["runs"] = Json::array();
    if (a.has("density") || a.has("tg")) {
      Json run = Json::object();
      run["s_sigma"] = N(a, "s_sigma", 1), run["s_eps"] = N(a, "s_eps", 1), run["density"] = N(a, "density", 0), run["tg"] = N(a, "tg", 0);
      h["runs"].push_back(run);
    }
    std::vector<CgCalibrationPoint> pts;
    for (const auto& e : h["runs"].items()) pts.push_back({e.num("s_sigma", 1), e.num("s_eps", 1), e.num("density", 0), e.num("tg", 0)});
    const double rho_t = N(a, "target_density", h.num("target_density", 0)), tg_t = N(a, "target_tg", h.num("target_tg", 0));
    h["target_density"] = rho_t, h["target_tg"] = tg_t;
    fs::create_directories(fs::path(hist).parent_path().empty() ? fs::path(".") : fs::path(hist).parent_path());
    write_text(hist, h.dump(1) + "\n");
    files.push_back(hist);
    CgCalibrationStep st;
    if (pts.empty()) { st.s_sigma = 1, st.s_eps = 1, st.note = "no run yet: the decks for s_σ = s_ε = 1 (the fitted set)"; }
    else st = calibration_step(pts, rho_t, tg_t);
    t << st.note << "\n";
    const std::string dir = S(a, "o", "calibration");
    CgTypes types;
    if (!S(a, "types").empty()) types = cg_types_from_json(read_json_file(S(a, "types")));
    else { std::set<std::string> bs; for (const auto& f : fits) { const auto d = f.key.find('-'); bs.insert(f.key.substr(0, d)), bs.insert(f.key.substr(d + 1)); } types.beads.assign(bs.begin(), bs.end()); }
    fs::create_directories(dir);
    write_text((fs::path(dir) / "pair.in").string(), lammps_pair_lines(scale_fits(fits, st.s_sigma, st.s_eps), types, fj.num("rc", 15) * st.s_sigma));
    files.push_back((fs::path(dir) / "pair.in").string());
    std::vector<std::string> deckfiles;
    write_cg_decks(dir, int(N(a, "exclude", 3)), fj.num("temperature", 300), N(a, "dt", 10), int64_t(N(a, "steps", 500000)), N(a, "t_hi", 500), N(a, "t_lo", 150), N(a, "t_step", 25), deckfiles);
    for (const auto& f : deckfiles) files.push_back(f);
    r["s_sigma"] = st.s_sigma, r["s_eps"] = st.s_eps, r["density_done"] = st.density_done, r["tg_done"] = st.tg_done;
    std::snprintf(b, sizeof b, "next: run %s/in.cg_run with -var ENS npt (density) and in.cg_tg (T_g) using -var PAIR %s, then report them with --s-sigma %.4f --s-eps %.4f\n",
                  dir.c_str(), dir.c_str(), st.s_sigma, st.s_eps);
    t << b;
  } else {
    throw std::invalid_argument("cgfit " + mode + ": bonded, refine, targets, ibi-start, ibi-step, ibi-run, fit, tg or calibrate");
  }
  t << r.text("command") << "\n";
  r["text"] = t.str();
}
}  // namespace

const std::vector<std::string>& cg_commands() {
  static std::vector<std::string> c;
  if (c.empty()) for (const auto& x : table()) c.push_back(x.name);
  return c;
}
bool is_cg_command(const std::string& c) { return std::find(cg_commands().begin(), cg_commands().end(), c) != cg_commands().end(); }
bool cg_is_switch(const std::string& c, std::string name) {
  while (!name.empty() && name[0] == '-') name.erase(0, 1);
  std::replace(name.begin(), name.end(), '-', '_');
  if (name == "json" || name == "help") return true;
  for (const auto& x : table())
    if (c == x.name) for (const auto& o : x.opts) if (name == o.name) return std::string(o.def) == "false";
  return false;
}

std::string cg_help(const std::string& c) {
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

Json cg_run(const std::string& c, const Json& a, const std::string& data_dir) {
  cmd_of(c);
  const auto in = inputs(a);
  Json r = Json::object();
  r["ok"] = true;
  r["command"] = dft_command_line(c, a);
  Json files = Json::array();
  if (c == "cgmap") {
    if (in.empty()) throw std::invalid_argument("give the all-atom structure(s) to map");
    const std::string dir = S(a, "o", ".");
    fs::create_directories(dir);
    std::string label;
    const bool from_map = !S(a, "map").empty();
    if (from_map && in.size() > 1) throw std::invalid_argument("--map belongs to one structure: map one at a time");
    CgRules rules;
    if (!from_map) rules = rules_of(a, data_dir, &label);
    // map every structure first, so the shared list covers them all
    std::vector<System> aa;
    std::vector<CgMapping> maps;
    for (const auto& p : in) {
      aa.push_back(open_file(p).frame(0));
      maps.push_back(from_map ? cg_mapping_from_json(read_json_file(S(a, "map")), aa.back()) : cg_mapping(aa.back(), rules));
    }
    CgTypes types;
    std::string types_file = S(a, "types");
    if (!types_file.empty()) {
      types = cg_types_from_json(read_json_file(types_file));
    } else {
      std::vector<CgTypes> lists;
      for (const auto& m : maps) lists.push_back(cg_types_of(m));
      types = cg_types_union(lists);
      types_file = (fs::path(dir) / "types.json").string();
      write_text(types_file, cg_types_json(types).dump(1) + "\n");
      files.push_back(types_file);
    }
    r["types"] = cg_types_json(types);
    Json systems = Json::array();
    const std::vector<std::string> cites;
    for (size_t k = 0; k < in.size(); ++k) {
      const auto& m = maps[k];
      const std::string stem = stem_of(in[k]);
      const std::string data = (fs::path(dir) / (stem + ".cg.data")).string();
      const std::string mapf = (fs::path(dir) / (stem + ".map.json")).string();
      std::vector<Vec3> pos;
      for (const auto& at : aa[k].atoms) pos.push_back(at.pos);
      const System cg = cg_structure(m, aa[k], cg_positions(m, aa[k], pos, aa[k].cell), aa[k].cell, types);
      write_cg_lammps_data(m, cg, types, data);
      write_text(mapf, cg_mapping_json(m, aa[k]).dump(1) + "\n");
      files.push_back(data);
      files.push_back(mapf);
      const auto sum = cg_mapping_summary(m);
      Json e = Json::object();
      e["input"] = in[k];
      e["atoms"] = int(aa[k].atoms.size());
      e["beads"] = int(m.beads());
      e["chains"] = sum.chains;
      e["beads_per_chain"] = std::to_string(sum.min_beads) + (sum.max_beads != sum.min_beads ? "–" + std::to_string(sum.max_beads) : "");
      // repeat units (the residues a grown cell numbers), when the atoms carry them
      std::set<std::pair<int64_t, int64_t>> units;
      for (const auto& at : aa[k].atoms) if (at.resid > 0) units.insert({at.mol, at.resid});
      if (!units.empty()) e["beads_per_unit"] = double(m.beads()) / double(units.size());
      Json kinds = Json::object(), frac = Json::object();
      for (const auto& [kk, n] : sum.kinds) kinds[kk] = n, frac[kk] = sum.fraction.at(kk);
      e["kinds"] = kinds;
      e["fractions"] = frac;
      e["mass_aa"] = m.mass_aa;
      e["mass_error"] = sum.mass_error;
      Json seq = Json::array();
      for (const auto& q : sum.sequences) seq.push_back(q);
      e["sequences"] = seq;
      Json cls = Json::array();
      for (const auto& cl : m.classes) {
        Json x = Json::object();
        x["name"] = cl.name;
        x["formula"] = cl.formula;
        x["count"] = cl.count;
        x["key"] = cl.key.substr(0, 8);
        cls.push_back(x);
      }
      e["classes"] = cls;
      Json notes = Json::array();
      for (const auto& nn : m.notes) notes.push_back(nn);
      e["notes"] = notes;
      e["data"] = data;
      e["map"] = mapf;
      // the trajectory (one structure)
      if (!S(a, "dump").empty()) {
        if (in.size() > 1) throw std::invalid_argument("--dump maps one structure's trajectory: give one structure");
        CgTrajectoryOptions to;
        to.types = types;
        to.stride = size_t(std::max(1.0, N(a, "stride", 1)));
        to.first = size_t(std::max(1.0, N(a, "first", 1))) - 1;
        if (a.has("last") && !S(a, "last").empty()) to.last = size_t(std::max(1.0, N(a, "last", 1))) - 1;
        const std::string outd = S(a, "out_dump", (fs::path(dir) / (stem + ".cg.lammpstrj")).string());
        const auto tr = map_lammps_dump(m, aa[k], S(a, "dump"), outd, to);
        Json t = Json::object();
        t["dump"] = S(a, "dump");
        t["out"] = outd;
        t["frames_read"] = double(tr.frames_read);
        t["frames_written"] = double(tr.frames_written);
        t["seconds"] = tr.seconds;
        t["unwrapped_input"] = tr.unwrapped_input;
        Json tn = Json::array();
        for (const auto& nn : tr.notes) tn.push_back(nn);
        t["notes"] = tn;
        e["trajectory"] = t;
        files.push_back(outd);
        provenance(outd, "cg.map_trajectory", std::to_string(tr.frames_written) + " frames mapped to " + std::to_string(m.beads()) + " beads",
                   {{"mapping", fs::path(mapf).filename().string()}, {"stride", std::to_string(to.stride)}, {"first", std::to_string(to.first + 1)}},
                   {{fs::path(in[k]).filename().string(), sha_of(in[k])}, {fs::path(S(a, "dump")).filename().string(), sha_of(S(a, "dump"))}}, cites);
      }
      KeyValues params = {{"mapping", from_map ? "map " + fs::path(S(a, "map")).filename().string() : label}, {"position", m.position},
                          {"beads", std::to_string(m.beads())}, {"bead types", std::to_string(types.beads.size())}};
      for (const auto& x : rules.cut) params.push_back({"cut", x});
      for (const auto& x : rules.names) params.push_back({"name " + x.name, x.smarts});
      provenance(data, "cg.map", std::to_string(m.atoms) + " atoms → " + std::to_string(m.beads()) + " beads", params,
                 {{fs::path(in[k]).filename().string(), sha_of(in[k])}}, cites,
                 {{"resolution", "coarse-grained beads at the " + std::string(m.position == "cog" ? "centre of geometry" : m.position == "com" ? "centre of mass" : "anchor atom") + " of their atoms"}});
      systems.push_back(e);
    }
    r["systems"] = systems;
    std::ostringstream t;
    char b[256];
    for (const auto& e : systems.items()) {
      std::snprintf(b, sizeof b, "%s: %d atoms → %d beads · %d chains × %s beads", e.text("input").c_str(), int(e.num("atoms", 0)), int(e.num("beads", 0)),
                    int(e.num("chains", 0)), e.text("beads_per_chain").c_str());
      t << b;
      if (e.has("beads_per_unit")) { std::snprintf(b, sizeof b, " · %.3g per repeat unit", e["beads_per_unit"].number()); t << b; }
      std::snprintf(b, sizeof b, " · mass conserved to %.1e g/mol\n   ", e.num("mass_error", 0));
      t << b;
      for (const auto& [k, v] : e["fractions"].members()) { std::snprintf(b, sizeof b, " %s %.1f %%", k.c_str(), 100 * v.number()); t << b; }
      t << "\n";
      for (const auto& q : e["sequences"].items()) t << "    " << (q.str().size() > 90 ? q.str().substr(0, 90) + " …" : q.str()) << "\n";
      for (const auto& nn : e["notes"].items()) t << "    note: " << nn.str() << "\n";
      if (e.has("trajectory")) {
        const Json& tr = e["trajectory"];
        std::snprintf(b, sizeof b, "    trajectory: %d frames → %s (%.1f s)\n", int(tr.num("frames_written", 0)), tr.text("out").c_str(), tr.num("seconds", 0));
        t << b;
      }
    }
    t << "types: " << types.beads.size() << " bead, " << types.bonds.size() << " bond, " << types.angles.size() << " angle, " << types.dihedrals.size() << " dihedral\n";
    t << "wrote " << files.size() << " files in " << dir << "\n" << r.text("command") << "\n";
    r["text"] = t.str();
  }
  if (c == "cgfit") {
    if (in.empty()) throw std::invalid_argument("cgfit bonded | refine, then the maps and their frames");
    const std::string mode = in[0];
    if (mode != "bonded" && mode != "refine") {
      cgfit_nonbonded(mode, a, in, r, files);
      r["files"] = files;
      return r;
    }
    // systems: each map.json followed by its frames
    struct SysIn { std::string map; std::vector<std::string> frames; Json mj; CgTopology t; };
    std::vector<SysIn> sys;
    for (size_t k = 1; k < in.size(); ++k) {
      const std::string& p = in[k];
      if (p.size() >= 9 && p.compare(p.size() - 9, 9, ".map.json") == 0) {
        SysIn x;
        x.map = p;
        x.mj = read_json_file(p);
        x.t = cg_topology_from_map(x.mj);
        sys.push_back(std::move(x));
      } else {
        if (sys.empty()) throw std::invalid_argument("give a mapping (STEM.map.json) before its frames: " + p);
        sys.back().frames.push_back(p);
      }
    }
    if (sys.empty()) throw std::invalid_argument("give at least one mapping (STEM.map.json) and its frames");
    CgTypes types;
    if (!S(a, "types").empty()) types = cg_types_from_json(read_json_file(S(a, "types")));
    else {
      std::vector<CgTypes> lists;
      for (const auto& x : sys) {
        CgTypes t;
        std::set<std::string> b, bo, an, di;
        for (const auto& k : x.t.kind) b.insert(k);
        for (const auto& [i, j] : x.t.bonds) bo.insert(cg_key({x.t.kind[size_t(i)], x.t.kind[size_t(j)]}));
        for (const auto& q : x.t.angles) an.insert(cg_key({x.t.kind[size_t(q[0])], x.t.kind[size_t(q[1])], x.t.kind[size_t(q[2])]}));
        for (const auto& q : x.t.dihedrals) di.insert(cg_key({x.t.kind[size_t(q[0])], x.t.kind[size_t(q[1])], x.t.kind[size_t(q[2])], x.t.kind[size_t(q[3])]}));
        lists.push_back({{b.begin(), b.end()}, {bo.begin(), bo.end()}, {an.begin(), an.end()}, {di.begin(), di.end()}});
      }
      types = cg_types_union(lists);
    }
    CgBondedOptions o;
    o.temperature = N(a, "T", 300);
    o.threshold = N(a, "threshold", o.threshold);
    o.smooth = N(a, "smooth", o.smooth);
    o.bond_bin = N(a, "bin_bond", o.bond_bin);
    o.angle_bin = N(a, "bin_angle", o.angle_bin);
    o.dihedral_bin = N(a, "bin_dihedral", o.dihedral_bin);
    o.bond_wall = N(a, "bond_wall", o.bond_wall);
    o.angle_wall = N(a, "angle_wall", o.angle_wall);
    if (a.has("aat")) {
      const auto v = L(a, "aat");
      if (v.size() != 3) throw std::invalid_argument("--aat takes K,θ1,θ2");
      o.aat_k = std::stod(v[0]), o.aat_theta1 = std::stod(v[1]), o.aat_theta2 = std::stod(v[2]);
    }
    std::vector<double> w(sys.size(), 1.0);
    if (a.has("weights")) {
      const auto v = L(a, "weights");
      if (v.size() != sys.size()) throw std::invalid_argument("--weights needs " + std::to_string(sys.size()) + " values, one per system");
      for (size_t k = 0; k < v.size(); ++k) w[k] = std::stod(v[k]);
    }
    const size_t stride = size_t(std::max(1.0, N(a, "stride", 1)));
    CgBondedAccumulator acc(types, o);
    Json sysrep = Json::array();
    KeyValues ins;
    for (size_t k = 0; k < sys.size(); ++k) {
      const auto& x = sys[k];
      const int si = acc.add_system(x.t, w[k], fs::path(x.map).filename().string());
      ins.push_back({fs::path(x.map).filename().string(), sha_of(x.map)});
      size_t nfr = 0;
      for (const auto& fpath : x.frames) {
        ins.push_back({fs::path(fpath).filename().string(), sha_of(fpath)});
        if (fs::path(fpath).extension() == ".data") {
          const System d = read_lammps_data(fpath);
          if (d.atoms.size() != x.t.beads()) throw std::invalid_argument(fpath + " has " + std::to_string(d.atoms.size()) + " beads, its map " + std::to_string(x.t.beads()));
          std::vector<Vec3> p(d.atoms.size());
          for (size_t b = 0; b < d.atoms.size(); ++b) p[size_t(d.atoms[b].id > 0 ? d.atoms[b].id - 1 : int64_t(b))] = d.atoms[b].pos;
          acc.add_frame(si, p, d.cell);
          ++nfr;
        } else {
          nfr += for_each_dump_frame(fpath, x.t.beads(), [&](size_t, int64_t, const std::vector<Vec3>& p, const Cell& cell, bool) {
            acc.add_frame(si, p, cell);
            return true;
          }, stride);
        }
      }
      if (!nfr) throw std::invalid_argument("no frames for " + x.map + ": give its CG dump(s) or STEM.cg.data after it");
      Json e = Json::object();
      e["map"] = x.map;
      e["frames"] = double(nfr);
      e["weight"] = w[k];
      sysrep.push_back(e);
    }
    const std::string dir = S(a, "o", mode == "bonded" ? "bonded" : "bonded_refined");
    CgBondedResult res;
    std::map<std::string, double> residual;
    if (mode == "bonded") {
      res = invert_bonded(acc);
    } else {
      if (S(a, "tables").empty()) throw std::invalid_argument("refine needs --tables (the bonded.json the CG run used)");
      res = refine_bonded(bonded_from_json(read_json_file(S(a, "tables"))), acc, N(a, "alpha", 0.5), &residual);
    }
    for (const auto& f : write_bonded(res, types, o, dir)) files.push_back(f);
    write_text((fs::path(dir) / "types.json").string(), cg_types_json(types).dump(1) + "\n");
    files.push_back((fs::path(dir) / "types.json").string());
    // the report: per table its minimum, curvature, sampled range and how well the halves agree
    Json tabs = Json::array();
    std::ostringstream t;
    char b[256];
    std::snprintf(b, sizeof b, "%s: %zu systems, %zu frames, %.0f K\n", mode == "bonded" ? "Boltzmann inversion" : "bonded IBI step", sys.size(), acc.frames(), o.temperature);
    t << b;
    for (const auto* list : {&res.bonds, &res.angles, &res.dihedrals})
      for (const auto& T : *list) {
        Json e = Json::object();
        e["kind"] = T.kind == 0 ? "bond" : T.kind == 1 ? "angle" : "dihedral";
        e["key"] = T.key;
        e["sampled"] = T.sampled;
        e["count"] = double(T.count);
        e["x0"] = T.x0;
        e["k_harmonic"] = T.k_harmonic;
        e["lo"] = T.lo;
        e["hi"] = T.hi;
        if (!std::isnan(T.half_diff)) e["half_diff_kT"] = T.half_diff;
        const std::string rk = std::string(e.text("kind")) + " " + T.key;
        if (residual.count(rk)) e["residual"] = residual[rk];
        tabs.push_back(e);
        if (!T.sampled) { t << "  " << e.text("kind") << " " << T.key << ": no samples\n"; continue; }
        const char* u = T.kind == 0 ? "Å" : "°";
        char hd[64];
        if (std::isnan(T.half_diff)) std::snprintf(hd, sizeof hd, "no split check");
        else std::snprintf(hd, sizeof hd, "halves differ by %.2f kT", T.half_diff);
        std::snprintf(b, sizeof b, "  %-8s %-10s min %.3g %s · range %.3g–%.3g %s · %s · %ld samples%s\n", e.text("kind").c_str(), T.key.c_str(), T.x0, u,
                      T.lo, T.hi, u, hd, T.count, residual.count(rk) ? (" · residual " + std::to_string(residual[rk]).substr(0, 5)).c_str() : "");
        t << b;
      }
    // a time step for these tables: a thirtieth of the fastest bond oscillation, the stiffer of the well and the wall
    // (masses: each bond type's lightest bead pair from the maps)
    {
      double dt = 1e30;
      std::string by;
      for (const auto& T : res.bonds) {
        if (!T.sampled) continue;
        const auto dash = T.key.find('-');
        double mu = 1e30;
        for (const auto& x : sys) {
          std::map<std::string, double> mk;
          for (const auto& bd : x.mj["beads"].items()) {
            const std::string k = bd["kind"].str();
            const double m = bd.num("mass", 0);
            if (m > 0 && (!mk.count(k) || m < mk[k])) mk[k] = m;
          }
          const std::string ka = T.key.substr(0, dash), kb = T.key.substr(dash + 1);
          if (mk.count(ka) && mk.count(kb)) mu = std::min(mu, mk[ka] * mk[kb] / (mk[ka] + mk[kb]));
        }
        if (mu >= 1e30) continue;
        const double k = std::max(T.k_harmonic, o.bond_wall);   // E = k x²: spring constant 2k
        const double period = 2 * 3.14159265358979 * std::sqrt(mu / (2 * k * 418.4)) * 1000;   // fs (kcal/mol/Å², g/mol)
        if (period / 30 < dt) dt = period / 30, by = T.key;
      }
      if (dt < 1e30) {
        const double step = std::max(0.5, std::floor(dt * 2) / 2);
        r["timestep_fs"] = step;
        std::snprintf(b, sizeof b, "  time step: about %.1f fs (a thirtieth of the %s bond's period with its wall); check the energy drift in a short NVE run with the non-bonded terms on\n", step, by.c_str());
        t << b;
      }
    }
    for (const auto& nn : res.notes) t << "  note: " << nn << "\n";
    t << "wrote " << dir << " (bonded.in: include after read_data with -var BONDED " << dir << ")\n" << r.text("command") << "\n";
    r["text"] = t.str();
    r["tables"] = tabs;
    r["systems"] = sysrep;
    r["dir"] = dir;
    provenance((fs::path(dir) / "bonded.in").string(), mode == "bonded" ? "cg.bonded_bi" : "cg.bonded_ibi",
               std::string(mode == "bonded" ? "bonded Boltzmann inversion" : "bonded IBI step") + " over " + std::to_string(acc.frames()) + " frames of " + std::to_string(sys.size()) + " systems",
               {{"temperature", std::to_string(o.temperature)}, {"threshold", std::to_string(o.threshold)}, {"smooth (bins)", std::to_string(o.smooth)},
                {"bins", std::to_string(o.bond_bin) + " Å, " + std::to_string(o.angle_bin) + "°, " + std::to_string(o.dihedral_bin) + "°"}},
               ins, {"reith2003"}, {{"bonded potentials", "tabulated, from mapped all-atom distributions at one temperature"}});
  }
  r["files"] = files;
  return r;
}

}  // namespace caps
