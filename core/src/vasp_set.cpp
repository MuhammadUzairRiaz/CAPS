// VASP calculation sets, job scripts and job tools (caps/vasp_set.hpp).
#include "caps/vasp_set.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <map>
#include <regex>
#include <set>
#include <sstream>
#include <stdexcept>

#include "caps/elements.hpp"
#include "caps/slab2d.hpp"
#include "caps/vasp_out.hpp"

namespace caps {

namespace fs = std::filesystem;

namespace {

std::string slurp(const std::string& p) {
  std::ifstream in(p, std::ios::binary);
  if (!in) return "";
  std::stringstream ss;
  ss << in.rdbuf();
  return ss.str();
}
void spit(const std::string& p, const std::string& t) {
  std::ofstream f(p, std::ios::binary);
  if (!f) throw std::runtime_error("cannot write " + p);
  f << t;
}
std::string join(const std::string& a, const std::string& b) { return (fs::path(a) / b).string(); }
std::string fmt(const char* f, double v) { char b[96]; std::snprintf(b, sizeof b, f, v); return b; }
Json read_json_file(const std::string& path) {
  const std::string t = slurp(path);
  if (t.empty()) throw std::runtime_error("cannot read " + path);
  return Json::parse(t);
}
// a tag's line replaced (or appended)
std::string set_tag(const std::string& incar, const std::string& key, const std::string& value) {
  const std::regex re("(^|\\n)" + key + "\\s*=[^\\n]*");
  if (std::regex_search(incar, re)) return std::regex_replace(incar, re, "$1" + key + std::string(std::max<size_t>(1, 8 - key.size()), ' ') + "= " + value, std::regex_constants::format_first_only);
  return incar + (incar.empty() || incar.back() == '\n' ? "" : "\n") + key + std::string(std::max<size_t>(1, 8 - key.size()), ' ') + "= " + value + "\n";
}
std::string remove_tag(const std::string& incar, const std::string& key) { return std::regex_replace(incar, std::regex("(^|\\n)" + key + "\\s*=[^\\n]*"), ""); }
std::string set_stages(const std::string& job, const std::string& stages) { return std::regex_replace(job, std::regex("STAGES=\"[^\"]*\""), "STAGES=\"" + stages + "\"", std::regex_constants::format_first_only); }
bool is_d_element(int z) { return (z >= 21 && z <= 30) || (z >= 39 && z <= 48) || (z >= 72 && z <= 80) || z == 57; }
void make_exec(const std::string& p) {
  std::error_code ec;
  fs::permissions(p, fs::perms::owner_exec | fs::perms::group_exec | fs::perms::others_exec, fs::perm_options::add, ec);
}
std::vector<std::string> species_of(const System& s) {
  std::vector<std::string> out;
  for (const auto& a : s.atoms) { const std::string e = element(a.element).symbol; if (out.empty() || out.back() != e) out.push_back(e); }
  return out;
}

}  // namespace

// ---------------------------------------------------------------- profiles
// data/dft/clusters.json (generic starting points), then the user's own ~/CAPS/dft-clusters.json (kept on their machine):
// a profile there replaces one of the same name
std::string user_clusters_file() {
  const char* home = std::getenv("HOME");
#ifdef _WIN32
  if (!home) home = std::getenv("USERPROFILE");
#endif
  return home ? std::string(home) + "/CAPS/dft-clusters.json" : std::string();
}
std::vector<ClusterProfile> cluster_profiles(const std::string& data_dir) {
  std::vector<ClusterProfile> out;
  std::vector<Json> files = {read_json_file(data_dir + "/dft/clusters.json")};
  if (const std::string u = user_clusters_file(); !u.empty() && std::filesystem::exists(u)) files.push_back(read_json_file(u));
  for (const Json& j : files) {
  if (!j.has("profiles")) continue;
  for (const auto& e : j["profiles"].items()) {
    ClusterProfile p;
    p.name = e.text("name"), p.about = e.text("about"), p.scheduler = e.text("scheduler", "slurm"), p.account = e.text("account");
    p.licenses = e.text("licenses"), p.modules = e.text("modules"), p.vasp = e.text("vasp", "vasp_std");
    p.cores_per_node = int(e.num("cores_per_node", 0)), p.mem_per_cpu = int(e.num("mem_per_cpu", 0));
    if (e.has("workspace")) p.ws_allocate = e["workspace"].text("allocate"), p.ws_find = e["workspace"].text("find"), p.ws_release = e["workspace"].text("release");
    for (const auto& r : e["parallel"].items()) p.parallel.push_back({int(r.num("max_atoms", 0)), int(r.num("ranks", 1)), int(r.num("kpar", 1)), int(r.num("ncore", 1)), r.text("time", "24:00:00")});
    if (e.has("gamma")) { const Json& g = e["gamma"]; p.gamma = {0, int(g.num("ranks", 1)), int(g.num("kpar", 1)), int(g.num("ncore", 1)), g.text("time", "24:00:00")}; }
    if (e.has("hang")) for (const auto& h : e["hang"].items()) p.hang.push_back({int(h.num("max_atoms", 0)), int(h.num("seconds", 1800))});
    const auto same = std::find_if(out.begin(), out.end(), [&](const ClusterProfile& q) { return q.name == p.name; });
    if (same != out.end()) *same = p;
    else out.push_back(p);
  }
  }
  return out;
}
ClusterProfile cluster_profile(const std::string& data_dir, const std::string& name) {
  for (const auto& p : cluster_profiles(data_dir)) if (p.name == name) return p;
  throw std::invalid_argument("no cluster profile " + name + " (data/dft/clusters.json, ~/CAPS/dft-clusters.json)");
}
ClusterRule parallel_rule(const ClusterProfile& p, int atoms, bool gamma_only) {
  if (gamma_only && p.gamma.ranks > 0) return p.gamma;
  for (const auto& r : p.parallel) if (r.max_atoms == 0 || atoms <= r.max_atoms) return r;
  return p.parallel.empty() ? ClusterRule{0, 1, 1, 1, "24:00:00"} : p.parallel.back();
}
int hang_seconds(const ClusterProfile& p, int atoms) {
  for (const auto& [m, s] : p.hang) if (m == 0 || atoms <= m) return s;
  return 1800;
}

// ---------------------------------------------------------------- INCAR
std::string incar_why(const std::string& k) {
  static const std::map<std::string, std::string> w = {
      {"ENCUT", "1.3 × the largest ENMAX of the POTCARs: high enough for a cell relaxation (Pulay stress)"},
      {"PREC", "Accurate: FFT grids fine enough for forces and stress"},
      {"ALGO", "Normal (blocked Davidson): robust for metals and molecules alike"},
      {"EDIFF", "1E-6 slab relaxation, 1E-5 large adsorption relaxations, 1E-7 single points"},
      {"NELM", "120 electronic steps before VASP gives up on a step"},
      {"ISMEAR", "0 (Gaussian): valid for a metallic slab AND a molecule; report E0 (σ → 0)"},
      {"SIGMA", "0.05 eV smearing width"},
      {"GGA", "PE: the PBE functional"},
      {"IVDW", "12: DFT-D3 with Becke–Johnson damping (dispersion binds a physisorbed molecule)"},
      {"LASPH", "non-spherical contributions inside the PAW spheres (d elements, GGA)"},
      {"LMAXMIX", "4 for d elements (charge-density mixing up to l = 4), else 2"},
      {"ISPIN", "1: non-magnetic (test ISPIN = 2 once for a new material)"},
      {"LREAL", ".FALSE. up to 40 atoms, else Auto — one value for every run of an adsorption set"},
      {"NCORE", "cores per orbital: from the cluster profile"},
      {"KPAR", "k-point groups: from the cluster profile"},
      {"IBRION", "2 conjugate gradient (relaxations), −1 no ionic step (single point), 0 MD, 5 finite differences"},
      {"ISIF", "3 relaxes ions and cell, 2 ions only"},
      {"LATTICE_CONSTRAINTS", ".TRUE. .TRUE. .FALSE.: a and b relax, c (the vacuum) stays (VASP ≥ 6.5)"},
      {"NSW", "the most ionic steps of the stage"},
      {"EDIFFG", "−0.01: stop when every force is below 0.01 eV/Å"},
      {"LWAVE", "no WAVECAR (large, never used)"},
      {"LCHARG", "CHGCAR only where Bader or a charge-density difference needs it"},
      {"LAECHG", "AECCAR0/2 for Bader (with LCHARG)"},
      {"LORBIT", "11: projected DOS per atom and orbital (DOSCAR, PROCAR)"},
      {"NEDOS", "3001 energies for a smooth DOS"},
      {"LVHAR", "LOCPOT with the electrostatic (Hartree + ionic) potential: the work function"},
      {"IDIPOL", "3: dipole correction along z, ENERGY only — no LDIPOL: on low-work-function surfaces (OH) LDIPOL lets "
                 "electrons leak into the vacuum and the energy runs away"},
      {"DIPOL", "the centre of the atoms (fractional): where the dipole layer sits"},
      {"MDALGO", "5: CSVR (Bussi) thermostat"},
      {"POTIM", "the MD step (fs), 1 fs with deuterium masses on H; 0.015 Å displacements for frequencies"},
      {"POMASS", "masses per species; H as deuterium (2.014) for a 1 fs MD step"},
      {"ISYM", "0: no symmetry (MD, frequencies)"},
      {"NFREE", "2: central differences"}};
  const auto it = w.find(k);
  return it == w.end() ? "" : it->second;
}

std::vector<IncarLine> incar_lines(const std::string& kind, int n, bool gamma, double encut, double dipz, const VaspSettings& s, const ClusterRule& par,
                                   const std::vector<int>& elements) {
  std::vector<IncarLine> L;
  auto add = [&](const std::string& k, const std::string& v) { L.push_back({k, v, incar_why(k)}); };
  add("ENCUT", fmt("%g", encut));
  add("PREC", s.prec);
  add("ALGO", s.algo);
  add("EDIFF", kind == "static" ? s.ediff_static : s.ediff_relax);
  add("NELM", std::to_string(s.nelm));
  add("ISMEAR", std::to_string(s.ismear));
  add("SIGMA", fmt("%g", s.sigma));
  add("GGA", s.gga);
  if (s.ivdw) add("IVDW", std::to_string(s.ivdw));
  add("LASPH", s.lasph ? ".TRUE." : ".FALSE.");
  bool d = false;
  for (int z : elements) d |= is_d_element(z);
  add("LMAXMIX", std::to_string(s.lmaxmix > 0 ? s.lmaxmix : (d ? 4 : 2)));
  add("ISPIN", std::to_string(s.ispin));
  add("LREAL", !s.lreal.empty() ? s.lreal : (n <= 40 ? ".FALSE." : "Auto"));
  add("NCORE", std::to_string(par.ncore));
  add("KPAR", std::to_string(par.kpar));
  if (kind == "cell") {
    add("IBRION", "2"); add("ISIF", "3"); add("LATTICE_CONSTRAINTS", ".TRUE. .TRUE. .FALSE.");
    add("NSW", std::to_string(s.nsw_cell)); add("EDIFFG", "-0.01"); add("LWAVE", ".FALSE."); add("LCHARG", ".FALSE.");
  } else if (kind == "relax") {
    add("IBRION", "2"); add("ISIF", "2");
    add("NSW", std::to_string(s.nsw_relax > 0 ? s.nsw_relax : (n <= 100 ? 300 : 600))); add("EDIFFG", "-0.01"); add("LWAVE", ".FALSE."); add("LCHARG", ".FALSE.");
  } else {
    const std::string ch = s.static_out == "all" ? ".TRUE." : ".FALSE.", pot = s.static_out == "all" || s.static_out == "pot" ? ".TRUE." : ".FALSE.";
    add("IBRION", "-1"); add("NSW", "0"); add("LCHARG", ch); add("LAECHG", ch); add("LORBIT", "11"); add("NEDOS", "3001"); add("LVHAR", pot); add("LWAVE", ".FALSE.");
  }
  if (dipz >= 0 && kind != "cell") { add("IDIPOL", "3"); add("DIPOL", "0.5 0.5 " + fmt("%.4f", dipz)); }
  (void)gamma;
  return L;
}

std::string incar_text(const std::string& system, const std::vector<IncarLine>& lines) {
  std::ostringstream o;
  o << "SYSTEM = " << system << "\n";
  for (const auto& l : lines) {
    std::string k = l.key;
    k.resize(std::max<size_t>(8, k.size() + 1), ' ');
    std::string v = l.value;
    v.resize(std::max<size_t>(12, v.size()), ' ');
    o << k << "= " << v << (l.why.empty() ? "" : "   # " + l.why) << "\n";
  }
  return o.str();
}

std::array<int, 3> kmesh(const System& s, double kdens, bool gamma) {
  if (gamma) return {1, 1, 1};
  return {std::max(1, int(std::ceil(kdens / norm(s.cell.a)))), std::max(1, int(std::ceil(kdens / norm(s.cell.b)))), 1};
}
std::string kpoints_text(const System& s, double kdens, bool gamma) {
  const auto k = kmesh(s, kdens, gamma);
  return "CAPS Gamma-centred mesh (kdens " + fmt("%g", kdens) + " A)\n0\nGamma\n  " + std::to_string(k[0]) + "  " + std::to_string(k[1]) + "  " + std::to_string(k[2]) + "\n  0  0  0\n";
}

// ---------------------------------------------------------------- POTCAR
std::map<std::string, std::string> potcar_map(const std::string& data_dir) {
  std::map<std::string, std::string> m;
  const Json j = read_json_file(data_dir + "/dft/potcar_map.json");
  for (const auto& [k, v] : j["map"].members()) m[k] = v.str();
  return m;
}

namespace {
PotcarInfo parse_potcar_text(const std::string& t) {
  PotcarInfo p;
  static const std::regex titel(R"(TITEL\s*=\s*\S+\s+(\S+))"), enmax(R"(ENMAX\s*=\s*([\d.]+))"), zval(R"(ZVAL\s*=\s*([\d.]+))");
  for (auto it = std::sregex_iterator(t.begin(), t.end(), titel); it != std::sregex_iterator(); ++it) p.titel.push_back((*it)[1]);
  for (auto it = std::sregex_iterator(t.begin(), t.end(), enmax); it != std::sregex_iterator(); ++it) p.enmax.push_back(std::stod((*it)[1]));
  for (auto it = std::sregex_iterator(t.begin(), t.end(), zval); it != std::sregex_iterator(); ++it) p.zval.push_back(std::stod((*it)[1]));
  double mx = 0;
  for (double e : p.enmax) mx = std::max(mx, e);
  p.encut = std::round(1.3 * mx);
  return p;
}
void check_order(const PotcarInfo& p, const std::vector<std::string>& species) {
  std::vector<std::string> got;
  for (const auto& t : p.titel) got.push_back(t.substr(0, t.find('_')));
  if (got != species) {
    std::string a, b;
    for (const auto& x : got) a += x + " ";
    for (const auto& x : species) b += x + " ";
    throw std::runtime_error("POTCAR order " + a + "does not match POSCAR " + b);
  }
}
}  // namespace

PotcarInfo assemble_potcar(const std::vector<std::string>& species, const std::string& pp_dir, const std::map<std::string, std::string>& map, const std::string& out) {
  std::string all;
  PotcarInfo info;
  for (const auto& s : species) {
    const auto it = map.find(s);
    const std::string ds = it == map.end() ? s : it->second;
    const std::string src = join(join(pp_dir, ds), "POTCAR");
    const std::string t = slurp(src);
    if (t.empty()) throw std::runtime_error("missing " + src);
    all += t;
    info.datasets.push_back(ds);
  }
  PotcarInfo p = parse_potcar_text(all);
  p.species = species;
  p.datasets = info.datasets;
  check_order(p, species);
  if (!out.empty()) spit(out, all);
  return p;
}
PotcarInfo verify_potcar(const std::string& potcar, const std::vector<std::string>& species) {
  const std::string t = slurp(potcar);
  if (t.empty()) throw std::runtime_error("cannot read " + potcar);
  PotcarInfo p = parse_potcar_text(t);
  p.species = species;
  check_order(p, species);
  return p;
}
std::string make_potcar_script(const std::map<std::string, std::string>& map) {
  std::string cases;
  for (const auto& [k, v] : map) cases += "  " + k + ") echo " + v + ";;\n";
  return "#!/bin/bash\n# Builds POTCAR from YOUR licensed PAW library. Order follows POSCAR line 6.\n# usage:  PP_DIR=/path/to/potpaw_PBE ./make_potcar.sh\nset -e\n"
         ": \"${PP_DIR:?set PP_DIR to the directory holding the PAW datasets (potpaw_PBE)}\"\npp() { case \"$1\" in\n" + cases +
         "  *) echo \"$1\";;\nesac; }\nrm -f POTCAR\nfor s in $(sed -n 6p POSCAR); do cat \"$PP_DIR/$(pp $s)/POTCAR\" >> POTCAR; done\n"
         "echo \"wrote POTCAR for: $(sed -n 6p POSCAR)\"; grep TITEL POTCAR\n";
}

// ---------------------------------------------------------------- job scripts
namespace {
std::string sbatch_head(const ClusterProfile& p, const std::string& name, int ranks, const std::string& time) {
  std::string h = "#!/bin/bash\n#SBATCH --job-name=" + name + "\n";
  if (!p.account.empty()) h += "#SBATCH --account=" + p.account + "\n";
  h += "#SBATCH --nodes=1\n#SBATCH --ntasks=" + std::to_string(ranks) + "\n#SBATCH --cpus-per-task=1\n";
  if (p.mem_per_cpu > 0) h += "#SBATCH --mem-per-cpu=" + std::to_string(p.mem_per_cpu) + "\n";
  h += "#SBATCH --time=" + time + "\n";
  if (!p.licenses.empty()) h += "#SBATCH --licenses=" + p.licenses + "\n";
  h += "#SBATCH --output=slurm-%j.out\n";
  return h;
}
std::string workspace_block(const ClusterProfile& p, const std::string& prefix) {
  std::string b = "CASE=$SLURM_SUBMIT_DIR\nWS=" + prefix + "-$SLURM_JOB_ID\n";
  if (!p.ws_allocate.empty()) b += p.ws_allocate + " > /dev/null\nWSDIR=$(" + p.ws_find + ")\n";
  else b += "WSDIR=\"\"\n";
  b += "if [ -z \"$WSDIR\" ] || [ ! -d \"$WSDIR\" ]; then echo \"workspace not found, running in $CASE\"; WSDIR=$CASE/scratch_$SLURM_JOB_ID; mkdir -p \"$WSDIR\"; fi\n"
       "echo \"workspace: $WSDIR\"\n\n"
       "## Results go to the long-term storage named in the project's .store file; home keeps only links.\n"
       "STORE=\"\"; r=$CASE\n"
       "while [ \"$r\" != \"/\" ]; do [ -f \"$r/.store\" ] && break; r=$(dirname \"$r\"); done\n"
       "if [ -f \"$r/.store\" ]; then STORE=$(head -1 \"$r/.store\")/${CASE#\"$r\"/}; echo \"results stored in: $STORE\"; fi\n"
       "link_stage() {   # $1 = stage folder: create it in $STORE and leave a link with the same name in $CASE\n"
       "  cd \"$CASE\"\n"
       "  if [ -z \"$STORE\" ]; then mkdir -p \"$1\"; return; fi\n"
       "  mkdir -p \"$STORE\"\n"
       "  if [ -d \"$1\" ] && [ ! -L \"$1\" ]; then\n"
       "    if [ -e \"$STORE/$1\" ]; then echo \"!! $1 exists in home and in $STORE: using the home copy\"; return; fi\n"
       "    mv \"$1\" \"$STORE/$1\"\n"
       "  fi\n"
       "  mkdir -p \"$STORE/$1\"\n"
       "  [ -L \"$1\" ] || ln -s \"$STORE/$1\" \"$1\"\n"
       "}\n";
  return b;
}
}  // namespace

std::string job_script(const ClusterProfile& p, const std::string& name, const std::vector<std::string>& stages, int atoms, bool gamma) {
  const ClusterRule r = parallel_rule(p, atoms, gamma);
  std::string st;
  for (const auto& s : stages) st += (st.empty() ? "" : " ") + s;
  std::string j = sbatch_head(p, name.substr(0, 40), r.ranks, r.time);
  j += "## Written by CAPS (profile " + p.name + "). Runs the stages below one after another in a scratch workspace.\n"
       "## Re-submitting the same script resumes: finished stages (file DONE) are skipped, an unfinished\n"
       "## relaxation continues from its CONTCAR (copied back every minute). If VASP writes nothing for HANG\n"
       "## seconds it is stopped and the stage is restarted from the last geometry (up to 3 attempts).\n";
  j += "STAGES=\"" + st + "\"\n\nexport OMP_NUM_THREADS=1\n";
  if (!p.modules.empty()) j += "module purge\nmodule load " + p.modules + "\n";
  j += "\nHANG=" + std::to_string(hang_seconds(p, atoms)) + "                     # seconds without any VASP output = hung -> stop and restart the stage\n";
  j += "POLL=${CAPS_POLL:-60}           # seconds between backups and watchdog checks\n";
  j += workspace_block(p, "vasp");
  j += "\nincar_for() { case $1 in 01_cell|02_cell2) echo INCAR.cell;; 03_relax) echo INCAR.relax;; *) echo INCAR.static;; esac; }\n\n"
       "prev=\"\"\n"
       "for st in $STAGES; do\n"
       "  cd \"$CASE\"\n"
       "  if [ -f \"$st/DONE\" ]; then echo \"== $st already done\"; prev=$st; continue; fi\n"
       "  link_stage \"$st\"\n"
       "  if [ -s \"$st/CONTCAR\" ] && [ \"$st\" != 04_static ]; then src=$st/CONTCAR        # resume an unfinished relaxation\n"
       "  elif [ -n \"$prev\" ]; then src=$prev/CONTCAR\n"
       "  else src=POSCAR; fi\n"
       "  cp \"$src\" \"$st/POSCAR.in\"\n"
       "  cp \"$(incar_for \"$st\")\" \"$st/INCAR\"\n"
       "  cp KPOINTS POTCAR \"$st/\"\n"
       "  mkdir -p \"${WSDIR:?}/${st:?}\"                       # fresh: the workspace is new for every job\n"
       "  cp \"$st/INCAR\" \"$st/KPOINTS\" \"$st/POTCAR\" \"$WSDIR/$st/\"\n"
       "  cp \"$st/POSCAR.in\" \"$WSDIR/$st/POSCAR\"\n"
       "  cd \"$WSDIR/$st\"\n"
       "  for attempt in 1 2 3; do\n"
       "    echo \"== $st  start $(date)$( [ $attempt -gt 1 ] && echo \"  (attempt $attempt, from last geometry)\")\"\n"
       "    ${CAPS_SRUN:-srun} " + p.vasp + " > vasp.out &\n"
       "    pid=$!\n"
       "    hung=0\n"
       "    while kill -0 $pid 2> /dev/null; do                # watchdog + live backup every minute\n"
       "      sleep $POLL\n"
       "      [ -s CONTCAR ] && cp CONTCAR OSZICAR \"$CASE/$st/\" 2> /dev/null\n"
       "      if [ $(( $(date +%s) - $(stat -c %Y vasp.out 2> /dev/null || stat -f %m vasp.out) )) -gt $HANG ]; then\n"
       "        echo \"== $st  no VASP output for $HANG s: VASP hung, stopping it $(date)\"\n"
       "        kill $pid 2> /dev/null; sleep 20; kill -9 $pid 2> /dev/null; hung=1\n"
       "      fi\n"
       "    done\n"
       "    wait $pid 2> /dev/null\n"
       "    if [ \"$st\" = 04_static ]; then grep -qs \"EDIFF is reached\" OUTCAR && break\n"
       "    else grep -qs \"reached required accuracy\" OUTCAR && break; fi\n"
       "    [ $hung = 1 ] || [ $attempt -lt 3 ] || break\n"
       "    [ \"$st\" != 04_static ] && [ -s CONTCAR ] && cp CONTCAR POSCAR    # continue from the last geometry\n"
       "    cp OSZICAR OSZICAR.attempt$attempt 2> /dev/null; cp OUTCAR OUTCAR.attempt$attempt 2> /dev/null\n"
       "  done\n"
       "  rm -f WAVECAR CHG AECCAR1 vaspout.h5                # large, never used\n"
       "  [ \"$st\" = 04_static ] || rm -f CHGCAR AECCAR*        # keep charge only from the static run\n"
       "  case $(basename \"$CASE\") in                        # charge density is needed only for complexes (Bader, CDD)\n"
       "    molecule) rm -f CHGCAR AECCAR* LOCPOT PROCAR;;\n"
       "    slab) rm -f CHGCAR AECCAR* PROCAR;;              # LOCPOT kept: work function\n"
       "  esac\n"
       "  cp -r ./* \"$CASE/$st/\"\n"
       "  cd \"$CASE/$st\"\n"
       "  if [ \"$st\" = 04_static ]; then\n"
       "    grep -q \"EDIFF is reached\" OUTCAR && touch DONE\n"
       "  else\n"
       "    grep -q \"reached required accuracy\" OUTCAR && touch DONE\n"
       "  fi\n"
       "  if [ ! -f DONE ]; then echo \"== $st NOT converged: re-submit this script to continue (sbatch job.slurm)\"; break; fi\n"
       "  echo \"== $st done  $(grep F= OSZICAR | tail -1)\"\n"
       "  prev=$st\n"
       "done\n";
  if (!p.ws_release.empty()) j += p.ws_release + "\n";
  return j;
}

std::string aimd_job_script(const ClusterProfile& p, const std::string& name, int segments, int nsw) {
  const ClusterRule r = parallel_rule(p, 1000, false);
  std::string j = sbatch_head(p, name.substr(0, 40), r.ranks, r.time);
  j += "## AIMD in segments: seg_001, seg_002, ... each starts from the previous CONTCAR (positions AND\n## velocities). Re-submit to continue. Total length = NSEG x NSW fs.\n";
  j += "NSEG=" + std::to_string(segments) + "\nNSW=" + std::to_string(nsw) + "\nexport OMP_NUM_THREADS=1\n";
  if (!p.modules.empty()) j += "module purge\nmodule load " + p.modules + "\n";
  j += workspace_block(p, "aimd");
  j += "prev=\"\"\n"
       "for i in $(seq -f \"%03g\" 1 $NSEG); do\n"
       "  cd \"$CASE\"\n"
       "  s=seg_$i\n"
       "  if [ -f \"$s/DONE\" ]; then prev=$s; continue; fi\n"
       "  link_stage \"$s\"\n"
       "  if [ -n \"$prev\" ]; then cp \"$prev/CONTCAR\" \"$s/POSCAR\"; else cp POSCAR \"$s/POSCAR\"; fi\n"
       "  cp INCAR KPOINTS POTCAR \"$s/\"\n"
       "  mkdir -p \"${WSDIR:?}/${s:?}\"\n"
       "  cp \"$s\"/* \"$WSDIR/$s/\"\n"
       "  cd \"$WSDIR/$s\"\n"
       "  echo \"== $s start $(date)\"\n"
       "  srun " + p.vasp + " > vasp.out\n"
       "  cp -r ./* \"$CASE/$s/\"\n"
       "  cd \"$CASE/$s\"\n"
       "  if [ \"$(grep -c 'T=' OSZICAR)\" -ge \"$NSW\" ]; then touch DONE; else echo \"$s incomplete: re-submit\"; break; fi\n"
       "  prev=$s\n"
       "done\n";
  if (!p.ws_release.empty()) j += p.ws_release + "\n";
  return j;
}

// ---------------------------------------------------------------- sets
const std::vector<std::string>& stages_slab() { static const std::vector<std::string> s = {"01_cell", "02_cell2", "03_relax", "04_static"}; return s; }
const std::vector<std::string>& stages_fixed() { static const std::vector<std::string> s = {"03_relax", "04_static"}; return s; }

VaspSetResult write_vasp_set(const System& s0, const std::string& dir, const VaspSetOptions& o) {
  fs::create_directories(dir);
  VaspSetResult r;
  r.dir = dir;
  const System s = s0;
  const int n = int(s.atoms.size());
  const auto stages = o.stages.empty() ? stages_fixed() : o.stages;
  const auto species = species_of(s);
  const std::string comment = formula_ordered(s) + " CAPS " + o.label;
  write_vasp_poscar(s, join(dir, "POSCAR"), comment);
  r.files.push_back("POSCAR");
  const auto map = potcar_map(o.data_dir);
  double encut = o.settings.encut;
  if (!o.pp_dir.empty()) {
    const auto p = assemble_potcar(species, o.pp_dir, map, join(dir, "POTCAR"));
    r.files.push_back("POTCAR");
    if (encut <= 0) encut = p.encut, r.notes.push_back("ENCUT " + fmt("%g", encut) + " eV = 1.3 × ENMAX " + fmt("%g", p.encut / 1.3) + " eV of the POTCARs");
  }
  if (encut <= 0) encut = 520, r.notes.push_back("ENCUT 520 eV (no POTCAR directory: 1.3 × 400 eV of C, N, O, F)");
  r.encut = encut;
  const bool dip = o.settings.dipole == 1 || (o.settings.dipole < 0 && o.asymmetric);
  double dz = -1;
  if (dip) {
    double zmin = 1e300, zmax = -1e300;
    for (const auto& a : s.atoms) { const double f = s.cell.to_fractional(a.pos)[2]; zmin = std::min(zmin, f), zmax = std::max(zmax, f); }
    dz = 0.5 * (zmin + zmax);
  }
  const auto prof = cluster_profile(o.data_dir, o.profile);
  r.parallel = parallel_rule(prof, n, o.gamma_only);
  std::vector<int> els;
  for (const auto& a : s.atoms) els.push_back(a.element);
  auto write_incar = [&](const std::string& kind, const std::string& file) {
    spit(join(dir, file), incar_text(o.label, incar_lines(kind, n, o.gamma_only, encut, dz, o.settings, r.parallel, els)));
    r.files.push_back(file);
  };
  if (std::find(stages.begin(), stages.end(), "01_cell") != stages.end()) write_incar("cell", "INCAR.cell");
  write_incar("relax", "INCAR.relax");
  write_incar("static", "INCAR.static");
  spit(join(dir, "KPOINTS"), kpoints_text(s, o.settings.kdens, o.gamma_only));
  r.kpoints = kmesh(s, o.settings.kdens, o.gamma_only);
  spit(join(dir, "job.slurm"), job_script(prof, o.label, stages, n, o.gamma_only));
  spit(join(dir, "make_potcar.sh"), make_potcar_script(map));
  make_exec(join(dir, "make_potcar.sh"));
  r.files.insert(r.files.end(), {"KPOINTS", "job.slurm", "make_potcar.sh"});
  if (!o.report_text.empty()) { spit(join(dir, "report.txt"), o.report_text); r.files.push_back("report.txt"); }
  if (o.report_json.kind() == Json::Object) { spit(join(dir, "report.json"), o.report_json.dump(2)); r.files.push_back("report.json"); }
  return r;
}

std::vector<VaspSetResult> write_adsorption_set(const AdsorbSet& set, const std::string& root, VaspSetOptions o, bool force) {
  std::vector<VaspSetResult> out;
  o.stages = stages_fixed();
  o.settings.lreal = "Auto";          // one value for every run of the set
  o.settings.ediff_relax = "1E-5";
  o.settings.dipole = 1;              // identical dipole settings everywhere: comparable energies
  auto one = [&](const System& s, const std::string& name, bool gamma, const std::string& out_level) {
    VaspSetOptions q = o;
    q.label = name;
    q.gamma_only = gamma;
    q.settings.static_out = out_level;
    out.push_back(write_vasp_set(s, join(root, name), q));
  };
  one(set.slab, "slab", false, "pot");
  one(set.molecule, "molecule", true, "none");
  for (const auto& c : set.complexes)
    if (c.status == "ok" || force) one(c.system, c.name, false, "none");
  return out;
}

// ---------------------------------------------------------------- convergence and lattice scan
std::vector<std::string> write_convergence_set(const std::string& cd, const std::vector<int>& encuts, const std::vector<int>& ks) {
  std::string base = slurp(join(cd, "INCAR.static"));
  if (base.empty()) throw std::runtime_error(cd + " has no INCAR.static");
  base = set_tag(set_tag(set_tag(base, "LAECHG", ".FALSE."), "LCHARG", ".FALSE."), "LORBIT", "0");
  const std::string job = set_stages(slurp(join(cd, "job.slurm")), "04_static");
  const std::string kp0 = slurp(join(cd, "KPOINTS"));
  std::vector<std::string> made;
  auto make = [&](const std::string& d, const std::string& incar, const std::string& kp) {
    fs::create_directories(d);
    spit(join(d, "INCAR.static"), incar);
    spit(join(d, "KPOINTS"), kp);
    spit(join(d, "job.slurm"), job);
    for (const char* f : {"POSCAR", "POTCAR"}) if (fs::exists(join(cd, f))) fs::copy_file(join(cd, f), join(d, f), fs::copy_options::overwrite_existing);
    made.push_back(d);
  };
  for (int e : encuts) make(join(cd, "conv_encut_" + std::to_string(e)), set_tag(base, "ENCUT", std::to_string(e)), kp0);
  for (int k : ks) make(join(cd, "conv_k_" + std::to_string(k)), base, "conv\n0\nGamma\n  " + std::to_string(k) + "  " + std::to_string(k) + "  1\n  0  0  0\n");
  return made;
}

std::string convergence_choice(std::vector<ConvRow>& rows, bool is_encut) {
  if (rows.empty()) return "";
  const ConvRow ref = rows.back();
  auto pass = [&](const ConvRow& r) {
    const double de = (r.e_atom - ref.e_atom) * 1000, ds = r.stress - ref.stress, df = r.fmax - ref.fmax;
    return is_encut ? std::fabs(ds) <= 1.0 && std::fabs(df) <= 0.01 : std::fabs(de) <= 1.0 && std::fabs(ds) <= 1.0;
  };
  std::string first;
  for (size_t i = 0; i < rows.size(); ++i) {
    auto& r = rows[i];
    r.de_mev = (r.e_atom - ref.e_atom) * 1000, r.ds = r.stress - ref.stress, r.df = r.fmax - ref.fmax;
    r.pass = true;
    for (size_t j = i; j < rows.size(); ++j) r.pass = r.pass && pass(rows[j]);   // this AND every more expensive one
    if (r.pass && first.empty()) first = r.value;
  }
  return first;
}

ConvResult collect_convergence(const std::string& cd) {
  ConvResult res;
  int natoms = 0;
  for (const auto& [n, c] : poscar_species(join(cd, "POSCAR"))) natoms += c;
  for (const std::string tag : {"encut", "k"}) {
    std::vector<std::pair<int, std::string>> dirs;
    std::error_code ec;
    for (const auto& e : fs::directory_iterator(cd, ec)) {
      const std::string n = e.path().filename().string(), pre = "conv_" + tag + "_";
      if (n.rfind(pre, 0) == 0) dirs.push_back({std::stoi(n.substr(pre.size())), e.path().string()});
    }
    std::sort(dirs.begin(), dirs.end());
    std::vector<ConvRow> rows;
    for (const auto& [v, d] : dirs) {
      const auto o = read_outcar(join(join(d, "04_static"), "OUTCAR"));
      if (o.e0.empty() || o.stress.empty() || o.fmax.empty()) continue;
      ConvRow r;
      r.value = std::to_string(v);
      r.e_atom = o.e0.back() / std::max(1, natoms);
      r.stress = 0.5 * (o.stress.back()[0] + o.stress.back()[1]);
      r.fmax = o.fmax.back();
      rows.push_back(r);
    }
    const std::string choice = convergence_choice(rows, tag == "encut");
    (tag == "encut" ? res.encut : res.k) = rows;
    (tag == "encut" ? res.encut_choice : res.k_choice) = choice;
  }
  return res;
}

std::vector<std::string> write_lattice_scan(const std::string& cd, double lo, double hi, double step) {
  std::string comment;
  const System at = read_vasp_poscar(join(cd, "POSCAR"), nullptr, &comment);
  const double a0 = norm(at.cell.a);
  const std::string job = set_stages(slurp(join(cd, "job.slurm")), "03_relax");
  std::vector<std::string> made;
  for (double a = lo; a <= hi + 1e-9; a += step) {
    System s = at;
    const double f = a / a0;
    for (auto& x : s.atoms) { Vec3 fr = at.cell.to_fractional(x.pos); s.cell.a = at.cell.a * f; s.cell.b = at.cell.b * f; x.pos = s.cell.to_cartesian(fr); }
    s.cell.a = at.cell.a * f, s.cell.b = at.cell.b * f;
    const std::string d = join(cd, "ascan_" + fmt("%.2f", a));
    fs::create_directories(d);
    write_vasp_poscar(s, join(d, "POSCAR"), comment);
    for (const char* f2 : {"INCAR.relax", "KPOINTS", "POTCAR"}) if (fs::exists(join(cd, f2))) fs::copy_file(join(cd, f2), join(d, f2), fs::copy_options::overwrite_existing);
    spit(join(d, "job.slurm"), job);
    made.push_back(d);
  }
  return made;
}

ScanFit fit_scan_points(const std::vector<std::pair<double, double>>& pts, const System& cell, double a0) {
  ScanFit f;
  f.points = pts;
  if (pts.size() < 3) throw std::runtime_error("need at least 3 finished points");
  // least squares: E = c0 + c1 a + c2 a²
  auto polyfit = [](const std::vector<double>& x, const std::vector<double>& y, int deg) {
    const int m = deg + 1;
    std::vector<std::vector<double>> A(size_t(m), std::vector<double>(size_t(m) + 1, 0));
    for (size_t k = 0; k < x.size(); ++k)
      for (int i = 0; i < m; ++i) {
        for (int j = 0; j < m; ++j) A[size_t(i)][size_t(j)] += std::pow(x[k], i + j);
        A[size_t(i)][size_t(m)] += y[k] * std::pow(x[k], i);
      }
    for (int i = 0; i < m; ++i) {
      int p = i;
      for (int r = i + 1; r < m; ++r) if (std::fabs(A[size_t(r)][size_t(i)]) > std::fabs(A[size_t(p)][size_t(i)])) p = r;
      std::swap(A[size_t(i)], A[size_t(p)]);
      for (int r = 0; r < m; ++r) {
        if (r == i) continue;
        const double fct = A[size_t(r)][size_t(i)] / A[size_t(i)][size_t(i)];
        for (int c = i; c <= m; ++c) A[size_t(r)][size_t(c)] -= fct * A[size_t(i)][size_t(c)];
      }
    }
    std::vector<double> c(static_cast<size_t>(m));
    for (int i = 0; i < m; ++i) c[size_t(i)] = A[size_t(i)][size_t(m)] / A[size_t(i)][size_t(i)];
    return c;
  };
  std::vector<double> a, e, v, x;
  for (const auto& [ai, ei] : pts) a.push_back(ai), e.push_back(ei);
  const auto c = polyfit(a, e, 2);
  f.a_parabola = -c[1] / (2 * c[2]);
  f.inside = f.a_parabola > *std::min_element(a.begin(), a.end()) && f.a_parabola < *std::max_element(a.begin(), a.end());
  // Birch–Murnaghan (3rd order) is linear in x = V^(−2/3): E = c0 + c1 x + c2 x² + c3 x³; V ∝ a² (c fixed)
  const double v0cell = std::fabs(dot(cell.cell.a, cross(cell.cell.b, cell.cell.c)));
  for (double ai : a) { const double vi = v0cell * (ai / a0) * (ai / a0); v.push_back(vi); x.push_back(std::pow(vi, -2.0 / 3)); }
  if (pts.size() >= 4) {
    const auto b = polyfit(x, e, 3);
    // dE/dx = 0 → c1 + 2 c2 x + 3 c3 x² = 0, the root inside the scanned range
    const double A = 3 * b[3], B = 2 * b[2], C = b[1];
    double xr = std::nan("");
    if (std::fabs(A) < 1e-30) xr = -C / B;
    else {
      const double disc = B * B - 4 * A * C;
      if (disc >= 0)
        for (double r : {(-B + std::sqrt(disc)) / (2 * A), (-B - std::sqrt(disc)) / (2 * A)})
          if (r > *std::min_element(x.begin(), x.end()) * 0.9 && r < *std::max_element(x.begin(), x.end()) * 1.1) xr = r;
    }
    if (std::isfinite(xr) && xr > 0) {
      f.v0 = std::pow(xr, -1.5);
      f.a_bm = a0 * std::sqrt(f.v0 / v0cell);
      // B0 = V d²E/dV² at V0; E(x(V)): d²E/dV² = E''(x) x'² + E'(x) x'', E'(x0) = 0
      const double xp = -2.0 / 3 * std::pow(f.v0, -5.0 / 3), e2 = 2 * b[2] + 6 * b[3] * xr;
      f.b0_gpa = f.v0 * e2 * xp * xp * 160.21766208;   // eV/Å³ → GPa
    }
  }
  return f;
}

ScanFit fit_lattice_scan(const std::string& cd) {
  std::vector<std::pair<double, double>> pts;
  std::error_code ec;
  for (const auto& e : fs::directory_iterator(cd, ec)) {
    const std::string n = e.path().filename().string();
    if (n.rfind("ascan_", 0) != 0) continue;
    const auto o = read_outcar(join(join(e.path().string(), "03_relax"), "OUTCAR"));
    if (!o.e0.empty()) pts.push_back({std::stod(n.substr(6)), o.e0.back()});
  }
  std::sort(pts.begin(), pts.end());
  const System s = read_vasp_poscar(join(cd, "POSCAR"));
  return fit_scan_points(pts, s, norm(s.cell.a));
}

// ---------------------------------------------------------------- derived sets
std::string write_charge_set(const std::string& cx) {
  const std::string src = join(join(cx, "04_static"), "POSCAR.in");
  if (!fs::exists(join(join(cx, "04_static"), "DONE"))) throw std::runtime_error(cx + "/04_static is not finished yet");
  const std::string d = join(cx, "charge");
  fs::create_directories(d);
  fs::copy_file(src, join(d, "POSCAR"), fs::copy_options::overwrite_existing);
  std::string inc = slurp(join(cx, "INCAR.static"));
  for (const char* t : {"LCHARG", "LAECHG", "LVHAR"}) inc = set_tag(inc, t, ".TRUE.");
  spit(join(d, "INCAR.static"), inc);
  for (const char* f : {"KPOINTS", "POTCAR"}) if (fs::exists(join(cx, f))) fs::copy_file(join(cx, f), join(d, f), fs::copy_options::overwrite_existing);
  std::string job = set_stages(slurp(join(cx, "job.slurm")), "04_static");
  job = std::regex_replace(job, std::regex("#SBATCH --job-name=[^\\n]*"), "#SBATCH --job-name=charge_" + fs::path(cx).filename().string().substr(0, 30));
  spit(join(d, "job.slurm"), job);
  return d;
}

std::vector<std::string> write_cdd_set(const std::string& cx, const System& slab_poscar, const std::string& pp_dir, const std::map<std::string, std::string>& map) {
  const std::string st = join(cx, "04_static");
  std::string src;
  for (const char* f : {"CONTCAR", "POSCAR.in", "POSCAR"}) if (fs::exists(join(st, f)) && fs::file_size(join(st, f)) > 0) { src = join(st, f); break; }
  if (src.empty()) throw std::runtime_error("no structure in " + st);
  const System at = read_vasp_poscar(src);
  std::map<int, int> sc, seen;
  for (const auto& a : slab_poscar.atoms) ++sc[a.element];
  System slab = at, mol = at;
  slab.atoms.clear(), mol.atoms.clear();
  for (const auto& a : at.atoms) (++seen[a.element] > sc[a.element] ? mol : slab).atoms.push_back(a);
  std::string inc = slurp(join(cx, "INCAR.static"));
  inc = set_tag(set_tag(set_tag(inc, "LCHARG", ".TRUE."), "LAECHG", ".FALSE."), "LVHAR", ".FALSE.");
  const std::string job = set_stages(slurp(join(cx, "job.slurm")), "04_static");
  std::vector<std::string> made;
  for (const auto& [name, part] : std::vector<std::pair<std::string, System>>{{"cdd_slab", slab}, {"cdd_mol", mol}}) {
    const std::string d = join(cx, name);
    fs::create_directories(d);
    write_vasp_poscar(part, join(d, "POSCAR"), formula_ordered(part) + " CAPS " + name + " (frozen complex geometry)");
    spit(join(d, "INCAR.static"), inc);
    fs::copy_file(join(cx, "KPOINTS"), join(d, "KPOINTS"), fs::copy_options::overwrite_existing);   // same cell and k → same FFT grid
    spit(join(d, "job.slurm"), job);
    if (!pp_dir.empty()) assemble_potcar(species_of(part), pp_dir, map, join(d, "POTCAR"));
    made.push_back(d);
  }
  return made;
}

std::string write_freq_set(const std::string& cd, std::vector<int> free) {
  const std::string st = join(cd, "04_static");
  std::string src;
  for (const char* f : {"CONTCAR", "POSCAR.in", "POSCAR"}) if (fs::exists(join(st, f)) && fs::file_size(join(st, f)) > 0) { src = join(st, f); break; }
  if (src.empty()) for (const char* f : {"POSCAR"}) if (fs::exists(join(cd, f))) src = join(cd, f);
  if (src.empty()) throw std::runtime_error("no structure in " + cd);
  const System at = read_vasp_poscar(src);
  if (free.empty()) {
    // the nitrile group: N (bonded to C < 1.25 Å), its C, the α-C and its H
    const auto nb = neighbours_pbc(at, 1.75);
    std::vector<std::vector<std::pair<int, double>>> adj(at.atoms.size());
    for (const auto& p : nb) adj[p.i].push_back({int(p.j), p.r});
    int N = -1;
    for (size_t i = 0; i < at.atoms.size() && N < 0; ++i)
      if (at.atoms[i].element == 7) for (const auto& [j, d] : adj[i]) if (at.atoms[size_t(j)].element == 6 && d < 1.25) { N = int(i); break; }
    if (N < 0) throw std::runtime_error("no nitrile group: give the atoms to free");
    int Cn = -1, Ca = -1;
    for (const auto& [j, d] : adj[size_t(N)]) if (at.atoms[size_t(j)].element == 6) { Cn = j; break; }
    for (const auto& [j, d] : adj[size_t(Cn)]) if (at.atoms[size_t(j)].element == 6 && j != N) { Ca = j; break; }
    free = {N, Cn, Ca};
    for (const auto& [j, d] : adj[size_t(Ca)]) if (at.atoms[size_t(j)].element == 1) free.push_back(j);
  }
  std::vector<bool> fixed(at.atoms.size(), true);
  for (int i : free) if (i >= 0 && size_t(i) < fixed.size()) fixed[size_t(i)] = false;
  const std::string d = join(cd, "freq");
  fs::create_directories(d);
  write_vasp_poscar(at, join(d, "POSCAR"), formula_ordered(at) + " CAPS frequencies (partial Hessian)", fixed);
  std::string inc = slurp(join(cd, "INCAR.static"));
  for (const char* t : {"IBRION", "NSW", "EDIFF", "LCHARG", "LAECHG", "LVHAR", "LORBIT", "NEDOS"}) inc = remove_tag(inc, t);
  inc += "## vibrational analysis (finite differences) of the freed atoms; compare SHIFTS, not absolute PBE values\n"
         "IBRION  = 5\nNFREE   = 2\nPOTIM   = 0.015\nNSW     = 1\nISYM    = 0\nEDIFF   = 1E-7\nLCHARG  = .FALSE.\nLAECHG  = .FALSE.\nLVHAR   = .FALSE.\nLORBIT  = 0\n";
  spit(join(d, "INCAR.static"), inc);
  for (const char* f : {"KPOINTS", "POTCAR"}) if (fs::exists(join(cd, f))) fs::copy_file(join(cd, f), join(d, f), fs::copy_options::overwrite_existing);
  spit(join(d, "job.slurm"), set_stages(slurp(join(cd, "job.slurm")), "04_static"));
  return d;
}

std::string write_aimd_set(const System& at, const std::string& out, const VaspSetOptions& o, double temp, int nsw, int seg) {
  fs::create_directories(out);
  write_vasp_poscar(at, join(out, "POSCAR"), formula_ordered(at) + " CAPS AIMD start");
  spit(join(out, "KPOINTS"), kpoints_text(at, o.settings.kdens));
  const auto k = kmesh(at, o.settings.kdens);
  double zmin = 1e300, zmax = -1e300;
  for (const auto& a : at.atoms) { const double f = at.cell.to_fractional(a.pos)[2]; zmin = std::min(zmin, f), zmax = std::max(zmax, f); }
  std::string pomass;
  for (const auto& s : species_of(at)) pomass += (pomass.empty() ? "" : " ") + fmt("%.3f", s == "H" ? 2.014 : element(element_from_symbol(s)).mass);
  const double encut = o.settings.encut > 0 ? o.settings.encut : 520;
  std::ostringstream in;
  in << "SYSTEM = AIMD " << formula_ordered(at) << "\n"
     << "## electronic: same functional, cutoff and dispersion as the adsorption set\n"
     << "ENCUT   = " << fmt("%g", encut) << "\nPREC    = Normal       # sufficient for MD forces; MD energies are not mixed with E_bind\n"
     << "ALGO    = Fast\nEDIFF   = 1E-5\nNELM    = 100\nISMEAR  = 0\nSIGMA   = 0.05\nGGA     = PE\nIVDW    = " << o.settings.ivdw
     << "\nLASPH   = .TRUE.\nLMAXMIX = 4\nISPIN   = 1\nLREAL   = Auto\nISYM    = 0            # no symmetry in MD\n"
     << "NCORE   = " << parallel_rule(cluster_profile(o.data_dir, o.profile), int(at.atoms.size()), false).ncore << "\nKPAR    = " << (k[0] * k[1] > 1 ? 2 : 1) << "\n"
     << "## molecular dynamics, NVT\nIBRION  = 0\nMDALGO  = 5            # CSVR (Bussi) thermostat\nCSVR_PERIOD = 20       # thermostat time constant in steps\n"
     << "TEBEG   = " << fmt("%g", temp) << "\nTEEND   = " << fmt("%g", temp) << "\nISIF    = 2\nPOTIM   = 1.0          # fs; safe because H carries the deuterium mass below\n"
     << "POMASS  = " << pomass << "\nNSW     = " << nsw << "        # one segment; job.slurm chains segments\nNBLOCK  = 1\nLWAVE   = .FALSE.\nLCHARG  = .FALSE.\n"
     << "IDIPOL  = 3            # energy-only dipole correction (no LDIPOL: stable on low-work-function surfaces)\nDIPOL   = 0.5 0.5 " << fmt("%.4f", 0.5 * (zmin + zmax)) << "\n";
  spit(join(out, "INCAR"), in.str());
  spit(join(out, "job.slurm"), aimd_job_script(cluster_profile(o.data_dir, o.profile), "aimd", seg, nsw));
  if (!o.pp_dir.empty()) assemble_potcar(species_of(at), o.pp_dir, potcar_map(o.data_dir), join(out, "POTCAR"));
  return out;
}

// ---------------------------------------------------------------- job tools
std::vector<std::string> submittable(const std::vector<std::string>& dirs, const std::vector<std::string>& queued) {
  std::vector<std::string> out;
  for (const auto& d : dirs) {
    if (!fs::is_directory(d) || !fs::exists(join(d, "job.slurm"))) continue;
    std::error_code ec;
    const std::string real = fs::weakly_canonical(d, ec).string();
    bool q = false;
    for (const auto& x : queued) q |= fs::weakly_canonical(x, ec).string() == real;
    if (q) continue;
    // the last stage of the job finished?
    std::string stages = "04_static";
    static const std::regex st("STAGES=\"([^\"]*)\"");
    std::smatch m;
    const std::string job = slurp(join(d, "job.slurm"));
    if (std::regex_search(job, m, st)) { std::stringstream ss(m[1]); for (std::string w; ss >> w;) stages = w; }
    if (job.find("NSEG=") != std::string::npos) {
      // AIMD: every segment done?
      static const std::regex ns("NSEG=(\\d+)");
      if (std::regex_search(job, m, ns)) { char b[16]; std::snprintf(b, sizeof b, "seg_%03d", std::stoi(m[1])); stages = b; }
    }
    if (fs::exists(join(join(d, stages), "DONE"))) continue;
    out.push_back(d);
  }
  return out;
}

int update_job_scripts(const std::vector<std::string>& dirs, const ClusterProfile& p) {
  int n = 0;
  for (const auto& d : dirs) {
    const std::string path = join(d, "job.slurm");
    const std::string old = slurp(path);
    if (old.empty() || old.find("NSEG=") != std::string::npos) continue;
    auto get = [&](const std::string& pat, const std::string& def) { std::smatch m; return std::regex_search(old, m, std::regex(pat)) ? m[1].str() : def; };
    const std::string name = get("#SBATCH --job-name=(\\S+)", fs::path(d).filename().string());
    const int ranks = std::stoi(get("#SBATCH --ntasks=(\\d+)", "104"));
    const std::string time = get("#SBATCH --time=(\\S+)", "48:00:00");
    const std::string stages = get("STAGES=\"([^\"]*)\"", "03_relax 04_static");
    const std::string hang = get("HANG=(\\d+)", "1800");
    std::vector<std::string> st;
    std::stringstream ss(stages);
    for (std::string w; ss >> w;) st.push_back(w);
    std::string fresh = job_script(p, name, st, 0, false);
    fresh = std::regex_replace(fresh, std::regex("#SBATCH --ntasks=\\d+"), "#SBATCH --ntasks=" + std::to_string(ranks));
    fresh = std::regex_replace(fresh, std::regex("#SBATCH --time=\\S+"), "#SBATCH --time=" + time);
    fresh = std::regex_replace(fresh, std::regex("HANG=\\d+"), "HANG=" + hang);
    if (fresh != old) { spit(path, fresh); ++n; }
  }
  return n;
}

std::string reset_stage(const std::string& cd, const std::string& stage) {
  const std::string p = join(cd, stage);
  std::error_code ec;
  if (!fs::exists(fs::symlink_status(p, ec))) throw std::runtime_error(p + " does not exist");
  const std::string target = fs::canonical(p).string();
  std::string nw = target + "_bad";
  for (int k = 2; fs::exists(nw); ++k) nw = target + "_bad" + std::to_string(k);
  fs::rename(target, nw);
  if (fs::is_symlink(p)) fs::remove(p);
  return nw;
}

std::vector<ToolAction> cleanup_plan(const std::string& root, bool keep_best) {
  std::vector<ToolAction> out;
  std::set<std::string> seen;
  auto add = [&](const std::string& p, const std::string& why) {
    std::error_code ec;
    if (!fs::exists(p, ec) || !seen.insert(p).second) return;
    long long b = 0;
    if (fs::is_directory(p)) { for (const auto& e : fs::recursive_directory_iterator(p, ec)) if (e.is_regular_file()) b += (long long)e.file_size(); }
    else b = (long long)fs::file_size(p, ec);
    out.push_back({p, why, b});
  };
  std::error_code ec;
  const std::vector<std::string> always = {"vaspout.h5", "AECCAR1", "WAVECAR", "CHG"}, charge = {"CHGCAR", "AECCAR0", "AECCAR2", "PROCAR"};
  for (auto it = fs::recursive_directory_iterator(root, fs::directory_options::follow_directory_symlink, ec); it != fs::recursive_directory_iterator(); ++it) {
    if (it->path().filename() != "DONE") continue;
    const std::string st = it->path().parent_path().string();
    const std::string case_name = it->path().parent_path().parent_path().filename().string();
    for (const auto& f : always) add(join(st, f), "not used");
    if (it->path().parent_path().filename() == "04_static") {
      if (case_name == "molecule") { for (const auto& f : charge) add(join(st, f), "molecule: only E0 and DOSCAR are used"); add(join(st, "LOCPOT"), "molecule: only E0 and DOSCAR are used"); }
      else if (case_name == "slab" || case_name.rfind("complex_", 0) != 0) for (const auto& f : charge) add(join(st, f), "slab: Bader not needed");
    }
  }
  if (keep_best) {
    for (auto it = fs::recursive_directory_iterator(root, ec); it != fs::recursive_directory_iterator(); ++it) {
      if (it->path().filename() != "binding_energies.csv") continue;
      const std::string set = it->path().parent_path().string();
      std::string best;
      double eb = 1e300;
      std::stringstream ss(slurp(it->path().string()));
      std::string l;
      std::getline(ss, l);
      while (std::getline(ss, l)) {
        std::stringstream ls(l);
        std::string c, ec2, e;
        std::getline(ls, c, ',');
        std::getline(ls, ec2, ',');
        std::getline(ls, e, ',');
        try { if (std::stod(e) < eb) eb = std::stod(e), best = c; } catch (...) {}
      }
      for (const auto& e : fs::directory_iterator(set, ec)) {
        const std::string n = e.path().filename().string();
        if (n.rfind("complex_", 0) == 0 && n != best && fs::exists(join(join(e.path().string(), "04_static"), "DONE")))
          for (const auto& f : {"CHGCAR", "AECCAR0", "AECCAR2", "PROCAR", "LOCPOT"}) add(join(join(e.path().string(), "04_static"), f), "not the best complex");
      }
    }
  }
  return out;
}

std::vector<ToolAction> store_plan(const std::string& root, const std::string& store) {
  std::vector<ToolAction> out;
  static const std::regex stage("^(0[1-4]_[a-z0-9]+|seg_\\d+)$");
  std::error_code ec;
  for (auto it = fs::recursive_directory_iterator(root, ec); it != fs::recursive_directory_iterator(); ++it) {
    if (!it->is_directory() || fs::is_symlink(it->path())) continue;
    if (!std::regex_match(it->path().filename().string(), stage)) continue;
    it.disable_recursion_pending();
    long long b = 0;
    for (const auto& e : fs::recursive_directory_iterator(it->path(), ec)) if (e.is_regular_file()) b += (long long)e.file_size();
    out.push_back({it->path().string(), join(store, fs::relative(it->path(), root).string()), b});
  }
  return out;
}

void apply_plan(const std::vector<ToolAction>& plan, const std::string& kind, const std::string& root, const std::string& store) {
  for (const auto& a : plan) {
    if (kind == "cleanup") {
      if (fs::is_symlink(a.path)) { const auto t = fs::read_symlink(a.path); fs::remove(a.path); fs::remove_all(t); }
      else fs::remove_all(a.path);
    } else if (kind == "store") {
      if (fs::exists(a.what)) continue;   // already in storage
      fs::create_directories(fs::path(a.what).parent_path());
      std::error_code ec;
      fs::rename(a.path, a.what, ec);
      if (ec) {   // another filesystem: copy, then remove
        fs::copy(a.path, a.what, fs::copy_options::recursive | fs::copy_options::copy_symlinks);
        fs::remove_all(a.path);
      }
      fs::create_directory_symlink(a.what, a.path);
    }
  }
  if (kind == "store" && !root.empty()) { fs::create_directories(store); spit(join(root, ".store"), store + "\n"); }
}

}  // namespace caps
