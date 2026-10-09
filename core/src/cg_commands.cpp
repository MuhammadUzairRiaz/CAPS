// The coarse-graining workflow as commands (caps/cg_commands.hpp).
#include "caps/cg_commands.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <set>
#include <sstream>
#include <stdexcept>

#include "caps/bundle.hpp"
#include "caps/cg_bonded.hpp"
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
      {"cgfit", "caps cgfit bonded MAP [FRAMES …] [MAP [FRAMES …] …] --types FILE -o DIR  |  caps cgfit refine MAP CG_FRAMES … --tables DIR/bonded.json -o DIR2",
       "Bonded coarse-grained potentials by tabulated Boltzmann inversion over many frames of many systems mapped alike: each MAP "
       "(STEM.map.json from cgmap) followed by its frames (CG dumps STEM.cg.lammpstrj or the single frame STEM.cg.data); each "
       "system's histograms normalised, then weighted. U(r) = −kT ln[P/r²], U(θ) = −kT ln[P/sin θ], U(φ) = −kT ln P, inverted only "
       "where P exceeds a share of its maximum, with slope-continuous walls outside. refine: one bonded IBI step from a CG run "
       "with the tables (U ← U + α kT ln(P_CG/P_target)), for the shift the non-bonded 1–3 and 1–4 terms cause.",
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
        {"o", "bonded", "the output folder (tables, bonded.in, gromacs/, bonded.json)"}},
       {"caps cgfit bonded cg/PBS.map.json cg/PBS.cg.lammpstrj cg/PBSA.map.json cg/PBSA.cg.lammpstrj cg/PBAT.map.json cg/PBAT.cg.lammpstrj --types cg/types.json -T 300 -o bonded",
        "caps cgfit refine cg/PBS.map.json run1/cg.lammpstrj --tables bonded/bonded.json --types cg/types.json -o bonded_it2"}},
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
    if (mode != "bonded" && mode != "refine") throw std::invalid_argument("cgfit " + mode + ": bonded or refine (the non-bonded fits follow)");
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
