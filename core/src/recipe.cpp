#include "caps/recipe.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <cctype>
#include <set>
#include <sstream>

#include "caps/analysis.hpp"
#include "caps/elements.hpp"
#include "caps/solvate.hpp"
#include "caps/nano.hpp"
#include "caps/lattice.hpp"
#include "caps/spacegroup.hpp"
#include "caps/crystal.hpp"
#include "caps/config.hpp"
#include "caps/dynamics.hpp"
#include "caps/cbmc.hpp"
#include "caps/dlpoly.hpp"
#include "caps/equilibrate.hpp"
#include "caps/polystats.hpp"
#include "caps/ffdef.hpp"
#include "caps/field.hpp"
#include "caps/grow.hpp"
#include "caps/io.hpp"
#include "caps/mechanics.hpp"
#include "caps/molecule.hpp"
#include "caps/pack.hpp"
#include "caps/polymer.hpp"
#include "caps/react.hpp"
#include "caps/relax.hpp"
#include "caps/typing.hpp"
#include "caps/uff.hpp"

namespace caps {

namespace {

const std::vector<std::string> kStages = {"build", "type", "grow", "react", "relax", "cbmc", "md", "equilibrate", "analyze", "export"};
const std::set<std::string> kTop = {"recipe", "name", "build", "type", "grow", "react", "relax", "cbmc", "md", "equilibrate", "analyze", "export", "electrostatics", "cutoff", "seed", "threads"};

std::string g6(double x) { char b[32]; std::snprintf(b, sizeof b, "%.6g", x); return b; }
double num(const Json& j, const char* k, double def) { return j.is_object() && j.has(k) && j[k].is_number() ? j[k].number() : def; }
std::string text(const Json& j, const char* k, const std::string& def) { return j.is_object() && j.has(k) && j[k].is_string() ? j[k].str() : def; }
bool flag(const Json& j, const char* k, bool def) {
  if (!j.is_object() || !j.has(k)) return def;
  return j[k].kind() == Json::Bool ? j[k].boolean() : j[k].is_number() ? j[k].number() != 0 : def;
}

KeyValues approx(const EnergyOptions& e, int threads) {
  const int t = threads > 0 ? threads : (max_threads() > 0 ? max_threads() : 1);
  return {{"van der Waals", "cut-off " + g6(e.cutoff) + " Å" + (e.tail ? " + tail correction" : "")},
          {"Electrostatics", !e.coulomb ? "off" : e.electrostatics == EnergyOptions::Electrostatics::PME ? "SPME · relative tolerance " + g6(e.ewald_rtol)
                                                                                                          : "damped shifted force · cut-off " + g6(e.cutoff) + " Å"},
          {"Constraints", "none"},
          {"Precision", "double · reproducible with " + std::to_string(t) + " threads"}};
}
void elec_cite(std::vector<std::string>& c, const EnergyOptions& e) {
  if (e.coulomb) c.push_back(e.electrostatics == EnergyOptions::Electrostatics::PME ? "essmann1995" : "fennell2006");
}
ProvStep step(const std::string& engine, const std::string& summary, KeyValues params, const std::string& rng, std::vector<std::string> cites, KeyValues ap = {}) {
  ProvStep s;
  s.engine = engine;
  s.summary = summary;
  s.params = std::move(params);
  s.rng = rng;
  s.cites = std::move(cites);
  s.approximations = std::move(ap);
  s.time = now_iso();
  return s;
}
std::string seeded(uint64_t s) { return "mt19937-64 · seed " + std::to_string(s); }

Trajectory as_trajectory(const System& s) {
  Trajectory t;
  t.topology = s;
  std::vector<Vec3> p;
  for (const auto& a : s.atoms) p.push_back(a.pos);
  t.positions.push_back(std::move(p));
  t.cells.push_back(s.cell);
  t.timesteps.push_back(0);
  return t;
}

}  // namespace

std::vector<std::string> recipe_stages(const Json& r) {
  if (!r.is_object()) throw RecipeError(2, "a recipe is a mapping of stages");
  for (const auto& [k, v] : r.members())
    if (!kTop.count(k)) throw RecipeError(2, "unknown key '" + k + "' (stages: build, type, grow, react, relax, md, equilibrate, analyze, export)");
  if (r.has("recipe") && r["recipe"].is_number() && r["recipe"].number() != 1) throw RecipeError(2, "this CAPS reads recipe version 1");
  if (!r.has("build")) throw RecipeError(2, "a recipe needs a build stage (polymer, molecule or file)");
  std::vector<std::string> out;
  for (const auto& s : kStages) if (r.has(s)) out.push_back(s);
  if (r.has("grow") && !(r["build"].is_object() && r["build"].has("polymer"))) throw RecipeError(2, "grow needs build: { polymer: … }");
  return out;
}

// react.templates: built-in names or template text (one string or a list)
static std::vector<ReactionTemplate> react_templates_of(const Json& J) {
  std::vector<std::string> items;
  if (J.has("templates") && J["templates"].is_array()) for (const auto& t : J["templates"].items()) items.push_back(t.str());
  else if (J.has("templates") && J["templates"].is_string()) items.push_back(J["templates"].str());
  if (items.empty()) throw RecipeError(2, "react needs templates (" + [] { std::string n; for (const auto& x : builtin_template_names()) n += (n.empty() ? "" : ", ") + x; return n; }() + " or template text)");
  std::string text;
  for (const auto& it : items) {
    const auto names = builtin_template_names();
    if (std::find(names.begin(), names.end(), it) != names.end()) text += builtin_template(it) + "\n";
    else if (it.find('\n') != std::string::npos) text += it + "\n";
    else throw RecipeError(2, "react: no built-in template '" + it + "'");
  }
  auto t = parse_templates(text);
  if (t.empty()) throw RecipeError(2, "react: the templates hold no reaction");
  return t;
}

RecipeCheck check_recipe(const Json& r) {
  RecipeCheck c;
  std::vector<std::string> stages;
  try { stages = recipe_stages(r); } catch (const RecipeError& e) { c.code = e.code; c.error = e.what(); return c; }
  auto list = [](const Json& a) { std::string out; if (a.is_array()) for (const auto& x : a.items()) out += (out.empty() ? "" : ", ") + (x.is_string() ? x.str() : g6(x.number())); return out; };
  for (const auto& st : stages) {
    const Json& J = r[st];
    RecipeStageInfo info{st, "", true};
    try {
      if (st == "build") {
        if (J.has("polymer")) {
          const Json& P = J["polymer"];
          std::vector<std::string> units;
          if (P.has("units") && P["units"].is_array()) for (const auto& u : P["units"].items()) units.push_back(u.is_string() ? u.str() : u.text("smiles"));
          else units.push_back(text(P, "smiles", ""));
          std::string formula;
          for (const auto& u : units) {
            if (u.empty()) throw RecipeError(2, "a repeat unit needs SMILES with two * points");
            try { formula += (formula.empty() ? "" : " + ") + repeat_unit_info(u).formula; } catch (const std::exception& e) { throw RecipeError(2, "repeat unit " + u + ": " + e.what()); }
          }
          info.summary = formula + " · DP " + g6(num(P, "dp", 20)) + " × " + g6(num(P, "chains", 10)) + " chains · " + text(P, "tacticity", "atactic") +
                         (units.size() > 1 ? " · " + text(P, "sequence", "homopolymer") : "");
        } else if (J.has("pack")) {
          const Json& P = J["pack"];
          info.summary = P.is_string() ? "pack " + P.str() : "pack " + g6(double(P.has("molecules") && P["molecules"].is_array() ? P["molecules"].size() : 0)) + " kinds of molecule";
        } else if (J.has("molecule")) {
          info.summary = "molecule " + (J["molecule"].is_string() ? J["molecule"].str() : text(J["molecule"], "smiles", ""));
        } else if (J.has("file")) {
          info.summary = "file " + J["file"].str();
        } else if (J.has("crystal")) {
          info.summary = "crystal " + (J["crystal"].has("cif") ? J["crystal"]["cif"].str() : text(J["crystal"], "group", "P 1"));
        } else if (J.has("surface")) {
          info.summary = "surface of " + text(J["surface"], "cif", "?");
        } else if (J.has("nano")) {
          info.summary = "nano " + text(J["nano"], "kind", "tube");
        } else if (J.has("solvate")) {
          info.summary = "solvate in " + text(J["solvate"], "solvent", "water");
        } else {
          throw RecipeError(2, "build needs polymer, molecule, pack, file, crystal, surface, nano or solvate");
        }
      } else if (st == "type") {
        info.summary = text(J, "forcefield", "default") + " · charges " + text(J, "charges", "auto");
      } else if (st == "grow") {
        info.summary = (J.has("box") ? "box " + g6(num(J, "box", 0)) + " Å" : "ρ " + g6(num(J, "density", 0.5)) + " g/cm³") + " · seed " + g6(num(J, "seed", 1)) +
                       " · best of " + g6(num(J, "trials", 120)) + " trials";
      } else if (st == "react") {
        const auto T = react_templates_of(J);
        std::string names;
        for (const auto& t : T) names += (names.empty() ? "" : ", ") + t.name;
        info.summary = names + " · target " + g6(num(J, "target", 1.0)) + (J.has("insert") ? " · insert " + g6(num(J["insert"], "count", 10)) + " × " + text(J["insert"], "smiles", "") : "");
      } else if (st == "relax") {
        const Minimiser m = minimiser_from_string(text(J, "method", "lbfgs"));
        info.summary = std::string(to_string(m)) + " · |F|max " + g6(num(J, "fmax", 0.5));
      } else if (st == "cbmc") {
        info.summary = g6(num(J, "moves", 1000)) + " end regrowths · k " + g6(num(J, "trials", 8)) + " · " + g6(num(J, "temperature", 300)) + " K";
      } else if (st == "md") {
        const std::string ens = text(J, "ensemble", "nvt");
        if (ens != "nve" && ens != "nvt" && ens != "npt" && ens != "nph") throw RecipeError(2, "md.ensemble: nve, nvt, npt or nph");
        info.summary = ens + " · " + g6(num(J, "temperature", 300)) + " K · " + g6(num(J, "ps", 10)) + " ps";
      } else if (st == "equilibrate") {
        ProtocolParams pp;
        pp.t_max = num(J, "t_max", pp.t_max);
        pp.t_final = num(J, "t_final", pp.t_final);
        if (J.has("p_max")) pp.p_max = J["p_max"].number() / 1.01325;
        pp.time_scale = num(J, "time_scale", 1.0);
        c.protocol = text(J, "protocol", "larsen21");
        try { c.schedule = protocol_by_name(c.protocol, pp); } catch (const std::exception& e) { throw RecipeError(2, e.what()); }
        double ps = 0;
        for (const auto& sg : c.schedule) ps += sg.ps;
        info.summary = c.protocol + " · " + std::to_string(c.schedule.size()) + " steps · " + g6(ps) + " ps";
      } else if (st == "analyze") {
        info.summary = J.has("properties") ? list(J["properties"]) : "density";
      } else if (st == "export") {
        std::vector<std::string> f;
        if (J.is_array()) for (const auto& x : J.items()) f.push_back(x.str());
        else if (J.is_string()) f.push_back(J.str());
        for (const auto& x : f)
          if (x != "lammps" && x != "moltemplate" && x != "dlpoly" && x != "gromacs" && x != "gro" && x != "pdb" && x != "xyz" && x != "mol2" && x != "sdf" && x != "cif")
            throw RecipeError(2, "export: unknown format '" + x + "'");
        info.summary = list(J.is_array() ? J : Json::array());
      }
    } catch (const RecipeError& e) {
      info.ok = false;
      info.summary = e.what();
      if (!c.code) { c.code = e.code; c.error = st + ": " + e.what(); }
    } catch (const std::exception& e) {
      info.ok = false;
      info.summary = e.what();
      if (!c.code) { c.code = 2; c.error = st + ": " + e.what(); }
    }
    c.stages.push_back(std::move(info));
  }
  return c;
}

RecipeResult run_recipe(const Json& r, const RecipeOptions& o) {
  const auto stages = recipe_stages(r);
  RecipeResult res;
  res.name = text(r, "name", "recipe");
  const int n = int(stages.size());
  const long long seed0 = o.seed >= 0 ? o.seed : (long long)num(r, "seed", 1);
  auto seed_of = [&](const Json& j) { return uint64_t(o.seed >= 0 ? o.seed : (long long)num(j, "seed", double(seed0))); };
  EnergyOptions energy;
  energy.cutoff = num(r, "cutoff", 10.0);
  if (text(r, "electrostatics", "dsf") == "pme") energy.electrostatics = EnergyOptions::Electrostatics::PME;
  energy.threads = o.threads > 0 ? o.threads : int(num(r, "threads", 0));
  const std::string by = "CAPS 0.1.0";
  auto report = [&](int k, const std::string& name, const std::string& detail, const std::string& status, double f) {
    if (o.progress) o.progress(RecipeEvent{k + 1, n, name, detail, status, f});
  };
  auto path_of = [&](const std::string& p) { return std::filesystem::path(p).is_absolute() ? p : (std::filesystem::path(o.base_dir) / p).string(); };

  System s;
  ChainSpec spec;
  bool polymer = false;
  int chains = 1;
  std::shared_ptr<const ForceField> ff;
  std::string filled_from;   // type.fill_from: the force fields gaps were filled from
  int filled_terms = 0;
  std::string ffname = "built-in default (GAFF for C and H, UFF otherwise)";
  if (!o.sha256.empty())
    res.manifest.steps.push_back(step("recipe.run", "recipe " + res.name, {{"recipe", res.name}, {"sha256", o.sha256}, {"stages", std::to_string(n)}}, "", {}));
  // Types a structure with the recipe's force field (throws RecipeError 3 for untyped atoms or missing parameters).
  std::string borrowed;   // a library force field typed with its family's rules
  auto type_now = [&](System& sys) {
    const Json T = r.has("type") ? r["type"] : Json::object();
    std::string name = text(T, "forcefield", "default");
    // library ids before the force fields got CAPS's own names ("opls2005-dlfield" is now "opls2005")
    if (name.size() > 8 && name.compare(name.size() - 8, 8, "-dlfield") == 0) name.resize(name.size() - 8);
    std::string charges = text(T, "charges", "auto");
    std::string companion_note;   // charges "increments": whose charges they are
    ff.reset();
    try {
      if (name == "default") {
        const bool ch = std::all_of(sys.atoms.begin(), sys.atoms.end(), [](const Atom& a) { return a.element == 1 || a.element == 6; });
        ff = std::make_shared<ForceField>(ch ? assign_gaff(sys) : assign_uff(sys));
        ffname = ch ? "the built-in GAFF (C and H)" : "UFF";
      } else if (is_uff(name)) {
        UffOptions uo;
        uo.qeq = charges == "qeq";
        uo.keep_charges = charges == "keep";
        ff = std::make_shared<ForceField>(assign_uff(sys, uo));
        ffname = "UFF";
      } else {
        std::string path = path_of(name), typing = T.has("typing") ? path_of(T["typing"].str()) : "";
        std::vector<std::pair<std::string, std::string>> family;
        if (!std::filesystem::exists(path) && !o.forcefield_dir.empty()) {   // a library id (or its first word): the catalogue's file and typing rules
          std::ifstream cf(std::filesystem::path(o.forcefield_dir) / "catalogue.json");
          if (cf) {
            const Json cat = Json::parse(std::string((std::istreambuf_iterator<char>(cf)), std::istreambuf_iterator<char>()));
            const Json* hit = nullptr;
            for (const auto& e : cat["forcefields"].items())
              if (e.text("id") == name) { hit = &e; break; }
            if (!hit)
              for (const auto& e : cat["forcefields"].items())
                if (e.text("id").rfind(name + "-", 0) == 0) { hit = &e; break; }
            auto rules = [](const Json& e) { return e.has("typing") && e["typing"].is_object() && e["typing"].has("rules") ? e["typing"]["rules"].str() : std::string(); };
            if (hit && !(hit->has("file") && (*hit)["file"].is_string()))   // an alias or a template: no parameter file of its own
              throw RecipeError(2, "'" + name + "' has no parameter file in the library: " + hit->text("notes", hit->text("status", "")));
            if (hit && hit->has("file")) {
              path = (std::filesystem::path(o.forcefield_dir) / (*hit)["file"].str()).string();
              std::string rel = rules(*hit);
              if (rel.empty()) {   // no rules of its own: the newest rules of its family (gaff2 → GAFF's; the type names are shared)
                std::string fam = hit->text("id");
                fam = fam.substr(0, fam.find('-'));
                while (!fam.empty() && std::isdigit((unsigned char)fam.back())) fam.pop_back();
                for (const auto& e : cat["forcefields"].items())
                  if (!rules(e).empty() && e.text("id").rfind(fam, 0) == 0) family.insert(family.begin(), {e.text("id"), (std::filesystem::path(o.forcefield_dir).parent_path() / rules(e)).string()});
              }
              if (typing.empty() && !rel.empty()) typing = (std::filesystem::path(o.forcefield_dir).parent_path() / rel).string();
            }
          }
        }
        if (!std::filesystem::exists(path)) throw RecipeError(2, "no force field '" + name + "' (a library id, uff, default or a path)");
        FFDef def = load_forcefield(path);
        if (!typing.empty()) load_typing(def, typing);
        // fill_from: another force field's bond, angle and dihedral rules, used only where this one has none (as Field's
        // "Fill gaps from…": OPLS-AA 2024 has no CM–CT–CT–CM torsion, OPLS 2005 does); what was borrowed is reported
        filled_from.clear();
        if (T.has("fill_from")) {
          std::vector<std::string> donors;
          if (T["fill_from"].is_array()) for (const auto& x : T["fill_from"].items()) donors.push_back(x.str());
          else donors.push_back(T["fill_from"].str());
          FFDef fill;
          for (const auto& dn : donors) {
            std::string dp = path_of(dn);
            if (!std::filesystem::exists(dp) && !o.forcefield_dir.empty()) {
              std::ifstream cf(std::filesystem::path(o.forcefield_dir) / "catalogue.json");
              if (cf) {
                const Json cat = Json::parse(std::string((std::istreambuf_iterator<char>(cf)), std::istreambuf_iterator<char>()));
                for (const auto& e : cat["forcefields"].items())
                  if (e.text("id") == dn && e.has("file") && e["file"].is_string()) dp = (std::filesystem::path(o.forcefield_dir) / e["file"].str()).string();
              }
            }
            if (!std::filesystem::exists(dp)) throw RecipeError(2, "type.fill_from: no force field '" + dn + "'");
            const FFDef donor = load_forcefield(dp);
            GapFill g = gap_fill_rules(def, donor, donor.name);
            if (!g.kept) throw RecipeError(2, "type.fill_from: none of " + donor.name + "'s rules can apply to " + def.name + " (no shared atom classes)");
            fill.bonds.insert(fill.bonds.end(), g.rules.bonds.begin(), g.rules.bonds.end());
            fill.angles.insert(fill.angles.end(), g.rules.angles.begin(), g.rules.angles.end());
            fill.dihedrals.insert(fill.dihedrals.end(), g.rules.dihedrals.begin(), g.rules.dihedrals.end());
            filled_from += (filled_from.empty() ? "" : ", ") + donor.name;
          }
          prepend_fill(def, fill);
        }
        for (const auto& [id, rules_path] : family) {   // newest first; the first whose types all exist in this force field
          if (!typing.empty()) break;
          try {
            FFDef trial = def;
            load_typing(trial, rules_path);
            def = std::move(trial);
            borrowed = id;
            break;
          } catch (const std::exception&) {}
        }
        // a united-atom force field: hydrogens on carbon fold into their carbons (charges computed first, then summed)
        const std::string ua = prepare_for_forcefield(sys, def, charges);
        std::vector<std::string> types;
        if (!def.typing.empty()) {
          const TypingResult tr = assign_types(sys, def);
          if (tr.untyped) throw RecipeError(3, untyped_message(def, sys, tr.untyped));
          types = tr.types;
        } else {
          for (const auto& a : sys.atoms) types.push_back(a.name);
        }
        ParamReport rep;
        if (charges == "increments") {   // another force field's bond-increment charges with these types
          const auto q = companion_charges(sys, def, path, &companion_note);
          for (size_t i = 0; i < sys.atoms.size(); ++i) sys.atoms[i].charge = q[i];
          sys.has_charges = true;
          charges = "keep";
        }
        const bool auto_charges = charges == "auto";
        if (charges == "auto") {   // the force field's charges, the file's, or (GAFF-like fields have none per type) Gasteiger, else QEq
          charges = sys.has_charges && !polymer ? "keep" : "types";
          try {
            ff = std::make_shared<ForceField>(parameterize(sys, def, types, charges, &rep, false));
          } catch (const FFError& e) {
            if (std::string(e.what()).find("has no charge for type") == std::string::npos) throw;
            charges = "gasteiger";
            rep = ParamReport{};
          }
        }
        if (!ff && charges == "gasteiger" && auto_charges) {
          try {
            ff = std::make_shared<ForceField>(parameterize(sys, def, types, charges, &rep, false));
          } catch (const std::exception&) {   // Gasteiger–Marsili has no parameters here (S=O, Si, metals): QEq covers every element
            charges = "qeq";
            rep = ParamReport{};
          }
        }
        if (!ff) ff = std::make_shared<ForceField>(parameterize(sys, def, types, charges, &rep, false));
        // automatic charges: fixed per-type charges that miss the formal charge give way to the companion force
        // field's bond-increment charges when it names one (OPLS-AA 2024 → OPLS 2005)
        if (auto_charges && charges == "types" && !def.charge_increments_from.empty()) {
          double net = 0;
          for (double q : ff->charge) net += q;
          int formal = 0;
          for (int fc : perceive(sys).charge) formal += fc;
          if (std::fabs(net - formal) > 1e-3) {
            try {
              std::string note;
              const auto q = companion_charges(sys, def, path, &note);
              for (size_t i = 0; i < sys.atoms.size(); ++i) sys.atoms[i].charge = q[i];
              sys.has_charges = true;
              ParamReport rep2;
              ff = std::make_shared<ForceField>(parameterize(sys, def, types, "keep", &rep2, false));
              rep = std::move(rep2);
              companion_note = note;
              charges = "keep";
            } catch (const FFError&) {}
          }
        }
        if (!rep.missing.empty()) throw RecipeError(3, std::to_string(rep.missing.size()) + " parameters missing in " + def.name + " (first: " + rep.missing.front() + ")");
        ffname = def.name + (ua.empty() ? "" : " (united-atom: hydrogens on carbon folded into their carbons)");
        filled_terms = 0;
        for (const auto& [rule, count] : rep.used)
          if (rule.find("filled:") != std::string::npos) filled_terms += count;
        if (!filled_from.empty()) ffname += " (gaps filled from " + filled_from + ": " + std::to_string(filled_terms) + " terms)";
      }
    } catch (const RecipeError&) { throw; } catch (const std::exception& e) { throw RecipeError(3, e.what()); }
    const std::string ch = !companion_note.empty() ? companion_note.substr(9) : charges == "qeq" ? "QEq" : charges == "gasteiger" ? "Gasteiger" : charges == "keep" ? "the file's" : charges == "auto" && ffname == "UFF" ? "no" : "from the force field";
    std::vector<std::string> c;
    if (ffname == "UFF") c.push_back("rappe1992");
    else if (ffname.find("GAFF") != std::string::npos) c.push_back("wang2004");
    if (charges == "qeq") c.push_back("rappe1991");
    if (charges == "gasteiger") c.push_back("gasteiger1980");
    KeyValues kv = {{"force field", ffname}, {"charges", ch}};
    if (!borrowed.empty()) kv.push_back({"typing rules", borrowed});
    if (!filled_from.empty()) kv.push_back({"gaps filled from", filled_from + " (" + std::to_string(filled_terms) + " terms where " + ffname.substr(0, ffname.find(" (gaps")) + " has none)"});
    bool found = false;
    for (auto& ps : res.manifest.steps)
      if (ps.engine == "field.assign") { ps = step("field.assign", ffname + " · charges " + ch, kv, "", c); found = true; }
    if (!found) res.manifest.steps.push_back(step("field.assign", ffname + " · charges " + ch, kv, "", c));
  };
  for (int k = 0; k < n; ++k) {
    const std::string& st = stages[size_t(k)];
    const Json& J = r[st];
    try {
      if (st == "build") {
        report(k, st, "", "running", 0);
        if (J.has("polymer")) {
          const Json& P = J["polymer"];
          polymer = true;
          if (P.has("units") && P["units"].is_array()) {
            char name = 'A';
            for (const auto& u : P["units"].items()) spec.units.push_back({std::string(1, name++), u.is_string() ? u.str() : u.text("smiles")});
          } else {
            spec.units.push_back({"A", text(P, "smiles", "")});
          }
          if (spec.units.empty() || spec.units[0].smiles.empty()) throw RecipeError(2, "build.polymer needs smiles (with two * points) or units");
          UnitInfo info;
          for (const auto& u : spec.units) {
            try { info = repeat_unit_info(u.smiles); } catch (const std::exception& e) { throw RecipeError(2, "repeat unit " + u.smiles + ": " + e.what()); }
          }
          spec.dp = int(num(P, "dp", 20));
          chains = int(num(P, "chains", 10));
          const std::string tac = text(P, "tacticity", "atactic");
          spec.tacticity = tac == "isotactic" ? Tacticity::Isotactic : tac == "syndiotactic" ? Tacticity::Syndiotactic : Tacticity::Atactic;
          const std::string seq = text(P, "sequence", "homopolymer");
          try { spec.sequence = sequence_from_string(seq); } catch (const std::exception& e) { throw RecipeError(2, std::string("build.polymer.sequence: ") + e.what()); }
          if (P.has("blocks") && P["blocks"].is_array()) for (const auto& b : P["blocks"].items()) spec.blocks.push_back(int(b.number()));
          if (P.has("weights") && P["weights"].is_array()) for (const auto& w : P["weights"].items()) spec.weights.push_back(w.number());
          spec.pattern = text(P, "pattern", "");
          if (P.has("pm")) spec.pm = P["pm"].number();
          if (P.has("p_mr") && P.has("p_rm")) spec.p_mr = P["p_mr"].number(), spec.p_rm = P["p_rm"].number();   // Markov tacticity
          if (P.has("r1")) spec.r1 = P["r1"].number();                                                        // terminal model
          if (P.has("r2")) spec.r2 = P["r2"].number();
          // end groups (a preset or a SMILES with *) and the linkage of the units
          spec.head_cap = text(P, "head_cap", "");
          spec.tail_cap = text(P, "tail_cap", "");
          try { spec.linkage = linkage_from_string(text(P, "linkage", "head-to-tail")); } catch (const std::exception& e) { throw RecipeError(2, std::string("build.polymer.linkage: ") + e.what()); }
          spec.inversion = num(P, "inversion", spec.inversion);
          // per-chain lengths: given (chain_dp: [..]) or drawn (lengths: {distribution, nn, pdi, seed})
          if (P.has("chain_dp") && P["chain_dp"].is_array()) {
            for (const auto& x : P["chain_dp"].items()) spec.chain_dp.push_back(int(x.number()));
          } else if (P.has("lengths") && P["lengths"].is_object()) {
            const Json& L = P["lengths"];
            try {
              if (L.has("histogram") && L["histogram"].is_array()) {   // [[length, weight], …]
                std::vector<std::pair<int, double>> h;
                for (const auto& b : L["histogram"].items()) if (b.is_array() && b.size() == 2) h.push_back({int(b[0].number()), b[1].number()});
                spec.chain_dp = draw_chain_lengths(h, chains, uint64_t(num(L, "seed", 1)));
              } else {
                spec.chain_dp = draw_chain_lengths(text(L, "distribution", "schulz-zimm"), num(L, "nn", spec.dp), num(L, "pdi", 1.1), chains, uint64_t(num(L, "seed", 1)));
              }
            } catch (const std::exception& e) { throw RecipeError(2, std::string("build.polymer.lengths: ") + e.what()); }
          }
          // architecture: linear (default), star {arms}, comb {arm_dp, spacing}, branched {arm_dp, branch_probability},
          // dendrimer {arms, arm_dp, generations}
          try { spec.architecture = architecture_from_string(text(P, "architecture", "linear")); } catch (const std::exception& e) { throw RecipeError(2, std::string("build.polymer.architecture: ") + e.what()); }
          if (P.has("arms")) spec.arms = int(P["arms"].number());
          if (P.has("arm_dp")) spec.arm_dp = int(P["arm_dp"].number());
          if (P.has("spacing")) spec.spacing = int(P["spacing"].number());
          if (P.has("branch_probability")) spec.branch_probability = P["branch_probability"].number();
          if (P.has("generations")) spec.generations = int(P["generations"].number());
          report(k, st, "DP " + std::to_string(spec.dp) + " × " + std::to_string(chains) + (spec.architecture == Architecture::Linear ? " chains" : std::string(" ") + to_string(spec.architecture) + " molecules") +
                            " · unit " + info.formula + " · " + tac, "done", 1);
        } else if (J.has("pack")) {
          // a packmol input file, or molecules (SMILES or files) with counts in a box:
          // pack: {box: [x, y, z] | edge, tolerance, seed, density (compress to), molecules: [{smiles | file, count}]}
          const Json& P = J["pack"];
          PackOptions po;
          std::vector<PackItem> items;
          std::string label;
          try {
            if (P.is_string()) {
              const std::string p = path_of(P.str());
              if (!std::filesystem::exists(p)) throw RecipeError(2, "pack: no file " + p);
              items = read_packmol_input(p, po, nullptr);
              label = std::filesystem::path(p).filename().string();
            } else {
              po.tolerance = num(P, "tolerance", 2.0);
              po.seed = uint64_t(num(P, "seed", 1));
              double L[3] = {30, 30, 30};
              if (P.has("box") && P["box"].is_array() && P["box"].size() == 3) for (size_t q = 0; q < 3; ++q) L[q] = P["box"][q].number();
              else if (P.has("box")) L[0] = L[1] = L[2] = P["box"].number();
              po.cell.a = {L[0], 0, 0}, po.cell.b = {0, L[1], 0}, po.cell.c = {0, 0, L[2]};
              if (P.has("density")) po.compress_to = P["density"].number();
              if (!P.has("molecules") || !P["molecules"].is_array() || P["molecules"].size() == 0) throw RecipeError(2, "pack needs molecules: [{smiles | file, count}]");
              for (const auto& m : P["molecules"].items()) {
                PackItem it;
                it.count = int(num(m, "count", 1));
                if (m.has("smiles")) {
                  BuildOptions bo;
                  bo.forcefield = "uff";
                  bo.seed = po.seed;
                  it.molecule = build_molecule(m["smiles"].str(), bo).system;
                  it.name = m["smiles"].str();
                } else if (m.has("file")) {
                  const std::string p = path_of(m["file"].str());
                  if (!std::filesystem::exists(p)) throw RecipeError(2, "pack: no file " + p);
                  it.molecule = open_file(p, "").frame(0);
                  it.name = std::filesystem::path(p).filename().string();
                } else {
                  throw RecipeError(2, "pack: each molecule needs smiles or file");
                }
                Region box;
                box.kind = Region::InsideBox;
                box.a = {0, 0, 0}, box.b = {L[0], L[1], L[2]};
                it.regions.push_back(box);
                label += (label.empty() ? "" : " + ") + std::to_string(it.count) + " × " + it.name;
                items.push_back(std::move(it));
              }
            }
            PackReport pr;
            s = pack(items, po, &pr);
            std::vector<std::string> c = {"martinez2009"};
            res.manifest.steps.push_back(step("pack.optimise", "molecules packed into a box: " + label,
                                              {{"tolerance", g6(po.tolerance) + " Å"}, {"molecules", label}}, seeded(po.seed), c));
            report(k, st, label + " · " + std::to_string(s.atoms.size()) + " atoms", "done", 1);
          } catch (const RecipeError&) { throw; }
          catch (const std::exception& e) { throw RecipeError(4, std::string("pack: ") + e.what()); }
        } else if (J.has("molecule")) {
          BuildOptions bo;
          bo.forcefield = "uff";
          bo.seed = seed_of(J["molecule"]);
          const std::string smi = J["molecule"].is_string() ? J["molecule"].str() : text(J["molecule"], "smiles", "");
          BuildResult br;
          try { br = build_molecule(smi, bo); } catch (const std::exception& e) { throw RecipeError(2, "molecule: " + std::string(e.what())); }
          s = br.system;
          res.manifest.steps.push_back(step("chem.build", "3D structure from SMILES · " + br.method, {{"smiles", br.graph.smiles}, {"clean-up", "UFF"}}, seeded(bo.seed), {"rappe1992"}));
          report(k, st, br.graph.smiles + " · " + std::to_string(s.atoms.size()) + " atoms", "done", 1);
        } else if (J.has("file")) {
          const std::string p = path_of(J["file"].str());
          if (!std::filesystem::exists(p)) throw RecipeError(2, "no file " + p);
          s = open_file(p, J.has("topology") ? path_of(J["topology"].str()) : "").frame(0);
          res.manifest.inputs.push_back({std::filesystem::path(p).filename().string(), ""});
          res.manifest.steps.push_back(step("io.read", "read " + std::filesystem::path(p).filename().string(), {{"file", std::filesystem::path(p).filename().string()}, {"atoms", std::to_string(s.atoms.size())}}, "", {}));
          report(k, st, std::filesystem::path(p).filename().string() + " · " + std::to_string(s.atoms.size()) + " atoms", "done", 1);
        } else if (J.has("crystal")) {
          // crystal: {cif} or {group, cell: [a, b, c, α, β, γ], sites: "Ti1 Ti 0 0 0; O1 O 0.3 0.3 0"}; supercell [na, nb, nc];
          // cell: primitive | niggli | conventional
          const Json& C = J["crystal"];
          std::string what;
          try {
            if (C.has("cif")) {
              const std::string p = path_of(C["cif"].str());
              s = read_cif(p);
              res.manifest.inputs.push_back({std::filesystem::path(p).filename().string(), ""});
              what = std::filesystem::path(p).filename().string();
            } else {
              CrystalSpec cs;
              cs.space_group = text(C, "group", "P 1");
              if (!C.has("cell") || !C["cell"].is_array() || C["cell"].size() < 3) throw std::invalid_argument("crystal needs cif, or group, cell [a, b, c(, α, β, γ)] and sites");
              const auto& cl = C["cell"];
              cs.a = cl[0].number(), cs.b = cl[1].number(), cs.c = cl[2].number();
              if (cl.size() >= 6) cs.alpha = cl[3].number(), cs.beta = cl[4].number(), cs.gamma = cl[5].number();
              std::string sites = text(C, "sites", "");
              std::replace(sites.begin(), sites.end(), ';', '\n');
              std::istringstream in(sites);
              for (std::string line; std::getline(in, line);) {
                std::istringstream ls(line);
                std::string lab, el;
                double x, y, z;
                if (ls >> lab >> el >> x >> y >> z) cs.sites.push_back({lab, element_from_symbol(el), {x, y, z}});
              }
              if (cs.sites.empty()) throw std::invalid_argument("crystal sites: 'LABEL EL x y z; …'");
              s = build_crystal(cs);
              what = cs.space_group;
            }
            const std::string cell = text(C, "cell_reduce", text(C, "reduce", ""));
            if (cell == "primitive") s = find_primitive_cell(s);
            else if (cell == "niggli") s = niggli_cell(s);
            else if (cell == "conventional") s = conventional_cell(s);
            else if (!cell.empty()) throw std::invalid_argument("reduce: primitive, niggli or conventional");
            if (C.has("supercell") && C["supercell"].is_array() && C["supercell"].size() == 3)
              s = supercell(s, int(C["supercell"][0].number()), int(C["supercell"][1].number()), int(C["supercell"][2].number()));
          } catch (const RecipeError&) { throw; } catch (const std::exception& e) { throw RecipeError(2, std::string("build.crystal: ") + e.what()); }
          res.manifest.steps.push_back(step("build.crystal", what, {{"atoms", std::to_string(s.atoms.size())}}, "", {}));
          report(k, st, what + " · " + std::to_string(s.atoms.size()) + " atoms", "done", 1);
        } else if (J.has("surface")) {
          // surface: {cif, hkl: [h, k, l], layers, termination, vacuum, supercell: [na, nb], passivate}
          const Json& C = J["surface"];
          SlabOptions so;
          if (C.has("hkl") && C["hkl"].is_array() && C["hkl"].size() == 3) so.h = int(C["hkl"][0].number()), so.k = int(C["hkl"][1].number()), so.l = int(C["hkl"][2].number());
          so.layers = int(num(C, "layers", so.layers));
          so.termination = int(num(C, "termination", so.termination));
          so.vacuum = num(C, "vacuum", so.vacuum);
          so.passivate = flag(C, "passivate", false);
          if (C.has("supercell") && C["supercell"].is_array() && C["supercell"].size() == 2) so.na = int(C["supercell"][0].number()), so.nb = int(C["supercell"][1].number());
          SlabReport sr;
          try {
            if (!C.has("cif")) throw std::invalid_argument("surface needs cif");
            s = cleave(read_cif(path_of(C["cif"].str())), so, &sr);
          } catch (const std::exception& e) { throw RecipeError(2, std::string("build.surface: ") + e.what()); }
          char b[160];
          std::snprintf(b, sizeof b, "(%d %d %d) · %d layers · %.1f Å slab · %zu atoms", so.h, so.k, so.l, so.layers, sr.thickness, s.atoms.size());
          res.manifest.steps.push_back(step("build.surface", b, {{"vacuum", g6(so.vacuum) + " Å"}}, "", {}));
          report(k, st, b, "done", 1);
        } else if (J.has("nano")) {
          // nano: {kind: tube | sheet | particle, …}: tube {n, m, length, material, walls, periodic}; sheet {lx, ly,
          // layers, material, periodic}; particle {cif, shape, radius, height, top_ratio, passivate}
          const Json& C = J["nano"];
          const std::string kind = text(C, "kind", "tube");
          try {
            if (kind == "tube") {
              NanotubeOptions o;
              o.material = text(C, "material", o.material);
              o.n = int(num(C, "n", o.n)), o.m = int(num(C, "m", o.m));
              o.length = num(C, "length", o.length);
              o.walls = int(num(C, "walls", o.walls));
              o.periodic = flag(C, "periodic", o.periodic);
              s = nanotube(o);
            } else if (kind == "sheet") {
              SheetOptions o;
              o.material = text(C, "material", o.material);
              o.lx = num(C, "lx", o.lx), o.ly = num(C, "ly", o.ly);
              o.layers = int(num(C, "layers", o.layers));
              o.periodic = flag(C, "periodic", o.periodic);
              s = graphene_sheet(o);
            } else if (kind == "particle") {
              ParticleOptions o;
              if (!C.has("cif")) throw std::invalid_argument("a particle needs cif");
              o.shape = particle_shape_from_string(text(C, "shape", "sphere"));
              o.radius = num(C, "radius", o.radius);
              o.height = num(C, "height", o.height);
              o.top_ratio = num(C, "top_ratio", o.top_ratio);
              o.passivate = flag(C, "passivate", false);
              s = nanoparticle(read_cif(path_of(C["cif"].str())), o);
            } else throw std::invalid_argument("kind: tube, sheet or particle");
          } catch (const std::exception& e) { throw RecipeError(2, "build.nano: " + std::string(e.what())); }
          res.manifest.steps.push_back(step("build.nano", kind, {{"atoms", std::to_string(s.atoms.size())}}, "", {}));
          report(k, st, kind + " · " + std::to_string(s.atoms.size()) + " atoms", "done", 1);
        } else if (J.has("solvate")) {
          // solvate: {solute: FILE (optional), solvent, water_model, padding | edge, salt, concentration, ions:
          // none | neutralise | concentration}
          const Json& C = J["solvate"];
          SolvateOptions o;
          o.solvent = text(C, "solvent", o.solvent);
          o.water_model = text(C, "water_model", o.water_model);
          o.tolerance = num(C, "tolerance", o.tolerance);
          o.salt = text(C, "salt", o.salt);
          o.concentration = num(C, "concentration", o.concentration);
          const std::string ions = text(C, "ions", "concentration");
          o.ion_mode = ions == "none" ? 0 : ions == "neutralise" || ions == "neutralize" ? 1 : 2;
          o.seed = uint64_t(seed_of(C));
          System solute;
          const bool has = C.has("solute");
          if (C.has("edge")) o.shape = 0, o.edge = num(C, "edge", o.edge);
          else o.shape = has ? 2 : 0, o.padding = num(C, "padding", o.padding);
          try {
            if (has) solute = open_file(path_of(C["solute"].str())).frame(0);
            SolvateReport sr;
            s = solvate(has ? &solute : nullptr, o, &sr);
            report(k, st, std::to_string(sr.plan.solvent) + " × " + sr.plan.solvent_name + (sr.plan.cations + sr.plan.anions ? " · " + std::to_string(sr.plan.cations) + " + " + std::to_string(sr.plan.anions) + " ions" : "") +
                              " · " + std::to_string(s.atoms.size()) + " atoms", "done", 1);
            res.manifest.steps.push_back(step("build.solvate", sr.plan.solvent_name, {{"molecules", std::to_string(sr.plan.solvent)}}, seeded(o.seed), {"martinez2009"}));
          } catch (const std::exception& e) { throw RecipeError(2, "build.solvate: " + std::string(e.what())); }
        } else {
          throw RecipeError(2, "build needs polymer, molecule, pack, file, crystal, surface, nano or solvate");
        }
      } else if (st == "type") {
        report(k, st, "", "running", 0);
        System probe = s;
        if (polymer) {   // the structure comes later (grow): type one chain of the spec now, so a gap in the force field stops here
          GrowOptions g;
          g.chains = 1;
          g.seed = uint64_t(seed0);
          g.density = 0.05;
          g.auto_scale = true;   // quaternary backbones (methacrylates, polyisobutylene) need a lower contact scale
          probe = grow_chains(spec, g, nullptr);
        }
        type_now(probe);
        report(k, st, ffname + (borrowed.empty() ? "" : " (rules of " + borrowed + ")") + " · " + std::to_string(probe.atoms.size()) + (polymer && r.has("grow") ? " typed (one chain) · " : " typed · ") + "0 missing parameters", "done", 1);
        if (polymer) ff.reset();   // typed again once the cell is grown
      } else if (st == "grow") {
        GrowOptions g;
        g.chains = chains;
        g.seed = seed_of(J);
        g.density = num(J, "density", 0.5);
        // contact_scale: a number (1.0 = full contact limits), or auto (step down only where a chain cannot be placed)
        const bool numeric = J.has("contact_scale") && J["contact_scale"].is_number();
        g.contact_scale = numeric ? J["contact_scale"].number() : 1.0;
        g.auto_scale = !numeric && text(J, "contact_scale", "auto") == "auto";
        if (J.has("box")) g.box = J["box"].number();
        if (J.has("region") && J["region"].is_object()) {   // {shape: slab, thickness, vacuum} | {shape: cylinder | around_cylinder, radius, length}
          const Json& R = J["region"];
          const std::string shape = text(R, "shape", "cubic");
          if (shape == "slab") g.slab_thickness = num(R, "thickness", 30), g.slab_vacuum = num(R, "vacuum", 30);
          else if (shape == "cylinder" || shape == "around_cylinder")
            g.cylinder_radius = num(R, "radius", 10), g.cylinder_length = num(R, "length", 0), g.cylinder_outside = shape == "around_cylinder";
          else if (shape != "cubic") throw RecipeError(2, "grow.region.shape: cubic, slab, cylinder or around_cylinder");
        }
        if (J.has("orientation") && J["orientation"].is_object()) {   // {axis: x | y | z | [x, y, z], strength: s (kT)}
          const Json& O = J["orientation"];
          if (O.has("axis") && O["axis"].is_array() && O["axis"].size() == 3)
            for (size_t k = 0; k < 3; ++k) g.orient_axis[k] = O["axis"][k].number();
          else {
            const std::string ax = text(O, "axis", "z");
            if (ax != "x" && ax != "y" && ax != "z") throw RecipeError(2, "grow.orientation.axis: x, y, z or [x, y, z]");
            g.orient_axis = {ax == "x" ? 1.0 : 0.0, ax == "y" ? 1.0 : 0.0, ax == "z" ? 1.0 : 0.0};
          }
          g.orient_strength = num(O, "strength", 4);
          if (g.orient_strength < 0) throw RecipeError(2, "grow.orientation.strength must be ≥ 0");
        }
        g.curve = flag(J, "curve", true);
        g.lookahead = std::clamp(int(num(J, "lookahead", 1)), 1, 4);
        const std::string method = text(J, "method", "trials");
        if (method != "trials" && method != "rosenbluth" && method != "rosenbluth_lj")
          throw RecipeError(2, "grow.method: trials (best of k by contact margin), rosenbluth (soft spheres) or rosenbluth_lj (UFF Lennard-Jones); "
                               "configurational-bias Monte Carlo with acceptance is not built");
        g.method = method == "rosenbluth" ? 1 : method == "rosenbluth_lj" ? 2 : 0;
        g.method_temperature = num(J, "temperature", 450);
        if (J.has("trials")) g.trials = int(J["trials"].number());
        g.progress = [&](int done, int total, int restarts) {
          report(k, st, "chain " + std::to_string(done) + "/" + std::to_string(total) + (restarts ? " · " + std::to_string(restarts) + " restarts" : ""), "running", total ? double(done) / total : 0);
          return true;
        };
        GrowReport gr;
        try { s = grow_chains(spec, g, &gr); } catch (const std::exception& e) { throw RecipeError(4, std::string("grow: ") + e.what()); }
        KeyValues gp = {{"chains", std::to_string(chains)}, {"DP", std::to_string(spec.dp)}, {"trials", std::to_string(g.trials)}, {"density", g6(g.density) + " g/cm³"}};
        if (!spec.chain_dp.empty()) {   // the sample actually built, never the target as if achieved
          double s1 = 0, s2 = 0;
          for (int n : spec.chain_dp) s1 += n, s2 += double(n) * n;
          const double nn = s1 / double(spec.chain_dp.size());
          gp[1] = {"DP", "per chain, sample Nn " + g6(nn) + ", Đ " + g6(s2 / s1 / nn) + ", " + std::to_string(*std::min_element(spec.chain_dp.begin(), spec.chain_dp.end())) + "–" +
                         std::to_string(*std::max_element(spec.chain_dp.begin(), spec.chain_dp.end()))};
        }
        if (spec.sequence == Sequence::Terminal) gp.push_back({"sequence", "terminal model, r1 " + g6(spec.r1) + ", r2 " + g6(spec.r2)});
        if (spec.p_mr >= 0 && spec.p_rm >= 0) gp.push_back({"tacticity", "Markov, P(r|m) " + g6(spec.p_mr) + ", P(m|r) " + g6(spec.p_rm)});
        if (spec.architecture == Architecture::Star) gp.push_back({"architecture", "star, " + std::to_string(spec.arms) + " arms of DP units on one core carbon"});
        else if (spec.architecture == Architecture::Comb)
          gp.push_back({"architecture", "comb, side chains of " + std::to_string(spec.arm_dp) + " units on every " + std::to_string(spec.spacing) + "th unit (or the one beside it with more room)"});
        else if (spec.architecture == Architecture::Branched)
          gp.push_back({"architecture", "branched, side chains of " + std::to_string(spec.arm_dp) + " units with probability " + g6(spec.branch_probability) + " per backbone unit"});
        else if (spec.architecture == Architecture::Dendrimer)
          gp.push_back({"architecture", "dendrimer, generation " + std::to_string(spec.generations) + ": " + std::to_string(spec.arms) +
                                            " core arms of DP units, each end splitting in two branches of " + std::to_string(spec.arm_dp) + " units"});
        if (!gr.notes.empty()) gp.push_back({"built", gr.notes.front()});
        if (g.method > 0) gp.push_back({"method", (g.method == 1 ? "Rosenbluth, soft spheres" : "Rosenbluth, UFF Lennard-Jones") + std::string(" at ") + g6(g.method_temperature) +
                                               " K; ln W per chain " + g6(gr.ln_rosenbluth)});
        res.manifest.steps.push_back(g.method == 0 ? step("grow.trials", std::to_string(chains) + " chains grown in a periodic cell, best-of-k trial placement by contact margin",
                                                          std::move(gp), seeded(g.seed), {"parsons2005", "matsumoto1998"})
                                                   : step("grow.rosenbluth", std::to_string(chains) + " chains grown in a periodic cell, trials drawn by Rosenbluth weight",
                                                          std::move(gp), seeded(g.seed),
                                                          g.method == 1 ? std::vector<std::string>{"rosenbluth1955", "theodorou1985", "jorgensen1984"}
                                                                        : std::vector<std::string>{"rosenbluth1955", "siepmann1992", "rappe1992", "jorgensen1984"}));
        report(k, st, "best of " + std::to_string(g.trials) + " trials · " + std::to_string(chains) + " chains · " + std::to_string(s.atoms.size()) + " atoms · box " + g6(gr.box) + " Å", "done", 1);
      } else if (st == "react") {
        // react: {templates: [sulfur_allylic] | text, insert: {smiles: SS, count: 12, tolerance: 2}, target, cycles, per_cycle,
        //         capture, relax, md_ps, temperature, during_md, seed} — curatives packed into the cell, then cure cycles; the
        // network is typed again with the recipe's force field (a term it lacks stops the run here)
        report(k, st, "", "running", 0);
        if (J.has("insert")) {
          const Json& I = J["insert"];
          BuildOptions bo;
          bo.forcefield = "uff";
          const BuildResult br = build_molecule(text(I, "smiles", ""), bo);
          PackOptions po;
          po.tolerance = num(I, "tolerance", po.tolerance);
          po.seed = seed_of(I);
          PackReport pr;
          const int n = std::max(1, int(num(I, "count", 10)));
          try { s = insert_molecules(s, br.system, n, po, &pr); } catch (const std::exception& e) { throw RecipeError(4, std::string("react.insert: ") + e.what()); }
          res.manifest.steps.push_back(step("pack.insert", std::to_string(n) + " × " + text(I, "smiles", "") + " inserted into the free space",
                                            {{"smiles", text(I, "smiles", "")}, {"count", std::to_string(n)}, {"tolerance", g6(po.tolerance) + " Å"}}, seeded(po.seed),
                                            {"martinez2009", "rappe1992"}));
        }
        ReactOptions ro;
        ro.templates = react_templates_of(J);
        if (J.has("capture")) for (auto& t : ro.templates) t.capture = J["capture"].number();
        ro.seed = seed_of(J);
        ro.max_cycles = int(num(J, "cycles", 50));
        ro.max_per_cycle = int(num(J, "per_cycle", 5));
        ro.target_conversion = num(J, "target", 1.0);
        ro.relax = flag(J, "relax", true);
        ro.md_ps = num(J, "md_ps", 0);
        ro.temperature = num(J, "temperature", 500);
        ro.during_md = flag(J, "during_md", false);
        ro.energy = energy;
        ro.progress = [&](const CycleRow& c) {
          report(k, st, "cycle " + std::to_string(c.cycle) + " · " + std::to_string(c.total) + " reactions · conversion " + g6(c.conversion), "running",
                 std::min(1.0, c.conversion / std::max(1e-9, ro.target_conversion)));
          return true;
        };
        if (!s.unwrapped && s.cell.valid()) make_molecules_whole(s);
        ReactReport rr;
        try { react(s, ro, &rr); } catch (const std::exception& e) { throw RecipeError(4, std::string("react: ") + e.what()); }
        std::string names;
        for (const auto& t : ro.templates) names += (names.empty() ? "" : ", ") + t.name;
        const double conv = rr.cycles.empty() ? 0 : rr.cycles.back().conversion;
        res.manifest.steps.push_back(step("react.templates", std::to_string(rr.reactions) + " reactions · conversion " + g6(conv),
                                          {{"templates", names}, {"target conversion", g6(ro.target_conversion)}, {"cycles", std::to_string(rr.cycles.size())},
                                           {"relax between cycles", ro.relax ? "yes" : "no"}, {"MD between cycles", g6(ro.md_ps) + " ps"}},
                                          seeded(ro.seed), {"matsumoto1998"}));
        ff.reset();
        if (r.has("type")) type_now(s);   // the network in the recipe's force field
        report(k, st, std::to_string(rr.reactions) + " reactions · conversion " + g6(conv) + (rr.gel_conversion >= 0 ? " · gel at " + g6(rr.gel_conversion) : "") +
                          (r.has("type") ? " · typed again: " + ffname : ""), "done", 1);
      } else if (st == "cbmc") {
        // configurational-bias regrowth of chain ends with the recipe's force field (Siepmann & Frenkel)
        if (!ff) type_now(s);
        CbmcOptions co;
        co.moves = int(num(J, "moves", 1000));
        co.trials = int(num(J, "trials", 8));
        co.max_torsions = int(num(J, "max_torsions", 4));
        co.temperature = num(J, "temperature", 300);
        co.cutoff = std::min(energy.cutoff, num(J, "cutoff", 9.0));
        co.coulomb = energy.coulomb;
        co.seed = seed_of(J);
        co.progress = [&](int done, int acc) {
          report(k, st, std::to_string(done) + " of " + std::to_string(co.moves) + " · " + std::to_string(acc) + " accepted", "running", co.moves ? double(done) / co.moves : 1);
          return true;
        };
        CbmcReport cr;
        try { cbmc_regrow(s, *ff, co, &cr); } catch (const std::exception& e) { throw RecipeError(4, std::string("cbmc: ") + e.what()); }
        res.manifest.steps.push_back(step("cbmc.regrow", "Configurational-bias Monte Carlo regrowth of chain ends",
                                          {{"moves", std::to_string(cr.attempted) + " end regrowths, " + std::to_string(cr.accepted) + " accepted (" + g6(std::round(1000 * cr.acceptance) / 10) + " %)"},
                                           {"trials", std::to_string(co.trials) + " torsions per bond, up to " + std::to_string(co.max_torsions) + " bonds per move"},
                                           {"temperature", g6(co.temperature) + " K"},
                                           {"trial energies", std::string("van der Waals") + (co.coulomb ? ", DSF electrostatics" : "") + ", torsions · cut-off " + g6(cr.cutoff) + " Å"},
                                           {"⟨R²⟩ end-to-end", g6(cr.r2_before) + " → " + g6(cr.r2_after) + " Å²"},
                                           {"force field", ffname}},
                                          seeded(co.seed), {"siepmann1992", "rosenbluth1955"}));
        report(k, st, std::to_string(cr.accepted) + " of " + std::to_string(cr.attempted) + " accepted · ⟨R²⟩ " + g6(cr.r2_before) + " → " + g6(cr.r2_after) + " Å²", "done", 1);
      } else if (st == "relax" || st == "md" || st == "equilibrate") {
        if (!ff) type_now(s);
        if (st == "relax") {
          RelaxOptions ro;
          ro.field = ff;
          ro.energy = energy;
          ro.method = minimiser_from_string(text(J, "method", "lbfgs"));
          ro.ftol = num(J, "fmax", 0.5);
          ro.max_iterations = int(num(J, "max_iterations", 5000));
          ro.pushoff = flag(J, "pushoff", true);
          ro.pushoff_ramp_ps = num(J, "pushoff_md_ps", 0);   // push-off by MD with the force cap ramped (Auhl et al.)
          ro.pushoff_cap = num(J, "pushoff_cap", 0);
          ro.pushoff_temperature = num(J, "pushoff_temperature", 300);
          ro.target_density = num(J, "target_density", 0);
          ro.progress = [&](const RelaxProgress& p) {
            report(k, st, p.stage + " · |F|max " + g6(p.fmax), "running", p.stages ? double(p.stage_index) / p.stages : 0);
            return true;
          };
          RelaxReport rr;
          try { relax(s, ro, &rr); } catch (const std::exception& e) { throw RecipeError(4, std::string("relax: ") + e.what()); }
          std::vector<std::string> c = {ro.method == Minimiser::LBFGS ? "liu1989" : ro.method == Minimiser::FIRE ? "bitzek2006" : "polak1969"};
          if (ro.pushoff) c.push_back("auhl2003");
          elec_cite(c, energy);
          res.manifest.steps.push_back(step(std::string("relax.") + (ro.method == Minimiser::LBFGS ? "lbfgs" : ro.method == Minimiser::FIRE ? "fire" : ro.method == Minimiser::ConjugateGradient ? "cg" : "sd"),
                                            std::string(to_string(ro.method)) + (rr.converged ? ", converged" : ", stopped before the tolerance"),
                                            {{"minimiser", to_string(ro.method)}, {"|F|max", g6(ro.ftol) + " kcal/mol/Å"}, {"push-off", !ro.pushoff ? "off" : ro.pushoff_ramp_ps > 0 ? "MD " + g6(ro.pushoff_ramp_ps) + " ps, cap ramped, then minimised" : "on"}, {"force field", ffname}}, "",
                                            c, approx(energy, o.threads)));
          report(k, st, std::string(to_string(ro.method)) + " · |F|max " + g6(ro.ftol) + (rr.converged ? "" : " · not reached"), "done", 1);
        } else if (st == "md") {
          DynamicsOptions m;
          m.field = ff;
          m.energy = energy;
          m.dt = num(J, "dt", 1.0);
          m.respa = int(num(J, "respa", 1));   // r-RESPA: bonded forces every dt / respa
          const double ps = num(J, "ps", 10);
          m.steps = int64_t(ps * 1000 / m.dt);
          m.temperature = num(J, "temperature", 300);
          const std::string ens = text(J, "ensemble", "nvt");
          // thermostat: bussi (default) | langevin | nose-hoover; barostat: crescale (default) | berendsen | mtk (with
          // nose-hoover or NPH). NPH without MTK: Berendsen (no thermostat)
          const Thermostat th = thermostat_from_string(text(J, "thermostat", "bussi"));
          m.thermostat = ens == "nve" || ens == "nph" ? Thermostat::None : th;
          const std::string bs = text(J, "barostat", th == Thermostat::NoseHoover ? "mtk" : ens == "nph" ? "berendsen" : "crescale");
          m.barostat = ens == "npt" || ens == "nph" ? barostat_from_string(bs) : Barostat::None;
          if (ens == "nph" && m.barostat == Barostat::CRescale) m.barostat = Barostat::Berendsen;
          if (m.tau_t = num(J, "tau_t", 100); m.tau_t <= 0) m.tau_t = 100;
          if (m.tau_p = num(J, "tau_p", 1000); m.tau_p <= 0) m.tau_p = 1000;
          m.pressure = num(J, "pressure", 1.0);
          // constraints: none (default) | h-bonds | all-bonds; constraint_solver: shake (default) | lincs
          try {
            m.constraints = constraints_from_string(text(J, "constraints", "none"));
            m.constraint_algorithm = constraint_algorithm_from_string(text(J, "constraint_solver", "shake"));
          } catch (const std::exception& ex) { throw RecipeError(2, std::string("md: ") + ex.what()); }
          m.seed = seed_of(J);
          m.new_velocities = true;
          m.frame_every = 0;
          m.thermo_every = int(std::max<int64_t>(1, m.steps / 200));
          m.progress = [&](const ThermoRow& t) {
            report(k, st, g6(t.time_ps) + " ps · T " + g6(std::round(t.temperature)) + " K · ρ " + g6(t.density), "running", m.steps ? double(t.step) / double(m.steps) : 1);
            return true;
          };
          try { run_dynamics(s, m); } catch (const std::exception& e) { throw RecipeError(4, std::string("md: ") + e.what()); }
          std::vector<std::string> c = {"swope1982"};
          if (m.respa > 1) c.push_back("tuckerman1992");
          if (m.thermostat == Thermostat::Bussi) c.push_back("bussi2007");
          if (m.thermostat == Thermostat::NoseHoover) c.push_back("martyna1992");
          if (m.barostat == Barostat::CRescale) c.push_back("bernetti2020");
          if (m.barostat == Barostat::Berendsen) c.push_back("berendsen1984");
          if (m.barostat == Barostat::MTK) c.push_back("martyna1994");
          if (m.constraints != ConstraintMode::None) c.push_back(m.constraint_algorithm == ConstraintAlgorithm::Lincs ? "hess1997" : "ryckaert1977"), c.push_back("andersen1983");
          elec_cite(c, energy);
          KeyValues pr = {{"length", g6(ps) + " ps · " + std::to_string(m.steps) + " steps of " + g6(m.dt) + " fs" +
                                         (m.respa > 1 ? " (r-RESPA: bonded forces every " + g6(m.dt / m.respa) + " fs)" : "")},
                          {"temperature", g6(m.temperature) + " K"},
                          {"thermostat", m.thermostat == Thermostat::None ? std::string("none") : std::string(to_string(m.thermostat)) + " · τ " + g6(m.tau_t) + " fs"}};
          if (m.barostat != Barostat::None)
            pr.push_back({"barostat", std::string(to_string(m.barostat)) + " · " + g6(m.pressure) + " atm · τ " + g6(m.tau_p) + " fs" + (m.thermostat == Thermostat::None ? " (no thermostat: NPH)" : "")});
          if (m.constraints != ConstraintMode::None) {
            pr.push_back({"constraints", to_string(m.constraints)});
            pr.push_back({"constraint solver", m.constraint_algorithm == ConstraintAlgorithm::Lincs ? "LINCS" : "SHAKE"});
          }
          pr.push_back({"force field", ffname});
          res.manifest.steps.push_back(step("dynamics." + ens, ens == "npt" ? "NPT molecular dynamics" : ens == "nvt" ? "NVT molecular dynamics" : ens == "nph" ? "NPH molecular dynamics" : "NVE molecular dynamics", pr,
                                            seeded(m.seed), c, approx(energy, o.threads)));
          report(k, st, ens + " · " + g6(ps) + " ps · " + g6(m.temperature) + " K", "done", 1);
        } else {
          EquilibrateOptions e;
          ProtocolParams pp;
          pp.t_max = num(J, "t_max", pp.t_max);
          pp.t_final = num(J, "t_final", pp.t_final);
          if (J.has("p_max")) pp.p_max = J["p_max"].number() / 1.01325;   // bar → atm
          pp.time_scale = num(J, "time_scale", 1.0);
          const std::string proto = text(J, "protocol", "larsen21");
          try {
            // protocol_text: the stages written out (the Equilibrate page's text: "nvt 50 ps T 600" …) in place of a name
            e.stages = J.has("protocol_text") && J["protocol_text"].is_string() ? parse_protocol(J["protocol_text"].str()) : protocol_by_name(proto, pp);
            if (e.stages.empty()) throw std::invalid_argument("the protocol has no stages");
          } catch (const std::exception& ex) { throw RecipeError(2, std::string("equilibrate: ") + ex.what()); }
          e.md.field = ff;
          e.md.energy = energy;
          e.md.seed = seed_of(J);
          e.md.thermostat = thermostat_from_string(text(J, "thermostat", "bussi"));
          if (e.md.thermostat == Thermostat::None) e.md.thermostat = Thermostat::Bussi;
          e.md.barostat = barostat_from_string(text(J, "barostat", e.md.thermostat == Thermostat::NoseHoover ? "mtk" : "crescale"));
          if (e.md.barostat == Barostat::None) e.md.barostat = Barostat::CRescale;
          if (e.md.barostat == Barostat::MTK) e.md.thermostat = Thermostat::NoseHoover;
          try {
            e.md.constraints = constraints_from_string(text(J, "constraints", "none"));
            e.md.constraint_algorithm = constraint_algorithm_from_string(text(J, "constraint_solver", "shake"));
          } catch (const std::exception& ex) { throw RecipeError(2, std::string("equilibrate: ") + ex.what()); }
          e.until_converged = flag(J, "until_converged", false);
          if (J.has("block_ps")) e.block_ps = num(J, "block_ps", e.block_ps);
          if (J.has("max_blocks")) e.max_blocks = int(num(J, "max_blocks", e.max_blocks));
          if (J.has("tol_internal")) e.tol_internal = num(J, "tol_internal", e.tol_internal);
          // internal_target: "ris-pe" (Flory's polyethylene C_n at t_final, for alkane cells) or [value at n = 0, 1, 2 …]
          if (J.has("internal_target")) {
            const Json& it = J["internal_target"];
            if (it.is_string() && it.str() == "ris-pe") {
              const auto c = ris_cn(RisModel{}, pp.t_final, 2000);
              e.internal_target.assign(1, 0.0);
              e.internal_target.insert(e.internal_target.end(), c.begin(), c.end());
            } else if (it.is_array()) {
              for (const auto& x : it.items()) e.internal_target.push_back(x.number());
            } else {
              throw RecipeError(2, "equilibrate: internal_target is \"ris-pe\" or an array of values indexed by n");
            }
          }
          e.progress = [&](int si, int sn, const std::string& label, const ThermoRow&) {
            report(k, st, proto + " · step " + std::to_string(si + 1) + "/" + std::to_string(sn) + " · " + label, "running", sn ? double(si) / sn : 0);
            return true;
          };
          EquilibrateReport er;
          try { equilibrate(s, e, &er); } catch (const std::exception& ex) { throw RecipeError(4, std::string("equilibrate: ") + ex.what()); }
          std::vector<std::string> c = {"swope1982", e.md.thermostat == Thermostat::NoseHoover ? "martyna1992" : e.md.thermostat == Thermostat::Langevin ? "leimkuhler2013" : "bussi2007",
                                        e.md.barostat == Barostat::MTK ? "martyna1994" : e.md.barostat == Barostat::Berendsen ? "berendsen1984" : "bernetti2020"};
          if (proto == "larsen21") c.insert(c.begin(), "larsen2011");
          elec_cite(c, energy);
          double pmax = 0;
          for (const auto& sg : e.stages) pmax = std::max(pmax, sg.pressure);
          res.manifest.steps.push_back(step(proto == "larsen21" ? "equilibrate.larsen21" : "equilibrate.protocol", proto, {{"stages", std::to_string(e.stages.size())}, {"Pmax", g6(pmax) + " atm"},
                                            {"length", g6(er.ps) + " ps"}, {"force field", ffname}}, seeded(e.md.seed), c, approx(energy, o.threads)));
          report(k, st, proto + " · " + std::to_string(e.stages.size()) + " stages · " + g6(er.ps) + " ps", "done", 1);
        }
      } else if (st == "analyze") {
        report(k, st, "", "running", 0);
        std::vector<std::string> ids;
        if (J.has("properties") && J["properties"].is_array()) for (const auto& p : J["properties"].items()) ids.push_back(p.str());
        if (ids.empty()) ids = {"density"};
        const bool tg = std::find(ids.begin(), ids.end(), "tg") != ids.end();
        ids.erase(std::remove(ids.begin(), ids.end(), "tg"), ids.end());
        AnalyzeOptions ao;
        try { if (!ids.empty()) res.properties = analyze(as_trajectory(s), ids, ao); } catch (const std::exception& e) { throw RecipeError(2, std::string("analyze: ") + e.what()); }
        if (tg) {   // analyze: {properties: [tg], tg: {t_start, t_end, t_step, ps_per_step, equilibrate_ps, seed}} — a stepwise NPT cooling scan
          if (!ff) type_now(s);
          const Json T = J.has("tg") ? J["tg"] : Json::object();
          CoolingOptions co;
          co.field = ff;
          co.energy = energy;
          co.t_start = num(T, "t_start", co.t_start);
          co.t_end = num(T, "t_end", co.t_end);
          co.t_step = num(T, "t_step", co.t_step);
          co.ps_per_step = num(T, "ps_per_step", co.ps_per_step);
          co.equilibrate_ps = num(T, "equilibrate_ps", co.equilibrate_ps);
          co.seed = seed_of(T);
          co.new_velocities = true;
          const int steps = int(std::floor(std::fabs(co.t_start - co.t_end) / std::max(1e-9, co.t_step))) + 1;
          co.progress = [&](const ThermoRow& row, int si, int sn) {
            report(k, st, "Tg · " + std::to_string(si + 1) + "/" + std::to_string(sn > 0 ? sn : steps) + " · " + g6(std::round(row.temperature)) + " K", "running",
                   sn > 0 ? double(si) / sn : 0);
            return true;
          };
          System copy = s;
          CoolingResult cr;
          try { cr = run_cooling(copy, co); } catch (const std::exception& e) { throw RecipeError(4, std::string("tg: ") + e.what()); }
          for (auto& p : cooling_properties(cr)) res.properties.push_back(std::move(p));
          res.manifest.steps.push_back(step("analysis.tg", "stepwise NPT cooling scan, two-line fit of specific volume",
                                            {{"from", g6(co.t_start) + " K"}, {"to", g6(co.t_end) + " K"}, {"step", g6(co.t_step) + " K"},
                                             {"hold", g6(co.ps_per_step) + " ps"}, {"force field", ffname}}, seeded(co.seed), {"soldera2006", "bussi2007", "bernetti2020"},
                                            approx(energy, o.threads)));
        }
        std::string d;
        for (const auto& p : res.properties) if (std::isfinite(p.value)) d += (d.empty() ? "" : " · ") + p.id + " " + g6(p.value) + (p.unit.empty() ? "" : " " + p.unit);
        report(k, st, d, "done", 1);
      } else if (st == "export") {
        report(k, st, "", "running", 0);
        if (!ff && r.has("type")) type_now(s);   // the recipe's force field for the engine files, even without a run before
        std::vector<std::string> formats;
        if (J.is_array()) for (const auto& f : J.items()) formats.push_back(f.str());
        else if (J.is_string()) formats.push_back(J.str());
        else if (J.has("formats")) for (const auto& f : J["formats"].items()) formats.push_back(f.str());
        std::filesystem::create_directories(o.out_dir);
        const std::string stem = (std::filesystem::path(o.out_dir) / text(J, "name", res.name)).string();
        std::string first;
        for (const auto& f : formats) {
          std::string path;
          if (f == "lammps") {
            path = stem + ".data";
            if (ff) {
              LammpsStyle ls;   // the force field's own styles unless the recipe says lammps_styles: exact
              ls.native = text(J, "lammps_styles", "native") != "exact";
              ls.hybrid = J.has("hybrid") && J["hybrid"].kind() == Json::Bool && J["hybrid"].boolean();
              ls.units = text(J, "units", "auto");
              write_lammps_data_ff(s, *ff, energy, path, false, ls);
              write_lammps_input(s, *ff, energy, std::filesystem::path(path).filename().string(), stem + ".in", 0, true, {}, ls);
              res.files.push_back(stem + ".in");
            } else {
              write_lammps_data(s, path);
            }
          } else if (f == "moltemplate") {   // the LAMMPS system as moltemplate input (the same energies through moltemplate.sh -overlay-all)
            if (!ff) throw RecipeError(3, "export moltemplate: type the structure first (a type stage)");
            LammpsStyle ls;
            ls.native = text(J, "lammps_styles", "native") != "exact";
            const std::string dp = stem + ".lt.data", ip = stem + ".lt.in";
            write_lammps_data_ff(s, *ff, energy, dp, false, ls);
            write_lammps_input(s, *ff, energy, std::filesystem::path(dp).filename().string(), ip, 0, true, {}, ls);
            path = stem + ".lt";
            std::ofstream(path) << lammps_to_moltemplate(dp, ip, ff->name);
            std::filesystem::remove(dp);
            std::filesystem::remove(ip);
          } else if (f == "dlpoly") {   // STEM_dlpoly/FIELD, CONFIG, CONTROL (DL_POLY 4)
            if (!ff) throw RecipeError(3, "export dlpoly: type the structure first (a type stage)");
            DlpolyOptions dop;
            dop.title = text(J, "name", res.name);
            dop.cutoff = energy.cutoff;
            const std::string dir = stem + "_dlpoly";
            write_dlpoly(s, *ff, dir, dop);
            res.files.push_back(dir + "/CONFIG");
            res.files.push_back(dir + "/CONTROL");
            path = dir + "/FIELD";
          } else if (f == "gromacs" && ff) {   // topology, coordinates and a single-point .mdp with the same force field
            write_gromacs(s, *ff, energy, stem);
            res.files.push_back(stem + ".top");
            res.files.push_back(stem + ".mdp");
            path = stem + ".gro";
          } else if (f == "gromacs" || f == "gro") { path = stem + ".gro"; write_gro(s, path); }
          else if (f == "pdb") { path = stem + ".pdb"; write_pdb(s, path); }
          else if (f == "xyz") { path = stem + ".xyz"; write_xyz(s, path); }
          else if (f == "mol2") { path = stem + ".mol2"; write_mol2(s, path); }
          else if (f == "sdf") { path = stem + ".sdf"; write_sdf(s, path); }
          else if (f == "cif") { path = stem + ".cif"; write_cif(s, path); }
          else throw RecipeError(2, "export: unknown format '" + f + "' (lammps, moltemplate, gromacs, pdb, xyz, mol2, sdf, cif)");
          res.files.push_back(path);
          if (first.empty()) first = path;
        }
        if (!first.empty()) write_manifest(res.manifest, first);
        std::string d;
        for (const auto& f : formats) d += (d.empty() ? "" : ", ") + f;
        report(k, st, d + (first.empty() ? "" : " · provenance beside " + std::filesystem::path(first).filename().string()), "done", 1);
      }
    } catch (const RecipeError& e) {
      report(k, st, e.what(), "failed", 0);
      for (int w = k + 1; w < n; ++w) report(w, stages[size_t(w)], "", "waiting", 0);
      throw;
    } catch (const std::exception& e) {
      report(k, st, e.what(), "failed", 0);
      for (int w = k + 1; w < n; ++w) report(w, stages[size_t(w)], "", "waiting", 0);
      throw RecipeError(4, st + ": " + e.what());
    }
  }
  res.system = s;
  res.field = ff;
  res.forcefield = ff ? ffname : "";
  res.manifest.generator = by;
  return res;
}

}  // namespace caps
