#include "caps/recipe.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <cctype>
#include <set>

#include "caps/analysis.hpp"
#include "caps/config.hpp"
#include "caps/dynamics.hpp"
#include "caps/equilibrate.hpp"
#include "caps/ffdef.hpp"
#include "caps/field.hpp"
#include "caps/grow.hpp"
#include "caps/io.hpp"
#include "caps/molecule.hpp"
#include "caps/polymer.hpp"
#include "caps/relax.hpp"
#include "caps/typing.hpp"
#include "caps/uff.hpp"

namespace caps {

namespace {

const std::vector<std::string> kStages = {"build", "type", "grow", "relax", "md", "equilibrate", "analyze", "export"};
const std::set<std::string> kTop = {"recipe", "name", "build", "type", "grow", "relax", "md", "equilibrate", "analyze", "export", "electrostatics", "cutoff", "seed", "threads"};

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
    if (!kTop.count(k)) throw RecipeError(2, "unknown key '" + k + "' (stages: build, type, grow, relax, md, equilibrate, analyze, export)");
  if (r.has("recipe") && r["recipe"].is_number() && r["recipe"].number() != 1) throw RecipeError(2, "this CAPS reads recipe version 1");
  if (!r.has("build")) throw RecipeError(2, "a recipe needs a build stage (polymer, molecule or file)");
  std::vector<std::string> out;
  for (const auto& s : kStages) if (r.has(s)) out.push_back(s);
  if (r.has("grow") && !(r["build"].is_object() && r["build"].has("polymer"))) throw RecipeError(2, "grow needs build: { polymer: … }");
  return out;
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
  std::string ffname = "built-in default (GAFF for C and H, UFF otherwise)";
  // Types a structure with the recipe's force field (throws RecipeError 3 for untyped atoms or missing parameters).
  std::string borrowed;   // a library force field typed with its family's rules
  auto type_now = [&](const System& sys) {
    const Json T = r.has("type") ? r["type"] : Json::object();
    const std::string name = text(T, "forcefield", "default");
    std::string charges = text(T, "charges", "auto");
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
        std::vector<std::string> types;
        if (!def.typing.empty()) {
          const TypingResult tr = assign_types(sys, def);
          if (tr.untyped) throw RecipeError(3, std::to_string(tr.untyped) + " atoms match no typing rule of " + def.name);
          types = tr.types;
        } else {
          for (const auto& a : sys.atoms) types.push_back(a.name);
        }
        ParamReport rep;
        if (charges == "auto") {   // the force field's charges, the file's, or (GAFF-like fields have none per type) Gasteiger
          charges = sys.has_charges && !polymer ? "keep" : "types";
          try {
            ff = std::make_shared<ForceField>(parameterize(sys, def, types, charges, &rep, false));
          } catch (const FFError& e) {
            if (std::string(e.what()).find("has no charge for type") == std::string::npos) throw;
            charges = "gasteiger";
            rep = ParamReport{};
          }
        }
        if (!ff || charges == "gasteiger") ff = std::make_shared<ForceField>(parameterize(sys, def, types, charges, &rep, false));
        if (!rep.missing.empty()) throw RecipeError(3, std::to_string(rep.missing.size()) + " parameters missing in " + def.name + " (first: " + rep.missing.front() + ")");
        ffname = def.name;
      }
    } catch (const RecipeError&) { throw; } catch (const std::exception& e) { throw RecipeError(3, e.what()); }
    const std::string ch = charges == "qeq" ? "QEq" : charges == "gasteiger" ? "Gasteiger" : charges == "keep" ? "the file's" : charges == "auto" && ffname == "UFF" ? "no" : "from the force field";
    std::vector<std::string> c;
    if (ffname == "UFF") c.push_back("rappe1992");
    else if (ffname.find("GAFF") != std::string::npos) c.push_back("wang2004");
    if (charges == "qeq") c.push_back("rappe1991");
    if (charges == "gasteiger") c.push_back("gasteiger1980");
    KeyValues kv = {{"force field", ffname}, {"charges", ch}};
    if (!borrowed.empty()) kv.push_back({"typing rules", borrowed});
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
          report(k, st, "DP " + std::to_string(spec.dp) + " × " + std::to_string(chains) + " chains · unit " + info.formula + " · " + tac, "done", 1);
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
        } else {
          throw RecipeError(2, "build needs polymer, molecule or file");
        }
      } else if (st == "type") {
        report(k, st, "", "running", 0);
        System probe = s;
        if (polymer) {   // the structure comes later (grow): type one chain of the spec now, so a gap in the force field stops here
          GrowOptions g;
          g.chains = 1;
          g.seed = uint64_t(seed0);
          g.density = 0.05;
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
        g.curve = flag(J, "curve", true);
        const std::string method = text(J, "method", "trials");
        if (method != "trials" && method != "cbmc") throw RecipeError(2, "grow.method: trials (best-of-k torsion trials)");
        if (J.has("trials")) g.trials = int(J["trials"].number());
        g.progress = [&](int done, int total, int restarts) {
          report(k, st, "chain " + std::to_string(done) + "/" + std::to_string(total) + (restarts ? " · " + std::to_string(restarts) + " restarts" : ""), "running", total ? double(done) / total : 0);
          return true;
        };
        GrowReport gr;
        try { s = grow_chains(spec, g, &gr); } catch (const std::exception& e) { throw RecipeError(4, std::string("grow: ") + e.what()); }
        res.manifest.steps.push_back(step("grow.trials", std::to_string(chains) + " chains grown in a periodic cell, best-of-k trial placement by contact margin",
                                          {{"chains", std::to_string(chains)}, {"DP", std::to_string(spec.dp)}, {"trials", std::to_string(g.trials)}, {"density", g6(g.density) + " g/cm³"}}, seeded(g.seed), {"parsons2005", "matsumoto1998"}));
        report(k, st, "best of " + std::to_string(g.trials) + " trials · " + std::to_string(chains) + " chains · " + std::to_string(s.atoms.size()) + " atoms · box " + g6(gr.box) + " Å", "done", 1);
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
                                            {{"minimiser", to_string(ro.method)}, {"|F|max", g6(ro.ftol) + " kcal/mol/Å"}, {"push-off", ro.pushoff ? "on" : "off"}, {"force field", ffname}}, "",
                                            c, approx(energy, o.threads)));
          report(k, st, std::string(to_string(ro.method)) + " · |F|max " + g6(ro.ftol) + (rr.converged ? "" : " · not reached"), "done", 1);
        } else if (st == "md") {
          DynamicsOptions m;
          m.field = ff;
          m.energy = energy;
          m.dt = num(J, "dt", 1.0);
          const double ps = num(J, "ps", 10);
          m.steps = int64_t(ps * 1000 / m.dt);
          m.temperature = num(J, "temperature", 300);
          const std::string ens = text(J, "ensemble", "nvt");
          m.thermostat = ens == "nve" ? Thermostat::None : Thermostat::Bussi;
          m.barostat = ens == "npt" ? Barostat::CRescale : Barostat::None;
          m.pressure = num(J, "pressure", 1.0);
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
          if (m.thermostat == Thermostat::Bussi) c.push_back("bussi2007");
          if (m.barostat == Barostat::CRescale) c.push_back("bernetti2020");
          elec_cite(c, energy);
          KeyValues pr = {{"length", g6(ps) + " ps · " + std::to_string(m.steps) + " steps of " + g6(m.dt) + " fs"}, {"temperature", g6(m.temperature) + " K"},
                          {"thermostat", m.thermostat == Thermostat::Bussi ? "Bussi velocity rescaling · τ 100 fs" : "none"}};
          if (m.barostat != Barostat::None) pr.push_back({"barostat", "stochastic cell rescaling · " + g6(m.pressure) + " atm · τ 1000 fs"});
          pr.push_back({"force field", ffname});
          res.manifest.steps.push_back(step("dynamics." + ens, ens == "npt" ? "NPT molecular dynamics" : ens == "nvt" ? "NVT molecular dynamics" : "NVE molecular dynamics", pr,
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
          try { e.stages = protocol_by_name(proto, pp); } catch (const std::exception& ex) { throw RecipeError(2, ex.what()); }
          e.md.field = ff;
          e.md.energy = energy;
          e.md.seed = seed_of(J);
          e.until_converged = flag(J, "until_converged", false);
          e.progress = [&](int si, int sn, const std::string& label, const ThermoRow&) {
            report(k, st, proto + " · step " + std::to_string(si + 1) + "/" + std::to_string(sn) + " · " + label, "running", sn ? double(si) / sn : 0);
            return true;
          };
          EquilibrateReport er;
          try { equilibrate(s, e, &er); } catch (const std::exception& ex) { throw RecipeError(4, std::string("equilibrate: ") + ex.what()); }
          std::vector<std::string> c = {"swope1982", "bussi2007", "bernetti2020"};
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
        AnalyzeOptions ao;
        try { res.properties = analyze(as_trajectory(s), ids, ao); } catch (const std::exception& e) { throw RecipeError(2, std::string("analyze: ") + e.what()); }
        std::string d;
        for (const auto& p : res.properties) if (std::isfinite(p.value)) d += (d.empty() ? "" : " · ") + p.id + " " + g6(p.value) + (p.unit.empty() ? "" : " " + p.unit);
        report(k, st, d, "done", 1);
      } else if (st == "export") {
        report(k, st, "", "running", 0);
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
              write_lammps_data_ff(s, *ff, energy, path);
              write_lammps_input(s, *ff, energy, std::filesystem::path(path).filename().string(), stem + ".in");
              res.files.push_back(stem + ".in");
            } else {
              write_lammps_data(s, path);
            }
          } else if (f == "gromacs" || f == "gro") { path = stem + ".gro"; write_gro(s, path); }
          else if (f == "pdb") { path = stem + ".pdb"; write_pdb(s, path); }
          else if (f == "xyz") { path = stem + ".xyz"; write_xyz(s, path); }
          else if (f == "mol2") { path = stem + ".mol2"; write_mol2(s, path); }
          else throw RecipeError(2, "export: unknown format '" + f + "' (lammps, gromacs, pdb, xyz, mol2)");
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
