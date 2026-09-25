#include "caps_c.h"

#include <cstdio>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <string>

#include "caps/analysis.hpp"
#include "caps/dynamics.hpp"
#include "caps/elements.hpp"
#include "caps/equilibrate.hpp"
#include "caps/ffdef.hpp"
#include "caps/grow.hpp"
#include "caps/io.hpp"
#include "caps/mechanics.hpp"
#include "caps/molecule.hpp"
#include "caps/config.hpp"
#include "caps/bench.hpp"
#include "caps/polymer.hpp"
#include "caps/pack.hpp"
#include "caps/properties.hpp"
#include "caps/react.hpp"
#include "caps/relax.hpp"
#include "caps/render.hpp"
#include "caps/typing.hpp"
#include "caps/uff.hpp"
#include "caps/crystal.hpp"
#include "caps/nano.hpp"
#include "caps/json.hpp"

#include <map>
#include <memory>
#include <set>
#include <sstream>

using caps::operator+;
using caps::operator-;
using caps::operator*;

// CAPS Field state of a document: the force field chosen from the library, parameters imported or entered by hand on
// top of it, per-atom type overrides, and the result (types, charges, parameters, what is missing).
struct FieldState {
  std::string ff_path;
  caps::FFDef base;                          // the library force field, with its typing rules
  caps::FFDef extra;                         // imported and hand-entered parameter rules (win over the base)
  std::vector<std::string> imported;         // files the imported rules came from
  std::map<int32_t, std::string> overrides;  // atom → type set by hand
  std::string charges = "types";             // types (force field), gasteiger, keep (from the file)
  caps::TypingResult typing;
  std::vector<std::string> types;
  caps::ParamReport rep;
  std::shared_ptr<const caps::ForceField> ff;   // null while atoms are untyped
  bool complete = false;
  std::string report;                        // JSON, see caps_field_report
  // the file's own types, restored by caps_field_clear
  std::vector<std::pair<int, std::string>> file_types;
  std::vector<caps::TypeInfo> file_type_table;
  std::vector<double> file_charges;
  bool file_has_charges = false;
};

struct caps_doc {
  caps::Trajectory traj;
  caps::System frame;
  caps::Renderer renderer;
  std::vector<double> dcom;   // distance of each atom to its own molecule's centre of mass
  size_t current = 0;
  bool wrap = false;
  std::unique_ptr<FieldState> field;
  std::string analysis;   // last caps_analyze result (JSON)
  std::string eq_checks;  // last caps_equilibrate convergence checks (JSON)
  int64_t held_mol = 0;   // molecule held in place by caps_relax (0: none)
};

namespace {

thread_local std::string g_error;

// Electrostatics chosen in the Studio's Settings (caps_set_electrostatics): applied to every energy evaluation here.
struct Elec {
  int mode = 0;   // 0 DSF, 1 PME
  double rtol = 1e-5, spacing = 1.0;
  int order = 5;
} g_elec;

caps::EnergyOptions elec(caps::EnergyOptions e = {}) {
  e.electrostatics = g_elec.mode == 1 ? caps::EnergyOptions::Electrostatics::PME : caps::EnergyOptions::Electrostatics::DSF;
  e.ewald_rtol = g_elec.rtol;
  e.pme_spacing = g_elec.spacing;
  e.pme_order = g_elec.order;
  return e;
}

template <class F>
int32_t guard(F&& f) {
  try {
    return f();
  } catch (const std::exception& e) {
    g_error = e.what();
  } catch (...) {
    g_error = "unknown error";
  }
  return -1;
}

void refresh(caps_doc* d) {
  d->frame = d->traj.frame(d->current);
  if (d->current + 1 != d->traj.frames()) d->frame.velocities.clear();   // velocities belong to the last frame
  if (!d->frame.unwrapped) caps::make_molecules_whole(d->frame);
  const auto shapes = caps::molecule_shapes(d->frame);   // always from whole molecules
  const auto mol = d->frame.molecules();
  d->dcom.resize(d->frame.atoms.size());
  for (size_t i = 0; i < d->frame.atoms.size(); ++i) d->dcom[i] = caps::norm(d->frame.atoms[i].pos - shapes[mol[i]].com);
  if (d->wrap && d->frame.cell.valid()) {
    for (auto& a : d->frame.atoms) a.pos = d->frame.cell.wrap(a.pos);
    d->frame.unwrapped = false;
  }
}

caps::Camera cam_of(const caps_camera* c) {
  caps::Camera k;
  if (!c) return k;
  k.yaw = c->yaw; k.pitch = c->pitch; k.zoom = c->zoom > 0 ? c->zoom : 1; k.pan_x = c->pan_x; k.pan_y = c->pan_y; k.perspective = c->perspective != 0;
  return k;
}

caps::RenderOptions opts_of(const caps_doc* d, const caps_render_opts* o) {
  caps::RenderOptions r;
  if (!o) return r;
  r.width = std::clamp(o->width, 16, 16384);
  r.height = std::clamp(o->height, 16, 16384);
  r.supersample = o->supersample > 0 ? o->supersample : 2;
  r.background = static_cast<caps::Background>(std::clamp(o->background, 0, 3));
  r.custom_rgb = o->custom_rgb;
  r.colour_by = static_cast<caps::ColourBy>(std::clamp(o->colour_by, 0, 3));
  r.style = static_cast<caps::Style>(std::clamp(o->style, 0, 4));
  r.outlines = o->outlines != 0;
  r.depth_cue = o->depth_cue != 0;
  r.show_cell = o->show_cell != 0;
  for (int k = 0; k < 4; ++k) if (o->highlight[k] >= 0) r.highlight.push_back(o->highlight[k]);
  if (r.colour_by == caps::ColourBy::Property) r.property = d->dcom;
  return r;
}

caps::ProtocolParams protocol_params(const caps_protocol_params* p) {
  caps::ProtocolParams q;
  if (!p) return q;
  if (p->t_final > 0) q.t_final = p->t_final;
  if (p->t_max > 0) q.t_max = p->t_max;
  q.p_final = p->p_final;
  if (p->p_max > 0) q.p_max = p->p_max;
  if (p->time_scale > 0) q.time_scale = p->time_scale;
  if (p->cycles > 0) q.cycles = p->cycles;
  if (p->t_low > 0) q.t_low = p->t_low;
  if (p->t_high > 0) q.t_high = p->t_high;
  if (p->ramp_ps > 0) q.ramp_ps = p->ramp_ps;
  if (p->hold_ps > 0) q.hold_ps = p->hold_ps;
  return q;
}

// The force field for Relax / Dynamics: the Field assignment when there is one (it must be complete: CAPS never guesses
// parameters), otherwise null (the built-in GAFF typing of C and H).
// The force field used when none is assigned in Field: the built-in GAFF for C/H structures, UFF for any other.
caps::ForceField default_ff(const caps::System& s) { return caps::default_forcefield(s); }

std::shared_ptr<const caps::ForceField> field_for_run(const caps_doc* d) {
  if (!d->field) {
    // no force field assigned: the built-in GAFF covers C and H; anything else runs with UFF (every element)
    const auto& atoms = d->traj.topology.atoms;
    if (std::any_of(atoms.begin(), atoms.end(), [](const caps::Atom& a) { return a.element != 1 && a.element != 6; })) {
      caps::System s = d->traj.frame(d->current);
      return std::make_shared<caps::ForceField>(caps::assign_uff(s));
    }
    return nullptr;
  }
  const FieldState& F = *d->field;
  if (!F.complete) {
    int untyped = 0;
    for (const auto& t : F.types) untyped += t.empty();
    throw caps::FieldError("CAPS Field: " + F.base.name + " is incomplete for this structure (" + std::to_string(untyped) + " atoms untyped, " +
                           std::to_string(F.rep.missing.size()) + " parameters missing); complete it in the Field panel or clear the assignment");
  }
  return F.ff;
}

std::string hex_colour(unsigned c) {
  char b[8];
  std::snprintf(b, sizeof b, "#%06X", c & 0xFFFFFF);
  return b;
}

// Types and parameterises the current structure with the field state, writes the types (and charges) into the
// document so the viewer colours by force-field type, and builds the JSON report.
void field_run(caps_doc* d) {
  FieldState& F = *d->field;
  caps::FFDef def = F.base;
  caps::merge_forcefield(def, F.extra);
  const caps::System& s = d->frame;
  const size_t n = s.atoms.size();
  const bool uff = caps::is_uff(F.ff_path);
  bool rules = !def.typing.empty();
  if (uff) {   // UFF: typed from elements, hybridisation and oxidation state
    F.typing = caps::TypingResult{};
    F.typing.types = caps::uff_types(s, &F.typing.why);
    F.typing.rule.assign(n, -1);
    F.typing.candidates.assign(n, {});
  } else if (rules) {
    F.typing = caps::assign_types(s, def);
  } else {   // no typing rules: the file's atom names are the types
    F.typing = caps::TypingResult{};
    F.typing.types.resize(n);
    F.typing.why.assign(n, "type name from the file");
    F.typing.rule.assign(n, -1);
    F.typing.candidates.assign(n, {});
    for (size_t i = 0; i < n; ++i) F.typing.types[i] = s.atoms[i].name;
  }
  F.types = F.typing.types;
  for (const auto& [i, t] : F.overrides)
    if (i >= 0 && size_t(i) < n) F.types[i] = t;
  std::set<std::string> known;
  for (const auto& t : def.types) known.insert(t.name);
  int untyped = 0;
  for (auto& t : F.types)
    if (t.empty() || !known.count(t)) { t.clear(); ++untyped; }
  F.rep = caps::ParamReport{};
  F.ff.reset();
  if (!untyped && uff) {
    caps::UffOptions uo;
    uo.keep_charges = F.charges == "keep";
    uo.labels = F.types;
    F.ff = std::make_shared<caps::ForceField>(caps::assign_uff(s, uo));
    for (const auto& note : F.ff->notes) F.rep.notes.push_back(note);
  } else if (!untyped) {
    F.ff = std::make_shared<caps::ForceField>(caps::parameterize(s, def, F.types, F.charges, &F.rep, true));
  }
  F.complete = F.ff && F.rep.missing.empty();

  // types into the document: colour by type shows the force-field types
  std::vector<std::string> names;
  std::map<std::string, int> tix;
  for (size_t i = 0; i < n; ++i) {
    const std::string t = F.types[i].empty() ? "?" : F.types[i];
    auto [it, fresh] = tix.emplace(t, int(names.size()) + 1);
    if (fresh) names.push_back(t);
    d->traj.topology.atoms[i].type = it->second;
    d->traj.topology.atoms[i].name = t;
    if (F.ff) d->traj.topology.atoms[i].charge = F.ff->charge[i];
  }
  d->traj.topology.types.clear();
  for (size_t k = 0; k < names.size(); ++k) {
    caps::TypeInfo ti;
    ti.type = int(k) + 1;
    ti.label = names[k];
    const caps::FFType* ft = def.type(names[k]);
    ti.mass = ft && ft->mass > 0 ? ft->mass : 0;
    d->traj.topology.types.push_back(ti);
  }
  if (F.ff) d->traj.topology.has_charges = true;
  refresh(d);

  // report
  caps::Json r = caps::Json::object();
  r["forcefield"] = def.name;
  r["version"] = def.version;
  r["source"] = def.source;
  r["file"] = F.ff_path;
  r["typing"] = uff ? std::string("UFF typer: element, hybridisation, conjugation, oxidation state")
                    : rules ? (def.typing_source.empty() ? F.ff_path : def.typing_source) : std::string("atom names in the file");
  r["rules"] = double(def.typing.size());
  r["charges"] = F.charges;
  caps::Json refs = caps::Json::array();
  for (const auto& x : def.references) refs.push_back(x);
  r["references"] = refs;
  const std::string rule_src = [&] {
    std::string p = def.typing_source.empty() ? F.ff_path : def.typing_source;
    const auto k = p.find_last_of("/\\");
    return k == std::string::npos ? p : p.substr(k + 1);
  }();
  caps::Json atoms = caps::Json::array();
  double qsum = 0;
  int overridden = 0;
  for (size_t i = 0; i < n; ++i) {
    caps::Json a = caps::Json::object();
    a["i"] = double(i + 1);
    a["el"] = std::string(caps::element(s.atoms[i].element).symbol);
    a["type"] = F.types[i];
    const bool ov = F.overrides.count(int32_t(i)) > 0;
    overridden += ov;
    a["ov"] = ov;
    const int ri = i < F.typing.rule.size() ? F.typing.rule[i] : -1;
    if (ov) {
      a["rule"] = "set by hand";
      a["src"] = "override";
    } else if (ri >= 0) {
      const auto& tr = def.typing[size_t(ri)];
      a["rule"] = tr.smarts;
      a["desc"] = tr.description;
      a["prio"] = double(tr.priority);
      a["src"] = rule_src;
    } else if (uff) {
      a["rule"] = "UFF typer";
      a["src"] = "UFF";
    } else {
      a["rule"] = rules ? "no rule matched" : "file";
      a["src"] = rules ? "" : "file";
    }
    a["why"] = i < F.typing.why.size() ? F.typing.why[i] : "";
    caps::Json c = caps::Json::array();
    if (i < F.typing.candidates.size())
      for (const auto& x : F.typing.candidates[i]) c.push_back(x);
    a["cands"] = c;
    if (F.ff) {
      a["q"] = F.ff->charge[i];
      qsum += F.ff->charge[i];
    }
    atoms.push_back(a);
  }
  r["atoms"] = atoms;
  r["typed"] = double(n - size_t(untyped));
  r["untyped"] = double(untyped);
  r["overridden"] = double(overridden);
  r["ambiguous"] = double(F.typing.ambiguous);
  r["net_charge"] = qsum;
  r["has_charges"] = F.ff != nullptr;
  caps::Json miss = caps::Json::array();
  for (const auto& m : F.rep.missing) miss.push_back(m);
  r["missing"] = miss;
  int estimated = 0, imported = 0;
  for (const auto& [k, v] : F.rep.used) {
    const auto sp = k.find(' ');
    const std::string name = sp == std::string::npos ? k : k.substr(sp + 1);
    if (name.rfind("user:", 0) == 0) estimated += v;
    if (name.rfind("imported:", 0) == 0) imported += v;
  }
  r["estimated"] = double(estimated);
  r["imported"] = double(imported);
  caps::Json hand = caps::Json::array();
  for (const auto* v : {&F.extra.pairs, &F.extra.bonds, &F.extra.angles, &F.extra.dihedrals, &F.extra.impropers})
    for (const auto& x : *v)
      if (x.name.rfind("user:", 0) == 0) {
        std::string t = x.name.substr(6) + " · " + (x.style.empty() ? "default" : x.style) + " ·";
        for (double p : x.params) { char b[32]; std::snprintf(b, sizeof b, " %g", p); t += b; }
        hand.push_back(t);
      }
  r["entered"] = hand;
  caps::Json imp = caps::Json::array();
  for (const auto& x : F.imported) imp.push_back(x);
  r["imported_files"] = imp;
  // types present, with the viewer's colour for each
  caps::Json used = caps::Json::array();
  std::map<std::string, int> count;
  for (const auto& t : F.types) ++count[t.empty() ? "?" : t];
  for (size_t k = 0; k < names.size(); ++k) {
    caps::Json u = caps::Json::object();
    u["name"] = names[k];
    u["count"] = double(count[names[k]]);
    u["colour"] = hex_colour(caps::molecule_colour(int(k)));
    const caps::FFType* ft = def.type(names[k]);
    if (ft) u["desc"] = ft->description;
    used.push_back(u);
  }
  r["used"] = used;
  caps::Json all = caps::Json::array();   // every type of the force field, for overrides
  for (const auto& t : def.types) {
    caps::Json u = caps::Json::object();
    u["name"] = t.name;
    u["el"] = std::string(caps::element(t.element).symbol);
    u["desc"] = t.description;
    all.push_back(u);
  }
  r["fftypes"] = all;
  caps::Json style = caps::Json::object();   // default styles, for entering parameters by hand
  style["pair"] = def.pair_style;
  style["bond"] = def.bond_style;
  style["angle"] = def.angle_style;
  style["dihedral"] = def.dihedral_style;
  style["improper"] = def.improper_style;
  r["styles"] = style;
  caps::Json notes = caps::Json::array();
  for (const auto& x : F.typing.notes) notes.push_back(x);
  for (const auto& x : F.rep.notes) notes.push_back(x);
  r["notes"] = notes;
  r["complete"] = F.complete;
  if (F.ff) {
    caps::Evaluator ev(*F.ff, elec());
    std::vector<double> x, f;
    for (const auto& a : s.atoms) x.insert(x.end(), a.pos.begin(), a.pos.end());
    const caps::EnergyTerms e = ev.compute(x, s.cell, f);
    caps::Json en = caps::Json::object();
    en["bond"] = e.bond; en["angle"] = e.angle; en["dihedral"] = e.dihedral; en["improper"] = e.improper;
    en["vdw"] = e.vdw; en["coulomb"] = e.coulomb; en["total"] = e.total();
    r["energy"] = en;
  }
  F.report = r.dump(0);
}

int32_t report_out(const std::string& t, char* buf, int32_t cap) {
  if (buf && cap > 0) {
    const size_t k = std::min(t.size(), size_t(cap) - 1);
    std::memcpy(buf, t.data(), k);
    buf[k] = 0;
  }
  return int32_t(t.size() + 1);
}

}  // namespace

extern "C" {

int32_t caps_abi_version(void) { return CAPS_ABI_VERSION; }
const char* caps_last_error(void) { return g_error.c_str(); }

caps_doc* caps_open(const char* path, const char* topology_path) {
  try {
    auto* d = new caps_doc;
    d->traj = caps::open_file(path, topology_path ? topology_path : "");
    refresh(d);
    return d;
  } catch (const std::exception& e) {
    g_error = e.what();
  } catch (...) {
    g_error = "unknown error";
  }
  return nullptr;
}

void caps_close(caps_doc* d) { delete d; }

caps_doc* caps_grow(const caps_grow_opts* o, caps_progress_fn progress, void* user, char* report, int32_t cap) {
  try {
    caps::GrowOptions g;
    g.chains = o->chains;
    g.dp = o->dp;
    g.tacticity = o->tacticity == 1 ? caps::Tacticity::Isotactic : o->tacticity == 2 ? caps::Tacticity::Syndiotactic : caps::Tacticity::Atactic;
    g.seed = o->seed;
    g.box = o->box;
    g.density = o->density;
    g.contact_scale = o->contact_scale > 0 ? o->contact_scale : 1.0;
    g.curve = o->curve != 0;
    if (progress) g.progress = [&](int done, int total, int restarts) { return progress(done, total, restarts, user) == 0; };
    caps::GrowReport rep;
    caps::System s = caps::grow(g, &rep);
    auto* d = new caps_doc;
    d->traj.topology = s;
    std::vector<caps::Vec3> p;
    for (const auto& a : s.atoms) p.push_back(a.pos);
    d->traj.positions.push_back(std::move(p));
    d->traj.cells.push_back(s.cell);
    d->traj.timesteps.push_back(0);
    refresh(d);
    if (report && cap > 0) {
      std::string t;
      for (const auto& n : rep.notes) t += n + "\n";
      std::strncpy(report, t.c_str(), size_t(cap) - 1);
      report[cap - 1] = 0;
    }
    return d;
  } catch (const std::exception& e) {
    g_error = e.what();
  } catch (...) {
    g_error = "unknown error";
  }
  return nullptr;
}

int32_t caps_save(caps_doc* d, const char* path) {
  return guard([&] {
    const std::string p = path;
    auto ends = [&](const char* e) { const std::string x = e; return p.size() >= x.size() && p.compare(p.size() - x.size(), x.size(), x) == 0; };
    if (ends(".pdb")) caps::write_pdb(d->frame, p);
    else if (ends(".xyz")) caps::write_xyz(d->frame, p);
    else if (ends(".mol2")) caps::write_mol2(d->frame, p);
    else if (d->field) {   // the Field assignment: its coefficients when complete, else the structure alone
      if (d->field->complete) caps::write_lammps_data_ff(d->frame, *d->field->ff, elec(), p);
      else caps::write_lammps_data(d->frame, p);
    } else {
      caps::ForceField ff;
      bool typed = true;
      try { ff = default_ff(d->frame); } catch (const caps::FieldError&) { typed = false; }
      if (typed) caps::write_lammps_data_ff(d->frame, ff, elec(), p);
      else caps::write_lammps_data(d->frame, p);
    }
    return 0;
  });
}

int32_t caps_lammps_input(caps_doc* d, const char* data_name, char* text, int32_t cap) {
  return guard([&] {
    caps::ForceField ff;
    if (d->field && d->field->complete) ff = *d->field->ff;
    else ff = default_ff(d->frame);
    const auto tmp = std::filesystem::temp_directory_path() / ("caps_input_" + std::to_string(reinterpret_cast<uintptr_t>(d)) + ".in");
    caps::write_lammps_input(d->frame, ff, elec(), data_name && *data_name ? data_name : "system.data", tmp.string());
    std::ifstream in(tmp);
    std::string s((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    in.close();
    std::filesystem::remove(tmp);
    // keep the setup only: the deck's thermo and run 0 lines are for single-point checks
    std::string out;
    std::istringstream ls(s);
    for (std::string line; std::getline(ls, line);)
      if (line.rfind("thermo_style", 0) != 0 && line.rfind("thermo_modify", 0) != 0 && line.rfind("run ", 0) != 0) out += line + "\n";
    const int32_t need = int32_t(out.size() + 1);
    if (text && cap > 0) {
      const size_t m = std::min<size_t>(size_t(cap - 1), out.size());
      std::memcpy(text, out.data(), m);
      text[m] = 0;
    }
    return need;
  });
}

int32_t caps_relax(caps_doc* d, const caps_relax_opts* o, caps_relax_progress_fn progress, void* user, char* report, int32_t cap) {
  return guard([&] {
    caps::RelaxOptions r;
    r.field = field_for_run(d);
    r.method = static_cast<caps::Minimiser>(std::clamp(o->method, 0, 3));
    if (o->ftol > 0) r.ftol = o->ftol;
    if (o->max_iterations > 0) r.max_iterations = o->max_iterations;
    r.target_density = o->target_density;
    if (o->compress_step > 0) r.compress_step = o->compress_step;
    r.pushoff = o->pushoff != 0;
    r.relax_box = o->relax_box != 0;
    r.pressure = o->pressure;
    if (o->cutoff > 0) r.energy.cutoff = o->cutoff;
    r.energy.coulomb = o->coulomb != 0;
    r.energy = elec(r.energy);
    r.energy.threads = o->threads;
    if (progress)
      r.progress = [&](const caps::RelaxProgress& p) {
        return progress(p.stage_index, p.stages, p.iteration, p.energy, p.fmax, p.density, user) == 0;
      };
    caps::System s = d->traj.frame(d->current);
    if (!s.unwrapped) caps::make_molecules_whole(s);
    if (d->held_mol > 0) {
      r.fixed.assign(s.atoms.size(), 0);
      for (size_t i = 0; i < s.atoms.size(); ++i) r.fixed[i] = s.atoms[i].mol == d->held_mol;
    }
    caps::Trajectory out;
    out.topology = s;
    auto push = [&](const std::vector<caps::Vec3>& p, const caps::Cell& c) {
      out.positions.push_back(p);
      out.cells.push_back(c);
      out.timesteps.push_back(int64_t(out.timesteps.size()));
    };
    {
      std::vector<caps::Vec3> p;
      for (const auto& a : s.atoms) p.push_back(a.pos);
      push(p, s.cell);
    }
    r.snapshot = [&](const std::vector<double>& x, const caps::Cell& c, const std::string&) {
      std::vector<caps::Vec3> p(x.size() / 3);
      for (size_t i = 0; i < p.size(); ++i) p[i] = {x[3 * i], x[3 * i + 1], x[3 * i + 2]};
      push(p, c);
    };
    caps::RelaxReport rep;
    caps::relax(s, r, &rep);
    // the final structure is the last snapshot; keep charges and the final cell on the topology
    out.topology.atoms = s.atoms;
    out.topology.cell = s.cell;
    out.topology.velocities.clear();   // minimised positions: old velocities no longer belong to them
    out.topology.has_charges = true;
    out.topology.unwrapped = true;
    out.topology.notes = rep.notes;
    d->traj = std::move(out);
    d->current = d->traj.frames() - 1;
    refresh(d);
    if (report && cap > 0) {
      std::string t = rep.field + "\n";
      for (const auto& n : rep.notes) t += n + "\n";
      for (const auto& st : rep.stages) {
        char b[200];
        std::snprintf(b, sizeof b, "%s: %d it, E %.1f, |F|max %.3f (%s)\n", st.name.c_str(), st.iterations, st.energy, st.fmax, st.stopped_by.c_str());
        t += b;
      }
      std::strncpy(report, t.c_str(), size_t(cap) - 1);
      report[cap - 1] = 0;
    }
    return rep.converged ? 0 : 1;
  });
}

int32_t caps_md(caps_doc* d, const caps_md_opts* o, caps_md_progress_fn progress, void* user, char* report, int32_t cap) {
  return guard([&] {
    caps::DynamicsOptions m;
    m.field = field_for_run(d);
    if (o->dt > 0) m.dt = o->dt;
    m.steps = std::max<int64_t>(0, o->steps);
    m.temperature = o->temperature;
    m.thermostat = static_cast<caps::Thermostat>(std::clamp(o->thermostat, 0, 2));
    if (o->tau_t > 0) m.tau_t = o->tau_t;
    m.barostat = static_cast<caps::Barostat>(std::clamp(o->barostat, 0, 2));
    m.pressure = o->pressure;
    if (o->tau_p > 0) m.tau_p = o->tau_p;
    m.new_velocities = o->new_velocities != 0;
    m.seed = o->seed;
    if (o->thermo_every > 0) m.thermo_every = o->thermo_every;
    m.frame_every = std::max(0, o->frame_every);
    if (o->cutoff > 0) m.energy.cutoff = o->cutoff;
    m.energy.coulomb = o->coulomb != 0;
    m.energy = elec(m.energy);
    m.energy.tail = o->tail != 0;
    m.energy.threads = o->threads;
    if (progress)
      m.progress = [&](const caps::ThermoRow& r) {
        caps_thermo t{r.step, r.time_ps, r.temperature, r.potential, r.kinetic, r.total, r.conserved, r.pressure, r.volume, r.density};
        return progress(&t, m.steps, user) == 0;
      };
    caps::System s = d->traj.frame(d->current);
    if (d->current + 1 != d->traj.frames()) s.velocities.clear();   // velocities belong to the last frame only
    if (d->held_mol > 0) {
      m.fixed.assign(s.atoms.size(), 0);
      for (size_t i = 0; i < s.atoms.size(); ++i) m.fixed[i] = s.atoms[i].mol == d->held_mol;
    }
    if (!s.unwrapped) caps::make_molecules_whole(s);
    caps::Trajectory out;
    out.topology = s;
    m.frame = [&](const std::vector<double>& x, const caps::Cell& c, int64_t step) {
      std::vector<caps::Vec3> p(x.size() / 3);
      for (size_t i = 0; i < p.size(); ++i) p[i] = {x[3 * i], x[3 * i + 1], x[3 * i + 2]};
      out.positions.push_back(std::move(p));
      out.cells.push_back(c);
      out.timesteps.push_back(step);
    };
    if (m.frame_every <= 0) m.frame_every = int(std::max<int64_t>(1, m.steps));   // at least the start and the end
    caps::DynamicsReport rep;
    caps::run_dynamics(s, m, &rep);
    out.topology.atoms = s.atoms;
    out.topology.cell = s.cell;
    out.topology.velocities = s.velocities;
    out.topology.has_charges = true;
    out.topology.unwrapped = true;
    out.topology.notes = rep.notes;
    d->traj = std::move(out);
    d->current = d->traj.frames() - 1;
    refresh(d);
    if (report && cap > 0) {
      std::string t;
      for (const auto& n : rep.notes) t += n + "\n";
      if (!rep.thermo.empty()) {
        const auto& r = rep.thermo.back();
        char b[256];
        std::snprintf(b, sizeof b, "end: T %.1f K, P %.0f atm, density %.4f g/cm³, E %.1f kcal/mol\n", r.temperature, r.pressure, r.density, r.total);
        t += b;
      }
      std::strncpy(report, t.c_str(), size_t(cap) - 1);
      report[cap - 1] = 0;
    }
    return 0;
  });
}

int32_t caps_save_trajectory(caps_doc* d, const char* path) {
  return guard([&] {
    caps::write_lammps_dump(d->traj, path);
    return 0;
  });
}


int32_t caps_protocol_text(const char* name, const caps_protocol_params* p, char* out, int32_t cap) {
  return guard([&] {
    const std::string t = caps::protocol_text(caps::protocol_by_name(name, protocol_params(p)));
    if (out && cap > 0) {
      std::strncpy(out, t.c_str(), size_t(cap) - 1);
      out[cap - 1] = 0;
    }
    return int32_t(t.size());
  });
}

int32_t caps_equilibrate(caps_doc* d, const char* protocol, const caps_equil_opts* o, caps_equil_progress_fn progress, void* user,
                         char* report, int32_t cap) {
  return guard([&] {
    caps::EquilibrateOptions e;
    e.md.field = field_for_run(d);
    e.stages = caps::parse_protocol(protocol ? protocol : "");
    if (o->dt > 0) e.md.dt = o->dt;
    e.md.thermostat = o->thermostat == 2 ? caps::Thermostat::Langevin : caps::Thermostat::Bussi;
    e.md.barostat = o->barostat == 2 ? caps::Barostat::Berendsen : caps::Barostat::CRescale;
    if (o->tau_t > 0) e.md.tau_t = o->tau_t;
    if (o->tau_p > 0) e.md.tau_p = o->tau_p;
    e.md.seed = o->seed;
    if (o->cutoff > 0) e.md.energy.cutoff = o->cutoff;
    e.md.energy.coulomb = o->coulomb != 0;
    e.md.energy = elec(e.md.energy);
    e.md.energy.tail = o->tail != 0;
    e.md.energy.threads = o->threads;
    if (o->frame_ps > 0) e.frame_ps = o->frame_ps;
    if (o->thermo_ps > 0) e.thermo_ps = o->thermo_ps;
    e.until_converged = o->until_converged != 0;
    if (o->block_ps > 0) e.block_ps = o->block_ps;
    if (o->max_blocks > 0) e.max_blocks = o->max_blocks;
    if (o->tol_density > 0) e.tol_density = o->tol_density;
    if (o->tol_energy > 0) e.tol_energy = o->tol_energy;
    if (o->tol_rg > 0) e.tol_rg = o->tol_rg;
    if (progress)
      e.progress = [&](int st, int n, const std::string& label, const caps::ThermoRow& r) {
        caps_thermo t{r.step, r.time_ps, r.temperature, r.potential, r.kinetic, r.total, r.conserved, r.pressure, r.volume, r.density};
        return progress(st, n, label.c_str(), &t, user) == 0;
      };
    caps::System s = d->traj.frame(d->current);
    if (d->current + 1 != d->traj.frames()) s.velocities.clear();
    if (!s.unwrapped) caps::make_molecules_whole(s);
    if (d->held_mol > 0) {
      e.md.fixed.assign(s.atoms.size(), 0);
      for (size_t i = 0; i < s.atoms.size(); ++i) e.md.fixed[i] = s.atoms[i].mol == d->held_mol;
    }
    caps::Trajectory out;
    out.topology = s;
    e.frame = [&](const std::vector<double>& x, const caps::Cell& c, int64_t step) {
      std::vector<caps::Vec3> p(x.size() / 3);
      for (size_t i = 0; i < p.size(); ++i) p[i] = {x[3 * i], x[3 * i + 1], x[3 * i + 2]};
      out.positions.push_back(std::move(p));
      out.cells.push_back(c);
      out.timesteps.push_back(step);
    };
    caps::EquilibrateReport rep;
    caps::equilibrate(s, e, &rep);
    out.topology.atoms = s.atoms;
    out.topology.cell = s.cell;
    out.topology.velocities = s.velocities;
    out.topology.has_charges = true;
    out.topology.unwrapped = true;
    out.topology.notes = rep.notes;
    d->traj = std::move(out);
    d->current = d->traj.frames() - 1;
    refresh(d);
    {
      std::ostringstream j;
      j << "{\"converged\":" << (rep.converged ? "true" : "false") << ",\"blocks\":" << rep.blocks << ",\"block_ps\":" << e.block_ps << ",\"checks\":[";
      for (size_t k = 0; k < rep.checks.size(); ++k) {
        const auto& c = rep.checks[k];
        j << (k ? "," : "") << "{\"quantity\":\"" << c.quantity << "\",\"ok\":" << (c.ok ? "true" : "false") << ",\"change\":" << c.change
          << ",\"tolerance\":" << c.tolerance << ",\"blocks\":[";
        for (size_t b = 0; b < c.blocks.size(); ++b) j << (b ? "," : "") << c.blocks[b];
        j << "]}";
      }
      j << "]}";
      d->eq_checks = j.str();
    }
    if (report && cap > 0) {
      std::string t;
      for (const auto& n : rep.notes) t += n + "\n";
      char b[256];
      t += "stage                          ps     T/K     P/atm   ρ/g·cm⁻³\n";
      for (const auto& st : rep.stages) {
        std::snprintf(b, sizeof b, "%-28s %6.1f %7.1f %9.0f %9.4f\n", st.label.c_str(), st.ps, st.temperature, st.pressure, st.density);
        t += b;
      }
      for (const auto& c : rep.checks) {
        std::snprintf(b, sizeof b, "check %s: %s (change %.3g, tolerance %.3g)\n", c.quantity.c_str(), c.ok ? "ok" : "not yet", c.change, c.tolerance);
        t += b;
      }
      std::strncpy(report, t.c_str(), size_t(cap) - 1);
      report[cap - 1] = 0;
    }
    return e.until_converged && !rep.converged ? 1 : 0;
  });
}

int32_t caps_equilibrate_checks(caps_doc* d, char* json, int32_t cap) {
  return guard([&] {
    const int32_t need = int32_t(d->eq_checks.size() + 1);
    if (json && cap > 0) {
      const size_t m = std::min<size_t>(size_t(cap - 1), d->eq_checks.size());
      std::memcpy(json, d->eq_checks.data(), m);
      json[m] = 0;
    }
    return need;
  });
}

int32_t caps_internal_distances(caps_doc* d, int32_t* n, double* ratio, int32_t cap, int32_t* chains, double* b2) {
  return guard([&] {
    const caps::InternalDistances r = caps::internal_distances(d->frame);
    if (chains) *chains = r.chains;
    if (b2) *b2 = r.b2;
    const int32_t m = std::min<int32_t>(cap, int32_t(r.n.size()));
    for (int32_t k = 0; k < m; ++k) {
      if (n) n[k] = r.n[k];
      if (ratio) ratio[k] = r.ratio[k];
    }
    return m;
  });
}

caps_doc* caps_pack(const char* text, const char* base_dir, int32_t threads, caps_pack_progress_fn progress, void* user, char* report, int32_t cap) {
  caps::PackReport rep;
  auto write_report = [&]() {
    if (!report || cap <= 0) return;
    std::string t;
    for (const auto& n : rep.notes) t += n + "\n";
    std::strncpy(report, t.c_str(), size_t(cap) - 1);
    report[cap - 1] = 0;
  };
  try {
    caps::PackOptions o;
    std::string output;
    auto items = caps::parse_packmol_input(text ? text : "", base_dir ? base_dir : ".", o, &output);
    o.threads = threads;
    if (progress)
      o.progress = [&](const caps::PackProgress& p) { return progress(p.loop, p.loops, p.penalty, p.bad, user) == 0; };
    caps::System s = caps::pack(items, o, &rep);
    auto* d = new caps_doc;
    d->traj.topology = s;
    std::vector<caps::Vec3> p;
    for (const auto& a : s.atoms) p.push_back(a.pos);
    d->traj.positions.push_back(std::move(p));
    d->traj.cells.push_back(s.cell);
    d->traj.timesteps.push_back(0);
    d->traj.topology.notes = rep.notes;
    refresh(d);
    write_report();
    return d;
  } catch (const std::exception& e) {
    g_error = e.what();
  } catch (...) {
    g_error = "unknown error";
  }
  write_report();
  return nullptr;
}

int32_t caps_reaction_template(const char* name, char* out, int32_t cap) {
  return guard([&] {
    std::string t;
    if (!name || !*name)
      for (const auto& n : caps::builtin_template_names()) t += n + "\n";
    else
      t = caps::builtin_template(name);
    if (out && cap > 0) {
      std::strncpy(out, t.c_str(), size_t(cap) - 1);
      out[cap - 1] = 0;
    }
    return int32_t(t.size());
  });
}

int32_t caps_react(caps_doc* d, const char* templates, const caps_react_opts* o, caps_react_progress_fn progress, void* user, char* report,
                   int32_t cap) {
  return guard([&] {
    caps::ReactOptions r;
    r.templates = caps::parse_templates(templates ? templates : "");
    if (o->capture > 0)
      for (auto& t : r.templates) t.capture = o->capture;
    r.seed = o->seed;
    if (o->max_cycles > 0) r.max_cycles = o->max_cycles;
    if (o->max_per_cycle > 0) r.max_per_cycle = o->max_per_cycle;
    if (o->target_conversion > 0) r.target_conversion = o->target_conversion;
    r.relax = o->relax != 0;
    if (o->relax_iterations > 0) r.relax_iterations = o->relax_iterations;
    r.md_ps = o->md_ps;
    if (o->temperature > 0) r.temperature = o->temperature;
    if (o->cutoff > 0) r.energy.cutoff = o->cutoff;
    r.energy.coulomb = o->coulomb != 0;
    r.energy = elec(r.energy);
    if (progress)
      r.progress = [&](const caps::CycleRow& c) {
        caps_react_cycle row{c.cycle, c.reactions, c.total, c.clusters.clusters, c.atoms, c.conversion, c.clusters.largest_fraction,
                             c.clusters.reduced_mw, c.energy};
        return progress(&row, user) == 0;
      };
    caps::System s = d->traj.frame(d->current);
    if (!s.unwrapped) caps::make_molecules_whole(s);
    // the atom count changes when atoms leave, so the record is one document per state: keep the start and the end
    std::vector<caps::System> frames{s};
    r.frame = [&](const caps::System& x, int) { frames.push_back(x); };
    caps::ReactReport rep;
    caps::react(s, r, &rep);
    caps::Trajectory out;
    const bool same = frames.front().atoms.size() == s.atoms.size();
    out.topology = s;
    out.topology.notes = rep.notes;
    // frames with the final topology: all cycles when no atoms left, else the final structure only
    for (const auto& f : frames) {
      if (f.atoms.size() != s.atoms.size()) continue;
      if (!same && &f != &frames.back()) continue;
      std::vector<caps::Vec3> p;
      for (const auto& a : f.atoms) p.push_back(a.pos);
      out.positions.push_back(std::move(p));
      out.cells.push_back(f.cell);
      out.timesteps.push_back(int64_t(out.timesteps.size()));
    }
    if (out.positions.empty()) {
      std::vector<caps::Vec3> p;
      for (const auto& a : s.atoms) p.push_back(a.pos);
      out.positions.push_back(std::move(p));
      out.cells.push_back(s.cell);
      out.timesteps.push_back(0);
    }
    d->field.reset();   // new topology: the Field assignment no longer applies
    d->traj = std::move(out);
    d->current = d->traj.frames() - 1;
    refresh(d);
    if (report && cap > 0) {
      std::string t;
      for (const auto& n : rep.notes) t += n + "\n";
      std::strncpy(report, t.c_str(), size_t(cap) - 1);
      report[cap - 1] = 0;
    }
    return 0;
  });
}

int32_t caps_field_assign(caps_doc* d, const char* ff_path, const char* rules_path, int32_t charges) {
  return guard([&] {
    auto F = std::make_unique<FieldState>();
    F->ff_path = ff_path ? ff_path : "";
    F->base = caps::is_uff(F->ff_path) ? caps::uff_definition() : caps::load_forcefield(F->ff_path);
    if (rules_path && *rules_path && !caps::is_uff(F->ff_path)) {
      F->base.typing.clear();
      caps::load_typing(F->base, rules_path);
    }
    F->charges = charges == 1 ? "gasteiger" : charges == 2 ? "keep" : "types";
    // keep the file's types to restore them on clear (and the previous assignment's, if any)
    if (d->field) {
      F->file_types = d->field->file_types;
      F->file_type_table = d->field->file_type_table;
      F->file_charges = d->field->file_charges;
      F->file_has_charges = d->field->file_has_charges;
      F->extra = d->field->extra.name.empty() && d->field->base.name == F->base.name ? d->field->extra : caps::FFDef{};
    } else {
      for (const auto& a : d->traj.topology.atoms) {
        F->file_types.push_back({a.type, a.name});
        F->file_charges.push_back(a.charge);
      }
      F->file_type_table = d->traj.topology.types;
      F->file_has_charges = d->traj.topology.has_charges;
    }
    if (F->charges == "keep" && !d->traj.topology.has_charges && !F->file_has_charges)
      throw caps::FFError("the structure has no charges to keep; use the force field's or Gasteiger charges");
    // "keep" means the file's charges, not the previous assignment's
    for (size_t i = 0; i < d->traj.topology.atoms.size() && i < F->file_charges.size(); ++i) d->traj.topology.atoms[i].charge = F->file_charges[i];
    d->field = std::move(F);
    refresh(d);
    field_run(d);
    return d->field->complete ? 0 : 1;
  });
}

int32_t caps_field_report(caps_doc* d, char* json, int32_t cap) {
  if (!d->field) return report_out("", json, cap);
  return report_out(d->field->report, json, cap);
}

int32_t caps_field_override(caps_doc* d, int32_t index, const char* type) {
  return guard([&] {
    if (!d->field) throw caps::FFError("assign a force field first");
    if (index < 0 || size_t(index) >= d->frame.atoms.size()) throw caps::FFError("atom " + std::to_string(index + 1) + " out of range");
    const std::string t = type ? type : "";
    if (t.empty()) d->field->overrides.erase(index);
    else {
      bool known = d->field->base.type(t) != nullptr || d->field->extra.type(t) != nullptr;
      if (!known) throw caps::FFError("type " + t + " is not in " + d->field->base.name);
      d->field->overrides[index] = t;
    }
    field_run(d);
    return d->field->complete ? 0 : 1;
  });
}

int32_t caps_field_add_rule(caps_doc* d, const char* kind, const char* types, const char* style, const char* params) {
  return guard([&] {
    if (!d->field) throw caps::FFError("assign a force field first");
    const std::string k = kind ? kind : "";
    std::vector<caps::FFRule>* dst = k == "pair" ? &d->field->extra.pairs : k == "bond" ? &d->field->extra.bonds
                                   : k == "angle" ? &d->field->extra.angles : k == "dihedral" ? &d->field->extra.dihedrals
                                   : k == "improper" ? &d->field->extra.impropers : nullptr;
    if (!dst) throw caps::FFError("kind must be pair, bond, angle, dihedral or improper");
    const size_t need = k == "pair" ? 1 : k == "bond" ? 2 : k == "angle" ? 3 : 4;
    caps::FFRule r;
    // names the rules can match: the types and their equivalence names
    std::set<std::string> known;
    for (const caps::FFDef* def : {&d->field->base, &d->field->extra})
      for (const auto& t : def->types) {
        known.insert(t.name);
        for (const auto& [kind, e] : t.equiv) known.insert(e);
      }
    std::istringstream ts(types ? types : "");
    std::string names;
    for (std::string w; ts >> w;) {
      if (w != "*" && !known.count(w)) throw caps::FFError(w + " is neither a type of " + d->field->base.name + " nor an equivalence name");
      r.match.push_back(w == "*" ? w : caps::glob_escape(w));
      names += " " + w;
    }
    if (r.match.size() != need) throw caps::FFError(k + " parameters need " + std::to_string(need) + " atom types");
    std::istringstream ps(params ? params : "");
    for (std::string w; ps >> w;) {
      char* end = nullptr;
      const double v = std::strtod(w.c_str(), &end);
      if (!end || *end) throw caps::FFError("not a number: " + w);
      r.params.push_back(v);
    }
    if (r.params.empty()) throw caps::FFError("give the parameters");
    r.style = style ? style : "";
    r.name = "user: " + k + names;
    r.comment = "entered by hand in CAPS Studio: estimated, not from the published force field";
    dst->push_back(r);
    try {
      field_run(d);
    } catch (...) {   // a rule the force field cannot use (wrong style or parameter count) is taken back
      dst->pop_back();
      field_run(d);
      throw;
    }
    return d->field->complete ? 0 : 1;
  });
}

int32_t caps_field_remove_rules(caps_doc* d) {
  return guard([&] {
    if (!d->field) throw caps::FFError("assign a force field first");
    d->field->extra = caps::FFDef{};
    d->field->imported.clear();
    field_run(d);
    return d->field->complete ? 0 : 1;
  });
}

int32_t caps_field_import(caps_doc* d, const char* path) {
  return guard([&] {
    if (!d->field) throw caps::FFError("assign a force field first");
    const std::string p = path ? path : "";
    const bool lt = p.size() > 3 && p.compare(p.size() - 3, 3, ".lt") == 0;
    caps::FFDef imp = lt ? caps::import_moltemplate(p) : caps::load_forcefield(p);
    // parameters only: types the force field already has keep their definitions, typing rules stay the library's
    std::vector<caps::FFType> fresh;
    for (auto& t : imp.types)
      if (!d->field->base.type(t.name)) fresh.push_back(t);
    imp.types = fresh;
    imp.typing.clear();
    const auto slash = p.find_last_of("/\\");
    const std::string file = slash == std::string::npos ? p : p.substr(slash + 1);
    for (auto* v : {&imp.pairs, &imp.bonds, &imp.angles, &imp.dihedrals, &imp.impropers, &imp.bond_increments})
      for (auto& r : *v) {
        r.name = "imported: " + r.name;
        r.comment = "imported from " + file + (r.comment.empty() ? "" : "; " + r.comment);
      }
    caps::merge_forcefield(d->field->extra, imp);
    d->field->imported.push_back(p);
    field_run(d);
    return d->field->complete ? 0 : 1;
  });
}

int32_t caps_field_clear(caps_doc* d) {
  return guard([&] {
    if (!d->field) return 0;
    auto& atoms = d->traj.topology.atoms;
    for (size_t i = 0; i < atoms.size() && i < d->field->file_types.size(); ++i) {
      atoms[i].type = d->field->file_types[i].first;
      atoms[i].name = d->field->file_types[i].second;
      atoms[i].charge = d->field->file_charges[i];
    }
    d->traj.topology.types = d->field->file_type_table;
    d->traj.topology.has_charges = d->field->file_has_charges;
    d->field.reset();
    refresh(d);
    return 0;
  });
}

int32_t caps_field_types_file(caps_doc* d, const char* path) {
  return guard([&] {
    if (!d->field) throw caps::FFError("assign a force field first");
    std::ofstream f(path);
    if (!f) throw caps::FFError(std::string("cannot write ") + path);
    f << "# " << d->field->base.name << " types from CAPS Field (" << d->field->overrides.size() << " set by hand); caps ff apply --types\n";
    for (const auto& t : d->field->types) f << (t.empty() ? "?" : t) << "\n";
    return 0;
  });
}

int32_t caps_field_info(caps_doc* d, char* text, int32_t cap) {
  return guard([&] {
    const auto assigned = field_for_run(d);
    const caps::ForceField ff = assigned ? *assigned : caps::default_forcefield(d->frame);
    caps::Evaluator ev(ff, elec());
    std::vector<double> x, f;
    for (const auto& a : d->frame.atoms) x.insert(x.end(), a.pos.begin(), a.pos.end());
    const caps::EnergyTerms e = ev.compute(x, d->frame.cell, f);
    std::string t = ff.name + "\n";
    for (const auto& n : ff.notes) t += n + "\n";
    char b[400];
    std::snprintf(b, sizeof b,
                  "energy (kcal/mol): bond %.1f, angle %.1f, dihedral %.1f, improper %.1f, vdW %.1f, Coulomb %.1f, total %.1f\n"
                  "largest force %.3f kcal/mol/Å",
                  e.bond, e.angle, e.dihedral, e.improper, e.vdw, e.coulomb, e.total(), caps::max_force(f));
    t += b;
    if (d->frame.cell.valid()) {
      std::snprintf(b, sizeof b, ", pressure %.0f atm (0 K virial)", caps::pressure_atm(e.virial, d->frame.cell.volume()));
      t += b;
    }
    t += "\n";
    if (text && cap > 0) {
      std::strncpy(text, t.c_str(), size_t(cap) - 1);
      text[cap - 1] = 0;
    }
    return 0;
  });
}

int32_t caps_summary_get(caps_doc* d, caps_summary* o) {
  return guard([&] {
    const auto& s = d->frame;
    std::memset(o, 0, sizeof *o);
    int nm = 0;
    s.molecules(&nm);
    o->atoms = int64_t(s.atoms.size());
    o->bonds = int64_t(s.bonds.size());
    o->molecules = nm;
    o->frames = int64_t(d->traj.frames());
    o->bonds_from_file = s.bonds_from_file;
    o->has_charges = s.has_charges;
    o->cell_valid = s.cell.valid();
    o->unwrapped = s.unwrapped;
    o->cell_a = caps::norm(s.cell.a); o->cell_b = caps::norm(s.cell.b); o->cell_c = caps::norm(s.cell.c);
    o->volume = s.cell.volume();
    o->density = s.density();
    o->total_mass = s.total_mass();
    for (const auto& a : s.atoms) o->total_charge += a.charge;
    std::strncpy(o->format, s.source_format.c_str(), sizeof o->format - 1);
    return 0;
  });
}

int32_t caps_set_frame(caps_doc* d, int64_t f) {
  return guard([&] {
    if (f < 0 || size_t(f) >= d->traj.frames()) throw std::out_of_range("frame out of range");
    d->current = size_t(f);
    refresh(d);
    return 0;
  });
}

int32_t caps_set_wrap(caps_doc* d, int32_t w) {
  return guard([&] {
    d->wrap = w != 0;
    refresh(d);
    return 0;
  });
}

int32_t caps_atom(caps_doc* d, int32_t i, caps_atom_info* o) {
  return guard([&] {
    if (i < 0 || size_t(i) >= d->frame.atoms.size()) throw std::out_of_range("atom index out of range");
    const auto& a = d->frame.atoms[size_t(i)];
    std::memset(o, 0, sizeof *o);
    o->id = a.id; o->mol = a.mol; o->type = a.type; o->element = a.element; o->index = i;
    o->charge = a.charge; o->x = a.pos[0]; o->y = a.pos[1]; o->z = a.pos[2];
    std::strncpy(o->element_symbol, caps::element(a.element).symbol, 3);
    std::strncpy(o->name, a.name.c_str(), 15);
    return 0;
  });
}

int32_t caps_note_count(caps_doc* d) { return int32_t(d->frame.notes.size()); }
const char* caps_note(caps_doc* d, int32_t k) { return (k >= 0 && size_t(k) < d->frame.notes.size()) ? d->frame.notes[size_t(k)].c_str() : ""; }

int32_t caps_render(caps_doc* d, const caps_camera* cam, const caps_render_opts* opt, uint8_t* rgba) {
  return guard([&] {
    const auto img = d->renderer.render(d->frame, cam_of(cam), opts_of(d, opt));
    std::memcpy(rgba, img.rgba.data(), img.rgba.size());
    return 0;
  });
}

int32_t caps_pick(caps_doc* d, int32_t x, int32_t y) { return d->renderer.pick(x, y); }

int32_t caps_export_png(caps_doc* d, const caps_camera* cam, const caps_render_opts* opt, const char* path) {
  return guard([&] {
    caps::Renderer r;   // separate renderer so the view's pick buffer is not replaced
    caps::write_png(r.render(d->frame, cam_of(cam), opts_of(d, opt)), path);
    return 0;
  });
}

int32_t caps_export_svg(caps_doc* d, const caps_camera* cam, const caps_render_opts* opt, const char* path) {
  return guard([&] {
    std::ofstream f(path);
    if (!f) throw std::runtime_error(std::string("cannot write ") + path);
    f << caps::render_svg(d->frame, cam_of(cam), opts_of(d, opt));
    return 0;
  });
}

int32_t caps_measure(caps_doc* d, const int32_t* idx, int32_t n, double* value) {
  return guard([&] {
    std::vector<uint32_t> v;
    for (int32_t k = 0; k < n; ++k) {
      if (idx[k] < 0) throw std::out_of_range("atom index out of range");
      v.push_back(uint32_t(idx[k]));
    }
    *value = caps::measure(d->frame, v);
    return 0;
  });
}

int32_t caps_analyze(caps_doc* d, const char* props, const caps_analyze_opts* p, caps_analyze_progress_fn progress, void* user) {
  return caps_analyze_ex(d, props, p, nullptr, progress, user);
}

int32_t caps_analyze_ex(caps_doc* d, const char* props, const caps_analyze_opts* p, const caps_mech_opts* m, caps_analyze_progress_fn progress,
                        void* user) {
  return guard([&] {
    std::vector<std::string> ids, protocols;
    {
      std::string cur;
      for (const char* c = props ? props : ""; ; ++c) {
        if (*c == ',' || *c == 0) {
          if (!cur.empty()) (cur == "cij_strain" || cur == "cij_run" || cur == "tensile" || cur == "tg" ? protocols : ids).push_back(cur);
          cur.clear();
          if (*c == 0) break;
        } else if (*c != ' ') cur += *c;
      }
    }
    if (ids.empty() && protocols.empty()) throw std::invalid_argument("no properties requested");
    caps::AnalyzeOptions o;
    if (p) {
      o.first = std::max<int64_t>(0, p->first);
      o.last = p->last > 0 ? p->last : -1;
      if (p->stride > 0) o.stride = p->stride;
      if (p->frame_ps > 0) o.frame_ps = p->frame_ps;
      if (p->timestep_fs > 0) o.timestep_fs = p->timestep_fs;
      if (p->blocks > 0) o.blocks = p->blocks;
      o.elem_a = std::max(0, p->elem_a);
      o.elem_b = std::max(0, p->elem_b);
      o.inter_only = p->inter_only != 0;
      if (p->rdf_rmax > 0) o.rdf_rmax = p->rdf_rmax;
      if (p->rdf_dr > 0) o.rdf_dr = p->rdf_dr;
      if (p->qmax > 0) o.qmax = p->qmax;
      if (p->dq > 0) o.dq = p->dq;
      if (p->q_direct > 0) o.q_direct = p->q_direct;
      if (p->fit_from > 0) o.fit_from = p->fit_from;
      if (p->fit_to > 0) o.fit_to = p->fit_to;
      if (p->probe > 0) o.probe = p->probe;
      if (p->grid > 0) o.grid = p->grid;
      if (p->cutoff > 0) o.energy.cutoff = p->cutoff;
      o.energy = elec(o.energy);
      if (p->threads > 0) o.threads = p->threads;
    }
    caps_mech_opts mo{};
    if (m) mo = *m;
    if (mo.temperature > 0) o.temperature = mo.temperature;
    if (progress) o.progress = [&](const std::string& what, double f) { return progress(what.c_str(), f, user) != 0; };
    auto has = [&](const char* k) { return std::find(ids.begin(), ids.end(), k) != ids.end(); };
    // force field: the Field assignment (must be complete), else GAFF of C and H
    std::shared_ptr<const caps::ForceField> ff = field_for_run(d);
    std::vector<std::string> extra_notes;
    if (!ff && (has("ced") || has("delta") || has("cij_fluct") || has("adhesion") || !protocols.empty())) {
      caps::System s0 = d->traj.frame(0);
      if (!s0.unwrapped) caps::make_molecules_whole(s0);
      ff = std::make_shared<caps::ForceField>(default_ff(s0));
      extra_notes.push_back("force field: " + ff->name + " (none assigned in Field)");
    }
    if (ff) o.ff = ff.get();
    std::vector<caps::Property> res = ids.empty() ? std::vector<caps::Property>{} : caps::analyze(d->traj, ids, o);
    // protocols run on a copy of the current frame; the document is not changed
    auto frame_copy = [&] {
      caps::System s = d->traj.frame(d->current);
      if (!s.unwrapped) caps::make_molecules_whole(s);
      return s;
    };
    auto cancelled = [&](const std::string& w, double f) { return progress && progress(w.c_str(), f, user) == 0; };
    for (const auto& id : protocols) {
      if (id == "cij_strain") {
        const auto fr = caps::analysis_frames(d->traj, o);
        const size_t nc = std::max<size_t>(1, std::min<size_t>(fr.size(), mo.configurations > 0 ? size_t(mo.configurations) : 1));
        std::vector<caps::System> cs;
        if (fr.size() <= 1 || nc == 1) cs.push_back(frame_copy());
        else
          for (size_t k = 0; k < nc; ++k) {
            caps::System s = d->traj.frame(fr[k * (fr.size() - 1) / (nc - 1)]);
            if (!s.unwrapped) caps::make_molecules_whole(s);
            cs.push_back(std::move(s));
          }
        caps::StaticElasticOptions so;
        so.field = ff;
        so.energy = o.energy;
        if (mo.strain > 0) so.strain = mo.strain;
        so.progress = [&](const std::string& w, double f) { return !cancelled(w, f); };
        for (auto& q : caps::elastic_properties(caps::static_elastic(cs, so), "")) res.push_back(std::move(q));
      } else if (id == "cij_run") {
        caps::System s = frame_copy();
        caps::FluctuationRunOptions fo;
        fo.field = ff;
        fo.energy = o.energy;
        if (mo.temperature > 0) fo.temperature = mo.temperature;
        if (mo.run_ps > 0) fo.ps = mo.run_ps;
        if (mo.dt > 0) fo.dt = mo.dt;
        if (mo.seed) fo.seed = mo.seed;
        fo.progress = [&](const std::string& w, double f) { return !cancelled(w, f); };
        for (auto& q : caps::elastic_properties(caps::fluctuation_run(s, fo), "_fluct")) res.push_back(std::move(q));
      } else if (id == "tensile") {
        caps::System s = frame_copy();
        caps::TensileOptions to;
        to.field = ff;
        to.energy = o.energy;
        to.axis = std::clamp(mo.axis, 0, 2);
        if (mo.rate > 0) to.rate = mo.rate;
        if (mo.max_strain > 0) to.max_strain = mo.max_strain;
        if (mo.temperature > 0) to.temperature = mo.temperature;
        if (mo.dt > 0) to.dt = mo.dt;
        if (mo.pressure > 0) to.pressure = mo.pressure;
        if (mo.fit_strain > 0) to.fit_strain = mo.fit_strain;
        if (mo.seed) to.seed = mo.seed;
        to.lateral_pressure = mo.lateral_fixed == 0;
        if (mo.equilibrate_ps != 0) to.equilibrate_ps = std::max(0.0, mo.equilibrate_ps);
        to.new_velocities = s.velocities.size() != s.atoms.size();
        to.progress = [&](const caps::TensilePoint& q) {
          char b[96];
          std::snprintf(b, sizeof b, "tensile: strain %.3f, stress %.1f MPa", q.strain, q.stress);
          return !cancelled(b, q.strain / to.max_strain);
        };
        for (auto& q : caps::tensile_properties(caps::run_tensile(s, to))) res.push_back(std::move(q));
      } else if (id == "tg") {
        caps::System s = frame_copy();
        caps::CoolingOptions co;
        co.field = ff;
        co.energy = o.energy;
        if (mo.t_start > 0) co.t_start = mo.t_start;
        if (mo.t_end > 0) co.t_end = mo.t_end;
        if (mo.t_step > 0) co.t_step = mo.t_step;
        if (mo.ps_per_step > 0) co.ps_per_step = mo.ps_per_step;
        if (mo.dt > 0) co.dt = mo.dt;
        if (mo.pressure > 0) co.pressure = mo.pressure;
        if (mo.seed) co.seed = mo.seed;
        co.new_velocities = s.velocities.size() != s.atoms.size();
        if (mo.equilibrate_ps != 0) co.equilibrate_ps = std::max(0.0, mo.equilibrate_ps);
        co.progress = [&](const caps::ThermoRow& r, int k, int n) {
          char b[128];
          std::snprintf(b, sizeof b, "cooling: %.0f K (%d of %d), density %.4f g/cm³", r.target_temperature, k + 1, n, r.density);
          if (k < 0) std::snprintf(b, sizeof b, "cooling: equilibrating at %.0f K, density %.4f g/cm³", r.target_temperature, r.density);
          return !cancelled(b, std::max(0, k) / double(n));
        };
        for (auto& q : caps::cooling_properties(caps::run_cooling(s, co))) res.push_back(std::move(q));
      }
    }
    for (auto& r : res)
      if (r.id == "ced" || r.id == "delta" || r.id.rfind("cij", 0) == 0 || r.id == "tensile_modulus" || r.id == "tg")
        for (const auto& n : extra_notes) r.notes.push_back(n);
    const auto fr = caps::analysis_frames(d->traj, o);
    d->analysis = "{\"frames\":" + std::to_string(fr.size()) + ",\"of\":" + std::to_string(d->traj.frames()) + ",\"atoms\":" +
                  std::to_string(d->traj.topology.atoms.size()) + ",\"properties\":" + caps::properties_json(res) + "}";
    return 0;
  });
}

int32_t caps_analyze_report(caps_doc* d, char* json, int32_t cap) {
  return guard([&] {
    const int32_t need = int32_t(d->analysis.size() + 1);
    if (json && cap > 0) {
      const size_t m = std::min<size_t>(size_t(cap - 1), d->analysis.size());
      std::memcpy(json, d->analysis.data(), m);
      json[m] = 0;
    }
    return need;
  });
}

int32_t caps_rdf(caps_doc* d, int32_t ea, int32_t eb, double rmax, double dr, int32_t inter, double* r, double* g, int32_t cap) {
  return guard([&] {
    if (!d->frame.cell.valid()) throw std::runtime_error("g(r) needs a periodic cell");
    const auto out = caps::rdf(d->frame, ea, eb, rmax, dr, inter != 0);
    const int32_t m = std::min<int32_t>(cap, int32_t(out.size()));
    for (int32_t k = 0; k < m; ++k) { r[k] = out[size_t(k)].first; g[k] = out[size_t(k)].second; }
    return m;
  });
}

int32_t caps_molecules(caps_doc* d, caps_molecule* out, int32_t cap) {
  return guard([&] {
    const auto m = caps::molecule_shapes(d->frame);
    const int32_t n = std::min<int32_t>(cap, int32_t(m.size()));
    for (int32_t k = 0; k < n; ++k) {
      const auto& x = m[size_t(k)];
      out[k] = {x.molecule + 1, x.atoms, x.mass, x.rg, x.kappa2, x.com[0], x.com[1], x.com[2]};
    }
    return int32_t(m.size());
  });
}

int32_t caps_property_range(caps_doc* d, double* lo, double* hi) {
  return guard([&] {
    if (d->dcom.empty()) { *lo = *hi = 0; return 0; }
    *lo = *std::min_element(d->dcom.begin(), d->dcom.end());
    *hi = *std::max_element(d->dcom.begin(), d->dcom.end());
    return 0;
  });
}

int32_t caps_neighbours(caps_doc* d, int32_t i, int32_t k, int32_t* idx, double* dist) {
  return guard([&] {
    const auto& s = d->frame;
    if (i < 0 || size_t(i) >= s.atoms.size()) throw std::out_of_range("atom index out of range");
    std::vector<std::pair<double, int32_t>> v;
    v.reserve(s.atoms.size());
    for (size_t j = 0; j < s.atoms.size(); ++j)
      if (j != size_t(i)) v.emplace_back(caps::norm(s.cell.minimum_image(s.atoms[j].pos - s.atoms[size_t(i)].pos)), int32_t(j));
    const size_t m = std::min<size_t>(size_t(std::max(0, k)), v.size());
    std::partial_sort(v.begin(), v.begin() + long(m), v.end());
    for (size_t q = 0; q < m; ++q) { idx[q] = v[q].second; dist[q] = v[q].first; }
    return int32_t(m);
  });
}

}  // extern "C"

namespace {
caps::Json mol_json(const caps::MolInfo& m) {
  caps::Json j = caps::Json::object();
  j["ok"] = m.problems.empty();
  j["formula"] = m.formula;
  j["mass"] = m.mass;
  j["atoms"] = double(m.atoms);
  j["heavy"] = double(m.heavy);
  j["bonds"] = double(m.bonds);
  j["rings"] = double(m.rings);
  j["stereocentres"] = double(m.stereocentres);
  j["stereo_bonds"] = double(m.stereo_bonds);
  j["charge"] = double(m.charge);
  caps::Json p = caps::Json::array();
  for (const auto& x : m.problems) p.push_back(x);
  j["problems"] = p;
  return j;
}
bool m_problems_empty(const caps::Json& j) { return j["problems"].size() == 0; }
}  // namespace

extern "C" int32_t caps_smiles_info(const char* smiles, char* json, int32_t cap) {
  caps::Json j = caps::Json::object();
  try {
    caps::MolGraph g = caps::parse_smiles(smiles ? smiles : "");
    caps::add_hydrogens(g);
    j = mol_json(caps::molecule_info(g));
    if (!m_problems_empty(j)) j["error"] = j["problems"][size_t(0)].str();
  } catch (const caps::SmilesError& e) {
    j["ok"] = false;
    j["error"] = std::string(e.what());
    j["position"] = double(e.position);
  } catch (const std::exception& e) {
    j["ok"] = false;
    j["error"] = std::string(e.what());
  }
  return report_out(j.dump(), json, cap);
}

extern "C" caps_doc* caps_build_smiles(const char* smiles, const char* ff_path, const caps_build_opts* o, char* report, int32_t cap) {
  try {
    caps::BuildOptions b;
    if (o) {
      b.conformers = o->conformers > 0 ? o->conformers : 1;
      b.seed = o->seed ? o->seed : 1;
    }
    b.forcefield = ff_path ? ff_path : "";
    const caps::BuildResult r = caps::build_molecule(smiles ? smiles : "", b);
    auto* d = new caps_doc;
    d->traj.topology = r.system;
    for (const auto& c : r.conformers) {
      d->traj.positions.push_back(c.pos);
      d->traj.cells.push_back(r.system.cell);
      d->traj.timesteps.push_back(int64_t(d->traj.timesteps.size()));
    }
    refresh(d);
    caps::Json j = mol_json(r.info);
    j["smiles"] = r.graph.smiles;
    j["method"] = r.method;
    caps::Json cs = caps::Json::array();
    const double e0 = r.conformers.front().energy;
    for (const auto& c : r.conformers) {
      caps::Json x = caps::Json::object();
      x["energy"] = c.energy;
      x["rel"] = c.energy - e0;
      x["minimised"] = c.minimised;
      cs.push_back(x);
    }
    j["conformers"] = cs;
    caps::Json notes = caps::Json::array();
    for (const auto& n : r.notes) notes.push_back(n);
    j["notes"] = notes;
    report_out(j.dump(), report, cap);
    return d;
  } catch (const std::exception& e) {
    g_error = e.what();
    return nullptr;
  }
}


extern "C" int32_t caps_smiles_depict(const char* smiles, char* json, int32_t cap) {
  caps::Json j = caps::Json::object();
  try {
    caps::MolGraph g = caps::parse_smiles(smiles ? smiles : "");
    caps::MolGraph gh = g;
    caps::add_hydrogens(gh);
    j = mol_json(caps::molecule_info(gh));
    if (!m_problems_empty(j)) j["error"] = j["problems"][size_t(0)].str();
    const auto p = caps::depict(g);
    caps::Json atoms = caps::Json::array();
    for (size_t i = 0; i < g.atoms.size(); ++i) {
      const auto& a = g.atoms[i];
      caps::Json x = caps::Json::object();
      x["z"] = double(a.element);
      x["symbol"] = std::string(a.element ? caps::element(a.element).symbol : "*");
      x["x"] = p[i][0];
      x["y"] = p[i][1];
      x["h"] = double(std::max(0, gh.atoms[i].hcount));
      x["charge"] = double(a.charge);
      x["isotope"] = double(a.isotope);
      x["hcount"] = double(a.hcount);
      x["aromatic"] = a.aromatic;
      x["bracket"] = a.bracket;
      x["chiral"] = double(a.chiral);
      x["map"] = double(a.map);
      caps::Json o = caps::Json::array();
      for (int v : a.order) o.push_back(double(v));
      x["order"] = o;
      atoms.push_back(x);
    }
    caps::Json bonds = caps::Json::array();
    for (const auto& b : g.bonds) {
      caps::Json x = caps::Json::object();
      x["a"] = double(b.a);
      x["b"] = double(b.b);
      x["order"] = double(b.order);
      x["dir"] = double(b.dir);
      bonds.push_back(x);
    }
    j["atoms"] = atoms;
    j["bonds"] = bonds;
  } catch (const caps::SmilesError& e) {
    j["ok"] = false;
    j["error"] = std::string(e.what());
    j["position"] = double(e.position);
  } catch (const std::exception& e) {
    j["ok"] = false;
    j["error"] = std::string(e.what());
  }
  return report_out(j.dump(), json, cap);
}

extern "C" int32_t caps_smiles_write(const char* graph_json, char* smiles, int32_t cap) {
  return guard([&] {
    const caps::Json j = caps::Json::parse(graph_json ? graph_json : "");
    caps::MolGraph g;
    for (const auto& x : j["atoms"].items()) {
      caps::MolAtom a;
      a.element = int(x.num("z", 6));
      a.charge = int(x.num("charge", 0));
      a.isotope = int(x.num("isotope", 0));
      a.hcount = int(x.num("hcount", -1));
      a.aromatic = x.has("aromatic") && x["aromatic"].boolean();
      a.bracket = x.has("bracket") && x["bracket"].boolean();
      a.chiral = int(x.num("chiral", 0));
      a.map = int(x.num("map", 0));
      if (x.has("order"))
        for (const auto& v : x["order"].items()) a.order.push_back(int(v.number()));
      g.atoms.push_back(a);
    }
    const int n = int(g.atoms.size());
    for (const auto& x : j["bonds"].items()) {
      caps::MolBond b;
      b.a = int(x.num("a", -1));
      b.b = int(x.num("b", -1));
      b.order = int(x.num("order", 1));
      b.dir = int(x.num("dir", 0));
      if (b.a < 0 || b.b < 0 || b.a >= n || b.b >= n || b.a == b.b) throw std::runtime_error("a bond refers to a missing atom");
      g.bonds.push_back(b);
    }
    g.heavy = n;
    return report_out(caps::write_smiles(g), smiles, cap);
  });
}

extern "C" void caps_set_palette(int32_t p) { caps::set_palette(p == 1 ? caps::Palette::OkabeIto : p == 2 ? caps::Palette::Monochrome : caps::Palette::Caps); }
extern "C" void caps_set_threads(int32_t n) { caps::set_max_threads(n); }

namespace {
caps::Json bench_json(const caps::BenchTable& t) {
  caps::Json j = caps::Json::object();
  j["id"] = t.id;
  j["title"] = t.title;
  j["scope"] = t.scope;
  caps::Json cols = caps::Json::array();
  for (const auto& c : t.columns) cols.push_back(c);
  j["columns"] = cols;
  caps::Json rows = caps::Json::array();
  for (const auto& r : t.rows) {
    caps::Json row = caps::Json::object(), cells = caps::Json::array();
    for (const auto& c : r.cells) cells.push_back(c);
    row["cells"] = cells;
    row["status"] = r.status;
    rows.push_back(row);
  }
  j["rows"] = rows;
  j["status"] = t.status;
  j["note"] = t.note;
  j["seconds"] = t.seconds;
  return j;
}

caps::BenchTable bench_from(const caps::Json& j) {
  caps::BenchTable t;
  t.id = j.text("id");
  t.title = j.text("title");
  t.scope = j.text("scope");
  t.status = j.text("status");
  t.note = j.text("note");
  t.seconds = j.num("seconds", 0);
  if (j.has("columns")) for (const auto& c : j["columns"].items()) t.columns.push_back(c.str());
  if (j.has("rows"))
    for (const auto& r : j["rows"].items()) {
      caps::BenchRow row;
      row.status = r.text("status");
      for (const auto& c : r["cells"].items()) row.cells.push_back(c.str());
      t.rows.push_back(row);
    }
  return t;
}
}  // namespace

extern "C" int32_t caps_bench_list(char* json, int32_t cap) {
  return guard([&] {
    caps::Json a = caps::Json::array();
    for (const auto& id : caps::bench_ids()) a.push_back(bench_json(caps::bench_describe(id)));
    return report_out(a.dump(), json, cap);
  });
}

extern "C" int32_t caps_bench_run(const char* id, const char* samples, const char* forcefields, int32_t repeats, int32_t quick,
                                  caps_bench_progress_fn progress, void* user, char* json, int32_t cap) {
  return guard([&] {
    caps::BenchOptions o;
    o.samples = samples ? samples : "";
    o.forcefields = forcefields ? forcefields : "";
    o.repeats = repeats > 0 ? repeats : 3;
    o.quick = quick != 0;
    if (progress) o.progress = [&](const std::string& t, const std::string& w, double f) { return progress(t.c_str(), w.c_str(), f, user) == 0; };
    return report_out(bench_json(caps::run_bench(id ? id : "", o)).dump(), json, cap);
  });
}

extern "C" int32_t caps_bench_write(const char* tables_json, const char* dir) {
  return guard([&] {
    const caps::Json a = caps::Json::parse(tables_json ? tables_json : "[]");
    std::vector<caps::BenchTable> ts;
    for (const auto& j : a.items()) ts.push_back(bench_from(j));
    const std::string d = dir ? dir : ".";
    std::filesystem::create_directories(d);
    std::ofstream(d + "/results.md") << caps::bench_markdown(ts);
    std::ofstream(d + "/results.tex") << caps::bench_latex(ts);
    for (const auto& t : ts)
      if (!t.rows.empty()) std::ofstream(d + "/" + t.id + ".csv") << caps::bench_csv(t);
    return 0;
  });
}

namespace {
caps::ChainSpec spec_from(const std::string& text) {
  const caps::Json j = caps::Json::parse(text);
  caps::ChainSpec c;
  if (j.has("units"))
    for (const auto& u : j["units"].items()) c.units.push_back({u.text("name"), u.text("smiles")});
  c.sequence = caps::sequence_from_string(j.text("sequence", "homopolymer"));
  c.dp = int(j.num("dp", 20));
  if (j.has("blocks")) for (const auto& b : j["blocks"].items()) c.blocks.push_back(int(b.number()));
  if (j.has("weights")) for (const auto& w : j["weights"].items()) c.weights.push_back(w.number());
  c.pattern = j.text("pattern");
  c.pm = j.num("pm", 0.5);
  c.forcefield = j.text("forcefield");
  const std::string tac = j.text("tacticity", "atactic");
  c.tacticity = caps::tacticity_from_string(tac);
  return c;
}
}  // namespace

extern "C" int32_t caps_unit_info(const char* smiles, char* json, int32_t cap) {
  caps::Json j = caps::Json::object();
  try {
    const caps::UnitInfo u = caps::repeat_unit_info(smiles ? smiles : "");
    j["ok"] = true;
    j["formula"] = u.formula;
    j["mass"] = u.mass;
    j["atoms"] = double(u.atoms);
    j["head_element"] = u.head_element;
    j["tail_element"] = u.tail_element;
    j["stereocentres"] = double(u.stereocentres);
  } catch (const std::exception& e) {
    j["ok"] = false;
    j["error"] = std::string(e.what());
  }
  return report_out(j.dump(), json, cap);
}

extern "C" int32_t caps_chain_preview(const char* spec_json, uint64_t seed, char* json, int32_t cap) {
  caps::Json j = caps::Json::object();
  try {
    const caps::ChainSpec c = spec_from(spec_json ? spec_json : "{}");
    if (c.units.empty()) throw std::runtime_error("no repeat unit");
    const auto seq = caps::chain_sequence(c, seed);
    const caps::MolGraph g = caps::chain_graph(c, seq);
    const caps::MolInfo m = caps::molecule_info(g);
    caps::Json s = caps::Json::array();
    for (int k : seq) s.push_back(double(k));
    j["ok"] = true;
    j["sequence"] = s;
    j["formula"] = m.formula;
    j["mass"] = m.mass;
    j["atoms"] = double(m.atoms);
    j["smiles"] = g.smiles;
  } catch (const std::exception& e) {
    j["ok"] = false;
    j["error"] = std::string(e.what());
  }
  return report_out(j.dump(), json, cap);
}

extern "C" caps_doc* caps_grow_chains(const char* spec_json, const caps_grow_opts* o, caps_progress_fn progress, void* user, char* report, int32_t cap) {
  try {
    caps::ChainSpec c = spec_from(spec_json ? spec_json : "{}");
    caps::GrowOptions g;
    g.chains = o->chains;
    if (o->dp > 0) c.dp = o->dp;
    c.tacticity = o->tacticity == 1 ? caps::Tacticity::Isotactic : o->tacticity == 2 ? caps::Tacticity::Syndiotactic : caps::Tacticity::Atactic;
    g.seed = o->seed;
    g.box = o->box;
    g.density = o->density;
    g.contact_scale = o->contact_scale > 0 ? o->contact_scale : o->contact_scale < 0 ? -o->contact_scale : 1.0;
    g.auto_scale = o->contact_scale < 0;   // negative: start there and step down when crowded
    g.curve = o->curve != 0;
    if (progress) g.progress = [&](int done, int total, int restarts) { return progress(done, total, restarts, user) == 0; };
    caps::GrowReport rep;
    caps::System s = caps::grow_chains(c, g, &rep);
    auto* d = new caps_doc;
    d->traj.topology = s;
    std::vector<caps::Vec3> p;
    for (const auto& a : s.atoms) p.push_back(a.pos);
    d->traj.positions.push_back(std::move(p));
    d->traj.cells.push_back(s.cell);
    d->traj.timesteps.push_back(0);
    refresh(d);
    std::string t;
    for (const auto& n : rep.notes) t += n + "\n";
    report_out(t, report, cap);
    return d;
  } catch (const std::exception& e) {
    g_error = e.what();
    return nullptr;
  }
}

namespace {
caps::SlabOptions slab_from(const caps::Json& j) {
  caps::SlabOptions o;
  o.h = int(j.num("h", 0)), o.k = int(j.num("k", 0)), o.l = int(j.num("l", 1));
  o.layers = int(j.num("layers", 3));
  o.termination = int(j.num("termination", 0));
  o.vacuum = j.num("vacuum", 15);
  o.orthogonal = j.num("orthogonal", 1) != 0;
  o.max_strain = j.num("max_strain", 0.02);
  o.na = int(j.num("na", 1)), o.nb = int(j.num("nb", 1));
  o.passivate = j.num("passivate", 0) != 0;
  return o;
}

caps_doc* doc_of(const caps::System& s) {
  auto* d = new caps_doc;
  d->traj.topology = s;
  std::vector<caps::Vec3> p;
  for (const auto& a : s.atoms) p.push_back(a.pos);
  d->traj.positions.push_back(std::move(p));
  d->traj.cells.push_back(s.cell);
  d->traj.timesteps.push_back(0);
  refresh(d);
  return d;
}
}  // namespace

extern "C" int32_t caps_surface_terminations(const char* cif_path, int32_t h, int32_t k, int32_t l, char* json, int32_t cap) {
  caps::Json j = caps::Json::object();
  try {
    const caps::System bulk = caps::read_cif(cif_path ? cif_path : "");
    double d = 0;
    const auto terms = caps::slab_terminations(bulk, h, k, l, &d);
    j["ok"] = true;
    j["d"] = d;
    std::vector<size_t> all(bulk.atoms.size());
    for (size_t i = 0; i < all.size(); ++i) all[i] = i;
    j["formula"] = caps::formula_of(bulk, all);
    j["atoms"] = double(bulk.atoms.size());
    j["density"] = bulk.density();
    const auto& c = bulk.cell;
    const double A = caps::norm(c.a), B = caps::norm(c.b), C = caps::norm(c.c);
    auto ang = [](const caps::Vec3& u, const caps::Vec3& v) { return std::acos(std::clamp(caps::dot(u, v) / (caps::norm(u) * caps::norm(v)), -1.0, 1.0)) * 57.29577951308232; };
    caps::Json cell = caps::Json::array();
    for (double x : {A, B, C, ang(c.b, c.c), ang(c.a, c.c), ang(c.a, c.b)}) cell.push_back(x);
    j["cell"] = cell;
    caps::Json notes = caps::Json::array();
    for (const auto& n : bulk.notes) notes.push_back(n);
    j["notes"] = notes;
    caps::Json arr = caps::Json::array();
    for (const auto& t : terms) {
      caps::Json x = caps::Json::object();
      x["label"] = t.label;
      x["top"] = t.top;
      x["bottom"] = t.bottom;
      x["gap"] = t.gap;
      x["bonds_per_nm2"] = t.bonds_per_nm2;
      arr.push_back(x);
    }
    j["terminations"] = arr;
  } catch (const std::exception& e) {
    j = caps::Json::object();
    j["ok"] = false;
    j["error"] = std::string(e.what());
  }
  return report_out(j.dump(), json, cap);
}

extern "C" caps_doc* caps_surface_build(const char* cif_path, const char* options_json, char* report, int32_t cap) {
  try {
    const caps::System bulk = caps::read_cif(cif_path ? cif_path : "");
    caps::SlabReport rep;
    const caps::System s = caps::cleave(bulk, slab_from(caps::Json::parse(options_json && *options_json ? options_json : "{}")), &rep);
    std::string t;
    for (const auto& n : rep.notes) t += n + "\n";
    report_out(t, report, cap);
    return doc_of(s);
  } catch (const std::exception& e) {
    g_error = e.what();
    return nullptr;
  }
}

extern "C" caps_doc* caps_interface_build(const char* options_json, const char* spec_json, const caps_grow_opts* o, caps_progress_fn progress, void* user, char* report,
                                          int32_t cap) {
  try {
    const caps::Json j = caps::Json::parse(options_json && *options_json ? options_json : "{}");
    caps::SlabOptions so = slab_from(j.has("slab") ? j["slab"] : caps::Json::object());
    so.vacuum = std::max(so.vacuum, 10.0);   // free surfaces; the interface sets the final cell
    caps::SlabReport sr;
    const caps::System slab = caps::cleave(caps::read_cif(j.text("crystal")), so, &sr);
    caps::ChainSpec c = spec_from(spec_json ? spec_json : "{}");
    caps::InterfaceOptions io;
    const caps::Json f = j.has("film") ? j["film"] : caps::Json::object();
    io.film = f.num("thickness", 30);
    io.density = f.num("density", 0.9);
    io.chains = int(f.num("chains", 0));
    io.gap = f.num("gap", 1.0);
    io.vacuum = f.num("vacuum", 0);
    if (o) {
      if (o->dp > 0) c.dp = o->dp;
      c.tacticity = o->tacticity == 1 ? caps::Tacticity::Isotactic : o->tacticity == 2 ? caps::Tacticity::Syndiotactic : caps::Tacticity::Atactic;
      io.grow.seed = o->seed;
      io.grow.contact_scale = o->contact_scale > 0 ? o->contact_scale : 1.0;
      io.grow.curve = o->curve != 0;
    }
    if (progress) io.grow.progress = [&](int done, int total, int restarts) { return progress(done, total, restarts, user) == 0; };
    caps::GrowReport rep;
    const caps::System s = caps::build_interface(slab, c, io, &rep);
    std::string t;
    for (const auto& n : sr.notes) t += n + "\n";
    for (const auto& n : rep.notes) t += n + "\n";
    report_out(t, report, cap);
    caps_doc* d = doc_of(s);
    d->held_mol = 1;
    return d;
  } catch (const std::exception& e) {
    g_error = e.what();
    return nullptr;
  }
}

namespace {
caps::System nano_from(const caps::Json& j, std::array<bool, 3>& keep, std::string& notes) {
  const std::string kind = j.text("kind", "tube");
  caps::NanoReport r;
  caps::System f;
  keep = {false, false, false};
  if (kind == "tube") {
    caps::NanotubeOptions t;
    t.n = int(j.num("n", 10)), t.m = int(j.num("m", 10));
    t.length = j.num("length", 25);
    t.periodic = j.num("periodic", 1) != 0;
    f = caps::nanotube(t, &r);
    keep = {false, false, t.periodic};
  } else if (kind == "sheet") {
    caps::SheetOptions sh;
    sh.lx = j.num("lx", 20), sh.ly = j.num("ly", 20);
    sh.layers = int(j.num("layers", 1));
    sh.periodic = j.num("periodic", 1) != 0;
    f = caps::graphene_sheet(sh, &r);
    keep = {sh.periodic, sh.periodic, false};
  } else if (kind == "particle") {
    caps::ParticleOptions po;
    po.shape = caps::particle_shape_from_string(j.text("shape", "sphere"));
    po.radius = j.num("radius", 12);
    po.on_atom = j.num("on_atom", 1) != 0;
    po.passivate = j.num("passivate", 0) != 0;
    f = caps::nanoparticle(caps::read_cif(j.text("crystal")), po, &r);
  } else {
    throw std::invalid_argument("kind must be tube, sheet or particle");
  }
  for (const auto& n : r.notes) notes += n + "\n";
  return f;
}
}  // namespace

extern "C" caps_doc* caps_nano_build(const char* options_json, char* report, int32_t cap) {
  try {
    std::array<bool, 3> keep;
    std::string notes;
    const caps::System f = nano_from(caps::Json::parse(options_json && *options_json ? options_json : "{}"), keep, notes);
    report_out(notes, report, cap);
    return doc_of(f);
  } catch (const std::exception& e) {
    g_error = e.what();
    return nullptr;
  }
}

extern "C" caps_doc* caps_nano_embed(const char* options_json, const char* spec_json, const caps_grow_opts* o, caps_progress_fn progress, void* user, char* report,
                                     int32_t cap) {
  try {
    const caps::Json j = caps::Json::parse(options_json && *options_json ? options_json : "{}");
    std::array<bool, 3> keep;
    std::string notes;
    const caps::System f = nano_from(j, keep, notes);
    caps::ChainSpec c = spec_from(spec_json ? spec_json : "{}");
    caps::FillerMatrixOptions fo;
    const caps::Json m = j.has("matrix") ? j["matrix"] : caps::Json::object();
    fo.chains = int(m.num("chains", 10));
    fo.density = m.num("density", 0.9);
    fo.keep_axis = keep;
    if (o) {
      if (o->dp > 0) c.dp = o->dp;
      c.tacticity = o->tacticity == 1 ? caps::Tacticity::Isotactic : o->tacticity == 2 ? caps::Tacticity::Syndiotactic : caps::Tacticity::Atactic;
      fo.grow.seed = o->seed;
      fo.grow.contact_scale = o->contact_scale > 0 ? o->contact_scale : 1.0;
      fo.grow.curve = o->curve != 0;
    }
    if (progress) fo.grow.progress = [&](int done, int total, int restarts) { return progress(done, total, restarts, user) == 0; };
    caps::FillerReport fr;
    const caps::System s = caps::embed_filler(f, c, fo, &fr);
    for (const auto& n : fr.notes) notes += n + "\n";
    report_out(notes, report, cap);
    caps_doc* d = doc_of(s);
    d->held_mol = 1;
    return d;
  } catch (const std::exception& e) {
    g_error = e.what();
    return nullptr;
  }
}

extern "C" int32_t caps_insert_molecules(caps_doc* d, const char* smiles, int32_t count, double tolerance, uint64_t seed, char* report, int32_t cap) {
  return guard([&] {
    caps::BuildOptions bo;
    bo.forcefield = "uff";
    const caps::BuildResult br = caps::build_molecule(smiles ? smiles : "", bo);
    caps::System guest = br.system;
    guest.title = smiles ? smiles : "molecule";
    caps::PackOptions po;
    if (tolerance > 0) po.tolerance = tolerance;
    po.seed = seed ? seed : 1;
    caps::PackReport pr;
    caps::System host = d->traj.frame(d->current);
    caps::System s = caps::insert_molecules(host, guest, std::max(1, count), po, &pr);
    d->field.reset();
    d->traj = caps::Trajectory{};
    d->traj.topology = s;
    std::vector<caps::Vec3> p;
    for (const auto& a : s.atoms) p.push_back(a.pos);
    d->traj.positions.push_back(std::move(p));
    d->traj.cells.push_back(s.cell);
    d->traj.timesteps.push_back(0);
    d->current = 0;
    refresh(d);
    std::string t;
    for (const auto& n : pr.notes) t += n + "\n";
    report_out(t, report, cap);
    return 0;
  });
}

extern "C" void caps_set_held_molecule(caps_doc* d, int64_t mol) {
  if (d) d->held_mol = std::max<int64_t>(0, mol);
}

extern "C" int64_t caps_held_molecule(const caps_doc* d) { return d ? d->held_mol : 0; }

extern "C" void caps_set_electrostatics(int32_t mode, double ewald_rtol, double pme_spacing, int32_t pme_order) {
  g_elec.mode = mode == 1 ? 1 : 0;
  if (ewald_rtol > 0) g_elec.rtol = ewald_rtol;
  if (pme_spacing > 0) g_elec.spacing = pme_spacing;
  if (pme_order >= 3 && pme_order <= 12) g_elec.order = pme_order;
}
