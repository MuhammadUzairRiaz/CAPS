#include "caps_c.h"

#include <cstdio>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <ctime>
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
#include "caps/checks.hpp"
#include "caps/pipeline.hpp"
#include "caps/bundle.hpp"
#include "caps/crystal.hpp"
#include "caps/spacegroup.hpp"
#include "caps/peptide.hpp"
#include "caps/solvate.hpp"
#include "caps/appearance.hpp"
#include "caps/trajectory.hpp"
#include "caps/torsion.hpp"
#include "caps/edit.hpp"
#include "caps/interactions.hpp"
#include "caps/import.hpp"
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

// Styles, colours, surfaces and polyhedra of the Studio view (caps_set_appearance), prepared for the current frame.
struct AppearanceState {
  bool active = false;
  std::vector<std::pair<std::string, int>> layers;   // expression ("" all) → caps::Style, applied in order
  int colour = -1;                 // −1 the view's; 0 element, 1 molecule, 2 type, 3 distance to own centre, 4 partial charge
  caps::Ramp ramp = caps::Ramp::BlueOrange;
  int surface = 0;                 // 0 none, 1 accessible, 2 van der Waals, 3 excluded
  double probe = 1.4, spacing = 0.6;
  float opacity = 0.6f;
  std::string surface_expr;        // atoms the surface wraps ("" all)
  int surface_colour = 1;          // 0 one colour, 1 electrostatic potential, 2 nearest atom
  // for the current frame
  std::vector<uint8_t> style;      // 255: the view's style
  std::unique_ptr<caps::Mesh> mesh, poly;
  std::vector<caps::Segment> ribbon;
  double phi_lo = 0, phi_hi = 0;
  std::string error;
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
  std::unique_ptr<caps::Pipeline> pipeline;        // caps_pipeline_set: steps run on every shown frame
  std::unique_ptr<caps::PipelineState> pstate;     // its result for the current frame
  std::vector<int32_t> shown_of;                   // frame index → first shown particle (−1: deleted)
  std::array<int, 3> cell_repeats{1, 1, 1};        // caps_crystal_build: a supercell of this many unit cells
  std::vector<caps::Segment> overlay;              // caps_peptide_build with ribbon: tubes drawn with the atoms (no pipeline)
  AppearanceState look;                            // caps_set_appearance
  int smooth_window = 1;                           // caps_set_smoothing: frames averaged for display
  std::vector<std::vector<caps::Vec3>> scan_frames; // caps_torsion_scan: the geometry of each point
  std::vector<caps::Vec3> scan_original;           // the frame before a scan point was shown
  std::vector<char> selection;                     // caps_select: the Studio's selection over the frame's atoms
  struct Snapshot { caps::System topology; std::vector<caps::Vec3> positions; std::string what; };
  std::vector<Snapshot> undo, redo;                // caps_edit history
  std::vector<caps::Segment> checks;               // caps_interactions: H-bonds, contacts, clashes drawn in the view
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

void run_doc_pipeline(caps_doc* d) {
  d->pstate.reset();
  d->shown_of.clear();
  if (!d->pipeline) return;
  const int64_t ts = d->current < d->traj.timesteps.size() ? d->traj.timesteps[d->current] : 0;
  d->pstate = std::make_unique<caps::PipelineState>(caps::run_pipeline(d->frame, *d->pipeline, int(d->current), ts, &d->traj));
  d->shown_of.assign(d->frame.atoms.size(), -1);
  for (size_t k = 0; k < d->pstate->origin.size(); ++k) {
    const int o = d->pstate->origin[k];
    if (o >= 0 && size_t(o) < d->shown_of.size() && d->shown_of[size_t(o)] < 0) d->shown_of[size_t(o)] = int32_t(k);
  }
}

// What the view draws: the pipeline's particles when there is a pipeline, else the frame.
const caps::System& shown(const caps_doc* d) { return d->pstate ? d->pstate->system : d->frame; }

void prepare_appearance(caps_doc* d) {
  AppearanceState& L = d->look;
  L.style.clear(), L.mesh.reset(), L.poly.reset(), L.ribbon.clear(), L.error.clear();
  L.phi_lo = L.phi_hi = 0;
  if (!L.active) return;
  const caps::System& f = d->frame;
  const size_t n = f.atoms.size();
  std::unique_ptr<caps::PipelineState> st;
  auto mask = [&](const std::string& expr) {
    std::vector<char> m(n, 1);
    if (expr.empty()) return m;
    if (!st) st = std::make_unique<caps::PipelineState>(caps::run_pipeline(f, caps::Pipeline{}));
    const auto v = caps::evaluate_expression(*st, expr);
    for (size_t i = 0; i < n && i < v.size(); ++i) m[i] = v[i] != 0;
    return m;
  };
  try {
    if (!L.layers.empty()) {
      L.style.assign(n, 255);
      for (const auto& [expr, style] : L.layers) {
        const auto m = mask(expr);
        for (size_t i = 0; i < n; ++i) if (m[i]) L.style[i] = uint8_t(style);
      }
      std::vector<char> centres(n, 0), ribbon(n, 0);
      bool any_poly = false, any_ribbon = false;
      for (size_t i = 0; i < n; ++i) {
        if (L.style[i] == uint8_t(caps::Style::Polyhedra)) centres[i] = 1, any_poly = true;
        if (L.style[i] == uint8_t(caps::Style::Ribbon)) ribbon[i] = 1, any_ribbon = true;
      }
      if (any_poly) L.poly = std::make_unique<caps::Mesh>(caps::polyhedra(f, centres));
      if (any_ribbon) {
        const auto mol = f.molecules();
        for (const auto& path : caps::ribbon_paths(f, ribbon)) {
          size_t nearest = 0;
          double best = 1e300;
          for (size_t i = 0; i < n; ++i) if (ribbon[i]) { const double dd = caps::norm(f.atoms[i].pos - path.front()); if (dd < best) best = dd, nearest = i; }
          const unsigned rgb = caps::molecule_colour(mol[nearest]);
          for (size_t k = 0; k + 1 < path.size(); ++k) {
            caps::Segment sg;
            sg.a = path[k], sg.b = path[k + 1], sg.rgb = rgb, sg.radius = 0.55;
            L.ribbon.push_back(sg);
          }
        }
      }
    }
    if (L.surface > 0) {
      caps::SurfaceOptions so;
      so.kind = L.surface == 2 ? caps::SurfaceKind::VanDerWaals : L.surface == 3 ? caps::SurfaceKind::Excluded : caps::SurfaceKind::Accessible;
      so.probe = L.probe;
      so.spacing = L.spacing;
      so.atoms = mask(L.surface_expr);
      L.mesh = std::make_unique<caps::Mesh>(caps::surface_mesh(f, so));
      auto& M = *L.mesh;
      if (L.surface_colour == 1) {
        const auto phi = caps::surface_potential(f, M);
        // the colour range: the 95th percentile of |φ|, so a few vertices beside formal charges do not wash out the map
        std::vector<double> a;
        a.reserve(phi.size());
        for (double x : phi) a.push_back(std::fabs(x));
        double m = 0;
        if (!a.empty()) {
          const size_t k = std::min(a.size() - 1, size_t(0.95 * double(a.size())));
          std::nth_element(a.begin(), a.begin() + long(k), a.end());
          m = a[k];
        }
        if (m < 1e-6) m = 1;
        L.phi_lo = -m, L.phi_hi = m;
        M.colours.resize(phi.size());
        for (size_t v = 0; v < phi.size(); ++v) M.colours[v] = caps::ramp_colour(L.ramp, 0.5 + 0.5 * phi[v] / m);
      } else if (L.surface_colour == 2) {
        M.colours.resize(M.vertices.size());
        for (size_t v = 0; v < M.vertices.size(); ++v) M.colours[v] = M.nearest[v] >= 0 ? caps::element_colour(f.atoms[size_t(M.nearest[v])].element) : 0x8FB8D8;
      }
    }
  } catch (const std::exception& e) {
    L.error = e.what();
  }
}

void refresh(caps_doc* d) {
  d->frame = d->traj.frame(d->current);
  if (d->smooth_window > 1 && d->traj.frames() > 1) {
    const auto p = caps::smoothed_positions(d->traj, d->current, d->smooth_window);
    for (size_t i = 0; i < d->frame.atoms.size() && i < p.size(); ++i) d->frame.atoms[i].pos = p[i];
  }
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
  run_doc_pipeline(d);
  prepare_appearance(d);
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
  r.cell_repeats = d->cell_repeats;
  // selection, focus and property colours refer to the frame's atoms; with a pipeline they map onto its particles
  auto to_shown = [&](int i) { return !d->pstate ? i : i >= 0 && size_t(i) < d->shown_of.size() ? d->shown_of[size_t(i)] : -1; };
  for (int k = 0; k < 4; ++k) if (o->highlight[k] >= 0 && to_shown(o->highlight[k]) >= 0) r.highlight.push_back(to_shown(o->highlight[k]));
  r.focus = to_shown(o->focus - 1);
  r.ambient_occlusion = o->ambient_occlusion != 0;
  if (!d->pstate && d->selection.size() == d->frame.atoms.size()) {   // the selection ringed (up to 50 000 atoms)
    for (size_t i = 0; i < d->selection.size() && r.highlight.size() < 50000; ++i) if (d->selection[i]) r.highlight.push_back(int(i));
  }
  if (!d->pstate) { r.segments = d->overlay; r.segments.insert(r.segments.end(), d->checks.begin(), d->checks.end()); }
  if (!d->pstate && d->look.active) {
    const AppearanceState& L = d->look;
    if (L.style.size() == d->frame.atoms.size()) {
      r.atom_style = L.style;
      for (auto& x : r.atom_style) if (x == 255) x = uint8_t(r.style);
    }
    if (L.colour == 4) {
      r.colour_by = caps::ColourBy::Property;
      r.property.clear();
      for (const auto& a : d->frame.atoms) r.property.push_back(a.charge);
      r.ramp = L.ramp;
      r.symmetric = true;
    } else if (L.colour >= 0) {
      r.colour_by = static_cast<caps::ColourBy>(std::clamp(L.colour, 0, 3));
    }
    if (L.mesh) r.meshes.push_back({L.mesh.get(), 0x8FB8D8, L.opacity});
    if (L.poly) r.meshes.push_back({L.poly.get(), 0x8FB8D8, 0.45f});
    r.segments.insert(r.segments.end(), L.ribbon.begin(), L.ribbon.end());
  }
  if (d->pstate) {
    const auto& st = *d->pstate;
    r.colours = st.colour;
    for (size_t i = 0; i < r.colours.size(); ++i) if (st.selected[i]) r.colours[i] = 0xE5484D;   // selected particles in red
    r.segments = st.segments;
    if (r.colour_by == caps::ColourBy::Property) caps::property_values(st, "DistanceToCOM", r.property);
  } else if (r.colour_by == caps::ColourBy::Property && r.property.size() != d->frame.atoms.size()) {
    r.property = d->dcom;
  }
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
    uo.qeq = F.charges == "qeq";
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

caps_doc* caps_open_staged(const char* path, const char* topology_path, int32_t max_frames, caps_open_progress_fn progress, void* user) {
  try {
    caps::OpenProgress p;
    p.max_frames = max_frames > 0 ? static_cast<size_t>(max_frames) : 0;
    if (progress) p.report = [&](int stage, double f, const std::string& detail) { return progress(stage, f, detail.c_str(), user) == 0; };
    auto* d = new caps_doc;
    d->traj = caps::open_file(path, topology_path ? topology_path : "", p);
    refresh(d);
    return d;
  } catch (const std::exception& e) {
    g_error = e.what();
  } catch (...) {
    g_error = "unknown error";
  }
  return nullptr;
}

int32_t caps_adopt_frames(caps_doc* dst, caps_doc* src) {
  if (!dst || !src || src->traj.topology.atoms.size() != dst->traj.topology.atoms.size()) return -1;
  dst->traj.positions = std::move(src->traj.positions);
  dst->traj.cells = std::move(src->traj.cells);
  dst->traj.timesteps = std::move(src->traj.timesteps);
  for (const auto& n : src->traj.topology.notes)
    if (std::find(dst->traj.topology.notes.begin(), dst->traj.topology.notes.end(), n) == dst->traj.topology.notes.end()) dst->traj.topology.notes.push_back(n);
  src->traj.positions.clear();
  src->traj.cells.clear();
  src->traj.timesteps.clear();
  if (dst->current >= dst->traj.frames()) dst->current = 0;
  refresh(dst);
  return static_cast<int32_t>(dst->traj.frames());
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
    caps::write_lammps_input(d->frame, ff, elec(), data_name && *data_name ? data_name : "system.data", tmp.string(), d->held_mol);
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
    F->charges = charges == 1 ? "gasteiger" : charges == 2 ? "keep" : charges == 3 ? "qeq" : "types";
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
    const auto img = d->renderer.render(shown(d), cam_of(cam), opts_of(d, opt));
    std::memcpy(rgba, img.rgba.data(), img.rgba.size());
    return 0;
  });
}

int32_t caps_pick(caps_doc* d, int32_t x, int32_t y) {
  const int k = d->renderer.pick(x, y);
  if (k < 0 || !d->pstate) return k;
  return size_t(k) < d->pstate->origin.size() ? d->pstate->origin[size_t(k)] : -1;   // the frame's atom under the pixel
}

int32_t caps_export_png(caps_doc* d, const caps_camera* cam, const caps_render_opts* opt, const char* path) {
  return guard([&] {
    caps::Renderer r;   // separate renderer so the view's pick buffer is not replaced
    caps::write_png(r.render(shown(d), cam_of(cam), opts_of(d, opt)), path);
    return 0;
  });
}

int32_t caps_export_svg(caps_doc* d, const caps_camera* cam, const caps_render_opts* opt, const char* path) {
  return guard([&] {
    std::ofstream f(path);
    if (!f) throw std::runtime_error(std::string("cannot write ") + path);
    f << caps::render_svg(shown(d), cam_of(cam), opts_of(d, opt));
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
          if (!cur.empty()) (cur == "cij_strain" || cur == "cij_run" || cur == "tensile" || cur == "tg" || cur == "pull_shear" || cur == "pull_normal" ? protocols : ids).push_back(cur);
          cur.clear();
          if (*c == 0) break;
        } else if (*c != ' ') cur += *c;
      }
    }
    if (ids.empty() && protocols.empty()) throw std::invalid_argument("no properties requested");
    caps::AnalyzeOptions o;
    o.exclude_mol = d->held_mol;   // a held surface is not a chain
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
      } else if (id == "pull_shear" || id == "pull_normal") {
        caps::System s = frame_copy();
        caps::PullOptions po;
        po.field = ff;
        po.energy = o.energy;
        po.normal = id == "pull_normal";
        po.axis = std::clamp(mo.axis, 0, 2);
        po.surface_mol = d->held_mol > 0 ? d->held_mol : 1;
        if (mo.temperature > 0) po.temperature = mo.temperature;
        if (mo.dt > 0) po.dt = mo.dt;
        if (mo.seed) po.seed = mo.seed;
        if (mo.equilibrate_ps != 0) po.equilibrate_ps = std::max(0.0, mo.equilibrate_ps);
        if (mo.max_strain > 0) po.distance = mo.max_strain;     // Å for the pull test
        if (mo.rate > 0) po.rate = mo.rate;                     // Å/ps for the pull test
        po.progress = [&](const caps::PullPoint& q) {
          char b[96];
          std::snprintf(b, sizeof b, "pull: %.2f Å, force %.2f kcal/mol/Å", q.displacement, q.force);
          return !cancelled(b, std::min(1.0, q.time_ps * po.rate / po.distance));
        };
        for (auto& q : caps::pull_properties(caps::run_pull(s, po), po.normal)) res.push_back(std::move(q));
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

double caps_view_scale(caps_doc* d, const caps_camera* cam, const caps_render_opts* opt) {
  try {
    return caps::view_scale(shown(d), cam_of(cam), opts_of(d, opt));
  } catch (...) {
    return 0;
  }
}

int32_t caps_bonded(caps_doc* d, int32_t i, int32_t* idx, int32_t cap) {
  return guard([&] {
    const auto& s = d->frame;
    if (i < 0 || size_t(i) >= s.atoms.size()) throw std::out_of_range("atom index out of range");
    int32_t n = 0;
    for (const auto& b : s.bonds) {
      const int32_t o = b.i == uint32_t(i) ? int32_t(b.j) : b.j == uint32_t(i) ? int32_t(b.i) : -1;
      if (o < 0) continue;
      if (n < cap) idx[n] = o;
      ++n;
    }
    return n;
  });
}

int32_t caps_molecule_index(caps_doc* d, int32_t* out, int32_t cap) {
  return guard([&] {
    int count = 0;
    const auto mol = d->frame.molecules(&count);
    for (size_t j = 0; j < mol.size() && j < size_t(std::max(0, cap)); ++j) out[j] = mol[j];
    return int32_t(count);
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
    po.length = j.num("length", 20);
    f = caps::nanoparticle(caps::read_cif(j.text("crystal")), po, &r);
    keep = {false, false, po.shape == caps::ParticleShape::Fibre};
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

extern "C" caps_doc* caps_grow_blend(const char* options_json, const caps_grow_opts* o, caps_progress_fn progress, void* user, char* report, int32_t cap) {
  try {
    const caps::Json j = caps::Json::parse(options_json && *options_json ? options_json : "{}");
    std::vector<caps::BlendComponent> comps;
    if (j.has("components"))
      for (const auto& c : j["components"].items()) {
        caps::BlendComponent b;
        b.spec = spec_from(c.has("spec") ? c["spec"].dump() : "{}");
        b.weight = c.num("weight", 1);
        b.chains = int(c.num("chains", 0));
        comps.push_back(b);
      }
    caps::BlendOptions bo;
    bo.chains = int(j.num("chains", 8));
    bo.density = j.num("density", 0.5);
    const std::string morph = j.text("morphology", "mixed");
    bo.morphology = morph == "slabs" ? caps::BlendMorphology::Slabs : morph == "droplet" ? caps::BlendMorphology::Droplet : caps::BlendMorphology::Mixed;
    if (o) {
      bo.grow.seed = o->seed;
      bo.grow.contact_scale = o->contact_scale > 0 ? o->contact_scale : 1.0;
      bo.grow.curve = o->curve != 0;
    }
    if (progress) bo.grow.progress = [&](int done, int total, int restarts) { return progress(done, total, restarts, user) == 0; };
    caps::BlendReport br;
    const caps::System s = caps::grow_blend(comps, bo, &br);
    std::string t;
    for (const auto& n : br.notes) t += n + "\n";
    report_out(t, report, cap);
    return doc_of(s);
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

extern "C" int32_t caps_pipeline_set(caps_doc* d, const char* json) {
  return guard([&] {
    if (!json || !*json) d->pipeline.reset();
    else {
      auto p = caps::pipeline_from_json(caps::Json::parse(json));
      if (p.steps.empty()) d->pipeline.reset();
      else d->pipeline = std::make_unique<caps::Pipeline>(std::move(p));
    }
    run_doc_pipeline(d);
    return 0;
  });
}

extern "C" int32_t caps_pipeline_result(caps_doc* d, char* json, int32_t cap) {
  try {
    if (!d->pstate) return report_out("", json, cap);
    return report_out(caps::pipeline_result_json(*d->pstate).dump(0), json, cap);
  } catch (const std::exception& e) {
    g_error = e.what();
    return -1;
  }
}

extern "C" int32_t caps_pipeline_particles(caps_doc* d, const char* filter, int32_t offset, int32_t count, char* json, int32_t cap) {
  try {
    caps::PipelineState plain;
    const caps::PipelineState* st = d->pstate.get();
    if (!st) {   // no pipeline: the frame itself
      plain = caps::run_pipeline(d->frame, caps::Pipeline{}, int(d->current), 0);
      st = &plain;
    }
    return report_out(caps::particles_json(*st, filter ? filter : "", size_t(std::max(0, offset)), size_t(std::max(0, count))).dump(0), json, cap);
  } catch (const std::exception& e) {
    g_error = e.what();
    return -1;
  }
}

extern "C" int32_t caps_pipeline_bonds(caps_doc* d, int32_t offset, int32_t count, char* json, int32_t cap) {
  try {
    caps::PipelineState plain;
    const caps::PipelineState* st = d->pstate.get();
    if (!st) {
      plain = caps::run_pipeline(d->frame, caps::Pipeline{}, int(d->current), 0);
      st = &plain;
    }
    return report_out(caps::bonds_json(*st, size_t(std::max(0, offset)), size_t(std::max(0, count))).dump(0), json, cap);
  } catch (const std::exception& e) {
    g_error = e.what();
    return -1;
  }
}

namespace {
// What an export writes, and the notes that say what was done.
void export_write(caps_doc* d, const std::string& fmt, const caps::Json& o, const std::string& path, std::vector<std::string>& notes) {
  const bool use_pipeline = o.has("pipeline") && o["pipeline"].kind() == caps::Json::Bool && o["pipeline"].boolean() && d->pstate;
  caps::System s = use_pipeline ? d->pstate->system : d->frame;
  if (use_pipeline) notes.push_back("the Visualize pipeline's particles (" + std::to_string(s.atoms.size()) + ")");
  const bool wrap = o.has("wrap") && o["wrap"].kind() == caps::Json::Bool && o["wrap"].boolean();
  if (wrap && s.cell.valid() && fmt != "lammps-data") {
    for (auto& a : s.atoms) a.pos = s.cell.wrap(a.pos);
    notes.push_back("positions wrapped into the cell");
  }
  const bool coeffs = !o.has("coeffs") || o["coeffs"].kind() != caps::Json::Bool || o["coeffs"].boolean();
  if (fmt == "pdb") caps::write_pdb(s, path);
  else if (fmt == "xyz") caps::write_xyz(s, path);
  else if (fmt == "mol2") caps::write_mol2(s, path);
  else if (fmt == "gro") caps::write_gro(s, path);
  else if (fmt == "lammps-dump") {
    caps::write_lammps_dump(d->traj, path);
    notes.push_back("every frame, unwrapped coordinates; the pipeline is not applied to a dump");
  } else if (fmt == "lammps-data") {
    notes.push_back("wrapped coordinates with image flags");
    const bool same = s.atoms.size() == d->frame.atoms.size() && !use_pipeline;
    if (coeffs && d->field && d->field->complete && same) {
      caps::write_lammps_data_ff(s, *d->field->ff, elec(), path);
      notes.push_back("coefficients from Field: " + d->field->ff->name);
    } else if (coeffs) {
      try {
        caps::write_lammps_data_ff(s, default_ff(s), elec(), path);
        notes.push_back("coefficients from the built-in force field (GAFF for C/H, else UFF)");
      } catch (const std::exception& e) {
        caps::write_lammps_data(s, path);
        notes.push_back(std::string("no coefficients: ") + e.what());
      }
    } else {
      caps::write_lammps_data(s, path);
      notes.push_back("structure only: no Coeffs sections");
    }
  } else {
    throw std::invalid_argument("unknown format " + fmt);
  }
}
}  // namespace

extern "C" int32_t caps_export_data(caps_doc* d, const char* path, const char* format, const char* options) {
  return guard([&] {
    std::vector<std::string> notes;
    export_write(d, format ? format : "lammps-data", caps::Json::parse(options && *options ? options : "{}"), path, notes);
    return 0;
  });
}

extern "C" int32_t caps_export_preview(caps_doc* d, const char* format, const char* options, int32_t lines, char* json, int32_t cap) {
  try {
    const std::string fmt = format ? format : "lammps-data";
    const auto tmp = std::filesystem::temp_directory_path() / ("caps_export_" + std::to_string(reinterpret_cast<uintptr_t>(d)) + "." + fmt);
    std::vector<std::string> notes;
    export_write(d, fmt, caps::Json::parse(options && *options ? options : "{}"), tmp.string(), notes);
    caps::Json j = caps::Json::object();
    j["bytes"] = double(std::filesystem::file_size(tmp));
    caps::Json out = caps::Json::array();
    std::ifstream in(tmp);
    std::string line;
    std::map<std::string, double> counts;
    int k = 0;
    while (std::getline(in, line)) {
      if (k < lines) out.push_back(line);
      ++k;
      if (fmt == "lammps-data" && k < 40) {   // header counts: "1300 atoms", "4 atom types"
        std::istringstream ls(line);
        double v;
        std::string a, b;
        if (ls >> v >> a) {
          if (ls >> b && b == "types") counts[a + "_types"] = v;
          else if (a == "atoms" || a == "bonds" || a == "angles" || a == "dihedrals" || a == "impropers") counts[a] = v;
        }
      }
    }
    in.close();
    std::filesystem::remove(tmp);
    j["lines"] = std::move(out);
    j["line_count"] = double(k);
    const caps::System& s = d->pstate && std::strstr(options ? options : "", "\"pipeline\":true") ? d->pstate->system : d->frame;
    if (counts.empty()) { counts["atoms"] = double(s.atoms.size()); counts["bonds"] = double(fmt == "xyz" || fmt == "gro" ? 0 : s.bonds.size()); }
    for (const auto& [key, v] : counts) j[key] = v;
    caps::Json nj = caps::Json::array();
    for (const auto& n : notes) nj.push_back(n);
    j["notes"] = std::move(nj);
    return report_out(j.dump(0), json, cap);
  } catch (const std::exception& e) {
    g_error = e.what();
    return -1;
  }
}

namespace {
std::pair<caps::BundleOptions, caps::Pipeline> bundle_options(const char* json) {
  const caps::Json o = caps::Json::parse(json && *json ? json : "{}");
  caps::BundleOptions b;
  b.name = o.text("name", "figure");
  b.input = o.text("input", "");
  b.topology = o.text("topology", "");
  b.frame = int(o.num("frame", 0));
  b.width = int(o.num("width", 1920));
  b.height = int(o.num("height", 1080));
  auto flag = [&](const char* k, bool def) { return o.has(k) && o[k].kind() == caps::Json::Bool ? o[k].boolean() : def; };
  b.include_input = flag("include_input", false);
  b.include_pipeline = flag("include_pipeline", true);
  b.include_data = flag("include_data", true);
  b.include_readme = flag("include_readme", true);
  caps::Pipeline p = o.has("pipeline") ? caps::pipeline_from_json(o["pipeline"]) : caps::Pipeline{};
  return {b, p};
}
}  // namespace

extern "C" int32_t caps_inspect_file(const char* path, const char* topology_path, char* json, int32_t cap) {
  try {
    const auto r = caps::inspect_file(path, topology_path ? topology_path : "");
    caps::Json j = caps::Json::object();
    j["format"] = r.format;
    j["format_name"] = r.format_name;
    j["atoms"] = double(r.atoms);
    j["frames"] = double(r.frames);
    j["bytes"] = double(r.bytes);
    j["bonds_from"] = r.bonds_from;
    j["units"] = r.units;
    caps::Json head = caps::Json::array();
    for (const auto& l : r.head) head.push_back(l);
    j["head"] = std::move(head);
    caps::Json cols = caps::Json::array();
    for (const auto& c : r.columns) {
      caps::Json o = caps::Json::object();
      o["name"] = c.name; o["maps_to"] = c.maps_to; o["kind"] = c.kind; o["used"] = c.used;
      cols.push_back(std::move(o));
    }
    j["columns"] = std::move(cols);
    caps::Json types = caps::Json::array();
    for (const auto& t : r.types) {
      caps::Json o = caps::Json::object();
      o["type"] = double(t.type); o["label"] = t.label; o["mass"] = t.mass; o["element"] = t.element;
      types.push_back(std::move(o));
    }
    j["types"] = std::move(types);
    caps::Json notes = caps::Json::array();
    for (const auto& n : r.notes) notes.push_back(n);
    j["notes"] = std::move(notes);
    return report_out(j.dump(0), json, cap);
  } catch (const std::exception& e) {
    g_error = e.what();
    return -1;
  }
}

extern "C" int32_t caps_bundle_preview(caps_doc* d, const char* options, char* json, int32_t cap) {
  try {
    auto [b, p] = bundle_options(options);
    b.include_figures = false;
    const auto files = caps::bundle_files(d->traj, p, b, caps::Camera{}, caps::RenderOptions{});
    caps::Json j = caps::Json::object();
    caps::Json fl = caps::Json::array();
    auto add = [&](const std::string& name, const std::string& note, double bytes, const std::string& hash) {
      caps::Json f = caps::Json::object();
      f["name"] = name;
      f["note"] = note;
      f["bytes"] = bytes;
      f["sha256"] = hash;
      fl.push_back(std::move(f));
    };
    add("figure.svg", "vector figure", -1, "");
    add("figure.png", std::to_string(b.width) + " × " + std::to_string(b.height) + " raster", -1, "");
    std::string first_csv, first_name;
    for (const auto& f : files) {
      add(f.name, f.content_note, double(f.bytes.size()), f.name == "provenance.json" ? "" : caps::sha256_hex(f.bytes));
      if (f.name == "provenance.json") j["provenance"] = f.bytes;
      if (first_csv.empty() && f.name.rfind("data/", 0) == 0) { first_csv = f.bytes; first_name = f.name; }
    }
    j["files"] = std::move(fl);
    j["csv"] = first_csv.substr(0, 2000);
    j["csv_name"] = first_name;
    return report_out(j.dump(0), json, cap);
  } catch (const std::exception& e) {
    g_error = e.what();
    return -1;
  }
}

extern "C" int32_t caps_bundle_write(caps_doc* d, const char* path, const char* options, const caps_camera* cam, const caps_render_opts* opt) {
  return guard([&] {
    auto [b, p] = bundle_options(options);
    caps::RenderOptions r = opts_of(d, opt);
    const auto files = caps::bundle_files(d->traj, p, b, cam_of(cam), r);
    caps::write_bundle(path, files);
    return int32_t(files.size());
  });
}

extern "C" int32_t caps_pipeline_series(caps_doc* d, int32_t stride, caps_analyze_progress_fn progress, void* user, char* json, int32_t cap) {
  try {
    const caps::Pipeline p = d->pipeline ? *d->pipeline : caps::Pipeline{};
    const auto t = caps::pipeline_series(d->traj, p, stride, d->wrap, [&](int done, int total) {
      return !progress || progress("frames", double(done) / std::max(1, total), user) == 0;
    });
    caps::Json j = caps::Json::object();
    caps::Json c = caps::Json::array();
    for (const auto& x : t.columns) c.push_back(x);
    j["columns"] = std::move(c);
    caps::Json rows = caps::Json::array();
    for (const auto& r : t.rows) {
      caps::Json row = caps::Json::array();
      for (double x : r) row.push_back(std::isfinite(x) ? caps::Json(x) : caps::Json());
      rows.push_back(std::move(row));
    }
    j["rows"] = std::move(rows);
    return report_out(j.dump(0), json, cap);
  } catch (const std::exception& e) {
    g_error = e.what();
    return -1;
  }
}

extern "C" int32_t caps_pipeline_to_yaml(const char* json, const char* name, const char* file, const char* topology, char* yaml, int32_t cap) {
  try {
    const auto p = caps::pipeline_from_json(caps::Json::parse(json && *json ? json : "[]"));
    return report_out(caps::pipeline_to_yaml(p, name ? name : "", file ? file : "", topology ? topology : ""), yaml, cap);
  } catch (const std::exception& e) {
    g_error = e.what();
    return -1;
  }
}

extern "C" int32_t caps_pipeline_from_yaml(const char* yaml, char* json, int32_t cap) {
  try {
    std::string name, file, topo;
    const auto p = caps::pipeline_from_yaml(yaml ? yaml : "", &name, &file, &topo);
    caps::Json j = caps::pipeline_to_json(p);
    j["name"] = name;
    j["file"] = file;
    j["topology"] = topo;
    return report_out(j.dump(0), json, cap);
  } catch (const std::exception& e) {
    g_error = e.what();
    return -1;
  }
}

extern "C" void caps_set_python(const char* package_dir, const char* interpreter) {
#ifdef _WIN32
  if (package_dir && *package_dir) _putenv_s("CAPS_PYTHON_PATH", package_dir);
  if (interpreter && *interpreter) _putenv_s("CAPS_PYTHON", interpreter);
#else
  if (package_dir && *package_dir) setenv("CAPS_PYTHON_PATH", package_dir, 1);
  if (interpreter && *interpreter) setenv("CAPS_PYTHON", interpreter, 1);
#endif
}

extern "C" int32_t caps_pipeline_catalogue(char* json, int32_t cap) {
  caps::Json a = caps::Json::array();
  for (const auto& [type, title, about] : caps::pipeline_step_catalogue()) {
    caps::Json o = caps::Json::object();
    o["type"] = type;
    o["title"] = title;
    o["about"] = about;
    a.push_back(std::move(o));
  }
  return report_out(a.dump(0), json, cap);
}

extern "C" int32_t caps_file_checks(caps_doc* d, char* json, int32_t cap) {
  try {
    return report_out(caps::file_checks_json(caps::file_checks(d->traj)), json, cap);
  } catch (const std::exception& e) {
    g_error = e.what();
    return -1;
  }
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

// ---------------------------------------------------------------- crystals from space groups (v20)

namespace {

caps::CrystalSpec crystal_spec_from(const caps::Json& j) {
  caps::CrystalSpec c;
  c.space_group = j.text("space_group", "P 1");
  c.a = j.num("a", 5), c.b = j.num("b", 5), c.c = j.num("c", 5);
  c.alpha = j.num("alpha", 90), c.beta = j.num("beta", 90), c.gamma = j.num("gamma", 90);
  c.tolerance = j.num("tolerance", 0.01);
  c.title = j.text("title");
  if (j.has("supercell") && j["supercell"].is_array())
    for (size_t k = 0; k < 3 && k < j["supercell"].size(); ++k) c.supercell[k] = std::max(1, int(j["supercell"][k].number()));
  if (j.has("sites") && j["sites"].is_array())
    for (const auto& x : j["sites"].items()) {
      caps::CrystalSite site;
      site.label = x.text("label");
      const std::string el = x.text("element");
      site.element = caps::element_from_symbol(el);
      if (site.element <= 0) throw std::invalid_argument("site " + site.label + ": unknown element '" + el + "'");
      site.frac = {x.num("x", 0), x.num("y", 0), x.num("z", 0)};
      site.occupancy = x.num("occupancy", 1);
      c.sites.push_back(site);
    }
  return c;
}

caps::Json crystal_spec_json(const caps::CrystalSpec& c) {
  caps::Json j = caps::Json::object();
  j["space_group"] = c.space_group;
  j["a"] = c.a, j["b"] = c.b, j["c"] = c.c, j["alpha"] = c.alpha, j["beta"] = c.beta, j["gamma"] = c.gamma;
  caps::Json sc = caps::Json::array();
  for (int n : c.supercell) sc.push_back(double(n));
  j["supercell"] = sc;
  j["tolerance"] = c.tolerance;
  caps::Json sites = caps::Json::array();
  for (const auto& s : c.sites) {
    caps::Json x = caps::Json::object();
    x["label"] = s.label;
    x["element"] = std::string(caps::element(s.element).symbol);
    x["x"] = s.frac[0], x["y"] = s.frac[1], x["z"] = s.frac[2];
    x["occupancy"] = s.occupancy;
    sites.push_back(x);
  }
  j["sites"] = sites;
  return j;
}

// The lattice centring of a setting: the first letter of its symbol; rhombohedral axes (":r") are primitive.
char crystal_centring(const caps::SpaceGroupSetting& sg) {
  if (sg.key.size() > 2 && sg.key.substr(sg.key.size() - 2) == ":r") return 'P';
  return sg.hm.empty() ? 'P' : sg.hm[0];
}

// The crystal of a spec; "primitive": true reduces a centred cell before the supercell repeats.
caps::System crystal_of(const caps::Json& j, caps::CrystalReport* rep) {
  caps::CrystalSpec c = crystal_spec_from(j);
  const bool primitive = j.num("primitive", 0) != 0 || (j.has("primitive") && j["primitive"].kind() == caps::Json::Bool && j["primitive"].boolean());
  if (!primitive) return caps::build_crystal(c, rep);
  const auto sc = c.supercell;
  c.supercell = {1, 1, 1};
  caps::System s = caps::build_crystal(c, rep);
  const auto* sg = caps::find_space_group(c.space_group);
  s = caps::primitive_cell(s, sg ? crystal_centring(*sg) : 'P');
  if (rep) rep->atoms_per_cell = s.atoms.size(), rep->volume = s.cell.volume();
  if (sc[0] * sc[1] * sc[2] > 1) s = caps::supercell(s, sc[0], sc[1], sc[2]);
  s.bonds = caps::crystal_bonds(s);
  return s;
}

}  // namespace

extern "C" int32_t caps_space_groups(char* json, int32_t cap) {
  caps::Json arr = caps::Json::array();
  for (const auto& sg : caps::space_group_settings()) {
    caps::Json x = caps::Json::object();
    x["key"] = sg.key;
    x["number"] = double(sg.number);
    x["hm"] = sg.hm;
    x["hall"] = sg.hall;
    x["system"] = caps::crystal_system(sg.number);
    arr.push_back(x);
  }
  return report_out(arr.dump(0), json, cap);
}

extern "C" int32_t caps_crystal_info(const char* spec_json, char* json, int32_t cap) {
  caps::Json j = caps::Json::object();
  try {
    const caps::Json spec = caps::Json::parse(spec_json && *spec_json ? spec_json : "{}");
    caps::CrystalReport rep;
    const caps::System s = crystal_of(spec, &rep);
    j["ok"] = true;
    j["key"] = rep.key, j["hm"] = rep.hm, j["hall"] = rep.hall, j["system"] = rep.system;
    j["number"] = double(rep.number), j["operations"] = double(rep.operations);
    j["atoms_per_cell"] = double(rep.atoms_per_cell), j["atoms"] = double(s.atoms.size()), j["bonds"] = double(s.bonds.size());
    caps::Json m = caps::Json::array();
    for (int k : rep.multiplicity) m.push_back(double(k));
    j["multiplicity"] = m;
    j["volume"] = rep.volume;
    j["density"] = s.density();
    std::vector<size_t> all(s.atoms.size());
    for (size_t i = 0; i < all.size(); ++i) all[i] = i;
    j["formula"] = caps::formula_of(s, all);
    const auto* sg = caps::find_space_group(rep.key);
    j["centring"] = std::string(1, sg ? crystal_centring(*sg) : 'P');
    caps::Json notes = caps::Json::array();
    for (const auto& n : rep.notes) notes.push_back(n);
    j["notes"] = notes;
  } catch (const std::exception& e) {
    j = caps::Json::object();
    j["ok"] = false;
    j["error"] = std::string(e.what());
  }
  return report_out(j.dump(0), json, cap);
}

extern "C" caps_doc* caps_crystal_build(const char* spec_json, char* report, int32_t cap) {
  try {
    caps::CrystalReport rep;
    const caps::System s = crystal_of(caps::Json::parse(spec_json && *spec_json ? spec_json : "{}"), &rep);
    std::string notes;
    for (const auto& n : s.notes) notes += n + "\n";
    report_out(notes, report, cap);
    caps_doc* d = doc_of(s);
    const caps::Json j = caps::Json::parse(spec_json && *spec_json ? spec_json : "{}");
    if (j.has("supercell") && j["supercell"].is_array())
      for (size_t k = 0; k < 3 && k < j["supercell"].size(); ++k) d->cell_repeats[k] = std::max(1, int(j["supercell"][k].number()));
    return d;
  } catch (const std::exception& e) {
    g_error = e.what();
    return nullptr;
  }
}

extern "C" int32_t caps_crystal_symmetrize(const char* spec_json, double snap, char* json, int32_t cap) {
  caps::Json j = caps::Json::object();
  try {
    const caps::Json in = caps::Json::parse(spec_json && *spec_json ? spec_json : "{}");
    int moved = 0;
    caps::Json out = crystal_spec_json(caps::symmetrize_sites(crystal_spec_from(in), snap > 0 ? snap : 0.3, &moved));
    if (in.has("primitive")) out["primitive"] = in["primitive"];
    j["ok"] = true;
    j["moved"] = double(moved);
    j["spec"] = out;
  } catch (const std::exception& e) {
    j["ok"] = false;
    j["error"] = std::string(e.what());
  }
  return report_out(j.dump(0), json, cap);
}

extern "C" int32_t caps_crystal_find_symmetry(const char* spec_json, const char* cif_path, double tolerance, char* json, int32_t cap) {
  caps::Json j = caps::Json::object();
  try {
    caps::System s;
    caps::CrystalSpec base;
    if (cif_path && *cif_path) {
      s = caps::read_cif(cif_path);
      base.title = s.title;
    } else {
      const caps::Json in = caps::Json::parse(spec_json && *spec_json ? spec_json : "{}");
      base = crystal_spec_from(in);
      caps::CrystalSpec one = base;
      one.supercell = {1, 1, 1};
      s = caps::build_crystal(one);
    }
    const caps::SymmetryFound f = caps::find_symmetry(s, tolerance > 0 ? tolerance : 0.1);
    caps::CrystalSpec out = base;
    const auto& c = s.cell;
    const double A = caps::norm(c.a), B = caps::norm(c.b), C = caps::norm(c.c);
    auto ang = [](const caps::Vec3& u, const caps::Vec3& v) { return std::acos(std::clamp(caps::dot(u, v) / (caps::norm(u) * caps::norm(v)), -1.0, 1.0)) * 57.29577951308232; };
    out.a = A, out.b = B, out.c = C, out.alpha = ang(c.b, c.c), out.beta = ang(c.a, c.c), out.gamma = ang(c.a, c.b);
    out.space_group = f.key;
    out.sites = f.sites;
    j["ok"] = true;
    j["key"] = f.key, j["hm"] = f.hm, j["system"] = f.system;
    j["number"] = double(f.number), j["operations"] = double(f.operations);
    j["atoms"] = double(s.atoms.size());
    j["spec"] = crystal_spec_json(out);
  } catch (const std::exception& e) {
    j["ok"] = false;
    j["error"] = std::string(e.what());
  }
  return report_out(j.dump(0), json, cap);
}

// ---------------------------------------------------------------- peptides (v20)

namespace {

caps::PeptideOptions peptide_options(const caps::Json& j) {
  caps::PeptideOptions o;
  o.sequence = j.text("sequence");
  o.structure = j.text("structure");
  auto tri = [&](const char* k, std::array<double, 3>& t) {
    if (j.has(k) && j[k].is_array() && j[k].size() >= 2)
      for (size_t q = 0; q < 3 && q < j[k].size(); ++q) t[q] = j[k][q].number();
  };
  tri("helix", o.helix);
  tri("strand", o.strand);
  tri("ppii", o.ppii);
  o.n_term = j.text("n_term", "NH3+");
  o.c_term = j.text("c_term", "COO-");
  o.ph = j.num("ph", 7.0);
  auto flag = [&](const char* k, bool def) {
    if (!j.has(k)) return def;
    return j[k].kind() == caps::Json::Bool ? j[k].boolean() : j[k].number() != 0;
  };
  o.neutral = flag("neutral", false);
  o.cleanup = flag("cleanup", true);
  o.seed = uint64_t(j.num("seed", 1));
  return o;
}

// A tube through the CA atoms (Catmull–Rom, six pieces per residue) coloured by secondary structure.
std::vector<caps::Segment> peptide_ribbon(const caps::System& s, const std::string& ss) {
  std::vector<caps::Vec3> ca;
  for (const auto& a : s.atoms) if (a.name == "CA" && a.resname != "ACE" && a.resname != "NME") ca.push_back(a.pos);
  std::vector<caps::Segment> out;
  if (ca.size() < 2) return out;
  auto colour = [&](size_t i) -> unsigned {
    const char c = i < ss.size() ? ss[i] : 'C';
    return c == 'H' ? 0xE07A5F : c == 'E' ? 0x6FA8DC : c == 'P' ? 0x9B7BD6 : 0x8A9097;
  };
  auto radius = [&](size_t i) { const char c = i < ss.size() ? ss[i] : 'C'; return c == 'H' ? 0.9 : c == 'E' ? 0.7 : 0.35; };
  for (size_t i = 0; i + 1 < ca.size(); ++i) {
    const caps::Vec3 &p0 = ca[i ? i - 1 : 0], &p1 = ca[i], &p2 = ca[i + 1], &p3 = ca[std::min(i + 2, ca.size() - 1)];
    caps::Vec3 prev = p1;
    for (int k = 1; k <= 6; ++k) {
      const double t = k / 6.0, t2 = t * t, t3 = t2 * t;
      caps::Vec3 q;
      for (int d = 0; d < 3; ++d)
        q[d] = 0.5 * (2 * p1[d] + (-p0[d] + p2[d]) * t + (2 * p0[d] - 5 * p1[d] + 4 * p2[d] - p3[d]) * t2 + (-p0[d] + 3 * p1[d] - 3 * p2[d] + p3[d]) * t3);
      caps::Segment seg;
      seg.a = prev, seg.b = q;
      const size_t owner = k <= 3 ? i : i + 1;
      seg.rgb = colour(owner);
      seg.radius = radius(owner);
      out.push_back(seg);
      prev = q;
    }
  }
  return out;
}

}  // namespace

extern "C" int32_t caps_peptide_info(const char* options_json, char* json, int32_t cap) {
  caps::Json j = caps::Json::object();
  try {
    caps::PeptideOptions o = peptide_options(caps::Json::parse(options_json && *options_json ? options_json : "{}"));
    o.cleanup = false;
    caps::PeptideReport rep;
    caps::build_peptide(o, &rep);
    j["ok"] = true;
    j["residues"] = double(rep.residues), j["atoms"] = double(rep.atoms), j["charge"] = double(rep.charge);
    j["formula"] = rep.formula, j["mass"] = rep.mass, j["smiles"] = rep.smiles, j["structure"] = rep.structure;
  } catch (const std::exception& e) {
    j["ok"] = false;
    j["error"] = std::string(e.what());
  }
  return report_out(j.dump(0), json, cap);
}

extern "C" caps_doc* caps_peptide_build(const char* options_json, char* report, int32_t cap) {
  try {
    const caps::Json j = caps::Json::parse(options_json && *options_json ? options_json : "{}");
    caps::PeptideReport rep;
    const caps::System s = caps::build_peptide(peptide_options(j), &rep);
    std::string notes;
    for (const auto& n : s.notes) notes += n + "\n";
    report_out(notes, report, cap);
    caps_doc* d = doc_of(s);
    if (j.has("ribbon") && (j["ribbon"].kind() == caps::Json::Bool ? j["ribbon"].boolean() : j["ribbon"].number() != 0)) d->overlay = peptide_ribbon(s, rep.structure);
    return d;
  } catch (const std::exception& e) {
    g_error = e.what();
    return nullptr;
  }
}

extern "C" int32_t caps_fasta_sequence(const char* text, char* seq, int32_t cap) {
  return report_out(caps::parse_fasta(text ? text : ""), seq, cap);
}

// ---------------------------------------------------------------- solvation (v20)

namespace {

caps::SolvateOptions solvate_options(const caps::Json& j) {
  caps::SolvateOptions o;
  o.shape = int(j.num("shape", 0));
  o.edge = j.num("edge", 30);
  if (j.has("edges") && j["edges"].is_array() && j["edges"].size() == 3) o.edges = {j["edges"][0].number(), j["edges"][1].number(), j["edges"][2].number()};
  o.padding = j.num("padding", 10);
  o.tolerance = j.num("tolerance", 2.0);
  o.solvent = j.text("solvent", "water");
  o.water_model = j.text("water_model", "TIP4P/2005");
  o.density = j.num("density", 0);
  o.molecules = int(j.num("molecules", 0));
  o.ion_mode = int(j.num("ion_mode", 2));
  o.salt = j.text("salt", "NaCl");
  o.concentration = j.num("concentration", 0.15);
  o.cations = int(j.num("cations", 0));
  o.anions = int(j.num("anions", 0));
  o.seed = uint64_t(j.num("seed", 1));
  return o;
}

caps::System solute_of(caps_doc* d) {
  caps::System s = d->frame;
  if (s.cell.valid() && !s.unwrapped) caps::make_molecules_whole(s);
  return s;
}

caps::Json plan_json(const caps::SolvatePlan& p) {
  caps::Json j = caps::Json::object();
  caps::Json box = caps::Json::array();
  for (int k = 0; k < 3; ++k) box.push_back(p.box[k]);
  j["box"] = box;
  j["box_volume"] = p.box_volume, j["solute_volume"] = p.solute_volume, j["free_volume"] = p.free_volume;
  j["solute_atoms"] = double(p.solute_atoms), j["solute_charge"] = p.solute_charge;
  j["solvent"] = double(p.solvent), j["cations"] = double(p.cations), j["anions"] = double(p.anions);
  j["solvent_name"] = p.solvent_name, j["cation"] = p.cation, j["anion"] = p.anion;
  j["solvent_mass"] = p.solvent_mass, j["concentration"] = p.concentration, j["density"] = p.density;
  caps::Json notes = caps::Json::array();
  for (const auto& n : p.notes) notes.push_back(n);
  j["notes"] = notes;
  return j;
}

}  // namespace

extern "C" int32_t caps_solvent_library(char* json, int32_t cap) {
  caps::Json j = caps::Json::object();
  caps::Json sv = caps::Json::array();
  for (const auto& s : caps::solvent_library()) {
    caps::Json x = caps::Json::object();
    x["id"] = s.id, x["name"] = s.name, x["smiles"] = s.smiles, x["density"] = s.density, x["use"] = s.use;
    sv.push_back(x);
  }
  j["solvents"] = sv;
  caps::Json sl = caps::Json::array();
  for (const auto& s : caps::salt_library()) {
    caps::Json x = caps::Json::object();
    x["id"] = s.id, x["cation"] = s.cation, x["anion"] = s.anion, x["zc"] = double(s.zc), x["za"] = double(s.za);
    sl.push_back(x);
  }
  j["salts"] = sl;
  return report_out(j.dump(0), json, cap);
}

extern "C" int32_t caps_solvate_plan(caps_doc* solute, const char* options_json, char* json, int32_t cap) {
  caps::Json j = caps::Json::object();
  try {
    const caps::SolvateOptions o = solvate_options(caps::Json::parse(options_json && *options_json ? options_json : "{}"));
    caps::System s;
    if (solute) s = solute_of(solute);
    j = plan_json(caps::solvate_plan(solute ? &s : nullptr, o));
    j["ok"] = true;
  } catch (const std::exception& e) {
    j = caps::Json::object();
    j["ok"] = false;
    j["error"] = std::string(e.what());
  }
  return report_out(j.dump(0), json, cap);
}

extern "C" caps_doc* caps_solvate(caps_doc* solute, const char* options_json, caps_stage_progress_fn progress, void* user, char* report, int32_t cap) {
  try {
    caps::SolvateOptions o = solvate_options(caps::Json::parse(options_json && *options_json ? options_json : "{}"));
    if (progress)
      o.progress = [&](const caps::PackProgress& p) {
        const int stage = p.stage == "insertion" ? 0 : p.stage == "optimisation" ? 1 : 2;
        return progress(stage, p.loop, p.loops, p.dmin, p.bad, user) == 0;
      };
    caps::System s;
    if (solute) s = solute_of(solute);
    caps::SolvateReport rep;
    const caps::System out = caps::solvate(solute ? &s : nullptr, o, &rep);
    std::string notes;
    for (const auto& n : out.notes) notes += n + "\n";
    char b[160];
    std::snprintf(b, sizeof b, "min. distance %.2f Å · %d molecules · %d atoms", rep.pack.dmin, rep.pack.molecules, rep.pack.atoms);
    notes += b;
    report_out(notes, report, cap);
    return doc_of(out);
  } catch (const std::exception& e) {
    g_error = e.what();
    return nullptr;
  }
}

// ---------------------------------------------------------------- appearance (v20)

extern "C" int32_t caps_set_appearance(caps_doc* d, const char* json) {
  return guard([&] {
    const caps::Json j = caps::Json::parse(json && *json ? json : "{}");
    AppearanceState& L = d->look;
    L.active = !j.has("active") || (j["active"].kind() == caps::Json::Bool ? j["active"].boolean() : j["active"].number() != 0);
    L.layers.clear();
    static const std::map<std::string, caps::Style> styles = {
        {"ball_and_stick", caps::Style::BallAndStick}, {"space_filling", caps::Style::SpaceFilling}, {"sticks", caps::Style::Sticks},
        {"no_hydrogens", caps::Style::NoHydrogens}, {"wireframe", caps::Style::Wireframe}, {"polyhedra", caps::Style::Polyhedra},
        {"ribbon", caps::Style::Ribbon}, {"hidden", caps::Style::Hidden}};
    if (j.has("layers") && j["layers"].is_array())
      for (const auto& l : j["layers"].items()) {
        const auto it = styles.find(l.text("style", "ball_and_stick"));
        if (it == styles.end()) throw std::invalid_argument("unknown style '" + l.text("style") + "'");
        L.layers.emplace_back(l.text("expression"), int(it->second));
      }
    const std::string colour = j.text("colour", "");
    L.colour = colour == "element" ? 0 : colour == "molecule" ? 1 : colour == "type" ? 2 : colour == "distance" ? 3 : colour == "charge" ? 4 : -1;
    const std::string ramp = j.text("ramp", "blue_orange");
    L.ramp = ramp == "viridis" ? caps::Ramp::Viridis : ramp == "red_white_blue" ? caps::Ramp::RedWhiteBlue : caps::Ramp::BlueOrange;
    const caps::Json sf = j.has("surface") ? j["surface"] : caps::Json::object();
    const std::string kind = sf.text("kind", "none");
    L.surface = kind == "accessible" ? 1 : kind == "vdw" ? 2 : kind == "excluded" ? 3 : 0;
    L.probe = sf.num("probe", 1.4);
    L.spacing = sf.num("spacing", 0.6);
    L.opacity = float(std::clamp(sf.num("opacity", 0.6), 0.05, 1.0));
    L.surface_expr = sf.text("expression");
    const std::string sc = sf.text("colour", "potential");
    L.surface_colour = sc == "uniform" ? 0 : sc == "atom" ? 2 : 1;
    prepare_appearance(d);
    if (!L.error.empty()) throw std::invalid_argument(L.error);
    return 0;
  });
}

extern "C" int32_t caps_appearance_info(caps_doc* d, char* json, int32_t cap) {
  caps::Json j = caps::Json::object();
  const AppearanceState& L = d->look;
  j["active"] = L.active;
  j["error"] = L.error;
  caps::Json counts = caps::Json::object();
  static const char* names[] = {"ball_and_stick", "space_filling", "sticks", "no_hydrogens", "backbone", "wireframe", "polyhedra", "ribbon", "hidden"};
  std::map<std::string, double> c;
  for (uint8_t x : L.style) if (x < 9) c[names[x]] += 1; else c["view"] += 1;
  for (const auto& [k, v] : c) counts[k] = v;
  j["styles"] = counts;
  double qlo = 0, qhi = 0;
  for (const auto& a : d->frame.atoms) qlo = std::min(qlo, a.charge), qhi = std::max(qhi, a.charge);
  caps::Json q = caps::Json::array();
  q.push_back(qlo), q.push_back(qhi);
  j["charge"] = q;
  if (L.mesh) {
    caps::Json sf = caps::Json::object();
    sf["area"] = L.mesh->area();
    sf["vertices"] = double(L.mesh->vertices.size());
    sf["triangles"] = double(L.mesh->triangles.size());
    caps::Json phi = caps::Json::array();
    phi.push_back(L.phi_lo), phi.push_back(L.phi_hi);
    sf["potential"] = phi;
    j["surface"] = sf;
  }
  if (L.poly) j["polyhedra"] = double(L.poly->triangles.size());
  return report_out(j.dump(0), json, cap);
}

extern "C" int32_t caps_atom_labels(caps_doc* d, const char* kind, char* json, int32_t cap) {
  caps::Json arr = caps::Json::array();
  const std::string k = kind ? kind : "element";
  const caps::System& f = d->frame;
  if (k == "rs") {
    for (const auto& x : caps::stereo_labels(f)) arr.push_back(x);
  } else {
    for (const auto& a : f.atoms) {
      if (k == "charge") {
        char b[24];
        std::snprintf(b, sizeof b, "%+.2f", a.charge);
        arr.push_back(std::string(b));
      } else if (k == "type") {
        std::string t;
        for (const auto& ti : f.types) if (ti.type == a.type) t = ti.label;
        arr.push_back(t.empty() ? std::to_string(a.type) : t);
      } else if (k == "name") {
        arr.push_back(a.name);
      } else {
        arr.push_back(std::string(caps::element(a.element).symbol));
      }
    }
  }
  return report_out(arr.dump(0), json, cap);
}

extern "C" int32_t caps_project_atoms(caps_doc* d, const caps_camera* cam, const caps_render_opts* opt, float* xyv, int32_t count) {
  return guard([&] {
    if (d->pstate) throw std::runtime_error("labels are drawn without a Visualize pipeline");
    const auto p = d->renderer.project(d->frame, cam_of(cam), opts_of(d, opt));
    const size_t n = std::min(p.size() / 3, size_t(std::max(0, count)));
    std::memcpy(xyv, p.data(), n * 3 * sizeof(float));
    return int32_t(n);
  });
}

// ---------------------------------------------------------------- trajectory player (v20)

extern "C" int32_t caps_trajectory_series(caps_doc* d, const char* options_json, caps_series_progress_fn progress, void* user, char* json, int32_t cap) {
  caps::Json j = caps::Json::object();
  try {
    const caps::Json o = caps::Json::parse(options_json && *options_json ? options_json : "{}");
    caps::TrajectorySeriesOptions so;
    so.molecule = int(o.num("molecule", 0));
    so.dt_fs = o.num("dt_fs", 1.0);
    so.stride = int(o.num("stride", 1));
    caps::ThermoLog log;
    const std::string log_path = o.text("log");
    if (!log_path.empty()) { log = caps::read_lammps_log(log_path); so.log = &log; }
    if (progress) so.progress = [&](int done, int total) { return progress(done, total, user) == 0; };
    const caps::DataTable T = caps::trajectory_series(d->traj, so);
    caps::System first = d->traj.frame(0);
    if (!first.unwrapped) caps::make_molecules_whole(first);
    const caps::ChainEnds e = caps::chain_ends(first, so.molecule);
    j["ok"] = true;
    caps::Json cols = caps::Json::array();
    for (const auto& c : T.columns) cols.push_back(c);
    j["columns"] = cols;
    caps::Json rows = caps::Json::array();
    for (const auto& r : T.rows) {
      caps::Json row = caps::Json::array();
      for (double v : r) row.push_back(std::isnan(v) ? caps::Json() : caps::Json(v));
      rows.push_back(row);
    }
    j["rows"] = rows;
    j["molecule"] = double(e.molecule);
    caps::Json ends = caps::Json::array();
    ends.push_back(double(e.first)), ends.push_back(double(e.last));
    j["ends"] = ends;
    // the frame nearest each run's first printed step (checkpoint markers on the timeline)
    caps::Json runs = caps::Json::array();
    if (so.log)
      for (size_t r : log.run_starts) {
        if (r >= log.rows.size()) continue;
        const double step = log.rows[r][0];
        size_t best = 0;
        double gap = 1e300;
        for (size_t k = 0; k < d->traj.timesteps.size(); ++k)
          if (std::fabs(double(d->traj.timesteps[k]) - step) < gap) gap = std::fabs(double(d->traj.timesteps[k]) - step), best = k;
        runs.push_back(double(best));
      }
    j["run_frames"] = runs;
    j["log_rows"] = double(so.log ? log.rows.size() : 0);
  } catch (const std::exception& e) {
    j = caps::Json::object();
    j["ok"] = false;
    j["error"] = std::string(e.what());
  }
  return report_out(j.dump(0), json, cap);
}

extern "C" void caps_set_smoothing(caps_doc* d, int32_t window) {
  d->smooth_window = std::clamp(int(window), 1, 101);
  refresh(d);
}

// ---------------------------------------------------------------- torsion scan (v20)

extern "C" int32_t caps_torsion_scan(caps_doc* d, const char* options_json, caps_series_progress_fn progress, void* user, char* json, int32_t cap) {
  caps::Json j = caps::Json::object();
  try {
    const caps::Json o = caps::Json::parse(options_json && *options_json ? options_json : "{}");
    caps::TorsionScanOptions so;
    if (!o.has("atoms") || !o["atoms"].is_array() || o["atoms"].size() != 4) throw std::invalid_argument("choose four atoms (a, b, c, d) for the torsion");
    for (size_t k = 0; k < 4; ++k) so.atoms[k] = int(o["atoms"][k].number());
    so.from = o.num("from", -180), so.to = o.num("to", 180), so.step = o.num("step", 15);
    so.relax = o.has("relax") && (o["relax"].kind() == caps::Json::Bool ? o["relax"].boolean() : o["relax"].number() != 0);
    so.ftol = o.num("ftol", 0.1);
    if (progress) so.progress = [&](int done, int total) { return progress(done, total, user) == 0; };
    // the document's force field when Field assigned one, else UFF
    caps::System s = d->frame;
    std::shared_ptr<const caps::ForceField> ff;
    std::string ffname;
    if (o.text("forcefield", "auto") != "uff" && d->field && d->field->ff && d->field->complete) ff = d->field->ff, ffname = d->field->ff->name;
    if (!ff) {
      ff = std::make_shared<caps::ForceField>(caps::assign_uff(s));
      ffname = "UFF";
    }
    const caps::TorsionScanResult R = caps::torsion_scan(s, *ff, so);
    d->scan_frames.clear();
    for (const auto& p : R.points) d->scan_frames.push_back(p.positions);
    d->scan_original.clear();
    for (const auto& a : d->frame.atoms) d->scan_original.push_back(a.pos);
    j["ok"] = true;
    j["forcefield"] = ffname;
    caps::Json pts = caps::Json::array();
    for (const auto& p : R.points) {
      caps::Json x = caps::Json::object();
      x["phi"] = p.phi, x["energy"] = p.energy - R.minimum, x["dihedral"] = p.terms.dihedral, x["vdw"] = p.terms.vdw, x["coulomb"] = p.terms.coulomb;
      x["angle"] = p.terms.angle, x["bond"] = p.terms.bond, x["total"] = p.energy;
      pts.push_back(x);
    }
    j["points"] = pts;
    caps::Json cf = caps::Json::array();
    for (const auto& c : R.conformers) {
      caps::Json x = caps::Json::object();
      x["phi"] = c.phi, x["energy"] = c.energy, x["state"] = c.state;
      cf.push_back(x);
    }
    j["conformers"] = cf;
    j["barrier"] = R.barrier, j["phi_start"] = R.phi_start, j["moving"] = double(R.moving.size());
    caps::Json notes = caps::Json::array();
    for (const auto& n : R.notes) notes.push_back(n);
    j["notes"] = notes;
  } catch (const std::exception& e) {
    j = caps::Json::object();
    j["ok"] = false;
    j["error"] = std::string(e.what());
  }
  return report_out(j.dump(0), json, cap);
}

extern "C" int32_t caps_torsion_show(caps_doc* d, int32_t index) {
  return guard([&] {
    const std::vector<caps::Vec3>* p = nullptr;
    if (index < 0) p = &d->scan_original;
    else if (size_t(index) < d->scan_frames.size()) p = &d->scan_frames[size_t(index)];
    else throw std::out_of_range("no such scan point");
    if (p->size() != d->traj.positions.at(d->current).size()) throw std::runtime_error("the scan belongs to another structure");
    d->traj.positions[d->current] = *p;
    refresh(d);
    return 0;
  });
}

extern "C" int32_t caps_default_torsion(caps_doc* d, int32_t* atoms) {
  const auto t = caps::default_torsion(d->frame);
  for (int k = 0; k < 4; ++k) atoms[k] = t[size_t(k)];
  return t[0] >= 0 ? 0 : -1;
}

// ---------------------------------------------------------------- editing, selection, tacticity (v20)

namespace {

std::vector<uint32_t> atoms_of(caps_doc* d, const caps::Json& j) {
  std::vector<uint32_t> out;
  const size_t n = d->frame.atoms.size();
  if (j.has("atoms") && j["atoms"].is_array()) {
    for (const auto& x : j["atoms"].items()) {
      const long v = long(x.number());
      if (v < 0 || size_t(v) >= n) throw std::out_of_range("atom " + std::to_string(v + 1) + " is not in the structure");
      out.push_back(uint32_t(v));
    }
  } else if (j.text("atoms") == "selection") {
    for (size_t i = 0; i < d->selection.size() && i < n; ++i) if (d->selection[i]) out.push_back(uint32_t(i));
  }
  return out;
}

void push_undo(caps_doc* d, const std::string& what) {
  caps_doc::Snapshot sn;
  sn.topology = d->traj.topology;
  sn.positions = d->traj.positions.at(d->current);
  sn.what = what;
  d->undo.push_back(std::move(sn));
  if (d->undo.size() > 100) d->undo.erase(d->undo.begin());
  d->redo.clear();
}

void store(caps_doc* d, const caps::System& s) {
  d->traj.topology = s;
  std::vector<caps::Vec3> p;
  for (const auto& a : s.atoms) p.push_back(a.pos);
  d->traj.positions[d->current] = std::move(p);
  if (d->current < d->traj.cells.size()) d->traj.cells[d->current] = s.cell;
  d->field.reset();   // the typing no longer matches
  d->scan_frames.clear();
  if (d->selection.size() != s.atoms.size()) d->selection.assign(s.atoms.size(), 0);
  refresh(d);
}

}  // namespace

extern "C" int32_t caps_edit(caps_doc* d, const char* json, char* out, int32_t cap) {
  caps::Json r = caps::Json::object();
  try {
    if (d->traj.frames() != 1) throw std::runtime_error("editing works on a single structure: this file has " + std::to_string(d->traj.frames()) +
                                                        " frames (save the frame you want as its own file)");
    const caps::Json j = caps::Json::parse(json && *json ? json : "{}");
    const std::string op = j.text("op");
    caps::System s = d->traj.frame(d->current);
    std::string what;
    caps::Json added = caps::Json::array();
    auto element_of = [&](const std::string& sym) {
      const int z = caps::element_from_symbol(sym);
      if (z <= 0) throw std::invalid_argument("unknown element '" + sym + "'");
      return z;
    };
    if (op == "element") {
      const auto at = atoms_of(d, j);
      if (at.empty()) throw std::invalid_argument("pick the atoms to change");
      const int z = element_of(j.text("element"));
      for (uint32_t a : at) caps::set_element(s, a, z);
      what = "Change " + std::to_string(at.size()) + " atom(s) to " + caps::element(z).symbol;
    } else if (op == "charge") {
      const auto at = atoms_of(d, j);
      for (uint32_t a : at) s.atoms[a].charge = j.num("charge", 0);
      s.has_charges = true;
      what = "Set the charge of " + std::to_string(at.size()) + " atom(s)";
    } else if (op == "add_atom") {
      const int to = int(j.num("to", -1));
      const int z = element_of(j.text("element", "C"));
      caps::Vec3 at{0, 0, 0};
      if (to < 0 && !s.atoms.empty()) {   // an isolated atom beside the structure
        caps::Vec3 hi{-1e300, -1e300, -1e300};
        for (const auto& a : s.atoms) for (int k = 0; k < 3; ++k) hi[k] = std::max(hi[k], a.pos[k]);
        at = {hi[0] + 2.5, hi[1], hi[2]};
      }
      const uint32_t id = caps::add_atom(s, to, z, int(j.num("order", 1)), int(j.num("geometry", 0)), at);
      if (j.has("charge")) s.atoms[id].charge = j.num("charge", 0);
      added.push_back(double(id));
      what = std::string("Place ") + caps::element(z).symbol + (to >= 0 ? " on atom " + std::to_string(to + 1) : "");
    } else if (op == "bond") {
      const uint32_t a = uint32_t(j.num("i", -1)), b = uint32_t(j.num("j", -1));
      caps::add_bond(s, a, b, int(j.num("order", 1)));
      what = "Bond " + std::to_string(a + 1) + "–" + std::to_string(b + 1);
    } else if (op == "unbond") {
      const uint32_t a = uint32_t(j.num("i", -1)), b = uint32_t(j.num("j", -1));
      if (!caps::remove_bond(s, a, b)) throw std::invalid_argument("those atoms are not bonded");
      what = "Break bond " + std::to_string(a + 1) + "–" + std::to_string(b + 1);
    } else if (op == "delete") {
      const auto at = atoms_of(d, j);
      if (at.empty()) throw std::invalid_argument("pick the atoms to delete");
      std::vector<char> m(s.atoms.size(), 0);
      for (uint32_t a : at) m[a] = 1;
      caps::delete_atoms(s, m);
      d->selection.assign(s.atoms.size(), 0);
      what = "Delete " + std::to_string(at.size()) + " atom(s)";
    } else if (op == "add_h") {
      const auto at = atoms_of(d, j);
      std::vector<char> m;
      if (!at.empty()) { m.assign(s.atoms.size(), 0); for (uint32_t a : at) m[a] = 1; }
      const size_t before = s.atoms.size();
      const int k = caps::add_hydrogens(s, m);
      if (k == 0) throw std::invalid_argument("no atom lacks hydrogens");
      for (size_t i = before; i < s.atoms.size(); ++i) added.push_back(double(i));
      what = "Add " + std::to_string(k) + " hydrogens";
    } else if (op == "invert") {
      const uint32_t c = uint32_t(j.num("centre", -1));
      caps::invert_centre(s, c);
      what = "Invert centre " + std::to_string(c + 1);
    } else if (op == "tacticity") {
      const bool iso = j.text("to", "isotactic") == "isotactic";
      const int k = caps::set_tacticity(s, iso);
      if (k > 0 && j.num("clean", 1) != 0) caps::clean_up(s);
      what = std::string("Make ") + (iso ? "isotactic" : "syndiotactic") + " (" + std::to_string(k) + " centres inverted)";
    } else if (op == "attach") {
      const uint32_t t = uint32_t(j.num("target", -1));
      const auto at = caps::attach_fragment(s, t, j.text("smiles"), int(j.num("which", 0)), j.num("replace_h", 1) != 0);
      for (uint32_t a : at) added.push_back(double(a));
      if (j.num("clean", 1) != 0) {
        std::vector<char> m(s.atoms.size(), 0);
        for (uint32_t a : at) m[a] = 1;
        caps::clean_up(s, m, 0.5);
      }
      what = "Attach " + j.text("name", "fragment") + " to atom " + std::to_string(t + 1);
    } else if (op == "place") {
      caps::BuildOptions bo;
      bo.forcefield = "uff";
      const caps::System m = caps::build_molecule(j.text("smiles"), bo).system;
      // beside the structure: its extent's far side along x, 3 Å of space
      caps::Vec3 hi{-1e300, -1e300, -1e300}, lo{1e300, 1e300, 1e300}, mlo{1e300, 1e300, 1e300}, mhi{-1e300, -1e300, -1e300};
      for (const auto& a : s.atoms) for (int k = 0; k < 3; ++k) hi[k] = std::max(hi[k], a.pos[k]), lo[k] = std::min(lo[k], a.pos[k]);
      for (const auto& a : m.atoms) for (int k = 0; k < 3; ++k) mhi[k] = std::max(mhi[k], a.pos[k]), mlo[k] = std::min(mlo[k], a.pos[k]);
      const caps::Vec3 shift = s.atoms.empty() ? caps::Vec3{0, 0, 0}
                                               : caps::Vec3{hi[0] + 3.0 - mlo[0], (lo[1] + hi[1]) * 0.5 - (mlo[1] + mhi[1]) * 0.5, (lo[2] + hi[2]) * 0.5 - (mlo[2] + mhi[2]) * 0.5};
      int64_t mol = 0;
      for (const auto& a : s.atoms) mol = std::max(mol, a.mol);
      const uint32_t base = uint32_t(s.atoms.size());
      for (const auto& a0 : m.atoms) {
        caps::Atom a = a0;
        a.pos = a0.pos + shift;
        a.mol = mol + 1;
        a.id = s.atoms.empty() ? 1 : s.atoms.back().id + 1;
        a.resname = j.text("resname", "MOL").substr(0, 3);
        const std::string sym = caps::element(a.element).symbol;
        int type = 0;
        for (const auto& t : s.types) if (t.label == sym) type = t.type;
        if (type == 0) {
          for (const auto& t : s.types) type = std::max(type, t.type);
          ++type;
          caps::TypeInfo ti;
          ti.type = type, ti.mass = caps::element(a.element).mass, ti.label = sym;
          s.types.push_back(ti);
        }
        a.type = type;
        added.push_back(double(s.atoms.size()));
        s.atoms.push_back(a);
      }
      for (const auto& b : m.bonds) s.bonds.push_back({b.i + base, b.j + base, b.order});
      s.has_mol = true;
      what = "Place " + j.text("name", "molecule");
    } else if (op == "clean") {
      const auto at = atoms_of(d, j);
      std::vector<char> m;
      if (!at.empty()) { m.assign(s.atoms.size(), 0); for (uint32_t a : at) m[a] = 1; }
      caps::clean_up(s, m, j.num("ftol", 0.5));
      what = "Clean up with UFF" + (at.empty() ? std::string() : " (" + std::to_string(at.size()) + " atoms)");
    } else {
      throw std::invalid_argument("unknown edit '" + op + "'");
    }
    push_undo(d, what);
    store(d, s);
    r["ok"] = true;
    r["what"] = what;
    r["atoms"] = double(s.atoms.size());
    r["added"] = added;
  } catch (const std::exception& e) {
    r = caps::Json::object();
    r["ok"] = false;
    r["error"] = std::string(e.what());
  }
  return report_out(r.dump(0), out, cap);
}

extern "C" int32_t caps_undo(caps_doc* d, int32_t redo) {
  return guard([&] {
    auto& from = redo ? d->redo : d->undo;
    auto& to = redo ? d->undo : d->redo;
    if (from.empty()) throw std::runtime_error(redo ? "nothing to redo" : "nothing to undo");
    caps_doc::Snapshot now;
    now.topology = d->traj.topology;
    now.positions = d->traj.positions.at(d->current);
    now.what = from.back().what;
    to.push_back(std::move(now));
    d->traj.topology = from.back().topology;
    d->traj.positions[d->current] = from.back().positions;
    from.pop_back();
    d->field.reset();
    d->selection.assign(d->traj.topology.atoms.size(), 0);
    refresh(d);
    return 0;
  });
}

extern "C" int32_t caps_history(caps_doc* d, char* json, int32_t cap) {
  caps::Json j = caps::Json::object();
  caps::Json u = caps::Json::array(), r = caps::Json::array();
  for (const auto& x : d->undo) u.push_back(x.what);
  for (const auto& x : d->redo) r.push_back(x.what);
  j["undo"] = u, j["redo"] = r;
  return report_out(j.dump(0), json, cap);
}

extern "C" int32_t caps_select(caps_doc* d, const char* json, char* out, int32_t cap) {
  caps::Json r = caps::Json::object();
  try {
    const caps::Json j = caps::Json::parse(json && *json ? json : "{}");
    const caps::System& s = d->frame;
    const size_t n = s.atoms.size();
    if (d->selection.size() != n) d->selection.assign(n, 0);
    const std::string mode = j.text("mode");
    std::vector<char> m(n, 0);
    if (mode == "smarts") m = caps::select_smarts(s, j.text("pattern"));
    else if (mode == "element") m = caps::select_element(s, j.text("pattern"));
    else if (mode == "type") m = caps::select_type(s, j.text("pattern"));
    else if (mode == "charge") m = caps::select_charge(s, j.num("lo", -1e9), j.num("hi", 1e9));
    else if (mode == "within") m = caps::select_within(s, d->selection, j.num("distance", 3.0));
    else if (mode == "grow") m = caps::select_grow(s, d->selection, int(j.num("steps", 1)));
    else if (mode == "all") m.assign(n, 1);
    else if (mode == "none") {}
    else if (mode == "molecule") {
      const auto mol = s.molecules();
      std::set<int> want;
      for (uint32_t a : atoms_of(d, j)) want.insert(mol[a]);
      for (size_t i = 0; i < n; ++i) m[i] = want.count(mol[i]) ? 1 : 0;
    } else if (mode == "indices") {
      for (uint32_t a : atoms_of(d, j)) m[a] = 1;
    } else if (mode == "expression") {
      const caps::PipelineState st = caps::run_pipeline(s, caps::Pipeline{});
      const auto v = caps::evaluate_expression(st, j.text("pattern"));
      for (size_t i = 0; i < n && i < v.size(); ++i) m[i] = v[i] != 0;
    } else {
      throw std::invalid_argument("unknown selection mode '" + mode + "'");
    }
    const std::string op = j.text("op", "replace");
    auto& sel = d->selection;
    for (size_t i = 0; i < n; ++i) {
      if (op == "replace" || mode == "within" || mode == "grow") sel[i] = m[i];
      else if (op == "add") sel[i] = sel[i] || m[i];
      else if (op == "subtract") sel[i] = sel[i] && !m[i];
      else if (op == "intersect") sel[i] = sel[i] && m[i];
      else if (op == "invert") sel[i] = !sel[i];
    }
    r["ok"] = true;
    r["count"] = double(std::count(sel.begin(), sel.end(), 1));
    r["matched"] = double(std::count(m.begin(), m.end(), 1));
  } catch (const std::exception& e) {
    r = caps::Json::object();
    r["ok"] = false;
    r["error"] = std::string(e.what());
  }
  return report_out(r.dump(0), out, cap);
}

extern "C" int32_t caps_selection(caps_doc* d, char* json, int32_t cap) {
  caps::Json r = caps::Json::object();
  caps::Json idx = caps::Json::array();
  size_t count = 0;
  for (size_t i = 0; i < d->selection.size(); ++i)
    if (d->selection[i]) { ++count; if (idx.size() < 200000) idx.push_back(double(i)); }
  r["count"] = double(count);
  r["indices"] = idx;
  return report_out(r.dump(0), json, cap);
}

extern "C" int32_t caps_tacticity(caps_doc* d, char* json, int32_t cap) {
  const auto T = caps::tacticity(d->frame);
  caps::Json r = caps::Json::object();
  r["label"] = T.label;
  r["m"] = double(T.m), r["r"] = double(T.r), r["mm"] = double(T.mm), r["mr"] = double(T.mr), r["rr"] = double(T.rr);
  caps::Json chains = caps::Json::array();
  size_t centres = 0;
  for (const auto& c : T.chains) {
    caps::Json x = caps::Json::object();
    caps::Json cs = caps::Json::array();
    for (uint32_t a : c.centres) cs.push_back(double(a));
    x["centres"] = cs;
    x["dyads"] = c.dyads;
    chains.push_back(x);
    centres += c.centres.size();
  }
  r["chains"] = chains;
  r["centres"] = double(centres);
  return report_out(r.dump(0), json, cap);
}

extern "C" int32_t caps_element_number(const char* symbol) { return symbol ? caps::element_from_symbol(symbol) : 0; }

extern "C" int32_t caps_element_info(int32_t z, double* mass, double* covalent, double* vdw, uint32_t* rgb) {
  if (z <= 0 || z > caps::max_element()) return -1;
  const auto& e = caps::element(z);
  if (mass) *mass = e.mass;
  if (covalent) *covalent = e.covalent;
  if (vdw) *vdw = e.vdw;
  if (rgb) *rgb = caps::element_colour(z);
  return 0;
}

// ---------------------------------------------------------------- interactions and checks (v20)

extern "C" int32_t caps_interactions(caps_doc* d, const char* options_json, char* json, int32_t cap) {
  caps::Json j = caps::Json::object();
  try {
    const caps::Json o = caps::Json::parse(options_json && *options_json ? options_json : "{}");
    caps::InteractionOptions io;
    io.hb_distance = o.num("hb_distance", 3.5);
    io.hb_angle = o.num("hb_angle", 30);
    io.contact_margin = o.num("contact_margin", 0.4);
    io.clash_factor = o.num("clash_factor", 0.75);
    auto flag = [&](const char* k, bool def) { return !o.has(k) ? def : o[k].kind() == caps::Json::Bool ? o[k].boolean() : o[k].number() != 0; };
    const caps::System& f = d->frame;
    const auto R = caps::find_interactions(f, io);
    // the lines in the view
    d->checks.clear();
    auto seg = [&](uint32_t a, uint32_t b, unsigned rgb, double radius, int dashes) {
      caps::Vec3 A = f.atoms[a].pos, B = f.atoms[b].pos;
      if (f.cell.valid()) B = A + f.cell.minimum_image(B - A);
      if (dashes <= 1) { caps::Segment sg; sg.a = A, sg.b = B, sg.rgb = rgb, sg.radius = radius; d->checks.push_back(sg); return; }
      for (int k = 0; k < dashes; ++k) {
        caps::Segment sg;
        sg.a = A + (B - A) * (double(2 * k) / (2 * dashes - 1));
        sg.b = A + (B - A) * (double(2 * k + 1) / (2 * dashes - 1));
        sg.rgb = rgb, sg.radius = radius;
        d->checks.push_back(sg);
      }
    };
    const size_t cap_lines = 20000;
    if (flag("show_hbonds", true)) for (const auto& h : R.hbonds) if (d->checks.size() < cap_lines) seg(h.hydrogen, h.acceptor, 0x6CC4D8, 0.07, 4);
    if (flag("show_contacts", false)) for (const auto& c : R.contacts) if (d->checks.size() < cap_lines) seg(c.i, c.j, 0xF5A524, 0.05, 1);
    if (flag("show_clashes", true)) for (const auto& c : R.clashes) if (d->checks.size() < cap_lines) seg(c.i, c.j, 0xFF7B72, 0.13, 1);
    j["ok"] = true;
    j["hbonds"] = double(R.n_hbonds), j["contacts"] = double(R.n_contacts), j["clashes"] = double(R.n_clashes);
    j["molecules"] = double(R.molecules), j["net_charge"] = R.net_charge;
    caps::Json hb = caps::Json::array();
    for (size_t k = 0; k < R.hbonds.size() && k < 2000; ++k) {
      caps::Json x = caps::Json::object();
      x["donor"] = double(R.hbonds[k].donor), x["hydrogen"] = double(R.hbonds[k].hydrogen), x["acceptor"] = double(R.hbonds[k].acceptor);
      x["distance"] = R.hbonds[k].distance, x["angle"] = R.hbonds[k].angle;
      hb.push_back(x);
    }
    j["hbond_list"] = hb;
    caps::Json cl = caps::Json::array();
    for (size_t k = 0; k < R.clashes.size() && k < 2000; ++k) {
      caps::Json x = caps::Json::object();
      x["i"] = double(R.clashes[k].i), x["j"] = double(R.clashes[k].j), x["distance"] = R.clashes[k].distance;
      cl.push_back(x);
    }
    j["clash_list"] = cl;
    caps::Json issues = caps::Json::array();
    for (const auto& is : R.issues) {
      caps::Json x = caps::Json::object();
      x["level"] = is.level, x["title"] = is.title, x["detail"] = is.detail, x["fix"] = is.fix;
      caps::Json at = caps::Json::array();
      for (uint32_t a : is.atoms) at.push_back(double(a));
      x["atoms"] = at;
      issues.push_back(x);
    }
    j["issues"] = issues;
  } catch (const std::exception& e) {
    j = caps::Json::object();
    j["ok"] = false;
    j["error"] = std::string(e.what());
  }
  return report_out(j.dump(0), json, cap);
}

extern "C" void caps_clear_checks(caps_doc* d) { d->checks.clear(); }

// ---------------------------------------------------------------- import (design/boards/ImportDialog)
namespace {
caps::ImportOptions import_options(const char* options_json) {
  const caps::Json o = caps::Json::parse(options_json && *options_json ? options_json : "{}");
  auto flag = [&](const char* k, bool def) { return !o.has(k) ? def : o[k].kind() == caps::Json::Bool ? o[k].boolean() : o[k].number() != 0; };
  caps::ImportOptions r;
  const std::string b = o.text("bonds", "perceive");
  r.bonds = b == "file" ? caps::ImportOptions::FromFile : b == "none" ? caps::ImportOptions::None : caps::ImportOptions::Perceive;
  r.tolerance = o.num("tolerance", 0.45);
  r.bond_orders = flag("bond_orders", true);
  r.split = flag("split", true);
  r.unwrap = flag("unwrap", true);
  r.use_cell = flag("use_cell", true);
  return r;
}
}  // namespace

extern "C" caps_doc* caps_import(const char* path, const char* topology_path, const char* options_json) {
  try {
    auto* d = new caps_doc;
    d->traj = caps::import_file(path, topology_path ? topology_path : "", import_options(options_json));
    refresh(d);
    return d;
  } catch (const std::exception& e) {
    g_error = e.what();
  } catch (...) {
    g_error = "unknown error";
  }
  return nullptr;
}

extern "C" int32_t caps_import_preview(const char* path, const char* options_json, char* json, int32_t cap) {
  caps::Json j = caps::Json::object();
  try {
    const auto r = caps::import_preview(path, import_options(options_json));
    j["ok"] = true;
    j["format"] = r.file.format;
    j["format_name"] = r.file.format_name;
    j["units"] = r.file.units;
    j["bytes"] = double(r.file.bytes);
    caps::Json head = caps::Json::array();
    for (const auto& l : r.file.head) head.push_back(l);
    j["head"] = std::move(head);
    j["atoms"] = double(r.atoms);
    j["bonds_in_file"] = double(r.bonds_in_file);
    j["bonds"] = double(r.bonds);
    j["molecules"] = double(r.molecules);
    j["single"] = r.single; j["double"] = r.dbl; j["triple"] = r.triple; j["aromatic"] = r.aromatic;
    j["cell"] = r.cell;
    j["fragment_heavy"] = r.fragment_heavy;
    j["fragment_atoms"] = double(r.fragment.atoms.size());
    j["fragment_bonds"] = double(r.fragment.bonds.size());
    caps::Json notes = caps::Json::array();
    for (const auto& n : r.notes) notes.push_back(n);
    j["notes"] = std::move(notes);
  } catch (const std::exception& e) {
    j = caps::Json::object();
    j["ok"] = false;
    j["error"] = std::string(e.what());
  }
  return report_out(j.dump(0), json, cap);
}

extern "C" caps_doc* caps_import_fragment(const char* path, const char* options_json) {
  try {
    return doc_of(caps::import_preview(path, import_options(options_json)).fragment);
  } catch (const std::exception& e) {
    g_error = e.what();
  } catch (...) {
    g_error = "unknown error";
  }
  return nullptr;
}

// ---------------------------------------------------------------- export image / movie (design/boards/ExportDialog)
namespace {
std::string iso_now() {
  const std::time_t t = std::time(nullptr);
  char buf[32];
  std::strftime(buf, sizeof buf, "%Y-%m-%dT%H:%M:%SZ", std::gmtime(&t));
  return buf;
}

// The provenance manifest embedded in exported images: what was shown, from which file (and its sha256), how.
std::string image_manifest(caps_doc* d, const caps_camera* cam, const caps_render_opts* opt, const caps::Json& o) {
  caps::Json m = caps::Json::object();
  m["schema"] = "caps-image/1.0";
  m["generator"] = std::string("CAPS 0.1.0");
  m["created"] = iso_now();
  const std::string src = o.text("source");
  if (!src.empty()) {
    m["source"] = std::filesystem::path(src).filename().string();
    std::ifstream f(src, std::ios::binary);
    if (f) {
      std::string bytes((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
      m["source_sha256"] = caps::sha256_hex(bytes);
    }
  }
  m["frame"] = double(d->current);
  m["frames"] = double(d->traj.frames());
  m["atoms"] = double(shown(d).atoms.size());
  caps::Json c = caps::Json::object();
  c["yaw"] = cam->yaw; c["pitch"] = cam->pitch; c["zoom"] = cam->zoom; c["pan_x"] = cam->pan_x; c["pan_y"] = cam->pan_y; c["perspective"] = cam->perspective != 0;
  m["camera"] = std::move(c);
  caps::Json r = caps::Json::object();
  r["width"] = opt->width; r["height"] = opt->height; r["supersample"] = opt->supersample;
  r["style"] = opt->style; r["colour_by"] = opt->colour_by; r["background"] = opt->background;
  m["render"] = std::move(r);
  return m.dump(0);
}
}  // namespace

extern "C" int32_t caps_export_image(caps_doc* d, const caps_camera* cam, const caps_render_opts* opt, const char* path, const char* options_json,
                                     const uint8_t* overlay) {
  return guard([&] {
    const caps::Json o = caps::Json::parse(options_json && *options_json ? options_json : "{}");
    auto flag = [&](const char* k, bool def) { return !o.has(k) ? def : o[k].kind() == caps::Json::Bool ? o[k].boolean() : o[k].number() != 0; };
    caps::PngOptions p;
    p.sixteen = o.num("bits", 8) >= 16;
    p.dpi = o.num("dpi", 0);
    p.srgb = o.text("colour_profile", "srgb") != "none";
    if (flag("provenance", true)) p.text.push_back({"caps:provenance", image_manifest(d, cam, opt, o)});
    caps::RenderOptions ro = opts_of(d, opt);
    ro.deep = p.sixteen;
    caps::Renderer r;
    caps::Image img = r.render(shown(d), cam_of(cam), ro);
    if (overlay) {
      // labels and measurements drawn by the caller (straight alpha), laid over the image at full precision
      for (size_t k = 0; k < size_t(img.width) * img.height; ++k) {
        const uint8_t* o8 = overlay + 4 * k;
        if (o8[3] == 0) continue;
        const float a = o8[3] / 255.f;
        uint8_t* p8 = &img.rgba[4 * k];
        const float da = p8[3] / 255.f, oa = a + da * (1 - a);
        for (int c = 0; c < 3; ++c) {
          const float v = oa > 0 ? (o8[c] / 255.f * a + p8[c] / 255.f * da * (1 - a)) / oa : 0;
          p8[c] = uint8_t(std::clamp(v, 0.f, 1.f) * 255 + .5f);
          if (!img.rgba16.empty()) {
            const float d16 = img.rgba16[4 * k + c] / 65535.f;
            const float v16 = oa > 0 ? (o8[c] / 255.f * a + d16 * da * (1 - a)) / oa : 0;
            img.rgba16[4 * k + c] = uint16_t(std::clamp(v16, 0.f, 1.f) * 65535 + .5f);
          }
        }
        p8[3] = uint8_t(oa * 255 + .5f);
        if (!img.rgba16.empty()) img.rgba16[4 * k + 3] = uint16_t(oa * 65535 + .5f);
      }
    }
    caps::write_png(img, path, p);
    return 0;
  });
}

extern "C" int32_t caps_png_text(const char* path, char* json, int32_t cap) {
  caps::Json j = caps::Json::object();
  try {
    for (const auto& [k, v] : caps::read_png_text(path)) j[k] = v;
  } catch (const std::exception& e) {
    g_error = e.what();
    return -1;
  }
  return report_out(j.dump(0), json, cap);
}

extern "C" int32_t caps_export_movie(caps_doc* d, const caps_camera* cam, const caps_render_opts* opt, const char* path, const char* options_json,
                                     caps_series_progress_fn progress, void* user) {
  return guard([&] {
    const caps::Json o = caps::Json::parse(options_json && *options_json ? options_json : "{}");
    const std::string format = o.text("format", "apng");
    const int nf = int(d->traj.frames());
    const int from = std::clamp(int(o.num("from", 0)), 0, nf - 1);
    const int to = std::clamp(int(o.num("to", nf - 1)), from, nf - 1);
    const int step = std::max(1, int(o.num("step", 1)));
    const double turn = o.num("turntable", 0);   // degrees of yaw over the whole movie
    const bool turntable = turn != 0;   // the current frame turned about the vertical
    int count = turntable ? std::max(2, int(o.num("turntable_frames", 120))) : (to - from) / step + 1;
    const int fps = std::clamp(int(o.num("fps", 24)), 1, 120);
    caps::PngOptions p;
    p.srgb = o.text("colour_profile", "srgb") != "none";
    caps::Json manifest_opts = o;
    std::unique_ptr<caps::ApngWriter> apng;
    if (format == "apng") {
      if (o.text("provenance", "1") != "0") p.text.push_back({"caps:provenance", image_manifest(d, cam, opt, o)});
      apng = std::make_unique<caps::ApngWriter>(path, fps, int(o.num("loops", 0)), p);
    } else {
      std::filesystem::create_directories(path);
    }
    const size_t keep = d->current;
    caps::Renderer r;
    const caps::RenderOptions ro = opts_of(d, opt);
    int written = 0;
    try {
      for (int k = 0; k < count; ++k) {
        const int f = turntable ? int(keep) : from + k * step;
        if (size_t(f) != d->current) { d->current = size_t(f); refresh(d); }
        caps::Camera c = cam_of(cam);
        if (turn != 0) c.yaw += turn * M_PI / 180.0 * double(k) / double(count);
        caps::Image img = r.render(shown(d), c, ro);
        if (apng) apng->add(img);
        else {
          char name[32];
          std::snprintf(name, sizeof name, "frame_%05d.png", k + 1);
          caps::write_png(img, (std::filesystem::path(path) / name).string(), p);
        }
        ++written;
        if (progress && progress(written, count, user) != 0) break;
      }
    } catch (...) {
      if (d->current != keep) { d->current = keep; refresh(d); }
      throw;
    }
    if (apng) apng->close();
    if (d->current != keep) { d->current = keep; refresh(d); }
    return written;
  });
}
