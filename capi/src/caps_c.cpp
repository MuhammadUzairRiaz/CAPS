#include "caps_c.h"

#include <cstdio>
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <filesystem>
#include <thread>
#include <fstream>
#include <string>

#include "caps/analysis.hpp"
#include "caps/dynamics.hpp"
#include "caps/superpose.hpp"
#include "caps/elements.hpp"
#include "caps/example_typing.hpp"
#include "caps/equilibrate.hpp"
#include "caps/ffdef.hpp"
#include "caps/grow.hpp"
#include "caps/io.hpp"
#include "caps/mechanics.hpp"
#include "caps/molecule.hpp"
#include "caps/config.hpp"
#include "caps/bench.hpp"
#include "caps/polymer.hpp"
#include "caps/chimd.hpp"
#include "caps/chipair.hpp"
#include "caps/polystats.hpp"
#include "caps/resolution.hpp"
#include "caps/kspace.hpp"
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
#include "caps/provenance.hpp"
#include "caps/recipe.hpp"
#include "caps/colourvision.hpp"
#include "caps/query.hpp"
#include "caps/charges.hpp"
#include "caps/molinfo.hpp"
#include "caps/yaml.hpp"
#include "caps/voids.hpp"
#include "caps/kremer_grest.hpp"
#include "caps/nano.hpp"
#include "caps/json.hpp"

#include <map>
#include <mutex>
#include <chrono>
#include <memory>
#include <set>
#include <numeric>
#include <functional>
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
  caps::FFDef fill;                          // borrowed rules used only where the base (and extra) define nothing
  std::vector<std::string> imported;         // files the imported rules came from
  std::map<int32_t, std::string> overrides;  // atom → type set by hand
  std::string charges = "types";             // types (force field), gasteiger, keep (from the file)
  bool auto_charges = false;                 // automatic: the force field's charges, Gasteiger when its types carry none
  caps::TypingResult typing;
  std::vector<std::string> types;
  caps::ParamReport rep;
  std::shared_ptr<const caps::ForceField> ff;   // null while atoms are untyped
  bool complete = false;
  std::vector<std::string> prep_notes;       // what was done to the structure for the force field (united-atom sites)
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
  std::vector<caps::RelaxOptions::Restraint> restraints;   // distance restraints for caps_relax
  std::vector<caps::RelaxOptions::DihedralRestraint> dihedral_restraints;   // and dihedral ones
  double ph = -1;         // Add hydrogens: residues protonated at this pH (< 0: neutral valences)
  std::unique_ptr<caps::Pipeline> pipeline;        // caps_pipeline_set: steps run on every shown frame
  std::unique_ptr<caps::PipelineState> pstate;     // its result for the current frame
  std::vector<int32_t> shown_of;                   // frame index → first shown particle (−1: deleted)
  std::array<int, 3> cell_repeats{1, 1, 1};        // caps_crystal_build: a supercell of this many unit cells
  std::vector<caps::Segment> overlay;              // caps_peptide_build with ribbon: tubes drawn with the atoms (no pipeline)
  caps::Scene scene;                               // caps_render_scene: the arrays the caller reads until the next call
  caps_live_fn live_fn = nullptr;             // caps_set_live: snapshots of a running caps_md / caps_equilibrate
  void* live_user = nullptr;
  AppearanceState look;                            // caps_set_appearance
  int smooth_window = 1;                           // caps_set_smoothing: frames averaged for display
  std::vector<std::vector<caps::Vec3>> scan_frames; // caps_torsion_scan: the geometry of each point
  std::vector<caps::Vec3> scan_original;           // the frame before a scan point was shown
  std::vector<char> selection;                     // caps_select: the Studio's selection over the frame's atoms
  struct Snapshot { caps::System topology; std::vector<caps::Vec3> positions; std::string what; };
  std::vector<Snapshot> undo, redo;                // caps_edit history
  std::vector<std::vector<Snapshot>> branches;     // redo steps set aside when an edit followed an undo (newest last)
  struct Named { std::string name; Snapshot state; int step = 0; bool on_branch = false; };
  std::vector<Named> snapshots;                    // caps_snapshot: named states to go back to or compare
  std::vector<caps::Segment> checks;               // caps_interactions: H-bonds, contacts, clashes drawn in the view
  caps::Manifest prov;                             // provenance: the steps that produced this structure
  std::vector<caps::VoidSphere> voids;             // caps_voids: the largest empty spheres of the frame
  struct Checkpoint {                              // the last checkpoint of an MD / equilibration run (caps_checkpoint)
    bool has = false;
    std::string kind, ended, error;
    std::vector<caps::Vec3> x, v;
    caps::Cell cell;
    int64_t step = 0, steps = 0;
    double dt = 1;
  } checkpoint;
  std::unique_ptr<caps::Mesh> void_mesh;           // … drawn translucent when shown
  std::array<int, 3> images{1, 1, 1};              // caps_set_images: periodic images drawn around the cell (faded)
  float image_fade = 0.7f;
  int save_wrap = 0;                               // caps_set_save_wrap: 0 as shown, 1 atoms into the cell, 2 molecule centres
  std::vector<double> atom_values;                 // caps_set_atom_values: a per-atom quantity coloured on a ramp (SASA …)
  int atom_values_ramp = 0;
  // Backbone style (design/boards/DisplayStyles): per atom 1 main chain, 2 heavy atom of a molecule without one, 0 hidden;
  // recomputed when the atom or bond count changes
  std::vector<uint8_t> bb_mask, bb_mask_p;   // the frame's, and the pipeline result's (keyed by the result object)
  const void* bb_p = nullptr;
  size_t bb_atoms = SIZE_MAX, bb_bonds = SIZE_MAX;
  int bb_chains = 0, bb_atoms_on = 0;
  // Display options (design/boards/DisplayStyles, LensView): view only, the structure is untouched
  struct Display {
    bool polar_h_only = false;     // hide hydrogens on carbon
    bool selection_full = false;   // selected atoms keep all their atoms (H) whatever the style
    bool lens = false;             // all-atom lens: inside one style, outside another
    int lens_centre = -1;          // atom index
    double lens_radius = 10;
    int lens_inside = 0, lens_outside = 4;
    bool lens_dim = false;
    // clip slab (design/boards/Appearance "Clip planes"): only atoms whose fractional coordinate along the axis (the cell's,
    // else the structure's extent) lies in [clip_from, clip_to] are drawn; clip_invert draws the rest instead
    bool clip = false;
    int clip_axis = 2;
    double clip_from = 0.0, clip_to = 0.5;
    bool clip_invert = false;
  } display;
  int vision = 0;                                  // caps_set_vision: the view as seen with a colour-vision deficiency
  double vision_severity = 1.0;
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

// A running MD's positions as a new document about four times a second (wall clock), for a live view: the callee owns
// it. Wrapped into the cell, as the Grow live view; stats: step, time, and the density of the cell.
std::function<void(const caps::EnergyTerms&, const std::vector<double>&, const caps::Cell&, int64_t)> live_hook(caps_doc* d, const caps::System& topo, double dt_fs) {
  if (!d->live_fn) return {};
  auto fn = d->live_fn;
  auto user = d->live_user;
  auto last = std::make_shared<std::chrono::steady_clock::time_point>();   // the first step shows at once
  auto base = std::make_shared<caps::System>(topo);
  double mass = 0;
  for (const auto& a : topo.atoms) mass += caps::element(a.element).mass;
  return [fn, user, last, base, mass, dt_fs](const caps::EnergyTerms&, const std::vector<double>& x, const caps::Cell& c, int64_t step) {
    const auto now = std::chrono::steady_clock::now();
    if (now - *last < std::chrono::milliseconds(250)) return;
    *last = now;
    auto* sd = new caps_doc;
    sd->traj.topology = *base;
    std::vector<caps::Vec3> pp(x.size() / 3);
    for (size_t i = 0; i < pp.size(); ++i) pp[i] = {x[3 * i], x[3 * i + 1], x[3 * i + 2]};
    sd->traj.topology.cell = c;
    sd->traj.positions.push_back(std::move(pp));
    sd->traj.cells.push_back(c);
    sd->traj.timesteps.push_back(step);
    sd->wrap = c.valid();
    refresh(sd);
    caps::Json j = caps::Json::object();
    j["step"] = double(step), j["time_ps"] = double(step) * dt_fs / 1000.0, j["atoms"] = double(base->atoms.size());
    if (c.valid()) j["density"] = mass / c.volume() * 1.66053906660;
    fn(sd, j.dump(0).c_str(), user);
  };
}

caps::Camera cam_of(const caps_camera* c) {
  caps::Camera k;
  if (!c) return k;
  k.yaw = c->yaw; k.pitch = c->pitch; k.zoom = c->zoom > 0 ? c->zoom : 1; k.pan_x = c->pan_x; k.pan_y = c->pan_y; k.perspective = c->perspective != 0;
  return k;
}

// The backbone mask for the frame (cached on the document; mutable cache behind a const document).
const std::vector<uint8_t>& backbone_mask(const caps_doc* dc) {
  auto* d = const_cast<caps_doc*>(dc);
  const auto& s = d->frame;
  if (d->bb_atoms == s.atoms.size() && d->bb_bonds == s.bonds.size()) return d->bb_mask;
  d->bb_atoms = s.atoms.size(), d->bb_bonds = s.bonds.size();
  d->bb_mask.assign(s.atoms.size(), 0);
  const auto bb = caps::backbones(s, 4);
  int nm = 0;
  const auto mol = s.molecules(&nm);
  std::vector<char> has(size_t(std::max(nm, 1)), 0);
  for (const auto& path : bb)
    for (uint32_t a : path) d->bb_mask[a] = 1, has[size_t(mol[a])] = 1;
  for (size_t i = 0; i < s.atoms.size(); ++i)
    if (!has[size_t(mol[i])] && s.atoms[i].element != 1) d->bb_mask[i] = 2;
  d->bb_chains = int(bb.size());
  d->bb_atoms_on = 0;
  for (const auto& path : bb) d->bb_atoms_on += int(path.size());
  return d->bb_mask;
}

// The same for the pipeline's particles (replicas and deletions change which atoms there are).
const std::vector<uint8_t>& backbone_mask_pipeline(const caps_doc* dc) {
  auto* d = const_cast<caps_doc*>(dc);
  const auto* key = d->pstate.get();
  if (d->bb_p == key && d->bb_mask_p.size() == d->pstate->system.atoms.size()) return d->bb_mask_p;
  const auto& s = d->pstate->system;
  d->bb_p = key;
  d->bb_mask_p.assign(s.atoms.size(), 0);
  int nm = 0;
  const auto mol = s.molecules(&nm);
  std::vector<char> has(size_t(std::max(nm, 1)), 0);
  for (const auto& path : caps::backbones(s, 4))
    for (uint32_t a : path) d->bb_mask_p[a] = 1, has[size_t(mol[a])] = 1;
  for (size_t i = 0; i < s.atoms.size(); ++i)
    if (!has[size_t(mol[i])] && s.atoms[i].element != 1) d->bb_mask_p[i] = 2;
  return d->bb_mask_p;
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
  r.lod_near = o->lod_near > 0 ? o->lod_near : 0;
  r.lod_far = o->lod_far > 0 ? o->lod_far : 0;
  if (!d->pstate && d->selection.size() == d->frame.atoms.size()) {   // the selection ringed (up to 50 000 atoms)
    for (size_t i = 0; i < d->selection.size() && r.highlight.size() < 50000; ++i) if (d->selection[i]) r.highlight.push_back(int(i));
  }
  if (!d->pstate) { r.segments = d->overlay; r.segments.insert(r.segments.end(), d->checks.begin(), d->checks.end()); }
  if (d->void_mesh) r.meshes.push_back({d->void_mesh.get(), 0x4FB3D9, 0.32f});
  if (d->pstate && r.style == caps::Style::Backbone) {   // tubes through the pipeline's chains too
    const auto& m = backbone_mask_pipeline(d);
    r.atom_style.resize(m.size());
    for (size_t i = 0; i < m.size(); ++i)
      r.atom_style[i] = uint8_t(m[i] == 1 ? caps::Style::Backbone : m[i] == 2 ? caps::Style::NoHydrogens : caps::Style::Hidden);
  }
  const auto& D = d->display;
  const bool appearance_styles = d->look.active && d->look.style.size() == d->frame.atoms.size();
  if (!d->pstate && !appearance_styles && (r.style == caps::Style::Backbone || D.polar_h_only || D.lens || (D.selection_full && !d->selection.empty()))) {
    const auto& A = d->frame.atoms;
    const size_t n = A.size();
    std::vector<std::vector<uint32_t>> host;   // hydrogens' heavy atom (polar H, selections)
    auto heavy_of = [&](size_t i) -> int {
      if (host.empty()) {
        host.assign(n, {});
        for (const auto& b : d->frame.bonds) host[b.i].push_back(b.j), host[b.j].push_back(b.i);
      }
      return host[i].empty() ? -1 : int(host[i][0]);
    };
    // one style for an atom under a global style: Backbone uses the main-chain mask, the others hide H themselves
    auto style_of = [&](caps::Style g, size_t i) -> uint8_t {
      if (g == caps::Style::Backbone) {
        const auto& m = backbone_mask(d);
        return uint8_t(m[i] == 1 ? caps::Style::Backbone : m[i] == 2 ? caps::Style::NoHydrogens : caps::Style::Hidden);
      }
      if (D.polar_h_only && A[i].element == 1 && g != caps::Style::NoHydrogens) {
        const int h = heavy_of(i);
        if (h >= 0 && A[size_t(h)].element == 6) return uint8_t(caps::Style::Hidden);
      }
      return uint8_t(g);
    };
    r.atom_style.resize(n);
    for (size_t i = 0; i < n; ++i) r.atom_style[i] = style_of(r.style, i);
    if (D.lens && D.lens_centre >= 0 && size_t(D.lens_centre) < n) {
      const caps::Vec3 c = A[size_t(D.lens_centre)].pos;
      const bool periodic = d->frame.cell.valid();
      const double r2 = D.lens_radius * D.lens_radius;
      if (D.lens_dim) r.faded.assign(n, 0), r.fade = 0.55f;
      for (size_t i = 0; i < n; ++i) {
        caps::Vec3 dv = A[i].pos - c;
        if (periodic) dv = d->frame.cell.minimum_image(dv);
        const bool in = caps::dot(dv, dv) <= r2;
        r.atom_style[i] = style_of(static_cast<caps::Style>(std::clamp(in ? D.lens_inside : D.lens_outside, 0, 4)), i);
        if (D.lens_dim && !in) r.faded[i] = 1;
      }
    }
    if (D.selection_full && d->selection.size() == n)   // selected atoms and their hydrogens drawn in full
      for (size_t i = 0; i < n; ++i) {
        const bool sel = d->selection[i] || (A[i].element == 1 && heavy_of(i) >= 0 && d->selection[size_t(heavy_of(i))]);
        if (sel) r.atom_style[i] = uint8_t(caps::Style::BallAndStick);
      }
  }
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
  if (!d->pstate && !d->atom_values.empty() && d->atom_values.size() == d->frame.atoms.size()) {
    r.colour_by = caps::ColourBy::Property;
    r.property = d->atom_values;
    r.ramp = static_cast<caps::Ramp>(std::clamp(d->atom_values_ramp, 0, 2));
    r.symmetric = false;
  }
  if (d->pstate) {
    const auto& st = *d->pstate;
    r.colours = st.colour;
    for (size_t i = 0; i < r.colours.size(); ++i) if (st.selected[i]) r.colours[i] = 0xE5484D;   // selected particles in red
    r.segments = st.segments;
    // the pipeline's Radius (Å) and Transparency (0 … 1) properties, drawn
    if (auto it = st.props.find("Radius"); it != st.props.end() && it->second.size() == st.system.atoms.size())
      r.radius.assign(it->second.begin(), it->second.end());
    if (auto it = st.props.find("Transparency"); it != st.props.end() && it->second.size() == st.system.atoms.size())
      r.transparency.assign(it->second.begin(), it->second.end());
    if (r.colour_by == caps::ColourBy::Property) caps::property_values(st, "DistanceToCOM", r.property);
  } else if (r.colour_by == caps::ColourBy::Property && r.property.size() != d->frame.atoms.size()) {
    r.property = d->dcom;
  }
  if (D.clip) {   // the clip slab hides what lies outside it, whatever the styles
    const caps::System& sys = d->pstate ? d->pstate->system : d->frame;
    const size_t n = sys.atoms.size();
    if (r.atom_style.size() != n) r.atom_style.assign(n, uint8_t(r.style));
    const int ax = std::clamp(D.clip_axis, 0, 2);
    double lo = 0, hi = 1;
    if (!sys.cell.valid()) {
      lo = 1e300, hi = -1e300;
      for (const auto& a : sys.atoms) lo = std::min(lo, a.pos[size_t(ax)]), hi = std::max(hi, a.pos[size_t(ax)]);
      if (hi - lo < 1e-9) hi = lo + 1;
    }
    for (size_t i = 0; i < n; ++i) {
      double f;
      if (sys.cell.valid()) { f = sys.cell.to_fractional(sys.atoms[i].pos)[size_t(ax)]; f -= std::floor(f); }
      else f = (sys.atoms[i].pos[size_t(ax)] - lo) / (hi - lo);
      const bool in = f >= D.clip_from && f <= D.clip_to;
      if (in == D.clip_invert) r.atom_style[i] = uint8_t(caps::Style::Hidden);
    }
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

// an OPLS file whose types are classes (CT, CA …) has no OPLS charges; the files with OPLS's own numbered types do
std::string opls_hint(const std::string& ff_name) {
  if (ff_name.find("OPLS") == std::string::npos) return "";
  return ". OPLS-AA's own charges come with its numbered types: choose OPLS-AA (2024 parameter file) or OPLS-AA (BOSS 4.8, 2008) to use them";
}

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
  // gap fillers first: the last matching rule wins, so the force field's own rules (and the imported ones) come later
  caps::prepend_fill(def, F.fill);
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
  // a Materials Studio .car carries its force-field types (IFF's inorganic ones, ClayFF's …): where they are types of this
  // force field they are kept, the rules type the rest
  int car_kept = 0;
  if (!uff && d->traj.topology.source_format == "car" && F.file_types.size() == n) {
    std::set<std::string> ffnames;
    for (const auto& t : def.types) ffnames.insert(t.name);
    for (size_t i = 0; i < n; ++i)
      if (ffnames.count(F.file_types[i].second)) {
        F.types[i] = F.file_types[i].second;
        if (i < F.typing.why.size()) F.typing.why[i] = "type from the .car file";
        if (i < F.typing.types.size()) F.typing.types[i] = F.types[i];
        ++car_kept;
      }
  }
  for (const auto& [i, t] : F.overrides)
    if (i >= 0 && size_t(i) < n) F.types[i] = t;
  std::set<std::string> known;
  for (const auto& t : def.types) known.insert(t.name);
  int untyped = 0;
  for (auto& t : F.types)
    if (t.empty() || !known.count(t)) { t.clear(); ++untyped; }
  F.rep = caps::ParamReport{};
  if (car_kept) F.rep.notes.push_back(std::to_string(car_kept) + " atoms keep the force-field types of the .car file; the typing rules gave the others");
  F.ff.reset();
  if (!untyped && uff) {
    caps::UffOptions uo;
    uo.keep_charges = F.charges == "keep";
    uo.qeq = F.charges == "qeq";
    uo.labels = F.types;
    F.ff = std::make_shared<caps::ForceField>(caps::assign_uff(s, uo));
    for (const auto& note : F.ff->notes) F.rep.notes.push_back(note);
  } else if (!untyped) {
    if (F.auto_charges) {
      // the force field's own charges when its types carry them; libraries that keep charges on molecule templates
      // (DL_FIELD's OPLS-AA, GAFF, CHARMM) have none per type, and then Gasteiger–Marsili, said in the report
      try {
        F.ff = std::make_shared<caps::ForceField>(caps::parameterize(s, def, F.types, "types", &F.rep, true));
        F.charges = "types";
        // a bond without an increment (pcff.frc has none for an alkoxysilane's o-sio): not the force field's charges
        for (const auto& m : F.rep.missing)
          if (m.rfind("bond increment", 0) == 0) throw caps::FFError(def.name + " has no charge for type pair " + m.substr(15) + " (no bond increment)");
        // fixed per-type charges that miss the formal charge (a group OPLS-AA's charges do not close): the companion
        // force field's bond-increment charges when it names one (OPLS-AA 2024 → OPLS 2005), said in the report
        if (!def.charge_increments_from.empty()) {
          double net = 0;
          for (double q : F.ff->charge) net += q;
          int formal = 0;
          for (int c : caps::perceive(s).charge) formal += c;
          if (std::fabs(net - formal) > 1e-3) {
            try {
              std::string note;
              caps::System sq = s;
              const auto q = caps::companion_charges(s, def, F.ff_path, &note);
              for (size_t i = 0; i < n; ++i) sq.atoms[i].charge = q[i];
              sq.has_charges = true;
              caps::ParamReport rep2;
              auto ff2 = std::make_shared<caps::ForceField>(caps::parameterize(sq, def, F.types, "keep", &rep2, true));
              char b[200];
              std::snprintf(b, sizeof b, "%s's own charges add up to %+.3f e here, not %+d e: ", def.name.c_str(), net, formal);
              rep2.notes.push_back(b + note + " (charges: automatic)");
              F.ff = std::move(ff2);
              F.rep = std::move(rep2);
              F.charges = "increments";
            } catch (const caps::FFError&) {}   // the companion cannot type this structure: the imbalance is reported below
          }
        }
      } catch (const caps::FFError& e) {
        if (std::string(e.what()).find("has no charge for type") == std::string::npos) throw;
        F.rep = caps::ParamReport{};
        const std::string why = std::string(e.what()).substr(std::string(e.what()).find("type"));
        try {
          F.ff = std::make_shared<caps::ForceField>(caps::parameterize(s, def, F.types, "gasteiger", &F.rep, true));
          F.charges = "gasteiger";
          F.rep.notes.push_back(def.name + " carries no charges on its atom types (" + why + "): Gasteiger–Marsili charges were used instead (charges: automatic)" +
                                opls_hint(def.name));
        } catch (const std::exception& g) {   // Gasteiger–Marsili has no parameters for some groups (S=O, most metals): QEq
          F.rep = caps::ParamReport{};
          F.ff = std::make_shared<caps::ForceField>(caps::parameterize(s, def, F.types, "qeq", &F.rep, true));
          F.charges = "qeq";
          F.rep.notes.push_back(def.name + " carries no charges on its atom types (" + why + ") and Gasteiger–Marsili has none for this structure (" + g.what() +
                                "): QEq charges were used instead (charges: automatic)" + opls_hint(def.name));
        }
      }
    } else if (F.charges == "increments") {
      // another force field's bond-increment charges (OPLS 2005's for OPLS-AA 2024) with this one's types
      std::string note;
      caps::System sq = s;
      const auto q = caps::companion_charges(s, def, F.ff_path, &note);
      for (size_t i = 0; i < n; ++i) sq.atoms[i].charge = q[i];
      sq.has_charges = true;
      F.ff = std::make_shared<caps::ForceField>(caps::parameterize(sq, def, F.types, "keep", &F.rep, true));
      F.rep.notes.push_back(note);
    } else {
      F.ff = std::make_shared<caps::ForceField>(caps::parameterize(s, def, F.types, F.charges, &F.rep, true));
    }
  }
  // the physics check: charges taken from the force field's own types must leave the structure at its formal charge
  // (OPLS-AA's group charges balance group by group; bond increments always do). A mismatch means a group was typed
  // only partly, so the assignment is incomplete, as a missing parameter would make it.
  if (F.ff && F.charges == "types" && !uff) {
    double net = 0;
    for (double q : F.ff->charge) net += q;
    int formal = 0;
    for (int c : caps::perceive(s).charge) formal += c;
    if (std::fabs(net - formal) > 1e-3) {
      char b[240];
      std::snprintf(b, sizeof b, "charges: the force field's charges add up to %+.3f e, not the formal charge %+d e (a group typed only partly)", net, formal);
      F.rep.missing.push_back(b);
    }
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
    if (i < F.rep.charge_keys.size()) a["ck"] = F.rep.charge_keys[i];   // the charge key its increments were looked up by
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
  int estimated = 0, imported = 0, filled = 0;
  caps::Json filled_terms = caps::Json::array();
  for (const auto& [k, v] : F.rep.used) {
    const auto sp = k.find(' ');
    const std::string name = sp == std::string::npos ? k : k.substr(sp + 1);
    if (name.rfind("user:", 0) == 0) estimated += v;
    if (name.rfind("imported:", 0) == 0) imported += v;
    if (name.rfind("filled:", 0) == 0) {
      filled += v;
      filled_terms.push_back(k + " × " + std::to_string(v));
    }
  }
  r["filled"] = double(filled);
  r["filled_terms"] = filled_terms;
  estimated += F.rep.estimated_terms;   // by analogy (the force field's "analogies"), each listed below
  r["estimated"] = double(estimated);
  r["imported"] = double(imported);
  caps::Json analog = caps::Json::array();
  for (const auto& x : F.rep.estimated) analog.push_back(x);
  r["by_analogy"] = analog;
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
    if (!std::isnan(t.charge)) u["q"] = t.charge;
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
  for (const auto& x : F.prep_notes) notes.push_back(x);
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

// ---------------------------------------------------------------- provenance helpers
std::string g6(double x) {
  char b[32];
  std::snprintf(b, sizeof b, "%.6g", x);
  return b;
}

// sha256 of a file (files over 1 GB are not hashed: their size stands in).
std::string file_sha256(const std::string& path) {
  std::error_code ec;
  const auto n = std::filesystem::file_size(path, ec);
  if (ec) return "";
  if (n > (1ull << 30)) return "not hashed (" + std::to_string(n) + " bytes)";
  std::ifstream f(path, std::ios::binary);
  std::string bytes((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
  return caps::sha256_hex(bytes);
}

void prov_input(caps_doc* d, const std::string& path) {
  if (path.empty()) return;
  const std::string name = std::filesystem::path(path).filename().string();
  for (const auto& [n, h] : d->prov.inputs) if (n == name) return;
  d->prov.inputs.push_back({name, file_sha256(path)});
}

void prov_step(caps_doc* d, const std::string& engine, const std::string& summary, caps::KeyValues params, const std::string& rng = "",
               std::vector<std::string> cites = {}, caps::KeyValues approx = {}) {
  caps::ProvStep s;
  s.engine = engine;
  s.summary = summary;
  s.params = std::move(params);
  s.rng = rng;
  s.cites = std::move(cites);
  s.approximations = std::move(approx);
  s.time = caps::now_iso();
  d->prov.steps.push_back(std::move(s));
}

std::string seeded(uint64_t seed) { return "mt19937-64 · seed " + std::to_string(seed); }

// A JSON options object as ordered parameters (nested values as compact JSON).
caps::KeyValues json_params(const caps::Json& j) {
  caps::KeyValues out;
  if (!j.is_object()) return out;
  for (const auto& [k, v] : j.members()) out.push_back({k, v.is_string() ? v.str() : v.is_number() ? g6(v.number()) : v.dump(0)});
  return out;
}
caps::KeyValues json_params(const char* text) {
  try { return json_params(caps::Json::parse(text && *text ? text : "{}")); } catch (...) { return {}; }
}

// The modelling choices of an energy evaluation, and the papers behind them.
// Worker threads an evaluation uses (0 asks for the default: the configured maximum, else the hardware's, at most 16).
int threads_used(int requested) {
  if (requested > 0) return requested;
  if (const int n = caps::max_threads(); n > 0) return n;
  const unsigned h = std::thread::hardware_concurrency();
  return int(std::clamp(h == 0 ? 1u : h, 1u, 16u));
}

caps::KeyValues energy_approx(double cutoff, bool coulomb, bool tail, int threads) {
  caps::KeyValues a;
  a.push_back({"van der Waals", "cut-off " + g6(cutoff) + " Å" + (tail ? " + tail correction" : "")});
  a.push_back({"Electrostatics", !coulomb ? "off" : g_elec.mode == 1 ? "SPME · relative tolerance " + g6(g_elec.rtol) : "damped shifted force · cut-off " + g6(cutoff) + " Å"});
  a.push_back({"Constraints", "none"});
  // forces are summed in worker order: the same thread count gives the same trajectory bit for bit
  a.push_back({"Precision", "double · reproducible with " + std::to_string(threads_used(threads)) + " threads"});
  return a;
}
void elec_cites(std::vector<std::string>& c, bool coulomb) {
  if (coulomb) c.push_back(g_elec.mode == 1 ? "essmann1995" : "fennell2006");
}

// A document just read from a file: its manifest from the sidecar, else the file as the first input.
std::string ff_label(const caps_doc* d) {
  if (!d->field) return "built-in default (GAFF for C and H, UFF otherwise)";
  const std::string name = d->field->base.name.empty() ? std::filesystem::path(d->field->ff_path).stem().string() : d->field->base.name;
  return name + (d->field->complete ? "" : " (incomplete)");
}

void prov_opened(caps_doc* d, const std::string& path, const std::string& topology, const std::string& how, caps::KeyValues params = {},
                 std::vector<std::string> cites = {}) {
  if (auto m = caps::read_manifest(path)) {
    d->prov = std::move(*m);
    return;
  }
  prov_input(d, path);
  if (!topology.empty()) prov_input(d, topology);
  const auto& t = d->traj.topology;
  params.insert(params.begin(), {"file", std::filesystem::path(path).filename().string()});
  params.push_back({"atoms", std::to_string(t.atoms.size())});
  params.push_back({"frames", std::to_string(d->traj.frames())});
  params.push_back({"bonds", std::to_string(t.bonds.size()) + (t.bonds_from_file ? " from the file" : " perceived")});
  prov_step(d, how, "read " + std::filesystem::path(path).filename().string(), std::move(params), "", std::move(cites));
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
    prov_opened(d, path, topology_path ? topology_path : "", "io.read");
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
    prov_opened(d, path, topology_path ? topology_path : "", "io.read");
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

caps_doc* caps_shadow(caps_doc* d) {
  if (!d) return nullptr;
  try {
    auto* sd = new caps_doc;
    sd->traj.topology = d->traj.topology;
    const size_t k = d->current;
    if (k < d->traj.positions.size()) sd->traj.positions.push_back(d->traj.positions[k]);
    if (k < d->traj.cells.size()) sd->traj.cells.push_back(d->traj.cells[k]);
    if (k < d->traj.timesteps.size()) sd->traj.timesteps.push_back(d->traj.timesteps[k]);
    sd->wrap = d->wrap;
    sd->look.active = d->look.active;   // the look's settings; its surface meshes are not carried
    sd->look.layers = d->look.layers;
    sd->look.colour = d->look.colour;
    sd->look.ramp = d->look.ramp;
    sd->selection = d->selection;
    sd->cell_repeats = d->cell_repeats;
    sd->ph = d->ph;
    sd->held_mol = d->held_mol;
    sd->restraints = d->restraints;
    sd->dihedral_restraints = d->dihedral_restraints;
    sd->analysis = d->analysis;
    sd->eq_checks = d->eq_checks;
    if (d->field) sd->field = std::make_unique<FieldState>(*d->field);
    refresh(sd);
    return sd;
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
    prov_step(d, "grow.trials", "polystyrene chains grown in a periodic cell, best-of-k trial placement by contact margin",
              {{"chains", std::to_string(o->chains)}, {"DP", std::to_string(o->dp)}, {"tacticity", o->tacticity == 1 ? "isotactic" : o->tacticity == 2 ? "syndiotactic" : "atactic"},
               {o->box > 0 ? "box" : "density", o->box > 0 ? g6(o->box) + " Å" : g6(o->density) + " g/cm³"}, {"contact scale", g6(g.contact_scale)}},
              seeded(o->seed), {"matsumoto1998"});
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

int32_t save_frame(caps_doc* d, const std::string& p);

int32_t caps_save(caps_doc* d, const char* path) {
  return guard([&] {
    const std::string p = path;
    if (d->save_wrap && d->frame.cell.valid()) {   // wrap on save: a copy, the document keeps its display
      caps::System w = d->frame;
      const auto mol = w.molecules();
      std::vector<caps::Vec3> shift(w.atoms.size(), caps::Vec3{0, 0, 0});
      if (d->save_wrap == 2) {   // molecule centres into the cell, molecules kept whole
        std::map<int, std::pair<caps::Vec3, int>> cen;
        for (size_t i = 0; i < w.atoms.size(); ++i) { auto& c = cen[mol[i]]; c.first = c.first + w.atoms[i].pos; c.second++; }
        std::map<int, caps::Vec3> sh;
        for (auto& [m, c] : cen) {
          const caps::Vec3 centre = c.first * (1.0 / c.second);
          caps::Vec3 f = w.cell.to_fractional(centre);
          for (int k = 0; k < 3; ++k) f[k] = std::floor(f[k]);
          sh[m] = w.cell.to_cartesian(f) - w.cell.to_cartesian(caps::Vec3{0, 0, 0});
        }
        for (size_t i = 0; i < w.atoms.size(); ++i) w.atoms[i].pos = w.atoms[i].pos - sh[mol[i]];
      } else {
        for (auto& a : w.atoms) {
          caps::Vec3 f = w.cell.to_fractional(a.pos);
          for (int k = 0; k < 3; ++k) f[k] -= std::floor(f[k]);
          a.pos = w.cell.to_cartesian(f);
        }
      }
      const caps::System keep = d->frame;
      d->frame = w;
      struct Restore { caps_doc* d; caps::System s; ~Restore() { d->frame = std::move(s); } } restore{d, keep};
      return save_frame(d, p);
    }
    return save_frame(d, p);
  });
}

int32_t save_frame(caps_doc* d, const std::string& p) {
  return guard([&] {
    auto ends = [&](const char* e) { const std::string x = e; return p.size() >= x.size() && p.compare(p.size() - x.size(), x.size(), x) == 0; };
    if (ends(".pdb")) caps::write_pdb(d->frame, p);
    else if (ends(".xyz")) caps::write_xyz(d->frame, p);
    else if (ends(".mol2")) caps::write_mol2(d->frame, p);
    else if (ends(".car")) caps::write_car(d->frame, p);   // Materials Studio, with its .mdf; the force-field types as names
    else if (ends(".gro")) caps::write_gro(d->frame, p);
    else if (ends(".sdf") || ends(".mol")) caps::write_sdf(d->frame, p);
    else if (ends(".cif")) caps::write_cif(d->frame, p);
    else if (d->field) {   // the Field assignment: its coefficients when complete, else the structure alone
      if (d->field->complete) caps::write_lammps_data_or_structure(d->frame, *d->field->ff, elec(), p);
      else caps::write_lammps_data(d->frame, p);
    } else {
      caps::ForceField ff;
      bool typed = true;
      try { ff = default_ff(d->frame); } catch (const caps::FieldError&) { typed = false; }
      if (typed) caps::write_lammps_data_ff(d->frame, ff, elec(), p);
      else caps::write_lammps_data(d->frame, p);
    }
    if (!d->prov.steps.empty()) {
      try { caps::write_manifest(d->prov, p); } catch (...) {}   // the structure is written; the sidecar is a bonus
    }
    return 0;
  });
}

int32_t caps_lammps_shake(caps_doc* d, int32_t mode, const char* group, char* text, int32_t cap) {
  return guard([&] {
    caps::ForceField ff;
    if (d->field && d->field->complete) ff = *d->field->ff;
    else ff = default_ff(d->frame);
    const std::string out = caps::lammps_shake_fix(d->frame, ff, elec(), static_cast<caps::ConstraintMode>(std::clamp(mode, 0, 2)),
                                                   group && *group ? group : "all");
    const int32_t need = int32_t(out.size() + 1);
    if (text && cap > 0) {
      const size_t m = std::min<size_t>(size_t(cap - 1), out.size());
      std::memcpy(text, out.data(), m);
      text[m] = 0;
    }
    return need;
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

int32_t caps_gromacs(caps_doc* d, const char* stem, char* text, int32_t cap) {
  return guard([&] {
    caps::ForceField ff;
    if (d->field && d->field->complete) ff = *d->field->ff;
    else ff = default_ff(d->frame);
    std::string out;
    std::vector<std::string> notes;
    notes = stem && *stem ? caps::write_gromacs(d->frame, ff, elec(), stem) : caps::gromacs_notes(d->frame, ff, elec());
    if (d->held_mol > 0) notes.push_back("the held molecule is not frozen here: give it an index group and freezegrps / freezedim");
    for (const auto& n : notes) out += "; note: " + n + "\n";
    out += caps::gromacs_mdp(d->frame, ff, elec());
    const int32_t need = int32_t(out.size() + 1);
    if (text && cap > 0) {
      const size_t m = std::min<size_t>(size_t(cap - 1), out.size());
      std::memcpy(text, out.data(), m);
      text[m] = 0;
    }
    return need;
  });
}

// Export center (ABI 23): the simulation files for LAMMPS and GROMACS in one call, from a complete force field.
extern "C" int32_t caps_export_engines(caps_doc* d, const char* dir, const char* options, char* out, int32_t cap) {
  caps::Json r = caps::Json::object();
  try {
    const caps::Json o = caps::Json::parse(options && *options ? options : "{}");
    auto flag = [&](const char* k, bool def) { return o.has(k) && o[k].kind() == caps::Json::Bool ? o[k].boolean() : def; };
    const bool lammps = flag("lammps", true), gromacs = flag("gromacs", true), preview = flag("preview", false);
    const std::string stem = o.text("stem", "system");
    if (stem.empty() || stem.find('/') != std::string::npos || stem.find('\\') != std::string::npos) throw std::runtime_error("the file stem must be a plain name");
    if (!d->field || !d->field->ff) throw std::runtime_error("assign a force field first (Force field step)");
    if (!d->field->complete) throw std::runtime_error("the force field is incomplete for this structure: see the Force field step for the untyped atoms and missing terms");
    const caps::ForceField& ff = *d->field->ff;
    const caps::System& s = d->frame;
    caps::LammpsRun run;
    const std::string kind = o.text("run", "check");
    run.kind = kind == "none" ? caps::LammpsRun::Kind::None : kind == "minimize" ? caps::LammpsRun::Kind::Minimize
             : kind == "nvt" ? caps::LammpsRun::Kind::NVT : kind == "npt" ? caps::LammpsRun::Kind::NPT : caps::LammpsRun::Kind::Check;
    run.minimize_first = flag("minimize_first", true);
    run.temperature = o.num("temperature", 300);
    run.pressure = o.num("pressure", 1.0);
    run.dt = o.num("dt", 0);   // 0: the force field's own
    run.steps = int64_t(o.num("steps", 100000));
    run.thermo_every = int(o.num("thermo_every", 1000));
    run.dump_every = int(o.num("dump_every", 5000));
    run.seed = uint64_t(o.num("seed", 4928459));
    run.constraints = caps::constraints_from_string(o.text("constraints", "none"));
    if (run.temperature <= 0 || run.dt < 0 || run.steps < 0) throw std::runtime_error("temperature and time step must be positive");
    namespace fs = std::filesystem;
    const fs::path folder = preview ? fs::temp_directory_path() / ("caps_export_" + std::to_string(reinterpret_cast<uintptr_t>(d))) : fs::path(dir ? dir : "");
    if (folder.empty()) throw std::runtime_error("choose a folder");
    fs::create_directories(folder);
    const auto base = (folder / stem).string();
    caps::Json files = caps::Json::array(), notes = caps::Json::array();
    std::vector<std::pair<std::string, std::string>> written;   // name, what
    const caps::EnergyOptions e = elec();
    // LAMMPS styles: the force field's own (default) or CAPS-exact; hybrid; long-range sum; cut-off
    caps::LammpsStyle ls;
    ls.native = o.text("lammps_styles", "native") != "exact";
    ls.hybrid = flag("hybrid", false);
    ls.coulomb = o.text("coulomb", "auto");
    ls.kspace_accuracy = o.num("kspace_accuracy", 1e-4);
    ls.cutoff = o.num("cutoff", 0);
    if (o.has("tail") && o["tail"].kind() == caps::Json::Bool) ls.tail = o["tail"].boolean() ? 1 : 0;
    if (ls.kspace_accuracy <= 0 || ls.cutoff < 0) throw std::runtime_error("the k-space accuracy and the cut-off must be positive");
    // a force field LAMMPS cannot express (GROMOS's reaction field …) refuses the LAMMPS files only: the GROMACS files
    // are still written, and the reason goes back as lammps_error
    std::string lammps_error;
    if (lammps) try {
      std::vector<std::string> lnotes;
      caps::write_lammps_data_ff(s, ff, e, base + ".data", false, ls);
      caps::write_lammps_input(s, ff, e, stem + ".data", base + ".in", d->held_mol, true, run, ls, &lnotes);
      for (const auto& n : lnotes) notes.push_back(caps::Json("LAMMPS: " + n));
      written.push_back({stem + ".data", "atoms, bonds, masses and bonded coefficients"});
      written.push_back({stem + ".in", "styles, every pair_coeff and the run"});
    } catch (const std::exception& ex) {
      lammps_error = ex.what();   // refused before its files (prepare checks the force field first)
      if (!gromacs) throw;
    }
    // a force field GROMACS cannot express (class II's 9-6 Lennard-Jones and cross terms …) refuses the GROMACS files
    // only: the LAMMPS files are still written, and the reason goes back as gromacs_error
    std::string gromacs_error;
    if (gromacs) try {
      for (const auto& n : caps::write_gromacs(s, ff, e, base)) notes.push_back("GROMACS: " + n);
      if (d->held_mol > 0) notes.push_back("GROMACS: the held molecule is not frozen: give it an index group and freezegrps / freezedim");
      if (run.constraints == caps::ConstraintMode::HBonds && (run.kind == caps::LammpsRun::Kind::NVT || run.kind == caps::LammpsRun::Kind::NPT)) {
        const auto nb = s.neighbours();
        size_t waters = 0;
        for (size_t o = 0; o < s.atoms.size(); ++o)
          if (s.atoms[o].element == 8 && nb[o].size() == 2 && s.atoms[nb[o][0]].element == 1 && s.atoms[nb[o][1]].element == 1) ++waters;
        if (waters)
          notes.push_back("GROMACS: constraints = h-bonds holds water's O–H bonds but not its angle; CAPS and LAMMPS hold " + std::to_string(waters) +
                          " water(s) rigid — use [ settles ] in the topology for the same model");
      }
      // the core's .mdp is a single point with the matching non-bonded settings; a protocol replaces its run lines
      if (run.kind != caps::LammpsRun::Kind::Check) {
        std::ifstream in(base + ".mdp");
        std::string mdp((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
        in.close();
        std::string kept;
        const bool md_run = run.kind == caps::LammpsRun::Kind::NVT || run.kind == caps::LammpsRun::Kind::NPT;
        std::istringstream ls(mdp);
        for (std::string line; std::getline(ls, line);) {
          const auto key = line.substr(0, line.find_first_of(" =\t"));
          if (key == "integrator" || key == "nsteps" || key == "dt" || key == "nstcalcenergy" || key == "nstenergy") continue;
          if (key == "constraints" && run.constraints != caps::ConstraintMode::None && md_run) continue;
          if (line.rfind("; GROMACS run parameters", 0) == 0) continue;
          kept += line + "\n";
        }
        char b[512];
        const bool md = run.kind == caps::LammpsRun::Kind::NVT || run.kind == caps::LammpsRun::Kind::NPT;
        std::string head = "; GROMACS run parameters written by CAPS: " + ff.name + "\n";
        if (!md) {
          std::snprintf(b, sizeof b, "integrator               = steep\nnsteps                   = %lld\nemtol                    = 10.0        ; kJ/mol/nm\n", static_cast<long long>(std::max<int64_t>(run.steps, 5000)));
          head += b;
        } else {
          std::snprintf(b, sizeof b, "integrator               = md\ndt                       = %.6g\nnsteps                   = %lld\nnstxout-compressed       = %d\nnstenergy                = %d\nnstlog                   = %d\n",
                        caps::lammps_timestep(run, ff) / 1000, static_cast<long long>(run.steps), run.dump_every, run.thermo_every, run.thermo_every);
          head += b;
        }
        std::string tail;
        if (md) {
          std::snprintf(b, sizeof b, "tcoupl                   = V-rescale   ; Bussi\ntc-grps                  = System\ntau-t                    = 0.1\nref-t                    = %.6g\n", run.temperature);
          tail += b;
          if (run.kind == caps::LammpsRun::Kind::NPT) {
            std::snprintf(b, sizeof b, "pcoupl                   = C-rescale\npcoupltype               = isotropic\ntau-p                    = 1.0\nref-p                    = %.6g\ncompressibility          = 4.5e-5\n", run.pressure * 1.01325);
            tail += b;
          }
          std::snprintf(b, sizeof b, "gen-vel                  = yes\ngen-temp                 = %.6g\ngen-seed                 = %llu\n", run.temperature, static_cast<unsigned long long>(run.seed % 2147483647));
          tail += b;
          if (run.constraints != caps::ConstraintMode::None)
            tail += run.constraints == caps::ConstraintMode::AllBonds ? "constraints              = all-bonds   ; as in CAPS (LINCS)\n"
                                                                     : "constraints              = h-bonds     ; bonds to hydrogen, as in CAPS (LINCS)\n";
        }
        std::ofstream mo(base + ".mdp");
        mo << head << kept << tail;
        if (md && run.minimize_first) {
          std::ofstream em(folder / (stem + "_em.mdp"));
          std::snprintf(b, sizeof b, "integrator               = steep\nnsteps                   = 5000\nemtol                    = 10.0        ; kJ/mol/nm\n");
          em << "; energy minimisation before " << (run.kind == caps::LammpsRun::Kind::NPT ? "NPT" : "NVT") << ", written by CAPS: " << ff.name << "\n" << b << kept;
        }
      }
      const std::string itp = stem + ".itp";
      written.push_back({stem + ".top", "defaults, atom types, every pair"});
      if (fs::exists(folder / itp)) written.push_back({itp, "molecule types"});
      written.push_back({stem + ".gro", "coordinates in nm, molecules whole"});
      written.push_back({stem + ".mdp", run.kind == caps::LammpsRun::Kind::Check ? "single point, matching cut-offs" : "the run, matching cut-offs"});
      if (fs::exists(folder / (stem + "_em.mdp")) && run.minimize_first && (run.kind == caps::LammpsRun::Kind::NVT || run.kind == caps::LammpsRun::Kind::NPT))
        written.push_back({stem + "_em.mdp", "minimisation first (gmx grompp -f " + stem + "_em.mdp)"});
    } catch (const std::exception& ex) {
      gromacs_error = ex.what();   // refused before any GROMACS file is written (write_gromacs checks first)
      if (!lammps_error.empty()) throw std::runtime_error("LAMMPS: " + lammps_error + "; GROMACS: " + gromacs_error);
      if (!lammps) throw;
    }
    if (d->field->rep.estimated_terms) {
      std::string ex;
      for (size_t k = 0; k < d->field->rep.estimated.size() && k < 3; ++k) ex += (k ? "; " : "") + d->field->rep.estimated[k];
      notes.push_back(std::to_string(d->field->rep.estimated_terms) + " terms use the parameters of analogous types (estimated): " + ex +
                      (d->field->rep.estimated.size() > 3 ? " …" : ""));
    }
    const int head_lines = int(o.num("head_lines", 0));
    for (const auto& [name, what] : written) {
      caps::Json f = caps::Json::object();
      f["name"] = name, f["what"] = what;
      f["bytes"] = double(fs::file_size(folder / name));
      if (head_lines > 0) {
        caps::Json lines = caps::Json::array();
        std::ifstream in(folder / name);
        std::string line;
        for (int k = 0; k < head_lines && std::getline(in, line); ++k) lines.push_back(line);
        f["head"] = std::move(lines);
      }
      files.push_back(std::move(f));
    }
    if (preview) fs::remove_all(folder);
    // what was checked before writing
    double net = 0;
    for (double q : ff.charge) net += q;
    std::set<std::string> types(d->field->types.begin(), d->field->types.end());
    caps::Json checks = caps::Json::object();
    checks["atoms"] = double(s.atoms.size());
    checks["typed"] = double(s.atoms.size() - size_t(d->field->typing.untyped));
    checks["types"] = double(types.size());
    checks["type_pairs"] = double(ff.type_names.size() * (ff.type_names.size() + 1) / 2);
    // every form counted (class II bonds, angles, torsions and out-of-plane terms too)
    checks["bonds"] = double(ff.bonds.size() + ff.bonds2.size() + ff.bonds_x.size());
    checks["angles"] = double(ff.angles.size() + ff.angles2.size() + ff.angles_x.size());
    checks["dihedrals"] = double(ff.dihedrals.size() + ff.dihedrals2.size() + ff.cbt.size());
    checks["impropers"] = double(ff.impropers.size() + ff.impropers_harmonic.size() + ff.inversions.size() + ff.impropers2.size());
    checks["missing"] = double(d->field->rep.missing.size());
    checks["net_charge"] = net;
    checks["charges"] = d->field->auto_charges ? "automatic" : d->field->charges;
    checks["forcefield"] = ff.name;
    checks["density"] = s.density();
    r["ok"] = true;
    if (!gromacs_error.empty()) r["gromacs_error"] = gromacs_error;
    if (!lammps_error.empty()) r["lammps_error"] = lammps_error;
    r["files"] = std::move(files);
    r["notes"] = std::move(notes);
    r["checks"] = std::move(checks);
    r["folder"] = preview ? std::string() : folder.string();
  } catch (const std::exception& ex) {
    r["ok"] = false;
    r["error"] = std::string(ex.what());
  }
  return report_out(r.dump(), out, cap);
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
    r.pushoff_ramp_ps = std::max(0.0, o->pushoff_ramp_ps);
    if (o->pushoff_cap > 0) r.pushoff_cap = o->pushoff_cap;
    if (o->pushoff_temperature > 0) r.pushoff_temperature = o->pushoff_temperature;
    r.relax_box = o->relax_box != 0;
    r.box_anisotropic = o->box_anisotropic != 0;
    if (o->box_axes & 7)
      for (int k = 0; k < 3; ++k) r.box_axes[k] = (o->box_axes >> k) & 1;
    r.dihedral_restraints = d->dihedral_restraints;
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
    for (const auto& rs : d->restraints)
      if (rs.i < s.atoms.size() && rs.j < s.atoms.size() && rs.i != rs.j) r.restraints.push_back(rs);
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
    {
      static const char* names[] = {"sd", "cg", "lbfgs", "fire"};
      static const char* titles[] = {"steepest descent", "Polak–Ribière conjugate gradients", "L-BFGS", "FIRE"};
      const int mth = std::clamp(o->method, 0, 3);
      std::vector<std::string> c;
      if (mth == 1) c.push_back("polak1969");
      if (mth == 2) c.push_back("liu1989");
      if (mth == 3) c.push_back("bitzek2006");
      if (r.pushoff) c.push_back("auhl2003");
      elec_cites(c, r.energy.coulomb);
      caps::KeyValues pr = {{"minimiser", titles[mth]}, {"|F|max", g6(r.ftol) + " kcal/mol/Å"}, {"max iterations", std::to_string(r.max_iterations)},
                            {"push-off", !r.pushoff ? "off" : r.pushoff_ramp_ps > 0 ? "MD " + g6(r.pushoff_ramp_ps) + " ps at " + g6(r.pushoff_temperature) +
                                                                  " K, cap ramped to " + g6(r.pushoff_cap > 0 ? r.pushoff_cap : 500) + " kcal/mol/Å, then minimised"
                                                                : "on"},
                            {"force field", ff_label(d)}};
      if (r.target_density > 0) pr.push_back({"target density", g6(r.target_density) + " g/cm³"});
      if (r.relax_box) pr.push_back({"box relaxation", g6(r.pressure) + " atm"});
      if (d->held_mol > 0) pr.push_back({"held molecule", std::to_string(d->held_mol)});
      if (!r.restraints.empty()) {
        std::string t;
        for (size_t k = 0; k < r.restraints.size() && k < 4; ++k)
          t += (k ? "; " : "") + std::to_string(r.restraints[k].i + 1) + "–" + std::to_string(r.restraints[k].j + 1) + " to " +
               std::to_string(r.restraints[k].r0).substr(0, 5) + " Å (k " + std::to_string(r.restraints[k].k).substr(0, 5) + ")";
        pr.push_back({"distance restraints", t + (r.restraints.size() > 4 ? " …" : "")});
      }
      prov_step(d, std::string("relax.") + names[mth], std::string(titles[mth]) + (rep.converged ? ", converged" : ", stopped before the tolerance"), std::move(pr), "",
                std::move(c), energy_approx(r.energy.cutoff, r.energy.coulomb, false, r.energy.threads));
    }
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
    m.thermostat = static_cast<caps::Thermostat>(std::clamp(o->thermostat, 0, 3));
    if (o->tau_t > 0) m.tau_t = o->tau_t;
    m.barostat = static_cast<caps::Barostat>(std::clamp(o->barostat, 0, 3));
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
    m.respa = std::clamp(o->respa, 1, 16);
    m.constraints = static_cast<caps::ConstraintMode>(std::clamp(o->constraints, 0, 2));
    m.step_offset = std::max<int64_t>(0, o->step_offset);
    m.checkpoint_every = o->checkpoint_every > 0 ? o->checkpoint_every : o->checkpoint_every < 0 ? 0 : std::max<int64_t>(100, m.steps / 50);
    d->checkpoint = {};
    d->checkpoint.kind = "md";
    d->checkpoint.steps = m.step_offset + m.steps;
    d->checkpoint.dt = m.dt;
    m.checkpoint = [d](const std::vector<double>& x, const std::vector<double>& v, const caps::Cell& c, int64_t step) {
      auto& k = d->checkpoint;
      k.x.resize(x.size() / 3);
      k.v.resize(v.size() / 3);
      for (size_t i = 0; i < k.x.size(); ++i) k.x[i] = {x[3 * i], x[3 * i + 1], x[3 * i + 2]}, k.v[i] = {v[3 * i], v[3 * i + 1], v[3 * i + 2]};
      k.cell = c;
      k.step = step;
      k.has = true;
    };
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
    m.each_step = live_hook(d, s, m.dt);
    caps::DynamicsReport rep;
    try {
      caps::run_dynamics(s, m, &rep);
    } catch (const caps::DynamicsCancelled&) {
      d->checkpoint.ended = "stopped";
      throw;
    } catch (const std::exception& e) {
      d->checkpoint.ended = "failed";
      d->checkpoint.error = e.what();
      throw;
    }
    d->checkpoint.ended = "finished";
    {
      const bool nvt = m.thermostat != caps::Thermostat::None, npt = nvt && m.barostat != caps::Barostat::None;
      std::vector<std::string> c = {"swope1982"};
      if (m.thermostat == caps::Thermostat::NoseHoover) c.push_back("martyna1992");
      if (m.barostat == caps::Barostat::MTK) c.push_back("martyna1994");
      if (m.constraints != caps::ConstraintMode::None) c.push_back("ryckaert1977"), c.push_back("andersen1983");
      if (m.thermostat == caps::Thermostat::Bussi) c.push_back("bussi2007");
      if (npt && m.barostat == caps::Barostat::CRescale) c.push_back("bernetti2020");
      if (npt && m.barostat == caps::Barostat::Berendsen) c.push_back("berendsen1984");
      elec_cites(c, m.energy.coulomb);
      caps::KeyValues pr = {{"length", g6(m.dt * double(m.steps) / 1000.0) + " ps · " + std::to_string(m.steps) + " steps of " + g6(m.dt) + " fs"},
                            {"temperature", g6(m.temperature) + " K"}, {"thermostat", std::string(caps::to_string(m.thermostat)) + (nvt ? " · τ " + g6(m.tau_t) + " fs" : "")}};
      if (npt) pr.push_back({"barostat", std::string(caps::to_string(m.barostat)) + " · " + g6(m.pressure) + " atm · τ " + g6(m.tau_p) + " fs"});
      if (m.constraints != caps::ConstraintMode::None) pr.push_back({"constraints", caps::to_string(m.constraints)});
      pr.push_back({"force field", ff_label(d)});
      const bool drew = m.new_velocities || s.velocities.empty();
      prov_step(d, npt ? "dynamics.npt" : nvt ? "dynamics.nvt" : "dynamics.nve", npt ? "NPT molecular dynamics" : nvt ? "NVT molecular dynamics" : "NVE molecular dynamics",
                std::move(pr), drew || nvt ? seeded(m.seed) : "", std::move(c), energy_approx(m.energy.cutoff, m.energy.coulomb, m.energy.tail, m.energy.threads));
    }
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
    e.md.thermostat = o->thermostat == 2 ? caps::Thermostat::Langevin : o->thermostat == 3 ? caps::Thermostat::NoseHoover : caps::Thermostat::Bussi;
    e.md.barostat = o->barostat == 2 ? caps::Barostat::Berendsen : o->barostat == 3 ? caps::Barostat::MTK : caps::Barostat::CRescale;
    if (e.md.barostat == caps::Barostat::MTK) e.md.thermostat = caps::Thermostat::NoseHoover;   // MTK is the NPT half of Nosé–Hoover
    if (o->tau_t > 0) e.md.tau_t = o->tau_t;
    if (o->tau_p > 0) e.md.tau_p = o->tau_p;
    e.md.constraints = static_cast<caps::ConstraintMode>(std::clamp(o->constraints, 0, 2));
    d->checkpoint = {};
    d->checkpoint.kind = "equilibrate";
    d->checkpoint.dt = e.md.dt;
    e.md.checkpoint_every = 1000;
    e.md.checkpoint = [d](const std::vector<double>& x, const std::vector<double>& v, const caps::Cell& c, int64_t step) {
      auto& k = d->checkpoint;
      k.x.resize(x.size() / 3);
      k.v.resize(v.size() / 3);
      for (size_t i = 0; i < k.x.size(); ++i) k.x[i] = {x[3 * i], x[3 * i + 1], x[3 * i + 2]}, k.v[i] = {v[3 * i], v[3 * i + 1], v[3 * i + 2]};
      k.cell = c;
      k.step = step;
      k.has = true;
    };
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
    if (o->tol_internal > 0) e.tol_internal = o->tol_internal;
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
    e.md.each_step = live_hook(d, s, e.md.dt);
    caps::EquilibrateReport rep;
    try {
      caps::equilibrate(s, e, &rep);
    } catch (const caps::DynamicsCancelled&) {
      d->checkpoint.ended = "stopped";
      throw;
    } catch (const std::exception& ex) {
      d->checkpoint.ended = "failed";
      d->checkpoint.error = ex.what();
      throw;
    }
    d->checkpoint.ended = "finished";
    {
      const bool l21 = e.stages.size() == 21 && e.stages.back().label.find("final") != std::string::npos;
      double pmax = 0;
      for (const auto& st : e.stages) pmax = std::max(pmax, st.pressure);
      std::vector<std::string> c = {"swope1982", e.md.thermostat == caps::Thermostat::Bussi ? "bussi2007" : e.md.thermostat == caps::Thermostat::NoseHoover ? "martyna1992" : ""};
      if (e.md.constraints != caps::ConstraintMode::None) c.push_back("ryckaert1977"), c.push_back("andersen1983");
      c.push_back(e.md.barostat == caps::Barostat::CRescale ? "bernetti2020" : e.md.barostat == caps::Barostat::MTK ? "martyna1994" : "berendsen1984");
      if (l21) c.insert(c.begin(), "larsen2011");
      elec_cites(c, e.md.energy.coulomb);
      c.erase(std::remove(c.begin(), c.end(), std::string()), c.end());
      const std::string text = protocol ? protocol : "";
      prov_step(d, l21 ? "equilibrate.larsen21" : "equilibrate.protocol",
                (l21 ? "21-step compression and decompression" : std::to_string(e.stages.size()) + "-stage protocol") + std::string(rep.converged ? ", converged" : ""),
                {{"stages", std::to_string(e.stages.size())}, {"Pmax", g6(pmax) + " atm"}, {"length", g6(rep.ps) + " ps"},
                 {"protocol sha256", caps::sha256_hex(text).substr(0, 12)}, {"force field", ff_label(d)}},
                seeded(e.md.seed), std::move(c), energy_approx(e.md.energy.cutoff, e.md.energy.coulomb, e.md.energy.tail, e.md.energy.threads));
    }
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
    prov_step(d, "pack.lbfgs", std::to_string(rep.molecules) + " molecules packed without overlaps",
              {{"molecules", std::to_string(rep.molecules)}, {"atoms", std::to_string(rep.atoms)}, {"closest contact", g6(rep.dmin) + " Å"},
               {"input sha256", caps::sha256_hex(text ? text : "").substr(0, 12)}},
              "", {"martinez2009", "liu1989"});
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
    r.during_md = o->during_md != 0;
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
    {
      std::string names;
      for (const auto& t : r.templates) names += (names.empty() ? "" : ", ") + t.name;
      const double conv = rep.cycles.empty() ? 0 : rep.cycles.back().conversion;
      prov_step(d, "react.templates", std::to_string(rep.reactions) + " reactions · conversion " + g6(conv),
                {{"templates", names}, {"target conversion", g6(r.target_conversion)}, {"cycles", std::to_string(rep.cycles.size())},
                 {"relax between cycles", r.relax ? "yes" : "no"}, {"MD between cycles", g6(r.md_ps) + " ps"}},
                seeded(r.seed), {"matsumoto1998"});
    }
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

namespace { void push_undo(caps_doc* d, const std::string& what); }   // below, with the edits

int32_t caps_field_assign(caps_doc* d, const char* ff_path, const char* rules_path, int32_t charges) {
  return guard([&] {
    auto F = std::make_unique<FieldState>();
    F->ff_path = ff_path ? ff_path : "";
    F->base = caps::is_uff(F->ff_path) ? caps::uff_definition() : caps::load_forcefield(F->ff_path);
    if (rules_path && *rules_path && !caps::is_uff(F->ff_path)) {
      F->base.typing.clear();
      caps::load_typing(F->base, rules_path);
    }
    F->charges = charges == 1 ? "gasteiger" : charges == 2 ? "keep" : charges == 3 ? "qeq" : charges == 5 ? "increments" : "types";
    F->auto_charges = charges == 4;
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
    // a united-atom force field: the structure's hydrogens on carbon fold into their carbons first (undoable)
    if (caps::needs_prepare(F->base)) {
      caps::System s = d->frame;
      std::string ch = F->auto_charges ? "auto" : F->charges;
      const std::string note = caps::prepare_for_forcefield(s, F->base, ch);
      if (!note.empty()) {
        push_undo(d, (F->base.united_atom ? "United-atom for " : "Prepared for ") + F->base.name);
        std::vector<caps::Vec3> pos;
        for (const auto& at : s.atoms) pos.push_back(at.pos);
        s.has_charges = true;
        d->traj.topology = s;
        d->traj.positions = {pos};
        d->traj.cells = {s.cell};
        d->traj.timesteps = {0};
        d->current = 0;
        refresh(d);
        F->file_types.clear();
        F->file_charges.clear();
        for (const auto& at : s.atoms) F->file_types.push_back({at.type, at.name}), F->file_charges.push_back(at.charge);
        F->file_type_table = s.types;
        F->file_has_charges = true;
        F->prep_notes.push_back(note);
      }
      // the model's own charges (a Martini protein's, a molecule template's), whether or not the structure changed
      if (ch == "keep") F->charges = "keep", F->auto_charges = false;
    }
    if (F->charges == "keep" && !d->traj.topology.has_charges && !F->file_has_charges)
      throw caps::FFError("the structure has no charges to keep; use the force field's or Gasteiger charges");
    // "keep" means the file's charges, not the previous assignment's
    for (size_t i = 0; i < d->traj.topology.atoms.size() && i < F->file_charges.size(); ++i) d->traj.topology.atoms[i].charge = F->file_charges[i];
    d->field = std::move(F);
    refresh(d);
    field_run(d);
    {
      const std::string ffp = d->field->ff_path;
      const bool uff = caps::is_uff(ffp);
      std::vector<std::string> c;
      if (uff) c.push_back("rappe1992");
      else if (ffp.find("gaff") != std::string::npos) c.push_back("wang2004");
      if (d->field->charges == "qeq") c.push_back("rappe1991");
      if (d->field->charges == "gasteiger") c.push_back("gasteiger1980");
      if (!uff) prov_input(d, ffp);
      const std::string ch = d->field->charges == "types" ? "from the force field" : d->field->charges == "keep" ? "kept from the file" : d->field->charges == "qeq" ? "QEq" : "Gasteiger–Marsili";
      prov_step(d, "field.assign", ff_label(d) + " · charges " + ch,
                {{"force field", uff ? "UFF" : std::filesystem::path(ffp).filename().string()}, {"charges", ch}, {"typed", d->field->complete ? "every atom" : "incomplete"}},
                "", std::move(c), {{"Estimated parameters", d->field->complete ? "none" : "some atoms untyped (see the typing report)"}});
    }
    return d->field->complete ? 0 : 1;
  });
}

// Which library force fields can describe the current structure (caps_field_coverage): each one's typing tried, then
// its parameters looked up, with the untyped atoms grouped by chemical environment and the net charge its own charges
// give. Force-field definitions are cached (the library does not change while the Studio runs).
namespace {
std::string atom_env(const caps::System& s, const caps::Perception& p, uint32_t i) {
  std::vector<std::string> nb;
  for (size_t k = 0; k < p.nb[i].size(); ++k) {
    const uint32_t j = p.nb[i][k];
    const int o = p.order[i][k];
    nb.push_back(std::string(o == 2 ? "=" : o == 3 ? "#" : "") + caps::element(s.atoms[j].element).symbol);
  }
  std::sort(nb.begin(), nb.end());
  std::string r = caps::element(s.atoms[i].element).symbol;
  if (p.aromatic[i]) r += " aromatic";
  r += nb.empty() ? " (no bonds)" : " bonded to ";
  for (size_t k = 0; k < nb.size(); ++k) r += (k ? " " : "") + nb[k];
  return r;
}
}  // namespace

extern "C" int32_t caps_field_coverage(caps_doc* d, const char* dir, caps_stage_fn progress, void* user, char* out, int32_t cap) {
  caps::Json r = caps::Json::object();
  try {
    static std::mutex cache_lock;
    static std::map<std::string, std::shared_ptr<caps::FFDef>> cache;
    const caps::System& s = d->frame;
    const std::string root = dir ? dir : "";
    std::ifstream cf(std::filesystem::path(root) / "catalogue.json");
    if (!cf) throw std::runtime_error("no force-field catalogue in " + root);
    const caps::Json cat = caps::Json::parse(std::string((std::istreambuf_iterator<char>(cf)), std::istreambuf_iterator<char>()));
    const caps::Perception per_all = caps::perceive(s);
    struct Entry { std::string id, name, path; };
    std::vector<Entry> entries{{"uff", "UFF (Rappé et al. 1992)", "uff"}};
    for (const auto& e : cat["forcefields"].items())
      if (!e.text("file").empty()) entries.push_back({e.text("id"), e.text("name"), (std::filesystem::path(root) / e.text("file")).string()});
    caps::Json list = caps::Json::array();
    int k = 0;
    for (const auto& en : entries) {
      if (progress && progress(en.name.c_str(), double(k++) / double(entries.size()), user) != 0) break;
      caps::Json x = caps::Json::object();
      x["id"] = en.id, x["name"] = en.name;
      try {
        caps::ForceField ff;
        caps::ParamReport rep;
        std::string charges = "types";
        if (en.path == "uff") {
          ff = caps::assign_uff(s);
          charges = "none";
          x["untyped"] = 0.0;
        } else {
          std::shared_ptr<caps::FFDef> def;
          {
            std::lock_guard<std::mutex> g(cache_lock);
            auto it = cache.find(en.path);
            if (it == cache.end()) it = cache.emplace(en.path, std::make_shared<caps::FFDef>(caps::load_forcefield(en.path))).first;
            def = it->second;
          }
          if (def->typing.empty()) { x["status"] = "no typing rules"; list.push_back(std::move(x)); continue; }
          // a united-atom force field sees the structure with its hydrogens on carbon folded in (charges summed)
          caps::System su;
          caps::Perception pu;
          std::string ua_charges = "auto";
          if (caps::needs_prepare(*def)) {
            su = s;
            const std::string note = caps::prepare_for_forcefield(su, *def, ua_charges);
            if (!note.empty()) x["united_atom"] = note;
            pu = caps::perceive(su);
          }
          const bool prepared = caps::needs_prepare(*def);
          const caps::System& s = prepared ? su : d->frame;   // NOLINT: shadows the all-atom structure on purpose
          const caps::Perception& per = prepared ? pu : per_all;
          caps::TypingResult tr = caps::assign_types(s, *def);
          {   // a rule may give a type this file lacks (rules shared with a larger version): untyped, as in field_run
            std::set<std::string> known;
            for (const auto& ty : def->types) known.insert(ty.name);
            for (const auto& ty : def->types) for (const auto& a : ty.aliases) known.insert(a);
            for (auto& ty : tr.types)
              if (!ty.empty() && !known.count(ty)) { ty.clear(); ++tr.untyped; }
          }
          x["untyped"] = double(tr.untyped);
          if (tr.untyped) {
            std::map<std::string, std::pair<int, std::vector<double>>> groups;
            for (uint32_t i = 0; i < tr.types.size(); ++i)
              if (tr.types[i].empty()) {
                auto& g = groups[atom_env(s, per, i)];
                ++g.first;
                if (g.second.size() < 200) g.second.push_back(double(i));
              }
            caps::Json gl = caps::Json::array();
            for (const auto& [env, g] : groups) {
              caps::Json ge = caps::Json::object();
              caps::Json at = caps::Json::array();
              for (double a : g.second) at.push_back(a);
              ge["environment"] = env, ge["count"] = double(g.first), ge["atoms"] = std::move(at);
              gl.push_back(std::move(ge));
            }
            x["untyped_groups"] = std::move(gl);
            x["status"] = "untyped atoms";
            list.push_back(std::move(x));
            continue;
          }
          if (ua_charges == "keep") {   // united-atom: the all-atom charges, summed into each site
            ff = caps::parameterize(s, *def, tr.types, "keep", &rep, true);
            charges = "gasteiger";
          } else try {
            ff = caps::parameterize(s, *def, tr.types, "types", &rep, true);
          } catch (const caps::FFError& e) {
            if (std::string(e.what()).find("has no charge for type") == std::string::npos) throw;
            rep = caps::ParamReport{};
            try {
              ff = caps::parameterize(s, *def, tr.types, "gasteiger", &rep, true);
              charges = "gasteiger";
            } catch (const std::exception&) {
              rep = caps::ParamReport{};
              ff = caps::parameterize(s, *def, tr.types, "qeq", &rep, true);
              charges = "qeq";
            }
          }
        }
        caps::Json miss = caps::Json::array();
        for (size_t m = 0; m < rep.missing.size() && m < 8; ++m) miss.push_back(rep.missing[m]);
        x["missing"] = std::move(miss);
        x["missing_count"] = double(rep.missing.size());
        x["estimated"] = double(rep.estimated_terms);
        x["charges"] = charges;
        double net = 0;
        for (double q : ff.charge) net += q;
        x["net_charge"] = net;
        // the physics check: a force field's own charges must leave every molecule neutral (the structure's formal
        // charge aside); with Gasteiger or no charges this says nothing about the force field
        int formal = 0;
        for (int c : per_all.charge) formal += c;
        const bool balanced = charges != "types" || std::fabs(net - formal) < 1e-3;
        x["balanced"] = balanced;
        x["complete"] = rep.missing.empty() && balanced;
        x["status"] = !rep.missing.empty() ? "missing parameters" : !balanced ? "charges do not balance" : "complete";
      } catch (const std::exception& e) {
        x["status"] = std::string("error: ") + e.what();
      }
      list.push_back(std::move(x));
    }
    r["ok"] = true;
    r["forcefields"] = std::move(list);
  } catch (const std::exception& e) {
    r["ok"] = false;
    r["error"] = std::string(e.what());
  }
  return report_out(r.dump(), out, cap);
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

int32_t caps_field_type_by_example(caps_doc* d, caps_doc* example, const char* types_json, char* report, int32_t cap) {
  return guard([&] {
    if (!d->field) throw caps::FFError("assign a force field first");
    if (!example) throw caps::FFError("no example");
    const caps::Json tj = caps::Json::parse(types_json ? types_json : "[]");
    std::vector<std::string> types;
    for (const auto& x : tj.items()) types.push_back(x.is_string() ? x.str() : "");
    const caps::System& ex = example->frame;
    if (types.size() != ex.atoms.size()) throw caps::FFError("one type per example atom (" + std::to_string(ex.atoms.size()) + "), " + std::to_string(types.size()) + " given");
    // types the force field does not have are not learned (and reported)
    caps::Json unknown = caps::Json::array();
    std::set<std::string> said;
    for (auto& t : types)
      if (!t.empty() && !d->field->base.type(t) && !d->field->extra.type(t)) {
        if (said.insert(t).second) unknown.push_back(t);
        t.clear();
      }
    const caps::ExampleTypes learned = caps::learn_types(ex, types);
    const caps::ExampleMatch m = caps::apply_types(d->frame, learned);
    size_t set = 0;
    for (size_t i = 0; i < m.types.size(); ++i)
      if (!m.types[i].empty()) d->field->overrides[int32_t(i)] = m.types[i], ++set;
    field_run(d);
    caps::Json r = caps::Json::object();
    r["radius"] = double(learned.radius);
    r["environments"] = double(learned.environments);
    r["exact"] = double(m.exact);
    r["shorter"] = double(m.shorter);
    r["unmatched"] = double(m.unmatched);
    r["set"] = double(set);
    caps::Json c = caps::Json::array();
    for (const auto& x : learned.conflicts) c.push_back(x);
    r["conflicts"] = c;
    r["unknown_types"] = unknown;
    report_out(r.dump(), report, cap);
    return d->field->complete ? 0 : 1;
  });
}

int32_t caps_equivalent_atoms(caps_doc* d, int32_t atom, int32_t radius, char* json, int32_t cap) {
  int32_t n = -1;
  guard([&] {
    caps::Json a = caps::Json::array();
    for (uint32_t i : caps::equivalent_atoms(d->frame, uint32_t(std::max(0, atom)), std::clamp(radius, 0, 8))) a.push_back(double(i));
    n = report_out(a.dump(), json, cap);
    return 0;
  });
  return n;
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
    d->field->fill = caps::FFDef{};
    d->field->imported.clear();
    field_run(d);
    return d->field->complete ? 0 : 1;
  });
}

int32_t caps_field_import(caps_doc* d, const char* path) { return caps_field_import_ex(d, path, nullptr); }


int32_t caps_field_import_ex(caps_doc* d, const char* path, const char* options) {
  return guard([&] {
    const caps::Json opt = caps::Json::parse(options && *options ? options : "{}");
    const bool fill = opt.text("mode", "override") == "fill";
    if (!d->field) throw caps::FFError("assign a force field first");
    const std::string p = path ? path : "";
    auto ends = [&](const char* e) { const std::string x = e; return p.size() > x.size() && p.compare(p.size() - x.size(), x.size(), x) == 0; };
    std::string lower_p = p;
    for (auto& ch : lower_p) ch = char(std::tolower(static_cast<unsigned char>(ch)));
    const bool frcmod = lower_p.find("frcmod") != std::string::npos || ends(".dat");
    const bool gmx = ends(".itp") || ends(".top");
    caps::FFDef imp = ends(".lt") ? caps::import_moltemplate(p) : frcmod ? caps::import_frcmod(p) : gmx ? caps::import_gromacs_params(p) : caps::load_forcefield(p);
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
    if (fill) {
      // borrowed only where the force field defines nothing (caps::gap_fill_rules)
      caps::GapFill g = caps::gap_fill_rules(d->field->base, imp, file);
      if (!g.kept) throw caps::FFError(file + ": none of its rules can apply to " + d->field->base.name + " (no shared atom classes)");
      auto& F = d->field->fill;
      for (const auto& p2 : {std::make_pair(&F.bonds, &g.rules.bonds), std::make_pair(&F.angles, &g.rules.angles), std::make_pair(&F.dihedrals, &g.rules.dihedrals)})
        p2.first->insert(p2.first->end(), p2.second->begin(), p2.second->end());
      d->field->imported.push_back(p + " (gaps only: " + std::to_string(g.kept) + " rules usable" +
                                   (g.dropped ? ", " + std::to_string(g.dropped) + " on classes " + d->field->base.name + " lacks" : "") + ")");
      field_run(d);
      return d->field->complete ? 0 : 1;
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

int32_t caps_set_live(caps_doc* d, caps_live_fn fn, void* user) {
  return guard([&] {
    d->live_fn = fn;
    d->live_user = user;
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

extern "C++" {
// What the view draws: the shown frame (with periodic images around the cell when asked, faded) and its options.
const caps::System& view_system(caps_doc* d, caps::RenderOptions& ro, caps::System& imaged) {
  const caps::System& base = shown(d);
  const bool images = (d->images[0] * d->images[1] * d->images[2] > 1) && base.cell.valid() && !d->pstate;
  if (!images) return base;
  // copies of the frame around the cell, faded; picks map back to the original atoms
  imaged = base;
  const size_t n = base.atoms.size();
  ro.faded.assign(n, 0);
  for (int a = 0; a < d->images[0]; ++a)
    for (int b = 0; b < d->images[1]; ++b)
      for (int c = 0; c < d->images[2]; ++c) {
        const int ia = a - (d->images[0] - 1) / 2, ib = b - (d->images[1] - 1) / 2, ic = c - (d->images[2] - 1) / 2;
        if (!ia && !ib && !ic) continue;
        const caps::Vec3 t = base.cell.a * double(ia) + base.cell.b * double(ib) + base.cell.c * double(ic);
        const uint32_t off = uint32_t(imaged.atoms.size());
        for (const auto& at : base.atoms) { imaged.atoms.push_back(at); imaged.atoms.back().pos = at.pos + t; }
        for (const auto& bd : base.bonds) imaged.bonds.push_back({bd.i + off, bd.j + off, bd.order});
        ro.faded.insert(ro.faded.end(), n, 1);
      }
  ro.fade = d->image_fade;
  return imaged;
}
}  // extern "C++"

int32_t caps_render(caps_doc* d, const caps_camera* cam, const caps_render_opts* opt, uint8_t* rgba) {
  return guard([&] {
    auto ro = opts_of(d, opt);
    caps::System imaged;
    const caps::System& sys = view_system(d, ro, imaged);
    auto img = d->renderer.render(sys, cam_of(cam), ro);
    if (d->vision) caps::simulate_vision(img, caps::Vision(d->vision), d->vision_severity);
    std::memcpy(rgba, img.rgba.data(), img.rgba.size());
    return 0;
  });
}

int32_t caps_render_scene(caps_doc* d, const caps_render_opts* opt, caps_scene* out) {
  return guard([&] {
    auto ro = opts_of(d, opt);
    caps::System imaged;
    const caps::System& sys = view_system(d, ro, imaged);
    d->scene = d->renderer.scene(sys, ro);
    const auto& sc = d->scene;
    std::memset(out, 0, sizeof *out);
    out->n_spheres = int32_t(sc.sphere_id.size());
    out->spheres = sc.spheres.data(), out->sphere_rgb = sc.sphere_rgb.data(), out->sphere_id = sc.sphere_id.data(), out->sphere_ring = sc.sphere_ring.data();
    out->n_capsules = int32_t(sc.capsule_rgb.size());
    out->capsules = sc.capsules.data(), out->capsule_rgb = sc.capsule_rgb.data();
    out->n_lines = int32_t(sc.line_rgb.size());
    out->lines = sc.lines.data(), out->line_rgb = sc.line_rgb.data(), out->line_width = sc.line_width.data();
    out->cpu_only = sc.has_meshes || d->vision ? 1 : 0;   // surfaces, polyhedra and the colour-vision preview are drawn on the CPU
    out->background = sc.background;
    out->transparent = sc.transparent, out->dark = sc.dark, out->depth_cue = sc.depth_cue, out->outlines = sc.outlines;
    out->fit_cx = sc.fit_centre[0], out->fit_cy = sc.fit_centre[1], out->fit_cz = sc.fit_centre[2];
    out->fit_pad = sc.fit_pad, out->fov_deg = sc.fov_deg;
    out->n_fit_corners = int32_t(sc.fit_corners.size() / 3), out->fit_corners = sc.fit_corners.data();
    out->n_fit_points = int32_t(sc.fit_points.size() / 3), out->fit_points = sc.fit_points.data();
    return 0;
  });
}

int32_t caps_view_fit(caps_doc* d, const caps_camera* cam, const caps_render_opts* opt, caps_view_fit_t* out) {
  return guard([&] {
    auto ro = opts_of(d, opt);
    caps::System imaged;
    const caps::System& sys = view_system(d, ro, imaged);
    const caps::ViewFit f = caps::view_fit(sys, cam_of(cam), ro);
    out->cos_yaw = f.cos_yaw, out->sin_yaw = f.sin_yaw, out->cos_pitch = f.cos_pitch, out->sin_pitch = f.sin_pitch;
    out->cx = f.centre[0], out->cy = f.centre[1], out->cz = f.centre[2];
    out->scale = f.scale, out->w = f.w, out->h = f.h, out->pan_x = f.pan_x, out->pan_y = f.pan_y;
    out->perspective = f.perspective ? 1 : 0;
    out->dist = f.dist, out->zmin = f.zmin, out->zmax = f.zmax;
    return 0;
  });
}

int32_t caps_pick(caps_doc* d, int32_t x, int32_t y) {
  int k = d->renderer.pick(x, y);
  if (k >= 0 && !d->pstate && !d->frame.atoms.empty()) k = int(size_t(k) % d->frame.atoms.size());   // an image atom picks its original
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
          if (!cur.empty()) (cur == "cij_strain" || cur == "cij_run" || cur == "viscosity" || cur == "tensile" || cur == "tg" || cur == "pull_shear" || cur == "pull_normal" ? protocols : ids).push_back(cur);
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
      o.deuterate = std::clamp(p->deuterate, 0, 4);
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
      } else if (id == "viscosity") {
        caps::System s = frame_copy();
        caps::ViscosityOptions vo;
        vo.field = ff;
        vo.energy = o.energy;
        if (mo.temperature > 0) vo.temperature = mo.temperature;
        if (mo.run_ps > 0) vo.ps = mo.run_ps;
        if (mo.dt > 0) vo.dt = mo.dt;
        if (mo.seed) vo.seed = mo.seed;
        vo.progress = [&](const std::string& w, double f) { return !cancelled(w, f); };
        for (auto& q : caps::viscosity_properties(caps::viscosity_green_kubo(s, vo))) res.push_back(std::move(q));
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

extern "C" caps_doc* caps_build_beads(const char* text, const char* ff_path, uint64_t seed, char* report, int32_t cap) {
  try {
    caps::FFDef def;
    if (ff_path && *ff_path) def = caps::load_forcefield(ff_path);
    const std::string t = text ? text : "";
    const bool tpl = caps::has_bead_template(def, t);
    const caps::System s = caps::build_bead_molecule(t, def, seed ? seed : 1);
    auto* d = new caps_doc;
    d->traj.topology = s;
    d->traj.positions.push_back({});
    for (const auto& a : s.atoms) d->traj.positions.back().push_back(a.pos);
    d->traj.cells.push_back(s.cell);
    d->traj.timesteps.push_back(0);
    refresh(d);
    prov_step(d, "cg.build", tpl ? "coarse-grained molecule from the " + def.name + " template " + t : "coarse-grained molecule from bead SMILES",
              {{"beads", tpl ? caps::bead_template_list(def).at(t) : t}, {"forcefield", def.name.empty() ? "none (4.7 Å bonds)" : def.name}}, seeded(seed ? seed : 1), {});
    double q = 0;
    for (const auto& a : s.atoms) q += a.charge;
    caps::Json j = caps::Json::object();
    j["beads"] = double(s.atoms.size());
    j["bonds"] = double(s.bonds.size());
    j["charge"] = q;
    j["template"] = tpl ? t : std::string();
    j["forcefield"] = def.name;
    report_out(j.dump(), report, cap);
    return d;
  } catch (const std::exception& e) {
    g_error = e.what();
    return nullptr;
  }
}

extern "C" int32_t caps_bead_templates(const char* ff_path, char* json, int32_t cap) {
  try {
    const caps::FFDef def = caps::load_forcefield(ff_path ? ff_path : "");
    caps::Json j = caps::Json::object();
    for (const auto& [k, v] : caps::bead_template_list(def)) j[k] = v;
    return report_out(j.dump(), json, cap);
  } catch (const std::exception& e) {
    g_error = e.what();
    return -1;
  }
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
    {
      const bool uff = caps::is_uff(b.forcefield);
      prov_step(d, "chem.build", "3D structure from SMILES · " + r.method,
                {{"smiles", r.graph.smiles}, {"conformers", std::to_string(b.conformers)}, {"clean-up", b.forcefield.empty() ? "embedding only" : uff ? "UFF" : std::filesystem::path(b.forcefield).filename().string()}},
                seeded(b.seed), uff ? std::vector<std::string>{"rappe1992"} : std::vector<std::string>{});
    }
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
  c.r1 = j.num("r1", 1), c.r2 = j.num("r2", 1);
  c.p_mr = j.num("p_mr", -1), c.p_rm = j.num("p_rm", -1);
  if (j.has("chain_dp")) for (const auto& x : j["chain_dp"].items()) c.chain_dp.push_back(int(x.number()));
  c.forcefield = j.text("forcefield");
  const std::string tac = j.text("tacticity", "atactic");
  c.tacticity = caps::tacticity_from_string(tac);
  c.architecture = caps::architecture_from_string(j.text("architecture", "linear"));
  c.arms = int(j.num("arms", 4));
  c.arm_dp = int(j.num("arm_dp", 5));
  c.spacing = int(j.num("spacing", 4));
  c.branch_probability = j.num("branch_probability", 0.1);
  c.generations = int(j.num("generations", 2));
  c.keep_configuration = j.num("keep_configuration", 0) != 0 || (j.has("keep_configuration") && j["keep_configuration"].kind() == caps::Json::Bool && j["keep_configuration"].boolean());
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
    j["ring_backbone"] = u.ring_backbone;
    if (u.ring_stereo_open)
      j["note"] = std::string("head and tail sit on one ring: give their configuration with @/@@ (e.g. 2,3-exo,exo for norbornenes); ") +
                  "without it the embedding picks one, and some (endo,endo) cannot form a chain";
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
    // a branched molecule: its arms (expected count for random branches), atoms and mass; each arm replaces a
    // hydrogen of the molecule and loses its own head cap
    if (c.architecture != caps::Architecture::Linear) {
      double arms = 0;
      caps::ChainSpec a = c;
      a.architecture = caps::Architecture::Linear;
      double core_arms = 0;   // a dendrimer: the core's arms of dp units besides its arm_dp branches
      if (c.architecture == caps::Architecture::Star) arms = c.arms - 1;
      else if (c.architecture == caps::Architecture::Dendrimer) {
        a.dp = std::max(1, c.arm_dp);
        core_arms = c.arms - 1;
        arms = c.arms * (std::pow(2.0, std::clamp(c.generations, 1, 6) + 1) - 2);
      } else {
        a.dp = std::max(1, c.arm_dp);
        arms = c.architecture == caps::Architecture::Comb ? double(c.dp / std::max(1, c.spacing)) : c.branch_probability * std::max(0, c.dp - 2);
      }
      const caps::MolInfo am = caps::molecule_info(caps::chain_graph(a, caps::chain_sequence(a, seed + 1)));
      caps::Json mo = caps::Json::object();
      mo["architecture"] = std::string(caps::to_string(c.architecture));
      mo["arms"] = arms;
      mo["atoms"] = double(m.atoms) + arms * (am.atoms - 2) + core_arms * (m.atoms - 2);
      mo["mass"] = m.mass + arms * (am.mass - 2 * 1.008) + core_arms * (m.mass - 2 * 1.008);
      if (core_arms > 0) mo["arms"] = arms + core_arms, mo["branches"] = arms;
      j["molecule"] = mo;
    }
  } catch (const std::exception& e) {
    j["ok"] = false;
    j["error"] = std::string(e.what());
  }
  return report_out(j.dump(), json, cap);
}

namespace {
caps_doc* grow_chains_impl(const char* spec_json, const caps_grow_opts* o, caps_progress_fn progress, caps_grow_live_fn live, void* user, char* report, int32_t cap);
}

extern "C" caps_doc* caps_grow_chains(const char* spec_json, const caps_grow_opts* o, caps_progress_fn progress, void* user, char* report, int32_t cap) {
  return grow_chains_impl(spec_json, o, progress, nullptr, user, report, cap);
}

extern "C" caps_doc* caps_grow_chains_live(const char* spec_json, const caps_grow_opts* o, caps_progress_fn progress, caps_grow_live_fn live, void* user, char* report,
                                           int32_t cap) {
  return grow_chains_impl(spec_json, o, progress, live, user, report, cap);
}

namespace {
caps_doc* grow_chains_impl(const char* spec_json, const caps_grow_opts* o, caps_progress_fn progress, caps_grow_live_fn live, void* user, char* report, int32_t cap) {
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
    try {   // optional growth settings carried in the spec: trial directions per step
      const caps::Json sj = caps::Json::parse(spec_json ? spec_json : "{}");
      if (sj.has("trials")) g.trials = std::clamp(int(sj["trials"].number()), 4, 5000);
      // method: trials (best by contact margin) | rosenbluth (soft spheres) | rosenbluth_lj (UFF Lennard-Jones); temperature K
      const std::string meth = sj.text("method", "trials");
      g.method = meth == "rosenbluth" ? 1 : meth == "rosenbluth_lj" ? 2 : 0;
      g.method_temperature = sj.num("temperature", 450);
      // region: {shape: slab, thickness, vacuum} | {shape: cylinder | around_cylinder, radius, length}
      if (sj.has("region") && sj["region"].is_object()) {
        const caps::Json& R = sj["region"];
        const std::string shape = R.text("shape", "cubic");
        if (shape == "slab") g.slab_thickness = R.num("thickness", 30), g.slab_vacuum = R.num("vacuum", 30);
        else if (shape == "cylinder" || shape == "around_cylinder")
          g.cylinder_radius = R.num("radius", 10), g.cylinder_length = R.num("length", 0), g.cylinder_outside = shape == "around_cylinder";
      }
    } catch (...) {}
    if (live)   // the chains so far as a new document each time (the callee closes it), with where the growth stands
      g.snapshot = [&](const caps::System& part, const caps::GrowOptions::Live& L) {
        auto* sd = new caps_doc;
        sd->traj.topology = part;
        std::vector<caps::Vec3> pp;
        for (const auto& a : part.atoms) pp.push_back(a.pos);
        sd->traj.positions.push_back(std::move(pp));
        sd->traj.cells.push_back(part.cell);
        sd->traj.timesteps.push_back(0);
        refresh(sd);
        caps::Json j = caps::Json::object();
        j["chains_done"] = double(L.chains_done), j["chains"] = double(L.chains), j["units"] = double(L.units), j["units_total"] = double(L.units_total);
        j["restarts"] = double(L.restarts), j["worst_margin"] = L.worst_margin, j["density"] = L.density, j["atoms"] = double(part.atoms.size());
        live(sd, j.dump(0).c_str(), user);
      };
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
    {
      caps::KeyValues pr = json_params(spec_json);
      pr.push_back({"chains", std::to_string(o->chains)});
      pr.push_back({o->box > 0 ? "box" : "density", o->box > 0 ? g6(o->box) + " Å" : g6(o->density) + " g/cm³"});
      pr.push_back({"contact scale", g6(g.contact_scale) + (g.auto_scale ? " (lowered when crowded)" : "")});
      prov_step(d, "grow.trials", std::to_string(o->chains) + " chains grown in a periodic cell, best-of-k trial placement by contact margin", std::move(pr),
                seeded(o->seed), {"matsumoto1998"});
    }
    std::string t;
    for (const auto& n : rep.notes) t += n + "\n";
    report_out(t, report, cap);
    return d;
  } catch (const std::exception& e) {
    g_error = e.what();
    return nullptr;
  }
}
}  // namespace

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
    caps_doc* d = doc_of(s);
    prov_input(d, cif_path ? cif_path : "");
    prov_step(d, "surface.build", "slab cut from a crystal", json_params(options_json));
    return d;
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
    {
      caps::KeyValues pr = json_params(options_json);
      for (auto& kv : json_params(spec_json)) pr.push_back({"chain " + kv.first, kv.second});
      prov_step(d, "interface.build", "polymer grown against a surface", std::move(pr), seeded(o ? o->seed : 0), {"matsumoto1998"});
    }
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
    t.cc = j.num("cc", 1.42);
    t.walls = int(j.num("walls", 1));
    f = caps::nanotube(t, &r);
    keep = {false, false, t.periodic};
  } else if (kind == "sheet") {
    caps::SheetOptions sh;
    sh.lx = j.num("lx", 20), sh.ly = j.num("ly", 20);
    sh.layers = int(j.num("layers", 1));
    sh.periodic = j.num("periodic", 1) != 0;
    sh.cc = j.num("cc", 1.42);
    f = caps::graphene_sheet(sh, &r);
    keep = {sh.periodic, sh.periodic, false};
  } else if (kind == "particle") {
    caps::ParticleOptions po;
    po.shape = caps::particle_shape_from_string(j.text("shape", "sphere"));
    po.radius = j.num("radius", 12);
    po.on_atom = j.num("on_atom", 1) != 0;
    po.passivate = j.num("passivate", 0) != 0;
    po.thiolate = j.text("thiolate", "");
    po.thiolate_fraction = j.num("thiolate_fraction", 1.0);
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
    caps_doc* d = doc_of(f);
    prov_step(d, "nano.build", "nanostructure", json_params(options_json));
    return d;
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
    {
      caps::KeyValues pr = json_params(options_json);
      for (auto& kv : json_params(spec_json)) pr.push_back({"chain " + kv.first, kv.second});
      prov_step(d, "nano.embed", "nanostructure in a grown polymer matrix", std::move(pr), seeded(o ? o->seed : 0), {"matsumoto1998"});
    }
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
    caps_doc* d = doc_of(s);
    prov_step(d, "grow.blend", "polymer blend grown in a periodic cell", json_params(options_json), seeded(o ? o->seed : 0), {"matsumoto1998"});
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
    prov_step(d, "pack.insert", std::to_string(std::max(1, count)) + " × " + std::string(smiles ? smiles : "") + " inserted",
              {{"smiles", smiles ? smiles : ""}, {"count", std::to_string(std::max(1, count))}, {"tolerance", g6(po.tolerance) + " Å"}}, seeded(po.seed),
              {"martinez2009", "rappe1992"});
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
  else if (fmt == "car") {
    caps::write_car(s, path);
    notes.push_back("the .mdf with the bonds beside it; atom types as the force field (or the file) names them");
  }
  else if (fmt == "gro") caps::write_gro(s, path);
  else if (fmt == "sdf" || fmt == "mol") {
    caps::write_sdf(s, path);
    notes.push_back(s.atoms.size() > 999 || s.bonds.size() > 999 ? "V3000 connection table (more than 999 atoms or bonds)" : "V2000 connection table");
    notes.push_back("formal charges from the structure's chemistry; partial charges are not part of the format");
  } else if (fmt == "cif") {
    caps::write_cif(s, path);
    notes.push_back("space group P 1: every atom at its fractional coordinates; no bonds");
  } else if (fmt == "dcd") {
    caps::write_dcd(d->traj, path);
    notes.push_back("every frame with its cell, single precision (as LAMMPS writes DCD); open it with a topology (the .data)");
  } else if (fmt == "lammps-dump") {
    caps::write_lammps_dump(d->traj, path);
    notes.push_back("every frame, unwrapped coordinates; the pipeline is not applied to a dump");
  } else if (fmt == "lammps-data") {
    notes.push_back("wrapped coordinates with image flags");
    const bool same = s.atoms.size() == d->frame.atoms.size() && !use_pipeline;
    if (coeffs && d->field && d->field->complete && same) {
      const std::string why = caps::write_lammps_data_or_structure(s, *d->field->ff, elec(), path);
      notes.push_back(why.empty() ? "coefficients from Field: " + d->field->ff->name : why);
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

extern "C" int32_t caps_pipeline_export_grid(caps_doc* d, const char* path) {
  return guard([&] {
    if (!d->pstate || !d->pstate->grid) throw std::invalid_argument("no grid: add a Density field step to the pipeline");
    caps::write_grid(*d->pstate->grid, d->pstate->system, path ? path : "");
    return 0;
  });
}

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

extern "C" void caps_set_ph(caps_doc* d, double ph) {
  if (d) d->ph = ph;
}

extern "C" int32_t caps_set_restraints(caps_doc* d, const char* json) {
  if (!d) return -1;
  int32_t n = 0;
  const int32_t rc = guard([&] {
    std::vector<caps::RelaxOptions::Restraint> v;
    std::vector<caps::RelaxOptions::DihedralRestraint> dv;
    if (json && *json) {
      const caps::Json j = caps::Json::parse(json);
      for (const auto& e : j.items()) {
        if (e.has("l")) {   // four atoms: a dihedral restraint {i, j, k, l, phi0 (°), kphi (kcal/mol/rad²)}
          caps::RelaxOptions::DihedralRestraint q;
          q.i = uint32_t(e.num("i", 0)), q.j = uint32_t(e.num("j", 0)), q.k = uint32_t(e.num("k", 0)), q.l = uint32_t(e.num("l", 0));
          q.phi0 = e.num("phi0", 180), q.kphi = e.num("kphi", 50);
          if (q.kphi < 0) throw std::invalid_argument("a dihedral restraint needs kphi ≥ 0");
          dv.push_back(q);
          continue;
        }
        caps::RelaxOptions::Restraint r;
        r.i = uint32_t(e.num("i", 0)), r.j = uint32_t(e.num("j", 0)), r.r0 = e.num("r0", 0), r.k = e.num("k", 10);
        if (r.r0 < 0 || r.k < 0) throw std::invalid_argument("a restraint needs r0 ≥ 0 and k ≥ 0");
        v.push_back(r);
      }
    }
    d->restraints = std::move(v);
    d->dihedral_restraints = std::move(dv);
    n = int32_t(d->restraints.size() + d->dihedral_restraints.size());
    return 0;
  });
  return rc < 0 ? -1 : n;
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
    prov_step(d, "crystal.build", "crystal from a space group and its asymmetric unit", json_params(spec_json), "", {"hall1981"});
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
    prov_step(d, "bio.peptide", "peptide from its sequence, backbone placed by NeRF", json_params(options_json), "", {"engh1991", "parsons2005"});
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
    caps_doc* d = doc_of(out);
    if (solute) { d->prov = solute->prov; }
    {
      std::vector<std::string> c = {"martinez2009"};
      const std::string wm = o.water_model;
      if (o.solvent == "water") c.push_back(wm == "TIP4P/2005" ? "abascal2005" : wm == "TIP3P" ? "jorgensen1983" : wm == "SPC/E" ? "berendsen1987" : "");
      c.erase(std::remove(c.begin(), c.end(), std::string()), c.end());
      prov_step(d, "solvate.pack", std::to_string(rep.pack.molecules) + " molecules packed · closest contact " + g6(rep.pack.dmin) + " Å", json_params(options_json), "", std::move(c));
    }
    return d;
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
  } else if (k == "ez" || k == "stereo") {   // E/Z on double-bond atoms; "stereo": R/S and E/Z together
    const auto ez = caps::ez_labels(f);
    const auto rs = k == "stereo" ? caps::stereo_labels(f) : std::vector<std::string>(f.atoms.size());
    for (size_t i = 0; i < f.atoms.size(); ++i) arr.push_back(rs[i].empty() ? ez[i] : rs[i]);
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
  // snapshots taken in the undone steps now belong to the branch those steps become
  if (!d->redo.empty())
    for (auto& n : d->snapshots)
      if (size_t(n.step) > d->undo.size() + 1) n.on_branch = true;
  caps_doc::Snapshot sn;
  sn.topology = d->traj.topology;
  sn.positions = d->traj.positions.at(d->current);
  sn.what = what;
  d->undo.push_back(std::move(sn));
  if (d->undo.size() > 100) d->undo.erase(d->undo.begin());
  // editing after an undo keeps the undone steps as a branch: nothing is lost until it is deleted
  if (!d->redo.empty()) {
    d->branches.push_back(std::move(d->redo));
    if (d->branches.size() > 10) d->branches.erase(d->branches.begin());
  }
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

namespace {
// No hydrogen and not one multiple bond: the orders were never assigned (a PDB opened directly), so read them from the
// geometry before counting what is missing.
// A structure whose bonds carry no orders (a LAMMPS data file, an XYZ): all single, so aromatic rings and double bonds
// would read as missing hydrogens.
bool needs_geometry_orders(const caps::System& s) {
  if (s.atoms.size() < 3 || s.bonds.empty()) return false;
  for (const auto& b : s.bonds) if (b.order >= 2) return false;
  return true;
}

// Bond orders from the geometry of the heavy atoms (orders_from_geometry works on them), mapped back; bonds to
// hydrogen stay single. Polystyrene from a LAMMPS data file: aromatic rings, sp³ backbone — no hydrogen missing.
void geometry_orders(caps::System& s) {
  caps::System heavy;
  heavy.cell = s.cell;
  std::vector<int64_t> to(s.atoms.size(), -1);
  for (size_t i = 0; i < s.atoms.size(); ++i)
    if (s.atoms[i].element != 1) { to[i] = int64_t(heavy.atoms.size()); heavy.atoms.push_back(s.atoms[i]); }
  std::vector<size_t> which;
  for (size_t k = 0; k < s.bonds.size(); ++k) {
    const auto& b = s.bonds[k];
    if (to[b.i] < 0 || to[b.j] < 0) continue;
    heavy.bonds.push_back({uint32_t(to[b.i]), uint32_t(to[b.j]), 1});
    which.push_back(k);
  }
  if (heavy.atoms.size() == s.atoms.size()) { caps::orders_from_geometry(s); return; }
  caps::orders_from_geometry(heavy);
  for (size_t m = 0; m < which.size(); ++m) s.bonds[which[m]].order = heavy.bonds[m].order;
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
      if (needs_geometry_orders(s)) geometry_orders(s);   // orders never assigned: from the heavy atoms' geometry
      const double ph = j.has("ph") ? j["ph"].number() : d->ph;
      const int k = ph >= 0 ? caps::add_hydrogens_at_ph(s, ph, m) : caps::add_hydrogens(s, m);
      if (k == 0) throw std::invalid_argument("no atom lacks hydrogens");
      for (size_t i = before; i < s.atoms.size(); ++i) added.push_back(double(i));
      what = "Add " + std::to_string(k) + " hydrogens" + (ph >= 0 ? " at pH " + std::to_string(ph).substr(0, 4) : "");
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
    } else if (op == "graft") {   // silane coupling agents on the surface silanols of silica
      caps::GraftOptions g;
      const std::string silane = j.text("silane", "TESPT");
      g.name = silane;
      g.smiles = silane.find('*') != std::string::npos ? silane : caps::silane_smiles(silane);
      if (silane.find('*') != std::string::npos) g.name = j.text("name", "silane");
      g.count = int(j.num("count", 0));
      g.fraction = j.num("fraction", 0.25);
      g.min_spacing = j.num("min_spacing", 5.0);
      g.seed = uint64_t(j.num("seed", 1));
      const auto gr = caps::graft_silanes(s, g);
      what = "Graft " + g.name + " on " + std::to_string(gr.grafted) + " of " + std::to_string(gr.silanols) + " silanols";
    } else if (op == "phosphate_ends") {   // P–H (a strand's 3′ cap) → P–OH
      const int k = caps::hydroxylate_phosphorus(s);
      if (k == 0) throw std::invalid_argument("no hydrogen on phosphorus");
      what = "Phosphate ends: " + std::to_string(k) + " P–H to P–OH";
    } else if (op == "translate") {   // {atoms | "selection", by: [dx, dy, dz]} Å: the atoms moved rigidly
      const auto at = atoms_of(d, j);
      if (at.empty()) throw std::invalid_argument("pick or select the atoms to move");
      if (!j.has("by") || !j["by"].is_array() || j["by"].size() != 3) throw std::invalid_argument("translate needs by: [dx, dy, dz]");
      const caps::Vec3 by{j["by"][0].number(), j["by"][1].number(), j["by"][2].number()};
      for (uint32_t a : at) s.atoms[a].pos = s.atoms[a].pos + by;
      char b[96];
      std::snprintf(b, sizeof b, "Move %zu atom(s) by %.2f Å", at.size(), caps::norm(by));
      what = b;
    } else if (op == "set_geometry") {   // {atoms: [i, j (, k (, l))], value}: the length, angle or dihedral made exact
      if (!j.has("atoms") || !j["atoms"].is_array()) throw std::invalid_argument("set_geometry needs atoms: [i, j …]");
      std::vector<uint32_t> at;
      for (const auto& x : j["atoms"].items()) at.push_back(uint32_t(x.number()));
      const double v = j.num("value", 0);
      char b[120];
      if (at.size() == 2) {
        caps::set_bond_length(s, at[0], at[1], v);
        std::snprintf(b, sizeof b, "Bond %u–%u to %.3f Å", at[0] + 1, at[1] + 1, v);
      } else if (at.size() == 3) {
        caps::set_bond_angle(s, at[0], at[1], at[2], v);
        std::snprintf(b, sizeof b, "Angle %u–%u–%u to %.2f°", at[0] + 1, at[1] + 1, at[2] + 1, v);
      } else if (at.size() == 4) {
        caps::set_torsion(s, at[0], at[1], at[2], at[3], v);
        std::snprintf(b, sizeof b, "Dihedral %u–%u–%u–%u to %.2f°", at[0] + 1, at[1] + 1, at[2] + 1, at[3] + 1, v);
      } else {
        throw std::invalid_argument("pick two, three or four atoms: a bond, an angle or a dihedral");
      }
      what = b;
    } else if (op == "rotate" || op == "mirror") {   // {atoms | "selection", axis: [x, y, z], degrees} / {…, normal: [x, y, z]}
      const auto at = atoms_of(d, j);
      if (at.empty()) throw std::invalid_argument("pick or select the atoms to " + op);
      const char* key = op == "rotate" ? "axis" : "normal";
      if (!j.has(key) || !j[key].is_array() || j[key].size() != 3) throw std::invalid_argument(op + " needs " + key + ": [x, y, z]");
      const caps::Vec3 ax{j[key][0].number(), j[key][1].number(), j[key][2].number()};
      char b[120];
      if (op == "rotate") {
        caps::rotate_atoms(s, at, ax, j.num("degrees", 90));
        std::snprintf(b, sizeof b, "Rotate %zu atom(s) by %.1f°", at.size(), j.num("degrees", 90));
      } else {
        caps::mirror_atoms(s, at, ax);
        std::snprintf(b, sizeof b, "Mirror %zu atom(s)", at.size());
      }
      what = b;
    } else if (op == "set_rs") {   // {centre, to: "R" | "S"}
      const uint32_t c = uint32_t(j.num("centre", -1));
      const std::string to = j.text("to", "R");
      if (!caps::set_configuration(s, c, to)) throw std::invalid_argument("atom " + std::to_string(c + 1) + " is not a stereocentre");
      what = "Centre " + std::to_string(c + 1) + " made " + to;
    } else if (op == "fuse_ring") {   // {i, j}: a benzene ring fused onto the bond i–j, cleaned with UFF
      const uint32_t a = uint32_t(j.num("i", -1)), b = uint32_t(j.num("j", -1));
      const auto at = caps::fuse_benzene(s, a, b);
      for (uint32_t x : at) added.push_back(double(x));
      if (j.num("clean", 1) != 0) {
        std::vector<char> m(s.atoms.size(), 0);
        for (uint32_t x : at) m[x] = 1;
        caps::clean_up(s, m, 0.5);
      }
      d->selection.assign(s.atoms.size(), 0);
      what = "Fuse a benzene ring onto " + std::to_string(a + 1) + "–" + std::to_string(b + 1);
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
    {
      // consecutive builder edits are one step listing its operations
      if (d->prov.steps.empty() || d->prov.steps.back().engine != "edit.builder") {
        prov_step(d, "edit.builder", "", {{"operations", "0"}}, "", op == "clean" ? std::vector<std::string>{"rappe1992"} : std::vector<std::string>{});
      }
      auto& st = d->prov.steps.back();
      const int n = std::stoi(st.params[0].second) + 1;
      st.params[0].second = std::to_string(n);
      st.params.push_back({"#" + std::to_string(n), what});
      if (st.params.size() > 41) st.params.erase(st.params.begin() + 1);   // the last 40 operations
      st.summary = std::to_string(n) + " edit" + (n == 1 ? "" : "s") + " in the builder";
      if (op == "clean" && std::find(st.cites.begin(), st.cites.end(), "rappe1992") == st.cites.end()) st.cites.push_back("rappe1992");
      st.time = caps::now_iso();
    }
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
    if (d->current < d->traj.cells.size()) d->traj.cells[d->current] = d->traj.topology.cell;   // a cell edit comes back too
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
  // the steps in order with the atoms after each: done (the last is current), then undone (next redo first)
  auto atoms = [](const caps_doc::Snapshot& sn) { return double(sn.topology.atoms.size()); };
  const double now = double(d->traj.topology.atoms.size());
  caps::Json steps = caps::Json::array();
  for (size_t i = 0; i < d->undo.size(); ++i) {
    caps::Json e = caps::Json::object();
    e["what"] = d->undo[i].what;
    e["atoms"] = i + 1 < d->undo.size() ? atoms(d->undo[i + 1]) : now;
    e["state"] = std::string(i + 1 == d->undo.size() ? "current" : "done");
    steps.push_back(std::move(e));
  }
  for (size_t i = d->redo.size(); i-- > 0;) {
    caps::Json e = caps::Json::object();
    e["what"] = d->redo[i].what;
    e["atoms"] = atoms(d->redo[i]);
    e["state"] = std::string("undone");
    steps.push_back(std::move(e));
  }
  j["steps"] = std::move(steps);
  j["start_atoms"] = d->undo.empty() ? now : atoms(d->undo.front());
  caps::Json br = caps::Json::array();
  for (const auto& b : d->branches) {
    caps::Json e = caps::Json::object(), st = caps::Json::array();
    for (size_t i = b.size(); i-- > 0;) st.push_back(b[i].what);
    e["steps"] = std::move(st);
    e["atoms"] = b.empty() ? 0.0 : atoms(b.front());
    br.push_back(std::move(e));
  }
  j["branches"] = std::move(br);
  caps::Json sn = caps::Json::array();
  for (const auto& x : d->snapshots) {
    caps::Json e = caps::Json::object();
    e["name"] = x.name;
    e["atoms"] = atoms(x.state);
    e["step"] = x.step;
    e["on_branch"] = x.on_branch;
    sn.push_back(std::move(e));
  }
  j["snapshots"] = std::move(sn);
  return report_out(j.dump(0), json, cap);
}

// Named snapshots and history branches (design/boards/History): {"op": "take", "name"} · {"op": "restore", "index"} ·
// {"op": "delete", "index"} · {"op": "save", "index", "path"} (a file to compare with) · {"op": "branch", "index"} (the
// branch becomes the redo steps; the current redo steps become a branch) · {"op": "drop_branch", "index"}.
extern "C" int32_t caps_snapshot(caps_doc* d, const char* json) {
  return guard([&] {
    const caps::Json j = caps::Json::parse(json && *json ? json : "{}");
    const std::string op = j.text("op");
    const long k = long(j.num("index", -1));
    auto current = [&] {
      caps_doc::Snapshot sn;
      sn.topology = d->traj.topology;
      sn.positions = d->traj.positions.at(d->current);
      return sn;
    };
    if (op == "take") {
      caps_doc::Named n;
      n.state = current();
      n.name = j.text("name", "snapshot " + std::to_string(d->snapshots.size() + 1));
      n.step = int(d->undo.size()) + 1;
      d->snapshots.push_back(std::move(n));
      return 0;
    }
    if (op == "branch" || op == "drop_branch") {
      if (k < 0 || size_t(k) >= d->branches.size()) throw std::out_of_range("no such branch");
      auto b = std::move(d->branches[size_t(k)]);
      d->branches.erase(d->branches.begin() + k);
      if (op == "branch") {
        if (!d->redo.empty()) d->branches.push_back(std::move(d->redo));
        d->redo = std::move(b);
      }
      return 0;
    }
    if (k < 0 || size_t(k) >= d->snapshots.size()) throw std::out_of_range("no such snapshot");
    auto& n = d->snapshots[size_t(k)];
    if (op == "delete") { d->snapshots.erase(d->snapshots.begin() + k); return 0; }
    if (op == "restore") {
      push_undo(d, "Restore snapshot · " + n.name);
      d->traj.topology = n.state.topology;
      d->traj.positions[d->current] = n.state.positions;
      d->field.reset();
      d->selection.assign(d->traj.topology.atoms.size(), 0);
      refresh(d);
      return 0;
    }
    if (op == "save") {
      caps::System s = n.state.topology;
      for (size_t i = 0; i < s.atoms.size() && i < n.state.positions.size(); ++i) s.atoms[i].pos = n.state.positions[i];
      const std::string path = j.text("path");
      if (path.empty()) throw std::invalid_argument("snapshot save: no path");
      auto ends = [&](const char* e) { const std::string x = e; return path.size() >= x.size() && path.compare(path.size() - x.size(), x.size(), x) == 0; };
      if (ends(".pdb")) caps::write_pdb(s, path);
      else if (ends(".xyz")) caps::write_xyz(s, path);
      else if (ends(".mol2")) caps::write_mol2(s, path);
      else caps::write_lammps_data(s, path);
      return 0;
    }
    throw std::invalid_argument("snapshot: unknown op '" + op + "'");
  });
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
    } else if (mode == "query") {
      const auto q = caps::select_query(s, j.text("pattern"), d->selection);
      m = q.atoms;
      r["rings"] = double(q.rings);
    } else if (mode == "expression") {
      const caps::PipelineState st = caps::run_pipeline(s, caps::Pipeline{});
      const auto v = caps::evaluate_expression(st, j.text("pattern"));
      for (size_t i = 0; i < n && i < v.size(); ++i) m[i] = v[i] != 0;
    } else {
      throw std::invalid_argument("unknown selection mode '" + mode + "'");
    }
    const std::string op = j.text("op", "replace");
    auto& sel = d->selection;
    if (op == "preview") {   // count only; the selection stays as it is
      r["ok"] = true;
      r["count"] = double(std::count(m.begin(), m.end(), 1));
      r["matched"] = r["count"];
      return report_out(r.dump(0), out, cap);
    }
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
    const caps::ImportOptions io = import_options(options_json);
    d->traj = caps::import_file(path, topology_path ? topology_path : "", io);
    refresh(d);
    caps::KeyValues ip = {{"bonds", io.bonds == 0 ? "perceived · covalent radii + " + g6(io.tolerance) + " Å" : io.bonds == 1 ? "from the file" : "none"},
                          {"bond orders", io.bond_orders ? "from valences" : "not assigned"}, {"molecules", io.split ? "by connectivity" : "as read"},
                          {"unwrap", io.unwrap ? "yes" : "no"}, {"cell", io.use_cell ? "from the file" : "dropped"}};
    prov_opened(d, path, topology_path ? topology_path : "", "io.import", std::move(ip), io.bonds == 0 ? std::vector<std::string>{"cordero2008"} : std::vector<std::string>{});
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

// ---------------------------------------------------------------- provenance (design/boards/Provenance)
extern "C" int32_t caps_provenance(caps_doc* d, char* json, int32_t cap) {
  return report_out(caps::manifest_json(d->prov).dump(0), json, cap);
}

extern "C" int32_t caps_provenance_file(const char* path, char* json, int32_t cap) {
  caps::Json j = caps::Json::object();
  if (auto m = caps::read_manifest(path ? path : "")) {
    j = caps::manifest_json(*m);
    j["ok"] = true;
  } else {
    j["ok"] = false;
    j["error"] = std::string("no provenance beside ") + (path ? path : "") + " (" + caps::sidecar_path(path ? path : "") + ")";
  }
  return report_out(j.dump(0), json, cap);
}

extern "C" int32_t caps_provenance_compare(const char* a_json, const char* b_json, char* json, int32_t cap) {
  caps::Json j = caps::Json::object();
  try {
    const auto a = caps::manifest_from_json(caps::Json::parse(a_json && *a_json ? a_json : "{}"));
    const auto b = caps::manifest_from_json(caps::Json::parse(b_json && *b_json ? b_json : "{}"));
    const auto d = caps::compare(a, b);
    j["ok"] = true;
    j["same_inputs"] = d.same_inputs;
    j["same_generator"] = d.same_generator;
    j["steps_a"] = d.steps_a;
    j["steps_b"] = d.steps_b;
    j["differing_steps"] = d.differing_steps;
    caps::Json rows = caps::Json::array();
    for (const auto& r : d.rows) {
      caps::Json o = caps::Json::object();
      o["step"] = r.step; o["engine"] = r.engine; o["key"] = r.key; o["a"] = r.a; o["b"] = r.b;
      rows.push_back(std::move(o));
    }
    j["rows"] = std::move(rows);
    caps::Json notes = caps::Json::array();
    for (const auto& n : d.notes) notes.push_back(n);
    j["notes"] = std::move(notes);
  } catch (const std::exception& e) {
    j = caps::Json::object();
    j["ok"] = false;
    j["error"] = std::string(e.what());
  }
  return report_out(j.dump(0), json, cap);
}

extern "C" int32_t caps_provenance_bibtex(const char* manifest_json, char* text, int32_t cap) {
  try {
    const auto m = caps::manifest_from_json(caps::Json::parse(manifest_json && *manifest_json ? manifest_json : "{}"));
    return report_out(caps::bibtex(caps::all_cites(m)), text, cap);
  } catch (const std::exception& e) {
    g_error = e.what();
    return -1;
  }
}

// Coherent neutron scattering length (fm) of an element (1001: ²H), NaN when CAPS has none.
extern "C" double caps_neutron_b(int32_t z) {
  try { return caps::neutron_b(z); } catch (...) { return std::numeric_limits<double>::quiet_NaN(); }
}

// ---------------------------------------------------------------- voids (design/boards/FreeVolume)
extern "C" int32_t caps_voids(caps_doc* d, const char* options_json, char* json, int32_t cap) {
  caps::Json j = caps::Json::object();
  try {
    const caps::Json o = caps::Json::parse(options_json && *options_json ? options_json : "{}");
    const bool show = !o.has("show") || (o["show"].kind() == caps::Json::Bool ? o["show"].boolean() : o["show"].number() != 0);
    if (o.has("clear") && o["clear"].kind() == caps::Json::Bool && o["clear"].boolean()) {
      d->voids.clear();
      d->void_mesh.reset();
      j["ok"] = true;
      return report_out(j.dump(0), json, cap);
    }
    caps::VoidOptions vo;
    vo.grid = std::clamp(o.num("grid", 0.5), 0.2, 2.0);
    vo.probe = std::max(0.0, o.num("probe", 1.4));
    vo.max_count = std::clamp(int(o.num("count", 40)), 1, 2000);
    vo.min_radius = std::max(0.2, o.num("min_radius", 1.0));
    const auto r = caps::largest_voids(d->frame, vo);
    d->voids = r.spheres;
    if (show) d->void_mesh = std::make_unique<caps::Mesh>(caps::void_mesh(r.spheres, 2));
    else d->void_mesh.reset();
    j["ok"] = true;
    j["accessible_point"] = r.accessible_point;
    j["accessible_probe"] = r.accessible_probe;
    j["largest"] = r.largest;
    caps::Json g = caps::Json::array();
    for (int k = 0; k < 3; ++k) g.push_back(r.grid[k]);
    j["grid"] = std::move(g);
    caps::Json sp = caps::Json::array();
    for (const auto& v : r.spheres) {
      caps::Json x = caps::Json::object();
      x["x"] = v.centre[0]; x["y"] = v.centre[1]; x["z"] = v.centre[2]; x["r"] = v.radius;
      sp.push_back(std::move(x));
    }
    j["spheres"] = std::move(sp);
  } catch (const std::exception& e) {
    j = caps::Json::object();
    j["ok"] = false;
    j["error"] = std::string(e.what());
  }
  return report_out(j.dump(0), json, cap);
}

extern "C" int32_t caps_voids_pdb(caps_doc* d, const char* path) {
  return guard([&] {
    if (d->voids.empty()) throw std::runtime_error("no voids computed yet");
    caps::write_voids_pdb(d->frame, d->voids, path);
    return 0;
  });
}

// ---------------------------------------------------------------- pores (design/boards/SlitPore)
extern "C" caps_doc* caps_pore_build(const char* options_json, char* report, int32_t cap) {
  try {
    const caps::Json j = caps::Json::parse(options_json && *options_json ? options_json : "{}");
    auto flag = [&](const char* k, bool def) { return !j.has(k) ? def : j[k].kind() == caps::Json::Bool ? j[k].boolean() : j[k].number() != 0; };
    caps::PoreOptions o;
    const std::string kind = j.text("kind", "slit");
    o.kind = kind == "cylinder" ? caps::PoreKind::Cylinder : kind == "framework" ? caps::PoreKind::Framework : caps::PoreKind::Slit;
    o.width = j.num("width", 10.0);
    o.layers = int(j.num("layers", 1));
    o.lx = j.num("lx", 26.0);
    o.ly = j.num("ly", 22.0);
    o.vacuum = flag("vacuum", false);
    o.vacuum_gap = j.num("vacuum_gap", 20.0);
    o.wall = j.num("wall", 6.0);
    o.length = j.num("length", 20.0);
    o.passivate = flag("passivate", false);
    if (j.has("repeat") && j["repeat"].is_array() && j["repeat"].size() == 3)
      for (int k = 0; k < 3; ++k) o.repeat[k] = int(j["repeat"][size_t(k)].number());
    o.tolerance = j.num("tolerance", 2.0);
    o.seed = uint64_t(j.num("seed", 1));
    caps::System crystal, fluid;
    const std::string cif = j.text("cif");
    if (o.kind != caps::PoreKind::Slit) {
      if (cif.empty()) throw std::invalid_argument("choose a crystal (CIF) for a cylindrical or framework pore");
      crystal = caps::read_cif(cif);
      o.crystal = &crystal;
    }
    const std::string smiles = j.text("fluid");
    o.count = int(j.num("count", 0));
    if (!smiles.empty() && o.count > 0) {
      caps::BuildOptions bo;
      bo.forcefield = "uff";
      fluid = caps::build_molecule(smiles, bo).system;
      fluid.title = j.text("fluid_name", smiles);
      o.fluid = &fluid;
    }
    caps::PoreReport r;
    const caps::System s = caps::build_pore(o, &r);
    caps::Json rj = caps::Json::object();
    rj["wall_atoms"] = r.wall_atoms;
    rj["fluid_molecules"] = r.fluid_molecules;
    rj["width"] = r.width;
    rj["pore_volume"] = r.pore_volume;
    rj["fluid_density"] = r.fluid_density;
    rj["dmin"] = r.dmin;
    caps::Json notes = caps::Json::array();
    for (const auto& n : r.notes) notes.push_back(n);
    rj["notes"] = std::move(notes);
    report_out(rj.dump(0), report, cap);
    caps_doc* d = doc_of(s);
    if (!cif.empty()) prov_input(d, cif);
    prov_step(d, "nano.pore", kind + " pore" + (o.fluid ? " with " + std::to_string(o.count) + " × " + smiles : std::string()), json_params(options_json),
              o.fluid ? seeded(o.seed) : "", o.fluid ? std::vector<std::string>{"martinez2009", "rappe1992"} : std::vector<std::string>{});
    return d;
  } catch (const std::exception& e) {
    g_error = e.what();
    return nullptr;
  }
}

// A readable reference for a built-in citation key (the theory manual).
extern "C" int32_t caps_citation_text(const char* key, char* text, int32_t cap) {
  return report_out(caps::citation_text(key ? key : ""), text, cap);
}

// A methods paragraph from a manifest (and manifests of replicas that differ only in their seeds): {text, refs[]}.
extern "C" int32_t caps_methods_text(const char* manifest_json, const char* replicas_json, char* json, int32_t cap) {
  caps::Json j = caps::Json::object();
  try {
    const auto m = caps::manifest_from_json(caps::Json::parse(manifest_json && *manifest_json ? manifest_json : "{}"));
    std::vector<caps::Manifest> reps;
    if (replicas_json && *replicas_json) {
      const caps::Json r = caps::Json::parse(replicas_json);
      if (r.is_array()) for (const auto& x : r.items()) reps.push_back(caps::manifest_from_json(x));
    }
    std::vector<std::string> refs;
    j["text"] = caps::methods_text(m, &refs, reps);
    caps::Json rl = caps::Json::array();
    for (const auto& x : refs) rl.push_back(x);
    j["refs"] = std::move(rl);
    j["ok"] = true;
  } catch (const std::exception& e) {
    j = caps::Json::object();
    j["ok"] = false;
    j["error"] = std::string(e.what());
  }
  return report_out(j.dump(0), json, cap);
}

// ---------------------------------------------------------------- coarse-grained melts (design/boards/CoarseGrained)
namespace {
caps::KgOptions kg_options(const char* options_json) {
  const caps::Json j = caps::Json::parse(options_json && *options_json ? options_json : "{}");
  caps::KgOptions o;
  o.chains = int(j.num("chains", 50));
  o.beads = int(j.num("beads", 100));
  o.density = j.num("density", 0.85);
  o.k_theta = j.num("k_theta", 0.0);
  o.seed = uint64_t(j.num("seed", 1));
  return o;
}
}  // namespace

extern "C" caps_doc* caps_kg_build(const char* options_json, char* report, int32_t cap) {
  try {
    const caps::KgOptions o = kg_options(options_json);
    caps::KgReport r;
    const caps::System s = caps::kremer_grest(o, &r);
    caps::Json j = caps::Json::object();
    j["box"] = r.box;
    j["closest"] = r.closest;
    j["r2_per_bond"] = r.mean_r2;
    report_out(j.dump(0), report, cap);
    caps_doc* d = doc_of(s);
    prov_step(d, "cg.kremer_grest", std::to_string(o.chains) + " × " + std::to_string(o.beads) + " bead-spring chains as random walks",
              {{"chains", std::to_string(o.chains)}, {"beads", std::to_string(o.beads)}, {"density", g6(o.density) + " σ⁻³"}, {"k_theta", g6(o.k_theta) + " ε"},
               {"box", g6(r.box) + " σ"}}, seeded(o.seed), {"kremer1990"});
    return d;
  } catch (const std::exception& e) {
    g_error = e.what();
    return nullptr;
  }
}

extern "C" int32_t caps_kg_lammps(caps_doc* d, const char* options_json, const char* stem, double pushoff_steps, double run_steps) {
  return guard([&] {
    caps::write_kg_lammps(d->frame, kg_options(options_json), stem, pushoff_steps > 0 ? pushoff_steps : 20000, run_steps > 0 ? run_steps : 100000);
    return 0;
  });
}

// ---------------------------------------------------------------- colour vision (design/boards/ColourVision)
extern "C" int32_t caps_vision_check(const char* palettes_json, double threshold, char* json, int32_t cap) {
  return guard([&] {
    const caps::Json in = caps::Json::parse(palettes_json ? palettes_json : "{}");
    std::vector<caps::NamedPalette> pals;
    for (const auto& [name, p] : in.members()) {
      caps::NamedPalette np;
      np.name = name;
      for (const auto& c : p["colours"].items()) {
        std::string h = c.str();
        if (!h.empty() && h[0] == '#') h = h.substr(1);
        np.colours.push_back(unsigned(std::stoul(h.substr(0, 6), nullptr, 16)));
      }
      if (p.has("labels")) for (const auto& l : p["labels"].items()) np.labels.push_back(l.str());
      while (np.labels.size() < np.colours.size()) np.labels.push_back(std::to_string(np.labels.size() + 1));
      pals.push_back(std::move(np));
    }
    auto hex = [](unsigned c) { char b[8]; std::snprintf(b, sizeof b, "#%06X", c & 0xFFFFFF); return std::string(b); };
    caps::Json out = caps::Json::object(), pa = caps::Json::array();
    for (const auto& p : pals) {
      caps::Json j = caps::Json::object(), labels = caps::Json::array();
      j["name"] = p.name;
      for (const auto& l : p.labels) labels.push_back(caps::Json(l));
      j["labels"] = std::move(labels);
      for (caps::Vision v : {caps::Vision::Normal, caps::Vision::Protan, caps::Vision::Deutan, caps::Vision::Tritan}) {
        caps::Json cs = caps::Json::array();
        for (unsigned c : p.colours) cs.push_back(caps::Json(hex(caps::simulate_vision(c, v))));
        j[caps::to_string(v)] = std::move(cs);
      }
      pa.push_back(std::move(j));
    }
    out["palettes"] = std::move(pa);
    caps::Json pairs = caps::Json::array();
    for (const auto& q : caps::confusable_pairs(pals, threshold > 0 ? threshold : 12.0)) {
      caps::Json j = caps::Json::object();
      const auto& p = *std::find_if(pals.begin(), pals.end(), [&](const caps::NamedPalette& x) { return x.name == q.palette; });
      j["palette"] = q.palette;
      j["vision"] = std::string(caps::to_string(q.vision));
      j["a"] = p.labels[size_t(q.a)];
      j["b"] = p.labels[size_t(q.b)];
      j["de"] = std::round(q.de * 10) / 10;
      pairs.push_back(std::move(j));
    }
    out["pairs"] = std::move(pairs);
    return report_out(out.dump(0), json, cap);
  });
}

extern "C" int32_t caps_periodic(caps_doc* d, const char* json, char* out, int32_t cap) {
  caps::Json r = caps::Json::object();
  try {
    const caps::Json j = caps::Json::parse(json && *json ? json : "{}");
    const caps::System& f = d->frame;
    if (!f.cell.valid()) throw std::runtime_error("the structure has no periodic cell");
    caps::System whole = f;
    caps::make_molecules_whole(whole);
    const auto mol = whole.molecules();
    const size_t n = whole.atoms.size();
    auto wrapped = [&](const caps::Vec3& p) {
      caps::Vec3 fr = whole.cell.to_fractional(p);
      for (int k = 0; k < 3; ++k) fr[k] -= std::floor(fr[k]);
      return whole.cell.to_cartesian(fr);
    };
    std::vector<caps::Vec3> w(n);
    for (size_t i = 0; i < n; ++i) w[i] = wrapped(whole.atoms[i].pos);
    // crossing: any atom of the whole molecule outside the cell; pieces: components of the wrapped frame whose bonds
    // are shorter than half the smallest cell width
    std::map<int, bool> crosses;
    for (size_t i = 0; i < n; ++i) {
      const caps::Vec3 fr = whole.cell.to_fractional(whole.atoms[i].pos);
      bool out = false;
      for (int k = 0; k < 3; ++k) out = out || fr[k] < -1e-9 || fr[k] >= 1.0 - 1e-12;
      crosses[mol[i]] = crosses[mol[i]] || out;
    }
    const double half = 0.5 * std::min({caps::norm(whole.cell.a), caps::norm(whole.cell.b), caps::norm(whole.cell.c)});
    std::vector<int> parent(n);
    std::iota(parent.begin(), parent.end(), 0);
    std::function<int(int)> find = [&](int x) { return parent[size_t(x)] == x ? x : parent[size_t(x)] = find(parent[size_t(x)]); };
    for (const auto& b : whole.bonds)
      if (caps::norm(w[b.i] - w[b.j]) < half) parent[size_t(find(int(b.i)))] = find(int(b.j));
    std::set<int> roots;
    for (size_t i = 0; i < n; ++i) roots.insert(find(int(i)));
    int crossing = 0;
    for (const auto& [m, c] : crosses) crossing += c ? 1 : 0;
    r["ok"] = true;
    r["molecules"] = double(crosses.size());
    r["crossing"] = double(crossing);
    r["pieces"] = double(roots.size());
    caps::Json box = caps::Json::array();
    for (const auto& v : {whole.cell.a, whole.cell.b, whole.cell.c}) box.push_back(caps::Json(caps::norm(v)));
    r["box"] = std::move(box);
    const double la = caps::norm(whole.cell.a), lb = caps::norm(whole.cell.b), lc = caps::norm(whole.cell.c);
    auto ang = [](const caps::Vec3& u, const caps::Vec3& v) { return std::acos(std::clamp(caps::dot(u, v) / (caps::norm(u) * caps::norm(v)), -1.0, 1.0)) * 180 / M_PI; };
    caps::Json angles = caps::Json::array();
    angles.push_back(caps::Json(ang(whole.cell.b, whole.cell.c)));
    angles.push_back(caps::Json(ang(whole.cell.a, whole.cell.c)));
    angles.push_back(caps::Json(ang(whole.cell.a, whole.cell.b)));
    r["angles"] = std::move(angles);
    r["cubic"] = std::fabs(la - lb) < 1e-3 && std::fabs(la - lc) < 1e-3 && std::fabs(caps::dot(whole.cell.a, whole.cell.b)) < 1e-6 && std::fabs(caps::dot(whole.cell.a, whole.cell.c)) < 1e-6;
    // one molecule, measured three ways
    long want = long(j.num("molecule", 0));
    std::vector<int> ids;
    for (const auto& [m, c] : crosses) ids.push_back(m);
    int m = ids.empty() ? -1 : ids.front();
    if (want > 0 && size_t(want) <= ids.size()) m = ids[size_t(want - 1)];
    else for (const auto& [mm, c] : crosses) if (c) { m = mm; break; }
    r["molecule"] = double(std::find(ids.begin(), ids.end(), m) - ids.begin() + 1);
    auto three = [&](uint32_t a, uint32_t b) {
      caps::Json e = caps::Json::object();
      e["i"] = double(a + 1), e["j"] = double(b + 1);
      e["wrapped"] = caps::norm(w[a] - w[b]);
      e["min_image"] = caps::norm(whole.cell.minimum_image(whole.atoms[a].pos - whole.atoms[b].pos));
      e["whole"] = caps::norm(whole.atoms[a].pos - whole.atoms[b].pos);
      return e;
    };
    for (const auto& b : whole.bonds)
      if (mol[b.i] == m && caps::norm(w[b.i] - w[b.j]) >= half) { r["bond"] = three(b.i, b.j); break; }
    uint32_t first = uint32_t(n), last = 0;
    for (uint32_t i = 0; i < n; ++i) if (mol[i] == m) { first = std::min(first, i); last = std::max(last, i); }
    // chain ends: the heavy atoms at the two ends of the molecule's index range
    while (first < last && whole.atoms[first].element == 1) ++first;
    while (last > first && whole.atoms[last].element == 1) --last;
    if (first < last) r["ends"] = three(first, last);
  } catch (const std::exception& e) {
    r = caps::Json::object();
    r["ok"] = false;
    r["error"] = std::string(e.what());
  }
  return report_out(r.dump(0), out, cap);
}

extern "C" int32_t caps_set_images(caps_doc* d, int32_t na, int32_t nb, int32_t nc, double fade) {
  if (!d) return -1;
  d->images = {std::clamp(na, 1, 5), std::clamp(nb, 1, 5), std::clamp(nc, 1, 5)};
  d->image_fade = float(std::clamp(fade, 0.0, 1.0));
  return 0;
}

extern "C" int32_t caps_set_save_wrap(caps_doc* d, int32_t mode) {
  if (!d || mode < 0 || mode > 2) return -1;
  d->save_wrap = mode;
  return 0;
}

extern "C" int32_t caps_centre_on(caps_doc* d, const int32_t* idx, int32_t n) {
  return guard([&] {
    const caps::System& f = d->frame;
    if (!f.cell.valid()) throw std::runtime_error("the structure has no periodic cell");
    caps::Vec3 c{0, 0, 0};
    int used = 0;
    for (int32_t k = 0; k < n; ++k)
      if (idx[k] >= 0 && size_t(idx[k]) < f.atoms.size()) { c = c + f.atoms[size_t(idx[k])].pos; ++used; }
    if (!used) throw std::invalid_argument("centre on: no atoms given");
    c = c * (1.0 / used);
    const caps::Vec3 mid = f.cell.origin + (f.cell.a + f.cell.b + f.cell.c) * 0.5;
    push_undo(d, "Centre on " + std::to_string(used) + " atoms");
    for (auto& p : d->traj.positions.at(d->current)) p = p + (mid - c);
    refresh(d);
    return 0;
  });
}

// Cell editor (design/boards/CellEditor): new parameters, with the atoms scaled (fractional coordinates kept) or left
// where they are; a supercell replicates atoms and bonds.
extern "C" int32_t caps_set_cell(caps_doc* d, const char* json) {
  return guard([&] {
    const caps::Json j = caps::Json::parse(json && *json ? json : "{}");
    const caps::Cell nc = caps::cell_parameters(j.num("a", 10), j.num("b", 10), j.num("c", 10), j.num("alpha", 90), j.num("beta", 90), j.num("gamma", 90));
    if (!(nc.volume() > 1e-6)) throw std::invalid_argument("these parameters give no volume");
    const bool scale = !j.has("scale") || j["scale"].kind() != caps::Json::Bool || j["scale"].boolean();
    const caps::Cell old = d->traj.topology.cell;
    push_undo(d, "Cell " + g6(j.num("a", 10)) + " × " + g6(j.num("b", 10)) + " × " + g6(j.num("c", 10)) + (scale ? " (atoms scaled)" : ""));
    caps::Cell c = nc;
    c.origin = old.valid() ? old.origin : caps::Vec3{0, 0, 0};
    if (scale && old.valid())
      for (auto& p : d->traj.positions.at(d->current)) p = c.to_cartesian(old.to_fractional(p));
    d->traj.topology.cell = c;
    for (auto& cc : d->traj.cells) cc = c;
    d->field.reset();
    refresh(d);
    return 0;
  });
}

extern "C" int32_t caps_supercell(caps_doc* d, int32_t na, int32_t nb, int32_t nc) {
  return guard([&] {
    if (na < 1 || nb < 1 || nc < 1 || na * nb * nc > 1000) throw std::invalid_argument("supercell: 1 … 1000 copies");
    const caps::Cell cell = d->traj.topology.cell;
    if (!cell.valid()) throw std::runtime_error("the structure has no cell to repeat");
    push_undo(d, "Supercell " + std::to_string(na) + " × " + std::to_string(nb) + " × " + std::to_string(nc));
    caps::System top = d->traj.topology;
    const auto pos = d->traj.positions.at(d->current);
    const size_t n = top.atoms.size();
    int64_t maxmol = 0;
    for (const auto& a : top.atoms) maxmol = std::max(maxmol, a.mol);
    std::vector<caps::Vec3> out;
    caps::System s = top;
    s.atoms.clear();
    s.bonds.clear();
    int copy = 0;
    for (int a = 0; a < na; ++a)
      for (int b = 0; b < nb; ++b)
        for (int c = 0; c < nc; ++c, ++copy) {
          const caps::Vec3 t = cell.a * double(a) + cell.b * double(b) + cell.c * double(c);
          const uint32_t off = uint32_t(copy * n);
          for (size_t i = 0; i < n; ++i) {
            caps::Atom at = top.atoms[i];
            at.id = int64_t(off + i + 1);
            if (at.mol > 0) at.mol += copy * maxmol;
            s.atoms.push_back(at);
            out.push_back(pos[i] + t);
          }
          for (const auto& bd : top.bonds) s.bonds.push_back({bd.i + off, bd.j + off, bd.order});
        }
    s.cell.a = cell.a * double(na);
    s.cell.b = cell.b * double(nb);
    s.cell.c = cell.c * double(nc);
    d->traj.topology = s;
    d->traj.positions = {out};
    d->traj.cells = {s.cell};
    d->traj.timesteps = {0};
    d->current = 0;
    d->field.reset();
    refresh(d);
    return 0;
  });
}

// Two states of the same atoms compared (design/boards/Compare): the moving one superposed on the reference (Horn 1987),
// RMSD over all, heavy and backbone atoms, the shift of each atom and the largest ones; colour: the view shows the shift.
extern "C" int32_t caps_provenance_note(caps_doc* d, const char* json) {
  return guard([&] {
    if (!d) throw std::runtime_error("no document");
    const caps::Json o = caps::Json::parse(json && *json ? json : "{}");
    const std::string engine = o.text("engine");
    if (engine.empty() || engine.find('.') == std::string::npos) throw std::runtime_error("a provenance note needs an engine such as equilibrate.accepted");
    caps::KeyValues kv;
    if (o.has("params") && o["params"].kind() == caps::Json::Object)
      for (const auto& [k, v] : o["params"].members()) kv.push_back({k, v.kind() == caps::Json::String ? v.str() : v.dump(0)});
    prov_step(d, engine, o.text("summary"), kv);
    return 0;
  });
}

extern "C" int32_t caps_checkpoint(caps_doc* d, const char* json, char* out, int32_t cap) {
  caps::Json r = caps::Json::object();
  try {
    if (!d) throw std::runtime_error("no document");
    const caps::Json o = caps::Json::parse(json && *json ? json : "{}");
    const std::string op = o.text("op", "info");
    auto& k = d->checkpoint;
    if (op == "clear") k = {};
    else if (op == "restore") {
      if (!k.has) throw std::runtime_error("no checkpoint to continue from");
      if (k.x.size() != d->traj.topology.atoms.size()) throw std::runtime_error("the structure changed since the checkpoint");
      d->traj.positions.push_back(k.x);
      d->traj.cells.push_back(k.cell);
      d->traj.timesteps.push_back(k.step);
      d->traj.topology.velocities = k.v;
      d->traj.topology.unwrapped = true;
      d->current = d->traj.frames() - 1;
      refresh(d);
      prov_step(d, "md.checkpoint", "continued from the checkpoint at step " + std::to_string(k.step) + " (the run " + k.ended + ")",
                {{"step", std::to_string(k.step) + " of " + std::to_string(k.steps)}, {"time", g6(double(k.step) * k.dt / 1000.0) + " ps"}});
    } else if (op != "info") throw std::runtime_error("checkpoint op info, restore or clear");
    r["ok"] = true;
    r["has"] = k.has;
    if (k.has) {
      r["kind"] = k.kind;
      r["step"] = double(k.step);
      r["steps"] = double(k.steps);
      r["time_ps"] = double(k.step) * k.dt / 1000.0;
      r["atoms"] = double(k.x.size());
      r["ended"] = k.ended.empty() ? std::string("running") : k.ended;
      if (!k.error.empty()) r["reason"] = k.error;
    }
  } catch (const std::exception& e) {
    r = caps::Json::object();
    r["ok"] = false;
    r["error"] = e.what();
  }
  return report_out(r.dump(0), out, cap);
}

// The Properties explorer's data (as Materials Studio's): the structure (formula, composition, mass, charge, cell with
// its angles, frames, force field) and one atom (types, charge, position, fractional position, bonded neighbours).
extern "C" int32_t caps_structure_info(caps_doc* d, char* out, int32_t cap) {
  caps::Json r = caps::Json::object();
  try {
    if (!d) throw std::runtime_error("no document");
    const caps::System& s = d->frame;
    std::map<int, size_t> count;
    double mass = 0, q = 0;
    for (const auto& a : s.atoms) {
      ++count[a.element];
      mass += s.mass_of(a);
      q += a.charge;
    }
    // Hill order: C, H, then the rest alphabetically
    std::vector<std::pair<std::string, size_t>> parts;
    for (const auto& [e, n] : count) parts.push_back({e > 0 ? caps::element(e).symbol : "X", n});
    auto rank = [&](const std::string& x) { return count.count(6) ? (x == "C" ? 0 : x == "H" ? 1 : 2) : 2; };
    std::sort(parts.begin(), parts.end(), [&](const auto& a, const auto& b) { return rank(a.first) != rank(b.first) ? rank(a.first) < rank(b.first) : a.first < b.first; });
    std::string formula;
    caps::Json comp = caps::Json::object();
    for (const auto& [sym, n] : parts) {
      formula += sym + (n > 1 ? std::to_string(n) : "");
      comp[sym] = double(n);
    }
    r["formula"] = formula;
    r["composition"] = comp;
    r["title"] = s.title;
    r["format"] = s.source_format;
    r["atoms"] = double(s.atoms.size());
    r["bonds"] = double(s.bonds.size());
    int nmol = 0;
    s.molecules(&nmol);
    r["molecules"] = double(nmol);
    r["frames"] = double(std::max<size_t>(1, d->traj.frames()));
    r["frame"] = double(d->current);
    r["timestep"] = double(s.timestep);
    r["mass"] = mass;
    r["charge"] = q;
    r["has_charges"] = s.has_charges;
    r["unwrapped"] = s.unwrapped;
    r["bonds_from_file"] = s.bonds_from_file;
    if (s.cell.valid()) {
      const auto& c = s.cell;
      const double la = caps::norm(c.a), lb = caps::norm(c.b), lc = caps::norm(c.c);
      auto ang = [](const caps::Vec3& u, const caps::Vec3& v) { return std::acos(std::clamp(caps::dot(u, v) / (caps::norm(u) * caps::norm(v)), -1.0, 1.0)) * 180.0 / M_PI; };
      caps::Json cell = caps::Json::object();
      cell["a"] = la, cell["b"] = lb, cell["c"] = lc;
      cell["alpha"] = ang(c.b, c.c), cell["beta"] = ang(c.a, c.c), cell["gamma"] = ang(c.a, c.b);
      cell["volume"] = c.volume();
      cell["density"] = s.density();
      caps::Json per = caps::Json::array();
      for (bool p : c.periodic) per.push_back(p);
      cell["periodic"] = per;
      r["cell"] = cell;
    }
    if (d->field) {
      caps::Json f = caps::Json::object();
      f["name"] = d->field->base.name.empty() ? std::filesystem::path(d->field->ff_path).stem().string() : d->field->base.name;
      f["complete"] = d->field->complete;
      f["missing"] = double(d->field->rep.missing.size());
      std::set<std::string> ts(d->field->types.begin(), d->field->types.end());
      f["types"] = double(ts.size());
      r["forcefield"] = f;
    }
    size_t sel = 0;
    for (char c : d->selection) sel += c != 0;
    r["selected"] = double(sel);
    r["ok"] = true;
  } catch (const std::exception& e) {
    r = caps::Json::object();
    r["ok"] = false;
    r["error"] = e.what();
  }
  return report_out(r.dump(0), out, cap);
}

extern "C" int32_t caps_atom_properties(caps_doc* d, int32_t index, char* out, int32_t cap) {
  caps::Json r = caps::Json::object();
  try {
    if (!d) throw std::runtime_error("no document");
    const caps::System& s = d->frame;
    if (index < 0 || size_t(index) >= s.atoms.size()) throw std::runtime_error("no atom " + std::to_string(index + 1));
    const auto& a = s.atoms[size_t(index)];
    r["index"] = double(index);
    r["id"] = double(a.id);
    r["element"] = a.element > 0 ? std::string(caps::element(a.element).symbol) : std::string("X");
    r["name"] = a.name;
    r["type"] = double(a.type);
    if (d->field && size_t(index) < d->field->types.size()) r["ff_type"] = d->field->types[size_t(index)];
    r["charge"] = a.charge;
    r["mass"] = s.mass_of(a);
    r["molecule"] = double(a.mol);
    if (a.resid) r["residue"] = double(a.resid);
    if (!a.resname.empty()) r["resname"] = a.resname;
    caps::Json xyz = caps::Json::array();
    for (double v : a.pos) xyz.push_back(v);
    r["xyz"] = xyz;
    if (s.cell.valid()) {
      const caps::Vec3 f = s.cell.to_fractional(a.pos);
      caps::Json fr = caps::Json::array();
      for (double v : f) fr.push_back(v);
      r["fractional"] = fr;
    }
    caps::Json nb = caps::Json::array();
    for (const auto& b : s.bonds) {
      if (b.i != uint32_t(index) && b.j != uint32_t(index)) continue;
      const uint32_t o = b.i == uint32_t(index) ? b.j : b.i;
      caps::Vec3 dv = s.atoms[o].pos - a.pos;
      if (s.cell.valid()) dv = s.cell.minimum_image(dv);
      caps::Json e = caps::Json::object();
      e["index"] = double(o);
      e["element"] = s.atoms[o].element > 0 ? std::string(caps::element(s.atoms[o].element).symbol) : std::string("X");
      e["distance"] = caps::norm(dv);
      e["order"] = double(b.order);
      nb.push_back(e);
    }
    r["neighbours"] = nb;
    r["ok"] = true;
  } catch (const std::exception& e) {
    r = caps::Json::object();
    r["ok"] = false;
    r["error"] = e.what();
  }
  return report_out(r.dump(0), out, cap);
}

extern "C" int32_t caps_compare_states(caps_doc* d, const char* json, char* out, int32_t cap) {
  caps::Json r = caps::Json::object();
  try {
    if (!d) throw std::runtime_error("no document");
    const caps::Json o = caps::Json::parse(json && *json ? json : "{}");
    if (o.text("op") == "clear") {
      d->atom_values.clear();
      r["ok"] = true;
      return report_out(r.dump(0), out, cap);
    }
    // a state: {kind: current | start | frame | snapshot, index}
    auto state = [&](const std::string& key, const std::string& def) -> std::pair<std::vector<caps::Vec3>, std::string> {
      caps::Json st = o.has(key) ? o[key] : caps::Json::object();
      const std::string kind = st.kind() == caps::Json::Object ? st.text("kind", def) : def;
      const int idx = st.kind() == caps::Json::Object ? int(st.num("index", 0)) : 0;
      std::vector<caps::Vec3> p;
      if (kind == "current") {
        for (const auto& a : d->frame.atoms) p.push_back(a.pos);
        return {p, "current"};
      }
      if (kind == "start") {
        if (!d->undo.empty()) return {d->undo.front().positions, "before the first edit"};
        if (d->traj.frames()) return {d->traj.positions.front(), "frame 1"};
        for (const auto& a : d->frame.atoms) p.push_back(a.pos);
        return {p, "current"};
      }
      if (kind == "frame") {
        if (idx < 0 || size_t(idx) >= d->traj.frames()) throw std::runtime_error("no frame " + std::to_string(idx + 1));
        return {d->traj.positions[size_t(idx)], "frame " + std::to_string(idx + 1)};
      }
      if (kind == "snapshot") {
        if (idx < 0 || size_t(idx) >= d->snapshots.size()) throw std::runtime_error("no snapshot " + std::to_string(idx + 1));
        return {d->snapshots[size_t(idx)].state.positions, d->snapshots[size_t(idx)].name};
      }
      throw std::runtime_error("unknown state '" + kind + "' (current, start, frame, snapshot)");
    };
    auto [ref, ref_name] = state("reference", "start");
    auto [mov, mov_name] = state("moving", "current");
    const auto& atoms = d->frame.atoms;
    if (ref.size() != mov.size()) throw std::runtime_error("the two states have different atoms (" + std::to_string(ref.size()) + " and " + std::to_string(mov.size()) + "): compare states of the same structure");
    if (ref.size() != atoms.size()) throw std::runtime_error("the states have " + std::to_string(ref.size()) + " atoms, the structure shown " + std::to_string(atoms.size()));
    const size_t n = ref.size();
    // periodic cells: each atom's nearest image to its reference position (frames may be wrapped)
    if (d->frame.cell.valid() && o.text("periodic", "yes") != "no")
      for (size_t i = 0; i < n; ++i) mov[i] = ref[i] + d->frame.cell.minimum_image(mov[i] - ref[i]);
    std::vector<char> heavy(n), back(n, 0);
    for (size_t i = 0; i < n; ++i) heavy[i] = atoms[i].element != 1;
    for (const auto& chain : caps::backbones(d->frame))
      for (uint32_t i : chain) back[i] = 1;
    size_t nback = 0;
    for (char b : back) nback += b;
    const std::string fit_on = o.text("fit", "all");
    std::vector<double> w(n, 1.0);
    if (fit_on == "heavy") for (size_t i = 0; i < n; ++i) w[i] = heavy[i] ? 1.0 : 0.0;
    else if (fit_on == "backbone") {
      if (!nback) throw std::runtime_error("no backbone found to fit on");
      for (size_t i = 0; i < n; ++i) w[i] = back[i] ? 1.0 : 0.0;
    } else if (fit_on == "selection") {
      if (d->selection.size() != n) throw std::runtime_error("select the atoms to fit on first");
      for (size_t i = 0; i < n; ++i) w[i] = d->selection[i] ? 1.0 : 0.0;
    } else if (fit_on == "none") {
      w.clear();
    } else if (fit_on != "all") throw std::runtime_error("fit on all, heavy, backbone, selection or none");
    caps::Superposition fit;
    if (!w.empty() || fit_on == "none") {
      if (fit_on == "none") fit = caps::Superposition{};
      else fit = caps::superpose(ref, mov, w);
    }
    std::vector<double> shift(n);
    for (size_t i = 0; i < n; ++i) shift[i] = caps::norm(fit.apply(mov[i]) - ref[i]);
    caps::Json rm = caps::Json::object();
    rm["all"] = caps::rmsd_after(fit, ref, mov);
    rm["heavy"] = caps::rmsd_after(fit, ref, mov, heavy);
    if (nback) rm["backbone"] = caps::rmsd_after(fit, ref, mov, back);
    r["rmsd"] = rm;
    r["fit"] = fit_on;
    r["fitted"] = double(fit.fitted);
    r["atoms"] = double(n);
    r["reference"] = ref_name;
    r["moving"] = mov_name;
    double mx = 0, mean = 0;
    for (double x : shift) mx = std::max(mx, x), mean += x;
    r["max"] = mx;
    r["mean"] = n ? mean / double(n) : 0.0;
    std::vector<size_t> order(n);
    for (size_t i = 0; i < n; ++i) order[i] = i;
    const size_t k = std::min(n, size_t(std::max(0.0, o.num("largest", 20))));
    std::partial_sort(order.begin(), order.begin() + k, order.end(), [&](size_t a, size_t b) { return shift[a] > shift[b]; });
    caps::Json big = caps::Json::array();
    for (size_t q = 0; q < k; ++q) {
      const size_t i = order[q];
      caps::Json e = caps::Json::object();
      e["atom"] = double(i);
      e["label"] = std::string(caps::element(atoms[i].element).symbol) + std::to_string(i + 1);
      e["element"] = std::string(caps::element(atoms[i].element).symbol);
      e["molecule"] = double(atoms[i].mol);
      e["backbone"] = bool(back[i]);
      e["shift"] = shift[i];
      big.push_back(e);
    }
    r["largest"] = big;
    if (o.num("per_atom", 0) > 0) {
      caps::Json all = caps::Json::array();
      for (double x : shift) all.push_back(x);
      r["shifts"] = all;
    }
    if (o.num("colour", 0) > 0) {
      d->atom_values = shift;
      d->atom_values_ramp = 0;
    }
    r["ok"] = true;
  } catch (const std::exception& e) {
    r = caps::Json::object();
    r["ok"] = false;
    r["error"] = e.what();
  }
  return report_out(r.dump(0), out, cap);
}

extern "C" int32_t caps_set_atom_values(caps_doc* d, const double* values, int32_t n, int32_t ramp) {
  if (!d) return -1;
  d->atom_values.assign(values && n > 0 ? values : nullptr, values && n > 0 ? values + n : nullptr);
  d->atom_values_ramp = ramp;
  return 0;
}

extern "C" int32_t caps_molecule_info(caps_doc* d, int32_t atom, char* out, int32_t cap) {
  caps::Json r = caps::Json::object();
  try {
    const auto m = caps::molecule_info(d->frame, uint32_t(std::max(0, atom)));
    auto num = [](double x) { return std::isfinite(x) ? caps::Json(x) : caps::Json(); };
    r["ok"] = true;
    r["molecule"] = double(m.molecule + 1);
    r["atoms"] = double(m.atoms), r["bonds"] = double(m.bonds), r["rings"] = double(m.rings);
    r["formula"] = m.formula, r["smiles"] = m.smiles;
    r["mass"] = m.mass, r["monoisotopic"] = num(m.monoisotopic), r["dbe"] = m.dbe;
    caps::Json in = caps::Json::array();
    for (double x : m.inertia) in.push_back(caps::Json(x));
    r["inertia"] = std::move(in);
    r["inertia_defect"] = m.inertia_defect, r["rg"] = m.rg;
    r["has_charges"] = m.has_charges, r["net_charge"] = m.net_charge, r["dipole"] = num(m.dipole);
    r["rotatable"] = double(m.rotatable);
  } catch (const std::exception& e) {
    r = caps::Json::object();
    r["ok"] = false;
    r["error"] = std::string(e.what());
  }
  return report_out(r.dump(0), out, cap);
}

extern "C" int32_t caps_sasa(caps_doc* d, const char* json, char* out, int32_t cap) {
  caps::Json r = caps::Json::object();
  try {
    const caps::Json j = caps::Json::parse(json && *json ? json : "{}");
    const double probe = j.num("probe", 1.4);
    const int points = int(j.num("points", 200));
    const caps::System& s = d->frame;
    const auto res = caps::sasa(s, probe, points);
    r["ok"] = true;
    r["total"] = res.total;
    caps::Json area = caps::Json::array();
    for (double a : res.area) area.push_back(caps::Json(a));
    r["area"] = std::move(area);
    // groups: aromatic rings with their hydrogens, polar atoms (N O S) with theirs, everything else
    const auto p = caps::perceive(s);
    double ring = 0, polar = 0, rest = 0;
    for (size_t i = 0; i < s.atoms.size(); ++i) {
      uint32_t host = uint32_t(i);
      if (s.atoms[i].element == 1 && !p.nb[i].empty()) host = p.nb[i][0];
      const int z = s.atoms[host].element;
      if (p.aromatic[host]) ring += res.area[i];
      else if (z == 7 || z == 8 || z == 16) polar += res.area[i];
      else rest += res.area[i];
    }
    caps::Json groups = caps::Json::array();
    auto group = [&](const char* name, double a) {
      if (a <= 0) return;
      caps::Json g = caps::Json::object();
      g["name"] = std::string(name), g["area"] = a, g["share"] = res.total > 0 ? a / res.total : 0.0;
      groups.push_back(std::move(g));
    };
    group("Aromatic rings (+ their H)", ring);
    group("N, O and S (+ their H)", polar);
    group(ring > 0 || polar > 0 ? "Everything else (+ its H)" : "All atoms", rest);
    r["groups"] = std::move(groups);
    if (j.has("convergence") && j["convergence"].kind() == caps::Json::Bool && j["convergence"].boolean()) {
      caps::Json conv = caps::Json::array();
      const double ref = caps::sasa(s, probe, 400).total;
      for (int pts : {30, 60, 120, 200, 400}) {
        const double t = pts == 400 ? ref : caps::sasa(s, probe, pts).total;
        caps::Json c = caps::Json::object();
        c["points"] = double(pts), c["total"] = t, c["delta"] = ref > 0 ? (t - ref) / ref : 0.0;
        conv.push_back(std::move(c));
      }
      r["convergence"] = std::move(conv);
    }
    if (j.has("colour") && j["colour"].kind() == caps::Json::Bool && j["colour"].boolean()) {   // exposure: the area over the atom's full sphere
      d->atom_values.assign(s.atoms.size(), 0.0);
      for (size_t i = 0; i < s.atoms.size(); ++i) {
        const double R = caps::element(s.atoms[i].element).vdw + probe;
        d->atom_values[i] = res.area[i] / (4 * M_PI * R * R);
      }
      d->atom_values_ramp = 0;
    }
  } catch (const std::exception& e) {
    r = caps::Json::object();
    r["ok"] = false;
    r["error"] = std::string(e.what());
  }
  return report_out(r.dump(0), out, cap);
}

extern "C" int32_t caps_charges(caps_doc* d, const char* json, char* out, int32_t cap) {
  caps::Json r = caps::Json::object();
  try {
    const caps::Json j = caps::Json::parse(json && *json ? json : "{}");
    const std::string method = j.text("method", "gasteiger");
    const caps::System& s = d->frame;
    caps::ChargeReport rep;
    if (method == "forcefield") {
      if (!d->field || !d->field->ff) throw std::runtime_error("no force field assigned: assign one in Field first (its library charges, or the charges it computed)");
      rep = caps::describe_charges(s, d->field->ff->charge, "forcefield");
      rep.notes.push_back("from the Field assignment (" + d->field->base.name + ")");
    } else {
      rep = caps::compute_charges(s, method, j.text("path"));
    }
    if (j.has("apply") && j["apply"].kind() == caps::Json::Bool && j["apply"].boolean()) {
      push_undo(d, "Charges · " + method);
      for (size_t i = 0; i < d->traj.topology.atoms.size() && i < rep.q.size(); ++i) d->traj.topology.atoms[i].charge = rep.q[i];
      d->traj.topology.has_charges = true;
      d->field.reset();
      refresh(d);
    }
    auto arr = [](const std::vector<double>& v) { caps::Json a = caps::Json::array(); for (double x : v) a.push_back(caps::Json(x)); return a; };
    r["ok"] = true;
    r["method"] = rep.method;
    r["q"] = arr(rep.q);
    r["net"] = rep.net;
    r["max_abs"] = rep.max_abs;
    r["formal"] = double(rep.formal);
    caps::Json g = caps::Json::array();
    for (const auto& x : rep.groups) {
      caps::Json e = caps::Json::object();
      e["name"] = x.name, e["n"] = double(x.n), e["mean"] = x.mean, e["lo"] = x.lo, e["hi"] = x.hi;
      g.push_back(std::move(e));
    }
    r["groups"] = std::move(g);
    r["edges"] = arr(rep.edges);
    r["counts"] = arr(rep.counts);
    caps::Json notes = caps::Json::array();
    for (const auto& n : rep.notes) notes.push_back(caps::Json(n));
    r["notes"] = std::move(notes);
  } catch (const std::exception& e) {
    r = caps::Json::object();
    r["ok"] = false;
    r["error"] = std::string(e.what());
  }
  return report_out(r.dump(0), out, cap);
}

extern "C" uint32_t caps_category_colour(int32_t k) { return caps::category_colour(k); }

extern "C" int32_t caps_camera_focus(caps_doc* d, const caps_camera* cam, const int32_t* idx, int32_t n, double fill, caps_camera* out) {
  return guard([&] {
    if (!cam || !out) throw std::invalid_argument("camera_focus: no camera");
    const caps::System& s = shown(d);
    *out = *cam;
    out->pan_x = out->pan_y = 0;
    if (s.atoms.empty() || n <= 0) { out->zoom = 1; return 0; }
    // the view's rotation about the frame centre, as the renderer fits it (cell centre, else the atoms' box)
    caps::Vec3 centre{0, 0, 0};
    if (s.cell.valid()) centre = s.cell.origin + (s.cell.a + s.cell.b + s.cell.c) * 0.5;
    else {
      caps::Vec3 lo{1e300, 1e300, 1e300}, hi{-1e300, -1e300, -1e300};
      for (const auto& a : s.atoms) for (int k = 0; k < 3; ++k) { lo[k] = std::min(lo[k], a.pos[k]); hi[k] = std::max(hi[k], a.pos[k]); }
      centre = (lo + hi) * 0.5;
    }
    const double cy = std::cos(cam->yaw), sy = std::sin(cam->yaw), cp = std::cos(cam->pitch), sp = std::sin(cam->pitch);
    auto rot = [&](const caps::Vec3& p) {
      const caps::Vec3 q = p - centre;
      const double x = q[0] * cy + q[2] * sy, z = -q[0] * sy + q[2] * cy, y = q[1];
      return std::array<double, 2>{x, y * cp - z * sp};
    };
    double ax = 1e-6, ay = 1e-6;   // the whole frame's half extent (zoom 1)
    if (s.cell.valid())
      for (int i = 0; i < 2; ++i) for (int j = 0; j < 2; ++j) for (int k = 0; k < 2; ++k) {
        const auto r = rot(s.cell.origin + s.cell.a * i + s.cell.b * j + s.cell.c * k);
        ax = std::max(ax, std::fabs(r[0])); ay = std::max(ay, std::fabs(r[1]));
      }
    for (const auto& a : s.atoms) { const auto r = rot(a.pos); ax = std::max(ax, std::fabs(r[0])); ay = std::max(ay, std::fabs(r[1])); }
    double lx = 1e300, hx = -1e300, ly = 1e300, hy = -1e300;
    int used = 0;
    for (int32_t k = 0; k < n; ++k) {
      if (idx[k] < 0 || size_t(idx[k]) >= s.atoms.size()) continue;
      const auto r = rot(s.atoms[size_t(idx[k])].pos);
      lx = std::min(lx, r[0]); hx = std::max(hx, r[0]); ly = std::min(ly, r[1]); hy = std::max(hy, r[1]);
      ++used;
    }
    if (!used) { out->zoom = 1; return 0; }
    out->pan_x = -(lx + hx) / 2;
    out->pan_y = -(ly + hy) / 2;
    const double f = fill > 0 && fill <= 1 ? fill : 0.6;
    const double hx2 = std::max(2.0, (hx - lx) / 2 + 2), hy2 = std::max(2.0, (hy - ly) / 2 + 2);   // 2 Å for the atoms' own size
    out->zoom = std::clamp(f * std::min(ax / hx2, ay / hy2), 1.0, 40.0);
    return 0;
  });
}

extern "C" int32_t caps_set_vision(caps_doc* d, int32_t vision, double severity) {
  if (!d || vision < 0 || vision > 3) return -1;
  d->vision = vision;
  d->vision_severity = severity > 0 ? std::min(severity, 1.0) : 1.0;
  return 0;
}

// ---------------------------------------------------------------- recipes (design/boards/CommandLine, JupyterNotebook)
extern "C" int32_t caps_yaml_to_json(const char* yaml, char* out, int32_t cap) {
  try {
    return report_out(caps::yaml_parse(yaml ? yaml : "").dump(0), out, cap);
  } catch (const std::exception& e) {
    caps::Json r = caps::Json::object();
    r["ok"] = false;
    r["error"] = std::string(e.what());
    return report_out(r.dump(0), out, cap);
  }
}

extern "C" int32_t caps_recipe_check(const char* recipe, char* out, int32_t cap) {
  caps::Json r = caps::Json::object();
  const std::string text = recipe ? recipe : "";
  r["sha256"] = caps::sha256_hex(text);
  try {
    const auto t0 = text.find_first_not_of(" \t\r\n");
    const caps::Json j = t0 != std::string::npos && text[t0] == '{' ? caps::Json::parse(text) : caps::yaml_parse(text);
    const auto c = caps::check_recipe(j);
    r["ok"] = c.code == 0;
    r["code"] = double(c.code);
    r["error"] = c.error;
    r["name"] = j.is_object() ? j.text("name", "recipe") : std::string("recipe");
    caps::Json st = caps::Json::array();
    for (const auto& x : c.stages) {
      caps::Json e = caps::Json::object();
      e["name"] = x.name, e["summary"] = x.summary, e["ok"] = x.ok;
      st.push_back(std::move(e));
    }
    r["stages"] = std::move(st);
    r["protocol"] = c.protocol;
    caps::Json sch = caps::Json::array();
    for (const auto& g : c.schedule) {
      caps::Json e = caps::Json::object();
      e["label"] = g.label;
      e["ensemble"] = std::string(g.ensemble == caps::Ensemble::NPT ? "NPT" : g.ensemble == caps::Ensemble::NVE ? "NVE" : "NVT");
      e["ps"] = g.ps, e["t_start"] = g.t_start, e["t_end"] = g.t_end >= 0 ? g.t_end : g.t_start;
      e["pressure_bar"] = g.ensemble == caps::Ensemble::NPT ? g.pressure * 1.01325 : 0.0;
      sch.push_back(std::move(e));
    }
    r["schedule"] = std::move(sch);
  } catch (const std::exception& e) {
    r["ok"] = false;
    r["code"] = 2.0;
    r["error"] = std::string(e.what());
    r["stages"] = caps::Json::array();
  }
  return report_out(r.dump(0), out, cap);
}

extern "C" caps_doc* caps_recipe_run(const char* recipe, const char* options_json, caps_recipe_progress_fn progress, void* user, char* report, int32_t cap) {
  caps::Json rep = caps::Json::object();
  try {
    const std::string text = recipe ? recipe : "";
    const auto t0 = text.find_first_not_of(" \t\r\n");
    caps::Json r;
    try {
      r = t0 != std::string::npos && text[t0] == '{' ? caps::Json::parse(text) : caps::yaml_parse(text);
    } catch (const std::exception& e) {
      throw caps::RecipeError(2, e.what());
    }
    const caps::Json o = options_json && *options_json ? caps::Json::parse(options_json) : caps::Json::object();
    caps::RecipeOptions ro;
    ro.base_dir = o.text("base_dir", ".");
    ro.out_dir = o.text("out_dir", ".");
    ro.forcefield_dir = o.text("forcefield_dir", "");
    ro.seed = (long long)o.num("seed", -1);
    ro.threads = int(o.num("threads", 0));
    ro.sha256 = caps::sha256_hex(text);
    if (progress)
      ro.progress = [&](const caps::RecipeEvent& e) {
        if (!progress(e.stage, e.stages, e.name.c_str(), e.status.c_str(), e.detail.c_str(), e.fraction, user)) throw caps::RecipeError(4, "cancelled");
      };
    auto res = caps::run_recipe(r, ro);
    caps_doc* d = doc_of(res.system);
    d->prov = res.manifest;
    if (res.field) {
      d->field = std::make_unique<FieldState>();
      d->field->ff_path = res.forcefield;
      d->field->base.name = res.forcefield;
      d->field->ff = res.field;
      d->field->complete = true;
    }
    rep["exit"] = 0;
    caps::Json files = caps::Json::array();
    for (const auto& f : res.files) files.push_back(caps::Json(f));
    rep["files"] = std::move(files);
    rep["properties"] = caps::Json::parse(caps::properties_json(res.properties));
    rep["forcefield"] = res.forcefield;
    report_out(rep.dump(0), report, cap);
    return d;
  } catch (const caps::RecipeError& e) {
    g_error = e.what();
    rep["exit"] = e.code;
    rep["error"] = std::string(e.what());
  } catch (const std::exception& e) {
    g_error = e.what();
    rep["exit"] = 4;
    rep["error"] = std::string(e.what());
  }
  report_out(rep.dump(0), report, cap);
  return nullptr;
}

extern "C" int32_t caps_scene_json(caps_doc* d, const char* options_json, char* json, int32_t cap) {
  return guard([&] {
    const caps::Json o = options_json && *options_json ? caps::Json::parse(options_json) : caps::Json::object();
    const size_t max_atoms = size_t(std::max(1.0, o.num("max_atoms", 60000)));
    const bool hyd = !o.has("hydrogens") || o["hydrogens"].kind() != caps::Json::Bool || o["hydrogens"].boolean();
    const caps::System& s = d->frame;
    const size_t n = s.atoms.size();
    size_t heavy = 0;
    for (const auto& a : s.atoms) heavy += hyd || a.element != 1;
    int64_t maxmol = 0;
    for (const auto& a : s.atoms) maxmol = std::max(maxmol, a.mol);
    const int64_t every = heavy > max_atoms ? int64_t(std::ceil(double(heavy) / double(max_atoms))) : 1;
    std::vector<int32_t> map(n, -1);
    std::string z = "[", xyz = "[", bonds = "[";
    char b[96];
    int32_t k = 0;
    std::set<int> elems;
    for (size_t i = 0; i < n; ++i) {
      const auto& a = s.atoms[i];
      if (!hyd && a.element == 1) continue;
      if (every > 1 && (maxmol > 0 ? (a.mol % every) != 0 : (int64_t(i) % every) != 0)) continue;
      map[i] = k++;
      elems.insert(a.element);
      z += (k > 1 ? "," : "") + std::to_string(a.element);
      std::snprintf(b, sizeof b, "%s%.3f,%.3f,%.3f", k > 1 ? "," : "", a.pos[0], a.pos[1], a.pos[2]);
      xyz += b;
    }
    bool first = true;
    for (const auto& bd : s.bonds) {
      if (bd.i >= n || bd.j >= n || map[bd.i] < 0 || map[bd.j] < 0) continue;
      std::snprintf(b, sizeof b, "%s%d,%d", first ? "" : ",", map[bd.i], map[bd.j]);
      bonds += b;
      first = false;
    }
    std::string colours = "{", radii = "{";
    for (int e : elems) {
      const auto& el = caps::element(e);
      std::snprintf(b, sizeof b, "%s\"%d\":\"#%06x\"", colours.size() > 1 ? "," : "", e, el.rgb);
      colours += b;
      std::snprintf(b, sizeof b, "%s\"%d\":%.2f", radii.size() > 1 ? "," : "", e, el.vdw);
      radii += b;
    }
    std::string cell = "null";
    if (s.cell.valid()) {
      cell = "[";
      const caps::Vec3* v[4] = {&s.cell.origin, &s.cell.a, &s.cell.b, &s.cell.c};
      for (int q = 0; q < 4; ++q) {
        std::snprintf(b, sizeof b, "%s%.3f,%.3f,%.3f", q ? "," : "", (*v[q])[0], (*v[q])[1], (*v[q])[2]);
        cell += b;
      }
      cell += "]";
    }
    const std::string out = "{\"atoms\":" + std::to_string(n) + ",\"shown\":" + std::to_string(k) + ",\"z\":" + z + "],\"xyz\":" + xyz + "],\"bonds\":" + bonds +
                            "],\"colours\":" + colours + "},\"radii\":" + radii + "},\"cell\":" + cell + "}";
    return report_out(out, json, cap);
  });
}

// ---------------------------------------------------------------- reaction template editor (design/boards/ReactionTemplate)
extern "C" int32_t caps_template_view(const char* text, char* json, int32_t cap) {
  caps::Json j = caps::Json::object();
  try {
    const auto ts = caps::parse_templates(text ? text : "");
    caps::Json list = caps::Json::array();
    for (const auto& t : ts) list.push_back(caps::Json::parse(caps::template_view(t)));
    j["ok"] = true;
    j["templates"] = std::move(list);
  } catch (const std::exception& e) {
    j = caps::Json::object();
    j["ok"] = false;
    j["error"] = std::string(e.what());
  }
  return report_out(j.dump(0), json, cap);
}

extern "C" int32_t caps_template_test(caps_doc* d, const char* text, char* json, int32_t cap) {
  caps::Json j = caps::Json::object();
  try {
    const auto ts = caps::parse_templates(text ? text : "");
    caps::System s = d->traj.frame(d->current);
    if (!s.unwrapped) caps::make_molecules_whole(s);
    caps::Json list = caps::Json::array();
    for (size_t k = 0; k < ts.size(); ++k) {
      caps::Json o = caps::Json::object();
      o["name"] = ts[k].name;
      o["sites"] = caps::count_sites(s, ts[k]);
      const auto m = caps::find_matches(s, ts[k], int(k));
      o["matches"] = double(m.size());
      o["closest"] = m.empty() ? -1.0 : m.front().distance;
      list.push_back(std::move(o));
    }
    j["ok"] = true;
    j["templates"] = std::move(list);
  } catch (const std::exception& e) {
    j = caps::Json::object();
    j["ok"] = false;
    j["error"] = std::string(e.what());
  }
  return report_out(j.dump(0), json, cap);
}

// ---------------------------------------------------------------- large systems (design/boards/MillionAtoms)
extern "C" int32_t caps_render_stats(caps_doc* d, int64_t* near, int64_t* mid, int64_t* far, int64_t* bonds) {
  const auto& st = d->renderer.stats;
  if (near) *near = int64_t(st.near);
  if (mid) *mid = int64_t(st.mid);
  if (far) *far = int64_t(st.far);
  if (bonds) *bonds = int64_t(st.bonds);
  return 0;
}

extern "C" int32_t caps_memory(caps_doc* d, char* json, int32_t cap) {
  const auto& t = d->traj;
  size_t topo = sizeof(caps::System) + t.topology.atoms.capacity() * sizeof(caps::Atom) + t.topology.bonds.capacity() * sizeof(caps::Bond);
  for (const auto& a : t.topology.atoms) topo += a.name.capacity() > 15 ? a.name.capacity() : 0;   // short strings live inside the atom
  size_t frames = 0;
  for (const auto& p : t.positions) frames += p.capacity() * sizeof(caps::Vec3);
  frames += t.cells.capacity() * sizeof(caps::Cell);
  const size_t shown = d->frame.atoms.capacity() * sizeof(caps::Atom);   // the current frame as a structure
  caps::Json j = caps::Json::object();
  j["atoms"] = double(t.topology.atoms.size());
  j["frames"] = double(t.frames());
  j["topology_bytes"] = double(topo + shown);
  j["frame_bytes"] = double(frames);
  j["per_atom_bytes"] = t.topology.atoms.empty() ? 0.0 : double(topo + shown + frames) / double(t.topology.atoms.size());
  j["atom_struct_bytes"] = double(sizeof(caps::Atom));
  return report_out(j.dump(0), json, cap);
}

// ---------------------------------------------------------------- v20 polymer statistics (design/boards row 18)

namespace {
caps::Json num_array(const std::vector<double>& v) {
  caps::Json a = caps::Json::array();
  for (double x : v) a.push_back(caps::Json(x));
  return a;
}
int32_t json_error(const std::exception& e, char* out, int32_t cap) {
  caps::Json r = caps::Json::object();
  r["ok"] = false;
  r["error"] = std::string(e.what());
  return report_out(r.dump(0), out, cap);
}
}  // namespace

extern "C" int32_t caps_chain_lengths(const char* json, char* out, int32_t cap) {
  try {
    const caps::Json j = caps::Json::parse(json && *json ? json : "{}");
    const std::string dist = j.text("distribution", "schulz-zimm");
    const double nn = std::max(2.0, j.num("nn", 40)), pdi = std::max(1.0, j.num("pdi", 1.1)), m0 = j.num("m0", 104.15);
    const int count = std::clamp(int(j.num("count", 20)), 1, 100000), best_of = std::clamp(int(j.num("best_of", 1)), 1, 1000);
    const uint64_t seed = uint64_t(std::max(0.0, j.num("seed", 2026)));
    auto stats = [&](const std::vector<int>& L, double& snn, double& swn) {
      double s1 = 0, s2 = 0;
      for (int n : L) s1 += n, s2 += double(n) * n;
      snn = s1 / double(L.size()), swn = s1 > 0 ? s2 / s1 : 0;
    };
    // best_of > 1: that many independent draws, the one whose Đ is closest to the target kept (said so in the report)
    std::vector<int> L;
    double best = 1e300;
    int kept = 0;
    for (int t = 0; t < best_of; ++t) {
      auto c = caps::draw_chain_lengths(dist, nn, pdi, count, seed + uint64_t(t) * 7919);
      double a, b;
      stats(c, a, b);
      const double target = dist == "flory" ? 2 - 1 / nn : dist == "poisson" ? 1 + (nn - 1) / (nn * nn) : dist == "monodisperse" ? 1.0 : pdi;
      const double miss = std::fabs(b / a - target) + 0.1 * std::fabs(a - nn) / nn;
      if (miss < best) best = miss, L = std::move(c), kept = t;
    }
    double snn, swn;
    stats(L, snn, swn);
    caps::Json r = caps::Json::object();
    r["ok"] = true;
    caps::Json lengths = caps::Json::array();
    for (int n : L) lengths.push_back(caps::Json(double(n)));
    r["lengths"] = std::move(lengths);
    auto sorted = L;
    std::sort(sorted.begin(), sorted.end());
    caps::Json s = caps::Json::object();
    s["nn"] = snn, s["mn"] = snn * m0, s["mw"] = swn * m0, s["pdi"] = snn > 0 ? swn / snn : 1.0;
    s["min"] = double(sorted.front()), s["max"] = double(sorted.back());
    double sum = 0;
    for (int n : L) sum += n;
    s["sum"] = sum;
    r["sample"] = std::move(s);
    // the distribution's own averages (discrete ones computed, the Gamma's are the inputs)
    const double tpdi = dist == "schulz-zimm" ? pdi : dist == "flory" ? 2 - 1 / nn : dist == "poisson" ? 1 + (nn - 1) / (nn * nn) : 1.0;
    caps::Json t = caps::Json::object();
    t["nn"] = nn, t["mn"] = nn * m0, t["mw"] = nn * m0 * tpdi, t["pdi"] = tpdi;
    r["target"] = std::move(t);
    r["k"] = dist == "schulz-zimm" ? (pdi > 1.0001 ? 1 / (pdi - 1) : 1e4) : 0.0;
    r["kept_draw"] = double(kept);
    r["draws"] = double(best_of);
    std::vector<double> xs, nf, wf;
    const double hi = std::max(double(sorted.back()) * 1.15, nn * (1 + 4 * std::sqrt(std::max(0.0, tpdi - 1)) + 0.3));
    for (int i = 0; i <= 160; ++i) {
      const double x = std::max(1.0, hi * i / 160);
      const double p = caps::chain_length_pdf(dist, nn, pdi, x);
      xs.push_back(x), nf.push_back(p), wf.push_back(x * p / nn);
    }
    caps::Json curve = caps::Json::object();
    curve["n"] = num_array(xs), curve["number"] = num_array(nf), curve["weight"] = num_array(wf);
    r["curve"] = std::move(curve);
    return report_out(r.dump(0), out, cap);
  } catch (const std::exception& e) {
    return json_error(e, out, cap);
  }
}

extern "C" int32_t caps_copolymer(const char* json, char* out, int32_t cap) {
  try {
    const caps::Json j = caps::Json::parse(json && *json ? json : "{}");
    const auto m = caps::copolymer_terminal(j.num("r1", 1), j.num("r2", 1), j.num("f1", 0.5));
    const int dp = std::clamp(int(j.num("dp", 80)), 2, 100000);
    const uint64_t seed = uint64_t(std::max(0.0, j.num("seed", 1)));
    caps::Json r = caps::Json::object();
    r["ok"] = true;
    r["F1"] = m.F1, r["paa"] = m.paa, r["pbb"] = m.pbb, r["run_a"] = m.run_a, r["run_b"] = m.run_b;
    r["azeotrope"] = m.azeotrope >= 0 ? caps::Json(m.azeotrope) : caps::Json();
    std::vector<double> f, F;
    for (int i = 0; i <= 100; ++i) f.push_back(i / 100.0), F.push_back(caps::mayo_lewis(m.r1, m.r2, i / 100.0));
    caps::Json curve = caps::Json::object();
    curve["f1"] = num_array(f), curve["F1"] = num_array(F);
    r["curve"] = std::move(curve);
    // the sequence Grow draws for chain 0 at this seed (caps_grow_chains uses seed + 101·chain)
    caps::ChainSpec spec;
    spec.units.resize(2);
    spec.sequence = caps::Sequence::Terminal;
    spec.dp = dp, spec.r1 = m.r1, spec.r2 = m.r2, spec.weights = {m.f1, 1 - m.f1};
    const auto seq = caps::chain_sequence(spec, seed);
    caps::Json sa = caps::Json::array();
    int na = 0, runs_a = 0, runs_b = 0, longest = 0, run = 0;
    for (size_t i = 0; i < seq.size(); ++i) {
      sa.push_back(caps::Json(double(seq[i])));
      na += seq[i] == 0;
      if (i == 0 || seq[i] != seq[i - 1]) (seq[i] == 0 ? runs_a : runs_b)++, run = 1;
      else ++run;
      longest = std::max(longest, run);
    }
    r["sequence"] = std::move(sa);
    caps::Json c = caps::Json::object();
    c["F1"] = double(na) / double(seq.size());
    c["run_a"] = runs_a ? double(na) / runs_a : 0.0, c["run_b"] = runs_b ? double(seq.size() - size_t(na)) / runs_b : 0.0;
    c["longest"] = double(longest), c["a"] = double(na), c["b"] = double(seq.size() - size_t(na));
    r["chain"] = std::move(c);
    return report_out(r.dump(0), out, cap);
  } catch (const std::exception& e) {
    return json_error(e, out, cap);
  }
}

extern "C" int32_t caps_stereo(const char* json, char* out, int32_t cap) {
  try {
    const caps::Json j = caps::Json::parse(json && *json ? json : "{}");
    auto model_json = [](const caps::StereoModel& m) {
      caps::Json o = caps::Json::object();
      o["kind"] = m.kind, o["pm"] = m.pm, o["p_mr"] = m.p_mr, o["p_rm"] = m.p_rm, o["mm"] = m.mm, o["mr"] = m.mr, o["rr"] = m.rr;
      o["pentads"] = num_array(std::vector<double>(m.pentads.begin(), m.pentads.end()));
      return o;
    };
    const caps::StereoModel m = j.text("model", "bernoulli") == "markov" ? caps::stereo_markov(j.num("p_mr", 0.5), j.num("p_rm", 0.5))
                                                                          : caps::stereo_bernoulli(j.num("pm", 0.5));
    caps::Json r = caps::Json::object();
    r["ok"] = true;
    r["model"] = model_json(m);
    caps::Json names = caps::Json::array();
    for (const char* n : caps::pentad_names()) names.push_back(caps::Json(std::string(n)));
    r["names"] = std::move(names);
    std::string dyads = j.text("dyads");
    if (dyads.empty() && j.has("dp")) dyads = caps::draw_dyads(m, int(j.num("dp", 200)), uint64_t(std::max(0.0, j.num("seed", 1))));
    if (!dyads.empty()) {
      const auto c = caps::count_stereo(dyads);
      caps::Json o = caps::Json::object();
      o["dyads"] = dyads, o["m"] = double(c.m), o["r"] = double(c.r);
      const double nt = std::max(1, c.mm + c.mr + c.rr), np = std::max(1, c.pentad_total);
      o["mm"] = c.mm / nt, o["mr"] = c.mr / nt, o["rr"] = c.rr / nt, o["triads"] = double(c.mm + c.mr + c.rr);
      std::vector<double> p;
      for (int k : c.pentads) p.push_back(k / np);
      o["pentads"] = num_array(p);
      o["pentad_count"] = double(c.pentad_total);
      r["chain"] = std::move(o);
    }
    if (j.has("measured")) {
      std::array<double, caps::kPentadCount> y{};
      size_t k = 0;
      for (const auto& x : j["measured"].items())
        if (k < y.size()) y[k++] = x.number();
      double rb = 0, rm = 0;
      auto fb = caps::fit_bernoulli(y, &rb), fm = caps::fit_markov(y, &rm);
      caps::Json o = caps::Json::object();
      o["bernoulli"] = model_json(fb), o["bernoulli_rms"] = rb;
      o["markov"] = model_json(fm), o["markov_rms"] = rm;
      // Bernoulli predicts mm·rr = (mr/2)²; the triads from the measured pentads (mm = mmmm + mmmr + rmmr, …)
      const double s = std::accumulate(y.begin(), y.end(), 0.0);
      const double mm = (y[0] + y[1] + y[2]) / s, rr = (y[7] + y[8] + y[9]) / s, mr = (y[3] + y[4] + y[5] + y[6]) / s;
      o["mm_rr"] = mm * rr, o["mr2_4"] = mr * mr / 4;
      r["fit"] = std::move(o);
    }
    return report_out(r.dump(0), out, cap);
  } catch (const std::exception& e) {
    return json_error(e, out, cap);
  }
}

extern "C" int32_t caps_ris_cn(double temperature, int32_t nmax, double* out) {
  try {
    const auto c = caps::ris_cn(caps::RisModel{}, temperature, std::max(1, int(nmax)));
    for (size_t k = 0; k < c.size(); ++k) out[k] = c[k];
    return int32_t(c.size());
  } catch (const std::exception& e) {
    g_error = e.what();
    return -1;
  }
}

extern "C" int32_t caps_chi_md(const char* json, caps_stage_fn progress, void* user, char* out, int32_t cap) {
  caps::Json r = caps::Json::object();
  try {
    const caps::Json j = caps::Json::parse(json && *json ? json : "{}");
    caps::ChiMdOptions o;
    o.polymer = spec_from(j["polymer"].dump());
    o.chains = int(j.num("chains", 6));
    if (j.has("polymer_b")) o.b_polymer = true, o.polymer_b = spec_from(j["polymer_b"].dump()), o.chains_b = int(j.num("chains_b", 6));
    o.solvent_smiles = j.text("solvent");
    o.solvent_molecules = int(j.num("solvent_molecules", 0));
    o.temperature = j.num("temperature", 300);
    o.pressure = j.num("pressure", 1);
    o.eq_ps = j.num("eq_ps", 20), o.prod_ps = j.num("prod_ps", 20);
    o.seed = uint64_t(j.num("seed", 1));
    if (progress) o.progress = [&](const std::string& st, double f) { return progress(st.c_str(), f, user) == 0; };
    const caps::ChiMdResult c = caps::chi_by_md(o);
    r["ok"] = true;
    r["chi"] = c.chi, r["chi_error"] = c.chi_error, r["phi_a"] = c.phi_a, r["v_ref"] = c.v_ref, r["de_mix"] = c.de_mix;
    caps::Json cells = caps::Json::array();
    for (const caps::ChiMdCell* x : {&c.a, &c.b, &c.mix}) {
      caps::Json e = caps::Json::object();
      e["name"] = x->name, e["atoms"] = double(x->atoms), e["molecules"] = double(x->molecules), e["density"] = x->density;
      e["ced"] = x->e_density / 1.4393e-4, e["ced_error"] = x->e_error / 1.4393e-4;
      cells.push_back(std::move(e));
    }
    r["cells"] = std::move(cells);
    caps::Json n = caps::Json::array();
    for (const auto& x : c.notes) n.push_back(x);
    r["notes"] = std::move(n);
  } catch (const std::exception& e) {
    r["ok"] = false;
    r["error"] = std::string(e.what());
  }
  return report_out(r.dump(), out, cap);
}

extern "C" int32_t caps_chi_contacts(const char* json, caps_stage_fn progress, void* user, char* out, int32_t cap) {
  caps::Json r = caps::Json::object();
  try {
    const caps::Json j = caps::Json::parse(json && *json ? json : "{}");
    caps::ChiPairOptions o;
    o.a_smiles = j.text("a"), o.b_smiles = j.text("b");
    o.forcefield = j.text("forcefield");
    o.samples = int(j.num("samples", o.samples));
    o.pack_trials = int(j.num("pack_trials", o.pack_trials));
    if (j.has("temperatures") && j["temperatures"].is_array()) {
      o.temperatures.clear();
      for (const auto& x : j["temperatures"].items()) o.temperatures.push_back(x.number());
    }
    o.report_temperature = j.num("t", o.report_temperature);
    o.seed = uint64_t(j.num("seed", 1));
    o.threads = int(j.num("threads", 0));
    if (progress) o.progress = [&](const std::string& st, double f) { return progress(st.c_str(), f, user) == 0; };
    const caps::ChiPairResult c = caps::chi_by_contacts(o);
    r["ok"] = true;
    r["chi"] = c.chi_report, r["chi_error"] = c.chi_report_error, r["fit_a"] = c.fit_a, r["fit_b"] = c.fit_b;
    r["temperatures"] = num_array(c.temperatures), r["chi_t"] = num_array(c.chi), r["chi_t_error"] = num_array(c.chi_error);
    caps::Json kinds = caps::Json::array();
    for (const caps::ChiPairKind* k : {&c.aa, &c.ab, &c.ba, &c.bb}) {
      caps::Json e = caps::Json::object();
      e["name"] = k->name, e["z"] = k->z, e["z_error"] = k->z_error, e["e_min"] = k->e_min, e["e_mean"] = k->e_mean;
      e["e_t"] = num_array(k->e_t), e["hist_e"] = num_array(k->hist_e), e["hist_p"] = num_array(k->hist_p);
      kinds.push_back(std::move(e));
    }
    r["kinds"] = std::move(kinds);
    r["forcefield"] = c.forcefield, r["a"] = c.a_name, r["b"] = c.b_name;
    caps::Json n = caps::Json::array();
    for (const auto& x : c.notes) n.push_back(x);
    r["notes"] = std::move(n);
  } catch (const std::exception& e) {
    r["ok"] = false;
    r["error"] = std::string(e.what());
  }
  return report_out(r.dump(), out, cap);
}

extern "C" int32_t caps_blend_phase(const char* json, char* out, int32_t cap) {
  try {
    const caps::Json j = caps::Json::parse(json && *json ? json : "{}");
    const double na = std::max(1.0, j.num("na", 100)), nb = std::max(1.0, j.num("nb", 200)), A = j.num("a", -0.02), B = j.num("b", 15),
                 T = j.num("t", 300);
    const auto crit = caps::blend_critical(na, nb);
    auto chi_at = [&](double t) { return A + B / t; };
    caps::Json r = caps::Json::object();
    r["ok"] = true;
    r["chi_c"] = crit.chi_c, r["phi_c"] = crit.phi_c;
    const bool has_tc = B != 0 && (crit.chi_c - A) * B > 0;
    const double tc = has_tc ? B / (crit.chi_c - A) : 0;
    r["tc"] = has_tc ? caps::Json(tc) : caps::Json();
    r["kind"] = B > 0 ? "ucst" : B < 0 ? "lcst" : "none";
    r["chi_t"] = chi_at(T);
    double lo, hi;
    if (caps::blend_binodal(na, nb, chi_at(T), lo, hi)) r["coexist"] = num_array({lo, hi});
    if (caps::blend_spinodal(na, nb, chi_at(T), lo, hi)) r["spinodal"] = num_array({lo, hi});
    // curves in T (two phases below Tc for UCST, above for LCST): χ from just above χc to χ at the far temperature
    if (has_tc) {
      const double tfar = B > 0 ? std::max(1.0, j.num("t_min", std::max(1.0, tc * 0.55))) : j.num("t_max", tc * 1.6);
      const double chi_far = chi_at(tfar);
      std::vector<double> bp, bt, sp, st;
      std::vector<std::array<double, 3>> left, right, sleft, sright;
      double glo = 0, ghi = 0;
      for (int i = 1; i <= 240; ++i) {
        const double u = double(i) / 240, chi = crit.chi_c + (chi_far - crit.chi_c) * u * u;   // dense near the critical point
        const double t = B / (chi - A);
        double b1, b2, s1, s2;
        if (caps::blend_binodal(na, nb, chi, b1, b2, glo, ghi)) glo = b1, ghi = b2, left.push_back({b1, t, 0}), right.push_back({b2, t, 0});
        if (caps::blend_spinodal(na, nb, chi, s1, s2)) sleft.push_back({s1, t, 0}), sright.push_back({s2, t, 0});
      }
      auto join = [](const std::vector<std::array<double, 3>>& L, const std::vector<std::array<double, 3>>& R, double pc, double tc, std::vector<double>& x, std::vector<double>& y) {
        for (auto it = L.rbegin(); it != L.rend(); ++it) x.push_back((*it)[0]), y.push_back((*it)[1]);
        x.push_back(pc), y.push_back(tc);
        for (const auto& p : R) x.push_back(p[0]), y.push_back(p[1]);
      };
      join(left, right, crit.phi_c, tc, bp, bt);
      join(sleft, sright, crit.phi_c, tc, sp, st);
      caps::Json b = caps::Json::object();
      b["phi"] = num_array(bp), b["t"] = num_array(bt);
      r["binodal"] = std::move(b);
      caps::Json s = caps::Json::object();
      s["phi"] = num_array(sp), s["t"] = num_array(st);
      r["spinodal_curve"] = std::move(s);
    }
    return report_out(r.dump(0), out, cap);
  } catch (const std::exception& e) {
    return json_error(e, out, cap);
  }
}

extern "C" int32_t caps_solvent_chi(const char* json, char* out, int32_t cap) {
  try {
    const caps::Json j = caps::Json::parse(json && *json ? json : "{}");
    const double dp = j.num("delta_polymer", 18.6), t = j.num("t", 298.15);
    caps::Json r = caps::Json::object(), list = caps::Json::array();
    r["ok"] = true;
    r["rt"] = 8.314462618 * t;
    if (j.has("solvents"))
      for (const auto& s : j["solvents"].items()) {
        const double chi = caps::hildebrand_chi(s.num("v", 100), s.num("delta", 18), dp, t);
        caps::Json o = caps::Json::object();
        o["name"] = s.text("name"), o["chi"] = chi;
        o["predicted"] = chi < 0.45 ? "solvent" : chi <= 0.55 ? "borderline" : "non-solvent";
        list.push_back(std::move(o));
      }
    r["solvents"] = std::move(list);
    return report_out(r.dump(0), out, cap);
  } catch (const std::exception& e) {
    return json_error(e, out, cap);
  }
}

extern "C" int32_t caps_ewald_params(caps_doc* d, const char* json, char* out, int32_t cap) {
  try {
    const caps::Json j = caps::Json::parse(json && *json ? json : "{}");
    const double rc = std::max(1.0, j.num("cutoff", 12)), tol = std::clamp(j.num("tolerance", 1e-5), 1e-12, 0.1),
                 spacing = std::max(0.3, j.num("spacing", 1.2));
    const int order = std::clamp(int(j.num("order", 4)), 3, 12);
    caps::Json r = caps::Json::object();
    r["ok"] = true;
    const double beta = caps::ewald_beta(rc, tol);
    r["beta"] = beta, r["beta_rc"] = beta * rc;
    caps::Json table = caps::Json::array();
    for (double t : {1e-4, 1e-5, 1e-6}) {
      caps::Json o = caps::Json::object();
      const double b = caps::ewald_beta(rc, t);
      o["tolerance"] = t, o["beta"] = b, o["beta_rc"] = b * rc;
      table.push_back(std::move(o));
    }
    r["table"] = std::move(table);
    std::array<double, 3> L = {0, 0, 0};
    if (j.has("edges")) {
      size_t k = 0;
      for (const auto& x : j["edges"].items())
        if (k < 3) L[k++] = x.number();
    } else if (d && d->frame.cell.valid()) {
      L = {caps::norm(d->frame.cell.a), caps::norm(d->frame.cell.b), caps::norm(d->frame.cell.c)};
    }
    if (L[0] > 0 && L[1] > 0 && L[2] > 0) {
      caps::Json mesh = caps::Json::array(), sp = caps::Json::array(), edges = caps::Json::array();
      for (double e : L) {
        const int n = caps::pme_mesh_size(e, spacing, order);
        mesh.push_back(caps::Json(double(n))), sp.push_back(caps::Json(e / n)), edges.push_back(caps::Json(e));
      }
      r["edges"] = std::move(edges), r["mesh"] = std::move(mesh), r["spacing"] = std::move(sp);
      r["fits"] = 2 * rc <= std::min({L[0], L[1], L[2]});
    }
    std::vector<double> x, y;
    for (int i = 0; i <= 200; ++i) {
      const double rr = 1 + (rc * 1.35 - 1) * i / 200.0;
      x.push_back(rr), y.push_back(std::erfc(beta * rr));
    }
    caps::Json curve = caps::Json::object();
    curve["r"] = num_array(x), curve["erfc"] = num_array(y);
    r["curve"] = std::move(curve);
    return report_out(r.dump(0), out, cap);
  } catch (const std::exception& e) {
    return json_error(e, out, cap);
  }
}

extern "C" int32_t caps_atom_residues(caps_doc* d, int32_t* out, int32_t cap) {
  const auto& atoms = d->frame.atoms;
  for (size_t i = 0; i < atoms.size() && out && int32_t(i) < cap; ++i) out[i] = int32_t(atoms[i].resid);
  return int32_t(atoms.size());
}

// ---------------------------------------------------------------- v20 display, lens, resolution, hydrogens (row 19)

namespace {
caps_doc* doc_of_system(caps::System sys, const caps_doc* from, const std::string& engine, const std::string& summary, caps::KeyValues params) {
  auto* d = new caps_doc;
  d->traj.topology = sys;
  std::vector<caps::Vec3> p;
  for (const auto& a : sys.atoms) p.push_back(a.pos);
  d->traj.positions.push_back(std::move(p));
  d->traj.cells.push_back(sys.cell);
  d->traj.timesteps.push_back(0);
  if (from) d->prov = from->prov;   // the new structure's history starts from the one it came from
  refresh(d);
  prov_step(d, engine, summary, std::move(params));
  return d;
}
}  // namespace

extern "C" caps_doc* caps_doc_copy(caps_doc* src) {
  try {
    caps_doc* d = doc_of_system(src->frame, src, "doc.copy", "a copy of the structure", {{"atoms", std::to_string(src->frame.atoms.size())}});
    // the same atoms: the force-field assignment (types set by hand, entered parameters, charges) holds for the copy
    if (src->field) d->field = std::make_unique<FieldState>(*src->field);
    d->held_mol = src->held_mol;
    return d;
  } catch (const std::exception& e) {
    g_error = e.what();
    return nullptr;
  }
}

extern "C" int32_t caps_set_display(caps_doc* d, const char* json) {
  try {
    const caps::Json j = caps::Json::parse(json && *json ? json : "{}");
    auto flag = [&](const caps::Json& o, const char* k, bool def) { return o.has(k) && o[k].kind() == caps::Json::Bool ? o[k].boolean() : def; };
    auto& D = d->display;
    D.polar_h_only = flag(j, "polar_h_only", D.polar_h_only);
    D.selection_full = flag(j, "selection_full", D.selection_full);
    if (j.has("clip")) {
      const caps::Json& C = j["clip"];
      D.clip = flag(C, "on", D.clip);
      D.clip_axis = std::clamp(int(C.num("axis", D.clip_axis)), 0, 2);
      D.clip_from = std::clamp(C.num("from", D.clip_from), 0.0, 1.0);
      D.clip_to = std::clamp(C.num("to", D.clip_to), 0.0, 1.0);
      if (D.clip_to < D.clip_from) std::swap(D.clip_from, D.clip_to);
      D.clip_invert = flag(C, "invert", D.clip_invert);
    }
    if (j.has("lens")) {
      const caps::Json& L = j["lens"];
      D.lens = flag(L, "on", D.lens);
      if (L.has("centre")) D.lens_centre = int(L["centre"].number());
      D.lens_radius = std::clamp(L.num("radius", D.lens_radius), 1.0, 500.0);
      D.lens_inside = std::clamp(int(L.num("inside", D.lens_inside)), 0, 4);
      D.lens_outside = std::clamp(int(L.num("outside", D.lens_outside)), 0, 4);
      D.lens_dim = flag(L, "dim", D.lens_dim);
      const auto& A = d->frame.atoms;
      if (D.lens && (D.lens_centre < 0 || size_t(D.lens_centre) >= A.size()) && !A.empty()) {   // no centre yet: the atom nearest the middle
        caps::Vec3 c{0, 0, 0};
        if (d->frame.cell.valid()) c = d->frame.cell.origin + (d->frame.cell.a + d->frame.cell.b + d->frame.cell.c) * 0.5;
        else for (const auto& a : A) c = c + a.pos * (1.0 / double(A.size()));
        double best = 1e300;
        for (size_t i = 0; i < A.size(); ++i) {
          const caps::Vec3 dv = A[i].pos - c;
          if (const double r = caps::dot(dv, dv); r < best && A[i].element != 1) best = r, D.lens_centre = int(i);
        }
      }
    }
    return 0;
  } catch (const std::exception& e) {
    g_error = e.what();
    return -1;
  }
}

extern "C" int32_t caps_lens_inside(caps_doc* d, int32_t atom) {
  const auto& D = d->display;
  const auto& A = d->frame.atoms;
  if (!D.lens || D.lens_centre < 0 || size_t(D.lens_centre) >= A.size() || atom < 0 || size_t(atom) >= A.size()) return 1;
  caps::Vec3 dv = A[size_t(atom)].pos - A[size_t(D.lens_centre)].pos;
  if (d->frame.cell.valid()) dv = d->frame.cell.minimum_image(dv);
  return caps::dot(dv, dv) <= D.lens_radius * D.lens_radius ? 1 : 0;
}

extern "C" int32_t caps_display_counts(caps_doc* d, char* out, int32_t cap) {
  const auto& A = d->frame.atoms;
  caps::Json r = caps::Json::object();
  long h = 0, polar = 0;
  std::vector<std::vector<uint32_t>> nb(A.size());
  for (const auto& b : d->frame.bonds) nb[b.i].push_back(b.j), nb[b.j].push_back(b.i);
  for (size_t i = 0; i < A.size(); ++i)
    if (A[i].element == 1) {
      ++h;
      if (nb[i].empty() || A[nb[i][0]].element != 6) ++polar;
    }
  backbone_mask(d);
  r["atoms"] = double(A.size()), r["h"] = double(h), r["heavy"] = double(long(A.size()) - h), r["polar_h"] = double(polar);
  r["chains"] = double(d->bb_chains), r["backbone_atoms"] = double(d->bb_atoms_on);
  const auto& D = d->display;
  if (D.lens && D.lens_centre >= 0 && size_t(D.lens_centre) < A.size()) {
    long in = 0, in_h = 0;
    for (size_t i = 0; i < A.size(); ++i)
      if (caps_lens_inside(d, int32_t(i))) ++in, in_h += A[i].element == 1 ? 1 : 0;
    r["lens_atoms"] = double(in), r["lens_h"] = double(in_h), r["lens_centre"] = double(D.lens_centre);
  }
  return report_out(r.dump(0), out, cap);
}


extern "C" int32_t caps_hydrogen_plan(caps_doc* d, char* out, int32_t cap) {
  try {
    caps::System perceived;
    const bool geo = needs_geometry_orders(d->frame);
    if (geo) perceived = d->frame, geometry_orders(perceived);
    if (d->ph >= 0) {   // residues protonated at the pH first: their formal charges change what each atom lacks
      if (!geo) perceived = d->frame;
      caps::protonate_residues(perceived, d->ph);
    }
    const auto& s = geo || d->ph >= 0 ? perceived : d->frame;
    std::vector<char> sel;
    if (d->selection.size() == s.atoms.size() && std::any_of(d->selection.begin(), d->selection.end(), [](char c) { return c != 0; })) sel = d->selection;
    const auto rows = caps::hydrogen_plan(s, sel);
    caps::Json r = caps::Json::object(), list = caps::Json::array();
    int heavy = 0, h = 0, add = 0;
    for (const auto& a : s.atoms) (a.element == 1 ? h : heavy)++;
    for (const auto& row : rows) {
      caps::Json o = caps::Json::object();
      o["label"] = row.label, o["atoms"] = double(row.atoms), o["hydrogens"] = double(row.hydrogens);
      add += row.hydrogens;
      list.push_back(std::move(o));
    }
    double q = 0;
    for (const auto& a : s.atoms) q += a.charge;
    r["ok"] = true;
    r["rows"] = std::move(list);
    r["heavy"] = double(heavy), r["h"] = double(h), r["add"] = double(add), r["net_charge"] = q;
    int aromatic = 0;
    for (const auto& b : s.bonds) aromatic += b.order == 4 ? 1 : 0;
    r["aromatic_bonds"] = double(aromatic);
    r["orders_from_geometry"] = geo;
    r["selection"] = !sel.empty();
    return report_out(r.dump(0), out, cap);
  } catch (const std::exception& e) {
    return json_error(e, out, cap);
  }
}

extern "C" int32_t caps_resolution_summary(caps_doc* d, const char* json, char* out, int32_t cap) {
  try {
    const caps::Json j = caps::Json::parse(json && *json ? json : "{}");
    const int per = std::clamp(int(j.num("per_bead", 5)), 1, 100);
    const auto aa = caps::all_atom_summary(d->frame);
    caps::ResolutionReport ua, cg;
    caps::united_atom(d->frame, &ua);
    caps::coarse_grain(d->frame, per, &cg);
    auto part = [](const caps::ResolutionReport& x) {
      caps::Json o = caps::Json::object();
      o["sites"] = double(x.sites), o["hydrogens"] = double(x.hydrogens), o["mass"] = x.mass;
      return o;
    };
    caps::Json r = caps::Json::object();
    r["ok"] = true;
    r["all_atom"] = part(aa), r["united_atom"] = part(ua), r["coarse_grained"] = part(cg);
    r["per_bead"] = double(per);
    return report_out(r.dump(0), out, cap);
  } catch (const std::exception& e) {
    return json_error(e, out, cap);
  }
}

extern "C" caps_doc* caps_resolution_convert(caps_doc* d, const char* json, char* report, int32_t cap) {
  try {
    const caps::Json j = caps::Json::parse(json && *json ? json : "{}");
    const std::string to = j.text("to", "united-atom");
    caps::ResolutionReport rep;
    caps::System s;
    if (to == "united-atom") s = caps::united_atom(d->frame, &rep);
    else if (to == "coarse-grained") s = caps::coarse_grain(d->frame, std::clamp(int(j.num("per_bead", 5)), 1, 100), &rep);
    else throw std::runtime_error("resolution: united-atom or coarse-grained");
    auto* nd = doc_of_system(std::move(s), d, "model.resolution", "converted to " + to + " (mass conserved: " + g6(rep.mass) + " g/mol)",
                             {{"to", to}, {"sites", std::to_string(rep.sites)}});
    std::string t;
    for (const auto& n : rep.notes) t += n + "\n";
    report_out(t, report, cap);
    return nd;
  } catch (const std::exception& e) {
    g_error = e.what();
    return nullptr;
  }
}

extern "C" caps_doc* caps_backmap(caps_doc* d, const char* beads_path, int32_t per_bead, int32_t relax_after, char* report, int32_t cap) {
  try {
    if (!beads_path || !*beads_path) throw std::runtime_error("choose the file with the moved beads");
    const caps::Trajectory bt = caps::open_file(beads_path, "");
    const caps::System beads = bt.frame(bt.frames() - 1);   // the last frame of a bead trajectory
    caps::BackmapReport rep;
    caps::System s = caps::backmap(d->frame, beads, std::clamp(int(per_bead), 1, 100), &rep);
    std::vector<std::string> notes = rep.notes;
    if (relax_after) {   // the bonds between beads settle: push-off, then minimisation
      caps::RelaxOptions ro;
      ro.ftol = 1.0;
      ro.max_iterations = 3000;
      ro.energy = elec(ro.energy);
      caps::RelaxReport rr;
      caps::relax(s, ro, &rr);
      notes.push_back("relaxed: " + (rr.notes.empty() ? std::string() : rr.notes.front()));
    }
    auto* nd = doc_of_system(std::move(s), d, "model.backmap", "backmapped onto the beads of " + std::string(beads_path),
                             {{"per_bead", std::to_string(per_bead)}, {"beads", std::to_string(rep.beads)}, {"relaxed", relax_after ? "yes" : "no"}});
    std::string t;
    for (const auto& n : notes) t += n + "\n";
    report_out(t, report, cap);
    return nd;
  } catch (const std::exception& e) {
    g_error = e.what();
    return nullptr;
  }
}

extern "C" caps_doc* caps_pipeline_materialize(caps_doc* d) {
  try {
    if (!d->pstate) throw std::runtime_error("no pipeline result to make real");
    caps::System sys = d->pstate->system;
    int64_t id = 0;
    for (auto& a : sys.atoms) a.id = ++id;   // unique identifiers in order
    return doc_of_system(std::move(sys), d, "pipeline.materialize", "the pipeline's particles made into a structure (" + std::to_string(d->pstate->system.atoms.size()) + " atoms)",
                         {{"steps", std::to_string(d->pipeline ? d->pipeline->steps.size() : 0)}});
  } catch (const std::exception& e) {
    g_error = e.what();
    return nullptr;
  }
}

extern "C" int32_t caps_expression_count(caps_doc* d, const char* expr, char* out, int32_t cap) {
  caps::Json r = caps::Json::object();
  try {
    caps::Pipeline p;
    caps::PipelineStep st;
    st.type = "select_expression";
    st.params["expression"] = std::string(expr ? expr : "");
    p.steps.push_back(st);
    const auto res = caps::run_pipeline(d->frame, p, int(d->current), 0, &d->traj);
    r["ok"] = true;
    r["count"] = double(res.selected_count());
    r["total"] = double(d->frame.atoms.size());
  } catch (const std::exception& e) {
    r["ok"] = false;
    r["error"] = std::string(e.what());
  }
  caps::Json types = caps::Json::object();   // type labels, to say what "Type == 2" means
  for (const auto& t : d->frame.types) if (!t.label.empty()) types[std::to_string(t.type)] = t.label;
  r["types"] = std::move(types);
  return report_out(r.dump(0), out, cap);
}
