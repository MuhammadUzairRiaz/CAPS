#include "caps_c.h"

#include <cstdio>
#include <algorithm>
#include <deque>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <filesystem>
#include <thread>
#include <fstream>
#include <string>

#include "caps/probe.hpp"
#include "caps/normal_modes.hpp"
#include "caps/conformers.hpp"
#include "caps/ibi.hpp"
#include "caps/free_energy.hpp"
#include "caps/bondrules.hpp"
#include "caps/piece.hpp"
#include "caps/tags.hpp"
#include "caps/amber.hpp"
#include "caps/cg_map.hpp"
#include "caps/analysis.hpp"
#include "caps/adsorption.hpp"
#include "caps/cbmc.hpp"
#include "caps/dpd.hpp"
#include "caps/functionalize.hpp"
#include "caps/layers.hpp"
#include "caps/ffmerge.hpp"
#include "caps/water.hpp"
#include "caps/ffio.hpp"
#include "caps/manybody.hpp"
#include "caps/molecule.hpp"
#include "caps/sorption.hpp"
#include "caps/dlpoly.hpp"
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
#include "caps/bond_react.hpp"
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
#include "caps/raytrace.hpp"
#include "caps/scene_export.hpp"
#include "caps/query.hpp"
#include "caps/charges.hpp"
#include "caps/molinfo.hpp"
#include "caps/yaml.hpp"
#include "caps/voids.hpp"
#include "caps/kremer_grest.hpp"
#include "caps/nano.hpp"
#include "caps/json.hpp"
#include "caps/labels.hpp"
#include "caps/lattice.hpp"

#include <map>
#include <unordered_map>
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
  bool ua_summed = false;                    // automatic on a united-atom structure: the fallback is the Gasteiger charges summed into its sites
  caps::TypingResult typing;
  std::vector<std::string> types;
  caps::ParamReport rep;
  std::shared_ptr<const caps::ForceField> ff;   // null while atoms are untyped
  std::vector<caps::LammpsStyle::Group> group_atoms;   // by-group force fields: each group's name and atoms (LAMMPS groups)
  bool complete = false;
  std::vector<std::string> prep_notes;       // what was done to the structure for the force field (united-atom sites)
  std::string report;                        // JSON, see caps_field_report
  // the file's own types, restored by caps_field_clear
  std::vector<std::pair<int, std::string>> file_types;
  std::vector<caps::TypeInfo> file_type_table;
  std::vector<double> file_charges;
  bool file_has_charges = false;
  std::string groups;                        // v36: a force field per group (caps_field_assign_groups), as JSON; "" one for all
  std::string model;                         // a model's own force field built with the structure (Kremer–Grest), as JSON
  std::string mixing;                        // v47: the mixing rule for unlike Lennard-Jones pairs in place of the force field's ("" its own)
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
  std::vector<double> preview_q;   // colour by these charges instead of the structure's (Charges page, before Apply)
  std::string column;              // colour 5: a per-atom column of the trajectory (a LAMMPS dump's c_pe, |f| …)
  // for the current frame
  std::vector<uint8_t> style;      // 255: the view's style
  std::unique_ptr<caps::Mesh> mesh, poly;
  std::vector<caps::Segment> ribbon;
  double phi_lo = 0, phi_hi = 0;
  std::string error;
};

struct caps_doc {
  // typing by hand learned for a by-group assignment (caps_field_groups_by_example): the example chain, its types and the
  // force field they belong to; every group assigned with that force field is typed from it, again at each assignment
  std::shared_ptr<caps::System> hand_example;
  std::vector<std::string> hand_types;
  std::string hand_ff;
  caps::Json hand_report;
  caps::Trajectory traj;
  caps::System frame;
  caps::Renderer renderer;
  std::vector<double> dcom;   // distance of each atom to its own molecule's centre of mass
  size_t current = 0;
  bool wrap = false;
  std::unique_ptr<FieldState> field;
  std::string analysis;   // last caps_analyze result (JSON)
  std::string eq_checks;  // last caps_equilibrate convergence checks (JSON)
  std::string react_json; // last caps_react network summary (JSON)
  std::vector<int64_t> react_chains;   // each atom's chain as the last caps_react left it (the next run starts from them)
  size_t react_bonds = 0;              // the bond count they belong to (an edit in between drops them)
  int64_t held_mol = 0;   // molecule held in place by caps_relax (0: none)
  std::string pack_items; // v51 the packing that made it: per input molecule its name, molecule ids and own force field (JSON)
  std::vector<uint32_t> fixed_atoms;   // v36: atoms held in place besides the held molecule (frame indices)
  int fixed_axes = 7;                  // v62: which coordinates of the fixed atoms are held (bits x 1, y 2, z 4)
  std::vector<int64_t> rigid_mols;     // v62: molecules written as rigid bodies in LAMMPS inputs
  std::vector<caps::RelaxOptions::Restraint> restraints;   // distance restraints for caps_relax
  std::vector<caps::RelaxOptions::DihedralRestraint> dihedral_restraints;   // and dihedral ones
  double ph = -1;         // Add hydrogens: residues protonated at this pH (< 0: neutral valences)
  std::unique_ptr<caps::Pipeline> pipeline;        // caps_pipeline_set: steps run on every shown frame
  std::unique_ptr<caps::PipelineState> pstate;     // its result for the current frame
  // results cached per frame (design/boards/PipelineSteps "results cached per frame"): kept while the pipeline and the
  // frame's contents (a fingerprint of positions, bonds, cell and the trajectory's length) are unchanged
  struct PipeCached { size_t frame; uint64_t fp, key; std::shared_ptr<const caps::PipelineState> st; };   // key: the steps it ran
  std::deque<PipeCached> pcache;
  bool pstate_cached = false;
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
  // Look (design/boards/Look): sizes for the view, and the selection's own size factor
  struct Sizes {
    double atom_scale = 0.28, bond_radius = 0.14, space_scale = 1.0, line_px = 1.4;
    bool bond_orders = false;
    std::vector<float> factor;   // per atom (empty: 1)
  } sizes;
  // Probes (design/boards/Probes): drawn in the view from their atoms in the frame shown
  struct ProbeDef { caps::ProbeKind kind; std::vector<size_t> atoms; unsigned rgb; };
  std::vector<ProbeDef> probes;
  // caps_set_atom_state (design/boards/SelectionBar, Layers): per atom 0 shown, 1 ghost (faint, not pickable), 2 hidden. View
  // only: the structure, its exports and calculations keep every atom. Ignored once the atom count no longer matches.
  std::vector<uint8_t> atom_state;
  int vision = 0;                                  // caps_set_vision: the view as seen with a colour-vision deficiency
  double vision_severity = 1.0;
  std::array<float, 9> scene_vision{};   // the scene's copy of the vision matrix (caps_scene.vision_matrix)
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

// what the pipeline's result on this frame depends on besides the steps: FNV-1a over positions, bonds, cell and the
// trajectory's length (steps read other frames)
uint64_t frame_fingerprint(const caps_doc* d) {
  uint64_t h = 1469598103934665603ull;
  auto mix = [&](const void* p, size_t n) {
    const auto* b = static_cast<const unsigned char*>(p);
    for (size_t k = 0; k < n; ++k) h = (h ^ b[k]) * 1099511628211ull;
  };
  const size_t n = d->frame.atoms.size(), nb = d->frame.bonds.size(), nf = d->traj.frames();
  mix(&n, sizeof n), mix(&nb, sizeof nb), mix(&nf, sizeof nf);
  for (const auto& a : d->frame.atoms) mix(a.pos.data(), sizeof(double) * 3), mix(&a.element, sizeof a.element), mix(&a.charge, sizeof a.charge);
  for (const auto& b : d->frame.bonds) mix(&b.i, sizeof b.i), mix(&b.j, sizeof b.j);
  for (const auto* v : {&d->frame.cell.a, &d->frame.cell.b, &d->frame.cell.c, &d->frame.cell.origin}) mix(v->data(), sizeof(double) * 3);
  return h;
}

// the cache's key of a pipeline's steps [lo, n) with its branch
uint64_t steps_key(const caps::Pipeline& p, size_t lo) {
  caps::Pipeline part;
  part.branch = lo == 0 ? p.branch : "";
  part.steps.assign(p.steps.begin() + long(lo), p.steps.end());
  return std::hash<std::string>{}(caps::pipeline_to_json(part).dump(0));
}

// Molecule numbers as CAPS counts them everywhere (groups, the held molecule, the LAMMPS molecule column): the file's
// own when it has them, else the bonded fragments in order (1-based); a PDB's chain ids are not molecules.
std::vector<int64_t> molecule_ids(const caps::System& s) {
  std::vector<int64_t> r(s.atoms.size());
  if (s.has_mol) {
    for (size_t i = 0; i < r.size(); ++i) r[i] = s.atoms[i].mol;
    return r;
  }
  const auto c = s.molecules();
  for (size_t i = 0; i < r.size(); ++i) r[i] = int64_t(c[i]) + 1;
  return r;
}

// the atoms held in place: the held molecule and the fixed atoms (empty: none)
std::vector<char> fixed_mask(const caps_doc* d, const caps::System& s) {
  std::vector<char> m;
  if (d->held_mol <= 0 && d->fixed_atoms.empty()) return m;
  m.assign(s.atoms.size(), 0);
  if (d->held_mol > 0) {
    const auto ids = molecule_ids(s);
    for (size_t i = 0; i < s.atoms.size(); ++i) m[i] = ids[i] == d->held_mol;
  }
  // the fixed atoms: every coordinate, or only some (bits 2, 4, 8 of the mask: x, y, z)
  const char code = d->fixed_axes == 7 ? 1 : char((d->fixed_axes & 7) << 1);
  for (uint32_t i : d->fixed_atoms) if (i < m.size() && m[i] != 1) m[i] = code;
  return m;
}

// CAPS does not evaluate a literature many-body potential (Tersoff, EAM …): its atoms must be held for CAPS's own runs
void require_manybody_held(const caps_doc* d, const std::vector<char>& m) {
  if (!d->field || !d->field->ff || !d->field->ff->manybody.on()) return;
  const caps::ForceField& ff = *d->field->ff;
  size_t loose = 0;
  for (size_t i = 0; i < ff.type_index.size(); ++i) {
    const size_t t = size_t(ff.type_index[i]);
    if (t < ff.manybody.element.size() && !ff.manybody.element[t].empty() && (i >= m.size() || !caps::holds_all(m[i]))) ++loose;
  }
  if (loose)
    throw std::runtime_error(std::to_string(loose) + " atoms are under the " + ff.manybody.style + " potential, which CAPS does not evaluate (LAMMPS does): hold them "
                             "(Hold the filler, or hold the selection) to run in CAPS, or run the LAMMPS files");
}

// GROMACS: the held atoms as an index group (STEM.ndx: System and Frozen) and freezegrps / freezedim in STEM.mdp
// Returns what was written ("" nothing): the held atoms as a freeze group, the tags as index groups.
std::string write_gromacs_freeze(const caps_doc* d, const caps::System& s, const std::string& stem) {
  const auto m = fixed_mask(d, s);
  size_t n = 0;
  for (char c : m) n += c != 0;
  if (n == 0 && s.tags.empty()) return "";
  const std::string file = std::filesystem::path(stem).filename().string() + ".ndx";
  const std::string tagged = s.tags.empty() ? "" : std::to_string(s.tags.size()) + " tag" + (s.tags.size() == 1 ? "" : "s") + " as index groups in " + file + " (gmx select, grompp -n)";
  std::ofstream ndx(stem + ".ndx");
  auto group = [&](const std::string& name, auto in) {
    ndx << "[ " << name << " ]\n";
    int col = 0;
    for (size_t i = 0; i < s.atoms.size(); ++i)
      if (in(i)) { ndx << (i + 1) << (++col % 15 == 0 ? "\n" : " "); }
    ndx << "\n";
  };
  group("System", [](size_t) { return true; });
  size_t nfull = 0, npart = 0;
  for (char c : m) nfull += c == 1, npart += c != 0 && c != 1;
  if (nfull > 0) group("Frozen", [&](size_t i) { return m[i] == 1; });
  if (npart > 0) group("FrozenAxes", [&](size_t i) { return m[i] != 0 && m[i] != 1; });
  // the tags as index groups (gmx select / grompp -n)
  for (size_t k = 0; k < s.tags.size(); ++k) group(caps::tag_group_name(s.tags[k].name), [&](size_t i) { return (s.atoms[i].tags >> k) & 1u; });
  if (n == 0) return tagged;
  std::ofstream mdp(stem + ".mdp", std::ios::app);
  const std::string axes = std::string(d->fixed_axes & 1 ? "Y" : "N") + " " + (d->fixed_axes & 2 ? "Y" : "N") + " " + (d->fixed_axes & 4 ? "Y" : "N");
  mdp << "\n; atoms held in place in CAPS (index groups in " << std::filesystem::path(stem).filename().string() << ".ndx: grompp -n)\nfreezegrps = "
      << (nfull ? "Frozen" : "") << (nfull && npart ? " " : "") << (npart ? "FrozenAxes" : "") << "\nfreezedim = " << (nfull ? "Y Y Y" : "")
      << (nfull && npart ? " " : "") << (npart ? axes : "") << "\n";
  return "the held atoms are a freeze group: " + file + " and freezegrps in the .mdp (grompp -n)" + (tagged.empty() ? "" : "; " + tagged);
}

void run_doc_pipeline(caps_doc* d) {
  d->pstate.reset();
  d->shown_of.clear();
  d->pstate_cached = false;
  if (!d->pipeline) { d->pcache.clear(); return; }
  const caps::Pipeline& P = *d->pipeline;
  const uint64_t fp = frame_fingerprint(d), full = steps_key(P, 0);
  auto find = [&](uint64_t key) -> std::shared_ptr<const caps::PipelineState> {
    for (const auto& e : d->pcache) if (e.frame == d->current && e.fp == fp && e.key == key) return e.st;
    return nullptr;
  };
  auto keep = [&](uint64_t key, const caps::PipelineState& st) {
    // an older entry of this frame and steps goes; up to 256 entries and 4 million particles, oldest out first
    for (auto it = d->pcache.begin(); it != d->pcache.end();) it = it->frame == d->current && it->key == key ? d->pcache.erase(it) : it + 1;
    d->pcache.push_back({d->current, fp, key, std::make_shared<const caps::PipelineState>(st)});
    size_t total = 0;
    for (const auto& e : d->pcache) total += e.st->system.atoms.size();
    while (d->pcache.size() > 1 && (d->pcache.size() > 256 || total > 4000000)) {
      total -= d->pcache.front().st->system.atoms.size();
      d->pcache.pop_front();
    }
  };
  if (auto hit = find(full)) {
    d->pstate = std::make_unique<caps::PipelineState>(*hit);
    d->pstate_cached = true;
  } else {
    const int64_t ts = d->current < d->traj.timesteps.size() ? d->traj.timesteps[d->current] : 0;
    const size_t n = P.steps.size(), trunk = caps::pipeline_trunk(P);
    caps::PipelineState st;
    if (trunk > 0 && trunk < n) {   // branches: the shared steps once per frame, the shown branch on top of them
      const uint64_t tkey = steps_key(P, trunk);
      if (auto t = find(tkey)) st = *t;
      else {
        st = caps::pipeline_begin(d->frame, P, int(d->current), ts, &d->traj);
        caps::pipeline_run_steps(st, P, n, trunk);
        st.pipeline = nullptr;
        keep(tkey, st);
      }
      st.traj = &d->traj;
      caps::pipeline_run_steps(st, P, trunk, 0);
      caps::pipeline_finish(st);
    } else {
      st = caps::run_pipeline(d->frame, P, int(d->current), ts, &d->traj);
    }
    keep(full, st);
    d->pstate = std::make_unique<caps::PipelineState>(std::move(st));
  }
  d->pstate->traj = &d->traj;
  d->pstate->pipeline = nullptr;
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
  k.fov_deg = c->fov_deg > 0 ? std::clamp(c->fov_deg, 10.0, 120.0) : 35.0;
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
  r.outline_strength = o->outlines >= 2 ? 2 : 1;   // 2: high-contrast outlines
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
  r.atom_scale = d->sizes.atom_scale, r.bond_radius = d->sizes.bond_radius, r.space_scale = d->sizes.space_scale, r.line_px = d->sizes.line_px;
  r.bond_orders = d->sizes.bond_orders;
  if (!d->pstate && d->sizes.factor.size() == d->frame.atoms.size()) r.size_factor = d->sizes.factor;
  if (!d->pstate && d->selection.size() == d->frame.atoms.size()) {   // the selection ringed (up to 50 000 atoms)
    for (size_t i = 0; i < d->selection.size() && r.highlight.size() < 50000; ++i) if (d->selection[i]) r.highlight.push_back(int(i));
  }
  if (!d->pstate) { r.segments = d->overlay; r.segments.insert(r.segments.end(), d->checks.begin(), d->checks.end()); }
  if (!d->pstate) for (const auto& pd : d->probes) {
      caps::Probe p;
      try { p = caps::make_probe(d->frame, pd.atoms, pd.kind); } catch (...) { continue; }
      auto tube = [&](const caps::Vec3& a, const caps::Vec3& b, double rad, bool arrow = false) { r.segments.push_back({a, b, pd.rgb, rad, arrow}); };
      const caps::Vec3 c = p.centre, e0 = p.axes[0], e1 = p.axes[1], e2 = p.axes[2];
      if (pd.kind == caps::ProbeKind::Point) {
        for (const auto& e : p.axes) tube(c - e * 0.7, c + e * 0.7, 0.1);
      } else if (pd.kind == caps::ProbeKind::Axis) {
        tube(c - e0 * (p.semi[0] + 1.0), c + e0 * (p.semi[0] + 1.0), 0.09, true);
      } else if (pd.kind == caps::ProbeKind::Plane) {
        const double u = p.semi[0] + 1.0, v = p.semi[1] + 1.0;
        const caps::Vec3 q[4] = {c + e0 * u + e1 * v, c - e0 * u + e1 * v, c - e0 * u - e1 * v, c + e0 * u - e1 * v};
        for (int k = 0; k < 4; ++k) tube(q[k], q[(k + 1) % 4], 0.06);
        tube(c, c + e2 * 3.0, 0.08, true);   // its normal: "above"
      } else {
        const caps::Vec3* ax[3] = {&e0, &e1, &e2};
        for (int k = 0; k < 3; ++k) {   // three great circles in the principal planes
          const int i = k, j = (k + 1) % 3;
          const double ai = p.semi[size_t(i)], aj = p.semi[size_t(j)];
          for (int t = 0; t < 36; ++t) {
            const double t0 = 2 * M_PI * t / 36, t1 = 2 * M_PI * (t + 1) / 36;
            tube(c + *ax[i] * (ai * std::cos(t0)) + *ax[j] * (aj * std::sin(t0)), c + *ax[i] * (ai * std::cos(t1)) + *ax[j] * (aj * std::sin(t1)), 0.05);
          }
        }
      }
    }
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
    if (L.colour == 5) {
      // a dump column in the frame shown, blue (low) to orange (high)
      r.colour_by = caps::ColourBy::Property;
      r.property.clear();
      auto it = d->traj.columns.find(L.column);
      if (it != d->traj.columns.end() && d->current < it->second.size() && it->second[d->current].size() == d->frame.atoms.size())
        for (float v : it->second[d->current]) r.property.push_back(double(v));
      else r.property.assign(d->frame.atoms.size(), 0.0);
      r.ramp = L.ramp;
      r.symmetric = false;
    } else if (L.colour == 4) {
      r.colour_by = caps::ColourBy::Property;
      r.property.clear();
      if (L.preview_q.size() == d->frame.atoms.size()) r.property = L.preview_q;
      else for (const auto& a : d->frame.atoms) r.property.push_back(a.charge);
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
    for (const auto& pm : st.meshes) if (pm.mesh) r.meshes.push_back({pm.mesh.get(), pm.rgb, pm.opacity});
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
  // hidden and ghosted atoms (after every style: they win)
  if (!d->pstate && d->atom_state.size() == d->frame.atoms.size() &&
      std::any_of(d->atom_state.begin(), d->atom_state.end(), [](uint8_t v) { return v != 0; })) {
    const size_t n = d->atom_state.size();
    if (r.atom_style.size() != n) r.atom_style.assign(n, uint8_t(r.style));
    r.unpickable.assign(n, 0);
    for (size_t i = 0; i < n; ++i) {
      const uint8_t v = d->atom_state[i] & 3;
      if (d->atom_state[i] & 4) r.unpickable[i] = 1;   // locked: drawn as it is, never picked
      if (v == 2) r.atom_style[i] = uint8_t(caps::Style::Hidden), r.unpickable[i] = 1;
      else if (v == 1) {
        if (r.transparency.size() != n) r.transparency.assign(n, 0.0f);
        r.transparency[i] = std::max(r.transparency[i], 0.78f);
        r.unpickable[i] = 1;
      }
    }
    r.highlight.erase(std::remove_if(r.highlight.begin(), r.highlight.end(),
                                     [&](int i) { return i >= 0 && size_t(i) < n && (d->atom_state[size_t(i)] & 3) != 0; }),
                      r.highlight.end());
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
void field_run_groups(caps_doc* d);
void field_run_model(caps_doc* d);
void install_file_field(caps_doc* d);

void field_run(caps_doc* d) {
  if (!d->field->groups.empty()) { field_run_groups(d); return; }
  if (!d->field->model.empty()) { field_run_model(d); return; }
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
          F.ff = std::make_shared<caps::ForceField>(caps::parameterize(s, def, F.types, F.ua_summed ? "keep" : "gasteiger", &F.rep, true));
          F.charges = F.ua_summed ? "keep" : "gasteiger";
          // molecules whose types all carry charges (TraPPE's CO2 beside its united-atom alkanes) keep the force field's own
          int nm = 0;
          const auto mol = s.molecules(&nm);
          std::vector<char> full(size_t(std::max(nm, 0)), 1);
          std::vector<double> qt(n, 0.0);
          for (size_t i = 0; i < n; ++i) {
            if (mol[i] < 0) continue;
            const caps::FFType* t = i < F.types.size() && !F.types[i].empty() ? def.type(F.types[i]) : nullptr;
            if (!t) continue;   // a hydrogen folded into its carbon (united atom): no site of its own
            if (std::isnan(t->charge)) full[size_t(mol[i])] = 0;
            else qt[i] = t->charge;
          }
          int own = 0;
          for (char f : full) own += f;
          if (own > 0) {
            caps::System sq = s;
            for (size_t i = 0; i < n; ++i) sq.atoms[i].charge = mol[i] >= 0 && full[size_t(mol[i])] ? qt[i] : F.ff->charge[i];
            sq.has_charges = true;
            caps::ParamReport rep2;
            F.ff = std::make_shared<caps::ForceField>(caps::parameterize(sq, def, F.types, "keep", &rep2, true));
            F.rep = std::move(rep2);
            F.charges = "mixed";
            F.rep.notes.push_back(def.name + "'s own charges on " + std::to_string(own) + " of " + std::to_string(nm) +
                                  " molecules (every type of theirs carries one); " + (F.ua_summed ? "Gasteiger–Marsili summed into the sites" : "Gasteiger–Marsili") +
                                  " on the others, whose types carry none (" + why.substr(0, why.find(';')) + ") (charges: automatic)");
          } else
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
  if (F.ff)
    for (const auto& note : caps::ring_angle_notes(s, *F.ff)) F.rep.notes.push_back(note);
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
  if (F.ff && !F.mixing.empty() && F.mixing != F.ff->mixing) {
    // the rule chosen in place of the force field's own: every run and export uses it; explicit pairs still win
    auto ff2 = std::make_shared<caps::ForceField>(*F.ff);
    F.rep.notes.push_back("mixing rule " + F.mixing + " in place of " + F.ff->name + "'s own " + F.ff->mixing +
                          " for unlike Lennard-Jones pairs (chosen in CAPS, not the published force field); its explicit pairs still apply");
    ff2->mixing = F.mixing;
    F.ff = std::move(ff2);
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
  r["mixing"] = def.pair_table.empty() ? (F.ff ? F.ff->mixing : def.mixing) : std::string("none: every pair from its table");
  r["mixing_own"] = def.mixing;
  r["mixing_override"] = F.mixing;
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
  // the structure's tags, saved beside it (NAME.tags.json), for the same atoms
  if (caps::read_tags(d->traj.topology, path) && d->frame.atoms.size() == d->traj.topology.atoms.size()) {
    d->frame.tags = d->traj.topology.tags;
    for (size_t i = 0; i < d->frame.atoms.size(); ++i) d->frame.atoms[i].tags = d->traj.topology.atoms[i].tags;
  }
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
    install_file_field(d);
    return d;
  } catch (const std::exception& e) {
    g_error = e.what();
  } catch (...) {
    g_error = "unknown error";
  }
  return nullptr;
}

caps_doc* caps_open_many(const char* paths_json, const char* topology_path) {
  try {
    const caps::Json j = caps::Json::parse(paths_json ? paths_json : "[]");
    std::vector<std::string> paths;
    for (const auto& x : j.items()) paths.push_back(x.str());
    if (paths.empty()) throw std::invalid_argument("no files");
    auto* d = new caps_doc;
    std::vector<std::string> notes;
    d->traj = caps::open_files(paths, topology_path ? topology_path : "", &notes);
    for (const auto& n : notes) d->traj.topology.notes.push_back(n);
    refresh(d);
    prov_opened(d, paths[0], topology_path ? topology_path : "", "io.read");
    install_file_field(d);
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
    install_file_field(d);
    return d;
  } catch (const std::exception& e) {
    g_error = e.what();
  } catch (...) {
    g_error = "unknown error";
  }
  return nullptr;
}

caps_doc* caps_open_frames(const char* path, const char* topology_path, int64_t first, int64_t last, int64_t stride, int32_t max_frames,
                           caps_open_progress_fn progress, void* user) {
  try {
    caps::OpenProgress p;
    p.max_frames = max_frames > 0 ? static_cast<size_t>(max_frames) : 0;
    p.frames.first = first > 0 ? static_cast<size_t>(first) : 0;
    p.frames.last = last >= 0 ? static_cast<size_t>(last) : SIZE_MAX;
    p.frames.stride = stride > 1 ? static_cast<size_t>(stride) : 1;
    if (progress) p.report = [&](int stage, double f, const std::string& detail) { return progress(stage, f, detail.c_str(), user) == 0; };
    auto* d = new caps_doc;
    try {
      d->traj = caps::open_file(path, topology_path ? topology_path : "", p);
      if (!p.frames.all() && d->traj.frames_read == 0 && p.frames.first >= d->traj.frames())
        throw std::runtime_error(std::string(path) + ": frame " + std::to_string(p.frames.first) + " is past the file's " + std::to_string(d->traj.frames()));
    } catch (...) {
      delete d;
      throw;
    }
    refresh(d);
    prov_opened(d, path, topology_path ? topology_path : "", "io.read");
    install_file_field(d);
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
    sd->look.preview_q = d->look.preview_q;
    sd->look.ramp = d->look.ramp;
    sd->selection = d->selection;
    sd->cell_repeats = d->cell_repeats;
    sd->ph = d->ph;
    sd->held_mol = d->held_mol;
    sd->fixed_atoms = d->fixed_atoms;
    sd->fixed_axes = d->fixed_axes;
    sd->rigid_mols = d->rigid_mols;
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

// POSCAR, CONTCAR (and POSCAR_…, CONTCAR-…) by the file's name, as VASP names them
bool is_poscar_name(const std::string& p) {
  const auto slash = p.find_last_of("/\\");
  std::string n = slash == std::string::npos ? p : p.substr(slash + 1);
  for (char& c : n) c = char(std::toupper(static_cast<unsigned char>(c)));
  return n.rfind("POSCAR", 0) == 0 || n.rfind("CONTCAR", 0) == 0;
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
    else if (ends(".vasp") || ends(".poscar") || is_poscar_name(p)) caps::write_poscar(d->frame, p, fixed_mask(d, d->frame));
    else if (ends(".caps.data")) caps::write_lammps_data(d->frame, p);   // CAPS to CAPS: every atom (a four-site water's M), no coefficients
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
    try { caps::write_tags(d->frame, p); } catch (...) {}   // the tags beside it (NAME.tags.json), as the provenance
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

// The LAMMPS groups of an input: the Field page's groups (a force field per group), else a composite's filler (the held
// molecule) and matrix; none otherwise.
static std::vector<caps::LammpsStyle::Group> lammps_groups(const caps_doc* d, const caps::System& s) {
  if (d->field) {
    const auto& ga = d->field->group_atoms;
    if (ga.size() > 1 && std::all_of(ga.begin(), ga.end(), [&](const auto& g) {
          return std::all_of(g.atoms.begin(), g.atoms.end(), [&](uint32_t i) { return i < s.atoms.size(); }); }))
      return ga;
  }
  if (d->held_mol > 0) {
    caps::LammpsStyle::Group filler{"filler", {}}, matrix{"matrix", {}};
    for (size_t i = 0; i < s.atoms.size(); ++i) (s.atoms[i].mol == d->held_mol ? filler : matrix).atoms.push_back(uint32_t(i));
    if (!filler.atoms.empty() && !matrix.atoms.empty()) return {filler, matrix};
  }
  return {};
}

int32_t caps_lammps_input(caps_doc* d, const char* data_name, char* text, int32_t cap) {
  return guard([&] {
    caps::ForceField ff;
    if (d->field && d->field->complete) ff = *d->field->ff;
    else ff = default_ff(d->frame);
    const auto tmp = std::filesystem::temp_directory_path() / ("caps_input_" + std::to_string(reinterpret_cast<uintptr_t>(d)) + ".in");
    caps::LammpsStyle style;
    style.groups = lammps_groups(d, d->frame);
    style.rigid_mols = d->rigid_mols;
    caps::write_lammps_input(d->frame, ff, elec(), data_name && *data_name ? data_name : "system.data", tmp.string(), d->held_mol, false, {}, style);
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
    if (stem && *stem) { if (auto w = write_gromacs_freeze(d, d->frame, stem); !w.empty()) notes.push_back(w); }
    else if (!(stem && *stem) && !fixed_mask(d, d->frame).empty()) notes.push_back("the held atoms go into STEM.ndx as a freeze group when the files are written");
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
extern "C" int32_t caps_field_save(caps_doc* d, const char* path) {
  return guard([&] {
    if (!d || !d->field || !d->field->ff) throw caps::FFError("no force field assigned");
    if (d->field->ff->type_index.size() != d->traj.topology.atoms.size()) throw caps::FFError("the force field is for another structure");
    std::ofstream f(path ? path : "");
    if (!f) throw std::runtime_error(std::string("cannot write ") + (path ? path : ""));
    caps::Json j = caps::Json::parse(caps::forcefield_to_json(*d->field->ff));
    caps::Json g = caps::Json::array();   // the groups it was assigned by (LAMMPS groups in the inputs)
    for (const auto& grp : d->field->group_atoms) {
      caps::Json o = caps::Json::object(), a = caps::Json::array();
      for (auto i : grp.atoms) a.push_back(caps::Json(double(i)));
      o["name"] = grp.name, o["atoms"] = a;
      g.push_back(o);
    }
    j["groups"] = g;
    f << j.dump_exact(0) << "\n";
    return 0;
  });
}

extern "C" int32_t caps_field_load(caps_doc* d, const char* path) {
  return guard([&] {
    std::ifstream f(path ? path : "");
    if (!f) throw std::runtime_error(std::string("cannot open ") + (path ? path : ""));
    std::stringstream ss;
    ss << f.rdbuf();
    auto ff = std::make_shared<caps::ForceField>(caps::forcefield_from_json(ss.str()));
    if (ff->type_index.size() != d->traj.topology.atoms.size())
      throw caps::FFError("the force field file is for " + std::to_string(ff->type_index.size()) + " atoms, the structure has " + std::to_string(d->traj.topology.atoms.size()));
    d->traj.topology.forcefield = ff;
    for (size_t i = 0; i < d->traj.topology.atoms.size(); ++i) d->traj.topology.atoms[i].charge = ff->charge[i];
    d->traj.topology.has_charges = true;
    install_file_field(d);
    if (!d->field || !d->field->ff) throw caps::FFError("the force field file could not be installed");
    {   // the groups it was assigned by
      const caps::Json j = caps::Json::parse(ss.str());
      if (j.has("groups") && j["groups"].is_array())
        for (const auto& g : j["groups"].items()) {
          caps::LammpsStyle::Group grp{g.text("name", "group"), {}};
          for (const auto& a : g["atoms"].items()) grp.atoms.push_back(uint32_t(a.number()));
          d->field->group_atoms.push_back(std::move(grp));
        }
    }
    prov_step(d, "field.assign", ff->name, {{"force field", "read whole from " + std::filesystem::path(path).filename().string()}}, "", {}, {});
    return 0;
  });
}

extern "C" int32_t caps_pack_items(caps_doc* d, char* json, int32_t cap) {
  if (!d) return -1;
  return report_out(d->pack_items.empty() ? "[]" : d->pack_items, json, cap);
}

extern "C" int32_t caps_pack_convert(const char* text, int32_t to_caps, char* out, int32_t cap) {
  try {
    const std::string t = text ? text : "";
    std::string r;
    if (to_caps) r = caps::is_caps_pack_input(t) ? t : caps::packmol_to_caps_pack(t);
    else {   // for packmol itself: CAPS's own lines (a molecule's force field or water model) left as comments
      r = caps::is_caps_pack_input(t) ? caps::caps_pack_to_packmol(t) : t;
      std::istringstream in(r);
      std::string out, line;
      while (std::getline(in, line)) {
        const auto p = line.find_first_not_of(" \t");
        const std::string w = p == std::string::npos ? "" : line.substr(p, line.find_first_of(" \t", p) - p);
        out += (w == "forcefield" || w == "water" || w == "compress" ? "# CAPS: " + line.substr(p) : line) + "\n";
      }
      r = out;
    }
    return report_out(r, out, cap);
  } catch (const std::exception& e) {
    g_error = e.what();
    return -1;
  }
}

extern "C" int32_t caps_water_models(char* json, int32_t cap) {
  try {
    caps::Json a = caps::Json::array();
    for (const auto& m : caps::water_models()) {
      caps::Json o = caps::Json::object();
      o["id"] = m.id, o["name"] = m.name, o["citation"] = m.citation, o["sites"] = double(m.sites), o["r_oh"] = m.r_oh, o["theta"] = m.theta;
      o["q_h"] = m.q_h, o["q_neg"] = m.q_neg, o["d_om"] = m.d_om, o["d_ol"] = m.d_ol, o["theta_l"] = m.theta_l, o["eps_o"] = m.eps_o, o["sigma_o"] = m.sigma_o, o["eps_h"] = m.eps_h, o["sigma_h"] = m.sigma_h;
      o["rigid"] = m.rigid, o["note"] = m.note;
      a.push_back(o);
    }
    return report_out(a.dump(0), json, cap);
  } catch (const std::exception& e) {
    g_error = e.what();
    return -1;
  }
}

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
             : kind == "nvt" ? caps::LammpsRun::Kind::NVT : kind == "npt" ? caps::LammpsRun::Kind::NPT
             : kind == "tensile" ? caps::LammpsRun::Kind::Tensile : kind == "creep" ? caps::LammpsRun::Kind::Creep
             : kind == "shear" ? caps::LammpsRun::Kind::Shear : kind == "protocol" ? caps::LammpsRun::Kind::Protocol : caps::LammpsRun::Kind::Check;
    {   // tensile / creep / shear (LAMMPS)
      const std::string a = o.text("axis", "x");
      run.axis = a == "y" ? 1 : a == "z" ? 2 : 0;
      run.strain_rate = o.num("strain_rate", 1e-3);
      run.max_strain = o.num("max_strain", 0);
      run.stress_mpa = o.num("stress_mpa", 50);
      run.shear_rate = o.num("shear_rate", 0.01);
    }
    if (kind == "protocol") {   // a CAPS protocol by name (with this run's temperature and pressure) or as text
      const std::string p = o.text("protocol", "larsen21");
      if (p == "larsen21" || p == "annealing" || p == "pushoff") {
        caps::ProtocolParams pp;
        pp.t_final = o.num("temperature", 300), pp.p_final = o.num("pressure", 1.0);
        if (o.has("t_max")) pp.t_max = o.num("t_max", pp.t_max);
        if (o.has("p_max")) pp.p_max = o.num("p_max", pp.p_max);
        if (o.has("time_scale")) pp.time_scale = o.num("time_scale", 1.0);
        run.protocol = caps::protocol_text(caps::protocol_by_name(p, pp));
      } else
        run.protocol = p;
      run.production_ps = o.num("production_ps", 0);
      if (caps::parse_protocol(run.protocol).empty()) throw std::runtime_error("the protocol has no stages");
    }
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
    if (gromacs && (kind == "tensile" || kind == "creep" || kind == "shear"))
      notes.push_back("the " + kind + " protocol is written for LAMMPS; the GROMACS files are a single point of the same system");
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
    ls.units = o.text("units", "auto");   // real | metal | auto (metal when a potential is read in metal units only)
    if (ls.kspace_accuracy <= 0 || ls.cutoff < 0) throw std::runtime_error("the k-space accuracy and the cut-off must be positive");
    // a force field LAMMPS cannot express (GROMOS's reaction field …) refuses the LAMMPS files only: the GROMACS files
    // are still written, and the reason goes back as lammps_error
    std::string lammps_error;
    const bool kg_model = !d->field->model.empty() && caps::Json::parse(d->field->model).text("model", "") == "kremer-grest";
    if (lammps && kg_model) try {
      // a Kremer–Grest melt: its own LAMMPS deck (units lj, or real as mapped): push-off, then FENE + WCA under Langevin
      const caps::Json m = caps::Json::parse(d->field->model);
      caps::KgOptions ko;
      ko.chains = int(m.num("chains", 0)), ko.beads = int(m.num("beads", 0)), ko.density = m.num("density", 0.85), ko.k_theta = m.num("k_theta", 0);
      ko.seed = uint64_t(m.num("seed", 1));
      ko.sigma = m.num("sigma", 0), ko.temperature = m.num("temperature", 0), ko.bead_mass = m.num("bead_mass", 0);
      caps::System ks = s;
      if (ko.sigma > 0) {
        for (auto& a : ks.atoms) a.pos = a.pos * (1.0 / ko.sigma);
        ks.cell.a = ks.cell.a * (1.0 / ko.sigma), ks.cell.b = ks.cell.b * (1.0 / ko.sigma), ks.cell.c = ks.cell.c * (1.0 / ko.sigma), ks.cell.origin = ks.cell.origin * (1.0 / ko.sigma);
      }
      caps::write_kg_lammps(ks, ko, base, o.num("pushoff_steps", 20000), double(run.steps > 0 ? run.steps : 100000));
      written.push_back({stem + ".data", "beads, bonds (and angles), in σ units" + std::string(ko.sigma > 0 ? " scaled to Å" : "")});
      written.push_back({stem + ".in", "Kremer–Grest deck: push-off, then FENE + WCA, Langevin thermostat"});
      notes.push_back(caps::Json(std::string("LAMMPS: the Kremer–Grest melt's own deck (") + (ko.sigma > 0 ? "units real, as mapped" : "units lj") +
                                 "); the run options of this page other than the steps do not apply to it"));
    } catch (const std::exception& ex) {
      lammps_error = ex.what();
      if (!gromacs) throw;
    }
    if (lammps && !kg_model) try {
      std::vector<std::string> lnotes;
      ls.groups = lammps_groups(d, s);
      ls.rigid_mols = d->rigid_mols;
      if (!d->rigid_mols.empty()) notes.push_back(caps::Json("LAMMPS: " + std::to_string(d->rigid_mols.size()) + " molecule(s) move as rigid bodies (fix rigid/nvt/small)"));
      caps::write_lammps_data_ff(s, ff, e, base + ".data", false, ls);
      caps::write_lammps_input(s, ff, e, stem + ".data", base + ".in", d->held_mol, true, run, ls, &lnotes);
      for (const auto& n : lnotes) notes.push_back(caps::Json("LAMMPS: " + n));
      if (!d->fixed_atoms.empty()) notes.push_back(caps::Json("LAMMPS: the " + std::to_string(d->fixed_atoms.size()) + " fixed atoms besides the held molecule are not held in this input (the GROMACS files freeze them)"));
      written.push_back({stem + ".data", "atoms, bonds, masses and bonded coefficients"});
      written.push_back({stem + ".in", "styles, every pair_coeff and the run"});
      if (ff.manybody.on()) {
        const bool metal_in = caps::lammps_metal_units(ff, ls);
        written.push_back({caps::manybody_file_name(ff.manybody, metal_in ? "metal" : "real"), "the " + ff.manybody.style + " potential file (" +
                           (!metal_in && caps::manybody_caps_converts(ff.manybody.style) ? std::string("converted by CAPS to real units")
                            : ff.manybody.tagged ? std::string("as given") : std::string("as given, with its units on the first line")) + ")"});
        if (!ff.manybody.file2.empty()) written.push_back({caps::manybody_file2_name(ff.manybody), "the MEAM parameter file (as given)"});
        const bool metal = caps::lammps_metal_units(ff, ls);
        notes.push_back(caps::Json("LAMMPS: " + ff.manybody.style + " overlays the pair terms for its elements (pair_style hybrid/overlay); " +
                                   (!metal && caps::manybody_caps_converts(ff.manybody.style)
                                        ? std::string("real units: CAPS wrote a converted copy of the file (A, B and the ε's × 23.060549, as LAMMPS cannot convert it)")
                                        : "the file in " + ff.manybody.units + " units" + (ff.manybody.units == (metal ? "metal" : "real") ? "" : ", which LAMMPS converts"))));
      }
      if (caps::lammps_metal_units(ff, ls)) {
        notes.push_back(caps::Json("LAMMPS: written in metal units (eV, ps, bar)" + std::string(ff.manybody.metal_only ? ", as " + ff.manybody.style + " requires" : ", as asked") +
                                   ": every energy parameter of the force field divided by 23.060549 (LAMMPS's own factor)"));
      }
      if (flag("moltemplate", false)) {   // the same as a moltemplate system (moltemplate.sh -overlay-all system.lt)
        std::ofstream lt(base + ".lt");
        lt << caps::lammps_to_moltemplate(base + ".data", base + ".in", ff.name);
        written.push_back({stem + ".lt", "the same system for moltemplate (moltemplate.sh -overlay-all " + stem + ".lt)"});
      }
    } catch (const std::exception& ex) {
      lammps_error = ex.what();   // refused before its files (prepare checks the force field first)
      if (!gromacs) throw;
    }
    // a force field GROMACS cannot express (class II's 9-6 Lennard-Jones and cross terms …) refuses the GROMACS files
    // only: the LAMMPS files are still written, and the reason goes back as gromacs_error
    std::string gromacs_error;
    if (gromacs && kg_model) gromacs_error = "Kremer–Grest's FENE bond with its WCA core has no GROMACS form (GROMACS's FENE bond carries no repulsive core, and a bonded pair's WCA cannot be cut off in [ pairs ]): export to LAMMPS";
    else if (gromacs) try {
      for (const auto& n : caps::write_gromacs(s, ff, e, base)) notes.push_back("GROMACS: " + n);
      if (auto w = write_gromacs_freeze(d, s, base); !w.empty()) notes.push_back("GROMACS: " + w);
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
    std::string dlpoly_error;
    if (flag("dlpoly", false)) try {   // STEM_dlpoly/FIELD, CONFIG, CONTROL (DL_POLY 4, units kcal)
      caps::DlpolyOptions dop;
      dop.title = stem;
      dop.cutoff = e.cutoff;
      dop.temperature = run.temperature;
      dop.steps = long(run.steps);
      dop.timestep_fs = caps::lammps_timestep(run, ff);
      for (const auto& n : caps::write_dlpoly(s, ff, (folder / (stem + "_dlpoly")).string(), dop)) notes.push_back("DL_POLY: " + n);
      written.push_back({stem + "_dlpoly/FIELD", "DL_POLY force field (units kcal): molecular types, every term, van der Waals pairs"});
      written.push_back({stem + "_dlpoly/CONFIG", "DL_POLY coordinates about the cell centre"});
      written.push_back({stem + "_dlpoly/CONTROL", "a generic DL_POLY NVT run with the same cut-off (edit before use)"});
    } catch (const std::exception& ex) {
      dlpoly_error = ex.what();
    }
    std::string amber_error;
    if (flag("amber", false)) try {   // STEM.prmtop and STEM.inpcrd (AMBER, OpenMM, ParmEd)
      if (kg_model) throw std::runtime_error("Kremer–Grest's FENE bond has no AMBER form: export to LAMMPS");
      for (const auto& n : caps::write_amber(s, ff, base)) notes.push_back("AMBER: " + n);
      written.push_back({stem + ".prmtop", "AMBER topology: every term, a Lennard-Jones pair table (OpenMM, ParmEd, AMBER)"});
      written.push_back({stem + ".inpcrd", "AMBER coordinates and box"});
    } catch (const std::exception& ex) {
      amber_error = ex.what();
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
    if (!dlpoly_error.empty()) r["dlpoly_error"] = dlpoly_error;
    if (!amber_error.empty()) r["amber_error"] = amber_error;
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

int32_t caps_cbmc(caps_doc* d, const caps_cbmc_opts* o, caps_cbmc_progress_fn progress, void* user, char* report, int32_t cap) {
  return guard([&] {
    if (!d || !o) throw std::invalid_argument("no document or options");
    caps::System s = d->traj.frame(d->current);
    const auto fp = field_for_run(d);
    const caps::ForceField ff = fp ? *fp : default_ff(s);
    caps::CbmcOptions c;
    if (o->moves > 0) c.moves = o->moves;
    if (o->trials > 0) c.trials = o->trials;
    if (o->max_torsions > 0) c.max_torsions = o->max_torsions;
    if (o->temperature > 0) c.temperature = o->temperature;
    if (o->cutoff > 0) c.cutoff = o->cutoff;
    c.coulomb = o->coulomb != 0;
    c.seed = o->seed ? o->seed : 1;
    bool cancelled = false;
    if (progress) c.progress = [&](int done, int acc) { cancelled = progress(done, c.moves, acc, user) != 0; return !cancelled; };
    caps::Trajectory out;
    int every = 0;
    c.snapshot = [&](const std::vector<caps::Vec3>& x) {   // about every 5 % of the run
      if (++every % 5 != 0) return;
      out.positions.push_back(x);
      out.cells.push_back(s.cell);
      out.timesteps.push_back(int64_t(out.timesteps.size()));
    };
    caps::make_molecules_whole(s);
    out.topology = s;
    {
      std::vector<caps::Vec3> p;
      for (const auto& a : s.atoms) p.push_back(a.pos);
      out.positions.push_back(p), out.cells.push_back(s.cell), out.timesteps.push_back(0);
    }
    caps::CbmcReport rep;
    caps::cbmc_regrow(s, ff, c, &rep);
    if (cancelled) throw std::runtime_error("CBMC cancelled");
    {
      std::vector<caps::Vec3> p;
      for (const auto& a : s.atoms) p.push_back(a.pos);
      out.positions.push_back(p), out.cells.push_back(s.cell), out.timesteps.push_back(int64_t(out.timesteps.size()));
    }
    caps::KeyValues pr = {{"moves", std::to_string(rep.attempted) + " end regrowths, " + std::to_string(rep.accepted) + " accepted (" + g6(std::round(1000 * rep.acceptance) / 10) + " %)"},
                          {"trials", std::to_string(c.trials) + " torsions per bond, up to " + std::to_string(c.max_torsions) + " bonds per move"},
                          {"temperature", g6(c.temperature) + " K"},
                          {"trial energies", std::string("van der Waals") + (c.coulomb ? ", DSF electrostatics" : "") + ", torsions · cut-off " + g6(rep.cutoff) + " Å"},
                          {"⟨R²⟩ end-to-end", g6(rep.r2_before) + " → " + g6(rep.r2_after) + " Å²"},
                          {"force field", ff_label(d)}};
    prov_step(d, "cbmc.regrow", "Configurational-bias Monte Carlo regrowth of chain ends", std::move(pr), seeded(c.seed), {"siepmann1992", "rosenbluth1955"},
              energy_approx(rep.cutoff, c.coulomb, false, 1));
    out.topology.atoms = s.atoms;
    out.topology.velocities.clear();
    out.topology.unwrapped = true;
    out.topology.notes = rep.notes;
    d->traj = std::move(out);
    d->current = d->traj.frames() - 1;
    refresh(d);
    caps::Json r = caps::Json::object();
    r["attempted"] = double(rep.attempted);
    r["accepted"] = double(rep.accepted);
    r["acceptance"] = rep.acceptance;
    r["chains"] = double(rep.chains);
    r["r2_before"] = rep.r2_before;
    r["r2_after"] = rep.r2_after;
    r["energy_change"] = rep.energy_change;
    r["cutoff"] = rep.cutoff;
    r["seconds"] = rep.seconds;
    caps::Json nt = caps::Json::array();
    for (const auto& x : rep.notes) nt.push_back(x);
    r["notes"] = std::move(nt);
    report_out(r.dump(), report, cap);
    return 0;
  });
}

namespace {
void push_undo(caps_doc* d, const std::string& what);   // with the edits, below
caps_doc* doc_of_system(caps::System sys, const caps_doc* from, const std::string& engine, const std::string& summary, caps::KeyValues params);

// The force field for s — the document's structure with molecules added after it: the document's Field assignment
// applied to s (on a scratch document; this one is not changed), else the built-in force field.
caps::ForceField field_for_extended(caps_doc* d, const caps::System& s) {
  if (!d->field) return default_ff(s);
  std::unique_ptr<caps_doc> t(doc_of_system(s, d, "scratch", "", {}));
  t->field = std::make_unique<FieldState>(*d->field);
  field_run(t.get());
  if (!t->field->complete || !t->field->ff) {
    int untyped = 0;
    for (const auto& x : t->field->types) untyped += x.empty();
    throw std::invalid_argument(d->field->base.name + " does not describe the added molecules (" + std::to_string(untyped) + " atoms untyped, " +
                                std::to_string(t->field->rep.missing.size()) + " parameters missing): choose a force field that covers them in Field, or clear the assignment");
  }
  return *t->field->ff;
}

// Molecules added after a prepared structure (sorbates, adsorbates) get the force field's own preparation (united-atom
// folding, a model's massless sites) on their own: the structure's atoms are prepared already, and folding them again
// would lose their sites' hydrogens.
void prepare_added(caps_doc* d, caps::System& s, size_t nh) {
  if (!d->field || !caps::needs_prepare(d->field->base) || nh >= s.atoms.size()) return;
  if (!s.has_mol) {   // molecule ids from the bonds first: an added site (unbonded) stays in its molecule
    int nm = 0;
    const auto comp = s.molecules(&nm);
    for (size_t i = 0; i < s.atoms.size(); ++i) s.atoms[i].mol = int64_t(comp[i]) + 1;
    s.has_mol = true;
  }
  caps::System sorb;
  sorb.cell = s.cell, sorb.has_mol = true;
  sorb.atoms.assign(s.atoms.begin() + long(nh), s.atoms.end());
  for (const auto& b : s.bonds)
    if (b.i >= nh && b.j >= nh) sorb.bonds.push_back({uint32_t(b.i - nh), uint32_t(b.j - nh), b.order});
  std::string ch = "keep";
  caps::prepare_for_forcefield(sorb, d->field->base, ch);
  s.atoms.resize(nh);
  std::vector<caps::Bond> hb;
  for (const auto& b : s.bonds) if (b.i < nh && b.j < nh) hb.push_back(b);
  s.bonds = std::move(hb);
  for (const auto& a : sorb.atoms) s.atoms.push_back(a);
  for (const auto& b : sorb.bonds) s.bonds.push_back({uint32_t(b.i + nh), uint32_t(b.j + nh), b.order});
  for (size_t i = 0; i < s.atoms.size(); ++i) s.atoms[i].id = int64_t(i) + 1;
}

// the structure's molecules from SMILES, count copies of each, added after its atoms (positions as built)
caps::System with_molecules(caps::System s, const caps::Json& list, int* added) {
  int64_t next_mol = 0;
  for (const auto& a : s.atoms) next_mol = std::max(next_mol, a.mol);
  *added = 0;
  for (const auto& a : list.items()) {
    const std::string smi = a.text("smiles");
    const int count = int(a.num("count", 1));
    if (smi.empty() || count < 1) continue;
    const caps::System m = caps::build_molecule(smi).system;
    for (int k = 0; k < count; ++k) {
      const uint32_t off = uint32_t(s.atoms.size());
      ++next_mol;
      for (auto at : m.atoms) { at.mol = next_mol; at.id = int64_t(s.atoms.size() + 1); s.atoms.push_back(at); }
      for (auto b : m.bonds) { b.i += off; b.j += off; s.bonds.push_back(b); }
      ++*added;
    }
  }
  return s;
}
}  // namespace

// Adsorption locator (ABI 34): adsorbates from SMILES added after the substrate, the Field assignment applied to the
// whole (or the built-in force field), Monte Carlo simulated annealing (adsorption.hpp).
extern "C" int32_t caps_adsorption(caps_doc* d, const char* json, caps_stage_fn progress, void* user, char* out, int32_t cap) {
  caps::Json r = caps::Json::object();
  try {
    if (!d) throw std::invalid_argument("no document");
    const caps::Json j = caps::Json::parse(json && *json ? json : "{}");
    caps::System s = d->traj.frame(d->current);
    const size_t substrate_atoms = s.atoms.size();
    int first = -1, added = 0;
    if (j.has("adsorbates") && j["adsorbates"].is_array()) {
      s = with_molecules(std::move(s), j["adsorbates"], &added);
      if (added) first = int(substrate_atoms);
    }
    if (j.has("first_atom")) first = int(j.num("first_atom", -1));
    if (added) prepare_added(d, s, substrate_atoms);   // the force field's own preparation of the adsorbates (sites, united atom)
    const caps::ForceField ff = added ? field_for_extended(d, s) : [&] { const auto fp = field_for_run(d); return fp ? *fp : default_ff(s); }();
    if (added) caps::apply_rigid_geometry(s, ff, substrate_atoms, s.atoms.size());   // rigid models: the force field's geometry
    caps::AdsorptionOptions o;
    o.first_mobile_atom = first;
    o.cycles = int(j.num("cycles", 3));
    o.steps = int(j.num("steps", 20000));
    o.t_high = j.num("t_high", 1e4);
    o.t_low = j.num("t_low", 100);
    o.cutoff = j.num("cutoff", 12.0);
    o.coulomb = !(j.has("coulomb") && j["coulomb"].kind() == caps::Json::Bool && !j["coulomb"].boolean());
    o.keep = int(j.num("keep", 10));
    o.seed = uint64_t(j.num("seed", 1));
    const std::string side = j.text("region", "cell");   // cell | above (over the substrate's top face)
    if (j.has("z_lo") || j.has("z_hi")) o.z_lo = j.num("z_lo", 0), o.z_hi = j.num("z_hi", 0);
    else if (side == "above" && substrate_atoms > 0 && s.cell.valid()) {
      double top = -1e30;
      for (size_t i = 0; i < substrate_atoms; ++i) top = std::max(top, s.atoms[i].pos[2]);
      o.z_lo = top + 1.0;
      o.z_hi = s.cell.origin[2] + s.cell.c[2];
      if (o.z_hi <= o.z_lo + 2) throw std::invalid_argument("no vacuum above the substrate: add vacuum along z (Surface builder) or search the whole cell");
    }
    bool cancelled = false;
    if (progress)
      o.progress = [&](int cyc, int step, double best) {
        char b[96];
        std::snprintf(b, sizeof b, "cycle %d · step %d · lowest %.2f kcal/mol", cyc + 1, step, best);
        cancelled = progress(b, (cyc + double(step) / o.steps) / o.cycles, user) != 0;
        return !cancelled;
      };
    caps::AdsorptionReport rep;
    caps::locate_adsorption(s, ff, o, &rep);
    if (cancelled) throw std::runtime_error("adsorption locator cancelled");
    push_undo(d, "Adsorption locator");
    // frames: the kept configurations, the lowest last (shown)
    caps::Trajectory t;
    t.topology = s;
    for (auto it = rep.configs.rbegin(); it != rep.configs.rend(); ++it) {
      t.positions.push_back(it->positions), t.cells.push_back(s.cell), t.timesteps.push_back(int64_t(t.timesteps.size()));
    }
    t.topology.notes = rep.notes;
    d->traj = std::move(t);
    d->current = d->traj.frames() - 1;
    refresh(d);
    if (added && d->field) field_run(d);   // the assignment now over the adsorbates too
    caps::KeyValues pr = {{"adsorbates", std::to_string(rep.components.size()) + " components"},
                          {"annealing", std::to_string(o.cycles) + " cycles of " + std::to_string(o.steps) + " steps, " + g6(o.t_high) + " → " + g6(o.t_low) + " K"},
                          {"adsorption energy", g6(rep.adsorption_energy) + " kcal/mol (rigid)"},
                          {"force field", ff_label(d)}};
    prov_step(d, "adsorption.locator", "Adsorption sites by Monte Carlo simulated annealing", std::move(pr), seeded(o.seed), {"kirkpatrick1983"},
              energy_approx(o.cutoff, o.coulomb, false, 1));
    r["ok"] = true;
    r["adsorption_energy"] = rep.adsorption_energy;
    r["adsorbate_substrate"] = rep.adsorbate_substrate;
    r["adsorbate_adsorbate"] = rep.adsorbate_adsorbate;
    caps::Json comps = caps::Json::array();
    for (const auto& c : rep.components) {
      caps::Json x = caps::Json::object();
      x["name"] = c.name, x["molecules"] = double(c.molecules), x["de_dn"] = c.de_dn;
      comps.push_back(std::move(x));
    }
    r["components"] = std::move(comps);
    caps::Json cf = caps::Json::array();
    for (const auto& c : rep.configs) { caps::Json x = caps::Json::object(); x["energy"] = c.energy, x["cycle"] = double(c.cycle); cf.push_back(std::move(x)); }
    r["configs"] = std::move(cf);
    caps::Json he = caps::Json::array(), hc = caps::Json::array();
    for (double x : rep.hist_edges) he.push_back(x);
    for (double x : rep.hist_counts) hc.push_back(x);
    r["hist_edges"] = std::move(he), r["hist_counts"] = std::move(hc);
    r["acceptance"] = rep.acceptance, r["steps"] = double(rep.steps), r["seconds"] = rep.seconds;
    r["z_lo"] = o.z_lo, r["z_hi"] = o.z_hi;
    r["forcefield"] = ff.name;
    caps::Json nt = caps::Json::array();
    for (const auto& x : rep.notes) nt.push_back(x);
    r["notes"] = std::move(nt);
  } catch (const std::exception& e) {
    r = caps::Json::object();
    r["ok"] = false;
    r["error"] = std::string(e.what());
  }
  return report_out(r.dump(0), out, cap);
}

// Sorption (ABI 34): Widom insertion and GCMC of a rigid sorbate in the document's structure held fixed (sorption.hpp).
extern "C" int32_t caps_sorption(caps_doc* d, const char* json, caps_stage_fn progress, void* user, char* out, int32_t cap) {
  caps::Json r = caps::Json::object();
  try {
    if (!d) throw std::invalid_argument("no document");
    const caps::Json j = caps::Json::parse(json && *json ? json : "{}");
    const caps::System host = d->traj.frame(d->current);
    // the sorbates: "sorbate" (one gas) or v63 "mixture": [{"smiles", "fraction"}, …] (an ideal gas mixture), one
    // template each after the host
    std::vector<std::string> smiles;
    std::vector<double> fractions;
    if (j.has("mixture") && j["mixture"].is_array() && !j["mixture"].items().empty()) {
      for (const auto& m : j["mixture"].items()) smiles.push_back(m.text("smiles", "")), fractions.push_back(m.num("fraction", 1.0));
    } else
      smiles.push_back(j.text("sorbate", "C"));
    caps::System s = host;
    std::vector<int> starts;
    for (const auto& smi : smiles) {
      if (smi.empty()) throw std::invalid_argument("give each sorbate as SMILES");
      caps::Json one = caps::Json::array();
      caps::Json sb = caps::Json::object();
      sb["smiles"] = smi;
      sb["count"] = 1.0;
      one.push_back(std::move(sb));
      int added = 0;
      const int first = int(s.atoms.size());
      s = with_molecules(std::move(s), one, &added);
      if (!added) throw std::invalid_argument("could not build the sorbate " + smi);
      starts.push_back(first);
    }
    if (d->field && caps::needs_prepare(d->field->base)) {   // the model's own preparation on the sorbates too
      if (!s.has_mol) {
        int nm = 0;
        const auto comp = s.molecules(&nm);
        for (size_t i = 0; i < s.atoms.size(); ++i) s.atoms[i].mol = int64_t(comp[i]) + 1;
        s.has_mol = true;
      }
      std::vector<int64_t> mols;
      for (int st : starts) mols.push_back(s.atoms[size_t(st)].mol);
      prepare_added(d, s, host.atoms.size());
      for (size_t k = 0; k < starts.size(); ++k)
        for (size_t i = host.atoms.size(); i < s.atoms.size(); ++i)
          if (s.atoms[i].mol == mols[k]) { starts[k] = int(i); break; }
    }
    const caps::ForceField ff = field_for_extended(d, s);
    // rigid models (TraPPE CO2, N2, O2) take their force field's bond lengths and angles
    const int reshaped = caps::apply_rigid_geometry(s, ff, host.atoms.size(), s.atoms.size());
    caps::SorptionOptions o;
    o.template_first_atom = int(host.atoms.size());
    if (smiles.size() > 1) o.species_first_atom = starts, o.mole_fractions = fractions;
    o.temperature = j.num("temperature", 300);
    o.insertions = int(j.num("insertions", 100000));
    if (j.has("pressures_kpa") && j["pressures_kpa"].is_array())
      for (const auto& x : j["pressures_kpa"].items()) if (x.number() > 0) o.pressures_kpa.push_back(x.number());
    o.steps = int(j.num("steps", 200000));
    o.cutoff = j.num("cutoff", 12.0);
    o.coulomb = !(j.has("coulomb") && j["coulomb"].kind() == caps::Json::Bool && !j["coulomb"].boolean());
    o.seed = uint64_t(j.num("seed", 1));
    o.map_grid = int(j.num("map_grid", 0));   // v62: sorbate density maps (C7)
    bool cancelled = false;
    if (progress) o.progress = [&](const std::string& st, double f) { cancelled = progress(st.c_str(), f, user) != 0; return !cancelled; };
    auto rep = caps::sorption(s, ff, o);
    if (cancelled) throw std::runtime_error("sorption cancelled");
    if (reshaped) rep.notes.push_back(std::to_string(reshaped) + " bond lengths and angles of the rigid sorbate set to " + ff.name + "'s own");
    std::string sorbates;
    for (size_t k = 0; k < smiles.size(); ++k)
      sorbates += (k ? " + " : "") + smiles[k] + (smiles.size() > 1 ? " (y " + g6(o.mole_fractions.empty() ? 1.0 : rep.species[k].fraction) + ")" : "");
    caps::KeyValues pr = {{"sorbate", sorbates}, {"temperature", g6(o.temperature) + " K"},
                          {"Widom", std::to_string(o.insertions) + " insertions · S " + g6(rep.solubility) + " cm³(STP)/(cm³ atm)"},
                          {"GCMC", std::to_string(o.pressures_kpa.size()) + " pressures × " + std::to_string(o.steps) + " steps"},
                          {"force field", ff.name}};
    prov_step(d, "sorption.widom_gcmc", "Sorption by test-particle insertion and grand-canonical Monte Carlo", std::move(pr), seeded(o.seed),
              {"widom1963", "frenkel2002"}, energy_approx(o.cutoff, o.coulomb, false, 1));
    r["ok"] = true;
    r["widom_w"] = rep.widom_w, r["widom_error"] = rep.widom_error, r["mu_ex"] = rep.mu_ex;
    r["henry_mol_kg_kpa"] = rep.henry_mol_kg_kpa, r["solubility"] = rep.solubility;
    r["host_mass"] = rep.host_mass, r["volume"] = rep.volume, r["forcefield"] = ff.name;
    caps::Json iso = caps::Json::array();
    for (const auto& p : rep.isotherm) {
      caps::Json x = caps::Json::object();
      x["pressure_kpa"] = p.pressure_kpa, x["loading"] = p.loading, x["loading_error"] = p.loading_error, x["mol_per_kg"] = p.mol_per_kg;
      x["cm3stp_per_cm3"] = p.cm3stp_per_cm3, x["heat"] = p.heat, x["acceptance_insert"] = p.acceptance_insert, x["acceptance_delete"] = p.acceptance_delete;
      auto arr = [](const std::vector<double>& v) { caps::Json a = caps::Json::array(); for (double q : v) a.push_back(q); return a; };
      x["species_loading"] = arr(p.species_loading), x["species_error"] = arr(p.species_error), x["species_mol_per_kg"] = arr(p.species_mol_per_kg);
      x["species_heat"] = arr(p.species_heat), x["selectivity"] = arr(p.selectivity);
      if (!p.density.empty()) {   // the density map: grid points along a, b, c and the values (per Å³), a slowest
        caps::Json m = caps::Json::object(), v = caps::Json::array();
        for (float q : p.density) v.push_back(double(q));
        m["grid"] = double(p.grid);
        m["density"] = std::move(v);
        x["map"] = std::move(m);
      }
      iso.push_back(std::move(x));
    }
    r["isotherm"] = std::move(iso);
    caps::Json sp = caps::Json::array();
    for (size_t k = 0; k < rep.species.size(); ++k) {
      const auto& w = rep.species[k];
      caps::Json x = caps::Json::object();
      x["smiles"] = k < smiles.size() ? smiles[k] : std::string();
      x["fraction"] = w.fraction, x["widom_w"] = w.widom_w, x["widom_error"] = w.widom_error, x["mu_ex"] = w.mu_ex;
      x["henry_mol_kg_kpa"] = w.henry_mol_kg_kpa, x["solubility"] = w.solubility;
      sp.push_back(std::move(x));
    }
    r["species"] = std::move(sp);
    caps::Json nt = caps::Json::array();
    for (const auto& x : rep.notes) nt.push_back(x);
    r["notes"] = std::move(nt);
  } catch (const std::exception& e) {
    r = caps::Json::object();
    r["ok"] = false;
    r["error"] = std::string(e.what());
  }
  return report_out(r.dump(0), out, cap);
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
    if (o->etol > 0) r.etol = o->etol;
    if (o->pressure_tol > 0) r.pressure_tol = o->pressure_tol;
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
    r.fixed = fixed_mask(d, s);
    require_manybody_held(d, r.fixed);
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
    if (o->box_anisotropic && m.barostat != caps::Barostat::None) {   // per axis: Berendsen scales each axis from P_kk
      m.anisotropic = true;
      m.barostat = caps::Barostat::Berendsen;
      const int ax = o->box_axes & 7 ? o->box_axes & 7 : 7;
      for (int k = 0; k < 3; ++k) m.couple_axis[k] = (ax >> k) & 1;
    }
    if (o->tau_p > 0) m.tau_p = o->tau_p;
    m.new_velocities = o->new_velocities != 0;
    m.efield = {o->efield[0], o->efield[1], o->efield[2]};
    m.full_shape = o->full_shape != 0 && m.anisotropic;
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
    m.constraint_algorithm = o->constraint_algorithm == 1 ? caps::ConstraintAlgorithm::Lincs : caps::ConstraintAlgorithm::Shake;
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
    m.fixed = fixed_mask(d, s);
    require_manybody_held(d, m.fixed);
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
      if (m.constraints != caps::ConstraintMode::None) c.push_back(m.constraint_algorithm == caps::ConstraintAlgorithm::Lincs ? "hess1997" : "ryckaert1977"), c.push_back("andersen1983");
      if (m.thermostat == caps::Thermostat::Bussi) c.push_back("bussi2007");
      if (npt && m.barostat == caps::Barostat::CRescale) c.push_back("bernetti2020");
      if (npt && m.barostat == caps::Barostat::Berendsen) c.push_back("berendsen1984");
      elec_cites(c, m.energy.coulomb);
      caps::KeyValues pr = {{"length", g6(m.dt * double(m.steps) / 1000.0) + " ps · " + std::to_string(m.steps) + " steps of " + g6(m.dt) + " fs"},
                            {"temperature", g6(m.temperature) + " K"}, {"thermostat", std::string(caps::to_string(m.thermostat)) + (nvt ? " · τ " + g6(m.tau_t) + " fs" : "")}};
      if (npt) pr.push_back({"barostat", std::string(caps::to_string(m.barostat)) + " · " + g6(m.pressure) + " atm · τ " + g6(m.tau_p) + " fs"});
      if (m.constraints != caps::ConstraintMode::None) {
        pr.push_back({"constraints", caps::to_string(m.constraints)});
        pr.push_back({"constraint solver", m.constraint_algorithm == caps::ConstraintAlgorithm::Lincs ? "LINCS" : "SHAKE"});
      }
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

int32_t caps_trajectory_columns(caps_doc* d, char* json, int32_t cap) {
  return guard([&] {
    caps::Json j = caps::Json::object();
    caps::Json cols = caps::Json::array();
    for (const auto& [name, frames] : d->traj.columns) {
      caps::Json c = caps::Json::object();
      c["name"] = name;
      double lo = 1e300, hi = -1e300;
      if (d->current < frames.size())
        for (float v : frames[d->current]) if (std::isfinite(v)) lo = std::min(lo, double(v)), hi = std::max(hi, double(v));
      c["min"] = lo <= hi ? lo : 0.0;
      c["max"] = lo <= hi ? hi : 0.0;
      cols.push_back(c);
    }
    j["columns"] = cols;
    j["velocities"] = !d->traj.velocities.empty();
    return report_out(j.dump(), json, cap);
  });
}

int32_t caps_save_trajectory(caps_doc* d, const char* path) {
  return guard([&] {
    caps::write_trajectory(d->traj, path ? path : "", 1.0);
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
    e.md.constraint_algorithm = o->constraint_algorithm == 1 ? caps::ConstraintAlgorithm::Lincs : caps::ConstraintAlgorithm::Shake;
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
    if (o->internal_target && o->internal_target_n > 0) e.internal_target.assign(o->internal_target, o->internal_target + o->internal_target_n);
    if (progress)
      e.progress = [&](int st, int n, const std::string& label, const caps::ThermoRow& r) {
        caps_thermo t{r.step, r.time_ps, r.temperature, r.potential, r.kinetic, r.total, r.conserved, r.pressure, r.volume, r.density};
        return progress(st, n, label.c_str(), &t, user) == 0;
      };
    caps::System s = d->traj.frame(d->current);
    if (d->current + 1 != d->traj.frames()) s.velocities.clear();
    if (!s.unwrapped) caps::make_molecules_whole(s);
    e.md.fixed = fixed_mask(d, s);
    require_manybody_held(d, e.md.fixed);
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
      if (e.md.constraints != caps::ConstraintMode::None) c.push_back(e.md.constraint_algorithm == caps::ConstraintAlgorithm::Lincs ? "hess1997" : "ryckaert1977"), c.push_back("andersen1983");
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
    {
      caps::Json a = caps::Json::array();
      for (const auto& it : rep.items) {
        caps::Json o = caps::Json::object();
        o["name"] = it.name, o["molecules"] = it.molecules, o["forcefield"] = it.forcefield;
        a.push_back(o);
      }
      d->pack_items = a.dump(0);
    }
    refresh(d);
    prov_step(d, "pack.lbfgs", std::to_string(rep.molecules) + " molecules packed without overlaps",
              {{"molecules", std::to_string(rep.molecules)}, {"atoms", std::to_string(rep.atoms)}, {"closest contact", g6(rep.dmin) + " Å"},
               {"input sha256", caps::sha256_hex(text ? text : "").substr(0, 12)},
               // the input itself (molecules, counts, regions, box): the packing can be set up again from it
               {"input", std::string(text ? text : "").substr(0, 20000)}},
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

}  // extern "C": C++ helpers

namespace {

// The Field assignment of d re-run on another structure (a reaction's product): the same force field, typing rules,
// imported parameters and charge choice; types set by hand are dropped (the atoms are renumbered). Throws with what is
// untyped or missing when the force field cannot describe the product.
std::unique_ptr<FieldState> field_on(const caps_doc* d, const caps::System& s) {
  caps_doc t;
  t.traj.topology = s;
  t.traj.topology.source_format.clear();
  std::vector<caps::Vec3> p;
  for (const auto& a : s.atoms) p.push_back(a.pos);
  t.traj.positions.push_back(std::move(p));
  t.traj.cells.push_back(s.cell);
  t.traj.timesteps.push_back(0);
  t.frame = s;
  t.field = std::make_unique<FieldState>(*d->field);
  t.field->overrides.clear();
  t.field->file_types.clear();
  t.field->file_charges.clear();
  field_run(&t);
  return std::move(t.field);
}

std::shared_ptr<const caps::ForceField> field_for_product(const caps_doc* d, const caps::System& s) {
  const auto F = field_on(d, s);
  if (!F->complete) {
    std::map<std::string, int> untyped;
    for (size_t i = 0; i < F->types.size() && i < s.atoms.size(); ++i)
      if (F->types[i].empty()) ++untyped[std::string(caps::element(s.atoms[i].element).symbol)];
    std::string why;
    for (const auto& [el, k] : untyped) why += (why.empty() ? "" : ", ") + std::to_string(k) + " " + el;
    if (!why.empty()) why = "untyped atoms: " + why;
    for (size_t k = 0; k < F->rep.missing.size() && k < 4; ++k) why += (why.empty() ? "missing " : "; missing ") + F->rep.missing[k];
    if (F->rep.missing.size() > 4) why += " (" + std::to_string(F->rep.missing.size()) + " terms missing in all)";
    throw caps::FieldError(why.empty() ? "the assignment is incomplete" : why);
  }
  return F->ff;
}

std::vector<double> number_list(const char* text) {
  std::vector<double> v;
  std::string t = text ? text : "";
  for (auto& c : t) if (c == ',' || c == ';') c = ' ';
  std::istringstream is(t);
  for (std::string w; is >> w;) {
    try { v.push_back(std::stod(w)); } catch (const std::exception&) { throw std::invalid_argument("weights: '" + w + "' is not a number"); }
  }
  return v;
}

}  // namespace

extern "C" {

int32_t caps_react_summary(caps_doc* d, char* json, int32_t cap) {
  return guard([&] {
    if (d->react_json.empty()) throw std::runtime_error("no reaction run on this structure yet");
    return report_out(d->react_json, json, cap);
  });
}

int32_t caps_bond_react_export(caps_doc* d, const char* templates, const char* dir, const char* options, char* report, int32_t cap) {
  return guard([&] {
    if (!d->field || !d->field->ff) throw std::runtime_error("assign a force field first (Force field step): the templates are typed with it");
    if (!d->field->complete) throw std::runtime_error("the force field is incomplete for this structure: complete it in the Force field step");
    const caps::Json o = caps::Json::parse(options && *options ? options : "{}");
    auto flag = [&](const char* k, bool def) { return o.has(k) && o[k].kind() == caps::Json::Bool ? o[k].boolean() : def; };
    const std::string stem = o.text("stem", "react");
    if (stem.empty() || stem.find('/') != std::string::npos || stem.find('\\') != std::string::npos) throw std::runtime_error("the file stem must be a plain name");
    caps::BondReactOptions b;
    b.radius = std::clamp(int(o.num("radius", 3)), 1, 8);
    b.max_variants = std::clamp(int(o.num("variants", 6)), 1, 50);
    b.keep_byproducts = flag("keep_byproducts", false);
    b.between_chains = flag("between_chains", false);
    if (o.has("weights") && o["weights"].is_array())
      for (const auto& w : o["weights"].items()) b.weights.push_back(w.number());
    b.nevery = std::max(1, int(o.num("nevery", 100)));
    b.rmax = std::max(0.0, o.num("rmax", 3.5));
    b.temperature = o.num("temperature", 300);
    b.steps = std::max<int64_t>(1, int64_t(o.num("steps", 100000)));
    b.seed = uint64_t(std::max(1.0, o.num("seed", 12345)));
    b.energy = elec(caps::EnergyOptions{});
    caps::System s = d->frame;
    const auto t = caps::parse_templates(templates ? templates : "");
    const auto r = caps::write_bond_react(s, t, [d](const caps::System& x) { return field_for_product(d, x); }, dir ? dir : ".", stem, b);
    caps::Json j = caps::Json::object();
    caps::Json files = caps::Json::array(), notes = caps::Json::array(), vars = caps::Json::array();
    for (const auto& f : r.files) files.push_back(f);
    for (const auto& n : r.notes) notes.push_back(n);
    for (const auto& v : r.variants) {
      caps::Json x = caps::Json::object();
      x["reaction"] = v.reaction;
      x["name"] = v.name;
      x["sites"] = double(v.sites);
      x["pre_atoms"] = double(v.pre_atoms);
      x["edge"] = double(v.edge);
      x["deleted"] = double(v.deleted);
      vars.push_back(x);
    }
    j["files"] = files;
    j["notes"] = notes;
    j["variants"] = vars;
    j["candidates"] = double(r.candidates);
    j["covered"] = double(r.covered);
    return report_out(j.dump(), report, cap);
  });
}

int32_t caps_reaction_library(const char* path, char* json, int32_t cap) {
  return guard([&] {
    std::ifstream in(path ? path : "");
    if (!in) throw std::runtime_error(std::string("cannot read the reaction library ") + (path ? path : ""));
    std::stringstream ss;
    ss << in.rdbuf();
    caps::Json lib = caps::Json::parse(ss.str());
    caps::Json out = caps::Json::array();
    for (auto e : lib["reactions"].items()) {
      std::vector<std::string> rs, ps;
      for (const auto& r : e["reactants"].items()) rs.push_back(r.text("smiles", ""));
      for (const auto& p : e["products"].items()) ps.push_back(p.text("smiles", ""));
      std::vector<std::string> notes;
      try {
        e["template"] = caps::template_text(caps::template_from_scheme(rs, ps, e.text("id", "reaction"), &notes));
      } catch (const std::exception& x) {
        e["error"] = std::string(x.what());
      }
      caps::Json n = caps::Json::array();
      for (const auto& x : notes) n.push_back(x);
      e["notes"] = n;
      out.push_back(e);
    }
    caps::Json j = caps::Json::object();
    j["reactions"] = out;
    return report_out(j.dump(), json, cap);
  });
}

int32_t caps_bond_react_import(const char* pre, const char* post, const char* map, const char* masses_from, const char* name, double capture,
                               char* json, int32_t cap) {
  return guard([&] {
    std::vector<std::string> notes;
    const auto t = caps::read_bond_react(pre ? pre : "", post ? post : "", map ? map : "", masses_from ? masses_from : "", name ? name : "", capture, &notes);
    caps::Json j = caps::Json::object();
    j["text"] = caps::template_text(t);
    caps::Json n = caps::Json::array();
    for (const auto& x : notes) n.push_back(x);
    j["notes"] = n;
    return report_out(j.dump(), json, cap);
  });
}

int32_t caps_react_sites(caps_doc* d, const char* templates_text, char* json, int32_t cap) {
  caps::Json r = caps::Json::object();
  guard([&] {
    const auto ts = caps::parse_templates(templates_text ? templates_text : "");
    // after a run: its chains (joined chains still counted one by one)
    const bool runs = d->react_chains.size() == d->frame.atoms.size() && d->react_bonds == d->frame.bonds.size();
    const auto cs = caps::chain_sites(d->frame, ts, runs ? d->react_chains : std::vector<int64_t>{});
    caps::Json a = caps::Json::array();
    int total = 0, units = 0;
    double mass = 0;
    for (const auto& c : cs) {
      caps::Json x = caps::Json::object();
      x["chain"] = double(c.chain), x["sites"] = double(c.sites), x["units"] = double(c.units), x["atoms"] = double(c.atoms), x["mass"] = c.mass;
      a.push_back(x);
      total += c.sites, units += c.units, mass += c.mass;
    }
    r["chains"] = a;
    r["chains_n"] = double(cs.size());
    r["total_sites"] = double(total);
    r["units"] = double(units);
    r["mass"] = mass;
    r["repeat_unit_mass"] = units > 0 ? mass / units : 0.0;
    return 0;
  });
  return report_out(r.dump(0), json, cap);
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
    r.between_chains = o->between_chains != 0;
    if (d->react_chains.size() == d->traj.topology.atoms.size() && d->react_bonds == d->traj.topology.bonds.size()) r.chains = d->react_chains;
    r.keep_byproducts = o->keep_byproducts != 0;
    r.selection = std::clamp(o->selection, 0, 1);
    r.weights = number_list(o->weights);
    if (!r.weights.empty() && r.weights.size() != r.templates.size())
      throw std::invalid_argument("weights: " + std::to_string(r.weights.size()) + " given for " + std::to_string(r.templates.size()) + " templates");
    for (double w : r.weights) if (w <= 0) throw std::invalid_argument("weights must be > 0");
    r.auto_capture = o->auto_capture != 0;
    if (o->capture_max > 0) r.capture_max = o->capture_max;
    if (o->capture_step > 0) r.capture_step = o->capture_step;
    if (o->target_kind < 0 || o->target_kind > 5) throw std::invalid_argument("target_kind: 0 conversion … 4 Mc, 5 degree of crosslinking %");
    r.target = caps::ReactTarget(o->target_kind);
    r.target_value = o->target_value;
    r.sites_per_chain = std::max(0, o->sites_per_chain);
    // the user's force field for every state of the network (a complete assignment is needed from the start)
    const bool use_field = o->field_mode == 0 && d->field;
    if (use_field) {
      if (!d->field->complete) throw caps::FieldError("the assigned force field is incomplete for this structure: complete it in the Force field step, or run with the built-in default");
      r.retype = [d](const caps::System& x) { return field_for_product(d, x); };
      r.field_name = ff_label(d);
    }
    if (progress)
      r.progress = [&](const caps::CycleRow& c) {
        caps_react_cycle row{c.cycle, c.reactions, c.total, c.clusters.clusters, c.atoms, c.conversion, c.clusters.largest_fraction,
                             c.clusters.reduced_mw, c.energy, c.max_force, c.crosslinks, c.capture, c.target, c.density, c.degree};
        return progress(&row, user) == 0;
      };
    if (d->live_fn) {
      auto fn = d->live_fn;
      auto lu = d->live_user;
      r.live = [fn, lu](const caps::System& x, const caps::CycleRow& c, const std::vector<int64_t>& chain, const std::vector<char>& linked) {
        auto* sd = new caps_doc;
        sd->traj.topology = x;
        int64_t top = 0;
        for (int64_t v : chain) top = std::max(top, v);
        for (size_t i = 0; i < x.atoms.size() && i < chain.size(); ++i)
          sd->traj.topology.atoms[i].mol = chain[i] > 0 ? chain[i] : top - chain[i];   // byproducts after the chains
        sd->traj.topology.has_mol = true;
        std::vector<caps::Vec3> p;
        for (const auto& a : x.atoms) p.push_back(a.pos);
        sd->traj.positions.push_back(std::move(p));
        sd->traj.cells.push_back(x.cell);
        sd->traj.timesteps.push_back(c.cycle);
        sd->wrap = x.cell.valid();
        refresh(sd);
        sd->selection.assign(x.atoms.size(), 0);
        for (size_t i = 0; i < x.atoms.size() && i < linked.size(); ++i) sd->selection[i] = linked[i];
        caps::Json j = caps::Json::object();
        j["cycle"] = double(c.cycle), j["reactions"] = double(c.total), j["crosslinks"] = double(c.crosslinks), j["target"] = double(c.target);
        j["density"] = c.density, j["degree"] = c.degree, j["conversion"] = c.conversion, j["atoms"] = double(x.atoms.size());
        fn(sd, j.dump(0).c_str(), lu);
      };
    }
    caps::System s = d->traj.frame(d->current);
    if (!s.unwrapped) caps::make_molecules_whole(s);
    // the atom count changes when atoms leave, so the record is one document per state: keep the start and the end
    std::vector<caps::System> frames{s};
    r.frame = [&](const caps::System& x, int) {
      frames.push_back(x);
      if (x.cell.valid()) caps::make_molecules_whole(frames.back());   // bonds new this cycle may cross the cell: whole again
    };
    caps::ReactReport rep;
    caps::react(s, r, &rep);
    {
      std::string names;
      for (const auto& t : r.templates) names += (names.empty() ? "" : ", ") + t.name;
      const double conv = rep.cycles.empty() ? 0 : rep.cycles.back().conversion;
      prov_step(d, "react.templates", std::to_string(rep.reactions) + " reactions · conversion " + g6(conv),
                {{"templates", names}, {"target conversion", g6(r.target_conversion)}, {"cycles", std::to_string(rep.cycles.size())},
                 {"relax between cycles", r.relax ? "yes" : "no"}, {"MD between cycles", g6(r.md_ps) + " ps"},
                 // the reaction templates as written: the reaction can be set up again from them
                 {"input", std::string(templates ? templates : "").substr(0, 20000)}},
                seeded(r.seed), {"matsumoto1998"},
                rep.failed_cycle > 0 ? caps::KeyValues{{"Stopped", "cycle " + std::to_string(rep.failed_cycle) + " failed (" + rep.failure + "); the structure after cycle " +
                                                                   std::to_string(rep.failed_cycle - 1) + " is kept"}}
                                     : caps::KeyValues{});
    }
    caps::Trajectory out;
    // the last frame as the run left it (made whole after the last cycle)
    if (!frames.empty() && frames.back().atoms.size() == s.atoms.size()) frames.back() = s;
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
    // the Field assignment follows the new topology: re-run on the product (types, charges and parameters of the new
    // bonds); without one, or with the built-in default chosen, there is none
    std::unique_ptr<FieldState> kept = use_field ? std::make_unique<FieldState>(*d->field) : nullptr;
    d->field.reset();
    d->traj = std::move(out);
    d->current = d->traj.frames() - 1;
    refresh(d);
    d->react_chains = rep.chains_after;
    d->react_bonds = d->traj.topology.bonds.size();
    std::string after = "none (assign one in the Force field step)";
    if (kept) {
      kept->overrides.clear();
      kept->file_types.clear();
      kept->file_charges.clear();
      d->field = std::move(kept);
      try {
        field_run(d);
        after = ff_label(d) + (d->field->complete ? " · complete" : " · incomplete: see the Force field step");
      } catch (const std::exception& e) {
        after = std::string("could not re-assign: ") + e.what();
        d->field.reset();
      }
      rep.notes.push_back("force field after the run: " + after);
      refresh(d);
    }
    {
      caps::Json j;
      j["chains"] = double(rep.chains);
      j["crosslinks"] = double(rep.crosslinks);
      j["intrachain"] = double(rep.intrachain);
      j["byproducts"] = double(rep.byproducts);
      j["target"] = double(rep.target_crosslinks);
      j["volume"] = rep.volume;
      j["chain_mass"] = rep.chain_mass;
      j["density"] = rep.density;
      j["per_chain"] = rep.per_chain;
      j["mc"] = rep.mc;
      j["monomers"] = double(rep.monomers);
      j["degree"] = rep.degree;
      j["reactions"] = double(rep.reactions);
      j["initial_sites"] = double(rep.initial_sites);
      j["conversion"] = rep.cycles.empty() ? 0.0 : rep.cycles.back().conversion;
      j["field"] = rep.field;
      j["field_after"] = after;
      caps::Json links = caps::Json::array();
      for (const auto& lr : rep.links) {
        caps::Json x = caps::Json::object();
        x["cycle"] = double(lr.cycle), x["reaction"] = lr.reaction;
        x["chain_a"] = double(lr.chain_a), x["unit_a"] = double(lr.unit_a), x["chain_b"] = double(lr.chain_b), x["unit_b"] = double(lr.unit_b);
        x["via"] = double(lr.via), x["via_name"] = lr.via_name;
        links.push_back(x);
      }
      j["links"] = links;
      caps::Json notes = caps::Json::array();
      for (const auto& n : rep.notes) notes.push_back(n);
      j["notes"] = notes;
      d->react_json = j.dump();
    }
    if (report && cap > 0) {
      std::string t;
      for (const auto& n : rep.notes) t += n + "\n";
      std::strncpy(report, t.c_str(), size_t(cap) - 1);
      report[cap - 1] = 0;
    }
    if (rep.failed_cycle > 0) {   // the completed cycles are kept; the caller hears which one failed
      g_error = "cycle " + std::to_string(rep.failed_cycle) + " failed: " + rep.failure;
      return 2;
    }
    return 0;
  });
}

namespace { void push_undo(caps_doc* d, const std::string& what); }   // below, with the edits

int32_t caps_field_assign(caps_doc* d, const char* ff_path, const char* rules_path, int32_t charges) {
  return guard([&] {
    if (ff_path && std::string(ff_path) == "file") {   // back to the force field the file carried
      if (!d->traj.topology.forcefield) throw caps::FFError("this structure's file carried no force field (an AMBER prmtop does)");
      if (d->field) {   // the file's own types and charges, as it was opened
        for (size_t i = 0; i < d->traj.topology.atoms.size() && i < d->field->file_types.size(); ++i)
          d->traj.topology.atoms[i].type = d->field->file_types[i].first, d->traj.topology.atoms[i].name = d->field->file_types[i].second;
        if (!d->field->file_type_table.empty()) d->traj.topology.types = d->field->file_type_table;
      }
      install_file_field(d);
      if (!d->field || d->field->model.empty()) throw caps::FFError("the structure no longer holds the atoms its topology file describes");
      prov_step(d, "field.assign", d->traj.topology.forcefield->name, {{"force field", "the topology file's"}, {"charges", "from the topology file"}}, "", {}, {});
      return 0;
    }
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
      if (ch == "keep" && F->auto_charges && F->base.united_atom) F->ua_summed = true;   // still automatic: types' own charges where they carry them
      else if (ch == "keep") F->charges = "keep", F->auto_charges = false;
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

int32_t caps_field_file_available(const caps_doc* d) {
  if (!d) return 0;
  const auto& t = d->traj.topology;
  return t.forcefield && t.forcefield->charge.size() == t.atoms.size() ? 1 : 0;
}

int32_t caps_field_assign_groups(caps_doc* d, const char* json) {
  return guard([&] {
    auto F = std::make_unique<FieldState>();
    F->groups = json && *json ? json : "{}";
    caps::Json::parse(F->groups);   // well-formed before anything changes
    if (d->field) {
      F->file_types = d->field->file_types, F->file_type_table = d->field->file_type_table;
      F->file_charges = d->field->file_charges, F->file_has_charges = d->field->file_has_charges;
    } else {
      for (const auto& a : d->traj.topology.atoms) F->file_types.push_back({a.type, a.name}), F->file_charges.push_back(a.charge);
      F->file_type_table = d->traj.topology.types;
      F->file_has_charges = d->traj.topology.has_charges;
    }
    auto old = std::move(d->field);
    d->field = std::move(F);
    try {
      field_run(d);
    } catch (...) {
      d->field = std::move(old);   // the previous assignment stays
      throw;
    }
    const caps::Json G = caps::Json::parse(d->field->groups);
    caps::KeyValues pr;
    for (const auto& g : G["groups"].items())
      pr.push_back({g.text("name", "group"), (g.has("potential") && g["potential"].is_object()
                                                   ? g["potential"].text("style", "") + " potential " + std::filesystem::path(g["potential"].text("file", "")).filename().string()
                                                   : std::filesystem::path(g.text("forcefield", "")).filename().string()) + " · molecules " + g.text("molecules", "")});
    pr.push_back({"between groups", "ε " + G.text("eps_rule", "auto") + ", σ " + G.text("sigma_rule", "auto")});
    if (G.text("scaling14", "own") == "first") pr.push_back({"1-4 scaling", "the first group's for all (asked)"});
    if (G.text("cross96", "refuse") == "rmin") pr.push_back({"9-6 sites in cross pairs", "12-6 with the same ε and r_min (asked)"});
    if (G.text("cross96", "refuse") == "area") pr.push_back({"9-6 sites in cross pairs", "12-6 with the same r_min, ε for the same ∫ U r² dr to the cut-off (asked)"});
    prov_step(d, "field.assign.groups", "force fields by group: " + d->field->ff->name, std::move(pr), "", {},
              {{"Cross interactions", "Lennard-Jones by the stated mixing rule between the groups' parameters"}});
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
    if (!d->field->groups.empty())
      throw caps::FFError("the force field is assigned by group: overrides, entered and imported parameters apply to one force field. Change the group's force field, or give its pairs in the group settings");
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

int32_t caps_field_groups_by_example(caps_doc* d, caps_doc* example, const char* types_json, const char* ff_path, const char* groups_json, char* report,
                                     int32_t cap) {
  return guard([&] {
    const bool given = groups_json && *groups_json;
    if (!given && (!d->field || d->field->groups.empty() || d->field->groups == "{}")) throw caps::FFError("assign the force fields by group first");
    const caps::Json tj = caps::Json::parse(types_json ? types_json : "[]");
    std::vector<std::string> types;
    if (example) {
      for (const auto& x : tj.items()) types.push_back(x.is_string() ? x.str() : "");
      if (types.size() != example->frame.atoms.size())
        throw caps::FFError("one type per example atom (" + std::to_string(example->frame.atoms.size()) + "), " + std::to_string(types.size()) + " given");
    }
    const std::string groups = given ? std::string(groups_json) : d->field->groups;
    auto keep_ex = d->hand_example;
    auto keep_types = d->hand_types;
    auto keep_ff = d->hand_ff;
    if (example) d->hand_example = std::make_shared<caps::System>(example->frame), d->hand_types = types, d->hand_ff = ff_path ? ff_path : "";
    else d->hand_example.reset(), d->hand_types.clear(), d->hand_ff.clear();   // forget what was taught
    d->hand_report = caps::Json::object();
    const int32_t rc = caps_field_assign_groups(d, groups.c_str());
    if (rc < 0) {   // nothing changed: the previous teaching stays
      const std::string err = g_error;
      d->hand_example = keep_ex, d->hand_types = keep_types, d->hand_ff = keep_ff;
      throw caps::FFError(err);
    }
    if (example && !d->hand_report.has("groups")) throw caps::FFError("no group is assigned with this force field: choose it for the polymer's group");
    report_out(d->hand_report.dump(), report, cap);
    return rc;
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

int32_t caps_field_set_mixing(caps_doc* d, const char* rule) {
  return guard([&] {
    if (!d->field) throw caps::FFError("assign a force field first");
    if (!d->field->groups.empty()) throw caps::FFError("the force field is assigned by group: its cross terms are set in the group settings");
    const std::string r = rule ? rule : "";
    if (!r.empty() && r != "arithmetic" && r != "geometric" && r != "sixthpower")
      throw caps::FFError("mixing rule: arithmetic (Lorentz–Berthelot), geometric or sixthpower (Waldman–Hagler); empty for the force field's own");
    d->field->mixing = r;
    field_run(d);
    refresh(d);
    return d->field->complete ? 0 : 1;
  });
}

int32_t caps_field_add_rule(caps_doc* d, const char* kind, const char* types, const char* style, const char* params) {
  return guard([&] {
    if (!d->field) throw caps::FFError("assign a force field first");
    if (!d->field->groups.empty())
      throw caps::FFError("the force field is assigned by group: overrides, entered and imported parameters apply to one force field. Change the group's force field, or give its pairs in the group settings");
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
    if (!d->field->groups.empty())
      throw caps::FFError("the force field is assigned by group: overrides, entered and imported parameters apply to one force field. Change the group's force field, or give its pairs in the group settings");
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
    out->cpu_only = 0;   // everything is drawn on the GPU (surfaces as triangles, the colour-vision preview per pixel)
    if (const double* m = caps::vision_matrix(caps::Vision(d->vision)); m && d->vision) {
      for (int k = 0; k < 9; ++k) d->scene_vision[size_t(k)] = float(m[k]);
      out->vision_matrix = d->scene_vision.data();
      out->vision_severity = d->vision_severity;
    }
    out->n_triangles = int32_t(sc.tri_rgb.size() / 3);
    out->tri_xyz = sc.tri_xyz.data(), out->tri_normal = sc.tri_normal.data(), out->tri_rgb = sc.tri_rgb.data();
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

int32_t caps_pick_at(caps_doc* d, const caps_camera* cam, const caps_render_opts* opt, int32_t x, int32_t y) {
  int32_t out = -1;
  guard([&] {
    auto ro = opts_of(d, opt);
    caps::System imaged;
    const caps::System& sys = view_system(d, ro, imaged);
    int k = caps::Renderer::pick_ray(sys, cam_of(cam), ro, x + 0.5, y + 0.5);   // the pixel's centre
    if (k >= 0 && !d->pstate && !d->frame.atoms.empty()) k = int(size_t(k) % d->frame.atoms.size());   // an image atom picks its original
    if (k >= 0 && d->pstate) k = size_t(k) < d->pstate->origin.size() ? d->pstate->origin[size_t(k)] : -1;
    out = k;
    return 0;
  });
  return out;
}

int32_t caps_set_atom_state(caps_doc* d, const int32_t* atoms, int32_t n, int32_t state) {
  int32_t changed = 0;
  guard([&] {
    const size_t na = d->frame.atoms.size();
    const uint8_t v = uint8_t(std::clamp(state, 0, 2));
    if (d->atom_state.size() != na) d->atom_state.assign(na, 0);
    auto set = [&](size_t i) { const uint8_t now = uint8_t((d->atom_state[i] & 4) | v); if (d->atom_state[i] != now) d->atom_state[i] = now, ++changed; };
    if (!atoms) for (size_t i = 0; i < na; ++i) set(i);
    else for (int32_t k = 0; k < n; ++k) if (atoms[k] >= 0 && size_t(atoms[k]) < na) set(size_t(atoms[k]));
    return 0;
  });
  return changed;
}

int32_t caps_set_atom_lock(caps_doc* d, const int32_t* atoms, int32_t n, int32_t locked) {
  int32_t changed = 0;
  guard([&] {
    const size_t na = d->frame.atoms.size();
    if (d->atom_state.size() != na) d->atom_state.assign(na, 0);
    auto set = [&](size_t i) {
      const uint8_t now = uint8_t(locked ? d->atom_state[i] | 4 : d->atom_state[i] & 3);
      if (d->atom_state[i] != now) d->atom_state[i] = now, ++changed;
    };
    if (!atoms) for (size_t i = 0; i < na; ++i) set(i);
    else for (int32_t k = 0; k < n; ++k) if (atoms[k] >= 0 && size_t(atoms[k]) < na) set(size_t(atoms[k]));
    return 0;
  });
  return changed;
}

// Layers (design/boards/Layers): the structure's molecules grouped by kind (residue name and formula), each with its atom
// count, a profile of its atoms along z (8 bins across the cell's c axis, else the structure's extent, peak 1), its view
// state (shown, ghost, hidden or mixed), whether it is locked and how many of its atoms are selected.
int32_t caps_layers(caps_doc* d, char* json, int32_t cap) {
  caps::Json r = caps::Json::object();
  guard([&] {
    const auto& S = d->frame;
    const size_t n = S.atoms.size();
    const auto ids = molecule_ids(S);
    const bool st = d->atom_state.size() == n, sel = d->selection.size() == n;
    // z of each atom as a fraction of the cell (or the extent)
    double lo = 1e300, hi = -1e300;
    if (!S.cell.valid()) for (const auto& a : S.atoms) lo = std::min(lo, a.pos[2]), hi = std::max(hi, a.pos[2]);
    auto zfrac = [&](const caps::Atom& a) {
      if (S.cell.valid()) { double f = S.cell.to_fractional(a.pos)[2]; return f - std::floor(f); }
      return hi > lo ? (a.pos[2] - lo) / (hi - lo) : 0.5;
    };
    struct Mol { int64_t id = 0; size_t atoms = 0, sel = 0; std::array<double, 8> z{}; size_t shown = 0, ghost = 0, hidden = 0, locked = 0; std::vector<std::pair<int, int>> el; std::string res; bool many_res = false; int64_t resid0 = 0; bool many_units = false; };
    std::vector<Mol> list;
    std::unordered_map<int64_t, size_t> slot;
    for (size_t i = 0; i < n; ++i) {
      auto [it, fresh] = slot.try_emplace(ids[i], list.size());
      if (fresh) list.emplace_back(), list.back().id = ids[i], list.back().res = S.atoms[i].resname;
      auto& m = list[it->second];
      ++m.atoms;
      if (!m.many_res && S.atoms[i].resname != m.res) m.many_res = true;   // a chain of several residues (a copolymer's units)
      if (m.atoms == 1) m.resid0 = S.atoms[i].resid;
      else if (S.atoms[i].resid != m.resid0) m.many_units = true;            // several repeat units: a chain of them
      const int e = S.atoms[i].element;
      auto el = std::find_if(m.el.begin(), m.el.end(), [e](const auto& p) { return p.first == e; });
      if (el == m.el.end()) m.el.push_back({e, 1}); else ++el->second;
      m.z[size_t(std::clamp(int(zfrac(S.atoms[i]) * 8), 0, 7))] += 1;
      const uint8_t v = st ? d->atom_state[i] : 0;
      ((v & 3) == 1 ? m.ghost : (v & 3) == 2 ? m.hidden : m.shown)++;
      if (v & 4) ++m.locked;
      if (sel && d->selection[i]) ++m.sel;
    }
    std::sort(list.begin(), list.end(), [](const Mol& a, const Mol& b) { return a.id < b.id; });
    auto formula = [](const std::vector<std::pair<int, int>>& v) {
      const std::map<int, int> el(v.begin(), v.end());
      std::string f;
      auto add = [&](int z, int c) { f += caps::element(z).symbol; if (c > 1) f += std::to_string(c); };
      if (el.count(6)) { add(6, el.at(6)); if (el.count(1)) add(1, el.at(1)); }
      std::vector<std::pair<std::string, int>> rest;
      for (const auto& [z, c] : el) if (!(el.count(6) && (z == 6 || z == 1))) rest.push_back({caps::element(z).symbol, c});
      std::sort(rest.begin(), rest.end());
      for (const auto& [sym, c] : rest) { f += sym; if (c > 1) f += std::to_string(c); }
      return f;
    };
    // kinds in order of first appearance
    std::vector<std::string> order;
    std::map<std::string, caps::Json> kinds;
    std::map<std::string, size_t> kind_atoms;
    for (const auto& m : list) {
      const int64_t id = m.id;
      const std::string f = formula(m.el);
      // a residue name names a molecule only when it is one residue: a chain of units is its formula
      // one unit name over several units: a homopolymer chain, poly(NAME) — PE's ETH units make poly(ETH)
      const std::string res = m.many_res || m.res.empty() ? std::string() : m.many_units ? "poly(" + m.res + ")" : m.res;
      const std::string key = res + "|" + f;
      if (!kinds.count(key)) {
        order.push_back(key);
        caps::Json k = caps::Json::object();
        k["name"] = res.empty() ? f : res;
        k["formula"] = f;
        k["elements"] = [&] {
          std::vector<int> zs;
          for (const auto& p : m.el) zs.push_back(p.first);
          std::sort(zs.begin(), zs.end());
          std::string e;
          for (int z : zs) { if (!e.empty()) e += " "; e += caps::element(z).symbol; }
          return e;
        }();
        k["molecules"] = caps::Json::array();
        kinds[key] = k;
      }
      caps::Json x = caps::Json::object();
      x["id"] = double(id);
      x["atoms"] = double(m.atoms);
      const double peak = *std::max_element(m.z.begin(), m.z.end());
      caps::Json z = caps::Json::array();
      for (double b : m.z) z.push_back(peak > 0 ? std::round(b / peak * 100) / 100 : 0.0);
      x["z"] = z;
      x["state"] = m.hidden == m.atoms ? "hidden" : m.ghost == m.atoms ? "ghost" : m.shown == m.atoms ? "shown" : "mixed";
      x["locked"] = m.locked == m.atoms && m.atoms > 0;
      x["selected"] = double(m.sel);
      x["held"] = d->held_mol > 0 && id == d->held_mol;
      kinds[key]["molecules"].push_back(x);
      kind_atoms[key] += m.atoms;
    }
    // a residue name names its kind when it is specific (not a builder's placeholder) and names no other kind
    std::map<std::string, int> res_uses;
    for (const auto& key : order) ++res_uses[key.substr(0, key.find('|'))];
    static const std::set<std::string> generic{"", "MOL", "UNK", "UNL", "LIG", "RES", "X", "SYS", "DUM"};
    caps::Json arr = caps::Json::array();
    for (const auto& key : order) {
      const std::string res = key.substr(0, key.find('|')), f = key.substr(key.find('|') + 1);
      // the unit's own name (poly(NAME) → NAME) for the placeholder and water checks
      std::string up = res.rfind("poly(", 0) == 0 && res.size() > 6 ? res.substr(5, res.size() - 6) : res;
      for (auto& c : up) c = char(std::toupper(static_cast<unsigned char>(c)));
      const bool water = f == "H2O" || up == "SOL" || up == "HOH" || up == "WAT" || up == "TIP3" || up == "SPC";
      const bool placeholder = generic.count(up) || up.size() <= 1;   // a builder's unit letter (A, B …) names nothing
      kinds[key]["name"] = water ? std::string("water") : !placeholder && res_uses[res] == 1 ? res : f;
      kinds[key]["atoms"] = double(kind_atoms[key]);
      arr.push_back(kinds[key]);
    }
    r["kinds"] = arr;
    r["atoms"] = double(n);
    r["z_axis"] = S.cell.valid() ? "cell c" : "extent";
    return 0;
  });
  return report_out(r.dump(0), json, cap);
}

// Tags (design/boards/Tags): named, coloured atom sets kept in the structure (and its topology, so every frame has them).
int32_t caps_tags(caps_doc* d, char* json, int32_t cap) {
  caps::Json r = caps::Json::object();
  guard([&] {
    caps::Json a = caps::Json::array();
    const auto& S = d->frame;
    for (size_t k = 0; k < S.tags.size(); ++k) {
      caps::Json t = caps::Json::object();
      const auto atoms = caps::tag_atoms(S, int(k));
      size_t sel = 0;
      for (size_t i : atoms) sel += i < d->selection.size() && d->selection[i];
      t["name"] = S.tags[k].name;
      t["colour"] = S.tags[k].colour;
      t["count"] = double(atoms.size());
      t["selected"] = double(sel);
      t["group"] = caps::tag_group_name(S.tags[k].name);
      a.push_back(t);
    }
    r["tags"] = a;
    return 0;
  });
  return report_out(r.dump(0), json, cap);
}

int32_t caps_tag_edit(caps_doc* d, const char* json) {
  return guard([&] {
    const auto j = caps::Json::parse(json ? json : "{}");
    const std::string op = j.text("op", "set"), name = j.text("name"), colour = j.text("colour");
    std::vector<size_t> atoms;
    if (j.has("atoms") && j["atoms"].is_array())
      for (const auto& x : j["atoms"].items()) atoms.push_back(size_t(x.number()));
    else if (j.text("atoms") == "selection")
      for (size_t i = 0; i < d->selection.size(); ++i) if (d->selection[i]) atoms.push_back(i);
    int result = -1;
    auto apply = [&](caps::System& s) {
      const int k = caps::tag_index(s, name);
      if (op == "delete") { caps::tag_delete(s, k); result = k; }
      else if (op == "rename") {
        const std::string to = j.text("new_name");
        if (k < 0 || to.empty()) throw std::runtime_error("no tag " + name + " to rename");
        if (caps::tag_index(s, to) >= 0 && to != name) throw std::runtime_error("a tag " + to + " is there already");
        s.tags[size_t(k)].name = to;
        result = k;
      } else if (op == "colour") {
        if (k < 0) throw std::runtime_error("no tag " + name);
        s.tags[size_t(k)].colour = colour;
        result = k;
      } else result = caps::tag_edit(s, name, colour, atoms, op);
    };
    apply(d->frame);
    if (d->traj.topology.atoms.size() == d->frame.atoms.size()) apply(d->traj.topology);
    return result;
  });
}

// Clipboard pieces (design/boards/Stamp): {atoms: [0-based …] | "selection" | "all", name} → caps-piece JSON.
int32_t caps_piece(caps_doc* d, const char* json, char* out, int32_t cap) {
  std::string text;
  const int32_t rc = guard([&] {
    const auto j = caps::Json::parse(json ? json : "{}");
    std::vector<size_t> atoms;
    const auto& S = d->frame;
    if (j.has("atoms") && j["atoms"].is_array()) for (const auto& x : j["atoms"].items()) atoms.push_back(size_t(x.number()));
    else if (j.text("atoms") == "selection") { for (size_t i = 0; i < d->selection.size(); ++i) if (d->selection[i]) atoms.push_back(i); }
    else for (size_t i = 0; i < S.atoms.size(); ++i) atoms.push_back(i);
    text = caps::piece_json(S, atoms, j.text("name"));
    return 0;
  });
  if (rc < 0) return -1;
  return report_out(text, out, cap);
}

// A structure file as a piece (every atom): a drop into the open structure.
int32_t caps_piece_file(const char* path, const char* name, char* out, int32_t cap) {
  std::string text;
  const int32_t rc = guard([&] {
    const auto t = caps::open_file(path ? path : "", "");
    std::vector<size_t> all(t.topology.atoms.size());
    for (size_t i = 0; i < all.size(); ++i) all[i] = i;
    text = caps::piece_json(t.frames() > 0 ? t.frame(0) : t.topology, all, name && *name ? name : std::filesystem::path(path ? path : "").stem().string());
    return 0;
  });
  if (rc < 0) return -1;
  return report_out(text, out, cap);
}

// Look (design/boards/Look): {atom_scale (× vdW, ball and stick), bond_radius (Å), space_scale (× vdW), line_px, bond_orders,
// factor: {atoms: [0-based …] | "selection", value} (those atoms' size, 1 back to the rest), clear_factors}.
int32_t caps_set_look(caps_doc* d, const char* json) {
  return guard([&] {
    const auto j = caps::Json::parse(json ? json : "{}");
    auto& z = d->sizes;
    z.atom_scale = std::clamp(j.num("atom_scale", z.atom_scale), 0.05, 1.5);
    z.bond_radius = std::clamp(j.num("bond_radius", z.bond_radius), 0.02, 0.6);
    z.space_scale = std::clamp(j.num("space_scale", z.space_scale), 0.3, 1.5);
    z.line_px = std::clamp(j.num("line_px", z.line_px), 0.5, 6.0);
    if (j.has("bond_orders")) z.bond_orders = j["bond_orders"].boolean();
    if (j.has("clear_factors") && j["clear_factors"].boolean()) z.factor.clear();
    if (j.has("factor") && j["factor"].is_object()) {
      const auto& f = j["factor"];
      const size_t n = d->frame.atoms.size();
      if (z.factor.size() != n) z.factor.assign(n, 1.0f);
      const float v = float(std::clamp(f.num("value", 1.0), 0.1, 4.0));
      if (f.has("atoms") && f["atoms"].is_array()) { for (const auto& x : f["atoms"].items()) if (size_t(x.number()) < n) z.factor[size_t(x.number())] = v; }
      else for (size_t i = 0; i < n && i < d->selection.size(); ++i) if (d->selection[i]) z.factor[i] = v;
    }
    return 0;
  });
}

// Bond rules (design/boards/BondRules): rules from JSON [{a, b, max, never}] (element symbols)
std::vector<caps::BondRule> rules_of(const caps::Json& j) {
  std::vector<caps::BondRule> r;
  if (!j.has("rules") || !j["rules"].is_array()) return r;
  for (const auto& x : j["rules"].items()) {
    caps::BondRule b;
    b.za = caps::element_from_symbol(x.text("a")), b.zb = caps::element_from_symbol(x.text("b"));
    b.max = x.num("max", 0);
    b.never = x.has("never") && x["never"].boolean();
    if (b.za > 0 && b.zb > 0) r.push_back(b);
  }
  return r;
}

// {pairs: [{a, b, counts, lo, bin, bonded, covalent, suggested, ionic}], bonds}
int32_t caps_pair_histograms(caps_doc* d, const char* json, char* out, int32_t cap) {
  caps::Json r = caps::Json::object();
  const int32_t rc = guard([&] {
    const auto j = caps::Json::parse(json && *json ? json : "{}");
    const auto hs = caps::pair_histograms(d->frame, j.num("lo", 0.8), j.num("hi", 3.2), j.num("bin", 0.04));
    caps::Json a = caps::Json::array();
    for (const auto& h : hs) {
      caps::Json o = caps::Json::object();
      o["a"] = caps::element(h.za).symbol, o["b"] = caps::element(h.zb).symbol;
      caps::Json c = caps::Json::array();
      for (int v : h.counts) c.push_back(double(v));
      o["counts"] = c;
      o["lo"] = h.lo, o["bin"] = h.bin, o["bonded"] = double(h.bonded), o["covalent"] = h.covalent, o["suggested"] = h.suggested, o["ionic"] = h.ionic;
      a.push_back(o);
    }
    r["pairs"] = a;
    r["bonds"] = double(d->frame.bonds.size());
    return 0;
  });
  if (rc < 0) return -1;
  return report_out(r.dump(0), out, cap);
}

// The rules' bonds compared with the structure's: {bonds, after, added, removed, pairs: {"B–N": +33}}; preview: the new
// bonds drawn in the view (thin accent sticks), else that preview cleared.
int32_t caps_bond_rules_preview(caps_doc* d, const char* json, char* out, int32_t cap) {
  caps::Json r = caps::Json::object();
  const int32_t rc = guard([&] {
    const auto j = caps::Json::parse(json && *json ? json : "{}");
    d->overlay.clear();
    if (!j.has("rules")) return 0;
    const auto& S = d->frame;
    const auto nb = caps::bonds_by_rules(S, rules_of(j));
    auto pk = [](uint32_t a, uint32_t b) { return a < b ? (uint64_t(a) << 32) | b : (uint64_t(b) << 32) | a; };
    std::set<uint64_t> now, after;
    for (const auto& b : S.bonds) now.insert(pk(b.i, b.j));
    for (const auto& b : nb) after.insert(pk(b.i, b.j));
    std::map<std::string, long> per;
    size_t added = 0, removed = 0;
    auto name = [&](uint32_t i, uint32_t k) {
      std::string a = caps::element(S.atoms[i].element).symbol, b = caps::element(S.atoms[k].element).symbol;
      if (S.atoms[i].element > S.atoms[k].element) std::swap(a, b);
      return a + "–" + b;
    };
    for (const auto& b : nb)
      if (!now.count(pk(b.i, b.j))) {
        ++added, ++per[name(b.i, b.j)];
        if (j.has("preview") && j["preview"].boolean() && d->overlay.size() < 20000) {
          const caps::Vec3 a = S.atoms[b.i].pos, dv = S.cell.valid() ? S.cell.minimum_image(S.atoms[b.j].pos - a) : S.atoms[b.j].pos - a;
          d->overlay.push_back({a, a + dv, 0xF0A83C, 0.07, false});
        }
      }
    for (const auto& b : S.bonds)
      if (!after.count(pk(b.i, b.j))) ++removed, --per[name(b.i, b.j)];
    r["bonds"] = double(S.bonds.size()), r["after"] = double(nb.size()), r["added"] = double(added), r["removed"] = double(removed);
    caps::Json p = caps::Json::object();
    for (const auto& [k, v] : per) if (v != 0) p[k] = double(v);
    r["pairs"] = p;
    return 0;
  });
  if (rc < 0) return -1;
  return report_out(r.dump(0), out, cap);
}

// Probes (design/boards/Probes): [{kind, atoms: [0-based …], rgb}] drawn in the view (replacing those set before).
int32_t caps_set_probes(caps_doc* d, const char* json) {
  return guard([&] {
    const auto j = caps::Json::parse(json && *json ? json : "[]");
    std::vector<caps_doc::ProbeDef> list;
    for (const auto& x : j.is_array() ? j.items() : std::vector<caps::Json>{}) {
      caps_doc::ProbeDef p;
      p.kind = caps::probe_kind(x.text("kind", "point"));
      if (x.has("atoms")) for (const auto& a : x["atoms"].items()) p.atoms.push_back(size_t(a.number()));
      p.rgb = unsigned(x.num("rgb", 0x6CC4D8));
      list.push_back(std::move(p));
    }
    d->probes = std::move(list);
    return 0;
  });
}

// A probe's geometry in the frame shown: {kind, atoms} → {centre, axes: [[…] largest first], semi, rms}.
int32_t caps_probe_geometry(caps_doc* d, const char* json, char* out, int32_t cap) {
  caps::Json r = caps::Json::object();
  const int32_t rc = guard([&] {
    const auto j = caps::Json::parse(json ? json : "{}");
    std::vector<size_t> at;
    if (j.has("atoms") && j["atoms"].is_array()) for (const auto& v : j["atoms"].items()) at.push_back(size_t(v.number()));
    else for (size_t i = 0; i < d->selection.size(); ++i) if (d->selection[i]) at.push_back(i);
    const auto p = caps::make_probe(d->frame, at, caps::probe_kind(j.text("kind", "ellipsoid")));
    auto vec = [](const caps::Vec3& v) { caps::Json a = caps::Json::array(); a.push_back(v[0]), a.push_back(v[1]), a.push_back(v[2]); return a; };
    r["centre"] = vec(p.centre);
    caps::Json ax = caps::Json::array(), se = caps::Json::array();
    for (int k = 0; k < 3; ++k) ax.push_back(vec(p.axes[size_t(k)])), se.push_back(p.semi[size_t(k)]);
    r["axes"] = ax, r["semi"] = se, r["rms"] = p.rms;
    return 0;
  });
  if (rc < 0) return -1;
  return report_out(r.dump(0), out, cap);
}

// A measurement between probes over every frame: {a: {kind, atoms}, b: {kind, atoms} | absent, measure} → {values}.
int32_t caps_probe_series(caps_doc* d, const char* json, char* out, int32_t cap) {
  caps::Json r = caps::Json::object();
  const int32_t rc = guard([&] {
    const auto j = caps::Json::parse(json ? json : "{}");
    auto atoms_of = [](const caps::Json& x) { std::vector<size_t> a; if (x.has("atoms")) for (const auto& v : x["atoms"].items()) a.push_back(size_t(v.number())); return a; };
    const auto ka = caps::probe_kind(j["a"].text("kind", "point"));
    const auto aa = atoms_of(j["a"]);
    const bool hasb = j.has("b") && j["b"].is_object();
    const auto kb = hasb ? caps::probe_kind(j["b"].text("kind", "point")) : caps::ProbeKind::Point;
    const auto ab = hasb ? atoms_of(j["b"]) : std::vector<size_t>{};
    const std::string what = j.text("measure", "distance");
    caps::Json v = caps::Json::array();
    const size_t nf = std::max<size_t>(1, d->traj.frames());
    for (size_t f = 0; f < nf; ++f) {
      const caps::System s = d->traj.frames() > 0 ? d->traj.frame(f) : d->frame;
      const auto pa = caps::make_probe(s, aa, ka);
      const caps::Probe pb = hasb ? caps::make_probe(s, ab, kb) : caps::Probe{};
      v.push_back(caps::probe_measure(s, pa, hasb ? &pb : nullptr, what));
    }
    r["values"] = v;
    r["unit"] = what == "angle" ? "°" : "Å";
    return 0;
  });
  if (rc < 0) return -1;
  return report_out(r.dump(0), out, cap);
}

// Conformers (C12, C13): the frame shown, or only its selected atoms (with their own share of the force field), in vacuum
std::pair<caps::System, std::shared_ptr<const caps::ForceField>> conformer_input(caps_doc* d, bool selection) {
  caps::System s = d->traj.frames() > 0 ? d->traj.frame(d->current) : d->frame;
  if (!s.unwrapped) caps::make_molecules_whole(s);
  auto ff = field_for_run(d);
  if (!ff) ff = std::make_shared<caps::ForceField>(default_ff(s));
  std::vector<uint32_t> keep;
  if (selection)
    for (size_t i = 0; i < s.atoms.size() && i < d->selection.size(); ++i)
      if (d->selection[i]) keep.push_back(uint32_t(i));
  if (!keep.empty()) {   // nothing selected: the whole frame
    std::vector<int64_t> newi(s.atoms.size(), -1);
    for (size_t k = 0; k < keep.size(); ++k) newi[keep[k]] = int64_t(k);
    caps::System t = s;
    t.atoms.clear(), t.bonds.clear(), t.velocities.clear();
    for (uint32_t i : keep) t.atoms.push_back(s.atoms[i]);
    for (const auto& b : s.bonds)
      if (newi[b.i] >= 0 && newi[b.j] >= 0) t.bonds.push_back({uint32_t(newi[b.i]), uint32_t(newi[b.j]), b.order});
    ff = std::make_shared<caps::ForceField>(caps::subset_forcefield(*ff, keep));
    s = std::move(t);
  }
  if (s.atoms.size() > 400) throw std::invalid_argument(std::to_string(s.atoms.size()) + " atoms: a conformer search is for one molecule — select it");
  s.cell = caps::Cell{};
  return {s, ff};
}

caps::Json conformers_json(const caps::ConformerSearchResult& r) {
  caps::Json j = caps::Json::object(), cs = caps::Json::array(), nt = caps::Json::array();
  for (const auto& c : r.conformers) {
    caps::Json x = caps::Json::object();
    x["energy"] = c.energy, x["relative"] = c.relative, x["population"] = c.population, x["found"] = double(c.found);
    cs.push_back(std::move(x));
  }
  for (const auto& n : r.notes) nt.push_back(n);
  j["conformers"] = std::move(cs);
  j["rotors"] = double(r.rotors);
  j["trials"] = double(r.minima.size());
  j["notes"] = std::move(nt);
  return j;
}

int32_t caps_conformer_frames(caps_doc* d, const char* json, char* out, int32_t cap) {
  caps::Json r = caps::Json::object();
  const int32_t rc = guard([&] {
    const auto j = caps::Json::parse(json ? json : "{}");
    auto [s, ff] = conformer_input(d, j.has("selection") && j["selection"].boolean());
    caps::ConformerSearchOptions co;
    co.field = ff;
    co.energy = elec(caps::EnergyOptions{});
    co.method = j.text("method", "torsions");
    co.trials = std::clamp(int(j.num("trials", 50)), 1, 5000);
    co.window = j.num("window", 10);
    co.rmsd = j.num("rmsd", 0.5);
    co.temperature = j.num("temperature", 298.15);
    co.seed = uint64_t(j.num("seed", 1));
    const auto res = caps::conformer_search(s, co);
    caps::Trajectory t;
    t.topology = s;
    t.topology.unwrapped = true;
    for (const auto& c : res.conformers) t.positions.push_back(c.pos), t.cells.push_back(caps::Cell{}), t.timesteps.push_back(int64_t(t.timesteps.size()));
    prov_step(d, "conformers", "Conformer search", {{"method", co.method + ", " + std::to_string(res.minima.size()) + " minimised starts"},
              {"rotatable bonds", std::to_string(res.rotors)}, {"kept", std::to_string(res.conformers.size()) + " within " + g6(co.window) + " kcal/mol, " + g6(co.rmsd) + " Å heavy-atom RMSD"},
              {"force field", ff->name}}, seeded(co.seed), {"horn1987"});
    d->traj = std::move(t);
    d->current = 0;
    refresh(d);
    r = conformers_json(res);
    return 0;
  });
  if (rc < 0) return -1;
  return report_out(r.dump(0), out, cap);
}

// Normal modes (C20): mode k of the frame shown as frames for the player — one period, the largest atom displacement
// `amplitude` Å. The document's frames become the animation (the structure itself is its first frame).
int32_t caps_mode_animate(caps_doc* d, const char* json, char* out, int32_t cap) {
  caps::Json r = caps::Json::object();
  const int32_t rc = guard([&] {
    const auto j = caps::Json::parse(json ? json : "{}");
    caps::System s = d->traj.frames() > 0 ? d->traj.frame(d->current) : d->frame;
    if (!s.unwrapped) caps::make_molecules_whole(s);
    caps::NormalModesOptions no;
    no.field = field_for_run(d);
    no.energy = elec(caps::EnergyOptions{});
    if (j.has("cutoff")) no.energy.cutoff = j.num("cutoff", 10);
    if (j.text("atoms", "") == "selection" && d->selection.size() == s.atoms.size()) no.moving.assign(d->selection.begin(), d->selection.end());
    const auto res = caps::normal_modes(s, no);
    const int k = int(j.num("mode", 1)) - 1;
    if (k < 0 || size_t(k) >= res.mode.size()) throw std::invalid_argument("mode " + std::to_string(k + 1) + " of " + std::to_string(res.mode.size()));
    const double amp = std::clamp(j.num("amplitude", 0.3), 0.01, 2.0);
    const int nfr = std::clamp(int(j.num("frames", 20)), 4, 200);
    caps::Trajectory t;
    t.topology = s;
    t.topology.velocities.clear();
    t.topology.unwrapped = true;
    for (const auto& p : caps::mode_frames(s, res, size_t(k), amp, nfr))
      t.positions.push_back(p), t.cells.push_back(s.cell), t.timesteps.push_back(int64_t(t.timesteps.size()));
    prov_step(d, "modes.animate", "Normal mode animation", {{"mode", std::to_string(k + 1) + " of " + std::to_string(res.mode.size())},
              {"wavenumber", g6(res.wavenumber[k]) + " cm⁻¹"}, {"amplitude", g6(amp) + " Å, " + std::to_string(nfr) + " frames per period"},
              {"force field", ff_label(d)}}, "", {"miller1980"});
    d->traj = std::move(t);
    d->current = 0;
    refresh(d);
    r["wavenumber"] = res.wavenumber[k];
    r["modes"] = double(res.mode.size());
    r["reduced_mass"] = res.reduced_mass[k];
    r["ir"] = res.ir[k];
    r["frames"] = double(nfr);
    caps::Json nt = caps::Json::array();
    for (const auto& x : res.notes) nt.push_back(x);
    r["notes"] = std::move(nt);
    return 0;
  });
  if (rc < 0) return -1;
  return report_out(r.dump(0), out, cap);
}

// Brush to select (design/boards/BrushSelect): a per-atom column for a histogram — atom_column's names, and
// "distance:N" (Å, minimum image, from atom N, 0-based) or "distance:tag:NAME" (from the nearest atom of a tag).
int32_t caps_atom_column(caps_doc* d, const char* name, double* out, int32_t cap) {
  return guard([&] {
    const auto& S = d->frame;
    const std::string w = name ? name : "";
    std::vector<double> v;
    if (w.rfind("distance:", 0) == 0) {
      std::vector<size_t> from;
      const std::string what = w.substr(9);
      if (what.rfind("tag:", 0) == 0) {
        from = caps::tag_atoms(S, caps::tag_index(S, what.substr(4)));
        if (from.empty()) throw std::invalid_argument("no atoms in tag " + what.substr(4));
      } else {
        const long k = std::stol(what);
        if (k < 0 || size_t(k) >= S.atoms.size()) throw std::invalid_argument("no atom " + what);
        from.push_back(size_t(k));
      }
      v.assign(S.atoms.size(), 1e300);
      for (size_t i = 0; i < S.atoms.size(); ++i)
        for (size_t f : from) {
          const caps::Vec3 dv = S.cell.valid() ? S.cell.minimum_image(S.atoms[i].pos - S.atoms[f].pos) : S.atoms[i].pos - S.atoms[f].pos;
          v[i] = std::min(v[i], std::sqrt(dv[0] * dv[0] + dv[1] * dv[1] + dv[2] * dv[2]));
        }
    } else {
      v = caps::atom_column(S, w);
    }
    if (out) for (size_t i = 0; i < v.size() && int32_t(i) < cap; ++i) out[i] = v[i];
    return int32_t(v.size());
  });
}

int32_t caps_tag_atoms(caps_doc* d, const char* name, int32_t* out, int32_t cap) {
  const auto atoms = caps::tag_atoms(d->frame, caps::tag_index(d->frame, name ? name : ""));
  if (out) for (size_t k = 0; k < atoms.size() && int32_t(k) < cap; ++k) out[k] = int32_t(atoms[k]);
  return int32_t(atoms.size());
}

int32_t caps_atom_states(caps_doc* d, uint8_t* out, int32_t n) {
  const size_t na = d->frame.atoms.size();
  if (out && n >= int32_t(na)) {
    if (d->atom_state.size() == na) std::copy(d->atom_state.begin(), d->atom_state.end(), out);
    else std::fill(out, out + na, uint8_t(0));
  }
  return int32_t(na);
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

// Molecule ids "1-4,7" (commas, semicolons or spaces between; a-b ranges inclusive).
static std::vector<int64_t> parse_mol_ranges(std::string t) {
  for (auto& c : t) if (c == ',' || c == ';') c = ' ';
  std::istringstream is(t);
  std::set<int64_t> want;
  for (std::string w; is >> w;) {
    const auto dash = w.find('-', 1);
    int64_t a = 0, b = 0;
    try {
      a = std::stoll(w.substr(0, dash));
      b = dash == std::string::npos ? a : std::stoll(w.substr(dash + 1));
    } catch (const std::exception&) { throw std::invalid_argument("molecule ids: cannot read " + w + " (write 1-4,7)"); }
    if (b < a || b - a > 10000000) throw std::invalid_argument("molecule ids: bad range " + w);
    for (int64_t m = a; m <= b; ++m) want.insert(m);
  }
  return {want.begin(), want.end()};
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
          if (!cur.empty()) (cur == "cij_strain" || cur == "cij_run" || cur == "viscosity" || cur == "nemd" || cur == "conformers" || cur == "creep" || cur == "friction" || cur == "solvation" || cur == "tensile" || cur == "tg" || cur == "pull_shear" || cur == "pull_normal" ? protocols : ids).push_back(cur);
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
      if (p->radii && *p->radii) o.radii = p->radii;
      if (p->zbin > 0) o.zbin = p->zbin;
      if (p->axis >= 1 && p->axis <= 3) o.axis = p->axis - 1;
      if (p->surface && *p->surface) o.surface_mols = parse_mol_ranges(p->surface);
    }
    caps_mech_opts mo{};
    if (m) mo = *m;
    if (mo.temperature > 0) o.temperature = mo.temperature;
    if (progress) o.progress = [&](const std::string& what, double f) { return progress(what.c_str(), f, user) != 0; };
    auto has = [&](const char* k) { return std::find(ids.begin(), ids.end(), k) != ids.end(); };
    // force field: the Field assignment (must be complete), else GAFF of C and H
    std::shared_ptr<const caps::ForceField> ff = field_for_run(d);
    std::vector<std::string> extra_notes;
    if (!ff && (has("ced") || has("delta") || has("cij_fluct") || has("fluct") || has("modes") || has("adhesion") || !protocols.empty())) {
      caps::System s0 = d->traj.frame(0);
      if (!s0.unwrapped) caps::make_molecules_whole(s0);
      ff = std::make_shared<caps::ForceField>(default_ff(s0));
      extra_notes.push_back("force field: " + ff->name + " (none assigned in Field)");
    }
    if (ff) o.ff = ff.get();
    // a group: the properties see a trajectory of those atoms alone (bonds among them, the force field's terms among them)
    const std::string group = p && p->group ? p->group : "";
    std::vector<caps::Property> res;
    if (!group.empty() && group != "all" && !ids.empty()) {
      const caps::System& top = d->traj.topology;
      const size_t n = top.atoms.size();
      std::vector<char> in(n, 0);
      if (group == "selection") {
        if (d->selection.size() == n) in.assign(d->selection.begin(), d->selection.end());
      } else if (group == "exclude-held") {
        for (size_t i = 0; i < n; ++i) in[i] = d->held_mol <= 0 || top.atoms[i].mol != d->held_mol;
      } else if (group.rfind("molecules:", 0) == 0) {
        const auto ids = parse_mol_ranges(group.substr(10));
        const std::set<int64_t> want(ids.begin(), ids.end());
        const auto mol = top.molecules();
        for (size_t i = 0; i < n; ++i) in[i] = want.count(top.has_mol ? top.atoms[i].mol : int64_t(mol[i]) + 1) > 0;
      } else if (group.rfind("tag:", 0) == 0) {   // a tag of the structure
        const int k = caps::tag_index(top, group.substr(4));
        if (k < 0) throw std::invalid_argument("no tag " + group.substr(4) + " in this structure");
        for (size_t i = 0; i < n; ++i) in[i] = (top.atoms[i].tags >> k) & 1u;
      } else {
        throw std::invalid_argument("group: \"selection\", \"molecules:1-4,7\", \"tag:NAME\" or \"exclude-held\"");
      }
      std::vector<uint32_t> keep;
      for (size_t i = 0; i < n; ++i) if (in[i]) keep.push_back(uint32_t(i));
      if (keep.empty()) throw std::invalid_argument("the group has no atoms");
      std::vector<int64_t> newi(n, -1);
      for (size_t k = 0; k < keep.size(); ++k) newi[keep[k]] = int64_t(k);
      caps::Trajectory sub;
      sub.topology = top;
      sub.topology.atoms.clear();
      sub.topology.bonds.clear();
      sub.topology.velocities.clear();
      for (uint32_t i : keep) sub.topology.atoms.push_back(top.atoms[i]);
      for (const auto& b : top.bonds)
        if (newi[b.i] >= 0 && newi[b.j] >= 0) sub.topology.bonds.push_back({uint32_t(newi[b.i]), uint32_t(newi[b.j]), b.order});
      sub.cells = d->traj.cells;
      sub.timesteps = d->traj.timesteps;
      for (const auto& fr : d->traj.positions) {
        std::vector<caps::Vec3> q;
        q.reserve(keep.size());
        for (uint32_t i : keep) q.push_back(fr[i]);
        sub.positions.push_back(std::move(q));
      }
      caps::ForceField subff;
      caps::AnalyzeOptions og = o;
      if (ff) { subff = caps::subset_forcefield(*ff, keep); og.ff = &subff; }
      og.exclude_mol = group == "exclude-held" ? 0 : o.exclude_mol;
      res = caps::analyze(sub, ids, og);
      for (auto& q : res) q.notes.insert(q.notes.begin(), "group " + group + ": " + std::to_string(keep.size()) + " of " + std::to_string(n) + " atoms");
    } else if (!ids.empty()) {
      res = caps::analyze(d->traj, ids, o);
    }
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
      } else if (id == "solvation") {
        caps::System s = frame_copy();
        caps::SolvationOptions so;
        so.field = ff;
        so.energy = o.energy;
        so.solute.assign(s.atoms.size(), 0);
        if (mo.solv_mol > 0) {
          const auto ids = molecule_ids(s);
          for (size_t i = 0; i < s.atoms.size(); ++i) so.solute[i] = ids[i] == mo.solv_mol;
        } else {
          for (size_t i = 0; i < s.atoms.size() && i < d->selection.size(); ++i) so.solute[i] = d->selection[i] ? 1 : 0;
        }
        auto spaced = [](int k) { std::vector<double> v; for (int i = 0; i < k; ++i) v.push_back(1.0 - double(i) / (k - 1)); return v; };
        if (mo.solv_coul_windows >= 2) so.coul_windows = spaced(mo.solv_coul_windows);
        if (mo.solv_lj_windows >= 2) {   // denser towards λ = 0, where the soft-core ⟨∂U/∂λ⟩ turns fastest
          so.lj_windows = spaced(mo.solv_lj_windows);
          for (auto& l : so.lj_windows) l = std::pow(l, 1.5);
        }
        if (mo.solv_ps > 0) so.ps = mo.solv_ps;
        if (mo.solv_eq_ps != 0) so.equilibrate_ps = std::max(0.0, mo.solv_eq_ps);
        so.temperature = mo.solv_t > 0 ? mo.solv_t : mo.temperature > 0 ? mo.temperature : 300;
        if (mo.dt > 0) so.dt = mo.dt;
        if (mo.seed) so.seed = mo.seed;
        so.progress = [&](const std::string& w, double f) { return !cancelled(w, f); };
        const auto sr = caps::solvation_free_energy(s, so);
        caps::Property q;
        q.id = "solvation";
        q.name = "Solvation free energy (TI)";
        q.unit = "kcal/mol";
        q.value = sr.dg, q.error = sr.dg_err;
        q.method = "thermodynamic integration: electrostatics then soft-core Lennard-Jones (Beutler et al. 1994) of the solute decoupled from the rest, trapezoid rule";
        q.extra["electrostatic part (kcal/mol)"] = sr.dg_coul;
        q.extra["Lennard-Jones part (kcal/mol)"] = sr.dg_lj;
        q.extra["ΔG (kJ/mol)"] = sr.dg * 4.184;
        q.extra["solute atoms"] = double(std::count(so.solute.begin(), so.solute.end(), 1));
        caps::Series sc{"⟨∂U/∂λ⟩ electrostatics", "λ_coul", "⟨∂U/∂λ⟩ (kcal/mol)", {}, {}}, sl{"⟨∂U/∂λ⟩ Lennard-Jones", "λ_lj", "⟨∂U/∂λ⟩ (kcal/mol)", {}, {}};
        for (const auto& w : sr.windows) (w.leg == "coulomb" ? sc : sl).x.push_back(w.lambda), (w.leg == "coulomb" ? sc : sl).y.push_back(w.dudl);
        if (!sc.x.empty()) q.series.push_back(sc);
        q.series.push_back(sl);
        q.notes = sr.notes;
        res.push_back(std::move(q));
      } else if (id == "friction") {
        caps::System s = frame_copy();
        caps::FrictionOptions fo;
        fo.field = ff;
        fo.energy = o.energy;
        if (mo.fr_moving > 0) fo.moving_mol = mo.fr_moving;
        fo.fixed_mol = std::max(0, mo.fr_fixed);
        if (mo.fr_velocity != 0) fo.velocity = mo.fr_velocity;
        if (mo.fr_ps > 0) fo.ps = mo.fr_ps;
        if (mo.fr_eq_ps != 0) fo.equilibrate_ps = std::max(0.0, mo.fr_eq_ps);
        fo.temperature = mo.fr_t > 0 ? mo.fr_t : mo.temperature > 0 ? mo.temperature : 300;
        if (mo.dt > 0) fo.dt = mo.dt;
        if (mo.seed) fo.seed = mo.seed;
        fo.progress = [&](const std::string& w, double f) { return !cancelled(w, f); };
        for (auto& q : caps::friction_properties(caps::run_friction(s, fo))) res.push_back(std::move(q));
      } else if (id == "creep") {
        caps::System s = frame_copy();
        caps::CreepOptions co;
        co.field = ff;
        co.energy = o.energy;
        co.axis = std::clamp(mo.creep_axis, 0, 2);
        if (mo.creep_stress != 0) co.stress = mo.creep_stress;
        if (mo.creep_t > 0) co.temperature = mo.creep_t;
        if (mo.creep_ps > 0) co.ps = mo.creep_ps;
        if (mo.creep_eq_ps != 0) co.equilibrate_ps = std::max(0.0, mo.creep_eq_ps);
        if (mo.pressure > 0) co.pressure = mo.pressure;
        if (mo.dt > 0) co.dt = mo.dt;
        if (mo.seed) co.seed = mo.seed;
        co.progress = [&](const std::string& w, double f) { return !cancelled(w, f); };
        for (auto& q : caps::creep_properties(caps::run_creep(s, co), co.stress)) res.push_back(std::move(q));
      } else if (id == "conformers") {
        auto [cs, cff] = conformer_input(d, mo.conf_selection != 0);
        caps::ConformerSearchOptions co;
        co.field = cff;
        co.energy = o.energy;
        co.method = mo.conf_method == 1 ? "anneal" : "torsions";
        if (mo.conf_trials > 0) co.trials = mo.conf_trials;
        if (mo.conf_window > 0) co.window = mo.conf_window;
        if (mo.conf_rmsd > 0) co.rmsd = mo.conf_rmsd;
        if (mo.temperature > 0) co.temperature = mo.temperature;
        if (mo.seed) co.seed = mo.seed;
        co.progress = [&](double f) { return !cancelled("conformer search", f); };
        const auto cr = caps::conformer_search(cs, co);
        caps::Property q;
        q.id = "conformers";
        q.name = "Conformers";
        q.unit = "";
        q.value = double(cr.conformers.size());
        q.method = (co.method == "anneal" ? "NVT at 1000 K, snapshots quenched" : "random staggered torsions") + std::string(", ") + std::to_string(cr.minima.size()) +
                   " minimised (" + cr.field + ", vacuum); clustered by heavy-atom RMSD < " + g6(co.rmsd) + " Å after superposition";
        q.extra["rotatable bonds"] = cr.rotors;
        q.extra["minimised starts"] = double(cr.minima.size());
        for (size_t k = 0; k < cr.conformers.size() && k < 12; ++k) {
          char key[80];
          std::snprintf(key, sizeof key, "#%zu ΔE (kcal/mol)", k + 1);
          q.extra[key] = cr.conformers[k].relative;
          std::snprintf(key, sizeof key, "#%zu population", k + 1);
          q.extra[key] = cr.conformers[k].population;
        }
        caps::Series se{"conformer energies", "conformer", "ΔE (kcal/mol)", {}, {}}, sm{"minima found", "start", "E (kcal/mol)", {}, {}};
        for (size_t k = 0; k < cr.conformers.size(); ++k) se.x.push_back(double(k + 1)), se.y.push_back(cr.conformers[k].relative);
        for (size_t k = 0; k < cr.minima.size(); ++k) sm.x.push_back(double(k + 1)), sm.y.push_back(cr.minima[k]);
        q.series = {se, sm};
        q.notes = cr.notes;
        res.push_back(std::move(q));
      } else if (id == "nemd") {
        caps::NemdOptions no;
        no.field = ff;
        no.energy = o.energy;
        if (mo.temperature > 0) no.temperature = mo.temperature;
        if (mo.shear_ps > 0) no.ps = mo.shear_ps;
        if (mo.shear_eq_ps != 0) no.equilibrate_ps = std::max(0.0, mo.shear_eq_ps);
        if (mo.dt > 0) no.dt = mo.dt;
        if (mo.seed) no.seed = mo.seed;
        if (mo.shear_lo > 0) {
          const double lo = mo.shear_lo, hi = std::max(mo.shear_hi, lo);
          const int np = hi > lo ? std::clamp(mo.shear_points, 2, 12) : 1;
          no.rates.clear();
          for (int k = 0; k < np; ++k) no.rates.push_back(np == 1 ? lo : lo * std::pow(hi / lo, double(k) / (np - 1)));
        }
        no.progress = [&](const std::string& w, double f) { return !cancelled(w, f); };
        for (auto& q : caps::nemd_properties(caps::nemd_viscosity(frame_copy(), no))) res.push_back(std::move(q));
      } else if (id == "tensile") {
        caps::System s = frame_copy();
        caps::TensileOptions to;
        to.field = ff;
        to.energy = o.energy;
        to.axis = std::clamp(mo.axis, 0, 3);   // 3: x, y and z averaged
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
        if (mo.barostat < 0 || mo.barostat > 2) throw std::invalid_argument("barostat: 0 c-rescale, 1 Berendsen, 2 MTK");
        co.barostat = mo.barostat == 1 ? caps::Barostat::Berendsen : mo.barostat == 2 ? caps::Barostat::MTK : caps::Barostat::CRescale;
        if (mo.tau_t > 0) co.tau_t = mo.tau_t;
        if (mo.tau_p > 0) co.tau_p = mo.tau_p;
        if (mo.average_from > 0) {
          if (mo.average_from >= 0.95) throw std::invalid_argument("average_from: the fraction of each hold discarded, below 0.95");
          co.average_from = mo.average_from;
        }
        co.property = std::clamp(mo.tg_property, 0, 1);
        co.fit = std::clamp(mo.tg_fit, 0, 1);
        co.glassy_max = std::max(0.0, mo.glassy_max);
        co.rubbery_min = std::max(0.0, mo.rubbery_min);
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
    // the results go into the structure's provenance (the Project table and the methods section read them from there)
    caps::KeyValues pr;
    std::string names;
    for (const auto& r : res) {
      if (!std::isfinite(r.value)) continue;
      pr.push_back({r.name, g6(r.value) + (std::isfinite(r.error) ? " ± " + g6(r.error) : "") + (r.unit.empty() ? "" : " " + r.unit)});
      names += (names.empty() ? "" : ", ") + r.id;
    }
    if (!pr.empty()) {
      pr.push_back({"frames", std::to_string(fr.size()) + " of " + std::to_string(d->traj.frames())});
      if (!group.empty()) pr.push_back({"group", group});
      prov_step(d, "analyze.properties", "properties: " + names, std::move(pr));
    }
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

int32_t caps_molecule_ids(caps_doc* d, int64_t* out, int32_t cap) {
  return guard([&] {
    const auto ids = molecule_ids(d->frame);
    for (size_t j = 0; j < ids.size() && j < size_t(std::max(0, cap)); ++j) out[j] = ids[j];
    return int32_t(ids.size());
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
      b.rotor_search = o->rotor_search != 0;
      b.implicit_hydrogens = o->heavy_only == 0;
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
    // T6 / T7: the reference ranges beside the force-field library; the user's equilibrated cells in CAPS_BENCH_CELLS,
    // else ~/.caps/bench/cells
    if (!o.forcefields.empty()) o.reference = (std::filesystem::path(o.forcefields).parent_path() / "reference" / "polymers.json").string();
    if (const char* c = std::getenv("CAPS_BENCH_CELLS"); c && *c) o.cells = c;
    else if (const char* h = std::getenv("HOME"); h && *h) o.cells = std::string(h) + "/.caps/bench/cells";
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
  c.dyads = j.text("dyads", "");
  for (char ch : c.dyads)
    if (ch != 'm' && ch != 'r' && ch != 'M' && ch != 'R') throw std::invalid_argument("dyads: a pattern of m (meso) and r (racemo), e.g. \"mr\"");
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
  c.head_cap = j.text("head_cap", "");   // end groups: a preset or a SMILES with *
  c.tail_cap = j.text("tail_cap", "");
  c.linkage = caps::linkage_from_string(j.text("linkage", "head-to-tail"));
  c.inversion = j.num("inversion", c.inversion);
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

namespace {
// Hill formulas as element counts and back (C, H first; alphabetical without carbon)
std::map<std::string, int> formula_counts(const std::string& f) {
  std::map<std::string, int> n;
  for (size_t i = 0; i < f.size();) {
    if (!std::isupper(static_cast<unsigned char>(f[i]))) { ++i; continue; }
    std::string e(1, f[i++]);
    while (i < f.size() && std::islower(static_cast<unsigned char>(f[i]))) e += f[i++];
    int k = 0;
    while (i < f.size() && std::isdigit(static_cast<unsigned char>(f[i]))) k = 10 * k + (f[i++] - '0');
    n[e] += k ? k : 1;
  }
  return n;
}
std::string hill(const std::map<std::string, int>& n) {
  std::string out;
  auto put = [&](const std::string& e) { const auto it = n.find(e); if (it != n.end() && it->second > 0) out += e + (it->second > 1 ? std::to_string(it->second) : ""); };
  const bool carbon = n.count("C") && n.at("C") > 0;
  if (carbon) { put("C"); put("H"); }
  for (const auto& [e, k] : n) if (!(carbon && (e == "C" || e == "H"))) put(e);
  return out;
}
}  // namespace

extern "C" int32_t caps_chain_preview(const char* spec_json, uint64_t seed, char* json, int32_t cap) {
  caps::Json j = caps::Json::object();
  try {
    const caps::ChainSpec c = spec_from(spec_json ? spec_json : "{}");
    if (c.units.empty()) throw std::runtime_error("no repeat unit");
    const auto seq = caps::chain_sequence(c, seed);
    const auto inv = caps::chain_inversions(c, seq.size(), seed);
    const caps::MolGraph g = caps::chain_graph(c, seq, inv);
    caps::MolInfo m = caps::molecule_info(g);
    // end groups: each replaces an end hydrogen, so the chain gains (the group with * as H) − H₂
    for (const std::string& end : {c.head_cap, c.tail_cap}) {
      std::string sm = caps::chain_end_smiles(end);
      if (sm.empty()) continue;
      // the group with its * made a hydrogen: the * removed, its atom takes the implicit H ("*C(=O)O" → formic acid)
      for (size_t at; (at = sm.find("(*)")) != std::string::npos;) sm.erase(at, 3);
      sm.erase(std::remove(sm.begin(), sm.end(), '*'), sm.end());
      caps::MolGraph eg = caps::parse_smiles(sm);
      caps::add_hydrogens(eg);
      const caps::MolInfo e = caps::molecule_info(eg);
      auto n = formula_counts(m.formula);
      for (const auto& [el, k] : formula_counts(e.formula)) n[el] += k;
      n["H"] -= 2;
      m.formula = hill(n);
      m.mass += e.mass - 2 * caps::element(1).mass;
      m.atoms += e.atoms - 2;
    }
    caps::Json s = caps::Json::array(), r = caps::Json::array();
    for (int k : seq) s.push_back(double(k));
    for (char x : inv) r.push_back(double(x));
    j["ok"] = true;
    j["sequence"] = s;
    j["inverted"] = r;   // units written backwards (head-to-head and random linkage)
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
      if (sj.has("lookahead")) g.lookahead = std::clamp(int(sj["lookahead"].number()), 1, 4);
      if (sj.has("threads")) g.threads = std::max(0, int(sj["threads"].number()));
      // method: trials (best by contact margin) | rosenbluth (soft spheres) | rosenbluth_lj (UFF Lennard-Jones); temperature K
      const std::string meth = sj.text("method", "trials");
      g.method = meth == "rosenbluth" ? 1 : meth == "rosenbluth_lj" ? 2 : 0;
      g.method_temperature = sj.num("temperature", 450);
      // orientation: {axis: "x" | "y" | "z" | [x, y, z], strength: s in kT} (an aligning field on the backbone)
      if (sj.has("orientation") && sj["orientation"].is_object()) {
        const caps::Json& O = sj["orientation"];
        if (O.has("axis") && O["axis"].is_array() && O["axis"].size() == 3)
          for (size_t k = 0; k < 3; ++k) g.orient_axis[k] = O["axis"][k].number();
        else {
          const std::string ax = O.text("axis", "z");
          g.orient_axis = ax == "x" ? std::array<double, 3>{1, 0, 0} : ax == "y" ? std::array<double, 3>{0, 1, 0} : std::array<double, 3>{0, 0, 1};
        }
        g.orient_strength = std::max(0.0, O.num("strength", 0));
      }
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
  o.whole_molecules = j.num("whole_molecules", 0) != 0;
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

// Force fields by group (v36, ffmerge.hpp): each group's atoms assigned on their own (typing, parameters and charges as
// caps_field_assign does), then merged into one force field with the cross pairs by the chosen rule; the report is the
// groups' reports joined, atoms renumbered to the structure.
// A model's own force field (Kremer–Grest: FENE + WCA, built from the melt's options), as a complete assignment.
void field_run_file(caps_doc* d);

void field_run_model(caps_doc* d) {
  FieldState& F = *d->field;
  const caps::Json J = caps::Json::parse(F.model);
  if (J.text("model", "") == "file") { field_run_file(d); return; }
  if (J.text("model", "") != "kremer-grest") throw caps::FFError("unknown model " + J.text("model", ""));
  caps::KgOptions o;
  o.k_theta = J.num("k_theta", 0);
  o.density = J.num("density", 0.85);
  o.sigma = J.num("sigma", 0), o.temperature = J.num("temperature", 0), o.bead_mass = J.num("bead_mass", 0);
  const caps::System& s = d->frame;
  auto M = std::make_shared<caps::ForceField>(caps::kremer_grest_forcefield(s, o));
  F.ff = M;
  F.types = M->atom_type;
  F.complete = true;
  F.ff_path = "kremer-grest";
  F.base = caps::FFDef{};
  F.base.name = M->name;
  const size_t n = s.atoms.size();
  for (size_t i = 0; i < n; ++i) {
    d->traj.topology.atoms[i].type = 1;
    d->traj.topology.atoms[i].name = "KG";
    d->traj.topology.atoms[i].charge = 0;
  }
  d->traj.topology.types.clear();
  caps::TypeInfo ti;
  ti.type = 1, ti.label = "KG", ti.mass = n ? M->mass[0] : 1.0;
  d->traj.topology.types.push_back(ti);
  d->traj.topology.has_charges = true;
  refresh(d);
  caps::Json r = caps::Json::object();
  r["forcefield"] = M->name;
  r["mixing"] = std::string("one bead type");
  r["version"] = std::string("1990");
  r["source"] = std::string("Kremer & Grest, J. Chem. Phys. 92, 5057 (1990)");
  r["file"] = std::string("kremer-grest");
  r["typing"] = std::string("every bead of the melt");
  r["charges"] = std::string("none");
  caps::Json atoms = caps::Json::array();
  for (size_t i = 0; i < n; ++i) {
    caps::Json a = caps::Json::object();
    a["i"] = double(i + 1), a["el"] = std::string("C"), a["type"] = std::string("KG"), a["ov"] = false;
    a["rule"] = std::string("Kremer–Grest bead"), a["src"] = std::string("model"), a["q"] = 0.0, a["cands"] = caps::Json::array();
    atoms.push_back(a);
  }
  r["atoms"] = atoms;
  r["typed"] = double(n), r["untyped"] = 0.0, r["overridden"] = 0.0, r["ambiguous"] = 0.0, r["rules"] = 0.0;
  r["net_charge"] = 0.0, r["has_charges"] = true;
  r["missing"] = caps::Json::array();
  r["filled"] = 0.0, r["filled_terms"] = caps::Json::array(), r["estimated"] = 0.0, r["imported"] = 0.0;
  r["by_analogy"] = caps::Json::array(), r["entered"] = caps::Json::array(), r["imported_files"] = caps::Json::array();
  caps::Json refs = caps::Json::array();
  refs.push_back(std::string("Kremer & Grest, J. Chem. Phys. 92, 5057 (1990)"));
  r["references"] = refs;
  caps::Json used = caps::Json::array(), u = caps::Json::object();
  u["name"] = std::string("KG"), u["count"] = double(n), u["colour"] = hex_colour(caps::molecule_colour(0));
  used.push_back(u);
  r["used"] = used;
  r["fftypes"] = caps::Json::array();
  r["styles"] = caps::Json::object();
  caps::Json notes = caps::Json::array();
  for (const auto& x : M->notes) notes.push_back(x);
  notes.push_back(std::string("FENE bonds (K = 30 ε/σ², R₀ = 1.5 σ) with their WCA core, the WCA pair cut at 2^(1/6) σ and shifted") +
                  (o.k_theta > 0 ? ", cosine bending k_θ = " + std::to_string(o.k_theta) + " ε" : std::string(", no bending term")) + "; special_bonds fene");
  r["notes"] = notes;
  r["complete"] = true;
  {
    caps::Evaluator ev(*M, elec());
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

// The force field the file carried (an AMBER prmtop): every term as the file gives it, nothing typed or estimated
void field_run_file(caps_doc* d) {
  FieldState& F = *d->field;
  const caps::System& s = d->frame;
  const auto M = d->traj.topology.forcefield;
  const size_t n = s.atoms.size();
  if (!M || M->charge.size() != n) throw caps::FFError("the structure no longer holds the atoms its topology file describes: assign a force field from the library");
  F.ff = M;
  F.types = M->atom_type;
  F.complete = true;
  F.ff_path = "file";
  F.base = caps::FFDef{};
  F.base.name = M->name;
  for (size_t i = 0; i < n; ++i) d->traj.topology.atoms[i].charge = M->charge[i];
  d->traj.topology.has_charges = true;
  refresh(d);
  caps::Json r = caps::Json::object();
  r["forcefield"] = M->name;
  r["mixing"] = std::string("arithmetic (Lorentz–Berthelot)") + (M->pair_override.empty() ? "" : ", " + std::to_string(M->pair_override.size()) + " type pairs with their own coefficients");
  r["version"] = std::string("");
  r["source"] = std::string("the topology file");
  r["file"] = std::string("file");
  r["typing"] = std::string("types from the topology file");
  r["charges"] = std::string("from the topology file");
  caps::Json atoms = caps::Json::array();
  std::map<std::string, int> count;
  double q = 0;
  for (size_t i = 0; i < n; ++i) {
    caps::Json a = caps::Json::object();
    a["i"] = double(i + 1), a["el"] = std::string(caps::element(s.atoms[i].element).symbol), a["type"] = M->atom_type[i], a["ov"] = false;
    a["rule"] = std::string("from the topology file"), a["src"] = std::string("file"), a["q"] = M->charge[i], a["cands"] = caps::Json::array();
    atoms.push_back(a);
    ++count[M->atom_type[i]];
    q += M->charge[i];
  }
  r["atoms"] = atoms;
  r["typed"] = double(n), r["untyped"] = 0.0, r["overridden"] = 0.0, r["ambiguous"] = 0.0, r["rules"] = 0.0;
  r["net_charge"] = q, r["has_charges"] = true;
  r["missing"] = caps::Json::array();
  r["filled"] = 0.0, r["filled_terms"] = caps::Json::array(), r["estimated"] = 0.0, r["imported"] = 0.0;
  r["by_analogy"] = caps::Json::array(), r["entered"] = caps::Json::array(), r["imported_files"] = caps::Json::array();
  r["references"] = caps::Json::array();
  caps::Json used = caps::Json::array();
  int k = 0;
  for (const auto& [t, c] : count) {
    caps::Json u = caps::Json::object();
    u["name"] = t, u["count"] = double(c), u["colour"] = hex_colour(caps::molecule_colour(k++));
    used.push_back(u);
  }
  r["used"] = used;
  r["fftypes"] = caps::Json::array();
  caps::Json st = caps::Json::object();
  st["pair"] = M->native_pair, st["bond"] = std::string("harmonic"), st["angle"] = std::string("harmonic"), st["dihedral"] = M->native_dihedral;
  st["improper"] = M->native_improper, st["special"] = M->native_special;
  r["styles"] = st;
  caps::Json notes = caps::Json::array();
  for (const auto& x : M->notes) notes.push_back(x);
  r["notes"] = notes;
  r["complete"] = true;
  {
    caps::Evaluator ev(*M, elec());
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

// A document whose file carried its force field starts with it assigned
void install_file_field(caps_doc* d) {
  if (!d->traj.topology.forcefield || d->traj.topology.forcefield->charge.size() != d->traj.topology.atoms.size()) return;
  auto F = std::make_unique<FieldState>();
  for (const auto& a : d->traj.topology.atoms) F->file_types.push_back({a.type, a.name}), F->file_charges.push_back(a.charge);
  F->file_type_table = d->traj.topology.types;
  F->file_has_charges = d->traj.topology.has_charges;
  F->charges = "keep";
  F->model = "{\"model\":\"file\"}";
  d->field = std::move(F);
  field_run(d);
}

void field_run_groups(caps_doc* d) {
  FieldState& F = *d->field;
  const caps::Json G = caps::Json::parse(F.groups);
  const caps::System& s = d->frame;
  const size_t n = s.atoms.size();
  const auto molx = s.molecules();
  auto mol_of = [&](size_t i) { return s.has_mol ? s.atoms[i].mol : int64_t(molx[i]) + 1; };
  if (!G.has("groups") || !G["groups"].is_array() || G["groups"].size() == 0) throw caps::FFError("groups: [{name, molecules, forcefield, charges}]");
  // which atoms each group holds: molecule ids ("1", "2-10, 12") or the rest
  std::vector<int> owner(n, -1);
  const size_t ng = G["groups"].size();
  int rest = -1;
  for (size_t g = 0; g < ng; ++g) {
    const caps::Json& J = G["groups"][g];
    std::string m = J.text("molecules", "");
    if (m == "rest" || m == "*") { rest = int(g); continue; }
    if (m == "water") {   // every water molecule (an O with two H, and its M site)
      for (const auto& w : caps::find_waters(s))
        for (auto a : w)
          if (a >= 0) {
            if (owner[size_t(a)] >= 0) throw caps::FFError("a water is in two groups");
            owner[size_t(a)] = int(g);
          }
      continue;
    }
    for (auto& c : m) if (c == ',' || c == ';') c = ' ';
    std::istringstream is(m);
    std::set<int64_t> want;
    for (std::string w; is >> w;) {
      const auto dash = w.find('-', 1);
      const int64_t a = std::stoll(w.substr(0, dash)), b = dash == std::string::npos ? a : std::stoll(w.substr(dash + 1));
      for (int64_t k = a; k <= b && k - a < 10000000; ++k) want.insert(k);
    }
    for (size_t i = 0; i < n; ++i)
      if (want.count(mol_of(i))) {
        if (owner[i] >= 0) throw caps::FFError("molecule " + std::to_string(mol_of(i)) + " is in two groups");
        owner[i] = int(g);
      }
  }
  for (size_t i = 0; i < n; ++i)
    if (owner[i] < 0) {
      if (rest < 0) throw caps::FFError("atom " + std::to_string(i + 1) + " (molecule " + std::to_string(mol_of(i)) + ") is in no group: add it, or a group for the rest");
      owner[i] = rest;
    }
  F.group_atoms.clear();
  for (size_t g = 0; g < ng; ++g) {
    caps::LammpsStyle::Group grp{G["groups"][g].text("name", "group" + std::to_string(g + 1)), {}};
    for (size_t i = 0; i < n; ++i) if (owner[i] == int(g)) grp.atoms.push_back(uint32_t(i));
    F.group_atoms.push_back(std::move(grp));
  }
  std::vector<std::shared_ptr<const caps::ForceField>> keep;
  std::vector<caps::FFPart> parts;
  std::vector<caps::Json> reports;
  std::vector<std::string> names, paths;
  bool complete = true, based = false;
  // the force-field groups first: a potential group's cross Lennard-Jones takes their form (12-6, or class II's 9-6)
  std::vector<std::shared_ptr<const caps::ForceField>> slot_ff(ng);
  std::vector<std::vector<uint32_t>> slot_atoms(ng);
  std::vector<caps::Json> slot_rep(ng);
  std::vector<std::string> slot_name(ng), slot_path(ng);
  std::string cross_form = "lj12-6";
  for (int pass = 0; pass < 2; ++pass)
  for (size_t g = 0; g < ng; ++g) {
    const caps::Json& J = G["groups"][g];
    if ((J.has("potential") && J["potential"].is_object()) != (pass == 1)) continue;
    std::vector<uint32_t> atoms;
    for (size_t i = 0; i < n; ++i) if (owner[i] == int(g)) atoms.push_back(uint32_t(i));
    if (atoms.empty()) continue;
    const std::string name = J.text("name", "group " + std::to_string(g + 1));
    if (J.has("water")) {   // a water model: its own geometry, charges, Lennard-Jones and M sites (caps/water.hpp)
      const caps::WaterModel& wm = caps::water_model(J.text("water", ""));
      auto wf = std::make_shared<caps::ForceField>(caps::water_forcefield(s, wm, atoms));
      slot_ff[g] = wf, slot_atoms[g] = atoms;
      caps::Json R = caps::Json::object();
      R["forcefield"] = wf->name;
      caps::Json ra = caps::Json::array(), rn = caps::Json::array(), rr = caps::Json::array();
      for (size_t k = 0; k < atoms.size(); ++k) {
        caps::Json a = caps::Json::object();
        a["i"] = double(k + 1), a["el"] = wf->atom_type[k], a["type"] = wf->atom_type[k], a["ov"] = false, a["rule"] = wf->why[k];
        a["src"] = std::string("water model"), a["q"] = wf->charge[k], a["cands"] = caps::Json::array();
        ra.push_back(a);
      }
      for (const auto& x : wf->notes) rn.push_back(x);
      rr.push_back(wm.citation);
      R["atoms"] = ra, R["notes"] = rn, R["references"] = rr, R["typed"] = double(atoms.size()), R["complete"] = true;
      slot_rep[g] = R, slot_name[g] = name, slot_path[g] = "water:" + wm.id;
      if (!based) { based = true; }
      continue;
    }
    const bool potential = J.has("potential") && J["potential"].is_object();
    const std::string path = potential ? J["potential"].text("file", "") : J.text("forcefield", "");
    if (path.empty()) throw caps::FFError(name + (potential ? ": no potential file" : ": no force field"));
    // the group alone: its atoms, the bonds among them, the cell, the charges the structure has
    caps::System sub;
    sub.cell = s.cell;
    sub.has_charges = s.has_charges;
    sub.has_mol = s.has_mol;
    sub.bonds_from_file = true;
    std::vector<int64_t> at(n, -1);
    for (size_t k = 0; k < atoms.size(); ++k) at[atoms[k]] = int64_t(k), sub.atoms.push_back(s.atoms[atoms[k]]);
    for (const auto& b : s.bonds)
      if (at[b.i] >= 0 && at[b.j] >= 0) sub.bonds.push_back({uint32_t(at[b.i]), uint32_t(at[b.j]), b.order});
    if (potential) {   // a literature many-body potential (Tersoff, EAM …) read by LAMMPS from its file
      caps::ManyBodySpec spec;
      spec.style = J["potential"].text("style", "");
      spec.file = path;
      spec.units = J["potential"].text("units", "");
      spec.args = J["potential"].text("args", "");
      spec.file2 = J["potential"].text("file2", "");   // MEAM's parameter file
      if (J["potential"].has("entries") && J["potential"]["entries"].is_object())   // MEAM: element → library entry
        for (const auto& [k, v] : J["potential"]["entries"].members()) spec.entries.push_back({k, v.str()});
      {   // class II's 9-6 when every force-field group is 9-6 (PCFF, COMPASS), else 12-6
        bool all96 = false, any = false;
        for (size_t h = 0; h < ng; ++h)
          if (slot_ff[h] && !slot_ff[h]->manybody.on()) all96 = (any ? all96 : true) && slot_ff[h]->pair_form == "lj9-6", any = true;
        cross_form = any && all96 ? "lj9-6" : "lj12-6";
      }
      spec.pair_form = cross_form;
      std::vector<std::string> mb_notes;
      auto mf = std::make_shared<caps::ForceField>(caps::manybody_part(sub, spec, &mb_notes));
      slot_ff[g] = mf, slot_atoms[g] = atoms;
      caps::Json R = caps::Json::object();
      R["forcefield"] = mf->name;
      caps::Json ra = caps::Json::array(), rn = caps::Json::array(), rr = caps::Json::array();
      for (size_t k = 0; k < atoms.size(); ++k) {
        caps::Json a = caps::Json::object();
        a["i"] = double(k + 1);
        a["el"] = mf->atom_type[k];
        a["type"] = mf->atom_type[k];
        a["ov"] = false;
        a["rule"] = mf->why[k];
        a["src"] = std::string("potential file");
        a["q"] = 0.0;
        a["cands"] = caps::Json::array();
        ra.push_back(a);
      }
      for (const auto& x : mb_notes) rn.push_back(x);
      if (!mf->manybody.citation.empty()) rr.push_back(mf->manybody.citation);
      R["atoms"] = ra, R["notes"] = rn, R["references"] = rr;
      R["typed"] = double(atoms.size());
      R["complete"] = true;
      slot_rep[g] = R, slot_name[g] = name, slot_path[g] = path;
      continue;
    }
    std::unique_ptr<caps_doc> tmp(doc_of(sub));
    const int rc = caps_field_assign(tmp.get(), path.c_str(), nullptr, int32_t(J.num("charges", 4)));
    if (rc < 0) throw caps::FFError(name + ": " + g_error);
    if (d->hand_example && path == d->hand_ff && tmp->field) {   // the types taught by hand, learned onto this group
      std::vector<std::string> types = d->hand_types;
      caps::Json unknown = caps::Json::array();
      std::set<std::string> said;
      for (auto& t : types)
        if (!t.empty() && !tmp->field->base.type(t) && !tmp->field->extra.type(t)) {
          if (said.insert(t).second) unknown.push_back(t);
          t.clear();
        }
      const caps::ExampleTypes learned = caps::learn_types(*d->hand_example, types);
      const caps::ExampleMatch m = caps::apply_types(tmp->frame, learned);
      size_t set = 0;
      for (size_t i = 0; i < m.types.size(); ++i)
        if (!m.types[i].empty()) tmp->field->overrides[int32_t(i)] = m.types[i], ++set;
      field_run(tmp.get());
      caps::Json& H = d->hand_report;
      if (!H.is_object()) H = caps::Json::object();
      H["radius"] = double(learned.radius), H["environments"] = double(learned.environments);
      H["set"] = H.num("set", 0) + double(set), H["unmatched"] = H.num("unmatched", 0) + double(m.unmatched);
      caps::Json gl = H.has("groups") ? H["groups"] : caps::Json::array();
      gl.push_back(name);
      H["groups"] = gl;
      caps::Json c = caps::Json::array();
      for (const auto& x : learned.conflicts) c.push_back(x);
      H["conflicts"] = c, H["unknown_types"] = unknown;
    }
    if (!tmp->field->ff) throw caps::FFError(name + ": " + std::to_string(size_t(caps::Json::parse(tmp->field->report)["untyped"].number())) + " atoms untyped by " + tmp->field->base.name);
    complete = complete && tmp->field->complete;
    slot_ff[g] = tmp->field->ff, slot_atoms[g] = atoms;
    slot_rep[g] = caps::Json::parse(tmp->field->report), slot_name[g] = name, slot_path[g] = path;
    if (!based) F.base = tmp->field->base, F.ff_path = path, based = true;
  }
  if (!based) throw caps::FFError("every group has a potential file: give the rest of the system a force field");
  for (size_t g = 0; g < ng; ++g)   // the groups in their order
    if (slot_ff[g]) {
      keep.push_back(slot_ff[g]);
      parts.push_back({keep.back().get(), slot_atoms[g], slot_name[g]});
      reports.push_back(slot_rep[g]);
      names.push_back(slot_name[g]), paths.push_back(slot_path[g]);
    }
  caps::MergeOptions mo;
  mo.eps_rule = G.text("eps_rule", mo.eps_rule);
  mo.sigma_rule = G.text("sigma_rule", mo.sigma_rule);
  mo.scaling14 = G.text("scaling14", mo.scaling14);
  mo.cross96 = G.text("cross96", mo.cross96);
  mo.refit_cutoff = G.num("refit_cutoff", elec().cutoff > 0 ? elec().cutoff : mo.refit_cutoff);
  if (G.has("pairs") && G["pairs"].is_array())
    for (const auto& x : G["pairs"].items()) mo.explicit_pairs.push_back({x.text("a", ""), x.text("b", ""), x.num("eps", 0), x.num("sigma", 0)});
  std::vector<std::string> merge_notes;
  auto M = std::make_shared<caps::ForceField>(caps::merge_forcefields(n, parts, mo, &merge_notes));
  F.ff = M;
  F.types = M->atom_type;
  F.complete = complete;
  F.rep = caps::ParamReport{};
  // types into the document
  for (size_t i = 0; i < n; ++i) {
    d->traj.topology.atoms[i].type = M->type_index[i] + 1;
    d->traj.topology.atoms[i].name = M->atom_type[i];
    d->traj.topology.atoms[i].charge = M->charge[i];
  }
  d->traj.topology.types.clear();
  for (size_t t = 0; t < M->type_names.size(); ++t) {
    caps::TypeInfo ti;
    ti.type = int(t) + 1;
    ti.label = M->type_names[t];
    for (size_t i = 0; i < n; ++i) if (M->type_index[i] == int(t)) { ti.mass = M->mass[i]; break; }
    d->traj.topology.types.push_back(ti);
  }
  d->traj.topology.has_charges = true;
  refresh(d);
  // the report: the groups' joined
  caps::Json r = caps::Json::object();
  r["forcefield"] = M->name;
  // the rule actually used between groups (an "auto" rule resolved by the merge, said in its note)
  std::string between = "ε " + (mo.sigma_rule == "sixthpower" ? std::string("and σ sixth-power") : mo.eps_rule + ", σ " + mo.sigma_rule);
  for (const auto& nt : M->notes)
    if (auto p = nt.find("cross type pairs by "); p != std::string::npos) {
      between = nt.substr(p + 20);
      if (auto q = between.find(", written out explicitly"); q != std::string::npos) between.erase(q, 24);
      break;
    }
  r["mixing"] = "by group · within each group its own force field's rule · between groups " + between;
  r["version"] = std::string("");
  r["source"] = std::string("");
  r["file"] = F.ff_path;
  r["typing"] = std::string("each group by its own force field");
  r["charges"] = std::string("by group");
  caps::Json atoms_out = caps::Json::array(), miss = caps::Json::array(), refs = caps::Json::array(), notes = caps::Json::array(), fftypes = caps::Json::array(), groups = caps::Json::array();
  std::set<std::string> ref_seen;
  double typed = 0, untyped = 0, overridden = 0, ambiguous = 0, filled = 0, estimated = 0, imported = 0, rules = 0;
  std::vector<caps::Json> atom_rows(n);
  for (size_t g = 0; g < parts.size(); ++g) {
    const caps::Json& R = reports[g];
    const std::string tag = names[g];
    typed += R.num("typed", 0), untyped += R.num("untyped", 0), overridden += R.num("overridden", 0), ambiguous += R.num("ambiguous", 0);
    filled += R.num("filled", 0), estimated += R.num("estimated", 0), imported += R.num("imported", 0), rules += R.num("rules", 0);
    if (R.has("atoms"))
      for (const auto& a : R["atoms"].items()) {
        const size_t k = size_t(a.num("i", 1)) - 1;
        if (k >= parts[g].atoms.size()) continue;
        const uint32_t gi = parts[g].atoms[k];
        caps::Json x = a;
        x["i"] = double(gi + 1);
        x["type"] = M->atom_type[gi];
        x["src"] = tag + " · " + a.text("src", "");
        atom_rows[gi] = x;
      }
    for (const char* key : {"missing", "notes"})
      if (R.has(key)) for (const auto& m : R[key].items()) (std::string(key) == "missing" ? miss : notes).push_back(caps::Json(tag + ": " + m.str()));
    if (R.has("references")) for (const auto& x : R["references"].items()) if (ref_seen.insert(x.dump(0)).second) refs.push_back(x);
    if (R.has("fftypes")) for (const auto& x : R["fftypes"].items()) fftypes.push_back(x);
    caps::Json gj = caps::Json::object();
    gj["name"] = tag, gj["forcefield"] = R.text("forcefield", ""), gj["file"] = paths[g], gj["atoms"] = double(parts[g].atoms.size());
    gj["complete"] = R.has("complete") && R["complete"].kind() == caps::Json::Bool && R["complete"].boolean();
    groups.push_back(gj);
  }
  for (auto& a : atom_rows) atoms_out.push_back(a);
  for (const auto& m : merge_notes) notes.push_back(caps::Json(m));
  double qsum = 0;
  for (double q : M->charge) qsum += q;
  r["atoms"] = atoms_out;
  r["typed"] = typed, r["untyped"] = untyped, r["overridden"] = overridden, r["ambiguous"] = ambiguous, r["rules"] = rules;
  r["net_charge"] = qsum;
  r["has_charges"] = true;
  r["missing"] = miss;
  r["filled"] = filled, r["filled_terms"] = caps::Json::array(), r["estimated"] = estimated, r["imported"] = imported;
  r["by_analogy"] = caps::Json::array(), r["entered"] = caps::Json::array(), r["imported_files"] = caps::Json::array();
  r["references"] = refs;
  caps::Json used = caps::Json::array();
  std::map<int, int> count;
  for (int t : M->type_index) ++count[t];
  for (size_t t = 0; t < M->type_names.size(); ++t) {
    caps::Json u = caps::Json::object();
    u["name"] = M->type_names[t];
    u["count"] = double(count[int(t)]);
    u["colour"] = hex_colour(caps::molecule_colour(int(t)));
    used.push_back(u);
  }
  r["used"] = used;
  r["fftypes"] = fftypes;
  r["styles"] = reports.empty() || !reports[0].has("styles") ? caps::Json::object() : reports[0]["styles"];
  r["notes"] = notes;
  r["groups"] = groups;
  r["complete"] = F.complete;
  {
    caps::Evaluator ev(*M, elec());
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
    caps::BrushReport brep;
    const bool brush = j.has("brush") && j["brush"].is_object();
    caps::System s;
    if (brush) {   // v62: chains grafted by one end to surface sites (a brush) instead of a free film
      const auto& b = j["brush"];
      caps::BrushOptions bo;
      bo.density = b.num("density", 0.3);
      bo.chains = int(b.num("chains", 0));
      bo.site = b.text("site", "O");
      bo.min_spacing = b.num("min_spacing", 4.0);
      bo.film = f.num("thickness", 40);
      bo.vacuum = std::max(5.0, f.num("vacuum", 20));
      bo.grow = io.grow;
      s = caps::build_brush(slab, c, bo, &brep);
    } else {
      s = caps::build_interface(slab, c, io, &rep);
    }
    std::string t;
    for (const auto& n : sr.notes) t += n + "\n";
    for (const auto& n : brush ? brep.notes : rep.notes) t += n + "\n";
    report_out(t, report, cap);
    caps_doc* d = doc_of(s);
    {
      caps::KeyValues pr = json_params(options_json);
      for (auto& kv : json_params(spec_json)) pr.push_back({"chain " + kv.first, kv.second});
      if (brush) prov_step(d, "brush.build", "polymer chains grafted to surface sites", std::move(pr), seeded(o ? o->seed : 0), {});
      else prov_step(d, "interface.build", "polymer grown against a surface", std::move(pr), seeded(o ? o->seed : 0), {"matsumoto1998"});
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
    t.material = j.text("material", "graphene");
    t.cc = j.num("cc", 0);
    t.walls = int(j.num("walls", 1));
    if (j.has("wall_chiralities") && j["wall_chiralities"].is_array())   // v62: each wall its own [n, m], innermost first
      for (const auto& w : j["wall_chiralities"].items())
        if (w.is_array() && w.size() == 2) t.wall_chiralities.push_back({int(w[0].number()), int(w[1].number())});
    t.bundle = int(j.num("bundle", 1));                                  // v62: ropes of 7, 19, 37 tubes
    t.bundle_lattice = j.num("bundle_lattice", 0) != 0;                 //      or the periodic triangular lattice
    t.tube_gap = j.num("tube_gap", 3.4);
    f = caps::nanotube(t, &r);
    keep = {t.bundle_lattice, t.bundle_lattice, t.periodic};
  } else if (kind == "sheet") {
    caps::SheetOptions sh;
    sh.lx = j.num("lx", 20), sh.ly = j.num("ly", 20);
    sh.layers = int(j.num("layers", 1));
    sh.periodic = j.num("periodic", 1) != 0;
    sh.material = j.text("material", "graphene");
    sh.cc = j.num("cc", 0);
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
    po.height = j.num("height", 0);
    po.top_ratio = j.num("top_ratio", 0.5);
    f = caps::nanoparticle(caps::read_cif(j.text("crystal")), po, &r);
    keep = {false, false, po.shape == caps::ParticleShape::Fibre};
    // the cut surface relaxed with UFF (design/boards/Nanoparticle): the surface layer (atoms within 4 Å of a hydrogen,
    // a ligand or an under-coordinated atom) free, the core held as the crystal made it. Not for metals: UFF is no model
    // of a metallic surface (EAM, in LAMMPS, is)
    if (j.num("relax_surface", 0) != 0) {
      const int metals[] = {26, 27, 28, 29, 44, 45, 46, 47, 76, 77, 78, 79};
      int heavy = 0, metal = 0;
      for (const auto& a : f.atoms)
        if (a.element > 1) ++heavy, metal += std::find(std::begin(metals), std::end(metals), a.element) != std::end(metals);
      if (metal * 2 > heavy)
        throw std::invalid_argument("a metal particle: UFF is no model of a metal surface — relax it with EAM in LAMMPS (Export) instead");
      const auto nb = f.neighbours();
      std::map<int, size_t> full;   // the most bonds each element has in the particle: its bulk coordination
      for (size_t i = 0; i < f.atoms.size(); ++i) full[f.atoms[i].element] = std::max(full[f.atoms[i].element], nb[i].size());
      std::vector<caps::Vec3> surf;
      for (size_t i = 0; i < f.atoms.size(); ++i) {
        const int el = f.atoms[i].element;
        bool s_atom = el == 1 || nb[i].size() < full[el] || f.atoms[i].mol != f.atoms[0].mol;
        for (uint32_t w : nb[i]) s_atom = s_atom || f.atoms[w].element == 1;
        if (s_atom) surf.push_back(f.atoms[i].pos);
      }
      caps::RelaxOptions ro;
      ro.fixed.assign(f.atoms.size(), 1);
      int free = 0;
      for (size_t i = 0; i < f.atoms.size(); ++i) {
        for (const auto& p : surf)
          if (caps::norm(f.atoms[i].pos - p) < 4.0) { ro.fixed[i] = 0; ++free; break; }
      }
      caps::UffOptions uo;
      ro.field = std::make_shared<caps::ForceField>(caps::assign_uff(f, uo));
      ro.pushoff = true;   // the cut leaves strained, close contacts at the surface: capped forces first
      ro.ftol = 1.0;
      ro.max_iterations = 20000;
      caps::RelaxReport rr;
      caps::relax(f, ro, &rr);
      char b[240];
      std::snprintf(b, sizeof b, "surface relaxed with UFF (no charges): %d surface atoms free, %zu core atoms held · E %.1f → %.1f kcal/mol · |F|max on the free atoms %.2f kcal/mol/Å",
                    free, f.atoms.size() - size_t(free), rr.initial.total(), rr.final.total(), rr.fmax_final);
      r.notes.push_back(b);
    }
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

// A polymer matrix grown around the document's structure (a filler built and functionalised in the Studio): the
// structure held at the centre, the directions it spans across its cell kept periodic.
extern "C" caps_doc* caps_embed_document(caps_doc* filler, const char* options_json, const char* spec_json, const caps_grow_opts* o, caps_progress_fn progress,
                                         void* user, char* report, int32_t cap) {
  try {
    if (!filler) throw std::invalid_argument("no document");
    const caps::Json j = caps::Json::parse(options_json && *options_json ? options_json : "{}");
    caps::System f = filler->traj.frame(filler->current);
    // the directions the filler fills across its cell stay periodic (a tube's axis, a sheet's plane): along each, the
    // largest empty gap between atoms (across the boundary too) under 2.5 Å; the others open for the matrix. Groups grafted
    // on a stacked sheet can close every gap: then the direction with the largest gap (the sheet's normal) opens
    std::array<bool, 3> keep{false, false, false};
    if (f.cell.valid()) {
      std::array<double, 3> gap{0, 0, 0};
      const caps::Vec3 edge[3] = {f.cell.a, f.cell.b, f.cell.c};
      for (int k = 0; k < 3; ++k) {
        if (!f.cell.periodic[size_t(k)]) { gap[size_t(k)] = 1e30; continue; }
        const double L = caps::norm(edge[k]);
        std::vector<double> u;
        for (const auto& a : f.atoms) {
          double x = f.cell.to_fractional(a.pos)[size_t(k)];
          x -= std::floor(x);
          u.push_back(x * L);
        }
        std::sort(u.begin(), u.end());
        double g = u.empty() ? L : L - u.back() + u.front();
        for (size_t i = 1; i < u.size(); ++i) g = std::max(g, u[i] - u[i - 1]);
        gap[size_t(k)] = g;
        keep[size_t(k)] = g < 2.5;
      }
      if (keep[0] && keep[1] && keep[2]) keep[size_t(std::max_element(gap.begin(), gap.end()) - gap.begin())] = false;
    }
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
    std::string notes;
    for (const auto& n : fr.notes) notes += n + "\n";
    report_out(notes, report, cap);
    caps_doc* d = doc_of(s);
    d->prov = filler->prov;   // the filler's history (built, functionalised) comes along
    caps::KeyValues pr = json_params(options_json);
    for (auto& kv : json_params(spec_json)) pr.push_back({"chain " + kv.first, kv.second});
    prov_step(d, "nano.embed", "the structure in a grown polymer matrix", std::move(pr), seeded(o ? o->seed : 0), {"matsumoto1998"});
    d->held_mol = 1;
    return d;
  } catch (const std::exception& e) {
    g_error = e.what();
    return nullptr;
  }
}

// Layer stacks (ABI 35, layers.hpp)
extern "C" caps_doc* caps_stack_documents(caps_doc* const* docs, int32_t n, const char* options_json, char* report, int32_t cap) {
  caps::Json r = caps::Json::object();
  try {
    if (!docs || n < 2) throw std::invalid_argument("a stack needs two documents or more");
    const caps::Json j = caps::Json::parse(options_json && *options_json ? options_json : "{}");
    std::vector<caps::System> frames;
    frames.reserve(size_t(n));
    for (int k = 0; k < n; ++k) {
      if (!docs[k]) throw std::invalid_argument("no document for layer " + std::to_string(k + 1));
      frames.push_back(docs[k]->traj.frame(docs[k]->current));
    }
    std::vector<caps::StackLayerInput> in;
    for (int k = 0; k < n; ++k) {
      std::string name = "layer " + std::to_string(k + 1);
      if (j.has("names") && j["names"].is_array() && size_t(k) < j["names"].size() && j["names"][size_t(k)].is_string()) name = j["names"][size_t(k)].str();
      caps::StackLayerInput li{name, &frames[size_t(k)]};
      // v62: per layer, flips: [bool] (180° about x) and shifts: [[dx, dy]] (Å in the plane)
      if (j.has("flips") && j["flips"].is_array() && size_t(k) < j["flips"].size()) li.flip = j["flips"][size_t(k)].kind() == caps::Json::Bool && j["flips"][size_t(k)].boolean();
      if (j.has("shifts") && j["shifts"].is_array() && size_t(k) < j["shifts"].size() && j["shifts"][size_t(k)].is_array() && j["shifts"][size_t(k)].size() == 2)
        li.shift_x = j["shifts"][size_t(k)][0].number(), li.shift_y = j["shifts"][size_t(k)][1].number();
      in.push_back(li);
    }
    caps::StackOptions so;
    so.gap = j.num("gap", so.gap);
    so.vacuum = j.num("vacuum", so.vacuum);
    so.match = j.text("match", so.match);
    so.max_repeat = int(j.num("max_repeat", so.max_repeat));
    caps::StackReport rep;
    const caps::System s = caps::stack_layers(in, so, &rep);
    caps_doc* d = doc_of(s);
    d->prov = docs[0]->prov;   // the first layer's history (a cleaved slab) comes along
    caps::KeyValues pr = json_params(options_json);
    std::string names;
    for (const auto& l : rep.layers) names += (names.empty() ? "" : " / ") + l.name;
    pr.push_back({"layers", names});
    prov_step(d, "build.stack", std::to_string(n) + " layers stacked along z", std::move(pr), "", {});
    int nm = 0;
    const auto mol = frames[0].molecules(&nm);
    (void)mol;
    if (nm == 1) d->held_mol = 1;
    r["ok"] = true;
    r["a"] = rep.a, r["b"] = rep.b, r["c"] = rep.c;
    caps::Json L = caps::Json::array(), N = caps::Json::array();
    for (const auto& l : rep.layers) {
      caps::Json o = caps::Json::object();
      o["name"] = l.name, o["na"] = double(l.na), o["nb"] = double(l.nb), o["strain_a"] = l.strain_a, o["strain_b"] = l.strain_b;
      o["z_lo"] = l.z_lo, o["z_hi"] = l.z_hi, o["atoms"] = double(l.atoms);
      L.push_back(o);
    }
    for (const auto& x : rep.notes) N.push_back(caps::Json(x));
    r["layers"] = L, r["notes"] = N;
    report_out(r.dump(0), report, cap);
    return d;
  } catch (const std::exception& e) {
    g_error = e.what();
    r["ok"] = false, r["error"] = std::string(e.what());
    report_out(r.dump(0), report, cap);
    return nullptr;
  }
}

extern "C" caps_doc* caps_frame_copy(caps_doc* d) {
  try {
    if (!d) throw std::invalid_argument("no document");
    caps_doc* c = doc_of(d->traj.frame(d->current));
    c->prov = d->prov;
    c->held_mol = d->held_mol;
    c->fixed_atoms = d->fixed_atoms;
    c->fixed_axes = d->fixed_axes;
    c->rigid_mols = d->rigid_mols;
    return c;
  } catch (const std::exception& e) {
    g_error = e.what();
    return nullptr;
  }
}

// DPD (ABI 34, dpd.hpp): a new document with the run's frames (beads as atoms) and a JSON report.
extern "C" caps_doc* caps_dpd(const char* json, caps_stage_fn progress, void* user, char* out, int32_t cap) {
  try {
    const caps::Json j = caps::Json::parse(json && *json ? json : "{}");
    caps::DpdOptions o;
    if (j.has("species") && j["species"].is_array())
      for (const auto& sp : j["species"].items()) o.species.push_back({sp.text("name", "molecule"), sp.text("sequence", "A"), int(sp.num("count", 0))});
    o.density = j.num("density", 3.0);
    if (j.has("chi") && j["chi"].is_object())
      for (const auto& [k, v] : j["chi"].members()) o.chi[k] = v.number();
    if (j.has("a") && j["a"].is_object())
      for (const auto& [k, v] : j["a"].members()) o.a[k] = v.number();
    o.gamma = j.num("gamma", 4.5);
    o.dt = j.num("dt", 0.04);
    o.bond_k = j.num("bond_k", 4.0);
    o.angle_k = j.num("angle_k", 0.0);   // v62: chain stiffness k_θ (1 + cos θ)
    o.start = j.text("start", "random");   // v62: random | lamellar | cylinders | spheres
    o.start_periods = int(j.num("start_periods", 1));
    o.steps = long(j.num("steps", 20000));
    o.equilibration = long(j.num("equilibration", double(o.steps) / 4));
    o.frame_every = int(j.num("frame_every", std::max(1.0, double(o.steps) / 40)));
    o.rc_angstrom = j.num("rc_angstrom", 6.46);
    o.seed = uint64_t(j.num("seed", 1));
    bool cancelled = false;
    if (progress)
      o.progress = [&](long step, double kT) {
        char b[80];
        std::snprintf(b, sizeof b, "step %ld · kT %.3f", step, kT);
        cancelled = progress(b, double(step) / double(o.steps), user) != 0;
        return !cancelled;
      };
    const auto rep = caps::run_dpd(o);
    if (cancelled) throw std::runtime_error("DPD cancelled");
    auto* d = new caps_doc;
    d->traj = rep.frames;
    d->current = d->traj.frames() - 1;
    refresh(d);
    caps::KeyValues pr = {{"beads", std::to_string(rep.beads) + " in " + std::to_string(rep.molecules) + " molecules, box " + g6(rep.box) + " r_c"},
                          {"steps", std::to_string(o.steps) + " of " + g6(o.dt) + " (equilibration " + std::to_string(o.equilibration) + ")"},
                          {"result", "ψ " + g6(rep.order) + " · domain spacing " + g6(rep.spacing) + " r_c"}};
    prov_step(d, "dpd.run", "Dissipative particle dynamics (Groot–Warren)", std::move(pr), seeded(o.seed), {"groot1997", "hoogerbrugge1992"});
    caps::Json r = caps::Json::object();
    r["ok"] = true;
    r["beads"] = double(rep.beads), r["molecules"] = double(rep.molecules), r["box"] = rep.box;
    r["kT"] = rep.kT, r["kT_error"] = rep.kT_error, r["pressure"] = rep.pressure, r["pressure_error"] = rep.pressure_error;
    r["order"] = rep.order, r["q_peak"] = rep.q_peak, r["spacing"] = rep.spacing, r["seconds"] = rep.seconds;
    r["domains_a"] = double(rep.domains_a), r["domains_b"] = double(rep.domains_b), r["domain_a_size"] = rep.domain_a_size, r["domain_b_size"] = rep.domain_b_size;
    r["largest_a"] = rep.largest_a, r["largest_b"] = rep.largest_b;
    caps::Json q = caps::Json::array(), sq = caps::Json::array(), os = caps::Json::array();
    for (double x : rep.q) q.push_back(x);
    for (double x : rep.sq) sq.push_back(x);
    for (const auto& [st, v] : rep.order_series) { caps::Json p = caps::Json::array(); p.push_back(st); p.push_back(v); os.push_back(std::move(p)); }
    r["q"] = std::move(q), r["sq"] = std::move(sq), r["order_series"] = std::move(os);
    r["types"] = rep.types;
    caps::Json nt = caps::Json::array();
    for (const auto& x : rep.notes) nt.push_back(x);
    r["notes"] = std::move(nt);
    report_out(r.dump(0), out, cap);
    return d;
  } catch (const std::exception& e) {
    g_error = e.what();
    caps::Json r = caps::Json::object();
    r["ok"] = false;
    r["error"] = std::string(e.what());
    report_out(r.dump(0), out, cap);
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
    // growth as Grow's: method trials | rosenbluth | rosenbluth_lj at a temperature (K), trial directions, look-ahead
    {
      const std::string meth = j.text("method", "trials");
      if (meth != "trials" && meth != "rosenbluth" && meth != "rosenbluth_lj") throw std::invalid_argument("blend method: trials, rosenbluth or rosenbluth_lj");
      bo.grow.method = meth == "rosenbluth" ? 1 : meth == "rosenbluth_lj" ? 2 : 0;
      bo.grow.method_temperature = j.num("temperature", 450);
      if (j.has("trials")) bo.grow.trials = std::clamp(int(j["trials"].number()), 4, 5000);
      if (j.has("lookahead")) bo.grow.lookahead = std::clamp(int(j["lookahead"].number()), 1, 4);
    }
    if (progress) bo.grow.progress = [&](int done, int total, int restarts) { return progress(done, total, restarts, user) == 0; };
    caps::BlendReport br;
    const caps::System s = caps::grow_blend(comps, bo, &br);
    std::string t;
    for (const auto& n : br.notes) t += n + "\n";
    // each component's molecules (for a force field per component: Field · by group)
    for (size_t k = 0; k < br.molecules.size(); ++k)
      t += "component " + std::to_string(k + 1) + " · molecules " + std::to_string(br.molecules[k].first) + "-" + std::to_string(br.molecules[k].second) + "\n";
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
    // the chains an earlier reaction run left carry on: the inserted molecules get ids of their own after them
    const bool chains = d->react_chains.size() == host.atoms.size() && d->react_bonds == host.bonds.size();
    if (chains && s.atoms.size() >= host.atoms.size()) {
      int64_t top = 0;
      for (int64_t v : d->react_chains) top = std::max(top, v);
      const size_t per = std::max<size_t>(1, guest.atoms.size());
      for (size_t i = host.atoms.size(); i < s.atoms.size(); ++i) d->react_chains.push_back(top + 1 + int64_t((i - host.atoms.size()) / per));
      d->react_bonds = s.bonds.size();
    } else {
      d->react_chains.clear();
    }
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
    // entries of other steps go (a branch's trunk stays while only the branch or the steps above it change)
    if (!d->pipeline) d->pcache.clear();
    else {
      std::set<uint64_t> live{steps_key(*d->pipeline, caps::pipeline_trunk(*d->pipeline))};
      std::set<std::string> branches{""};
      for (const auto& st : d->pipeline->steps) branches.insert(caps::step_branch(st));
      for (const auto& b : branches) {   // every branch of these steps (switching back shows it at once)
        caps::Pipeline q = *d->pipeline;
        q.branch = b;
        live.insert(steps_key(q, 0));
      }
      for (auto it = d->pcache.begin(); it != d->pcache.end();) it = live.count(it->key) ? it + 1 : d->pcache.erase(it);
    }
    run_doc_pipeline(d);
    return 0;
  });
}

extern "C" int32_t caps_pipeline_result(caps_doc* d, char* json, int32_t cap) {
  try {
    if (!d->pstate) return report_out("", json, cap);
    caps::Json j = caps::pipeline_result_json(*d->pstate);
    caps::Json c = caps::Json::object();
    c["hit"] = d->pstate_cached;
    c["frames"] = double(d->pcache.size());
    j["cache"] = c;
    return report_out(j.dump(0), json, cap);
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
  } else if (fmt == "poscar") {
    const auto held = use_pipeline ? std::vector<char>{} : fixed_mask(d, s);
    caps::write_poscar(s, path, held);
    notes.push_back("VASP 5 POSCAR: species and counts lines, Direct (fractional) coordinates in the cell");
    notes.push_back("atoms grouped by element in the order each first appears (VASP needs them grouped); no bonds, charges or types");
    if (std::any_of(held.begin(), held.end(), [](char f) { return f != 0; })) notes.push_back("Selective dynamics: F on the held coordinates, T on the rest");
    notes.push_back("name it POSCAR for VASP; the POTCAR must list the same species in the same order");
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

extern "C" int32_t caps_pipeline_write_outputs(caps_doc* d, const char* pipeline_json, const char* dir, char* report, int32_t cap) {
  caps::Json r = caps::Json::object();
  try {
    if (!d) throw std::invalid_argument("no document");
    const caps::Pipeline p = caps::pipeline_from_json(caps::Json::parse(pipeline_json && *pipeline_json ? pipeline_json : "{\"steps\": []}"));
    const int64_t ts = d->current < d->traj.timesteps.size() ? d->traj.timesteps[d->current] : 0;
    const caps::PipelineState st = caps::run_pipeline(d->frame, p, int(d->current), ts, &d->traj);
    caps::Json lines = caps::Json::array();
    for (const auto& l : caps::write_pipeline_outputs(st, p, dir ? dir : ".")) lines.push_back(caps::Json(l));
    r["ok"] = true;
    r["lines"] = lines;
    report_out(r.dump(0), report, cap);
    return 0;
  } catch (const std::exception& e) {
    g_error = e.what();
    r["ok"] = false, r["error"] = std::string(e.what());
    report_out(r.dump(0), report, cap);
    return -1;
  }
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

extern "C" int32_t caps_energy_terms(caps_doc* d, char* json, int32_t cap) {
  try {
    caps::ForceField ff;
    if (d->field && d->field->complete) ff = *d->field->ff;
    else ff = default_ff(d->frame);
    const caps::EnergyOptions eo = elec();
    caps::Evaluator ev(ff, eo);
    std::vector<double> x, f;
    for (const auto& a : d->frame.atoms) x.insert(x.end(), a.pos.begin(), a.pos.end());
    const caps::EnergyTerms e = ev.compute(x, d->frame.cell, f);
    caps::Json en = caps::Json::object();
    en["bond"] = e.bond, en["angle"] = e.angle, en["dihedral"] = e.dihedral, en["improper"] = e.improper;
    en["vdw"] = e.vdw, en["coulomb"] = e.coulomb, en["total"] = e.total();
    en["electrostatics"] = std::string(eo.electrostatics == caps::EnergyOptions::Electrostatics::PME ? "pme" : "dsf");
    en["forcefield"] = ff.name;
    return report_out(en.dump(0), json, cap);
  } catch (const std::exception& e) {
    g_error = e.what();
    return -1;
  }
}

extern "C" int32_t caps_set_fixed_atoms(caps_doc* d, const int32_t* atoms, int32_t n) {
  return guard([&] {
    d->fixed_atoms.clear();
    const size_t na = d->frame.atoms.size();
    for (int32_t k = 0; k < n && atoms; ++k)
      if (atoms[k] >= 0 && size_t(atoms[k]) < na) d->fixed_atoms.push_back(uint32_t(atoms[k]));
    std::sort(d->fixed_atoms.begin(), d->fixed_atoms.end());
    d->fixed_atoms.erase(std::unique(d->fixed_atoms.begin(), d->fixed_atoms.end()), d->fixed_atoms.end());
    return int32_t(d->fixed_atoms.size());
  });
}

// v62: molecules written as rigid bodies in the LAMMPS inputs ("1-3,7"; "" none). Returns how many.
extern "C" int32_t caps_set_rigid_molecules(caps_doc* d, const char* ranges) {
  return guard([&] {
    d->rigid_mols = ranges && *ranges ? parse_mol_ranges(ranges) : std::vector<int64_t>{};
    std::sort(d->rigid_mols.begin(), d->rigid_mols.end());
    d->rigid_mols.erase(std::unique(d->rigid_mols.begin(), d->rigid_mols.end()), d->rigid_mols.end());
    return int32_t(d->rigid_mols.size());
  });
}

// v62 analogs (R groups): JSON {core, groups: [{r, subs: [smiles …]}], max} → {analogs: [{smiles, name}]}.
extern "C" int32_t caps_enumerate_analogs(const char* json, char* out, int32_t cap) {
  caps::Json r = caps::Json::object();
  const int32_t rc = guard([&] {
    const auto j = caps::Json::parse(json ? json : "{}");
    std::vector<std::pair<int, std::vector<std::string>>> groups;
    if (j.has("groups"))
      for (const auto& g : j["groups"].items()) {
        std::vector<std::string> subs;
        if (g.has("subs")) for (const auto& x : g["subs"].items()) subs.push_back(x.str());
        groups.push_back({int(g.num("r", 1)), subs});
      }
    const auto an = caps::enumerate_analogs(j.text("core", ""), groups, size_t(std::clamp(j.num("max", 200), 1.0, 5000.0)));
    caps::Json a = caps::Json::array();
    for (const auto& x : an) {
      caps::Json o = caps::Json::object();
      o["smiles"] = x.smiles, o["name"] = x.name;
      a.push_back(std::move(o));
    }
    r["analogs"] = std::move(a);
    return 0;
  });
  if (rc < 0) return -1;
  return report_out(r.dump(0), out, cap);
}

// v62: every frame's timestep (as the file gives it); returns the frame count.
extern "C" int32_t caps_frame_timesteps(caps_doc* d, int64_t* out, int32_t cap) {
  const auto& t = d->traj.timesteps;
  const int32_t n = int32_t(d->traj.frames());
  for (int32_t k = 0; k < cap && k < n; ++k) out[k] = size_t(k) < t.size() ? t[size_t(k)] : int64_t(k);
  return n;
}

// v62: the rigid molecules as ranges ("1-3,7"), "" none.
extern "C" int32_t caps_rigid_molecules(caps_doc* d, char* out, int32_t cap) {
  std::string t;
  const auto& m = d->rigid_mols;
  for (size_t k = 0; k < m.size();) {
    size_t e = k;
    while (e + 1 < m.size() && m[e + 1] == m[e] + 1) ++e;
    if (!t.empty()) t += ",";
    t += std::to_string(m[k]) + (e > k ? "-" + std::to_string(m[e]) : "");
    k = e + 1;
  }
  return report_out(t, out, cap);
}

// v62: which coordinates of the fixed atoms are held: bits x 1, y 2, z 4 (7: all, the default). Returns the axes set.
extern "C" int32_t caps_set_fixed_axes(caps_doc* d, int32_t axes) {
  d->fixed_axes = (axes & 7) ? (axes & 7) : 7;
  return d->fixed_axes;
}
extern "C" int32_t caps_fixed_axes(const caps_doc* d) { return d->fixed_axes; }

extern "C" int32_t caps_fixed_atoms(const caps_doc* d, int32_t* atoms, int32_t cap) {
  if (!d) return 0;
  for (int32_t k = 0; k < cap && size_t(k) < d->fixed_atoms.size(); ++k) atoms[k] = int32_t(d->fixed_atoms[size_t(k)]);
  return int32_t(d->fixed_atoms.size());
}

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
    L.column.clear();
    if (colour.rfind("column:", 0) == 0) {
      L.column = colour.substr(7);
      if (!d->traj.columns.count(L.column)) throw std::invalid_argument("no per-atom column '" + L.column + "' in this trajectory");
      L.colour = 5;
    }
    L.preview_q.clear();
    if (j.has("charges") && j["charges"].is_array())
      for (const auto& x : j["charges"].items()) L.preview_q.push_back(x.number());
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
  if (L.preview_q.size() == d->frame.atoms.size()) for (double x : L.preview_q) qlo = std::min(qlo, x), qhi = std::max(qhi, x);
  else for (const auto& a : d->frame.atoms) qlo = std::min(qlo, a.charge), qhi = std::max(qhi, a.charge);
  j["preview"] = L.preview_q.size() == d->frame.atoms.size();
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
  try {
    if (k == "rs") {
      for (const auto& x : caps::stereo_labels(f)) arr.push_back(x);
    } else if (k == "ez") {
      for (const auto& x : caps::ez_labels(f)) arr.push_back(x);
    } else {
      // the assignment's types when the Field has typed this structure
      const std::vector<std::string>* types = d->field && d->field->types.size() == f.atoms.size() ? &d->field->types : nullptr;
      for (auto& x : caps::atom_labels(f, k, types)) arr.push_back(std::move(x));
    }
  } catch (const std::exception& e) {
    g_error = e.what();
    return -1;
  }
  return report_out(arr.dump(0), json, cap);
}

extern "C" int32_t caps_bond_labels(caps_doc* d, const char* kind, char* json, int32_t cap) {
  caps::Json j = caps::Json::object(), pairs = caps::Json::array(), labels = caps::Json::array(), crossing = caps::Json::array();
  try {
    const caps::System& f = d->frame;
    const bool typed = d->field && d->field->types.size() == f.atoms.size();
    const auto ls = caps::bond_labels(f, kind ? kind : "length", typed ? &d->field->types : nullptr, typed && d->field->ff ? d->field->ff.get() : nullptr);
    for (const auto& l : ls) {
      pairs.push_back(double(l.i));
      pairs.push_back(double(l.j));
      labels.push_back(l.text);
      crossing.push_back(l.crossing);
    }
  } catch (const std::exception& e) {
    g_error = e.what();
    return -1;
  }
  j["pairs"] = pairs;
  j["labels"] = labels;
  j["crossing"] = crossing;
  return report_out(j.dump(0), json, cap);
}

extern "C" int32_t caps_type_table(caps_doc* d, char* json, int32_t cap) {
  caps::Json arr = caps::Json::array();
  const caps::System& f = d->frame;
  std::map<int, std::map<int, size_t>> el;   // type → element → atoms
  for (const auto& a : f.atoms) ++el[a.type][a.element];
  for (const auto& [t, m] : el) {
    caps::Json r = caps::Json::object();
    r["type"] = double(t);
    std::string label;
    double mass = 0;
    for (const auto& ti : f.types) if (ti.type == t) label = ti.label, mass = ti.mass;
    r["label"] = label;
    r["mass"] = mass;
    size_t n = 0;
    int best = 0;
    size_t bestn = 0;
    caps::Json els = caps::Json::array();
    for (const auto& [z, c] : m) {
      n += c;
      if (c > bestn) best = z, bestn = c;
      els.push_back(std::string(caps::element(z).symbol));
    }
    r["count"] = double(n);
    r["element"] = std::string(caps::element(best).symbol);
    r["elements"] = els;
    // the element the mass points to, to show a mismatch
    r["from_mass"] = mass > 0 ? std::string(caps::element(caps::element_from_mass(mass, 0.5)).symbol) : std::string();
    arr.push_back(r);
  }
  return report_out(arr.dump(0), json, cap);
}

extern "C" int32_t caps_repeat_unit_smiles(caps_doc* d, int32_t head, int32_t tail, char* json, int32_t cap) {
  caps::Json j = caps::Json::object();
  try {
    if (head < 0 || tail < 0) throw std::invalid_argument("pick the head atom, then the tail atom");
    j["smiles"] = caps::repeat_unit_smiles(d->frame, uint32_t(head), uint32_t(tail));
    j["ok"] = true;
  } catch (const std::exception& e) {
    j["ok"] = false;
    j["error"] = std::string(e.what());
  }
  return report_out(j.dump(0), json, cap);
}

extern "C" int32_t caps_label_kinds(char* json, int32_t cap) {
  caps::Json j = caps::Json::object();
  for (const auto& [name, list] : {std::make_pair("atom", &caps::atom_label_kinds()), std::make_pair("bond", &caps::bond_label_kinds())}) {
    caps::Json a = caps::Json::array();
    for (const auto& k : *list) {
      caps::Json o = caps::Json::object();
      o["id"] = k.id;
      o["title"] = k.title;
      o["group"] = k.group;
      a.push_back(o);
    }
    j[name] = a;
  }
  return report_out(j.dump(0), json, cap);
}

extern "C" int32_t caps_project_indices(caps_doc* d, const caps_camera* cam, const caps_render_opts* opt, const int32_t* atoms, int32_t n, float* xy) {
  return guard([&] {
    std::vector<int> idx(atoms, atoms + std::max(0, n));
    const auto p = caps::Renderer::project_some(d->frame, cam_of(cam), opts_of(d, opt), idx);
    std::memcpy(xy, p.data(), p.size() * sizeof(float));
    return int32_t(idx.size());
  });
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
  // the force field an AMBER topology carried holds while the atoms (elements) and bonds are unchanged
  bool same = s.forcefield && s.atoms.size() == d->traj.topology.atoms.size() && s.bonds.size() == d->traj.topology.bonds.size();
  for (size_t i = 0; same && i < s.atoms.size(); ++i) same = s.atoms[i].element == d->traj.topology.atoms[i].element;
  for (size_t i = 0; same && i < s.bonds.size(); ++i) same = s.bonds[i].i == d->traj.topology.bonds[i].i && s.bonds[i].j == d->traj.topology.bonds[i].j;
  const bool file_ff = same && d->field && !d->field->model.empty() && caps::Json::parse(d->field->model).text("model", "") == "file";
  d->traj.topology = s;
  if (!same) d->traj.topology.forcefield.reset();
  std::vector<caps::Vec3> p;
  for (const auto& a : s.atoms) p.push_back(a.pos);
  d->traj.positions[d->current] = std::move(p);
  if (d->current < d->traj.cells.size()) d->traj.cells[d->current] = s.cell;
  d->field.reset();   // the typing no longer matches
  d->scan_frames.clear();
  if (d->selection.size() != s.atoms.size()) d->selection.assign(s.atoms.size(), 0);
  refresh(d);
  if (file_ff) install_file_field(d);   // moved atoms, the same molecule: its topology's force field still applies
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
    } else if (op == "fix_h") {   // H autopilot: surplus hydrogens of the atoms removed, missing ones added
      const auto at = atoms_of(d, j);
      std::vector<char> m;
      if (!at.empty()) { m.assign(s.atoms.size(), 0); for (uint32_t a : at) m[a] = 1; }
      if (needs_geometry_orders(s)) geometry_orders(s);
      const auto [added_h, removed_h] = caps::fix_hydrogens(s, m);
      if (added_h == 0 && removed_h == 0) throw std::invalid_argument("the hydrogens are right");
      r["h_added"] = double(added_h), r["h_removed"] = double(removed_h);
      what = "H autopilot · +" + std::to_string(added_h) + " H, −" + std::to_string(removed_h) + " H";
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
    } else if (op == "functionalize") {   // groups on sidewalls, ends, edges or chosen atoms of a filler (functionalize.hpp)
      caps::FunctionalizeOptions f;
      f.group = j.text("group", "hydroxyl");
      f.pattern = j.text("pattern", "random");
      f.elements = j.text("elements", "");
      f.fraction = j.num("fraction", 0.05);
      f.count = int(j.num("count", 0));
      f.min_spacing = j.num("min_spacing", 3.0);
      f.from = j.num("from", 0.0), f.to = j.num("to", 1.0);
      f.pitch = j.num("pitch", 20.0), f.phase = j.num("phase", 0.0);
      f.side = j.text("side", "outer");
      f.seed = uint64_t(j.num("seed", 1));
      if (f.pattern == "atoms") for (uint32_t a : atoms_of(d, j)) f.atoms.push_back(a);
      const auto fr = caps::functionalize(s, f);
      what = fr.notes.empty() ? "Functionalised" : fr.notes.front();
    } else if (op == "phosphate_ends") {   // P–H (a strand's 3′ cap) → P–OH
      const int k = caps::hydroxylate_phosphorus(s);
      if (k == 0) throw std::invalid_argument("no hydrogen on phosphorus");
      what = "Phosphate ends: " + std::to_string(k) + " P–H to P–OH";
    } else if (op == "type_element") {   // {type, element, mass?}: every atom of a type made this element; the type is kept
      const int t = int(j.num("type", 0));
      const int z = element_of(j.text("element"));
      if (z <= 0) throw std::invalid_argument("unknown element");
      size_t n = 0;
      for (auto& a : s.atoms) if (a.type == t) a.element = z, ++n;
      if (n == 0) throw std::invalid_argument("no atom has type " + std::to_string(t));
      if (j.has("mass")) {
        const double m = j.num("mass", 0);
        if (m <= 0) throw std::invalid_argument("a mass must be positive");
        for (auto& ti : s.types) if (ti.type == t) ti.mass = m;
      }
      what = "Type " + std::to_string(t) + ": " + std::to_string(n) + " atom(s) made " + caps::element(z).symbol + (j.has("mass") ? " (mass set)" : "");
    } else if (op == "redefine_lattice") {   // {matrix: [9 numbers, row-major; column j = the new vector j in the old ones]}
      if (!j.has("matrix") || !j["matrix"].is_array() || j["matrix"].size() != 9) throw std::invalid_argument("redefine_lattice needs matrix: [9 numbers]");
      caps::Mat3 m{};
      for (int q = 0; q < 9; ++q) m[size_t(q / 3)][size_t(q % 3)] = j["matrix"][size_t(q)].number();
      const size_t n0 = s.atoms.size();
      s = caps::transform_cell(s, m, j.num("tolerance", 0.05));
      what = "Redefine lattice: " + std::to_string(n0) + " → " + std::to_string(s.atoms.size()) + " atoms";
    } else if (op == "water_model") {   // {model: "tip4p2005" | "TIP4P/2005" | …}: geometry, charges, M sites of every water
      const caps::WaterModel& wm = caps::water_model(j.text("model", ""));
      std::vector<std::string> wn;
      const size_t nw = caps::apply_water_model(s, wm, &wn);
      if (nw == 0) throw std::invalid_argument("no water molecules (an O bonded to two H) in this structure");
      what = std::to_string(nw) + " waters as " + wm.name + (wm.sites == 4 ? " (with their M sites)" : wm.sites == 5 ? " (with their lone pairs)" : "");
    } else if (op == "niggli") {
      caps::NiggliResult r;
      s = caps::niggli_cell(s, &r);
      char b[160];
      std::snprintf(b, sizeof b, "Niggli cell: a %.4f b %.4f c %.4f Å, α %.2f β %.2f γ %.2f°", r.a, r.b, r.c, r.alpha, r.beta, r.gamma);
      what = b;
    } else if (op == "find_primitive") {   // {tolerance}
      int k = 1;
      const size_t n0 = s.atoms.size();
      s = caps::find_primitive_cell(s, j.num("tolerance", 0.1), &k);
      what = k == 1 ? "Already primitive" : "Primitive cell: " + std::to_string(n0) + " → " + std::to_string(s.atoms.size()) + " atoms";
      if (k == 1) throw std::invalid_argument("the cell is already primitive (no translation maps the structure onto itself)");
    } else if (op == "conventional") {   // {tolerance}
      caps::ConventionalResult r;
      s = caps::conventional_cell(s, j.num("tolerance", 0.1), &r);
      what = "Conventional cell: " + r.hm + " (No. " + std::to_string(r.number) + "), " + std::to_string(s.atoms.size()) + " atoms";
    } else if (op == "supercell") {   // {n: [na, nb, nc]}
      if (!j.has("n") || !j["n"].is_array() || j["n"].size() != 3) throw std::invalid_argument("supercell needs n: [na, nb, nc]");
      const int na = int(j["n"][0].number()), nb = int(j["n"][1].number()), nc = int(j["n"][2].number());
      if (na < 1 || nb < 1 || nc < 1) throw std::invalid_argument("supercell repeats must be 1 or more");
      if (double(na) * nb * nc * double(s.atoms.size()) > 2e6) throw std::invalid_argument("the supercell would hold more than two million atoms");
      s = caps::supercell(s, na, nb, nc);
      what = "Supercell " + std::to_string(na) + " × " + std::to_string(nb) + " × " + std::to_string(nc);
    } else if (op == "vacuum_slab") {   // {vacuum Å, centre: bool}
      caps::SlabResult r;
      s = caps::vacuum_slab(s, j.num("vacuum", 15), !(j.has("centre") && j["centre"].kind() == caps::Json::Bool && !j["centre"].boolean()), &r);
      char b[128];
      std::snprintf(b, sizeof b, "Vacuum slab: %.2f Å slab, %.2f Å vacuum", r.thickness, r.vacuum);
      what = b;
    } else if (op == "defects") {   // {from, to ("" vacancies), fraction | count, min_spacing, seed, atoms (the region)}
      caps::DefectOptions o;
      o.from = element_of(j.text("from", "C"));
      const std::string to = j.text("to", "");
      o.to = to.empty() ? 0 : element_of(to);
      o.fraction = j.num("fraction", 0.05);
      o.count = int(j.num("count", 0));
      o.min_spacing = j.num("min_spacing", 0);
      o.seed = uint64_t(j.num("seed", 1));
      if (j.has("atoms")) {
        o.region.assign(s.atoms.size(), 0);
        for (uint32_t a : atoms_of(d, j)) o.region[a] = 1;
      }
      std::vector<std::string> nt;
      caps::point_defects(s, o, &nt);
      what = nt.front();
      for (size_t k = 1; k < nt.size(); ++k) what += " · " + nt[k];
    } else if (op == "backbone_torsions") {   // {pattern: [deg …], atoms: [...] (their chains; none: every chain)}
      std::vector<double> pat;
      if (j.has("pattern") && j["pattern"].is_array())
        for (const auto& x : j["pattern"].items()) pat.push_back(x.number());
      std::vector<uint32_t> at;
      if (j.has("atoms")) at = atoms_of(d, j);
      std::vector<std::string> nt;
      const int k = caps::set_backbone_torsions(s, pat, at, &nt);
      if (k == 0) throw std::invalid_argument("no backbone dihedral to set (a chain needs four backbone atoms)");
      what = "Backbone torsions: " + nt.front();
    } else if (op == "cluster") {   // {radius (≤ 0: all), centre: [x, y, z] | "selection" | "cell", any_atom: bool}
      caps::Vec3 c = s.cell.valid() ? s.cell.origin + (s.cell.a + s.cell.b + s.cell.c) * 0.5 : caps::Vec3{0, 0, 0};
      if (j.has("centre") && j["centre"].is_array() && j["centre"].size() == 3) c = {j["centre"][0].number(), j["centre"][1].number(), j["centre"][2].number()};
      else if (j.text("centre", "") == "selection") {
        const auto at = atoms_of(d, j);
        if (at.empty()) throw std::invalid_argument("select the atoms to centre the cluster on");
        caps::Vec3 m{0, 0, 0};
        for (uint32_t a : at) m = m + s.atoms[a].pos;
        c = m * (1.0 / double(at.size()));
      }
      caps::ClusterResult cr;
      s = caps::cut_cluster(s, c, j.num("radius", 0), j.has("any_atom") && j["any_atom"].boolean(), &cr);
      char b[200];
      if (j.num("radius", 0) > 0)
        std::snprintf(b, sizeof b, "Cluster: %zu whole molecules within %.1f Å (%zu atoms, %zu molecules left out), no cell", cr.molecules, j.num("radius", 0), cr.atoms, cr.dropped);
      else
        std::snprintf(b, sizeof b, "Periodicity removed: %zu whole molecules, %zu atoms, no cell", cr.molecules, cr.atoms);
      what = b;
    } else if (op == "nanowire") {   // {uvw: [u, v, w], radius, repeats, shape, vacuum}
      caps::WireOptions o;
      if (j.has("uvw") && j["uvw"].is_array() && j["uvw"].size() == 3)
        o.uvw = {int(j["uvw"][0].number()), int(j["uvw"][1].number()), int(j["uvw"][2].number())};
      o.radius = j.num("radius", o.radius);
      o.repeats = int(j.num("repeats", o.repeats));
      o.shape = j.text("shape", o.shape);
      o.vacuum = j.num("vacuum", o.vacuum);
      caps::WireResult r;
      s = caps::nanowire(s, o, &r);
      what = "Nanowire: " + std::to_string(r.atoms) + " atoms";
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
    } else if (op == "set_coordination") {   // {atom, geometry}: its neighbours at the ideal directions of the geometry
      const uint32_t c = uint32_t(j.num("atom", -1));
      const std::string g = j.text("geometry", "");
      const double moved = caps::set_coordination(s, c, g);
      char b[160];
      std::snprintf(b, sizeof b, "Atom %u made %s (ligands moved up to %.1f°)", c + 1, g.c_str(), moved);
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
    } else if (op == "bond_rules") {   // {rules: [{a, b, max, never}]}: the bonds from the rules (orders of kept bonds kept)
      const auto nb = caps::bonds_by_rules(s, rules_of(j));
      std::map<uint64_t, int> order;
      for (const auto& b : s.bonds) order[b.i < b.j ? (uint64_t(b.i) << 32) | b.j : (uint64_t(b.j) << 32) | b.i] = b.order;
      std::vector<caps::Bond> next;
      for (auto b : nb) {
        const auto it = order.find((uint64_t(b.i) << 32) | b.j);
        if (it != order.end()) b.order = it->second;
        next.push_back(b);
      }
      const long delta = long(next.size()) - long(s.bonds.size());
      s.bonds = std::move(next);
      s.bonds_from_file = false;
      d->overlay.clear();
      what = "Bond rules · " + std::to_string(s.bonds.size()) + " bonds (" + (delta >= 0 ? "+" : "") + std::to_string(delta) + ")";
    } else if (op == "stamp") {   // {piece: caps-piece JSON text, at: [x, y, z], axis, degrees, clear: Å} (design/boards/Stamp)
      caps::Vec3 at{0, 0, 0}, ax{0, 0, 1};
      if (j.has("at") && j["at"].size() == 3) at = {j["at"][0].number(), j["at"][1].number(), j["at"][2].number()};
      if (j.has("axis") && j["axis"].size() == 3) ax = {j["axis"][0].number(), j["axis"][1].number(), j["axis"][2].number()};
      const auto r = caps::stamp_piece(s, j.text("piece"), at, ax, j.num("degrees", 0), j.num("clear", 1.5));
      for (size_t i : r.added) added.push_back(double(i));
      char b[160];
      const double want = j.num("clear", 1.5);
      std::snprintf(b, sizeof b, " · %zu atoms%s%s", r.added.size(), r.shift > 0 ? (", moved " + std::to_string(r.shift).substr(0, 4) + " Å").c_str() : "",
                    r.closest < want && r.closest > 0 ? (" · closest contact " + std::to_string(r.closest).substr(0, 4) + " Å: no room for " + std::to_string(want).substr(0, 3) + " Å here, Minimise clears it").c_str() : "");
      what = "Stamp " + j.text("name", "piece") + b;
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
    caps::Image img;
    if (o.text("engine", "raster") == "raytrace") {
      // CAPS's ray tracer on the same scene and camera (raytrace.hpp): occlusion, soft shadows, depth of field
      caps::RenderOptions rs = ro;
      rs.supersample = 1;
      const caps::System& sys = shown(d);
      const caps::Scene sc = r.scene(sys, rs);
      const caps::ViewFit fit = caps::view_fit(sys, cam_of(cam), rs);
      caps::TraceOptions to;
      to.samples = int(std::clamp(o.num("samples", 64), 1.0, 4096.0));
      to.ambient_occlusion = flag("occlusion", true);
      to.shadows = flag("shadows", true);
      to.ao_strength = std::clamp(o.num("occlusion_strength", to.ao_strength), 0.0, 1.0);
      to.aperture = std::max(0.0, o.num("aperture", 0.0));
      // or relative to the view: a fraction of its half-width in Å (the same blur on a molecule or a large cell)
      if (o.has("aperture_fraction")) to.aperture = std::max(0.0, o.num("aperture_fraction", 0.0)) * (rs.width * 0.5 / std::max(1e-9, fit.scale));
      to.focus = o.num("focus", 0.0);
      to.outlines = flag("outlines", false);
      img = caps::raytrace(sc, fit, rs.width, rs.height, to);
      if (!p.sixteen) img.rgba16.clear();
    } else {
      img = r.render(shown(d), cam_of(cam), ro);
    }
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

// The view's scene for other renderers and 3D tools (scene_export.hpp): format "pov" | "glb" | "obj", with the camera
// of `cam` for opt's size (POV-Ray). report: {spheres, cylinders, triangles}.
extern "C" int32_t caps_export_scene(caps_doc* d, const caps_camera* cam, const caps_render_opts* opt, const char* path, const char* format,
                                     char* report, int32_t cap) {
  return guard([&] {
    caps::RenderOptions ro = opts_of(d, opt);
    ro.supersample = 1;
    const caps::System& sys = shown(d);
    caps::Renderer r;
    const caps::Scene sc = r.scene(sys, ro);
    const std::string f = format ? format : "glb";
    caps::SceneExportReport rep;
    if (f == "pov") rep = caps::write_povray(sc, caps::view_fit(sys, cam_of(cam), ro), ro.width, ro.height, path ? path : "");
    else if (f == "obj") rep = caps::write_obj(sc, path ? path : "");
    else if (f == "glb" || f == "gltf") rep = caps::write_gltf(sc, path ? path : "");
    else throw std::invalid_argument("scene format: pov, glb or obj");
    caps::Json j = caps::Json::object();
    j["spheres"] = double(rep.spheres), j["cylinders"] = double(rep.cylinders), j["triangles"] = double(rep.triangles);
    report_out(j.dump(), report, cap);
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
    if (const std::string rk = o.text("radii", "bondi"); rk != "bondi")   // uff | forcefield (the Field assignment)
      vo.radii = caps::free_volume_radii(d->frame, rk, d->field && d->field->ff ? d->field->ff.get() : nullptr);
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
      // united-atom fluid (TraPPE-UA, GROMOS): hydrogens on carbon folded into their carbons before packing — CH4 one site
      if (flag("united_atom", false)) fluid = caps::united_atom(fluid);
      fluid.title = j.text("fluid_name", smiles);
      o.fluid = &fluid;
    }
    caps::PoreReport r;
    const caps::System s = caps::build_pore(o, &r);
    if (o.fluid && flag("united_atom", false))
      r.notes.push_back("the fluid in united atoms (hydrogens folded into their carbons): assign TraPPE-UA to it in Field · by group, the walls their own");
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
  o.sigma = j.num("sigma", 0);             // real units: σ (Å), T (K) for ε = k_B T, bead mass (g/mol) — all three, else reduced
  o.temperature = j.num("temperature", 0);
  o.bead_mass = j.num("bead_mass", 0);
  return o;
}
}  // namespace

extern "C" caps_doc* caps_kg_build(const char* options_json, char* report, int32_t cap) {
  try {
    const caps::KgOptions o = kg_options(options_json);
    caps::KgReport r;
    caps::System s = caps::kremer_grest(o, &r);
    const bool mapped = caps::kg_units(o).eps > 0;
    if (mapped) {   // the melt in Å: σ as mapped (its density then in real g/cm³)
      for (auto& a : s.atoms) a.pos = a.pos * o.sigma;
      s.cell.a = s.cell.a * o.sigma, s.cell.b = s.cell.b * o.sigma, s.cell.c = s.cell.c * o.sigma, s.cell.origin = s.cell.origin * o.sigma;
      for (auto& t : s.types) t.mass = o.bead_mass;
    }
    caps::Json j = caps::Json::object();
    j["box"] = r.box;
    j["closest"] = r.closest;
    j["r2_per_bond"] = r.mean_r2;
    report_out(j.dump(0), report, cap);
    caps_doc* d = doc_of(s);
    prov_step(d, "cg.kremer_grest", std::to_string(o.chains) + " × " + std::to_string(o.beads) + " bead-spring chains as random walks",
              {{"chains", std::to_string(o.chains)}, {"beads", std::to_string(o.beads)}, {"density", g6(o.density) + " σ⁻³"}, {"k_theta", g6(o.k_theta) + " ε"},
               {"box", g6(r.box) + " σ"}}, seeded(o.seed), {"kremer1990"});
    // its own force field: the melt is ready to run, analyse and export without a library assignment
    caps::Json m = caps::Json::object();
    m["model"] = std::string("kremer-grest");
    m["k_theta"] = o.k_theta, m["density"] = o.density, m["chains"] = double(o.chains), m["beads"] = double(o.beads), m["seed"] = double(o.seed);
    if (mapped) m["sigma"] = o.sigma, m["temperature"] = o.temperature, m["bead_mass"] = o.bead_mass;
    d->field = std::make_unique<FieldState>();
    for (const auto& a : d->traj.topology.atoms) d->field->file_types.push_back({a.type, a.name}), d->field->file_charges.push_back(a.charge);
    d->field->file_type_table = d->traj.topology.types;
    d->field->model = m.dump(0);
    field_run(d);
    return d;
  } catch (const std::exception& e) {
    g_error = e.what();
    return nullptr;
  }
}

extern "C" caps_doc* caps_martini_melt(const char* options_json, caps_progress_fn progress, void* user, char* report, int32_t cap) {
  try {
    const caps::Json j = caps::Json::parse(options_json && *options_json ? options_json : "{}");
    const std::string ff_path = j.text("forcefield", "");
    const std::string repeat = j.text("repeat", "");
    const int n = int(j.num("repeats", 20)), chains = int(j.num("chains", 20));
    const double density = j.num("density", 1.0), tol = j.num("tolerance", 3.0);
    const uint64_t seed = uint64_t(j.num("seed", 1));
    if (ff_path.empty()) throw std::runtime_error("choose the MARTINI force field");
    if (repeat.empty()) throw std::runtime_error("give the repeat unit as bead SMILES, e.g. [SN0] (PEO) or [C1] (four CH2 per bead)");
    if (n < 1 || chains < 1 || density <= 0) throw std::runtime_error("repeats, chains and the density must be positive");
    const caps::FFDef def = caps::load_forcefield(ff_path);
    {   // every bead of the repeat must be a type of this force field ([SN0], [Qd+1]: the name before a charge)
      std::string missing;
      for (size_t a = repeat.find('['); a != std::string::npos; a = repeat.find('[', a + 1)) {
        const size_t b = repeat.find(']', a);
        if (b == std::string::npos) throw std::runtime_error("unclosed [ in the repeat unit " + repeat);
        std::string name = repeat.substr(a + 1, b - a - 1);
        const size_t q = name.find_first_of("+-", 1);
        if (q != std::string::npos) name.resize(q);
        bool known = def.type(name) != nullptr;   // or a moltemplate-style name built on it (SN0_bSN0_aSN0_…)
        for (size_t t = 0; !known && t < def.types.size(); ++t) known = def.types[t].name.rfind(name + "_", 0) == 0;
        if (!known && missing.find(" " + name + ",") == std::string::npos) missing += " " + name + ",";
      }
      if (!missing.empty()) {
        missing.pop_back();
        throw std::runtime_error(def.name + " has no bead type" + missing + " (MARTINI 2 names such as C1, SN0 are not Martini 3's; choose the force field the mapping is written for)");
      }
    }
    std::string chain_text;
    for (int k = 0; k < n; ++k) chain_text += repeat;
    caps::System chain = caps::build_bead_molecule(chain_text, def, seed);
    chain.title = "chain";
    for (auto& a : chain.atoms) a.mol = 1;
    chain.has_mol = true;
    // pack loosely (a quarter of the target density), then compress with the MARTINI force field itself
    double mass = 0;
    for (const auto& a : chain.atoms) {
      const caps::FFType* t = def.type(a.name);
      mass += t && t->mass > 0 ? t->mass : 72.0;
    }
    const double loose = 0.25 * density;
    const double L = std::cbrt(chains * mass / 6.02214076e23 / loose * 1e24);
    caps::PackItem it;
    it.name = "chain";
    it.molecule = chain;
    it.count = chains;
    caps::PackOptions po;
    po.cell.a = {L, 0, 0}, po.cell.b = {0, L, 0}, po.cell.c = {0, 0, L};
    po.periodic = true;
    po.tolerance = tol;
    po.seed = seed;
    caps::PackReport pr;
    caps::System packed = caps::pack({it}, po, &pr);
    packed.title = "MARTINI melt";
    // the force field on the loose cell, for the compression
    std::unique_ptr<caps_doc> tmp(doc_of(packed));
    if (caps_field_assign(tmp.get(), ff_path.c_str(), nullptr, 2) < 0)   // the charges of the bead SMILES throw std::runtime_error(g_error);
    if (!tmp->field->complete) throw std::runtime_error(def.name + " does not describe every bead of " + repeat + ": see Field for the untyped beads and missing terms");
    caps::RelaxOptions ro;
    ro.field = tmp->field->ff;
    ro.energy = elec();
    ro.target_density = density;
    ro.compress_step = 0.06;
    ro.ftol = 2.0;
    ro.max_iterations = 2000;
    if (progress) ro.progress = [&](const caps::RelaxProgress& p) { return progress(p.stage_index, p.stages, 0, user) == 0; };
    caps::RelaxReport rr;
    caps::System melt = packed;
    caps::relax(melt, ro, &rr);
    melt.velocities.clear();
    caps_doc* d = doc_of(melt);
    if (caps_field_assign(d, ff_path.c_str(), nullptr, 2) < 0) { const std::string e = g_error; delete d; throw std::runtime_error(e); }
    char b[256];
    std::snprintf(b, sizeof b, "%d chains of %d × %s packed at %.3f g/cm³, compressed with %s to %.3f g/cm³", chains, n, repeat.c_str(), packed.density(), def.name.c_str(),
                  melt.density());
    prov_step(d, "cg.martini_melt", b, {{"forcefield", def.name}, {"repeat", repeat}, {"repeats", std::to_string(n)}, {"chains", std::to_string(chains)},
                                        {"density", g6(density) + " g/cm³"}}, seeded(seed), {});
    caps::Json rj = caps::Json::object();
    rj["beads"] = double(melt.atoms.size());
    rj["beads_per_chain"] = double(chain.atoms.size());
    rj["chain_mass"] = mass;
    rj["density"] = melt.density();
    rj["loose_density"] = packed.density();
    rj["forcefield"] = def.name;
    rj["complete"] = d->field && d->field->complete;
    report_out(rj.dump(0), report, cap);
    return d;
  } catch (const std::exception& e) {
    g_error = e.what();
    return nullptr;
  }
}

extern "C" int32_t caps_kg_lammps(caps_doc* d, const char* options_json, const char* stem, double pushoff_steps, double run_steps) {
  return guard([&] {
    // the structure in σ units (a mapped melt is held in Å)
    caps::System s = d->frame;
    if (d->field && !d->field->model.empty()) {
      const double sg = caps::Json::parse(d->field->model).num("sigma", 0);
      if (sg > 0) {
        for (auto& a : s.atoms) a.pos = a.pos * (1.0 / sg);
        s.cell.a = s.cell.a * (1.0 / sg), s.cell.b = s.cell.b * (1.0 / sg), s.cell.c = s.cell.c * (1.0 / sg), s.cell.origin = s.cell.origin * (1.0 / sg);
      }
    }
    caps::write_kg_lammps(s, kg_options(options_json), stem, pushoff_steps > 0 ? pushoff_steps : 20000, run_steps > 0 ? run_steps : 100000);
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

namespace {
// A state of the document as positions, with its name: {kind: current | start | frame | snapshot, index}
std::pair<std::vector<caps::Vec3>, std::string> state_positions(caps_doc* d, const caps::Json& st, const std::string& def) {
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
}
}  // namespace

// One state of the document as a new document (the split view shows two states side by side).
extern "C" caps_doc* caps_state_document(caps_doc* d, const char* json) {
  try {
    if (!d) throw std::runtime_error("no document");
    const auto [pos, name] = state_positions(d, caps::Json::parse(json && *json ? json : "{}"), "current");
    caps::System s = d->frame;
    if (pos.size() != s.atoms.size()) throw std::runtime_error("that state has other atoms than the structure shown");
    for (size_t i = 0; i < pos.size(); ++i) s.atoms[i].pos = pos[i];
    return doc_of_system(s, d, "state.copy", "the structure at " + name, {{"state", name}});
  } catch (const std::exception& e) {
    g_error = e.what();
    return nullptr;
  }
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
    auto state = [&](const std::string& key, const std::string& def) {
      return state_positions(d, o.has(key) ? o[key] : caps::Json::object(), def);
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
    // over the trajectory: {"frames": {first, last (−1: the last), every}} → the total per frame, its mean ± sd
    if (j.has("frames") && d->traj.frames() > 1) {
      const auto& fr = j["frames"];
      const long nf = long(d->traj.frames());
      long first = long(fr.num("first", 0)), last = long(fr.num("last", -1)), every = std::max(1L, long(fr.num("every", 1)));
      if (last < 0 || last >= nf) last = nf - 1;
      first = std::clamp(first, 0L, last);
      caps::Json series = caps::Json::array();
      std::vector<double> tot;
      for (long f = first; f <= last; f += every) {
        const caps::System sf = d->traj.frame(size_t(f));
        const double t = caps::sasa(sf, probe, points).total;
        tot.push_back(t);
        caps::Json row = caps::Json::object();
        row["frame"] = double(f), row["total"] = t;
        if (d->traj.timesteps.size() > size_t(f)) row["timestep"] = double(d->traj.timesteps[size_t(f)]);
        series.push_back(std::move(row));
      }
      double m = 0, v = 0;
      for (double t : tot) m += t;
      m /= double(tot.size());
      for (double t : tot) v += (t - m) * (t - m);
      r["series"] = std::move(series);
      r["mean"] = m;
      r["sd"] = tot.size() > 1 ? std::sqrt(v / double(tot.size() - 1)) : 0.0;
      r["frames_used"] = double(tot.size());
    }
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
    if (j.has("adjust") && j["adjust"].is_object()) {   // average equivalent atoms, scale, neutralise
      const auto& a = j["adjust"];
      caps::ChargeAdjust ad;
      ad.average = a.has("average") && a["average"].boolean();
      ad.scale = a.num("scale", 1.0);
      ad.neutralise = a.text("neutralise", "none");
      if (a.has("target") && a["target"].kind() == caps::Json::Number) ad.target = a.num("target", 0);
      std::vector<std::string> adj;
      caps::adjust_charges(s, rep.q, ad, rep.formal, &adj);
      auto notes = rep.notes;
      rep = caps::describe_charges(s, rep.q, rep.method);
      rep.notes = notes;
      for (auto& x : adj) rep.notes.push_back(x);
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
    // a histogram: [[length, number weight], …] (a measured distribution binned by chain length)
    std::vector<std::pair<int, double>> hist;
    if (dist == "histogram") {
      if (!j.has("histogram") || !j["histogram"].is_array()) throw std::invalid_argument("histogram: [[length, weight], …]");
      for (const auto& b : j["histogram"].items())
        if (b.is_array() && b.size() == 2) hist.push_back({int(b[0].number()), b[1].number()});
    }
    double hist_nn = 0, hist_pdi = 1;
    if (!hist.empty()) {
      double w0 = 0, w1 = 0, w2 = 0;
      for (const auto& [n, w] : hist) w0 += w, w1 += w * n, w2 += w * double(n) * n;
      if (w0 > 0 && w1 > 0) hist_nn = w1 / w0, hist_pdi = (w2 / w1) / hist_nn;
    }
    for (int t = 0; t < best_of; ++t) {
      auto c = !hist.empty() ? caps::draw_chain_lengths(hist, count, seed + uint64_t(t) * 7919) : caps::draw_chain_lengths(dist, nn, pdi, count, seed + uint64_t(t) * 7919);
      double a, b;
      stats(c, a, b);
      const double target = !hist.empty() ? hist_pdi : dist == "flory" ? 2 - 1 / nn : dist == "poisson" ? 1 + (nn - 1) / (nn * nn) : dist == "monodisperse" ? 1.0 : pdi;
      const double nref = !hist.empty() ? hist_nn : nn;
      const double miss = std::fabs(b / a - target) + 0.1 * std::fabs(a - nref) / nref;
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
    const double tpdi = !hist.empty() ? hist_pdi : dist == "schulz-zimm" || dist == "log-normal" ? pdi : dist == "flory" ? 2 - 1 / nn : dist == "poisson" ? 1 + (nn - 1) / (nn * nn) : 1.0;
    const double tnn = !hist.empty() ? hist_nn : nn;
    caps::Json t = caps::Json::object();
    t["nn"] = tnn, t["mn"] = tnn * m0, t["mw"] = tnn * m0 * tpdi, t["pdi"] = tpdi;
    r["target"] = std::move(t);
    r["k"] = dist == "schulz-zimm" ? (pdi > 1.0001 ? 1 / (pdi - 1) : 1e4) : 0.0;
    r["kept_draw"] = double(kept);
    r["draws"] = double(best_of);
    std::vector<double> xs, nf, wf;
    const double hi = std::max(double(sorted.back()) * 1.15, nn * (1 + 4 * std::sqrt(std::max(0.0, tpdi - 1)) + 0.3));
    if (!hist.empty()) {   // the histogram itself, as number and weight fractions at its lengths
      double w0 = 0;
      for (const auto& [n, w] : hist) w0 += w;
      auto hs = hist;
      std::sort(hs.begin(), hs.end());
      for (const auto& [n, w] : hs) xs.push_back(n), nf.push_back(w / w0), wf.push_back(n * (w / w0) / tnn);
    } else
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
    d->fixed_atoms = src->fixed_atoms;
    d->fixed_axes = src->fixed_axes;
    d->rigid_mols = src->rigid_mols;
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

// Structure-based coarse-graining (cg_map.hpp): the document's all-atom structure — every frame of it — mapped to beads,
// with Boltzmann-inverted bonds and angles and a repulsive WCA from the intermolecular bead g(r). A new document carrying
// its bead model. options: {scheme: unit | backbone_side | backbone_n, per_bead, temperature}. Report JSON: notes, bond
// and angle types, sigma, cut, epsilon, g(r).
namespace {
caps::CgMapOptions cg_map_options(const caps::Json& j) {
  caps::CgMapOptions o;
  o.scheme = j.text("scheme", "unit");
  o.per_bead = int(j.num("per_bead", 3));
  o.temperature = j.num("temperature", 300);
  return o;
}
std::string cg_map_report(const caps::CgMapResult& r) {
  caps::Json j = caps::Json::object();
  caps::Json notes = caps::Json::array(), bonds = caps::Json::array(), angles = caps::Json::array(), gr = caps::Json::array();
  for (const auto& n : r.notes) notes.push_back(n);
  for (const auto& b : r.bonds) {
    caps::Json o = caps::Json::object();
    o["types"] = b.a + "–" + b.b, o["r0"] = b.r0, o["k"] = b.k, o["sd"] = b.sd, o["count"] = double(b.count);
    bonds.push_back(o);
  }
  for (const auto& a : r.angles) {
    caps::Json o = caps::Json::object();
    o["types"] = a.a + "–" + a.b + "–" + a.c, o["theta0"] = a.theta0 * 180 / 3.14159265358979323846, o["k"] = a.k, o["sd"] = a.sd * 180 / 3.14159265358979323846,
    o["count"] = double(a.count);
    angles.push_back(o);
  }
  for (const auto& [x, y] : r.gr) { caps::Json p = caps::Json::array(); p.push_back(x); p.push_back(y); gr.push_back(p); }
  j["ok"] = true, j["notes"] = notes, j["bonds"] = bonds, j["angles"] = angles, j["gr"] = gr;
  j["sigma"] = r.sigma, j["cut"] = r.cut, j["epsilon"] = r.epsilon, j["beads"] = double(r.beads.atoms.size());
  return j.dump(0);
}
caps_doc* cg_doc(caps::CgMapResult& r, const caps_doc* from, const std::string& summary, caps::KeyValues params) {
  caps::System beads = r.beads;
  caps::Trajectory t;
  t.topology = beads;
  t.positions = r.frames;
  t.cells = r.cells;
  for (size_t k = 0; k < r.frames.size(); ++k) t.timesteps.push_back(int64_t(k));
  auto* d = new caps_doc;
  d->traj = std::move(t);
  d->current = d->traj.frames() - 1;
  refresh(d);
  if (from) d->prov = from->prov;
  prov_step(d, "cg.map", summary, std::move(params), "", {"tschop1998", "reith2003"});
  install_file_field(d);   // the bead model goes with the beads
  return d;
}
}  // namespace

// v62: the bead model's non-bonded part refined by iterative Boltzmann inversion against the mapped g(r) (options "ibi":
// {iterations, run_ps, temperature, pressure_correction}); the table replaces the WCA, the bonded terms stay
void apply_ibi(caps::CgMapResult& r, const caps::Json& j, double temperature) {
  if (!j.has("ibi") || !j["ibi"].is_object()) return;
  const auto& b = j["ibi"];
  if (r.gr.size() < 10) throw std::invalid_argument("IBI needs the mapped g(r) (a periodic cell of several chains)");
  std::vector<double> rr, gg;
  for (const auto& [x, y] : r.gr) rr.push_back(x), gg.push_back(y);
  caps::IbiOptions io;
  io.temperature = b.num("temperature", temperature);
  io.iterations = std::clamp(int(b.num("iterations", 6)), 1, 50);
  io.run_ps = b.num("run_ps", 20);
  io.first_equilibrate_ps = b.num("equilibrate_ps", 5);
  io.pressure_correction = b.has("pressure_correction") && b["pressure_correction"].boolean();
  io.seed = uint64_t(b.num("seed", 1));
  const auto res = caps::run_ibi(r.beads, *r.ff, rr, gg, io);
  r.ff = res.ff;
  r.beads.forcefield = res.ff;   // the document's model is the structure's own
  for (const auto& n : res.notes) r.notes.push_back(n);
  std::string hist;
  for (const auto& h : res.history) hist += (hist.empty() ? "" : ", ") + std::string(std::to_string(h.residual).substr(0, 6));
  r.notes.push_back("IBI residual per iteration: " + hist + " (one table for every bead pair; LAMMPS: pair_style table, caps_pairs.table)");
}

extern "C" caps_doc* caps_cg_map(caps_doc* d, const char* options_json, char* report, int32_t cap) {
  try {
    const caps::Json j = caps::Json::parse(options_json && *options_json ? options_json : "{}");
    const caps::CgMapOptions o = cg_map_options(j);
    caps::System aa = d->traj.topology;
    caps::CgMapResult r = caps::cg_map(aa, o, d->traj.positions, d->traj.cells);
    apply_ibi(r, j, o.temperature);
    report_out(cg_map_report(r), report, cap);
    return cg_doc(r, d, "all-atom structure mapped to beads (" + o.scheme + ")",
                  {{"scheme", o.scheme}, {"per_bead", std::to_string(o.per_bead)}, {"temperature", g6(o.temperature)}, {"frames", std::to_string(r.frames.size())}});
  } catch (const std::exception& e) {
    g_error = e.what();
    return nullptr;
  }
}

// A polymer (spec JSON as Grow takes it) coarse-grained from an all-atom reference melt: `chains` × `dp` grown, compressed
// to `density` (g/cm³) by the built-in force field (push-off, then minimisation), then mapped (options as caps_cg_map).
extern "C" caps_doc* caps_cg_from_polymer(const char* spec_json, const char* options_json, caps_progress_fn progress, void* user, char* report, int32_t cap) {
  try {
    const caps::Json j = caps::Json::parse(options_json && *options_json ? options_json : "{}");
    caps::ChainSpec c = spec_from(spec_json ? spec_json : "{}");
    caps::GrowOptions g;
    g.chains = int(j.num("chains", 8));
    c.dp = int(j.num("dp", 20));
    g.density = 0.3;
    g.auto_scale = true;
    g.curve = true;
    g.seed = uint64_t(j.num("seed", 1));
    if (progress) g.progress = [&](int done, int total, int restarts) { return progress(done, total, restarts, user) == 0; };
    caps::GrowReport gr;
    caps::System aa = caps::grow_chains(c, g, &gr);
    caps::RelaxOptions ro;
    ro.target_density = j.num("density", 1.0);
    ro.ftol = 2.0;
    ro.max_iterations = 3000;
    ro.energy = elec(ro.energy);
    caps::RelaxReport rr;
    caps::relax(aa, ro, &rr);
    caps::make_molecules_whole(aa);
    caps::CgMapResult r = caps::cg_map(aa, cg_map_options(j));
    apply_ibi(r, j, cg_map_options(j).temperature);
    r.notes.insert(r.notes.begin(), "reference: " + std::to_string(g.chains) + " chains × " + std::to_string(c.dp) + " units grown and compressed to " +
                                        g6(aa.density()) + " g/cm³ (" + rr.field + ")");
    report_out(cg_map_report(r), report, cap);
    return cg_doc(r, nullptr, "polymer coarse-grained from an all-atom reference melt (" + cg_map_options(j).scheme + ")",
                  {{"chains", std::to_string(g.chains)}, {"dp", std::to_string(c.dp)}, {"density", g6(ro.target_density)}, {"scheme", cg_map_options(j).scheme},
                   {"temperature", g6(cg_map_options(j).temperature)}});
  } catch (const std::exception& e) {
    g_error = e.what();
    return nullptr;
  }
}

// Kremer–Grest → all atoms (design/boards/CoarseGrained "Backmap to all atoms"): one bead per repeat unit. Each KG chain
// gets an all-atom chain of as many units (grown sparse with the polymer builder), whose units — its coarse-graining by
// the unit's backbone atoms per bead — are carried rigidly onto the KG beads (backmap: translated to the bead, turned by
// the fit to the neighbouring beads), in the KG cell; then relaxed (push-off, minimisation) with the default force field.
extern "C" caps_doc* caps_kg_backmap(caps_doc* d, const char* spec_json, const char* options_json, char* report, int32_t cap) {
  try {
    const caps::System kg = d->frame;
    if (kg.atoms.empty()) throw std::runtime_error("no beads");
    const bool kg_model = d->field && !d->field->model.empty() && caps::Json::parse(d->field->model).text("model", "") == "kremer-grest";
    bool named = true;
    for (const auto& a : kg.atoms) named = named && a.name == "KG";
    if (!kg_model && !named) throw std::runtime_error("not a Kremer–Grest melt: build one on the Coarse-grained page (MARTINI beads backmap from Model resolution)");
    if (kg_model && caps::Json::parse(d->field->model).num("sigma", 0) <= 0)
      throw std::runtime_error("the melt is in reduced units: map it to the polymer first (σ in Å, T, bead mass), so the beads sit where the units will");
    const caps::Json o = caps::Json::parse(options_json && *options_json ? options_json : "{}");
    caps::ChainSpec c = spec_from(spec_json ? spec_json : "{}");
    if (c.units.empty()) throw std::runtime_error("choose the repeat unit each bead stands for");
    // the unit's backbone atoms (head to tail along bonds): the atoms per bead of the coarse-graining
    int per_bead = 0;
    {
      const caps::UnitInfo ui = caps::repeat_unit_info(c.units.front().smiles);
      caps::MolGraph g = caps::parse_smiles(c.units.front().smiles);
      if (ui.head < 0 || ui.tail < 0) throw std::runtime_error("the repeat unit needs a head and a tail (two * attachment points)");
      std::vector<int> dist(g.atoms.size(), -1);
      std::vector<int> q{ui.head};
      dist[size_t(ui.head)] = 0;
      for (size_t h = 0; h < q.size(); ++h)
        for (const auto& b : g.bonds) {
          const int u = q[h], w = b.a == u ? b.b : b.b == u ? b.a : -1;
          if (w >= 0 && dist[size_t(w)] < 0) dist[size_t(w)] = dist[size_t(u)] + 1, q.push_back(w);
        }
      per_bead = dist[size_t(ui.tail)] + 1;
      if (per_bead < 1) throw std::runtime_error("the unit's head and tail are not connected");
    }
    // the KG chains: molecules in order, beads in file order within each
    int nm = 0;
    caps::System kgm = kg;
    kgm.has_mol = false;
    const auto mol = kgm.molecules(&nm);
    std::vector<std::vector<uint32_t>> chains(static_cast<size_t>(nm));
    for (uint32_t i = 0; i < kg.atoms.size(); ++i) chains[size_t(mol[i])].push_back(i);
    caps::GrowOptions g;
    g.chains = nm;
    g.seed = uint64_t(o.num("seed", 1));
    g.density = o.num("density", 0.1);   // sparse: the units are placed on the beads afterwards
    g.auto_scale = true;
    g.curve = true;
    c.chain_dp.clear();
    for (const auto& ch : chains) c.chain_dp.push_back(int(ch.size()));
    c.dp = c.chain_dp.front();
    caps::GrowReport gr;
    const caps::System aa = caps::grow_chains(c, g, &gr);
    // the all-atom chains coarse-grained the same way must give the KG chains' bead counts
    caps::ResolutionReport rr0;
    const caps::System cg = caps::coarse_grain(aa, per_bead, &rr0);
    if (cg.atoms.size() != kg.atoms.size())
      throw std::runtime_error("the grown chains coarse-grain to " + std::to_string(cg.atoms.size()) + " beads, the melt has " + std::to_string(kg.atoms.size()) +
                               " (end groups add to the backbone): choose the unit each bead stands for");
    // the KG beads in the order of cg's beads: molecule by molecule
    int cgm = 0;
    caps::System cgs = cg;
    cgs.has_mol = false;
    const auto cmol = cgs.molecules(&cgm);
    std::vector<std::vector<uint32_t>> cg_chains(static_cast<size_t>(cgm));
    for (uint32_t i = 0; i < cg.atoms.size(); ++i) cg_chains[size_t(cmol[i])].push_back(i);
    if (cgm != nm) throw std::runtime_error("the grown chains do not match the melt's chains");
    caps::System beads = cg;
    for (size_t k = 0; k < cg_chains.size(); ++k) {
      if (cg_chains[k].size() != chains[k].size()) throw std::runtime_error("chain " + std::to_string(k + 1) + ": bead counts differ after coarse-graining");
      for (size_t j = 0; j < chains[k].size(); ++j) beads.atoms[cg_chains[k][j]].pos = kg.atoms[chains[k][j]].pos;
    }
    beads.cell = kg.cell;
    caps::BackmapReport rep;
    caps::System s = caps::backmap(aa, beads, per_bead, &rep);
    std::vector<std::string> notes = {std::to_string(nm) + " chains · " + std::to_string(kg.atoms.size()) + " beads → " + std::to_string(s.atoms.size()) +
                                      " atoms (one " + c.units.front().name + " unit per bead, " + std::to_string(per_bead) + " backbone atoms each)"};
    for (const auto& n : rep.notes) notes.push_back(n);
    if (o.num("relax", 1) != 0) {
      caps::RelaxOptions ro;
      ro.ftol = 1.0;
      ro.max_iterations = 5000;
      ro.energy = elec(ro.energy);
      caps::RelaxReport rr;
      caps::relax(s, ro, &rr);
      char b[200];
      std::snprintf(b, sizeof b, "relaxed (%s): E %.1f → %.1f kcal/mol, |F|max %.2f kcal/mol/Å, %.3f g/cm³", rr.field.c_str(), rr.initial.total(), rr.final.total(),
                    rr.fmax_final, s.density());
      notes.push_back(b);
    }
    auto* nd = doc_of_system(std::move(s), d, "model.backmap_kg", "Kremer–Grest beads backmapped to all atoms, one repeat unit per bead",
                             {{"unit", c.units.front().smiles}, {"per_bead", std::to_string(per_bead)}, {"chains", std::to_string(nm)}, {"beads", std::to_string(kg.atoms.size())}});
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
