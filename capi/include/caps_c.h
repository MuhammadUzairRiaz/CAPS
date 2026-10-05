/* CAPS C ABI v19 — the stable boundary used by the Studio (P/Invoke) and other languages.
   Every function is exception-safe: errors are returned as codes and caps_last_error() explains them. */
#ifndef CAPS_C_H
#define CAPS_C_H
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define CAPS_ABI_VERSION 60  /* v60 tags: caps_tags, caps_tag_edit, caps_tag_atoms; v59 caps_react_opts.sites_per_chain, caps_react_sites; v58 layers: caps_layers, caps_set_atom_lock (atom state bit 4); v57 caps_set_atom_state, caps_atom_states (hide and ghost atoms in the view); v56 caps_project_indices; v55 caps_export_scene (POV-Ray, glTF, OBJ), caps_export_image engine raytrace; v54 caps_scene triangles (surfaces and polyhedra on the GPU; cpu_only only for the colour-vision preview); v53 caps_pick_at (picking by the view ray, no render); v52 caps_field_groups_by_example (typing by hand inside a by-group assignment), LAMMPS inputs list each group's types with element and mass; v51 water models (caps_water_models; edit op water_model; field groups "molecules": "water" and "water": model); v50 caps_open_frames (read a frame selection); v49 caps_camera.fov_deg; v48 labels: caps_atom_labels kinds of caps_label_kinds, caps_bond_labels; v2 relax, field; v3 md, trajectory; v4 equilibrate, chains; v5 pack; v6 react; v7 CAPS Field; v8 Analyze; v9 mechanics, Tg; v10 LAMMPS input; v11 convergence checks; v12 molecule builder; v13 palette, threads; v14 bench; v15 polymer builder; v16 electrostatics; v17 surfaces, interfaces, held molecule, inserted curatives; v18 progressive open, keyboard focus; v19 ambient occlusion, view scale; v20 space groups, crystal builder, peptides, solvation, appearance, trajectory player, torsion scan, editing, selections; v21 r-RESPA (caps_md_opts.respa), reactions during MD (caps_react_opts.during_md), restraints; v22 GROMACS export (caps_gromacs), χ from pair contacts (caps_chi_contacts); v23 export center (caps_export_engines); v24 coarse-grained beads (caps_build_beads, caps_bead_templates); v25 live view of MD and equilibration (caps_set_live); v26 GPU view (caps_render_scene, caps_view_fit); v27 the scene carries its camera-fit inputs (a view turns while a run holds the document); v28 caps_shadow (a copy of the shown frame the window reads while a run holds the document); v29 bond constraints (caps_md_opts / caps_equil_opts .constraints: SHAKE/RATTLE), typing by example; v30 relax push-off by MD with a ramped force cap (caps_relax_opts.pushoff_ramp_ps …); v31 caps_equil_opts.tol_internal (the internal-distance convergence check), caps_pipeline_export_grid; v32 LINCS (caps_md_opts / caps_equil_opts .constraint_algorithm), an internal-distance target curve (caps_equil_opts.internal_target); v33 CBMC regrowth (caps_cbmc); v34 adsorption locator (caps_adsorption), sorption (caps_sorption); v35 layer stacks (caps_stack_documents), caps_frame_copy, pipeline outputs (caps_pipeline_write_outputs); v36 relax etol / pressure_tol, MD per-axis pressure coupling, fixed atoms (caps_set_fixed_atoms), caps_energy_terms, caps_analyze_opts.group, caps_field_assign_groups; v37 caps_molecule_ids; v38 caps_analyze_opts.radii, electron scattering (analyze id electron); v39 caps_martini_melt, Kremer–Grest melts carry their force field; v40 caps_build_opts.rotor_search, heavy_only; v41 AMBER prmtop topologies carry their force field (caps_field_assign "file", caps_field_file_available); v42 caps_react_cycle.max_force, caps_react returns 2 when a failed cycle kept the completed ones, caps_kg_backmap; v43 caps_cg_map, caps_cg_from_polymer (structure-based CG); v44 caps_analyze_opts.zbin, axis, surface; v47 caps_field_set_mixing; v46 caps_mech_opts tg barostat, tau_t / tau_p, average_from, tg_property, tg_fit, glassy_max, rubbery_min; v45 React with the assigned force field, between chains, crosslink targets, byproducts, weights, auto capture (caps_react_opts), caps_react_summary */

typedef struct caps_doc caps_doc;   /* an opened file: trajectory + current frame + renderer */

typedef struct {
  double yaw, pitch, zoom, pan_x, pan_y;
  int32_t perspective;
  double fov_deg;   /* ABI 49: the perspective view's field of view, degrees (0: 35; 10 … 120) */
} caps_camera;

typedef struct {
  int32_t width, height, supersample;
  int32_t background;     /* 0 dark, 1 white, 2 transparent, 3 custom */
  uint32_t custom_rgb;
  int32_t colour_by;      /* 0 element, 1 molecule, 2 type, 3 property (distance to own molecule's centre) */
  int32_t style;          /* 0 ball & stick, 1 space filling, 2 sticks, 3 no hydrogens, 4 backbone */
  int32_t outlines, depth_cue, show_cell;   /* outlines: 0 none, 1 light, 2 high contrast */
  int32_t highlight[4];   /* up to four selected atom indices, -1 for unused */
  int32_t focus;          /* v18: atom index + 1 drawn with the keyboard-focus ring, 0 for none */
  int32_t ambient_occlusion;   /* v19: darken atoms by the open sky they see */
  double lod_near, lod_far;    /* v20 level of detail: full within lod_near Å of the focus, spheres to lod_far, points beyond; 0 off */
} caps_render_opts;

typedef struct {
  int64_t id, mol;
  int32_t type, element, index;
  double charge, x, y, z;
  char element_symbol[4];
  char name[16];
} caps_atom_info;

typedef struct {
  int64_t atoms, bonds, molecules, frames;
  int32_t bonds_from_file, has_charges, cell_valid, unwrapped;
  double cell_a, cell_b, cell_c, volume, density, total_mass, total_charge;
  char format[24];
} caps_summary;

typedef struct {
  int32_t chains, dp, tacticity;   /* tacticity: 0 atactic, 1 isotactic, 2 syndiotactic */
  uint64_t seed;
  double box, density;             /* box edge in Å, or 0 to use density (g/cm^3) */
  double contact_scale;            /* 1 = full contact limits; negative (caps_grow_chains): start at |value| and lower it when crowded */
  int32_t curve;                   /* allow gauche backbone torsions */
} caps_grow_opts;

/* Progress: (chains finished, chains, restarts, user) -> 0 to continue, non-zero to cancel. */
typedef int32_t (*caps_progress_fn)(int32_t done, int32_t total, int32_t restarts, void* user);

int32_t caps_abi_version(void);
const char* caps_last_error(void);

caps_doc* caps_open(const char* path, const char* topology_path);   /* NULL on error */
/* v48 several trajectory files of one run as one trajectory (JSON array of paths; open_files in io.hpp): parts in time
   order, a frame repeated at a boundary kept once, the atoms checked to agree. NULL on error. */
caps_doc* caps_open_many(const char* paths_json, const char* topology_path);
/* v19: what a file holds before opening it (format, first lines, dump columns, types, frames, bonds) as JSON. */
int32_t caps_inspect_file(const char* path, const char* topology_path, char* json, int32_t cap);
/* Staged open (Studio progressive open): stage 0 format detected, 1 frame 0 read, 2 topology joined, 3 frames read
   (fraction of the file); return non-zero to stop — the frames read so far are kept. max_frames > 0 stops early. */
typedef int32_t (*caps_open_progress_fn)(int32_t stage, double fraction, const char* detail, void* user);
caps_doc* caps_open_staged(const char* path, const char* topology_path, int32_t max_frames, caps_open_progress_fn progress, void* user);
/* v50 the same, keeping only file frames first, first + stride, … up to last (inclusive; last < 0: to the end). A LAMMPS
   dump passes over the other frames unread and stops after last; binary trajectories read and drop them; multi-frame
   .gro/.pdb/.xyz are thinned once read. The first kept frame is frame 0; a note says what was kept. NULL on error
   (a selection past the file's end included). */
caps_doc* caps_open_frames(const char* path, const char* topology_path, int64_t first, int64_t last, int64_t stride, int32_t max_frames,
                           caps_open_progress_fn progress, void* user);
/* v20 import (design/boards/ImportDialog). options: {"bonds":"perceive"|"file"|"none","tolerance":0.45,"bond_orders":true,
   "split":true,"unwrap":true,"use_cell":true}. caps_import_preview reads frame 0 only and returns {ok, format, format_name,
   units, bytes, head[], atoms, bonds_in_file, bonds, molecules, single, double, triple, aromatic, cell, fragment_heavy,
   fragment_atoms, fragment_bonds, notes[]}; caps_import_fragment is that preview's fragment (the first ten connected
   heavy atoms, carbons first, with their hydrogens) as a document. */
caps_doc* caps_import(const char* path, const char* topology_path, const char* options_json);
int32_t caps_import_preview(const char* path, const char* options_json, char* json, int32_t cap);
caps_doc* caps_import_fragment(const char* path, const char* options_json);
/* Moves the frames of src (the same file read in full) into dst, keeping dst's selection, field and current frame.
   Returns dst's frame count, or -1 when the atom counts differ. src is left empty. */
int32_t caps_adopt_frames(caps_doc* dst, caps_doc* src);
/* A new document holding d's shown frame (topology, positions, cell), its wrap, appearance, selection, supercell and
   field: what the window reads while a run holds d. NULL on error. */
caps_doc* caps_shadow(caps_doc* d);
void caps_close(caps_doc* d);

/* Grow polystyrene chains in a periodic cell; the result is a new document. NULL on error or cancel.
   The report (restarts, backtracks, worst contact margin) is copied into `report` when given. */
caps_doc* caps_grow(const caps_grow_opts* o, caps_progress_fn progress, void* user, char* report, int32_t report_cap);

typedef struct {
  int32_t method;                  /* 0 steepest descent, 1 Polak–Ribière CG, 2 L-BFGS, 3 FIRE */
  double ftol;                     /* largest atomic force to stop at, kcal/mol/Å */
  int32_t max_iterations;          /* per minimisation */
  double target_density;           /* g/cm^3 to compress (or expand) to; 0 keeps the cell */
  double compress_step;            /* relative density change per stage, e.g. 0.06 */
  int32_t pushoff;                 /* capped-force push-off before minimising */
  int32_t relax_box;               /* isotropic box relaxation towards `pressure` */
  double pressure;                 /* atm */
  double cutoff;                   /* Å, LJ and Coulomb */
  int32_t coulomb;                 /* damped shifted force electrostatics */
  int32_t threads;                 /* 0 = automatic */
  int32_t box_anisotropic;         /* ABI 29: 1 each axis follows its own diagonal pressure (orthorhombic cells) */
  int32_t box_axes;                /* ABI 29: axes that move, bits 1 x, 2 y, 4 z (0 = all) */
  double pushoff_ramp_ps;          /* ABI 30: > 0 push-off by NVT MD first, the LJ force cap raised over this time (λ ramp) */
  double pushoff_cap;              /* ABI 30: final force cap, kcal/mol/Å (0 = 500) */
  double pushoff_temperature;      /* ABI 30: K (0 = 300) */
  double etol;                     /* ABI 36: stop when the relative energy change per step is below this (0 = 1e-8) */
  double pressure_tol;             /* ABI 36: box relaxation stops when |P − pressure| is below this, atm (0 = 100) */
} caps_relax_opts;

/* Relax progress: (stage, stages, iteration, energy kcal/mol, largest force, density, user) -> non-zero cancels. */
typedef int32_t (*caps_relax_progress_fn)(int32_t stage, int32_t stages, int32_t iteration, double energy, double fmax, double density,
                                          void* user);

/* Configurational-bias Monte Carlo regrowth of chain ends (v33; Siepmann & Frenkel 1992): moves attempted, torsion trials
   per bond, rotatable backbone bonds regrown per move at most, temperature (K), pair cut-off (Å), damped shifted force
   electrostatics, seed. Trial energies come from the Field assignment (or the built-in force field). */
typedef struct caps_cbmc_opts {
  int32_t moves, trials, max_torsions;
  double temperature, cutoff;
  int32_t coulomb;
  uint64_t seed;
} caps_cbmc_opts;

/* CBMC progress: (moves done, moves, accepted, user) -> non-zero cancels. */
typedef int32_t (*caps_cbmc_progress_fn)(int32_t done, int32_t moves, int32_t accepted, void* user);

/* Regrow chain ends of the current frame. On success the document holds the start and snapshots through the run (the
   last is the result) and report gets JSON {attempted, accepted, acceptance, chains, r2_before, r2_after, energy_change,
   cutoff, seconds, notes[]}. Returns 0, or -1 on error / cancel (the document is unchanged). */
int32_t caps_cbmc(caps_doc* d, const caps_cbmc_opts* o, caps_cbmc_progress_fn progress, void* user, char* report, int32_t report_cap);

/* Relax the current frame with GAFF (C/H in this version). On success the document holds one frame per stage
   (the start, each push-off / compression stage, the final structure) and shows the last one. Returns 0, 1 when
   the run finished without meeting the force tolerance, or -1 on error / cancel (the document is unchanged). */
int32_t caps_relax(caps_doc* d, const caps_relax_opts* o, caps_relax_progress_fn progress, void* user, char* report, int32_t report_cap);

typedef struct {
  double dt;                       /* fs */
  int64_t steps;
  double temperature;              /* K */
  int32_t thermostat;              /* 0 none (NVE), 1 Bussi, 2 Langevin, 3 Nosé–Hoover chain (LAMMPS fix nvt) */
  double tau_t;                    /* fs */
  int32_t barostat;                /* 0 none, 1 stochastic cell rescaling, 2 Berendsen, 3 MTK (isotropic, LAMMPS fix npt iso;
                                      with the Nosé–Hoover thermostat or none, no constraints) */
  double pressure, tau_p;          /* atm, fs */
  int32_t new_velocities;          /* draw velocities even when the document has them */
  uint64_t seed;
  int32_t thermo_every, frame_every;
  double cutoff;
  int32_t coulomb, tail, threads;
  int32_t respa;                   /* r-RESPA inner steps: bonded forces every dt / respa (0 or 1: off) */
  int32_t constraints;             /* 0 none, 1 bonds to hydrogen and rigid water, 2 all bonds (SHAKE/RATTLE) */
  int64_t step_offset;             /* added to reported steps and times (a run continued from a checkpoint) */
  int64_t checkpoint_every;        /* steps between checkpoints (0: automatic, about 50 per run; < 0: none) */
  int32_t constraint_algorithm;    /* ABI 32: 0 SHAKE (positions) / RATTLE, 1 LINCS (positions) / RATTLE */
  int32_t box_anisotropic;         /* ABI 36: 1 each axis in box_axes scaled on its own from P_kk (Berendsen barostat) */
  int32_t box_axes;                /* ABI 36: bits 1 x, 2 y, 4 z (0 = all) */
} caps_md_opts;

typedef struct {
  int64_t step;
  double time_ps, temperature, potential, kinetic, total, conserved, pressure, volume, density;
} caps_thermo;

/* MD progress: one thermo row per call -> non-zero cancels. */
typedef int32_t (*caps_md_progress_fn)(const caps_thermo* row, int64_t steps, void* user);

/* Molecular dynamics from the current frame. Velocities carry over when the current frame is the end of a previous
   run in this document, otherwise they are drawn at the temperature. On success the document holds the recorded
   frames (start, every frame_every steps, end) and shows the last one. Returns 0, or -1 on error / cancel (the
   document is unchanged). */
int32_t caps_md(caps_doc* d, const caps_md_opts* o, caps_md_progress_fn progress, void* user, char* report, int32_t report_cap);
/* The last checkpoint of a caps_md / caps_equilibrate run in this document, kept when the run fails or is stopped
   (ABI 29): JSON {op: "info" | "restore" | "clear"} → {ok, has, kind: "md" | "equilibrate", step, steps, time_ps,
   atoms, ended: "failed" | "stopped" | "finished", error}. restore appends the checkpoint's positions and cell as the
   last frame with its velocities, so a caps_md with new_velocities 0 continues the run exactly. */
int32_t caps_checkpoint(caps_doc* d, const char* json, char* out, int32_t cap);

/* Write every frame of the document, the format by the extension: .lammpstrj / .dump (LAMMPS text dump, id mol type
   xu yu zu), .dcd, .xyz (extended XYZ), .pdb (MODEL per frame), .gro (frames one after another), .trr (GROMACS). */
int32_t caps_save_trajectory(caps_doc* d, const char* path);
/* v47: the per-atom columns the trajectory carries beyond positions (a LAMMPS dump's fx, c_pe, v_…; |f| and |v| from
   the components) with their range in the frame shown: JSON {columns: [{name, min, max}], velocities}. Colour by one
   with caps_set_appearance {"colour": "column:NAME"}. Returns the length needed or -1. */
int32_t caps_trajectory_columns(caps_doc* d, char* json, int32_t cap);

typedef struct {
  double t_final, t_max, p_final, p_max;   /* K, atm */
  double time_scale;                       /* multiplies every stage duration */
  int32_t cycles;                          /* annealing */
  double t_low, t_high, ramp_ps, hold_ps;
} caps_protocol_params;

/* Text of a named protocol ("larsen21", "annealing", "pushoff"), one stage per line; returns its length or -1. */
int32_t caps_protocol_text(const char* name, const caps_protocol_params* p, char* out, int32_t cap);

typedef struct {
  double dt;                               /* fs */
  int32_t thermostat;                      /* 1 Bussi, 2 Langevin, 3 Nosé–Hoover chain */
  int32_t barostat;                        /* 1 stochastic cell rescaling, 2 Berendsen, 3 MTK (with Nosé–Hoover) */
  double tau_t, tau_p;                     /* fs */
  uint64_t seed;
  double cutoff;
  int32_t coulomb, tail, threads;
  double frame_ps, thermo_ps;
  int32_t until_converged;                 /* production NPT blocks after the stages until the checks pass */
  double block_ps;
  int32_t max_blocks;
  double tol_density, tol_energy, tol_rg;  /* relative, kcal/mol per atom, relative */
  int32_t constraints;                     /* 0 none, 1 bonds to hydrogen and rigid water, 2 all bonds (SHAKE/RATTLE) */
  double tol_internal;                     /* ABI 31: internal distances ⟨R²(n)⟩/(n⟨b²⟩), relative change at any n (0: 5 %) */
  int32_t constraint_algorithm;            /* ABI 32: 0 SHAKE, 1 LINCS */
  const double* internal_target;           /* ABI 32: a target ⟨R²(n)⟩/(n⟨b²⟩) curve indexed by n (e.g. RIS C_n; 0 entries skipped), or NULL */
  int32_t internal_target_n;               /* its length (n = 0 … internal_target_n − 1) */
} caps_equil_opts;

/* Equilibrate progress: (stage, stages, stage label, thermo row, user) -> non-zero cancels. */
typedef int32_t (*caps_equil_progress_fn)(int32_t stage, int32_t stages, const char* label, const caps_thermo* row, void* user);

/* Run a protocol (text as from caps_protocol_text, or written by hand) from the current frame. On success the document
   holds the recorded frames and shows the last one; returns 0, 1 when the convergence checks did not pass, -1 on
   error / cancel (the document is unchanged). */
/* The convergence checks of the last caps_equilibrate run with until_converged, as JSON:
   {"converged", "blocks", "block_ps", "checks": [{"quantity", "ok", "change", "tolerance", "blocks": [block means]}]}.
   Returns the length needed including the final NUL (json = NULL to size); empty before a run. */
int32_t caps_equilibrate_checks(caps_doc* d, char* json, int32_t cap);
int32_t caps_equilibrate(caps_doc* d, const char* protocol, const caps_equil_opts* o, caps_equil_progress_fn progress, void* user,
                         char* report, int32_t report_cap);

/* Mean-square internal distances of the backbones in the current frame: n and <R²(n)>/(n <b²>). Returns the count
   written; chains and <b²> are returned through the pointers when given. */
int32_t caps_internal_distances(caps_doc* d, int32_t* n, double* ratio, int32_t cap, int32_t* chains, double* b2);

/* Pack progress: (round, rounds, penalty, molecules in violation, user) -> non-zero cancels. */
typedef int32_t (*caps_pack_progress_fn)(int32_t loop, int32_t loops, double penalty, int32_t bad, void* user);

/* Pack from packmol-format text (structure paths relative to base_dir). Returns a new document, or NULL when the
   input is wrong, the run is cancelled or the tolerance cannot be met (nothing below tolerance is ever returned);
   the report (distances, rounds, time) is filled in both cases. */
caps_doc* caps_pack(const char* text, const char* base_dir, int32_t threads, caps_pack_progress_fn progress, void* user, char* report,
                    int32_t report_cap);

/* Text of a built-in reaction template ("cc_crosslink", "epoxy_amine_primary", "epoxy_amine_secondary");
   an empty name gives the list of names, one per line. Returns the length or -1. */
int32_t caps_reaction_template(const char* name, char* out, int32_t cap);

typedef struct {
  uint64_t seed;
  int32_t max_cycles, max_per_cycle;
  double target_conversion;
  double capture;                  /* > 0 overrides every template's capture distance, Å */
  int32_t relax;                   /* minimise after each cycle (needs typed atoms) */
  int32_t relax_iterations;
  double md_ps, temperature;       /* NVT between cycles */
  double cutoff;
  int32_t coulomb;
  int32_t during_md;               /* REACTER-style: continuous NVT, reactions checked every md_ps, sites stabilised locally */
  /* v45 */
  int32_t field_mode;              /* 0: the Field assignment (when there is one) types and parameterises the structure after
                                      every cycle and is re-assigned to the product; 1: the built-in default during the run, the
                                      assignment dropped (the old behaviour) */
  int32_t between_chains;          /* bonds only between different chains (a crosslinker belongs to the chains it has joined) */
  int32_t keep_byproducts;         /* a template's byproduct atoms (H2, H2O) stay as molecules; 0: removed */
  int32_t selection;               /* several templates: 0 closest pairs first; 1 by relative weights */
  const char* weights;             /* "2,1,0.5": a weight per template in order ("" / NULL: the templates' own) */
  int32_t auto_capture;            /* no pair found: every capture grows by capture_step up to capture_max */
  double capture_max, capture_step; /* Å (0: 8 and 0.5) */
  int32_t target_kind;             /* 0 conversion (target_conversion), 1 links between chains, 2 links per chain, 3 crosslink
                                      density (mol/m³), 4 molecular weight between crosslinks Mc (g/mol), 5 degree of
                                      crosslinking DC = 2 links / monomers × 100 (%) */
  double target_value;
  /* v59 */
  int32_t sites_per_chain;         /* at most this many of each chain's reactive sites react (0: no limit) */
} caps_react_opts;

/* v59: each chain's reactive sites for the templates (text, or built-in names): {chains: [{chain, sites, units, atoms, mass}],
   total_sites, chains_n, units, mass, repeat_unit_mass (g/mol: the chains' mass over their repeat units)}. */
int32_t caps_react_sites(caps_doc* d, const char* templates_text, char* json, int32_t cap);

typedef struct {
  int32_t cycle, reactions, total, clusters, atoms;
  double conversion, largest_fraction, reduced_mw, energy;
  double max_force;                /* ABI 42: largest force after the cycle's relaxation, kcal/mol/Å (0: not relaxed) */
  int32_t crosslinks;              /* v45: links between chains so far */
  double capture;                  /* v45: the capture distance the cycle used, Å */
  int32_t target;                  /* v45: the crosslink target as links (0: a conversion target) */
  double density, degree;          /* v45: crosslink density so far (mol/m³), degree of crosslinking so far (%) */
} caps_react_cycle;

typedef int32_t (*caps_react_progress_fn)(const caps_react_cycle* row, void* user);

/* React the current frame with the templates in `templates` (text). On success the document holds one frame per
   cycle and shows the last; returns 0, or -1 on error / cancel (the document is unchanged). ABI 42: 2 when a cycle
   failed after others completed — the document holds the structure after the last completed cycle, the report says
   which cycle failed and why ("failed at cycle N: …"), and caps_last_error gives the reason. */
int32_t caps_react(caps_doc* d, const char* templates, const caps_react_opts* o, caps_react_progress_fn progress, void* user, char* report,
                   int32_t report_cap);
/* With caps_set_live on d, every cycle sends the structure as a snapshot: each atom's molecule id is the chain it started
   in (colour by molecule keeps every chain its colour as they join; byproducts get their own ids), the atoms of the links
   formed between chains are selected, and the stats JSON is {cycle, reactions, crosslinks, target, density, degree,
   conversion, atoms}. */
/* v45: the last caps_react run on this document as JSON {chains, crosslinks, intrachain, byproducts, target, volume,
   chain_mass, density (mol/m³), per_chain, mc (g/mol), reactions, initial_sites, conversion, field, field_after,
   notes[]}. Returns the length needed or -1. */
int32_t caps_react_summary(caps_doc* d, char* json, int32_t cap);
/* v45: the reactions as a LAMMPS fix bond/react set in dir: STEM.data (every type the reactions create, with its
   coefficients), STEM.in (the force field's styles, the molecule templates, fix bond/react with stabilisation, NVT, the
   reaction counts in thermo), and per template and chemical environment STEM_<name>_<k>_pre.mol, _post.mol, _map.txt —
   cut from real reaction sites of the current frame and typed with the assigned force field before and after the
   reaction (a complete assignment is needed). options JSON: {stem, radius (3 bonds), variants (6), keep_byproducts,
   between_chains (molecule inter), weights [..], nevery (100), rmax (3.5 Å; 0: the template's capture), temperature (300), steps (100000), seed}. Report JSON
   {files[], notes[], variants[{reaction, name, sites, pre_atoms, edge, deleted}], candidates, covered}. Returns the
   length needed or -1. */
int32_t caps_bond_react_export(caps_doc* d, const char* templates, const char* dir, const char* options, char* report, int32_t cap);
/* v45: a LAMMPS fix bond/react set read back as a CAPS reaction template: pre- and post-reaction molecule files and the
   map file (CAPS's or anyone's); masses_from: a data file whose Masses give the elements of numbered types (NULL when
   the molecule files carry Masses or element-like type labels). JSON {text, notes[]}: text is the template for
   caps_react (form, break, move, delete and byproducts from comparing the two templates). Returns the length or -1. */
/* v45: the reaction library (data/reactions/library.json at path): every scheme with its CAPS template. JSON {reactions:
   [{id, name, category, description, reactants[{label, smiles}], products[...], tags[{map, role}], note, steps[ids],
   template (text, when it converts), error (why not), notes[]}]}. Returns the length needed or -1. */
int32_t caps_reaction_library(const char* path, char* json, int32_t cap);
int32_t caps_bond_react_import(const char* pre, const char* post, const char* map, const char* masses_from, const char* name, double capture,
                               char* json, int32_t cap);

/* Force-field summary of the current frame (the Field assignment, else the built-in GAFF typing of C and H): types,
   term counts and energy terms, as text. Returns 0 or -1. */
int32_t caps_field_info(caps_doc* d, char* text, int32_t cap);
/* The current frame's energy by term with the force field runs use (v36): JSON {bond, angle, dihedral, improper, vdw
   (with the tail term), coulomb, total} in kcal/mol, and the electrostatics used. Returns the length needed with NUL. */
int32_t caps_energy_terms(caps_doc* d, char* json, int32_t cap);

/* CAPS Field. Assigns a force field (a caps-forcefield JSON file; typed by its own rules, or by rules_path when given)
   to the document: types every atom, sets charges (charges: 0 from the force field, 1 Gasteiger, 2 keep the file's, 3 QEq, 4 automatic: the force
   field's when its types carry charges, else Gasteiger, noted in the report; 5 bond increments keyed by each type's
   number, from the force field's own table or the one its "charge_increments_from" names — OPLS-AA 2024 → OPLS 2005)
   and looks up every parameter. The types are written into the document (colour by type shows them). Relax,
   Dynamics, Equilibrate and LAMMPS data then use this force field; while atoms are untyped or parameters missing they
   refuse to run (CAPS never guesses parameters). Returns 0 when complete, 1 when something is missing, -1 on error. */
int32_t caps_field_assign(caps_doc* d, const char* ff_path, const char* rules_path, int32_t charges);
/* A force field per group of molecules (v36; ffmerge.hpp): {groups: [{name, molecules: "1" | "2-10, 12" | "rest",
   forcefield: path, charges: as caps_field_assign}], eps_rule: geometric | arithmetic, sigma_rule: arithmetic | geometric
   | sixthpower, scaling14: refuse | first, cross96: refuse | rmin, pairs: [{a, b, eps, sigma}] (merged type names)}.
   Each group is typed and parameterised on its own, the parts merged with the cross pairs written out explicitly (the
   LAMMPS, GROMACS and DL_POLY files carry them). The report is as caps_field_report's, with "groups". Returns as
   caps_field_assign; what one simulation cannot hold (different 1-4 scalings, 9-6 with 12-6) is refused with the reason.
   A group may take a literature many-body potential instead of a force field (manybody.hpp): potential: {style: tersoff |
   tersoff/mod | tersoff/mod/c | tersoff/zbl | sw | vashishta | gw | gw/zbl | eam/alloy | eam/fs, file: path, units: metal |
   real (only when the file's first line does not say)} — one type per element, standard masses, no charge, UFF
   Lennard-Jones across; LAMMPS files only (pair_style hybrid/overlay, the file written beside them), GROMACS and DL_POLY
   refused; CAPS's relax / MD / equilibrate refuse until those atoms are held. */
int32_t caps_field_assign_groups(caps_doc* d, const char* json);
/* ABI 43: structure-based coarse-graining. caps_cg_map: the document (every frame) mapped to beads — options
   {"scheme": "unit" | "backbone_side" | "backbone_n", "per_bead", "temperature"}; caps_cg_from_polymer: an all-atom
   reference melt of the polymer (spec JSON; options also "chains", "dp", "density", "seed") grown, compressed and mapped.
   Bonds and angles Boltzmann-inverted, a repulsive WCA from the intermolecular bead g(r); the new document carries the
   bead model. Report JSON {ok, notes, bonds, angles, sigma, cut, epsilon, gr}. NULL on error (caps_last_error). */
caps_doc* caps_cg_map(caps_doc* d, const char* options_json, char* report, int32_t cap);
caps_doc* caps_cg_from_polymer(const char* spec_json, const char* options_json, caps_progress_fn progress, void* user, char* report, int32_t cap);
/* ABI 42: Kremer–Grest beads → all atoms, one repeat unit (spec_json: a polymer spec, its first unit) per bead: all-atom
   chains grown with the melt's bead counts, carried rigidly onto the beads in the melt's cell, relaxed unless
   options {"relax": 0}; options also {"seed", "density"}. A new document, or NULL (caps_last_error). */
caps_doc* caps_kg_backmap(caps_doc* d, const char* spec_json, const char* options_json, char* report, int32_t cap);
/* ABI 41: 1 when the structure's own file carried a force field (an AMBER prmtop) that still fits its atoms — caps_field_assign(d, "file", …) assigns it; else 0 */
int32_t caps_field_file_available(const caps_doc* d);
/* The assignment as JSON (atoms with type, rule, candidates, charge; counts; missing terms; types present with the
   viewer's colours; all types of the force field; notes; energy). Returns the length needed including the final NUL;
   call with json = NULL to size the buffer. Empty when nothing is assigned. */
int32_t caps_field_report(caps_doc* d, char* json, int32_t cap);
/* Which library force fields describe the current structure (ABI 22): dir is the library folder (catalogue.json).
   Each force field's typing is tried, then its parameters: {ok, forcefields: [{id, name, status (complete | untyped atoms
   | missing parameters | charges do not balance | no typing rules | error: …), untyped, untyped_groups: [{environment,
   count, atoms[]}], missing[], missing_count, charges (types | gasteiger | none), net_charge, balanced, complete}]}.
   balanced: the force field's own charges leave the structure at its formal charge. Seconds for a few thousand atoms;
   progress gets (force field, fraction) and returns non-zero to stop. */
typedef int32_t (*caps_stage_fn)(const char* stage, double fraction, void* user);
int32_t caps_field_coverage(caps_doc* d, const char* dir, caps_stage_fn progress, void* user, char* out, int32_t cap);

/* Adsorption locator (v34, adsorption.hpp): JSON {adsorbates: [{smiles, count}] (added after the substrate; or
   first_atom: the atom from which the existing molecules move), cycles, steps, t_high, t_low (K), cutoff, coulomb,
   region: "cell" | "above" (over the substrate's top face) or z_lo/z_hi, keep, seed}. The Field assignment is applied
   to the whole (else the built-in force field). The document gets the kept configurations as frames, the lowest last.
   → {ok, error, adsorption_energy, adsorbate_substrate, adsorbate_adsorbate, components: [{name, molecules, de_dn}],
   configs: [{energy, cycle}], hist_edges, hist_counts, acceptance, steps, seconds, z_lo, z_hi, forcefield, notes}. */
int32_t caps_adsorption(caps_doc* d, const char* json, caps_stage_fn progress, void* user, char* out, int32_t cap);
/* Sorption (v34, sorption.hpp): JSON {sorbate (SMILES), temperature, insertions (Widom), pressures_kpa: [..] (GCMC),
   steps, cutoff, coulomb, seed}; the structure is the fixed host (typed with its Field assignment plus the sorbate, else
   the built-in force field); the document is not changed. → {ok, error, widom_w, widom_error, mu_ex, henry_mol_kg_kpa,
   solubility (cm³(STP)/(cm³ atm)), host_mass, volume, forcefield, isotherm: [{pressure_kpa, loading (per cell),
   loading_error, mol_per_kg, cm3stp_per_cm3, heat (kcal/mol), acceptance_insert, acceptance_delete}], notes}. */
int32_t caps_sorption(caps_doc* d, const char* json, caps_stage_fn progress, void* user, char* out, int32_t cap);
/* DPD (v34, dpd.hpp): JSON {species: [{name, sequence ("AAAAABBBBB"), count}], density, chi: {"AB": χ}, a: {"AB": a_ij},
   gamma, dt, bond_k, steps, equilibration, frame_every, rc_angstrom, seed} → a new document with the frames (beads as
   atoms: A C, B O, C N, D S, others P) and report JSON {ok, error, beads, molecules, box, kT, kT_error, pressure,
   pressure_error, order (ψ), q_peak, spacing, q[], sq[], order_series [[step, ψ]], types, notes, seconds}; NULL on error. */
caps_doc* caps_dpd(const char* json, caps_stage_fn progress, void* user, char* out, int32_t cap);
/* Sets atom `index` (0-based) to a type by hand; type "" or NULL removes the override. Returns as caps_field_assign. */
int32_t caps_field_override(caps_doc* d, int32_t index, const char* type);
/* Typing by example (ABI 29): types given on an example document (a repeat unit's head, body and tail in a short chain,
   a copolymer's units and junctions; types_json: one type name per example atom, "" = not given) are learned by each
   atom's chemical environment and set as overrides on every atom of `d` whose environment the example has (d's Field
   must be assigned with the force field the types belong to). Report: {"radius", "environments", "exact", "shorter",
   "unmatched", "set", "conflicts": [...], "unknown_types": [...]}. Returns as caps_field_assign, -1 on error. */
int32_t caps_field_type_by_example(caps_doc* d, caps_doc* example, const char* types_json, char* report, int32_t cap);
/* v52: the types taught by hand on an example chain (as caps_field_type_by_example) learned onto every group of a by-group
   assignment whose force field is ff_path; the other groups keep theirs (a potential file, a water model), the cross terms
   are merged again, and the teaching is kept for later group assignments of this document (example NULL forgets it).
   report: {"set", "unmatched", "radius", "environments", "conflicts", "unknown_types", "groups"}. Returns 0 complete, 1 some
   parameters missing, -1 error. groups_json: the groups to assign (as caps_field_assign_groups); NULL or "" the assigned ones. */
int32_t caps_field_groups_by_example(caps_doc* d, caps_doc* example, const char* types_json, const char* ff_path, const char* groups_json, char* report,
                                     int32_t cap);
/* Atoms of the current frame with the same chemical environment as `atom` to `radius` bonds (JSON array of indices,
   `atom` included): where one hand-assigned type applies. Returns the length needed or -1. */
int32_t caps_equivalent_atoms(caps_doc* d, int32_t atom, int32_t radius, char* json, int32_t cap);
/* Adds a parameter rule entered by hand: kind pair | bond | angle | dihedral | improper, the atom types (space-
   separated, as the missing-term list names them), a style ("" for the force field's default) and the parameters in
   that style's order. Such terms are reported as estimated. Returns as caps_field_assign. */
int32_t caps_field_add_rule(caps_doc* d, const char* kind, const char* types, const char* style, const char* params);
/* Merges parameters from another file over the force field: caps-forcefield .json, moltemplate .lt, AMBER frcmod (any
   name containing frcmod, or .dat) or the [ *types ] sections of a GROMACS .itp / .top. */
int32_t caps_field_import(caps_doc* d, const char* path);
/* The same with options (ABI 29): {mode: "override" (the default: the file's rules win) | "fill" (the file's rules are
   used only where the force field defines nothing — borrowed parameters, each listed in the report as filled_terms;
   for a moltemplate OPLS-AA force field the file's atom classes become its type patterns)}. */
int32_t caps_field_import_ex(caps_doc* d, const char* path, const char* options_json);
/* v47: the mixing rule for unlike Lennard-Jones pairs in place of the force field's own: "arithmetic" (Lorentz–
   Berthelot), "geometric" or "sixthpower" (Waldman–Hagler); "" back to the force field's. Every run and export uses it
   (explicit pair parameters still win) and the report says it is not the published rule. Returns as caps_field_assign. */
int32_t caps_field_set_mixing(caps_doc* d, const char* rule);
/* Removes every imported and hand-entered parameter. */
int32_t caps_field_remove_rules(caps_doc* d);
/* Ends the assignment and restores the file's own types and charges. */
int32_t caps_field_clear(caps_doc* d);
/* Writes the assigned types, one per line (for caps ff apply --types). */
int32_t caps_field_types_file(caps_doc* d, const char* path);

/* Analyze: properties over the frames of the document (ids as for caps analyze: density, rdf, sq, xray, electron, neutron, rg,
   ree, cn, persistence, msd, diffusion, relaxation, ced, delta, ffv, psd; comma-separated). Zero or negative fields
   of the options take the defaults. Cohesive energy uses the Field assignment when there is one (it must be complete),
   else the built-in GAFF typing of C and H. The progress callback gets what is being computed and the fraction done;
   return 0 to cancel. Returns 0, or -1 on error / cancel; the result is read with caps_analyze_report. */
typedef struct {
  int64_t first, last, stride;       /* frames (last: -1 = the last; 0 is taken as -1) */
  double frame_ps, timestep_fs;      /* time between frames; else timesteps × timestep_fs */
  int32_t blocks;                    /* block averages for the errors (5) */
  int32_t elem_a, elem_b, inter_only;   /* g(r) pair, 0 = any */
  double rdf_rmax, rdf_dr, qmax, dq, q_direct;
  double fit_from, fit_to;           /* Einstein fit window, fractions of the run (0.2, 0.5) */
  double probe, grid;                /* free volume probe radius (0 = point) and grid (0.4 Å) */
  double cutoff;                     /* cohesive energy cutoff (0: default) */
  int32_t threads;
  int32_t deuterate;                 /* v20 neutron contrast: 0 none, 1 every H, 2 H on aliphatic C (d-backbone), 3 H on aromatic C (d-ring), 4 H on O/N */
  const char* group;                 /* v36: the atoms analysed — NULL or "" all; "selection" the document's selection;
                                        "molecules:1-4,7" molecule ids; "tag:NAME" a tag (v60); "exclude-held" all but the held molecule. The
                                        properties see those atoms alone (density: theirs in the whole cell); the
                                        mechanics and Tg protocols always use the whole cell */
  const char* radii;                 /* v38 free volume and pores: NULL / "bondi" (Bondi 1964), "uff" (x/2), "forcefield" (half
                                        the assigned Lennard-Jones minimum) */
  double zbin;                       /* v44 density profile bin width along the axis (Å; 0: 0.5) */
  int32_t axis;                      /* v44 interfaces and Herman's f: 1 a (x), 2 b (y), 3 c (z); 0: c */
  const char* surface;               /* v44 the surface or filler of zprofile, adhesion, interaction: molecule ids
                                        "1", "1-3,7"; NULL / "": molecule 1 */
} caps_analyze_opts;
typedef int32_t (*caps_analyze_progress_fn)(const char* what, double fraction, void* user);
/* v20: coherent neutron scattering length (fm) of element z (1001 = ²H); NaN when unknown. */
double caps_neutron_b(int32_t z);
int32_t caps_analyze(caps_doc* d, const char* props, const caps_analyze_opts* o, caps_analyze_progress_fn progress, void* user);
/* Mechanics and thermal transitions (ABI 9). caps_analyze_ex takes the ids of caps_analyze plus:
     cij_fluct   elastic constants from stress fluctuations of the saved frames (an NVT run at `temperature`); for
                 bonded systems far too few frames are saved: prefer cij_run
     cij_run     NVT run of a copy of the current frame at `temperature` for `run_ps`, stress sampled every step
     cij_strain  static elastic constants (minimise, ±strain, re-minimise) of the current frame, or of `configurations`
                 frames spread over the selection
     tensile     uniaxial deformation MD of a copy of the current frame (axis, rate 1/ps, max_strain, temperature,
                 lateral_fixed: 0 lateral faces at `pressure`, 1 fixed)
     tg          stepwise cooling of a copy of the current frame (t_start → t_end by t_step K, ps_per_step each)
   The document is never changed; zero fields take the defaults. m may be NULL. */
typedef struct {
  int32_t configurations;
  double strain;
  double temperature;
  int32_t axis;
  double rate, max_strain, fit_strain;
  int32_t lateral_fixed;
  double t_start, t_end, t_step, ps_per_step;
  double dt, pressure;
  uint64_t seed;
  double run_ps;                     /* cij_run: sampled time (100 ps) */
  double equilibrate_ps;             /* tensile and tg: unsampled NPT run first (0: default, < 0: none) */
  /* v46 tg */
  int32_t barostat;                  /* 0 stochastic cell rescaling (default), 1 Berendsen, 2 MTK (Nosé–Hoover chains) */
  double tau_t, tau_p;               /* fs (0: 100 and 1000) */
  double average_from;               /* fraction of each hold discarded before averaging (0: 0.5) */
  int32_t tg_property;               /* 0 specific volume, 1 potential energy per atom */
  int32_t tg_fit;                    /* 0 continuous two-line fit, hinge free; 1 two separate lines through glassy and rubbery ranges */
  double glassy_max, rubbery_min;    /* K, for tg_fit 1 (0: the lowest and highest third of the scan) */
} caps_mech_opts;
int32_t caps_analyze_ex(caps_doc* d, const char* props, const caps_analyze_opts* o, const caps_mech_opts* m, caps_analyze_progress_fn progress,
                        void* user);
/* The last analysis as JSON, as caps analyze --json writes it, plus the frames used:
   {"frames": n, "of": total, "atoms": n, "properties": [...]}. Returns the length needed including the final NUL (call
   with json = NULL to size the buffer); empty before the first analysis. */
int32_t caps_analyze_report(caps_doc* d, char* json, int32_t cap);

/* Save the current frame: .data (LAMMPS full; with GAFF coefficients, angles and dihedrals when every atom can be
   typed), .pdb or .xyz by extension. */
int32_t caps_save(caps_doc* d, const char* path);
/* The LAMMPS input that goes with the LAMMPS data file caps_save writes (same force field: the Field assignment when
   complete, else GAFF of C and H): units, styles, special bonds, read_data <data_name>, neighbour settings, ending
   before any run command. Returns the length needed including the final NUL (text = NULL to size the buffer). */
int32_t caps_lammps_input(caps_doc* d, const char* data_name, char* text, int32_t cap);
/* The LAMMPS fix shake line holding the constraints of mode (1 bonds to hydrogen and rigid water, 2 all bonds) on
   group, with the type numbers of the data file caps_lammps_input reads; "" when nothing is held (ABI 29). */
int32_t caps_lammps_shake(caps_doc* d, int32_t mode, const char* group, char* text, int32_t cap);
/* GROMACS files with the same force field (ABI 22): when stem is non-empty, writes stem.top, stem.gro and stem.mdp (a
   single-point run). Returns the .mdp non-bonded settings (cut-offs, modifiers, dispersion correction, electrostatics)
   preceded by "; note: " lines where GROMACS cannot compute exactly what CAPS does. Fails for forms GROMACS lacks
   (class II, inversions, 9-6 or Buckingham pairs). Returns the length needed including the final NUL. */
int32_t caps_gromacs(caps_doc* d, const char* stem, char* text, int32_t cap);

/* Export center (ABI 23): LAMMPS (STEM.data without pair coefficients, STEM.in with every pair_coeff and the run) and
   GROMACS (STEM.top, STEM.itp, STEM.gro, STEM.mdp; STEM_em.mdp when a run minimises first) written into dir from the
   complete Field assignment (refused otherwise). options {lammps: true, gromacs: true, stem: "system", run: "check" |
   "none" | "minimize" | "nvt" | "npt", minimize_first, temperature (K), pressure (atm), dt (fs; 0 or absent: the force field's own, Martini 20, else 0.5), steps, thermo_every,
   dump_every, seed, lammps_styles: "native" (the force field's own styles: OPLS dihedrals, PPPM …) | "exact" (CAPS's
   energy exactly), hybrid: bool (every style as hybrid), coulomb: "auto" | "pppm" | "ewald" | "dsf" | "cut",
   kspace_accuracy (1e-4), cutoff (Å, 0 the force field's), tail: bool, head_lines: N (each file's first N lines in the
   reply), preview: true (written to a scratch folder, read and removed), dlpoly: true (also STEM_dlpoly/FIELD, CONFIG,
   CONTROL for DL_POLY 4)} → {ok, error, gromacs_error, dlpoly_error, folder, files: [{name, what, bytes, head[]}], notes[], checks: {atoms, typed, types,
   type_pairs, bonds, angles, dihedrals, impropers, missing, net_charge, charges, forcefield, density}}. */
/* v51 the water models (SPC, SPC/E, SPC/Fw, TIP3P, TIP3P CHARMM, TIP3P-Ew, TIP4P, TIP4P-Ew, TIP4P/2005, TIP4P/Ice, OPC) as JSON:
   [{id, name, citation, sites, r_oh, theta, q_h, q_neg, d_om, eps_o, sigma_o, eps_h, sigma_h, rigid, note}]. Apply one with caps_edit
   {"op":"water_model","model":"tip4p2005"}; type the waters with it by a field group {"molecules":"water","water":"tip4p2005"}. */
int32_t caps_water_models(char* json, int32_t cap);
/* v51 Pack input in the CAPS form (cell, distance, molecule … end; pack.hpp) or packmol's: to_caps 1 gives the CAPS form,
   0 packmol's (either given). Returns the length needed with NUL, or -1 (the error names the line). */
int32_t caps_pack_convert(const char* text, int32_t to_caps, char* out, int32_t cap);
/* v51 the packing that made a document: [{name, molecules: "11-110", forcefield: "" | library id | "water:MODEL"}] per input
   molecule, in the input's order ([] for a document not packed). For a force field per group of molecules. */
int32_t caps_pack_items(caps_doc* d, char* json, int32_t cap);
/* v51 the assigned force field written whole (every type, charge, term and setting, 17 digits; caps/ffio.hpp) and read
   back onto the same structure exactly (its atom count checked): a remote job or a saved recipe takes it instead of typing
   again — the same parameters, groups and water models kept. 0, or -1 (caps_last_error). */
int32_t caps_field_save(caps_doc* d, const char* path);
int32_t caps_field_load(caps_doc* d, const char* path);
int32_t caps_export_engines(caps_doc* d, const char* dir, const char* options, char* out, int32_t cap);
int32_t caps_summary_get(caps_doc* d, caps_summary* out);
int32_t caps_set_frame(caps_doc* d, int64_t frame);
/* 1: fold atoms into the cell (bonds across faces are hidden); 0: keep molecules whole (default). */
int32_t caps_set_wrap(caps_doc* d, int32_t wrap);
/* v25 live view of a run: while caps_md or caps_equilibrate runs on d, about four times a second `live` gets the
   positions so far as a new document (wrapped into the cell; the callee closes it) and JSON {step, time_ps, atoms,
   density}. NULL stops it. The callback runs on the run's thread. */
typedef void (*caps_live_fn)(caps_doc* snapshot, const char* stats_json, void* user);
int32_t caps_set_live(caps_doc* d, caps_live_fn live, void* user);
int32_t caps_atom(caps_doc* d, int32_t index, caps_atom_info* out);
int32_t caps_note_count(caps_doc* d);
const char* caps_note(caps_doc* d, int32_t k);

/* Renders into caller memory of width*height*4 bytes, RGBA with straight alpha. */
int32_t caps_render(caps_doc* d, const caps_camera* cam, const caps_render_opts* opt, uint8_t* rgba);
/* v26 GPU view: the scene caps_render draws, in world coordinates, for a view that uploads it once and turns it on the
   GPU (the arrays belong to the document until the next call): spheres x y z r, half-bond and segment capsules
   ax ay az bx by bz r, lines ax ay az bx by bz with a width in pixels; colours 0xRRGGBB with fading and ambient occlusion
   applied; ring 1 selection, 2 keyboard focus. cpu_only: surfaces, polyhedra or the colour-vision preview, which only
   caps_render draws. caps_view_fit: the camera exactly as caps_render fits it for opt's size (× supersample):
   r = R(p − c) + pan, R yaw about y then pitch about x; x = w/2 + r.x·scale·k, y = h/2 − r.y·scale·k,
   k = dist/(dist − r.z) in perspective, else 1; zmin, zmax the shown atoms' view depths (depth cue). */
typedef struct {
  int32_t n_spheres;
  const float* spheres;
  const uint32_t* sphere_rgb;
  const int32_t* sphere_id;
  const uint8_t* sphere_ring;
  int32_t n_capsules;
  const float* capsules;
  const uint32_t* capsule_rgb;
  int32_t n_lines;
  const float* lines;
  const uint32_t* line_rgb;
  const float* line_width;
  int32_t cpu_only;
  uint32_t background;
  int32_t transparent, dark, depth_cue, outlines;
  /* v27: the camera fit's inputs, so a view can fit (caps_view_fit's numbers) without the document, e.g. while a run
     holds it: centre, cell corners (8 × xyz, or 0 when the atoms frame the view), shown atoms' positions, the pad (Å)
     and the perspective field of view. The fit: e = max over corners and points of |R(p − c)| per axis;
     ex, ey = max(e + pad, 2.5); scale = min(0.45 w / ex, 0.45 h / ey) · zoom; dist = ez / tan(fov/2) + ez. */
  double fit_cx, fit_cy, fit_cz, fit_pad, fov_deg;
  int32_t n_fit_corners;
  const float* fit_corners;
  int32_t n_fit_points;
  const float* fit_points;
  /* v54: meshes (surfaces, polyhedra) as triangles: n_triangles; xyz of each corner (9 floats per triangle), the corner
     normals (9 per triangle), corner colours 0xAARRGGBB with AA the opacity. Drawn over the atoms as caps_render does
     (the nearest surface of each pixel blended once), so a view with surfaces stays on the GPU. */
  int32_t n_triangles;
  const float* tri_xyz;
  const float* tri_normal;
  const uint32_t* tri_rgb;
  /* the colour-vision preview (caps_set_vision): the linear-RGB matrix (9, row-major) or NULL for normal vision, and the
     severity it is mixed with the identity by — a GPU view applies it per pixel in linear light, as caps_render does */
  const float* vision_matrix;
  double vision_severity;
} caps_scene;
typedef struct {
  double cos_yaw, sin_yaw, cos_pitch, sin_pitch, cx, cy, cz, scale, w, h, pan_x, pan_y;
  int32_t perspective;
  double dist, zmin, zmax;
} caps_view_fit_t;
int32_t caps_render_scene(caps_doc* d, const caps_render_opts* opt, caps_scene* out);
int32_t caps_view_fit(caps_doc* d, const caps_camera* cam, const caps_render_opts* opt, caps_view_fit_t* out);
int32_t caps_pick(caps_doc* d, int32_t x, int32_t y);
/* v53: the atom under output pixel (x, y) for this camera and these options, by casting the view ray against the drawn
   atoms (the render's styles, radii, hidden and clipped atoms, periodic images, the Visualize pipeline) — nothing is
   rendered, so it answers at once for any size of system and needs no render first. -1: none. */
/* v55: the view's scene for other renderers and 3D tools: format "pov" (POV-Ray 3.7: exact spheres and cylinders, the
   view's camera and key light), "glb" (glTF 2.0 binary, vertex colours: Blender, ParaView …) or "obj" (+ .mtl beside it).
   report: {spheres, cylinders, triangles}. Returns 0, −1 on error. caps_export_image options also take engine: "raytrace"
   with samples, shadows, occlusion, occlusion_strength, aperture (Å) or aperture_fraction (of the view's half-width),
   focus (view depth, Å) and outlines. */
int32_t caps_export_scene(caps_doc* d, const caps_camera* cam, const caps_render_opts* opt, const char* path, const char* format, char* report,
                          int32_t cap);
/* v56: screen x, y (output pixels, 2 floats each) of the listed atoms only — overlays that follow a few atoms without
   projecting the whole structure (caps_project_atoms projects all, with visibility). Returns n, −1 on error. */
int32_t caps_project_indices(caps_doc* d, const caps_camera* cam, const caps_render_opts* opt, const int32_t* atoms, int32_t n, float* xy);
/* The view state of atoms: 0 shown, 1 ghost (faint, never picked), 2 hidden. View only — the structure, its exports and
   calculations keep every atom. atoms NULL: every atom. Returns how many atoms changed. */
int32_t caps_set_atom_state(caps_doc* d, const int32_t* atoms, int32_t n, int32_t state);
/* Lock atoms (locked 1) against picking and edits, or free them (0); atoms NULL: all. Their view state is kept. */
int32_t caps_set_atom_lock(caps_doc* d, const int32_t* atoms, int32_t n, int32_t locked);
/* The layers JSON: {kinds: [{name, formula, elements, atoms, molecules: [{id, atoms, z[8], state, locked, selected, held}]}],
   atoms, z_axis}. Returns the length needed (as the other JSON calls). */
int32_t caps_layers(caps_doc* d, char* json, int32_t cap);
/* The atom count; with out (n ≥ the count) each atom's state: bits 0-1 the view (0 shown, 1 ghost, 2 hidden), bit 2 locked. */
int32_t caps_atom_states(caps_doc* d, uint8_t* out, int32_t n);
/* v60 tags (design/boards/Tags): named, coloured atom sets saved with the structure (NAME.tags.json beside it) and written
   as groups (LAMMPS input "group NAME id …", GROMACS STEM.ndx). caps_tags: {tags: [{name, colour, count, selected, group}]};
   caps_tag_edit JSON {op: set|add|remove|delete|rename|colour, name, colour, new_name, atoms: [0-based …] | "selection"}
   returns the tag's index (-1 with caps_last_error on failure; at most 32 tags); caps_tag_atoms: its atoms, the count returned. */
int32_t caps_tags(caps_doc* d, char* json, int32_t cap);
int32_t caps_tag_edit(caps_doc* d, const char* json);
int32_t caps_tag_atoms(caps_doc* d, const char* name, int32_t* out, int32_t cap);
int32_t caps_pick_at(caps_doc* d, const caps_camera* cam, const caps_render_opts* opt, int32_t x, int32_t y);   /* atom index under pixel of last render, -1 none */
/* v19: pixels per Å at the focal plane for a width × height image of the current frame (exact when orthographic). */
double caps_view_scale(caps_doc* d, const caps_camera* cam, const caps_render_opts* opt);

/* v20: atoms drawn in each level of detail by the last caps_render (near, mid, far) and bond halves drawn. */
int32_t caps_render_stats(caps_doc* d, int64_t* near, int64_t* mid, int64_t* far, int64_t* bonds);
/* v20: memory the document holds: {atoms, frames, topology_bytes, frame_bytes, per_atom_bytes} as JSON. */
int32_t caps_memory(caps_doc* d, char* json, int32_t cap);
int32_t caps_export_png(caps_doc* d, const caps_camera* cam, const caps_render_opts* opt, const char* path);
int32_t caps_export_svg(caps_doc* d, const caps_camera* cam, const caps_render_opts* opt, const char* path);
/* Progress for long series: (done, total, user) -> 0 to continue. */
typedef int32_t (*caps_series_progress_fn)(int32_t done, int32_t total, void* user);
/* v20 provenance (design/boards/Provenance): the steps that produced the document — {schema "caps-manifest/1.0",
   generator, deterministic, inputs [{name, sha256}], steps [{engine, summary, params{}, rng, cites[], approximations{},
   time}], approximations{}}. caps_save writes it beside the file as <file>.provenance.json and opening a file reads it
   back. caps_provenance_file reads a file's (ok false without one); caps_provenance_compare gives {same_inputs,
   same_generator, steps_a, steps_b, differing_steps, rows [{step, engine, key, a, b}], notes[]}; caps_provenance_bibtex
   the BibTeX of every method cited (-1 on error). */
int32_t caps_provenance(caps_doc* d, char* json, int32_t cap);
/* A decision taken on the structure recorded as a provenance step (ABI 29): JSON {engine: "equilibrate.accepted" …,
   summary, params: {key: value …}}. 0, or -1 with caps_last_error. */
int32_t caps_provenance_note(caps_doc* d, const char* json);
int32_t caps_provenance_file(const char* path, char* json, int32_t cap);
int32_t caps_provenance_compare(const char* a_json, const char* b_json, char* json, int32_t cap);
int32_t caps_provenance_bibtex(const char* manifest_json, char* text, int32_t cap);
/* A readable reference ("Authors, \"Title\", Journal Volume, Pages (Year). doi:…") for a citation key; "" if unknown. */
int32_t caps_citation_text(const char* key, char* text, int32_t cap);
/* A methods paragraph for a paper from a manifest (replicas_json: an array of manifests that differ only in their seeds,
   or NULL): {ok, text, refs[]} with numbered references in order of first use. */
int32_t caps_methods_text(const char* manifest_json, const char* replicas_json, char* json, int32_t cap);
/* v20 voids (design/boards/FreeVolume) of the current frame: options {"grid":0.5, "probe":1.4, "count":40, "min_radius":1.0,
   "show":true (drawn as translucent spheres), "clear":true (remove them)} → {ok, accessible_point, accessible_probe,
   largest, grid[3], spheres [{x, y, z, r}]} — the largest non-overlapping empty spheres (radius to the nearest Bondi
   surface), biggest first. caps_voids_pdb writes them as HETATM VOI records with the radius as the B-factor. */
int32_t caps_voids(caps_doc* d, const char* options_json, char* json, int32_t cap);
int32_t caps_voids_pdb(caps_doc* d, const char* path);
/* v20 pores (design/boards/SlitPore): {"kind":"slit"|"cylinder"|"framework", "width":10 (slit H, C–C centres, or channel
   diameter), "layers":1, "lx":26, "ly":22, "vacuum":false, "vacuum_gap":20, "cif":"<crystal for cylinder/framework>",
   "wall":6, "length":20, "repeat":[2,2,2], "passivate":false, "fluid":"<SMILES>", "fluid_name", "count", "tolerance":2,
   "seed":1}. The walls are molecule 1; the fluid is packed inside the pore. Report JSON {wall_atoms, fluid_molecules,
   width, pore_volume, fluid_density, dmin, notes[]}. NULL on error. */
caps_doc* caps_pore_build(const char* options_json, char* report, int32_t cap);
/* v20 Kremer–Grest melts (design/boards/CoarseGrained): {"chains":50, "beads":100, "density":0.85, "k_theta":0, "seed":1};
   report {box, closest, r2_per_bond}. caps_kg_lammps writes STEM.data and STEM.in (push-off, then FENE + WCA; steps ≤ 0:
   20 000 and 100 000). */
caps_doc* caps_kg_build(const char* options_json, char* report, int32_t cap);
/* v39 MARTINI polymer melt: {forcefield: path (MARTINI), repeat: bead SMILES of the repeat unit ([SN0], [C1] …), repeats,
   chains, density (g/cm³), tolerance (Å, 3), seed}. One chain built from the repeat, `chains` copies packed at a quarter
   of the density, the MARTINI force field assigned, the cell compressed to the density with it; the melt comes back
   assigned (ready for Relax, Dynamics, GROMACS and LAMMPS). progress gets (stage, stages, 0). Report: {beads,
   beads_per_chain, chain_mass, density, loose_density, forcefield, complete}. A Kremer–Grest melt (caps_kg_build) also
   comes back with its own force field (FENE + WCA; v39). */
caps_doc* caps_martini_melt(const char* options_json, caps_progress_fn progress, void* user, char* report, int32_t cap);
int32_t caps_kg_lammps(caps_doc* d, const char* options_json, const char* stem, double pushoff_steps, double run_steps);
/* v20 colour vision (design/boards/ColourVision): palettes {"Elements": {"labels": [...], "colours": ["#909090", ...]}, ...}
   → {"palettes": [{name, labels, normal[], protanopia[], deuteranopia[], tritanopia[]}], "pairs": [{palette, vision,
   a, b, de}]} with every pair closer than threshold ΔE*ab (12) under a deficiency, closest first (Machado et al. 2009,
   severity 1). caps_set_vision previews the view as seen with one (0 normal, 1 protanopia, 2 deuteranopia,
   3 tritanopia); exports are never simulated. */
int32_t caps_vision_check(const char* palettes_json, double threshold, char* json, int32_t cap);
int32_t caps_set_vision(caps_doc* d, int32_t vision, double severity);
/* v20 motion (design/boards/Motion): the camera that frames atoms idx[0..n) — same yaw and pitch as cam, the pan that
   centres them and the zoom at which their extent fills `fill` (0.6) of the view (1 … 40); n = 0 frames everything. */
int32_t caps_camera_focus(caps_doc* d, const caps_camera* cam, const int32_t* idx, int32_t n, double fill, caps_camera* out);
/* v20 periodic box (design/boards/PeriodicBox). caps_periodic JSON {"molecule": k (1-based; 0: the first that crosses a
   face)} → {ok, molecules, crossing (molecules that cross a face when whole), pieces (fragments when every atom is
   wrapped into the cell), box: [a, b, c], cubic, molecule, bond: {i, j, wrapped, min_image, whole} (a bond of it that
   crosses a face), ends: {i, j, wrapped, min_image, whole} (its first and last atom)}; distances in Å.
   caps_set_images draws na × nb × nc copies around the cell, faded (1 1 1: off); caps_set_save_wrap(d, 0 as shown,
   1 atoms into the cell, 2 molecule centres into the cell) applies to caps_save; caps_centre_on moves every atom so the
   centre of atoms idx[0..n) sits at the cell centre (undoable). */
int32_t caps_periodic(caps_doc* d, const char* json, char* out, int32_t cap);
int32_t caps_set_images(caps_doc* d, int32_t na, int32_t nb, int32_t nc, double fade);
int32_t caps_set_save_wrap(caps_doc* d, int32_t mode);
int32_t caps_centre_on(caps_doc* d, const int32_t* idx, int32_t n);
/* v20 one molecule (design/boards/MoleculeInspector), the one containing `atom`: {ok, molecule, atoms, bonds, rings,
   formula, smiles (written from the perceived bonds, not canonical), mass, monoisotopic, dbe, inertia [3] (amu·Å²,
   ascending), inertia_defect, rg, has_charges, net_charge, dipole (D; null without charges), rotatable}. */
int32_t caps_molecule_info(caps_doc* d, int32_t atom, char* out, int32_t cap);
/* v20 solvent-accessible surface area (design/boards/SurfaceArea), Shrake–Rupley with Bondi radii: {"probe": 1.4,
   "points": 200, "convergence": false, "colour": false, "frames": {first, last, every}} → {ok, total, area[] (Å² per atom),
   groups [{name, area, share}], convergence [{points, total, delta}]; with frames (a trajectory): series [{frame, timestep,
   total}], mean, sd, frames_used}; colour: true colours the view by exposure (area / full sphere).
   caps_set_atom_values colours the view by any per-atom quantity on a ramp (0 viridis, 1 blue–orange, 2 red–white–blue;
   n = 0 clears). */
int32_t caps_sasa(caps_doc* d, const char* json, char* out, int32_t cap);
int32_t caps_set_atom_values(caps_doc* d, const double* values, int32_t n, int32_t ramp);
/* Two states of the same atoms compared (ABI 29, design/boards/Compare): JSON {reference: {kind: "start" | "current" |
   "frame" | "snapshot", index}, moving: {…} (default current), fit: "all" | "heavy" | "backbone" | "selection" | "none",
   periodic: "yes" | "no" (each atom's nearest image to its reference), largest: 20, per_atom: 1 (the shifts array),
   colour: 1 (the view shows each atom's shift on a ramp)} or {op: "clear"} (colouring off) → {ok, error, rmsd: {all,
   heavy, backbone}, fit, fitted, atoms, reference, moving, max, mean, largest: [{atom, label, element, molecule,
   backbone, shift}], shifts}. The moving state is superposed on the reference (Horn 1987 quaternions). */
int32_t caps_compare_states(caps_doc* d, const char* json, char* out, int32_t cap);
/* One state of the document ({kind: current | start | frame | snapshot, index}) as a new document (v34), for the split view. */
caps_doc* caps_state_document(caps_doc* d, const char* json);
/* The Properties explorer (ABI 29). caps_structure_info: {ok, formula (Hill), composition {symbol: count}, title, format,
   atoms, bonds, molecules, frames, frame, timestep, mass (g/mol), charge, has_charges, unwrapped, bonds_from_file,
   cell {a, b, c, alpha, beta, gamma, volume, density, periodic [3]}, forcefield {name, complete, missing, types},
   selected}. caps_atom_properties: {ok, index, id, element, name, type, ff_type, charge, mass, molecule, residue,
   resname, xyz [3], fractional [3], neighbours [{index, element, distance, order}]}. */
int32_t caps_structure_info(caps_doc* d, char* out, int32_t cap);
int32_t caps_atom_properties(caps_doc* d, int32_t index, char* out, int32_t cap);
/* v20 cell editor (design/boards/CellEditor): caps_set_cell JSON {a, b, c, alpha, beta, gamma, scale: true} (a along x,
   b in xy; scale keeps fractional coordinates, false leaves the atoms where they are); caps_supercell replicates the
   frame na × nb × nc (atoms, bonds, molecules). Both undoable. */
int32_t caps_set_cell(caps_doc* d, const char* json);
int32_t caps_supercell(caps_doc* d, int32_t na, int32_t nb, int32_t nc);
/* v20 partial charges (design/boards/Charges): {"method": "gasteiger" | "qeq" | "forcefield" (the Field assignment) |
   "file" {path: .chg} | "keep", "apply": false} → {ok, error, method, q[], net, max_abs, groups: [{name, n, mean, lo, hi}],
   edges[], counts[], notes[]}. apply: true sets the charges on the structure (undoable, clears the Field assignment). */
int32_t caps_charges(caps_doc* d, const char* json, char* out, int32_t cap);
/* v20 the k-th molecule / chain colour of the current palette (0xRRGGBB). */
uint32_t caps_category_colour(int32_t k);
/* v20 recipes (design/boards/CommandLine, JupyterNotebook): runs a recipe (YAML or JSON text, see recipe.hpp) — build,
   type, grow, relax, md, equilibrate, analyze, export — and returns the structure with its force field and provenance.
   options: {"base_dir", "out_dir", "forcefield_dir", "seed": -1 (≥ 0 overrides), "threads": 0}. progress (may be NULL)
   gets every stage event; return 0 to cancel. report: {exit, error, files[], properties[], forcefield}. NULL on error
   (exit 2 input, 3 missing parameters, 4 failed run). */
typedef int32_t (*caps_recipe_progress_fn)(int32_t stage, int32_t stages, const char* name, const char* status, const char* detail, double fraction, void* user);
caps_doc* caps_recipe_run(const char* recipe, const char* options_json, caps_recipe_progress_fn progress, void* user, char* report, int32_t cap);
/* v20 a recipe checked without running it (design/boards/RecipeEditor): {ok, code, error, name, sha256, stages: [{name,
   summary, ok}], protocol, schedule: [{label, ensemble, ps, t_start, t_end, pressure_bar}]}. Runs record the recipe's
   SHA-256 as the first step of their provenance. */
int32_t caps_recipe_check(const char* recipe, char* out, int32_t cap);
/* v20 the YAML subset recipes use, as JSON ({"ok": false, "error"} when it does not parse). */
int32_t caps_yaml_to_json(const char* yaml, char* out, int32_t cap);
/* v20 the current frame for a viewer outside the Studio (the notebook's caps.View): {"atoms", "shown", "z": [...],
   "xyz": [x0, y0, z0, …] (Å, 3 decimals), "bonds": [i0, j0, …], "colours": {"6": "#909090", …}, "radii": {"6": 1.7, …}
   (van der Waals), "cell": [ox, oy, oz, ax, ay, az, bx, …] or null}. options: {"max_atoms": 60000, "hydrogens": true}
   (beyond max_atoms every k-th molecule is kept). */
int32_t caps_scene_json(caps_doc* d, const char* options_json, char* json, int32_t cap);
/* v20 reaction template editor (design/boards/ReactionTemplate): caps_template_view parses template text and returns
   {ok, error, templates: [view]} (see template_view in react.hpp); caps_template_test counts, on the document's current
   frame, each template's reactive sites and the matches within its capture distance: {ok, templates: [{name, sites,
   matches, closest}]}. */
int32_t caps_template_view(const char* text, char* json, int32_t cap);
int32_t caps_template_test(caps_doc* d, const char* text, char* json, int32_t cap);
/* v20 export dialog (design/boards/ExportDialog). options: {"bits":8|16, "dpi":600, "colour_profile":"srgb"|"none",
   "provenance":true, "source":"<the structure's path, hashed into the manifest>"}. 16-bit keeps the supersampled
   average at full precision. The manifest (caps-image/1.0: generator, created, source + sha256, frame, atoms, camera,
   render) goes into an iTXt chunk "caps:provenance"; caps_png_text reads a PNG's text chunks back as a JSON object. */
int32_t caps_export_image(caps_doc* d, const caps_camera* cam, const caps_render_opts* opt, const char* path, const char* options_json,
                          const uint8_t* overlay);   /* overlay: width*height*4 straight-alpha RGBA laid over the image, or NULL */
int32_t caps_png_text(const char* path, char* json, int32_t cap);
/* Movie: {"format":"apng"|"png_sequence", "from", "to", "step", "fps":24, "loops":0, "turntable":<degrees of yaw over the
   movie, 0 none: the current frame turned, "turntable_frames":120 of it>, "colour_profile", "provenance", "source"}.
   png_sequence writes path/frame_00001.png …; returns the frames written (progress may stop it early), -1 on error.
   The document's current frame is restored. */
int32_t caps_export_movie(caps_doc* d, const caps_camera* cam, const caps_render_opts* opt, const char* path, const char* options_json,
                          caps_series_progress_fn progress, void* user);

/* Geometry of 2–4 atoms with minimum-image separations: distance (Å), angle (deg) or dihedral (deg). */
int32_t caps_measure(caps_doc* d, const int32_t* idx, int32_t n, double* value);

/* g(r) between elements (0 = any), optionally only between different molecules. Returns bins written. */
int32_t caps_rdf(caps_doc* d, int32_t elem_a, int32_t elem_b, double rmax, double dr, int32_t inter_only, double* r, double* g, int32_t cap);

typedef struct {
  int32_t molecule, atoms;
  double mass, rg, kappa2, com_x, com_y, com_z;
} caps_molecule;

/* Per-molecule size and shape (gyration tensor, whole molecules). Returns count written. */
int32_t caps_molecules(caps_doc* d, caps_molecule* out, int32_t cap);

/* Range of the per-atom property used by colour_by = 3 (distance to own molecule centre, Å). */
int32_t caps_property_range(caps_doc* d, double* lo, double* hi);

/* Nearest neighbours of an atom (minimum image): fills up to k indices and distances, returns count. */
int32_t caps_neighbours(caps_doc* d, int32_t index, int32_t k, int32_t* idx, double* dist);
/* v19 visualize pipeline (caps/pipeline.hpp): steps as JSON ({"steps":[{"type":…,"enabled":…,…}]}), run on every
   shown frame; the view then draws its particles and colours (selected in red), and picks return the frame's atom.
   NULL or "" clears it. Results, particle rows (filter: an expression) and bond rows are JSON, sized like the others. */
int32_t caps_pipeline_set(caps_doc* d, const char* json);
int32_t caps_pipeline_result(caps_doc* d, char* json, int32_t cap);
int32_t caps_pipeline_particles(caps_doc* d, const char* filter, int32_t offset, int32_t count, char* json, int32_t cap);
int32_t caps_pipeline_bonds(caps_doc* d, int32_t offset, int32_t count, char* json, int32_t cap);
int32_t caps_pipeline_catalogue(char* json, int32_t cap);
/* v19: where the caps Python package for Python steps is (data/python), and the interpreter (NULL keeps python3). */
void caps_set_python(const char* package_dir, const char* interpreter);
/* v19 pipeline YAML (design/boards/SavePipeline): steps JSON to YAML text (name, source file and topology optional),
   and YAML back to steps JSON ({"steps":[…],"name":…,"file":…,"topology":…}). */
int32_t caps_pipeline_to_yaml(const char* json, const char* name, const char* file, const char* topology, char* yaml, int32_t cap);
int32_t caps_pipeline_from_yaml(const char* yaml, char* json, int32_t cap);
/* v19 export (design/boards/ExportData): format lammps-data | lammps-dump (all frames) | dcd (all frames) | gro | pdb | xyz |
   mol2 | sdf | cif | car;
   options JSON {"pipeline": bool (the Visualize result instead of the frame), "wrap": bool, "coeffs": bool (force-field
   sections in LAMMPS data)}. The preview writes to a scratch file and returns {"lines":[first n],"bytes","atoms","bonds",
   "angles","dihedrals","atom_types","bond_types","angle_types","dihedral_types","notes":[…]}. */
int32_t caps_export_data(caps_doc* d, const char* path, const char* format, const char* options);
/* v30 the Visualize pipeline's last grid (the density field step's) written by extension: .cube (Gaussian cube, with the
   atoms), .vtk (legacy structured grid), .npy (NumPy float64 n0 × n1 × n2). Returns 0, or −1 (no grid; caps_last_error). */
int32_t caps_pipeline_export_grid(caps_doc* d, const char* path);
/* The outputs of a pipeline (v35, design/boards/SavePipeline "outputs:"): pipeline JSON {steps: [...], outputs: [{kind:
   table | plot | attributes | render | grid, what, path, width, height}]} run on the shown frame and its outputs written
   under dir; the report JSON {ok, lines: ["wrote …" | "not written: …"]}. */
int32_t caps_pipeline_write_outputs(caps_doc* d, const char* pipeline_json, const char* dir, char* report, int32_t cap);
int32_t caps_export_preview(caps_doc* d, const char* format, const char* options, int32_t lines, char* json, int32_t cap);
/* v19 figure bundles (design/boards/FigureBundle): options JSON {"name","input","topology","frame","width","height",
   "include_input","include_pipeline","include_data","include_readme","pipeline":{"steps":[…]}}. The preview lists
   the files (name, note, bytes, sha256; figures only on write) with provenance.json and the first data file. */
int32_t caps_bundle_preview(caps_doc* d, const char* options, char* json, int32_t cap);
int32_t caps_bundle_write(caps_doc* d, const char* path, const char* options, const caps_camera* cam, const caps_render_opts* opt);
/* The pipeline's global attributes on every stride-th frame: {"columns":[…],"rows":[[…]]}. progress may be NULL. */
int32_t caps_pipeline_series(caps_doc* d, int32_t stride, caps_analyze_progress_fn progress, void* user, char* json, int32_t cap);
/* Atoms bonded to atom index (up to cap written); returns the count. */
int32_t caps_bonded(caps_doc* d, int32_t index, int32_t* idx, int32_t cap);
/* Each atom's molecule number as CAPS counts molecules everywhere (Field groups, the held molecule, the LAMMPS molecule
   column): the file's own when it has them, else the bonded fragments in order from 1 (a PDB's chain ids are not
   molecules). Returns the atom count (out = NULL to size). v37 */
int32_t caps_molecule_ids(caps_doc* d, int64_t* out, int32_t cap);
/* The molecule (connected by bonds, 0-based) of every atom, up to cap written; returns the number of molecules. */
int32_t caps_molecule_index(caps_doc* d, int32_t* mol, int32_t cap);

/* Molecule builder (ABI 12). caps_smiles_info parses a SMILES and returns JSON without building: {ok, error, position,
   formula, mass, atoms, heavy, bonds, rings, stereocentres, stereo_bonds, charge, problems[]}; the length needed
   including the final NUL is returned (json = NULL to size). caps_build_smiles embeds `conformers` conformers in 3D and,
   when ff_path names a caps-forcefield JSON with typing rules, minimises each with it; the document holds one frame
   per conformer, lowest energy first. The report is JSON: the fields of caps_smiles_info plus method, conformers
   [{energy, rel, minimised}] and notes[]. NULL on error (caps_last_error explains). */
typedef struct {
  int32_t conformers;              /* 1 */
  uint64_t seed;                   /* 1 */
  int32_t rotor_search;            /* v40: 1 each minimised conformer's rotatable bonds tried at their staggered positions */
  int32_t heavy_only;              /* v40: 1 no implicit hydrogens (the SMILES as written: united-atom models) */
} caps_build_opts;
int32_t caps_smiles_info(const char* smiles, char* json, int32_t cap);
/* The 2D drawing of a SMILES, as JSON: the caps_smiles_info fields plus atoms [{z, symbol, x, y, h (hydrogens drawn as a
   label), charge, isotope, hcount (-1 implicit), aromatic, bracket, chiral, map, order[]}] and bonds [{a, b, order, dir}];
   bond length 1. Returns the length needed including the final NUL (json = NULL to size). */
int32_t caps_smiles_depict(const char* smiles, char* json, int32_t cap);
/* SMILES from a graph in the JSON form caps_smiles_depict gives (x, y, h and symbol are not needed). Returns the length
   needed including the final NUL, or -1 when the graph is not valid (caps_last_error explains). */
int32_t caps_smiles_write(const char* graph_json, char* smiles, int32_t cap);
caps_doc* caps_build_smiles(const char* smiles, const char* ff_path, const caps_build_opts* o, char* report, int32_t cap);

/* Coarse-grained molecules (v24). caps_build_beads: `text` is a bead template of the force field at ff_path (its
   "bead_templates", e.g. MARTINI's DPPC) or bead SMILES ([TYPE] and [TYPE±q] beads, branches, ring closures); one site per
   bead named by its type (element 0), bond lengths and masses from the force field when ff_path is given; a start for a
   relax. Report (JSON): {beads, bonds, charge, template, forcefield}. caps_bead_templates writes {"NAME": "bead SMILES",
   …} for ff_path; returns the size needed including the final NUL (json = NULL to size), or −1 on error. */
caps_doc* caps_build_beads(const char* text, const char* ff_path, uint64_t seed, char* report, int32_t cap);
int32_t caps_bead_templates(const char* ff_path, char* json, int32_t cap);

/* Settings (ABI 13), process-wide: the colour palette of elements and molecules/chains (0 CAPS, 1 Okabe–Ito, 2
   monochrome) and the worker threads of the parallel loops (0 = one per hardware thread, at most 16). */
void caps_set_palette(int32_t palette);
void caps_set_threads(int32_t threads);
/* Electrostatics for every later calculation (ABI 16): mode 0 damped shifted force, 1 particle-mesh Ewald (periodic
   cells; others keep DSF); ewald_rtol sets β by erfc(β rc) = rtol, pme_spacing the grid (Å), pme_order the B-splines. */
void caps_set_electrostatics(int32_t mode, double ewald_rtol, double pme_spacing, int32_t pme_order);

/* Bench (ABI 14): the built-in validation suite. caps_bench_list gives the tables as JSON [{id, title, scope, columns[],
   status: "not run", note}]. caps_bench_run runs one table (samples: the directory with ps_melt.data and water.pdb;
   forcefields: data/forcefields, or NULL) and returns it as JSON {id, title, scope, columns[], rows[{cells[], status}],
   status, note, seconds}. caps_bench_write writes results.md, results.tex and one CSV per table from a JSON array of
   tables into `dir`. Sizes as caps_field_report: the length needed including the final NUL, or -1 on error. T6 / T7 measure the user's equilibrated cells (CAPS_BENCH_CELLS, else ~/.caps/bench/cells; named
   by reference id, e.g. ps-atactic.data, with a trajectory of the same name for frames) against the ranges in
   reference/polymers.json beside forcefields; without cells they are not run and say so. */
typedef int32_t (*caps_bench_progress_fn)(const char* table, const char* what, double fraction, void* user);
int32_t caps_bench_list(char* json, int32_t cap);
int32_t caps_bench_run(const char* id, const char* samples, const char* forcefields, int32_t repeats, int32_t quick,
                       caps_bench_progress_fn progress, void* user, char* json, int32_t cap);
int32_t caps_bench_write(const char* tables_json, const char* dir);

/* Polymer builder (ABI 15). A chain spec is JSON: {units: [{name, smiles}], sequence: "homopolymer" | "alternating" |
   "block" | "random" | "gradient" | "pattern", dp, blocks: [..], weights: [..], pattern: "AAB", pm, forcefield}; the
   repeat-unit SMILES carry two attachment points, head first. caps_unit_info: {ok, error, formula, mass, atoms,
   head_element, tail_element, stereocentres}. caps_chain_preview: {ok, error, sequence: [unit index …], formula, mass,
   atoms, smiles} of one chain drawn with `seed`. caps_grow_chains grows o->chains chains of the spec (o->dp, when > 0,
   overrides the spec's; o->tacticity applies) into a periodic cell; NULL on error or cancel. */
int32_t caps_unit_info(const char* smiles, char* json, int32_t cap);
int32_t caps_chain_preview(const char* spec_json, uint64_t seed, char* json, int32_t cap);
caps_doc* caps_grow_chains(const char* spec_json, const caps_grow_opts* o, caps_progress_fn progress, void* user, char* report, int32_t cap);
/* v20 the same with a live view (design/boards/GrowAllAtom): about four times a second `live` gets the chains so far as a
   new document (atoms and bonds, molecule = chain, no end caps; the callee closes it with caps_close) and
   {chains_done, chains, units, units_total, restarts, worst_margin, density, atoms}. Called on the growing thread. */
typedef void (*caps_grow_live_fn)(caps_doc* snapshot, const char* stats_json, void* user);
caps_doc* caps_grow_chains_live(const char* spec_json, const caps_grow_opts* o, caps_progress_fn progress, caps_grow_live_fn live, void* user, char* report, int32_t cap);

/* Surfaces and interfaces (v17). caps_surface_terminations: {ok, error, d, formula, atoms, density, cell: [a, b, c,
   alpha, beta, gamma], notes, terminations: [{label, top, bottom, gap, bonds_per_nm2}]} for (hkl) of a CIF file, fewest
   bonds cut first. caps_surface_build: a slab document; options JSON {h, k, l, layers, termination (0-based), vacuum,
   orthogonal, max_strain, na, nb, passivate}. caps_interface_build: a slab with a polymer film grown on it; options JSON
   {crystal: CIF path, slab: {as caps_surface_build}, film: {thickness, density, chains, gap, vacuum}}, the chain spec as
   caps_grow_chains, grow options (seed, contact_scale; chains and density come from film). The slab is molecule 1. */
int32_t caps_surface_terminations(const char* cif_path, int32_t h, int32_t k, int32_t l, char* json, int32_t cap);
caps_doc* caps_surface_build(const char* cif_path, const char* options_json, char* report, int32_t cap);
caps_doc* caps_interface_build(const char* options_json, const char* spec_json, const caps_grow_opts* o, caps_progress_fn progress, void* user, char* report,
                               int32_t cap);
/* Nanostructures (v17). Options JSON {kind: "tube" | "sheet" | "particle"; tube: n, m, length, periodic; sheet: lx, ly,
   layers, periodic; particle: crystal (CIF path), shape (sphere | cube | octahedron | cuboctahedron), radius, on_atom,
   passivate}. caps_nano_build: the filler alone. caps_nano_embed: the filler held at the centre of a periodic cell with
   polymer chains grown around it; options add matrix: {chains, density}; the chain spec as caps_grow_chains; the filler
   is molecule 1 and is held in Relax / Dynamics. */
caps_doc* caps_nano_build(const char* options_json, char* report, int32_t cap);
/* The same around the document's structure (v34: a filler built and functionalised in the Studio): options {matrix:
   {chains, density}}; the directions the structure spans across its cell stay periodic; its history comes along. */
caps_doc* caps_embed_document(caps_doc* filler, const char* options_json, const char* spec_json, const caps_grow_opts* o, caps_progress_fn progress, void* user,
                              char* report, int32_t cap);
caps_doc* caps_nano_embed(const char* options_json, const char* spec_json, const caps_grow_opts* o, caps_progress_fn progress, void* user, char* report,
                          int32_t cap);
/* Layer stacks (v35, layers.hpp): the documents' shown frames piled along z, the first at the bottom; options {names:
   [...], gap (Å), vacuum (Å, 0: periodic in z), match: "both" | "first", max_repeat}. Every cell rectangular. The report
   is JSON {ok, a, b, c, layers: [{name, na, nb, strain_a, strain_b, z_lo, z_hi, atoms}], notes: [...]}; the first
   layer is molecule 1 … and is held when it is one molecule (a slab). */
caps_doc* caps_stack_documents(caps_doc* const* docs, int32_t n, const char* options_json, char* report, int32_t cap);
/* A new document holding the shown frame of d with its history (v35): a layer or a copy that outlives d. */
caps_doc* caps_frame_copy(caps_doc* d);
/* Polymer blends (v17): options JSON {components: [{spec: {as caps_grow_chains}, weight, chains}], chains (of the first
   component), density (growth), morphology: "mixed" | "slabs"}; grow options for seed and contact scale. */
caps_doc* caps_grow_blend(const char* options_json, const caps_grow_opts* o, caps_progress_fn progress, void* user, char* report, int32_t cap);
/* Inserts `count` copies of a molecule (SMILES; hydrogens added, cleaned with UFF) into the free space of the current
   frame (the structure stays where it is), e.g. H–S–S–H sulfur donors for the sulfur_allylic cure. The document becomes
   that one frame; a Field assignment is cleared. 0 on success. */
int32_t caps_insert_molecules(caps_doc* d, const char* smiles, int32_t count, double tolerance, uint64_t seed, char* report, int32_t cap);
/* File checks of the open document (v17): JSON [{level: pass|note|warn|error, title, detail, action}]; the action is
   one suggested fix id (wrap, relax, field) or empty. Returns the length needed including the final NUL. */
int32_t caps_file_checks(caps_doc* d, char* json, int32_t cap);
/* Holds every atom of molecule `mol` in place in caps_relax, caps_md and caps_equilibrate (0: none), e.g. the substrate of
   an interface. A document from caps_interface_build holds molecule 1. */
void caps_set_held_molecule(caps_doc* d, int64_t mol);
/* Distance restraints for caps_relax, as JSON: [{"i": 0, "j": 5, "r0": 3.0, "k": 10}, …] (atom indices from 0, Å,
   kcal/mol/Å²); entries with "l" are dihedral restraints (ABI 29): {"i", "j", "k", "l", "phi0": degrees, "kphi":
   kcal/mol/rad²}; "[]" or NULL clears them. Returns the number set, or -1 on a malformed list. */
int32_t caps_set_restraints(caps_doc* d, const char* json);
/* Add hydrogens by pH (v21): amino-acid residues protonated at `ph` (model pKa values) before the "add_h" edit and
   in caps_hydrogen_plan; a negative pH goes back to neutral valences. "add_h" also takes {ph} for one edit. */
void caps_set_ph(caps_doc* d, double ph);

/* χ by MD (v21, chimd.hpp): JSON {polymer: {units: [{name, smiles}], dp …} (a chain spec), chains, solvent (SMILES) and
   solvent_molecules, or polymer_b and chains_b, temperature, pressure, eq_ps, prod_ps, seed} → {ok, error, chi,
   chi_error, phi_a, v_ref, de_mix, cells: [{name, atoms, molecules, density, ced, ced_error}], notes}. Minutes of
   CPU: progress gets (stage, fraction) and returns non-zero to cancel. */
int32_t caps_chi_md(const char* json, caps_stage_fn progress, void* user, char* out, int32_t cap);
/* χ(T) from pair contacts (v22; Fan, Olafson, Blanco & Hsu 1992): json {a, b (SMILES; * ends capped with H),
   forcefield (a caps-forcefield JSON with typing rules; empty: built-in GAFF subset or UFF), samples, pack_trials,
   temperatures [K], t (report temperature), seed}. Returns {ok, chi, chi_error (at t), fit_a, fit_b (χ = A + B/T),
   temperatures[], chi_t[], chi_t_error[], kinds: [{name, z, z_error, e_min, e_mean, e_t[], hist_e[], hist_p[]}],
   forcefield, a, b, notes}. Seconds for small molecules. */
int32_t caps_chi_contacts(const char* json, caps_stage_fn progress, void* user, char* out, int32_t cap);
/* RIS reference (v21): C_n of polyethylene (Flory's three-state model) at temperature T for n = 1 … nmax into out;
   returns nmax. */
int32_t caps_ris_cn(double temperature, int32_t nmax, double* out);
/* Backmap (v21): the current structure's coarse-grained beads (per_bead backbone atoms each, as caps_resolution_convert
   makes them) moved to those of the file (the last frame; same count and order, e.g. a LAMMPS dump of the beads); each
   bead's atoms carried along and turned with it, then relaxed when relax is set. A new document, or NULL. */
caps_doc* caps_backmap(caps_doc* d, const char* beads_path, int32_t per_bead, int32_t relax, char* report, int32_t cap);
int64_t caps_held_molecule(const caps_doc* d);
/* Atoms held in place besides the held molecule (v36; indices from 0 in the current frame): no force or motion in Relax,
   Dynamics and Equilibrate, a freeze group in the GROMACS files. n = 0 clears them. Returns the number set, -1 on error.
   caps_fixed_atoms copies up to cap indices and returns how many there are. */
int32_t caps_set_fixed_atoms(caps_doc* d, const int32_t* atoms, int32_t n);
int32_t caps_fixed_atoms(const caps_doc* d, int32_t* atoms, int32_t cap);

/* Crystals from space groups (v20, design/boards/CrystalBuilder). A spec is JSON {space_group (key "227:2", number or
   Hermann–Mauguin symbol), a, b, c (Å), alpha, beta, gamma (°), sites: [{label, element (symbol), x, y, z (fractional),
   occupancy}], supercell: [na, nb, nc], tolerance (Å: images closer are one atom), primitive (bool: the primitive cell of
   a centred lattice), title}.
   caps_space_groups: [{key, number, hm, hall, system}] of the 530 settings of the 230 groups.
   caps_crystal_info: {ok, error, key, hm, hall, system, number, operations, atoms_per_cell, atoms, bonds,
   multiplicity: [per site], volume (Å³ of the cell), density, formula, centring, notes}.
   caps_crystal_build: the crystal as a document; report lines in `report`.
   caps_crystal_symmetrize: {ok, error, moved, spec} with every site moved onto its special position (images within
   `snap` Å of the site averaged).
   caps_crystal_find_symmetry: the highest-symmetry setting that maps the crystal of a spec (cif_path NULL or empty) or of
   a CIF file onto itself within `tolerance` Å: {ok, error, key, hm, system, number, operations, atoms, spec (that
   setting's asymmetric unit and the cell)}. All return the length needed including the final NUL. */
int32_t caps_space_groups(char* json, int32_t cap);
int32_t caps_crystal_info(const char* spec_json, char* json, int32_t cap);
caps_doc* caps_crystal_build(const char* spec_json, char* report, int32_t cap);
int32_t caps_crystal_symmetrize(const char* spec_json, double snap, char* json, int32_t cap);
int32_t caps_crystal_find_symmetry(const char* spec_json, const char* cif_path, double tolerance, char* json, int32_t cap);

/* Peptides (v20, design/boards/BioBuilder). Options JSON {sequence (one-letter codes), structure (per residue H α-helix,
   E β-strand, P PPII, C coil), helix / strand / ppii: [φ, ψ, ω] in degrees, n_term: "NH3+" | "NH2" | "ACE", c_term:
   "COO-" | "COOH" | "NME", ph, neutral (bool), cleanup (bool: UFF), seed, ribbon (bool: the document draws a tube
   through the CA atoms coloured by structure until a pipeline replaces it)}.
   caps_peptide_info: {ok, error, residues, atoms, charge, formula, mass, smiles, structure} without building the clean-up.
   caps_peptide_build: the peptide as a document (atoms named as in PDB files, residue names and numbers).
   caps_fasta_sequence: the sequence of the first record of a FASTA text. */
int32_t caps_peptide_info(const char* options_json, char* json, int32_t cap);
caps_doc* caps_peptide_build(const char* options_json, char* report, int32_t cap);
int32_t caps_fasta_sequence(const char* text, char* seq, int32_t cap);

/* Solvation (v20, design/boards/SolvationBuilder). Options JSON {shape: 0 cubic | 1 rectangular | 2 solute plus padding,
   edge, edges: [a, b, c], padding (Å), tolerance (Å), solvent (an id of caps_solvent_library), water_model: "TIP3P" |
   "SPC/E" | "TIP4P/2005", density (g/cm³, 0: the solvent's), molecules (0: from the free volume), ion_mode: 0 none |
   1 neutralise | 2 concentration | 3 custom, salt, concentration (mol/L), cations, anions, seed}. The solute is the
   document's current frame (NULL: pure solvent), held at the centre.
   caps_solvent_library: {solvents: [{id, name, smiles, density, use}], salts: [{id, cation, anion, zc, za}]}.
   caps_solvate_plan: {ok, error, box: [a, b, c], box_volume, solute_volume, free_volume, solute_atoms, solute_charge,
   solvent, cations, anions, solvent_name, cation, anion, solvent_mass, concentration, density, notes} without packing.
   caps_solvate: the packed box as a new document; progress gets the stage (0 insertion, 1 optimisation, 2 verification),
   the loop, the smallest distance so far and the molecules still too close; return non-zero to cancel. */
typedef int32_t (*caps_stage_progress_fn)(int32_t stage, int32_t loop, int32_t loops, double dmin, int32_t bad, void* user);
int32_t caps_solvent_library(char* json, int32_t cap);
int32_t caps_solvate_plan(caps_doc* solute, const char* options_json, char* json, int32_t cap);
caps_doc* caps_solvate(caps_doc* solute, const char* options_json, caps_stage_progress_fn progress, void* user, char* report, int32_t cap);

/* Appearance of the Studio view (v20, design/boards/Appearance), kept by the document and prepared for each frame (not
   drawn while a Visualize pipeline runs). caps_set_appearance JSON {active, layers: [{expression ("" all atoms, else a
   pipeline expression such as "Molecule <= 3"), style: "ball_and_stick" | "sticks" | "wireframe" | "space_filling" |
   "polyhedra" | "ribbon" | "hidden" | "no_hydrogens"}] applied in order (atoms no layer names keep the view's style),
   colour: "" (the view's) | "element" | "molecule" | "type" | "distance" | "charge", ramp: "blue_orange" | "viridis" |
   "red_white_blue", surface: {kind: "none" | "accessible" | "vdw" | "excluded", probe, spacing (Å), opacity,
   expression (the atoms it wraps), colour: "potential" | "atom" | "uniform"}, charges: [one per atom] (optional: colour
   "charge" shows these instead of the structure's, a preview)}; 0 or −1 (caps_last_error).
   caps_appearance_info: {active, error, styles: {style: count}, charge: [min, max], preview (charges shown are the
   preview's), surface: {area, vertices,
   triangles, potential: [min, max] kcal/mol/e}, polyhedra: triangles}.
   caps_atom_labels: a JSON array with one label per atom of the frame; kind "element" | "rs" (CIP) | "type" | "charge" |
   "name".
   caps_project_atoms: after caps_render with the same camera and options, x, y (output pixels) and visibility (1 when
   the atom shows at its centre) of each frame atom into xyv (3 floats per atom); returns the atoms written. */
/* v48 labels (labels.hpp): caps_label_kinds {atom: [{id, title, group}], bond: [..]}; caps_atom_labels takes any atom kind
   (and "rs", "ez"); caps_bond_labels {pairs: [i, j, …], labels: [..], crossing: [bool: across the cell]} of a bond kind,
   the Field assignment's types and terms when the structure is typed. */
int32_t caps_label_kinds(char* json, int32_t cap);
/* v49 a repeat unit from a molecule in 3D (molinfo.hpp repeat_unit_smiles): head and tail atoms (0-based; a hydrogen, or
   a heavy atom that gives one of its hydrogens) → {ok, smiles (head's * first), error}. */
int32_t caps_repeat_unit_smiles(caps_doc* d, int32_t head, int32_t tail, char* json, int32_t cap);
/* v48 caps_type_table: [{type, label, mass, count, element (the most common), elements, from_mass (the element the mass
   points to)}] of the frame; caps_edit {"op": "type_element", type, element, mass?} sets every atom of a type (the type
   and its label kept). */
int32_t caps_type_table(caps_doc* d, char* json, int32_t cap);
int32_t caps_bond_labels(caps_doc* d, const char* kind, char* json, int32_t cap);
int32_t caps_set_appearance(caps_doc* d, const char* json);
int32_t caps_appearance_info(caps_doc* d, char* json, int32_t cap);
int32_t caps_atom_labels(caps_doc* d, const char* kind, char* json, int32_t cap);
int32_t caps_project_atoms(caps_doc* d, const caps_camera* cam, const caps_render_opts* opt, float* xyv, int32_t count);

/* Trajectory player (v20, design/boards/Trajectory). caps_trajectory_series: options JSON {molecule (0: the largest),
   dt_fs (fs per timestep), stride, log (a LAMMPS log path, optional)} → {ok, error, columns: [Frame, Timestep,
   Time (ps), Density, Volume, Rg (Å), Ree (Å), the log's thermo columns …], rows (null where a value is missing),
   molecule, ends: [first, last] (the chain's backbone ends), run_frames (the frame nearest each run's start),
   log_rows}. progress(done, total) returns non-zero to stop. caps_set_smoothing: positions shown averaged over
   `window` frames (1: off). */
int32_t caps_trajectory_series(caps_doc* d, const char* options_json, caps_series_progress_fn progress, void* user, char* json, int32_t cap);
void caps_set_smoothing(caps_doc* d, int32_t window);

/* Torsion scan (v20, design/boards/TorsionScan) of the current frame. Options JSON {atoms: [a, b, c, d] (0-based),
   from, to, step (degrees), relax (bool: the four atoms held and the rest minimised at each step), ftol, forcefield:
   "auto" (the Field assignment when complete, else UFF) | "uff"} → {ok, error, forcefield, points: [{phi, energy
   (kcal/mol above the lowest), dihedral, vdw, coulomb, angle, bond, total}], conformers: [{phi, energy, state}],
   barrier, phi_start, moving (atoms rotated), notes}. The geometry of each point is kept: caps_torsion_show(d, k) puts
   point k into the current frame, −1 restores the frame as it was before the scan. */
int32_t caps_torsion_scan(caps_doc* d, const char* options_json, caps_series_progress_fn progress, void* user, char* json, int32_t cap);
int32_t caps_torsion_show(caps_doc* d, int32_t index);
/* A torsion to scan when none is picked: four heavy atoms around the middle bond of the longest backbone (not in a
   ring); −1 when there is none. */
int32_t caps_default_torsion(caps_doc* d, int32_t* atoms);

/* Structure editing (v20, the Studio's builder tools, design/boards/ElementPicker, SelectionStereo, AddHydrogens) of a
   single-frame document, with undo. caps_edit JSON {op, …}: "element" {atoms: [i…] | "selection", element: "N"},
   "charge" {atoms, charge}, "add_atom" {to (−1: beside the structure), element, order, geometry (0 auto, 3 sp³, 2 sp²,
   1 sp), charge}, "bond" {i, j, order (9: coordinate / dative)}, "unbond" {i, j}, "set_coordination" {atom, geometry: linear |
   trigonal | tetrahedral | square_planar | trigonal_bipyramidal | square_pyramidal | octahedral}, "delete" {atoms}, "add_h" {atoms (absent: all)}, "invert"
   {centre}, "tacticity" {to: "isotactic" | "syndiotactic", clean}, "clean" {atoms, ftol}, "attach" {target, smiles (with *
   attachment points), which, clean, name}, "place" {smiles, name, resname}, "translate" {atoms | "selection", by: [dx, dy,
   dz]}, "fuse_ring" {i, j, clean} (a benzene ring onto the bond i–j), "set_geometry" {atoms: [i, j (, k (, l))], value (Å or
   degrees): the side of the last atom moves}, "rotate" {atoms | "selection", axis, degrees}, "mirror" {atoms | "selection",
   normal}, "set_rs" {centre, to: "R" | "S"}, "functionalize" {group (preset or SMILES with *), pattern: "random" | "all" |
   "band" | "helix" | "ends" | "edges" | "atoms" (with atoms | "selection"), elements, fraction, count, min_spacing, from,
   to, pitch, phase, side: "outer" | "inner" | "both", seed} (functionalize.hpp) → {ok, error, what, atoms,
   added: [new atom indices]}. A Field assignment is cleared by any edit. caps_undo(d, 0) undoes, (d, 1) redoes;
   caps_history: {undo: [what…], redo: […]}.
   caps_select JSON {mode: "smarts" | "element" | "type" | "charge" {lo, hi} | "within" {distance} | "grow" {steps} |
   "molecule" {atoms} | "indices" {atoms} | "expression" | "query" (v20, query.hpp: smarts "…" and chain 1-4 …; the reply
   adds rings) | "all" | "none", pattern, op: "replace" | "add" | "subtract" |
   "intersect" | "invert"} → {ok, error, count, matched}; the selection is ringed in the view. caps_selection:
   {count, indices}. caps_tacticity: {label, m, r, mm, mr, rr, centres, chains: [{centres, dyads}]}. */
int32_t caps_edit(caps_doc* d, const char* json, char* out, int32_t cap);
int32_t caps_undo(caps_doc* d, int32_t redo);
int32_t caps_history(caps_doc* d, char* json, int32_t cap);
/* v20 history (design/boards/History): caps_history also gives {steps: [{what, atoms, state: done|current|undone}],
   start_atoms, branches: [{steps: [what…], atoms}], snapshots: [{name, atoms, step, on_branch}]}; an edit after an undo
   keeps the undone steps as a branch. caps_snapshot JSON {op: "take" {name} | "restore" {index} | "delete" {index} |
   "save" {index, path} | "branch" {index} (make that branch the redo steps) | "drop_branch" {index}}. */
int32_t caps_snapshot(caps_doc* d, const char* json);
int32_t caps_select(caps_doc* d, const char* json, char* out, int32_t cap);
int32_t caps_selection(caps_doc* d, char* json, int32_t cap);
int32_t caps_tacticity(caps_doc* d, char* json, int32_t cap);

/* Elements (v20): the atomic number of a symbol (0 unknown); mass (g/mol), covalent radius (Cordero 2008), van der
   Waals radius (Bondi 1964) and display colour of an element; −1 for an unknown number. */
int32_t caps_element_number(const char* symbol);
int32_t caps_element_info(int32_t z, double* mass, double* covalent, double* vdw, uint32_t* rgb);

/* Interactions and checks (v20, design/boards/Interactions) of the current frame: options JSON {hb_distance (3.5 Å),
   hb_angle (30°), contact_margin (0.4 Å), clash_factor (0.75), show_hbonds, show_contacts, show_clashes} →
   {ok, error, hbonds, contacts, clashes, molecules, net_charge, hbond_list: [{donor, hydrogen, acceptor, distance,
   angle}], clash_list: [{i, j, distance}], issues: [{level, title, detail, fix: "push_apart" | "wrap" | "add_h" | "",
   atoms}]}; the H-bonds (dashed), contacts and clashes are drawn in the view until caps_clear_checks. */
int32_t caps_interactions(caps_doc* d, const char* options_json, char* json, int32_t cap);
void caps_clear_checks(caps_doc* d);

/* v20 polymer statistics (design/boards Polydispersity, Copolymer, Tacticity, BlendPhase, SolventScreen, Electrostatics).
   JSON in, JSON out ({ok, error} on failure).
   caps_chain_lengths {distribution: monodisperse | schulz-zimm | flory | poisson, nn, pdi, count, seed, m0, best_of} →
     {lengths[], sample {nn, mn, mw, pdi, min, max, sum}, target {nn, mn, mw, pdi}, k, curve {n[], number[], weight[]}}
   caps_copolymer {r1, r2, f1, dp, seed} → {F1, paa, pbb, run_a, run_b, azeotrope|null, curve {f1[], F1[]}, sequence[]
     (Grow's chain 0 at this seed), chain {F1, run_a, run_b, longest, a, b}}
   caps_stereo {model: bernoulli | markov, pm | p_mr, p_rm, dyads | dp + seed, measured[10]} → {model {pm, p_mr, p_rm, mm, mr,
     rr, pentads[]}, names[], chain {dyads, m, r, mm, mr, rr, pentads[]}, fit {bernoulli, bernoulli_rms, markov, markov_rms,
     mm_rr, mr2_4}}
   caps_blend_phase {na, nb, a, b, t, t_min | t_max} (χ = a + b/T) → {chi_c, phi_c, tc|null, kind, chi_t, coexist[2],
     spinodal[2], binodal {phi[], t[]}, spinodal_curve {phi[], t[]}}
   caps_solvent_chi {delta_polymer, t, solvents: [{name, v, delta}]} → {rt, solvents: [{name, chi, predicted}]}
   caps_ewald_params (doc may be NULL) {cutoff, tolerance, spacing, order, edges[3]} → {beta, beta_rc, table[], edges[],
     mesh[], spacing[], fits, curve {r[], erfc[]}} — edges from the document's cell when not given. */
int32_t caps_chain_lengths(const char* json, char* out, int32_t cap);
int32_t caps_copolymer(const char* json, char* out, int32_t cap);
int32_t caps_stereo(const char* json, char* out, int32_t cap);
int32_t caps_blend_phase(const char* json, char* out, int32_t cap);
int32_t caps_solvent_chi(const char* json, char* out, int32_t cap);
int32_t caps_ewald_params(caps_doc* d, const char* json, char* out, int32_t cap);
/* v20 each atom's residue number (0 = none; Grow numbers the repeat units along each chain from 1); returns the atom
   count, filling at most cap entries. */
int32_t caps_atom_residues(caps_doc* d, int32_t* out, int32_t cap);
/* v20 display (design/boards DisplayStyles, LensView, AddHydrogens, ModelResolution). The Backbone style draws tubes
   through each chain's main-chain atoms (side groups hidden, small molecules without H).
   caps_set_display {polar_h_only, selection_full, lens: {on, centre (atom), radius Å, inside, outside (styles 0–4), dim},
   clip: {on, axis (0 x, 1 y, 2 z), from, to (fractions of the cell along the axis, else of the structure's extent), invert}}:
     view only, the structure is untouched. caps_lens_inside: 1 when the atom is inside the lens (or no lens).
   caps_display_counts → {atoms, h, heavy, polar_h, chains, backbone_atoms, lens_atoms, lens_h}.
   caps_hydrogen_plan → {rows: [{label, atoms, hydrogens}], heavy, h, add, net_charge, aromatic_bonds, selection}: what
     the "add_h" edit would add (the selection when there is one).
   caps_doc_copy: a new document holding the current frame and its force-field assignment (its provenance carried over).
   caps_resolution_summary {per_bead} → {all_atom, united_atom, coarse_grained: {sites, hydrogens, mass}}.
   caps_resolution_convert {to: united-atom | coarse-grained, per_bead} → a new document (NULL on error). */
int32_t caps_set_display(caps_doc* d, const char* json);
int32_t caps_lens_inside(caps_doc* d, int32_t atom);
int32_t caps_display_counts(caps_doc* d, char* out, int32_t cap);
int32_t caps_hydrogen_plan(caps_doc* d, char* out, int32_t cap);
caps_doc* caps_doc_copy(caps_doc* d);
int32_t caps_resolution_summary(caps_doc* d, const char* json, char* out, int32_t cap);
caps_doc* caps_resolution_convert(caps_doc* d, const char* json, char* report, int32_t cap);
/* v20 "Make real" (design/boards/Replicate): the pipeline's current particles (replicas, deletions …) as a new document
   with unique identifiers; NULL without a pipeline result. */
caps_doc* caps_pipeline_materialize(caps_doc* d);
/* v20 how many particles of the current frame an expression selects (design/boards/ExpressionSelect): {ok, error,
   count, total, types: {"2": "ca", …}}. */
int32_t caps_expression_count(caps_doc* d, const char* expr, char* out, int32_t cap);

#ifdef __cplusplus
}
#endif
#endif
