/* CAPS C ABI v19 — the stable boundary used by the Studio (P/Invoke) and other languages.
   Every function is exception-safe: errors are returned as codes and caps_last_error() explains them. */
#ifndef CAPS_C_H
#define CAPS_C_H
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define CAPS_ABI_VERSION 20  /* v2 relax, field; v3 md, trajectory; v4 equilibrate, chains; v5 pack; v6 react; v7 CAPS Field; v8 Analyze; v9 mechanics, Tg; v10 LAMMPS input; v11 convergence checks; v12 molecule builder; v13 palette, threads; v14 bench; v15 polymer builder; v16 electrostatics; v17 surfaces, interfaces, held molecule, inserted curatives; v18 progressive open, keyboard focus; v19 ambient occlusion, view scale; v20 space groups, crystal builder, peptides, solvation, appearance, trajectory player, torsion scan, editing, selections */

typedef struct caps_doc caps_doc;   /* an opened file: trajectory + current frame + renderer */

typedef struct {
  double yaw, pitch, zoom, pan_x, pan_y;
  int32_t perspective;
} caps_camera;

typedef struct {
  int32_t width, height, supersample;
  int32_t background;     /* 0 dark, 1 white, 2 transparent, 3 custom */
  uint32_t custom_rgb;
  int32_t colour_by;      /* 0 element, 1 molecule, 2 type, 3 property (distance to own molecule's centre) */
  int32_t style;          /* 0 ball & stick, 1 space filling, 2 sticks, 3 no hydrogens, 4 backbone */
  int32_t outlines, depth_cue, show_cell;
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
/* v19: what a file holds before opening it (format, first lines, dump columns, types, frames, bonds) as JSON. */
int32_t caps_inspect_file(const char* path, const char* topology_path, char* json, int32_t cap);
/* Staged open (Studio progressive open): stage 0 format detected, 1 frame 0 read, 2 topology joined, 3 frames read
   (fraction of the file); return non-zero to stop — the frames read so far are kept. max_frames > 0 stops early. */
typedef int32_t (*caps_open_progress_fn)(int32_t stage, double fraction, const char* detail, void* user);
caps_doc* caps_open_staged(const char* path, const char* topology_path, int32_t max_frames, caps_open_progress_fn progress, void* user);
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
} caps_relax_opts;

/* Relax progress: (stage, stages, iteration, energy kcal/mol, largest force, density, user) -> non-zero cancels. */
typedef int32_t (*caps_relax_progress_fn)(int32_t stage, int32_t stages, int32_t iteration, double energy, double fmax, double density,
                                          void* user);

/* Relax the current frame with GAFF (C/H in this version). On success the document holds one frame per stage
   (the start, each push-off / compression stage, the final structure) and shows the last one. Returns 0, 1 when
   the run finished without meeting the force tolerance, or -1 on error / cancel (the document is unchanged). */
int32_t caps_relax(caps_doc* d, const caps_relax_opts* o, caps_relax_progress_fn progress, void* user, char* report, int32_t report_cap);

typedef struct {
  double dt;                       /* fs */
  int64_t steps;
  double temperature;              /* K */
  int32_t thermostat;              /* 0 none (NVE), 1 Bussi, 2 Langevin */
  double tau_t;                    /* fs */
  int32_t barostat;                /* 0 none, 1 stochastic cell rescaling, 2 Berendsen */
  double pressure, tau_p;          /* atm, fs */
  int32_t new_velocities;          /* draw velocities even when the document has them */
  uint64_t seed;
  int32_t thermo_every, frame_every;
  double cutoff;
  int32_t coulomb, tail, threads;
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

/* Write every frame of the document as a LAMMPS text dump (id mol type xu yu zu). */
int32_t caps_save_trajectory(caps_doc* d, const char* path);

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
  int32_t thermostat;                      /* 1 Bussi, 2 Langevin */
  int32_t barostat;                        /* 1 stochastic cell rescaling, 2 Berendsen */
  double tau_t, tau_p;                     /* fs */
  uint64_t seed;
  double cutoff;
  int32_t coulomb, tail, threads;
  double frame_ps, thermo_ps;
  int32_t until_converged;                 /* production NPT blocks after the stages until the checks pass */
  double block_ps;
  int32_t max_blocks;
  double tol_density, tol_energy, tol_rg;  /* relative, kcal/mol per atom, relative */
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
} caps_react_opts;

typedef struct {
  int32_t cycle, reactions, total, clusters, atoms;
  double conversion, largest_fraction, reduced_mw, energy;
} caps_react_cycle;

typedef int32_t (*caps_react_progress_fn)(const caps_react_cycle* row, void* user);

/* React the current frame with the templates in `templates` (text). On success the document holds one frame per
   cycle and shows the last; returns 0, or -1 on error / cancel (the document is unchanged). */
int32_t caps_react(caps_doc* d, const char* templates, const caps_react_opts* o, caps_react_progress_fn progress, void* user, char* report,
                   int32_t report_cap);

/* Force-field summary of the current frame (the Field assignment, else the built-in GAFF typing of C and H): types,
   term counts and energy terms, as text. Returns 0 or -1. */
int32_t caps_field_info(caps_doc* d, char* text, int32_t cap);

/* CAPS Field. Assigns a force field (a caps-forcefield JSON file; typed by its own rules, or by rules_path when given)
   to the document: types every atom, sets charges (charges: 0 from the force field, 1 Gasteiger, 2 keep the file's, 3 QEq)
   and looks up every parameter. The types are written into the document (colour by type shows them). Relax,
   Dynamics, Equilibrate and LAMMPS data then use this force field; while atoms are untyped or parameters missing they
   refuse to run (CAPS never guesses parameters). Returns 0 when complete, 1 when something is missing, -1 on error. */
int32_t caps_field_assign(caps_doc* d, const char* ff_path, const char* rules_path, int32_t charges);
/* The assignment as JSON (atoms with type, rule, candidates, charge; counts; missing terms; types present with the
   viewer's colours; all types of the force field; notes; energy). Returns the length needed including the final NUL;
   call with json = NULL to size the buffer. Empty when nothing is assigned. */
int32_t caps_field_report(caps_doc* d, char* json, int32_t cap);
/* Sets atom `index` (0-based) to a type by hand; type "" or NULL removes the override. Returns as caps_field_assign. */
int32_t caps_field_override(caps_doc* d, int32_t index, const char* type);
/* Adds a parameter rule entered by hand: kind pair | bond | angle | dihedral | improper, the atom types (space-
   separated, as the missing-term list names them), a style ("" for the force field's default) and the parameters in
   that style's order. Such terms are reported as estimated. Returns as caps_field_assign. */
int32_t caps_field_add_rule(caps_doc* d, const char* kind, const char* types, const char* style, const char* params);
/* Merges parameters from another file (caps-forcefield .json or moltemplate .lt) over the force field. */
int32_t caps_field_import(caps_doc* d, const char* path);
/* Removes every imported and hand-entered parameter. */
int32_t caps_field_remove_rules(caps_doc* d);
/* Ends the assignment and restores the file's own types and charges. */
int32_t caps_field_clear(caps_doc* d);
/* Writes the assigned types, one per line (for caps ff apply --types). */
int32_t caps_field_types_file(caps_doc* d, const char* path);

/* Analyze: properties over the frames of the document (ids as for caps analyze: density, rdf, sq, xray, neutron, rg,
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

int32_t caps_summary_get(caps_doc* d, caps_summary* out);
int32_t caps_set_frame(caps_doc* d, int64_t frame);
/* 1: fold atoms into the cell (bonds across faces are hidden); 0: keep molecules whole (default). */
int32_t caps_set_wrap(caps_doc* d, int32_t wrap);
int32_t caps_atom(caps_doc* d, int32_t index, caps_atom_info* out);
int32_t caps_note_count(caps_doc* d);
const char* caps_note(caps_doc* d, int32_t k);

/* Renders into caller memory of width*height*4 bytes, RGBA with straight alpha. */
int32_t caps_render(caps_doc* d, const caps_camera* cam, const caps_render_opts* opt, uint8_t* rgba);
int32_t caps_pick(caps_doc* d, int32_t x, int32_t y);   /* atom index under pixel of last render, -1 none */
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
   "points": 200, "convergence": false, "colour": false} → {ok, total, area[] (Å² per atom), groups [{name, area, share}],
   convergence [{points, total, delta}]}; colour: true colours the view by exposure (area / full sphere).
   caps_set_atom_values colours the view by any per-atom quantity on a ramp (0 viridis, 1 blue–orange, 2 red–white–blue;
   n = 0 clears). */
int32_t caps_sasa(caps_doc* d, const char* json, char* out, int32_t cap);
int32_t caps_set_atom_values(caps_doc* d, const double* values, int32_t n, int32_t ramp);
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
/* v19 export (design/boards/ExportData): format lammps-data | lammps-dump (all frames) | gro | pdb | xyz | mol2;
   options JSON {"pipeline": bool (the Visualize result instead of the frame), "wrap": bool, "coeffs": bool (force-field
   sections in LAMMPS data)}. The preview writes to a scratch file and returns {"lines":[first n],"bytes","atoms","bonds",
   "angles","dihedrals","atom_types","bond_types","angle_types","dihedral_types","notes":[…]}. */
int32_t caps_export_data(caps_doc* d, const char* path, const char* format, const char* options);
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
   tables into `dir`. Sizes as caps_field_report: the length needed including the final NUL, or -1 on error. */
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
caps_doc* caps_nano_embed(const char* options_json, const char* spec_json, const caps_grow_opts* o, caps_progress_fn progress, void* user, char* report,
                          int32_t cap);
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
int64_t caps_held_molecule(const caps_doc* d);

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
   expression (the atoms it wraps), colour: "potential" | "atom" | "uniform"}}; 0 or −1 (caps_last_error).
   caps_appearance_info: {active, error, styles: {style: count}, charge: [min, max], surface: {area, vertices,
   triangles, potential: [min, max] kcal/mol/e}, polyhedra: triangles}.
   caps_atom_labels: a JSON array with one label per atom of the frame; kind "element" | "rs" (CIP) | "type" | "charge" |
   "name".
   caps_project_atoms: after caps_render with the same camera and options, x, y (output pixels) and visibility (1 when
   the atom shows at its centre) of each frame atom into xyv (3 floats per atom); returns the atoms written. */
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
   1 sp), charge}, "bond" {i, j, order}, "unbond" {i, j}, "delete" {atoms}, "add_h" {atoms (absent: all)}, "invert"
   {centre}, "tacticity" {to: "isotactic" | "syndiotactic", clean}, "clean" {atoms, ftol}, "attach" {target, smiles (with *
   attachment points), which, clean, name}, "place" {smiles, name, resname} → {ok, error, what, atoms,
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
   caps_set_display {polar_h_only, selection_full, lens: {on, centre (atom), radius Å, inside, outside (styles 0–4), dim}}:
     view only, the structure is untouched. caps_lens_inside: 1 when the atom is inside the lens (or no lens).
   caps_display_counts → {atoms, h, heavy, polar_h, chains, backbone_atoms, lens_atoms, lens_h}.
   caps_hydrogen_plan → {rows: [{label, atoms, hydrogens}], heavy, h, add, net_charge, aromatic_bonds, selection}: what
     the "add_h" edit would add (the selection when there is one).
   caps_doc_copy: a new document holding the current frame (its provenance carried over).
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
