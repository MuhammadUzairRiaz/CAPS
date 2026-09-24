/* CAPS C ABI v12 — the stable boundary used by the Studio (P/Invoke) and other languages.
   Every function is exception-safe: errors are returned as codes and caps_last_error() explains them. */
#ifndef CAPS_C_H
#define CAPS_C_H
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

#define CAPS_ABI_VERSION 12  /* v2 relax, field; v3 md, trajectory; v4 equilibrate, chains; v5 pack; v6 react; v7 CAPS Field; v8 Analyze; v9 mechanics, Tg; v10 LAMMPS input; v11 convergence checks; v12 molecule builder */

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
  double contact_scale;            /* 1 = full contact limits */
  int32_t curve;                   /* allow gauche backbone torsions */
} caps_grow_opts;

/* Progress: (chains finished, chains, restarts, user) -> 0 to continue, non-zero to cancel. */
typedef int32_t (*caps_progress_fn)(int32_t done, int32_t total, int32_t restarts, void* user);

int32_t caps_abi_version(void);
const char* caps_last_error(void);

caps_doc* caps_open(const char* path, const char* topology_path);   /* NULL on error */
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
   to the document: types every atom, sets charges (charges: 0 from the force field, 1 Gasteiger, 2 keep the file's)
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
} caps_analyze_opts;
typedef int32_t (*caps_analyze_progress_fn)(const char* what, double fraction, void* user);
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

int32_t caps_export_png(caps_doc* d, const caps_camera* cam, const caps_render_opts* opt, const char* path);
int32_t caps_export_svg(caps_doc* d, const caps_camera* cam, const caps_render_opts* opt, const char* path);

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
caps_doc* caps_build_smiles(const char* smiles, const char* ff_path, const caps_build_opts* o, char* report, int32_t cap);

#ifdef __cplusplus
}
#endif
#endif
