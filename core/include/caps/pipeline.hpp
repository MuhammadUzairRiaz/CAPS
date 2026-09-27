// CAPS visualize pipeline (design/boards/VisPipeline, PipelineSteps, ColourBy, DataInspector): non-destructive steps
// applied to one frame, bottom to top as they are listed. Each step reads and writes per-particle properties, the
// selection and colours, adds global attributes and data tables, and may delete or replicate particles. The source
// file is never changed; the result is what the view draws and the data inspector lists.
//
// Steps (type, parameters):
//   select_expression   expression                       Type == 2 && Position.Z > 13 · Element == "O"
//   invert_selection · clear_selection
//   expand_selection    mode bonds|cutoff, iterations, cutoff
//   delete_selected
//   slice               normal [x,y,z], distance (default: through the cell centre), width, invert, select_only
//   colour_coding       property, mode auto|categorical|continuous, map viridis|diverging, start, end, lighten_h, only_selected
//   assign_colour       colour "#RRGGBB", keep_selection
//   cluster             mode bonds|cutoff, cutoff, heavy_only, unit atoms|molecules, only_selected, sort_by_size, colour, sweep
//                       → Cluster, tables clusters (with the molecules each holds) and cluster_sweep (cutoff mode)
//   coordination        cutoff, rmax, bins, element_a, element_b (0: any), inter_only, only_selected, average_frames,
//                       every → Coordination (within cutoff), table rdf (to rmax; averaged over frames when asked)
//   compute_property    name, expression, only_selected (Position.X/Y/Z, Charge and Selection write through)
//   wrap                mode atoms (each atom folded into the cell) | molecules (each molecule whole, moved by
//                       lattice vectors so its centre of mass is inside)   → table outside (image shifts)
//   replicate           nx, ny, nz, adjust_cell
//   histogram           property, bins, start, end, only_selected, stack_by (Type, Element …) → table histogram
//   binning             property, axis 0|1|2, bins, reduction mean|sum|density             → table binning;
//                       axis2 (0|1|2), bins2: a map over two axes → table binning2d (x, y, value, count)
//   create_bonds        mode perceive|cutoff|pairs, cutoff, pairs {"C-C": 1.7, …}, tolerance, inter_only, only_selected,
//                       keep_file (default: a file's bonds stay and the new ones are compared with them), replace
//   transparency        mode selected (value) | property (property, start, end → low … high), only_selected → Transparency
//   particle_radius     the same → Radius (Å); the view draws both
//   freeze_property     property, frame (reference, default 0), output (default "<property> frozen"): the values it had
//                       at the reference frame, after the steps below this one, matched by particle identifier
//   orientation         axis director|x|y|z, radius, angle, neighbours → Orientation (P₂ per atom), Crystalline, Orientation.S
//   affine_transform    strain [εxx εyy εzz] or matrix (9, row-major), translation, target all|particles|cell, only_selected
//   cna                 adaptive common neighbour analysis → Structure Type (0 other, 1 FCC, 2 HCP, 3 BCC, 4 ICO), table structures
//   centrosymmetry      neighbours (even; 12 FCC, 8 BCC) → Centrosymmetry (Å²)
//   ptm                 rmsd_max, fcc|hcp|bcc|ico|sc     polyhedral template matching (Larsen 2016) → Structure Type (as cna,
//                       5 SC), RMSD, Interatomic Distance, Shear Strain, Orientation.W/X/Y/Z, tables structures, rmsd
//   wigner_seitz        reference frame|file, frame, path, output particles|sites → Site Index, Occupancy; attributes
//                       WignerSeitz.vacancy_count, .interstitial_count; table defects
//   combine             path, frame                      the particles and bonds of a second file appended
//   unwrap              method bonds (whole along bonds) | images (the file's image flags) | nojump (each atom
//                       followed through the frames by its shortest step)  → table images, Unwrap.bonds_split
//   molecule_shape      per molecule: Rg, κ², asphericity, mass → properties and table molecules
//   topology            bond lengths, bond angles and dihedrals as histograms          → tables bonds, angles, dihedrals
//   displacements       reference first|previous|frame, frame, subtract_drift   → Displacement(.X .Y .Z), MSD, table
//                       displacements (per molecule: mean and max |d| of heavy atoms, COM shift)
//   smooth              window (frames, centred): positions averaged over the trajectory; C–C lengths stored and
//                       averaged and the RMS shift as attributes; the frame is marked Smoothed
//   vectors             property end_to_end|dipole|displacement|velocity, scale, radius       → arrows, table vectors
//   trajectory_lines    particles centres|selected, from, to, stride, radius   → paths over the frames, table paths
//                       (path length, net shift, net / path)
//   primitive_paths     radius, max_steps   chains pulled tight between fixed ends without crossing (Everaers 2004) → lines,
//                       attributes PrimitivePath.Ne (modified S-coil), .Ne_coil, .Z, .Lpp, .a_pp, table primitive_paths
//   voids               probe, grid, show (points)       accessible volume, voids by volume, probe sweep (needs a cell)
//   voronoi             method exact|exact_radical|grid|radical → AtomicVolume (+ exact: Coordination, Max Face Order, Voronoi Index.3–6, table voronoi_indices); cells sum to the box
//   density_field       grid, sigma, axis, position       Gaussian mass density: mean, empty share, profile, slice points
//   msd                 heavy_only, every, max_lag, fit_from, fit_to (lags, frames), timestep_fs
//                       MSD(τ) of atoms and chain centres, D from the centres over the fit lags
//   scatter             x, y, only_selected                         → table scatter (points), Pearson r
//   python              file (a script with an @step function) or code (its text), timeout — run in a Python process through the caps
//                       package in data/python ($CAPS_PYTHON_PATH; interpreter $CAPS_PYTHON, else python3): the frame's
//                       particles, bonds and attributes go in; attributes, tables, properties and a selection come back
//
// Expressions: numbers, "C" (an element, for Element comparisons), particle properties (Identifier, Index, Molecule,
// Type, Element, Mass, Charge, Position.X/Y/Z, Selection, DistanceToCOM, Monomer, MoleculeCOM.X/Y/Z, any computed
// property), + - * / % ^, == != < <= > >=, && || ! (also and, or, not), abs sqrt exp log min max floor ceil round,
// parentheses. Vectors: a name with .X .Y .Z parts (Position, MoleculeCOM(MoleculeIdentifier)) adds, subtracts and
// scales, and norm(v), dot(a, b) and the components make numbers of them: norm(Position − MoleculeCOM(MoleculeIdentifier)).
#pragma once
#include <array>
#include <cstdint>
#include <functional>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "caps/json.hpp"
#include "caps/render.hpp"
#include "caps/system.hpp"

namespace caps {

constexpr unsigned kNoColour = 0xFFFFFFFFu;

struct DataTable {
  std::string name, title;
  std::vector<std::string> columns;
  std::vector<std::vector<double>> rows;
  bool points = false;   // a scatter: plot the rows as points, not a line
  // an optional text column (one entry per row), shown after the first column: which molecules a cluster holds …
  std::string label_column;
  std::vector<std::string> labels;
};

struct PipelineLegend {
  std::string property;
  bool continuous = false;
  double lo = 0, hi = 0;
  std::string map = "viridis";
  std::vector<std::pair<std::string, unsigned>> entries;   // categorical: label, colour
};

struct StepStatus {
  std::string type, title, summary;
  std::string level = "ok";   // ok | warning | error | off
  std::string output;         // what the step printed (a Python step's console)
};

struct PipelineState {
  System system;                                       // the particles after the steps
  std::vector<int> origin;                             // index in the source frame of each particle
  std::vector<char> selected;
  std::vector<unsigned> colour;                        // 0xRRGGBB, or kNoColour for the view's own colouring
  std::map<std::string, std::vector<double>> props;    // computed per-particle properties
  std::vector<std::pair<std::string, double>> attributes;
  std::vector<DataTable> tables;
  std::vector<StepStatus> steps;                       // one per step, in the listed order
  PipelineLegend legend;
  bool has_legend = false;
  std::vector<Segment> segments;                       // arrows and paths the view draws with the particles
  int frame = 0;
  int64_t timestep = 0;
  const Trajectory* traj = nullptr;                    // the whole trajectory, for steps that read other frames
  const struct Pipeline* pipeline = nullptr;           // the steps being run and the one running (freeze re-runs those
  size_t step_index = 0;                               // below it on the reference frame)

  double attribute(const std::string& name, double def = 0) const;
  void set_attribute(const std::string& name, double v);
  size_t selected_count() const;
};

struct PipelineStep {
  std::string type;
  bool enabled = true;
  Json params = Json::object();
};

struct Pipeline {
  std::vector<PipelineStep> steps;   // top first, as listed; evaluated bottom to top
};

// {"steps": [{"type": …, "enabled": …, <parameters>}]} or the bare array.
Pipeline pipeline_from_json(const Json& j);
Json pipeline_to_json(const Pipeline& p);
// YAML (design/boards/SavePipeline): steps in the order they run, one flow mapping each; name and source optional.
std::string pipeline_to_yaml(const Pipeline& p, const std::string& name = "", const std::string& file = "", const std::string& topology = "");
Pipeline pipeline_from_yaml(const std::string& text, std::string* name = nullptr, std::string* file = nullptr, std::string* topology = nullptr);
// Types known to run_pipeline, with a title and a one-line description each.
std::vector<std::array<std::string, 3>> pipeline_step_catalogue();
std::string step_title(const std::string& type);

// traj (optional) lets displacements and smoothing read other frames; frame is traj's frame frame_index as shown.
PipelineState run_pipeline(const System& frame, const Pipeline& p, int frame_index = 0, int64_t timestep = 0, const Trajectory* traj = nullptr);

// Time series (design/boards/TimeSeries): the pipeline on every stride-th frame, one row per frame with every numeric
// global attribute (Frame, Timestep, then the attributes in order). progress(done, total) returns false to stop.
DataTable pipeline_series(const Trajectory& traj, const Pipeline& p, int stride = 1, bool wrap = false,
                          const std::function<bool(int, int)>& progress = {});

// Per-particle property names available to expressions and colour coding.
std::vector<std::string> property_names(const PipelineState& st);
// Values of a named per-particle property; false if unknown.
bool property_values(const PipelineState& st, const std::string& name, std::vector<double>& out);
// One value per particle; throws std::invalid_argument with the position of a syntax error or an unknown name.
std::vector<double> evaluate_expression(const PipelineState& st, const std::string& expr);

// Result for the Studio: attributes, step status, tables, legend, property names, counts.
Json pipeline_result_json(const PipelineState& st);
// Particles matching filter (an expression; empty: all), rows offset … offset + count, with their columns.
Json particles_json(const PipelineState& st, const std::string& filter, size_t offset, size_t count);
Json bonds_json(const PipelineState& st, size_t offset, size_t count);

}  // namespace caps
