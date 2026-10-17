#pragma once
// caps ff complete: the bonded terms a force field lacks for a structure, found and fitted by CAPS itself.
//
// The structure is typed and parameterised with the force field; every bond, angle and proper dihedral it has no rule
// for is located in the structure, a fragment is cut round one occurrence (whole small rings, no bond between two
// unsaturated atoms cut, one H cap per cut bond, periodic images made whole; widened until its atoms keep their types),
// the coordinate is scanned with CAPS's machine-learned potential (each point relaxed with the coordinate held), and the
// term is fitted to the potential's energy minus CAPS's MM energy without it: harmonic bonds and angles, OPLS Fourier
// V1-V4 torsions (every missing class about the scanned bond together). A fit outside the physical range or with a poor
// RMS is not written; a fragment with an element the potential does not cover is reported, never guessed.
// The result is an overlay (FFDef with the fitted rules, each with its provenance in its comment) and a report.
#include <functional>
#include <string>
#include <vector>

#include "caps/ffdef.hpp"
#include "caps/system.hpp"

namespace caps {

struct CompleteOptions {
  std::string model;           // the machine-learned potential (data/ml/caps-mace-v1.model)
  int step = 15;               // degrees between torsion scan points
  int radius = 2;              // bonds round the term's atoms the fragment starts from
  double fmax = 0.01;          // eV/Å, relaxation convergence
  int max_steps = 400;
  double torsion_rms_max = 1.0;   // kcal/mol: a poorer torsion fit is reported, not written
  double stretch_rms_max = 0.5;   // kcal/mol, bonds and angles
  double window = 30.0;           // kcal/mol: a torsion point this far above its scan's lowest, in the model or in MM
                                  // without the term, is a clash rather than the torsion's profile and is not fitted
};

struct CompleteTerm {
  std::string kind;                  // bond, angle, dihedral
  std::vector<std::string> classes;  // the class names the force field looks the term up with
  std::string status;                // fitted, rejected, not_fitted
  std::vector<double> params;        // harmonic: K, r0 (Å) or theta0 (deg); opls: V1..V4
  double rms = 0;                    // kcal/mol
  int points = 0, fragment_atoms = 0;
  std::string why;                   // for rejected / not_fitted
  // the scan, point by point: the scanned value (Å or degrees), the model's energy and the MM energy without the term
  // (kcal/mol, each relative to its own lowest point), the relaxation's convergence
  struct Point { double value, e_ml, e_mm; bool converged; bool used = true; };   // used: in the fit (see window)
  std::vector<Point> scan;
};

struct CompleteReport {
  std::vector<std::string> missing_before, missing_after;
  std::vector<CompleteTerm> terms;
  bool complete = false;
  FFDef overlay;                     // the fitted rules only (merge_forcefield onto the force field)
};

CompleteReport complete_forcefield(const System& s, const FFDef& ff, const CompleteOptions& opt,
                                   const std::function<void(const std::string&)>& log = {});
// with the structure's types given (an assignment typed by hand, or a .car file's types kept)
CompleteReport complete_forcefield(const System& s, const FFDef& ff, const std::vector<std::string>& types, const CompleteOptions& opt,
                                   const std::function<void(const std::string&)>& log = {});

// ------------------------------------------------------------------ completion during an assignment
// Opt-in ("caps ff apply --complete ml", a recipe's type: { complete: ml }, the Field step's "Fill missing terms with
// CAPS ML"): after typing, the bonded terms the force field lacks are fitted by caps ff complete and merged onto it, and
// the assignment goes on. Each fitted rule carries its source (FFRule::source, ml_source), which follows its terms into
// the exported files. A term the model cannot fit, or whose fit fails the quality gate, stays missing with its reason;
// charges and untyped atoms are never ML's to fill.

// whether this CAPS has its ML engine (a build with -DCAPS_WITH_TORCH=ON)
bool ml_engine_available();
// what a build without it says when completion is asked for
std::string ml_engine_missing_message();

// CAPS's model, data/ml/caps-mace-v1.model: $CAPS_ML_MODEL, else under $CAPS_HOME, near the given paths (a force field's
// file in data/forcefields: its library's data/ml), beside the CAPS library or program (and its ../Resources: the macOS
// app), or in the source tree. "" when not found.
std::string find_ml_model(const std::vector<std::string>& near = {});
// The model to use: `given` when not empty, else find_ml_model(near); throws std::runtime_error naming where to put it
// when there is none. Shared by caps ff complete (the CLI, caps_ff_complete) and completion during an assignment.
std::string resolve_ml_model(const std::string& given, const std::vector<std::string>& near = {});

struct MlFill {
  int fitted = 0;                          // rules CAPS ML added
  std::vector<std::string> fitted_terms;   // "dihedral CA-CT-CT-CM · CAPS ML (...), RMS 0.31 kcal/mol"
  std::vector<std::string> left;           // what is still missing, with the reason ("dihedral X-Y-Z-W: a poor fit (RMS 1.40 kcal/mol)")
  CompleteReport report;
  std::string model;                       // the model file used
};

// Fits what `def` lacks for `s` (typed `types`) with CAPS ML and merges the fitted rules onto `def`. Throws
// std::runtime_error with ml_engine_missing_message() in a build without the ML engine, and when the model is not found.
MlFill fill_missing_with_ml(const System& s, FFDef& def, const std::vector<std::string>& types, const CompleteOptions& opt,
                            const std::function<void(const std::string&)>& log = {});
// The completion option's values: "ml" (fill with CAPS ML) or "none" / "off" / "" (report what is missing). Throws
// std::invalid_argument for anything else.
bool parse_complete_option(const std::string& value);

// The report as JSON, the same for caps ff complete --report and the C API (caps_ff_complete): {structure, forcefield,
// model, missing_before[], missing_after[], terms[{kind, classes[], status, params[], rms_kcal, points, fragment_atoms,
// why, scan[[value, E_model, E_MM without the term, converged, used]]}], complete}.
std::string complete_report_json(const CompleteReport& rep, const std::string& structure, const std::string& forcefield,
                                 const std::string& model);

}  // namespace caps
