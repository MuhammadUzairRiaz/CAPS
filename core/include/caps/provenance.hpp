// CAPS provenance (design/boards/Provenance): everything that produced a structure, in order — each step's engine,
// parameters, random-number generator and seed, and the papers its method comes from — with the inputs' sha256.
// A document keeps its chain; saving writes it beside the file (<file>.provenance.json, caps-manifest/1.0) and
// opening a file reads it back, so the chain follows the file through later sessions.
//
//   approximations()  the modelling choices in force (van der Waals, electrostatics, constraints, precision,
//                     estimated parameters), the latest step that set each winning
//   compare()         two manifests step by step: the parameters that differ, and whether inputs and versions agree
//   bibtex()          BibTeX of the cited methods (a built-in table of the methods CAPS implements)
#pragma once
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "caps/json.hpp"

namespace caps {

using KeyValues = std::vector<std::pair<std::string, std::string>>;

struct ProvStep {
  std::string engine;        // "relax.lbfgs", "dynamics.npt", "field.assign" …
  std::string summary;       // one line for people
  KeyValues params;          // in the order they matter
  std::string rng;           // "mt19937-64 · seed 20260923", or ""
  std::vector<std::string> cites;   // keys of the built-in bibliography
  KeyValues approximations;  // modelling choices this step puts in force (vdw, electrostatics, constraints, precision, estimated)
  std::string time;          // ISO 8601 UTC
};

struct Manifest {
  std::string generator = "CAPS 0.1.0";
  std::vector<ProvStep> steps;
  KeyValues inputs;          // file name → sha256
  bool deterministic = true; // every step seeded or deterministic
};

std::string now_iso();
Json manifest_json(const Manifest& m);
Manifest manifest_from_json(const Json& j);

std::string sidecar_path(const std::string& data_path);   // data_path + ".provenance.json"
void write_manifest(const Manifest& m, const std::string& data_path);
std::optional<Manifest> read_manifest(const std::string& data_path);

KeyValues approximations(const Manifest& m);

struct ManifestDiff {
  struct Row { int step; std::string engine, key, a, b; };
  std::vector<Row> rows;           // parameters (and seeds) that differ, by step
  std::vector<std::string> notes;  // steps only in one, different engines, different inputs or versions
  bool same_inputs = true, same_generator = true;
  int steps_a = 0, steps_b = 0, differing_steps = 0;
};
ManifestDiff compare(const Manifest& a, const Manifest& b);

// BibTeX entries for the keys (unknown keys are skipped); all_cites gathers the manifest's keys in order.
std::string bibtex(const std::vector<std::string>& keys);
std::vector<std::string> all_cites(const Manifest& m);
bool known_citation(const std::string& key);
// A reference for reading: "Authors, "Title", Journal Volume, Pages (Year). doi:…" from the built-in entry ("" if unknown).
std::string citation_text(const std::string& key);

}  // namespace caps
