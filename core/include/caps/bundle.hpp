// CAPS figure bundles (design/boards/FigureBundle): a zip that carries a figure with the data behind it, the pipeline
// as run, the input structure and a provenance record with the sha256 of every file — and reproduce(), which rebuilds
// the data from the bundle alone and compares the hashes.
//
//   figure.png, figure.svg        the view (no overlays), at the recorded size
//   data/<table>.csv              every data table of the pipeline, columns as named, values %.10g
//   pipeline.json                 the steps as run
//   input/<file> [+ topology]     the structure, when included
//   provenance.json               CAPS version, input (name, sha256, atoms, frame), box, pipeline sha256, files and hashes
//   README.txt                    what is inside and how to reproduce it
#pragma once
#include <cstdint>
#include <map>
#include <string>
#include <vector>

#include "caps/pipeline.hpp"
#include "caps/render.hpp"
#include "caps/system.hpp"

namespace caps {

std::string sha256_hex(const std::string& bytes);

// Minimal zip (deflate, CRC-32) for bundles: names to contents.
void write_zip(const std::string& path, const std::vector<std::pair<std::string, std::string>>& files);
std::map<std::string, std::string> read_zip(const std::string& path);

struct BundleOptions {
  std::string name = "figure";          // bundle stem
  std::string input, topology;          // paths of the structure (and its data file for a dump)
  int frame = 0;
  bool include_input = false, include_pipeline = true, include_data = true, include_readme = true, include_figures = true;
  int width = 1920, height = 1080;
};

struct BundleFile {
  std::string name, content_note, bytes;   // bytes: the file's content
};

// The bundle's files for a trajectory, its pipeline and a camera (figures only when opt.include_figures).
std::vector<BundleFile> bundle_files(const Trajectory& traj, const Pipeline& pipeline, const BundleOptions& opt, const Camera& cam,
                                     const RenderOptions& render);
void write_bundle(const std::string& path, const std::vector<BundleFile>& files);

// Rebuilds the data files from the bundle's input and pipeline and compares their sha256 with the provenance record.
// Returns true when every data file matches; report lists each file as "match" or "differs".
bool reproduce_bundle(const std::string& path, std::vector<std::string>& report);

// The CSV of a data table, as the bundle writes it.
std::string table_csv(const DataTable& t);

}  // namespace caps
