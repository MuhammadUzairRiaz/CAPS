// CAPS's own ray tracer for publication images: the same scene the views draw (spheres, half-bond capsules, lines as thin
// tubes, surfaces and polyhedra as triangles, with their colours and transparency) traced through a bounding-volume
// hierarchy, with ambient occlusion, soft shadows from an area light, depth of field and true transparency. Any size:
// the image is traced in tiles on every core and kept at 16 bits per channel. The camera is the views' own (view_fit),
// so a traced image lines up with what the view shows.
#pragma once

#include <atomic>
#include <functional>

#include "caps/render.hpp"

namespace caps {

struct TraceOptions {
  int samples = 64;            // per pixel (antialiasing, soft shadows, occlusion and depth of field converge with them)
  bool ambient_occlusion = true;
  double ao_distance = 6.0;    // Å: occluders farther than this do not darken
  double ao_strength = 0.85;   // 0 … 1
  bool shadows = true;
  double light_size = 0.12;    // the key light's angular radius (radians): larger, softer shadows
  double aperture = 0.0;       // Å: the lens radius for depth of field (0: all in focus)
  double focus = 0.0;          // view depth (Å, 0 the view centre) that is sharp
  bool outlines = false;       // a thin ink rim at silhouettes (the CAPS view look)
  int max_layers = 8;          // transparent surfaces a ray passes through
  int threads = 0;             // 0: every core
  uint64_t seed = 1;
  // done (0 … 1) after each tile; return false to cancel
  std::function<bool(double)> progress;
};

// The scene (Renderer::scene) traced for a camera fitted to w × h (view_fit with those sizes, supersample 1). The image
// has straight alpha; rgba16 holds the full 16-bit result. Throws std::runtime_error when cancelled.
Image raytrace(const Scene& scene, const ViewFit& fit, int w, int h, const TraceOptions& o);

}  // namespace caps
