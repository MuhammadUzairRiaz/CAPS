// The views' scene written for other renderers and 3D tools, with the views' colours: POV-Ray (exact spheres and
// cylinders, the camera and key light of the view), glTF 2.0 binary (.glb: one mesh with vertex colours, for Blender,
// ParaView, three.js …) and Wavefront OBJ (the same triangles, colours per vertex). Transparent atoms keep their
// transparency in POV-Ray and glTF.
#pragma once

#include <string>

#include "caps/render.hpp"

namespace caps {

struct SceneExportReport {
  size_t spheres = 0, cylinders = 0, triangles = 0;   // what was written (glTF / OBJ: triangles of the tessellation)
};

// POV-Ray 3.7 scene: the camera of `fit` for a w × h image (orthographic or perspective), the key light, the background.
SceneExportReport write_povray(const Scene& sc, const ViewFit& fit, int w, int h, const std::string& path);
// glTF 2.0 binary. detail: sphere subdivisions (0 … 3; the default follows the number of atoms so large systems stay light).
SceneExportReport write_gltf(const Scene& sc, const std::string& path, int detail = -1);
SceneExportReport write_obj(const Scene& sc, const std::string& path, int detail = -1);

}  // namespace caps
