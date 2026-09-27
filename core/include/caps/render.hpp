// CPU ball-and-stick renderer: z-buffered sphere and capsule impostors, depth outlines,
// supersampled antialiasing with a true alpha channel (for transparent figures).
#pragma once
#include <cstdint>
#include <string>
#include <array>
#include <vector>

#include "caps/system.hpp"

namespace caps {

enum class Background { Dark, White, Transparent, Custom };
enum class ColourBy { Element, Molecule, Type, Property };
enum class Style { BallAndStick, SpaceFilling, Sticks, NoHydrogens, Backbone, Wireframe, Polyhedra, Ribbon, Hidden };
// Per-atom colour ramps for ColourBy::Property: viridis; blue–orange and red–white–blue diverging about zero.
enum class Ramp { Viridis, BlueOrange, RedWhiteBlue };

struct Mesh;   // caps/appearance.hpp
// A triangle mesh drawn with the atoms: one colour or per-vertex colours, translucent.
struct MeshDraw {
  const Mesh* mesh = nullptr;
  unsigned rgb = 0x8FB8D8;
  float opacity = 0.6f;
};

struct Camera {
  double yaw = 0.55, pitch = 0.40;     // radians
  double zoom = 1.0;                   // 1 = cell fills the view
  double pan_x = 0.0, pan_y = 0.0;     // Å in view plane
  bool perspective = false;
  double fov_deg = 35.0;
};

// A world-space tube (pipeline vectors, trajectory lines); an arrow ends in a cone.
struct Segment {
  Vec3 a{0, 0, 0}, b{0, 0, 0};
  unsigned rgb = 0xF5A524;
  double radius = 0.25;   // Å
  bool arrow = false;
};

struct RenderOptions {
  int width = 1280, height = 800;
  int supersample = 2;
  bool deep = false;                   // also keep 16-bit channels (Image::rgba16) from the supersampled average
  // Level of detail (design/boards/MillionAtoms), by distance from the focus (the view centre in its focal plane): atoms
  // nearer than lod_near are drawn in full with their bonds, those nearer than lod_far as spheres without bonds, the rest
  // as points. 0 turns it off.
  double lod_near = 0, lod_far = 0;   // Å
  Background background = Background::Dark;
  unsigned custom_rgb = 0x0F1113;
  ColourBy colour_by = ColourBy::Molecule;
  Style style = Style::BallAndStick;
  std::vector<double> property;        // per-atom values for ColourBy::Property
  double atom_scale = 0.28;            // × vdW radius for ball-and-stick
  double bond_radius = 0.14;           // Å
  bool outlines = true;
  bool depth_cue = true;
  bool show_cell = true;
  std::array<int, 3> cell_repeats{1, 1, 1};   // the cell is a supercell of these unit cells: its box dashed, one unit cell in the accent
  std::vector<char> faded;             // per atom: drawn faded toward the background (periodic images)
  float fade = 0.7f;                   // how far faded atoms move toward the background colour (0 … 1)
  std::vector<int> highlight;          // atom indices drawn with a selection ring
  int focus = -1;                      // atom drawn with the keyboard-focus ring (accent, outside any selection ring)
  bool ambient_occlusion = false;      // darken atoms by how little open sky they see (object space, per atom)
  std::vector<unsigned> colours;       // per atom 0xRRGGBB overriding colour_by (a pipeline's colours); 0xFFFFFFFF keeps it
  std::vector<Segment> segments;       // tubes and arrows drawn with the atoms
  std::vector<uint8_t> atom_style;     // per atom (Style values) overriding `style`: mixed styles; Ribbon and Hidden hide the atom
  Ramp ramp = Ramp::Viridis;           // ColourBy::Property
  bool symmetric = false;              // ColourBy::Property: the range ±max|value| (charges)
  double range_min = 0, range_max = 0; // ColourBy::Property: a fixed range when range_max > range_min
  std::vector<MeshDraw> meshes;        // surfaces and coordination polyhedra
};

struct Image {
  int width = 0, height = 0;
  std::vector<uint8_t> rgba;           // straight (non-premultiplied) alpha, row-major, top row first
  std::vector<uint16_t> rgba16;        // the same at 16 bits per channel when rendered with deep = true
};

struct RenderStats {
  size_t near = 0, mid = 0, far = 0;   // atoms drawn in each level of detail (all near when LOD is off)
  size_t bonds = 0;                    // bond halves drawn
};

// What the renderer draws, in world coordinates, for a GPU view (the view uploads it once and turns it cheaply):
// spheres, half-bond capsules and lines with their final colours (fading and ambient occlusion applied), the same
// atoms, styles, colours and radii as render(). Level of detail is not applied; meshes are not included (has_meshes).
struct Scene {
  std::vector<float> spheres;          // x y z r per sphere (Å)
  std::vector<uint32_t> sphere_rgb;    // 0xRRGGBB
  std::vector<int32_t> sphere_id;      // atom index
  std::vector<uint8_t> sphere_ring;    // 1 selection ring, 2 keyboard-focus ring, 3 both
  std::vector<float> capsules;         // ax ay az bx by bz r per capsule (Å)
  std::vector<uint32_t> capsule_rgb;
  std::vector<float> lines;            // ax ay az bx by bz per line
  std::vector<uint32_t> line_rgb;
  std::vector<float> line_width;       // pixels at 1×
  bool has_meshes = false;             // surfaces or polyhedra: a view that needs them draws on the CPU
  // What the camera fit uses (view_fit from these alone, without the structure): its centre, the cell corners when the
  // cell frames the view, the shown atoms' positions, the pad around them (Å), the perspective field of view.
  Vec3 fit_centre{0, 0, 0};
  std::vector<float> fit_corners;      // 8 × xyz, or none
  std::vector<float> fit_points;       // xyz per shown atom
  double fit_pad = 1.0, fov_deg = 35.0;
  uint32_t background = 0;             // 0xRRGGBB
  bool transparent = false, dark = true, depth_cue = true, outlines = true;
};

// The camera as render() fits it for an image of opt.width × opt.height × opt.supersample pixels: a world point p maps to
// r = R (p − centre) + pan (R: yaw about y, then pitch about x), screen x = w/2 + r.x·scale·k, y = h/2 − r.y·scale·k,
// k = dist / (dist − r.z) in perspective (else 1); zmin, zmax: the view depths of the shown atoms (depth cue).
struct ViewFit {
  double cos_yaw = 1, sin_yaw = 0, cos_pitch = 1, sin_pitch = 0;
  Vec3 centre{0, 0, 0};
  double scale = 1, w = 1, h = 1, pan_x = 0, pan_y = 0;
  bool perspective = false;
  double dist = 100, zmin = 0, zmax = 0;
};
ViewFit view_fit(const System& s, const Camera& cam, const RenderOptions& opt);

struct Renderer {
  RenderStats stats;                   // of the last render
  // The drawable scene (see Scene); ambient occlusion from the renderer's cache.
  Scene scene(const System& s, const RenderOptions& opt);
  // Picks the atom under a pixel of the last render (-1 if none).
  int pick(int x, int y) const;
  Image render(const System& s, const Camera& cam, const RenderOptions& opt);
  // For every atom: output-pixel x, y and 1 when the last render (same camera and options) shows it at its centre, else 0.
  std::vector<float> project(const System& s, const Camera& cam, const RenderOptions& opt) const;

 private:
  std::vector<int32_t> id_buffer_;
  int id_w_ = 0, id_h_ = 0;
  std::vector<float> ao_;              // per-atom accessibility of the last frame (camera independent), and its key
  double ao_key_ = 0;
  std::vector<float> scene_ao_;        // the same for scene()
  double scene_ao_key_ = 0;
};

// Per-atom ambient accessibility in [0, 1]: the share of 32 directions from each atom's surface that leave a 5 Å shell
// without meeting another atom (radius = max(r, 0.7 Å)); atoms with show[i] == 0 neither occlude nor are shaded.
std::vector<float> ambient_accessibility(const System& s, const std::vector<double>& radius, const std::vector<char>& show);

unsigned molecule_colour(int k);     // design palette, 10 entries cycled
unsigned viridis(double t);          // t in [0, 1]
unsigned ramp_colour(Ramp r, double t);   // t in [0, 1]
unsigned background_rgb(Background b, unsigned custom);

void write_png(const Image& img, const std::string& path);
std::vector<uint8_t> encode_png(const Image& img);

// PNG with 16-bit channels (from rgba16, else the 8-bit values widened), the dpi (pHYs), the sRGB chunk and text
// chunks (iTXt, UTF-8) such as a provenance manifest.
struct PngOptions {
  bool sixteen = false;
  double dpi = 0;
  bool srgb = true;
  std::vector<std::pair<std::string, std::string>> text;
};
std::vector<uint8_t> encode_png(const Image& img, const PngOptions& o);
void write_png(const Image& img, const std::string& path, const PngOptions& o);
// The text chunks (tEXt, iTXt without compression) of a PNG file.
std::vector<std::pair<std::string, std::string>> read_png_text(const std::string& path);

// Animated PNG, streamed frame by frame (8-bit RGBA; every frame the size of the first). Plays in browsers and most
// image viewers; the first frame is the still image for viewers without APNG.
class ApngWriter {
 public:
  ApngWriter(const std::string& path, int fps, int loops = 0, const PngOptions& o = {});
  ~ApngWriter();
  void add(const Image& frame);
  void close();   // writes IEND and the frame count
  int frames() const { return frames_; }

 private:
  struct Impl;
  Impl* impl_;
  int frames_ = 0;
};
// Vector figure: painter-sorted spheres with shading gradients; no background shape when transparent.
std::string render_svg(const System& s, const Camera& cam, const RenderOptions& opt);
// Pixels per Å at the focal plane for an opt.width × opt.height image (exact in orthographic views): scale bars.
double view_scale(const System& s, const Camera& cam, const RenderOptions& opt);

}  // namespace caps
