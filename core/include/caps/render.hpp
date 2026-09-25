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
};

struct Renderer {
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
// Vector figure: painter-sorted spheres with shading gradients; no background shape when transparent.
std::string render_svg(const System& s, const Camera& cam, const RenderOptions& opt);
// Pixels per Å at the focal plane for an opt.width × opt.height image (exact in orthographic views): scale bars.
double view_scale(const System& s, const Camera& cam, const RenderOptions& opt);

}  // namespace caps
