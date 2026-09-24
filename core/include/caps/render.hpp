// CPU ball-and-stick renderer: z-buffered sphere and capsule impostors, depth outlines,
// supersampled antialiasing with a true alpha channel (for transparent figures).
#pragma once
#include <cstdint>
#include <string>
#include <vector>

#include "caps/system.hpp"

namespace caps {

enum class Background { Dark, White, Transparent, Custom };
enum class ColourBy { Element, Molecule, Type, Property };
enum class Style { BallAndStick, SpaceFilling, Sticks, NoHydrogens, Backbone };

struct Camera {
  double yaw = 0.55, pitch = 0.40;     // radians
  double zoom = 1.0;                   // 1 = cell fills the view
  double pan_x = 0.0, pan_y = 0.0;     // Å in view plane
  bool perspective = false;
  double fov_deg = 35.0;
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
  std::vector<int> highlight;          // atom indices drawn with a selection ring
};

struct Image {
  int width = 0, height = 0;
  std::vector<uint8_t> rgba;           // straight (non-premultiplied) alpha, row-major, top row first
};

struct Renderer {
  // Picks the atom under a pixel of the last render (-1 if none).
  int pick(int x, int y) const;
  Image render(const System& s, const Camera& cam, const RenderOptions& opt);

 private:
  std::vector<int32_t> id_buffer_;
  int id_w_ = 0, id_h_ = 0;
};

unsigned molecule_colour(int k);     // design palette, 10 entries cycled
unsigned viridis(double t);          // t in [0, 1]
unsigned background_rgb(Background b, unsigned custom);

void write_png(const Image& img, const std::string& path);
std::vector<uint8_t> encode_png(const Image& img);
// Vector figure: painter-sorted spheres with shading gradients; no background shape when transparent.
std::string render_svg(const System& s, const Camera& cam, const RenderOptions& opt);

}  // namespace caps
