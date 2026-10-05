// CAPS layer stacks (design/boards/SurfaceBuilder "Layer stack", "Add layer", "From CAPS Grow cell"): periodic cells
// piled along z into one cell — a crystal slab under a grown polymer cell, a polymer film between two slabs, a bilayer
// of two polymers, any rectangular cells from files.
//
// Every layer needs a rectangular cell (a along x, b along y, c along z; the Surface builder makes its slabs so). The
// first layer sets the lateral cell: the others are repeated na × nb times to come closest to it and then stretched
// (x and y scaled, positions affinely) to fit exactly; the strains are reported. With match "both", the first layer is
// repeated too (up to max_repeat) where that brings the others' strains lower. Crystals should come first: an amorphous
// polymer cell takes a few per cent of lateral strain and relaxes it away; a crystal does not. With match "average" the
// strain is shared (two crystals of similar spacing). A layer can be flipped (180° about x) and shifted in the plane.
//
// Inside each layer, molecules are made whole along z (bonds followed by minimum image) and put back by their centres,
// so a polymer cell's chains stay intact and stick out of its box a little; the layer's thickness is then its atoms'
// extent in z, and the next layer starts gap Å above the previous one's top atom. Along x and y atoms are wrapped
// into the new cell (the layer stays periodic there). vacuum > 0 leaves that much empty space on top; 0 makes the
// stack periodic in z with the same gap between the top layer and the image of the first.
//
// Molecules are renumbered layer by layer (the first layer's molecules first); atom types are merged by label and
// mass. Charges, names and residues are kept; velocities are dropped. Relax before dynamics: the gaps are a guess.
#pragma once

#include <string>
#include <vector>

#include "caps/system.hpp"

namespace caps {

struct StackLayerInput {
  std::string name;                // for the report
  const System* system = nullptr;  // a rectangular periodic cell
  bool flip = false;               // turned upside down (180° about x: the other face down, handedness kept)
  double shift_x = 0, shift_y = 0; // Å, moved in the plane (registry against the layer below)
};

struct StackOptions {
  double gap = 2.0;                // Å between one layer's top atom and the next one's lowest
  double vacuum = 0.0;             // Å above the stack; 0: periodic in z
  std::string match = "both";      // first: the first layer's cell as it is; both: it may be repeated too; average: the
                                   // lateral cell the mean of the layers' (repeated) cells, every layer strained a little
  int max_repeat = 6;              // repeats per direction considered
  double max_cell = 150.0;         // Å, the largest lateral edge considered
};

struct StackLayerReport {
  std::string name;
  int na = 1, nb = 1;              // repeats of the layer's own cell
  double strain_a = 0, strain_b = 0;   // fractional (x, y), 0 for the first layer
  double z_lo = 0, z_hi = 0;       // Å, where its atoms are in the stack
  size_t atoms = 0;
};

struct StackReport {
  std::vector<StackLayerReport> layers;
  double a = 0, b = 0, c = 0;      // the stacked cell, Å
  std::vector<std::string> notes;
};

System stack_layers(const std::vector<StackLayerInput>& layers, const StackOptions& o, StackReport* report = nullptr);

}  // namespace caps
