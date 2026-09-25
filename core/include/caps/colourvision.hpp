// Colour-vision deficiency (design/boards/ColourVision): how a colour looks with protanopia, deuteranopia or
// tritanopia, after Machado, Oliveira & Fernandes (2009) — their severity-1.0 matrices applied to linear sRGB; lower
// severities blend linearly with normal vision (an approximation of their tabulated matrices) — and CIE 1976 ΔE*ab
// to find colours that become hard to tell apart.
#pragma once
#include <string>
#include <vector>

#include "caps/render.hpp"

namespace caps {

enum class Vision { Normal = 0, Protan = 1, Deutan = 2, Tritan = 3 };
const char* to_string(Vision v);   // normal, protanopia, deuteranopia, tritanopia

unsigned simulate_vision(unsigned rgb, Vision v, double severity = 1.0);   // 0xRRGGBB
void simulate_vision(Image& img, Vision v, double severity = 1.0);          // in place (alpha kept)
double delta_e76(unsigned a, unsigned b);                                   // CIE L*a*b* (D65) distance

struct NamedPalette {
  std::string name;
  std::vector<std::string> labels;
  std::vector<unsigned> colours;
};
struct VisionPair {
  std::string palette;
  Vision vision = Vision::Normal;
  int a = 0, b = 0;   // indices into the palette
  double de = 0;
};
// Every pair of each palette closer than threshold under any of the three deficiencies (or normal vision), closest first.
std::vector<VisionPair> confusable_pairs(const std::vector<NamedPalette>& palettes, double threshold = 12.0);

}  // namespace caps
