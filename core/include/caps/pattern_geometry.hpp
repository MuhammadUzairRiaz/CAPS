// CAPS geometry measured at every match of a SMARTS pattern, over the frames of a trajectory.
//
//  Points     the pattern's atoms carrying atom-map numbers ([S:1][S:2]); atoms sharing a number form one point at
//             their centroid ([c:1]1[c:1][c:1][c:1][c:1][c:1]1 is a ring's centre), each group made whole by minimum
//             images from its first atom. 2 points: a distance (Å); 3: the angle at the middle one (°); 4: the
//             torsion about the middle pair (°, IUPAC sign, −180 … 180).
//  Matches    every embedding of the pattern; two that measure the same points in reverse order (S–S read either way,
//             a torsion read from either end) count once.
//  Result     every value over matches and frames: mean, standard deviation, a histogram, and the mean per frame.
#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "caps/system.hpp"

namespace caps {

struct PatternGeometry {
  std::string kind;                                        // distance | angle | torsion
  std::string unit;                                        // Å | °
  std::vector<std::vector<std::vector<uint32_t>>> matches; // per match: per point, its atoms
  std::vector<double> values;                              // frame-major: frame f, match m at f · matches + m
  double mean = 0, sd = 0, min = 0, max = 0;
  std::vector<double> hist_x, hist_y;                      // bin centres, fraction per unit
  std::vector<double> frame_mean;
  std::vector<std::string> notes;
};

// bin: histogram width (0: 0.01 Å for distances, 1° for angles and torsions)
PatternGeometry pattern_geometry(const Trajectory& t, const std::vector<size_t>& frames, const std::string& smarts, double bin = 0);

// The atom-map number of each pattern atom in the SMARTS text (0: none), in pattern order, and the text without them.
std::vector<int> smarts_map_numbers(const std::string& smarts, std::string* stripped = nullptr);

}  // namespace caps
