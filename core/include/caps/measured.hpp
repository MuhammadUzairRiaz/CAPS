// CAPS readers for measured data and grids: the files instruments and other programs write, read into one curve (x, y,
// σ) or one grid.
//
//  Curves (read_measured, by content first, then by extension)
//    JCAMP-DX     IR and other spectra (McDonald & Wilks, Appl. Spectrosc. 42, 151 (1988)): ##XYDATA=(X++(Y..Y)) in
//                 AFFN, PAC or the compressed ASDF forms (SQZ, DIF, DUP, with the Y check value that repeats a DIF line's
//                 last ordinate), and ##XYPOINTS=(XY..XY); x from FIRSTX to LASTX over NPOINTS, y × YFACTOR.
//    XRDML        PANalytical XML: the first scan's 2Theta positions (start and end, or listed) and its counts (2.x) or
//                 intensities (1.x), as written; the counting time noted.
//    GSAS         constant-wavelength banks (BINTYP CONS): STD (σ = √(Y/NCTR)), ESD (Y with its σ), FXY and FXYE (x y [σ]);
//                 x in centidegrees, given in degrees. The first bank is read.
//    UXD          Bruker / Siemens text: a range's _START and _STEPSIZE with single-column _COUNTS / _CPS data, or the
//                 two-column _2THETACOUNTS / _2THETACPS / _2THETACOUNTSTIME data. The first range is read.
//    Rigaku       the *TYPE text format (.asc, .dat): a group's *START, *STEP and *COUNT values. The first group is read.
//    columns      delimited text. A header line whose names match the columns picks x (2θ, x, ...), y (intensity,
//                 counts, y, ...) and σ (sigma, esd, ...) by name; a q column is never the intensity, and a header with
//                 no intensity column is refused. Without a header, two or three columns (x, y[, σ]) in order. Lines
//                 starting with # or ! and other non-numeric lines skipped.
//  Grids (read_cube)
//    .cube        Gaussian cube: origin, three voxel vectors (bohr, or Å for a negative count), the atoms, the values
//                 with the third index fastest; an orbital list (negative atom count) skipped, the first value of each
//                 point kept. The grid's cell spans n voxel vectors from the first point less half a voxel (as CAPS
//                 writes its own grids).
#pragma once
#include <string>
#include <vector>

#include "caps/json.hpp"
#include "caps/pipeline.hpp"
#include "caps/system.hpp"

namespace caps {

struct MeasuredCurve {
  std::vector<double> x, y, sigma;   // sigma empty: none given
  std::string format, title, xunits, yunits;
  std::vector<std::string> notes;
};

MeasuredCurve read_measured(const std::string& path);
// The same from text; `name` (a file name) breaks ties between formats by its extension.
MeasuredCurve parse_measured(const std::string& text, const std::string& name = "");
// A Gaussian cube file: the grid and, when asked, its atoms (Å).
GridField read_cube(const std::string& path, System* atoms = nullptr);
GridField parse_cube(const std::string& text, System* atoms = nullptr);
Json measured_json(const MeasuredCurve& c);

}  // namespace caps
