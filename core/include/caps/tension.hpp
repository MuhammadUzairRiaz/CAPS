#pragma once
// Surface and interfacial tension of a slab cell from its pressure tensor: for planar interfaces normal to an axis n, the
// Kirkwood–Buff integral γ = ½ ∫ [P_N(z) − P_T(z)] dz (J. Chem. Phys. 17, 338 (1949)) over a periodic cell with
// `interfaces` interfaces (a film with vacuum or a liquid slab has two) is γ = (L_n / interfaces) [P_nn − (P_t1 + P_t2)/2]
// with the cell's averaged pressure components. Each sample gives a γ; the mean and its error come from block averages.
// An isotropic long-range tail term (CAPS's or LAMMPS's homogeneous correction) cancels in P_N − P_T; the truncation's
// effect on γ itself is not corrected, so use a cut-off of 12 Å or more (or compare two cut-offs).
#include <string>
#include <vector>

#include "caps/dynamics.hpp"
#include "caps/trajectory.hpp"

namespace caps {

struct TensionSample {
  double pn = 0, pt1 = 0, pt2 = 0;   // atm: normal and the two tangential diagonal components
  double ln = 0;                     // Å: the cell's length along the normal
};
struct TensionResult {
  double gamma = 0, error = 0;       // mN/m, the error the standard error of the block means
  double pn = 0, pt = 0, ln = 0;     // means: atm, atm, Å
  int samples = 0, blocks = 0;
  std::vector<double> series;        // γ of each sample used, mN/m
  std::vector<std::string> notes;
};
// skip: the first share of the samples left out (equilibration); blocks: for the error (at least 2 samples each)
TensionResult surface_tension(const std::vector<TensionSample>& s, int interfaces = 2, double skip = 0.2, int blocks = 5);
// axis 0, 1, 2 (x, y, z): the normal
std::vector<TensionSample> tension_samples(const std::vector<ThermoRow>& rows, int axis = 2);
// A LAMMPS log with columns Pxx Pyy Pzz and Lz (or Zhi and Zlo; any case); the box length from cell_length (Å) when the
// log has none (NVT: constant).
std::vector<TensionSample> tension_samples(const ThermoLog& log, int axis = 2, double cell_length = 0);
// A file of either kind: CAPS's thermo CSV (caps md --log; columns pxx_atm … lz_A) or a LAMMPS log.
std::vector<TensionSample> tension_samples_from_file(const std::string& path, int axis = 2, double cell_length = 0);

}  // namespace caps
