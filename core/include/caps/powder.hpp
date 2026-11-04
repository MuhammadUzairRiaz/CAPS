// CAPS powder pattern: the diffraction pattern a powder of the cell would give on a laboratory diffractometer (or with
// neutrons), from the cell's reciprocal lattice. Works for a crystal's unit cell (sharp Bragg reflections) and for any
// periodic MD cell (a semicrystalline stack, a filled rubber), whose frames are averaged.
//
// Reflections: every reciprocal-lattice vector G = h b1 + k b2 + l b3 of the cell with |G| ≤ 4π sin θ_max / λ, its
// structure factor F(G) = Σ_j w_j(|G|) e^(−B s²) e^(iG·r_j) (w: Cromer–Mann X-ray form factors or neutron coherent
// lengths, s = |G|/4π, B a global isotropic temperature factor in Å²), |F|² averaged over the frames (the thermal disorder
// of an MD run gives its own Debye–Waller damping and diffuse scattering). Vectors of equal |G| make one reflection: its
// multiplicity is their number (systematic absences come out of F itself). Bragg's law sin θ = λ|G|/4π, per line.
// Integrated intensity per line: I = m |F|² P_MD LP, with
//   Lorentz–polarisation LP = (1 + cos² 2θ_m cos² 2θ) / ((1 + cos² 2θ_m) sin² θ cos θ) (X-ray; θ_m the monochromator's
//   Bragg angle, 0 none; a polarised synchrotron beam: numerator 1), neutrons 1 / (sin² θ cos θ);
//   March–Dollase preferred orientation P_MD = (r² cos² α + sin² α / r)^(−3/2) per vector, α its angle to the preferred
//   direction's reciprocal vector (Dollase, J. Appl. Cryst. 19, 267 (1986));
//   the Kα doublet: Kα2 at its own λ with ka2_ratio of Kα1's intensity.
// Peak shapes (area 1 in degrees 2θ):
//   Gaussian width H_G² = U tan² θ + V tan θ + W + (strain)² (Caglioti et al., Nucl. Instrum. 3, 223 (1958)), strain
//   4ε tan θ (Stokes & Wilson 1944); Lorentzian width H_L = X tan θ + Y / cos θ + size, size Kλ / (L cos θ) (Scherrer,
//   with L the crystallite size along G: one length, or one per cell axis as an ellipsoid);
//   "tch": the pseudo-Voigt of Thompson, Cox & Hastings (J. Appl. Cryst. 20, 79 (1987)), H and η from H_G and H_L;
//   "pseudo-voigt": the same H with a fixed η; "gaussian": all widths as one Gaussian (in quadrature);
//   "lorentzian": all widths as one Lorentzian (summed).
// 2θ is shifted by zero_shift. No axial-divergence asymmetry, absorption or anomalous dispersion.
#pragma once

#include <array>
#include <functional>
#include <string>
#include <vector>

#include "caps/system.hpp"

namespace caps {

struct PowderOptions {
  std::string radiation = "xray";   // xray | neutron
  std::string anode = "Cu";         // Cu Mo Co Fe Cr Ag (Kα1, Kα2); ignored when wavelength > 0
  double wavelength = 0;            // Å, one line (no doublet)
  bool kalpha2 = true;
  double ka2_ratio = 0.5;
  double monochromator_2theta = 0;  // degrees (0: none)
  bool polarised = false;           // a polarised (synchrotron) beam
  double two_theta_min = 5, two_theta_max = 90, step = 0.02;   // degrees
  std::string profile = "tch";      // tch | pseudo-voigt | gaussian | lorentzian
  double U = 0, V = 0, W = 0.0025;  // Caglioti, degrees²
  double X = 0, Y = 0;              // Lorentzian, degrees
  double eta = 0.5;                 // pseudo-voigt mixing (fixed)
  double size = 0;                  // crystallite size, Å (0: no size broadening)
  std::array<double, 3> size_abc{0, 0, 0};   // per cell axis, Å (all > 0: used instead of size)
  double scherrer_k = 0.9;
  double strain = 0;                // microstrain ε (0.001 = 0.1 %)
  double b_iso = 0;                 // global temperature factor B, Å²
  std::array<int, 3> preferred{0, 0, 0};   // preferred-orientation direction (hkl of the cell); 0 0 0: none
  double march_r = 1;               // March–Dollase r (< 1 platelets along it, > 1 needles)
  double zero_shift = 0;            // degrees 2θ
  int max_peaks = 200;              // reflections listed (strongest first, then by angle)
  int threads = 0;                  // 0: all cores
  std::function<bool(double)> progress;   // fraction; false stops
};

struct PowderPeak {
  int h = 0, k = 0, l = 0;          // one of the equivalent vectors (indices of the cell given)
  int multiplicity = 0;
  double d = 0, q = 0, two_theta = 0;   // Å, Å⁻¹, degrees (Kα1 or the single line, without the zero shift)
  double f = 0;                     // √⟨|F|²⟩ of one vector
  double intensity = 0;             // integrated, the strongest 100 (line 1)
  double fwhm = 0;                  // degrees 2θ
};

struct PowderReport {
  std::vector<double> wavelengths, weights;   // the lines, Å
  std::vector<PowderPeak> peaks;    // strongest first
  std::vector<double> two_theta, intensity, q;   // the pattern (intensity: the strongest point 100)
  int reflections = 0, vectors = 0, frames = 0;
  double seconds = 0;
  std::vector<std::string> notes;
};

// The wavelengths of an anode's Kα1 and Kα2 (Å); throws std::invalid_argument for an unknown one.
std::array<double, 2> anode_kalpha(const std::string& anode);

// t: a periodic structure; frames: those to average (empty: the first). Throws std::invalid_argument.
PowderReport powder_pattern(const Trajectory& t, const std::vector<size_t>& frames, const PowderOptions& o);

}  // namespace caps
