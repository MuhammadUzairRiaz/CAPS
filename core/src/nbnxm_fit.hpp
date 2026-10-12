// Ewald real space for the cluster kernels: erf(κr)/r = κ gV(u) and the force part 2 κ³ dgV/du = κ³ gF(u), u = (κr)² ≤ 16,
// as P(t)/Q(t) with t = u/16 (coefficients from t⁰ up; Q(0) = 1). Rational fits made for CAPS (least squares reweighted
// towards minimax): maximum absolute error 3e-8 (gV) and 5e-8 (gF), below single-precision rounding.
#pragma once

namespace caps::nbnxm {

constexpr float kGvP[7] = {1.1283791379e+00f, 3.1401995010e+00f, 1.2811100603e+01f, 1.4467038445e+01f, 8.9639355099e+00f, -2.0071435048e+00f, 3.2804742414e-01f};
constexpr float kGvQ[5] = {1.0f, 8.1162471304e+00f, 2.9041673956e+01f, 5.7404020415e+01f, 5.9764272617e+01f};
constexpr float kGfP[7] = {-7.5225272615e-01f, 1.2926832369e+00f, -5.1479940366e+00f, 5.0524604436e+00f, -4.3528567753e+00f, 2.0075107894e+00f, -3.9860971314e-01f};
constexpr float kGfQ[5] = {1.0f, 7.8815427017e+00f, 2.7653014088e+01f, 5.3807163983e+01f, 5.6797690753e+01f};
constexpr float kFitUMax = 16.0f;   // the fits hold for (κ r_c)² up to this

// erfc by Abramowitz & Stegun 7.1.26 (|error| < 1.5e-7), for κ r_c beyond the fits: t = 1/(1 + p κr)
constexpr float kErfcP = 0.3275911f;
constexpr float kErfcA[5] = {0.254829592f, -0.284496736f, 1.421413741f, -1.453152027f, 1.061405429f};

}  // namespace caps::nbnxm
