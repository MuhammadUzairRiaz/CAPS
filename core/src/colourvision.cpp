#include "caps/colourvision.hpp"

#include <algorithm>
#include <array>
#include <cmath>

namespace caps {

namespace {

using M3 = std::array<double, 9>;
// Machado, Oliveira & Fernandes, IEEE TVCG 15, 1291 (2009), severity 1.0.
const M3 kProtan = {0.152286, 1.052583, -0.204868, 0.114503, 0.786281, 0.099216, -0.003882, -0.048116, 1.051998};
const M3 kDeutan = {0.367322, 0.860646, -0.227968, 0.280085, 0.672501, 0.047413, -0.011820, 0.042940, 0.968881};
const M3 kTritan = {1.255528, -0.076749, -0.178779, -0.078411, 0.930809, 0.147602, 0.004733, 0.691367, 0.303900};

double to_linear(double c) { return c <= 0.04045 ? c / 12.92 : std::pow((c + 0.055) / 1.055, 2.4); }
double to_srgb(double c) {
  c = std::clamp(c, 0.0, 1.0);
  return c <= 0.0031308 ? 12.92 * c : 1.055 * std::pow(c, 1 / 2.4) - 0.055;
}
const M3* matrix(Vision v) { return v == Vision::Protan ? &kProtan : v == Vision::Deutan ? &kDeutan : v == Vision::Tritan ? &kTritan : nullptr; }

std::array<double, 3> apply(const M3& m, double severity, const std::array<double, 3>& lin) {
  std::array<double, 3> o{};
  for (int r = 0; r < 3; ++r) {
    const double s = m[size_t(3 * r)] * lin[0] + m[size_t(3 * r + 1)] * lin[1] + m[size_t(3 * r + 2)] * lin[2];
    o[size_t(r)] = severity * s + (1 - severity) * lin[size_t(r)];
  }
  return o;
}

std::array<double, 3> lab(unsigned rgb) {
  const double r = to_linear(((rgb >> 16) & 255) / 255.0), g = to_linear(((rgb >> 8) & 255) / 255.0), b = to_linear((rgb & 255) / 255.0);
  const double X = (0.4124564 * r + 0.3575761 * g + 0.1804375 * b) / 0.95047;
  const double Y = 0.2126729 * r + 0.7151522 * g + 0.0721750 * b;
  const double Z = (0.0193339 * r + 0.1191920 * g + 0.9503041 * b) / 1.08883;
  auto f = [](double t) { return t > 216.0 / 24389 ? std::cbrt(t) : (24389.0 / 27 * t + 16) / 116; };
  return {116 * f(Y) - 16, 500 * (f(X) - f(Y)), 200 * (f(Y) - f(Z))};
}

}  // namespace

const char* to_string(Vision v) {
  switch (v) {
    case Vision::Protan: return "protanopia";
    case Vision::Deutan: return "deuteranopia";
    case Vision::Tritan: return "tritanopia";
    default: return "normal";
  }
}

unsigned simulate_vision(unsigned rgb, Vision v, double severity) {
  const M3* m = matrix(v);
  if (!m || severity <= 0) return rgb & 0xFFFFFF;
  severity = std::min(severity, 1.0);
  const auto o = apply(*m, severity, {to_linear(((rgb >> 16) & 255) / 255.0), to_linear(((rgb >> 8) & 255) / 255.0), to_linear((rgb & 255) / 255.0)});
  auto q = [](double c) { return unsigned(std::lround(to_srgb(c) * 255)); };
  return (q(o[0]) << 16) | (q(o[1]) << 8) | q(o[2]);
}

void simulate_vision(Image& img, Vision v, double severity) {
  const M3* m = matrix(v);
  if (!m || severity <= 0) return;
  severity = std::min(severity, 1.0);
  std::array<double, 256> lin{};
  for (int i = 0; i < 256; ++i) lin[size_t(i)] = to_linear(i / 255.0);
  std::array<uint8_t, 4096> back{};   // linear → sRGB byte, 12-bit table
  for (int i = 0; i < 4096; ++i) back[size_t(i)] = uint8_t(std::lround(to_srgb(i / 4095.0) * 255));
  for (size_t p = 0; p + 3 < img.rgba.size(); p += 4) {
    const auto o = apply(*m, severity, {lin[img.rgba[p]], lin[img.rgba[p + 1]], lin[img.rgba[p + 2]]});
    for (int c = 0; c < 3; ++c) img.rgba[p + size_t(c)] = back[size_t(std::lround(std::clamp(o[size_t(c)], 0.0, 1.0) * 4095))];
  }
  img.rgba16.clear();   // the 16-bit copy would no longer match
}

double delta_e76(unsigned a, unsigned b) {
  const auto x = lab(a), y = lab(b);
  return std::sqrt((x[0] - y[0]) * (x[0] - y[0]) + (x[1] - y[1]) * (x[1] - y[1]) + (x[2] - y[2]) * (x[2] - y[2]));
}

std::vector<VisionPair> confusable_pairs(const std::vector<NamedPalette>& palettes, double threshold) {
  std::vector<VisionPair> out;
  for (const auto& p : palettes)
    for (Vision v : {Vision::Normal, Vision::Protan, Vision::Deutan, Vision::Tritan}) {
      std::vector<unsigned> c;
      for (unsigned x : p.colours) c.push_back(simulate_vision(x, v));
      for (size_t i = 0; i < c.size(); ++i)
        for (size_t j = i + 1; j < c.size(); ++j) {
          const double de = delta_e76(c[i], c[j]);
          if (de >= threshold) continue;
          // a pair already close in normal vision is reported once, as normal
          if (v != Vision::Normal && delta_e76(p.colours[i], p.colours[j]) < threshold) continue;
          out.push_back({p.name, v, int(i), int(j), de});
        }
    }
  std::sort(out.begin(), out.end(), [](const VisionPair& a, const VisionPair& b) { return a.de < b.de; });
  return out;
}

}  // namespace caps
