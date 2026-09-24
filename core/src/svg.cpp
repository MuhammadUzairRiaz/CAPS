// Vector figure export: the same camera as the raster renderer, drawn back to front.
// Atoms are circles with a radial shading gradient; bonds are two half-lines in the atoms' colours.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <sstream>

#include "caps/elements.hpp"
#include "caps/render.hpp"

namespace caps {
namespace {

std::string hex(unsigned c) {
  char b[8];
  std::snprintf(b, sizeof b, "#%06X", c & 0xFFFFFF);
  return b;
}

unsigned mixu(unsigned a, unsigned b, double t) {
  auto ch = [&](int s) { return unsigned(std::lround(((a >> s) & 255) + (((b >> s) & 255) - double((a >> s) & 255)) * t)); };
  return (ch(16) << 16) | (ch(8) << 8) | ch(0);
}

}  // namespace

std::string render_svg(const System& s, const Camera& cam, const RenderOptions& opt) {
  const size_t n = s.atoms.size();
  const double W = opt.width, H = opt.height;
  const bool transparent = opt.background == Background::Transparent;
  const unsigned bg = background_rgb(opt.background, opt.custom_rgb);
  const bool dark_bg = opt.background == Background::Dark;

  int nmol = 0;
  const auto mol = s.molecules(&nmol);
  std::vector<unsigned> col(n);
  double pmin = 0, pmax = 1;
  if (opt.colour_by == ColourBy::Property && opt.property.size() == n && n) {
    pmin = *std::min_element(opt.property.begin(), opt.property.end());
    pmax = std::max(pmin + 1e-12, *std::max_element(opt.property.begin(), opt.property.end()));
  }
  for (size_t i = 0; i < n; ++i) {
    const Atom& a = s.atoms[i];
    unsigned c = element(a.element).rgb;
    if (opt.colour_by == ColourBy::Molecule) c = a.element == 1 ? mixu(molecule_colour(mol[i]), 0xFFFFFF, 0.55) : molecule_colour(mol[i]);
    else if (opt.colour_by == ColourBy::Type) c = molecule_colour(a.type - 1);
    else if (opt.colour_by == ColourBy::Property && opt.property.size() == n) c = viridis((opt.property[i] - pmin) / (pmax - pmin));
    col[i] = c;
  }
  std::vector<char> show(n, 1);
  if (opt.style == Style::NoHydrogens || opt.style == Style::Backbone)
    for (size_t i = 0; i < n; ++i) show[i] = s.atoms[i].element != 1;

  // Same camera as Renderer::render (orthographic here; perspective is approximated by scale).
  const double cy = std::cos(cam.yaw), sy = std::sin(cam.yaw), cp = std::cos(cam.pitch), sp = std::sin(cam.pitch);
  Vec3 centre{0, 0, 0};
  std::vector<Vec3> corners;
  if (s.cell.valid()) {
    for (int i = 0; i < 2; ++i)
      for (int j = 0; j < 2; ++j)
        for (int k = 0; k < 2; ++k) corners.push_back(s.cell.origin + s.cell.a * i + s.cell.b * j + s.cell.c * k);
    centre = s.cell.origin + (s.cell.a + s.cell.b + s.cell.c) * 0.5;
  } else if (n) {
    Vec3 lo = s.atoms[0].pos, hi = lo;
    for (const auto& a : s.atoms)
      for (int k = 0; k < 3; ++k) { lo[k] = std::min(lo[k], a.pos[k]); hi[k] = std::max(hi[k], a.pos[k]); }
    centre = (lo + hi) * 0.5;
  }
  auto rot = [&](const Vec3& p) {
    const Vec3 d = p - centre;
    const double x = d[0] * cy + d[2] * sy, z = -d[0] * sy + d[2] * cy, y = d[1];
    return Vec3{x + cam.pan_x, y * cp - z * sp + cam.pan_y, y * sp + z * cp};
  };
  double ex = 1e-6, ey = 1e-6;
  for (const auto& c : corners) { const Vec3 r = rot(c); ex = std::max(ex, std::fabs(r[0] - cam.pan_x)); ey = std::max(ey, std::fabs(r[1] - cam.pan_y)); }
  for (size_t i = 0; i < n; ++i)
    if (show[i]) { const Vec3 r = rot(s.atoms[i].pos); ex = std::max(ex, std::fabs(r[0] - cam.pan_x)); ey = std::max(ey, std::fabs(r[1] - cam.pan_y)); }
  const double scale = std::min(W * 0.45 / ex, H * 0.45 / ey) * cam.zoom;
  auto P = [&](const Vec3& p, double& x, double& y, double& z) { const Vec3 r = rot(p); x = W / 2 + r[0] * scale; y = H / 2 - r[1] * scale; z = r[2]; };

  auto radius = [&](size_t i) {
    const double vdw = element(s.atoms[i].element).vdw;
    if (opt.style == Style::SpaceFilling) return vdw;
    if (opt.style == Style::Sticks) return opt.bond_radius;
    if (opt.style == Style::Backbone) return opt.bond_radius * 2.2;
    return std::max(opt.bond_radius * 1.25, vdw * opt.atom_scale);
  };
  const double bond_r = opt.style == Style::Backbone ? opt.bond_radius * 2.2 : opt.bond_radius;
  const unsigned ink = dark_bg ? 0x0B0D0F : 0x141413;

  struct Item { double z; std::string svg; };
  std::vector<Item> items;
  std::map<unsigned, int> grads;
  char buf[512];
  std::vector<double> X(n), Y(n), Z(n);
  for (size_t i = 0; i < n; ++i) P(s.atoms[i].pos, X[i], Y[i], Z[i]);
  if (opt.style != Style::SpaceFilling) {
    const double half_cell = s.cell.valid() ? 0.5 * std::min({norm(s.cell.a), norm(s.cell.b), norm(s.cell.c)}) : 1e300;
    for (const auto& b : s.bonds) {
      if (!show[b.i] || !show[b.j] || norm(s.atoms[b.i].pos - s.atoms[b.j].pos) > half_cell) continue;
      const double mx = (X[b.i] + X[b.j]) / 2, my = (Y[b.i] + Y[b.j]) / 2;
      const double w = 2 * bond_r * scale;
      for (int h = 0; h < 2; ++h) {
        const uint32_t a = h ? b.j : b.i;
        std::snprintf(buf, sizeof buf,
                      "<line x1=\"%.2f\" y1=\"%.2f\" x2=\"%.2f\" y2=\"%.2f\" stroke=\"%s\" stroke-width=\"%.2f\" stroke-linecap=\"round\"/>"
                      "<line x1=\"%.2f\" y1=\"%.2f\" x2=\"%.2f\" y2=\"%.2f\" stroke=\"%s\" stroke-width=\"%.2f\" stroke-linecap=\"round\"/>",
                      X[a], Y[a], mx, my, hex(ink).c_str(), w + 1.2, X[a], Y[a], mx, my, hex(col[a]).c_str(), w);
        items.push_back({(Z[a] + (Z[b.i] + Z[b.j]) / 2) / 2 - bond_r, buf});
      }
    }
  }
  for (size_t i = 0; i < n; ++i) {
    if (!show[i]) continue;
    auto [it, fresh] = grads.try_emplace(col[i], int(grads.size()));
    std::snprintf(buf, sizeof buf, "<circle cx=\"%.2f\" cy=\"%.2f\" r=\"%.2f\" fill=\"url(#g%d)\" stroke=\"%s\" stroke-width=\"0.8\"/>", X[i], Y[i], radius(i) * scale,
                  it->second, hex(ink).c_str());
    items.push_back({Z[i], buf});
  }
  std::stable_sort(items.begin(), items.end(), [](const Item& a, const Item& b) { return a.z < b.z; });

  std::ostringstream o;
  o << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << opt.width << "\" height=\"" << opt.height << "\" viewBox=\"0 0 " << opt.width << " " << opt.height << "\">\n";
  o << "<defs>";
  for (const auto& [c, id] : grads)
    o << "<radialGradient id=\"g" << id << "\" cx=\"0.36\" cy=\"0.32\" r=\"0.75\"><stop offset=\"0\" stop-color=\"" << hex(mixu(c, 0xFFFFFF, 0.35))
      << "\"/><stop offset=\"0.55\" stop-color=\"" << hex(c) << "\"/><stop offset=\"1\" stop-color=\"" << hex(mixu(c, 0x000000, 0.45)) << "\"/></radialGradient>";
  o << "</defs>\n";
  if (!transparent) o << "<rect width=\"100%\" height=\"100%\" fill=\"" << hex(bg) << "\"/>\n";
  if (opt.show_cell && s.cell.valid()) {
    const int E[12][2] = {{0, 1}, {0, 2}, {0, 4}, {1, 3}, {1, 5}, {2, 3}, {2, 6}, {3, 7}, {4, 5}, {4, 6}, {5, 7}, {6, 7}};
    o << "<g stroke=\"" << (dark_bg ? "#A5ABB1" : "#6B7178") << "\" stroke-width=\"1\" fill=\"none\">";
    for (auto& e : E) {
      double x0, y0, z0, x1, y1, z1;
      P(corners[e[0]], x0, y0, z0);
      P(corners[e[1]], x1, y1, z1);
      std::snprintf(buf, sizeof buf, "<line x1=\"%.2f\" y1=\"%.2f\" x2=\"%.2f\" y2=\"%.2f\"/>", x0, y0, x1, y1);
      o << buf;
    }
    o << "</g>\n";
  }
  for (const auto& it : items) o << it.svg << "\n";
  o << "</svg>\n";
  return o.str();
}

}  // namespace caps
