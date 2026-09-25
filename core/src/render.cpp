// CPU ball-and-stick renderer.
// Spheres and bond capsules are rasterised as depth-correct impostors into a supersampled buffer,
// then outlined from depth discontinuities and averaged down with a real alpha channel.
#include "caps/render.hpp"
#include "caps/config.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <limits>

#include "caps/elements.hpp"

namespace caps {
namespace {

struct RGB { float r, g, b; };

RGB rgb(unsigned c) { return {((c >> 16) & 255) / 255.f, ((c >> 8) & 255) / 255.f, (c & 255) / 255.f}; }
RGB mixc(RGB a, RGB b, float t) { return {a.r + (b.r - a.r) * t, a.g + (b.g - a.g) * t, a.b + (b.b - a.b) * t}; }

// Design palette (CAPS "Paper"/dark themes): ten molecule colours, distinct in lightness as well as hue.
// chain colours: design palette v2 (safe for protan, deutan and tritan vision), then two extra
const unsigned kMol[10] = {0xF0A83C, 0x6CC4D8, 0xDE775D, 0x9B7AD5, 0x7DC884, 0xD6AC5C, 0xE9ECEF, 0x2271DB, 0xC77DBA, 0x8FB8A8};
const unsigned kViridis[9] = {0x440154, 0x472D7B, 0x3B528B, 0x2C728E, 0x21918C, 0x28AE80, 0x5EC962, 0xADDC30, 0xFDE725};

struct View {
  double cy, sy, cp, sp;
  Vec3 centre;
  double scale;        // px per Å at the focal plane (supersampled)
  double w, h;         // supersampled size
  double pan_x, pan_y;
  bool persp;
  double dist;         // eye distance for perspective, Å
  Vec3 rot(const Vec3& p) const {
    Vec3 d = p - centre;
    double x = d[0] * cy + d[2] * sy, z = -d[0] * sy + d[2] * cy, y = d[1];
    double y2 = y * cp - z * sp, z2 = y * sp + z * cp;
    return {x + pan_x, y2 + pan_y, z2};
  }
  // Screen x, y (px), view z (Å, larger = nearer), perspective factor.
  void project(const Vec3& p, double& sx, double& sy_, double& z, double& k) const {
    const Vec3 r = rot(p);
    k = persp ? dist / std::max(1e-3, dist - r[2]) : 1.0;
    sx = w / 2 + r[0] * scale * k;
    sy_ = h / 2 - r[1] * scale * k;
    z = r[2];
  }
};

// Camera: centre on the cell (or the atoms), fit the rotated extent into W × H pixels.
View fit_view(const System& s, const Camera& cam, const RenderOptions& opt, const std::vector<char>& show, int W, int H) {
  const size_t n = s.atoms.size();
  View v;
  v.cy = std::cos(cam.yaw); v.sy = std::sin(cam.yaw); v.cp = std::cos(cam.pitch); v.sp = std::sin(cam.pitch);
  std::vector<Vec3> corners;
  if (s.cell.valid()) {
    for (int i = 0; i < 2; ++i)
      for (int j = 0; j < 2; ++j)
        for (int k = 0; k < 2; ++k) corners.push_back(s.cell.origin + s.cell.a * i + s.cell.b * j + s.cell.c * k);
    v.centre = s.cell.origin + (s.cell.a + s.cell.b + s.cell.c) * 0.5;
  } else {
    Vec3 lo{1e300, 1e300, 1e300}, hi{-1e300, -1e300, -1e300};
    for (const auto& a : s.atoms)
      for (int k = 0; k < 3; ++k) { lo[k] = std::min(lo[k], a.pos[k]); hi[k] = std::max(hi[k], a.pos[k]); }
    if (!n) lo = hi = {0, 0, 0};
    v.centre = (lo + hi) * 0.5;
  }
  v.pan_x = cam.pan_x; v.pan_y = cam.pan_y; v.persp = false; v.scale = 1; v.w = W; v.h = H;
  double ex = 1e-6, ey = 1e-6, ez = 1e-6;
  auto extend = [&](const Vec3& p) { const Vec3 r = v.rot(p); ex = std::max(ex, std::fabs(r[0] - v.pan_x)); ey = std::max(ey, std::fabs(r[1] - v.pan_y)); ez = std::max(ez, std::fabs(r[2])); };
  for (const auto& c : corners) extend(c);
  for (size_t i = 0; i < n; ++i) if (show[i]) extend(s.atoms[i].pos);
  // the atoms' own size, and a smallest frame so a small molecule is not blown up to fill the view
  const double pad = opt.style == Style::SpaceFilling ? 2.0 : 1.0;
  ex = std::max(ex + pad, 2.5); ey = std::max(ey + pad, 2.5);
  v.scale = std::min(W * 0.45 / ex, H * 0.45 / ey) * cam.zoom;
  v.persp = cam.perspective;
  v.dist = ez / std::tan(cam.fov_deg * M_PI / 360.0) + ez;
  return v;
}

const Vec3 kLight = [] { Vec3 l{-0.45, 0.6, 1.0}; return l * (1.0 / norm(l)); }();
const Vec3 kHalf = [] { Vec3 h = kLight + Vec3{0, 0, 1}; return h * (1.0 / norm(h)); }();

RGB shade(RGB base, const Vec3& n) {
  const float diff = static_cast<float>(std::max(0.0, dot(n, kLight)));
  const float spec = static_cast<float>(std::pow(std::max(0.0, dot(n, kHalf)), 48.0));
  const float k = 0.38f + 0.62f * diff;
  return {std::min(1.f, base.r * k + 0.28f * spec), std::min(1.f, base.g * k + 0.28f * spec), std::min(1.f, base.b * k + 0.28f * spec)};
}

struct Buffers {
  int w, h;
  std::vector<RGB> col;
  std::vector<float> z;
  std::vector<int32_t> id;
  Buffers(int w_, int h_) : w(w_), h(h_), col(size_t(w_) * h_), z(size_t(w_) * h_, -std::numeric_limits<float>::infinity()), id(size_t(w_) * h_, -1) {}
};

void sphere(Buffers& B, double sx, double sy, double zc, double R, double r_world, RGB base, int32_t id) {
  const int x0 = std::max(0, int(std::floor(sx - R))), x1 = std::min(B.w - 1, int(std::ceil(sx + R)));
  const int y0 = std::max(0, int(std::floor(sy - R))), y1 = std::min(B.h - 1, int(std::ceil(sy + R)));
  const double R2 = R * R;
  for (int y = y0; y <= y1; ++y)
    for (int x = x0; x <= x1; ++x) {
      const double dx = x + 0.5 - sx, dy = y + 0.5 - sy, d2 = dx * dx + dy * dy;
      if (d2 > R2) continue;
      const double nz = std::sqrt(1 - d2 / R2);
      const float z = static_cast<float>(zc + r_world * nz);
      const size_t k = size_t(y) * B.w + x;
      if (z <= B.z[k]) continue;
      B.z[k] = z;
      B.id[k] = id;
      B.col[k] = shade(base, {dx / R, -dy / R, nz});
    }
}

// Capsule from (ax,ay,az) to (bx,by,bz) with screen radius R: a cylinder impostor with rounded ends.
void capsule(Buffers& B, double ax, double ay, double az, double bx, double by, double bz, double R, double r_world, RGB base, int32_t id) {
  const int x0 = std::max(0, int(std::floor(std::min(ax, bx) - R))), x1 = std::min(B.w - 1, int(std::ceil(std::max(ax, bx) + R)));
  const int y0 = std::max(0, int(std::floor(std::min(ay, by) - R))), y1 = std::min(B.h - 1, int(std::ceil(std::max(ay, by) + R)));
  const double ex = bx - ax, ey = by - ay, L2 = ex * ex + ey * ey;
  for (int y = y0; y <= y1; ++y)
    for (int x = x0; x <= x1; ++x) {
      const double px = x + 0.5 - ax, py = y + 0.5 - ay;
      const double t = L2 > 1e-9 ? std::clamp((px * ex + py * ey) / L2, 0.0, 1.0) : 0.0;
      const double dx = px - t * ex, dy = py - t * ey, d2 = dx * dx + dy * dy;
      if (d2 > R * R) continue;
      const double nz = std::sqrt(1 - d2 / (R * R));
      const float z = static_cast<float>(az + t * (bz - az) + r_world * nz);
      const size_t k = size_t(y) * B.w + x;
      if (z <= B.z[k]) continue;
      B.z[k] = z;
      B.id[k] = id;
      B.col[k] = shade(base, {dx / R, -dy / R, nz});
    }
}

void line(Buffers& B, double ax, double ay, double az, double bx, double by, double bz, double width, RGB c) {
  const double len = std::hypot(bx - ax, by - ay);
  const int n = std::max(1, int(len));
  const int hw = std::max(0, int(width / 2));
  for (int s = 0; s <= n; ++s) {
    const double t = double(s) / n;
    const int x = int(ax + t * (bx - ax)), y = int(ay + t * (by - ay));
    const float z = static_cast<float>(az + t * (bz - az)) + 0.05f;   // small bias: edges stay visible on contact
    for (int oy = -hw; oy <= hw; ++oy)
      for (int ox = -hw; ox <= hw; ++ox) {
        const int xx = x + ox, yy = y + oy;
        if (xx < 0 || yy < 0 || xx >= B.w || yy >= B.h) continue;
        const size_t k = size_t(yy) * B.w + xx;
        if (z <= B.z[k]) continue;
        B.z[k] = z;
        B.id[k] = -2;
        B.col[k] = c;
      }
  }
}

}  // namespace

unsigned molecule_colour(int k) {
  if (palette() != Palette::Caps) return category_colour(k);
  return kMol[((k % 10) + 10) % 10];
}

namespace {
std::atomic<int> g_palette{0}, g_threads{0};
}

void set_palette(Palette p) { g_palette = int(p); }
Palette palette() { return Palette(g_palette.load()); }
void set_max_threads(int n) { g_threads = n < 0 ? 0 : n; }
int max_threads() { return g_threads.load(); }

unsigned element_colour(int z) {
  switch (palette()) {
    case Palette::Caps: return element(z).rgb;
    case Palette::OkabeIto:   // Okabe & Ito (2008) colour-universal design; carbon stays neutral
      switch (z) {
        case 1: return 0xE9ECEF;
        case 6: return 0x9AA1A8;
        case 7: return 0x56B4E9;
        case 8: return 0xD55E00;
        case 16: return 0xF0E442;
        case 15: return 0xE69F00;
        case 9: case 17: case 35: case 53: return 0x009E73;
        case 14: return 0x0072B2;
        default: return 0xCC79A7;
      }
    case Palette::Monochrome:   // lightness only, so categories read in greyscale
      switch (z) {
        case 1: return 0xF0F0F0;
        case 6: return 0x7C7C7C;
        case 7: return 0xB4B4B4;
        case 8: return 0x4A4A4A;
        case 16: return 0xD8D8D8;
        default: return 0x9C9C9C;
      }
  }
  return element(z).rgb;
}

unsigned category_colour(int k) {
  static const unsigned okabe[] = {0xE69F00, 0x56B4E9, 0x009E73, 0xF0E442, 0x0072B2, 0xD55E00, 0xCC79A7, 0x9AA1A8};
  static const unsigned mono[] = {0xE8E8E8, 0x6E6E6E, 0xB8B8B8, 0x4A4A4A, 0xD0D0D0, 0x8E8E8E};
  const int i = ((k % 8) + 8) % 8;
  switch (palette()) {
    case Palette::OkabeIto: return okabe[i];
    case Palette::Monochrome: return mono[((k % 6) + 6) % 6];
    default: return kMol[((k % 10) + 10) % 10];
  }
}

unsigned viridis(double t) {
  t = std::clamp(t, 0.0, 1.0) * 8;
  const int i = std::min(7, int(t));
  const RGB c = mixc(rgb(kViridis[i]), rgb(kViridis[i + 1]), float(t - i));
  return (unsigned(c.r * 255 + .5) << 16) | (unsigned(c.g * 255 + .5) << 8) | unsigned(c.b * 255 + .5);
}

unsigned background_rgb(Background b, unsigned custom) {
  switch (b) {
    case Background::Dark: return 0x0F1113;
    case Background::White: return 0xFFFFFF;
    case Background::Transparent: return 0xFFFFFF;
    case Background::Custom: return custom;
  }
  return 0x0F1113;
}

std::vector<float> ambient_accessibility(const System& s, const std::vector<double>& radius, const std::vector<char>& show) {
  const size_t n = s.atoms.size();
  std::vector<float> acc(n, 1.f);
  if (!n) return acc;
  constexpr int kDirs = 32;
  constexpr double kShell = 5.0;
  static const std::vector<Vec3> dirs = [] {
    std::vector<Vec3> d;
    const double ga = M_PI * (3 - std::sqrt(5.0));
    for (int k = 0; k < kDirs; ++k) {
      const double z = 1 - (k + 0.5) * 2.0 / kDirs, r = std::sqrt(1 - z * z);
      d.push_back({r * std::cos(ga * k), r * std::sin(ga * k), z});
    }
    return d;
  }();
  double rmax = 0.7;
  Vec3 lo{1e300, 1e300, 1e300}, hi{-1e300, -1e300, -1e300};
  for (size_t i = 0; i < n; ++i) {
    if (!show[i]) continue;
    rmax = std::max(rmax, radius[i]);
    for (int k = 0; k < 3; ++k) { lo[k] = std::min(lo[k], s.atoms[i].pos[k]); hi[k] = std::max(hi[k], s.atoms[i].pos[k]); }
  }
  if (lo[0] > hi[0]) return acc;
  // a plain bin grid over the drawn atoms (the view shows no periodic images, so neither does the occlusion)
  const double bin = kShell + 2 * rmax;
  int nb[3];
  for (int k = 0; k < 3; ++k) nb[k] = std::clamp(int((hi[k] - lo[k]) / bin) + 1, 1, 256);
  std::vector<std::vector<uint32_t>> bins(size_t(nb[0]) * nb[1] * nb[2]);
  auto bidx = [&](const Vec3& p, int k) { return std::clamp(int((p[k] - lo[k]) / bin), 0, nb[k] - 1); };
  for (size_t i = 0; i < n; ++i)
    if (show[i]) bins[(size_t(bidx(s.atoms[i].pos, 0)) * nb[1] + bidx(s.atoms[i].pos, 1)) * nb[2] + bidx(s.atoms[i].pos, 2)].push_back(uint32_t(i));
  std::vector<uint32_t> near;
  for (size_t i = 0; i < n; ++i) {
    if (!show[i]) continue;
    const Vec3 pi = s.atoms[i].pos;
    const double ri = std::max(radius[i], 0.7);
    near.clear();
    const int b0 = bidx(pi, 0), b1 = bidx(pi, 1), b2 = bidx(pi, 2);
    for (int x = std::max(0, b0 - 1); x <= std::min(nb[0] - 1, b0 + 1); ++x)
      for (int y = std::max(0, b1 - 1); y <= std::min(nb[1] - 1, b1 + 1); ++y)
        for (int z = std::max(0, b2 - 1); z <= std::min(nb[2] - 1, b2 + 1); ++z)
          for (uint32_t j : bins[(size_t(x) * nb[1] + y) * nb[2] + z]) {
            if (j == i) continue;
            const Vec3 d = s.atoms[j].pos - pi;
            const double lim = ri + kShell + std::max(radius[j], 0.7);
            if (dot(d, d) < lim * lim) near.push_back(j);
          }
    int open = 0;
    for (const Vec3& u : dirs) {
      bool hit = false;
      for (uint32_t j : near) {
        const Vec3 d = s.atoms[j].pos - pi;
        const double t = dot(d, u);
        if (t <= 0) continue;
        const double rj = std::max(radius[j], 0.7);
        if (t > ri + kShell + rj) continue;
        const double perp2 = dot(d, d) - t * t;
        if (perp2 < rj * rj) { hit = true; break; }
      }
      if (!hit) ++open;
    }
    acc[i] = float(open) / kDirs;
  }
  return acc;
}

double view_scale(const System& s, const Camera& cam, const RenderOptions& opt) {
  std::vector<char> show(s.atoms.size(), 1);
  if (opt.style == Style::NoHydrogens || opt.style == Style::Backbone)
    for (size_t i = 0; i < s.atoms.size(); ++i) show[i] = s.atoms[i].element != 1;
  return fit_view(s, cam, opt, show, opt.width, opt.height).scale;
}

int Renderer::pick(int x, int y) const {
  if (x < 0 || y < 0 || x >= id_w_ || y >= id_h_) return -1;
  const int v = id_buffer_[size_t(y) * id_w_ + x];
  return v >= 0 ? v : -1;
}

Image Renderer::render(const System& s, const Camera& cam, const RenderOptions& opt) {
  const int ss = std::clamp(opt.supersample, 1, 4);
  const int W = opt.width * ss, H = opt.height * ss;
  Buffers B(W, H);
  const bool dark_bg = opt.background == Background::Dark ||
                       (opt.background == Background::Custom && (((opt.custom_rgb >> 16) & 255) + ((opt.custom_rgb >> 8) & 255) + (opt.custom_rgb & 255)) < 384);
  const RGB bg = rgb(background_rgb(opt.background, opt.custom_rgb));

  // Which atoms are drawn.
  const size_t n = s.atoms.size();
  std::vector<char> show(n, 1);
  if (opt.style == Style::NoHydrogens || opt.style == Style::Backbone)
    for (size_t i = 0; i < n; ++i) show[i] = s.atoms[i].element != 1;

  // Colours.
  int nmol = 0;
  const auto mol = s.molecules(&nmol);
  double pmin = 0, pmax = 1;
  if (opt.colour_by == ColourBy::Property && opt.property.size() == n && n) {
    pmin = *std::min_element(opt.property.begin(), opt.property.end());
    pmax = *std::max_element(opt.property.begin(), opt.property.end());
    if (pmax - pmin < 1e-12) pmax = pmin + 1;
  }
  std::vector<RGB> colour(n);
  const bool overrides = opt.colours.size() == n;
  for (size_t i = 0; i < n; ++i) {
    const Atom& a = s.atoms[i];
    if (overrides && opt.colours[i] != 0xFFFFFFFFu) { colour[i] = rgb(opt.colours[i]); continue; }
    unsigned c = element_colour(a.element);
    switch (opt.colour_by) {
      case ColourBy::Element: break;
      case ColourBy::Molecule: {
        c = molecule_colour(mol[i]);
        if (a.element == 1) { const RGB m = mixc(rgb(c), {1, 1, 1}, 0.55f); colour[i] = m; continue; }
        break;
      }
      case ColourBy::Type: c = molecule_colour(a.type - 1); break;
      case ColourBy::Property:
        if (opt.property.size() == n) c = viridis((opt.property[i] - pmin) / (pmax - pmin));
        break;
    }
    colour[i] = rgb(c);
  }

  // Camera: centre on the cell (or the atoms), fit the rotated extent.
  const View v = fit_view(s, cam, opt, show, W, H);

  // Radii.
  auto radius = [&](size_t i) {
    const double vdw = element(s.atoms[i].element).vdw;
    switch (opt.style) {
      case Style::SpaceFilling: return vdw;
      case Style::Sticks: return opt.bond_radius;
      case Style::Backbone: return opt.bond_radius * 2.2;
      default: return std::max(opt.bond_radius * 1.25, vdw * opt.atom_scale);
    }
  };
  const double bond_r = opt.style == Style::Backbone ? opt.bond_radius * 2.2 : opt.bond_radius;

  std::vector<double> px(n), py(n), pz(n), pk(n);
  double zmin = 1e300, zmax = -1e300;
  for (size_t i = 0; i < n; ++i) {
    v.project(s.atoms[i].pos, px[i], py[i], pz[i], pk[i]);
    if (show[i]) { zmin = std::min(zmin, pz[i]); zmax = std::max(zmax, pz[i]); }
  }

  // Bonds (two half-capsules, each in its atom's colour). Bonds longer than half the cell are not drawn.
  if (opt.style != Style::SpaceFilling) {
    double half_cell = 1e300;
    if (s.cell.valid()) half_cell = 0.5 * std::min({norm(s.cell.a), norm(s.cell.b), norm(s.cell.c)});
    for (const auto& b : s.bonds) {
      if (!show[b.i] || !show[b.j]) continue;
      if (norm(s.atoms[b.i].pos - s.atoms[b.j].pos) > half_cell) continue;
      const double mx = (px[b.i] + px[b.j]) / 2, my = (py[b.i] + py[b.j]) / 2, mz = (pz[b.i] + pz[b.j]) / 2;
      const double R = bond_r * v.scale * (pk[b.i] + pk[b.j]) / 2;
      capsule(B, px[b.i], py[b.i], pz[b.i], mx, my, mz, R, bond_r, colour[b.i], int32_t(b.i));
      capsule(B, mx, my, mz, px[b.j], py[b.j], pz[b.j], R, bond_r, colour[b.j], int32_t(b.j));
    }
  }
  // Atoms.
  for (size_t i = 0; i < n; ++i) {
    if (!show[i]) continue;
    const double r = radius(i);
    sphere(B, px[i], py[i], pz[i], r * v.scale * pk[i], r, colour[i], int32_t(i));
  }
  // Segments: tubes, and arrows whose last part is a stepped cone (not pickable).
  for (const auto& sg : opt.segments) {
    double ax, ay, az, ak, bx, by, bz, bk;
    v.project(sg.a, ax, ay, az, ak);
    v.project(sg.b, bx, by, bz, bk);
    const RGB c = rgb(sg.rgb);
    const double R = std::max(0.6 * ss, sg.radius * v.scale * (ak + bk) / 2);
    if (!sg.arrow) { capsule(B, ax, ay, az, bx, by, bz, R, sg.radius, c, -3); continue; }
    const double len = norm(sg.b - sg.a);
    const double head = std::min(0.35 * len, 6.0 * sg.radius);
    const double f = len > 1e-9 ? 1 - head / len : 0;
    const double hx = ax + (bx - ax) * f, hy = ay + (by - ay) * f, hz = az + (bz - az) * f;
    capsule(B, ax, ay, az, hx, hy, hz, R, sg.radius, c, -3);
    constexpr int steps = 6;
    for (int k = 0; k < steps; ++k) {
      const double t0 = double(k) / steps, t1 = double(k + 1) / steps;
      const double r = 2.4 * R * (1 - t0) + 0.3 * R;
      capsule(B, hx + (bx - hx) * t0, hy + (by - hy) * t0, hz + (bz - hz) * t0, hx + (bx - hx) * t1, hy + (by - hy) * t1, hz + (bz - hz) * t1, r,
              sg.radius * 2.4 * (1 - t0), c, -3);
    }
  }

  // Ambient occlusion: each atom's pixels (and its half-bonds) darkened by how little open sky the atom sees.
  if (opt.ambient_occlusion && n) {
    double key = double(n) * 1e-3 + double(opt.style) * 7;
    for (size_t i = 0; i < n; i += std::max<size_t>(1, n / 512)) key += s.atoms[i].pos[0] * 1.3 + s.atoms[i].pos[1] * 1.7 + s.atoms[i].pos[2] * 2.9;
    if (ao_.size() != n || ao_key_ != key) {
      std::vector<double> rad(n);
      for (size_t i = 0; i < n; ++i) rad[i] = radius(i);
      ao_ = ambient_accessibility(s, rad, show);
      ao_key_ = key;
    }
    for (size_t k = 0; k < B.col.size(); ++k) {
      const int32_t id = B.id[k];
      if (id < 0 || size_t(id) >= n) continue;
      const float f = 0.3f + 0.7f * ao_[size_t(id)];
      B.col[k] = {B.col[k].r * f, B.col[k].g * f, B.col[k].b * f};
    }
  }
  // Cell edges.
  if (opt.show_cell && s.cell.valid()) {
    const RGB ec = dark_bg ? rgb(0xA5ABB1) : rgb(0x6B7178);
    const int E[12][2] = {{0, 1}, {0, 2}, {0, 4}, {1, 3}, {1, 5}, {2, 3}, {2, 6}, {3, 7}, {4, 5}, {4, 6}, {5, 7}, {6, 7}};
    double cx[8], cyy[8], cz[8], ck[8];
    for (int k = 0; k < 8; ++k)
      v.project(s.cell.origin + s.cell.a * (k >> 2) + s.cell.b * ((k >> 1) & 1) + s.cell.c * (k & 1), cx[k], cyy[k], cz[k], ck[k]);
    for (auto& e : E) line(B, cx[e[0]], cyy[e[0]], cz[e[0]], cx[e[1]], cyy[e[1]], cz[e[1]], std::max(1.0, 1.1 * ss), ec);
  }

  // Depth cue: fade far atoms toward the background; with no background, darken them instead.
  if (opt.depth_cue && zmax > zmin) {
    for (size_t k = 0; k < B.col.size(); ++k) {
      if (B.id[k] < 0) continue;
      const float t = static_cast<float>(0.35 * (zmax - B.z[k]) / (zmax - zmin));
      B.col[k] = opt.background == Background::Transparent ? mixc(B.col[k], {0, 0, 0}, t * 0.6f) : mixc(B.col[k], bg, t);
    }
  }

  // Outlines where depth jumps or coverage ends.
  if (opt.outlines) {
    const float jump = 0.6f;
    const int rpx = std::max(1, ss);
    std::vector<char> edge(B.col.size(), 0);
    for (int y = 0; y < H; ++y)
      for (int x = 0; x < W; ++x) {
        const size_t k = size_t(y) * W + x;
        if (B.id[k] < 0) continue;
        bool e = false;
        for (int d = 1; d <= rpx && !e; ++d) {
          const int nx[4] = {x + d, x - d, x, x}, ny[4] = {y, y, y + d, y - d};
          for (int q = 0; q < 4 && !e; ++q) {
            if (nx[q] < 0 || ny[q] < 0 || nx[q] >= W || ny[q] >= H) { e = true; break; }
            const size_t m = size_t(ny[q]) * W + nx[q];
            if (B.id[m] == -1 || (B.z[m] < B.z[k] - jump && B.id[m] != B.id[k])) e = true;
          }
        }
        edge[k] = e;
      }
    const RGB ink = dark_bg ? RGB{0.04f, 0.045f, 0.05f} : RGB{0.08f, 0.08f, 0.08f};
    for (size_t k = 0; k < edge.size(); ++k)
      if (edge[k] && B.id[k] >= 0) B.col[k] = mixc(B.col[k], ink, 0.8f);
  }

  // Selection rings.
  for (int hi : opt.highlight) {
    if (hi < 0 || size_t(hi) >= n || !show[size_t(hi)]) continue;
    const size_t i = size_t(hi);
    const double R = radius(i) * v.scale * pk[i] + 3.0 * ss;
    const RGB sel = rgb(0x6CC4D8);
    for (int y = int(py[i] - R - 3 * ss); y <= int(py[i] + R + 3 * ss); ++y)
      for (int x = int(px[i] - R - 3 * ss); x <= int(px[i] + R + 3 * ss); ++x) {
        if (x < 0 || y < 0 || x >= W || y >= H) continue;
        const double d = std::hypot(x + 0.5 - px[i], y + 0.5 - py[i]);
        if (std::fabs(d - R) <= 1.2 * ss) {
          const size_t k = size_t(y) * W + x;
          B.col[k] = sel; B.id[k] = int32_t(i); B.z[k] = 1e30f;
        }
      }
  }

  // Keyboard-focus ring: accent, a little outside the selection ring, drawn in two arcs so it reads as focus.
  if (opt.focus >= 0 && size_t(opt.focus) < n && show[size_t(opt.focus)]) {
    const size_t i = size_t(opt.focus);
    const double R = radius(i) * v.scale * pk[i] + 6.5 * ss;
    const RGB acc = rgb(0xF5A524);
    for (int y = int(py[i] - R - 3 * ss); y <= int(py[i] + R + 3 * ss); ++y)
      for (int x = int(px[i] - R - 3 * ss); x <= int(px[i] + R + 3 * ss); ++x) {
        if (x < 0 || y < 0 || x >= W || y >= H) continue;
        const double dx = x + 0.5 - px[i], dy = y + 0.5 - py[i];
        const double d = std::hypot(dx, dy);
        const double a = std::atan2(dy, dx);
        if (std::fabs(d - R) <= 1.6 * ss && std::fabs(std::sin(a)) > 0.26) {   // gaps left and right
          const size_t k = size_t(y) * W + x;
          B.col[k] = acc; B.id[k] = int32_t(i); B.z[k] = 1e30f;
        }
      }
  }

  // Downsample with coverage as alpha.
  Image img;
  img.width = opt.width;
  img.height = opt.height;
  img.rgba.resize(size_t(img.width) * img.height * 4);
  id_w_ = img.width; id_h_ = img.height;
  id_buffer_.assign(size_t(id_w_) * id_h_, -1);
  const float inv = 1.f / (ss * ss);
  for (int y = 0; y < img.height; ++y)
    for (int x = 0; x < img.width; ++x) {
      float r = 0, g = 0, b = 0, a = 0;
      for (int j = 0; j < ss; ++j)
        for (int i = 0; i < ss; ++i) {
          const size_t k = size_t(y * ss + j) * W + (x * ss + i);
          if (B.id[k] == -1) continue;
          r += B.col[k].r; g += B.col[k].g; b += B.col[k].b; a += 1;
        }
      r *= inv; g *= inv; b *= inv; a *= inv;   // premultiplied
      float R_, G_, Bc, A_;
      if (opt.background == Background::Transparent) {
        A_ = a;
        R_ = a > 0 ? r / a : 0; G_ = a > 0 ? g / a : 0; Bc = a > 0 ? b / a : 0;
      } else {
        R_ = r + (1 - a) * bg.r; G_ = g + (1 - a) * bg.g; Bc = b + (1 - a) * bg.b; A_ = 1;
      }
      uint8_t* p = &img.rgba[(size_t(y) * img.width + x) * 4];
      p[0] = uint8_t(std::clamp(R_, 0.f, 1.f) * 255 + .5f);
      p[1] = uint8_t(std::clamp(G_, 0.f, 1.f) * 255 + .5f);
      p[2] = uint8_t(std::clamp(Bc, 0.f, 1.f) * 255 + .5f);
      p[3] = uint8_t(std::clamp(A_, 0.f, 1.f) * 255 + .5f);
      const size_t c = size_t(y * ss + ss / 2) * W + (x * ss + ss / 2);
      id_buffer_[size_t(y) * id_w_ + x] = B.id[c];
    }
  return img;
}

}  // namespace caps
