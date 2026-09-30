#include "caps/raytrace.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <thread>
#include <vector>

namespace caps {
namespace {

struct V {
  float x = 0, y = 0, z = 0;
};
inline V operator+(V a, V b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
inline V operator-(V a, V b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
inline V operator*(V a, float s) { return {a.x * s, a.y * s, a.z * s}; }
inline float dotv(V a, V b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline V crossv(V a, V b) { return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x}; }
inline V unit(V a) {
  const float l = std::sqrt(dotv(a, a));
  return l > 1e-20f ? a * (1.f / l) : V{0, 0, 1};
}
inline float comp(V a, int k) { return k == 0 ? a.x : k == 1 ? a.y : a.z; }

// A small fast generator per sample stream (splitmix64).
struct Rng {
  uint64_t s;
  explicit Rng(uint64_t seed) : s(seed) {}
  uint64_t next() {
    uint64_t z = (s += 0x9E3779B97F4A7C15ull);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ull;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBull;
    return z ^ (z >> 31);
  }
  float uni() { return float(next() >> 40) * (1.f / 16777216.f); }
};

struct Col {
  float r = 0, g = 0, b = 0;
};
Col unpack(uint32_t c) { return {((c >> 16) & 255) / 255.f, ((c >> 8) & 255) / 255.f, (c & 255) / 255.f}; }

// Everything in view space (x right, y up, z toward the viewer).
enum Kind : uint8_t { Sphere, Capsule, Triangle };
struct Prim {
  Kind kind;
  uint32_t index;   // into its own array
};
struct SphereP {
  V c;
  float r;
  Col col;
  float opacity;
};
struct CapsuleP {
  V a, b;
  float r;
  Col col;
  float opacity;
};
struct TriP {
  V p[3], n[3];
  Col col[3];
  float opacity;
};

struct Box {
  V lo{1e30f, 1e30f, 1e30f}, hi{-1e30f, -1e30f, -1e30f};
  void grow(V p) {
    lo = {std::min(lo.x, p.x), std::min(lo.y, p.y), std::min(lo.z, p.z)};
    hi = {std::max(hi.x, p.x), std::max(hi.y, p.y), std::max(hi.z, p.z)};
  }
  void grow(const Box& b) { grow(b.lo); grow(b.hi); }
  V centre() const { return (lo + hi) * 0.5f; }
};

struct Node {
  Box box;
  uint32_t left = 0, count = 0;   // leaf: first primitive and count; inner: left child (right = left + 1), count 0
};

struct Hit {
  float t = std::numeric_limits<float>::infinity();
  V n;
  Col col;
  float opacity = 1;
};

class Tracer {
 public:
  std::vector<SphereP> spheres;
  std::vector<CapsuleP> capsules;
  std::vector<TriP> tris;
  std::vector<Prim> prims;
  std::vector<Box> pbox;
  std::vector<Node> nodes;

  void build() {
    for (const auto& s : spheres) {
      Box b;
      b.grow(s.c - V{s.r, s.r, s.r});
      b.grow(s.c + V{s.r, s.r, s.r});
      prims.push_back({Sphere, uint32_t(&s - spheres.data())});
      pbox.push_back(b);
    }
    for (const auto& c : capsules) {
      Box b;
      const V r{c.r, c.r, c.r};
      b.grow(c.a - r); b.grow(c.a + r); b.grow(c.b - r); b.grow(c.b + r);
      prims.push_back({Capsule, uint32_t(&c - capsules.data())});
      pbox.push_back(b);
    }
    for (const auto& t : tris) {
      Box b;
      for (const auto& p : t.p) b.grow(p);
      prims.push_back({Triangle, uint32_t(&t - tris.data())});
      pbox.push_back(b);
    }
    order.resize(prims.size());
    for (size_t i = 0; i < order.size(); ++i) order[i] = uint32_t(i);
    nodes.clear();
    nodes.reserve(prims.size() * 2 + 1);
    nodes.push_back({});
    if (!prims.empty()) split(0, 0, uint32_t(prims.size()));
    // primitives in leaf order
    std::vector<Prim> p2(prims.size());
    std::vector<Box> b2(prims.size());
    for (size_t i = 0; i < order.size(); ++i) p2[i] = prims[order[i]], b2[i] = pbox[order[i]];
    prims.swap(p2);
    pbox.swap(b2);
  }

  // The nearest surface along the ray within (tmin, tmax).
  bool nearest(V o, V d, float tmin, float tmax, Hit& h) const {
    if (nodes.empty() || prims.empty()) return false;
    const V inv{1.f / (std::fabs(d.x) > 1e-12f ? d.x : 1e-12f), 1.f / (std::fabs(d.y) > 1e-12f ? d.y : 1e-12f), 1.f / (std::fabs(d.z) > 1e-12f ? d.z : 1e-12f)};
    uint32_t stack[96];
    int sp = 0;
    stack[sp++] = 0;
    bool any = false;
    h.t = tmax;
    while (sp) {
      const Node& nd = nodes[stack[--sp]];
      if (!slab(nd.box, o, inv, tmin, h.t)) continue;
      if (nd.count) {
        for (uint32_t k = nd.left; k < nd.left + nd.count; ++k) any |= intersect(prims[k], o, d, tmin, h);
      } else if (sp + 2 <= 96) {
        stack[sp++] = nd.left + 1;
        stack[sp++] = nd.left;
      }
    }
    return any;
  }

  // How much light passes along a shadow or occlusion ray (1 clear, 0 blocked; translucent surfaces let some through).
  float transmit(V o, V d, float tmin, float tmax) const {
    float pass = 1;
    float t0 = tmin;
    for (int layer = 0; layer < 8 && pass > 0.01f; ++layer) {
      Hit h;
      if (!nearest(o, d, t0, tmax, h)) break;
      pass *= 1 - h.opacity;
      t0 = h.t + 1e-3f;
    }
    return pass;
  }

 private:
  std::vector<uint32_t> order;

  static bool slab(const Box& b, V o, V inv, float tmin, float tmax) {
    float t0 = tmin, t1 = tmax;
    for (int k = 0; k < 3; ++k) {
      float ta = (comp(b.lo, k) - comp(o, k)) * comp(inv, k), tb = (comp(b.hi, k) - comp(o, k)) * comp(inv, k);
      if (ta > tb) std::swap(ta, tb);
      t0 = std::max(t0, ta);
      t1 = std::min(t1, tb);
      if (t0 > t1) return false;
    }
    return true;
  }

  void split(uint32_t node, uint32_t first, uint32_t count) {
    Box b, cb;
    for (uint32_t i = first; i < first + count; ++i) { b.grow(pbox[order[i]]); cb.grow(pbox[order[i]].centre()); }
    nodes[node].box = b;
    if (count <= 4) { nodes[node].left = first, nodes[node].count = count; return; }
    // the longest axis of the centres, split at the median
    const V ext = cb.hi - cb.lo;
    const int ax = ext.x > ext.y && ext.x > ext.z ? 0 : ext.y > ext.z ? 1 : 2;
    const uint32_t mid = first + count / 2;
    std::nth_element(order.begin() + first, order.begin() + mid, order.begin() + first + count,
                     [&](uint32_t a, uint32_t c) { return comp(pbox[a].centre(), ax) < comp(pbox[c].centre(), ax); });
    const uint32_t left = uint32_t(nodes.size());
    nodes.push_back({});
    nodes.push_back({});
    nodes[node].left = left, nodes[node].count = 0;
    split(left, first, mid - first);
    split(left + 1, mid, first + count - mid);
  }

  static bool hit_sphere(V c, float r, V o, V d, float tmin, float& t) {
    const V oc = o - c;
    const float b = dotv(oc, d), cc = dotv(oc, oc) - r * r, disc = b * b - cc;
    if (disc < 0) return false;
    const float s = std::sqrt(disc);
    float x = -b - s;
    if (x <= tmin) x = -b + s;
    if (x <= tmin || x >= t) return false;
    t = x;
    return true;
  }

  bool intersect(const Prim& p, V o, V d, float tmin, Hit& h) const {
    switch (p.kind) {
      case Sphere: {
        const auto& s = spheres[p.index];
        float t = h.t;
        if (!hit_sphere(s.c, s.r, o, d, tmin, t)) return false;
        h.t = t;
        h.n = unit(o + d * t - s.c);
        h.col = s.col, h.opacity = s.opacity;
        return true;
      }
      case Capsule: {
        const auto& c = capsules[p.index];
        // the cylinder's side, then its two end spheres
        bool got = false;
        const V ba = c.b - c.a, oa = o - c.a;
        const float baba = dotv(ba, ba), bard = dotv(ba, d), baoa = dotv(ba, oa);
        const float k2 = baba - bard * bard, k1 = baba * dotv(oa, d) - baoa * bard, k0 = baba * dotv(oa, oa) - baoa * baoa - c.r * c.r * baba;
        if (std::fabs(k2) > 1e-12f) {
          const float disc = k1 * k1 - k2 * k0;
          if (disc >= 0) {
            const float s = std::sqrt(disc);
            for (float t : {(-k1 - s) / k2, (-k1 + s) / k2}) {
              if (t <= tmin || t >= h.t) continue;
              const float y = baoa + t * bard;
              if (y <= 0 || y >= baba) continue;
              const V q = oa + d * t - ba * (y / baba);
              h.t = t, h.n = unit(q), h.col = c.col, h.opacity = c.opacity;
              got = true;
              break;
            }
          }
        }
        for (const V& e : {c.a, c.b}) {
          float t = h.t;
          if (hit_sphere(e, c.r, o, d, tmin, t)) { h.t = t, h.n = unit(o + d * t - e), h.col = c.col, h.opacity = c.opacity; got = true; }
        }
        return got;
      }
      case Triangle: {
        const auto& tr = tris[p.index];
        const V e1 = tr.p[1] - tr.p[0], e2 = tr.p[2] - tr.p[0], pv = crossv(d, e2);
        const float det = dotv(e1, pv);
        if (std::fabs(det) < 1e-12f) return false;
        const float inv = 1.f / det;
        const V tv = o - tr.p[0];
        const float u = dotv(tv, pv) * inv;
        if (u < 0 || u > 1) return false;
        const V qv = crossv(tv, e1);
        const float v = dotv(d, qv) * inv;
        if (v < 0 || u + v > 1) return false;
        const float t = dotv(e2, qv) * inv;
        if (t <= tmin || t >= h.t) return false;
        const float w = 1 - u - v;
        h.t = t;
        h.n = unit(tr.n[0] * w + tr.n[1] * u + tr.n[2] * v);
        if (dotv(h.n, h.n) < 1e-12f) h.n = unit(crossv(e1, e2));
        h.col = {tr.col[0].r * w + tr.col[1].r * u + tr.col[2].r * v, tr.col[0].g * w + tr.col[1].g * u + tr.col[2].g * v,
                 tr.col[0].b * w + tr.col[1].b * u + tr.col[2].b * v};
        h.opacity = tr.opacity;
        return true;
      }
    }
    return false;
  }
};

// Cosine-weighted direction about n.
V cosine_dir(V n, Rng& g) {
  const float u1 = g.uni(), u2 = g.uni();
  const float r = std::sqrt(u1), phi = 6.2831853f * u2;
  const V t = unit(std::fabs(n.x) > 0.9f ? crossv(n, V{0, 1, 0}) : crossv(n, V{1, 0, 0}));
  const V b = crossv(n, t);
  return unit(t * (r * std::cos(phi)) + b * (r * std::sin(phi)) + n * std::sqrt(std::max(0.f, 1 - u1)));
}

}  // namespace

Image raytrace(const Scene& sc, const ViewFit& f, int w, int h, const TraceOptions& o) {
  if (w <= 0 || h <= 0) throw std::invalid_argument("raytrace: the image needs a size");
  Tracer T;
  const float cy = float(f.cos_yaw), sy = float(f.sin_yaw), cp = float(f.cos_pitch), spn = float(f.sin_pitch);
  const V centre{float(f.centre[0]), float(f.centre[1]), float(f.centre[2])};
  auto view = [&](float x, float y, float z) {   // world → view (the views' rotation and pan)
    const V dd = V{x, y, z} - centre;
    const float xr = dd.x * cy + dd.z * sy, zr = -dd.x * sy + dd.z * cy;
    return V{xr + float(f.pan_x), dd.y * cp - zr * spn + float(f.pan_y), dd.y * spn + zr * cp};
  };
  auto viewd = [&](float x, float y, float z) {   // world → view for a direction
    const float xr = x * cy + z * sy, zr = -x * sy + z * cy;
    return V{xr, y * cp - zr * spn, y * spn + zr * cp};
  };
  // the scene's colours carry transparency (spheres, capsules: top byte, 0 opaque) or opacity (triangles)
  for (size_t i = 0; i + 3 < sc.spheres.size() && i / 4 < sc.sphere_rgb.size(); i += 4) {
    const uint32_t c = sc.sphere_rgb[i / 4];
    if (sc.spheres[i + 3] <= 0) continue;
    T.spheres.push_back({view(sc.spheres[i], sc.spheres[i + 1], sc.spheres[i + 2]), sc.spheres[i + 3], unpack(c), 1 - ((c >> 24) & 255) / 255.f});
  }
  for (size_t i = 0; i + 6 < sc.capsules.size() && i / 7 < sc.capsule_rgb.size(); i += 7) {
    const uint32_t c = sc.capsule_rgb[i / 7];
    T.capsules.push_back({view(sc.capsules[i], sc.capsules[i + 1], sc.capsules[i + 2]), view(sc.capsules[i + 3], sc.capsules[i + 4], sc.capsules[i + 5]),
                          sc.capsules[i + 6], unpack(c), 1 - ((c >> 24) & 255) / 255.f});
  }
  const float scale = float(f.scale);
  for (size_t i = 0; i + 5 < sc.lines.size() && i / 6 < sc.line_rgb.size(); i += 6) {   // lines: thin tubes of their pixel width
    const float wpx = i / 6 < sc.line_width.size() ? sc.line_width[i / 6] : 1.f;
    T.capsules.push_back({view(sc.lines[i], sc.lines[i + 1], sc.lines[i + 2]), view(sc.lines[i + 3], sc.lines[i + 4], sc.lines[i + 5]),
                          std::max(0.5f * wpx / std::max(scale, 1e-6f), 0.02f), unpack(sc.line_rgb[i / 6]), 1.f});
  }
  for (size_t t = 0; 9 * t + 8 < sc.tri_xyz.size() && 3 * t + 2 < sc.tri_rgb.size(); ++t) {
    TriP tp;
    for (int c = 0; c < 3; ++c) {
      const size_t k = 9 * t + 3 * size_t(c);
      tp.p[c] = view(sc.tri_xyz[k], sc.tri_xyz[k + 1], sc.tri_xyz[k + 2]);
      tp.n[c] = k + 2 < sc.tri_normal.size() ? viewd(sc.tri_normal[k], sc.tri_normal[k + 1], sc.tri_normal[k + 2]) : V{0, 0, 1};
      tp.col[c] = unpack(sc.tri_rgb[3 * t + size_t(c)]);
    }
    tp.opacity = ((sc.tri_rgb[3 * t] >> 24) & 255) / 255.f;
    T.tris.push_back(tp);
  }
  T.build();

  // the scene's depth extent: where rays start (orthographic) and how far they go
  float zhi = 1, zlo = -1;
  if (!T.nodes.empty() && !T.prims.empty()) zhi = T.nodes[0].box.hi.z + 1, zlo = T.nodes[0].box.lo.z - 1;
  const float span = std::max(zhi - zlo, 1.f);
  const bool persp = f.perspective;
  const float dist = float(f.dist);
  const Col bg = unpack(sc.background);
  const V L = unit(V{-0.45f, 0.6f, 1.0f}), Hh = unit(L + V{0, 0, 1});
  const V Lu = unit(crossv(L, std::fabs(L.y) < 0.9f ? V{0, 1, 0} : V{1, 0, 0})), Lv = crossv(L, Lu);
  const Col ink = sc.dark ? Col{0.04f, 0.045f, 0.05f} : Col{0.08f, 0.08f, 0.08f};
  const int spp = std::max(1, o.samples);
  const float ao_d = float(std::max(0.1, o.ao_distance));
  const float scale_px = scale;   // px per Å at the focal plane (the fit was made for w × h)

  Image img;
  img.width = w, img.height = h;
  img.rgba.assign(size_t(w) * h * 4, 0);
  img.rgba16.assign(size_t(w) * h * 4, 0);

  // one shaded sample: the colour and coverage along a primary ray (front to back through translucent surfaces)
  auto trace = [&](V ro, V rd, Rng& g, float& A) {
    Col acc{0, 0, 0};
    A = 0;
    float t0 = 1e-4f;
    const float tfar = persp ? dist + span * 4 : span * 4 + 10;
    for (int layer = 0; layer < std::max(1, o.max_layers) && A < 0.995f; ++layer) {
      Hit hit;
      if (!T.nearest(ro, rd, t0, tfar, hit)) break;
      const V p = ro + rd * hit.t;
      V n = hit.n;
      if (dotv(n, rd) > 0) n = n * -1.f;   // two-sided (surfaces)
      // the key light: a small disc about L, blocked (partly, by translucent ones) by what lies toward it
      float light = 1;
      if (o.shadows) {
        const float a = 6.2831853f * g.uni(), rr = float(o.light_size) * std::sqrt(g.uni());
        const V ldir = unit(L + Lu * (rr * std::cos(a)) + Lv * (rr * std::sin(a)));
        light = dotv(n, ldir) > 0 ? T.transmit(p + n * 1e-3f, ldir, 1e-3f, 1e6f) : 0.f;
      }
      float ao = 1;
      if (o.ambient_occlusion) {
        const V od = cosine_dir(n, g);
        const float pass = T.transmit(p + n * 1e-3f, od, 1e-3f, ao_d);
        ao = 1 - float(o.ao_strength) * (1 - pass);
      }
      const float diff = std::max(0.f, dotv(n, L)), spec = std::pow(std::max(0.f, dotv(n, Hh)), 40.f);
      Col c{std::min(1.f, hit.col.r * (0.32f * ao + 0.72f * diff * light) + 0.5f * spec * light),
            std::min(1.f, hit.col.g * (0.32f * ao + 0.72f * diff * light) + 0.5f * spec * light),
            std::min(1.f, hit.col.b * (0.32f * ao + 0.72f * diff * light) + 0.5f * spec * light)};
      if (o.outlines) {   // a thin rim where the surface turns away (about a pixel at this size)
        const float facing = std::fabs(dotv(n, rd));
        const float rim = std::clamp((0.22f - facing) / 0.22f, 0.f, 1.f) * 0.55f;
        c = {c.r + (ink.r - c.r) * rim, c.g + (ink.g - c.g) * rim, c.b + (ink.b - c.b) * rim};
      }
      const float a = hit.opacity * (1 - A);
      acc = {acc.r + c.r * a, acc.g + c.g * a, acc.b + c.b * a};
      A += a;
      t0 = hit.t + 1e-3f;
    }
    if (!sc.transparent) {   // the background behind whatever is left uncovered
      const float rest = 1 - A;
      acc = {acc.r + bg.r * rest, acc.g + bg.g * rest, acc.b + bg.b * rest};
      A = 1;
    }
    return acc;
  };

  // the primary ray through a point of the image (pixels), with the lens for depth of field
  auto primary = [&](float px, float py, Rng& g, V& ro, V& rd) {
    const float X = (px - w * 0.5f) / scale_px, Y = -(py - h * 0.5f) / scale_px;
    const float zf = float(o.focus);
    if (persp) {
      const V eye{0, 0, dist};
      rd = unit(V{X, Y, 0} - eye);
      ro = eye;
    } else {
      ro = V{X, Y, zhi + 1};
      rd = V{0, 0, -1};
    }
    if (o.aperture > 0) {   // a thin lens: every ray through the same point of the focal plane
      const float tf = (zf - ro.z) / (std::fabs(rd.z) > 1e-6f ? rd.z : -1e-6f);
      const V focal = ro + rd * tf;
      const float a = 6.2831853f * g.uni(), rr = float(o.aperture) * std::sqrt(g.uni());
      ro = ro + V{rr * std::cos(a), rr * std::sin(a), 0};
      rd = unit(focal - ro);
    }
  };

  // tiles over the threads; each tile is finished with all its samples
  const int tile = 32;
  const int tx = (w + tile - 1) / tile, ty = (h + tile - 1) / tile, ntiles = tx * ty;
  std::atomic<int> next{0}, done{0};
  std::atomic<bool> cancel{false};
  int nthreads = o.threads > 0 ? o.threads : int(std::max(1u, std::thread::hardware_concurrency()));
  nthreads = std::min(nthreads, ntiles);
  const int side = int(std::ceil(std::sqrt(double(spp))));   // stratified: a side × side grid per pixel
  auto work = [&](int thread) {
    (void)thread;
    for (;;) {
      const int k = next.fetch_add(1);
      if (k >= ntiles || cancel.load()) return;
      const int x0 = (k % tx) * tile, y0 = (k / tx) * tile;
      for (int y = y0; y < std::min(h, y0 + tile); ++y)
        for (int x = x0; x < std::min(w, x0 + tile); ++x) {
          Rng g(o.seed * 0x9E3779B97F4A7C15ull ^ (uint64_t(y) * 73856093ull + uint64_t(x) * 19349663ull + 1));
          float R = 0, G = 0, B = 0, AA = 0;
          for (int s = 0; s < spp; ++s) {
            const int sxi = s % side, syi = (s / side) % side;
            const float jx = (sxi + g.uni()) / side, jy = (syi + g.uni()) / side;
            V ro, rd;
            primary(x + jx, y + jy, g, ro, rd);
            float a;
            const Col c = trace(ro, rd, g, a);
            R += c.r, G += c.g, B += c.b, AA += a;
          }
          const float inv = 1.f / spp;
          float a = AA * inv;
          // straight alpha: the colour of what covers the pixel
          const float r = a > 1e-6f ? R * inv / a : 0, gg = a > 1e-6f ? G * inv / a : 0, b = a > 1e-6f ? B * inv / a : 0;
          const size_t q = (size_t(y) * w + x) * 4;
          const float v[4] = {r, gg, b, a};
          for (int c = 0; c < 4; ++c) {
            const float cl = std::clamp(v[c], 0.f, 1.f);
            img.rgba16[q + size_t(c)] = uint16_t(std::lround(cl * 65535));
            img.rgba[q + size_t(c)] = uint8_t(std::lround(cl * 255));
          }
        }
      const int d = done.fetch_add(1) + 1;
      if (o.progress && thread == 0 && !o.progress(double(d) / ntiles)) cancel = true;
    }
  };
  std::vector<std::thread> pool;
  for (int t = 1; t < nthreads; ++t) pool.emplace_back(work, t);
  work(0);
  for (auto& th : pool) th.join();
  if (cancel.load()) throw std::runtime_error("raytrace: cancelled");
  if (o.progress) o.progress(1.0);
  return img;
}

}  // namespace caps
