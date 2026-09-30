#include "caps/scene_export.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <map>
#include <stdexcept>
#include <vector>

namespace caps {
namespace {

std::array<float, 3> rgbf(uint32_t c) { return {((c >> 16) & 255) / 255.f, ((c >> 8) & 255) / 255.f, (c & 255) / 255.f}; }

// One coloured triangle soup: positions, normals, colours (RGBA 0 … 1) per vertex.
struct Soup {
  std::vector<float> pos, nrm, col;
  std::vector<uint32_t> idx;
  uint32_t vertex(const std::array<double, 3>& p, const std::array<double, 3>& n, const std::array<float, 4>& c) {
    pos.insert(pos.end(), {float(p[0]), float(p[1]), float(p[2])});
    nrm.insert(nrm.end(), {float(n[0]), float(n[1]), float(n[2])});
    col.insert(col.end(), c.begin(), c.end());
    return uint32_t(pos.size() / 3 - 1);
  }
};

// A unit icosphere, subdivided `level` times.
void icosphere(int level, std::vector<std::array<double, 3>>& v, std::vector<std::array<uint32_t, 3>>& f) {
  const double t = (1 + std::sqrt(5.0)) / 2;
  v = {{-1, t, 0}, {1, t, 0}, {-1, -t, 0}, {1, -t, 0}, {0, -1, t}, {0, 1, t}, {0, -1, -t}, {0, 1, -t}, {t, 0, -1}, {t, 0, 1}, {-t, 0, -1}, {-t, 0, 1}};
  for (auto& p : v) { const double l = std::sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]); p = {p[0] / l, p[1] / l, p[2] / l}; }
  f = {{0, 11, 5}, {0, 5, 1}, {0, 1, 7}, {0, 7, 10}, {0, 10, 11}, {1, 5, 9}, {5, 11, 4}, {11, 10, 2}, {10, 7, 6}, {7, 1, 8},
       {3, 9, 4}, {3, 4, 2}, {3, 2, 6}, {3, 6, 8}, {3, 8, 9}, {4, 9, 5}, {2, 4, 11}, {6, 2, 10}, {8, 6, 7}, {9, 8, 1}};
  for (int l = 0; l < level; ++l) {
    std::map<std::pair<uint32_t, uint32_t>, uint32_t> mid;
    auto m = [&](uint32_t a, uint32_t b) {
      const auto key = std::minmax(a, b);
      if (auto it = mid.find(key); it != mid.end()) return it->second;
      std::array<double, 3> p{(v[a][0] + v[b][0]) / 2, (v[a][1] + v[b][1]) / 2, (v[a][2] + v[b][2]) / 2};
      const double n = std::sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
      v.push_back({p[0] / n, p[1] / n, p[2] / n});
      return mid[key] = uint32_t(v.size() - 1);
    };
    std::vector<std::array<uint32_t, 3>> g;
    for (const auto& t3 : f) {
      const uint32_t a = m(t3[0], t3[1]), b = m(t3[1], t3[2]), c = m(t3[2], t3[0]);
      g.push_back({t3[0], a, c}); g.push_back({t3[1], b, a}); g.push_back({t3[2], c, b}); g.push_back({a, b, c});
    }
    f.swap(g);
  }
}

// The scene tessellated: spheres as icospheres, capsules and lines as open tubes, mesh triangles as they are.
Soup tessellate(const Scene& sc, int detail, SceneExportReport& rep) {
  const size_t n = sc.sphere_rgb.size();
  if (detail < 0) detail = n > 200000 ? 0 : n > 20000 ? 1 : 2;
  detail = std::clamp(detail, 0, 3);
  std::vector<std::array<double, 3>> sv;
  std::vector<std::array<uint32_t, 3>> sf;
  icosphere(detail, sv, sf);
  const int sides = detail >= 2 ? 12 : detail == 1 ? 8 : 6;
  Soup s;
  auto rgba = [](uint32_t c, bool transparency) {
    const auto f = rgbf(c);
    const float a = ((c >> 24) & 255) / 255.f;
    return std::array<float, 4>{f[0], f[1], f[2], transparency ? 1 - a : a};
  };
  for (size_t i = 0; i < n && 4 * i + 3 < sc.spheres.size(); ++i) {
    const double r = sc.spheres[4 * i + 3];
    if (r <= 0) continue;
    const std::array<double, 3> c{sc.spheres[4 * i], sc.spheres[4 * i + 1], sc.spheres[4 * i + 2]};
    const auto col = rgba(sc.sphere_rgb[i], true);
    const uint32_t base = uint32_t(s.pos.size() / 3);
    for (const auto& p : sv) s.vertex({c[0] + r * p[0], c[1] + r * p[1], c[2] + r * p[2]}, p, col);
    for (const auto& t : sf) s.idx.insert(s.idx.end(), {base + t[0], base + t[1], base + t[2]});
    ++rep.spheres;
  }
  auto tube = [&](std::array<double, 3> a, std::array<double, 3> b, double r, const std::array<float, 4>& col) {
    std::array<double, 3> d{b[0] - a[0], b[1] - a[1], b[2] - a[2]};
    const double L = std::sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
    if (L < 1e-9 || r <= 0) return;
    d = {d[0] / L, d[1] / L, d[2] / L};
    std::array<double, 3> u = std::fabs(d[0]) < 0.9 ? std::array<double, 3>{1, 0, 0} : std::array<double, 3>{0, 1, 0};
    std::array<double, 3> e{d[1] * u[2] - d[2] * u[1], d[2] * u[0] - d[0] * u[2], d[0] * u[1] - d[1] * u[0]};
    const double el = std::sqrt(e[0] * e[0] + e[1] * e[1] + e[2] * e[2]);
    e = {e[0] / el, e[1] / el, e[2] / el};
    const std::array<double, 3> g{d[1] * e[2] - d[2] * e[1], d[2] * e[0] - d[0] * e[2], d[0] * e[1] - d[1] * e[0]};
    const uint32_t base = uint32_t(s.pos.size() / 3);
    for (int k = 0; k < sides; ++k) {
      const double ang = 2 * M_PI * k / sides, cs = std::cos(ang), sn = std::sin(ang);
      const std::array<double, 3> nn{e[0] * cs + g[0] * sn, e[1] * cs + g[1] * sn, e[2] * cs + g[2] * sn};
      s.vertex({a[0] + r * nn[0], a[1] + r * nn[1], a[2] + r * nn[2]}, nn, col);
      s.vertex({b[0] + r * nn[0], b[1] + r * nn[1], b[2] + r * nn[2]}, nn, col);
    }
    for (int k = 0; k < sides; ++k) {
      const uint32_t a0 = base + 2 * k, b0 = a0 + 1, a1 = base + 2 * ((k + 1) % sides), b1 = a1 + 1;
      s.idx.insert(s.idx.end(), {a0, a1, b0, b0, a1, b1});
    }
    ++rep.cylinders;
  };
  for (size_t i = 0; i < sc.capsule_rgb.size() && 7 * i + 6 < sc.capsules.size(); ++i) {
    const float* c = &sc.capsules[7 * i];
    tube({c[0], c[1], c[2]}, {c[3], c[4], c[5]}, c[6], rgba(sc.capsule_rgb[i], true));
  }
  for (size_t i = 0; i < sc.line_rgb.size() && 6 * i + 5 < sc.lines.size(); ++i) {
    const float* c = &sc.lines[6 * i];
    tube({c[0], c[1], c[2]}, {c[3], c[4], c[5]}, 0.04, rgba(sc.line_rgb[i] & 0xFFFFFFu, true));
  }
  for (size_t t = 0; 9 * t + 8 < sc.tri_xyz.size() && 3 * t + 2 < sc.tri_rgb.size(); ++t) {
    const uint32_t base = uint32_t(s.pos.size() / 3);
    for (int c = 0; c < 3; ++c) {
      const size_t k = 9 * t + 3 * size_t(c);
      const std::array<double, 3> nn = k + 2 < sc.tri_normal.size() ? std::array<double, 3>{sc.tri_normal[k], sc.tri_normal[k + 1], sc.tri_normal[k + 2]}
                                                                     : std::array<double, 3>{0, 0, 1};
      s.vertex({sc.tri_xyz[k], sc.tri_xyz[k + 1], sc.tri_xyz[k + 2]}, nn, rgba(sc.tri_rgb[3 * t + size_t(c)], false));
    }
    s.idx.insert(s.idx.end(), {base, base + 1, base + 2});
  }
  rep.triangles = s.idx.size() / 3;
  return s;
}

std::ofstream open_out(const std::string& path, bool binary) {
  std::ofstream f(path, binary ? std::ios::binary : std::ios::out);
  if (!f) throw std::runtime_error("cannot write " + path);
  return f;
}

}  // namespace

SceneExportReport write_povray(const Scene& sc, const ViewFit& fit, int w, int h, const std::string& path) {
  SceneExportReport rep;
  auto f = open_out(path, false);
  char b[512];
  // the view's axes in world space: the rotation's rows (world → view is R), so right, up and toward-the-viewer are
  // R's rows; the view centre sits at centre − Rᵀ·pan
  const double cy = fit.cos_yaw, sy = fit.sin_yaw, cp = fit.cos_pitch, spn = fit.sin_pitch;
  const std::array<double, 3> right{cy, 0, sy}, up{sy * spn, cp, -cy * spn}, back{-sy * cp, spn, cy * cp};
  const std::array<double, 3> focus{fit.centre[0] - right[0] * fit.pan_x - up[0] * fit.pan_y, fit.centre[1] - right[1] * fit.pan_x - up[1] * fit.pan_y,
                                    fit.centre[2] - right[2] * fit.pan_x - up[2] * fit.pan_y};
  const double hw = w / (2 * fit.scale), hh = h / (2 * fit.scale);   // half the view at the focal plane (Å)
  f << "// POV-Ray 3.7 scene written by CAPS: the view's camera, key light and colours. Render: povray +W" << w << " +H" << h << " +A +Q11 " << path << "\n";
  f << "#version 3.7;\nglobal_settings { assumed_gamma 1.0 }\n";
  const auto bg = rgbf(sc.background);
  if (sc.transparent) f << "background { color rgbt <1, 1, 1, 1> }\n";
  else { std::snprintf(b, sizeof b, "background { color srgb <%.4f, %.4f, %.4f> }\n", bg[0], bg[1], bg[2]); f << b; }
  // explicit camera vectors (no look_at): a ray is direction + x·right + y·up, x to the image's right, so right is the
  // view's right and up its up — the image is the view's, unmirrored
  if (fit.perspective) {
    const double d = fit.dist;
    std::snprintf(b, sizeof b,
                  "camera {\n  perspective\n  location <%.5f, %.5f, %.5f>\n  direction <%.5f, %.5f, %.5f>\n  right <%.5f, %.5f, %.5f>\n  up <%.5f, %.5f, %.5f>\n}\n",
                  focus[0] + back[0] * d, focus[1] + back[1] * d, focus[2] + back[2] * d, -back[0] * d, -back[1] * d, -back[2] * d,
                  right[0] * 2 * hw, right[1] * 2 * hw, right[2] * 2 * hw, up[0] * 2 * hh, up[1] * 2 * hh, up[2] * 2 * hh);
  } else {
    const double d = 1000;
    std::snprintf(b, sizeof b,
                  "camera {\n  orthographic\n  location <%.5f, %.5f, %.5f>\n  direction <%.5f, %.5f, %.5f>\n  right <%.5f, %.5f, %.5f>\n  up <%.5f, %.5f, %.5f>\n}\n",
                  focus[0] + back[0] * d, focus[1] + back[1] * d, focus[2] + back[2] * d, -back[0], -back[1], -back[2],
                  right[0] * 2 * hw, right[1] * 2 * hw, right[2] * 2 * hw, up[0] * 2 * hh, up[1] * 2 * hh, up[2] * 2 * hh);
  }
  f << b;
  // the key light of the views (upper left, in front), soft; a dim fill from the camera
  const double lx = -0.45, ly = 0.6, lz = 1.0, ln = std::sqrt(lx * lx + ly * ly + lz * lz);
  const double far = 4 * std::max(hw, hh) + 200;
  const std::array<double, 3> L{(right[0] * lx + up[0] * ly + back[0] * lz) / ln, (right[1] * lx + up[1] * ly + back[1] * lz) / ln,
                                (right[2] * lx + up[2] * ly + back[2] * lz) / ln};
  std::snprintf(b, sizeof b, "light_source { <%.3f, %.3f, %.3f> color rgb 1.0 area_light <%.2f, 0, 0>, <0, 0, %.2f>, 5, 5 adaptive 1 jitter }\n",
                focus[0] + L[0] * far, focus[1] + L[1] * far, focus[2] + L[2] * far, far * 0.08, far * 0.08);
  f << b;
  std::snprintf(b, sizeof b, "light_source { <%.3f, %.3f, %.3f> color rgb 0.35 shadowless }\n", focus[0] + back[0] * far, focus[1] + back[1] * far,
                focus[2] + back[2] * far);
  f << b;
  f << "#declare CapsFinish = finish { ambient 0.18 diffuse 0.75 specular 0.45 roughness 0.02 }\n";
  auto pig = [&](uint32_t c, bool transparency) {
    const auto q = rgbf(c);
    const double t = transparency ? ((c >> 24) & 255) / 255.0 : 1 - ((c >> 24) & 255) / 255.0;
    char p[128];
    std::snprintf(p, sizeof p, "pigment { color srgbt <%.4f, %.4f, %.4f, %.3f> } finish { CapsFinish }", q[0], q[1], q[2], t);
    return std::string(p);
  };
  f << "union {\n";
  for (size_t i = 0; i < sc.sphere_rgb.size() && 4 * i + 3 < sc.spheres.size(); ++i) {
    const float* s = &sc.spheres[4 * i];
    if (s[3] <= 0) continue;
    std::snprintf(b, sizeof b, "sphere { <%.4f, %.4f, %.4f>, %.4f ", s[0], s[1], s[2], s[3]);
    f << b << pig(sc.sphere_rgb[i], true) << " }\n";
    ++rep.spheres;
  }
  for (size_t i = 0; i < sc.capsule_rgb.size() && 7 * i + 6 < sc.capsules.size(); ++i) {
    const float* c = &sc.capsules[7 * i];
    if (std::fabs(c[0] - c[3]) + std::fabs(c[1] - c[4]) + std::fabs(c[2] - c[5]) < 1e-6) continue;
    std::snprintf(b, sizeof b, "cylinder { <%.4f, %.4f, %.4f>, <%.4f, %.4f, %.4f>, %.4f ", c[0], c[1], c[2], c[3], c[4], c[5], c[6]);
    f << b << pig(sc.capsule_rgb[i], true) << " }\n";
    ++rep.cylinders;
  }
  for (size_t i = 0; i < sc.line_rgb.size() && 6 * i + 5 < sc.lines.size(); ++i) {
    const float* c = &sc.lines[6 * i];
    if (std::fabs(c[0] - c[3]) + std::fabs(c[1] - c[4]) + std::fabs(c[2] - c[5]) < 1e-6) continue;
    std::snprintf(b, sizeof b, "cylinder { <%.4f, %.4f, %.4f>, <%.4f, %.4f, %.4f>, %.4f ", c[0], c[1], c[2], c[3], c[4], c[5], 0.5 * (i < sc.line_width.size() ? sc.line_width[i] : 1.f) / fit.scale);
    f << b << pig(sc.line_rgb[i] & 0xFFFFFFu, true) << " }\n";
    ++rep.cylinders;
  }
  f << "}\n";
  if (!sc.tri_rgb.empty()) {
    const size_t nt = sc.tri_rgb.size() / 3;
    f << "mesh2 {\n  vertex_vectors { " << 3 * nt;
    for (size_t k = 0; k < 3 * nt; ++k) { std::snprintf(b, sizeof b, ", <%.4f, %.4f, %.4f>", sc.tri_xyz[3 * k], sc.tri_xyz[3 * k + 1], sc.tri_xyz[3 * k + 2]); f << b; }
    f << " }\n  normal_vectors { " << 3 * nt;
    for (size_t k = 0; k < 3 * nt; ++k) {
      const bool has = 3 * k + 2 < sc.tri_normal.size();
      std::snprintf(b, sizeof b, ", <%.4f, %.4f, %.4f>", has ? sc.tri_normal[3 * k] : 0.f, has ? sc.tri_normal[3 * k + 1] : 0.f, has ? sc.tri_normal[3 * k + 2] : 1.f);
      f << b;
    }
    f << " }\n  texture_list { " << 3 * nt;
    for (size_t k = 0; k < 3 * nt; ++k) f << ", texture { " << pig(sc.tri_rgb[k], false) << " }";
    f << " }\n  face_indices { " << nt;
    for (size_t t = 0; t < nt; ++t) f << ", <" << 3 * t << ", " << 3 * t + 1 << ", " << 3 * t + 2 << ">, " << 3 * t << ", " << 3 * t + 1 << ", " << 3 * t + 2;
    f << " }\n}\n";
    rep.triangles = nt;
  }
  if (!f) throw std::runtime_error("could not write " + path);
  return rep;
}

SceneExportReport write_obj(const Scene& sc, const std::string& path, int detail) {
  SceneExportReport rep;
  const Soup s = tessellate(sc, detail, rep);
  // colours as materials (a companion .mtl beside it, one material per colour), the portable way every OBJ reader knows
  std::string mtl = path;
  const auto dot = mtl.find_last_of('.');
  mtl = (dot == std::string::npos ? mtl : mtl.substr(0, dot)) + ".mtl";
  const std::string mtl_name = mtl.substr(mtl.find_last_of("/\\") == std::string::npos ? 0 : mtl.find_last_of("/\\") + 1);
  std::map<uint32_t, int> mats;   // packed 8-bit RGB + opacity → material index
  auto key = [&](uint32_t v) {
    const float* c = &s.col[4 * v];
    return (uint32_t(std::lround(c[0] * 255)) << 24) | (uint32_t(std::lround(c[1] * 255)) << 16) | (uint32_t(std::lround(c[2] * 255)) << 8) |
           uint32_t(std::lround(c[3] * 255));
  };
  std::vector<int> face_mat(s.idx.size() / 3);
  for (size_t t = 0; t < face_mat.size(); ++t) {
    const uint32_t k = key(s.idx[3 * t]);
    auto it = mats.find(k);
    if (it == mats.end()) it = mats.emplace(k, int(mats.size())).first;
    face_mat[t] = it->second;
  }
  {
    auto m = open_out(mtl, false);
    m << "# materials written by CAPS (one per colour of the view)\n";
    for (const auto& [k, idx] : mats) {
      char b[200];
      std::snprintf(b, sizeof b, "newmtl caps%d\nKd %.4f %.4f %.4f\nKa %.4f %.4f %.4f\nKs 0.4 0.4 0.4\nNs 60\nd %.3f\n\n", idx, ((k >> 24) & 255) / 255.0,
                    ((k >> 16) & 255) / 255.0, ((k >> 8) & 255) / 255.0, 0.2 * ((k >> 24) & 255) / 255.0, 0.2 * ((k >> 16) & 255) / 255.0,
                    0.2 * ((k >> 8) & 255) / 255.0, (k & 255) / 255.0);
      m << b;
    }
  }
  auto f = open_out(path, false);
  f << "# Wavefront OBJ written by CAPS: the view's atoms, bonds and surfaces as triangles; colours in " << mtl_name << "\n";
  f << "mtllib " << mtl_name << "\n";
  char b[160];
  for (size_t k = 0; k < s.pos.size() / 3; ++k) { std::snprintf(b, sizeof b, "v %.4f %.4f %.4f\n", s.pos[3 * k], s.pos[3 * k + 1], s.pos[3 * k + 2]); f << b; }
  for (size_t k = 0; k < s.nrm.size() / 3; ++k) { std::snprintf(b, sizeof b, "vn %.4f %.4f %.4f\n", s.nrm[3 * k], s.nrm[3 * k + 1], s.nrm[3 * k + 2]); f << b; }
  // faces grouped by material
  std::vector<std::vector<size_t>> by(mats.size());
  for (size_t t = 0; t < face_mat.size(); ++t) by[size_t(face_mat[t])].push_back(t);
  for (size_t mi = 0; mi < by.size(); ++mi) {
    if (by[mi].empty()) continue;
    f << "usemtl caps" << mi << '\n';
    for (size_t t : by[mi]) {
      const uint32_t a = s.idx[3 * t] + 1, c = s.idx[3 * t + 1] + 1, d = s.idx[3 * t + 2] + 1;
      f << "f " << a << "//" << a << ' ' << c << "//" << c << ' ' << d << "//" << d << '\n';
    }
  }
  if (!f) throw std::runtime_error("could not write " + path);
  return rep;
}

SceneExportReport write_gltf(const Scene& sc, const std::string& path, int detail) {
  SceneExportReport rep;
  const Soup s = tessellate(sc, detail, rep);
  const size_t nv = s.pos.size() / 3;
  // one binary buffer: positions, normals, colours (RGBA float), indices (uint32)
  std::vector<uint8_t> bin;
  auto put = [&](const void* p, size_t bytes) { const auto* c = static_cast<const uint8_t*>(p); bin.insert(bin.end(), c, c + bytes); while (bin.size() % 4) bin.push_back(0); };
  const size_t o_pos = bin.size(); put(s.pos.data(), s.pos.size() * 4);
  const size_t o_nrm = bin.size(); put(s.nrm.data(), s.nrm.size() * 4);
  const size_t o_col = bin.size(); put(s.col.data(), s.col.size() * 4);
  const size_t o_idx = bin.size(); put(s.idx.data(), s.idx.size() * 4);
  float lo[3] = {1e30f, 1e30f, 1e30f}, hi[3] = {-1e30f, -1e30f, -1e30f};
  for (size_t k = 0; k < nv; ++k)
    for (int c = 0; c < 3; ++c) lo[c] = std::min(lo[c], s.pos[3 * k + size_t(c)]), hi[c] = std::max(hi[c], s.pos[3 * k + size_t(c)]);
  if (!nv) lo[0] = lo[1] = lo[2] = hi[0] = hi[1] = hi[2] = 0;
  bool alpha = false;
  for (size_t k = 3; k < s.col.size(); k += 4) if (s.col[k] < 0.999f) { alpha = true; break; }
  char j[4096];
  std::snprintf(j, sizeof j,
    "{\"asset\":{\"version\":\"2.0\",\"generator\":\"CAPS\"},\"scene\":0,\"scenes\":[{\"nodes\":[0]}],"
    "\"nodes\":[{\"mesh\":0,\"name\":\"CAPS structure\"}],"
    "\"materials\":[{\"name\":\"CAPS\",\"pbrMetallicRoughness\":{\"baseColorFactor\":[1,1,1,1],\"metallicFactor\":0.0,\"roughnessFactor\":0.45},\"alphaMode\":\"%s\",\"doubleSided\":true}],"
    "\"meshes\":[{\"primitives\":[{\"attributes\":{\"POSITION\":0,\"NORMAL\":1,\"COLOR_0\":2},\"indices\":3,\"material\":0}]}],"
    "\"buffers\":[{\"byteLength\":%zu}],"
    "\"bufferViews\":[{\"buffer\":0,\"byteOffset\":%zu,\"byteLength\":%zu,\"target\":34962},{\"buffer\":0,\"byteOffset\":%zu,\"byteLength\":%zu,\"target\":34962},"
    "{\"buffer\":0,\"byteOffset\":%zu,\"byteLength\":%zu,\"target\":34962},{\"buffer\":0,\"byteOffset\":%zu,\"byteLength\":%zu,\"target\":34963}],"
    "\"accessors\":[{\"bufferView\":0,\"componentType\":5126,\"count\":%zu,\"type\":\"VEC3\",\"min\":[%g,%g,%g],\"max\":[%g,%g,%g]},"
    "{\"bufferView\":1,\"componentType\":5126,\"count\":%zu,\"type\":\"VEC3\"},{\"bufferView\":2,\"componentType\":5126,\"count\":%zu,\"type\":\"VEC4\"},"
    "{\"bufferView\":3,\"componentType\":5125,\"count\":%zu,\"type\":\"SCALAR\"}]}",
    alpha ? "BLEND" : "OPAQUE", bin.size(), o_pos, s.pos.size() * 4, o_nrm, s.nrm.size() * 4, o_col, s.col.size() * 4, o_idx, s.idx.size() * 4,
    nv, lo[0], lo[1], lo[2], hi[0], hi[1], hi[2], nv, nv, s.idx.size());
  std::string json(j);
  while (json.size() % 4) json += ' ';
  auto f = open_out(path, true);
  auto u32 = [&](uint32_t v) { f.write(reinterpret_cast<const char*>(&v), 4); };
  u32(0x46546C67); u32(2); u32(uint32_t(12 + 8 + json.size() + 8 + bin.size()));
  u32(uint32_t(json.size())); u32(0x4E4F534A); f.write(json.data(), std::streamsize(json.size()));
  u32(uint32_t(bin.size())); u32(0x004E4942); f.write(reinterpret_cast<const char*>(bin.data()), std::streamsize(bin.size()));
  if (!f) throw std::runtime_error("could not write " + path);
  return rep;
}

}  // namespace caps
