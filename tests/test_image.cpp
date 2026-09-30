#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>

#include "caps/io.hpp"
#include "caps/json.hpp"
#include "caps/raytrace.hpp"
#include "caps/scene_export.hpp"
#include <cstring>
#include "caps/render.hpp"

using namespace caps;

namespace {

std::string slurp(const std::string& path) {
  std::ifstream f(path, std::ios::binary);
  return std::string((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
}
uint32_t be32(const std::string& b, size_t k) {
  return uint32_t(uint8_t(b[k])) << 24 | uint32_t(uint8_t(b[k + 1])) << 16 | uint32_t(uint8_t(b[k + 2])) << 8 | uint8_t(b[k + 3]);
}
Image gradient(int w, int h, int shift = 0) {
  Image img;
  img.width = w;
  img.height = h;
  img.rgba.resize(size_t(w) * h * 4);
  img.rgba16.resize(img.rgba.size());
  for (int y = 0; y < h; ++y)
    for (int x = 0; x < w; ++x)
      for (int c = 0; c < 4; ++c) {
        const size_t k = (size_t(y) * w + x) * 4 + c;
        img.rgba16[k] = c == 3 ? 65535 : uint16_t((x + shift) * 65535 / (w - 1 + shift));
        img.rgba[k] = uint8_t(img.rgba16[k] >> 8);
      }
  return img;
}

}  // namespace

TEST(Image, SixteenBitPngWithDpiSrgbAndText) {
  const auto path = (std::filesystem::temp_directory_path() / "caps_image16.png").string();
  PngOptions o;
  o.sixteen = true;
  o.dpi = 600;
  o.text = {{"caps:provenance", "{\"schema\":\"caps-image/1.0\",\"note\":\"ρ = 1.05 g/cm³\"}"}};
  write_png(gradient(64, 8), path, o);
  const std::string b = slurp(path);
  ASSERT_GT(b.size(), 33u);
  EXPECT_EQ(b.substr(12, 4), "IHDR");
  EXPECT_EQ(uint8_t(b[24]), 16);   // bit depth
  EXPECT_NE(b.find("sRGB"), std::string::npos);
  const size_t phys = b.find("pHYs");
  ASSERT_NE(phys, std::string::npos);
  EXPECT_EQ(be32(b, phys + 4), 23622u);   // 600 dpi in pixels per metre
  const auto text = read_png_text(path);
  ASSERT_EQ(text.size(), 1u);
  EXPECT_EQ(text[0].first, "caps:provenance");
  EXPECT_EQ(text[0].second, o.text[0].second);
}

TEST(Image, AnimatedPngCountsItsFrames) {
  const auto path = (std::filesystem::temp_directory_path() / "caps_movie.png").string();
  {
    ApngWriter w(path, 12, 0);
    for (int k = 0; k < 5; ++k) w.add(gradient(32, 16, k));
    w.close();
    EXPECT_EQ(w.frames(), 5);
  }
  const std::string b = slurp(path);
  const size_t actl = b.find("acTL");
  ASSERT_NE(actl, std::string::npos);
  EXPECT_LT(actl, b.find("IDAT"));          // before the image data, as APNG requires
  EXPECT_EQ(be32(b, actl + 4), 5u);         // frame count patched on close
  size_t fctl = 0, fdat = 0;
  for (size_t at = b.find("fcTL"); at != std::string::npos; at = b.find("fcTL", at + 4)) ++fctl;
  for (size_t at = b.find("fdAT"); at != std::string::npos; at = b.find("fdAT", at + 4)) ++fdat;
  EXPECT_EQ(fctl, 5u);
  EXPECT_EQ(fdat, 4u);                      // the first frame is the IDAT
  EXPECT_EQ(b.substr(b.size() - 8, 4), "IEND");
}

TEST(Image, DeepRenderKeepsSixteenBits) {
  System s;
  Atom a;
  a.element = 6;
  s.atoms.push_back(a);
  RenderOptions o;
  o.width = 40;
  o.height = 40;
  o.supersample = 4;
  o.deep = true;
  o.background = Background::White;
  Renderer r;
  const Image img = r.render(s, Camera{}, o);
  ASSERT_EQ(img.rgba16.size(), img.rgba.size());
  int finer = 0;   // antialiased edge pixels whose 16-bit value is not an 8-bit value widened
  for (size_t k = 0; k < img.rgba16.size(); ++k) {
    EXPECT_NEAR(img.rgba16[k] / 257.0, img.rgba[k], 0.51);
    finer += img.rgba16[k] % 257 != 0;
  }
  EXPECT_GT(finer, 0);
}

// Picking by the view ray (no image) names the same atom as the rendered id buffer at each visible atom's centre.
TEST(Image, RayPickAgreesWithTheRenderedIds) {
  const caps::System s = caps::open_file(std::string(CAPS_SAMPLES) + "/ps_melt.data").frame(0);
  caps::Renderer r;
  caps::Camera cam;
  cam.yaw = 0.7, cam.pitch = 0.35, cam.zoom = 1.3;
  for (bool persp : {false, true}) {
    cam.perspective = persp;
    caps::RenderOptions opt;
    opt.width = 640, opt.height = 480, opt.supersample = 2;
    r.render(s, cam, opt);
    const auto p = r.project(s, cam, opt);
    size_t tested = 0, same = 0;
    for (size_t i = 0; i < s.atoms.size(); ++i) {
      if (p[3 * i + 2] < 0.5f) continue;   // hidden behind others
      const int x = int(p[3 * i]), y = int(p[3 * i + 1]);
      if (x < 0 || y < 0 || x >= 640 || y >= 480) continue;
      const int buf = r.pick(x, y);
      if (buf < 0) continue;
      ++tested;
      same += caps::Renderer::pick_ray(s, cam, opt, x + 0.75, y + 0.75) == buf;   // where the id buffer samples (the supersample right of centre)
    }
    EXPECT_GT(tested, 200u);
    EXPECT_GE(double(same) / double(tested), 0.97) << (persp ? "perspective" : "orthographic") << ": " << same << " of " << tested;
  }
  caps::RenderOptions opt;
  opt.width = 640, opt.height = 480;
  EXPECT_EQ(caps::Renderer::pick_ray(s, cam, opt, 2, 2), -1);   // a corner of empty space
}

// The ray tracer sees the same geometry through the same camera as the CPU renderer: their silhouettes (coverage on a
// transparent background) agree, and the traced image keeps 16 bits.
TEST(Image, RayTracedSilhouetteMatchesTheRender) {
  const caps::System s = caps::open_file(std::string(CAPS_SAMPLES) + "/ps_melt.data").frame(0);
  caps::Renderer r;
  caps::Camera cam;
  cam.yaw = 0.6, cam.pitch = 0.3;
  for (bool persp : {false, true}) {
    cam.perspective = persp;
    caps::RenderOptions opt;
    opt.width = 320, opt.height = 240, opt.supersample = 1, opt.background = caps::Background::Transparent, opt.show_cell = false;
    const caps::Image cpu = r.render(s, cam, opt);
    const caps::Scene sc = r.scene(s, opt);
    const caps::ViewFit fit = caps::view_fit(s, cam, opt);
    caps::TraceOptions to;
    to.samples = 4;
    const caps::Image rt = caps::raytrace(sc, fit, 320, 240, to);
    ASSERT_EQ(rt.rgba16.size(), size_t(320 * 240 * 4));
    size_t both = 0, either = 0;
    for (size_t p = 0; p < size_t(320 * 240); ++p) {
      const bool a = cpu.rgba[4 * p + 3] > 127, b = rt.rgba[4 * p + 3] > 127;
      both += a && b, either += a || b;
    }
    ASSERT_GT(either, 2000u);
    EXPECT_GT(double(both) / double(either), 0.88) << (persp ? "perspective" : "orthographic");   // edges of 1–2 px sticks sample differently
  }
}

// The scene written for other tools: POV-Ray (a sphere per atom, a cylinder per half-bond), glTF binary (a valid
// container whose JSON counts match its binary chunk) and OBJ (vertices with colours, triangles).
TEST(Image, SceneExportsPovGltfObj) {
  const caps::System s = caps::open_file(std::string(CAPS_SAMPLES) + "/ps_melt.data").frame(0);
  caps::Renderer r;
  caps::Camera cam;
  caps::RenderOptions opt;
  opt.width = 800, opt.height = 600, opt.supersample = 1, opt.show_cell = false;
  const caps::Scene sc = r.scene(s, opt);
  const caps::ViewFit fit = caps::view_fit(s, cam, opt);
  const auto dir = std::filesystem::temp_directory_path();
  const auto pov = caps::write_povray(sc, fit, 800, 600, (dir / "caps_scene.pov").string());
  EXPECT_EQ(pov.spheres, sc.sphere_rgb.size());
  EXPECT_EQ(pov.cylinders, sc.capsule_rgb.size());
  std::ifstream pf(dir / "caps_scene.pov");
  const std::string ptxt((std::istreambuf_iterator<char>(pf)), std::istreambuf_iterator<char>());
  EXPECT_NE(ptxt.find("camera {"), std::string::npos);
  EXPECT_NE(ptxt.find("light_source"), std::string::npos);
  const auto glb = caps::write_gltf(sc, (dir / "caps_scene.glb").string());
  std::ifstream gf(dir / "caps_scene.glb", std::ios::binary);
  const std::string g((std::istreambuf_iterator<char>(gf)), std::istreambuf_iterator<char>());
  ASSERT_GT(g.size(), 28u);
  uint32_t magic, version, total, jlen, jtype;
  std::memcpy(&magic, g.data(), 4); std::memcpy(&version, g.data() + 4, 4); std::memcpy(&total, g.data() + 8, 4);
  std::memcpy(&jlen, g.data() + 12, 4); std::memcpy(&jtype, g.data() + 16, 4);
  EXPECT_EQ(magic, 0x46546C67u); EXPECT_EQ(version, 2u); EXPECT_EQ(total, g.size()); EXPECT_EQ(jtype, 0x4E4F534Au);
  const caps::Json j = caps::Json::parse(g.substr(20, jlen));
  uint32_t blen;
  std::memcpy(&blen, g.data() + 20 + jlen, 4);
  EXPECT_EQ(size_t(j["buffers"][0]["byteLength"].number()), size_t(blen));
  EXPECT_EQ(size_t(j["accessors"][3]["count"].number()), glb.triangles * 3);
  EXPECT_EQ(glb.spheres, sc.sphere_rgb.size());
  const auto obj = caps::write_obj(sc, (dir / "caps_scene.obj").string());
  std::ifstream of(dir / "caps_scene.obj");
  size_t nv = 0, nf = 0;
  for (std::string line; std::getline(of, line);) { nv += line.rfind("v ", 0) == 0; nf += line.rfind("f ", 0) == 0; }
  EXPECT_EQ(nf, obj.triangles);
  EXPECT_GT(nv, sc.sphere_rgb.size() * 12);
}
