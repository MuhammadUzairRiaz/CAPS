#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>

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
