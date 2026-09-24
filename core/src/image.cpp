// PNG encoder (RGBA 8-bit, zlib via the system library). Keeps the alpha channel for transparent figures.
#include <zlib.h>

#include <cstring>
#include <fstream>
#include <stdexcept>

#include "caps/render.hpp"

namespace caps {
namespace {

void put32(std::vector<uint8_t>& v, uint32_t x) {
  v.push_back(uint8_t(x >> 24)); v.push_back(uint8_t(x >> 16)); v.push_back(uint8_t(x >> 8)); v.push_back(uint8_t(x));
}

void chunk(std::vector<uint8_t>& out, const char* type, const std::vector<uint8_t>& data) {
  put32(out, uint32_t(data.size()));
  const size_t start = out.size();
  out.insert(out.end(), type, type + 4);
  out.insert(out.end(), data.begin(), data.end());
  put32(out, uint32_t(crc32(0L, out.data() + start, uInt(out.size() - start))));
}

}  // namespace

std::vector<uint8_t> encode_png(const Image& img) {
  std::vector<uint8_t> out = {0x89, 'P', 'N', 'G', 0x0D, 0x0A, 0x1A, 0x0A};
  std::vector<uint8_t> ihdr;
  put32(ihdr, uint32_t(img.width));
  put32(ihdr, uint32_t(img.height));
  ihdr.insert(ihdr.end(), {8, 6, 0, 0, 0});   // 8-bit, RGBA, deflate, adaptive filter, no interlace
  chunk(out, "IHDR", ihdr);

  // Filter type 0 (none) per row; zlib does the rest well enough for flat-shaded figures.
  const size_t stride = size_t(img.width) * 4;
  std::vector<uint8_t> raw((stride + 1) * img.height);
  for (int y = 0; y < img.height; ++y) {
    raw[y * (stride + 1)] = 0;
    std::memcpy(&raw[y * (stride + 1) + 1], &img.rgba[y * stride], stride);
  }
  uLongf clen = compressBound(uLong(raw.size()));
  std::vector<uint8_t> z(clen);
  if (compress2(z.data(), &clen, raw.data(), uLong(raw.size()), 6) != Z_OK) throw std::runtime_error("PNG compression failed");
  z.resize(clen);
  chunk(out, "IDAT", z);
  chunk(out, "IEND", {});
  return out;
}

void write_png(const Image& img, const std::string& path) {
  const auto bytes = encode_png(img);
  std::ofstream f(path, std::ios::binary);
  if (!f) throw std::runtime_error("cannot write " + path);
  f.write(reinterpret_cast<const char*>(bytes.data()), std::streamsize(bytes.size()));
}

}  // namespace caps
