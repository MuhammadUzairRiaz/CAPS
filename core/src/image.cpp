// PNG encoder (RGBA 8- or 16-bit, zlib via the system library), with dpi, sRGB and text chunks, and a streaming
// animated-PNG writer. Keeps the alpha channel for transparent figures.
#include <cstdio>
#include <zlib.h>

#include <cmath>
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

namespace {

std::vector<uint8_t> deflate_bytes(const std::vector<uint8_t>& raw) {
  uLongf clen = compressBound(uLong(raw.size()));
  std::vector<uint8_t> z(clen);
  if (compress2(z.data(), &clen, raw.data(), uLong(raw.size()), 6) != Z_OK) throw std::runtime_error("PNG compression failed");
  z.resize(clen);
  return z;
}

// Rows with filter 0 (none), 8 or 16 bits per channel (big-endian samples).
std::vector<uint8_t> scanlines(const Image& img, bool sixteen) {
  const size_t px = size_t(img.width) * img.height;
  const size_t stride = size_t(img.width) * 4 * (sixteen ? 2 : 1);
  std::vector<uint8_t> raw((stride + 1) * img.height);
  const bool deep = img.rgba16.size() == px * 4;
  for (int y = 0; y < img.height; ++y) {
    uint8_t* row = &raw[y * (stride + 1)];
    row[0] = 0;
    if (!sixteen) { std::memcpy(row + 1, &img.rgba[size_t(y) * img.width * 4], stride); continue; }
    for (size_t k = 0; k < size_t(img.width) * 4; ++k) {
      const size_t at = size_t(y) * img.width * 4 + k;
      const uint16_t v = deep ? img.rgba16[at] : uint16_t(img.rgba[at] * 257);
      row[1 + 2 * k] = uint8_t(v >> 8);
      row[2 + 2 * k] = uint8_t(v);
    }
  }
  return raw;
}

void header_chunks(std::vector<uint8_t>& out, const Image& img, const PngOptions& o) {
  std::vector<uint8_t> ihdr;
  put32(ihdr, uint32_t(img.width));
  put32(ihdr, uint32_t(img.height));
  ihdr.insert(ihdr.end(), {uint8_t(o.sixteen ? 16 : 8), 6, 0, 0, 0});
  chunk(out, "IHDR", ihdr);
  if (o.srgb) chunk(out, "sRGB", {0});   // perceptual rendering intent
  if (o.dpi > 0) {
    std::vector<uint8_t> phys;
    const uint32_t ppm = uint32_t(std::lround(o.dpi / 0.0254));
    put32(phys, ppm);
    put32(phys, ppm);
    phys.push_back(1);   // metre
    chunk(out, "pHYs", phys);
  }
  for (const auto& [key, text] : o.text) {
    std::vector<uint8_t> t(key.begin(), key.begin() + std::min<size_t>(key.size(), 79));
    t.insert(t.end(), {0, 0, 0, 0, 0});   // keyword NUL, no compression, method 0, empty language NUL, empty translation NUL
    t.insert(t.end(), text.begin(), text.end());
    chunk(out, "iTXt", t);
  }
}

}  // namespace

std::vector<uint8_t> encode_png(const Image& img, const PngOptions& o) {
  std::vector<uint8_t> out = {0x89, 'P', 'N', 'G', 0x0D, 0x0A, 0x1A, 0x0A};
  header_chunks(out, img, o);
  chunk(out, "IDAT", deflate_bytes(scanlines(img, o.sixteen)));
  chunk(out, "IEND", {});
  return out;
}

void write_png(const Image& img, const std::string& path, const PngOptions& o) {
  const auto bytes = encode_png(img, o);
  std::ofstream f(path, std::ios::binary);
  if (!f) throw std::runtime_error("cannot write " + path);
  f.write(reinterpret_cast<const char*>(bytes.data()), std::streamsize(bytes.size()));
}

std::vector<std::pair<std::string, std::string>> read_png_text(const std::string& path) {
  std::ifstream f(path, std::ios::binary);
  if (!f) throw std::runtime_error("cannot read " + path);
  std::string b((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
  std::vector<std::pair<std::string, std::string>> out;
  if (b.size() < 8 || b.compare(1, 3, "PNG") != 0) throw std::runtime_error(path + " is not a PNG file");
  size_t at = 8;
  auto be = [&](size_t k) { return uint32_t(uint8_t(b[k])) << 24 | uint32_t(uint8_t(b[k + 1])) << 16 | uint32_t(uint8_t(b[k + 2])) << 8 | uint8_t(b[k + 3]); };
  while (at + 12 <= b.size()) {
    const uint32_t n = be(at);
    const std::string type = b.substr(at + 4, 4);
    if (at + 12 + n > b.size()) break;
    const std::string data = b.substr(at + 8, n);
    if (type == "tEXt") {
      const size_t z = data.find('\0');
      if (z != std::string::npos) out.push_back({data.substr(0, z), data.substr(z + 1)});
    } else if (type == "iTXt") {
      const size_t z = data.find('\0');
      if (z != std::string::npos && z + 2 < data.size() && data[z + 1] == 0) {
        const size_t lang = data.find('\0', z + 3);
        const size_t trans = lang == std::string::npos ? lang : data.find('\0', lang + 1);
        if (trans != std::string::npos) out.push_back({data.substr(0, z), data.substr(trans + 1)});
      }
    }
    if (type == "IEND") break;
    at += 12 + n;
  }
  return out;
}

// ---------------------------------------------------------------- APNG
struct ApngWriter::Impl {
  std::ofstream f;
  std::string path;
  int fps, loops, w = 0, h = 0;
  uint32_t seq = 0;
  std::streampos actl_at = 0;
  PngOptions o;
  bool closed = false;
};

ApngWriter::ApngWriter(const std::string& path, int fps, int loops, const PngOptions& o) : impl_(new Impl) {
  impl_->path = path;
  impl_->fps = std::max(1, fps);
  impl_->loops = std::max(0, loops);
  impl_->o = o;
  impl_->o.sixteen = false;
  impl_->f.open(path, std::ios::binary);
  if (!impl_->f) { delete impl_; impl_ = nullptr; throw std::runtime_error("cannot write " + path); }
}

ApngWriter::~ApngWriter() {
  try { close(); } catch (...) {}
  delete impl_;
}

void ApngWriter::add(const Image& img) {
  Impl& I = *impl_;
  std::vector<uint8_t> out;
  if (frames_ == 0) {
    I.w = img.width;
    I.h = img.height;
    out = {0x89, 'P', 'N', 'G', 0x0D, 0x0A, 0x1A, 0x0A};
    std::vector<uint8_t> head;
    header_chunks(head, img, I.o);
    // IHDR first, then acTL (frame count patched on close), then the rest of the header chunks
    const size_t ihdr_len = 8 + 13 + 4;
    out.insert(out.end(), head.begin(), head.begin() + ihdr_len);
    I.f.write(reinterpret_cast<const char*>(out.data()), std::streamsize(out.size()));
    out.clear();
    I.actl_at = I.f.tellp();
    std::vector<uint8_t> actl;
    put32(actl, 0);
    put32(actl, uint32_t(I.loops));
    chunk(out, "acTL", actl);
    out.insert(out.end(), head.begin() + ihdr_len, head.end());
  } else if (img.width != I.w || img.height != I.h) {
    throw std::runtime_error("every frame of an animated PNG must be the size of the first");
  }
  std::vector<uint8_t> fctl;
  put32(fctl, I.seq++);
  put32(fctl, uint32_t(I.w));
  put32(fctl, uint32_t(I.h));
  put32(fctl, 0);
  put32(fctl, 0);
  fctl.push_back(0); fctl.push_back(1);                       // delay 1 / fps s
  fctl.push_back(uint8_t(I.fps >> 8)); fctl.push_back(uint8_t(I.fps));
  fctl.push_back(0);                                          // dispose none
  fctl.push_back(0);                                          // blend source
  chunk(out, "fcTL", fctl);
  const auto z = deflate_bytes(scanlines(img, false));
  if (frames_ == 0) chunk(out, "IDAT", z);
  else {
    std::vector<uint8_t> d;
    put32(d, I.seq++);
    d.insert(d.end(), z.begin(), z.end());
    chunk(out, "fdAT", d);
  }
  I.f.write(reinterpret_cast<const char*>(out.data()), std::streamsize(out.size()));
  ++frames_;
}

void ApngWriter::close() {
  if (!impl_ || impl_->closed) return;
  Impl& I = *impl_;
  I.closed = true;
  if (frames_ == 0) { I.f.close(); std::remove(I.path.c_str()); return; }
  std::vector<uint8_t> end;
  chunk(end, "IEND", {});
  I.f.write(reinterpret_cast<const char*>(end.data()), std::streamsize(end.size()));
  std::vector<uint8_t> actl;
  put32(actl, uint32_t(frames_));
  put32(actl, uint32_t(I.loops));
  std::vector<uint8_t> c;
  chunk(c, "acTL", actl);
  I.f.seekp(I.actl_at);
  I.f.write(reinterpret_cast<const char*>(c.data()), std::streamsize(c.size()));
  I.f.close();
}

}  // namespace caps
