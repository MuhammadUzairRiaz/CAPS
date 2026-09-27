// CAPS figure bundles (see caps/bundle.hpp).
#include "caps/bundle.hpp"
#include "caps/provenance.hpp"

#include <zlib.h>

#include <algorithm>
#include <array>
#include <cctype>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <functional>
#include <sstream>
#include <stdexcept>

#include "caps/analysis.hpp"
#include "caps/io.hpp"
#include "caps/json.hpp"

#ifndef CAPS_VERSION_STRING
#define CAPS_VERSION_STRING "0.1.0"
#endif

namespace caps {

// ---------------------------------------------------------------- SHA-256 (FIPS 180-4)

std::string sha256_hex(const std::string& bytes) {
  static const uint32_t K[64] = {
      0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5, 0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74,
      0x80deb1fe, 0x9bdc06a7, 0xc19bf174, 0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da, 0x983e5152, 0xa831c66d,
      0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967, 0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e,
      0x92722c85, 0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070, 0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5,
      0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3, 0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2};
  uint32_t h[8] = {0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19};
  std::string m = bytes;
  const uint64_t bits = uint64_t(bytes.size()) * 8;
  m.push_back(char(0x80));
  while (m.size() % 64 != 56) m.push_back(0);
  for (int k = 7; k >= 0; --k) m.push_back(char((bits >> (8 * k)) & 0xFF));
  auto rotr = [](uint32_t x, int n) { return (x >> n) | (x << (32 - n)); };
  for (size_t off = 0; off < m.size(); off += 64) {
    uint32_t w[64];
    for (int t = 0; t < 16; ++t)
      w[t] = (uint32_t(uint8_t(m[off + 4 * t])) << 24) | (uint32_t(uint8_t(m[off + 4 * t + 1])) << 16) | (uint32_t(uint8_t(m[off + 4 * t + 2])) << 8) |
             uint32_t(uint8_t(m[off + 4 * t + 3]));
    for (int t = 16; t < 64; ++t) {
      const uint32_t s0 = rotr(w[t - 15], 7) ^ rotr(w[t - 15], 18) ^ (w[t - 15] >> 3);
      const uint32_t s1 = rotr(w[t - 2], 17) ^ rotr(w[t - 2], 19) ^ (w[t - 2] >> 10);
      w[t] = w[t - 16] + s0 + w[t - 7] + s1;
    }
    uint32_t a = h[0], b = h[1], c = h[2], d = h[3], e = h[4], f = h[5], g = h[6], hh = h[7];
    for (int t = 0; t < 64; ++t) {
      const uint32_t S1 = rotr(e, 6) ^ rotr(e, 11) ^ rotr(e, 25), ch = (e & f) ^ (~e & g), t1 = hh + S1 + ch + K[t] + w[t];
      const uint32_t S0 = rotr(a, 2) ^ rotr(a, 13) ^ rotr(a, 22), mj = (a & b) ^ (a & c) ^ (b & c), t2 = S0 + mj;
      hh = g; g = f; f = e; e = d + t1; d = c; c = b; b = a; a = t1 + t2;
    }
    h[0] += a; h[1] += b; h[2] += c; h[3] += d; h[4] += e; h[5] += f; h[6] += g; h[7] += hh;
  }
  char out[65];
  for (int k = 0; k < 8; ++k) std::snprintf(out + 8 * k, 9, "%08x", h[k]);
  return std::string(out, 64);
}

// ---------------------------------------------------------------- zip

namespace {
void put16(std::string& s, uint16_t v) { s.push_back(char(v & 255)); s.push_back(char(v >> 8)); }
void put32(std::string& s, uint32_t v) { for (int k = 0; k < 4; ++k) s.push_back(char((v >> (8 * k)) & 255)); }
uint16_t get16(const std::string& s, size_t p) { return uint16_t(uint8_t(s[p]) | (uint8_t(s[p + 1]) << 8)); }
uint32_t get32(const std::string& s, size_t p) { return uint32_t(get16(s, p)) | (uint32_t(get16(s, p + 2)) << 16); }

std::string deflate_raw(const std::string& in) {
  z_stream z{};
  if (deflateInit2(&z, 6, Z_DEFLATED, -15, 8, Z_DEFAULT_STRATEGY) != Z_OK) throw std::runtime_error("deflate");
  std::string out(deflateBound(&z, uLong(in.size())) + 16, '\0');
  z.next_in = reinterpret_cast<Bytef*>(const_cast<char*>(in.data()));
  z.avail_in = uInt(in.size());
  z.next_out = reinterpret_cast<Bytef*>(out.data());
  z.avail_out = uInt(out.size());
  deflate(&z, Z_FINISH);
  out.resize(z.total_out);
  deflateEnd(&z);
  return out;
}

std::string inflate_raw(const std::string& in, size_t size) {
  z_stream z{};
  if (inflateInit2(&z, -15) != Z_OK) throw std::runtime_error("inflate");
  std::string out(size, '\0');
  z.next_in = reinterpret_cast<Bytef*>(const_cast<char*>(in.data()));
  z.avail_in = uInt(in.size());
  z.next_out = reinterpret_cast<Bytef*>(out.data());
  z.avail_out = uInt(out.size());
  const int r = inflate(&z, Z_FINISH);
  inflateEnd(&z);
  if (r != Z_STREAM_END) throw std::runtime_error("corrupt zip entry");
  return out;
}
}  // namespace

void write_zip(const std::string& path, const std::vector<std::pair<std::string, std::string>>& files) {
  std::string out, central;
  for (const auto& [name, data] : files) {
    const uint32_t crc = uint32_t(crc32(0, reinterpret_cast<const Bytef*>(data.data()), uInt(data.size())));
    const std::string comp = deflate_raw(data);
    const uint32_t offset = uint32_t(out.size());
    put32(out, 0x04034b50); put16(out, 20); put16(out, 0x0800); put16(out, 8); put16(out, 0); put16(out, 0x21);   // UTF-8 names, deflate, 1980-01-01
    put32(out, crc); put32(out, uint32_t(comp.size())); put32(out, uint32_t(data.size())); put16(out, uint16_t(name.size())); put16(out, 0);
    out += name;
    out += comp;
    put32(central, 0x02014b50); put16(central, 20); put16(central, 20); put16(central, 0x0800); put16(central, 8); put16(central, 0); put16(central, 0x21);
    put32(central, crc); put32(central, uint32_t(comp.size())); put32(central, uint32_t(data.size())); put16(central, uint16_t(name.size()));
    put16(central, 0); put16(central, 0); put16(central, 0); put16(central, 0); put32(central, 0); put32(central, offset);
    central += name;
  }
  const uint32_t cd = uint32_t(out.size());
  out += central;
  put32(out, 0x06054b50); put16(out, 0); put16(out, 0); put16(out, uint16_t(files.size())); put16(out, uint16_t(files.size()));
  put32(out, uint32_t(central.size())); put32(out, cd); put16(out, 0);
  std::ofstream f(path, std::ios::binary);
  if (!f) throw std::runtime_error("cannot write " + path);
  f.write(out.data(), std::streamsize(out.size()));
}

std::map<std::string, std::string> read_zip(const std::string& path) {
  std::ifstream f(path, std::ios::binary);
  if (!f) throw std::runtime_error("cannot open " + path);
  const std::string s((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
  size_t eocd = std::string::npos;
  for (size_t p = s.size() >= 22 ? s.size() - 22 : 0; p != std::string::npos && s.size() - p < 70000; --p) {
    if (get32(s, p) == 0x06054b50) { eocd = p; break; }
    if (p == 0) break;
  }
  if (eocd == std::string::npos) throw std::runtime_error(path + " is not a zip file");
  const uint16_t count = get16(s, eocd + 10);
  size_t p = get32(s, eocd + 16);
  std::map<std::string, std::string> files;
  for (int k = 0; k < count; ++k) {
    if (p + 46 > s.size() || get32(s, p) != 0x02014b50) throw std::runtime_error("corrupt zip directory");
    const uint16_t method = get16(s, p + 10);
    const uint32_t csize = get32(s, p + 20), usize = get32(s, p + 24);
    const uint16_t nlen = get16(s, p + 28), xlen = get16(s, p + 30), clen = get16(s, p + 32);
    const uint32_t local = get32(s, p + 42);
    const std::string name = s.substr(p + 46, nlen);
    const size_t data = local + 30 + get16(s, local + 26) + get16(s, local + 28);
    const std::string comp = s.substr(data, csize);
    files[name] = method == 0 ? comp : inflate_raw(comp, usize);
    p += 46 + nlen + xlen + clen;
  }
  return files;
}

// ---------------------------------------------------------------- bundle

std::string table_csv(const DataTable& t) {
  std::string out;
  for (size_t c = 0; c < t.columns.size(); ++c) {
    std::string col = t.columns[c];
    if (col.find(',') != std::string::npos) col = "\"" + col + "\"";
    out += (c ? "," : "") + col;
  }
  out += "\n";
  char b[40];
  for (const auto& r : t.rows) {
    for (size_t c = 0; c < r.size(); ++c) {
      std::snprintf(b, sizeof b, "%.10g", r[c]);
      out += (c ? "," : "") + std::string(b);
    }
    out += "\n";
  }
  return out;
}

namespace {
std::string read_file(const std::string& p) {
  std::ifstream f(p, std::ios::binary);
  if (!f) throw std::runtime_error("cannot read " + p);
  return std::string((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
}

// The pipeline's result on the recorded frame, as the bundle's data files see it.
PipelineState run_for_bundle(const Trajectory& traj, const Pipeline& pipeline, int frame) {
  frame = std::clamp(frame, 0, int(traj.frames()) - 1);
  System s = traj.frame(size_t(frame));
  if (!s.unwrapped && s.cell.valid()) make_molecules_whole(s);
  return run_pipeline(s, pipeline, frame, frame < int(traj.timesteps.size()) ? traj.timesteps[size_t(frame)] : 0, &traj);
}

std::string safe(std::string n) {
  for (char& c : n) if (!(std::isalnum(static_cast<unsigned char>(c)) || c == '_' || c == '-' || c == '.')) c = '_';
  return n;
}
}  // namespace

std::vector<BundleFile> bundle_files(const Trajectory& traj, const Pipeline& pipeline, const BundleOptions& opt, const Camera& cam,
                                     const RenderOptions& render) {
  std::vector<BundleFile> files;
  const auto st = run_for_bundle(traj, pipeline, opt.frame);
  if (opt.include_figures) {
    RenderOptions r = render;
    r.width = opt.width;
    r.height = opt.height;
    r.colours = st.colour;
    for (size_t i = 0; i < r.colours.size(); ++i) if (st.selected[i]) r.colours[i] = 0xE5484D;
    r.segments = st.segments;
    r.highlight.clear();
    r.focus = -1;
    if (r.colour_by == ColourBy::Property) property_values(st, "DistanceToCOM", r.property);
    Renderer renderer;
    const auto png = encode_png(renderer.render(st.system, cam, r));
    files.push_back({"figure.png", std::to_string(opt.width) + " × " + std::to_string(opt.height) + " raster", std::string(png.begin(), png.end())});
    files.push_back({"figure.svg", "vector figure", render_svg(st.system, cam, r)});
  }
  if (opt.include_data)
    for (const auto& t : st.tables)
      files.push_back({"data/" + safe(t.name) + ".csv", std::to_string(t.rows.size()) + " rows · " + t.title, table_csv(t)});
  const std::string pj = pipeline_to_json(pipeline).dump(2);
  if (opt.include_pipeline) files.push_back({"pipeline.json", "the pipeline as run (" + std::to_string(pipeline.steps.size()) + " steps)", pj});
  std::string input_sha;
  if (!opt.input.empty() && std::filesystem::exists(opt.input)) {
    const std::string bytes = read_file(opt.input);
    input_sha = sha256_hex(bytes);
    if (opt.include_input) {
      files.push_back({"input/" + std::filesystem::path(opt.input).filename().string(), "the input structure", bytes});
      if (!opt.topology.empty() && std::filesystem::exists(opt.topology))
        files.push_back({"input/" + std::filesystem::path(opt.topology).filename().string(), "its topology", read_file(opt.topology)});
    }
  }
  // provenance: what went in, what came out, every file's hash
  Json prov = Json::object();
  prov["caps"] = CAPS_VERSION_STRING;
  {
    char b[32];
    const std::time_t now = std::time(nullptr);
    std::strftime(b, sizeof b, "%Y-%m-%dT%H:%M:%SZ", std::gmtime(&now));
    prov["created"] = std::string(b);
  }
  Json in = Json::object();
  in["file"] = opt.input.empty() ? std::string("(built in CAPS)") : std::filesystem::path(opt.input).filename().string();
  if (!opt.topology.empty()) in["topology"] = std::filesystem::path(opt.topology).filename().string();
  in["sha256"] = input_sha;
  in["atoms"] = double(traj.topology.atoms.size());
  in["frames"] = double(traj.frames());
  in["frame"] = double(opt.frame);
  in["included"] = opt.include_input;
  prov["input"] = std::move(in);
  const Cell& c = traj.cells.empty() ? traj.topology.cell : traj.cells[size_t(std::clamp(opt.frame, 0, int(traj.cells.size()) - 1))];
  Json box = Json::array();
  for (const Vec3& e : {c.a, c.b, c.c}) box.push_back(norm(e));
  prov["box_A"] = std::move(box);
  Json pp = Json::object();
  pp["steps"] = double(pipeline.steps.size());
  pp["sha256"] = sha256_hex(pj);
  // the methods the steps use: their references beside the data
  const auto refs = pipeline_citations(pipeline);
  if (!refs.empty()) {
    Json ra = Json::array();
    for (const auto& k : refs) ra.push_back(k);
    pp["references"] = std::move(ra);
    files.push_back({"references.bib", "the references of the methods the pipeline uses", bibtex(refs)});
  }
  prov["pipeline"] = std::move(pp);
  Json fig = Json::array();
  fig.push_back(double(opt.width));
  fig.push_back(double(opt.height));
  prov["figure_px"] = std::move(fig);
  Json hashes = Json::object();
  for (const auto& f : files) hashes[f.name] = sha256_hex(f.bytes);
  prov["files"] = std::move(hashes);
  files.push_back({"provenance.json", "inputs, pipeline and every file's sha256", prov.dump(2)});
  if (opt.include_readme) {
    std::string r = "CAPS figure bundle · " + opt.name + "\n\n";
    r += "figure.png, figure.svg   the view of frame " + std::to_string(opt.frame) + " with the pipeline applied\n";
    r += "data/*.csv               the data behind each plot, one file per table\n";
    r += "pipeline.json            the steps as run (caps pipeline FILE --steps pipeline.json)\n";
    r += "input/                   the structure, when included\n";
    r += "provenance.json          CAPS version, input hash, pipeline hash and the sha256 of every file\n";
    if (!refs.empty()) r += "references.bib           the references of the methods the pipeline uses\n";
    r += "\n";
    r += "Reproduce: caps reproduce " + opt.name + ".caps-bundle.zip\n";
    r += "  rebuilds the data files from input/ and pipeline.json and compares their sha256 with provenance.json.\n";
    files.push_back({"README.txt", "what is inside, how to reproduce it", r});
  }
  return files;
}

void write_bundle(const std::string& path, const std::vector<BundleFile>& files) {
  std::vector<std::pair<std::string, std::string>> z;
  for (const auto& f : files) z.push_back({f.name, f.bytes});
  write_zip(path, z);
}

bool reproduce_bundle(const std::string& path, std::vector<std::string>& report) {
  const auto files = read_zip(path);
  auto it = files.find("provenance.json");
  if (it == files.end()) throw std::runtime_error("no provenance.json in " + path);
  const Json prov = Json::parse(it->second);
  const Json& in = prov["input"];
  if (!in.has("included") || !in["included"].boolean()) throw std::runtime_error("the bundle does not include its input structure");
  auto pit = files.find("pipeline.json");
  if (pit == files.end()) throw std::runtime_error("the bundle does not include its pipeline");
  // the input back on disk, beside its topology, so the readers see the same files
  const auto dir = std::filesystem::temp_directory_path() / ("caps_reproduce_" + std::to_string(std::hash<std::string>{}(path)));
  std::filesystem::create_directories(dir);
  auto extract = [&](const std::string& name) {
    const auto f = files.find("input/" + name);
    if (f == files.end()) throw std::runtime_error("input/" + name + " missing from the bundle");
    const auto p = dir / name;
    std::ofstream o(p, std::ios::binary);
    o.write(f->second.data(), std::streamsize(f->second.size()));
    return p.string();
  };
  const std::string input = extract(in["file"].str());
  if (sha256_hex(files.at("input/" + in["file"].str())) != in["sha256"].str()) report.push_back("input differs from its recorded sha256");
  const std::string topo = in.has("topology") ? extract(in["topology"].str()) : "";
  const Trajectory traj = open_file(input, topo);
  const auto st = run_for_bundle(traj, pipeline_from_json(Json::parse(pit->second)), int(in.num("frame", 0)));
  bool ok = true;
  size_t checked = 0;
  for (const auto& [name, hash] : prov["files"].members()) {
    if (name.rfind("data/", 0) != 0) continue;
    std::string now;
    for (const auto& t : st.tables) if ("data/" + safe(t.name) + ".csv" == name) now = table_csv(t);
    const bool same = !now.empty() && sha256_hex(now) == hash.str();
    report.push_back(name + (same ? "  match" : "  differs"));
    ok &= same;
    ++checked;
  }
  std::filesystem::remove_all(dir);
  if (!checked) report.push_back("no data files to compare");
  return ok && checked > 0;
}

}  // namespace caps
