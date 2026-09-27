// Binary trajectories: GROMACS .xtc (compressed, the xdr3dfcoord scheme of GROMACS's libxdrf) and .trr (XDR, single or
// double precision), and CHARMM/NAMD/LAMMPS .dcd (Fortran records, either byte order). Coordinates only: the atoms,
// bonds and types come from the structure or topology given with them.
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>

#include "caps/io.hpp"

namespace caps {

namespace {

// ---------------------------------------------------------------------------------------------------------------- XDR

struct Xdr {
  std::ifstream in;
  size_t size = 0;
  explicit Xdr(const std::string& path) : in(path, std::ios::binary) {
    if (!in) throw ReadError("cannot open " + path);
    in.seekg(0, std::ios::end);
    size = size_t(in.tellg());
    in.seekg(0);
  }
  bool eof() { return in.peek() == std::char_traits<char>::eof(); }
  double where() { return size ? double(in.tellg()) / double(size) : 1.0; }
  uint32_t u32() {
    unsigned char b[4];
    if (!in.read(reinterpret_cast<char*>(b), 4)) throw ReadError("unexpected end of file");
    return (uint32_t(b[0]) << 24) | (uint32_t(b[1]) << 16) | (uint32_t(b[2]) << 8) | uint32_t(b[3]);
  }
  int32_t i32() { return int32_t(u32()); }
  float f32() {
    const uint32_t u = u32();
    float f;
    std::memcpy(&f, &u, 4);
    return f;
  }
  double f64() {
    const uint64_t hi = u32(), lo = u32();
    const uint64_t u = (hi << 32) | lo;
    double d;
    std::memcpy(&d, &u, 8);
    return d;
  }
  double real(bool dbl) { return dbl ? f64() : double(f32()); }
  void bytes(std::vector<unsigned char>& out, size_t n) {
    out.resize(n);
    if (n && !in.read(reinterpret_cast<char*>(out.data()), std::streamsize(n))) throw ReadError("unexpected end of file");
    const size_t pad = (4 - n % 4) % 4;
    if (pad) in.ignore(std::streamsize(pad));
  }
};

// GROMACS box rows (nm) as a cell (Å)
Cell gmx_cell(const double box[9]) {
  Cell c;
  c.a = {box[0] * 10, box[1] * 10, box[2] * 10};
  c.b = {box[3] * 10, box[4] * 10, box[5] * 10};
  c.c = {box[6] * 10, box[7] * 10, box[8] * 10};
  return c;
}

// ---------------------------------------------------------------------------------------------------------------- XTC

const int kMagicInts[] = {0,       0,       0,       0,       0,       0,        0,        0,        0,       8,       10,      12,      16,
                          20,      25,      32,      40,      50,      64,       80,       101,      128,     161,     203,     256,     322,
                          406,     512,     645,     812,     1024,    1290,     1625,     2048,     2580,    3250,    4096,    5060,    6501,
                          8192,    10321,   13003,   16384,   20642,   26007,    32768,    41285,    52015,   65536,   82570,   104031,  131072,
                          165140,  208063,  262144,  330280,  416127,  524287,   660561,   832255,   1048576, 1321122, 1664510, 2097152, 2642245,
                          3329021, 4194304, 5284491, 6658042, 8388607, 10568983, 13316085, 16777216};
constexpr int kFirstIdx = 9;
constexpr int kLastIdx = int(sizeof(kMagicInts) / sizeof(kMagicInts[0]));

struct BitReader {
  const unsigned char* p;
  size_t n, cnt = 0;
  unsigned int lastbits = 0, lastbyte = 0;
  unsigned char next() {
    if (cnt >= n) throw ReadError("xtc: compressed coordinates end early");
    return p[cnt++];
  }
  int bits(int num_of_bits) {
    const int mask = num_of_bits >= 32 ? -1 : (1 << num_of_bits) - 1;
    int num = 0;
    while (num_of_bits >= 8) {
      lastbyte = (lastbyte << 8) | next();
      num |= int((lastbyte >> lastbits) << (num_of_bits - 8));
      num_of_bits -= 8;
    }
    if (num_of_bits > 0) {
      if (int(lastbits) < num_of_bits) {
        lastbits += 8;
        lastbyte = (lastbyte << 8) | next();
      }
      lastbits -= unsigned(num_of_bits);
      num |= int((lastbyte >> lastbits) & ((1u << num_of_bits) - 1));
    }
    return num & mask;
  }
  void ints(int num_of_ints, int num_of_bits, const unsigned int sizes[], int nums[]) {
    int bytes[32];
    bytes[1] = bytes[2] = bytes[3] = 0;
    int nb = 0;
    while (num_of_bits > 8) {
      bytes[nb++] = bits(8);
      num_of_bits -= 8;
    }
    if (num_of_bits > 0) bytes[nb++] = bits(num_of_bits);
    for (int i = num_of_ints - 1; i > 0; --i) {
      unsigned int num = 0;
      for (int j = nb - 1; j >= 0; --j) {
        num = (num << 8) | unsigned(bytes[j]);
        const unsigned int q = num / sizes[i];
        bytes[j] = int(q);
        num = num - q * sizes[i];
      }
      nums[i] = int(num);
    }
    nums[0] = bytes[0] | (bytes[1] << 8) | (bytes[2] << 16) | (bytes[3] << 24);
  }
};

int size_of_int(unsigned int size) {
  unsigned int num = 1;
  int bits = 0;
  while (size >= num && bits < 32) {
    ++bits;
    num <<= 1;
  }
  return bits;
}

int size_of_ints(int num_of_ints, const unsigned int sizes[]) {
  unsigned int bytes[32];
  int num_of_bytes = 1, bits = 0;
  bytes[0] = 1;
  for (int i = 0; i < num_of_ints; ++i) {
    unsigned int tmp = 0;
    int bytecnt;
    for (bytecnt = 0; bytecnt < num_of_bytes; ++bytecnt) {
      tmp = bytes[bytecnt] * sizes[i] + tmp;
      bytes[bytecnt] = tmp & 0xff;
      tmp >>= 8;
    }
    while (tmp != 0) {
      bytes[bytecnt++] = tmp & 0xff;
      tmp >>= 8;
    }
    num_of_bytes = bytecnt;
  }
  unsigned int num = 1;
  --num_of_bytes;
  while (bytes[num_of_bytes] >= num) {
    ++bits;
    num *= 2;
  }
  return bits + num_of_bytes * 8;
}

// One frame's coordinates (nm), as GROMACS's xdr3dfcoord decodes them.
void xtc_coords(Xdr& x, size_t natoms, std::vector<float>& out) {
  const int lsize = x.i32();
  if (size_t(lsize) != natoms) throw ReadError("xtc: a frame holds " + std::to_string(lsize) + " atoms, the header " + std::to_string(natoms));
  out.assign(3 * natoms, 0.0f);
  if (natoms <= 9) {   // small systems are stored uncompressed
    for (auto& v : out) v = x.f32();
    return;
  }
  const float precision = x.f32();
  int minint[3], maxint[3];
  for (int& v : minint) v = x.i32();
  for (int& v : maxint) v = x.i32();
  unsigned int sizeint[3], bitsizeint[3] = {0, 0, 0};
  int bitsize = 0;
  bool large = false;
  for (int k = 0; k < 3; ++k) {
    sizeint[k] = unsigned(maxint[k] - minint[k] + 1);
    if (sizeint[k] > 0xffffff) large = true;
  }
  if (large) {
    for (int k = 0; k < 3; ++k) bitsizeint[k] = unsigned(size_of_int(sizeint[k]));
    bitsize = 0;
  } else {
    bitsize = size_of_ints(3, sizeint);
  }
  int smallidx = x.i32();
  if (smallidx < kFirstIdx || smallidx >= kLastIdx) throw ReadError("xtc: bad compression index");
  int smaller = kMagicInts[std::max(kFirstIdx, smallidx - 1)] / 2;
  int smallnum = kMagicInts[smallidx] / 2;
  unsigned int sizesmall[3];
  sizesmall[0] = sizesmall[1] = sizesmall[2] = unsigned(kMagicInts[smallidx]);
  const int byte_cnt = x.i32();
  if (byte_cnt < 0) throw ReadError("xtc: bad byte count");
  std::vector<unsigned char> buf;
  x.bytes(buf, size_t(byte_cnt));
  BitReader br{buf.data(), buf.size()};
  const float inv = 1.0f / precision;
  size_t i = 0, o = 0;
  int run = 0;
  int prev[3], cur[3];
  auto put = [&](const int* c) {
    if (o + 3 > out.size()) throw ReadError("xtc: more coordinates than atoms");
    out[o++] = float(c[0]) * inv;
    out[o++] = float(c[1]) * inv;
    out[o++] = float(c[2]) * inv;
  };
  while (i < natoms) {
    if (bitsize == 0) {
      cur[0] = br.bits(int(bitsizeint[0]));
      cur[1] = br.bits(int(bitsizeint[1]));
      cur[2] = br.bits(int(bitsizeint[2]));
    } else {
      br.ints(3, bitsize, sizeint, cur);
    }
    ++i;
    for (int k = 0; k < 3; ++k) cur[k] += minint[k], prev[k] = cur[k];
    const int flag = br.bits(1);
    int is_smaller = 0;
    if (flag == 1) {
      run = br.bits(5);
      is_smaller = run % 3;
      run -= is_smaller;
      --is_smaller;
    }
    if (run > 0) {
      for (int k = 0; k < run; k += 3) {
        br.ints(3, smallidx, sizesmall, cur);
        ++i;
        for (int j = 0; j < 3; ++j) cur[j] += prev[j] - smallnum;
        if (k == 0) {   // the first two atoms of a run are swapped (water compresses better that way)
          for (int j = 0; j < 3; ++j) std::swap(cur[j], prev[j]);
          put(prev);
        } else {
          for (int j = 0; j < 3; ++j) prev[j] = cur[j];
        }
        put(cur);
      }
    } else {
      put(cur);
    }
    smallidx += is_smaller;
    if (smallidx < kFirstIdx || smallidx >= kLastIdx) throw ReadError("xtc: bad compression index");
    if (is_smaller < 0) {
      smallnum = smaller;
      smaller = smallidx > kFirstIdx ? kMagicInts[smallidx - 1] / 2 : 0;
    } else if (is_smaller > 0) {
      smaller = smallnum;
      smallnum = kMagicInts[smallidx] / 2;
    }
    sizesmall[0] = sizesmall[1] = sizesmall[2] = unsigned(kMagicInts[smallidx]);
  }
  if (o != out.size()) throw ReadError("xtc: fewer coordinates than atoms");
}

void check_atoms(const System& top, size_t n, const std::string& path) {
  if (top.atoms.size() != n)
    throw ReadError(path + " has " + std::to_string(n) + " atoms, the structure given with it " + std::to_string(top.atoms.size()));
}

}  // namespace

Trajectory read_xtc(const std::string& path, const System& topology, size_t max_frames, const std::function<bool(double, const Trajectory&)>& progress) {
  Xdr x(path);
  Trajectory t;
  t.topology = topology;
  std::vector<float> c;
  while (!x.eof()) {
    if (x.i32() != 1995) throw ReadError(path + ": not an .xtc file (magic number)");
    const size_t n = size_t(x.i32());
    if (t.frames() == 0) check_atoms(topology, n, path);
    const int step = x.i32();
    (void)x.f32();   // time, ps
    double box[9];
    for (double& b : box) b = x.f32();
    xtc_coords(x, n, c);
    std::vector<Vec3> p(n);
    for (size_t i = 0; i < n; ++i) p[i] = {c[3 * i] * 10.0, c[3 * i + 1] * 10.0, c[3 * i + 2] * 10.0};
    t.positions.push_back(std::move(p));
    t.cells.push_back(gmx_cell(box));
    t.timesteps.push_back(step);
    if (progress && !progress(x.where(), t)) {
      t.topology.notes.push_back("reading stopped after " + std::to_string(t.frames()) + " frames");
      break;
    }
    if (max_frames && t.frames() >= max_frames) break;
  }
  if (t.frames() == 0) throw ReadError(path + ": no frames");
  t.topology.cell = t.cells.front();
  t.topology.unwrapped = false;
  t.topology.source_format = "xtc";
  t.topology.notes.push_back(std::to_string(t.frames()) + " frames from " + std::filesystem::path(path).filename().string() + " (GROMACS xtc, 0.001 nm precision by default)");
  return t;
}

Trajectory read_trr(const std::string& path, const System& topology, size_t max_frames, const std::function<bool(double, const Trajectory&)>& progress) {
  Xdr x(path);
  Trajectory t;
  t.topology = topology;
  bool velocities = false;
  while (!x.eof()) {
    if (x.i32() != 1993) throw ReadError(path + ": not a .trr file (magic number)");
    const int slen = x.i32();
    std::vector<unsigned char> version;
    const int vlen = x.i32();
    if (slen <= 0 || vlen < 0 || vlen > 256) throw ReadError(path + ": bad .trr header");
    x.bytes(version, size_t(vlen));
    int sz[13];
    for (int& v : sz) v = x.i32();
    const int box_size = sz[2], vir_size = sz[3], pres_size = sz[4], x_size = sz[7], v_size = sz[8], f_size = sz[9], natoms = sz[10], step = sz[11];
    if (natoms <= 0) throw ReadError(path + ": no atoms in a .trr frame");
    if (t.frames() == 0) check_atoms(topology, size_t(natoms), path);
    int real = 4;
    if (box_size) real = box_size / 9;
    else if (x_size) real = x_size / (3 * natoms);
    else if (v_size) real = v_size / (3 * natoms);
    else if (f_size) real = f_size / (3 * natoms);
    const bool dbl = real == 8;
    (void)x.real(dbl);   // time
    (void)x.real(dbl);   // lambda
    double box[9] = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    if (box_size)
      for (double& b : box) b = x.real(dbl);
    for (int k = 0; k < (vir_size ? 9 : 0) + (pres_size ? 9 : 0); ++k) (void)x.real(dbl);
    std::vector<Vec3> p;
    if (x_size) {
      p.resize(size_t(natoms));
      for (auto& q : p) q = {x.real(dbl) * 10, x.real(dbl) * 10, x.real(dbl) * 10};
    }
    if (v_size) {
      std::vector<Vec3> v(static_cast<size_t>(natoms));
      for (auto& q : v) q = {x.real(dbl) * 1e-2, x.real(dbl) * 1e-2, x.real(dbl) * 1e-2};   // nm/ps → Å/fs
      velocities = true;
      if (t.frames() == 0 || p.empty()) t.topology.velocities = v;
    }
    if (f_size)
      for (int k = 0; k < 3 * natoms; ++k) (void)x.real(dbl);
    if (p.empty()) continue;   // a frame of velocities or forces only
    t.positions.push_back(std::move(p));
    t.cells.push_back(box_size ? gmx_cell(box) : topology.cell);
    t.timesteps.push_back(step);
    if (progress && !progress(x.where(), t)) {
      t.topology.notes.push_back("reading stopped after " + std::to_string(t.frames()) + " frames");
      break;
    }
    if (max_frames && t.frames() >= max_frames) break;
  }
  if (t.frames() == 0) throw ReadError(path + ": no frames with coordinates");
  t.topology.cell = t.cells.front();
  t.topology.unwrapped = false;
  t.topology.source_format = "trr";
  t.topology.notes.push_back(std::to_string(t.frames()) + " frames from " + std::filesystem::path(path).filename().string() + " (GROMACS trr" +
                             (velocities ? ", velocities of the first frame kept" : "") + ")");
  return t;
}

namespace {

// Fortran unformatted records, either byte order (the first record is 84 bytes long)
struct Fortran {
  std::ifstream in;
  bool swap = false;
  size_t size = 0;
  explicit Fortran(const std::string& path) : in(path, std::ios::binary) {
    if (!in) throw ReadError("cannot open " + path);
    in.seekg(0, std::ios::end);
    size = size_t(in.tellg());
    in.seekg(0);
  }
  uint32_t raw() {
    unsigned char b[4];
    if (!in.read(reinterpret_cast<char*>(b), 4)) throw ReadError("unexpected end of file");
    return swap ? (uint32_t(b[0]) << 24) | (uint32_t(b[1]) << 16) | (uint32_t(b[2]) << 8) | b[3]
                : (uint32_t(b[3]) << 24) | (uint32_t(b[2]) << 16) | (uint32_t(b[1]) << 8) | b[0];
  }
  std::vector<unsigned char> record() {
    const uint32_t n = raw();
    std::vector<unsigned char> r(n);
    if (n && !in.read(reinterpret_cast<char*>(r.data()), n)) throw ReadError("unexpected end of file");
    if (raw() != n) throw ReadError("a Fortran record's end marker does not match its start");
    return r;
  }
  int32_t i32(const unsigned char* p) const {
    uint32_t u = swap ? (uint32_t(p[0]) << 24) | (uint32_t(p[1]) << 16) | (uint32_t(p[2]) << 8) | p[3]
                      : (uint32_t(p[3]) << 24) | (uint32_t(p[2]) << 16) | (uint32_t(p[1]) << 8) | p[0];
    return int32_t(u);
  }
  float f32(const unsigned char* p) const {
    const uint32_t u = uint32_t(i32(p));
    float f;
    std::memcpy(&f, &u, 4);
    return f;
  }
  double f64(const unsigned char* p) const {
    unsigned char b[8];
    for (int k = 0; k < 8; ++k) b[k] = swap ? p[7 - k] : p[k];
    double d;
    std::memcpy(&d, b, 8);
    return d;
  }
};

}  // namespace

Trajectory read_dcd(const std::string& path, const System& topology, size_t max_frames, const std::function<bool(double, const Trajectory&)>& progress) {
  Fortran f(path);
  {
    unsigned char b[4];
    if (!f.in.read(reinterpret_cast<char*>(b), 4)) throw ReadError(path + ": empty");
    const uint32_t le = (uint32_t(b[3]) << 24) | (uint32_t(b[2]) << 16) | (uint32_t(b[1]) << 8) | b[0];
    const uint32_t be = (uint32_t(b[0]) << 24) | (uint32_t(b[1]) << 16) | (uint32_t(b[2]) << 8) | b[3];
    if (le == 84) f.swap = false;
    else if (be == 84) f.swap = true;
    else throw ReadError(path + ": not a .dcd file (first record)");
    f.in.seekg(0);
  }
  const auto h = f.record();
  if (h.size() != 84 || std::memcmp(h.data(), "CORD", 4) != 0) throw ReadError(path + ": not a coordinate .dcd file");
  int icntrl[20];
  for (int k = 0; k < 20; ++k) icntrl[k] = f.i32(h.data() + 4 + 4 * k);
  const bool charmm = icntrl[19] != 0;
  const bool has_cell = charmm && icntrl[10] != 0;
  const bool four_d = charmm && icntrl[11] != 0;
  const int nfixed = icntrl[8];
  const int nsavc = std::max(1, icntrl[2]);
  const int istart = icntrl[1];
  if (nfixed != 0) throw ReadError(path + ": fixed atoms in a .dcd are not read");
  (void)f.record();   // titles
  const auto na = f.record();
  if (na.size() != 4) throw ReadError(path + ": bad atom-count record");
  const size_t n = size_t(f.i32(na.data()));
  check_atoms(topology, n, path);
  Trajectory t;
  t.topology = topology;
  int64_t frame = 0;
  while (f.in.peek() != std::char_traits<char>::eof()) {
    Cell cell = topology.cell;
    if (has_cell) {
      const auto u = f.record();
      if (u.size() != 48) throw ReadError(path + ": bad unit-cell record");
      double d[6];
      for (int k = 0; k < 6; ++k) d[k] = f.f64(u.data() + 8 * k);
      // A, γ, B, β, α, C: angles as cosines (LAMMPS, newer CHARMM/NAMD) or in degrees
      const double la = d[0], lb = d[2], lc = d[5];
      auto ang = [](double v) { return std::abs(v) <= 1.0 ? std::acos(v) : v * M_PI / 180.0; };
      const double gam = ang(d[1]), bet = ang(d[3]), alp = ang(d[4]);
      Cell c;
      c.a = {la, 0, 0};
      c.b = {lb * std::cos(gam), lb * std::sin(gam), 0};
      const double cx = lc * std::cos(bet), cy = lc * (std::cos(alp) - std::cos(bet) * std::cos(gam)) / std::sin(gam);
      c.c = {cx, cy, std::sqrt(std::max(0.0, lc * lc - cx * cx - cy * cy))};
      if (la > 0 && lb > 0 && lc > 0) cell = c;
    }
    std::vector<Vec3> p(n);
    for (int k = 0; k < 3; ++k) {
      const auto r = f.record();
      if (r.size() != 4 * n) throw ReadError(path + ": a coordinate record does not match the atom count");
      for (size_t i = 0; i < n; ++i) p[i][k] = f.f32(r.data() + 4 * i);
    }
    if (four_d) (void)f.record();
    t.positions.push_back(std::move(p));
    t.cells.push_back(cell);
    t.timesteps.push_back(int64_t(istart) + frame * nsavc);
    ++frame;
    if (progress && !progress(f.size ? double(f.in.tellg()) / double(f.size) : 1.0, t)) {
      t.topology.notes.push_back("reading stopped after " + std::to_string(t.frames()) + " frames");
      break;
    }
    if (max_frames && t.frames() >= max_frames) break;
  }
  if (t.frames() == 0) throw ReadError(path + ": no frames");
  t.topology.cell = t.cells.front();
  t.topology.unwrapped = false;
  t.topology.source_format = "dcd";
  t.topology.notes.push_back(std::to_string(t.frames()) + " frames from " + std::filesystem::path(path).filename().string() + " (DCD" +
                             (has_cell ? ", with the cell of each frame" : ", no cell in the file: the structure's") + ")");
  return t;
}

}  // namespace caps
