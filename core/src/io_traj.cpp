// Trajectories in the common formats: see caps/io.hpp (write_trajectory).
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <sstream>

#include "caps/io.hpp"

namespace caps {

namespace {

std::string lower_ext(const std::string& path) {
  const auto dot = path.find_last_of('.');
  std::string e = dot == std::string::npos ? "" : path.substr(dot + 1);
  for (auto& c : e) c = char(std::tolower(static_cast<unsigned char>(c)));
  return e;
}

// XDR: big-endian 32-bit words
struct Xdr {
  std::ofstream& out;
  void i32(int32_t v) {
    const uint32_t u = uint32_t(v);
    const unsigned char b[4] = {static_cast<unsigned char>(u >> 24), static_cast<unsigned char>(u >> 16), static_cast<unsigned char>(u >> 8),
                                static_cast<unsigned char>(u)};
    out.write(reinterpret_cast<const char*>(b), 4);
  }
  void f32(float f) {
    int32_t v;
    std::memcpy(&v, &f, 4);
    i32(v);
  }
  void str(const char* s) {   // xdr string: length, then the bytes padded to a multiple of four
    const int32_t n = int32_t(std::strlen(s));
    i32(n);
    out.write(s, n);
    for (int32_t k = n; k % 4; ++k) out.put('\0');
  }
};

}  // namespace

void write_trajectory(const Trajectory& t, const std::string& path, double dt_fs) {
  const std::string e = lower_ext(path);
  if (e == "lammpstrj" || e == "dump") { write_lammps_dump(t, path); return; }
  if (e == "dcd") { write_dcd(t, path, dt_fs); return; }
  const size_t nf = t.frames();
  if (nf == 0) throw ReadError("no frames to write");
  auto time_ps = [&](size_t k) { return k < t.timesteps.size() ? double(t.timesteps[k]) * dt_fs / 1000.0 : double(k); };
  char b[128];
  if (e == "xyz") {
    std::ofstream out(path);
    if (!out) throw ReadError("cannot write " + path);
    for (size_t k = 0; k < nf; ++k) {
      std::snprintf(b, sizeof b, "Time=%.6g", time_ps(k));   // ps
      write_xyz_frame(out, t.frame(k), b);
    }
    return;
  }
  if (e == "gro") {
    std::ofstream out(path);
    if (!out) throw ReadError("cannot write " + path);
    for (size_t k = 0; k < nf; ++k) {
      std::snprintf(b, sizeof b, "CAPS trajectory t= %.5f step= %lld", time_ps(k), static_cast<long long>(k < t.timesteps.size() ? t.timesteps[k] : int64_t(k)));
      write_gro_frame(out, t.frame(k), b);
    }
    return;
  }
  if (e == "pdb") {
    std::FILE* f = std::fopen(path.c_str(), "w");
    if (!f) throw ReadError("cannot write " + path);
    for (size_t k = 0; k < nf; ++k) write_pdb_frame(f, t.frame(k), int(k + 1), k + 1 == nf);
    std::fprintf(f, "END\n");
    std::fclose(f);
    return;
  }
  if (e == "trr") {
    // GROMACS .trr (single precision): per frame the header — magic 1993, the version string, the block sizes, natoms,
    // step, nre, time, lambda — then the box (3 × 3, nm) and the positions (nm)
    std::ofstream out(path, std::ios::binary);
    if (!out) throw ReadError("cannot write " + path);
    Xdr x{out};
    const int32_t n = int32_t(t.topology.atoms.size());
    for (size_t k = 0; k < nf; ++k) {
      const System f = t.frame(k);
      const bool box = f.cell.valid();
      x.i32(1993);
      x.i32(13);   // the string's length plus its terminator, as GROMACS writes it
      x.str("GMX_trn_file");
      const int32_t sizes[] = {0 /*ir*/, 0 /*e*/, box ? 9 * 4 : 0 /*box*/, 0 /*vir*/, 0 /*pres*/, 0 /*top*/, 0 /*sym*/, n * 3 * 4 /*x*/, 0 /*v*/, 0 /*f*/};
      for (int32_t s : sizes) x.i32(s);
      x.i32(n);
      x.i32(int32_t(k < t.timesteps.size() ? t.timesteps[k] : int64_t(k)));
      x.i32(0);   // nre
      x.f32(float(time_ps(k)));
      x.f32(0.0f);   // lambda
      if (box)
        for (const Vec3* v : {&f.cell.a, &f.cell.b, &f.cell.c})
          for (int d = 0; d < 3; ++d) x.f32(float((*v)[d] / 10.0));
      for (const auto& a : f.atoms)
        for (int d = 0; d < 3; ++d) x.f32(float(a.pos[d] / 10.0));
    }
    return;
  }
  throw ReadError("trajectory format '" + e + "': write .lammpstrj, .dcd, .xyz, .pdb, .gro or .trr");
}

}  // namespace caps
