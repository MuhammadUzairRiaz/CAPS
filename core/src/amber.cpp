// AMBER prmtop and inpcrd / rst7 readers (see amber.hpp). The prmtop is a list of %FLAG sections, each with a Fortran
// %FORMAT (20a4, 10I8, 5E16.8 …) read by its fixed field widths. Units as the file: charges × 18.2223 (e·√(kcal Å/mol)),
// bond k in kcal/(mol Å²) for k (r − r0)², angle k in kcal/(mol rad²) for k (θ − θ0)², torsions V [1 + cos(nφ − δ)],
// Lennard-Jones A/r¹² − B/r⁶ per type pair.
#include "caps/amber.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <map>
#include <set>
#include <sstream>

#include "caps/elements.hpp"
#include "caps/io.hpp"

namespace caps {

namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr double kChargeUnit = 18.2223;          // AMBER's internal charge unit per e (√332.0522173)
constexpr double kVelUnit = 20.455;              // AMBER velocity unit: Å per (1/20.455 ps)

std::string lower(std::string s) {
  for (auto& c : s) c = char(std::tolower(static_cast<unsigned char>(c)));
  return s;
}

std::string trim(const std::string& s) {
  const auto a = s.find_first_not_of(" \t\r\n");
  if (a == std::string::npos) return "";
  return s.substr(a, s.find_last_not_of(" \t\r\n") - a + 1);
}

struct Section {
  char kind = 'a';   // a text, I integer, E / F real
  int width = 80;
  std::vector<std::string> fields;
};

// %FORMAT(10I8) → I, 8; (5E16.8) → E, 16; (20a4) → a, 4; (1a80) → a, 80
void parse_format(const std::string& f, Section& s, const std::string& where) {
  const auto a = f.find('('), b = f.find(')');
  if (a == std::string::npos || b == std::string::npos) throw ReadError(where + ": bad %FORMAT " + f);
  const std::string in = f.substr(a + 1, b - a - 1);
  size_t k = 0;
  while (k < in.size() && std::isdigit(static_cast<unsigned char>(in[k]))) ++k;
  if (k >= in.size()) throw ReadError(where + ": bad %FORMAT " + f);
  s.kind = char(std::toupper(static_cast<unsigned char>(in[k])));
  if (s.kind == 'A') s.kind = 'a';
  size_t e = k + 1;
  while (e < in.size() && std::isdigit(static_cast<unsigned char>(in[e]))) ++e;
  s.width = std::stoi(in.substr(k + 1, e - k - 1));
  if (s.width <= 0) throw ReadError(where + ": bad %FORMAT " + f);
}

std::map<std::string, Section> read_sections(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw ReadError("cannot open " + path);
  std::map<std::string, Section> out;
  Section* cur = nullptr;
  std::string line, flag;
  size_t ln = 0;
  bool any = false;
  while (std::getline(in, line)) {
    ++ln;
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (line.rfind("%VERSION", 0) == 0 || line.rfind("%COMMENT", 0) == 0) continue;
    if (line.rfind("%FLAG", 0) == 0) {
      flag = trim(line.substr(5));
      cur = &out[flag];
      *cur = Section{};
      any = true;
      continue;
    }
    if (line.rfind("%FORMAT", 0) == 0) {
      if (!cur) throw ReadError(path + ":" + std::to_string(ln) + ": %FORMAT before any %FLAG");
      parse_format(line, *cur, path + ":" + std::to_string(ln));
      continue;
    }
    if (!cur) continue;
    // fixed-width fields; a text field keeps its spaces trimmed, a blank numeric tail ends the line
    for (size_t p = 0; p < line.size(); p += size_t(cur->width)) {
      const std::string f = line.substr(p, size_t(cur->width));
      if (cur->kind == 'a') cur->fields.push_back(trim(f));
      else {
        const std::string t = trim(f);
        if (t.empty()) break;
        cur->fields.push_back(t);
      }
    }
    if (cur->kind == 'a' && cur->width >= 80) {   // a title line: one field, whatever its length
      if (!cur->fields.empty()) cur->fields.back() = trim(line);
    }
  }
  if (!any) throw ReadError(path + ": not an AMBER topology (no %FLAG sections)");
  return out;
}

struct Reader {
  const std::map<std::string, Section>& S;
  std::string path;
  bool has(const std::string& f) const { return S.count(f) > 0; }
  const Section& sec(const std::string& f) const {
    auto it = S.find(f);
    if (it == S.end()) throw ReadError(path + ": no %FLAG " + f);
    return it->second;
  }
  std::vector<long> ints(const std::string& f, size_t want) const {
    const auto& s = sec(f);
    std::vector<long> v;
    v.reserve(s.fields.size());
    for (const auto& x : s.fields) {
      try { v.push_back(std::stol(x)); } catch (...) { throw ReadError(path + ": %FLAG " + f + ": '" + x + "' is not an integer"); }
    }
    if (want != size_t(-1) && v.size() < want)
      throw ReadError(path + ": %FLAG " + f + " has " + std::to_string(v.size()) + " values, " + std::to_string(want) + " expected");
    return v;
  }
  std::vector<double> reals(const std::string& f, size_t want) const {
    const auto& s = sec(f);
    std::vector<double> v;
    v.reserve(s.fields.size());
    for (auto x : s.fields) {
      std::replace(x.begin(), x.end(), 'D', 'E');   // Fortran double exponents
      try { v.push_back(std::stod(x)); } catch (...) { throw ReadError(path + ": %FLAG " + f + ": '" + x + "' is not a number"); }
    }
    if (want != size_t(-1) && v.size() < want)
      throw ReadError(path + ": %FLAG " + f + " has " + std::to_string(v.size()) + " values, " + std::to_string(want) + " expected");
    return v;
  }
  std::vector<std::string> texts(const std::string& f, size_t want) const {
    const auto& s = sec(f);
    if (want != size_t(-1) && s.fields.size() < want)
      throw ReadError(path + ": %FLAG " + f + " has " + std::to_string(s.fields.size()) + " entries, " + std::to_string(want) + " expected");
    return s.fields;
  }
};

Cell box_cell(double a, double b, double c, double al, double be, double ga) {
  const double ca = std::cos(al * kPi / 180), cb = std::cos(be * kPi / 180), cg = std::cos(ga * kPi / 180), sg = std::sin(ga * kPi / 180);
  Cell cell;
  cell.a = {a, 0, 0};
  cell.b = {b * cg, b * sg, 0};
  const double cx = c * cb, cy = c * (ca - cb * cg) / sg;
  cell.c = {cx, cy, std::sqrt(std::max(0.0, c * c - cx * cx - cy * cy))};
  return cell;
}

// ---- NetCDF-3 (classic "CDF\1" and 64-bit offset "CDF\2"), big-endian: the header, then variables by name
// (Unidata, "The NetCDF Classic Format Specification")
struct NcAtt {
  int type = 0;
  std::string text;
  std::vector<double> values;
};
struct NcVar {
  std::string name;
  std::vector<int> dims;
  std::map<std::string, NcAtt> atts;
  int type = 0;
  uint64_t begin = 0, vsize = 0;
  bool record = false;
};
struct NcFile {
  std::string path;
  int version = 0;
  uint64_t numrecs = 0, recsize = 0;
  std::vector<std::pair<std::string, uint64_t>> dims;
  int recdim = -1;
  std::map<std::string, NcAtt> atts;
  std::vector<NcVar> vars;
  std::ifstream in;

  const NcVar* var(const std::string& n) const {
    for (const auto& v : vars)
      if (v.name == n) return &v;
    return nullptr;
  }
  uint64_t dim(const std::string& n) const {
    for (const auto& [k, v] : dims)
      if (k == n) return v;
    return 0;
  }
};

int nc_size(int type) {
  switch (type) {
    case 1: case 2: return 1;
    case 3: return 2;
    case 4: case 5: return 4;
    case 6: return 8;
    default: return 0;
  }
}

class NcHeader {
 public:
  NcHeader(std::ifstream& in, const std::string& path) : in_(in), path_(path) {}
  uint32_t u32() { unsigned char b[4]; get(b, 4); return uint32_t(b[0]) << 24 | uint32_t(b[1]) << 16 | uint32_t(b[2]) << 8 | b[3]; }
  uint64_t u64() { const uint64_t hi = u32(); return hi << 32 | u32(); }
  std::string name() {
    const uint32_t n = u32();
    if (n > 1 << 20) throw ReadError(path_ + ": NetCDF header damaged (a name of " + std::to_string(n) + " bytes)");
    std::string s(n, '\0');
    get(reinterpret_cast<unsigned char*>(s.data()), n);
    skip_pad(n);
    return s;
  }
  NcAtt att_values(int type, uint32_t n) {
    NcAtt a;
    a.type = type;
    const int sz = nc_size(type);
    if (!sz) throw ReadError(path_ + ": NetCDF type " + std::to_string(type) + " is not a NetCDF-3 type");
    std::vector<unsigned char> b(size_t(n) * size_t(sz));
    if (!b.empty()) get(b.data(), b.size());
    skip_pad(uint32_t(b.size()));
    if (type == 2) a.text.assign(b.begin(), b.end());
    else
      for (uint32_t k = 0; k < n; ++k) a.values.push_back(decode(b.data() + size_t(k) * size_t(sz), type));
    while (!a.text.empty() && a.text.back() == '\0') a.text.pop_back();
    return a;
  }
  std::map<std::string, NcAtt> atts() {
    std::map<std::string, NcAtt> out;
    const uint32_t tag = u32(), n = u32();
    if (tag == 0 && n == 0) return out;
    if (tag != 0x0C) throw ReadError(path_ + ": NetCDF header damaged (attribute list)");
    for (uint32_t k = 0; k < n; ++k) {
      const std::string nm = name();
      const int type = int(u32());
      const uint32_t cnt = u32();
      out[nm] = att_values(type, cnt);
    }
    return out;
  }
  static double decode(const unsigned char* b, int type) {
    switch (type) {
      case 1: return double(static_cast<signed char>(b[0]));
      case 2: return double(b[0]);
      case 3: return double(int16_t(uint16_t(b[0]) << 8 | b[1]));
      case 4: return double(int32_t(uint32_t(b[0]) << 24 | uint32_t(b[1]) << 16 | uint32_t(b[2]) << 8 | b[3]));
      case 5: {
        const uint32_t u = uint32_t(b[0]) << 24 | uint32_t(b[1]) << 16 | uint32_t(b[2]) << 8 | b[3];
        float f;
        std::memcpy(&f, &u, 4);
        return double(f);
      }
      case 6: {
        uint64_t u = 0;
        for (int k = 0; k < 8; ++k) u = u << 8 | b[k];
        double d;
        std::memcpy(&d, &u, 8);
        return d;
      }
      default: return 0;
    }
  }

 private:
  void get(unsigned char* b, size_t n) {
    in_.read(reinterpret_cast<char*>(b), std::streamsize(n));
    if (size_t(in_.gcount()) != n) throw ReadError(path_ + ": NetCDF file ends inside its header");
  }
  void skip_pad(uint32_t n) {
    const uint32_t pad = (4 - n % 4) % 4;
    unsigned char b[4];
    if (pad) get(b, pad);
  }
  std::ifstream& in_;
  std::string path_;
};

void open_netcdf(NcFile& f, const std::string& path) {
  f.path = path;
  f.in.open(path, std::ios::binary);
  if (!f.in) throw ReadError("cannot open " + path);
  unsigned char m[4] = {0, 0, 0, 0};
  f.in.read(reinterpret_cast<char*>(m), 4);
  if (m[0] == 0x89 && m[1] == 'H' && m[2] == 'D' && m[3] == 'F')
    throw ReadError(path + ": a NetCDF-4 (HDF5) file — write it as NetCDF-3 (cpptraj: trajout name.nc netcdf; or nccopy -k nc6 in.nc out.nc) to open it");
  if (m[0] != 'C' || m[1] != 'D' || m[2] != 'F') throw ReadError(path + ": not a NetCDF file");
  if (m[3] == 5) throw ReadError(path + ": a CDF-5 (64-bit data) NetCDF file — write it as NetCDF-3 64-bit offset (nccopy -k nc6) to open it");
  if (m[3] != 1 && m[3] != 2) throw ReadError(path + ": unknown NetCDF version " + std::to_string(int(m[3])));
  f.version = m[3];
  NcHeader h(f.in, path);
  f.numrecs = h.u32();
  {   // dimensions
    const uint32_t tag = h.u32(), n = h.u32();
    if (!(tag == 0 && n == 0)) {
      if (tag != 0x0A) throw ReadError(path + ": NetCDF header damaged (dimension list)");
      for (uint32_t k = 0; k < n; ++k) {
        const std::string nm = h.name();
        const uint64_t len = h.u32();
        if (len == 0) f.recdim = int(k);
        f.dims.push_back({nm, len});
      }
    }
  }
  f.atts = h.atts();
  {   // variables
    const uint32_t tag = h.u32(), n = h.u32();
    if (!(tag == 0 && n == 0)) {
      if (tag != 0x0B) throw ReadError(path + ": NetCDF header damaged (variable list)");
      for (uint32_t k = 0; k < n; ++k) {
        NcVar v;
        v.name = h.name();
        const uint32_t nd = h.u32();
        for (uint32_t d = 0; d < nd; ++d) {
          const uint32_t id = h.u32();
          if (id >= f.dims.size()) throw ReadError(path + ": NetCDF variable " + v.name + " refers to dimension " + std::to_string(id));
          v.dims.push_back(int(id));
        }
        v.atts = h.atts();
        v.type = int(h.u32());
        v.vsize = h.u32();
        v.begin = f.version == 2 ? h.u64() : h.u32();
        v.record = !v.dims.empty() && v.dims.front() == f.recdim;
        // the size per record (or of the whole variable) from the dimensions: the header's vsize saturates past 4 GiB
        uint64_t cnt = 1;
        for (size_t d = v.record ? 1 : 0; d < v.dims.size(); ++d) cnt *= f.dims[size_t(v.dims[d])].second;
        const uint64_t bytes = cnt * uint64_t(nc_size(v.type));
        v.vsize = (bytes + 3) / 4 * 4;
        f.vars.push_back(std::move(v));
      }
    }
  }
  int nrec = 0;
  for (const auto& v : f.vars)
    if (v.record) f.recsize += v.vsize, ++nrec;
  if (nrec == 1)   // one record variable: its records are not padded
    for (const auto& v : f.vars)
      if (v.record) {
        uint64_t cnt = 1;
        for (size_t d = 1; d < v.dims.size(); ++d) cnt *= f.dims[size_t(v.dims[d])].second;
        f.recsize = cnt * uint64_t(nc_size(v.type));
      }
  if (f.numrecs == 0xFFFFFFFFu && f.recsize > 0) {   // streaming: the records the file holds
    uint64_t first = UINT64_MAX;
    for (const auto& v : f.vars)
      if (v.record) first = std::min(first, v.begin);
    const uint64_t size = uint64_t(std::filesystem::file_size(path));
    f.numrecs = size > first ? (size - first) / f.recsize : 0;
  }
}

// one record (or the whole of a non-record variable) as doubles
std::vector<double> nc_read(NcFile& f, const NcVar& v, uint64_t rec) {
  uint64_t cnt = 1;
  for (size_t d = v.record ? 1 : 0; d < v.dims.size(); ++d) cnt *= f.dims[size_t(v.dims[d])].second;
  const int sz = nc_size(v.type);
  if (!sz) throw ReadError(f.path + ": variable " + v.name + " has NetCDF type " + std::to_string(v.type));
  const uint64_t off = v.begin + (v.record ? rec * f.recsize : 0);
  std::vector<unsigned char> b(size_t(cnt) * size_t(sz));
  f.in.clear();
  f.in.seekg(std::streamoff(off));
  f.in.read(reinterpret_cast<char*>(b.data()), std::streamsize(b.size()));
  if (size_t(f.in.gcount()) != b.size()) throw ReadError(f.path + ": the file ends inside " + v.name + (v.record ? " of frame " + std::to_string(rec + 1) : ""));
  std::vector<double> out(static_cast<size_t>(cnt));
  for (size_t k = 0; k < size_t(cnt); ++k) out[k] = NcHeader::decode(b.data() + k * size_t(sz), v.type);
  auto it = v.atts.find("scale_factor");
  if (it != v.atts.end() && !it->second.values.empty() && it->second.values[0] != 1)
    for (auto& x : out) x *= it->second.values[0];
  return out;
}

std::string nc_convention(const NcFile& f) {
  auto it = f.atts.find("Conventions");
  return it == f.atts.end() ? std::string() : it->second.text;
}

}  // namespace

bool is_amber_topology_path(const std::string& path) {
  std::string p = lower(path);
  if (p.size() > 3 && p.substr(p.size() - 3) == ".gz") p = p.substr(0, p.size() - 3);
  const std::string e = std::filesystem::path(p).extension().string();
  return e == ".prmtop" || e == ".parm7" || e == ".prmtp";
}

bool is_amber_coordinates_path(const std::string& path) {
  const std::string e = lower(std::filesystem::path(path).extension().string());
  return e == ".inpcrd" || e == ".rst7" || e == ".restrt" || e == ".crd" || e == ".rst" || e == ".ncrst" || e == ".mdcrd_rst";
}

AmberTopology read_amber_prmtop(const std::string& path) {
  const auto S = read_sections(path);
  const Reader R{S, path};
  const std::string name = std::filesystem::path(path).filename().string();
  auto refuse = [&](const std::string& what) { throw ReadError(path + ": " + what); };
  // what the force field needs and CAPS cannot hold: the structure still opens, without the file's force field
  std::string ff_off;
  auto ff_stop = [&](const std::string& why) { if (ff_off.empty()) ff_off = why; };
  if (R.has("CTITLE") || R.has("FORCE_FIELD_TYPE") || R.has("CHARMM_UREY_BRADLEY"))
    ff_stop("a CHARMM topology converted by chamber (Urey–Bradley, harmonic impropers, CHARMM 1-4 parameters) — open the CHARMM files for its force field");
  for (const char* f : {"CHARMM_CMAP_COUNT", "CMAP_COUNT"})
    if (R.has(f)) ff_stop("CMAP correction maps (ff19SB, CHARMM): CAPS has no CMAP term, and leaving them out would change the energy");
  if (R.has("POLARIZABILITY") || R.has("DIPOLE_DAMP_FACTOR")) ff_stop("a polarisable topology (induced dipoles)");
  if (R.has("LENNARD_JONES_CCOEF")) {   // 12-6-4 ions (Li & Merz): C/r⁴ on the type pairs listed
    for (double c : R.reals("LENNARD_JONES_CCOEF", 0))
      if (c != 0) { ff_stop("12-6-4 ion terms (LENNARD_JONES_CCOEF, C/r⁴): CAPS has no r⁻⁴ pair term, and leaving them out would change the energy"); break; }
  }

  const auto P = R.ints("POINTERS", 30);
  const size_t natom = size_t(P[0]), ntypes = size_t(P[1]), nbonh = size_t(P[2]), mbona = size_t(P[3]), ntheth = size_t(P[4]),
               mtheta = size_t(P[5]), nphih = size_t(P[6]), mphia = size_t(P[7]), nnb = size_t(P[10]), nres = size_t(P[11]),
               numbnd = size_t(P[15]), numang = size_t(P[16]), nptra = size_t(P[17]), nphb = size_t(P[19]);
  const long ifpert = P[20], ifbox = P[27], ifcap = P.size() > 29 ? P[29] : 0, numextra = P.size() > 30 ? P[30] : 0;
  if (ifpert) refuse("a perturbation topology (IFPERT) is not read");
  if (ifcap) refuse("a solvent-cap topology (IFCAP) is not read");
  if (numextra > 0) ff_stop(std::to_string(numextra) + " extra points (TIP4P / TIP5P virtual sites)");
  if (natom == 0) refuse("no atoms");

  const auto names = R.texts("ATOM_NAME", natom);
  auto charge = R.reals("CHARGE", natom);
  const auto mass = R.reals("MASS", natom);
  const auto tindex = R.ints("ATOM_TYPE_INDEX", natom);
  const auto nexcl = R.ints("NUMBER_EXCLUDED_ATOMS", natom);
  const auto nbindex = R.ints("NONBONDED_PARM_INDEX", ntypes * ntypes);
  const auto reslab = R.texts("RESIDUE_LABEL", nres);
  const auto resptr = R.ints("RESIDUE_POINTER", nres);
  const auto bk = R.reals("BOND_FORCE_CONSTANT", numbnd), br = R.reals("BOND_EQUIL_VALUE", numbnd);
  const auto ak = R.reals("ANGLE_FORCE_CONSTANT", numang), at = R.reals("ANGLE_EQUIL_VALUE", numang);
  const auto dk = R.reals("DIHEDRAL_FORCE_CONSTANT", nptra), dn = R.reals("DIHEDRAL_PERIODICITY", nptra), dp = R.reals("DIHEDRAL_PHASE", nptra);
  std::vector<double> scee(nptra, 1.2), scnb(nptra, 2.0);
  const bool has_scale = R.has("SCEE_SCALE_FACTOR") && R.has("SCNB_SCALE_FACTOR");
  if (has_scale) scee = R.reals("SCEE_SCALE_FACTOR", nptra), scnb = R.reals("SCNB_SCALE_FACTOR", nptra);
  const size_t npair = ntypes * (ntypes + 1) / 2;
  const auto lja = R.reals("LENNARD_JONES_ACOEF", npair), ljb = R.reals("LENNARD_JONES_BCOEF", npair);
  if (nphb > 0 && R.has("HBOND_ACOEF")) {
    const auto ha = R.reals("HBOND_ACOEF", 0), hb = R.reals("HBOND_BCOEF", 0);
    for (size_t k = 0; k < std::min(ha.size(), hb.size()); ++k)
      if (ha[k] != 0 || hb[k] != 0) { ff_stop("10-12 hydrogen-bond terms with non-zero coefficients"); break; }
  }
  const auto amber_type = R.has("AMBER_ATOM_TYPE") ? R.texts("AMBER_ATOM_TYPE", natom) : names;
  const auto excl = R.ints("EXCLUDED_ATOMS_LIST", nnb);
  std::vector<long> atomic_number;
  if (R.has("ATOMIC_NUMBER")) atomic_number = R.ints("ATOMIC_NUMBER", natom);

  AmberTopology out;
  auto& notes = out.notes;
  System& s = out.system;
  s.source_format = "amber-prmtop";
  if (R.has("TITLE") && !R.sec("TITLE").fields.empty()) s.title = R.sec("TITLE").fields.front();
  else if (R.has("CTITLE") && !R.sec("CTITLE").fields.empty()) s.title = R.sec("CTITLE").fields.front();

  auto ff = std::make_shared<ForceField>();
  ForceField& F = *ff;
  F.name = "AMBER topology " + name;
  F.pair_form = "lj12-6";
  F.mixing = "arithmetic";
  F.native_pair = "lj/cut/coul/long";
  F.native_dihedral = "fourier";
  F.native_improper = "cvff";
  F.native_special = "amber";
  F.native_cutoff = 10;

  // types: an AMBER atom type with its Lennard-Jones index and mass (hydrogen mass repartitioning gives one type two masses)
  std::map<std::tuple<std::string, long, double>, int> tkey;
  std::map<std::string, int> name_uses;
  std::vector<long> type_lj;   // per CAPS type: the file's LJ type (0-based)
  F.atom_type.resize(natom);
  F.type_index.resize(natom);
  F.charge.resize(natom);
  F.mass.resize(natom);
  F.why.assign(natom, "type from the AMBER topology");
  for (size_t i = 0; i < natom; ++i) {
    const long lt = tindex[i] - 1;
    if (lt < 0 || size_t(lt) >= ntypes) refuse("atom " + std::to_string(i + 1) + " has Lennard-Jones type " + std::to_string(tindex[i]) + " of " + std::to_string(ntypes));
    const auto key = std::make_tuple(amber_type[i], lt, mass[i]);
    auto it = tkey.find(key);
    if (it == tkey.end()) {
      const int uses = ++name_uses[amber_type[i]];
      const std::string tn = uses == 1 ? amber_type[i] : amber_type[i] + "_" + std::to_string(uses);
      it = tkey.emplace(key, int(F.type_names.size())).first;
      F.type_names.push_back(tn);
      type_lj.push_back(lt);
    }
    F.type_index[i] = it->second;
    F.atom_type[i] = F.type_names[size_t(it->second)];
    F.charge[i] = charge[i] / kChargeUnit;
    F.mass[i] = mass[i];
  }
  for (const auto& [n, u] : name_uses)
    if (u > 1) notes.push_back("AMBER type " + n + " comes with " + std::to_string(u) + " masses or Lennard-Jones types: kept apart as " + n + ", " + n + "_2 …");

  // Lennard-Jones: A = 4εσ¹², B = 4εσ⁶ for each type with itself; a pair the arithmetic (Lorentz–Berthelot) rule does
  // not reproduce (NBFIX, off-diagonal edits) keeps its own coefficients
  auto ab = [&](long a, long b, double& A, double& B) {
    const long ix = nbindex[size_t(a) * ntypes + size_t(b)];
    if (ix < 0) { A = B = 0; return false; }   // a 10-12 pair (checked zero above)
    if (ix == 0 || size_t(ix) > npair) refuse("NONBONDED_PARM_INDEX " + std::to_string(ix) + " out of range");
    A = lja[size_t(ix - 1)], B = ljb[size_t(ix - 1)];
    return true;
  };
  auto eps_sigma = [](double A, double B) -> PairType {
    if (A <= 0 || B <= 0) return {0, 0};
    return {B * B / (4 * A), std::pow(A / B, 1.0 / 6)};
  };
  const size_t nt = F.type_names.size();
  F.lj.resize(nt);
  for (size_t t = 0; t < nt; ++t) {
    double A, B;
    ab(type_lj[t], type_lj[t], A, B);
    if ((A > 0) != (B > 0)) ff_stop("type " + F.type_names[t] + ": a purely repulsive or purely attractive Lennard-Jones pair (A " + std::to_string(A) + ", B " + std::to_string(B) + ") has no ε, σ");
    F.lj[t] = eps_sigma(A, B);
  }
  int overrides = 0;
  for (size_t a = 0; a < nt; ++a)
    for (size_t b = a + 1; b < nt; ++b) {
      double A, B;
      ab(type_lj[a], type_lj[b], A, B);
      const PairType &p = F.lj[a], &q = F.lj[b];
      const double e = std::sqrt(p.eps * q.eps), sg = 0.5 * (p.sigma + q.sigma);
      const double Am = 4 * e * std::pow(sg, 12), Bm = 4 * e * std::pow(sg, 6);
      auto close = [](double x, double y) { return std::fabs(x - y) <= 2e-6 * std::max({std::fabs(x), std::fabs(y), 1e-30}); };
      if (close(A, Am) && close(B, Bm)) continue;
      if ((A > 0) != (B > 0)) ff_stop("types " + F.type_names[a] + " and " + F.type_names[b] + ": a Lennard-Jones pair with only A or only B has no ε, σ");
      F.pair_override[{int(a), int(b)}] = A > 0 ? eps_sigma(A, B) : PairType{0, 1};
      ++overrides;
    }
  if (overrides) notes.push_back(std::to_string(overrides) + " type pairs with their own Lennard-Jones coefficients (not the arithmetic rule)");

  // atoms: residues, molecules
  s.atoms.resize(natom);
  std::vector<long> resstart(resptr.begin(), resptr.begin() + long(nres));
  for (size_t r = 0; r < nres; ++r) {
    const size_t a0 = size_t(resstart[r] - 1), a1 = r + 1 < nres ? size_t(resstart[r + 1] - 1) : natom;
    if (a0 > a1 || a1 > natom) refuse("RESIDUE_POINTER out of order");
    for (size_t i = a0; i < a1; ++i) s.atoms[i].resname = reslab[r], s.atoms[i].resid = long(r + 1);
  }
  std::vector<long> mol(natom, 0);
  if (ifbox > 0 && R.has("ATOMS_PER_MOLECULE")) {
    const auto sp = R.ints("SOLVENT_POINTERS", 3);
    const auto apm = R.ints("ATOMS_PER_MOLECULE", size_t(sp[1]));
    size_t i = 0;
    for (size_t m = 0; m < size_t(sp[1]); ++m)
      for (long k = 0; k < apm[m] && i < natom; ++k) mol[i++] = long(m + 1);
    if (i != natom) refuse("ATOMS_PER_MOLECULE covers " + std::to_string(i) + " of " + std::to_string(natom) + " atoms");
  }

  // bonded terms (coordinate indices: 3 × atom)
  const size_t bad = size_t(-1);
  auto atom = [&](long c) -> uint32_t {
    const long a = std::labs(c);
    if (a % 3 != 0 || size_t(a / 3) >= natom) return uint32_t(bad);
    return uint32_t(a / 3);
  };
  auto check = [&](uint32_t a, const char* what) {
    if (a == uint32_t(bad)) refuse(std::string("a ") + what + " refers to an atom out of range");
    return a;
  };
  std::set<std::pair<uint32_t, uint32_t>> bonded;
  auto bonds = [&](const char* flag, size_t n) {
    const auto v = R.ints(flag, 3 * n);
    for (size_t k = 0; k < n; ++k) {
      const uint32_t i = check(atom(v[3 * k]), "bond"), j = check(atom(v[3 * k + 1]), "bond");
      const long t = v[3 * k + 2] - 1;
      if (t < 0 || size_t(t) >= numbnd) refuse("bond type out of range");
      F.bonds.push_back({i, j, bk[size_t(t)], br[size_t(t)]});
      s.bonds.push_back({i, j, 1});
      bonded.insert({std::min(i, j), std::max(i, j)});
    }
  };
  bonds("BONDS_INC_HYDROGEN", nbonh);
  bonds("BONDS_WITHOUT_HYDROGEN", mbona);
  auto angles = [&](const char* flag, size_t n) {
    const auto v = R.ints(flag, 4 * n);
    for (size_t k = 0; k < n; ++k) {
      const long t = v[4 * k + 3] - 1;
      if (t < 0 || size_t(t) >= numang) refuse("angle type out of range");
      F.angles.push_back({check(atom(v[4 * k]), "angle"), check(atom(v[4 * k + 1]), "angle"), check(atom(v[4 * k + 2]), "angle"), ak[size_t(t)], at[size_t(t)]});
    }
  };
  angles("ANGLES_INC_HYDROGEN", ntheth);
  angles("ANGLES_WITHOUT_HYDROGEN", mtheta);
  // torsions: a negative third index — the 1-4 pair is not computed here (another term or a ring already counts it); a
  // negative fourth — an improper. AMBER impropers are periodic torsions too; one whose phase is not 0 or 180° is kept
  // with the proper torsions (the same energy, and an exact LAMMPS form)
  std::set<std::pair<uint32_t, uint32_t>> p14;
  std::set<std::pair<double, double>> scales;   // (SCEE, SCNB) of the torsions that carry a 1-4 pair
  int improper_as_proper = 0;
  auto torsions = [&](const char* flag, size_t n) {
    const auto v = R.ints(flag, 5 * n);
    for (size_t k = 0; k < n; ++k) {
      const uint32_t i = check(atom(v[5 * k]), "torsion"), j = check(atom(v[5 * k + 1]), "torsion"), c = check(atom(v[5 * k + 2]), "torsion"),
                     l = check(atom(v[5 * k + 3]), "torsion");
      const long t = v[5 * k + 4] - 1;
      if (t < 0 || size_t(t) >= nptra) refuse("torsion type out of range");
      const bool no14 = v[5 * k + 2] < 0, improper = v[5 * k + 3] < 0;
      const double per = dn[size_t(t)];
      if (std::fabs(per - std::round(per)) > 1e-6) ff_stop("a torsion with periodicity " + std::to_string(per) + " (not whole)");
      const TorsionTerm term{i, j, c, l, dk[size_t(t)], int(std::lround(std::fabs(per))), dp[size_t(t)]};
      if (improper) {
        const double cs = std::cos(term.delta);
        if (std::fabs(std::fabs(cs) - 1) < 1e-9) F.impropers.push_back(term);
        else { F.dihedrals.push_back(term); ++improper_as_proper; }
      } else {
        F.dihedrals.push_back(term);
        if (!no14) {
          p14.insert({std::min(i, l), std::max(i, l)});
          scales.insert({scee[size_t(t)], scnb[size_t(t)]});
        }
      }
    }
  };
  torsions("DIHEDRALS_INC_HYDROGEN", nphih);
  torsions("DIHEDRALS_WITHOUT_HYDROGEN", mphia);
  if (improper_as_proper) notes.push_back(std::to_string(improper_as_proper) + " improper torsions with a phase other than 0 or 180° kept with the proper torsions (the same periodic form)");
  // zero-barrier torsion terms carry only their 1-4 pair: drop the term, keep the pair
  F.dihedrals.erase(std::remove_if(F.dihedrals.begin(), F.dihedrals.end(), [](const TorsionTerm& t) { return t.v == 0; }), F.dihedrals.end());
  if (scales.size() > 1) {
    std::string list;
    for (const auto& [e, v] : scales) list += (list.empty() ? "" : ", ") + std::to_string(e) + "/" + std::to_string(v);
    ff_stop("1-4 scaling differs between torsions (SCEE/SCNB " + list + ", e.g. GLYCAM with a protein force field): CAPS applies one 1-4 scaling to the whole structure");
  }
  const double sc_e = scales.empty() ? 1.2 : scales.begin()->first, sc_n = scales.empty() ? 2.0 : scales.begin()->second;
  if (sc_e <= 0 || sc_n <= 0) ff_stop("a 1-4 scale factor of zero");
  F.coul14 = sc_e > 0 ? 1 / sc_e : 0;
  F.lj14 = sc_n > 0 ? 1 / sc_n : 0;
  if (!has_scale) notes.push_back("no SCEE / SCNB sections (an older topology): the AMBER defaults 1.2 and 2.0");
  for (const auto& [a, b] : p14) F.pairs14.push_back({a, b});

  // exclusions: the file's list (1-2, 1-3 and 1-4 partners; a lone 0 marks an atom with none), made symmetric
  F.excluded.assign(natom, {});
  size_t pos = 0;
  for (size_t i = 0; i < natom; ++i) {
    const long ne = nexcl[i];
    if (ne < 0 || pos + size_t(ne) > excl.size()) refuse("EXCLUDED_ATOMS_LIST is shorter than NUMBER_EXCLUDED_ATOMS says");
    for (long k = 0; k < ne; ++k) {
      const long j = excl[pos + size_t(k)];
      if (j <= 0) continue;
      if (size_t(j) > natom) refuse("an excluded atom out of range");
      F.excluded[i].push_back(uint32_t(j - 1));
      F.excluded[size_t(j - 1)].push_back(uint32_t(i));
    }
    pos += size_t(ne);
  }
  auto add_ex = [&](uint32_t a, uint32_t b) { F.excluded[a].push_back(b), F.excluded[b].push_back(a); };
  int added = 0;
  auto excluded = [&](uint32_t a, uint32_t b) { return std::find(F.excluded[a].begin(), F.excluded[a].end(), b) != F.excluded[a].end(); };
  for (const auto& [a, b] : bonded) if (!excluded(a, b)) add_ex(a, b), ++added;
  for (const auto& pr : F.pairs14) if (!excluded(pr[0], pr[1])) add_ex(pr[0], pr[1]), ++added;
  for (auto& e : F.excluded) {
    std::sort(e.begin(), e.end());
    e.erase(std::unique(e.begin(), e.end()), e.end());
  }
  if (added) notes.push_back(std::to_string(added) + " bonded or 1-4 pairs missing from EXCLUDED_ATOMS_LIST added to the exclusions");

  // atoms of the structure: named by their types (as a GROMACS topology), elements from the file or the mass
  int from_mass = 0;
  for (size_t i = 0; i < natom; ++i) {
    Atom& a = s.atoms[i];
    a.id = int64_t(i + 1);
    a.mol = mol[i];
    a.type = F.type_index[i] + 1;
    a.name = F.atom_type[i];
    a.charge = F.charge[i];
    int z = atomic_number.empty() ? 0 : int(atomic_number[i]);
    if (z <= 0 || z > max_element()) {
      // the atom name's element when its mass fits it (hydrogens repartitioned to ~3 u stay H), else the nearest mass
      const int zn = element_from_name(names[i]);
      if (zn > 0 && (std::fabs(element(zn).mass - mass[i]) < 0.6 || (zn == 1 && mass[i] < 4.1))) z = zn;
      else z = element_from_mass(mass[i], 0.6);
      ++from_mass;
    }
    a.element = z;
  }
  if (from_mass && atomic_number.empty()) notes.push_back("no ATOMIC_NUMBER section: elements from the atom names and masses");
  s.types.clear();
  for (size_t t = 0; t < nt; ++t) {
    TypeInfo ti;
    ti.type = int(t + 1);
    ti.label = F.type_names[t];
    for (size_t i = 0; i < natom; ++i)
      if (size_t(F.type_index[i]) == t) { ti.mass = mass[i]; break; }
    s.types.push_back(ti);
  }
  s.bonds_from_file = true;
  s.has_charges = true;
  if (ifbox > 0) {
    if (R.has("BOX_DIMENSIONS")) {
      const auto b = R.reals("BOX_DIMENSIONS", 4);
      s.cell = box_cell(b[1], b[2], b[3], ifbox == 2 ? 109.4712206 : 90, b[0], ifbox == 2 ? 109.4712206 : 90);
    }
    if (ifbox == 2) notes.push_back("a truncated octahedron (IFBOX 2): the box as its triclinic cell");
  }
  if (mol.front() == 0) {   // no molecule list: bonded fragments
    int count = 0;
    const auto m = s.molecules(&count);
    for (size_t i = 0; i < natom; ++i) s.atoms[i].mol = m[i] + 1;
  }
  s.has_mol = true;

  // net charge as the file holds it (AMBER charges are × 18.2223, so the sum is rarely exactly whole)
  double q = 0;
  for (double x : F.charge) q += x;
  F.notes.push_back(std::string("from ") + name + ": " + std::to_string(natom) + " atoms, " + std::to_string(nt) + " types, " + std::to_string(F.bonds.size()) +
                    " bonds, " + std::to_string(F.angles.size()) + " angles, " + std::to_string(F.dihedrals.size()) + " torsion terms, " +
                    std::to_string(F.impropers.size()) + " impropers, " + std::to_string(F.pairs14.size()) + " 1-4 pairs (Coulomb × 1/" +
                    std::to_string(sc_e).substr(0, 4) + ", LJ × 1/" + std::to_string(sc_n).substr(0, 4) + ")");
  char buf[64];
  std::snprintf(buf, sizeof buf, "%.4f", q);
  F.notes.push_back(std::string("net charge ") + buf + " e");
  for (const auto& n : notes) F.notes.push_back(n);
  if (!ff_off.empty()) {   // the atoms, residues, bonds, charges and masses stand; the force field is not taken
    notes.push_back("the file's force field is not taken — " + ff_off + "; assign one from the library");
    s.notes = notes;
    return out;
  }
  s.notes = notes;
  out.ff = ff;
  return out;
}

AmberCoordinates read_amber_coordinates(const std::string& path) {
  std::ifstream in(path, std::ios::binary);
  if (!in) throw ReadError("cannot open " + path);
  char magic[4] = {0, 0, 0, 0};
  in.read(magic, 4);
  if (std::string(magic, 3) == "CDF" || (magic[0] == '\x89' && magic[1] == 'H')) {   // a NetCDF restart (AMBERRESTART)
    in.close();
    NcFile f;
    open_netcdf(f, path);
    const std::string conv = nc_convention(f);
    if (conv.find("AMBERRESTART") == std::string::npos)
      throw ReadError(path + ": a NetCDF file without the AMBERRESTART convention (Conventions: '" + conv + "') — open a trajectory with its topology instead");
    AmberCoordinates c;
    const NcVar* xv = f.var("coordinates");
    if (!xv) throw ReadError(path + ": no coordinates variable");
    const auto x = nc_read(f, *xv, 0);
    const size_t n = size_t(f.dim("atom"));
    if (x.size() != 3 * n) throw ReadError(path + ": coordinates are not atom × 3");
    for (size_t i = 0; i < n; ++i) c.positions.push_back({x[3 * i], x[3 * i + 1], x[3 * i + 2]});
    if (const NcVar* vv = f.var("velocities")) {   // Å/ps after the file's scale factor (20.455)
      const auto v = nc_read(f, *vv, 0);
      for (size_t i = 0; i < n && v.size() == 3 * n; ++i) c.velocities.push_back(Vec3{v[3 * i], v[3 * i + 1], v[3 * i + 2]} * 1e-3);
    }
    const NcVar *lv = f.var("cell_lengths"), *av = f.var("cell_angles");
    if (lv && av) {
      const auto l = nc_read(f, *lv, 0), a = nc_read(f, *av, 0);
      if (l.size() == 3 && a.size() == 3 && l[0] > 0) c.cell = box_cell(l[0], l[1], l[2], a[0], a[1], a[2]), c.has_box = true;
    }
    if (auto it = f.atts.find("title"); it != f.atts.end()) c.title = it->second.text;
    return c;
  }
  in.clear();
  in.seekg(0);
  AmberCoordinates c;
  std::string line;
  std::getline(in, line);
  if (!line.empty() && line.back() == '\r') line.pop_back();
  c.title = trim(line);
  if (!std::getline(in, line)) throw ReadError(path + ": no atom count line");
  size_t natom = 0;
  {
    std::istringstream is(line);
    long n = 0;
    if (!(is >> n) || n <= 0) throw ReadError(path + ": no atom count on line 2");
    natom = size_t(n);
  }
  // 6F12.7 records; a line may be shorter at the end of a block
  std::vector<double> v;
  v.reserve(natom * 3 + 6);
  size_t ln = 2;
  while (std::getline(in, line)) {
    ++ln;
    if (!line.empty() && line.back() == '\r') line.pop_back();
    for (size_t p = 0; p + 1 <= line.size(); p += 12) {
      const std::string f = trim(line.substr(p, 12));
      if (f.empty()) continue;
      try { v.push_back(std::stod(f)); } catch (...) { throw ReadError(path + ":" + std::to_string(ln) + ": '" + f + "' is not a number"); }
    }
  }
  const size_t n3 = 3 * natom;
  if (v.size() < n3) throw ReadError(path + ": " + std::to_string(v.size() / 3) + " coordinates for " + std::to_string(natom) + " atoms");
  const bool vel = v.size() >= 2 * n3, box = v.size() == n3 + 6 || v.size() == 2 * n3 + 6;
  if (!(v.size() == n3 || v.size() == n3 + 6 || v.size() == 2 * n3 || v.size() == 2 * n3 + 6))
    throw ReadError(path + ": " + std::to_string(v.size()) + " numbers after the header fit neither coordinates, velocities nor a box for " + std::to_string(natom) + " atoms");
  c.positions.resize(natom);
  for (size_t i = 0; i < natom; ++i) c.positions[i] = {v[3 * i], v[3 * i + 1], v[3 * i + 2]};
  if (vel) {
    c.velocities.resize(natom);
    for (size_t i = 0; i < natom; ++i) c.velocities[i] = Vec3{v[n3 + 3 * i], v[n3 + 3 * i + 1], v[n3 + 3 * i + 2]} * (kVelUnit / 1000.0);
  }
  if (box) {
    const double* b = &v[v.size() - 6];
    c.cell = box_cell(b[0], b[1], b[2], b[3], b[4], b[5]);
    c.has_box = true;
  }
  return c;
}

// ---------------------------------------------------------------------------------------------------------------- writer

namespace {

class PrmtopWriter {
 public:
  explicit PrmtopWriter(std::ostream& o) : o_(o) {}
  void flag(const std::string& name, const std::string& fmt) { o_ << "%FLAG " << name << "\n%FORMAT(" << fmt << ")\n"; }
  void ints(const std::string& name, const std::vector<long>& v) {
    flag(name, "10I8");
    for (size_t k = 0; k < v.size(); ++k) {
      char b[16];
      std::snprintf(b, sizeof b, "%8ld", v[k]);
      o_ << b << ((k + 1) % 10 == 0 || k + 1 == v.size() ? "\n" : "");
    }
    if (v.empty()) o_ << "\n";
  }
  void reals(const std::string& name, const std::vector<double>& v) {
    flag(name, "5E16.8");
    for (size_t k = 0; k < v.size(); ++k) {
      char b[32];
      std::snprintf(b, sizeof b, "%16.8E", v[k]);
      o_ << b << ((k + 1) % 5 == 0 || k + 1 == v.size() ? "\n" : "");
    }
    if (v.empty()) o_ << "\n";
  }
  void texts(const std::string& name, const std::vector<std::string>& v) {
    flag(name, "20a4");
    for (size_t k = 0; k < v.size(); ++k) {
      std::string t = v[k].substr(0, 4);
      t.resize(4, ' ');
      o_ << t << ((k + 1) % 20 == 0 || k + 1 == v.size() ? "\n" : "");
    }
    if (v.empty()) o_ << "\n";
  }

 private:
  std::ostream& o_;
};

}  // namespace

std::vector<std::string> write_amber(const System& s, const ForceField& ff, const std::string& stem) {
  std::vector<std::string> notes;
  const size_t n = s.atoms.size();
  auto refuse = [&](const std::string& why) { throw FieldError("AMBER prmtop: " + why); };
  if (ff.charge.size() != n) refuse("the force field is for another structure");
  if (ff.pair_form != "lj12-6") refuse("the " + ff.pair_form + " pair form (class II 9-6) has no AMBER form — AMBER holds 12-6 Lennard-Jones only");
  if (!ff.bonds2.empty() || !ff.angles2.empty() || !ff.dihedrals2.empty() || !ff.impropers2.empty()) refuse("class II cross terms have no AMBER form");
  if (!ff.pair_func.empty()) refuse("Buckingham / Morse pairs have no AMBER form");
  if (!ff.lj14_types.empty()) refuse("separate 1-4 Lennard-Jones parameters (CHARMM, GROMOS) need a chamber topology; export GROMACS or LAMMPS instead");
  if (!ff.lj_pairs.empty()) refuse("explicit Lennard-Jones atom pairs have no AMBER form");
  if (!ff.urey_bradley.empty()) refuse("Urey–Bradley terms need a chamber topology");
  if (!ff.impropers_harmonic.empty()) refuse("harmonic impropers need a chamber topology");
  if (!ff.inversions.empty()) refuse("inversion (umbrella) impropers have no AMBER form");
  if (!ff.bonds_x.empty() || !ff.angles_x.empty() || !ff.cbt.empty()) refuse("Morse / cosine / bending–torsion terms have no AMBER form");
  if (!ff.vsites.empty()) refuse("virtual sites (extra points) are not written");
  if (ff.sw.on || ff.manybody.on() || ff.hbond.on()) refuse("many-body or hydrogen-bond terms have no AMBER form");
  if (ff.keep13 || !ff.excluded_type_pairs.empty()) refuse("1-3 pairs in full or excluded type pairs have no AMBER form");
  if (ff.dielectric != 1 || ff.coul_gromacs || ff.coul_rf) refuse("a relative permittivity or reaction-field Coulomb has no AMBER topology form");
  if ((ff.lj14 == 0) != (ff.coul14 == 0)) refuse("1-4 scaling with LJ or Coulomb alone at zero has no AMBER form");
  const bool no14 = ff.lj14 == 0 && ff.coul14 == 0;
  if (ff.lj_shift || ff.lj_fsw) notes.push_back("the force field's own Lennard-Jones shift / switch is not in a topology: set the engine's cut-off treatment to match");

  // types used, in order of first use; names of at most 4 characters (the file's field), else T1, T2 …
  std::vector<int> tmap(ff.type_names.size(), -1), tused;
  for (size_t i = 0; i < n; ++i)
    if (tmap[size_t(ff.type_index[i])] < 0) tmap[size_t(ff.type_index[i])] = int(tused.size()), tused.push_back(ff.type_index[i]);
  const size_t nt = tused.size();
  std::vector<std::string> tname(nt);
  {
    std::set<std::string> seen;
    bool fit = true;
    for (size_t t = 0; t < nt; ++t) {
      const std::string& nm = ff.type_names[size_t(tused[t])];
      if (nm.empty() || nm.size() > 4 || nm.find(' ') != std::string::npos || !seen.insert(nm).second) fit = false;
    }
    for (size_t t = 0; t < nt; ++t) tname[t] = fit ? ff.type_names[size_t(tused[t])] : "T" + std::to_string(t + 1);
    if (!fit) {
      std::string m;
      for (size_t t = 0; t < nt && t < 40; ++t) m += (m.empty() ? "" : ", ") + tname[t] + " = " + ff.type_names[size_t(tused[t])];
      notes.push_back("AMBER atom types hold 4 characters: the types are written as " + m + (nt > 40 ? " …" : ""));
    }
  }
  // Lennard-Jones A/B for every pair of used types, as the force field mixes them
  std::vector<double> A, B;
  std::vector<long> nbidx(nt * nt);
  for (size_t j = 0; j < nt; ++j)
    for (size_t i = 0; i <= j; ++i) {
      const PairType p = mixed_pair(ff, tused[i], tused[j]);
      const double s6 = std::pow(p.sigma, 6);
      A.push_back(4 * p.eps * s6 * s6);
      B.push_back(4 * p.eps * s6);
      const long ix = long(A.size());
      nbidx[i * nt + j] = ix, nbidx[j * nt + i] = ix;
    }

  // residues: runs of atoms with the same molecule, residue number and name
  std::vector<std::string> reslab;
  std::vector<long> resptr;
  std::vector<int> resof(n);
  for (size_t i = 0; i < n; ++i) {
    const Atom& a = s.atoms[i];
    const bool fresh = i == 0 || a.mol != s.atoms[i - 1].mol || a.resid != s.atoms[i - 1].resid || a.resname != s.atoms[i - 1].resname;
    if (fresh) {
      reslab.push_back(a.resname.empty() ? "MOL" : a.resname);
      resptr.push_back(long(i + 1));
    }
    resof[i] = int(reslab.size() - 1);
  }
  long nmxrs = 0;
  for (size_t r = 0; r < resptr.size(); ++r)
    nmxrs = std::max(nmxrs, (r + 1 < resptr.size() ? resptr[r + 1] : long(n) + 1) - resptr[r]);
  // atom names: the file's when short, else element and a count within the residue
  std::vector<std::string> names(n);
  {
    std::map<std::pair<int, std::string>, int> count;
    for (size_t i = 0; i < n; ++i) {
      const Atom& a = s.atoms[i];
      if (!a.name.empty() && a.name.size() <= 4 && a.name.find(' ') == std::string::npos && a.name != ff.atom_type[i]) { names[i] = a.name; continue; }
      const std::string el = a.element > 0 ? element(a.element).symbol : "X";
      const int k = ++count[{resof[i], el}];
      names[i] = (el + std::to_string(k)).substr(0, 4);
    }
  }

  // bonded terms by type (exact values), split by whether a hydrogen takes part
  auto is_h = [&](uint32_t a) { return s.atoms[a].element == 1; };
  auto bits = [](double x) { uint64_t u; std::memcpy(&u, &x, 8); return u; };
  std::map<std::pair<uint64_t, uint64_t>, long> btype;
  std::vector<double> bk, br;
  std::vector<long> bh, bn;
  for (const auto& b : ff.bonds) {
    auto [it, fresh] = btype.emplace(std::make_pair(bits(b.k), bits(b.r0)), long(bk.size() + 1));
    if (fresh) bk.push_back(b.k), br.push_back(b.r0);
    auto& L = is_h(b.i) || is_h(b.j) ? bh : bn;
    L.insert(L.end(), {long(3 * b.i), long(3 * b.j), it->second});
  }
  std::map<std::pair<uint64_t, uint64_t>, long> atype;
  std::vector<double> ak, at;
  std::vector<long> ah, an;
  for (const auto& a : ff.angles) {
    auto [it, fresh] = atype.emplace(std::make_pair(bits(a.kt), bits(a.theta0)), long(ak.size() + 1));
    if (fresh) ak.push_back(a.kt), at.push_back(a.theta0);
    auto& L = is_h(a.i) || is_h(a.j) || is_h(a.k) ? ah : an;
    L.insert(L.end(), {long(3 * a.i), long(3 * a.j), long(3 * a.k), it->second});
  }
  // torsions: the first term reaching a 1-4 pair carries it; the other terms of the pair, and impropers, do not
  std::set<std::pair<uint32_t, uint32_t>> p14;
  for (const auto& p : ff.pairs14) p14.insert({std::min(p[0], p[1]), std::max(p[0], p[1])});
  if (no14 && !p14.empty()) notes.push_back("1-4 scaling of zero: the 1-4 pairs are excluded outright");
  std::set<std::pair<uint32_t, uint32_t>> carried;
  std::map<std::tuple<uint64_t, int, uint64_t>, long> dtype;
  std::vector<double> dk, dn, dp;
  std::vector<long> dh, dnn;
  double constant = 0;
  int dropped0 = 0;
  auto torsion = [&](TorsionTerm t, bool improper) {
    if (t.n == 0) { constant += t.v * (1 + std::cos(t.delta)); ++dropped0; return; }
    if (t.n < 0) refuse("a torsion with negative periodicity");
    const auto key = std::make_pair(std::min(t.i, t.l), std::max(t.i, t.l));
    bool carry = false;
    if (!improper && !no14 && p14.count(key) && !carried.count(key)) carry = true, carried.insert(key);
    auto [it, fresh] = dtype.emplace(std::make_tuple(bits(t.v), t.n, bits(t.delta)), long(dk.size() + 1));
    if (fresh) dk.push_back(t.v), dn.push_back(double(t.n)), dp.push_back(t.delta);
    uint32_t a = t.i, b = t.j, c = t.k, d = t.l;
    if (c == 0 || d == 0) std::swap(a, d), std::swap(b, c);   // a negative flag cannot mark atom 0: the same angle read backwards
    const long lc = long(3 * c) * (carry ? 1 : -1), ld = long(3 * d) * (improper ? -1 : 1);
    auto& L = is_h(a) || is_h(b) || is_h(c) || is_h(d) ? dh : dnn;
    L.insert(L.end(), {long(3 * a), long(3 * b), lc, ld, it->second});
  };
  for (const auto& t : ff.dihedrals) torsion(t, false);
  for (const auto& t : ff.impropers) torsion(t, true);
  if (dropped0) {
    char b[160];
    std::snprintf(b, sizeof b, "%d constant torsion terms (n = 0) left out: the energy is lower by %.6f kcal/mol, no force changes", dropped0, constant);
    notes.push_back(b);
  }
  // 1-4 pairs no torsion term reaches: a zero-barrier term along the bonds carries each
  if (!no14) {
    const auto nb = s.neighbours();
    int added = 0;
    for (const auto& key : p14) {
      if (carried.count(key)) continue;
      const uint32_t i = key.first, l = key.second;
      bool found = false;
      for (uint32_t j : nb[i]) {
        for (uint32_t k : nb[j])
          if (k != i && std::find(nb[k].begin(), nb[k].end(), l) != nb[k].end() && l != j) {
            torsion(TorsionTerm{i, j, k, l, 0.0, 1, 0.0}, false);
            found = true;
            break;
          }
        if (found) break;
      }
      if (!found) refuse("a 1-4 pair (atoms " + std::to_string(i + 1) + ", " + std::to_string(l + 1) + ") three bonds apart along no path of the structure's bonds");
      ++added;
    }
    if (added) notes.push_back(std::to_string(added) + " 1-4 pairs carried by zero-barrier torsion terms (no torsion of the force field reaches them)");
  }

  // exclusions: each atom's higher-numbered partners (1-2, 1-3, 1-4), a 0 for an atom with none
  std::vector<long> nex(n), exl;
  for (size_t i = 0; i < n; ++i) {
    std::vector<uint32_t> e;
    if (i < ff.excluded.size())
      for (uint32_t j : ff.excluded[i])
        if (j > i) e.push_back(j);
    for (const auto& [a, b] : p14)   // the 1-4 pairs are excluded from the ordinary pairs too
      if (a == i && std::find(e.begin(), e.end(), b) == e.end()) e.push_back(b);
    std::sort(e.begin(), e.end());
    if (e.empty()) { nex[i] = 1; exl.push_back(0); continue; }
    nex[i] = long(e.size());
    for (uint32_t j : e) exl.push_back(long(j + 1));
  }

  // molecules (a periodic topology lists them in order: each must be a run of atoms)
  const bool box = s.cell.valid();
  std::vector<long> apm;
  if (box) {
    int nm = 0;
    const auto m = s.molecules(&nm);
    for (size_t i = 0; i < n; ++i) {
      if (i == 0 || m[i] != m[i - 1]) {
        if (m[i] != int(apm.size())) refuse("the molecules are not runs of atoms in order (AMBER's molecule list needs them so): reorder the atoms by molecule first");
        apm.push_back(0);
      }
      ++apm.back();
    }
  }
  auto ang = [](const Vec3& u, const Vec3& v) { return std::acos(std::clamp(dot(u, v) / (norm(u) * norm(v)), -1.0, 1.0)) * 180 / kPi; };
  const double la = norm(s.cell.a), lb = norm(s.cell.b), lc = norm(s.cell.c);
  const double al = box ? ang(s.cell.b, s.cell.c) : 90, be = box ? ang(s.cell.a, s.cell.c) : 90, ga = box ? ang(s.cell.a, s.cell.b) : 90;
  const bool octahedron = box && std::fabs(al - 109.4712206) < 1e-4 && std::fabs(be - 109.4712206) < 1e-4 && std::fabs(ga - 109.4712206) < 1e-4;

  std::ofstream o(stem + ".prmtop");
  if (!o) throw std::runtime_error("cannot write " + stem + ".prmtop");
  PrmtopWriter w(o);
  o << "%VERSION  VERSION_STAMP = V0001.000  DATE = 00/00/00  00:00:00  (CAPS)\n";
  w.flag("TITLE", "20a4");   // one 80-column line
  {
    std::string t = (s.title.empty() ? std::string("CAPS ") + ff.name : s.title).substr(0, 80);
    t.resize(80, ' ');
    o << t << "\n";
  }
  const long nbonh = long(bh.size() / 3), mbona = long(bn.size() / 3), ntheth = long(ah.size() / 4), mtheta = long(an.size() / 4);
  const long nphih = long(dh.size() / 5), mphia = long(dnn.size() / 5);
  w.ints("POINTERS", {long(n), long(nt), nbonh, mbona, ntheth, mtheta, nphih, mphia, 0, 0, long(exl.size()), long(reslab.size()), mbona, mtheta, mphia,
                      long(bk.size()), long(ak.size()), long(dk.size()), long(nt), 0, 0, 0, 0, 0, 0, 0, 0, box ? (octahedron ? 2 : 1) : 0, nmxrs, 0, 0, 0});
  w.texts("ATOM_NAME", names);
  {
    std::vector<double> q(n);
    for (size_t i = 0; i < n; ++i) q[i] = ff.charge[i] * kChargeUnit;
    w.reals("CHARGE", q);
  }
  {
    std::vector<long> z(n);
    for (size_t i = 0; i < n; ++i) z[i] = s.atoms[i].element > 0 ? s.atoms[i].element : -1;
    w.ints("ATOMIC_NUMBER", z);
  }
  w.reals("MASS", ff.mass);
  {
    std::vector<long> ti(n);
    for (size_t i = 0; i < n; ++i) ti[i] = tmap[size_t(ff.type_index[i])] + 1;
    w.ints("ATOM_TYPE_INDEX", ti);
  }
  w.ints("NUMBER_EXCLUDED_ATOMS", nex);
  w.ints("NONBONDED_PARM_INDEX", nbidx);
  w.texts("RESIDUE_LABEL", reslab);
  w.ints("RESIDUE_POINTER", resptr);
  w.reals("BOND_FORCE_CONSTANT", bk);
  w.reals("BOND_EQUIL_VALUE", br);
  w.reals("ANGLE_FORCE_CONSTANT", ak);
  w.reals("ANGLE_EQUIL_VALUE", at);
  w.reals("DIHEDRAL_FORCE_CONSTANT", dk);
  w.reals("DIHEDRAL_PERIODICITY", dn);
  w.reals("DIHEDRAL_PHASE", dp);
  w.reals("SCEE_SCALE_FACTOR", std::vector<double>(dk.size(), no14 ? 1.2 : 1 / ff.coul14));
  w.reals("SCNB_SCALE_FACTOR", std::vector<double>(dk.size(), no14 ? 2.0 : 1 / ff.lj14));
  w.reals("SOLTY", std::vector<double>(nt, 0.0));
  w.reals("LENNARD_JONES_ACOEF", A);
  w.reals("LENNARD_JONES_BCOEF", B);
  w.ints("BONDS_INC_HYDROGEN", bh);
  w.ints("BONDS_WITHOUT_HYDROGEN", bn);
  w.ints("ANGLES_INC_HYDROGEN", ah);
  w.ints("ANGLES_WITHOUT_HYDROGEN", an);
  w.ints("DIHEDRALS_INC_HYDROGEN", dh);
  w.ints("DIHEDRALS_WITHOUT_HYDROGEN", dnn);
  w.ints("EXCLUDED_ATOMS_LIST", exl);
  w.reals("HBOND_ACOEF", {});
  w.reals("HBOND_BCOEF", {});
  w.reals("HBCUT", {});
  {
    std::vector<std::string> ty(n);
    for (size_t i = 0; i < n; ++i) ty[i] = tname[size_t(tmap[size_t(ff.type_index[i])])];
    w.texts("AMBER_ATOM_TYPE", ty);
  }
  w.texts("TREE_CHAIN_CLASSIFICATION", std::vector<std::string>(n, "BLA"));
  w.ints("JOIN_ARRAY", std::vector<long>(n, 0));
  w.ints("IROTAT", std::vector<long>(n, 0));
  if (box) {
    w.ints("SOLVENT_POINTERS", {long(reslab.size()), long(apm.size()), long(apm.size()) + 1});
    w.ints("ATOMS_PER_MOLECULE", apm);
    w.reals("BOX_DIMENSIONS", {be, la, lb, lc});
  }
  w.ints("IPOL", {0});
  if (!o) throw std::runtime_error("could not write " + stem + ".prmtop");
  o.close();

  // the restart: 6F12.7, the box as lengths and angles
  std::ofstream c(stem + ".inpcrd");
  if (!c) throw std::runtime_error("cannot write " + stem + ".inpcrd");
  {
    std::string t = (s.title.empty() ? std::string("CAPS") : s.title).substr(0, 80);
    c << t << "\n";
    char b[64];
    std::snprintf(b, sizeof b, "%6zu\n", n);
    c << b;
    size_t k = 0;
    auto put = [&](double x) {
      if (std::fabs(x) >= 9999.99999995) refuse("a coordinate beyond ±9999.9999999 Å does not fit the restart's F12.7");
      std::snprintf(b, sizeof b, "%12.7f", x);
      c << b << (++k % 6 == 0 ? "\n" : "");
    };
    for (const auto& a : s.atoms) put(a.pos[0]), put(a.pos[1]), put(a.pos[2]);
    if (k % 6) c << "\n";
    if (box) {
      k = 0;
      for (double x : {la, lb, lc, al, be, ga}) put(x);
    }
  }
  if (!c) throw std::runtime_error("could not write " + stem + ".inpcrd");
  if (box && (std::fabs(al - 90) > 1e-6 || std::fabs(ga - 90) > 1e-6) && !octahedron)
    notes.push_back("a triclinic cell: its three angles are in the restart (the prmtop's BOX_DIMENSIONS holds β only)");
  notes.push_back("no generalised-Born radii: set them (ParmEd changeRadii) before an implicit-solvent run");
  return notes;
}

Trajectory read_amber_netcdf(const std::string& path, const System& topology, size_t max_frames,
                             const std::function<bool(double, const Trajectory&)>& progress) {
  NcFile f;
  open_netcdf(f, path);
  const std::string conv = nc_convention(f);
  if (conv.find("AMBERRESTART") != std::string::npos) {   // one frame: a restart
    const AmberCoordinates c = read_amber_coordinates(path);
    if (c.positions.size() != topology.atoms.size())
      throw ReadError(path + " has " + std::to_string(c.positions.size()) + " atoms, the topology " + std::to_string(topology.atoms.size()));
    Trajectory t;
    t.topology = topology;
    t.positions.push_back(c.positions);
    t.cells.push_back(c.has_box ? c.cell : topology.cell);
    t.timesteps.push_back(0);
    if (!c.velocities.empty()) t.topology.velocities = c.velocities;
    t.topology.cell = t.cells.front();
    return t;
  }
  if (conv.find("AMBER") == std::string::npos) throw ReadError(path + ": a NetCDF file without the AMBER convention (Conventions: '" + conv + "')");
  const NcVar* xv = f.var("coordinates");
  if (!xv || !xv->record) throw ReadError(path + ": no coordinates per frame");
  const size_t n = size_t(f.dim("atom"));
  if (n != topology.atoms.size())
    throw ReadError(path + " has " + std::to_string(n) + " atoms, the topology " + std::to_string(topology.atoms.size()));
  const NcVar *lv = f.var("cell_lengths"), *av = f.var("cell_angles"), *tv = f.var("time"), *vv = f.var("velocities");
  Trajectory t;
  t.topology = topology;
  std::vector<double> times;
  for (uint64_t r = 0; r < f.numrecs; ++r) {
    const auto x = nc_read(f, *xv, r);
    std::vector<Vec3> p(n);
    for (size_t i = 0; i < n; ++i) p[i] = {x[3 * i], x[3 * i + 1], x[3 * i + 2]};
    t.positions.push_back(std::move(p));
    Cell cell = topology.cell;
    if (lv && av) {
      const auto l = nc_read(f, *lv, lv->record ? r : 0), a = nc_read(f, *av, av->record ? r : 0);
      if (l.size() == 3 && a.size() == 3 && l[0] > 0) cell = box_cell(l[0], l[1], l[2], a[0], a[1], a[2]);
    }
    t.cells.push_back(cell);
    t.timesteps.push_back(int64_t(r));
    if (tv && tv->record) times.push_back(nc_read(f, *tv, r).at(0));
    if (r == 0 && vv && vv->record) {   // the first frame's velocities (Å/ps after the file's scale factor) as the structure's
      const auto v = nc_read(f, *vv, 0);
      t.topology.velocities.clear();
      for (size_t i = 0; i < n && v.size() == 3 * n; ++i) t.topology.velocities.push_back(Vec3{v[3 * i], v[3 * i + 1], v[3 * i + 2]} * 1e-3);
    }
    if (progress && !progress(double(r + 1) / double(std::max<uint64_t>(1, f.numrecs)), t)) {
      t.topology.notes.push_back("reading stopped after " + std::to_string(t.frames()) + " frames");
      break;
    }
    if (max_frames && t.frames() >= max_frames) break;
  }
  if (t.frames() == 0) throw ReadError(path + ": no frames");
  t.topology.cell = t.cells.front();
  t.topology.unwrapped = false;
  std::string note = std::to_string(t.frames()) + " frames from " + std::filesystem::path(path).filename().string() + " (AMBER NetCDF)";
  if (times.size() > 1) {
    char b[96];
    std::snprintf(b, sizeof b, ", %.4g to %.4g ps", times.front(), times.back());
    note += b;
  }
  t.topology.notes.push_back(note);
  return t;
}

Trajectory read_amber_mdcrd(const std::string& path, const System& topology, size_t max_frames,
                            const std::function<bool(double, const Trajectory&)>& progress) {
  std::ifstream in(path);
  if (!in) throw ReadError("cannot open " + path);
  const size_t n = topology.atoms.size(), n3 = 3 * n;
  if (n == 0) throw ReadError(path + ": the topology has no atoms");
  const bool box = topology.cell.valid();   // AMBER writes a box line after each frame of a periodic run
  double al = 90, be = 90, ga = 90;
  if (box) {
    auto ang = [](const Vec3& u, const Vec3& v) { return std::acos(std::clamp(dot(u, v) / (norm(u) * norm(v)), -1.0, 1.0)) * 180 / kPi; };
    al = ang(topology.cell.b, topology.cell.c), be = ang(topology.cell.a, topology.cell.c), ga = ang(topology.cell.a, topology.cell.b);
  }
  const auto size = double(std::filesystem::file_size(path));
  std::string line;
  std::getline(in, line);   // title
  Trajectory t;
  t.topology = topology;
  std::vector<double> v;
  v.reserve(n3);
  size_t ln = 1;
  auto fields = [&](const std::string& l, std::vector<double>& out) {
    for (size_t p = 0; p < l.size(); p += 8) {
      const std::string f = trim(l.substr(p, 8));
      if (f.empty()) continue;
      if (f.find('*') != std::string::npos) throw ReadError(path + ":" + std::to_string(ln) + ": a coordinate overflowed its 8 columns (********)");
      try { out.push_back(std::stod(f)); } catch (...) { throw ReadError(path + ":" + std::to_string(ln) + ": '" + f + "' is not a number"); }
    }
  };
  while (std::getline(in, line)) {
    ++ln;
    if (!line.empty() && line.back() == '\r') line.pop_back();
    if (trim(line).empty()) continue;
    fields(line, v);
    if (v.size() < n3) continue;
    if (v.size() > n3) throw ReadError(path + ":" + std::to_string(ln) + ": a frame's coordinates run past " + std::to_string(n) + " atoms (another topology?)");
    std::vector<Vec3> p(n);
    for (size_t i = 0; i < n; ++i) p[i] = {v[3 * i], v[3 * i + 1], v[3 * i + 2]};
    v.clear();
    Cell cell = topology.cell;
    if (box) {
      if (!std::getline(in, line)) throw ReadError(path + ": the last frame has no box line (the topology is periodic)");
      ++ln;
      std::vector<double> b;
      {   // 3F8.3 as AMBER writes it; some writers separate wider columns by spaces (mdtraj: %8.3f %8.3f %8.3f)
        std::istringstream ws(line);
        for (double x; ws >> x;) b.push_back(x);
        if (!ws.eof() || (b.size() != 3 && b.size() != 6)) b.clear(), fields(line, b);
      }
      if (b.size() != 3 && b.size() != 6) throw ReadError(path + ":" + std::to_string(ln) + ": a box line of 3 lengths expected after the frame, found " + std::to_string(b.size()) + " values");
      cell = b.size() == 6 ? box_cell(b[0], b[1], b[2], b[3], b[4], b[5]) : box_cell(b[0], b[1], b[2], al, be, ga);
    }
    t.positions.push_back(std::move(p));
    t.cells.push_back(cell);
    t.timesteps.push_back(int64_t(t.frames() - 1));
    if (progress && !progress(size > 0 ? double(in.tellg()) / size : 1.0, t)) {
      t.topology.notes.push_back("reading stopped after " + std::to_string(t.frames()) + " frames");
      break;
    }
    if (max_frames && t.frames() >= max_frames) break;
  }
  if (!v.empty()) throw ReadError(path + ": the file ends inside a frame (" + std::to_string(v.size() / 3) + " of " + std::to_string(n) + " atoms)");
  if (t.frames() == 0) throw ReadError(path + ": no frames");
  t.topology.cell = t.cells.front();
  t.topology.unwrapped = false;
  t.topology.notes.push_back(std::to_string(t.frames()) + " frames from " + std::filesystem::path(path).filename().string() + " (AMBER mdcrd, 0.001 Å)" +
                             (box ? ", a box per frame" : ""));
  return t;
}

}  // namespace caps
