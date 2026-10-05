// LAMMPS data (atom_style full / molecular / charge / atomic) and text dump readers, data writer.
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <limits>
#include <numeric>
#include <unordered_map>

#include "caps/analysis.hpp"
#include "caps/elements.hpp"
#include "caps/io.hpp"
#include "io_util.hpp"

namespace caps {

std::string residue_comment(const Atom& a) {
  if (a.resid <= 0) return "";
  std::string name = a.resname;
  for (auto& c : name) if (std::isspace(static_cast<unsigned char>(c))) c = '_';
  return "  # res " + std::to_string(a.resid) + (name.empty() ? "" : " " + name);
}

void read_residue_comment(const std::string& comment, Atom& a) {
  std::istringstream is(comment);
  std::string w;
  is >> w;
  if (w != "res") return;
  long long r = 0;
  if (!(is >> r) || r <= 0) return;
  a.resid = r;
  std::string name;
  if (is >> name) a.resname = name;
}

namespace {

Cell lammps_cell(double xlo, double xhi, double ylo, double yhi, double zlo, double zhi, double xy, double xz, double yz) {
  Cell c;
  c.origin = {xlo, ylo, zlo};
  c.a = {xhi - xlo, 0, 0};
  c.b = {xy, yhi - ylo, 0};
  c.c = {xz, yz, zhi - zlo};
  return c;
}

void assign_elements(System& s) {
  int guessed = 0, sites = 0;
  for (auto& at : s.atoms) {
    if (at.element) continue;
    for (const auto& t : s.types)
      if (t.type == at.type) {
        // A label such as "c3" or "CT" names the element first; the mass decides when it does not.
        int z = element_from_mass(t.mass);
        int ua_h = 0;   // hydrogens a united-atom site holds (named for them: CH2, so typing knows its hydrogens)
        if (!z && !t.label.empty()) z = element_from_name(t.label);
        // a labelled united-atom site whose mass falls near another element ("CH2" 14.027 reads as N by mass): the
        // label's element when the mass is that element's plus 1–4 hydrogens
        // without a label: a united-atom reading clearly closer than the plain element (CH2 14.027 against N 14.007)
        if (z && (z == 7 || z == 8 || z == 16 || z == 9) && t.label.empty()) {
          const double d0 = std::fabs(element(z).mass - t.mass);
          for (int host : {6, 7, 8, 16})
            for (int nh = 1; nh <= 4; ++nh)
              if (std::fabs(element(host).mass + nh * 1.008 - t.mass) < d0 - 0.003) { z = host, ua_h = nh, ++sites; nh = 5; host = 99; }
        }
        if (z && !t.label.empty()) {
          const int zl = element_from_name(t.label);
          if (zl > 0 && zl != z && (zl == 6 || zl == 7 || zl == 8 || zl == 16))
            for (int nh = 1; nh <= 4; ++nh)
              if (std::fabs(element(zl).mass + nh * 1.008 - t.mass) < 0.02) { z = zl, ua_h = nh, ++sites; break; }
        }
        // a united-atom site (mW water's 18.02, a CH2's 14.03 under another label): its heavy atom plus hydrogens
        // (C, N, O or S hosts only: a coarse-grained bead's 54 or 72 must not read as a metal hydride)
        for (int nh = 1; !z && nh <= 4 && t.mass > 1.5 + nh * 1.008; ++nh) {
          const int h = element_from_mass(t.mass - nh * 1.008, 0.02);
          if (h == 6 || h == 7 || h == 8 || h == 16) z = h, ua_h = nh, ++sites;
        }
        at.element = z;
        if (!at.name.size() && !t.label.empty()) at.name = t.label;
        if (!at.name.size() && ua_h > 0) at.name = std::string(element(z).symbol) + "H" + (ua_h > 1 ? std::to_string(ua_h) : "");
        break;
      }
    if (at.element) ++guessed;
  }
  if (guessed) s.notes.push_back("elements guessed from type masses for " + std::to_string(guessed) + " atoms");
  if (sites) s.notes.push_back(std::to_string(sites) + " united-atom sites: element from the mass less its hydrogens");
  // Bead-spring models (Kremer–Grest, reduced units) give beads mass 1, which reads as hydrogen. A type whose atoms
  // bond to two or more others is beads, not hydrogens: they are drawn as carbon, as CAPS's own bead-spring melts are.
  std::vector<int> nbond(s.atoms.size(), 0);
  for (const auto& b : s.bonds) ++nbond[b.i], ++nbond[b.j];
  std::set<int> bead_types;
  for (size_t i = 0; i < s.atoms.size(); ++i)
    if (s.atoms[i].element == 1 && nbond[i] >= 2) bead_types.insert(s.atoms[i].type);
  if (!bead_types.empty()) {
    size_t nb = 0;
    for (auto& at : s.atoms)
      if (at.element == 1 && bead_types.count(at.type)) at.element = 6, ++nb;
    s.notes.push_back("bead-spring model: " + std::to_string(nb) + " mass-1 beads bonded to two or more others are beads, not hydrogens (drawn as carbon)");
  }
}

}  // namespace

System read_lammps_data(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw ReadError("cannot open " + path);
  System s;
  s.source_format = "lammps-data";
  std::string line;
  std::getline(in, line);
  s.title = trim(line);

  size_t natoms = 0, nbonds = 0;
  bool bonds_declared = false;   // a "N bonds" header line: the file says what is bonded, even when N is 0
  double xlo = 0, xhi = 0, ylo = 0, yhi = 0, zlo = 0, zhi = 0, xy = 0, xz = 0, yz = 0;
  std::string section, style_hint;
  std::unordered_map<int64_t, uint32_t> index;
  std::vector<std::pair<int64_t, Vec3>> vel;
  size_t lineno = 1;
  size_t skipped_bonds = 0;

  while (std::getline(in, line)) {
    ++lineno;
    std::string comment;
    const std::string body = trim(strip_comment(line, &comment));
    if (body.empty()) continue;
    auto tok = split(body);

    // Section headers are a single keyword (optionally followed by a style comment).
    static const char* kSections[] = {"Masses", "Atoms", "Velocities", "Bonds", "Angles", "Dihedrals", "Impropers",
                                      "Pair Coeffs", "Bond Coeffs", "Angle Coeffs", "Dihedral Coeffs", "Improper Coeffs",
                                      "PairIJ Coeffs", "Ellipsoids", "Lines", "Triangles", "Bodies"};
    bool is_header = false;
    for (const char* h : kSections)
      if (body == h) { is_header = true; section = h; style_hint = trim(comment); break; }
    if (is_header) continue;

    if (section.empty()) {
      if (tok.size() >= 2 && tok[1] == "atoms") natoms = std::stoull(tok[0]);
      else if (tok.size() >= 2 && tok[1] == "bonds") { nbonds = std::stoull(tok[0]); bonds_declared = true; }
      else if (tok.size() >= 4 && tok[2] == "xlo") { xlo = std::stod(tok[0]); xhi = std::stod(tok[1]); }
      else if (tok.size() >= 4 && tok[2] == "ylo") { ylo = std::stod(tok[0]); yhi = std::stod(tok[1]); }
      else if (tok.size() >= 4 && tok[2] == "zlo") { zlo = std::stod(tok[0]); zhi = std::stod(tok[1]); }
      else if (tok.size() >= 6 && tok[3] == "xy") { xy = std::stod(tok[0]); xz = std::stod(tok[1]); yz = std::stod(tok[2]); }
      continue;
    }

    if (section == "Masses" && tok.size() >= 2) {
      TypeInfo t;
      t.type = std::stoi(tok[0]);
      t.mass = std::stod(tok[1]);
      t.label = trim(comment);
      s.types.push_back(t);
    } else if (section == "Atoms") {
      // Column layout from the style comment, else from the column count.
      std::string st = style_hint;
      const size_t n = tok.size();
      if (st.empty()) {
        if (n == 7 || n == 10) st = "full";
        else if (n == 6 || n == 9) st = "molecular";   // ambiguous with charge; molecular is more common here
        else if (n == 5 || n == 8) st = "atomic";
        else throw ReadError(path + ":" + std::to_string(lineno) + ": cannot infer atom_style from " + std::to_string(n) + " columns");
      }
      Atom a;
      size_t p = 0;
      a.id = std::stoll(tok[p++]);
      if (st == "full" || st == "molecular" || st == "bond" || st == "angle") { a.mol = std::stoll(tok[p++]); s.has_mol = true; }
      a.type = std::stoi(tok[p++]);
      if (st == "full" || st == "charge") { a.charge = std::stod(tok[p++]); s.has_charges = true; }
      if (p + 3 > n) throw ReadError(path + ":" + std::to_string(lineno) + ": too few columns for atom_style " + st);
      a.pos = {std::stod(tok[p]), std::stod(tok[p + 1]), std::stod(tok[p + 2])};
      p += 3;
      if (p + 3 <= n) a.image = {std::stoi(tok[p]), std::stoi(tok[p + 1]), std::stoi(tok[p + 2])};
      read_residue_comment(comment, a);
      index[a.id] = static_cast<uint32_t>(s.atoms.size());
      s.atoms.push_back(a);
    } else if (section == "Velocities" && tok.size() >= 4) {
      vel.push_back({std::stoll(tok[0]), Vec3{std::stod(tok[1]), std::stod(tok[2]), std::stod(tok[3])}});
    } else if (section == "Bonds" && tok.size() >= 4) {
      auto i = index.find(std::stoll(tok[2]));
      auto j = index.find(std::stoll(tok[3]));
      if (i == index.end() || j == index.end()) { ++skipped_bonds; continue; }
      s.bonds.push_back({i->second, j->second});
    }
  }
  if (!vel.empty()) {
    s.velocities.assign(s.atoms.size(), Vec3{0, 0, 0});
    size_t found = 0;
    for (const auto& [id, v] : vel)
      if (auto it = index.find(id); it != index.end()) { s.velocities[it->second] = v; ++found; }
    if (found != s.atoms.size()) s.notes.push_back("velocities given for " + std::to_string(found) + " of " + std::to_string(s.atoms.size()) + " atoms");
  }
  if (natoms && s.atoms.size() != natoms)
    s.notes.push_back("header says " + std::to_string(natoms) + " atoms, read " + std::to_string(s.atoms.size()));
  if (nbonds && s.bonds.size() + skipped_bonds != nbonds)
    s.notes.push_back("header says " + std::to_string(nbonds) + " bonds, read " + std::to_string(s.bonds.size()));
  if (skipped_bonds) s.notes.push_back(std::to_string(skipped_bonds) + " bonds refer to atoms that are not defined; skipped");
  s.bonds_from_file = !s.bonds.empty() || bonds_declared;
  s.cell = lammps_cell(xlo, xhi, ylo, yhi, zlo, zhi, xy, xz, yz);
  // Image flags make molecules whole: store unwrapped positions, as a dump with xu/yu/zu would.
  size_t flagged = 0;
  for (auto& a : s.atoms)
    if (a.image[0] || a.image[1] || a.image[2]) {
      a.pos = a.pos + s.cell.a * a.image[0] + s.cell.b * a.image[1] + s.cell.c * a.image[2];
      ++flagged;
    }
  if (flagged) {
    s.unwrapped = true;
    s.notes.push_back("image flags applied to " + std::to_string(flagged) + " atoms (positions unwrapped)");
  }
  assign_elements(s);
  return s;
}

Trajectory read_lammps_dump(const std::string& path, const System* topology) { return read_lammps_dump(path, topology, 0, {}); }

Trajectory read_lammps_dump(const std::string& path, const System* topology, size_t max_frames, const std::function<bool(double, const Trajectory&)>& progress,
                            const FrameSelection& frames) {
  std::ifstream in(path, std::ios::binary);
  if (!in) throw ReadError("cannot open " + path);
  double file_size = 0;
  if (progress) {
    in.seekg(0, std::ios::end);
    file_size = static_cast<double>(in.tellg());
    in.seekg(0, std::ios::beg);
  }
  Trajectory tr;
  tr.selection = frames;
  std::string line;
  int64_t step = -1;
  size_t natoms = 0;
  Cell cell;
  bool first = true;
  std::vector<std::string> cols;
  size_t lineno = 0;

  while (std::getline(in, line)) {
    ++lineno;
    if (line.rfind("ITEM: TIMESTEP", 0) == 0) {
      std::getline(in, line); ++lineno;
      step = std::stoll(trim(line));
    } else if (line.rfind("ITEM: NUMBER OF ATOMS", 0) == 0) {
      std::getline(in, line); ++lineno;
      natoms = std::stoull(trim(line));
    } else if (line.rfind("ITEM: BOX BOUNDS", 0) == 0) {
      auto flags = split(line.substr(16));
      const bool tric = !flags.empty() && flags[0] == "xy";
      double lo[3], hi[3], tilt[3] = {0, 0, 0};
      for (int k = 0; k < 3; ++k) {
        std::getline(in, line); ++lineno;
        auto t = split(line);
        lo[k] = std::stod(t[0]); hi[k] = std::stod(t[1]);
        if (tric && t.size() >= 3) tilt[k] = std::stod(t[2]);
      }
      if (tric) {  // bounding box -> lo/hi of the parallelepiped (LAMMPS docs, Howto triclinic)
        const double xy = tilt[0], xz = tilt[1], yz = tilt[2];
        lo[0] -= std::min({0.0, xy, xz, xy + xz}); hi[0] -= std::max({0.0, xy, xz, xy + xz});
        lo[1] -= std::min(0.0, yz); hi[1] -= std::max(0.0, yz);
      }
      cell = lammps_cell(lo[0], hi[0], lo[1], hi[1], lo[2], hi[2], tilt[0], tilt[1], tilt[2]);
      const size_t pb = tric ? 3 : 0;
      for (int k = 0; k < 3; ++k) cell.periodic[k] = flags.size() > pb + k ? flags[pb + k][0] == 'p' : true;
    } else if (line.rfind("ITEM: ATOMS", 0) == 0) {
      if (!tr.selection.wants(tr.frames_read)) {   // a frame not kept: its atom lines passed over unread
        for (size_t i = 0; i < natoms; ++i) {
          if (!std::getline(in, line)) throw ReadError(path + ": file ends inside frame at timestep " + std::to_string(step));
          ++lineno;
        }
        if (tr.selection.past(tr.frames_read++)) break;
        continue;
      }
      cols = split(line.substr(11));
      std::map<std::string, int> ci;
      for (size_t k = 0; k < cols.size(); ++k) ci[cols[k]] = static_cast<int>(k);
      auto col = [&](std::initializer_list<const char*> names) {
        for (const char* nm : names) if (auto it = ci.find(nm); it != ci.end()) return it->second;
        return -1;
      };
      const int cid = col({"id"}), cmol = col({"mol"}), ctype = col({"type"}), cel = col({"element"}), cq = col({"q"});
      int cx = col({"xu", "x", "xs", "xsu"}), cy = col({"yu", "y", "ys", "ysu"}), cz = col({"zu", "z", "zs", "zsu"});
      const bool scaled = cx >= 0 && (cols[cx] == "xs" || cols[cx] == "xsu");
      const bool unwrapped = cx >= 0 && (cols[cx] == "xu" || cols[cx] == "xsu");
      const int cix = col({"ix"}), ciy = col({"iy"}), ciz = col({"iz"});
      if (cx < 0 || cy < 0 || cz < 0) throw ReadError(path + ":" + std::to_string(lineno) + ": no position columns (x/xu/xs)");
      const int cvx = col({"vx"}), cvy = col({"vy"}), cvz = col({"vz"});
      const bool vel = cvx >= 0 && cvy >= 0 && cvz >= 0;
      // every other column is kept by its name (numbers only)
      std::vector<int> extra;
      {
        static const std::set<std::string> known = {"id", "mol", "type", "element", "q", "x", "y", "z", "xu", "yu", "zu", "xs", "ys", "zs",
                                                     "xsu", "ysu", "zsu", "ix", "iy", "iz", "vx", "vy", "vz"};
        for (size_t k = 0; k < cols.size(); ++k) if (!known.count(cols[k])) extra.push_back(int(k));
      }
      std::vector<std::vector<float>> xvals(extra.size(), std::vector<float>(natoms, 0.0f));
      std::vector<Vec3> vvals(vel ? natoms : 0);

      std::vector<Atom> atoms(natoms);
      for (size_t i = 0; i < natoms; ++i) {
        if (!std::getline(in, line)) throw ReadError(path + ": file ends inside frame at timestep " + std::to_string(step));
        ++lineno;
        auto t = split(line);
        if (t.size() < cols.size()) throw ReadError(path + ":" + std::to_string(lineno) + ": expected " + std::to_string(cols.size()) + " columns");
        Atom& a = atoms[i];
        a.id = cid >= 0 ? std::stoll(t[cid]) : static_cast<int64_t>(i + 1);
        if (cmol >= 0) a.mol = std::stoll(t[cmol]);
        if (ctype >= 0) a.type = std::stoi(t[ctype]);
        if (cel >= 0) { a.element = element_from_symbol(t[cel]); a.name = t[cel]; }
        if (cq >= 0) a.charge = std::stod(t[cq]);
        Vec3 r{std::stod(t[cx]), std::stod(t[cy]), std::stod(t[cz])};
        if (scaled) r = cell.to_cartesian(r);
        if (cix >= 0 && ciy >= 0 && ciz >= 0) {
          a.image = {std::stoi(t[cix]), std::stoi(t[ciy]), std::stoi(t[ciz])};
          if (!unwrapped) r = r + cell.a * a.image[0] + cell.b * a.image[1] + cell.c * a.image[2];
        }
        a.pos = r;
        if (vel) vvals[i] = {std::stod(t[cvx]), std::stod(t[cvy]), std::stod(t[cvz])};
        for (size_t e = 0; e < extra.size(); ++e) {
          char* end = nullptr;
          const double v = std::strtod(t[size_t(extra[e])].c_str(), &end);
          xvals[e][i] = end && *end == 0 ? float(v) : std::numeric_limits<float>::quiet_NaN();
        }
      }
      // Frames are stored in id order so topology indices line up (the extra columns follow the same order).
      std::vector<size_t> ord(atoms.size());
      std::iota(ord.begin(), ord.end(), size_t(0));
      std::sort(ord.begin(), ord.end(), [&](size_t x, size_t y) { return atoms[x].id < atoms[y].id; });
      {
        std::vector<Atom> sorted(atoms.size());
        for (size_t i = 0; i < ord.size(); ++i) sorted[i] = atoms[ord[i]];
        atoms = std::move(sorted);
      }
      auto reorder = [&](const auto& v) {
        std::decay_t<decltype(v)> o(v.size());
        for (size_t i = 0; i < ord.size(); ++i) o[i] = v[ord[i]];
        return o;
      };
      if (vel) tr.velocities.push_back(reorder(vvals));
      for (size_t e = 0; e < extra.size(); ++e) tr.columns[cols[size_t(extra[e])]].push_back(reorder(xvals[e]));
      {   // magnitudes of vector columns: |f| from fx fy fz, |v| from the velocities
        auto mag = [&](const std::string& n, const std::vector<float>& a, const std::vector<float>& b, const std::vector<float>& c) {
          std::vector<float> m(a.size());
          for (size_t i = 0; i < a.size(); ++i) m[i] = std::sqrt(a[i] * a[i] + b[i] * b[i] + c[i] * c[i]);
          tr.columns[n].push_back(std::move(m));
        };
        if (tr.columns.count("fx") && tr.columns.count("fy") && tr.columns.count("fz"))
          mag("|f|", tr.columns["fx"].back(), tr.columns["fy"].back(), tr.columns["fz"].back());
        if (vel) {
          std::vector<float> m(tr.velocities.back().size());
          for (size_t i = 0; i < m.size(); ++i) m[i] = float(norm(tr.velocities.back()[i]));
          tr.columns["|v|"].push_back(std::move(m));
        }
      }
      if (first) {
        System s;
        if (topology) {
          s = *topology;
          if (s.atoms.size() != atoms.size())
            throw ReadError("dump has " + std::to_string(atoms.size()) + " atoms but the topology has " + std::to_string(s.atoms.size()));
          for (size_t i = 0; i < atoms.size(); ++i) {
            s.atoms[i].pos = atoms[i].pos;
            if (cq >= 0) s.atoms[i].charge = atoms[i].charge;
          }
        } else {
          s.atoms = atoms;
          s.has_mol = cmol >= 0;
          s.has_charges = cq >= 0;
          if (cel < 0) s.notes.push_back("no element column; elements guessed from type numbers are unknown — open with a data file for masses");
        }
        s.source_format = "lammps-dump";
        s.cell = cell;
        s.timestep = step;
        s.unwrapped = unwrapped || cix >= 0;
        tr.topology = std::move(s);
        first = false;
      }
      std::vector<Vec3> pos(atoms.size());
      for (size_t i = 0; i < atoms.size(); ++i) pos[i] = atoms[i].pos;
      tr.positions.push_back(std::move(pos));
      tr.cells.push_back(cell);
      tr.timesteps.push_back(step);
      const bool more = !tr.selection.past(tr.frames_read++);
      if (max_frames && tr.positions.size() >= max_frames) break;
      if (progress && !progress(file_size > 0 ? std::min(1.0, static_cast<double>(in.tellg()) / file_size) : 1.0, tr)) {
        tr.topology.notes.push_back("reading stopped after " + std::to_string(tr.positions.size()) + " frames");
        break;
      }
      if (!more) break;
    }
  }
  if (tr.positions.empty())
    throw ReadError(path + (tr.selection.all() ? ": no frames found" : ": no frames in the selection (the file has " + std::to_string(tr.frames_read) + ")"));
  return tr;
}

void write_lammps_data(const System& s_in, const std::string& path) {
  // atoms without a numeric type (read from a PDB, .xyz or .gro, or packed from such files): every atom gets a type by
  // its label (its type's, else its element or name) and mass, so the file keeps what the atoms are
  System retyped;
  const bool untyped = std::any_of(s_in.atoms.begin(), s_in.atoms.end(), [](const Atom& a) { return a.type <= 0; });
  if (untyped) {
    retyped = s_in;
    std::map<std::pair<std::string, long>, int> tid;
    std::vector<TypeInfo> types;
    for (auto& a : retyped.atoms) {
      std::string label;
      double mass = 0;
      for (const auto& t : s_in.types)
        if (a.type > 0 && t.type == a.type) label = t.label, mass = t.mass;
      if (label.empty()) label = a.element > 0 ? std::string(element(a.element).symbol) : (a.name.empty() ? "X" : a.name);
      if (mass == 0) mass = element(a.element).mass;
      auto [it, fresh] = tid.emplace(std::make_pair(label, std::lround(mass * 1e4)), int(tid.size()) + 1);
      a.type = it->second;
      if (fresh) types.push_back({it->second, mass, label});
    }
    retyped.types = types;
  }
  const System& s = untyped ? retyped : s_in;
  std::ofstream out(path);
  if (!out) throw ReadError("cannot write " + path);
  out << "CAPS 0.1 · " << (s.title.empty() ? "structure" : s.title) << " · atom_style full · units real\n\n";
  out << s.atoms.size() << " atoms\n" << s.bonds.size() << " bonds\n\n";
  int ntypes = 0;
  for (const auto& a : s.atoms) ntypes = std::max(ntypes, a.type);
  out << ntypes << " atom types\n" << (s.bonds.empty() ? 0 : 1) << " bond types\n\n";
  Cell c = s.cell;
  if (!c.valid()) {   // no cell: a box around the atoms with 25 Å of vacuum on every side, at least 100 Å
    Vec3 mn{1e30, 1e30, 1e30}, mx{-1e30, -1e30, -1e30};
    for (const auto& at : s.atoms)
      for (int k = 0; k < 3; ++k) mn[k] = std::min(mn[k], at.pos[k]), mx[k] = std::max(mx[k], at.pos[k]);
    if (s.atoms.empty()) mn = {0, 0, 0}, mx = {0, 0, 0};
    Vec3 w;
    for (int k = 0; k < 3; ++k) w[k] = std::max(100.0, mx[k] - mn[k] + 50.0), c.origin[k] = 0.5 * (mn[k] + mx[k]) - 0.5 * w[k];
    c.a = {w[0], 0, 0}, c.b = {0, w[1], 0}, c.c = {0, 0, w[2]};
  }
  char buf[256];
  std::snprintf(buf, sizeof buf, "%.10f %.10f xlo xhi\n%.10f %.10f ylo yhi\n%.10f %.10f zlo zhi\n", c.origin[0], c.origin[0] + c.a[0],
                c.origin[1], c.origin[1] + c.b[1], c.origin[2], c.origin[2] + c.c[2]);
  out << buf;
  if (std::fabs(c.b[0]) + std::fabs(c.c[0]) + std::fabs(c.c[1]) > 0) {
    std::snprintf(buf, sizeof buf, "%.10f %.10f %.10f xy xz yz\n", c.b[0], c.c[0], c.c[1]);
    out << buf;
  }
  out << "\nMasses\n\n";
  for (int t = 1; t <= ntypes; ++t) {
    double m = 0;
    std::string label;
    for (const auto& ti : s.types) if (ti.type == t) { m = ti.mass; label = ti.label; }
    if (m == 0)
      for (const auto& a : s.atoms) if (a.type == t) { m = element(a.element).mass; break; }
    std::snprintf(buf, sizeof buf, "%d %.4f", t, m);
    out << buf << (label.empty() ? "" : "  # " + label) << "\n";
  }
  out << "\nAtoms  # full\n\n";
  const std::vector<Vec3> whole = whole_positions(s);   // consistent image flags along every bond
  for (size_t ai = 0; ai < s.atoms.size(); ++ai) {
    const auto& a = s.atoms[ai];
    Vec3 f = c.valid() ? c.to_fractional(whole[ai]) : Vec3{0, 0, 0};
    int im[3] = {0, 0, 0};
    for (int k = 0; k < 3; ++k) im[k] = c.valid() && c.periodic[k] ? static_cast<int>(std::floor(f[k])) : 0;
    const Vec3 w = whole[ai] - (c.a * im[0] + c.b * im[1] + c.c * im[2]);
    std::snprintf(buf, sizeof buf, "%lld %lld %d %.8f %.10f %.10f %.10f %d %d %d", static_cast<long long>(a.id), static_cast<long long>(a.mol), a.type,
                  a.charge, w[0], w[1], w[2], im[0], im[1], im[2]);
    out << buf << residue_comment(a) << "\n";
  }
  if (!s.bonds.empty()) {
    out << "\nBonds\n\n";
    for (size_t k = 0; k < s.bonds.size(); ++k)
      out << k + 1 << " 1 " << s.atoms[s.bonds[k].i].id << " " << s.atoms[s.bonds[k].j].id << "\n";
  }
}

void write_lammps_dump(const Trajectory& t, const std::string& path) {
  std::ofstream out(path);
  if (!out) throw ReadError("cannot write " + path);
  const System& top = t.topology;
  char buf[256];
  for (size_t k = 0; k < t.frames(); ++k) {
    const Cell& c = k < t.cells.size() ? t.cells[k] : top.cell;
    const bool tri = std::fabs(c.b[0]) + std::fabs(c.c[0]) + std::fabs(c.c[1]) > 0;
    out << "ITEM: TIMESTEP\n" << (k < t.timesteps.size() ? t.timesteps[k] : int64_t(k)) << "\nITEM: NUMBER OF ATOMS\n" << top.atoms.size() << "\n";
    // LAMMPS bounding-box form of the cell
    const double xlo = c.origin[0], ylo = c.origin[1], zlo = c.origin[2];
    const double xhi = xlo + c.a[0], yhi = ylo + c.b[1], zhi = zlo + c.c[2];
    if (tri) {
      const double xy = c.b[0], xz = c.c[0], yz = c.c[1];
      std::snprintf(buf, sizeof buf, "ITEM: BOX BOUNDS xy xz yz pp pp pp\n%.6f %.6f %.6f\n%.6f %.6f %.6f\n%.6f %.6f %.6f\n",
                    xlo + std::min({0.0, xy, xz, xy + xz}), xhi + std::max({0.0, xy, xz, xy + xz}), xy, ylo + std::min(0.0, yz),
                    yhi + std::max(0.0, yz), xz, zlo, zhi, yz);
    } else {
      std::snprintf(buf, sizeof buf, "ITEM: BOX BOUNDS pp pp pp\n%.6f %.6f\n%.6f %.6f\n%.6f %.6f\n", xlo, xhi, ylo, yhi, zlo, zhi);
    }
    out << buf << "ITEM: ATOMS id mol type xu yu zu\n";
    const auto& pos = t.positions[k];
    for (size_t i = 0; i < top.atoms.size(); ++i) {
      const auto& a = top.atoms[i];
      std::snprintf(buf, sizeof buf, "%lld %lld %d %.5f %.5f %.5f\n", static_cast<long long>(a.id), static_cast<long long>(a.mol), a.type,
                    pos[i][0], pos[i][1], pos[i][2]);
      out << buf;
    }
  }
}

}  // namespace caps
