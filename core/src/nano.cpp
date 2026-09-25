// CAPS nanostructures (see caps/nano.hpp).
#include "caps/nano.hpp"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <numeric>
#include <stdexcept>

#include "caps/crystal.hpp"
#include "caps/elements.hpp"

namespace caps {

namespace {

constexpr double kPi = 3.14159265358979323846;

void add(System& s, int z, const Vec3& p, int64_t mol = 1) {
  int type = 0;
  for (const auto& t : s.types)
    if (t.label == element(z).symbol) type = t.type;
  if (!type) {
    TypeInfo ti;
    ti.type = int(s.types.size()) + 1;
    ti.mass = element(z).mass;
    ti.label = element(z).symbol;
    s.types.push_back(ti);
    type = ti.type;
  }
  Atom a;
  a.id = int64_t(s.atoms.size() + 1);
  a.mol = mol;
  a.element = z;
  a.type = type;
  a.name = element(z).symbol;
  a.pos = p;
  s.atoms.push_back(a);
}

// Hydrogen on every carbon with two carbon neighbours (graphene edges, open tube ends), in the plane of its bonds.
int cap_edges(System& s) {
  // carbons hanging on by one bond are trimmed (again, until none is left), the rest of the edge is capped
  for (bool again = true; again;) {
    again = false;
    std::vector<int> deg(s.atoms.size(), 0);
    for (const auto& b : crystal_bonds(s)) ++deg[b.i], ++deg[b.j];
    System t;
    t.cell = s.cell;
    for (size_t i = 0; i < s.atoms.size(); ++i) {
      if (s.atoms[i].element == 6 && deg[i] <= 1) { again = true; continue; }
      add(t, s.atoms[i].element, s.atoms[i].pos, s.atoms[i].mol);
    }
    if (again) s = std::move(t);
  }
  const auto bonds = crystal_bonds(s);
  std::vector<std::vector<Vec3>> dirs(s.atoms.size());
  for (const auto& b : bonds) {
    const Vec3 d = s.cell.minimum_image(s.atoms[b.j].pos - s.atoms[b.i].pos);
    dirs[b.i].push_back(d * (1 / norm(d)));
    dirs[b.j].push_back(d * (-1 / norm(d)));
  }
  const size_t n0 = s.atoms.size();
  int added = 0;
  for (size_t i = 0; i < n0; ++i) {
    if (s.atoms[i].element != 6 || dirs[i].size() != 2) continue;
    Vec3 u = (dirs[i][0] + dirs[i][1]) * -1.0;
    const double l = norm(u);
    if (l < 1e-6) continue;
    add(s, 1, s.atoms[i].pos + u * (1.09 / l), s.atoms[i].mol);
    ++added;
  }
  s.bonds = crystal_bonds(s);
  return added;
}

void finish(System& s, const std::string& title) {
  s.title = title;
  s.has_mol = true;
  s.unwrapped = true;
  s.source_format = "caps-nano";
  if (s.bonds.empty()) s.bonds = crystal_bonds(s);
}

bool inside(ParticleShape sh, const Vec3& r, double R) {
  const double ax = std::fabs(r[0]), ay = std::fabs(r[1]), az = std::fabs(r[2]);
  switch (sh) {
    case ParticleShape::Sphere: return norm(r) <= R;
    case ParticleShape::Cube: return std::max({ax, ay, az}) <= R / std::sqrt(3.0);
    case ParticleShape::Octahedron: return ax + ay + az <= R;
    case ParticleShape::Cuboctahedron: {
      const double a = R / std::sqrt(2.0);   // vertices at (±a, ±a, 0) and permutations
      return std::max({ax, ay, az}) <= a && ax + ay + az <= 2 * a;
    }
  }
  return false;
}

}  // namespace

const char* to_string(ParticleShape s) {
  switch (s) {
    case ParticleShape::Sphere: return "sphere";
    case ParticleShape::Cube: return "cube";
    case ParticleShape::Octahedron: return "octahedron";
    case ParticleShape::Cuboctahedron: return "cuboctahedron";
  }
  return "sphere";
}

ParticleShape particle_shape_from_string(const std::string& s) {
  if (s == "cube") return ParticleShape::Cube;
  if (s == "octahedron") return ParticleShape::Octahedron;
  if (s == "cuboctahedron") return ParticleShape::Cuboctahedron;
  if (s == "sphere") return ParticleShape::Sphere;
  throw std::invalid_argument("shape must be sphere, cube, octahedron or cuboctahedron");
}

// ---------------------------------------------------------------- graphene

System graphene_sheet(const SheetOptions& o, NanoReport* rep) {
  const double cc = o.cc > 0 ? o.cc : 1.42, a = std::sqrt(3.0) * cc, h = 3 * cc;   // rectangular cell a × 3cc, 4 atoms
  const int nx = std::max(3, int(std::lround(o.lx / a))), ny = std::max(2, int(std::lround(o.ly / h)));
  const int layers = std::max(1, o.layers);
  const double gap = 3.35;
  const double Lx = nx * a, Ly = ny * h;
  System s;
  const double vac = std::max(0.0, o.vacuum);
  const double off_xy = o.periodic ? 0.0 : vac / 2;
  s.cell.a = {Lx + (o.periodic ? 0 : vac), 0, 0};
  s.cell.b = {0, Ly + (o.periodic ? 0 : vac), 0};
  s.cell.c = {0, 0, (layers - 1) * gap + std::max(vac, 2.0 * gap)};
  const double z0 = std::max(vac, 2.0 * gap) / 2;
  const Vec3 basis[4] = {{0, 0, 0}, {a / 2, cc / 2, 0}, {a / 2, 1.5 * cc, 0}, {0, 2 * cc, 0}};
  for (int l = 0; l < layers; ++l) {
    const double sy = (l % 2) ? cc : 0.0;   // AB stacking: every other layer shifted by one C–C along y
    for (int i = 0; i < nx; ++i)
      for (int j = 0; j < ny; ++j)
        for (const auto& b : basis) {
          double y = j * h + b[1] + sy;
          if (o.periodic) y -= Ly * std::floor(y / Ly);
          add(s, 6, {off_xy + i * a + b[0], off_xy + y, z0 + l * gap});
        }
  }
  finish(s, "graphene " + std::to_string(layers) + (layers > 1 ? " layers" : " layer"));
  NanoReport r;
  if (!o.periodic) r.capped = cap_edges(s);
  char b[160];
  std::snprintf(b, sizeof b, "graphene %.2f × %.2f Å (%d × %d cells), %d layer%s, %zu atoms%s", Lx, Ly, nx, ny, layers, layers > 1 ? "s" : "", s.atoms.size(),
                o.periodic ? " · periodic in the plane" : " · flake, edges capped with H");
  r.notes.push_back(b);
  s.notes = r.notes;
  if (rep) *rep = r;
  return s;
}

// ---------------------------------------------------------------- nanotubes

std::array<double, 3> nanotube_geometry(int n, int m, double cc) {
  const double a = std::sqrt(3.0) * cc;
  const double ch = a * std::sqrt(double(n) * n + double(n) * m + double(m) * m);
  const int dr = std::gcd(2 * m + n, 2 * n + m);
  const double T = std::sqrt(3.0) * ch / dr;
  const double theta = std::atan2(std::sqrt(3.0) * m, 2.0 * n + m) * 180 / kPi;
  return {ch / kPi, theta, T};
}

System nanotube(const NanotubeOptions& o, NanoReport* rep) {
  const int n = o.n, m = o.m;
  if (n < 1 || m < 0 || m > n) throw std::invalid_argument("chirality (n, m) needs n ≥ 1 and 0 ≤ m ≤ n");
  const double cc = o.cc > 0 ? o.cc : 1.42, a = std::sqrt(3.0) * cc;
  const Vec3 a1{a, 0, 0}, a2{a / 2, a * std::sqrt(3.0) / 2, 0};
  const Vec3 Ch = a1 * n + a2 * m;
  const int dr = std::gcd(2 * m + n, 2 * n + m);
  const int t1 = (2 * m + n) / dr, t2 = -(2 * n + m) / dr;
  const Vec3 T = a1 * t1 + a2 * t2;
  const double lc = norm(Ch), lt = norm(T), R = lc / (2 * kPi);
  const Vec3 basis[2] = {{0, 0, 0}, (a1 + a2) * (1.0 / 3)};
  // atoms of one period: graphene points inside the parallelogram of Ch and T
  std::vector<std::pair<double, double>> uv;
  const int span = 2 * (std::abs(n) + std::abs(m) + std::abs(t1) + std::abs(t2)) + 2;
  for (int i = -span; i <= span; ++i)
    for (int j = -span; j <= span; ++j)
      for (const auto& b : basis) {
        const Vec3 r = a1 * i + a2 * j + b;
        const double u = dot(r, Ch) / (lc * lc), v = dot(r, T) / (lt * lt);
        if (u >= -1e-9 && u < 1 - 1e-9 && v >= -1e-9 && v < 1 - 1e-9) uv.push_back({u, v});
      }
  const int expect = 4 * (n * n + n * m + m * m) / dr;
  if (int(uv.size()) != expect) throw std::runtime_error("nanotube: " + std::to_string(uv.size()) + " atoms per period, expected " + std::to_string(expect));
  int periods = std::max(1, int(std::lround(o.length / lt)));
  if (o.periodic) periods = std::max(periods, int(std::ceil(6.0 / lt)));   // at least 6 Å along the axis
  const double Lz = periods * lt, vac = std::max(0.0, o.vacuum);
  System s;
  const double box = 2 * R + 2 * vac;
  s.cell.a = {box, 0, 0};
  s.cell.b = {0, box, 0};
  s.cell.c = {0, 0, o.periodic ? Lz : Lz + 2 * vac};
  const double zoff = o.periodic ? 0 : vac;
  for (int p = 0; p < periods; ++p)
    for (const auto& [u, v] : uv) {
      const double th = 2 * kPi * u;
      add(s, 6, {box / 2 + R * std::cos(th), box / 2 + R * std::sin(th), zoff + (v + p) * lt});
    }
  NanoReport r;
  const auto g = nanotube_geometry(n, m, cc);
  r.diameter = g[0], r.chiral_angle = g[1], r.translation = g[2], r.atoms_per_period = expect;
  finish(s, "(" + std::to_string(n) + "," + std::to_string(m) + ") nanotube");
  if (!o.periodic) r.capped = cap_edges(s);
  const char* kind = m == n ? "armchair" : m == 0 ? "zigzag" : "chiral";
  char b[200];
  std::snprintf(b, sizeof b, "(%d,%d) %s nanotube · d = %.2f Å · chiral angle %.2f° · |T| = %.3f Å · %d atoms per period × %d · %s", n, m, kind, r.diameter,
                r.chiral_angle, r.translation, expect, periods, o.periodic ? "periodic along z" : "finite, ends capped with H");
  r.notes.push_back(b);
  s.notes = r.notes;
  if (rep) *rep = r;
  return s;
}

// ---------------------------------------------------------------- particles

System nanoparticle(const System& bulk, const ParticleOptions& o, NanoReport* rep) {
  if (!bulk.cell.valid() || bulk.atoms.empty()) throw std::invalid_argument("the crystal has no cell or no atoms");
  if (o.radius <= 0) throw std::invalid_argument("the particle radius must be positive");
  const Cell& c = bulk.cell;
  const double V = std::fabs(dot(c.a, cross(c.b, c.c)));
  const double w[3] = {V / norm(cross(c.b, c.c)), V / norm(cross(c.c, c.a)), V / norm(cross(c.a, c.b))};
  int nrep[3];
  for (int k = 0; k < 3; ++k) nrep[k] = int(std::ceil((o.radius + 3) / w[k])) + 1;
  // the centre: the cell centre, or the atom nearest it
  Vec3 centre = c.origin + (c.a + c.b + c.c) * 0.5;
  if (o.on_atom) {
    double best = 1e300;
    Vec3 pick = centre;
    for (const auto& at : bulk.atoms) {
      const double d = norm(c.minimum_image(at.pos - centre));
      if (d < best) best = d, pick = centre + c.minimum_image(at.pos - centre);
    }
    centre = pick;
  }
  const double vac = std::max(0.0, o.vacuum), box = 2 * o.radius + 2 * vac + 4;
  System s;
  s.cell.a = {box, 0, 0};
  s.cell.b = {0, box, 0};
  s.cell.c = {0, 0, box};
  const Vec3 mid{box / 2, box / 2, box / 2};
  for (int i = -nrep[0]; i <= nrep[0]; ++i)
    for (int j = -nrep[1]; j <= nrep[1]; ++j)
      for (int k = -nrep[2]; k <= nrep[2]; ++k)
        for (const auto& at : bulk.atoms) {
          const Vec3 r = at.pos + c.a * i + c.b * j + c.c * k - centre;
          if (inside(o.shape, r, o.radius)) add(s, at.element, mid + r);
        }
  if (s.atoms.empty()) throw std::invalid_argument("no atom inside the particle: make it larger");
  s.bonds = crystal_bonds(s);
  // atoms left without a bonded partner (a lone O at the surface of an oxide) are dropped, unless the crystal is a metal
  {
    std::vector<int> deg(s.atoms.size(), 0);
    for (const auto& b : s.bonds) ++deg[b.i], ++deg[b.j];
    const bool any_bond = std::any_of(deg.begin(), deg.end(), [](int d) { return d > 0; });
    if (any_bond) {
      System t;
      t.cell = s.cell;
      for (size_t i = 0; i < s.atoms.size(); ++i)
        if (deg[i] > 0) add(t, s.atoms[i].element, s.atoms[i].pos);
      s = std::move(t);
      s.bonds = crystal_bonds(s);
    }
  }
  NanoReport r;
  if (o.passivate) {
    const auto [nh, noh] = passivate_surface(s, bulk, [&](const Vec3& p) {
      const Vec3 d = p - mid;
      const double l = norm(d);
      return l > 1e-9 ? d * (1 / l) : Vec3{0, 0, 1};
    });
    r.added_h = nh, r.added_oh = noh;
  }
  double rmax = 0;
  for (const auto& at : s.atoms) rmax = std::max(rmax, norm(at.pos - mid));
  r.diameter = 2 * rmax;
  const std::string name = bulk.title.empty() ? std::string("crystal") : bulk.title;
  finish(s, name + " " + to_string(o.shape));
  std::vector<size_t> all(s.atoms.size());
  std::iota(all.begin(), all.end(), size_t(0));
  char b[220];
  std::snprintf(b, sizeof b, "%s %s · circumscribed radius %.1f Å · %zu atoms (%s) · diameter %.1f Å", name.c_str(), to_string(o.shape), o.radius,
                s.atoms.size(), formula_of(s, all).c_str(), r.diameter);
  r.notes.push_back(b);
  if (r.added_h || r.added_oh) r.notes.push_back("passivated: " + std::to_string(r.added_oh) + " OH, " + std::to_string(r.added_h) + " H");
  s.notes = r.notes;
  if (rep) *rep = r;
  return s;
}

// ---------------------------------------------------------------- filler in a matrix

double occupied_volume(const System& s) {
  if (s.atoms.empty()) return 0;
  Vec3 lo{1e300, 1e300, 1e300}, hi{-1e300, -1e300, -1e300};
  for (const auto& a : s.atoms)
    for (int k = 0; k < 3; ++k) lo[k] = std::min(lo[k], a.pos[k] - 2.5), hi[k] = std::max(hi[k], a.pos[k] + 2.5);
  const double g = 0.4;
  int n[3];
  for (int k = 0; k < 3; ++k) n[k] = std::max(1, int(std::ceil((hi[k] - lo[k]) / g)));
  std::vector<char> in(size_t(n[0]) * n[1] * n[2], 0);
  for (const auto& a : s.atoms) {
    const double r = element(a.element).vdw;
    int b0[3], b1[3];
    for (int k = 0; k < 3; ++k) b0[k] = std::max(0, int((a.pos[k] - r - lo[k]) / g)), b1[k] = std::min(n[k] - 1, int((a.pos[k] + r - lo[k]) / g));
    for (int i = b0[0]; i <= b1[0]; ++i)
      for (int j = b0[1]; j <= b1[1]; ++j)
        for (int k = b0[2]; k <= b1[2]; ++k) {
          const Vec3 p{lo[0] + (i + 0.5) * g, lo[1] + (j + 0.5) * g, lo[2] + (k + 0.5) * g};
          if (norm(p - a.pos) <= r) in[(size_t(i) * n[1] + j) * n[2] + k] = 1;
        }
  }
  return double(std::count(in.begin(), in.end(), 1)) * g * g * g;
}

System embed_filler(const System& filler, const ChainSpec& spec, const FillerMatrixOptions& o, FillerReport* rep) {
  if (filler.atoms.empty()) throw GrowError("the filler has no atoms");
  if (o.chains < 1) throw GrowError("give at least one chain");
  if (o.density <= 0) throw GrowError("give the matrix density");
  const Cell& fc = filler.cell;
  const std::array<double, 3> fl{norm(fc.a), norm(fc.b), norm(fc.c)};
  int kept = 0;
  for (int k = 0; k < 3; ++k) kept += o.keep_axis[size_t(k)];
  if (kept == 3) throw GrowError("at least one cell edge must be free to hold the matrix");
  if (kept > 0 && (std::fabs(fc.a[1]) + std::fabs(fc.a[2]) + std::fabs(fc.b[0]) + std::fabs(fc.b[2]) + std::fabs(fc.c[0]) + std::fabs(fc.c[1]) > 1e-6))
    throw GrowError("a filler with periodic axes needs a rectangular cell");
  // the matrix volume from its mass and density, the filler's from its van der Waals spheres
  double mchain = 0;
  for (int k = 0; k < 8; ++k) mchain += chain_mass(spec, chain_sequence(spec, o.grow.seed + uint64_t(k) * 101));
  mchain /= 8;
  const double vm = o.chains * mchain / (o.density * 0.602214076);
  const double vf = occupied_volume(filler);
  double fixed = 1;
  for (int k = 0; k < 3; ++k)
    if (o.keep_axis[size_t(k)]) fixed *= fl[size_t(k)];
  const double free_len = std::pow((vm + vf) / fixed, 1.0 / (3 - kept));
  std::array<double, 3> L;
  // extent of the filler along the free axes: the cell must leave room for chains around it
  Vec3 lo{1e300, 1e300, 1e300}, hi{-1e300, -1e300, -1e300};
  for (const auto& a : filler.atoms)
    for (int k = 0; k < 3; ++k) lo[k] = std::min(lo[k], a.pos[k]), hi[k] = std::max(hi[k], a.pos[k]);
  FillerReport R;
  std::string widened;
  for (int k = 0; k < 3; ++k) {
    if (o.keep_axis[size_t(k)]) { L[size_t(k)] = fl[size_t(k)]; continue; }
    L[size_t(k)] = free_len;
    const double need = hi[k] - lo[k] + 7.0;
    if (L[size_t(k)] < need) {
      L[size_t(k)] = need;
      widened += std::string(widened.empty() ? "" : ", ") + "xyz"[k];
    }
  }
  if (!widened.empty()) {
    const double room = L[0] * L[1] * L[2] - vf;
    const int want = int(std::ceil(room * o.density * 0.602214076 / mchain));
    R.notes.push_back("the cell was widened along " + widened + " to fit the filler, so the matrix is below " + std::to_string(o.density).substr(0, 4) +
                      " g/cm³: about " + std::to_string(want) + " chains would fill it");
  }
  // the filler centred (free axes) in the new cell
  System sub = filler;
  for (auto& a : sub.atoms)
    for (int k = 0; k < 3; ++k) {
      if (o.keep_axis[size_t(k)]) a.pos[k] -= fc.origin[k];
      else a.pos[k] += L[size_t(k)] / 2 - (lo[k] + hi[k]) / 2;
      if (o.keep_axis[size_t(k)]) a.pos[k] -= L[size_t(k)] * std::floor(a.pos[k] / L[size_t(k)]);
    }
  for (auto& a : sub.atoms) a.mol = 1;
  GrowOptions g = o.grow;
  g.cell = L;
  g.substrate = &sub;
  g.chains = o.chains;
  g.z_lo = g.z_hi = 0;
  g.auto_scale = true;
  GrowReport gr;
  System s = grow_chains(spec, g, &gr);
  double mf = 0, mt = 0;
  for (size_t i = 0; i < s.atoms.size(); ++i) {
    const double m = s.mass_of(s.atoms[i]);
    mt += m;
    if (i < sub.atoms.size()) mf += m;
  }
  const double vol = L[0] * L[1] * L[2];
  R.filler_mass_fraction = mt > 0 ? mf / mt : 0;
  R.filler_volume_fraction = vf / vol;
  R.density = mt / (0.602214076 * vol);
  R.cell = L;
  char b[240];
  std::snprintf(b, sizeof b, "filler %.1f %% by mass, %.1f %% by volume · cell %.2f × %.2f × %.2f Å · %.3f g/cm³ overall", 100 * R.filler_mass_fraction,
                100 * R.filler_volume_fraction, L[0], L[1], L[2], R.density);
  R.notes.insert(R.notes.begin(), b);
  for (const auto& n : gr.notes) R.notes.push_back(n);
  s.title = "CAPS composite: " + (filler.title.empty() ? std::string("filler") : filler.title) + " + " + std::to_string(o.chains) + " chains";
  s.notes = R.notes;
  if (rep) *rep = R;
  return s;
}

}  // namespace caps
