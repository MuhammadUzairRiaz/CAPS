#include <gtest/gtest.h>

#include <cmath>
#include <set>

#include "caps/functionalize.hpp"
#include "caps/nano.hpp"

using namespace caps;

namespace {
double radius(const System& s, const Vec3& p, const Vec3& c) { return std::hypot(p[0] - c[0], p[1] - c[1]); }
Vec3 axis_centre(const System& s) { return s.cell.a * 0.5 + s.cell.b * 0.5; }
// the atom a group attaches through: the new heavy atom bonded to the site
uint32_t attached(const System& s, uint32_t site, size_t first_new) {
  for (const auto& b : s.bonds) {
    if (b.i == site && b.j >= first_new && s.atoms[b.j].element != 1) return b.j;
    if (b.j == site && b.i >= first_new && s.atoms[b.i].element != 1) return b.i;
  }
  return uint32_t(-1);
}
}  // namespace

// Carboxyl groups on a (10,10) tube: a twentieth of the sidewall carbons, spaced; each outside the wall (or inside with
// side "inner"), bonded to its carbon
TEST(Functionalize, TubeSidewallOutsideAndInside) {
  NanotubeOptions t;
  t.length = 20;
  NanoReport tr;
  System s = nanotube(t, &tr);
  const size_t n0 = s.atoms.size();
  const double R = tr.diameter / 2;
  const Vec3 c = axis_centre(s);
  FunctionalizeOptions o;
  o.group = "carboxyl";
  o.fraction = 0.05;
  const auto r = functionalize(s, o);
  EXPECT_EQ(r.eligible, n0);
  EXPECT_EQ(r.grafted, size_t(std::lround(0.05 * double(n0))));
  EXPECT_EQ(r.added_atoms, 4 * r.grafted);   // C, O, O, H each
  for (uint32_t site : r.sites) {
    const uint32_t x = attached(s, site, n0);
    ASSERT_NE(x, uint32_t(-1));
    EXPECT_GT(radius(s, s.atoms[x].pos, c), R + 1.0);
  }
  System s2 = nanotube(t);
  o.side = "inner";
  o.fraction = 0.02;
  const auto r2 = functionalize(s2, o);
  for (uint32_t site : r2.sites) EXPECT_LT(radius(s2, s2.atoms[attached(s2, site, n0)].pos, c), R - 1.0);
}

// A flat sheet: hydroxyls on top (+z), every site at least the spacing apart; "both" puts some below; boron nitride with
// elements B takes only borons
TEST(Functionalize, SheetSidesSpacingAndElements) {
  SheetOptions so;
  System s = graphene_sheet(so);
  const double z0 = s.atoms[0].pos[2];
  const size_t n0 = s.atoms.size();
  FunctionalizeOptions o;
  o.group = "hydroxyl";
  o.pattern = "all";
  o.min_spacing = 5.0;
  const auto r = functionalize(s, o);
  EXPECT_GT(r.grafted, 5u);
  for (size_t a = 0; a < r.sites.size(); ++a) {
    EXPECT_GT(s.atoms[attached(s, r.sites[a], n0)].pos[2], z0 + 1.0);
    for (size_t b = a + 1; b < r.sites.size(); ++b) EXPECT_GE(norm(s.cell.minimum_image(s.atoms[r.sites[a]].pos - s.atoms[r.sites[b]].pos)), 5.0 - 1e-9);
  }
  System s2 = graphene_sheet(so);
  o.side = "both";
  const auto r2 = functionalize(s2, o);
  int below = 0;
  for (uint32_t site : r2.sites) below += s2.atoms[attached(s2, site, n0)].pos[2] < z0;
  EXPECT_GT(below, 0);
  EXPECT_LT(below, int(r2.grafted));
  so.material = "h-BN";
  System bn = graphene_sheet(so);
  o.elements = "B";
  o.side = "outer";
  const auto r3 = functionalize(bn, o);
  for (uint32_t site : r3.sites) EXPECT_EQ(bn.atoms[site].element, 5);
}

// A finite tube's rim: the hydrogens replaced by amines; a helix follows 2πz/pitch
TEST(Functionalize, EndsAndHelix) {
  NanotubeOptions t;
  t.length = 20;
  t.periodic = false;
  System s = nanotube(t);
  int h0 = 0;
  for (const auto& a : s.atoms) h0 += a.element == 1;
  FunctionalizeOptions o;
  o.group = "amine";
  o.pattern = "ends";
  o.min_spacing = 0;
  const auto r = functionalize(s, o);
  EXPECT_EQ(int(r.grafted), h0);   // every rim hydrogen became an NH2
  int nN = 0, h1 = 0;
  for (const auto& a : s.atoms) nN += a.element == 7, h1 += a.element == 1;
  EXPECT_EQ(nN, h0);
  EXPECT_EQ(h1, 2 * h0);
  NanotubeOptions tp;
  tp.length = 40;
  System hx = nanotube(tp);
  FunctionalizeOptions oh;
  oh.group = "fluoro";
  oh.pattern = "helix";
  oh.pitch = 25;
  oh.min_spacing = 2.6;
  const auto rh = functionalize(hx, oh);
  EXPECT_GT(rh.grafted, 5u);
  const Vec3 c = axis_centre(hx);
  double zlo = 1e9;
  for (const auto& a : hx.atoms) if (a.element == 6) zlo = std::min(zlo, a.pos[2]);
  for (uint32_t site : rh.sites) {
    const auto& p = hx.atoms[site].pos;
    const double d = std::remainder(std::atan2(p[1] - c[1], p[0] - c[0]) - 2 * 3.14159265358979323846 * (p[2] - zlo) / 25, 2 * 3.14159265358979323846);
    EXPECT_LE(std::fabs(d), 20 * 3.14159265358979323846 / 180 + 1e-9);
  }
  EXPECT_THROW(functionalize(hx, FunctionalizeOptions{"CO"}), std::invalid_argument);
}
