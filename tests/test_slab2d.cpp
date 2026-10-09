// 2D sheets and terminations (caps/slab2d.hpp): the reference Ti3C2Tx slabs rebuilt exactly, the height from the bond,
// the validator's PASS and FAIL cases, and NumPy's generator bit for bit.
#include <gtest/gtest.h>

#include <cmath>
#include <filesystem>
#include <fstream>

#include "../core/src/npy_random.hpp"
#include "caps/slab2d.hpp"
#include "caps/elements.hpp"

using namespace caps;

namespace {
const std::string kData = std::string(CAPS_SOURCE_DIR) + "/data";
const std::string kRef = std::string(CAPS_SOURCE_DIR) + "/tests/data/slab2d/";

System build(const std::string& spec, int n, uint64_t seed = 1) {
  const System sheet = sheet_from_layers(sheet_preset(kData, "Ti3C2"));
  TerminateOptions o;
  o.top.fractions = parse_fractions(spec);
  o.na = o.nb = n;
  o.seed = seed;
  o.data_dir = kData;
  return terminate_slab(sheet, o);
}

void expect_same(const System& a, const System& b, const std::string& what) {
  ASSERT_EQ(a.atoms.size(), b.atoms.size()) << what;
  for (int k = 0; k < 3; ++k) {
    const Vec3 ca = k == 0 ? a.cell.a : k == 1 ? a.cell.b : a.cell.c, cb = k == 0 ? b.cell.a : k == 1 ? b.cell.b : b.cell.c;
    for (int q = 0; q < 3; ++q) EXPECT_NEAR(ca[q], cb[q], 1e-6) << what << " cell";
  }
  double worst = 0;
  for (size_t i = 0; i < a.atoms.size(); ++i) {
    ASSERT_EQ(a.atoms[i].element, b.atoms[i].element) << what << " atom " << i;
    for (int q = 0; q < 3; ++q) worst = std::max(worst, std::fabs(a.atoms[i].pos[q] - b.atoms[i].pos[q]));
  }
  EXPECT_LT(worst, 1e-6) << what;
}
}  // namespace

TEST(Slab2D, NumpyGeneratorBitForBit) {
  // np.random.default_rng(12345).random() == 0.22733602246716966 (NumPy's documentation)
  npy::Pcg64 g(12345);
  EXPECT_DOUBLE_EQ(double(g.next64() >> 11) * (1.0 / 9007199254740992.0), 0.22733602246716966);
}

TEST(Slab2D, ReferenceSlabsRebuiltExactly) {
  expect_same(build("O", 1), read_vasp_poscar(kRef + "Ti3C2_O.POSCAR"), "O");
  expect_same(build("OH", 1), read_vasp_poscar(kRef + "Ti3C2_OH.POSCAR"), "OH");
  expect_same(build("F", 1), read_vasp_poscar(kRef + "Ti3C2_F.POSCAR"), "F");
  expect_same(build("O:0.5,OH:0.25,F:0.25", 3), read_vasp_poscar(kRef + "Ti3C2_mixed.POSCAR"), "mixed 3x3 seed 1");
  // written as the reference was (comment "<formula> CAPS <label>")
  const auto dir = std::filesystem::temp_directory_path();
  const System m = build("O:0.5,OH:0.25,F:0.25", 3);
  write_vasp_poscar(m, (dir / "caps_mixed.POSCAR").string(), formula_ordered(m) + " CAPS Ti3C2_mixed");
  auto slurp = [](const std::string& p) { std::ifstream f(p); return std::string((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>()); };
  // the header (comment, cell, species and counts, Cartesian) the same text; coordinates agree within 1e-6 above
  auto head = [](const std::string& t) { size_t p = 0; for (int k = 0; k < 8; ++k) p = t.find('\n', p) + 1; return t.substr(0, p); };
  EXPECT_EQ(head(slurp((dir / "caps_mixed.POSCAR").string())), head(slurp(kRef + "Ti3C2_mixed.POSCAR")));
}

TEST(Slab2D, HeightFromBondLength) {
  // Ti–O 2.05 Å over a hollow of a = 3.0814 Å: h = √(d² − (a/√3)²) = 1.0186 Å; O–H 0.97 Å straight out on both faces
  const System o = build("O", 1);
  double ti_top = -1e9, o_top = -1e9;
  for (const auto& a : o.atoms) { if (a.element == 22) ti_top = std::max(ti_top, a.pos[2]); if (a.element == 8) o_top = std::max(o_top, a.pos[2]); }
  EXPECT_NEAR(o_top - ti_top, 1.0186, 1e-4);
  const System oh = build("OH", 1);
  for (const auto& h : oh.atoms) {
    if (h.element != 1) continue;
    double best = 1e9;
    Vec3 d{};
    for (const auto& x : oh.atoms) if (x.element == 8 && norm(x.pos - h.pos) < best) best = norm(x.pos - h.pos), d = h.pos - x.pos;
    EXPECT_NEAR(best, 0.97, 1e-9);
    EXPECT_NEAR(std::fabs(d[2]), 0.97, 1e-9);
  }
}

TEST(Slab2D, ValidatorPassAndFail) {
  ValidateOptions v;
  v.data_dir = kData;
  for (const auto& [f, expect] : std::vector<std::pair<std::string, std::string>>{{"O", "Ti3C2O2"}, {"OH", "Ti3C2O2H2"}, {"F", "Ti3C2F2"}, {"mixed", ""}}) {
    std::string comment;
    const System s = read_vasp_poscar(kRef + "Ti3C2_" + f + ".POSCAR", nullptr, &comment);
    v.expect = expect;
    v.core = "Ti3C2";
    v.comment = comment;
    const auto r = validate_2d(s, v);
    EXPECT_EQ(r.status, "PASS") << f << "\n" << r.text();
  }
  v.expect = "";
  const auto dir = std::filesystem::temp_directory_path();
  // the digit 0 as a species label
  {
    std::ifstream in(kRef + "Ti3C2_O.POSCAR");
    std::string t((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    t.replace(t.find("  Ti  C  O\n"), 11, "  Ti  C  0\n");
    std::ofstream((dir / "caps_zero.POSCAR").string()) << t;
    std::map<std::string, std::string> bad;
    const System s = read_vasp_poscar((dir / "caps_zero.POSCAR").string(), &bad);
    v.bad_labels = bad;
    EXPECT_EQ(validate_2d(s, v).status, "FAIL");
    v.bad_labels.clear();
  }
  // two stacked sheets
  {
    System s = read_vasp_poscar(kRef + "Ti3C2_O.POSCAR");
    const size_t n = s.atoms.size();
    for (size_t i = 0; i < n; ++i) { Atom a = s.atoms[i]; a.pos[2] += 9.0; s.atoms.push_back(a); }
    const auto r = validate_2d(s, v);
    EXPECT_EQ(r.status, "FAIL");
    EXPECT_NE(r.text().find("complete sheets"), std::string::npos) << r.text();
  }
  // an H not bonded to O
  {
    System s = read_vasp_poscar(kRef + "Ti3C2_OH.POSCAR");
    for (auto& a : s.atoms) if (a.element == 1) { a.pos[2] += 0.8; break; }
    const auto r = validate_2d(s, v);
    EXPECT_EQ(r.status, "FAIL");
    EXPECT_NE(r.text().find("not bonded to any O"), std::string::npos) << r.text();
  }
  // 8 Å of vacuum
  {
    System s = read_vasp_poscar(kRef + "Ti3C2_O.POSCAR");
    double z0 = 1e9, z1 = -1e9;
    for (const auto& a : s.atoms) z0 = std::min(z0, a.pos[2]), z1 = std::max(z1, a.pos[2]);
    for (auto& a : s.atoms) a.pos[2] += 4.0 - z0;
    s.cell.c = {0, 0, z1 - z0 + 8.0};
    const auto r = validate_2d(s, v);
    EXPECT_EQ(r.status, "FAIL");
    EXPECT_NE(r.text().find("vacuum gap only"), std::string::npos) << r.text();
  }
}

TEST(Slab2D, SitesAndIsolate) {
  const System sheet = sheet_from_layers(sheet_preset(kData, "Ti3C2"));
  const auto sites = surface_sites(sheet);
  int fcc = 0, hcp = 0, top = 0;
  for (const auto& s : sites) if (s.face == "top") fcc += s.kind == "fcc", hcp += s.kind == "hcp", top += s.kind == "top";
  EXPECT_EQ(fcc, 1);
  EXPECT_EQ(hcp, 1);
  EXPECT_EQ(top, 1);
  // a stack of two sheets 9 Å apart in a 2 × 9 Å cell: either sheet comes out complete
  System two = sheet;
  for (const auto& a : sheet.atoms) { Atom b = a; b.pos[2] += 9.0; two.atoms.push_back(b); }
  two.cell.c = {0, 0, 18.0};
  IsolateOptions io;
  io.formula = "Ti3C2";
  IsolateReport rep;
  const System one = isolate_layer(two, io, &rep);
  EXPECT_EQ(rep.complete, 2);
  EXPECT_EQ(one.atoms.size(), 5u);
  // a MAX-type stack (an Al layer between the sheets) gives the sheet once Al is removed
  System max = two;
  Atom al;
  al.element = 13;
  al.pos = {0, 0, 7.0};
  max.atoms.push_back(al);
  io.remove = {"Al"};
  EXPECT_EQ(isolate_layer(max, io).atoms.size(), 5u);
}
