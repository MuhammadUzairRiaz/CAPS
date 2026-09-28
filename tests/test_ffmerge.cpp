#include <gtest/gtest.h>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <sstream>

#include "caps/analysis.hpp"
#include "caps/ffdef.hpp"
#include "caps/ffmerge.hpp"
#include "caps/io.hpp"
#include "caps/json.hpp"
#include "caps/elements.hpp"
#include "caps/manybody.hpp"
#include "caps/relax.hpp"
#include "caps/uff.hpp"
#include "caps/typing.hpp"

using namespace caps;

namespace {

const std::string kFF = std::string(CAPS_SOURCE_DIR) + "/data/forcefields/";

// the part of s made of the atoms given (bonds among them), and its force field
System part_of(const System& s, const std::vector<uint32_t>& atoms) {
  System p;
  p.cell = s.cell;
  std::vector<int> at(s.atoms.size(), -1);
  for (size_t k = 0; k < atoms.size(); ++k) { at[atoms[k]] = int(k); p.atoms.push_back(s.atoms[atoms[k]]); }
  for (const auto& b : s.bonds) if (at[b.i] >= 0 && at[b.j] >= 0) p.bonds.push_back({uint32_t(at[b.i]), uint32_t(at[b.j]), b.order});
  p.bonds_from_file = true;
  return p;
}

ForceField typed(const System& s, const FFDef& def) {
  const TypingResult tr = assign_types(s, def);
  EXPECT_EQ(tr.untyped, 0);
  ParamReport rep;
  return parameterize(s, def, tr.types, "gasteiger", &rep, false);
}

EnergyTerms energy(const ForceField& ff, const System& s) {
  EnergyOptions o;
  Evaluator ev(ff, o);
  std::vector<double> x, f;
  for (const auto& a : s.atoms) x.insert(x.end(), a.pos.begin(), a.pos.end());
  return ev.compute(x, s.cell, f);
}

}  // namespace

TEST(FFMerge, SameForceFieldInTwoPartsIsExact) {
  // the polystyrene melt split into chains 1-5 and 6-10, both GAFF2, merged with GAFF's own cross rule (ε geometric,
  // σ arithmetic): every term equals the single assignment's
  System s = open_file(std::string(CAPS_SOURCE_DIR) + "/samples/ps_melt.data").frame(0);
  make_molecules_whole(s);
  const FFDef def = load_forcefield(kFF + "gaff-amber25.json");
  const auto mol = s.molecules();
  std::vector<uint32_t> a, b;
  for (size_t i = 0; i < s.atoms.size(); ++i) (mol[i] < 5 ? a : b).push_back(uint32_t(i));
  // one charge model for both: Gasteiger–Marsili per chain equals per cell (it is per molecule)
  const ForceField whole = typed(s, def);
  const System sa = part_of(s, a), sb = part_of(s, b);
  const ForceField fa = typed(sa, def), fb = typed(sb, def);
  std::vector<std::string> notes;
  const ForceField m = merge_forcefields(s.atoms.size(), {{&fa, a, "A"}, {&fb, b, "B"}}, MergeOptions{}, &notes);
  const EnergyTerms e1 = energy(whole, s), e2 = energy(m, s);
  EXPECT_NEAR(e2.bond, e1.bond, 1e-9);
  EXPECT_NEAR(e2.angle, e1.angle, 1e-9);
  EXPECT_NEAR(e2.dihedral, e1.dihedral, 1e-9);
  EXPECT_NEAR(e2.improper, e1.improper, 1e-9);
  EXPECT_NEAR(e2.vdw, e1.vdw, 1e-8);
  EXPECT_NEAR(e2.coulomb, e1.coulomb, 1e-8);
  EXPECT_FALSE(notes.empty());
  EXPECT_EQ(m.type_names.size(), fa.type_names.size() + fb.type_names.size());   // B's names tagged
}

TEST(FFMerge, DifferentFamiliesAreRefusedOrMergedAsAsked) {
  System s = open_file(std::string(CAPS_SOURCE_DIR) + "/samples/ps_melt.data").frame(0);
  make_molecules_whole(s);
  const auto mol = s.molecules();
  std::vector<uint32_t> a, b;
  for (size_t i = 0; i < s.atoms.size(); ++i) (mol[i] < 5 ? a : b).push_back(uint32_t(i));
  const FFDef gaff = load_forcefield(kFF + "gaff-amber25.json"), opls = load_forcefield(kFF + "opls2005.json");
  const ForceField fa = typed(part_of(s, a), gaff), fb = typed(part_of(s, b), opls);
  // GAFF 1-4 Coulomb 1/1.2, OPLS 1/2: refused unless asked
  EXPECT_THROW(merge_forcefields(s.atoms.size(), {{&fa, a, "GAFF"}, {&fb, b, "OPLS"}}, MergeOptions{}), FieldError);
  MergeOptions o;
  o.scaling14 = "first";
  o.eps_rule = "arithmetic";
  std::vector<std::string> notes;
  const ForceField m = merge_forcefields(s.atoms.size(), {{&fa, a, "GAFF"}, {&fb, b, "OPLS"}}, o, &notes);
  // a cross pair: ε arithmetic, σ arithmetic of the two sites
  const int ta = m.type_index[a[0]], tb = m.type_index[b[0]];
  const PairType x = mixed_pair(m, ta, tb);
  EXPECT_NEAR(x.eps, 0.5 * (m.lj[size_t(ta)].eps + m.lj[size_t(tb)].eps), 1e-12);
  EXPECT_NEAR(x.sigma, 0.5 * (m.lj[size_t(ta)].sigma + m.lj[size_t(tb)].sigma), 1e-12);
  // OPLS's own unlike pairs keep OPLS's geometric σ
  if (fb.type_names.size() > 1) {
    const int t0 = m.type_index[b[0]];
    int t1 = t0;
    for (uint32_t i : b) if (m.type_index[i] != t0) { t1 = m.type_index[i]; break; }
    const PairType y = mixed_pair(m, t0, t1);
    EXPECT_NEAR(y.sigma, std::sqrt(m.lj[size_t(t0)].sigma * m.lj[size_t(t1)].sigma), 1e-12);
  }
  EXPECT_TRUE(std::any_of(notes.begin(), notes.end(), [](const std::string& n) { return n.find("1-4 scaling") != std::string::npos; }));
  EXPECT_TRUE(std::isfinite(energy(m, s).total()));
  // an atom left out
  std::vector<uint32_t> short_b(b.begin(), b.end() - 1);
  EXPECT_THROW(merge_forcefields(s.atoms.size(), {{&fa, a, "GAFF"}, {&fb, short_b, "OPLS"}}, o), FieldError);
}

// A crystal group under a literature many-body potential (Tersoff's silicon, Phys. Rev. B 37, 6991 (1988), as LAMMPS's
// potentials/Si.tersoff gives it): one type per element, no charge, zero Lennard-Jones inside, UFF across; the file's units
// checked; the LAMMPS files overlay the style and carry a copy that says its units.
TEST(FFMerge, ManyBodyGroup) {
  namespace fs = std::filesystem;
  const fs::path dir = fs::temp_directory_path() / "caps_manybody_test";
  fs::create_directories(dir);
  const std::string entry = "Si Si Si 3.0 1.0 1.3258 4.8381 2.0417 0.0000 22.956\n  0.33675 1.3258 95.373 3.0 0.2 3.2394 3264.7\n";
  { std::ofstream(dir / "tagged.tersoff") << "# DATE: 2007-10-25 UNITS: metal CITATION: Tersoff, Phys Rev B, 37, 6991 (1988)\n" << entry; }
  { std::ofstream(dir / "bare.tersoff") << "# Tersoff silicon\n" << entry; }
  // a silicon dimer below a methane
  System s;
  s.cell.a = {30, 0, 0}, s.cell.b = {0, 30, 0}, s.cell.c = {0, 0, 30};
  auto add = [&](int z, double x, double y, double w, int64_t mol) {
    Atom a;
    a.element = z, a.pos = {x, y, w}, a.mol = mol;
    s.atoms.push_back(a);
  };
  add(14, 10, 10, 10, 1), add(14, 12.35, 10, 10, 1);
  add(6, 11, 10, 14, 2), add(1, 11, 10, 15.09, 2), add(1, 12.03, 10, 13.64, 2), add(1, 10.49, 10.89, 13.64, 2), add(1, 10.49, 9.11, 13.64, 2);
  s.has_mol = true;
  s.bonds = {{0, 1, 1}, {2, 3, 1}, {2, 4, 1}, {2, 5, 1}, {2, 6, 1}};
  const std::vector<uint32_t> si = {0, 1}, me = {2, 3, 4, 5, 6};
  EXPECT_THROW(manybody_part(part_of(s, si), {"tersoff", (dir / "bare.tersoff").string(), ""}), FieldError);        // units unsaid
  EXPECT_THROW(manybody_part(part_of(s, si), {"tersoff", (dir / "tagged.tersoff").string(), "real"}), FieldError);   // units contradicted
  EXPECT_THROW(manybody_part(part_of(s, si), {"airebo", (dir / "tagged.tersoff").string(), ""}), FieldError);        // metal units only
  EXPECT_THROW(manybody_part(part_of(s, si), {"sw", (dir / "tagged.tersoff").string(), ""}), FieldError);            // not an sw file
  EXPECT_THROW(manybody_part(part_of(s, me), {"tersoff", (dir / "tagged.tersoff").string(), ""}), FieldError);       // no C or H entries
  const ForceField fs_ = manybody_part(part_of(s, si), {"tersoff", (dir / "bare.tersoff").string(), "metal"});
  EXPECT_EQ(fs_.type_names, std::vector<std::string>{"Si"});
  EXPECT_NEAR(fs_.mass[0], 28.085, 0.01);
  EXPECT_EQ(fs_.charge[0], 0.0);
  const FFDef gaff = load_forcefield(kFF + "gaff-amber25.json");
  const ForceField fm = typed(part_of(s, me), gaff);
  const ForceField m = merge_forcefields(s.atoms.size(), {{&fs_, si, "Si"}, {&fm, me, "methane"}}, MergeOptions{});
  ASSERT_TRUE(m.manybody.on());
  EXPECT_EQ(m.manybody.element[size_t(m.type_index[0])], "Si");
  EXPECT_EQ(m.manybody.element[size_t(m.type_index[2])], "");
  EXPECT_EQ(mixed_pair(m, m.type_index[0], m.type_index[0]).eps, 0.0);   // the potential does Si-Si
  double x = 0, d = 0;
  ASSERT_TRUE(uff_vdw(14, x, d));
  const PairType sc = mixed_pair(m, m.type_index[0], m.type_index[2]);
  EXPECT_NEAR(sc.eps, std::sqrt(d * m.lj[size_t(m.type_index[2])].eps), 1e-12);
  EXPECT_NEAR(sc.sigma, 0.5 * (x / std::pow(2.0, 1.0 / 6) + m.lj[size_t(m.type_index[2])].sigma), 1e-12);
  EXPECT_EQ(m.lj14, fm.lj14);   // the force field's settings, not the potential group's
  {   // beside a class II group the cross Lennard-Jones is 9-6: σ is UFF's minimum x itself
    ManyBodySpec sp{"tersoff", (dir / "tagged.tersoff").string(), ""};
    sp.pair_form = "lj9-6";
    const ForceField f96 = manybody_part(part_of(s, si), sp);
    EXPECT_EQ(f96.pair_form, "lj9-6");
    EXPECT_NEAR(f96.lj[0].sigma, x, 1e-12);
    EXPECT_NEAR(f96.lj[0].eps, d, 1e-12);
  }
  // LAMMPS: overlay, the element map, no Si-Si bond (special_bonds would hide the pair from Tersoff), the copy says its units
  EnergyOptions e;
  write_lammps_data_ff(s, m, e, (dir / "sys.data").string(), false);
  write_lammps_input(s, m, e, "sys.data", (dir / "sys.in").string());
  std::ifstream in(dir / "sys.in"), data(dir / "sys.data"), copy(dir / "bare-metal.tersoff");
  std::stringstream a, b, c;
  a << in.rdbuf(), b << data.rdbuf(), c << copy.rdbuf();
  EXPECT_NE(a.str().find(" tersoff\n"), std::string::npos);
  EXPECT_NE(a.str().find("* * tersoff bare-metal.tersoff Si NULL NULL"), std::string::npos);
  EXPECT_NE(b.str().find("4 bonds"), std::string::npos);
  EXPECT_EQ(c.str().rfind("# Tersoff silicon UNITS: metal\n", 0), 0u);
  EXPECT_THROW(write_gromacs(s, m, e, (dir / "sys").string()), FieldError);
  fs::remove_all(dir);
}

// Every potential of the library (data/potentials) reads as its style says and covers the elements it lists.
TEST(FFMerge, PotentialLibrary) {
  const std::string dir = std::string(CAPS_SOURCE_DIR) + "/data/potentials/";
  std::ifstream f(dir + "catalogue.json");
  std::stringstream ss;
  ss << f.rdbuf();
  const Json cat = Json::parse(ss.str());
  ASSERT_GE(cat["potentials"].size(), 10u);
  for (const auto& p : cat["potentials"].items()) {
    System g;
    for (const auto& e : p["elements"].items()) {
      Atom a;
      a.element = element_from_symbol(e.str());
      a.pos = {double(g.atoms.size()) * 3.0, 0, 0};
      g.atoms.push_back(a);
    }
    std::vector<std::string> notes;
    ForceField mf;
    ManyBodySpec sp{p.text("style", ""), dir + p.text("file", ""), ""};
    if (!p.text("file2", "").empty()) sp.file2 = dir + p.text("file2", "");
    if (p.has("entries")) for (const auto& [k, v] : p["entries"].members()) sp.entries.push_back({k, v.str()});
    EXPECT_NO_THROW(mf = manybody_part(g, sp, &notes)) << p.text("id", "");
    EXPECT_EQ(mf.manybody.units, "metal") << p.text("id", "");
    EXPECT_FALSE(p.text("citation", "").empty()) << p.text("id", "");
    EXPECT_EQ(mf.manybody.metal_only, p.text("style", "") == "airebo" || p.text("style", "") == "airebo/morse" || p.text("style", "") == "rebo" ||
                                          p.text("style", "") == "meam") << p.text("id", "");
  }
}

// AIREBO (read by LAMMPS in metal units only) on a carbon group: the LAMMPS files in metal units, every energy parameter
// divided by 23.060549, the times in ps; real units refused with the reason.
TEST(FFMerge, AireboWritesMetalUnits) {
  namespace fs = std::filesystem;
  const fs::path dir = fs::temp_directory_path() / "caps_airebo_test";
  fs::create_directories(dir);
  System s;
  s.cell.a = {30, 0, 0}, s.cell.b = {0, 30, 0}, s.cell.c = {0, 0, 30};
  auto add = [&](int z, double x, double y, double w, int64_t mol) {
    Atom a;
    a.element = z, a.pos = {x, y, w}, a.mol = mol;
    s.atoms.push_back(a);
  };
  add(6, 10, 10, 10, 1), add(6, 11.42, 10, 10, 1);
  add(6, 11, 10, 14, 2), add(1, 11, 10, 15.09, 2), add(1, 12.03, 10, 13.64, 2), add(1, 10.49, 10.89, 13.64, 2), add(1, 10.49, 9.11, 13.64, 2);
  s.has_mol = true;
  s.bonds = {{0, 1, 1}, {2, 3, 1}, {2, 4, 1}, {2, 5, 1}, {2, 6, 1}};
  const std::vector<uint32_t> cc = {0, 1}, me = {2, 3, 4, 5, 6};
  const std::string file = std::string(CAPS_SOURCE_DIR) + "/data/potentials/lammps/CH.airebo";
  EXPECT_THROW(manybody_part(part_of(s, cc), {"airebo", std::string(CAPS_SOURCE_DIR) + "/data/potentials/lammps/Si.tersoff", ""}), FieldError);   // not a Brenner file
  const ForceField fc = manybody_part(part_of(s, cc), {"airebo", file, ""});
  EXPECT_TRUE(fc.manybody.metal_only);
  EXPECT_EQ(fc.manybody.args, "3.0 1 1");
  const ForceField fm = typed(part_of(s, me), load_forcefield(kFF + "gaff-amber25.json"));
  const ForceField m = merge_forcefields(s.atoms.size(), {{&fc, cc, "CNT"}, {&fm, me, "methane"}}, MergeOptions{});
  EXPECT_TRUE(lammps_metal_units(m));   // automatic: metal, the published file as it is
  const ForceField mm = forcefield_in_metal_units(m);
  for (size_t t = 0; t < m.lj.size(); ++t) EXPECT_NEAR(mm.lj[t].eps, m.lj[t].eps / 23.060549, 1e-15);
  for (size_t k = 0; k < m.bonds.size(); ++k) EXPECT_NEAR(mm.bonds[k].k, m.bonds[k].k / 23.060549, 1e-12);
  EnergyOptions e;
  write_lammps_data_ff(s, m, e, (dir / "sys.data").string(), false);
  LammpsRun run;
  run.kind = LammpsRun::Kind::NPT;
  write_lammps_input(s, m, e, "sys.data", (dir / "sys.in").string(), 0, true, run);
  std::ifstream in(dir / "sys.in");
  std::stringstream a;
  a << in.rdbuf();
  const std::string t = a.str();
  EXPECT_NE(t.find("units           metal"), std::string::npos);
  EXPECT_NE(t.find(" airebo 3.0 1 1\n"), std::string::npos);
  EXPECT_NE(t.find("* * airebo CH.airebo C NULL NULL"), std::string::npos);
  EXPECT_NE(t.find("timestep        0.0005"), std::string::npos);
  EXPECT_NE(t.find("300 300 0.1 iso 1.01325 1.01325 1"), std::string::npos);   // Tdamp 100 fs, 1 atm, Pdamp 1000 fs in ps and bar
  EXPECT_TRUE(fs::exists(dir / "CH.airebo"));
  // real units, asked for: CAPS's converted copy (A_CC × 23.060549, a spline value untouched, the first line says real)
  LammpsStyle rs;
  rs.units = "real";
  EXPECT_FALSE(lammps_metal_units(m, rs));
  write_lammps_data_ff(s, m, e, (dir / "r.data").string(), false, rs);
  write_lammps_input(s, m, e, "r.data", (dir / "r.in").string(), 0, true, {}, rs);
  std::ifstream rin(dir / "r.in"), conv(dir / "CH-real.airebo");
  std::stringstream ri, rc;
  ri << rin.rdbuf(), rc << conv.rdbuf();
  EXPECT_NE(ri.str().find("units           real"), std::string::npos);
  EXPECT_NE(ri.str().find("* * airebo CH-real.airebo C NULL NULL"), std::string::npos);
  EXPECT_NE(rc.str().find("UNITS: real"), std::string::npos);
  EXPECT_NE(rc.str().find("252594.7418753"), std::string::npos);   // A_CC 10953.54416216992 eV × 23.060549
  EXPECT_NE(rc.str().find("1.7\t     rcmin_CC"), std::string::npos);   // a distance as it was
  fs::remove_all(dir);
}

// MEAM: each element mapped to a library entry (by default the first of its atomic number; a wrong element refused),
// the parameter file beside the inputs, pair_coeff with the extracted entries and a mapping per type; metal units only.
TEST(FFMerge, MeamMapsLibraryEntries) {
  namespace fs = std::filesystem;
  const fs::path dir = fs::temp_directory_path() / "caps_meam_test";
  fs::create_directories(dir);
  const std::string lib = std::string(CAPS_SOURCE_DIR) + "/data/potentials/lammps/library.meam", par = std::string(CAPS_SOURCE_DIR) + "/data/potentials/lammps/SiC.meam";
  const auto entries = meam_library(lib);
  EXPECT_EQ(entries.front().name, "AlS");
  EXPECT_EQ(std::count_if(entries.begin(), entries.end(), [](const MeamEntry& e) { return e.name == "C"; }), 1);   // the first of a repeated name
  System s;
  s.cell.a = {30, 0, 0}, s.cell.b = {0, 30, 0}, s.cell.c = {0, 0, 30};
  auto add = [&](int z, double x, double y, double w, int64_t mol) {
    Atom a;
    a.element = z, a.pos = {x, y, w}, a.mol = mol;
    s.atoms.push_back(a);
  };
  add(14, 10, 10, 10, 1), add(6, 11.9, 10, 10, 1);
  add(6, 11, 10, 14, 2), add(1, 11, 10, 15.09, 2), add(1, 12.03, 10, 13.64, 2), add(1, 10.49, 10.89, 13.64, 2), add(1, 10.49, 9.11, 13.64, 2);
  s.has_mol = true;
  s.bonds = {{2, 3, 1}, {2, 4, 1}, {2, 5, 1}, {2, 6, 1}};
  const std::vector<uint32_t> sc = {0, 1}, me = {2, 3, 4, 5, 6};
  ManyBodySpec sp{"meam", lib, ""};
  const ForceField byz = manybody_part(part_of(s, sc), sp);
  EXPECT_EQ(byz.manybody.entry, (std::vector<std::string>{"SiS", "C"}));   // the first entries of Z 14 and 6
  sp.entries = {{"Si", "C"}};
  EXPECT_THROW(manybody_part(part_of(s, sc), sp), FieldError);   // 'C' is carbon, not silicon
  sp.file2 = par;
  {   // a silicon-only group with the SiC set: SiC.meam names element 2, so C must be read too, in the file's order
    std::vector<uint32_t> si_only = {0};
    ManyBodySpec s2{"meam", lib, ""};
    s2.file2 = par;
    s2.entries = {{"Si", "Si"}};
    EXPECT_THROW(manybody_part(part_of(s, si_only), s2), FieldError);
    s2.entries = {{"Si", "Si"}, {"C", "C"}};
    const ForceField f1 = manybody_part(part_of(s, si_only), s2);
    EXPECT_EQ(f1.manybody.extract, (std::vector<std::string>{"Si", "C"}));
  }
  sp.entries = {{"Si", "Si"}, {"C", "C"}};
  const ForceField fc = manybody_part(part_of(s, sc), sp);
  const ForceField fm = typed(part_of(s, me), load_forcefield(kFF + "gaff-amber25.json"));
  const ForceField m = merge_forcefields(s.atoms.size(), {{&fc, sc, "SiC"}, {&fm, me, "methane"}}, MergeOptions{});
  LammpsStyle rs;
  rs.units = "real";
  EXPECT_THROW(lammps_metal_units(m, rs), FieldError);
  EnergyOptions e;
  write_lammps_data_ff(s, m, e, (dir / "m.data").string(), false);
  write_lammps_input(s, m, e, "m.data", (dir / "m.in").string(), 0, true);
  std::ifstream in(dir / "m.in");
  std::stringstream a;
  a << in.rdbuf();
  EXPECT_NE(a.str().find("units           metal"), std::string::npos);
  EXPECT_NE(a.str().find("* * meam library.meam Si C SiC.meam Si C NULL NULL"), std::string::npos) << a.str();
  EXPECT_TRUE(fs::exists(dir / "library.meam") && fs::exists(dir / "SiC.meam"));
  fs::remove_all(dir);
}
