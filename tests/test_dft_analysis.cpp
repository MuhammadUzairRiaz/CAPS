// The DFT analyses on small synthetic VASP files (written here in VASP's formats, the lines CAPS reads): Bader transfer,
// grouped DOS, the work function from LOCPOT, binding energies, AIMD stability.
#include <gtest/gtest.h>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <sstream>

#include "caps/dft_analysis.hpp"
#include "caps/slab2d.hpp"

using namespace caps;
namespace fs = std::filesystem;

namespace {
const std::string kRef = std::string(CAPS_SOURCE_DIR) + "/tests/data/slab2d/";
void spit(const fs::path& p, const std::string& t) { fs::create_directories(p.parent_path()); std::ofstream(p) << t; }
fs::path tmpdir(const std::string& n) { const auto d = fs::temp_directory_path() / n; fs::remove_all(d); fs::create_directories(d); return d; }
// slab + an N2 molecule 3 Å above its top
System with_n2(const System& slab) {
  System s = slab;
  double top = -1e9;
  for (const auto& a : slab.atoms) top = std::max(top, a.pos[2]);
  for (double x : {0.0, 1.10}) { Atom a; a.element = 7; a.pos = {1.0 + x, 1.0, top + 3.0}; s.atoms.push_back(a); }
  return s;
}
std::string outcar(double e0, double ef = 0) {
  std::ostringstream o;
  o << " NIONS =      2\n NELM   =    120;\n E-fermi :  " << ef << "   XC(G=0): 0\n energy  without entropy=  0  energy(sigma->0) =  " << e0
    << "\n EDIFF is reached\n General timing and accounting\n";
  return o.str();
}
}  // namespace

TEST(DftAnalysis, BaderTransfer) {
  const auto d = tmpdir("caps_bader");
  const System slab = read_vasp_poscar(kRef + "Ti3C2_OH.POSCAR");
  const System cx = with_n2(slab);
  write_vasp_poscar(slab, (d / "slab.vasp").string(), "slab");
  write_vasp_poscar(cx, (d / "cx.vasp").string(), "complex");
  // ZVAL Ti 10, C 4, O 6, H 1, N 5: the molecule gains 0.2 e from one O of the slab
  std::ostringstream acf;
  acf << "    #         X           Y           Z       CHARGE      MIN DIST   ATOMIC VOL\n ------\n";
  int k = 1;
  for (const auto& a : cx.atoms) {
    double q = a.element == 22 ? 10 : a.element == 6 ? 4 : a.element == 8 ? 6 : a.element == 1 ? 1 : 5.1;
    if (k == 6) q -= 0.2;   // the first O
    acf << "    " << k++ << "    0.0 0.0 0.0    " << q << "    1.0    10.0\n";
  }
  spit(d / "ACF.dat", acf.str());
  const Json r = bader_transfer((d / "slab.vasp").string(), (d / "cx.vasp").string(), (d / "ACF.dat").string());
  EXPECT_NEAR(r["molecule_gained"].number(), 0.2, 1e-9);
  EXPECT_NEAR(r["slab_gained"].number(), -0.2, 1e-9);
  EXPECT_NEAR(r["total"].number(), 0.0, 1e-9);
}

TEST(DftAnalysis, WorkFunctionFromLocpot) {
  // V(z) = 5 eV in the vacuum, 0 inside the slab; E_F = −1 eV → φ = 6 eV on both faces, flat plateau
  const auto d = tmpdir("caps_wf");
  System s = read_vasp_poscar(kRef + "Ti3C2_O.POSCAR");
  write_vasp_poscar(s, (d / "POSCAR").string(), "x");
  std::ifstream pin((d / "POSCAR").string());
  std::string head, l;
  for (int i = 0; i < 8 && std::getline(pin, l); ++i) head += (i == 7 ? "Direct" : l) + "\n";
  std::ostringstream g;
  g << head;
  for (const auto& a : s.atoms) { const Vec3 f = s.cell.to_fractional(a.pos); g << f[0] << " " << f[1] << " " << f[2] << "\n"; }
  const int nz = 60;
  const double c = s.cell.c[2];
  g << "\n 2 2 " << nz << "\n";
  double zmin = 1e9, zmax = -1e9;
  for (const auto& a : s.atoms) zmin = std::min(zmin, a.pos[2]), zmax = std::max(zmax, a.pos[2]);
  for (int k = 0; k < nz; ++k) {
    const double z = c * k / nz;
    const double v = (z > zmax + 2 || z < zmin - 2) ? 5.0 : 0.0;
    for (int q = 0; q < 4; ++q) g << " " << v;
    g << "\n";
  }
  spit(d / "LOCPOT", g.str());
  spit(d / "OUTCAR", outcar(-50, -1.0));
  const Json r = work_function(d.string());
  ASSERT_FALSE(r.has("error")) << r.dump(1);
  EXPECT_NEAR(r["phi_top"].number(), 6.0, 1e-9);
  EXPECT_NEAR(r["phi_bottom"].number(), 6.0, 1e-9);
  EXPECT_NEAR(r["plateau_flatness_top"].number(), 0.0, 1e-12);
}

TEST(DftAnalysis, GroupedDos) {
  const auto d = tmpdir("caps_dos") / "slab" / "04_static";
  fs::create_directories(d);
  System s = read_vasp_poscar(kRef + "Ti3C2_O.POSCAR");
  write_vasp_poscar(s, (d / "POSCAR").string(), "x");
  const int n = int(s.atoms.size()), nedos = 5;
  std::ostringstream o;
  o << n << " " << n << " 1 0\n x\n x\n CAR\n x\n 5.0 -5.0 " << nedos << " 1.0 1.0\n";
  for (int e = 0; e < nedos; ++e) o << (-1.0 + e) << " " << (e == 2 ? 3.0 : 1.0) << " 0\n";   // E_F = 1 → E − E_F = −2 … 2
  for (int a = 0; a < n; ++a) {
    o << " 5.0 -5.0 " << nedos << " 1.0 1.0\n";
    for (int e = 0; e < nedos; ++e) o << (-1.0 + e) << " 0.1 0.1 0.1\n";
  }
  spit(d / "DOSCAR", o.str());
  const Json r = dos_analysis(d.string());
  ASSERT_FALSE(r.has("error")) << r.dump(1);
  EXPECT_DOUBLE_EQ(r["DOS_at_EF_total"].number(), 3.0);   // E − E_F = 0 at the third energy
  EXPECT_TRUE(r["groups"].has("Ti"));
  EXPECT_NEAR(r["groups"]["Ti"][0].number(), 3 * 0.3, 1e-12);   // three Ti, three orbitals of 0.1
}

TEST(DftAnalysis, BindingEnergies) {
  const auto root = tmpdir("caps_bind");
  const std::string incar = "ENCUT = 520\nIVDW = 12\nISMEAR = 0\nSIGMA = 0.05\nPREC = Accurate\nIDIPOL = 3\nLREAL = Auto\nGGA = PE\n";
  for (const auto& [name, e] : std::vector<std::pair<std::string, double>>{{"slab", -100.0}, {"molecule", -6.8}, {"complex_a_az0", -107.3}, {"complex_b_az0", -107.0}}) {
    spit(root / name / "04_static" / "OUTCAR", outcar(e));
    spit(root / name / "INCAR.static", incar);
  }
  spit(root / "complex_b_az0" / "INCAR.static", incar + "ENCUT = 400\n");   // a different setting is reported
  const Json r = binding_energies(root.string());
  ASSERT_FALSE(r.has("error")) << r.dump(1);
  EXPECT_EQ(r.text("best"), "complex_a_az0");
  EXPECT_NEAR(r["best_E_bind"].number(), -0.5, 1e-9);
  EXPECT_NEAR(r["complexes"][0]["E_bind_kJmol"].number(), -0.5 * 96.485, 1e-6);
  EXPECT_TRUE(r["complexes"][0]["final"].boolean());
  EXPECT_EQ(r["settings_warnings"].items().size(), 1u);
  EXPECT_TRUE(fs::exists(root / "binding_energies.csv"));
}

TEST(DftAnalysis, AimdStaysAdsorbed) {
  // twelve frames of a slab with N2 3 Å above it, sitting still at 300 K: no broken bonds, height constant → STAYS ADSORBED
  const auto d = tmpdir("caps_md");
  const System s = with_n2(read_vasp_poscar(kRef + "Ti3C2_OH.POSCAR"));
  write_vasp_poscar(s, (d / "P").string(), "x");
  std::ifstream pin((d / "P").string());
  std::string head, l;
  for (int i = 0; i < 7 && std::getline(pin, l); ++i) head += l + "\n";
  std::ostringstream x, z;
  x << head;
  for (int f = 1; f <= 12; ++f) {
    x << "Direct configuration=" << f << "\n";
    for (const auto& a : s.atoms) { const Vec3 fr = s.cell.to_fractional(a.pos); x << " " << fr[0] << " " << fr[1] << " " << fr[2] << "\n"; }
    z << "     " << f << " T=   " << (300 + (f % 3) - 1) << ". E= -.10000000E+03 F= -.10500000E+03 E0= -.10500000E+03  EK= 0.1E+00 SP= 0 SK= 0\n";
  }
  spit(d / "seg_001" / "XDATCAR", x.str());
  spit(d / "seg_001" / "OSZICAR", z.str());
  const Json r = md_analysis(d.string(), 0.0, 1.0);
  EXPECT_EQ(r.text("verdict"), "STAYS ADSORBED") << r.dump(1);
  EXPECT_EQ(int(r["frames"].number()), 12);
  EXPECT_NEAR(r["T_mean"].number(), 300.0, 1.0);
  EXPECT_EQ(int(r["broken_bonds_max"].number()), 0);
}
