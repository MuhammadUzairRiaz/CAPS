// VASP sets, readers and checks (caps/vasp_set.hpp, caps/vasp_out.hpp, caps/dft_analysis.hpp): k-meshes, POTCAR order,
// the job script (syntax; a hung first attempt restarted by the watchdog; stages stored with links; resubmission skips
// DONE stages), OUTCAR reading, the health check's flags, and the convergence choice.
#include <gtest/gtest.h>

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <regex>
#include <set>

#include "caps/dft_analysis.hpp"
#include "caps/slab2d.hpp"
#include "caps/vasp_out.hpp"
#include "caps/vasp_set.hpp"

using namespace caps;
namespace fs = std::filesystem;

namespace {
const std::string kData = std::string(CAPS_SOURCE_DIR) + "/data";
const std::string kRef = std::string(CAPS_SOURCE_DIR) + "/tests/data/slab2d/";
const std::string kVasp = std::string(CAPS_SOURCE_DIR) + "/tests/data/vasp/";
void spit(const fs::path& p, const std::string& t) { fs::create_directories(p.parent_path()); std::ofstream(p) << t; }
std::string slurp(const fs::path& p) { std::ifstream f(p); return std::string((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>()); }
fs::path tmpdir(const std::string& name) { const auto d = fs::temp_directory_path() / name; fs::remove_all(d); fs::create_directories(d); return d; }
System hex(double a) {
  System s;
  s.cell.a = {a, 0, 0};
  s.cell.b = {-a / 2, a * std::sqrt(3.0) / 2, 0};
  s.cell.c = {0, 0, 30};
  return s;
}
}  // namespace

TEST(DftSets, KpointMeshes) {
  for (double a : {3.0166, 3.0589, 3.0513}) {
    const auto k = kmesh(hex(a), 45.0);
    EXPECT_EQ(k[0], 15) << a;
    EXPECT_EQ(k[1], 15) << a;
    EXPECT_EQ(k[2], 1);
  }
  EXPECT_EQ(kmesh(hex(3 * 3.0813748837), 45.0)[0], 5);   // the 3 × 3 mixed cell
  EXPECT_EQ(kmesh(hex(5 * 3.0513), 45.0)[0], 3);          // 5 × 5 adsorption cells
  EXPECT_EQ(kmesh(hex(5 * 3.0513), 45.0, true)[0], 1);    // Γ only
}

TEST(DftSets, PotcarOrderChecked) {
  // fake PAW files (only the lines CAPS reads; never real POTCARs)
  const auto pp = tmpdir("caps_fake_paw");
  for (const auto& [ds, el, enmax, z] : std::vector<std::tuple<std::string, std::string, double, double>>{{"Ti_pv", "Ti_pv", 222.335, 10}, {"C", "C", 400, 4}, {"O", "O", 400, 6}})
    spit(pp / ds / "POTCAR", "  PAW_PBE " + el + " 07Sep2000\n   TITEL  = PAW_PBE " + el + " 07Sep2000\n   ENMAX  =  " + std::to_string(enmax) + "; ENMIN = 1.0 eV\n   POMASS = 1.0; ZVAL   =   " + std::to_string(z) + "    mass and valenz\n End of Dataset\n");
  const auto map = potcar_map(kData);
  const auto out = (pp / "POTCAR").string();
  const auto p = assemble_potcar({"Ti", "C", "O"}, pp.string(), map, out);
  EXPECT_DOUBLE_EQ(p.encut, 520);   // 1.3 × 400
  EXPECT_NO_THROW(verify_potcar(out, {"Ti", "C", "O"}));
  EXPECT_THROW(verify_potcar(out, {"Ti", "O", "C"}), std::runtime_error);   // permuted
}

TEST(DftSets, SlabSetFiles) {
  const auto dir = tmpdir("caps_slab_set");
  std::map<std::string, std::string> bad;
  const System s = read_vasp_poscar(kRef + "Ti3C2_O.POSCAR");
  VaspSetOptions o;
  o.label = "Ti3C2_O";
  o.stages = stages_slab();
  o.data_dir = kData;
  o.settings.static_out = "pot";
  const auto r = write_vasp_set(s, dir.string(), o);
  EXPECT_DOUBLE_EQ(r.encut, 520);
  EXPECT_EQ(r.kpoints[0], 15);
  const std::string cell = slurp(dir / "INCAR.cell"), stat = slurp(dir / "INCAR.static");
  EXPECT_NE(cell.find("LATTICE_CONSTRAINTS = .TRUE. .TRUE. .FALSE."), std::string::npos) << cell;
  EXPECT_NE(cell.find("ISIF    = 3"), std::string::npos);
  EXPECT_NE(stat.find("LVHAR   = .TRUE."), std::string::npos);
  EXPECT_NE(stat.find("LCHARG  = .FALSE."), std::string::npos);
  EXPECT_EQ(stat.find("LDIPOL"), std::string::npos);
  EXPECT_NE(stat.find("KPAR    = 4"), std::string::npos);   // ≤ 20 atoms, slurm-workspace profile: 52 ranks, KPAR 4, NCORE 13
  EXPECT_NE(stat.find("NCORE   = 13"), std::string::npos);
  EXPECT_EQ(std::system(("bash -n " + (dir / "job.slurm").string()).c_str()), 0);
  EXPECT_EQ(std::system(("bash -n " + (dir / "make_potcar.sh").string()).c_str()), 0);
  EXPECT_EQ(read_vasp_poscar((dir / "POSCAR").string()).atoms.size(), 7u);
}

TEST(DftSets, JobScriptWatchdogStoreAndResume) {
  if (std::system("command -v bash > /dev/null 2>&1") != 0) GTEST_SKIP() << "no bash";
  const auto root = tmpdir("caps_job_test");
  const auto proj = root / "proj", store = root / "store", cs = proj / "slab";
  spit(proj / ".store", store.string() + "\n");
  const System s = read_vasp_poscar(kRef + "Ti3C2_O.POSCAR");
  VaspSetOptions o;
  o.label = "slab";
  o.data_dir = kData;
  o.profile = "generic-slurm";
  write_vasp_set(s, cs.string(), o);
  spit(cs / "POTCAR", "fake\n");
  // a fake srun: the first 03_relax attempt hangs (writes vasp.out once, then nothing); every other run finishes
  const auto fake = root / "fake_srun.sh";
  spit(fake, "#!/bin/bash\n"
             "echo running > vasp.out\n"
             "if grep -q 'ISIF    = 2' INCAR && [ ! -f \"$MARK\" ]; then touch \"$MARK\"; sleep 30; exit 1; fi\n"
             "cp POSCAR CONTCAR\n"
             "echo '   1 F= -.1E+03 E0= -.1E+03  d E =-.1E+03' > OSZICAR\n"
             "if grep -q 'NSW     = 0' INCAR; then echo ' EDIFF is reached' > OUTCAR; else echo ' reached required accuracy' > OUTCAR; fi\n");
  fs::permissions(fake, fs::perms::owner_all);
  std::string job = slurp(cs / "job.slurm");
  job = std::regex_replace(job, std::regex("HANG=\\d+"), "HANG=2");
  spit(cs / "job.slurm", job);
  const std::string env = "cd " + cs.string() + " && SLURM_SUBMIT_DIR=" + cs.string() + " SLURM_JOB_ID=7 CAPS_POLL=1 CAPS_SRUN=" + fake.string() + " MARK=" + (root / "hung").string() +
                          " bash job.slurm > log1.txt 2>&1";
  ASSERT_EQ(std::system(env.c_str()), 0) << slurp(cs / "log1.txt");
  const std::string log = slurp(cs / "log1.txt");
  EXPECT_NE(log.find("VASP hung"), std::string::npos) << log;
  EXPECT_NE(log.find("attempt 2"), std::string::npos) << log;
  EXPECT_TRUE(fs::exists(cs / "03_relax" / "DONE")) << log;
  EXPECT_TRUE(fs::exists(cs / "04_static" / "DONE")) << log;
  EXPECT_TRUE(fs::is_symlink(cs / "03_relax"));
  EXPECT_TRUE(fs::exists(store / "slab" / "04_static" / "OUTCAR"));
  // resubmitted: both stages skipped
  const std::string env2 = "cd " + cs.string() + " && SLURM_SUBMIT_DIR=" + cs.string() + " SLURM_JOB_ID=8 CAPS_POLL=1 CAPS_SRUN=" + fake.string() + " MARK=" + (root / "hung").string() +
                           " bash job.slurm > log2.txt 2>&1";
  ASSERT_EQ(std::system(env2.c_str()), 0);
  const std::string log2 = slurp(cs / "log2.txt");
  EXPECT_NE(log2.find("03_relax already done"), std::string::npos) << log2;
  EXPECT_NE(log2.find("04_static already done"), std::string::npos) << log2;
  EXPECT_TRUE(submittable({cs.string()}, {}).empty());
}

TEST(DftSets, OutcarRead) {
  const Outcar o = read_outcar(kVasp + "OUTCAR_static");
  EXPECT_TRUE(o.finished);
  EXPECT_TRUE(o.ediff_reached);
  EXPECT_FALSE(o.scf_failed);
  EXPECT_EQ(o.nions, 3);
  ASSERT_EQ(o.e0.size(), 1u);
  EXPECT_DOUBLE_EQ(o.e0[0], -51.232);
  ASSERT_EQ(o.fmax.size(), 1u);
  EXPECT_NEAR(o.fmax[0], 0.005, 1e-12);
  ASSERT_EQ(o.stress.size(), 1u);
  EXPECT_DOUBLE_EQ(o.stress[0][1], 1.3);
  EXPECT_DOUBLE_EQ(o.efermi.back(), -1.2345);
  EXPECT_DOUBLE_EQ(o.dipole.back(), -0.012345);
  EXPECT_DOUBLE_EQ(o.encut, 520);
  EXPECT_EQ(o.nelm, 120);
  EXPECT_EQ(o.nsw, 0);
}

TEST(DftSets, HealthCheckFlags) {
  // a set with a normal complex and one whose energy jumped by 245 eV with a 1.6 e·Å dipole
  const auto root = tmpdir("caps_health");
  const System slab = read_vasp_poscar(kRef + "Ti3C2_OH.POSCAR");
  write_vasp_poscar(slab, (root / "slab" / "POSCAR").string(), "slab");
  System mol;
  mol.cell = slab.cell;
  Atom h1, h2;
  h1.element = h2.element = 7;
  h1.pos = {0, 0, 0}, h2.pos = {1.10, 0, 0};
  mol.atoms = {h1, h2};
  write_vasp_poscar(mol, (root / "molecule" / "POSCAR").string(), "N2");
  auto outcar = [](double e0, double dip) {
    return " NIONS =      2\n NELM   =    120;\n energy  without entropy=  0  energy(sigma->0) =  " + std::to_string(e0) + "\n dipolmoment  0.0 0.0 " + std::to_string(dip) +
           " electrons x Angstroem\n EDIFF is reached\n General timing and accounting\n";
  };
  spit(root / "slab" / "04_static" / "OUTCAR", outcar(-100.0, 0.0));
  spit(root / "molecule" / "04_static" / "OUTCAR", outcar(-6.8, 0.0));
  for (const auto& [name, e, jump, dip] : std::vector<std::tuple<std::string, double, double, double>>{{"complex_ok_az0", -107.3, 0.0, 0.2}, {"complex_bad_az0", 138.0, 245.0, 1.6}}) {
    System cx = slab;
    for (auto a : mol.atoms) { a.pos = a.pos + Vec3{1.0, 1.0, 22.5}; cx.atoms.push_back(a); }
    write_vasp_poscar(cx, (root / name / "POSCAR").string(), name);
    spit(root / name / "03_relax" / "OUTCAR", outcar(e, dip));
    spit(root / name / "03_relax" / "OSZICAR", "   1 F= " + std::to_string(e - jump) + " E0= " + std::to_string(e - jump) + "\n   2 F= " + std::to_string(e) + " E0= " + std::to_string(e) + "\n");
  }
  const Json h = health_check(root.string());
  ASSERT_FALSE(h.has("error")) << h.dump(1);
  std::map<std::string, std::set<std::string>> kinds;
  for (const auto& r : h["complexes"].items()) for (const auto& f : r["flags"].items()) kinds[r.text("complex")].insert(f.text("kind"));
  EXPECT_TRUE(kinds["complex_bad_az0"].count("ENERGY")) << h.dump(1);
  EXPECT_TRUE(kinds["complex_bad_az0"].count("JUMP"));
  EXPECT_TRUE(kinds["complex_bad_az0"].count("DIPOLE"));
  EXPECT_TRUE(kinds["complex_ok_az0"].empty()) << h.dump(1);
}

TEST(DftSets, ConvergenceChoice) {
  // ENCUT: stress settles (≤ 1 kB of 800 eV) from 520 on, forces everywhere within 0.01 eV/Å
  std::vector<ConvRow> en;
  for (const auto& [v, e, s, f] : std::vector<std::tuple<int, double, double, double>>{
           {400, -9.10, -12.0, 0.020}, {450, -9.13, -4.5, 0.015}, {520, -9.140, -1.6, 0.012}, {600, -9.141, -1.1, 0.011}, {700, -9.1412, -0.9, 0.011}, {800, -9.1413, -0.8, 0.011}}) {
    ConvRow r;
    r.value = std::to_string(v), r.e_atom = e, r.stress = s, r.fmax = f;
    en.push_back(r);
  }
  EXPECT_EQ(convergence_choice(en, true), "520");
  // k: E0 within 1 meV/atom and stress within 1 kB from 15 × 15
  std::vector<ConvRow> k;
  for (const auto& [v, e, s] : std::vector<std::tuple<int, double, double>>{{6, -9.120, 3.0}, {9, -9.136, 1.8}, {12, -9.1395, 1.2}, {15, -9.1408, 0.6}, {18, -9.1410, 0.5}}) {
    ConvRow r;
    r.value = std::to_string(v), r.e_atom = e, r.stress = s, r.fmax = 0.01;
    k.push_back(r);
  }
  EXPECT_EQ(convergence_choice(k, false), "15");
}

TEST(DftSets, LatticeScanFit) {
  // E(a) = (a − 3.05)² × 20 − 50 sampled at five points: the parabola and Birch–Murnaghan minima at 3.05 Å
  const System s = hex(3.0);
  std::vector<std::pair<double, double>> pts;
  for (double a : {2.98, 3.02, 3.06, 3.10, 3.14}) pts.push_back({a, 20 * (a - 3.05) * (a - 3.05) - 50});
  const auto f = fit_scan_points(pts, s, 3.0);
  EXPECT_NEAR(f.a_parabola, 3.05, 1e-9);
  EXPECT_TRUE(f.inside);
  EXPECT_NEAR(f.a_bm, 3.05, 0.005);
}
