// The VASP relaxation chain and the phono3py workflow (caps/vasp_set.hpp, caps/phono3py.hpp): the guide's ENCUT / ISIF
// schedule, k-meshes, 2D handling (in-plane relaxation, 1 along c, dipole), the supercell rule, phono3py's displacement
// count (the guide's PE: 8367 fc3 and 15 fc2 supercells), the step scripts (syntax, and the whole chain run with fake
// VASP and phono3py), the check and the κ analysis on a synthetic extract.
#include <gtest/gtest.h>

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <regex>

#include "caps/crystal.hpp"
#include "caps/dft_commands.hpp"
#include "caps/phono3py.hpp"
#include "caps/provenance.hpp"
#include "caps/recipe.hpp"
#include "caps/slab2d.hpp"
#include "caps/vasp_set.hpp"

using namespace caps;
namespace fs = std::filesystem;

namespace {
const std::string kData = std::string(CAPS_SOURCE_DIR) + "/data";
const std::string kPoly = kData + "/crystals/polymers/";
const std::string kSheet = std::string(CAPS_SOURCE_DIR) + "/tests/data/slab2d/Ti3C2_O.POSCAR";
void spit(const fs::path& p, const std::string& t) { fs::create_directories(p.parent_path()); std::ofstream(p) << t; }
std::string slurp(const fs::path& p) { std::ifstream f(p); return std::string((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>()); }
fs::path tmpdir(const std::string& name) { const auto d = fs::temp_directory_path() / name; fs::remove_all(d); fs::create_directories(d); return d; }
System pe() { return group_by_species(read_cif(kPoly + "pe-orthorhombic.cif"), {"C", "H"}); }
std::string tag(const std::string& incar, const std::string& key) {
  std::smatch m;
  return std::regex_search(incar, m, std::regex("(^|\\n)" + key + "\\s*=\\s*([^#\\n]*?)\\s*(#|\\n)")) ? m[2].str() : "";
}
bool bash_ok(const fs::path& p) { return std::system(("bash -n '" + p.string() + "'").c_str()) == 0; }
}  // namespace

TEST(RelaxChain, GuideScheduleForPE) {
  const auto dir = tmpdir("caps_relax_pe");
  RelaxChainOptions o;
  o.label = "pe";
  o.data_dir = kData;
  o.dimensionality = "3d";
  const auto r = write_relax_chain(pe(), dir.string(), o);
  // ENMAX 400 (C) → 2 × and 3 × rounded to 10 eV; relax3 and relax4 back at relax1's
  ASSERT_EQ(r.steps.size(), 4u);
  const double enc[] = {800, 1200, 800, 800};
  const int isif[] = {3, 3, 4, 2};
  for (int k = 0; k < 4; ++k) {
    EXPECT_DOUBLE_EQ(r.steps[size_t(k)].encut, enc[k]);
    EXPECT_EQ(r.steps[size_t(k)].isif, isif[k]);
    const std::string in = slurp(dir / r.steps[size_t(k)].name / "INCAR");
    EXPECT_EQ(tag(in, "ENCUT"), std::to_string(int(enc[k]))) << in;
    EXPECT_EQ(tag(in, "ISIF"), std::to_string(isif[k]));
    EXPECT_EQ(tag(in, "IVDW"), "12");
    EXPECT_EQ(tag(in, "IBRION"), "2");
    EXPECT_EQ(tag(in, "LATTICE_CONSTRAINTS"), "");   // 3D: the whole cell relaxes
    EXPECT_EQ(tag(in, "IDIPOL"), "");
  }
  // k: ceil(50 / |a_i|) for a 7.40 × 4.93 × 2.534 Å cell
  EXPECT_EQ(r.kpoints, (std::array<int, 3>{7, 11, 20}));
  EXPECT_NE(slurp(dir / "KPOINTS").find("7  11  20"), std::string::npos);
  EXPECT_TRUE(bash_ok(dir / "job.slurm"));
  EXPECT_EQ(read_vasp_poscar((dir / "POSCAR").string()).atoms.size(), 12u);
  // the guide's own PE example: 900 / 1350 / 900 / 900 eV
  o.encuts = {900, 1350, 900, 900};
  const auto g = write_relax_chain(pe(), (dir / "guide").string(), o);
  EXPECT_DOUBLE_EQ(g.steps[1].encut, 1350);
  EXPECT_DOUBLE_EQ(g.steps[3].encut, 900);
}

TEST(RelaxChain, PotcarEnmaxGivesTheGuideCutoffs) {
  // a fake PAW directory whose C ENMAX is 450: 2 × and 3 × give the guide's 900 / 1350
  const auto pp = tmpdir("caps_relax_paw");
  for (const auto& [el, enmax, z] : std::vector<std::tuple<std::string, double, double>>{{"C", 450, 4}, {"H", 250, 1}})
    spit(pp / el / "POTCAR", "  PAW_PBE " + el + " 08Apr2002\n   TITEL  = PAW_PBE " + el + " 08Apr2002\n   ENMAX  =  " + std::to_string(enmax) + "; ENMIN = 1.0 eV\n   POMASS = 1.0; ZVAL   =   " + std::to_string(z) + "    mass and valenz\n End of Dataset\n");
  RelaxChainOptions o;
  o.label = "pe";
  o.data_dir = kData;
  o.pp_dir = pp.string();
  const auto r = write_relax_chain(pe(), (pp / "case").string(), o);
  EXPECT_DOUBLE_EQ(r.enmax, 450);
  EXPECT_DOUBLE_EQ(r.steps[0].encut, 900);
  EXPECT_DOUBLE_EQ(r.steps[1].encut, 1350);
  EXPECT_TRUE(fs::exists(pp / "case" / "POTCAR"));
}

TEST(RelaxChain, ChainRunsInOrderAndStopsOnFailure) {
  if (std::system("command -v bash > /dev/null 2>&1") != 0) GTEST_SKIP() << "no bash";
  const auto root = tmpdir("caps_relax_run");
  const auto cs = root / "pe";
  RelaxChainOptions o;
  o.label = "pe";
  o.data_dir = kData;
  o.profile = "generic-slurm";
  write_relax_chain(pe(), cs.string(), o);
  spit(cs / "POTCAR", "fake\n");
  // fake VASP: the cell shrinks 1 % per step (CONTCAR = POSCAR scaled), relax3 fails while FAIL3 exists
  const auto fake = root / "fake_srun.sh";
  spit(fake, "#!/bin/bash\n"
             "if grep -q 'ISIF    = 4' INCAR && [ -f \"$FAIL3\" ]; then echo 'not converged' > OUT; cp POSCAR CONTCAR; exit 0; fi\n"
             "awk 'NR == 2 {print 0.99 * $1; next} {print}' POSCAR > CONTCAR\n"
             "echo '   1 F= -.7E+02 E0= -.7E+02  d E =-.1E+00' > OSZICAR\n"
             "echo ' reached required accuracy - stopping structural energy minimisation' > OUT\n");
  fs::permissions(fake, fs::perms::owner_all);
  spit(root / "fail3", "x");
  const std::string env = "cd " + cs.string() + " && SLURM_SUBMIT_DIR=" + cs.string() + " SLURM_JOB_ID=1 CAPS_SRUN=" + fake.string() + " FAIL3=" + (root / "fail3").string();
  EXPECT_NE(std::system((env + " bash job.slurm > log1.txt 2>&1").c_str()), 0);
  const std::string log1 = slurp(cs / "log1.txt");
  EXPECT_TRUE(fs::exists(cs / "relax1" / "DONE") && fs::exists(cs / "relax2" / "DONE")) << log1;
  EXPECT_FALSE(fs::exists(cs / "relax3" / "DONE"));
  EXPECT_FALSE(fs::exists(cs / "relax4"/ "OUT"));
  EXPECT_NE(log1.find("relax3 did NOT reach the required accuracy"), std::string::npos) << log1;
  EXPECT_NE(log1.find("volume"), std::string::npos) << log1;   // the cell change of each finished step
  // relax2 started from relax1's CONTCAR, unchanged
  EXPECT_EQ(slurp(cs / "relax2" / "POSCAR.start"), slurp(cs / "relax1" / "CONTCAR"));
  // fixed, resubmitted: relax1 and relax2 skipped, the chain finishes and POSCAR-unitcell is relax4's CONTCAR
  fs::remove(root / "fail3");
  fs::remove(cs / "relax3" / "CONTCAR");
  ASSERT_EQ(std::system((env + " bash job.slurm > log2.txt 2>&1").c_str()), 0) << slurp(cs / "log2.txt");
  const std::string log2 = slurp(cs / "log2.txt");
  EXPECT_NE(log2.find("relax1 already done"), std::string::npos) << log2;
  EXPECT_EQ(slurp(cs / "POSCAR-unitcell"), slurp(cs / "relax4" / "CONTCAR"));
  // the report: four steps done, each 1 % shorter (volume −2.97 %)
  const Json rep = relax_chain_report(cs.string());
  EXPECT_TRUE(rep["finished"].boolean());
  EXPECT_NEAR(rep["steps"].items()[0].num("volume_change_percent", 0), 100 * (0.99 * 0.99 * 0.99 - 1), 1e-6);
}

TEST(RelaxChain, TwoDimensionalSheet) {
  const System s = read_vasp_poscar(kSheet);   // Ti3C2O2 with vacuum along c
  // no choice: the vacuum gap is a hint, never a silent switch
  EXPECT_EQ(dimensionality_hint(s).guess, "2d");
  EXPECT_THROW(resolve_dimensionality(s, ""), std::invalid_argument);
  EXPECT_EQ(resolve_dimensionality(pe(), ""), "3d");
  EXPECT_THROW(resolve_dimensionality(pe(), "2d"), std::invalid_argument);   // no vacuum
  const auto dir = tmpdir("caps_relax_2d");
  RelaxChainOptions o;
  o.label = "Ti3C2O2";
  o.data_dir = kData;
  o.dimensionality = "2d";
  const auto r = write_relax_chain(s, dir.string(), o);
  EXPECT_EQ(r.dimensionality, "2d");
  EXPECT_EQ(r.kpoints[2], 1);
  EXPECT_FALSE(r.polar);   // the same O face on both sides
  for (const auto& st : r.steps) {
    const std::string in = slurp(dir / st.name / "INCAR");
    EXPECT_EQ(tag(in, "LATTICE_CONSTRAINTS"), st.isif >= 3 ? ".TRUE. .TRUE. .FALSE." : "") << st.name;
    EXPECT_EQ(tag(in, "IDIPOL"), "");
  }
  // a polar sheet (the top O removed): the dipole correction with LDIPOL
  System p = s;
  size_t top = 0;
  for (size_t i = 0; i < p.atoms.size(); ++i) if (p.atoms[i].pos[2] > p.atoms[top].pos[2]) top = i;
  p.atoms.erase(p.atoms.begin() + long(top));
  EXPECT_TRUE(polar_along_c(p));
  const auto rp = write_relax_chain(p, (dir / "polar").string(), o);
  const std::string in = slurp(dir / "polar" / "relax1" / "INCAR");
  EXPECT_EQ(tag(in, "IDIPOL"), "3");
  EXPECT_EQ(tag(in, "LDIPOL"), ".TRUE.");
}

TEST(Phono3py, ChainAxesAndSupercellRule) {
  EXPECT_EQ(chain_axes(pe()), std::vector<int>{2});
  // fc3: at least 10 Å (5 % short accepted), fc2: 15 Å; the latter is the guide's PE fc2 supercell (its 6 3 2, chain along a)
  EXPECT_EQ(supercell_for_length(pe(), 10.0), (std::array<int, 3>{2, 2, 4}));
  EXPECT_EQ(supercell_for_length(pe(), 15.0), (std::array<int, 3>{2, 3, 6}));
  const System nr = read_cif(kPoly + "natural-rubber.cif");
  EXPECT_EQ(chain_axes(nr), std::vector<int>{2});
  // PTFE phase IV (one 15/7 helix per hexagonal cell, c 19.50 A): the chain along c, one cell along it is already 10 A
  const System ptfe = read_cif(kPoly + "ptfe-iv.cif");
  EXPECT_EQ(chain_axes(ptfe), std::vector<int>{2});
  EXPECT_EQ(supercell_for_length(ptfe, 10.0), (std::array<int, 3>{2, 2, 1}));
  EXPECT_EQ(supercell_for_length(ptfe, 15.0), (std::array<int, 3>{3, 3, 1}));
  EXPECT_EQ(crystal_operations(ptfe).size(), 6u);   // P3112: P31 and the helix's twofold axes
}

TEST(RelaxChain, PtfeAccepted) {
  // the generated PTFE cells go through the relaxation chain like the papers' cells: C and F at ENMAX 400, k from the cell
  for (const char* id : {"ptfe-iv", "ptfe-ii"}) {
    const auto dir = tmpdir(std::string("caps_relax_") + id);
    RelaxChainOptions o;
    o.label = id;
    o.data_dir = kData;
    o.dimensionality = "3d";
    const System s = group_by_species(read_cif(kPoly + id + ".cif"), {"C", "F"});
    const auto r = write_relax_chain(s, dir.string(), o);
    ASSERT_EQ(r.steps.size(), 4u);
    EXPECT_DOUBLE_EQ(r.steps[0].encut, 800);
    EXPECT_EQ(r.kpoints, (std::array<int, 3>{9, 9, 3}));
    EXPECT_EQ(read_vasp_poscar((dir / "POSCAR").string()).atoms.size(), s.atoms.size());
  }
}

TEST(Phono3py, GuidePolyethyleneDisplacementCount) {
  // the guide's PE (7.2 / 7.3): --dim "3 2 2" --dim-fc2 "6 3 2" with the chain along a gave 8367 fc3 and 15 fc2
  // supercells; in CAPS's cell the chain is c, so the same supercells are 2 2 3 and 2 3 6
  const auto ops = crystal_operations(pe());
  EXPECT_EQ(ops.size(), 8u);   // Pnam
  const auto c = count_displacements(pe(), {2, 2, 3}, {2, 3, 6});
  EXPECT_EQ(c.independent_atoms, 3);
  EXPECT_EQ(c.fc2, 15);
  EXPECT_EQ(c.first, 15);
  EXPECT_EQ(c.fc3(), 8367);
}

TEST(Phono3py, PlanAndScripts) {
  const auto dir = tmpdir("caps_p3_pe");
  Phono3pyOptions o;
  o.label = "pe";
  o.data_dir = kData;
  o.mode = "both";
  const auto r = write_phono3py_set(pe(), dir.string(), o);
  ASSERT_EQ(r.plans.size(), 2u);
  const auto& full = r.plans[0];
  const auto& rnd = r.plans[1];
  EXPECT_EQ(full.dim_fc3, (std::array<int, 3>{2, 2, 4}));
  EXPECT_EQ(full.dim_fc2, (std::array<int, 3>{2, 3, 6}));
  EXPECT_EQ(full.mesh, (std::array<int, 3>{10, 15, 160}));   // the guide's 160 15 10 in CAPS's axes
  EXPECT_EQ(full.chain_axis, 2);
  EXPECT_EQ(full.fc3_supercells, full.full.fc3());
  EXPECT_EQ(rnd.fc3_supercells, (full.full.fc3() + 9) / 10);
  EXPECT_LT(rnd.core_h_total, 0.2 * full.core_h_total);
  // cost model: 288 valence electrons (a 144-atom PE supercell) = 4.33 core-h
  EXPECT_NEAR(single_point_core_hours(288), 4.33, 1e-12);
  for (const char* m : {"full", "random"}) {
    const fs::path d = dir / m;
    for (const char* f : {"run_step1.sh", "run_step2.sh", "run_step3.sh", "run_step4.sh", "run_step5.sh", "run_step6.sh", "run_step7.sh", "run_dryrun.sh",
                          "cleanup.sh", "vasp_task.sh", "phono3py_env.sh", "make_potcar.sh"})
      EXPECT_TRUE(bash_ok(d / f)) << m << " " << f;
    const std::string inc = slurp(d / "INCAR");
    EXPECT_EQ(tag(inc, "IBRION"), "-1");
    EXPECT_EQ(tag(inc, "ISIF"), "0");
    EXPECT_EQ(tag(inc, "EDIFF"), "1E-8");
    EXPECT_EQ(tag(inc, "PREC"), "Accurate");
    EXPECT_EQ(tag(inc, "LREAL"), ".FALSE.");
    EXPECT_EQ(tag(inc, "ENCUT"), "800");
    const std::string s1 = slurp(d / "run_step1.sh");
    EXPECT_NE(s1.find("$INIT -d --dim 2 2 4 --dim-fc2 2 3 6 --pa auto -c POSCAR-unitcell"), std::string::npos) << s1;
    EXPECT_EQ(s1.find("--rd") != std::string::npos, std::string(m) == "random");
    EXPECT_EQ(slurp(d / "run_step6.sh").find("--symfc") != std::string::npos, std::string(m) == "random");
    EXPECT_NE(slurp(d / "run_step7.sh").find("--mesh 10 15 160 --br --ts 100 200 300 400"), std::string::npos);
  }
}

TEST(Phono3py, WholeChainWithFakes) {
  // fake phono3py (version 4 names) and fake VASP: run_dryrun → step2/3 tasks → step4 → step5 → step6 → step7
  if (std::system("command -v bash > /dev/null 2>&1") != 0) GTEST_SKIP() << "no bash";
  const auto root = tmpdir("caps_p3_run");
  const auto d = root / "set";
  Phono3pyOptions o;
  o.label = "pe";
  o.data_dir = kData;
  o.mode = "random";
  o.rd = 4;
  o.profile = "generic-slurm";
  write_phono3py_set(pe(), d.string(), o);
  spit(d / "POTCAR", "fake\n");
  const auto bin = root / "bin";
  spit(bin / "phono3py-init",
       "#!/bin/bash\n"
       "case \" $* \" in\n"
       "  *\" -d \"*) for i in 00001 00002 00003 00004; do cp POSCAR-unitcell POSCAR-$i; done; for i in 00001 00002; do cp POSCAR-unitcell POSCAR_FC2-$i; done; echo disp > phono3py_disp.yaml;;\n"
       "  *--cf3*) echo forces > FORCES_FC3;;\n"
       "  *--cf2*) echo forces > FORCES_FC2;;\n"
       "esac\n");
  spit(bin / "phono3py",
       "#!/bin/bash\n"
       "case \" $* \" in *--mesh*) echo kappa > kappa-m1015160.hdf5;; *) echo fc > fc3.hdf5; echo fc > fc2.hdf5;; esac\n");
  spit(bin / "srun",
       "#!/bin/bash\n"
       "printf '<?xml version=\"1.0\"?>\\n<modeling>\\n<varray name=\"forces\" >\\n%s\\n</varray>\\n</modeling>\\n' \"$(head -c 30000 /dev/zero | tr '\\0' x)\" > vasprun.xml\n"
       "echo ' EDIFF is reached' > OUTCAR; echo 'F=' > OSZICAR; echo ok\n");
  spit(bin / "python3", "#!/bin/bash\nexit 0\n");
  for (const char* f : {"phono3py-init", "phono3py", "srun", "python3"}) fs::permissions(bin / f, fs::perms::owner_all);
  const std::string env = "cd " + d.string() + " && PATH=" + bin.string() + ":$PATH HOME=" + root.string() + " SLURM_SUBMIT_DIR=" + d.string() + " SLURM_JOB_ID=5 ";
  ASSERT_EQ(std::system((env + "bash run_dryrun.sh > dry.log 2>&1").c_str()), 0) << slurp(d / "dry.log");
  EXPECT_NE(slurp(d / "run_step2.sh").find("#SBATCH --array=1-4\n"), std::string::npos);   // the real counts
  EXPECT_NE(slurp(d / "run_step3.sh").find("#SBATCH --array=1-2\n"), std::string::npos);
  EXPECT_NE(slurp(d / "dryrun_report.txt").find("projected"), std::string::npos);
  EXPECT_TRUE(fs::exists(d / "disp-00001" / "DONE"));
  // the check before the arrays ran: 2 of 4 fc3 done
  Json c = phono3py_check(d.string());
  EXPECT_EQ(int(c["fc3"].num("done", -1)), 2);
  EXPECT_EQ(int(c["fc3"].num("not_run", -1)), 2);
  EXPECT_FALSE(c["ready_for_step5"].boolean());
  for (int i = 1; i <= 4; ++i) ASSERT_EQ(std::system((env + "SLURM_ARRAY_TASK_ID=" + std::to_string(i) + " bash run_step2.sh >> arr.log 2>&1").c_str()), 0) << slurp(d / "arr.log");
  for (int i = 1; i <= 2; ++i) ASSERT_EQ(std::system((env + "SLURM_ARRAY_TASK_ID=" + std::to_string(i) + " bash run_step3.sh >> arr.log 2>&1").c_str()), 0);
  // a broken vasprun is caught and listed for a rerun
  spit(d / "disp-00003" / "vasprun.xml", "<modeling>\n");
  EXPECT_NE(std::system((env + "bash run_step4.sh > s4.log 2>&1").c_str()), 0);
  EXPECT_NE(slurp(d / "s4.log").find("sbatch --array=3 run_step2.sh"), std::string::npos) << slurp(d / "s4.log");
  ASSERT_EQ(std::system((env + "SLURM_ARRAY_TASK_ID=3 bash run_step2.sh >> arr.log 2>&1").c_str()), 0);
  c = phono3py_check(d.string());
  EXPECT_TRUE(c["ready_for_step5"].boolean()) << c.dump(1);
  for (const char* s : {"run_step5.sh", "run_step6.sh", "run_step7.sh"}) ASSERT_EQ(std::system((env + "bash " + s + " > " + s + ".log 2>&1").c_str()), 0) << s << slurp(d / (std::string(s) + ".log"));
  EXPECT_TRUE(fs::exists(d / "FORCES_FC3") && fs::exists(d / "fc3.hdf5") && fs::exists(d / "kappa-m1015160.hdf5"));
  // quota cleanup dry run lists, --yes deletes
  spit(d / "disp-00001" / "WAVECAR", std::string(1000, 'w'));
  ASSERT_EQ(std::system((env + "bash cleanup.sh > cl.log 2>&1").c_str()), 0);
  EXPECT_TRUE(fs::exists(d / "disp-00001" / "WAVECAR"));
  ASSERT_EQ(std::system((env + "bash cleanup.sh --yes > cl.log 2>&1").c_str()), 0);
  EXPECT_FALSE(fs::exists(d / "disp-00001" / "WAVECAR"));
}

TEST(Phono3py, TwoDimensionalPlan) {
  const System s = read_vasp_poscar(kSheet);
  Phono3pyOptions o;
  o.label = "Ti3C2O2";
  o.data_dir = kData;
  EXPECT_THROW(plan_phono3py(s, o, "full"), std::invalid_argument);   // a vacuum gap and no choice
  o.dimensionality = "2d";
  o.thickness = 9.8;
  o.dim_fc3 = {3, 3, 2};
  EXPECT_THROW(plan_phono3py(s, o, "full"), std::invalid_argument);   // 2 along the vacuum
  o.dim_fc3 = {0, 0, 0};
  const auto p = plan_phono3py(s, o, "full");
  EXPECT_EQ(p.dimensionality, "2d");
  EXPECT_EQ(p.dim_fc3[2], 1);
  EXPECT_EQ(p.dim_fc2[2], 1);
  EXPECT_EQ(p.kpoints_fc3[2], 1);
  EXPECT_EQ(p.mesh[2], 1);
  EXPECT_GT(p.dim_fc3[0], 2);
  EXPECT_NEAR(p.cell_height, norm(s.cell.c), 1e-9);
}

TEST(Phono3py, KappaAnalysisOnSyntheticExtract) {
  // synthetic extract: two temperatures, four modes; kappa = sum of the mode contributions
  const auto ex = tmpdir("caps_kappa_ex");
  spit(ex / "kappa_summary.json",
       "{\"temperatures\": [100, 300], \"kappa\": [[30, 2, 60, 0, 0, 0], [10, 1, 20, 0, 0, 0]], \"modes\": [\"modes_T100.csv\", \"modes_T300.csv\"]}");
  const std::string hdr = "frequency_THz,weight,k_xx,k_yy,k_zz,k_yz,k_xz,k_xy,gamma_THz,v_x,v_y,v_z\n";
  spit(ex / "modes_T100.csv", hdr + "1,1,15,1,30,0,0,0,0.01,10,0,40\n2,1,15,1,30,0,0,0,0.02,5,0,20\n");
  spit(ex / "modes_T300.csv", hdr + "1,1,4,0.5,8,0,0,0,0.0795774715,10,0,40\n2,1,3,0.25,6,0,0,0,0.2,5,0,20\n3,1,2,0.25,4,0,0,0,0.5,2,0,5\n4,1,1,0,2,0,0,0,0,0,0,0\n");
  const auto out = ex / "analysis";
  const Json r = analyze_kappa_extract(ex.string(), out.string(), 2);
  const auto& t = r["temperatures"].items();
  ASSERT_EQ(t.size(), 2u);
  EXPECT_DOUBLE_EQ(t[1].num("along_chain", 0), 20);
  EXPECT_DOUBLE_EQ(t[1].num("transverse_mean", 0), 5.5);
  EXPECT_DOUBLE_EQ(t[1]["sum_of_modes"].items()[2].number(), 20);
  // half of kappa_zz (10) is reached at the second mode by frequency (8 + 6)
  EXPECT_DOUBLE_EQ(t[1]["half_kappa_below_frequency_THz"].items()[2].number(), 2);
  // tau = 1 / (4 pi gamma): gamma 0.0795774715 THz → 1.0 ps; |v| = sqrt(10² + 40²) THz·Å → mfp = 41.23 Å = 4.123 nm
  const std::string modes = slurp(out / "modes_T300.csv");
  EXPECT_NE(modes.find("1,1,4.12311"), std::string::npos) << modes;
  EXPECT_TRUE(fs::exists(out / "cumulative_mfp_T300.csv"));
  EXPECT_TRUE(fs::exists(out / "kappa_tensor.csv"));
  // a 2D sheet: phono3py's kappa per cell volume → sheet conductance kappa × h, and kappa × h / t with a thickness
  const Json s2 = analyze_kappa_extract(ex.string(), (ex / "a2").string(), -1, "2d", 20.0, 5.0);
  const auto& r2 = s2["temperatures"].items()[1];
  EXPECT_DOUBLE_EQ(r2["kappa_sheet_xx_yy_xy"].items()[0].number(), 40);           // 10 × 20 / 5
  EXPECT_NEAR(r2["sheet_conductance_nW_per_K"].items()[0].number(), 20, 1e-9);    // 10 W/m/K × 20 Å = 2e-8 W/K = 20 nW/K
}

TEST(Phono3py, CommandLayerAndRecipe) {
  // the same commands the CLI, the C API, Python and the Studio call; a dft recipe runs them in order
  const auto dir = tmpdir("caps_p3_cmd");
  Json a = Json::object();
  Json in = Json::array();
  in.push_back("estimate"), in.push_back("pe-orthorhombic");
  a["inputs"] = in;
  const Json est = dft_run("phono3py", a, kData);
  EXPECT_EQ(est["plans"].items().size(), 2u);   // auto: the full and the random plan side by side
  EXPECT_NE(est.text("text").find("15 fc2 supercells"), std::string::npos);
  EXPECT_FALSE(fs::exists(dir / "x"));
  // a sheet without --dimensionality is refused, with it written
  Json b = Json::object();
  Json bi = Json::array();
  bi.push_back(kSheet);
  b["inputs"] = bi, b["out"] = (dir / "sheet").string();
  EXPECT_THROW(dft_run("vasp-relax", b, kData), std::invalid_argument);
  b["dimensionality"] = "2d";
  EXPECT_EQ(dft_run("vasp-relax", b, kData)["chain"].text("dimensionality"), "2d");
  EXPECT_NE(dft_help("phono3py").find("--dimensionality"), std::string::npos);
  EXPECT_NE(dft_help("vasp-relax").find("--encuts"), std::string::npos);
  // recipe: relax chain, then the phono3py set from the relaxation folder (its relax4 ENCUT)
  RecipeOptions ro;
  ro.out_dir = dir.string();
  ro.forcefield_dir = kData + "/forcefields";
  const Json recipe = Json::parse(R"({"dft": [
    {"command": "vasp-relax", "inputs": ["pe-orthorhombic"], "out": "pe_relax", "dimensionality": "3d", "encuts": "900,1350,900,900"},
    {"command": "phono3py", "inputs": ["setup", "pe_relax"], "mode": "random", "out": "pe_kappa"}]})");
  run_recipe(recipe, ro);
  EXPECT_TRUE(fs::exists(dir / "pe_relax" / "relax4" / "INCAR"));
  EXPECT_EQ(tag(slurp(dir / "pe_kappa" / "INCAR"), "ENCUT"), "900");
  EXPECT_NE(slurp(dir / "pe_kappa" / "run_step1.sh").find("--rd"), std::string::npos);
}

TEST(Phono3py, SetsCiteTheirMethods) {
  // phono3py and VASP run on the cluster, not in CAPS; the folders CAPS writes and their methods text cite them
  const auto dir = tmpdir("caps_p3_cite");
  auto run = [&](const char* cmd, std::vector<std::string> inputs, const std::string& out, std::vector<std::pair<const char*, const char*>> extra) {
    Json a = Json::object();
    Json in = Json::array();
    for (const auto& i : inputs) in.push_back(i);
    a["inputs"] = in, a["out"] = (dir / out).string();
    for (const auto& [k, v] : extra) a[k] = v;
    return dft_run(cmd, a, kData);
  };
  run("vasp-relax", {"pe-orthorhombic"}, "relax", {{"dimensionality", "3d"}});
  const std::string rb = slurp(dir / "relax" / "references.bib");
  for (const char* k : {"kresse1996", "blochl1994", "kresse1999", "perdew1996", "grimme2010", "grimme2011"}) EXPECT_NE(rb.find(k), std::string::npos) << k;
  std::vector<std::string> refs;
  EXPECT_NE(methods_text(*read_manifest((dir / "relax" / "POSCAR").string()), &refs).find("relaxed with VASP"), std::string::npos);
  EXPECT_EQ(refs.size(), 6u);
  run("phono3py", {"setup", "pe-orthorhombic"}, "random", {{"mode", "random"}});
  const std::string pr = slurp(dir / "random" / "references.bib");
  EXPECT_NE(pr.find("togo2015"), std::string::npos);
  EXPECT_NE(pr.find("togo2023"), std::string::npos);
  EXPECT_NE(pr.find("seko2024"), std::string::npos);   // symfc fits the random displacements
  EXPECT_EQ(pr.find("chaput2013"), std::string::npos);
  EXPECT_NE(slurp(dir / "random" / "README.txt").find("Phys. Rev. B 91, 094306 (2015)"), std::string::npos);
  run("phono3py", {"setup", "pe-orthorhombic"}, "full", {{"mode", "full"}, {"method", "lbte"}});
  const std::string fb = slurp(dir / "full" / "references.bib");
  EXPECT_EQ(fb.find("seko2024"), std::string::npos);
  EXPECT_NE(fb.find("chaput2013"), std::string::npos);   // the direct solution of the LBTE
  const std::string mt = methods_text(*read_manifest((dir / "full" / "POSCAR-unitcell").string()), &refs);
  EXPECT_NE(mt.find("with phono3py"), std::string::npos) << mt;
  EXPECT_NE(mt.find("solved directly"), std::string::npos) << mt;
  EXPECT_EQ(refs.size(), 9u);
}

TEST(Phono3py, AnalyzeKappaPyOnSyntheticHdf5) {
  // analyze_kappa.py (written by setup) on a synthetic kappa-m*.hdf5 with phono3py's dataset names and shapes; needs a
  // Python with h5py and numpy: CAPS_PYTHON (skipped without it)
  const char* py = std::getenv("CAPS_PYTHON");
  if (!py || std::system(("'" + std::string(py) + "' -c 'import h5py, numpy' > /dev/null 2>&1").c_str()) != 0) GTEST_SKIP() << "CAPS_PYTHON with h5py not set";
  const auto dir = tmpdir("caps_p3_h5");
  Phono3pyOptions o;
  o.label = "pe";
  o.data_dir = kData;
  o.mode = "full";
  write_phono3py_set(pe(), (dir / "set").string(), o);
  // 3 temperatures, 2 irreducible q-points (weights 1 and 3), 2 bands; mode_kappa sums / weight sum = kappa
  spit(dir / "make.py",
       "import h5py, numpy as np, sys\n"
       "T = np.array([100., 200., 300.]); w = np.array([1, 3]); f = np.array([[1.0, 2.0], [3.0, 4.0]])\n"
       "mk = np.zeros((3, 2, 2, 6)); mk[:, :, :, 2] = np.arange(1, 5).reshape(2, 2) * 4.0; mk[:, :, :, 0] = 4.0\n"
       "mk = mk / T[:, None, None, None] * 300\n"
       "kappa = mk.sum(axis=(1, 2)) / w.sum()\n"
       "g = np.full((3, 2, 2), 0.1); gv = np.zeros((2, 2, 3)); gv[:, :, 2] = 50.0\n"
       "with h5py.File(sys.argv[1], 'w') as h:\n"
       "    h['temperature'] = T; h['kappa'] = kappa; h['mode_kappa'] = mk; h['frequency'] = f; h['weight'] = w\n"
       "    h['gamma'] = g; h['group_velocity'] = gv; h['mesh'] = np.array([10, 15, 160]); h['heat_capacity'] = np.ones((3, 2, 2))\n");
  const std::string h5 = (dir / "set" / "kappa-m1015160.hdf5").string();
  ASSERT_EQ(std::system(("'" + std::string(py) + "' " + (dir / "make.py").string() + " " + h5).c_str()), 0);
  Json a = Json::object();
  Json in = Json::array();
  in.push_back("analyze"), in.push_back((dir / "set").string());
  a["inputs"] = in;
  a["python"] = py;
  const Json r = dft_run("phono3py", a, kData);
  const auto& rows = r["analysis"]["temperatures"].items();
  ASSERT_EQ(rows.size(), 3u);
  // at 300 K: zz = (4+8+12+16)/4 = 10, xx = 16/4 = 4; along the chain (c) 10, across (4 + 0) / 2 = 2
  EXPECT_NEAR(rows[2].num("along_chain", 0), 10, 1e-9);
  EXPECT_NEAR(rows[2].num("transverse_mean", 0), 2, 1e-9);
  EXPECT_NEAR(rows[2]["sum_of_modes"].items()[2].number(), 10, 1e-9);
  // tau = 1 / (4 pi 0.1) = 0.7958 ps, |v| 50 THz·Å = 5 km/s
  EXPECT_NEAR(rows[2].num("kappa_weighted_lifetime_ps", 0), 1 / (4 * M_PI * 0.1), 1e-6);
  EXPECT_NEAR(rows[2].num("max_group_velocity_km_s", 0), 5, 1e-9);
}
