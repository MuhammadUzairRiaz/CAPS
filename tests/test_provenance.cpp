#include <gtest/gtest.h>

#include <filesystem>
#include <fstream>

#include "caps/provenance.hpp"

using namespace caps;

namespace {

Manifest sample(uint64_t seed) {
  Manifest m;
  m.inputs = {{"PS_atactic_DP40.caps", "9f2c"}, {"gaff2.json", "71be"}};
  ProvStep g;
  g.engine = "grow.trials";
  g.params = {{"chains", "20"}, {"DP", "40"}};
  g.rng = "mt19937-64 · seed " + std::to_string(seed);
  g.cites = {"matsumoto1998"};
  ProvStep r;
  r.engine = "relax.lbfgs";
  r.params = {{"|F|max", "0.02 kcal/mol/Å"}};
  r.cites = {"liu1989", "auhl2003"};
  r.approximations = {{"van der Waals", "cut-off 12 Å"}, {"Electrostatics", "off"}};
  ProvStep md;
  md.engine = "dynamics.npt";
  md.cites = {"swope1982", "bussi2007", "bernetti2020", "essmann1995"};
  md.approximations = {{"van der Waals", "cut-off 12 Å + tail correction"}, {"Electrostatics", "SPME · relative tolerance 1e-05"}, {"Precision", "double · reproducible with 10 threads"}};
  m.steps = {g, r, md};
  return m;
}

}  // namespace

TEST(Provenance, JsonRoundTripAndSidecar) {
  const Manifest m = sample(20260923);
  const Manifest back = manifest_from_json(Json::parse(manifest_json(m).dump(2)));
  ASSERT_EQ(back.steps.size(), 3u);
  EXPECT_EQ(back.steps[0].rng, "mt19937-64 · seed 20260923");
  EXPECT_EQ(back.steps[1].cites, (std::vector<std::string>{"liu1989", "auhl2003"}));
  EXPECT_EQ(back.inputs, m.inputs);

  const auto data = (std::filesystem::temp_directory_path() / "caps_prov_cell.data").string();
  write_manifest(m, data);
  EXPECT_TRUE(std::filesystem::exists(data + ".provenance.json"));
  const auto r = read_manifest(data);
  ASSERT_TRUE(r.has_value());
  EXPECT_EQ(r->steps[2].engine, "dynamics.npt");
  EXPECT_FALSE(read_manifest(data + ".missing").has_value());
}

TEST(Provenance, ApproximationsTakeTheLatestStep) {
  const auto a = approximations(sample(1));
  ASSERT_GE(a.size(), 3u);
  EXPECT_EQ(a[0].first, "van der Waals");
  EXPECT_EQ(a[0].second, "cut-off 12 Å + tail correction");
  EXPECT_EQ(a[1].second, "SPME · relative tolerance 1e-05");
  EXPECT_EQ(a[2].first, "Precision");
}

TEST(Provenance, CompareFindsOnlyTheSeed) {
  const auto d = compare(sample(20260923), sample(20260924));
  ASSERT_EQ(d.rows.size(), 1u);
  EXPECT_EQ(d.rows[0].engine, "grow.trials");
  EXPECT_EQ(d.rows[0].key, "rng");
  EXPECT_EQ(d.differing_steps, 1);
  EXPECT_TRUE(d.same_inputs);
  EXPECT_TRUE(d.notes.empty());

  Manifest b = sample(20260923);
  b.steps.pop_back();
  b.inputs[1].second = "0000";
  const auto e = compare(sample(20260923), b);
  EXPECT_FALSE(e.same_inputs);
  EXPECT_EQ(e.notes.size(), 2u);   // inputs differ; dynamics only in the first
}

TEST(Provenance, BibtexForEveryCitedMethod) {
  const std::string bib = bibtex(all_cites(sample(1)));
  for (const char* k : {"matsumoto1998", "liu1989", "auhl2003", "swope1982", "bussi2007", "bernetti2020", "essmann1995"})
    EXPECT_NE(bib.find(std::string("{") + k + ","), std::string::npos) << k;
  EXPECT_NE(bib.find("doi = {10.1063/1.2408420}"), std::string::npos);
  // every key the C API cites is in the table
  for (const char* k : {"abascal2005", "auhl2003", "berendsen1984", "berendsen1987", "bernetti2020", "bitzek2006", "bussi2007", "cordero2008",
                        "engh1991", "essmann1995", "fennell2006", "gasteiger1980", "hall1981", "jorgensen1983", "larsen2011", "liu1989",
                        "martinez2009", "matsumoto1998", "parsons2005", "polak1969", "rappe1991", "rappe1992", "swope1982", "wang2004"})
    EXPECT_TRUE(known_citation(k)) << k;
}
