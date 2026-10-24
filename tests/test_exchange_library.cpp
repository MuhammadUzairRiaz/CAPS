#include <gtest/gtest.h>

#include <cmath>
#include <fstream>
#include <set>
#include <sstream>
#include <string>

#include "caps/exchange_library.hpp"
#include "caps/json.hpp"

using namespace caps;

namespace {

const std::vector<ExchangeChemistry>& library() {
  static const auto lib = load_exchange_library(std::string(CAPS_SOURCE_DIR) + "/data/exchange/library.json");
  return lib;
}

}  // namespace

// Every family of the study's table is there; every scheme named is in the reaction library; every barrier carries a
// DOI, a status and, for one number in kJ/mol, the same value in kcal/mol.
TEST(ExchangeLibrary, ChemistriesSchemesAndSources) {
  const auto& lib = library();
  std::set<std::string> ids;
  for (const auto& c : lib) EXPECT_TRUE(ids.insert(c.id).second) << "duplicate id " << c.id;
  for (const char* id : {"transesterification", "transcarbamoylation", "vinylogous_urethane", "vinylogous_urea", "transimination", "imine_metathesis",
                         "disulfide_aromatic", "thiol_disulfide", "polysulfide_rubber", "sulfenamide", "diselenide", "boronic_ester", "dioxaborolane",
                         "siloxane", "sulfonium_transalkylation", "enamine_isocyanide", "acetal", "diels_alder", "hindered_urea", "thiol_michael",
                         "hydrogen_bond", "bead_spring_swap", "three_body_swap"})
    EXPECT_TRUE(ids.count(id)) << id;
  std::ifstream in(std::string(CAPS_SOURCE_DIR) + "/data/reactions/library.json");
  std::stringstream ss;
  ss << in.rdbuf();
  std::set<std::string> schemes;
  const Json reactions = Json::parse(ss.str());
  for (const auto& r : reactions["reactions"].items()) schemes.insert(r.text("id"));
  for (const auto& c : lib) {
    if (!c.scheme.empty()) EXPECT_TRUE(schemes.count(c.scheme)) << c.id << ": scheme " << c.scheme;
    EXPECT_TRUE(c.move == "swap" || c.move == "shift" || c.move == "bind" || c.move == "none") << c.id;
    for (const auto& b : c.barriers) {
      EXPECT_EQ(b.source.rfind("10.", 0), 0u) << c.id << ": source " << b.source;
      EXPECT_TRUE(b.status == "read" || b.status == "quoted") << c.id;
      EXPECT_FALSE(b.kind.empty()) << c.id;
      if (b.unit == "kJ/mol" && std::isfinite(b.kcal_per_mol)) EXPECT_NEAR(b.kcal_per_mol, std::stod(b.value) / 4.184, 1e-3) << c.id;
      if (b.unit == "kcal/mol" && std::isfinite(b.kcal_per_mol)) EXPECT_NEAR(b.kcal_per_mol, std::stod(b.value), 1e-9) << c.id;
    }
    if (std::isfinite(c.default_kcal_per_mol)) EXPECT_FALSE(c.default_why.empty()) << c.id;
    if (c.barriers.empty()) EXPECT_FALSE(c.notes.empty() && c.to_read.empty()) << c.id << ": no barrier and no word on what to read";
  }
}

// A chemistry gives the exchange engine its move, barrier and capture distance.
TEST(ExchangeLibrary, OptionsForTheEngine) {
  const auto& lib = library();
  const auto te = exchange_options_for(find_exchange_chemistry(lib, "transesterification"));
  EXPECT_EQ(te.kind, ExchangeKind::Shift);
  EXPECT_NEAR(te.barrier, 14.6, 1e-12);
  EXPECT_NEAR(te.capture, 4.25, 1e-12);
  const auto ds = exchange_options_for(find_exchange_chemistry(lib, "disulfide_aromatic"));
  EXPECT_EQ(ds.kind, ExchangeKind::Swap);
  EXPECT_NEAR(ds.barrier, 12.02, 1e-12);
  EXPECT_NEAR(ds.capture, 4.12, 1e-12);
  const auto da = exchange_options_for(find_exchange_chemistry(lib, "diels_alder"));
  EXPECT_EQ(da.kind, ExchangeKind::Bind);
  EXPECT_EQ(da.barrier, 0.0);   // no barrier read yet
  EXPECT_THROW(exchange_options_for(find_exchange_chemistry(lib, "hydrogen_bond")), std::invalid_argument);
  EXPECT_THROW(find_exchange_chemistry(lib, "no_such_chemistry"), std::invalid_argument);
}

// Yang et al. 2016: with E_a = 68.18 kJ/mol the rate at 420 K is 27 % of the rate at 450 K.
TEST(ExchangeLibrary, ArrheniusProbability) {
  EXPECT_NEAR(arrhenius_probability(68.18 / 4.184, 420, 450), 0.27, 0.005);
  EXPECT_DOUBLE_EQ(arrhenius_probability(10, 500, 450), 1.0);   // capped above the reference
  EXPECT_THROW(arrhenius_probability(10, 0, 450), std::invalid_argument);
}
