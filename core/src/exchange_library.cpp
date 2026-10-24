// CAPS Exchange, chemistry library (see exchange_library.hpp).
#include "caps/exchange_library.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>

#include "caps/json.hpp"

namespace caps {

namespace {

constexpr double kGasKcal = 1.98720425864083e-3;   // kcal/(mol K)

std::string text_of(const Json& j, const char* k) {
  if (!j.has(k) || j[k].is_null()) return "";
  return j[k].is_string() ? j[k].str() : j[k].dump(0);
}

}  // namespace

std::vector<ExchangeChemistry> load_exchange_library(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("cannot read the exchange library " + path);
  std::stringstream ss;
  ss << in.rdbuf();
  const Json j = Json::parse(ss.str());
  if (j.text("format") != "caps-exchange") throw std::runtime_error(path + " is not a caps-exchange library");
  std::vector<ExchangeChemistry> lib;
  for (const auto& e : j["chemistries"].items()) {
    ExchangeChemistry c;
    c.id = e.text("id"), c.name = e.text("name"), c.family = e.text("family"), c.mechanism = e.text("mechanism");
    c.move = e.text("move"), c.scheme = text_of(e, "scheme"), c.notes = e.text("notes");
    for (const auto& b : e["barriers"].items()) {
      ExchangeBarrier x;
      x.value = text_of(b, "value");
      if (b.has("kcal_per_mol")) x.kcal_per_mol = b["kcal_per_mol"].number();
      x.unit = b.text("unit"), x.kind = b.text("kind"), x.method = b.text("method"), x.conditions = b.text("conditions");
      x.source = b.text("source"), x.status = b.text("status"), x.note = b.text("note");
      c.barriers.push_back(x);
    }
    if (e.has("default") && e["default"].is_object()) {
      const Json& d = e["default"];
      if (d.has("barrier_kcal_per_mol")) c.default_kcal_per_mol = d["barrier_kcal_per_mol"].number();
      c.default_why = d.text("why");
    }
    if (e.has("capture") && e["capture"].is_object() && e["capture"].text("unit") == "Å") {
      c.capture_angstrom = e["capture"]["value"].number();
      c.capture_pair = e["capture"].text("pair"), c.capture_source = e["capture"].text("source");
    }
    for (const auto& s : e["systems"].items()) c.systems.push_back(s.str());
    for (const auto& s : e["to_read"].items()) c.to_read.push_back(s.str());
    lib.push_back(std::move(c));
  }
  return lib;
}

const ExchangeChemistry& find_exchange_chemistry(const std::vector<ExchangeChemistry>& lib, const std::string& id) {
  std::string ids;
  for (const auto& c : lib) {
    if (c.id == id) return c;
    ids += (ids.empty() ? "" : ", ") + c.id;
  }
  throw std::invalid_argument("no exchange chemistry '" + id + "' (there are: " + ids + ")");
}

ExchangeOptions exchange_options_for(const ExchangeChemistry& c) {
  if (c.move != "swap" && c.move != "shift" && c.move != "bind")
    throw std::invalid_argument(c.id + ": its exchange is not a move (" + (c.move.empty() ? "none" : c.move) + ")");
  ExchangeOptions o;
  o.kind = exchange_kind_from_string(c.move);
  o.barrier = std::isfinite(c.default_kcal_per_mol) ? c.default_kcal_per_mol : 0.0;
  if (std::isfinite(c.capture_angstrom)) o.capture = c.capture_angstrom;
  return o;
}

double arrhenius_probability(double ea_kcal_per_mol, double temperature, double t_ref) {
  if (!(temperature > 0) || !(t_ref > 0)) throw std::invalid_argument("arrhenius_probability: temperatures must be positive");
  return std::min(1.0, std::exp(-ea_kcal_per_mol / kGasKcal * (1.0 / temperature - 1.0 / t_ref)));
}

}  // namespace caps
