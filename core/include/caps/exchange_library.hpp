// CAPS Exchange, chemistry library (stage E4 of docs/studies/vitrimers/caps_module_proposal.md): the exchange
// chemistries of data/exchange/library.json (format caps-exchange), each with its mechanism, the move the exchange engine
// makes, its reaction scheme (data/reactions/library.json) and every barrier read in the literature with its source.
#pragma once
#include <limits>
#include <string>
#include <vector>

#include "caps/exchange_moves.hpp"

namespace caps {

struct ExchangeBarrier {
  std::string value;            // as written ("12.02", "150 to 180")
  double kcal_per_mol = std::numeric_limits<double>::quiet_NaN();   // when the value is one number in an energy unit
  std::string unit, kind, method, conditions, source, status, note;  // source: a DOI; status: read | quoted
};

struct ExchangeChemistry {
  std::string id, name, family, mechanism, move, scheme;   // move: swap | shift | bind | none; scheme: "" when none yet
  std::string template_name;    // the built-in React template a network runs in CAPS's engine (its scheme is a model compound's)
  std::vector<ExchangeBarrier> barriers;
  double default_kcal_per_mol = std::numeric_limits<double>::quiet_NaN();   // the barrier a run uses (NaN: none)
  std::string default_why;
  double capture_angstrom = std::numeric_limits<double>::quiet_NaN();       // all-atom capture distance (NaN: none)
  std::string capture_pair, capture_source;
  std::vector<std::string> systems, to_read;
  std::string notes;
};

// Reads the library (throws std::runtime_error when the file cannot be read or is not caps-exchange).
std::vector<ExchangeChemistry> load_exchange_library(const std::string& path);
// The chemistry with this id (throws std::invalid_argument naming the ids there are).
const ExchangeChemistry& find_exchange_chemistry(const std::vector<ExchangeChemistry>& lib, const std::string& id);
// Exchange engine settings from a chemistry: its move, its default barrier (0 when it has none) and its capture
// distance (the engine's default when it has none). Groups and labels are the caller's. Throws for a chemistry whose
// exchange is not a move (physical bonds, the three-body potential).
ExchangeOptions exchange_options_for(const ExchangeChemistry& c);
// The probability of an exchange attempt at T relative to T_ref for an Arrhenius rate, exp(−E_a/R (1/T − 1/T_ref)),
// capped at 1 (the acceptance Yang et al. 2016 normalise to 1 at 450 K; a LAMMPS fix bond/react prob). E_a in kcal/mol.
double arrhenius_probability(double ea_kcal_per_mol, double temperature, double t_ref);

}  // namespace caps
