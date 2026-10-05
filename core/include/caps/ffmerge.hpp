// CAPS force fields by group: parts of one structure typed and parameterised by different force fields (a crystal or
// filler with one, the polymer with another; each component of a blend with its own), merged into one ForceField that
// the engine, LAMMPS, GROMACS and DL_POLY exports use.
//
// Within each part every term is that force field's own, its unlike Lennard-Jones pairs by its own mixing rule. Between
// parts the pairs follow the cross rule chosen here — ε geometric or arithmetic, σ arithmetic or geometric, or the
// sixth-power rule for both (as DL_FIELD's "multiple potentials" asks for the rule between different force fields,
// with no default) — or explicit ε, σ given for a pair of types. Every cross pair is written out explicitly, so the
// engine files carry exactly what CAPS computes.
//
// What one simulation cannot hold is refused with the reason, never approximated silently:
//   · parts with different 1-4 scaling (GAFF 0.5 / 0.8333 against OPLS-AA 0.5 / 0.5) keep their own (scaling14 "own",
//     the default): each 1-4 pair scaled by its part's force field — LAMMPS a pair sub-style per part with its own
//     special weights, GROMACS each 1-4 pair with its own fudge; "first" takes the first part's for all, "refuse" stops
//   · 9-6 (class II) and 12-6 parts: the cross pairs need one form; cross96 "rmin" gives the 9-6 site a 12-6 form with
//     the same well depth ε and minimum r_min; "area" keeps r_min and takes the 12-6 depth with the same ∫ U r² dr from
//     r_min to the cut-off (both said in the notes)
//   · coarse-grained settings (dielectric, reaction field, force switches), Stillinger–Weber in more than one part,
//     DREIDING hydrogen bonds, CHARMM 1-4 types, virtual sites across parts
#pragma once
#include <cstdint>

#include <string>
#include <vector>

#include "caps/field.hpp"

namespace caps {

struct FFPart {
  const ForceField* ff = nullptr;       // parameterised for the part's atoms alone, in their order
  std::vector<uint32_t> atoms;          // the part's atoms in the whole structure (ff's atom k is atoms[k])
  std::string tag;                      // a short name ("filler", "PS"): added to type names that collide
};

struct CrossPair { std::string a, b; double eps = 0, sigma = 0; };   // merged type names (kcal/mol, Å)

struct MergeOptions {
  // auto: as the parts' force fields mix when they all mix alike (OPLS geometric ε and σ; AMBER/GAFF/CHARMM Lorentz–Berthelot;
  // class II sixth power); parts that mix differently: sixth power with a class II part, else Lorentz–Berthelot
  std::string eps_rule = "auto";     // auto | geometric | arithmetic
  std::string sigma_rule = "auto";   // auto | arithmetic | geometric | sixthpower (ε then by the sixth-power rule too)
  std::string scaling14 = "own";        // own (each part's 1-4 pairs by its own scaling) | first | refuse
  std::string cross96 = "refuse";       // refuse | rmin | area
  double refit_cutoff = 12.0;           // Å, cross96 "area": the upper limit of the matched integral (the run's cut-off)
  std::vector<CrossPair> explicit_pairs;
};

// cross96 "area": the 12-6 well depth over the 9-6 one for a site of minimum r0 (Å), equal ∫ U r² dr from r0 to rc
double lj96_to_126_depth_ratio(double r0, double rc);

// Throws FieldError with the reason when the parts cannot share one simulation (above) or an atom is in no part / two.
ForceField merge_forcefields(size_t natoms, const std::vector<FFPart>& parts, const MergeOptions& o, std::vector<std::string>* notes = nullptr);

}  // namespace caps
