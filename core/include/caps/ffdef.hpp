// CAPS force-field definitions: a data format for force fields (types, typing rules, parameters, functional forms),
// importers from other tools' files, and assignment of parameters to a structure.
//
// The format is JSON ("caps-forcefield", version 1). Parameter rules match atom-type names with glob patterns
// (* and ?). Within a section the LAST matching rule wins, as in moltemplate's "By Type" sections, so a converted
// moltemplate force field assigns exactly what moltemplate assigns; user overlays appended later win over the base.
#pragma once
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include "caps/field.hpp"
#include "caps/system.hpp"

namespace caps {

struct FFError : std::runtime_error {
  using std::runtime_error::runtime_error;
};

struct FFType {
  std::string name;
  int element = 0;
  double mass = 0;
  double charge = std::numeric_limits<double>::quiet_NaN();   // type charge (OPLS), NaN when the force field has none
  std::string description;
  std::string smarts;                  // typing rule (empty: cannot be assigned automatically)
  std::vector<std::string> overrides;  // types this rule takes precedence over when both match
  int priority = 0;                    // among matching rules not overridden, the highest priority wins
  std::string source;                  // literature reference or origin of the parameters
  // Equivalences: the type name used when looking up one kind of parameter ("vdw", "bond", "angle",
  // "dihedral", "improper", "increment"), as DL_FIELD's EQUIVALENCE and the .frc #equivalence tables.
  std::map<std::string, std::string> equiv;
  std::vector<std::string> aliases;    // other names for this type (DL_FIELD atom types, moltemplate short names)
};

// A typing rule: atoms matching the SMARTS (first pattern atom = the atom typed) get `type`. Among the rules that match
// an atom, rules whose type another match overrides drop out, then the highest priority wins, then the earliest rule.
struct TypingRule {
  std::string type, smarts, description;
  std::vector<std::string> overrides;
  int priority = 0;
};

// One parameter rule: glob patterns on type names (2 for bonds, 3 angles, 4 dihedrals / impropers), a style and its
// parameters in LAMMPS order and units (real).
//
// Class II rules (style class2) keep their cross terms as named groups, each in LAMMPS order and units:
//   angles     params θ0 K2 K3 K4;  bb: M r1 r2;  ba: N1 N2 r1 r2
//   dihedrals  params K1 φ1 K2 φ2 K3 φ3;  mbt: A1 A2 A3 r2;  ebt: B1 B2 B3 C1 C2 C3 r1 r3;
//              at: D1 D2 D3 E1 E2 E3 θ1 θ2;  aat: M θ1 θ2;  bb13: N r1 r3
//   impropers  params K χ0;  aa: M1 M2 M3 θ1 θ2 θ3
// Cross terms are direction dependent: a rule that matches an interaction backwards assigns it with the atoms
// reversed, so r1 / B / D always belong to the end matched by the rule's first pattern.
struct FFRule {
  std::string name;
  std::vector<std::string> match;
  std::string style;
  std::vector<double> params;
  std::map<std::string, std::vector<double>> cross;
  std::string comment;
};

struct FFDef {
  std::string name, version, source;
  std::vector<std::string> references;
  std::string units = "real";
  // default styles (a rule may name its own style, as in LAMMPS hybrid styles)
  std::string pair_style = "lj/cut/coul/long", bond_style = "harmonic", angle_style = "harmonic", dihedral_style = "fourier",
              improper_style = "cvff";
  std::string mixing = "arithmetic";   // arithmetic (Lorentz–Berthelot), geometric (OPLS), sixthpower (class II)
  double special_lj[3] = {0, 0, 0.5};  // 1-2, 1-3, 1-4 scaling
  double special_coul[3] = {0, 0, 0.8333333333};
  double cutoff = 10.0;
  // How improper quadruples are formed and ordered (moltemplate symmetry plugins):
  //   "center3_sorted"  centre in position 3, the others sorted by atom index (AMBER / GAFF, gaff_imp.py)
  //   "center1_sorted"  centre in position 1, the others sorted by atom index (OPLS, cenIsortJKL.py)
  //   "center3_type_sorted"  centre third, the others by type name then bond order (AMBER tleap, DL_FIELD)
  //   "center2_sorted"  centre in position 2 (class II, cenJsortIKL); the outer atoms keep the order of the
  //                     matching rule, since the angle-angle terms depend on it
  std::string improper_order = "center3_sorted";
  // How equivalences are used: "replace" (DL_FIELD: look up by the equivalent name only) or "fallback"
  // (msi2lmp: the type's own name first, then the equivalent).
  std::string equivalence = "replace";
  // Impropers keep the outer atoms in the order the rule matched them (DL_FIELD, msi2lmp) rather than sorted.
  bool improper_matched_order = false;
  // Improper rules may also match their atom pattern read backwards (DL_FIELD's torsion-type impropers).
  bool improper_reversible = false;
  // Every explicit (wildcard-free) improper rule that matches a centre applies, not just the best one (CHARMM).
  bool improper_all_explicit = false;
  // Only centres with at most this many neighbours get impropers (3: planar centres); 0 = any.
  int improper_max_neighbours = 0;
  // Torsions whose rule has a wildcard end get K divided by (connections − 1) of the neighbouring central atom,
  // per wildcard end (msi2lmp and DL_FIELD for CVFF): "none" or "msi2lmp".
  std::string wildcard_torsion_scaling = "none";
  std::vector<FFType> types;
  std::vector<FFRule> pairs;           // match: [a] (self) or [a, b] (explicit pair); params: epsilon, sigma
  std::vector<FFRule> bonds, angles, dihedrals, impropers;
  // Bond increments (class II charges): match [a, b], params [δa, δb]; an atom's charge is the sum over its bonds
  // (plus its type charge, if any). The last matching rule wins, in either direction.
  std::vector<FFRule> bond_increments;
  // Automatic parameters (Accelrys "cff91_auto" / "cvff_auto"): used only when no rule above matches. They are
  // looked up by the types' auto equivalences, which depend on the position in the interaction:
  //   auto_bond; auto_angle_end / auto_angle_apex; auto_torsion_end / auto_torsion_center; auto_oop_end / _center.
  std::vector<FFRule> auto_bonds, auto_angles, auto_dihedrals, auto_impropers;
  // Class II cross terms given as their own rule sets (Accelrys .frc): looked up separately from the main term, each by
  // its own patterns, with the reference lengths / angles taken from the structure's assigned bonds and angles
  // (as msi2lmp / Discover do). Groups and parameters (LAMMPS order):
  //   bb: M (3 types) · ba: N1 N2 · mbt: A1 A2 A3 (4) · ebt: B1 B2 B3 C1 C2 C3 · at: D1 D2 D3 E1 E2 E3 · aat: M · bb13: N
  //   aa: K for the angle pair i-j-k / k-j-l (4 types, j the centre)
  std::map<std::string, std::vector<FFRule>> cross_rules;
  // Out-of-plane terms of class II force fields as msi2lmp builds them: "msi2lmp": every three-connected centre gets a
  // Wilson term and angle-angle terms, every centre with more neighbours angle-angle terms for each triple.
  std::string oop_scheme;
  // Typing rules (caps/typing.hpp). In the JSON file "typing" is either the rules or the name of a "caps-typing" file
  // next to the force field; types' own "smarts" become rules after them.
  std::vector<TypingRule> typing;
  std::vector<std::pair<std::string, std::string>> typing_pairs;   // conjugated pairs (GAFF cc/cd, ...), see typing.hpp
  bool typing_ordered = false;
  bool typing_unknown_untyped = false;   // rules may name types this file lacks: their atoms end up untyped
  bool typing_pairs_double_same = false;   // pairs keep one type across a double bond (CGenFF CG2DC1/2), not GAFF's   // rules are an ordered list (antechamber): the first match is intended, not ambiguous
  std::string typing_source;
  std::vector<std::string> notes;

  const FFType* type(const std::string& name) const;
};

FFDef load_forcefield(const std::string& path);
// Reads a "caps-typing" rules file and appends its rules to the force field (types must exist in it).
void load_typing(FFDef& ff, const std::string& path);
void save_forcefield(const FFDef& ff, const std::string& path);
// Appends the overlay's types and rules after the base's (so they win), replacing types of the same name.
void merge_forcefield(FFDef& base, const FFDef& overlay);

// moltemplate force-field files (.lt): "In Init" styles and special_bonds, "Data Masses", "In Charges",
// pair / bond / angle / dihedral / improper coefficients, the "By Type" rules and their symmetry plugins, "replace"
// aliases, and the type descriptions from the comments.
FFDef import_moltemplate(const std::string& path);

// DL_FIELD force-field libraries: NAME.par (parameters), NAME.sf (atom types; templates are not converted yet)
// and NAME.bci (bond charge increments), as shipped in DL_FIELD's lib/ directory. Parameters are converted to
// CAPS units (kcal/mol, Å, degrees) and LAMMPS functional forms as DL_FIELD itself writes them for LAMMPS.
// Families not yet interpreted throw FFError naming the family.
FFDef import_dlfield(const std::string& par, const std::string& sf = "", const std::string& bci = "");

bool glob_match(const std::string& pattern, const std::string& text);
std::string glob_escape(const std::string& name);   // a pattern that matches exactly this name

struct ParamReport {
  std::vector<std::string> missing;    // interactions with no matching rule (each once, with an example)
  std::map<std::string, int> used;     // rule name → interactions
  std::vector<std::string> notes;
  bool complete() const { return missing.empty(); }
};

// Build the evaluator force field for a structure whose atoms carry force-field type names (one per atom).
// Charges: `charges` = "types" (from the force field: type charges and / or bond increments; error if neither
// gives a charge), "keep" (the structure's) or
// "gasteiger" (C/H/N/O only). Bonded interactions are generated from the bonds as moltemplate does. Throws FFError
// listing every missing parameter unless `allow_missing` (then the report lists them and those terms are left out).
ForceField parameterize(const System& s, const FFDef& ff, const std::vector<std::string>& types, const std::string& charges,
                        ParamReport* report = nullptr, bool allow_missing = false);

}  // namespace caps
