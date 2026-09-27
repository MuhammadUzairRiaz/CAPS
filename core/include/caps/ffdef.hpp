// CAPS force-field definitions: a data format for force fields (types, typing rules, parameters, functional forms),
// importers from other tools' files, and assignment of parameters to a structure.
//
// The format is JSON ("caps-forcefield", version 1). Parameter rules match atom-type names with glob patterns
// (* and ?). Within a section the LAST matching rule wins, as in moltemplate's "By Type" sections, so a converted
// moltemplate force field assigns exactly what moltemplate assigns; user overlays appended later win over the base.
#pragma once
#include <limits>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "caps/field.hpp"
#include "caps/resolution.hpp"
#include "caps/json.hpp"
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
  // Conditions on the whole structure (compound-specific potential sets of ionic solids: NaCl's Na is not NaF's Na):
  // every element of `requires` present, none of `excludes`. `atom_name`: only atoms of that name (the shells CAPS adds).
  std::vector<int> needs_elements, no_elements;
  std::string atom_name;
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
  double timestep = 0;   // fs: the model's usual MD time step for engine inputs (Martini 20); 0: 0.5 fs
  // Coarse-grained pair settings (JSON "pair_settings"): lj/gromacs inner radius, coul/gromacs inner radius and
  // relative permittivity (MARTINI 9 Å, 1e-6 Å, 15); model_cutoff: the cut-off belongs to the model (MARTINI 12 Å)
  double lj_inner = 0, coul_inner = 0, dielectric = 1;
  bool model_cutoff = false;
  // Martini 3 ("coulomb": "reaction-field", "eps_rf", "lj_modifier": "potential-shift" in pair_settings)
  bool coul_rf = false, lj_shift = false;
  // CHARMM (pair_settings "lj_modifier": "charmm-force-switch", "lj_inner" 10 with the model cut-off 12): LAMMPS lj/charmmfsw
  bool lj_fsw = false;
  double eps_rf = 0;
  // Torsions only where the file defines them ("torsion_terms": "if_defined"; MARTINI, SDK: TORSION IGNORE in their
  // sources): a dihedral with no term is not missing
  bool torsions_if_defined = false;
  bool angles_if_defined = false;   // "angle_terms": "if_defined" (MARTINI's ANGLE WARN): a missing angle is a note, not an error
  // How improper quadruples are formed and ordered (moltemplate symmetry plugins):
  //   "center3_sorted"  centre in position 3, the others sorted by atom index (AMBER / GAFF, gaff_imp.py)
  //   "center1_sorted"  centre in position 1, the others sorted by atom index (OPLS, cenIsortJKL.py)
  //   "center3_type_sorted"  centre third, the others by type name then bond order (AMBER tleap, DL_FIELD)
  //   "center2_sorted"  centre in position 2 (class II, cenJsortIKL); the outer atoms keep the order of the
  //                     matching rule, since the angle-angle terms depend on it
  std::string improper_order = "center3_sorted";
  // How a matched improper is written (and evaluated): "" as matched; "center2": the matched quartet i j C l written
  // i C j l (the DL-derived OPLS files: their rules put the centre third, their LAMMPS terms second; their DL_POLY
  // FIELD terms are j i C l, which write_dlpoly reproduces)
  std::string improper_written;
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
  // per wildcard end (msi2lmp): "msi2lmp"; "torsions": a rule with both ends wild over the torsions that exist about
  // the bond (a three-membered ring has fewer); "none".
  std::string wildcard_torsion_scaling = "none";
  // Charges "increments": the force field named here (a file beside this one, e.g. "opls2005.json") types the
  // structure and gives its bond-increment charges, used with this force field's types (see companion_charges).
  std::string charge_increments_from;
  // "dreiding1990": a torsion no rule lists takes DREIDING's own rule by the hybridisation of its central atoms (Mayo,
  // Olafson, Goddard 1990, cases a-j), its barrier divided over the torsions about the bond
  std::string torsion_rules;
  // ClayFF ("angle_contacts": 2.6): a bend i-j-k whose i has no bonds of its own (a metal of an octahedral sheet) and lies
  // within this distance of j, j bonded to k (the hydroxyl O-H): its rule applies without a bond i-j, so i and k keep
  // their non-bonded terms (M-O-H bends of ClayFF, Cygan, Liang, Kalinichev 2004)
  double angle_contacts = 0;
  // the special_bonds keyword LAMMPS has for this force field ("amber", "dreiding"; styles "special" in the file)
  std::string special_style;
  // GROMACS topology form of the van der Waals parameters: "" σ/ε (comb-rule 2), "c6c12" C6/C12 (GROMOS, comb-rule 1)
  std::string gromacs_lj;
  // DREIDING's hydrogen bond (JSON "hbonds"; LAMMPS hbond/dreiding/lj): E = S(r) ε [5 (σ/r)¹² − 6 (σ/r)¹⁰] cos^n θ between a
  // donor D and an acceptor A (r = D–A), θ the D–H···A angle at a hydrogen H bonded to D, only for θ > the angle cut-off;
  // S switches from inner to outer. After typing, N/O/F atoms take their "_hd" variant when they carry a hydrogen (which
  // becomes the hydrogen type), else their "_ha" variant, where the force field has one.
  struct HBondDef {
    int power = 4;
    double inner = 6.0, outer = 6.5, angle = 90;
    struct Term { std::string donor, acceptor, hydrogen; double eps = 0, sigma = 0; int n = 4; };
    std::vector<Term> terms;
  } hbonds;
  // Parameters by analogy for terms the file lacks (the typing file's "analogies", as AmberTools' parmchk2 fills
  // missing GAFF terms): type → types whose parameters stand in, most similar first. Every term found this way is
  // listed in ParamReport::estimated; nothing is filled silently.
  std::map<std::string, std::vector<std::string>> analogies;
  std::string analogy_source;
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
  // Charge keys (typing file "charge_rules"): where a force field's charges come from bond increments between keys finer
  // than its types (OPLS 2005: the ether C and the alcohol C are both CT, but their increments are keyed 181 and 157),
  // each atom gets a key by these rules (same matching as the typing rules) and bond_increments are looked up by keys
  std::vector<TypingRule> charge_typing;
  std::vector<std::pair<std::string, std::string>> typing_pairs;   // conjugated pairs (GAFF cc/cd, ...), see typing.hpp
  bool typing_ordered = false;
  bool typing_unknown_untyped = false;   // rules may name types this file lacks: their atoms end up untyped
  // United-atom force field (typing file "united_atom": true; GROMOS, TraPPE-UA, CHARMM19): hydrogens on carbon are
  // part of their carbon's site, so an all-atom structure is converted before typing (prepare_for_forcefield)
  bool united_atom = false;
  std::vector<int> united_atom_hosts{6};   // elements whose hydrogens fold in ("united_atom_hosts"; O for mW water)
  // Shell models (typing file "shells": core type → shell type): each core gets a shell particle on it, bonded by the
  // file's core-shell spring and named for its shell type (prepare_for_forcefield)
  std::map<std::string, std::string> shells;
  // Ionic solids (typing file "bonds": "defined"): the builder's neighbour bonds are not bonds of the model; after
  // typing, only bonds the force field has a term for stay (core-shell springs, O-H of water and hydroxyls)
  bool keep_defined_bonds = false;
  // Types that are never bonded (typing file "unbonded_types"): IFF's metals, 12-6 Lennard-Jones atoms; the neighbour
  // bonds perception gives them are dropped after typing (prepare_for_forcefield)
  std::vector<std::string> unbonded_types;
  // Type pairs that never interact (typing file "exclude_pairs": [["C", "C"]]; LAMMPS neigh_modify exclude type): a
  // graphene sheet held in place, as the source excludes its carbons from each other
  std::vector<std::pair<std::string, std::string>> exclude_type_pairs;
  // Coarse-grained force fields (typing file "coarse_grained": true): sites are beads, typed by name. "beads" and
  // "bead_groups" map an all-atom structure onto beads first (map_to_beads; SDK). bead_templates: named bead SMILES
  // from the force field's sources (force-field file "bead_templates"), for build_beads.
  bool coarse_grained = false;
  std::vector<BeadRule> bead_rules;
  std::vector<BeadGroup> bead_groups;
  std::map<std::string, std::string> bead_templates;
  // Molecules given term by term (force-field file "molecule_templates": a caps-martini-molecules file, GROMACS units;
  // Martini 3's solvents, ions, small molecules ...): built with their explicit topology, and recognised again in a
  // structure without one (same bead types in the same order, same bonds). Loaded with the force field.
  std::string molecule_templates_path;
  std::shared_ptr<const Json> molecule_templates;
  // Every pair's Lennard-Jones terms from a table (force-field file "pair_table": gzip text, "T name mass" and
  // "P a b sigma(nm) epsilon(kJ/mol)"; Martini 3's 355,746 pairs): no mixing rule, a pair not listed has no LJ
  std::string pair_table;
  // stiff bonds for constraints of molecule templates (kJ/mol/nm², "constraint_k")
  double constraint_kj = 1e6;
  // Martini proteins (typing file "martini_protein": the model's JSON): an all-atom protein becomes beads with the
  // model's explicit topology (martini22_protein; DSSP for the secondary structure)
  std::string martini_protein;
  // Martini 3 small molecules (typing file "martini_small_molecules"): all-atom molecules matched by graph to the model's
  // CHARMM residues become their beads (martini3_small_molecules)
  std::string martini_small;
  // Bond-order variants (DREIDING): a base type may have variants that differ only in which bonds get which force
  // constant (moltemplate's C_2 / C_2_b1 / C_2_b2, C_R / C_R_b1; the other file's C_2 / C_2S, C_R / C_RS). After the rules,
  // each conjugated system takes the variants that make every bond's constant equal bond_k_per_order x its bond order
  // (DREIDING: k = 700 n kcal/mol/Å², 350 n in LAMMPS's harmonic K).
  std::map<std::string, std::vector<std::string>> type_variants;
  double bond_k_per_order = 0;
  double bond_conjugated_single = 1.0;   // the order a force field gives a single bond between two conjugated atoms
  bool typing_pairs_double_same = false;   // pairs keep one type across a double bond (CGenFF CG2DC1/2), not GAFF's   // rules are an ordered list (antechamber): the first match is intended, not ambiguous
  std::string typing_source;
  std::vector<std::string> typing_files;   // rule files loaded (each once: a recipe may name the file the force field already names)
  std::vector<std::string> notes;

  const FFType* type(const std::string& name) const;
};

FFDef load_forcefield(const std::string& path);
// Charges "increments" for a force field whose charge_increments_from names another (a file beside ff_path): that
// force field types the structure by its own rules and gives its charges from its bond increments, which balance bond
// by bond; they are used with this force field's types and parameters (OPLS-AA 2024 with OPLS 2005's charges). Throws
// FFError with the reason when it cannot (atoms it cannot type, a bond without an increment).
std::vector<double> companion_charges(const System& s, const FFDef& def, const std::string& ff_path, std::string* note = nullptr);
// Reads a "caps-typing" rules file and appends its rules to the force field (types must exist in it).
void load_typing(FFDef& ff, const std::string& path);
void save_forcefield(const FFDef& ff, const std::string& path);
// Appends the overlay's types and rules after the base's (so they win), replacing types of the same name.
void merge_forcefield(FFDef& base, const FFDef& overlay);

// Gap filling: the bond, angle and dihedral rules of `donor` usable where `base` defines nothing (impropers are never
// borrowed — an absent improper is not a gap). For a moltemplate OPLS-AA base (numbered types 136_bCT_aCT_dCT_iCT) the
// donor's class names become the base's patterns; a rule on a class the base lacks is dropped (nothing is guessed).
// Rules are renamed "filled: …" and commented with `source`. Prepend them with prepend_fill: the last matching rule
// wins, so the base's own rules still decide wherever they apply.
struct GapFill {
  FFDef rules;
  size_t kept = 0, dropped = 0;
};
GapFill gap_fill_rules(const FFDef& base, const FFDef& donor, const std::string& source);
void prepend_fill(FFDef& def, const FFDef& fill);
// A moltemplate OPLS-AA type name, 136_bCT_aCT_dCT_iCT: its bond, angle, dihedral and improper classes
bool opls_classes(const std::string& name, std::array<std::string, 4>& classes);

// moltemplate force-field files (.lt): "In Init" styles and special_bonds, "Data Masses", "In Charges",
// pair / bond / angle / dihedral / improper coefficients, the "By Type" rules and their symmetry plugins, "replace"
// aliases, and the type descriptions from the comments.
FFDef import_moltemplate(const std::string& path);
// AMBER frcmod (MASS, BOND, ANGL, DIHE with multi-term dihedrals, IMPR, NONB as R*/2 and ε): rules in CAPS's styles,
// X as the wildcard. Throws FFError when nothing is found.
FFDef import_frcmod(const std::string& path);
// GROMACS [ atomtypes ] (σ ε, or C6 C12 with combination rule 1), [ bondtypes ] 1, [ angletypes ] 1 and 5 (the
// Urey–Bradley part noted), [ dihedraltypes ] 1, 9 (fourier), 4 (periodic improper) and 2 (harmonic improper) of a
// .itp or .top, in kcal/mol and Å; other functions are listed in references as not imported.
FFDef import_gromacs_params(const std::string& path);

// DL_FIELD force-field libraries: NAME.par (parameters), NAME.sf (atom types; templates are not converted yet)
// and NAME.bci (bond charge increments), as shipped in DL_FIELD's lib/ directory. Parameters are converted to
// CAPS units (kcal/mol, Å, degrees) and LAMMPS functional forms as DL_FIELD itself writes them for LAMMPS.
// Families not yet interpreted throw FFError naming the family.
FFDef import_dlfield(const std::string& par, const std::string& sf = "", const std::string& bci = "");

bool glob_match(const std::string& pattern, const std::string& text);
std::string glob_escape(const std::string& name);   // a pattern that matches exactly this name

struct ParamReport {
  std::vector<std::string> missing;    // interactions with no matching rule (each once, with an example)
  std::vector<std::string> estimated;  // interactions given the parameters of analogous types (FFDef::analogies), each once
  int estimated_terms = 0;             // how many interactions (all of them, not each kind once)
  std::map<std::string, int> used;     // rule name → interactions
  std::vector<std::string> charge_keys;   // per atom, when the charges came from bond increments between charge keys
  std::vector<std::string> notes;
  bool complete() const { return missing.empty(); }
};

// A structure as the force field describes it. A shell-model force field gets a shell on every core it types (at the
// core's place, bonded to it, named for the shell type so its rule types it). For a united-atom force field an all-atom structure (hydrogens on
// carbon) becomes united-atom: each such hydrogen folds into its carbon (CH, CH2, CH3 sites, polar hydrogens kept,
// resolution.hpp united_atom). Charges "auto" or "gasteiger" are computed on the all-atom structure first (Gasteiger,
// else QEq) and summed into the sites, so the caller then uses charges "keep" (charges is updated). Returns a note
// when it converted, "" otherwise.
std::string prepare_for_forcefield(System& s, const FFDef& ff, std::string& charges);
// Does the force field change the structure before typing (united atom, shells, ionic bonds)?
inline bool needs_prepare(const FFDef& ff) {
  return ff.united_atom || !ff.shells.empty() || ff.keep_defined_bonds || !ff.unbonded_types.empty() || ff.coarse_grained || !ff.bead_rules.empty() || !ff.bead_groups.empty() ||
         !ff.martini_protein.empty();
}

// "N atoms match no typing rule of FF", or for a coarse-grained force field that maps atoms onto beads, that they are in
// molecules its bead fragments cannot cover.
std::string untyped_message(const FFDef& ff, const System& s, int untyped);

// A bead structure for a coarse-grained force field: `text` is one of its bead templates by name, or bead SMILES.
// Bond lengths come from the force field's bond terms (typed by bead name), masses from its types; the result is a
// start for a relax.
System build_bead_molecule(const std::string& text, const FFDef& ff, uint64_t seed = 1);

// The force field's named molecules: bead templates (bead SMILES) and molecule templates (a short description:
// "N beads, source file"), by name.
std::map<std::string, std::string> bead_template_list(const FFDef& ff);
bool has_bead_template(const FFDef& ff, const std::string& name);

// Give a structure without an explicit topology its molecule templates' topology: every molecule (by molecule id, else
// bonded component) whose bead types in order and bonds are a template's. Returns how many molecules were recognised;
// `unmatched` gets the others' first bead types. Templates with the same beads and bonds (Martini 3's ions) are told apart by
// the structure's charges when it has them, else the first by name is taken. Does nothing without molecule templates.
int recognise_molecule_templates(System& s, const FFDef& ff, std::vector<std::string>* unmatched = nullptr);

// Build the evaluator force field for a structure whose atoms carry force-field type names (one per atom).
// Charges: `charges` = "types" (from the force field: type charges and / or bond increments; error if neither
// gives a charge), "keep" (the structure's) or
// "gasteiger" (C/H/N/O only). Bonded interactions are generated from the bonds as moltemplate does. Throws FFError
// listing every missing parameter unless `allow_missing` (then the report lists them and those terms are left out).
ForceField parameterize(const System& s, const FFDef& ff, const std::vector<std::string>& types, const std::string& charges,
                        ParamReport* report = nullptr, bool allow_missing = false);

}  // namespace caps
