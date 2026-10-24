// CAPS Exchange, topology: bonds broken and made in a structure and its force field while a run goes on (stage E1 of
// docs/studies/vitrimers/caps_module_proposal.md). One change can break and make any bonds: a swap of partners
// (associative exchange: transesterification, disulfide metathesis, transimination, boronic ester and siloxane
// exchange, bead-spring swaps), a bond broken alone or made alone (dissociative exchange: Diels–Alder, hindered urea,
// reversible crosslinks and stickers). The force field follows:
//
//   bond terms        the broken bond's own parameters when the same pair is made again or a swap moves them to the
//                     new partner; otherwise the force field's own bond between the same two types
//   angles, torsions  every angle and proper torsion through a made bond, from the force field's own terms between the
//                     same types (a combination of types the structure had with no term gets none, as the force field
//                     gave it none); angles and torsions through a broken bond go
//   impropers         at each atom whose neighbours changed, from an atom of the same type with neighbours of the same
//                     types (harmonic, Fourier and inversion forms)
//   1-4 pairs         from the bond graph, with each pair's own scales when it had them
//   exclusions        1-2, and 1-3 and 1-4 when the force field excludes those (read from its own list), plus any pair it
//                     excluded for another reason
//   the structure     its bond list
//
// A change the force field's own terms cannot describe (a pair of types it never bonded, class II cross terms, combined
// bending–torsion terms, virtual sites) throws TopologyError before anything is changed: such chemistry needs the
// force field assigned again to the new structure (the run's force field replaced: DynamicsOptions::topology).
#pragma once
#include <array>
#include <cstdint>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include "caps/field.hpp"
#include "caps/system.hpp"

namespace caps {

struct TopologyError : std::runtime_error {
  explicit TopologyError(const std::string& m) : std::runtime_error(m) {}
};

class TopologyEditor {
 public:
  using Pair = std::array<uint32_t, 2>;
  // s and ff are edited in place and must outlive the editor; ff must be for s's atoms. Reads the force field's
  // exclusion rule and its terms as tables by type.
  TopologyEditor(System& s, ForceField& ff);
  bool bonded(uint32_t i, uint32_t j) const;
  const std::vector<std::vector<uint32_t>>& neighbours() const { return nb_; }
  // One change: every bond in `broken` (they must exist) broken and every bond in `made` (they must not) made, then the
  // terms, 1-4 pairs and exclusions set again. Throws TopologyError, with s and ff as they were, when the force field
  // cannot describe the result.
  void change(const std::vector<Pair>& broken, const std::vector<Pair>& made);
  // The same with atoms retyped first (an exchange whose atoms trade roles: an ester oxygen that becomes an alcohol's):
  // each takes its new type and charge, and every other bond of a retyped atom is made again from the table of its new
  // types, with the angles, torsions and impropers through it.
  struct Retype { uint32_t atom = 0; int type = 0; std::string name; double charge = 0; };
  void change(const std::vector<Pair>& broken, const std::vector<Pair>& made, const std::vector<Retype>& retype);
  // a1–a2 and b1–b2 become a1–b2 and b1–a2; a1 and b1 each keep their own bond parameters.
  void swap(uint32_t a1, uint32_t a2, uint32_t b1, uint32_t b2);
  // The bond made between the types of `like.i` and `like.j` from now on (stickers whose types the structure has not bonded
  // yet; replaces the force field's own bond for those types). Its atoms only name the types.
  void define_bond(const BondTerm& like);
  void define_bond(const BondX& like);
  // Terms by type learned from another structure typed by the same force field (a change's neighbourhood typed again: a
  // combination of types the structure never had, such as C–S–S–C once polysulfides exchange into a disulfide). Only
  // terms whose atoms are all marked complete (their neighbourhood whole in `frag`) are read, and only for type names
  // this force field has; combinations already known stay as they are. Bonds, angles and torsions; impropers are not
  // learned (a centre without one is counted in missing_impropers). Returns how many combinations were added.
  size_t learn(const System& frag, const ForceField& f, const std::vector<char>& complete);
  // Changes done, and impropers that could not be placed (a new environment of a centre the force field had no improper
  // for; none placed there).
  size_t changes() const { return changes_; }
  size_t missing_impropers() const { return missing_impropers_; }
  bool excludes_13() const { return x13_; }
  bool excludes_14() const { return x14_; }

 private:
  struct BondParams { int kind = -1; BondTerm b{}; BondX x{}; Class2Bond c{}; };   // kind 0 harmonic, 1 extended, 2 class II
  using Key2 = std::pair<int, int>;
  using Key3 = std::array<int, 3>;
  using Key4 = std::array<int, 4>;
  struct AngleSet { std::vector<AngleTerm> h; std::vector<AngleX> x; std::vector<UreyBradley> ub; };
  struct ImproperPattern { int kind; std::array<int, 4> slot; TorsionTerm t{}; HarmonicTorsion h{}; InversionTerm v{}; };   // slot: −1 the centre, else a neighbour's type
  using ImproperKey = std::pair<int, std::vector<int>>;

  // fresh: bonds made from the type table even when their own parameters are at hand (a retyped atom's)
  void change_impl(const std::vector<Pair>& broken, const std::vector<Pair>& made, const std::map<Pair, BondParams>& hint,
                   const std::set<Pair>& fresh = {});
  int type(uint32_t i) const { return ff_.type_index[i]; }
  Key2 key2(uint32_t i, uint32_t j) const;
  Key3 key3(uint32_t i, uint32_t j, uint32_t k) const;
  Key4 key4(uint32_t i, uint32_t j, uint32_t k, uint32_t l) const;
  ImproperKey improper_key(uint32_t c, const std::vector<std::vector<uint32_t>>& nb) const;
  BondParams bond_params(uint32_t i, uint32_t j) const;
  void graph_pairs(const std::vector<std::vector<uint32_t>>& nb, std::set<Pair>& p12, std::set<Pair>& p13, std::set<Pair>& p14) const;

  System& s_;
  ForceField& ff_;
  size_t n_;
  std::vector<std::vector<uint32_t>> nb_;
  bool x13_ = true, x14_ = true;
  bool p14_from_graph_ = false;                 // ff.pairs14 are the graph's 1-4 pairs (rebuilt on a change)
  std::set<Pair> extra_excluded_;               // excluded pairs the bond graph does not explain (kept)
  std::map<Pair, std::array<double, 2>> p14_scale_;     // per-pair 1-4 scales (lj, coul), merged force fields
  std::map<Key2, std::array<double, 2>> p14_type_scale_;
  std::map<Key2, BondParams> bond_tab_;
  std::set<Key2> bond_seen_;
  std::map<Key3, AngleSet> angle_tab_;
  std::set<Key3> angle_seen_;
  std::map<Key4, std::vector<TorsionTerm>> dihedral_tab_;
  std::set<Key4> dihedral_seen_;
  std::map<ImproperKey, std::vector<ImproperPattern>> improper_tab_;
  std::set<ImproperKey> improper_seen_;
  std::map<Pair, BondParams> stash_;            // parameters of bonds broken, for when they are made again
  size_t changes_ = 0, missing_impropers_ = 0;
};

}  // namespace caps
