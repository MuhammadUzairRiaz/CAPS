#!/usr/bin/env python3
"""MARTINI 2 polymers for CAPS: PEO / PEG and polystyrene, an overlay on martini2-marrink2007 built from the papers (PDFs
in ~/project/ff-sources/papers/cg, never in the repository):

  PEO   H. Lee, A. H. de Vries, S. J. Marrink, R. W. Pastor, J. Phys. Chem. B 113, 13186 (2009), Table 1 and section 2:
        one bead per -CH2-O-CH2- (the monomer's C-O-C), named PEO; PEO-PEO sigma 4.3 A, epsilon 3.375 kJ/mol (as SNda),
        PEO-W 4.7 A, 4.0 kJ/mol; bond 3.30 A, K_b 17000 kJ mol-1 nm-2; cosine angle 130 deg, K 85 kJ/mol; dihedral
        sum K (1 + cos(n phi - phi0)): 1.96 kJ/mol n 1 at 180 deg, 0.18 n 2, 0.33 n 3, 0.12 n 4 at 0 deg. With the other
        bead types PEO takes SNa's interactions, as the paper proposes ("the SNa type will be more appropriate").
  PS    G. Rossi, L. Monticelli, S. R. Puisto, I. Vattulainen, T. Ala-Nissila, Soft Matter 7, 698 (2011), Tables 2 and 4,
        the B-mapping (a backbone bead per styrene unit, bonded along the chain; three ring beads): backbone SC1 (named SC1PS, so cyclohexane's SC1 keeps its own bonds), ring
        beads SC4PS (SC4 except SC4PS-SC4PS: sigma 0.41 nm, epsilon 2.4 kJ/mol); bonds B-B 0.25 nm / 8000 kJ mol-1 nm-2,
        B-R 0.217 nm and R-R 0.27 nm constrained in the paper (stiff bonds here, 1e6 kJ mol-1 nm-2: CAPS has no
        constraints in a force field), harmonic angles B-B-B 170 deg / 25, R-B-B 125 / 45, B-R-R 156 / 200 kJ/mol;
        the ring's own triangle has no angle terms.

Mapping all-atom structures (typing file "beads", an exact cover of the heavy atoms): PEO's -CH2-O-CH2- monomers; for
polystyrene the backbone CH and one CH2 per unit as the backbone bead and the phenyl's six carbons as three pairs. The
paper's mappings share backbone atoms between beads (each shared atom counted half); CAPS's exact cover cannot, so the
backbone bead's centre sits a little off the paper's: relax the mapped structure before a run.

usage: make_martini2_polymers.py [DATA_DIR]
"""
import json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.dirname(HERE)), "data")
FID = "martini2-polymers"
KJ = 4.184
LEE = "H. Lee, A. H. de Vries, S. J. Marrink, R. W. Pastor, J. Phys. Chem. B 113, 13186 (2009)"
ROSSI = "G. Rossi, L. Monticelli, S. R. Puisto, I. Vattulainen, T. Ala-Nissila, Soft Matter 7, 698 (2011)"
STIFF = 1.0e6   # kJ mol-1 nm-2, a constraint written as a stiff bond


def main():
    base = json.load(open(os.path.join(DATA, "forcefields", "martini2-marrink2007.json")))
    bp = {}
    for p in base["pairs"]:
        a, b = (p["match"][0], p["match"][0]) if len(p["match"]) == 1 else p["match"]
        bp[(a, b)] = bp[(b, a)] = p["params"]
    names = [t["name"] for t in base["atom_types"]]
    types = [{"name": "PEO", "element": "", "mass": 72.0, "charge": 0.0, "description": "PEO monomer -CH2-O-CH2- (Lee et al. 2009)"},
             {"name": "SC4PS", "element": "", "mass": 72.0, "charge": 0.0, "description": "polystyrene ring bead (Rossi et al. 2011, SC4 PS)"},
             {"name": "SC1PS", "element": "", "mass": 72.0, "charge": 0.0, "equivalence": {"vdw": "SC1"},
              "description": "polystyrene backbone bead: SC1's non-bonded terms (Rossi et al. 2011), its own bonded terms (cyclohexane's SC1 keeps MARTINI 2.0's)"}]
    pairs = [{"name": "PEO-PEO", "match": ["PEO"], "params": [3.375 / KJ, 4.3], "comment": "Lee 2009 Table 1: 4.3 A, 3.375 kJ/mol"}]
    for t in names:
        if t in ("P4", "BP4"):
            pairs.append({"name": f"PEO-{t}", "match": ["PEO", t], "params": [4.0 / KJ, 4.7], "comment": "Lee 2009 Table 1: PEO-W 4.7 A, 4.0 kJ/mol"})
        else:
            eps, sig = bp[("SNa", t)]
            pairs.append({"name": f"PEO-{t}", "match": ["PEO", t], "params": [eps, sig], "comment": "as SNa (Lee 2009, section 2)"})
    pairs.append({"name": "PEO-SC4PS", "match": ["PEO", "SC4PS"], "params": bp[("SNa", "SC4")], "comment": "as SNa-SC4"})
    pairs.append({"name": "SC4PS-SC4PS", "match": ["SC4PS"], "params": [2.4 / KJ, 4.1], "comment": "Rossi 2011 Table 4: 0.41 nm, 2.4 kJ/mol"})
    for t in names:
        eps, sig = bp[("SC4", t)]
        pairs.append({"name": f"SC4PS-{t}", "match": ["SC4PS", t], "params": [eps, sig], "comment": "as SC4 (Rossi 2011, section 2.4)"})
    hb = lambda k: k / KJ / 100 / 2   # GROMACS (1/2) k (r - r0)^2 in kJ mol-1 nm-2 -> LAMMPS K, kcal/mol/A^2
    bonds = [{"name": "PEO-PEO", "match": ["PEO", "PEO"], "style": "harmonic", "params": [hb(17000), 3.30], "comment": "Lee 2009 Table 1"},
             {"name": "SC1PS-SC1PS", "match": ["SC1PS", "SC1PS"], "style": "harmonic", "params": [hb(8000), 2.5], "comment": "Rossi 2011 Table 2, B-B (B-mapping)"},
             {"name": "SC1PS-SC4PS", "match": ["SC1PS", "SC4PS"], "style": "harmonic", "params": [hb(STIFF), 2.17],
              "comment": "Rossi 2011 Table 2, B-R 0.217 nm, constrained (a stiff bond here)"},
             {"name": "SC4PS-SC4PS", "match": ["SC4PS", "SC4PS"], "style": "harmonic", "params": [hb(STIFF), 2.7],
              "comment": "Rossi 2011 Table 2, R-R 0.27 nm, constrained (a stiff bond here)"}]
    angles = [{"name": "PEO-PEO-PEO", "match": ["PEO", "PEO", "PEO"], "style": "cosine/squared", "params": [85 / KJ / 2, 130.0],
               "comment": "Lee 2009 Table 1: 130 deg, 85 kJ/mol, (1/2) K (cos - cos0)^2"},
              {"name": "SC1PS-SC1PS-SC1PS", "match": ["SC1PS", "SC1PS", "SC1PS"], "style": "harmonic", "params": [25 / KJ / 2, 170.0], "comment": "Rossi 2011 Table 2, B-B-B"},
              {"name": "SC4PS-SC1PS-SC1PS", "match": ["SC4PS", "SC1PS", "SC1PS"], "style": "harmonic", "params": [45 / KJ / 2, 125.0], "comment": "Rossi 2011 Table 2, R-B-B"},
              {"name": "SC1PS-SC4PS-SC4PS", "match": ["SC1PS", "SC4PS", "SC4PS"], "style": "harmonic", "params": [200 / KJ / 2, 156.0], "comment": "Rossi 2011 Table 2, B-R-R"},
              {"name": "SC4PS-SC4PS-SC4PS", "match": ["SC4PS", "SC4PS", "SC4PS"], "style": "harmonic", "params": [0.0, 60.0],
               "comment": "the ring triangle: no angle term (constrained bonds)"}]
    dihedrals = [{"name": "PEO-PEO-PEO-PEO", "match": ["PEO", "PEO", "PEO", "PEO"], "style": "fourier",
                  "params": [4, 1.96 / KJ, 1, 180.0, 0.18 / KJ, 2, 0.0, 0.33 / KJ, 3, 0.0, 0.12 / KJ, 4, 0.0],
                  "comment": "Lee 2009 Table 1: K (1 + cos(n phi - phi0)), 1.96 / 0.18 / 0.33 / 0.12 kJ/mol for n 1-4"}]
    ps = lambda n: "".join(f"[SC1PS]([SC4PS]{k % 9 + 1}[SC4PS][SC4PS]{k % 9 + 1})" for k in range(n))
    doc = {"format": "caps-forcefield", "format_version": 1, "name": "MARTINI 2 polymers (PEO, polystyrene)",
           "source": LEE + ", Table 1; " + ROSSI + ", Tables 2 and 4 (B-mapping)", "references": [LEE, ROSSI], "extends": "martini2-marrink2007.json",
           "atom_types": types, "pairs": pairs, "bonds": bonds, "angles": angles, "dihedrals": dihedrals, "impropers": [],
           "bead_templates": {**base["bead_templates"], "PEO9": "[PEO]" * 9, "PEO18": "[PEO]" * 18, "PEO36": "[PEO]" * 36, "PS10": ps(10)},
           "notes": ["built by CAPS (bench/ff/make_martini2_polymers.py) from the papers' tables as an overlay on martini2-marrink2007",
                     "PEO interacts with other beads as SNa (the paper's proposal); polystyrene's ring beads as SC4 except among themselves",
                     "polystyrene's constrained bonds are stiff bonds (1e6 kJ mol-1 nm-2): use a short time step or constrain them",
                     "angle styles: cosine/squared (MARTINI, PEO) and harmonic (polystyrene, as the paper)"],
           "version": "2009 / 2011", "status": "validated", "typing": f"../typing/{FID}.typing.json"}
    json.dump(doc, open(os.path.join(DATA, "forcefields", FID + ".json"), "w"), indent=1, ensure_ascii=False)
    btyp = json.load(open(os.path.join(DATA, "typing", "martini2-marrink2007.typing.json")))
    typing = {"format": "caps-typing", "version": 1, "forcefield": doc["name"], "coarse_grained": True, "unknown_types": "untyped",
              "description": "CAPS rules for the MARTINI 2 polymers (bench/ff/make_martini2_polymers.py): beads by name, and all-atom PEO and "
                             "polystyrene mapped onto them",
              "beads": [{"type": "PEO", "smarts": "[CH2X4][OX2][CH2X4]", "description": "-CH2-O-CH2- (Lee 2009: one bead per C-O-C)"},
                        {"type": "PEO", "smarts": "[CX4;H2,H3][CH2X4][OX2]", "description": "a chain end's C-C-O (CAPS's choice; the paper maps C-O-C)"},
                        {"type": "SC1PS", "smarts": "[CX4;H2,H3][CX4;H1,H2;$(*c)]", "description": "polystyrene backbone CH2-CH (B-mapping's backbone bead)"},
                        {"type": "SC4PS", "smarts": "cc", "description": "two phenyl carbons (a ring bead)"}],
              "rules": btyp["rules"] + [{"type": t, "smarts": "*", "atom_name": t, "priority": 0, "description": f"a bead named {t}"}
                                        for t in ("PEO", "SC4PS", "SC1PS")]}
    json.dump(typing, open(os.path.join(DATA, "typing", FID + ".typing.json"), "w"), indent=1, ensure_ascii=False)
    e = {"id": FID, "name": doc["name"], "version": doc["version"], "references": doc["references"], "origin": "built from source",
         "source": doc["source"] + "; overlay on martini2-marrink2007; bench/ff/make_martini2_polymers.py",
         "license": "published parameters of the papers", "status": "validated",
         "notes": "PEO (Lee 2009) and polystyrene (Rossi 2011, B-mapping) on MARTINI 2.0; all-atom PEO and PS mapped onto them",
         "file": FID + ".json", "typing": {"rules": f"typing/{FID}.typing.json", "evidence": "beads by name; all-atom PEO and polystyrene by fragment"},
         "counts": {"atom_types": len(types), "pairs": len(pairs), "bonds": len(bonds), "angles": len(angles), "dihedrals": len(dihedrals), "impropers": 0},
         "group": "cg", "aliases": ["martini-polymers", "martini2-polymers"]}
    cat_p = os.path.join(DATA, "forcefields", "catalogue.json")
    cat = json.load(open(cat_p))
    cat["forcefields"] = [x for x in cat["forcefields"] if x["id"] != FID] + [e]
    json.dump(cat, open(cat_p, "w"), indent=1, ensure_ascii=False)
    print(f"{FID} written: {len(pairs)} pairs")


if __name__ == "__main__":
    main()
