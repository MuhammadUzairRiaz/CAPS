#!/usr/bin/env python3
"""CAPS typing rules and extensions for the library's miscellaneous force field (misc.json): four molecule-specific
parameter sets, each typed only where its charges add up.

  alkanes      O'Malley & Catlow, Phys. Chem. Chem. Phys. 15, 19024 (2013) (alkanes in silicalite): sp3 carbons with
               only carbon and hydrogen neighbours, CH3 / CH2 / CH as Cp1 / Cs1 / Ct1 (-0.3 / -0.2 / -0.1) with H1
               (+0.1), so every CHn group is neutral; a quaternary carbon is the generic C1 (CAPS extension: charge 0,
               the same -0.1 e per hydrogen). Methane has no type (its carbon would need -0.4) and stays untyped
  methanol     Plant, Maurin & Bell, J. Phys. Chem. B 111, 2836 (2007): CH3OH only
  chloroform   Ramsahye & Bell, J. Phys. Chem. B 109, 4738 (2005): CHCl3 only
  HFA-134a     Peguin, Kamath, Potoff & da Rocha, J. Phys. Chem. B 113, 178 (2009) (force field 2, OPLS-AA based):
               CF3-CH2F only. The source gives both carbons the type CT with different charges (0.534 in CF3, 0.002 in
               CH2F); CAPS splits them as CTf3 / CTf1 with CT's parameters

Force-field extensions (in misc.json's notes): type charges from the source's molecule templates; CTf3 / CTf1; zero
torsions for the alkanes (the source's alkane template switches torsions off) with Cp1 / Cs1 / Ct1 taking C1's
torsions as they take its other terms. The source scales HFA-134a's 1-4 pairs by 0.5 and methanol's by 1.0; one
factor applies to a force-field file, and misc.json keeps 1.0.

usage: make_misc_rules.py [SOURCE_LIB_DIR]
"""
import json, os, re, sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA = os.path.join(ROOT, "data")
SRC = sys.argv[1] if len(sys.argv) > 1 else "~/project/dl_f_4.13/lib"
FID = "misc"
REFS = ["A. J. O'Malley, C. R. A. Catlow, Phys. Chem. Chem. Phys. 15, 19024 (2013)",
        "D. F. Plant, G. Maurin, R. G. Bell, J. Phys. Chem. B 111, 2836 (2007)",
        "N. A. Ramsahye, R. G. Bell, J. Phys. Chem. B 109, 4738 (2005)",
        "R. P. S. Peguin, G. Kamath, J. J. Potoff, S. R. P. da Rocha, J. Phys. Chem. B 113, 178 (2009)"]

HCARB = "!$(*~[!#6;!#1])"                 # a carbon whose neighbours are carbons and hydrogens only
HFA_C3 = "[CX4](-[#9])(-[#9])(-[#9])-[CX4H2]-[#9]"
HFA_C1 = "[CX4H2](-[#9])-[CX4](-[#9])(-[#9])-[#9]"
RULES = [
    ("Cp1", "[CX4H3;%s]" % HCARB, "alkane CH3 carbon (O'Malley & Catlow)"),
    ("Cs1", "[CX4H2;%s]" % HCARB, "alkane CH2 carbon"),
    ("Ct1", "[CX4H1;%s]" % HCARB, "alkane CH carbon"),
    ("C1", "[CX4H0;%s]" % HCARB, "alkane quaternary carbon (generic C1)"),
    ("H1", "[#1][CX4;H1,H2,H3;%s]" % HCARB, "H on an alkane carbon"),
    ("C2", "[CX4H3]-[OX2H1]", "methanol carbon (Plant, Maurin & Bell)"),
    ("H2", "[#1][CX4H3]-[OX2H1]", "methanol methyl H"),
    ("O2", "[OX2H1]-[CX4H3]", "methanol oxygen"),
    ("HO2", "[#1][OX2H1]-[CX4H3]", "methanol hydroxyl H"),
    ("C3", "[CX4H1](-[#17])(-[#17])-[#17]", "chloroform carbon (Ramsahye & Bell)"),
    ("H3", "[#1][CX4H1](-[#17])(-[#17])-[#17]", "chloroform H"),
    ("Cl3", "[#17]-[CX4H1](-[#17])-[#17]", "chloroform Cl"),
    ("CTf3", HFA_C3, "HFA-134a CF3 carbon (Peguin et al., force field 2)"),
    ("CTf1", HFA_C1, "HFA-134a CH2F carbon"),
    ("F", "[#9;$([#9]-[CX4](-[#9])(-[#9])-[CX4H2]-[#9]),$([#9]-[CX4H2]-[CX4](-[#9])(-[#9])-[#9])]", "HFA-134a F"),
    ("HC", "[#1]%s" % HFA_C1, "HFA-134a H"),
]


def template_charges():
    """{(molecule, atom label): (type alias, charge)} from the source's molecule templates"""
    txt = open(os.path.join(SRC, "MISC_FF.sf")).read()
    out = {}
    for m in re.finditer(r"^MOLECULE (\S+).*?\n(.*?)^END MOLECULE", txt, re.S | re.M):
        for l in m.group(2).splitlines():
            w = l.split()
            if len(w) >= 3 and not w[0].startswith(("CONNECT", "DIHEDRAL", "#")):
                try:
                    out[(m.group(1), w[0])] = (w[1], float(w[2]))
                except ValueError:
                    pass
    return out


def note(ff, text):
    ff.setdefault("notes", [])
    if text not in ff["notes"]:
        ff["notes"].append(text)


if __name__ == "__main__":
    p = os.path.join(DATA, "forcefields", FID + ".json")
    ff = json.load(open(p))
    types = {t["name"]: t for t in ff["atom_types"]}
    by_alias = {a: t for t in ff["atom_types"] for a in t.get("aliases", [])}
    ch = template_charges()
    for (mol, label), (alias, q) in ch.items():
        t = by_alias.get(alias)
        if t is not None and t["name"] != "CT":   # CT: two charges, split below
            t["charge"] = q
    # HFA-134a's carbons: one type in the source, two charges
    for name, label, what in (("CTf3", "C1", "CF3"), ("CTf1", "C2", "CH2F")):
        q = ch[("HFA_134a", label)][1]
        if name not in types:
            t = {"name": name, "element": "C", "mass": types["CT"]["mass"], "charge": q,
                 "description": f"HFA-134a {what} carbon: CT's parameters, its template charge",
                 "equivalence": {"vdw": "CT", "bond": "CT", "angle": "CT", "dihedral": "CT"}}
            ff["atom_types"].append(t)
            types[name] = t
        types[name]["charge"] = q
    note(ff, "CAPS extension: CTf3 / CTf1 are HFA-134a's CF3 and CH2F carbons (charges 0.534 / 0.002 in the source's template, "
             "which gives both the type CT), with CT's parameters")
    # the generic alkane carbon: a quaternary centre, no hydrogens
    types["C1"]["charge"] = 0.0
    note(ff, "CAPS extension: C1 (a quaternary alkane carbon) has charge 0, the alkane template's -0.1 e per hydrogen that keeps each CHn group neutral")
    note(ff, "CAPS extension: type charges from the force field's molecule templates")
    # alkane torsions are off in the source: zero torsions, so CAPS sees them as given
    for n in ("Cp1", "Cs1", "Ct1"):
        types[n].setdefault("equivalence", {})["dihedral"] = "C1"
    have = {tuple(d["match"]) for d in ff["dihedrals"]}
    for m in (("C1", "C1", "C1", "C1"), ("H1", "C1", "C1", "C1"), ("H1", "C1", "C1", "H1")):
        if m not in have:
            ff["dihedrals"].append({"name": "-".join(m), "match": list(m), "style": "fourier", "params": [1, 0.0, 3, 0.0],
                                    "comment": "CAPS extension: zero (the source's alkane template has DIHEDRAL OFF)"})
    note(ff, "CAPS extension: zero alkane torsions (the source's alkane template switches torsions off)")
    note(ff, "1-4 pairs: the source scales HFA-134a's by 0.5 (OPLS-AA) and methanol's by 1.0; CAPS applies one factor per "
             "force-field file, 1.0 here, so HFA-134a's F...F and F...H 1-4 pairs count in full")
    ff["references"] = REFS
    ff["typing"] = f"../typing/{FID}.typing.json"
    json.dump(ff, open(p, "w"), ensure_ascii=False, indent=1)

    doc = {"format": "caps-typing", "version": 1, "forcefield": ff["name"], "unknown_types": "untyped",
           "description": "CAPS rules for the miscellaneous set (bench/typing/make_misc_rules.py): alkanes (hydrocarbon sp3 carbons), "
                          "methanol, chloroform and HFA-134a, each only where its template charges add up; other atoms are untyped.",
           "rules": [{"type": t, "smarts": s, "priority": 0, "description": d} for t, s, d in RULES]}
    out = os.path.join(DATA, "typing", FID + ".typing.json")
    json.dump(doc, open(out, "w"), ensure_ascii=False, indent=1)

    cat_p = os.path.join(DATA, "forcefields", "catalogue.json")
    cat = json.load(open(cat_p))
    for e in cat["forcefields"]:
        if e["id"] == FID:
            e["references"] = REFS
            e["notes"] = "alkanes (silicalite), methanol and chloroform (faujasite) and HFA-134a: molecule-specific sets, typed per molecule"
            e["typing"] = {"rules": f"typing/{FID}.typing.json", "evidence": f"CAPS rules (bench/typing/make_misc_rules.py): {len(RULES)} rules"}
    json.dump(cat, open(cat_p, "w"), ensure_ascii=False, indent=1)
    print(f"{out}: {len(RULES)} rules; charges: " + ", ".join(f"{t['name']} {t['charge']:+g}" for t in ff["atom_types"] if "charge" in t))
