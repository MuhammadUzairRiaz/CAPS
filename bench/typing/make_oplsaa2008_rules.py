#!/usr/bin/env python3
"""CAPS typing rules for OPLS-AA as BOSS 4.8's 2008 table numbers it (moltemplate oplsaa2008.lt).

The chemistry is the same as the 2024 rules (bench/typing/make_oplsaa_rules.py); only the numbers differ (the 2008
table's alkane CH3 is 80, CH2 81, H 85, where 2024 has 135, 136, 140). Each 2024 type the rules use is paired with the
2008 type describing the same group: same element, the same charge where the two tables agree, the bonded class that
corresponds (2024 CT / HC / CA / OS ... with 2008 013 / 046 / 048 / 020 ...) and the same description. Where the tables
give a group different charges (2008's phenol, sulfides, esters, silanes), the whole group takes the 2008 types, so its
charges still add up (the Field report's charge-balance check verifies this). Groups the 2008 table has no numbers for
(OPLS/2020 tertiary alcohols, epoxides, silanols, silyl ethers, dialkyl carbonates, phenylsilane) keep a rule naming no
type, so their atoms are reported untyped instead of falling to a more general rule with the wrong charges.

usage: make_oplsaa2008_rules.py [DATA_DIR]
"""
import json, os, re, sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "data")

# 2024 bonded class → 2008 bonded class
CLASS = {"CT": "013", "HC": "046", "CA": "048", "HA": "049", "C~": "003", "O~": "004", "OH": "005", "HO": "007", "OS": "020", "N~": "024",
         "H~": "045", "CM": "047", "S~": "016", "SH": "015", "HS": "017", "NT": "044", "Cl": "021", "F~": "001", "Br": "065", "I~": "066",
         "CZ": "019", "NZ": "018", "SY": "079", "OY": "023", "O2": "052", "N3": "053", "H3": "054", "NO": "102", "ON": "103", "CW": "084",
         "CS": "087", "NA": "057", "C!": "086", "CB": "048", "OA": "020", "tipH": "032", "tipO": "031", "CF": "013"}

# pairings the automatic match gets wrong or cannot make, read from both tables' descriptions (2024 → 2008; None: no 2008 type)
MANUAL = {
    "135": "80", "139": "84", "140": "85",                      # alkane CH3, >C<, alkane H
    "157": "99",                                                # primary alcohol CH2 / CH3 (not the diol type)
    "213": "155", "214": "156", "215": "157", "216": "158",     # disulfide carbons (same charges as the sulfide ones)
    "546": "487",                                               # pyrrole H2
    "756": "697", "757": "698", "758": "699",                   # alkyl nitrile carbons
    "903": "733", "904": "734", "907": "737", "909": "739", "910": "740", "912": "742",   # amines
    "956": "786",                                               # alkyl fluoride F
    "158": "100", "159": "101", "160": "96",                    # secondary / tertiary alcohols (2008's R3COH with its alcohol O)
    "166": "108", "167": "109", "168": "110",                   # phenol (2008 charges 0.15 / -0.585 / 0.435)
    "199": "141", "202": "144", "208": "150",                   # anisole C, sulfide S, tertiary thiol C
    "209": "151", "210": "152", "211": "153", "212": "154",     # sulfide carbons
    "152": "801", "153": "802", "972": "802", "973": "803", "974": "804",   # alkyl chlorides
    "976": "806", "978": "808",                                 # alkyl bromides
    "957": "787", "958": "788", "959": "789", "960": "790",     # alkyl fluorides
    "718": "659", "719": "660", "729": "670", "730": "671",     # fluoro- and bromobenzene
    "465": "406", "466": "407", "468": "409", "469": "410", "471": "412",   # esters (2008's 0.51 / -0.43, methyl ester, aryl ester)
    "906": "736", "908": "738", "911": "741",                   # amine CH2 carbons, H on carbons next to amine N
    "1060": "866", "1061": "867", "1062": "868", "1063": "869", "1064": "870",   # alkylsilanes (870 is H on Si)
    "1065": "871", "1066": "872", "1067": "873", "1068": "874",
    "1069": None, "1070": None, "1071": None, "1072": None, "1073": None, "1074": None, "1075": None, "1076": None, "1077": None, "1078": None,
    "1025": None, "1026": None, "1027": None, "1028": None, "1029": None,   # epoxides
    "789": None,                                                # methoxy carbon of dialkyl carbonates
}


def short(n):
    return n.split("_b")[0]


def desc(t):
    s = t.get("description", "")
    s = re.sub(r"moltemplate @atom:\S+ · ", "", s)
    return re.sub(r"\s*·\s*\S+$", "", s).replace('"', "").strip()


def words(s):
    return {w for w in re.findall(r"[a-z0-9]+", s.lower().replace("-", " ")) if len(w) > 1}


def bclass(t):
    return re.search(r"_b([^_]+)_", t["name"]).group(1)


if __name__ == "__main__":
    d4 = json.load(open(os.path.join(DATA, "forcefields", "oplsaa2024-moltemplate.json")))
    d8 = json.load(open(os.path.join(DATA, "forcefields", "oplsaa2008-moltemplate.json")))
    src = json.load(open(os.path.join(DATA, "typing", "oplsaa2024-moltemplate.typing.json")))
    t4 = {short(t["name"]): t for t in d4["atom_types"]}
    t8 = {short(t["name"]): t for t in d8["atom_types"]}
    pair, unmatched = {}, []
    for u in sorted({r["type"] for r in src["rules"]}, key=lambda x: int(re.match(r"\d+", x).group())):
        if u in MANUAL:
            pair[u] = MANUAL[u]
            continue
        a = t4[u]
        best = None
        for k, b in t8.items():
            if b.get("element") != a.get("element") or a.get("charge") is None or b.get("charge") is None or "DON'T USE" in desc(b):
                continue
            if abs(b["charge"] - a["charge"]) > 1e-4 or CLASS.get(bclass(a)) != bclass(b):
                continue
            score = len(words(desc(a)) & words(desc(b)))
            if best is None or score > best[0]:
                best = (score, k)
        pair[u] = best[1] if best else None
        if best is None:
            unmatched.append(u)
    rules = []
    for r in src["rules"]:
        m = pair[r["type"]]
        rules.append(dict(r, type=m if m else "none-in-2008", description=r.get("description", "") + ("" if m else " (no 2008 type: untyped)")))
    doc = {"format": "caps-typing", "version": 1, "forcefield": "OPLS-AA (BOSS 4.8, 2008)",
           "description": "The OPLS-AA 2024 rules' chemistry with the 2008 table's numbers (bench/typing/make_oplsaa2008_rules.py); groups the 2008 "
                          "table has no numbers for are reported untyped.",
           "unknown_types": "untyped", "rules": rules}
    out = os.path.join(DATA, "typing", "oplsaa2008-moltemplate.typing.json")
    json.dump(doc, open(out, "w"), ensure_ascii=False, indent=1)
    d8["typing"] = "../typing/oplsaa2008-moltemplate.typing.json"
    json.dump(d8, open(os.path.join(DATA, "forcefields", "oplsaa2008-moltemplate.json"), "w"), ensure_ascii=False, indent=1)
    cat_p = os.path.join(DATA, "forcefields", "catalogue.json")
    cat = json.load(open(cat_p))
    for e in cat["forcefields"]:
        if e["id"] == "oplsaa2008-moltemplate":
            e["typing"] = {"rules": "typing/oplsaa2008-moltemplate.typing.json", "evidence": "the OPLS-AA 2024 rules mapped to the 2008 numbering by element, charge, class and description (bench/typing/make_oplsaa2008_rules.py)"}
    json.dump(cat, open(cat_p, "w"), ensure_ascii=False, indent=1)
    print(f"{out}: {len(rules)} rules; {sum(1 for v in pair.values() if v)} of {len(pair)} types paired; no 2008 type: {', '.join(sorted(k for k, v in pair.items() if not v))}")
    for u, v in sorted(pair.items(), key=lambda kv: int(re.match(r'\d+', kv[0]).group())):
        if v and u not in MANUAL:
            print(f"  {u:>5} {desc(t4[u])[:44]:44} -> {v:>4} {desc(t8[v])[:44]}")
