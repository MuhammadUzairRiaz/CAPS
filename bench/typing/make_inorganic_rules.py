#!/usr/bin/env python3
"""CAPS typing rules for the library's inorganic and water force fields.

Compound-specific potential sets (binary halides and oxides, ternary oxides, glasses, zeolites, GaN) are typed from the
force field's own compound templates: each atom of a compound gets that compound's type when the structure holds the
compound's other elements (the rule's "requires"); where compounds nest (sodium silicate in soda-lime silicate) the rule
naming more elements wins. A composition with several potential sets takes the first reference set, except TiO2, which
takes the general set (rutile and anatase cannot be told apart by composition); the others are chosen by hand in the
Field step. Shell models declare their core → shell pairs ("shells"): CAPS puts a shell on every core it types.

CLAYFF (Cygan, Liang & Kalinichev, J. Phys. Chem. B 108, 1255 (2004)) and the Hill-Sauer zeolite field (Hill & Sauer,
J. Phys. Chem. 98, 1238 (1994); 99, 9536 (1995)) type by coordination and neighbours. Water models: O and H of water.

usage: make_inorganic_rules.py [SOURCE_LIB_DIR]
"""
import json, os, re, sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA = os.path.join(ROOT, "data")
SRC = sys.argv[1] if len(sys.argv) > 1 else "~/project/dl_f_4.13/lib"
FILES = {"binary_halides": "inorganic-binary-halides", "binary_misc": "inorganic-binary-misc", "binary_oxides": "inorganic-binary-oxides",
         "glass": "inorganic-glass", "ternary_oxides": "inorganic-ternary-oxides", "zeolite": "inorganic-zeolite"}
PREFER = {"inorganic-binary-oxides": {("O", "Ti"): "TiO2_1"}}   # composition → template, where the first would not do
RE_GLASS = {"Sc", "Y", "Lu", "La"}                               # the rare-earth aluminosilicate lists every rare earth; a glass holds one


def rule(t, prio, smarts, desc, **kw):
    r = {"type": t, "smarts": smarts, "priority": prio, "description": desc}
    r.update({k: v for k, v in kw.items() if v})
    return r


def templates(sf):
    txt = open(os.path.join(SRC, sf)).read()
    at = re.search(r"^ATOM_TYPE(.*?)^END ATOM_TYPE", txt, re.S | re.M).group(1)
    keys = {l.split()[1] for l in at.splitlines() if len(l.split()) >= 3 and not l.startswith("#")}
    out = []
    for m in re.finditer(r"^MOLECULE (\S+).*?\n(.*?)^END MOLECULE", txt, re.S | re.M):
        labels = []
        for l in m.group(2).splitlines():
            w = l.split()
            if len(w) >= 3 and w[0] in keys and w[0] not in labels:
                labels.append(w[0])
        if labels:
            out.append((m.group(1), labels))
    return out


def compound_rules(short, fid):
    ff = json.load(open(os.path.join(DATA, "forcefields", fid + ".json")))
    types = {t["name"]: t for t in ff["atom_types"]}
    springs = {tuple(sorted(b["match"])) for b in ff.get("bonds", [])}
    rules, shells, chosen = [], {}, {}
    everything = {types[l]["element"] for _, labels in templates(f"INORGANIC_{short}.sf") for l in labels if l in types and types[l].get("element")}
    for name, labels in templates(f"INORGANIC_{short}.sf"):
        if name.startswith("water"):
            continue
        if any(l not in types or "charge" not in types[l] for l in labels):
            continue   # a template naming a type without parameters or charges
        el = {l: types[l]["element"] for l in labels}
        comp = tuple(sorted(set(el.values())))
        pref = PREFER.get(fid, {}).get(comp)
        if comp in chosen or (pref and name != pref):
            continue   # the first reference set (or the preferred one) for each composition
        chosen[comp] = name
        paired = lambda c: c.endswith("c") and c[:-1] + "s" in labels and tuple(sorted((c, c[:-1] + "s"))) in springs
        cores = {l for l in labels if not (l.endswith("s") and paired(l[:-1] + "c"))}
        for l in labels:
            if l in cores and paired(l):
                shells[l] = l[:-1] + "s"
        for l in sorted(cores):
            e = el[l]
            others = sorted(set(el.values()) - {e})
            if name.startswith("rare-earth"):   # Al, Si, O and the one rare earth present
                others = sorted(({"Al", "Si", "O"} | ({e} if e in RE_GLASS else set())) - {e})
                if e == "O":
                    others = ["Al", "Si"]
            z = f"[{e}]"
            # a compound's own set only: none of the elements of the file's other compounds (BaTiO3 is not BaO + TiO2)
            own = set(others) | {e} | (RE_GLASS if name.startswith("rare-earth") else set())
            rules.append(rule(l, len(others), z, f"{e} in {name.replace('_', ' ')}", requires=others,
                              excludes=sorted(everything - own - {"H"})))
        for c, s in shells.items():
            if c in cores:
                rules.append(rule(s, 99, f"[{el[c]}]", f"shell of {c}", atom_name=s))
    if fid == "inorganic-ternary-oxides":   # LaNiO3: its own oxygen charge split (types O6cN / O6sN)
        rules.append(rule("O6cN", 5, "[O]", "O core in LaNiO3", requires=["La", "Ni"]))
        rules.append(rule("O6sN", 99, "[O]", "shell of O6cN", atom_name="O6sN"))
        shells["O6cN"] = "O6sN"
    if fid == "inorganic-zeolite":         # the acid faujasite's bridging hydroxyl: by its H, not by composition alone
        rules = [r for r in rules if r["type"] not in ("Ob2", "Hb2")]
        rules.append(rule("Ob2", 5, "[O;$(*[H])]", "bridging hydroxyl O of an acid zeolite", requires=["Al", "Si", "H"]))
        rules.append(rule("Hb2", 5, "[H][O]", "bridging hydroxyl H of an acid zeolite", requires=["Al", "Si"]))
    return rules, shells, chosen


def clayff():
    R = []
    AL_OCT, AL_TET, SUB = "[Al;!X4]", "[Al;X4]", "[Mg,Fe,Li]"
    R.append(rule("st", 0, "[Si]", "tetrahedral silicon"))
    R.append(rule("ao", 0, "[Al]", "octahedral aluminium")); R.append(rule("at", 1, "[Al;X4]", "tetrahedral aluminium"))
    R.append(rule("mgo", 0, "[Mg;!X0]", "octahedral magnesium")); R.append(rule("mgh", 1, "[Mg;!X0;!$(*~[O;H0])]", "magnesium of a hydroxide (all O are OH)"))
    R.append(rule("cao", 0, "[Ca;!X0]", "octahedral calcium")); R.append(rule("cah", 1, "[Ca;!X0;!$(*~[O;H0])]", "calcium of a hydroxide"))
    R.append(rule("feo", 0, "[Fe]", "octahedral iron")); R.append(rule("lio", 0, "[Li]", "octahedral lithium"))
    for ion, s in (("Na", "[Na;X0]"), ("K", "[K;X0]"), ("Cs", "[Cs;X0]"), ("Ca", "[Ca;X0]"), ("Ba", "[Ba;X0]"), ("Cl", "[Cl;X0]")):
        R.append(rule(ion, 2, s, f"aqueous {ion} ion"))
    R.append(rule("ob", 0, "[O;H0]", "bridging oxygen"))
    R.append(rule("obb", 1, f"[O;H0;$(*(~{AL_OCT})~{AL_OCT});$(*~[Si])]", "bridging O with two octahedral Al and one Si"))
    R.append(rule("obt", 2, f"[O;H0;$(*~{AL_TET})]", "bridging O with tetrahedral substitution (Al for Si)"))
    R.append(rule("obc", 3, f"[O;H0;$(*~{AL_TET});$(*(~{AL_OCT})~{AL_OCT})]", "bridging O with two octahedral Al and a tetrahedral Al"))
    R.append(rule("obo", 2, f"[O;H0;$(*~{SUB});$(*~[Si,Al])]", "bridging O with octahedral substitution"))
    R.append(rule("obs", 4, f"[O;H0;$(*~{AL_TET});$(*~{SUB})]", "bridging O with double substitution"))
    R.append(rule("oh", 1, "[O;H1]", "hydroxyl oxygen")); R.append(rule("ohs", 2, f"[O;H1;$(*~{SUB});$(*~[Al,Si])]", "hydroxyl O with substitution"))
    R.append(rule("ho", 1, "[H][O;H1]", "hydroxyl hydrogen"))
    R.append(rule("o*", 5, "[OX2H2]", "SPC water O")); R.append(rule("h*", 5, "[H][OX2H2]", "SPC water H"))
    return R


def hill_sauer():
    return [rule("Si", 0, "[Si]", "silicon"), rule("Al", 0, "[Al]", "aluminium"),
            rule("O_SS", 1, "[O;$(*(~[Si])~[Si])]", "Si-O-Si oxygen"), rule("O_AS", 1, "[O;H0;$(*(~[Si])~[Al])]", "Si-O-Al oxygen"),
            rule("Ob", 2, "[O;$(*(~[Si])~[Al]);$(*[H])]", "bridging hydroxyl oxygen Si-O(H)-Al"),
            rule("O_SH", 1, "[O;H1;$(*[Si]);!$(*[Al])]", "silanol oxygen"), rule("O_AH", 1, "[O;H1;$(*[Al]);!$(*[Si])]", "aluminol oxygen"),
            rule("Hb", 2, "[H][O;$(*~[Si]);$(*~[Al])]", "bridging hydroxyl hydrogen"),
            rule("H_OS", 1, "[H][O;$(*[Si]);!$(*[Al])]", "silanol hydrogen"), rule("H_OA", 1, "[H][O;$(*[Al]);!$(*[Si])]", "aluminol hydrogen")]


def water(o, h):
    return [rule(o, 0, "[OX2H2]", "water oxygen"), rule(h, 0, "[H][OX2H2]", "water hydrogen")]


def write(fid, rules, label, extra):
    ff_p = os.path.join(DATA, "forcefields", fid + ".json")
    ff = json.load(open(ff_p))
    have = {t["name"] for t in ff["atom_types"]}
    absent = sorted({r["type"] for r in rules if r["type"] not in have})
    doc = {"format": "caps-typing", "version": 1, "forcefield": label, "unknown_types": "untyped", "rules": rules,
           **({"bonds": "defined"} if fid.startswith("inorganic") and fid != "inorganic-zeolite-hill-sauer" else {}),
           "description": f"CAPS rules for {label} (bench/typing/make_inorganic_rules.py)." + (" " + extra if extra else "")
                          + (" Types this file lacks: " + ", ".join(absent) + "." if absent else "")}
    out = os.path.join(DATA, "typing", fid + ".typing.json")
    json.dump(doc, open(out, "w"), ensure_ascii=False, indent=1)
    ff["typing"] = f"../typing/{fid}.typing.json"
    json.dump(ff, open(ff_p, "w"), ensure_ascii=False, indent=1)
    return out, doc


if __name__ == "__main__":
    cat_p = os.path.join(DATA, "forcefields", "catalogue.json")
    cat = json.load(open(cat_p))
    made = {}
    for short, fid in FILES.items():
        rules, shells, chosen = compound_rules(short, fid)
        label = next(e["name"] for e in cat["forcefields"] if e["id"] == fid)
        out, doc = write(fid, rules, label, "Compound potential sets by composition (the first reference set of each; others by hand)."
                         + (" Shell model: a shell is added on each core." if shells else ""))
        if shells:
            doc["shells"] = shells
            json.dump(doc, open(out, "w"), ensure_ascii=False, indent=1)
        made[fid] = f"{len(rules)} rules, {len(chosen)} compounds" + (f", {len(shells)} shelled cores" if shells else "")
    for fid, rules, label in (("inorganic-clay", clayff(), "CLAYFF (metal-O-H bends need metal-oxygen bonds, which are not kept: left out)"), ("inorganic-zeolite-hill-sauer", hill_sauer(), "Hill-Sauer zeolite field"),
                              ("spce-moltemplate", water("O", "H"), "SPC/E water"), ("tip3p-1983-moltemplate", water("O", "H"), "TIP3P water (1983)"),
                              ("tip3p-2004-moltemplate", water("O", "H"), "TIP3P water (2004)")):
        write(fid, rules, label, "")
        made[fid] = f"{len(rules)} rules"
    for e in cat["forcefields"]:
        if e["id"] in made:
            e["typing"] = {"rules": f"typing/{e['id']}.typing.json", "evidence": f"CAPS rules (bench/typing/make_inorganic_rules.py): {made[e['id']]}"}
    json.dump(cat, open(cat_p, "w"), ensure_ascii=False, indent=1)
    for k, v in made.items():
        print(f"{k}: {v}")
