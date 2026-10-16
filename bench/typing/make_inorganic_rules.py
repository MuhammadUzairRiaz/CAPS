#!/usr/bin/env python3
"""CAPS typing rules for the library's inorganic force fields.

Compound-specific potential sets (binary halides and oxides, ternary oxides, glasses, zeolites, GaN) are typed from the
force field's own compound templates: each atom of a compound gets that compound's type when the structure holds the
compound's other elements (the rule's "requires"); where compounds nest (sodium silicate in soda-lime silicate) the rule
naming more elements wins. A composition with several potential sets takes the first reference set, except TiO2, which
takes the general set (rutile and anatase cannot be told apart by composition); the others are chosen by hand in the
Field step. Shell models declare their core → shell pairs ("shells"): CAPS puts a shell on every core it types.

The Hill-Sauer zeolite field (Hill & Sauer, J. Phys. Chem. 98, 1238 (1994); 99, 9536 (1995)) types by coordination and
neighbours. ClayFF is not here: its library and rules come from the paper
(bench/ff/make_clayff_cygan2004.py, data/typing/clayff-cygan2004.typing.json).

usage: make_inorganic_rules.py [SOURCE_LIB_DIR]
"""
import json, os, re, sys

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA = os.path.join(ROOT, "data")
SRC = sys.argv[1] if len(sys.argv) > 1 else os.path.expanduser("~/project/dl_f_4.13/lib")
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
    key = {}
    for l in at.splitlines():   # the first key of a type the table lists twice, as DL_FIELD takes it
        if len(l.split()) >= 3 and not l.startswith("#"):
            key.setdefault(l.split()[0], l.split()[1])
    keys = set(key.values())
    out = []
    for m in re.finditer(r"^MOLECULE (\S+).*?\n(.*?)^END MOLECULE", txt, re.S | re.M):
        labels = []
        for l in m.group(2).splitlines():
            w = l.split()
            # an atom line: label, DL_FIELD type (its key is the force-field type; the label need not be)
            k = key.get(w[1]) if len(w) >= 3 and w[1] in key else (w[0] if len(w) >= 3 and w[0] in keys else None)
            if k and k not in labels:
                labels.append(k)
        if labels:
            out.append((m.group(1), labels))
    return out


def set_rules(name, labels, el, types, springs, everything):
    """one compound template's rules: each core typed by composition, each shell by its name; and its core -> shell pairs"""
    rules, shells = [], {}
    # a core and its shell, as the template pairs them (Xc with Xs); the spring may be missing from the file (Fe in
    # alpha-Fe2O3 has none): CAPS then reports the missing core-shell term instead of dropping the shell
    paired = lambda c: c.endswith("c") and c[:-1] + "s" in labels
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
    for c, sh in shells.items():
        rules.append(rule(sh, 99, f"[{el[c]}]", f"shell of {c}", atom_name=sh))
    return rules, shells


def compound_rules(short, fid):
    ff = json.load(open(os.path.join(DATA, "forcefields", fid + ".json")))
    types = {t["name"]: t for t in ff["atom_types"]}
    springs = {tuple(sorted(b["match"])) for b in ff.get("bonds", [])}
    rules, shells, chosen, alternatives = [], {}, {}, {}
    everything = {types[l]["element"] for _, labels in templates(f"INORGANIC_{short}.sf") for l in labels if l in types and types[l].get("element")}
    # the file's water model, by composition (this file keeps only the bonds it defines, so water's O-H is not seen)
    for name, labels in templates(f"INORGANIC_{short}.sf"):
        if name.startswith("water") and all(l in types for l in labels):
            for l in labels:
                e = types[l]["element"]
                rules.append(rule(l, 2, f"[{e}]", f"{e} of water ({name.replace('_', ' ')})", requires=["H" if e == "O" else "O"],
                                  excludes=sorted(everything - {"O", "H"})))
            break
    for name, labels in templates(f"INORGANIC_{short}.sf"):
        if name.startswith("water"):
            continue
        if any(l not in types or "charge" not in types[l] for l in labels):
            continue   # a template naming a type without parameters or charges
        el = {l: types[l]["element"] for l in labels}
        comp = tuple(sorted(set(el.values())))
        pref = PREFER.get(fid, {}).get(comp)
        r, sh = set_rules(name, labels, el, types, springs, everything)
        if comp in chosen or (pref and name != pref):
            alternatives[name] = (r, sh)   # another reference set for this composition: an overlay the user picks
            continue
        chosen[comp] = name
        rules += r
        shells.update(sh)
    if fid == "inorganic-ternary-oxides":   # LaNiO3: its own oxygen core charge (O6cN), the shell the library's O6s
        rules.append(rule("O6cN", 5, "[O]", "O core in LaNiO3", requires=["La", "Ni"]))
        shells["O6cN"] = "O6s"
    if fid == "inorganic-zeolite":         # the acid faujasite's bridging hydroxyl: by its H, not by composition alone
        rules = [r for r in rules if r["type"] not in ("Ob2", "Hb2")]
        rules.append(rule("Ob2", 5, "[O;$(*[H])]", "bridging hydroxyl O of an acid zeolite", requires=["Al", "Si", "H"]))
        rules.append(rule("Hb2", 5, "[H][O]", "bridging hydroxyl H of an acid zeolite", requires=["Al", "Si"]))
    return rules, shells, chosen, alternatives


def hill_sauer():
    return [rule("Si", 0, "[Si]", "silicon"), rule("Al", 0, "[Al]", "aluminium"),
            rule("O_SS", 1, "[O;$(*(~[Si])~[Si])]", "Si-O-Si oxygen"), rule("O_AS", 1, "[O;H0;$(*(~[Si])~[Al])]", "Si-O-Al oxygen"),
            rule("Ob", 2, "[O;$(*(~[Si])~[Al]);$(*[H])]", "bridging hydroxyl oxygen Si-O(H)-Al"),
            rule("O_SH", 1, "[O;H1;$(*[Si]);!$(*[Al])]", "silanol oxygen"), rule("O_AH", 1, "[O;H1;$(*[Al]);!$(*[Si])]", "aluminol oxygen"),
            rule("Hb", 2, "[H][O;$(*~[Si]);$(*~[Al])]", "bridging hydroxyl hydrogen"),
            rule("H_OS", 1, "[H][O;$(*[Si]);!$(*[Al])]", "silanol hydrogen"), rule("H_OA", 1, "[H][O;$(*[Al]);!$(*[Si])]", "aluminol hydrogen")]


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
    have = lambda fid: os.path.exists(os.path.join(DATA, "forcefields", fid + ".json"))
    for short, fid in FILES.items():
        if not have(fid):   # the libraries converted from DL_FIELD's files are kept on disk only, never committed
            print(f"{fid}: no data/forcefields/{fid}.json here, skipped")
            continue
        rules, shells, chosen, alternatives = compound_rules(short, fid)
        # every other reference set as an overlay: caps ff type ... --typing <fid>.typing.json,sets/<fid>/<template>.typing.json
        sdir = os.path.join(DATA, "typing", "sets", fid)
        os.makedirs(sdir, exist_ok=True)
        for name, (r, sh) in alternatives.items():
            ov = {"format": "caps-typing", "version": 1, "forcefield": f"{fid}: {name}", "unknown_types": "untyped",
                  "description": f"Overlay on the {fid} rules: the {name.replace('_', ' ')} potential set in place of the composition's default "
                                 "(bench/typing/make_inorganic_rules.py).",
                  "rules": [dict(x, priority=x["priority"] + 50) if x["priority"] < 99 else x for x in r]}
            if sh:
                ov["shells"] = sh
            json.dump(ov, open(os.path.join(sdir, name + ".typing.json"), "w"), ensure_ascii=False, indent=1)
        label = next((e["name"] for e in cat["forcefields"] if e["id"] == fid), None) or \
            json.load(open(os.path.join(DATA, "forcefields", fid + ".json"))).get("name", fid)
        out, doc = write(fid, rules, label, "Compound potential sets by composition (the first reference set of each; others by hand)."
                         + (" Shell model: a shell is added on each core." if shells else ""))
        if shells:
            doc["shells"] = shells
            json.dump(doc, open(out, "w"), ensure_ascii=False, indent=1)
        made[fid] = f"{len(rules)} rules, {len(chosen)} compounds" + (f", {len(shells)} shelled cores" if shells else "") + \
                    (f", {len(alternatives)} other potential sets as overlays (typing/sets/{fid})" if alternatives else "")
    for fid, rules, label in (("inorganic-zeolite-hill-sauer", hill_sauer(), "Hill-Sauer zeolite field"),):
        if not have(fid):
            print(f"{fid}: no data/forcefields/{fid}.json here, skipped")
            continue
        write(fid, rules, label, "")
        made[fid] = f"{len(rules)} rules"
    for e in cat["forcefields"]:
        if e["id"] in made:
            e["typing"] = {"rules": f"typing/{e['id']}.typing.json", "evidence": f"CAPS rules (bench/typing/make_inorganic_rules.py): {made[e['id']]}"}
    json.dump(cat, open(cat_p, "w"), ensure_ascii=False, indent=1)
    for k, v in made.items():
        print(f"{k}: {v}")
