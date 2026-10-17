#!/usr/bin/env python3
"""MARTINI 2.0 and Dry Martini for CAPS from their papers (PDFs in ~/project/ff-sources/papers/cg, never in the
repository; the cgmartini.nl .itp files are not used):

  martini2-marrink2007      S. J. Marrink, H. J. Risselada, S. Yefimov, D. P. Tieleman, A. H. de Vries, J. Phys. Chem. B 111,
                            7812 (2007): Table 1 (the interaction matrix of the 18 bead types), section 2.2 (the levels'
                            epsilon, sigma 0.47 nm, 0.62 nm for level IX, shifted LJ and Coulomb, epsilon_r 15), 2.3 (bonds
                            0.47 nm / 1250 kJ mol-1 nm-2, cosine angles 25 kJ mol-1 at 180 deg, cis 45 at 120, trans 45 at 180;
                            bonded neighbours excluded, second neighbours not), 2.4 (ring "S" beads: sigma 0.43 nm and 75 % of
                            epsilon between ring beads), 2.5 (antifreeze BP4: BP4-P4 at level O with sigma 0.57 nm, otherwise
                            P4), 2.6 (cut-off 1.2 nm, LJ shifted from 0.9 nm, Coulomb from 0), 2.7 (masses 72 amu), 2.8
                            (topologies: water P4; dodecane and hexadecane three and four C1; octanol C1-P1; chloroform C4;
                            ether N0; cyclohexane three SC1, 0.3 nm, 5000 kJ mol-1 nm-2; PC head Q0 + Qa, PE head Qd + Qa,
                            glycerol two Na 0.37 nm apart, tails of four C1 per 16 carbons, oleoyl C1 C1 C3 C1 C1; ions Q
                            beads with full charges)
  drymartini-arnarez2015    C. Arnarez et al., J. Chem. Theory Comput. 11, 260 (2015): Table 1 (the implicit-solvent
                            interaction matrix and levels X-XII, epsilon 1.5 / 1.0 / 0.5 kJ/mol; Q-Q sigma 0.60 nm, level IX
                            0.62 nm, ring beads as in MARTINI), section 3.1 (head bond 0.45 nm, tail bonds 0.48 nm, saturated
                            tail angles 35 kJ mol-1)

Single-atom ions: the paper models them as Q beads carrying a full charge; CAPS takes Qd for the cation and Qa for the
anion (as the paper's PE and PC head groups use them for donor and acceptor charges). Units: kcal/mol and A (LAMMPS
real; harmonic K = k/2, cosine/squared K = k/2 from GROMACS's (1/2) k forms).

usage: make_martini2.py [DATA_DIR]
"""
import json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.dirname(HERE)), "data")
KJ = 4.184
TYPES = ["Qda", "Qd", "Qa", "Q0", "P5", "P4", "P3", "P2", "P1", "Nda", "Nd", "Na", "N0", "C5", "C4", "C3", "C2", "C1"]
LEVEL = {"O": 5.6, "I": 5.0, "II": 4.5, "III": 4.0, "IV": 3.5, "V": 3.1, "VI": 2.7, "VII": 2.3, "VIII": 2.0, "IX": 2.0,
         "X": 1.5, "XI": 1.0, "XII": 0.5}
MARRINK = "S. J. Marrink, H. J. Risselada, S. Yefimov, D. P. Tieleman, A. H. de Vries, J. Phys. Chem. B 111, 7812 (2007)"
ARNAREZ = ("C. Arnarez, J. J. Uusitalo, M. F. Masman, H. I. Ingolfsson, D. H. de Jong, M. N. Melo, X. Periole, A. H. de Vries, "
           "S. J. Marrink, J. Chem. Theory Comput. 11, 260 (2015)")
# Marrink 2007, Table 1, row by row in TYPES order (read from the page)
M2 = """O O O II O O O I I I I I IV V VI VII IX IX
O I O II O O O I I I III I IV V VI VII IX IX
O O I II O O O I I I I III IV V VI VII IX IX
II II II IV I O I II III III III III IV V VI VII IX IX
O O O I O O O O O I I I IV V VI VI VII VIII
O O O O O I I II II III III III IV V VI VI VII VIII
O O O I O I I II II II II II IV IV V V VI VII
I I I II O II II II II II II II III IV IV V VI VII
I I I III O II II II II II II II III IV IV IV V VI
I I I III I III II II II II II II IV IV V VI VI VI
I III I III I III II II II II III II IV IV V VI VI VI
I I III III I III II II II II II III IV IV V VI VI VI
IV IV IV IV IV IV IV III III IV IV IV IV IV IV IV V VI
V V V V V V IV IV IV IV IV IV IV IV IV IV V V
VI VI VI VI VI VI V IV IV V V V IV IV IV IV V V
VII VII VII VII VI VI V V IV VI VI VI IV IV IV IV IV IV
IX IX IX IX VII VII VI VI V VI VI VI V V V IV IV IV
IX IX IX IX VIII VIII VII VII VI VI VI VI VI V V IV IV IV"""
# Arnarez 2015, Table 1 (lower triangle, rows in TYPES order)
DRY = """VI
VI VII
VI VI VII
VIII VIII VIII VIII
XII XII XII XII V
XII XII XII XII V VI
XII XII XII XII V VI VI
XII XII XII XII V VII VII VII
XII XII XII XII X X X VIII VII
XII XII XII XII XII XI XI XI VII VI
XII XII XII XII XII XI XI XI VII VI VII
XII XII XII XII XII XI XI XI VII VI VI VII
XII XII XII XII XII XI X IX VII VII VII VII V
XII XII XII XII XII XI X VI VI VI VI VI V III
XII XII XII XII XII XI X VI VI VI VI VI IV III II
XII XII XII XII XII XI X VII VI VI VI VI IV III II II
IX IX IX IX XII XII X VIII VI VI VI VI IV III III II II
IX IX IX IX XII XII XI X VII VI VI VI IV III III II II II"""


def matrix_full(text):
    rows = [r.split() for r in text.strip().splitlines()]
    n = len(TYPES)
    m = [[None] * n for _ in range(n)]
    for i, r in enumerate(rows):
        for j, lv in enumerate(r):
            m[i][j] = lv
            if m[j][i] is None:
                m[j][i] = lv
    for i in range(n):
        for j in range(n):
            assert m[i][j] is not None and m[i][j] == m[j][i], (TYPES[i], TYPES[j], m[i][j], m[j][i])
    return m


def bond(match, r_nm, k_kj, why):
    return {"name": "-".join(match), "match": match, "style": "harmonic", "params": [k_kj / KJ / 100 / 2, r_nm * 10], "comment": why}


def angle(match, th, k_kj, why):
    return {"name": "-".join(match), "match": match, "style": "cosine/squared", "params": [k_kj / KJ / 2, th], "comment": why}


# templates: name -> bead SMILES (2007 section 2.8; ions as noted above)
TAIL16 = "[C1][C1][C1][C1]"
OLEOYL = "[C1][C1][C3][C1][C1]"
TEMPLATES = {
    "W": "[P4]", "WF": "[BP4]", "NA+": "[Qd+1]", "CL-": "[Qa-1]",
    "DOD": "[C1][C1][C1]", "HD": "[C1][C1][C1][C1]", "OCO": "[C1][P1]", "CLF": "[C4]", "ETH": "[N0]", "CHEX": "[SC1]1[SC1][SC1]1",
    "DPPC": f"[Q0+1][Qa-1][Na]({TAIL16})[Na]{TAIL16}", "DLPC": "[Q0+1][Qa-1][Na]([C1][C1][C1])[Na][C1][C1][C1]",
    "DCPC": "[Q0+1][Qa-1][Na]([C1][C1])[Na][C1][C1]", "DOPC": f"[Q0+1][Qa-1][Na]({OLEOYL})[Na]{OLEOYL}",
    "POPC": f"[Q0+1][Qa-1][Na]({OLEOYL})[Na]{TAIL16}",
    "DPPE": f"[Qd+1][Qa-1][Na]({TAIL16})[Na]{TAIL16}", "DOPE": f"[Qd+1][Qa-1][Na]({OLEOYL})[Na]{OLEOYL}",
    "POPE": f"[Qd+1][Qa-1][Na]({OLEOYL})[Na]{TAIL16}",
}


def library(fid, name, ref, source, levels, dry):
    m = matrix_full(levels)
    ring = ["S" + t for t in TYPES]
    names = TYPES + ring + ([] if dry else ["BP4"])
    base = lambda t: t[1:] if t.startswith("S") and t[1:] in TYPES else ("P4" if t == "BP4" else t)
    types = [{"name": t, "element": "", "mass": 72.0, "charge": 0.0,
              "description": ("ring bead " if t.startswith("S") and t[1:] in TYPES else "antifreeze bead " if t == "BP4" else "bead ") + t}
             for t in names]
    pairs = []
    for a in range(len(names)):
        for b in range(a, len(names)):
            ta, tb = names[a], names[b]
            lv = m[TYPES.index(base(ta))][TYPES.index(base(tb))]
            eps, sig = LEVEL[lv], 0.47
            both_ring = ta in ring and tb in ring
            if lv == "IX":
                sig = 0.62
            elif dry and base(ta).startswith("Q") and base(tb).startswith("Q"):
                sig = 0.60
            elif both_ring:
                sig = 0.43
            if both_ring:
                eps *= 0.75
            why = f"level {lv}" + (", ring pair: sigma 0.43 nm, 75 % epsilon" if both_ring and lv != "IX" else ", ring pair: 75 % epsilon" if both_ring else "")
            if {ta, tb} == {"BP4", "P4"}:
                eps, sig, why = LEVEL["O"], 0.57, "antifreeze BP4-P4: level O, sigma 0.57 nm (section 2.5)"
            pairs.append({"name": f"{ta}-{tb}", "match": [ta, tb] if ta != tb else [ta], "params": [eps / KJ, sig * 10],
                          "comment": f"{why} ({eps:.4f} kJ/mol, {sig} nm)"})
    if dry:
        bonds = [bond(["*", "*"], 0.47, 1250, "standard bond (Marrink 2007, section 2.3)"),
                 bond(["Q?", "Qa"], 0.45, 1250, "head-group bond 0.45 nm (section 3.1)"),
                 bond(["Na", "Na"], 0.37, 1250, "glycerol backbone 0.37 nm (Marrink 2007, section 2.8.3)"),
                 bond(["C?", "C?"], 0.48, 1250, "tail bonds 0.48 nm (section 3.1)"),
                 bond(["SC1", "SC1"], 0.30, 5000, "cyclohexane ring (Marrink 2007, section 2.8.1)")]
        angles = [angle(["*", "*", "*"], 180.0, 25, "standard angle (Marrink 2007, section 2.3)"),
                  angle(["*", "C?", "*"], 180.0, 35, "saturated tails 35 kJ/mol (section 3.1)"),
                  angle(["C1", "C3", "C1"], 120.0, 45, "cis double bond (Marrink 2007, section 2.3)"),
                  angle(["SC1", "SC1", "SC1"], 60.0, 0, "the cyclohexane triangle: no angle term (its beads are all bonded)")]
    else:
        bonds = [bond(["*", "*"], 0.47, 1250, "standard bond (section 2.3)"),
                 bond(["Na", "Na"], 0.37, 1250, "glycerol backbone 0.37 nm (section 2.8.3)"),
                 bond(["SC1", "SC1"], 0.30, 5000, "cyclohexane ring (section 2.8.1)")]
        angles = [angle(["*", "*", "*"], 180.0, 25, "standard angle (section 2.3)"),
                  angle(["C1", "C3", "C1"], 120.0, 45, "cis double bond (section 2.3)"),
                  angle(["SC1", "SC1", "SC1"], 60.0, 0, "the cyclohexane triangle: no angle term (its beads are all bonded)")]
    doc = {"format": "caps-forcefield", "format_version": 1, "name": name, "source": source, "references": [ref] if not dry else [ref, MARRINK],
           "units": "real",
           "styles": {"pair": "lj/gromacs/coul/gromacs", "bond": "harmonic", "angle": "cosine/squared", "dihedral": "charmm", "improper": "cvff"},
           "mixing": "arithmetic", "special_lj": [0, 1, 1], "special_coul": [0, 1, 1], "cutoff": 12.0, "timestep": 20,
           "pair_settings": {"lj_inner": 9.0, "coul_inner": 1e-06, "dielectric": 15.0, "model_cutoff": True},
           "angle_terms": "if_defined", "torsion_terms": "if_defined",
           "atom_types": types, "pairs": pairs, "bonds": bonds, "angles": angles, "dihedrals": [], "impropers": [],
           "bead_templates": {k: v for k, v in TEMPLATES.items() if not (dry and k in ("W", "WF"))},
           "notes": [f"built by CAPS (bench/ff/make_martini2.py) from {source}",
                     "shifted LJ (from 0.9 nm) and Coulomb (from 0, epsilon_r 15) to 1.2 nm, as GROMACS's shift functions; masses 72 amu",
                     "single-atom ions: Qd (cation) and Qa (anion) with full charges (CAPS's choice of Q subtype; the paper says Q beads)"],
           "version": "2007" if not dry else "2015", "status": "validated", "typing": f"../typing/{fid}.typing.json"}
    typing = {"format": "caps-typing", "version": 1, "forcefield": name, "coarse_grained": True, "unknown_types": "untyped",
              "description": f"CAPS rules for {name} (bench/ff/make_martini2.py): beads typed by name (bead templates or files that name "
                             "their sites by bead type)",
              "rules": [{"type": t, "smarts": "*", "atom_name": t, "priority": 0, "description": f"a bead named {t}"} for t in names]}
    return doc, typing


def main():
    libs = [("martini2-marrink2007", "MARTINI 2.0", MARRINK, MARRINK + ", Table 1 and sections 2.2-2.8", M2, False,
             ["martini", "martini2"], {"label": "MARTINI 2", "key": "martini2", "order": 37}),
            ("drymartini-arnarez2015", "Dry Martini (implicit solvent)", ARNAREZ, ARNAREZ + ", Table 1 and section 3.1; MARTINI 2.0's bonded "
             "terms otherwise", DRY, True, ["drymartini", "dry-martini"], {"label": "Dry Martini", "key": "drymartini", "order": 38})]
    cat_p = os.path.join(DATA, "forcefields", "catalogue.json")
    cat = json.load(open(cat_p))
    for fid, name, ref, source, levels, dry, aliases, listing in libs:
        doc, typing = library(fid, name, ref, source, levels, dry)
        json.dump(doc, open(os.path.join(DATA, "forcefields", fid + ".json"), "w"), indent=1, ensure_ascii=False)
        json.dump(typing, open(os.path.join(DATA, "typing", fid + ".typing.json"), "w"), indent=1, ensure_ascii=False)
        e = {"id": fid, "name": name, "version": doc["version"], "references": doc["references"], "origin": "built from source",
             "source": source + "; bench/ff/make_martini2.py", "license": "published parameters of the paper(s)", "status": "validated",
             "notes": f"{len(doc['atom_types'])} bead types, every pair from the paper's interaction matrix; templates "
                      + ", ".join(doc["bead_templates"]), "file": fid + ".json",
             "typing": {"rules": f"typing/{fid}.typing.json", "evidence": "beads typed by name"},
             "counts": {"atom_types": len(doc["atom_types"]), **{k: len(doc.get(k, [])) for k in ("pairs", "bonds", "angles", "dihedrals", "impropers")}},
             "group": "cg", "aliases": aliases, "list": listing}
        cat["forcefields"] = [x for x in cat["forcefields"] if x["id"] != fid] + [e]
        print(fid, "written:", len(doc["pairs"]), "pairs")
    json.dump(cat, open(cat_p, "w"), indent=1, ensure_ascii=False)


if __name__ == "__main__":
    main()
