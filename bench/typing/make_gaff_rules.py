#!/usr/bin/env python3
"""Write CAPS typing rules for GAFF and GAFF2 (data/typing/gaff*.typing.json).

The rules follow antechamber's definition files (ATOMTYPE_GFF.DEF for GAFF, ATOMTYPE_GFF2.DEF for GAFF2 in
AmberTools) line by line and in the same order: antechamber takes the first definition that matches, and so does
CAPS among rules of equal priority. antechamber's ring classes AR1..AR5 are the CAPS SMARTS extension {AR1}..{AR5};
its conjugated pairs (cc/cd, ce/cf, ...) are the rules file's "pairs". Water keeps GAFF's own ow / hw.

usage: make_gaff_rules.py OUTDIR
"""
import json, os, sys

EW = "[#7,#8,#9,#16,#17,#35,#53]"          # antechamber's electron-withdrawing atoms (ring.c: N O S F Cl Br I), for h1-h5
RING = "{AR1},{AR2},{AR3}"                  # planar rings
CONJ_SB = "$(*-[#6X3,#6X2,#7X2,#15X2]),$(*-[#16X3,#16X4,#15X3,#15X4;$(*=*)])"        # ce: C3 C2 XB2 XD3 XD4 by a pure single bond
CONJ_SB_XA = "$(*-[#6X3,#6X2,#8X1,#16X1,#7X2,#15X2]),$(*-[#16X3,#16X4,#15X3,#15X4;$(*=*)])"   # ne / pe also XA1
CONJ_sb = "$(*-,:[#7X2,#15X2,#6X3]),$(*-,:[#16X3,#16X4,#15X3,#15X4;$(*=*)])"         # px / py / sx / sy: single incl. aromatic


def has(k):
    """at least k electron-withdrawing neighbours"""
    return "$(*" + "".join("(" + EW + ")" for _ in range(k - 1)) + EW + ")"


def ew(n):
    """a hydrogen's carbon with exactly n electron-withdrawing neighbours"""
    return has(n) + ";!" + has(n + 1)


def rules(v2):
    R = []
    add = lambda t, *sm, d="": [R.append({"type": t, "smarts": x, **({"description": d} if d else {})}) for x in sm]
    add("ow", "[OX2H2]", d="water O (GAFF / TIP3P)")
    add("hw", "[H][OX2H2]", d="water H")
    # carbon
    add("cx", "[#6X4;r3]", d="sp3 C in a 3-ring")
    add("cy", "[#6X4;r4]", d="sp3 C in a 4-ring")
    if v2:
        add("c5", "[#6X4;r5]", d="sp3 C in a 5-ring")
        add("c6", "[#6X4;r6]", d="sp3 C in a 6-ring")
    add("c3", "[#6X4]", d="sp3 C")
    if v2:
        add("cs", "[#6X3]~[#16X1]", d="C=S")
    add("c", "[#6X3]~[#8X1]" if v2 else "[#6X3]~[#8X1,#16X1]", d="C=O" + ("" if v2 else " / C=S"))
    add("cz", "[#6X3](~[#7X3])(~[#7X3])~[#7X3]", d="guanidinium C")
    add("cp", "[#6X3;{AR1};R1](-!@[*;{AR1}])", d="aromatic C at a single bond between aromatic rings (biphenyl)")
    add("ca", "[#6X3;{AR1}]", d="pure aromatic C")
    add("cc", "[#6X3;{AR2},{AR3}]", d="sp2 C in a conjugated ring (with cd)")
    add("ce", f"[#6X3;$(*=*);{CONJ_SB}]", d="sp2 C in a conjugated chain (with cf)")
    add("cu", "[#6X3;r3]", d="sp2 C in a 3-ring")
    add("cv", "[#6X3;r4]", d="sp2 C in a 4-ring")
    add("c2", "[#6X3]", d="sp2 C")
    add("cg", "[#6X2;$(*#*);$(*-[#6X2,#6X3,#7X1,#7X2,#15X2])]", d="sp C conjugated (with ch)")
    add("c1", "[#6X2]", "[#6X1]", d="sp C")
    # hydrogen
    add("hn", "[H][#7]"); add("ho", "[H][#8]"); add("hs", "[H][#16]"); add("hp", "[H][#15]")
    add("hx", "[H][#6][#7X4]", d="H on a C next to N+")
    add("h3", f"[H][#6X4;{ew(3)}]", d="H on sp3 C with 3 electron-withdrawing neighbours")
    add("h2", f"[H][#6X4;{ew(2)}]")
    add("h1", f"[H][#6X4;{ew(1)}]")
    add("hc", "[H][#6X4]", d="H on sp3 C")
    add("h5", f"[H][#6X3;{ew(2)}]", d="H on sp2 C with 2 electron-withdrawing neighbours")
    add("h4", f"[H][#6X3;{ew(1)}]")
    add("ha", "[H]", d="H on sp2 / aromatic C")
    add("f", "[F]"); add("cl", "[Cl]"); add("br", "[Br]"); add("i", "[I]")
    # phosphorus
    add("pb", "[#15X2;{AR1}]")
    add("pc", "[#15X2;{AR2},{AR3}]")
    add("pe", f"[#15X2;$(*=*);{CONJ_SB_XA}]")
    add("p2", "[#15X2]", "[#15X1]")
    add("px", f"[#15X3;$(*=*);{CONJ_sb}]")
    add("p4", "[#15X3;$(*=*);$(*~[#8X1,#16X1])]")
    add("p3", "[#15X3]")
    add("py", f"[#15X4;$(*=*);{CONJ_sb}]")
    add("p5", "[#15X4]", "[#15X5]", "[#15X6]")
    # nitrogen
    amide = "$(*[#6X3]~[#8X1,#16X1])"
    nh_env = [f"[#7X3;$(*[#6,#7,#8,#16,#15;{RING}])", "[#7X3;$(*[#6X3;$(*=*)])", "[#7X3;$(*[#7X2;$(*=*)])", "[#7X3;$(*[#15X2;$(*=*)])"]
    if v2:
        add("ns", f"[#7X3H1;{amide}]", d="amide N with 1 H")
        add("nt", f"[#7X3H2;{amide}]", d="amide N with 2 H")
        add("ni", f"[#7X3;r3;{amide}]"); add("nj", f"[#7X3;r4;{amide}]")
    add("n", f"[#7X3;{amide}]", d="amide N")
    if v2:
        add("nk", "[#7X4;r3]"); add("nl", "[#7X4;r4]")
        add("nx", "[#7X4H1]"); add("ny", "[#7X4H2]"); add("nz", "[#7X4H3]"); add("n+", "[#7X4H4]")
    add("n4", "[#7X4]", d="sp3 N with 4 connections")
    add("no", "[#7X3](~[#8X1])~[#8X1]", d="nitro N")
    add("na", f"[#7X3;{RING}]", d="N with 3 connections in a planar ring")
    if v2:
        add("nu", *[x + ";H1]" for x in nh_env], d="amine N on an aromatic / conjugated atom, 1 H")
        add("nv", *[x + ";H2]" for x in nh_env], d="same, 2 H")
        add("nm", *[x + ";r3]" for x in nh_env]); add("nn", *[x + ";r4]" for x in nh_env])
    add("nh", *[x + "]" for x in nh_env], d="amine N on an aromatic / conjugated atom")
    if v2:
        add("np", "[#7X3H0;r3]"); add("nq", "[#7X3H0;r4]"); add("n5", "[#7X3H1;r3]"); add("n6", "[#7X3H1;r4]")
        add("n7", "[#7X3H1]", d="sp3 amine N with 1 H"); add("n8", "[#7X3H2]", d="with 2 H"); add("n9", "[#7X3H3]", d="ammonia")
    add("n3", "[#7X3]", d="sp3 N")
    add("nb", "[#7X2;{AR1}]", d="pure aromatic N")
    add("nc", "[#7X2;{AR2},{AR3}]", d="sp2 N in a conjugated ring (with nd)")
    add("ne", f"[#7X2;$(*=*);{CONJ_SB_XA}]", d="sp2 N in a conjugated chain (with nf)")
    add("n1", "[#7X2;$(*(=*)=*)]", "[#7X2;$(*#*)]", "[#7X1]", d="sp N")
    add("n2", "[#7X2]", d="sp2 N")
    # oxygen
    add("o", "[#8X1]", d="O with one connection")
    add("oh", "[#8X2;H1,H2]", "[#8X3;H1,H2,H3]", d="hydroxyl O")
    if v2:
        add("op", "[#8X2;r3]"); add("oq", "[#8X2;r4]")
    add("os", "[#8X2]", "[#8X3]", "[#8]", d="ether / ester O")
    # sulfur
    add("s", "[#16X1]")
    add("s2", "[#16X2;$(*=*),$(*#*)]")
    add("sh", "[#16X2;H1,H2]")
    if v2:
        add("sp", "[#16X2;r3]"); add("sq", "[#16X2;r4]")
    add("ss", "[#16X2]")
    add("sx", f"[#16X3;$(*=*);{CONJ_sb}]")
    add("s4", "[#16X3]")
    add("sy", f"[#16X4;$(*=*);{CONJ_sb}]")
    add("s6", "[#16X4]", "[#16X5]", "[#16X6]")
    return R


PAIRS = [["cc", "cd"], ["ce", "cf"], ["cp", "cq"], ["nc", "nd"], ["ne", "nf"], ["pc", "pd"], ["pe", "pf"]]
out = sys.argv[1]
for name, v2, ff, src in [("gaff-amber16", False, "GAFF (AmberTools 16, DL_FIELD)", "ATOMTYPE_GFF.DEF"),
                          ("gaff-amber25", True, "GAFF2 (AmberTools 25, DL_FIELD)", "ATOMTYPE_GFF2.DEF")]:
    json.dump({"format": "caps-typing", "version": 1, "forcefield": ff,
               "description": f"GAFF{'2' if v2 else ''} atom types following antechamber's {src} in order (first matching rule wins); "
                              "{AR1}..{AR5} are antechamber's ring classes, pairs its conjugation pattern.",
               "ordered": True, "pairs": PAIRS, "rules": rules(v2)},
              open(os.path.join(out, name + ".typing.json"), "w"), indent=1)
    print("wrote", name)
