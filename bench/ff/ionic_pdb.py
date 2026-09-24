#!/usr/bin/env python3
"""Validation structures for DL_FIELD's ionic templates (every ion CONNECTs only to itself): formula units on a
simple cubic grid, ions ≥ 4.5 Å apart (no bonds), residue = the template's MOLECULE_KEY. Core-shell pairs
(names ending c / s) are put 0.1 Å apart.

usage: ionic_pdb.py LIB.sf MOLECULE OUT.pdb [--units 8]
       ionic_pdb.py LIB.sf --list      (ionic templates)
"""
import re, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import template_pdb

def ionic(m):
    labels = [a for a, _ in m["atoms"]]
    return m["atoms"] and all(js == [a] for a, js in m["connect"].items()) and set(m["connect"]) >= set(labels)

def build(sf, name, out, units=8, gap=4.5):
    types, keys, mols = template_pdb.parse_sf(sf)
    m = mols[name]
    key = keys.get(name, name[:4].upper())
    n = round(units ** (1 / 3))
    lines, k = [], 0
    for u in range(units):
        ox, oy, oz = (u % n) * gap * 3, (u // n % n) * gap * 3, (u // (n * n)) * gap * 3
        for q, (lab, ty) in enumerate(m["atoms"]):
            el = types.get(ty, re.sub(r"[^A-Za-z]", "", lab)[:2])
            shell = lab.endswith("s") and any(l == lab[:-1] + "c" for l, _ in m["atoms"])
            base = [l for l, _ in m["atoms"]].index(lab[:-1] + "c") if shell else q
            x = [ox + (base % 2) * gap, oy + (base // 2 % 2) * gap, oz + (base // 4) * gap]
            if shell: x[0] += 0.1
            k += 1
            nm = lab if len(lab) == 4 else f" {lab:<3s}"
            lines.append(f"HETATM{k:5d} {nm:4s} {key:<4s}{u + 1:5d}    {x[0]:8.3f}{x[1]:8.3f}{x[2]:8.3f}  1.00  0.00          {el[:2]:>2s}")
    open(out, "w").write("REMARK   CAPS ionic validation structure\n" + "\n".join(lines) + "\nEND\n")

if __name__ == "__main__":
    if "--list" in sys.argv:
        types, keys, mols = template_pdb.parse_sf(sys.argv[1])
        print("\n".join(n for n, m in mols.items() if n in keys and ionic(m)))
    else:
        units = int(sys.argv[sys.argv.index("--units") + 1]) if "--units" in sys.argv else 8
        build(sys.argv[1], sys.argv[2], sys.argv[3], units)
