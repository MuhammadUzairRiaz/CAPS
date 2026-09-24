#!/usr/bin/env python3
"""Validation molecules for COMPASS (moltemplate compass_published.lt): a siloxane, an ester-ether and an
alkyl benzene. usage: compass_molecules.py OUTDIR  (writes NAME.lt for each; run moltemplate on them)."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from zmat import build, methyl, write_lt

def siloxane():
    # Me3Si-O-SiMe2-O-SiMe3 (octamethyltrisiloxane)
    r = [("Si1", "si4c", 0, None, 0, None, 0, None, 0),
         ("O1", "o2z", 0, "Si1", 1.64, None, 0, None, 0),
         ("Si2", "si4c", 0, "O1", 1.64, "Si1", 145, None, 0),
         ("O2", "o2z", 0, "Si2", 1.64, "O1", 109.5, "Si1", 180),
         ("Si3", "si4c", 0, "O2", 1.64, "Si2", 145, "O1", 170)]
    for si, x, y, n in (("Si1", "O1", "Si2", 3), ("Si2", "O1", "Si1", 2), ("Si3", "O2", "Si2", 3)):
        dih = [60, 180, 300] if n == 3 else [60, 300]
        for k in range(n):
            c = f"C{si[2]}{k+1}"
            r.append((c, "c4", 0, si, 1.87, x, 110.0, y, dih[k] + (0 if si != "Si2" else 0)))
            r += methyl(c, si, x)
    return r

def ester_ether():
    # CH3-CH2-C(=O)-O-CH2-CH2-O-CH3  (2-methoxyethyl propanoate)
    r = [("C1", "c4", 0, None, 0, None, 0, None, 0),
         ("C2", "c4", 0, "C1", 1.53, None, 0, None, 0),
         ("C3", "c3prime", 0, "C2", 1.51, "C1", 112, None, 0),
         ("O1", "o1=", 0, "C3", 1.21, "C2", 124, "C1", 0),
         ("O2", "o2s", 0, "C3", 1.35, "C2", 111, "C1", 180),
         ("C4", "c4", 0, "O2", 1.44, "C3", 116, "C2", 180),
         ("C5", "c4", 0, "C4", 1.52, "O2", 108, "C3", 175),
         ("O3", "o2e", 0, "C5", 1.42, "C4", 109, "O2", 65),
         ("C6", "c4", 0, "O3", 1.42, "C5", 112, "C4", 178)]
    r += methyl("C1", "C2", "C3")
    r += [("H21", "h1", 0, "C2", 1.09, "C1", 110, "C3", 120), ("H22", "h1", 0, "C2", 1.09, "C1", 110, "C3", -120),
          ("H41", "h1", 0, "C4", 1.09, "O2", 109, "C5", 120), ("H42", "h1", 0, "C4", 1.09, "O2", 109, "C5", -120),
          ("H51", "h1", 0, "C5", 1.09, "C4", 110, "O3", 120), ("H52", "h1", 0, "C5", 1.09, "C4", 110, "O3", -120)]
    r += methyl("C6", "O3", "C5")
    return r

def propylbenzene():
    r = [("A1", "c3a", 0, None, 0, None, 0, None, 0),
         ("A2", "c3a", 0, "A1", 1.39, None, 0, None, 0),
         ("A3", "c3a", 0, "A2", 1.39, "A1", 120, None, 0),
         ("A4", "c3a", 0, "A3", 1.39, "A2", 120, "A1", 0),
         ("A5", "c3a", 0, "A4", 1.39, "A3", 120, "A2", 0),
         ("A6", "c3a", 0, "A5", 1.39, "A4", 120, "A3", 0)]
    r += [(f"HA{k}", "h1", 0, f"A{k}", 1.08, f"A{k-1}", 120, f"A{k-2}", 180) for k in (3, 4, 5)]
    r += [("HA6", "h1", 0, "A6", 1.08, "A5", 120, "A4", 180), ("HA2", "h1", 0, "A2", 1.08, "A3", 120, "A4", 180)]
    r += [("C1", "c4", 0, "A1", 1.51, "A2", 120, "A3", 180),
          ("C2", "c4", 0, "C1", 1.53, "A1", 113, "A2", 90),
          ("C3", "c4", 0, "C2", 1.53, "C1", 112, "A1", 178),
          ("H11", "h1", 0, "C1", 1.09, "A1", 109, "C2", 120), ("H12", "h1", 0, "C1", 1.09, "A1", 109, "C2", -120),
          ("H21", "h1", 0, "C2", 1.09, "C1", 109, "C3", 120), ("H22", "h1", 0, "C2", 1.09, "C1", 109, "C3", -120)]
    r += methyl("C3", "C2", "C1")
    return r, [("A6", "A1")]

if __name__ == "__main__":
    out = sys.argv[1]
    os.makedirs(out, exist_ok=True)
    write_lt(build(siloxane()), "compass_published", "COMPASS", "Siloxane", os.path.join(out, "siloxane.lt"))
    write_lt(build(ester_ether()), "compass_published", "COMPASS", "EsterEther", os.path.join(out, "ester.lt"))
    rows, extra = propylbenzene()
    write_lt(build(rows), "compass_published", "COMPASS", "PropylBenzene", os.path.join(out, "pbz.lt"), extra_bonds=extra)
