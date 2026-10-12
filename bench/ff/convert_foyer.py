#!/usr/bin/env python3
"""A CAPS force-field file from a foyer force field (OpenMM XML; mosdef-hub/foyer, MIT): foyer's trappe-ua.xml.

  convert_foyer.py FILE.xml OUT.json --name NAME [--reference "..."] [--version "..."] [--typing ../typing/ID.typing.json]
                   [--cutoff 14] [--tail]

OpenMM's units and forms (kJ/mol, nm, rad) to CAPS's (real units):
  atoms      each Type's class names it in the bonded terms; NonbondedForce charge, sigma, epsilon (epsilon / 4.184, sigma x 10)
  bonds      HarmonicBondForce 1/2 k (r - r0)^2: harmonic K = k / 2 / 4.184 / 100, r0 x 10
  angles     HarmonicAngleForce 1/2 k (theta - theta0)^2: harmonic K = k / 2 / 4.184, theta0 in degrees
  torsions   RBTorsionForce sum c_n cos^n(psi), psi = phi - 180 deg: written exactly as fourier terms of phi,
             a_n = (-1)^n c_n, E = a0 + a2/2 + (a1 + 3 a3 / 4) cos phi + a2 / 2 cos 2phi + a3 / 4 cos 3phi (c4 = c5 = 0),
             each b_n cos n phi as b_n (1 + cos n phi) and the constant left over as an n = 0 term (K (1 + cos 90 deg) = K)
  1-4        NonbondedForce coulomb14scale / lj14scale; combining_rule lorentz: arithmetic mixing
Atom types are typed by CAPS's own rules (the typing file), not foyer's SMARTS.
"""
import json, math, sys
import xml.etree.ElementTree as ET

src = sys.argv[1]
OUT = sys.argv[2]
arg = lambda k, d=None: sys.argv[sys.argv.index(k) + 1] if k in sys.argv else d
many = lambda k: [sys.argv[i + 1] for i, a in enumerate(sys.argv) if a == k]
KJ = 4.184
root = ET.parse(src).getroot()
cls = {}
types, pairs = [], []
el_of = {"_CH4": "C", "_CH3": "C", "_CH2": "C", "_HC": "C", "_CH": "C"}
for t in root.find("AtomTypes"):
    cls[t.get("name")] = t.get("class")
nb = root.find("NonbondedForce")
q = {a.get("type"): a for a in nb}
for t in root.find("AtomTypes"):
    n, c = t.get("name"), t.get("class")
    a = q[n]
    e = t.get("element")
    types.append({"name": n, "element": el_of.get(e, e), "mass": float(t.get("mass")), "charge": float(a.get("charge")),
                  "equivalence": {k: c for k in ("bond", "angle", "dihedral", "improper")},
                  "description": f"{t.get('desc')} (doi:{t.get('doi')})"})
    pairs.append({"name": n, "match": [n], "params": [float(a.get("epsilon")) / KJ, float(a.get("sigma")) * 10]})
bonds = [{"name": f"{b.get('class1')}-{b.get('class2')}", "match": [b.get("class1"), b.get("class2")], "style": "harmonic",
          "params": [float(b.get("k")) / 2 / KJ / 100, float(b.get("length")) * 10]} for b in root.find("HarmonicBondForce")]
angles = [{"name": "-".join(a.get(f"class{i}") for i in (1, 2, 3)), "match": [a.get(f"class{i}") for i in (1, 2, 3)], "style": "harmonic",
           "params": [float(a.get("k")) / 2 / KJ, math.degrees(float(a.get("angle")))]} for a in root.find("HarmonicAngleForce")]
dihedrals = []
for p in root.find("RBTorsionForce"):
    c = [float(p.get(f"c{i}")) / KJ for i in range(6)]
    if c[4] or c[5]:
        raise SystemExit("RB terms beyond cos^3 are not written")
    a = [c[n] * (-1) ** n for n in range(4)]
    b = {1: a[1] + 0.75 * a[3], 2: a[2] / 2, 3: a[3] / 4}
    const = a[0] + a[2] / 2 - sum(b.values())
    terms = [(v, n, 0.0) for n, v in b.items() if abs(v) > 1e-15]
    if abs(const) > 1e-15:
        terms.append((const, 0, 90.0))
    m = [p.get(f"class{i}") for i in (1, 2, 3, 4)]
    dihedrals.append({"name": "-".join(m), "match": m, "style": "fourier", "params": [len(terms)] + [x for t in terms for x in t],
                      "comment": "RB " + " ".join(f"c{i} {float(p.get(f'c{i}')):g}" for i in range(4)) + " kJ/mol (foyer)"})
s14 = [0, 0, float(nb.get("lj14scale", 0.5))]
doc = {"format": "caps-forcefield", "format_version": 1, "name": arg("--name"), "source": src.split("/")[-1] + " (foyer, mosdef-hub)",
       "references": many("--reference"), "units": "real",
       "styles": {"pair": "lj/cut/coul/long", "bond": "harmonic", "angle": "harmonic", "dihedral": "fourier", "improper": "harmonic"},
       "mixing": "arithmetic" if root.get("combining_rule", "lorentz") == "lorentz" else "geometric",
       "special_lj": s14, "special_coul": [0, 0, float(nb.get("coulomb14scale", 0.5))], "cutoff": float(arg("--cutoff", 14)),
       "equivalence": "replace", "atom_types": types, "pairs": pairs, "bonds": bonds, "angles": angles, "dihedrals": dihedrals, "impropers": [],
       "notes": [f"converted by CAPS (bench/ff/convert_foyer.py) from {src.split('/')[-1]} (foyer, MIT)",
                 "bonds harmonic with foyer's stiff constants (TraPPE fixes bond lengths)"],
       "version": arg("--version", root.get("version", ""))}
if "--tail" in sys.argv:
    doc["tail"] = True
if arg("--typing"):
    doc["typing"] = arg("--typing")
json.dump(doc, open(OUT, "w"), indent=1, ensure_ascii=False)
print(f"{OUT}: {len(types)} types, {len(bonds)} bonds, {len(angles)} angles, {len(dihedrals)} torsions")
