#!/usr/bin/env python3
"""Oxide glass potentials for CAPS from their papers (each term with its table):

teter-glass        rigid ions with Teter's partial charges (O -1.2, Si +2.4, Al/Y/RE +1.8, Na +0.6), Buckingham terms with O
                   only (O-O and cation-O; cation-cation Coulomb only):
                     A. N. Cormack, J. Du, T. R. Zeitler, Phys. Chem. Chem. Phys. 4, 3193 (2002), Table 1 (O-O, Si-O, Na-O)
                     J. Du, J. Am. Ceram. Soc. 92, 87 (2009), Table II (Al-O, Y-O)
                     K. Okhotnikov, B. Stevensson, M. Eden, Phys. Chem. Chem. Phys. 15, 15041 (2013), Table 1 (La, Lu, Sc and
                     its own Y-O, the overlay okhotnikov2013-y)
                   overlay tilocca2006-teter: A. Tilocca, N. H. de Leeuw, A. N. Cormack, Phys. Rev. B 73, 104209 (2006),
                   Table III (its O-O 2029.2204 eV, C 192.28; Ca +1.2). The papers splice a short-range repulsion onto the
                   Buckingham term for melts (Tilocca's b/r^n has no switching radius in the paper): not written.
tilocca2006-sm     the shell model of Tilocca et al. 2006, Table I: formal charges, O core +0.8482 / shell -2.8482, spring
                   74.92 eV/A^2, shell 0.2 amu; O-Si-O bends on the O shells: SM1 (default) the screened harmonic
                   1/2 k (theta - theta0)^2 exp(-(r1 + r2)/rho), k 100 eV/rad^2, rho 1.0 A, r_max 2.5 A; SM0 (overlay sm0)
                   1/2 k (theta - theta0)^2, k 2.097, r_max 1.9 A
vogel2003-lipo3    M. Vogel, arXiv:cond-mat/0307738v2 (2003), Table I: Li +0.6, P +3.0, O -1.2, A exp(-r/0.29) for every pair

usage: make_glass.py [DATA_DIR]
"""
import json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.dirname(HERE)), "data")
EV = 23.060548
MASS = {"O": 15.999, "Si": 28.0855, "Al": 26.982, "Y": 88.906, "La": 138.905, "Lu": 174.967, "Sc": 44.956, "Na": 22.98977,
        "Ca": 40.078, "Li": 6.94, "P": 30.973762}
Z = {"O": 8, "Si": 14, "Al": 13, "Y": 39, "La": 57, "Lu": 71, "Sc": 21, "Na": 11, "Ca": 20, "Li": 3, "P": 15}
CORMACK = "A. N. Cormack, J. Du, T. R. Zeitler, Phys. Chem. Chem. Phys. 4, 3193 (2002)"
DU = "J. Du, J. Am. Ceram. Soc. 92, 87 (2009)"
OKH = "K. Okhotnikov, B. Stevensson, M. Eden, Phys. Chem. Chem. Phys. 15, 15041 (2013)"
TIL = "A. Tilocca, N. H. de Leeuw, A. N. Cormack, Phys. Rev. B 73, 104209 (2006)"
VOGEL = "M. Vogel, arXiv:cond-mat/0307738v2 (2003)"


def write(fid, doc, typing, overlays):
    doc.update({"format": "caps-forcefield", "format_version": 1, "units": "real", "typing": f"../typing/{fid}.typing.json"})
    doc.setdefault("styles", {"pair": "buck/coul/long", "bond": "harmonic", "angle": "harmonic", "dihedral": "fourier", "improper": "fourier"})
    doc.setdefault("mixing", "geometric"); doc.setdefault("special_lj", [0, 0, 1]); doc.setdefault("special_coul", [0, 0, 1])
    for k in ("bonds", "angles", "dihedrals", "impropers"):
        doc.setdefault(k, [])
    json.dump(doc, open(os.path.join(DATA, "forcefields", fid + ".json"), "w"), indent=1)
    json.dump(dict({"format": "caps-typing", "version": 1, "forcefield": doc["name"], "unknown_types": "untyped", "bonds": "defined"}, **typing),
              open(os.path.join(DATA, "typing", fid + ".typing.json"), "w"), indent=1)
    if overlays:
        os.makedirs(os.path.join(DATA, "typing", "sets", fid), exist_ok=True)
        for name, (desc, rs) in overlays.items():
            json.dump({"format": "caps-typing", "version": 1, "forcefield": f"{fid}: {name}", "unknown_types": "untyped", "description": desc,
                       "rules": [dict(r, priority=50) for r in rs]}, open(os.path.join(DATA, "typing", "sets", fid, name + ".typing.json"), "w"), indent=1)
    print(f"{fid}: {len(doc['atom_types'])} types, {len(doc['pairs'])} pair entries, {len(overlays)} overlays")


def rule(t, el, others, d):
    return {"type": t, "smarts": f"[#{Z[el]}]", "priority": 0, "requires": others, "description": d}


# ---- teter-glass
T_CH = {"O": -1.2, "Si": 2.4, "Al": 1.8, "Y": 1.8, "La": 1.8, "Lu": 1.8, "Sc": 1.8, "Na": 0.6}
# (cation, A eV, rho A, C eV A^6, source)
T_PAIRS = [("O", 1844.7458, 0.343645, 192.58, CORMACK + ", Table 1; " + DU + ", Table II; " + OKH + ", Table 1"),
           ("Si", 13702.905, 0.193817, 54.681, CORMACK + ", Table 1; " + DU + ", Table II; " + OKH + ", Table 1"),
           ("Na", 4383.7555, 0.243838, 30.700, CORMACK + ", Table 1"),
           ("Al", 12201.417, 0.195628, 31.997, DU + ", Table II; " + OKH + ", Table 1"),
           ("Y", 29526.977, 0.211377, 50.477, DU + ", Table II; " + OKH + ", Table 1"),
           ("La", 9537.95, 0.235917, 0.0, OKH + ", Table 1"),
           ("Lu", 6657.09, 0.230360, 0.0, OKH + ", Table 1"),
           ("Sc", 3095.04, 0.244357, 0.0, OKH + ", Table 1")]
types = [{"name": f"{el}_te", "element": el, "mass": MASS[el], "charge": q, "description": f"{el} ({q:+g}, Teter's partial charge)"} for el, q in T_CH.items()]
types.append({"name": "Y_ok", "element": "Y", "mass": MASS["Y"], "charge": 1.8, "description": "Y (+1.8) with Okhotnikov's own Y-O (Table 1)"})
# Tilocca 2006, Table III: (A, rho, C) with O; charges
TIL_T = {"O": (2029.22040, 0.343645, 192.28, -1.2), "Si": (13702.9050, 0.193817, 54.681, 2.4), "Na": (4383.75550, 0.243838, 30.700, 0.6),
         "Ca": (7747.18340, 0.252623, 93.109, 1.2)}
types += [{"name": f"{el}_tt", "element": el, "mass": MASS[el], "charge": v[3], "description": f"{el} ({v[3]:+g}) in Tilocca's Teter set (Table III)"}
          for el, v in TIL_T.items()]
pairs = [{"name": t["name"], "match": [t["name"]], "params": [0.0, 0.0]} for t in types]
for el, A, rho, C, src in T_PAIRS:
    a = "O_te" if el == "O" else f"{el}_te"
    pairs.append({"name": f"{a}-O_te", "match": [a, "O_te"], "style": "buck", "params": [A * EV, rho, C * EV], "comment": f"A {A} eV, rho {rho} A, C {C} eV A^6 ({src})"})
pairs.append({"name": "Y_ok-O_te", "match": ["Y_ok", "O_te"], "style": "buck", "params": [11019.6 * EV, 0.224173, 0.0], "comment": f"A 11019.6 eV, rho 0.224173 A ({OKH}, Table 1: derived there)"})
for el, (A, rho, C, q) in TIL_T.items():
    a = f"{el}_tt"
    pairs.append({"name": f"{a}-O_tt", "match": [a, "O_tt"], "style": "buck", "params": [A * EV, rho, C * EV], "comment": f"A {A} eV, rho {rho} A, C {C} eV A^6 ({TIL}, Table III)"})
els = list(T_CH)
rules = [rule(f"{el}_te", el, ["O"] if el != "O" else [], f"{el} in an oxide glass") for el in els]
overlays = {"okhotnikov2013-y": ("Overlay: Okhotnikov's Y-O (Born-Mayer, derived in the paper) in place of Du's", [rule("Y_ok", "Y", ["O"], "Y, Okhotnikov's Y-O")]),
            "tilocca2006-teter": ("Overlay: Tilocca et al.'s Teter set (Table III) for O, Si, Na, Ca", [rule(f"{el}_tt", el, [], f"{el} (Tilocca's Teter set)") for el in TIL_T])}
write("teter-glass", {"name": "Oxide glasses, Teter partial charges (Cormack 2002; Du 2009; Okhotnikov 2013)",
                      "source": f"{CORMACK}, Table 1; {DU}, Table II; {OKH}, Table 1; overlay: {TIL}, Table III",
                      "references": [CORMACK, DU, OKH, TIL], "cutoff": 8, "atom_types": types, "pairs": pairs,
                      "notes": ["built by CAPS (bench/ff/make_glass.py) from the papers' tables",
                                "Buckingham terms with O only; cation-cation pairs interact by Coulomb alone",
                                "short-range terms cut at 8 A as Okhotnikov and Tilocca; the short-range splice the papers use for melts is not written",
                                "Ca is in the tilocca2006-teter overlay only (Table III)"],
                      "version": "2002-2013"},
      {"description": "each element by itself (bench/ff/make_glass.py); Okhotnikov's Y and Tilocca's Teter set are overlays", "rules": rules}, overlays)

# ---- tilocca2006-sm
SM_CH = {"Si": 4.0, "Na": 1.0, "Ca": 2.0}
types = [{"name": f"{el}_sm", "element": el, "mass": MASS[el], "charge": q, "description": f"{el} ({q:+g}, rigid ion; Table I)"} for el, q in SM_CH.items()]
types += [{"name": "Si_sm0", "element": "Si", "mass": MASS["Si"], "charge": 4.0, "description": "Si (+4), SM0 bends"},
          {"name": "O_smc", "element": "O", "mass": MASS["O"] - 0.2, "charge": 0.8482, "description": "O core (Table I: Y -2.8482)"},
          {"name": "O_sms", "element": "O", "mass": 0.2, "charge": -2.8482, "description": "O shell (0.2 amu, as the paper)"},
          {"name": "O_sm0c", "element": "O", "mass": MASS["O"] - 0.2, "charge": 0.8482, "description": "O core, SM0"},
          {"name": "O_sm0s", "element": "O", "mass": 0.2, "charge": -2.8482, "description": "O shell, SM0"}]
pairs = [{"name": t["name"], "match": [t["name"]], "params": [0.0, 0.0]} for t in types]
SM_P = [("O", 22764.30, 0.14900, 27.88), ("Si", 1283.91, 0.320520, 10.661580), ("Na", 56465.3453, 0.193931, 0.0), ("Ca", 2152.3566, 0.309227, 0.099440)]
for shell, si in (("O_sms", "Si_sm"), ("O_sm0s", "Si_sm0")):
    for el, A, rho, C in SM_P:
        a = shell if el == "O" else (si if el == "Si" else f"{el}_sm")
        pairs.append({"name": f"{a}-{shell}", "match": [a, shell], "style": "buck", "params": [A * EV, rho, C * EV],
                      "comment": f"{el}-Os A {A} eV, rho {rho} A, C {C} eV A^6 (Table I)"})
bonds = [{"name": f"{c}-{s}", "match": [c, s], "style": "class2", "params": [0.0, 74.92 / 2 * EV, 0.0, 0.0], "comment": "Oc-Os 74.92 eV/A^2 (Table I)"}
         for c, s in (("O_smc", "O_sms"), ("O_sm0c", "O_sm0s"))]
angles = [{"name": "O_sms-Si_sm-O_sms", "match": ["O_sms", "Si_sm", "O_sms"], "style": "screened", "params": [100.0 / 2 * EV, 109.47, 1.0, 1.0], "contact": [2.5],
           "comment": "SM1: 1/2 k3b (theta - theta0)^2 exp(-(r1/rho + r2/rho)), k3b 100.0 eV/rad^2, rho 1.0 A, r_max 2.5 A (Table I, eq. 4)"},
          {"name": "O_sm0s-Si_sm0-O_sm0s", "match": ["O_sm0s", "Si_sm0", "O_sm0s"], "style": "harmonic", "params": [2.097 / 2 * EV, 109.47], "contact": [1.9],
           "comment": "SM0: 1/2 k3b (theta - theta0)^2, k3b 2.097 eV/rad^2, r_max 1.9 A (Table I, eq. 3)"}]
shells = {"O_smc": "O_sms", "O_sm0c": "O_sm0s"}
rules = [rule("Si_sm", "Si", ["O"], "Si4+"), rule("Na_sm", "Na", ["O"], "Na+"), rule("Ca_sm", "Ca", ["O"], "Ca2+"), rule("O_smc", "O", [], "O2- core")]
rules += [{"type": s, "smarts": "[#8]", "priority": 99, "atom_name": s, "description": f"shell of {c}"} for c, s in shells.items()]
overlays = {"sm0": ("Overlay: SM0, the raw harmonic O-Si-O bend (Table I)", [rule("Si_sm0", "Si", ["O"], "Si4+, SM0"), rule("O_sm0c", "O", [], "O2- core, SM0")])}
write("tilocca2006-sm", {"name": "Silicate glasses, shell model (Tilocca, de Leeuw & Cormack 2006)", "source": f"{TIL}, Table I", "references": [TIL],
                         "styles": {"pair": "buck/coul/long", "bond": "class2", "angle": "harmonic", "dihedral": "fourier", "improper": "fourier"},
                         "cutoff": 8, "atom_types": types, "pairs": pairs, "bonds": bonds, "angles": angles,
                         "notes": ["built by CAPS (bench/ff/make_glass.py) from the paper's Table I",
                                   "SM1 by default (the screened O-Si-O bend: CAPS's own engine and DL_POLY; LAMMPS and GROMACS have no such style); SM0 is the overlay 'sm0'",
                                   "short-range terms cut at 8 A as the paper"],
                         "version": "2006 (Table I)"},
      {"shells": shells, "description": "Si, Na, Ca and shell-model O (bench/ff/make_glass.py); SM0 is an overlay", "rules": rules}, overlays)

# ---- vogel2003-lipo3
V_CH = {"Li": 0.6, "P": 3.0, "O": -1.2}
V_A = {("Li", "Li"): 167.47, ("Li", "P"): 158.09, ("Li", "O"): 300.00, ("P", "P"): 148.93, ("P", "O"): 694.66, ("O", "O"): 2644.8}
types = [{"name": f"{el}_v", "element": el, "mass": MASS[el], "charge": q, "description": f"{el} ({q:+g}; the paper's effective charge)"} for el, q in V_CH.items()]
pairs = [{"name": t["name"], "match": [t["name"]], "params": [0.0, 0.0]} for t in types]
pairs += [{"name": f"{a}_v-{b}_v", "match": [f"{a}_v", f"{b}_v"], "style": "buck", "params": [A * EV, 0.29, 0.0], "comment": f"A {A} eV, rho 0.29 A (Table I)"}
          for (a, b), A in V_A.items()]
write("vogel2003-lipo3", {"name": "LiPO3 glass (Vogel 2003)", "source": f"{VOGEL}, Table I and eq. 1", "references": [VOGEL], "cutoff": 10,
                          "atom_types": types, "pairs": pairs,
                          "notes": ["built by CAPS (bench/ff/make_glass.py) from the paper's Table I (A_LiO reduced by the paper; the rest after Karthikeyan et al.)"],
                          "version": "2003"},
      {"description": "Li, P and O of lithium metaphosphate (bench/ff/make_glass.py)",
       "rules": [rule("Li_v", "Li", ["P", "O"], "Li+"), rule("P_v", "P", ["O"], "P"), rule("O_v", "O", ["P"], "O")]}, {})
