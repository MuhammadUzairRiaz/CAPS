#!/usr/bin/env python3
"""Molecular gas models for CAPS from their papers (PDFs in ~/project/ff-sources/papers/gases, never in the repository):

  epm2-harris1995   CO2, EPM2: J. G. Harris, K. H. Yung, J. Phys. Chem. 99, 12021 (1995), Table 1 (EPM2 rows): epsilon_CC
                    28.129 K, sigma_CC 2.757 A, epsilon_OO 80.507 K, sigma_OO 3.033 A, epsilon_CO 47.588 K, sigma_CO 2.892 A
                    (given, not mixed), l_CO 1.149 A, q_C +0.6512 e (each O half of it, negative), k_theta 1236 kJ/mol/rad^2
                    for (1/2) k (theta - 180)^2 (the flexible variant; bonds rigid)
  h2-buch1994       H2, one isotropic Lennard-Jones site at the centre: V. Buch, J. Chem. Phys. 100, 7610 (1994), section
                    II A: epsilon 34.2 K, sigma 2.96 A. The paper's anisotropic quadrupole-quadrupole term is left out.
  h2-darkrim1998    H2 with its quadrupole: F. Darkrim, D. Levesque, J. Chem. Phys. 109, 4981 (1998), p 4982: Lennard-
                    Jones sigma 0.2958 nm, epsilon 36.7 K at the centre of mass; charges q on the protons 0.0741 nm apart
                    and -2q at the centre. The paper prints q = 0.615e-26 esu, the units of a quadrupole moment
                    (esu cm^2): taken as H2's quadrupole Q, q = 2 Q / d^2 = 0.466378 e (Q = q d^2 / 2 for this charge set)

The H2 models keep the molecule's two hydrogens (their mass and, for Darkrim-Levesque, their charges) and put the centre
as a massless virtual site midway between them (as TraPPE's N2 centre): GROMACS writes it as a virtual site, LAMMPS
has no form for it. Molecules are rigid: run them as rigid bodies (bond constants 0).

usage: make_gas_models.py [DATA_DIR]
"""
import json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.dirname(HERE)), "data")
KB = 0.0019872043          # kcal/mol/K (CODATA 2018 k_B N_A)
KJ = 4.184
E_ESU = 4.803204712570263e-10   # elementary charge in esu (CODATA 2018)
M_C, M_O, M_H = 12.011, 15.9994, 1.008
STY = {"pair": "lj/cut/coul/long", "bond": "harmonic", "angle": "harmonic", "dihedral": "opls", "improper": "cvff"}

HARRIS = "J. G. Harris, K. H. Yung, J. Phys. Chem. 99, 12021 (1995)"
BUCH = "V. Buch, J. Chem. Phys. 100, 7610 (1994)"
DARKRIM = "F. Darkrim, D. Levesque, J. Chem. Phys. 109, 4981 (1998)"
D_HH = 0.741   # A (Darkrim & Levesque: protons 0.0741 nm apart)
Q_DARKRIM = 2 * 0.615e-26 / (D_HH * 1e-8) ** 2 / E_ESU


def base(name, source, refs, fid, aliases, mixing="arithmetic"):
    return {"format": "caps-forcefield", "format_version": 1, "name": name, "source": source, "references": refs, "units": "real",
            "styles": dict(STY), "mixing": mixing, "special_lj": [0, 0, 0], "special_coul": [0, 0, 0], "cutoff": 12,
            "atom_types": [], "pairs": [], "bonds": [], "angles": [], "dihedrals": [], "impropers": [],
            "notes": [f"built by CAPS (bench/ff/make_gas_models.py) from {source}"], "typing": f"../typing/{fid}.typing.json"}


def epm2():
    fid = "epm2-harris1995"
    d = base("CO2 (EPM2, Harris & Yung 1995)", HARRIS + ", Table 1 (EPM2)", [HARRIS], fid, ["epm2", "co2-epm2"])
    qc = 0.6512
    d["atom_types"] = [{"name": "C_co2", "element": "C", "mass": M_C, "charge": qc, "description": "CO2 carbon (EPM2)"},
                       {"name": "O_co2", "element": "O", "mass": M_O, "charge": -qc / 2, "description": "CO2 oxygen (EPM2)"}]
    d["pairs"] = [{"name": "C_co2", "match": ["C_co2"], "params": [28.129 * KB, 2.757], "comment": "epsilon_CC/k 28.129 K, sigma_CC 2.757 A"},
                  {"name": "O_co2", "match": ["O_co2"], "params": [80.507 * KB, 3.033], "comment": "epsilon_OO/k 80.507 K, sigma_OO 3.033 A"},
                  {"name": "C_co2-O_co2", "match": ["C_co2", "O_co2"], "params": [47.588 * KB, 2.892],
                   "comment": "epsilon_CO/k 47.588 K, sigma_CO 2.892 A (the paper's own cross term)"}]
    d["bonds"] = [{"name": "C_co2-O_co2", "match": ["C_co2", "O_co2"], "style": "harmonic", "params": [0.0, 1.149],
                   "comment": "l_CO 1.149 A; rigid bonds: run the molecule rigid"}]
    d["angles"] = [{"name": "O_co2-C_co2-O_co2", "match": ["O_co2", "C_co2", "O_co2"], "style": "harmonic", "params": [1236 / KJ / 2, 180.0],
                    "comment": "k_theta 1236 kJ/mol/rad^2 for (1/2) k (theta - 180)^2 (the flexible EPM2), LAMMPS K = k/2"}]
    d["notes"].append("bonds rigid (constant 0: constrain them or run the molecule as a rigid body), the bend flexible as the paper's flexible model")
    d.update(version="1995", status="validated")
    rules = [{"type": "C_co2", "smarts": "[CX2](=[OX1])=[OX1]", "priority": 0, "description": "O=[C]=O"},
             {"type": "O_co2", "smarts": "[OX1]=[CX2]=[OX1]", "priority": 0, "description": "[O]=C=O"}]
    return fid, d, rules, ["epm2", "co2-epm2"], "CO2 three sites, charges and Lennard-Jones of Table 1 (EPM2)"


def h2(fid, name, ref, source, eps_k, sigma, q, aliases):
    d = base(name, source, [ref], fid, aliases)
    d["atom_types"] = [{"name": "H_h2", "element": "H", "mass": M_H, "charge": q, "description": "H2 proton"},
                       {"name": "M_h2", "element": "", "mass": 0.0, "charge": -2 * q, "description": "H2 centre (massless)"}]
    d["pairs"] = [{"name": "H_h2", "match": ["H_h2"], "params": [0.0, 0.0], "comment": "charge only" if q else "no Lennard-Jones"},
                  {"name": "M_h2", "match": ["M_h2"], "params": [eps_k * KB, sigma], "comment": f"epsilon/k {eps_k} K, sigma {sigma} A at the centre"}]
    d["bonds"] = [{"name": "H_h2-H_h2", "match": ["H_h2", "H_h2"], "style": "harmonic", "params": [0.0, D_HH],
                   "comment": "H-H 0.741 A (Darkrim & Levesque: 0.0741 nm); rigid: run the molecule as a rigid body"}]
    d["virtual_sites"] = [{"smarts": "[#1][#1]", "type": "M_h2", "from": [0, 1], "w": [0.5, 0.5], "comment": "H2's centre, midway between the protons"}]
    d["notes"].append("the centre is a massless virtual site: GROMACS writes it as one, LAMMPS has no form for it")
    d.update(version=ref[-5:-1], status="validated")
    rules = [{"type": "H_h2", "smarts": "[#1][#1]", "priority": 0, "description": "a hydrogen of H2"},
             {"type": "M_h2", "smarts": "[#0]", "atom_name": "M_h2", "priority": 0, "description": "H2's centre (added by CAPS)"}]
    return fid, d, rules, aliases, "H2 by its two hydrogens; the centre added by the force field"


def main():
    libs = [epm2(),
            h2("h2-buch1994", "H2 (one Lennard-Jones site, Buch 1994)", BUCH, BUCH + ", section II A", 34.2, 2.96, 0.0, ["h2-buch"]),
            h2("h2-darkrim1998", "H2 with its quadrupole (Darkrim & Levesque 1998)", DARKRIM, DARKRIM + ", p 4982", 36.7, 2.958, Q_DARKRIM,
               ["h2-darkrim"])]
    cat_p = os.path.join(DATA, "forcefields", "catalogue.json")
    cat = json.load(open(cat_p))
    for fid, d, rules, aliases, evidence in libs:
        json.dump(d, open(os.path.join(DATA, "forcefields", fid + ".json"), "w"), indent=1, ensure_ascii=False)
        json.dump({"format": "caps-typing", "version": 1, "forcefield": d["name"], "unknown_types": "untyped",
                   "description": f"CAPS rules for {d['name']} (bench/ff/make_gas_models.py)", "rules": rules},
                  open(os.path.join(DATA, "typing", fid + ".typing.json"), "w"), indent=1, ensure_ascii=False)
        e = {"id": fid, "name": d["name"], "version": d["version"], "references": d["references"], "origin": "built from source",
             "source": d["source"] + "; bench/ff/make_gas_models.py", "license": "published parameters of the paper", "status": "validated",
             "notes": d["notes"][-1], "file": fid + ".json", "typing": {"rules": f"typing/{fid}.typing.json", "evidence": evidence},
             "counts": {"atom_types": len(d["atom_types"]), **{k: len(d.get(k, [])) for k in ("pairs", "bonds", "angles", "dihedrals", "impropers")}},
             "group": "gas", "aliases": aliases}
        cat["forcefields"] = [x for x in cat["forcefields"] if x["id"] != fid] + [e]
        print(fid, "written")
    json.dump(cat, open(cat_p, "w"), indent=1, ensure_ascii=False)


if __name__ == "__main__":
    main()
