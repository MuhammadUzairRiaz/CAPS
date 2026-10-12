"""Builds zeolite frameworks for data/crystals from the IZA-SC Database of Zeolite Structures (Ch. Baerlocher, L. B.
McCusker, http://www.iza-structure.org/databases/): each framework type's idealised SiO2 cell, space group and asymmetric
unit (refined with DLS-76) as the database gives it, built with CAPS's own crystal builder and checked atom by atom
against the CIF's own symmetry operators (the setting CAPS needs, e.g. Fd-3m's origin choice, is the one that agrees).
The CIFs stay local (~/project/ff-sources/zeolites, one per code); CAPS keeps its own P1 cells.

usage: python3 tools/make_zeolites.py [CIF_DIR] [--write]
"""
import itertools, json, math, os, re, subprocess, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CAPS = os.path.join(ROOT, "build", "cli", "caps")
OUT = os.path.join(ROOT, "data", "crystals")
SRC = next((a for a in sys.argv[1:] if not a.startswith("--")), os.path.expanduser("~/project/ff-sources/zeolites"))
REF = "Ch. Baerlocher, L. B. McCusker, Database of Zeolite Structures, IZA Structure Commission (http://www.iza-structure.org/databases/)"
# code: name, typical use (the materials CAPS users build them for)
NAMES = {"SOD": "sodalite", "LTA": "zeolite A (Linde Type A)", "FAU": "faujasite (zeolite X / Y)", "EMT": "EMC-2 (hexagonal faujasite)",
         "MFI": "ZSM-5 / silicalite-1", "MEL": "ZSM-11 / silicalite-2", "MTW": "ZSM-12", "TON": "theta-1 / ZSM-22", "FER": "ferrierite",
         "MOR": "mordenite", "CHA": "chabazite (SSZ-13, SAPO-34)", "AEI": "AlPO-18 / SSZ-39", "AFI": "AlPO-5", "AEL": "AlPO-11",
         "LTL": "zeolite L", "MWW": "MCM-22", "GIS": "gismondine (zeolite P, MAP)", "ANA": "analcime", "NAT": "natrolite",
         "HEU": "heulandite / clinoptilolite", "RHO": "zeolite rho", "KFI": "ZK-5", "ERI": "erionite", "OFF": "offretite",
         "MAZ": "mazzite (omega)", "DDR": "deca-dodecasil 3R", "ITE": "ITQ-3", "STT": "SSZ-23", "IFR": "ITQ-4"}


def num(x):
    return float(re.sub(r"\(.*\)", "", x))


def parse(path):
    t = open(path).read()
    cell = [num(re.search(r"_cell_" + k + r"\s+(\S+)", t).group(1)) for k in ("length_a", "length_b", "length_c", "angle_alpha", "angle_beta", "angle_gamma")]
    hm = re.search(r"_symmetry_space_group_name_H-M\s+'([^']+)'", t).group(1)
    no = int(re.search(r"_symmetry_Int_Tables_number\s+(\d+)", t).group(1))
    ops = [l.strip().strip("'").replace(" ", "") for l in t.split("_symmetry_equiv_pos_as_xyz", 1)[1].split("loop_", 1)[0].splitlines()
           if l.strip().count(",") == 2]   # quoted ('+x,+y,+z') or bare (X,Y,Z) operators
    sites = []
    for l in t.split("_atom_site_fract_z", 1)[1].splitlines():
        w = l.split()
        if len(w) >= 5:
            sites.append((w[0], w[1], num(w[2]), num(w[3]), num(w[4])))
    return cell, hm, no, ops, sites


def expand(ops, sites):
    out = []
    for lab, el, x, y, z in sites:
        for op in ops:
            p = [eval(e.replace("x", f"({x})").replace("y", f"({y})").replace("z", f"({z})"), {}) % 1.0 for e in op.lower().split(",")]
            p = [v % 1.0 if abs(v % 1.0 - 1) > 1e-4 else 0.0 for v in p]
            if not any(e == el and all(min(abs(a - b), 1 - abs(a - b)) < 1e-3 for a, b in zip(p, q)) for e, q in out):
                out.append((el, p))
    return out


def built(path):
    t = open(path).read()
    atoms = []
    for l in t.split("_atom_site_fract_z", 1)[1].splitlines():
        w = l.split()
        if len(w) >= 5:
            atoms.append((w[1], [float(w[2]) % 1.0, float(w[3]) % 1.0, float(w[4]) % 1.0]))
    return atoms


def same(a, b):
    if len(a) != len(b):
        return False
    left = list(b)
    for el, p in a:
        k = next((i for i, (e, q) in enumerate(left) if e == el and all(min(abs(u - v), 1 - abs(u - v)) < 2e-3 for u, v in zip(p, q))), None)
        if k is None:
            return False
        left.pop(k)
    return True


def main():
    write = "--write" in sys.argv
    tmp = os.path.join(ROOT, "build", "zeolites")
    os.makedirs(tmp, exist_ok=True)
    made, bad = [], 0
    for code in sorted(NAMES):
        path = os.path.join(SRC, code + ".cif")
        if not os.path.exists(path):
            print(f"{code:4s} no CIF in {SRC}")
            continue
        cell, hm, no, ops, sites = parse(path)
        want = expand(ops, sites)
        site_txt = "; ".join(f"{lab} {el} {x} {y} {z}" for lab, el, x, y, z in sites)
        ok_grp = None
        for grp in (f"{no}", f"{no}:2", f"{no}:1", f"{no}:h", hm):
            out = os.path.join(tmp, code + ".cif")
            r = subprocess.run([CAPS, "crystal", "--group", grp, "--cell", ",".join(str(x) for x in cell), "--sites", site_txt, "-o", out],
                               capture_output=True, text=True)
            if r.returncode == 0 and os.path.exists(out) and same(want, built(out)):
                ok_grp = grp
                break
        nt = sum(1 for e, _ in want if e == "Si")
        a, b, c, al, be, ga = cell
        ca, cb, cg = (math.cos(math.radians(x)) for x in (al, be, ga))
        vol = a * b * c * math.sqrt(1 - ca * ca - cb * cb - cg * cg + 2 * ca * cb * cg)
        fd = 1000 * nt / vol
        print(f"{code:4s} {hm:12s} {len(want):4d} atoms ({nt} T)  FD {fd:5.1f} T/1000 A^3  {'ok (' + ok_grp + ')' if ok_grp else 'CHECK: no setting reproduces the CIF'}")
        if not ok_grp:
            bad += 1
            continue
        made.append((code, hm, ok_grp, cell, site_txt, fd, nt))
        if write:
            name = NAMES[code]
            head = (f"# {code} framework, {name} (pure SiO2) — {hm} ({ok_grp}) a b c α β γ = {' '.join(str(x) for x in cell)}; sites {site_txt}\n"
                    f"# {REF}: idealised SiO2 framework refined with DLS-76; built by tools/make_zeolites.py (framework density {fd:.1f} T/1000 A^3)\n")
            open(os.path.join(OUT, f"zeolite-{code.lower()}.cif"), "w").write(head + open(os.path.join(tmp, code + ".cif")).read())
    if write and made:
        cat_path = os.path.join(OUT, "catalogue.json")
        text = open(cat_path).read()
        have = {c["id"] for c in json.loads(text)["crystals"]}
        lines = [json.dumps({"id": f"zeolite-{code.lower()}", "name": f"{code} framework ({NAMES[code]}), SiO2", "file": f"zeolite-{code.lower()}.cif",
                             "use": "zeolite framework (IZA-SC idealised SiO2): adsorption, sorption, catalysis", "hkl": [0, 0, 1]}, ensure_ascii=False)
                 for code, *_ in made if f"zeolite-{code.lower()}" not in have]
        if lines:
            end = text.rindex("}", 0, text.rindex("]"))
            text = text[:end + 1] + "".join(",\n    " + l for l in lines) + text[end + 1:]
            json.loads(text)
            open(cat_path, "w").write(text)
    print(f"\n{len(made)} frameworks agree with their CIF atom by atom; {bad} need a look")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
