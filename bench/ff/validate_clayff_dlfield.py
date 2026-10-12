#!/usr/bin/env python3
"""CAPS's ClayFF typing (data/typing/clayff-cygan2004.typing.json) against DL_FIELD's automatic ClayFF typing, atom by atom.

DL_FIELD is the accuracy reference only: nothing of it is read into the library. Each structure is written as an xyz
file with one molecular group of DL_FIELD's general_CLAYFF template (MOLECULE_KEY CLYF) and its periodic cell, and
DL_FIELD types it by element, coordination and neighbours (its own neighbour distances, minimum image); its types
(dlf_notation.output) are mapped to the paper's species by their ATOM_TYPE table (INORGANIC_clay.sf), where some
carry other names for the same species:
  o2h, o3h -> oh      obb (two octahedral Al and a Si) -> ob      obc, obt -> obts      obo -> obos      obs -> obss
CAPS types the same xyz file with `caps ff type`. Every atom is compared. One difference is the paper's: an O next to a
tetrahedral Al and an octahedral Mg / Fe / Li is obss in the paper (p 1259: "two octahedral substitutions or
tetrahedral-octahedral substitutions") and obo (obos) in DL_FIELD; such atoms are counted apart, every other atom must
agree. The default set adds IFF's muscovite with every third octahedral Al made Mg, which holds both kinds of obss.

DL_FIELD stops on structures outside its clay rules: brucite (an OH with three Mg) and LAMMPS's msi2lmp pyrophyllite
(PyAC_bulk-clayff.car, whose hydroxyl H 889 lies 1.31 A from Al 880, inside DL_FIELD's Al-H distance).

usage: validate_clayff_dlfield.py [STRUCTURE ...]   (default: the clays of the INTERFACE force field's model database:
       pyrophyllite, kaolinite, muscovite, montmorillonites with Na or K)
Environment: DLFIELD (default ~/project/dl_f_4.13).
"""
import os, re, shutil, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
CAPS = os.path.join(ROOT, "build", "cli", "caps")
FF = os.path.join(ROOT, "data", "forcefields", "clayff-cygan2004.json")
SRC = os.environ.get("DLFIELD", os.path.expanduser("~/project/dl_f_4.13"))
WORK = os.path.join(ROOT, "build", "dlfield_work_clay")
IFF = os.path.expanduser("~/iff-ref/INTERFACE_FF_1_5/MODEL_DATABASE/CLAY_MINERALS")
DEFAULT = [os.path.join(IFF, f + ".car") for f in ("pyrophyllite30_cell", "kaolinite30_cell", "mica15_cell", "mont0_333_Na_15_cell",
                                                   "mont0_333_K_15_cell", "mont0_400_Na_15_cell", "mont0_533_Na_15_cell", "mont0_321_Na_28_cell")]
SAME = {"o2h": "oh", "o3h": "oh", "obb": "ob", "obc": "obts", "obt": "obts", "obo": "obos", "obs": "obss"}


def dl_keys():
    """DL_FIELD ATOM_TYPE -> ATOM_KEY of its clay library"""
    txt = open(os.path.join(SRC, "lib", "INORGANIC_clay.sf")).read()
    block = re.search(r"^ATOM_TYPE(.*?)^END ATOM_TYPE", txt, re.S | re.M).group(1)
    return {w[0]: w[1] for w in (l.split() for l in block.splitlines()) if len(w) >= 3 and not w[0].startswith("#")}


def read_xyz(path):
    lines = open(path).read().splitlines()
    n = int(lines[0])
    lat = re.search(r'Lattice="([^"]+)"', lines[1])
    cell = [float(v) for v in lat.group(1).split()] if lat else None
    atoms = [l.split()[:4] for l in lines[2:2 + n]]
    return atoms, cell


def run_dlfield(atoms, cell):
    if not os.path.exists(os.path.join(WORK, "dl_field")):
        shutil.copytree(SRC, WORK, ignore=shutil.ignore_patterns("output", "Examples", "source", "*.zip"))
    out = os.path.join(WORK, "out")
    shutil.rmtree(out, ignore_errors=True)
    os.makedirs(out)
    with open(os.path.join(WORK, "caps_clay.xyz"), "w") as f:
        f.write(f"{len(atoms)}\nCAPS ClayFF typing check\n# MOLECULE_KEY CLYF MOLECULE_GROUP A\n")
        for el, x, y, z in atoms:
            f.write(f"{el} {x} {y} {z}\n")
    with open(os.path.join(WORK, "dl_f_path"), "w") as f:
        f.write("library = lib/\nsolvent = solvent/\noutput = out/\ncontrol = caps_clay.control\n")
    lines = open(os.path.join(SRC, "control_files", "ethanol.control")).read().splitlines()
    lines[0] = "CAPS ClayFF typing check"
    lines[2] = "none     * Seconday output files (gromacs, chemshell or none)."
    lines[4] = "inorganic_clay * Type of force field require (see list below for choices)."
    lines[9] = "none     * Include user-defined information. Put 'none' or a .udff filename"
    lines[11] = "caps_clay.xyz  * Configuration file."
    for k, l in enumerate(lines):
        if "Periodic condition" in l:
            lines[k] = "3        * Periodic condition ? 0=no, other number = type of box (see below)"
            for v in range(3):
                lines[k + 1 + v] = " ".join(f"{c:.6f}" for c in cell[3 * v:3 * v + 3]) + f"  * Cell vector {'abc'[v]} (x, y, z)"
            break
    open(os.path.join(WORK, "caps_clay.control"), "w").write("\n".join(lines) + "\n")
    r = subprocess.run(["./dl_field"], cwd=WORK, capture_output=True, text=True)
    found = [os.path.join(d, "dlf_notation.output") for d, _, fs in os.walk(out) if "dlf_notation.output" in fs]
    path = found[0] if found else ""
    if not path:
        raise RuntimeError("DL_FIELD wrote no dlf_notation.output: " + " | ".join(r.stdout.strip().splitlines()[-4:]))
    types, nbs = [], []
    for l in open(path):
        w = l.split()
        if len(w) >= 4 and w[0].isdigit():
            types.append(w[1])
            nbs.append([int(x) - 1 for x in w[4:4 + int(w[3])]])
    if not types:
        raise RuntimeError(" | ".join(r.stdout.strip().splitlines()[-4:]))
    return types, nbs


def paper_obss(i, atoms, nbs):
    """an O next to a tetrahedral Al and an octahedral substitution: obss in the paper ("a bridging oxygen coordinated to
    two octahedral substitutions or tetrahedral-octahedral substitutions", p 1259), obo (obos) in DL_FIELD"""
    el = [atoms[j][0] for j in nbs[i]]
    return any(atoms[j][0] == "Al" and len(nbs[j]) == 4 for j in nbs[i]) and any(e in ("Mg", "Fe", "Li") for e in el)


def substituted_mica(work):
    """IFF's muscovite with every third octahedral Al (as CAPS types them) made Mg: oxygens next to a tetrahedral and an
    octahedral substitution, and next to two octahedral ones"""
    src = os.path.join(IFF, "mica15_cell.car")
    if not os.path.exists(src):
        return None
    xyz, tf, out = (os.path.join(work, n) for n in ("mica15.xyz", "mica15.types", "mica15_mg.xyz"))
    subprocess.run([CAPS, "convert", src, xyz], check=True, capture_output=True)
    subprocess.run([CAPS, "ff", "type", xyz, "--ff", FF, "-o", tf], check=True, capture_output=True)
    ao = [k + 1 for k, l in enumerate(open(tf)) if l.split() and l.split()[-1] == "ao"]
    subprocess.run([CAPS, "edit", xyz, "--ops", "element " + ",".join(map(str, ao[2::3])) + " Mg", "-o", out], check=True, capture_output=True)
    return out


def main():
    work = tempfile.mkdtemp()
    files = [a for a in sys.argv[1:] if not a.startswith("-")]
    if not files:
        files = DEFAULT + [substituted_mica(work)]
    key = dl_keys()
    total_ok = total = total_dev = 0
    failed = False
    for path in files:
        if not path or not os.path.exists(path):
            print(f"{os.path.basename(path or 'mica15_mg')}: missing, skipped")
            continue
        name = os.path.splitext(os.path.basename(path))[0]
        xyz = os.path.join(work, name + ".in.xyz")
        subprocess.run([CAPS, "convert", path, xyz], check=True, capture_output=True)
        atoms, cell = read_xyz(xyz)
        tf = os.path.join(work, name + ".in.types")
        subprocess.run([CAPS, "ff", "type", xyz, "--ff", FF, "-o", tf], check=True, capture_output=True)
        caps = [l.split()[-1] for l in open(tf) if l.strip()]
        try:
            dlt, nbs = run_dlfield(atoms, cell)
            dl = [SAME.get(key[t], key[t]) for t in dlt]
        except (RuntimeError, KeyError) as e:
            print(f"{name}: DL_FIELD could not type it ({e})")
            failed = True
            continue
        if len(dl) != len(caps) or len(dl) != len(atoms):
            print(f"{name}: atom counts differ (CAPS {len(caps)}, DL_FIELD {len(dl)}, file {len(atoms)})")
            failed = True
            continue
        dev = [i for i, (c, d) in enumerate(zip(caps, dl)) if c != d and (c, d) == ("obss", "obos") and paper_obss(i, atoms, nbs)]
        bad = [(i, atoms[i][0], c, d) for i, (c, d) in enumerate(zip(caps, dl)) if c != d and i not in set(dev)]
        comp = {}
        for t in caps:
            comp[t] = comp.get(t, 0) + 1
        ok = len(dl) - len(bad) - len(dev)
        total_ok += ok
        total_dev += len(dev)
        total += len(dl)
        print(f"{name}: {ok} / {len(dl)} atoms as DL_FIELD types them" + (f", {len(dev)} O next to a tetrahedral and an "
              "octahedral substitution obss as the paper names them (DL_FIELD obos)" if dev else "")
              + f"  ({' '.join(f'{t} {n}' for t, n in sorted(comp.items()))})")
        for i, el, c, d in bad[:20]:
            print(f"   atom {i + 1} {el}: CAPS {c}, DL_FIELD {d}")
        failed = failed or bool(bad)
    print(f"total: {total_ok} / {total} as DL_FIELD" + (f", {total_dev} where the paper's obss differs from DL_FIELD's obos" if total_dev else "")
          + (", every other atom as DL_FIELD" if not failed else ""))
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
