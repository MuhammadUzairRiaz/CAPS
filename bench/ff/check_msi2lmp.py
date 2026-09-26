#!/usr/bin/env python3
"""CVFF, PCFF and COMPASS from their .frc files (data/forcefields/*-frc.json) against LAMMPS's own msi2lmp tool.

msi2lmp's reference data files (lammps/tools/msi2lmp/test/reference: class1 = cvff.frc, class2a = compass_published.frc,
class2b = pcff.frc) give each structure with the types Materials Studio assigned and the parameters msi2lmp found for
them. CAPS reads the same structure with the same types (--types), assigns its parameters from the .frc conversion
and writes its data file and input in the force field's own LAMMPS styles; LAMMPS then evaluates both with the same
settings (Coulomb and van der Waals cut at 15 Å, special_bonds 0 0 1, no tail) and the energy terms are compared.

Known msi2lmp behaviour CAPS does not copy (reported, not counted as differences):
  · bond-bond_1_3 terms only for dihedrals with a type named cp… (COMPASS's aromatic c3a gets none from msi2lmp);
  · angle-angle terms of three-connected centres looked up by equivalence only when msi2lmp runs with -p 3.

usage: check_msi2lmp.py [--only NAME] [--keep DIR]
"""
import os, re, subprocess, sys, tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
CAPS = os.path.join(ROOT, "build", "cli", "caps")
LMP = os.environ.get("LMP", os.path.expanduser("~/lammps/build-class2/lmp"))
TEST = os.path.expanduser("~/lammps/tools/msi2lmp/test")
arg = lambda k, d: sys.argv[sys.argv.index(k) + 1] if k in sys.argv else d
only = arg("--only", "")
work = arg("--keep", "") or tempfile.mkdtemp()
FF = {"class1": "cvff-frc.json", "class2a": "compass-frc.json", "class2b": "pcff-frc.json"}
TERMS = ["ebond", "eangle", "edihed", "eimp", "evdwl", "ecoul"]


def lammps(inp, cwd):
    r = subprocess.run([LMP, "-in", inp, "-log", "none", "-screen", "none"], cwd=cwd, capture_output=True, text=True)
    out = os.path.join(cwd, "thermo.txt")
    if not os.path.exists(out):
        raise RuntimeError((r.stdout + r.stderr).strip().splitlines()[-1] if (r.stdout + r.stderr).strip() else "LAMMPS failed")
    vals = open(out).read().split()
    os.remove(out)
    return dict(zip(TERMS, map(float, vals)))


def section(text, name):
    out, on = [], False
    for line in text.splitlines():
        if line.startswith(name):
            on = True
            continue
        if on and line and line[0].isalpha():
            break
        if on and line.split():
            out.append(line.split("#")[0].split())
    return out


def explained(ref_text, caps_text, types):
    """For CVFF, whether the differences are msi2lmp's own shortcuts only: (1) an out-of-plane term's atoms in the order
    the .mdf file lists the centre's neighbours (the same parameters for every centre), (2) a wildcard torsion's barrier
    divided by the neighbour count of the first atom of each central type rather than of the atoms themselves.
    Returns (True, note) or (False, what differs)."""
    nb = {}
    for w in section(ref_text, "Bonds"):
        a, b = int(w[2]), int(w[3])
        nb.setdefault(a, set()).add(b)
        nb.setdefault(b, set()).add(a)
    first_nc = {}
    for i, t in enumerate(types, 1):
        first_nc.setdefault(t, len(nb.get(i, ())))

    def terms(text, sec, key):
        co = {int(w[0]): w[1:] for w in section(text, sec + " Coeffs")}
        m = {}
        for w in section(text, sec + "s"):
            a = tuple(map(int, w[2:6]))
            m.setdefault(key(a), []).append(co[int(w[1])])
        return m

    def dnorm(v):
        out = []
        for x in v or []:
            if len(x) == 3:   # harmonic K d n
                out.append((float(x[0]), int(x[2]), int(x[1])))
            else:             # fourier m (K n phase)…
                for q in range(int(x[0])):
                    out.append((float(x[1 + 3 * q]), int(x[2 + 3 * q]), 1 if float(x[3 + 3 * q]) == 0 else -1))
        return sorted(o for o in out if o[0] != 0)
    dk = lambda a: a if a[1] < a[2] else a[::-1]
    R, C = terms(ref_text, "Dihedral", dk), terms(caps_text, "Dihedral", dk)
    by_count = 0
    for k in set(R) | set(C):
        r, c = dnorm(R.get(k)), dnorm(C.get(k))
        if r == c:
            continue
        # the barrier spread over the torsions about j-k: by the atoms' own neighbours (CAPS) or by type (msi2lmp)
        j, kk = k[1], k[2]
        own = (len(nb[j]) - 1) * (len(nb[kk]) - 1)
        typ = (first_nc[types[j - 1]] - 1) * (first_nc[types[kk - 1]] - 1)
        if len(r) == len(c) and own != typ and all(abs(a[0] * typ - b[0] * own) < 2e-3 * max(1, own, typ) and a[1:] == b[1:] for a, b in zip(r, c)):
            by_count += 1
            continue
        return False, f"torsion {[types[i - 1] for i in k]}: msi2lmp {r}, CAPS {c}"
    ck = lambda a: a[1]
    R, C = terms(ref_text, "Improper", ck), terms(caps_text, "Improper", ck)
    reordered = 0
    for k in set(R) | set(C):
        r = sorted((round(float(x[0]), 4), int(x[1]), int(x[2])) for x in R.get(k, []) if float(x[0]) != 0)
        c = sorted((round(float(x[0]), 4), int(x[1]), int(x[2])) for x in C.get(k, []) if float(x[0]) != 0)
        if r != c:
            return False, f"out-of-plane at {types[k - 1]} {k}: msi2lmp {r}, CAPS {c}"
        reordered += 1
    notes = []
    if by_count:
        notes.append(f"{by_count} wildcard torsions spread over the atoms' own neighbours (msi2lmp counts by type)")
    notes.append("out-of-plane parameters equal at every centre (atom order: msi2lmp's from the .mdf file)")
    return True, "; ".join(notes)


FOOT = ("thermo_style custom step pe ebond eangle edihed eimp evdwl ecoul\nrun 0\n"
        "print \"$(ebond:%.10f) $(eangle:%.10f) $(edihed:%.10f) $(eimp:%.10f) $(evdwl:%.10f) $(ecoul:%.10f)\" file thermo.txt\n")


def main():
    rows, fails = [], 0
    for fn in sorted(os.listdir(os.path.join(TEST, "reference"))):
        m = re.match(r"(.*)-(class1|class2a|class2b)\.data$", fn)
        if not m or (only and only not in fn):
            continue
        name, cls = m.groups()
        ref = os.path.join(TEST, "reference", fn)
        d = os.path.join(work, name + "-" + cls)
        os.makedirs(d, exist_ok=True)
        text = open(ref).read()
        # types by the Masses comments, atoms in id order
        # the types Materials Studio assigned, from the .car file msi2lmp read (atoms in the same order)
        atoms = []
        for line in open(os.path.join(TEST, name + "-" + cls + ".car")).read().splitlines()[4:]:
            w = line.split()
            if len(w) >= 9 and w[0] != "end":
                atoms.append((len(atoms) + 1, w[6][:4]))
        open(os.path.join(d, "types.txt"), "w").write("\n".join(t for _, t in sorted(atoms)) + "\n")
        refin = open(os.path.join(TEST, "in." + name + "-" + cls)).read()
        boundary = (re.search(r"^boundary.*$", refin, re.M) or [None])[0] if re.search(r"^boundary.*$", refin, re.M) else "boundary p p p"
        styles = [l for l in refin.splitlines() if re.match(r"(bond|angle|dihedral|improper)_style", l)]
        pair = "lj/class2/coul/cut 15.0" if cls != "class1" else "lj/cut/coul/cut 15.0"
        open(os.path.join(d, "ref.in"), "w").write(
            f"units real\n{boundary}\natom_style full\npair_style {pair}\n" + "\n".join(styles) +
            f"\nspecial_bonds lj/coul 0.0 0.0 1.0\nread_data {ref}\n" + FOOT)
        r = subprocess.run([CAPS, "ff", "apply", ref, "--ff", os.path.join(ROOT, "data", "forcefields", FF[cls]), "--types",
                            os.path.join(d, "types.txt"), "--charges", "keep", "-o", os.path.join(d, "caps.data"),
                            "--lammps-input", os.path.join(d, "caps.in"), "--kspace", "cut", "--lammps-cutoff", "15"],
                           capture_output=True, text=True)
        if r.returncode:
            rows.append((fn, "CAPS failed: " + (r.stderr.strip().splitlines() or ["?"])[-1][:200], ""))
            fails += 1
            continue
        auto = re.findall(r"automatic[^\n]*", r.stdout)
        # compared as msi2lmp writes: bonded coefficients rounded to its 4 decimals, and bond-bond-1-3 terms only
        # when the structure has a type named cp… (msi2lmp's rule; the .frc file gives them to any matching torsion)
        has_cp = any(t.startswith("cp") for _, t in atoms)
        data, sec, out, zeroed = open(os.path.join(d, "caps.data")).read().splitlines(), "", [], 0
        for line in data:
            w = line.split("#")[0].split()
            if line and line[0].isalpha():
                sec = line.split("#")[0].strip()
            elif w and sec.endswith("Coeffs") and not sec.startswith("Pair") and not sec.startswith("PairIJ"):
                vals = [w[0]] + [f"{round(float(x), 4):.4f}" if re.match(r"^-?\d*\.\d*(e-?\d+)?$|^-?\d+e-?\d+$", x) else x for x in w[1:]]
                if sec.startswith("BondBond13") and not has_cp and float(vals[1]) != 0:
                    vals[1] = "0.0000"
                    zeroed += 1
                line = " ".join(vals)
            out.append(line)
        open(os.path.join(d, "caps.data"), "w").write("\n".join(out) + "\n")
        cin = open(os.path.join(d, "caps.in")).read()
        cin = re.sub(r"(?m)^boundary.*$", boundary, cin)
        cin = re.sub(r"(?m)^pair_modify.*\n", "", cin)
        cin = re.sub(r"(?ms)^thermo_style.*", "", cin) + FOOT
        open(os.path.join(d, "caps.in"), "w").write(cin)
        try:
            ge, ce = lammps("ref.in", d), lammps("caps.in", d)
        except RuntimeError as ex:
            rows.append((fn, "LAMMPS failed: " + str(ex)[:200], ""))
            fails += 1
            continue
        diff = {k: ce[k] - ge[k] for k in TERMS}
        worst = max(abs(diff[k]) / max(1.0, abs(ge[k])) for k in TERMS)
        ok = worst < 1e-6
        why = ""
        if not ok and cls == "class1" and max(abs(diff[k]) / max(1.0, abs(ge[k])) for k in TERMS if k not in ("edihed", "eimp")) < 1e-6:
            ok, why = explained(text, open(os.path.join(d, "caps.data")).read(), [t for _, t in sorted(atoms)])
            why = (" · " if ok else " · ") + why
        fails += 0 if ok else 1
        styl = " ".join(l.split(None, 1)[1] for l in open(os.path.join(d, "caps.in")).read().splitlines()
                        if re.match(r"(bond|angle|dihedral|improper)_style", l))
        quirk = f" · {zeroed} bond-bond-1-3 types left out as msi2lmp does (no cp type)" if zeroed else ""
        rows.append((fn, f"{'ok' if ok else 'DIFFERS'} · {len(atoms)} atoms · largest relative difference {worst:.1e} · CAPS styles: {styl}{quirk}{why}",
                     "" if ok else "   msi2lmp/CAPS: " + "  ".join(f"{k} {ge[k]:.6f}/{ce[k]:.6f}" for k in TERMS)))
    for fn, res, more in rows:
        print(f"{fn}\n   {res}" + (f"\n{more}" if more else ""))
    print(f"\n{sum(1 for r in rows if ' ok ' in ' ' + r[1][:3] + ' ')} of {len(rows)} msi2lmp reference structures agree with CAPS (work: {work})")
    sys.exit(1 if fails else 0)


if __name__ == "__main__":
    main()
