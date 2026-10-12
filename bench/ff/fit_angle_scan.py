#!/usr/bin/env python3
"""A bend a force field lacks, fitted to a quantum-chemical relaxed scan of that angle.

    fit_angle_scan.py SCANDIR MOL2 FF.json I J K [--overlay OV.json] [--json out.json]
    fit_angle_scan.py SCANDIR MOL2 FF.json I J --bond [...]   (a relaxed scan of a bond length: E = K (r - r0)^2)
                      [--class CA-SY-CA] [--molecule "diphenyl sulfone"] [--level "B3LYP-D3BJ/def2-TZVP"]

SCANDIR holds an ORCA relaxed scan (%geom Scan A I J K = lo, hi, n end end): the optimised geometry of every point
(*.allxyz, or the numbered *.NNN.xyz) and their energies. MOL2 is the model molecule with its bonds, atoms in the
order of the scan (I J K counted from 0, as ORCA counts). Each point's geometry is evaluated with CAPS (FF.json, the
missing bend at zero; --overlay adds terms fitted already, such as the molecule's missing torsions), and

    E_QM(theta) - E_MM(theta) = K (theta - theta0)^2 + c

is fitted by least squares (a theta^2 + b theta + c: K = a, theta0 = -b / 2a), K in kcal/mol/rad^2 as OPLS-AA's
harmonic bend E = K (theta - theta0)^2. --json keeps the whole profile (QM, MM without the bend, MM with the fitted
bend) for the paper's figures (docs/paper/scripts/plots/ff_derive.py).
"""
import glob, json, math, os, re, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "data", "python"))
from caps import torsionfit as T   # noqa: E402  (mol2 writing, CAPS energies)

HARTREE = 627.509474


def frames(path):
    """geometries of a multi-frame xyz file: [(comment, elements, positions)]"""
    out, lines, k = [], open(path).read().splitlines(), 0
    while k < len(lines):
        if not lines[k].strip():
            k += 1
            continue
        n = int(lines[k].split()[0])
        el, pos = [], []
        for l in lines[k + 2: k + 2 + n]:
            w = l.split()
            el.append(re.sub(r"\d", "", w[0]))
            pos.append([float(x) for x in w[1:4]])
        out.append((lines[k + 1], el, pos))
        k += 2 + n
    return out


def angle(p, i, j, k):
    u, v = T.sub(p[i], p[j]), T.sub(p[k], p[j])
    return math.degrees(math.acos(max(-1.0, min(1.0, T.dot(u, v) / T.norm(u) / T.norm(v)))))


def read_mol2_bonds(path):
    bonds, sec = [], None
    for l in open(path):
        if l.startswith("@<TRIPOS>"):
            sec = l.strip()[9:]
            continue
        w = l.split()
        if sec == "BOND" and len(w) >= 3:
            bonds.append((int(w[1]) - 1, int(w[2]) - 1))
    return bonds


def lstsq3(x, y):
    """a x^2 + b x + c by least squares (normal equations, 3 x 3)"""
    rows = [[t * t, t, 1.0] for t in x]
    A = [[sum(r[i] * r[j] for r in rows) for j in range(3)] for i in range(3)]
    v = [sum(r[i] * yy for r, yy in zip(rows, y)) for i in range(3)]
    for c in range(3):   # Gauss-Jordan
        piv = max(range(c, 3), key=lambda r: abs(A[r][c]))
        A[c], A[piv], v[c], v[piv] = A[piv], A[c], v[piv], v[c]
        for r in range(3):
            if r != c:
                f = A[r][c] / A[c][c]
                A[r] = [a - f * b for a, b in zip(A[r], A[c])]
                v[r] -= f * v[c]
    return [v[i] / A[i][i] for i in range(3)]


def main():
    a = sys.argv[1:]
    bond = "--bond" in a
    if len(a) < (5 if bond else 6):
        raise SystemExit(__doc__)
    scan, mol2, ff = a[0], a[1], os.path.abspath(a[2])
    if bond:
        i, j = map(int, a[3:5])
        k = None
    else:
        i, j, k = map(int, a[3:6])
    overlay = a[a.index("--overlay") + 1] if "--overlay" in a else None
    allxyz = sorted(glob.glob(os.path.join(scan, "*.allxyz")))
    if allxyz:   # the scan job's own trajectory (the last .allxyz: the scan is the second job)
        pts = frames(allxyz[-1])
        E = [float(re.search(r"E\s+(-?\d+\.\d+)", c).group(1)) for c, _, _ in pts]
    else:
        files = sorted(f for f in glob.glob(os.path.join(scan, "*.[0-9][0-9][0-9].xyz")))
        pts = [frames(f)[0] for f in files]
        E = [float(re.search(r"E\s+(-?\d+\.\d+)", c).group(1)) for c, _, _ in pts]
    bonds = read_mol2_bonds(mol2)
    th, eqm, emm = [], [], []
    tmp = os.path.join(scan, "fit")
    os.makedirs(tmp, exist_ok=True)
    for n, ((_, el, pos), e) in enumerate(zip(pts, E)):
        m2 = os.path.join(tmp, f"p{n:02d}.mol2")
        T.write_mol2(m2, el, pos, bonds)
        emm.append(T.caps_apply(m2, ff, os.path.join(tmp, f"p{n:02d}.data"), overlay))
        eqm.append(e * HARTREE)
        th.append(T.norm(T.sub(pos[i], pos[j])) if bond else angle(pos, i, j, k))
    y = [q - m for q, m in zip(eqm, emm)]
    A, B, C = lstsq3(th, y)
    K = A if bond else A * (180.0 / math.pi) ** 2   # a bond: per A^2; an angle: per deg^2 -> per rad^2
    th0 = -B / (2 * A)
    fit = [A * t * t + B * t + C for t in th]
    rms = math.sqrt(sum((f - v) ** 2 for f, v in zip(fit, y)) / len(y))
    lo = min(eqm)
    print(f"{'r' if bond else 'theta':>8} {'E_QM':>9} {'E_MM(no term)':>14} {'fit':>9}")
    for t, q, m, f, v in zip(th, eqm, emm, fit, y):
        print(f"{t:8.2f} {q - lo:9.3f} {m:14.3f} {f - min(fit):9.3f}")
    print(f"K = {K:.2f} kcal/mol/{'A' if bond else 'rad'}^2, {'r0' if bond else 'theta0'} = {th0:.4f} "
          f"{'A' if bond else 'deg'}, RMS {rms:.3f} kcal/mol over {len(th)} points")
    if "--json" in a:
        opt = lambda k, d="": a[a.index(k) + 1] if k in a else d
        bend = [K * ((t - th0) if bond else math.radians(t - th0)) ** 2 for t in th]
        mlo = min(m + b for m, b in zip(emm, bend))
        x = "r" if bond else "theta"
        json.dump({"kind": "bond" if bond else "angle", "class": opt("--class"), "molecule": opt("--molecule"),
                   "level": opt("--level"), "atoms": [i, j] if bond else [i, j, k], "K": round(K, 2),
                   ("r0" if bond else "theta0"): round(th0, 4 if bond else 2), "rms": round(rms, 3), "points": len(th),
                   "scan": [{x: round(t, 4), "E_qm": round(q - lo, 4), "E_mm_without": round(m - min(emm), 4),
                             "E_mm_with": round(m + b - mlo, 4)} for t, q, m, b in zip(th, eqm, emm, bend)]},
                  open(opt("--json"), "w"), indent=1)


if __name__ == "__main__":
    main()
