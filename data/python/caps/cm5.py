"""CM5 charges from Hirshfeld charges, and OPLS-AA charges from CM5 (1.27 x CM5).

    python -m caps.cm5 ORCA.out [--scale 1.27] [--no-average] [--json out.json]
    python -m caps.cm5 --check CM5PAC_TESTS_DIR

CM5 (Marenich, Jerome, Cramer & Truhlar, J. Chem. Theory Comput. 8, 527 (2012)) maps Hirshfeld charges to charges
that reproduce gas-phase dipoles:

    q_k(CM5) = q_k(Hirshfeld) + sum_{k' != k} T_{Z_k Z_k'} B_{kk'},   B_{kk'} = exp[-alpha (r_kk' - R_Z_k - R_Z_k')]

alpha = 2.474 1/A, R the covalent radii the paper takes (CRC Handbook, 91st ed.), T_ZZ' = D_Z - D_Z' except the six
pairs among H, C, N and O the paper fits apart. The numbers below are the paper's (Table 1, eq. 4-5); the radii past
the paper's elements are those CM5PAC (Apache 2.0) uses.

OPLS-AA takes 1.27 x CM5 for neutral molecules (Dodda, Vilseck, Tirado-Rives & Jorgensen, J. Phys. Chem. B 121, 3864
(2017); Vilseck, Tirado-Rives & Jorgensen, J. Chem. Theory Comput. 10, 2802 (2014)), computed from the Hirshfeld
charges of a density at a good level (CAPS uses M06-2X/6-311+G(2df,2p)), with symmetry-equivalent atoms given their
mean and the total kept exact.

ORCA prints Hirshfeld charges (%output Print[P_Hirshfeld] 1 end) but not CM5; this module adds it.
"""
import json, math, os, re, sys

ALPHA = 2.474
# covalent radii (A) by atomic number, 1..118
RAD = [0.32, 0.37, 1.3, 0.99, 0.84, 0.75, 0.71, 0.64, 0.6, 0.62, 1.6, 1.4, 1.24, 1.14, 1.09, 1.04, 1.0, 1.01, 2.0,
       1.74, 1.59, 1.48, 1.44, 1.3, 1.29, 1.24, 1.18, 1.17, 1.22, 1.2, 1.23, 1.2, 1.2, 1.18, 1.17, 1.16, 2.15, 1.9,
       1.76, 1.64, 1.56, 1.46, 1.38, 1.36, 1.34, 1.3, 1.36, 1.4, 1.42, 1.4, 1.4, 1.37, 1.36, 1.36, 2.38, 2.06, 1.94,
       1.84, 1.9, 1.88, 1.86, 1.85, 1.83, 1.82, 1.81, 1.8, 1.79, 1.77, 1.77, 1.78, 1.74, 1.64, 1.58, 1.5, 1.41, 1.36,
       1.32, 1.3, 1.3, 1.32, 1.44, 1.45, 1.5, 1.42, 1.48, 1.46, 2.42, 2.11, 2.01, 1.9, 1.84, 1.83, 1.8, 1.8, 1.73,
       1.68, 1.68, 1.68, 1.65, 1.67, 1.73, 1.76, 1.61, 1.57, 1.49, 1.43, 1.41, 1.34, 1.29, 1.28, 1.21, 1.22, 1.36,
       1.43, 1.62, 1.75, 1.65, 1.57]
# D_Z by atomic number, 1..118 (Table 1; zero where the paper fits none)
D = [0.0056, -0.1543, 0.0, 0.0333, -0.103, -0.0446, -0.1072, -0.0802, -0.0629, -0.1088, 0.0184, 0.0, -0.0726,
     -0.079, -0.0756, -0.0565, -0.0444, -0.0767, 0.013] + [0.0] * 11 + [
     -0.0512, -0.0557, -0.0533, -0.0399, -0.0313, -0.0541, 0.0092] + [0.0] * 11 + [
     -0.0361, -0.0393, -0.0376, -0.0281, -0.022, -0.0381, 0.0065] + [0.0] * 25 + [
     -0.0255, -0.0277, -0.0265, -0.0198, -0.0155, -0.0269, 0.0046] + [0.0] * 25 + [
     -0.0179, -0.0195, -0.0187, -0.014, -0.011, -0.0189]
PAIR = {(1, 6): 0.0502, (1, 7): 0.1747, (1, 8): 0.1671, (6, 7): 0.0556, (6, 8): 0.0234, (7, 8): -0.0346}
SYMBOL = ("H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr "
          "Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe").split()
ZOF = {s: i + 1 for i, s in enumerate(SYMBOL)}


def T(z1, z2):
    if (z1, z2) in PAIR:
        return PAIR[(z1, z2)]
    if (z2, z1) in PAIR:
        return -PAIR[(z2, z1)]
    return D[z1 - 1] - D[z2 - 1]


def dist(a, b):
    return math.sqrt(sum((a[k] - b[k]) ** 2 for k in range(3)))


def cm5(z, xyz, q):
    """CM5 charges of atoms with atomic numbers z, positions xyz (A) and Hirshfeld charges q"""
    out = list(q)
    for k in range(len(z)):
        for j in range(len(z)):
            if z[j] != z[k]:
                out[k] += T(z[k], z[j]) * math.exp(-ALPHA * (dist(xyz[k], xyz[j]) - RAD[z[k] - 1] - RAD[z[j] - 1]))
    return out


# ------------------------------------------------------------------ symmetry classes
def bonds_of(z, xyz, tol=0.45):
    return [(i, j) for i in range(len(z)) for j in range(i + 1, len(z))
            if dist(xyz[i], xyz[j]) < RAD[z[i] - 1] + RAD[z[j] - 1] + tol]


def classes(z, xyz):
    """topologically equivalent atoms (element and bonded surroundings, refined to convergence)"""
    nb = [[] for _ in z]
    for i, j in bonds_of(z, xyz):
        nb[i].append(j)
        nb[j].append(i)
    lab = list(z)
    while True:
        keys = [(lab[i], tuple(sorted(lab[j] for j in nb[i]))) for i in range(len(z))]
        ids = {k: n for n, k in enumerate(sorted(set(keys)))}
        new = [ids[k] for k in keys]
        if len(set(new)) == len(set(lab)):
            return new
        lab = new


def opls(z, xyz, q_cm5, scale=1.27, average=True, total=None):
    """scale x CM5, equivalent atoms averaged, total kept (the rounded charge of the molecule unless given)"""
    q = [scale * c for c in q_cm5]
    if average:
        cl = classes(z, xyz)
        mean = {c: sum(q[i] for i in range(len(q)) if cl[i] == c) / cl.count(c) for c in set(cl)}
        q = [mean[cl[i]] for i in range(len(q))]
    want = round(sum(q_cm5)) if total is None else total
    shift = (want - sum(q)) / len(q)
    return [c + shift for c in q]


# ------------------------------------------------------------------ readers
def read_orca(path):
    """(atomic numbers, last geometry in A, Hirshfeld charges) of an ORCA output"""
    text = open(path, errors="replace").read()
    geo = text.rfind("CARTESIAN COORDINATES (ANGSTROEM)")
    hir = text.rfind("HIRSHFELD ANALYSIS")
    if geo < 0 or hir < 0:
        raise SystemExit(f"{path}: no geometry or no Hirshfeld analysis (add %output Print[P_Hirshfeld] 1 end)")
    z, xyz = [], []
    for line in text[geo:].splitlines()[2:]:
        w = line.split()
        if len(w) != 4:
            break
        z.append(ZOF[re.sub(r"\d", "", w[0])])
        xyz.append([float(v) for v in w[1:]])
    q = []
    for line in text[hir:].splitlines():
        m = re.match(r"\s*(\d+)\s+([A-Z][a-z]?)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s*$", line)
        if m:
            q.append(float(m.group(3)))
        elif q:
            break
    if len(q) != len(z):
        raise SystemExit(f"{path}: {len(z)} atoms but {len(q)} Hirshfeld charges")
    return z, xyz, q


def check(tests):
    """CM5 against CM5PAC's own test outputs (geometry from testN.inp, Hirshfeld and CM5 from testN_RevA02.out)"""
    worst = 0.0
    for n in range(1, 20):
        inp, out = os.path.join(tests, f"test{n}.inp"), os.path.join(tests, f"test{n}_RevA02.out")
        if not (os.path.exists(inp) and os.path.exists(out)):
            continue
        if "CM5M" in open(out).read():
            print(f"test{n}: CM5M (the metal-cluster variant), not CM5; skipped")
            continue
        z, xyz = [], []
        for line in open(inp):
            w = line.split()
            if len(w) == 4 and re.match(r"^\d+$", w[0]) and all(re.match(r"^-?\d+\.\d+$", v) for v in w[1:]):
                z.append(int(w[0]))
                xyz.append([float(v) for v in w[1:]])
        ref, hir = [], []
        for line in open(out):
            w = line.split()
            if len(w) == 4 and re.match(r"^\d+$", w[0]) and re.match(r"^-?\d+\.\d+$", w[2]):
                ref.append(float(w[2]))
                hir.append(float(w[3]))
        if len(ref) != len(z):
            print(f"test{n}: {len(z)} atoms in the input, {len(ref)} in the output; skipped")
            continue
        got = cm5(z, xyz, hir)
        dev = max(abs(a - b) for a, b in zip(got, ref))
        worst = max(worst, dev)
        print(f"test{n}: {len(z)} atoms, max |CM5 - CM5PAC| = {dev:.6f}")
    print(f"worst {worst:.6f}")
    return worst


def main():
    a = sys.argv[1:]
    if not a:
        raise SystemExit(__doc__)
    if a[0] == "--check":
        sys.exit(0 if check(a[1]) < 5e-6 else 1)
    z, xyz, q = read_orca(a[0])
    scale = float(a[a.index("--scale") + 1]) if "--scale" in a else 1.27
    c = cm5(z, xyz, q)
    o = opls(z, xyz, c, scale, "--no-average" not in a)
    cl = classes(z, xyz)
    print(f"{'atom':>5} {'el':>3} {'class':>5} {'Hirshfeld':>10} {'CM5':>10} {f'{scale}xCM5':>10}")
    for i in range(len(z)):
        print(f"{i + 1:5d} {SYMBOL[z[i] - 1]:>3} {cl[i]:5d} {q[i]:10.5f} {c[i]:10.5f} {o[i]:10.5f}")
    print(f"total {sum(q):.5f} {sum(c):.5f} {sum(o):.5f}")
    if "--json" in a:
        json.dump({"source": os.path.basename(a[0]), "elements": [SYMBOL[v - 1] for v in z], "hirshfeld": q,
                   "cm5": c, "opls": o, "scale": scale, "classes": cl},
                  open(a[a.index("--json") + 1], "w"), indent=1)


if __name__ == "__main__":
    main()
