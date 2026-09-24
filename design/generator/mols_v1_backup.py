"""Procedural structures for the mockups; geometry uses standard bond lengths."""
import math, random
from lib import add, sub, mul, norm, cross, dot3, dist


def polystyrene(n_units=8, tact="atactic", seed=3, hydrogens=True, curve=True):
    """Polystyrene chain: C-C 1.53 Å, aromatic C-C 1.39 Å, C-H 1.09 Å.
    Returns atoms, and indices of backbone carbons and of CH (stereo) carbons."""
    rnd = random.Random(seed)
    atoms, backbone, stereo = [], [], []
    nb = 2 * n_units
    dx, dy = 1.27, 0.43  # zigzag giving 1.53 Å bonds, ~112° angle
    def bb(k):
        x = (k - nb / 2) * dx
        y = dy * (1 if k % 2 == 0 else -1)
        z = 0.0
        if curve:
            z = 2.2 * math.sin(x / 5.5)
            y += 0.018 * x * x - 1.2
        return (x, y, z)
    for k in range(nb):
        atoms.append({"e": "C", "p": bb(k)})
        backbone.append(k)
    for k in range(nb):
        c = atoms[k]["p"]
        sy = 1 if k % 2 == 0 else -1
        # local tangent
        t = norm(sub(bb(k + 1), bb(k - 1)))
        up = norm(sub((c[0], c[1] + sy, c[2]), c))
        side = norm(cross(t, up))
        if k % 2 == 1:  # CH carrying phenyl
            stereo.append(k)
            if tact == "isotactic":
                s = 1
            elif tact == "syndiotactic":
                s = 1 if (k // 2) % 2 == 0 else -1
            else:
                s = rnd.choice([1, -1])
            d = norm(add(mul(up, 0.55), mul(side, s * 0.83)))
            ipso = add(c, mul(d, 1.51))
            cen = add(ipso, mul(d, 1.39))
            ring = []
            for m in range(6):
                th = math.pi + m * math.pi / 3
                p = add(cen, add(mul(d, 1.39 * math.cos(th)), mul(t, 1.39 * math.sin(th))))
                atoms.append({"e": "C", "p": p, "ar": True})
                ring.append(len(atoms) - 1)
            if hydrogens:
                for m in range(1, 6):
                    th = math.pi + m * math.pi / 3
                    p = add(cen, add(mul(d, 2.48 * math.cos(th)), mul(t, 2.48 * math.sin(th))))
                    atoms.append({"e": "H", "p": p})
                hd = norm(add(mul(up, 0.5), mul(side, -s * 0.87)))
                atoms.append({"e": "H", "p": add(c, mul(hd, 1.09))})
        elif hydrogens:
            for s in (1, -1):
                hd = norm(add(mul(up, 0.58), mul(side, s * 0.81)))
                atoms.append({"e": "H", "p": add(c, mul(hd, 1.09))})
    # end caps
    if hydrogens:
        for k, sgn in ((0, -1), (nb - 1, 1)):
            c = atoms[k]["p"]
            atoms.append({"e": "H", "p": add(c, (sgn * 1.0, 0.2 if k == 0 else -0.2, 0.35))})
    return atoms, backbone, stereo


def pe_crystal(na=2, nb=3, nc=5, a=7.40, b=4.93, c=2.534, hydrogens=False):
    """Orthorhombic PE (Bunn 1939): chains along c at (0,0) and (a/2,b/2), setting angles ±~42°."""
    atoms = []
    lat = 0.428  # half lateral zigzag width for 1.53 Å bonds
    for i in range(na):
        for j in range(nb):
            for (fx, fy, ang) in ((0.0, 0.0, 42.0), (0.5, 0.5, -42.0)):
                x0, y0 = (i + fx) * a, (j + fy) * b
                ca, sa = math.cos(math.radians(ang)), math.sin(math.radians(ang))
                for k in range(2 * nc):
                    s = 1 if k % 2 == 0 else -1
                    atoms.append({"e": "C", "p": (x0 + s * lat * sa, y0 + s * lat * ca, (k + 0.5) * c / 2)})
    return atoms


def silica_slab(nx=7, ny=7, nz=3, s=2.6):
    atoms = []
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                off = (s / 2 if k % 2 else 0.0)
                atoms.append({"e": "Si", "p": (i * s + off, j * s + off, k * s * 0.9)})
                if i < nx - 1:
                    atoms.append({"e": "O", "p": (i * s + off + s / 2, j * s + off, k * s * 0.9 + 0.35)})
                if j < ny - 1:
                    atoms.append({"e": "O", "p": (i * s + off, j * s + off + s / 2, k * s * 0.9 - 0.35)})
    return atoms


def random_chain(start, n, seed, bond=1.54, box=None, persistence=0.72):
    rnd = random.Random(seed)
    pts = [start]
    d = norm((rnd.uniform(-1, 1), rnd.uniform(-1, 1), rnd.uniform(-1, 1)))
    for _ in range(n - 1):
        r = norm((rnd.uniform(-1, 1), rnd.uniform(-1, 1), rnd.uniform(-1, 1)))
        d = norm(add(mul(d, persistence), mul(r, 1 - persistence)))
        p = add(pts[-1], mul(d, bond))
        if box:
            lo, hi = box
            p = tuple(min(max(p[k], lo[k]), hi[k]) for k in range(3))
            for k in range(3):
                if p[k] in (lo[k], hi[k]):
                    d = tuple(-d[m] if m == k else d[m] for m in range(3))
        pts.append(p)
    return pts


def nanotube(n=10, m=10, length=26.0, acc=1.42):
    a = acc * math.sqrt(3)
    a1 = (a, 0.0)
    a2 = (a / 2, a * math.sqrt(3) / 2)
    Ch = (n * a1[0] + m * a2[0], n * a1[1] + m * a2[1])
    L = math.hypot(*Ch)
    R = L / (2 * math.pi)
    u = (Ch[0] / L, Ch[1] / L)
    v = (-u[1], u[0])
    basis = [(0.0, 0.0), ((a1[0] + a2[0]) / 3, (a1[1] + a2[1]) / 3)]
    pts = []
    rng = int(max(n, m) * 3 + length / a * 2 + 10)
    for i in range(-rng, rng):
        for j in range(-rng, rng):
            for bx, by in basis:
                x = i * a1[0] + j * a2[0] + bx
                y = i * a1[1] + j * a2[1] + by
                s = x * u[0] + y * u[1]
                t = x * v[0] + y * v[1]
                if 0 <= s < L - 1e-6 and 0 <= t < length:
                    th = 2 * math.pi * s / L
                    pts.append((R * math.cos(th), R * math.sin(th), t))
    atoms = [{"e": "C", "p": p} for p in pts]
    bonds = []
    for i in range(len(pts)):
        for j in range(i + 1, len(pts)):
            if dist(pts[i], pts[j]) < acc + 0.12:
                bonds.append((i, j))
    return atoms, bonds, R


def helix_ca(nres=22, rise=1.5, radius=2.3, per_turn=3.6):
    pts = []
    for i in range(nres):
        th = 2 * math.pi * i / per_turn
        pts.append((radius * math.cos(th), radius * math.sin(th), i * rise))
    return pts


def water(o, rnd):
    d1 = norm((rnd.uniform(-1, 1), rnd.uniform(-1, 1), rnd.uniform(-1, 1)))
    tmp = norm((rnd.uniform(-1, 1), rnd.uniform(-1, 1), rnd.uniform(-1, 1)))
    perp = norm(cross(d1, tmp))
    ang = math.radians(104.52 / 2)
    h1 = add(o, mul(norm(add(mul(d1, math.cos(ang)), mul(perp, math.sin(ang)))), 0.9572))
    h2 = add(o, mul(norm(add(mul(d1, math.cos(ang)), mul(perp, -math.sin(ang)))), 0.9572))
    return [{"e": "O", "p": o}, {"e": "H", "p": h1}, {"e": "H", "p": h2}], [(0, 1), (0, 2)]
