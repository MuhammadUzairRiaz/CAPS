"""Procedural structures for the mockups; geometry uses standard bond lengths."""
import math, random
from lib import add, sub, mul, norm, cross, dot3, dist


def _place(a, b, c, bond, ang, tor):
    """NeRF: position of d given a-b-c, |cd| = bond, angle bcd = ang (deg), torsion abcd = tor (deg)."""
    bc = norm(sub(c, b))
    n = norm(cross(sub(b, a), bc))
    m = cross(n, bc)
    A, T = math.radians(ang), math.radians(tor)
    d2 = (-bond * math.cos(A), bond * math.sin(A) * math.cos(T), bond * math.sin(A) * math.sin(T))
    return add(c, add(add(mul(bc, d2[0]), mul(m, d2[1])), mul(n, d2[2])))


def _subst(c, u1, u2, s):
    """Unit vector for a tetrahedral substituent on atom c with neighbour directions u1, u2; s = +1 / -1 picks the side."""
    bis = norm(add(u1, u2))
    nn = norm(cross(u1, u2))
    return norm(add(mul(bis, -0.575), mul(nn, s * 0.818)))


def polystyrene(n_units=8, tact="atactic", seed=3, hydrogens=True, curve=True):
    """All-atom polystyrene built from internal coordinates: backbone C–C 1.53 Å, C–C–C 114°,
    C(H)–C(ipso) 1.51 Å, planar phenyl ring C–C 1.39 Å, aromatic C–H 1.08 Å, aliphatic C–H 1.09 Å.
    Each step places the next backbone atom and the finished atom's substituents (H included) together;
    backbone torsion and ring rotation are chosen by scoring clashes against everything placed before,
    like a small configurational-bias growth. curve=False restricts torsions to near-trans.
    Atom order: backbone first, then per backbone atom its substituents (CH: 6 ring C, 5 ring H, 1 H;
    CH2: 2 H), then the two end-cap H. Returns atoms, backbone indices, stereo (CH) indices."""
    rnd = random.Random(seed)
    nb = 2 * n_units
    sides = []
    for k in range(nb):
        if k % 2 == 1:
            sides.append(1 if tact == "isotactic" else ((1 if (k // 2) % 2 == 0 else -1) if tact == "syndiotactic" else rnd.choice([1, -1])))
        else:
            sides.append(0)
    LIM = {("C", "C"): 3.0, ("C", "H"): 2.45, ("H", "H"): 2.0}

    def lim(a, b):
        return LIM[(a, b) if a <= b else (b, a)]

    def ring(bk, d, phi):
        cen = add(bk, mul(d, 1.51 + 1.39))
        ref = norm(cross(d, (0.0, 0.0, 1.0) if abs(d[2]) < 0.9 else (1.0, 0.0, 0.0)))
        w = norm(add(mul(ref, math.cos(phi)), mul(cross(d, ref), math.sin(phi))))
        C = [add(cen, add(mul(d, 1.39 * math.cos(math.pi + m * math.pi / 3)), mul(w, 1.39 * math.sin(math.pi + m * math.pi / 3)))) for m in range(6)]
        H = [add(cen, add(mul(d, 2.47 * math.cos(math.pi + m * math.pi / 3)), mul(w, 2.47 * math.sin(math.pi + m * math.pi / 3)))) for m in range(1, 6)]
        return C, H

    def subs(j, prev, bj, nxt, phi):
        """Substituents of backbone atom j: list of (elem, pos, tag)."""
        u1, u2 = norm(sub(prev, bj)), norm(sub(nxt, bj))
        out = []
        if j % 2 == 1:
            d = _subst(bj, u1, u2, sides[j])
            C, H = ring(bj, d, phi)
            out += [("C", p, "ring") for p in C] + [("H", p, "ringH") for p in H]
            out.append(("H", add(bj, mul(_subst(bj, u1, u2, -sides[j]), 1.09)), "H"))
        else:
            out += [("H", add(bj, mul(_subst(bj, u1, u2, sv), 1.09)), "H") for sv in (1, -1)]
        return out

    # placed atoms: id -> (elem, pos); adjacency for topological exclusion; spatial grid for speed
    E, P, adj, grid = [], [], [], {}
    CELL = 3.0

    def gkey(p):
        return (int(math.floor(p[0] / CELL)), int(math.floor(p[1] / CELL)), int(math.floor(p[2] / CELL)))

    def add_atom(e, p, bonded):
        idx = len(E)
        E.append(e); P.append(p); adj.append(set(bonded))
        for b in bonded:
            adj[b].add(idx)
        grid.setdefault(gkey(p), []).append(idx)
        return idx

    def within3(a_nbrs, target):
        """Is target within 3 bonds, given the candidate atom's (placed) neighbours a_nbrs?"""
        frontier, seen = set(a_nbrs), set(a_nbrs)
        for _ in range(2):
            if target in seen:
                return True
            nf = set()
            for x in frontier:
                nf |= adj[x]
            frontier = nf - seen
            seen |= nf
        return target in seen

    def score(cand):
        """cand: list of (elem, pos, placed-neighbour ids, local-index-of-candidate-neighbours)."""
        worst = 9.0
        for ci, (e1, p1, pn, ln) in enumerate(cand):
            # neighbours of this candidate among placed atoms, including via candidate partners
            reach = set(pn)
            for l in ln:
                reach |= set(cand[l][2])
            kx, ky, kz = gkey(p1)
            for dx in (-1, 0, 1):
                for dy in (-1, 0, 1):
                    for dz in (-1, 0, 1):
                        for t in grid.get((kx + dx, ky + dy, kz + dz), ()):
                            v = dist(p1, P[t]) - lim(e1, E[t])
                            if v < worst and not within3(reach, t):
                                worst = v
        return worst

    B = [(0.0, 0.0, 0.0), (1.53, 0.0, 0.0)]
    B.append(add(B[1], (1.53 * math.cos(math.radians(66)), 1.53 * math.sin(math.radians(66)), 0.0)))
    first_virt = _place(B[2], B[1], B[0], 1.53, 114.0, 180.0)
    bid = [add_atom("C", B[0], []), None, None]
    bid[1] = add_atom("C", B[1], [bid[0]])
    bid[2] = add_atom("C", B[2], [bid[1]])
    S = {}
    snaps = {}
    backtracks = 0

    def rollback(to_j):
        n_e, n_b = snaps[to_j]
        for idx in range(n_e, len(E)):
            k_ = gkey(P[idx])
            grid[k_] = [t for t in grid[k_] if t < n_e]
        del E[n_e:]; del P[n_e:]; del adj[n_e:]
        for s_ in adj:
            s_.difference_update([t for t in s_ if t >= n_e])
        del B[n_b:]; del bid[n_b:]
        for k_ in [k_ for k_ in S if k_ >= to_j]:
            del S[k_]

    jj = 0
    while jj < nb:
        snaps[jj] = (len(E), len(B))
        prev = B[jj - 1] if jj > 0 else first_virt
        need_next = len(B) <= jj + 1
        best = None
        for trial in range(120):
            if need_next:
                t = (180.0 if (not curve or rnd.random() < 0.5) else rnd.choice([60.0, -60.0])) + rnd.uniform(-15, 15)
                nxt = _place(B[jj - 2], B[jj - 1], B[jj], 1.53, 114.0, t)
            else:
                nxt = B[jj + 1]
            phi = rnd.uniform(0, math.pi)
            sb = subs(jj, prev, B[jj], nxt, phi)
            # candidate list with bonding: ring C0 (ipso) bonded to backbone jj; ring cyclic; ring H m -> ring C m; H -> backbone jj
            cand = []
            if jj % 2 == 1:
                for m in range(6):
                    ln = [(m - 1) % 6, (m + 1) % 6]
                    cand.append(("C", sb[m][1], [bid[jj]] if m == 0 else [], ln))
                for m in range(1, 6):
                    cand.append(("H", sb[5 + m][1], [], [m]))
                cand.append(("H", sb[11][1], [bid[jj]], []))
            else:
                for e, pos, _ in sb:
                    cand.append(("H", pos, [bid[jj]], []))
            if need_next:
                cand.append(("C", nxt, [bid[jj]], []))
            sc = score(cand)
            if best is None or sc > best[0]:
                best = (sc, nxt, sb, cand)
            if sc >= 0:
                break
        sc, nxt, sb, cand = best
        if sc < -0.05 and jj >= 4 and backtracks < 400:
            backtracks += 1
            back = min(jj - 2, 2 + rnd.randrange(4))
            back -= back % 1
            rollback(jj - back)
            jj -= back
            continue
        ids = []
        for ci, (e, pos, pn, ln) in enumerate(cand):
            nbrs = list(pn) + [ids[l] for l in ln if l < len(ids)]
            ids.append(add_atom(e, pos, nbrs))
        if need_next:
            B.append(nxt)
            bid.append(ids[-1])
        S[jj] = sb
        jj += 1
    virt = B[nb] if len(B) > nb else _place(B[nb - 3], B[nb - 2], B[nb - 1], 1.53, 114.0, 180.0)
    B = B[:nb]
    atoms, backbone, stereo = [], [], []
    for k in range(nb):
        atoms.append({"e": "C", "p": B[k]})
        backbone.append(k)
    for k in range(nb):
        if k % 2 == 1:
            stereo.append(k)
        for e, p, tag in S[k]:
            if e == "C":
                atoms.append({"e": "C", "p": p, "ar": True})
            elif hydrogens:
                atoms.append({"e": "H", "p": p})
    if hydrogens:
        for k, ext in ((0, first_virt), (nb - 1, virt)):
            atoms.append({"e": "H", "p": add(B[k], mul(norm(sub(ext, B[k])), 1.09))})
    c = mul(tuple(map(sum, zip(*[a["p"] for a in atoms]))), 1 / len(atoms))
    for a in atoms:
        a["p"] = sub(a["p"], c)
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
