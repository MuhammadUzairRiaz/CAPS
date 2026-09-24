"""Shared display-style pieces: all-atom PS cells, the viewport style switcher and a zoom inset."""
import math, random
from lib import *
from mols import polystyrene
from screen_row16b import ps_bonds

STYLES = ["All atoms", "No H", "Backbone"]


def _rotm(rnd):
    a, b, c = (rnd.uniform(0, 2 * math.pi) for _ in range(3))
    def R(p):
        return rot(p, a, b, c)
    return R


def ps_cell(n_chains, dp, L, seed=1, min_d=3.0, tries=3000):
    """Rigid all-atom PS chains (with H) placed so no heavy atoms of different chains come closer than
    min_d under periodic boundaries (minimum image, cubic box L). Returns [(atoms, bonds, backbone)]."""
    rnd = random.Random(seed)
    cs = min_d
    nc = max(1, int(L // cs))
    cs = L / nc
    grid = {}

    def key(p):
        return tuple(int((p[k] % L) // cs) % nc for k in range(3))

    def clash(p):
        kx, ky, kz = key(p)
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for q in grid.get(((kx + dx) % nc, (ky + dy) % nc, (kz + dz) % nc), ()):
                        d2 = 0.0
                        for k in range(3):
                            x = q[k] - p[k]
                            x -= L * round(x / L)
                            d2 += x * x
                        if d2 < min_d * min_d:
                            return True
        return False

    placed = []
    for k in range(n_chains):
        atoms0, backbone, stereo = polystyrene(dp, "atactic", seed=seed * 100 + k)
        bonds = ps_bonds(atoms0, backbone, stereo)
        c0 = mul(tuple(map(sum, zip(*[a["p"] for a in atoms0]))), 1 / len(atoms0))
        ok = None
        for _ in range(tries):
            R = _rotm(rnd)
            cen = (rnd.uniform(0, L), rnd.uniform(0, L), rnd.uniform(0, L))
            pts = [add(R(sub(a["p"], c0)), cen) for a in atoms0]
            hv = [p for p, a in zip(pts, atoms0) if a["e"] != "H"]
            if not any(clash(p) for p in hv):
                ok = pts
                break
        if ok is None:
            raise RuntimeError(f"could not place chain {k + 1} without overlap; lower n_chains or raise L")
        atoms = [dict(a, p=p) for a, p in zip(atoms0, ok)]
        placed.append((atoms, bonds, backbone))
        for p, a in zip(ok, atoms0):
            if a["e"] != "H":
                grid.setdefault(key(p), []).append(p)
    return placed


def render_cell(sc, chains, style, colours=None, tube_w=0.55):
    """Add a PS cell to a Scene in one of the three display styles."""
    for i, (atoms, bonds, backbone) in enumerate(chains):
        col_ = colours[i % len(colours)] if colours else None
        if style == "Backbone":
            sc.add_tube([atoms[j]["p"] for j in backbone], col_ or ACC, tube_w)
            continue
        keep = [j for j, a in enumerate(atoms) if style == "All atoms" or a["e"] != "H"]
        idx = {j: n for n, j in enumerate(keep)}
        sel = [dict(atoms[j], c=col_) if (col_ and atoms[j]["e"] != "H") else atoms[j] for j in keep]
        sc.add_atoms(sel, [(idx[a], idx[b]) for a, b in bonds if a in idx and b in idx])


def count(chains):
    n = sum(len(a) for a, _, _ in chains)
    h = sum(1 for a, _, _ in chains for x in a if x["e"] == "H")
    return n, h


def style_bar(active, note, pos="right: 12px; top: 10px"):
    return (f'<div style="position: absolute; {pos}; display: flex; flex-direction: column; align-items: flex-end; gap: 6px">'
            f'<div style="display: flex; align-items: center; gap: 8px"><span style="font-size: 11.5px; color: {MUTED}">Display</span>{seg(STYLES, active)}</div>'
            f'{chip(note, TEXT, "#16191Ccc")}</div>')


def zoom_inset(sid, w=230, h=170, dp=3, seed=7, pos="left: 12px; bottom: 12px", label="Zoom · all atoms incl. H"):
    atoms, backbone, stereo = polystyrene(dp, "atactic", seed=seed)
    bonds = ps_bonds(atoms, backbone, stereo)
    sc = Scene(sid, w, h - 26, yaw=-0.1, pitch=0.5, persp=0.15, fog=0.2, atom_k=1.05)
    sc.add_atoms(atoms, bonds)
    nH = sum(1 for a in atoms if a["e"] == "H")
    return (f'<div style="position: absolute; {pos}; width: {w}px; background: {BG0}; border: 1px solid {ACC}; border-radius: 8px; overflow: hidden; box-shadow: 0 8px 24px #0009">'
            f'<div style="height: 26px; display: flex; align-items: center; justify-content: space-between; gap: 6px; padding: 0 8px; background: {BG1}; border-bottom: 1px solid {LINE}; font-size: 10.5px; white-space: nowrap">'
            f'<span>{label}</span><span style="font-family: {MONO}; color: {DIM}">{len(atoms)} · {nH} H</span></div>{sc.svg()}</div>')
