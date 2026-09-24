"""Tiny Z-matrix builder for validation molecules (no external chemistry toolkit needed).

rows: (label, type, charge, a, r, b, theta, c, phi) with a/b/c earlier labels (None for the first atoms).
Writes a moltemplate .lt molecule with its bond list (bonds are the a-references).
"""
import math

def _place(pa, pb, pc, r, theta, phi):
    th, ph = math.radians(theta), math.radians(phi)
    bc = [pa[i] - pb[i] for i in range(3)]
    n = math.sqrt(sum(x * x for x in bc)); bc = [x / n for x in bc]
    ab = [pb[i] - pc[i] for i in range(3)]
    nv = [ab[1] * bc[2] - ab[2] * bc[1], ab[2] * bc[0] - ab[0] * bc[2], ab[0] * bc[1] - ab[1] * bc[0]]
    n = math.sqrt(sum(x * x for x in nv)); nv = [x / n for x in nv]
    m = [nv[1] * bc[2] - nv[2] * bc[1], nv[2] * bc[0] - nv[0] * bc[2], nv[0] * bc[1] - nv[1] * bc[0]]
    d2 = [-r * math.cos(th), r * math.sin(th) * math.cos(ph), r * math.sin(th) * math.sin(ph)]
    return [pa[i] + d2[0] * bc[i] + d2[1] * m[i] + d2[2] * nv[i] for i in range(3)]

def build(rows):
    pos, out = {}, []
    for k, (lab, ty, q, a, r, b, th, c, ph) in enumerate(rows):
        if k == 0: p = [0.0, 0.0, 0.0]
        elif k == 1: p = [r, 0.0, 0.0]
        elif k == 2:
            pa, pb = pos[a], pos[b]
            d = [pb[i] - pa[i] for i in range(3)]; n = math.sqrt(sum(x * x for x in d)); d = [x / n for x in d]
            t = math.radians(th)
            p = [pa[0] + r * (d[0] * math.cos(t) - d[1] * math.sin(t)), pa[1] + r * (d[0] * math.sin(t) + d[1] * math.cos(t)), 0.0]
        else: p = _place(pos[a], pos[b], pos[c], r, th, ph)
        pos[lab] = p
        out.append((lab, ty, q, p, a))
    return out

def methyl(c, x, y, h="h1", r=1.09, prefix=None, q=0.0):
    """three H on carbon c, staggered about the x-c bond (y fixes the dihedral origin)."""
    prefix = prefix or c + "H"
    return [(f"{prefix}{i+1}", h, q, c, r, x, 110.5, y, 60 + 120 * i) for i in range(3)]

def write_lt(atoms, ffname, ffclass, molname, path, box=20.0, extra_bonds=()):
    with open(path, "w") as f:
        f.write(f'import "{ffname}.lt"\n{molname} inherits {ffclass} {{\n  write("Data Atoms") {{\n')
        for lab, ty, q, p, a in atoms:
            f.write(f"    $atom:{lab} $mol:. @atom:{ty} {q:.4f} {p[0]:.5f} {p[1]:.5f} {p[2]:.5f}\n")
        f.write('  }\n  write("Data Bond List") {\n')
        for k, (lab, ty, q, p, a) in enumerate(atoms):
            if a: f.write(f"    $bond:b{k} $atom:{a} $atom:{lab}\n")
        for k, (a, b) in enumerate(extra_bonds):
            f.write(f"    $bond:x{k} $atom:{a} $atom:{b}\n")
        f.write("  }\n}\n")
        f.write(f'write_once("Data Boundary") {{\n  {-box} {box} xlo xhi\n  {-box} {box} ylo yhi\n  {-box} {box} zlo zhi\n}}\n')
        f.write(f"m = new {molname}\n")

def element_of(ty):
    """Element from a force-field type name of the validation molecules (h1, c3a, o2e, si4c, ...)."""
    t = ty.lower()
    for e in ("si", "cl", "br"):
        if t.startswith(e): return e.capitalize()
    return t[0].upper()

def write_xyz(atoms, path, title="CAPS validation molecule"):
    with open(path, "w") as f:
        f.write(f"# {title}\n")
        for lab, ty, q, p, a in atoms:
            f.write(f"{element_of(ty):3s} {p[0]:10.5f} {p[1]:10.5f} {p[2]:10.5f}\n")

def write_mol2(atoms, types, bonds, path, name="CAPS"):
    """Tripos mol2 with force-field atom types (DL_FIELD reads the type column). types: label -> type."""
    with open(path, "w") as f:
        f.write(f"@<TRIPOS>MOLECULE\n{name}\n{len(atoms):5d} {len(bonds):5d}     1     0     0\nSMALL\nbcc\n\n\n@<TRIPOS>ATOM\n")
        idx, count = {}, {}
        for k, (lab, ty, q, p, a) in enumerate(atoms):
            idx[lab] = k + 1
            el = element_of(ty)
            count[el] = count.get(el, 0) + 1
            name = f"{el}{count[el]}"   # DL_FIELD takes the element from the atom name
            f.write(f"{k+1:7d} {name:<8s} {p[0]:10.4f} {p[1]:10.4f} {p[2]:10.4f} {types[lab]:<8s} 1 MOL  0.000000\n")
        f.write("@<TRIPOS>BOND\n")
        for k, (a, b) in enumerate(bonds):
            f.write(f"{k+1:6d} {idx[a]:5d} {idx[b]:5d} 1\n")

def bonds_of(atoms, extra=()):
    return [(a, lab) for lab, ty, q, p, a in atoms if a] + list(extra)
