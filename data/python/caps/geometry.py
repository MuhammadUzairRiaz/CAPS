"""Small geometry helpers for steps (numpy, imported when used)."""


def dihedrals(pos):
    """Dihedral angles (degrees) along a chain of points: one per four consecutive points."""
    import numpy as np
    p = np.asarray(pos, dtype=float)
    b0, b1, b2 = p[1:-2] - p[:-3], p[2:-1] - p[1:-2], p[3:] - p[2:-1]
    n1, n2 = np.cross(b0, b1), np.cross(b1, b2)
    m1 = np.cross(n1, b1 / np.linalg.norm(b1, axis=1)[:, None])
    return np.degrees(np.arctan2(np.einsum("ij,ij->i", m1, n2), np.einsum("ij,ij->i", n1, n2)))


def angles(pos):
    """Bond angles (degrees) along a chain of points."""
    import numpy as np
    p = np.asarray(pos, dtype=float)
    u, v = p[:-2] - p[1:-1], p[2:] - p[1:-1]
    c = np.einsum("ij,ij->i", u, v) / (np.linalg.norm(u, axis=1) * np.linalg.norm(v, axis=1))
    return np.degrees(np.arccos(np.clip(c, -1, 1)))
