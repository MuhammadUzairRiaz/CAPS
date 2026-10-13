"""A periodic box of a molecular mixture for the export checks (check_data_lammps.py, check_gromacs.py): each SMILES built
by CAPS (caps build, cleaned up with the force field given, else UFF), the copies packed by CAPS Pack 2.6 Å apart (so no
contact is read as a bond) and 1 Å from the faces of a cube of edge L, written as extended XYZ with the cell (bonds are
perceived from the positions)."""
import os, subprocess


def build(caps, parts, edge, out, ff=None, seed=1, tolerance=2.6):
    """parts: [(smiles, count)]; returns the .xyz path"""
    base = os.path.splitext(out)[0]
    inp = base + ".inp"
    with open(inp, "w") as f:
        f.write(f"tolerance {tolerance}\nseed {seed}\nfiletype xyz\noutput {out}\n")
        for k, (smi, n) in enumerate(parts):
            one = f"{base}.part{k}.xyz"
            subprocess.run([caps, "build", smi, "-o", one] + (["--ff", ff] if ff else []), capture_output=True, check=True)
            f.write(f"structure {one}\n  number {n}\n  inside box 1. 1. 1. {edge - 1} {edge - 1} {edge - 1}\nend structure\n")
    subprocess.run([caps, "pack", inp, "--quiet"], capture_output=True, check=True)
    lines = open(out).read().splitlines()
    lines[1] = f'Lattice="{edge} 0 0 0 {edge} 0 0 0 {edge}" Properties=species:S:1:pos:R:3'
    with open(out, "w") as f:
        f.write("\n".join(lines) + "\n")
    return out
