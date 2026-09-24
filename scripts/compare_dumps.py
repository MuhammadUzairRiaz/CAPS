#!/usr/bin/env python3
"""Largest per-atom position difference between two LAMMPS dumps (id xu yu zu), frame by frame. Exit 1 above a tolerance."""
import sys
def frames(p):
    out={}
    for f in open(p).read().split("ITEM: TIMESTEP\n")[1:]:
        lines=f.splitlines(); step=int(lines[0]); n=int(lines[2])
        hdr=lines[7].split()[2:]
        ix=[hdr.index(c) for c in ("id","xu","yu","zu")]
        out[step]={}
        for l in lines[8:8+n]:
            t=l.split(); out[step][int(t[ix[0]])]=[float(t[k]) for k in ix[1:]]
    return out
a=frames(sys.argv[1]); b=frames(sys.argv[2])
tol=float(sys.argv[3]) if len(sys.argv)>3 else 1e-3
worst=0.0
for st in sorted(set(a)&set(b)):
    d=max(sum((p-q)**2 for p,q in zip(a[st][i],b[st][i]))**0.5 for i in a[st])
    worst=max(worst,d)
    print(f"step {st}: largest position difference {d:.2e} Å")
sys.exit(0 if worst<tol else 1)
