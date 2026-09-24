#!/usr/bin/env python3
"""LAMMPS data files written by CAPS reproduce CAPS's energy and forces in LAMMPS, for every force-field family.

For each case: CAPS assigns the force field (caps ff apply), writes the data file and the matching LAMMPS input
(-o, --lammps-input) and its own forces; LAMMPS runs that input (run 0); every energy term and every atomic force are
compared, and the six components of the virial tensor (LAMMPS compute pressure NULL virial). The cases cover class I (CVFF, OPLS-AA, GAFF), class II with class I torsions in one file (DL_FIELD's PCFF
and COMPASS: hybrid styles with skip lines in the class II sections), DREIDING (umbrella inversions), ionic crystals
(Buckingham pairs, periodic; Gasteiger charges, as DL_FIELD keeps ionic charges in its templates), a periodic polymer melt (tail corrections), and CHARMM-type force fields, whose separate
1-4 Lennard-Jones parameters LAMMPS cannot reproduce without switching (the writer refuses them; reported as such).

With --pme, CAPS uses particle-mesh Ewald on a fine grid (β from ewald-rtol 1e-7, spacing 0.5 Å, order 6) and LAMMPS
its Ewald sum (kspace_style ewald): both converge to the same electrostatics; only periodic cases are run.

usage: check_data_lammps.py [--only substring] [--keep DIR] [--pme]
Needs LMP (default ~/lammps/build-class2/lmp) with CLASS2, MOLECULE, EXTRA-MOLECULE, EXTRA-PAIR and MOFFF.
"""
import math, os, re, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE)
import template_pdb, ionic_pdb

CAPS = os.path.join(ROOT, "build", "cli", "caps")
LMP = os.environ.get("LMP", os.path.expanduser("~/lammps/build-class2/lmp"))
LIB = os.environ.get("DLFIELD", "~/project/dl_f_4.13") + "/lib"
FF = os.path.join(ROOT, "data", "forcefields")
arg = lambda k, d: sys.argv[sys.argv.index(k) + 1] if k in sys.argv else d
only = arg("--only", "")
PME = "--pme" in sys.argv
work = arg("--keep", "") or tempfile.mkdtemp()
os.makedirs(work, exist_ok=True)

# (label, structure source, force field, charges, typing: "rules" or "keys")
CASES = [
    ("PCFF favipiravir (class II + class I)", ("template", "PCFF", "favipiravir"), "pcff-dlfield", "types", "rules"),
    ("PCFF acetylcholine (N+)", ("template", "PCFF", "acetylcholine"), "pcff-dlfield", "types", "rules"),
    ("PCFF dimethyl carbonate", ("template", "PCFF", "dimethyl_carbonate"), "pcff-dlfield", "types", "rules"),
    ("COMPASS (DL_FIELD) first template", ("template-first", "COMPASS", ""), "compass-dlfield", "types", "keys"),
    ("CVFF phenol (cvff impropers)", ("template", "CVFF", "phenol"), "cvff-dlfield", "types", "rules"),
    ("CVFF ethyl acetate", ("template", "CVFF", "ethylacetate"), "cvff-dlfield", "types", "rules"),
    ("OPLS-AA methyl vinyl ketone", ("template", "OPLS2005", "methyl_vinyl_ketone"), "opls2005-dlfield", "gasteiger", "rules"),
    ("OPLS-AA phthalimide", ("template", "OPLS2005", "phthalimide"), "opls2005-dlfield", "gasteiger", "rules"),
    ("GAFF toluene", ("template", "AMBER16_gaff", "toluene"), "gaff-amber16-dlfield", "gasteiger", "rules"),
    ("GAFF2 phenol", ("template", "AMBER25_gaff", "phenol"), "gaff-amber25-dlfield", "gasteiger", "rules"),
    ("DREIDING first template (umbrella inversions)", ("template-first", "DREIDING", ""), "dreiding-dlfield", "gasteiger", "keys"),
    ("Ionic halides (Buckingham, periodic)", ("ionic-first", "INORGANIC_binary_halides", ""), "inorganic-binary-halides-dlfield", "gasteiger", "keys"),
    ("Binary oxides (Buckingham, periodic)", ("ionic-first", "INORGANIC_binary_oxides", ""), "inorganic-binary-oxides-dlfield", "gasteiger", "keys"),
    ("Polystyrene melt, GAFF2 (periodic, 1300 atoms)", ("file", os.path.join(ROOT, "samples", "ps_melt.data")), "gaff-amber25-dlfield", "gasteiger", "rules"),
    ("Polystyrene melt, PCFF (periodic, class II + class I)", ("file", os.path.join(ROOT, "samples", "ps_melt.data")), "pcff-dlfield", "types", "rules"),
    # every bonded kind hybrid: COMPASS (class II with all cross terms) plus an overlay turning one angle, one torsion
    # and one improper type into class I forms, so each class II section carries skip lines. (The improper is cvff,
    # K[1 − cos 2φ]: smooth. A harmonic improper with χ0 = 0 on a class II centre-second quadruple sits at χ ≈ 180°,
    # a cusp of K(χ − χ0)² where the two programs take different one-sided slopes.)
    ("COMPASS polystyrene + class I overlay (hybrid in every kind)", ("compass-ps",), "compass-published-moltemplate", "types", "keys"),
    ("CGenFF toluene (separate 1-4 LJ)", ("template", "CHARMM36_cgenff", "toluene"), "cgenff-dlfield", "gasteiger", "rules"),
    # UFF (every element): Fourier and periodic angles (linear, trigonal, square planar, octahedral, trigonal
    # bipyramid), sp2 and pyramidal-P inversions (improper fourier), group-16 torsions, full 1-4 van der Waals
    ("UFF mixed elements (P, S, Si, Pt, F, Cl)", ("smiles", "CC#CC(=O)Oc1ccc(cc1)P(C)C.F[S](F)(F)(F)(F)F.N[Pt](N)(Cl)Cl."
                                               "FP(F)(F)(F)F.C[Si](C)(C)O[Si](C)(C)C.CSSC"), "uff", "types", "rules"),
    ("Polystyrene melt, UFF (periodic, 1300 atoms)", ("file", os.path.join(ROOT, "samples", "ps_melt.data")), "uff", "types", "rules"),
]


OVERLAY = """{"format": "caps-forcefield", "version": 1, "name": "class I test overlay", "atom_types": [],
 "angles": [{"name": "test c4-c4-c4", "match": ["c4~*", "c4~*", "c4~*"], "style": "harmonic", "params": [60.0, 112.0]}],
 "dihedrals": [{"name": "test c4-c4-c4-c4", "match": ["c4~*", "c4~*", "c4~*", "c4~*"], "style": "fourier", "params": [2, 0.4, 1, 0.0, 0.2, 3, 0.0]}],
 "impropers": [{"name": "test ipso", "match": ["c3a~*", "c3a~*", "c3a~*", "c4~*"], "style": "cvff", "params": [8.0, -1, 2]}]}
"""


def ff_file(fid):
    if fid == "uff":   # built in: typed from elements and bonds
        return "uff"
    p = os.path.join(FF, fid + ".json")
    return p if os.path.exists(p) else None


def structure(src, base):
    """(structure path, types file or None)"""
    kind = src[0]
    if kind == "file":
        return src[1], None
    if kind == "smiles":   # built and cleaned up by CAPS with UFF
        m = os.path.join(work, base + ".mol2")
        subprocess.run([CAPS, "build", src[1], "--ff", "uff", "-o", m], capture_output=True, check=True)
        return m, None
    if kind == "compass-ps":   # polystyrene melt with COMPASS types (from the GAFF2 typing: c3 → c4, ca → c3a, H → h1)
        ps = os.path.join(ROOT, "samples", "ps_melt.data")
        gt = os.path.join(work, base + ".gaff")
        subprocess.run([CAPS, "ff", "type", ps, "--ff", ff_file("gaff-amber25-dlfield"), "-o", gt], capture_output=True, check=True)
        tfile = os.path.join(work, base + ".types")
        with open(tfile, "w") as f:
            for t in open(gt):
                f.write({"c3": "c4", "ca": "c3a", "hc": "h1", "ha": "h1"}[t.strip()] + "\n")
        return ps, tfile
    lib = src[1]
    sf = os.path.join(LIB, lib + ".sf")
    if not os.path.exists(sf):
        return None, None
    if kind == "template-first":
        names = template_pdb.candidates(sf, 30)
        if not names:
            return None, None
        name = names[0]
    elif kind == "ionic-first":   # rigid-ion templates (no core-shell pairs, whose springs this check does not cover)
        names = [n for n, m in template_pdb.parse_sf(sf)[2].items()
                 if ionic_pdb.ionic(m) and not any(re.search(r"[A-Za-z][0-9]*[cs]$", a) for a, _ in m["atoms"])]
        if not names:
            return None, None
        name = names[0]
    else:
        name = src[2]
        mols0 = template_pdb.parse_sf(sf)[2]
        if name not in mols0:   # this library names it differently: its first small template instead
            names = template_pdb.candidates(sf, 30)
            if not names:
                return None, None
            name = names[0]
    pdb = os.path.join(work, base + ".pdb")
    if kind == "ionic-first":
        ionic_pdb.build(sf, name, pdb, units=4)
    else:
        template_pdb.build(sf, name, pdb)
    types, keys, mols = template_pdb.parse_sf(sf)
    key_of = dict(template_pdb.ATOM_KEY)
    tfile = os.path.join(work, base + ".types")
    atoms = mols[name]["atoms"]
    if kind == "ionic-first":   # the builder repeats the formula unit
        n = sum(1 for l in open(pdb) if l.startswith(("ATOM", "HETATM")))
        atoms = [atoms[i % len(atoms)] for i in range(n)]
    with open(tfile, "w") as f:
        for _, t in atoms:
            f.write(key_of.get(t, t) + "\n")
    return pdb, tfile


def lammps(infile, dump):
    txt = open(infile).read().replace("run 0", f"run 0\nwrite_dump all custom {dump} id fx fy fz modify sort id format float %.10f")
    # virial pressure tensor (no kinetic part) and the volume, for the virial tensor W = P V / nktv2p
    txt = re.sub(r"thermo_style custom (.*)", r"compute pv all pressure NULL virial\nthermo_style custom \1 vol c_pv[1] c_pv[2] c_pv[3] c_pv[4] c_pv[5] c_pv[6]", txt)
    with open(infile, "w") as f:
        f.write(txt)
    r = subprocess.run([LMP, "-in", os.path.basename(infile), "-log", "none"], capture_output=True, text=True, cwd=os.path.dirname(infile))
    if r.returncode:
        raise RuntimeError((r.stdout + r.stderr).strip().splitlines()[-1][:200])
    lines = r.stdout.splitlines()
    k = next(i for i, l in enumerate(lines) if l.split()[:2] == ["Step", "PotEng"])
    head, vals = lines[k].split(), [float(x) for x in lines[k + 1].split()]
    e = dict(zip(head, vals))
    f = {}
    for l in open(dump).read().split("ITEM: ATOMS")[1].splitlines()[1:]:
        w = l.split()
        if len(w) == 4:
            f[int(w[0])] = tuple(map(float, w[1:]))
    return e, f


rows, fails = [], 0
for label, src, fid, charges, typing in CASES:
    if only and only.lower() not in label.lower():
        continue
    base = re.sub(r"[^a-z0-9]+", "_", label.lower()).strip("_")[:40]
    ffj = ff_file(fid)
    s, tfile = structure(src, base)
    if not ffj or not s:
        rows.append((label, "skipped (source not available)", "", ""))
        continue
    d = os.path.join(work, base)
    os.makedirs(d, exist_ok=True)
    cmd = [CAPS, "ff", "apply", s, "--ff", ffj, "--charges", charges, "-o", os.path.join(d, "case.data"),
           "--lammps-input", os.path.join(d, "case.in"), "--forces", os.path.join(d, "caps_f.txt")]
    if typing == "keys" and tfile:
        cmd += ["--types", tfile]
    if PME:
        if src[0] not in ("file", "ionic-first"):
            continue
        cmd += ["--pme", "--ewald-rtol", "1e-7", "--pme-spacing", "0.5", "--pme-order", "6"]
    if src[0] == "compass-ps":
        ov = os.path.join(work, base + ".overlay.json")
        with open(ov, "w") as f:
            f.write(OVERLAY)
        cmd += ["--overlay", ov]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode or "missing parameters" in r.stdout:
        why = (r.stderr.strip().splitlines() or r.stdout.strip().splitlines() or ["?"])[-1]
        expected = "1-4 Lennard-Jones" in why
        rows.append((label, ("refused, as intended: " if expected else "CAPS failed: ") + why[:150], "", ""))
        fails += 0 if expected else 1
        continue
    m = re.search(r"energy \(kcal/mol\): bond (\S+)  angle (\S+)  dihedral (\S+)  improper (\S+)  vdW (\S+)  Coulomb (\S+)", r.stdout)
    ce = dict(zip(["bond", "angle", "dihedral", "improper", "vdw", "coulomb"], map(float, m.groups())))
    cf = {int(w[0]): tuple(map(float, w[1:4])) for w in (l.split() for l in open(os.path.join(d, "caps_f.txt"))) if len(w) == 4}
    try:
        le, lf = lammps(os.path.join(d, "case.in"), os.path.join(d, "f.dump"))
    except Exception as ex:
        rows.append((label, "LAMMPS failed: " + str(ex), "", ""))
        fails += 1
        continue
    lm = {"bond": le["E_bond"], "angle": le["E_angle"], "dihedral": le["E_dihed"], "improper": le["E_impro"], "vdw": le["E_vdwl"],
          "coulomb": le["E_coul"] + le["E_long"]}
    # with PME both codes reach the Ewald limit only to their discretisation and the real-space erfc approximation
    # (Abramowitz–Stegun, 1.5e-7 per pair, as LAMMPS): the Coulomb term is compared to 1e-3 kcal/mol absolute there
    de = max(abs(ce[k] - lm[k]) / (max(1.0, abs(ce[k])) if not (PME and k == "coulomb") else 100.0) for k in ce)
    mw = re.search(r"virial tensor \(kcal/mol\): xx (\S+)  yy (\S+)  zz (\S+)  xy (\S+)  xz (\S+)  yz (\S+)", r.stdout)
    cw = list(map(float, mw.groups()))
    lw = [le[f"c_pv[{k}]"] * le["Volume"] / 68568.415 for k in range(1, 7)]
    scale = max(1.0, max(abs(w) for w in cw))
    dw = max(abs(a - b) for a, b in zip(cw, lw)) / scale
    df = max(math.dist(cf[i], lf[i]) for i in cf)
    styles = " · ".join(l.strip() for l in open(os.path.join(d, "case.in")) if re.match(r"(bond|angle|dihedral|improper|pair)_style", l))
    ok = de < 1e-5 and df < 1e-3 and dw < 1e-5
    fails += 0 if ok else 1
    rows.append((label, f"{'ok' if ok else 'DIFFERS'} · energy terms {de:.1e} (relative) · forces {df:.1e} kcal/mol/Å · virial tensor {dw:.1e} (relative)", styles,
                 " ".join(f"{k} {ce[k]:.4f}/{lm[k]:.4f}" for k in ce if abs(ce[k]) > 0 or abs(lm[k]) > 0)))

for label, res, styles, terms in rows:
    print(f"{label}\n   {res}")
    if styles:
        print(f"   styles: {styles}")
    if res.startswith("DIFFERS"):
        print(f"   CAPS/LAMMPS: {terms}")
print(f"\n{sum(1 for r in rows if r[1].startswith('ok'))} of {len(rows)} cases agree with LAMMPS"
      f"{'' if not fails else f'; {fails} failed'} (work: {work})")
sys.exit(1 if fails else 0)
