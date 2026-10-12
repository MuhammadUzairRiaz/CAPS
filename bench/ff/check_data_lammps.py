#!/usr/bin/env python3
"""LAMMPS data files written by CAPS reproduce CAPS's energy and forces in LAMMPS, for every force-field family.

For each case: CAPS assigns the force field (caps ff apply), writes the data file and the matching LAMMPS input
(-o, --lammps-input) and its own forces; LAMMPS runs that input (run 0); every energy term and every atomic force are
compared, and the six components of the virial tensor (LAMMPS compute pressure NULL virial). The cases cover class I (CVFF, OPLS-AA, GAFF), class II with class I torsions in one file (DL_FIELD's PCFF
and COMPASS: hybrid styles with skip lines in the class II sections), DREIDING (umbrella inversions), ionic crystals
(Buckingham pairs, periodic; QEq charges, as DL_FIELD keeps ionic charges in its templates and Gasteiger–Marsili has no parameters for the metals), a periodic polymer melt (tail corrections), mW water (Stillinger–Weber three-body term, pair_style sw), and CHARMM-type force fields, whose separate
1-4 Lennard-Jones parameters LAMMPS cannot reproduce without switching (the writer refuses them; reported as such).

With --pme, CAPS uses particle-mesh Ewald on a fine grid (β from ewald-rtol 1e-9, spacing 0.3 Å, order 8) and LAMMPS
its Ewald sum (kspace_style ewald): both converge to the same electrostatics; only periodic cases are run.

usage: check_data_lammps.py [--only substring] [--keep DIR] [--pme]
Needs LMP (default ~/lammps/build-class2/lmp) with CLASS2, MOLECULE, EXTRA-MOLECULE, EXTRA-PAIR, MANYBODY and MOFFF.
"""
import collections, math, os, re, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, HERE)
import template_pdb, ionic_pdb

CAPS = os.path.join(ROOT, "build", "cli", "caps")
LMP = os.environ.get("LMP", os.path.expanduser("~/lammps/build-class2/lmp"))
LIB = os.environ.get("DLFIELD", os.path.expanduser("~/project/dl_f_4.13")) + "/lib"
FF = os.path.join(ROOT, "data", "forcefields")
arg = lambda k, d: sys.argv[sys.argv.index(k) + 1] if k in sys.argv else d
only = arg("--only", "")
PME = "--pme" in sys.argv
# --native: the force field's own LAMMPS styles (OPLS dihedrals, long-range Coulomb by PPPM ...), as exports write them;
# --hybrid with it: every style hybrid. The bonded and van der Waals terms must equal CAPS's; Coulomb is PPPM there,
# CAPS's damped shifted force here (different methods), so it and the forces are not compared.
NATIVE = "--native" in sys.argv
HYBRID = "--hybrid" in sys.argv
work = arg("--keep", "") or tempfile.mkdtemp()
os.makedirs(work, exist_ok=True)

# (label, structure source, force field, charges, typing: "rules" or "keys")
CASES = [
    ("PCFF favipiravir (class II + class I)", ("template", "PCFF", "favipiravir"), "pcff", "types", "rules"),
    ("PCFF acetylcholine (N+)", ("template", "PCFF", "acetylcholine"), "pcff", "types", "rules"),
    ("PCFF dimethyl carbonate", ("template", "PCFF", "dimethyl_carbonate"), "pcff", "types", "rules"),
    ("COMPASS (DL_FIELD) first template", ("template-first", "COMPASS", ""), "compass", "types", "keys"),
    ("CVFF phenol (cvff impropers)", ("template", "CVFF", "phenol"), "cvff", "types", "rules"),
    ("CVFF ethyl acetate", ("template", "CVFF", "ethylacetate"), "cvff", "types", "rules"),
    ("OPLS-AA methyl vinyl ketone", ("template", "OPLS2005", "methyl_vinyl_ketone"), "oplsaa2024", "auto", "rules"),
    ("OPLS-AA phthalimide", ("template", "OPLS2005", "phthalimide"), "oplsaa2024", "auto", "rules"),
    ("GAFF toluene", ("template", "AMBER16_gaff", "toluene"), "gaff-1.81-ambertools", "gasteiger", "rules"),
    ("GAFF2 phenol", ("template", "AMBER25_gaff", "phenol"), "gaff2-2.2.20-ambertools", "gasteiger", "rules"),
    ("DREIDING first template (umbrella inversions)", ("template-first", "DREIDING", ""), "dreiding", "gasteiger", "keys"),
    # DREIDING's own torsion rules (Mayo 1990, cases a-j) where the file lists none: C=C double and conjugated single
    # bonds (polybutadiene, butadiene), sp3 next to sp2 (b / j), biphenyl's ring link (f), the aromatic ring (d)
    ("DREIDING ethanol liquid: hydrogen bonds (periodic, 270 molecules)", ("solvent-box", "ethanol", 26), "dreiding1990", "gasteiger", "rules"),
    ("Polystyrene melt, DREIDING (periodic; sp3-aromatic torsions by DREIDING's rules)", ("file", os.path.join(ROOT, "samples", "ps_melt.data")), "dreiding1990", "gasteiger", "rules"),
    ("DREIDING rules: cis-polybutadiene fragment, butadiene, biphenyl, methyl acetate", ("smiles", "C/C=C\\CC/C=C\\CC/C=C\\C.C=CC=C.c1ccccc1-c1ccccc1.CC(=O)OC"), "dreiding1990", "gasteiger", "rules"),
    ("Ionic halides (Buckingham, periodic)", ("ionic-first", "INORGANIC_binary_halides", ""), "inorganic-binary-halides", "qeq", "keys"),
    # a shell model (core and shell on each ion, a spring between them): CORESHELL's buck/coul/long/cs (born/coul/dsf/cs
    # in CAPS's own styles, which LAMMPS gives no tail term: compared without one). LAMMPS's /cs styles evaluate each
    # core-shell pair at r = 1e-10 Å, where 1 - erfc cancels: its Coulomb carries ~1e-5 of that noise; CAPS the exact limit
    ("SrTiO3 perovskite, core-shell (CORESHELL styles)", ("crystal", ["--group", "P m -3 m", "--cell", "3.905,3.905,3.905", "--sites",
     "Sr1 Sr 0 0 0; Ti1 Ti 0.5 0.5 0.5; O1 O 0.5 0.5 0", "--supercell", "3,3,3"]), "inorganic-ternary-oxides", "types", "rules", ["--no-tail"]),
    ("Binary oxides (Buckingham, periodic)", ("ionic-first", "INORGANIC_binary_oxides", ""), "inorganic-binary-oxides", "qeq", "keys"),
    ("Alkali halides (Catlow, Diller & Norgett 1977, set 1): NaCl rock salt, rigid ions", ("crystal", ["--group", "F m -3 m", "--cell", "5.64,5.64,5.64",
     "--sites", "Na1 Na 0 0 0; Cl1 Cl 0.5 0.5 0.5", "--supercell", "3,3,3"]), "catlow1977-halides", "types", "rules", ["--no-tail"]),
    ("Alkali halides (Catlow, Diller & Norgett 1977, set 1): LiF rock salt, rigid ions", ("crystal", ["--group", "F m -3 m", "--cell", "4.03,4.03,4.03",
     "--sites", "Li1 Li 0 0 0; F1 F 0.5 0.5 0.5", "--supercell", "5,5,5"]), "catlow1977-halides", "types", "rules", ["--no-tail"]),
    ("Lewis & Catlow 1985 (Table 3): MgO rigid Mg, O shell", ("crystal", ["--group", "F m -3 m", "--cell", "4.212,4.212,4.212", "--sites", "Mg1 Mg 0 0 0; O1 O 0.5 0.5 0.5",
     "--supercell", "3,3,3"]), "lewis-catlow1985-oxides", "types", "rules", ["--no-tail"]),
    ("Lewis & Catlow 1985 (Table 3): CaO, cation and anion shells", ("crystal", ["--group", "F m -3 m", "--cell", "4.81,4.81,4.81", "--sites", "Ca1 Ca 0 0 0; O1 O 0.5 0.5 0.5",
     "--supercell", "3,3,3"]), "lewis-catlow1985-oxides", "types", "rules", ["--no-tail"]),
    ("Catlow 1977 UO2 (Table 3, potential 1): fluorite, U and O shells", ("crystal", ["--group", "F m -3 m", "--cell", "5.47,5.47,5.47", "--sites",
     "U1 U 0 0 0; O1 O 0.25 0.25 0.25", "--supercell", "3,3,3"]), "catlow1977-uo2", "types", "rules", ["--no-tail"]),
    ("Matsui & Akaogi 1991 (Table 1): rutile, rigid ions", ("crystal", ["--group", "P 42/m n m", "--cell", "4.594,4.594,2.959", "--sites",
     "Ti1 Ti 0 0 0; O1 O 0.3053 0.3053 0", "--supercell", "3,3,5"]), "matsui1991-tio2", "types", "rules", ["--no-tail"]),
    ("Sanders, Leslie & Catlow 1984 (Table 1): alpha-quartz, O shells, O-Si-O bends by distance", ("file", os.path.join(ROOT, "data", "crystals",
     "alpha-quartz.cif")), "sanders1984-sio2", "types", "rules", ["--no-tail"]),
    # Xiao, Edwards & Grater 2011 CaCO3: carbonate bends and plane terms by distance, no C-O bond term, no exclusions
    ("Xiao 2011 CaCO3 (Table 1): calcite", ("file", os.path.join(ROOT, "data", "crystals", "calcite.cif")), "xiao2011-caco3", "types", "rules", ["--no-tail"]),
    ("Xiao 2011 CaCO3 (Table 1): aragonite (De Villiers 1971)", ("crystal", ["--group", "P m c n", "--cell", "4.9614,7.9671,5.7404", "--sites",
     "Ca1 Ca 0.25 0.41500 0.75970; C1 C 0.25 0.76220 -0.0852; O1 O 0.25 0.92250 -0.0961; O2 O 0.4736 0.68065 -0.08620"]), "xiao2011-caco3", "types", "rules", ["--no-tail"]),
    ("Gunn 2012 (Table 1): Gd2Ti2O7 pyrochlore", ("crystal", ["--group", "F d -3 m:2", "--cell", "10.185,10.185,10.185", "--sites",
     "Gd1 Gd 0.5 0.5 0.5; Ti1 Ti 0 0 0; O1 O 0.3300 0.125 0.125; O2 O 0.375 0.375 0.375"]), "gunn2012-pyrochlore", "types", "rules", ["--no-tail"]),
    ("Gunn 2012 (Table 1): Gd2Zr2O7 pyrochlore", ("crystal", ["--group", "F d -3 m:2", "--cell", "10.53,10.53,10.53", "--sites",
     "Gd1 Gd 0.5 0.5 0.5; Zr1 Zr 0 0 0; O1 O 0.3450 0.125 0.125; O2 O 0.375 0.375 0.375"]), "gunn2012-pyrochlore", "types", "rules", ["--no-tail"]),
    ("La perovskites, Cherry 1995 (Table 1): cubic LaFeO3, all shells", ("crystal", ["--group", "P m -3 m", "--cell", "3.93,3.93,3.93", "--sites",
     "La1 La 0.5 0.5 0.5; Fe1 Fe 0 0 0; O1 O 0.5 0 0", "--supercell", "3,3,3"]), "la-perovskites", "types", "rules", ["--no-tail"]),
    ("La perovskites, Cherry 1995 (Table 1): cubic LaCoO3", ("crystal", ["--group", "P m -3 m", "--cell", "3.83,3.83,3.83", "--sites",
     "La1 La 0.5 0.5 0.5; Co1 Co 0 0 0; O1 O 0.5 0 0", "--supercell", "3,3,3"]), "la-perovskites", "types", "rules", ["--no-tail"]),
    ("La perovskites, Read 1999 (Table 1): La2NiO4 (I4/mmm), rigid Ni", ("crystal", ["--group", "I 4/m m m", "--cell", "3.86,3.86,12.65", "--sites",
     "La1 La 0 0 0.361; Ni1 Ni 0 0 0; O1 O 0 0.5 0; O2 O 0 0 0.176", "--supercell", "3,3,1"]), "la-perovskites", "types", "rules", ["--no-tail"]),
    # zeolites: Schroder et al. 1992 (framework O shells, hydroxyl O-H Morse pair within 1.4 A with its Coulomb term kept)
    ("Schroder 1992 zeolite set: alpha-quartz, O shells", ("file", os.path.join(ROOT, "data", "crystals", "alpha-quartz.cif")), "schroder1992-zeolite", "types", "rules", ["--no-tail"]),
    ("Schroder 1992 zeolite set: Si(OH)4 and Al(OH)3, hydroxyl Morse pairs", ("smiles", "O[Si](O)(O)O.O[Al](O)O"), "schroder1992-zeolite", "types", "rules"),
    ("Ramsahye 2005 NaY set: alpha-quartz", ("file", os.path.join(ROOT, "data", "crystals", "alpha-quartz.cif")), "ramsahye2005-nay", "types", "rules", ["--no-tail"]),
    ("Ramsahye 2005 NaY set: chloroform, Na+, AlO4 and SiO4 (every cross term)", ("smiles", "ClC(Cl)Cl.ClC(Cl)Cl.[Na+].[O-][Al]([O-])([O-])[O-].[O-][Si]([O-])([O-])[O-]"),
     "ramsahye2005-nay", "types", "rules"),
    ("Forester 1997 silicalite set: alpha-quartz, O shells", ("file", os.path.join(ROOT, "data", "crystals", "alpha-quartz.cif")), "forester1997-silicalite", "types", "rules", ["--no-tail"]),
    ("Forester 1997 silicalite set: benzene with SiO4 (benzene-O shell pairs)", ("smiles", "c1ccccc1.[O-][Si]([O-])([O-])[O-]"), "forester1997-silicalite", "types", "rules"),
    # TraPPE-EH (bench/ff/make_trappe_eh.py): per-model charges, C-H bond-centre sites as the H atoms, 1-4 excluded
    ("TraPPE-EH (TraPPE database): octane, methylamine, acetamide, acetonitrile, nitroethane (flexible)", ("smiles",
     "CCCCCCCC.CN.CC(N)=O.CC#N.CC[N+](=O)[O-]", os.path.join(FF, "trappe-eh-umn.json")), "trappe-eh-umn", "types", "rules"),
    ("TraPPE-EH (TraPPE database): benzene, pyridine, phenol, 2-chlorothiophene, naphthalene", ("smiles",
     "c1ccccc1.c1ccncc1.Oc1ccccc1.Clc1cccs1.c1ccc2ccccc2c1"), "trappe-eh-umn", "types", "rules"),
    # zeolite sorbates (bench/ff/make_zeolite_sorbates.py)
    ("O'Malley 2013 silicalite set: alpha-quartz, formal charges, O-Si-O bends by distance", ("file", os.path.join(ROOT, "data", "crystals", "alpha-quartz.cif")),
     "omalley2013-alkanes", "types", "rules", ["--no-tail"]),
    ("O'Malley 2013 all-atom alkanes: octane, isobutane, neopentane, methane with SiO4 (A/B pairs)", ("smiles",
     "CCCCCCCC.CC(C)C.CC(C)(C)C.C.[O-][Si]([O-])([O-])[O-]"), "omalley2013-alkanes", "types", "rules"),
    ("Plant 2007 methanol in NaY: two methanols, Na+, AlO4 and SiO4 (every cross term)", ("smiles",
     "CO.CO.[Na+].[O-][Al]([O-])([O-])[O-].[O-][Si]([O-])([O-])[O-]"), "plant2007-methanol", "types", "rules"),
    ("Peguin 2009 HFA134a (FF no. 2): two molecules", ("smiles", "FCC(F)(F)F.FCC(F)(F)F"), "peguin2009-hfa134a", "types", "rules"),
    # Hill & Sauer 1995 (CFF class II, every cross term; charges from the paper's bond increments)
    ("Hill & Sauer 1995: disilicic acid and a bridging-hydroxyl Si-O(H)-Al cluster", ("smiles",
     "O[Si](O)(O)O[Si](O)(O)O.O[Si](O)(O)[OH+][Al-](O)(O)O"), "hill-sauer1995", "types", "rules"),
    ("Hill & Sauer 1995: (HO)3Si-O-Al(OH)3 anion and a four-ring", ("smiles", "O[Si](O)(O)O[Al-](O)(O)O.O[Si]1(O)O[Si](O)(O)O[Si](O)(O)O[Si](O)(O)O1"),
     "hill-sauer1995", "types", "rules"),
    # oxide glasses (bench/ff/make_glass.py): the Teter partial-charge set, Tilocca's shell model (SM0; SM1's screened bend has
    # no LAMMPS style and is refused), Vogel's LiPO3; overlays by --typing
    ("Teter glass set (Cormack 2002 / Du 2009): alpha-quartz", ("file", os.path.join(ROOT, "data", "crystals", "alpha-quartz.cif")), "teter-glass", "types", "rules", ["--no-tail"]),
    ("Teter glass set (Du 2009): corundum, Al 1.8", ("file", os.path.join(ROOT, "data", "crystals", "corundum.cif")), "teter-glass", "types", "rules", ["--no-tail"]),
    ("Teter glass set (Du 2009): Y2O3 in the corundum structure", ("crystal", ["--group", "R -3 c", "--cell", "5.6,5.6,15.3", "--sites",
     "Y1 Y 0 0 0.352; O1 O 0.306 0 0.25"]), "teter-glass", "types", "rules", ["--no-tail"]),
    ("Teter glass set, Okhotnikov 2013 Y-O overlay: Y2O3 in the corundum structure", ("crystal", ["--group", "R -3 c", "--cell", "5.6,5.6,15.3", "--sites",
     "Y1 Y 0 0 0.352; O1 O 0.306 0 0.25"]), "teter-glass", "types", "rules",
     ["--no-tail", "--typing", os.path.join(ROOT, "data", "typing", "sets", "teter-glass", "okhotnikov2013-y.typing.json")]),
    ("Teter glass set, Tilocca 2006 Table III overlay: alpha-quartz", ("file", os.path.join(ROOT, "data", "crystals", "alpha-quartz.cif")), "teter-glass", "types", "rules",
     ["--no-tail", "--typing", os.path.join(ROOT, "data", "typing", "sets", "teter-glass", "tilocca2006-teter.typing.json")]),
    ("Tilocca 2006 shell model SM0 (overlay): alpha-quartz", ("file", os.path.join(ROOT, "data", "crystals", "alpha-quartz.cif")), "tilocca2006-sm", "types", "rules",
     ["--no-tail", "--typing", os.path.join(ROOT, "data", "typing", "sets", "tilocca2006-sm", "sm0.typing.json")]),
    ("Vogel 2003 LiPO3: a cubic Li-P-O3 test cell", ("crystal", ["--group", "P m -3 m", "--cell", "3.6,3.6,3.6", "--sites",
     "Li1 Li 0 0 0; P1 P 0.5 0.5 0.5; O1 O 0.5 0.5 0", "--supercell", "3,3,3"]), "vogel2003-lipo3", "types", "rules", ["--no-tail"]),
    ("GaN (Wonglakhon & Zahn 2020): wurtzite, formal charges, LJ + Born-Mayer", ("crystal", ["--group", "P 63 m c", "--cell", "3.189,3.189,5.185,90,90,120",
     "--sites", "Ga1 Ga 0.33333333 0.66666667 0; N1 N 0.33333333 0.66666667 0.377", "--supercell", "4,4,3"]), "gan2020", "types", "rules", ["--no-tail"]),
    # Bush, Gale, Catlow & Battle 1994 (one shell-model set for the oxides): MgO, SrTiO3, corundum
    ("Bush 1994: MgO rock salt, core-shell", ("crystal", ["--group", "F m -3 m", "--cell", "4.212,4.212,4.212", "--sites", "Mg1 Mg 0 0 0; O1 O 0.5 0.5 0.5",
     "--supercell", "3,3,3"]), "bush1994-oxides", "types", "rules", ["--no-tail"]),
    ("Bush 1994: SrTiO3 perovskite, core-shell", ("crystal", ["--group", "P m -3 m", "--cell", "3.905,3.905,3.905", "--sites",
     "Sr1 Sr 0 0 0; Ti1 Ti 0.5 0.5 0.5; O1 O 0.5 0.5 0", "--supercell", "3,3,3"]), "bush1994-oxides", "types", "rules", ["--no-tail"]),
    ("Bush 1994: corundum (alpha-Al2O3), core-shell", ("file", os.path.join(ROOT, "data", "crystals", "corundum.cif")), "bush1994-oxides", "types", "rules", ["--no-tail"]),
    ("Polystyrene melt, GAFF2 (periodic, 1300 atoms)", ("file", os.path.join(ROOT, "samples", "ps_melt.data")), "gaff2-2.2.20-ambertools", "gasteiger", "rules"),
    ("Polystyrene melt, PCFF (periodic, class II + class I)", ("file", os.path.join(ROOT, "samples", "ps_melt.data")), "pcff", "types", "rules"),
    # CVFF, PCFF and COMPASS from their .frc files (msi2lmp's assignment; check_msi2lmp.py compares with msi2lmp itself)
    ("Polystyrene melt, PCFF from pcff.frc (periodic, full class II)", ("file", os.path.join(ROOT, "samples", "ps_melt.data")), "pcff-frc", "types", "rules"),
    ("PCFF from pcff.frc: ester, ether, amide, siloxane, amine", ("smiles", "CCOC(=O)C.COCC.CC(=O)NC.C[Si](C)(C)O[Si](C)(C)C.CCN(C)C"), "pcff-frc", "types", "rules"),
    ("Polystyrene melt, CVFF from cvff.frc (periodic)", ("file", os.path.join(ROOT, "samples", "ps_melt.data")), "cvff-frc", "types", "rules"),
    ("CVFF from cvff.frc: phenol, ethyl acetate, N-methylacetamide", ("smiles", "Oc1ccccc1.CCOC(=O)C.CC(=O)NC"), "cvff-frc", "types", "rules"),
    ("Polystyrene melt, COMPASS from compass_published.frc (periodic)", ("file", os.path.join(ROOT, "samples", "ps_melt.data")), "compass-frc", "types", "rules"),
    ("COMPASS from compass_published.frc: propylbenzene, ester-ether", ("smiles", "CCCc1ccccc1.COCCOC(=O)C"), "compass-frc", "types", "rules"),
    # every bonded kind hybrid: COMPASS (class II with all cross terms) plus an overlay turning one angle, one torsion
    # and one improper type into class I forms, so each class II section carries skip lines. (The improper is cvff,
    # K[1 − cos 2φ]: smooth. A harmonic improper with χ0 = 0 on a class II centre-second quadruple sits at χ ≈ 180°,
    # a cusp of K(χ − χ0)² where the two programs take different one-sided slopes.)
    ("COMPASS polystyrene + class I overlay (hybrid in every kind)", ("compass-ps",), "compass-published-moltemplate", "types", "keys"),
    # OPLS-AA with its numbered types and their own charges (the 2024 and BOSS 2008 files), L-OPLS for long alkyl chains:
    # rubber-relevant groups — trisubstituted alkene (isoprene), styrene, nitrile, ester, chloroalkene, disulfide
    ("OPLS-AA 2024: isoprene, ethylbenzene, nitrile, ester, disulfide", ("smiles", "CC=C(C)CC.CCc1ccccc1.CCC#N.CCOC(=O)C.CSSC"), "oplsaa2024-moltemplate", "types", "rules"),
    ("OPLS-AA 2024 polystyrene melt (periodic)", ("file", os.path.join(ROOT, "samples", "ps_melt.data")), "oplsaa2024-moltemplate", "types", "rules"),
    ("OPLS-AA 2008 (BOSS 4.8): isoprene, chloroethene, NBR nitrile (R2CH-CN); nudged off linear", ("smiles-nudged", "CC=C(C)CC.CC=CCl.CCC(C#N)CC"), "oplsaa2008-moltemplate", "types", "rules"),
    ("L-OPLS 2024: hexadecane", ("smiles", "CCCCCCCCCCCCCCCC"), "loplsaa2024-moltemplate", "types", "rules"),
    ("CGenFF methane template (separate 1-4 LJ)", ("template", "CHARMM36_cgenff", "toluene"), "cgenff", "gasteiger", "rules"),
    # CHARMM in its own styles: lj/charmmfsw (force switch 10-12 Å, pairs in the switching shell), the 1-4 pairs with
    # their own ε14 / σ14 through dihedral charmmfsw weights (a phenyl ring's para pairs reached by two torsions once)
    ("Polystyrene melt, CGenFF (periodic; lj/charmmfsw, 1-4 by dihedral weights)", ("file", os.path.join(ROOT, "samples", "ps_melt.data")), "cgenff", "gasteiger", "rules"),
    # the libraries CAPS builds from the force fields' own files (bench/ff/open: CHARMM's toppar, AmberTools' dat/leap/parm):
    # CHARMM in its own styles (Urey-Bradley angles, multi-term torsions, periodic impropers, separate 1-4 LJ; run with
    # --native), GAFF / GAFF2 typed by antechamber's rules
    ("CHARMM36 toppar: peptide GAFKW (Urey-Bradley, 1-4 LJ)", ("peptide", "GAFKW"), "charmm36-prot-toppar", "gasteiger", "rules"),
    ("CHARMM22 toppar: peptide GAFKW", ("peptide", "GAFKW"), "charmm22-prot-toppar", "gasteiger", "rules"),
    ("CGenFF toppar: toluene, N-methylacetamide, ethyl acetate", ("smiles", "Cc1ccccc1.CC(=O)NC.CCOC(=O)C"), "cgenff-toppar", "gasteiger", "rules"),
    ("CHARMM36 carbohydrates toppar: glucose", ("smiles", "OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O"), "charmm36-carb-toppar", "gasteiger", "rules"),
    ("CHARMM36 nucleic acids toppar: thymidine", ("smiles", "CC1=CN([C@H]2C[C@H](O)[C@@H](CO)O2)C(=O)NC1=O"), "charmm36-nucl-toppar", "gasteiger", "rules"),
    ("CHARMM36 lipids toppar: phosphatidylcholine", ("smiles", "CCCCCCCC(=O)OC[C@H](COP(=O)([O-])OCC[N+](C)(C)C)OC(=O)CCCCCCC"), "charmm36-lipid-toppar", "qeq", "rules"),
    # CHARMM19 (param19 / toph19, united atom: only polar H and formyl H explicit; 1-4 Coulomb x 0.4)
    ("CHARMM19 toppar (united atom): peptide GAFKW", ("peptide", "GAFKW"), "charmm19-toppar", "gasteiger", "rules"),
    ("AMBER ff03: peptide GAFKW", ("peptide", "GAFKW"), "amber-ff03-ambertools", "gasteiger", "rules"),
    ("AMBER ff03: thymidine (N*, CM, CT / H2 sugar)", ("smiles-nudged", "CC1=CN([C@H]2C[C@H](O)[C@@H](CO)O2)C(=O)NC1=O"), "amber-ff03-ambertools", "gasteiger", "rules"),
    # GLYCAM06j (its 1-4 terms unscaled): a methyl glucoside and a beta-1,4 disaccharide
    ("GLYCAM06j: methyl alpha-D-glucopyranoside, cellobiose", ("smiles", "CO[C@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O.OC[C@H]1O[C@@H](O[C@H]2[C@H](O)[C@@H](O)[C@H](O)O[C@@H]2CO)[C@H](O)[C@@H](O)[C@@H]1O"),
     "glycam06j-ambertools", "gasteiger", "rules"),
    ("AMBER ff14SB: peptide GAFKWDR (CX, 2C / 3C, C8, CO)", ("peptide", "GAFKWDR"), "amber-ff14sb-ambertools", "gasteiger", "rules"),
    # OPLS-AA from GROMACS's oplsaa.ff (type charges): a peptide with its own C alpha, Arg, His, Trp types; small organics
    ("OPLS-AA (GROMACS oplsaa.ff): peptide GHWRK", ("peptide", "GHWRK"), "oplsaa-gromacs", "types", "rules"),
    ("OPLS-AA (GROMACS oplsaa.ff): isoprene, ethylbenzene, ester, disulfide, phenol", ("smiles", "CC=C(C)CC.CCc1ccccc1.CCOC(=O)C.CSSC.Oc1ccccc1"),
     "oplsaa-gromacs", "types", "rules"),
    # OPLS-AA 2024 from Jorgensen's own files (type charges)
    ("OPLS-AA 2024 (Jorgensen): isoprene, ethylbenzene, ester, disulfide, PVA unit, perfluorobutane", ("smiles", "CC=C(C)CC.CCc1ccccc1.CCOC(=O)C.CSSC.CC(O)CC(O)C.FC(F)(F)C(F)(F)C(F)(F)C(F)(F)F"),
     "oplsaa2024", "types", "rules"),
    ("OPLS-AA 2024 (Jorgensen): polystyrene melt (periodic)", ("file", os.path.join(ROOT, "samples", "ps_melt.data")), "oplsaa2024", "types", "rules"),
    # DREIDING from its paper (bench/ff/make_dreiding.py): bond constants by perceived order, cosine-harmonic angles,
    # the paper's torsion rules, planar inversions, hydrogen bonds
    ("DREIDING (Mayo 1990): butadienyl amide, benzene ring, phenol, nitrile", ("smiles", "C=CC=CCC(=O)NCc1ccccc1O.CCC#N"), "dreiding1990", "gasteiger", "rules"),
    ("DREIDING (Mayo 1990) ethanol liquid: hydrogen bonds (periodic)", ("solvent-box", "ethanol", 26), "dreiding1990", "gasteiger", "rules"),
    ("OPLS-AA/M (Jorgensen lab): peptide GHWRKD", ("peptide", "GHWRKD"), "oplsaam-jorgensen", "types", "rules"),
    # CL&P ionic liquids from il.ff: an imidazolium / NTf2 pair, pyrrolidinium dicyanamide, ethylsulfate (charges of the types)
    ("CL&P (il.ff): [C4C1im][NTf2], [C4C1pyrr][dca], ethylsulfate, BF4; nudged off linear", ("smiles-nudged", "CCCCn1cc[n+](C)c1.FC(F)(F)S(=O)(=O)[N-]S(=O)(=O)C(F)(F)F.CCCC[N+]1(C)CCCC1.N#C[N-]C#N.CCOS(=O)(=O)[O-].F[B-](F)(F)F"),
     "clandp", "types", "rules"),
    ("GAFF 1.81 (gaff.dat): ester, epoxide, amide, thioether", ("smiles", "CCOC(=O)C.C1CO1.CC(=O)NC.CSC"), "gaff-1.81-ambertools", "gasteiger", "rules"),
    ("GAFF2 2.11 (gaff211.dat): isoprene, styrene, nitrile; nudged off linear", ("smiles-nudged", "CC=C(C)CC.C=Cc1ccccc1.CCC#N"), "gaff2-2.11-ambertools", "gasteiger", "rules"),
    ("GAFF2 2.2.20 (gaff2.dat): polystyrene melt (periodic)", ("file", os.path.join(ROOT, "samples", "ps_melt.data")), "gaff2-2.2.20-ambertools", "gasteiger", "rules"),
    # UFF (every element): Fourier and periodic angles (linear, trigonal, square planar, octahedral, trigonal
    # bipyramid), sp2 and pyramidal-P inversions (improper fourier), group-16 torsions, full 1-4 van der Waals
    ("UFF mixed elements (P, S, Si, Pt, F, Cl)", ("smiles", "CC#CC(=O)Oc1ccc(cc1)P(C)C.F[S](F)(F)(F)(F)F.N[Pt](N)(Cl)Cl."
                                               "FP(F)(F)(F)F.C[Si](C)(C)O[Si](C)(C)C.CSSC"), "uff", "types", "rules"),
    # ClayFF (Cygan 2004) on LAMMPS's own ClayFF test structure (msi2lmp's pyrophyllite, its types and charges): O-H
    # bonds, the M-O-H bends by contact, Lorentz-Berthelot mixing. Ewald: this LAMMPS build's PPPM fails on that cell
    # (its own in.PyAC_bulk-clayff as well).
    ("ClayFF pyrophyllite (LAMMPS's ClayFF test, triclinic, M-O-H bends)", ("file", os.path.expanduser("~/lammps/tools/msi2lmp/test/PyAC_bulk-clayff.car")),
     "inorganic-clay", "keep", "names", ["--kspace", "ewald"]),
    ("ClayFF (Cygan 2004, Tables 1 and 2) pyrophyllite: LAMMPS's ClayFF test structure, its own types", ("file", os.path.expanduser("~/lammps/tools/msi2lmp/test/PyAC_bulk-clayff.car")),
     "clayff2004", "types", "names", ["--kspace", "ewald"]),
    # IFF typed by rule on slabs CAPS builds: hydroxylated quartz (sc4, oc23, oc24, hoy; the charges IFF states from
    # bond increments), a copper (111) slab (Lennard-Jones atoms, perceived bonds dropped)
    ("IFF (PCFF) hydroxylated quartz (001) slab, typed by rule", ("surface", os.path.join(ROOT, "data", "crystals", "alpha-quartz.cif"),
     ["--hkl", "0,0,1", "--layers", "3", "--supercell", "2,2", "--passivate"]), "iff-pcff", "types", "rules"),
    ("IFF (CVFF) hydroxylated quartz (001) slab, typed by rule", ("surface", os.path.join(ROOT, "data", "crystals", "alpha-quartz.cif"),
     ["--hkl", "0,0,1", "--layers", "3", "--supercell", "2,2", "--passivate"]), "iff-cvff", "types", "rules"),
    ("IFF copper (111) slab, typed by rule (no bonds)", ("surface", os.path.join(ROOT, "data", "crystals", "copper.cif"),
     ["--hkl", "1,1,1", "--layers", "4", "--supercell", "3,3"]), "iff-pcff", "types", "rules"),
    ("Polystyrene melt, UFF (periodic, 1300 atoms)", ("file", os.path.join(ROOT, "samples", "ps_melt.data")), "uff", "types", "rules"),
    # TraPPE in its own styles: the c1-c3 torsions as dihedral opls, 14 Å with TraPPE's 1-4 exclusions (united atom;
    # explicit-hydrogen aromatics)
    ("TraPPE-UA: dodecane, 2-methylpentane (opls torsions)", ("smiles", "CCCCCCCCCCCC.CC(C)CCC"), "trappe-ua-umn", "keep", "rules"),
    ("TraPPE-UA (foyer): butanol, 2-propanol, tert-butanol, 2,3-dimethylbutane, neopentane (RB torsions as fourier)",
     ("smiles", "CCCCO.CC(C)O.CC(C)(C)O.CC(C)C(C)C.CC(C)(C)C"), "trappe-ua-foyer", "types", "rules"),
    ("TraPPE-UA (TraPPE database): DME alone (1-4 Coulomb x 0.5 for every molecule)", ("smiles", "COCCOC"), "trappe-ua-umn", "types", "rules"),
    ("TraPPE-UA (TraPPE database): DME with dimethyl ether (1-4 Coulomb x 0.5 for DME only: a sub-style per part)", ("smiles", "COCCOC.COC"), "trappe-ua-umn", "types", "rules"),
    ("TraPPE (TraPPE database): CO2, ethylene oxide, decane (the models with charge sites go to GROMACS); nudged off linear", ("smiles-nudged", "O=C=O.C1CO1.CCCCCCCCCC"),
     "trappe-ua-umn", "types", "rules"),
    ("TraPPE-UA (TraPPE database): hexane, DME (1-4 Coulomb x 0.5), methyl methacrylate, ethanol, propionitrile, cis-2-butene", ("smiles", "CCCCCC.COCCOC.CC(=C)C(=O)OC.CCO.CCC#N.C/C=C\\C"),
     "trappe-ua-umn", "types", "rules"),
    ("TraPPE-UA alkanes (moltemplate): dodecane", ("smiles", "CCCCCCCCCCCC"), "trappe1998-moltemplate", "keep", "rules"),
    # coarse-grained: MARTINI (lj/gromacs + coul/gromacs, dielectric 15, cosine/squared angles, 1-3 and 1-4 pairs kept)
    # and SDK (lj/sdk 9-6 / 12-4, angle sdk with its 1-3 repulsion), bead molecules and water beads in periodic boxes
    ("MARTINI DPPC + POPE + ions + water (periodic)", ("cg-box", "martini-moltemplate", [("template", "DPPC", 6), ("template", "POPE", 4),
                                                   ("template", "NA+", 5), ("template", "CL-", 5), ("template", "W", 150)], 42.0),
     "martini-moltemplate", "keep", "rules"),
    ("MARTINI polymers: PEO (torsions) + water (periodic)", ("cg-box", "martini-polymers", [("template", "PEO", 3), ("beads", "[P4]", 120)], 44.0),
     "martini-polymers", "keep", "rules"),
    ("MARTINI sugars: sucrose, maltose, glucose + water (periodic)", ("cg-box", "martini-sugars", [("template", "SUCR", 4), ("template", "MALT", 4),
                                                                   ("template", "GLUC", 4), ("beads", "[P4]", 120)], 40.0),
     "martini-sugars", "keep", "rules"),
    ("MARTINI amino acids: HIS, PHE, TYR, TRP (impropers), ARG, ASP, LYS + water (periodic)",
     ("cg-box", "martini-aminoacids", [("template", t, 3) for t in ("HIS", "PHE", "TYR", "TRP", "ARG", "ASP", "LYS")] + [("beads", "[P4]", 100)], 40.0),
     "martini-aminoacids", "keep", "rules"),
    ("Martini 2.2 protein: AK helix (DSSP; helix constraints, BBB angles, backbone dihedrals)", ("peptide", "AEAAAKEAAAKEAAAKA", "--helix"),
     "martini22-proteins", "keep", "rules"),
    ("Martini 2.2 protein: aromatic helix (W, Y, F, H ring constraints and impropers, charged termini)", ("peptide", "AEAWAKEAYAKEAFAKHA", "--helix"),
     "martini22-proteins", "keep", "rules"),
    ("Martini 2.2 protein: 1ICO beta-hairpin (VAL, ILE as AC1 / AC2; elastic bonds, disulfide)", ("file", os.path.join(ROOT, "tests", "data", "vermouth", "1ico_aa.pdb")),
     "martini22-proteins", "keep", "rules"),
    ("SDK DMPC + DMPE + water (periodic)", ("cg-box", "sdk-moltemplate", [("beads", "[NC][PH][GL]([EST1][CM][CM][CM][CT2])[EST2][CM][CM][CM][CT2]", 6),
                                                                     ("beads", "[NH][PHE][GL]([EST1][CM][CM][CM][CT2])[EST2][CM][CM][CM][CT2]", 4),
                                                                     ("beads", "[W]", 120)], 36.0),
     "sdk-moltemplate", "types", "rules"),
    ("SDK C12E8 + water (periodic)", ("cg-box", "sdk-moltemplate", [("beads", "[OA]" + "[EO]" * 8 + "[CM][CM][CM][CT2]", 6), ("beads", "[W]", 120)], 40.0),
     "sdk-moltemplate", "types", "rules"),
    ("SDK SDS + Na+ + water (periodic, relative permittivity 80)", ("cg-box", "sdk-moltemplate", [("beads", "[SO4][CM][CM][CM][CT]", 8),
                                                                   ("beads", "[SOD]", 8), ("beads", "[W]", 120)], 38.0),
     "sdk-moltemplate", "types", "rules"),
    ("Cooke-Deserno lipids (cosine/squared, FENE; periodic)", ("cg-box", "cooke-deserno-moltemplate", [("template", "lipid", 60)], 16.0),
     "cooke-deserno-moltemplate", "types", "rules"),
    # DREIDING's hydrogen bond (hbond/dreiding/lj, moltemplate's DREIDING): water typed O_3_hd / H_HB, every O-H···O within
    # 6.5 Å and past 90°, switched from 6 Å
    ("DREIDING (moltemplate) water: hydrogen bonds (periodic, 480 molecules)", ("water-box", 480, 24.84), "dreiding-moltemplate", "gasteiger", "rules"),
    # mW water: all-atom water packed by CAPS, one Stillinger–Weber site per molecule (pair_style sw with a .sw file)
    ("mW water (Stillinger-Weber, periodic, 480 sites)", ("water-box", 480, 24.84), "mw-moltemplate", "types", "rules"),
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
    if kind == "crystal":   # CAPS's space-group builder: caps crystal ARGS
        out = os.path.join(work, base + ".data")
        subprocess.run([CAPS, "crystal"] + list(src[1]) + ["-o", out], capture_output=True, check=True)
        return out, None
    if kind == "surface":   # CAPS's slab builder: caps surface CIF ARGS
        out = os.path.join(work, base + ".data")
        subprocess.run([CAPS, "surface", src[1], "-o", out] + list(src[2]), capture_output=True, check=True)
        return out, None
    if kind == "solvent-box":   # a liquid from CAPS's solvent packing (caps solvate), bonds from the packed molecules
        out = os.path.join(work, base + ".pdb")
        subprocess.run([CAPS, "solvate", "-o", out, "--edge", str(src[2]), "--solvent", src[1], "--no-ions"], capture_output=True, check=True)
        return out, None
    if kind == "water-box":   # N waters packed 2 Å apart (and 2 Å from the cell faces), in a periodic cube of edge L
        n, edge = src[1], src[2]
        w, box = os.path.join(work, base + ".w.xyz"), os.path.join(work, base + ".xyz")
        with open(w, "w") as f:
            f.write("3\nwater\nO 0 0 0\nH 0.757 0 0.586\nH -0.757 0 0.586\n")
        inp = os.path.join(work, base + ".inp")
        with open(inp, "w") as f:
            f.write(f"tolerance 2.0\nfiletype xyz\noutput {box}\nstructure {w}\n  number {n}\n  inside box 1. 1. 1. {edge - 1} {edge - 1} {edge - 1}\n"
                    "end structure\n")
        subprocess.run([CAPS, "pack", inp, "--quiet"], capture_output=True, check=True)
        lines = open(box).read().splitlines()
        lines[1] = f'Lattice="{edge} 0 0 0 {edge} 0 0 0 {edge}" Properties=species:S:1:pos:R:3'
        with open(box, "w") as f:
            f.write("\n".join(lines) + "\n")
        return box, None
    if kind == "cg-box":   # copies of bead molecules (CAPS's builder), random orientations and places, periodic
        import random
        rng = random.Random(7)
        ffj, parts, edge = ff_file(src[1]), src[2], src[3]
        atoms, bonds, labels = [], [], {}
        for how, text, count in parts:
            tmp = os.path.join(work, base + ".one.data")
            subprocess.run([CAPS, "build", "--" + how, text, "--ff", ffj, "-o", tmp], capture_output=True, check=True)
            lines = open(tmp).read().splitlines()
            sec, mass, mol_at, mol_b = None, {}, [], []
            for l in lines:
                w = l.split()
                if l.strip() in ("Masses", "Atoms  # full", "Bonds"):
                    sec = l.split()[0]
                    continue
                if not w or len(w) < 2:
                    continue
                if sec == "Masses":
                    mass[int(w[0])] = (float(w[1]), w[3])
                elif sec == "Atoms" and len(w) >= 7:
                    mol_at.append((mass[int(w[2])], float(w[3]), float(w[4]), float(w[5]), float(w[6])))
                elif sec == "Bonds" and len(w) >= 4:
                    mol_b.append((int(w[2]) - 1, int(w[3]) - 1))
            cx = [sum(a[k] for a in mol_at) / len(mol_at) for k in (2, 3, 4)]
            for c in range(count):
                for _ in range(2000):
                    # a random rotation (QR of a Gaussian matrix) and place, kept 3.5 A from the others
                    import math
                    q = [rng.gauss(0, 1) for _ in range(4)]
                    nq = math.sqrt(sum(x * x for x in q))
                    a_, b_, c_, d_ = (x / nq for x in q)
                    R = [[a_*a_+b_*b_-c_*c_-d_*d_, 2*(b_*c_-a_*d_), 2*(b_*d_+a_*c_)], [2*(b_*c_+a_*d_), a_*a_-b_*b_+c_*c_-d_*d_, 2*(c_*d_-a_*b_)],
                         [2*(b_*d_-a_*c_), 2*(c_*d_+a_*b_), a_*a_-b_*b_-c_*c_+d_*d_]]
                    t = [rng.uniform(0, edge) for _ in range(3)]
                    new = []
                    for (m, lab), qq, x, y, z in mol_at:
                        d = (x - cx[0], y - cx[1], z - cx[2])
                        new.append(((m, lab), qq, *[t[k] + sum(R[k][j] * d[j] for j in range(3)) for k in range(3)]))
                    def close(p1, p2):
                        return math.sqrt(sum(((p1[k] - p2[k] + edge / 2) % edge - edge / 2) ** 2 for k in range(3))) < (1.2 if "cooke" in src[1] else 3.5)
                    if all(not close(n1[2:], a2[2:]) for n1 in new for a2 in atoms) or _ == 1999:
                        break
                off = len(atoms)
                atoms.extend(new)
                bonds.extend((i + off, j + off) for i, j in mol_b)
        types = sorted({a[0] for a in atoms}, key=lambda t: t[1])
        tid = {t: k + 1 for k, t in enumerate(types)}
        out = os.path.join(work, base + ".data")
        with open(out, "w") as f:
            f.write(f"CAPS · {base}\n\n{len(atoms)} atoms\n{len(bonds)} bonds\n\n{len(types)} atom types\n1 bond types\n\n"
                    f"0 {edge} xlo xhi\n0 {edge} ylo yhi\n0 {edge} zlo zhi\n\nMasses\n\n")
            for t in types:
                f.write(f"{tid[t]} {t[0]}  # {t[1]}\n")
            f.write("\nAtoms  # full\n\n")
            mols = 0
            for k, a in enumerate(atoms):
                f.write(f"{k + 1} 1 {tid[a[0]]} {a[1]} {a[2] % edge:.6f} {a[3] % edge:.6f} {a[4] % edge:.6f}\n")
            f.write("\nBonds\n\n")
            for k, (i, j) in enumerate(bonds):
                f.write(f"{k + 1} 1 {i + 1} {j + 1}\n")
        return out, None
    if kind == "peptide":   # an all-atom peptide from CAPS's builder
        out = os.path.join(work, base + ".pdb")
        subprocess.run([CAPS, "peptide", src[1], "-o", out] + list(src[2:]), capture_output=True, check=True)
        return out, None
    if kind in ("smiles", "smiles-nudged"):   # built and cleaned up by CAPS with UFF (or the force field a third item names)
        m = os.path.join(work, base + ".mol2")
        subprocess.run([CAPS, "build", src[1], "--ff", src[2] if len(src) > 2 else "uff", "-o", m], capture_output=True, check=True)
        if kind == "smiles-nudged":
            # every atom moved by up to 0.03 Å (fixed seed): a planar centre is no longer within 0.06° of flat, where
            # LAMMPS's improper harmonic clamps 1/sin χ (SMALL = 0.001) and its forces there are not the analytic ones
            import random
            rng = random.Random(7)
            lines, atoms = open(m).read().split("\n"), False
            for i, l in enumerate(lines):
                if l.startswith("@<TRIPOS>"):
                    atoms = l.strip() == "@<TRIPOS>ATOM"
                    continue
                w = l.split()
                if atoms and len(w) >= 6:
                    for c in (2, 3, 4):
                        w[c] = "%.5f" % (float(w[c]) + rng.uniform(-0.03, 0.03))
                    lines[i] = " ".join(w)
            open(m, "w").write("\n".join(lines))
        return m, None
    if kind == "compass-ps":   # polystyrene melt with COMPASS types (from the GAFF2 typing: c3 → c4, ca → c3a, H → h1)
        ps = os.path.join(ROOT, "samples", "ps_melt.data")
        gt = os.path.join(work, base + ".gaff")
        subprocess.run([CAPS, "ff", "type", ps, "--ff", ff_file("gaff2-2.2.20-ambertools"), "-o", gt], capture_output=True, check=True)
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



def harmonic_improper_clamped(d):
    """Atoms of improper-harmonic quadruplets inside LAMMPS's clamp window (see improper_harmonic.cpp)."""
    inp = open(os.path.join(d, "case.in")).read()
    if not re.search(r"improper_style\s+harmonic", inp):
        return set()
    lines = open(os.path.join(d, "case.data")).read().split("\n")
    def section(name):
        i = next((k for k, l in enumerate(lines) if l.strip().split("#")[0].strip() == name), None)
        out = []
        if i is None:
            return out
        i += 2
        while i < len(lines) and lines[i].strip():
            out.append(lines[i].split())
            i += 1
        return out
    pos = {int(w[0]): tuple(float(x) for x in w[4:7]) for w in section("Atoms")}
    sub = lambda a, b: tuple(a[k] - b[k] for k in range(3))
    dot = lambda a, b: sum(a[k] * b[k] for k in range(3))
    out = set()
    for w in section("Impropers"):
        i1, i2, i3, i4 = (int(x) for x in w[2:6])
        v1, v2, v3 = sub(pos[i1], pos[i2]), sub(pos[i3], pos[i2]), sub(pos[i4], pos[i3])
        r1, r2, r3 = (math.sqrt(dot(v, v)) for v in (v1, v2, v3))
        c0 = dot(v1, v3) / (r1 * r3)
        c1 = dot(v1, v2) / (r1 * r2)
        c2 = -dot(v3, v2) / (r3 * r2)
        s1, s2 = 1 - c1 * c1, 1 - c2 * c2
        c = max(-1.0, min(1.0, (c1 * c2 + c0) / math.sqrt(max(s1, 0.001) * max(s2, 0.001))))
        if s1 < 0.001 or s2 < 0.001 or math.sqrt(1 - c * c) < 0.001:
            out |= {i1, i2, i3, i4}
    return out


def lammps_hbond_offset(d):
    """LAMMPS's hbond/dreiding/lj force in its switching region (r_in < r < r_out) takes the switch's derivative without
    cos^n(theta) (pair_hbond_dreiding_lj.cpp: force_switch = eng_lj*switch2/rsq), so its forces are not the derivative of
    its energy there; CAPS's are. This gives, per atom, LAMMPS's force minus the exact one, and the virial of it (kcal/mol),
    from the files CAPS wrote: -(1 - cos^n) U(r) S'(r) along D-A on the donor, the opposite on the acceptor."""
    inp = open(os.path.join(d, "case.in")).read()
    m = re.search(r"hbond/dreiding/lj (\d+) (\S+) (\S+) (\S+)", inp)
    if not m:
        return {}, [0.0] * 6
    rin, rout, cut = float(m.group(2)), float(m.group(3)), math.cos(math.radians(float(m.group(4))))
    terms = {}
    for w in re.findall(r"^pair_coeff\s+(\d+) (\d+) hbond/dreiding/lj (\d+) ([ij]) (\S+) (\S+) (\d+)", inp, re.M):
        i, j, h, flag, eps, sig, n = int(w[0]), int(w[1]), int(w[2]), w[3], float(w[4]), float(w[5]), int(w[6])
        dt, at = (i, j) if flag == "i" else (j, i)
        terms[(dt, at)] = (h, eps, sig, n)
    txt = open(os.path.join(d, "case.data")).read()
    box = []
    for ax in "xyz":
        mm = re.search(r"(\S+)\s+(\S+)\s+" + ax + "lo " + ax + "hi", txt)
        box.append(float(mm.group(2)) - float(mm.group(1)))
    at = {}
    for l in txt.split("Atoms", 1)[1].split("\n\n", 2)[1].splitlines():
        w = l.split()
        if len(w) >= 7:
            at[int(w[0])] = (int(w[2]), tuple(map(float, w[4:7])))
    nb = {}
    for l in txt.split("\nBonds", 1)[1].split("\n\n", 2)[1].splitlines() if "\nBonds" in txt else []:
        w = l.split()
        if len(w) >= 4:
            a, b = int(w[2]), int(w[3])
            nb.setdefault(a, []).append(b)
            nb.setdefault(b, []).append(a)
    mi = lambda u, v: [(u[k] - v[k]) - box[k] * round((u[k] - v[k]) / box[k]) for k in range(3)]
    donors = {t for t, _ in terms}
    accs = {a for _, a in terms}
    off, wv = {}, [0.0] * 6
    ro2, ri2 = rout * rout, rin * rin
    den = (ro2 - ri2) ** 3
    ids = sorted(at)
    for D in ids:
        tD = at[D][0]
        if tD not in donors:
            continue
        for A in ids:
            if A == D or (tD, at[A][0]) not in terms:
                continue
            h, eps, sig, n = terms[(tD, at[A][0])]
            dA = mi(at[A][1], at[D][1])   # A - D
            r2 = sum(x * x for x in dA)
            if not (ri2 < r2 < ro2):
                continue
            r = math.sqrt(r2)
            U = eps * (5 * (sig / r) ** 12 - 6 * (sig / r) ** 10)
            dS = -12 * r * (ro2 - r2) * (r2 - ri2) / den
            for H in nb.get(D, []):
                if at[H][0] != h:
                    continue
                d1 = mi(at[D][1], at[H][1])
                d2 = [dA[k] + d1[k] for k in range(3)]
                c = sum(d1[k] * d2[k] for k in range(3)) / math.sqrt(sum(x * x for x in d1) * sum(x * x for x in d2))
                if not c < cut:
                    continue
                g = -(1 - c ** n) * U * dS / r   # along D - A = -dA, on the donor
                fD = [-g * dA[k] for k in range(3)]
                for k in range(3):
                    off.setdefault(D, [0.0] * 3)[k] += fD[k]
                    off.setdefault(A, [0.0] * 3)[k] -= fD[k]
                dm = [-x for x in dA]   # D - A
                wv[0] += dm[0] * fD[0]; wv[1] += dm[1] * fD[1]; wv[2] += dm[2] * fD[2]
                wv[3] += dm[0] * fD[1]; wv[4] += dm[0] * fD[2]; wv[5] += dm[1] * fD[2]
    return off, wv

DSF = re.compile(r"(?:^|\s)(?:lj/cut/)?coul/dsf (\S+) (\S+)")
LONG = re.compile(r"(?:^|\s)(?:lj/cut/|lj/class2/)?coul/long (\S+)")

def dsf_self_repeat(d):
    """LAMMPS's coul/dsf adds each atom's self energy -(e_shift/2 + a/sqrt(pi)) q^2 in every coul/dsf sub-style
    whose pairs it takes part in: with a sub-style per part (force fields or molecules with their own 1-4 scaling),
    an atom of a type in two of them is counted twice. A constant (no force); this gives LAMMPS's surplus (kcal/mol)."""
    inp = open(os.path.join(d, "case.in")).read()
    ps = next((l for l in inp.splitlines() if l.startswith("pair_style")), "")
    m = re.findall(r"coul/dsf\s+(\S+)\s+(\S+)", ps)
    if len(m) < 2:
        return 0.0
    a, rc = float(m[0][0]), float(m[0][1])
    subs = collections.defaultdict(set)   # type -> the coul/dsf sub-styles it takes part in
    for l in inp.splitlines():
        w = l.split()
        if len(w) >= 5 and w[0] == "pair_coeff" and w[3] == "coul/dsf":
            subs[int(w[1])].add(w[4])
            subs[int(w[2])].add(w[4])
    q, ty, sec = {}, {}, None
    for l in open(os.path.join(d, "case.data")):
        w = l.split("#")[0].split()
        if l.startswith("Atoms"):
            sec = "a"
        elif l[:1].isalpha():
            sec = None
        elif sec == "a" and len(w) >= 7:
            q[int(w[0])], ty[int(w[0])] = float(w[3]), int(w[2])
    # LAMMPS's own shifts (pair_coul_dsf.cpp init_style): f_shift = -(erfc(a rc)/rc^2 + 2a/sqrt(pi) exp(-a^2 rc^2)/rc),
    # e_shift = erfc(a rc)/rc - f_shift rc
    f_shift = -(math.erfc(a * rc) / rc ** 2 + 2 * a / math.sqrt(math.pi) * math.exp(-a * a * rc * rc) / rc)
    e_shift = math.erfc(a * rc) / rc - f_shift * rc
    e1 = -(e_shift / 2 + a / math.sqrt(math.pi)) * 332.06371
    return sum(e1 * q[i] ** 2 * (len(subs[ty[i]]) - 1) for i in q if len(subs[ty[i]]) > 1)


def erfc_offset(d, pairs, g_ewald):
    """LAMMPS's coul/dsf and coul/long styles (the latter with pair_modify table 0) take erfc from Abramowitz and Stegun
    7.1.26 (|error| <= 1.5e-7), CAPS the exact erfc. This gives LAMMPS's Coulomb energy minus the exact one (kcal/mol,
    or eV in metal units) over the pairs of its own neighbour list (pairs: atom ids and distance, from the run):
    q_i q_j [erfc_AS(ar) - erfc(ar)] / r with a the DSF damping or LAMMPS's Ewald g, the same for a pair at full weight
    and an excluded one (LAMMPS takes the whole erfc term there and removes the bare 1/r)."""
    inp = open(os.path.join(d, "case.in")).read()
    ps = next((l for l in inp.splitlines() if l.startswith("pair_style")), "")
    m, ml = DSF.search(ps), LONG.search(ps)
    if not pairs or not (m or (ml and g_ewald)):
        return 0.0
    a, rc = (float(m.group(1)), float(m.group(2))) if m else (g_ewald, float(ml.group(1)))
    diel = re.search(r"^dielectric\s+(\S+)", inp, re.M)
    qqrd2e = (14.399645 if re.search(r"^units\s+metal", inp, re.M) else 332.06371) / (float(diel.group(1)) if diel else 1.0)
    q, sec = {}, None
    for l in open(os.path.join(d, "case.data")):
        w = l.split("#")[0].split()
        if l.startswith("Atoms"):
            sec = "a"
        elif l[:1].isalpha():
            sec = None
        elif sec == "a" and len(w) >= 7:
            q[int(w[0])] = float(w[3])
    off = 0.0
    for i, j, r in pairs:
        if r < rc and q[i] and q[j]:
            t = 1 / (1 + 0.3275911 * a * r)
            erfc_as = t * (0.254829592 + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429)))) * math.exp(-a * a * r * r)
            if r < 1e-6:   # a shell on its core (CORESHELL's /cs styles take r = 1e-10 there: no erfc difference to speak of)
                continue
            off += q[i] * q[j] * (erfc_as - math.erfc(a * r)) / r
    return qqrd2e * off

def lammps(infile, dump):
    txt = open(infile).read()
    # coul/dsf and coul/long: every pair of LAMMPS's list with its distance, for its erfc (erfc_offset); coul/long
    # without LAMMPS's 12-bit Coulomb table, whose interpolation is ~1e-7 off. The core-shell styles (born/coul/dsf/cs,
    # coul/long/cs) and TIP4P keep theirs, as their core-shell pairs sit at r = 1e-10 A
    pdump = dump + ".pairs"
    ps = next((l for l in txt.splitlines() if l.startswith("pair_style")), "")
    dsf = (DSF.search(ps) or LONG.search(ps)) is not None
    if dsf:
        if LONG.search(ps):
            txt = txt.replace("run 0", "pair_modify table 0\nrun 0", 1)
        txt = txt.replace("run 0", f"compute caps_pp all property/local patom1 patom2\ncompute caps_pd all pair/local dist\n"
                          f"dump caps_pl all local 1 {pdump} c_caps_pp[1] c_caps_pp[2] c_caps_pd\ndump_modify caps_pl format float %.12g\nrun 0", 1)
    txt = txt.replace("run 0", f"run 0\nwrite_dump all custom {dump} id fx fy fz modify sort id format float %.10f")
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
    pairs = []
    if dsf:
        for l in open(pdump).read().split("ITEM: ENTRIES")[1].splitlines()[1:]:
            w = l.split()
            if len(w) == 3:
                pairs.append((int(w[0]), int(w[1]), float(w[2])))
    g = re.search(r"G vector \(1/distance\) = (\S+)", r.stdout)
    return e, f, pairs, float(g.group(1)) if g else 0.0


rows, fails = [], 0
for case in CASES:
    label, src, fid, charges, typing = case[:5]
    case_extra = case[5] if len(case) > 5 else []   # options of this case only (e.g. Ewald where LAMMPS's PPPM fails)
    if only and only.lower() not in label.lower():
        continue
    base = re.sub(r"[^a-z0-9]+", "_", label.lower()).strip("_")[:40]
    ffj = ff_file(fid)
    s, tfile = structure(src, base) if ffj else (None, None)
    if not ffj or not s:
        rows.append((label, "skipped (source not available)", "", ""))
        continue
    d = os.path.join(work, base)
    os.makedirs(d, exist_ok=True)
    cmd = [CAPS, "ff", "apply", s, "--ff", ffj, "--charges", charges, "-o", os.path.join(d, "case.data"),
           "--lammps-input", os.path.join(d, "case.in"), "--forces", os.path.join(d, "caps_f.txt")]
    cmd += ["--lammps-style", "native"] if NATIVE else ["--lammps-style", "exact"]   # native: the force field's own cut-off, as CAPS evaluates
    if HYBRID:
        cmd += ["--hybrid"]
    if typing == "names":   # the file's own force-field types (a Materials Studio .car)
        cmd += ["--names"]
    cmd += case_extra
    if typing == "keys" and tfile:
        cmd += ["--types", tfile]
    if PME:
        if src[0] not in ("file", "ionic-first", "cg-box", "crystal"):
            continue
        cmd += ["--pme", "--ewald-rtol", "1e-9", "--pme-spacing", "0.3", "--pme-order", "8"]
    if src[0] == "compass-ps":
        ov = os.path.join(work, base + ".overlay.json")
        with open(ov, "w") as f:
            f.write(OVERLAY)
        cmd += ["--overlay", ov]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode or "missing parameters" in r.stdout:
        why = (r.stderr.strip().splitlines() or r.stdout.strip().splitlines() or ["?"])[-1]
        # CHARMM's force switch and 1-4 terms exist in CHARMM's own styles only (lj/charmmfsw): refused in the exact mode
        expected = "1-4 Lennard-Jones" in why or ("CHARMM's own LAMMPS styles only" in why and not NATIVE)
        rows.append((label, ("refused, as intended: " if expected else "CAPS failed: ") + why[:150], "", ""))
        fails += 0 if expected else 1
        continue
    m = re.search(r"energy \(kcal/mol\): bond (\S+)  angle (\S+)  dihedral (\S+)  improper (\S+)  vdW (\S+)  Coulomb (\S+)", r.stdout)
    ce = dict(zip(["bond", "angle", "dihedral", "improper", "vdw", "coulomb"], map(float, m.groups())))
    cf = {int(w[0]): tuple(map(float, w[1:4])) for w in (l.split() for l in open(os.path.join(d, "caps_f.txt"))) if len(w) == 4}
    try:
        le, lf, lpairs, lg = lammps(os.path.join(d, "case.in"), os.path.join(d, "f.dump"))
    except Exception as ex:
        rows.append((label, "LAMMPS failed: " + str(ex), "", ""))
        fails += 1
        continue
    lm = {"bond": le["E_bond"], "angle": le["E_angle"], "dihedral": le["E_dihed"], "improper": le["E_impro"], "vdw": le["E_vdwl"],
          "coulomb": le["E_coul"] + le["E_long"]}
    # torsion constants the opls style cannot hold (TraPPE's c0) are left out of the LAMMPS file; CAPS reports their total
    c0 = re.search(r"the dihedral energy lower by (\S+) kcal/mol", r.stdout)
    if c0:
        lm["dihedral"] += float(c0.group(1))
    # LAMMPS's coul/dsf and coul/long real-space energy with the exact erfc, as CAPS's
    lm["coulomb"] -= erfc_offset(d, lpairs, lg)
    lm["coulomb"] -= dsf_self_repeat(d)
    # with PME both codes reach the Ewald limit only to their discretisation (CAPS's grid, LAMMPS's k-space cut-off; the
    # real-space erfc is exact on both sides after erfc_offset): the Coulomb term is compared to 1e-3 kcal/mol absolute there
    de = max(abs(ce[k] - lm[k]) / (max(1.0, abs(ce[k])) if not (PME and k == "coulomb") else 100.0) for k in ce)
    mw = re.search(r"virial tensor \(kcal/mol\): xx (\S+)  yy (\S+)  zz (\S+)  xy (\S+)  xz (\S+)  yz (\S+)", r.stdout)
    cw = list(map(float, mw.groups()))
    # DREIDING hydrogen bonds: LAMMPS's own (inexact) force in the switching region, rebuilt, added to CAPS's exact one
    hoff, hw = lammps_hbond_offset(d)
    for a, v in hoff.items():
        cf[a] = tuple(cf[a][k] + v[k] for k in range(3))
    cw = [cw[k] + hw[k] for k in range(6)]
    lw = [le[f"c_pv[{k}]"] * le["Volume"] / 68568.415 for k in range(1, 7)]
    # relative to the size of the terms it sums: an ionic crystal's Coulomb and repulsive virials (each of the order of
    # their energies) nearly cancel, so the net virial alone would magnify LAMMPS's 1e-7 erfc in Coulomb (formal charges)
    scale = max(1.0, max(abs(w) for w in cw), abs(ce.get("coulomb", 0.0)), abs(ce.get("vdw", 0.0)))
    dw = max(abs(a - b) for a, b in zip(cw, lw)) / scale
    # improper harmonic: LAMMPS clamps 1 − cos² of both bond angles and sin χ at 0.001 (improper_harmonic.cpp SMALL), so
    # at a (near-)flat centre its forces are not the derivative of its energy: those atoms are left out of the force
    # comparison (their energy is still compared) and counted
    clamped = harmonic_improper_clamped(d)
    df = max((math.dist(cf[i], lf[i]) for i in cf if i not in clamped), default=0.0)
    styles = " · ".join(l.strip() for l in open(os.path.join(d, "case.in")) if re.match(r"(bond|angle|dihedral|improper|pair)_style", l))
    if NATIVE:   # Coulomb by another method: the other terms only
        de = max(abs(ce[k] - lm[k]) / max(1.0, abs(ce[k])) for k in ce if k != "coulomb")
        ok = de < 1e-5
        rows.append((label, f"{'ok' if ok else 'DIFFERS'} · bonded and van der Waals terms {de:.1e} (relative); Coulomb CAPS DSF {ce['coulomb']:.3f} / LAMMPS "
                     f"{lm['coulomb']:.3f} kcal/mol", " · ".join(l.strip() for l in open(os.path.join(d, "case.in")) if re.match(r"(bond|angle|dihedral|improper|pair|kspace)_style", l)),
                     " ".join(f"{k} {ce[k]:.4f}/{lm[k]:.4f}" for k in ce if abs(ce[k]) > 0 or abs(lm[k]) > 0)))
        fails += 0 if ok else 1
        continue
    ok = de < 1e-5 and df < 1e-3 and dw < 1e-5
    fails += 0 if ok else 1
    rows.append((label, f"{'ok' if ok else 'DIFFERS'} · energy terms {de:.1e} (relative) · forces {df:.1e} kcal/mol/Å · virial tensor {dw:.1e} (relative)" +
                 (f" · {len(clamped)} atoms of flat harmonic impropers left out of the forces (LAMMPS clamps there)" if clamped else ""), styles,
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
