"""Builds the polymer crystals of data/crystals/polymers from their original crystal-structure papers.

Each entry below holds the paper's space group (its equivalent positions), cell and atomic coordinates as printed, with
the page they are on. Papers that give only the heavy atoms get their hydrogens here, placed geometrically:

  CH2 (sp3, two heavy neighbours)  C-H 1.09 A, H-C-H 109.47 deg, in the plane that bisects the X-C-Y angle
  CH  (sp3, three heavy neighbours) C-H 1.09 A, opposite the sum of the three bond directions
  CH3 (sp3, one heavy neighbour A)  C-H 1.09 A, A-C-H 109.47 deg, staggered (H anti to A's nearest other heavy neighbour)
  sp2 CH (aromatic, two neighbours) C-H 1.08 A, in the plane, opposite the bisector
  NH  (amide N, two neighbours)     N-H 1.01 A, in the C-N-C plane, opposite the bisector

A paper that gives only a lattice and the helix symmetry (PTFE) gets its chain generated from the helix; see ptfe_sites.

The cell is the paper's (the DFT start: VASP relaxes it), the heavy atoms are the paper's (corrections of misprints are
listed per entry, with the reason), and the cell is written in P1 with every atom, so any CIF reader gets the same
structure. polymers.json holds what the tests check against: the paper's cell, space group, Z and density.

usage: python3 tools/make_polymer_crystals.py [--write]   (without --write: build and check only)
"""
import json
import math
import os
import subprocess
import sys
import tempfile
from fractions import Fraction

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "data", "crystals", "polymers")
CAPS = os.environ.get("CAPS_EXE", os.path.join(ROOT, "build", "cli", "caps"))
PRERELAX_FF = "oplsaa2024.json"
MASS = {"H": 1.008, "C": 12.011, "N": 14.007, "O": 15.999, "F": 18.998}
NA = 6.02214076e23

# ---------------------------------------------------------------------------------------------------------- the papers
CRYSTALS = [
    {
        "id": "pe-orthorhombic", "name": "Polyethylene, orthorhombic", "polymer": "polyethylene (PE)", "repeat": "C2H4",
        "source": "C. W. Bunn, Trans. Faraday Soc. 35, 482 (1939)", "doi": "10.1039/tf9393500482",
        "space_group": "Pnam", "number": 62,
        "cell": (7.40, 4.93, 2.534, 90, 90, 90), "chain_axis": "c", "Z": 2, "chains": 2,
        "density": None,
        "density_note": "the paper gives no density; 4 CH2 in the paper's cell give 1.008 g/cm3",
        "ops": ["x,y,z", "-x,-y,-z", "x,y,1/2-z", "-x,-y,1/2+z", "1/2+x,1/2-y,z", "1/2-x,1/2+y,-z", "1/2-x,1/2+y,1/2+z", "1/2+x,1/2-y,1/2-z"],
        "sites": [("C1", "C", 0.038, -0.065, 0.25)],
        "where": "cell and space group p. 483; carbon positions p. 490 (x = 0.038 a, y = 0.065 b; 0.038, .935, .250 ...)",
        "hydrogens": {"C1": "CH2"},
        "use": "PE crystal: VASP and phono3py thermal conductivity; the crystalline layers of semicrystalline PE", "hkl": (1, 1, 0),
    },
    {
        "id": "pe-monoclinic", "name": "Polyethylene, monoclinic", "polymer": "polyethylene (PE)", "repeat": "C2H4",
        "source": "K. Seto, T. Hara & K. Tanaka, Jpn. J. Appl. Phys. 7, 31 (1968)", "doi": "10.1143/JJAP.7.31",
        "space_group": "C2/m (b unique, b the chain axis)", "number": 12,
        "cell": (8.09, 2.53, 4.79, 90, 107.9, 90), "chain_axis": "b", "Z": 2, "chains": 2,
        "density": 0.997,
        "ops": ["x,y,z", "-x,y,-z", "-x,-y,-z", "x,-y,z", "1/2+x,1/2+y,z", "1/2-x,1/2+y,-z", "1/2-x,1/2-y,-z", "1/2+x,1/2-y,z"],
        # the paper prints no coordinates: the zigzag plane is perpendicular to a (phi = 90 deg, p. 34) and the chain is
        # "the same as the usual OR form" (p. 33), so C is on 4i (1/4 + dx, 0, dz) with Bunn's C-C 1.53 A and the
        # 2.53 A repeat: the zigzag half-width sqrt(1.53^2 - 1.265^2)/2 = 0.4303 A along z (perpendicular to a)
        "sites": "pe-monoclinic",
        "where": "cell, space group and density p. 32; phi = 90 deg p. 34; no coordinate table in the paper",
        "hydrogens": {"C1": "CH2"},
        "corrections": ["the paper prints no atomic coordinates: C placed from its chain setting (zigzag plane perpendicular to a, "
                        "phi = 90 deg, p. 34) with the orthorhombic chain geometry it assumes (Bunn 1939: C-C 1.53 A)"],
        "use": "PE crystal polymorph (strained or drawn PE)", "hkl": (0, 0, 1),
    },
    {
        "id": "ipp-alpha", "name": "Isotactic polypropylene, α (α2, ordered)", "polymer": "isotactic polypropylene (iPP)", "repeat": "C3H6",
        "source": "M. Hikosaka & T. Seto, Polym. J. 5, 111 (1973)", "doi": "10.1295/polymj.5.111",
        "space_group": "P21/c (b unique, c the chain axis)", "number": 14,
        "cell": (6.65, 20.73, 6.50, 90, 98.6667, 90), "chain_axis": "c", "Z": 12, "chains": 4,
        "density": None,
        "density_note": "the paper gives no density; 12 C3H6 in the specimen B cell give 0.947 g/cm3",
        "ops": ["x,y,z", "-x,1/2+y,1/2-z", "-x,-y,-z", "x,1/2-y,1/2+z"],
        "sites": [("C1", "C", -0.0423, -0.1759, 0.4699), ("C2", "C", 0.1669, -0.1527, 0.5805), ("C3", "C", 0.1680, -0.1478, 0.8171),
                  ("C4", "C", 0.5423, -0.1759, 0.8968), ("C5", "C", 0.3744, -0.1300, 0.9470), ("C6", "C", 0.3606, -0.1321, 0.1813),
                  ("C7", "C", 0.2500, -0.0152, 0.1833), ("C8", "C", 0.2087, -0.0843, 0.2538), ("C9", "C", 0.2214, -0.0870, 0.4923)],
        "where": "cells p. 112 (specimen B: a 6.65, b 20.73, c 6.50 A, beta 98 deg 40 min); Table Ia p. 115 (carbon only)",
        "hydrogens": {"C1": "CH3", "C4": "CH3", "C7": "CH3", "C2": "CH", "C5": "CH", "C8": "CH", "C3": "CH2", "C6": "CH2", "C9": "CH2"},
        "corrections": ["the appendix's equivalent positions (pp. 126-127) do not fit Table Ia; the standard P21/c positions "
                        "(x,y,z; -x,1/2+y,1/2-z; -x,-y,-z; x,1/2-y,1/2+z) do: every C-C 1.54 A and the paper's structure factors "
                        "of Table Ib are reproduced",
                        "the ordered limiting structure alpha2 (all chains up/down as in P21/c); the alpha1 form is the statistical "
                        "C2/c up/down model of the same coordinates"],
        "use": "iPP crystal", "hkl": (0, 1, 0),
    },
    {
        "id": "peo", "name": "Poly(ethylene oxide), monoclinic", "polymer": "poly(ethylene oxide) (PEO)", "repeat": "C2H4O",
        "source": "Y. Takahashi & H. Tadokoro, Macromolecules 6, 672 (1973)", "doi": "10.1021/ma60035a005",
        "space_group": "P21/a (b unique, c the chain axis)", "number": 14,
        "cell": (8.05, 13.04, 19.48, 90, 125.4, 90), "chain_axis": "c", "Z": 28, "chains": 4,
        "density": None,
        "density_note": "the paper gives no density; 28 C2H4O in the paper's cell give 1.229 g/cm3",
        "ops": ["x,y,z", "1/2-x,1/2+y,-z", "-x,-y,-z", "1/2+x,1/2-y,z"],
        "sites": "peo",
        "where": "cell and space group p. 672; Table II p. 674 (all atoms, hydrogens included)",
        "hydrogens": "paper",
        "use": "PEO crystal (7/2 helix)", "hkl": (1, 2, 0),
    },
    {
        "id": "pvdf-beta", "name": "Poly(vinylidene fluoride), β (form I)", "polymer": "poly(vinylidene fluoride) (PVDF)", "repeat": "C2H2F2",
        "source": "R. Hasegawa, Y. Takahashi, Y. Chatani & H. Tadokoro, Polym. J. 3, 600 (1972)", "doi": "10.1295/polymj.3.600",
        "space_group": "Cm2m (c the chain axis)", "number": 38,
        "cell": (8.58, 4.91, 2.56, 90, 90, 90), "chain_axis": "c", "Z": 2, "chains": 2,
        "density": 1.973,
        "ops": ["x,y,z", "-x,y,z", "x,y,-z", "-x,y,-z", "1/2+x,1/2+y,z", "1/2-x,1/2+y,z", "1/2+x,1/2+y,-z", "1/2-x,1/2+y,-z"],
        "sites": [("C1", "C", 0.0, 0.0, 0.0), ("C2", "C", 0.0, 0.174, 0.5), ("F2", "F", 0.126, 0.335, 0.5), ("H1", "H", 0.105, -0.124, 0.0)],
        "where": "cell p. 601; Table I(a) p. 603 (planar zigzag model, hydrogens included); density p. 604",
        "hydrogens": "paper",
        "corrections": ["F y/b = 0.335 instead of the printed 0.355 (Table I(a)): 0.355 gives C-F 1.40 A and F-C-F 101 deg against the "
                        "paper's stated C-F 1.34 A and F-C-F 108 deg (p. 603), which give 0.335; the printed form III value (0.330) "
                        "and the disordered model's F pair (mean 0.332) agree with 0.335",
                        "the ordered planar zigzag model (a), not the statistical +-7 deg model (b)"],
        "use": "piezoelectric PVDF crystal", "hkl": (1, 1, 0),
    },
    {
        "id": "pvdf-alpha", "name": "Poly(vinylidene fluoride), α (form II)", "polymer": "poly(vinylidene fluoride) (PVDF)", "repeat": "C2H2F2",
        "source": "R. Hasegawa, Y. Takahashi, Y. Chatani & H. Tadokoro, Polym. J. 3, 600 (1972)", "doi": "10.1295/polymj.3.600",
        "space_group": "P21/c (b unique, c the chain axis)", "number": 14,
        "cell": (4.96, 9.64, 4.62, 90, 90, 90), "chain_axis": "c", "Z": 4, "chains": 2,
        "density": 1.925,
        "ops": ["x,y,z", "-x,1/2+y,1/2-z", "-x,-y,-z", "x,1/2-y,1/2+z"],
        "sites": [("C1", "C", 0.238, 0.183, 0.825), ("C2", "C", 0.308, 0.186, 0.150), ("F1", "F", 0.218, 0.070, 0.276),
                  ("F2", "F", 0.576, 0.178, 0.185), ("H1", "H", 0.366, 0.096, 0.760), ("H2", "H", 0.021, 0.166, 0.806)],
        "where": "cell and density p. 604; Table III p. 605 (hydrogens included)",
        "hydrogens": "paper",
        "use": "PVDF crystal (TGTG', the common melt-crystallised form)", "hkl": (0, 1, 0),
    },
    {
        "id": "pom-orthorhombic", "name": "Polyoxymethylene, orthorhombic", "polymer": "polyoxymethylene (POM)", "repeat": "CH2O",
        "source": "G. Carazzolo & M. Mammi, J. Polym. Sci. A 1, 965 (1963)", "doi": "10.1002/pol.1963.100010311",
        "space_group": "P212121 (c the chain axis)", "number": 19,
        "cell": (4.77, 7.65, 3.56, 90, 90, 90), "chain_axis": "c", "Z": 4, "chains": 2,
        "density": 1.54,
        "ops": ["x,y,z", "1/2-x,-y,1/2+z", "-x,1/2+y,1/2-z", "1/2+x,1/2-y,-z"],
        "sites": [("C1", "C", 0.367, 0.073, 0.875), ("O1", "O", 0.133, 0.073, 0.625)],
        "where": "cell and density p. 968; space group p. 969; Table I p. 974 (C and O)",
        "hydrogens": {"C1": "CH2"},
        "corrections": ["this paper's structure is the orthorhombic (2/1 helix) form; it gives no coordinates for the common trigonal "
                        "9/5 helix form"],
        "use": "POM crystal (orthorhombic form)", "hkl": (1, 0, 0),
    },
    {
        "id": "nylon6-alpha", "name": "Nylon 6, α", "polymer": "polyamide 6 (nylon 6)", "repeat": "C6H11NO",
        "source": "D. R. Holmes, C. W. Bunn & D. J. Smith, J. Polym. Sci. 17, 159 (1955)", "doi": "10.1002/pol.1955.120178401",
        "space_group": "P21 (b unique, b the chain axis)", "number": 4,
        "cell": (9.56, 17.24, 8.01, 90, 67.5, 90), "chain_axis": "b", "Z": 8, "chains": 4,
        "density": 1.23,
        "ops": ["x,y,z", "-x,1/2+y,-z"],
        "sites": "nylon6",
        "where": "cell, space group and density p. 161; Table II p. 172 (heavy atoms)",
        "hydrogens": {"CH2": "CH2", "N": "NH"},
        "corrections": ["molecule A's O z = -0.027 instead of the printed -0.007 (C=O 1.15 A as printed): molecule B prints -0.027, "
                        "and molecule C, which the paper says is A displaced by 3/14 b and c/2 (p. 163), prints 0.473 = 0.5 - 0.027",
                        "molecules B and D: the CH2 zigzag phase swapped (x 0.544 <-> 0.460, z +-0.007 <-> -+0.007) so it fits their "
                        "own N, C and O (as printed N-CH2 is 1.26 A and C(O)-CH2 1.19 A; swapped 1.48 and 1.44 A)"],
        "prerelax": True,
        "use": "nylon 6 crystal (hydrogen-bonded sheets)", "hkl": (0, 0, 1),
    },
    {
        "id": "nylon66-alpha", "name": "Nylon 66, α", "polymer": "polyamide 66 (nylon 66)", "repeat": "C12H22N2O2",
        "source": "C. W. Bunn & E. V. Garner, Proc. R. Soc. Lond. A 189, 39 (1947)", "doi": "10.1098/rspa.1947.0028",
        "space_group": "P-1 (c the chain axis)", "number": 2,
        "cell": (4.9, 5.4, 17.2, 48.5, 77, 63.5), "chain_axis": "c", "Z": 1, "chains": 1,
        "density": 1.24,
        "ops": ["x,y,z", "-x,-y,-z"],
        "sites": [("CH2a", "C", -0.083, 0.041, 0.037), ("CH2b", "C", 0.120, -0.040, 0.115), ("C", "C", -0.088, 0.010, 0.186),
                  ("O", "O", -0.371, 0.010, 0.201), ("N", "N", 0.058, -0.019, 0.253), ("CH2c", "C", -0.131, 0.095, 0.313),
                  ("CH2d", "C", 0.089, -0.039, 0.392), ("CH2e", "C", -0.106, 0.065, 0.458)],
        "where": "cell p. 39 and p. 49; density p. 42; coordinates p. 49 (heavy atoms)",
        "hydrogens": {"CH2a": "CH2", "CH2b": "CH2", "CH2c": "CH2", "CH2d": "CH2", "CH2e": "CH2", "N": "NH"},
        "use": "nylon 66 crystal", "hkl": (0, 1, 0),
    },
    {
        "id": "pet", "name": "Poly(ethylene terephthalate), triclinic", "polymer": "poly(ethylene terephthalate) (PET)", "repeat": "C10H8O4",
        "source": "R. de P. Daubeny, C. W. Bunn & C. J. Brown, Proc. R. Soc. Lond. A 226, 531 (1954)", "doi": "10.1098/rspa.1954.0273",
        "space_group": "P-1 (c near the chain axis)", "number": 2,
        "cell": (4.56, 5.94, 10.75, 98.5, 118, 112), "chain_axis": "c", "Z": 1, "chains": 1,
        "density": 1.455,
        "ops": ["x,y,z", "-x,-y,-z"],
        "sites": [("C1", "C", 0.000, -0.078, 0.371), ("C2", "C", 0.000, -0.181, 0.237), ("C3", "C", 0.100, 0.170, 0.438),
                  ("C4", "C", -0.117, -0.257, 0.433), ("C5", "C", 0.045, -0.085, 0.040), ("O6", "O", -0.098, -0.415, 0.170),
                  ("O7", "O", 0.052, 0.0, 0.175)],
        "where": "cell p. 531 and p. 533; density p. 536; Table 1 p. 538 (heavy atoms)",
        "hydrogens": {"C3": "CHar", "C4": "CHar", "C5": "CH2"},
        "use": "PET crystal", "hkl": (1, 0, 0),
    },
    {
        "id": "plla-alpha", "name": "Poly(L-lactide), α", "polymer": "poly(L-lactic acid) (PLLA)", "repeat": "C3H4O2",
        "source": "S. Sasaki & T. Asakura, Macromolecules 36, 8385 (2003)", "doi": "10.1021/ma0348674",
        "space_group": "P212121 (c the chain axis)", "number": 19,
        "cell": (10.66, 6.16, 28.88, 90, 90, 90), "chain_axis": "c", "Z": 20, "chains": 2,
        "density": 1.26,
        "ops": ["x,y,z", "-x,-y,1/2+z", "1/2+x,1/2-y,-z", "1/2-x,1/2+y,1/2-z"],
        "sites": "plla",
        "where": "cell, space group and density p. 8388; Table 1 p. 8389 (heavy atoms and the methine H)",
        "hydrogens": {"Cb": "CH3"},
        "use": "PLA crystal (10/3 helix, antiparallel chains)", "hkl": (1, 1, 0),
    },
    {
        "id": "natural-rubber", "name": "Natural rubber (cis-1,4-polyisoprene)", "polymer": "natural rubber (cis-1,4-polyisoprene)", "repeat": "C5H8",
        "source": "Y. Takahashi & T. Kumano, Macromolecules 37, 4860 (2004)", "doi": "10.1021/ma0305683",
        "space_group": "P21/a (b unique, c the chain axis)", "number": 14,
        "cell": (12.41, 8.81, 8.23, 90, 93.1, 90), "chain_axis": "c", "Z": 8, "chains": 4,
        "density": None,
        "density_note": "the paper gives no density; 8 C5H8 in the paper's cell give 1.007 g/cm3",
        "ops": ["x,y,z", "-x,-y,-z", "1/2+x,1/2-y,z", "1/2-x,1/2+y,-z"],
        "sites": "nr",
        "where": "cell p. 4860 and Fig. 5 p. 4862 (beta 93.1 deg; the abstract prints 94.6); Table 2 p. 4862 (all atoms)",
        "hydrogens": "paper",
        "corrections": ["beta = 93.1 deg (text and Fig. 5) rather than the abstract's 94.6 deg: the coordinates give regular bond lengths "
                        "with 93.1 deg only",
                        "ordered model: every site holds the Table 2 molecule. The paper refines a statistical structure in which the "
                        "mirror image about y = 1/8 occupies a site with 33 % weight; a DFT cell needs one orientation"],
        "use": "strain-crystallised natural rubber", "hkl": (0, 1, 0),
    },
    {
        "id": "ptfe-iv", "name": "Polytetrafluoroethylene, phase IV (15/7 helix)", "polymer": "polytetrafluoroethylene (PTFE)", "repeat": "CF2",
        "source": "E. S. Clark & L. T. Muus, Z. Kristallogr. 117, 119 (1962)", "doi": "10.1524/zkri.1962.117.2-3.119",
        "space_group": "P31 (c the chain axis)", "number": 144,
        "cell": (5.66, 5.66, 19.50, 90, 90, 120), "chain_axis": "c", "Z": 15, "chains": 1,
        "density": None,
        "density_note": "the paper gives no density; 15 CF2 in the paper's hexagonal cell give 2.302 g/cm3",
        "ops": ["x,y,z", "-y,x-y,1/3+z", "y-x,-x,2/3+z"],
        "sites": "ptfe-iv", "helix": (15, 7),
        "where": "hexagonal cell (25 C) p. 121; 15/7 helix packed in P31, one chain per cell, p. 127; no atomic coordinates in the paper",
        "hydrogens": "none",
        "symmetry_note": ["at this setting angle the chain's own twofold axes fall on the lattice's, so the cell also has them: "
                          "P3112 (No. 151, origin at z = 11/30), which contains the paper's P31"],
        "use": "PTFE crystal between 19 and 30 C (room temperature), the form most used in simulations; VASP and phono3py", "hkl": (1, 0, 0),
    },
    {
        "id": "ptfe-ii", "name": "Polytetrafluoroethylene, phase II (13/6 helix)", "polymer": "polytetrafluoroethylene (PTFE)", "repeat": "CF2",
        "source": "E. S. Clark & L. T. Muus, Z. Kristallogr. 117, 119 (1962)", "doi": "10.1524/zkri.1962.117.2-3.119",
        "space_group": "P1 (c the chain axis)", "number": 1,
        "cell": (5.59, 5.59, 16.88, 90, 90, 119.3), "chain_axis": "c", "Z": 13, "chains": 1,
        "density": None,
        "density_note": "the paper gives no density; 13 CF2 in the paper's oblique cell give 2.347 g/cm3",
        "ops": ["x,y,z"],
        "sites": "ptfe-ii", "helix": (13, 6),
        "where": "oblique (probably triclinic) cell at 0 C, a' = b' = 5.59 A, c = 16.88 A, gamma = 119.3 deg, p. 121; 13/6 helix on "
                 "that lattice pp. 123 and 126; no atomic coordinates in the paper",
        "hydrogens": "none",
        "corrections": ["the paper treats the lattice as the periodic arrangement of the chain axes (p. 121) and gives no space group: "
                        "one chain per cell, written as P1 with alpha = beta = 90 deg (c along the chain axis)"],
        "use": "PTFE crystal below 19 C", "hkl": (1, 0, 0),
    },
]

# ------------------------------------------------------------------------------------------------- PTFE helices
# Clark & Muus give the lattices and the helix symmetry but no coordinates, so the chain is generated: N CF2 in T turns
# span exactly c, carbon k sits at (r cos(psi + k phi), r sin(psi + k phi), k c / N) with phi = 2 pi T / N, and r follows
# from the C-C bond. Bond lengths and F-C-F are those Hasegawa, Takahashi, Chatani & Tadokoro used for their fluoropolymer
# chain (Polym. J. 3, 600 (1972), p. 603): C-C 1.54 A, C-F 1.34 A, F-C-F 108 deg. The C-C-C angle then follows from the
# helix (115.6 deg for 15/7; Bunn & Howells, Nature 174, 549 (1954), p. 551, found 116 deg). Each F pair lies in the
# plane that bisects C-C-C, symmetric about it (the helix's own twofold axis through every C). The setting angle psi of
# the chain about its axis is not in the paper: it is the one, in 1 deg steps, that gives the longest shortest F...F
# contact between neighbouring chains.
HELIX_CC, HELIX_CF, HELIX_FCF = 1.54, 1.34, 108.0
HELIX_REF = "Hasegawa, Takahashi, Chatani & Tadokoro, Polym. J. 3, 600 (1972), p. 603"


def helix_chain(cell, units, turns, psi_deg):
    L = lattice(*cell)
    h = cell[2] / units
    phi = 2 * math.pi * turns / units
    r = math.sqrt((HELIX_CC ** 2 - h ** 2) / (2 * (1 - math.cos(phi))))
    psi = math.radians(psi_deg)
    carbons = [[r * math.cos(psi + k * phi), r * math.sin(psi + k * phi), k * h] for k in range(-1, units + 1)]
    half = math.radians(HELIX_FCF) / 2
    sites = []
    for k in range(1, units + 1):
        c = carbons[k]
        u1, u2 = unit(sub(carbons[k - 1], c)), unit(sub(carbons[k + 1], c))
        bis = unit(mul(add(u1, u2), -1))
        n = unit(cross(u1, u2))
        sites.append((f"C{k}", "C", c))
        for s in (1, -1):
            sites.append((f"F{k}", "F", add(c, mul(add(mul(bis, math.cos(half)), mul(n, s * math.sin(half))), HELIX_CF))))
    return L, r, h, phi, [(lab, el, *to_frac(L, p)) for lab, el, p in sites]    # unwrapped: the chain stays whole


def interchain_ff(L, sites):
    """Shortest F...F between the chain and its neighbours (the chain itself is one cell wide)."""
    fs = [to_cart(L, s[2:]) for s in sites if s[1] == "F"]
    best = 9.0
    for i in (-1, 0, 1):
        for j in (-1, 0, 1):
            if i or j:
                for k in (-1, 0, 1):
                    t = to_cart(L, [i, j, k])
                    moved = [add(q, t) for q in fs]
                    best = min(best, min(math.dist(p, q) for p in fs for q in moved))
    return best


HELIX = {}


def ptfe_sites(x):
    units, turns = x["helix"]
    best = None
    for i in range(360):
        psi = float(i)
        L, r, h, phi, sites = helix_chain(x["cell"], units, turns, psi)
        d = interchain_ff(L, sites)
        if best is None or d > best[0] + 1e-9:
            best = (d, psi)
    L, r, h, phi, sites = helix_chain(x["cell"], units, turns, best[1])
    s13 = math.sqrt(2 * r * r * (1 - math.cos(2 * phi)) + 4 * h * h)     # C(k-1)...C(k+1)
    ccc = math.degrees(math.acos(1 - s13 ** 2 / (2 * HELIX_CC ** 2)))
    fr = [math.hypot(*to_cart(L, s[2:])[:2]) for s in sites if s[1] == "F"]
    HELIX[x["id"]] = {"units": units, "turns": turns, "rise": round(h, 6), "twist_deg": round(math.degrees(phi), 4),
                      "C_radius": round(r, 4), "F_radius": round(sum(fr) / len(fr), 4), "C-C": HELIX_CC, "C-F": HELIX_CF,
                      "F-C-F": HELIX_FCF, "C-C-C": round(ccc, 3), "setting_angle_deg": best[1],
                      "interchain_FF": round(best[0], 3), "bonds_from": HELIX_REF}
    x["corrections"] = list(x.get("corrections", [])) + [
        f"no atomic coordinates in the paper: the chain is generated from the {units}/{turns} helix and the cited lattice, not refined: "
        f"{units} CF2 in {turns} turns span c exactly (rise {h:.4f} A, twist {math.degrees(phi):.2f} deg per CF2), C-C {HELIX_CC} A, "
        f"C-F {HELIX_CF} A, F-C-F {HELIX_FCF:g} deg ({HELIX_REF}); this gives C-C-C {ccc:.2f} deg, carbon helix radius {r:.3f} A "
        f"and fluorine radius {HELIX[x['id']]['F_radius']:.3f} A (Bunn & Howells, Nature 174, 549 (1954): 116 deg, about 0.42 A "
        f"and 1.64 +- 0.05 A); one chain on the cell axis, setting angle {best[1]:g} deg (the paper gives none: the one with "
        f"the longest shortest F...F between chains, {best[0]:.3f} A)"] + x.get("symmetry_note", [])
    return [(lab, el, *wrap([fx, fy, fz])) for lab, el, fx, fy, fz in sites]


def sites_named(key):
    if key == "pe-monoclinic":
        a, b, c, al, be, ga = 8.09, 2.53, 4.79, 90, 107.9, 90
        half = math.sqrt(1.53 ** 2 - (b / 2) ** 2) / 2         # zigzag half-width, A
        bet = math.radians(be)
        dz = half / (c * math.sin(bet))
        dx = -dz * c * math.cos(bet) / a                         # keeps the offset perpendicular to a
        return [("C1", "C", 0.25 + dx, 0.0, dz)]
    if key == "peo":
        rows = """O1 0.064 0.366 -0.042|C1 -0.070 0.313 -0.028|C2 0.055 0.250 0.053|O2 0.194 0.315 0.123|C3 0.345 0.257 0.195
C4 0.469 0.329 0.272|O3 0.343 0.367 0.296|C5 0.449 0.438 0.364|C6 0.302 0.485 0.382|O4 0.254 0.412 0.423|C7 0.063 0.433 0.409
C8 0.025 0.355 0.458|O5 0.201 0.347 0.543|C9 0.172 0.275 0.591|C10 0.340 0.289 0.685|O6 0.274 0.359 0.722|C11 0.433 0.382 0.807
C12 0.355 0.457 0.843|O7 0.188 0.413 0.840|C13 0.105 0.481 0.870|C14 -0.045 0.423 0.881
H1 -0.162 0.369 -0.023|H2 -0.169 0.263 -0.081|H3 -0.048 0.214 0.065|H4 0.140 0.191 0.046|H5 0.270 0.199 0.209|H6 0.446 0.220 0.183
H7 0.596 0.287 0.324|H8 0.526 0.393 0.256|H9 0.571 0.398 0.420|H10 0.514 0.499 0.348|H11 0.373 0.552 0.422|H12 0.162 0.508 0.322
H13 0.066 0.510 0.431|H14 -0.058 0.427 0.342|H15 -0.105 0.380 0.457|H16 -0.008 0.280 0.428|H17 0.024 0.288 0.579|H18 0.179 0.198 0.572
H19 0.374 0.215 0.717|H20 0.477 0.319 0.694|H21 0.484 0.312 0.844|H22 0.561 0.417 0.810|H23 0.477 0.475 0.909|H24 0.304 0.528 0.806
H25 0.228 0.512 0.931|H26 0.026 0.543 0.825|H27 -0.138 0.371 0.828|H28 -0.144 0.478 0.883"""
        return table(rows)
    if key == "nr":
        rows = """C1 0.248 0.229 0.180|C2 0.156 0.193 0.291|C3 0.048 0.222 0.199|C4 0.161 0.142 0.446|C5 0.270 0.113 0.537
C6 0.289 0.221 0.683|C7 0.370 0.151 0.809|C8 0.480 0.138 0.734|C9 0.356 0.104 0.961|C10 0.246 0.117 1.035
H1 0.029 0.336 0.152|H2 0.041 0.151 0.091|H3 -0.018 0.197 0.277|H4 0.080 0.129 0.487|H5 0.275 -0.006 0.574|H6 0.335 0.138 0.459
H7 0.325 0.220 0.249|H8 0.237 0.344 0.133|H9 0.521 0.240 0.693|H10 0.473 0.066 0.626|H11 0.538 0.088 0.823|H12 0.432 0.060 1.015
H13 0.218 0.005 1.069|H14 0.188 0.167 0.946|H15 0.213 0.242 0.740|H16 0.321 0.328 0.639"""
        return table(rows)
    if key == "plla":
        rows = """O''1 0.0036 0.1115 0.0569|Ca1 0.0891 0.0603 0.0948|Cb1 0.1350 0.2722 0.1175|Ha1 0.1686 -0.0290 0.0816|C'1 0.0206 -0.0776 0.1309|O'1 -0.0912 -0.0977 0.1335
O''2 0.1018 -0.1705 0.1586|Ca2 0.0549 -0.2709 0.2010|Cb2 0.1641 -0.3802 0.2269|Ha2 -0.0149 -0.3914 0.1923|C'2 -0.0041 -0.0970 0.2319|O'2 0.0226 0.0922 0.2313
O''3 -0.0881 -0.1851 0.2592|Ca3 -0.1321 -0.0589 0.2988|Cb3 -0.2262 -0.1948 0.3270|Ha3 -0.1774 0.0873 0.2867|C'3 -0.0204 0.0016 0.3295|O'3 0.0733 -0.1046 0.3337
O''4 -0.0421 0.1863 0.3504|Ca4 0.0423 0.2510 0.3878|Cb4 0.0086 0.4817 0.4041|Ha4 0.1379 0.2489 0.3754|C'4 0.0291 0.0925 0.4283|O'4 -0.0600 -0.0222 0.4351
O''5 0.1294 0.1005 0.4544|Ca5 0.1324 -0.0352 0.4959|Cb5 0.2643 -0.0269 0.5175|Ha5 0.1100 -0.2007 0.4867|C'5 0.0362 0.0486 0.5309|O'5 -0.0002 0.2326 0.5333"""
        return table(rows)
    if key == "nylon6":
        # Table II as printed, with the two corrections of the entry: A's O z -0.027, and B/D's CH2 zigzag phase swapped
        mol = {
            "A": [("N", 0.037, 0.285, 0.006), ("C", -0.040, 0.219, -0.007), ("O", -0.160, 0.219, -0.027), ("CH2", 0.044, 0.146, 0.007),
                  ("CH2", -0.040, 0.073, -0.007), ("CH2", 0.044, 0.0, 0.007), ("CH2", -0.040, -0.073, -0.007), ("CH2", 0.044, -0.146, 0.007)],
            "B": [("N", 0.537, 0.219, 0.006), ("C", 0.460, 0.285, -0.007), ("O", 0.340, 0.285, -0.027), ("CH2", 0.460, 0.146, -0.007),
                  ("CH2", 0.544, 0.073, 0.007), ("CH2", 0.460, 0.0, -0.007), ("CH2", 0.544, -0.073, 0.007), ("CH2", 0.460, -0.146, -0.007)],
            "C": [("N", 0.037, 0.499, 0.506), ("C", -0.040, 0.433, 0.493), ("O", -0.160, 0.433, 0.473), ("CH2", 0.044, 0.360, 0.507),
                  ("CH2", -0.040, 0.287, 0.493), ("CH2", 0.044, 0.214, 0.507), ("CH2", -0.040, 0.141, 0.493), ("CH2", 0.044, 0.068, 0.507)],
            "D": [("N", 0.537, 0.433, 0.506), ("C", 0.460, 0.499, 0.493), ("O", 0.340, 0.499, 0.473), ("CH2", 0.460, 0.360, 0.493),
                  ("CH2", 0.544, 0.287, 0.507), ("CH2", 0.460, 0.214, 0.493), ("CH2", 0.544, 0.141, 0.507), ("CH2", 0.460, 0.068, 0.493)],
        }
        out = []
        for m, rows in mol.items():
            for lab, x, y, z in rows:
                out.append((lab, "C" if lab in ("C", "CH2") else lab, x, y, z))
        return out
    raise KeyError(key)


def table(rows):
    out = []
    for item in rows.replace("\n", "|").split("|"):
        lab, x, y, z = item.split()
        el = lab[0] if lab[0] in "CHNOF" else "C"
        if lab.startswith("Ha"):
            el = "H"
        out.append((lab, el, float(x), float(y), float(z)))
    return out


# ------------------------------------------------------------------------------------------------------ vector helpers
def add(a, b): return [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
def sub(a, b): return [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
def mul(a, s): return [a[0] * s, a[1] * s, a[2] * s]
def dot(a, b): return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
def cross(a, b): return [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
def norm(a): return math.sqrt(dot(a, a))
def unit(a): n = norm(a); return [a[0] / n, a[1] / n, a[2] / n]


def lattice(a, b, c, al, be, ga):
    """a along x, b in the xy plane (as caps::cell_parameters)."""
    al, be, ga = (math.radians(v) for v in (al, be, ga))
    ax = [a, 0, 0]
    bx = [b * math.cos(ga), b * math.sin(ga), 0]
    cx = c * math.cos(be)
    cy = c * (math.cos(al) - math.cos(be) * math.cos(ga)) / math.sin(ga)
    cz = math.sqrt(max(c * c - cx * cx - cy * cy, 0))
    return [ax, bx, [cx, cy, cz]]


def to_cart(L, f): return [f[0] * L[0][k] + f[1] * L[1][k] + f[2] * L[2][k] for k in range(3)]


def to_frac(L, r):
    # solve r = f0 a + f1 b + f2 c
    a, b, c = L
    v = dot(a, cross(b, c))
    return [dot(r, cross(b, c)) / v, dot(a, cross(r, c)) / v, dot(a, cross(b, r)) / v]


def parse_op(s):
    """'1/2-x,1/2+y,-z' -> (R, t) with t a Fraction triple."""
    R, t = [], []
    for part in s.replace(" ", "").split(","):
        row, const, sign, num = [0, 0, 0], Fraction(0), 1, ""
        part += "+"
        for ch in part:
            if ch in "+-":
                if num:
                    const += sign * Fraction(num)
                    num = ""
                sign = 1 if ch == "+" else -1
            elif ch in "xyz":
                row["xyz".index(ch)] += sign
                num = ""
            else:
                num += ch
        R.append(row)
        t.append(const)
    return R, t


def apply(op, f):
    R, t = op
    return [sum(R[i][k] * f[k] for k in range(3)) + float(t[i]) for i in range(3)]


def close_group(ops):
    group = [parse_op(o) for o in ops]
    key = lambda op: (tuple(map(tuple, op[0])), tuple(x % 1 for x in op[1]))
    seen = {key(o): o for o in group}
    changed = True
    while changed:
        changed = False
        for a in list(seen.values()):
            for b in list(seen.values()):
                R = [[sum(a[0][i][k] * b[0][k][j] for k in range(3)) for j in range(3)] for i in range(3)]
                t = [sum(a[0][i][k] * b[1][k] for k in range(3)) + a[1][i] for i in range(3)]
                op = (R, [x % 1 for x in t])
                if key(op) not in seen:
                    seen[key(op)] = op
                    changed = True
    return list(seen.values())


def wrap(f): return [x - math.floor(x) for x in f]


def frac_dist(L, f, g):
    d = [f[k] - g[k] for k in range(3)]
    d = [x - round(x) for x in d]
    return norm(to_cart(L, d))


def expand(L, ops, sites, tol=0.01):
    atoms = []
    for lab, el, x, y, z in sites:
        for op in ops:
            f = wrap(apply(op, [x, y, z]))
            if any(a["el"] == el and frac_dist(L, a["f"], f) < tol for a in atoms):
                continue
            atoms.append({"label": lab, "el": el, "f": f})
    return atoms


# ------------------------------------------------------------------------------------------------- neighbours in PBC
def images(L, cut):
    # enough images along each axis for the cutoff (thin cells such as PE's 2.5 A repeat need several)
    a, b, c = L
    v = abs(dot(a, cross(b, c)))
    h = [v / norm(cross(b, c)), v / norm(cross(c, a)), v / norm(cross(a, b))]
    n = [int(math.ceil(cut / h[k])) + 1 for k in range(3)]
    return [(i, j, k) for i in range(-n[0], n[0] + 1) for j in range(-n[1], n[1] + 1) for k in range(-n[2], n[2] + 1)]


def neighbours(L, atoms, cut):
    """[(i, j, vector j-i in A)] for every pair (and image) closer than cut."""
    out = [[] for _ in atoms]
    ims = images(L, cut)
    cart = [to_cart(L, a["f"]) for a in atoms]
    for i in range(len(atoms)):
        for j in range(len(atoms)):
            for im in ims:
                if i == j and im == (0, 0, 0):
                    continue
                r = add(sub(cart[j], cart[i]), to_cart(L, im))
                d = norm(r)
                if d < cut:
                    out[i].append((j, r, d))
    return out


BOND = {("C", "C"): 1.75, ("C", "O"): 1.65, ("C", "N"): 1.65, ("C", "F"): 1.6, ("C", "H"): 1.25, ("N", "H"): 1.2, ("O", "H"): 1.2}


def bonded(e1, e2, d):
    k = (e1, e2) if (e1, e2) in BOND else (e2, e1)
    return k in BOND and d < BOND[k]


def heavy_bonds(L, atoms):
    nb = neighbours(L, atoms, 2.0)
    return [[(j, r) for j, r, d in nb[i] if bonded(atoms[i]["el"], atoms[j]["el"], d)] for i in range(len(atoms))]


# ------------------------------------------------------------------------------------------------------- hydrogens
TET = math.radians(109.4712)


def place_h(L, atoms, rules):
    """Adds hydrogens to the expanded heavy atoms by the rule of each atom's site label."""
    bonds = heavy_bonds(L, atoms)
    new = []
    for i, a in enumerate(atoms):
        kind = rules.get(a["label"]) or rules.get(a["label"].rstrip("0123456789"))
        if not kind:
            continue
        c = to_cart(L, a["f"])
        heavy = [(j, r) for j, r in bonds[i] if atoms[j]["el"] != "H"]
        us = [unit(r) for _, r in heavy]
        hs = []
        if kind == "CH2":
            assert len(us) == 2, (a["label"], len(us))
            bis = unit(mul(add(us[0], us[1]), -1))
            n = unit(cross(us[0], us[1]))
            for s in (1, -1):
                hs.append(add(c, mul(add(mul(bis, math.cos(TET / 2)), mul(n, s * math.sin(TET / 2))), 1.09)))
        elif kind == "CH":
            assert len(us) == 3, (a["label"], len(us))
            hs.append(add(c, mul(unit(mul(add(add(us[0], us[1]), us[2]), -1)), 1.09)))
        elif kind in ("CHar", "NH"):
            assert len(us) == 2, (a["label"], len(us))
            d = 1.08 if kind == "CHar" else 1.01
            hs.append(add(c, mul(unit(mul(add(us[0], us[1]), -1)), d)))
        elif kind == "CH3":
            assert len(us) == 1, (a["label"], len(us))
            j, r = heavy[0]
            ua = us[0]                                      # C -> A
            # A's nearest other heavy neighbour B sets the staggering (H anti to B)
            others = sorted(((norm(r2), r2) for k, r2 in bonds[j] if atoms[k]["el"] != "H" and norm(add(r, r2)) > 0.1), key=lambda x: x[0])
            rb = others[0][1]                               # A -> B
            # a frame on the C-A axis: e1 along -ua, e2 the direction of B projected off the axis, anti = -e2
            e1 = mul(ua, -1)
            pb = sub(rb, mul(ua, dot(rb, ua)))
            e2 = unit(pb)
            e3 = cross(e1, e2)
            for phi in (180.0, 60.0, -60.0):
                p = math.radians(phi)
                # angle A-C-H = 109.47 deg: the H direction makes 180 - 109.47 deg with -ua
                ct, st = math.cos(math.pi - TET), math.sin(math.pi - TET)
                v = add(mul(e1, ct), add(mul(e2, st * math.cos(p)), mul(e3, st * math.sin(p))))
                hs.append(add(c, mul(v, 1.09)))
        else:
            raise ValueError(kind)
        for h in hs:
            new.append({"label": "H" + a["label"], "el": "H", "f": wrap(to_frac(L, h))})
    return atoms + new


# ---------------------------------------------------------------------------------------------------------- checks
VALENCE = {"C": 4, "N": 3, "O": 2, "F": 1, "H": 1}


def check(L, atoms):
    """Valences (every atom has its neighbour count: C 4, N 3, O 2, F 1, H 1), the shortest contacts beyond 1-3."""
    nb = neighbours(L, atoms, 3.2)
    bonds = [[(j, r) for j, r, d in nb[i] if bonded(atoms[i]["el"], atoms[j]["el"], d)] for i in range(len(atoms))]
    issues = []
    for i, a in enumerate(atoms):
        n = len(bonds[i])
        # sp2 C (C=O, C=C, aromatic) has three neighbours, a carbonyl O one
        want = {VALENCE[a["el"]]}
        if a["el"] == "C":
            want = {3, 4}
        if a["el"] == "O":
            want = {1, 2}
        if n not in want:
            issues.append(f"{a['label']} ({a['el']}) has {n} neighbours")
    # contacts: pairs not within two bonds (1-2 and 1-3 pairs are geometry, not contacts); hydrogen bonds X-H...Y
    # (X, Y = N, O) are listed apart, they are short by nature
    near = {}
    for i in range(len(atoms)):
        s = []
        for j, r in bonds[i]:
            s.append((j, r))
            for k, r2 in bonds[j]:
                s.append((k, add(r, r2)))
        near[i] = s
    polar = {i for i, a in enumerate(atoms) if a["el"] == "H" and any(atoms[j]["el"] in "NO" for j, _ in bonds[i])}
    mins = {"HH": 9.0, "HX": 9.0, "XX": 9.0, "hbond": 9.0}
    for i in range(len(atoms)):
        for j, r, d in nb[i]:
            if any(k == j and norm(sub(r, rk)) < 0.05 for k, rk in near[i]):
                continue
            hi, hj = atoms[i]["el"] == "H", atoms[j]["el"] == "H"
            if (i in polar and atoms[j]["el"] in "NO") or (j in polar and atoms[i]["el"] in "NO"):
                k = "hbond"
            else:
                k = "HH" if hi and hj else "XX" if not hi and not hj else "HX"
            mins[k] = min(mins[k], d)
    return issues, mins


def density(L, atoms):
    v = abs(dot(L[0], cross(L[1], L[2])))
    m = sum(MASS[a["el"]] for a in atoms)
    return m / NA / (v * 1e-24), v


def formula(atoms):
    n = {}
    for a in atoms:
        n[a["el"]] = n.get(a["el"], 0) + 1
    order = ["C", "H"] + sorted(e for e in n if e not in ("C", "H"))
    return "".join(e + (str(n[e]) if n[e] > 1 else "") for e in order if e in n), n


def repeat_counts(f):
    out, el, num = {}, "", ""
    for ch in f + "X":
        if ch.isalpha() and ch.isupper():
            if el:
                out[el] = out.get(el, 0) + (int(num) if num else 1)
            el, num = ch, ""
        elif ch.isdigit():
            num += ch
    return out


# --------------------------------------------------------------------------------------------------------- writing
def cif_text(x, L, atoms, mins, rho, vol, order=True):
    a, b, c, al, be, ga = x["cell"]
    lines = [f"# {x['name']}: {x['source']}, doi:{x['doi']}",
             f"# space group {x['space_group']} (No. {x['number']}), Z = {x['Z']} {x['repeat']} per cell; {x['where']}"]
    if "helix" in x:
        lines.append("# coordinates: generated from the helix symmetry and the cited lattice (tools/make_polymer_crystals.py), not refined")
    for s in x.get("corrections", []):
        lines.append("# correction / model: " + s)
    h = x["hydrogens"]
    if h == "none":
        lines.append("# hydrogens: none (perfluorinated chain)")
    elif h == "paper":
        lines.append("# hydrogens: the paper's")
    else:
        kinds = sorted(set(h.values()))
        lines.append("# hydrogens: none in the paper; placed by tools/make_polymer_crystals.py (" + ", ".join(kinds) +
                     "; C-H 1.09 A sp3, 1.08 A sp2, N-H 1.01 A, tetrahedral / in-plane angles)")
    if x.get("shift"):
        big, rms, before = x["shift"]
        lines.append(f"# pre-relaxed at the paper's fixed cell with OPLS-AA 2024 (caps relax, positions only) to remove clashes "
                     f"(before: H...H {before['HH']:.2f} A): heavy atoms moved {rms:.3f} A RMS, {big:.3f} A at most")
    def m(k): return f"{mins[k]:.2f} A" if mins[k] < 9 else "none below 3.2 A"
    lines.append(f"# written in P1 with every atom of the paper's cell; {rho:.4f} g/cm3, shortest contacts beyond 1-3: "
                 f"H...H {m('HH')}, H...X {m('HX')}, X...X {m('XX')}"
                 + (f", hydrogen bond H...O/N {mins['hbond']:.2f} A" if mins['hbond'] < 9 else ""))
    lines += [f"data_{x['id'].replace('-', '_')}", f"_chemical_name_common '{x['name']}'",
              "_symmetry_space_group_name_H-M 'P 1'", "_symmetry_Int_Tables_number 1",
              f"_cell_length_a {a}", f"_cell_length_b {b}", f"_cell_length_c {c}",
              f"_cell_angle_alpha {al}", f"_cell_angle_beta {be}", f"_cell_angle_gamma {ga}", f"_cell_volume {vol:.4f}",
              "loop_", "_symmetry_equiv_pos_as_xyz", "'x,y,z'",
              "loop_", "_atom_site_label", "_atom_site_type_symbol", "_atom_site_fract_x", "_atom_site_fract_y", "_atom_site_fract_z"]
    order = sorted(range(len(atoms)), key=lambda i: (["C", "N", "O", "F", "H"].index(atoms[i]["el"]), i)) if order else range(len(atoms))
    count = {}
    for i in order:
        e = atoms[i]["el"]
        count[e] = count.get(e, 0) + 1
        f = atoms[i]["f"]
        lines.append(f"{e}{count[e]} {e} {f[0]:.6f} {f[1]:.6f} {f[2]:.6f}")
    return "\n".join(lines) + "\n"


def read_p1_cif(path):
    atoms, loop = [], False
    for line in open(path):
        w = line.split()
        if w and w[0] == "_atom_site_fract_z":
            loop = True
            continue
        if loop and len(w) >= 5 and not w[0].startswith("_"):
            atoms.append({"el": w[1], "f": [float(v) for v in w[2:5]]})
    return atoms


def prerelax(x, L, atoms):
    """Atom positions minimised at the paper's fixed cell with OPLS-AA 2024 (caps relax): only to remove the clashes of the
    paper's coordinates. Returns the moved atoms and (largest, RMS) heavy-atom shift in A."""
    with tempfile.TemporaryDirectory() as td:
        src, out = os.path.join(td, "in.cif"), os.path.join(td, "out.xyz")
        open(src, "w").write(cif_text(x, L, atoms, {"HH": 9, "HX": 9, "XX": 9, "hbond": 9}, 0, 0, order=False))
        r = subprocess.run([CAPS, "relax", src, "-o", out, "--ff", os.path.join(ROOT, "data", "forcefields", PRERELAX_FF), "--iterations", "2000", "--quiet"],
                           capture_output=True, text=True)
        if r.returncode:
            raise RuntimeError("caps relax: " + r.stderr[-400:])
        moved = []
        for line in open(out).read().splitlines()[2:]:     # XYZ in the cell's frame (a along x, b in the xy plane)
            w = line.split()
            if len(w) >= 4:
                moved.append({"el": w[0], "f": to_frac(L, [float(v) for v in w[1:4]])})
    if len(moved) != len(atoms):
        raise RuntimeError("caps relax changed the number of atoms")
    # the writer groups atoms by element: match each atom to the nearest moved one of its element
    taken, match = set(), []
    for a in atoms:
        d, k = min((frac_dist(L, a["f"], m["f"]), k) for k, m in enumerate(moved) if m["el"] == a["el"] and k not in taken)
        if d > (1.0 if a["el"] == "H" else 0.6):
            raise RuntimeError(f"caps relax moved an atom {d:.2f} A: not a clash removal")
        taken.add(k)
        match.append(moved[k])
    moved = match
    shifts = [frac_dist(L, a["f"], m["f"]) for a, m in zip(atoms, moved) if a["el"] != "H"]
    out = [dict(a, f=wrap(m["f"])) for a, m in zip(atoms, moved)]
    return out, max(shifts), math.sqrt(sum(d * d for d in shifts) / len(shifts))


def build(x):
    L = lattice(*x["cell"])
    ops = close_group(x["ops"])
    if "helix" in x:
        sites = ptfe_sites(x)
    else:
        sites = sites_named(x["sites"]) if isinstance(x["sites"], str) else x["sites"]
    atoms = expand(L, ops, sites)
    if x["hydrogens"] not in ("paper", "none"):
        atoms = place_h(L, atoms, x["hydrogens"])
    x["shift"] = None
    if x.get("prerelax"):
        before = check(L, atoms)[1]
        atoms, big, rms = prerelax(x, L, atoms)
        x["shift"] = (big, rms, before)
    issues, mins = check(L, atoms)
    rho, vol = density(L, atoms)
    return L, ops, atoms, issues, mins, rho, vol


def main():
    write = "--write" in sys.argv
    only = [a for a in sys.argv[1:] if not a.startswith("--")]
    meta = []
    bad = 0
    for x in CRYSTALS:
        if only and x["id"] not in only:
            continue
        L, ops, atoms, issues, mins, rho, vol = build(x)
        f, n = formula(atoms)
        rep = repeat_counts(x["repeat"])
        z_ok = all(n.get(e, 0) == rep[e] * x["Z"] for e in rep) and set(n) == set(rep)
        ref = x["density"]
        rho_ok = ref is None or abs(rho - ref) / ref < 0.01
        clash = mins["HH"] < 1.8 or mins["HX"] < 2.2 or mins["XX"] < 2.5
        ok = z_ok and rho_ok and not issues
        bad += not ok
        print(f"{x['id']:16s} {len(ops):2d} ops  {len(atoms):4d} atoms  {f:14s} {rho:.4f} g/cm3 (paper {ref if ref else '-'})  "
              f"H..H {mins['HH']:.2f} H..X {mins['HX']:.2f} X..X {mins['XX']:.2f} Hb {mins['hbond']:.2f}  {'CLASH ' if clash else ''}{'ok' if ok else 'CHECK ' + '; '.join(issues[:4])}")
        if not z_ok:
            print("   composition", n, "expected", {e: v * x["Z"] for e, v in rep.items()})
        if write and ok:
            os.makedirs(OUT, exist_ok=True)
            open(os.path.join(OUT, x["id"] + ".cif"), "w").write(cif_text(x, L, atoms, mins, rho, vol))
        meta.append({
            "id": x["id"], "name": x["name"], "polymer": x["polymer"], "file": x["id"] + ".cif", "source": x["source"], "doi": x["doi"],
            "where": x["where"], "space_group": x["space_group"], "space_group_number": x["number"],
            "cell": {"a": x["cell"][0], "b": x["cell"][1], "c": x["cell"][2], "alpha": x["cell"][3], "beta": x["cell"][4], "gamma": x["cell"][5]},
            "operations": x["ops"], "chain_axis": x["chain_axis"], "repeat_unit": x["repeat"], "Z": x["Z"], "chains": x["chains"], "atoms": len(atoms), "formula": f,
            "density_paper": ref, "density": round(rho, 4), **({"density_note": x["density_note"]} if "density_note" in x else {}),
            "hydrogens": x["hydrogens"] if x["hydrogens"] in ("paper", "none") else "added (tools/make_polymer_crystals.py)",
            **({"coordinates": "generated from the helix symmetry and the cited lattice (tools/make_polymer_crystals.py), not refined",
                "helix": HELIX[x["id"]]} if "helix" in x else {}),
            "corrections": x.get("corrections", []),
            "min_contacts": {k: round(v, 3) for k, v in mins.items() if v < 9},
            **({"prerelax": {"forcefield": "OPLS-AA 2024", "fixed_cell": True, "heavy_shift_max": round(x["shift"][0], 3),
                             "heavy_shift_rms": round(x["shift"][1], 3)}} if x.get("shift") else {}),
            "use": x["use"], "hkl": list(x["hkl"]),
        })
    if write and not only:
        doc = {"format": "caps-polymer-crystals",
               "about": "Polymer crystal unit cells from their original crystal-structure papers (one entry per paper and form): the paper's "
                        "cell and heavy atoms, hydrogens from the paper or placed geometrically (tools/make_polymer_crystals.py), written in P1. "
                        "PTFE's paper gives lattices and helix symmetry only: its chains are generated from the helix and the cited "
                        "lattice, not refined (see each entry's helix block). These are DFT unit cells (VASP relaxation, phono3py); "
                        "large MD cells are built from them.",
               "crystals": meta}
        open(os.path.join(OUT, "polymers.json"), "w").write(json.dumps(doc, indent=1, ensure_ascii=False) + "\n")
        print("wrote", OUT)
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
