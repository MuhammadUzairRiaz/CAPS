#!/usr/bin/env python3
"""Build the curated CAPS force-field library (data/forcefields) from moltemplate and DL_FIELD.

usage: build_library.py [--moltemplate DIR] [--dlfield LIBDIR] [--out data/forcefields]

Every source file gets a catalogue entry, whether or not it converts: its name, version / year, primary
reference, origin file, what CAPS can do with it, and the evidence behind that status. Entries whose status is
"validated" were checked interaction by interaction (and energies / forces) against their reference implementations and LAMMPS
with bench/ff/compare_moltemplate.py, compare_class2.py or compare_dlfield.py.
"""
import json, os, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
CAPS = os.path.join(ROOT, "build", "cli", "caps")
arg = lambda k, d: sys.argv[sys.argv.index(k) + 1] if k in sys.argv else d
MT = arg("--moltemplate", os.path.expanduser("~/moltemplate/moltemplate/force_fields"))
DLF = arg("--dlfield", os.path.expanduser("~/project/dl_f_4.13/lib"))
OUT = arg("--out", os.path.join(ROOT, "data", "forcefields"))

V = "validated"
C = "converted"          # parsed into CAPS rules; not yet checked against the source tool
P = "pending"            # the family's functional forms are not interpreted yet
T = "template"           # a molecule built from another force field's types, not a force field

INORG = {
    "INORGANIC_binary_halides": ("validated", "14 of 14 ionic templates: every Buckingham pair identical to DL_FIELD's FIELD file."),
    "INORGANIC_binary_misc": ("validated", "GaN: Buckingham pairs identical to DL_FIELD's FIELD file."),
    "INORGANIC_binary_oxides": ("validated", "17 of 17 buildable templates plus the alpha-Al2O3 and MgO shell-model examples: Buckingham pairs and core-shell springs identical."),
    "INORGANIC_ternary_oxides": ("validated", "12 of 12 perovskite / pyrochlore templates: Buckingham pairs and quartic core-shell springs identical (DL_FIELD's FIELD file runs large spring constants together)."),
    "INORGANIC_zeolite": ("validated", "4 of 4 silicalite / faujasite templates identical."),
}
DLF_NOTE = "DL_FIELD 4.13 library (C. W. Yong, STFC Daresbury Laboratory)"
GAFF_REF = "J. Wang, R. M. Wolf, J. W. Caldwell, P. A. Kollman, D. A. Case, J. Comput. Chem. 25, 1157 (2004)"

# id, source, file, name, version / year, references, status, evidence / notes
ENTRIES = [
    # ---- DL_FIELD ----
    ("pcff", "dlfield", "PCFF", "PCFF (DL_FIELD, diagonal terms)", "cff91 family, 1994; DL_FIELD library v4.5 (2019)",
     ["H. Sun, S. J. Mumby, J. R. Maple, A. T. Hagler, J. Am. Chem. Soc. 116, 2978 (1994)"], V,
     "matches DL_FIELD's LAMMPS output term by term and in energies / forces (≤6e-5 kcal/mol/Å) for an ester-ether, propylbenzene "
     "and a siloxane, including distorted geometries and auto (cff91_auto) parameters. No class II cross terms: DL_FIELD keeps "
     "the diagonal part only."),
    ("compass", "dlfield", "COMPASS", "COMPASS (DL_FIELD, diagonal terms)", "1998",
     ["H. Sun, J. Phys. Chem. B 102, 7338 (1998)"], V,
     "matches DL_FIELD for propylbenzene and a siloxane; DL_FIELD and CAPS both find no c4-c4-c3' angle for the ester. No cross terms."),
    ("cvff", "dlfield", "CVFF", "CVFF (DL_FIELD)", "1988",
     ["P. Dauber-Osguthorpe et al., Proteins 4, 31 (1988)"], V,
     "matches DL_FIELD for ester, propylbenzene and siloxane, including wildcard-torsion scaling and cvff_auto parameters."),
    ("opls2005", "dlfield", "OPLS2005", "OPLS-AA / OPLS 2005 (DL_FIELD)", "2005",
     ["J. L. Banks et al., J. Comput. Chem. 26, 1752 (2005)", "W. L. Jorgensen, D. S. Maxwell, J. Tirado-Rives, J. Am. Chem. Soc. 118, 11225 (1996)"], V,
     "matches DL_FIELD for ester, propylbenzene and siloxane including impropers on distorted geometries."),
    ("opls2020", "dlfield", "OPLS2020", "OPLS 2020 bond / angle supplement (DL_FIELD)", "as distributed with DL_FIELD 4.13",
     [DLF_NOTE], C, "bonds and angles only; used together with OPLS 2005."),
    ("opls-aam", "dlfield", "OPLS_AAM", "OPLS-AA/M (DL_FIELD)", "2015",
     ["M. J. Robertson, J. Tirado-Rives, W. L. Jorgensen, J. Chem. Theory Comput. 11, 3499 (2015)"], C,
     "same functional forms as OPLS 2005 (validated); DL_FIELD cannot type the test molecules with this set, so not yet compared."),
    ("opls-clp", "dlfield", "OPLS_CL_P", "CL&P ionic liquids (OPLS-AA based, DL_FIELD)", "2004 onwards",
     ["J. N. Canongia Lopes, A. A. H. Pádua, J. Phys. Chem. B 108, 2038 (2004)"], V,
     "checked term by term against DL_FIELD's DL_POLY FIELD file (bench/ff/validate_family.py): 24 of 24 identical (kJ/mol constants confirmed as DL_POLY's k)."),
    ("opls-des", "dlfield", "OPLS_DES", "OPLS-DES deep eutectic solvents (DL_FIELD)", "2018",
     ["B. Doherty, O. Acevedo, J. Phys. Chem. B 122, 9982 (2018)"], C, "same forms as OPLS 2005; not yet compared with DL_FIELD."),
    ("amber", "dlfield", "AMBER", "AMBER (DL_FIELD)", "Cornell et al. 1995 family",
     ["W. D. Cornell et al., J. Am. Chem. Soc. 117, 5179 (1995)"], V,
     "checked term by term against DL_FIELD's DL_POLY FIELD file (bench/ff/validate_family.py): BPTI (6PTI) and a methyl glucoside identical."),
    ("gaff-amber16", "dlfield", "AMBER16_gaff", "GAFF (AmberTools 16, DL_FIELD)", "GAFF 1.8x, 2016",
     [GAFF_REF], V, "matches DL_FIELD for the ester and propylbenzene from mol2 input."),
    ("gaff-amber25", "dlfield", "AMBER25_gaff", "GAFF (AmberTools 25, DL_FIELD)", "2025",
     [GAFF_REF], V,
     "matches DL_FIELD for the ester and propylbenzene; one carbonyl improper (X-X-c-o) is listed by DL_FIELD in an order that "
     "moves the key's end atom — CAPS follows the key."),
    ("charmm", "dlfield", "CHARMM", "CHARMM (DL_FIELD)", "", ["B. R. Brooks et al., J. Comput. Chem. 30, 1545 (2009)"], V,
     "checked term by term against DL_FIELD's DL_POLY FIELD file (bench/ff/validate_family.py): 48 of 48 buildable templates identical; CAPS adds impropers only at centres DL_FIELD's templates leave out."),
    ("charmm19", "dlfield", "CHARMM19", "CHARMM19 united atom (DL_FIELD)", "", ["E. Neria, S. Fischer, M. Karplus, J. Chem. Phys. 105, 1902 (1996)"], V,
     "checked term by term against DL_FIELD's DL_POLY FIELD file (bench/ff/validate_family.py): butane, octane, ethanol and DL_FIELD's C18 alkane example identical."),
    ("charmm22-prot", "dlfield", "CHARMM22_prot", "CHARMM22 proteins (DL_FIELD)", "1998", ["A. D. MacKerell Jr. et al., J. Phys. Chem. B 102, 3586 (1998)"], V,
     "checked term by term against DL_FIELD's DL_POLY FIELD file (bench/ff/validate_family.py): 35 of 35 buildable templates and the whole SOD1 protein (4412 bonds, 7986 angles, 8968 torsion terms, 3556 Urey-Bradley, 11472 1-4 pairs) identical."),
    ("charmm36-carb", "dlfield", "CHARMM36_carb", "CHARMM36 carbohydrates (DL_FIELD)", "", ["O. Guvench et al., J. Chem. Theory Comput. 7, 3162 (2011)"], C,
     "same forms as CHARMM36 (validated); DL_FIELD could not build its own carbohydrate examples."),
    ("cgenff", "dlfield", "CHARMM36_cgenff", "CGenFF (CHARMM36, DL_FIELD)", "", ["K. Vanommeslaeghe et al., J. Comput. Chem. 31, 671 (2010)"], V,
     "checked term by term against DL_FIELD's DL_POLY FIELD file (bench/ff/validate_family.py): 60 of 60 buildable templates identical, including Urey-Bradley, separate 1-4 LJ and NBFIX."),
    ("charmm36-lipid", "dlfield", "CHARMM36_lipid", "CHARMM36 lipids (DL_FIELD)", "2010", ["J. B. Klauda et al., J. Phys. Chem. B 114, 7830 (2010)"], C,
     "same forms as CHARMM36 (validated); the library holds no buildable lipid templates."),
    ("charmm36-nucl", "dlfield", "CHARMM36_nucl", "CHARMM36 nucleic acids (DL_FIELD)", "", ["K. Hart et al., J. Chem. Theory Comput. 8, 348 (2012)"], V,
     "checked term by term against DL_FIELD's DL_POLY FIELD file (bench/ff/validate_family.py): DL_FIELD's nucleic-acid example identical (DL_FIELD's nucleic templates carry no impropers; CAPS adds them from the rules)."),
    ("charmm36-prot", "dlfield", "CHARMM36_prot", "CHARMM36 proteins (DL_FIELD)", "2012", ["R. B. Best et al., J. Chem. Theory Comput. 8, 3257 (2012)"], V,
     "checked term by term against DL_FIELD's DL_POLY FIELD file (bench/ff/validate_family.py): BPTI (6PTI) identical: 896 bonds, 1608 angles, 1843 torsion terms, 700 Urey-Bradley, 2322 1-4 pairs."),
    ("dreiding", "dlfield", "DREIDING", "DREIDING (DL_FIELD)", "1990", ["S. L. Mayo, B. D. Olafson, W. A. Goddard III, J. Phys. Chem. 94, 8897 (1990)"], V,
     "checked term by term against DL_FIELD's DL_POLY FIELD file (bench/ff/validate_family.py): 29 of 29 buildable templates identical (planar inversions, linear angles, conjugated and ammonium torsions)."),
    ("gromos-54a7", "dlfield", "GROMOS_G54A7", "GROMOS 54A7 (DL_FIELD)", "2011", ["N. Schmid et al., Eur. Biophys. J. 40, 843 (2011)"], V,
     "checked term by term against DL_FIELD's DL_POLY FIELD file (bench/ff/validate_family.py): octane identical; bonds kept as the GROMOS quartic form (DL_FIELD writes a harmonic approximation)."),
    ("trappe-ua", "dlfield", "TRAPPE_UA", "TraPPE-UA (DL_FIELD)", "1998", ["M. G. Martin, J. I. Siepmann, J. Phys. Chem. B 102, 2569 (1998)"], V,
     "checked term by term against DL_FIELD's DL_POLY FIELD file (bench/ff/validate_family.py): 23 of 25 identical; the other two are DL_FIELD's C-C-O-H torsion, which takes terms from the O-C-C-O entry - CAPS keeps the published TraPPE alcohol torsion."),
    ("trappe-eh", "dlfield", "TRAPPE_EH", "TraPPE-EH (DL_FIELD)", "2007", ["N. Rai, J. I. Siepmann, J. Phys. Chem. B 111, 10790 (2007)"], V,
     "checked term by term against DL_FIELD's DL_POLY FIELD file (bench/ff/validate_family.py): 9 of 9 identical."),
    ("misc", "dlfield", "MISC_FF", "Miscellaneous (DL_FIELD)", "", [DLF_NOTE], C,
     "converted; DL_FIELD cannot write this scheme for LAMMPS and has no buildable examples."),
] + [
    (f"inorganic-{n.split('_', 1)[1].lower().replace('_', '-')}-dlfield", "dlfield", n, f"Inorganic: {n.split('_', 1)[1].replace('_', ' ')} (DL_FIELD)", "",
     [DLF_NOTE], INORG.get(n, (C, "Buckingham / Morse pairs and core-shell springs converted; no buildable DL_FIELD example."))[0],
     INORG.get(n, (C, "Buckingham / Morse pairs and core-shell springs converted; no buildable DL_FIELD example."))[1])
    for n in ["INORGANIC_binary_halides", "INORGANIC_binary_misc", "INORGANIC_binary_oxides", "INORGANIC_clay", "INORGANIC_glass",
              "INORGANIC_ternary_oxides", "INORGANIC_zeolite", "INORGANIC_zeolite_Hill_Sauer"]
] + [
    # ---- moltemplate ----
    ("compass-published-moltemplate", "moltemplate", "compass_published", "COMPASS (published subset, full class II)", "1998",
     ["H. Sun, J. Phys. Chem. B 102, 7338 (1998)"], V,
     "all cross terms; energies / forces equal LAMMPS CLASS2 (≤7e-8 kcal/mol/Å) for an ester-ether, propylbenzene and a siloxane. "
     "Every cross-term reference length / angle lies on its own bond / angle in CAPS (moltemplate's canonical sort misplaces "
     "up to 55 %, and it doubles bond increments with 'Data Bond List')."),
    ("gaff-moltemplate", "moltemplate", "gaff", "GAFF (moltemplate)", "GAFF 1.x", [GAFF_REF], C,
     "same forms as GAFF2 (validated). moltemplate uses lj/charmm/coul/long; CAPS evaluates plain LJ."),
    ("gaff2-moltemplate", "moltemplate", "gaff2", "GAFF2 (moltemplate)", "GAFF 2", [GAFF_REF + " (GAFF2: AmberTools)"], V,
     "methyl acetate: every coefficient identical, forces to 8.8e-8 kcal/mol/Å against LAMMPS."),
    ("oplsaa2024-moltemplate", "moltemplate", "oplsaa2024", "OPLS-AA (2024 parameter file)", "2024",
     ["W. L. Jorgensen et al., J. Phys. Chem. B (2024), doi:10.1021/acs.jpcb.3c06602"], V,
     "methyl benzoate: forces to 4.7e-8 kcal/mol/Å against LAMMPS."),
    ("oplsaa2008-moltemplate", "moltemplate", "oplsaa2008", "OPLS-AA (BOSS 4.8, 2008)", "2008",
     ["W. L. Jorgensen, D. S. Maxwell, J. Tirado-Rives, J. Am. Chem. Soc. 118, 11225 (1996)"], C, "same forms as the 2024 file."),
    ("loplsaa2024-moltemplate", "moltemplate", "loplsaa2024", "L-OPLS overlay (on OPLS-AA 2024)", "2012 / 2015",
     ["S. W. I. Siu, K. Pluhackova, R. A. Böckmann, J. Chem. Theory Comput. 8, 1459 (2012)",
      "K. Pluhackova et al., J. Phys. Chem. B 119, 15287 (2015)"], C, "overlay: merge on top of oplsaa2024 (caps ff apply --overlay)."),
    ("loplsaa2008-moltemplate", "moltemplate", "loplsaa2008", "L-OPLS overlay (on OPLS-AA 2008)", "2012 / 2015",
     ["S. W. I. Siu, K. Pluhackova, R. A. Böckmann, J. Chem. Theory Comput. 8, 1459 (2012)"], C, "overlay on oplsaa2008."),
    ("dreiding-moltemplate", "moltemplate", "dreiding", "DREIDING (moltemplate)", "1990",
     ["S. L. Mayo, B. D. Olafson, W. A. Goddard III, J. Phys. Chem. 94, 8897 (1990)"], C,
     "not yet compared; its hydrogen-bond term is not evaluated by CAPS."),
    ("trappe1998-moltemplate", "moltemplate", "trappe1998", "TraPPE-UA alkanes", "1998",
     ["M. G. Martin, J. I. Siepmann, J. Phys. Chem. B 102, 2569 (1998)"], C, "united atoms: element guesses from mass are not meaningful."),
    ("sdk-moltemplate", "moltemplate", "sdk", "SDK coarse-grained", "2007 / 2010",
     ["W. Shinoda, R. DeVane, M. L. Klein, Mol. Simul. 33, 27 (2007)"], C, "lj/sdk 9-6 / 12-4 forms are not evaluated by CAPS yet."),
    ("martini-moltemplate", "moltemplate", "martini", "MARTINI 2.0", "2007",
     ["S. J. Marrink et al., J. Phys. Chem. B 111, 7812 (2007)"], C, "lj/gromacs shifted forms are not evaluated by CAPS yet."),
    ("drymartini-moltemplate", "moltemplate", "drymartini", "Dry MARTINI", "2015",
     ["C. Arnarez et al., J. Chem. Theory Comput. 11, 260 (2015)"], C, "as MARTINI."),
    ("cooke-deserno-moltemplate", "moltemplate", "cooke_deserno_lipid", "Cooke–Deserno lipid model", "2005",
     ["I. R. Cooke, M. Deserno, J. Chem. Phys. 123, 224710 (2005)"], C, "cosine-squared attraction not evaluated yet."),
    ("graphene-moltemplate", "moltemplate", "graphene", "Graphene (LJ carbon)", "", ["moltemplate example"], C, ""),
    ("spce-moltemplate", "moltemplate", "spce", "SPC/E water", "1987",
     ["H. J. C. Berendsen, J. R. Grigera, T. P. Straatsma, J. Phys. Chem. 91, 6269 (1987)"], C, "rigid model: needs SHAKE / RATTLE."),
    ("tip3p-1983-moltemplate", "moltemplate", "tip3p_1983", "TIP3P water (1983)", "1983",
     ["W. L. Jorgensen et al., J. Chem. Phys. 79, 926 (1983)"], C, "rigid model: needs SHAKE / RATTLE."),
    ("tip3p-2004-moltemplate", "moltemplate", "tip3p_2004", "TIP3P water (Ewald, 2004)", "2004",
     ["D. J. Price, C. L. Brooks III, J. Chem. Phys. 121, 10096 (2004)"], C, "rigid model: needs SHAKE / RATTLE."),
    ("mw-moltemplate", "moltemplate", "watmw", "mW water", "2009",
     ["V. Molinero, E. B. Moore, J. Phys. Chem. B 113, 4008 (2009)"], C, "Stillinger–Weber three-body form not evaluated by CAPS."),
]
# moltemplate files that are molecules built from another force field (recorded, not converted)
TEMPLATES = ["graphite", "spc_oplsaa", "spc_oplsaa2008", "spc_oplsaa2024", "spce_oplsaa", "spce_oplsaa2024", "spce_ice_rect8",
             "spce_ice_rect16", "spce_ice_rect32", "tip3p_1983_oplsaa", "tip3p_1983_oplsaa2008", "tip3p_1983_oplsaa2024",
             "tip3p_2004_oplsaa", "tip3p_2004_oplsaa2024", "tip5p_oplsaa", "tip5p_oplsaa2008", "tip5p_oplsaa2024"]
# variants that are the same parameters under another name
ALIASES = {"oplsaa": "oplsaa2024", "loplsaa": "loplsaa2024", "graphene_rectangular": "graphene", "spce_amber": "spce",
           "spce_dreiding": "spce", "spce_hybrid": "spce", "spce_more_comments": "spce", "tip3p_1983_charmm": "tip3p_1983",
           "tip3p_1983_charmm_hybrid": "tip3p_1983", "tip3p_1983_hybrid": "tip3p_1983", "tip3p_2004_hybrid": "tip3p_2004"}

# automatic typing rules (data/typing) and how they were checked (bench/ff/validate_typing.py, compare_msi_types.py)
TYPING = {
    "pcff": ("pcff.typing.json",
                     "92 DL_FIELD PCFF molecule templates typed from connectivity alone: 1387 / 1400 atoms (99.1 %) as DL_FIELD's "
                     "templates (the rest: DL_FIELD's templates give furan, oxazole and indole ring carbons cp but pyrrole's c5, and "
                     "neutral histidine carbons ci)."),
    "cvff": ("cvff.typing.json",
                     "23 DL_FIELD CVFF templates: 272 / 278 atoms (97.8 %; the rest are TIP3P / SPC water, chosen with the "
                     "cvff-tip3p / -spc overlays); Materials Studio's own typing of the msi2lmp examples (crambin, nylon, "
                     "aromatics): 811 / 823 atoms (98.5 %)."),
    "gaff-amber16": ("gaff-amber16.typing.json",
                             "Rules follow antechamber's ATOMTYPE_GFF.DEF in order, with its ring classes and cc/cd-type conjugation "
                             "pairs. 32 DL_FIELD GAFF templates (water models aside): 365 / 373 atoms (97.9 %); the rest are template "
                             "deviations from antechamber (DMSO methyl H hc, azobenzene N n2)."),
    "gaff-amber25": ("gaff-amber25.typing.json",
                             "GAFF2 rules after antechamber's ATOMTYPE_GFF2.DEF (c5 / c6, n7-n9, ns / nt, nu / nv ...). 27 DL_FIELD "
                             "templates: 304 / 316 atoms (96.2 %; the rest: templates type H on C-O / C-S as hc where antechamber "
                             "gives h1). A 6440-atom polystyrene melt: every atom as CAPS's GAFF builder typed it."),
    "opls2005": ("opls2005.typing.json",
                         "99 DL_FIELD OPLS2005 templates (water models aside): 1368 / 1398 atoms (97.9 %). Residue-specific keys "
                         "(amino-acid side-chain CB, sugars, nucleic bases) are not assigned automatically."),
    "cgenff": ("cgenff.typing.json",
                       "Rules from CGenFF's own type definitions. The 494 model compounds of top_all36_cgenff.rtf (CGenFF 3.0.1), "
                       "typed from their bonds alone: 9024 / 9174 atoms (98.4 %) as CGenFF's developers typed them, 417 residues "
                       "fully right."),
}
os.makedirs(OUT, exist_ok=True)
catalogue = []
for fid, src, fname, name, version, refs, status, notes in ENTRIES:
    path = os.path.join(MT, fname + ".lt") if src == "moltemplate" else os.path.join(DLF, fname + ".par")
    entry = {"id": fid, "name": name, "version": version, "references": refs, "origin": src,
             "source_file": os.path.basename(path), "status": status, "notes": notes, "file": None}
    if status != P and os.path.exists(path):
        out = os.path.join(OUT, fid + ".json")
        cmd = [CAPS, "ff", "import-lt" if src == "moltemplate" else "import-dlf", path, "-o", out]
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode != 0:
            entry["status"], entry["notes"] = "failed", (r.stderr or r.stdout).strip()
        else:
            ff = json.load(open(out))
            ff["name"], ff["version"] = name, version
            ff["references"] = refs + ff.get("references", [])
            ff["source"] = f"{'moltemplate' if src == 'moltemplate' else 'DL_FIELD 4.13'}: {os.path.basename(path)}"
            ff["status"] = status
            ff["validation"] = notes
            if fid in TYPING:
                ff["typing"] = "../typing/" + TYPING[fid][0]
                entry["typing"] = {"rules": "typing/" + TYPING[fid][0], "evidence": TYPING[fid][1]}
            json.dump(ff, open(out, "w"), indent=1, ensure_ascii=False)
            entry["file"] = fid + ".json"
            entry["counts"] = {k: len(ff.get(k, [])) for k in ("atom_types", "pairs", "bonds", "angles", "dihedrals", "impropers")}
    elif not os.path.exists(path):
        entry["status"], entry["notes"] = "missing", f"{path} not found"
    catalogue.append(entry)
for t in TEMPLATES:
    catalogue.append({"id": t + "-moltemplate", "name": t, "origin": "moltemplate", "source_file": t + ".lt", "status": T,
                      "notes": "a molecule written with another force field's atom types; build it with that force field", "file": None})
for a, target in ALIASES.items():
    catalogue.append({"id": a + "-moltemplate", "name": a, "origin": "moltemplate", "source_file": a + ".lt", "status": "alias",
                      "notes": f"the same parameters as {target} (style or comment variant)", "file": None})
json.dump({"format": "caps-forcefield-catalogue", "format_version": 1, "forcefields": catalogue},
          open(os.path.join(OUT, "catalogue.json"), "w"), indent=1, ensure_ascii=False)

# README table
rows = ["| id | force field | version / year | from | status | automatic typing |", "|---|---|---|---|---|---|"]
for e in catalogue:
    rows.append(f"| {e['id']} | {e['name']} | {e.get('version', '')} | {e['origin']}: {e['source_file']} | {e['status']} | "
                f"{'yes' if e.get('typing') else ''} |")
open(os.path.join(OUT, "README.md"), "w").write(
    "# CAPS force-field library\n\nBuilt by `bench/ff/build_library.py` from the published parameter files. Status: **validated** — "
    "checked interaction by interaction and in energies / forces against their reference implementations and LAMMPS (evidence in each "
    "file's `validation` field); **converted** — parsed into CAPS rules, not yet compared; **pending** — the family's "
    "functional forms are not interpreted yet; **template** / **alias** — not a separate force field.\n\n" + "\n".join(rows) + "\n")
counts = {}
for e in catalogue:
    counts[e["status"]] = counts.get(e["status"], 0) + 1
print(f"{len(catalogue)} catalogue entries in {OUT}: " + ", ".join(f"{v} {k}" for k, v in sorted(counts.items())))

# CAPS's own names in the library: no mention of the conversion tool (bench/ff/library_names.py)
import library_names  # noqa: E402
library_names.run()
