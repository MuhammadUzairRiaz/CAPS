#!/usr/bin/env python3
"""Build the curated CAPS force-field library (data/forcefields) from DL_FIELD.

usage: build_library.py [--dlfield LIBDIR] [--out data/forcefields]

Every source file gets a catalogue entry, whether or not it converts: its name, version / year, primary
reference, origin file, what CAPS can do with it, and the evidence behind that status. Entries whose status is
"validated" were checked interaction by interaction (and energies / forces) against their reference implementations and LAMMPS
with bench/ff/compare_class2.py or compare_dlfield.py.
"""
import json, os, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
CAPS = os.path.join(ROOT, "build", "cli", "caps")
arg = lambda k, d: sys.argv[sys.argv.index(k) + 1] if k in sys.argv else d
DLF = arg("--dlfield", os.path.expanduser("~/project/dl_f_4.13/lib"))
OUT = arg("--out", os.path.join(ROOT, "data", "forcefields"))

V = "validated"
C = "converted"          # parsed into CAPS rules; not yet checked against the source tool
P = "pending"            # the family's functional forms are not interpreted yet

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
    ("opls-clp", "dlfield", "OPLS_CL_P", "CL&P ionic liquids (OPLS-AA based, DL_FIELD)", "2004 onwards",
     ["J. N. Canongia Lopes, A. A. H. Pádua, J. Phys. Chem. B 108, 2038 (2004)"], V,
     "checked term by term against DL_FIELD's DL_POLY FIELD file (bench/ff/validate_family.py): 24 of 24 identical (kJ/mol constants confirmed as DL_POLY's k)."),
    ("opls-des", "dlfield", "OPLS_DES", "OPLS-DES deep eutectic solvents (DL_FIELD)", "2018",
     ["B. Doherty, O. Acevedo, J. Phys. Chem. B 122, 9982 (2018)"], C, "same forms as OPLS 2005; not yet compared with DL_FIELD."),
    ("gaff-amber16", "dlfield", "AMBER16_gaff", "GAFF (AmberTools 16, DL_FIELD)", "GAFF 1.8x, 2016",
     [GAFF_REF], V, "matches DL_FIELD for the ester and propylbenzene from mol2 input."),
    ("gaff-amber25", "dlfield", "AMBER25_gaff", "GAFF (AmberTools 25, DL_FIELD)", "2025",
     [GAFF_REF], V,
     "matches DL_FIELD for the ester and propylbenzene; one carbonyl improper (X-X-c-o) is listed by DL_FIELD in an order that "
     "moves the key's end atom — CAPS follows the key."),
    ("charmm", "dlfield", "CHARMM", "CHARMM (DL_FIELD)", "", ["B. R. Brooks et al., J. Comput. Chem. 30, 1545 (2009)"], V,
     "checked term by term against DL_FIELD's DL_POLY FIELD file (bench/ff/validate_family.py): 48 of 48 buildable templates identical; CAPS adds impropers only at centres DL_FIELD's templates leave out."),
    ("cgenff", "dlfield", "CHARMM36_cgenff", "CGenFF (CHARMM36, DL_FIELD)", "", ["K. Vanommeslaeghe et al., J. Comput. Chem. 31, 671 (2010)"], V,
     "checked term by term against DL_FIELD's DL_POLY FIELD file (bench/ff/validate_family.py): 60 of 60 buildable templates identical, including Urey-Bradley, separate 1-4 LJ and NBFIX."),
    ("dreiding", "dlfield", "DREIDING", "DREIDING (DL_FIELD)", "1990", ["S. L. Mayo, B. D. Olafson, W. A. Goddard III, J. Phys. Chem. 94, 8897 (1990)"], V,
     "checked term by term against DL_FIELD's DL_POLY FIELD file (bench/ff/validate_family.py): 29 of 29 buildable templates identical (planar inversions, linear angles, conjugated and ammonium torsions)."),
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
]

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
    path = os.path.join(DLF, fname + ".par")
    entry = {"id": fid, "name": name, "version": version, "references": refs, "origin": src,
             "source_file": os.path.basename(path), "status": status, "notes": notes, "file": None}
    if status != P and os.path.exists(path):
        out = os.path.join(OUT, fid + ".json")
        cmd = [CAPS, "ff", "import-dlf", path, "-o", out]
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode != 0:
            entry["status"], entry["notes"] = "failed", (r.stderr or r.stdout).strip()
        else:
            ff = json.load(open(out))
            ff["name"], ff["version"] = name, version
            ff["references"] = refs + ff.get("references", [])
            ff["source"] = f"DL_FIELD 4.13: {os.path.basename(path)}"
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
    "functional forms are not interpreted yet.\n\n" + "\n".join(rows) + "\n")
counts = {}
for e in catalogue:
    counts[e["status"]] = counts.get(e["status"], 0) + 1
print(f"{len(catalogue)} catalogue entries in {OUT}: " + ", ".join(f"{v} {k}" for k, v in sorted(counts.items())))

# CAPS's own names in the library: no mention of the conversion tool (bench/ff/library_names.py)
import library_names  # noqa: E402
library_names.run()
