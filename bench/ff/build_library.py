#!/usr/bin/env python3
"""Build the curated CAPS force-field library (data/forcefields) from moltemplate's force-field files.

usage: build_library.py [--moltemplate DIR] [--out data/forcefields]

Every source file gets a catalogue entry, whether or not it converts: its name, version / year, primary
reference, origin file, what CAPS can do with it, and the evidence behind that status. Entries whose status is
"validated" were checked interaction by interaction (and energies / forces) against their reference implementations and LAMMPS
with bench/ff/compare_moltemplate.py or compare_class2.py.
"""
import json, os, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
CAPS = os.path.join(ROOT, "build", "cli", "caps")
arg = lambda k, d: sys.argv[sys.argv.index(k) + 1] if k in sys.argv else d
MT = arg("--moltemplate", os.path.expanduser("~/moltemplate/moltemplate/force_fields"))
OUT = arg("--out", os.path.join(ROOT, "data", "forcefields"))

V = "validated"
C = "converted"          # parsed into CAPS rules; not yet checked against the source tool
P = "pending"            # the family's functional forms are not interpreted yet
T = "template"           # a molecule built from another force field's types, not a force field

GAFF_REF = "J. Wang, R. M. Wolf, J. W. Caldwell, P. A. Kollman, D. A. Case, J. Comput. Chem. 25, 1157 (2004)"

# id, source, file, name, version / year, references, status, evidence / notes
ENTRIES = [
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

os.makedirs(OUT, exist_ok=True)
catalogue = []
for fid, src, fname, name, version, refs, status, notes in ENTRIES:
    path = os.path.join(MT, fname + ".lt")
    entry = {"id": fid, "name": name, "version": version, "references": refs, "origin": src,
             "source_file": os.path.basename(path), "status": status, "notes": notes, "file": None}
    if status != P and os.path.exists(path):
        out = os.path.join(OUT, fid + ".json")
        cmd = [CAPS, "ff", "import-lt", path, "-o", out]
        r = subprocess.run(cmd, capture_output=True, text=True)
        if r.returncode != 0:
            entry["status"], entry["notes"] = "failed", (r.stderr or r.stdout).strip()
        else:
            ff = json.load(open(out))
            ff["name"], ff["version"] = name, version
            ff["references"] = refs + ff.get("references", [])
            ff["source"] = f"moltemplate: {os.path.basename(path)}"
            ff["status"] = status
            ff["validation"] = notes
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

