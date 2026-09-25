"""The Python bindings (data/python/caps/core.py) against the built library: run by ctest as python_bindings."""
import os
import sys
import tempfile

import caps


def check(cond, what):
    if not cond:
        print("FAIL", what)
        sys.exit(1)
    print("ok  ", what)


check(caps.abi_version() >= 20, f"ABI {caps.abi_version()}")
samples = os.environ["CAPS_SAMPLES"]
melt = caps.open(os.path.join(samples, "ps_melt.data"))
s = melt.summary()
check(s["atoms"] == 1300 and s["molecules"] == 10, f"open: {s['atoms']} atoms, {s['molecules']} molecules")
check(melt.tacticity()["centres"] > 10, "tacticity")
xyz = caps.import_file(os.path.join(samples, "ps_melt.xyz"))
check(xyz.summary()["bonds"] == 1370 and xyz.summary()["molecules"] == 10, "import_file: bonds perceived, molecules split")
check(caps.import_preview(os.path.join(samples, "ps_melt.xyz"), bonds="none")["bonds"] == 0, "import_preview: bonds none")
pep = caps.build.peptide("AEAAAKA", structure="HHHHHHH", cleanup=False)
check(pep.atoms == 90, f"peptide: {pep.atoms} atoms")
check(pep.field.assign("uff")["complete"], "UFF assigned")
check(pep.relax(ftol=3.0) in (0, 1), "relax")
prov = pep.provenance()
engines = [st["engine"] for st in prov["steps"]]
check(engines[:3] == ["bio.peptide", "field.assign", "relax.lbfgs"], f"provenance: {engines}")
saved = os.path.join(tempfile.gettempdir(), "caps_prov_peptide.data")
pep.save(saved)
side = caps.provenance_file(saved)
reopened = caps.open(saved).provenance()
check(side["ok"] and [st["engine"] for st in reopened["steps"]] == engines, "provenance saved beside the file and read back")
md = caps.methods(prov)
check("minimised with L-BFGS" in md["text"] and len(md["refs"]) >= 3, "methods text from provenance")
check("@article{rappe1992," in caps.bibtex(prov) and "@article{engh1991," in caps.bibtex(prov), "bibtex of the cited methods")
with tempfile.TemporaryDirectory() as tmp:
    out = os.path.join(tmp, "pep.data")
    pep.save(out)
    back = caps.open(out)
    check(back.atoms == 90, "save and reopen")
    png = os.path.join(tmp, "pep.png")
    pep.render(png, 160, 120)
    check(os.path.getsize(png) > 1000, "render PNG")
rutile = caps.build.crystal("P 42/m n m", 4.5937, 4.5937, 2.9587, [("Ti1", "Ti", 0, 0, 0), ("O1", "O", 0.3048, 0.3048, 0)])
check(rutile.atoms == 6 and abs(rutile.summary()["density"] - 4.248) < 0.01, "crystal: rutile")
box = caps.build.solvate(edge=12, ion_mode=0)
check(box.interactions()["hbonds"] > 0, "solvate + interactions")
eth = caps.build.smiles("CC")
r = eth.edit(op="add_atom", to=0, element="O")
check(r["ok"] and eth.atoms == 9, "edit: place O")
eth.undo()
check(eth.atoms == 8, "undo")
check(len(caps.space_groups()) == 530, "space groups")
try:
    caps.open("/no/such/file.data")
    check(False, "a missing file raises")
except caps.CapsError:
    check(True, "a missing file raises CapsError")
# recipes, the notebook API (design/boards/CommandLine, JupyterNotebook)
with tempfile.TemporaryDirectory() as tmp:
    events = []
    cell = caps.run({"recipe": 1, "name": "pe", "build": {"polymer": {"smiles": "*CC*", "dp": 4, "chains": 2}}, "type": {"forcefield": "default"},
                     "grow": {"density": 0.5, "seed": 1}, "relax": {"fmax": 5}, "analyze": {"properties": ["density"]}, "export": ["pdb"]},
                    out_dir=tmp, progress=lambda e: events.append((e["name"], e["status"])))
    check(cell.atoms == 2 * 26 and os.path.exists(os.path.join(tmp, "pe.pdb")) and ("relax", "done") in events, "run: recipe with progress")
    check(abs(cell.properties[0]["value"] - 0.5) < 0.05, "run: analysed density")
    check([st["engine"] for st in cell.provenance()["steps"]][:2] == ["field.assign", "grow.trials"], "run: provenance")
    try:
        caps.run("build: {molecule: CCO}\nbogus: 1\n")
        check(False, "a wrong recipe raises")
    except caps.CapsError as e:
        check(e.exit_code == 2, "a wrong recipe raises CapsError with exit code 2")
    ps = caps.polymer("*CC(*)c1ccccc1", dp=6, tacticity="isotactic", seed=5)
    check(ps.atoms == 6 * 16 + 2, f"polymer: {ps.atoms} atoms")
    sc = ps.scene()
    check(sc["shown"] == ps.atoms and len(sc["xyz"]) == 3 * ps.atoms and len(sc["bonds"]) == 2 * (ps.atoms - 1 + 6), "scene")
    html = ps.view()._repr_html_()
    check("<canvas" in html and "atoms shown" in html and len(ps.view()._repr_png_()) > 1000, "view: HTML and PNG")
    bib = ps.provenance.citations(fmt="bibtex")
    check(len(bib) >= 2 and all(b.startswith("@") for b in bib), "provenance.citations(bibtex)")
    sw = caps.sweep.run("*CC*", tacticities=("atactic",), dps=(4,), seeds=(1, 2), chains=2, folder=os.path.join(tmp, "sw"))
    res = caps.sweep.result("atactic", folder=str(sw.folder))
    m, sd, n = res["density"]
    check(n == 2 and abs(m - 0.5) < 0.05, "sweep: result pooled over seeds")
    t = caps.table([res], ["density", "tg", "c_inf"])
    check("—" in repr(t) and "<table" in t._repr_html_() and t.to_csv().startswith("condition,density"), "table")
print("all python checks passed")
