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


check(caps.abi_version() >= 22, f"ABI {caps.abi_version()}")
samples = os.environ["CAPS_SAMPLES"]
root = os.path.dirname(os.path.abspath(samples))
melt = caps.open(os.path.join(samples, "ps_melt.data"))
s = melt.summary()
check(s["atoms"] == 1300 and s["molecules"] == 10, f"open: {s['atoms']} atoms, {s['molecules']} molecules")
check(melt.tacticity()["centres"] > 10, "tacticity")
frc = os.path.join(tempfile.gettempdir(), "caps_test_override.frcmod")
open(frc, "w").write("override\nBOND\nc3-c3  222.2   1.5555\n\n")
melt_ff = caps.open(os.path.join(samples, "ps_melt.data"))
melt_ff.field.assign("gaff")
rep_ff = melt_ff.field.import_params(frc)
check(rep_ff["complete"] and rep_ff["imported"] > 100, f"frcmod import: {rep_ff['imported']} terms from the imported rule")
br = caps.polymer("[*]C/C=C\\C[*]", dp=6)
before = br.field.assign(os.path.join(root, "data", "forcefields", "oplsaa2024-moltemplate.json"))["complete"]
filled = br.field.import_params(os.path.join(root, "data", "forcefields", "opls2005.json"), fill_gaps=True)
check(not before and filled["complete"] and filled["filled_terms"] == ["dihedral filled: CM-CT-CT-CM × 5"],
      f"OPLS-AA 2024 on polybutadiene: gaps filled from OPLS 2005 only where missing: {filled['filled_terms']}")
xyz = caps.import_file(os.path.join(samples, "ps_melt.xyz"))
check(xyz.summary()["bonds"] == 1370 and xyz.summary()["molecules"] == 10, "import_file: bonds perceived, molecules split")
check(caps.import_preview(os.path.join(samples, "ps_melt.xyz"), bonds="none")["bonds"] == 0, "import_preview: bonds none")
pep = caps.build.peptide("AEAAAKA", structure="HHHHHHH", cleanup=False)
check(pep.atoms == 90, f"peptide: {pep.atoms} atoms")
check(pep.field.assign("uff")["complete"], "UFF assigned")
check(pep.relax(ftol=3.0) in (0, 1), "relax")
rep = pep.md(steps=50, dt=2.0, thermo_every=10, frame_every=50, constraints="h-bonds")
check("SHAKE/RATTLE" in rep, "md with bonds to hydrogen constrained")
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
    check([st["engine"] for st in cell.provenance()["steps"]][:3] == ["recipe.run", "field.assign", "grow.trials"], "run: provenance names the recipe")
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
melt2 = caps.open(os.path.join(samples, "ps_melt.data"))
check(melt2.query('smarts "c1ccccc1" and chain 1', op="preview") == 48 and melt2.query("element H") == 660, "query grammar")
# polymer statistics (row 18): the models, and a polydisperse terminal-model copolymer grown from them
cl = caps.chain_lengths(nn=40, pdi=1.1, count=20, seed=2026)
check(len(cl["lengths"]) == 20 and abs(cl["target"]["mn"] - 4166) < 1 and cl["sample"]["pdi"] > 1, "chain_lengths")
co = caps.copolymer_model(0.52, 0.46, 0.5)
check(abs(co["F1"] - 0.510) < 5e-4 and abs(co["azeotrope"] - 0.529) < 5e-4, "copolymer_model: Mayo–Lewis")
st = caps.stereo(pm=0.7)
fit = caps.stereo(measured=st["model"]["pentads"])["fit"]
check(abs(fit["bernoulli"]["pm"] - 0.7) < 1e-6 and abs(caps.stereo(dyads="mmrr")["chain"]["mm"] - 1 / 3) < 1e-9, "stereo: model, counts, fit")
bp = caps.blend_phase(100, 200, -0.02, 15, t=300)
check(abs(bp["tc"] - 433.9) < 0.1 and abs(bp["coexist"][0] - 0.0385) < 1e-4, "blend_phase")
check(caps.solvent_chi(18.6, [{"name": "toluene", "v": 106.3, "delta": 18.2}])[0]["predicted"] == "solvent", "solvent_chi")
check(caps.ewald_params(12, 1e-5, 1.2, edges=[45.3] * 3)["mesh"] == [40, 40, 40], "ewald_params")
sbr = caps.polymer(["*CC=CC*", "*CC(*)c1ccccc1"], sequence="terminal", r1=1.35, r2=0.58, weights=[0.75, 0.25], chains=3, density=0.3,
                   lengths={"distribution": "schulz-zimm", "nn": 12, "pdi": 1.2, "seed": 4})
grow = [st for st in sbr.provenance()["steps"] if st["engine"] == "grow.trials"][0]
check("sample Nn" in grow["params"]["DP"] and "terminal model" in grow["params"]["sequence"], f"polymer: SBR, terminal model, drawn lengths in provenance · {sbr.atoms} atoms")
# display and resolution (row 19): the heavy-atom PS 4-mer gets its 34 H; C20H42 keeps its mass at every resolution
frag = caps.open(os.path.join(samples, "ps_frag.pdb"))
check(frag.hydrogen_plan()["add"] == 34 and frag.add_hydrogens() == 34 and frag.atoms == 66, "hydrogen_plan / add_hydrogens")
eic = caps.build.smiles("CCCCCCCCCCCCCCCCCCCC")
r = eic.resolution()
check(r["all_atom"]["sites"] == 62 and r["united_atom"]["sites"] == 20 and r["coarse_grained"]["sites"] == 4
      and abs(r["coarse_grained"]["mass"] - 282.556) < 1e-3, "resolution: 62 / 20 / 4 sites, mass conserved")
check(eic.convert("coarse-grained").atoms == 4 and eic.atoms == 62, "convert: a new document, the original kept")
# branched molecules: a 4-arm polybutadiene star (star BR, as coupled on silicon) is one molecule of 4 × 8 units
star = caps.polymer("*CC=CC*", dp=8, chains=2, density=0.3, architecture="star", arms=4, seed=3)
arch = [st for st in star.provenance()["steps"] if st["engine"] == "grow.trials"][0]["params"]["architecture"]
check(star.summary()["molecules"] == 2 and star.atoms == 2 * (4 * 8 * 10 + 4 - 2) and arch.startswith("star, 4 arms"), f"polymer: star BR, {star.atoms} atoms")
film = caps.polymer("*CC=CC*", dp=10, chains=4, density=0.6, region={"shape": "slab", "thickness": 20, "vacuum": 30}, seed=2)
check(abs(film.summary()["cell_c"] - 50.0) < 1e-6 and film.atoms == 4 * 102, f"polymer: a slab with vacuum · {film.atoms} atoms, c {film.summary()['cell_c']:.1f} Å")
rb = caps.polymer("*CC*", dp=20, chains=3, density=0.4, method="rosenbluth", seed=4)
st = [x for x in rb.provenance()["steps"] if x["engine"] == "grow.rosenbluth"]
check(len(st) == 1 and "ln W" in st[0]["params"]["method"], "polymer: Rosenbluth growth recorded in the provenance")
pair = caps.polymer("*CC*", dp=6, chains=2, density=0.3, seed=3)
pair.relax(ftol=1.0, restraints=[(0, 30, 4.5, 50)])
check("restraint 1–31: target 4.500 Å" in pair.report, "relax: a distance restraint is applied and reported")
pair.relax(ftol=1.0, restraints=[])
check("restraint" not in pair.report, "relax: restraints=[] clears them")
# χ by MD runs (a very short run: the call and the result's shape, not the number)
chi = caps.chi_by_md("*CC*", polymer_b="*CC*", dp=4, chains=2, eq_ps=0.2, prod_ps=0.2)
check(len(chi["cells"]) == 3 and "chi" in chi and chi["phi_a"] > 0.3, f"chi_by_md: χ {chi['chi']:.2f} ± {chi['chi_error']:.2f} (a smoke test)")
# GROMACS export: three files and the non-bonded settings (PME for a periodic cell)
with tempfile.TemporaryDirectory() as tmp:
    stem = os.path.join(tmp, "melt")
    mdp = caps.open(os.path.join(samples, "ps_melt.data")).save_gromacs(stem)
    have = all(os.path.exists(stem + x) for x in (".top", ".gro", ".mdp"))
    check(have and "coulombtype              = PME" in mdp and "[ nonbond_params ]" in open(stem + ".top").read(), "save_gromacs: top, gro, mdp")
# Export center: refused without a complete force field; with one, LAMMPS and GROMACS files and the run protocol
with tempfile.TemporaryDirectory() as tmp:
    melt = caps.open(os.path.join(samples, "ps_melt.data"))
    try:
        melt.export_engines(tmp)
        refused = False
    except caps.CapsError as e:
        refused = "force field" in str(e)
    melt.field.assign("gaff2")
    pkg = melt.export_engines(tmp, run="npt", steps=1000)
    names = {f["name"] for f in pkg["files"]}
    deck = open(os.path.join(tmp, "system.in")).read()
    check(refused and {"system.data", "system.in", "system.top", "system.itp", "system.gro", "system.mdp", "system_em.mdp"} <= names
          and "pair_coeff" in deck and "fix             integrate all npt" in deck and "Pair Coeffs" not in open(os.path.join(tmp, "system.data")).read()
          and "pcoupl                   = C-rescale" in open(os.path.join(tmp, "system.mdp")).read() and pkg["checks"]["missing"] == 0,
          "export_engines: refused before a force field; LAMMPS data + in and GROMACS top/itp/gro/mdp with the NPT protocol")
# parameters by analogy (parmchk2's approach): GAFF2 lacks torsions through a cyclopropyl carbon (cx); they are taken
# from c3 and listed as estimated, never filled silently
cp = caps.build.smiles("C1CC1CCO", forcefield="uff")
rep = cp.field.assign("gaff2")
check(rep["complete"] and rep["estimated"] > 0 and all(" as " in x for x in rep["by_analogy"]) and not rep["missing"],
      f"field: {rep['estimated']:.0f} terms by analogy, e.g. {rep['by_analogy'][0] if rep['by_analogy'] else '-'}")
# united-atom force fields: TraPPE-UA folds the hydrogens on carbon into their carbons, then types CH3 / CH2 sites
pe = caps.polymer("*CC*", dp=6, chains=1, density=0.3, seed=3)
ua = pe.field.assign("trappe-ua")
types = {a["type"] for a in ua["atoms"]}
check(ua["complete"] and types <= {"CH2", "CH3"} and any("united-atom" in n for n in ua["notes"]),
      f"field: TraPPE-UA on polyethylene · {len(ua['atoms'])} sites, types {sorted(types)}")
# OPLS descendants: reline (choline chloride + 2 urea) in OPLS-DES, an imidazolium NTf2 in CL&P, polyethylene in L-OPLS
reline = caps.build.smiles("C[N+](C)(C)CCO.NC(N)=O.NC(N)=O.[Cl-]", forcefield="uff").field.assign("opls-des")
il = caps.build.smiles("CCCCn1cc[n+](C)c1.FC(F)(F)S(=O)(=O)[N-]S(=O)(=O)C(F)(F)F", forcefield="uff").field.assign("opls-clp")
lopls = caps.polymer("*CC*", dp=8, chains=1, density=0.3, seed=3).field.assign("loplsaa2024-moltemplate")
ltypes = {a["type"].split("_")[0] for a in lopls["atoms"]}
check(reline["complete"] and il["complete"] and lopls["complete"] and {"54L", "57L", "60LCH2"} <= ltypes,
      f"field: OPLS-DES reline, CL&P [C4mim][NTf2], L-OPLS polyethylene ({sorted(ltypes)})")
# OPLS 2020: the bond / angle supplement on OPLS 2005, typed with OPLS 2005's rules; the 2020 values replace 2005's
o20 = caps.build.smiles("CC(=O)OCc1ccccc1", forcefield="uff")
r20 = o20.field.assign("opls2020", charges="gasteiger")
check(r20["complete"] and "2020" in r20.get("forcefield", r20.get("name", "OPLS 2020")), f"field: OPLS 2020 on benzyl acetate · {r20.get('forcefield', r20.get('name'))}")
# inorganic and water: NaCl typed as NaCl's potential set (neighbour bonds dropped, the set's charges); MgO's shell
# model gets a shell on every core; SPC/E water its charges
nacl = caps.build.crystal("F m -3 m", 5.64, 5.64, 5.64, [("Na1", "Na", 0, 0, 0), ("Cl1", "Cl", 0.5, 0.5, 0.5)], supercell=(2, 2, 2))
rn = nacl.field.assign("inorganic-binary-halides", charges="forcefield")
mgo = caps.build.crystal("F m -3 m", 4.212, 4.212, 4.212, [("Mg1", "Mg", 0, 0, 0), ("O1", "O", 0.5, 0.5, 0.5)], supercell=(2, 2, 2))
rm = mgo.field.assign("inorganic-binary-oxides", charges="forcefield")
mtypes = sorted({a["type"] for a in rm["atoms"]})
w = caps.build.smiles("O", forcefield="uff").field.assign("spce-moltemplate", charges="forcefield")
check(rn["complete"] and {a["type"] for a in rn["atoms"]} == {"Na1", "Cl1"} and abs(rn["atoms"][0]["q"]) == 0.988
      and rm["complete"] and mtypes == ["Mg1c", "Mg1s", "O2c", "O2s"] and len(rm["atoms"]) == 128
      and w["complete"] and abs(w["atoms"][0]["q"] + 0.8476) < 1e-9,
      f"field: NaCl potential set, MgO shell model ({len(rm['atoms'])} sites: {mtypes}), SPC/E water")
# miscellaneous set: each molecule-specific set where its charges add up (HFA-134a's two carbons split by charge)
hfa = caps.build.smiles("FCC(F)(F)F", forcefield="uff").field.assign("misc", charges="forcefield")
ipn = caps.build.smiles("CCC(C)(C)C", forcefield="uff").field.assign("misc", charges="forcefield")
check(hfa["complete"] and ipn["complete"] and {a["type"] for a in hfa["atoms"]} == {"CTf3", "CTf1", "F", "HC"}
      and {a["type"] for a in ipn["atoms"]} == {"Cp1", "Cs1", "C1", "H1"}
      and abs(sum(a["q"] for a in hfa["atoms"])) < 1e-9 and abs(sum(a["q"] for a in ipn["atoms"])) < 1e-9,
      "field: misc (HFA-134a, 2,2-dimethylbutane)")
# mW: a water box becomes one Stillinger–Weber site per molecule (the document converted, undoable), which relaxes
wb = caps.build.solvate(None, shape="cube", edge=20, water_model="SPC/E", ion_mode=0)
nw = wb.summary()["molecules"]
rw = wb.field.assign("mw-moltemplate")
check(rw["complete"] and wb.atoms == nw and {a["type"] for a in rw["atoms"]} == {"MW"} and wb.relax(ftol=1.0) in (0, 1),
      f"field: mW water ({nw} molecules -> {wb.atoms} Stillinger-Weber sites), relaxed")
# coarse-grained: MARTINI's DPPC template (typed by bead name, charges kept); SDK maps all-atom DMPC onto its beads
tpl = caps.bead_templates("martini-moltemplate")
dppc = caps.build.beads("DPPC", forcefield="martini-moltemplate")
rd = dppc.field.assign("martini-moltemplate", charges="keep")
dmpc = caps.build.smiles("CCCCCCCCCCCCCC(=O)OCC(COP(=O)([O-])OCC[N+](C)(C)C)OC(=O)CCCCCCCCCCCCC", forcefield="")
rs = dmpc.field.assign("sdk-moltemplate", charges="forcefield")
check("DPPC" in tpl and dppc.atoms == 12 and rd["complete"] and {a["type"].split("_b")[0] for a in rd["atoms"]} == {"Q0", "Qa", "Na", "C1"}
      and rs["complete"] and dmpc.atoms == 13 and abs(sum(a["q"] for a in rs["atoms"])) < 1e-9,
      f"coarse-grained: MARTINI DPPC template ({len(tpl)} templates), SDK DMPC mapped ({dmpc.atoms} beads)")
# MARTINI overlays from the source's other parameter files: PEO (polymers.prm, torsions), sucrose (sugars.prm)
peo = caps.build.beads("PEO", forcefield="martini-polymers")
suc = caps.build.beads("SUCR", forcefield="martini-sugars")
rp, ru = peo.field.assign("martini-polymers", charges="keep"), suc.field.assign("martini-sugars", charges="keep")
check(rp["complete"] and ru["complete"] and peo.atoms == 37, f"MARTINI overlays: PEO {peo.atoms} beads, sucrose {suc.atoms} beads")
trp = caps.build.beads("TRP", forcefield="martini-aminoacids")
rt = trp.field.assign("martini-aminoacids", charges="keep")
check(rt["complete"] and trp.atoms == 5 and "ILE" not in caps.bead_templates("martini-aminoacids"),
      f"MARTINI amino acids: TRP {trp.atoms} beads with its impropers; ILE / LEU / PRO / VAL left out (no AC1 / AC2 parameters)")
# Martini 2.2 proteins: an all-atom helical peptide mapped as martinize does (DSSP, explicit topology)
ak = caps.build.peptide("AEAAAKEAAAKEAAAKA", structure="H" * 17)
ra = ak.field.assign("martini22-proteins", charges="keep")
check(ra["complete"] and ak.atoms == 26, f"Martini 2.2 protein: AK peptide -> {ak.atoms} beads with its topology")
# χ from pair contacts: the self-mixing control is 0 within its error; a hydrocarbon against water is far above ½
ctl = caps.chi_by_contacts("*CC*", "*CC*", samples=200000, pack_trials=1000)
wat = caps.chi_by_contacts("*CC*", "O", samples=200000, pack_trials=1000)
check(abs(ctl["chi"]) < 3 * ctl["chi_error"] + 0.05 and wat["chi"] > 2 and len(ctl["kinds"]) == 4,
      f"chi_by_contacts: control {ctl['chi']:.3f} ± {ctl['chi_error']:.3f}, water {wat['chi']:.2f}")
# a sulfur cure: H–S–S–H donors inserted into a natural-rubber cell, cured, then the network typed with PCFF
nr = caps.polymer("[*]C/C=C(C)\\C[*]", dp=8, chains=3, density=0.5, seed=2)
n0 = nr.atoms
nr.insert("SS", 6)
cure = nr.react("sulfur_allylic", relax=False, seed=3)
rn = nr.field.assign("pcff-frc")
check("sulfur_allylic" in caps.reaction_templates() and "reactions" in cure and nr.atoms < n0 + 24 and rn["complete"],
      f"react: sulfur cure {n0} -> {nr.atoms} atoms · {cure.splitlines()[0] if cure else ''} · PCFF complete {rn['complete']}")
print("all python checks passed")
