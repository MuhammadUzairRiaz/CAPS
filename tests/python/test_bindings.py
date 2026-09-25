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
print("all python checks passed")
