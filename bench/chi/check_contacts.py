#!/usr/bin/env python3
"""χ by pair contacts (Fan, Olafson, Blanco & Hsu 1992; core chipair.hpp) against known behaviour.

Every polymer of data/reference/solvents.json against every solvent there (GAFF2, Gasteiger charges), scored with the
Solvent screen's rule (χ < 0.45 solvent, 0.45–0.55 borderline, above that a non-solvent) against the file's known
behaviour, beside the Hildebrand estimate from the same file; then blends with a known answer and the self-mixing
controls (χ must be 0). Nothing here is tuned to the answers: the force field and the settings are the Studio's.

usage: PYTHONPATH=data/python CAPS_LIB=build/capi/libcaps_c.dylib python3 bench/chi/check_contacts.py [--t 298.15]
"""
import json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(ROOT, "data", "python"))
os.environ.setdefault("CAPS_LIB", os.path.join(ROOT, "build", "capi", "libcaps_c.dylib" if sys.platform == "darwin" else "libcaps_c.so"))
import caps

T = float(sys.argv[sys.argv.index("--t") + 1]) if "--t" in sys.argv else 298.15
data = json.load(open(os.path.join(ROOT, "data", "reference", "solvents.json")))


def verdict(chi):
    return "solvent" if chi < 0.45 else "borderline" if chi <= 0.55 else "non-solvent"


def score(pred, known):
    good, bad = "good" in known, "non-solvent" in known
    if (pred == "solvent" and good) or (pred == "non-solvent" and bad) or (pred == "borderline" and "Θ" in known):
        return "agrees"
    if (pred == "solvent" and bad) or (pred == "non-solvent" and good):
        return "WRONG"
    return "—"


tally = {"contacts": [0, 0], "hildebrand": [0, 0]}
print(f"Polymer–solvent at {T:.2f} K (χ contacts, GAFF2 · χ Hildebrand + 0.34 · known)")
for p in data["polymers"]:
    for s in data["solvents"]:
        known = p["known"].get(s["name"], "")
        r = caps.chi_by_contacts(p["unit"], s["smiles"], t=T)
        c = r["fit_a"] + r["fit_b"] / T
        h = s["v"] * (s["delta"] - p["delta"]) ** 2 / (8.314462618 * T) + 0.34
        sc, sh = score(verdict(c), known), score(verdict(h), known)
        for k, v in (("contacts", sc), ("hildebrand", sh)):
            tally[k][0] += v == "agrees"
            tally[k][1] += v != "—" or True
        print(f"  {p['short']:3s} {s['name']:12s} contacts {c:7.2f} ± {r['chi_error']:.2f} {verdict(c):11s} {sc:6s} · "
              f"Hildebrand {h:5.2f} {verdict(h):11s} {sh:6s} · known: {known}")
print(f"agree with known behaviour: contacts {tally['contacts'][0]} of {tally['contacts'][1]}, "
      f"Hildebrand {tally['hildebrand'][0]} of {tally['hildebrand'][1]}")

print("\nBlends (repeat units capped with H)")
BLENDS = [("PS / PVME", "*CC(*)c1ccccc1", "*CC(*)OC", "miscible (LCST near 420 K): χ ≤ 0 expected"),
          ("PS / PI", "*CC(*)c1ccccc1", "*C/C=C(C)\\C*", "immiscible: χ ≫ χ_c"),
          ("NR / BR", "*C/C=C(C)\\C*", "*C/C=C\\C*", "near-miscible: small χ"),
          ("PS / PS (control)", "*CC(*)c1ccccc1", "*CC(*)c1ccccc1", "χ = 0"),
          ("NR / NR (control)", "*C/C=C(C)\\C*", "*C/C=C(C)\\C*", "χ = 0")]
for name, a, b, expect in BLENDS:
    r = caps.chi_by_contacts(a, b, t=T)
    print(f"  {name:18s} χ({T:.0f} K) = {r['fit_a'] + r['fit_b'] / T:6.3f} ± {r['chi_error']:.3f}; χ(T) = {r['fit_a']:.3f} + {r['fit_b']:.1f}/T · expected {expect}")
