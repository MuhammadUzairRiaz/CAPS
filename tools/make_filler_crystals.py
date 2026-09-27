"""Builds the filler crystals of data/crystals from their space groups, cells and asymmetric units (the values below,
from the references named) with CAPS's own crystal builder, and checks each against its known density.

usage: python3 tools/make_filler_crystals.py [--write]   (without --write: check only)
"""
import itertools, json, math, os, re, subprocess, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CAPS = os.path.join(ROOT, "build", "cli", "caps")
OUT = os.path.join(ROOT, "data", "crystals")

# id, name, space-group key, cell (a b c α β γ), sites, known density (g/cm³), use, surface hkl, reference
CRYSTALS = [
    ("corundum", "Corundum α-Al₂O₃", "167:h", (4.7589, 4.7589, 12.9910, 90, 90, 120), "Al1 Al 0 0 0.35216; O1 O 0.30624 0 0.25", 3.987,
     "alumina filler, thermally conductive composites", (0, 0, 1), "Lewis, Schwarzenbach & Flack, Acta Cryst. A38, 733 (1982)"),
    ("hematite", "Hematite α-Fe₂O₃", "167:h", (5.0356, 5.0356, 13.7489, 90, 90, 120), "Fe1 Fe 0 0 0.35530; O1 O 0.3059 0 0.25", 5.26,
     "iron oxide pigment and filler", (0, 0, 1), "Blake, Hessevick, Zoltai & Finger, Am. Mineral. 51, 123 (1966)"),
    ("calcite", "Calcite CaCO₃", "167:h", (4.9896, 4.9896, 17.0610, 90, 90, 120), "Ca1 Ca 0 0 0; C1 C 0 0 0.25; O1 O 0.25706 0 0.25", 2.71,
     "calcium carbonate filler", (1, 0, 4), "Maslen, Streltsov, Streltsova & Ishizawa, Acta Cryst. B49, 636 (1993)"),
    ("periclase", "Periclase MgO", "225", (4.2112, 4.2112, 4.2112, 90, 90, 90), "Mg1 Mg 0 0 0; O1 O 0.5 0.5 0.5", 3.58,
     "magnesia filler", (1, 0, 0), "Hazen, Am. Mineral. 61, 266 (1976)"),
    ("lime", "Lime CaO", "225", (4.8105, 4.8105, 4.8105, 90, 90, 90), "Ca1 Ca 0 0 0; O1 O 0.5 0.5 0.5", 3.34,
     "activator in sulfur cure, desiccant filler", (1, 0, 0), "Swanson & Fuyat, NBS Circular 539 (1953)"),
    ("anatase", "Anatase TiO₂", "141:1", (3.7845, 3.7845, 9.5143, 90, 90, 90), "Ti1 Ti 0 0 0; O1 O 0 0 0.2081", 3.89,
     "TiO₂ pigment, photocatalytic filler", (1, 0, 1), "Horn, Schwerdtfeger & Meagher, Z. Kristallogr. 136, 273 (1972)"),
    ("brucite", "Brucite Mg(OH)₂", "164", (3.142, 3.142, 4.766, 90, 90, 120), "Mg1 Mg 0 0 0; O1 O 0.33333333 0.66666667 0.2203; H1 H 0.33333333 0.66666667 0.4213", 2.39,
     "magnesium hydroxide flame retardant", (0, 0, 1), "Zigan & Rothbauer, Neues Jahrb. Mineral. Monatsh. 1967, 137 (Mg, O); H at 0.958 Å from O along c"),
    ("fluorite", "Fluorite CaF₂", "225", (5.4626, 5.4626, 5.4626, 90, 90, 90), "Ca1 Ca 0 0 0; F1 F 0.25 0.25 0.25", 3.18,
     "fluoride filler", (1, 1, 1), "Swanson & Tatge, NBS Circular 539 (1953)"),
    ("ceria", "Ceria CeO₂", "225", (5.411, 5.411, 5.411, 90, 90, 90), "Ce1 Ce 0 0 0; O1 O 0.25 0.25 0.25", 7.22,
     "UV-absorbing and radical-scavenging filler", (1, 1, 1), "Swanson & Tatge, NBS Circular 539 (1953)"),
    ("h-bn", "Hexagonal boron nitride", "194", (2.504, 2.504, 6.656, 90, 90, 120), "B1 B 0.33333333 0.66666667 0.25; N1 N 0.33333333 0.66666667 0.75", 2.28,
     "thermally conductive, electrically insulating filler", (0, 0, 1), "Pease, Acta Cryst. 5, 356 (1952)"),
    ("molybdenite", "2H-MoS₂", "194", (3.160, 3.160, 12.294, 90, 90, 120), "Mo1 Mo 0.33333333 0.66666667 0.25; S1 S 0.33333333 0.66666667 0.621", 5.06,
     "solid lubricant filler", (0, 0, 1), "Dickinson & Pauling, J. Am. Chem. Soc. 45, 1466 (1923); Bronsema et al., Z. Anorg. Allg. Chem. 540, 15 (1986)"),
    ("silicon-carbide", "3C-SiC (β, zinc blende)", "216", (4.3596, 4.3596, 4.3596, 90, 90, 90), "Si1 Si 0 0 0; C1 C 0.25 0.25 0.25", 3.21,
     "SiC particle and whisker reinforcement", (1, 1, 1), "Taylor & Jones, Phys. Rev. 1960; a = 4.3596 Å"),
    ("silicon", "Silicon (diamond)", "227:1", (5.4307, 5.4307, 5.4307, 90, 90, 90), "Si1 Si 0 0 0", 2.33,
     "silicon particles, wafers", (1, 1, 1), "Okada & Tokumaru, J. Appl. Phys. 56, 314 (1984)"),
    ("sphalerite", "Sphalerite ZnS", "216", (5.4093, 5.4093, 5.4093, 90, 90, 90), "Zn1 Zn 0 0 0; S1 S 0.25 0.25 0.25", 4.10,
     "ZnS pigment", (1, 1, 0), "Swanson & Fuyat, NBS Circular 539 (1953)"),
    ("barium-titanate", "Barium titanate BaTiO₃ (cubic, above 120 °C)", "221", (4.0094, 4.0094, 4.0094, 90, 90, 90), "Ba1 Ba 0 0 0; Ti1 Ti 0.5 0.5 0.5; O1 O 0.5 0.5 0", 6.01,
     "high-permittivity dielectric filler", (1, 0, 0), "Kwei, Lawson, Billinge & Cheong, J. Phys. Chem. 97, 2368 (1993)"),
    ("strontium-titanate", "Strontium titanate SrTiO₃", "221", (3.905, 3.905, 3.905, 90, 90, 90), "Sr1 Sr 0 0 0; Ti1 Ti 0.5 0.5 0.5; O1 O 0.5 0.5 0", 5.12,
     "dielectric filler", (1, 0, 0), "Mitchell et al., Phys. Rev. B 61, 11434 (2000)"),
    ("magnetite", "Magnetite Fe₃O₄", "227:2", (8.396, 8.396, 8.396, 90, 90, 90), "Fe1 Fe 0.125 0.125 0.125; Fe2 Fe 0.5 0.5 0.5; O1 O 0.2549 0.2549 0.2549", 5.20,
     "magnetic filler", (1, 1, 1), "Fleet, Acta Cryst. B37, 917 (1981)"),
    ("spinel", "Spinel MgAl₂O₄", "227:2", (8.0831, 8.0831, 8.0831, 90, 90, 90), "Mg1 Mg 0.125 0.125 0.125; Al1 Al 0.5 0.5 0.5; O1 O 0.2624 0.2624 0.2624", 3.58,
     "thermally conductive oxide filler", (1, 1, 1), "Yamanaka & Takéuchi, Z. Kristallogr. 165, 65 (1983)"),
    ("aluminium", "Aluminium (fcc)", "225", (4.0495, 4.0495, 4.0495, 90, 90, 90), "Al1 Al 0 0 0", 2.70,
     "metal flake and powder filler", (1, 1, 1), "Swanson & Tatge, NBS Circular 539 (1953)"),
    ("silver", "Silver (fcc)", "225", (4.0862, 4.0862, 4.0862, 90, 90, 90), "Ag1 Ag 0 0 0", 10.50,
     "conductive and antimicrobial filler", (1, 1, 1), "Swanson & Tatge, NBS Circular 539 (1953)"),
    ("nickel", "Nickel (fcc)", "225", (3.5240, 3.5240, 3.5240, 90, 90, 90), "Ni1 Ni 0 0 0", 8.91,
     "conductive filler", (1, 1, 1), "Swanson & Tatge, NBS Circular 539 (1953)"),
    ("titanium", "Titanium (hcp)", "194", (2.9508, 2.9508, 4.6855, 90, 90, 120), "Ti1 Ti 0.33333333 0.66666667 0.25", 4.51,
     "titanium powder, implant composites", (0, 0, 1), "Wood, Proc. Phys. Soc. 80, 783 (1962)"),
]


def main():
    write = "--write" in sys.argv
    tmp = os.path.join(ROOT, "build", "fillers")
    os.makedirs(tmp, exist_ok=True)
    bad = 0
    for cid, name, grp, cell, sites, rho, use, hkl, ref in CRYSTALS:
        out = os.path.join(tmp, cid + ".cif")
        r = subprocess.run([CAPS, "crystal", "--group", grp, "--cell", ",".join(str(x) for x in cell), "--sites", sites, "-o", out], capture_output=True, text=True)
        m = re.search(r"([\d.]+) g/cm³", r.stdout)
        got = float(m.group(1)) if m else -1
        ok = m is not None and abs(got - rho) / rho < 0.015
        bad += not ok
        print(f"{cid:20s} {grp:6s} {got:8.3f} g/cm³ (known {rho:.3f})  {'ok' if ok else 'CHECK'}  {r.stdout.strip().splitlines()[-1][:60] if r.stdout.strip() else r.stderr.strip()[-120:]}")
        if write and ok:
            text = open(out).read()
            head = f"# {name} — {grp} a b c α β γ = {' '.join(str(x) for x in cell)}; sites {sites}\n# {ref}; built by tools/make_filler_crystals.py\n"
            open(os.path.join(OUT, cid + ".cif"), "w").write(head + text)
    if write:
        # the catalogue keeps one line per crystal: new ones are appended in that style
        cat_path = os.path.join(OUT, "catalogue.json")
        text = open(cat_path).read()
        have = {c["id"] for c in json.loads(text)["crystals"]}
        lines = [json.dumps({"id": cid, "name": name, "file": cid + ".cif", "use": use, "hkl": list(hkl)}, ensure_ascii=False)
                 for cid, name, grp, cell, sites, rho, use, hkl, ref in CRYSTALS if cid not in have]
        if lines:
            end = text.rindex("}", 0, text.rindex("]"))
            text = text[:end + 1] + "".join(",\n    " + l for l in lines) + text[end + 1:]
            json.loads(text)
            open(cat_path, "w").write(text)
    print(f"\n{len(CRYSTALS) - bad} of {len(CRYSTALS)} at their known density")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
