import math, random
from lib import *
from mols import *
from screen_app import shell, pagehead, card
from screen_studio2 import studio_shell, bar, spacer

NA = 6.02214076e23


# ---------------------------------------------------------------- Chain statistics
def chain_set(n_chains=24, n=80, bond=1.54, p=0.62, seed=40):
    return [random_chain((0.0, 0.0, 0.0), n + 1, seed + k, bond=bond, box=None, persistence=p) for k in range(n_chains)]


def chain_stats():
    bond = 1.54
    chains = chain_set()
    N = len(chains[0]) - 1
    # C_n from all internal sub-chains of n bonds
    cn = []
    for n in range(1, N + 1):
        acc, m = 0.0, 0
        for ch in chains:
            for i in range(0, N - n + 1):
                acc += dist(ch[i], ch[i + n]) ** 2; m += 1
        cn.append((n, acc / m / (n * bond * bond)))
    ree = [dist(ch[0], ch[-1]) for ch in chains]
    rg = []
    for ch in chains:
        c = mul(tuple(map(sum, zip(*ch))), 1 / len(ch))
        rg.append(math.sqrt(sum(dist(p, c) ** 2 for p in ch) / len(ch)))
    R2 = sum(r * r for r in ree) / len(ree)
    Rg2 = sum(r * r for r in rg) / len(rg)
    ratio = R2 / Rg2
    cN = cn[-1][1]
    cmax = math.ceil(max(v for _, v in cn))
    pl = plot(560, 250, [(cn, ACC, 2, None), ([(0, cN), (N, cN)], DIM, 1, "4 4")], (0, N), (0, cmax), [0, 20, 40, 60, 80], list(range(0, cmax + 1, max(1, cmax // 4))),
              "n (bonds in sub-chain)", "C<tspan>n</tspan> = ⟨R²ₙ⟩ / n l²", pad=(46, 12, 16, 32))
    pl = pl.replace("C<tspan>n</tspan>", "Cₙ")
    bins = 10
    lo, hi = 0, math.ceil(max(ree) / 10) * 10
    h = [0] * bins
    for r in ree:
        h[min(bins - 1, int((r - lo) / (hi - lo) * bins))] += 1
    hm = max(h)
    hw, hh = 560, 130
    hs = "".join(f'<rect x="{40 + i * (hw - 50) / bins + 1:.1f}" y="{hh - 24 - v / hm * (hh - 36):.1f}" width="{(hw - 50) / bins - 3:.1f}" height="{v / hm * (hh - 36):.1f}" rx="2" fill="{SEL}" opacity="0.8"></rect>' for i, v in enumerate(h))
    hist = (f'<svg width="{hw}" height="{hh}" viewBox="0 0 {hw} {hh}" role="img" aria-label="Histogram of end-to-end distance" style="display: block">{hs}<line x1="40" y1="{hh - 24}" x2="{hw - 10}" y2="{hh - 24}" stroke="{LINE}"></line>'
            + "".join(f'<text x="{40 + (v - lo) / (hi - lo) * (hw - 50):.1f}" y="{hh - 8}" text-anchor="middle" font-size="10" font-family="IBM Plex Mono" fill="{DIM}">{v}</text>' for v in range(lo, hi + 1, hi // 5))
            + f'<text x="{hw - 10}" y="12" text-anchor="end" font-size="10" fill="{DIM}">R<tspan font-size="8">ee</tspan> (Å) · {len(ree)} chains</text></svg>')
    sc = Scene("cs", 390, 300, yaw=0.5, pitch=0.35, persp=0.2, fog=0.35)
    for k, ch in enumerate(chains[:6]):
        sc.add_tube(ch, ACC if k == 0 else mix("#8B969E", BG0, 0.35), 0.5 if k == 0 else 0.35)
    sc.add_line(chains[0][0], chains[0][-1], SEL, 2, "6 4", 1)
    view = sc.svg()
    res = table(["Quantity", "Drawn chains", "Polystyrene"],
                [["Bonds per chain, n", str(N), "—"], ["⟨R²ₑₑ⟩<sup>½</sup>", f"{math.sqrt(R2):.1f} Å", "[result]"], ["⟨R²<sub>g</sub>⟩<sup>½</sup>", f"{math.sqrt(Rg2):.1f} Å", "[result]"],
                 ["⟨R²ₑₑ⟩ / ⟨R²<sub>g</sub>⟩", f"{ratio:.2f}", "[result]"], ["C<sub>n</sub> at n = N", f"{cN:.2f}", "[result]"], ["C<sub>∞</sub>", "extrapolate", "lit. ≈ 9.5–10"]],
                ["44%", "28%", "28%"], mono_cols=(1, 2), align_right=(1, 2), fs=12, rowh=25)
    left = (f'<div style="width: 420px; flex-shrink: 0; display: flex; flex-direction: column; gap: 12px">'
            + card("Sample", f'<div style="background: {BG0}; border-radius: 8px; overflow: hidden">{view}</div><span style="font-size: 11.5px; color: {DIM}">First 6 of {len(chains)} chains; highlighted chain with its end-to-end vector. Display: backbone trace — these are model chains of bond vectors, with no atoms to show.</span>', "", 14)
            + card("Values", res, "", 12) + '</div>')
    right = (f'<div style="flex-grow: 1; min-width: 0; display: flex; flex-direction: column; gap: 12px">'
             + card("Characteristic ratio", col(pl, f'<div style="font-size: 11.5px; color: {DIM}">Averaged over every sub-chain of each length. An ideal Gaussian chain has ⟨R²ₑₑ⟩/⟨R²<sub>g</sub>⟩ = 6; {len(chains)} chains is a small sample, so expect scatter around that.</div>', gap=6), chip("computed", OK, BG2), 14)
             + card("End-to-end distance", hist, "", 14) + '</div>')
    content = (pagehead("Chain statistics", "Size and stiffness of chains: R<sub>ee</sub>, R<sub>g</sub>, C<sub>n</sub>", row(select("Chains", "all · 24"), btn("Export CSV", ic="download"), gap=8))
               + f'<div style="flex-grow: 1; display: flex; gap: 14px; padding: 16px 20px; min-height: 0; overflow: hidden">{left}{right}</div>')
    return shell("Analyze", "CAPS Analyze — chain statistics", "Analyze › Chain statistics", content,
                 (f"<span>{len(chains)} drawn chains · {N} bonds · l = {bond} Å</span>", "<span>Middle column computed from the drawn chains · right column filled from a trajectory</span>"))


# ---------------------------------------------------------------- Density calculator
def density_calc():
    L = 40.0
    V_A3 = L ** 3
    V_cm3 = V_A3 * 1e-24
    comps = [("Water", "H₂O", 18.015, 0.99705, "25 °C"), ("Toluene", "C₇H₈", 92.14, 0.862, "25 °C"), ("Ethanol", "C₂H₆O", 46.07, 0.785, "25 °C")]
    rows = []
    for name, f, M, rho, T in comps:
        n_exact = rho * V_cm3 * NA / M
        n = round(n_exact)
        rho_got = n * M / NA / V_cm3
        rows.append([f"{name} · {f}", f"{M:.3f}", f"{rho}", f"{n_exact:,.1f}".replace(",", " "), f"{n:,}".replace(",", " "), f"{rho_got:.4f}"])
    tb = table(["Molecule", "M (g/mol)", "ρ target", "N exact", "N", "ρ achieved"], rows, ["26%", "14%", "14%", "16%", "12%", "18%"], mono_cols=(1, 2, 3, 4, 5), align_right=(1, 2, 3, 4, 5), fs=12, rowh=30, hl={0})
    # mixture: 20 PS chains DP 40 + toluene to 0.95 g/cm3 overall in same box
    M_ps = 104.15 * 40 + 2 * 1.008
    n_ps = 4
    m_ps = n_ps * M_ps / NA
    rho_mix = 0.95
    m_tot = rho_mix * V_cm3
    n_tol = round((m_tot - m_ps) / (92.14 / NA))
    assert n_tol > 0
    wt = m_ps / (m_ps + n_tol * 92.14 / NA)
    mix_ = col(kv("Polymer", f"{n_ps} × PS DP 40 · M = {M_ps:,.1f} g/mol".replace(",", " ")), kv("Overall target", f"{rho_mix} g/cm³"), kv("Toluene to add", f"{n_tol:,}".replace(",", " ")),
               kv("Polymer weight fraction", f"{wt * 100:.1f} %"), gap=7)
    # inverse: box edge for 1000 waters
    n_w = 1000
    L_w = (n_w * 18.015 / NA / 0.99705 * 1e24) ** (1 / 3)
    inv = col(row(field("Molecules", str(n_w)), select("Species", "Water"), field("ρ", "0.99705", "g/cm³"), gap=8), kv("Cubic box edge", f"{L_w:.3f} Å"), gap=8)
    sc = Scene("dc", 420, 300, yaw=0.6, pitch=0.35, persp=0.25, fog=0.45)
    rnd = random.Random(3)
    for _ in range(90):
        o = (rnd.uniform(1.5, L - 1.5), rnd.uniform(1.5, L - 1.5), rnd.uniform(1.5, L - 1.5))
        a, b = water(o, rnd)
        sc.add_atoms(a, b)
    sc.add_box((0, 0, 0), (L, 0, 0), (0, L, 0), (0, 0, L), MUTED, 1.2, None, 0.9)
    left = (f'<div style="width: 460px; flex-shrink: 0; display: flex; flex-direction: column; gap: 12px">'
            + card("Box", col(f'<div style="background: {BG0}; border-radius: 8px; overflow: hidden">{sc.svg()}</div>',
                              row(field("a", f"{L:.1f}", "Å"), field("b", f"{L:.1f}", "Å"), field("c", f"{L:.1f}", "Å"), gap=8),
                              kv("Volume", f"{V_A3:,.0f} Å³ = {V_cm3:.3e} cm³".replace(",", " ").replace("e-20", " × 10⁻²⁰")), gap=10), chip("preview: 90 of the molecules", DIM, BG2), 14)
            + '</div>')
    right = (f'<div style="flex-grow: 1; min-width: 0; display: flex; flex-direction: column; gap: 12px">'
             + card("How many molecules?", col(tb, f'<div style="font-family: {MONO}; font-size: 12px; color: {MUTED}">N = ρ · V · N<sub>A</sub> / M, rounded to an integer; ρ achieved uses the rounded N.</div>'), chip("computed", OK, BG2), 14)
             + '<div style="display: flex; gap: 12px">'
             + card("Polymer + solvent", mix_, "", 14, "flex: 1 1 0; min-width: 0")
             + card("Box for a given N", inv, "", 14, "flex: 1 1 0; min-width: 0") + '</div>'
             + f'<div style="font-size: 11.5px; color: {DIM}; line-height: 1.5">Target densities are experimental values at 25 °C used as inputs. Packing at the experimental density is a starting point; the equilibrated density comes from NPT and is reported separately.</div>'
             + '</div>')
    content = (pagehead("Pack · density calculator", "Molecule counts from box size and target density, and back", row(btn("Send to Pack", True, "pack"), gap=8))
               + f'<div style="flex-grow: 1; display: flex; gap: 14px; padding: 16px 20px; min-height: 0; overflow: hidden">{left}{right}</div>')
    return shell("Pack", "CAPS Pack — density calculator", "Pack › Density calculator", content, ("<span>N<sub>A</sub> = 6.02214076 × 10²³ mol⁻¹ (exact, SI 2019)</span>", "<span>40.0 Å cubic box</span>"))


# ---------------------------------------------------------------- Units & number formats
def units():
    kcal = 4.184
    eV_kJ = 96.485332
    kB = 1.380649e-23
    R_kJ = kB * NA / 1000
    T = 298.15
    conv = [("Energy", "1 kcal/mol", f"{kcal:.3f} kJ/mol", "exact (thermochemical calorie)"),
            ("Energy", "1 eV per particle", f"{eV_kJ:.3f} kJ/mol · {eV_kJ / kcal:.4f} kcal/mol", "e · N<sub>A</sub>"),
            ("Thermal energy", "k<sub>B</sub>T at 298.15 K", f"{R_kJ * T:.4f} kJ/mol · {R_kJ * T / kcal:.4f} kcal/mol", "k<sub>B</sub> exact"),
            ("Pressure", "1 atm", "1.01325 bar", "exact"), ("Pressure", "1 GPa", f"{1e4:,.0f} bar".replace(",", " "), "exact"),
            ("Density", "1 amu/Å³", f"{1.66053907:.5f} g/cm³", "1 u = 1.66053907 × 10⁻²⁴ g"),
            ("Diffusion", "1 Å²/ps", "1 × 10⁻⁸ m²/s", "exact"), ("Dipole", "1 e·Å", f"{1.602176634e-29 / 3.33564e-30:.4f} D", "1 D = 3.33564 × 10⁻³⁰ C·m")]
    ct = table(["Quantity", "From", "Equals", "Basis"], [[a, b, c, d] for a, b, c, d in conv], ["18%", "20%", "38%", "24%"], mono_cols=(1, 2), fs=12, rowh=30)
    systems = [("CAPS default", "Å · ps · kcal/mol · K · atm · g/cm³", True), ("SI-derived", "nm · ps · kJ/mol · K · bar · kg/m³", False), ("LAMMPS metal", "Å · ps · eV · K · bar", False)]
    sl = "".join(f'<div style="display: flex; align-items: center; gap: 10px; padding: 9px 12px; border-radius: 7px; background: {BG3 if on else "transparent"}; border: 1px solid {ACC if on else LINE}">'
                 f'<span style="width: 14px; height: 14px; border-radius: 50%; border: 2px solid {ACC if on else DIM}; flex-shrink: 0"></span><span style="display: flex; flex-direction: column; gap: 2px"><span style="font-size: 12.5px">{a}</span><span style="font-family: {MONO}; font-size: 11px; color: {DIM}">{b}</span></span></div>' for a, b, on in systems)
    fmt = [("Thousands", "12 840 atoms", "thin space, never a comma"), ("Uncertainty", "1.041 ± 0.003 g/cm³", "or 1.041(3) in exports for papers"),
           ("Significant figures", "from the uncertainty", "no more digits than the error supports"), ("Exponents", "2.3 × 10⁻¹⁰ m²/s", "never 2.3e-10 in the UI; e-notation in CSV"),
           ("Missing value", "[result] · —", "placeholder vs not applicable"), ("Minus sign", "−0.062", "U+2212 in the UI, ASCII in files")]
    fl = "".join(f'<div style="display: flex; gap: 12px; padding: 8px 0; border-bottom: 1px solid {BG2}; font-size: 12.5px"><span style="width: 140px; flex-shrink: 0; color: {MUTED}">{a}</span><span style="width: 190px; flex-shrink: 0; font-family: {MONO}">{b}</span><span style="color: {DIM}">{c}</span></div>' for a, b, c in fmt)
    left = (f'<div style="width: 420px; flex-shrink: 0; display: flex; flex-direction: column; gap: 12px">'
            + card("Unit system", col(sl, f'<span style="font-size: 11.5px; color: {DIM}; line-height: 1.5">Display only. Files keep the units of their format, and every export states them in its header.</span>', gap=8), "", 14)
            + card("Per quantity", col(row(select("Energy", "kcal/mol"), select("Length", "Å"), gap=8), row(select("Pressure", "atm"), select("Time", "ps"), gap=8), gap=8), "", 14)
            + '</div>')
    right = (f'<div style="flex-grow: 1; min-width: 0; display: flex; flex-direction: column; gap: 12px">'
             + card("Conversions used by CAPS", ct, chip("computed from SI 2019 constants", OK, BG2), 14)
             + card("Number formats", f'<div>{fl}</div>', "", 14) + '</div>')
    content = (pagehead("Settings · units & number formats", "How quantities are shown, converted and written", "")
               + f'<div style="flex-grow: 1; display: flex; gap: 14px; padding: 16px 20px; min-height: 0; overflow: hidden">{left}{right}</div>')
    return shell("Studio", "CAPS — units & number formats", "Settings › Units", content, ("<span>Settings</span>", "<span>changes apply to display only</span>"))


if __name__ == "__main__":
    import os
    os.makedirs("stage10/project", exist_ok=True)
    for name, fn in (("ChainStats", chain_stats), ("DensityCalc", density_calc), ("Units", units)):
        open(f"stage10/project/{name}.dc.html", "w").write(fn())
        print(name, "ok")
