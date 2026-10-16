#!/usr/bin/env python3
"""The six presentation figures of CAPS's ML parameter engine, from the results files and collect.py's outputs.

    plot.py DATA_DIR OUT_DIR

DATA_DIR holds parity_v1.npz, profiles_all.json (collect.py) and the caps ff complete reports and timings
(sbr/tgddm/lactic .report.json, .log with /usr/bin/time -p). The rest is read from bench/ff: the census, the v1
training summary and torsion test, the DLPNO overlay and its audit (overlay_audit.json). Writes fig1..fig6 as PNG (300 dpi) and PDF.
"""
import json, os, re, sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyBboxPatch

HERE = os.path.dirname(os.path.abspath(__file__))
BENCH = os.path.abspath(os.path.join(HERE, "..", ".."))
EV_KCAL = 23.0605478

INK, MUTED, GRID = "#1f2933", "#6b7785", "#e3e8ee"
BLUE, ORANGE, GREEN, RED, PURPLE = "#2563eb", "#ea8c00", "#16a34a", "#dc2626", "#7c3aed"
plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 11, "axes.edgecolor": MUTED, "axes.labelcolor": INK,
    "xtick.color": INK, "ytick.color": INK, "axes.spines.top": False, "axes.spines.right": False,
    "axes.titleweight": "bold", "axes.titlesize": 12, "legend.frameon": False, "mathtext.fontset": "dejavusans",
    "savefig.facecolor": "white", "figure.facecolor": "white",
})


def save(fig, out, name):
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(out, f"{name}.{ext}"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(name)


def box(ax, x, y, w, h, title, body, color):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.08", fc=color + "14",
                                ec=color, lw=1.6))
    ax.text(x + w / 2, y + h - 0.17, title, ha="center", va="top", fontsize=11.5, fontweight="bold", color=color)
    ax.text(x + w / 2, y + h - 0.52, body, ha="center", va="top", fontsize=9.3, color=INK, linespacing=1.45)


def arrow(ax, x0, y0, x1, y1):
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0), arrowprops=dict(arrowstyle="-|>", color=MUTED, lw=1.6))


def fig1_workflow(out, census, summary, n_dlpno, n_classes):
    fig, ax = plt.subplots(figsize=(14, 5.6))
    ax.set_xlim(0, 14)
    ax.set_ylim(0, 5.6)
    ax.axis("off")
    ax.text(0.1, 5.45, "Built once, on the ZIH clusters", fontsize=12.5, fontweight="bold", color=INK, va="top")
    ax.text(0.1, 2.55, "At the user's desk: caps ff complete (no QM, no Python)", fontsize=12.5, fontweight="bold",
            color=INK, va="top")
    w, h, y = 2.55, 1.85, 3.1
    top = [
        ("1  Census", f"{census['structures']} polymers, copolymers,\ncrosslink products\n"
                      f"{census['envs']:,} rotor environments\n{census['classes']} torsion classes", BLUE),
        ("2  Fragments", "cut round each rotor\nrings kept whole, H caps\nrelaxed torsion scans\n+ distorted conformers",
         BLUE),
        ("3  DFT reference", f"$\\omega$B97M-D3(BJ)/def2-TZVPPD\nGPU4PySCF on A100 / H100\n"
                             f"{summary['frames']:,} conformers\n(Alpha, Capella)", ORANGE),
        ("4  CAPS model", f"MACE, trained on Alpha\n{summary['train_frames']:,} train / "
                          f"{summary['test_frames']:,} test\nH C N O F Si S Cl\nshipped as one file", PURPLE),
        ("5  Gold standard", f"DLPNO-CCSD(T)/cc-pVTZ\n{n_dlpno:,} points on Barnard\nchecks the model,\n"
                             f"adds {n_classes} torsion classes", GREEN),
    ]
    for k, (t, b, c) in enumerate(top):
        x = 0.15 + k * 2.77
        box(ax, x, y, w, h, t, b, c)
        if k:
            arrow(ax, x - 0.2, y + h / 2, x - 0.02, y + h / 2)
    bottom = [
        ("Detect", "type the structure\nlist every bond, angle,\ntorsion with no parameter"),
        ("Cut", "fragment round one\noccurrence, same\nrules as the census"),
        ("Scan", "CAPS model in C++\n(libtorch), each point\nrelaxed, held exactly"),
        ("Fit", "harmonic bonds/angles\nOPLS $V_1$-$V_4$ torsions\nRMS gate, else reported"),
        ("Overlay", "rules with provenance\nmerged onto OPLS-AA 2024\nready for LAMMPS/GROMACS"),
    ]
    for k, (t, b) in enumerate(bottom):
        x = 0.15 + k * 2.77
        box(ax, x, 0.15, w, 1.75, t, b, INK)
        if k:
            arrow(ax, x - 0.2, 1.02, x - 0.02, 1.02)
    # the shipped model drives the scans
    ax.annotate("", xy=(2 * 2.77 + 0.15 + w / 2, 1.92), xytext=(3 * 2.77 + 0.15 + w / 2, y - 0.02),
                arrowprops=dict(arrowstyle="-|>", color=PURPLE, lw=1.8))
    ax.text(2.5 * 2.77 + 0.15 + w / 2 + 0.25, 2.55, "model file ships inside CAPS", color=PURPLE, fontsize=9.5,
            ha="left", va="center")
    save(fig, out, "fig1_workflow")


def fig2_parity(out, npz):
    d = np.load(npz)
    n = d["n"].astype(float)
    # energies per atom; plotted with each fragment's mean removed (an absolute energy mostly shows the composition),
    # the RMSE quoted on the absolute ones, as MACE's test table
    e_ml, e_ref0 = d["e_ml"] / n, d["e_ref"] / n
    rel_ref, rel_ml = np.empty_like(e_ml), np.empty_like(e_ml)
    for env in np.unique(d["env"]):
        m = d["env"] == env
        rel_ref[m] = (e_ref0[m] - e_ref0[m].mean()) * 1000
        rel_ml[m] = (e_ml[m] - e_ml[m].mean()) * 1000
    rmse_e_abs = np.sqrt(np.mean((e_ml - e_ref0) ** 2)) * 1000
    f_ref, f_ml = d["f_ref"].ravel() * 1000, d["f_ml"].ravel() * 1000
    rmse_f = np.sqrt(np.mean((f_ml - f_ref) ** 2))
    fig, axs = plt.subplots(1, 2, figsize=(11.5, 5.2))
    for ax, x, y, lab, unit, txt in (
            (axs[0], rel_ref, rel_ml, "Energy per atom", "meV/atom",
             f"RMSE {rmse_e_abs:.1f} meV/atom\n{len(n):,} conformers"),
            (axs[1], f_ref, f_ml, "Force components", "meV/$\\AA$",
             f"RMSE {rmse_f:.1f} meV/$\\AA$\n{len(f_ref):,} components")):
        lo, hi = np.percentile(np.concatenate([x, y]), [0.05, 99.95])
        pad = 0.05 * (hi - lo)
        lo, hi = lo - pad, hi + pad
        ax.hexbin(x, y, gridsize=120, bins="log", cmap="Blues", mincnt=1, extent=(lo, hi, lo, hi), linewidths=0)
        ax.plot([lo, hi], [lo, hi], color=RED, lw=1, ls="--")
        ax.set_xlim(lo, hi)
        ax.set_ylim(lo, hi)
        ax.set_aspect("equal")
        ax.set_xlabel(f"DFT ({unit})")
        ax.set_ylabel(f"CAPS model ({unit})")
        ax.set_title(lab)
        ax.text(0.04, 0.96, txt, transform=ax.transAxes, va="top", fontsize=10.5,
                bbox=dict(fc="white", ec=GRID, boxstyle="round,pad=0.4"))
    axs[0].text(0.96, 0.04, "relative to each\nfragment's mean", transform=axs[0].transAxes, ha="right",
                fontsize=8.5, color=MUTED)
    fig.suptitle(f"CAPS model on {len(np.unique(d['env']))} fragments it never saw in training",
                 fontsize=13.5, fontweight="bold", y=1.0)
    fig.tight_layout()
    save(fig, out, "fig2_parity")
    return rmse_e_abs, rmse_f


def opls(phi, v):
    p = np.radians(phi)
    return 0.5 * (v[0] * (1 + np.cos(p)) + v[1] * (1 - np.cos(2 * p)) + v[2] * (1 + np.cos(3 * p))
                  + v[3] * (1 - np.cos(4 * p)))


def fig3_profiles(out, prof, labels):
    keys = [k for k in labels if k in prof]
    fig, axs = plt.subplots(2, 3, figsize=(14, 7.8), sharex=True)
    for ax, k in zip(axs.ravel(), keys):
        p = prof[k]
        phi = np.array(p["phi"], float)
        o = np.argsort(phi)
        phi = phi[o]
        g = lambda n: np.array(p[n])[o]
        rms = lambda a, b: float(np.sqrt(np.mean((a - b - np.mean(a - b)) ** 2)))
        ax.plot(phi, g("opls_without"), color=MUTED, ls=":", lw=1.6, label="OPLS-AA 2024 in CAPS before (missing, or an MP2 fit)")
        ax.plot(phi, g("dlpno"), "o", color=INK, ms=5.5, label="DLPNO-CCSD(T)", zorder=3)
        ax.plot(phi, g("caps_ml"), "-", color=PURPLE, lw=2.2, label="CAPS model")
        ax.plot(phi, g("opls_fitted"), "--", color=GREEN, lw=2, label="OPLS-AA 2024 + CAPS fit")
        ax.set_title(f"{k}\n{labels[k]}", fontsize=11)
        ax.text(0.03, 0.97, f"model vs DLPNO {rms(g('caps_ml'), g('dlpno')):.2f}\n"
                            f"fit vs DLPNO {rms(g('opls_fitted'), g('dlpno')):.2f} kcal/mol",
                transform=ax.transAxes, va="top", fontsize=8.8, bbox=dict(fc="white", ec=GRID, boxstyle="round"))
        top = max(max(g("dlpno")), max(g("caps_ml")), max(g("opls_fitted")), max(g("opls_without")))
        ax.set_ylim(-0.04 * top, 1.32 * top)   # room for the RMS box above the highest barrier
        ax.set_xticks([-180, -90, 0, 90, 180])
        ax.grid(color=GRID, lw=0.8)
    for ax in axs[1]:
        ax.set_xlabel("dihedral angle (degrees)")
    for ax in axs[:, 0]:
        ax.set_ylabel("relative energy (kcal/mol)")
    h, l = axs[0, 0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=4, bbox_to_anchor=(0.5, -0.03), fontsize=11)
    fig.suptitle("Torsions OPLS-AA 2024 lacks, refitted by CAPS to DLPNO-CCSD(T)", fontsize=13.5, fontweight="bold")
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    save(fig, out, "fig3_profiles")


def fig4_accuracy(out, tors, audit):
    s = tors["summary"]
    fig, axs = plt.subplots(1, 2, figsize=(13, 5.2), gridspec_kw={"width_ratios": [1, 1.25]})
    ax = axs[0]
    metrics = [("median", "rmse_median"), ("90th percentile", "rmse_p90")]
    x = np.arange(len(metrics))
    for k, (name, col) in enumerate((("CAPS model", PURPLE), ("AIMNet2", ORANGE))):
        key = "caps_mace" if k == 0 else "aimnet2"
        vals = [s[key][m] for _, m in metrics]
        b = ax.bar(x + (k - 0.5) * 0.36, vals, 0.34, color=col, label=name)
        ax.bar_label(b, fmt="%.2f", padding=3, fontsize=10)
    ax.set_xticks(x, [m for m, _ in metrics])
    ax.set_ylabel("torsion profile RMSE vs DFT (kcal/mol)")
    n = s["caps_mace"]["environments"]
    ax.set_title(f"Held-out torsion scans ({n} environments)")
    ax.legend(loc="upper left")
    ax.text(0.98, 0.97, f"within 0.5 kcal/mol:\nCAPS {s['caps_mace']['within_0.5']}/{n}, "
                        f"AIMNet2 {s['aimnet2']['within_0.5']}/{n}", transform=ax.transAxes, ha="right", va="top",
            fontsize=9.5, color=INK)
    ax.set_ylim(0, max(s["aimnet2"]["rmse_p90"], s["caps_mace"]["rmse_p90"]) * 1.45)
    ax.grid(axis="y", color=GRID)
    ax.set_axisbelow(True)
    ax = axs[1]
    # every usable DLPNO-CCSD(T) scan of the library: OPLS-AA 2024 alone, with CAPS's DLPNO overlay, and CAPS's model
    sc = [x for x in audit["scans"] if x["in_overlay"]]
    ml = np.array([x["rms_caps_model"] for x in audit["scans"]])
    alone = np.array([x["rms_opls_alone"] for x in sc])
    fitted = np.array([x["rms_opls_with_overlay"] for x in sc])
    bins = np.arange(0, 3.01, 0.15)
    ax.hist(np.clip(alone, 0, 2.99), bins=bins, color=MUTED, alpha=0.55, label=f"OPLS-AA 2024 alone ({len(sc)} classes)")
    ax.hist(fitted, bins=bins, color=GREEN, alpha=0.85, label="OPLS-AA 2024 + CAPS DLPNO fit")
    ax.hist(ml, bins=bins, histtype="step", color=PURPLE, lw=2.2, label=f"CAPS model ({len(ml)} scans)")
    ax.axvline(1.0, color=RED, ls="--", lw=1.2)
    ax.legend(loc="upper right", fontsize=9.5)
    ax.text(0.98, 0.55, f"within 1 kcal/mol:\nOPLS-AA 2024 alone {int((alone <= 1).sum())}/{len(sc)}\n"
                        f"with CAPS fit {int((fitted <= 1).sum())}/{len(sc)}\nCAPS model {int((ml <= 1).sum())}/{len(ml)}",
            transform=ax.transAxes, ha="right", va="top", fontsize=9.5, color=INK)
    ax.set_xlabel("RMSE vs DLPNO-CCSD(T) per torsion scan (kcal/mol)")
    ax.set_ylabel("scans")
    ax.set_title("Against the gold standard, DLPNO-CCSD(T)")
    ax.grid(axis="y", color=GRID)
    ax.set_axisbelow(True)
    fig.tight_layout()
    save(fig, out, "fig4_accuracy")


def seconds(log):
    m = re.search(r"^real\s+([\d.]+)", open(log).read(), re.M)
    if m:
        return float(m.group(1))
    m = re.search(r"(\d+):([\d.]+) total", open(log).read())
    return 60 * float(m.group(1)) + float(m.group(2))


def fig5_speed(out, py_s, cpp_s, fc):
    fig, axs = plt.subplots(1, 2, figsize=(12.5, 4.6), gridspec_kw={"width_ratios": [1, 1]})
    ax = axs[0]
    b = ax.barh(["Python (MACE + ASE)", "CAPS C++ engine"], [py_s / 60, cpp_s / 60], color=[MUTED, PURPLE],
                height=0.55)
    ax.bar_label(b, labels=[f"{py_s / 60:.1f} min", f"{cpp_s / 60:.1f} min"], padding=4, fontsize=10.5)
    ax.invert_yaxis()
    ax.set_xlabel("wall time (minutes, Apple M5 laptop CPU)")
    ax.set_title(f"Same model, same SBR scan: {py_s / cpp_s:.0f}x faster")
    ax.set_xlim(0, py_s / 60 * 1.25)
    ax = axs[1]
    names = list(fc)
    b = ax.barh(names, [fc[k]["seconds"] / 60 for k in names], color=PURPLE, height=0.55)
    ax.bar_label(b, labels=[f"{fc[k]['seconds'] / 60:.1f} min, {fc[k]['fitted']} term{'' if fc[k]['fitted'] == 1 else 's'} fitted" for k in names],
                 padding=4, fontsize=10.5)
    ax.invert_yaxis()
    ax.set_xlabel("wall time (minutes)")
    ax.set_title("caps ff complete, start to finished overlay")
    ax.set_xlim(0, max(fc[k]["seconds"] for k in names) / 60 * 1.6)
    for a in axs:
        a.grid(axis="x", color=GRID)
        a.set_axisbelow(True)
    fig.tight_layout()
    save(fig, out, "fig5_speed")


def fig6_coverage(out, census, fc):
    fig, axs = plt.subplots(1, 2, figsize=(13, 4.9), gridspec_kw={"width_ratios": [1.15, 1]})
    ax = axs[0]
    rows = [("torsion classes", census["torsion"]), ("rotor environments", census["rotor"])]
    for i, (name, (pub, fit, miss)) in enumerate(rows):
        tot = pub + fit + miss
        left = 0
        for v, c, lab in ((pub, BLUE, "published OPLS-AA 2024"), (fit, GREEN, "CAPS-fitted"),
                          (miss, RED, "missing")):
            ax.barh(i, v / tot * 100, left=left, color=c, height=0.55, label=lab if i == 0 else None)
            if v / tot > 0.04:
                ax.text(left + v / tot * 50, i, f"{v}", ha="center", va="center", color="white", fontsize=10,
                        fontweight="bold")
            left += v / tot * 100
        ax.text(101, i, f"{tot:,}", va="center", fontsize=10)
    ax.set_yticks([0, 1], [r[0] for r in rows])
    ax.invert_yaxis()
    ax.set_xlim(0, 112)
    ax.set_xlabel("% of what the census found")
    ax.set_title(f"Census of {census['structures']} structures")
    ax.legend(loc="lower center", bbox_to_anchor=(0.45, -0.42), ncol=3, fontsize=9.5)
    ax = axs[1]
    names = list(fc)
    x = np.arange(len(names))
    before = [fc[k]["before"] for k in names]
    after = [fc[k]["after"] for k in names]
    b1 = ax.bar(x - 0.19, before, 0.36, color=RED, label="missing before")
    b2 = ax.bar(x + 0.19, after, 0.36, color=GREEN, label="missing after caps ff complete")
    ax.bar_label(b1, padding=2)
    ax.bar_label(b2, padding=2)
    ax.set_xticks(x, [n + ("\n" + fc[n]["note"] if fc[n].get("note") else "") for n in names])
    ax.set_ylabel("missing bonded terms")
    ax.set_ylim(0, max(before) * 1.35)
    ax.legend(loc="upper right", fontsize=9.5)
    ax.set_title("Per structure")
    for a in axs:
        a.set_axisbelow(True)
    axs[1].grid(axis="y", color=GRID)
    fig.tight_layout()
    save(fig, out, "fig6_coverage")


def main():
    data, out = sys.argv[1], sys.argv[2]
    os.makedirs(out, exist_ok=True)
    md = open(os.path.join(BENCH, "census", "oplsaa2024", "summary.md")).read()
    row = lambda name: [int(v) for v in re.search(rf"\| {name} \| (\d+) \| (\d+) \| (\d+) \|", md).groups()]
    census = {"structures": int(re.search(r"Structures typed: (\d+)", md).group(1)),
              "torsion": row("torsion classes"), "rotor": row("rotor environments")}
    census["envs"] = sum(census["rotor"])
    census["classes"] = sum(census["torsion"])
    summary = json.load(open(os.path.join(BENCH, "mltrain", "mace_v1", "summary.json")))
    plan = json.load(open(os.path.join(BENCH, "mltrain", "dlpno_library", "plan.json")))
    n_dlpno = sum(len(s["points"]) for s in plan["scans"])
    ov = json.load(open(os.path.join(BENCH, "..", "..", "data", "forcefields", "oplsaa2024-dlpno.json")))
    fig1_workflow(out, census, summary, n_dlpno, sum(1 for d in ov["dihedrals"] if not d.get("fit", {}).get("zeroed")))
    fig2_parity(out, os.path.join(data, "parity_v1.npz"))
    labels = json.load(open(os.path.join(data, "profile_labels.json")))
    fig3_profiles(out, json.load(open(os.path.join(data, "profiles_all.json"))), labels)
    fig4_accuracy(out, json.load(open(os.path.join(BENCH, "mltrain", "mace_v1", "torsions.json"))),
                  json.load(open(os.path.join(BENCH, "mltrain", "dlpno_library", "overlay_audit.json"))))
    fc = json.load(open(os.path.join(data, "complete.json")))
    fig5_speed(out, fc.pop("_python_s"), fc.pop("_native_s"), fc)
    fig6_coverage(out, census, fc)


if __name__ == "__main__":
    main()
