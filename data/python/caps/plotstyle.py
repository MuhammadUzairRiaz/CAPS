"""One journal look for every figure CAPS writes with matplotlib.

    from caps import plotstyle as ps
    fig, ax = ps.figure("single")                      # 3.35 in wide; "double" is 7 in
    ax.plot(t, T, **ps.line(0, len(t)))                # points on a sparse series, a plain line on a dense one
    ax.set_xlabel(ps.label("t", "ps")); ax.set_ylabel(ps.label("T", "K"))
    ps.save(fig, "temperature")                        # temperature.png (300 dpi) and temperature.pdf (vector)

The look is the one of the Studio's plots: a closed frame with inward major and minor ticks on all four sides, axis titles
as "quantity / unit" with the symbol in italics and the unit upright (mathtext), Okabe and Ito's colour-blind safe
colours, no legend frame, 8 pt text at print size. matplotlib is imported only when a figure is made, so ``import caps``
does not need it.

Scripts CAPS writes for a cluster (``caps phono3py`` writes analyze_kappa.py) run where CAPS is not installed and carry
the same settings; ``RC`` is the source the tests compare them with.
"""
import os
import re
import sys

__all__ = ["SINGLE", "DOUBLE", "COLOURS", "MARKERS", "MARKER_LIMIT", "RC", "apply", "figure", "symbol", "label", "line",
           "band", "panel_labels", "save"]

SINGLE = 3.35          # in: one column of a two-column journal page
DOUBLE = 7.0           # in: the full page width
WIDTHS = {"single": SINGLE, "onehalf": 5.0, "double": DOUBLE}

# Okabe & Ito (2008): blue, vermillion, bluish green, reddish purple, orange, sky blue, black, yellow
COLOURS = ("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#000000", "#F0E442")
MARKERS = ("o", "s", "^", "D", "v", "<", ">", "p")
LINESTYLES = ("-", "--", "-.", ":")
# up to this many samples a measured series shows each sample as a point (as the Studio's plots do)
MARKER_LIMIT = 150

RC = {
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans", "DejaVu Sans"],
    "font.size": 8.0,
    "axes.labelsize": 8.0,
    "axes.titlesize": 8.0,
    "xtick.labelsize": 7.5,
    "ytick.labelsize": 7.5,
    "legend.fontsize": 7.5,
    "mathtext.default": "regular",
    "mathtext.fontset": "custom",
    "mathtext.rm": "sans",
    "mathtext.it": "sans:italic",
    "mathtext.bf": "sans:bold",
    "axes.linewidth": 0.6,
    "axes.edgecolor": "black",
    "axes.labelcolor": "black",
    "axes.grid": False,
    "axes.unicode_minus": True,
    "axes.formatter.use_mathtext": True,
    "axes.formatter.limits": [-3, 4],
    "axes.labelpad": 3.0,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "xtick.minor.visible": True,
    "ytick.minor.visible": True,
    "xtick.major.size": 3.5,
    "ytick.major.size": 3.5,
    "xtick.minor.size": 2.0,
    "ytick.minor.size": 2.0,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.minor.width": 0.45,
    "ytick.minor.width": 0.45,
    "xtick.major.pad": 3.0,
    "ytick.major.pad": 3.0,
    "xtick.color": "black",
    "ytick.color": "black",
    "lines.linewidth": 1.0,
    "lines.markersize": 3.0,
    "lines.markeredgewidth": 0.6,
    "errorbar.capsize": 1.5,
    "legend.frameon": False,
    "legend.handlelength": 1.8,
    "legend.borderaxespad": 0.4,
    "legend.labelspacing": 0.3,
    "figure.dpi": 150,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "savefig.dpi": 300,
    "savefig.facecolor": "white",
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.02,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
}


def apply():
    """Sets the journal look in matplotlib (the Agg backend unless pyplot is already in use) and returns pyplot."""
    import matplotlib
    if "matplotlib.pyplot" not in sys.modules:
        matplotlib.use("Agg")
    from cycler import cycler
    rc = dict(RC)
    rc["axes.prop_cycle"] = cycler(color=list(COLOURS))
    matplotlib.rcParams.update(rc)
    import matplotlib.pyplot as plt
    return plt


def figure(width="single", aspect=0.75, nrows=1, ncols=1, height=None, **kw):
    """A figure the width of a journal column ("single", 3.35 in), one and a half columns ("onehalf") or the page
    ("double", 7 in), or a width in inches; height = width × aspect unless given. Returns (fig, ax or axes)."""
    plt = apply()
    w = WIDTHS[width] if isinstance(width, str) else float(width)
    return plt.subplots(nrows, ncols, figsize=(w, height if height is not None else w * aspect), **kw)


GREEK = {"alpha", "beta", "gamma", "delta", "Delta", "epsilon", "eta", "theta", "kappa", "lambda", "mu", "nu", "rho",
         "sigma", "tau", "phi", "chi", "psi", "omega"}


def symbol(q):
    """A quantity's symbol as italic mathtext: "T" gives T in italics, "kappa" gives κ, "E_pot" gives E with an upright
    subscript "pot", "v_x" an italic x. Words (more than two letters, not Greek: "MSD", "frequency") and text holding "$" come back as given."""
    if "$" in q:
        return q
    base, _, sub = q.partition("_")
    if base in GREEK:
        body = "\\mathit{\\" + base + "}"
    elif 0 < len(base) <= 2 and base.isalpha():
        body = "\\mathit{" + base + "}"
    else:
        return q
    if sub:
        body += "_{\\mathrm{" + sub + "}}" if len(sub) > 1 else "_{\\mathit{" + sub + "}}"
    return "$" + body + "$"


def _unit(s):
    """A plain-text unit as mathtext with raised powers: "W m-1 K-1" → W m⁻¹ K⁻¹ (as $^{-1}$), "Å2" → Å²."""
    return re.sub(r"(?<=[A-Za-zÅ°])(-?\d+)(?![\d.])", lambda m: "$^{" + m.group(1).replace("-", "−") + "}$", s)


def label(quantity, unit=None):
    """An axis title "quantity / unit" (IUPAC): the symbol in italics, the unit upright with its powers raised.
    label("kappa", "W m-1 K-1") is κ / W m⁻¹ K⁻¹; label("t", "ps") is t / ps; label("frequency", "THz") stays upright."""
    q = symbol(quantity)
    return q if not unit else q + " / " + _unit(unit)


def line(k=0, n=None, dashed=False, **kw):
    """Keyword arguments for ax.plot of series k with n samples: its colour, and an open marker on each sample while the
    series has at most MARKER_LIMIT points (a line between sparse samples is only a guide to the eye)."""
    d = {"color": COLOURS[k % len(COLOURS)], "lw": 1.0, "ls": LINESTYLES[k % len(LINESTYLES)] if dashed else "-"}
    if n is not None and n <= MARKER_LIMIT:
        d.update(marker=MARKERS[k % len(MARKERS)], ms=3.0, mfc="white", mew=0.7)
    d.update(kw)
    return d


def band(ax, x, lo, hi, k=0, **kw):
    """An uncertainty band (mean ± error) under series k."""
    return ax.fill_between(x, lo, hi, color=COLOURS[k % len(COLOURS)], alpha=0.2, lw=0, **kw)


def panel_labels(axes, x=-0.2, y=1.0):
    """(a), (b) … in bold at the top left of each panel."""
    for k, a in enumerate(getattr(axes, "flat", axes)):
        a.text(x, y, "(" + "abcdefghijklmnopqrstuvwxyz"[k] + ")", transform=a.transAxes, fontweight="bold", va="bottom", ha="left")


def save(fig, stem, formats=("png", "pdf"), dpi=300):
    """Writes stem.png at dpi (300) and stem.pdf (vector, TrueType fonts embedded), or the given formats ("svg" too);
    a stem ending in one of those extensions puts that format first. Closes the figure and returns the paths."""
    import matplotlib.pyplot as plt
    stem = os.fspath(stem)
    base, ext = os.path.splitext(stem)
    ext = ext.lower().lstrip(".")
    if ext in ("png", "pdf", "svg", "eps", "tif", "tiff"):
        stem, formats = base, (ext,) + tuple(f for f in formats if f != ext)
    d = os.path.dirname(stem)
    if d:
        os.makedirs(d, exist_ok=True)
    out = []
    for f in formats:
        p = stem + "." + f
        if f in ("png", "tif", "tiff"):
            fig.savefig(p, dpi=dpi)
        else:
            fig.savefig(p)
        out.append(p)
    plt.close(fig)
    return out
