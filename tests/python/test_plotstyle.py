"""caps.plotstyle (the journal look of CAPS's matplotlib figures): run by ctest as python_plotstyle.

Checks that the settings written into analyze_kappa.py (core/src/phono3py.cpp, a script for clusters where CAPS is not
installed) are caps.plotstyle.RC, that axis titles come out as "quantity / unit", and, where matplotlib is installed, that
a figure is written as a 300 dpi PNG and a vector PDF of the column width.
"""
import ast
import os
import re
import sys
import tempfile

from caps import plotstyle as ps


def check(cond, what):
    if not cond:
        print("FAIL", what)
        sys.exit(1)
    print("ok  ", what)


root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
src = open(os.path.join(root, "core", "src", "phono3py.cpp"), encoding="utf-8").read()
m = re.search(r"\nRC = (\{.*?\n\})\nCOLOURS = (\[.*?\])\nSINGLE, DOUBLE = ([\d.]+), ([\d.]+)\n", src, re.S)
check(m is not None, "analyze_kappa.py carries RC, COLOURS and the column widths")
check(ast.literal_eval(m.group(1)) == ps.RC, "analyze_kappa.py's RC is caps.plotstyle.RC")
check(ast.literal_eval(m.group(2)) == list(ps.COLOURS), "analyze_kappa.py's colours are caps.plotstyle.COLOURS")
check((float(m.group(3)), float(m.group(4))) == (ps.SINGLE, ps.DOUBLE), "column widths 3.35 and 7 in")

check(ps.RC["xtick.direction"] == "in" and ps.RC["ytick.direction"] == "in" and ps.RC["xtick.top"] and ps.RC["ytick.right"],
      "inward ticks on all four sides")
check(ps.RC["xtick.minor.visible"] and ps.RC["ytick.minor.visible"], "minor ticks")
check(not ps.RC["legend.frameon"] and ps.RC["savefig.dpi"] == 300 and ps.RC["pdf.fonttype"] == 42, "no legend frame, 300 dpi, fonts embedded")
check(8 <= ps.RC["axes.labelsize"] <= 9, "8-9 pt axis titles")
check(ps.label("T", "K") == "$\\mathit{T}$ / K", "label T / K")
check(ps.label("kappa", "W m-1 K-1") == "$\\mathit{\\kappa}$ / W m$^{−1}$ K$^{−1}$", "label kappa with raised powers")
check(ps.label("E_pot", "kcal mol-1") == "$\\mathit{E}_{\\mathrm{pot}}$ / kcal mol$^{−1}$", "label with an upright subscript")
check(ps.label("frequency", "THz") == "frequency / THz", "a word stays upright")
check(ps.label("MSD", "Å2") == "MSD / Å$^{2}$", "an abbreviation upright, Å² raised")
check("marker" in ps.line(0, 40) and "marker" not in ps.line(0, 400), "points on a sparse series only")

try:
    import matplotlib  # noqa: F401
except ImportError:
    print("skip  figure writing: no matplotlib in this Python")
    sys.exit(0)

fig, ax = ps.figure("single")
x = [0, 1, 2, 3]
ax.plot(x, [1, 3, 2, 4], **ps.line(0, len(x)))
ax.set_xlabel(ps.label("t", "ps"))
ax.set_ylabel(ps.label("T", "K"))
d = tempfile.mkdtemp(prefix="caps_plotstyle_")
paths = ps.save(fig, os.path.join(d, "check"))
check([os.path.basename(p) for p in paths] == ["check.png", "check.pdf"], "PNG and PDF written")
png = open(paths[0], "rb").read()
check(png[:8] == b"\x89PNG\r\n\x1a\n", "a PNG")
w = int.from_bytes(png[16:20], "big")
check(abs(w - 3.35 * 300) < 0.12 * 3.35 * 300, f"PNG {w} px wide: one column at 300 dpi (tight bounding box)")
check(open(paths[1], "rb").read(5) == b"%PDF-", "a PDF")
print("all plotstyle checks passed")
