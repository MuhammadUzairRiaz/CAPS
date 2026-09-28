"""Python render overlays (design/boards/RenderOverlays "Python overlay · draw your own").

A script draws on the rendered image with a small canvas; CAPS Studio paints what it drew over the render (and over the
render-frame guide in the view), one call per frame::

    from caps.overlay import overlay

    @overlay
    def draw(canvas, data):
        rho = data.attributes.get("Density", 0)
        canvas.text(24, canvas.height - 60, f"ρ = {rho:.3f} g/cm³", size=28, colour="#FFFFFF")
        rdf = data.tables.get("rdf")
        if rdf:
            canvas.plot(rdf.column(0), rdf.column(1), box=(canvas.width - 520, 40, 480, 300), title="g(r)")

Coordinates are pixels of the output image, x to the right and y down from the top-left corner. data holds the frame
(number, count), the title, the global attributes of the Visualize pipeline on that frame (Particles, Density,
CoordinationAnalysis.mean …) and its data tables. Colours are "#RRGGBB" or "#RRGGBBAA".

The Studio runs ``python3 -m caps.overlay SCRIPT`` with the frame's data as JSON on stdin; the drawing comes back as
JSON commands on stdout. Anything the script prints goes to the Render page's console.
"""
import json
import math
import runpy
import sys

__all__ = ["overlay", "Canvas", "Data", "Table"]

_drawers = []


def overlay(fn):
    """Marks the function that draws: fn(canvas, data)."""
    _drawers.append(fn)
    return fn


class Table:
    """A data table of the pipeline: columns (names) and rows (lists of numbers)."""

    def __init__(self, name, columns, rows):
        self.name, self.columns, self.rows = name, list(columns), [list(r) for r in rows]

    def column(self, k):
        """Column k (an index or a name) as a list."""
        if isinstance(k, str):
            k = self.columns.index(k)
        return [r[k] for r in self.rows if k < len(r)]

    def __len__(self):
        return len(self.rows)


class Data:
    """What the overlay knows about the frame."""

    def __init__(self, d):
        self.frame = int(d.get("frame", 0))
        self.frames = int(d.get("frames", 1))
        self.title = d.get("title", "")
        self.attributes = dict(d.get("attributes", {}))
        self.tables = {name: Table(name, t.get("columns", []), t.get("rows", [])) for name, t in d.get("tables", {}).items()}


class Canvas:
    """Draws in the output image's pixels; every call is recorded as a command for the Studio to paint."""

    def __init__(self, width, height, dark=True):
        self.width, self.height, self.dark = float(width), float(height), bool(dark)
        self.commands = []
        # the ink that reads on this background
        self.ink = "#E9ECEF" if dark else "#141413"

    def _c(self, colour):
        return self.ink if colour is None else str(colour)

    def text(self, x, y, s, size=16, colour=None, align="left", bold=False, mono=False):
        """Text with its top-left (align left), top-centre (centre) or top-right (right) corner at x, y."""
        self.commands.append({"op": "text", "x": float(x), "y": float(y), "s": str(s), "size": float(size), "colour": self._c(colour),
                              "align": align, "bold": bool(bold), "mono": bool(mono)})

    def line(self, x1, y1, x2, y2, colour=None, width=2.0):
        self.commands.append({"op": "line", "x1": float(x1), "y1": float(y1), "x2": float(x2), "y2": float(y2), "colour": self._c(colour), "width": float(width)})

    def polyline(self, points, colour=None, width=2.0, closed=False):
        """A path through [(x, y), …]."""
        pts = [[float(x), float(y)] for x, y in points]
        if len(pts) >= 2:
            self.commands.append({"op": "polyline", "points": pts, "colour": self._c(colour), "width": float(width), "closed": bool(closed)})

    def rect(self, x, y, w, h, fill=None, stroke=None, width=1.0, radius=0.0):
        self.commands.append({"op": "rect", "x": float(x), "y": float(y), "w": float(w), "h": float(h), "fill": fill, "stroke": stroke,
                              "width": float(width), "radius": float(radius)})

    def circle(self, x, y, r, fill=None, stroke=None, width=1.0):
        self.commands.append({"op": "circle", "x": float(x), "y": float(y), "r": float(r), "fill": fill, "stroke": stroke, "width": float(width)})

    def plot(self, xs, ys, box, colour="#D98A1E", title="", xlabel="", ylabel="", background=True, marker=None):
        """An inset line plot of ys against xs in box = (x, y, w, h), with its axes, ticks and title."""
        x0, y0, w, h = (float(v) for v in box)
        pts = [(float(a), float(b)) for a, b in zip(xs, ys) if math.isfinite(float(a)) and math.isfinite(float(b))]
        if background:
            self.rect(x0, y0, w, h, fill="#0B0D0FCC" if self.dark else "#FFFFFFCC", stroke=None, radius=6)
        if not pts:
            return
        ml, mr, mt, mb = 46.0, 10.0, 24.0 if title else 10.0, 30.0 if xlabel else 22.0
        px, py, pw, ph = x0 + ml, y0 + mt, max(1.0, w - ml - mr), max(1.0, h - mt - mb)
        lo_x, hi_x = min(p[0] for p in pts), max(p[0] for p in pts)
        lo_y, hi_y = min(p[1] for p in pts), max(p[1] for p in pts)
        if hi_x == lo_x:
            lo_x, hi_x = lo_x - 0.5, hi_x + 0.5
        if hi_y == lo_y:
            lo_y, hi_y = lo_y - 0.5, hi_y + 0.5
        pad = 0.05 * (hi_y - lo_y)
        lo_y, hi_y = lo_y - pad, hi_y + pad
        X = lambda v: px + (v - lo_x) / (hi_x - lo_x) * pw
        Y = lambda v: py + (1 - (v - lo_y) / (hi_y - lo_y)) * ph
        muted = "#A5ABB1" if self.dark else "#5A6168"
        self.rect(px, py, pw, ph, stroke=muted, width=1)
        for t in _ticks(lo_x, hi_x):
            self.line(X(t), py + ph, X(t), py + ph + 4, colour=muted, width=1)
            self.text(X(t), py + ph + 5, _fmt(t), size=11, colour=muted, align="centre", mono=True)
        for t in _ticks(lo_y, hi_y):
            self.line(px - 4, Y(t), px, Y(t), colour=muted, width=1)
            self.text(px - 6, Y(t) - 7, _fmt(t), size=11, colour=muted, align="right", mono=True)
        self.polyline([(X(a), Y(b)) for a, b in pts], colour=colour, width=2)
        if marker is not None:   # a marker at x = marker (the current time, say)
            self.line(X(marker), py, X(marker), py + ph, colour=self.ink, width=1)
        if title:
            self.text(x0 + ml, y0 + 5, title, size=13, bold=True)
        if xlabel:
            self.text(px + pw / 2, y0 + h - 16, xlabel, size=11, colour=muted, align="centre")
        if ylabel:
            self.text(x0 + 4, py - 16, ylabel, size=11, colour=muted)


def _ticks(lo, hi):
    raw = (hi - lo) / 4
    mag = 10 ** math.floor(math.log10(raw))
    r = raw / mag
    step = (1 if r < 1.5 else 2 if r < 3 else 5 if r < 7 else 10) * mag
    t, out = math.ceil(lo / step) * step, []
    while t <= hi + 1e-9 * step:
        out.append(0.0 if abs(t) < 1e-12 * step else t)
        t += step
    return out


def _fmt(v):
    return f"{v:.4g}"


def main(argv):
    if len(argv) < 2:
        print("usage: python3 -m caps.overlay SCRIPT < frame.json", file=sys.stderr)
        return 2
    d = json.loads(sys.stdin.read() or "{}")
    canvas = Canvas(d.get("width", 1920), d.get("height", 1080), d.get("dark", True))
    data = Data(d)
    real_stdout = sys.stdout
    sys.stdout = sys.stderr   # what the script prints goes to the console, not into the commands
    try:
        g = runpy.run_path(argv[1], run_name="__caps_overlay__")
        # run as python -m caps.overlay, this module is __main__: the script's "from caps.overlay import overlay" marks
        # its functions in the imported copy
        other = sys.modules.get("caps.overlay")
        marked = list(_drawers) + (list(other._drawers) if other is not None and other is not sys.modules.get(__name__) else [])
        drawers = marked or [g[n] for n in ("draw", "render") if callable(g.get(n))]
        if not drawers:
            raise SystemExit("no overlay: mark a function draw(canvas, data) with @overlay")
        for fn in drawers:
            fn(canvas, data)
    finally:
        sys.stdout = real_stdout
    json.dump({"commands": canvas.commands}, sys.stdout)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
