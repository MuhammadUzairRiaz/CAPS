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
            canvas.plot(rdf.column(0), rdf.column(1), box=(canvas.width - 520, 40, 480, 300), xlabel="r / Å", ylabel="g(r)")

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

    def plot(self, xs, ys, box, colour=None, title="", xlabel="", ylabel="", background=True, marker=None, markers=None):
        """An inset plot of ys against xs in box = (x, y, w, h), drawn the way the Studio's plots are: a closed frame,
        inward major and minor ticks on all four sides at round values (1, 2, 2.5 or 5 times a power of ten), the axis
        titles as "quantity / unit" (ylabel above the frame, as an inset has no room to turn it), every finite (x, y) drawn
        as it is. A series of up to 150 points also gets an open circle on each point (markers=True/False to choose).
        colour: Okabe and Ito's blue unless given; marker: a vertical line at that x (the current time, say)."""
        x0, y0, w, h = (float(v) for v in box)
        pts = [(float(a), float(b)) for a, b in zip(xs, ys) if math.isfinite(float(a)) and math.isfinite(float(b))]
        if background:
            self.rect(x0, y0, w, h, fill="#0B0D0FCC" if self.dark else "#FFFFFFCC", stroke=None, radius=6)
        if not pts:
            return
        colour = colour or ("#56B4E9" if self.dark else "#0072B2")
        fs = max(9.0, min(14.0, h / 22.0))            # tick labels; titles one size up
        tick = max(3.0, fs * 0.45)
        lo_x, hi_x = min(p[0] for p in pts), max(p[0] for p in pts)
        lo_y, hi_y = min(p[1] for p in pts), max(p[1] for p in pts)
        if hi_x == lo_x:
            lo_x, hi_x = lo_x - 0.5, hi_x + 0.5
        if hi_y == lo_y:
            lo_y, hi_y = lo_y - 0.5, hi_y + 0.5
        pad = 0.05 * (hi_y - lo_y)
        lo_y, hi_y = lo_y - pad, hi_y + pad
        xstep, ystep = _step(lo_x, hi_x), _step(lo_y, hi_y)
        # the axis ends on a minor tick, so the frame does not stop at an odd value
        lo_x, hi_x = _snap(lo_x, hi_x, xstep)
        lo_y, hi_y = _snap(lo_y, hi_y, ystep)
        ylab_w = max(len(_fmt(t, ystep)) for t in _ticks(lo_y, hi_y, ystep)) * fs * 0.6
        ml, mr = ylab_w + fs * 0.9 + 4, 10.0
        mt = (fs * 1.3 if title else 4.0) + (fs * 1.4 if ylabel else 6.0)
        mb = fs * 1.5 + (fs * 1.5 if xlabel else 2.0)
        px, py, pw, ph = x0 + ml, y0 + mt, max(1.0, w - ml - mr), max(1.0, h - mt - mb)
        X = lambda v: px + (v - lo_x) / (hi_x - lo_x) * pw
        Y = lambda v: py + (1 - (v - lo_y) / (hi_y - lo_y)) * ph
        ink = self.ink
        lw = max(1.0, fs / 11.0)
        for t in _ticks(lo_x, hi_x, xstep / _minor(xstep)):
            self.line(X(t), py + ph, X(t), py + ph - tick * 0.55, colour=ink, width=lw)
            self.line(X(t), py, X(t), py + tick * 0.55, colour=ink, width=lw)
        for t in _ticks(lo_y, hi_y, ystep / _minor(ystep)):
            self.line(px, Y(t), px + tick * 0.55, Y(t), colour=ink, width=lw)
            self.line(px + pw, Y(t), px + pw - tick * 0.55, Y(t), colour=ink, width=lw)
        for t in _ticks(lo_x, hi_x, xstep):
            self.line(X(t), py + ph, X(t), py + ph - tick, colour=ink, width=lw)
            self.line(X(t), py, X(t), py + tick, colour=ink, width=lw)
            self.text(X(t), py + ph + fs * 0.35, _fmt(t, xstep), size=fs, colour=ink, align="centre")
        for t in _ticks(lo_y, hi_y, ystep):
            self.line(px, Y(t), px + tick, Y(t), colour=ink, width=lw)
            self.line(px + pw, Y(t), px + pw - tick, Y(t), colour=ink, width=lw)
            self.text(px - fs * 0.4, Y(t) - fs * 0.62, _fmt(t, ystep), size=fs, colour=ink, align="right")
        self.polyline([(X(a), Y(b)) for a, b in pts], colour=colour, width=max(1.5, fs / 7.0))
        if markers if markers is not None else len(pts) <= 150:
            r = max(2.5, fs * 0.28)
            bg = "#0B0D0F" if self.dark else "#FFFFFF"
            for a, b in pts:
                self.circle(X(a), Y(b), r, fill=bg, stroke=colour, width=lw)
        if marker is not None and lo_x <= marker <= hi_x:   # a vertical line at x = marker (the current time, say)
            self.line(X(marker), py, X(marker), py + ph, colour=ink, width=lw)
        self.rect(px, py, pw, ph, stroke=ink, width=lw)
        if title:
            self.text(x0 + 6, y0 + 3, title, size=fs + 1, bold=True)
        if ylabel:
            self.text(px, py - fs * 1.45, ylabel, size=fs + 1)
        if xlabel:
            self.text(px + pw / 2, py + ph + fs * 1.55, xlabel, size=fs + 1, align="centre")


def _step(lo, hi, target=4):
    """A round major step (1, 2, 2.5 or 5 times a power of ten) giving about target intervals over [lo, hi]."""
    raw = (hi - lo) / target
    mag = 10 ** math.floor(math.log10(raw))
    r = raw / mag
    return (1 if r < 1.5 else 2 if r < 2.25 else 2.5 if r < 3.5 else 5 if r < 7.5 else 10) * mag


def _minor(step):
    """Minor ticks per major step: 4 for a step of 2, else 5."""
    f = step / 10 ** math.floor(math.log10(step))
    return 4 if abs(f - 2) < 1e-9 else 5


def _snap(lo, hi, step):
    m = step / _minor(step)
    return math.floor(lo / m + 1e-9) * m, math.ceil(hi / m - 1e-9) * m


def _ticks(lo, hi, step=None):
    step = step or _step(lo, hi)
    t, out = math.ceil(lo / step - 1e-9) * step, []
    while t <= hi + 1e-9 * step:
        out.append(0.0 if abs(t) < 1e-9 * step else t)
        t += step
    return out


def _fmt(v, step=None):
    """A tick value with as many decimals as the step needs (0, 0.5, 1.0 …), a true minus sign."""
    if step is None:
        return f"{v:.4g}".replace("-", "\u2212")
    a = abs(v)
    if a and (a >= 1e5 or a < 1e-3):
        return f"{v:.1e}".replace("-", "\u2212")
    dec = 0
    while dec < 6 and abs(step * 10 ** dec - round(step * 10 ** dec)) > 1e-6 * 10 ** dec:
        dec += 1
    return f"{v:.{dec}f}".replace("-", "\u2212")


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
