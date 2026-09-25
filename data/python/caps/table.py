"""caps.table: results side by side — one row per result (a sweep condition, or a Document with .properties), one
column per property, mean ± standard deviation over seeds. Shown as HTML in notebooks; .to_csv(), .rows for code."""
from __future__ import annotations

import html
import math

from .core import _PROPERTY_ALIASES

_HEAD = {"density": ("density", "g/cm³"), "rg": ("Rg", "Å"), "cn": ("c_inf", ""), "tg": ("tg", "K"), "ree": ("⟨R²⟩½", "Å")}


def _fmt(m: float, sd: float, prop: str) -> str:
    if not math.isfinite(m):
        return "—"
    digits = 3 if prop == "density" else 0 if prop == "tg" else 2
    s = f"{m:.{digits}f}"
    return s + (f" ± {sd:.{digits}f}" if math.isfinite(sd) else "")


def _stat(cell, prop: str):
    pid = _PROPERTY_ALIASES.get(prop, prop)
    if hasattr(cell, "runs"):   # a sweep result
        return cell[pid]
    for p in getattr(cell, "properties", []) or []:   # a Document from caps.run / caps.polymer
        if p.get("id") == pid and p.get("value") is not None:
            e = p.get("error")
            return (p["value"], e if e is not None else math.nan, 1)
    return (math.nan, math.nan, 0)


class Table:
    def __init__(self, cells, properties):
        self.cells = list(cells)
        self.properties = [properties] if isinstance(properties, str) else list(properties)

    def _label(self, c) -> str:
        return getattr(c, "label", None) or str(c)

    def _unit(self, prop: str) -> str:
        pid = _PROPERTY_ALIASES.get(prop, prop)
        for c in self.cells:
            if hasattr(c, "unit"):
                u = c.unit(pid)
                if u:
                    return u
        return _HEAD.get(pid, ("", ""))[1]

    @property
    def rows(self) -> list:
        """[(label, {prop: (mean, sd, n)})]"""
        return [(self._label(c), {p: _stat(c, p) for p in self.properties}) for c in self.cells]

    def _header(self, p: str) -> str:
        u = self._unit(p)
        return f"{p} ({u})" if u else p

    def __repr__(self) -> str:
        head = [""] + [self._header(p) for p in self.properties]
        body = [[lbl] + [_fmt(*v[p][:2], _PROPERTY_ALIASES.get(p, p)) for p in self.properties] for lbl, v in self.rows]
        w = [max(len(r[i]) for r in [head] + body) for i in range(len(head))]
        line = lambda r: "  ".join(x.ljust(w[i]) if i == 0 else x.rjust(w[i]) for i, x in enumerate(r))
        return "\n".join([line(head)] + [line(r) for r in body])

    def _repr_html_(self) -> str:
        th = "".join(f"<th style='text-align:right;padding:4px 12px'>{html.escape(self._header(p))}</th>" for p in self.properties)
        trs = ""
        for lbl, v in self.rows:
            tds = "".join(f"<td style='text-align:right;padding:4px 12px;font-variant-numeric:tabular-nums'>{html.escape(_fmt(*v[p][:2], _PROPERTY_ALIASES.get(p, p)))}</td>"
                          for p in self.properties)
            trs += f"<tr><td style='padding:4px 12px'>{html.escape(lbl)}</td>{tds}</tr>"
        return f"<table style='border-collapse:collapse;font-size:13px'><thead><tr><th></th>{th}</tr></thead><tbody>{trs}</tbody></table>"

    def to_csv(self, path: str = "") -> str:
        cols = ["condition"] + [c for p in self.properties for c in (p, p + "_sd", p + "_n")]
        lines = [",".join(cols)]
        for lbl, v in self.rows:
            lines.append(",".join([lbl] + [f"{x:.10g}" if isinstance(x, float) and math.isfinite(x) else ("" if isinstance(x, float) else str(x))
                                           for p in self.properties for x in v[p]]))
        text = "\n".join(lines) + "\n"
        if path:
            with open(path, "w", encoding="utf-8") as f:
                f.write(text)
        return text


def table(cells, properties) -> Table:
    return Table(cells, properties)
