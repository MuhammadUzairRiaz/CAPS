"""caps.View: the current frame of a Document in a notebook — drag to rotate, wheel to zoom, double-click to reset.

The page gets the atoms once (positions, elements, bonds, the cell) and draws them itself in a <canvas>: no widget
extension, no server round trip, so it works in Jupyter, JupyterLab, VS Code and exported HTML. Where HTML is not
shown (a terminal, GitHub's preview of a saved notebook without scripts) the PNG from CAPS's renderer stands in.
"""
from __future__ import annotations

import json
import os
import tempfile
import uuid
from typing import Optional

_STYLES = {"ball-and-stick", "space-filling", "sticks", "no-hydrogens", "lines"}

_JS = r"""
(function(){
const D = %(data)s, W = %(w)d, H = %(h)d, root = document.getElementById("%(id)s");
const cv = root.querySelector("canvas"), g = cv.getContext("2d"), dpr = window.devicePixelRatio || 1;
cv.width = W * dpr; cv.height = H * dpr; cv.style.width = W + "px"; cv.style.height = H + "px";
const n = D.z.length, P = D.xyz, B = D.bonds, style = D.style;
let cx = 0, cy = 0, cz = 0;
for (let i = 0; i < n; i++) { cx += P[3*i]; cy += P[3*i+1]; cz += P[3*i+2]; }
cx /= Math.max(n, 1); cy /= Math.max(n, 1); cz /= Math.max(n, 1);
let ext = 1;
for (let i = 0; i < n; i++) ext = Math.max(ext, Math.hypot(P[3*i]-cx, P[3*i+1]-cy, P[3*i+2]-cz));
const scale0 = 0.46 * Math.min(W, H) / (ext + 2);
let R = [1,0,0, 0,1,0, 0,0,1], zoom = 1;
function rot(ax, ay) {   // turn about the screen axes
  const c1 = Math.cos(ay), s1 = Math.sin(ay), c2 = Math.cos(ax), s2 = Math.sin(ax);
  const Y = [c1,0,s1, 0,1,0, -s1,0,c1], X = [1,0,0, 0,c2,-s2, 0,s2,c2];
  const m = (a, b) => { const o = new Array(9); for (let r = 0; r < 3; r++) for (let c = 0; c < 3; c++) o[3*r+c] = a[3*r]*b[c] + a[3*r+1]*b[3+c] + a[3*r+2]*b[6+c]; return o; };
  R = m(X, m(Y, R));
}
rot(-0.40, 0.55);
const sprites = {};
function sprite(z) {
  if (sprites[z]) return sprites[z];
  const s = document.createElement("canvas"), k = 96; s.width = s.height = k;
  const c = s.getContext("2d"), col = D.colours[z] || "#ff1493";
  const gr = c.createRadialGradient(k*0.36, k*0.34, k*0.04, k*0.5, k*0.5, k*0.5);
  gr.addColorStop(0, "#ffffff"); gr.addColorStop(0.22, col); gr.addColorStop(1, shade(col, 0.45));
  c.fillStyle = gr; c.beginPath(); c.arc(k/2, k/2, k/2 - 2, 0, 2*Math.PI); c.fill();
  c.lineWidth = 2.5; c.strokeStyle = "rgba(0,0,0,0.55)"; c.stroke();
  return sprites[z] = s;
}
function shade(hex, f) {
  const v = parseInt(hex.slice(1), 16);
  const r = Math.round(((v >> 16) & 255) * f), gg = Math.round(((v >> 8) & 255) * f), b = Math.round((v & 255) * f);
  return "rgb(" + r + "," + gg + "," + b + ")";
}
function radius(z) {
  const r = D.radii[z] || 1.7;
  return style === "space-filling" ? r : style === "sticks" || style === "lines" ? 0.18 : 0.26 * r;
}
const X = new Float32Array(n), Yp = new Float32Array(n), Z = new Float32Array(n);
function draw() {
  const s = scale0 * zoom;
  for (let i = 0; i < n; i++) {
    const x = P[3*i]-cx, y = P[3*i+1]-cy, z = P[3*i+2]-cz;
    X[i] = W/2 + s * (R[0]*x + R[1]*y + R[2]*z);
    Yp[i] = H/2 - s * (R[3]*x + R[4]*y + R[5]*z);
    Z[i] = R[6]*x + R[7]*y + R[8]*z;
  }
  g.setTransform(dpr, 0, 0, dpr, 0, 0);
  g.fillStyle = D.background; g.fillRect(0, 0, W, H);
  if (D.cell) {   // the cell's twelve edges
    const C = D.cell, o = [C[0]-cx, C[1]-cy, C[2]-cz], a = C.slice(3, 6), b = C.slice(6, 9), c = C.slice(9, 12);
    const corner = (i, j, k) => { const p = [o[0]+i*a[0]+j*b[0]+k*c[0], o[1]+i*a[1]+j*b[1]+k*c[1], o[2]+i*a[2]+j*b[2]+k*c[2]];
      return [W/2 + s*(R[0]*p[0]+R[1]*p[1]+R[2]*p[2]), H/2 - s*(R[3]*p[0]+R[4]*p[1]+R[5]*p[2])]; };
    g.strokeStyle = "rgba(120,120,120,0.7)"; g.lineWidth = 1; g.beginPath();
    for (const [p, q] of [[[0,0,0],[1,0,0]],[[0,0,0],[0,1,0]],[[0,0,0],[0,0,1]],[[1,1,0],[0,1,0]],[[1,1,0],[1,0,0]],[[1,1,0],[1,1,1]],
                          [[1,0,1],[0,0,1]],[[1,0,1],[1,1,1]],[[1,0,1],[1,0,0]],[[0,1,1],[0,1,0]],[[0,1,1],[0,0,1]],[[0,1,1],[1,1,1]]]) {
      const A = corner(...p), Bq = corner(...q); g.moveTo(A[0], A[1]); g.lineTo(Bq[0], Bq[1]);
    }
    g.stroke();
  }
  // bonds and atoms back to front: each bond half is drawn with the atom it starts from
  const order = Array.from({length: n}, (_, i) => i).sort((i, j) => Z[i] - Z[j]);
  const nb = B.length / 2, half = {};
  if (style !== "space-filling")
    for (let k = 0; k < nb; k++) { const i = B[2*k], j = B[2*k+1];
      if (Math.hypot(P[3*i]-P[3*j], P[3*i+1]-P[3*j+1], P[3*i+2]-P[3*j+2]) > 3.2) continue;   // bonds across the cell
      (half[i] = half[i] || []).push(j); (half[j] = half[j] || []).push(i); }
  const bw = style === "lines" ? 1.6 : Math.max(1.2, (style === "sticks" ? 0.30 : 0.16) * s);
  g.lineCap = "round";
  for (const i of order) {
    const zi = D.z[i];
    if (half[i]) {
      g.lineWidth = bw + 1.2; g.strokeStyle = "rgba(0,0,0,0.5)";
      for (const j of half[i]) { g.beginPath(); g.moveTo(X[i], Yp[i]); g.lineTo((X[i]+X[j])/2, (Yp[i]+Yp[j])/2); g.stroke(); }
      g.lineWidth = bw; g.strokeStyle = D.colours[zi] || "#ff1493";
      for (const j of half[i]) { g.beginPath(); g.moveTo(X[i], Yp[i]); g.lineTo((X[i]+X[j])/2, (Yp[i]+Yp[j])/2); g.stroke(); }
    }
    if (style !== "lines") { const r = Math.max(1, radius(zi) * s); g.drawImage(sprite(zi), X[i]-r, Yp[i]-r, 2*r, 2*r); }
  }
}
let drag = null;
cv.addEventListener("pointerdown", e => { drag = [e.clientX, e.clientY]; cv.setPointerCapture(e.pointerId); });
cv.addEventListener("pointermove", e => { if (!drag) return; rot(-(e.clientY-drag[1])*0.01, (e.clientX-drag[0])*0.01); drag = [e.clientX, e.clientY]; draw(); });
cv.addEventListener("pointerup", () => drag = null);
cv.addEventListener("wheel", e => { e.preventDefault(); zoom *= Math.exp(-e.deltaY * 0.0015); zoom = Math.min(20, Math.max(0.1, zoom)); draw(); }, {passive: false});
cv.addEventListener("dblclick", () => { R = [1,0,0, 0,1,0, 0,0,1]; rot(-0.40, 0.55); zoom = 1; draw(); });
draw();
})();
"""


class View:
    """What Document.view returns; shown by the notebook (HTML with a canvas, or PNG)."""

    def __init__(self, doc, style: str = "ball-and-stick", width: int = 640, height: int = 400, hydrogens: bool = True,
                 max_atoms: int = 60000, background: str = "#ffffff", cell: Optional[bool] = None):
        style = style.replace("_", "-")
        if style not in _STYLES:
            raise ValueError("style: " + ", ".join(sorted(_STYLES)))
        self.doc, self.style, self.width, self.height, self.background = doc, style, width, height, background
        self.hydrogens = hydrogens and style != "no-hydrogens"
        self.max_atoms = max_atoms
        self.cell = cell   # None: drawn when there is more than one molecule

    def _show_cell(self) -> bool:
        return self.cell if self.cell is not None else self.doc.summary()["molecules"] > 1

    def scene(self) -> dict:
        return self.doc.scene(hydrogens=self.hydrogens, max_atoms=self.max_atoms)

    def __repr__(self) -> str:
        s = self.scene()
        return f"<caps.View · {s['shown']} atoms shown · {self.style}>"

    def _repr_html_(self) -> str:
        s = self.scene()
        s["style"] = self.style
        s["background"] = self.background
        if not self._show_cell():
            s["cell"] = None
        vid = "caps-view-" + uuid.uuid4().hex[:10]
        caption = f"drag to rotate · CAPS View widget · {s['shown']:,} atoms shown".replace(",", " ")
        if s["shown"] < s["atoms"]:
            caption += f" of {s['atoms']:,}".replace(",", " ")
        js = _JS % {"data": json.dumps(s, separators=(",", ":")), "w": self.width, "h": self.height, "id": vid}
        return (f'<div id="{vid}" style="display:inline-block;border:1px solid #d8d8d8;border-radius:6px;overflow:hidden">'
                f'<canvas style="display:block;cursor:grab;touch-action:none"></canvas>'
                f'<div style="font:12px system-ui,sans-serif;color:#666;padding:4px 8px;border-top:1px solid #e6e6e6">{caption}</div>'
                f'</div><script>{js}</script>')

    def _repr_png_(self) -> bytes:
        fd, path = tempfile.mkstemp(suffix=".png")
        os.close(fd)
        try:
            self.doc.render(path, width=self.width, height=self.height, style={"ball-and-stick": "ball_and_stick", "space-filling": "space_filling",
                                                                                "sticks": "sticks", "no-hydrogens": "no_hydrogens", "lines": "sticks"}[self.style])
            with open(path, "rb") as f:
                return f.read()
        finally:
            os.remove(path)

    def save_html(self, path: str) -> None:
        """A standalone page with the view."""
        with open(path, "w", encoding="utf-8") as f:
            f.write("<!doctype html><meta charset='utf-8'><title>CAPS view</title><body style='margin:16px'>" + self._repr_html_() + "</body>")
