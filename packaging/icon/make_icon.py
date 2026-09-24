#!/usr/bin/env python3
"""Draws the CAPS icon (the Studio's hexagon mark on a dark rounded square) at 1024 px, then writes caps.png,
caps.ico (Windows) and, on macOS, caps.icns (iconutil). Needs Pillow."""
import math, os, shutil, subprocess, sys
from PIL import Image, ImageDraw

HERE = os.path.dirname(os.path.abspath(__file__))
S = 1024
BG, ACC, TXT = (22, 25, 28, 255), (240, 168, 60, 255), (232, 230, 225, 255)


def draw(size):
    k = 4   # supersample
    n = size * k
    im = Image.new("RGBA", (n, n), (0, 0, 0, 0))
    d = ImageDraw.Draw(im)
    pad = int(n * 0.08)
    d.rounded_rectangle([pad, pad, n - pad, n - pad], radius=int(n * 0.2), fill=BG)
    c, r = n / 2, n * 0.3
    def hexagon(rad):
        return [(c + rad * math.cos(math.radians(90 + 60 * i)), c - rad * math.sin(math.radians(90 + 60 * i))) for i in range(6)]
    w = n * 0.045   # ring: outer hexagon filled, inner one cut out (clean corners)
    d.polygon(hexagon(r + w / 2 / math.cos(math.radians(30))), fill=ACC)
    d.polygon(hexagon(r - w / 2 / math.cos(math.radians(30))), fill=BG)
    rc = n * 0.085
    d.ellipse([c - rc, c - rc, c + rc, c + rc], fill=ACC)
    rs = n * 0.05
    for ang in (90, 210, 330):   # three atoms around the centre
        x, y = c + r * 0.62 * math.cos(math.radians(ang)), c - r * 0.62 * math.sin(math.radians(ang))
        d.ellipse([x - rs, y - rs, x + rs, y + rs], fill=TXT)
    return im.resize((size, size), Image.LANCZOS)


big = draw(S)
big.save(os.path.join(HERE, "caps.png"))
big.save(os.path.join(HERE, "caps.ico"), sizes=[(16, 16), (24, 24), (32, 32), (48, 48), (64, 64), (128, 128), (256, 256)])
for s in (16, 32, 48, 64, 128, 256, 512):
    draw(s).save(os.path.join(HERE, f"caps_{s}.png"))
if sys.platform == "darwin" and shutil.which("iconutil"):
    iset = os.path.join(HERE, "caps.iconset")
    os.makedirs(iset, exist_ok=True)
    for s in (16, 32, 128, 256, 512):
        draw(s).save(os.path.join(iset, f"icon_{s}x{s}.png"))
        draw(2 * s).save(os.path.join(iset, f"icon_{s}x{s}@2x.png"))
    subprocess.run(["iconutil", "-c", "icns", iset, "-o", os.path.join(HERE, "caps.icns")], check=True)
    shutil.rmtree(iset)
print("wrote", ", ".join(f for f in sorted(os.listdir(HERE)) if f.startswith("caps")))
