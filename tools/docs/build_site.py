#!/usr/bin/env python3
"""The CAPS documentation site (GitHub Pages, served from docs/): generated, so it never drifts from the program.

    python3 tools/docs/build_site.py [--caps build/cli/caps] [--lib build/capi/libcaps.dylib] [--out docs]

Sources
  tools/docs/src/*.md          the guide and the tutorials (a small Markdown: headings, paragraphs, lists, tables,
                               fenced code, > notes, $$ display maths $$, `code`, **bold**, *italic*, [links](url))
  caps (no arguments)          the usage text: one entry per command
  caps COMMAND --help          the options, defaults and reasons of the commands that have them
  data/manual/manual.json      the theory pages: equations (LaTeX, typeset by KaTeX), symbols, when to use, tests
  libcaps                      the references (caps_citation_text) and the Python package's docstrings
Writes docs/index.html, docs/cli/, docs/tutorials/, docs/theory/, docs/python.html, docs/assets/ and docs/.nojekyll;
docs/manual/ (the Studio tour) is left as it is.
"""
from __future__ import annotations

import argparse
import ctypes
import html
import inspect
import json
import os
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = Path(__file__).resolve().parent / "src"
KATEX = "https://cdn.jsdelivr.net/npm/katex@0.16.11/dist"
REPO = "https://github.com/MuhammadUzairRiaz/CAPS"

# ---------------------------------------------------------------------------------------------------------- commands
# (group, blurb, commands) in the order a project uses them
GROUPS = [
    ("Build structures", "Molecules, chains, cells, crystals, surfaces and fillers.",
     ["build", "grow", "blend", "peptide", "crystal", "surface", "interface", "nano", "pore", "solvate", "pack", "edit"]),
    ("Force fields", "Types, charges and parameters; exports for LAMMPS, GROMACS and DL_POLY.", ["ff", "field"]),
    ("Relax and run", "Minimisation, molecular dynamics, equilibration protocols, regrowth, crosslinking and recipes.",
     ["relax", "md", "equilibrate", "cbmc", "react", "run"]),
    ("Properties", "Mechanics, glass transition, elastic constants and interfaces.", ["tensile", "tg", "elastic", "pull"]),
    ("Inspect and analyse", "Files, structure, trajectories and figures.",
     ["info", "check", "contacts", "chains", "shape", "rdf", "analyze", "frames", "convert", "render", "pipeline"]),
    ("Reproducibility", "Provenance, bundles and the validation suite.", ["provenance", "bundle", "reproduce", "bench"]),
    ("Coarse-graining", "Chemistry-aware mapping, bonded and non-bonded potentials, melts, entanglements, backmapping.",
     ["cgmap", "cgfit", "cgbuild", "ppa", "mech", "cgdyn", "backmap"]),
    ("Biomolecular coarse-graining", "Secondary structure and Martini proteins.", ["dssp", "martini"]),
    ("DFT surfaces", "2D sheets, terminations, adsorption complexes and VASP case folders.",
     ["sheet", "terminate", "validate", "adsorb-dft", "vasp-set", "vasp-conv", "vasp-scan", "vasp-derived", "vasp-jobs",
      "vasp-check", "vasp-progress", "vasp-health", "vasp-bind", "vasp-analyze"]),
]
# one line per command for the guide's tables (the reference shows the program's own text)
SUMMARY = {
    "build": "a 3D molecule from SMILES, with conformers and a force-field clean-up",
    "grow": "polymer chains of any repeat units grown in a periodic cell (tacticity, copolymer sequences)",
    "blend": "two or more polymers in one cell: mixed, as slabs or as a droplet",
    "peptide": "an all-atom peptide from a sequence: helix, strand or your secondary structure, at a pH",
    "crystal": "a crystal from a space group and sites, or a CIF's symmetry found",
    "surface": "a slab cut from a crystal at (hkl), terminations listed",
    "interface": "a polymer film grown on a crystal surface",
    "nano": "nanotubes, graphene sheets and nanoparticles, alone or in a polymer matrix",
    "pore": "a fluid in a slit, cylindrical or framework pore",
    "solvate": "solvent and ions packed around a solute",
    "pack": "molecules packed by a packmol-style input",
    "edit": "scripted edits: elements, bonds, hydrogens, geometry, stereo, supercells, vacuum",
    "ff": "force fields: import, inspect, type atoms by SMARTS rules, apply parameters and export",
    "field": "the default force field's types, terms and energy for a structure",
    "relax": "energy minimisation (L-BFGS, CG, FIRE) with push-off and box relaxation",
    "md": "molecular dynamics: NVE, NVT, NPT, constraints, checkpoints",
    "equilibrate": "equilibration protocols: Larsen's 21 steps, annealing cycles, push-off, until converged",
    "cbmc": "configurational-bias Monte Carlo regrowth of chain ends",
    "react": "crosslinking and curing by reaction templates, cycle by cycle",
    "run": "a recipe (build to export in one file) or a saved pipeline over many inputs",
    "tensile": "uniaxial tension: stress–strain curve, modulus, yield",
    "tg": "glass transition from a cooling scan of density",
    "elastic": "elastic constants by static strain or stress fluctuations",
    "pull": "pull-out of a film from a surface: interfacial shear strength, work of separation",
    "info": "what a file holds: atoms, molecules, bonds, cell, frames",
    "check": "file checks: counts, bonds, contacts, charges, cell",
    "contacts": "the closest distance between atoms of different molecules",
    "chains": "backbones and mean-square internal distances",
    "shape": "per-molecule radius of gyration and shape",
    "rdf": "radial distribution function of one structure",
    "analyze": "properties over a trajectory: density, g(r), S(q), Rg, C∞, MSD, CED, free volume, interfaces",
    "frames": "a trajectory thinned or converted",
    "convert": "a structure in another format",
    "render": "a picture (PNG or SVG) of a structure without the Studio",
    "pipeline": "an analysis pipeline on one frame: attributes and tables",
    "provenance": "how a file was made: steps, methods paragraph, BibTeX, comparison of two runs",
    "bundle": "a figure's data with its pipeline, provenance and hashes",
    "reproduce": "a bundle rebuilt from its input, checked by SHA-256",
    "bench": "the built-in validation suite against reference values",
    "cgmap": "chemistry-aware mapping to beads by SMARTS rules; trajectories mapped frame by frame",
    "cgfit": "CG potentials: bonded Boltzmann inversion, per-pair IBI, fits, T_g, σ/ε calibration",
    "cgbuild": "a CG melt of repeat units with real sequence statistics",
    "ppa": "entanglements: primitive paths and the N_e estimators",
    "mech": "tension decks for CG melts and their analysis (modulus, yield, strain hardening)",
    "cgdyn": "chain dynamics: g₁, g₂, g₃, D, τ_R, τ_e; CG-to-AA time mapping",
    "backmap": "CG beads back to all atoms, fragment by fragment, with a relaxation deck",
    "dssp": "DSSP secondary structure of a protein",
    "martini": "a Martini 2.2 or 3 protein model",
    "sheet": "a 2D sheet from a preset or cut from a bulk structure",
    "terminate": "terminations (O, OH, F …) on both faces of a sheet",
    "validate": "PASS/FAIL checks of a 2D slab",
    "adsorb-dft": "a molecule on a slab in several placements, for DFT",
    "vasp-set": "one VASP case folder: POSCAR, INCARs with reasons, KPOINTS, job script",
    "vasp-conv": "ENCUT and k-mesh convergence",
    "vasp-scan": "scans: make the cases, fit the curve",
    "vasp-derived": "follow-up runs: charges, charge-density difference, frequencies, AIMD",
    "vasp-jobs": "submit, update, reset, clean up and store many VASP jobs",
    "vasp-check": "check a case before or after running",
    "vasp-progress": "progress of running cases",
    "vasp-health": "stuck or failing runs",
    "vasp-bind": "adsorption energies from slab, molecule and complex",
    "vasp-analyze": "geometry, work function, DOS, Δρ, Bader, frequencies, AIMD, the summary table",
}
# words the public pages leave out (CAPS's own names for what it does)
SCRUB = [re.compile(r"^\s*caps ff import-dlf.*$", re.M)]


def run(cmd: list[str]) -> str:
    r = subprocess.run(cmd, capture_output=True, text=True, env=dict(os.environ, CAPS_HOME=str(ROOT)))
    return (r.stdout or "") + (r.stderr or "")


def usage_blocks(caps: str) -> dict[str, list[str]]:
    """The usage text cut into entries: a line '  caps NAME …' and its indented continuation lines."""
    text = run([caps])
    for rx in SCRUB:
        text = rx.sub("", text)
    blocks: dict[str, list[str]] = {}
    cur = None
    for line in text.splitlines():
        m = re.match(r"^  caps ([a-z][a-z0-9-]*)", line)
        if m:
            cur = m.group(1)
            blocks.setdefault(cur, []).append(line[2:])
        elif cur and line.startswith("    "):
            blocks[cur].append(line[2:])
        else:
            cur = None
    # "vasp-conv … · vasp-scan …" lines name several commands
    for name, lines in list(blocks.items()):
        for other in re.findall(r"· (vasp-[a-z]+)", " ".join(lines)):
            blocks.setdefault(other, lines)
    return blocks


def command_help(caps: str, name: str) -> str:
    out = run([caps, name, "--help"])
    return out if out.startswith("usage:") else ""


# ---------------------------------------------------------------------------------------------------------- markdown
def inline(s: str) -> str:
    """`code`, **bold**, *italic*, [text](url) and \\( maths \\) (left for KaTeX) on escaped text."""
    links: list[str] = []

    def keep(m):
        links.append(f'<a href="{html.escape(m.group(2))}">{inline(m.group(1))}</a>')
        return f"\x00{len(links) - 1}\x00"
    s = re.sub(r"\[((?:[^\]`]|`[^`]*`)+)\]\(([^)\s]+)\)", keep, s)
    out = _inline(s)
    return re.sub(r"\x00(\d+)\x00", lambda m: links[int(m.group(1))], out)


def _inline(s: str) -> str:
    parts = re.split(r"(`[^`]+`)", s)
    out = []
    for p in parts:
        if p.startswith("`") and p.endswith("`") and len(p) > 1:
            out.append("<code>" + html.escape(p[1:-1]) + "</code>")
            continue
        p = html.escape(p, quote=False)
        p = re.sub(r"\*\*(.+?)\*\*", r"<strong>\1</strong>", p)
        p = re.sub(r"(?<![\w*])\*(?!\s)(.+?)(?<!\s)\*(?![\w*])", r"<em>\1</em>", p)
        p = re.sub(r"\[([^\]]+)\]\(([^)\s]+)\)", lambda m: f'<a href="{m.group(2)}">{m.group(1)}</a>', p)
        out.append(p)
    return "".join(out)


def slug(s: str) -> str:
    s = re.sub(r"<[^>]+>", "", s)
    s = re.sub(r"[^a-z0-9]+", "-", s.lower()).strip("-")
    return s or "section"


def markdown(text: str) -> tuple[str, list[tuple[int, str, str]], dict]:
    """HTML, the headings (level, id, text) and the front matter (--- key: value --- at the top)."""
    meta: dict = {}
    lines = text.splitlines()
    if lines and lines[0].strip() == "---":
        end = lines.index("---", 1)
        for l in lines[1:end]:
            k, _, v = l.partition(":")
            meta[k.strip()] = v.strip()
        lines = lines[end + 1:]
    out, toc, i = [], [], 0
    used: set[str] = set()

    def para(buf):
        if buf:
            out.append("<p>" + inline(" ".join(x.strip() for x in buf)) + "</p>")
            buf.clear()

    buf: list[str] = []
    while i < len(lines):
        l = lines[i]
        if l.startswith("```"):
            para(buf)
            lang = l[3:].strip()
            j = i + 1
            code = []
            while j < len(lines) and not lines[j].startswith("```"):
                code.append(lines[j])
                j += 1
            label = {"bash": "shell", "sh": "shell", "python": "Python", "yaml": "recipe (YAML)", "text": "output"}.get(lang, lang)
            out.append(f'<div class="code" data-lang="{html.escape(label)}"><pre><code>' + html.escape("\n".join(code)) + "</code></pre></div>")
            i = j + 1
            continue
        if l.strip() == "$$":
            para(buf)
            j = i + 1
            tex = []
            while j < len(lines) and lines[j].strip() != "$$":
                tex.append(lines[j])
                j += 1
            out.append('<div class="eq">$$' + html.escape("\n".join(tex), quote=False) + "$$</div>")
            i = j + 1
            continue
        m = re.match(r"^(#{1,4}) (.+)$", l)
        if m:
            para(buf)
            level, title = len(m.group(1)), m.group(2).strip()
            hid = slug(title)
            while hid in used:
                hid += "-x"
            used.add(hid)
            if level >= 2:
                toc.append((level, hid, title))
            out.append(f'<h{level} id="{hid}">' + inline(title) + f"</h{level}>")
            i += 1
            continue
        if l.startswith("> "):
            para(buf)
            note = []
            while i < len(lines) and lines[i].startswith(">"):
                note.append(lines[i][1:].strip())
                i += 1
            kind = "note"
            if note and re.match(r"^\*\*(Warning|Careful)", note[0]):
                kind = "warn"
            out.append(f'<aside class="{kind}">' + inline(" ".join(note)) + "</aside>")
            continue
        if l.startswith("|"):
            para(buf)
            rows = []
            while i < len(lines) and lines[i].startswith("|"):
                cells = re.split(r"(?<!\\)\|", lines[i].strip()[1:-1] if lines[i].strip().endswith("|") else lines[i].strip()[1:])
                rows.append([c.strip().replace("\\|", "|") for c in cells])
                i += 1
            head, body = rows[0], [r for r in rows[1:] if not all(re.match(r"^:?-+:?$", c) for c in r)]
            t = ['<div class="table"><table><thead><tr>' + "".join("<th>" + inline(c) + "</th>" for c in head) + "</tr></thead><tbody>"]
            for r in body:
                t.append("<tr>" + "".join("<td>" + inline(c) + "</td>" for c in r) + "</tr>")
            t.append("</tbody></table></div>")
            out.append("".join(t))
            continue
        m = re.match(r"^(\s*)([-*]|\d+\.) (.+)$", l)
        if m:
            para(buf)
            ordered = m.group(2)[0].isdigit()
            items = []
            while i < len(lines):
                mm = re.match(r"^(\s*)([-*]|\d+\.) (.+)$", lines[i])
                if mm:
                    items.append(mm.group(3))
                elif lines[i].startswith("   ") and items:
                    items[-1] += " " + lines[i].strip()
                else:
                    break
                i += 1
            tag = "ol" if ordered else "ul"
            out.append(f"<{tag}>" + "".join("<li>" + inline(x) + "</li>" for x in items) + f"</{tag}>")
            continue
        if not l.strip():
            para(buf)
        else:
            buf.append(l)
        i += 1
    para(buf)
    return "\n".join(out), toc, meta


# ---------------------------------------------------------------------------------------------------------- page shell
NAV = [("index.html", "Home"), ("cli/index.html", "CLI guide"), ("cli/reference.html", "Command reference"),
       ("tutorials/index.html", "Tutorials"), ("theory/index.html", "Theory"), ("python.html", "Python"),
       ("manual/index.html", "Studio tour")]


def page(path: str, title: str, body: str, toc_html: str = "", maths: bool = False, wide: bool = False, desc: str = "") -> str:
    depth = path.count("/")
    up = "../" * depth
    section = {"cli/reference.html": "cli/reference.html"}.get(path, path.split("/")[0] + "/index.html" if "/" in path else path)
    nav = "".join(f'<a href="{up}{href}"' + (' aria-current="page"' if href == section else "") + f">{name}</a>" for href, name in NAV)
    head_maths = (f'<link rel="stylesheet" href="{KATEX}/katex.min.css">'
                  f'<script defer src="{KATEX}/katex.min.js"></script>'
                  f'<script defer src="{KATEX}/contrib/auto-render.min.js"></script>') if maths else ""
    layout = "doc wide" if wide else "doc"
    side = f'<nav class="toc" aria-label="On this page">{toc_html}</nav>' if toc_html else ""
    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">
<title>{html.escape(title)}</title>
<meta name="description" content="{html.escape(desc or 'CAPS: a polymer and materials simulation workbench — command line, Python and Studio.')}">
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Archivo:wght@500;600;700&family=IBM+Plex+Sans:ital,wght@0,400;0,500;0,600;1,400&family=IBM+Plex+Mono:wght@400;500&display=swap">
<link rel="stylesheet" href="{up}assets/site.css">
{head_maths}
<script defer src="{up}assets/site.js"></script>
</head>
<body>
<header class="top">
  <a class="brand" href="{up}index.html"><span class="mark" aria-hidden="true"></span>CAPS</a>
  <button class="menu" aria-expanded="false" aria-controls="nav">Menu</button>
  <nav id="nav" class="links">{nav}<a href="{REPO}">GitHub</a></nav>
</header>
<div class="{layout}">
{side}
<main>
{body}
</main>
</div>
<footer class="foot"><span>CAPS 0.1.0 · generated by <code>tools/docs/build_site.py</code> from the program itself</span><a href="{REPO}">{REPO.replace('https://', '')}</a></footer>
</body>
</html>
"""


def toc_from(headings, title="On this page") -> str:
    items = "".join(f'<a class="l{lvl}" href="#{hid}">{inline(t)}</a>' for lvl, hid, t in headings if lvl <= 3)
    return f"<h2>{title}</h2>{items}" if items else ""


# ---------------------------------------------------------------------------------------------------------- pages
def build_reference(caps: str, out: Path) -> dict[str, str]:
    blocks = usage_blocks(caps)
    listed = {c for _, _, cs in GROUPS for c in cs}
    missing = sorted(set(blocks) - listed)
    groups = GROUPS + ([("Other", "", missing)] if missing else [])
    summary: dict[str, str] = {}
    body = ['<header class="hero"><p class="eyebrow">Command line · every command</p><h1>Command reference</h1>'
            '<p>Every <code>caps</code> command with its synopsis, as the program prints it. Commands marked '
            '<span class="pill">--help</span> print each option with its default and the reason for it; that text is '
            'reproduced in full. Shared conventions (units, files, seeds, exit codes) are in the '
            '<a href="index.html">CLI guide</a>.</p>'
            '<label class="filter"><span>Filter</span><input id="cmdfilter" type="search" placeholder="grow, rdf, tensile …" autocomplete="off"></label></header>']
    toc = []
    for g, blurb, cs in groups:
        gid = slug(g)
        toc.append(f'<h2>{html.escape(g)}</h2>' + "".join(f'<a href="#{c}"><code>{c}</code></a>' for c in cs if c in blocks or command_help(caps, c)))
        body.append(f'<section class="group" id="{gid}"><h2>{html.escape(g)}</h2><p class="lede">{html.escape(blurb)}</p>')
        for c in cs:
            usage = "\n".join(blocks.get(c, []))
            hp = command_help(caps, c)
            if not usage and not hp:
                continue
            lead = ""
            if hp:
                paras = hp.split("\n\n")
                lead = paras[1].strip() if len(paras) > 1 and not paras[1].startswith("options") else ""
            else:
                m = re.search(r"\s{3,}([a-z][^\n]*)$", blocks[c][0])
                lead = m.group(1) if m else ""
            summary[c] = (lead.split(". ")[0].split(": ")[0] if lead else "").strip().rstrip(".")
            body.append(f'<article class="cmd" id="{c}" data-name="{c}"><h3><code>caps {c}</code>'
                        + ('<span class="pill">--help</span>' if hp else "") + "</h3>")
            if lead and not hp:   # the help text below starts with the same description
                body.append("<p>" + html.escape(lead) + "</p>")
            text = hp if hp else usage
            body.append('<div class="code" data-lang="' + ("caps " + c + " --help" if hp else "synopsis") + '"><pre><code>'
                        + html.escape(text.rstrip()) + "</code></pre></div></article>")
        body.append("</section>")
    (out / "cli").mkdir(parents=True, exist_ok=True)
    (out / "cli" / "reference.html").write_text(page("cli/reference.html", "CAPS command reference", "\n".join(body), "".join(toc),
                                                     desc="Every caps command with its options, defaults and reasons."))
    return summary


def build_guide(out: Path, summary: dict[str, str]) -> None:
    md = (SRC / "cli-guide.md").read_text()
    # the generated "which command" tables
    rows = []
    for g, blurb, cs in GROUPS:
        rows.append(f"### {g}\n\n{blurb}\n\n| Command | What it does |\n|---|---|")
        for c in cs:
            s = SUMMARY.get(c) or summary.get(c, "")
            rows.append(f"| [`caps {c}`](reference.html#{c}) | {s.replace('|', '/')} |")
        rows.append("")
    md = md.replace("<!-- COMMAND TABLES -->", "\n".join(rows))
    body, toc, meta = markdown(md)
    hero = (f'<header class="hero"><p class="eyebrow">{html.escape(meta.get("eyebrow", "Command line"))}</p>'
            f'<h1>{html.escape(meta.get("title", "CLI guide"))}</h1><p>{inline(meta.get("lede", ""))}</p></header>')
    (out / "cli" / "index.html").write_text(page("cli/index.html", "CAPS CLI guide", hero + body, toc_from(toc), maths=True,
                                                 desc=meta.get("lede", "")))


def build_tutorials(out: Path) -> None:
    d = out / "tutorials"
    d.mkdir(parents=True, exist_ok=True)
    files = sorted((SRC / "tutorials").glob("*.md"))
    cards = []
    metas = []
    for f in files:
        body, toc, meta = markdown(f.read_text())
        name = f.stem.split("-", 1)[1] + ".html"
        metas.append((name, meta))
    for k, f in enumerate(files):
        body, toc, meta = markdown(f.read_text())
        name, _ = metas[k]
        prev = metas[k - 1] if k > 0 else None
        nxt = metas[k + 1] if k + 1 < len(metas) else None
        pager = '<nav class="pager">' + (f'<a href="{prev[0]}"><small>Previous</small>{html.escape(prev[1]["title"])}</a>' if prev else "<span></span>") \
            + (f'<a class="next" href="{nxt[0]}"><small>Next</small>{html.escape(nxt[1]["title"])}</a>' if nxt else "<span></span>") + "</nav>"
        facts = "".join(f"<span><i>{html.escape(a)}</i> {html.escape(meta[b])}</span>" for a, b in (("Time", "time"), ("Level", "level"), ("Uses", "uses")) if meta.get(b))
        hero = (f'<header class="hero"><p class="eyebrow">Tutorial {k + 1:02d}</p><h1>{html.escape(meta["title"])}</h1>'
                f'<p>{inline(meta.get("lede", ""))}</p><div class="strip">{facts}</div></header>')
        (d / name).write_text(page(f"tutorials/{name}", meta["title"] + " · CAPS tutorial", hero + body + pager, toc_from(toc), maths=True,
                                   desc=meta.get("lede", "")))
        cards.append(f'<a class="card" href="{name}"><span class="num">{k + 1:02d}</span><b>{html.escape(meta["title"])}</b>'
                     f'<span>{inline(meta.get("lede", ""))}</span><span class="meta">{html.escape(meta.get("time", ""))} · {html.escape(meta.get("uses", ""))}</span></a>')
    hero = ('<header class="hero"><p class="eyebrow">Learn by doing</p><h1>Tutorials</h1><p>Each tutorial is a complete route '
            'through one kind of project, with commands that were run as written. Runs are kept short so they finish on a '
            'laptop; the text says where a production study needs longer ones.</p></header>')
    (d / "index.html").write_text(page("tutorials/index.html", "CAPS tutorials", hero + '<div class="cards">' + "".join(cards) + "</div>",
                                       desc="Step-by-step CAPS tutorials: cells, force fields, crosslinking, coarse-graining, DFT surfaces, Python."))


class Refs:
    def __init__(self, lib: str):
        self.lib = None
        if lib and Path(lib).exists():
            self.lib = ctypes.CDLL(lib)
            self.lib.caps_citation_text.argtypes = [ctypes.c_char_p, ctypes.c_char_p, ctypes.c_int32]
            self.lib.caps_citation_text.restype = ctypes.c_int32

    def text(self, key: str) -> str:
        if not self.lib:
            return ""
        buf = ctypes.create_string_buffer(4096)
        n = self.lib.caps_citation_text(key.encode(), buf, 4096)
        return buf.value.decode() if n >= 0 else ""


def ref_html(t: str) -> str:
    e = html.escape(t)
    return re.sub(r"doi:(\S+?)(\.?)$", lambda m: f'<a href="https://doi.org/{m.group(1)}">doi:{m.group(1)}</a>{m.group(2)}', e)


def build_theory(out: Path, refs: Refs) -> None:
    pages = json.loads((ROOT / "data/manual/manual.json").read_text())["pages"]
    d = out / "theory"
    d.mkdir(parents=True, exist_ok=True)
    groups: dict[str, list] = {}
    for p in pages:
        groups.setdefault(p["group"], []).append(p)
    order = ["Integrators", "Thermostats", "Barostats", "Minimisers", "Electrostatics", "Force fields", "Charges", "Growth",
             "Packing", "Protocols", "Analysis", "Coarse-grained", "DFT surfaces"]
    order += [g for g in groups if g not in order]
    cards = []
    for g in order:
        if g not in groups:
            continue
        gid = slug(g)
        body = [f'<header class="hero"><p class="eyebrow">Theory · {len(groups[g])} methods</p><h1>{html.escape(g)}</h1>'
                '<p>What CAPS computes, as it computes it: the equation, its symbols with CAPS\'s defaults, when to use the '
                'method, the source file, the tests that check it, and any departure from the cited method.</p></header>']
        toc = [f'<h2>{html.escape(g)}</h2>']
        for p in groups[g]:
            toc.append(f'<a href="#{p["id"]}">{html.escape(p["title"])}</a>')
            tex = p["tex"]
            if "\\\\" in tex and "\\begin" not in tex:
                tex = "\\begin{gathered}" + tex + "\\end{gathered}"
            sym = "".join(f"<tr><td>{html.escape(a)}</td><td>{html.escape(b)}</td><td>{html.escape(c)}</td></tr>" for a, b, c in p["symbols"])
            cites = "".join(f'<li id="{p["id"]}-{k}">{ref_html(refs.text(k)) or html.escape(k)}</li>' for k in p["cites"])
            body.append(f'<section class="method" id="{p["id"]}"><h2>{html.escape(p["title"])}</h2>'
                        f'<p>{html.escape(p["summary"])}</p>'
                        f'<div class="eq">$$' + html.escape(tex, quote=False) + '$$</div>'
                        + (f'<p class="eqnote">{html.escape(p["equation_note"])}</p>' if p.get("equation_note") else "")
                        + (f'<div class="table"><table class="sym"><thead><tr><th>Symbol</th><th>Meaning</th><th>In CAPS</th></tr></thead><tbody>{sym}</tbody></table></div>' if sym else "")
                        + f'<h3>When to use it</h3><p>{html.escape(p["when"])}</p>'
                        + '<dl class="facts">'
                        + f'<dt>Source</dt><dd><a href="{REPO}/blob/main/{html.escape(p["source"].split(" ")[0])}"><code>{html.escape(p["source"])}</code></a></dd>'
                        + f'<dt>Tested by</dt><dd>{html.escape(p["tested"])}</dd>'
                        + f'<dt>Departure from the reference</dt><dd>{html.escape(p["deviation"])}</dd>'
                        + "</dl>"
                        + (f'<h3>References</h3><ol class="refs">{cites}</ol>' if cites else "")
                        + "</section>")
        (d / f"{gid}.html").write_text(page(f"theory/{gid}.html", f"{g} · CAPS theory", "\n".join(body), "".join(toc), maths=True,
                                            desc=f"The theory behind CAPS's {g.lower()} methods."))
        cards.append(f'<a class="card" href="{gid}.html"><b>{html.escape(g)}</b><span>'
                     + html.escape(" · ".join(p["title"] for p in groups[g][:6]) + (" …" if len(groups[g]) > 6 else ""))
                     + f'</span><span class="meta">{len(groups[g])} methods</span></a>')
    hero = ('<header class="hero"><p class="eyebrow">Theory</p><h1>The methods behind CAPS</h1><p>One page per family of methods: '
            f'{len(pages)} methods, each with its equation, the symbols and their defaults in CAPS, when to use it, where it '
            'lives in the source, which tests check it, and the papers it follows. The same pages open from the Studio\'s '
            'Manual.</p></header>')
    (d / "index.html").write_text(page("theory/index.html", "CAPS theory", hero + '<div class="cards">' + "".join(cards) + "</div>",
                                       desc="The equations and references behind every CAPS method."))


def build_python(out: Path, lib: str) -> None:
    os.environ["CAPS_LIB"] = lib
    sys.path.insert(0, str(ROOT / "data/python"))
    import caps  # noqa: E402

    def sig(f, name):
        try:
            s = str(inspect.signature(f))
        except (TypeError, ValueError):
            s = "(…)"
        return name + s

    def entry(name, obj, kind):
        doc = inspect.getdoc(obj) or ""
        return (f'<article class="api" id="{html.escape(name)}"><h3><span class="kind">{kind}</span><code>{html.escape(sig(obj, name))}</code></h3>'
                + (f'<div class="doc">{html.escape(doc)}</div>' if doc else "") + "</article>")

    intro, toc_h, _ = markdown((SRC / "python.md").read_text())
    body = [intro]
    toc = [toc_from(toc_h, "Guide")]
    funcs = [(n, getattr(caps, n)) for n in caps.__all__ if inspect.isfunction(getattr(caps, n, None))]
    body.append('<h2 id="functions">Functions</h2>')
    toc.append('<h2>Functions</h2>' + "".join(f'<a href="#{n}"><code>{n}</code></a>' for n, _ in funcs))
    body += [entry(f"caps.{n}", f, "function") for n, f in funcs]
    doc_methods = [(n, m) for n, m in inspect.getmembers(caps.Document) if not n.startswith("_") and (inspect.isfunction(m) or isinstance(m, property))]
    body.append('<h2 id="document">caps.Document</h2><p>' + html.escape(inspect.getdoc(caps.Document) or "") + "</p>")
    toc.append('<h2>Document</h2>' + "".join(f'<a href="#Document.{n}"><code>.{n}</code></a>' for n, _ in doc_methods))
    for n, m in doc_methods:
        if isinstance(m, property):
            d = inspect.getdoc(m) or ""
            body.append(f'<article class="api" id="Document.{n}"><h3><span class="kind">property</span><code>Document.{n}</code></h3>'
                        + (f'<div class="doc">{html.escape(d)}</div>' if d else "") + "</article>")
        else:
            body.append(entry(f"Document.{n}", m, "method").replace(f'id="Document.{n}"', f'id="Document.{n}"'))
    builders = [(n, m) for n, m in inspect.getmembers(caps.build) if not n.startswith("_") and callable(m)]
    body.append('<h2 id="build">caps.build</h2><p>' + html.escape(inspect.getdoc(caps.build) or "") + "</p>")
    toc.append('<h2>build</h2>' + "".join(f'<a href="#build.{n}"><code>.{n}</code></a>' for n, _ in builders))
    body += [entry(f"build.{n}", m, "builder") for n, m in builders]
    for modname in ("pipeline", "sweep", "geometry"):
        mod = getattr(caps, modname)
        members = [(n, m) for n, m in inspect.getmembers(mod) if not n.startswith("_") and (inspect.isfunction(m) or inspect.isclass(m)) and getattr(m, "__module__", "") == mod.__name__]
        body.append(f'<h2 id="{modname}">caps.{modname}</h2><div class="doc">' + html.escape(inspect.getdoc(mod) or "") + "</div>")
        toc.append(f'<h2>{modname}</h2>' + "".join(f'<a href="#{modname}.{n}"><code>.{n}</code></a>' for n, _ in members))
        body += [entry(f"{modname}.{n}", m, "class" if inspect.isclass(m) else "function") for n, m in members]
    hero = ('<header class="hero"><p class="eyebrow">Python · <code>import caps</code></p><h1>Python API</h1><p>The same core as the '
            'command line and the Studio, through its C library: open and build structures, assign force fields, relax, run, '
            'crosslink, analyse, export, and write pipeline steps the Studio runs. Every entry below is the docstring of the '
            f'installed package (ABI {caps.abi_version()}).</p></header>')
    (out / "python.html").write_text(page("python.html", "CAPS Python API", hero + "\n".join(body), "".join(toc),
                                          desc="The caps Python package: functions, Document methods, builders, pipeline steps."))


def build_home(out: Path, n_methods: int, n_cmds: int) -> None:
    md = (SRC / "home.html").read_text()
    n_ff = len([f for f in (ROOT / "data/forcefields").glob("*.json") if f.name != "catalogue.json"])
    n_tut = len(list((SRC / "tutorials").glob("*.md")))
    md = md.replace("{N_METHODS}", str(n_methods)).replace("{N_CMDS}", str(n_cmds)).replace("{N_FF}", str(n_ff)).replace("{N_TUT}", str(n_tut))
    (out / "index.html").write_text(page("index.html", "CAPS documentation", md, wide=True,
                                         desc="CAPS: build, type, run and analyse polymers and materials from the Studio, the command line or Python."))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--caps", default=str(ROOT / "build/cli/caps"))
    lib_default = next((str(p) for p in (ROOT / "build/capi").glob("libcaps.*") if p.suffix in (".dylib", ".so", ".dll")), "")
    ap.add_argument("--lib", default=lib_default)
    ap.add_argument("--out", default=str(ROOT / "docs"))
    a = ap.parse_args()
    out = Path(a.out)
    (out / "assets").mkdir(parents=True, exist_ok=True)
    for f in ("site.css", "site.js"):
        (out / "assets" / f).write_text((SRC / f).read_text())
    (out / ".nojekyll").write_text("")
    summary = build_reference(a.caps, out)
    build_guide(out, summary)
    build_tutorials(out)
    build_theory(out, Refs(a.lib))
    build_python(out, a.lib)
    n_methods = len(json.loads((ROOT / "data/manual/manual.json").read_text())["pages"])
    build_home(out, n_methods, len(summary))
    print(f"wrote {out}: {len(summary)} commands, {n_methods} methods")


if __name__ == "__main__":
    main()
