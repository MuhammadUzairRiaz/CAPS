"""caps.sweep: the Studio's Sweep from Python — one polymer over tacticities, chain lengths and seeds, each cell grown,
relaxed (optionally run NPT) and analysed, saved with its provenance; and the results of a finished sweep.

    runs = caps.sweep.run("*CC(*)c1ccccc1", tacticities=("isotactic", "syndiotactic", "atactic"), dps=(20,), seeds=(1, 2, 3))
    cells = [caps.sweep.result(t) for t in ("isotactic", "syndiotactic", "atactic")]
    caps.table(cells, ["density", "rg", "c_inf"])

A sweep folder holds TACTICITY_dpN_seedS.data (+ .provenance.json) for every run and results.json with the values;
the Studio writes the same layout under ~/CAPS/sweeps.
"""
from __future__ import annotations

import glob
import json
import math
import os
import re
import time
from pathlib import Path
from typing import Iterable, Optional

from .core import CapsError, _PROPERTY_ALIASES, open as _open, run as _run

SWEEPS = Path.home() / "CAPS" / "sweeps"
_DEFAULT = ("density", "rg", "cn")


def _pid(p: str) -> str:
    return _PROPERTY_ALIASES.get(p, p)


def run(polymer, tacticities: Iterable[str] = ("isotactic", "syndiotactic", "atactic"), dps: Iterable[int] = (20,),
        seeds: Iterable[int] = (1, 2, 3), chains: int = 10, density: float = 0.5, forcefield: str = "default", relax: bool = True,
        npt_ps: float = 0.0, temperature: float = 300.0, properties: Iterable[str] = _DEFAULT, tg: Optional[dict] = None,
        folder: Optional[str] = None, progress=None) -> "Sweep":
    """Runs every (tacticity, DP, seed); failures are recorded, not raised. tg={"t_start": 500, "t_end": 200, "t_step": 25,
    "ps_per_step": 20} adds a stepwise-cooling Tg to each run (slow). progress(run, status) is told as runs finish."""
    units = [polymer] if isinstance(polymer, str) else list(polymer)
    folder = Path(folder) if folder else SWEEPS / ("sweep-" + time.strftime("%Y%m%d-%H%M"))
    folder.mkdir(parents=True, exist_ok=True)
    props = [_pid(p) for p in properties]
    runs = []
    for tac in tacticities:
        for dp in dps:
            for seed in seeds:
                name = f"{tac}_dp{dp}_seed{seed}"
                rec = {"recipe": 1, "name": name,
                       "build": {"polymer": {"units": units, "dp": dp, "chains": chains, "tacticity": tac}},
                       "type": {"forcefield": forcefield}, "grow": {"density": density, "seed": seed},
                       "analyze": {"properties": props}, "export": ["lammps"]}
                if relax:
                    rec["relax"] = {"method": "lbfgs", "fmax": 1.0}
                if npt_ps > 0:
                    rec["md"] = {"ps": npt_ps, "temperature": temperature, "ensemble": "npt", "seed": seed}
                entry = {"tacticity": tac, "dp": dp, "seed": seed, "file": name + ".data", "status": "running", "properties": {}}
                try:
                    doc = _run(rec, out_dir=str(folder))
                    for p in doc.properties:
                        if p.get("value") is not None and math.isfinite(p["value"]):
                            entry["properties"][p["id"]] = {"value": p["value"], "unit": p.get("unit", ""), "name": p.get("name", p["id"])}
                    if tg:
                        for p in doc.analyze("tg", **tg):
                            if p.get("value") is not None and math.isfinite(p["value"]):
                                entry["properties"]["tg"] = {"value": p["value"], "unit": p.get("unit", "K"), "name": "Tg"}
                    entry["status"] = "done"
                except CapsError as e:
                    entry["status"], entry["error"] = "failed", str(e)
                runs.append(entry)
                _write(folder, units, runs)
                if progress:
                    progress(entry, entry["status"])
    return Sweep(folder)


def _write(folder: Path, units, runs) -> None:
    (folder / "results.json").write_text(json.dumps({"format": "caps-sweep/1", "units": units, "runs": runs}, indent=1))


def latest(root: Optional[str] = None) -> Path:
    """The newest sweep folder (under ~/CAPS/sweeps)."""
    base = Path(root) if root else SWEEPS
    dirs = [d for d in base.glob("*") if d.is_dir() and ((d / "results.json").exists() or any(d.glob("*_dp*_seed*.data")))]
    if not dirs:
        raise CapsError(f"no sweep in {base}")
    return max(dirs, key=lambda d: d.stat().st_mtime)


class Sweep:
    """A sweep folder: its runs and the result of each condition."""

    def __init__(self, folder):
        self.folder = Path(folder)

    @property
    def runs(self) -> list:
        f = self.folder / "results.json"
        if f.exists():
            return json.loads(f.read_text())["runs"]
        runs = []   # a folder without results.json (an older Studio sweep): the files' names
        for p in sorted(self.folder.glob("*_dp*_seed*.data")):
            m = re.match(r"(.+)_dp(\d+)_seed(\d+)\.data$", p.name)
            if m:
                runs.append({"tacticity": m.group(1), "dp": int(m.group(2)), "seed": int(m.group(3)), "file": p.name, "status": "done", "properties": {}})
        return runs

    @property
    def conditions(self) -> list:
        seen = []
        for r in self.runs:
            if r["tacticity"] not in seen:
                seen.append(r["tacticity"])
        return seen

    def result(self, condition: str, dp: Optional[int] = None) -> "Result":
        runs = [r for r in self.runs if r["tacticity"] == condition and (dp is None or r["dp"] == dp)]
        if not runs:
            raise CapsError(f"no runs for {condition!r} in {self.folder} (conditions: {', '.join(self.conditions)})")
        return Result(condition if dp is None else f"{condition} · DP {dp}", runs, self.folder)

    def __repr__(self) -> str:
        runs = self.runs
        done = sum(r["status"] == "done" for r in runs)
        return f"<caps.Sweep {self.folder} · {done}/{len(runs)} done>"


def result(condition: str, folder: Optional[str] = None, dp: Optional[int] = None) -> "Result":
    """One condition (a tacticity) of a sweep, its seeds pooled: the newest sweep unless folder is given."""
    return Sweep(folder or latest()).result(condition, dp)


class Result:
    """The runs of one condition: r["density"] → (mean, standard deviation, n) over the runs that have it."""

    def __init__(self, label: str, runs: list, folder: Path):
        self.label, self.runs, self.folder = label, runs, Path(folder)

    def values(self, prop: str) -> list:
        pid = _pid(prop)
        out = []
        for r in self.runs:
            if r.get("status") != "done":
                continue
            p = r.get("properties", {}).get(pid)
            if p is None and pid != "tg" and not r.get("properties"):   # a Studio folder: analyse the saved cell
                try:
                    doc = _open(str(self.folder / r["file"]))
                    for q in doc.analyze([x for x in _DEFAULT]):
                        if q.get("value") is not None and math.isfinite(q["value"]):
                            r.setdefault("properties", {})[q["id"]] = {"value": q["value"], "unit": q.get("unit", ""), "name": q.get("name", q["id"])}
                    p = r["properties"].get(pid)
                except CapsError:
                    p = None
            if p is not None:
                out.append(p["value"])
        return out

    def unit(self, prop: str) -> str:
        pid = _pid(prop)
        for r in self.runs:
            p = r.get("properties", {}).get(pid)
            if p:
                return p.get("unit", "")
        return {"density": "g/cm³", "rg": "Å", "tg": "K", "cn": ""}.get(pid, "")

    def __getitem__(self, prop: str):
        v = self.values(prop)
        if not v:
            return (math.nan, math.nan, 0)
        m = sum(v) / len(v)
        sd = math.sqrt(sum((x - m) ** 2 for x in v) / (len(v) - 1)) if len(v) > 1 else math.nan
        return (m, sd, len(v))

    def __repr__(self) -> str:
        done = sum(r.get("status") == "done" for r in self.runs)
        return f"<caps.sweep.Result {self.label} · {done}/{len(self.runs)} runs>"
