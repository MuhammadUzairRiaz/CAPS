"""python -m caps.runner SCRIPT IN.json OUT.json — what CAPS runs for a Python step."""
import json
import runpy
import sys
import traceback


def _plain(v):
    if hasattr(v, "tolist"):
        return v.tolist()
    if isinstance(v, dict):
        return {str(k): _plain(x) for k, x in v.items()}
    if isinstance(v, (list, tuple)):
        return [_plain(x) for x in v]
    return v


def _table(name, t):
    """A table as columns and rows: a dict of key → value (or → list), or {"columns": [...], "rows": [[...]]}."""
    t = _plain(t)
    if isinstance(t, dict) and "columns" in t and "rows" in t:
        return {"name": name, "columns": [str(c) for c in t["columns"]], "rows": [[float(x) for x in r] for r in t["rows"]]}
    if isinstance(t, dict):
        rows = []
        for k, v in t.items():
            try:
                key = float(k)
            except ValueError:
                key = float(len(rows) + 1)
            vals = v if isinstance(v, list) else [v]
            rows.append([key] + [float(x) for x in vals])
        width = max((len(r) for r in rows), default=2)
        cols = ["key"] + [name if width == 2 else f"{name}.{i}" for i in range(1, width)]
        return {"name": name, "columns": cols, "rows": rows}
    raise TypeError(f"table {name}: a dict or {{'columns', 'rows'}}")


def main():
    script, src, dst = sys.argv[1:4]
    from caps import pipeline
    from caps.pipeline import Data
    with open(src) as f:
        raw = json.load(f)
    out = {"ok": False, "log": []}
    try:
        runpy.run_path(script, run_name="caps_step")
        if not pipeline._registered:
            raise RuntimeError("the script defines no @step function")
        name, fn = pipeline._registered[-1]
        data = Data(raw)
        fn(raw.get("frame", 0), data)
        out.update(ok=True, name=name,
                   attributes={k: float(v) for k, v in _plain(data.attributes).items() if isinstance(v, (int, float)) or hasattr(v, "__float__")},
                   properties={k: [float(x) for x in _plain(v)] for k, v in data.particles.new.items()},
                   tables=[_table(k, v) for k, v in data.tables.items()])
        if data.selection is not None:
            out["selection"] = [1 if x else 0 for x in _plain(data.selection)]
    except Exception as e:  # reported on the step, with the line
        out["error"] = f"{type(e).__name__}: {e}"
        out["trace"] = traceback.format_exc()
    with open(dst, "w") as f:
        json.dump(out, f)


if __name__ == "__main__":
    main()
