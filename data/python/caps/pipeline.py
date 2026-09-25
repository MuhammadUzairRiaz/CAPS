"""The @step decorator and the data a step sees (particles, bonds, cell, attributes, tables)."""

_registered = []


def step(name=None):
    """Marks the function CAPS calls: modify(frame, data)."""
    def wrap(fn):
        _registered.append((name or fn.__name__, fn))
        return fn
    return wrap


def _array(values):
    try:
        import numpy as np
        return np.asarray(values)
    except ImportError:  # numpy is optional: plain lists then
        return values


class Particles:
    """Per-particle columns by name: "Particle Identifier", "Molecule Identifier", "Particle Type", "Element", "Charge",
    "Selection", "Backbone" (1 on backbone atoms), "Position" (N × 3) and any property earlier steps made."""

    def __init__(self, columns, count):
        self._columns = {k: _array(v) for k, v in columns.items()}
        self.count = count
        self.new = {}

    def __len__(self):
        return self.count

    def __contains__(self, key):
        return key in self._columns or key in self.new

    def __getitem__(self, key):
        if key in self.new:
            return self.new[key]
        if key not in self._columns:
            raise KeyError(f"no particle property {key!r}; have {sorted(self._columns)}")
        return self._columns[key]

    def __setitem__(self, key, values):
        values = list(values.tolist() if hasattr(values, "tolist") else values)
        if len(values) != self.count:
            raise ValueError(f"{key}: {len(values)} values for {self.count} particles")
        self.new[key] = _array(values)

    def keys(self):
        return list(self._columns) + [k for k in self.new if k not in self._columns]

    @property
    def positions(self):
        return self._columns["Position"]

    positions_unwrapped = positions   # CAPS hands over whole molecules (unwrapped)


class Data:
    def __init__(self, raw):
        self.particles = Particles(raw["particles"], raw["count"])
        self.bonds = _array(raw.get("bonds", []))
        self.cell = raw.get("cell")
        self.attributes = dict(raw.get("attributes", {}))
        self.tables = {}
        self.selection = None   # set to a list of 0/1 to select particles
