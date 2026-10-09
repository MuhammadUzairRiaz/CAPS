## Setup

The package is plain Python 3 (standard library only; NumPy for pipeline steps) and ships with the program. Put its folder on `PYTHONPATH`; it finds the CAPS library beside it.

| System | `PYTHONPATH` |
|---|---|
| macOS | `/Applications/CAPS Studio.app/Contents/Resources/data/python` |
| Windows | `<install folder>\data\python` |
| Linux (.deb) | `/opt/caps/data/python` |
| Linux (tarball) | `<unpacked folder>/data/python` |

```bash
export PYTHONPATH="/Applications/CAPS Studio.app/Contents/Resources/data/python"
python3 -c "import caps; print(caps.abi_version())"
```

`CAPS_LIB` overrides the library's path (a development build, for instance).

## A first script

Grow a natural-rubber cell, type it with OPLS-AA, cure it with a peroxide template, and write LAMMPS and GROMACS inputs. Every call is a method of a `Document` (one structure with its frames, force field and provenance).

```python
import caps

d = caps.polymer("*CC(C)=CC*", dp=10, chains=6, density=0.6, seed=1)   # cis-1,4-polyisoprene
d.relax(max_iterations=300)
d.field.assign("opls2005")                     # types and charges from the library
print(d.react("peroxide_allylic", cycles=5, per_cycle=4, target=0.3, seed=1))
d.save("nr_xl.data")
d.export_engines("engines", stem="nr", run="npt", temperature=300, steps=5000)
```

The [crosslinking tutorial](tutorials/crosslinking.html) goes through each step, and the LAMMPS fix bond/react export.

## Pipeline steps in Python

The Studio's analysis pipeline (and `caps pipeline`) runs Python steps in a separate process. A step is one function decorated with `@step`; it reads the frame's particles and writes attributes, tables or new per-particle properties.

```python
import numpy as np
from caps.pipeline import step

@step(name="Backbone height")
def modify(frame, data):
    pos = data.particles.positions_unwrapped
    data.attributes["MeanZ"] = float(np.mean(pos[:, 2]))
```

The reference below is generated from the installed package's docstrings.
