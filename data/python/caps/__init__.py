"""CAPS Python steps: a pipeline step written in Python.

CAPS runs the script in a separate Python process with the frame's particles, bonds, cell and attributes. The script
defines one function decorated with @step; it reads data.particles and writes data.attributes, data.tables and new
per-particle properties, which CAPS shows in the data inspector and passes to the next steps.

    import numpy as np
    import caps
    from caps.pipeline import step

    @step(name="Backbone conformation")
    def modify(frame, data):
        mol = data.particles["Molecule Identifier"]
        pos = data.particles.positions_unwrapped
        data.attributes["MeanZ"] = float(np.mean(pos[:, 2]))
"""
from . import geometry, pipeline  # noqa: F401

__all__ = ["geometry", "pipeline"]
