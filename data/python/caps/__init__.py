"""CAPS from Python: the core through its C library (caps.open, caps.build.*, Document.relax / md / edit / save …),
and pipeline steps written in Python.

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
from .core import (CapsError, Document, Provenance, abi_version, bibtex, build, compare_provenance, current, hand_back, import_file, import_preview, library, methods,  # noqa: F401
                   open, pack, polymer, potentials, provenance_file, run, space_groups, chain_lengths, bead_templates, copolymer_model, stereo, blend_phase, solvent_chi,
                   ewald_params, chi_by_md, chi_by_contacts, reaction_templates, reaction_template)
from . import sweep  # noqa: F401,E402
from .table import Table, table  # noqa: F401,E402
from .view import View  # noqa: F401,E402

__version__ = "0.1.0"
__all__ = ["geometry", "pipeline", "sweep", "Document", "Provenance", "View", "Table", "open", "pack", "potentials", "current", "hand_back", "import_file", "import_preview", "build", "run", "polymer",
           "table", "library", "abi_version", "space_groups", "provenance_file", "compare_provenance", "bibtex", "methods", "CapsError",
           "chain_lengths", "bead_templates", "copolymer_model", "stereo", "blend_phase", "solvent_chi", "ewald_params", "chi_by_md", "chi_by_contacts", "reaction_templates", "reaction_template"]
