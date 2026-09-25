"""CAPS from Python: the C library (libcaps) through ctypes, for scripts, macros and notebooks.

    import caps
    doc = caps.open("PS_melt.data")
    print(doc.summary())
    pep = caps.build.peptide("AEAAAKEAAAKA", structure="HHHHHHHHHHHH")
    pep.field.assign("uff")            # or a force field id from the library, e.g. "gaff2"
    pep.relax(ftol=0.5)
    pep.save("peptide.data")

The library is found through CAPS_LIB (a path to libcaps), next to the Studio, or in the build tree.
"""
from __future__ import annotations

import ctypes as C
import json
import os
import sys
from pathlib import Path
from typing import Any, Iterable, Optional

__all__ = ["Document", "open", "build", "library", "abi_version", "CapsError"]


class CapsError(RuntimeError):
    """An error from the CAPS core (the message is caps_last_error)."""


# ---------------------------------------------------------------------------------------------------------------- library

def _candidates() -> Iterable[Path]:
    env = os.environ.get("CAPS_LIB")
    if env:
        yield Path(env)
    names = {"darwin": ["libcaps.dylib"], "win32": ["caps.dll", "libcaps.dll"]}.get(sys.platform, ["libcaps.so"])
    here = Path(__file__).resolve().parent            # …/data/python/caps
    roots = [here.parents[2] / "build" / "capi",      # the source tree: build/capi
             here.parents[2],                          # an installed app: data/ sits beside the library
             here.parents[3] / "MacOS",                # the macOS bundle: Contents/Resources/data → Contents/MacOS
             here.parents[2] / "studio" / "CapsStudio" / "bin" / "Release" / "net10.0"]
    for r in roots:
        for n in names:
            yield r / n


_lib: Optional[C.CDLL] = None


def library() -> C.CDLL:
    """The loaded libcaps (loaded on first use)."""
    global _lib
    if _lib is not None:
        return _lib
    tried = []
    for p in _candidates():
        tried.append(str(p))
        if p.exists():
            _lib = C.CDLL(str(p))
            _declare(_lib)
            return _lib
    raise CapsError("libcaps not found; set CAPS_LIB to its path (tried: " + ", ".join(tried) + ")")


class _Summary(C.Structure):
    _fields_ = [("atoms", C.c_int64), ("bonds", C.c_int64), ("molecules", C.c_int64), ("frames", C.c_int64),
                ("bonds_from_file", C.c_int32), ("has_charges", C.c_int32), ("cell_valid", C.c_int32), ("unwrapped", C.c_int32),
                ("cell_a", C.c_double), ("cell_b", C.c_double), ("cell_c", C.c_double), ("volume", C.c_double), ("density", C.c_double),
                ("total_mass", C.c_double), ("total_charge", C.c_double), ("format", C.c_char * 24)]


class _Atom(C.Structure):
    _fields_ = [("id", C.c_int64), ("mol", C.c_int64), ("type", C.c_int32), ("element", C.c_int32), ("index", C.c_int32),
                ("charge", C.c_double), ("x", C.c_double), ("y", C.c_double), ("z", C.c_double),
                ("element_symbol", C.c_char * 4), ("name", C.c_char * 16)]


class _Camera(C.Structure):
    _fields_ = [("yaw", C.c_double), ("pitch", C.c_double), ("zoom", C.c_double), ("pan_x", C.c_double), ("pan_y", C.c_double),
                ("perspective", C.c_int32)]


class _RenderOpts(C.Structure):
    _fields_ = [("width", C.c_int32), ("height", C.c_int32), ("supersample", C.c_int32), ("background", C.c_int32),
                ("custom_rgb", C.c_uint32), ("colour_by", C.c_int32), ("style", C.c_int32), ("outlines", C.c_int32),
                ("depth_cue", C.c_int32), ("show_cell", C.c_int32), ("highlight", C.c_int32 * 4), ("focus", C.c_int32),
                ("ambient_occlusion", C.c_int32)]


class _RelaxOpts(C.Structure):
    _fields_ = [("method", C.c_int32), ("ftol", C.c_double), ("max_iterations", C.c_int32), ("target_density", C.c_double),
                ("compress_step", C.c_double), ("pushoff", C.c_int32), ("relax_box", C.c_int32), ("pressure", C.c_double),
                ("cutoff", C.c_double), ("coulomb", C.c_int32), ("threads", C.c_int32)]


class _MdOpts(C.Structure):
    _fields_ = [("dt", C.c_double), ("steps", C.c_int64), ("temperature", C.c_double), ("thermostat", C.c_int32), ("tau_t", C.c_double),
                ("barostat", C.c_int32), ("pressure", C.c_double), ("tau_p", C.c_double), ("new_velocities", C.c_int32),
                ("seed", C.c_uint64), ("thermo_every", C.c_int32), ("frame_every", C.c_int32), ("cutoff", C.c_double),
                ("coulomb", C.c_int32), ("tail", C.c_int32), ("threads", C.c_int32)]


class _BuildOpts(C.Structure):
    _fields_ = [("conformers", C.c_int32), ("seed", C.c_uint64)]


def _declare(L: C.CDLL) -> None:
    P, S, I, D, B = C.c_void_p, C.c_char_p, C.c_int32, C.c_double, C.c_char_p
    sig = {
        "caps_abi_version": ([], I), "caps_last_error": ([], S),
        "caps_open": ([S, S], P), "caps_close": ([P], None), "caps_import": ([S, S, S], P), "caps_import_preview": ([S, S, B, I], I),
        "caps_summary_get": ([P, C.POINTER(_Summary)], I), "caps_set_frame": ([P, C.c_int64], I),
        "caps_atom": ([P, I, C.POINTER(_Atom)], I), "caps_save": ([P, S], I), "caps_save_trajectory": ([P, S], I),
        "caps_export_png": ([P, C.POINTER(_Camera), C.POINTER(_RenderOpts), S], I),
        "caps_relax": ([P, C.POINTER(_RelaxOpts), P, P, B, I], I), "caps_md": ([P, C.POINTER(_MdOpts), P, P, B, I], I),
        "caps_field_assign": ([P, S, S, I], I), "caps_field_report": ([P, B, I], I),
        "caps_build_smiles": ([S, S, C.POINTER(_BuildOpts), B, I], P),
        "caps_peptide_build": ([S, B, I], P), "caps_crystal_build": ([S, B, I], P), "caps_nano_build": ([S, B, I], P),
        "caps_solvate": ([P, S, P, P, B, I], P),
        "caps_edit": ([P, S, B, I], I), "caps_undo": ([P, I], I), "caps_select": ([P, S, B, I], I), "caps_selection": ([P, B, I], I),
        "caps_tacticity": ([P, B, I], I), "caps_interactions": ([P, S, B, I], I), "caps_torsion_scan": ([P, S, P, P, B, I], I),
        "caps_trajectory_series": ([P, S, P, P, B, I], I), "caps_file_checks": ([P, B, I], I),
        "caps_space_groups": ([B, I], I), "caps_crystal_find_symmetry": ([S, S, D, B, I], I),
    }
    for name, (args, res) in sig.items():
        f = getattr(L, name, None)
        if f is None:
            continue
        f.argtypes = args
        f.restype = res


def abi_version() -> int:
    return library().caps_abi_version()


def _error() -> CapsError:
    msg = library().caps_last_error()
    return CapsError(msg.decode() if msg else "unknown error")


def _enc(s: Optional[str]) -> Optional[bytes]:
    return None if s is None else str(s).encode()


def _json_call(fn, *args) -> Any:
    """Calls a caps_… that writes JSON into (buffer, capacity) and returns the length needed."""
    cap = 1 << 16
    for _ in range(3):
        buf = C.create_string_buffer(cap)
        n = fn(*args, buf, cap)
        if n < 0:
            raise _error()
        if n <= cap:
            text = buf.value.decode()
            return json.loads(text) if text else {}
        cap = n + 1
    raise CapsError("reply too large")


def _report() -> C.Array:
    return C.create_string_buffer(8192)


# ---------------------------------------------------------------------------------------------------------------- document

class _Field:
    def __init__(self, doc: "Document"):
        self._doc = doc

    def assign(self, forcefield: str = "uff", charges: str = "forcefield") -> dict:
        """Types every atom and looks up every parameter: a force-field id from the library (gaff2, opls2005 …), a path
        to a caps-forcefield JSON, or "uff". charges: forcefield | gasteiger | keep | qeq."""
        path = _forcefield_path(forcefield)
        code = {"forcefield": 0, "gasteiger": 1, "keep": 2, "qeq": 3}[charges]
        rc = library().caps_field_assign(self._doc._h, _enc(path), None, code)
        if rc < 0:
            raise _error()
        rep = _json_call(library().caps_field_report, self._doc._h)
        rep["complete"] = rc == 0
        return rep


def _forcefield_path(ff: str) -> str:
    if ff.lower() == "uff" or os.path.exists(ff):
        return ff if ff.lower() != "uff" else "uff"
    lib = Path(__file__).resolve().parents[1].parent / "forcefields"
    aliases = {"gaff2": "gaff-amber25-dlfield", "gaff": "gaff-amber16-dlfield", "opls": "opls2005-dlfield", "pcff": "pcff-dlfield",
               "compass": "compass-dlfield", "cvff": "cvff-dlfield", "amber": "amber-dlfield", "charmm": "charmm-dlfield"}
    name = aliases.get(ff.lower(), ff)
    p = lib / (name + ".json")
    if p.exists():
        return str(p)
    raise CapsError(f"unknown force field {ff!r} (a library id such as gaff2, a caps-forcefield JSON path, or uff)")


class Document:
    """A structure or trajectory: frames, a current frame, and everything CAPS does to it."""

    def __init__(self, handle: int, label: str = ""):
        if not handle:
            raise _error()
        self._h = C.c_void_p(handle)
        self.label = label
        self.field = _Field(self)
        self.report = ""

    # life
    def close(self) -> None:
        if self._h:
            library().caps_close(self._h)
            self._h = C.c_void_p()

    def __del__(self):
        try:
            self.close()
        except Exception:
            pass

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    def __repr__(self) -> str:
        s = self.summary()
        return f"<caps.Document {self.label or ''} · {s['atoms']} atoms · {s['frames']} frame(s)>"

    # contents
    def summary(self) -> dict:
        s = _Summary()
        if library().caps_summary_get(self._h, C.byref(s)) != 0:
            raise _error()
        return {k: (getattr(s, k).decode() if k == "format" else getattr(s, k)) for k, _ in _Summary._fields_}

    @property
    def atoms(self) -> int:
        return self.summary()["atoms"]

    @property
    def frames(self) -> int:
        return self.summary()["frames"]

    def set_frame(self, k: int) -> None:
        if library().caps_set_frame(self._h, k) != 0:
            raise _error()

    def atom(self, i: int) -> dict:
        a = _Atom()
        if library().caps_atom(self._h, i, C.byref(a)) != 0:
            raise _error()
        return {"id": a.id, "molecule": a.mol, "type": a.type, "element": a.element_symbol.decode(), "name": a.name.decode(),
                "charge": a.charge, "position": (a.x, a.y, a.z)}

    def positions(self) -> list:
        return [self.atom(i)["position"] for i in range(self.atoms)]

    # engines
    def relax(self, ftol: float = 0.5, method: str = "lbfgs", max_iterations: int = 5000, density: float = 0.0, pushoff: bool = True,
              box: bool = False, pressure: float = 1.0, cutoff: float = 10.0, coulomb: bool = True, threads: int = 0) -> int:
        """Minimises the current frame with the Field assignment (else the built-in GAFF): 0 converged, 1 not quite."""
        o = _RelaxOpts({"sd": 0, "cg": 1, "lbfgs": 2, "fire": 3}[method], ftol, max_iterations, density, 0.06, int(pushoff), int(box),
                       pressure, cutoff, int(coulomb), threads)
        rep = _report()
        rc = library().caps_relax(self._h, C.byref(o), None, None, rep, len(rep))
        if rc < 0:
            raise _error()
        self.report = rep.value.decode()
        return rc

    def md(self, steps: int = 10000, dt: float = 1.0, temperature: float = 300.0, thermostat: str = "bussi", barostat: str = "none",
           pressure: float = 1.0, seed: int = 1, frame_every: int = 1000, thermo_every: int = 100, cutoff: float = 10.0) -> str:
        """Molecular dynamics from the current frame; the frames recorded become the document's frames."""
        o = _MdOpts(dt, steps, temperature, {"none": 0, "bussi": 1, "langevin": 2}[thermostat], 100.0,
                    {"none": 0, "crescale": 1, "berendsen": 2}[barostat], pressure, 1000.0, 0, seed, thermo_every, frame_every,
                    cutoff, 1, 1, 0)
        rep = _report()
        if library().caps_md(self._h, C.byref(o), None, None, rep, len(rep)) < 0:
            raise _error()
        self.report = rep.value.decode()
        return self.report

    # editing, selection, analysis
    def edit(self, **op) -> dict:
        """One structure edit (the Studio's builder tools): edit(op="add_h"), edit(op="attach", target=5, smiles="*C(=O)O*") …"""
        r = _json_call(library().caps_edit, self._h, _enc(json.dumps(op)))
        if not r.get("ok"):
            raise CapsError(r.get("error", "edit failed"))
        return r

    def undo(self, redo: bool = False) -> None:
        if library().caps_undo(self._h, int(redo)) != 0:
            raise _error()

    def select(self, mode: str, pattern: str = "", op: str = "replace", **kw) -> int:
        r = _json_call(library().caps_select, self._h, _enc(json.dumps({"mode": mode, "pattern": pattern, "op": op, **kw})))
        if not r.get("ok"):
            raise CapsError(r.get("error", "selection failed"))
        return int(r["count"])

    def selection(self) -> list:
        return [int(i) for i in _json_call(library().caps_selection, self._h)["indices"]]

    def tacticity(self) -> dict:
        return _json_call(library().caps_tacticity, self._h)

    def interactions(self, **options) -> dict:
        return _json_call(library().caps_interactions, self._h, _enc(json.dumps(options)))

    def checks(self) -> list:
        return _json_call(library().caps_file_checks, self._h)

    def torsion_scan(self, atoms, step: float = 15, relax: bool = False, forcefield: str = "auto") -> dict:
        return _json_call(library().caps_torsion_scan, self._h, _enc(json.dumps({"atoms": list(atoms), "step": step, "relax": relax,
                                                                                  "forcefield": forcefield})), None, None)

    def series(self, molecule: int = 0, dt_fs: float = 1.0, log: str = "") -> dict:
        return _json_call(library().caps_trajectory_series, self._h, _enc(json.dumps({"molecule": molecule, "dt_fs": dt_fs, "log": log})),
                          None, None)

    # files
    def save(self, path: str) -> None:
        """Writes the current frame: .data (LAMMPS, with force-field sections when assigned), .pdb, .xyz, .mol2, .gro …"""
        if library().caps_save(self._h, _enc(str(path))) != 0:
            raise _error()

    def save_trajectory(self, path: str) -> None:
        if library().caps_save_trajectory(self._h, _enc(str(path))) != 0:
            raise _error()

    def render(self, path: str, width: int = 1280, height: int = 800, background: str = "white", style: str = "ball_and_stick",
               colour: str = "molecule", yaw: float = 0.55, pitch: float = 0.40, zoom: float = 1.0) -> None:
        """A PNG of the current frame."""
        cam = _Camera(yaw, pitch, zoom, 0, 0, 0)
        hl = (C.c_int32 * 4)(-1, -1, -1, -1)
        opt = _RenderOpts(width, height, 2, {"dark": 0, "white": 1, "transparent": 2}[background], 0,
                          {"element": 0, "molecule": 1, "type": 2}[colour],
                          {"ball_and_stick": 0, "space_filling": 1, "sticks": 2, "no_hydrogens": 3, "backbone": 4}[style], 1, 1, 1, hl, 0, 1)
        if library().caps_export_png(self._h, C.byref(cam), C.byref(opt), _enc(str(path))) != 0:
            raise _error()


def open(path: str, topology: Optional[str] = None) -> Document:  # noqa: A001 (the natural name)
    """Opens a structure or trajectory (LAMMPS data/dump, GROMACS .gro, PDB, XYZ, mol2, CIF …)."""
    return Document(library().caps_open(_enc(str(path)), _enc(topology)), str(path))


def import_file(path: str, bonds: str = "perceive", tolerance: float = 0.45, bond_orders: bool = True, split: bool = True,
                unwrap: bool = True, use_cell: bool = True, topology: Optional[str] = None) -> Document:
    """Opens a file that has no topology (XYZ, PDB, CIF …) with the Import choices: bonds "perceive" (covalent radii +
    tolerance Å), "file" or "none"; bond orders and aromaticity; molecules by connectivity; unwrapping; the cell."""
    opts = json.dumps({"bonds": bonds, "tolerance": tolerance, "bond_orders": bond_orders, "split": split, "unwrap": unwrap, "use_cell": use_cell})
    return Document(library().caps_import(_enc(str(path)), _enc(topology), _enc(opts)), str(path))


def import_preview(path: str, **options) -> dict:
    """What import_file would make of frame 0: counts, cell, bond orders and the first lines."""
    return _json_call(library().caps_import_preview, _enc(str(path)), _enc(json.dumps(options)))


# ---------------------------------------------------------------------------------------------------------------- builders

class build:
    """The builders: each returns a new Document."""

    @staticmethod
    def smiles(smiles: str, forcefield: str = "uff", conformers: int = 1, seed: int = 1) -> Document:
        o = _BuildOpts(conformers, seed)
        rep = _report()
        ff = None if forcefield in ("", None) else ("uff" if forcefield == "uff" else _forcefield_path(forcefield))
        d = Document(library().caps_build_smiles(_enc(smiles), _enc(ff), C.byref(o), rep, len(rep)), smiles)
        d.report = rep.value.decode()
        return d

    @staticmethod
    def peptide(sequence: str, structure: str = "", n_term: str = "NH3+", c_term: str = "COO-", ph: float = 7.0, cleanup: bool = True) -> Document:
        rep = _report()
        spec = {"sequence": sequence, "structure": structure, "n_term": n_term, "c_term": c_term, "ph": ph, "cleanup": cleanup}
        d = Document(library().caps_peptide_build(_enc(json.dumps(spec)), rep, len(rep)), "peptide")
        d.report = rep.value.decode()
        return d

    @staticmethod
    def crystal(space_group: str, a: float, b: float, c: float, sites, alpha: float = 90, beta: float = 90, gamma: float = 90,
                supercell=(1, 1, 1), primitive: bool = False) -> Document:
        """sites: [("Ti1", "Ti", x, y, z), …] (fractional)."""
        spec = {"space_group": space_group, "a": a, "b": b, "c": c, "alpha": alpha, "beta": beta, "gamma": gamma,
                "sites": [{"label": s[0], "element": s[1], "x": s[2], "y": s[3], "z": s[4]} for s in sites],
                "supercell": list(supercell), "primitive": primitive}
        rep = _report()
        d = Document(library().caps_crystal_build(_enc(json.dumps(spec)), rep, len(rep)), space_group)
        d.report = rep.value.decode()
        return d

    @staticmethod
    def nano(**options) -> Document:
        """kind="tube" (n, m, length, periodic) | "sheet" (lx, ly, layers) | "particle" (crystal CIF, shape, radius)."""
        rep = _report()
        d = Document(library().caps_nano_build(_enc(json.dumps(options)), rep, len(rep)), options.get("kind", "nano"))
        d.report = rep.value.decode()
        return d

    @staticmethod
    def solvate(solute: Optional[Document] = None, **options) -> Document:
        """Solvent and ions around a solute (or a box of solvent): shape, edge, padding, solvent, water_model, ion_mode,
        salt, concentration …"""
        rep = _report()
        h = solute._h if solute is not None else None
        d = Document(library().caps_solvate(h, _enc(json.dumps(options)), None, None, rep, len(rep)), "solvated")
        d.report = rep.value.decode()
        return d


def space_groups() -> list:
    return _json_call(library().caps_space_groups)
