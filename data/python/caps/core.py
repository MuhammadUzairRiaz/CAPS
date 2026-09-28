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
import builtins as _builtins
import json
import os
import os as _os
from os.path import abspath as _os_path_abspath, exists as _os_path_exists, join as _os_path_join
import re
import sys
from pathlib import Path
from typing import Any, Iterable, Optional

__all__ = ["Document", "Provenance", "open", "build", "run", "polymer", "library", "abi_version", "CapsError"]


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
                ("ambient_occlusion", C.c_int32), ("lod_near", C.c_double), ("lod_far", C.c_double)]


class _RelaxOpts(C.Structure):
    _fields_ = [("method", C.c_int32), ("ftol", C.c_double), ("max_iterations", C.c_int32), ("target_density", C.c_double),
                ("compress_step", C.c_double), ("pushoff", C.c_int32), ("relax_box", C.c_int32), ("pressure", C.c_double),
                ("cutoff", C.c_double), ("coulomb", C.c_int32), ("threads", C.c_int32), ("box_anisotropic", C.c_int32), ("box_axes", C.c_int32),
                ("pushoff_ramp_ps", C.c_double), ("pushoff_cap", C.c_double), ("pushoff_temperature", C.c_double),
                ("etol", C.c_double), ("pressure_tol", C.c_double)]


class _MdOpts(C.Structure):
    _fields_ = [("dt", C.c_double), ("steps", C.c_int64), ("temperature", C.c_double), ("thermostat", C.c_int32), ("tau_t", C.c_double),
                ("barostat", C.c_int32), ("pressure", C.c_double), ("tau_p", C.c_double), ("new_velocities", C.c_int32),
                ("seed", C.c_uint64), ("thermo_every", C.c_int32), ("frame_every", C.c_int32), ("cutoff", C.c_double),
                ("coulomb", C.c_int32), ("tail", C.c_int32), ("threads", C.c_int32), ("respa", C.c_int32), ("constraints", C.c_int32),
                ("step_offset", C.c_int64), ("checkpoint_every", C.c_int64), ("constraint_algorithm", C.c_int32),
                ("box_anisotropic", C.c_int32), ("box_axes", C.c_int32)]


class _ReactOpts(C.Structure):
    _fields_ = [("seed", C.c_uint64), ("max_cycles", C.c_int32), ("max_per_cycle", C.c_int32), ("target_conversion", C.c_double),
                ("capture", C.c_double), ("relax", C.c_int32), ("relax_iterations", C.c_int32), ("md_ps", C.c_double),
                ("temperature", C.c_double), ("cutoff", C.c_double), ("coulomb", C.c_int32), ("during_md", C.c_int32)]


class _BuildOpts(C.Structure):
    _fields_ = [("conformers", C.c_int32), ("seed", C.c_uint64)]


class _AnalyzeOpts(C.Structure):
    _fields_ = [("first", C.c_int64), ("last", C.c_int64), ("stride", C.c_int64), ("frame_ps", C.c_double), ("timestep_fs", C.c_double),
                ("blocks", C.c_int32), ("elem_a", C.c_int32), ("elem_b", C.c_int32), ("inter_only", C.c_int32),
                ("rdf_rmax", C.c_double), ("rdf_dr", C.c_double), ("qmax", C.c_double), ("dq", C.c_double), ("q_direct", C.c_double),
                ("fit_from", C.c_double), ("fit_to", C.c_double), ("probe", C.c_double), ("grid", C.c_double), ("cutoff", C.c_double),
                ("threads", C.c_int32), ("deuterate", C.c_int32), ("group", C.c_char_p)]


class _MechOpts(C.Structure):
    _fields_ = [("configurations", C.c_int32), ("strain", C.c_double), ("temperature", C.c_double), ("axis", C.c_int32),
                ("rate", C.c_double), ("max_strain", C.c_double), ("fit_strain", C.c_double), ("lateral_fixed", C.c_int32),
                ("t_start", C.c_double), ("t_end", C.c_double), ("t_step", C.c_double), ("ps_per_step", C.c_double),
                ("dt", C.c_double), ("pressure", C.c_double), ("seed", C.c_uint64), ("run_ps", C.c_double), ("equilibrate_ps", C.c_double)]


_RecipeProgress = C.CFUNCTYPE(C.c_int32, C.c_int32, C.c_int32, C.c_char_p, C.c_char_p, C.c_char_p, C.c_double, C.c_void_p)


def _declare(L: C.CDLL) -> None:
    P, S, I, D, B = C.c_void_p, C.c_char_p, C.c_int32, C.c_double, C.c_char_p
    sig = {
        "caps_abi_version": ([], I), "caps_last_error": ([], S), "caps_set_restraints": ([P, C.c_char_p], I),
        "caps_pack": ([S, S, I, P, P, B, I], P),
        "caps_set_held_molecule": ([P, C.c_int64], None), "caps_set_fixed_atoms": ([P, C.POINTER(C.c_int32), I], I),
        "caps_chi_md": ([C.c_char_p, P, P, B, I], I), "caps_chi_contacts": ([C.c_char_p, P, P, B, I], I),
        "caps_open": ([S, S], P), "caps_close": ([P], None), "caps_import": ([S, S, S], P), "caps_provenance": ([P, B, I], I), "caps_provenance_file": ([S, B, I], I), "caps_provenance_compare": ([S, S, B, I], I), "caps_provenance_bibtex": ([S, B, I], I), "caps_methods_text": ([S, S, B, I], I), "caps_import_preview": ([S, S, B, I], I),
        "caps_summary_get": ([P, C.POINTER(_Summary)], I), "caps_set_frame": ([P, C.c_int64], I),
        "caps_atom": ([P, I, C.POINTER(_Atom)], I), "caps_save": ([P, S], I), "caps_save_trajectory": ([P, S], I), "caps_gromacs": ([P, S, B, I], I), "caps_export_engines": ([P, S, S, B, I], I),
        "caps_export_png": ([P, C.POINTER(_Camera), C.POINTER(_RenderOpts), S], I),
        "caps_relax": ([P, C.POINTER(_RelaxOpts), P, P, B, I], I), "caps_md": ([P, C.POINTER(_MdOpts), P, P, B, I], I),
        "caps_field_assign": ([P, S, S, I], I), "caps_field_assign_groups": ([P, S], I), "caps_field_report": ([P, B, I], I), "caps_field_import": ([P, S], I), "caps_field_import_ex": ([P, S, S], I),
        "caps_build_smiles": ([S, S, C.POINTER(_BuildOpts), B, I], P),
        "caps_build_beads": ([S, S, C.c_uint64, B, I], P), "caps_bead_templates": ([S, B, I], I),
        "caps_peptide_build": ([S, B, I], P), "caps_crystal_build": ([S, B, I], P), "caps_nano_build": ([S, B, I], P),
        "caps_solvate": ([P, S, P, P, B, I], P),
        "caps_edit": ([P, S, B, I], I), "caps_undo": ([P, I], I), "caps_select": ([P, S, B, I], I), "caps_selection": ([P, B, I], I),
        "caps_tacticity": ([P, B, I], I), "caps_interactions": ([P, S, B, I], I), "caps_torsion_scan": ([P, S, P, P, B, I], I),
        "caps_trajectory_series": ([P, S, P, P, B, I], I), "caps_file_checks": ([P, B, I], I),
        "caps_space_groups": ([B, I], I), "caps_crystal_find_symmetry": ([S, S, D, B, I], I),
        "caps_recipe_run": ([S, S, _RecipeProgress, P, B, I], P), "caps_scene_json": ([P, S, B, I], I),
        "caps_analyze": ([P, S, C.POINTER(_AnalyzeOpts), P, P], I), "caps_analyze_ex": ([P, S, C.POINTER(_AnalyzeOpts), C.POINTER(_MechOpts), P, P], I),
        "caps_analyze_report": ([P, B, I], I),
        "caps_hydrogen_plan": ([P, B, I], I), "caps_resolution_summary": ([P, S, B, I], I), "caps_resolution_convert": ([P, S, B, I], P),
        "caps_chain_lengths": ([S, B, I], I), "caps_copolymer": ([S, B, I], I), "caps_stereo": ([S, B, I], I),
        "caps_react": ([P, S, C.POINTER(_ReactOpts), P, P, B, I], I), "caps_reaction_template": ([S, B, I], I),
        "caps_insert_molecules": ([P, S, I, D, C.c_uint64, B, I], I),
        "caps_blend_phase": ([S, B, I], I), "caps_solvent_chi": ([S, B, I], I), "caps_ewald_params": ([P, S, B, I], I),
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

    def assign(self, forcefield: str = "uff", charges: str = "auto") -> dict:
        """Types every atom and looks up every parameter: a force-field id from the library (gaff2, opls2005 …), a path
        to a caps-forcefield JSON, or "uff". charges: auto (the force field's, else Gasteiger) | forcefield | gasteiger | keep | qeq | increments (bond increments by the types' numbers: OPLS-AA 2024 with OPLS 2005's)."""
        path = _forcefield_path(forcefield)
        code = {"forcefield": 0, "gasteiger": 1, "keep": 2, "qeq": 3, "auto": 4, "increments": 5}[charges]
        rc = library().caps_field_assign(self._doc._h, _enc(path), None, code)
        if rc < 0:
            raise _error()
        rep = _json_call(library().caps_field_report, self._doc._h)
        rep["complete"] = rc == 0
        return rep

    def assign_groups(self, groups, eps_rule: str = "geometric", sigma_rule: str = "arithmetic", scaling14: str = "refuse",
                      cross96: str = "refuse", pairs: Optional[list] = None) -> dict:
        """A force field per group of molecules: groups = [{"name": "filler", "molecules": "1", "forcefield": "iff-cvff"},
        {"name": "matrix", "molecules": "rest", "forcefield": "gaff2", "charges": "auto"}]. A crystal group may instead take
        a literature many-body potential that LAMMPS reads from its file: {"name": "Si", "molecules": "1", "potential":
        {"style": "tersoff", "file": "Si.tersoff", "units": "metal"}} (units only when the file does not say them); its
        atoms get one type per element, standard masses, no charge, UFF Lennard-Jones for the cross pairs, and CAPS runs
        need them held (LAMMPS evaluates the potential). Between groups the Lennard-Jones
        pairs follow eps_rule (geometric | arithmetic) and sigma_rule (arithmetic | geometric | sixthpower), or pairs =
        [{"a": type, "b": type, "eps": kcal/mol, "sigma": Å}]. Different 1-4 scalings are refused unless scaling14="first";
        9-6 with 12-6 unless cross96="rmin" (the 9-6 sites keep ε and r_min)."""
        codes = {"forcefield": 0, "gasteiger": 1, "keep": 2, "qeq": 3, "auto": 4, "increments": 5}
        gs = []
        for g in groups:
            g = dict(g)
            if "potential" in g:   # a literature many-body potential: {"style": "tersoff", "file": PATH, "units": "metal"}
                g["potential"] = dict(g["potential"], file=os.path.abspath(os.path.expanduser(g["potential"]["file"])))
            else:
                g["forcefield"] = _forcefield_path(g["forcefield"])
                g["charges"] = codes[g.get("charges", "auto")]
            gs.append(g)
        spec = {"groups": gs, "eps_rule": eps_rule, "sigma_rule": sigma_rule, "scaling14": scaling14, "cross96": cross96, "pairs": pairs or []}
        rc = library().caps_field_assign_groups(self._doc._h, _enc(json.dumps(spec)))
        if rc < 0:
            raise _error()
        rep = _json_call(library().caps_field_report, self._doc._h)
        rep["complete"] = rc == 0
        return rep

    def import_params(self, path: str, fill_gaps: bool = False) -> dict:
        """Adds parameters from a file over the assigned force field: caps-forcefield JSON, moltemplate .lt, AMBER frcmod
        or the [ *types ] sections of a GROMACS .itp / .top. fill_gaps: used only where the force field defines nothing
        (borrowed terms are listed in the report as filled_terms)."""
        opts = json.dumps({"mode": "fill" if fill_gaps else "override"})
        rc = library().caps_field_import_ex(self._doc._h, _enc(str(path)), _enc(opts))
        if rc < 0:
            raise _error()
        rep = _json_call(library().caps_field_report, self._doc._h)
        rep["complete"] = rc == 0
        return rep


def _forcefield_path(ff: str) -> str:
    if ff.lower() == "uff" or os.path.exists(ff):
        return ff if ff.lower() != "uff" else "uff"
    lib = Path(__file__).resolve().parents[1].parent / "forcefields"
    if ff.endswith("-dlfield"):   # library ids before the force fields got CAPS's own names
        ff = ff[: -len("-dlfield")]
    aliases = {"gaff2": "gaff-amber25", "gaff": "gaff-amber16", "opls": "opls2005", "pcff": "pcff",
               "compass": "compass", "cvff": "cvff", "amber": "amber", "charmm": "charmm"}
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
              box: bool = False, pressure: float = 1.0, cutoff: float = 10.0, coulomb: bool = True, threads: int = 0,
              restraints: Optional[list] = None, box_axes: str = "", pushoff_md_ps: float = 0.0, pushoff_cap: float = 0.0,
              pushoff_temperature: float = 0.0, etol: float = 0.0, pressure_tol: float = 0.0) -> int:
        """Minimises the current frame with the Field assignment (else the built-in GAFF): 0 converged, 1 not quite.
        restraints: [(i, j, r0), …] or [(i, j, r0, k), …] — k (r − r0)² between atoms i and j (indices from 0, Å,
        k kcal/mol/Å², default 10); dihedral ones as {"i", "j", "k", "l", "phi0", "kphi"} dicts; they stay set for later
        relaxations until relax(restraints=[]) clears them. box_axes: "z", "xy" … relaxes the box axis by axis (only
        those axes move); "" keeps the isotropic box relaxation (box=True). pushoff_md_ps > 0: the push-off by NVT MD
        first (Auhl et al. 2003), the LJ force cap raised to pushoff_cap (default 500 kcal/mol/Å) over that time."""
        if restraints is not None:
            rs = [dict(r) if isinstance(r, dict) else {"i": int(r[0]), "j": int(r[1]), "r0": float(r[2]), "k": float(r[3]) if len(r) > 3 else 10.0}
                  for r in restraints]
            if library().caps_set_restraints(self._h, _enc(json.dumps(rs))) < 0:
                raise _error()
        o = _RelaxOpts({"sd": 0, "cg": 1, "lbfgs": 2, "fire": 3}[method], ftol, max_iterations, density, 0.06, int(pushoff), int(box),
                       pressure, cutoff, int(coulomb), threads, int(bool(box_axes)),
                       sum({"x": 1, "y": 2, "z": 4}[a] for a in set(box_axes)), pushoff_md_ps, pushoff_cap, pushoff_temperature)
        if box_axes:
            o.relax_box = 1
        o.etol, o.pressure_tol = etol, pressure_tol   # 0: 1e-8 and 100 atm
        rep = _report()
        rc = library().caps_relax(self._h, C.byref(o), None, None, rep, len(rep))
        if rc < 0:
            raise _error()
        self.report = rep.value.decode()
        return rc

    def hold(self, molecule: int = 0, atoms: Optional[list] = None) -> None:
        """Holds atoms in place in relax, md and equilibration (no force, no motion; a freeze group in GROMACS files):
        molecule > 0 holds that molecule (an interface's surface is 1), atoms (indices from 0) any others; hold() frees them."""
        library().caps_set_held_molecule(self._h, int(molecule))
        idx = [int(a) for a in (atoms or [])]
        arr = (C.c_int32 * max(1, len(idx)))(*idx)
        if library().caps_set_fixed_atoms(self._h, arr, len(idx)) < 0:
            raise _error()

    def md(self, steps: int = 10000, dt: float = 1.0, temperature: float = 300.0, thermostat: str = "bussi", barostat: str = "none",
           pressure: float = 1.0, seed: int = 1, frame_every: int = 1000, thermo_every: int = 100, cutoff: float = 10.0,
           respa: int = 1, constraints: str = "none", constraint_solver: str = "shake", couple_axes: str = "") -> str:
        """Molecular dynamics from the current frame; the frames recorded become the document's frames. respa > 1: r-RESPA,
        the bonded forces every dt / respa (e.g. dt=2, respa=4 with hydrogens). constraints "h-bonds" (bonds to hydrogen,
        rigid water) or "all-bonds": SHAKE/RATTLE, for dt=2 (the alternative to respa). thermostat "nose-hoover" with
        barostat "mtk" runs as LAMMPS's fix nvt / fix npt iso. constraint_solver "lincs" puts the positions back on the
        constraints with LINCS (as GROMACS), "shake" with SHAKE (as LAMMPS); the same result to the tolerance."""
        o = _MdOpts(dt, steps, temperature, {"none": 0, "bussi": 1, "langevin": 2, "nose-hoover": 3}[thermostat], 100.0,
                    {"none": 0, "crescale": 1, "berendsen": 2, "mtk": 3}[barostat], pressure, 1000.0, 0, seed, thermo_every, frame_every,
                    cutoff, 1, 1, 0, respa, {"none": 0, "h-bonds": 1, "all-bonds": 2}[constraints], 0, 0,
                    {"shake": 0, "lincs": 1}[constraint_solver])
        if couple_axes:   # "z", "xy" …: those axes coupled on their own to the pressure (Berendsen), the others fixed
            o.box_anisotropic, o.box_axes = 1, sum({"x": 1, "y": 2, "z": 4}[a] for a in set(couple_axes))
        rep = _report()
        if library().caps_md(self._h, C.byref(o), None, None, rep, len(rep)) < 0:
            raise _error()
        self.report = rep.value.decode()
        return self.report

    def insert(self, smiles: str, count: int, tolerance: float = 2.0, seed: int = 1) -> str:
        """Inserts count copies of a molecule (SMILES; hydrogens added, UFF-cleaned) into the free space of the current
        frame — curatives before a cure, e.g. insert("SS", 40) for the sulfur_allylic template. The document becomes that
        frame and its force-field assignment is cleared (assign again after)."""
        rep = _report()
        if library().caps_insert_molecules(self._h, _enc(smiles), int(count), float(tolerance), int(seed), rep, len(rep)) < 0:
            raise _error()
        self.report = rep.value.decode()
        return self.report

    def react(self, templates, cycles: int = 50, per_cycle: int = 5, target: float = 1.0, capture: float = 0.0, relax: bool = True,
              relax_iterations: int = 500, md_ps: float = 0.0, temperature: float = 500.0, cutoff: float = 10.0, seed: int = 1,
              during_md: bool = False) -> str:
        """Crosslinks the current frame cycle by cycle (Polymatic-style; REACTER-style capture and probability) with
        reaction templates: built-in names ("cc_crosslink", "sulfur_allylic", "peroxide_allylic", "polysulfide_allylic",
        "epoxy_amine_primary", …; see reaction_templates()) or template text, one or a list. target is the conversion to
        stop at (0 … 1); relax minimises after each cycle, md_ps runs NVT at temperature between cycles. The force-field
        assignment is cleared (the topology changed): assign again to type the network. Returns the report."""
        names = [templates] if isinstance(templates, str) else list(templates)
        text = "\n".join(reaction_template(t) if "\n" not in t and t.strip() in reaction_templates() else t for t in names)
        o = _ReactOpts(int(seed), int(cycles), int(per_cycle), float(target), float(capture), int(relax), int(relax_iterations),
                       float(md_ps), float(temperature), float(cutoff), 1, int(during_md))
        rep = _report()
        if library().caps_react(self._h, _enc(text), C.byref(o), None, None, rep, len(rep)) < 0:
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

    def query(self, query: str, op: str = "replace") -> int:
        """Selects by the query grammar (design/boards/SmartSelect): smarts "c1ccccc1", element C N O, type c3,
        chain 1-4, index 1-20, ring 5, stereo R|S|*, within 5.0 of <query>, sel, combined with and / or / not / ( ).
        op "preview" counts without selecting. Returns the number of atoms."""
        return self.select("query", query, op=op)

    def selection(self) -> list:
        return [int(i) for i in _json_call(library().caps_selection, self._h)["indices"]]

    def tacticity(self) -> dict:
        return _json_call(library().caps_tacticity, self._h)

    def hydrogen_plan(self) -> dict:
        """What add_h would add, by kind of atom: {rows: [{label, atoms, hydrogens}], heavy, h, add, …} (bond orders
        from the geometry when the structure has no H and no multiple bonds)."""
        return _json_call(library().caps_hydrogen_plan, self._h)

    def add_hydrogens(self) -> int:
        """Adds the hydrogens every atom lacks (valence rules); returns how many."""
        before = self.atoms
        _json_call(library().caps_edit, self._h, _enc('{"op": "add_h"}'))
        return self.atoms - before

    def resolution(self, per_bead: int = 5) -> dict:
        """Sites, hydrogens and mass all-atom, united-atom and coarse-grained (mass conserved)."""
        return _json_call(library().caps_resolution_summary, self._h, _enc(json.dumps({"per_bead": per_bead})))

    def convert(self, to: str = "united-atom", per_bead: int = 5) -> "Document":
        """A new document at another resolution: "united-atom" (H on carbon folded in) or "coarse-grained" (beads of
        per_bead backbone atoms at their centre of mass); this one is unchanged."""
        rep = _report()
        h = library().caps_resolution_convert(self._h, _enc(json.dumps({"to": to, "per_bead": per_bead})), rep, len(rep))
        d = Document(h, f"{self.label} ({to})")
        d.report = rep.value.decode()
        return d

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
    @property
    def provenance(self) -> "Provenance":
        """The steps that produced this structure: doc.provenance() is the manifest (caps-manifest/1.0, written beside
        the file on save); doc.provenance.methods(), .citations(fmt="bibtex"|"text"), .bibtex()."""
        return Provenance(lambda: _json_call(library().caps_provenance, self._h))

    # analysis and viewing
    def analyze(self, properties="density", first: int = 0, last: int = -1, stride: int = 1, blocks: int = 5, threads: int = 0,
                **options) -> list:
        """Properties of the frames (the ids of caps analyze: density, rdf, rg, ree, cn, persistence, msd, diffusion, ced,
        ffv …; tg runs a stepwise cooling of a copy, t_start/t_end/t_step/ps_per_step in options). A list of
        {id, name, value, error, unit, ...}."""
        ids = ",".join(_PROPERTY_ALIASES.get(p, p) for p in ([properties] if isinstance(properties, str) else properties))
        o = _AnalyzeOpts(first, last, stride, 0, 0, blocks, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, threads, int(options.pop("deuterate", 0)),
                         _enc(options.pop("group", "")))   # group: "selection", "molecules:1-4,7", "exclude-held"
        if "tg" in ids.split(","):
            m = _MechOpts()
            for k, v in options.items():
                setattr(m, k, v)
            rc = library().caps_analyze_ex(self._h, _enc(ids), C.byref(o), C.byref(m), None, None)
        else:
            rc = library().caps_analyze(self._h, _enc(ids), C.byref(o), None, None)
        if rc < 0:
            raise _error()
        return _json_call(library().caps_analyze_report, self._h)["properties"]

    def view(self, style: str = "ball-and-stick", width: int = 640, height: int = 400, hydrogens: bool = True, max_atoms: int = 60000,
             cell: Optional[bool] = None):
        """An interactive view for notebooks (drag to rotate, wheel to zoom, double-click to reset); a PNG where HTML is not
        shown. style: ball-and-stick, space-filling, sticks, no-hydrogens, lines; cell: draw the periodic cell (default:
        when there is more than one molecule)."""
        from .view import View
        return View(self, style=style, width=width, height=height, hydrogens=hydrogens, max_atoms=max_atoms, cell=cell)

    def scene(self, hydrogens: bool = True, max_atoms: int = 60000) -> dict:
        """The current frame for a viewer: z, xyz (flat), bonds (flat pairs), colours, radii, cell."""
        return _json_call(library().caps_scene_json, self._h, _enc(json.dumps({"hydrogens": hydrogens, "max_atoms": max_atoms})))

    def save(self, path: str) -> None:
        """Writes the current frame: .data (LAMMPS, with force-field sections when assigned), .pdb, .xyz, .mol2, .gro …"""
        if library().caps_save(self._h, _enc(str(path))) != 0:
            raise _error()

    def save_gromacs(self, stem: str) -> str:
        """Writes stem.top, stem.gro and stem.mdp for GROMACS with the same force field (a single-point run; energies and
        forces checked against GROMACS by bench/ff/check_gromacs.py). Returns the .mdp non-bonded settings, with
        "; note:" lines where GROMACS computes differently (DSF becomes PME; no cell)."""
        need = library().caps_gromacs(self._h, None, None, 0)
        if need < 0:
            raise _error()
        buf = C.create_string_buffer(need)
        if library().caps_gromacs(self._h, _enc(str(stem)), buf, need) < 0:
            raise _error()
        return buf.value.decode()

    def export_engines(self, folder: str, stem: str = "system", lammps: bool = True, gromacs: bool = True, run: str = "check", **opts) -> dict:
        """The simulation files for LAMMPS (stem.data, stem.in with every pair_coeff and the run) and GROMACS (stem.top,
        stem.itp, stem.gro, stem.mdp; stem_em.mdp when a run minimises first) from the assigned force field, which must be
        complete. run: check | none | minimize | nvt | npt; opts: minimize_first, temperature (K), pressure (atm), dt (fs),
        steps, thermo_every, dump_every, seed. Returns {folder, files, notes, checks}; raises CapsError when refused."""
        o = dict(opts, stem=stem, lammps=lammps, gromacs=gromacs, run=run)
        r = _json_call(library().caps_export_engines, self._h, _enc(str(folder)), _enc(json.dumps(o)))
        if not r.get("ok"):
            raise CapsError(r.get("error", "export failed"))
        return r

    def save_trajectory(self, path: str) -> None:
        if library().caps_save_trajectory(self._h, _enc(str(path))) != 0:
            raise _error()

    def render(self, path: str, width: int = 1280, height: int = 800, background: str = "white", style: str = "ball_and_stick",
               colour: str = "molecule", yaw: float = 0.55, pitch: float = 0.40, zoom: float = 1.0) -> None:
        """A PNG of the current frame."""
        style = style.replace("-", "_")
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


def current() -> Document:
    """The structure a CAPS Studio macro runs on (Macro › Target: the open structure): the Studio saves it and names the
    file in CAPS_DOC. Outside the Studio, set CAPS_DOC to a structure file."""
    path = os.environ.get("CAPS_DOC", "")
    if not path:
        raise CapsError("no current structure: run from CAPS Studio with Target: the open structure, or set CAPS_DOC")
    return open(path, os.environ.get("CAPS_DOC_TOPOLOGY") or None)


def hand_back(doc: "Document") -> str:
    """Gives the macro's result back to CAPS Studio (it opens it when the macro ends): saved to CAPS_OUT, else to
    result.data beside the script. Returns the path."""
    path = os.environ.get("CAPS_OUT") or "result.data"
    doc.save(path)
    return path


def import_file(path: str, bonds: str = "perceive", tolerance: float = 0.45, bond_orders: bool = True, split: bool = True,
                unwrap: bool = True, use_cell: bool = True, topology: Optional[str] = None) -> Document:
    """Opens a file that has no topology (XYZ, PDB, CIF …) with the Import choices: bonds "perceive" (covalent radii +
    tolerance Å), "file" or "none"; bond orders and aromaticity; molecules by connectivity; unwrapping; the cell."""
    opts = json.dumps({"bonds": bonds, "tolerance": tolerance, "bond_orders": bond_orders, "split": split, "unwrap": unwrap, "use_cell": use_cell})
    return Document(library().caps_import(_enc(str(path)), _enc(topology), _enc(opts)), str(path))


def import_preview(path: str, **options) -> dict:
    """What import_file would make of frame 0: counts, cell, bond orders and the first lines."""
    return _json_call(library().caps_import_preview, _enc(str(path)), _enc(json.dumps(options)))


def provenance_file(path: str) -> dict:
    """The provenance saved beside a file (<file>.provenance.json); {"ok": False, ...} without one."""
    return _json_call(library().caps_provenance_file, _enc(str(path)))


def compare_provenance(a: dict, b: dict) -> dict:
    """Two manifests step by step: the parameters and seeds that differ."""
    return _json_call(library().caps_provenance_compare, _enc(json.dumps(a)), _enc(json.dumps(b)))


def methods(manifest: dict, replicas: Optional[list] = None) -> dict:
    """A methods paragraph for a paper from a manifest: {"text": ..., "refs": [...]} with numbered references."""
    return _json_call(library().caps_methods_text, _enc(json.dumps(manifest)), _enc(json.dumps(replicas) if replicas else None))


def bibtex(manifest: dict) -> str:
    """BibTeX of every method a manifest cites."""
    cap = 1 << 16
    for _ in range(3):
        buf = C.create_string_buffer(cap)
        n = library().caps_provenance_bibtex(_enc(json.dumps(manifest)), buf, cap)
        if n < 0:
            raise _error()
        if n <= cap:
            return buf.value.decode()
        cap = n + 1
    raise CapsError("reply too large")


_PROPERTY_ALIASES = {"c_inf": "cn", "cinf": "cn", "C_inf": "cn", "end_to_end": "ree", "free_volume": "ffv"}


class Provenance:
    """How a structure was made. Call it for the manifest; methods(), citations(), bibtex() for a paper."""

    def __init__(self, source):
        self._source = source

    def __call__(self) -> dict:
        return self._source()

    @property
    def manifest(self) -> dict:
        return self._source()

    @property
    def steps(self) -> list:
        return self._source().get("steps", [])

    def methods(self, replicas: Optional[list] = None) -> str:
        """The methods paragraph with numbered references."""
        return methods(self._source(), replicas)["text"]

    def bibtex(self) -> str:
        return bibtex(self._source())

    def citations(self, fmt: str = "text") -> list:
        """Every method cited: fmt="bibtex" gives one BibTeX entry per item, "text" one formatted reference."""
        if fmt == "bibtex":
            text = self.bibtex().strip()
            return [e.strip() for e in re.split(r"\n(?=@)", text) if e.strip()] if text else []
        if fmt == "text":
            return list(methods(self._source())["refs"])
        raise ValueError('fmt: "bibtex" or "text"')

    def __repr__(self) -> str:
        m = self._source()
        return f"<caps.Provenance · {len(m.get('steps', []))} steps>"

    def _repr_html_(self) -> str:
        import html as _h
        rows = "".join(f"<tr><td style='padding:2px 10px 2px 0;color:#666'>{i + 1}</td><td style='padding:2px 10px 2px 0'><code>{_h.escape(s.get('engine', ''))}</code></td>"
                       f"<td style='padding:2px 0'>{_h.escape(s.get('summary', ''))}</td></tr>" for i, s in enumerate(self.steps))
        return f"<table style='font-size:12px;border-collapse:collapse'>{rows}</table>"


# ---------------------------------------------------------------------------------------------------------------- recipes

def _data_dir() -> Path:
    return Path(__file__).resolve().parents[1].parent


def run(recipe, out_dir: str = ".", seed: Optional[int] = None, threads: int = 0, progress=None, base_dir: Optional[str] = None) -> Document:
    """Runs a recipe (a dict, YAML/JSON text, or a path to one) as caps run does: build → type → grow → relax → md →
    equilibrate → analyze → export. Returns the structure; .properties and .files hold what the recipe analysed and
    wrote. progress(event) gets {stage, stages, name, status, detail, fraction}; return False to stop."""
    label = "recipe"
    if isinstance(recipe, dict):
        text = json.dumps(recipe)
        label = recipe.get("name", label)
    else:
        text = str(recipe)
        if "\n" not in text and os.path.exists(text):
            base_dir = base_dir or str(Path(text).resolve().parent)
            label = Path(text).stem
            text = Path(text).read_text()
    opts = {"base_dir": base_dir or os.getcwd(), "out_dir": str(out_dir), "forcefield_dir": str(_data_dir() / "forcefields"),
            "seed": -1 if seed is None else int(seed), "threads": threads}

    def _cb(stage, stages, name, status, detail, fraction, _user):
        if progress is None:
            return 1
        try:
            r = progress({"stage": stage, "stages": stages, "name": name.decode(), "status": status.decode(), "detail": detail.decode(), "fraction": fraction})
        except Exception:
            return 0
        return 0 if r is False else 1

    cb = _RecipeProgress(_cb)
    buf = C.create_string_buffer(1 << 20)
    h = library().caps_recipe_run(_enc(text), _enc(json.dumps(opts)), cb, None, buf, len(buf))
    rep = json.loads(buf.value.decode() or "{}")
    if not h:
        err = CapsError(rep.get("error") or _error().args[0])
        err.exit_code = rep.get("exit", 4)
        raise err
    d = Document(h, label)
    d.properties = rep.get("properties", [])
    d.files = rep.get("files", [])
    d.forcefield = rep.get("forcefield", "")
    return d


def polymer(smiles, dp: int = 20, chains: int = 1, tacticity: str = "atactic", seed: int = 1, density: Optional[float] = None,
            forcefield: Optional[str] = None, relax: bool = False, sequence: str = "homopolymer", trials: int = 120,
            blocks: Optional[list] = None, weights: Optional[list] = None, pattern: str = "", r1: Optional[float] = None,
            r2: Optional[float] = None, pm: Optional[float] = None, p_mr: Optional[float] = None, p_rm: Optional[float] = None,
            lengths: Optional[dict] = None, chain_dp: Optional[list] = None, architecture: str = "linear",
            arms: Optional[int] = None, arm_dp: Optional[int] = None, spacing: Optional[int] = None,
            branch_probability: Optional[float] = None, generations: Optional[int] = None, region: Optional[dict] = None,
            method: str = "trials", method_temperature: float = 450.0, orientation: Optional[dict] = None, lookahead: int = 1) -> Document:
    """Chains of a repeat unit (SMILES with two * points, or a list of them for copolymers — sequence alternating, block
    with blocks=[…], random with weights=[…], gradient, pattern="AAB", terminal with r1, r2 and weights=[f1, f2]) grown
    in a periodic cell: one chain in a roomy cell by default (0.1 g/cm³), a melt with chains=… density=…. Atactic
    chains: pm (Bernoulli) or p_mr, p_rm (first-order Markov). Polydisperse: lengths={"distribution": "schulz-zimm",
    "nn": 40, "pdi": 1.1, "seed": 1} or chain_dp=[…]; the provenance records the sample drawn. Branched molecules:
    architecture="star" with arms=3|4 (each arm dp units on one core carbon), "comb" with arm_dp and spacing, or
    "branched" with arm_dp and branch_probability, "dendrimer" with arms, arm_dp and generations (each end splits in two
    branches of arm_dp units per generation); chains counts molecules. forcefield types it (default: the built-in
    GAFF for C and H, else UFF); relax=True minimises. region={"shape": "slab", "thickness": 30, "vacuum": 30} grows a
    film, {"shape": "cylinder" | "around_cylinder", "radius": 10} chains in or around a cylinder along z. method: "trials"
    (the roomiest of k trials), "rosenbluth" (a trial drawn by its Boltzmann weight: soft spheres and butane torsions) or
    "rosenbluth_lj" (the same with UFF Lennard-Jones), at method_temperature (K). orientation={"axis": "z", "strength": 4}
    grows oriented chains (an aligning field −s P₂ in kT on each unit's backbone chord; the report gives ⟨P₂⟩)."""
    units = [smiles] if isinstance(smiles, str) else list(smiles)
    r = {"recipe": 1, "name": "polymer",
         "build": {"polymer": {"units": units, "dp": dp, "chains": chains, "tacticity": tacticity, "sequence": sequence}},
         "grow": {"density": density if density is not None else (0.1 if chains == 1 else 0.5), "seed": seed, "trials": trials}}
    for k, v in (("blocks", blocks), ("weights", weights), ("pattern", pattern), ("lengths", lengths), ("chain_dp", chain_dp)):
        if v:
            r["build"]["polymer"][k] = v
    for k, v in (("r1", r1), ("r2", r2), ("pm", pm), ("p_mr", p_mr), ("p_rm", p_rm), ("arms", arms), ("arm_dp", arm_dp),
                 ("spacing", spacing), ("branch_probability", branch_probability), ("generations", generations)):
        if v is not None:
            r["build"]["polymer"][k] = v
    if architecture != "linear":
        r["build"]["polymer"]["architecture"] = architecture
    if region:
        r["grow"]["region"] = region
    if orientation:
        r["grow"]["orientation"] = orientation
    if lookahead > 1:   # bonds ahead that need room (fewer dead ends in dense cells)
        r["grow"]["lookahead"] = lookahead
    if method != "trials":
        r["grow"]["method"] = method
        r["grow"]["temperature"] = method_temperature
    if forcefield or relax:
        r["type"] = {"forcefield": forcefield or "default"}
    if relax:
        r["relax"] = {"method": "lbfgs", "fmax": 1.0}
    d = run(r)
    d.label = units[0] if len(units) == 1 else "copolymer"
    return d


def pack(molecules=None, box=30.0, tolerance: float = 2.0, seed: int = 1, density: Optional[float] = None, inp: Optional[str] = None,
         forcefield: Optional[str] = None, relax: bool = False, base_dir: Optional[str] = None) -> Document:
    """Molecules packed into a periodic box with no two atoms of different molecules closer than tolerance (Å), as
    packmol does: molecules = [("CCO", 50), ("mol.pdb", 3), (doc, 10) …] — SMILES, structure files or documents with
    their counts — in box = edge or (x, y, z) Å; density compresses the packed cell to that g/cm³. Or inp = a packmol
    input file. forcefield types the cell, relax minimises it."""
    import tempfile
    tmp = tempfile.mkdtemp(prefix="caps-pack-")
    if inp is not None and not forcefield and not relax:
        # packmol input: its structure files are found in base_dir (else beside the input)
        with _builtins.open(str(inp)) as f:
            text = f.read()
        rep = _report()
        h = library().caps_pack(_enc(text), _enc(base_dir or _os_path_abspath(_os.path.dirname(str(inp)) or ".")), 0, None, None, rep, len(rep))
        if not h:
            raise _error()
        d = Document(h, str(inp))
        d.report = rep.value.decode()
        return d
    if inp is not None:
        build = {"pack": _os_path_abspath(str(inp))}
    else:
        mols = []
        for k, (m, count) in enumerate(molecules or []):
            if isinstance(m, Document):
                path = _os_path_join(tmp, f"molecule_{k}.pdb")
                m.save(path)
                mols.append({"file": path, "count": int(count)})
            elif isinstance(m, str) and _os_path_exists(m):
                mols.append({"file": _os_path_abspath(m), "count": int(count)})
            else:
                mols.append({"smiles": str(m), "count": int(count)})
        p = {"molecules": mols, "tolerance": tolerance, "seed": seed, "box": list(box) if isinstance(box, (list, tuple)) else float(box)}
        if density:
            p["density"] = density
        build = {"pack": p}
    r = {"recipe": 1, "name": "packed", "build": build}
    if forcefield or relax:
        r["type"] = {"forcefield": forcefield or "default"}
    if relax:
        r["relax"] = {"method": "lbfgs", "fmax": 1.0}
    return run(r, out_dir=tmp, base_dir=base_dir)


def chi_by_md(polymer, solvent: Optional[str] = None, polymer_b=None, dp: int = 10, chains: int = 6, temperature: float = 300.0,
              eq_ps: float = 200.0, prod_ps: float = 300.0, seed: int = 1) -> dict:
    """Flory–Huggins χ from the energy of mixing by MD (core chimd.hpp): A alone, B alone and a mixture, each relaxed
    and run NPT, χ = V_ref (φ_A CED_A + φ_B CED_B − CED_mix) / (RT φ_A φ_B). Enthalpic only and noisy: natural rubber
    against itself (χ must be 0) gives −1.2 ± 0.5 after 300 ps per cell, so run long, use larger cells and run the
    self-mixing control (polymer_b = polymer) beside it. polymer / polymer_b: repeat-unit SMILES; solvent: SMILES.
    Takes minutes to hours."""
    spec = lambda s: {"units": [{"name": "A", "smiles": s}], "dp": dp}
    o = {"polymer": spec(polymer), "chains": chains, "temperature": temperature, "eq_ps": eq_ps, "prod_ps": prod_ps, "seed": seed}
    if polymer_b is not None:
        o["polymer_b"] = spec(polymer_b)
        o["chains_b"] = chains
    elif solvent:
        o["solvent"] = solvent
    else:
        raise CapsError("give a solvent SMILES or polymer_b")
    lib = library()
    buf = C.create_string_buffer(1 << 16)
    lib.caps_chi_md(_enc(json.dumps(o)), None, None, buf, len(buf))
    r = json.loads(buf.value.decode())
    if not r.get("ok"):
        raise CapsError(r.get("error", "chi_by_md failed"))
    return r


def chi_by_contacts(a: str, b: str, forcefield: Optional[str] = "gaff2", samples: int = 1000000, pack_trials: int = 5000,
                    temperatures=(250, 275, 300, 325, 350, 375, 400), t: float = 298.15, seed: int = 1) -> dict:
    """Flory–Huggins χ(T) from pair contacts (Fan, Olafson, Blanco & Hsu, Macromolecules 1992; core chipair.hpp): pair
    energies of rigid molecules at van der Waals contact, Boltzmann-averaged, and coordination numbers from packing;
    χ = ½(Z_AB E_AB + Z_BA E_BA − Z_AA E_AA − Z_BB E_BB)/RT, fitted to A + B/T. a, b: SMILES of molecules or repeat
    units (* ends capped with H). forcefield: a library id (gaff2 by default) or a path; None for the built-in GAFF
    subset / UFF. A screen, measured on known cases in the README (wrong for PS/THF and PS/PVME). Seconds."""
    o = {"a": a, "b": b, "forcefield": _forcefield_path(forcefield) if forcefield else "", "samples": samples,
         "pack_trials": pack_trials, "temperatures": list(temperatures), "t": t, "seed": seed}
    buf = C.create_string_buffer(1 << 18)
    library().caps_chi_contacts(_enc(json.dumps(o)), None, None, buf, len(buf))
    r = json.loads(buf.value.decode())
    if not r.get("ok"):
        raise CapsError(r.get("error", "chi_by_contacts failed"))
    return r


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
    def beads(text: str, forcefield: Optional[str] = None, seed: int = 1) -> Document:
        """A coarse-grained molecule: a bead template of the force field (MARTINI's DPPC, W, NA+ …; bead_templates()) or
        bead SMILES ("[Q0+1][Qa-1][Na]…"); one site per bead named by its type, bonds at the force field's lengths."""
        rep = _report()
        ff = _forcefield_path(forcefield) if forcefield else None
        d = Document(library().caps_build_beads(_enc(text), _enc(ff), seed, rep, len(rep)), text)
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
        """kind="tube" (n, m, length, periodic) | "sheet" (lx, ly, layers) | "particle" (crystal CIF, shape, radius; a metal
        particle takes thiolate="C6" | "C12" | "C18" | "MPA" | "MUA" | "MHA" | "*S…" SMILES and thiolate_fraction)."""
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


# ---------------------------------------------------------------------------------------------------------------- polymer statistics

def bead_templates(forcefield: str) -> dict:
    """{name: bead SMILES} for a coarse-grained force field (its sources' molecule templates)."""
    return _json_call(library().caps_bead_templates, _enc(_forcefield_path(forcefield)))


def chain_lengths(distribution: str = "schulz-zimm", nn: float = 40, pdi: float = 1.1, count: int = 20, seed: int = 2026,
                  m0: float = 104.15, best_of: int = 1) -> dict:
    """Chain lengths drawn from a distribution (monodisperse, schulz-zimm, flory, poisson): {lengths, sample {nn, mn,
    mw, pdi, min, max, sum}, target {…}, curve {n, number, weight}}. m0 is the repeat unit's molar mass (g/mol)."""
    return _json_call(library().caps_chain_lengths, _enc(json.dumps({"distribution": distribution, "nn": nn, "pdi": pdi, "count": count,
                                                                     "seed": seed, "m0": m0, "best_of": best_of})))


def copolymer_model(r1: float, r2: float, f1: float, dp: int = 80, seed: int = 1) -> dict:
    """The terminal (Mayo–Lewis) model: {F1, paa, pbb, run_a, run_b, azeotrope, curve {f1, F1}, sequence, chain {…}};
    sequence is what polymer(…, sequence="terminal") grows for its first chain at this seed."""
    return _json_call(library().caps_copolymer, _enc(json.dumps({"r1": r1, "r2": r2, "f1": f1, "dp": dp, "seed": seed})))


def stereo(pm: float = 0.5, p_mr: Optional[float] = None, p_rm: Optional[float] = None, dyads: str = "",
           measured: Optional[list] = None) -> dict:
    """Bernoulli (pm) or first-order Markov (p_mr, p_rm) triads and pentads; with dyads ("mrrm…", e.g. from
    Document.tacticity()) the counted values; with measured (ten pentad fractions, NMR) Bernoulli and Markov fits."""
    q: dict = {"model": "markov", "p_mr": p_mr, "p_rm": p_rm} if p_mr is not None and p_rm is not None else {"model": "bernoulli", "pm": pm}
    if dyads:
        q["dyads"] = dyads
    if measured is not None:
        q["measured"] = list(measured)
    return _json_call(library().caps_stereo, _enc(json.dumps(q)))


def blend_phase(na: float, nb: float, a: float, b: float, t: float = 300.0) -> dict:
    """Flory–Huggins binary blend with χ = a + b/T: {chi_c, phi_c, tc, chi_t, coexist, spinodal, binodal {phi, t},
    spinodal_curve {phi, t}}."""
    return _json_call(library().caps_blend_phase, _enc(json.dumps({"na": na, "nb": nb, "a": a, "b": b, "t": t})))


def solvent_chi(delta_polymer: float, solvents: list, t: float = 298.15) -> list:
    """Hildebrand χ ≈ V(δs − δp)²/RT + 0.34 for [{name, v (cm³/mol), delta (MPa½)}]: [{name, chi, predicted}]. It ignores
    polarity and hydrogen bonding: check it against known behaviour."""
    return _json_call(library().caps_solvent_chi, _enc(json.dumps({"delta_polymer": delta_polymer, "t": t, "solvents": solvents})))["solvents"]


def ewald_params(cutoff: float = 12.0, tolerance: float = 1e-5, spacing: float = 1.2, order: int = 4, edges: Optional[list] = None,
                 doc: Optional["Document"] = None) -> dict:
    """β from erfc(β rc) = tolerance and the PME mesh (FFT sizes with factors 2, 3, 5, 7) for edges or doc's cell."""
    q: dict = {"cutoff": cutoff, "tolerance": tolerance, "spacing": spacing, "order": order}
    if edges is not None:
        q["edges"] = list(edges)
    return _json_call(library().caps_ewald_params, doc._h if doc is not None else None, _enc(json.dumps(q)))


def reaction_templates() -> list:
    """The names of the built-in reaction templates (Document.react)."""
    n = library().caps_reaction_template(b"", None, 0)
    buf = C.create_string_buffer(max(1, n + 1))
    library().caps_reaction_template(b"", buf, len(buf))
    return [x for x in buf.value.decode().split("\n") if x.strip()]


def reaction_template(name: str) -> str:
    """The text of a built-in reaction template (edit it and pass the text to Document.react)."""
    n = library().caps_reaction_template(_enc(name), None, 0)
    if n < 0:
        raise _error()
    buf = C.create_string_buffer(n + 1)
    library().caps_reaction_template(_enc(name), buf, len(buf))
    return buf.value.decode()


def space_groups() -> list:
    return _json_call(library().caps_space_groups)
