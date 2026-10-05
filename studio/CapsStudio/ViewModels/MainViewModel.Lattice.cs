using System;
using System.Globalization;
using System.Linq;
using System.Text.Json.Nodes;

namespace CapsStudio.ViewModels;

/// <summary>Cell tools on the open structure (Build › Crystal tools): find the primitive cell, Niggli-reduce, the
/// conventional cell, a redefined lattice, a vacuum slab, a nanowire (lattice parameters and supercells: the Cell editor). Each is one
/// undoable edit (core lattice.hpp).</summary>
public sealed partial class MainViewModel
{
    private string _cellMsg = "";
    /// <summary>What the last cell tool did (or why it could not).</summary>
    public string CellToolText { get => _cellMsg; private set => Set(ref _cellMsg, value); }

    private bool CellEdit(object spec)
    {
        if (_doc == null) { CellToolText = "Open or build a periodic structure first"; return false; }
        var r = RunEdit(spec);
        CellToolText = r == null ? EditError : r["what"]?.GetValue<string>() ?? "Done";
        return r != null;
    }

    private decimal _cellTol = 0.1m;
    public decimal? CellTolerance { get => _cellTol; set { if (value != null) Set(ref _cellTol, Math.Clamp(value.Value, 0.001m, 1m)); } }

    public void FindPrimitiveCell() => CellEdit(new { op = "find_primitive", tolerance = (double)_cellTol });
    public void NiggliCell() => CellEdit(new { op = "niggli" });
    public void ConventionalCell() => CellEdit(new { op = "conventional", tolerance = (double)_cellTol });

    // the new lattice: column j is the new vector j in the old a, b, c (fractions such as 0.5 allowed)
    private string _redefine = "1 0 0\n0 1 0\n0 0 1";
    /// <summary>The matrix as three rows of three numbers (fractions like 1/2 are read).</summary>
    public string RedefineMatrix { get => _redefine; set => Set(ref _redefine, value ?? ""); }
    public static readonly (string Name, string Matrix)[] RedefinePresets =
    [
        ("Identity", "1 0 0\n0 1 0\n0 0 1"),
        ("F-centred → primitive", "0 1/2 1/2\n1/2 0 1/2\n1/2 1/2 0"),
        ("I-centred → primitive", "-1/2 1/2 1/2\n1/2 -1/2 1/2\n1/2 1/2 -1/2"),
        ("C-centred → primitive", "1/2 -1/2 0\n1/2 1/2 0\n0 0 1"),
        ("Hexagonal → orthohexagonal (2× volume)", "1 1 0\n0 2 0\n0 0 1"),
        ("√2 × √2 rotated (2× volume)", "1 -1 0\n1 1 0\n0 0 1"),
        ("Swap a and c", "0 0 1\n0 -1 0\n1 0 0"),
    ];
    public static string[] RedefinePresetNames => RedefinePresets.Select(p => p.Name).ToArray();
    private int _redefinePreset;
    public int RedefinePreset { get => _redefinePreset; set { if (Set(ref _redefinePreset, Math.Clamp(value, 0, RedefinePresets.Length - 1))) RedefineMatrix = RedefinePresets[_redefinePreset].Matrix; } }
    public void RedefineLattice()
    {
        var m = new double[9];
        var parts = _redefine.Split((char[])[' ', '\t', '\n', '\r', ','], StringSplitOptions.RemoveEmptyEntries);
        if (parts.Length != 9) { CellToolText = $"The matrix needs 9 numbers ({parts.Length} given)"; return; }
        for (var i = 0; i < 9; ++i)
        {
            var p = parts[i];
            double v;
            var slash = p.IndexOf('/');
            if (slash > 0 && double.TryParse(p[..slash], NumberStyles.Float, CultureInfo.InvariantCulture, out var num)
                && double.TryParse(p[(slash + 1)..], NumberStyles.Float, CultureInfo.InvariantCulture, out var den) && den != 0) v = num / den;
            else if (!double.TryParse(p, NumberStyles.Float, CultureInfo.InvariantCulture, out v)) { CellToolText = $"Not a number: {p}"; return; }
            m[i] = v;
        }
        CellEdit(new { op = "redefine_lattice", matrix = m });
    }

    private decimal _slabVac = 15;
    private bool _slabCentre = true;
    public decimal? SlabVacuum { get => _slabVac; set { if (value != null) Set(ref _slabVac, Math.Clamp(value.Value, 0, 500)); } }
    public bool SlabCentre { get => _slabCentre; set => Set(ref _slabCentre, value); }
    public void MakeVacuumSlab() => CellEdit(new { op = "vacuum_slab", vacuum = (double)_slabVac, centre = _slabCentre });

    // a cluster from the periodic structure (core lattice.hpp cut_cluster): whole molecules about a centre, no cell
    private decimal _clusterR = 12;
    private bool _clusterSel, _clusterAny;
    public decimal? ClusterRadius { get => _clusterR; set { if (value != null) Set(ref _clusterR, Math.Clamp(value.Value, 0, 500)); } }
    public bool ClusterOnSelection { get => _clusterSel; set => Set(ref _clusterSel, value); }
    public bool ClusterAnyAtom { get => _clusterAny; set => Set(ref _clusterAny, value); }
    public void CutCluster()
    {
        if (_clusterSel) CellEdit(new { op = "cluster", radius = (double)_clusterR, centre = "selection", any_atom = _clusterAny, atoms = SelectionAtoms() });
        else CellEdit(new { op = "cluster", radius = (double)_clusterR, any_atom = _clusterAny });
    }

    public static readonly string[] WireShapes = ["Cylinder", "Hexagonal prism", "Square prism"];
    private static readonly string[] WireShapeIds = ["cylinder", "hexagonal", "square"];
    private string _wireUvw = "0 0 1";
    private decimal _wireR = 10, _wireRep = 2, _wireVac = 10;
    private int _wireShape;
    public string WireUvw { get => _wireUvw; set => Set(ref _wireUvw, value ?? ""); }
    public decimal? WireRadius { get => _wireR; set { if (value != null) Set(ref _wireR, Math.Clamp(value.Value, 1, 200)); } }
    public decimal? WireRepeats { get => _wireRep; set { if (value != null) Set(ref _wireRep, Math.Clamp(value.Value, 1, 100)); } }
    public decimal? WireVacuum { get => _wireVac; set { if (value != null) Set(ref _wireVac, Math.Clamp(value.Value, 0, 200)); } }
    public int WireShape { get => _wireShape; set => Set(ref _wireShape, Math.Clamp(value, 0, 2)); }
    public void MakeNanowire()
    {
        var uvw = _wireUvw.Split((char[])[' ', ',', '[', ']'], StringSplitOptions.RemoveEmptyEntries).Select(x => int.TryParse(x, out var k) ? k : int.MinValue).ToArray();
        if (uvw.Length != 3 || uvw.Contains(int.MinValue)) { CellToolText = "The direction is three integers, e.g. 1 1 0"; return; }
        CellEdit(new { op = "nanowire", uvw, radius = (double)_wireR, repeats = (int)_wireRep, shape = WireShapeIds[_wireShape], vacuum = (double)_wireVac });
    }
}
