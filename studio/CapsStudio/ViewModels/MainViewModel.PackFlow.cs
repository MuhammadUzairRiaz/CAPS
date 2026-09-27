using System;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Threading.Tasks;

namespace CapsStudio.ViewModels;

/// <summary>Pack as its own workflow (as Materials Studio keeps Amorphous Cell's Packing apart from Construction):
/// molecules packed into an empty box, or around the current structure (a grown cell, a surface), which stays fixed;
/// its own force field, assigned when packing finishes; then the cell goes straight to LAMMPS or GROMACS, with
/// minimisation and dynamics in CAPS as optional steps.</summary>
public partial class MainViewModel
{
    public static readonly string[] PackStartModes = ["An empty box", "Around the current structure (kept fixed; molecules fill its cell)"];
    private int _packStart;
    public int PackStart
    {
        get => _packStart;
        set
        {
            var v = value == 1 && !PackCanUseCurrent ? 0 : Math.Clamp(value, 0, 1);
            if (!Set(ref _packStart, v)) return;
            if (v == 1 && _doc != null)
            {
                // the box is the structure's cell: new structure blocks go inside it
                var s = _doc.Summary();
                PackXD = (decimal)s.CellA; PackYD = (decimal)s.CellB; PackZD = (decimal)s.CellC;
            }
            Raise(nameof(PackStartNote));
        }
    }
    /// <summary>Packing around a structure needs its cell (orthorhombic: packmol's boxes are axis-aligned).</summary>
    public bool PackCanUseCurrent
    {
        get
        {
            if (_doc == null) return false;
            var s = _doc.Summary();
            return s.CellValid != 0 && Math.Abs(s.Volume - s.CellA * s.CellB * s.CellC) < 1e-6 * s.Volume;   // V = abc only with right angles
        }
    }
    public string PackStartNote => _packStart == 1 && _doc != null
        ? string.Format(CultureInfo.InvariantCulture, "{0} stays where it is; the molecules below fill the free space of its {1:0.#} × {2:0.#} × {3:0.#} Å cell (periodic)",
            Title.Replace(" (unsaved)", ""), _packX, _packY, _packZ)
        : "The molecules below are packed into the box of the input";

    // ---- the Pack's own force field (Grow keeps its own too)
    private int _packFf = -1, _packCharges;
    private bool _packAssign = true;
    public int PackFfIndex { get => _packFf < 0 ? Field.FfIndex : _packFf; set { _packFf = value; Raise(); } }
    public int PackChargeMode { get => _packCharges; set => Set(ref _packCharges, value); }
    public bool PackAssignField { get => _packAssign; set => Set(ref _packAssign, value); }

    // ---- after packing: export straight away, or minimise / run dynamics first
    private bool _packDone;
    public bool PackDone { get => _packDone; private set => Set(ref _packDone, value); }

    /// <summary>The packmol text actually run: with the current structure as a fixed block and its cell as the periodic box.</summary>
    private string PackTextToRun()
    {
        if (_packStart != 1 || _doc == null) return _packText;
        if (!PackCanUseCurrent) throw new InvalidOperationException("the current structure has no orthorhombic cell to pack into");
        var dir = Path.Combine(Path.GetTempPath(), "caps-pack");
        Directory.CreateDirectory(dir);
        var host = Path.Combine(dir, "host_" + Environment.ProcessId + ".data");
        _doc.Save(host);
        var s = _doc.Summary();
        var inv = CultureInfo.InvariantCulture;
        var pbc = string.Format(inv, "pbc 0. 0. 0. {0:0.####} {1:0.####} {2:0.####}", s.CellA, s.CellB, s.CellC);
        var lines = _packText.Split('\n').Where(l => !l.TrimStart().StartsWith("pbc ", StringComparison.OrdinalIgnoreCase)).ToList();
        var first = lines.FindIndex(l => l.TrimStart().StartsWith("structure ", StringComparison.OrdinalIgnoreCase));
        if (first < 0) first = lines.Count;
        lines.Insert(first, $"{pbc}\n\nstructure {host}   # {Title.Replace(" (unsaved)", "")}, kept where it is\n  number 1\n  fixed 0. 0. 0. 0. 0. 0.\nend structure\n");
        return string.Join('\n', lines);
    }

    /// <summary>After Grow or Pack made a cell: that builder's own force field, typed and checked.</summary>
    private async Task AssignForBuilder(bool pack)
    {
        if (_doc == null) return;
        if (pack ? !_packAssign : !_growAssignFf) return;
        var idx = pack ? PackFfIndex : GrowFfIndex;
        if (idx < 0 || idx >= Field.Library.Count) return;
        Field.FfIndex = idx;
        Field.ChargeMode = pack ? _packCharges : _growCharges;
        await Field.Assign();
        _pipeAutoFf = true;
        RaiseGrowField();
        RefreshSteps();
    }

    // ---- Grow's own force field
    private int _growFf = -1, _growCharges;
    public int GrowFfIndex { get => _growFf < 0 ? Field.FfIndex : _growFf; set { _growFf = value; Raise(); } }
    public int GrowChargeMode { get => _growCharges; set => Set(ref _growCharges, value); }

    /// <summary>Grow › Assign now: Grow's own force field on the structure in the viewer.</summary>
    public async Task AssignGrowFieldNow()
    {
        if (_doc == null) return;
        var idx = GrowFfIndex;
        if (idx >= 0 && idx < Field.Library.Count) Field.FfIndex = idx;
        Field.ChargeMode = _growCharges;
        await Field.Assign();
        RaiseGrowField();
        RefreshSteps();
    }

    /// <summary>Pack › Export: the packed cell to LAMMPS or GROMACS (the Export center, with its checks).</summary>
    public void PackExport() => OpenExportCenter();
}
