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
            Raise(nameof(PackIntoCurrent));
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
    // stage 3: pack loosely, then compress the cell to a density (push-off and minimisation between affine steps)
    /// <summary>Packing around the current structure: it stays fixed, so the cell cannot be compressed.</summary>
    public bool PackIntoCurrent => _packStart == 1;
    private bool _packCompress;
    private decimal _packCompressTo = 0.9m;
    public bool PackCompress { get => _packCompress; set => Set(ref _packCompress, value); }
    public decimal PackCompressTo { get => _packCompressTo; set => Set(ref _packCompressTo, Math.Clamp(value, 0.05m, 5m)); }
    private string PackTextToRun()
    {
        var text = PackTextToRunCore();
        if (_packCompress && _packStart != 1 && !text.Split('\n').Any(l => l.TrimStart().StartsWith("compress ", StringComparison.OrdinalIgnoreCase)))
            text = string.Format(CultureInfo.InvariantCulture, "compress {0:0.###}\n", _packCompressTo) + text;
        return text;
    }
    private string PackTextToRunCore()
    {
        if (_packStart != 1 || _doc == null) return _packText;
        if (!PackCanUseCurrent) throw new InvalidOperationException("the current structure has no orthorhombic cell to pack into");
        var dir = Path.Combine(Path.GetTempPath(), "caps-pack");
        Directory.CreateDirectory(dir);
        var host = Path.Combine(dir, "host_" + Environment.ProcessId + ".data");
        _doc.Save(host);
        var s = _doc.Summary();
        var inv = CultureInfo.InvariantCulture;
        var caps = IsCapsPack(_packText);
        var cell = caps ? string.Format(inv, "cell      {0:0.####} {1:0.####} {2:0.####}", s.CellA, s.CellB, s.CellC)
                        : string.Format(inv, "pbc 0. 0. 0. {0:0.####} {1:0.####} {2:0.####}", s.CellA, s.CellB, s.CellC);
        var lines = _packText.Split('\n').Where(l => !l.TrimStart().StartsWith(caps ? "cell " : "pbc ", StringComparison.OrdinalIgnoreCase)).ToList();
        var first = lines.FindIndex(l => l.TrimStart().StartsWith(caps ? "molecule " : "structure ", StringComparison.OrdinalIgnoreCase));
        if (first < 0) first = lines.Count;
        var title = Title.Replace(" (unsaved)", "");
        lines.Insert(first, caps ? $"{cell}\n\nmolecule  {host}   # {title}, kept where it is\n  count   1\n  fixed   at 0 0 0\nend\n"
                                 : $"{cell}\n\nstructure {host}   # {title}, kept where it is\n  number 1\n  fixed 0. 0. 0. 0. 0. 0.\nend structure\n");
        return string.Join('\n', lines);
    }

    // ---- fill to a density: the number of the last molecule block from the target density of the whole cell
    private decimal _fillDensity = 0.9m;
    private string _fillText = "";
    public decimal? FillDensity { get => _fillDensity; set { if (value != null) Set(ref _fillDensity, Math.Clamp(value.Value, 0.01m, 5m)); } }
    public string FillText { get => _fillText; private set => Set(ref _fillText, value); }
    /// <summary>Sets the last structure block's number so the cell (the current structure plus the molecules) reaches the
    /// density: N = (ρ V N_A − m_host) / m_molecule, rounded down (ρ in g/cm³, V the cell's volume).</summary>
    public void FillToDensity()
    {
        var inv = CultureInfo.InvariantCulture;
        if (_doc == null || _packStart != 1) { FillText = "Pack around the current structure first (Pack into › Around the current structure)"; return; }
        var lines = _packText.Split('\n').ToList();
        var caps = IsCapsPack(_packText);
        string open = caps ? "molecule " : "structure ", countKey = caps ? "count " : "number ", endKey = caps ? "end" : "end structure";
        var k = lines.FindLastIndex(l => l.TrimStart().StartsWith(open, StringComparison.OrdinalIgnoreCase));
        if (k < 0) { FillText = "Add the molecule to fill with first"; return; }
        var path = lines[k].Trim()[open.Length..].Split('#')[0].Trim();
        if (!Path.IsPathRooted(path)) path = Path.Combine(_packBaseDir, path);
        double mMol;
        try { using var mol = Interop.CapsDocument.Open(path); mMol = mol.Summary().TotalMass; }
        catch (Exception e) { FillText = "Cannot read the molecule: " + e.Message; return; }
        var s = _doc.Summary();
        const double avogadroPerA3 = 0.602214076;   // g/cm³ × Å³ → g/mol
        var target = (double)_fillDensity * s.Volume * avogadroPerA3;
        var room = target - s.TotalMass;
        if (mMol <= 0) { FillText = "The molecule has no mass"; return; }
        if (room < mMol) { FillText = string.Format(inv, "The cell is already at {0:0.000} g/cm³: nothing to add for {1:0.000}", s.Density, _fillDensity); return; }
        var n = (int)Math.Floor(room / mMol);
        var j = k + 1;
        while (j < lines.Count && !lines[j].TrimStart().StartsWith(endKey, StringComparison.OrdinalIgnoreCase) && !lines[j].TrimStart().StartsWith(countKey, StringComparison.OrdinalIgnoreCase)) ++j;
        if (j < lines.Count && lines[j].TrimStart().StartsWith(countKey, StringComparison.OrdinalIgnoreCase)) lines[j] = $"  {countKey}{n}";
        else lines.Insert(k + 1, $"  {countKey}{n}");
        PackText = string.Join('\n', lines);
        var reached = (s.TotalMass + n * mMol) / (s.Volume * avogadroPerA3);
        FillText = string.Format(inv, "{0} × {1} ({2:0.0} g/mol): {3:0.000} → {4:0.000} g/cm³ (target {5:0.000}; whole molecules)", n, Path.GetFileNameWithoutExtension(path), mMol, s.Density, reached, _fillDensity);
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
