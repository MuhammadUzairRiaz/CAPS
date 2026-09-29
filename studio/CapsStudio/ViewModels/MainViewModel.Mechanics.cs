using System.Collections.ObjectModel;
using System.Globalization;

namespace CapsStudio.ViewModels;

/// <summary>A labelled result cell of a focused Analyze page.</summary>
public sealed class ResultCell : ObservableObject
{
    public ResultCell(string label, string caption) { _label = label; _caption = caption; }
    private string _label;
    public string Label { get => _label; set => Set(ref _label, value); }
    private string _value = "—", _caption;
    public string Value { get => _value; set { if (Set(ref _value, value)) Raise(nameof(HasValue)); } }
    public string Caption { get => _caption; set => Set(ref _caption, value); }
    public bool HasValue => _value != "—";
}

/// <summary>One stiffness entry (C11 … C66) of the matrix card.</summary>
public sealed record CijCell(string Label, string Value, bool Blank = false);

/// <summary>Analyze › Mechanics (design/boards/Mechanics): elastic constants by constant strain (or from stress
/// fluctuations) and uniaxial stress–strain at a stated rate — E, ν, K and yield, the curve with its fit window, the
/// stiffness matrix and the protocol. Runs through the Analyze engine with only these calculations switched on.</summary>
public sealed partial class MainViewModel
{
    public bool IsMechanics => _module == 38;
    public static int FocusOf(string id) => id switch
    {
        "cij_strain" or "cij_run" or "tensile" => 38,
        "sq" or "xray" or "electron" or "neutron" => 39,
        "ffv" or "psd" => 40,
        _ => 1,
    };

    /// <summary>A calculation picked in the side list: its focused page, or the Analyze properties with it switched on.</summary>
    public void FocusChip(CalcChip c)
    {
        var m = FocusOf(c.Id);
        if (m == 1) { c.IsOn = true; SetModule(1); return; }
        if (m == 38) OpenMechanics();
        else if (m == 39) OpenScattering();
        else if (m == 40) OpenFreeVolume();
    }

    private void SyncFocusChips()
    {
        foreach (var chip in Analyze.Groups.SelectMany(g => g.Chips))
            chip.Active = _module switch
            {
                38 => chip.Id == "tensile" && _mechTensile || chip.Id == (_mechMethod == 0 ? "cij_strain" : "cij_run"),
                39 => chip.Id == ScatterBeamId && _scatterXray || chip.Id == "neutron" && _scatterNeutron,
                40 => chip.Id is "ffv" or "psd",
                _ => false,
            };
    }

    public ObservableCollection<ResultCell> MechCells { get; } =
    [
        new("Young's modulus E", "fit 0–2 % strain"), new("Poisson ratio ν", "from lateral strain"),
        new("Bulk modulus K", "(C₁₁ + 2C₁₂)/3, Voigt"), new("Yield stress", "0.2 % offset"),
    ];
    public ObservableCollection<CijCell> CijCells { get; } = new();
    public static readonly string[] MechMethods = ["Constant strain", "Stress fluctuations"];
    public static readonly string[] MechRates = ["1e8 s⁻¹", "1e9 s⁻¹", "1e10 s⁻¹"];
    public static readonly string[] MechLaterals = ["NPT · 1 atm", "Fixed (uniaxial strain)"];
    private int _mechMethod, _mechRate = 1;
    private bool _mechTensile = true;
    public int MechMethod { get => _mechMethod; set { if (Set(ref _mechMethod, value)) { Raise(nameof(MechCite)); Raise(nameof(MechIsStrain)); SyncFocusChips(); } } }
    public bool MechIsStrain => _mechMethod == 0;
    public bool MechTensile { get => _mechTensile; set { if (Set(ref _mechTensile, value)) SyncFocusChips(); } }
    public int MechRate { get => _mechRate; set { if (Set(ref _mechRate, value)) Analyze.TensRateD = value switch { 0 => 1e-4m, 2 => 1e-2m, _ => 1e-3m }; } }
    public int MechLateral { get => Analyze.TensFixed ? 1 : 0; set { Analyze.TensFixed = value == 1; Raise(); } }
    public decimal MechStrainPercent { get => Analyze.CijStrainD * 100; set { Analyze.CijStrainD = value / 100; Raise(); } }
    public string MechCite => _mechMethod == 0
        ? "Theodorou & Suter, Macromolecules 19, 139 (1986) — static elastic constants of glassy polymers"
        : "Lutsko, J. Appl. Phys. 65, 2991 (1989); Clavier et al., Mol. Simul. 43, 1413 (2017) — elastic constants from stress fluctuations";
    public string MechNote => _mechMethod == 0
        ? $"±{Analyze.CijStrainD * 100:0.###} % strain per component, stress from the virial, averaged over {Analyze.CijConfigsD:0} relaxed configuration{(Analyze.CijConfigsD == 1 ? "" : "s")}."
        : $"Stress fluctuations over {Analyze.FluctPsD:0} ps of NVT at {Analyze.FluctTD:0} K; bonded polymers need long runs (errors stay GPa-level over 100 ps).";
    public string MechCurveChip => $"tension along {MainViewModel.AxisName(Analyze.TensAxis)} · {MechRates[_mechRate]}";
    private static string AxisName(int a) => a switch { 1 => "y", 2 => "z", _ => "x" };
    public string MechFooter => _doc == null ? "" : $"{Title} · {_doc.Summary().Atoms:N0} atoms";
    public event Action? MechCurveChanged;
    public (double X, double Y)[] MechCurve { get; private set; } = [];
    public (double X, double Y)[] MechFit { get; private set; } = [];
    public double MechFitTo => 2.0;   // % strain (core fit_strain 0.02)
    public bool MechHasCurve => MechCurve.Length > 1;

    public void OpenMechanics()
    {
        SetModule(38);
        FillMechanics();
    }

    /// <summary>Runs the chosen elastic-constant method (and the tensile test) with only these calculations on.</summary>
    public async Task RunMechanics()
    {
        if (_doc == null || Analyze.Working) return;
        var chips = Analyze.Groups.SelectMany(g => g.Chips).ToList();
        var was = chips.ToDictionary(c => c, c => c.IsOn);
        foreach (var c in chips) c.IsOn = false;
        (_mechMethod == 0 ? Analyze.StrainChip : Analyze.FluctChip).IsOn = true;
        Analyze.TensileChip.IsOn = _mechTensile;
        try { await Analyze.Run(); }
        finally { foreach (var (c, on) in was) c.IsOn = on; }
        FillMechanics();
    }

    private void FillMechanics()
    {
        var inv = CultureInfo.InvariantCulture;
        ResultCard? R(params string[] ids) => ids.Select(id => Analyze.Results.FirstOrDefault(r => r.Id == id)).FirstOrDefault(r => r != null && r.HasValue);
        double X(ResultCard? c, string key) => c?.Extra.FirstOrDefault(e => e.Key == key) is { Key: not null } kv ? kv.Value : double.NaN;
        string V(double v, string f, string unit = "") => double.IsNaN(v) ? "—" : v.ToString(f, inv) + unit;
        var tens = R("tensile_modulus");
        var suffix = _mechMethod == 0 ? "" : "_fluct";
        var cij = R("cij" + suffix);
        var E = tens?.Value ?? R("youngs" + suffix)?.Value ?? double.NaN;
        MechCells[0].Value = V(E, "0.00", " GPa");
        MechCells[0].Caption = tens != null ? "fit 0–2 % strain" : "Hill average of Cij";
        var nu = tens != null ? X(tens, "Poisson's ratio") : double.NaN;
        if (double.IsNaN(nu)) nu = R("poisson" + suffix)?.Value ?? double.NaN;
        MechCells[1].Value = V(nu, "0.000");
        MechCells[1].Caption = tens != null && !double.IsNaN(X(tens, "Poisson's ratio")) ? "from lateral strain" : "Hill average of Cij";
        MechCells[2].Value = V(R("bulk" + suffix)?.Value ?? double.NaN, "0.00", " GPa");
        MechCells[2].Caption = "Hill average of Cij";
        MechCells[3].Value = V(R("yield")?.Value ?? double.NaN, "0", " MPa");
        CijCells.Clear();
        string C(int i, int j) => V(X(cij, $"C{i}{j} (GPa)"), "0.00");
        foreach (var (i, j) in new[] { (1, 1), (1, 2), (1, 3), (0, 0), (2, 2), (2, 3), (0, 0), (0, 0), (3, 3), (4, 4), (5, 5), (6, 6) })
            CijCells.Add(i == 0 ? new CijCell("", "", true) : new CijCell($"C{Sub(i)}{Sub(j)}", cij == null ? "—" : C(i, j)));
        // the curve: strain % against stress GPa, with the linear fit over the window
        var series = Analyze.Curves.FirstOrDefault(c => c.Property == tens?.Name && c.Label.Contains("smoothed"));
        MechCurve = series == null ? [] : series.X.Zip(series.Y, (x, y) => (x * 100, y / 1000)).ToArray();
        MechFit = double.IsNaN(E) || tens == null ? [] : [(0, 0), (MechFitTo * 1.5, E * MechFitTo * 1.5 / 100)];
        foreach (var n in new[] { nameof(MechHasCurve), nameof(MechCurveChip), nameof(MechNote), nameof(MechFooter) }) Raise(n);
        MechCurveChanged?.Invoke();
    }

    private static string Sub(int k) => ((char)('₀' + k)).ToString();
}
