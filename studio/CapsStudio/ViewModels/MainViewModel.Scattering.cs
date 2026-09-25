using System.Collections.ObjectModel;
using System.Globalization;
using CapsStudio.Interop;

namespace CapsStudio.ViewModels;

/// <summary>An X-ray source: name and wavelength (Å).</summary>
public sealed record XraySource(string Name, double Lambda)
{
    public override string ToString() => $"{Name} · λ {Lambda.ToString("0.0000", CultureInfo.InvariantCulture)} Å";
}

/// <summary>Analyze › Scattering (design/boards/Scattering): X-ray and neutron patterns of the cell — the direct
/// reciprocal-lattice sum at low q and the g(r) transform above, Cromer–Mann form factors, coherent scattering lengths
/// with a deuteration pattern for contrast — peaks in q, d and 2θ for a chosen source, and a measured curve overlaid.</summary>
public sealed partial class MainViewModel
{
    public bool IsScattering => _module == 39;
    public static readonly XraySource[] XraySources =
    [
        new("Cu Kα", 1.5406), new("Mo Kα", 0.7107), new("Co Kα", 1.7890), new("Cr Kα", 2.2910), new("Ag Kα", 0.5594),
    ];
    public static readonly string[] IsotopePatterns = ["Natural (all ¹H)", "Fully deuterated", "d-backbone, h-ring", "h-backbone, d-ring", "Exchangeable H → D (O–H, N–H)"];
    public ObservableCollection<ProvRow> ScatterLengths { get; } = new();
    public ObservableCollection<ProvRow> ScatterPeaks { get; } = new();

    private int _xraySource, _isotope;
    private bool _scatterXray = true, _scatterNeutron = true;
    public int XraySourceIndex { get => _xraySource; set { if (Set(ref _xraySource, value)) FillScattering(); } }
    public int IsotopePattern { get => _isotope; set { if (Set(ref _isotope, value)) { Analyze.Deuterate = value; Raise(nameof(NeutronChip)); } } }
    public bool ScatterXray { get => _scatterXray; set { if (Set(ref _scatterXray, value)) { SyncFocusChips(); ScatterChanged?.Invoke(); } } }
    public bool ScatterNeutron { get => _scatterNeutron; set { if (Set(ref _scatterNeutron, value)) { SyncFocusChips(); ScatterChanged?.Invoke(); } } }
    private decimal _scatterPlotTo = 6;
    /// <summary>The plotted q range ends here (the computed range may reach further).</summary>
    public decimal ScatterPlotTo { get => _scatterPlotTo; set { if (Set(ref _scatterPlotTo, Math.Clamp(value, 1, 40))) FillScattering(); } }
    public string XrayChip => "X-ray (" + XraySources[_xraySource].Name + ")";
    public string NeutronChip => "neutron · " + IsotopePatterns[_isotope].ToLowerInvariant();
    public string ScatterTitle => "I(q) · " + (_doc != null ? Title : "no structure");
    private string _scatterText = "", _expName = "";
    public string ScatterText { get => _scatterText; private set => Set(ref _scatterText, value); }
    public string ExperimentName { get => _expName; private set { if (Set(ref _expName, value)) Raise(nameof(HasExperiment)); } }
    public bool HasExperiment => _expName.Length > 0;
    public event Action? ScatterChanged;
    public (double X, double Y)[] ScatterXrayCurve { get; private set; } = [];
    public (double X, double Y)[] ScatterNeutronCurve { get; private set; } = [];
    public (double X, double Y)[] ScatterExperiment { get; private set; } = [];
    private (double X, double Y)[] _expRaw = [];
    public bool ScatterHasCurve => ScatterXrayCurve.Length > 1 || ScatterNeutronCurve.Length > 1;
    public string ScatterFooter => _doc == null ? "" : $"{Analyze.FramesText} · {_doc.Summary().Atoms:N0} atoms";

    public void OpenScattering()
    {
        SetModule(39);
        FillLengths();
        FillScattering();
    }

    public async Task RunScattering()
    {
        if (_doc == null || Analyze.Working) return;
        var chips = Analyze.Groups.SelectMany(g => g.Chips).ToList();
        var was = chips.ToDictionary(c => c, c => c.IsOn);
        foreach (var c in chips) c.IsOn = c.Id == "xray" && _scatterXray || c.Id == "neutron" && _scatterNeutron;
        try { await Analyze.Run(); }
        finally { foreach (var (c, on) in was) c.IsOn = on; }
        FillScattering();
    }

    /// <summary>The coherent scattering lengths of the elements in the structure (and ²H when deuterating).</summary>
    private void FillLengths()
    {
        ScatterLengths.Clear();
        if (_doc == null) return;
        var inv = CultureInfo.InvariantCulture;
        var n = (int)Math.Min(_doc.Summary().Atoms, 200000);
        var zs = new SortedDictionary<int, string>();
        for (int i = 0; i < n; i++) { var a = _doc.Atom(i); zs.TryAdd(a.Element, a.ElementSymbol); }
        foreach (var (z, symbol) in zs)
        {
            var b = Native.NeutronB(z);
            var sym = z == 1 ? "¹H" : symbol;
            ScatterLengths.Add(new ProvRow(sym, double.IsNaN(b) ? "—" : b.ToString("0.000", inv)));
            if (z == 1) ScatterLengths.Add(new ProvRow("²H (D)", Native.NeutronB(1001).ToString("0.000", inv)));
        }
    }

    private void FillScattering()
    {
        var inv = CultureInfo.InvariantCulture;
        (double X, double Y)[] Curve(string id)
        {
            var card = Analyze.Results.FirstOrDefault(r => r.Id == id);
            var s = card == null ? null : Analyze.Curves.FirstOrDefault(c => c.Property == card.Name);
            if (s == null) return [];
            var to = (double)_scatterPlotTo;
            var pts = s.X.Zip(s.Y).Where(p => double.IsFinite(p.First) && double.IsFinite(p.Second) && p.First >= 0.2 && p.First <= to).ToArray();
            // normalised to the highest point above 0.4 Å⁻¹: a low-density cell's small-angle upturn would flatten the halo
            var max = pts.Where(p => p.First >= 0.4).Select(p => p.Second).DefaultIfEmpty(1).Max();
            return max > 0 ? pts.Select(p => (p.First, Math.Min(1.6, p.Second / max))).ToArray() : pts;
        }
        ScatterXrayCurve = Curve("xray");
        ScatterNeutronCurve = Curve("neutron");
        ScatterPeaks.Clear();
        var lambda = XraySources[_xraySource].Lambda;
        var basis = ScatterXrayCurve.Length > 1 ? ScatterXrayCurve : ScatterNeutronCurve;
        foreach (var (q, _) in Peaks(basis).Take(4))
        {
            var s = q * lambda / (4 * Math.PI);
            var tth = s < 1 ? 2 * Math.Asin(s) * 180 / Math.PI : double.NaN;
            ScatterPeaks.Add(new ProvRow(q.ToString("0.00", inv) + " Å⁻¹", (2 * Math.PI / q).ToString("0.00", inv) + " Å",
                                         double.IsNaN(tth) ? "—" : tth.ToString("0.0", inv) + "°"));
        }
        ScatterText = ScatterPeaks.Count == 0 && basis.Length > 1
            ? "No maxima stand above the small-angle rise in the plotted range: a low-density or porous cell scatters most strongly near the box's first shell (2π/L). Check the density before comparing with a measured glass."
            : ScatterPeaks.Count == 0
            ? "Run to compute the patterns. The low-q part is the exact sum over the cell's reciprocal lattice (every k-vector the box allows); above it, the partial g(r) are transformed."
            : $"Maxima at q ≈ {string.Join(", ", ScatterPeaks.Select(p => p.Key.Split(' ')[0]))} Å⁻¹ " +
              $"(d = 2π/q ≈ {string.Join(", ", ScatterPeaks.Select(p => p.Value.Split(' ')[0]))} Å; 2θ for {XraySources[_xraySource].Name}: {string.Join(", ", ScatterPeaks.Select(p => p.Other))}). " +
              "Below 2π/L the cell holds no k-vectors: the curve starts at the box's first shell.";
        ScaleExperiment();
        foreach (var n in new[] { nameof(ScatterHasCurve), nameof(XrayChip), nameof(NeutronChip), nameof(ScatterTitle), nameof(ScatterFooter) }) Raise(n);
        ScatterChanged?.Invoke();
    }

    /// <summary>Local maxima of a normalised curve (smoothed over 5 points), highest first then by q.</summary>
    private static IEnumerable<(double Q, double I)> Peaks((double X, double Y)[] c)
    {
        if (c.Length < 7) return [];
        var sm = new double[c.Length];
        for (int i = 0; i < c.Length; i++)
        {
            int a = Math.Max(0, i - 2), b = Math.Min(c.Length - 1, i + 2);
            double s = 0;
            for (int k = a; k <= b; k++) s += c[k].Y;
            sm[i] = s / (b - a + 1);
        }
        var list = new List<(double, double)>();
        for (int i = 3; i < c.Length - 3; i++)
            if (sm[i] > sm[i - 1] && sm[i] >= sm[i + 1] && sm[i] > sm[i - 3] + 0.01 && sm[i] > sm[i + 3] + 0.01 && sm[i] > 0.3 && c[i].X >= 0.4) list.Add((c[i].X, sm[i]));
        return list.OrderByDescending(p => p.Item2).Take(4).OrderBy(p => p.Item1);
    }

    /// <summary>Reads a measured pattern: two columns (q in Å⁻¹, or 2θ in degrees when the header says so or x runs past 20).</summary>
    public string? LoadExperiment(string path)
    {
        var inv = CultureInfo.InvariantCulture;
        var pts = new List<(double, double)>();
        var twoTheta = false;
        foreach (var raw in File.ReadLines(path))
        {
            var line = raw.Trim();
            if (line.Length == 0) continue;
            if (line[0] is '#' or '"' || char.IsLetter(line[0]))
            {
                var l = line.ToLowerInvariant();
                if (l.Contains("2theta") || l.Contains("2θ") || l.Contains("tth") || l.Contains("two_theta")) twoTheta = true;
                continue;
            }
            var f = line.Split([',', ';', '\t', ' '], StringSplitOptions.RemoveEmptyEntries);
            if (f.Length >= 2 && double.TryParse(f[0], NumberStyles.Float, inv, out var x) && double.TryParse(f[1], NumberStyles.Float, inv, out var y)) pts.Add((x, y));
        }
        if (pts.Count < 3) return "No two-column data found";
        if (!twoTheta && pts.Max(p => p.Item1) > 20) twoTheta = true;
        var lambda = XraySources[_xraySource].Lambda;
        _expRaw = pts.Select(p => twoTheta ? (4 * Math.PI * Math.Sin(p.Item1 * Math.PI / 360) / lambda, p.Item2) : p).OrderBy(p => p.Item1).ToArray();
        ExperimentName = Path.GetFileName(path) + (twoTheta ? $" · 2θ with {XraySources[_xraySource].Name}" : " · q");
        ScaleExperiment();
        ScatterChanged?.Invoke();
        return null;
    }

    public void ClearExperiment() { _expRaw = []; ScatterExperiment = []; ExperimentName = ""; ScatterChanged?.Invoke(); }

    private void ScaleExperiment()
    {
        if (_expRaw.Length < 3) { ScatterExperiment = []; return; }
        var max = _expRaw.Where(p => p.X >= 0.2).Select(p => p.Y).DefaultIfEmpty(1).Max();
        ScatterExperiment = max > 0 ? _expRaw.Select(p => (p.X, p.Y / max)).ToArray() : _expRaw;
    }
}
