using System;
using System.Collections.ObjectModel;
using System.IO;
using System.Text.Json.Nodes;
using CapsStudio.Interop;
using System.Globalization;
using System.Linq;
using System.Threading.Tasks;

namespace CapsStudio.ViewModels;

/// <summary>Analyze › Diffusion (design/boards/Diffusion): the molecule-centre MSD on log–log axes with its fit window,
/// the local slope (1 where the motion is diffusive), the Einstein D and the Yeh–Hummer finite-size correction
/// D∞ = D_PBC + k_B T ξ / (6π η L), ξ = 2.837297, for a cubic box of edge L.</summary>
public sealed partial class MainViewModel
{
    public bool IsDiffusion => _module == 49;
    public ResultCell DfPbc { get; } = new("D_PBC (fit)", "Einstein, slope/6");
    public ResultCell DfInf { get; } = new("D∞ (corrected)", "Yeh & Hummer 2004");
    public ResultCell DfSlope { get; } = new("Slope in the fit window", "d ln MSD / d ln t");
    private decimal _dfEta = 0.89m, _dfT = 300m, _dfL;
    public decimal DfEta { get => _dfEta; set { if (Set(ref _dfEta, Math.Clamp(value, 0.001m, 1e6m))) DfCorrect(); } }
    public decimal DfT { get => _dfT; set { if (Set(ref _dfT, Math.Clamp(value, 1m, 5000m))) DfCorrect(); } }
    public decimal DfL { get => _dfL; set { if (Set(ref _dfL, Math.Clamp(value, 0m, 1000m))) DfCorrect(); } }
    public decimal DfFitFrom { get => Analyze.FitFromD * 100; set { Analyze.FitFromD = value / 100; Raise(); } }
    public decimal DfFitTo { get => Analyze.FitToD * 100; set { Analyze.FitToD = value / 100; Raise(); } }
    private string _dfStatus = "", _dfCorrection = "", _dfWarning = "";
    public string DfStatus { get => _dfStatus; private set => Set(ref _dfStatus, value); }
    public string DfCorrectionText { get => _dfCorrection; private set => Set(ref _dfCorrection, value); }
    public string DfWarning { get => _dfWarning; private set { Set(ref _dfWarning, value); Raise(nameof(DfHasWarning)); } }
    public bool DfHasWarning => _dfWarning.Length > 0;
    public (double X, double Y)[] DfMsd { get; private set; } = [];
    public (double X, double Y)[] DfSlopeCurve { get; private set; } = [];
    public (double From, double To)? DfWindow { get; private set; }
    private double _dfPbc = double.NaN;
    public event Action? DfChanged;

    public void OpenDiffusion()
    {
        SetModule(49);
        if (_doc != null && _dfL == 0 && _doc.Summary() is { CellValid: not 0 } s) DfL = (decimal)Math.Round(Math.Pow(s.Volume, 1.0 / 3) / 10, 3);   // Å → nm, the cube of the same volume
        DfCorrect();
        DfChanged?.Invoke();
    }

    /// <summary>Yeh & Hummer (2004): k_B T ξ / (6π η L), in m²/s.</summary>
    public static double YehHummer(double temperatureK, double etaMPaS, double boxNm) =>
        boxNm > 0 && etaMPaS > 0 ? 1.380649e-23 * temperatureK * 2.837297 / (6 * Math.PI * etaMPaS * 1e-3 * boxNm * 1e-9) : double.NaN;

    private void DfCorrect()
    {
        var inv = CultureInfo.InvariantCulture;
        var c = YehHummer((double)_dfT, (double)_dfEta, (double)_dfL);
        DfCorrectionText = double.IsFinite(c) ? c.ToString("0.00e0", inv) + " m²/s" : "—";
        if (double.IsFinite(_dfPbc) && double.IsFinite(c))
        {
            DfInf.Value = (_dfPbc + c).ToString("0.00e0", inv) + " m²/s";
            DfInf.Caption = $"+ {(c / _dfPbc * 100).ToString("0", inv)} % for L = {_dfL.ToString("0.###", inv)} nm, η = {_dfEta.ToString("0.###", inv)} mPa·s, {_dfT.ToString("0", inv)} K";
        }
    }

    // ---- replicas (design/boards/Diffusion "Replicas · spread"): other runs of the same system, each fitted as this one,
    // D given as mean ± standard deviation over all runs
    public ObservableCollection<DfReplica> DfReplicas { get; } = new();
    public ResultCell DfSpread { get; } = new("Replicas · spread", "add other runs of the same system");
    public bool HasDfReplicas => DfReplicas.Count > 0;
    public void AddDfReplicas(IEnumerable<string> paths)
    {
        foreach (var p in paths)
            if (File.Exists(p) && DfReplicas.All(r => r.Path != p)) DfReplicas.Add(new DfReplica(p, TopologyFor(p)));
        Raise(nameof(HasDfReplicas));
        DfSpread.Caption = $"{DfReplicas.Count + 1} runs · Run to fit them all";
    }
    public void RemoveDfReplica(DfReplica r) { DfReplicas.Remove(r); Raise(nameof(HasDfReplicas)); DfSpreadUpdate(); }

    /// <summary>D of one run (m²/s): its file opened on its own and fitted with this page's options.</summary>
    internal static double DiffusionOf(string path, string topology, CapsAnalyzeOpts o)
    {
        using var doc = CapsDocument.Open(path, topology.Length > 0 ? topology : null);
        var j = JsonNode.Parse(doc.Analyze("diffusion", o, null))!;
        var p = (j["properties"] as JsonArray)?.OfType<JsonObject>().FirstOrDefault(x => (string?)x["id"] == "diffusion");
        if (p?["extra"] is JsonObject ex)
            foreach (var kv in ex)
                if (kv.Key.StartsWith("D (m²/s)", StringComparison.Ordinal) && kv.Value is JsonValue v && v.TryGetValue<double>(out var dv)) return dv;
        return double.NaN;
    }

    private async Task RunDfReplicas()
    {
        var o = Analyze.Options();
        var inv = CultureInfo.InvariantCulture;
        foreach (var r in DfReplicas.ToList())
        {
            DfStatus = $"Replica {r.Name}…";
            try { r.D = await Task.Run(() => DiffusionOf(r.Path, r.Topology, o)); r.Fitted = true; }
            catch (Exception e) { r.D = double.NaN; r.Note = e.Message; }
        }
        DfSpreadUpdate();
    }

    private void DfSpreadUpdate()
    {
        var inv = CultureInfo.InvariantCulture;
        var (mean, sd, n) = Spread(new[] { _dfPbc }.Concat(DfReplicas.Select(r => r.D)));
        if (n < 2)
        {
            DfSpread.Value = "—";
            DfSpread.Caption = DfReplicas.Count == 0 ? "add other runs of the same system" : $"needs two fitted runs ({n} fitted: a run too short for its fit window has no D)";
            return;
        }
        DfSpread.Value = $"{mean.ToString("0.00e0", inv)} ± {sd.ToString("0.0e0", inv)} m²/s";
        DfSpread.Caption = $"mean ± s.d. of {n} runs ({(sd / mean * 100).ToString("0", inv)} %); standard error {(sd / Math.Sqrt(n)).ToString("0.0e0", inv)} m²/s";
    }
    /// <summary>Mean and sample standard deviation of the finite values, and how many there were.</summary>
    public static (double Mean, double Sd, int N) Spread(IEnumerable<double> values)
    {
        var ds = values.Where(double.IsFinite).ToList();
        if (ds.Count == 0) return (double.NaN, double.NaN, 0);
        var mean = ds.Average();
        return (mean, ds.Count > 1 ? Math.Sqrt(ds.Sum(x => (x - mean) * (x - mean)) / (ds.Count - 1)) : double.NaN, ds.Count);
    }

    public async Task RunDiffusion()
    {
        if (_doc == null || Analyze.Working) return;
        var inv = CultureInfo.InvariantCulture;
        DfStatus = "MSD over all time origins…";
        await RunChips("msd", "diffusion");
        var msd = Analyze.Results.FirstOrDefault(r => r.Id == "msd");
        var d = Analyze.Results.FirstOrDefault(r => r.Id == "diffusion");
        var centres = msd == null ? null : Analyze.Curves.FirstOrDefault(c => c.Property == msd.Name && c.Label == "molecule centres");
        var slope = msd == null ? null : Analyze.Curves.FirstOrDefault(c => c.Property == msd.Name && c.Label == "log-log slope");
        DfMsd = centres == null ? [] : centres.X.Zip(centres.Y).Where(p => p.First > 0 && p.Second > 0).Select(p => (Math.Log10(p.First), Math.Log10(p.Second))).ToArray();
        DfSlopeCurve = slope == null ? [] : slope.X.Zip(slope.Y).Where(p => p.First > 0 && double.IsFinite(p.Second)).Select(p => (Math.Log10(p.First), p.Second)).ToArray();
        double X(string k) => d == null ? double.NaN : d.Extra.Where(e => e.Key.StartsWith(k, StringComparison.Ordinal)).Select(e => e.Value).DefaultIfEmpty(double.NaN).First();
        var from = X("fit from"); var to = X("fit to");
        DfWindow = from > 0 && to > from ? (Math.Log10(from), Math.Log10(to)) : null;
        _dfPbc = X("D (m²/s)");
        var beta = X("log-log slope");
        DfPbc.Value = double.IsFinite(_dfPbc) ? _dfPbc.ToString("0.00e0", inv) + " m²/s" : "—";
        DfPbc.Caption = d != null && double.IsFinite(d.Error) ? $"± {(d.Error * 1e-9).ToString("0.0e0", inv)} m²/s (window halves)" : "Einstein, slope/6";
        DfSlope.Value = double.IsFinite(beta) ? beta.ToString("0.00", inv) : "—";
        DfWarning = !double.IsFinite(beta) ? (d?.Notes.FirstOrDefault() ?? "")
                  : beta < 0.9 ? $"The fit window is not diffusive yet (slope {beta.ToString("0.00", inv)} < 0.9): sub-diffusive, so D is an upper bound. Run longer or move the window later."
                  : beta > 1.1 ? $"The fit window includes ballistic motion (slope {beta.ToString("0.00", inv)} > 1.1): move it later."
                  : "";
        if (!double.IsFinite(_dfPbc)) DfInf.Value = "—";
        DfCorrect();
        if (DfReplicas.Count > 0) await RunDfReplicas();
        DfStatus = Analyze.Log;
        DfChanged?.Invoke();
    }
}

/// <summary>Another run of the same system (Diffusion replicas): its file, topology and fitted D.</summary>
public sealed class DfReplica : ObservableObject
{
    public DfReplica(string path, string topology) { Path = path; Topology = topology; }
    public string Path { get; }
    public string Topology { get; }
    public string Name => System.IO.Path.GetFileName(Path);
    private double _d = double.NaN;
    public double D { get => _d; set { if (Set(ref _d, value)) Raise(nameof(DText)); } }
    private string _note = "";
    public string Note { get => _note; set { if (Set(ref _note, value)) Raise(nameof(DText)); } }
    /// <summary>Its analysis ran (D may still be NaN: a run too short for the fit window).</summary>
    private bool _fitted;
    public bool Fitted { get => _fitted; set { if (Set(ref _fitted, value)) Raise(nameof(DText)); } }
    public string DText => double.IsFinite(_d) ? _d.ToString("0.00e0", CultureInfo.InvariantCulture) + " m²/s" : _note.Length > 0 ? "could not open" : Fitted ? "no D (too short)" : "not run";
}
