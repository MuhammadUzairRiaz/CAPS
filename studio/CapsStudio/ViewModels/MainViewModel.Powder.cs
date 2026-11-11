using System.Collections.ObjectModel;
using System.Globalization;
using System.Text.Json.Nodes;

namespace CapsStudio.ViewModels;

/// <summary>One reflection of the powder pattern.</summary>
public sealed record PowderRow(string Hkl, string D, string TwoTheta, string M, string I);

/// <summary>Analyze › Scattering › powder pattern (core powder.hpp): the cell's reflections with their hkl, d, 2θ and
/// intensity, and I(2θ) as a diffractometer records it — Kα doublet, Lorentz–polarisation, peak shapes with crystallite
/// size and strain, preferred orientation — beside a measured pattern.</summary>
public sealed partial class MainViewModel
{
    public static readonly string[] PowderAnodes = ["Cu Kα", "Mo Kα", "Co Kα", "Fe Kα", "Cr Kα", "Ag Kα"];
    private static readonly string[] PowderAnodeIds = ["Cu", "Mo", "Co", "Fe", "Cr", "Ag"];
    public static readonly string[] PowderProfiles = ["Pseudo-Voigt (Thompson–Cox–Hastings)", "Gaussian", "Lorentzian"];
    /// <summary>The full wording of each profile, shown as the dropdown's tooltip (the labels above stay short enough for the card).</summary>
    private static readonly string[] PowderProfileTips =
    [
        "Pseudo-Voigt whose Gaussian and Lorentzian widths come from the crystallite size, the strain and the instrument (Thompson, Cox and Hastings)",
        "Gaussian peaks from the crystallite size, the strain and the instrument",
        "Lorentzian peaks from the crystallite size, the strain and the instrument",
    ];
    private static readonly string[] PowderProfileIds = ["tch", "gaussian", "lorentzian"];
    private bool _powderView, _powderRunning, _powKa2 = true, _powAllFrames;
    private int _powAnode, _powProfile;
    private double _powFrom = 5, _powTo = 60, _powSizeNm, _powStrainPct, _powB, _powR = 1, _powW = 0.0025;
    private string _powPreferred = "", _powText = "", _powError = "";
    /// <summary>The plot shows the powder pattern in 2θ (else I(q)).</summary>
    public bool PowderView { get => _powderView; set { if (Set(ref _powderView, value)) { Raise(nameof(QView)); Raise(nameof(ScatterTitle)); ScatterChanged?.Invoke(); } } }
    public bool QView { get => !_powderView; set => PowderView = !value; }
    public bool PowderRunning { get => _powderRunning; private set { if (Set(ref _powderRunning, value)) RaiseBusy(); } }
    public int PowAnode { get => _powAnode; set => Set(ref _powAnode, Math.Clamp(value, 0, PowderAnodes.Length - 1)); }
    public int PowProfile { get => _powProfile; set { if (Set(ref _powProfile, Math.Clamp(value, 0, PowderProfiles.Length - 1))) Raise(nameof(PowProfileTip)); } }
    public string PowProfileTip => PowderProfileTips[_powProfile];
    public bool PowKa2 { get => _powKa2; set => Set(ref _powKa2, value); }
    public bool PowAllFrames { get => _powAllFrames; set => Set(ref _powAllFrames, value); }
    public decimal PowFromD { get => (decimal)_powFrom; set => Set(ref _powFrom, Math.Clamp((double)value, 0.5, 170)); }
    public decimal PowToD { get => (decimal)_powTo; set => Set(ref _powTo, Math.Clamp((double)value, 1, 175)); }
    public decimal PowSizeD { get => (decimal)_powSizeNm; set => Set(ref _powSizeNm, Math.Clamp((double)value, 0, 10000)); }
    public decimal PowStrainD { get => (decimal)_powStrainPct; set => Set(ref _powStrainPct, Math.Clamp((double)value, 0, 10)); }
    public decimal PowBD { get => (decimal)_powB; set => Set(ref _powB, Math.Clamp((double)value, 0, 50)); }
    public decimal PowRD { get => (decimal)_powR; set => Set(ref _powR, Math.Clamp((double)value, 0.05, 20)); }
    public decimal PowWD { get => (decimal)_powW; set => Set(ref _powW, Math.Clamp((double)value, 0, 1)); }
    public string PowPreferred { get => _powPreferred; set => Set(ref _powPreferred, value ?? ""); }
    public string PowText { get => _powText; private set => Set(ref _powText, value); }
    public string PowError { get => _powError; private set { if (Set(ref _powError, value)) Raise(nameof(PowHasError)); } }
    public bool PowHasError => _powError.Length > 0;
    public ObservableCollection<PowderRow> PowderRows { get; } = new();
    public (double X, double Y)[] PowderCurve { get; private set; } = [];
    private double _powLambda = 1.54056;
    /// <summary>The measured curve (held in q) redrawn in 2θ at the pattern's Kα1, scaled to the pattern's 100.</summary>
    public (double X, double Y)[] PowderExperiment
    {
        get
        {
            if (_expRaw.Length < 3 || PowderCurve.Length < 2) return [];
            var pts = _expRaw.Select(p => (S: p.X * _powLambda / (4 * Math.PI), p.Y)).Where(p => p.S > 0 && p.S < 1)
                             .Select(p => (X: 2 * Math.Asin(p.S) * 180 / Math.PI, p.Y)).Where(p => p.X >= _powFrom && p.X <= _powTo).ToArray();
            var max = pts.Select(p => p.Y).DefaultIfEmpty(0).Max();
            return max > 0 ? pts.Select(p => (p.X, 100 * p.Y / max)).ToArray() : [];
        }
    }

    public async Task RunPowder()
    {
        if (_doc == null || !Idle) return;
        var inv = CultureInfo.InvariantCulture;
        var jo = new JsonObject
        {
            ["anode"] = PowderAnodeIds[_powAnode], ["kalpha2"] = _powKa2, ["two_theta_min"] = Math.Min(_powFrom, _powTo - 1), ["two_theta_max"] = _powTo,
            ["step"] = 0.02, ["profile"] = PowderProfileIds[_powProfile], ["W"] = _powW, ["size"] = _powSizeNm * 10, ["strain"] = _powStrainPct / 100,
            ["b_iso"] = _powB, ["frames"] = _powAllFrames ? "all" : "current", ["max_peaks"] = 60,
        };
        var pref = _powPreferred.Split([' ', ',', ';'], StringSplitOptions.RemoveEmptyEntries);
        if (pref.Length == 3 && pref.All(x => int.TryParse(x, out _)) && pref.Any(x => x != "0"))
        {
            jo["preferred"] = new JsonArray(pref.Select(x => (JsonNode)int.Parse(x, inv)).ToArray());
            jo["march_r"] = _powR;
        }
        else if (_powPreferred.Trim().Length > 0) { PowError = "Preferred orientation: three integers h k l, e.g. 0 0 1"; return; }
        var doc = _doc!;
        PowderRunning = true;
        PowError = "";
        RefineHas = false;   // a new pattern replaces the last refinement on the plot
        Status = "Powder pattern…";
        try
        {
            var res = await Task.Run(() => doc.Powder(jo.ToJsonString()));
            var r = JsonNode.Parse(res)!;
            if (r["ok"]?.GetValue<bool>() != true) { PowError = (string?)r["error"] ?? "the powder pattern failed"; Status = "Powder pattern: " + PowError; return; }
            var x = ((JsonArray)r["two_theta"]!).Select(v => (double?)v ?? 0).ToArray();
            var y = ((JsonArray)r["intensity"]!).Select(v => (double?)v ?? 0).ToArray();
            PowderCurve = x.Zip(y).Select(p => (p.First, p.Second)).ToArray();
            _powLambda = r["wavelengths"] is JsonArray lw && lw.Count > 0 ? (double?)lw[0] ?? 1.54056 : 1.54056;
            PowderRows.Clear();
            foreach (var p in ((JsonArray)r["peaks"]!).Take(16))
            {
                double P(string k) => (double?)p![k] ?? 0;
                PowderRows.Add(new PowderRow($"{P("h"):0} {P("k"):0} {P("l"):0}", P("d").ToString("0.000", inv), P("two_theta").ToString("0.00", inv),
                                             P("multiplicity").ToString("0", inv), P("intensity").ToString("0.0", inv)));
            }
            PowText = string.Join(" · ", ((JsonArray)r["notes"]!).Select(v => (string?)v ?? ""));
            PowderView = true;
            Status = $"Powder pattern: {(int)((double?)r["reflections"] ?? 0)} reflections";
            ScatterChanged?.Invoke();
        }
        catch (Exception e) { PowError = e.Message; }
        finally { PowderRunning = false; }
    }
}
