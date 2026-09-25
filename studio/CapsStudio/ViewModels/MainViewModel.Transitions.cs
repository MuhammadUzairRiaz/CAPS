using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Threading.Tasks;

namespace CapsStudio.ViewModels;

/// <summary>Analyze › Glass transition (design/boards/GlassTransition): stepwise NPT cooling from the open structure,
/// replicas that differ only in their seed, specific volume against temperature with the replica spread, and the
/// two-line fit (Soldera &amp; Metatla 2006). Simulated cooling rates are reported beside every Tg.</summary>
public sealed partial class MainViewModel
{
    public bool IsGlass => _module == 47;
    public ResultCell GtTg { get; } = new("Tg", "two lines, free break");
    public ResultCell GtAlphaGlass { get; } = new("α glass", "volumetric expansion below Tg");
    public ResultCell GtAlphaMelt { get; } = new("α melt", "volumetric expansion above Tg");
    public ResultCell GtResidual { get; } = new("Fit residual", "rms of the two-line fit");
    private int _gtReplicas = 1;
    public decimal GtReplicas { get => _gtReplicas; set { if (Set(ref _gtReplicas, (int)Math.Clamp(value, 1, 10))) RaiseGlassPlan(); } }
    private string _gtStatus = "", _gtReplicaText = "";
    public string GtStatus { get => _gtStatus; private set => Set(ref _gtStatus, value); }
    public string GtReplicaText { get => _gtReplicaText; private set => Set(ref _gtReplicaText, value); }
    public (double X, double Y)[] GtPoints { get; private set; } = [];
    public double[] GtErrors { get; private set; } = [];
    public (double X, double Y)[] GtFit { get; private set; } = [];
    public bool GtHasPoints => GtPoints.Length > 1;
    public event Action? GtChanged;

    // the schedule from the settings (Analyze's cooling fields, shared with the Analyze page)
    private double GtFrom => (double)Analyze.TgFromD;
    private double GtTo => (double)Analyze.TgToD;
    private double GtStepK => (double)Analyze.TgStepD;
    private double GtHoldPs => (double)Analyze.TgPsD;
    private double GtAnnealPs => (double)Analyze.EqPsD;
    private int GtTemperatures => (int)Math.Floor(Math.Abs(GtFrom - GtTo) / Math.Max(1e-9, GtStepK) + 1e-9) + 1;
    public string GtTemperaturesText => GtTemperatures.ToString(CultureInfo.InvariantCulture);
    public string GtRateText => Sci(GtStepK / (GtHoldPs * 1e-12)) + " K/s";

    /// <summary>2.5e13 → "2.5 × 10¹³".</summary>
    public static string Sci(double v)
    {
        if (!(v > 0) || !double.IsFinite(v)) return "—";
        var e = (int)Math.Floor(Math.Log10(v));
        var m = v / Math.Pow(10, e);
        const string sup = "⁰¹²³⁴⁵⁶⁷⁸⁹";
        var exp = string.Concat(Math.Abs(e).ToString(CultureInfo.InvariantCulture).Select(c => sup[c - '0']));
        return m.ToString("0.0#", CultureInfo.InvariantCulture) + " × 10" + (e < 0 ? "⁻" : "") + exp;
    }
    public string GtTotalText
    {
        get
        {
            var ps = GtAnnealPs + GtTemperatures * GtHoldPs;
            var t = ps >= 1000 ? (ps / 1000).ToString("0.#", CultureInfo.InvariantCulture) + " ns" : ps.ToString("0", CultureInfo.InvariantCulture) + " ps";
            return _gtReplicas > 1 ? $"{t} × {_gtReplicas}" : t;
        }
    }
    /// <summary>Temperature against time for one replica: the anneal at the start, then each hold.</summary>
    public (double X, double Y)[] GtSchedule
    {
        get
        {
            var pts = new List<(double, double)>();
            double t = 0;
            var dir = GtTo < GtFrom ? -1 : 1;
            var unit = GtAnnealPs + GtTemperatures * GtHoldPs >= 1000 ? 1000.0 : 1.0;   // ns, or ps for short scans
            if (GtAnnealPs > 0) { pts.Add((0, GtFrom)); t = GtAnnealPs / unit; pts.Add((t, GtFrom)); }
            for (var k = 0; k < GtTemperatures; ++k)
            {
                var T = GtFrom + dir * k * GtStepK;
                pts.Add((t, T));
                t += GtHoldPs / unit;
                pts.Add((t, T));
            }
            return pts.ToArray();
        }
    }
    public string GtScheduleAxis => GtAnnealPs + GtTemperatures * GtHoldPs >= 1000 ? "time (ns) · one replica" : "time (ps) · one replica";
    public void RaiseGlassPlan()
    {
        foreach (var n in new[] { nameof(GtTemperaturesText), nameof(GtRateText), nameof(GtTotalText), nameof(GtSchedule), nameof(GtScheduleAxis) }) Raise(n);
        GtChanged?.Invoke();
    }

    public void OpenGlass() { SetModule(47); RaiseGlassPlan(); }

    /// <summary>Runs Analyze's cooling scan once per replica (seed 1, 2, …) and pools them.</summary>
    public async Task RunGlass()
    {
        if (_doc == null || Analyze.Working) return;
        var inv = CultureInfo.InvariantCulture;
        var runs = new List<(double Tg, double Err, double AlphaLow, double AlphaHigh, double Vtg, double Rms, (double T, double V)[] Points)>();
        var seed0 = Analyze.MechSeed;
        try
        {
            for (var r = 1; r <= _gtReplicas; ++r)
            {
                GtStatus = $"Replica {r} of {_gtReplicas}: cooling {GtFrom.ToString(inv)} → {GtTo.ToString(inv)} K";
                Analyze.MechSeed = (ulong)r;
                await RunChips("tg");
                var card = Analyze.Results.FirstOrDefault(c => c.Id == "tg");
                var curve = card == null ? null : Analyze.Curves.FirstOrDefault(c => c.Property == card.Name && c.Markers);
                if (card == null || curve == null) { GtStatus = Analyze.Log; if (runs.Count == 0) return; break; }
                double X(string k) => card.Extra.Where(e => e.Key.StartsWith(k, StringComparison.Ordinal)).Select(e => e.Value).DefaultIfEmpty(double.NaN).First();
                runs.Add((card.Value, card.Error, X("expansion below"), X("expansion above"), X("specific volume at Tg"), X("fit residual"),
                          curve.X.Zip(curve.Y).ToArray()));
            }
        }
        finally { Analyze.MechSeed = seed0; }
        // the replicas pooled: specific volume per temperature (mean ± SD), Tg and α (mean ± SD, or the fit's bootstrap for one)
        var temps = runs[0].Points.Select(p => p.T).ToArray();
        GtPoints = temps.Select((T, i) => (T, runs.Average(r => r.Points.Length > i ? r.Points[i].V : double.NaN))).ToArray();
        GtErrors = runs.Count > 1 ? temps.Select((_, i) => Sd(runs.Select(r => r.Points.Length > i ? r.Points[i].V : double.NaN).ToArray())).ToArray() : [];
        var ok = runs.Where(r => double.IsFinite(r.Tg)).ToList();
        if (ok.Count == 0) { GtTg.Value = "—"; GtStatus = "No break in the specific volume: widen the temperature range"; Raise(nameof(GtHasPoints)); GtChanged?.Invoke(); return; }
        var tg = ok.Average(r => r.Tg);
        var tgErr = ok.Count > 1 ? Sd(ok.Select(r => r.Tg).ToArray()) : ok[0].Err;
        var aLow = ok.Average(r => r.AlphaLow);
        var aHigh = ok.Average(r => r.AlphaHigh);
        var vtg = ok.Average(r => r.Vtg);
        GtTg.Value = $"{tg.ToString("0", inv)} ± {tgErr.ToString("0", inv)} K";
        GtTg.Caption = ok.Count > 1 ? $"mean ± SD of {ok.Count} replicas · at {GtRateText}" : $"bootstrap of the fit · at {GtRateText}";
        GtAlphaGlass.Value = $"{aLow.ToString("0.0e0", inv)} K⁻¹";
        GtAlphaMelt.Value = $"{aHigh.ToString("0.0e0", inv)} K⁻¹";
        GtResidual.Value = $"{ok.Average(r => r.Rms).ToString("0.0e0", inv)} cm³/g";
        // the pooled two-line fit: v(T) = v(Tg) + s_low·min(T − Tg, 0) + s_high·max(T − Tg, 0), s = α·v(Tg)
        var (lo, hi) = (temps.Min(), temps.Max());
        GtFit = Enumerable.Range(0, 101).Select(k => lo + (hi - lo) * k / 100.0)
                          .Select(T => (T, vtg + aLow * vtg * Math.Min(T - tg, 0) + aHigh * vtg * Math.Max(T - tg, 0))).ToArray();
        GtReplicaText = ok.Count > 1 ? "Replicas: " + string.Join(" · ", ok.Select(r => r.Tg.ToString("0", inv) + " K")) : "";
        GtStatus = $"{runs.Count} replica{(runs.Count == 1 ? "" : "s")} · {temps.Length} temperatures · {Analyze.Log}";
        Raise(nameof(GtHasPoints));
        GtChanged?.Invoke();
    }

    private static double Sd(double[] v)
    {
        var x = v.Where(double.IsFinite).ToArray();
        if (x.Length < 2) return double.NaN;
        var m = x.Average();
        return Math.Sqrt(x.Sum(a => (a - m) * (a - m)) / (x.Length - 1));
    }

    /// <summary>Runs only these Analyze calculations, then puts the chips back as they were.</summary>
    private async Task RunChips(params string[] ids)
    {
        var chips = Analyze.Groups.SelectMany(c => c.Chips).ToList();
        var was = chips.ToDictionary(c => c, c => c.IsOn);
        foreach (var c in chips) c.IsOn = ids.Contains(c.Id);
        try { await Analyze.Run(); }
        finally { foreach (var (c, on) in was) c.IsOn = on; }
    }

    /// <summary>The cooling scan as a recipe (analyze: tg) for this structure, saved beside it or anywhere.</summary>
    public string GlassRecipe(string structurePath)
    {
        var inv = CultureInfo.InvariantCulture;
        return "# Glass transition from CAPS Studio · run: caps run " + System.IO.Path.GetFileNameWithoutExtension(structurePath) + "_tg.yaml\n" +
               "recipe: 1\nname: " + System.IO.Path.GetFileNameWithoutExtension(structurePath) + "_tg\n" +
               "build: { file: \"" + structurePath.Replace("\\", "/") + "\" }\ntype: { forcefield: default }\n" +
               string.Format(inv, "analyze:\n  properties: [tg]\n  tg: {{ t_start: {0}, t_end: {1}, t_step: {2}, ps_per_step: {3}, equilibrate_ps: {4}, seed: 1 }}\n",
                             GtFrom, GtTo, GtStepK, GtHoldPs, GtAnnealPs > 0 ? GtAnnealPs : -1);
    }
}
