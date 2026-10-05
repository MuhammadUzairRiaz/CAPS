using System.Collections.ObjectModel;
using System.Globalization;
using System.Text.Json.Nodes;

namespace CapsStudio.ViewModels;

public sealed record SorbRow(string Pressure, string Loading, string MolKg, string Cm3, string Heat);
/// <summary>One gas of a mixture: its Widom solubility and, at the highest pressure, its loading, selectivity and heat.</summary>
public sealed record SorbSpeciesRow(string Smiles, string Y, string S, string Loading, string Selectivity, string Heat);

/// <summary>Analyze › Sorption: how much of a gas or small molecule the open structure takes up (core sorption.hpp):
/// Widom insertion (excess chemical potential, Henry constant, solubility coefficient) and grand-canonical Monte Carlo
/// at each pressure (the isotherm, the isosteric heat). The structure is held fixed and is not changed.</summary>
public sealed partial class MainViewModel
{
    public bool IsSorption => _module == 71;
    private bool _sorbRunning;
    private string _sorbate = "O=C=O", _sorbPressures = "10, 50, 100, 500, 1000", _sorbError = "", _sorbProgress = "", _sorbNote = "";
    private double _sorbT = 300;
    private int _sorbInsert = 100000, _sorbSteps = 200000;
    private string _sorbW = "—", _sorbMu = "—", _sorbHenry = "—", _sorbS = "—";
    private CancellationTokenSource? _sorbCts;
    public static readonly string[] SorbPresets = ["O=C=O · CO₂", "C · CH₄", "N#N · N₂", "O=O · O₂", "O · H₂O", "[He] · He", "[H][H] · H₂"];
    public ObservableCollection<SorbRow> SorbRows { get; } = new();
    // ---- gas mixtures (C6): one "SMILES fraction" per line, each gas at y_i p
    private bool _sorbMix;
    private string _sorbMixText = "O=C=O 0.15\nN#N 0.85";
    public bool SorbMix { get => _sorbMix; set => Set(ref _sorbMix, value); }
    public string SorbMixText { get => _sorbMixText; set { Set(ref _sorbMixText, value ?? ""); Raise(nameof(SorbMixHint)); } }
    public ObservableCollection<SorbSpeciesRow> SorbSpeciesRows { get; } = new();
    public bool SorbHasSpecies => SorbSpeciesRows.Count > 1;
    public string SorbSpeciesTitle { get; private set; } = "";
    /// <summary>Per species' isotherms (mol/kg against kPa) for the plot, in the mixture's order.</summary>
    public List<(double X, double Y)[]> SorbSpeciesIsotherms { get; } = new();
    /// <summary>The mixture lines as (SMILES, fraction); an error text when a line does not read.</summary>
    public static (List<(string Smiles, double Y)> Gases, string Error) ParseSorbMixture(string text)
    {
        var list = new List<(string, double)>();
        foreach (var raw in (text ?? "").Split(['\n', ';'], StringSplitOptions.RemoveEmptyEntries))
        {
            var line = raw.Trim();
            if (line.Length == 0 || line.StartsWith('#')) continue;
            var parts = line.Split([' ', '\t', ','], StringSplitOptions.RemoveEmptyEntries);
            if (parts.Length != 2 || !double.TryParse(parts[1].TrimEnd('%'), NumberStyles.Float, CultureInfo.InvariantCulture, out var y) || y < 0)
                return (list, $"'{line}': write a SMILES and its mole fraction, e.g. O=C=O 0.15");
            list.Add((parts[0], parts[1].EndsWith('%') ? y / 100 : y));
        }
        if (list.Count < 2) return (list, "a mixture needs two gases or more");
        if (list.Sum(g => g.Item2) <= 0) return (list, "the mole fractions add up to zero");
        return (list, "");
    }
    public string SorbMixHint
    {
        get
        {
            var (g, err) = ParseSorbMixture(_sorbMixText);
            if (err.Length > 0) return err;
            var t = g.Sum(x => x.Y);
            return string.Join(" · ", g.Select(x => $"{x.Smiles} y = {(x.Y / t).ToString("0.###", CultureInfo.InvariantCulture)}")) + (Math.Abs(t - 1) > 1e-6 ? " (normalised)" : "");
        }
    }
    public (double X, double Y)[] SorbIsotherm { get; private set; } = [];
    public string Sorbate { get => _sorbate; set => Set(ref _sorbate, value ?? ""); }
    public string SorbPreset { get => ""; set { if (value?.Split(" · ")[0] is { Length: > 0 } smi) Sorbate = smi; } }
    public string SorbPressures { get => _sorbPressures; set => Set(ref _sorbPressures, value ?? ""); }
    public decimal SorbTD { get => (decimal)_sorbT; set => Set(ref _sorbT, Math.Clamp((double)value, 1, 5000)); }
    public decimal SorbInsertD { get => _sorbInsert; set => Set(ref _sorbInsert, (int)Math.Clamp(value, 0, 100_000_000)); }
    public decimal SorbStepsD { get => _sorbSteps; set => Set(ref _sorbSteps, (int)Math.Clamp(value, 1000, 1_000_000_000)); }
    public bool SorbRunning { get => _sorbRunning; private set { if (Set(ref _sorbRunning, value)) RaiseBusy(); } }
    public string SorbError { get => _sorbError; private set { if (Set(ref _sorbError, value)) Raise(nameof(SorbHasError)); } }
    public bool SorbHasError => _sorbError.Length > 0;
    public string SorbProgress { get => _sorbProgress; private set => Set(ref _sorbProgress, value); }
    public string SorbNote { get => _sorbNote; private set => Set(ref _sorbNote, value); }
    public string SorbW { get => _sorbW; private set => Set(ref _sorbW, value); }
    public string SorbMu { get => _sorbMu; private set => Set(ref _sorbMu, value); }
    public string SorbHenry { get => _sorbHenry; private set => Set(ref _sorbHenry, value); }
    public string SorbS { get => _sorbS; private set => Set(ref _sorbS, value); }
    public event Action? SorptionChanged;

    public void OpenSorption() => SetModule(71);
    public void CancelSorption() => _sorbCts?.Cancel();

    public async Task RunSorption()
    {
        if (_doc == null || !Idle || BlockedByField("Sorption")) return;
        var inv = CultureInfo.InvariantCulture;
        var pressures = new JsonArray();
        foreach (var p in SorbPressureList(_sorbPressures)) pressures.Add(p);
        var jo = new JsonObject
        {
            ["sorbate"] = _sorbate.Trim(), ["temperature"] = _sorbT, ["insertions"] = _sorbInsert, ["pressures_kpa"] = pressures,
            ["steps"] = _sorbSteps, ["cutoff"] = Math.Min(12.0, _relaxCutoff), ["coulomb"] = _relaxCoulomb, ["seed"] = SorbSeed.Take(),
            ["map_grid"] = _sorbMapOn ? 24 : 0,
        };
        if (_sorbMix)
        {
            var (gases, err) = ParseSorbMixture(_sorbMixText);
            if (err.Length > 0) { SorbError = err; return; }
            var mix = new JsonArray();
            foreach (var (smi, y) in gases) mix.Add(new JsonObject { ["smiles"] = smi, ["fraction"] = y });
            jo["mixture"] = mix;
        }
        var json = jo.ToJsonString();
        var doc = _doc!;
        var cellA = doc.Summary();
        SorbRunning = true;
        SorbError = "";
        _sorbCts = new CancellationTokenSource();
        var ct = _sorbCts.Token;
        Status = $"Sorption of {_sorbate} in {Title}…";
        try
        {
            var res = await Task.Run(() => doc.Sorption(json, (st, f) =>
            {
                Avalonia.Threading.Dispatcher.UIThread.Post(() => SorbProgress = $"{st} · {f * 100:0} %");
                return !ct.IsCancellationRequested;
            }));
            var r = JsonNode.Parse(res)!;
            if (r["ok"]?.GetValue<bool>() != true) { SorbError = (string?)r["error"] ?? "sorption failed"; Status = "Sorption: " + SorbError; return; }
            double D(string k) => r[k]?.GetValue<double>() ?? 0;
            SorbW = $"{D("widom_w").ToString("0.####E+0", inv)} ± {D("widom_error").ToString("0.#E+0", inv)}";
            SorbMu = D("mu_ex").ToString("0.00", inv) + " kcal/mol";
            SorbHenry = D("henry_mol_kg_kpa").ToString("0.###E+0", inv) + " mol/(kg·kPa)";
            SorbS = D("solubility").ToString("0.###", inv) + " cm³(STP)/(cm³·atm)";
            SorbRows.Clear();
            _sorbMaps.Clear();
            var pts = new List<(double, double)>();
            foreach (var p in (JsonArray)r["isotherm"]!)
            {
                double P(string k) => p![k]!.GetValue<double>();
                SorbRows.Add(new SorbRow(P("pressure_kpa").ToString("0.##", inv), $"{P("loading").ToString("0.##", inv)} ± {P("loading_error").ToString("0.##", inv)}",
                    P("mol_per_kg").ToString("0.###", inv), P("cm3stp_per_cm3").ToString("0.##", inv), P("heat").ToString("0.0", inv)));
                pts.Add((P("pressure_kpa"), P("mol_per_kg")));
                if (p!["map"] is JsonObject mp && mp["density"] is JsonArray da)
                    _sorbMaps.Add((P("pressure_kpa"), (int)((double?)mp["grid"] ?? 0), da.Select(v => (double?)v ?? 0).ToArray()));
            }
            _sorbCell = (cellA.CellA, cellA.CellB, cellA.CellC);
            SorbMapPressures.Clear();
            foreach (var m in _sorbMaps) SorbMapPressures.Add(m.P.ToString("0.##", inv) + " kPa");
            _sorbMapIndex = Math.Max(0, _sorbMaps.Count - 1);
            Raise(nameof(SorbMapIndex)); Raise(nameof(SorbHasMap));
            SorptionMapChanged?.Invoke();
            SorbIsotherm = pts.ToArray();
            SorbSpeciesRows.Clear();
            SorbSpeciesIsotherms.Clear();
            if (r["species"] is JsonArray spj && spj.Count > 1 && r["isotherm"] is JsonArray isoj)
            {
                var lastPt = isoj.Count > 0 ? isoj[^1] : null;
                double At(JsonNode? p, string k, int i) => p?[k] is JsonArray a && i < a.Count ? (double?)a[i] ?? double.NaN : double.NaN;
                for (var i = 0; i < spj.Count; ++i)
                {
                    var sp = spj[i]!;
                    SorbSpeciesRows.Add(new SorbSpeciesRow((string?)sp["smiles"] ?? "", ((double?)sp["fraction"] ?? 0).ToString("0.###", inv),
                        ((double?)sp["solubility"] ?? 0).ToString("0.###", inv),
                        lastPt == null ? "—" : $"{At(lastPt, "species_loading", i).ToString("0.##", inv)} ± {At(lastPt, "species_error", i).ToString("0.##", inv)}",
                        lastPt == null ? "—" : At(lastPt, "selectivity", i).ToString("0.###", inv),
                        lastPt == null ? "—" : At(lastPt, "species_heat", i).ToString("0.0", inv)));
                    SorbSpeciesIsotherms.Add(isoj.Select(p => ((double?)p!["pressure_kpa"] ?? 0, At(p, "species_mol_per_kg", i))).Where(q => double.IsFinite(q.Item2)).ToArray());
                }
                SorbSpeciesTitle = lastPt == null ? "Gases (Widom)" : $"Gases at {((double?)lastPt["pressure_kpa"] ?? 0).ToString("0.##", inv)} kPa · selectivity over {(string?)spj[0]!["smiles"]}";
            }
            Raise(nameof(SorbHasSpecies)); Raise(nameof(SorbSpeciesTitle));
            SorbNote = string.Join(" · ", ((JsonArray)r["notes"]!).Select(x => (string?)x ?? "")) + " · " + ((string?)r["forcefield"] ?? "");
            Status = $"Sorption of {(_sorbMix ? "the mixture" : _sorbate)}: S {SorbS}" + (SorbRows.Count > 0 ? $" · {SorbRows.Count} isotherm points" : "");
            SorptionChanged?.Invoke();
        }
        catch (Exception e) { SorbError = e.Message; }
        finally { SorbRunning = false; SorbProgress = ""; }
    }

    // ---- the sorbate's density map (C7): per pressure, the 3D grid projected onto a cell face as molecules per nm²
    private bool _sorbMapOn = true;
    private readonly List<(double P, int Grid, double[] Density)> _sorbMaps = new();
    private (double A, double B, double C) _sorbCell;
    private int _sorbMapIndex, _sorbMapFace;
    public bool SorbMapOn { get => _sorbMapOn; set => Set(ref _sorbMapOn, value); }
    public System.Collections.ObjectModel.ObservableCollection<string> SorbMapPressures { get; } = new();
    public static readonly string[] SorbMapFaces = ["Onto the a–b face (along c)", "Onto the a–c face (along b)", "Onto the b–c face (along a)"];
    public bool SorbHasMap => _sorbMaps.Count > 0;
    public int SorbMapIndex { get => _sorbMapIndex; set { if (Set(ref _sorbMapIndex, Math.Clamp(value, 0, Math.Max(0, _sorbMaps.Count - 1)))) SorptionMapChanged?.Invoke(); } }
    public int SorbMapFace { get => _sorbMapFace; set { if (Set(ref _sorbMapFace, Math.Clamp(value, 0, 2))) SorptionMapChanged?.Invoke(); } }
    public event Action? SorptionMapChanged;
    private string _sorbMapText = "";
    public string SorbMapText { get => _sorbMapText; private set => Set(ref _sorbMapText, value); }
    /// <summary>The chosen pressure's map summed along one cell edge: (x, y, areal density per nm²) for the heat plot.</summary>
    public (double[] X, double[] Y, double[] Z, string XLabel, string YLabel) SorbMapProjection()
    {
        if (_sorbMapIndex < 0 || _sorbMapIndex >= _sorbMaps.Count) return ([], [], [], "", "");
        var (p, g, d) = _sorbMaps[_sorbMapIndex];
        if (g <= 0 || d.Length != g * g * g) return ([], [], [], "", "");
        var L = new[] { _sorbCell.A, _sorbCell.B, _sorbCell.C };
        var (u, v, w) = _sorbMapFace switch { 1 => (0, 2, 1), 2 => (1, 2, 0), _ => (0, 1, 2) };
        var xs = new List<double>(); var ys = new List<double>(); var zs = new List<double>();
        var peak = 0.0;
        for (var a = 0; a < g; ++a)
            for (var b = 0; b < g; ++b)
            {
                var sum = 0.0;
                for (var c = 0; c < g; ++c)
                {
                    var ix = new int[3];
                    ix[u] = a; ix[v] = b; ix[w] = c;
                    sum += d[(ix[0] * g + ix[1]) * g + ix[2]];
                }
                var areal = sum * L[w] / g * 100;   // per Å³ × Å along the edge → per Å² → per nm²
                xs.Add((a + 0.5) * L[u] / g); ys.Add((b + 0.5) * L[v] / g); zs.Add(areal);
                peak = Math.Max(peak, areal);
            }
        var names = new[] { "a", "b", "c" };
        SorbMapText = string.Format(CultureInfo.InvariantCulture, "{0:0.##} kPa · summed along {1} · up to {2:0.###} molecules/nm² · {3} × {3} × {3} grid; empty regions are where the host leaves no room",
            p, names[w], peak, g);
        return (xs.ToArray(), ys.ToArray(), zs.ToArray(), names[u] + " (Å)", names[v] + " (Å)");
    }

    /// <summary>The GCMC pressures: a list ("10, 50, 100") and/or log-spaced sweeps "a..b xN" (N points from a to b, evenly in
    /// log p — isotherms span decades).</summary>
    public static List<double> SorbPressureList(string text)
    {
        var inv = CultureInfo.InvariantCulture;
        var outp = new List<double>();
        var rest = System.Text.RegularExpressions.Regex.Replace(text, @"([0-9.eE+-]+)\s*\.\.\s*([0-9.eE+-]+)\s*[xX×]\s*(\d+)", m =>
        {
            if (double.TryParse(m.Groups[1].Value, NumberStyles.Float, inv, out var a) && double.TryParse(m.Groups[2].Value, NumberStyles.Float, inv, out var b) &&
                int.TryParse(m.Groups[3].Value, out var n) && a > 0 && b > 0 && n >= 2)
                for (var k = 0; k < Math.Min(n, 50); ++k) outp.Add(a * Math.Pow(b / a, k / (double)(n - 1)));
            return " ";
        });
        foreach (var t in rest.Split([',', ' ', ';'], StringSplitOptions.RemoveEmptyEntries))
            if (double.TryParse(t, NumberStyles.Float, inv, out var p) && p > 0) outp.Add(p);
        return outp.Distinct().OrderBy(v => v).ToList();
    }
}
